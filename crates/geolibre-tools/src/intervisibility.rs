//! GeoLibre tool: line-of-sight visibility between observer/target pairs
//! against 3D vector obstructions.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Intervisibility* (3D Analyst).
//!
//! ## Why `line_of_sight` does not cover this
//!
//! `line_of_sight` (round 2) tests against a **surface raster** and returns a
//! profile with obstruction points. `visibility_index`, `skyline_analysis` and
//! `viewshed` are likewise raster-DEM tools. None of them can test visibility
//! against **3D vector obstructions** — buildings as multipatches — which is
//! the whole point of Intervisibility: a DEM does not contain buildings, and
//! burning them in loses their overhangs, arcades and interior voids.
//!
//! Batch pair testing against solids has no counterpart in either registry.
//!
//! ## Why a segment test rather than the ray machinery
//!
//! `inside_3d` casts an infinite ray and counts crossings for parity, which is
//! why it needs the shared-edge dedup fix from round 16. Here the question is
//! "does anything lie between these two points", so the test is a **segment**
//! against each triangle with the parameter clamped to `[0, 1]`. There is no
//! parity to break, so the diagonal double-hit hazard does not apply at all.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::args_common::{bool_or, f64_or, req_str};
use crate::inside_3d::{collect_triangles, Solid};
use crate::mesh3d::segment_triangle;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct IntervisibilityTool;

impl Tool for IntervisibilityTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "intervisibility",
            display_name: "Intervisibility",
            summary: "Tests each sight line against 3D vector obstructions and flags it visible or blocked, reporting which obstruction blocked it and where (ArcGIS Intervisibility). line_of_sight, visibility_index, viewshed and skyline_analysis all work against a surface raster; none can test against 3D multipatch obstructions such as buildings, whose overhangs and arcades are lost when burnt into a DEM.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Sight lines: 3D polylines whose first vertex is the observer and last vertex the target.",
                    required: true,
                },
                ToolParamSpec {
                    name: "obstructions",
                    description: "One or more comma-separated 3D obstruction layers (triangle-mesh MultiPolygons with Z).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output sight lines with the visibility fields appended. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "visible_field",
                    description: "Name of the boolean visibility field (default 'VISIBLE'), matching the ArcGIS parameter.",
                    required: false,
                },
                ToolParamSpec {
                    name: "observer_offset",
                    description: "Height added to the observer endpoint before testing, in CRS units (default 0). Models eye or instrument height.",
                    required: false,
                },
                ToolParamSpec {
                    name: "target_offset",
                    description: "Height added to the target endpoint before testing, in CRS units (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "visible_only",
                    description: "Emit only the sight lines that are unobstructed (default false).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        req_str(args, "obstructions")?;
        f64_or(args, "observer_offset", 0.0)?;
        f64_or(args, "target_offset", 0.0)?;
        bool_or(args, "visible_only", false)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?;
        let obstruction_spec = req_str(args, "obstructions")?;
        let output = parse_optional_str(args, "output")?;
        let visible_field = parse_optional_str(args, "visible_field")?
            .unwrap_or("VISIBLE")
            .to_string();
        let observer_offset = f64_or(args, "observer_offset", 0.0)?;
        let target_offset = f64_or(args, "target_offset", 0.0)?;
        let visible_only = bool_or(args, "visible_only", false)?;

        // Load every obstruction layer into one flat list, tagged by source so
        // BLOCKED_BY identifies the feature rather than just "something".
        let mut blockers: Vec<(usize, Solid)> = Vec::new();
        for (li, path) in obstruction_spec
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .enumerate()
        {
            let layer = load_input_layer(path)?;
            for (fid, feature) in layer.iter().enumerate() {
                let Some(geom) = feature.geometry.as_ref() else {
                    continue;
                };
                let tris = collect_triangles(geom);
                if tris.is_empty() {
                    continue;
                }
                // Obstructions need not be closed: an open wall blocks sight
                // just as well as a solid one, so no closedness filter here.
                blockers.push((li, Solid::new(fid, tris)));
            }
        }
        if blockers.is_empty() {
            return Err(ToolError::Execution(
                "'obstructions' holds no triangle-mesh geometry".to_string(),
            ));
        }

        let lines = load_input_layer(input)?;
        ctx.progress.info(&format!(
            "{} obstruction mesh(es), {} sight line(s)",
            blockers.len(),
            lines.iter().count()
        ));

        let mut out = Layer::new("intervisibility");
        out.geom_type = lines.geom_type;
        out.crs = lines.crs.clone();
        for f in lines.schema.fields() {
            out.add_field(f.clone());
        }
        out.add_field(FieldDef::new(&visible_field, FieldType::Boolean));
        out.add_field(FieldDef::new("BLOCKED_BY", FieldType::Integer));
        out.add_field(FieldDef::new("BLOCKED_LAYER", FieldType::Integer));
        out.add_field(FieldDef::new("BLOCK_FRACTION", FieldType::Float));

        let names: Vec<String> = lines
            .schema
            .fields()
            .iter()
            .map(|f| f.name.clone())
            .collect();

        let mut visible_count = 0_u64;
        let mut blocked_count = 0_u64;
        let mut skipped = 0_u64;
        let total = lines.iter().count().max(1);

        for (i, feature) in lines.iter().enumerate() {
            let Some((mut a, mut b)) = endpoints(feature.geometry.as_ref()) else {
                skipped += 1;
                continue;
            };
            a[2] += observer_offset;
            b[2] += target_offset;

            // Nearest blocker along the segment, so BLOCK_FRACTION is the
            // first obstruction rather than an arbitrary one.
            let mut best: Option<(f64, usize, usize)> = None;
            for (li, solid) in &blockers {
                if !segment_bbox_overlap(a, b, solid) {
                    continue;
                }
                for tri in &solid.tris {
                    if let Some(t) = segment_triangle(a, b, tri) {
                        // Ignore hits pinned exactly at the endpoints: an
                        // observer standing on a rooftop is not blocked by it.
                        if !(1e-9..=1.0 - 1e-9).contains(&t) {
                            continue;
                        }
                        if best.map_or(true, |(bt, _, _)| t < bt) {
                            best = Some((t, solid.fid, *li));
                        }
                    }
                }
            }

            let visible = best.is_none();
            if visible {
                visible_count += 1;
            } else {
                blocked_count += 1;
            }
            if visible_only && !visible {
                continue;
            }

            let mut attrs: Vec<(&str, FieldValue)> = names
                .iter()
                .enumerate()
                .filter_map(|(k, n)| {
                    feature.attributes.get(k).map(|v| (n.as_str(), v.clone()))
                })
                .collect();
            attrs.push((visible_field.as_str(), FieldValue::Boolean(visible)));
            attrs.push((
                "BLOCKED_BY",
                match best {
                    Some((_, fid, _)) => FieldValue::Integer(fid as i64),
                    None => FieldValue::Integer(-1),
                },
            ));
            attrs.push((
                "BLOCKED_LAYER",
                match best {
                    Some((_, _, li)) => FieldValue::Integer(li as i64),
                    None => FieldValue::Integer(-1),
                },
            ));
            attrs.push((
                "BLOCK_FRACTION",
                FieldValue::Float(best.map(|(t, _, _)| t).unwrap_or(-1.0)),
            ));
            out.add_feature(feature.geometry.clone(), &attrs)
                .map_err(|e| ToolError::Execution(e.to_string()))?;
            ctx.progress.progress((i as f64 + 1.0) / total as f64);
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("visible_count".to_string(), json!(visible_count));
        outputs.insert("blocked_count".to_string(), json!(blocked_count));
        outputs.insert("obstruction_count".to_string(), json!(blockers.len()));
        outputs.insert("skipped_count".to_string(), json!(skipped));
        Ok(ToolRunResult { outputs })
    }
}

/// First and last vertex of a sight line, as 3D points.
fn endpoints(geom: Option<&Geometry>) -> Option<([f64; 3], [f64; 3])> {
    let coords: &Vec<Coord> = match geom? {
        Geometry::LineString(cs) => cs,
        Geometry::MultiLineString(parts) => parts.first()?,
        _ => return None,
    };
    let a = coords.first()?;
    let b = coords.last()?;
    if coords.len() < 2 {
        return None;
    }
    Some((
        [a.x, a.y, a.z.unwrap_or(0.0)],
        [b.x, b.y, b.z.unwrap_or(0.0)],
    ))
}

/// Cheap rejection: does the segment's bounding box meet the solid's?
fn segment_bbox_overlap(a: [f64; 3], b: [f64; 3], solid: &Solid) -> bool {
    (0..3).all(|k| {
        let lo = a[k].min(b[k]);
        let hi = a[k].max(b[k]);
        lo <= solid.max[k] && hi >= solid.min[k]
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, GeometryType};

    use crate::mesh3d::box_mesh;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn solids(geoms: Vec<Geometry>) -> String {
        let mut l = Layer::new("obs");
        l.geom_type = Some(GeometryType::MultiPolygon);
        for g in geoms {
            l.add_feature(Some(g), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    /// One sight line per (observer, target) pair.
    fn sightlines(pairs: Vec<([f64; 3], [f64; 3])>) -> String {
        let mut l = Layer::new("los");
        l.geom_type = Some(GeometryType::LineString);
        for (a, b) in pairs {
            l.add_feature(
                Some(Geometry::LineString(vec![
                    Coord::xyz(a[0], a[1], a[2]),
                    Coord::xyz(b[0], b[1], b[2]),
                ])),
                &[],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = IntervisibilityTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    fn val(layer: &Layer, fid: usize, name: &str) -> FieldValue {
        let i = layer.schema.field_index(name).unwrap();
        layer.iter().nth(fid).unwrap().attributes[i].clone()
    }

    /// A wall spanning y and z, one unit thick in x at x = 5.
    fn wall() -> Geometry {
        box_mesh([5.0, -10.0, 0.0], [6.0, 10.0, 10.0])
    }

    #[test]
    fn a_wall_between_observer_and_target_blocks_the_sight_line() {
        let (out, res) = run(json!({
            "input": sightlines(vec![([0.0, 0.0, 5.0], [10.0, 0.0, 5.0])]),
            "obstructions": solids(vec![wall()]),
        }));
        assert_eq!(val(&out, 0, "VISIBLE"), FieldValue::Boolean(false));
        assert_eq!(res.outputs["blocked_count"], json!(1));
        assert_eq!(res.outputs["visible_count"], json!(0));
    }

    #[test]
    fn a_sight_line_passing_over_the_wall_stays_visible() {
        // The wall tops out at z = 10; this line runs at z = 20.
        let (out, res) = run(json!({
            "input": sightlines(vec![([0.0, 0.0, 20.0], [10.0, 0.0, 20.0])]),
            "obstructions": solids(vec![wall()]),
        }));
        assert_eq!(val(&out, 0, "VISIBLE"), FieldValue::Boolean(true));
        assert_eq!(res.outputs["visible_count"], json!(1));
        assert_eq!(val(&out, 0, "BLOCKED_BY"), FieldValue::Integer(-1));
    }

    #[test]
    fn a_sight_line_that_stops_short_of_the_wall_is_not_blocked() {
        // The segment ends at x = 3, well before the wall at x = 5. An
        // infinite-ray test would wrongly call this blocked.
        let (out, _) = run(json!({
            "input": sightlines(vec![([0.0, 0.0, 5.0], [3.0, 0.0, 5.0])]),
            "obstructions": solids(vec![wall()]),
        }));
        assert_eq!(val(&out, 0, "VISIBLE"), FieldValue::Boolean(true));
    }

    #[test]
    fn the_blocking_feature_and_position_are_reported() {
        let (out, _) = run(json!({
            "input": sightlines(vec![([0.0, 0.0, 5.0], [10.0, 0.0, 5.0])]),
            "obstructions": solids(vec![wall()]),
        }));
        assert_eq!(val(&out, 0, "BLOCKED_BY"), FieldValue::Integer(0));
        // First contact is the wall's near face at x = 5 along a 10-unit run.
        match val(&out, 0, "BLOCK_FRACTION") {
            FieldValue::Float(t) => assert!((t - 0.5).abs() < 1e-6, "got {t}"),
            other => panic!("expected a float, got {other:?}"),
        }
    }

    #[test]
    fn the_nearest_obstruction_wins_when_several_intervene() {
        // Two walls; the reported blocker must be the closer one (fid 0 at
        // x = 2), not whichever happened to be tested first.
        let near = box_mesh([2.0, -10.0, 0.0], [3.0, 10.0, 10.0]);
        let far = box_mesh([7.0, -10.0, 0.0], [8.0, 10.0, 10.0]);
        let (out, _) = run(json!({
            "input": sightlines(vec![([0.0, 0.0, 5.0], [10.0, 0.0, 5.0])]),
            "obstructions": solids(vec![near, far]),
        }));
        assert_eq!(val(&out, 0, "BLOCKED_BY"), FieldValue::Integer(0));
        match val(&out, 0, "BLOCK_FRACTION") {
            FieldValue::Float(t) => assert!((t - 0.2).abs() < 1e-6, "got {t}"),
            other => panic!("expected a float, got {other:?}"),
        }
    }

    #[test]
    fn observer_offset_can_raise_a_line_over_an_obstruction() {
        // A low wall (to z = 4). At eye level 1 the line is blocked; raising
        // the observer to 9 clears it.
        let low = box_mesh([5.0, -10.0, 0.0], [6.0, 10.0, 4.0]);
        let obs = solids(vec![low]);
        // Both ends at z = 1, well under the wall's z = 4 top.
        let lines = sightlines(vec![([0.0, 0.0, 1.0], [10.0, 0.0, 1.0])]);
        let (blocked, _) = run(json!({"input": lines, "obstructions": obs}));
        assert_eq!(val(&blocked, 0, "VISIBLE"), FieldValue::Boolean(false));
        let (clear, _) = run(json!({
            "input": lines, "obstructions": obs, "observer_offset": 8.0,
        }));
        assert_eq!(val(&clear, 0, "VISIBLE"), FieldValue::Boolean(true));
    }

    #[test]
    fn visible_only_filters_the_output_but_not_the_counts() {
        let (out, res) = run(json!({
            "input": sightlines(vec![
                ([0.0, 0.0, 5.0], [10.0, 0.0, 5.0]),
                ([0.0, 0.0, 20.0], [10.0, 0.0, 20.0]),
            ]),
            "obstructions": solids(vec![wall()]),
            "visible_only": true,
        }));
        assert_eq!(out.iter().count(), 1);
        assert_eq!(res.outputs["visible_count"], json!(1));
        assert_eq!(res.outputs["blocked_count"], json!(1));
    }

    #[test]
    fn a_custom_visible_field_name_is_honoured() {
        let (out, _) = run(json!({
            "input": sightlines(vec![([0.0, 0.0, 20.0], [10.0, 0.0, 20.0])]),
            "obstructions": solids(vec![wall()]),
            "visible_field": "CAN_SEE",
        }));
        assert!(out.schema.field_index("CAN_SEE").is_some());
        assert_eq!(val(&out, 0, "CAN_SEE"), FieldValue::Boolean(true));
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            IntervisibilityTool.validate(&args).is_err()
        };
        assert!(bad(json!({"obstructions": "x"})));
        assert!(bad(json!({"input": "x"})));
    }
}
