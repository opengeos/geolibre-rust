//! GeoLibre tool: 2D footprint polygons from 3D multipatch features.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Multipatch Footprint* (3D Analyst).
//!
//! ## Why this matters for the repo specifically
//!
//! GeoLibre's identity is the raster-to-vector cleanup pipeline
//! (`polygonize` → `regularize_building_footprints` → `smooth_natural_features`
//! → web formats). That pipeline currently has **no 3D entry point**: solids
//! produced by `buffer_3d`, `voxel_isosurface` or `minimum_bounding_volume`, and
//! building meshes arriving from outside, cannot be fed into it at all.
//!
//! `layer_footprint_vector` returns one extent rectangle for an entire layer,
//! which is a different thing entirely — it cannot give per-building outlines.
//!
//! ## Why the triangles are unioned rather than hulled
//!
//! Projecting a mesh gives a pile of overlapping triangles. Their **union** is
//! the true footprint; a convex hull would fill in courtyards and L-shaped
//! plans, which is precisely the detail `regularize_building_footprints` exists
//! to work on. Vertical faces project to zero-area slivers and are dropped
//! before the union so they cannot inject degenerate rings.

use std::collections::BTreeMap;

use geo::{Area, BooleanOps, Coord as GeoCoord, LineString, MultiPolygon, Polygon, Simplify};
use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::args_common::{opt_positive_f64, req_str};
use crate::inside_3d::collect_triangles;
use crate::mesh3d::tri_area;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct MultipatchFootprintTool;

impl Tool for MultipatchFootprintTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "multipatch_footprint",
            display_name: "Multipatch Footprint",
            summary: "Projects 3D multipatch features onto the XY plane and dissolves them into 2D footprint polygons, one per feature or per group field (ArcGIS Multipatch Footprint). layer_footprint_vector returns a single extent rectangle for a whole layer, not per-feature outlines, so there is currently no way to feed buffer_3d / voxel_isosurface / imported building solids into the regularize_building_footprints pipeline.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "3D multipatch features (triangle-mesh MultiPolygons with Z).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output 2D polygon footprints. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "group_field",
                    description: "Optional field whose values dissolve footprints together (one output polygon per distinct value). Default: one output per input feature.",
                    required: false,
                },
                ToolParamSpec {
                    name: "simplify_tolerance",
                    description: "Optional Douglas-Peucker tolerance (CRS units) applied to the dissolved footprint. Omitted by default so the raw outline is preserved for regularize_building_footprints.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        opt_positive_f64(args, "simplify_tolerance")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let group_field = parse_optional_str(args, "group_field")?;
        let tolerance = opt_positive_f64(args, "simplify_tolerance")?;

        let layer = load_input_layer(input)?;
        let group_idx = match group_field {
            Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                ToolError::Validation(format!("group_field '{f}' not found on the input layer"))
            })?),
            None => None,
        };

        // Group the projected triangles. Insertion order is preserved for the
        // ungrouped case so outputs line up with input features.
        let mut groups: Vec<(String, Vec<Polygon>)> = Vec::new();
        let mut index: BTreeMap<String, usize> = BTreeMap::new();
        let mut degenerate = 0_u64;
        let mut skipped = 0_u64;

        for (fid, feature) in layer.iter().enumerate() {
            let tris = feature
                .geometry
                .as_ref()
                .map(collect_triangles)
                .unwrap_or_default();
            if tris.is_empty() {
                skipped += 1;
                continue;
            }
            let key = match group_idx {
                Some(i) => field_to_string(&feature.attributes[i]),
                None => fid.to_string(),
            };
            let slot = *index.entry(key.clone()).or_insert_with(|| {
                groups.push((key.clone(), Vec::new()));
                groups.len() - 1
            });
            for t in &tris {
                let poly = Polygon::new(
                    LineString::from(vec![
                        GeoCoord { x: t[0][0], y: t[0][1] },
                        GeoCoord { x: t[1][0], y: t[1][1] },
                        GeoCoord { x: t[2][0], y: t[2][1] },
                        GeoCoord { x: t[0][0], y: t[0][1] },
                    ]),
                    vec![],
                );
                // A vertical wall projects to a zero-area sliver; unioning
                // those in adds degenerate rings for no benefit.
                if poly.unsigned_area() <= f64::EPSILON {
                    if tri_area(t) > 0.0 {
                        degenerate += 1;
                    }
                    continue;
                }
                groups[slot].1.push(poly);
            }
        }

        ctx.progress
            .info(&format!("{} footprint group(s)", groups.len()));

        let mut out = Layer::new("multipatch_footprint");
        out.geom_type = Some(GeometryType::MultiPolygon);
        out.crs = layer.crs.clone();
        out.add_field(FieldDef::new("GROUP_ID", FieldType::Text));
        out.add_field(FieldDef::new("TRI_COUNT", FieldType::Integer));
        out.add_field(FieldDef::new("FOOTPRINT_AREA", FieldType::Float));

        let total = groups.len().max(1);
        let mut emitted = 0_u64;
        for (gi, (key, polys)) in groups.iter().enumerate() {
            if polys.is_empty() {
                continue;
            }
            // Accumulate with BooleanOps rather than hulling: courtyards and
            // re-entrant corners must survive.
            let mut acc = MultiPolygon::new(Vec::new());
            for p in polys {
                acc = acc.union(&MultiPolygon::new(vec![p.clone()]));
            }
            if let Some(tol) = tolerance {
                acc = MultiPolygon::new(
                    acc.0.iter().map(|p| p.simplify(tol)).collect::<Vec<_>>(),
                );
            }
            let area = acc.unsigned_area();
            if acc.0.is_empty() || area <= 0.0 {
                continue;
            }
            out.add_feature(
                Some(multipolygon_to_geometry(&acc)),
                &[
                    ("GROUP_ID", FieldValue::Text(key.clone())),
                    ("TRI_COUNT", FieldValue::Integer(polys.len() as i64)),
                    ("FOOTPRINT_AREA", FieldValue::Float(area)),
                ],
            )
            .map_err(|e| ToolError::Execution(e.to_string()))?;
            emitted += 1;
            ctx.progress.progress((gi as f64 + 1.0) / total as f64);
        }

        if emitted == 0 {
            return Err(ToolError::Execution(format!(
                "input holds no projectable triangle meshes ({skipped} non-mesh feature(s) skipped)"
            )));
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("footprint_count".to_string(), json!(emitted));
        outputs.insert("vertical_triangles".to_string(), json!(degenerate));
        outputs.insert("non_mesh_count".to_string(), json!(skipped));
        Ok(ToolRunResult { outputs })
    }
}

fn multipolygon_to_geometry(mp: &MultiPolygon) -> Geometry {
    Geometry::MultiPolygon(
        mp.0.iter()
            .map(|p| {
                (
                    linestring_to_ring(p.exterior()),
                    p.interiors().iter().map(linestring_to_ring).collect(),
                )
            })
            .collect(),
    )
}

fn linestring_to_ring(ls: &LineString) -> Ring {
    let mut coords: Vec<Coord> = ls.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
    if coords.len() >= 2 && coords.first() == coords.last() {
        coords.pop();
    }
    Ring::new(coords)
}

fn field_to_string(v: &FieldValue) -> String {
    match v {
        FieldValue::Text(s) => s.clone(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) => f.to_string(),
        FieldValue::Boolean(b) => b.to_string(),
        other => format!("{other:?}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::memory_store;

    use crate::inside_3d::box_mesh;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn layer_of(geoms: Vec<Geometry>) -> String {
        let mut l = Layer::new("in");
        l.geom_type = Some(GeometryType::MultiPolygon);
        for g in geoms {
            l.add_feature(Some(g), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = MultipatchFootprintTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    fn area(layer: &Layer, fid: usize) -> f64 {
        let i = layer.schema.field_index("FOOTPRINT_AREA").unwrap();
        match &layer.iter().nth(fid).unwrap().attributes[i] {
            FieldValue::Float(v) => *v,
            other => panic!("expected a float, got {other:?}"),
        }
    }

    #[test]
    fn a_box_projects_to_its_base_rectangle() {
        let (out, res) = run(json!({
            "input": layer_of(vec![box_mesh([0.0, 0.0, 0.0], [3.0, 5.0, 9.0])]),
        }));
        assert_eq!(res.outputs["footprint_count"], json!(1));
        // 3 x 5 = 15, regardless of the 9-unit height.
        assert!((area(&out, 0) - 15.0).abs() < 1e-6, "got {}", area(&out, 0));
    }

    #[test]
    fn one_footprint_per_feature_by_default() {
        let (out, res) = run(json!({
            "input": layer_of(vec![
                box_mesh([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]),
                box_mesh([10.0, 10.0, 0.0], [12.0, 12.0, 2.0]),
            ]),
        }));
        assert_eq!(res.outputs["footprint_count"], json!(2));
        assert_eq!(out.iter().count(), 2);
        assert!((area(&out, 0) - 4.0).abs() < 1e-6);
        assert!((area(&out, 1) - 4.0).abs() < 1e-6);
    }

    #[test]
    fn a_group_field_dissolves_overlapping_solids_without_double_counting() {
        // Two boxes overlapping by 1x2: the dissolved footprint must be
        // 4 + 4 - 2 = 6, not 8. This is what the union buys over summing.
        let mut l = Layer::new("in");
        l.geom_type = Some(GeometryType::MultiPolygon);
        l.add_field(FieldDef::new("blk", FieldType::Text));
        for g in [
            box_mesh([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]),
            box_mesh([1.0, 0.0, 0.0], [3.0, 2.0, 2.0]),
        ] {
            l.add_feature(Some(g), &[("blk", FieldValue::Text("a".into()))])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        let path = memory_store::make_vector_memory_path(&id);

        let (out, res) = run(json!({"input": path, "group_field": "blk"}));
        assert_eq!(res.outputs["footprint_count"], json!(1));
        assert!((area(&out, 0) - 6.0).abs() < 1e-6, "got {}", area(&out, 0));
    }

    #[test]
    fn a_courtyard_survives_rather_than_being_hulled_over() {
        // Four bars forming a hollow square ring. A convex hull would report
        // 9; the union must report 9 - 1 = 8, keeping the courtyard.
        let bars = vec![
            box_mesh([0.0, 0.0, 0.0], [3.0, 1.0, 1.0]),
            box_mesh([0.0, 2.0, 0.0], [3.0, 3.0, 1.0]),
            box_mesh([0.0, 1.0, 0.0], [1.0, 2.0, 1.0]),
            box_mesh([2.0, 1.0, 0.0], [3.0, 2.0, 1.0]),
        ];
        let mut l = Layer::new("in");
        l.geom_type = Some(GeometryType::MultiPolygon);
        l.add_field(FieldDef::new("blk", FieldType::Text));
        for g in bars {
            l.add_feature(Some(g), &[("blk", FieldValue::Text("ring".into()))])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        let path = memory_store::make_vector_memory_path(&id);

        let (out, _) = run(json!({"input": path, "group_field": "blk"}));
        assert!(
            (area(&out, 0) - 8.0).abs() < 1e-6,
            "courtyard was filled in: area {}",
            area(&out, 0)
        );
    }

    #[test]
    fn the_footprint_feeds_regularize_building_footprints() {
        // The composition this tool exists for: its output must be a 2D areal
        // geometry the cleanup pipeline accepts.
        let (out, _) = run(json!({
            "input": layer_of(vec![box_mesh([0.0, 0.0, 0.0], [4.0, 4.0, 4.0])]),
        }));
        let geom = out.iter().next().unwrap().geometry.clone().unwrap();
        assert!(matches!(geom, Geometry::MultiPolygon(_)));
        // And it is genuinely 2D: no Z survives the projection.
        if let Geometry::MultiPolygon(parts) = geom {
            for (ext, _) in &parts {
                assert!(ext.0.iter().all(|c| c.z.is_none() || c.z == Some(0.0)));
            }
        }
    }

    #[test]
    fn errors_when_nothing_is_projectable() {
        let args: ToolArgs = serde_json::from_value(json!({
            "input": layer_of(vec![Geometry::Point(Coord::xyz(0.0, 0.0, 0.0))]),
        }))
        .unwrap();
        assert!(MultipatchFootprintTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn rejects_bad_parameters() {
        let path = layer_of(vec![box_mesh([0.0; 3], [1.0, 1.0, 1.0])]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            MultipatchFootprintTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"input": path, "simplify_tolerance": -1})));

        // A missing group field is caught at run time against the real schema.
        let args: ToolArgs =
            serde_json::from_value(json!({"input": path, "group_field": "nope"})).unwrap();
        assert!(MultipatchFootprintTool.run(&args, &ctx()).is_err());
    }
}
