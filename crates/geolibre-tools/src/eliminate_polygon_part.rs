//! GeoLibre tool: eliminate small polygon parts and holes below a size threshold.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Eliminate Polygon Part* (Data
//! Management). Where the bundled whitebox `remove_polygon_holes` strips **all**
//! interior rings unconditionally, and GeoLibre's `eliminate_polygons` merges
//! whole sliver polygons into neighbours, this tool performs the finer cleanup:
//! it drops the *small* holes and (optionally) small outer parts of each
//! (multi)polygon while leaving the geometry otherwise intact and preserving
//! every attribute. It is the natural speckle-removal step after `polygonize`
//! turns a raster classification into vectors.
//!
//! Two selection modes, mirroring ArcGIS:
//!
//! - `condition = AREA` — a part or hole is removed when its area (in CRS
//!   units²) is below `min_area`.
//! - `condition = PERCENT` — removed when its area, as a percentage of the
//!   feature's total outer area (sum of the exterior-ring areas of all parts),
//!   is below `percentage`.
//!
//! `part_option` controls what is eligible:
//!
//! - `CONTAINED_ONLY` (default) — only interior rings (holes) are removed;
//!   exterior parts are always kept.
//! - `ANY` — small exterior parts are removed too, except the single largest
//!   part, which is always retained so the feature never disappears.
//!
//! The work is pure ring-area (shoelace) arithmetic plus ring reassembly — no
//! overlay engine required. Non-polygon features pass through untouched.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Feature, Geometry, GeometryType, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct EliminatePolygonPartTool;

impl Tool for EliminatePolygonPartTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "eliminate_polygon_part",
            display_name: "Eliminate Polygon Part",
            summary: "Remove interior holes (and optionally small outer parts) of polygons whose area is below an absolute-area or percentage threshold, preserving attributes.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon vector file path, format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "condition",
                    description: "Size test: 'AREA' (default) compares against 'min_area'; 'PERCENT' compares the part/hole area as a percentage of the feature's total outer area against 'percentage'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_area",
                    description: "AREA condition: remove parts/holes whose area (in CRS units squared) is below this value. Required when condition = AREA.",
                    required: false,
                },
                ToolParamSpec {
                    name: "percentage",
                    description: "PERCENT condition: remove parts/holes whose area, as a percentage (0-100) of the feature's total outer area, is below this value. Required when condition = PERCENT.",
                    required: false,
                },
                ToolParamSpec {
                    name: "part_option",
                    description: "What is eligible for removal: 'CONTAINED_ONLY' (default) removes only interior rings (holes); 'ANY' also removes small exterior parts except the single largest part.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args
            .get("input")
            .and_then(Value::as_str)
            .map(str::trim)
            .unwrap_or("")
            .is_empty()
        {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args
            .get("input")
            .and_then(Value::as_str)
            .filter(|s| !s.trim().is_empty())
            .ok_or_else(|| {
                ToolError::Validation("missing required parameter 'input'".to_string())
            })?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let schema = layer.schema.clone();
        let layer_name = layer.name.clone();
        let layer_crs = layer.crs.clone();
        let input_count = layer.len();

        let mut holes_removed = 0usize;
        let mut parts_removed = 0usize;
        let mut has_multipolygon = false;

        let mut out_features: Vec<Feature> = Vec::with_capacity(input_count);
        for feature in layer.features.into_iter() {
            let mut feature = feature;
            if let Some(geom) = feature.geometry.as_ref() {
                if let Some(parts) = to_parts(geom) {
                    let (kept, hr, pr) = clean_parts(parts, &prm);
                    holes_removed += hr;
                    parts_removed += pr;
                    let new_geom = parts_to_geometry(kept);
                    has_multipolygon |= matches!(new_geom, Geometry::MultiPolygon(_));
                    feature.geometry = Some(new_geom);
                }
            }
            feature.fid = out_features.len() as u64;
            out_features.push(feature);
        }

        ctx.progress.info(&format!(
            "{input_count} feature(s): removed {holes_removed} hole(s) and {parts_removed} outer part(s)",
        ));

        let mut out_layer = wbvector::Layer::new(layer_name);
        out_layer.schema = schema;
        out_layer.crs = layer_crs;
        out_layer.features = out_features;
        out_layer.geom_type = if has_multipolygon {
            Some(GeometryType::MultiPolygon)
        } else {
            Some(GeometryType::Polygon)
        };

        let feature_count = out_layer.len();
        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_count".to_string(), json!(input_count));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("holes_removed".to_string(), json!(holes_removed));
        outputs.insert("parts_removed".to_string(), json!(parts_removed));
        Ok(ToolRunResult { outputs })
    }
}

// ── Core cleanup ────────────────────────────────────────────────────────────

/// One polygon part: an exterior ring and its interior rings (holes).
type Part = (Ring, Vec<Ring>);

/// Filters holes and (in `ANY` mode) small outer parts from a feature's parts.
/// Returns the surviving parts plus counts of removed holes and outer parts.
fn clean_parts(parts: Vec<Part>, prm: &Params) -> (Vec<Part>, usize, usize) {
    // Reference for the PERCENT condition: the feature's total outer area, i.e.
    // the sum of every part's exterior-ring area.
    let ext_areas: Vec<f64> = parts.iter().map(|(ext, _)| ring_area(ext)).collect();
    let total_outer: f64 = ext_areas.iter().sum();

    let mut holes_removed = 0usize;
    let mut parts_removed = 0usize;

    // In ANY mode, protect the single largest part (by exterior area) so a
    // feature can never be reduced to nothing, matching ArcGIS.
    let largest_idx = ext_areas
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i);

    let mut kept: Vec<Part> = Vec::with_capacity(parts.len());
    for (idx, (ext, interiors)) in parts.into_iter().enumerate() {
        let ext_area = ext_areas[idx];

        // Drop the whole outer part when eligible and below threshold. Its holes
        // vanish with it but are not counted as hole removals.
        if prm.part_option == PartOption::Any
            && Some(idx) != largest_idx
            && prm.below_threshold(ext_area, total_outer)
        {
            parts_removed += 1;
            continue;
        }

        // Filter this part's holes.
        let mut new_interiors: Vec<Ring> = Vec::with_capacity(interiors.len());
        for hole in interiors {
            let hole_area = ring_area(&hole);
            if prm.below_threshold(hole_area, total_outer) {
                holes_removed += 1;
            } else {
                new_interiors.push(hole);
            }
        }
        kept.push((ext, new_interiors));
    }

    (kept, holes_removed, parts_removed)
}

// ── Parameters ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Condition {
    Area,
    Percent,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum PartOption {
    ContainedOnly,
    Any,
}

struct Params {
    condition: Condition,
    min_area: f64,
    percentage: f64,
    part_option: PartOption,
}

impl Params {
    /// True when a ring/part of area `area` falls below the active threshold.
    /// `total_outer` is the feature's total outer area (for PERCENT mode).
    fn below_threshold(&self, area: f64, total_outer: f64) -> bool {
        match self.condition {
            Condition::Area => area < self.min_area,
            Condition::Percent => {
                if total_outer <= 0.0 {
                    // Degenerate feature: nothing to compare against; keep it.
                    false
                } else {
                    (area / total_outer) * 100.0 < self.percentage
                }
            }
        }
    }
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let condition = match parse_optional_str(args, "condition")?
        .map(|s| s.trim().to_ascii_uppercase())
        .as_deref()
    {
        None | Some("AREA") => Condition::Area,
        Some("PERCENT") => Condition::Percent,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown condition '{other}' (expected AREA or PERCENT)"
            )))
        }
    };

    let min_area = parse_optional_f64(args, "min_area")?;
    let percentage = parse_optional_f64(args, "percentage")?;

    match condition {
        Condition::Area => {
            let m = min_area.ok_or_else(|| {
                ToolError::Validation(
                    "condition AREA requires 'min_area' (a positive number)".to_string(),
                )
            })?;
            if !(m > 0.0 && m.is_finite()) {
                return Err(ToolError::Validation(
                    "parameter 'min_area' must be a positive number".to_string(),
                ));
            }
        }
        Condition::Percent => {
            let p = percentage.ok_or_else(|| {
                ToolError::Validation(
                    "condition PERCENT requires 'percentage' (a number in 0..100)".to_string(),
                )
            })?;
            if !(p > 0.0 && p <= 100.0) {
                return Err(ToolError::Validation(
                    "parameter 'percentage' must be greater than 0 and at most 100".to_string(),
                ));
            }
        }
    }

    let part_option = match parse_optional_str(args, "part_option")?
        .map(|s| s.trim().to_ascii_uppercase())
        .as_deref()
    {
        None | Some("CONTAINED_ONLY") => PartOption::ContainedOnly,
        // Accept ALL as a friendly alias for ArcGIS's ANY.
        Some("ANY") | Some("ALL") => PartOption::Any,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown part_option '{other}' (expected CONTAINED_ONLY or ANY)"
            )))
        }
    };

    Ok(Params {
        condition,
        min_area: min_area.unwrap_or(0.0),
        percentage: percentage.unwrap_or(0.0),
        part_option,
    })
}

/// Parses an optional numeric parameter, accepting a JSON number or a numeric
/// string (host UIs often post form values as strings).
fn parse_optional_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be a number"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a number"
        ))),
    }
}

// ── Geometry helpers ────────────────────────────────────────────────────────

/// One polygon part reassembly type shared with the reader/writer helpers.
///
/// Decomposes a polygonal `wbvector` geometry into its parts (exterior + holes).
/// Returns `None` for non-polygon geometries (passed through untouched).
fn to_parts(geom: &Geometry) -> Option<Vec<Part>> {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => Some(vec![(exterior.clone(), interiors.clone())]),
        Geometry::MultiPolygon(parts) => Some(
            parts
                .iter()
                .map(|(ext, ints)| (ext.clone(), ints.clone()))
                .collect(),
        ),
        _ => None,
    }
}

/// Reassembles surviving parts into a `Polygon` (single part) or `MultiPolygon`.
fn parts_to_geometry(mut parts: Vec<Part>) -> Geometry {
    if parts.len() == 1 {
        let (exterior, interiors) = parts.pop().unwrap();
        Geometry::Polygon {
            exterior,
            interiors,
        }
    } else {
        Geometry::MultiPolygon(parts)
    }
}

/// Unsigned area of a ring via the shoelace formula. `Ring` stores the closing
/// vertex implicitly, so the last→first edge is added explicitly.
fn ring_area(ring: &Ring) -> f64 {
    let pts = ring.coords();
    let n = pts.len();
    if n < 3 {
        return 0.0;
    }
    let mut sum = 0.0;
    for i in 0..n {
        let a = &pts[i];
        let b = &pts[(i + 1) % n];
        sum += a.x * b.y - b.x * a.y;
    }
    (sum * 0.5).abs()
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Coord, FieldDef, FieldType, FieldValue, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Coord> {
        vec![
            Coord::xy(x0, y0),
            Coord::xy(x1, y0),
            Coord::xy(x1, y1),
            Coord::xy(x0, y1),
        ]
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = EliminatePolygonPartTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn geom_area(geom: &Geometry) -> f64 {
        // Net area (exterior minus holes) summed over parts.
        to_parts(geom)
            .map(|parts| {
                parts
                    .iter()
                    .map(|(ext, ints)| ring_area(ext) - ints.iter().map(ring_area).sum::<f64>())
                    .sum()
            })
            .unwrap_or(0.0)
    }

    /// A 20x20 square (area 400) with a tiny 1x1 hole (area 1) and a big 5x5
    /// hole (area 25). AREA=10 should drop the tiny hole only: net area rises
    /// from 400-1-25=374 to 400-25=375.
    #[test]
    fn drops_small_hole_keeps_big_hole() {
        let mut layer = Layer::new("shapes");
        layer.add_field(FieldDef::new("id", FieldType::Integer));
        let small_hole = rect(2.0, 2.0, 3.0, 3.0);
        let big_hole = rect(10.0, 10.0, 15.0, 15.0);
        layer
            .add_feature(
                Some(Geometry::polygon(
                    rect(0.0, 0.0, 20.0, 20.0),
                    vec![small_hole, big_hole],
                )),
                &[("id", 1i64.into())],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) =
            run_tool(json!({ "input": input, "condition": "AREA", "min_area": 10.0 }));
        assert_eq!(out.outputs["holes_removed"], json!(1));
        assert_eq!(out.outputs["parts_removed"], json!(0));
        assert_eq!(out.outputs["feature_count"], json!(1));
        // id preserved
        assert_eq!(
            layer.features[0].get(&layer.schema, "id").unwrap(),
            &FieldValue::Integer(1)
        );
        // net area 374 -> 375 (tiny hole filled)
        assert!((geom_area(layer.features[0].geometry.as_ref().unwrap()) - 375.0).abs() < 1e-6);
    }

    /// PERCENT mode: hole area / total outer area * 100 must be below the value.
    /// Outer 400; hole 1 -> 0.25%; hole 25 -> 6.25%. percentage=1 drops only the
    /// 0.25% hole.
    #[test]
    fn percent_condition_drops_relative_small_hole() {
        let mut layer = Layer::new("shapes");
        layer
            .add_feature(
                Some(Geometry::polygon(
                    rect(0.0, 0.0, 20.0, 20.0),
                    vec![rect(2.0, 2.0, 3.0, 3.0), rect(10.0, 10.0, 15.0, 15.0)],
                )),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, _) =
            run_tool(json!({ "input": input, "condition": "PERCENT", "percentage": 1.0 }));
        assert_eq!(out.outputs["holes_removed"], json!(1));
    }

    /// ANY mode drops a small outer part but never the largest part. A big 10x10
    /// part (area 100) and a tiny 1x1 part (area 1); AREA=10, part_option=ANY
    /// removes the tiny part, leaving a single-part polygon of area 100.
    #[test]
    fn any_mode_drops_small_outer_part_but_keeps_largest() {
        let mut layer = Layer::new("shapes");
        let big = (Ring::new(rect(0.0, 0.0, 10.0, 10.0)), vec![]);
        let tiny = (Ring::new(rect(50.0, 50.0, 51.0, 51.0)), vec![]);
        layer
            .add_feature(Some(Geometry::MultiPolygon(vec![big, tiny])), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({
            "input": input, "condition": "AREA", "min_area": 10.0, "part_option": "ANY"
        }));
        assert_eq!(out.outputs["parts_removed"], json!(1));
        assert_eq!(out.outputs["feature_count"], json!(1));
        assert!((geom_area(layer.features[0].geometry.as_ref().unwrap()) - 100.0).abs() < 1e-6);
    }

    /// CONTAINED_ONLY (default) never removes outer parts, even tiny ones.
    #[test]
    fn contained_only_leaves_outer_parts() {
        let mut layer = Layer::new("shapes");
        layer
            .add_feature(
                Some(Geometry::MultiPolygon(vec![
                    (Ring::new(rect(0.0, 0.0, 10.0, 10.0)), vec![]),
                    (Ring::new(rect(50.0, 50.0, 51.0, 51.0)), vec![]),
                ])),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, _) = run_tool(json!({ "input": input, "condition": "AREA", "min_area": 10.0 }));
        assert_eq!(out.outputs["parts_removed"], json!(0));
    }

    /// A polygon with only large holes is unchanged (pass-through property), and
    /// non-polygon features pass through untouched.
    #[test]
    fn non_matching_and_non_polygon_pass_through() {
        let mut layer = Layer::new("mixed");
        layer
            .add_feature(Some(Geometry::point(1.0, 2.0)), &[])
            .unwrap();
        // Polygon with a big hole (25) that stays under AREA=10.
        layer
            .add_feature(
                Some(Geometry::polygon(
                    rect(0.0, 0.0, 20.0, 20.0),
                    vec![rect(10.0, 10.0, 15.0, 15.0)],
                )),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) =
            run_tool(json!({ "input": input, "condition": "AREA", "min_area": 10.0 }));
        assert_eq!(out.outputs["holes_removed"], json!(0));
        assert_eq!(out.outputs["feature_count"], json!(2));
        assert!(layer
            .features
            .iter()
            .any(|f| matches!(f.geometry, Some(Geometry::Point(_)))));
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = EliminatePolygonPartTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input must fail");
        assert!(
            bad(json!({ "input": "x.geojson" })).is_err(),
            "AREA condition needs min_area"
        );
        assert!(bad(json!({ "input": "x.geojson", "min_area": 0 })).is_err());
        assert!(
            bad(json!({ "input": "x.geojson", "condition": "PERCENT" })).is_err(),
            "PERCENT needs percentage"
        );
        assert!(
            bad(json!({ "input": "x.geojson", "condition": "PERCENT", "percentage": 150 }))
                .is_err()
        );
        assert!(bad(json!({ "input": "x.geojson", "condition": "BOGUS", "min_area": 5 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "min_area": 5, "part_option": "wat" })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "min_area": 5 })).is_ok());
        assert!(
            bad(json!({ "input": "x.geojson", "min_area": "5.0" })).is_ok(),
            "numeric strings ok"
        );
        assert!(
            bad(json!({ "input": "x.geojson", "condition": "PERCENT", "percentage": "10" }))
                .is_ok()
        );
    }
}
