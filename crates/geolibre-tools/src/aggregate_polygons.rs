//! GeoLibre tool: aggregate nearby polygons into larger polygons.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Aggregate Polygons* (Cartography):
//! combine polygons that lie within an aggregation distance of one another into
//! single larger polygons, for smaller-scale display or for merging fragmented
//! patches of one class (e.g. regularized building footprints, or land-cover
//! slivers of the same type). It pairs with `regularize_building_footprints`
//! (regularize, then aggregate) and complements `delineate_built_up_areas`
//! (which is building-count-aware and produces settlement extents); this tool is
//! the general-purpose "merge polygons within a distance" primitive.
//!
//! The core is a morphological *closing*: dilate every polygon by
//! `aggregation_distance / 2`, union, then erode back by the same radius, so
//! polygons whose gap is at most `aggregation_distance` fuse while the outer
//! boundary returns to roughly the input extent (with rounded joins). The
//! original polygons are unioned back so aggregation only ever adds area.
//!
//! Options mirror the ArcGIS tool:
//!
//! - `min_area` drops aggregated polygons smaller than the threshold.
//! - `min_hole_size` fills interior holes smaller than the threshold (the gaps
//!   left between aggregated polygons; larger courtyards are kept).
//! - `barrier` is an optional line/polygon layer that aggregation may not cross:
//!   any bridge that would span a barrier is severed, so polygons on opposite
//!   sides of a road or river are not merged. Existing input polygons are never
//!   split — a barrier only blocks *new* aggregation.
//!
//! Output is a new polygon layer; each feature carries `part_count` (how many
//! input polygons it aggregates) and `area`. Non-polygon input features are
//! ignored.
//!
//! Scope for v1: the ArcGIS `orthogonal` option (preserving right angles for
//! buildings — regularize first with `regularize_building_footprints`) and the
//! one-to-many source→aggregate link table are not implemented.

use std::collections::BTreeMap;

use geo::{
    Area, BooleanOps, Buffer, Contains, Coord as GeoCoord, LineString, MultiPolygon, Point, Polygon,
};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct AggregatePolygonsTool;

impl Tool for AggregatePolygonsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "aggregate_polygons",
            display_name: "Aggregate Polygons",
            summary: "Combine polygons within an aggregation distance into larger polygons (morphological closing), with minimum-area and minimum-hole-size filters and an optional barrier layer aggregation may not cross.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon vector layer, format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "aggregation_distance",
                    description: "Polygons within this distance (CRS units) of each other are aggregated. Sets the morphological-closing radius to half this value.",
                    required: true,
                },
                ToolParamSpec {
                    name: "min_area",
                    description: "Drop aggregated polygons whose area (CRS units squared) is below this value. Default: keep all.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_hole_size",
                    description: "Fill interior holes whose area (CRS units squared) is below this value. Default: keep all holes.",
                    required: false,
                },
                ToolParamSpec {
                    name: "barrier",
                    description: "Optional barrier vector layer (lines or polygons) that aggregation may not cross; bridges spanning a barrier are severed. Existing input polygons are never split.",
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
        let layer_name = layer.name.clone();
        let layer_crs = layer.crs.clone();
        let input_count = layer.len();

        // Collect input polygons as `geo` polygons; keep each one's
        // representative point so an aggregate can count its source polygons.
        let mut polys: Vec<Polygon> = Vec::new();
        let mut reps: Vec<Point> = Vec::new();
        for feature in &layer.features {
            let Some(mp) = feature.geometry.as_ref().and_then(to_multipolygon) else {
                continue;
            };
            for poly in &mp.0 {
                if let Some(p) = representative_point(poly) {
                    reps.push(p);
                }
            }
            polys.extend(mp.0);
        }
        let part_total = polys.len();

        let radius = prm.aggregation_distance / 2.0;
        ctx.progress.info(&format!(
            "{input_count} feature(s): aggregating {part_total} polygon(s) (radius {radius:.4})"
        ));

        // Morphological closing, then restore the originals so aggregation only
        // ever adds area (closing is extensive; this guards the discrete buffer).
        let originals = MultiPolygon(polys);
        let mut closed = originals.buffer(radius).buffer(-radius).union(&originals);

        // A barrier severs any bridge that spans it: subtract a thin buffer of
        // the barrier, then re-union the originals so no input polygon is split.
        if let Some(barrier_path) = prm.barrier.as_deref() {
            let barrier = load_barrier(barrier_path, radius)?;
            if !barrier.0.is_empty() {
                closed = closed.difference(&barrier).union(&originals);
            }
        }

        let mut out_layer = Layer::new(layer_name);
        out_layer.crs = layer_crs;
        out_layer.add_field(FieldDef::new("part_count", FieldType::Integer));
        out_layer.add_field(FieldDef::new("area", FieldType::Float));
        out_layer.geom_type = Some(GeometryType::Polygon);

        let mut kept = 0usize;
        let mut dropped = 0usize;
        for part in closed.0 {
            let part = match prm.min_hole_size {
                Some(min) => fill_small_holes(part, min),
                None => part,
            };
            let area = part.unsigned_area();
            if prm.min_area.is_some_and(|m| area < m) {
                dropped += 1;
                continue;
            }
            let count = reps.iter().filter(|p| part.contains(*p)).count();
            out_layer
                .add_feature(
                    Some(polygon_to_geometry(&part)),
                    &[
                        ("part_count", FieldValue::Integer(count as i64)),
                        ("area", FieldValue::Float(area)),
                    ],
                )
                .map_err(|e| {
                    ToolError::Execution(format!("failed building output feature: {e}"))
                })?;
            kept += 1;
        }

        ctx.progress.info(&format!(
            "aggregated into {kept} polygon(s), dropped {dropped} below min_area"
        ));

        let feature_count = out_layer.len();
        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_count".to_string(), json!(input_count));
        outputs.insert("part_count".to_string(), json!(part_total));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("aggregate_count".to_string(), json!(kept));
        outputs.insert("dropped_count".to_string(), json!(dropped));
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    aggregation_distance: f64,
    min_area: Option<f64>,
    min_hole_size: Option<f64>,
    barrier: Option<String>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let aggregation_distance =
        parse_optional_f64(args, "aggregation_distance")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'aggregation_distance'".to_string())
        })?;
    if !(aggregation_distance > 0.0 && aggregation_distance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'aggregation_distance' must be a positive number".to_string(),
        ));
    }
    let min_area = parse_positive_opt(args, "min_area")?;
    let min_hole_size = parse_positive_opt(args, "min_hole_size")?;
    let barrier = parse_optional_str(args, "barrier")?.map(str::to_string);
    Ok(Params {
        aggregation_distance,
        min_area,
        min_hole_size,
        barrier,
    })
}

fn parse_positive_opt(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    let v = parse_optional_f64(args, key)?;
    if let Some(m) = v {
        if !(m > 0.0 && m.is_finite()) {
            return Err(ToolError::Validation(format!(
                "parameter '{key}' must be a positive number"
            )));
        }
    }
    Ok(v)
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

// ── Geometry helpers ──────────────────────────────────────────────────────────

/// Loads a barrier layer and returns a thin buffered `MultiPolygon` to subtract
/// from the aggregation. Lines gain a hairline width so a `difference` severs a
/// crossed bridge; polygon barriers are used directly. The width is a tiny
/// fraction of the closing radius so it does not erode meaningful area.
fn load_barrier(path: &str, radius: f64) -> Result<MultiPolygon, ToolError> {
    let layer = load_input_layer(path)?;
    let width = (radius * 1e-3).max(1e-6);
    let mut lines: Vec<LineString> = Vec::new();
    let mut polys: Vec<Polygon> = Vec::new();
    for feature in &layer.features {
        match feature.geometry.as_ref() {
            Some(Geometry::LineString(coords)) => lines.push(coords_to_linestring(coords)),
            Some(Geometry::MultiLineString(parts)) => {
                lines.extend(parts.iter().map(|c| coords_to_linestring(c)))
            }
            Some(g) => {
                if let Some(mp) = to_multipolygon(g) {
                    polys.extend(mp.0);
                }
            }
            None => {}
        }
    }
    let mut barrier = MultiPolygon(polys);
    if !lines.is_empty() {
        let line_buf = geo::MultiLineString(lines).buffer(width);
        barrier = barrier.union(&line_buf);
    }
    Ok(barrier)
}

fn coords_to_linestring(coords: &[Coord]) -> LineString {
    LineString::new(coords.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect())
}

/// A point guaranteed to lie inside `poly` (its centroid if that is interior,
/// otherwise a representative point), used to attribute a source polygon to the
/// aggregate that contains it.
fn representative_point(poly: &Polygon) -> Option<Point> {
    use geo::{Centroid, InteriorPoint};
    poly.centroid()
        .filter(|c| poly.contains(c))
        .or_else(|| poly.interior_point())
}

/// Drops interior rings (holes) whose area is below `min_hole_size`.
fn fill_small_holes(poly: Polygon, min_hole_size: f64) -> Polygon {
    let (exterior, interiors) = poly.into_inner();
    let kept: Vec<LineString> = interiors
        .into_iter()
        .filter(|ring| Polygon::new(ring.clone(), vec![]).unsigned_area() >= min_hole_size)
        .collect();
    Polygon::new(exterior, kept)
}

/// Converts a polygonal `wbvector` geometry to a `geo` `MultiPolygon`. Returns
/// `None` for non-polygon geometries (which the tool ignores).
fn to_multipolygon(geom: &Geometry) -> Option<MultiPolygon> {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => Some(MultiPolygon(vec![rings_to_polygon(exterior, interiors)])),
        Geometry::MultiPolygon(parts) => Some(MultiPolygon(
            parts
                .iter()
                .map(|(ext, ints)| rings_to_polygon(ext, ints))
                .collect(),
        )),
        _ => None,
    }
}

fn rings_to_polygon(exterior: &Ring, interiors: &[Ring]) -> Polygon {
    Polygon::new(
        ring_to_linestring(exterior),
        interiors.iter().map(ring_to_linestring).collect(),
    )
}

fn ring_to_linestring(ring: &Ring) -> LineString {
    LineString::new(
        ring.coords()
            .iter()
            .map(|c| GeoCoord { x: c.x, y: c.y })
            .collect(),
    )
}

/// Converts a single `geo` `Polygon` to a `wbvector` `Geometry::Polygon`.
fn polygon_to_geometry(poly: &Polygon) -> Geometry {
    Geometry::Polygon {
        exterior: linestring_to_ring(poly.exterior()),
        interiors: poly.interiors().iter().map(linestring_to_ring).collect(),
    }
}

fn linestring_to_ring(ls: &LineString) -> Ring {
    let mut coords: Vec<Coord> = ls.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
    if coords.len() >= 2 && coords.first() == coords.last() {
        coords.pop();
    }
    Ring::new(coords)
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn square(x0: f64, y0: f64, s: f64) -> Geometry {
        rect(x0, y0, x0 + s, y0 + s)
    }

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Geometry {
        Geometry::polygon(
            vec![
                Coord::xy(x0, y0),
                Coord::xy(x1, y0),
                Coord::xy(x1, y1),
                Coord::xy(x0, y1),
            ],
            vec![],
        )
    }

    fn layer_path(geoms: Vec<Geometry>) -> String {
        let mut layer = Layer::new("polys");
        for g in geoms {
            layer.add_feature(Some(g), &[]).unwrap();
        }
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = AggregatePolygonsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn part_count(layer: &Layer, idx: usize) -> i64 {
        match layer.features[idx]
            .get(&layer.schema, "part_count")
            .unwrap()
        {
            FieldValue::Integer(i) => *i,
            other => panic!("part_count should be integer, got {other:?}"),
        }
    }

    #[test]
    fn aggregates_close_polygons_and_isolates_distant() {
        let input = layer_path(vec![
            square(0.0, 0.0, 4.0),
            square(6.0, 0.0, 4.0),  // gap 2 -> aggregate with the first
            square(30.0, 0.0, 4.0), // far -> its own aggregate
        ]);
        let (out, layer) = run_tool(json!({ "input": input, "aggregation_distance": 6.0 }));
        assert_eq!(out.outputs["part_count"], json!(3));
        assert_eq!(out.outputs["feature_count"], json!(2));
        let mut counts: Vec<i64> = (0..layer.len()).map(|i| part_count(&layer, i)).collect();
        counts.sort_unstable();
        assert_eq!(counts, vec![1, 2]);
    }

    #[test]
    fn min_area_drops_small_aggregates() {
        let input = layer_path(vec![
            square(0.0, 0.0, 4.0),
            square(6.0, 0.0, 4.0),
            square(40.0, 0.0, 1.0), // tiny loner
        ]);
        let (out, _) =
            run_tool(json!({ "input": input, "aggregation_distance": 6.0, "min_area": 20.0 }));
        assert_eq!(out.outputs["feature_count"], json!(1));
        assert_eq!(out.outputs["dropped_count"], json!(1));
    }

    #[test]
    fn min_hole_size_fills_small_holes() {
        // Four bars overlapping at the corners form a frame around a 14x14
        // courtyard (area 196); min_hole_size above that should fill it.
        let input = layer_path(vec![
            rect(0.0, 0.0, 30.0, 8.0),   // bottom
            rect(0.0, 22.0, 30.0, 30.0), // top
            rect(0.0, 0.0, 8.0, 30.0),   // left
            rect(22.0, 0.0, 30.0, 30.0), // right
        ]);
        let hole_count = |layer: &Layer| {
            layer
                .features
                .iter()
                .map(|f| match f.geometry.as_ref().unwrap() {
                    Geometry::Polygon { interiors, .. } => interiors.len(),
                    _ => 0,
                })
                .sum::<usize>()
        };
        let (_, kept) = run_tool(json!({ "input": input.clone(), "aggregation_distance": 2.0 }));
        let (_, filled) = run_tool(
            json!({ "input": input, "aggregation_distance": 2.0, "min_hole_size": 400.0 }),
        );
        assert!(
            hole_count(&kept) >= 1,
            "expected a courtyard hole before filling"
        );
        assert_eq!(hole_count(&filled), 0, "courtyard hole should be filled");
    }

    #[test]
    fn barrier_prevents_aggregation_across_it() {
        // Two blocks a short gap apart would aggregate; a barrier line running
        // through the gap keeps them separate.
        let input = layer_path(vec![square(0.0, 0.0, 4.0), square(6.0, 0.0, 4.0)]);
        let mut barrier_layer = Layer::new("barrier");
        barrier_layer
            .add_feature(
                Some(Geometry::line_string(vec![
                    Coord::xy(5.0, -5.0),
                    Coord::xy(5.0, 15.0),
                ])),
                &[],
            )
            .unwrap();
        let bid = memory_store::put_vector(barrier_layer);
        let barrier = memory_store::make_vector_memory_path(&bid);

        let (with_barrier, _) = run_tool(
            json!({ "input": input.clone(), "aggregation_distance": 6.0, "barrier": barrier }),
        );
        let (no_barrier, _) = run_tool(json!({ "input": input, "aggregation_distance": 6.0 }));
        assert_eq!(
            no_barrier.outputs["feature_count"],
            json!(1),
            "should merge without a barrier"
        );
        assert_eq!(
            with_barrier.outputs["feature_count"],
            json!(2),
            "barrier must block the aggregation"
        );
    }

    #[test]
    fn ignores_non_polygons() {
        let mut layer = Layer::new("mixed");
        layer
            .add_feature(Some(Geometry::point(1.0, 1.0)), &[])
            .unwrap();
        layer.add_feature(Some(square(0.0, 0.0, 4.0)), &[]).unwrap();
        layer.add_feature(Some(square(6.0, 0.0, 4.0)), &[]).unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);
        let (out, layer) = run_tool(json!({ "input": input, "aggregation_distance": 6.0 }));
        assert_eq!(out.outputs["part_count"], json!(2));
        assert_eq!(out.outputs["feature_count"], json!(1));
        assert!(layer
            .features
            .iter()
            .all(|f| matches!(f.geometry, Some(Geometry::Polygon { .. }))));
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = AggregatePolygonsTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "input": "x.geojson" })).is_err(),
            "missing distance"
        );
        assert!(bad(json!({ "input": "x.geojson", "aggregation_distance": 0 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "aggregation_distance": -1 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "aggregation_distance": "wide" })).is_err());
        assert!(
            bad(json!({ "input": "x.geojson", "aggregation_distance": 5, "min_area": -1 }))
                .is_err()
        );
        assert!(bad(
            json!({ "input": "x.geojson", "aggregation_distance": "5", "min_hole_size": "10" })
        )
        .is_ok());
    }
}
