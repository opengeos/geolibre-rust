//! GeoLibre tool: orthogonality-preserving building footprint simplification.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Simplify Building* (Cartography).
//! This is the missing middle step of the repo's raster-to-vector cleanup
//! pipeline (`polygonize` → `regularize_building_footprints` →
//! `smooth_natural_features`). The two neighbouring tools solve different
//! problems:
//!
//! * `regularize_building_footprints` snaps edges to dominant angles, fixing
//!   *orientation*, but never reduces vertex count — a footprint traced from a
//!   raster keeps its stair-stepped vertex density after regularization.
//! * The bundled `simplify_features` is Douglas-Peucker / Visvalingam: it
//!   reduces vertices but is angle-agnostic, so it rounds off exactly the right
//!   angles that make a building read as a building.
//!
//! The removal metric here is **corner-protected**: each vertex is scored by
//! how far the boundary would move if it were dropped, then penalized by how
//! close its interior angle is to a right angle. A near-90° corner therefore
//! survives far past the tolerance that removes an intermediate vertex on a
//! straight run, which is what preserves building character.
//!
//! `minimum_area` drops footprints too small to render at the target scale,
//! optionally emitting a point in their place so the feature is not lost.

use std::collections::BTreeMap;

use geo::{Area, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct SimplifyBuildingTool;

impl Tool for SimplifyBuildingTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "simplify_building",
            display_name: "Simplify Building",
            summary: "Reduce building footprint vertices while protecting right-angle corners, with a minimum-area drop threshold, like ArcGIS Simplify Building.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Building footprint polygons.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Simplification tolerance in map units: the maximum distance the boundary may move.",
                    required: true,
                },
                ToolParamSpec {
                    name: "minimum_area",
                    description: "Footprints below this area are dropped (default 0, keep all).",
                    required: false,
                },
                ToolParamSpec {
                    name: "keep_collapsed_points",
                    description: "Emit a point at the centroid of each footprint dropped by minimum_area (default false).",
                    required: false,
                },
                ToolParamSpec {
                    name: "corner_tolerance",
                    description: "Degrees within which an interior angle counts as a protected right angle (default 20).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        let tol = parse_optional_f64(args, "tolerance")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'tolerance'".to_string())
        })?;
        if !tol.is_finite() || tol <= 0.0 {
            return Err(ToolError::Validation(
                "'tolerance' must be greater than 0".to_string(),
            ));
        }
        if let Some(a) = parse_optional_f64(args, "minimum_area")? {
            if !a.is_finite() || a < 0.0 {
                return Err(ToolError::Validation(
                    "'minimum_area' must be zero or greater".to_string(),
                ));
            }
        }
        if let Some(c) = parse_optional_f64(args, "corner_tolerance")? {
            if !c.is_finite() || !(0.0..90.0).contains(&c) {
                return Err(ToolError::Validation(
                    "'corner_tolerance' must be in [0, 90)".to_string(),
                ));
            }
        }
        parse_optional_bool(args, "keep_collapsed_points")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let tolerance = parse_optional_f64(args, "tolerance")?.unwrap();
        let min_area = parse_optional_f64(args, "minimum_area")?.unwrap_or(0.0);
        let keep_points = parse_optional_bool(args, "keep_collapsed_points")?.unwrap_or(false);
        let corner_tol = parse_optional_f64(args, "corner_tolerance")?.unwrap_or(20.0);

        let layer = load_input_layer(input)?;

        let mut out = Layer::new(layer.name.clone());
        out.crs = layer.crs.clone();
        // Every simplified footprint is written as a MultiPolygon, so declare
        // that rather than inheriting the input's (usually Polygon) type --
        // writers that key off geom_type would otherwise mislabel the layer.
        out.geom_type = Some(GeometryType::MultiPolygon);
        for fd in layer.schema.fields().iter() {
            out.add_field(fd.clone());
        }
        out.add_field(FieldDef::new("collapsed", FieldType::Integer));

        ctx.progress
            .info(&format!("simplifying {} footprint(s)", layer.len()));

        let mut vertices_before = 0usize;
        let mut vertices_after = 0usize;
        let mut dropped = 0usize;
        let mut kept = 0usize;

        for (fi, feat) in layer.features.iter().enumerate() {
            let Some(geom) = feat.geometry.as_ref() else {
                continue;
            };
            let Some(mp) = to_multipolygon(geom) else {
                // Non-areal geometry passes straight through, untouched.
                let mut fields = base_fields(&layer, feat);
                fields.push(("collapsed".to_string(), FieldValue::Integer(0)));
                let refs: Vec<(&str, FieldValue)> = fields
                    .iter()
                    .map(|(k, v)| (k.as_str(), v.clone()))
                    .collect();
                out.add_feature(Some(geom.clone()), &refs)
                    .map_err(|e| ToolError::Execution(format!("failed writing feature: {e}")))?;
                kept += 1;
                continue;
            };

            vertices_before += count_vertices(geom);

            // Simplify every ring with the corner-protected metric.
            let simplified = MultiPolygon(
                mp.0.iter()
                    .map(|p| {
                        Polygon::new(
                            simplify_ring(p.exterior(), tolerance, corner_tol),
                            p.interiors()
                                .iter()
                                .map(|r| simplify_ring(r, tolerance, corner_tol))
                                .collect(),
                        )
                    })
                    .collect(),
            );

            let area = simplified.unsigned_area();
            if min_area > 0.0 && area < min_area {
                dropped += 1;
                if keep_points {
                    let (cx, cy) = centroid_of(&mp);
                    let mut fields = base_fields(&layer, feat);
                    fields.push(("collapsed".to_string(), FieldValue::Integer(1)));
                    let refs: Vec<(&str, FieldValue)> = fields
                        .iter()
                        .map(|(k, v)| (k.as_str(), v.clone()))
                        .collect();
                    out.add_feature(Some(Geometry::Point(Coord::xy(cx, cy))), &refs)
                        .map_err(|e| {
                            ToolError::Execution(format!("failed writing collapsed point: {e}"))
                        })?;
                }
                continue;
            }

            let g = multipolygon_to_geometry(&simplified);
            vertices_after += count_vertices(&g);
            let mut fields = base_fields(&layer, feat);
            fields.push(("collapsed".to_string(), FieldValue::Integer(0)));
            let refs: Vec<(&str, FieldValue)> = fields
                .iter()
                .map(|(k, v)| (k.as_str(), v.clone()))
                .collect();
            out.add_feature(Some(g), &refs)
                .map_err(|e| ToolError::Execution(format!("failed writing feature: {e}")))?;
            kept += 1;
            ctx.progress
                .progress((fi as f64 + 1.0) / layer.len().max(1) as f64);
        }

        // Emitting points changes the geometry type, so widen it rather than
        // leaving a polygon declaration that no longer matches the contents.
        if keep_points && dropped > 0 {
            out.geom_type = Some(GeometryType::GeometryCollection);
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(kept));
        outputs.insert("dropped".to_string(), json!(dropped));
        outputs.insert("vertices_before".to_string(), json!(vertices_before));
        outputs.insert("vertices_after".to_string(), json!(vertices_after));
        outputs.insert(
            "vertices_removed".to_string(),
            json!(vertices_before.saturating_sub(vertices_after)),
        );
        Ok(ToolRunResult { outputs })
    }
}

fn base_fields(layer: &Layer, feat: &wbvector::Feature) -> Vec<(String, FieldValue)> {
    layer
        .schema
        .fields()
        .iter()
        .enumerate()
        .map(|(i, fd)| (fd.name.clone(), feat.attributes[i].clone()))
        .collect()
}

/// Greedy corner-protected vertex removal.
///
/// Repeatedly drops the vertex whose removal displaces the boundary least,
/// where the displacement is scaled up for vertices near a right angle so
/// corners survive. Stops when the cheapest remaining removal would exceed
/// `tolerance`, or the ring is down to a triangle.
fn simplify_ring(ring: &LineString, tolerance: f64, corner_tol_deg: f64) -> LineString {
    // Work on the open ring (no closing duplicate).
    let mut pts: Vec<GeoCoord<f64>> = ring.0.clone();
    if pts.len() >= 2 && pts.first() == pts.last() {
        pts.pop();
    }
    if pts.len() <= 3 {
        return close(pts);
    }

    loop {
        if pts.len() <= 3 {
            break;
        }
        let n = pts.len();
        let mut best: Option<(f64, usize)> = None;
        for i in 0..n {
            let prev = pts[(i + n - 1) % n];
            let cur = pts[i];
            let next = pts[(i + 1) % n];
            let d = perpendicular_distance(cur, prev, next);
            let cost = d * corner_penalty(prev, cur, next, corner_tol_deg);
            if best.is_none_or(|(bc, _)| cost < bc) {
                best = Some((cost, i));
            }
        }
        match best {
            Some((cost, i)) if cost <= tolerance => {
                pts.remove(i);
            }
            _ => break,
        }
    }
    close(pts)
}

/// Multiplier applied to a vertex's displacement cost. A vertex whose interior
/// angle is within `corner_tol_deg` of 90° (or 270°) is a building corner and
/// is made expensive to remove; a vertex on a straight run is cheap.
fn corner_penalty(
    prev: GeoCoord<f64>,
    cur: GeoCoord<f64>,
    next: GeoCoord<f64>,
    corner_tol_deg: f64,
) -> f64 {
    let a = (prev.y - cur.y).atan2(prev.x - cur.x);
    let b = (next.y - cur.y).atan2(next.x - cur.x);
    let mut angle = (b - a).abs().to_degrees();
    if angle > 180.0 {
        angle = 360.0 - angle;
    }
    // Distance from a right angle, in degrees.
    let from_right = (angle - 90.0).abs();
    if from_right <= corner_tol_deg {
        // Scale smoothly from a large penalty at exactly 90 down to 1 at the
        // edge of the tolerance band, so the protection is not a hard cliff.
        let t = from_right / corner_tol_deg.max(f64::EPSILON);
        1.0 + (CORNER_WEIGHT - 1.0) * (1.0 - t)
    } else {
        1.0
    }
}

/// How much more expensive an exact right-angle corner is than a straight-run
/// vertex. Large enough that corners survive tolerances that flatten runs.
const CORNER_WEIGHT: f64 = 1000.0;

fn perpendicular_distance(p: GeoCoord<f64>, a: GeoCoord<f64>, b: GeoCoord<f64>) -> f64 {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let len2 = dx * dx + dy * dy;
    if len2 <= f64::EPSILON {
        return ((p.x - a.x).powi(2) + (p.y - a.y).powi(2)).sqrt();
    }
    ((dx * (a.y - p.y) - (a.x - p.x) * dy).abs()) / len2.sqrt()
}

fn close(mut pts: Vec<GeoCoord<f64>>) -> LineString {
    if let Some(first) = pts.first().copied() {
        if pts.last() != Some(&first) {
            pts.push(first);
        }
    }
    LineString::new(pts)
}

fn centroid_of(mp: &MultiPolygon) -> (f64, f64) {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0usize;
    for p in &mp.0 {
        for c in p.exterior().0.iter() {
            sx += c.x;
            sy += c.y;
            n += 1;
        }
    }
    if n == 0 {
        (0.0, 0.0)
    } else {
        (sx / n as f64, sy / n as f64)
    }
}

fn count_vertices(geom: &Geometry) -> usize {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => exterior.coords().len() + interiors.iter().map(|r| r.coords().len()).sum::<usize>(),
        Geometry::MultiPolygon(parts) => parts
            .iter()
            .map(|(e, hs)| e.coords().len() + hs.iter().map(|r| r.coords().len()).sum::<usize>())
            .sum(),
        _ => 0,
    }
}

fn to_multipolygon(geom: &Geometry) -> Option<MultiPolygon> {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => Some(MultiPolygon(vec![rings_to_polygon(exterior, interiors)])),
        Geometry::MultiPolygon(parts) => Some(MultiPolygon(
            parts.iter().map(|(e, i)| rings_to_polygon(e, i)).collect(),
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
    if coords.len() >= 2 && coords.first().map(|c| (c.x, c.y)) == coords.last().map(|c| (c.x, c.y))
    {
        coords.pop();
    }
    Ring::new(coords)
}

// ── parameter parsing ────────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

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

fn parse_optional_bool(args: &ToolArgs, k: &str) -> Result<Option<bool>, ToolError> {
    match args.get(k) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{k}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{k}' must be a boolean"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::memory_store;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn poly_layer(rings: Vec<Vec<(f64, f64)>>) -> String {
        let mut l = Layer::new("b")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        for r in rings {
            l.add_feature(
                Some(Geometry::polygon(
                    r.into_iter().map(|(x, y)| Coord::xy(x, y)).collect(),
                    vec![],
                )),
                &[],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SimplifyBuildingTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn ring_of(layer: &Layer, i: usize) -> Vec<(f64, f64)> {
        match layer.features[i].geometry.as_ref().unwrap() {
            Geometry::MultiPolygon(parts) => {
                parts[0].0.coords().iter().map(|c| (c.x, c.y)).collect()
            }
            Geometry::Polygon { exterior, .. } => {
                exterior.coords().iter().map(|c| (c.x, c.y)).collect()
            }
            other => panic!("expected polygon, got {other:?}"),
        }
    }

    /// THE property: a stair-stepped edge collapses to a straight run, but the
    /// building's four right-angle corners survive.
    #[test]
    fn removes_stairsteps_but_keeps_corners() {
        // A 100x100 square whose bottom edge has been traced as small steps.
        let mut r = vec![(0.0, 0.0)];
        for i in 1..=10 {
            let x = i as f64 * 10.0;
            r.push((x, 0.4)); // tiny stair-step jitter
            r.push((x, 0.0));
        }
        r.push((100.0, 100.0));
        r.push((0.0, 100.0));
        let (out, layer) = run(json!({
            "input": poly_layer(vec![r]), "tolerance": 2.0
        }));

        let simplified = ring_of(&layer, 0);
        assert!(
            out.outputs["vertices_removed"].as_f64().unwrap() > 0.0,
            "stair-steps should be removed"
        );
        // The four corners of the square must still be present.
        for corner in [(0.0, 0.0), (100.0, 100.0), (0.0, 100.0)] {
            assert!(
                simplified
                    .iter()
                    .any(|(x, y)| (x - corner.0).abs() < 1e-6 && (y - corner.1).abs() < 1e-6),
                "corner {corner:?} was removed; ring is {simplified:?}"
            );
        }
    }

    /// A vertex mid-way along a straight edge is redundant and goes.
    #[test]
    fn removes_collinear_vertex() {
        let r = vec![
            (0.0, 0.0),
            (50.0, 0.0), // collinear on the bottom edge
            (100.0, 0.0),
            (100.0, 100.0),
            (0.0, 100.0),
        ];
        let (out, _) = run(json!({ "input": poly_layer(vec![r]), "tolerance": 1.0 }));
        assert!(out.outputs["vertices_removed"].as_f64().unwrap() >= 1.0);
    }

    /// A clean rectangle is already minimal and must be left alone.
    #[test]
    fn leaves_a_clean_rectangle_untouched() {
        let r = vec![(0.0, 0.0), (100.0, 0.0), (100.0, 50.0), (0.0, 50.0)];
        let (out, layer) = run(json!({ "input": poly_layer(vec![r]), "tolerance": 5.0 }));
        assert_eq!(out.outputs["vertices_removed"], json!(0));
        assert_eq!(ring_of(&layer, 0).len(), 4);
    }

    /// minimum_area drops small footprints.
    #[test]
    fn minimum_area_drops_small_buildings() {
        let big = vec![(0.0, 0.0), (100.0, 0.0), (100.0, 100.0), (0.0, 100.0)];
        let small = vec![(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)];
        let (out, layer) = run(json!({
            "input": poly_layer(vec![big, small]),
            "tolerance": 1.0, "minimum_area": 100.0
        }));
        assert_eq!(out.outputs["dropped"], json!(1));
        assert_eq!(layer.len(), 1);
    }

    /// keep_collapsed_points replaces a dropped footprint with a point.
    #[test]
    fn keep_collapsed_points_emits_a_point() {
        let small = vec![(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)];
        let (out, layer) = run(json!({
            "input": poly_layer(vec![small]),
            "tolerance": 1.0, "minimum_area": 100.0,
            "keep_collapsed_points": true
        }));
        assert_eq!(out.outputs["dropped"], json!(1));
        assert_eq!(layer.len(), 1);
        assert!(matches!(
            layer.features[0].geometry.as_ref().unwrap(),
            Geometry::Point(_)
        ));
        let ci = layer.schema.field_index("collapsed").unwrap();
        assert_eq!(layer.features[0].attributes[ci].as_f64(), Some(1.0));
    }

    /// A ring can never be reduced below a triangle.
    #[test]
    fn never_degenerates_below_a_triangle() {
        let r = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        let (_, layer) = run(json!({
            "input": poly_layer(vec![r]), "tolerance": 1e9
        }));
        assert!(ring_of(&layer, 0).len() >= 3);
    }

    /// Non-polygon geometry passes through untouched.
    #[test]
    fn passes_through_non_polygons() {
        let mut l = Layer::new("x")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_feature(Some(Geometry::Point(Coord::xy(5.0, 5.0))), &[])
            .unwrap();
        let id = memory_store::put_vector(l);
        let (_, layer) = run(json!({
            "input": memory_store::make_vector_memory_path(&id), "tolerance": 1.0
        }));
        assert!(matches!(
            layer.features[0].geometry.as_ref().unwrap(),
            Geometry::Point(_)
        ));
    }

    #[test]
    fn rejects_bad_parameters() {
        let p = poly_layer(vec![vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]]);
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SimplifyBuildingTool.validate(&args).is_err()
        };
        assert!(bad(json!({ "input": p })));
        assert!(bad(json!({ "input": p, "tolerance": 0 })));
        assert!(bad(
            json!({ "input": p, "tolerance": 1, "minimum_area": -5 })
        ));
        assert!(bad(
            json!({ "input": p, "tolerance": 1, "corner_tolerance": 95 })
        ));
        assert!(bad(json!({ "tolerance": 1 })));
    }

    /// The corner penalty itself: a right angle costs far more than a straight run.
    #[test]
    fn corner_penalty_protects_right_angles() {
        let right = corner_penalty(
            GeoCoord { x: -1.0, y: 0.0 },
            GeoCoord { x: 0.0, y: 0.0 },
            GeoCoord { x: 0.0, y: 1.0 },
            20.0,
        );
        let straight = corner_penalty(
            GeoCoord { x: -1.0, y: 0.0 },
            GeoCoord { x: 0.0, y: 0.0 },
            GeoCoord { x: 1.0, y: 0.0 },
            20.0,
        );
        assert!(
            right > straight * 100.0,
            "right={right} straight={straight}"
        );
        assert!((straight - 1.0).abs() < 1e-9);
    }
}
