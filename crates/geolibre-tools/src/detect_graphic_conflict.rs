//! GeoLibre tool: find where the symbolized extents of two feature layers
//! graphically conflict.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Detect Graphic Conflict*
//! (Cartography). The repo already ships the *resolution* half of
//! cartographic conflict handling (`resolve_road_conflicts`,
//! `resolve_building_conflicts`), which move geometry to clear conflicts; this
//! tool is the QA/review counterpart that only *reports* them, leaving the
//! source data untouched. Nothing in the bundled whitebox suite reasons about
//! drawn symbol width — `count_overlapping_features` only sees raw geometry.
//!
//! ## Approach
//!
//! At a target scale a feature is drawn as a shape some *symbol width* wide
//! (a line's stroke width, a point's marker diameter, a polygon's outline
//! weight), not just its bare geometry. Two features graphically conflict when
//! those drawn extents overlap or crowd closer than an extra required
//! clearance. Each input geometry is buffered outward by half its symbol
//! width plus half `conflict_distance` using `geo`'s `Buffer` trait (the
//! `i_overlay` backend already used by `multiple_ring_buffer` and
//! `resolve_road_conflicts`'s sibling tools — pure Rust, no GEOS). Candidate
//! pairs are pruned with a bounding-box test (as `count_overlapping_features`
//! does), then the buffered polygons of a surviving pair are intersected with
//! `geo`'s `BooleanOps::intersection`; every non-empty result becomes one
//! output conflict polygon tagged with both source feature ids and the
//! overlap area. This buffers the *true* symbol footprint and intersects it
//! exactly, rather than approximating conflicts from centreline distance the
//! way `resolve_road_conflicts` does — real polygon buffering is available
//! here via `geo::Buffer`, so the "acceptable v1" distance-only fallback
//! described in the issue was not needed.
//!
//! `input` is always buffered by `symbol_width`. When `conflict` names a
//! second layer it is buffered by `conflict_symbol_width` (default: the same
//! as `symbol_width`) and every `input` × `conflict` pair is tested. When
//! `conflict` is omitted, `input` is checked against itself: every unordered
//! pair of its own features is tested once (`conflict_symbol_width` is
//! ignored in this mode — there is only one role to buffer).
//!
//! For a self-comparison of line features, two lines that are legitimately
//! joined end-to-end (a road continuing at a junction) always overlap right at
//! the shared vertex — that is normal connectivity, not a conflict. When
//! `line_connection_allowance > 0`, any pair of `LineString`/`MultiLineString`
//! features that share an endpoint (coincident within a tight snap tolerance)
//! has a disc of that radius, centred on each shared endpoint, subtracted from
//! their conflict polygon; if nothing survives the pair is dropped entirely.
//!
//! ## Scope for v1
//!
//! Symbol widths are fixed scalars (`symbol_width` / `conflict_symbol_width`),
//! not per-feature fields — unlike `resolve_road_conflicts`'s
//! `symbol_width_field`. `input` and `conflict` are assumed to share one CRS
//! and unit system (no on-the-fly reprojection). Shared-endpoint detection
//! uses exact (tightly-snapped) vertex coincidence, so line networks that are
//! not topologically snapped will not get connection-allowance relief at
//! near-but-not-exact junctions.

use std::collections::BTreeMap;

use geo::{
    Area, BooleanOps, BoundingRect, Buffer, Geometry as GeoGeometry, LineString, MultiLineString,
    MultiPoint, MultiPolygon, Point, Polygon,
};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Intersection/difference results smaller than this area (CRS units^2) are
/// treated as numerical slivers, not conflicts.
const SLIVER_AREA_EPS: f64 = 1e-9;
/// Coordinate tolerance for treating two line endpoints as the same vertex.
const ENDPOINT_SNAP_EPS: f64 = 1e-9;

pub struct DetectGraphicConflictTool;

impl Tool for DetectGraphicConflictTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "detect_graphic_conflict",
            display_name: "Detect Graphic Conflict",
            summary: "Find where the symbolized (drawn-width) extents of two feature layers — or one layer against itself — overlap or crowd within a minimum separation, emitting a conflict polygon per conflicting pair tagged with both source feature ids and the overlap area, like ArcGIS Detect Graphic Conflict.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (points, lines, or polygons).",
                    required: true,
                },
                ToolParamSpec {
                    name: "conflict",
                    description: "Second vector layer to check 'input' against. Omit to check 'input' against itself.",
                    required: false,
                },
                ToolParamSpec {
                    name: "symbol_width",
                    description: "Symbolized width/diameter of 'input' features, in map units.",
                    required: true,
                },
                ToolParamSpec {
                    name: "conflict_symbol_width",
                    description: "Symbolized width/diameter of 'conflict' features, in map units. Default: same as 'symbol_width'. Ignored when 'conflict' is omitted.",
                    required: false,
                },
                ToolParamSpec {
                    name: "conflict_distance",
                    description: "Extra required separation beyond touching, in map units. Default 0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "line_connection_allowance",
                    description: "Radius around a shared endpoint of two same-layer line features within which their overlap is treated as normal end-to-end connectivity, not a conflict. Default 0 (no allowance).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon vector path with 'input_fid', 'conflict_fid', and 'overlap_area' fields. If omitted, stored in memory.",
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

        let input_layer = load_input_layer(input)?;
        let half_input = prm.symbol_width / 2.0 + prm.conflict_distance / 2.0;
        let feats_a = build_features(&input_layer, half_input);

        let (feats_b, same_layer) = match prm.conflict.as_deref() {
            Some(path) => {
                let half_conflict = prm.conflict_symbol_width.unwrap_or(prm.symbol_width) / 2.0
                    + prm.conflict_distance / 2.0;
                let conflict_layer = load_input_layer(path)?;
                (build_features(&conflict_layer, half_conflict), false)
            }
            None => (feats_a.clone(), true),
        };

        if feats_a.is_empty() || feats_b.is_empty() {
            return Err(ToolError::Execution(
                "no usable geometry in input/conflict layer(s)".to_string(),
            ));
        }

        ctx.progress.info(&format!(
            "{} input feature(s), {} conflict feature(s){}",
            feats_a.len(),
            feats_b.len(),
            if same_layer { " (self-comparison)" } else { "" }
        ));

        let mut out = Layer::new("graphic_conflicts").with_geom_type(GeometryType::MultiPolygon);
        if let Some(epsg) = input_layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("input_fid", FieldType::Integer));
        out.add_field(FieldDef::new("conflict_fid", FieldType::Integer));
        out.add_field(FieldDef::new("overlap_area", FieldType::Float));

        let mut conflict_count = 0usize;
        let mut suppressed_count = 0usize;
        let mut total_overlap_area = 0.0f64;

        let mut emit = |a: &FeatureGeom,
                        b: &FeatureGeom,
                        out: &mut Layer|
         -> Result<(), ToolError> {
            let Some((geom, area)) = evaluate_pair(
                a,
                b,
                same_layer,
                prm.line_connection_allowance,
                &mut suppressed_count,
            ) else {
                return Ok(());
            };
            out.add_feature(
                Some(geom),
                &[
                    ("input_fid", (a.fid as i64).into()),
                    ("conflict_fid", (b.fid as i64).into()),
                    ("overlap_area", area.into()),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed writing conflict polygon: {e}")))?;
            conflict_count += 1;
            total_overlap_area += area;
            Ok(())
        };

        if same_layer {
            for i in 0..feats_a.len() {
                for j in (i + 1)..feats_a.len() {
                    emit(&feats_a[i], &feats_a[j], &mut out)?;
                }
            }
        } else {
            for a in &feats_a {
                for b in &feats_b {
                    emit(a, b, &mut out)?;
                }
            }
        }

        ctx.progress.info(&format!(
            "{conflict_count} conflict polygon(s); {suppressed_count} end-connection(s) suppressed"
        ));

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_feature_count".to_string(), json!(feats_a.len()));
        outputs.insert("conflict_feature_count".to_string(), json!(feats_b.len()));
        outputs.insert("conflict_polygon_count".to_string(), json!(conflict_count));
        outputs.insert(
            "suppressed_connections".to_string(),
            json!(suppressed_count),
        );
        outputs.insert("total_overlap_area".to_string(), json!(total_overlap_area));
        outputs.insert("same_layer".to_string(), json!(same_layer));
        Ok(ToolRunResult { outputs })
    }
}

// ── Feature geometry ────────────────────────────────────────────────────────

#[derive(Clone)]
struct FeatureGeom {
    fid: usize,
    buffered: MultiPolygon,
    bbox: [f64; 4],
    is_line: bool,
    endpoints: Vec<(f64, f64)>,
}

/// Converts every usable feature of `layer` into a symbol-buffered geometry
/// (buffered by `half`) plus its bounding box and, for line features, its
/// endpoints (used for connection-allowance suppression).
fn build_features(layer: &Layer, half: f64) -> Vec<FeatureGeom> {
    let mut out = Vec::new();
    for (fidx, f) in layer.features.iter().enumerate() {
        let Some(geom) = f.geometry.as_ref() else {
            continue;
        };
        let Some(g) = to_geo_geometry(geom) else {
            continue;
        };
        let is_line = matches!(
            g,
            GeoGeometry::LineString(_) | GeoGeometry::MultiLineString(_)
        );
        let endpoints = if is_line {
            line_endpoints(&g)
        } else {
            Vec::new()
        };
        let buffered = g.buffer(half);
        if buffered.0.is_empty() {
            continue;
        }
        let bbox = bbox(&buffered);
        out.push(FeatureGeom {
            fid: fidx,
            buffered,
            bbox,
            is_line,
            endpoints,
        });
    }
    out
}

/// Tests one candidate pair. Returns the conflict polygon and its area, or
/// `None` if the pair does not conflict (or the whole conflict is suppressed
/// as normal end-to-end line connectivity, in which case `suppressed` is
/// bumped).
fn evaluate_pair(
    a: &FeatureGeom,
    b: &FeatureGeom,
    same_layer: bool,
    line_connection_allowance: f64,
    suppressed: &mut usize,
) -> Option<(Geometry, f64)> {
    if !bbox_overlap(&a.bbox, &b.bbox) {
        return None;
    }
    let mut inter = a.buffered.intersection(&b.buffered);
    if inter.unsigned_area() <= SLIVER_AREA_EPS {
        return None;
    }
    if same_layer && line_connection_allowance > 0.0 && a.is_line && b.is_line {
        let shared = shared_endpoints(&a.endpoints, &b.endpoints);
        if !shared.is_empty() {
            let discs = union_discs(&shared, line_connection_allowance);
            inter = inter.difference(&discs);
            if inter.unsigned_area() <= SLIVER_AREA_EPS {
                *suppressed += 1;
                return None;
            }
        }
    }
    let area = inter.unsigned_area();
    Some((multipolygon_to_geometry(&inter), area))
}

/// The first/last coordinate of every part of a line geometry.
fn line_endpoints(g: &GeoGeometry) -> Vec<(f64, f64)> {
    let push_ends = |ls: &LineString, out: &mut Vec<(f64, f64)>| {
        if let Some(first) = ls.0.first() {
            out.push((first.x, first.y));
        }
        if ls.0.len() > 1 {
            if let Some(last) = ls.0.last() {
                out.push((last.x, last.y));
            }
        }
    };
    let mut out = Vec::new();
    match g {
        GeoGeometry::LineString(ls) => push_ends(ls, &mut out),
        GeoGeometry::MultiLineString(mls) => {
            for ls in &mls.0 {
                push_ends(ls, &mut out);
            }
        }
        _ => {}
    }
    out
}

/// Endpoints of `a` that coincide (within a tight snap tolerance) with an
/// endpoint of `b`.
fn shared_endpoints(a: &[(f64, f64)], b: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for &(ax, ay) in a {
        for &(bx, by) in b {
            if (ax - bx).abs() < ENDPOINT_SNAP_EPS && (ay - by).abs() < ENDPOINT_SNAP_EPS {
                out.push((ax, ay));
            }
        }
    }
    out
}

/// Union of discs of `radius` centred on each of `points`.
fn union_discs(points: &[(f64, f64)], radius: f64) -> MultiPolygon {
    let mut acc = MultiPolygon::<f64>::new(vec![]);
    for &(x, y) in points {
        let disc = GeoGeometry::Point(Point::new(x, y)).buffer(radius);
        acc = acc.union(&disc);
    }
    acc
}

fn bbox(mp: &MultiPolygon) -> [f64; 4] {
    match mp.bounding_rect() {
        Some(r) => [r.min().x, r.min().y, r.max().x, r.max().y],
        None => [
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        ],
    }
}

fn bbox_overlap(a: &[f64; 4], b: &[f64; 4]) -> bool {
    a[0] <= b[2] && b[0] <= a[2] && a[1] <= b[3] && b[1] <= a[3]
}

// ── geo <-> wbvector geometry conversion ───────────────────────────────────

fn to_geo_geometry(g: &Geometry) -> Option<GeoGeometry> {
    let geom = match g {
        Geometry::Point(c) => GeoGeometry::Point(Point::new(c.x, c.y)),
        Geometry::MultiPoint(cs) => GeoGeometry::MultiPoint(MultiPoint(
            cs.iter().map(|c| Point::new(c.x, c.y)).collect(),
        )),
        Geometry::LineString(cs) => GeoGeometry::LineString(coords_to_linestring(cs)),
        Geometry::MultiLineString(ls) => GeoGeometry::MultiLineString(MultiLineString(
            ls.iter().map(|c| coords_to_linestring(c)).collect(),
        )),
        Geometry::Polygon {
            exterior,
            interiors,
        } => GeoGeometry::Polygon(rings_to_polygon(exterior, interiors)),
        Geometry::MultiPolygon(parts) => GeoGeometry::MultiPolygon(MultiPolygon(
            parts.iter().map(|(e, i)| rings_to_polygon(e, i)).collect(),
        )),
        Geometry::GeometryCollection(_) => return None,
    };
    Some(geom)
}

fn coords_to_linestring(coords: &[Coord]) -> LineString {
    LineString::new(
        coords
            .iter()
            .map(|c| geo::Coord { x: c.x, y: c.y })
            .collect(),
    )
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
            .map(|c| geo::Coord { x: c.x, y: c.y })
            .collect(),
    )
}

fn multipolygon_to_geometry(mp: &MultiPolygon) -> Geometry {
    if mp.0.len() == 1 {
        let (exterior, interiors) = polygon_to_rings(&mp.0[0]);
        Geometry::Polygon {
            exterior,
            interiors,
        }
    } else {
        Geometry::MultiPolygon(mp.0.iter().map(polygon_to_rings).collect())
    }
}

fn polygon_to_rings(poly: &Polygon) -> (Ring, Vec<Ring>) {
    (
        linestring_to_ring(poly.exterior()),
        poly.interiors().iter().map(linestring_to_ring).collect(),
    )
}

fn linestring_to_ring(ls: &LineString) -> Ring {
    let mut coords: Vec<Coord> = ls.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
    if coords.len() >= 2 && coords.first() == coords.last() {
        coords.pop();
    }
    Ring::new(coords)
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    symbol_width: f64,
    conflict_symbol_width: Option<f64>,
    conflict_distance: f64,
    line_connection_allowance: f64,
    conflict: Option<String>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let symbol_width = opt_f64(args, "symbol_width")?.ok_or_else(|| {
        ToolError::Validation("missing required parameter 'symbol_width'".to_string())
    })?;
    if !(symbol_width.is_finite() && symbol_width > 0.0) {
        return Err(ToolError::Validation(
            "'symbol_width' must be a positive number".to_string(),
        ));
    }
    let conflict_symbol_width = opt_f64(args, "conflict_symbol_width")?;
    if let Some(w) = conflict_symbol_width {
        if !(w.is_finite() && w > 0.0) {
            return Err(ToolError::Validation(
                "'conflict_symbol_width' must be a positive number".to_string(),
            ));
        }
    }
    let conflict_distance = opt_f64(args, "conflict_distance")?.unwrap_or(0.0);
    if !(conflict_distance.is_finite() && conflict_distance >= 0.0) {
        return Err(ToolError::Validation(
            "'conflict_distance' must be a non-negative number".to_string(),
        ));
    }
    let line_connection_allowance = opt_f64(args, "line_connection_allowance")?.unwrap_or(0.0);
    if !(line_connection_allowance.is_finite() && line_connection_allowance >= 0.0) {
        return Err(ToolError::Validation(
            "'line_connection_allowance' must be a non-negative number".to_string(),
        ));
    }
    let conflict = parse_optional_str(args, "conflict")?.map(str::to_string);
    Ok(Params {
        symbol_width,
        conflict_symbol_width,
        conflict_distance,
        line_connection_allowance,
        conflict,
    })
}

fn opt_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
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

    fn line_layer(lines: &[Vec<(f64, f64)>]) -> String {
        let mut l = Layer::new("lines")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        for pts in lines {
            l.add_feature(
                Some(Geometry::LineString(
                    pts.iter().map(|(x, y)| Coord::xy(*x, *y)).collect(),
                )),
                &[],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        DetectGraphicConflictTool.run(&args, &ctx()).unwrap()
    }

    /// Two parallel lines closer together than the sum of their half symbol
    /// widths produce at least one conflict polygon with positive area.
    #[test]
    fn flags_close_parallel_lines() {
        let a: Vec<(f64, f64)> = (0..=10).map(|i| (i as f64 * 5.0, 0.0)).collect();
        let b: Vec<(f64, f64)> = (0..=10).map(|i| (i as f64 * 5.0, 3.0)).collect();
        let path = line_layer(&[a, b]);
        // symbol_width 10 -> half 5 each; lines are only 3 apart, so buffers overlap.
        let out = run(json!({ "input": path, "symbol_width": 10.0 }));
        assert_eq!(out.outputs["same_layer"], json!(true));
        let count = out.outputs["conflict_polygon_count"].as_u64().unwrap();
        assert!(count >= 1, "expected at least one conflict polygon");
        assert!(out.outputs["total_overlap_area"].as_f64().unwrap() > 0.0);
    }

    /// The same lines, far enough apart, produce no conflict.
    #[test]
    fn clears_far_apart_lines() {
        let a: Vec<(f64, f64)> = (0..=10).map(|i| (i as f64 * 5.0, 0.0)).collect();
        let b: Vec<(f64, f64)> = (0..=10).map(|i| (i as f64 * 5.0, 100.0)).collect();
        let path = line_layer(&[a, b]);
        let out = run(json!({ "input": path, "symbol_width": 10.0 }));
        assert_eq!(out.outputs["conflict_polygon_count"], json!(0));
    }

    /// Two lines joined end-to-end always overlap right at the shared vertex;
    /// a large enough connection allowance suppresses that as normal
    /// connectivity, not a conflict.
    #[test]
    fn connection_allowance_suppresses_shared_endpoint() {
        let a: Vec<(f64, f64)> = vec![(0.0, 0.0), (10.0, 0.0)];
        let b: Vec<(f64, f64)> = vec![(10.0, 0.0), (10.0, 10.0)];
        let path = line_layer(&[a, b]);
        // symbol_width 6 -> half 3 each; the L-joint buffers overlap near (10, 0).
        let without = run(json!({ "input": path.clone(), "symbol_width": 6.0 }));
        assert!(
            without.outputs["conflict_polygon_count"].as_u64().unwrap() >= 1,
            "expected the shared-endpoint overlap to be reported without an allowance"
        );

        let with = run(json!({
            "input": path, "symbol_width": 6.0, "line_connection_allowance": 5.0
        }));
        assert_eq!(with.outputs["conflict_polygon_count"], json!(0));
        assert!(
            with.outputs["suppressed_connections"].as_u64().unwrap() >= 1,
            "expected the connection to be recorded as suppressed"
        );
    }

    /// A separate 'conflict' layer is compared against 'input' (not against
    /// itself), and every input x conflict pair is a candidate.
    #[test]
    fn detects_conflict_between_two_layers() {
        let input_path = line_layer(&[vec![(0.0, 0.0), (10.0, 0.0)]]);
        let conflict_path = line_layer(&[vec![(0.0, 2.0), (10.0, 2.0)]]);
        let out = run(json!({
            "input": input_path, "conflict": conflict_path,
            "symbol_width": 6.0, "conflict_symbol_width": 6.0
        }));
        assert_eq!(out.outputs["same_layer"], json!(false));
        assert_eq!(out.outputs["conflict_polygon_count"], json!(1));
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = DetectGraphicConflictTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        let path = line_layer(&[vec![(0.0, 0.0), (10.0, 0.0)]]);
        assert!(bad(json!({})).is_err(), "missing input");
        assert!(
            bad(json!({ "input": path.clone() })).is_err(),
            "missing symbol_width"
        );
        assert!(
            bad(json!({ "input": path.clone(), "symbol_width": 0.0 })).is_err(),
            "zero symbol_width"
        );
        assert!(
            bad(json!({ "input": path.clone(), "symbol_width": 5.0, "conflict_distance": -1.0 }))
                .is_err(),
            "negative conflict_distance"
        );
        assert!(
            bad(json!({
                "input": path.clone(), "symbol_width": 5.0, "line_connection_allowance": -1.0
            }))
            .is_err(),
            "negative line_connection_allowance"
        );
        assert!(
            bad(json!({ "input": path, "symbol_width": 5.0 })).is_ok(),
            "valid parameters"
        );
    }
}
