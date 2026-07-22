//! GeoLibre tool: densify vector geometries along the true geodesic (or rhumb)
//! path between existing vertices.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Geodetic Densify*. The bundled
//! `densify_features`-style tools interpolate linearly in coordinate space,
//! which is correct for a projected CRS but wrong for long segments in
//! geographic (lon/lat) coordinates: a straight chord between two far-apart
//! lon/lat vertices does **not** follow the shortest path over the ellipsoid
//! (the geodesic) or a constant-bearing course (the rhumb line/loxodrome) —
//! it cuts across meridians at a different, misleading angle. This tool
//! assumes the input is geographic (EPSG:4326-style lon/lat) and replaces
//! each original segment with extra vertices placed along the true geodesic
//! or rhumb path, spaced so consecutive vertices are no more than
//! `max_segment_length` meters apart (or an exact `vertices_per_segment`
//! count of evenly-spaced insertions).
//!
//! # Geodesic math
//!
//! `geo` 0.33's geodesic algorithms (`Geodesic`/`Rhumb` + the `Bearing`,
//! `Destination`, `Distance`, `InterpolatePoint` traits) are unconditional
//! dependencies of the crate — not gated behind any Cargo feature — so they
//! compile with `default-features = false`. `Geodesic` wraps
//! `geographiclib-rs`'s pure-Rust port of Karney (2013)'s algorithms for
//! geodesics on the WGS84 ellipsoid; `Rhumb` implements loxodrome
//! bearing/destination/distance on a spherical approximation. Both expose
//! `points_along_line(start, end, max_distance, include_ends)`, which is
//! exactly the "insert vertices every N meters" operation this tool needs,
//! and `point_at_ratio_between` for the fixed-vertex-count mode.
//!
//! # Antimeridian handling
//!
//! For **line** geometries, if densifying a segment produces two consecutive
//! points whose longitude jumps by more than 180° (a wrap around ±180°), the
//! output is split at the antimeridian: a synthetic vertex is inserted at
//! exactly ±180° (latitude linearly interpolated) and the line continues as a
//! new part, so the result is a `MultiLineString` with no part crossing the
//! dateline. **Scope cut:** this splitting is only applied to line output.
//! Polygon rings are densified in place without being re-split into a
//! `MultiPolygon` at the dateline — true dateline-aware polygon clipping
//! needs general polygon boolean algebra, which duplicates GEOS/PROJ
//! functionality and is out of scope for a pure-Rust/WASM tool this size.
//! Polygons that do not cross the antimeridian (the overwhelming majority)
//! are unaffected by this cut.

use std::collections::BTreeMap;

use geo::{Geodesic, InterpolatePoint, Point, Rhumb};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, Geometry, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct GeodeticDensifyTool;

impl Tool for GeodeticDensifyTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "geodetic_densify",
            display_name: "Geodetic Densify",
            summary: "Insert vertices along the true geodesic or rhumb (constant-bearing) path between existing vertices of a geographic (lon/lat) line or polygon layer, like ArcGIS's Geodetic Densify — so long segments follow the earth's curvature instead of a planar chord.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input line or polygon vector layer, assumed geographic (lon/lat, EPSG:4326-style).",
                    required: true,
                },
                ToolParamSpec {
                    name: "geodetic_type",
                    description: "Path type between vertices: 'geodesic' (default; shortest path over the WGS84 ellipsoid) or 'rhumb'/'loxodrome' (constant true bearing).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_segment_length",
                    description: "Target maximum spacing, in meters, between consecutive vertices after densifying. Exactly one of max_segment_length / vertices_per_segment must be given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "vertices_per_segment",
                    description: "Alternative to max_segment_length: insert exactly this many evenly-spaced vertices into every original segment, regardless of its length. Exactly one of max_segment_length / vertices_per_segment must be given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output path; if omitted, stored in memory.",
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
                "missing required string parameter 'input'".into(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;

        let mut out = wbvector::Layer::new("geodetic_densified");
        out.geom_type = layer.geom_type;
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for f in layer.schema.fields() {
            out.add_field(f.clone());
        }

        let mut vertices_in = 0usize;
        let mut vertices_out = 0usize;
        let mut features_densified = 0usize;
        let mut features_passthrough = 0usize;

        for feature in layer.iter() {
            let geometry = match &feature.geometry {
                None => None,
                Some(g) => {
                    let (new_geom, n_in, n_out, touched) = densify_geometry(g, &prm);
                    vertices_in += n_in;
                    vertices_out += n_out;
                    if touched {
                        features_densified += 1;
                    } else {
                        features_passthrough += 1;
                    }
                    Some(new_geom)
                }
            };
            out.push(Feature {
                fid: 0,
                geometry,
                attributes: feature.attributes.clone(),
            });
        }

        ctx.progress.info(&format!(
            "geodetic_densify: {features_densified} feature(s) densified, {features_passthrough} passed through unchanged, {vertices_in} -> {vertices_out} vertices ({})",
            prm.geodetic_type.as_str()
        ));

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("features_densified".to_string(), json!(features_densified));
        outputs.insert(
            "features_passthrough".to_string(),
            json!(features_passthrough),
        );
        outputs.insert("vertices_in".to_string(), json!(vertices_in));
        outputs.insert("vertices_out".to_string(), json!(vertices_out));
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum GeodeticType {
    Geodesic,
    Rhumb,
}

impl GeodeticType {
    fn parse(s: &str) -> Option<GeodeticType> {
        match s.trim().to_ascii_lowercase().as_str() {
            "geodesic" | "geodetic" => Some(GeodeticType::Geodesic),
            "rhumb" | "loxodrome" | "rhumb_line" | "rhumbline" => Some(GeodeticType::Rhumb),
            _ => None,
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            GeodeticType::Geodesic => "geodesic",
            GeodeticType::Rhumb => "rhumb",
        }
    }
}

struct Params {
    geodetic_type: GeodeticType,
    max_segment_length: Option<f64>,
    vertices_per_segment: Option<u64>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let geodetic_type = match parse_optional_str(args, "geodetic_type")? {
        None => GeodeticType::Geodesic,
        Some(s) => GeodeticType::parse(s).ok_or_else(|| {
            ToolError::Validation(format!(
                "parameter 'geodetic_type' must be 'geodesic' or 'rhumb' (got '{s}')"
            ))
        })?,
    };
    let max_segment_length = parse_optional_f64(args, "max_segment_length")?;
    if let Some(v) = max_segment_length {
        if !(v.is_finite() && v > 0.0) {
            return Err(ToolError::Validation(
                "parameter 'max_segment_length' must be a positive number of meters".into(),
            ));
        }
    }
    let vertices_per_segment = parse_optional_u64(args, "vertices_per_segment")?;
    if let Some(v) = vertices_per_segment {
        if v == 0 {
            return Err(ToolError::Validation(
                "parameter 'vertices_per_segment' must be at least 1".into(),
            ));
        }
    }
    match (max_segment_length, vertices_per_segment) {
        (Some(_), None) | (None, Some(_)) => {}
        (None, None) => {
            return Err(ToolError::Validation(
                "exactly one of 'max_segment_length' or 'vertices_per_segment' is required".into(),
            ));
        }
        (Some(_), Some(_)) => {
            return Err(ToolError::Validation(
                "provide only one of 'max_segment_length' or 'vertices_per_segment', not both"
                    .into(),
            ));
        }
    }
    Ok(Params {
        geodetic_type,
        max_segment_length,
        vertices_per_segment,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
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

fn parse_optional_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64().or_else(|| n.as_f64().map(|f| f as u64))),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be an integer"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be an integer"
        ))),
    }
}

// ── Core densify ─────────────────────────────────────────────────────────────

/// Points strictly between `a` and `b` (exclusive), ordered from `a` to `b`.
fn intermediate_points(a: (f64, f64), b: (f64, f64), prm: &Params) -> Vec<(f64, f64)> {
    if a == b {
        return Vec::new();
    }
    let start = Point::new(a.0, a.1);
    let end = Point::new(b.0, b.1);
    match (prm.max_segment_length, prm.vertices_per_segment) {
        (Some(max_len), None) => {
            let pts: Vec<Point<f64>> = match prm.geodetic_type {
                GeodeticType::Geodesic => Geodesic
                    .points_along_line(start, end, max_len, false)
                    .collect(),
                GeodeticType::Rhumb => Rhumb
                    .points_along_line(start, end, max_len, false)
                    .collect(),
            };
            pts.into_iter().map(|p| (p.x(), p.y())).collect()
        }
        (None, Some(n)) => {
            let n = n.max(1);
            let mut out = Vec::with_capacity(n as usize);
            for k in 1..=n {
                let ratio = k as f64 / (n as f64 + 1.0);
                let p = match prm.geodetic_type {
                    GeodeticType::Geodesic => Geodesic.point_at_ratio_between(start, end, ratio),
                    GeodeticType::Rhumb => Rhumb.point_at_ratio_between(start, end, ratio),
                };
                out.push((p.x(), p.y()));
            }
            out
        }
        _ => unreachable!("validate() enforces exactly one of the two params"),
    }
}

/// Densifies an open polyline (no implicit closure): every original vertex is
/// kept, with intermediate points inserted between each consecutive pair.
fn densify_open(coords: &[(f64, f64)], prm: &Params) -> Vec<(f64, f64)> {
    if coords.len() < 2 {
        return coords.to_vec();
    }
    let mut out = Vec::with_capacity(coords.len() * 2);
    for w in coords.windows(2) {
        out.push(w[0]);
        out.extend(intermediate_points(w[0], w[1], prm));
    }
    out.push(*coords.last().unwrap());
    out
}

/// Densifies a closed ring (no duplicated closing vertex, per `wbvector::Ring`
/// convention): every original vertex is kept, with intermediate points
/// inserted along each edge including the closing edge (last -> first).
fn densify_ring(coords: &[(f64, f64)], prm: &Params) -> Vec<(f64, f64)> {
    let n = coords.len();
    if n < 2 {
        return coords.to_vec();
    }
    let mut out = Vec::with_capacity(n * 2);
    for i in 0..n {
        let a = coords[i];
        let b = coords[(i + 1) % n];
        out.push(a);
        out.extend(intermediate_points(a, b, prm));
    }
    out
}

/// Splits a densified open polyline into parts wherever consecutive vertices
/// jump more than 180° in longitude (an antimeridian crossing), inserting a
/// synthetic vertex at exactly ±180° so no output part straddles the dateline.
fn split_at_antimeridian(points: &[(f64, f64)]) -> Vec<Vec<(f64, f64)>> {
    if points.len() < 2 {
        return vec![points.to_vec()];
    }
    let mut parts: Vec<Vec<(f64, f64)>> = vec![vec![points[0]]];
    for w in points.windows(2) {
        let (lon1, lat1) = w[0];
        let (lon2, lat2) = w[1];
        let d = lon2 - lon1;
        if d > 180.0 {
            // Wrapped westbound through -180/+180 (unwrapped lon2 = lon2 - 360).
            let unwrapped = lon2 - 360.0;
            let t = (-180.0 - lon1) / (unwrapped - lon1);
            let lat_cross = lat1 + t * (lat2 - lat1);
            parts.last_mut().unwrap().push((-180.0, lat_cross));
            parts.push(vec![(180.0, lat_cross)]);
            parts.last_mut().unwrap().push((lon2, lat2));
        } else if d < -180.0 {
            // Wrapped eastbound through +180/-180 (unwrapped lon2 = lon2 + 360).
            let unwrapped = lon2 + 360.0;
            let t = (180.0 - lon1) / (unwrapped - lon1);
            let lat_cross = lat1 + t * (lat2 - lat1);
            parts.last_mut().unwrap().push((180.0, lat_cross));
            parts.push(vec![(-180.0, lat_cross)]);
            parts.last_mut().unwrap().push((lon2, lat2));
        } else {
            parts.last_mut().unwrap().push((lon2, lat2));
        }
    }
    parts
}

fn coords_to_ring(coords: &[(f64, f64)]) -> Ring {
    Ring::new(coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect())
}

fn coords_to_vec(coords: &[(f64, f64)]) -> Vec<Coord> {
    coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect()
}

/// Densifies one geometry. Returns (new geometry, input vertex count, output
/// vertex count, whether anything was actually inserted).
fn densify_geometry(geom: &Geometry, prm: &Params) -> (Geometry, usize, usize, bool) {
    match geom {
        Geometry::Point(_) | Geometry::MultiPoint(_) => {
            let n = geom_vertex_count(geom);
            (geom.clone(), n, n, false)
        }
        Geometry::LineString(coords) => {
            let pts: Vec<(f64, f64)> = coords.iter().map(|c| (c.x, c.y)).collect();
            let n_in = pts.len();
            let densified = densify_open(&pts, prm);
            let parts = split_at_antimeridian(&densified);
            let n_out: usize = parts.iter().map(|p| p.len()).sum();
            let touched = n_out != n_in || parts.len() > 1;
            let new_geom = if parts.len() == 1 {
                Geometry::LineString(coords_to_vec(&parts[0]))
            } else {
                Geometry::MultiLineString(parts.iter().map(|p| coords_to_vec(p)).collect())
            };
            (new_geom, n_in, n_out, touched)
        }
        Geometry::MultiLineString(lines) => {
            let mut n_in = 0usize;
            let mut n_out = 0usize;
            let mut all_parts: Vec<Vec<Coord>> = Vec::new();
            for coords in lines {
                let pts: Vec<(f64, f64)> = coords.iter().map(|c| (c.x, c.y)).collect();
                n_in += pts.len();
                let densified = densify_open(&pts, prm);
                let parts = split_at_antimeridian(&densified);
                for p in &parts {
                    n_out += p.len();
                    all_parts.push(coords_to_vec(p));
                }
            }
            let touched = n_out != n_in || all_parts.len() != lines.len();
            (Geometry::MultiLineString(all_parts), n_in, n_out, touched)
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            let (new_ext, n_in_e, n_out_e) = densify_ring_geom(exterior, prm);
            let mut n_in = n_in_e;
            let mut n_out = n_out_e;
            let mut new_interiors = Vec::with_capacity(interiors.len());
            for ring in interiors {
                let (r, ni, no) = densify_ring_geom(ring, prm);
                n_in += ni;
                n_out += no;
                new_interiors.push(r);
            }
            let touched = n_out != n_in;
            (
                Geometry::Polygon {
                    exterior: new_ext,
                    interiors: new_interiors,
                },
                n_in,
                n_out,
                touched,
            )
        }
        Geometry::MultiPolygon(parts) => {
            let mut n_in = 0usize;
            let mut n_out = 0usize;
            let mut new_parts = Vec::with_capacity(parts.len());
            for (exterior, interiors) in parts {
                let (new_ext, ni, no) = densify_ring_geom(exterior, prm);
                n_in += ni;
                n_out += no;
                let mut new_interiors = Vec::with_capacity(interiors.len());
                for ring in interiors {
                    let (r, ni2, no2) = densify_ring_geom(ring, prm);
                    n_in += ni2;
                    n_out += no2;
                    new_interiors.push(r);
                }
                new_parts.push((new_ext, new_interiors));
            }
            let touched = n_out != n_in;
            (Geometry::MultiPolygon(new_parts), n_in, n_out, touched)
        }
        Geometry::GeometryCollection(geoms) => {
            let mut n_in = 0usize;
            let mut n_out = 0usize;
            let mut touched = false;
            let mut new_geoms = Vec::with_capacity(geoms.len());
            for g in geoms {
                let (ng, ni, no, t) = densify_geometry(g, prm);
                n_in += ni;
                n_out += no;
                touched |= t;
                new_geoms.push(ng);
            }
            (
                Geometry::GeometryCollection(new_geoms),
                n_in,
                n_out,
                touched,
            )
        }
    }
}

fn densify_ring_geom(ring: &Ring, prm: &Params) -> (Ring, usize, usize) {
    let pts: Vec<(f64, f64)> = ring.coords().iter().map(|c| (c.x, c.y)).collect();
    let n_in = pts.len();
    let densified = densify_ring(&pts, prm);
    let n_out = densified.len();
    (coords_to_ring(&densified), n_in, n_out)
}

fn geom_vertex_count(geom: &Geometry) -> usize {
    match geom {
        Geometry::Point(_) => 1,
        Geometry::MultiPoint(pts) => pts.len(),
        Geometry::LineString(coords) => coords.len(),
        Geometry::MultiLineString(lines) => lines.iter().map(|l| l.len()).sum(),
        Geometry::Polygon {
            exterior,
            interiors,
        } => exterior.len() + interiors.iter().map(Ring::len).sum::<usize>(),
        Geometry::MultiPolygon(parts) => parts
            .iter()
            .map(|(e, ints)| e.len() + ints.iter().map(Ring::len).sum::<usize>())
            .sum(),
        Geometry::GeometryCollection(geoms) => geoms.iter().map(geom_vertex_count).sum(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use geo::{Bearing, Distance};
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn line_layer(coords: &[(f64, f64)]) -> String {
        let mut l = Layer::new("lines")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.add_feature(
            Some(Geometry::line_string(
                coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect(),
            )),
            &[("name", "seg".into())],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = GeodeticDensifyTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// A long east-west segment at high latitude, densified along the true
    /// geodesic, must bow poleward of the planar-chord midpoint: at 70°N the
    /// geodesic between two points 100° of longitude apart passes noticeably
    /// closer to the pole than the straight average of the two latitudes
    /// (which is just 70° here, since both endpoints share the same lat).
    #[test]
    fn geodesic_bows_poleward_at_high_latitude() {
        let input = line_layer(&[(-40.0, 70.0), (60.0, 70.0)]);
        let (_out, layer) = run(json!({
            "input": input,
            "geodetic_type": "geodesic",
            "max_segment_length": 200_000.0,
        }));
        let Geometry::LineString(coords) = layer.features[0].geometry.as_ref().unwrap() else {
            panic!("expected LineString");
        };
        assert!(coords.len() > 2, "expected inserted vertices");
        let planar_lat = 70.0; // both endpoints share this latitude
        let max_lat = coords.iter().map(|c| c.y).fold(f64::MIN, f64::max);
        assert!(
            max_lat > planar_lat + 1.0,
            "geodesic midsection should bow well north of {planar_lat}, got max {max_lat}"
        );
    }

    /// With max_segment_length, the number of output segments on one input
    /// segment is within 1 of ceil(geodesic_distance / max_len) (the extra
    /// slack absorbs floating-point drift in repeated-addition step
    /// accumulation), and — the property that actually matters — every
    /// output segment's geodesic length is at most max_segment_length.
    #[test]
    fn segment_count_matches_ceil_distance_over_max_len() {
        let a = Point::new(-40.0, 70.0);
        let b = Point::new(60.0, 70.0);
        let dist = Geodesic.distance(a, b);
        let max_len = 250_000.0;
        let expected_segments = (dist / max_len).ceil() as usize;
        let input = line_layer(&[(-40.0, 70.0), (60.0, 70.0)]);
        let (_out, layer) = run(json!({
            "input": input,
            "geodetic_type": "geodesic",
            "max_segment_length": max_len,
        }));
        let Geometry::LineString(coords) = layer.features[0].geometry.as_ref().unwrap() else {
            panic!("expected LineString");
        };
        let actual_segments = coords.len() - 1;
        assert!(
            actual_segments.abs_diff(expected_segments) <= 1,
            "expected ~{expected_segments} output segments (ceil({dist} / {max_len})), got {actual_segments}"
        );
        for w in coords.windows(2) {
            let p1 = Point::new(w[0].x, w[0].y);
            let p2 = Point::new(w[1].x, w[1].y);
            let seg_len = Geodesic.distance(p1, p2);
            assert!(
                seg_len <= max_len + 1.0,
                "output segment {seg_len}m exceeds max_segment_length {max_len}m"
            );
        }
    }

    /// Rhumb-line densification keeps a constant bearing between every
    /// consecutive pair of output vertices (within numerical tolerance),
    /// unlike the geodesic, whose bearing continuously changes.
    #[test]
    fn rhumb_keeps_constant_bearing() {
        let input = line_layer(&[(-40.0, 10.0), (60.0, 55.0)]);
        let (_out, layer) = run(json!({
            "input": input,
            "geodetic_type": "rhumb",
            "vertices_per_segment": 6u64,
        }));
        let Geometry::LineString(coords) = layer.features[0].geometry.as_ref().unwrap() else {
            panic!("expected LineString");
        };
        assert_eq!(coords.len(), 8); // 2 endpoints + 6 inserted
        let bearings: Vec<f64> = coords
            .windows(2)
            .map(|w| {
                let p1 = Point::new(w[0].x, w[0].y);
                let p2 = Point::new(w[1].x, w[1].y);
                Rhumb.bearing(p1, p2)
            })
            .collect();
        let first = bearings[0];
        for b in &bearings {
            assert!(
                (b - first).abs() < 1e-6,
                "rhumb bearing should stay constant: {first} vs {b}"
            );
        }
    }

    /// Point geometries pass through untouched (no vertices to insert).
    #[test]
    fn points_pass_through_unchanged() {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.add_feature(Some(Geometry::point(10.0, 20.0)), &[("name", "p".into())])
            .unwrap();
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);
        let (out, layer) = run(json!({
            "input": input,
            "geodetic_type": "geodesic",
            "max_segment_length": 1000.0,
        }));
        assert_eq!(out.outputs["features_passthrough"], json!(1));
        assert_eq!(out.outputs["features_densified"], json!(0));
        match layer.features[0].geometry.as_ref().unwrap() {
            Geometry::Point(c) => {
                assert!((c.x - 10.0).abs() < 1e-9);
                assert!((c.y - 20.0).abs() < 1e-9);
            }
            _ => panic!("expected Point"),
        }
    }

    /// A segment crossing the antimeridian is split into a MultiLineString
    /// with no part straddling ±180°.
    #[test]
    fn splits_at_antimeridian() {
        let input = line_layer(&[(170.0, 10.0), (-170.0, 12.0)]);
        let (_out, layer) = run(json!({
            "input": input,
            "geodetic_type": "geodesic",
            "max_segment_length": 100_000.0,
        }));
        match layer.features[0].geometry.as_ref().unwrap() {
            Geometry::MultiLineString(parts) => {
                assert_eq!(parts.len(), 2, "expected a split into two parts");
                for part in parts {
                    for w in part.windows(2) {
                        assert!(
                            (w[1].x - w[0].x).abs() <= 180.0,
                            "part should not itself straddle the antimeridian"
                        );
                    }
                }
            }
            other => panic!("expected MultiLineString after antimeridian split, got {other:?}"),
        }
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            GeodeticDensifyTool.validate(&args)
        };
        // Missing input.
        assert!(bad(json!({ "max_segment_length": 1000.0 })).is_err());
        // Neither max_segment_length nor vertices_per_segment.
        assert!(bad(json!({ "input": "a.geojson" })).is_err());
        // Both given.
        assert!(bad(json!({
            "input": "a.geojson",
            "max_segment_length": 1000.0,
            "vertices_per_segment": 3,
        }))
        .is_err());
        // Non-positive max_segment_length.
        assert!(bad(json!({
            "input": "a.geojson", "max_segment_length": -5.0,
        }))
        .is_err());
        // Zero vertices_per_segment.
        assert!(bad(json!({
            "input": "a.geojson", "vertices_per_segment": 0,
        }))
        .is_err());
        // Unknown geodetic_type.
        assert!(bad(json!({
            "input": "a.geojson", "max_segment_length": 1000.0, "geodetic_type": "flat",
        }))
        .is_err());
        // Valid.
        assert!(bad(json!({
            "input": "a.geojson", "max_segment_length": 1000.0,
        }))
        .is_ok());
        assert!(bad(json!({
            "input": "a.geojson", "vertices_per_segment": 4, "geodetic_type": "rhumb",
        }))
        .is_ok());
    }
}
