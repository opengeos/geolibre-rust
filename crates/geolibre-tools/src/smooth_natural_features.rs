//! GeoLibre tool: smooth natural feature boundaries.
//!
//! Pure-Rust counterpart of [Smoothify](https://github.com/DPIRD-DMA/Smoothify):
//! it turns pixelated polygons and jagged lines — typically raster-to-vector
//! outputs such as classified land cover, water bodies, or vegetation masks —
//! into smooth, natural-looking curves while preserving each polygon's area.
//!
//! The pipeline per ring/line (mirroring Smoothify's two passes):
//! 1. De-noise with Douglas–Peucker at `segment_length` tolerance (removes the
//!    staircase artifacts of raster boundaries). Closed rings are opened at
//!    the vertex farthest from their centroid so the split point is stable and
//!    start-position bias is avoided (Smoothify instead merges several
//!    start-offset variants; the stable split achieves the same end cheaply).
//! 2. Densify so no segment exceeds 4x`segment_length`, bounding the scale of
//!    Chaikin's quarter-length corner cuts, then apply one Chaikin pass.
//! 3. Re-simplify at `segment_length / 5`, re-densify, and apply the final
//!    Chaikin corner-cutting passes (`iterations` of them).
//! 4. For polygons with `preserve_area`, restore the original area (Chaikin
//!    shrinks convex boundaries) by offsetting the smoothed rings along their
//!    miter vertex normals — exterior outward, holes inward — with the offset
//!    distance solved by Newton iteration (dA/dd = boundary length), to within
//!    0.01% of the original area.
//!
//! `segment_length` should match the source raster's cell size; when omitted
//! it is auto-detected as the layer's median edge length (raster-traced
//! boundaries step one cell at a time, so the median is the resolution). Per
//! feature the effective value is capped at 1/12 of the ring perimeter (or
//! line length) so clean, coarse geometries are not simplified away.
//!
//! Not implemented from Smoothify: dissolving adjacent geometries, joining
//! nearly-touching holes, and the morphological artifact repair (all need
//! polygon boolean ops). Point features pass through unchanged.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata,
    ToolParamSpec, ToolRunResult,
};
use wbvector::{Coord, Geometry, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Densification bound: Chaikin's corner cuts span a quarter of each segment,
/// so segments no longer than 4x`segment_length` keep cuts within one cell.
const SEGMENT_FACTOR: f64 = 4.0;
/// Second-pass simplification runs at `segment_length / 5` (Smoothify's ratio):
/// enough to erase first-pass Chaikin detail without moving the boundary.
const POST_SIMPLIFY_DIVISOR: f64 = 5.0;
/// Area preservation converges when within this fraction of the original area.
const AREA_REL_TOL: f64 = 1e-4;
/// Chaikin doubles the vertex count per pass; stop before exceeding this.
const MAX_VERTICES: usize = 16_384;
/// A feature's effective segment length is capped at `perimeter / 12` so an
/// auto-detected (layer-median) value cannot flatten a clean coarse shape.
const FEATURE_SIZE_DIVISOR: f64 = 12.0;

pub struct SmoothNaturalFeaturesTool;

impl Tool for SmoothNaturalFeaturesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "smooth_natural_features",
            display_name: "Smooth Natural Features",
            summary: "Smooth pixelated polygons and jagged lines (raster-to-vector outputs) into natural-looking curves via Chaikin corner cutting, preserving polygon area.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector file path with polygon or line features, format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "segment_length",
                    description: "Scale of the staircase noise to remove, in CRS units — ideally the source raster's cell size. Auto-detected from the layer's median edge length if omitted.",
                    required: false,
                },
                ToolParamSpec {
                    name: "iterations",
                    description: "Chaikin smoothing iterations (1-8); more iterations give smoother output with more vertices. Default 3.",
                    required: false,
                },
                ToolParamSpec {
                    name: "preserve_area",
                    description: "Restore each polygon's original area after smoothing (corner cutting shrinks convex shapes). Default true.",
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
            .ok_or_else(|| ToolError::Validation("missing required parameter 'input'".to_string()))?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let segment_length = match prm.segment_length {
            Some(v) => v,
            None => detect_segment_length(&layer).ok_or_else(|| {
                ToolError::Validation(
                    "could not auto-detect 'segment_length' (no line or polygon edges in input); pass it explicitly".to_string(),
                )
            })?,
        };
        ctx.progress.info(&format!(
            "smoothing {} feature(s) (segment length {segment_length:.4}, {} iteration(s))",
            layer.len(),
            prm.iterations
        ));

        layer.extent = None; // geometries change; invalidate the cached bbox
        let (mut smoothed, mut skipped) = (0usize, 0usize);
        for feature in layer.iter_mut() {
            match feature
                .geometry
                .as_ref()
                .and_then(|g| smooth_geometry(g, &prm, segment_length))
            {
                Some(geom) => {
                    feature.geometry = Some(geom);
                    smoothed += 1;
                }
                // Point (or missing) geometry: pass through unchanged.
                None => skipped += 1,
            }
        }
        ctx.progress.info(&format!(
            "{smoothed} smoothed, {skipped} passed through unchanged"
        ));

        let feature_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("smoothed_count".to_string(), json!(smoothed));
        outputs.insert("skipped_count".to_string(), json!(skipped));
        outputs.insert("segment_length".to_string(), json!(segment_length));
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    segment_length: Option<f64>,
    iterations: usize,
    preserve_area: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let segment_length = parse_optional_f64(args, "segment_length")?;
    if let Some(v) = segment_length {
        if !(v > 0.0 && v.is_finite()) {
            return Err(ToolError::Validation(
                "parameter 'segment_length' must be a positive number".to_string(),
            ));
        }
    }
    let iterations = match parse_optional_f64(args, "iterations")? {
        None => 3,
        Some(v) if v.fract() == 0.0 && (1.0..=8.0).contains(&v) => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "parameter 'iterations' must be an integer between 1 and 8".to_string(),
            ))
        }
    };
    let preserve_area = parse_optional_bool(args, "preserve_area")?.unwrap_or(true);
    Ok(Params { segment_length, iterations, preserve_area })
}

/// Parses an optional numeric parameter, accepting a JSON number or a numeric
/// string (host UIs often post form values as strings).
fn parse_optional_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s.trim().parse::<f64>().map(Some).map_err(|_| {
            ToolError::Validation(format!("parameter '{key}' must be a number"))
        }),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a number"
        ))),
    }
}

/// Parses an optional boolean parameter, accepting a JSON bool or a
/// "true"/"false" string (host UIs often post form values as strings).
fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "true" | "yes" | "1" => Ok(Some(true)),
            "false" | "no" | "0" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
}

// ── Geometry dispatch ─────────────────────────────────────────────────────────

/// Smooths a polygonal or linear geometry. Returns `None` for point (or other
/// unsupported) geometries, which pass through unchanged.
fn smooth_geometry(geom: &Geometry, prm: &Params, seg: f64) -> Option<Geometry> {
    match geom {
        Geometry::Polygon { exterior, interiors } => {
            let (exterior, interiors) = smooth_part(exterior, interiors, prm, seg);
            Some(Geometry::Polygon { exterior, interiors })
        }
        Geometry::MultiPolygon(parts) => Some(Geometry::MultiPolygon(
            parts
                .iter()
                .map(|(ext, holes)| smooth_part(ext, holes, prm, seg))
                .collect(),
        )),
        Geometry::LineString(coords) => {
            Some(Geometry::LineString(smooth_line(coords, prm, seg)))
        }
        Geometry::MultiLineString(lines) => Some(Geometry::MultiLineString(
            lines.iter().map(|l| smooth_line(l, prm, seg)).collect(),
        )),
        _ => None,
    }
}

/// Smooths one polygon part (exterior + holes), then restores its original
/// area if requested. Area is preserved per part: each part's exterior grows
/// and its holes shrink by the same offset distance.
fn smooth_part(exterior: &Ring, interiors: &[Ring], prm: &Params, seg: f64) -> (Ring, Vec<Ring>) {
    let ext_pts = ring_points(exterior);
    let hole_pts: Vec<Vec<P>> = interiors.iter().map(ring_points).collect();
    if ext_pts.len() < 3 {
        return (exterior.clone(), interiors.to_vec());
    }
    let target = part_area(&ext_pts, &hole_pts);
    let seg = seg.min(perimeter(&ext_pts, true) / FEATURE_SIZE_DIVISOR).max(f64::MIN_POSITIVE);

    let mut smooth_ext = smooth_closed_ring(&ext_pts, seg, prm.iterations);
    let mut smooth_holes: Vec<Vec<P>> = hole_pts
        .iter()
        .map(|h| smooth_closed_ring(h, seg, prm.iterations))
        .collect();
    if prm.preserve_area {
        (smooth_ext, smooth_holes) = preserve_part_area(target, smooth_ext, smooth_holes, seg);
    }
    (
        pts_to_ring(&smooth_ext),
        smooth_holes.iter().map(|h| pts_to_ring(h)).collect(),
    )
}

fn smooth_line(coords: &[Coord], prm: &Params, seg: f64) -> Vec<Coord> {
    let pts = dedup_points(coords);
    if pts.len() < 3 {
        return coords.to_vec();
    }
    let seg = seg.min(perimeter(&pts, false) / FEATURE_SIZE_DIVISOR).max(f64::MIN_POSITIVE);
    let mut cur = rdp(&pts, seg);
    cur = densify(&cur, SEGMENT_FACTOR * seg, false);
    cur = chaikin(&cur, 1, false);
    cur = rdp(&cur, seg / POST_SIMPLIFY_DIVISOR);
    cur = densify(&cur, SEGMENT_FACTOR * seg, false);
    cur = chaikin(&cur, prm.iterations, false);
    cur.iter().map(|p| Coord::xy(p.x, p.y)).collect()
}

// ── Smoothing pipeline ────────────────────────────────────────────────────────

/// Simplify → densify → Chaikin, twice (a light first pass, then the full
/// `iterations`), on an unclosed ring. Degenerate rings are returned as-is.
fn smooth_closed_ring(pts: &[P], seg: f64, iterations: usize) -> Vec<P> {
    if pts.len() < 4 {
        return pts.to_vec();
    }
    let mut cur = rdp_ring(pts, seg);
    cur = densify(&cur, SEGMENT_FACTOR * seg, true);
    cur = chaikin(&cur, 1, true);
    cur = rdp_ring(&cur, seg / POST_SIMPLIFY_DIVISOR);
    cur = densify(&cur, SEGMENT_FACTOR * seg, true);
    chaikin(&cur, iterations, true)
}

/// Chaikin corner cutting: each segment is replaced by its 1/4 and 3/4 points.
/// Open lines keep their endpoints fixed. Stops early rather than exceed
/// `MAX_VERTICES`.
fn chaikin(pts: &[P], iterations: usize, closed: bool) -> Vec<P> {
    let mut cur = pts.to_vec();
    for _ in 0..iterations {
        let n = cur.len();
        if n < 3 || n * 2 > MAX_VERTICES {
            break;
        }
        let mut next = Vec::with_capacity(n * 2 + 2);
        if closed {
            for i in 0..n {
                let (a, b) = (cur[i], cur[(i + 1) % n]);
                next.push(lerp(a, b, 0.25));
                next.push(lerp(a, b, 0.75));
            }
        } else {
            next.push(cur[0]);
            for i in 0..n - 1 {
                let (a, b) = (cur[i], cur[i + 1]);
                next.push(lerp(a, b, 0.25));
                next.push(lerp(a, b, 0.75));
            }
            next.push(cur[n - 1]);
        }
        cur = next;
    }
    cur
}

/// Inserts vertices so no segment is longer than `max_len`.
fn densify(pts: &[P], max_len: f64, closed: bool) -> Vec<P> {
    let n = pts.len();
    if n < 2 || max_len <= 0.0 {
        return pts.to_vec();
    }
    let edges = if closed { n } else { n - 1 };
    let mut out = Vec::with_capacity(n * 2);
    for i in 0..edges {
        let (a, b) = (pts[i], pts[(i + 1) % n]);
        out.push(a);
        let pieces = (dist(a, b) / max_len).ceil() as usize;
        for j in 1..pieces {
            out.push(lerp(a, b, j as f64 / pieces as f64));
        }
    }
    if !closed {
        out.push(pts[n - 1]);
    }
    out
}

// ── Area preservation ─────────────────────────────────────────────────────────

/// Restores a smoothed polygon part to its original (`target`) area by
/// offsetting the exterior outward and the holes inward along miter vertex
/// normals, solving the offset distance with Newton iteration (dA/dd equals
/// the boundary length). Falls back to the closest candidate seen — including
/// the unoffset rings — if the iteration fails to converge.
fn preserve_part_area(
    target: f64,
    ext: Vec<P>,
    holes: Vec<Vec<P>>,
    seg: f64,
) -> (Vec<P>, Vec<Vec<P>>) {
    if target.is_nan() || target <= 0.0 || ext.len() < 3 {
        return (ext, holes);
    }
    let tol = target * AREA_REL_TOL;
    let ext_sign = signed_area(&ext).signum();
    let hole_signs: Vec<f64> = holes.iter().map(|h| signed_area(h).signum()).collect();
    let max_d = 10.0 * seg;

    let mut d = 0.0;
    let mut best: Option<(f64, Vec<P>, Vec<Vec<P>>)> = None;
    for _ in 0..12 {
        let cand_ext = offset_ring(&ext, d);
        if signed_area(&cand_ext).signum() != ext_sign {
            // Overshot into self-inversion; retreat toward the smoothed rings.
            d *= 0.5;
            if d.abs() < 1e-12 {
                break;
            }
            continue;
        }
        // Holes shrink as the exterior grows; a hole that inverts has been
        // consumed by the offset and is dropped.
        let cand_holes: Vec<Vec<P>> = holes
            .iter()
            .zip(&hole_signs)
            .map(|(h, sign)| (offset_ring(h, -d), sign))
            .filter(|(h, sign)| h.len() >= 3 && signed_area(h).signum() == **sign)
            .map(|(h, _)| h)
            .collect();
        let area = part_area(&cand_ext, &cand_holes);
        let err = area - target;
        if best.as_ref().is_none_or(|(e, _, _)| err.abs() < *e) {
            best = Some((err.abs(), cand_ext.clone(), cand_holes.clone()));
        }
        if err.abs() <= tol {
            break;
        }
        let slope = perimeter(&cand_ext, true)
            + cand_holes.iter().map(|h| perimeter(h, true)).sum::<f64>();
        if slope.is_nan() || slope <= 1e-12 {
            break;
        }
        d = (d - err / slope).clamp(-max_d, max_d);
    }
    match best {
        Some((_, e, h)) => (e, h),
        None => (ext, holes),
    }
}

/// Offsets a ring's vertices along their miter normals so the ring's enclosed
/// area grows for positive `d` (regardless of winding). The miter length is
/// capped at 2x`d` to keep sharp corners bounded.
fn offset_ring(pts: &[P], d: f64) -> Vec<P> {
    let n = pts.len();
    if n < 3 || d == 0.0 {
        return pts.to_vec();
    }
    let orient = signed_area(pts).signum();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let prev = pts[(i + n - 1) % n];
        let cur = pts[i];
        let next = pts[(i + 1) % n];
        let n1 = outward_normal(prev, cur, orient);
        let n2 = outward_normal(cur, next, orient);
        let (mut bx, mut by) = (n1.0 + n2.0, n1.1 + n2.1);
        let blen = bx.hypot(by);
        if blen < 1e-12 {
            (bx, by) = n1;
        } else {
            bx /= blen;
            by /= blen;
        }
        // Moving both edges by d displaces their intersection by d / cos(α/2);
        // the cap at cos = 0.5 limits the miter to 2d at sharp corners.
        let cos_half = ((1.0 + (n1.0 * n2.0 + n1.1 * n2.1)) * 0.5).max(0.0).sqrt().max(0.5);
        out.push(P { x: cur.x + bx * d / cos_half, y: cur.y + by * d / cos_half });
    }
    out
}

/// Unit normal of the edge `a`→`b` pointing away from the ring's interior
/// (`orient` is the sign of the ring's shoelace area; CCW keeps the interior
/// on the left, so outward is the right-hand normal).
fn outward_normal(a: P, b: P, orient: f64) -> (f64, f64) {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let len = dx.hypot(dy);
    if len < 1e-12 {
        return (0.0, 0.0);
    }
    if orient >= 0.0 {
        (dy / len, -dx / len)
    } else {
        (-dy / len, dx / len)
    }
}

/// Area of a polygon part: |exterior| minus the |holes|.
fn part_area(ext: &[P], holes: &[Vec<P>]) -> f64 {
    signed_area(ext).abs() - holes.iter().map(|h| signed_area(h).abs()).sum::<f64>()
}

// ── Geometry primitives ───────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
struct P {
    x: f64,
    y: f64,
}

fn lerp(a: P, b: P, t: f64) -> P {
    P { x: a.x + (b.x - a.x) * t, y: a.y + (b.y - a.y) * t }
}

fn dist(a: P, b: P) -> f64 {
    (a.x - b.x).hypot(a.y - b.y)
}

fn perimeter(pts: &[P], closed: bool) -> f64 {
    let n = pts.len();
    if n < 2 {
        return 0.0;
    }
    let edges = if closed { n } else { n - 1 };
    (0..edges).map(|i| dist(pts[i], pts[(i + 1) % n])).sum()
}

/// Extracts a ring's vertices, dropping consecutive duplicates and the closing
/// duplicate if present.
fn ring_points(ring: &Ring) -> Vec<P> {
    let mut pts = dedup_points(ring.coords());
    while pts.len() >= 2 && dist(pts[0], *pts.last().unwrap()) <= 1e-12 {
        pts.pop();
    }
    pts
}

/// Extracts vertices, dropping consecutive duplicates.
fn dedup_points(coords: &[Coord]) -> Vec<P> {
    let mut pts: Vec<P> = Vec::with_capacity(coords.len());
    for c in coords {
        let p = P { x: c.x, y: c.y };
        if pts.last().is_none_or(|last| dist(*last, p) > 1e-12) {
            pts.push(p);
        }
    }
    pts
}

fn pts_to_ring(pts: &[P]) -> Ring {
    Ring::new(pts.iter().map(|p| Coord::xy(p.x, p.y)).collect())
}

/// Shoelace signed area of an unclosed ring (positive = CCW).
fn signed_area(pts: &[P]) -> f64 {
    let n = pts.len();
    if n < 3 {
        return 0.0;
    }
    let mut a = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        a += pts[i].x * pts[j].y - pts[j].x * pts[i].y;
    }
    a * 0.5
}

/// Distance from `p` to the segment `a`-`b`.
fn point_seg_dist(p: P, a: P, b: P) -> f64 {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let len2 = dx * dx + dy * dy;
    if len2 <= 0.0 {
        return dist(p, a);
    }
    let t = (((p.x - a.x) * dx + (p.y - a.y) * dy) / len2).clamp(0.0, 1.0);
    dist(p, P { x: a.x + t * dx, y: a.y + t * dy })
}

/// Douglas–Peucker on an open polyline; endpoints are always kept.
fn rdp(points: &[P], tol: f64) -> Vec<P> {
    let n = points.len();
    if n < 3 {
        return points.to_vec();
    }
    let mut keep = vec![false; n];
    keep[0] = true;
    keep[n - 1] = true;
    let mut stack = vec![(0usize, n - 1)];
    while let Some((i, j)) = stack.pop() {
        if j <= i + 1 {
            continue;
        }
        let (mut best, mut best_d) = (i + 1, -1.0);
        for (k, p) in points.iter().enumerate().take(j).skip(i + 1) {
            let d = point_seg_dist(*p, points[i], points[j]);
            if d > best_d {
                best_d = d;
                best = k;
            }
        }
        if best_d > tol {
            keep[best] = true;
            stack.push((i, best));
            stack.push((best, j));
        }
    }
    points
        .iter()
        .zip(&keep)
        .filter_map(|(p, k)| k.then_some(*p))
        .collect()
}

/// Douglas–Peucker for a closed ring: the ring is opened at the vertex farthest
/// from its centroid (almost certainly a real corner) so the split point is
/// stable, then simplified as a polyline.
fn rdp_ring(pts: &[P], tol: f64) -> Vec<P> {
    let n = pts.len();
    if n < 4 {
        return pts.to_vec();
    }
    let cx = pts.iter().map(|p| p.x).sum::<f64>() / n as f64;
    let cy = pts.iter().map(|p| p.y).sum::<f64>() / n as f64;
    let centroid = P { x: cx, y: cy };
    let start = (0..n)
        .max_by(|&a, &b| dist(pts[a], centroid).total_cmp(&dist(pts[b], centroid)))
        .unwrap_or(0);
    let closed: Vec<P> = (0..=n).map(|i| pts[(start + i) % n]).collect();
    let mut out = rdp(&closed, tol);
    out.pop(); // drop the duplicated closing vertex
    if out.len() < 3 {
        pts.to_vec()
    } else {
        out
    }
}

// ── Segment-length auto-detection ─────────────────────────────────────────────

/// Median edge length across the layer's line and polygon geometries. For
/// raster-traced boundaries (which step one cell at a time) this recovers the
/// source raster's resolution, matching Smoothify's `segment_length` default.
fn detect_segment_length(layer: &Layer) -> Option<f64> {
    let mut lens: Vec<f64> = Vec::new();
    for feature in layer.features.iter() {
        if let Some(geom) = &feature.geometry {
            collect_edge_lengths(geom, &mut lens);
        }
    }
    lens.retain(|l| *l > 1e-12 && l.is_finite());
    if lens.is_empty() {
        return None;
    }
    lens.sort_by(f64::total_cmp);
    Some(lens[lens.len() / 2])
}

fn collect_edge_lengths(geom: &Geometry, lens: &mut Vec<f64>) {
    let mut ring_edges = |ring: &Ring| {
        let pts = ring_points(ring);
        lens.extend((0..pts.len()).map(|i| dist(pts[i], pts[(i + 1) % pts.len()])));
    };
    match geom {
        Geometry::Polygon { exterior, interiors } => {
            ring_edges(exterior);
            interiors.iter().for_each(&mut ring_edges);
        }
        Geometry::MultiPolygon(parts) => {
            for (ext, holes) in parts {
                ring_edges(ext);
                holes.iter().for_each(&mut ring_edges);
            }
        }
        Geometry::LineString(coords) => {
            let pts = dedup_points(coords);
            lens.extend(pts.windows(2).map(|w| dist(w[0], w[1])));
        }
        Geometry::MultiLineString(lines) => {
            for coords in lines {
                let pts = dedup_points(coords);
                lens.extend(pts.windows(2).map(|w| dist(w[0], w[1])));
            }
        }
        Geometry::GeometryCollection(geoms) => {
            geoms.iter().for_each(|g| collect_edge_lengths(g, lens));
        }
        _ => {}
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, FieldValue, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext { progress: &NullProgress, capabilities: &AllowAllCapabilities }
    }

    /// Densifies the ring `corners` every ~`step` units, jittering each vertex
    /// perpendicular to its edge by an alternating ±`amp` zigzag (a stand-in
    /// for raster staircase noise).
    fn zigzag_ring(corners: &[(f64, f64)], step: f64, amp: f64) -> Vec<Coord> {
        let mut out = Vec::new();
        let mut idx = 0usize;
        let n = corners.len();
        for i in 0..n {
            let (ax, ay) = corners[i];
            let (bx, by) = corners[(i + 1) % n];
            let len = (bx - ax).hypot(by - ay);
            let (nx, ny) = (-(by - ay) / len, (bx - ax) / len);
            let count = (len / step).ceil().max(1.0) as usize;
            for k in 0..count {
                let t = k as f64 / count as f64;
                let j = if idx.is_multiple_of(2) { amp } else { -amp };
                idx += 1;
                out.push(Coord::xy(ax + t * (bx - ax) + j * nx, ay + t * (by - ay) + j * ny));
            }
        }
        out
    }

    fn layer_with_geometry(geom: Geometry) -> String {
        let mut layer = Layer::new("features");
        layer.add_field(FieldDef::new("name", FieldType::Text));
        layer
            .add_feature(Some(geom), &[("name", FieldValue::Text("feat".into()))])
            .unwrap();
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SmoothNaturalFeaturesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn geometry_of(layer: &Layer, idx: usize) -> &Geometry {
        layer.features[idx].geometry.as_ref().unwrap()
    }

    /// Sharpest turn (degrees) between consecutive edge directions.
    fn max_turn_degrees(pts: &[P], closed: bool) -> f64 {
        let n = pts.len();
        let edges: Vec<(f64, f64)> = (0..if closed { n } else { n - 1 })
            .map(|i| (pts[(i + 1) % n].x - pts[i].x, pts[(i + 1) % n].y - pts[i].y))
            .collect();
        let m = edges.len();
        let turns = if closed { m } else { m - 1 };
        (0..turns)
            .map(|i| {
                let (ax, ay) = edges[i];
                let (bx, by) = edges[(i + 1) % m];
                (ax * by - ay * bx).atan2(ax * bx + ay * by).abs().to_degrees()
            })
            .fold(0.0, f64::max)
    }

    fn max_dist_to_rect(pts: &[P], w: f64, h: f64) -> f64 {
        let corners = [
            P { x: 0.0, y: 0.0 },
            P { x: w, y: 0.0 },
            P { x: w, y: h },
            P { x: 0.0, y: h },
        ];
        pts.iter()
            .map(|p| {
                (0..4)
                    .map(|i| point_seg_dist(*p, corners[i], corners[(i + 1) % 4]))
                    .fold(f64::INFINITY, f64::min)
            })
            .fold(0.0, f64::max)
    }

    #[test]
    fn smooths_a_zigzag_polygon_and_preserves_area() {
        let corners = [(0.0, 0.0), (20.0, 0.0), (20.0, 10.0), (0.0, 10.0)];
        let input_ring = zigzag_ring(&corners, 1.0, 0.5);
        let input_area = signed_area(&dedup_points(&input_ring)).abs();
        let input = layer_with_geometry(Geometry::polygon(input_ring, vec![]));

        // segment_length is auto-detected (zigzag edges are ~1.1 units).
        let (out, layer) = run_tool(json!({ "input": input }));
        assert_eq!(out.outputs["smoothed_count"], json!(1));
        let seg = out.outputs["segment_length"].as_f64().unwrap();
        assert!((0.5..2.0).contains(&seg), "unexpected auto segment length {seg}");

        let pts = match geometry_of(&layer, 0) {
            Geometry::Polygon { exterior, .. } => ring_points(exterior),
            other => panic!("expected polygon, got {other:?}"),
        };
        assert!(pts.len() > 8, "expected a smoothed curve, got {} vertices", pts.len());
        // The zigzag (±0.5) must be gone: the boundary hugs the clean rectangle.
        let dev = max_dist_to_rect(&pts, 20.0, 10.0);
        assert!(dev <= 1.0, "smoothed boundary strays {dev} from the rectangle");
        // No sharp corners survive: the zigzag turns ~90°, the curve stays gentle.
        let turn = max_turn_degrees(&pts, true);
        assert!(turn <= 45.0, "sharpest turn is still {turn} degrees");
        // Area restored to within 0.1% (Newton targets 0.01%).
        let area = signed_area(&pts).abs();
        assert!(
            (area - input_area).abs() <= 1e-3 * input_area,
            "area {area} drifted from original {input_area}"
        );
        // Attributes pass through.
        assert_eq!(
            layer.features[0].get(&layer.schema, "name").unwrap(),
            &FieldValue::Text("feat".into())
        );
    }

    #[test]
    fn preserve_area_false_lets_corner_cutting_shrink_the_polygon() {
        let corners = [(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0)];
        let square = Geometry::polygon(
            corners.iter().map(|&(x, y)| Coord::xy(x, y)).collect(),
            vec![],
        );
        let args = |preserve: bool| {
            json!({
                "input": layer_with_geometry(square.clone()),
                "segment_length": 1.0,
                "preserve_area": preserve,
            })
        };
        let area_of = |layer: &Layer| match geometry_of(layer, 0) {
            Geometry::Polygon { exterior, .. } => signed_area(&ring_points(exterior)).abs(),
            other => panic!("expected polygon, got {other:?}"),
        };

        let (_, shrunk) = run_tool(args(false));
        let (_, restored) = run_tool(args(true));
        let (shrunk_area, restored_area) = (area_of(&shrunk), area_of(&restored));
        assert!(shrunk_area < 399.0, "corner cutting should shrink the square, got {shrunk_area}");
        assert!(
            (restored_area - 400.0).abs() <= 0.4,
            "restored area {restored_area} should be within 0.1% of 400"
        );
    }

    #[test]
    fn preserves_area_of_a_polygon_with_a_hole() {
        let ext = [(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0)];
        // Hole wound opposite to the exterior, as writers emit it.
        let hole = [(6.0, 6.0), (6.0, 14.0), (14.0, 14.0), (14.0, 6.0)];
        let to_coords = |pts: &[(f64, f64)]| pts.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
        let input = layer_with_geometry(Geometry::polygon(to_coords(&ext), vec![to_coords(&hole)]));

        let (_, layer) = run_tool(json!({ "input": input, "segment_length": 1.0 }));
        let (ext_pts, holes) = match geometry_of(&layer, 0) {
            Geometry::Polygon { exterior, interiors } => (
                ring_points(exterior),
                interiors.iter().map(ring_points).collect::<Vec<_>>(),
            ),
            other => panic!("expected polygon, got {other:?}"),
        };
        assert_eq!(holes.len(), 1, "the hole must survive smoothing");
        let area = part_area(&ext_pts, &holes);
        let target = 400.0 - 64.0;
        assert!(
            (area - target).abs() <= 1e-3 * target,
            "part area {area} drifted from original {target}"
        );
    }

    #[test]
    fn smooths_a_zigzag_line_keeping_its_endpoints() {
        let mut coords = vec![Coord::xy(0.0, 0.0)];
        coords.extend((1..20).map(|i| Coord::xy(i as f64, if i % 2 == 0 { 0.5 } else { -0.5 })));
        coords.push(Coord::xy(20.0, 0.0));
        let input = layer_with_geometry(Geometry::line_string(coords));

        let (out, layer) = run_tool(json!({ "input": input }));
        assert_eq!(out.outputs["smoothed_count"], json!(1));
        let pts = match geometry_of(&layer, 0) {
            Geometry::LineString(coords) => dedup_points(coords),
            other => panic!("expected linestring, got {other:?}"),
        };
        let (first, last) = (pts[0], *pts.last().unwrap());
        assert!(dist(first, P { x: 0.0, y: 0.0 }) <= 1e-9, "start endpoint moved to {first:?}");
        assert!(dist(last, P { x: 20.0, y: 0.0 }) <= 1e-9, "end endpoint moved to {last:?}");
        let worst_y = pts.iter().map(|p| p.y.abs()).fold(0.0, f64::max);
        assert!(worst_y <= 0.5, "zigzag not damped: |y| up to {worst_y}");
        let turn = max_turn_degrees(&pts, false);
        assert!(turn <= 45.0, "sharpest turn is still {turn} degrees");
    }

    #[test]
    fn passes_points_through_unchanged() {
        let mut layer = Layer::new("mixed");
        layer.add_feature(Some(Geometry::point(3.0, 4.0)), &[]).unwrap();
        layer
            .add_feature(
                Some(Geometry::polygon(zigzag_ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)], 1.0, 0.4), vec![])),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input }));
        assert_eq!(out.outputs["smoothed_count"], json!(1));
        assert_eq!(out.outputs["skipped_count"], json!(1));
        match geometry_of(&layer, 0) {
            Geometry::Point(c) => assert_eq!((c.x, c.y), (3.0, 4.0)),
            other => panic!("point should pass through, got {other:?}"),
        }
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = SmoothNaturalFeaturesTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args).unwrap_err()
        };
        bad(json!({}));
        bad(json!({ "input": "a.geojson", "segment_length": -1.0 }));
        bad(json!({ "input": "a.geojson", "segment_length": "zero" }));
        bad(json!({ "input": "a.geojson", "iterations": 0 }));
        bad(json!({ "input": "a.geojson", "iterations": 2.5 }));
        bad(json!({ "input": "a.geojson", "iterations": 99 }));
        bad(json!({ "input": "a.geojson", "preserve_area": "banana" }));
    }
}
