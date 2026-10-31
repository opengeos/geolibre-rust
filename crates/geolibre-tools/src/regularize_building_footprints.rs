//! GeoLibre tool: regularize building footprint polygons.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Regularize Building Footprint*
//! (3D Analyst): it normalizes noisy footprint polygons — typically digitized
//! from imagery or extracted from lidar/segmentation rasters — into clean
//! shapes whose walls follow a small set of directions.
//!
//! Methods (matching the ArcGIS tool where practical):
//! - `right_angles` — every wall snaps to the footprint's dominant direction or
//!   its perpendicular (90° corners only).
//! - `right_angles_and_diagonals` — additionally allows 45° walls; the
//!   `diagonal_penalty` factor controls how strongly right angles are preferred.
//! - `any_angle` — straightens walls (removes noise vertices within the
//!   tolerance) without constraining their directions.
//! - `circle` — replaces the footprint with its least-squares best-fit circle,
//!   subject to `min_radius`/`max_radius`.
//!
//! The pipeline per polygon part: de-noise the ring (Douglas–Peucker within
//! `tolerance`; retried at half and quarter tolerance when the fit fails,
//! since a deep slit narrower than 2x the simplification tolerance would
//! otherwise collapse), then try several candidate wall directions (length-weighted
//! circular mean of edge orientations mod 90°, the longest wall, and the
//! minimum-area rotated rectangle's orientation). For each candidate: classify
//! each wall to its nearest allowed direction, merge consecutive walls of the
//! same direction, prune walls shorter than the tolerance (noise), fit one
//! line per wall run (length-weighted least squares at the fixed direction),
//! and rebuild corners by intersecting consecutive lines. The footprint's
//! minimum rotated rectangle, scaled to the original area, competes as one
//! more candidate (it usually wins for small blob-like masks). The
//! best-fitting candidate that passes the acceptance gate is the result.
//! Interior rings (courtyards) are regularized against the exterior's chosen
//! direction so holes stay aligned with the building.
//!
//! Like the ArcGIS tool, a feature that cannot be regularized within the
//! tolerance keeps its original geometry, and every feature gets a `status`
//! attribute: 0 = regularized, 1 = original geometry retained (Null for
//! non-polygon features, which pass through unchanged). The acceptance gate
//! rejects a result (status 1) when it self-intersects, flips orientation,
//! changes area by more than a factor of 2, has any corner farther than
//! 3.5x`tolerance` from the original boundary (square corners over rounded
//! mask corners overshoot by design), or leaves the original boundary farther
//! than 2.5x`tolerance` away for more than 10% of its vertices (a few
//! isolated segmentation spikes up to 4x`tolerance` are tolerated).
//! `tolerance` is in the layer's CRS units.
//!
//! Not implemented from the ArcGIS tool: the alignment-feature option and the
//! densification/precision knobs of its internal solver (our line-fit
//! formulation does not need them).

use std::collections::BTreeMap;
use std::f64::consts::{FRAC_PI_4, PI};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata,
    ToolParamSpec, ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct RegularizeBuildingFootprintsTool;

impl Tool for RegularizeBuildingFootprintsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "regularize_building_footprints",
            display_name: "Regularize Building Footprints",
            summary: "Normalize noisy building footprint polygons into regular shapes (right-angle, diagonal, any-angle, or circular walls).",
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
                    name: "method",
                    description: "Regularization method: 'right_angles' (default), 'right_angles_and_diagonals', 'any_angle', or 'circle'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Maximum distance (in CRS units) the regularized footprint may deviate from the original boundary. Default 1.0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "diagonal_penalty",
                    description: "For 'right_angles_and_diagonals': how strongly right angles are preferred over 45-degree walls (a wall becomes diagonal only when its angular deviation from a diagonal, times this factor, is smaller than its deviation from an axis). Default 1.5; larger values yield fewer diagonals.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_radius",
                    description: "For 'circle': smallest allowed fitted radius; smaller fits keep the original footprint. Default 0.1.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_radius",
                    description: "For 'circle': largest allowed fitted radius; larger fits keep the original footprint. Unlimited if omitted.",
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
        ctx.progress.info(&format!(
            "regularizing {} feature(s) with method '{}'",
            layer.len(),
            prm.method.as_str()
        ));

        // Record per-feature outcome like ArcGIS's STATUS field. If the input
        // already has an incompatible 'status' field, fall back to 'reg_status'.
        let status_name = match layer.schema.field("status") {
            None => "status",
            Some(def) if def.field_type == FieldType::Integer => "status",
            Some(_) => "reg_status",
        };
        layer.add_field(FieldDef::new(status_name, FieldType::Integer));
        let schema = layer.schema.clone();
        layer.extent = None; // geometries change; invalidate the cached bbox

        let (mut regularized, mut retained, mut skipped) = (0usize, 0usize, 0usize);
        for feature in layer.iter_mut() {
            let status = match feature
                .geometry
                .as_ref()
                .and_then(|g| regularize_geometry(g, &prm))
            {
                Some((geom, ok)) => {
                    feature.geometry = Some(geom);
                    if ok {
                        regularized += 1;
                        FieldValue::Integer(0)
                    } else {
                        retained += 1;
                        FieldValue::Integer(1)
                    }
                }
                // Non-polygon (or missing) geometry: pass through unchanged.
                None => {
                    skipped += 1;
                    FieldValue::Null
                }
            };
            feature
                .set(&schema, status_name, status)
                .map_err(|e| ToolError::Execution(format!("failed writing status field: {e}")))?;
        }
        ctx.progress.info(&format!(
            "{regularized} regularized, {retained} kept original, {skipped} non-polygon"
        ));

        let feature_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("regularized_count".to_string(), json!(regularized));
        outputs.insert("retained_count".to_string(), json!(retained));
        outputs.insert("skipped_count".to_string(), json!(skipped));
        outputs.insert("method".to_string(), json!(prm.method.as_str()));
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Method {
    RightAngles,
    RightAnglesAndDiagonals,
    AnyAngle,
    Circle,
}

impl Method {
    fn as_str(self) -> &'static str {
        match self {
            Self::RightAngles => "right_angles",
            Self::RightAnglesAndDiagonals => "right_angles_and_diagonals",
            Self::AnyAngle => "any_angle",
            Self::Circle => "circle",
        }
    }
}

struct Params {
    method: Method,
    tolerance: f64,
    diagonal_penalty: f64,
    min_radius: f64,
    max_radius: Option<f64>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let method = match parse_optional_str(args, "method")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("right_angles") => Method::RightAngles,
        Some("right_angles_and_diagonals") => Method::RightAnglesAndDiagonals,
        Some("any_angle") => Method::AnyAngle,
        Some("circle") => Method::Circle,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown method '{other}' (expected right_angles, right_angles_and_diagonals, any_angle, or circle)"
            )))
        }
    };
    let tolerance = parse_optional_f64(args, "tolerance")?.unwrap_or(1.0);
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'tolerance' must be a positive number".to_string(),
        ));
    }
    let diagonal_penalty = parse_optional_f64(args, "diagonal_penalty")?.unwrap_or(1.5);
    if !(diagonal_penalty >= 0.0 && diagonal_penalty.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'diagonal_penalty' must be a non-negative number".to_string(),
        ));
    }
    let min_radius = parse_optional_f64(args, "min_radius")?.unwrap_or(0.1);
    if !(min_radius >= 0.0 && min_radius.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'min_radius' must be a non-negative number".to_string(),
        ));
    }
    let max_radius = parse_optional_f64(args, "max_radius")?;
    if let Some(mx) = max_radius {
        if !(mx > min_radius && mx.is_finite()) {
            return Err(ToolError::Validation(
                "parameter 'max_radius' must be greater than 'min_radius'".to_string(),
            ));
        }
    }
    Ok(Params { method, tolerance, diagonal_penalty, min_radius, max_radius })
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

// ── Geometry dispatch ─────────────────────────────────────────────────────────

/// Regularizes a polygonal geometry. Returns `None` for non-polygon geometries;
/// otherwise `(new_geometry, fully_regularized)` — parts that fail keep their
/// original rings and clear the flag.
fn regularize_geometry(geom: &Geometry, prm: &Params) -> Option<(Geometry, bool)> {
    match geom {
        Geometry::Polygon { exterior, interiors } => {
            let ((exterior, interiors), ok) = regularize_part(exterior, interiors, prm);
            Some((Geometry::Polygon { exterior, interiors }, ok))
        }
        Geometry::MultiPolygon(parts) => {
            let mut all_ok = true;
            let new_parts = parts
                .iter()
                .map(|(ext, ints)| {
                    let (part, ok) = regularize_part(ext, ints, prm);
                    all_ok &= ok;
                    part
                })
                .collect();
            Some((Geometry::MultiPolygon(new_parts), all_ok))
        }
        _ => None,
    }
}

/// Regularizes one polygon part (exterior + holes). On failure the original
/// rings are returned with `false`.
fn regularize_part(exterior: &Ring, interiors: &[Ring], prm: &Params) -> ((Ring, Vec<Ring>), bool) {
    let original = || ((exterior.clone(), interiors.to_vec()), false);
    let raw = ring_points(exterior);
    if raw.len() < 3 {
        return original();
    }

    match prm.method {
        Method::Circle => match regularize_circle(&raw, prm) {
            // A circular footprint has no courtyards.
            Some(new) => ((pts_to_ring(&new), Vec::new()), true),
            None => original(),
        },
        Method::AnyAngle => {
            let simplified = rdp_ring(&raw, prm.tolerance);
            if simplified.len() >= 3 && acceptance_score(&raw, &simplified, prm.tolerance).is_some()
            {
                let holes = interiors
                    .iter()
                    .map(|hole| {
                        let hraw = ring_points(hole);
                        if hraw.len() < 3 {
                            return hole.clone();
                        }
                        let hsimp = rdp_ring(&hraw, prm.tolerance);
                        if hsimp.len() >= 3
                            && acceptance_score(&hraw, &hsimp, prm.tolerance).is_some()
                        {
                            pts_to_ring(&hsimp)
                        } else {
                            hole.clone()
                        }
                    })
                    .collect();
                ((pts_to_ring(&simplified), holes), true)
            } else {
                original()
            }
        }
        Method::RightAngles | Method::RightAnglesAndDiagonals => {
            match fit_rectilinear_exterior(&raw, prm) {
                Some((new_ext, theta)) => {
                    // Holes align to the exterior's chosen direction; a hole
                    // that cannot be regularized keeps its original shape.
                    let holes = interiors
                        .iter()
                        .map(|hole| {
                            let hraw = ring_points(hole);
                            if hraw.len() < 3 {
                                return hole.clone();
                            }
                            let hsimp = rdp_ring(&hraw, prm.tolerance);
                            if hsimp.len() < 3 {
                                return hole.clone();
                            }
                            regularize_ring_rectilinear(&hsimp, theta, prm)
                                .filter(|pts| {
                                    acceptance_score(&hraw, pts, prm.tolerance).is_some()
                                })
                                .map(|pts| pts_to_ring(&pts))
                                .unwrap_or_else(|| hole.clone())
                        })
                        .collect();
                    ((pts_to_ring(&new_ext), holes), true)
                }
                None => original(),
            }
        }
    }
}

// ── Core algorithm ────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Debug)]
struct P {
    x: f64,
    y: f64,
}

/// Extracts a ring's vertices, dropping consecutive duplicates and the closing
/// duplicate if present.
fn ring_points(ring: &Ring) -> Vec<P> {
    let mut pts: Vec<P> = Vec::with_capacity(ring.len());
    for c in ring.coords() {
        let p = P { x: c.x, y: c.y };
        if pts.last().is_none_or(|last| dist(*last, p) > 1e-12) {
            pts.push(p);
        }
    }
    while pts.len() >= 2 && dist(pts[0], *pts.last().unwrap()) <= 1e-12 {
        pts.pop();
    }
    pts
}

fn pts_to_ring(pts: &[P]) -> Ring {
    Ring::new(pts.iter().map(|p| Coord::xy(p.x, p.y)).collect())
}

fn dist(a: P, b: P) -> f64 {
    (a.x - b.x).hypot(a.y - b.y)
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

/// Minimum distance from `p` to the boundary of the ring `pts` (unclosed).
fn dist_to_ring(p: P, pts: &[P]) -> f64 {
    let n = pts.len();
    (0..n)
        .map(|i| point_seg_dist(p, pts[i], pts[(i + 1) % n]))
        .fold(f64::INFINITY, f64::min)
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

/// Dominant wall direction in radians, in (-45°, 45°]: the length-weighted
/// circular mean of edge orientations folded modulo 90° (the classic 4-alpha
/// doubling trick, so a wall and its perpendicular vote for the same answer).
fn dominant_angle(pts: &[P]) -> f64 {
    let n = pts.len();
    let (mut sx, mut sy) = (0.0f64, 0.0f64);
    for i in 0..n {
        let (a, b) = (pts[i], pts[(i + 1) % n]);
        let (dx, dy) = (b.x - a.x, b.y - a.y);
        let len = dx.hypot(dy);
        if len <= 0.0 {
            continue;
        }
        let alpha = dy.atan2(dx);
        sx += len * (4.0 * alpha).cos();
        sy += len * (4.0 * alpha).sin();
    }
    0.25 * sy.atan2(sx)
}

/// Angular distance between two undirected line orientations (mod 180°).
fn line_angle_dist(a: f64, b: f64) -> f64 {
    let d = (a - b).rem_euclid(PI);
    d.min(PI - d)
}

/// Classifies an edge orientation to a direction class `k` in 0..4, meaning
/// `theta + k*45°`: even k = axis directions, odd k = diagonals.
fn classify_edge(alpha: f64, theta: f64, allow_diagonals: bool, penalty: f64) -> usize {
    let d: Vec<f64> = (0..4)
        .map(|k| line_angle_dist(alpha, theta + k as f64 * FRAC_PI_4))
        .collect();
    let (k_axis, d_axis) = if d[0] <= d[2] { (0, d[0]) } else { (2, d[2]) };
    if allow_diagonals {
        let (k_diag, d_diag) = if d[1] <= d[3] { (1, d[1]) } else { (3, d[3]) };
        if d_diag * penalty < d_axis {
            return k_diag;
        }
    }
    k_axis
}

/// Fits a rectilinear footprint to an exterior ring, trying several candidate
/// wall directions and keeping the best result that passes the acceptance
/// gate. The footprint's minimum rotated rectangle, scaled to the original
/// area, competes as one more candidate (the right answer for small blob-like
/// masks). Douglas–Peucker cannot see a ring doubling back on itself, so a
/// deep slit narrower than 2x the simplification tolerance collapses and the
/// fit degenerates; when no candidate passes, the fit retries with the
/// simplification tolerance halved, then quartered (the acceptance gate keeps
/// judging against the user's tolerance). Returns the new ring and the wall
/// direction used, so holes can be aligned to it.
fn fit_rectilinear_exterior(raw: &[P], prm: &Params) -> Option<(Vec<P>, f64)> {
    let rect = min_area_rect(raw);
    let rect_candidate = rect.as_ref().and_then(|rect| {
        let scale = (signed_area(raw).abs() / rect.area.max(1e-300)).sqrt();
        let (cx, cy) = (
            rect.corners.iter().map(|p| p.x).sum::<f64>() / 4.0,
            rect.corners.iter().map(|p| p.y).sum::<f64>() / 4.0,
        );
        let mut pts: Vec<P> = rect
            .corners
            .iter()
            .map(|p| P { x: cx + (p.x - cx) * scale, y: cy + (p.y - cy) * scale })
            .collect();
        if signed_area(raw) < 0.0 {
            pts.reverse();
        }
        acceptance_score(raw, &pts, prm.tolerance).map(|score| (pts, rect.theta, score))
    });

    for divisor in [1.0, 2.0, 4.0] {
        let simplified = rdp_ring(raw, prm.tolerance / divisor);
        if simplified.len() < 3 {
            continue;
        }
        let mut candidates = vec![dominant_angle(&simplified)];
        if let Some(t) = longest_edge_angle(&simplified) {
            candidates.push(t);
        }
        if let Some(r) = &rect {
            candidates.push(r.theta);
        }
        // Directions equivalent mod 45° produce identical classifications.
        let mut unique: Vec<f64> = Vec::new();
        for t in candidates {
            if unique
                .iter()
                .all(|u| line_angle_dist(4.0 * t, 4.0 * u) > 4.0f64.to_radians())
            {
                unique.push(t);
            }
        }

        let mut best: Option<(Vec<P>, f64, f64)> = None; // (ring, theta, score)
        for theta in unique {
            if let Some(pts) = regularize_ring_rectilinear(&simplified, theta, prm) {
                if let Some(score) = acceptance_score(raw, &pts, prm.tolerance) {
                    if best.as_ref().is_none_or(|(_, _, s)| score < *s) {
                        best = Some((pts, theta, score));
                    }
                }
            }
        }
        // The rect candidate competes at this level too.
        if let Some((pts, theta, score)) = &rect_candidate {
            if best.as_ref().is_none_or(|(_, _, s)| score < s) {
                best = Some((pts.clone(), *theta, *score));
            }
        }
        if let Some((pts, theta, _)) = best {
            return Some((pts, theta));
        }
    }
    None
}

/// Orientation of the ring's longest edge (an alternative wall-direction
/// estimate: robust when one long true wall dominates a noisy outline).
fn longest_edge_angle(pts: &[P]) -> Option<f64> {
    let n = pts.len();
    let mut best: Option<(f64, f64)> = None; // (length, angle)
    for i in 0..n {
        let (a, b) = (pts[i], pts[(i + 1) % n]);
        let (dx, dy) = (b.x - a.x, b.y - a.y);
        let len = dx.hypot(dy);
        if best.is_none_or(|(l, _)| len > l) {
            best = Some((len, dy.atan2(dx)));
        }
    }
    best.map(|(_, a)| a)
}

/// Andrew's monotone-chain convex hull, counter-clockwise.
fn convex_hull(pts: &[P]) -> Vec<P> {
    let mut sorted: Vec<P> = pts.to_vec();
    sorted.sort_by(|a, b| a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y)));
    sorted.dedup_by(|a, b| dist(*a, *b) <= 1e-12);
    let n = sorted.len();
    if n < 3 {
        return sorted;
    }
    let cross = |o: P, a: P, b: P| (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
    let mut lower: Vec<P> = Vec::new();
    for &p in &sorted {
        while lower.len() >= 2 && cross(lower[lower.len() - 2], lower[lower.len() - 1], p) <= 0.0 {
            lower.pop();
        }
        lower.push(p);
    }
    let mut upper: Vec<P> = Vec::new();
    for &p in sorted.iter().rev() {
        while upper.len() >= 2 && cross(upper[upper.len() - 2], upper[upper.len() - 1], p) <= 0.0 {
            upper.pop();
        }
        upper.push(p);
    }
    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

struct MinRect {
    corners: [P; 4],
    theta: f64,
    area: f64,
}

/// Minimum-area rotated bounding rectangle via rotating calipers over the
/// convex hull.
fn min_area_rect(pts: &[P]) -> Option<MinRect> {
    let hull = convex_hull(pts);
    let n = hull.len();
    if n < 3 {
        return None;
    }
    let mut best: Option<MinRect> = None;
    for i in 0..n {
        let (a, b) = (hull[i], hull[(i + 1) % n]);
        let (dx, dy) = (b.x - a.x, b.y - a.y);
        let len = dx.hypot(dy);
        if len <= 0.0 {
            continue;
        }
        let (ux, uy) = (dx / len, dy / len); // edge direction
        let (mut u0, mut u1, mut v0, mut v1) =
            (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
        for p in &hull {
            let u = p.x * ux + p.y * uy;
            let v = -p.x * uy + p.y * ux;
            u0 = u0.min(u);
            u1 = u1.max(u);
            v0 = v0.min(v);
            v1 = v1.max(v);
        }
        let area = (u1 - u0) * (v1 - v0);
        if best.as_ref().is_none_or(|r| area < r.area) {
            let corner = |u: f64, v: f64| P { x: u * ux - v * uy, y: u * uy + v * ux };
            best = Some(MinRect {
                corners: [corner(u0, v0), corner(u1, v0), corner(u1, v1), corner(u0, v1)],
                theta: uy.atan2(ux),
                area,
            });
        }
    }
    best
}

/// Snaps a simplified ring to walls at `theta + k*45°`: classify each edge,
/// merge consecutive edges of the same direction into runs, fit one line per
/// run (length-weighted, direction fixed), and intersect consecutive lines to
/// rebuild the corners. Returns `None` when the ring degenerates (fewer than
/// three wall runs).
fn regularize_ring_rectilinear(pts: &[P], theta: f64, prm: &Params) -> Option<Vec<P>> {
    let n = pts.len();
    if n < 3 {
        return None;
    }
    let allow_diag = prm.method == Method::RightAnglesAndDiagonals;

    // Per-edge direction class, length, and midpoint.
    let mut classes = Vec::with_capacity(n);
    let mut weights = Vec::with_capacity(n);
    let mut midpoints = Vec::with_capacity(n);
    for i in 0..n {
        let (a, b) = (pts[i], pts[(i + 1) % n]);
        let (dx, dy) = (b.x - a.x, b.y - a.y);
        let alpha = dy.atan2(dx);
        classes.push(classify_edge(alpha, theta, allow_diag, prm.diagonal_penalty));
        weights.push(dx.hypot(dy));
        midpoints.push(P { x: (a.x + b.x) * 0.5, y: (a.y + b.y) * 0.5 });
    }

    // Start scanning at a class boundary so runs never wrap.
    let start = (0..n).find(|&i| classes[i] != classes[(i + n - 1) % n])?;

    // Merge consecutive same-class edges; accumulate the length-weighted mean
    // offset `c` of the line `normal . p = c` at the run's fixed direction.
    struct Run {
        class: usize,
        offset_sum: f64,
        weight_sum: f64,
    }
    let mut runs: Vec<Run> = Vec::new();
    for off in 0..n {
        let i = (start + off) % n;
        let phi = theta + classes[i] as f64 * FRAC_PI_4;
        let normal = (-phi.sin(), phi.cos());
        let offset = normal.0 * midpoints[i].x + normal.1 * midpoints[i].y;
        match runs.last_mut() {
            Some(run) if run.class == classes[i] => {
                run.offset_sum += weights[i] * offset;
                run.weight_sum += weights[i];
            }
            _ => runs.push(Run {
                class: classes[i],
                offset_sum: weights[i] * offset,
                weight_sum: weights[i],
            }),
        }
    }
    // Prune walls shorter than the tolerance: they are noise the simplifier
    // could not remove (a jag wider than the tolerance but insignificant as a
    // wall). Dropping one may leave equal-direction neighbours adjacent, which
    // then merge into a single wall.
    loop {
        if runs.len() <= 3 {
            break;
        }
        let (weakest, weight) = runs
            .iter()
            .enumerate()
            .map(|(i, r)| (i, r.weight_sum))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .expect("runs is non-empty");
        if weight >= prm.tolerance {
            break;
        }
        runs.remove(weakest);
        let len = runs.len();
        let next = weakest % len;
        let prev = (weakest + len - 1) % len;
        if prev != next && runs[prev].class == runs[next].class {
            let merged = runs.remove(next);
            let prev = if next < prev { prev - 1 } else { prev };
            runs[prev].offset_sum += merged.offset_sum;
            runs[prev].weight_sum += merged.weight_sum;
        }
    }
    if runs.len() < 3 {
        return None;
    }

    // One fitted line per run.
    let lines: Vec<(f64, f64, f64)> = runs
        .iter()
        .map(|run| {
            let phi = theta + run.class as f64 * FRAC_PI_4;
            (-phi.sin(), phi.cos(), run.offset_sum / run.weight_sum.max(1e-300))
        })
        .collect();

    // Corner between run i and run i+1 = intersection of their lines. Adjacent
    // runs differ by at least 45°, so the intersection is well-conditioned.
    let m = lines.len();
    let mut out: Vec<P> = Vec::with_capacity(m);
    for i in 0..m {
        let (a1, b1, c1) = lines[i];
        let (a2, b2, c2) = lines[(i + 1) % m];
        let det = a1 * b2 - a2 * b1;
        if det.abs() < 1e-12 {
            return None;
        }
        let p = P { x: (c1 * b2 - c2 * b1) / det, y: (a1 * c2 - a2 * c1) / det };
        if out.last().is_none_or(|last| dist(*last, p) > 1e-9) {
            out.push(p);
        }
    }
    while out.len() >= 2 && dist(out[0], *out.last().unwrap()) <= 1e-9 {
        out.pop();
    }
    (out.len() >= 3).then_some(out)
}

/// Accepts a regularized ring only when it stays faithful to the original:
/// same orientation, area within a factor of two, no self-intersections,
/// every new corner within 3.5x tolerance of the original boundary (a square
/// corner over a rounded mask corner necessarily overshoots it), and the
/// original boundary within 2.5x tolerance of the result for at least 90% of
/// its vertices (isolated segmentation spikes up to 4x tolerance are
/// tolerated). Returns a fit score (lower is better) on acceptance, `None`
/// on rejection.
fn acceptance_score(original: &[P], new_pts: &[P], tolerance: f64) -> Option<f64> {
    if new_pts.len() < 3 {
        return None;
    }
    let a_old = signed_area(original);
    let a_new = signed_area(new_pts);
    if a_old == 0.0 || a_new == 0.0 || a_old.signum() != a_new.signum() {
        return None;
    }
    let ratio = a_new.abs() / a_old.abs();
    if !(0.5..=2.0).contains(&ratio) {
        return None;
    }
    if self_intersects(new_pts) {
        return None;
    }
    let corner_max = new_pts
        .iter()
        .map(|p| dist_to_ring(*p, original))
        .fold(0.0f64, f64::max);
    if corner_max > 3.5 * tolerance {
        return None;
    }
    let mut devs: Vec<f64> = original.iter().map(|p| dist_to_ring(*p, new_pts)).collect();
    devs.sort_unstable_by(f64::total_cmp);
    let max = *devs.last().expect("original ring is non-empty");
    let p90 = devs[((devs.len() - 1) as f64 * 0.9).ceil() as usize];
    if p90 > 2.5 * tolerance || max > 4.0 * tolerance {
        return None;
    }
    Some(p90 + corner_max)
}

/// True when any two non-adjacent edges of the ring properly cross.
fn self_intersects(pts: &[P]) -> bool {
    let n = pts.len();
    if n < 4 {
        return false;
    }
    let cross = |o: P, a: P, b: P| (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
    for i in 0..n {
        let (a1, a2) = (pts[i], pts[(i + 1) % n]);
        for j in (i + 2)..n {
            // Skip adjacent edges (they share a vertex), including the wrap.
            if i == 0 && j == n - 1 {
                continue;
            }
            let (b1, b2) = (pts[j], pts[(j + 1) % n]);
            let d1 = cross(b1, b2, a1);
            let d2 = cross(b1, b2, a2);
            let d3 = cross(a1, a2, b1);
            let d4 = cross(a1, a2, b2);
            if d1 * d2 < 0.0 && d3 * d4 < 0.0 {
                return true;
            }
        }
    }
    false
}

// ── Circle method ─────────────────────────────────────────────────────────────

/// Least-squares (Kåsa) circle fit through the ring's vertices.
fn fit_circle(pts: &[P]) -> Option<(P, f64)> {
    let n = pts.len();
    if n < 3 {
        return None;
    }
    let inv_n = 1.0 / n as f64;
    let mx = pts.iter().map(|p| p.x).sum::<f64>() * inv_n;
    let my = pts.iter().map(|p| p.y).sum::<f64>() * inv_n;
    let (mut suu, mut suv, mut svv, mut suuu, mut svvv, mut suvv, mut svuu) =
        (0.0f64, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for p in pts {
        let (u, v) = (p.x - mx, p.y - my);
        suu += u * u;
        suv += u * v;
        svv += v * v;
        suuu += u * u * u;
        svvv += v * v * v;
        suvv += u * v * v;
        svuu += v * u * u;
    }
    let det = suu * svv - suv * suv;
    if det.abs() < 1e-12 {
        return None;
    }
    let rhs_u = 0.5 * (suuu + suvv);
    let rhs_v = 0.5 * (svvv + svuu);
    let uc = (rhs_u * svv - rhs_v * suv) / det;
    let vc = (rhs_v * suu - rhs_u * suv) / det;
    let r2 = uc * uc + vc * vc + (suu + svv) * inv_n;
    if !(r2.is_finite() && r2 > 0.0) {
        return None;
    }
    Some((P { x: mx + uc, y: my + vc }, r2.sqrt()))
}

/// Replaces a footprint ring with its best-fit circle when the fit stays
/// within 2x tolerance of every vertex and the radius satisfies the bounds.
fn regularize_circle(raw: &[P], prm: &Params) -> Option<Vec<P>> {
    let (center, radius) = fit_circle(raw)?;
    if radius < prm.min_radius {
        return None;
    }
    if let Some(max_r) = prm.max_radius {
        if radius > max_r {
            return None;
        }
    }
    if raw.iter().any(|p| (dist(*p, center) - radius).abs() > 2.0 * prm.tolerance) {
        return None;
    }

    // Segment count so the chord sagitta stays under tolerance/2.
    let frac = (prm.tolerance / (2.0 * radius)).clamp(1e-6, 0.5);
    let step = 2.0 * (1.0 - frac).acos();
    let segments = ((2.0 * PI / step).ceil() as usize).clamp(16, 256);
    // Preserve the original winding.
    let sign = if signed_area(raw) >= 0.0 { 1.0 } else { -1.0 };
    Some(
        (0..segments)
            .map(|i| {
                let ang = sign * 2.0 * PI * i as f64 / segments as f64;
                P { x: center.x + radius * ang.cos(), y: center.y + radius * ang.sin() }
            })
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext { progress: &NullProgress, capabilities: &AllowAllCapabilities }
    }

    /// Deterministic pseudo-noise in (-1, 1).
    fn noise(i: usize) -> f64 {
        ((i as f64 * 12.9898).sin() * 43758.5453).fract()
    }

    fn rotate(x: f64, y: f64, deg: f64) -> (f64, f64) {
        let (s, c) = deg.to_radians().sin_cos();
        (x * c - y * s, x * s + y * c)
    }

    /// Densifies the ring `corners` every ~`step` units, jittering each vertex
    /// perpendicular to its edge by up to `amp`.
    fn noisy_ring(corners: &[(f64, f64)], step: f64, amp: f64) -> Vec<Coord> {
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
                let j = amp * noise(idx);
                idx += 1;
                out.push(Coord::xy(ax + t * (bx - ax) + j * nx, ay + t * (by - ay) + j * ny));
            }
        }
        out
    }

    fn layer_with_polygon(exterior: Vec<Coord>) -> String {
        let mut layer = Layer::new("buildings");
        layer.add_field(FieldDef::new("name", FieldType::Text));
        layer
            .add_feature(Some(Geometry::polygon(exterior, vec![])), &[("name", FieldValue::Text("bldg".into()))])
            .unwrap();
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = RegularizeBuildingFootprintsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn exterior_pts(layer: &Layer, idx: usize) -> Vec<P> {
        match layer.features[idx].geometry.as_ref().unwrap() {
            Geometry::Polygon { exterior, .. } => ring_points(exterior),
            other => panic!("expected polygon, got {other:?}"),
        }
    }

    fn status_of(layer: &Layer, idx: usize) -> FieldValue {
        layer.features[idx].get(&layer.schema, "status").unwrap().clone()
    }

    fn assert_has_corner_near(pts: &[P], x: f64, y: f64, tol: f64) {
        let target = P { x, y };
        let best = pts.iter().map(|p| dist(*p, target)).fold(f64::INFINITY, f64::min);
        assert!(best <= tol, "no vertex within {tol} of ({x}, {y}); closest was {best}");
    }

    #[test]
    fn squares_a_noisy_rotated_rectangle() {
        let corners: Vec<(f64, f64)> =
            [(0.0, 0.0), (20.0, 0.0), (20.0, 10.0), (0.0, 10.0)]
                .iter()
                .map(|&(x, y)| rotate(x, y, 25.0))
                .collect();
        let input = layer_with_polygon(noisy_ring(&corners, 1.0, 0.15));

        let (out, layer) = run_tool(json!({ "input": input, "tolerance": 0.5 }));
        assert_eq!(out.outputs["regularized_count"], json!(1));

        let pts = exterior_pts(&layer, 0);
        assert_eq!(pts.len(), 4, "expected a clean 4-corner rectangle, got {pts:?}");
        for &(x, y) in &corners {
            assert_has_corner_near(&pts, x, y, 0.5);
        }
        // Right angles: consecutive edges must be perpendicular.
        for i in 0..4 {
            let (a, b, c) = (pts[i], pts[(i + 1) % 4], pts[(i + 2) % 4]);
            let dot = (b.x - a.x) * (c.x - b.x) + (b.y - a.y) * (c.y - b.y);
            assert!(dot.abs() < 1e-6, "corner {i} is not square (dot={dot})");
        }
        assert_eq!(status_of(&layer, 0), FieldValue::Integer(0));
        assert_eq!(
            layer.features[0].get(&layer.schema, "name").unwrap(),
            &FieldValue::Text("bldg".into())
        );
    }

    #[test]
    fn squares_a_noisy_rotated_l_shape() {
        let corners: Vec<(f64, f64)> =
            [(0.0, 0.0), (20.0, 0.0), (20.0, 10.0), (12.0, 10.0), (12.0, 16.0), (0.0, 16.0)]
                .iter()
                .map(|&(x, y)| rotate(x, y, 15.0))
                .collect();
        let input = layer_with_polygon(noisy_ring(&corners, 1.0, 0.1));

        let (_, layer) = run_tool(json!({ "input": input, "tolerance": 0.4 }));
        let pts = exterior_pts(&layer, 0);
        assert_eq!(pts.len(), 6, "expected 6 corners for an L-shape, got {pts:?}");
        for &(x, y) in &corners {
            assert_has_corner_near(&pts, x, y, 0.4);
        }
        assert_eq!(status_of(&layer, 0), FieldValue::Integer(0));
    }

    #[test]
    fn diagonals_keep_chamfered_corners_and_right_angles_reject_them() {
        // A 20x10 rectangle with all four corners chamfered by 3 units (45°).
        let corners = [
            (3.0, 0.0), (17.0, 0.0), (20.0, 3.0), (20.0, 7.0),
            (17.0, 10.0), (3.0, 10.0), (0.0, 7.0), (0.0, 3.0),
        ];
        let rotated: Vec<(f64, f64)> =
            corners.iter().map(|&(x, y)| rotate(x, y, 10.0)).collect();

        let input = layer_with_polygon(noisy_ring(&rotated, 0.5, 0.08));
        let (out, layer) = run_tool(json!({
            "input": input,
            "method": "right_angles_and_diagonals",
            "tolerance": 0.4,
        }));
        assert_eq!(out.outputs["regularized_count"], json!(1));
        let pts = exterior_pts(&layer, 0);
        assert_eq!(pts.len(), 8, "expected 8 corners with chamfers kept, got {pts:?}");
        for &(x, y) in &rotated {
            assert_has_corner_near(&pts, x, y, 0.5);
        }

        // Right-angles-only cannot express the 45° chamfers within tolerance,
        // so the original geometry must be retained (status 1).
        let input2 = layer_with_polygon(noisy_ring(&rotated, 0.5, 0.08));
        let (out2, layer2) = run_tool(json!({
            "input": input2,
            "method": "right_angles",
            "tolerance": 0.4,
        }));
        assert_eq!(out2.outputs["retained_count"], json!(1));
        assert_eq!(status_of(&layer2, 0), FieldValue::Integer(1));
    }

    #[test]
    fn any_angle_straightens_without_snapping_directions() {
        // A triangle: nothing rectilinear about it, but any_angle should still
        // strip the noise vertices down to (roughly) the three corners.
        let corners = [(0.0, 0.0), (18.0, 2.0), (7.0, 14.0)];
        let input = layer_with_polygon(noisy_ring(&corners, 1.0, 0.1));
        let (_, layer) = run_tool(json!({
            "input": input,
            "method": "any_angle",
            "tolerance": 0.5,
        }));
        let pts = exterior_pts(&layer, 0);
        assert!(pts.len() <= 5, "expected heavy simplification, got {} vertices", pts.len());
        for &(x, y) in &corners {
            assert_has_corner_near(&pts, x, y, 0.5);
        }
        assert_eq!(status_of(&layer, 0), FieldValue::Integer(0));
    }

    #[test]
    fn fits_a_circle_and_respects_radius_bounds() {
        let (cx, cy, r) = (5.0, -3.0, 8.0);
        let ring: Vec<Coord> = (0..40)
            .map(|i| {
                let ang = 2.0 * PI * i as f64 / 40.0;
                let rr = r + 0.1 * noise(i);
                Coord::xy(cx + rr * ang.cos(), cy + rr * ang.sin())
            })
            .collect();

        let input = layer_with_polygon(ring.clone());
        let (out, layer) = run_tool(json!({ "input": input, "method": "circle", "tolerance": 0.5 }));
        assert_eq!(out.outputs["regularized_count"], json!(1));
        let pts = exterior_pts(&layer, 0);
        assert!(pts.len() >= 16);
        for p in &pts {
            let d = (dist(*p, P { x: cx, y: cy }) - r).abs();
            assert!(d < 0.1, "output vertex {d} off the fitted circle");
        }

        // A min_radius larger than the footprint keeps the original.
        let input2 = layer_with_polygon(ring);
        let (out2, layer2) = run_tool(json!({
            "input": input2,
            "method": "circle",
            "tolerance": 0.5,
            "min_radius": 20,
        }));
        assert_eq!(out2.outputs["retained_count"], json!(1));
        assert_eq!(exterior_pts(&layer2, 0).len(), 40);
    }

    #[test]
    fn regularizes_holes_against_the_exterior_direction() {
        let outer: Vec<(f64, f64)> =
            [(0.0, 0.0), (30.0, 0.0), (30.0, 20.0), (0.0, 20.0)]
                .iter()
                .map(|&(x, y)| rotate(x, y, 20.0))
                .collect();
        // Courtyard, wound opposite to the exterior.
        let inner: Vec<(f64, f64)> =
            [(10.0, 5.0), (10.0, 12.0), (20.0, 12.0), (20.0, 5.0)]
                .iter()
                .map(|&(x, y)| rotate(x, y, 20.0))
                .collect();
        let mut layer = Layer::new("buildings");
        let ext = noisy_ring(&outer, 1.0, 0.1);
        let hole = noisy_ring(&inner, 1.0, 0.1);
        layer
            .add_feature(Some(Geometry::polygon(ext, vec![hole])), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (_, out_layer) = run_tool(json!({ "input": input, "tolerance": 0.4 }));
        match out_layer.features[0].geometry.as_ref().unwrap() {
            Geometry::Polygon { exterior, interiors } => {
                assert_eq!(ring_points(exterior).len(), 4);
                assert_eq!(interiors.len(), 1);
                let hole_pts = ring_points(&interiors[0]);
                assert_eq!(hole_pts.len(), 4, "hole should regularize to 4 corners");
                for &(x, y) in &inner {
                    assert_has_corner_near(&hole_pts, x, y, 0.4);
                }
            }
            other => panic!("expected polygon, got {other:?}"),
        }
        assert_eq!(status_of(&out_layer, 0), FieldValue::Integer(0));
    }

    #[test]
    fn deep_narrow_notch_survives_regularization() {
        // Real-world regression (NAIP masks, Spokane): an E-shaped building
        // whose notch is only ~1 px wide. Douglas–Peucker at the full
        // tolerance collapses the slit (every slit vertex is within tolerance
        // of a chord down its middle), so the fit must retry at a finer
        // simplification instead of giving up.
        let coords: &[(f64, f64)] = &[
            (7.8, 21.0), (7.8, 20.4), (2.4, 20.4), (2.4, 19.8), (1.2, 19.8), (1.2, 19.2),
            (0.6, 19.2), (0.6, 17.4), (0.0, 17.4), (0.0, 12.0), (0.6, 12.0), (0.6, 10.2),
            (0.0, 10.2), (0.0, 5.4), (0.6, 5.4), (0.6, 2.4), (1.2, 2.4), (1.2, 0.6),
            (1.8, 0.6), (1.8, 0.0), (16.8, 0.0), (16.8, 0.6), (17.4, 0.6), (17.4, 7.2),
            (16.8, 7.2), (16.8, 8.4), (16.2, 8.4), (16.2, 9.0), (15.0, 9.0), (15.0, 9.6),
            (13.8, 9.6), (13.8, 10.2), (6.6, 10.2), (6.6, 10.8), (6.0, 10.8), (6.0, 11.4),
            (7.8, 11.4), (7.8, 12.0), (15.0, 12.0), (15.0, 12.6), (16.2, 12.6), (16.2, 13.2),
            (16.8, 13.2), (16.8, 19.8), (16.2, 19.8), (16.2, 20.4), (15.0, 20.4), (15.0, 21.0),
            (9.0, 21.0), (9.0, 20.4), (8.4, 20.4), (8.4, 21.0),
        ];
        let ring: Vec<Coord> = coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
        let orig_area = {
            let pts: Vec<P> = coords.iter().map(|&(x, y)| P { x, y }).collect();
            signed_area(&pts).abs()
        };
        let input = layer_with_polygon(ring);
        let (out, layer) = run_tool(json!({
            "input": input,
            "method": "right_angles_and_diagonals",
            "tolerance": 1.0,
        }));
        assert_eq!(out.outputs["regularized_count"], json!(1));
        let pts = exterior_pts(&layer, 0);
        let new_area = signed_area(&pts).abs();
        assert!(
            (new_area / orig_area - 1.0).abs() < 0.15,
            "area not preserved: {orig_area} -> {new_area}"
        );
        // The notch must survive: some vertex well inside the bbox interior.
        assert!(
            pts.iter().any(|p| p.x < 9.0 && p.x > 3.0 && p.y > 8.0 && p.y < 13.0),
            "notch collapsed: {pts:?}"
        );
    }

    #[test]
    fn min_area_rect_recovers_a_rotated_rectangle() {
        let pts: Vec<P> = [(0.0, 0.0), (20.0, 0.0), (20.0, 10.0), (0.0, 10.0)]
            .iter()
            .map(|&(x, y)| {
                let (rx, ry) = rotate(x, y, 33.0);
                P { x: rx, y: ry }
            })
            .collect();
        let rect = min_area_rect(&pts).unwrap();
        assert!((rect.area - 200.0).abs() < 1e-6, "area {}", rect.area);
        let d = line_angle_dist(rect.theta, 33f64.to_radians());
        // Orientation is meaningful mod 90 deg (either rectangle axis wins).
        assert!(
            d < 1e-6 || (d - std::f64::consts::FRAC_PI_2).abs() < 1e-6,
            "theta off by {d}"
        );
    }

    #[test]
    fn blob_masks_fall_back_to_an_area_matched_rectangle() {
        // A rounded elliptical blob: with no straight walls, the fitter (or
        // the min-rotated-rect fallback) should still produce a clean 4-corner
        // footprint of roughly the same area rather than giving up.
        let ring: Vec<Coord> = (0..24)
            .map(|i| {
                let ang = 2.0 * PI * i as f64 / 24.0;
                let (x, y) = (5.0 * ang.cos() * (1.0 + 0.05 * noise(i)),
                              3.5 * ang.sin() * (1.0 + 0.05 * noise(i + 100)));
                let (rx, ry) = rotate(x, y, 33.0);
                Coord::xy(rx, ry)
            })
            .collect();
        let orig_area = {
            let pts: Vec<P> = ring.iter().map(|c| P { x: c.x, y: c.y }).collect();
            signed_area(&pts).abs()
        };
        let input = layer_with_polygon(ring);
        let (out, layer) = run_tool(json!({ "input": input, "tolerance": 1.0 }));
        assert_eq!(out.outputs["regularized_count"], json!(1));
        let pts = exterior_pts(&layer, 0);
        assert_eq!(pts.len(), 4, "expected the min-rect fallback, got {pts:?}");
        let new_area = signed_area(&pts).abs();
        assert!(
            (new_area / orig_area - 1.0).abs() < 0.15,
            "area not preserved: {orig_area} -> {new_area}"
        );
    }

    #[test]
    fn passes_non_polygons_through_unchanged() {
        let mut layer = Layer::new("mixed");
        layer.add_feature(Some(Geometry::point(1.0, 2.0)), &[]).unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, out_layer) = run_tool(json!({ "input": input }));
        assert_eq!(out.outputs["skipped_count"], json!(1));
        assert_eq!(
            out_layer.features[0].geometry,
            Some(Geometry::point(1.0, 2.0))
        );
        assert_eq!(status_of(&out_layer, 0), FieldValue::Null);
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = RegularizeBuildingFootprintsTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input must fail");
        assert!(bad(json!({ "input": "x.geojson", "method": "bogus" })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "tolerance": 0 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "tolerance": -1.0 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "min_radius": 5, "max_radius": 2 })).is_err());
        // Numeric strings are accepted (form-posted values).
        assert!(bad(json!({ "input": "x.geojson", "tolerance": "0.5" })).is_ok());
    }
}
