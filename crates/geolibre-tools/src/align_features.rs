//! GeoLibre tool: conform source edges to a nearby trusted target layer.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Align Features* (Editing). The
//! bundled suite only offers `snap_endnodes` (polyline endpoints only) and
//! `snap_points_to_network` (points only); nothing moves the *interior*
//! vertices of a line or polygon boundary onto a better-accuracy reference —
//! the everyday task of conforming parcel edges to a freshly resurveyed road
//! casing, or hydrography to an updated DEM-derived stream network. Unlike
//! `rubbersheet_features` / `edgematch_features`, this is **link-free**: there
//! is no correspondence step, just "pull anything close enough onto the
//! nearest target edge".
//!
//! Algorithm, per source vertex:
//!
//! 1. **Nearest target edge.** Target line/polygon-boundary edges are
//!    segmented and grid-indexed (cell size = `search_distance`, mirroring the
//!    bbox-bucket index in `integrate`'s `insert_on_edges`). For each source
//!    vertex the nearest point on any target segment within `search_distance`
//!    is found (optionally restricted to segments whose owning target feature's
//!    `target_match_field` equals the source feature's `match_field`, like
//!    `edgematch_features`'s attribute gate). A vertex with no candidate is
//!    left untouched.
//! 2. **Tapered blend.** The raw per-vertex "move fully onto the projection"
//!    vector is smoothed with an arc-length-windowed moving average (a
//!    triangular kernel of radius `search_distance`) along its own line or
//!    ring, exactly the *rubbersheet*-style displacement idiom but applied to
//!    a step function instead of TIN barycentric weights. This is what turns a
//!    hard cut between "moved" and "fixed" vertices into a smooth taper — no
//!    kinks at the transition.
//! 3. **Shared-edge safety.** Every distinct vertex (by exact coordinate, plus
//!    the match-field grouping) is displaced exactly once: the first ring or
//!    line to visit it computes and caches the final displacement, and every
//!    other feature that shares that vertex — the standard way adjoining
//!    polygons reference a common boundary — reuses the cached value instead
//!    of recomputing it. This is the same "process a shared edge once" idea as
//!    `simplify_shared_edges`'s undirected arc hashing, applied at the vertex
//!    level (sufficient here because the final quantity is per-vertex, not
//!    per-arc), so coincident boundaries stay coincident after alignment.
//!
//! Every source feature gets `align_max_disp` / `align_mean_disp` fields (the
//! per-feature max/mean vertex displacement actually applied). Point/MultiPoint
//! features pass through unchanged (Align Features only moves edges).
//! Distances are in the layer CRS units.

use std::collections::{BTreeMap, HashMap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct AlignFeaturesTool;

impl Tool for AlignFeaturesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "align_features",
            display_name: "Align Features",
            summary: "Conform source line/polygon edges to a nearby trusted target layer within a search distance, pulling whole edge runs (not just endpoints) with a smooth taper and coincident shared boundaries, like ArcGIS Align Features.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Source vector layer to align (lines or polygons).",
                    required: true,
                },
                ToolParamSpec {
                    name: "target",
                    description: "Trusted target layer (lines, or polygon boundaries) that source vertices are pulled toward.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "search_distance",
                    description: "Max distance (CRS units) a source vertex may be pulled toward the target. Required.",
                    required: true,
                },
                ToolParamSpec {
                    name: "match_field",
                    description: "Optional source attribute field; only target edges whose target_match_field value matches are candidates for that feature's vertices.",
                    required: false,
                },
                ToolParamSpec {
                    name: "target_match_field",
                    description: "Target attribute field compared against match_field. Defaults to match_field's name if omitted; requires match_field.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "target")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let target_layer = load_input_layer(&prm.target)?;

        let match_field_idx = match &prm.match_field {
            Some(name) => Some(layer.schema.field_index(name).ok_or_else(|| {
                ToolError::Validation(format!("match_field '{name}' not found in input schema"))
            })?),
            None => None,
        };
        // target_match_field defaults to match_field's name; only meaningful
        // when match_field is set (parse_params already rejects the reverse).
        let target_field_name: Option<&str> = prm
            .target_match_field
            .as_deref()
            .or(prm.match_field.as_deref());
        if let (Some(name), true) = (target_field_name, prm.match_field.is_some()) {
            if target_layer.schema.field_index(name).is_none() {
                return Err(ToolError::Validation(format!(
                    "target_match_field '{name}' not found in target schema"
                )));
            }
        }
        let target_field_name = if prm.match_field.is_some() {
            target_field_name
        } else {
            None
        };

        let segments = build_target_segments(&target_layer, target_field_name);
        ctx.progress.info(&format!(
            "{} target edge segment(s); search_distance={}",
            segments.len(),
            prm.search_distance
        ));
        let grid = build_grid(&segments, prm.search_distance);

        layer.add_field(FieldDef::new("align_max_disp", FieldType::Float));
        layer.add_field(FieldDef::new("align_mean_disp", FieldType::Float));
        let idx_max = layer
            .schema
            .field_index("align_max_disp")
            .expect("just added");
        let idx_mean = layer
            .schema
            .field_index("align_mean_disp")
            .expect("just added");

        let mut raw_cache: HashMap<(Key, Option<String>), (f64, f64)> = HashMap::new();
        let mut final_cache: HashMap<(Key, Option<String>), (f64, f64)> = HashMap::new();
        let mut total_disp = 0.0f64;
        let mut total_vertices = 0usize;
        let mut aligned_vertices = 0usize;
        let mut max_disp_overall = 0.0f64;
        let mut aligned_features = 0usize;

        for feature in layer.features.iter_mut() {
            let match_val: Option<String> = match match_field_idx {
                Some(idx) => feature.get_by_index(idx).and_then(field_value_to_string),
                None => None,
            };
            let Some(geom) = feature.geometry.take() else {
                continue;
            };
            let mut feature_max = 0.0f64;
            let mut feature_sum = 0.0f64;
            let mut feature_n = 0usize;
            let new_geom = align_geometry(&geom, &mut |points: &[P], closed: bool| -> Vec<P> {
                let raw: Vec<(f64, f64)> = points
                    .iter()
                    .map(|p| {
                        *raw_cache
                            .entry((key(*p), match_val.clone()))
                            .or_insert_with(|| {
                                nearest_disp(
                                    *p,
                                    match_val.as_deref(),
                                    &segments,
                                    &grid,
                                    prm.search_distance,
                                )
                            })
                    })
                    .collect();
                let smoothed = smooth_chain(points, &raw, closed, prm.search_distance);
                points
                    .iter()
                    .zip(smoothed.iter())
                    .map(|(p, &(dx, dy))| {
                        let k = (key(*p), match_val.clone());
                        let (fdx, fdy) = *final_cache.entry(k).or_insert((dx, dy));
                        let d = (fdx * fdx + fdy * fdy).sqrt();
                        feature_max = feature_max.max(d);
                        feature_sum += d;
                        feature_n += 1;
                        total_disp += d;
                        total_vertices += 1;
                        if d > 1e-9 {
                            aligned_vertices += 1;
                        }
                        max_disp_overall = max_disp_overall.max(d);
                        P {
                            x: p.x + fdx,
                            y: p.y + fdy,
                        }
                    })
                    .collect()
            });
            feature.geometry = Some(new_geom);
            let feature_mean = if feature_n > 0 {
                feature_sum / feature_n as f64
            } else {
                0.0
            };
            if feature_max > 1e-9 {
                aligned_features += 1;
            }
            feature.set_by_index(idx_max, feature_max.into());
            feature.set_by_index(idx_mean, feature_mean.into());
        }
        layer.extent = None;

        let mean_disp = if total_vertices > 0 {
            total_disp / total_vertices as f64
        } else {
            0.0
        };
        ctx.progress.info(&format!(
            "{aligned_vertices}/{total_vertices} vertex(es) aligned; mean displacement {mean_disp:.4}, max {max_disp_overall:.4}"
        ));

        let feature_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("aligned_feature_count".to_string(), json!(aligned_features));
        outputs.insert("target_segment_count".to_string(), json!(segments.len()));
        outputs.insert("vertices_aligned".to_string(), json!(aligned_vertices));
        outputs.insert("vertices_total".to_string(), json!(total_vertices));
        outputs.insert("mean_displacement".to_string(), json!(mean_disp));
        outputs.insert("max_displacement".to_string(), json!(max_disp_overall));
        Ok(ToolRunResult { outputs })
    }
}

// ── Points, keys ─────────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
struct P {
    x: f64,
    y: f64,
}

type Key = (u64, u64);

fn key(p: P) -> Key {
    (p.x.to_bits(), p.y.to_bits())
}

fn dist(a: P, b: P) -> f64 {
    (a.x - b.x).hypot(a.y - b.y)
}

fn field_value_to_string(v: &FieldValue) -> Option<String> {
    match v {
        FieldValue::Null => None,
        FieldValue::Text(s) | FieldValue::Date(s) | FieldValue::DateTime(s) => Some(s.clone()),
        FieldValue::Integer(i) => Some(i.to_string()),
        FieldValue::Float(f) => Some(f.to_string()),
        FieldValue::Boolean(b) => Some(b.to_string()),
        FieldValue::Blob(_) => None,
    }
}

// ── Target edge index ────────────────────────────────────────────────────────

struct Segment {
    a: P,
    b: P,
    match_val: Option<String>,
}

fn ring_edges(ring: &Ring, f: &mut impl FnMut(P, P)) {
    let cs = ring.coords();
    let n = cs.len();
    if n < 2 {
        return;
    }
    for i in 0..n {
        let a = P {
            x: cs[i].x,
            y: cs[i].y,
        };
        let b = P {
            x: cs[(i + 1) % n].x,
            y: cs[(i + 1) % n].y,
        };
        if key(a) != key(b) {
            f(a, b);
        }
    }
}

fn for_each_edge(geom: &Geometry, f: &mut impl FnMut(P, P)) {
    match geom {
        Geometry::LineString(cs) => {
            for w in cs.windows(2) {
                let a = P {
                    x: w[0].x,
                    y: w[0].y,
                };
                let b = P {
                    x: w[1].x,
                    y: w[1].y,
                };
                if key(a) != key(b) {
                    f(a, b);
                }
            }
        }
        Geometry::MultiLineString(lines) => {
            for l in lines {
                for w in l.windows(2) {
                    let a = P {
                        x: w[0].x,
                        y: w[0].y,
                    };
                    let b = P {
                        x: w[1].x,
                        y: w[1].y,
                    };
                    if key(a) != key(b) {
                        f(a, b);
                    }
                }
            }
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            ring_edges(exterior, f);
            for r in interiors {
                ring_edges(r, f);
            }
        }
        Geometry::MultiPolygon(parts) => {
            for (e, holes) in parts {
                ring_edges(e, f);
                for r in holes {
                    ring_edges(r, f);
                }
            }
        }
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                for_each_edge(g, f);
            }
        }
        Geometry::Point(_) | Geometry::MultiPoint(_) => {}
    }
}

fn build_target_segments(target: &Layer, target_match_field: Option<&str>) -> Vec<Segment> {
    let midx = target_match_field.and_then(|f| target.schema.field_index(f));
    let mut segs = Vec::new();
    for feature in target.features.iter() {
        let match_val = midx
            .and_then(|i| feature.get_by_index(i))
            .and_then(field_value_to_string);
        if let Some(g) = feature.geometry.as_ref() {
            for_each_edge(g, &mut |a, b| {
                segs.push(Segment {
                    a,
                    b,
                    match_val: match_val.clone(),
                });
            });
        }
    }
    segs
}

/// Grid-buckets every segment into cells its (tolerance-expanded) bbox spans,
/// so a single-cell lookup at a query point's own cell finds every candidate
/// within `cell` distance (mirrors `integrate.rs`'s `insert_on_edges` index).
fn build_grid(segments: &[Segment], cell: f64) -> HashMap<(i64, i64), Vec<usize>> {
    let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    let cell = cell.max(1e-9);
    for (i, s) in segments.iter().enumerate() {
        let minx = s.a.x.min(s.b.x) - cell;
        let maxx = s.a.x.max(s.b.x) + cell;
        let miny = s.a.y.min(s.b.y) - cell;
        let maxy = s.a.y.max(s.b.y) + cell;
        let c0x = (minx / cell).floor() as i64;
        let c1x = (maxx / cell).floor() as i64;
        let c0y = (miny / cell).floor() as i64;
        let c1y = (maxy / cell).floor() as i64;
        for gx in c0x..=c1x {
            for gy in c0y..=c1y {
                grid.entry((gx, gy)).or_default().push(i);
            }
        }
    }
    grid
}

/// Nearest point on segment `a`-`b` to `p`, clamped to the segment (t in
/// [0,1]), and the distance to it.
fn project(p: P, a: P, b: P) -> (P, f64) {
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let len2 = dx * dx + dy * dy;
    if len2 <= 1e-18 {
        return (a, dist(p, a));
    }
    let t = (((p.x - a.x) * dx + (p.y - a.y) * dy) / len2).clamp(0.0, 1.0);
    let proj = P {
        x: a.x + t * dx,
        y: a.y + t * dy,
    };
    (proj, dist(p, proj))
}

/// Raw "move fully onto the nearest target edge" displacement at `p`, or
/// `(0, 0)` if nothing matches within `search_distance`.
fn nearest_disp(
    p: P,
    match_val: Option<&str>,
    segments: &[Segment],
    grid: &HashMap<(i64, i64), Vec<usize>>,
    search_distance: f64,
) -> (f64, f64) {
    let cell = search_distance.max(1e-9);
    let cx = (p.x / cell).floor() as i64;
    let cy = (p.y / cell).floor() as i64;
    let mut best: Option<(f64, P)> = None;
    if let Some(bucket) = grid.get(&(cx, cy)) {
        for &si in bucket {
            let seg = &segments[si];
            let ok = match (match_val, seg.match_val.as_deref()) {
                (Some(a), Some(b)) => a == b,
                (Some(_), None) => false,
                (None, _) => true,
            };
            if !ok {
                continue;
            }
            let (proj, d) = project(p, seg.a, seg.b);
            if d <= search_distance && best.is_none_or(|(bd, _)| d < bd) {
                best = Some((d, proj));
            }
        }
    }
    match best {
        Some((_, proj)) => (proj.x - p.x, proj.y - p.y),
        None => (0.0, 0.0),
    }
}

// ── Moving-average taper along arc length ───────────────────────────────────

/// Triangular-kernel moving average of `raw` over an arc-length window,
/// wrapping around for closed rings so the taper crosses the seam smoothly.
fn smooth_chain(points: &[P], raw: &[(f64, f64)], closed: bool, window: f64) -> Vec<(f64, f64)> {
    let n = points.len();
    if n == 0 {
        return Vec::new();
    }
    if window <= 0.0 {
        return raw.to_vec();
    }
    if !closed {
        return smooth_open(points, raw, window);
    }

    // Unroll the ring 3x so a window straddling the seam sees both sides.
    let mut cum = vec![0.0; n];
    for i in 1..n {
        cum[i] = cum[i - 1] + dist(points[i - 1], points[i]);
    }
    let total = cum[n - 1] + dist(points[n - 1], points[0]);
    if total <= 0.0 {
        return raw.to_vec();
    }
    let ext_cum: Vec<f64> = (0..3 * n)
        .map(|k| cum[k % n] + (k / n) as f64 * total)
        .collect();
    let ext_raw: Vec<(f64, f64)> = (0..3 * n).map(|k| raw[k % n]).collect();

    let mut out = Vec::with_capacity(n);
    let mut lo = 0usize;
    let mut hi = 0usize;
    for i in n..2 * n {
        while lo < i && ext_cum[i] - ext_cum[lo] > window {
            lo += 1;
        }
        while hi + 1 < 3 * n && ext_cum[hi + 1] - ext_cum[i] <= window {
            hi += 1;
        }
        out.push(weighted_avg(&ext_cum, &ext_raw, i, lo, hi, window));
    }
    out
}

fn smooth_open(points: &[P], raw: &[(f64, f64)], window: f64) -> Vec<(f64, f64)> {
    let n = points.len();
    let mut cum = vec![0.0; n];
    for i in 1..n {
        cum[i] = cum[i - 1] + dist(points[i - 1], points[i]);
    }
    let mut out = Vec::with_capacity(n);
    let mut lo = 0usize;
    let mut hi = 0usize;
    for i in 0..n {
        while lo < i && cum[i] - cum[lo] > window {
            lo += 1;
        }
        while hi + 1 < n && cum[hi + 1] - cum[i] <= window {
            hi += 1;
        }
        out.push(weighted_avg(&cum, raw, i, lo, hi, window));
    }
    out
}

fn weighted_avg(
    cum: &[f64],
    raw: &[(f64, f64)],
    i: usize,
    lo: usize,
    hi: usize,
    window: f64,
) -> (f64, f64) {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sw = 0.0;
    for j in lo..=hi {
        let d = (cum[i] - cum[j]).abs();
        let w = 1.0 - d / window;
        if w > 0.0 {
            sx += w * raw[j].0;
            sy += w * raw[j].1;
            sw += w;
        }
    }
    if sw > 0.0 {
        (sx / sw, sy / sw)
    } else {
        (0.0, 0.0)
    }
}

// ── Geometry traversal (per line/ring chain) ────────────────────────────────

/// Applies `f` to every line/ring vertex chain of `geom` (open for lines,
/// closed for polygon rings), rebuilding the geometry from the results.
/// Point/MultiPoint geometries pass through unchanged.
fn align_geometry(geom: &Geometry, f: &mut impl FnMut(&[P], bool) -> Vec<P>) -> Geometry {
    let open = |cs: &[Coord], f: &mut dyn FnMut(&[P], bool) -> Vec<P>| -> Vec<Coord> {
        let pts: Vec<P> = cs.iter().map(|c| P { x: c.x, y: c.y }).collect();
        f(&pts, false)
            .into_iter()
            .map(|p| Coord::xy(p.x, p.y))
            .collect()
    };
    let closed = |ring: &Ring, f: &mut dyn FnMut(&[P], bool) -> Vec<P>| -> Ring {
        let pts: Vec<P> = ring.coords().iter().map(|c| P { x: c.x, y: c.y }).collect();
        Ring::new(
            f(&pts, true)
                .into_iter()
                .map(|p| Coord::xy(p.x, p.y))
                .collect(),
        )
    };
    match geom {
        Geometry::Point(c) => Geometry::Point(c.clone()),
        Geometry::MultiPoint(cs) => Geometry::MultiPoint(cs.clone()),
        Geometry::LineString(cs) => Geometry::LineString(open(cs, f)),
        Geometry::MultiLineString(lines) => {
            Geometry::MultiLineString(lines.iter().map(|l| open(l, f)).collect())
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => Geometry::Polygon {
            exterior: closed(exterior, f),
            interiors: interiors.iter().map(|r| closed(r, f)).collect(),
        },
        Geometry::MultiPolygon(parts) => Geometry::MultiPolygon(
            parts
                .iter()
                .map(|(e, holes)| (closed(e, f), holes.iter().map(|r| closed(r, f)).collect()))
                .collect(),
        ),
        Geometry::GeometryCollection(gs) => {
            Geometry::GeometryCollection(gs.iter().map(|g| align_geometry(g, f)).collect())
        }
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    target: String,
    search_distance: f64,
    match_field: Option<String>,
    target_match_field: Option<String>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let target = require_str(args, "target")?.to_string();
    let search_distance = parse_optional_f64(args, "search_distance")?.ok_or_else(|| {
        ToolError::Validation("required parameter 'search_distance' is missing".to_string())
    })?;
    if !(search_distance > 0.0 && search_distance.is_finite()) {
        return Err(ToolError::Validation(
            "'search_distance' must be a positive number".to_string(),
        ));
    }
    let match_field = parse_optional_str(args, "match_field")?.map(str::to_string);
    let target_match_field = parse_optional_str(args, "target_match_field")?.map(str::to_string);
    if target_match_field.is_some() && match_field.is_none() {
        return Err(ToolError::Validation(
            "'target_match_field' requires 'match_field'".to_string(),
        ));
    }
    Ok(Params {
        target,
        search_distance,
        match_field,
        target_match_field,
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

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn line_layer(chains: &[&[(f64, f64)]]) -> String {
        let mut l = Layer::new("l")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        for chain in chains {
            let cs: Vec<Coord> = chain.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
            l.add_feature(Some(Geometry::line_string(cs)), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn poly(coords: &[(f64, f64)]) -> Geometry {
        Geometry::polygon(
            coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect(),
            vec![],
        )
    }

    fn mixed_layer(geoms: Vec<Geometry>) -> String {
        let mut l = Layer::new("m").with_crs_epsg(3857);
        for g in geoms {
            l.add_feature(Some(g), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = AlignFeaturesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn linestring_coords(layer: &Layer, i: usize) -> Vec<(f64, f64)> {
        match layer.features[i].geometry.as_ref().unwrap() {
            Geometry::LineString(cs) => cs.iter().map(|c| (c.x, c.y)).collect(),
            other => panic!("expected line, got {other:?}"),
        }
    }

    fn exterior(layer: &Layer, i: usize) -> Vec<(f64, f64)> {
        match layer.features[i].geometry.as_ref().unwrap() {
            Geometry::Polygon { exterior, .. } => {
                exterior.coords().iter().map(|c| (c.x, c.y)).collect()
            }
            other => panic!("expected polygon, got {other:?}"),
        }
    }

    /// A source line offset by less than search_distance from a straight
    /// target line snaps onto it (distance to target ~0 afterward).
    #[test]
    fn near_line_snaps_onto_target() {
        let source = line_layer(&[&[(0.0, 2.0), (50.0, 2.0), (100.0, 2.0)]]);
        let target = line_layer(&[&[(0.0, 0.0), (100.0, 0.0)]]);
        let (out, layer) = run(json!({
            "input": source, "target": target, "search_distance": 5.0,
        }));
        assert!(out.outputs["vertices_aligned"].as_u64().unwrap() >= 3);
        for &(_x, y) in &linestring_coords(&layer, 0) {
            assert!(y.abs() < 1e-6, "vertex not snapped onto target: y={y}");
        }
    }

    /// A source vertex farther than search_distance from the target does not
    /// move at all.
    #[test]
    fn far_vertex_stays_put() {
        let source = line_layer(&[&[(0.0, 100.0), (50.0, 100.0), (100.0, 100.0)]]);
        let target = line_layer(&[&[(0.0, 0.0), (100.0, 0.0)]]);
        let (out, layer) = run(json!({
            "input": source, "target": target, "search_distance": 5.0,
        }));
        assert_eq!(out.outputs["vertices_aligned"], json!(0));
        assert_eq!(
            linestring_coords(&layer, 0),
            vec![(0.0, 100.0), (50.0, 100.0), (100.0, 100.0)]
        );
    }

    /// Two polygons sharing a boundary stay coincident on that boundary after
    /// alignment (both feature's copies of the shared vertices move
    /// identically).
    #[test]
    fn polygon_shared_edge_stays_coincident() {
        let left = poly(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]);
        let right = poly(&[(10.0, 0.0), (20.0, 0.0), (20.0, 10.0), (10.0, 10.0)]);
        let input = mixed_layer(vec![left, right]);
        let target = line_layer(&[&[(10.5, -5.0), (10.5, 15.0)]]);
        let (out, layer) = run(json!({
            "input": input, "target": target, "search_distance": 2.0,
        }));
        assert!(out.outputs["vertices_aligned"].as_u64().unwrap() >= 4);
        let l = exterior(&layer, 0);
        let r = exterior(&layer, 1);
        let on_shared = |v: &[(f64, f64)]| -> Vec<(u64, u64)> {
            let mut s: Vec<(u64, u64)> = v
                .iter()
                .filter(|(x, _)| (*x - 10.0).abs() < 2.0)
                .map(|(x, y)| (x.to_bits(), y.to_bits()))
                .collect();
            s.sort_unstable();
            s
        };
        assert_eq!(
            on_shared(&l),
            on_shared(&r),
            "shared boundary diverged: {l:?} vs {r:?}"
        );
        // And it actually moved toward the target (x ~10.5, not 10.0).
        assert!(l.iter().any(|(x, _)| (*x - 10.5).abs() < 1e-6));
    }

    /// Point geometries have no edges to align and pass through untouched.
    #[test]
    fn passes_non_edge_geometry_through() {
        let mut l = Layer::new("mixed").with_crs_epsg(3857);
        l.add_feature(Some(Geometry::point(1.0, 2.0)), &[]).unwrap();
        l.add_feature(
            Some(Geometry::line_string(vec![
                Coord::xy(0.0, 100.0),
                Coord::xy(100.0, 100.0),
            ])),
            &[],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);
        let target = line_layer(&[&[(0.0, 0.0), (100.0, 0.0)]]);
        let (_out, layer) = run(json!({
            "input": input, "target": target, "search_distance": 5.0,
        }));
        assert_eq!(layer.features[0].geometry, Some(Geometry::point(1.0, 2.0)));
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            AlignFeaturesTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err()); // no target
        assert!(bad(json!({ "input": "a.geojson", "target": "b.geojson" })).is_err()); // no search_distance
        assert!(bad(json!({
            "input": "a.geojson", "target": "b.geojson", "search_distance": 0
        }))
        .is_err());
        assert!(bad(json!({
            "input": "a.geojson", "target": "b.geojson", "search_distance": -1.0
        }))
        .is_err());
        assert!(bad(json!({
            "input": "a.geojson", "target": "b.geojson",
            "search_distance": 1.0, "target_match_field": "id",
        }))
        .is_err()); // target_match_field without match_field
        assert!(bad(json!({
            "input": "a.geojson", "target": "b.geojson", "search_distance": 1.0
        }))
        .is_ok());
    }
}
