//! GeoLibre tool: sequence-aware map matching of GPS tracks to a network.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Snap Tracks* (GeoAnalytics). The
//! bundled `snap_points_to_network` snaps each point independently to the
//! nearest edge — with parallel carriageways or a dense grid that produces
//! physically impossible tracks that hop between streets. This does true
//! **sequence-aware** map matching: each fix is assigned to a network edge
//! trading off snap distance against route continuity, so the matched track
//! follows one plausible path. Continues the movement arc of `reconstruct_tracks`
//! and pairs with `thin_road_network`.
//!
//! Per track (grouped by `track_field`, time-sorted), a Viterbi dynamic program
//! runs over the candidate edges of each fix:
//!
//! * **emission** cost = the snap distance from the fix to the candidate edge;
//! * **transition** cost (Newson–Krumm style) = `|route_estimate − straight
//!   distance between the two fixes|`, where the route estimate is the
//!   along-network distance for candidates on the same or an adjacent (node-
//!   sharing) edge, and a large disconnect penalty otherwise — so hopping to an
//!   unconnected parallel street is strongly discouraged.
//!
//! Output is the input points snapped onto their matched edge, carrying the
//! matched `edge_id` and `snap_dist`. Use a projected CRS (distances in its
//! units).

use std::collections::{BTreeMap, HashMap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct SnapTracksTool;

impl Tool for SnapTracksTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "snap_tracks",
            display_name: "Snap Tracks",
            summary: "Map-match timestamped GPS tracks onto a road network with a Viterbi dynamic program (emission = snap distance, transition = route continuity), so matched tracks follow one plausible path instead of zig-zagging between parallel streets — like ArcGIS Snap Tracks.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input timestamped GPS point layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "network",
                    description: "Road network line layer to match onto.",
                    required: true,
                },
                ToolParamSpec {
                    name: "track_field",
                    description: "Field identifying each track/mover.",
                    required: true,
                },
                ToolParamSpec {
                    name: "time_field",
                    description: "Field holding each point's time (numeric seconds or ISO-8601).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output snapped point layer (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "search_distance",
                    description: "Candidate radius: only network edges within this distance of a fix are considered. Required.",
                    required: true,
                },
                ToolParamSpec {
                    name: "max_candidates",
                    description: "Maximum candidate edges kept per fix (default 5).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "network")?;
        require_str(args, "track_field")?;
        require_str(args, "time_field")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let network_path = require_str(args, "network")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let network = load_input_layer(network_path)?;
        let track_idx = layer.schema.field_index(&prm.track_field).ok_or_else(|| {
            ToolError::Validation(format!("track_field '{}' not found", prm.track_field))
        })?;
        let time_idx = layer.schema.field_index(&prm.time_field).ok_or_else(|| {
            ToolError::Validation(format!("time_field '{}' not found", prm.time_field))
        })?;

        // ── Build the network (edges + node adjacency + spatial grid) ─────────
        let net = Network::build(&network);
        if net.edges.is_empty() {
            return Err(ToolError::Execution(
                "network has no line edges".to_string(),
            ));
        }
        ctx.progress
            .info(&format!("network: {} edge(s)", net.edges.len()));

        // Group fixes by track.
        let mut tracks: BTreeMap<String, Vec<Fix>> = BTreeMap::new();
        for (fi, feature) in layer.features.iter().enumerate() {
            let Some((x, y)) = feature.geometry.as_ref().and_then(point_xy) else {
                continue;
            };
            let Some(t) = feature.attributes.get(time_idx).and_then(parse_time_value) else {
                continue;
            };
            let id = feature
                .attributes
                .get(track_idx)
                .map(value_string)
                .unwrap_or_default();
            tracks
                .entry(id)
                .or_default()
                .push(Fix { feat: fi, x, y, t });
        }

        // ── Output layer: input schema + matched fields ───────────────────────
        let mut out = Layer::new("snapped").with_geom_type(GeometryType::Point);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for field in layer.schema.fields() {
            out.add_field(field.clone());
        }
        out.add_field(FieldDef::new("edge_id", FieldType::Integer));
        out.add_field(FieldDef::new("snap_dist", FieldType::Float));

        let mut matched = 0usize;
        let mut unmatched = 0usize;
        for (_id, mut fixes) in tracks {
            fixes.sort_by(|a, b| a.t.total_cmp(&b.t));
            let snapped = match_track(&fixes, &net, &prm);
            for (fi, res) in fixes.iter().zip(snapped) {
                let src = &layer.features[fi.feat];
                let mut attrs = src.attributes.clone();
                match res {
                    Some((edge_id, px, py, d)) => {
                        attrs.push(FieldValue::Integer(edge_id as i64));
                        attrs.push(FieldValue::Float(d));
                        out.push(Feature {
                            fid: 0,
                            geometry: Some(Geometry::point(px, py)),
                            attributes: attrs,
                        });
                        matched += 1;
                    }
                    None => {
                        attrs.push(FieldValue::Integer(-1));
                        attrs.push(FieldValue::Float(f64::NAN));
                        out.push(Feature {
                            fid: 0,
                            geometry: src.geometry.clone(),
                            attributes: attrs,
                        });
                        unmatched += 1;
                    }
                }
            }
        }

        ctx.progress
            .info(&format!("{matched} fix(es) matched, {unmatched} unmatched"));

        let feature_count = out.len();
        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("matched".to_string(), json!(matched));
        outputs.insert("unmatched".to_string(), json!(unmatched));
        Ok(ToolRunResult { outputs })
    }
}

// ── Map matching (Viterbi) ───────────────────────────────────────────────────

struct Fix {
    feat: usize,
    x: f64,
    y: f64,
    t: f64,
}

/// A candidate match: which edge, the projection point, along-edge distance, and
/// the snap distance.
#[derive(Clone, Copy)]
struct Cand {
    edge: usize,
    px: f64,
    py: f64,
    along: f64,
    dist: f64,
}

/// Runs the per-track Viterbi. Returns, per fix, `Some((edge, x, y, dist))` or
/// `None` when no candidate was in range.
fn match_track(fixes: &[Fix], net: &Network, prm: &Params) -> Vec<Option<(usize, f64, f64, f64)>> {
    let n = fixes.len();
    // Candidates per fix.
    let cands: Vec<Vec<Cand>> = fixes
        .iter()
        .map(|f| net.candidates(f.x, f.y, prm.search_distance, prm.max_candidates))
        .collect();

    let mut result = vec![None; n];
    if n == 0 {
        return result;
    }

    // Viterbi over fixes that have at least one candidate. Runs are split by gaps
    // of no candidates (each contiguous run is matched independently).
    let disconnect = prm.search_distance * 8.0;
    let mut i = 0;
    while i < n {
        if cands[i].is_empty() {
            i += 1;
            continue;
        }
        // Extend the run while candidates exist.
        let start = i;
        let mut end = i;
        while end + 1 < n && !cands[end + 1].is_empty() {
            end += 1;
        }
        // DP.
        let mut cost: Vec<Vec<f64>> = Vec::new();
        let mut back: Vec<Vec<usize>> = Vec::new();
        // First column: emission only.
        cost.push(cands[start].iter().map(|c| c.dist).collect());
        back.push(vec![usize::MAX; cands[start].len()]);
        for k in (start + 1)..=end {
            let prev = k - 1;
            let gcd = dist(fixes[prev].x, fixes[prev].y, fixes[k].x, fixes[k].y);
            let mut col = vec![f64::INFINITY; cands[k].len()];
            let mut bk = vec![0usize; cands[k].len()];
            for (si, s) in cands[k].iter().enumerate() {
                for (pi, p) in cands[prev].iter().enumerate() {
                    let route = net.route_estimate(*p, *s, disconnect);
                    let trans = (route - gcd).abs();
                    let c = cost[cost.len() - 1][pi] + trans;
                    if c < col[si] {
                        col[si] = c;
                        bk[si] = pi;
                    }
                }
                col[si] += s.dist; // emission
            }
            cost.push(col);
            back.push(bk);
        }
        // Backtrack from the best final state.
        let last = cost.len() - 1;
        let mut si = (0..cost[last].len())
            .min_by(|&a, &b| cost[last][a].total_cmp(&cost[last][b]))
            .unwrap_or(0);
        for k in (start..=end).rev() {
            let ci = k - start;
            let c = cands[k][si];
            result[k] = Some((c.edge, c.px, c.py, c.dist));
            if ci > 0 {
                si = back[ci][si];
            }
        }
        i = end + 1;
    }
    result
}

// ── Network model ────────────────────────────────────────────────────────────

struct NetEdge {
    verts: Vec<(f64, f64)>,
    cum: Vec<f64>, // cumulative length at each vertex
    start_node: NodeKey,
    end_node: NodeKey,
}

type NodeKey = (i64, i64);

struct Network {
    edges: Vec<NetEdge>,
    /// grid cell -> edges present in it (for candidate search).
    grid: HashMap<(i64, i64), Vec<usize>>,
    cell: f64,
}

impl Network {
    fn build(layer: &Layer) -> Network {
        // Choose a grid/node scale from the data extent.
        let mut span = 0.0f64;
        let mut edges: Vec<NetEdge> = Vec::new();
        for f in layer.features.iter() {
            let Some(g) = f.geometry.as_ref() else {
                continue;
            };
            for chain in line_chains(g) {
                if chain.len() < 2 {
                    continue;
                }
                let mut cum = vec![0.0; chain.len()];
                for i in 1..chain.len() {
                    cum[i] =
                        cum[i - 1] + dist(chain[i - 1].0, chain[i - 1].1, chain[i].0, chain[i].1);
                }
                let (mut minx, mut miny, mut maxx, mut maxy) = (
                    f64::INFINITY,
                    f64::INFINITY,
                    f64::NEG_INFINITY,
                    f64::NEG_INFINITY,
                );
                for &(x, y) in &chain {
                    minx = minx.min(x);
                    miny = miny.min(y);
                    maxx = maxx.max(x);
                    maxy = maxy.max(y);
                }
                span = span.max(maxx - minx).max(maxy - miny);
                edges.push(NetEdge {
                    verts: chain,
                    cum,
                    start_node: (0, 0),
                    end_node: (0, 0),
                });
            }
        }
        let node_snap = (span / 100_000.0).max(1e-6);
        let cell = (span / 100.0).max(node_snap * 10.0);
        let nkey = |x: f64, y: f64| -> NodeKey {
            (
                (x / node_snap).round() as i64,
                (y / node_snap).round() as i64,
            )
        };
        let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
        for (ei, e) in edges.iter_mut().enumerate() {
            let a = e.verts[0];
            let b = e.verts[e.verts.len() - 1];
            e.start_node = nkey(a.0, a.1);
            e.end_node = nkey(b.0, b.1);
            // Rasterize densified vertices into grid cells.
            for w in e.verts.windows(2) {
                let steps = ((dist(w[0].0, w[0].1, w[1].0, w[1].1) / cell).ceil() as usize).max(1);
                for k in 0..=steps {
                    let t = k as f64 / steps as f64;
                    let x = w[0].0 + (w[1].0 - w[0].0) * t;
                    let y = w[0].1 + (w[1].1 - w[0].1) * t;
                    let c = ((x / cell).floor() as i64, (y / cell).floor() as i64);
                    let bucket = grid.entry(c).or_default();
                    if bucket.last() != Some(&ei) {
                        bucket.push(ei);
                    }
                }
            }
        }
        Network { edges, grid, cell }
    }

    /// Nearest candidate matches within `radius`, keeping the closest `max`.
    fn candidates(&self, x: f64, y: f64, radius: f64, max: usize) -> Vec<Cand> {
        let reach = (radius / self.cell).ceil() as i64 + 1;
        let (cx, cy) = (
            (x / self.cell).floor() as i64,
            (y / self.cell).floor() as i64,
        );
        let mut seen: HashMap<usize, ()> = HashMap::new();
        let mut cands: Vec<Cand> = Vec::new();
        for dx in -reach..=reach {
            for dy in -reach..=reach {
                if let Some(bucket) = self.grid.get(&(cx + dx, cy + dy)) {
                    for &ei in bucket {
                        if seen.insert(ei, ()).is_none() {
                            if let Some((px, py, along, d)) = project_onto(&self.edges[ei], x, y) {
                                if d <= radius {
                                    cands.push(Cand {
                                        edge: ei,
                                        px,
                                        py,
                                        along,
                                        dist: d,
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
        cands.sort_by(|a, b| a.dist.total_cmp(&b.dist));
        cands.truncate(max);
        cands
    }

    /// Estimated along-network distance between two candidate projections; a big
    /// `disconnect` penalty when their edges are neither identical nor adjacent.
    fn route_estimate(&self, p: Cand, s: Cand, disconnect: f64) -> f64 {
        if p.edge == s.edge {
            return (p.along - s.along).abs();
        }
        // Adjacent if the two edges share a node.
        let ep = &self.edges[p.edge];
        let es = &self.edges[s.edge];
        let shared = [ep.start_node, ep.end_node]
            .into_iter()
            .find(|n| *n == es.start_node || *n == es.end_node);
        if let Some(node) = shared {
            let np = if node == ep.start_node {
                0.0
            } else {
                ep.cum[ep.cum.len() - 1]
            };
            let ns = if node == es.start_node {
                0.0
            } else {
                es.cum[es.cum.len() - 1]
            };
            return (p.along - np).abs() + (s.along - ns).abs();
        }
        dist(p.px, p.py, s.px, s.py) + disconnect
    }
}

/// Projects (x, y) onto an edge: returns (proj_x, proj_y, along-distance, dist).
fn project_onto(e: &NetEdge, x: f64, y: f64) -> Option<(f64, f64, f64, f64)> {
    // Quick bbox reject handled by the grid; still guard empty.
    if e.verts.len() < 2 {
        return None;
    }
    let mut best = f64::INFINITY;
    let mut out = (0.0, 0.0, 0.0);
    for (i, w) in e.verts.windows(2).enumerate() {
        let (a, b) = (w[0], w[1]);
        let (dx, dy) = (b.0 - a.0, b.1 - a.1);
        let len2 = dx * dx + dy * dy;
        let t = if len2 <= 0.0 {
            0.0
        } else {
            (((x - a.0) * dx + (y - a.1) * dy) / len2).clamp(0.0, 1.0)
        };
        let (px, py) = (a.0 + t * dx, a.1 + t * dy);
        let d = dist(x, y, px, py);
        if d < best {
            best = d;
            let seg_len = len2.sqrt();
            out = (px, py, e.cum[i] + t * seg_len);
        }
    }
    Some((out.0, out.1, out.2, best))
}

// ── Helpers ──────────────────────────────────────────────────────────────────

fn dist(x0: f64, y0: f64, x1: f64, y1: f64) -> f64 {
    (x1 - x0).hypot(y1 - y0)
}

fn line_chains(geom: &Geometry) -> Vec<Vec<(f64, f64)>> {
    let to = |cs: &[Coord]| cs.iter().map(|c| (c.x, c.y)).collect::<Vec<_>>();
    match geom {
        Geometry::LineString(cs) => vec![to(cs)],
        Geometry::MultiLineString(lines) => lines.iter().map(|l| to(l)).collect(),
        _ => Vec::new(),
    }
}

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        _ => None,
    }
}

fn value_string(fv: &FieldValue) -> String {
    if let Some(i) = fv.as_i64() {
        i.to_string()
    } else if let Some(f) = fv.as_f64() {
        format!("{f}")
    } else {
        fv.as_str().unwrap_or("").to_string()
    }
}

// Time parsing (shared shape with reconstruct_tracks).
fn parse_time_value(fv: &FieldValue) -> Option<f64> {
    if let Some(n) = fv.as_f64() {
        return Some(n);
    }
    fv.as_str().and_then(parse_iso8601_seconds)
}

fn parse_iso8601_seconds(s: &str) -> Option<f64> {
    let s = s.trim();
    if s.len() < 10 {
        return None;
    }
    let b = s.as_bytes();
    let year: i64 = s.get(0..4)?.parse().ok()?;
    if b[4] != b'-' {
        return None;
    }
    let month: i64 = s.get(5..7)?.parse().ok()?;
    let day: i64 = s.get(8..10)?.parse().ok()?;
    let (mut hh, mut mm, mut ss) = (0i64, 0i64, 0i64);
    if s.len() >= 19 && (b[10] == b'T' || b[10] == b' ') {
        hh = s.get(11..13)?.parse().ok()?;
        mm = s.get(14..16)?.parse().ok()?;
        ss = s.get(17..19)?.parse().ok()?;
    }
    Some((days_from_civil(year, month, day) * 86400 + hh * 3600 + mm * 60 + ss) as f64)
}

fn days_from_civil(y: i64, m: i64, d: i64) -> i64 {
    let y = if m <= 2 { y - 1 } else { y };
    let era = if y >= 0 { y } else { y - 399 } / 400;
    let yoe = y - era * 400;
    let doy = (153 * (if m > 2 { m - 3 } else { m + 9 }) + 2) / 5 + d - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146097 + doe - 719468
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    track_field: String,
    time_field: String,
    search_distance: f64,
    max_candidates: usize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let track_field = require_str(args, "track_field")?.to_string();
    let time_field = require_str(args, "time_field")?.to_string();
    let search_distance = opt_f64(args, "search_distance")?.ok_or_else(|| {
        ToolError::Validation("required parameter 'search_distance' is missing".to_string())
    })?;
    if !(search_distance > 0.0 && search_distance.is_finite()) {
        return Err(ToolError::Validation(
            "'search_distance' must be a positive number".to_string(),
        ));
    }
    let max_candidates = match opt_f64(args, "max_candidates")? {
        None => 5,
        Some(v) if v >= 1.0 => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "'max_candidates' must be >= 1".to_string(),
            ))
        }
    };
    Ok(Params {
        track_field,
        time_field,
        search_distance,
        max_candidates,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
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

    fn net_layer(lines: &[&[(f64, f64)]]) -> String {
        let mut l = Layer::new("net")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        for coords in lines {
            let cs = coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
            l.add_feature(Some(Geometry::line_string(cs)), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn pt_layer(pts: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("gps")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("trk", FieldType::Text));
        l.add_field(FieldDef::new("t", FieldType::Float));
        for (x, y, t) in pts {
            l.add_feature(
                Some(Geometry::point(*x, *y)),
                &[("trk", "A".into()), ("t", (*t).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SnapTracksTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn edge_ids(layer: &Layer) -> Vec<i64> {
        let idx = layer.schema.field_index("edge_id").unwrap();
        layer
            .iter()
            .map(|f| f.attributes[idx].as_i64().unwrap())
            .collect()
    }

    /// Two parallel roads; a noisy track running along road 0 must stay on road 0
    /// even where individual fixes drift closer to road 1.
    #[test]
    fn stays_on_one_of_two_parallel_roads() {
        // Road 0 at y=0, road 1 at y=8 (both x 0..100). Fixes hug y≈1 but one
        // fix jumps to y=6 (closer to road 1) — sequence must keep it on road 0.
        let net = net_layer(&[
            &[(0.0, 0.0), (100.0, 0.0)][..],
            &[(0.0, 8.0), (100.0, 8.0)][..],
        ]);
        let gps = pt_layer(&[
            (10.0, 1.0, 0.0),
            (30.0, 1.0, 1.0),
            (50.0, 6.0, 2.0), // drifts toward road 1
            (70.0, 1.0, 3.0),
            (90.0, 1.0, 4.0),
        ]);
        let (out, layer) = run(json!({
            "input": gps, "network": net, "track_field": "trk", "time_field": "t",
            "search_distance": 10.0,
        }));
        assert_eq!(out.outputs["matched"], json!(5));
        let ids = edge_ids(&layer);
        assert!(
            ids.iter().all(|&e| e == 0),
            "track should stay on road 0, got {ids:?}"
        );
    }

    /// Snapped points land on the network (y≈0 for a road at y=0).
    #[test]
    fn snaps_points_onto_the_road() {
        let net = net_layer(&[&[(0.0, 0.0), (100.0, 0.0)][..]]);
        let gps = pt_layer(&[(10.0, 3.0, 0.0), (50.0, 2.5, 1.0), (90.0, 3.5, 2.0)]);
        let (_o, layer) = run(json!({
            "input": gps, "network": net, "track_field": "trk", "time_field": "t",
            "search_distance": 10.0,
        }));
        for f in layer.iter() {
            if let Some(Geometry::Point(c)) = f.geometry.as_ref() {
                assert!(c.y.abs() < 1e-6, "snapped point off the road: y={}", c.y);
            }
        }
    }

    /// A fix beyond search_distance from any edge is left unmatched.
    #[test]
    fn far_fix_is_unmatched() {
        let net = net_layer(&[&[(0.0, 0.0), (100.0, 0.0)][..]]);
        let gps = pt_layer(&[(10.0, 1.0, 0.0), (50.0, 500.0, 1.0), (90.0, 1.0, 2.0)]);
        let (out, _l) = run(json!({
            "input": gps, "network": net, "track_field": "trk", "time_field": "t",
            "search_distance": 10.0,
        }));
        assert_eq!(out.outputs["unmatched"], json!(1));
        assert_eq!(out.outputs["matched"], json!(2));
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SnapTracksTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(
            json!({ "input": "g", "network": "n", "track_field": "t", "time_field": "tt" })
        )
        .is_err());
        assert!(bad(json!({ "input": "g", "network": "n", "track_field": "t", "time_field": "tt", "search_distance": 5 })).is_ok());
    }
}
