//! GeoLibre tool: collapse small road detail (loops/jogs) for generalization.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Collapse Road Detail* (Cartography):
//! at small map scales, small interchange loops, roundabouts, dogbones, and
//! offset jogs in a road network are graphically indistinct and should be
//! replaced by a single through-connection while **connectivity is preserved**.
//! This completes the road-generalization suite alongside `thin_road_network`
//! (density thinning) and `collapse_dual_lines_to_centerline` (paired lines).
//!
//! Algorithm — reuses the road-graph assembly of `thin_road_network` (junction
//! nodes from snapped segment endpoints, one edge per line feature):
//!
//! 1. **Self-loop roundabouts** — a road drawn as a single closed ring
//!    (start node == end node) whose bounding diameter is below
//!    `collapse_distance` is deleted; its incident spokes already meet at that
//!    junction, so removing the ring leaves a clean through-junction.
//! 2. **Small cycles** (multi-edge roundabouts, dogbones, offset jogs) — for
//!    every edge, the shortest alternate path between its endpoints closes the
//!    *smallest* cycle through it. If the bounding diameter of that cycle's
//!    boundary geometry is below `collapse_distance`, all of its junction nodes
//!    are merged into one node at the cycle's centroid; the cycle's own edges
//!    (now self-loops) are dropped and every incident road is reconnected to the
//!    centroid — a single through-connection.
//!
//! Merging nodes and dropping only cycle-internal edges never changes the number
//! of connected components, so the network stays as routable as it started.
//! Roads not touching any collapsed loop keep their exact geometry. An optional
//! `road_class_field` restricts collapse to loops whose boundary roads all share
//! the same class, so genuine junctions between different road classes survive.
//!
//! The `collapse_distance` diameter is measured in metres for a geographic
//! (lon/lat) CRS via a haversine bounding-box diagonal, and in CRS units
//! otherwise.
//!
//! Scope for v1: each loop collapses to its **centroid** with incident roads
//! snapped to that point; a full medial-axis / skeleton *face collapse* (routing
//! a curved centerline through the collapsed detail) is deferred to a later
//! revision, as is collapsing loops larger than the diameter threshold.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldValue, Geometry};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CollapseRoadDetailTool;

impl Tool for CollapseRoadDetailTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "collapse_road_detail",
            display_name: "Collapse Road Detail",
            summary: "Collapse small road loops and jogs (roundabouts, dogbones, offset jogs) below a diameter threshold to single through-connections for small-scale display, preserving network connectivity — like ArcGIS Collapse Road Detail.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input road (line) vector layer, format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "collapse_distance",
                    description: "Maximum detail diameter to collapse. Loops/jogs whose bounding diameter is below this (metres for a geographic CRS, CRS units otherwise) are collapsed; larger ones are kept.",
                    required: true,
                },
                ToolParamSpec {
                    name: "road_class_field",
                    description: "Optional road-class field. When set, only loops whose boundary roads all share the same class value are collapsed (junctions between different classes are preserved).",
                    required: false,
                },
                ToolParamSpec {
                    name: "snap_tolerance",
                    description: "Distance within which road endpoints are treated as the same junction. Default 1e-6.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let schema = layer.schema.clone();
        let geographic = layer.crs_epsg().map(|e| e == 4326).unwrap_or(true);

        // Build the road graph: one edge per line feature, snapping endpoints to
        // shared junction nodes. Non-line features carry `edge_of_feature = None`
        // and pass through unchanged.
        let inv_tol = 1.0 / prm.snap_tolerance;
        let mut node_key: BTreeMap<(i64, i64), usize> = BTreeMap::new();
        let mut positions: Vec<(f64, f64)> = Vec::new();
        let mut key = |x: f64, y: f64| -> usize {
            let k = ((x * inv_tol).round() as i64, (y * inv_tol).round() as i64);
            if let Some(&id) = node_key.get(&k) {
                id
            } else {
                let id = positions.len();
                node_key.insert(k, id);
                positions.push((x, y));
                id
            }
        };

        let mut edges: Vec<Edge> = Vec::new();
        let mut edge_of_feature: Vec<Option<usize>> = vec![None; layer.features.len()];
        for (fi, feature) in layer.features.iter().enumerate() {
            let Some(coords) = feature.geometry.as_ref().and_then(line_coords) else {
                continue; // non-line: passes through unchanged
            };
            if coords.len() < 2 {
                continue;
            }
            let a = key(coords[0].0, coords[0].1);
            let b = key(coords[coords.len() - 1].0, coords[coords.len() - 1].1);
            let class = prm.road_class_field.as_ref().and_then(|f| {
                feature
                    .get(&schema, f)
                    .ok()
                    .filter(|v| !matches!(v, FieldValue::Null))
                    .map(value_string)
            });
            edge_of_feature[fi] = Some(edges.len());
            edges.push(Edge {
                a,
                b,
                coords,
                class,
                alive: true,
            });
        }
        let n_nodes = positions.len();
        let n_roads = edges.len();

        // Union-find over junction nodes (loop collapse merges nodes).
        let mut parent: Vec<usize> = (0..n_nodes).collect();

        ctx.progress
            .info(&format!("{n_roads} road(s), {n_nodes} junction(s)"));

        let mut loops_collapsed = 0usize;

        // Pass 1: single-feature closed rings (roundabouts drawn as one way).
        for e in edges.iter_mut() {
            if e.a == e.b && bbox_diameter(&e.coords, geographic) < prm.collapse_distance {
                e.alive = false;
                loops_collapsed += 1;
            }
        }

        // Pass 2: small multi-edge cycles (roundabouts, dogbones, offset jogs).
        // Each iteration collapses at most one loop, then rebuilds and rescans.
        let max_iter = n_roads + 5;
        for _ in 0..max_iter {
            // Adjacency over current representative nodes and alive edges.
            let mut adj: Vec<Vec<(usize, usize, f64)>> = vec![Vec::new(); n_nodes];
            for (ei, e) in edges.iter().enumerate() {
                if !e.alive {
                    continue;
                }
                let (ra, rb) = (find(&mut parent, e.a), find(&mut parent, e.b));
                if ra == rb {
                    continue;
                }
                let w = e.length(geographic);
                adj[ra].push((rb, ei, w));
                adj[rb].push((ra, ei, w));
            }

            // Find the first edge that closes a small, same-class cycle.
            let mut collapse: Option<Vec<usize>> = None;
            for (ei, e) in edges.iter().enumerate() {
                if !e.alive {
                    continue;
                }
                let (ra, rb) = (find(&mut parent, e.a), find(&mut parent, e.b));
                if ra == rb {
                    continue;
                }
                let Some(path_edges) = shortest_path(&adj, ra, rb, ei) else {
                    continue;
                };
                let mut cycle: Vec<usize> = path_edges;
                cycle.push(ei);
                // Bounding diameter of the whole loop boundary geometry.
                let mut pts: Vec<(f64, f64)> = Vec::new();
                for &ce in &cycle {
                    pts.extend_from_slice(&edges[ce].coords);
                }
                if bbox_diameter(&pts, geographic) >= prm.collapse_distance {
                    continue;
                }
                if !same_class(&cycle, &edges) {
                    continue;
                }
                collapse = Some(cycle);
                break;
            }

            let Some(cycle) = collapse else {
                break; // no more small loops
            };

            // Merge every junction on the cycle into one node.
            let mut nodes: Vec<usize> = Vec::new();
            for &ce in &cycle {
                nodes.push(edges[ce].a);
                nodes.push(edges[ce].b);
            }
            let anchor = nodes[0];
            for &nd in &nodes[1..] {
                union(&mut parent, anchor, nd);
            }
            // Drop every edge that is now internal to the merged blob.
            for e in edges.iter_mut() {
                if e.alive && find(&mut parent, e.a) == find(&mut parent, e.b) {
                    e.alive = false;
                }
            }
            loops_collapsed += 1;
        }

        // Per-group centroid and size (a node in a collapsed loop has size > 1).
        let mut sum: Vec<(f64, f64)> = vec![(0.0, 0.0); n_nodes];
        let mut cnt: Vec<usize> = vec![0; n_nodes];
        for (i, &(x, y)) in positions.iter().enumerate() {
            let r = find(&mut parent, i);
            sum[r].0 += x;
            sum[r].1 += y;
            cnt[r] += 1;
        }
        let centroid = |node: usize, parent: &mut [usize]| -> Option<(f64, f64)> {
            let r = find(parent, node);
            if cnt[r] > 1 {
                Some((sum[r].0 / cnt[r] as f64, sum[r].1 / cnt[r] as f64))
            } else {
                None
            }
        };

        // Rebuild the output: drop collapsed loop edges, snap surviving roads'
        // endpoints that meet a collapsed junction to that junction's centroid.
        let mut roads_removed = 0usize;
        let mut roads_unchanged = 0usize;
        let features = std::mem::take(&mut layer.features);
        let mut out_features: Vec<Feature> = Vec::with_capacity(features.len());
        for (fi, mut feature) in features.into_iter().enumerate() {
            match edge_of_feature[fi] {
                Some(ei) if !edges[ei].alive => {
                    roads_removed += 1;
                    continue;
                }
                Some(ei) => {
                    let start = centroid(edges[ei].a, &mut parent);
                    let end = centroid(edges[ei].b, &mut parent);
                    let changed = rewrite_endpoints(&mut feature.geometry, start, end);
                    if !changed {
                        roads_unchanged += 1;
                    }
                }
                None => {
                    roads_unchanged += 1; // non-line feature, untouched
                }
            }
            out_features.push(feature);
        }
        layer.features = out_features;

        ctx.progress.info(&format!(
            "collapsed {loops_collapsed} loop(s); removed {roads_removed} road(s)"
        ));

        let feature_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("road_count".to_string(), json!(n_roads));
        outputs.insert("junction_count".to_string(), json!(n_nodes));
        outputs.insert("loops_collapsed".to_string(), json!(loops_collapsed));
        outputs.insert("roads_removed".to_string(), json!(roads_removed));
        outputs.insert("roads_unchanged".to_string(), json!(roads_unchanged));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        Ok(ToolRunResult { outputs })
    }
}

struct Edge {
    a: usize,
    b: usize,
    coords: Vec<(f64, f64)>,
    class: Option<String>,
    alive: bool,
}

impl Edge {
    fn length(&self, geographic: bool) -> f64 {
        let mut len = 0.0;
        for w in self.coords.windows(2) {
            len += distance(w[0].0, w[0].1, w[1].0, w[1].1, geographic);
        }
        len
    }
}

// ── Graph helpers ──────────────────────────────────────────────────────────────

fn find(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

fn union(parent: &mut [usize], a: usize, b: usize) {
    let (ra, rb) = (find(parent, a), find(parent, b));
    if ra != rb {
        parent[rb] = ra;
    }
}

/// Shortest path (by edge length) between rep-nodes `src` and `dst` over `adj`,
/// excluding edge index `exclude`. Returns the edge indices on the path, or None
/// if `dst` is unreachable without `exclude`.
fn shortest_path(
    adj: &[Vec<(usize, usize, f64)>],
    src: usize,
    dst: usize,
    exclude: usize,
) -> Option<Vec<usize>> {
    let n = adj.len();
    let mut dist = vec![f64::INFINITY; n];
    let mut prev_edge: Vec<Option<(usize, usize)>> = vec![None; n]; // (prev_node, edge_idx)
    let mut visited = vec![false; n];
    dist[src] = 0.0;
    // Small networks: a linear-scan Dijkstra is plenty.
    loop {
        let mut u = usize::MAX;
        let mut best = f64::INFINITY;
        for (i, &d) in dist.iter().enumerate() {
            if !visited[i] && d < best {
                best = d;
                u = i;
            }
        }
        if u == usize::MAX {
            break;
        }
        if u == dst {
            break;
        }
        visited[u] = true;
        for &(v, ei, w) in &adj[u] {
            if ei == exclude || visited[v] {
                continue;
            }
            let nd = dist[u] + w;
            if nd < dist[v] {
                dist[v] = nd;
                prev_edge[v] = Some((u, ei));
            }
        }
    }
    if !dist[dst].is_finite() {
        return None;
    }
    let mut edges_on_path = Vec::new();
    let mut cur = dst;
    while cur != src {
        let (p, ei) = prev_edge[cur]?;
        edges_on_path.push(ei);
        cur = p;
    }
    Some(edges_on_path)
}

/// True when every boundary edge of the cycle shares the same road class. Edges
/// with no class field configured are always considered matching.
fn same_class(cycle: &[usize], edges: &[Edge]) -> bool {
    let mut first: Option<&Option<String>> = None;
    for &ce in cycle {
        let c = &edges[ce].class;
        match first {
            None => first = Some(c),
            Some(f) => {
                if f != c {
                    return false;
                }
            }
        }
    }
    true
}

// ── Geometry helpers ───────────────────────────────────────────────────────────

/// Flattens a line geometry's vertices to (x, y) pairs; None for non-line.
fn line_coords(geom: &Geometry) -> Option<Vec<(f64, f64)>> {
    match geom {
        Geometry::LineString(cs) => Some(cs.iter().map(|c| (c.x, c.y)).collect()),
        Geometry::MultiLineString(parts) => {
            let mut out = Vec::new();
            for p in parts {
                out.extend(p.iter().map(|c| (c.x, c.y)));
            }
            if out.is_empty() {
                None
            } else {
                Some(out)
            }
        }
        _ => None,
    }
}

/// Diameter of the bounding box of a point set (metres for geographic CRS).
fn bbox_diameter(pts: &[(f64, f64)], geographic: bool) -> f64 {
    if pts.is_empty() {
        return 0.0;
    }
    let (mut minx, mut miny, mut maxx, mut maxy) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
    for &(x, y) in pts {
        minx = minx.min(x);
        miny = miny.min(y);
        maxx = maxx.max(x);
        maxy = maxy.max(y);
    }
    distance(minx, miny, maxx, maxy, geographic)
}

/// Moves the first and/or last vertex of a line geometry to new positions.
/// Returns true when the geometry was actually changed.
fn rewrite_endpoints(
    geom: &mut Option<Geometry>,
    start: Option<(f64, f64)>,
    end: Option<(f64, f64)>,
) -> bool {
    if start.is_none() && end.is_none() {
        return false;
    }
    let Some(g) = geom.as_mut() else {
        return false;
    };
    let mut changed = false;
    match g {
        Geometry::LineString(cs) => {
            if let (Some((x, y)), Some(f)) = (start, cs.first_mut()) {
                *f = Coord::xy(x, y);
                changed = true;
            }
            if let (Some((x, y)), Some(l)) = (end, cs.last_mut()) {
                *l = Coord::xy(x, y);
                changed = true;
            }
        }
        Geometry::MultiLineString(parts) => {
            if let Some((x, y)) = start {
                if let Some(f) = parts.iter_mut().flat_map(|p| p.first_mut()).next() {
                    *f = Coord::xy(x, y);
                    changed = true;
                }
            }
            if let Some((x, y)) = end {
                if let Some(l) = parts.iter_mut().rev().flat_map(|p| p.last_mut()).next() {
                    *l = Coord::xy(x, y);
                    changed = true;
                }
            }
        }
        _ => {}
    }
    changed
}

/// Distance between two points: haversine metres for a geographic CRS, planar
/// CRS units otherwise.
fn distance(x0: f64, y0: f64, x1: f64, y1: f64, geographic: bool) -> f64 {
    if geographic {
        haversine(y0, x0, y1, x1)
    } else {
        (x1 - x0).hypot(y1 - y0)
    }
}

fn haversine(lat0: f64, lon0: f64, lat1: f64, lon1: f64) -> f64 {
    const R: f64 = 6_371_000.0;
    let (p0, p1) = (lat0.to_radians(), lat1.to_radians());
    let dphi = (lat1 - lat0).to_radians();
    let dlmb = (lon1 - lon0).to_radians();
    let a = (dphi / 2.0).sin().powi(2) + p0.cos() * p1.cos() * (dlmb / 2.0).sin().powi(2);
    2.0 * R * a.sqrt().asin()
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

// ── Parameters ──────────────────────────────────────────────────────────────────

struct Params {
    collapse_distance: f64,
    road_class_field: Option<String>,
    snap_tolerance: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let collapse_distance = opt_f64(args, "collapse_distance")?.ok_or_else(|| {
        ToolError::Validation("missing required parameter 'collapse_distance'".to_string())
    })?;
    if !(collapse_distance > 0.0 && collapse_distance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'collapse_distance' must be a positive number".to_string(),
        ));
    }
    let road_class_field = parse_optional_str(args, "road_class_field")?.map(str::to_string);
    let snap_tolerance = opt_f64(args, "snap_tolerance")?.unwrap_or(1e-6);
    if !(snap_tolerance > 0.0 && snap_tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'snap_tolerance' must be a positive number".to_string(),
        ));
    }
    Ok(Params {
        collapse_distance,
        road_class_field,
        snap_tolerance,
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
    use wbvector::{memory_store, FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn ls(pts: &[(f64, f64)]) -> Geometry {
        Geometry::line_string(pts.iter().map(|&(x, y)| Coord::xy(x, y)).collect())
    }

    fn run(build: impl FnOnce(&mut Layer), args: serde_json::Value) -> (ToolRunResult, Layer) {
        let mut layer = Layer::new("roads")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        build(&mut layer);
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);
        let mut v = args;
        v["input"] = json!(input);
        let args: ToolArgs = serde_json::from_value(v).unwrap();
        let out = CollapseRoadDetailTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// A small square roundabout (four edges) below the threshold collapses to a
    /// single node; its two through-spokes are kept and reconnect to the centre.
    #[test]
    fn collapses_small_roundabout_to_centroid() {
        let (out, layer) = run(
            |l| {
                // 10x10 square loop centred at (5,5): four edges A-B-C-D-A.
                l.add_feature(Some(ls(&[(0.0, 0.0), (10.0, 0.0)])), &[])
                    .unwrap();
                l.add_feature(Some(ls(&[(10.0, 0.0), (10.0, 10.0)])), &[])
                    .unwrap();
                l.add_feature(Some(ls(&[(10.0, 10.0), (0.0, 10.0)])), &[])
                    .unwrap();
                l.add_feature(Some(ls(&[(0.0, 10.0), (0.0, 0.0)])), &[])
                    .unwrap();
                // Two spokes into the loop corners.
                l.add_feature(Some(ls(&[(-50.0, 0.0), (0.0, 0.0)])), &[])
                    .unwrap();
                l.add_feature(Some(ls(&[(10.0, 10.0), (60.0, 10.0)])), &[])
                    .unwrap();
            },
            // Loop bbox diagonal ~14.1 < 20; spokes are far longer.
            json!({ "collapse_distance": 20.0 }),
        );
        assert_eq!(out.outputs["loops_collapsed"], json!(1));
        assert_eq!(
            out.outputs["roads_removed"],
            json!(4),
            "four ring edges gone"
        );
        assert_eq!(layer.len(), 2, "only the two spokes remain");
        // Both surviving spokes now touch the loop centroid (5,5).
        for f in &layer.features {
            if let Some(Geometry::LineString(cs)) = &f.geometry {
                let touches = cs
                    .iter()
                    .any(|c| (c.x - 5.0).abs() < 1e-6 && (c.y - 5.0).abs() < 1e-6);
                assert!(touches, "spoke should reconnect to centroid");
            }
        }
    }

    /// A dogbone: two curved edges between the same pair of nodes. Below the
    /// threshold they collapse; through-roads reconnect (connectivity kept).
    #[test]
    fn collapses_dogbone_and_keeps_through_roads() {
        let (out, layer) = run(
            |l| {
                l.add_feature(Some(ls(&[(0.0, 0.0), (5.0, 4.0), (10.0, 0.0)])), &[])
                    .unwrap(); // upper arc
                l.add_feature(Some(ls(&[(0.0, 0.0), (5.0, -4.0), (10.0, 0.0)])), &[])
                    .unwrap(); // lower arc
                l.add_feature(Some(ls(&[(-40.0, 0.0), (0.0, 0.0)])), &[])
                    .unwrap(); // west approach
                l.add_feature(Some(ls(&[(10.0, 0.0), (50.0, 0.0)])), &[])
                    .unwrap(); // east approach
            },
            json!({ "collapse_distance": 15.0 }),
        );
        assert_eq!(out.outputs["loops_collapsed"], json!(1));
        assert_eq!(out.outputs["roads_removed"], json!(2), "both arcs removed");
        assert_eq!(layer.len(), 2, "the two approaches remain and reconnect");
    }

    /// A loop larger than the threshold is left untouched (no over-collapse).
    #[test]
    fn keeps_loop_above_threshold() {
        let (out, layer) = run(
            |l| {
                l.add_feature(Some(ls(&[(0.0, 0.0), (100.0, 0.0)])), &[])
                    .unwrap();
                l.add_feature(Some(ls(&[(100.0, 0.0), (100.0, 100.0)])), &[])
                    .unwrap();
                l.add_feature(Some(ls(&[(100.0, 100.0), (0.0, 0.0)])), &[])
                    .unwrap();
            },
            json!({ "collapse_distance": 20.0 }),
        );
        assert_eq!(out.outputs["loops_collapsed"], json!(0));
        assert_eq!(layer.len(), 3, "large loop kept intact");
        assert_eq!(out.outputs["roads_unchanged"], json!(3));
    }

    /// A single-feature closed ring (start == end) below the threshold is deleted.
    #[test]
    fn collapses_single_feature_closed_ring() {
        let (out, _layer) = run(
            |l| {
                l.add_feature(
                    Some(ls(&[
                        (0.0, 0.0),
                        (8.0, 0.0),
                        (8.0, 8.0),
                        (0.0, 8.0),
                        (0.0, 0.0),
                    ])),
                    &[],
                )
                .unwrap();
                l.add_feature(Some(ls(&[(0.0, 0.0), (-50.0, 0.0)])), &[])
                    .unwrap(); // spoke
            },
            json!({ "collapse_distance": 20.0 }),
        );
        assert_eq!(out.outputs["loops_collapsed"], json!(1));
        assert_eq!(
            out.outputs["roads_removed"],
            json!(1),
            "the ring is removed"
        );
    }

    /// A non-line feature passes through unchanged.
    #[test]
    fn passes_through_non_line() {
        let (out, layer) = run(
            |l| {
                l.add_feature(Some(Geometry::point(1.0, 1.0)), &[]).unwrap();
                l.add_feature(Some(ls(&[(0.0, 0.0), (100.0, 0.0)])), &[])
                    .unwrap();
            },
            json!({ "collapse_distance": 20.0 }),
        );
        assert_eq!(out.outputs["loops_collapsed"], json!(0));
        assert_eq!(layer.len(), 2);
    }

    /// The road-class field guards collapse: a mixed-class loop is preserved.
    #[test]
    fn class_field_preserves_mixed_class_loop() {
        let build = |l: &mut Layer| {
            l.add_field(FieldDef::new("class", FieldType::Text));
            // Triangle loop, one edge of a different class.
            l.add_feature(
                Some(ls(&[(0.0, 0.0), (10.0, 0.0)])),
                &[("class", "res".into())],
            )
            .unwrap();
            l.add_feature(
                Some(ls(&[(10.0, 0.0), (5.0, 8.0)])),
                &[("class", "res".into())],
            )
            .unwrap();
            l.add_feature(
                Some(ls(&[(5.0, 8.0), (0.0, 0.0)])),
                &[("class", "primary".into())],
            )
            .unwrap();
        };
        // With the class field the mixed loop is kept.
        let (out_guarded, _) = run(
            build,
            json!({ "collapse_distance": 30.0, "road_class_field": "class" }),
        );
        assert_eq!(out_guarded.outputs["loops_collapsed"], json!(0));
        // Without it, the small loop collapses.
        let (out_free, _) = run(build, json!({ "collapse_distance": 30.0 }));
        assert_eq!(out_free.outputs["loops_collapsed"], json!(1));
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = CollapseRoadDetailTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "input": "r.geojson" })).is_err(),
            "missing collapse_distance"
        );
        assert!(bad(json!({ "input": "r.geojson", "collapse_distance": 0 })).is_err());
        assert!(bad(json!({ "input": "r.geojson", "collapse_distance": -5 })).is_err());
        assert!(bad(json!({ "input": "r.geojson", "collapse_distance": 25 })).is_ok());
    }
}
