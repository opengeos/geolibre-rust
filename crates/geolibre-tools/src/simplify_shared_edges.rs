//! GeoLibre tool: coverage-safe polygon simplification.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Simplify Shared Edges* (Cartography).
//! Ordinary per-feature simplification (`simplify_features`) tears a polygon
//! *coverage* apart: each polygon is simplified on its own, so a boundary shared
//! with a neighbor is simplified twice — independently — and the two results
//! diverge, opening gaps and slivers. This tool simplifies the coverage as a
//! whole so every shared boundary stays coincident.
//!
//! The approach mirrors GEOS 3.12's `CoverageSimplifier` (which geolibre-rust
//! cannot use — it is GEOS-free by design):
//!
//! 1. Build a planar **arc–node topology** from every polygon boundary. A
//!    *node* is a vertex where the topology branches (its undirected-edge degree
//!    is not 2 — a shared edge meeting the coverage boundary, or three polygons
//!    meeting). An *arc* is a maximal chain of edges between two nodes; a shared
//!    boundary is one arc referenced by both neighbors, and a ring with no node
//!    (an island's boundary) is a single closed arc.
//! 2. Simplify **each arc once** (Douglas–Peucker), pinning its node endpoints.
//! 3. Reassemble every polygon from its (now-simplified) arcs.
//!
//! Because a shared arc is simplified a single time and both neighbors reference
//! the same simplified vertices, the shared boundary is byte-identical on both
//! sides — no gaps, no slivers. `simplify_boundary` controls whether arcs on the
//! coverage's outer edge (degree-1, belonging to a single polygon) are
//! simplified too. `snap_tolerance` quantizes vertices onto a grid before
//! building the topology, so a coverage whose shared vertices are only *nearly*
//! coincident (a common artifact of imported data) still shares arcs.
//!
//! A polygon whose ring would collapse below three vertices keeps its original
//! geometry and is reported; everything else, and every non-polygon feature,
//! passes through. `tolerance`/`snap_tolerance` are in the layer's CRS units.

use std::collections::{BTreeMap, HashMap, HashSet};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Geometry, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct SimplifySharedEdgesTool;

impl Tool for SimplifySharedEdgesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "simplify_shared_edges",
            display_name: "Simplify Shared Edges",
            summary: "Simplify a polygon coverage while keeping boundaries shared between adjacent polygons coincident (no gaps or slivers), like ArcGIS Simplify Shared Edges.",
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
                    name: "tolerance",
                    description: "Douglas-Peucker simplification tolerance: the maximum distance (in CRS units) a simplified arc may deviate from the original. Default 1.0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "simplify_boundary",
                    description: "Also simplify arcs on the coverage's outer boundary (edges belonging to a single polygon). Default true; set false to keep the outer boundary unchanged.",
                    required: false,
                },
                ToolParamSpec {
                    name: "snap_tolerance",
                    description: "Quantize vertices onto a grid of this size (CRS units) before building the topology, so nearly-coincident shared vertices are treated as one. Default 0 (use coordinates exactly).",
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

        let mut layer = load_input_layer(input)?;

        // Collect every polygon ring (exterior + interiors) as a canonical
        // vertex chain, remembering where it came from so the simplified rings
        // can be written back into the right feature. Non-polygon features are
        // left untouched.
        let mut rings: Vec<Vec<P>> = Vec::new();
        let mut origins: Vec<RingOrigin> = Vec::new();
        for (fidx, feature) in layer.features.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            collect_rings(geom, fidx, prm.snap_tolerance, &mut rings, &mut origins);
        }

        ctx.progress.info(&format!(
            "{} feature(s): {} polygon ring(s) to simplify",
            layer.len(),
            rings.len()
        ));

        // Build the arc-node topology and simplify each arc once.
        let topo = Topology::build(&rings);
        let arc_simplified: Vec<Vec<P>> = topo
            .arcs
            .iter()
            .map(|arc| {
                let simplify = prm.simplify_boundary || arc.shared;
                if simplify {
                    simplify_arc(&arc.pts, arc.closed, prm.tolerance)
                } else {
                    arc.pts.clone()
                }
            })
            .collect();

        // Reassemble each ring from its (shared) simplified arcs.
        let new_rings: Vec<Vec<P>> = (0..rings.len())
            .map(|ri| topo.reassemble_ring(ri, &arc_simplified))
            .collect();

        // Write the simplified rings back, grouped by feature. A ring that
        // collapsed below three vertices is dropped back to its original shape
        // (rare; only tiny polygons whose every arc degenerated).
        let mut by_feature: HashMap<usize, Vec<(RingSlot, Vec<P>)>> = HashMap::new();
        let mut collapsed = 0usize;
        for (ri, origin) in origins.iter().enumerate() {
            let simplified = &new_rings[ri];
            let ring = if simplified.len() >= 3 {
                simplified.clone()
            } else {
                collapsed += 1;
                rings[ri].clone()
            };
            by_feature
                .entry(origin.feature)
                .or_default()
                .push((origin.slot, ring));
        }

        let mut simplified_features = 0usize;
        for (fidx, parts) in by_feature {
            let Some(geom) = layer.features[fidx].geometry.as_ref() else {
                continue;
            };
            if let Some(new_geom) = rebuild_geometry(geom, &parts) {
                layer.features[fidx].geometry = Some(new_geom);
                simplified_features += 1;
            }
        }
        layer.extent = None; // geometries changed; drop the cached bbox

        let arcs_total = topo.arcs.len();
        let shared_arcs = topo.arcs.iter().filter(|a| a.shared).count();
        ctx.progress.info(&format!(
            "{arcs_total} arc(s) ({shared_arcs} shared); simplified {simplified_features} feature(s), {collapsed} ring(s) kept original"
        ));

        let feature_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("ring_count".to_string(), json!(rings.len()));
        outputs.insert("arc_count".to_string(), json!(arcs_total));
        outputs.insert("shared_arc_count".to_string(), json!(shared_arcs));
        outputs.insert("collapsed_ring_count".to_string(), json!(collapsed));
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    tolerance: f64,
    simplify_boundary: bool,
    snap_tolerance: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let tolerance = parse_optional_f64(args, "tolerance")?.unwrap_or(1.0);
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'tolerance' must be a positive number".to_string(),
        ));
    }
    let simplify_boundary = parse_optional_bool(args, "simplify_boundary")?.unwrap_or(true);
    let snap_tolerance = parse_optional_f64(args, "snap_tolerance")?.unwrap_or(0.0);
    if !(snap_tolerance >= 0.0 && snap_tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'snap_tolerance' must be a non-negative number".to_string(),
        ));
    }
    Ok(Params {
        tolerance,
        simplify_boundary,
        snap_tolerance,
    })
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

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
}

// ── Points and canonical keys ───────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
struct P {
    x: f64,
    y: f64,
}

/// Bit-exact key for a canonical (already snapped) vertex, so topology lookups
/// treat two coincident vertices as the same node.
type Key = (u64, u64);

fn key(p: P) -> Key {
    (p.x.to_bits(), p.y.to_bits())
}

/// Canonical undirected edge key (endpoints ordered), so a shared edge hashes
/// the same regardless of which polygon's winding produced it.
fn edge_key(a: P, b: P) -> (Key, Key) {
    let (ka, kb) = (key(a), key(b));
    if ka <= kb {
        (ka, kb)
    } else {
        (kb, ka)
    }
}

/// Snaps a coordinate onto a grid of `snap` (no-op when `snap` is 0), giving
/// nearly-coincident vertices an identical canonical value.
fn canonical(x: f64, y: f64, snap: f64) -> P {
    if snap > 0.0 {
        P {
            x: (x / snap).round() * snap,
            y: (y / snap).round() * snap,
        }
    } else {
        P { x, y }
    }
}

// ── Ring collection ─────────────────────────────────────────────────────────

/// Which feature and which ring slot a collected ring belongs to, so the
/// simplified result is written back into the correct place.
#[derive(Clone, Copy)]
struct RingOrigin {
    feature: usize,
    slot: RingSlot,
}

/// Address of a ring inside a (possibly multi-part) polygon geometry.
#[derive(Clone, Copy, PartialEq, Eq)]
enum RingSlot {
    /// Exterior of polygon part `part` (0 for a single Polygon).
    Exterior { part: usize },
    /// Interior ring `hole` of polygon part `part`.
    Interior { part: usize, hole: usize },
}

fn collect_rings(
    geom: &Geometry,
    feature: usize,
    snap: f64,
    rings: &mut Vec<Vec<P>>,
    origins: &mut Vec<RingOrigin>,
) {
    let push =
        |ring: &Ring, slot: RingSlot, rings: &mut Vec<Vec<P>>, origins: &mut Vec<RingOrigin>| {
            let pts = ring_points(ring, snap);
            if pts.len() >= 3 {
                rings.push(pts);
                origins.push(RingOrigin { feature, slot });
            }
        };
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            push(exterior, RingSlot::Exterior { part: 0 }, rings, origins);
            for (h, hole) in interiors.iter().enumerate() {
                push(
                    hole,
                    RingSlot::Interior { part: 0, hole: h },
                    rings,
                    origins,
                );
            }
        }
        Geometry::MultiPolygon(parts) => {
            for (p, (exterior, interiors)) in parts.iter().enumerate() {
                push(exterior, RingSlot::Exterior { part: p }, rings, origins);
                for (h, hole) in interiors.iter().enumerate() {
                    push(
                        hole,
                        RingSlot::Interior { part: p, hole: h },
                        rings,
                        origins,
                    );
                }
            }
        }
        _ => {}
    }
}

/// Extracts a ring's vertices (canonicalized), dropping consecutive duplicates
/// and the closing duplicate, so the ring is an open cyclic chain.
fn ring_points(ring: &Ring, snap: f64) -> Vec<P> {
    let mut pts: Vec<P> = Vec::with_capacity(ring.len());
    for c in ring.coords() {
        let p = canonical(c.x, c.y, snap);
        if pts.last().is_none_or(|last| key(*last) != key(p)) {
            pts.push(p);
        }
    }
    while pts.len() >= 2 && key(pts[0]) == key(*pts.last().unwrap()) {
        pts.pop();
    }
    pts
}

// ── Arc–node topology ────────────────────────────────────────────────────────

struct Arc {
    /// Vertex chain node→node (inclusive of both endpoints); for a closed arc
    /// the endpoints are the same node and the chain is stored unclosed.
    pts: Vec<P>,
    closed: bool,
    /// True when the arc is referenced by more than one ring (an interior,
    /// shared boundary) — as opposed to a coverage-outer-boundary arc.
    shared: bool,
}

/// A ring expressed as an ordered list of arc references (arc id + direction).
struct RingArcs {
    refs: Vec<(usize, bool)>, // (arc index, forward)
}

struct Topology {
    arcs: Vec<Arc>,
    ring_arcs: Vec<RingArcs>,
}

impl Topology {
    fn build(rings: &[Vec<P>]) -> Topology {
        // Distinct undirected edges and per-vertex distinct neighbors.
        let mut edges: HashSet<(Key, Key)> = HashSet::new();
        let mut neighbors: HashMap<Key, HashSet<Key>> = HashMap::new();
        let mut coord: HashMap<Key, P> = HashMap::new();
        for ring in rings {
            let n = ring.len();
            for i in 0..n {
                let (a, b) = (ring[i], ring[(i + 1) % n]);
                coord.entry(key(a)).or_insert(a);
                coord.entry(key(b)).or_insert(b);
                if key(a) == key(b) {
                    continue;
                }
                edges.insert(edge_key(a, b));
                neighbors.entry(key(a)).or_default().insert(key(b));
                neighbors.entry(key(b)).or_default().insert(key(a));
            }
        }

        // A vertex is a node (junction) when its undirected-edge degree != 2.
        // For a valid coverage this marks exactly the topological branch points.
        let is_node = |k: &Key| neighbors.get(k).map(|s| s.len()).unwrap_or(0) != 2;

        // How many distinct rings reference each undirected edge, so an arc can
        // be flagged shared (interior) vs boundary.
        let mut edge_ring_count: HashMap<(Key, Key), usize> = HashMap::new();
        for ring in rings {
            let n = ring.len();
            let mut seen: HashSet<(Key, Key)> = HashSet::new();
            for i in 0..n {
                let e = edge_key(ring[i], ring[(i + 1) % n]);
                if e.0 != e.1 && seen.insert(e) {
                    *edge_ring_count.entry(e).or_default() += 1;
                }
            }
        }

        // Extract arcs by walking degree-2 chains between nodes; leftover
        // all-degree-2 components are closed-loop arcs.
        let mut arcs: Vec<Arc> = Vec::new();
        let mut edge_arc: HashMap<(Key, Key), usize> = HashMap::new();
        let mut used: HashSet<(Key, Key)> = HashSet::new();

        let neighbor_pts = |k: &Key| -> Vec<P> {
            neighbors
                .get(k)
                .map(|s| s.iter().map(|nk| coord[nk]).collect())
                .unwrap_or_default()
        };

        // Arcs anchored at nodes.
        let node_keys: Vec<Key> = coord.keys().copied().filter(|k| is_node(k)).collect();
        for nk in &node_keys {
            let start = coord[nk];
            for nbr in neighbor_pts(nk) {
                let e0 = edge_key(start, nbr);
                if used.contains(&e0) {
                    continue;
                }
                let arc = walk_arc(start, nbr, &is_node, &neighbors, &coord, &mut used);
                record_arc(arc, false, &edge_ring_count, &mut arcs, &mut edge_arc);
            }
        }
        // Closed loops: any edge not yet consumed lies on a node-free ring.
        for e in &edges {
            if used.contains(e) {
                continue;
            }
            let a = coord[&e.0];
            let b = coord[&e.1];
            let arc = walk_arc(a, b, &is_node, &neighbors, &coord, &mut used);
            record_arc(arc, true, &edge_ring_count, &mut arcs, &mut edge_arc);
        }

        // Map every ring to its ordered arc references.
        let node_set: HashSet<Key> = node_keys.into_iter().collect();
        let ring_arcs = rings
            .iter()
            .map(|ring| decompose_ring(ring, &node_set, &edge_arc, &arcs))
            .collect();

        Topology { arcs, ring_arcs }
    }

    /// Rebuilds ring `ri` by concatenating its arcs' simplified vertex chains,
    /// following each arc in the ring's original direction and dropping the
    /// shared node shared between consecutive arcs.
    fn reassemble_ring(&self, ri: usize, arc_simplified: &[Vec<P>]) -> Vec<P> {
        let mut out: Vec<P> = Vec::new();
        for &(aid, forward) in &self.ring_arcs[ri].refs {
            let pts = &arc_simplified[aid];
            let mut seq: Vec<P> = if forward {
                pts.clone()
            } else {
                pts.iter().rev().copied().collect()
            };
            if out.is_empty() {
                out.append(&mut seq);
            } else {
                // Drop the leading node; it repeats the previous arc's tail.
                out.extend(seq.into_iter().skip(1));
            }
        }
        while out.len() >= 2 && key(out[0]) == key(*out.last().unwrap()) {
            out.pop();
        }
        out
    }
}

/// Walks a degree-2 chain from node `start` into `first`, stopping at the next
/// node (or back at `start` for a closed loop). Marks traversed edges used.
fn walk_arc(
    start: P,
    first: P,
    is_node: &impl Fn(&Key) -> bool,
    neighbors: &HashMap<Key, HashSet<Key>>,
    coord: &HashMap<Key, P>,
    used: &mut HashSet<(Key, Key)>,
) -> Vec<P> {
    let mut arc = vec![start, first];
    used.insert(edge_key(start, first));
    let mut prev = start;
    let mut cur = first;
    while !is_node(&key(cur)) {
        // Degree-2: exactly two neighbors; step to the one we did not come from.
        let nbrs = &neighbors[&key(cur)];
        let next_key = nbrs.iter().find(|&&nk| nk != key(prev));
        let Some(&nk) = next_key else { break };
        let next = coord[&nk];
        let e = edge_key(cur, next);
        if used.contains(&e) {
            break; // closed loop returned to the start edge
        }
        used.insert(e);
        arc.push(next);
        prev = cur;
        cur = next;
        if key(cur) == key(start) {
            break; // closed loop
        }
    }
    arc
}

/// Stores an arc and indexes each of its edges back to it. `shared` is derived
/// from how many rings reference the arc's edges.
fn record_arc(
    pts: Vec<P>,
    closed: bool,
    edge_ring_count: &HashMap<(Key, Key), usize>,
    arcs: &mut Vec<Arc>,
    edge_arc: &mut HashMap<(Key, Key), usize>,
) {
    if pts.len() < 2 {
        return;
    }
    let id = arcs.len();
    let mut shared = false;
    let seg_count = if closed { pts.len() } else { pts.len() - 1 };
    for i in 0..seg_count {
        let e = edge_key(pts[i], pts[(i + 1) % pts.len()]);
        edge_arc.insert(e, id);
        if edge_ring_count.get(&e).copied().unwrap_or(0) > 1 {
            shared = true;
        }
    }
    arcs.push(Arc {
        pts,
        closed,
        shared,
    });
}

/// Expresses a ring as an ordered list of (arc, direction). A ring with no node
/// is a single closed arc.
fn decompose_ring(
    ring: &[P],
    nodes: &HashSet<Key>,
    edge_arc: &HashMap<(Key, Key), usize>,
    arcs: &[Arc],
) -> RingArcs {
    let n = ring.len();
    let start = (0..n).find(|&i| nodes.contains(&key(ring[i])));
    let Some(s0) = start else {
        // Closed loop: one arc covering the whole ring.
        let e = edge_key(ring[0], ring[1]);
        let refs = edge_arc
            .get(&e)
            .map(|&aid| vec![(aid, closed_forward(ring, &arcs[aid]))])
            .unwrap_or_default();
        return RingArcs { refs };
    };

    let mut refs = Vec::new();
    let mut i = s0;
    let mut guard = 0;
    loop {
        let e = edge_key(ring[i], ring[(i + 1) % n]);
        let Some(&aid) = edge_arc.get(&e) else { break };
        let arc = &arcs[aid];
        let forward = key(ring[i]) == key(arc.pts[0]);
        refs.push((aid, forward));
        // The arc spans (len-1) ring edges to the next node.
        i = (i + arc.pts.len().saturating_sub(1)) % n;
        guard += 1;
        if i == s0 || guard > n {
            break;
        }
    }
    RingArcs { refs }
}

/// For a closed-loop ring, whether the ring runs in the arc's stored order.
fn closed_forward(ring: &[P], arc: &Arc) -> bool {
    if arc.pts.len() < 2 || ring.len() < 2 {
        return true;
    }
    // Find where the ring starts within the arc and compare the next step.
    let start_k = key(ring[0]);
    if let Some(pos) = arc.pts.iter().position(|p| key(*p) == start_k) {
        let next_in_arc = key(arc.pts[(pos + 1) % arc.pts.len()]);
        return key(ring[1 % ring.len()]) == next_in_arc;
    }
    true
}

// ── Simplification ──────────────────────────────────────────────────────────

/// Simplifies one arc. Open arcs keep their node endpoints; closed loops keep a
/// single pinned start vertex so both referencing rings stay identical.
fn simplify_arc(pts: &[P], closed: bool, tol: f64) -> Vec<P> {
    if closed {
        rdp_closed(pts, tol)
    } else {
        rdp(pts, tol)
    }
}

/// Douglas–Peucker on an open polyline; the endpoints are always kept.
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

/// Douglas–Peucker for a closed loop: pin the start vertex (and the farthest
/// vertex from it, a stable second anchor) and simplify the two halves. Keeping
/// the start pinned means both rings that reference this loop simplify to the
/// same vertices.
fn rdp_closed(pts: &[P], tol: f64) -> Vec<P> {
    let n = pts.len();
    if n < 4 {
        return pts.to_vec();
    }
    // Anchor the start vertex and the vertex farthest from it (a stable second
    // pin), then simplify the two halves as open polylines and stitch. Pinning
    // the start keeps both rings that reference this loop identical.
    let far = (1..n)
        .max_by(|&a, &b| dist(pts[0], pts[a]).total_cmp(&dist(pts[0], pts[b])))
        .unwrap_or(n / 2);
    let first_half: Vec<P> = pts[0..=far].to_vec();
    let mut second_half: Vec<P> = pts[far..n].to_vec();
    second_half.push(pts[0]); // close the loop back to the start
    let a = rdp(&first_half, tol); // [start .. far]
    let b = rdp(&second_half, tol); // [far .. start]
                                    // a already ends at `far` and b already starts at `far` and ends at the
                                    // start duplicate; keep only b's interior so the loop stays unclosed.
    let mut out = a;
    if b.len() > 2 {
        out.extend_from_slice(&b[1..b.len() - 1]);
    }
    if out.len() < 3 {
        pts.to_vec()
    } else {
        out
    }
}

fn dist(a: P, b: P) -> f64 {
    (a.x - b.x).hypot(a.y - b.y)
}

fn point_seg_dist(p: P, a: P, b: P) -> f64 {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let len2 = dx * dx + dy * dy;
    if len2 <= 0.0 {
        return dist(p, a);
    }
    let t = (((p.x - a.x) * dx + (p.y - a.y) * dy) / len2).clamp(0.0, 1.0);
    dist(
        p,
        P {
            x: a.x + t * dx,
            y: a.y + t * dy,
        },
    )
}

// ── Geometry rebuild ─────────────────────────────────────────────────────────

fn pts_to_ring(pts: &[P]) -> Ring {
    Ring::new(pts.iter().map(|p| Coord::xy(p.x, p.y)).collect())
}

/// Rewrites `geom`'s rings from the per-slot simplified chains in `parts`.
fn rebuild_geometry(geom: &Geometry, parts: &[(RingSlot, Vec<P>)]) -> Option<Geometry> {
    let lookup = |slot: RingSlot, parts: &[(RingSlot, Vec<P>)]| -> Option<Ring> {
        parts
            .iter()
            .find(|(s, _)| *s == slot)
            .map(|(_, pts)| pts_to_ring(pts))
    };
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            let ext =
                lookup(RingSlot::Exterior { part: 0 }, parts).unwrap_or_else(|| exterior.clone());
            let holes = interiors
                .iter()
                .enumerate()
                .map(|(h, hole)| {
                    lookup(RingSlot::Interior { part: 0, hole: h }, parts)
                        .unwrap_or_else(|| hole.clone())
                })
                .collect();
            Some(Geometry::Polygon {
                exterior: ext,
                interiors: holes,
            })
        }
        Geometry::MultiPolygon(orig_parts) => {
            let new_parts = orig_parts
                .iter()
                .enumerate()
                .map(|(p, (exterior, interiors))| {
                    let ext = lookup(RingSlot::Exterior { part: p }, parts)
                        .unwrap_or_else(|| exterior.clone());
                    let holes = interiors
                        .iter()
                        .enumerate()
                        .map(|(h, hole)| {
                            lookup(RingSlot::Interior { part: p, hole: h }, parts)
                                .unwrap_or_else(|| hole.clone())
                        })
                        .collect();
                    (ext, holes)
                })
                .collect();
            Some(Geometry::MultiPolygon(new_parts))
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn ring(coords: &[(f64, f64)]) -> Vec<Coord> {
        coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect()
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SimplifySharedEdgesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn exterior(layer: &Layer, idx: usize) -> Vec<(f64, f64)> {
        match layer.features[idx].geometry.as_ref().unwrap() {
            Geometry::Polygon { exterior, .. } => {
                exterior.coords().iter().map(|c| (c.x, c.y)).collect()
            }
            other => panic!("expected polygon, got {other:?}"),
        }
    }

    /// The shared boundary between two polygons that both carry the identical
    /// jagged edge must be simplified identically in both — no divergence.
    #[test]
    fn shared_edge_stays_coincident() {
        // A jagged vertical boundary at x≈10 with a spike the tolerance removes.
        let shared: Vec<(f64, f64)> = vec![
            (10.0, 0.0),
            (10.0, 3.0),
            (10.2, 5.0),
            (10.0, 7.0),
            (10.0, 10.0),
        ];
        // Left polygon: 0..10 with the shared edge on its right (top-to-bottom).
        let mut left = vec![(0.0, 0.0), (0.0, 10.0)];
        left.extend(shared.iter().rev().copied()); // (10,10)->...->(10,0)
                                                   // Right polygon: 10..20, its left edge following the same shared chain.
        let right_ring: Vec<(f64, f64)> = {
            let mut r = vec![(10.0, 0.0), (20.0, 0.0), (20.0, 10.0), (10.0, 10.0)];
            // then follow shared from (10,10) down to (10,0) excluding endpoints
            for &pt in shared
                .iter()
                .rev()
                .skip(1)
                .take(shared.len().saturating_sub(2))
            {
                r.push(pt);
            }
            r
        };

        let mut layer = Layer::new("cov");
        layer.add_field(FieldDef::new("name", FieldType::Text));
        layer
            .add_feature(
                Some(Geometry::polygon(ring(&left), vec![])),
                &[("name", "L".into())],
            )
            .unwrap();
        layer
            .add_feature(
                Some(Geometry::polygon(ring(&right_ring), vec![])),
                &[("name", "R".into())],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input, "tolerance": 0.5 }));
        // The x=10.2 spike is within 0.5 of the straight edge, so it is removed.
        assert!(out.outputs["shared_arc_count"].as_u64().unwrap() >= 1);

        // Extract each polygon's copy of the shared boundary (the x≈10 vertices)
        // and assert they are identical point sets.
        let l = exterior(&layer, 0);
        let r = exterior(&layer, 1);
        let on_shared = |v: &Vec<(f64, f64)>| -> Vec<(u64, u64)> {
            let mut s: Vec<(u64, u64)> = v
                .iter()
                .filter(|(x, _)| (*x - 10.0).abs() < 1.0)
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
        // And the spike vertex (10.2, 5.0) must be gone from both.
        assert!(!l.iter().any(|(x, _)| (*x - 10.2).abs() < 1e-9));
        assert!(!r.iter().any(|(x, _)| (*x - 10.2).abs() < 1e-9));
    }

    /// A vertex where three polygons meet is a node and must be preserved.
    #[test]
    fn triple_junction_is_preserved() {
        // Three polygons meeting at (10,10): quadrant-ish layout.
        let a = ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]);
        let b = ring(&[(10.0, 0.0), (20.0, 0.0), (20.0, 10.0), (10.0, 10.0)]);
        let c = ring(&[
            (0.0, 10.0),
            (10.0, 10.0),
            (20.0, 10.0),
            (20.0, 20.0),
            (0.0, 20.0),
        ]);
        let mut layer = Layer::new("cov");
        for (i, r) in [a, b, c].into_iter().enumerate() {
            layer
                .add_feature(Some(Geometry::polygon(r, vec![])), &[])
                .unwrap();
            let _ = i;
        }
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (_out, layer) = run_tool(json!({ "input": input, "tolerance": 2.0 }));
        // (10,10) is a triple junction; it must survive in every polygon touching it.
        for idx in 0..3 {
            let ext = exterior(&layer, idx);
            assert!(
                ext.iter()
                    .any(|(x, y)| (*x - 10.0).abs() < 1e-9 && (*y - 10.0).abs() < 1e-9),
                "polygon {idx} lost the triple junction (10,10): {ext:?}"
            );
        }
    }

    /// simplify_boundary=false keeps the coverage's outer boundary vertices.
    #[test]
    fn boundary_preserved_when_requested() {
        // Two polygons; the outer boundary has a removable jog on the left edge.
        let left = ring(&[
            (0.0, 0.0),
            (0.1, 5.0),
            (0.0, 10.0),
            (10.0, 10.0),
            (10.0, 0.0),
        ]);
        let right = ring(&[(10.0, 0.0), (10.0, 10.0), (20.0, 10.0), (20.0, 0.0)]);
        let build = || {
            let mut layer = Layer::new("cov");
            layer
                .add_feature(Some(Geometry::polygon(left.clone(), vec![])), &[])
                .unwrap();
            layer
                .add_feature(Some(Geometry::polygon(right.clone(), vec![])), &[])
                .unwrap();
            let id = memory_store::put_vector(layer);
            memory_store::make_vector_memory_path(&id)
        };

        let (_o1, kept) =
            run_tool(json!({ "input": build(), "tolerance": 0.5, "simplify_boundary": false }));
        assert!(
            exterior(&kept, 0)
                .iter()
                .any(|(x, _)| (*x - 0.1).abs() < 1e-9),
            "boundary jog was removed despite simplify_boundary=false"
        );
        let (_o2, simplified) =
            run_tool(json!({ "input": build(), "tolerance": 0.5, "simplify_boundary": true }));
        assert!(
            !exterior(&simplified, 0)
                .iter()
                .any(|(x, _)| (*x - 0.1).abs() < 1e-9),
            "boundary jog should be removed when simplify_boundary=true"
        );
    }

    /// snap_tolerance merges near-coincident shared vertices so the arc is
    /// detected as shared and simplified once.
    #[test]
    fn snap_tolerance_merges_near_coincident_vertices() {
        // Two polygons whose shared edge vertices differ by 0.001.
        let left = ring(&[(0.0, 0.0), (0.0, 10.0), (5.0, 10.0), (5.0, 5.0), (5.0, 0.0)]);
        let right = ring(&[
            (5.001, 0.0),
            (5.001, 5.0),
            (5.001, 10.0),
            (10.0, 10.0),
            (10.0, 0.0),
        ]);
        let mut layer = Layer::new("cov");
        layer
            .add_feature(Some(Geometry::polygon(left, vec![])), &[])
            .unwrap();
        layer
            .add_feature(Some(Geometry::polygon(right, vec![])), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        // Without snapping the shared edge is not recognized (0 shared arcs)...
        let (o0, _) = run_tool(json!({ "input": input, "tolerance": 0.5 }));
        assert_eq!(o0.outputs["shared_arc_count"], json!(0));

        // ...with a 0.01 snap grid the two edges collapse to one shared arc.
        let id2 = {
            let left = ring(&[(0.0, 0.0), (0.0, 10.0), (5.0, 10.0), (5.0, 5.0), (5.0, 0.0)]);
            let right = ring(&[
                (5.001, 0.0),
                (5.001, 5.0),
                (5.001, 10.0),
                (10.0, 10.0),
                (10.0, 0.0),
            ]);
            let mut layer = Layer::new("cov");
            layer
                .add_feature(Some(Geometry::polygon(left, vec![])), &[])
                .unwrap();
            layer
                .add_feature(Some(Geometry::polygon(right, vec![])), &[])
                .unwrap();
            memory_store::put_vector(layer)
        };
        let input2 = memory_store::make_vector_memory_path(&id2);
        let (o1, _) =
            run_tool(json!({ "input": input2, "tolerance": 0.5, "snap_tolerance": 0.01 }));
        assert!(o1.outputs["shared_arc_count"].as_u64().unwrap() >= 1);
    }

    #[test]
    fn passes_non_polygons_through() {
        let mut layer = Layer::new("mixed");
        layer
            .add_feature(Some(Geometry::point(1.0, 2.0)), &[])
            .unwrap();
        layer
            .add_feature(
                Some(Geometry::polygon(
                    ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]),
                    vec![],
                )),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (_out, layer) = run_tool(json!({ "input": input, "tolerance": 1.0 }));
        assert_eq!(layer.features[0].geometry, Some(Geometry::point(1.0, 2.0)));
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = SimplifySharedEdgesTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "x.geojson", "tolerance": 0 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "tolerance": -1 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "snap_tolerance": -1 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "tolerance": 1.0 })).is_ok());
        assert!(bad(json!({ "input": "x.geojson", "simplify_boundary": "false" })).is_ok());
    }
}
