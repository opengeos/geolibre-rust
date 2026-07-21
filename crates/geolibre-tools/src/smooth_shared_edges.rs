//! GeoLibre tool: coverage-safe polygon smoothing.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Smooth Shared Edges* (Cartography).
//! Ordinary per-feature smoothing (`smooth_natural_features`, the bundled
//! `feature_preserving_smoothing`/`smooth_vectors`) tears a polygon *coverage*
//! apart: each polygon is smoothed on its own, so a boundary shared with a
//! neighbor is smoothed twice — independently — and the two results diverge,
//! opening gaps and slivers. This tool smooths the coverage as a whole so every
//! shared boundary stays coincident.
//!
//! It is the smoothing twin of `simplify_shared_edges` and reuses the same
//! **arc–node topology**:
//!
//! 1. Build a planar arc–node topology from every polygon boundary. A *node* is
//!    a vertex where the topology branches (undirected-edge degree != 2 — a
//!    shared edge meeting the coverage boundary, or three polygons meeting). An
//!    *arc* is a maximal chain of edges between two nodes; a shared boundary is
//!    one arc referenced by both neighbors, and a node-free ring (an island) is
//!    a single closed arc.
//! 2. Smooth **each arc once**, pinning its node endpoints so junctions never
//!    move.
//! 3. Reassemble every polygon from its (now-smoothed) arcs.
//!
//! Because a shared arc is smoothed a single time and both neighbors reference
//! the same smoothed vertices, the shared boundary is byte-identical on both
//! sides — no gaps, no slivers.
//!
//! Two smoothing algorithms are offered, matching ArcGIS:
//!
//! * **PAEK** (Polynomial Approximation with Exponential Kernel, the default):
//!   the arc is resampled at a fine step and each interior vertex is replaced by
//!   a Gaussian-weighted average of the vertices within `tolerance` arc-length
//!   of it. `tolerance` is the smoothing length in CRS units — larger values
//!   give smoother, rounder curves. Endpoints are pinned.
//! * **Bezier** (`bezier`): Chaikin corner-cutting subdivision, endpoints
//!   pinned. Like ArcGIS's Bezier interpolation it ignores `tolerance` and fits
//!   a smooth curve through the arc's vertices.
//!
//! `smooth_boundary` controls whether arcs on the coverage's outer edge
//! (degree-1, a single polygon) are smoothed too. `snap_tolerance` quantizes
//! vertices onto a grid before building the topology, so a coverage whose shared
//! vertices are only *nearly* coincident (a common artifact of imported data)
//! still shares arcs. Any non-polygon feature passes through unchanged.

use std::collections::{BTreeMap, HashMap, HashSet};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Geometry, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Chaikin doubles the vertex count per pass; stop before exceeding this.
const MAX_VERTICES: usize = 16_384;
/// PAEK resamples an arc at `tolerance / RESAMPLE_DIVISOR` before averaging, so
/// the Gaussian window spans enough samples to smooth cleanly.
const RESAMPLE_DIVISOR: f64 = 4.0;
/// Bezier (Chaikin) smoothing passes.
const CHAIKIN_ITERATIONS: usize = 3;

pub struct SmoothSharedEdgesTool;

impl Tool for SmoothSharedEdgesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "smooth_shared_edges",
            display_name: "Smooth Shared Edges",
            summary: "Smooth a polygon coverage while keeping boundaries shared between adjacent polygons coincident (no gaps or slivers), like ArcGIS Smooth Shared Edges.",
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
                    name: "algorithm",
                    description: "Smoothing algorithm: 'paek' (Gaussian kernel, tolerance-controlled; default) or 'bezier' (Chaikin corner-cutting, ignores tolerance).",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "PAEK smoothing length in CRS units: the arc-length window over which vertices are averaged. Larger values give smoother curves. Default 1.0. Ignored by the 'bezier' algorithm.",
                    required: false,
                },
                ToolParamSpec {
                    name: "smooth_boundary",
                    description: "Also smooth arcs on the coverage's outer boundary (edges belonging to a single polygon). Default true; set false to keep the outer boundary unchanged.",
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
        // vertex chain, remembering where it came from so the smoothed rings can
        // be written back into the right feature. Non-polygon features are left
        // untouched.
        let mut rings: Vec<Vec<P>> = Vec::new();
        let mut origins: Vec<RingOrigin> = Vec::new();
        for (fidx, feature) in layer.features.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            collect_rings(geom, fidx, prm.snap_tolerance, &mut rings, &mut origins);
        }

        ctx.progress.info(&format!(
            "{} feature(s): {} polygon ring(s) to smooth",
            layer.len(),
            rings.len()
        ));

        // Build the arc-node topology and smooth each arc once.
        let topo = Topology::build(&rings);
        let arc_smoothed: Vec<Vec<P>> = topo
            .arcs
            .iter()
            .map(|arc| {
                if prm.smooth_boundary || arc.shared {
                    smooth_arc(&arc.pts, arc.closed, prm.algorithm, prm.tolerance)
                } else {
                    arc.pts.clone()
                }
            })
            .collect();

        // Reassemble each ring from its (shared) smoothed arcs.
        let new_rings: Vec<Vec<P>> = (0..rings.len())
            .map(|ri| topo.reassemble_ring(ri, &arc_smoothed))
            .collect();

        // Write the smoothed rings back, grouped by feature. A ring that
        // collapsed below three vertices is dropped back to its original shape.
        let mut by_feature: HashMap<usize, Vec<(RingSlot, Vec<P>)>> = HashMap::new();
        let mut collapsed = 0usize;
        for (ri, origin) in origins.iter().enumerate() {
            let smoothed = &new_rings[ri];
            let ring = if smoothed.len() >= 3 {
                smoothed.clone()
            } else {
                collapsed += 1;
                rings[ri].clone()
            };
            by_feature
                .entry(origin.feature)
                .or_default()
                .push((origin.slot, ring));
        }

        let mut smoothed_features = 0usize;
        for (fidx, parts) in by_feature {
            let Some(geom) = layer.features[fidx].geometry.as_ref() else {
                continue;
            };
            if let Some(new_geom) = rebuild_geometry(geom, &parts) {
                layer.features[fidx].geometry = Some(new_geom);
                smoothed_features += 1;
            }
        }
        layer.extent = None; // geometries changed; drop the cached bbox

        let arcs_total = topo.arcs.len();
        let shared_arcs = topo.arcs.iter().filter(|a| a.shared).count();
        ctx.progress.info(&format!(
            "{arcs_total} arc(s) ({shared_arcs} shared); smoothed {smoothed_features} feature(s), {collapsed} ring(s) kept original"
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

#[derive(Clone, Copy, PartialEq, Eq)]
enum Algorithm {
    Paek,
    Bezier,
}

struct Params {
    algorithm: Algorithm,
    tolerance: f64,
    smooth_boundary: bool,
    snap_tolerance: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let algorithm = match parse_optional_str(args, "algorithm")? {
        None => Algorithm::Paek,
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "paek" => Algorithm::Paek,
            "bezier" => Algorithm::Bezier,
            other => {
                return Err(ToolError::Validation(format!(
                    "parameter 'algorithm' must be 'paek' or 'bezier', got '{other}'"
                )))
            }
        },
    };
    let tolerance = parse_optional_f64(args, "tolerance")?.unwrap_or(1.0);
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'tolerance' must be a positive number".to_string(),
        ));
    }
    let smooth_boundary = parse_optional_bool(args, "smooth_boundary")?.unwrap_or(true);
    let snap_tolerance = parse_optional_f64(args, "snap_tolerance")?.unwrap_or(0.0);
    if !(snap_tolerance >= 0.0 && snap_tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'snap_tolerance' must be a non-negative number".to_string(),
        ));
    }
    Ok(Params {
        algorithm,
        tolerance,
        smooth_boundary,
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
/// smoothed result is written back into the correct place.
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

    /// Rebuilds ring `ri` by concatenating its arcs' smoothed vertex chains,
    /// following each arc in the ring's original direction and dropping the
    /// shared node between consecutive arcs.
    fn reassemble_ring(&self, ri: usize, arc_smoothed: &[Vec<P>]) -> Vec<P> {
        let mut out: Vec<P> = Vec::new();
        for &(aid, forward) in &self.ring_arcs[ri].refs {
            let pts = &arc_smoothed[aid];
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
    let start_k = key(ring[0]);
    if let Some(pos) = arc.pts.iter().position(|p| key(*p) == start_k) {
        let next_in_arc = key(arc.pts[(pos + 1) % arc.pts.len()]);
        return key(ring[1 % ring.len()]) == next_in_arc;
    }
    true
}

// ── Smoothing ────────────────────────────────────────────────────────────────

/// Smooths one arc with the requested algorithm. Node endpoints are always
/// preserved (open arcs keep both; closed loops keep a single pinned start) so
/// both referencing rings stay identical.
fn smooth_arc(pts: &[P], closed: bool, algo: Algorithm, tol: f64) -> Vec<P> {
    match algo {
        Algorithm::Paek => {
            if closed {
                paek_closed(pts, tol)
            } else {
                paek_open(pts, tol)
            }
        }
        Algorithm::Bezier => {
            if closed {
                chaikin(pts, CHAIKIN_ITERATIONS, true)
            } else {
                chaikin(pts, CHAIKIN_ITERATIONS, false)
            }
        }
    }
}

/// PAEK on an open polyline: resample at a fine step, then replace each interior
/// vertex with a Gaussian-weighted average of the samples within `tol`
/// arc-length. Endpoints are pinned.
fn paek_open(pts: &[P], tol: f64) -> Vec<P> {
    if pts.len() < 3 || tol <= 0.0 {
        return pts.to_vec();
    }
    let step = resample_step(pts, tol, false);
    let samples = densify(pts, step, false);
    gaussian_smooth(&samples, tol, false)
}

/// PAEK on a closed loop: rotate to pin the start vertex (so both rings that
/// reference the loop smooth identically), smooth as an open chain start→…→start,
/// then drop the closing duplicate.
fn paek_closed(pts: &[P], tol: f64) -> Vec<P> {
    let n = pts.len();
    if n < 4 || tol <= 0.0 {
        return pts.to_vec();
    }
    let mut closed: Vec<P> = pts.to_vec();
    closed.push(pts[0]); // reclose so the smoother pins the start at both ends
    let step = resample_step(pts, tol, true);
    let samples = densify(&closed, step, false);
    let mut out = gaussian_smooth(&samples, tol, false);
    while out.len() >= 2 && dist(out[0], *out.last().unwrap()) <= 1e-9 {
        out.pop();
    }
    if out.len() < 3 {
        pts.to_vec()
    } else {
        out
    }
}

/// Chooses a resample step so the Gaussian window spans several samples while
/// keeping the total vertex count bounded.
fn resample_step(pts: &[P], tol: f64, closed: bool) -> f64 {
    let len = perimeter(pts, closed);
    let fine = tol / RESAMPLE_DIVISOR;
    let by_budget = len / (MAX_VERTICES as f64);
    fine.max(by_budget).max(f64::MIN_POSITIVE)
}

/// Gaussian-weighted moving average along arc length. `tol` is the window
/// half-... actually the full smoothing length: sigma = tol/2, samples beyond
/// `tol` of a vertex are ignored. Endpoints stay fixed (open) — for the closed
/// case the caller pins the start by reclosing the chain.
fn gaussian_smooth(pts: &[P], tol: f64, _closed: bool) -> Vec<P> {
    let n = pts.len();
    if n < 3 {
        return pts.to_vec();
    }
    // Cumulative arc-length position of each sample.
    let mut s = vec![0.0f64; n];
    for i in 1..n {
        s[i] = s[i - 1] + dist(pts[i - 1], pts[i]);
    }
    let sigma = (tol * 0.5).max(f64::MIN_POSITIVE);
    let window = tol;
    let mut out = Vec::with_capacity(n);
    out.push(pts[0]); // pinned endpoint
    for i in 1..n - 1 {
        let (mut sx, mut sy, mut sw) = (0.0, 0.0, 0.0);
        // Walk outward from i until beyond the window on both sides.
        let mut j = i;
        while j < n && s[j] - s[i] <= window {
            let w = gauss((s[j] - s[i]) / sigma);
            sx += pts[j].x * w;
            sy += pts[j].y * w;
            sw += w;
            j += 1;
        }
        let mut k = i;
        while k > 0 {
            k -= 1;
            if s[i] - s[k] > window {
                break;
            }
            let w = gauss((s[i] - s[k]) / sigma);
            sx += pts[k].x * w;
            sy += pts[k].y * w;
            sw += w;
        }
        if sw > 0.0 {
            out.push(P {
                x: sx / sw,
                y: sy / sw,
            });
        } else {
            out.push(pts[i]);
        }
    }
    out.push(pts[n - 1]); // pinned endpoint
    dedup_p(&out)
}

fn gauss(z: f64) -> f64 {
    (-0.5 * z * z).exp()
}

/// Chaikin corner cutting: each segment is replaced by its 1/4 and 3/4 points.
/// Open arcs keep their endpoints fixed; closed arcs are cut cyclically but the
/// caller relies on the pinned start staying put — Chaikin keeps a vertex only
/// approximately, so for closed loops the first vertex is re-pinned afterward.
fn chaikin(pts: &[P], iterations: usize, closed: bool) -> Vec<P> {
    if pts.len() < 3 {
        return pts.to_vec();
    }
    let pinned_start = pts[0];
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
    if closed {
        // Re-pin the node so both referencing rings stay coincident, and keep
        // the chain unclosed.
        cur.insert(0, pinned_start);
        while cur.len() >= 2 && dist(cur[0], *cur.last().unwrap()) <= 1e-9 {
            cur.pop();
        }
    }
    dedup_p(&cur)
}

// ── Geometry primitives ───────────────────────────────────────────────────────

fn lerp(a: P, b: P, t: f64) -> P {
    P {
        x: a.x + (b.x - a.x) * t,
        y: a.y + (b.y - a.y) * t,
    }
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
        let pieces = (dist(a, b) / max_len).ceil().max(1.0) as usize;
        for j in 1..pieces {
            out.push(lerp(a, b, j as f64 / pieces as f64));
        }
    }
    if !closed {
        out.push(pts[n - 1]);
    }
    out
}

/// Drops consecutive near-duplicate vertices.
fn dedup_p(pts: &[P]) -> Vec<P> {
    let mut out: Vec<P> = Vec::with_capacity(pts.len());
    for &p in pts {
        if out.last().is_none_or(|last| dist(*last, p) > 1e-12) {
            out.push(p);
        }
    }
    out
}

// ── Geometry rebuild ─────────────────────────────────────────────────────────

fn pts_to_ring(pts: &[P]) -> Ring {
    Ring::new(pts.iter().map(|p| Coord::xy(p.x, p.y)).collect())
}

/// Rewrites `geom`'s rings from the per-slot smoothed chains in `parts`.
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
        let out = SmoothSharedEdgesTool.run(&args, &ctx()).unwrap();
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

    /// The shared-boundary vertices of a polygon: those near x=10 (the jagged
    /// interior edge), as a sorted bit-exact set. Only reliable when the outer
    /// boundary is left unsmoothed (`smooth_boundary=false`) so it contributes
    /// no vertices in the band.
    fn shared_band(v: &[(f64, f64)]) -> Vec<(u64, u64)> {
        let mut s: Vec<(u64, u64)> = v
            .iter()
            .filter(|(x, _)| (*x - 10.0).abs() < 3.0)
            .map(|(x, y)| (x.to_bits(), y.to_bits()))
            .collect();
        s.sort_unstable();
        s
    }

    /// A boundary shared by two polygons must be smoothed identically in both —
    /// the two copies of the shared arc must be point-for-point equal.
    #[test]
    fn shared_edge_stays_coincident() {
        // A jagged vertical boundary at x≈10 shared by a left and right polygon.
        let shared: Vec<(f64, f64)> = vec![
            (10.0, 0.0),
            (11.0, 2.5),
            (9.0, 5.0),
            (11.0, 7.5),
            (10.0, 10.0),
        ];
        // Left polygon: 0..10 with the shared edge on its right (top-to-bottom).
        let mut left = vec![(0.0, 0.0), (0.0, 10.0)];
        left.extend(shared.iter().rev().copied());
        // Right polygon: its left edge follows the same shared chain.
        let right_ring: Vec<(f64, f64)> = {
            let mut r = vec![(10.0, 0.0), (20.0, 0.0), (20.0, 10.0), (10.0, 10.0)];
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

        // Leave the outer boundary untouched so the x≈10 band is exactly the
        // shared arc, then assert both polygons' copies of it are identical.
        let (out, layer) =
            run_tool(json!({ "input": input, "tolerance": 3.0, "smooth_boundary": false }));
        assert!(out.outputs["shared_arc_count"].as_u64().unwrap() >= 1);

        let l = exterior(&layer, 0);
        let r = exterior(&layer, 1);
        assert_eq!(
            shared_band(&l),
            shared_band(&r),
            "shared boundary diverged: {l:?} vs {r:?}"
        );
        // And smoothing actually happened: more vertices than the 5-point input.
        assert!(
            l.iter().filter(|(x, _)| (*x - 10.0).abs() < 3.0).count() > 5,
            "shared edge was not smoothed: {l:?}"
        );
    }

    /// A vertex where three polygons meet is a node and must be preserved.
    #[test]
    fn triple_junction_is_preserved() {
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
        for r in [a, b, c] {
            layer
                .add_feature(Some(Geometry::polygon(r, vec![])), &[])
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (_out, layer) = run_tool(json!({ "input": input, "tolerance": 4.0 }));
        // (10,10) is a triple junction; it must survive in every polygon.
        for idx in 0..3 {
            let ext = exterior(&layer, idx);
            assert!(
                ext.iter()
                    .any(|(x, y)| (*x - 10.0).abs() < 1e-9 && (*y - 10.0).abs() < 1e-9),
                "polygon {idx} lost the triple junction (10,10): {ext:?}"
            );
        }
    }

    /// PAEK smoothing pulls in a sharp protruding spike on the shared edge.
    #[test]
    fn paek_rounds_sharp_corners() {
        // Two polygons sharing a boundary with a hard spike out to x=14.
        let shared: Vec<(f64, f64)> = vec![
            (10.0, 0.0),
            (10.0, 4.0),
            (14.0, 5.0),
            (10.0, 6.0),
            (10.0, 10.0),
        ];
        let mut left = vec![(0.0, 0.0), (0.0, 10.0)];
        left.extend(shared.iter().rev().copied());
        let mut right = vec![(10.0, 0.0), (20.0, 0.0), (20.0, 10.0), (10.0, 10.0)];
        for &pt in shared.iter().rev().skip(1).take(3) {
            right.push(pt);
        }
        let mut layer = Layer::new("cov");
        layer
            .add_feature(Some(Geometry::polygon(ring(&left), vec![])), &[])
            .unwrap();
        layer
            .add_feature(Some(Geometry::polygon(ring(&right), vec![])), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        // A generous window pulls the x=14 spike well back toward the edge.
        let (_o, layer) =
            run_tool(json!({ "input": input, "tolerance": 8.0, "smooth_boundary": false }));
        let max_x = exterior(&layer, 0)
            .iter()
            .filter(|(x, _)| (*x - 10.0).abs() < 6.0)
            .map(|(x, _)| *x)
            .fold(0.0, f64::max);
        assert!(
            max_x < 13.0,
            "sharp spike not pulled in, max x still {max_x}"
        );
    }

    /// The 'bezier' algorithm also smooths and stays coverage-safe.
    #[test]
    fn bezier_algorithm_smooths_and_stays_shared() {
        let shared: Vec<(f64, f64)> = vec![(10.0, 0.0), (12.0, 3.0), (8.0, 6.0), (10.0, 10.0)];
        let mut left = vec![(0.0, 0.0), (0.0, 10.0)];
        left.extend(shared.iter().rev().copied());
        let mut right = vec![(10.0, 0.0), (20.0, 0.0), (20.0, 10.0), (10.0, 10.0)];
        for &pt in shared.iter().rev().skip(1).take(2) {
            right.push(pt);
        }
        let mut layer = Layer::new("cov");
        layer
            .add_feature(Some(Geometry::polygon(ring(&left), vec![])), &[])
            .unwrap();
        layer
            .add_feature(Some(Geometry::polygon(ring(&right), vec![])), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) =
            run_tool(json!({ "input": input, "algorithm": "bezier", "smooth_boundary": false }));
        assert!(out.outputs["shared_arc_count"].as_u64().unwrap() >= 1);
        let l = exterior(&layer, 0);
        let r = exterior(&layer, 1);
        assert_eq!(
            shared_band(&l),
            shared_band(&r),
            "bezier shared boundary diverged"
        );
        assert!(
            l.iter().filter(|(x, _)| (*x - 10.0).abs() < 3.0).count() > 4,
            "bezier did not add smoothing vertices to the shared edge"
        );
    }

    /// smooth_boundary=false keeps the coverage's outer boundary vertices,
    /// while the shared interior edge is still smoothed.
    #[test]
    fn boundary_preserved_when_requested() {
        let shared: Vec<(f64, f64)> = vec![(10.0, 0.0), (12.0, 5.0), (10.0, 10.0)];
        let mut left = vec![(0.0, 0.0), (0.0, 10.0)];
        left.extend(shared.iter().rev().copied());
        let mut right = vec![(10.0, 0.0), (20.0, 0.0), (20.0, 10.0), (10.0, 10.0)];
        right.push(shared[1]);
        let build = || {
            let mut layer = Layer::new("cov");
            layer
                .add_feature(Some(Geometry::polygon(ring(&left), vec![])), &[])
                .unwrap();
            layer
                .add_feature(Some(Geometry::polygon(ring(&right), vec![])), &[])
                .unwrap();
            let id = memory_store::put_vector(layer);
            memory_store::make_vector_memory_path(&id)
        };

        let (_o, kept) =
            run_tool(json!({ "input": build(), "tolerance": 3.0, "smooth_boundary": false }));
        // The outer corner (20,0) must be untouched.
        assert!(
            exterior(&kept, 1)
                .iter()
                .any(|(x, y)| (*x - 20.0).abs() < 1e-9 && *y < 1e-9),
            "outer boundary corner was moved despite smooth_boundary=false"
        );
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
        let tool = SmoothSharedEdgesTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "x.geojson", "tolerance": 0 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "tolerance": -1 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "snap_tolerance": -1 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "algorithm": "spline" })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "tolerance": 1.0 })).is_ok());
        assert!(bad(json!({ "input": "x.geojson", "algorithm": "bezier" })).is_ok());
        assert!(bad(json!({ "input": "x.geojson", "smooth_boundary": "false" })).is_ok());
    }
}
