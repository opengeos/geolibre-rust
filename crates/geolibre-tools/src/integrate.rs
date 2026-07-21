//! GeoLibre tool: cluster-tolerance vertex snapping across features (Integrate).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Integrate* (Data Management), with the
//! editing *Snap* tool's vertex and edge modes. It snaps **every vertex across
//! every feature** that falls within a cluster tolerance to a single shared
//! location, so nearly-coincident shared boundaries become exactly coincident —
//! the topology-cleaning step that the repo's coverage-safe family
//! (`simplify_shared_edges`, `smooth_shared_edges`, `polygon_neighbors`,
//! `eliminate_polygons`) all assume has already happened.
//!
//! The bundled suite only has `snap_endnodes` (polyline endpoints); nothing
//! integrates all vertices, polygon borders included.
//!
//! Two passes:
//!
//! 1. **Vertex clustering.** All distinct vertices are grid-hashed at the
//!    tolerance scale and union-found into clusters (any two within `tolerance`
//!    join). Every vertex is moved to its cluster's centroid, so vertices that
//!    were merely *near*-coincident become bit-identical.
//! 2. **Edge snapping** (`snap_to_edges`, default on). A clustered vertex that
//!    lands within `tolerance` of another feature's *segment* — a T-junction —
//!    is inserted into that segment as a shared vertex, so the junction is
//!    topologically closed.
//!
//! Works on any geometry type; distances are in the layer CRS units. Reports the
//! number of vertices moved and inserted.

use std::collections::{BTreeMap, HashMap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Geometry, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct IntegrateTool;

impl Tool for IntegrateTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "integrate",
            display_name: "Integrate",
            summary: "Snap all vertices across all features within a cluster tolerance to a shared location so nearly-coincident shared boundaries become exactly coincident, optionally inserting T-junction vertices onto nearby edges — like ArcGIS Integrate.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (any geometry type).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Cluster distance in CRS units: vertices within this distance snap to a shared location. Required.",
                    required: true,
                },
                ToolParamSpec {
                    name: "snap_to_edges",
                    description: "Also insert a vertex where one feature's vertex lands within tolerance of another feature's segment (T-junctions). Default true.",
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

        // ── Pass 0: collect every distinct vertex ─────────────────────────────
        let mut distinct: HashMap<Key, usize> = HashMap::new();
        let mut coords: Vec<P> = Vec::new();
        for feature in layer.features.iter() {
            if let Some(g) = feature.geometry.as_ref() {
                for_each_vertex(g, &mut |x, y| {
                    let k = key(x, y);
                    distinct.entry(k).or_insert_with(|| {
                        coords.push(P { x, y });
                        coords.len() - 1
                    });
                });
            }
        }
        let n = coords.len();
        ctx.progress.info(&format!(
            "{n} distinct vertex(es); clustering at {}",
            prm.tolerance
        ));

        // ── Pass 1: cluster within tolerance (grid-hashed union-find) ─────────
        let tol = prm.tolerance;
        let cell_of =
            |p: P| -> (i64, i64) { ((p.x / tol).floor() as i64, (p.y / tol).floor() as i64) };
        let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
        for (i, &p) in coords.iter().enumerate() {
            grid.entry(cell_of(p)).or_default().push(i);
        }
        let mut uf = UnionFind::new(n);
        for (i, &p) in coords.iter().enumerate() {
            let (cx, cy) = cell_of(p);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    if let Some(bucket) = grid.get(&(cx + dx, cy + dy)) {
                        for &j in bucket {
                            if j > i && dist(p, coords[j]) <= tol {
                                uf.union(i, j);
                            }
                        }
                    }
                }
            }
        }
        // Cluster centroid: mean of the distinct member coordinates.
        let mut sum: HashMap<usize, (f64, f64, usize)> = HashMap::new();
        for (i, &p) in coords.iter().enumerate() {
            let r = uf.find(i);
            let e = sum.entry(r).or_insert((0.0, 0.0, 0));
            e.0 += p.x;
            e.1 += p.y;
            e.2 += 1;
        }
        let centroid_of = |i: usize, uf: &mut UnionFind| -> P {
            let (sx, sy, c) = sum[&uf.find(i)];
            P {
                x: sx / c as f64,
                y: sy / c as f64,
            }
        };
        // Map every distinct coord key to its cluster centroid.
        let mut snap: HashMap<Key, P> = HashMap::with_capacity(n);
        for (&k, &idx) in &distinct {
            snap.insert(k, centroid_of(idx, &mut uf));
        }
        let snap_fn = |x: f64, y: f64| -> P { snap.get(&key(x, y)).copied().unwrap_or(P { x, y }) };

        // Unique cluster centroids, for the edge-snap phase.
        let mut centroids: Vec<P> = Vec::new();
        {
            let mut seen: HashMap<Key, ()> = HashMap::new();
            for p in snap.values() {
                if seen.insert(key(p.x, p.y), ()).is_none() {
                    centroids.push(*p);
                }
            }
        }

        // ── Rebuild geometries with snapped vertices; count moves ─────────────
        let mut moved = 0usize;
        for feature in layer.features.iter_mut() {
            if let Some(g) = feature.geometry.as_ref() {
                let new = snap_geometry(g, &snap_fn, &mut moved);
                feature.geometry = Some(new);
            }
        }

        // ── Pass 2: edge snapping (insert T-junction vertices) ────────────────
        let mut inserted = 0usize;
        if prm.snap_to_edges && !centroids.is_empty() {
            let mut cgrid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
            for (i, &p) in centroids.iter().enumerate() {
                cgrid.entry(cell_of(p)).or_default().push(i);
            }
            for feature in layer.features.iter_mut() {
                if let Some(g) = feature.geometry.as_ref() {
                    let new = insert_on_edges(g, &centroids, &cgrid, tol, &mut inserted);
                    feature.geometry = Some(new);
                }
            }
        }

        layer.extent = None;
        ctx.progress.info(&format!(
            "{} cluster(s); {moved} vertex(es) moved, {inserted} inserted on edges",
            sum.len()
        ));

        let feature_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("cluster_count".to_string(), json!(sum.len()));
        outputs.insert("vertices_moved".to_string(), json!(moved));
        outputs.insert("vertices_inserted".to_string(), json!(inserted));
        Ok(ToolRunResult { outputs })
    }
}

// ── Geometry traversal / rebuild ─────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
struct P {
    x: f64,
    y: f64,
}

type Key = (u64, u64);

fn key(x: f64, y: f64) -> Key {
    (x.to_bits(), y.to_bits())
}

fn dist(a: P, b: P) -> f64 {
    (a.x - b.x).hypot(a.y - b.y)
}

fn for_each_vertex(geom: &Geometry, f: &mut impl FnMut(f64, f64)) {
    match geom {
        Geometry::Point(c) => f(c.x, c.y),
        Geometry::LineString(cs) | Geometry::MultiPoint(cs) => cs.iter().for_each(|c| f(c.x, c.y)),
        Geometry::MultiLineString(lines) => lines.iter().flatten().for_each(|c| f(c.x, c.y)),
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            exterior.coords().iter().for_each(|c| f(c.x, c.y));
            interiors
                .iter()
                .for_each(|r| r.coords().iter().for_each(|c| f(c.x, c.y)));
        }
        Geometry::MultiPolygon(parts) => {
            for (e, holes) in parts {
                e.coords().iter().for_each(|c| f(c.x, c.y));
                holes
                    .iter()
                    .for_each(|r| r.coords().iter().for_each(|c| f(c.x, c.y)));
            }
        }
        Geometry::GeometryCollection(gs) => gs.iter().for_each(|g| for_each_vertex(g, f)),
    }
}

/// Snaps every vertex of a geometry to its cluster centroid, dropping vertices
/// that collapse onto their neighbour. Counts vertices whose position changed.
fn snap_geometry(geom: &Geometry, snap: &impl Fn(f64, f64) -> P, moved: &mut usize) -> Geometry {
    let snap_chain = |cs: &[Coord], closed: bool, moved: &mut usize| -> Vec<Coord> {
        let mut out: Vec<Coord> = Vec::with_capacity(cs.len());
        for c in cs {
            let s = snap(c.x, c.y);
            if s.x != c.x || s.y != c.y {
                *moved += 1;
            }
            if out.last().is_none_or(|l| key(l.x, l.y) != key(s.x, s.y)) {
                out.push(Coord::xy(s.x, s.y));
            }
        }
        if closed {
            // Keep the ring closed if the writer expects a duplicate; Ring stores
            // it unclosed, so drop a trailing duplicate of the first vertex.
            while out.len() >= 2
                && key(out[0].x, out[0].y) == key(out.last().unwrap().x, out.last().unwrap().y)
            {
                out.pop();
            }
        }
        out
    };
    match geom {
        Geometry::Point(c) => {
            let s = snap(c.x, c.y);
            if s.x != c.x || s.y != c.y {
                *moved += 1;
            }
            Geometry::Point(Coord::xy(s.x, s.y))
        }
        Geometry::MultiPoint(cs) => Geometry::MultiPoint(
            cs.iter()
                .map(|c| {
                    let s = snap(c.x, c.y);
                    if s.x != c.x || s.y != c.y {
                        *moved += 1;
                    }
                    Coord::xy(s.x, s.y)
                })
                .collect(),
        ),
        Geometry::LineString(cs) => Geometry::LineString(snap_chain(cs, false, moved)),
        Geometry::MultiLineString(lines) => {
            Geometry::MultiLineString(lines.iter().map(|l| snap_chain(l, false, moved)).collect())
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => Geometry::Polygon {
            exterior: Ring::new(snap_chain(exterior.coords(), true, moved)),
            interiors: interiors
                .iter()
                .map(|r| Ring::new(snap_chain(r.coords(), true, moved)))
                .collect(),
        },
        Geometry::MultiPolygon(parts) => Geometry::MultiPolygon(
            parts
                .iter()
                .map(|(e, holes)| {
                    (
                        Ring::new(snap_chain(e.coords(), true, moved)),
                        holes
                            .iter()
                            .map(|r| Ring::new(snap_chain(r.coords(), true, moved)))
                            .collect(),
                    )
                })
                .collect(),
        ),
        Geometry::GeometryCollection(gs) => {
            Geometry::GeometryCollection(gs.iter().map(|g| snap_geometry(g, snap, moved)).collect())
        }
    }
}

/// Inserts cluster centroids that lie within `tol` of a line/ring segment (a
/// T-junction) into that segment, ordered along it.
fn insert_on_edges(
    geom: &Geometry,
    centroids: &[P],
    cgrid: &HashMap<(i64, i64), Vec<usize>>,
    tol: f64,
    inserted: &mut usize,
) -> Geometry {
    let densify_chain = |cs: &[Coord], inserted: &mut usize| -> Vec<Coord> {
        if cs.len() < 2 {
            return cs.to_vec();
        }
        let mut out: Vec<Coord> = Vec::with_capacity(cs.len());
        for w in cs.windows(2) {
            let a = P {
                x: w[0].x,
                y: w[0].y,
            };
            let b = P {
                x: w[1].x,
                y: w[1].y,
            };
            out.push(Coord::xy(a.x, a.y));
            // Candidate centroids near this segment (grid cells spanned by its bbox).
            let mut inserts: Vec<(f64, P)> = Vec::new();
            let (minx, maxx) = (a.x.min(b.x), a.x.max(b.x));
            let (miny, maxy) = (a.y.min(b.y), a.y.max(b.y));
            let (c0x, c1x) = (
                ((minx - tol) / tol).floor() as i64,
                ((maxx + tol) / tol).floor() as i64,
            );
            let (c0y, c1y) = (
                ((miny - tol) / tol).floor() as i64,
                ((maxy + tol) / tol).floor() as i64,
            );
            for gx in c0x..=c1x {
                for gy in c0y..=c1y {
                    if let Some(bucket) = cgrid.get(&(gx, gy)) {
                        for &ci in bucket {
                            let c = centroids[ci];
                            // Skip the segment's own endpoints.
                            if key(c.x, c.y) == key(a.x, a.y) || key(c.x, c.y) == key(b.x, b.y) {
                                continue;
                            }
                            if let Some(t) = interior_projection(c, a, b, tol) {
                                inserts.push((t, c));
                            }
                        }
                    }
                }
            }
            inserts.sort_by(|p, q| p.0.total_cmp(&q.0));
            for (_, c) in inserts {
                if out.last().is_none_or(|l| key(l.x, l.y) != key(c.x, c.y)) {
                    out.push(Coord::xy(c.x, c.y));
                    *inserted += 1;
                }
            }
        }
        out.push(cs.last().unwrap().clone());
        out
    };
    match geom {
        Geometry::LineString(cs) => Geometry::LineString(densify_chain(cs, inserted)),
        Geometry::MultiLineString(lines) => {
            Geometry::MultiLineString(lines.iter().map(|l| densify_chain(l, inserted)).collect())
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => Geometry::Polygon {
            exterior: ring_insert(exterior, &densify_chain, inserted),
            interiors: interiors
                .iter()
                .map(|r| ring_insert(r, &densify_chain, inserted))
                .collect(),
        },
        Geometry::MultiPolygon(parts) => Geometry::MultiPolygon(
            parts
                .iter()
                .map(|(e, holes)| {
                    (
                        ring_insert(e, &densify_chain, inserted),
                        holes
                            .iter()
                            .map(|r| ring_insert(r, &densify_chain, inserted))
                            .collect(),
                    )
                })
                .collect(),
        ),
        other => other.clone(),
    }
}

/// Applies segment insertion to a closed ring (temporarily closed for the walk).
fn ring_insert(
    ring: &Ring,
    densify: &impl Fn(&[Coord], &mut usize) -> Vec<Coord>,
    inserted: &mut usize,
) -> Ring {
    let mut closed = ring.coords().to_vec();
    if let (Some(first), Some(last)) = (closed.first().cloned(), closed.last().cloned()) {
        if key(first.x, first.y) != key(last.x, last.y) {
            closed.push(first);
        }
    }
    let mut out = densify(&closed, inserted);
    // Drop the trailing closing duplicate the ring stores unclosed.
    while out.len() >= 2
        && key(out[0].x, out[0].y) == key(out.last().unwrap().x, out.last().unwrap().y)
    {
        out.pop();
    }
    Ring::new(out)
}

/// If `p` projects onto the interior of segment `a`-`b` within `tol`, returns the
/// projection parameter t ∈ (0, 1).
fn interior_projection(p: P, a: P, b: P, tol: f64) -> Option<f64> {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let len2 = dx * dx + dy * dy;
    if len2 <= 0.0 {
        return None;
    }
    let t = ((p.x - a.x) * dx + (p.y - a.y) * dy) / len2;
    if t <= 1e-9 || t >= 1.0 - 1e-9 {
        return None; // at or beyond an endpoint
    }
    let proj = P {
        x: a.x + t * dx,
        y: a.y + t * dy,
    };
    (dist(p, proj) <= tol).then_some(t)
}

// ── Union-find ───────────────────────────────────────────────────────────────

struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}
impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }
    fn find(&mut self, x: usize) -> usize {
        let mut r = x;
        while self.parent[r] != r {
            r = self.parent[r];
        }
        let mut c = x;
        while self.parent[c] != r {
            let next = self.parent[c];
            self.parent[c] = r;
            c = next;
        }
        r
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        if self.rank[ra] < self.rank[rb] {
            self.parent[ra] = rb;
        } else if self.rank[ra] > self.rank[rb] {
            self.parent[rb] = ra;
        } else {
            self.parent[rb] = ra;
            self.rank[ra] += 1;
        }
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    tolerance: f64,
    snap_to_edges: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let tolerance = parse_optional_f64(args, "tolerance")?.ok_or_else(|| {
        ToolError::Validation("required parameter 'tolerance' is missing".to_string())
    })?;
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "'tolerance' must be a positive number".to_string(),
        ));
    }
    let snap_to_edges = parse_optional_bool(args, "snap_to_edges")?.unwrap_or(true);
    Ok(Params {
        tolerance,
        snap_to_edges,
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

    fn poly(coords: &[(f64, f64)]) -> Geometry {
        Geometry::polygon(
            coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect(),
            vec![],
        )
    }

    fn layer_of(geoms: Vec<Geometry>) -> String {
        let mut l = Layer::new("v")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        for g in geoms {
            l.add_feature(Some(g), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = IntegrateTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn exterior(layer: &Layer, i: usize) -> Vec<(f64, f64)> {
        match layer.features[i].geometry.as_ref().unwrap() {
            Geometry::Polygon { exterior, .. } => {
                exterior.coords().iter().map(|c| (c.x, c.y)).collect()
            }
            other => panic!("expected polygon, got {other:?}"),
        }
    }

    /// Two polygons whose shared border vertices differ by < tolerance become
    /// bit-identical on that border after integrate.
    #[test]
    fn snaps_near_coincident_borders() {
        // Left shares x=10 with right, but right's border is offset by 0.3.
        let left = poly(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]);
        let right = poly(&[(10.3, 0.0), (20.0, 0.0), (20.0, 10.0), (10.3, 10.0)]);
        let (out, layer) = run(json!({ "input": layer_of(vec![left, right]), "tolerance": 1.0 }));
        assert!(out.outputs["vertices_moved"].as_u64().unwrap() >= 2);
        // The shared border vertices are now identical between the two polygons.
        let l = exterior(&layer, 0);
        let r = exterior(&layer, 1);
        let on_border = |v: &[(f64, f64)]| -> Vec<(u64, u64)> {
            let mut s: Vec<(u64, u64)> = v
                .iter()
                .filter(|(x, _)| (*x - 10.15).abs() < 1.0)
                .map(|(x, y)| (x.to_bits(), y.to_bits()))
                .collect();
            s.sort_unstable();
            s
        };
        assert_eq!(
            on_border(&l),
            on_border(&r),
            "borders not coincident: {l:?} vs {r:?}"
        );
    }

    /// A vertex landing on another feature's edge is inserted (T-junction).
    #[test]
    fn inserts_t_junction_vertex() {
        // A long horizontal edge, and a second polygon whose corner touches its
        // middle (at x=50) but the long edge has no vertex there.
        let bar = poly(&[(0.0, 0.0), (100.0, 0.0), (100.0, 5.0), (0.0, 5.0)]);
        let tee = poly(&[(50.0, 5.0), (60.0, 5.0), (60.0, 30.0), (50.0, 30.0)]);
        let (out, layer) = run(json!({
            "input": layer_of(vec![bar, tee]), "tolerance": 1.0, "snap_to_edges": true,
        }));
        assert!(
            out.outputs["vertices_inserted"].as_u64().unwrap() >= 1,
            "expected a T-junction insert"
        );
        // The bar's top edge (y=5) now has vertices at x=50 and/or x=60.
        let bar_v = exterior(&layer, 0);
        assert!(
            bar_v.iter().any(|(x, y)| (*y - 5.0).abs() < 1e-9
                && ((*x - 50.0).abs() < 1e-9 || (*x - 60.0).abs() < 1e-9)),
            "T-junction vertex not inserted into the bar edge: {bar_v:?}"
        );
    }

    /// snap_to_edges=false leaves segments un-split.
    #[test]
    fn edge_snap_can_be_disabled() {
        let bar = poly(&[(0.0, 0.0), (100.0, 0.0), (100.0, 5.0), (0.0, 5.0)]);
        let tee = poly(&[(50.0, 5.0), (60.0, 5.0), (60.0, 30.0), (50.0, 30.0)]);
        let (out, _l) = run(json!({
            "input": layer_of(vec![bar, tee]), "tolerance": 1.0, "snap_to_edges": false,
        }));
        assert_eq!(out.outputs["vertices_inserted"], json!(0));
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            IntegrateTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "x.geojson" })).is_err()); // no tolerance
        assert!(bad(json!({ "input": "x.geojson", "tolerance": 0 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "tolerance": 1.0 })).is_ok());
    }
}
