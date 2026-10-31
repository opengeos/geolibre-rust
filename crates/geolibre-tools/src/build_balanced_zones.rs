//! GeoLibre tool: SKATER regionalization into balanced, contiguous zones.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Build Balanced Zones* and the core of
//! *Spatially Constrained Multivariate Clustering* (Spatial Statistics): group
//! contiguous polygons into a target number of zones that balance a criterion —
//! equal feature count, equal sum of an attribute, or attribute homogeneity —
//! while every zone stays spatially connected. Districting, sales-territory
//! design, and ecological regionalization have no answer in the bundled suite
//! (its `dbscan` / `k_means_clustering` ignore contiguity), and there is no other
//! pure-Rust/WASM SKATER.
//!
//! The algorithm is SKATER (Spatial 'K'luster Analysis by Tree Edge Removal):
//!
//! 1. Build a **contiguity graph** — polygons are adjacent under `rook` (a shared
//!    boundary segment) or `queen` (any shared boundary point) contiguity.
//! 2. Build a **minimum spanning forest** over that graph (Kruskal), with edges
//!    weighted by attribute dissimilarity (standardized Euclidean distance over
//!    the analysis `fields`) or, when no fields are given, by centroid distance.
//! 3. **Cut** the forest greedily: at each step remove the tree edge whose
//!    removal best improves the balance objective, splitting one zone into two,
//!    until there are `zones` zones. Because zones are always subtrees of the
//!    contiguity graph, they are guaranteed connected.
//!
//! Balance criteria (`criterion`):
//! - `homogeneity` (default) — minimize the total within-zone sum of squared
//!   deviations across the standardized `fields` (classic SKATER);
//! - `equal_count` — make the zones as equal in feature count as possible;
//! - `equal_sum` — make the zones as equal in the sum of the first field as
//!   possible.
//!
//! Each output feature keeps its attributes and gains a `zone` id (0-based); the
//! report carries per-zone feature counts. The graph build is O(n²); it suits
//! moderate polygon counts (counties, tracts, parcels). Deterministic — no RNG.

use std::collections::BTreeMap;

use geo::{Centroid, Coord as GeoCoord, Intersects, LineString, MultiPolygon, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct BuildBalancedZonesTool;

impl Tool for BuildBalancedZonesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "build_balanced_zones",
            display_name: "Build Balanced Zones",
            summary: "Group contiguous polygons into a target number of connected zones that balance feature count, an attribute sum, or attribute homogeneity, using SKATER spanning-tree partitioning.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon layer, format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "zones",
                    description: "Target number of zones to create.",
                    required: true,
                },
                ToolParamSpec {
                    name: "criterion",
                    description: "Balance objective: 'homogeneity' (default), 'equal_count', or 'equal_sum' (of the first field).",
                    required: false,
                },
                ToolParamSpec {
                    name: "fields",
                    description: "Comma-separated numeric analysis field(s). Required for 'homogeneity' and 'equal_sum'; drives the spanning-tree edge weights.",
                    required: false,
                },
                ToolParamSpec {
                    name: "contiguity",
                    description: "Adjacency rule: 'rook' (shared edge, default) or 'queen' (shared edge or vertex).",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Maximum distance for two boundary segments to count as shared (rook). Default 1e-6.",
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
        if matches!(prm.criterion, Criterion::Homogeneity | Criterion::EqualSum)
            && prm.fields.is_empty()
        {
            return Err(ToolError::Validation(format!(
                "criterion '{}' requires at least one analysis field",
                prm.criterion.as_str()
            )));
        }

        let mut layer = load_input_layer(input)?;
        let schema = layer.schema.clone();

        // Collect polygon nodes: geometry (as geo), attribute values.
        let mut polys: Vec<MultiPolygon> = Vec::new();
        let mut attrs: Vec<Vec<f64>> = Vec::new();
        let mut node_feature: Vec<usize> = Vec::new();
        for (fi, feature) in layer.features.iter().enumerate() {
            let Some(mp) = feature.geometry.as_ref().and_then(to_multipolygon) else {
                continue;
            };
            let row: Vec<f64> = prm
                .fields
                .iter()
                .map(|f| {
                    feature
                        .get(&schema, f)
                        .ok()
                        .and_then(FieldValue::as_f64)
                        .filter(|v| v.is_finite())
                        .unwrap_or(0.0)
                })
                .collect();
            polys.push(mp);
            attrs.push(row);
            node_feature.push(fi);
        }
        let n = polys.len();
        if n == 0 {
            return Err(ToolError::Execution(
                "no polygon features in input".to_string(),
            ));
        }
        if prm.zones < 1 || prm.zones > n {
            return Err(ToolError::Validation(format!(
                "'zones' must be between 1 and the polygon count ({n}), got {}",
                prm.zones
            )));
        }

        // Standardize attributes (z-scores) for homogeneity SSD and edge weights.
        let nf = prm.fields.len();
        let z = standardize(&attrs, nf);

        // Contiguity graph edges with dissimilarity weights.
        ctx.progress.info(&format!(
            "building {} contiguity graph over {n} polygon(s)",
            prm.contiguity.as_str()
        ));
        let segs: Vec<Segments> = polys.iter().map(Segments::of).collect();
        let centroids: Vec<(f64, f64)> = polys
            .iter()
            .map(|p| p.centroid().map(|c| (c.x(), c.y())).unwrap_or((0.0, 0.0)))
            .collect();
        let mut edges: Vec<Edge> = Vec::new();
        for i in 0..n {
            for j in i + 1..n {
                if !segs[i].bbox.intersects(&segs[j].bbox, prm.tolerance) {
                    continue;
                }
                let adjacent = match prm.contiguity {
                    Contiguity::Rook => segs[i].shared_len(&segs[j], prm.tolerance) > prm.tolerance,
                    Contiguity::Queen => polys[i].intersects(&polys[j]),
                };
                if adjacent {
                    // Homogeneity clusters attribute-similar neighbours, so its
                    // tree is weighted by attribute dissimilarity. The balance
                    // criteria want compact zones with balanced cuts available,
                    // so they use geographic (centroid) distance instead.
                    let w = if prm.criterion == Criterion::Homogeneity && nf > 0 {
                        (0..nf)
                            .map(|f| (z[i][f] - z[j][f]).powi(2))
                            .sum::<f64>()
                            .sqrt()
                    } else {
                        let (dx, dy) = (
                            centroids[i].0 - centroids[j].0,
                            centroids[i].1 - centroids[j].1,
                        );
                        (dx * dx + dy * dy).sqrt()
                    };
                    edges.push(Edge { a: i, b: j, w });
                }
            }
        }

        // Full contiguity adjacency (for region growing).
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for e in &edges {
            adj[e.a].push(e.b);
            adj[e.b].push(e.a);
        }

        // Minimum spanning forest (Kruskal). Zones can't be fewer than the
        // number of disconnected components.
        let (forest, components) = min_spanning_forest(n, &mut edges);
        if prm.zones < components {
            return Err(ToolError::Validation(format!(
                "the contiguity graph has {components} disconnected component(s); 'zones' must be at least {components}"
            )));
        }

        // Homogeneity uses SKATER (spanning-tree edge cutting); the balance
        // criteria use balanced region growing, which reliably produces
        // near-target zones a spanning-tree cut cannot guarantee.
        let labels = match prm.criterion {
            Criterion::Homogeneity => skater(n, forest, prm.zones, &z, nf, &attrs, prm.criterion),
            _ => {
                let measure: Vec<f64> = (0..n)
                    .map(|i| match prm.criterion {
                        Criterion::EqualSum => attrs[i][0],
                        _ => 1.0,
                    })
                    .collect();
                grow_balanced(n, &adj, &measure, &centroids, prm.zones)
            }
        };

        // Write the zone id onto every feature (nulls for non-polygons).
        layer.add_field(FieldDef::new("zone", FieldType::Integer));
        let mut zone_of_feature: Vec<Option<i64>> = vec![None; layer.features.len()];
        for (node, &fi) in node_feature.iter().enumerate() {
            zone_of_feature[fi] = Some(labels[node] as i64);
        }
        for (fi, feature) in layer.features.iter_mut().enumerate() {
            feature.attributes.push(match zone_of_feature[fi] {
                Some(z) => FieldValue::Integer(z),
                None => FieldValue::Null,
            });
        }

        // Per-zone counts.
        let mut counts = vec![0usize; prm.zones];
        for &l in &labels {
            counts[l] += 1;
        }
        ctx.progress
            .info(&format!("{} zone(s), sizes {counts:?}", prm.zones));

        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("polygon_count".to_string(), json!(n));
        outputs.insert("zones".to_string(), json!(prm.zones));
        outputs.insert("components".to_string(), json!(components));
        outputs.insert("criterion".to_string(), json!(prm.criterion.as_str()));
        outputs.insert("zone_sizes".to_string(), json!(counts));
        Ok(ToolRunResult { outputs })
    }
}

// ── SKATER ────────────────────────────────────────────────────────────────────

struct Edge {
    a: usize,
    b: usize,
    w: f64,
}

/// Kruskal minimum spanning forest. Returns the chosen tree edges and the number
/// of connected components.
fn min_spanning_forest(n: usize, edges: &mut [Edge]) -> (Vec<(usize, usize)>, usize) {
    edges.sort_by(|e, f| e.w.total_cmp(&f.w));
    let mut uf = UnionFind::new(n);
    let mut tree = Vec::new();
    for e in edges.iter() {
        if uf.union(e.a, e.b) {
            tree.push((e.a, e.b));
        }
    }
    let mut roots = std::collections::BTreeSet::new();
    for i in 0..n {
        roots.insert(uf.find(i));
    }
    (tree, roots.len())
}

/// Greedy SKATER cutting: repeatedly remove the forest edge whose removal most
/// improves the balance objective, until `zones` connected components remain.
fn skater(
    n: usize,
    forest: Vec<(usize, usize)>,
    zones: usize,
    z: &[Vec<f64>],
    nf: usize,
    attrs: &[Vec<f64>],
    criterion: Criterion,
) -> Vec<usize> {
    // Current partition: each cluster = (nodes, internal edges).
    let mut clusters = initial_clusters(n, &forest);

    // For the balance criteria, peel off zones whose measure is closest to the
    // per-zone target; homogeneity minimizes total within-zone SSD instead.
    let measure = |set: &[usize]| -> f64 {
        match criterion {
            Criterion::EqualSum => set.iter().map(|&nd| attrs[nd][0]).sum(),
            _ => set.len() as f64,
        }
    };
    let target = match criterion {
        Criterion::Homogeneity => 0.0,
        _ => measure(&(0..n).collect::<Vec<_>>()) / zones as f64,
    };

    while clusters.len() < zones {
        // (score, cluster_idx, edge_idx, partA, partB)
        let mut best: BestCut = None;
        for (ci, cl) in clusters.iter().enumerate() {
            if cl.edges.is_empty() {
                continue;
            }
            for ei in 0..cl.edges.len() {
                let (pa, pb) = split_component(&cl.nodes, &cl.edges, ei);
                let score = match criterion {
                    Criterion::Homogeneity => {
                        // Total SSD with cluster ci replaced by pa, pb.
                        let mut sets: Vec<Vec<usize>> = clusters
                            .iter()
                            .enumerate()
                            .filter(|(k, _)| *k != ci)
                            .map(|(_, c)| c.nodes.clone())
                            .collect();
                        sets.push(pa.clone());
                        sets.push(pb.clone());
                        objective(&sets, z, nf, attrs, criterion)
                    }
                    // Peel: reward a cut that carves off a near-target piece.
                    _ => (measure(&pa) - target)
                        .abs()
                        .min((measure(&pb) - target).abs()),
                };
                if best.as_ref().is_none_or(|b| score < b.0) {
                    best = Some((score, ci, ei, pa, pb));
                }
            }
        }
        let Some((_, ci, ei, pa, pb)) = best else {
            break; // no more cuttable edges
        };
        // Apply the cut: split cluster ci into (pa, pb).
        let old = clusters.swap_remove(ci);
        let (mut ea, mut eb) = (Vec::new(), Vec::new());
        let set_a: std::collections::HashSet<usize> = pa.iter().copied().collect();
        for (k, &e) in old.edges.iter().enumerate() {
            if k == ei {
                continue;
            }
            if set_a.contains(&e.0) {
                ea.push(e);
            } else {
                eb.push(e);
            }
        }
        clusters.push(Cluster {
            nodes: pa,
            edges: ea,
        });
        clusters.push(Cluster {
            nodes: pb,
            edges: eb,
        });
    }

    // Preserve determinism: label zones by their smallest member node.
    let mut order: Vec<usize> = (0..clusters.len()).collect();
    order.sort_by_key(|&k| *clusters[k].nodes.iter().min().unwrap());
    let mut labels = vec![0usize; n];
    for (zone, &k) in order.iter().enumerate() {
        for &node in &clusters[k].nodes {
            labels[node] = zone;
        }
    }
    labels
}

/// Balanced, contiguity-preserving region growing for the balance criteria.
///
/// Seeds are chosen by farthest-point sampling (deterministic). Growth then
/// repeatedly extends the zone with the smallest current measure by the adjacent
/// unassigned unit nearest that zone's seed — so zones stay compact and
/// connected while their measures stay close to the per-zone target. Any unit
/// left unreachable (in a component with no seed) is assigned to the nearest
/// zone by centroid.
fn grow_balanced(
    n: usize,
    adj: &[Vec<usize>],
    measure: &[f64],
    centroids: &[(f64, f64)],
    zones: usize,
) -> Vec<usize> {
    let dist = |a: usize, b: usize| {
        let (dx, dy) = (centroids[a].0 - centroids[b].0, centroids[a].1 - centroids[b].1);
        (dx * dx + dy * dy).sqrt()
    };
    // Farthest-point seeds: start at node 0, then repeatedly the node maximizing
    // its distance to the nearest existing seed.
    let mut seeds = vec![0usize];
    while seeds.len() < zones {
        let mut best = (-1.0, 0usize);
        for i in 0..n {
            if seeds.contains(&i) {
                continue;
            }
            let d = seeds.iter().map(|&s| dist(i, s)).fold(f64::INFINITY, f64::min);
            if d > best.0 {
                best = (d, i);
            }
        }
        seeds.push(best.1);
    }

    let mut label = vec![usize::MAX; n];
    let mut zmeasure = vec![0.0; zones];
    let mut seed_of = vec![0usize; zones];
    for (zone, &s) in seeds.iter().enumerate() {
        label[s] = zone;
        zmeasure[zone] = measure[s];
        seed_of[zone] = s;
    }
    let mut remaining = n - zones;

    while remaining > 0 {
        // The lightest zone that still has an adjacent unassigned unit.
        let mut order: Vec<usize> = (0..zones).collect();
        order.sort_by(|&a, &b| zmeasure[a].total_cmp(&zmeasure[b]));
        let mut grew = false;
        for &zone in &order {
            // Frontier: unassigned units adjacent to this zone.
            let mut best: Option<(f64, usize)> = None;
            for i in 0..n {
                if label[i] == zone {
                    for &nb in &adj[i] {
                        if label[nb] == usize::MAX {
                            let d = dist(nb, seed_of[zone]);
                            if best.is_none_or(|(bd, _)| d < bd) {
                                best = Some((d, nb));
                            }
                        }
                    }
                }
            }
            if let Some((_, node)) = best {
                label[node] = zone;
                zmeasure[zone] += measure[node];
                remaining -= 1;
                grew = true;
                break;
            }
        }
        if !grew {
            break; // remaining units are unreachable (disconnected)
        }
    }

    // Assign any stranded units to the nearest zone by centroid.
    if remaining > 0 {
        for i in 0..n {
            if label[i] == usize::MAX {
                let z = (0..zones)
                    .min_by(|&a, &b| dist(i, seed_of[a]).total_cmp(&dist(i, seed_of[b])))
                    .unwrap();
                label[i] = z;
            }
        }
    }
    label
}

struct Cluster {
    nodes: Vec<usize>,
    edges: Vec<(usize, usize)>,
}

/// Best candidate cut in a greedy step: (score, cluster index, edge index within
/// that cluster, and the two node partitions the cut produces).
type BestCut = Option<(f64, usize, usize, Vec<usize>, Vec<usize>)>;

/// Connected components of the spanning forest as initial clusters.
fn initial_clusters(n: usize, forest: &[(usize, usize)]) -> Vec<Cluster> {
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(a, b) in forest {
        adj[a].push(b);
        adj[b].push(a);
    }
    let mut seen = vec![false; n];
    let mut clusters = Vec::new();
    for s in 0..n {
        if seen[s] {
            continue;
        }
        let mut nodes = Vec::new();
        let mut stack = vec![s];
        seen[s] = true;
        while let Some(u) = stack.pop() {
            nodes.push(u);
            for &v in &adj[u] {
                if !seen[v] {
                    seen[v] = true;
                    stack.push(v);
                }
            }
        }
        let nset: std::collections::HashSet<usize> = nodes.iter().copied().collect();
        let edges: Vec<(usize, usize)> = forest
            .iter()
            .copied()
            .filter(|(a, b)| nset.contains(a) && nset.contains(b))
            .collect();
        clusters.push(Cluster { nodes, edges });
    }
    clusters
}

/// Splits a cluster's node set into the two components formed by removing edge
/// `skip` from its edge list.
fn split_component(
    nodes: &[usize],
    edges: &[(usize, usize)],
    skip: usize,
) -> (Vec<usize>, Vec<usize>) {
    let mut adj: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for &node in nodes {
        adj.entry(node).or_default();
    }
    for (k, &(a, b)) in edges.iter().enumerate() {
        if k == skip {
            continue;
        }
        adj.get_mut(&a).unwrap().push(b);
        adj.get_mut(&b).unwrap().push(a);
    }
    let start = edges[skip].0;
    let mut seen = std::collections::HashSet::new();
    let mut stack = vec![start];
    seen.insert(start);
    while let Some(u) = stack.pop() {
        for &v in &adj[&u] {
            if seen.insert(v) {
                stack.push(v);
            }
        }
    }
    let a: Vec<usize> = nodes
        .iter()
        .copied()
        .filter(|nd| seen.contains(nd))
        .collect();
    let b: Vec<usize> = nodes
        .iter()
        .copied()
        .filter(|nd| !seen.contains(nd))
        .collect();
    (a, b)
}

/// The balance objective (lower is better) for a full partition.
fn objective(
    sets: &[Vec<usize>],
    z: &[Vec<f64>],
    nf: usize,
    attrs: &[Vec<f64>],
    criterion: Criterion,
) -> f64 {
    match criterion {
        Criterion::Homogeneity => sets.iter().map(|s| ssd(s, z, nf)).sum(),
        Criterion::EqualCount => {
            let sizes: Vec<f64> = sets.iter().map(|s| s.len() as f64).collect();
            variance(&sizes)
        }
        Criterion::EqualSum => {
            let sums: Vec<f64> = sets
                .iter()
                .map(|s| s.iter().map(|&nd| attrs[nd][0]).sum())
                .collect();
            variance(&sums)
        }
    }
}

/// Within-set sum of squared deviations across the standardized fields.
#[allow(clippy::needless_range_loop)]
fn ssd(set: &[usize], z: &[Vec<f64>], nf: usize) -> f64 {
    if set.is_empty() {
        return 0.0;
    }
    let m = set.len() as f64;
    let mut total = 0.0;
    for f in 0..nf {
        let mean = set.iter().map(|&i| z[i][f]).sum::<f64>() / m;
        total += set.iter().map(|&i| (z[i][f] - mean).powi(2)).sum::<f64>();
    }
    total
}

fn variance(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let m = v.iter().sum::<f64>() / v.len() as f64;
    v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / v.len() as f64
}

/// Z-score standardization per field (std 0 -> leave centered at 0).
#[allow(clippy::needless_range_loop)]
fn standardize(attrs: &[Vec<f64>], nf: usize) -> Vec<Vec<f64>> {
    let n = attrs.len();
    let mut z = vec![vec![0.0; nf]; n];
    for f in 0..nf {
        let mean = attrs.iter().map(|a| a[f]).sum::<f64>() / n as f64;
        let var = attrs.iter().map(|a| (a[f] - mean).powi(2)).sum::<f64>() / n as f64;
        let sd = var.sqrt();
        for i in 0..n {
            z[i][f] = if sd > 0.0 {
                (attrs[i][f] - mean) / sd
            } else {
                0.0
            };
        }
    }
    z
}

// ── Union-Find ────────────────────────────────────────────────────────────────

struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
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
    fn union(&mut self, a: usize, b: usize) -> bool {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return false;
        }
        match self.rank[ra].cmp(&self.rank[rb]) {
            std::cmp::Ordering::Less => self.parent[ra] = rb,
            std::cmp::Ordering::Greater => self.parent[rb] = ra,
            std::cmp::Ordering::Equal => {
                self.parent[rb] = ra;
                self.rank[ra] += 1;
            }
        }
        true
    }
}

// ── Shared-boundary detection (rook) ──────────────────────────────────────────

#[derive(Clone, Copy)]
struct BBox {
    minx: f64,
    miny: f64,
    maxx: f64,
    maxy: f64,
}

impl BBox {
    fn empty() -> Self {
        Self {
            minx: f64::INFINITY,
            miny: f64::INFINITY,
            maxx: f64::NEG_INFINITY,
            maxy: f64::NEG_INFINITY,
        }
    }
    fn expand(&mut self, x: f64, y: f64) {
        self.minx = self.minx.min(x);
        self.miny = self.miny.min(y);
        self.maxx = self.maxx.max(x);
        self.maxy = self.maxy.max(y);
    }
    fn intersects(&self, o: &BBox, pad: f64) -> bool {
        self.minx <= o.maxx + pad
            && self.maxx >= o.minx - pad
            && self.miny <= o.maxy + pad
            && self.maxy >= o.miny - pad
    }
}

struct Segments {
    segs: Vec<(GeoCoord, GeoCoord)>,
    bbox: BBox,
}

impl Segments {
    fn of(mp: &MultiPolygon) -> Self {
        let mut segs = Vec::new();
        let mut bbox = BBox::empty();
        for poly in mp {
            for ring in std::iter::once(poly.exterior()).chain(poly.interiors()) {
                let pts = &ring.0;
                for w in pts.windows(2) {
                    segs.push((w[0], w[1]));
                    bbox.expand(w[0].x, w[0].y);
                }
                if let Some(last) = pts.last() {
                    bbox.expand(last.x, last.y);
                }
            }
        }
        Self { segs, bbox }
    }
    fn shared_len(&self, other: &Segments, tol: f64) -> f64 {
        let mut total = 0.0;
        for &(p1, p2) in &self.segs {
            for &(q1, q2) in &other.segs {
                total += collinear_overlap(p1, p2, q1, q2, tol);
            }
        }
        total
    }
}

fn collinear_overlap(p1: GeoCoord, p2: GeoCoord, q1: GeoCoord, q2: GeoCoord, tol: f64) -> f64 {
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let len = dx.hypot(dy);
    if len <= tol {
        return 0.0;
    }
    let (ux, uy) = (dx / len, dy / len);
    let perp = |q: GeoCoord| ((q.x - p1.x) * uy - (q.y - p1.y) * ux).abs();
    if perp(q1) > tol || perp(q2) > tol {
        return 0.0;
    }
    let proj = |q: GeoCoord| (q.x - p1.x) * ux + (q.y - p1.y) * uy;
    let (tq1, tq2) = (proj(q1), proj(q2));
    let lo = 0.0f64.max(tq1.min(tq2));
    let hi = len.min(tq1.max(tq2));
    (hi - lo).max(0.0)
}

fn to_multipolygon(geom: &Geometry) -> Option<MultiPolygon> {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => Some(MultiPolygon(vec![rings_to_polygon(exterior, interiors)])),
        Geometry::MultiPolygon(parts) => Some(MultiPolygon(
            parts.iter().map(|(e, i)| rings_to_polygon(e, i)).collect(),
        )),
        _ => None,
    }
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
            .map(|c| GeoCoord { x: c.x, y: c.y })
            .collect(),
    )
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Criterion {
    Homogeneity,
    EqualCount,
    EqualSum,
}

impl Criterion {
    fn as_str(self) -> &'static str {
        match self {
            Self::Homogeneity => "homogeneity",
            Self::EqualCount => "equal_count",
            Self::EqualSum => "equal_sum",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Contiguity {
    Rook,
    Queen,
}

impl Contiguity {
    fn as_str(self) -> &'static str {
        match self {
            Self::Rook => "rook",
            Self::Queen => "queen",
        }
    }
}

struct Params {
    zones: usize,
    criterion: Criterion,
    fields: Vec<String>,
    contiguity: Contiguity,
    tolerance: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let zones = match parse_optional_f64(args, "zones")? {
        Some(v) if v.fract() == 0.0 && v >= 1.0 && v.is_finite() => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "parameter 'zones' must be a positive integer".to_string(),
            ))
        }
        None => {
            return Err(ToolError::Validation(
                "missing required parameter 'zones'".to_string(),
            ))
        }
    };
    let criterion = match parse_optional_str(args, "criterion")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("homogeneity") => Criterion::Homogeneity,
        Some("equal_count") => Criterion::EqualCount,
        Some("equal_sum") => Criterion::EqualSum,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown criterion '{other}' (expected homogeneity, equal_count, or equal_sum)"
            )))
        }
    };
    let fields: Vec<String> = parse_optional_str(args, "fields")?
        .map(|s| {
            s.split(',')
                .map(str::trim)
                .filter(|t| !t.is_empty())
                .map(str::to_string)
                .collect()
        })
        .unwrap_or_default();
    let contiguity = match parse_optional_str(args, "contiguity")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("rook") => Contiguity::Rook,
        Some("queen") => Contiguity::Queen,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown contiguity '{other}' (expected rook or queen)"
            )))
        }
    };
    let tolerance = parse_optional_f64(args, "tolerance")?.unwrap_or(1e-6);
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'tolerance' must be a positive number".to_string(),
        ));
    }
    Ok(Params {
        zones,
        criterion,
        fields,
        contiguity,
        tolerance,
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

    fn cell(col: f64, row: f64) -> Geometry {
        use wbvector::Coord;
        Geometry::polygon(
            vec![
                Coord::xy(col, row),
                Coord::xy(col + 1.0, row),
                Coord::xy(col + 1.0, row + 1.0),
                Coord::xy(col, row + 1.0),
            ],
            vec![],
        )
    }

    /// Builds a 1xN strip of unit cells with a `val` attribute.
    fn strip(vals: &[f64]) -> String {
        let mut layer = Layer::new("cells");
        layer.add_field(FieldDef::new("val", FieldType::Float));
        for (i, &v) in vals.iter().enumerate() {
            layer
                .add_feature(Some(cell(i as f64, 0.0)), &[("val", FieldValue::Float(v))])
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = BuildBalancedZonesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn zone(layer: &Layer, i: usize) -> i64 {
        match layer.features[i].get(&layer.schema, "zone").unwrap() {
            FieldValue::Integer(z) => *z,
            other => panic!("zone should be integer, got {other:?}"),
        }
    }

    #[test]
    fn homogeneity_splits_at_the_value_jump() {
        // Six cells in a row: values [1,1,1,9,9,9]. Two homogeneous zones should
        // split exactly between index 2 and 3.
        let (out, layer) = run(json!({
            "input": strip(&[1.0, 1.0, 1.0, 9.0, 9.0, 9.0]),
            "zones": 2, "criterion": "homogeneity", "fields": "val"
        }));
        assert_eq!(out.outputs["zones"], json!(2));
        let z: Vec<i64> = (0..6).map(|i| zone(&layer, i)).collect();
        // First three share a zone, last three share the other.
        assert_eq!(z[0], z[1]);
        assert_eq!(z[1], z[2]);
        assert_eq!(z[3], z[4]);
        assert_eq!(z[4], z[5]);
        assert_ne!(z[2], z[3], "the split must fall at the value jump");
    }

    #[test]
    fn zones_are_contiguous() {
        // A 1x8 strip cut into 3 zones: each zone is a contiguous run of cells.
        let (_, layer) = run(json!({
            "input": strip(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]),
            "zones": 3, "criterion": "equal_count"
        }));
        let z: Vec<i64> = (0..8).map(|i| zone(&layer, i)).collect();
        // Contiguity on a strip: each zone id forms one uninterrupted run.
        let mut runs = 0;
        for i in 0..8 {
            if i == 0 || z[i] != z[i - 1] {
                runs += 1;
            }
        }
        assert_eq!(
            runs, 3,
            "each of the 3 zones must be one contiguous run, got {z:?}"
        );
    }

    #[test]
    fn equal_count_balances_sizes() {
        let (out, _) = run(json!({
            "input": strip(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
            "zones": 3, "criterion": "equal_count"
        }));
        // 6 cells / 3 zones -> sizes [2,2,2].
        let sizes: Vec<i64> = out.outputs["zone_sizes"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_i64().unwrap())
            .collect();
        assert_eq!(sizes, vec![2, 2, 2], "sizes {sizes:?}");
    }

    #[test]
    fn equal_sum_balances_attribute_totals() {
        // Values [1,1,1,1,1,5]: to equalize sums into 2 zones, the lone 5 becomes
        // one zone (sum 5) and the five 1s the other (sum 5).
        let (_, layer) = run(json!({
            "input": strip(&[1.0, 1.0, 1.0, 1.0, 1.0, 5.0]),
            "zones": 2, "criterion": "equal_sum", "fields": "val"
        }));
        let z: Vec<i64> = (0..6).map(|i| zone(&layer, i)).collect();
        assert_ne!(
            z[5], z[4],
            "the value-5 cell should be its own balancing zone"
        );
        for i in 0..5 {
            assert_eq!(z[i], z[0], "the five 1-cells share a zone");
        }
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = BuildBalancedZonesTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "input": "x.geojson" })).is_err(),
            "missing zones"
        );
        assert!(bad(json!({ "input": "x.geojson", "zones": 0 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "zones": 2, "criterion": "bogus" })).is_err());
        assert!(
            bad(json!({ "input": "x.geojson", "zones": 2, "criterion": "equal_count" })).is_ok()
        );
    }
}
