//! GeoLibre tool: HDBSCAN* hierarchical density-based clustering.
//!
//! Pure-Rust counterpart of the HDBSCAN option in ArcGIS Pro's *Density-based
//! Clustering* (Spatial Statistics). The bundled suite has `dbscan` and k-means
//! only; HDBSCAN is the current default for point clustering because real point
//! data (crime, wildlife, POIs) is variable-density and DBSCAN's single epsilon
//! is brittle.
//!
//! The full HDBSCAN* pipeline:
//!
//! 1. **Core distance** of each point — the distance to its `min_samples`-th
//!    nearest neighbour (a local density estimate).
//! 2. **Mutual reachability** distance `max(core_a, core_b, d(a,b))` — pushes
//!    sparse points apart.
//! 3. **Minimum spanning tree** of the mutual-reachability graph (Prim's).
//! 4. **Single-linkage dendrogram** from the sorted MST edges (union-find).
//! 5. **Condensed tree**: walking the dendrogram top-down, a split where both
//!    sides keep at least `min_cluster_size` points is a real split; smaller
//!    sides "fall out" as their points leave the parent cluster.
//! 6. **Excess-of-mass** cluster selection — pick the set of non-nested clusters
//!    with the greatest total stability.
//!
//! Output copies the input points and adds `cluster_id` (−1 = noise),
//! `probability` (membership strength, 0-1), and `outlier_score` (GLOSH-style,
//! higher = more outlying). Deterministic; O(n²) so it suits moderate point
//! counts. Use a projected CRS (distances in its units).

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct HdbscanTool;

impl Tool for HdbscanTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "hdbscan",
            display_name: "HDBSCAN Clustering",
            summary: "Hierarchical density-based clustering (HDBSCAN*): mutual-reachability MST, condensed cluster tree, and excess-of-mass selection, with per-point membership probability and outlier score — handles variable density with no epsilon, like the HDBSCAN option of ArcGIS Density-based Clustering.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer (projected CRS; distances in its units).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (a copy of the input with cluster fields). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_cluster_size",
                    description: "Smallest grouping considered a cluster (>= 2). Default 5.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_samples",
                    description: "Neighbourhood size for the core-distance density estimate (higher = more points called noise). Default: min_cluster_size.",
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

        // Representative point per feature; remember which feature it came from.
        let mut pts: Vec<(f64, f64)> = Vec::new();
        let mut feat_of: Vec<usize> = Vec::new();
        for (fi, feature) in layer.features.iter().enumerate() {
            if let Some((x, y)) = feature.geometry.as_ref().and_then(point_xy) {
                pts.push((x, y));
                feat_of.push(fi);
            }
        }
        let n = pts.len();
        let min_samples = prm.min_samples.unwrap_or(prm.min_cluster_size).max(1);
        if n < prm.min_cluster_size.max(min_samples + 1) {
            return Err(ToolError::Execution(format!(
                "need more than min_cluster_size/min_samples points, found {n}"
            )));
        }

        ctx.progress
            .info(&format!("clustering {n} point(s) with HDBSCAN*"));

        let result = hdbscan(&pts, prm.min_cluster_size, min_samples);

        // Write cluster fields onto a copy of the input.
        layer.add_field(FieldDef::new("cluster_id", FieldType::Integer));
        layer.add_field(FieldDef::new("probability", FieldType::Float));
        layer.add_field(FieldDef::new("outlier_score", FieldType::Float));
        // Default values for features without geometry (not clustered).
        let mut per_feat: Vec<(i64, f64, f64)> = vec![(-1, 0.0, 0.0); layer.len()];
        for i in 0..n {
            per_feat[feat_of[i]] = (result.labels[i], result.probabilities[i], result.outlier[i]);
        }
        for (fi, feature) in layer.features.iter_mut().enumerate() {
            let (c, p, o) = per_feat[fi];
            feature.attributes.push(FieldValue::Integer(c));
            feature.attributes.push(FieldValue::Float(p));
            feature.attributes.push(FieldValue::Float(o));
        }

        let n_clusters = result
            .labels
            .iter()
            .filter(|&&l| l >= 0)
            .copied()
            .max()
            .map(|m| (m + 1) as usize)
            .unwrap_or(0);
        let n_noise = result.labels.iter().filter(|&&l| l < 0).count();
        ctx.progress.info(&format!(
            "{n_clusters} cluster(s), {n_noise} noise point(s)"
        ));

        let feature_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("cluster_count".to_string(), json!(n_clusters));
        outputs.insert("noise_count".to_string(), json!(n_noise));
        Ok(ToolRunResult { outputs })
    }
}

// ── HDBSCAN* ─────────────────────────────────────────────────────────────────

struct ClusterResult {
    labels: Vec<i64>,
    probabilities: Vec<f64>,
    outlier: Vec<f64>,
}

fn hdbscan(pts: &[(f64, f64)], min_cluster_size: usize, min_samples: usize) -> ClusterResult {
    let n = pts.len();
    // Pairwise distances.
    let dist = |a: usize, b: usize| -> f64 {
        let (dx, dy) = (pts[a].0 - pts[b].0, pts[a].1 - pts[b].1);
        dx.hypot(dy)
    };

    // Core distance: distance to the min_samples-th nearest neighbour.
    let mut core = vec![0.0f64; n];
    for (i, core_i) in core.iter_mut().enumerate() {
        let mut ds: Vec<f64> = (0..n).filter(|&j| j != i).map(|j| dist(i, j)).collect();
        ds.sort_by(f64::total_cmp);
        let k = (min_samples - 1).min(ds.len().saturating_sub(1));
        *core_i = ds.get(k).copied().unwrap_or(0.0);
    }

    let mreach = |a: usize, b: usize| -> f64 { dist(a, b).max(core[a]).max(core[b]) };

    // Prim's MST on the mutual-reachability graph.
    let mut in_tree = vec![false; n];
    let mut best = vec![f64::INFINITY; n];
    let mut parent = vec![usize::MAX; n];
    best[0] = 0.0;
    let mut edges: Vec<(f64, usize, usize)> = Vec::with_capacity(n.saturating_sub(1));
    for _ in 0..n {
        let mut u = usize::MAX;
        let mut bu = f64::INFINITY;
        for v in 0..n {
            if !in_tree[v] && best[v] < bu {
                bu = best[v];
                u = v;
            }
        }
        if u == usize::MAX {
            break;
        }
        in_tree[u] = true;
        if parent[u] != usize::MAX {
            edges.push((bu, parent[u], u));
        }
        for v in 0..n {
            if !in_tree[v] {
                let w = mreach(u, v);
                if w < best[v] {
                    best[v] = w;
                    parent[v] = u;
                }
            }
        }
    }
    edges.sort_by(|a, b| a.0.total_cmp(&b.0));

    // Single-linkage dendrogram via union-find. Node ids: 0..n leaves, then
    // internal nodes n.. ; each MST edge merges two components into a new node.
    let mut uf = UnionFind::new(n);
    let mut node_of = (0..n).collect::<Vec<usize>>(); // component root leaf/internal id
    let mut size = vec![1usize; 2 * n];
    let mut children: Vec<(usize, usize)> = vec![(usize::MAX, usize::MAX); 2 * n];
    let mut lambda_birth = vec![0.0f64; 2 * n]; // birth lambda of the node as a child
    let mut next = n;
    // For each merge, record the new node's formation distance.
    let mut node_dist = vec![0.0f64; 2 * n];
    for &(w, a, b) in &edges {
        let ra = uf.find(a);
        let rb = uf.find(b);
        let na = node_of[ra];
        let nb = node_of[rb];
        let new_id = next;
        next += 1;
        children[new_id] = (na, nb);
        node_dist[new_id] = w;
        size[new_id] = size[na] + size[nb];
        let lam = if w > 0.0 { 1.0 / w } else { f64::INFINITY };
        lambda_birth[na] = lam;
        lambda_birth[nb] = lam;
        let r = uf.union(ra, rb);
        node_of[r] = new_id;
    }
    let root = next - 1;

    // ── Condense the tree ────────────────────────────────────────────────────
    // Walk from the root; a child bigger than min_cluster_size stays a cluster,
    // otherwise its points fall out at the split's lambda.
    // condensed cluster nodes are the internal dendrogram nodes that survive.
    let mut point_lambda = vec![0.0f64; n]; // lambda at which each point leaves its cluster
    let mut point_cluster = vec![root; n]; // the condensed cluster each point ends in
                                           // cluster stability accumulation and structure
    let mut cluster_lambda_birth: BTreeMap<usize, f64> = BTreeMap::new();
    let mut cluster_children: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    let mut cluster_parent: BTreeMap<usize, usize> = BTreeMap::new();
    let mut stability: BTreeMap<usize, f64> = BTreeMap::new();
    cluster_lambda_birth.insert(root, 0.0);
    stability.insert(root, 0.0);

    // Iterative descent. Each item: (dendro node, condensed cluster it belongs to).
    let mut stack = vec![(root, root)];
    while let Some((node, cluster)) = stack.pop() {
        if node < n {
            // Leaf: it stays in `cluster` down to λ = ∞ (its own birth lambda).
            let lam = lambda_birth[node].min(1e18);
            point_lambda[node] = lam;
            point_cluster[node] = cluster;
            *stability.get_mut(&cluster).unwrap() += lam - cluster_lambda_birth[&cluster];
            continue;
        }
        let (a, b) = children[node];
        let lam_split = if node_dist[node] > 0.0 {
            1.0 / node_dist[node]
        } else {
            f64::INFINITY
        };
        let a_big = size[a] >= min_cluster_size;
        let b_big = size[b] >= min_cluster_size;
        match (a_big, b_big) {
            (true, true) => {
                // Real split: both children become new condensed clusters.
                for &child in &[a, b] {
                    cluster_lambda_birth.insert(child, lam_split);
                    cluster_parent.insert(child, cluster);
                    cluster_children.entry(cluster).or_default().push(child);
                    stability.insert(child, 0.0);
                    stack.push((child, child));
                }
            }
            (true, false) => {
                // b's points fall out of `cluster` at lam_split; a continues it.
                fall_out(
                    b,
                    lam_split,
                    cluster,
                    &children,
                    n,
                    &mut point_lambda,
                    &mut point_cluster,
                    &mut stability,
                    &cluster_lambda_birth,
                );
                stack.push((a, cluster));
            }
            (false, true) => {
                fall_out(
                    a,
                    lam_split,
                    cluster,
                    &children,
                    n,
                    &mut point_lambda,
                    &mut point_cluster,
                    &mut stability,
                    &cluster_lambda_birth,
                );
                stack.push((b, cluster));
            }
            (false, false) => {
                // Both small: all points fall out of `cluster` here.
                fall_out(
                    a,
                    lam_split,
                    cluster,
                    &children,
                    n,
                    &mut point_lambda,
                    &mut point_cluster,
                    &mut stability,
                    &cluster_lambda_birth,
                );
                fall_out(
                    b,
                    lam_split,
                    cluster,
                    &children,
                    n,
                    &mut point_lambda,
                    &mut point_cluster,
                    &mut stability,
                    &cluster_lambda_birth,
                );
            }
        }
    }

    // ── Excess of mass selection ─────────────────────────────────────────────
    // Process clusters bottom-up: select a cluster if its own stability exceeds
    // the sum of its selected descendants; else propagate the descendants' sum.
    let mut cluster_ids: Vec<usize> = stability.keys().copied().collect();
    // Deeper clusters (smaller size) first — larger dendro id ≈ higher up, so
    // sort by size ascending processes children before parents.
    cluster_ids.sort_by_key(|&c| size[c]);
    let mut selected: std::collections::HashSet<usize> = Default::default();
    let mut prop_stability: BTreeMap<usize, f64> = BTreeMap::new();
    for &c in &cluster_ids {
        if c == root {
            continue; // the root is never itself a cluster
        }
        let own = stability[&c];
        let child_sum: f64 = cluster_children
            .get(&c)
            .map(|ch| {
                ch.iter()
                    .map(|k| prop_stability.get(k).copied().unwrap_or(0.0))
                    .sum()
            })
            .unwrap_or(0.0);
        if own >= child_sum {
            prop_stability.insert(c, own);
            // Deselect descendants, select c.
            if let Some(ch) = cluster_children.get(&c) {
                let mut dstack = ch.clone();
                while let Some(d) = dstack.pop() {
                    selected.remove(&d);
                    if let Some(gc) = cluster_children.get(&d) {
                        dstack.extend(gc.iter().copied());
                    }
                }
            }
            selected.insert(c);
        } else {
            prop_stability.insert(c, child_sum);
        }
    }

    // ── Labels, probabilities, outlier scores ────────────────────────────────
    // Map each point to the nearest selected ancestor cluster (or noise).
    let mut sel_index: BTreeMap<usize, i64> = BTreeMap::new();
    for (i, &c) in selected.iter().enumerate() {
        sel_index.insert(c, i as i64);
    }
    let mut labels = vec![-1i64; n];
    let mut probabilities = vec![0.0f64; n];
    let mut lambda_max: BTreeMap<usize, f64> = BTreeMap::new();
    for i in 0..n {
        // Walk up from the point's condensed cluster to a selected one.
        let mut c = point_cluster[i];
        loop {
            if let Some(&lab) = sel_index.get(&c) {
                labels[i] = lab;
                let lm = lambda_max.entry(c).or_insert(0.0);
                if point_lambda[i] > *lm {
                    *lm = point_lambda[i];
                }
                break;
            }
            match cluster_parent.get(&c) {
                Some(&p) => c = p,
                None => break,
            }
        }
    }
    // Probability = point_lambda / cluster_lambda_max.
    for i in 0..n {
        if labels[i] >= 0 {
            let c = {
                // Find the selected cluster again.
                let mut cc = point_cluster[i];
                while !sel_index.contains_key(&cc) {
                    match cluster_parent.get(&cc) {
                        Some(&p) => cc = p,
                        None => break,
                    }
                }
                cc
            };
            let lm = lambda_max.get(&c).copied().unwrap_or(1.0).max(1e-12);
            probabilities[i] = (point_lambda[i] / lm).clamp(0.0, 1.0);
        }
    }
    // Outlier score (GLOSH-lite): 1 - lambda_point / lambda_max_of_its_cluster.
    let outlier: Vec<f64> = probabilities
        .iter()
        .zip(&labels)
        .map(|(p, &l)| if l >= 0 { 1.0 - p } else { 1.0 })
        .collect();

    ClusterResult {
        labels,
        probabilities,
        outlier,
    }
}

/// A whole subtree of points leaves `cluster` at `lam`.
#[allow(clippy::too_many_arguments)]
fn fall_out(
    node: usize,
    lam: f64,
    cluster: usize,
    children: &[(usize, usize)],
    n: usize,
    point_lambda: &mut [f64],
    point_cluster: &mut [usize],
    stability: &mut BTreeMap<usize, f64>,
    cluster_lambda_birth: &BTreeMap<usize, f64>,
) {
    let birth = cluster_lambda_birth[&cluster];
    let mut stack = vec![node];
    while let Some(nd) = stack.pop() {
        if nd < n {
            point_lambda[nd] = lam;
            point_cluster[nd] = cluster;
            *stability.get_mut(&cluster).unwrap() += lam - birth;
        } else {
            let (a, b) = children[nd];
            stack.push(a);
            stack.push(b);
        }
    }
}

struct UnionFind {
    parent: Vec<usize>,
}
impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
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
    fn union(&mut self, a: usize, b: usize) -> usize {
        let ra = self.find(a);
        let rb = self.find(b);
        self.parent[rb] = ra;
        ra
    }
}

// ── Geometry / parameters ────────────────────────────────────────────────────

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => {
            let n = cs.len() as f64;
            Some((
                cs.iter().map(|c| c.x).sum::<f64>() / n,
                cs.iter().map(|c| c.y).sum::<f64>() / n,
            ))
        }
        _ => None,
    }
}

struct Params {
    min_cluster_size: usize,
    min_samples: Option<usize>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let min_cluster_size = match parse_opt_usize(args, "min_cluster_size")? {
        None => 5,
        Some(v) if v >= 2 => v,
        Some(_) => {
            return Err(ToolError::Validation(
                "'min_cluster_size' must be >= 2".to_string(),
            ))
        }
    };
    let min_samples = parse_opt_usize(args, "min_samples")?;
    if let Some(v) = min_samples {
        if v < 1 {
            return Err(ToolError::Validation(
                "'min_samples' must be >= 1".to_string(),
            ));
        }
    }
    Ok(Params {
        min_cluster_size,
        min_samples,
    })
}

fn parse_opt_usize(args: &ToolArgs, key: &str) -> Result<Option<usize>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64().map(|v| v as usize)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be an integer"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be an integer"
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

    fn layer_of(pts: &[(f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        for &(x, y) in pts {
            l.add_feature(Some(Geometry::point(x, y)), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = HdbscanTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn labels(layer: &Layer) -> Vec<i64> {
        let idx = layer.schema.field_index("cluster_id").unwrap();
        layer
            .iter()
            .map(|f| f.attributes[idx].as_i64().unwrap())
            .collect()
    }

    /// Two tight, well-separated blobs -> two clusters, no noise between them.
    #[test]
    fn separates_two_blobs() {
        let mut pts = Vec::new();
        // Blob A around (0,0).
        for i in 0..10 {
            pts.push((0.1 * (i % 3) as f64, 0.1 * (i / 3) as f64));
        }
        // Blob B around (100,100).
        for i in 0..10 {
            pts.push((100.0 + 0.1 * (i % 3) as f64, 100.0 + 0.1 * (i / 3) as f64));
        }
        let (out, layer) = run(json!({ "input": layer_of(&pts), "min_cluster_size": 5 }));
        assert_eq!(
            out.outputs["cluster_count"],
            json!(2),
            "should find 2 clusters"
        );
        let lab = labels(&layer);
        // The two blobs get different labels.
        assert_ne!(lab[0], lab[15], "the two blobs must be different clusters");
        // Points within a blob share a label.
        assert!(lab[0..10].iter().all(|&l| l == lab[0] && l >= 0));
    }

    /// A far-away lone point is labelled noise.
    #[test]
    fn isolates_noise() {
        let mut pts = Vec::new();
        for i in 0..12 {
            pts.push((0.1 * (i % 4) as f64, 0.1 * (i / 4) as f64));
        }
        pts.push((1000.0, 1000.0)); // outlier
        let (_o, layer) = run(json!({ "input": layer_of(&pts), "min_cluster_size": 5 }));
        let lab = labels(&layer);
        assert_eq!(lab[12], -1, "the far outlier should be noise");
        // The outlier score of noise is 1.
        let oidx = layer.schema.field_index("outlier_score").unwrap();
        assert!((layer.features[12].attributes[oidx].as_f64().unwrap() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            HdbscanTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "min_cluster_size": 1 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "min_cluster_size": 5 })).is_ok());
    }
}
