//! GeoLibre tool: SKATER contiguity-constrained regionalization.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Spatially Constrained Multivariate
//! Clustering* (Spatial Statistics). Neither registry had contiguity-constrained
//! clustering:
//!
//!   * GeoLibre's `multivariate_clustering` is unconstrained k-means over the
//!     attribute space — its clusters are routinely scattered across the map,
//!     which makes them useless as districts or management regions.
//!   * GeoLibre's `build_balanced_zones` optimises for *balance* on a target
//!     variable, not attribute homogeneity, under a different objective.
//!   * The bundled `dbscan` / `hdbscan` / `k_means_clustering` are density- or
//!     distance-based on point geometry, not regionalization over a contiguity
//!     graph.
//!
//! SKATER is the standard method: build the contiguity graph, build a minimum
//! spanning tree whose edge weights are attribute distance, then repeatedly cut
//! the tree edge whose removal buys the largest reduction in within-cluster sum
//! of squared deviations. Every cluster is therefore connected **by
//! construction**, which is the property k-means cannot give you.
//!
//! **Scope for v1:** `contiguity_edges`, `contiguity_edges_corners` and `knn`
//! neighbourhoods are implemented. Trimmed Delaunay is not — there is no
//! Delaunay triangulator in the dependency set, and `knn` covers the case it
//! was wanted for (linking islands and other non-touching features).

use std::collections::BTreeMap;
use std::collections::HashMap;

use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Coordinate quantum for matching shared vertices between polygons.
///
/// 1e-6, matching `merge_divided_roads`. At 1e-9 the match is effectively
/// exact, so two polygons sharing a boundary whose vertices differ by ordinary
/// storage rounding (~1e-6) would not be neighbours; the graph would then split
/// and the tool would silently report extra disconnected clusters.
const SNAP: f64 = 1e-6;

pub struct SpatiallyConstrainedMultivariateClusteringTool;

impl Tool for SpatiallyConstrainedMultivariateClusteringTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "spatially_constrained_multivariate_clustering",
            display_name: "Spatially Constrained Multivariate Clustering",
            summary: "Partition features into spatially contiguous, attribute-homogeneous regions with the SKATER algorithm (minimum spanning tree over a contiguity graph, cut to maximise between-cluster variance), with optional size/balance constraints. Every cluster is connected by construction, unlike multivariate_clustering. Like ArcGIS Spatially Constrained Multivariate Clustering.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon (or point) features.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output path for the input features plus a CLUSTER_ID field. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "analysis_fields",
                    description: "Comma/semicolon-separated numeric field names to cluster on.",
                    required: true,
                },
                ToolParamSpec {
                    name: "number_of_clusters",
                    description: "Target number of clusters (default 5).",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighborhood",
                    description: "'contiguity_edges' (rook, default), 'contiguity_edges_corners' (queen) or 'knn'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "number_of_neighbors",
                    description: "Neighbour count for the 'knn' neighbourhood (default 8).",
                    required: false,
                },
                ToolParamSpec {
                    name: "constraint",
                    description: "'none' (default), 'feature_count' or 'attribute_value'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "constraint_field",
                    description: "Field summed per cluster when 'constraint' is 'attribute_value'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_constraint",
                    description: "Minimum per-cluster feature count or attribute total; cuts violating it are skipped.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_constraint",
                    description: "Maximum per-cluster feature count or attribute total; clusters above it are split further even past 'number_of_clusters'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "scale_data",
                    description: "Z-standardise the analysis fields before clustering (default true).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_table",
                    description: "Optional path for the per-cluster summary (means, R-squared per variable, pseudo-F).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        if split_list(require_str(args, "analysis_fields")?).is_empty() {
            return Err(ToolError::Validation(
                "'analysis_fields' must name at least one field".to_string(),
            ));
        }
        let nb = parse_neighborhood(args)?;
        let constraint = parse_constraint(args)?;
        if constraint == Constraint::AttributeValue
            && parse_optional_str(args, "constraint_field")?.is_none()
        {
            return Err(ToolError::Validation(
                "'constraint_field' is required when 'constraint' is 'attribute_value'".to_string(),
            ));
        }
        if let Some(k) = parse_optional_f64(args, "number_of_clusters")? {
            if k < 2.0 {
                return Err(ToolError::Validation(
                    "'number_of_clusters' must be at least 2".to_string(),
                ));
            }
        }
        if nb == Neighborhood::Knn {
            if let Some(k) = parse_optional_f64(args, "number_of_neighbors")? {
                if k < 1.0 {
                    return Err(ToolError::Validation(
                        "'number_of_neighbors' must be at least 1".to_string(),
                    ));
                }
            }
        }
        parse_optional_bool(args, "scale_data")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let field_names = split_list(require_str(args, "analysis_fields")?);
        let k_target = parse_optional_f64(args, "number_of_clusters")?.unwrap_or(5.0) as usize;
        let neighborhood = parse_neighborhood(args)?;
        let knn_k = parse_optional_f64(args, "number_of_neighbors")?.unwrap_or(8.0) as usize;
        let constraint = parse_constraint(args)?;
        let min_c = parse_optional_f64(args, "min_constraint")?;
        let max_c = parse_optional_f64(args, "max_constraint")?;
        let scale = parse_optional_bool(args, "scale_data")?.unwrap_or(true);

        let layer = load_input_layer(input)?;
        let n = layer.features.len();
        if n < 2 {
            return Err(ToolError::Execution(
                "at least 2 features are required".to_string(),
            ));
        }
        let idx: Vec<usize> = field_names
            .iter()
            .map(|f| {
                layer.schema.field_index(f).ok_or_else(|| {
                    ToolError::Validation(format!("analysis field '{f}' not found"))
                })
            })
            .collect::<Result<_, _>>()?;
        let c_idx = match (constraint, parse_optional_str(args, "constraint_field")?) {
            (Constraint::AttributeValue, Some(f)) => Some(
                layer.schema.field_index(f).ok_or_else(|| {
                    ToolError::Validation(format!("constraint_field '{f}' not found"))
                })?,
            ),
            _ => None,
        };

        // Attribute matrix (rows = features, cols = analysis fields).
        let p = idx.len();
        let mut x = vec![0.0_f64; n * p];
        // A substituted 0 is not neutral once the data is z-standardised — it
        // becomes exactly the mean, so a feature with null attributes would be
        // quietly clustered as "average". Substitute, but surface the count.
        let mut unusable_values = 0usize;
        let mut features_with_gaps = 0usize;
        for (i, f) in layer.iter().enumerate() {
            let mut gap = false;
            for (j, &fi) in idx.iter().enumerate() {
                match f.attributes.get(fi).and_then(FieldValue::as_f64) {
                    Some(v) if v.is_finite() => x[i * p + j] = v,
                    _ => {
                        unusable_values += 1;
                        gap = true;
                    }
                }
            }
            if gap {
                features_with_gaps += 1;
            }
        }
        if scale {
            for j in 0..p {
                let col: Vec<f64> = (0..n).map(|i| x[i * p + j]).collect();
                let mean = col.iter().sum::<f64>() / n as f64;
                let var = col.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n as f64;
                let sd = var.sqrt();
                if sd > 0.0 {
                    for i in 0..n {
                        x[i * p + j] = (x[i * p + j] - mean) / sd;
                    }
                }
            }
        }

        // Per-feature constraint weight.
        let weight: Vec<f64> = (0..n)
            .map(|i| match constraint {
                Constraint::None | Constraint::FeatureCount => 1.0,
                Constraint::AttributeValue => layer.features[i]
                    .attributes
                    .get(c_idx.expect("resolved above"))
                    .and_then(FieldValue::as_f64)
                    .unwrap_or(0.0),
            })
            .collect();

        ctx.progress.info("building contiguity graph");
        let adjacency = build_graph(&layer, neighborhood, knn_k)?;
        let edge_count: usize = adjacency.iter().map(|a| a.len()).sum::<usize>() / 2;
        if edge_count == 0 {
            return Err(ToolError::Execution(
                "the contiguity graph has no edges; try 'contiguity_edges_corners' or 'knn'"
                    .to_string(),
            ));
        }

        // Minimum spanning forest (a disconnected input yields one tree per
        // component, which already respects contiguity).
        let mst = minimum_spanning_forest(&adjacency, &x, p, n);
        ctx.progress
            .info(&format!("MST has {} edge(s); cutting", mst.len()));

        let labels = skater(&mst, n, &x, p, &weight, k_target, min_c, max_c, constraint);
        let clusters: Vec<i64> = labels.clone();
        let n_clusters = clusters.iter().collect::<std::collections::BTreeSet<_>>().len();

        // Output layer: input attributes + CLUSTER_ID.
        let mut out = Layer::new("skater_clusters");
        if let Some(gt) = layer.geom_type {
            out = out.with_geom_type(gt);
        }
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for fd in layer.schema.fields() {
            out.add_field(fd.clone());
        }
        out.add_field(FieldDef::new("CLUSTER_ID", FieldType::Integer));
        let names: Vec<String> = layer
            .schema
            .fields()
            .iter()
            .map(|f| f.name.clone())
            .collect();
        for (i, feat) in layer.iter().enumerate() {
            let mut attrs: Vec<(&str, FieldValue)> = names
                .iter()
                .enumerate()
                .map(|(fi, nm)| {
                    (
                        nm.as_str(),
                        feat.attributes.get(fi).cloned().unwrap_or(FieldValue::Null),
                    )
                })
                .collect();
            attrs.push(("CLUSTER_ID", FieldValue::Integer(clusters[i])));
            out.add_feature(feat.geometry.clone(), &attrs)
                .map_err(|e| ToolError::Execution(format!("failed adding feature: {e}")))?;
        }

        // Variance decomposition: how much of each variable's spread the
        // clustering explains.
        let (total_ss, within_ss) = variance_split(&x, p, n, &clusters);
        let r2: Vec<f64> = (0..p)
            .map(|j| {
                if total_ss[j] > 0.0 {
                    1.0 - within_ss[j] / total_ss[j]
                } else {
                    0.0
                }
            })
            .collect();
        let tot: f64 = total_ss.iter().sum();
        let wit: f64 = within_ss.iter().sum();
        let pseudo_f = if n_clusters > 1 && n > n_clusters && wit > 0.0 {
            ((tot - wit) / (n_clusters as f64 - 1.0)) / (wit / (n as f64 - n_clusters as f64))
        } else {
            f64::NAN
        };

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("cluster_count".to_string(), json!(n_clusters));
        outputs.insert("edge_count".to_string(), json!(edge_count));
        outputs.insert("unusable_value_count".to_string(), json!(unusable_values));
        outputs.insert(
            "features_with_missing_values".to_string(),
            json!(features_with_gaps),
        );
        outputs.insert("mst_edge_count".to_string(), json!(mst.len()));
        outputs.insert("neighborhood".to_string(), json!(neighborhood.name()));
        outputs.insert("r_squared".to_string(), json!(r2));
        if pseudo_f.is_finite() {
            outputs.insert("pseudo_f".to_string(), json!(pseudo_f));
        }

        if matches!(args.get("output_table"), Some(v) if !v.is_null()) {
            let mut table = Layer::new("skater_summary");
            table.add_field(FieldDef::new("CLUSTER_ID", FieldType::Integer));
            table.add_field(FieldDef::new("FEATURE_COUNT", FieldType::Integer));
            table.add_field(FieldDef::new("CONSTRAINT_TOTAL", FieldType::Float));
            for nm in &field_names {
                table.add_field(FieldDef::new(format!("MEAN_{nm}"), FieldType::Float));
            }
            let mut per: BTreeMap<i64, (usize, f64, Vec<f64>)> = BTreeMap::new();
            for i in 0..n {
                let e = per
                    .entry(clusters[i])
                    .or_insert_with(|| (0, 0.0, vec![0.0; p]));
                e.0 += 1;
                e.1 += weight[i];
                for j in 0..p {
                    e.2[j] += x[i * p + j];
                }
            }
            for (cid, (count, wsum, sums)) in per {
                let mut attrs = vec![
                    ("CLUSTER_ID", FieldValue::Integer(cid)),
                    ("FEATURE_COUNT", FieldValue::Integer(count as i64)),
                    ("CONSTRAINT_TOTAL", FieldValue::Float(wsum)),
                ];
                let means: Vec<String> = field_names.iter().map(|nm| format!("MEAN_{nm}")).collect();
                for (j, nm) in means.iter().enumerate() {
                    attrs.push((nm.as_str(), FieldValue::Float(sums[j] / count as f64)));
                }
                table
                    .add_feature(None, &attrs)
                    .map_err(|e| ToolError::Execution(format!("failed adding summary row: {e}")))?;
            }
            let p = parse_optional_str(args, "output_table")?;
            outputs.insert("output_table".to_string(), json!(write_or_store_layer(table, p)?));
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── Graph construction ──────────────────────────────────────────────────────

/// Adjacency lists over the chosen neighbourhood definition.
fn build_graph(
    layer: &Layer,
    nb: Neighborhood,
    knn_k: usize,
) -> Result<Vec<Vec<usize>>, ToolError> {
    let n = layer.features.len();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];

    match nb {
        Neighborhood::Knn => {
            let cents: Vec<(f64, f64)> = layer
                .iter()
                .map(|f| f.geometry.as_ref().and_then(centroid).unwrap_or((0.0, 0.0)))
                .collect();
            let mut tree: KdTree<f64, usize, [f64; 2]> = KdTree::new(2);
            for (i, c) in cents.iter().enumerate() {
                tree.add([c.0, c.1], i).ok();
            }
            for (i, c) in cents.iter().enumerate() {
                let found = tree
                    .nearest(&[c.0, c.1], (knn_k + 1).min(n), &squared_euclidean)
                    .map_err(|e| ToolError::Execution(format!("knn search failed: {e}")))?;
                for (_, j) in found {
                    let j = *j;
                    if j != i {
                        push_unique(&mut adj, i, j);
                    }
                }
            }
        }
        // Rook/queen contiguity via shared boundary vertices: two features are
        // queen-neighbours when they share >= 1 vertex, rook-neighbours when
        // they share >= 2 (i.e. an edge rather than only a corner).
        _ => {
            let need = if nb == Neighborhood::ContiguityEdges { 2 } else { 1 };
            let mut by_vertex: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
            for (i, f) in layer.iter().enumerate() {
                let Some(g) = &f.geometry else { continue };
                let mut seen: std::collections::HashSet<(i64, i64)> =
                    std::collections::HashSet::new();
                for c in g.all_coords() {
                    // A HashSet, not a Vec: densified geometry has thousands of
                    // vertices per feature, and `contains` on a Vec makes graph
                    // construction quadratic.
                    let Some(key) = snap_key(c.x, c.y) else { continue };
                    if seen.insert(key) {
                        by_vertex.entry(key).or_default().push(i);
                    }
                }
            }
            let mut shared: HashMap<(usize, usize), usize> = HashMap::new();
            for owners in by_vertex.values() {
                for (a, &i) in owners.iter().enumerate() {
                    for &j in &owners[a + 1..] {
                        let key = if i < j { (i, j) } else { (j, i) };
                        *shared.entry(key).or_insert(0) += 1;
                    }
                }
            }
            for ((i, j), count) in shared {
                if count >= need {
                    push_unique(&mut adj, i, j);
                    push_unique(&mut adj, j, i);
                }
            }
        }
    }
    Ok(adj)
}

/// Quantised vertex key, or `None` when the coordinate is non-finite or so
/// large that `as i64` would saturate — at saturation unrelated features would
/// appear to share a vertex.
fn snap_key(x: f64, y: f64) -> Option<(i64, i64)> {
    let (qx, qy) = (x / SNAP, y / SNAP);
    if !qx.is_finite() || !qy.is_finite() {
        return None;
    }
    const LIMIT: f64 = 9.0e18;
    if qx.abs() >= LIMIT || qy.abs() >= LIMIT {
        return None;
    }
    Some((qx.round() as i64, qy.round() as i64))
}

fn push_unique(adj: &mut [Vec<usize>], i: usize, j: usize) {
    if !adj[i].contains(&j) {
        adj[i].push(j);
    }
    if !adj[j].contains(&i) {
        adj[j].push(i);
    }
}

fn centroid(g: &Geometry) -> Option<(f64, f64)> {
    let coords = g.all_coords();
    if coords.is_empty() {
        return None;
    }
    let n = coords.len() as f64;
    Some((
        coords.iter().map(|c| c.x).sum::<f64>() / n,
        coords.iter().map(|c| c.y).sum::<f64>() / n,
    ))
}

/// Attribute-space distance between two features.
fn attr_dist(x: &[f64], p: usize, i: usize, j: usize) -> f64 {
    (0..p)
        .map(|k| (x[i * p + k] - x[j * p + k]).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Prim's algorithm per connected component, returning the retained edges.
///
/// The frontier is a binary heap with a lazy `in_tree` check on pop. A linearly
/// scanned Vec that is also re-filtered after every accepted edge is O(E^2),
/// which the default 8-neighbour graph reaches on a few thousand features.
fn minimum_spanning_forest(
    adj: &[Vec<usize>],
    x: &[f64],
    p: usize,
    n: usize,
) -> Vec<(usize, usize)> {
    let mut in_tree = vec![false; n];
    let mut edges = Vec::new();
    for seed in 0..n {
        if in_tree[seed] {
            continue;
        }
        in_tree[seed] = true;
        let mut heap: std::collections::BinaryHeap<Candidate> = adj[seed]
            .iter()
            .map(|&j| Candidate {
                dist: attr_dist(x, p, seed, j),
                from: seed,
                to: j,
            })
            .collect();
        while let Some(Candidate { from, to, .. }) = heap.pop() {
            if in_tree[to] {
                continue;
            }
            in_tree[to] = true;
            edges.push((from, to));
            for &j in &adj[to] {
                if !in_tree[j] {
                    heap.push(Candidate {
                        dist: attr_dist(x, p, to, j),
                        from: to,
                        to: j,
                    });
                }
            }
        }
    }
    edges
}

/// Min-heap entry for Prim (BinaryHeap is a max-heap, so ordering is reversed).
struct Candidate {
    dist: f64,
    from: usize,
    to: usize,
}
impl PartialEq for Candidate {
    fn eq(&self, other: &Self) -> bool {
        self.dist == other.dist && self.to == other.to
    }
}
impl Eq for Candidate {}
impl Ord for Candidate {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other
            .dist
            .total_cmp(&self.dist)
            .then_with(|| other.to.cmp(&self.to))
    }
}
impl PartialOrd for Candidate {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

// ── SKATER ──────────────────────────────────────────────────────────────────

/// Recursively cuts the MST, each time removing the edge whose removal buys the
/// largest reduction in total within-cluster sum of squared deviations.
#[allow(clippy::too_many_arguments)]
fn skater(
    mst: &[(usize, usize)],
    n: usize,
    x: &[f64],
    p: usize,
    weight: &[f64],
    k_target: usize,
    min_c: Option<f64>,
    max_c: Option<f64>,
    constraint: Constraint,
) -> Vec<i64> {
    // Start from the MST's connected components (usually one).
    let mut active: Vec<Vec<(usize, usize)>> = Vec::new();
    let mut members: Vec<Vec<usize>> = Vec::new();
    {
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for &(a, b) in mst {
            adj[a].push(b);
            adj[b].push(a);
        }
        let mut seen = vec![false; n];
        for s in 0..n {
            if seen[s] {
                continue;
            }
            let mut comp = Vec::new();
            let mut stack = vec![s];
            seen[s] = true;
            while let Some(u) = stack.pop() {
                comp.push(u);
                for &v in &adj[u] {
                    if !seen[v] {
                        seen[v] = true;
                        stack.push(v);
                    }
                }
            }
            let set: std::collections::HashSet<usize> = comp.iter().copied().collect();
            let comp_edges: Vec<(usize, usize)> = mst
                .iter()
                .filter(|(a, b)| set.contains(a) && set.contains(b))
                .copied()
                .collect();
            members.push(comp);
            active.push(comp_edges);
        }
    }

    // Keep cutting while below the target, or while any cluster busts max.
    loop {
        let over_max = max_c.is_some_and(|m| {
            members
                .iter()
                .any(|c| cluster_total(c, weight, constraint) > m && c.len() > 1)
        });
        if members.len() >= k_target && !over_max {
            break;
        }
        // Pick the (cluster, edge) whose removal reduces within-SS the most.
        let mut best: Option<(usize, usize, f64)> = None; // cluster, edge index, gain
        for (ci, edges) in active.iter().enumerate() {
            if members[ci].len() < 2 {
                continue;
            }
            let base = sum_sq(&members[ci], x, p);
            // Built once per cluster, then reused for every candidate edge:
            // rebuilding it inside the edge loop made the cut search allocate a
            // fresh adjacency map per candidate, per cluster, per cut.
            let adj = cluster_adjacency(edges);
            for (ei, _) in edges.iter().enumerate() {
                let (left, right) = split_components(&members[ci], edges, ei, &adj);
                if left.is_empty() || right.is_empty() {
                    continue;
                }
                if let Some(m) = min_c {
                    if cluster_total(&left, weight, constraint) < m
                        || cluster_total(&right, weight, constraint) < m
                    {
                        continue;
                    }
                }
                let gain = base - sum_sq(&left, x, p) - sum_sq(&right, x, p);
                if best.is_none_or(|(_, _, bg)| gain > bg) {
                    best = Some((ci, ei, gain));
                }
            }
        }
        let Some((ci, ei, _)) = best else { break };
        let edges = active[ci].clone();
        let adj = cluster_adjacency(&edges);
        let (left, right) = split_components(&members[ci], &edges, ei, &adj);
        let set_l: std::collections::HashSet<usize> = left.iter().copied().collect();
        let el: Vec<(usize, usize)> = edges
            .iter()
            .enumerate()
            .filter(|&(i, (a, b))| i != ei && set_l.contains(a) && set_l.contains(b))
            .map(|(_, e)| *e)
            .collect();
        let er: Vec<(usize, usize)> = edges
            .iter()
            .enumerate()
            .filter(|&(i, (a, b))| i != ei && !set_l.contains(a) && !set_l.contains(b))
            .map(|(_, e)| *e)
            .collect();
        members[ci] = left;
        active[ci] = el;
        members.push(right);
        active.push(er);
    }

    let mut labels = vec![0_i64; n];
    for (cid, comp) in members.iter().enumerate() {
        for &i in comp {
            labels[i] = cid as i64;
        }
    }
    labels
}

/// Adjacency for one cluster's spanning tree, as `node -> [(neighbour, edge)]`.
///
/// Built once per cluster and shared across every candidate cut, which skips an
/// edge by index rather than by rebuilding the map without it.
fn cluster_adjacency(edges: &[(usize, usize)]) -> HashMap<usize, Vec<(usize, usize)>> {
    let mut adj: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
    for (i, &(a, b)) in edges.iter().enumerate() {
        adj.entry(a).or_default().push((b, i));
        adj.entry(b).or_default().push((a, i));
    }
    adj
}

/// Splits a cluster's members by removing edge `skip` from its spanning tree.
fn split_components(
    members: &[usize],
    edges: &[(usize, usize)],
    skip: usize,
    adj: &HashMap<usize, Vec<(usize, usize)>>,
) -> (Vec<usize>, Vec<usize>) {
    let start = edges[skip].0;
    let mut seen: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut stack = vec![start];
    seen.insert(start);
    while let Some(u) = stack.pop() {
        if let Some(ns) = adj.get(&u) {
            for &(v, ei) in ns {
                // Skipping by edge index is what lets the adjacency be shared.
                if ei != skip && seen.insert(v) {
                    stack.push(v);
                }
            }
        }
    }
    let left: Vec<usize> = members.iter().copied().filter(|i| seen.contains(i)).collect();
    let right: Vec<usize> = members.iter().copied().filter(|i| !seen.contains(i)).collect();
    (left, right)
}

fn sum_sq(members: &[usize], x: &[f64], p: usize) -> f64 {
    if members.is_empty() {
        return 0.0;
    }
    let n = members.len() as f64;
    let mut total = 0.0;
    for j in 0..p {
        let mean: f64 = members.iter().map(|&i| x[i * p + j]).sum::<f64>() / n;
        total += members
            .iter()
            .map(|&i| (x[i * p + j] - mean).powi(2))
            .sum::<f64>();
    }
    total
}

fn cluster_total(members: &[usize], weight: &[f64], constraint: Constraint) -> f64 {
    match constraint {
        Constraint::None => f64::INFINITY,
        Constraint::FeatureCount => members.len() as f64,
        Constraint::AttributeValue => members.iter().map(|&i| weight[i]).sum(),
    }
}

/// Per-variable total and within-cluster sums of squares.
fn variance_split(x: &[f64], p: usize, n: usize, labels: &[i64]) -> (Vec<f64>, Vec<f64>) {
    let mut total = vec![0.0; p];
    for j in 0..p {
        let mean: f64 = (0..n).map(|i| x[i * p + j]).sum::<f64>() / n as f64;
        total[j] = (0..n).map(|i| (x[i * p + j] - mean).powi(2)).sum();
    }
    let mut groups: BTreeMap<i64, Vec<usize>> = BTreeMap::new();
    for (i, &l) in labels.iter().enumerate() {
        groups.entry(l).or_default().push(i);
    }
    let mut within = vec![0.0; p];
    for members in groups.values() {
        let m = members.len() as f64;
        for j in 0..p {
            let mean: f64 = members.iter().map(|&i| x[i * p + j]).sum::<f64>() / m;
            within[j] += members
                .iter()
                .map(|&i| (x[i * p + j] - mean).powi(2))
                .sum::<f64>();
        }
    }
    (total, within)
}

// ── Params ──────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Neighborhood {
    ContiguityEdges,
    ContiguityEdgesCorners,
    Knn,
}

impl Neighborhood {
    fn name(self) -> &'static str {
        match self {
            Neighborhood::ContiguityEdges => "contiguity_edges",
            Neighborhood::ContiguityEdgesCorners => "contiguity_edges_corners",
            Neighborhood::Knn => "knn",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Constraint {
    None,
    FeatureCount,
    AttributeValue,
}

fn parse_neighborhood(args: &ToolArgs) -> Result<Neighborhood, ToolError> {
    match args
        .get("neighborhood")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        None | Some("") | Some("contiguity_edges") => Ok(Neighborhood::ContiguityEdges),
        Some("contiguity_edges_corners") => Ok(Neighborhood::ContiguityEdgesCorners),
        Some("knn") => Ok(Neighborhood::Knn),
        Some("trimmed_delaunay") => Err(ToolError::Validation(
            "'trimmed_delaunay' is not implemented in this version; use 'knn' to link \
             features that do not share a boundary"
                .to_string(),
        )),
        Some(o) => Err(ToolError::Validation(format!(
            "'neighborhood' must be 'contiguity_edges', 'contiguity_edges_corners' or 'knn', got '{o}'"
        ))),
    }
}

fn parse_constraint(args: &ToolArgs) -> Result<Constraint, ToolError> {
    match args
        .get("constraint")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        None | Some("") | Some("none") => Ok(Constraint::None),
        Some("feature_count") => Ok(Constraint::FeatureCount),
        Some("attribute_value") => Ok(Constraint::AttributeValue),
        Some(o) => Err(ToolError::Validation(format!(
            "'constraint' must be 'none', 'feature_count' or 'attribute_value', got '{o}'"
        ))),
    }
}

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
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

fn split_list(s: &str) -> Vec<String> {
    s.split([',', ';'])
        .map(str::trim)
        .filter(|p| !p.is_empty())
        .map(str::to_string)
        .collect()
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{Coord, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// A 1 x n strip of unit squares, each carrying value `vals[i]`.
    fn strip(vals: &[f64]) -> String {
        let mut l = Layer::new("s")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (i, v) in vals.iter().enumerate() {
            let x = i as f64;
            l.add_feature(
                Some(Geometry::polygon(
                    vec![
                        Coord::xy(x, 0.0),
                        Coord::xy(x + 1.0, 0.0),
                        Coord::xy(x + 1.0, 1.0),
                        Coord::xy(x, 1.0),
                        Coord::xy(x, 0.0),
                    ],
                    vec![],
                )),
                &[("v", FieldValue::Float(*v))],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SpatiallyConstrainedMultivariateClusteringTool
            .run(&args, &ctx())
            .unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn ids(l: &Layer) -> Vec<i64> {
        let i = l.schema.field_index("CLUSTER_ID").unwrap();
        l.iter().map(|f| f.attributes[i].as_i64().unwrap()).collect()
    }

    #[test]
    fn clusters_are_spatially_contiguous_runs() {
        // Values 1,1,1,9,9,9 along a strip: the natural cut is between index
        // 2 and 3, and every cluster must be a contiguous run.
        let (out, layer) = run(json!({
            "input": strip(&[1.0, 1.0, 1.0, 9.0, 9.0, 9.0]),
            "analysis_fields": "v", "number_of_clusters": 2
        }));
        assert_eq!(out.outputs["cluster_count"], json!(2));
        let c = ids(&layer);
        assert_eq!(c[0], c[1]);
        assert_eq!(c[1], c[2]);
        assert_eq!(c[3], c[4]);
        assert_eq!(c[4], c[5]);
        assert_ne!(c[2], c[3]);
    }

    #[test]
    fn a_high_value_in_the_middle_cannot_join_a_distant_twin() {
        // Values 1,9,1,1,1,9: an unconstrained k-means would put the two 9s
        // together. Contiguity forbids that — the defining property.
        let (_o, layer) = run(json!({
            "input": strip(&[1.0, 9.0, 1.0, 1.0, 1.0, 9.0]),
            "analysis_fields": "v", "number_of_clusters": 3
        }));
        let c = ids(&layer);
        assert_ne!(c[1], c[5], "non-adjacent features shared a cluster");
    }

    #[test]
    fn every_cluster_is_connected_in_the_contiguity_graph() {
        let (_o, layer) = run(json!({
            "input": strip(&[1.0, 5.0, 2.0, 8.0, 3.0, 9.0, 4.0, 7.0]),
            "analysis_fields": "v", "number_of_clusters": 4
        }));
        let c = ids(&layer);
        // On a strip, contiguity means each cluster's indices form one run.
        let mut seen: BTreeMap<i64, (usize, usize, usize)> = BTreeMap::new();
        for (i, &cid) in c.iter().enumerate() {
            let e = seen.entry(cid).or_insert((i, i, 0));
            e.0 = e.0.min(i);
            e.1 = e.1.max(i);
            e.2 += 1;
        }
        for (cid, (lo, hi, count)) in seen {
            assert_eq!(hi - lo + 1, count, "cluster {cid} is not a contiguous run");
        }
    }

    #[test]
    fn r_squared_is_reported_and_bounded() {
        let (out, _l) = run(json!({
            "input": strip(&[1.0, 1.0, 1.0, 9.0, 9.0, 9.0]),
            "analysis_fields": "v", "number_of_clusters": 2
        }));
        let r2 = out.outputs["r_squared"].as_array().unwrap()[0]
            .as_f64()
            .unwrap();
        assert!(r2 > 0.9, "clean two-group split should explain most variance, got {r2}");
        assert!(r2 <= 1.0 + 1e-9);
    }

    #[test]
    fn min_constraint_blocks_tiny_clusters() {
        // Without a constraint the best first cut isolates the lone 9.
        let vals = [1.0, 1.0, 1.0, 1.0, 9.0];
        let (_o, free) = run(json!({
            "input": strip(&vals), "analysis_fields": "v", "number_of_clusters": 2
        }));
        let sizes = |l: &Layer| -> Vec<usize> {
            let c = ids(l);
            let mut m: BTreeMap<i64, usize> = BTreeMap::new();
            for cid in c {
                *m.entry(cid).or_default() += 1;
            }
            let mut v: Vec<usize> = m.into_values().collect();
            v.sort_unstable();
            v
        };
        assert_eq!(sizes(&free), vec![1, 4]);
        let (_o, held) = run(json!({
            "input": strip(&vals), "analysis_fields": "v", "number_of_clusters": 2,
            "constraint": "feature_count", "min_constraint": 2
        }));
        assert!(sizes(&held).iter().all(|&s| s >= 2), "min constraint ignored");
    }

    #[test]
    fn max_constraint_forces_extra_splits() {
        let (out, _l) = run(json!({
            "input": strip(&[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
            "analysis_fields": "v", "number_of_clusters": 2,
            "constraint": "feature_count", "max_constraint": 2
        }));
        // 6 features capped at 2 each needs at least 3 clusters, not 2.
        assert!(out.outputs["cluster_count"].as_f64().unwrap() >= 3.0);
    }

    #[test]
    fn knn_neighborhood_links_disjoint_features() {
        // Two separated squares share no vertices, so rook contiguity finds no
        // edges at all; knn must still connect them.
        let mut l = Layer::new("d")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (i, x) in [0.0, 100.0, 200.0].iter().enumerate() {
            l.add_feature(
                Some(Geometry::polygon(
                    vec![
                        Coord::xy(*x, 0.0),
                        Coord::xy(x + 1.0, 0.0),
                        Coord::xy(x + 1.0, 1.0),
                        Coord::xy(*x, 1.0),
                        Coord::xy(*x, 0.0),
                    ],
                    vec![],
                )),
                &[("v", FieldValue::Float(i as f64))],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        let input = wbvector::memory_store::make_vector_memory_path(&id);

        let rook: ToolArgs = serde_json::from_value(json!({
            "input": input.clone(), "analysis_fields": "v", "number_of_clusters": 2
        }))
        .unwrap();
        assert!(SpatiallyConstrainedMultivariateClusteringTool
            .run(&rook, &ctx())
            .is_err());

        let (out, _l) = run(json!({
            "input": input, "analysis_fields": "v", "number_of_clusters": 2,
            "neighborhood": "knn", "number_of_neighbors": 2
        }));
        assert_eq!(out.outputs["cluster_count"], json!(2));
    }

    #[test]
    fn summary_table_is_emitted_on_request() {
        let args: ToolArgs = serde_json::from_value(json!({
            "input": strip(&[1.0, 1.0, 9.0, 9.0]),
            "analysis_fields": "v", "number_of_clusters": 2, "output_table": ""
        }))
        .unwrap();
        let out = SpatiallyConstrainedMultivariateClusteringTool
            .run(&args, &ctx())
            .unwrap();
        let t = load_input_layer(out.outputs["output_table"].as_str().unwrap()).unwrap();
        assert_eq!(t.features.len(), 2);
        assert!(t.schema.field_index("MEAN_v").is_some());
        assert!(t.schema.field_index("FEATURE_COUNT").is_some());
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SpatiallyConstrainedMultivariateClusteringTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "p.shp" })).is_err());
        assert!(bad(json!({
            "input": "p.shp", "analysis_fields": "v", "number_of_clusters": 1
        }))
        .is_err());
        assert!(bad(json!({
            "input": "p.shp", "analysis_fields": "v", "neighborhood": "trimmed_delaunay"
        }))
        .is_err());
        assert!(bad(json!({
            "input": "p.shp", "analysis_fields": "v", "constraint": "attribute_value"
        }))
        .is_err());
        assert!(bad(json!({ "input": "p.shp", "analysis_fields": "v" })).is_ok());
    }
}
