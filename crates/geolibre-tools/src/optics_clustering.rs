//! GeoLibre tool: OPTICS multi-scale density-based clustering.
//!
//! Pure-Rust counterpart of the OPTICS (Multi-scale) option in ArcGIS Pro's
//! *Density-based Clustering* (Spatial Statistics). The authored suite ships
//! `dbscan` and `hdbscan`; OPTICS is the third member ArcGIS bundles and the one
//! that recovers clusters of **varying density** without a single `eps`: it
//! builds a reachability ordering of the points and extracts clusters from the
//! shape of that plot.
//!
//! The implementation mirrors scikit-learn's OPTICS:
//!
//! 1. **Core distance** of each point — the distance to its `min_features_cluster`-th
//!    nearest neighbour (including itself). Beyond an optional `search_distance`
//!    (max ε) the core distance is undefined and the point cannot expand.
//! 2. **Reachability ordering** — process the point with the smallest current
//!    reachability, relaxing every neighbour to `max(core_dist, point_dist)`.
//! 3. **ξ-steep extraction** — walk the reachability plot, pairing steep-down and
//!    steep-up areas into clusters (Ankerst et al., the same method sklearn's
//!    `cluster_method='xi'` uses). `cluster_sensitivity` (0-100, higher → more
//!    clusters) maps to ξ = 1 − sensitivity/100.
//!
//! Output copies the input points and adds `cluster_id` (−1 = noise) and
//! `reachability` (the point's reachability-plot value, a density diagnostic).
//! Deterministic; O(n²), so it suits moderate point counts. Use a projected CRS.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct OpticsClusteringTool;

impl Tool for OpticsClusteringTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "optics_clustering",
            display_name: "OPTICS Clustering",
            summary: "Multi-scale density-based clustering via OPTICS (like the OPTICS option of ArcGIS Density-based Clustering): builds a reachability ordering and extracts clusters of varying density with the ξ-steep method — the multi-scale member of the DBSCAN/HDBSCAN/OPTICS trio the authored dbscan and hdbscan don't cover.",
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
                    name: "min_features_cluster",
                    description: "Minimum points for the core-distance density estimate and the smallest cluster (>= 2). Default 5.",
                    required: false,
                },
                ToolParamSpec {
                    name: "search_distance",
                    description: "Optional maximum reachability distance (max epsilon). Default: unbounded (all scales).",
                    required: false,
                },
                ToolParamSpec {
                    name: "cluster_sensitivity",
                    description: "0-100; higher detects more/tighter clusters. Maps to xi = 1 - sensitivity/100. Default 95 (xi = 0.05).",
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

        // Representative point per feature; remember the source feature.
        let mut pts: Vec<(f64, f64)> = Vec::new();
        let mut feat_of: Vec<usize> = Vec::new();
        for (fi, feature) in layer.features.iter().enumerate() {
            if let Some((x, y)) = feature.geometry.as_ref().and_then(point_xy) {
                pts.push((x, y));
                feat_of.push(fi);
            }
        }
        let n = pts.len();
        if n < prm.min_samples {
            return Err(ToolError::Execution(format!(
                "need at least min_features_cluster ({}) point features, found {n}",
                prm.min_samples
            )));
        }

        ctx.progress
            .info(&format!("OPTICS ordering of {n} point(s)"));
        let graph = optics_graph(&pts, prm.min_samples, prm.max_eps);
        let reach_plot: Vec<f64> = graph
            .ordering
            .iter()
            .map(|&i| graph.reachability[i])
            .collect();
        let pred_plot: Vec<i64> = graph
            .ordering
            .iter()
            .map(|&i| graph.predecessor[i])
            .collect();

        ctx.progress
            .info(&format!("xi-extraction (xi = {:.3})", prm.xi));
        let clusters = xi_cluster(
            &reach_plot,
            &pred_plot,
            &graph.ordering,
            prm.xi,
            prm.min_samples,
            prm.min_samples,
        );
        let labels_by_sample = extract_labels(&graph.ordering, &clusters, n);

        // Write cluster fields onto a copy of the input.
        layer.add_field(FieldDef::new("cluster_id", FieldType::Integer));
        layer.add_field(FieldDef::new("reachability", FieldType::Float));
        // reachability per sample (undefined/first point -> its stored value).
        let mut per_feat: Vec<(i64, f64)> = vec![(-1, -1.0); layer.len()];
        for i in 0..n {
            let r = graph.reachability[i];
            let r_out = if r.is_finite() { r } else { -1.0 };
            per_feat[feat_of[i]] = (labels_by_sample[i], r_out);
        }
        for (fi, feature) in layer.features.iter_mut().enumerate() {
            let (c, r) = per_feat[fi];
            feature.attributes.push(FieldValue::Integer(c));
            feature.attributes.push(FieldValue::Float(r));
        }

        let n_clusters = labels_by_sample
            .iter()
            .filter(|&&l| l >= 0)
            .copied()
            .max()
            .map(|m| (m + 1) as usize)
            .unwrap_or(0);
        let n_noise = labels_by_sample.iter().filter(|&&l| l < 0).count();
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

// ── OPTICS reachability graph (mirrors sklearn compute_optics_graph) ──────────

struct Graph {
    ordering: Vec<usize>,
    reachability: Vec<f64>,
    predecessor: Vec<i64>,
}

/// Rounds to 15 decimal places, matching sklearn's `np.around(x, precision)` on
/// core/reachability distances (float64 precision = 15).
fn round15(x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    (x * 1e15).round() / 1e15
}

fn optics_graph(pts: &[(f64, f64)], min_samples: usize, max_eps: f64) -> Graph {
    let n = pts.len();
    let dist = |a: usize, b: usize| -> f64 {
        let (dx, dy) = (pts[a].0 - pts[b].0, pts[a].1 - pts[b].1);
        dx.hypot(dy)
    };

    // Core distance: the min_samples-th smallest distance (self included at 0),
    // i.e. index min_samples-1 of the ascending distance list. Beyond max_eps it
    // is undefined (inf).
    let mut core = vec![f64::INFINITY; n];
    for (i, core_i) in core.iter_mut().enumerate() {
        let mut ds: Vec<f64> = (0..n).map(|j| dist(i, j)).collect();
        ds.sort_by(f64::total_cmp);
        let k = (min_samples - 1).min(n - 1);
        let cd = round15(ds[k]);
        *core_i = if cd > max_eps { f64::INFINITY } else { cd };
    }

    let mut reachability = vec![f64::INFINITY; n];
    let mut predecessor = vec![-1i64; n];
    let mut processed = vec![false; n];
    let mut ordering = vec![0usize; n];

    for slot in ordering.iter_mut() {
        // Choose the unprocessed point with the smallest reachability; ties go to
        // the smallest index (np.argmin over ascending `where(~processed)`).
        let mut point = usize::MAX;
        let mut best = f64::INFINITY;
        for j in 0..n {
            if processed[j] && point != usize::MAX {
                continue;
            }
            if !processed[j] {
                let r = reachability[j];
                if point == usize::MAX || r < best {
                    best = r;
                    point = j;
                }
            }
        }
        processed[point] = true;
        *slot = point;

        if core[point].is_finite() {
            // Relax every unprocessed neighbour within max_eps.
            for j in 0..n {
                if processed[j] {
                    continue;
                }
                let d = dist(point, j);
                if d > max_eps {
                    continue;
                }
                let rdist = round15(d.max(core[point]));
                if rdist < reachability[j] {
                    reachability[j] = rdist;
                    predecessor[j] = point as i64;
                }
            }
        }
    }

    Graph {
        ordering,
        reachability,
        predecessor,
    }
}

// ── ξ-steep cluster extraction (mirrors sklearn _xi_cluster) ──────────────────

#[derive(Clone)]
struct Sda {
    start: usize,
    mib: f64,
}

/// Extend a steep region until maximal. `steep_point`/`xward_point` are the
/// steep and monotone masks for the direction being extended.
fn extend_region(
    steep_point: &[bool],
    xward_point: &[bool],
    start: usize,
    min_samples: usize,
) -> usize {
    let n = steep_point.len();
    let mut non_xward = 0usize;
    let mut index = start;
    let mut end = start;
    while index < n {
        if steep_point[index] {
            non_xward = 0;
            end = index;
        } else if !xward_point[index] {
            non_xward += 1;
            if non_xward > min_samples {
                break;
            }
        } else {
            return end;
        }
        index += 1;
    }
    end
}

/// Drop steep-down areas no longer reachable given the new maximum-in-between.
fn update_filter_sdas(sdas: &[Sda], mib: f64, xi_complement: f64, rplot: &[f64]) -> Vec<Sda> {
    if mib.is_infinite() {
        return Vec::new();
    }
    sdas.iter()
        .filter(|sda| mib <= rplot[sda.start] * xi_complement)
        .map(|sda| Sda {
            start: sda.start,
            mib: sda.mib.max(mib),
        })
        .collect()
}

/// Predecessor correction (Schubert & Gertz 2018, Algorithm 2).
fn correct_predecessor(
    rplot: &[f64],
    pred_plot: &[i64],
    ordering: &[usize],
    s: usize,
    mut e: usize,
) -> Option<(usize, usize)> {
    while s < e {
        if rplot[s] > rplot[e] {
            return Some((s, e));
        }
        let p_e = pred_plot[e];
        if ordering[s..e].iter().any(|&o| o as i64 == p_e) {
            return Some((s, e));
        }
        e -= 1;
    }
    None
}

/// The Xi-steep method. Returns clusters as inclusive `(start, end)` ranges in
/// ordering-position space, smaller clusters before the larger ones enclosing
/// them.
fn xi_cluster(
    reach_plot: &[f64],
    pred_plot: &[i64],
    ordering: &[usize],
    xi: f64,
    min_samples: usize,
    min_cluster_size: usize,
) -> Vec<(usize, usize)> {
    let n = reach_plot.len();
    // Append inf so clusters can end at the plot's end.
    let mut rplot = Vec::with_capacity(n + 1);
    rplot.extend_from_slice(reach_plot);
    rplot.push(f64::INFINITY);

    let xi_complement = 1.0 - xi;

    // ratio[i] = rplot[i] / rplot[i+1], length n.
    let mut steep_upward = vec![false; n];
    let mut steep_downward = vec![false; n];
    let mut upward = vec![false; n];
    let mut downward = vec![false; n];
    for i in 0..n {
        let ratio = rplot[i] / rplot[i + 1];
        // NaN (inf/inf, 0/0) yields all-false, matching numpy's errstate ignore.
        steep_upward[i] = ratio <= xi_complement;
        steep_downward[i] = ratio >= 1.0 / xi_complement;
        downward[i] = ratio > 1.0;
        upward[i] = ratio < 1.0;
    }

    let mut sdas: Vec<Sda> = Vec::new();
    let mut clusters: Vec<(usize, usize)> = Vec::new();
    let mut index = 0usize;
    let mut mib = 0.0f64;

    let steep_indices: Vec<usize> = (0..n)
        .filter(|&i| steep_upward[i] || steep_downward[i])
        .collect();

    for &steep_index in &steep_indices {
        if steep_index < index {
            continue;
        }
        // mib = max(mib, max(rplot[index ..= steep_index]))
        for &v in &rplot[index..=steep_index] {
            if v > mib {
                mib = v;
            }
        }

        if steep_downward[steep_index] {
            sdas = update_filter_sdas(&sdas, mib, xi_complement, &rplot);
            let d_start = steep_index;
            let d_end = extend_region(&steep_downward, &upward, d_start, min_samples);
            sdas.push(Sda {
                start: d_start,
                mib: 0.0,
            });
            index = d_end + 1;
            mib = rplot[index];
        } else {
            sdas = update_filter_sdas(&sdas, mib, xi_complement, &rplot);
            let u_start = steep_index;
            let u_end = extend_region(&steep_upward, &downward, u_start, min_samples);
            index = u_end + 1;
            mib = rplot[index];

            let mut u_clusters: Vec<(usize, usize)> = Vec::new();
            for d in &sdas {
                let mut c_start = d.start;
                let mut c_end = u_end;

                // line (**), sc2*
                if rplot[c_end + 1] * xi_complement < d.mib {
                    continue;
                }
                // Definition 11: criterion 4
                let d_max = rplot[d.start];
                let d_end_region = extend_region(&steep_downward, &upward, d.start, min_samples);
                if d_max * xi_complement >= rplot[c_end + 1] {
                    while rplot[c_start + 1] > rplot[c_end + 1] && c_start < d_end_region {
                        c_start += 1;
                    }
                } else if rplot[c_end + 1] * xi_complement >= d_max {
                    while c_end > u_start && rplot[c_end - 1] > d_max {
                        c_end -= 1;
                    }
                }

                // predecessor correction
                match correct_predecessor(&rplot, pred_plot, ordering, c_start, c_end) {
                    Some((s, e)) => {
                        c_start = s;
                        c_end = e;
                    }
                    None => continue,
                }

                // criterion 3.a
                if c_end + 1 - c_start < min_cluster_size {
                    continue;
                }
                // criterion 1
                if c_start > d_end_region {
                    continue;
                }
                // criterion 2
                if c_end < u_start {
                    continue;
                }
                u_clusters.push((c_start, c_end));
            }
            u_clusters.reverse();
            clusters.extend(u_clusters);
        }
    }

    clusters
}

/// Assign labels from the extracted clusters (smaller clusters first win their
/// points), then map back to sample-id order.
fn extract_labels(ordering: &[usize], clusters: &[(usize, usize)], n: usize) -> Vec<i64> {
    let mut labels_ord = vec![-1i64; n];
    let mut label = 0i64;
    for &(s, e) in clusters {
        if (s..=e).all(|i| labels_ord[i] == -1) {
            for l in labels_ord.iter_mut().take(e + 1).skip(s) {
                *l = label;
            }
            label += 1;
        }
    }
    // labels[ordering] = labels: sample ordering[pos] gets the position's label.
    let mut by_sample = vec![-1i64; n];
    for (pos, &sample) in ordering.iter().enumerate() {
        by_sample[sample] = labels_ord[pos];
    }
    by_sample
}

// ── Geometry / parameters ─────────────────────────────────────────────────────

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
    min_samples: usize,
    max_eps: f64,
    xi: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let min_samples = match parse_opt_usize(args, "min_features_cluster")? {
        None => 5,
        Some(v) if v >= 2 => v,
        Some(_) => {
            return Err(ToolError::Validation(
                "'min_features_cluster' must be >= 2".to_string(),
            ))
        }
    };
    let max_eps = match parse_opt_f64(args, "search_distance")? {
        None => f64::INFINITY,
        Some(v) if v > 0.0 => v,
        Some(_) => {
            return Err(ToolError::Validation(
                "'search_distance' must be > 0".to_string(),
            ))
        }
    };
    let sensitivity = match parse_opt_f64(args, "cluster_sensitivity")? {
        None => 95.0,
        Some(v) if (0.0..=100.0).contains(&v) => v,
        Some(_) => {
            return Err(ToolError::Validation(
                "'cluster_sensitivity' must be in [0, 100]".to_string(),
            ))
        }
    };
    let xi = (1.0 - sensitivity / 100.0).clamp(0.001, 0.999);
    Ok(Params {
        min_samples,
        max_eps,
        xi,
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

fn parse_opt_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(Some(n.as_f64().unwrap_or(f64::NAN))),
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
        let out = OpticsClusteringTool.run(&args, &ctx()).unwrap();
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

    /// A grid blob and a second grid blob far away become two clusters.
    #[test]
    fn separates_two_blobs() {
        let mut pts = Vec::new();
        for i in 0..16 {
            pts.push((0.1 * (i % 4) as f64, 0.1 * (i / 4) as f64));
        }
        for i in 0..16 {
            pts.push((50.0 + 0.1 * (i % 4) as f64, 50.0 + 0.1 * (i / 4) as f64));
        }
        let (out, layer) = run(json!({ "input": layer_of(&pts), "min_features_cluster": 4 }));
        assert_eq!(out.outputs["cluster_count"], json!(2));
        let lab = labels(&layer);
        assert!(lab[0..16].iter().all(|&l| l == lab[0] && l >= 0));
        assert!(lab[16..32].iter().all(|&l| l == lab[16] && l >= 0));
        assert_ne!(lab[0], lab[16]);
    }

    /// Multi-scale: a dense blob embedded near a sparse blob — OPTICS recovers
    /// both, which a single-eps DBSCAN cannot.
    #[test]
    fn recovers_varying_density() {
        let mut pts = Vec::new();
        // Dense blob: 20 pts in a tight 0.5-wide cluster.
        for i in 0..20 {
            let a = i as f64 * 0.31;
            pts.push((10.0 + 0.25 * a.cos(), 10.0 + 0.25 * a.sin()));
        }
        // Sparse blob: 20 pts spread over a 6-wide cluster, far away.
        for i in 0..20 {
            let a = i as f64 * 0.31;
            pts.push((60.0 + 3.0 * a.cos(), 60.0 + 3.0 * a.sin()));
        }
        let (out, _l) = run(json!({ "input": layer_of(&pts), "min_features_cluster": 5 }));
        assert!(
            out.outputs["cluster_count"].as_u64().unwrap() >= 2,
            "should recover both the dense and sparse clusters"
        );
    }

    /// Scattered points between two dense blobs are labelled noise. This exact
    /// configuration yields 2 clusters and 1 noise point (index 21) under
    /// scikit-learn OPTICS(min_samples=4, xi=0.05); the tool reproduces it.
    #[test]
    fn isolates_noise() {
        let mut pts = Vec::new();
        for i in 0..10 {
            pts.push((0.15 * (i % 5) as f64, 0.15 * (i / 5) as f64)); // blob A
        }
        for i in 0..10 {
            pts.push((30.0 + 0.15 * (i % 5) as f64, 0.15 * (i / 5) as f64)); // blob B
        }
        pts.push((8.0, 10.0)); // 20
        pts.push((14.0, -9.0)); // 21 -> noise
        pts.push((20.0, 12.0)); // 22
        let (out, layer) = run(json!({ "input": layer_of(&pts), "min_features_cluster": 4 }));
        assert_eq!(out.outputs["cluster_count"], json!(2));
        let lab = labels(&layer);
        assert_eq!(lab[21], -1, "the isolated scatter point is noise");
        assert!(out.outputs["noise_count"].as_u64().unwrap() >= 1);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            OpticsClusteringTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "min_features_cluster": 1 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "cluster_sensitivity": 150 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "min_features_cluster": 5 })).is_ok());
    }
}
