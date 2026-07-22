//! GeoLibre tool: attribute-space (multivariate) clustering of vector features.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Multivariate Clustering* (Spatial
//! Statistics). The bundled `k_means_clustering` clusters raster bands only;
//! `dbscan`/`hdbscan` group points by spatial density; `build_balanced_zones`
//! partitions contiguous polygons. None of them group arbitrary vector
//! features purely by their *attribute profile*, ignoring location — e.g.
//! segmenting census tracts by demographics regardless of where they sit.
//!
//! Pipeline:
//!
//! 1. Read the analysis `fields` for every feature and z-score standardize
//!    each field across all features (mean 0, std 1) so fields on different
//!    scales contribute equally.
//! 2. Cluster the standardized rows with either:
//!    - `kmeans` (default) — seeded k-means++ initialization followed by
//!      Lloyd iterations, with several seeded restarts; the lowest-WCSS
//!      result wins.
//!    - `kmedoids` — deterministic k-medoids (PAM) on the full pairwise
//!      distance matrix, suited to smaller feature counts.
//! 3. If `num_clusters` is omitted or `0`, evaluate every k in `2..=15`
//!    (capped at `n - 1`) and keep the k with the highest
//!    Calinski-Harabasz pseudo-F (the classic ArcGIS "evaluate optimal
//!    number of clusters" heuristic).
//!
//! All randomness is a tiny inline seeded splitmix64 generator — no
//! `Date::now`, no unseeded RNG — so the same `seed` reproduces identical
//! cluster assignments run to run, which matters for WASM determinism.
//!
//! Output copies the input features and adds `cluster_id` (0-based; `-1` for
//! features missing one or more analysis field values) and `silhouette` (a
//! per-feature separation score in `[-1, 1]`). The result also reports
//! per-cluster field means and, in optimal-k mode, the full pseudo-F curve.

use std::collections::BTreeMap;

use serde_json::{json, Map, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

const MIN_K: usize = 2;
const MAX_K: usize = 15;
const KMEANS_MAX_ITERS: usize = 100;
const KMEANS_RESTARTS: usize = 8;
const PAM_MAX_ITERS: usize = 50;
const PAM_RESTARTS: usize = 8;

#[derive(Clone, Copy, PartialEq)]
enum Method {
    Kmeans,
    Kmedoids,
}

impl Method {
    fn as_str(self) -> &'static str {
        match self {
            Method::Kmeans => "kmeans",
            Method::Kmedoids => "kmedoids",
        }
    }
}

pub struct MultivariateClusteringTool;

impl Tool for MultivariateClusteringTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "multivariate_clustering",
            display_name: "Multivariate Clustering",
            summary: "Cluster vector features by the similarity of their standardized attribute values (seeded k-means or k-medoids), automatically choosing the number of clusters by Calinski-Harabasz pseudo-F when not given, like ArcGIS Multivariate Clustering.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer, format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "fields",
                    description: "Comma-separated numeric analysis fields (attribute space to cluster on).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (a copy of the input with cluster_id/silhouette). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "num_clusters",
                    description: "Number of clusters k. Omit or pass 0 to evaluate k = 2..15 automatically by Calinski-Harabasz pseudo-F.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "Clustering method: 'kmeans' (default) or 'kmedoids' (PAM; suited to smaller feature counts).",
                    required: false,
                },
                ToolParamSpec {
                    name: "seed",
                    description: "Random seed for deterministic k-means++ / medoid initialization (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "fields")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let schema = layer.schema.clone();

        let field_indices: Vec<usize> = prm
            .fields
            .iter()
            .map(|f| {
                schema
                    .field_index(f)
                    .ok_or_else(|| ToolError::Validation(format!("field '{f}' not found in input")))
            })
            .collect::<Result<_, _>>()?;
        let nf = field_indices.len();

        // Rows with a value for every analysis field are clustered; the rest
        // get cluster_id = -1 (like hdbscan's noise convention).
        let mut attrs: Vec<Vec<f64>> = Vec::new();
        let mut included: Vec<usize> = Vec::new();
        for (fi, feature) in layer.features.iter().enumerate() {
            let row: Option<Vec<f64>> = field_indices
                .iter()
                .map(|&idx| feature.attributes.get(idx).and_then(FieldValue::as_f64))
                .collect();
            if let Some(row) = row {
                if row.iter().all(|v| v.is_finite()) {
                    attrs.push(row);
                    included.push(fi);
                }
            }
        }
        let n = attrs.len();
        if n < MIN_K + 1 {
            return Err(ToolError::Execution(format!(
                "need at least {} features with valid values for all analysis fields, found {n}",
                MIN_K + 1
            )));
        }

        let z = standardize(&attrs, nf);

        let auto = prm.num_clusters.is_none();
        // Cap k so every candidate cluster can average at least 2 members;
        // otherwise Calinski-Harabasz over-fits toward near-singleton
        // clusters on small samples.
        let k_max = MAX_K.min(n - 1).min((n / 2).max(MIN_K));
        let (k, assign, pseudo_f, k_curve) = if auto {
            ctx.progress.info(&format!(
                "evaluating k = {MIN_K}..={k_max} by Calinski-Harabasz pseudo-F"
            ));
            let mut best: Option<(usize, Vec<usize>, f64)> = None;
            let mut curve = Vec::new();
            for k in MIN_K..=k_max {
                let assign = cluster(&z, k, prm.method, prm.seed);
                let ch = calinski_harabasz(&z, &assign, k);
                curve.push(json!({ "k": k, "pseudo_f": ch }));
                if best.as_ref().is_none_or(|(_, _, bch)| ch > *bch) {
                    best = Some((k, assign, ch));
                }
            }
            let (k, assign, ch) = best.expect("MIN_K..=k_max is non-empty");
            (k, assign, ch, Some(curve))
        } else {
            let k = prm.num_clusters.unwrap();
            if !(MIN_K..=n - 1).contains(&k) {
                return Err(ToolError::Validation(format!(
                    "'num_clusters' must be between {MIN_K} and {} ({} feature(s) with valid data), got {k}",
                    n - 1,
                    n
                )));
            }
            ctx.progress
                .info(&format!("clustering {n} feature(s) into {k} cluster(s)"));
            let assign = cluster(&z, k, prm.method, prm.seed);
            let ch = calinski_harabasz(&z, &assign, k);
            (k, assign, ch, None)
        };

        let sil = silhouette_scores(&z, &assign, k);

        // Per-cluster means in original (non-standardized) units.
        let mut counts = vec![0usize; k];
        let mut sums = vec![vec![0.0f64; nf]; k];
        for (i, &c) in assign.iter().enumerate() {
            counts[c] += 1;
            for d in 0..nf {
                sums[c][d] += attrs[i][d];
            }
        }
        let cluster_means: Vec<Value> = (0..k)
            .map(|c| {
                let mut m = Map::new();
                m.insert("cluster_id".to_string(), json!(c));
                m.insert("count".to_string(), json!(counts[c]));
                for (d, name) in prm.fields.iter().enumerate() {
                    let mean = if counts[c] > 0 {
                        sums[c][d] / counts[c] as f64
                    } else {
                        0.0
                    };
                    m.insert(name.clone(), json!(mean));
                }
                Value::Object(m)
            })
            .collect();

        // Write cluster_id / silhouette onto a copy of the input.
        layer.add_field(FieldDef::new("cluster_id", FieldType::Integer));
        layer.add_field(FieldDef::new("silhouette", FieldType::Float));
        let mut per_feature: Vec<(i64, f64)> = vec![(-1, 0.0); layer.features.len()];
        for (row_i, &fi) in included.iter().enumerate() {
            per_feature[fi] = (assign[row_i] as i64, sil[row_i]);
        }
        for (fi, feature) in layer.features.iter_mut().enumerate() {
            let (c, s) = per_feature[fi];
            feature.attributes.push(FieldValue::Integer(c));
            feature.attributes.push(FieldValue::Float(s));
        }

        let feature_count = layer.len();
        let unclustered_count = feature_count - n;
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("unclustered_count".to_string(), json!(unclustered_count));
        outputs.insert("num_clusters".to_string(), json!(k));
        outputs.insert("method".to_string(), json!(prm.method.as_str()));
        outputs.insert("seed".to_string(), json!(prm.seed));
        outputs.insert("pseudo_f".to_string(), json!(pseudo_f));
        outputs.insert("cluster_sizes".to_string(), json!(counts));
        outputs.insert("cluster_means".to_string(), json!(cluster_means));
        if let Some(curve) = k_curve {
            outputs.insert("pseudo_f_curve".to_string(), json!(curve));
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── Clustering dispatch ────────────────────────────────────────────────────────

fn cluster(z: &[Vec<f64>], k: usize, method: Method, seed: u64) -> Vec<usize> {
    match method {
        Method::Kmeans => kmeans(z, k, seed),
        Method::Kmedoids => k_medoids(z, k, seed),
    }
}

// ── k-means (seeded k-means++ + Lloyd) ──────────────────────────────────────────

/// Several seeded restarts of k-means++ + Lloyd; the lowest-WCSS assignment
/// wins. Deterministic: the same seed always walks the same sequence of
/// restarts.
fn kmeans(z: &[Vec<f64>], k: usize, seed: u64) -> Vec<usize> {
    let mut state = seed ^ 0x9E37_79B9_7F4A_7C15;
    let mut best: Option<(f64, Vec<usize>)> = None;
    for _ in 0..KMEANS_RESTARTS {
        let (assign, wcss) = kmeans_once(z, k, &mut state);
        if best.as_ref().is_none_or(|(bw, _)| wcss < *bw) {
            best = Some((wcss, assign));
        }
    }
    best.expect("at least one restart").1
}

/// One k-means run: seeded k-means++ init, then Lloyd iterations until
/// assignments stop changing (or `KMEANS_MAX_ITERS`). Ties in nearest-centroid
/// assignment go to the lowest cluster index (deterministic scan order).
fn kmeans_once(z: &[Vec<f64>], k: usize, state: &mut u64) -> (Vec<usize>, f64) {
    let n = z.len();
    let nf = z.first().map(|r| r.len()).unwrap_or(0);
    let mut centroids: Vec<Vec<f64>> = kmeans_plus_plus_init(z, k, state)
        .into_iter()
        .map(|i| z[i].clone())
        .collect();

    let mut assign = vec![0usize; n];
    for _ in 0..KMEANS_MAX_ITERS {
        let mut changed = false;
        for (i, row) in z.iter().enumerate() {
            let c = nearest(row, &centroids);
            if assign[i] != c {
                assign[i] = c;
                changed = true;
            }
        }

        let mut counts = vec![0usize; k];
        let mut sums = vec![vec![0.0f64; nf]; k];
        for (i, row) in z.iter().enumerate() {
            let c = assign[i];
            counts[c] += 1;
            for d in 0..nf {
                sums[c][d] += row[d];
            }
        }
        let mut new_centroids = centroids.clone();
        for c in 0..k {
            if counts[c] == 0 {
                // Empty cluster: reseed deterministically at the point
                // currently farthest from its own centroid.
                let far = (0..n)
                    .max_by(|&a, &b| {
                        let da = dist2(&z[a], &centroids[assign[a]]);
                        let db = dist2(&z[b], &centroids[assign[b]]);
                        da.total_cmp(&db)
                    })
                    .unwrap_or(0);
                new_centroids[c] = z[far].clone();
            } else {
                new_centroids[c] = sums[c].iter().map(|s| s / counts[c] as f64).collect();
            }
        }
        let stable = !changed
            && new_centroids
                .iter()
                .zip(&centroids)
                .all(|(a, b)| dist2(a, b) < 1e-24);
        centroids = new_centroids;
        if stable {
            break;
        }
    }

    let wcss: f64 = z
        .iter()
        .enumerate()
        .map(|(i, row)| dist2(row, &centroids[assign[i]]))
        .sum();
    (assign, wcss)
}

/// Seeded k-means++ initialization: the first center is picked uniformly (via
/// the seeded generator), then each subsequent center is picked with
/// probability proportional to its squared distance to the nearest existing
/// center — all draws come from the seeded splitmix64 stream, never an
/// unseeded RNG or wall clock.
fn kmeans_plus_plus_init(z: &[Vec<f64>], k: usize, state: &mut u64) -> Vec<usize> {
    let n = z.len();
    let mut centers = Vec::with_capacity(k);
    centers.push((splitmix(state) % n as u64) as usize);
    let mut d2 = vec![f64::INFINITY; n];
    while centers.len() < k {
        let last = *centers.last().unwrap();
        for (i, row) in z.iter().enumerate() {
            let d = dist2(row, &z[last]);
            if d < d2[i] {
                d2[i] = d;
            }
        }
        let total: f64 = d2.iter().sum();
        let next = if total <= 0.0 {
            (0..n)
                .find(|i| !centers.contains(i))
                .unwrap_or(centers.len() % n)
        } else {
            let r = next_unit_f64(state) * total;
            let mut cum = 0.0;
            let mut chosen = n - 1;
            for (i, &d) in d2.iter().enumerate() {
                cum += d;
                if cum >= r {
                    chosen = i;
                    break;
                }
            }
            chosen
        };
        centers.push(next);
    }
    centers
}

fn nearest(row: &[f64], centroids: &[Vec<f64>]) -> usize {
    let mut best = 0usize;
    let mut best_d = f64::INFINITY;
    for (c, cent) in centroids.iter().enumerate() {
        let d = dist2(row, cent);
        if d < best_d {
            best_d = d;
            best = c;
        }
    }
    best
}

// ── k-medoids (PAM) ──────────────────────────────────────────────────────────

/// Deterministic k-medoids with several seeded restarts on a precomputed
/// distance matrix; the partition with the lowest total point-to-medoid cost
/// wins. O(n^2) — intended for moderate feature counts.
fn k_medoids(z: &[Vec<f64>], k: usize, seed: u64) -> Vec<usize> {
    let n = z.len();
    let mut d = vec![0.0f64; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            let v = dist2(&z[i], &z[j]).sqrt();
            d[i * n + j] = v;
            d[j * n + i] = v;
        }
    }

    let mut state = seed ^ 0x2545_F491_4F6C_DD1D;
    let mut best: Option<(f64, Vec<usize>)> = None;
    for _ in 0..PAM_RESTARTS {
        let (assign, medoids) = pam_once(&d, n, k, &mut state);
        let cost: f64 = (0..n).map(|i| d[i * n + medoids[assign[i]]]).sum();
        if best.as_ref().is_none_or(|(bc, _)| cost < *bc) {
            best = Some((cost, assign));
        }
    }
    best.expect("at least one restart").1
}

fn pam_once(d: &[f64], n: usize, k: usize, state: &mut u64) -> (Vec<usize>, Vec<usize>) {
    let mut medoids: Vec<usize> = Vec::with_capacity(k);
    while medoids.len() < k {
        let cand = (splitmix(state) % n as u64) as usize;
        if !medoids.contains(&cand) {
            medoids.push(cand);
        }
    }

    let mut assign = vec![0usize; n];
    for _ in 0..PAM_MAX_ITERS {
        let mut changed = false;
        for i in 0..n {
            let mut best = 0usize;
            let mut best_d = f64::INFINITY;
            for (m, &med) in medoids.iter().enumerate() {
                let dd = d[i * n + med];
                if dd < best_d {
                    best_d = dd;
                    best = m;
                }
            }
            if assign[i] != best {
                assign[i] = best;
                changed = true;
            }
        }
        let mut new_medoids = medoids.clone();
        for (m, med) in medoids.iter().enumerate() {
            let members: Vec<usize> = (0..n).filter(|&i| assign[i] == m).collect();
            if members.is_empty() {
                continue;
            }
            let mut best = *med;
            let mut best_cost = f64::INFINITY;
            for &cand in &members {
                let cost: f64 = members.iter().map(|&i| d[cand * n + i]).sum();
                if cost < best_cost {
                    best_cost = cost;
                    best = cand;
                }
            }
            new_medoids[m] = best;
        }
        if new_medoids == medoids && !changed {
            break;
        }
        medoids = new_medoids;
    }
    (assign, medoids)
}

// ── Seeded RNG (splitmix64) ─────────────────────────────────────────────────────

fn splitmix(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn next_unit_f64(state: &mut u64) -> f64 {
    (splitmix(state) >> 11) as f64 / (1u64 << 53) as f64
}

// ── Statistics ────────────────────────────────────────────────────────────────

fn dist2(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).powi(2)).sum()
}

/// Z-score standardization per field (std <= 0 collapses to 0, not NaN).
fn standardize(attrs: &[Vec<f64>], nf: usize) -> Vec<Vec<f64>> {
    let n = attrs.len() as f64;
    let mut mean = vec![0.0f64; nf];
    for row in attrs {
        for (d, v) in row.iter().enumerate() {
            mean[d] += v;
        }
    }
    for m in mean.iter_mut() {
        *m /= n;
    }
    let mut var = vec![0.0f64; nf];
    for row in attrs {
        for (d, v) in row.iter().enumerate() {
            var[d] += (v - mean[d]).powi(2);
        }
    }
    for v in var.iter_mut() {
        *v /= n;
    }
    let sd: Vec<f64> = var.iter().map(|v| v.sqrt()).collect();
    attrs
        .iter()
        .map(|row| {
            row.iter()
                .enumerate()
                .map(|(d, v)| {
                    if sd[d] <= 1e-12 {
                        0.0
                    } else {
                        (v - mean[d]) / sd[d]
                    }
                })
                .collect()
        })
        .collect()
}

/// Calinski-Harabasz pseudo-F: (between-cluster SS / (k-1)) / (within-cluster
/// SS / (n-k)). Higher is better; used to pick k automatically.
fn calinski_harabasz(z: &[Vec<f64>], assign: &[usize], k: usize) -> f64 {
    let n = z.len();
    if k < 2 || n <= k {
        return 0.0;
    }
    let nf = z.first().map(|r| r.len()).unwrap_or(0);

    let mut overall = vec![0.0f64; nf];
    for row in z {
        for (d, v) in row.iter().enumerate() {
            overall[d] += v;
        }
    }
    for v in overall.iter_mut() {
        *v /= n as f64;
    }

    let mut counts = vec![0usize; k];
    let mut sums = vec![vec![0.0f64; nf]; k];
    for (i, row) in z.iter().enumerate() {
        let c = assign[i];
        counts[c] += 1;
        for (d, v) in row.iter().enumerate() {
            sums[c][d] += v;
        }
    }
    let centroids: Vec<Vec<f64>> = (0..k)
        .map(|c| {
            if counts[c] == 0 {
                vec![0.0; nf]
            } else {
                sums[c].iter().map(|s| s / counts[c] as f64).collect()
            }
        })
        .collect();

    let between: f64 = (0..k)
        .map(|c| counts[c] as f64 * dist2(&centroids[c], &overall))
        .sum();
    let within: f64 = z
        .iter()
        .enumerate()
        .map(|(i, row)| dist2(row, &centroids[assign[i]]))
        .sum();

    if within <= 0.0 {
        return 0.0;
    }
    (between / (k as f64 - 1.0)) / (within / (n - k) as f64)
}

/// Per-point silhouette score in `[-1, 1]`: `(b - a) / max(a, b)` where `a` is
/// the mean distance to same-cluster points and `b` is the mean distance to
/// the nearest other cluster. O(n^2).
fn silhouette_scores(z: &[Vec<f64>], assign: &[usize], k: usize) -> Vec<f64> {
    let n = z.len();
    let mut sil = vec![0.0f64; n];
    if k < 2 || n < 3 {
        return sil;
    }
    for i in 0..n {
        let ci = assign[i];
        let mut same_sum = 0.0;
        let mut same_count = 0usize;
        let mut other_sum = vec![0.0f64; k];
        let mut other_count = vec![0usize; k];
        for j in 0..n {
            if i == j {
                continue;
            }
            let d = dist2(&z[i], &z[j]).sqrt();
            let cj = assign[j];
            if cj == ci {
                same_sum += d;
                same_count += 1;
            } else {
                other_sum[cj] += d;
                other_count[cj] += 1;
            }
        }
        let a = if same_count > 0 {
            same_sum / same_count as f64
        } else {
            0.0
        };
        let b = (0..k)
            .filter(|&c| c != ci && other_count[c] > 0)
            .map(|c| other_sum[c] / other_count[c] as f64)
            .fold(f64::INFINITY, f64::min);
        sil[i] = if b.is_finite() && a.max(b) > 0.0 {
            (b - a) / a.max(b)
        } else {
            0.0
        };
    }
    sil
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    fields: Vec<String>,
    num_clusters: Option<usize>,
    method: Method,
    seed: u64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let fields: Vec<String> = require_str(args, "fields")?
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_string)
        .collect();
    if fields.is_empty() {
        return Err(ToolError::Validation(
            "'fields' must list at least one analysis field".to_string(),
        ));
    }

    let num_clusters = match parse_optional_u64(args, "num_clusters")? {
        None | Some(0) => None,
        Some(v) => {
            if v < MIN_K as u64 {
                return Err(ToolError::Validation(format!(
                    "'num_clusters' must be 0 (automatic) or at least {MIN_K}, got {v}"
                )));
            }
            Some(v as usize)
        }
    };

    let method = match parse_optional_str(args, "method")?.map(|s| s.trim().to_lowercase()) {
        None => Method::Kmeans,
        Some(s) if s.is_empty() || s == "kmeans" => Method::Kmeans,
        Some(s) if s == "kmedoids" => Method::Kmedoids,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "'method' must be kmeans|kmedoids, got '{other}'"
            )))
        }
    };

    let seed = parse_optional_u64(args, "seed")?.unwrap_or(1);

    Ok(Params {
        fields,
        num_clusters,
        method,
        seed,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn parse_optional_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64().or_else(|| n.as_f64().map(|f| f as u64))),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
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
    use wbvector::{memory_store, Coord, Geometry, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Points at (x, y, a, b) where (a, b) is the attribute-space location.
    fn layer_of(pts: &[(f64, f64, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("a", FieldType::Float));
        l.add_field(FieldDef::new("b", FieldType::Float));
        for (x, y, a, b) in pts {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(*x, *y))),
                &[("a", (*a).into()), ("b", (*b).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = MultivariateClusteringTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// Points per blob in [`three_blobs`]; large enough that Calinski-Harabasz
    /// reliably prefers the true k=3 over over-fitting to near-singleton
    /// sub-clusters.
    const BLOB_SIZE: usize = 25;

    /// Three well-separated Gaussian blobs in attribute space (via Box-Muller
    /// over the tool's own seeded splitmix64 stream, so this is deterministic
    /// without pulling in a `rand` dependency). Centers are ~85-120 units
    /// apart with an internal std of ~1.5, so the true k=3 structure recovers
    /// with the right memberships (chunk `i` -> blob `i / BLOB_SIZE`).
    fn three_blobs() -> Vec<(f64, f64, f64, f64)> {
        let mut state = 42u64;
        let centers = [(0.0, 0.0), (60.0, 60.0), (-60.0, 60.0)];
        let mut pts = Vec::new();
        for (bi, &(ca, cb)) in centers.iter().enumerate() {
            for i in 0..BLOB_SIZE {
                let u1 = next_unit_f64(&mut state).max(1e-9);
                let u2 = next_unit_f64(&mut state);
                let r = (-2.0 * u1.ln()).sqrt();
                let theta = 2.0 * std::f64::consts::PI * u2;
                let da = r * theta.cos() * 1.5;
                let db = r * theta.sin() * 1.5;
                pts.push((bi as f64 * 200.0 + i as f64, bi as f64, ca + da, cb + db));
            }
        }
        pts
    }

    #[test]
    fn recovers_three_blobs_with_kmeans() {
        let input = layer_of(&three_blobs());
        let (out, layer) = run(json!({
            "input": input, "fields": "a,b", "num_clusters": 3, "seed": 7
        }));
        assert_eq!(out.outputs["num_clusters"], json!(3));
        let gi = layer.schema.field_index("cluster_id").unwrap();
        let clusters: Vec<i64> = layer
            .iter()
            .map(|f| f.attributes[gi].as_i64().unwrap())
            .collect();
        // Each blob (BLOB_SIZE consecutive rows) stays in a single cluster,
        // and the three blobs land in three distinct clusters.
        for chunk in clusters.chunks(BLOB_SIZE) {
            assert!(chunk.iter().all(|&c| c == chunk[0]));
        }
        assert_ne!(clusters[0], clusters[BLOB_SIZE]);
        assert_ne!(clusters[0], clusters[2 * BLOB_SIZE]);
        assert_ne!(clusters[BLOB_SIZE], clusters[2 * BLOB_SIZE]);
    }

    #[test]
    fn recovers_three_blobs_with_kmedoids() {
        let input = layer_of(&three_blobs());
        let (_out, layer) = run(json!({
            "input": input, "fields": "a,b", "num_clusters": 3, "method": "kmedoids", "seed": 3
        }));
        let gi = layer.schema.field_index("cluster_id").unwrap();
        let clusters: Vec<i64> = layer
            .iter()
            .map(|f| f.attributes[gi].as_i64().unwrap())
            .collect();
        for chunk in clusters.chunks(BLOB_SIZE) {
            assert!(chunk.iter().all(|&c| c == chunk[0]));
        }
        assert_ne!(clusters[0], clusters[BLOB_SIZE]);
        assert_ne!(clusters[0], clusters[2 * BLOB_SIZE]);
        assert_ne!(clusters[BLOB_SIZE], clusters[2 * BLOB_SIZE]);
    }

    /// Same seed -> identical labels across independent runs (determinism).
    #[test]
    fn deterministic_by_seed() {
        let input = layer_of(&three_blobs());
        let labels = |seed: u64| -> Vec<i64> {
            let (_o, l) = run(json!({
                "input": input, "fields": "a,b", "num_clusters": 3, "seed": seed
            }));
            let gi = l.schema.field_index("cluster_id").unwrap();
            l.iter()
                .map(|f| f.attributes[gi].as_i64().unwrap())
                .collect()
        };
        assert_eq!(labels(11), labels(11));
    }

    /// Optimal-k mode (num_clusters omitted) picks 3 on well-separated blobs.
    #[test]
    fn optimal_k_picks_three() {
        let input = layer_of(&three_blobs());
        let (out, _layer) = run(json!({ "input": input, "fields": "a,b", "seed": 1 }));
        assert_eq!(out.outputs["num_clusters"], json!(3));
        assert!(out.outputs["pseudo_f_curve"].as_array().unwrap().len() >= 2);
    }

    #[test]
    fn rejects_bad_parameters() {
        let input = layer_of(&three_blobs());

        let args: ToolArgs = serde_json::from_value(json!({ "fields": "a,b" })).unwrap();
        assert!(
            MultivariateClusteringTool.validate(&args).is_err(),
            "missing input"
        );

        let args: ToolArgs = serde_json::from_value(json!({ "input": input })).unwrap();
        assert!(
            MultivariateClusteringTool.validate(&args).is_err(),
            "missing fields"
        );

        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "fields": "a,b", "method": "nope" }))
                .unwrap();
        assert!(
            MultivariateClusteringTool.validate(&args).is_err(),
            "bad method"
        );

        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "fields": "a,b", "num_clusters": 1 }))
                .unwrap();
        assert!(
            MultivariateClusteringTool.validate(&args).is_err(),
            "num_clusters below MIN_K"
        );
    }

    #[test]
    fn rejects_unknown_field() {
        let input = layer_of(&three_blobs());
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "fields": "a,nonexistent" })).unwrap();
        // Unknown field is only caught at run-time (needs the loaded schema).
        assert!(MultivariateClusteringTool.validate(&args).is_ok());
        let out = MultivariateClusteringTool.run(&args, &ctx());
        assert!(out.is_err(), "unknown field should fail at run time");
    }
}
