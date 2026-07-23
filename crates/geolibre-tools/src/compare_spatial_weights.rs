//! GeoLibre tool: rank spatial-weights conceptualizations by clustering strength.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Compare Neighborhood Conceptualizations*
//! (Spatial Statistics). Choosing the spatial-weights matrix **W** (KNN vs
//! distance band vs contiguity vs Delaunay, and the bandwidth/k) is the single
//! most consequential and least-guided decision in every clustering statistic
//! the repo already ships (`global_morans_i`, `local_morans_i_lisa`,
//! `getis_ord_gi_star`, `incremental_spatial_autocorrelation`). This tool takes
//! the guesswork out of it: for each analysis field × each candidate
//! conceptualization it builds the neighbour structure (row-standardized),
//! computes Global Moran's *I* with its normal z-score and two-sided p-value,
//! and emits a ranked comparison table plus the single best method — the
//! conceptualization whose neighborhood definition captures the strongest,
//! most statistically significant spatial clustering.
//!
//! Candidate conceptualizations:
//!
//! * `knn` — the k nearest features (kd-tree),
//! * `fixed_distance_band` — every feature within `threshold_distance`,
//! * `contiguity_edges` — polygons sharing a border segment (rook),
//! * `contiguity_edges_corners` — polygons sharing a border **or** corner (queen),
//! * `delaunay` — features adjacent in the Delaunay triangulation.
//!
//! Contiguity methods are dropped automatically when the layer has no polygons;
//! point-based methods are dropped when the layer has no usable point geometry;
//! any method that errors is skipped rather than failing the whole run. Rows
//! with fewer than three features or zero variance are emitted with null
//! statistics. Output is a geometry-less attribute table (or a CSV when the
//! path ends in `.csv`). Deterministic.

use std::collections::{BTreeMap, HashMap, HashSet};

use geo::{Centroid, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Feature, FieldDef, FieldType, FieldValue, Geometry, Layer, Ring};

use crate::common::write_text_output;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Per-origin adjacency: `adj[i]` lists `(neighbour index, weight)` pairs.
type Adjacency = Vec<Vec<(usize, f64)>>;

pub struct CompareSpatialWeightsTool;

impl Tool for CompareSpatialWeightsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "compare_spatial_weights",
            display_name: "Compare Spatial Weights",
            summary: "Rank candidate spatial-weights conceptualizations (KNN, fixed distance band, contiguity edges / edges+corners, Delaunay) by how strongly each captures spatial clustering of one or more analysis fields, using row-standardized Global Moran's I with its z-score and two-sided p-value — like ArcGIS Compare Neighborhood Conceptualizations. Emits a ranked comparison table and the single best method.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (points/lines/polygons; contiguity requires polygons).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output comparison table path — a CSV (extension .csv) or a geometry-less vector table. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "input_fields",
                    description: "Comma-separated list of numeric field names to analyse for spatial clustering.",
                    required: true,
                },
                ToolParamSpec {
                    name: "id_field",
                    description: "Field holding each feature's unique id (accepted for compatibility; not required to compute Moran's I).",
                    required: false,
                },
                ToolParamSpec {
                    name: "methods",
                    description: "Comma-separated candidate conceptualizations from: knn, fixed_distance_band, contiguity_edges, contiguity_edges_corners, delaunay. Default 'knn,contiguity_edges,contiguity_edges_corners'. Methods that do not apply to the geometry (or error) are skipped.",
                    required: false,
                },
                ToolParamSpec {
                    name: "number_of_neighbors",
                    description: "Number of nearest neighbours k for the knn method. Default 8.",
                    required: false,
                },
                ToolParamSpec {
                    name: "threshold_distance",
                    description: "Distance band cutoff in map units. Required only when fixed_distance_band is among the methods.",
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

        let layer = load_input_layer(input)?;
        let n = layer.features.len();

        // Optional id_field — validated for existence, unused in the statistic.
        if let Some(f) = &prm.id_field {
            if layer.schema.field_index(f).is_none() {
                return Err(ToolError::Validation(format!("id_field '{f}' not found")));
            }
        }

        // Resolve each analysis field's values once. `None` where the field is
        // absent or any feature's value is non-numeric (that field is skipped).
        let mut fields: Vec<(String, Option<Vec<f64>>)> = Vec::new();
        for name in &prm.input_fields {
            let idx = layer
                .schema
                .field_index(name)
                .ok_or_else(|| ToolError::Validation(format!("input field '{name}' not found")))?;
            let mut vals = Vec::with_capacity(n);
            let mut ok = true;
            for feat in &layer.features {
                match feat.attributes.get(idx).and_then(FieldValue::as_f64) {
                    Some(v) if v.is_finite() => vals.push(v),
                    _ => {
                        ok = false;
                        break;
                    }
                }
            }
            fields.push((name.clone(), if ok { Some(vals) } else { None }));
        }

        // Build each requested method's row-standardized adjacency, skipping
        // methods that do not apply to this geometry (or otherwise error).
        let mut built: Vec<(Method, Adjacency, f64)> = Vec::new();
        for &method in &prm.methods {
            match build_adjacency(&layer, method, &prm) {
                Ok(mut adj) => {
                    let link_count: usize = adj.iter().map(Vec::len).sum();
                    let mean_neighbors = if n > 0 {
                        link_count as f64 / n as f64
                    } else {
                        0.0
                    };
                    row_standardize(&mut adj);
                    built.push((method, adj, mean_neighbors));
                }
                Err(e) => {
                    ctx.progress
                        .info(&format!("skipping method '{}': {e}", method.label()));
                }
            }
        }
        if built.is_empty() {
            return Err(ToolError::Execution(
                "no candidate conceptualization applied to this layer's geometry".to_string(),
            ));
        }

        // ── Score every (method, field) pair ──────────────────────────────────
        let mut results: Vec<Row> = Vec::new();
        for (method, adj, mean_neighbors) in &built {
            for (fname, fvals) in &fields {
                let stats = fvals.as_ref().and_then(|vals| morans_i(adj, vals, n));
                results.push(Row {
                    method: method.label().to_string(),
                    field: fname.clone(),
                    stats,
                    mean_neighbors: *mean_neighbors,
                });
            }
        }

        // Rank by descending z-score (null stats sink to the bottom).
        results.sort_by(|a, b| {
            let za = a.stats.map(|s| s.z_score).unwrap_or(f64::NEG_INFINITY);
            let zb = b.stats.map(|s| s.z_score).unwrap_or(f64::NEG_INFINITY);
            zb.partial_cmp(&za)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.method.cmp(&b.method))
                .then_with(|| a.field.cmp(&b.field))
        });

        // Best method: highest mean z across fields (ties broken by count).
        let mut per_method: BTreeMap<String, (f64, usize)> = BTreeMap::new();
        for r in &results {
            if let Some(s) = r.stats {
                let e = per_method.entry(r.method.clone()).or_insert((0.0, 0));
                e.0 += s.z_score;
                e.1 += 1;
            }
        }
        let (best_method, best_z_score) = per_method
            .iter()
            .map(|(m, (sum, cnt))| (m.clone(), sum / *cnt as f64))
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap_or_default();

        ctx.progress.info(&format!(
            "{n} feature(s); {} method(s) scored; best = '{best_method}' (mean z {best_z_score:.4})",
            built.len()
        ));

        // ── Emit the ranked comparison table ──────────────────────────────────
        let mut table = Layer::new("weights_comparison");
        table.add_field(FieldDef::new("method", FieldType::Text));
        table.add_field(FieldDef::new("field", FieldType::Text));
        table.add_field(FieldDef::new("morans_i", FieldType::Float));
        table.add_field(FieldDef::new("expected_i", FieldType::Float));
        table.add_field(FieldDef::new("z_score", FieldType::Float));
        table.add_field(FieldDef::new("p_value", FieldType::Float));
        table.add_field(FieldDef::new("mean_neighbors", FieldType::Float));

        let mut csv =
            String::from("method,field,morans_i,expected_i,z_score,p_value,mean_neighbors\n");
        let fmt = |v: Option<f64>| match v {
            Some(x) => format!("{x}"),
            None => String::new(),
        };
        for r in &results {
            let (mi, ei, z, p) = match r.stats {
                Some(s) => (
                    FieldValue::Float(s.morans_i),
                    FieldValue::Float(s.expected_i),
                    FieldValue::Float(s.z_score),
                    FieldValue::Float(s.p_value),
                ),
                None => (
                    FieldValue::Null,
                    FieldValue::Null,
                    FieldValue::Null,
                    FieldValue::Null,
                ),
            };
            table.push(Feature {
                fid: 0,
                geometry: None,
                attributes: vec![
                    FieldValue::Text(r.method.clone()),
                    FieldValue::Text(r.field.clone()),
                    mi,
                    ei,
                    z,
                    p,
                    FieldValue::Float(r.mean_neighbors),
                ],
            });
            let (smi, sei, sz, sp) = match r.stats {
                Some(s) => (
                    fmt(Some(s.morans_i)),
                    fmt(Some(s.expected_i)),
                    fmt(Some(s.z_score)),
                    fmt(Some(s.p_value)),
                ),
                None => (String::new(), String::new(), String::new(), String::new()),
            };
            csv.push_str(&format!(
                "{},{},{},{},{},{},{}\n",
                r.method, r.field, smi, sei, sz, sp, r.mean_neighbors
            ));
        }

        let out_path = match output {
            Some(path) if path.to_ascii_lowercase().ends_with(".csv") => {
                write_text_output(&csv, path)?;
                path.to_string()
            }
            other => write_or_store_layer(table, other)?,
        };

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("rows".to_string(), json!(results.len()));
        outputs.insert("best_method".to_string(), json!(best_method));
        outputs.insert("best_z_score".to_string(), json!(best_z_score));
        Ok(ToolRunResult { outputs })
    }
}

/// One row of the comparison table: a (method, field) pairing.
struct Row {
    method: String,
    field: String,
    stats: Option<MoranStats>,
    mean_neighbors: f64,
}

#[derive(Clone, Copy)]
struct MoranStats {
    morans_i: f64,
    expected_i: f64,
    z_score: f64,
    p_value: f64,
}

// ── Global Moran's I (row-standardized weights) ───────────────────────────────

/// Global Moran's *I* with normal-approximation z-score and two-sided p-value.
///
/// `adj[i]` holds `(j, w_ij)` pairs (already row-standardized). Returns `None`
/// when there are fewer than three features, zero variance, or degenerate
/// weights (no links / zero variance of I).
fn morans_i(adj: &[Vec<(usize, f64)>], values: &[f64], n: usize) -> Option<MoranStats> {
    if n < 3 || values.len() != n {
        return None;
    }
    let nf = n as f64;
    let mean = values.iter().sum::<f64>() / nf;
    let z: Vec<f64> = values.iter().map(|v| v - mean).collect();
    let denom: f64 = z.iter().map(|v| v * v).sum();
    if denom <= 0.0 {
        return None; // zero variance
    }

    // Sparse weight map plus row / column sums.
    let mut w: HashMap<(usize, usize), f64> = HashMap::new();
    let mut row_sum = vec![0.0f64; n];
    let mut col_sum = vec![0.0f64; n];
    for (i, row) in adj.iter().enumerate() {
        for &(j, wij) in row {
            if i == j || wij == 0.0 {
                continue;
            }
            *w.entry((i, j)).or_insert(0.0) += wij;
        }
    }
    for (&(i, j), &wij) in &w {
        row_sum[i] += wij;
        col_sum[j] += wij;
    }

    let s0: f64 = w.values().sum();
    if s0 <= 0.0 {
        return None;
    }

    // Cross-product numerator Σ_i Σ_j w_ij z_i z_j.
    let mut cross = 0.0;
    for (&(i, j), &wij) in &w {
        cross += wij * z[i] * z[j];
    }
    let morans = (nf / s0) * (cross / denom);
    let expected = -1.0 / (nf - 1.0);

    // S1 = 0.5 Σ_i Σ_j (w_ij + w_ji)^2 over the union of directed pairs.
    let mut pairs: HashSet<(usize, usize)> = HashSet::new();
    for &(i, j) in w.keys() {
        pairs.insert((i, j));
        pairs.insert((j, i));
    }
    let mut s1 = 0.0;
    for &(i, j) in &pairs {
        let wij = w.get(&(i, j)).copied().unwrap_or(0.0);
        let wji = w.get(&(j, i)).copied().unwrap_or(0.0);
        s1 += (wij + wji).powi(2);
    }
    s1 *= 0.5;

    // S2 = Σ_i ( Σ_j w_ij + Σ_j w_ji )^2.
    let mut s2 = 0.0;
    for i in 0..n {
        s2 += (row_sum[i] + col_sum[i]).powi(2);
    }

    let var = (nf * nf * s1 - nf * s2 + 3.0 * s0 * s0) / ((nf * nf - 1.0) * s0 * s0)
        - expected * expected;
    if !(var.is_finite() && var > 0.0) {
        return None;
    }
    let z_score = (morans - expected) / var.sqrt();
    let p_value = two_sided_p(z_score);

    Some(MoranStats {
        morans_i: morans,
        expected_i: expected,
        z_score,
        p_value,
    })
}

/// Two-sided p-value from a z-score: `2·(1 − Φ(|z|))`.
fn two_sided_p(z: f64) -> f64 {
    2.0 * (1.0 - normal_cdf(z.abs()))
}

/// Standard normal CDF `Φ(x) = 0.5·erfc(−x/√2)`.
fn normal_cdf(x: f64) -> f64 {
    0.5 * erfc(-x / std::f64::consts::SQRT_2)
}

fn erfc(x: f64) -> f64 {
    1.0 - erf(x)
}

/// Abramowitz & Stegun 7.1.26 approximation of the error function.
fn erf(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let y = 1.0
        - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t
            + 0.254829592)
            * t
            * (-x * x).exp();
    sign * y
}

// ── Neighbour builders (self-contained; copied from generate_spatial_weights) ──

fn build_adjacency(
    layer: &Layer,
    method: Method,
    prm: &Params,
) -> Result<Vec<Vec<(usize, f64)>>, ToolError> {
    match method {
        Method::Knn => neighbors_knn(layer, prm.number_of_neighbors),
        Method::FixedDistanceBand => {
            let t = prm.threshold_distance.ok_or_else(|| {
                ToolError::Validation(
                    "'threshold_distance' is required for fixed_distance_band".to_string(),
                )
            })?;
            neighbors_distance_band(layer, t)
        }
        Method::ContiguityEdges => neighbors_contiguity(layer, false),
        Method::ContiguityEdgesCorners => neighbors_contiguity(layer, true),
        Method::Delaunay => neighbors_delaunay(layer),
    }
}

/// Divide each origin's weights by their row sum so every origin sums to 1.
fn row_standardize(adj: &mut [Vec<(usize, f64)>]) {
    for row in adj.iter_mut() {
        let sum: f64 = row.iter().map(|(_, w)| *w).sum();
        if sum > 0.0 {
            for (_, w) in row.iter_mut() {
                *w /= sum;
            }
        }
    }
}

/// Representative point per feature (`None` when no usable geometry).
fn representative_points(layer: &Layer) -> Vec<Option<(f64, f64)>> {
    layer
        .features
        .iter()
        .map(|f| f.geometry.as_ref().and_then(representative_point))
        .collect()
}

/// Build a kd-tree over the representative points that exist.
fn build_tree(reps: &[Option<(f64, f64)>]) -> Result<KdTree<f64, usize, [f64; 2]>, ToolError> {
    let mut tree: KdTree<f64, usize, [f64; 2]> = KdTree::new(2);
    for (i, r) in reps.iter().enumerate() {
        if let Some((x, y)) = r {
            tree.add([*x, *y], i)
                .map_err(|e| ToolError::Execution(format!("kd-tree insert failed: {e:?}")))?;
        }
    }
    Ok(tree)
}

fn neighbors_knn(layer: &Layer, k: usize) -> Result<Vec<Vec<(usize, f64)>>, ToolError> {
    let reps = representative_points(layer);
    if reps.iter().all(Option::is_none) {
        return Err(ToolError::Execution(
            "input has no usable point geometry for knn".to_string(),
        ));
    }
    let tree = build_tree(&reps)?;
    let n = reps.len();
    let mut out = vec![Vec::new(); n];
    for (i, rep) in reps.iter().enumerate() {
        let Some((x, y)) = rep else { continue };
        // Ask for k+1 so a feature can drop its own self-match.
        let found = tree
            .nearest(&[*x, *y], k + 1, &squared_euclidean)
            .map_err(|e| ToolError::Execution(format!("kd-tree query failed: {e:?}")))?;
        for (_d2, &j) in found.into_iter().filter(|(_, &j)| j != i).take(k) {
            out[i].push((j, 1.0));
        }
    }
    Ok(out)
}

fn neighbors_distance_band(
    layer: &Layer,
    threshold: f64,
) -> Result<Vec<Vec<(usize, f64)>>, ToolError> {
    let reps = representative_points(layer);
    if reps.iter().all(Option::is_none) {
        return Err(ToolError::Execution(
            "input has no usable point geometry".to_string(),
        ));
    }
    let tree = build_tree(&reps)?;
    let r2 = threshold * threshold;
    let n = reps.len();
    let mut out = vec![Vec::new(); n];
    for (i, rep) in reps.iter().enumerate() {
        let Some((x, y)) = rep else { continue };
        let found = tree
            .within(&[*x, *y], r2, &squared_euclidean)
            .map_err(|e| ToolError::Execution(format!("kd-tree query failed: {e:?}")))?;
        for (_d2, &j) in found {
            if j == i {
                continue;
            }
            out[i].push((j, 1.0));
        }
    }
    Ok(out)
}

fn neighbors_delaunay(layer: &Layer) -> Result<Vec<Vec<(usize, f64)>>, ToolError> {
    let reps = representative_points(layer);
    if reps.iter().all(Option::is_none) {
        return Err(ToolError::Execution(
            "input has no usable point geometry for delaunay".to_string(),
        ));
    }
    // Collect the indices that have a representative point, keeping a mapping.
    let present: Vec<usize> = reps
        .iter()
        .enumerate()
        .filter_map(|(i, r)| r.map(|_| i))
        .collect();
    let pts: Vec<(f64, f64)> = present.iter().map(|&i| reps[i].unwrap()).collect();
    let n = reps.len();
    let mut sets: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    for t in delaunay(&pts) {
        for &(a, b) in &[(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let (ga, gb) = (present[a], present[b]);
            sets[ga].insert(gb);
            sets[gb].insert(ga);
        }
    }
    Ok(sets
        .into_iter()
        .map(|s| s.into_iter().map(|j| (j, 1.0)).collect())
        .collect())
}

/// Polygon contiguity via shared edges (rook) and, when `include_corners`,
/// shared vertices too (queen). Symmetric by construction.
fn neighbors_contiguity(
    layer: &Layer,
    include_corners: bool,
) -> Result<Vec<Vec<(usize, f64)>>, ToolError> {
    let n = layer.features.len();
    let mut edge_feats: HashMap<(Key, Key), HashSet<usize>> = HashMap::new();
    let mut vert_feats: HashMap<Key, HashSet<usize>> = HashMap::new();
    let mut any_poly = false;

    for (fidx, feature) in layer.features.iter().enumerate() {
        let Some(geom) = feature.geometry.as_ref() else {
            continue;
        };
        let rings = polygon_rings(geom);
        if rings.is_empty() {
            continue;
        }
        any_poly = true;
        for ring in &rings {
            let m = ring.len();
            for i in 0..m {
                let a = ring[i];
                let b = ring[(i + 1) % m];
                vert_feats.entry(key(a)).or_default().insert(fidx);
                if key(a) == key(b) {
                    continue;
                }
                edge_feats.entry(edge_key(a, b)).or_default().insert(fidx);
            }
        }
    }
    if !any_poly {
        return Err(ToolError::Execution(
            "contiguity requires polygon geometry, but none was found".to_string(),
        ));
    }

    let mut sets: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    let link = |feats: &HashSet<usize>, sets: &mut Vec<HashSet<usize>>| {
        if feats.len() < 2 {
            return;
        }
        let list: Vec<usize> = feats.iter().copied().collect();
        for i in 0..list.len() {
            for j in (i + 1)..list.len() {
                sets[list[i]].insert(list[j]);
                sets[list[j]].insert(list[i]);
            }
        }
    };
    for feats in edge_feats.values() {
        link(feats, &mut sets);
    }
    if include_corners {
        for feats in vert_feats.values() {
            link(feats, &mut sets);
        }
    }

    Ok(sets
        .into_iter()
        .map(|s| s.into_iter().map(|j| (j, 1.0)).collect())
        .collect())
}

// ── Representative points ─────────────────────────────────────────────────────

fn representative_point(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => {
            let (sx, sy) = cs
                .iter()
                .fold((0.0, 0.0), |(ax, ay), c| (ax + c.x, ay + c.y));
            let k = cs.len() as f64;
            Some((sx / k, sy / k))
        }
        Geometry::LineString(cs) if !cs.is_empty() => {
            let ls = LineString::new(cs.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect());
            ls.centroid().map(|p| (p.x(), p.y()))
        }
        Geometry::MultiLineString(parts) => {
            let mls = geo::MultiLineString(
                parts
                    .iter()
                    .map(|cs| {
                        LineString::new(cs.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect())
                    })
                    .collect(),
            );
            mls.centroid().map(|p| (p.x(), p.y()))
        }
        Geometry::Polygon { .. } | Geometry::MultiPolygon(_) => to_multipolygon(geom)
            .and_then(|mp| mp.centroid())
            .map(|p| (p.x(), p.y())),
        _ => None,
    }
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

// ── Contiguity keys and rings ─────────────────────────────────────────────────

#[derive(Clone, Copy)]
struct P {
    x: f64,
    y: f64,
}

type Key = (u64, u64);

fn key(p: P) -> Key {
    (p.x.to_bits(), p.y.to_bits())
}

fn edge_key(a: P, b: P) -> (Key, Key) {
    let (ka, kb) = (key(a), key(b));
    if ka <= kb {
        (ka, kb)
    } else {
        (kb, ka)
    }
}

fn canonical(x: f64, y: f64) -> P {
    P { x, y }
}

/// All rings (exterior + interiors) of a polygon geometry as vertex chains
/// without the closing duplicate.
fn polygon_rings(geom: &Geometry) -> Vec<Vec<P>> {
    let ring_pts = |ring: &Ring| -> Vec<P> {
        let mut pts: Vec<P> = Vec::with_capacity(ring.len());
        for c in ring.coords() {
            let p = canonical(c.x, c.y);
            if pts.last().is_none_or(|l| key(*l) != key(p)) {
                pts.push(p);
            }
        }
        while pts.len() >= 2 && key(pts[0]) == key(*pts.last().unwrap()) {
            pts.pop();
        }
        pts
    };
    let mut out = Vec::new();
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            out.push(ring_pts(exterior));
            out.extend(interiors.iter().map(&ring_pts));
        }
        Geometry::MultiPolygon(parts) => {
            for (ext, holes) in parts {
                out.push(ring_pts(ext));
                out.extend(holes.iter().map(&ring_pts));
            }
        }
        _ => {}
    }
    out.retain(|r| r.len() >= 3);
    out
}

// ── Delaunay triangulation (Bowyer–Watson) ────────────────────────────────────

/// Returns triangles as index triples into `pts` (super-triangle removed).
fn delaunay(pts: &[(f64, f64)]) -> Vec<[usize; 3]> {
    let n = pts.len();
    if n < 3 {
        return Vec::new();
    }
    let (mut minx, mut miny, mut maxx, mut maxy) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
    for &(x, y) in pts {
        minx = minx.min(x);
        miny = miny.min(y);
        maxx = maxx.max(x);
        maxy = maxy.max(y);
    }
    let dmax = (maxx - minx).max(maxy - miny).max(1.0) * 20.0;
    let (mx, my) = ((minx + maxx) / 2.0, (miny + maxy) / 2.0);
    let mut verts: Vec<(f64, f64)> = pts.to_vec();
    let s0 = verts.len();
    verts.push((mx - dmax, my - dmax));
    verts.push((mx + dmax, my - dmax));
    verts.push((mx, my + dmax));
    let mut tris: Vec<[usize; 3]> = vec![[s0, s0 + 1, s0 + 2]];

    for i in 0..n {
        let p = verts[i];
        let mut bad: Vec<usize> = Vec::new();
        for (ti, t) in tris.iter().enumerate() {
            if in_circumcircle(p, verts[t[0]], verts[t[1]], verts[t[2]]) {
                bad.push(ti);
            }
        }
        let mut boundary: Vec<(usize, usize)> = Vec::new();
        for &ti in &bad {
            let t = tris[ti];
            for e in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                let shared = bad.iter().any(|&oj| {
                    oj != ti && {
                        let o = tris[oj];
                        let edges = [(o[0], o[1]), (o[1], o[2]), (o[2], o[0])];
                        edges
                            .iter()
                            .any(|&(u, v)| (u == e.0 && v == e.1) || (u == e.1 && v == e.0))
                    }
                });
                if !shared {
                    boundary.push(e);
                }
            }
        }
        bad.sort_unstable_by(|a, b| b.cmp(a));
        for &ti in &bad {
            tris.swap_remove(ti);
        }
        for (u, v) in boundary {
            tris.push([u, v, i]);
        }
    }
    tris.retain(|t| t.iter().all(|&v| v < s0));
    tris
}

fn in_circumcircle(p: (f64, f64), a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> bool {
    let ax = a.0 - p.0;
    let ay = a.1 - p.1;
    let bx = b.0 - p.0;
    let by = b.1 - p.1;
    let cx = c.0 - p.0;
    let cy = c.1 - p.1;
    let det = (ax * ax + ay * ay) * (bx * cy - cx * by) - (bx * bx + by * by) * (ax * cy - cx * ay)
        + (cx * cx + cy * cy) * (ax * by - bx * ay);
    let area2 = (b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1);
    if area2 > 0.0 {
        det > 0.0
    } else {
        det < 0.0
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Knn,
    FixedDistanceBand,
    ContiguityEdges,
    ContiguityEdgesCorners,
    Delaunay,
}

impl Method {
    fn label(self) -> &'static str {
        match self {
            Method::Knn => "knn",
            Method::FixedDistanceBand => "fixed_distance_band",
            Method::ContiguityEdges => "contiguity_edges",
            Method::ContiguityEdgesCorners => "contiguity_edges_corners",
            Method::Delaunay => "delaunay",
        }
    }
}

struct Params {
    input_fields: Vec<String>,
    id_field: Option<String>,
    methods: Vec<Method>,
    number_of_neighbors: usize,
    threshold_distance: Option<f64>,
}

fn parse_method(s: &str) -> Result<Method, ToolError> {
    match s.trim().to_ascii_lowercase().as_str() {
        "knn" | "k_nearest_neighbors" => Ok(Method::Knn),
        "fixed_distance_band" | "distance_band" => Ok(Method::FixedDistanceBand),
        "contiguity_edges" | "edges" | "rook" => Ok(Method::ContiguityEdges),
        "contiguity_edges_corners" | "edges_corners" | "queen" => Ok(Method::ContiguityEdgesCorners),
        "delaunay" | "triangulation" => Ok(Method::Delaunay),
        other => Err(ToolError::Validation(format!(
            "'methods' entry must be one of knn, fixed_distance_band, contiguity_edges, contiguity_edges_corners, delaunay (got '{other}')"
        ))),
    }
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let input_fields: Vec<String> = parse_optional_str(args, "input_fields")?
        .map(|s| {
            s.split(',')
                .map(str::trim)
                .filter(|t| !t.is_empty())
                .map(str::to_string)
                .collect::<Vec<_>>()
        })
        .unwrap_or_default();
    if input_fields.is_empty() {
        return Err(ToolError::Validation(
            "missing required parameter 'input_fields' (comma-separated numeric field names)"
                .to_string(),
        ));
    }

    let id_field = parse_optional_str(args, "id_field")?.map(str::to_string);

    let methods = match parse_optional_str(args, "methods")? {
        Some(s) => {
            let mut out = Vec::new();
            for tok in s.split(',').map(str::trim).filter(|t| !t.is_empty()) {
                let m = parse_method(tok)?;
                if !out.contains(&m) {
                    out.push(m);
                }
            }
            if out.is_empty() {
                return Err(ToolError::Validation(
                    "'methods' must list at least one conceptualization".to_string(),
                ));
            }
            out
        }
        None => vec![
            Method::Knn,
            Method::ContiguityEdges,
            Method::ContiguityEdgesCorners,
        ],
    };

    let number_of_neighbors = match args.get("number_of_neighbors") {
        None | Some(Value::Null) => 8,
        Some(Value::Number(n)) => n.as_u64().ok_or_else(|| {
            ToolError::Validation("'number_of_neighbors' must be a positive integer".into())
        })? as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 8,
        Some(Value::String(s)) => s.trim().parse::<usize>().map_err(|_| {
            ToolError::Validation("'number_of_neighbors' must be an integer".into())
        })?,
        Some(_) => {
            return Err(ToolError::Validation(
                "'number_of_neighbors' must be a number".into(),
            ))
        }
    };
    if methods.contains(&Method::Knn) && number_of_neighbors == 0 {
        return Err(ToolError::Validation(
            "'number_of_neighbors' must be >= 1 for knn".to_string(),
        ));
    }

    let threshold_distance = parse_optional_f64(args, "threshold_distance")?;
    if let Some(t) = threshold_distance {
        if !(t.is_finite() && t > 0.0) {
            return Err(ToolError::Validation(
                "'threshold_distance' must be positive".to_string(),
            ));
        }
    }
    if methods.contains(&Method::FixedDistanceBand) && threshold_distance.is_none() {
        return Err(ToolError::Validation(
            "'threshold_distance' is required when fixed_distance_band is among the methods"
                .to_string(),
        ));
    }

    Ok(Params {
        input_fields,
        id_field,
        methods,
        number_of_neighbors,
        threshold_distance,
    })
}

fn parse_optional_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => {
            Ok(Some(s.trim().parse::<f64>().map_err(|_| {
                ToolError::Validation(format!("'{key}' must be a number"))
            })?))
        }
        Some(_) => Err(ToolError::Validation(format!("'{key}' must be a number"))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Coord, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Points on a grid with an attribute value increasing with x -> strong
    /// positive spatial autocorrelation.
    fn clustered_point_layer() -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("id", FieldType::Integer));
        l.add_field(FieldDef::new("val", FieldType::Float));
        let mut id = 0i64;
        for gx in 0..6 {
            for gy in 0..6 {
                l.add_feature(
                    Some(Geometry::Point(Coord::xy(gx as f64, gy as f64))),
                    &[("id", id.into()), ("val", (gx as f64).into())],
                )
                .unwrap();
                id += 1;
            }
        }
        let vid = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&vid)
    }

    fn square(x0: f64, y0: f64, s: f64) -> Geometry {
        Geometry::polygon(
            vec![
                Coord::xy(x0, y0),
                Coord::xy(x0 + s, y0),
                Coord::xy(x0 + s, y0 + s),
                Coord::xy(x0, y0 + s),
            ],
            vec![],
        )
    }

    /// 3x3 grid of squares with an attribute increasing with column index.
    fn poly_grid_layer() -> String {
        let mut l = Layer::new("polys")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("val", FieldType::Float));
        for cx in 0..3 {
            for cy in 0..3 {
                l.add_feature(
                    Some(square(cx as f64 * 10.0, cy as f64 * 10.0, 10.0)),
                    &[("val", (cx as f64).into())],
                )
                .unwrap();
            }
        }
        let vid = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&vid)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CompareSpatialWeightsTool.run(&args, &ctx()).unwrap();
        let table = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, table)
    }

    type TestRow = (String, String, Option<f64>, Option<f64>, Option<f64>);

    fn rows(table: &Layer) -> Vec<TestRow> {
        let mi = table.schema.field_index("method").unwrap();
        let fi = table.schema.field_index("field").unwrap();
        let ii = table.schema.field_index("morans_i").unwrap();
        let zi = table.schema.field_index("z_score").unwrap();
        let ei = table.schema.field_index("expected_i").unwrap();
        table
            .iter()
            .map(|f| {
                (
                    f.attributes[mi].as_str().unwrap().to_string(),
                    f.attributes[fi].as_str().unwrap().to_string(),
                    f.attributes[ii].as_f64(),
                    f.attributes[zi].as_f64(),
                    f.attributes[ei].as_f64(),
                )
            })
            .collect()
    }

    /// A clustered point field yields positive Moran's I with z-score above the
    /// expected value under knn.
    #[test]
    fn clustered_field_is_positively_autocorrelated() {
        let input = clustered_point_layer();
        let (out, table) = run(json!({
            "input": input, "input_fields": "val", "methods": "knn",
            "number_of_neighbors": 4
        }));
        let r = rows(&table);
        let knn = r
            .iter()
            .find(|(m, f, ..)| m == "knn" && f == "val")
            .unwrap();
        let morans = knn.2.unwrap();
        let z = knn.3.unwrap();
        let expected = knn.4.unwrap();
        assert!(morans > 0.0, "Moran's I should be positive, got {morans}");
        assert!(
            z > expected,
            "z-score {z} should exceed expected I {expected}"
        );
        assert!(z > 0.0, "clustered field should have positive z, got {z}");
        // Best method selection returns a non-empty string.
        assert!(!out.outputs["best_method"].as_str().unwrap().is_empty());
    }

    /// A polygon layer supports contiguity; the default method set works.
    #[test]
    fn polygon_layer_supports_contiguity() {
        let input = poly_grid_layer();
        let (out, table) = run(json!({
            "input": input, "input_fields": "val"
        }));
        let methods: HashSet<String> = rows(&table).into_iter().map(|(m, ..)| m).collect();
        assert!(methods.contains("contiguity_edges"));
        assert!(methods.contains("contiguity_edges_corners"));
        assert!(!out.outputs["best_method"].as_str().unwrap().is_empty());
    }

    /// A point-only layer silently omits contiguity methods.
    #[test]
    fn point_layer_omits_contiguity() {
        let input = clustered_point_layer();
        let (_out, table) = run(json!({
            "input": input, "input_fields": "val",
            "methods": "knn,contiguity_edges,contiguity_edges_corners",
            "number_of_neighbors": 4
        }));
        let methods: HashSet<String> = rows(&table).into_iter().map(|(m, ..)| m).collect();
        assert!(methods.contains("knn"));
        assert!(
            !methods.contains("contiguity_edges"),
            "contiguity should be dropped for point-only input"
        );
    }

    #[test]
    fn rejects_missing_input() {
        let args: ToolArgs = serde_json::from_value(json!({ "input_fields": "val" })).unwrap();
        assert!(CompareSpatialWeightsTool.validate(&args).is_err());
    }

    #[test]
    fn rejects_missing_input_fields() {
        let input = clustered_point_layer();
        let args: ToolArgs = serde_json::from_value(json!({ "input": input })).unwrap();
        assert!(CompareSpatialWeightsTool.validate(&args).is_err());
    }

    #[test]
    fn rejects_unknown_method() {
        let input = clustered_point_layer();
        let args: ToolArgs = serde_json::from_value(
            json!({ "input": input, "input_fields": "val", "methods": "nonsense" }),
        )
        .unwrap();
        assert!(CompareSpatialWeightsTool.validate(&args).is_err());
    }
}
