//! GeoLibre tool: Moran eigenvector spatial filtering (MESF).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Decompose Spatial Structure (MESF)*
//! and *Create Spatial Component Explanatory Variables* (Spatial Statistics).
//! It manufactures **Moran eigenvectors** — the synthetic, mutually orthogonal
//! spatial predictors ("spatial filters") that let an ordinary OLS/GLR model
//! absorb residual spatial autocorrelation without moving to a full
//! spatial-lag/GWR model.
//!
//! The repo already ships `generate_spatial_weights_matrix`, which emits the
//! neighbour/weights structure **W** but never decomposes it. This tool takes
//! the same conceptualizations, builds a symmetric binary connectivity matrix
//! **C**, double-centres it, eigendecomposes it, and keeps the eigenvectors
//! carrying meaningful positive spatial autocorrelation.
//!
//! ## How it works
//!
//! 1. Build a symmetric binary connectivity matrix **C** (`C_ij ∈ {0,1}`,
//!    `C_ii = 0`) from `method`: `contiguity_edges` (rook) /
//!    `contiguity_edges_corners` (queen) — already symmetric — or `knn`,
//!    symmetrized so `C_ij = C_ji = 1` when *either* feature ranks the other a
//!    neighbour.
//! 2. Double-centre with `M = I − 11ᵀ/n`: `B = M C M`. Because C is symmetric
//!    this is the closed form `B_ij = C_ij − r_i − r_j + g`, where `r_i` is
//!    row-sum(C)/n and `g` is grand-sum(C)/n² — an O(n²) build, no matrix
//!    products.
//! 3. Symmetric eigendecomposition of **B** by **cyclic Jacobi rotations** (pure
//!    Rust, no LAPACK/`nalgebra`): sweep the off-diagonal entries, applying the
//!    plane rotation that zeros each `B_pq` to B (both sides) and to an
//!    accumulating identity V (one side), until the off-diagonal Frobenius norm
//!    drops below `1e-10·‖B‖` or 100 sweeps elapse. Eigenvalues are `diag(B)`;
//!    eigenvector *k* is column *k* of V, already orthonormal.
//! 4. Each eigenvector `E_k` has Moran's I `= (n / S0)·λ_k` where `S0` is the
//!    grand sum of C. Keep those with `I ≥ min_autocorrelation` (positive
//!    autocorrelation), rank them descending by Moran's I, and cap at
//!    `max_components`.
//! 5. Append the kept eigenvectors as numeric fields `MEV1, MEV2, …` (ranked
//!    order). Each column is already unit length and mean-zero (orthogonal to
//!    the constant vector, which B sends to eigenvalue 0).
//!
//! Deterministic — no RNG, no clocks. The Jacobi solver is O(n³); inputs above
//! 1500 features are rejected with a suggestion to subset.

use std::collections::{BTreeMap, HashMap, HashSet};

use geo::{Centroid, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Hard cap on feature count — the Jacobi eigensolver is O(n³) and a dense n×n
/// matrix is materialized.
const MAX_FEATURES: usize = 1500;

pub struct SpatialEigenvectorFilteringTool;

impl Tool for SpatialEigenvectorFilteringTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "spatial_eigenvector_filtering",
            display_name: "Spatial Eigenvector Filtering (MESF)",
            summary: "Generate Moran eigenvectors (spatial filters) from a spatial weights matrix and append the autocorrelated components as explanatory-variable fields (MEV1, MEV2, …), like ArcGIS Decompose Spatial Structure / Create Spatial Component Explanatory Variables. Builds a symmetric connectivity matrix (contiguity or symmetrized KNN), double-centres it (M C M), and eigendecomposes it with a pure-Rust cyclic-Jacobi solver; keeps eigenvectors whose Moran's I ≥ min_autocorrelation, ranked descending and capped at max_components.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (points/lines/polygons; contiguity methods require polygons).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path — the input features with the kept eigenvectors appended as MEV1, MEV2, …. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "Conceptualization of spatial relationships: 'contiguity_edges' (rook), 'contiguity_edges_corners' (queen), or 'knn' (symmetrized). Default 'contiguity_edges_corners'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "number_of_neighbors",
                    description: "Number of nearest neighbours k for method = knn. Default 8.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_autocorrelation",
                    description: "Minimum Moran's I an eigenvector must have to be kept (positive spatial autocorrelation). Default 0.25.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_components",
                    description: "Maximum number of eigenvectors to emit, taken in descending Moran's I order. Default 15.",
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
        let n = layer.features.len();
        if n < 3 {
            return Err(ToolError::Validation(format!(
                "spatial_eigenvector_filtering needs at least 3 features (got {n})"
            )));
        }
        if n > MAX_FEATURES {
            return Err(ToolError::Validation(format!(
                "input has {n} features, above the {MAX_FEATURES}-feature limit for the O(n³) Jacobi eigensolver; subset the input to a smaller area first"
            )));
        }

        // ── 1. Symmetric binary connectivity matrix C (n×n, C_ii = 0) ──────────
        let adjacency: Vec<Vec<(usize, f64)>> = match prm.method {
            Method::Knn => neighbors_knn(&layer, prm.number_of_neighbors)?,
            Method::ContiguityEdges => neighbors_contiguity(&layer, false)?,
            Method::ContiguityEdgesCorners => neighbors_contiguity(&layer, true)?,
        };
        let mut c = vec![0.0f64; n * n];
        for (i, row) in adjacency.iter().enumerate() {
            for &(j, _w) in row {
                if i == j {
                    continue;
                }
                c[i * n + j] = 1.0;
                c[j * n + i] = 1.0; // symmetrize (no-op for contiguity)
            }
        }

        // Grand sum S0 and per-row sums of C.
        let row_sum: Vec<f64> = (0..n).map(|i| (0..n).map(|j| c[i * n + j]).sum()).collect();
        let s0: f64 = row_sum.iter().sum();
        if s0 <= 0.0 {
            return Err(ToolError::Execution(
                "connectivity matrix is empty — no neighbour links were found (check method/geometry)".to_string(),
            ));
        }
        let nf = n as f64;
        let g = s0 / (nf * nf); // grand sum / n²
        let r: Vec<f64> = row_sum.iter().map(|s| s / nf).collect(); // row sum / n

        // ── 2. Double-centre: B = M C M = C_ij − r_i − r_j + g ─────────────────
        let mut b = vec![0.0f64; n * n];
        for i in 0..n {
            for j in 0..n {
                b[i * n + j] = c[i * n + j] - r[i] - r[j] + g;
            }
        }

        // ── 3. Symmetric eigendecomposition via cyclic Jacobi rotations ────────
        ctx.progress
            .info(&format!("{n} features; running Jacobi eigensolver"));
        let (eigenvalues, vectors) = jacobi_eigen(b, n);

        // ── 4. Moran's I per eigenvector; keep the autocorrelated ones ─────────
        let scale = nf / s0; // Moran's I = (n/S0)·λ
        let mut ranked: Vec<(usize, f64, f64)> = (0..n)
            .map(|k| {
                let lambda = eigenvalues[k];
                (k, lambda, scale * lambda)
            })
            .filter(|&(_, _, i)| i >= prm.min_autocorrelation)
            .collect();
        // Descending by Moran's I; index as a deterministic tie-breaker.
        ranked.sort_by(|a, b| {
            b.2.partial_cmp(&a.2)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.0.cmp(&b.0))
        });
        ranked.truncate(prm.max_components);
        let kept = ranked.len();
        ctx.progress
            .info(&format!("keeping {kept} eigenvector(s) as spatial filters"));

        // ── 5. Append kept eigenvectors as fields MEV1, MEV2, … ────────────────
        let base = layer.schema.len();
        for rank in 0..kept {
            layer.add_field(FieldDef::new(format!("MEV{}", rank + 1), FieldType::Float));
        }
        for (i, feature) in layer.features.iter_mut().enumerate() {
            // Keep positional attribute array aligned with the schema.
            if feature.attributes.len() < base {
                feature.attributes.resize(base, FieldValue::Null);
            }
            for &(k, _lambda, _i) in &ranked {
                let v = vectors[i * n + k];
                if v.is_finite() {
                    feature.attributes.push(FieldValue::Float(v));
                } else {
                    feature.attributes.push(FieldValue::Null);
                }
            }
        }

        let eig_kept: Vec<f64> = ranked.iter().map(|&(_, l, _)| l).collect();
        let morans_kept: Vec<f64> = ranked.iter().map(|&(_, _, i)| i).collect();

        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("components_kept".to_string(), json!(kept));
        outputs.insert("eigenvalues".to_string(), json!(eig_kept));
        outputs.insert("morans_i".to_string(), json!(morans_kept));
        Ok(ToolRunResult { outputs })
    }
}

// ── Symmetric eigensolver (cyclic Jacobi) ─────────────────────────────────────

/// Eigendecompose a symmetric `n×n` matrix `a` (row-major) by cyclic Jacobi
/// rotations. Returns `(eigenvalues, vectors)` where `eigenvalues[k]` is the
/// k-th eigenvalue and `vectors[i*n + k]` is component `i` of the (orthonormal)
/// eigenvector for `eigenvalues[k]` — i.e. eigenvector `k` is column `k` of the
/// returned row-major matrix.
fn jacobi_eigen(mut a: Vec<f64>, n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut v = vec![0.0f64; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    if n <= 1 {
        let eig = (0..n).map(|i| a[i * n + i]).collect();
        return (eig, v);
    }

    // Off-diagonal Frobenius norm.
    let off = |a: &[f64]| -> f64 {
        let mut s = 0.0;
        for p in 0..n {
            for q in (p + 1)..n {
                s += a[p * n + q] * a[p * n + q];
            }
        }
        s.sqrt()
    };
    let frob: f64 = a.iter().map(|x| x * x).sum::<f64>().sqrt();
    let tol = 1e-10 * frob.max(1e-300);

    for _sweep in 0..100 {
        if off(&a) <= tol {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p * n + q];
                if apq.abs() <= f64::MIN_POSITIVE {
                    continue;
                }
                let app = a[p * n + p];
                let aqq = a[q * n + q];
                // Rotation that zeros a_pq: cot(2θ) = (a_qq − a_pp)/(2 a_pq).
                let theta = (aqq - app) / (2.0 * apq);
                let t = if theta == 0.0 {
                    1.0
                } else {
                    theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt())
                };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;

                // A ← Jᵀ A J. First right-multiply (columns p,q)…
                for k in 0..n {
                    let akp = a[k * n + p];
                    let akq = a[k * n + q];
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                // …then left-multiply (rows p,q), reading the updated columns.
                for k in 0..n {
                    let apk = a[p * n + k];
                    let aqk = a[q * n + k];
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
                // V ← V J (accumulate eigenvectors as columns).
                for k in 0..n {
                    let vkp = v[k * n + p];
                    let vkq = v[k * n + q];
                    v[k * n + p] = c * vkp - s * vkq;
                    v[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }

    let eig = (0..n).map(|i| a[i * n + i]).collect();
    (eig, v)
}

// ── Neighbour builders (vendored from generate_spatial_weights_matrix) ─────────

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
    if reps.iter().any(Option::is_none) {
        return Err(ToolError::Execution(
            "every feature must have a usable point/centroid geometry for knn".to_string(),
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

/// All rings (exterior + interiors) of a polygon geometry as vertex chains
/// without the closing duplicate.
fn polygon_rings(geom: &Geometry) -> Vec<Vec<P>> {
    let ring_pts = |ring: &Ring| -> Vec<P> {
        let mut pts: Vec<P> = Vec::with_capacity(ring.len());
        for c in ring.coords() {
            let p = P { x: c.x, y: c.y };
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

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Knn,
    ContiguityEdges,
    ContiguityEdgesCorners,
}

struct Params {
    method: Method,
    number_of_neighbors: usize,
    min_autocorrelation: f64,
    max_components: usize,
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

fn parse_optional_usize(args: &ToolArgs, key: &str) -> Result<Option<usize>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(Some(n.as_u64().ok_or_else(|| {
            ToolError::Validation(format!("'{key}' must be a non-negative integer"))
        })? as usize)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => {
            Ok(Some(s.trim().parse::<usize>().map_err(|_| {
                ToolError::Validation(format!("'{key}' must be an integer"))
            })?))
        }
        Some(_) => Err(ToolError::Validation(format!("'{key}' must be an integer"))),
    }
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let method = match parse_optional_str(args, "method")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("contiguity_edges_corners") | Some("edges_corners") | Some("queen") => {
            Method::ContiguityEdgesCorners
        }
        Some("contiguity_edges") | Some("edges") | Some("rook") => Method::ContiguityEdges,
        Some("knn") | Some("k_nearest_neighbors") => Method::Knn,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "'method' must be one of contiguity_edges, contiguity_edges_corners, knn (got '{other}')"
            )))
        }
    };

    let number_of_neighbors = parse_optional_usize(args, "number_of_neighbors")?.unwrap_or(8);
    if method == Method::Knn && number_of_neighbors == 0 {
        return Err(ToolError::Validation(
            "'number_of_neighbors' must be >= 1 for knn".to_string(),
        ));
    }

    let min_autocorrelation = parse_optional_f64(args, "min_autocorrelation")?.unwrap_or(0.25);
    if !min_autocorrelation.is_finite() {
        return Err(ToolError::Validation(
            "'min_autocorrelation' must be a finite number".to_string(),
        ));
    }

    let max_components = parse_optional_usize(args, "max_components")?.unwrap_or(15);
    if max_components == 0 {
        return Err(ToolError::Validation(
            "'max_components' must be >= 1".to_string(),
        ));
    }

    Ok(Params {
        method,
        number_of_neighbors,
        min_autocorrelation,
        max_components,
    })
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

    /// Build a rows×cols grid of unit squares as a polygon layer.
    fn grid_layer(rows: usize, cols: usize) -> String {
        let mut l = Layer::new("grid")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("id", FieldType::Integer));
        let mut idx = 0i64;
        for r in 0..rows {
            for c in 0..cols {
                l.add_feature(Some(square(c as f64, r as f64, 1.0)), &[("id", idx.into())])
                    .unwrap();
                idx += 1;
            }
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SpatialEigenvectorFilteringTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn mev_columns(layer: &Layer) -> Vec<Vec<f64>> {
        let mut cols = Vec::new();
        let mut k = 1;
        loop {
            let name = format!("MEV{k}");
            let Some(fi) = layer.schema.field_index(&name) else {
                break;
            };
            let col: Vec<f64> = layer
                .iter()
                .map(|f| f.attributes[fi].as_f64().unwrap())
                .collect();
            cols.push(col);
            k += 1;
        }
        cols
    }

    fn dot(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b).map(|(x, y)| x * y).sum()
    }

    /// A regular grid has strong positive autocorrelation: at least one
    /// component is kept, Moran's I values respect the threshold and descend.
    #[test]
    fn grid_keeps_descending_components() {
        let input = grid_layer(4, 4);
        let (out, layer) = run(json!({
            "input": input,
            "method": "contiguity_edges_corners",
            "min_autocorrelation": 0.25
        }));
        let kept = out.outputs["components_kept"].as_u64().unwrap() as usize;
        assert!(kept >= 1, "expected at least one component kept");

        let morans: Vec<f64> = out.outputs["morans_i"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap())
            .collect();
        assert_eq!(morans.len(), kept);
        for w in morans.windows(2) {
            assert!(w[0] >= w[1] - 1e-12, "Moran's I not descending: {morans:?}");
        }
        for &i in &morans {
            assert!(i >= 0.25 - 1e-12, "kept component below threshold: {i}");
        }

        // One MEV field per kept component, all finite.
        let cols = mev_columns(&layer);
        assert_eq!(
            cols.len(),
            kept,
            "MEV field count must match components_kept"
        );
        for col in &cols {
            assert!(col.iter().all(|v| v.is_finite()));
        }
    }

    /// KNN path also yields autocorrelated components on a grid of point-like
    /// polygons (centroids) and appends fields.
    #[test]
    fn knn_path_keeps_components() {
        let input = grid_layer(5, 5);
        let (out, layer) = run(json!({
            "input": input,
            "method": "knn",
            "number_of_neighbors": 4,
            "min_autocorrelation": 0.1
        }));
        let kept = out.outputs["components_kept"].as_u64().unwrap() as usize;
        assert!(kept >= 1);
        assert_eq!(mev_columns(&layer).len(), kept);
    }

    /// Eigenvectors are orthonormal and double-centred: MEV1·MEV2 ≈ 0 and each
    /// column has ≈ zero mean.
    #[test]
    fn components_orthogonal_and_mean_zero() {
        let input = grid_layer(5, 5);
        let (out, layer) = run(json!({
            "input": input,
            "method": "contiguity_edges_corners",
            "min_autocorrelation": 0.1,
            "max_components": 5
        }));
        let kept = out.outputs["components_kept"].as_u64().unwrap() as usize;
        assert!(kept >= 2, "need at least two components for this test");
        let cols = mev_columns(&layer);

        // Each column is mean-zero (orthogonal to the constant vector).
        for col in &cols {
            let mean: f64 = col.iter().sum::<f64>() / col.len() as f64;
            assert!(mean.abs() < 1e-6, "column mean not ≈0: {mean}");
        }
        // Distinct columns are orthogonal; each is unit length.
        for a in 0..cols.len() {
            assert!(
                (dot(&cols[a], &cols[a]) - 1.0).abs() < 1e-6,
                "not unit length"
            );
            for b in (a + 1)..cols.len() {
                assert!(
                    dot(&cols[a], &cols[b]).abs() < 1e-6,
                    "MEV{}·MEV{} not orthogonal",
                    a + 1,
                    b + 1
                );
            }
        }
    }

    /// Standalone check of the Jacobi eigensolver on a known symmetric matrix.
    /// [[2,1],[1,2]] has eigenvalues 3 and 1 with eigenvectors (1,1)/√2 and
    /// (1,−1)/√2.
    #[test]
    fn jacobi_matches_analytic_2x2() {
        let (eig, vec) = jacobi_eigen(vec![2.0, 1.0, 1.0, 2.0], 2);
        let mut vals = eig.clone();
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((vals[0] - 1.0).abs() < 1e-9, "λ0 = {}", vals[0]);
        assert!((vals[1] - 3.0).abs() < 1e-9, "λ1 = {}", vals[1]);
        // Columns are orthonormal.
        let c0 = [vec[0], vec[2]];
        let c1 = [vec[1], vec[3]];
        assert!((dot(&c0, &c0) - 1.0).abs() < 1e-9);
        assert!((dot(&c1, &c1) - 1.0).abs() < 1e-9);
        assert!(dot(&c0, &c1).abs() < 1e-9);
    }

    /// A 3×3 symmetric matrix with analytic eigenvalues (diagonal-dominant
    /// arrowhead): verify the eigenvalues sum/product invariants.
    #[test]
    fn jacobi_diagonal_3x3() {
        // Diagonal matrix: eigenvalues are the diagonal entries.
        let (eig, _v) = jacobi_eigen(vec![5.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, -2.0], 3);
        let mut vals = eig.clone();
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((vals[0] + 2.0).abs() < 1e-9);
        assert!((vals[1] - 3.0).abs() < 1e-9);
        assert!((vals[2] - 5.0).abs() < 1e-9);
    }

    #[test]
    fn rejects_missing_input() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(SpatialEigenvectorFilteringTool.validate(&args).is_err());
    }

    #[test]
    fn rejects_bad_method() {
        let input = grid_layer(2, 2);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "method": "nonsense" })).unwrap();
        assert!(SpatialEigenvectorFilteringTool.validate(&args).is_err());
    }

    #[test]
    fn rejects_zero_max_components() {
        let input = grid_layer(2, 2);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "max_components": 0 })).unwrap();
        assert!(SpatialEigenvectorFilteringTool.validate(&args).is_err());
    }
}
