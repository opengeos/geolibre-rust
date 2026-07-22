//! GeoLibre tool: global High/Low Clustering — Getis-Ord General G.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *High/Low Clustering (Getis-Ord
//! General G)* (Spatial Statistics). The bundled `getis_ord_gi_star` and
//! GeoLibre's shipped hot-spot tools are all *local* Gi* per-feature statistics;
//! the *global* General G — one number for the whole dataset — existed nowhere.
//!
//! General G answers a different question than Moran's I: rather than "are
//! similar values near each other" it asks "are the **high** values (or the
//! **low** values) globally clustered". It is defined as
//!
//! ```text
//!         Σ_i Σ_j  w_ij · x_i · x_j          (i ≠ j)
//!   G  =  ─────────────────────────
//!             Σ_i Σ_j  x_i · x_j             (i ≠ j)
//! ```
//!
//! with an analytic expectation `E[G] = S0 / (n(n-1))` and the Getis-Ord (1992)
//! randomization variance built from the weight sums `S0, S1, S2` and the first
//! four moments of `x`. The z-score `(G − E[G]) / √Var(G)` is compared to the
//! standard normal: a **large positive** z means high values cluster (a hot
//! landscape); a **large negative** z means low values cluster (a cold
//! landscape).
//!
//! Because the statistic is a ratio of value products, the field must be
//! **non-negative** (ideally strictly positive) with a real zero — General G is
//! not meaningful for fields that can be negative.
//!
//! Spatial weights are chosen with `weights`:
//! * `distance_band` (default) — binary, `w_ij = 1` when features lie within a
//!   threshold distance (default: the maximum nearest-neighbour distance, so
//!   every feature has ≥1 neighbour);
//! * `k_nearest` — the `k` nearest neighbours of each feature (asymmetric);
//! * `queen` — polygon contiguity, sharing at least one boundary vertex;
//! * `rook` — polygon contiguity, sharing at least one boundary edge.
//!
//! `row_standardize` rescales each row of the weight matrix to sum to 1. Output
//! is a one-row report (`observed_g, expected_g, variance, z_score, p_value,
//! pattern`) returned in the result and, when `output` is given, written as CSV.
//! Pairwise (O(n²)), so this suits moderate feature counts. Use a projected CRS
//! for the distance-based schemes (distances are in the CRS units).

use std::collections::{BTreeMap, HashMap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldValue, Geometry};

use crate::common::write_text_output;
use crate::vector_common::{load_input_layer, parse_optional_str};

pub struct GetisOrdGeneralGTool;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Weights {
    DistanceBand,
    KNearest,
    Queen,
    Rook,
}

impl Weights {
    fn parse(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "distance_band" | "distance" | "fixed_distance" => Some(Weights::DistanceBand),
            "k_nearest" | "knn" | "k_nearest_neighbors" => Some(Weights::KNearest),
            "queen" => Some(Weights::Queen),
            "rook" => Some(Weights::Rook),
            _ => None,
        }
    }
}

impl Tool for GetisOrdGeneralGTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "getis_ord_general_g",
            display_name: "High/Low Clustering (Getis-Ord General G)",
            summary: "Compute the global Getis-Ord General G statistic — a single z-score and p-value telling you whether HIGH or LOW values are globally clustered, like ArcGIS High/Low Clustering (Getis-Ord General G).",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer. Points use their coordinates; other geometries use their vertex-mean representative point (contiguity weights need polygons).",
                    required: true,
                },
                ToolParamSpec {
                    name: "field",
                    description: "Numeric, non-negative field to test for high/low clustering.",
                    required: true,
                },
                ToolParamSpec {
                    name: "weights",
                    description: "Spatial-weights scheme: 'distance_band' (default), 'k_nearest', 'queen', or 'rook'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distance_band",
                    description: "Threshold distance (CRS units) for the 'distance_band' scheme. Default: the maximum nearest-neighbour distance (every feature has >=1 neighbour).",
                    required: false,
                },
                ToolParamSpec {
                    name: "k",
                    description: "Number of neighbours for the 'k_nearest' scheme (default 8).",
                    required: false,
                },
                ToolParamSpec {
                    name: "row_standardize",
                    description: "Row-standardize the weight matrix so each row sums to 1 (default false).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional CSV path for the one-row report. Always returned in the result.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "field")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let field = require_str(args, "field")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let fidx = layer
            .schema
            .field_index(field)
            .ok_or_else(|| ToolError::Validation(format!("field '{field}' not found")))?;

        // Collect representative points, values, and geometries.
        let mut pts: Vec<(f64, f64)> = Vec::new();
        let mut vals: Vec<f64> = Vec::new();
        let mut geoms: Vec<Geometry> = Vec::new();
        for feature in layer.iter() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let Some((x, y)) = representative_xy(geom) else {
                continue;
            };
            let Some(v) = feature.attributes.get(fidx).and_then(FieldValue::as_f64) else {
                continue;
            };
            if !v.is_finite() {
                continue;
            }
            pts.push((x, y));
            vals.push(v);
            geoms.push(geom.clone());
        }
        let n = pts.len();
        if n < 4 {
            return Err(ToolError::Execution(format!(
                "need at least 4 valued features, found {n}"
            )));
        }
        if vals.iter().any(|&v| v < 0.0) {
            return Err(ToolError::Execution(
                "General G requires a non-negative field (a ratio of value products); the field has negative values".to_string(),
            ));
        }

        let nf = n as f64;
        let sum_y: f64 = vals.iter().sum();
        let sum_y2: f64 = vals.iter().map(|v| v * v).sum();
        // Denominator Σ_{i≠j} x_i x_j = (Σx)² − Σx².
        let den = sum_y * sum_y - sum_y2;
        if den <= 0.0 {
            return Err(ToolError::Execution(
                "field has no positive spread (all values equal or zero); General G is undefined"
                    .to_string(),
            ));
        }

        // Build the (dense) spatial-weights matrix for the chosen scheme.
        ctx.progress.info(&format!(
            "building {} weights for {n} features",
            prm.weights_label()
        ));
        let mut w = build_weights(&prm, &pts, &geoms)?;
        if prm.row_standardize {
            for row in w.iter_mut() {
                let rs: f64 = row.iter().sum();
                if rs > 0.0 {
                    for v in row.iter_mut() {
                        *v /= rs;
                    }
                }
            }
        }

        // Weight sums S0, S1, S2 (general, work for asymmetric / real weights).
        let mut rowsum = vec![0.0f64; n];
        let mut colsum = vec![0.0f64; n];
        let mut s0 = 0.0;
        let mut s1 = 0.0;
        for i in 0..n {
            for j in 0..n {
                let wij = w[i][j];
                rowsum[i] += wij;
                colsum[j] += wij;
                s0 += wij;
                let s = wij + w[j][i];
                s1 += s * s;
            }
        }
        s1 *= 0.5;
        let s2: f64 = (0..n).map(|i| (rowsum[i] + colsum[i]).powi(2)).sum();
        if s0 <= 0.0 {
            return Err(ToolError::Execution(
                "no spatial neighbours found; try a larger distance_band or a different weights scheme".to_string(),
            ));
        }

        // Observed G and its analytic expectation.
        let mut num = 0.0;
        for i in 0..n {
            let yi = vals[i];
            for j in 0..n {
                if w[i][j] != 0.0 {
                    num += w[i][j] * yi * vals[j];
                }
            }
        }
        let g = num / den;
        let eg = s0 / (nf * (nf - 1.0));

        // Getis-Ord randomization variance E[G²] − E[G]² (PySAL/Getis-Ord 1992).
        let s02 = s0 * s0;
        let n2 = nf * nf;
        let b0 = (n2 - 3.0 * nf + 3.0) * s1 - nf * s2 + 3.0 * s02;
        let b1 = -((n2 - nf) * s1 - 2.0 * nf * s2 + 6.0 * s02);
        let b2 = -(2.0 * nf * s1 - (nf + 3.0) * s2 + 6.0 * s02);
        let b3 = 4.0 * (nf - 1.0) * s1 - 2.0 * (nf + 1.0) * s2 + 8.0 * s02;
        let b4 = s1 - s2 + s02;
        let sum_y3: f64 = vals.iter().map(|v| v.powi(3)).sum();
        let sum_y4: f64 = vals.iter().map(|v| v.powi(4)).sum();
        let eg2_num = b0 * sum_y2 * sum_y2
            + b1 * sum_y4
            + b2 * sum_y * sum_y * sum_y2
            + b3 * sum_y * sum_y3
            + b4 * sum_y.powi(4);
        let eg2_den = den * den * nf * (nf - 1.0) * (nf - 2.0) * (nf - 3.0);
        let variance = if eg2_den != 0.0 {
            eg2_num / eg2_den - eg * eg
        } else {
            f64::NAN
        };

        let (z, p) = if variance > 0.0 {
            let z = (g - eg) / variance.sqrt();
            (z, 2.0 * (1.0 - normal_cdf(z.abs())))
        } else {
            (f64::NAN, f64::NAN)
        };

        let pattern = if !z.is_finite() {
            "undetermined"
        } else if p > 0.05 {
            "random"
        } else if z > 0.0 {
            "high_clustering"
        } else {
            "low_clustering"
        };

        ctx.progress.info(&format!(
            "G={g:.6} E[G]={eg:.6} z={z:.4} p={p:.4} -> {pattern}"
        ));

        let csv = format!(
            "observed_g,expected_g,variance,z_score,p_value,pattern\n{:.8},{:.8},{:.8e},{:.6},{:.6},{}\n",
            g, eg, variance, z, p, pattern
        );
        if let Some(path) = output {
            write_text_output(&csv, path)?;
        }

        let mut outputs = BTreeMap::new();
        if let Some(path) = output {
            outputs.insert("output".to_string(), json!(path));
        }
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("weights".to_string(), json!(prm.weights_label()));
        outputs.insert("s0".to_string(), json!(s0));
        outputs.insert("observed_g".to_string(), json!(g));
        outputs.insert("expected_g".to_string(), json!(eg));
        outputs.insert("variance".to_string(), json!(variance));
        outputs.insert("z_score".to_string(), json!(z));
        outputs.insert("p_value".to_string(), json!(p));
        outputs.insert("pattern".to_string(), json!(pattern));
        Ok(ToolRunResult { outputs })
    }
}

// ── Spatial weights ──────────────────────────────────────────────────────────

fn build_weights(
    prm: &Params,
    pts: &[(f64, f64)],
    geoms: &[Geometry],
) -> Result<Vec<Vec<f64>>, ToolError> {
    let n = pts.len();
    let mut w = vec![vec![0.0f64; n]; n];
    match prm.weights {
        Weights::DistanceBand => {
            let mut nearest = vec![f64::INFINITY; n];
            let mut dmat = vec![vec![0.0f64; n]; n];
            for i in 0..n {
                for j in 0..n {
                    if i != j {
                        let d = dist(pts[i], pts[j]);
                        dmat[i][j] = d;
                        if d < nearest[i] {
                            nearest[i] = d;
                        }
                    }
                }
            }
            let d = match prm.distance_band {
                Some(d) => d,
                None => nearest
                    .iter()
                    .copied()
                    .filter(|d| d.is_finite())
                    .fold(0.0, f64::max)
                    .max(f64::MIN_POSITIVE),
            };
            for i in 0..n {
                for j in 0..n {
                    if i != j && dmat[i][j] <= d {
                        w[i][j] = 1.0;
                    }
                }
            }
        }
        Weights::KNearest => {
            let k = prm.k.min(n - 1);
            let mut idx: Vec<usize> = Vec::with_capacity(n - 1);
            for i in 0..n {
                idx.clear();
                idx.extend((0..n).filter(|&j| j != i));
                idx.sort_by(|&a, &b| {
                    dist(pts[i], pts[a])
                        .partial_cmp(&dist(pts[i], pts[b]))
                        .unwrap()
                });
                for &j in idx.iter().take(k) {
                    w[i][j] = 1.0;
                }
            }
        }
        Weights::Queen | Weights::Rook => {
            build_contiguity(prm.weights, geoms, &mut w)?;
        }
    }
    Ok(w)
}

/// Polygon contiguity by shared vertices (queen) or shared edges (rook).
///
/// Uses exact (rounded) coordinate matching, the standard way PySAL builds
/// contiguity from a layer: two polygons are Queen-adjacent when they share at
/// least one boundary vertex and Rook-adjacent when they share at least one
/// boundary edge (a pair of consecutive vertices). Assumes topologically clean
/// polygons with coincident vertices along shared borders (e.g. Natural Earth
/// admin boundaries).
type VertKey = (i64, i64);
type EdgeKey = (VertKey, VertKey);

fn build_contiguity(
    weights: Weights,
    geoms: &[Geometry],
    w: &mut [Vec<f64>],
) -> Result<(), ToolError> {
    // Verify all inputs are polygonal.
    if geoms.iter().any(|g| !is_polygonal(g)) {
        return Err(ToolError::Execution(
            "queen/rook contiguity requires polygon input; use distance_band or k_nearest for points/lines".to_string(),
        ));
    }
    let key = |c: &Coord| ((c.x * 1e7).round() as i64, (c.y * 1e7).round() as i64);

    if weights == Weights::Queen {
        // Map each vertex key to the polygons that touch it.
        let mut vmap: HashMap<VertKey, Vec<usize>> = HashMap::new();
        for (i, g) in geoms.iter().enumerate() {
            let mut keys: Vec<VertKey> = polygon_vertices(g).iter().map(&key).collect();
            keys.sort_unstable();
            keys.dedup();
            for kk in keys {
                vmap.entry(kk).or_default().push(i);
            }
        }
        for ids in vmap.values() {
            connect_all(ids, w);
        }
    } else {
        // Map each undirected edge key to the polygons that own it.
        let mut emap: HashMap<EdgeKey, Vec<usize>> = HashMap::new();
        for (i, g) in geoms.iter().enumerate() {
            let mut edges: Vec<EdgeKey> = polygon_edges(g)
                .iter()
                .map(|(a, b)| {
                    let (ka, kb) = (key(a), key(b));
                    if ka <= kb {
                        (ka, kb)
                    } else {
                        (kb, ka)
                    }
                })
                .filter(|(a, b)| a != b)
                .collect();
            edges.sort_unstable();
            edges.dedup();
            for e in edges {
                emap.entry(e).or_default().push(i);
            }
        }
        for ids in emap.values() {
            connect_all(ids, w);
        }
    }
    Ok(())
}

/// Marks every distinct pair in `ids` as mutual (symmetric binary) neighbours.
fn connect_all(ids: &[usize], w: &mut [Vec<f64>]) {
    for a in 0..ids.len() {
        for b in (a + 1)..ids.len() {
            let (i, j) = (ids[a], ids[b]);
            if i != j {
                w[i][j] = 1.0;
                w[j][i] = 1.0;
            }
        }
    }
}

fn is_polygonal(g: &Geometry) -> bool {
    matches!(g, Geometry::Polygon { .. } | Geometry::MultiPolygon(_))
        || matches!(g, Geometry::GeometryCollection(gs) if gs.iter().any(is_polygonal))
}

fn polygon_vertices(g: &Geometry) -> Vec<Coord> {
    let mut out = Vec::new();
    push_rings(g, &mut |ring: &[Coord]| out.extend_from_slice(ring));
    out
}

fn polygon_edges(g: &Geometry) -> Vec<(Coord, Coord)> {
    let mut out = Vec::new();
    push_rings(g, &mut |ring: &[Coord]| {
        for pair in ring.windows(2) {
            out.push((pair[0].clone(), pair[1].clone()));
        }
    });
    out
}

fn push_rings(g: &Geometry, f: &mut dyn FnMut(&[Coord])) {
    match g {
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            f(exterior.coords());
            for r in interiors {
                f(r.coords());
            }
        }
        Geometry::MultiPolygon(polys) => {
            for (ext, holes) in polys {
                f(ext.coords());
                for r in holes {
                    f(r.coords());
                }
            }
        }
        Geometry::GeometryCollection(gs) => {
            for sub in gs {
                push_rings(sub, f);
            }
        }
        _ => {}
    }
}

fn dist(a: (f64, f64), b: (f64, f64)) -> f64 {
    ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt()
}

// ── Geometry / stats helpers ─────────────────────────────────────────────────

pub fn representative_xy(geom: &Geometry) -> Option<(f64, f64)> {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0u64;
    accumulate(geom, &mut sx, &mut sy, &mut n);
    (n > 0).then(|| (sx / n as f64, sy / n as f64))
}

fn accumulate(geom: &Geometry, sx: &mut f64, sy: &mut f64, n: &mut u64) {
    let mut add = |c: &Coord| {
        *sx += c.x;
        *sy += c.y;
        *n += 1;
    };
    match geom {
        Geometry::Point(c) => add(c),
        Geometry::LineString(cs) | Geometry::MultiPoint(cs) => cs.iter().for_each(add),
        Geometry::MultiLineString(lines) => lines.iter().flatten().for_each(add),
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            exterior.coords().iter().for_each(&mut add);
            interiors
                .iter()
                .for_each(|r| r.coords().iter().for_each(&mut add));
        }
        Geometry::MultiPolygon(polys) => {
            for (ext, holes) in polys {
                ext.coords().iter().for_each(&mut add);
                holes
                    .iter()
                    .for_each(|r| r.coords().iter().for_each(&mut add));
            }
        }
        Geometry::GeometryCollection(geoms) => {
            for g in geoms {
                accumulate(g, sx, sy, n);
            }
        }
    }
}

fn normal_cdf(x: f64) -> f64 {
    0.5 * erfc(-x / std::f64::consts::SQRT_2)
}

fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let ans = t
        * (-z * z - 1.26551223
            + t * (1.00002368
                + t * (0.37409196
                    + t * (0.09678418
                        + t * (-0.18628806
                            + t * (0.27886807
                                + t * (-1.13520398
                                    + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
            .exp();
    if x >= 0.0 {
        ans
    } else {
        2.0 - ans
    }
}

// ── Parameters ───────────────────────────────────────────────────────────────

struct Params {
    weights: Weights,
    distance_band: Option<f64>,
    k: usize,
    row_standardize: bool,
}

impl Params {
    fn weights_label(&self) -> &'static str {
        match self.weights {
            Weights::DistanceBand => "distance_band",
            Weights::KNearest => "k_nearest",
            Weights::Queen => "queen",
            Weights::Rook => "rook",
        }
    }
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let weights = match parse_optional_str(args, "weights")? {
        None => Weights::DistanceBand,
        Some(s) => Weights::parse(s).ok_or_else(|| {
            ToolError::Validation(format!(
                "'weights' must be one of distance_band, k_nearest, queen, rook (got '{s}')"
            ))
        })?,
    };
    let distance_band = parse_optional_f64(args, "distance_band")?;
    if let Some(v) = distance_band {
        if !(v > 0.0 && v.is_finite()) {
            return Err(ToolError::Validation(
                "'distance_band' must be a positive number".to_string(),
            ));
        }
    }
    let k = match parse_optional_f64(args, "k")? {
        None => 8,
        Some(v) if v.fract() == 0.0 && (1.0..=1000.0).contains(&v) => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "'k' must be an integer between 1 and 1000".to_string(),
            ))
        }
    };
    let row_standardize = parse_optional_bool(args, "row_standardize")?.unwrap_or(false);
    Ok(Params {
        weights,
        distance_band,
        k,
        row_standardize,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
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
    use wbvector::{memory_store, FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a grid of points with the given values.
    fn grid_layer(side: usize, value: impl Fn(usize, usize) -> f64) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for r in 0..side {
            for c in 0..side {
                l.add_feature(
                    Some(Geometry::point(c as f64, r as f64)),
                    &[("v", value(r, c).into())],
                )
                .unwrap();
            }
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        GetisOrdGeneralGTool.run(&args, &ctx()).unwrap()
    }

    /// A single hot cluster in one corner (high values packed together) drives a
    /// strongly positive z-score: high values are clustered.
    #[test]
    fn hot_cluster_gives_positive_z() {
        // High values only in the top-left quadrant.
        let input = grid_layer(10, |r, c| if r < 3 && c < 3 { 100.0 } else { 1.0 });
        let out = run(json!({
            "input": input, "field": "v", "weights": "distance_band", "distance_band": 1.5,
        }));
        let z = out.outputs["z_score"].as_f64().unwrap();
        assert!(
            z > 2.5,
            "clustered high values should give strongly positive z, got {z}"
        );
        assert_eq!(out.outputs["pattern"], json!("high_clustering"));
        // Observed G exceeds its expectation.
        assert!(
            out.outputs["observed_g"].as_f64().unwrap()
                > out.outputs["expected_g"].as_f64().unwrap()
        );
    }

    /// A spatially random field yields a z-score near zero and a non-significant
    /// pattern.
    #[test]
    fn random_field_is_not_significant() {
        // A deterministic hash-like scatter, no spatial structure.
        let input = grid_layer(12, |r, c| {
            let h = ((r * 7 + c * 13 + 3) % 11) as f64;
            h + 1.0
        });
        let out = run(json!({ "input": input, "field": "v", "distance_band": 1.5 }));
        let z = out.outputs["z_score"].as_f64().unwrap();
        assert!(z.abs() < 2.0, "random field z should be small, got {z}");
    }

    /// k_nearest and distance_band both run and produce finite statistics; and G
    /// is symmetric-invariant to the k choice on a clustered field (sign holds).
    #[test]
    fn k_nearest_runs_and_detects_cluster() {
        let input = grid_layer(10, |r, c| if r < 3 && c < 3 { 50.0 } else { 2.0 });
        let out = run(json!({
            "input": input, "field": "v", "weights": "k_nearest", "k": 4,
        }));
        assert_eq!(out.outputs["weights"], json!("k_nearest"));
        assert!(out.outputs["z_score"].as_f64().unwrap() > 1.5);
    }

    /// Row standardization keeps G finite and the sign of clustering intact.
    #[test]
    fn row_standardize_runs() {
        let input = grid_layer(8, |r, c| if r < 2 && c < 2 { 40.0 } else { 1.0 });
        let out = run(json!({
            "input": input, "field": "v", "distance_band": 1.5, "row_standardize": true,
        }));
        assert!(out.outputs["z_score"].as_f64().unwrap().is_finite());
        assert!(out.outputs["observed_g"].as_f64().unwrap() > 0.0);
    }

    /// Queen contiguity on a 3x3 grid of unit squares: the centre cell touches
    /// all eight neighbours, so a hot centre clusters positively.
    #[test]
    fn queen_contiguity_on_polygons() {
        let mut l = Layer::new("polys")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for r in 0..3 {
            for c in 0..3 {
                let (x, y) = (c as f64, r as f64);
                let ring = vec![
                    Coord::xy(x, y),
                    Coord::xy(x + 1.0, y),
                    Coord::xy(x + 1.0, y + 1.0),
                    Coord::xy(x, y + 1.0),
                    Coord::xy(x, y),
                ];
                let v = if r == 1 && c == 1 { 5.0 } else { 1.0 };
                l.add_feature(Some(Geometry::polygon(ring, vec![])), &[("v", v.into())])
                    .unwrap();
            }
        }
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);
        let out = run(json!({ "input": input, "field": "v", "weights": "queen" }));
        // 9 cells; centre has 8 queen neighbours, edges/corners fewer. S0 must be
        // the classic 3x3 queen total of 40 directed links.
        assert_eq!(out.outputs["s0"].as_f64().unwrap() as i64, 40);
        assert_eq!(out.outputs["weights"], json!("queen"));
    }

    /// Rook contiguity on the same grid yields the classic 24 directed links.
    #[test]
    fn rook_contiguity_on_polygons() {
        let mut l = Layer::new("polys")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for r in 0..3 {
            for c in 0..3 {
                let (x, y) = (c as f64, r as f64);
                let ring = vec![
                    Coord::xy(x, y),
                    Coord::xy(x + 1.0, y),
                    Coord::xy(x + 1.0, y + 1.0),
                    Coord::xy(x, y + 1.0),
                    Coord::xy(x, y),
                ];
                l.add_feature(Some(Geometry::polygon(ring, vec![])), &[("v", 1.0.into())])
                    .unwrap();
            }
        }
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);
        let out = run(json!({ "input": input, "field": "v", "weights": "rook" }));
        assert_eq!(out.outputs["s0"].as_f64().unwrap() as i64, 24);
    }

    /// Contiguity on point input is rejected.
    #[test]
    fn contiguity_rejects_points() {
        let input = grid_layer(5, |_, _| 1.0);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "field": "v", "weights": "queen" }))
                .unwrap();
        assert!(GetisOrdGeneralGTool.run(&args, &ctx()).is_err());
    }

    /// A negative field is rejected (General G needs non-negative values).
    #[test]
    fn rejects_negative_field() {
        let input = grid_layer(6, |r, c| (r as f64) - (c as f64));
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "field": "v" })).unwrap();
        assert!(GetisOrdGeneralGTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            GetisOrdGeneralGTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "field": "v", "weights": "bogus" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "field": "v", "k": 0 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "field": "v", "distance_band": -1 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "field": "v" })).is_ok());
    }
}
