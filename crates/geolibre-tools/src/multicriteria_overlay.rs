//! GeoLibre tool: rank-based multicriteria overlay (MCDA) over a raster stack.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Multicriteria Overlay* (Spatial
//! Analyst). The bundled suitability tools are crisp linear combiners
//! (`weighted_overlay`/`weighted_sum`/`sum_overlay`) or fuzzy logic
//! (`fuzzy_overlay`); none implements the rank-based **TOPSIS** or **OWA**
//! methods that Multicriteria Overlay adds. This tool takes N co-registered
//! criterion rasters plus per-raster weights and combines them cell-wise with
//! one of four closed-form, deterministic methods:
//!
//! * **`weighted_sum`** — Σ wᵢ·nᵢ (the classic linear suitability score).
//! * **`weighted_geometric_mean`** — Π nᵢ^wᵢ (a compensatory-but-conjunctive
//!   combiner; a low score on any criterion drags the result down harder than
//!   the weighted sum does).
//! * **`owa`** — Ordered Weighted Averaging (Yager / Malczewski). The criteria
//!   are re-ordered per cell by value; a set of *order weights* then controls
//!   the degree of ORness/ANDness (how optimistic vs. pessimistic the trade-off
//!   is). With all-equal order weights OWA reduces exactly to the weighted sum.
//! * **`topsis`** — Technique for Order Preference by Similarity to Ideal
//!   Solution. Each cell's Euclidean distance to the positive and negative
//!   ideal points is measured; the output is the *closeness* ratio
//!   D⁻/(D⁺+D⁻) ∈ [0, 1] (1 = coincident with the ideal).
//!
//! Every criterion is min-max rescaled to a 0..1 benefit surface using its own
//! valid-cell range before combination, so rasters in different native units
//! are directly comparable. `from_scale`/`to_scale` linearly rescale the final
//! score (default 0..1; e.g. set to 1..10 for an ArcGIS-style evaluation scale).
//!
//! No-data propagates: a cell that is no-data (or non-finite) in *any* input is
//! no-data in the output. The result is a single F32 band with a distinct
//! no-data value.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::common::{load_input_raster, parse_optional_output, raster_like_with_data};

/// No-data value for the output score (outside every plausible score range).
const OUT_NODATA: f64 = -9999.0;

#[derive(Clone, Copy, PartialEq)]
enum Method {
    WeightedSum,
    WeightedGeometricMean,
    Owa,
    Topsis,
}

pub struct MulticriteriaOverlayTool;

impl Tool for MulticriteriaOverlayTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "multicriteria_overlay",
            display_name: "Multicriteria Overlay",
            summary: "Combine a stack of criterion rasters into one suitability surface with rank-based MCDA: weighted_sum, weighted_geometric_mean, OWA (ordered weighted averaging), or TOPSIS (closeness to the ideal point), with per-raster weights and from_scale/to_scale rescaling, like ArcGIS Multicriteria Overlay.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "inputs",
                    description: "Comma-separated list of ≥2 co-registered criterion rasters (higher value = more suitable; each is min-max rescaled to 0..1 before combination).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output raster path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "Combination method: 'weighted_sum' (default), 'weighted_geometric_mean', 'owa', 'topsis'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "weights",
                    description: "Comma-separated per-raster importance weights, one per input (default equal). Normalized to sum to 1.",
                    required: false,
                },
                ToolParamSpec {
                    name: "order_weights",
                    description: "OWA only: comma-separated order weights, one per input, applied to the criteria sorted high→low per cell (default equal, which makes OWA == weighted mean). Normalized to sum to 1.",
                    required: false,
                },
                ToolParamSpec {
                    name: "from_scale",
                    description: "Low bound of the output evaluation scale (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "to_scale",
                    description: "High bound of the output evaluation scale (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        let list = args
            .get("inputs")
            .and_then(Value::as_str)
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .ok_or_else(|| {
                ToolError::Validation("missing required parameter 'inputs'".to_string())
            })?;
        let n = split_paths(list).len();
        if n < 2 {
            return Err(ToolError::Validation(
                "'inputs' must list at least 2 criterion rasters".to_string(),
            ));
        }
        parse_method(args)?;
        // Validate weight-vector lengths eagerly where present.
        if let Some(w) = parse_weight_vec(args, "weights")? {
            if w.len() != n {
                return Err(ToolError::Validation(format!(
                    "'weights' has {} entries but there are {n} inputs",
                    w.len()
                )));
            }
        }
        if let Some(w) = parse_weight_vec(args, "order_weights")? {
            if w.len() != n {
                return Err(ToolError::Validation(format!(
                    "'order_weights' has {} entries but there are {n} inputs",
                    w.len()
                )));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let output = parse_optional_output(args, "output")?;
        let method = parse_method(args)?;

        let list = args
            .get("inputs")
            .and_then(Value::as_str)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .ok_or_else(|| {
                ToolError::Validation("missing required parameter 'inputs'".to_string())
            })?;
        let paths = split_paths(list);
        if paths.len() < 2 {
            return Err(ToolError::Validation(
                "'inputs' must list at least 2 criterion rasters".to_string(),
            ));
        }
        let n = paths.len();

        let rasters: Vec<Raster> = paths
            .iter()
            .map(|p| load_input_raster(p))
            .collect::<Result<_, _>>()?;
        let rows = rasters[0].rows;
        let cols = rasters[0].cols;
        for (i, r) in rasters.iter().enumerate() {
            if r.rows != rows || r.cols != cols {
                return Err(ToolError::Validation(format!(
                    "input {} is {}x{}, expected {rows}x{cols} — all criterion rasters must align",
                    i, r.rows, r.cols
                )));
            }
        }

        // Importance weights (normalized to sum 1); default equal.
        let weights = normalize_weights(
            parse_weight_vec(args, "weights")?.unwrap_or_else(|| vec![1.0; n]),
            n,
            "weights",
        )?;
        // OWA order weights (normalized); default equal → OWA == weighted sum.
        let order_weights = normalize_weights(
            parse_weight_vec(args, "order_weights")?.unwrap_or_else(|| vec![1.0; n]),
            n,
            "order_weights",
        )?;

        let from_scale = parse_f64(args, "from_scale")?.unwrap_or(0.0);
        let to_scale = parse_f64(args, "to_scale")?.unwrap_or(1.0);

        // Pass 1: per-criterion min/max over fully-valid cells (a cell counts
        // only when every criterion is present there, matching the no-data mask
        // used for the actual combination).
        let mut cmin = vec![f64::INFINITY; n];
        let mut cmax = vec![f64::NEG_INFINITY; n];
        let mut valid_cells = 0u64;
        for row in 0..rows as isize {
            for col in 0..cols as isize {
                if let Some(vals) = read_cell(&rasters, row, col) {
                    valid_cells += 1;
                    for i in 0..n {
                        cmin[i] = cmin[i].min(vals[i]);
                        cmax[i] = cmax[i].max(vals[i]);
                    }
                }
            }
        }
        if valid_cells == 0 {
            return Err(ToolError::Execution(
                "no cell is valid (non-nodata) across all criterion rasters".to_string(),
            ));
        }

        // Normalized extremes per criterion (0/1 when the range is non-zero,
        // 0.5/0.5 for a constant criterion). Drive the TOPSIS ideal points.
        let norm_min: Vec<f64> = (0..n)
            .map(|i| normalize(cmin[i], cmin[i], cmax[i]))
            .collect();
        let norm_max: Vec<f64> = (0..n)
            .map(|i| normalize(cmax[i], cmin[i], cmax[i]))
            .collect();
        // Positive/negative ideal points in weighted-normalized space.
        let ideal_pos: Vec<f64> = (0..n).map(|i| weights[i] * norm_max[i]).collect();
        let ideal_neg: Vec<f64> = (0..n).map(|i| weights[i] * norm_min[i]).collect();

        ctx.progress.info(&format!(
            "multicriteria overlay of {n} rasters via {}",
            method_name(method)
        ));

        // Pass 2: combine.
        let mut data = vec![OUT_NODATA; rows * cols];
        let mut smin = f64::INFINITY;
        let mut smax = f64::NEG_INFINITY;
        let mut ssum = 0.0;
        for row in 0..rows {
            for col in 0..cols {
                let Some(raw) = read_cell(&rasters, row as isize, col as isize) else {
                    continue;
                };
                let norm: Vec<f64> = (0..n)
                    .map(|i| normalize(raw[i], cmin[i], cmax[i]))
                    .collect();
                let score = match method {
                    Method::WeightedSum => weighted_sum(&norm, &weights),
                    Method::WeightedGeometricMean => weighted_geometric_mean(&norm, &weights),
                    Method::Owa => owa(&norm, &weights, &order_weights),
                    Method::Topsis => topsis(&norm, &weights, &ideal_pos, &ideal_neg),
                };
                let scaled = from_scale + score.clamp(0.0, 1.0) * (to_scale - from_scale);
                data[row * cols + col] = scaled;
                smin = smin.min(scaled);
                smax = smax.max(scaled);
                ssum += scaled;
            }
        }

        let out = raster_like_with_data(&rasters[0], data, OUT_NODATA, DataType::F32)?;
        let out_path = crate::common::write_or_store_output(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("method".to_string(), json!(method_name(method)));
        outputs.insert("input_count".to_string(), json!(n));
        outputs.insert("valid_cells".to_string(), json!(valid_cells));
        outputs.insert("weights".to_string(), json!(weights));
        outputs.insert("score_min".to_string(), json!(smin));
        outputs.insert("score_max".to_string(), json!(smax));
        outputs.insert("score_mean".to_string(), json!(ssum / valid_cells as f64));
        Ok(ToolRunResult { outputs })
    }
}

// ── Cell access ────────────────────────────────────────────────────────────────

/// Reads all criteria at one cell, returning `None` when any criterion is
/// no-data or non-finite (the shared no-data mask).
fn read_cell(rasters: &[Raster], row: isize, col: isize) -> Option<Vec<f64>> {
    let mut vals = Vec::with_capacity(rasters.len());
    for r in rasters {
        let v = r.get(0, row, col);
        if v == r.nodata || !v.is_finite() {
            return None;
        }
        vals.push(v);
    }
    Some(vals)
}

/// Min-max rescale to [0, 1]; a constant criterion (range 0) maps to the
/// neutral 0.5 so it neither helps nor hurts the score.
fn normalize(x: f64, min: f64, max: f64) -> f64 {
    if max <= min {
        0.5
    } else {
        ((x - min) / (max - min)).clamp(0.0, 1.0)
    }
}

// ── Combination methods ─────────────────────────────────────────────────────────

fn weighted_sum(norm: &[f64], w: &[f64]) -> f64 {
    norm.iter().zip(w).map(|(n, wi)| wi * n).sum()
}

fn weighted_geometric_mean(norm: &[f64], w: &[f64]) -> f64 {
    // Π nᵢ^wᵢ with Σwᵢ = 1; a zero on any criterion zeroes the product.
    let mut prod = 1.0;
    for (n, wi) in norm.iter().zip(w) {
        if *n <= 0.0 {
            return 0.0;
        }
        prod *= n.powf(*wi);
    }
    prod
}

/// Ordered Weighted Averaging (Malczewski's importance-weighted form). The
/// (importance-weight, value) pairs are sorted by value high→low; the order
/// weight for position j is combined with the reordered importance weight, then
/// renormalized. With all-equal order weights this collapses to Σwᵢnᵢ.
fn owa(norm: &[f64], w: &[f64], order: &[f64]) -> f64 {
    let n = norm.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| {
        norm[b]
            .partial_cmp(&norm[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    // combined_j = u_j * v_j, where u_j is the importance weight of the j-th
    // largest value and v_j the j-th order weight.
    let mut combined = vec![0.0; n];
    let mut denom = 0.0;
    for (j, &i) in idx.iter().enumerate() {
        combined[j] = w[i] * order[j];
        denom += combined[j];
    }
    if denom <= 0.0 {
        return 0.0;
    }
    let mut acc = 0.0;
    for (j, &i) in idx.iter().enumerate() {
        acc += (combined[j] / denom) * norm[i];
    }
    acc
}

/// TOPSIS closeness coefficient. Distances are taken in weighted-normalized
/// space to the positive/negative ideal points; the result is D⁻/(D⁺+D⁻).
fn topsis(norm: &[f64], w: &[f64], ideal_pos: &[f64], ideal_neg: &[f64]) -> f64 {
    let mut dp = 0.0;
    let mut dn = 0.0;
    for i in 0..norm.len() {
        let v = w[i] * norm[i];
        dp += (v - ideal_pos[i]).powi(2);
        dn += (v - ideal_neg[i]).powi(2);
    }
    let dp = dp.sqrt();
    let dn = dn.sqrt();
    if dp + dn <= 0.0 {
        0.5
    } else {
        dn / (dp + dn)
    }
}

// ── Parameters ──────────────────────────────────────────────────────────────────

fn split_paths(list: &str) -> Vec<&str> {
    list.split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .collect()
}

fn parse_method(args: &ToolArgs) -> Result<Method, ToolError> {
    match args
        .get("method")
        .and_then(Value::as_str)
        .map(|s| s.trim().to_lowercase())
    {
        None => Ok(Method::WeightedSum),
        Some(s) if s.is_empty() || s == "weighted_sum" || s == "weighted sum" => {
            Ok(Method::WeightedSum)
        }
        Some(s)
            if s == "weighted_geometric_mean"
                || s == "weighted geometric mean"
                || s == "geometric_mean"
                || s == "geometric" =>
        {
            Ok(Method::WeightedGeometricMean)
        }
        Some(s) if s == "owa" => Ok(Method::Owa),
        Some(s) if s == "topsis" => Ok(Method::Topsis),
        Some(other) => Err(ToolError::Validation(format!(
            "'method' must be weighted_sum|weighted_geometric_mean|owa|topsis, got '{other}'"
        ))),
    }
}

/// Parses a comma-separated numeric vector; `None` when the param is absent.
fn parse_weight_vec(args: &ToolArgs, key: &str) -> Result<Option<Vec<f64>>, ToolError> {
    let s = match args.get(key) {
        None | Some(Value::Null) => return Ok(None),
        Some(Value::String(s)) if s.trim().is_empty() => return Ok(None),
        Some(Value::String(s)) => s.clone(),
        Some(_) => {
            return Err(ToolError::Validation(format!(
                "'{key}' must be a comma-separated string of numbers"
            )))
        }
    };
    let mut out = Vec::new();
    for tok in s.split(',').map(str::trim).filter(|t| !t.is_empty()) {
        let v = tok.parse::<f64>().map_err(|_| {
            ToolError::Validation(format!("'{key}' has a non-numeric entry '{tok}'"))
        })?;
        if !v.is_finite() || v < 0.0 {
            return Err(ToolError::Validation(format!(
                "'{key}' entries must be finite and non-negative"
            )));
        }
        out.push(v);
    }
    if out.is_empty() {
        return Ok(None);
    }
    Ok(Some(out))
}

fn normalize_weights(w: Vec<f64>, n: usize, key: &str) -> Result<Vec<f64>, ToolError> {
    if w.len() != n {
        return Err(ToolError::Validation(format!(
            "'{key}' has {} entries but there are {n} inputs",
            w.len()
        )));
    }
    let sum: f64 = w.iter().sum();
    if sum <= 0.0 {
        return Err(ToolError::Validation(format!(
            "'{key}' must have a positive sum"
        )));
    }
    Ok(w.into_iter().map(|x| x / sum).collect())
}

fn parse_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(num)) => Ok(num.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("'{key}' must be a number"))),
        Some(_) => Err(ToolError::Validation(format!("'{key}' must be a number"))),
    }
}

fn method_name(m: Method) -> &'static str {
    match m {
        Method::WeightedSum => "weighted_sum",
        Method::WeightedGeometricMean => "weighted_geometric_mean",
        Method::Owa => "owa",
        Method::Topsis => "topsis",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{memory_store, CrsInfo, DataType, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn raster_of(rows: usize, cols: usize, vals: &[f64], nodata: f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, vals[row * cols + col])
                    .unwrap();
            }
        }
        let id = memory_store::put_raster(r);
        memory_store::make_raster_memory_path(&id)
    }

    fn read_all(path: &str) -> (Vec<f64>, f64) {
        let r = load_input_raster(path).unwrap();
        let mut v = Vec::new();
        for row in 0..r.rows as isize {
            for col in 0..r.cols as isize {
                v.push(r.get(0, row, col));
            }
        }
        (v, r.nodata)
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        MulticriteriaOverlayTool.run(&args, &ctx()).unwrap()
    }

    /// Weighted sum matches a hand computation on normalized criteria.
    /// A: [0,5,10] -> norm [0,0.5,1]; B: [10,5,0] -> norm [1,0.5,0].
    /// weights 0.75/0.25 -> 0.75*normA + 0.25*normB.
    #[test]
    fn weighted_sum_matches_hand_computation() {
        let a = raster_of(1, 3, &[0.0, 5.0, 10.0], -9999.0);
        let b = raster_of(1, 3, &[10.0, 5.0, 0.0], -9999.0);
        let out = run(json!({
            "inputs": format!("{a},{b}"),
            "method": "weighted_sum",
            "weights": "0.75,0.25",
        }));
        let (v, _) = read_all(out.outputs["output"].as_str().unwrap());
        assert!(
            (v[0] - 0.25).abs() < 1e-6,
            "0.75*0+0.25*1=0.25, got {}",
            v[0]
        );
        assert!((v[1] - 0.50).abs() < 1e-6, "0.75*0.5+0.25*0.5=0.5");
        assert!((v[2] - 0.75).abs() < 1e-6, "0.75*1+0.25*0=0.75");
    }

    /// OWA with all-equal order weights equals the weighted mean.
    #[test]
    fn owa_equal_order_weights_is_weighted_mean() {
        let a = raster_of(1, 3, &[0.0, 5.0, 10.0], -9999.0);
        let b = raster_of(1, 3, &[10.0, 5.0, 0.0], -9999.0);
        let ws = run(json!({
            "inputs": format!("{a},{b}"), "method": "weighted_sum", "weights": "0.75,0.25",
        }));
        let owa = run(json!({
            "inputs": format!("{a},{b}"), "method": "owa", "weights": "0.75,0.25",
            "order_weights": "0.5,0.5",
        }));
        let (vw, _) = read_all(ws.outputs["output"].as_str().unwrap());
        let (vo, _) = read_all(owa.outputs["output"].as_str().unwrap());
        for (x, y) in vw.iter().zip(&vo) {
            assert!(
                (x - y).abs() < 1e-9,
                "OWA equal order weights == weighted mean"
            );
        }
    }

    /// OWA order weights biased to the largest value = fuzzy OR (optimistic):
    /// in a cell where the criteria disagree it returns the larger norm.
    #[test]
    fn owa_or_like_takes_the_max() {
        let a = raster_of(1, 3, &[0.0, 10.0, 10.0], -9999.0); // norm [0,1,1]
        let b = raster_of(1, 3, &[0.0, 0.0, 10.0], -9999.0); //  norm [0,0,1]
                                                             // Equal importance; all order weight on the top-ranked criterion.
        let out = run(json!({
            "inputs": format!("{a},{b}"), "method": "owa",
            "weights": "0.5,0.5", "order_weights": "1,0",
        }));
        let (v, _) = read_all(out.outputs["output"].as_str().unwrap());
        // cell1 norms are {1, 0}; OR-like OWA picks the max = 1 (mean would be 0.5).
        assert!(
            (v[1] - 1.0).abs() < 1e-9,
            "OR-like OWA returns the max, got {}",
            v[1]
        );
    }

    /// TOPSIS closeness is in [0,1]; the all-max cell -> 1, the all-min -> 0.
    #[test]
    fn topsis_closeness_bounds() {
        // Two cells: cell0 = criterion maxima, cell1 = criterion minima.
        let a = raster_of(1, 2, &[10.0, 0.0], -9999.0);
        let b = raster_of(1, 2, &[20.0, 5.0], -9999.0);
        let out = run(json!({ "inputs": format!("{a},{b}"), "method": "topsis" }));
        let (v, _) = read_all(out.outputs["output"].as_str().unwrap());
        assert!(
            (v[0] - 1.0).abs() < 1e-9,
            "best cell closeness 1, got {}",
            v[0]
        );
        assert!(
            (v[1] - 0.0).abs() < 1e-9,
            "worst cell closeness 0, got {}",
            v[1]
        );
        assert!(out.outputs["method"].as_str() == Some("topsis"));
    }

    /// Weighted geometric mean zeroes when any criterion is zero, and equals the
    /// plain geometric mean of the norms under equal weights otherwise.
    #[test]
    fn weighted_geometric_mean_behaviour() {
        let a = raster_of(1, 3, &[0.0, 10.0, 20.0], -9999.0); // norm [0,0.5,1]
        let b = raster_of(1, 3, &[0.0, 10.0, 20.0], -9999.0); // norm [0,0.5,1]
        let out = run(json!({ "inputs": format!("{a},{b}"), "method": "weighted_geometric_mean" }));
        let (v, _) = read_all(out.outputs["output"].as_str().unwrap());
        // cell0: normA=0 -> product 0.
        assert!((v[0] - 0.0).abs() < 1e-9, "zero criterion zeroes GM");
        // cell1: 0.5^0.5 * 0.5^0.5 = 0.5.
        assert!(
            (v[1] - 0.5).abs() < 1e-9,
            "GM of {{0.5,0.5}} = 0.5, got {}",
            v[1]
        );
        // cell2: normA=1, normB=1 -> 1.
        assert!((v[2] - 1.0).abs() < 1e-9);
    }

    /// from_scale/to_scale linearly rescale the score.
    #[test]
    fn from_to_scale_rescales() {
        let a = raster_of(1, 3, &[0.0, 5.0, 10.0], -9999.0);
        let b = raster_of(1, 3, &[0.0, 5.0, 10.0], -9999.0);
        let out = run(json!({
            "inputs": format!("{a},{b}"), "method": "weighted_sum",
            "from_scale": 1, "to_scale": 10,
        }));
        let (v, _) = read_all(out.outputs["output"].as_str().unwrap());
        // scores 0,0.5,1 -> 1, 5.5, 10.
        assert!((v[0] - 1.0).abs() < 1e-6);
        assert!((v[1] - 5.5).abs() < 1e-6);
        assert!((v[2] - 10.0).abs() < 1e-6);
    }

    /// No-data in any input propagates to no-data out.
    #[test]
    fn propagates_nodata() {
        let a = raster_of(1, 2, &[5.0, -9999.0], -9999.0);
        let b = raster_of(1, 2, &[5.0, 5.0], -9999.0);
        let out = run(json!({ "inputs": format!("{a},{b}"), "method": "weighted_sum" }));
        let (v, nd) = read_all(out.outputs["output"].as_str().unwrap());
        assert_eq!(v[1], nd, "nodata propagated");
    }

    #[test]
    fn rejects_bad_parameters() {
        // Missing inputs.
        let args: ToolArgs = serde_json::from_value(json!({ "method": "topsis" })).unwrap();
        assert!(MulticriteriaOverlayTool.validate(&args).is_err());
        // Single input.
        let a = raster_of(1, 1, &[1.0], -9999.0);
        let args: ToolArgs = serde_json::from_value(json!({ "inputs": a.clone() })).unwrap();
        assert!(MulticriteriaOverlayTool.validate(&args).is_err());
        // Bad method.
        let args: ToolArgs =
            serde_json::from_value(json!({ "inputs": format!("{a},{a}"), "method": "bogus" }))
                .unwrap();
        assert!(MulticriteriaOverlayTool.validate(&args).is_err());
        // Weight-count mismatch.
        let args: ToolArgs =
            serde_json::from_value(json!({ "inputs": format!("{a},{a}"), "weights": "1,2,3" }))
                .unwrap();
        assert!(MulticriteriaOverlayTool.validate(&args).is_err());
    }
}
