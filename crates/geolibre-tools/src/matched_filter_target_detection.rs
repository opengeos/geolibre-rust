//! GeoLibre tool: statistical matched-filter target detection (CEM / ACE).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Detect Target Using Spectra*
//! (Image Analyst, CEM/ACE mode). Where `detect_image_anomalies` is the
//! *unsupervised* RX detector (no target signature) and the bundled
//! `spectral_angle_mapper`/`spectral_library_matching` tools score the *angle*
//! between a pixel and a reference spectrum (ignoring the scene background),
//! this tool scores every pixel against a **known target spectrum** while using
//! the scene statistics to suppress the background.
//!
//! Two classic detectors are offered:
//!
//! * `cem` — **Constrained Energy Minimization**. Builds the matched filter
//!   `w = R⁻¹ d / (dᵀ R⁻¹ d)` from the scene autocorrelation matrix
//!   `R = (1/N) Σ xₖ xₖᵀ` and reports `CEM(x) = wᵀ x = (dᵀ R⁻¹ x)/(dᵀ R⁻¹ d)`.
//!   A pixel whose spectrum equals the target scores exactly `1.0`; background
//!   is driven toward `0`.
//! * `ace` — **Adaptive Coherence Estimator**. Mean-centres both the target and
//!   each pixel against the scene mean `μ` and normalizes by the pixel energy:
//!   `ACE(x) = (δᵀ Σ⁻¹ ξ)² / [(δᵀ Σ⁻¹ δ)(ξᵀ Σ⁻¹ ξ)]` with `δ = d−μ`,
//!   `ξ = x−μ`. Output is a bounded coherence in `[0, 1]` (a CFAR detector,
//!   invariant to per-pixel scaling), with the exact-target pixel at `1.0`.
//!
//! The `R`/`Σ` matrix is `bands × bands` (small) and inverted once with the same
//! ridge-stabilized Gauss–Jordan solve `detect_image_anomalies` uses. Output is
//! a single-band F32 detection-score raster; with `threshold` (a percentile in
//! 0..100) a second binary mask marks the top scores.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};

const OUT_NODATA: f64 = -1.0;
/// Ridge added to the R/Σ diagonal (× its mean) for invertibility.
const RIDGE: f64 = 1e-6;

#[derive(Clone, Copy, PartialEq)]
enum Method {
    Cem,
    Ace,
}

pub struct MatchedFilterTargetDetectionTool;

impl Tool for MatchedFilterTargetDetectionTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "matched_filter_target_detection",
            display_name: "Matched Filter Target Detection (CEM/ACE)",
            summary: "Score every pixel of a multiband image against a known target spectrum using the scene statistics to suppress the background — Constrained Energy Minimization (CEM) or the Adaptive Coherence Estimator (ACE) — with an optional percentile threshold mask, like ArcGIS Detect Target Using Spectra.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input multiband raster (2+ bands).",
                    required: true,
                },
                ToolParamSpec {
                    name: "target_spectrum",
                    description: "Target spectrum: one value per band, as a JSON array or a comma/space-separated list (length must equal the band count).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output detection-score raster (F32). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'cem' (Constrained Energy Minimization, default) or 'ace' (Adaptive Coherence Estimator, scores in 0..1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "threshold",
                    description: "Optional percentile (0-100) of the score distribution above which a binary detection mask is written.",
                    required: false,
                },
                ToolParamSpec {
                    name: "mask_output",
                    description: "Output path for the binary detection mask (requires 'threshold'). If omitted with a threshold, stored in memory.",
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
        if parse_target_spectrum(args)?.is_none() {
            return Err(ToolError::Validation(
                "missing required parameter 'target_spectrum'".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args.get("input").and_then(Value::as_str).unwrap();
        let output = parse_optional_output(args, "output")?;
        let prm = parse_params(args)?;
        let target = parse_target_spectrum(args)?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'target_spectrum'".to_string())
        })?;

        let raster = load_input_raster(input)?;
        let b = raster.bands;
        if b < 2 {
            return Err(ToolError::Validation(format!(
                "matched-filter detection needs a multiband image; input has {b} band(s)"
            )));
        }
        if target.len() != b {
            return Err(ToolError::Validation(format!(
                "'target_spectrum' has {} values but the image has {b} band(s)",
                target.len()
            )));
        }
        let rows = raster.rows;
        let cols = raster.cols;
        let n = rows * cols;
        let nodata = raster.nodata;

        // Materialize the band stack and a per-pixel validity mask.
        let mut stack: Vec<Vec<f64>> = vec![vec![0.0; n]; b];
        let mut valid = vec![true; n];
        for (band, stack_band) in stack.iter_mut().enumerate() {
            for row in 0..rows {
                for col in 0..cols {
                    let v = raster.get(band as isize, row as isize, col as isize);
                    let i = row * cols + col;
                    if v == nodata || !v.is_finite() {
                        valid[i] = false;
                    }
                    stack_band[i] = v;
                }
            }
        }

        ctx.progress.info(&format!(
            "matched-filter target detection ({b}-band, {} method)",
            method_name(prm.method)
        ));

        let scores = match prm.method {
            Method::Cem => cem(&stack, &valid, &target, b, n)?,
            Method::Ace => ace(&stack, &valid, &target, b, n)?,
        };

        let data: Vec<f64> = (0..n)
            .map(|i| if valid[i] { scores[i] } else { OUT_NODATA })
            .collect();
        let out = raster_like_with_data(&raster, data, OUT_NODATA, DataType::F32)?;
        let out_path = write_or_store_output(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("bands".to_string(), json!(b));
        outputs.insert("method".to_string(), json!(method_name(prm.method)));
        let valid_count = valid.iter().filter(|v| **v).count();
        outputs.insert("valid_pixels".to_string(), json!(valid_count));
        let max_score = (0..n)
            .filter(|&i| valid[i])
            .map(|i| scores[i])
            .fold(f64::NEG_INFINITY, f64::max);
        if max_score.is_finite() {
            outputs.insert("max_score".to_string(), json!(max_score));
        }

        // Optional percentile threshold mask.
        if let Some(pct) = prm.threshold {
            let mut vals: Vec<f64> = (0..n).filter(|&i| valid[i]).map(|i| scores[i]).collect();
            vals.sort_by(f64::total_cmp);
            let cut = if vals.is_empty() {
                0.0
            } else {
                let idx = ((pct / 100.0) * (vals.len() - 1) as f64).round() as usize;
                vals[idx.min(vals.len() - 1)]
            };
            let mask: Vec<f64> = (0..n)
                .map(|i| {
                    if !valid[i] {
                        OUT_NODATA
                    } else if scores[i] >= cut {
                        1.0
                    } else {
                        0.0
                    }
                })
                .collect();
            let mask_raster = raster_like_with_data(&raster, mask, OUT_NODATA, DataType::F32)?;
            let mask_path = write_or_store_output(mask_raster, prm.mask_output.as_deref())?;
            let detections = (0..n).filter(|&i| valid[i] && scores[i] >= cut).count();
            outputs.insert("mask_output".to_string(), json!(mask_path));
            outputs.insert("threshold_score".to_string(), json!(cut));
            outputs.insert("detection_pixels".to_string(), json!(detections));
        }

        Ok(ToolRunResult { outputs })
    }
}

// ── Detectors ───────────────────────────────────────────────────────────────

/// Constrained Energy Minimization.
///
/// Uses the scene autocorrelation matrix `R = (1/N) Σ xₖ xₖᵀ` (no mean
/// subtraction) and reports `CEM(x) = (dᵀ R⁻¹ x)/(dᵀ R⁻¹ d)`. The exact-target
/// pixel scores `1.0`.
fn cem(
    stack: &[Vec<f64>],
    valid: &[bool],
    target: &[f64],
    b: usize,
    n: usize,
) -> Result<Vec<f64>, ToolError> {
    let r = autocorrelation(stack, valid, b, n);
    let inv = invert(&r, b)
        .ok_or_else(|| ToolError::Execution("autocorrelation matrix is singular".to_string()))?;
    // g = R⁻¹ d  and  denom = dᵀ R⁻¹ d.
    let g = mat_vec(&inv, target, b);
    let denom = dot(target, &g, b);
    if denom.abs() < 1e-30 {
        return Err(ToolError::Execution(
            "degenerate target spectrum (dᵀR⁻¹d ≈ 0)".to_string(),
        ));
    }
    let mut scores = vec![0.0; n];
    let mut x = vec![0.0; b];
    for (i, score) in scores.iter_mut().enumerate() {
        if !valid[i] {
            continue;
        }
        for (k, xk) in x.iter_mut().enumerate() {
            *xk = stack[k][i];
        }
        // CEM(x) = (R⁻¹ d)ᵀ x / denom  =  gᵀ x / denom.
        *score = dot(&g, &x, b) / denom;
    }
    Ok(scores)
}

/// Adaptive Coherence Estimator.
///
/// Mean-centres against the scene mean `μ` and uses the covariance `Σ`:
/// `ACE(x) = (δᵀ Σ⁻¹ ξ)² / [(δᵀ Σ⁻¹ δ)(ξᵀ Σ⁻¹ ξ)]` with `δ = d−μ`, `ξ = x−μ`.
/// Bounded in `[0, 1]`; the exact-target pixel scores `1.0`.
fn ace(
    stack: &[Vec<f64>],
    valid: &[bool],
    target: &[f64],
    b: usize,
    n: usize,
) -> Result<Vec<f64>, ToolError> {
    let (mean, cov) = mean_cov(stack, valid, b, n);
    let inv = invert(&cov, b)
        .ok_or_else(|| ToolError::Execution("covariance matrix is singular".to_string()))?;
    let delta: Vec<f64> = (0..b).map(|k| target[k] - mean[k]).collect();
    let ginv_delta = mat_vec(&inv, &delta, b); // Σ⁻¹ δ
    let dt = dot(&delta, &ginv_delta, b); // δᵀ Σ⁻¹ δ
    if dt.abs() < 1e-30 {
        return Err(ToolError::Execution(
            "degenerate target spectrum (δᵀΣ⁻¹δ ≈ 0)".to_string(),
        ));
    }
    let mut scores = vec![0.0; n];
    let mut xi = vec![0.0; b];
    for (i, score) in scores.iter_mut().enumerate() {
        if !valid[i] {
            continue;
        }
        for (k, x) in xi.iter_mut().enumerate() {
            *x = stack[k][i] - mean[k];
        }
        let num = dot(&ginv_delta, &xi, b); // δᵀ Σ⁻¹ ξ
        let xx = mahalanobis(&xi, &inv, b); // ξᵀ Σ⁻¹ ξ
        *score = if xx > 1e-30 {
            (num * num / (dt * xx)).clamp(0.0, 1.0)
        } else {
            0.0
        };
    }
    Ok(scores)
}

/// Scene autocorrelation matrix `R = (1/N) Σ xₖ xₖᵀ`, ridge-stabilized.
fn autocorrelation(stack: &[Vec<f64>], valid: &[bool], b: usize, n: usize) -> Vec<f64> {
    let mut r = vec![0.0; b * b];
    let mut count = 0.0_f64;
    for i in 0..n {
        if !valid[i] {
            continue;
        }
        for k in 0..b {
            let xk = stack[k][i];
            for l in k..b {
                r[k * b + l] += xk * stack[l][i];
            }
        }
        count += 1.0;
    }
    let denom = count.max(1.0);
    for k in 0..b {
        for l in k..b {
            let v = r[k * b + l] / denom;
            r[k * b + l] = v;
            r[l * b + k] = v;
        }
    }
    add_ridge(&mut r, b);
    r
}

/// Scene mean vector and (population) covariance over the valid pixels,
/// ridge-stabilized.
fn mean_cov(stack: &[Vec<f64>], valid: &[bool], b: usize, n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut mean = vec![0.0; b];
    let mut count = 0.0;
    for i in 0..n {
        if !valid[i] {
            continue;
        }
        for (k, m) in mean.iter_mut().enumerate() {
            *m += stack[k][i];
        }
        count += 1.0;
    }
    if count > 0.0 {
        for m in mean.iter_mut() {
            *m /= count;
        }
    }
    let mut cov = vec![0.0; b * b];
    for i in 0..n {
        if !valid[i] {
            continue;
        }
        for k in 0..b {
            let dk = stack[k][i] - mean[k];
            for l in k..b {
                cov[k * b + l] += dk * (stack[l][i] - mean[l]);
            }
        }
    }
    let denom = count.max(1.0);
    for k in 0..b {
        for l in k..b {
            let v = cov[k * b + l] / denom;
            cov[k * b + l] = v;
            cov[l * b + k] = v;
        }
    }
    add_ridge(&mut cov, b);
    (mean, cov)
}

/// Adds a small ridge (× the mean diagonal) to the diagonal for invertibility.
fn add_ridge(m: &mut [f64], b: usize) {
    let diag_mean: f64 = (0..b).map(|k| m[k * b + k]).sum::<f64>() / b as f64;
    let ridge = RIDGE * diag_mean.abs().max(1e-12);
    for k in 0..b {
        m[k * b + k] += ridge;
    }
}

fn dot(a: &[f64], c: &[f64], b: usize) -> f64 {
    (0..b).map(|k| a[k] * c[k]).sum()
}

fn mat_vec(m: &[f64], v: &[f64], b: usize) -> Vec<f64> {
    (0..b)
        .map(|k| (0..b).map(|l| m[k * b + l] * v[l]).sum())
        .collect()
}

fn mahalanobis(x: &[f64], inv: &[f64], b: usize) -> f64 {
    let mut s = 0.0;
    for k in 0..b {
        let mut row = 0.0;
        for l in 0..b {
            row += inv[k * b + l] * x[l];
        }
        s += x[k] * row;
    }
    s.max(0.0)
}

/// Gauss–Jordan inverse of a b×b matrix (row-major); None if singular.
fn invert(m: &[f64], b: usize) -> Option<Vec<f64>> {
    let mut a = m.to_vec();
    let mut inv = vec![0.0; b * b];
    for k in 0..b {
        inv[k * b + k] = 1.0;
    }
    for col in 0..b {
        let mut piv = col;
        let mut best = a[col * b + col].abs();
        for r in (col + 1)..b {
            let v = a[r * b + col].abs();
            if v > best {
                best = v;
                piv = r;
            }
        }
        if best < 1e-15 {
            return None;
        }
        if piv != col {
            for c in 0..b {
                a.swap(col * b + c, piv * b + c);
                inv.swap(col * b + c, piv * b + c);
            }
        }
        let d = a[col * b + col];
        for c in 0..b {
            a[col * b + c] /= d;
            inv[col * b + c] /= d;
        }
        for r in 0..b {
            if r == col {
                continue;
            }
            let f = a[r * b + col];
            if f == 0.0 {
                continue;
            }
            for c in 0..b {
                a[r * b + c] -= f * a[col * b + c];
                inv[r * b + c] -= f * inv[col * b + c];
            }
        }
    }
    Some(inv)
}

fn method_name(m: Method) -> &'static str {
    match m {
        Method::Cem => "cem",
        Method::Ace => "ace",
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    method: Method,
    threshold: Option<f64>,
    mask_output: Option<String>,
}

/// Parses `target_spectrum` from a JSON array of numbers or a comma/space/
/// semicolon-separated string. Returns `None` when the parameter is absent.
fn parse_target_spectrum(args: &ToolArgs) -> Result<Option<Vec<f64>>, ToolError> {
    match args.get("target_spectrum") {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Array(arr)) => {
            if arr.is_empty() {
                return Ok(None);
            }
            let mut out = Vec::with_capacity(arr.len());
            for v in arr {
                let f = v.as_f64().ok_or_else(|| {
                    ToolError::Validation(
                        "'target_spectrum' array must contain only numbers".into(),
                    )
                })?;
                if !f.is_finite() {
                    return Err(ToolError::Validation(
                        "'target_spectrum' values must be finite".into(),
                    ));
                }
                out.push(f);
            }
            Ok(Some(out))
        }
        Some(Value::String(s)) => {
            let s = s.trim();
            if s.is_empty() {
                return Ok(None);
            }
            let mut out = Vec::new();
            for tok in s.split(|c: char| c == ',' || c == ';' || c.is_whitespace()) {
                let tok = tok.trim();
                if tok.is_empty() {
                    continue;
                }
                let f = tok.parse::<f64>().map_err(|_| {
                    ToolError::Validation(format!(
                        "'target_spectrum' contains a non-numeric value '{tok}'"
                    ))
                })?;
                if !f.is_finite() {
                    return Err(ToolError::Validation(
                        "'target_spectrum' values must be finite".into(),
                    ));
                }
                out.push(f);
            }
            if out.is_empty() {
                return Ok(None);
            }
            Ok(Some(out))
        }
        Some(Value::Number(num)) => {
            // A single scalar is not a usable multiband spectrum.
            let f = num.as_f64().unwrap_or(f64::NAN);
            Ok(Some(vec![f]))
        }
        Some(_) => Err(ToolError::Validation(
            "'target_spectrum' must be a JSON array or a comma-separated string of numbers".into(),
        )),
    }
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let method = match args
        .get("method")
        .and_then(Value::as_str)
        .map(|s| s.trim().to_lowercase())
    {
        None => Method::Cem,
        Some(s) if s.is_empty() || s == "cem" => Method::Cem,
        Some(s) if s == "ace" => Method::Ace,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "'method' must be 'cem' or 'ace', got '{other}'"
            )))
        }
    };
    let threshold = match args.get("threshold") {
        None | Some(Value::Null) => None,
        Some(Value::Number(n)) => n.as_f64(),
        Some(Value::String(s)) if s.trim().is_empty() => None,
        Some(Value::String(s)) => Some(
            s.trim()
                .parse::<f64>()
                .map_err(|_| ToolError::Validation("'threshold' must be a number".into()))?,
        ),
        Some(_) => return Err(ToolError::Validation("'threshold' must be a number".into())),
    };
    if let Some(p) = threshold {
        if !(0.0..=100.0).contains(&p) {
            return Err(ToolError::Validation(
                "'threshold' must be a percentile in 0..100".into(),
            ));
        }
    }
    let mask_output = match args.get("mask_output") {
        None | Some(Value::Null) => None,
        Some(Value::String(s)) if s.trim().is_empty() => None,
        Some(Value::String(s)) => Some(s.clone()),
        Some(_) => {
            return Err(ToolError::Validation(
                "'mask_output' must be a string".into(),
            ))
        }
    };
    Ok(Params {
        method,
        threshold,
        mask_output,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{memory_store, CrsInfo, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Build a `bands`-band raster from per-band row-major buffers.
    fn raster_of(rows: usize, cols: usize, bands: &[Vec<f64>]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: bands.len(),
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for (bi, band) in bands.iter().enumerate() {
            for row in 0..rows {
                for col in 0..cols {
                    r.set(
                        bi as isize,
                        row as isize,
                        col as isize,
                        band[row * cols + col],
                    )
                    .unwrap();
                }
            }
        }
        let id = memory_store::put_raster(r);
        memory_store::make_raster_memory_path(&id)
    }

    fn read(path: &str) -> Vec<f64> {
        let r = load_input_raster(path).unwrap();
        let mut v = Vec::new();
        for row in 0..r.rows as isize {
            for col in 0..r.cols as isize {
                v.push(r.get(0, row, col));
            }
        }
        v
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Vec<f64>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = MatchedFilterTargetDetectionTool.run(&args, &ctx()).unwrap();
        (out.clone(), read(out.outputs["output"].as_str().unwrap()))
    }

    /// A 3-band scene with a distinct background and one target pixel. Using
    /// that pixel's exact spectrum as the target, CEM scores it at ~1.0 and far
    /// above the background.
    #[test]
    fn cem_scores_exact_target_at_one() {
        let (rows, cols) = (4, 4);
        let n = rows * cols;
        // Correlated background across 3 bands.
        let mut b1 = vec![0.0; n];
        let mut b2 = vec![0.0; n];
        let mut b3 = vec![0.0; n];
        for i in 0..n {
            let t = 0.2 + i as f64 * 0.03;
            b1[i] = t;
            b2[i] = 0.8 * t;
            b3[i] = 0.5 * t;
        }
        // A spectrally distinct target planted at pixel 9.
        let (t1, t2, t3) = (0.05, 0.9, 0.4);
        b1[9] = t1;
        b2[9] = t2;
        b3[9] = t3;
        let input = raster_of(rows, cols, &[b1, b2, b3]);
        let (_o, scores) = run(json!({
            "input": input,
            "target_spectrum": [t1, t2, t3],
            "method": "cem",
        }));
        let target_score = scores[9];
        assert!(
            (target_score - 1.0).abs() < 1e-6,
            "exact-target CEM score should be 1.0, got {target_score}"
        );
        let bg_max = scores
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != 9)
            .map(|(_, s)| s.abs())
            .fold(0.0, f64::max);
        assert!(
            target_score > bg_max * 2.0,
            "target {target_score} vs background max {bg_max}"
        );
    }

    /// ACE scores lie in [0,1], the exact-target pixel scores ~1.0, and it is
    /// the strongest detection.
    #[test]
    fn ace_bounded_and_peaks_at_target() {
        let (rows, cols) = (5, 5);
        let n = rows * cols;
        let mut b1 = vec![0.0; n];
        let mut b2 = vec![0.0; n];
        let mut b3 = vec![0.0; n];
        for i in 0..n {
            let t = 0.1 + i as f64 * 0.02;
            b1[i] = t;
            b2[i] = 0.7 * t + 0.05;
            b3[i] = 0.3 * t;
        }
        let (t1, t2, t3) = (0.9, 0.1, 0.6);
        b1[12] = t1;
        b2[12] = t2;
        b3[12] = t3;
        let input = raster_of(rows, cols, &[b1, b2, b3]);
        let (_o, scores) = run(json!({
            "input": input,
            "target_spectrum": format!("{t1}, {t2}, {t3}"),
            "method": "ace",
        }));
        for (i, s) in scores.iter().enumerate() {
            assert!(
                (0.0..=1.0).contains(s),
                "ACE score at {i} out of [0,1]: {s}"
            );
        }
        let target_score = scores[12];
        assert!(
            (target_score - 1.0).abs() < 1e-6,
            "exact-target ACE score should be 1.0, got {target_score}"
        );
        let bg_max = scores
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != 12)
            .map(|(_, s)| *s)
            .fold(0.0, f64::max);
        assert!(
            target_score > bg_max,
            "target {target_score} should exceed background max {bg_max}"
        );
    }

    /// A percentile threshold produces a binary mask flagging the top scores,
    /// including the planted target.
    #[test]
    fn threshold_mask_flags_target() {
        let (rows, cols) = (5, 5);
        let n = rows * cols;
        let mut b1 = vec![0.0; n];
        let mut b2 = vec![0.0; n];
        for i in 0..n {
            let t = 1.0 + i as f64 * 0.05;
            b1[i] = t;
            b2[i] = 0.6 * t;
        }
        b1[7] = 0.2;
        b2[7] = 3.0; // distinct target
        let input = raster_of(rows, cols, &[b1, b2]);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input,
            "target_spectrum": [0.2, 3.0],
            "method": "cem",
            "threshold": 90.0,
        }))
        .unwrap();
        let out = MatchedFilterTargetDetectionTool.run(&args, &ctx()).unwrap();
        assert!(out.outputs.contains_key("mask_output"));
        let mask = read(out.outputs["mask_output"].as_str().unwrap());
        assert_eq!(mask[7], 1.0, "the target is flagged");
        let flagged: usize = mask.iter().filter(|v| **v == 1.0).count();
        assert!(
            (1..=4).contains(&flagged),
            "≈top 10% flagged, got {flagged}"
        );
    }

    /// A flat, spectrally-uniform scene with a target that does not match its
    /// direction yields no dominant detection: no pixel scores near 1.
    #[test]
    fn non_matching_scene_has_no_dominant_detection() {
        let (rows, cols) = (4, 4);
        let n = rows * cols;
        // Every pixel is a scaled copy of the same spectrum (a flat scene).
        let mut b1 = vec![0.0; n];
        let mut b2 = vec![0.0; n];
        for i in 0..n {
            let s = 1.0 + i as f64 * 0.1;
            b1[i] = s * 1.0;
            b2[i] = s * 2.0;
        }
        let input = raster_of(rows, cols, &[b1, b2]);
        // Target orthogonal-ish to the scene direction: no pixel should peak.
        let (_o, scores) = run(json!({
            "input": input,
            "target_spectrum": [2.0, 1.0],
            "method": "ace",
        }));
        let max = scores.iter().cloned().fold(0.0_f64, f64::max);
        assert!(max < 0.5, "flat mismatched scene peaked at {max}");
    }

    #[test]
    fn rejects_single_band() {
        let input = raster_of(2, 2, &[vec![1.0, 2.0, 3.0, 4.0]]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "target_spectrum": [1.0, 2.0] }))
                .unwrap();
        assert!(MatchedFilterTargetDetectionTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn rejects_spectrum_length_mismatch() {
        let input = raster_of(2, 2, &[vec![1.0; 4], vec![2.0; 4], vec![3.0; 4]]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "target_spectrum": [1.0, 2.0] }))
                .unwrap();
        assert!(MatchedFilterTargetDetectionTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn rejects_bad_parameters() {
        let input = raster_of(2, 2, &[vec![1.0; 4], vec![2.0; 4]]);
        // Missing target_spectrum.
        let a: ToolArgs = serde_json::from_value(json!({ "input": &input })).unwrap();
        assert!(MatchedFilterTargetDetectionTool.validate(&a).is_err());
        // Bad method.
        let a: ToolArgs = serde_json::from_value(
            json!({ "input": &input, "target_spectrum": [1.0, 2.0], "method": "xyz" }),
        )
        .unwrap();
        assert!(MatchedFilterTargetDetectionTool.validate(&a).is_err());
        // Non-numeric spectrum token.
        let a: ToolArgs =
            serde_json::from_value(json!({ "input": &input, "target_spectrum": "1.0, foo" }))
                .unwrap();
        assert!(MatchedFilterTargetDetectionTool.validate(&a).is_err());
        // Threshold out of range.
        let a: ToolArgs = serde_json::from_value(
            json!({ "input": &input, "target_spectrum": [1.0, 2.0], "threshold": 150.0 }),
        )
        .unwrap();
        assert!(MatchedFilterTargetDetectionTool.validate(&a).is_err());
    }
}
