//! GeoLibre tool: Reed–Xiaoli (RX) anomaly detector for multiband imagery.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Detect Image Anomalies* (Image
//! Analyst). The bundled remote-sensing suite *transforms* imagery
//! (`principal_component_analysis`, `minimum_noise_fraction`,
//! `linear_spectral_unmixing`) but never *scores anomalies*. The RX detector is
//! the standard unsupervised anomaly score: each pixel's squared Mahalanobis
//! distance to the scene's band statistics, `(x − μ)ᵀ Σ⁻¹ (x − μ)`, needs no
//! training data and complements `spectral_index`.
//!
//! * `global` (default) — μ and Σ over all valid pixels of the whole scene.
//! * `local` — μ and Σ recomputed in a moving `window` around each pixel
//!   (RX-LRX), catching anomalies against their local background.
//!
//! The covariance is inverted with a hand-rolled Gauss–Jordan solve (plus a tiny
//! ridge for numerical stability). Output is a single-band F32 anomaly-score
//! raster; with `threshold` (a percentile in 0..100) a second binary mask marks
//! the top scores.

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
/// Ridge added to the covariance diagonal (× its mean) for invertibility.
const RIDGE: f64 = 1e-6;

#[derive(Clone, Copy, PartialEq)]
enum Mode {
    Global,
    Local,
}

pub struct DetectImageAnomaliesTool;

impl Tool for DetectImageAnomaliesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "detect_image_anomalies",
            display_name: "Detect Image Anomalies",
            summary: "Score every pixel of a multiband image by its Mahalanobis distance to the scene (or a moving-window) band statistics — the Reed–Xiaoli (RX) anomaly detector — with an optional percentile threshold mask, like ArcGIS Detect Image Anomalies.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input multiband raster (2+ bands).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output anomaly-score raster (F32). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "mode",
                    description: "'global' (whole-scene statistics, default) or 'local' (moving-window statistics).",
                    required: false,
                },
                ToolParamSpec {
                    name: "window",
                    description: "Odd window size (pixels) for local mode (default 15).",
                    required: false,
                },
                ToolParamSpec {
                    name: "threshold",
                    description: "Optional percentile (0-100) of the score distribution above which a binary anomaly mask is written.",
                    required: false,
                },
                ToolParamSpec {
                    name: "mask_output",
                    description: "Output path for the binary anomaly mask (requires 'threshold'). If omitted with a threshold, stored in memory.",
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
        let input = args.get("input").and_then(Value::as_str).unwrap();
        let output = parse_optional_output(args, "output")?;
        let prm = parse_params(args)?;

        let raster = load_input_raster(input)?;
        let b = raster.bands;
        if b < 2 {
            return Err(ToolError::Validation(format!(
                "RX needs a multiband image; input has {b} band(s)"
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
            "RX anomaly detection ({b}-band, {} mode)",
            mode_name(prm.mode)
        ));

        let scores = match prm.mode {
            Mode::Global => rx_global(&stack, &valid, b, n)?,
            Mode::Local => rx_local(&stack, &valid, b, rows, cols, prm.window),
        };

        let data: Vec<f64> = (0..n)
            .map(|i| if valid[i] { scores[i] } else { OUT_NODATA })
            .collect();
        let out = raster_like_with_data(&raster, data.clone(), OUT_NODATA, DataType::F32)?;
        let out_path = write_or_store_output(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("bands".to_string(), json!(b));
        outputs.insert("mode".to_string(), json!(mode_name(prm.mode)));
        let valid_count = valid.iter().filter(|v| **v).count();
        outputs.insert("valid_pixels".to_string(), json!(valid_count));

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
            let anomalies = (0..n).filter(|&i| valid[i] && scores[i] >= cut).count();
            outputs.insert("mask_output".to_string(), json!(mask_path));
            outputs.insert("threshold_score".to_string(), json!(cut));
            outputs.insert("anomaly_pixels".to_string(), json!(anomalies));
        }

        Ok(ToolRunResult { outputs })
    }
}

// ── RX detector ────────────────────────────────────────────────────────────────

/// Global RX: one mean vector and covariance over all valid pixels.
fn rx_global(
    stack: &[Vec<f64>],
    valid: &[bool],
    b: usize,
    n: usize,
) -> Result<Vec<f64>, ToolError> {
    let (mean, cov) = mean_cov(stack, valid, b, (0..n).filter(|&i| valid[i]));
    let inv = invert(&cov, b)
        .ok_or_else(|| ToolError::Execution("covariance matrix is singular".to_string()))?;
    let mut scores = vec![0.0; n];
    let mut x = vec![0.0; b];
    for (i, score) in scores.iter_mut().enumerate() {
        if !valid[i] {
            continue;
        }
        for (k, xk) in x.iter_mut().enumerate() {
            *xk = stack[k][i] - mean[k];
        }
        *score = mahalanobis(&x, &inv, b);
    }
    Ok(scores)
}

/// Local RX: mean/covariance recomputed in a moving window around each pixel.
fn rx_local(
    stack: &[Vec<f64>],
    valid: &[bool],
    b: usize,
    rows: usize,
    cols: usize,
    window: usize,
) -> Vec<f64> {
    let half = (window / 2) as isize;
    let mut scores = vec![0.0; rows * cols];
    let mut x = vec![0.0; b];
    for row in 0..rows as isize {
        for col in 0..cols as isize {
            let i = row as usize * cols + col as usize;
            if !valid[i] {
                continue;
            }
            let idxs = (row - half..=row + half).flat_map(|r| {
                (col - half..=col + half).filter_map(move |c| {
                    if r >= 0 && c >= 0 && r < rows as isize && c < cols as isize {
                        Some(r as usize * cols + c as usize)
                    } else {
                        None
                    }
                })
            });
            let win: Vec<usize> = idxs.filter(|&j| valid[j]).collect();
            if win.len() <= b {
                continue; // not enough samples for a covariance
            }
            let (mean, cov) = mean_cov(stack, valid, b, win.iter().copied());
            let Some(inv) = invert(&cov, b) else {
                continue;
            };
            for (k, xk) in x.iter_mut().enumerate() {
                *xk = stack[k][i] - mean[k];
            }
            scores[i] = mahalanobis(&x, &inv, b);
        }
    }
    scores
}

/// Mean vector and (population) covariance over the given pixel indices.
fn mean_cov(
    stack: &[Vec<f64>],
    _valid: &[bool],
    b: usize,
    idxs: impl Iterator<Item = usize> + Clone,
) -> (Vec<f64>, Vec<f64>) {
    let mut mean = vec![0.0; b];
    let mut count = 0.0;
    for i in idxs.clone() {
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
    for i in idxs {
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
    // Ridge for stability.
    let diag_mean: f64 = (0..b).map(|k| cov[k * b + k]).sum::<f64>() / b as f64;
    let ridge = RIDGE * diag_mean.max(1e-12);
    for k in 0..b {
        cov[k * b + k] += ridge;
    }
    (mean, cov)
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
        // Partial pivot.
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

fn mode_name(m: Mode) -> &'static str {
    match m {
        Mode::Global => "global",
        Mode::Local => "local",
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    mode: Mode,
    window: usize,
    threshold: Option<f64>,
    mask_output: Option<String>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let mode = match args
        .get("mode")
        .and_then(Value::as_str)
        .map(|s| s.trim().to_lowercase())
    {
        None => Mode::Global,
        Some(s) if s.is_empty() || s == "global" => Mode::Global,
        Some(s) if s == "local" => Mode::Local,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "'mode' must be 'global' or 'local', got '{other}'"
            )))
        }
    };
    let window = match args.get("window") {
        None | Some(Value::Null) => 15,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(15).max(3) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 15,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'window' must be an integer".into()))?
            .max(3),
        Some(_) => return Err(ToolError::Validation("'window' must be a number".into())),
    };
    let window = if window % 2 == 0 { window + 1 } else { window };
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
        mode,
        window,
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
        let out = DetectImageAnomaliesTool.run(&args, &ctx()).unwrap();
        (out.clone(), read(out.outputs["output"].as_str().unwrap()))
    }

    /// The single off-trend pixel scores far higher than the correlated background.
    #[test]
    fn flags_the_spectral_outlier() {
        // 4x4, 2 bands. Background: band2 = band1 (correlated). One pixel breaks it.
        let n = 16;
        let mut b1 = vec![0.0; n];
        let mut b2 = vec![0.0; n];
        for i in 0..n {
            b1[i] = (i as f64) * 0.1;
            b2[i] = b1[i]; // perfectly correlated background
        }
        // pixel 5 is an anomaly: high band1, low band2.
        b1[5] = 5.0;
        b2[5] = 0.0;
        let input = raster_of(4, 4, &[b1, b2]);
        let (_o, scores) = run(json!({ "input": input, "mode": "global" }));
        let anom = scores[5];
        let others_max = scores
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != 5)
            .map(|(_, s)| *s)
            .fold(0.0, f64::max);
        assert!(
            anom > others_max * 5.0,
            "anomaly {anom} vs others max {others_max}"
        );
    }

    /// A percentile threshold produces a binary mask flagging the top scores.
    #[test]
    fn threshold_mask_flags_top_scores() {
        let n = 25;
        let mut b1 = vec![0.0; n];
        let mut b2 = vec![0.0; n];
        // Distinct correlated background (a gradient, so scores are distinct).
        for i in 0..n {
            b1[i] = 1.0 + i as f64 * 0.05;
            b2[i] = b1[i];
        }
        b1[12] = 9.0;
        b2[12] = 0.0; // strong anomaly
        let input = raster_of(5, 5, &[b1, b2]);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input, "threshold": 90.0
        }))
        .unwrap();
        let out = DetectImageAnomaliesTool.run(&args, &ctx()).unwrap();
        assert!(out.outputs.contains_key("mask_output"));
        let mask = read(out.outputs["mask_output"].as_str().unwrap());
        assert_eq!(mask[12], 1.0, "the anomaly is flagged");
        let flagged: usize = mask.iter().filter(|v| **v == 1.0).count();
        assert!(
            (1..=4).contains(&flagged),
            "≈top 10% flagged, got {flagged}"
        );
    }

    /// Local mode runs and scores the local outlier.
    #[test]
    fn local_mode_runs() {
        let n = 49;
        let mut b1 = vec![2.0; n];
        let mut b2 = vec![2.0; n];
        b1[24] = 8.0;
        b2[24] = 1.0;
        let input = raster_of(7, 7, &[b1, b2]);
        let (_o, scores) = run(json!({ "input": input, "mode": "local", "window": 5 }));
        let anom = scores[24];
        assert!(anom > 0.0, "local anomaly scored");
    }

    #[test]
    fn rejects_single_band() {
        let input = raster_of(2, 2, &[vec![1.0, 2.0, 3.0, 4.0]]);
        let args: ToolArgs = serde_json::from_value(json!({ "input": input })).unwrap();
        assert!(DetectImageAnomaliesTool.run(&args, &ctx()).is_err());
    }
}
