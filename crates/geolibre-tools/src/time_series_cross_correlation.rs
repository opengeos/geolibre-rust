//! GeoLibre tool: per-cell lagged cross-correlation between two space-time
//! raster stacks.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Time Series Cross Correlation*
//! (Space Time Pattern Mining) — and of the Image Analyst *Multidimensional
//! Raster Correlation* workflow. The space-time cube family already covers a
//! single variable (`generate_trend_raster`, `change_point_detection`,
//! `time_series_forecast`, `local_outlier_analysis`), but nothing relates *two*
//! variables through time. This tool does: for every cell it correlates two
//! aligned time series (e.g. rainfall vs. NDVI) across a range of time lags and
//! reports where they line up best.
//!
//! Each stack is either a single multiband raster (one band per time step) or a
//! comma-separated list of co-registered single-band rasters in time order. The
//! two stacks must share dimensions and length.
//!
//! Per cell, over the two series `A` (`input`) and `B` (`secondary`):
//! * optionally remove a per-series linear trend (`detrend`) and/or an additive
//!   seasonal cycle of period `season_length` (`deseasonalize`);
//! * for each integer lag `L` in `min_lag..=max_lag`, pair `A[t]` with
//!   `B[t+L]` over the valid overlapping window and compute Pearson `r_L`;
//! * pick the lag with the largest `|r_L|` (needing at least `min_valid`
//!   pairs); `output` is that best lag, `corr_output` its (signed) correlation.
//!
//! Optional extra outputs: `corr0_output` (correlation at lag 0) and
//! `pvalue_output`, a pseudo two-sided p-value that corrects for serial
//! autocorrelation via the Pyper & Peterman (1998) effective sample size
//! combined with a Fisher z-transform — so autocorrelated series are not
//! declared significant on the strength of inflated degrees of freedom.
//!
//! Cells whose two series never share `min_valid` valid overlapping
//! observations at any lag become no-data.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::common::{load_input_raster, parse_optional_output, raster_like_with_data};

pub struct TimeSeriesCrossCorrelationTool;

impl Tool for TimeSeriesCrossCorrelationTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "time_series_cross_correlation",
            display_name: "Time Series Cross Correlation",
            summary: "Per-cell lagged cross-correlation between two aligned space-time raster stacks (like ArcGIS's Time Series Cross Correlation / Multidimensional Raster Correlation): for each cell, optionally detrend/deseasonalize, then correlate series A against series B across min_lag..=max_lag and report the best lag and its Pearson correlation, plus optional lag-0 correlation and a Pyper-Peterman effective-sample-size pseudo p-value. The two-variable, lag-aware complement to the single-series generate_trend_raster / change_point_detection cube tools — nothing bundled relates two variables through time.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Stack A: a multiband raster (one band per time step) OR a comma-separated list of co-registered single-band raster paths in time order.",
                    required: true,
                },
                ToolParamSpec {
                    name: "secondary",
                    description: "Stack B: same shape and length as 'input' (multiband raster or comma-separated single-band list).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output best-lag raster (integer lag maximizing |correlation|). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "corr_output",
                    description: "Optional output raster of the (signed) Pearson correlation at the best lag.",
                    required: false,
                },
                ToolParamSpec {
                    name: "corr0_output",
                    description: "Optional output raster of the correlation at lag 0 (contemporaneous).",
                    required: false,
                },
                ToolParamSpec {
                    name: "pvalue_output",
                    description: "Optional output raster of a pseudo two-sided p-value at the best lag (Pyper-Peterman effective sample size + Fisher z).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_lag",
                    description: "Minimum lag in time steps (default -6). A[t] is paired with B[t+lag].",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_lag",
                    description: "Maximum lag in time steps (default 6). Must be >= min_lag.",
                    required: false,
                },
                ToolParamSpec {
                    name: "detrend",
                    description: "If true, subtract a per-series OLS linear trend before correlating (default false).",
                    required: false,
                },
                ToolParamSpec {
                    name: "deseasonalize",
                    description: "If true, subtract a per-series additive seasonal cycle of period 'season_length' before correlating (default false).",
                    required: false,
                },
                ToolParamSpec {
                    name: "season_length",
                    description: "Seasonal period (time steps) used when 'deseasonalize' is true (default 12).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_valid",
                    description: "Minimum valid overlapping pairs required to correlate at a lag (default 5).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read when a stack is given as a comma-separated single-band list (default 1). Ignored for multiband stacks.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "secondary")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let prm = parse_params(args)?;
        let output = parse_optional_output(args, "output")?;
        let corr_out = parse_optional_output(args, "corr_output")?;
        let corr0_out = parse_optional_output(args, "corr0_output")?;
        let pval_out = parse_optional_output(args, "pvalue_output")?;

        let stack_a = Stack::load(require_str(args, "input")?, prm.band)?;
        let stack_b = Stack::load(require_str(args, "secondary")?, prm.band)?;

        if stack_a.n != stack_b.n {
            return Err(ToolError::Validation(format!(
                "stacks differ in length: 'input' has {} time step(s), 'secondary' has {}",
                stack_a.n, stack_b.n
            )));
        }
        if stack_a.rows != stack_b.rows || stack_a.cols != stack_b.cols {
            return Err(ToolError::Validation(format!(
                "stacks differ in size: {}x{} vs {}x{}",
                stack_a.rows, stack_a.cols, stack_b.rows, stack_b.cols
            )));
        }
        let n = stack_a.n;
        if prm.max_lag >= n as i64 || prm.min_lag <= -(n as i64) {
            return Err(ToolError::Validation(format!(
                "lag range [{},{}] exceeds the {n}-step series",
                prm.min_lag, prm.max_lag
            )));
        }

        let (rows, cols) = (stack_a.rows, stack_a.cols);
        let nodata = -9999.0_f64;
        ctx.progress.info(&format!(
            "{n} time step(s), {rows}x{cols}, lags {}..={}{}{}",
            prm.min_lag,
            prm.max_lag,
            if prm.detrend { ", detrend" } else { "" },
            if prm.deseasonalize {
                ", deseasonalize"
            } else {
                ""
            },
        ));

        let mut best_lag = vec![nodata; rows * cols];
        let mut best_corr = vec![nodata; rows * cols];
        let mut corr0 = vec![nodata; rows * cols];
        let mut pval = vec![nodata; rows * cols];

        // Per-cell scratch buffers (valid flag + value).
        let mut a: Vec<Option<f64>> = vec![None; n];
        let mut b: Vec<Option<f64>> = vec![None; n];
        for r in 0..rows {
            for c in 0..cols {
                for t in 0..n {
                    a[t] = stack_a.get(t, r as isize, c as isize);
                    b[t] = stack_b.get(t, r as isize, c as isize);
                }
                if prm.detrend {
                    detrend(&mut a);
                    detrend(&mut b);
                }
                if prm.deseasonalize {
                    deseasonalize(&mut a, prm.season_length);
                    deseasonalize(&mut b, prm.season_length);
                }

                let mut best: Option<(i64, f64, usize)> = None; // (lag, r, npairs)
                let mut r0: Option<f64> = None;
                for lag in prm.min_lag..=prm.max_lag {
                    if let Some((corr, npairs)) = pearson_at_lag(&a, &b, lag, prm.min_valid) {
                        if lag == 0 {
                            r0 = Some(corr);
                        }
                        let better = match best {
                            None => true,
                            Some((_, bc, _)) => corr.abs() > bc.abs(),
                        };
                        if better {
                            best = Some((lag, corr, npairs));
                        }
                    }
                }

                let idx = r * cols + c;
                if let Some((lag, corr, _)) = best {
                    best_lag[idx] = lag as f64;
                    best_corr[idx] = corr;
                    if let Some(v) = r0 {
                        corr0[idx] = v;
                    }
                    if pval_out.is_some() {
                        pval[idx] = pseudo_pvalue(&a, &b, lag, prm.min_valid).unwrap_or(nodata);
                    }
                }
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let template = &stack_a.template;
        let best_lag_r = raster_like_with_data(template, best_lag, nodata, DataType::F32)?;
        let out_path = crate::common::write_or_store_output(best_lag_r, output)?;
        let best_corr_r = raster_like_with_data(template, best_corr, nodata, DataType::F32)?;
        let corr_path = crate::common::write_or_store_output(best_corr_r, corr_out)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        // Best-lag correlation is a first-class result; always emit its handle.
        outputs.insert("corr_output".to_string(), json!(corr_path));
        if let Some(p) = corr0_out {
            let r = raster_like_with_data(template, corr0, nodata, DataType::F32)?;
            outputs.insert(
                "corr0_output".to_string(),
                json!(crate::common::write_or_store_output(r, Some(p))?),
            );
        }
        if let Some(p) = pval_out {
            let r = raster_like_with_data(template, pval, nodata, DataType::F32)?;
            outputs.insert(
                "pvalue_output".to_string(),
                json!(crate::common::write_or_store_output(r, Some(p))?),
            );
        }
        outputs.insert("time_steps".to_string(), json!(n));
        outputs.insert("min_lag".to_string(), json!(prm.min_lag));
        outputs.insert("max_lag".to_string(), json!(prm.max_lag));
        Ok(ToolRunResult { outputs })
    }
}

// ── Stack abstraction ────────────────────────────────────────────────────────

/// A space-time stack: either one multiband raster (band = time) or a list of
/// co-registered single-band rasters (one per time step).
struct Stack {
    rasters: Vec<Raster>,
    multiband: bool,
    band: isize,
    n: usize,
    rows: usize,
    cols: usize,
    template: Raster,
}

impl Stack {
    fn load(spec: &str, band_1based: usize) -> Result<Stack, ToolError> {
        let paths: Vec<&str> = spec
            .split(',')
            .map(str::trim)
            .filter(|p| !p.is_empty())
            .collect();
        if paths.is_empty() {
            return Err(ToolError::Validation("empty stack specification".into()));
        }
        let band = (band_1based.max(1) - 1) as isize;

        if paths.len() == 1 {
            let ras = load_input_raster(paths[0])?;
            if ras.bands >= 2 {
                // Multiband stack: each band is a time step.
                let (rows, cols, n) = (ras.rows, ras.cols, ras.bands);
                let template = ras.clone();
                return Ok(Stack {
                    rasters: vec![ras],
                    multiband: true,
                    band: 0,
                    n,
                    rows,
                    cols,
                    template,
                });
            }
            // Single single-band raster is a degenerate 1-step stack.
            let (rows, cols) = (ras.rows, ras.cols);
            let template = ras.clone();
            return Ok(Stack {
                rasters: vec![ras],
                multiband: false,
                band,
                n: 1,
                rows,
                cols,
                template,
            });
        }

        // Comma-separated list of single-band rasters.
        let rasters: Vec<Raster> = paths
            .iter()
            .map(|p| load_input_raster(p))
            .collect::<Result<_, _>>()?;
        let (rows, cols) = (rasters[0].rows, rasters[0].cols);
        for (i, r) in rasters.iter().enumerate() {
            if r.rows != rows || r.cols != cols {
                return Err(ToolError::Validation(format!(
                    "raster {i} in stack is {}x{}, expected {rows}x{cols}",
                    r.rows, r.cols
                )));
            }
            if band < 0 || band as usize >= r.bands {
                return Err(ToolError::Validation(format!(
                    "band {} out of range for raster {i}",
                    band_1based
                )));
            }
        }
        let template = rasters[0].clone();
        let n = rasters.len();
        Ok(Stack {
            rasters,
            multiband: false,
            band,
            n,
            rows,
            cols,
            template,
        })
    }

    /// Value at time step `t` for cell (row, col), or `None` if no-data / non-finite.
    fn get(&self, t: usize, row: isize, col: isize) -> Option<f64> {
        let (ras, band) = if self.multiband {
            (&self.rasters[0], t as isize)
        } else {
            (&self.rasters[t], self.band)
        };
        let v = ras.get(band, row, col);
        if v != ras.nodata && v.is_finite() {
            Some(v)
        } else {
            None
        }
    }
}

// ── Numerics ─────────────────────────────────────────────────────────────────

/// Pearson correlation of `A[t]` with `B[t+lag]` over valid overlapping `t`.
/// Returns `(r, n_pairs)`, or `None` if fewer than `min_valid` pairs or the
/// correlation is undefined (a constant series).
fn pearson_at_lag(
    a: &[Option<f64>],
    b: &[Option<f64>],
    lag: i64,
    min_valid: usize,
) -> Option<(f64, usize)> {
    let n = a.len() as i64;
    let (mut sx, mut sy, mut sxx, mut syy, mut sxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let mut m = 0usize;
    for t in 0..n {
        let tb = t + lag;
        if tb < 0 || tb >= n {
            continue;
        }
        if let (Some(x), Some(y)) = (a[t as usize], b[tb as usize]) {
            sx += x;
            sy += y;
            sxx += x * x;
            syy += y * y;
            sxy += x * y;
            m += 1;
        }
    }
    if m < min_valid || m < 2 {
        return None;
    }
    let mf = m as f64;
    let cov = sxy - sx * sy / mf;
    let vx = sxx - sx * sx / mf;
    let vy = syy - sy * sy / mf;
    if vx <= 0.0 || vy <= 0.0 {
        return None;
    }
    let r = (cov / (vx * vy).sqrt()).clamp(-1.0, 1.0);
    Some((r, m))
}

/// Pyper & Peterman (1998) effective-sample-size corrected pseudo two-sided
/// p-value for the correlation at `lag`, via a Fisher z-transform. `None` if
/// undefined.
fn pseudo_pvalue(a: &[Option<f64>], b: &[Option<f64>], lag: i64, min_valid: usize) -> Option<f64> {
    // Build the aligned, valid pairs at this lag.
    let n = a.len() as i64;
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for t in 0..n {
        let tb = t + lag;
        if tb < 0 || tb >= n {
            continue;
        }
        if let (Some(x), Some(y)) = (a[t as usize], b[tb as usize]) {
            xs.push(x);
            ys.push(y);
        }
    }
    let m = xs.len();
    if m < min_valid || m < 4 {
        return None;
    }
    let r = pearson(&xs, &ys)?;
    // Effective sample size from serial autocorrelation.
    let jmax = (m / 5).max(1);
    let mut s = 0.0;
    for j in 1..=jmax {
        if j >= m {
            break;
        }
        let rxx = autocorr(&xs, j);
        let ryy = autocorr(&ys, j);
        s += rxx * ryy;
    }
    let mf = m as f64;
    let inv_neff = 1.0 / mf + (2.0 / mf) * s;
    let n_eff = if inv_neff > 0.0 {
        (1.0 / inv_neff).clamp(2.0, mf)
    } else {
        mf
    };
    if n_eff <= 3.0 || r.abs() >= 1.0 {
        return Some(if r.abs() >= 1.0 { 0.0 } else { 1.0 });
    }
    let z = 0.5 * ((1.0 + r) / (1.0 - r)).ln() * (n_eff - 3.0).sqrt();
    let p = 2.0 * (1.0 - normal_cdf(z.abs()));
    Some(p.clamp(0.0, 1.0))
}

fn pearson(x: &[f64], y: &[f64]) -> Option<f64> {
    let n = x.len() as f64;
    if n < 2.0 {
        return None;
    }
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let (mut cov, mut vx, mut vy) = (0.0, 0.0, 0.0);
    for i in 0..x.len() {
        let dx = x[i] - mx;
        let dy = y[i] - my;
        cov += dx * dy;
        vx += dx * dx;
        vy += dy * dy;
    }
    if vx <= 0.0 || vy <= 0.0 {
        return None;
    }
    Some((cov / (vx * vy).sqrt()).clamp(-1.0, 1.0))
}

/// Modified (Pyper-Peterman) lag-`j` autocorrelation: the lag-j cross-product is
/// averaged over its `m-j` terms, normalized by the full-series variance.
fn autocorr(x: &[f64], j: usize) -> f64 {
    let m = x.len();
    if j >= m {
        return 0.0;
    }
    let mean = x.iter().sum::<f64>() / m as f64;
    let mut num = 0.0;
    let mut den = 0.0;
    for &v in x.iter() {
        den += (v - mean).powi(2);
    }
    for i in 0..(m - j) {
        num += (x[i] - mean) * (x[i + j] - mean);
    }
    if den <= 0.0 {
        return 0.0;
    }
    // num averaged over (m-j), den over m -> the modified estimator.
    (num / (m - j) as f64) / (den / m as f64)
}

/// Subtract a per-series OLS linear trend in place, over valid observations.
fn detrend(series: &mut [Option<f64>]) {
    let n = series.len();
    let (mut st, mut sv, mut n_valid) = (0.0, 0.0, 0.0);
    for (t, s) in series.iter().enumerate() {
        if let Some(v) = s {
            st += t as f64;
            sv += *v;
            n_valid += 1.0;
        }
    }
    if n_valid < 2.0 {
        return;
    }
    let mt = st / n_valid;
    let mv = sv / n_valid;
    let (mut sxx, mut sxy) = (0.0, 0.0);
    for (t, s) in series.iter().enumerate() {
        if let Some(v) = s {
            let dt = t as f64 - mt;
            sxx += dt * dt;
            sxy += dt * (*v - mv);
        }
    }
    if sxx <= 0.0 {
        return;
    }
    let slope = sxy / sxx;
    let intercept = mv - slope * mt;
    for (t, s) in series.iter_mut().enumerate().take(n) {
        if let Some(v) = s {
            *s = Some(*v - (intercept + slope * t as f64));
        }
    }
}

/// Subtract an additive seasonal cycle of `period` in place: for each phase
/// `t % period`, subtract the mean of that phase's valid observations.
fn deseasonalize(series: &mut [Option<f64>], period: usize) {
    if period < 2 {
        return;
    }
    let mut sums = vec![0.0; period];
    let mut counts = vec![0.0; period];
    for (t, s) in series.iter().enumerate() {
        if let Some(v) = s {
            let ph = t % period;
            sums[ph] += *v;
            counts[ph] += 1.0;
        }
    }
    for (t, s) in series.iter_mut().enumerate() {
        if let Some(v) = s {
            let ph = t % period;
            if counts[ph] > 0.0 {
                *s = Some(*v - sums[ph] / counts[ph]);
            }
        }
    }
}

fn normal_cdf(x: f64) -> f64 {
    0.5 * (1.0 + erf(x / std::f64::consts::SQRT_2))
}

fn erf(x: f64) -> f64 {
    // Abramowitz & Stegun 7.1.26.
    let t = 1.0 / (1.0 + 0.3275911 * x.abs());
    let y = 1.0
        - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t
            + 0.254829592)
            * t
            * (-x * x).exp();
    if x >= 0.0 {
        y
    } else {
        -y
    }
}

// ── Parameters ───────────────────────────────────────────────────────────────

struct Params {
    min_lag: i64,
    max_lag: i64,
    detrend: bool,
    deseasonalize: bool,
    season_length: usize,
    min_valid: usize,
    band: usize,
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    match args.get(key).and_then(Value::as_str) {
        Some(s) if !s.trim().is_empty() => Ok(s),
        _ => Err(ToolError::Validation(format!(
            "missing required parameter '{key}'"
        ))),
    }
}

/// Parses an optional integer that may arrive as a JSON number or a string.
fn opt_i64(args: &ToolArgs, key: &str, default: i64) -> Result<i64, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(default),
        Some(Value::Number(n)) => n
            .as_i64()
            .or_else(|| n.as_f64().map(|f| f.round() as i64))
            .ok_or_else(|| ToolError::Validation(format!("'{key}' must be an integer"))),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(default),
        Some(Value::String(s)) => s
            .trim()
            .parse::<i64>()
            .map_err(|_| ToolError::Validation(format!("'{key}' must be an integer"))),
        _ => Err(ToolError::Validation(format!("'{key}' must be an integer"))),
    }
}

fn opt_usize(args: &ToolArgs, key: &str, default: usize) -> Result<usize, ToolError> {
    let v = opt_i64(args, key, default as i64)?;
    if v < 0 {
        return Err(ToolError::Validation(format!("'{key}' must be >= 0")));
    }
    Ok(v as usize)
}

fn opt_bool(args: &ToolArgs, key: &str, default: bool) -> Result<bool, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(default),
        Some(Value::Bool(b)) => Ok(*b),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(default),
            "true" | "1" | "yes" | "t" => Ok(true),
            "false" | "0" | "no" | "f" => Ok(false),
            other => Err(ToolError::Validation(format!(
                "'{key}' must be a boolean, got '{other}'"
            ))),
        },
        Some(Value::Number(n)) => Ok(n.as_f64().map(|f| f != 0.0).unwrap_or(default)),
        _ => Err(ToolError::Validation(format!("'{key}' must be a boolean"))),
    }
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let min_lag = opt_i64(args, "min_lag", -6)?;
    let max_lag = opt_i64(args, "max_lag", 6)?;
    if max_lag < min_lag {
        return Err(ToolError::Validation(format!(
            "'max_lag' ({max_lag}) must be >= 'min_lag' ({min_lag})"
        )));
    }
    let detrend = opt_bool(args, "detrend", false)?;
    let deseasonalize = opt_bool(args, "deseasonalize", false)?;
    let season_length = opt_usize(args, "season_length", 12)?;
    if deseasonalize && season_length < 2 {
        return Err(ToolError::Validation(
            "'season_length' must be >= 2 when 'deseasonalize' is true".into(),
        ));
    }
    let min_valid = opt_usize(args, "min_valid", 5)?.max(2);
    let band = opt_usize(args, "band", 1)?.max(1);
    Ok(Params {
        min_lag,
        max_lag,
        detrend,
        deseasonalize,
        season_length,
        min_valid,
        band,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a multiband raster of `cols*rows` cells and `bands` time steps
    /// from a per-cell series generator, returning a memory path.
    fn stack_from<F>(cols: usize, rows: usize, bands: usize, mut f: F) -> String
    where
        F: FnMut(usize, usize, usize) -> f64,
    {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands,
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
        for band in 0..bands {
            for row in 0..rows {
                for col in 0..cols {
                    r.set(band as isize, row as isize, col as isize, f(row, col, band))
                        .unwrap();
                }
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (Raster, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = TimeSeriesCrossCorrelationTool.run(&args, &ctx()).unwrap();
        let lag = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        let corr = load_input_raster(out.outputs["corr_output"].as_str().unwrap()).unwrap();
        (lag, corr)
    }

    /// A single cell where B is A delayed by a known lag recovers that lag with
    /// a near-perfect correlation.
    #[test]
    fn recovers_planted_lag() {
        let n = 24;
        let planted = 3i64;
        // A[t] = sin(t) + t*0.1 ; B[t] = A[t - planted]. B before the shift is 0.
        let a = |t: usize| (t as f64 * 0.7).sin() + 0.1 * t as f64;
        let stack_a = stack_from(1, 1, n, |_, _, t| a(t));
        let stack_b = stack_from(1, 1, n, |_, _, t| {
            if (t as i64) - planted >= 0 {
                a((t as i64 - planted) as usize)
            } else {
                0.0
            }
        });
        let (lag, corr) = run(json!({
            "input": stack_a,
            "secondary": stack_b,
            "min_lag": -6,
            "max_lag": 6,
        }));
        assert_eq!(
            lag.get(0, 0, 0),
            planted as f64,
            "should recover planted lag"
        );
        assert!(
            corr.get(0, 0, 0) > 0.95,
            "best correlation should be high, got {}",
            corr.get(0, 0, 0)
        );
    }

    /// Lag 0: perfectly correlated identical series -> best lag 0, r = 1.
    #[test]
    fn identical_series_lag_zero() {
        let n = 16;
        let a = |t: usize| (t as f64).cos() + 0.3 * t as f64;
        let sa = stack_from(1, 1, n, |_, _, t| a(t));
        let sb = stack_from(1, 1, n, |_, _, t| a(t));
        let (lag, corr) = run(json!({ "input": sa, "secondary": sb, "min_lag": -4, "max_lag": 4 }));
        assert_eq!(lag.get(0, 0, 0), 0.0);
        assert!((corr.get(0, 0, 0) - 1.0).abs() < 1e-6);
    }

    /// Cells that never reach `min_valid` overlapping pairs become no-data.
    #[test]
    fn insufficient_data_is_nodata() {
        let n = 6;
        // Mostly nodata: only one valid step in each stack -> never enough pairs.
        let sa = stack_from(1, 1, n, |_, _, t| if t == 0 { 1.0 } else { -9999.0 });
        let sb = stack_from(1, 1, n, |_, _, t| if t == 0 { 2.0 } else { -9999.0 });
        let (lag, _corr) = run(
            json!({ "input": sa, "secondary": sb, "min_lag": -2, "max_lag": 2, "min_valid": 5 }),
        );
        assert_eq!(lag.get(0, 0, 0), -9999.0);
    }

    /// Detrending removes a shared linear ramp so the residual (anti-correlated)
    /// signal drives the best correlation negative.
    #[test]
    fn detrend_exposes_negative_residual() {
        let n = 24;
        // Both share a strong upward ramp; residuals are opposite-sign sinusoids.
        let a = |t: usize| 5.0 * t as f64 + (t as f64 * 0.9).sin();
        let b = |t: usize| 5.0 * t as f64 - (t as f64 * 0.9).sin();
        let sa = stack_from(1, 1, n, |_, _, t| a(t));
        let sb = stack_from(1, 1, n, |_, _, t| b(t));
        // Without detrend the shared ramp dominates -> strong positive at lag 0.
        let (_lag0, corr0) = run(
            json!({ "input": sa.clone(), "secondary": sb.clone(), "min_lag": 0, "max_lag": 0 }),
        );
        assert!(corr0.get(0, 0, 0) > 0.9, "ramp dominates without detrend");
        // With detrend the anti-correlated residual shows through at lag 0.
        let (_lagd, corrd) = run(json!({
            "input": sa,
            "secondary": sb,
            "min_lag": 0,
            "max_lag": 0,
            "detrend": true,
        }));
        assert!(
            corrd.get(0, 0, 0) < -0.5,
            "detrended residuals are anti-correlated, got {}",
            corrd.get(0, 0, 0)
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            TimeSeriesCrossCorrelationTool.validate(&args)
        };
        assert!(bad(json!({})).is_err()); // missing input/secondary
        assert!(bad(json!({ "input": "a.tif" })).is_err()); // missing secondary
        assert!(
            bad(json!({ "input": "a.tif", "secondary": "b.tif", "min_lag": 5, "max_lag": 2 }))
                .is_err()
        ); // max < min
        assert!(bad(
            json!({ "input": "a.tif", "secondary": "b.tif", "deseasonalize": true, "season_length": 1 })
        )
        .is_err()); // season_length too small
        assert!(bad(json!({ "input": "a.tif", "secondary": "b.tif" })).is_ok());
    }
}
