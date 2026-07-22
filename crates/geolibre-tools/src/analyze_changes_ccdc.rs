//! GeoLibre tool: CCDC continuous change detection on a dense image time series.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Analyze Changes Using CCDC* (Image
//! Analyst) — Continuous Change Detection and Classification (Zhu & Woodcock
//! 2014). CCDC is the fourth member of the GeoLibre change-detection family
//! (`generate_trend_raster`, `change_point_detection`, `landtrendr`) and the
//! only one that models intra-annual **seasonality**: LandTrendr segments a
//! yearly series and cannot handle sub-annual observations, while the bundled
//! `change_vector_analysis`/`pca_based_change_detection` are two-date methods.
//!
//! Each pixel's ordered observations are fit with a harmonic (seasonal)
//! regression — a constant, a linear trend, and `harmonic_order` sine/cosine
//! pairs at the annual frequency and its harmonics — by ordinary least squares
//! (closed-form normal equations solved with Gaussian elimination; no
//! linear-algebra crate). The tool then slides through time flagging
//! observations whose residual exceeds `change_threshold × RMSE`; after
//! `min_consecutive` flagged observations it declares a **break** at the first
//! flagged date and refits a fresh segment from there.
//!
//! Outputs: a per-pixel **break-count** raster (primary), plus optional rasters
//! for the **first break date**, and the initial-segment model's **RMSE**,
//! **slope** (trend) and first-harmonic **amplitude**. Pixels with fewer than
//! `min_observations` valid observations become no-data.

use std::collections::BTreeMap;
use std::f64::consts::TAU;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::common::{band_to_vec, load_input_raster, parse_optional_output, raster_like_with_data};

pub struct AnalyzeChangesCcdcTool;

impl Tool for AnalyzeChangesCcdcTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "analyze_changes_ccdc",
            display_name: "Analyze Changes Using CCDC",
            summary: "CCDC continuous change detection on a dense image time series (like ArcGIS Analyze Changes Using CCDC): per-pixel harmonic (seasonal) OLS regression, then a sliding residual test (k×RMSE, min consecutive anomalies) that declares breaks and refits segments. Emits break-count, first-break-date, and initial-segment RMSE/slope/amplitude rasters — the seasonal, multi-break change history the yearly landtrendr and two-date change_vector_analysis can't produce.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "A single multiband raster (each band is a time slice) OR a comma-separated list of co-registered single-band raster paths, in time order.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output break-count raster (breaks detected per pixel). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "dates",
                    description: "Comma-separated numeric dates (decimal years or day numbers) matching the slices. Default 0,1,2,...",
                    required: false,
                },
                ToolParamSpec {
                    name: "period",
                    description: "Length of one seasonal cycle in the units of 'dates' (default 1.0 for decimal years; use 365 for day numbers).",
                    required: false,
                },
                ToolParamSpec {
                    name: "harmonic_order",
                    description: "Number of sine/cosine harmonic pairs at the annual frequency (default 1 = annual only; 2 adds the semi-annual term).",
                    required: false,
                },
                ToolParamSpec {
                    name: "change_threshold",
                    description: "Anomaly threshold as a multiple k of the segment RMSE (default 3.0). A residual > k×RMSE flags an observation.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_consecutive",
                    description: "Number of consecutive flagged observations required to declare a break (default 3).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_observations",
                    description: "Minimum valid observations used to initialize a segment's harmonic model (default 12). Pixels with fewer valid observations become no-data.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "For a comma-separated list input: 1-based band to read from each raster (default 1). Ignored for a single multiband input.",
                    required: false,
                },
                ToolParamSpec {
                    name: "break_date_output",
                    description: "Optional raster of the date of the first detected break per pixel (no-data where no break).",
                    required: false,
                },
                ToolParamSpec {
                    name: "rmse_output",
                    description: "Optional raster of the initial-segment harmonic model RMSE.",
                    required: false,
                },
                ToolParamSpec {
                    name: "slope_output",
                    description: "Optional raster of the initial-segment linear trend (slope) coefficient.",
                    required: false,
                },
                ToolParamSpec {
                    name: "amplitude_output",
                    description: "Optional raster of the initial-segment first-harmonic amplitude sqrt(a²+b²).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        parse_inputs(args)?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let paths = parse_inputs(args)?;
        let prm = parse_params(args)?;
        let output = parse_optional_output(args, "output")?;
        let break_date_out = parse_optional_output(args, "break_date_output")?;
        let rmse_out = parse_optional_output(args, "rmse_output")?;
        let slope_out = parse_optional_output(args, "slope_output")?;
        let amp_out = parse_optional_output(args, "amplitude_output")?;

        // Load the stack into per-slice row-major buffers sharing one geometry.
        let (template, slices, nodata) = load_stack(&paths, prm.band)?;
        let n_slices = slices.len();

        let dates = match &prm.dates {
            Some(d) if d.len() != n_slices => {
                return Err(ToolError::Validation(format!(
                    "{} dates for {} slices",
                    d.len(),
                    n_slices
                )))
            }
            Some(d) => d.clone(),
            None => (0..n_slices).map(|i| i as f64).collect(),
        };

        let (rows, cols) = (template.rows, template.cols);
        let p = 2 + 2 * prm.harmonic_order;
        let min_obs = prm.min_observations.max(p + 1);

        ctx.progress.info(&format!(
            "CCDC over {n_slices} slice(s), {rows}x{cols}, harmonic order {}, k={}, min_consec={}, min_obs={min_obs}",
            prm.harmonic_order, prm.change_threshold, prm.min_consecutive
        ));

        let mut count_r = vec![nodata; rows * cols];
        let mut date_r = vec![nodata; rows * cols];
        let mut rmse_r = vec![nodata; rows * cols];
        let mut slope_r = vec![nodata; rows * cols];
        let mut amp_r = vec![nodata; rows * cols];

        let mut ts: Vec<f64> = Vec::with_capacity(n_slices);
        let mut ys: Vec<f64> = Vec::with_capacity(n_slices);
        for r in 0..rows {
            for c in 0..cols {
                let idx = r * cols + c;
                ts.clear();
                ys.clear();
                for (k, band) in slices.iter().enumerate() {
                    let v = band[idx];
                    if v != nodata && v.is_finite() {
                        ts.push(dates[k]);
                        ys.push(v);
                    }
                }
                if ys.len() < min_obs {
                    continue; // insufficient data -> no-data
                }
                let seg = ccdc_pixel(&ts, &ys, &prm, min_obs);
                count_r[idx] = seg.breaks.len() as f64;
                if let Some(d) = seg.breaks.first() {
                    date_r[idx] = *d;
                }
                if let Some(m) = seg.initial {
                    rmse_r[idx] = m.rmse;
                    slope_r[idx] = m.slope;
                    amp_r[idx] = m.amplitude;
                }
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let count = raster_like_with_data(&template, count_r, nodata, DataType::F32)?;
        let out_path = crate::common::write_or_store_output(count, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));

        for (buf, path, key) in [
            (date_r, break_date_out, "break_date_output"),
            (rmse_r, rmse_out, "rmse_output"),
            (slope_r, slope_out, "slope_output"),
            (amp_r, amp_out, "amplitude_output"),
        ] {
            if let Some(p) = path {
                let r = raster_like_with_data(&template, buf, nodata, DataType::F32)?;
                outputs.insert(
                    key.to_string(),
                    json!(crate::common::write_or_store_output(r, Some(p))?),
                );
            }
        }
        Ok(ToolRunResult { outputs })
    }
}

/// Per-pixel CCDC result: the detected break dates (in input date units) and the
/// initial segment's fitted model summary.
struct PixelResult {
    breaks: Vec<f64>,
    initial: Option<ModelSummary>,
}

struct ModelSummary {
    rmse: f64,
    slope: f64,
    amplitude: f64,
}

/// Absolute floor on a segment's RMSE, so a near-perfect fit does not turn
/// floating-point round-off into spurious anomalies. Real data always sits well
/// above this.
const RMSE_FLOOR: f64 = 1e-6;

/// Runs the CCDC sliding-window detection for one pixel's valid observations
/// `(t, y)` (already time-ordered). A segment's harmonic model is initialized on
/// its first `min_obs` observations and then **continuously refit** as each
/// subsequent clean observation is accepted (so RMSE tracks the true noise).
/// Observations whose residual exceeds `change_threshold × RMSE` are flagged;
/// `min_consecutive` consecutive flags declare a break at the first flagged
/// date, after which a fresh segment is fit from the break onward.
fn ccdc_pixel(t: &[f64], y: &[f64], prm: &Params, min_obs: usize) -> PixelResult {
    let n = t.len();
    let order = prm.harmonic_order;
    let period = prm.period;
    let mut breaks: Vec<f64> = Vec::new();
    let mut initial: Option<ModelSummary> = None;
    let mut seg_start = 0usize;

    while n - seg_start >= min_obs {
        let t0 = t[seg_start];
        // Training set for the current segment, grown as clean observations are
        // accepted; seeded with the first `min_obs` of the segment.
        let mut train: Vec<usize> = (seg_start..seg_start + min_obs).collect();
        let mut fit = match fit_indices(t, y, &train, t0, period, order) {
            Some(f) => f,
            None => break, // singular design; stop
        };
        if initial.is_none() {
            let amplitude = (fit.beta[2] * fit.beta[2] + fit.beta[3] * fit.beta[3]).sqrt();
            initial = Some(ModelSummary {
                rmse: fit.rmse,
                slope: fit.beta[1],
                amplitude,
            });
        }

        let mut consec = 0usize;
        let mut first_anom: Option<usize> = None;
        let mut break_idx: Option<usize> = None;
        for i in (seg_start + min_obs)..n {
            let thresh = prm.change_threshold * fit.rmse.max(RMSE_FLOOR);
            let pred = predict(t[i], t0, period, order, &fit.beta);
            if (y[i] - pred).abs() > thresh {
                if consec == 0 {
                    first_anom = Some(i);
                }
                consec += 1;
                if consec >= prm.min_consecutive {
                    break_idx = first_anom;
                    break;
                }
            } else {
                consec = 0;
                first_anom = None;
                // Accept this clean observation and refit the segment model.
                train.push(i);
                if let Some(f) = fit_indices(t, y, &train, t0, period, order) {
                    fit = f;
                }
            }
        }
        match break_idx {
            Some(bi) => {
                breaks.push(t[bi]);
                seg_start = bi; // start a fresh segment at the break
            }
            None => break, // no further breaks
        }
    }

    PixelResult { breaks, initial }
}

/// Fits the harmonic model on the observations at the given `indices`.
fn fit_indices(
    t: &[f64],
    y: &[f64],
    indices: &[usize],
    t0: f64,
    period: f64,
    order: usize,
) -> Option<HarmonicFit> {
    let tt: Vec<f64> = indices.iter().map(|&i| t[i]).collect();
    let yy: Vec<f64> = indices.iter().map(|&i| y[i]).collect();
    fit_harmonic(&tt, &yy, t0, period, order)
}

/// A fitted harmonic model: coefficient vector and RMSE.
struct HarmonicFit {
    beta: Vec<f64>,
    rmse: f64,
}

/// Design-matrix row for date `t`: [1, (t - t0), sin(w·t), cos(w·t), ...] with
/// one sine/cosine pair per harmonic `k = 1..=order` at angular frequency
/// `2π·k / period`. The trend term is centered at the segment start `t0` for
/// numerical conditioning; harmonics use absolute `t` so seasonality is
/// phase-anchored to the calendar.
fn design(t: f64, t0: f64, period: f64, order: usize) -> Vec<f64> {
    let mut row = Vec::with_capacity(2 + 2 * order);
    row.push(1.0);
    row.push(t - t0);
    for k in 1..=order {
        let w = TAU * k as f64 * t / period;
        row.push(w.sin());
        row.push(w.cos());
    }
    row
}

/// Predicts the modeled value at date `t` from coefficient vector `beta`.
fn predict(t: f64, t0: f64, period: f64, order: usize, beta: &[f64]) -> f64 {
    design(t, t0, period, order)
        .iter()
        .zip(beta)
        .map(|(x, b)| x * b)
        .sum()
}

/// Ordinary least squares fit of the harmonic basis via the normal equations
/// `(XᵀX) β = Xᵀy`, solved by Gaussian elimination with partial pivoting.
/// Returns `None` if the system is singular. RMSE uses `n - p` degrees of
/// freedom (floored at 1).
fn fit_harmonic(t: &[f64], y: &[f64], t0: f64, period: f64, order: usize) -> Option<HarmonicFit> {
    let p = 2 + 2 * order;
    let n = t.len();
    if n < p {
        return None;
    }
    let mut ata = vec![vec![0.0f64; p]; p];
    let mut aty = vec![0.0f64; p];
    for i in 0..n {
        let x = design(t[i], t0, period, order);
        for a in 0..p {
            aty[a] += x[a] * y[i];
            for b in 0..p {
                ata[a][b] += x[a] * x[b];
            }
        }
    }
    let beta = solve(ata, aty)?;
    let mut sse = 0.0;
    for i in 0..n {
        let e = y[i] - predict(t[i], t0, period, order, &beta);
        sse += e * e;
    }
    let dof = (n as f64 - p as f64).max(1.0);
    Some(HarmonicFit {
        beta,
        rmse: (sse / dof).sqrt(),
    })
}

/// Solves a dense linear system `A x = b` (A is `n×n`, row-major as `Vec<Vec>`)
/// by Gaussian elimination with partial pivoting. Returns `None` if singular.
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        let mut piv = col;
        for r in (col + 1)..n {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-12 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        let d = a[col][col];
        let pivot_row = a[col].clone();
        let bcol = b[col];
        for r in (col + 1)..n {
            let f = a[r][col] / d;
            if f != 0.0 {
                for (ac, pc) in a[r][col..].iter_mut().zip(&pivot_row[col..]) {
                    *ac -= f * pc;
                }
                b[r] -= f * bcol;
            }
        }
    }
    let mut x = vec![0.0f64; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for c in (i + 1)..n {
            s -= a[i][c] * x[c];
        }
        x[i] = s / a[i][i];
    }
    Some(x)
}

/// Loads the time-series stack into per-slice row-major `f64` buffers that share
/// a single template geometry. A comma-separated `paths` list reads `band` from
/// each raster; a single path reads every band of one multiband raster.
fn load_stack(
    paths: &[String],
    band_1based: usize,
) -> Result<(Raster, Vec<Vec<f64>>, f64), ToolError> {
    if paths.len() == 1 {
        let r = load_input_raster(&paths[0])?;
        if r.bands < 2 {
            return Err(ToolError::Validation(
                "a single-raster input must be multiband (each band is a time slice); pass a comma-separated list otherwise".to_string(),
            ));
        }
        let nodata = r.nodata;
        let slices: Vec<Vec<f64>> = (0..r.bands).map(|b| band_to_vec(&r, b as isize)).collect();
        Ok((r, slices, nodata))
    } else {
        let band = (band_1based - 1) as isize;
        let rasters: Vec<Raster> = paths
            .iter()
            .map(|p| load_input_raster(p))
            .collect::<Result<_, _>>()?;
        let (rows, cols) = (rasters[0].rows, rasters[0].cols);
        for (i, r) in rasters.iter().enumerate() {
            if r.rows != rows || r.cols != cols {
                return Err(ToolError::Validation(format!(
                    "raster {i} is {}x{}, expected {rows}x{cols}",
                    r.rows, r.cols
                )));
            }
            if band < 0 || band as usize >= r.bands {
                return Err(ToolError::Validation(format!(
                    "band {band_1based} out of range for raster {i}"
                )));
            }
        }
        let nodata = rasters[0].nodata;
        let slices: Vec<Vec<f64>> = rasters.iter().map(|r| band_to_vec(r, band)).collect();
        let template = rasters.into_iter().next().unwrap();
        Ok((template, slices, nodata))
    }
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    dates: Option<Vec<f64>>,
    period: f64,
    harmonic_order: usize,
    change_threshold: f64,
    min_consecutive: usize,
    min_observations: usize,
    band: usize,
}

fn parse_inputs(args: &ToolArgs) -> Result<Vec<String>, ToolError> {
    let s = args
        .get("input")
        .and_then(Value::as_str)
        .ok_or_else(|| ToolError::Validation("missing required parameter 'input'".to_string()))?;
    let paths: Vec<String> = s
        .split(',')
        .map(str::trim)
        .filter(|p| !p.is_empty())
        .map(String::from)
        .collect();
    if paths.is_empty() {
        return Err(ToolError::Validation("'input' is empty".to_string()));
    }
    Ok(paths)
}

fn opt_f64(args: &ToolArgs, key: &str, default: f64) -> Result<f64, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(default),
        Some(Value::Number(n)) => Ok(n.as_f64().unwrap_or(default)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(default),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation(format!("'{key}' must be a number"))),
        Some(_) => Err(ToolError::Validation(format!("'{key}' must be a number"))),
    }
}

fn opt_usize(args: &ToolArgs, key: &str, default: usize, min: usize) -> Result<usize, ToolError> {
    let v = match args.get(key) {
        None | Some(Value::Null) => default,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(default as u64) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => default,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation(format!("'{key}' must be an integer")))?,
        Some(_) => return Err(ToolError::Validation(format!("'{key}' must be an integer"))),
    };
    Ok(v.max(min))
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let dates = match args.get("dates").and_then(Value::as_str) {
        None => None,
        Some(s) if s.trim().is_empty() => None,
        Some(s) => Some(
            s.split(',')
                .map(str::trim)
                .filter(|x| !x.is_empty())
                .map(|x| {
                    x.parse::<f64>()
                        .map_err(|_| ToolError::Validation(format!("date '{x}' is not a number")))
                })
                .collect::<Result<Vec<_>, _>>()?,
        ),
    };
    let period = opt_f64(args, "period", 1.0)?;
    if period <= 0.0 {
        return Err(ToolError::Validation(
            "'period' must be positive".to_string(),
        ));
    }
    let change_threshold = opt_f64(args, "change_threshold", 3.0)?;
    if change_threshold <= 0.0 {
        return Err(ToolError::Validation(
            "'change_threshold' must be positive".to_string(),
        ));
    }
    Ok(Params {
        dates,
        period,
        harmonic_order: opt_usize(args, "harmonic_order", 1, 1)?,
        change_threshold,
        min_consecutive: opt_usize(args, "min_consecutive", 3, 1)?,
        min_observations: opt_usize(args, "min_observations", 12, 4)?,
        band: opt_usize(args, "band", 1, 1)?,
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

    /// Builds an in-memory multiband raster from `cols×rows` per-band buffers.
    fn multiband(cols: usize, rows: usize, bands: &[Vec<f64>]) -> String {
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
        for (b, buf) in bands.iter().enumerate() {
            for row in 0..rows {
                for col in 0..cols {
                    r.set(
                        b as isize,
                        row as isize,
                        col as isize,
                        buf[row * cols + col],
                    )
                    .unwrap();
                }
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    /// A synthetic monthly harmonic series (trend + annual season + small
    /// deterministic noise) with an optional level shift at a known date.
    fn harmonic_series(
        n: usize,
        start: f64,
        shift_at: Option<usize>,
        shift: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        let mut dates = Vec::new();
        let mut vals = Vec::new();
        for k in 0..n {
            let t = start + k as f64 / 12.0;
            // Deterministic pseudo-noise in [-0.02, 0.02] (real EO series are noisy).
            let h = (k as f64 * 12.9898 + 4.0).sin() * 43758.5453;
            let noise = ((h - h.floor()) - 0.5) * 0.04;
            let mut v =
                0.5 + 0.02 * (t - start) + 0.2 * (TAU * t).sin() + 0.1 * (TAU * t).cos() + noise;
            if let Some(s) = shift_at {
                if k >= s {
                    v += shift;
                }
            }
            dates.push(t);
            vals.push(v);
        }
        (dates, vals)
    }

    /// A pixel with a planted level shift at a known date reports exactly one
    /// break, dated within one monthly time-step of the plant.
    #[test]
    fn detects_planted_break() {
        let n = 48;
        let shift_k = 30; // 2018.0 + 30/12 = 2020.5
        let (dates, vals) = harmonic_series(n, 2018.0, Some(shift_k), 0.4);
        // 1x1 multiband raster: one band per observation.
        let bands: Vec<Vec<f64>> = vals.iter().map(|&v| vec![v]).collect();
        let stack = multiband(1, 1, &bands);
        let date_str = dates
            .iter()
            .map(|d| format!("{d}"))
            .collect::<Vec<_>>()
            .join(",");
        let args: ToolArgs = serde_json::from_value(json!({
            "input": stack, "dates": date_str,
        }))
        .unwrap();
        let out = AnalyzeChangesCcdcTool.run(&args, &ctx()).unwrap();
        let count = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        assert_eq!(count.get(0, 0, 0), 1.0, "exactly one break expected");

        // The detected break date is within one monthly step of the plant.
        let expected = 2018.0 + shift_k as f64 / 12.0;
        let prm = parse_params(&serde_json::from_value(json!({})).unwrap()).unwrap();
        let res = ccdc_pixel(&dates, &vals, &prm, prm.min_observations);
        assert_eq!(res.breaks.len(), 1);
        assert!(
            (res.breaks[0] - expected).abs() <= 1.0 / 12.0 + 1e-9,
            "break {} should be within one month of {expected}",
            res.breaks[0]
        );
    }

    /// A stable harmonic series (no shift) reports zero breaks.
    #[test]
    fn stable_pixel_zero_breaks() {
        let (dates, vals) = harmonic_series(48, 2018.0, None, 0.0);
        let bands: Vec<Vec<f64>> = vals.iter().map(|&v| vec![v]).collect();
        let stack = multiband(1, 1, &bands);
        let date_str = dates
            .iter()
            .map(|d| format!("{d}"))
            .collect::<Vec<_>>()
            .join(",");
        let args: ToolArgs = serde_json::from_value(json!({
            "input": stack, "dates": date_str,
        }))
        .unwrap();
        let out = AnalyzeChangesCcdcTool.run(&args, &ctx()).unwrap();
        let count = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        assert_eq!(count.get(0, 0, 0), 0.0, "stable series -> zero breaks");
    }

    /// The harmonic OLS solver recovers a known amplitude on a clean signal.
    #[test]
    fn harmonic_fit_recovers_amplitude() {
        let n = 24;
        let mut t = Vec::new();
        let mut y = Vec::new();
        for k in 0..n {
            let tt = 2000.0 + k as f64 / 12.0;
            t.push(tt);
            y.push(1.0 + 0.3 * (TAU * tt).sin()); // amplitude sqrt(0.3^2)=0.3
        }
        let fit = fit_harmonic(&t, &y, t[0], 1.0, 1).unwrap();
        let amp = (fit.beta[2] * fit.beta[2] + fit.beta[3] * fit.beta[3]).sqrt();
        assert!(
            (amp - 0.3).abs() < 1e-6,
            "amplitude should be ~0.3, got {amp}"
        );
        assert!(fit.rmse < 1e-6, "clean signal -> ~0 RMSE, got {}", fit.rmse);
    }

    /// Pixels with too few valid observations become no-data.
    #[test]
    fn insufficient_data_is_nodata() {
        // 8 slices but min_observations default 12 -> no-data.
        let (dates, vals) = harmonic_series(8, 2018.0, None, 0.0);
        let bands: Vec<Vec<f64>> = vals.iter().map(|&v| vec![v]).collect();
        let stack = multiband(1, 1, &bands);
        let date_str = dates
            .iter()
            .map(|d| format!("{d}"))
            .collect::<Vec<_>>()
            .join(",");
        let args: ToolArgs = serde_json::from_value(json!({
            "input": stack, "dates": date_str,
        }))
        .unwrap();
        let out = AnalyzeChangesCcdcTool.run(&args, &ctx()).unwrap();
        let count = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        assert_eq!(count.get(0, 0, 0), -9999.0, "too few obs -> no-data");
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            AnalyzeChangesCcdcTool.validate(&args)
        };
        assert!(bad(json!({})).is_err()); // missing input
        assert!(bad(json!({ "input": "a.tif,b.tif", "period": 0 })).is_err());
        assert!(bad(json!({ "input": "a.tif,b.tif", "change_threshold": -1 })).is_err());
        assert!(bad(json!({ "input": "a.tif,b.tif", "dates": "x,y" })).is_err());
        assert!(bad(json!({ "input": "a.tif,b.tif,c.tif" })).is_ok());
    }
}
