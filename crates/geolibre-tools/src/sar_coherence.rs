//! GeoLibre tool: interferometric coherence and phase between SAR acquisitions.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Compute Coherence* (Image Analyst),
//! also covering *Generate Interferogram* through the optional phase output.
//!
//! Coherence is the workhorse SAR change-detection product: it measures how
//! similar the *scattering geometry* stayed between two passes, independent of
//! brightness. Because it is a phase-stability measure rather than an intensity
//! measure, it detects changes intensity differencing cannot see at all —
//! building collapse, flooding under vegetation, deforestation, landslides,
//! crop harvest, ice motion — in all weather, day or night.
//!
//! Nothing in either registry computes it. GeoLibre's change-detection tools
//! (`detect_feature_changes`, `detect_image_anomalies`, `landtrendr`,
//! `analyze_changes_ccdc`, `image_regression`) are all optical/intensity-domain
//! and none uses phase; the bundled speckle filters work on detected intensity,
//! which has already discarded the phase this depends on.
//!
//! Complex input follows the two-band I/Q convention established by `multilook`.
//!
//! ## Scope, deliberately
//!
//! The inputs are assumed **already co-registered**. Co-registration, orbit
//! correction and phase unwrapping are separate tools with much heavier
//! requirements (orbit state vectors, external metadata); bundling them here
//! would make this unshippable.
//!
//! ## The correctness trap
//!
//! The numerator is the **magnitude of the complex sum**, not the sum of
//! magnitudes. Those differ by exactly the phase stability being measured, so
//! confusing them yields a coherence map that is ≈1 everywhere.
//! `random_phase_decorrelates` and `stable_phase_stays_coherent` pin this down.

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

/// Computes normalised interferometric coherence between two co-registered
/// complex SAR acquisitions.
pub struct SarCoherenceTool;

impl Tool for SarCoherenceTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "sar_coherence",
            display_name: "SAR Coherence",
            summary: "Computes normalised interferometric coherence magnitude (and optionally the interferometric phase) between two co-registered complex SAR acquisitions over a moving window (ArcGIS Compute Coherence / Generate Interferogram). Coherence measures phase stability rather than brightness, so it detects change that intensity differencing cannot see — collapse, flooding under vegetation, harvest, ice motion — in all weather. No tool in either registry uses SAR phase. Inputs are two-band I/Q pairs and must already be co-registered.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "reference",
                    description: "Complex SAR raster (two bands: I and Q) from the first acquisition.",
                    required: true,
                },
                ToolParamSpec {
                    name: "secondary",
                    description: "Co-registered complex SAR raster (two bands: I and Q) from the second acquisition.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output coherence raster path (magnitude, 0..1). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_phase",
                    description: "Optional output path for the interferometric phase raster in radians, wrapped to (-pi, pi].",
                    required: false,
                },
                ToolParamSpec {
                    name: "window_size",
                    description: "Estimation window as 'range,azimuth' in cells (default '5,5'). The estimator is centred, so an even value is rounded UP to the next odd size; the effective window is reported back as window_range/window_azimuth.",
                    required: false,
                },
                ToolParamSpec {
                    name: "bias_correction",
                    description: "Correct the upward bias of the coherence estimator at small window sizes (default true).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["reference", "secondary"] {
            if args
                .get(key)
                .and_then(Value::as_str)
                .map(str::trim)
                .unwrap_or("")
                .is_empty()
            {
                return Err(ToolError::Validation(format!(
                    "missing required string parameter '{key}'"
                )));
            }
        }
        parse_window(args)?;
        opt_bool(args, "bias_correction")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let ref_path = required_str(args, "reference")?;
        let sec_path = required_str(args, "secondary")?;
        let output = parse_optional_output(args, "output")?;
        let output_phase = parse_optional_output(args, "output_phase")?;
        let (req_r, req_a) = parse_window(args)?;
        // The estimator uses an inclusive -half..=half range, so the realised
        // window is always odd. Normalise up front and report the EFFECTIVE
        // size rather than echoing an even request the run never honoured.
        let (win_r, win_a) = (req_r | 1, req_a | 1);
        let bias_correction = opt_bool(args, "bias_correction")?.unwrap_or(true);

        let reference = load_input_raster(ref_path)?;
        let secondary = load_input_raster(sec_path)?;

        for (name, r) in [("reference", &reference), ("secondary", &secondary)] {
            if r.bands < 2 {
                return Err(ToolError::Validation(format!(
                    "'{name}' must be a complex two-band I/Q raster; got {} band(s)",
                    r.bands
                )));
            }
        }
        if reference.rows != secondary.rows || reference.cols != secondary.cols {
            return Err(ToolError::Validation(format!(
                "reference is {}x{} but secondary is {}x{}; inputs must be co-registered onto the same grid",
                reference.rows, reference.cols, secondary.rows, secondary.cols
            )));
        }

        let rows = reference.rows;
        let cols = reference.cols;

        ctx.progress.info("reading complex samples");
        // Pull both scenes into flat I/Q buffers, marking invalid samples NaN.
        let (ri, rq) = complex_bands(&reference);
        let (si, sq) = complex_bands(&secondary);

        ctx.progress.info(&format!(
            "estimating coherence over a {win_r} x {win_a} window"
        ));

        let out_nodata = -1.0_f64;
        let est = estimate_coherence(
            &ri,
            &rq,
            &si,
            &sq,
            rows,
            cols,
            win_r,
            win_a,
            bias_correction,
        );
        ctx.progress.progress(1.0);
        let (valid, coh_sum) = (est.valid, est.sum);
        // NaN marks "not estimated"; the raster wants the -1.0 sentinel.
        let coh: Vec<f64> = est
            .coherence
            .iter()
            .map(|v| if v.is_nan() { out_nodata } else { *v })
            .collect();
        let phase = est.phase;

        let coh_raster = raster_like_with_data(&reference, coh, out_nodata, DataType::F32)?;
        let out_path = write_or_store_output(coh_raster, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("valid_cells".to_string(), json!(valid));
        outputs.insert(
            "mean_coherence".to_string(),
            json!(if valid > 0 {
                coh_sum / valid as f64
            } else {
                0.0
            }),
        );
        outputs.insert("window_range".to_string(), json!(win_r));
        outputs.insert("window_azimuth".to_string(), json!(win_a));
        outputs.insert("window_range_requested".to_string(), json!(req_r));
        outputs.insert("window_azimuth_requested".to_string(), json!(req_a));

        // Phase is always produced (stored in memory when no path is given),
        // matching how create_overpass handles its secondary output. Phase is
        // wrapped to (-pi, pi], so -1 is a legitimate value; use a sentinel
        // outside that range for no-data.
        let phase_nodata = -9999.0_f64;
        let buf: Vec<f64> = phase
            .iter()
            .map(|v| if v.is_nan() { phase_nodata } else { *v })
            .collect();
        let phase_raster = raster_like_with_data(&reference, buf, phase_nodata, DataType::F32)?;
        let phase_path = write_or_store_output(phase_raster, output_phase)?;
        outputs.insert("output_phase".to_string(), json!(phase_path));

        Ok(ToolRunResult { outputs })
    }
}

/// Reads bands 1 and 2 as I and Q, mapping no-data to NaN.
pub(crate) fn complex_bands(r: &wbraster::Raster) -> (Vec<f64>, Vec<f64>) {
    let rows = r.rows;
    let cols = r.cols;
    let nd = r.nodata;
    let mut i = vec![f64::NAN; rows * cols];
    let mut q = vec![f64::NAN; rows * cols];
    for row in 0..rows {
        for col in 0..cols {
            let a = r.get(0, row as isize, col as isize);
            let b = r.get(1, row as isize, col as isize);
            if a != nd && a.is_finite() && b != nd && b.is_finite() {
                i[row * cols + col] = a;
                q[row * cols + col] = b;
            }
        }
    }
    (i, q)
}

/// The result of [`estimate_coherence`].
///
/// Cells that could not be estimated are `NaN` in both `coherence` and `phase`;
/// callers map that onto whatever no-data sentinel their raster uses.
pub(crate) struct CoherenceEstimate {
    pub(crate) coherence: Vec<f64>,
    pub(crate) phase: Vec<f64>,
    pub(crate) valid: u64,
    pub(crate) sum: f64,
}

/// Per-cell coherence and interferometric phase over a moving window.
///
/// Factored out of `run` so `multitemporal_coherence` estimates every pair with
/// the *same* code rather than a copy. The correctness trap lives here: the
/// numerator is the **magnitude of the complex sum**, not the sum of
/// magnitudes — averaging magnitudes discards the phase alignment that
/// coherence measures and reports a decorrelated pair as fully coherent.
///
/// `ri`, `rq`, `si` and `sq` are row-major buffers of exactly `rows * cols`
/// samples; invalid samples are `NaN`. Cells that could not be estimated are
/// `NaN` in both outputs.
///
/// # Panics
///
/// Panics if any input buffer is shorter than `rows * cols`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn estimate_coherence(
    ri: &[f64],
    rq: &[f64],
    si: &[f64],
    sq: &[f64],
    rows: usize,
    cols: usize,
    win_r: usize,
    win_a: usize,
    bias_correction: bool,
) -> CoherenceEstimate {
    // Four same-typed slices and nine positional arguments make a transposed
    // call easy to write; without this the failure surfaces as an index panic
    // in the inner loop rather than at the boundary.
    let cells = rows * cols;
    for (name, buf) in [("ri", ri), ("rq", rq), ("si", si), ("sq", sq)] {
        assert!(
            buf.len() >= cells,
            "estimate_coherence: '{name}' holds {} samples but {rows}x{cols} needs {cells}",
            buf.len()
        );
    }

    let mut coherence = vec![f64::NAN; cells];
    let mut phase = vec![f64::NAN; cells];
    let mut valid = 0_u64;
    let mut sum = 0.0_f64;

    let half_r = (win_r / 2) as isize;
    let half_a = (win_a / 2) as isize;

    for row in 0..rows {
        for col in 0..cols {
            // Accumulate the complex cross-product and both intensities.
            let (mut cr, mut ci) = (0.0_f64, 0.0_f64);
            let (mut p_ref, mut p_sec) = (0.0_f64, 0.0_f64);
            let mut n = 0_usize;
            let mut poisoned = false;

            for dr in -half_a..=half_a {
                let r = row as isize + dr;
                if r < 0 || r >= rows as isize {
                    continue;
                }
                for dc in -half_r..=half_r {
                    let c = col as isize + dc;
                    if c < 0 || c >= cols as isize {
                        continue;
                    }
                    let k = r as usize * cols + c as usize;
                    let (a_i, a_q) = (ri[k], rq[k]);
                    let (b_i, b_q) = (si[k], sq[k]);
                    if !a_i.is_finite() || !a_q.is_finite() || !b_i.is_finite() || !b_q.is_finite()
                    {
                        poisoned = true;
                        continue;
                    }
                    // reference * conj(secondary)
                    cr += a_i * b_i + a_q * b_q;
                    ci += a_q * b_i - a_i * b_q;
                    p_ref += a_i * a_i + a_q * a_q;
                    p_sec += b_i * b_i + b_q * b_q;
                    n += 1;
                }
            }

            if n == 0 || poisoned {
                continue;
            }
            let denom = (p_ref * p_sec).sqrt();
            if denom <= 0.0 {
                continue;
            }
            // Magnitude of the complex sum — NOT the sum of magnitudes.
            let mut gamma = (cr * cr + ci * ci).sqrt() / denom;
            if bias_correction {
                gamma = correct_bias(gamma, n);
            }
            let gamma = gamma.clamp(0.0, 1.0);
            coherence[row * cols + col] = gamma;
            phase[row * cols + col] = ci.atan2(cr);
            valid += 1;
            sum += gamma;
        }
    }

    CoherenceEstimate {
        coherence,
        phase,
        valid,
        sum,
    }
}

/// Removes the well-known upward bias of the coherence magnitude estimator at
/// small sample counts. The estimator's expectation exceeds the true value by
/// roughly `(1 - gamma^2) / (2n)`; inverting that first-order relation is the
/// standard correction. Values stay clamped to `[0, 1]`.
fn correct_bias(gamma: f64, n: usize) -> f64 {
    if n <= 1 {
        return gamma;
    }
    let corrected = gamma - (1.0 - gamma * gamma) / (2.0 * n as f64);
    corrected.clamp(0.0, 1.0)
}

fn parse_window(args: &ToolArgs) -> Result<(usize, usize), ToolError> {
    let Some(s) = opt_str(args, "window_size")? else {
        return Ok((5, 5));
    };
    let parts: Vec<&str> = s
        .split(',')
        .map(str::trim)
        .filter(|p| !p.is_empty())
        .collect();
    let nums: Result<Vec<usize>, ToolError> = parts
        .iter()
        .map(|p| {
            p.parse::<usize>().map_err(|_| {
                ToolError::Validation(format!(
                    "parameter 'window_size' has non-integer component '{p}'"
                ))
            })
        })
        .collect();
    let nums = nums?;
    let (r, a) = match nums.len() {
        0 => (5, 5),
        1 => (nums[0], nums[0]),
        _ => (nums[0], nums[1]),
    };
    if r == 0 || a == 0 {
        return Err(ToolError::Validation(
            "'window_size' components must be at least 1".to_string(),
        ));
    }
    Ok((r, a))
}

fn opt_str<'a>(args: &'a ToolArgs, key: &str) -> Result<Option<&'a str>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => Ok(Some(s.as_str())),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a string when provided"
        ))),
    }
}

fn opt_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
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

fn required_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a complex raster from per-pixel (I, Q) pairs in row-major order.
    fn complex_raster(cols: usize, rows: usize, iq: &[(f64, f64)]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 2,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: -9999.0,
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
                let (i, q) = iq[row * cols + col];
                r.set(0, row as isize, col as isize, i).unwrap();
                r.set(1, row as isize, col as isize, q).unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    /// A deterministic pseudo-random phase generator — no RNG, so WASM and
    /// native agree exactly.
    fn pseudo_phase(i: usize) -> f64 {
        let x = ((i as f64) * 12.9898).sin() * 43758.5453;
        (x - x.floor()) * std::f64::consts::TAU
    }

    fn run(reference: String, secondary: String, extra: Value) -> (Raster, ToolRunResult) {
        let mut obj = serde_json::Map::new();
        obj.insert("reference".to_string(), json!(reference));
        obj.insert("secondary".to_string(), json!(secondary));
        if let Value::Object(m) = extra {
            for (k, v) in m {
                obj.insert(k, v);
            }
        }
        let args: ToolArgs = serde_json::from_value(Value::Object(obj)).unwrap();
        let res = SarCoherenceTool.run(&args, &ctx()).unwrap();
        let r = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        (r, res)
    }

    /// Identical scenes are perfectly coherent.
    #[test]
    fn identical_scenes_are_fully_coherent() {
        let n = 49;
        let iq: Vec<(f64, f64)> = (0..n)
            .map(|i| {
                let p = pseudo_phase(i);
                (p.cos(), p.sin())
            })
            .collect();
        let a = complex_raster(7, 7, &iq);
        let b = complex_raster(7, 7, &iq);
        let (out, _) = run(a, b, json!({ "bias_correction": false }));
        assert!(
            (out.get(0, 3, 3) - 1.0).abs() < 1e-9,
            "identical scenes must give coherence 1, got {}",
            out.get(0, 3, 3)
        );
    }

    /// A constant phase *shift* between the scenes preserves coherence — the
    /// measure is phase-stability, not phase-equality. This is what makes it a
    /// deformation-insensitive change detector.
    #[test]
    fn stable_phase_stays_coherent() {
        let n = 49;
        let shift = 0.7_f64;
        let a_iq: Vec<(f64, f64)> = (0..n)
            .map(|i| {
                let p = pseudo_phase(i);
                (p.cos(), p.sin())
            })
            .collect();
        let b_iq: Vec<(f64, f64)> = (0..n)
            .map(|i| {
                let p = pseudo_phase(i) + shift;
                (p.cos(), p.sin())
            })
            .collect();
        let a = complex_raster(7, 7, &a_iq);
        let b = complex_raster(7, 7, &b_iq);
        let (out, _) = run(a, b, json!({ "bias_correction": false }));
        assert!(
            (out.get(0, 3, 3) - 1.0).abs() < 1e-9,
            "a constant phase shift must not decorrelate, got {}",
            out.get(0, 3, 3)
        );
    }

    /// Random relative phase decorrelates. If the implementation summed
    /// magnitudes instead of the complex values, this would still read ~1 —
    /// which is precisely the bug this test exists to catch.
    #[test]
    fn random_phase_decorrelates() {
        let n = 81;
        let a_iq: Vec<(f64, f64)> = (0..n)
            .map(|i| {
                let p = pseudo_phase(i);
                (p.cos(), p.sin())
            })
            .collect();
        // Independent phase in the secondary scene.
        let b_iq: Vec<(f64, f64)> = (0..n)
            .map(|i| {
                let p = pseudo_phase(i + 5000);
                (p.cos(), p.sin())
            })
            .collect();
        let a = complex_raster(9, 9, &a_iq);
        let b = complex_raster(9, 9, &b_iq);
        let (out, res) = run(
            a,
            b,
            json!({ "window_size": "9,9", "bias_correction": false }),
        );
        let g = out.get(0, 4, 4);
        assert!(
            g < 0.6,
            "random relative phase must decorrelate; got {g} (a value near 1 means magnitudes were summed instead of complex values)"
        );
        assert!(res.outputs["mean_coherence"].as_f64().unwrap() < 0.8);
    }

    /// The phase output is the interferogram, wrapped to (-pi, pi].
    #[test]
    fn phase_output_recovers_the_shift() {
        let n = 25;
        let shift = 1.1_f64;
        let a_iq: Vec<(f64, f64)> = (0..n).map(|_| (1.0, 0.0)).collect();
        let b_iq: Vec<(f64, f64)> = (0..n).map(|_| (shift.cos(), shift.sin())).collect();
        let a = complex_raster(5, 5, &a_iq);
        let b = complex_raster(5, 5, &b_iq);

        let args: ToolArgs = serde_json::from_value(json!({
            "reference": a,
            "secondary": b,
            "output_phase": ""
        }))
        .unwrap();
        let res = SarCoherenceTool.run(&args, &ctx()).unwrap();
        let ph = load_input_raster(res.outputs["output_phase"].as_str().unwrap()).unwrap();
        // reference * conj(secondary) has phase -shift.
        assert!(
            (ph.get(0, 2, 2) + shift).abs() < 1e-6,
            "expected phase {}, got {}",
            -shift,
            ph.get(0, 2, 2)
        );
        assert!(ph.get(0, 2, 2) > -std::f64::consts::PI);
        assert!(ph.get(0, 2, 2) <= std::f64::consts::PI);
    }

    /// An even window request is rounded up to the next odd size and the
    /// effective value is reported, rather than echoing a size never used.
    #[test]
    fn even_window_is_normalised_and_reported() {
        let iq: Vec<(f64, f64)> = (0..49)
            .map(|i| {
                let p = pseudo_phase(i);
                (p.cos(), p.sin())
            })
            .collect();
        let a = complex_raster(7, 7, &iq);
        let b = complex_raster(7, 7, &iq);
        let (_, res) = run(a, b, json!({ "window_size": "4,4" }));
        assert_eq!(res.outputs["window_range"], json!(5));
        assert_eq!(res.outputs["window_azimuth"], json!(5));
        assert_eq!(res.outputs["window_range_requested"], json!(4));
    }

    /// A phase of exactly -1.0 rad is a real value, not a no-data marker.
    #[test]
    fn phase_of_minus_one_radian_is_preserved() {
        let shift = 1.0_f64; // reference * conj(secondary) then has phase -1.0
        let a_iq: Vec<(f64, f64)> = (0..25).map(|_| (1.0, 0.0)).collect();
        let b_iq: Vec<(f64, f64)> = (0..25).map(|_| (shift.cos(), shift.sin())).collect();
        let a = complex_raster(5, 5, &a_iq);
        let b = complex_raster(5, 5, &b_iq);
        let args: ToolArgs =
            serde_json::from_value(json!({ "reference": a, "secondary": b })).unwrap();
        let res = SarCoherenceTool.run(&args, &ctx()).unwrap();
        let ph = load_input_raster(res.outputs["output_phase"].as_str().unwrap()).unwrap();
        let v = ph.get(0, 2, 2);
        assert!(
            (v + 1.0).abs() < 1e-6 && v != ph.nodata,
            "-1.0 rad must survive as a real phase, got {v} (nodata {})",
            ph.nodata
        );
    }

    /// Bias correction lowers the estimate at small window sizes, where the
    /// magnitude estimator is known to read high.
    #[test]
    fn bias_correction_lowers_small_window_estimates() {
        let n = 25;
        let a_iq: Vec<(f64, f64)> = (0..n)
            .map(|i| {
                let p = pseudo_phase(i);
                (p.cos(), p.sin())
            })
            .collect();
        let b_iq: Vec<(f64, f64)> = (0..n)
            .map(|i| {
                let p = pseudo_phase(i) + pseudo_phase(i + 99) * 0.8;
                (p.cos(), p.sin())
            })
            .collect();
        let a = complex_raster(5, 5, &a_iq);
        let b = complex_raster(5, 5, &b_iq);
        let (raw, _) = run(a.clone(), b.clone(), json!({ "bias_correction": false }));
        let (corrected, _) = run(a, b, json!({ "bias_correction": true }));
        assert!(
            corrected.get(0, 2, 2) <= raw.get(0, 2, 2),
            "bias correction must not increase coherence"
        );
    }

    #[test]
    fn rejects_bad_parameters() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(SarCoherenceTool.validate(&args).is_err());

        let a = complex_raster(2, 2, &[(1.0, 0.0); 4]);
        let b = complex_raster(2, 2, &[(1.0, 0.0); 4]);
        for bad in [
            json!({ "reference": a.clone(), "secondary": b.clone(), "window_size": "x" }),
            json!({ "reference": a.clone(), "secondary": b.clone(), "window_size": "nan" }),
            json!({ "reference": a.clone(), "secondary": b.clone(), "window_size": "0,0" }),
        ] {
            let args: ToolArgs = serde_json::from_value(bad).unwrap();
            assert!(SarCoherenceTool.validate(&args).is_err());
        }

        // Mismatched grids are rejected: the inputs must be co-registered.
        let big = complex_raster(3, 3, &[(1.0, 0.0); 9]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "reference": a, "secondary": big })).unwrap();
        assert!(SarCoherenceTool.run(&args, &ctx()).is_err());
    }
}
