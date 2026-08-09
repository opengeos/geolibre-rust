//! GeoLibre tool: moving-window statistics along a raster cube's dimension.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Dimensional Moving Statistics*
//! (Image Analyst / Spatial Analyst).
//!
//! ## Why the catalog needs it
//!
//! Smoothing a per-pixel time series in place is how a noisy satellite record
//! becomes usable: a 5-slice rolling mean removes sensor noise from an NDVI
//! series without shortening it, and a rolling median removes undetected cloud
//! outright while leaving genuine phenology alone.
//!
//! Nothing in either registry does it on rasters. `time_series_smoothing`
//! (round 17) is vector-only — it groups *features* by an id field. The bundled
//! focal filters (`gaussian_filter`, `percentile_filter`, `total_filter`) move
//! a window *spatially within one slice*, which is a different axis entirely.
//! `aggregate_multidimensional_raster` reduces slices into bins and so
//! shortens the cube; this preserves one output slice per input slice.
//!
//! ## Circular mean
//!
//! `circular_mean` exists because averaging directional data — aspect, wind
//! direction, wave heading, interferometric phase — with an arithmetic mean is
//! simply wrong: 350 degrees and 10 degrees average to 180, the exact opposite
//! of the correct 0. The circular form averages the unit vectors instead. Use
//! `period` to declare the wrap (360 for degrees, the default; 2*pi for
//! radians).
//!
//! ## Window and edges
//!
//! The window spans `backward_window` slices before and `forward_window` after
//! the slice being written, inclusive of that slice. Near the ends of the cube
//! the window is truncated rather than padded, so the statistic is always
//! computed from real observations; `min_valid` sets how many are required
//! before a value is written at all.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::args_common::{choice_or, f64_or, usize_or};
use crate::common::{parse_optional_output, write_or_store_output};
use crate::cube::load_cube;
use crate::raster_stack::raster_like_multiband;

pub struct DimensionalMovingStatisticsTool;

impl Tool for DimensionalMovingStatisticsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "dimensional_moving_statistics",
            display_name: "Dimensional Moving Statistics",
            summary: "Runs a moving window along a raster cube's dimension, writing one smoothed slice per input slice — the in-place rolling mean or median that turns a noisy satellite time series into a usable one (ArcGIS Dimensional Moving Statistics). Round 17's time_series_smoothing is vector-only, the bundled focal filters move a window spatially within a single slice, and aggregate_multidimensional_raster shortens the cube into bins. Includes a circular mean, without which averaging aspect or phase across the wrap point is meaningless.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "One multiband raster (each band is a slice) or a comma-separated list of co-registered rasters.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output multiband raster with one band per input slice. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "dimension",
                    description: "Name of the dimension, used only in the report (default 'slice').",
                    required: false,
                },
                ToolParamSpec {
                    name: "backward_window",
                    description: "Slices before the current one to include (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "forward_window",
                    description: "Slices after the current one to include (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "statistic",
                    description: "'mean' (default), 'circular_mean', 'majority', 'maximum', 'median', 'minimum', 'minority', 'percentile', 'std', or 'sum'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "percentile_value",
                    description: "Percentile in [0, 100] for statistic 'percentile' (default 50).",
                    required: false,
                },
                ToolParamSpec {
                    name: "period",
                    description: "Wrap period for 'circular_mean' (default 360, i.e. degrees; use 6.283185307179586 for radians).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_valid",
                    description: "Minimum non-no-data observations in the window before a value is written (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "nodata_handling",
                    description: "'data' (default; skip no-data and compute from the rest), 'nodata' (any no-data in the window makes the output no-data), or 'fill_nodata' (also write a value where the centre slice itself is no-data).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        crate::raster_stack::parse_input_paths(args, "input")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let prm = parse_params(args)?;
        let output = parse_optional_output(args, "output")?;
        let cube = load_cube(args, "input", "dimension_values", "dimension", 1)?;

        let (rows, cols, n) = (cube.rows, cube.cols, cube.len());
        ctx.progress.info(&format!(
            "{n} slice(s) over '{}', window -{}/+{}, statistic {}",
            cube.dimension,
            prm.backward,
            prm.forward,
            prm.statistic.label()
        ));

        let nodata = -9999.0_f64;
        let mut bands: Vec<Vec<f64>> = vec![vec![nodata; rows * cols]; n];
        let mut series: Vec<Option<f64>> = Vec::with_capacity(n);
        let mut window: Vec<f64> = Vec::with_capacity(n);
        let mut written = 0usize;

        for r in 0..rows {
            for c in 0..cols {
                cube.series(r, c, &mut series);
                for s in 0..n {
                    // `fill_nodata` is the only mode that writes where the
                    // centre slice itself has no observation; the other two
                    // leave a gap where the input had one.
                    if series[s].is_none() && prm.fill != NoDataHandling::Fill {
                        continue;
                    }
                    let lo = s.saturating_sub(prm.backward);
                    let hi = (s + prm.forward).min(n - 1);

                    window.clear();
                    let mut had_nodata = false;
                    for item in series.iter().take(hi + 1).skip(lo) {
                        match item {
                            Some(v) => window.push(*v),
                            None => had_nodata = true,
                        }
                    }
                    if prm.fill == NoDataHandling::Strict && had_nodata {
                        continue;
                    }
                    if window.len() < prm.min_valid {
                        continue;
                    }
                    bands[s][r * cols + c] = reduce(&mut window, &prm);
                    written += 1;
                }
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let out_raster = raster_like_multiband(cube.template(), &bands, nodata, DataType::F32)?;
        let out_path = write_or_store_output(out_raster, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("dimension".to_string(), json!(cube.dimension));
        outputs.insert("slices".to_string(), json!(n));
        outputs.insert("statistic".to_string(), json!(prm.statistic.label()));
        outputs.insert("backward_window".to_string(), json!(prm.backward));
        outputs.insert("forward_window".to_string(), json!(prm.forward));
        outputs.insert("values_written".to_string(), json!(written));
        Ok(ToolRunResult { outputs })
    }
}

/// Reduces one window. `vals` may be reordered.
fn reduce(vals: &mut [f64], prm: &Params) -> f64 {
    match prm.statistic {
        Statistic::Mean => vals.iter().sum::<f64>() / vals.len() as f64,
        Statistic::CircularMean => circular_mean(vals, prm.period),
        Statistic::Sum => vals.iter().sum(),
        Statistic::Maximum => vals.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        Statistic::Minimum => vals.iter().copied().fold(f64::INFINITY, f64::min),
        Statistic::Std => {
            let n = vals.len() as f64;
            let mean = vals.iter().sum::<f64>() / n;
            (vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n).sqrt()
        }
        Statistic::Median => percentile(vals, 50.0),
        Statistic::Percentile => percentile(vals, prm.percentile_value),
        Statistic::Majority => extreme_frequency(vals, true),
        Statistic::Minority => extreme_frequency(vals, false),
    }
}

/// Mean direction of angular data, returned in `[0, period)`.
///
/// Averages the unit vectors rather than the raw numbers, so values straddling
/// the wrap point combine correctly.
fn circular_mean(vals: &[f64], period: f64) -> f64 {
    let k = 2.0 * PI / period;
    let (mut sx, mut sy) = (0.0, 0.0);
    for &v in vals {
        let a = v * k;
        sx += a.cos();
        sy += a.sin();
    }
    // A perfectly opposed pair has no mean direction; report the first value
    // rather than an arbitrary angle from atan2(0, 0).
    if sx.abs() < 1e-12 && sy.abs() < 1e-12 {
        return vals[0].rem_euclid(period);
    }
    (sy.atan2(sx) / k).rem_euclid(period)
}

/// Linear-interpolated percentile, matching `cell_statistics`.
fn percentile(vals: &mut [f64], p: f64) -> f64 {
    vals.sort_by(f64::total_cmp);
    let n = vals.len();
    if n == 1 {
        return vals[0];
    }
    let rank = (p / 100.0) * (n - 1) as f64;
    let lo = rank.floor() as usize;
    let hi = rank.ceil() as usize;
    if lo == hi {
        vals[lo]
    } else {
        let frac = rank - lo as f64;
        vals[lo] * (1.0 - frac) + vals[hi] * frac
    }
}

/// Most or least frequent value; ties resolve to the smaller value.
fn extreme_frequency(vals: &[f64], want_max: bool) -> f64 {
    let mut counts: BTreeMap<u64, (usize, f64)> = BTreeMap::new();
    for &v in vals {
        let e = counts.entry(v.to_bits()).or_insert((0, v));
        e.0 += 1;
    }
    let mut ordered: Vec<(f64, usize)> = counts.values().map(|(c, v)| (*v, *c)).collect();
    ordered.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut best_count = if want_max { 0usize } else { usize::MAX };
    let mut best_val = f64::NAN;
    for (v, count) in ordered {
        let better = if want_max {
            count > best_count
        } else {
            count < best_count
        };
        if better {
            best_count = count;
            best_val = v;
        }
    }
    best_val
}

// ── Parameters ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Statistic {
    Mean,
    CircularMean,
    Majority,
    Maximum,
    Median,
    Minimum,
    Minority,
    Percentile,
    Std,
    Sum,
}

impl Statistic {
    fn label(self) -> &'static str {
        match self {
            Statistic::Mean => "mean",
            Statistic::CircularMean => "circular_mean",
            Statistic::Majority => "majority",
            Statistic::Maximum => "maximum",
            Statistic::Median => "median",
            Statistic::Minimum => "minimum",
            Statistic::Minority => "minority",
            Statistic::Percentile => "percentile",
            Statistic::Std => "std",
            Statistic::Sum => "sum",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum NoDataHandling {
    /// Skip no-data in the window; leave gaps where the centre slice is empty.
    Data,
    /// Any no-data in the window makes the output no-data.
    Strict,
    /// Also write where the centre slice itself is no-data.
    Fill,
}

struct Params {
    backward: usize,
    forward: usize,
    statistic: Statistic,
    percentile_value: f64,
    period: f64,
    min_valid: usize,
    fill: NoDataHandling,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let backward = usize_or(args, "backward_window", 1)?;
    let forward = usize_or(args, "forward_window", 1)?;
    if backward == 0 && forward == 0 {
        return Err(ToolError::Validation(
            "'backward_window' and 'forward_window' cannot both be 0; the window would be the \
             single centre slice and the tool a no-op"
                .to_string(),
        ));
    }

    let statistic = match choice_or(
        args,
        "statistic",
        &[
            "mean",
            "circular_mean",
            "majority",
            "maximum",
            "median",
            "minimum",
            "minority",
            "percentile",
            "std",
            "sum",
        ],
        "mean",
    )? {
        "circular_mean" => Statistic::CircularMean,
        "majority" => Statistic::Majority,
        "maximum" => Statistic::Maximum,
        "median" => Statistic::Median,
        "minimum" => Statistic::Minimum,
        "minority" => Statistic::Minority,
        "percentile" => Statistic::Percentile,
        "std" => Statistic::Std,
        "sum" => Statistic::Sum,
        _ => Statistic::Mean,
    };

    let percentile_value = f64_or(args, "percentile_value", 50.0)?;
    if !(0.0..=100.0).contains(&percentile_value) {
        return Err(ToolError::Validation(format!(
            "'percentile_value' must be in [0, 100], got {percentile_value}"
        )));
    }

    let period = f64_or(args, "period", 360.0)?;
    if !period.is_finite() || period <= 0.0 {
        return Err(ToolError::Validation(format!(
            "'period' must be positive, got {period}"
        )));
    }

    let min_valid = usize_or(args, "min_valid", 1)?.max(1);
    let fill = match choice_or(
        args,
        "nodata_handling",
        &["data", "nodata", "fill_nodata"],
        "data",
    )? {
        "nodata" => NoDataHandling::Strict,
        "fill_nodata" => NoDataHandling::Fill,
        _ => NoDataHandling::Data,
    };

    Ok(Params {
        backward,
        forward,
        statistic,
        percentile_value,
        period,
        min_valid,
        fill,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::load_input_raster;
    use crate::cube::test_support::cube_raster;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::Raster;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn run(args: Value) -> (Raster, BTreeMap<String, Value>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = DimensionalMovingStatisticsTool.run(&args, &ctx()).unwrap();
        let r = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (r, out.outputs)
    }

    fn series_cube(vals: &[f64]) -> String {
        let slices: Vec<Vec<f64>> = vals.iter().map(|&v| vec![v]).collect();
        cube_raster(1, 1, &slices)
    }

    /// A centred 3-slice mean, with the window truncated at both ends.
    #[test]
    fn centred_rolling_mean() {
        let (out, outputs) = run(json!({ "input": series_cube(&[1.0, 2.0, 3.0, 4.0, 5.0]) }));
        assert_eq!(out.bands, 5, "the cube length must be preserved");
        assert_eq!(outputs["slices"].as_u64().unwrap(), 5);
        // Edges use a truncated window rather than padding.
        assert!((out.get(0, 0, 0) - 1.5).abs() < 1e-6, "mean(1,2)");
        assert!((out.get(1, 0, 0) - 2.0).abs() < 1e-6, "mean(1,2,3)");
        assert!((out.get(2, 0, 0) - 3.0).abs() < 1e-6, "mean(2,3,4)");
        assert!((out.get(3, 0, 0) - 4.0).abs() < 1e-6, "mean(3,4,5)");
        assert!((out.get(4, 0, 0) - 4.5).abs() < 1e-6, "mean(4,5)");
    }

    /// Asymmetric windows: a trailing average uses only the past.
    #[test]
    fn backward_only_window_is_causal() {
        let (out, _) = run(json!({
            "input": series_cube(&[1.0, 2.0, 3.0, 4.0]),
            "backward_window": 2, "forward_window": 0
        }));
        assert!((out.get(0, 0, 0) - 1.0).abs() < 1e-6);
        assert!((out.get(1, 0, 0) - 1.5).abs() < 1e-6);
        assert!((out.get(2, 0, 0) - 2.0).abs() < 1e-6, "mean(1,2,3)");
        assert!((out.get(3, 0, 0) - 3.0).abs() < 1e-6, "mean(2,3,4)");
    }

    /// The reason the tool is useful: a rolling median removes a single-slice
    /// spike (an undetected cloud) that a mean would only smear.
    #[test]
    fn median_removes_a_spike_the_mean_smears() {
        // A flat series with one cloud spike at slice 3.
        let vals = [0.5, 0.5, 0.5, 9.0, 0.5, 0.5, 0.5];
        let median = run(json!({
            "input": series_cube(&vals), "statistic": "median"
        }))
        .0;
        for s in 0..7 {
            let v = median.get(s as isize, 0, 0);
            assert!(
                (v - 0.5).abs() < 1e-6,
                "median slice {s} should be 0.5, got {v}"
            );
        }
        let mean = run(json!({ "input": series_cube(&vals) })).0;
        // The mean spreads the spike over three slices instead of removing it.
        for s in [2, 3, 4] {
            assert!(
                mean.get(s, 0, 0) > 1.0,
                "mean should be contaminated at slice {s}"
            );
        }
    }

    /// The circular mean is why this parameter exists: an arithmetic mean of
    /// 350 and 10 degrees gives 180, the exact opposite of the right answer.
    #[test]
    fn circular_mean_handles_the_wrap() {
        // A symmetric pair straddling the wrap has an exact circular mean of 0.
        // Use a trailing window so slice 1's window is exactly {350, 10}.
        let (out, _) = run(json!({
            "input": series_cube(&[350.0, 10.0]),
            "statistic": "circular_mean", "backward_window": 1, "forward_window": 0
        }));
        let mid = out.get(1, 0, 0);
        let from_zero = mid.min(360.0 - mid);
        assert!(
            from_zero < 1e-3,
            "circular mean of 350 and 10 is exactly 0/360, got {mid}"
        );

        // The same window under the arithmetic mean gives 180 — the exact
        // opposite direction. That is the error this statistic exists to avoid.
        let plain = run(json!({
            "input": series_cube(&[350.0, 10.0]),
            "backward_window": 1, "forward_window": 0
        }))
        .0;
        assert!(
            (plain.get(1, 0, 0) - 180.0).abs() < 1e-3,
            "the arithmetic mean should give the wrong answer this test contrasts with, got {}",
            plain.get(1, 0, 0)
        );
    }

    /// Radians work through the `period` parameter.
    #[test]
    fn circular_mean_in_radians() {
        let two_pi = 2.0 * PI;
        let (out, _) = run(json!({
            "input": series_cube(&[two_pi - 0.1, 0.1]),
            "statistic": "circular_mean", "period": two_pi,
            "backward_window": 1, "forward_window": 0
        }));
        let mid = out.get(1, 0, 0);
        let from_zero = mid.min(two_pi - mid);
        assert!(from_zero < 1e-4, "expected exactly 0, got {mid}");
    }

    /// No-data modes behave differently at a gap.
    #[test]
    fn nodata_modes() {
        // Slice 1 is missing.
        let src = series_cube(&[1.0, -9999.0, 3.0]);

        // 'data': the gap stays a gap, neighbours are computed from what exists.
        let d = run(json!({ "input": src.clone() })).0;
        assert_eq!(d.get(1, 0, 0), -9999.0, "'data' must not fill the gap");
        assert!((d.get(0, 0, 0) - 1.0).abs() < 1e-6);

        // 'nodata': the gap poisons every window that touches it.
        let strict = run(json!({ "input": src.clone(), "nodata_handling": "nodata" })).0;
        for s in 0..3 {
            assert_eq!(
                strict.get(s, 0, 0),
                -9999.0,
                "slice {s} window touches the gap"
            );
        }

        // 'fill_nodata': the gap is filled from its neighbours.
        let fill = run(json!({ "input": src, "nodata_handling": "fill_nodata" })).0;
        assert!(
            (fill.get(1, 0, 0) - 2.0).abs() < 1e-6,
            "gap should be filled with mean(1,3) = 2"
        );
    }

    /// `min_valid` suppresses values computed from too little data.
    #[test]
    fn min_valid_suppresses_thin_windows() {
        let (out, _) = run(json!({
            "input": series_cube(&[1.0, 2.0, 3.0]), "min_valid": 3
        }));
        // Only the middle slice has a full 3-slice window.
        assert_eq!(out.get(0, 0, 0), -9999.0);
        assert!((out.get(1, 0, 0) - 2.0).abs() < 1e-6);
        assert_eq!(out.get(2, 0, 0), -9999.0);
    }

    /// Cells are smoothed independently.
    #[test]
    fn cells_are_independent() {
        let src = cube_raster(
            2,
            1,
            &[vec![1.0, 100.0], vec![2.0, 200.0], vec![3.0, 300.0]],
        );
        let (out, _) = run(json!({ "input": src }));
        assert!((out.get(1, 0, 0) - 2.0).abs() < 1e-6);
        assert!((out.get(1, 0, 1) - 200.0).abs() < 1e-6);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            DimensionalMovingStatisticsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        // A zero-width window would be a no-op rather than a smoother.
        assert!(
            bad(json!({"input": "a.tif", "backward_window": 0, "forward_window": 0})).is_err()
        );
        assert!(bad(json!({"input": "a.tif", "statistic": "mode"})).is_err());
        assert!(bad(json!({"input": "a.tif", "percentile_value": -1})).is_err());
        assert!(bad(json!({"input": "a.tif", "period": 0})).is_err());
        assert!(bad(json!({"input": "a.tif", "nodata_handling": "ignore"})).is_err());
        assert!(bad(json!({"input": "a.tif", "statistic": "circular_mean"})).is_ok());
    }
}
