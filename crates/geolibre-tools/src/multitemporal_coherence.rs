//! GeoLibre tool: InSAR coherence statistics across an N-image stack.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Generate Multitemporal Coherence*
//! (Image Analyst).
//!
//! ## Why pairwise coherence was not enough
//!
//! `sar_coherence` (round 16) is scoped in its own module doc to **two**
//! co-registered acquisitions. That is the right primitive but not how
//! coherence is used: real InSAR change detection runs over a time series, and
//! the products that matter are stack-level.
//!
//! * **Mean coherence** separates persistent scatterers (buildings, rock) from
//!   decorrelating surfaces (vegetation, water) — the standard first pass for
//!   PS/DS candidate selection.
//! * **Coherence change** between consecutive pairs localizes an event *in
//!   time*: building collapse, flooding under canopy, harvest, landslide.
//! * **Decorrelation slope** — how fast coherence decays with temporal baseline
//!   — is a land-cover discriminator in its own right.
//!
//! None of these was reachable by running `sar_coherence` pairwise by hand,
//! because nothing reduced the resulting pairs into a stack statistic.
//!
//! ## Shared estimator, not a copy
//!
//! Every pair goes through `sar_coherence::estimate_coherence`, the same code
//! the pairwise tool uses. Copying that loop would have duplicated its central
//! correctness trap — the numerator is the magnitude of the complex sum, not
//! the sum of magnitudes.
//!
//! ## Scope
//!
//! Inputs are assumed **already co-registered**, exactly as `sar_coherence`
//! assumes; `coregister_rasters` (round 18) is the upstream tool. Each input is
//! a two-band I/Q complex raster, the convention `multilook` and `sar_coherence`
//! established.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::args_common::{bool_or, choice_or, usize_or};
use crate::common::{load_input_raster, parse_optional_output, write_or_store_output};
use crate::raster_stack::{check_alignment_refs, raster_like_multiband};
use crate::sar_coherence::{complex_bands, estimate_coherence};

const PAIRS: [&str; 3] = ["consecutive", "all", "reference"];
const STATS: [&str; 5] = ["mean", "min", "max", "std", "decorrelation_slope"];

/// Guard on the `all` pairing: N images give N(N-1)/2 pairs, and each pair is a
/// full windowed pass over the grid. Without a cap a 60-image stack would
/// quietly queue 1,770 passes.
const MAX_PAIRS: usize = 512;

pub struct MultitemporalCoherenceTool;

impl Tool for MultitemporalCoherenceTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "multitemporal_coherence",
            display_name: "Multitemporal Coherence",
            summary: "Computes interferometric coherence across a stack of co-registered SAR acquisitions and reduces it to mean/min/max/std and temporal-decorrelation slope, optionally emitting the per-pair coherence cube (ArcGIS Generate Multitemporal Coherence). sar_coherence handles exactly two scenes, so persistent-scatterer selection and coherence-change timing were out of reach.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "inputs",
                    description: "Comma- or semicolon-separated co-registered complex SLC rasters (two bands each: I then Q). At least two.",
                    required: true,
                },
                ToolParamSpec {
                    name: "dates",
                    description: "Optional acquisition times, one per input, as numbers in consistent units (e.g. days). Required by statistic 'decorrelation_slope'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "pairs",
                    description: "'consecutive' (default, i with i+1), 'all' (every unordered pair), or 'reference' (every scene against reference_index).",
                    required: false,
                },
                ToolParamSpec {
                    name: "reference_index",
                    description: "0-based index of the reference scene for pairs='reference' (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "window_size",
                    description: "Estimation window in cells; a single value or 'range,azimuth' (default 5). Even values are rounded up to odd.",
                    required: false,
                },
                ToolParamSpec {
                    name: "bias_correction",
                    description: "Apply the small-sample coherence bias correction (default true).",
                    required: false,
                },
                ToolParamSpec {
                    name: "statistics",
                    description: "Comma-separated: 'mean' (default), 'min', 'max', 'std', 'decorrelation_slope'. One output band per statistic, in this order.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output raster, one band per requested statistic. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_pairs",
                    description: "Optional raster receiving the per-pair coherence cube, one band per pair.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        let inputs = parse_list(args, "inputs")?;
        if inputs.len() < 2 {
            return Err(ToolError::Validation(format!(
                "'inputs' must list at least 2 acquisitions, got {}",
                inputs.len()
            )));
        }
        choice_or(args, "pairs", &PAIRS, "consecutive")?;
        bool_or(args, "bias_correction", true)?;
        let stats = parse_stats(args)?;
        let dates = parse_dates(args, inputs.len())?;
        if stats.iter().any(|s| s == "decorrelation_slope") && dates.is_none() {
            // Substituting the slice index for time would silently report a
            // slope per acquisition rather than per unit of temporal baseline.
            return Err(ToolError::Validation(
                "statistic 'decorrelation_slope' is a fit against temporal baseline, so 'dates' \
                 is required"
                    .to_string(),
            ));
        }
        if let Some(i) = args.get("reference_index").and_then(Value::as_u64) {
            if i as usize >= inputs.len() {
                return Err(ToolError::Validation(format!(
                    "'reference_index' {i} is out of range for {} input(s)",
                    inputs.len()
                )));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let inputs = parse_list(args, "inputs")?;
        if inputs.len() < 2 {
            return Err(ToolError::Validation(
                "'inputs' must list at least 2 acquisitions".to_string(),
            ));
        }
        let pairing = choice_or(args, "pairs", &PAIRS, "consecutive")?;
        let reference_index = usize_or(args, "reference_index", 0)?;
        let bias_correction = bool_or(args, "bias_correction", true)?;
        let stats = parse_stats(args)?;
        let dates = parse_dates(args, inputs.len())?;
        let (win_r, win_a) = parse_window(args)?;

        let rasters: Vec<Raster> = inputs
            .iter()
            .map(|p| load_input_raster(p))
            .collect::<Result<_, _>>()?;
        for (i, r) in rasters.iter().enumerate() {
            if r.bands < 2 {
                return Err(ToolError::Validation(format!(
                    "input {i} ('{}') must be a complex two-band I/Q raster; got {} band(s)",
                    inputs[i], r.bands
                )));
            }
        }
        let refs: Vec<&Raster> = rasters.iter().collect();
        check_alignment_refs(&refs)?;

        let n = rasters.len();
        if reference_index >= n {
            return Err(ToolError::Validation(format!(
                "'reference_index' {reference_index} is out of range for {n} input(s)"
            )));
        }

        // Build the pair list, capped BEFORE any estimation runs. N(N-1)/2
        // grows quadratically and each pair is a full windowed pass.
        let pair_list: Vec<(usize, usize)> = match pairing {
            "all" => {
                let count = n * (n - 1) / 2;
                if count > MAX_PAIRS {
                    return Err(ToolError::Validation(format!(
                        "pairs='all' over {n} acquisitions needs {count} pair(s), over the \
                         {MAX_PAIRS} cap; use pairs='consecutive' or 'reference', or split the \
                         stack"
                    )));
                }
                (0..n).flat_map(|a| ((a + 1)..n).map(move |b| (a, b))).collect()
            }
            "reference" => (0..n).filter(|&i| i != reference_index).map(|i| (reference_index, i)).collect(),
            _ => (0..n - 1).map(|i| (i, i + 1)).collect(),
        };
        ctx.progress.info(&format!(
            "{n} acquisition(s), {} pair(s), {win_r}x{win_a} window",
            pair_list.len()
        ));

        let rows = rasters[0].rows;
        let cols = rasters[0].cols;
        let cells = rows * cols;

        // Per-pair coherence, NaN where not estimable.
        let mut per_pair: Vec<Vec<f64>> = Vec::with_capacity(pair_list.len());
        let iq: Vec<(Vec<f64>, Vec<f64>)> = rasters.iter().map(complex_bands).collect();
        for (k, &(a, b)) in pair_list.iter().enumerate() {
            let est = estimate_coherence(
                &iq[a].0,
                &iq[a].1,
                &iq[b].0,
                &iq[b].1,
                rows,
                cols,
                win_r,
                win_a,
                bias_correction,
            );
            per_pair.push(est.coherence);
            ctx.progress
                .progress((k as f64 + 1.0) / pair_list.len() as f64);
        }

        // Temporal baseline per pair, for the slope fit.
        let baselines: Option<Vec<f64>> = dates
            .as_ref()
            .map(|d| pair_list.iter().map(|&(a, b)| (d[b] - d[a]).abs()).collect());

        let nodata = -1.0_f64;
        let mut bands: Vec<Vec<f64>> = Vec::with_capacity(stats.len());
        for stat in &stats {
            let mut buf = vec![nodata; cells];
            for cell in 0..cells {
                let vals: Vec<f64> = per_pair
                    .iter()
                    .map(|p| p[cell])
                    .filter(|v| v.is_finite())
                    .collect();
                if vals.is_empty() {
                    continue;
                }
                buf[cell] = match stat.as_str() {
                    "min" => vals.iter().cloned().fold(f64::INFINITY, f64::min),
                    "max" => vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
                    "std" => {
                        if vals.len() < 2 {
                            0.0
                        } else {
                            let m = vals.iter().sum::<f64>() / vals.len() as f64;
                            (vals.iter().map(|v| (v - m).powi(2)).sum::<f64>()
                                / (vals.len() - 1) as f64)
                                .sqrt()
                        }
                    }
                    "decorrelation_slope" => {
                        // Least-squares slope of coherence against temporal
                        // baseline. Needs the pairs that survived at THIS cell,
                        // so the baselines are re-gathered with the same filter.
                        let Some(bl) = &baselines else { continue };
                        let pts: Vec<(f64, f64)> = per_pair
                            .iter()
                            .zip(bl.iter())
                            .filter(|(p, _)| p[cell].is_finite())
                            .map(|(p, t)| (*t, p[cell]))
                            .collect();
                        match slope(&pts) {
                            Some(s) => s,
                            None => continue,
                        }
                    }
                    _ => vals.iter().sum::<f64>() / vals.len() as f64,
                };
            }
            bands.push(buf);
        }

        let template = &rasters[0];
        let out = raster_like_multiband(template, &bands, nodata, DataType::F32)?;
        let out_path = write_or_store_output(out, parse_optional_output(args, "output")?)?;

        // The per-pair cube is emitted unconditionally so a caller with no
        // scratch path still gets it back (the round-16 lesson).
        let pair_bands: Vec<Vec<f64>> = per_pair
            .iter()
            .map(|p| {
                p.iter()
                    .map(|v| if v.is_finite() { *v } else { nodata })
                    .collect()
            })
            .collect();
        let pair_raster = raster_like_multiband(template, &pair_bands, nodata, DataType::F32)?;
        let pair_path =
            write_or_store_output(pair_raster, parse_optional_output(args, "output_pairs")?)?;

        let pair_json: Vec<Value> = pair_list
            .iter()
            .enumerate()
            .map(|(k, &(a, b))| {
                json!({
                    "band": k + 1,
                    "reference": a,
                    "secondary": b,
                    "temporal_baseline": baselines.as_ref().map(|bl| bl[k]),
                })
            })
            .collect();

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("output_pairs".to_string(), json!(pair_path));
        outputs.insert("acquisition_count".to_string(), json!(n));
        outputs.insert("pair_count".to_string(), json!(pair_list.len()));
        outputs.insert("statistics".to_string(), json!(stats));
        outputs.insert("pairs".to_string(), json!(pair_json));
        outputs.insert("window_range".to_string(), json!(win_r));
        outputs.insert("window_azimuth".to_string(), json!(win_a));
        Ok(ToolRunResult { outputs })
    }
}

/// Ordinary least-squares slope of y against x. `None` when the x values do not
/// vary, where the slope is undefined rather than zero.
fn slope(pts: &[(f64, f64)]) -> Option<f64> {
    if pts.len() < 2 {
        return None;
    }
    let n = pts.len() as f64;
    let mx = pts.iter().map(|p| p.0).sum::<f64>() / n;
    let my = pts.iter().map(|p| p.1).sum::<f64>() / n;
    let mut num = 0.0;
    let mut den = 0.0;
    for (x, y) in pts {
        num += (x - mx) * (y - my);
        den += (x - mx) * (x - mx);
    }
    (den > 0.0).then(|| num / den)
}

fn parse_list(args: &ToolArgs, key: &str) -> Result<Vec<String>, ToolError> {
    match args.get(key) {
        Some(Value::String(s)) => Ok(s
            .split([',', ';'])
            .map(str::trim)
            .filter(|p| !p.is_empty())
            .map(str::to_string)
            .collect()),
        Some(Value::Array(items)) => items
            .iter()
            .map(|v| {
                v.as_str().map(str::to_string).ok_or_else(|| {
                    ToolError::Validation(format!("every entry of '{key}' must be a string"))
                })
            })
            .collect(),
        Some(_) => Err(ToolError::Validation(format!(
            "'{key}' must be a delimited string or an array of strings"
        ))),
        None => Err(ToolError::Validation(format!(
            "missing required parameter '{key}'"
        ))),
    }
}

fn parse_stats(args: &ToolArgs) -> Result<Vec<String>, ToolError> {
    let raw = match args.get("statistics") {
        None | Some(Value::Null) => return Ok(vec!["mean".to_string()]),
        Some(Value::String(s)) if s.trim().is_empty() => return Ok(vec!["mean".to_string()]),
        Some(Value::String(s)) => s
            .split([',', ';'])
            .map(str::trim)
            .filter(|p| !p.is_empty())
            .map(|p| p.to_ascii_lowercase())
            .collect::<Vec<_>>(),
        Some(Value::Array(items)) => items
            .iter()
            .filter_map(|v| v.as_str())
            .map(|s| s.trim().to_ascii_lowercase())
            .collect(),
        Some(_) => {
            return Err(ToolError::Validation(
                "'statistics' must be a delimited string or an array of strings".to_string(),
            ))
        }
    };
    let mut out = Vec::new();
    for s in raw {
        if !STATS.contains(&s.as_str()) {
            return Err(ToolError::Validation(format!(
                "'statistics' entry '{s}' must be one of {}",
                STATS.join("|")
            )));
        }
        if !out.contains(&s) {
            out.push(s);
        }
    }
    if out.is_empty() {
        out.push("mean".to_string());
    }
    Ok(out)
}

fn parse_dates(args: &ToolArgs, n: usize) -> Result<Option<Vec<f64>>, ToolError> {
    let Some(s) = args.get("dates").and_then(Value::as_str) else {
        return Ok(None);
    };
    if s.trim().is_empty() {
        return Ok(None);
    }
    let vals: Vec<f64> = s
        .split([',', ';'])
        .map(str::trim)
        .filter(|t| !t.is_empty())
        .map(|t| {
            t.parse::<f64>()
                .map_err(|_| ToolError::Validation(format!("'dates' entry '{t}' is not a number")))
        })
        .collect::<Result<_, _>>()?;
    if vals.len() != n {
        return Err(ToolError::Validation(format!(
            "'dates' has {} value(s) but there are {n} input(s)",
            vals.len()
        )));
    }
    if vals.iter().any(|v| !v.is_finite()) {
        return Err(ToolError::Validation(
            "'dates' contains a non-finite value".to_string(),
        ));
    }
    Ok(Some(vals))
}

/// Same convention as `sar_coherence`: a single value or `range,azimuth`, and
/// the realised window is always odd because the estimator uses an inclusive
/// `-half..=half` range.
fn parse_window(args: &ToolArgs) -> Result<(usize, usize), ToolError> {
    let Some(s) = args.get("window_size").and_then(Value::as_str) else {
        let w = usize_or(args, "window_size", 5)?.max(1);
        return Ok((w | 1, w | 1));
    };
    let parts: Vec<&str> = s.split([',', ';']).map(str::trim).filter(|p| !p.is_empty()).collect();
    let nums: Vec<usize> = parts
        .iter()
        .map(|p| {
            p.parse::<usize>().map_err(|_| {
                ToolError::Validation(format!("'window_size' entry '{p}' is not an integer"))
            })
        })
        .collect::<Result<_, _>>()?;
    match nums.len() {
        1 => Ok((nums[0].max(1) | 1, nums[0].max(1) | 1)),
        2 => Ok((nums[0].max(1) | 1, nums[1].max(1) | 1)),
        _ => Err(ToolError::Validation(
            "'window_size' must be one value or 'range,azimuth'".to_string(),
        )),
    }
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

    /// A 4x4 complex scene whose phase at every cell is `phase`, amplitude 1.
    fn scene(phase: f64) -> String {
        complex(4, 4, &(0..16).map(|_| (phase.cos(), phase.sin())).collect::<Vec<_>>())
    }

    /// A 4x4 scene whose per-cell phase varies pseudo-randomly with `seed`,
    /// which decorrelates it from anything else.
    fn noisy(seed: u64) -> String {
        let iq: Vec<(f64, f64)> = (0..16)
            .map(|i| {
                let mut z = seed
                    .wrapping_mul(0x9E37_79B9_7F4A_7C15)
                    .wrapping_add((i as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9));
                z ^= z >> 31;
                let p = (z % 6283) as f64 / 1000.0;
                (p.cos(), p.sin())
            })
            .collect();
        complex(4, 4, &iq)
    }

    fn complex(cols: usize, rows: usize, iq: &[(f64, f64)]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 2,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F64,
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

    fn run(args: Value) -> (Raster, Raster, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = MultitemporalCoherenceTool.run(&args, &ctx()).unwrap();
        let out = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        let pairs = load_input_raster(res.outputs["output_pairs"].as_str().unwrap()).unwrap();
        (out, pairs, res)
    }

    #[test]
    fn consecutive_pairing_gives_n_minus_one_pairs() {
        let a = scene(0.0);
        let b = scene(0.5);
        let c = scene(1.0);
        let (_, pairs, res) = run(json!({"inputs": format!("{a},{b},{c}")}));
        assert_eq!(res.outputs["pair_count"], json!(2));
        assert_eq!(pairs.bands, 2);
    }

    #[test]
    fn all_pairing_gives_every_unordered_pair() {
        let a = scene(0.0);
        let b = scene(0.5);
        let c = scene(1.0);
        let (_, _, res) = run(json!({
            "inputs": format!("{a},{b},{c}"), "pairs": "all",
        }));
        assert_eq!(res.outputs["pair_count"], json!(3));
    }

    #[test]
    fn reference_pairing_compares_everything_against_one_scene() {
        let a = scene(0.0);
        let b = scene(0.5);
        let c = scene(1.0);
        let (_, _, res) = run(json!({
            "inputs": format!("{a},{b},{c}"), "pairs": "reference", "reference_index": 1,
        }));
        assert_eq!(res.outputs["pair_count"], json!(2));
        let pairs = res.outputs["pairs"].as_array().unwrap();
        assert!(pairs.iter().all(|p| p["reference"] == json!(1)));
    }

    #[test]
    fn a_stack_of_stable_scenes_is_fully_coherent() {
        // Constant phase offsets preserve the scattering geometry, so mean
        // coherence must sit at 1.
        let a = scene(0.0);
        let b = scene(0.5);
        let c = scene(2.0);
        let (out, _, _) = run(json!({"inputs": format!("{a},{b},{c}")}));
        for r in 0..4 {
            for c2 in 0..4 {
                let v = out.get(0, r, c2);
                assert!((v - 1.0).abs() < 1e-6, "expected 1.0, got {v}");
            }
        }
    }

    #[test]
    fn a_decorrelating_stack_scores_far_below_one() {
        let a = noisy(1);
        let b = noisy(2);
        let c = noisy(3);
        let (out, _, _) = run(json!({"inputs": format!("{a},{b},{c}")}));
        let v = out.get(0, 2, 2);
        assert!(v < 0.7, "random phase should decorrelate, got {v}");
    }

    #[test]
    fn mean_over_the_stack_lies_between_min_and_max() {
        let a = scene(0.0);
        let b = noisy(7);
        let c = scene(1.0);
        let (out, _, res) = run(json!({
            "inputs": format!("{a},{b},{c}"), "pairs": "all",
            "statistics": "mean,min,max",
        }));
        assert_eq!(res.outputs["statistics"], json!(["mean", "min", "max"]));
        let (mean, min, max) = (out.get(0, 2, 2), out.get(1, 2, 2), out.get(2, 2, 2));
        assert!(min <= mean + 1e-9 && mean <= max + 1e-9, "{min} {mean} {max}");
        assert!(min < max, "a mixed stack should have spread");
    }

    #[test]
    fn one_band_is_emitted_per_requested_statistic_in_order() {
        let a = scene(0.0);
        let b = scene(0.5);
        let (out, _, _) = run(json!({
            "inputs": format!("{a},{b}"), "statistics": "max,mean",
        }));
        assert_eq!(out.bands, 2);
    }

    #[test]
    fn duplicate_statistics_are_collapsed() {
        let a = scene(0.0);
        let b = scene(0.5);
        let (out, _, res) = run(json!({
            "inputs": format!("{a},{b}"), "statistics": "mean,mean,std",
        }));
        assert_eq!(res.outputs["statistics"], json!(["mean", "std"]));
        assert_eq!(out.bands, 2);
    }

    #[test]
    fn std_of_a_uniformly_coherent_stack_is_zero() {
        let a = scene(0.0);
        let b = scene(0.5);
        let c = scene(1.0);
        let (out, _, _) = run(json!({
            "inputs": format!("{a},{b},{c}"), "statistics": "std",
        }));
        assert!(out.get(0, 2, 2).abs() < 1e-9);
    }

    #[test]
    fn decorrelation_slope_is_negative_when_coherence_decays_with_baseline() {
        // Scene 0 and 1 are identical (short baseline, coherence 1); scene 2 is
        // noise at a long baseline. Against the reference, coherence falls as
        // the baseline grows, so the fitted slope must be negative.
        let a = scene(0.0);
        let b = scene(0.0);
        let c = noisy(11);
        let (out, _, _) = run(json!({
            "inputs": format!("{a},{b},{c}"), "pairs": "reference",
            "dates": "0,1,100", "statistics": "decorrelation_slope",
        }));
        let s = out.get(0, 2, 2);
        assert!(s < 0.0, "expected decay, got slope {s}");
    }

    #[test]
    fn temporal_baselines_are_reported_per_pair() {
        let a = scene(0.0);
        let b = scene(0.5);
        let c = scene(1.0);
        let (_, _, res) = run(json!({
            "inputs": format!("{a},{b},{c}"), "dates": "0,12,36",
        }));
        let pairs = res.outputs["pairs"].as_array().unwrap();
        assert_eq!(pairs[0]["temporal_baseline"], json!(12.0));
        assert_eq!(pairs[1]["temporal_baseline"], json!(24.0));
    }

    #[test]
    fn the_per_pair_cube_matches_what_sar_coherence_computes_alone() {
        // The shared-estimator guarantee: a pair inside the stack must equal
        // the pairwise tool's answer for the same two scenes and window.
        let a = scene(0.0);
        let b = noisy(5);
        let (_, pairs, _) = run(json!({
            "inputs": format!("{a},{b}"), "window_size": 3,
        }));
        let args: ToolArgs = serde_json::from_value(json!({
            "reference": a, "secondary": b, "window_size": "3",
        }))
        .unwrap();
        let res = crate::sar_coherence::SarCoherenceTool.run(&args, &ctx()).unwrap();
        let single = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        for r in 0..4 {
            for c in 0..4 {
                let (x, y) = (pairs.get(0, r, c), single.get(0, r, c));
                assert!((x - y).abs() < 1e-9, "cell ({r},{c}): {x} vs {y}");
            }
        }
    }

    #[test]
    fn the_pair_cube_is_produced_without_a_path() {
        let a = scene(0.0);
        let b = scene(0.5);
        let args: ToolArgs =
            serde_json::from_value(json!({"inputs": format!("{a},{b}")})).unwrap();
        let res = MultitemporalCoherenceTool.run(&args, &ctx()).unwrap();
        let p = res.outputs["output_pairs"].as_str().unwrap();
        assert!(load_input_raster(p).is_ok());
    }

    #[test]
    fn a_single_band_input_is_refused() {
        let mut r = Raster::new(RasterConfig {
            cols: 4,
            rows: 4,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F64,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        r.set(0, 0, 0, 1.0).unwrap();
        let id = wbraster::memory_store::put_raster(r);
        let bad = wbraster::memory_store::make_raster_memory_path(&id);
        let good = scene(0.0);
        let args: ToolArgs =
            serde_json::from_value(json!({"inputs": format!("{good},{bad}")})).unwrap();
        let err = MultitemporalCoherenceTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("two-band"), "{err}");
    }

    #[test]
    fn an_oversized_all_pairing_is_capped_before_estimating() {
        // 40 acquisitions is 780 pairs, over the cap. The guard has to fire
        // before any windowed pass runs.
        let a = scene(0.0);
        let inputs: Vec<String> = (0..40).map(|_| a.clone()).collect();
        let args: ToolArgs = serde_json::from_value(json!({
            "inputs": inputs.join(","), "pairs": "all",
        }))
        .unwrap();
        let err = MultitemporalCoherenceTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("cap"), "{err}");
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            MultitemporalCoherenceTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"inputs": "a.tif"})));
        assert!(bad(json!({"inputs": "a.tif,b.tif", "pairs": "sequential"})));
        assert!(bad(json!({"inputs": "a.tif,b.tif", "statistics": "median"})));
        // A slope against acquisition index is not a slope against time.
        assert!(bad(json!({
            "inputs": "a.tif,b.tif", "statistics": "decorrelation_slope",
        })));
        assert!(bad(json!({
            "inputs": "a.tif,b.tif", "dates": "0,1,2",
        })));
        assert!(bad(json!({"inputs": "a.tif,b.tif", "reference_index": 5})));
    }
}
