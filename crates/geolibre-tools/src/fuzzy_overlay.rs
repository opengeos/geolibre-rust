//! GeoLibre tool: fuzzy membership transforms and fuzzy overlay combination.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Fuzzy Membership* and *Fuzzy Overlay*
//! (Spatial Analyst). The bundled `weighted_overlay`/`weighted_sum` are crisp
//! reclass-and-add, and the only fuzzy tool is `fuzzy_knn_classification` (a
//! classifier); the standard multi-criteria suitability workflow — rescale each
//! criterion raster to a 0..1 membership surface, then combine the surfaces —
//! has no home. This single tool covers both halves, all closed-form cell math
//! with no dependencies.
//!
//! Two modes, chosen by which input parameter is present:
//!
//! * **Membership** (`input`): transform one band into a 0..1 membership surface
//!   with `function` = `linear` | `gaussian` | `small` | `large` | `ms_small` |
//!   `ms_large`. Function parameters (`midpoint`, `spread`, `min`, `max`)
//!   default to values derived from the band's statistics.
//! * **Overlay** (`inputs`, a comma-separated list of ≥2 membership rasters):
//!   combine them cell-wise with `overlay` = `and` (min) | `or` (max) |
//!   `product` | `sum` (algebraic sum 1−∏(1−xᵢ)) | `gamma`
//!   (sumᵞ · product¹⁻ᵞ).
//!
//! No-data propagates: a cell that is no-data in the input (or in any overlay
//! input) is no-data in the output. Output is a single F32 band clamped to
//! [0, 1] with a distinct no-data value.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::common::{load_input_raster, parse_optional_output, raster_like_with_data};

/// No-data value for the 0..1 membership output (outside the valid range).
const OUT_NODATA: f64 = -9999.0;

#[derive(Clone, Copy, PartialEq)]
enum MembershipFn {
    Linear,
    Gaussian,
    Small,
    Large,
    MsSmall,
    MsLarge,
}

#[derive(Clone, Copy, PartialEq)]
enum Overlay {
    And,
    Or,
    Product,
    Sum,
    Gamma,
}

pub struct FuzzyOverlayTool;

impl Tool for FuzzyOverlayTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "fuzzy_overlay",
            display_name: "Fuzzy Overlay",
            summary: "Rescale a raster to a 0..1 fuzzy membership surface (linear/gaussian/small/large/ms_small/ms_large), or combine several membership rasters with fuzzy AND/OR/PRODUCT/SUM/GAMMA, like ArcGIS Fuzzy Membership + Fuzzy Overlay.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Membership mode: input raster to transform into a 0..1 surface.",
                    required: false,
                },
                ToolParamSpec {
                    name: "inputs",
                    description: "Overlay mode: comma-separated list of ≥2 membership rasters to combine. Takes precedence over 'input'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output raster path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "function",
                    description: "Membership function: 'linear' (default), 'gaussian', 'small', 'large', 'ms_small', 'ms_large'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "overlay",
                    description: "Overlay operator: 'and' (min, default), 'or' (max), 'product', 'sum', 'gamma'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "midpoint",
                    description: "linear: value mapping to 0.5 is ignored (use min/max); gaussian/small/large: the function midpoint; ms_small/ms_large: the mean multiplier (default from band stats).",
                    required: false,
                },
                ToolParamSpec {
                    name: "spread",
                    description: "gaussian/small/large: spread (default gaussian 0.1, small/large 5); ms_small/ms_large: the standard-deviation multiplier (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min",
                    description: "linear: value mapping to membership 0 (default band minimum). If min>max the ramp is decreasing.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max",
                    description: "linear: value mapping to membership 1 (default band maximum).",
                    required: false,
                },
                ToolParamSpec {
                    name: "gamma",
                    description: "gamma overlay exponent in [0,1] (default 0.9).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to transform in membership mode (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        let has_input = args
            .get("input")
            .and_then(Value::as_str)
            .map(|s| !s.trim().is_empty())
            .unwrap_or(false);
        let has_inputs = args
            .get("inputs")
            .and_then(Value::as_str)
            .map(|s| !s.trim().is_empty())
            .unwrap_or(false);
        if !has_input && !has_inputs {
            return Err(ToolError::Validation(
                "provide 'input' (membership mode) or 'inputs' (overlay mode)".to_string(),
            ));
        }
        parse_membership_fn(args)?;
        parse_overlay(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let output = parse_optional_output(args, "output")?;

        // Overlay mode takes precedence when 'inputs' is present.
        if let Some(list) = args
            .get("inputs")
            .and_then(Value::as_str)
            .filter(|s| !s.trim().is_empty())
        {
            return self.run_overlay(list, args, output, ctx);
        }

        let input = args
            .get("input")
            .and_then(Value::as_str)
            .filter(|s| !s.trim().is_empty())
            .ok_or_else(|| {
                ToolError::Validation("missing required parameter 'input'".to_string())
            })?;
        self.run_membership(input, args, output, ctx)
    }
}

impl FuzzyOverlayTool {
    fn run_membership(
        &self,
        input: &str,
        args: &ToolArgs,
        output: Option<&str>,
        ctx: &ToolContext,
    ) -> Result<ToolRunResult, ToolError> {
        let func = parse_membership_fn(args)?;
        let band_1based = args.get("band").and_then(Value::as_u64).unwrap_or(1).max(1);
        let band = (band_1based - 1) as isize;

        let raster = load_input_raster(input)?;
        if band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "band {band_1based} out of range (raster has {} band(s))",
                raster.bands
            )));
        }
        let nodata = raster.nodata;
        let rows = raster.rows;
        let cols = raster.cols;

        // Band statistics (min/max/mean/std) drive the parameter defaults.
        let mut min = f64::INFINITY;
        let mut max = f64::NEG_INFINITY;
        let mut sum = 0.0;
        let mut sumsq = 0.0;
        let mut n = 0u64;
        for row in 0..rows as isize {
            for col in 0..cols as isize {
                let v = raster.get(band, row, col);
                if v != nodata && v.is_finite() {
                    min = min.min(v);
                    max = max.max(v);
                    sum += v;
                    sumsq += v * v;
                    n += 1;
                }
            }
        }
        if n == 0 {
            return Err(ToolError::Execution(
                "raster band contains no valid (non-nodata) values".to_string(),
            ));
        }
        let mean = sum / n as f64;
        let var = (sumsq / n as f64 - mean * mean).max(0.0);
        let std = var.sqrt();

        let p = FnParams::resolve(func, args, min, max, mean, std)?;
        ctx.progress.info("applying fuzzy membership");

        let mut data = vec![OUT_NODATA; rows * cols];
        for row in 0..rows {
            for col in 0..cols {
                let v = raster.get(band, row as isize, col as isize);
                if v != nodata && v.is_finite() {
                    data[row * cols + col] = membership(func, v, &p).clamp(0.0, 1.0);
                }
            }
        }

        let out = raster_like_with_data(&raster, data, OUT_NODATA, DataType::F32)?;
        let out_path = crate::common::write_or_store_output(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("mode".to_string(), json!("membership"));
        outputs.insert("function".to_string(), json!(func_name(func)));
        outputs.insert("band_min".to_string(), json!(min));
        outputs.insert("band_max".to_string(), json!(max));
        outputs.insert("band_mean".to_string(), json!(mean));
        outputs.insert("band_std".to_string(), json!(std));
        Ok(ToolRunResult { outputs })
    }

    fn run_overlay(
        &self,
        list: &str,
        args: &ToolArgs,
        output: Option<&str>,
        ctx: &ToolContext,
    ) -> Result<ToolRunResult, ToolError> {
        let op = parse_overlay(args)?;
        let gamma = parse_f64(args, "gamma")?.unwrap_or(0.9);
        if op == Overlay::Gamma && !(0.0..=1.0).contains(&gamma) {
            return Err(ToolError::Validation(
                "'gamma' must be in [0, 1]".to_string(),
            ));
        }

        let paths: Vec<&str> = list
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .collect();
        if paths.len() < 2 {
            return Err(ToolError::Validation(
                "overlay mode needs at least 2 rasters in 'inputs'".to_string(),
            ));
        }

        let rasters: Vec<Raster> = paths
            .iter()
            .map(|p| load_input_raster(p))
            .collect::<Result<_, _>>()?;
        let rows = rasters[0].rows;
        let cols = rasters[0].cols;
        for (i, r) in rasters.iter().enumerate() {
            if r.rows != rows || r.cols != cols {
                return Err(ToolError::Validation(format!(
                    "input {} is {}x{}, expected {rows}x{cols} — all overlay rasters must align",
                    i, r.rows, r.cols
                )));
            }
        }

        ctx.progress
            .info(&format!("fuzzy overlay of {} rasters", rasters.len()));
        let mut data = vec![OUT_NODATA; rows * cols];
        for row in 0..rows {
            for col in 0..cols {
                let mut vals: Vec<f64> = Vec::with_capacity(rasters.len());
                let mut ok = true;
                for r in &rasters {
                    let v = r.get(0, row as isize, col as isize);
                    if v == r.nodata || !v.is_finite() {
                        ok = false;
                        break;
                    }
                    vals.push(v.clamp(0.0, 1.0));
                }
                if ok {
                    data[row * cols + col] = combine(op, &vals, gamma).clamp(0.0, 1.0);
                }
            }
        }

        let out = raster_like_with_data(&rasters[0], data, OUT_NODATA, DataType::F32)?;
        let out_path = crate::common::write_or_store_output(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("mode".to_string(), json!("overlay"));
        outputs.insert("overlay".to_string(), json!(overlay_name(op)));
        outputs.insert("input_count".to_string(), json!(rasters.len()));
        Ok(ToolRunResult { outputs })
    }
}

// ── Membership functions ──────────────────────────────────────────────────────

struct FnParams {
    midpoint: f64,
    spread: f64,
    lo: f64,
    hi: f64,
    mean: f64,
    std: f64,
}

impl FnParams {
    fn resolve(
        func: MembershipFn,
        args: &ToolArgs,
        min: f64,
        max: f64,
        mean: f64,
        std: f64,
    ) -> Result<FnParams, ToolError> {
        let midpoint = parse_f64(args, "midpoint")?;
        let spread = parse_f64(args, "spread")?;
        let lo = parse_f64(args, "min")?.unwrap_or(min);
        let hi = parse_f64(args, "max")?.unwrap_or(max);
        let default_spread = match func {
            MembershipFn::Gaussian => 0.1,
            MembershipFn::Small | MembershipFn::Large => 5.0,
            MembershipFn::MsSmall | MembershipFn::MsLarge => 1.0, // std multiplier
            MembershipFn::Linear => 1.0,
        };
        let default_mid = match func {
            MembershipFn::MsSmall | MembershipFn::MsLarge => 1.0, // mean multiplier
            _ => mean,
        };
        Ok(FnParams {
            midpoint: midpoint.unwrap_or(default_mid),
            spread: spread.unwrap_or(default_spread),
            lo,
            hi,
            mean,
            std,
        })
    }
}

fn membership(func: MembershipFn, x: f64, p: &FnParams) -> f64 {
    match func {
        MembershipFn::Linear => {
            let (lo, hi) = (p.lo, p.hi);
            if hi == lo {
                return if x >= hi { 1.0 } else { 0.0 };
            }
            ((x - lo) / (hi - lo)).clamp(0.0, 1.0)
        }
        // Near/Gaussian: peak 1 at the midpoint, falling off with spread.
        MembershipFn::Gaussian => (-p.spread * (x - p.midpoint).powi(2)).exp(),
        // Small: high membership for small x. 1 / (1 + (x/mid)^spread).
        MembershipFn::Small => {
            if p.midpoint == 0.0 {
                return 0.0;
            }
            1.0 / (1.0 + (x / p.midpoint).powf(p.spread))
        }
        // Large: high membership for large x. 1 / (1 + (x/mid)^-spread).
        MembershipFn::Large => {
            if p.midpoint == 0.0 || x <= 0.0 {
                return 0.0;
            }
            1.0 / (1.0 + (x / p.midpoint).powf(-p.spread))
        }
        // MSSmall: mean/std based; a = mean multiplier, b = std multiplier.
        MembershipFn::MsSmall => {
            let a = p.midpoint;
            let b = p.spread;
            let t = x - a * p.mean;
            let bs = b * p.std;
            if t <= 0.0 {
                1.0
            } else {
                bs / (t + bs)
            }
        }
        // MSLarge: high membership for large x, mean/std based.
        MembershipFn::MsLarge => {
            let a = p.midpoint;
            let b = p.spread;
            let t = x - a * p.mean;
            let bs = b * p.std;
            if t <= 0.0 {
                0.0
            } else {
                t / (t + bs)
            }
        }
    }
}

// ── Overlay operators ─────────────────────────────────────────────────────────

fn combine(op: Overlay, vals: &[f64], gamma: f64) -> f64 {
    match op {
        Overlay::And => vals.iter().cloned().fold(f64::INFINITY, f64::min),
        Overlay::Or => vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        Overlay::Product => vals.iter().product(),
        Overlay::Sum => 1.0 - vals.iter().map(|v| 1.0 - v).product::<f64>(),
        Overlay::Gamma => {
            let product: f64 = vals.iter().product();
            let sum = 1.0 - vals.iter().map(|v| 1.0 - v).product::<f64>();
            sum.powf(gamma) * product.powf(1.0 - gamma)
        }
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

fn parse_membership_fn(args: &ToolArgs) -> Result<MembershipFn, ToolError> {
    match args
        .get("function")
        .and_then(Value::as_str)
        .map(|s| s.trim().to_lowercase())
    {
        None => Ok(MembershipFn::Linear),
        Some(s) if s.is_empty() || s == "linear" => Ok(MembershipFn::Linear),
        Some(s) if s == "gaussian" || s == "near" => Ok(MembershipFn::Gaussian),
        Some(s) if s == "small" => Ok(MembershipFn::Small),
        Some(s) if s == "large" => Ok(MembershipFn::Large),
        Some(s) if s == "ms_small" || s == "mssmall" => Ok(MembershipFn::MsSmall),
        Some(s) if s == "ms_large" || s == "mslarge" => Ok(MembershipFn::MsLarge),
        Some(other) => Err(ToolError::Validation(format!(
            "'function' must be linear|gaussian|small|large|ms_small|ms_large, got '{other}'"
        ))),
    }
}

fn parse_overlay(args: &ToolArgs) -> Result<Overlay, ToolError> {
    match args
        .get("overlay")
        .and_then(Value::as_str)
        .map(|s| s.trim().to_lowercase())
    {
        None => Ok(Overlay::And),
        Some(s) if s.is_empty() || s == "and" => Ok(Overlay::And),
        Some(s) if s == "or" => Ok(Overlay::Or),
        Some(s) if s == "product" => Ok(Overlay::Product),
        Some(s) if s == "sum" => Ok(Overlay::Sum),
        Some(s) if s == "gamma" => Ok(Overlay::Gamma),
        Some(other) => Err(ToolError::Validation(format!(
            "'overlay' must be and|or|product|sum|gamma, got '{other}'"
        ))),
    }
}

fn parse_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("'{key}' must be a number"))),
        Some(_) => Err(ToolError::Validation(format!("'{key}' must be a number"))),
    }
}

fn func_name(f: MembershipFn) -> &'static str {
    match f {
        MembershipFn::Linear => "linear",
        MembershipFn::Gaussian => "gaussian",
        MembershipFn::Small => "small",
        MembershipFn::Large => "large",
        MembershipFn::MsSmall => "ms_small",
        MembershipFn::MsLarge => "ms_large",
    }
}

fn overlay_name(o: Overlay) -> &'static str {
    match o {
        Overlay::And => "and",
        Overlay::Or => "or",
        Overlay::Product => "product",
        Overlay::Sum => "sum",
        Overlay::Gamma => "gamma",
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

    /// Build a 1-band raster from a row-major values buffer.
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
        FuzzyOverlayTool.run(&args, &ctx()).unwrap()
    }

    /// Linear membership maps band min->0, max->1, midpoint->0.5.
    #[test]
    fn linear_membership_ramps_zero_to_one() {
        let input = raster_of(1, 3, &[0.0, 5.0, 10.0], -9999.0);
        let out = run(json!({ "input": input, "function": "linear" }));
        let (v, _) = read_all(out.outputs["output"].as_str().unwrap());
        assert!((v[0] - 0.0).abs() < 1e-6);
        assert!((v[1] - 0.5).abs() < 1e-6);
        assert!((v[2] - 1.0).abs() < 1e-6);
    }

    /// Gaussian peaks at the midpoint and falls off symmetrically.
    #[test]
    fn gaussian_peaks_at_midpoint() {
        let input = raster_of(1, 3, &[0.0, 5.0, 10.0], -9999.0);
        let out =
            run(json!({ "input": input, "function": "gaussian", "midpoint": 5.0, "spread": 0.1 }));
        let (v, _) = read_all(out.outputs["output"].as_str().unwrap());
        assert!((v[1] - 1.0).abs() < 1e-9, "peak at midpoint");
        assert!(v[0] < v[1] && v[2] < v[1]);
        assert!((v[0] - v[2]).abs() < 1e-9, "symmetric");
    }

    /// 'large' is monotonically increasing; 'small' monotonically decreasing.
    #[test]
    fn small_and_large_monotonic() {
        let input = raster_of(1, 4, &[1.0, 2.0, 4.0, 8.0], -9999.0);
        let large =
            run(json!({ "input": input, "function": "large", "midpoint": 4.0, "spread": 5.0 }));
        let (lv, _) = read_all(large.outputs["output"].as_str().unwrap());
        assert!(
            lv[0] < lv[1] && lv[1] < lv[2] && lv[2] < lv[3],
            "large increasing"
        );
        let small =
            run(json!({ "input": input, "function": "small", "midpoint": 4.0, "spread": 5.0 }));
        let (sv, _) = read_all(small.outputs["output"].as_str().unwrap());
        assert!(
            sv[0] > sv[1] && sv[1] > sv[2] && sv[2] > sv[3],
            "small decreasing"
        );
    }

    /// Fuzzy AND is the cell-wise minimum; OR the maximum.
    #[test]
    fn overlay_and_or() {
        let a = raster_of(1, 2, &[0.2, 0.9], -9999.0);
        let b = raster_of(1, 2, &[0.6, 0.4], -9999.0);
        let and = run(json!({ "inputs": format!("{a},{b}"), "overlay": "and" }));
        let (v, _) = read_all(and.outputs["output"].as_str().unwrap());
        assert!((v[0] - 0.2).abs() < 1e-6 && (v[1] - 0.4).abs() < 1e-6);
        let or = run(json!({ "inputs": format!("{a},{b}"), "overlay": "or" }));
        let (v, _) = read_all(or.outputs["output"].as_str().unwrap());
        assert!((v[0] - 0.6).abs() < 1e-6 && (v[1] - 0.9).abs() < 1e-6);
    }

    /// Product and algebraic sum match their formulas.
    #[test]
    fn overlay_product_and_sum() {
        let a = raster_of(1, 1, &[0.5], -9999.0);
        let b = raster_of(1, 1, &[0.5], -9999.0);
        let prod = run(json!({ "inputs": format!("{a},{b}"), "overlay": "product" }));
        let (v, _) = read_all(prod.outputs["output"].as_str().unwrap());
        assert!((v[0] - 0.25).abs() < 1e-6, "0.5*0.5=0.25");
        let sum = run(json!({ "inputs": format!("{a},{b}"), "overlay": "sum" }));
        let (v, _) = read_all(sum.outputs["output"].as_str().unwrap());
        assert!((v[0] - 0.75).abs() < 1e-6, "1-(0.5)(0.5)=0.75");
    }

    /// No-data in any overlay input propagates to no-data out.
    #[test]
    fn overlay_propagates_nodata() {
        let a = raster_of(1, 2, &[0.5, -9999.0], -9999.0);
        let b = raster_of(1, 2, &[0.5, 0.5], -9999.0);
        let out = run(json!({ "inputs": format!("{a},{b}"), "overlay": "and" }));
        let (v, nd) = read_all(out.outputs["output"].as_str().unwrap());
        assert!((v[0] - 0.5).abs() < 1e-6);
        assert_eq!(v[1], nd, "nodata propagated");
    }

    #[test]
    fn rejects_no_input() {
        let args: ToolArgs = serde_json::from_value(json!({ "function": "linear" })).unwrap();
        assert!(FuzzyOverlayTool.validate(&args).is_err());
    }

    #[test]
    fn rejects_bad_function() {
        let input = raster_of(1, 1, &[1.0], -9999.0);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "function": "bogus" })).unwrap();
        assert!(FuzzyOverlayTool.validate(&args).is_err());
    }
}
