//! GeoLibre tool: rescale a continuous raster onto a common suitability scale
//! through a choice of transfer functions.
//!
//! Pure-Rust counterpart of ArcGIS Spatial Analyst's *Rescale By Function*. Each
//! transfer function maps the input values (normalized against the band range)
//! monotonically or unimodally onto `[from_scale, to_scale]`, so disparate
//! criteria layers can be combined in a suitability model.
//!
//! `fuzzy_overlay` already ships several fuzzy-membership curves as part of its
//! overlay pipeline; this tool exposes the rescale step on its own, with an
//! arbitrary output range and optional below/above-threshold cut-offs, so a
//! single criterion can be rescaled without running an overlay combination.
//!
//! No-data cells are preserved. The band min/max used for normalization are
//! computed once up front; a degenerate (single-value) band maps every valid
//! cell to `from_scale`.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};

use crate::common::{load_input_raster, parse_optional_output, write_or_store_output};

pub struct RescaleByFunctionTool;

impl Tool for RescaleByFunctionTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "rescale_by_function",
            display_name: "Rescale By Function",
            summary: "Rescale a continuous raster onto a common suitability scale via a transfer function (linear, inverse_linear, power, exponential, logarithmic, logistic, gaussian, near, small, large, symmetric_linear), like ArcGIS Rescale By Function.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input continuous raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output raster path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "function",
                    description: "Transfer function: linear (default), inverse_linear, power, exponential, logarithmic, logistic, gaussian, near, small, large, symmetric_linear.",
                    required: false,
                },
                ToolParamSpec {
                    name: "from_scale",
                    description: "Low end of the output scale (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "to_scale",
                    description: "High end of the output scale (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "param1",
                    description: "First shape parameter (exponent for power; steepness for exponential/logarithmic/logistic/small/large; center for gaussian/near). Function-specific default.",
                    required: false,
                },
                ToolParamSpec {
                    name: "param2",
                    description: "Second shape parameter (spread for gaussian/near). Function-specific default.",
                    required: false,
                },
                ToolParamSpec {
                    name: "low_threshold",
                    description: "Optional lower cut-off: values below map to 'value_below'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "high_threshold",
                    description: "Optional upper cut-off: values above map to 'value_above'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "value_below",
                    description: "Output value for cells below 'low_threshold' (default = from_scale).",
                    required: false,
                },
                ToolParamSpec {
                    name: "value_above",
                    description: "Output value for cells above 'high_threshold' (default = to_scale).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to rescale (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_function(args)?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_output(args, "output")?;
        let func = parse_function(args)?;
        let p = parse_params(args)?;

        let mut raster = load_input_raster(input)?;
        if p.band_0 as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "band {} out of range (raster has {} band(s))",
                p.band_0 + 1,
                raster.bands
            )));
        }
        let band = p.band_0;
        let nodata = raster.nodata;
        let rows = raster.rows as isize;
        let cols = raster.cols as isize;

        // Band range for normalization.
        let (mut vmin, mut vmax) = (f64::INFINITY, f64::NEG_INFINITY);
        for r in 0..rows {
            for c in 0..cols {
                let v = raster.get(band, r, c);
                if v != nodata && v.is_finite() {
                    vmin = vmin.min(v);
                    vmax = vmax.max(v);
                }
            }
        }
        if !vmin.is_finite() || !vmax.is_finite() {
            return Err(ToolError::Execution(
                "raster band contains no valid (non-nodata) values".to_string(),
            ));
        }
        let range = vmax - vmin;

        ctx.progress
            .info(&format!("rescaling with '{}'", func.label()));
        let value_below = p.value_below.unwrap_or(p.from_scale);
        let value_above = p.value_above.unwrap_or(p.to_scale);

        for r in 0..rows {
            for c in 0..cols {
                let v = raster.get(band, r, c);
                if v == nodata || !v.is_finite() {
                    continue;
                }
                let out = if p.low_threshold.map(|t| v < t).unwrap_or(false) {
                    value_below
                } else if p.high_threshold.map(|t| v > t).unwrap_or(false) {
                    value_above
                } else {
                    let t = if range == 0.0 {
                        0.0
                    } else {
                        (v - vmin) / range
                    };
                    let f = func.apply(t, v, vmin, vmax, &p).clamp(0.0, 1.0);
                    p.from_scale + f * (p.to_scale - p.from_scale)
                };
                raster
                    .set(band, r, c, out)
                    .map_err(|e| ToolError::Execution(format!("failed writing cell: {e}")))?;
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let out_path = write_or_store_output(raster, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("function".to_string(), json!(func.label()));
        outputs.insert("input_min".to_string(), json!(vmin));
        outputs.insert("input_max".to_string(), json!(vmax));
        Ok(ToolRunResult { outputs })
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Function {
    Linear,
    InverseLinear,
    Power,
    Exponential,
    Logarithmic,
    Logistic,
    Gaussian,
    Near,
    Small,
    Large,
    SymmetricLinear,
}

impl Function {
    fn label(self) -> &'static str {
        match self {
            Function::Linear => "linear",
            Function::InverseLinear => "inverse_linear",
            Function::Power => "power",
            Function::Exponential => "exponential",
            Function::Logarithmic => "logarithmic",
            Function::Logistic => "logistic",
            Function::Gaussian => "gaussian",
            Function::Near => "near",
            Function::Small => "small",
            Function::Large => "large",
            Function::SymmetricLinear => "symmetric_linear",
        }
    }

    /// Evaluates the transfer function, returning a value in `[0, 1]` (clamped by
    /// the caller). `t` is the range-normalized value in `[0, 1]`; `v` is the raw
    /// value with band bounds `vmin`/`vmax` for the center-based functions.
    fn apply(self, t: f64, v: f64, vmin: f64, vmax: f64, p: &Params) -> f64 {
        match self {
            Function::Linear => t,
            Function::InverseLinear => 1.0 - t,
            Function::Power => {
                let e = p.param1.unwrap_or(2.0).max(1e-6);
                t.powf(e)
            }
            Function::Exponential => {
                let k = p.param1.unwrap_or(5.0);
                if k.abs() < 1e-9 {
                    t
                } else {
                    ((k * t).exp() - 1.0) / (k.exp() - 1.0)
                }
            }
            Function::Logarithmic => {
                let k = p.param1.unwrap_or(9.0).max(1e-6);
                (1.0 + k * t).ln() / (1.0 + k).ln()
            }
            Function::Logistic | Function::Large => {
                let k = p.param1.unwrap_or(10.0);
                logistic(t, 0.5, k)
            }
            Function::Small => {
                let k = p.param1.unwrap_or(10.0);
                1.0 - logistic(t, 0.5, k)
            }
            Function::Gaussian => {
                let mid = p.param1.unwrap_or(0.5 * (vmin + vmax));
                let width = p.param2.unwrap_or((vmax - vmin) / 6.0).abs().max(1e-9);
                let d = (v - mid) / width;
                (-0.5 * d * d).exp()
            }
            Function::Near => {
                let mid = p.param1.unwrap_or(0.5 * (vmin + vmax));
                let spread = p.param2.unwrap_or((vmax - vmin) / 6.0).abs().max(1e-9);
                let d = (v - mid) / spread;
                1.0 / (1.0 + d * d)
            }
            Function::SymmetricLinear => {
                // Peaks at the midpoint, falls linearly to the ends.
                1.0 - (2.0 * t - 1.0).abs()
            }
        }
    }
}

/// Logistic sigmoid rescaled so `f(0)≈0` and `f(1)≈1` around center `c` with
/// steepness `k`.
fn logistic(t: f64, c: f64, k: f64) -> f64 {
    let raw = |x: f64| 1.0 / (1.0 + (-k * (x - c)).exp());
    let lo = raw(0.0);
    let hi = raw(1.0);
    if (hi - lo).abs() < 1e-12 {
        t
    } else {
        (raw(t) - lo) / (hi - lo)
    }
}

struct Params {
    from_scale: f64,
    to_scale: f64,
    param1: Option<f64>,
    param2: Option<f64>,
    low_threshold: Option<f64>,
    high_threshold: Option<f64>,
    value_below: Option<f64>,
    value_above: Option<f64>,
    band_0: isize,
}

fn parse_function(args: &ToolArgs) -> Result<Function, ToolError> {
    Ok(
        match args.get("function").and_then(Value::as_str).map(str::trim) {
            None | Some("") | Some("linear") => Function::Linear,
            Some("inverse_linear") => Function::InverseLinear,
            Some("power") => Function::Power,
            Some("exponential") => Function::Exponential,
            Some("logarithmic") => Function::Logarithmic,
            Some("logistic") => Function::Logistic,
            Some("gaussian") => Function::Gaussian,
            Some("near") => Function::Near,
            Some("small") => Function::Small,
            Some("large") => Function::Large,
            Some("symmetric_linear") => Function::SymmetricLinear,
            Some(o) => {
                return Err(ToolError::Validation(format!(
                    "'function' must be one of linear|inverse_linear|power|exponential|logarithmic|logistic|gaussian|near|small|large|symmetric_linear, got '{o}'"
                )))
            }
        },
    )
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let from_scale = parse_optional_f64(args, "from_scale")?.unwrap_or(0.0);
    let to_scale = parse_optional_f64(args, "to_scale")?.unwrap_or(1.0);
    let band_1 = match parse_optional_f64(args, "band")? {
        None => 1,
        Some(v) if v.fract() == 0.0 && v >= 1.0 => v as isize,
        Some(_) => {
            return Err(ToolError::Validation(
                "parameter 'band' must be a positive integer".to_string(),
            ))
        }
    };
    Ok(Params {
        from_scale,
        to_scale,
        param1: parse_optional_f64(args, "param1")?,
        param2: parse_optional_f64(args, "param2")?,
        low_threshold: parse_optional_f64(args, "low_threshold")?,
        high_threshold: parse_optional_f64(args, "high_threshold")?,
        value_below: parse_optional_f64(args, "value_below")?,
        value_above: parse_optional_f64(args, "value_above")?,
        band_0: band_1 - 1,
    })
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

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{DataType, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn raster_path(rows: usize, cols: usize, vals: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: Default::default(),
            metadata: Default::default(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, vals[row * cols + col])
                    .unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> Raster {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = RescaleByFunctionTool.run(&args, &ctx()).unwrap();
        load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    #[test]
    fn linear_maps_min_to_from_and_max_to_to() {
        let input = raster_path(1, 3, &[0.0, 5.0, 10.0]);
        let r = run(
            json!({ "input": input, "function": "linear", "from_scale": 1.0, "to_scale": 10.0 }),
        );
        assert!((r.get(0, 0, 0) - 1.0).abs() < 1e-9);
        assert!((r.get(0, 0, 1) - 5.5).abs() < 1e-9);
        assert!((r.get(0, 0, 2) - 10.0).abs() < 1e-9);
    }

    #[test]
    fn inverse_linear_flips() {
        let input = raster_path(1, 3, &[0.0, 5.0, 10.0]);
        let r = run(json!({ "input": input, "function": "inverse_linear" }));
        assert!((r.get(0, 0, 0) - 1.0).abs() < 1e-9);
        assert!((r.get(0, 0, 2) - 0.0).abs() < 1e-9);
    }

    #[test]
    fn thresholds_clamp_to_below_and_above() {
        let input = raster_path(1, 5, &[0.0, 2.0, 5.0, 8.0, 10.0]);
        let r = run(json!({
            "input": input, "function": "linear",
            "low_threshold": 2.5, "high_threshold": 7.5,
            "value_below": -1.0, "value_above": 2.0
        }));
        assert_eq!(r.get(0, 0, 0), -1.0);
        assert_eq!(r.get(0, 0, 1), -1.0);
        assert_eq!(r.get(0, 0, 4), 2.0);
        assert_eq!(r.get(0, 0, 3), 2.0);
        // Middle value untouched by clamps, in [0,1].
        assert!((0.0..=1.0).contains(&r.get(0, 0, 2)));
    }

    #[test]
    fn gaussian_peaks_at_center() {
        let input = raster_path(1, 3, &[0.0, 5.0, 10.0]);
        let r = run(json!({ "input": input, "function": "gaussian" }));
        // Midpoint (5) should have the highest membership.
        assert!(r.get(0, 0, 1) > r.get(0, 0, 0));
        assert!(r.get(0, 0, 1) > r.get(0, 0, 2));
    }

    #[test]
    fn nodata_is_preserved() {
        let input = raster_path(1, 3, &[0.0, -9999.0, 10.0]);
        let r = run(json!({ "input": input, "function": "linear" }));
        assert_eq!(r.get(0, 0, 1), -9999.0);
    }

    #[test]
    fn monotone_functions_stay_ordered() {
        let input = raster_path(1, 5, &[0.0, 1.0, 2.0, 3.0, 4.0]);
        for f in ["power", "exponential", "logarithmic", "logistic", "large"] {
            let r = run(json!({ "input": input.clone(), "function": f }));
            let vals: Vec<f64> = (0..5).map(|c| r.get(0, 0, c)).collect();
            for w in vals.windows(2) {
                assert!(w[1] >= w[0] - 1e-9, "function {f} not monotone: {vals:?}");
            }
        }
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            RescaleByFunctionTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "x.tif", "function": "bogus" })).is_err());
        assert!(bad(json!({ "input": "x.tif", "band": 0 })).is_err());
        assert!(bad(json!({ "input": "x.tif", "function": "gaussian" })).is_ok());
    }
}
