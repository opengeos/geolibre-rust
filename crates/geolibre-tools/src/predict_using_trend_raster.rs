//! GeoLibre tool: evaluate a fitted per-pixel temporal trend at arbitrary times.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Predict Using Trend Raster* (Image
//! Analyst). It is the consumer half of the shipped `generate_trend_raster`,
//! which fits a per-pixel model across a raster time series and writes the
//! fitted **slope** and **intercept** rasters — but nothing in the registry
//! evaluates those coefficients, so the fit currently cannot be turned back into
//! a predicted surface. The bundled `trend_surface` / `trend_surface_vector_points`
//! fit a *spatial* polynomial across x/y and never read a coefficient stack.
//!
//! Given the slope raster `m` and intercept raster `b` from
//! `generate_trend_raster`, the predicted value at time `t` is the same linear
//! model that tool fitted:
//!
//! ```text
//!   v(t) = b + m · t
//! ```
//!
//! Times are supplied either explicitly (`times = 12,18,24`) or as a range
//! (`start` / `end` / `interval`), and must use the same units the fit used — if
//! `generate_trend_raster` was given `times`, pass values on that scale; if it
//! defaulted to raster indices `0..n-1`, pass indices. A cell is predicted only
//! where **both** coefficients are valid, so no-data propagates from either
//! input.
//!
//! Output is one raster per requested time. With a single time the result is
//! written to `output`; with several, `output` receives the first and the rest
//! are written alongside it with a `_t<value>` suffix (or all are stored in
//! memory and returned in the `outputs` list when no path is given).

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

pub struct PredictUsingTrendRasterTool;

impl Tool for PredictUsingTrendRasterTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "predict_using_trend_raster",
            display_name: "Predict Using Trend Raster",
            summary: "Evaluate the slope/intercept rasters fitted by generate_trend_raster at one or more times, producing predicted surfaces, like ArcGIS Predict Using Trend Raster.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Slope (trend) raster written by generate_trend_raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "intercept",
                    description: "Intercept raster written by generate_trend_raster (its 'intercept_output'), co-registered to the slope raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output raster path for the first predicted time. If omitted, results are stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "times",
                    description: "Comma-separated times to predict at, on the same scale used for the fit (e.g. '12,18,24'). Overrides start/end/interval.",
                    required: false,
                },
                ToolParamSpec {
                    name: "start",
                    description: "First time of a generated range (used when 'times' is absent).",
                    required: false,
                },
                ToolParamSpec {
                    name: "end",
                    description: "Last time of a generated range (inclusive; used when 'times' is absent).",
                    required: false,
                },
                ToolParamSpec {
                    name: "interval",
                    description: "Step between generated times (default 1; used when 'times' is absent).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read from each coefficient raster (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "intercept")?;
        parse_times(args)?;
        parse_band(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let slope_path = require_str(args, "input")?;
        let intercept_path = require_str(args, "intercept")?;
        let output = parse_optional_output(args, "output")?;
        let times = parse_times(args)?;
        let band = parse_band(args)?;

        let slope_r = load_input_raster(slope_path)?;
        let intercept_r = load_input_raster(intercept_path)?;

        if slope_r.rows != intercept_r.rows || slope_r.cols != intercept_r.cols {
            return Err(ToolError::Validation(format!(
                "intercept raster is {}x{}, expected {}x{} to match the slope raster",
                intercept_r.rows, intercept_r.cols, slope_r.rows, slope_r.cols
            )));
        }
        for (label, r) in [("input", &slope_r), ("intercept", &intercept_r)] {
            if band < 0 || band as usize >= r.bands {
                return Err(ToolError::Validation(format!(
                    "band {} out of range for the '{label}' raster ({} band(s))",
                    band + 1,
                    r.bands
                )));
            }
        }

        let (rows, cols) = (slope_r.rows, slope_r.cols);
        let nodata = slope_r.nodata;
        ctx.progress.info(&format!(
            "predicting {} time(s) over {rows}x{cols}",
            times.len()
        ));

        // Read both coefficient bands once; every predicted time reuses them.
        let mut slope = vec![0.0_f64; rows * cols];
        let mut intercept = vec![0.0_f64; rows * cols];
        let mut valid = vec![false; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let idx = r * cols + c;
                let m = slope_r.get(band, r as isize, c as isize);
                let b = intercept_r.get(band, r as isize, c as isize);
                // A cell is predictable only where BOTH coefficients are valid.
                if m != slope_r.nodata && m.is_finite() && b != intercept_r.nodata && b.is_finite()
                {
                    slope[idx] = m;
                    intercept[idx] = b;
                    valid[idx] = true;
                }
            }
        }
        let valid_cells = valid.iter().filter(|v| **v).count();

        let mut written: Vec<String> = Vec::with_capacity(times.len());
        for (i, &t) in times.iter().enumerate() {
            let mut data = vec![nodata; rows * cols];
            for idx in 0..rows * cols {
                if valid[idx] {
                    data[idx] = intercept[idx] + slope[idx] * t;
                }
            }
            let raster = raster_like_with_data(&slope_r, data, nodata, DataType::F32)?;
            let target = match (output, i) {
                (Some(path), 0) => Some(path.to_string()),
                (Some(path), _) => Some(suffixed_path(path, t)),
                (None, _) => None,
            };
            written.push(write_or_store_output(raster, target.as_deref())?);
            ctx.progress.progress((i as f64 + 1.0) / times.len() as f64);
        }

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(written[0]));
        outputs.insert("outputs".to_string(), json!(written));
        outputs.insert("times".to_string(), json!(times));
        outputs.insert("time_count".to_string(), json!(times.len()));
        outputs.insert("predicted_cells".to_string(), json!(valid_cells));
        Ok(ToolRunResult { outputs })
    }
}

/// Inserts a `_t<value>` tag before the extension so multi-time runs do not
/// overwrite one another (`ndvi.tif` -> `ndvi_t18.tif`).
fn suffixed_path(path: &str, t: f64) -> String {
    // Format the tag without a trailing ".0" for whole numbers, and with '.'/'-'
    // replaced so the filename stays portable.
    let tag = if t.fract() == 0.0 {
        format!("{}", t as i64)
    } else {
        format!("{t}")
    }
    .replace(['.', '-'], "_");

    match path.rfind('.') {
        Some(dot) if !path[dot..].contains('/') && !path[dot..].contains('\\') => {
            format!("{}_t{}{}", &path[..dot], tag, &path[dot..])
        }
        _ => format!("{path}_t{tag}"),
    }
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    match args.get(key).and_then(Value::as_str) {
        Some(s) if !s.trim().is_empty() => Ok(s),
        _ => Err(ToolError::Validation(format!(
            "missing required string parameter '{key}'"
        ))),
    }
}

fn parse_band(args: &ToolArgs) -> Result<isize, ToolError> {
    let band_1based = match args.get("band") {
        None | Some(Value::Null) => 1_u64,
        Some(Value::Number(n)) => n
            .as_u64()
            .ok_or_else(|| ToolError::Validation("'band' must be a positive integer".into()))?,
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
            .map_err(|_| ToolError::Validation("'band' must be a positive integer".into()))?,
        Some(_) => {
            return Err(ToolError::Validation(
                "'band' must be a positive integer".into(),
            ))
        }
    };
    Ok((band_1based.max(1) - 1) as isize)
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

/// Resolves the requested prediction times from either an explicit `times` list
/// or a `start`/`end`/`interval` range.
fn parse_times(args: &ToolArgs) -> Result<Vec<f64>, ToolError> {
    if let Some(raw) = args.get("times").and_then(Value::as_str) {
        if !raw.trim().is_empty() {
            let mut out = Vec::new();
            for tok in raw.split(',') {
                let tok = tok.trim();
                if tok.is_empty() {
                    continue;
                }
                out.push(tok.parse::<f64>().map_err(|_| {
                    ToolError::Validation(format!("'times' entry '{tok}' is not a number"))
                })?);
            }
            if out.is_empty() {
                return Err(ToolError::Validation(
                    "'times' did not contain any values".into(),
                ));
            }
            return Ok(out);
        }
    }

    let start = parse_optional_f64(args, "start")?;
    let end = parse_optional_f64(args, "end")?;
    let interval = parse_optional_f64(args, "interval")?.unwrap_or(1.0);

    match (start, end) {
        (Some(s), Some(e)) => {
            if interval <= 0.0 {
                return Err(ToolError::Validation(
                    "'interval' must be greater than 0".into(),
                ));
            }
            if e < s {
                return Err(ToolError::Validation(
                    "'end' must be greater than or equal to 'start'".into(),
                ));
            }
            // Guard against a pathologically small interval producing a huge stack.
            let steps = ((e - s) / interval).floor() as usize;
            if steps > 10_000 {
                return Err(ToolError::Validation(format!(
                    "start/end/interval would generate {} rasters; supply a coarser 'interval'",
                    steps + 1
                )));
            }
            let mut out = Vec::with_capacity(steps + 1);
            for k in 0..=steps {
                out.push(s + interval * k as f64);
            }
            Ok(out)
        }
        (Some(s), None) => Ok(vec![s]),
        _ => Err(ToolError::Validation(
            "supply either 'times' or both 'start' and 'end'".into(),
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{memory_store, DataType, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    const ND: f64 = -9999.0;

    /// Builds a 1-band raster from a row-major buffer and returns a memory path.
    fn raster_of(rows: usize, cols: usize, data: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: ND,
            data_type: DataType::F32,
            crs: wbraster::CrsInfo::default(),
            metadata: Default::default(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, data[row * cols + col])
                    .unwrap();
            }
        }
        let id = memory_store::put_raster(r);
        memory_store::make_raster_memory_path(&id)
    }

    fn read_all(path: &str) -> Vec<f64> {
        let r = load_input_raster(path).unwrap();
        let mut out = Vec::new();
        for row in 0..r.rows {
            for col in 0..r.cols {
                out.push(r.get(0, row as isize, col as isize));
            }
        }
        out
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        PredictUsingTrendRasterTool.run(&args, &ctx()).unwrap()
    }

    /// v(t) = b + m*t evaluated at an explicit time.
    #[test]
    fn predicts_linear_model_at_explicit_time() {
        let slope = raster_of(1, 3, &[2.0, -1.0, 0.5]);
        let intercept = raster_of(1, 3, &[10.0, 100.0, 0.0]);
        let out = run(json!({ "input": slope, "intercept": intercept, "times": "4" }));

        let vals = read_all(out.outputs["output"].as_str().unwrap());
        // 10 + 2*4 = 18 ; 100 + (-1)*4 = 96 ; 0 + 0.5*4 = 2
        assert_eq!(vals, vec![18.0, 96.0, 2.0]);
        assert_eq!(out.outputs["time_count"], json!(1));
        assert_eq!(out.outputs["predicted_cells"], json!(3));
    }

    /// t = 0 must reproduce the intercept exactly.
    #[test]
    fn time_zero_returns_intercept() {
        let slope = raster_of(1, 2, &[3.0, -7.5]);
        let intercept = raster_of(1, 2, &[1.0, 2.0]);
        let out = run(json!({ "input": slope, "intercept": intercept, "times": "0" }));
        assert_eq!(
            read_all(out.outputs["output"].as_str().unwrap()),
            vec![1.0, 2.0]
        );
    }

    /// start/end/interval generates an inclusive, evenly spaced stack.
    #[test]
    fn range_generates_inclusive_series() {
        let slope = raster_of(1, 1, &[1.0]);
        let intercept = raster_of(1, 1, &[0.0]);
        let out = run(json!({
            "input": slope, "intercept": intercept,
            "start": 2, "end": 8, "interval": 2
        }));
        assert_eq!(out.outputs["times"], json!([2.0, 4.0, 6.0, 8.0]));
        let paths = out.outputs["outputs"].as_array().unwrap();
        assert_eq!(paths.len(), 4);
        // slope 1, intercept 0 -> prediction equals the time itself.
        for (i, t) in [2.0, 4.0, 6.0, 8.0].iter().enumerate() {
            assert_eq!(read_all(paths[i].as_str().unwrap()), vec![*t]);
        }
    }

    /// No-data in EITHER coefficient must suppress the cell.
    #[test]
    fn nodata_in_either_coefficient_propagates() {
        let slope = raster_of(1, 3, &[1.0, ND, 1.0]);
        let intercept = raster_of(1, 3, &[0.0, 0.0, ND]);
        let out = run(json!({ "input": slope, "intercept": intercept, "times": "5" }));
        let vals = read_all(out.outputs["output"].as_str().unwrap());
        assert_eq!(vals[0], 5.0);
        assert_eq!(vals[1], ND, "nodata slope must not predict");
        assert_eq!(vals[2], ND, "nodata intercept must not predict");
        assert_eq!(out.outputs["predicted_cells"], json!(1));
    }

    /// Mismatched coefficient grids are a validation error, not a silent crop.
    #[test]
    fn rejects_mismatched_grids() {
        let slope = raster_of(2, 2, &[1.0, 1.0, 1.0, 1.0]);
        let intercept = raster_of(1, 2, &[0.0, 0.0]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": slope, "intercept": intercept, "times": "1" }))
                .unwrap();
        let err = PredictUsingTrendRasterTool.run(&args, &ctx()).unwrap_err();
        assert!(matches!(err, ToolError::Validation(_)));
    }

    #[test]
    fn rejects_bad_parameters() {
        let slope = raster_of(1, 1, &[1.0]);
        let intercept = raster_of(1, 1, &[0.0]);
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            PredictUsingTrendRasterTool.validate(&args).is_err()
        };
        // no times and no start/end
        assert!(bad(json!({ "input": slope, "intercept": intercept })));
        // non-numeric time
        assert!(bad(
            json!({ "input": slope, "intercept": intercept, "times": "abc" })
        ));
        // end before start
        assert!(bad(
            json!({ "input": slope, "intercept": intercept, "start": 5, "end": 1 })
        ));
        // non-positive interval
        assert!(bad(
            json!({ "input": slope, "intercept": intercept, "start": 1, "end": 5, "interval": 0 })
        ));
        // missing intercept
        assert!(bad(json!({ "input": slope, "times": "1" })));
    }

    /// Multi-time file output must not collide on one path.
    #[test]
    fn suffixed_path_disambiguates() {
        assert_eq!(suffixed_path("/tmp/ndvi.tif", 18.0), "/tmp/ndvi_t18.tif");
        assert_eq!(suffixed_path("/tmp/ndvi.tif", -2.5), "/tmp/ndvi_t_2_5.tif");
        assert_eq!(suffixed_path("/tmp/ndvi", 3.0), "/tmp/ndvi_t3");
    }
}
