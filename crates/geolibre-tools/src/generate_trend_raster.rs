//! GeoLibre tool: per-pixel temporal trend across a time series of rasters.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Generate Trend Raster* (Image
//! Analyst) — the workhorse of Earth-observation change monitoring (NDVI
//! greening/browning, temperature trends). The bundled suite has nothing
//! temporal-per-pixel: `trend_surface` is a spatial polynomial,
//! `change_vector_analysis` compares two dates, and `image_stack_profile` only
//! extracts a series at sample points. Pairs with `spectral_index` and
//! `detect_image_anomalies` for an EO monitoring pipeline.
//!
//! Given an ordered stack of co-registered single-band rasters and their times,
//! each pixel's valid observations are fit two ways:
//! * **linear** — ordinary least squares; `output` is the slope, with optional
//!   `intercept_output` and `significance_output` (the R² of the fit).
//! * **mann_kendall** — the non-parametric Mann-Kendall trend test with Sen's
//!   slope; `output` is Sen's slope and `significance_output` is the two-sided
//!   p-value (normal approximation with tie and continuity correction).
//!
//! Pixels with fewer than `min_valid` non-no-data observations become no-data.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::common::{load_input_raster, parse_optional_output, raster_like_with_data};

pub struct GenerateTrendRasterTool;

impl Tool for GenerateTrendRasterTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "generate_trend_raster",
            display_name: "Generate Trend Raster",
            summary: "Fit a per-pixel temporal trend across a time series of co-registered rasters (like ArcGIS Generate Trend Raster): OLS slope/intercept/R² or non-parametric Mann-Kendall + Sen's slope with a p-value — the per-pixel temporal trend the bundled trend_surface (spatial) and change_vector_analysis (two-date) can't produce.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "inputs",
                    description: "Comma-separated list of co-registered single-band raster paths, in time order (>= 3).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output slope raster (Sen's slope for mann_kendall). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "times",
                    description: "Comma-separated numeric times (e.g. years) matching 'inputs'. Default 0,1,2,...",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'linear' (OLS; default) or 'mann_kendall' (non-parametric + Sen's slope).",
                    required: false,
                },
                ToolParamSpec {
                    name: "intercept_output",
                    description: "linear only: optional output raster of the OLS intercept.",
                    required: false,
                },
                ToolParamSpec {
                    name: "significance_output",
                    description: "Optional output raster: R² (linear) or the Mann-Kendall two-sided p-value (mann_kendall).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_valid",
                    description: "Minimum non-no-data observations per pixel (default 3).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read from each input (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        let paths = parse_inputs(args)?;
        if paths.len() < 3 {
            return Err(ToolError::Validation(
                "'inputs' needs at least 3 rasters for a trend".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let paths = parse_inputs(args)?;
        let prm = parse_params(args)?;
        let output = parse_optional_output(args, "output")?;
        let intercept_out = parse_optional_output(args, "intercept_output")?;
        let signif_out = parse_optional_output(args, "significance_output")?;

        let times = match &prm.times {
            Some(t) if t.len() != paths.len() => {
                return Err(ToolError::Validation(format!(
                    "{} times for {} rasters",
                    t.len(),
                    paths.len()
                )))
            }
            Some(t) => t.clone(),
            None => (0..paths.len()).map(|i| i as f64).collect(),
        };

        // Load the stack; all rasters must share dimensions.
        let rasters: Vec<Raster> = paths
            .iter()
            .map(|p| load_input_raster(p))
            .collect::<Result<_, _>>()?;
        let (rows, cols) = (rasters[0].rows, rasters[0].cols);
        let band = prm.band;
        for (i, r) in rasters.iter().enumerate() {
            if r.rows != rows || r.cols != cols {
                return Err(ToolError::Validation(format!(
                    "raster {i} is {}x{}, expected {rows}x{cols}",
                    r.rows, r.cols
                )));
            }
            if band < 0 || band as usize >= r.bands {
                return Err(ToolError::Validation(format!(
                    "band {} out of range for raster {i}",
                    band + 1
                )));
            }
        }

        ctx.progress.info(&format!(
            "{} raster(s), {}x{}, method {}",
            rasters.len(),
            rows,
            cols,
            prm.method.label()
        ));

        let nodata = rasters[0].nodata;
        let mut slope = vec![nodata; rows * cols];
        let mut intercept = vec![nodata; rows * cols];
        let mut signif = vec![nodata; rows * cols];

        let mut ts: Vec<f64> = Vec::with_capacity(rasters.len());
        let mut vs: Vec<f64> = Vec::with_capacity(rasters.len());
        for r in 0..rows {
            for c in 0..cols {
                ts.clear();
                vs.clear();
                for (k, ras) in rasters.iter().enumerate() {
                    let v = ras.get(band, r as isize, c as isize);
                    if v != ras.nodata && v.is_finite() {
                        ts.push(times[k]);
                        vs.push(v);
                    }
                }
                if vs.len() < prm.min_valid {
                    continue;
                }
                let idx = r * cols + c;
                match prm.method {
                    Method::Linear => {
                        if let Some((m, b, r2)) = ols(&ts, &vs) {
                            slope[idx] = m;
                            intercept[idx] = b;
                            signif[idx] = r2;
                        }
                    }
                    Method::MannKendall => {
                        let (sen, p) = mann_kendall(&ts, &vs);
                        slope[idx] = sen;
                        signif[idx] = p;
                    }
                }
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let slope_r = raster_like_with_data(&rasters[0], slope, nodata, DataType::F32)?;
        let out_path = crate::common::write_or_store_output(slope_r, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        if let Some(p) = intercept_out {
            let r = raster_like_with_data(&rasters[0], intercept, nodata, DataType::F32)?;
            outputs.insert(
                "intercept_output".to_string(),
                json!(crate::common::write_or_store_output(r, Some(p))?),
            );
        }
        if let Some(p) = signif_out {
            let r = raster_like_with_data(&rasters[0], signif, nodata, DataType::F32)?;
            outputs.insert(
                "significance_output".to_string(),
                json!(crate::common::write_or_store_output(r, Some(p))?),
            );
        }
        Ok(ToolRunResult { outputs })
    }
}

/// OLS fit: returns (slope, intercept, R²).
fn ols(t: &[f64], v: &[f64]) -> Option<(f64, f64, f64)> {
    let n = t.len() as f64;
    if n < 2.0 {
        return None;
    }
    let mt = t.iter().sum::<f64>() / n;
    let mv = v.iter().sum::<f64>() / n;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut syy = 0.0;
    for i in 0..t.len() {
        let dt = t[i] - mt;
        let dv = v[i] - mv;
        sxx += dt * dt;
        sxy += dt * dv;
        syy += dv * dv;
    }
    if sxx <= 0.0 {
        return None;
    }
    let slope = sxy / sxx;
    let intercept = mv - slope * mt;
    let r2 = if syy > 0.0 {
        (sxy * sxy) / (sxx * syy)
    } else {
        1.0
    };
    Some((slope, intercept, r2))
}

/// Mann-Kendall test with Sen's slope. Returns (sen_slope, two_sided_p).
fn mann_kendall(t: &[f64], v: &[f64]) -> (f64, f64) {
    let n = v.len();
    // S statistic.
    let mut s = 0i64;
    for i in 0..n {
        for j in (i + 1)..n {
            s += (v[j] - v[i]).signum() as i64;
        }
    }
    // Variance with tie correction.
    let mut ties: BTreeMap<u64, i64> = BTreeMap::new();
    for &x in v {
        *ties.entry(x.to_bits()).or_insert(0) += 1;
    }
    let n_f = n as f64;
    let mut var = n_f * (n_f - 1.0) * (2.0 * n_f + 5.0);
    for &tie in ties.values() {
        let tf = tie as f64;
        var -= tf * (tf - 1.0) * (2.0 * tf + 5.0);
    }
    var /= 18.0;
    let z = if var <= 0.0 {
        0.0
    } else if s > 0 {
        (s as f64 - 1.0) / var.sqrt()
    } else if s < 0 {
        (s as f64 + 1.0) / var.sqrt()
    } else {
        0.0
    };
    let p = 2.0 * (1.0 - normal_cdf(z.abs()));

    // Sen's slope: median of pairwise slopes.
    let mut slopes: Vec<f64> = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            if t[j] != t[i] {
                slopes.push((v[j] - v[i]) / (t[j] - t[i]));
            }
        }
    }
    let sen = median(&mut slopes);
    (sen, p.clamp(0.0, 1.0))
}

fn median(v: &mut [f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(|a, b| a.total_cmp(b));
    let m = v.len() / 2;
    if v.len() % 2 == 1 {
        v[m]
    } else {
        (v[m - 1] + v[m]) / 2.0
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

// ── Parameters ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Linear,
    MannKendall,
}

impl Method {
    fn label(&self) -> &'static str {
        match self {
            Method::Linear => "linear",
            Method::MannKendall => "mann_kendall",
        }
    }
}

struct Params {
    method: Method,
    times: Option<Vec<f64>>,
    min_valid: usize,
    band: isize,
}

fn parse_inputs(args: &ToolArgs) -> Result<Vec<String>, ToolError> {
    let s = args
        .get("inputs")
        .and_then(Value::as_str)
        .ok_or_else(|| ToolError::Validation("missing required parameter 'inputs'".to_string()))?;
    let paths: Vec<String> = s
        .split(',')
        .map(str::trim)
        .filter(|p| !p.is_empty())
        .map(String::from)
        .collect();
    if paths.is_empty() {
        return Err(ToolError::Validation("'inputs' is empty".to_string()));
    }
    Ok(paths)
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let method = match args.get("method").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("linear") => Method::Linear,
        Some("mann_kendall") => Method::MannKendall,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'method' must be 'linear' or 'mann_kendall', got '{o}'"
            )))
        }
    };
    let times = match args.get("times").and_then(Value::as_str) {
        None => None,
        Some(s) if s.trim().is_empty() => None,
        Some(s) => Some(
            s.split(',')
                .map(str::trim)
                .filter(|x| !x.is_empty())
                .map(|x| {
                    x.parse::<f64>()
                        .map_err(|_| ToolError::Validation(format!("time '{x}' is not a number")))
                })
                .collect::<Result<Vec<_>, _>>()?,
        ),
    };
    let min_valid = match args.get("min_valid") {
        None | Some(Value::Null) => 3,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(3).max(2) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 3,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'min_valid' must be an integer".into()))?
            .max(2),
        _ => 3,
    };
    let band_1based = match args.get("band") {
        None | Some(Value::Null) => 1,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(1).max(1) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'band' must be an integer".into()))?
            .max(1),
        _ => 1,
    };
    Ok(Params {
        method,
        times,
        min_valid,
        band: (band_1based - 1) as isize,
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

    fn raster_from(cols: usize, rows: usize, data: Vec<f64>) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
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
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, data[row * cols + col])
                    .unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> Raster {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = GenerateTrendRasterTool.run(&args, &ctx()).unwrap();
        load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    /// A pixel with a perfect linear increase recovers the exact slope.
    #[test]
    fn ols_recovers_known_slope() {
        // 1x1 pixel; value = 10 + 2*t over t=0..4.
        let rasters: Vec<String> = (0..5)
            .map(|t| raster_from(1, 1, vec![10.0 + 2.0 * t as f64]))
            .collect();
        let r = run(json!({ "inputs": rasters.join(","), "method": "linear" }));
        assert!((r.get(0, 0, 0) - 2.0).abs() < 1e-9, "slope should be 2");
    }

    /// Mann-Kendall gives Sen's slope; a monotone increase is positive and
    /// significant.
    #[test]
    fn mann_kendall_positive_trend() {
        let vals = [1.0, 2.0, 4.0, 7.0, 11.0, 16.0];
        let rasters: Vec<String> = vals.iter().map(|&v| raster_from(1, 1, vec![v])).collect();
        let r = run(json!({ "inputs": rasters.join(","), "method": "mann_kendall" }));
        assert!(
            r.get(0, 0, 0) > 0.0,
            "monotone increase -> positive Sen slope"
        );
        // The p-value path is exercised directly on the statistic.
        let t: Vec<f64> = (0..6).map(|i| i as f64).collect();
        let (sen, p) = mann_kendall(&t, &vals);
        assert!(sen > 0.0);
        assert!(
            p < 0.1,
            "strong monotone trend should be significant, got p={p}"
        );
    }

    /// Pixels with too few valid observations become no-data.
    #[test]
    fn insufficient_data_is_nodata() {
        // Two rasters mostly nodata; min_valid default 3 -> all nodata.
        let r1 = raster_from(1, 1, vec![5.0]);
        let r2 = raster_from(1, 1, vec![6.0]);
        let r3 = raster_from(1, 1, vec![-9999.0]);
        let r = run(json!({ "inputs": format!("{r1},{r2},{r3}"), "method": "linear" }));
        assert_eq!(r.get(0, 0, 0), -9999.0, "only 2 valid obs < min_valid 3");
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            GenerateTrendRasterTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "inputs": "a.tif,b.tif" })).is_err()); // < 3
        assert!(bad(json!({ "inputs": "a.tif,b.tif,c.tif", "method": "quadratic" })).is_err());
        assert!(bad(json!({ "inputs": "a.tif,b.tif,c.tif" })).is_ok());
    }
}
