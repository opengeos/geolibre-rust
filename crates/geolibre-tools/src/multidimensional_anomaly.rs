//! GeoLibre tool: per-slice temporal anomaly cube from a raster time stack.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Generate Multidimensional Anomaly*
//! (Image Analyst). For every pixel of a multiband/temporal raster stack, this
//! computes a per-pixel temporal **baseline** (mean/std or median across the
//! slice dimension) and then rewrites every time slice as its deviation from
//! that baseline — emitting an anomaly stack with one output slice per input
//! slice.
//!
//! This is distinct from the neighbouring stack tools:
//! * `detect_image_anomalies` is spatial/spectral (Mahalanobis across bands in
//!   *one* scene), not temporal.
//! * `find_argument_statistics` returns the *position* (slice/date) of an
//!   extreme, not the deviation value.
//! * `generate_trend_raster` fits a trend line; the bundled overlay ops
//!   (`average_overlay`, `standard_deviation_overlay`) collapse the stack to a
//!   single layer and cannot emit an anomaly *per slice*.
//!
//! Inputs are one multiband raster (each band is a time slice) or a
//! comma-separated list of co-registered rasters (every band of every raster is
//! a slice, in order). Calculation methods:
//! * **difference_from_mean** — `value - mean`.
//! * **percent_of_mean** — `100 * value / mean`.
//! * **z_score** — `(value - mean) / std` (population std; 0 when std is 0).
//! * **difference_from_median** — `value - median`.
//!
//! An optional `reference_range` ("start-end", 1-based inclusive) restricts the
//! slices used to build the baseline (e.g. a climatological normal period);
//! anomalies are still emitted for every input slice. Pixels with fewer than
//! `min_valid` non-no-data observations inside the reference range become
//! no-data across all output slices.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster, RasterConfig};

use crate::common::{load_input_raster, parse_optional_output, write_or_store_output};

const OUT_NODATA: f64 = -9999.0;

pub struct MultidimensionalAnomalyTool;

impl Tool for MultidimensionalAnomalyTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "multidimensional_anomaly",
            display_name: "Generate Multidimensional Anomaly",
            summary: "Per-slice temporal anomaly cube from a raster time stack (like ArcGIS Generate Multidimensional Anomaly): for each pixel compute a temporal baseline (mean/std or median across the slice dimension) and rewrite every time slice as its deviation from that baseline — difference-from-mean, percent-of-mean, z-score, or difference-from-median — emitting an anomaly stack with one slice per input slice.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "One multiband raster (each band is a time slice) or a comma-separated list of co-registered rasters, in time order.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output multiband anomaly raster (one band per input slice). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "Anomaly calculation: 'difference_from_mean' (default), 'percent_of_mean', 'z_score', or 'difference_from_median'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "reference_range",
                    description: "Optional 1-based inclusive slice range 'start-end' used to compute the per-pixel baseline (default: all slices).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_valid",
                    description: "Minimum non-no-data observations per pixel inside the reference range (default 1); below this the pixel is no-data.",
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

        // Load the stack; flatten every band of every input into an ordered
        // list of (raster_index, band) slices.
        let rasters: Vec<Raster> = paths
            .iter()
            .map(|p| load_input_raster(p))
            .collect::<Result<_, _>>()?;
        let (rows, cols) = (rasters[0].rows, rasters[0].cols);
        let mut slices: Vec<(usize, isize)> = Vec::new();
        for (i, r) in rasters.iter().enumerate() {
            if r.rows != rows || r.cols != cols {
                return Err(ToolError::Validation(format!(
                    "raster {i} is {}x{}, expected {rows}x{cols}",
                    r.rows, r.cols
                )));
            }
            for b in 0..r.bands {
                slices.push((i, b as isize));
            }
        }
        let n_slices = slices.len();
        if n_slices < 2 {
            return Err(ToolError::Validation(format!(
                "need at least 2 slices (bands across inputs), got {n_slices}"
            )));
        }

        // Resolve the reference range (0-based, exclusive end) against the stack.
        let (ref_start, ref_end) = match prm.reference_range {
            Some((s, e)) => {
                if s < 1 || e < s || e > n_slices {
                    return Err(ToolError::Validation(format!(
                        "reference_range {s}-{e} out of bounds for {n_slices} slices"
                    )));
                }
                (s - 1, e)
            }
            None => (0, n_slices),
        };

        ctx.progress.info(&format!(
            "{} input(s), {n_slices} slice(s), {}x{}, method {}, reference slices {}..{}",
            rasters.len(),
            rows,
            cols,
            prm.method.label(),
            ref_start + 1,
            ref_end
        ));

        // One output band per input slice.
        let mut out = Raster::new(RasterConfig {
            cols,
            rows,
            bands: n_slices,
            x_min: rasters[0].x_min,
            y_min: rasters[0].y_min,
            cell_size: rasters[0].cell_size_x,
            cell_size_y: Some(rasters[0].cell_size_y),
            nodata: OUT_NODATA,
            data_type: DataType::F32,
            crs: rasters[0].crs.clone(),
            metadata: rasters[0].metadata.clone(),
        });

        // Reusable per-pixel buffer of valid reference values.
        let mut ref_vals: Vec<f64> = Vec::with_capacity(n_slices);
        for r in 0..rows {
            for c in 0..cols {
                // Collect valid observations inside the reference range.
                ref_vals.clear();
                for &(ri, band) in &slices[ref_start..ref_end] {
                    let ras = &rasters[ri];
                    let v = ras.get(band, r as isize, c as isize);
                    if v != ras.nodata && v.is_finite() {
                        ref_vals.push(v);
                    }
                }

                if ref_vals.len() < prm.min_valid {
                    // Leave every output band as no-data for this pixel.
                    for b in 0..n_slices {
                        out.set(b as isize, r as isize, c as isize, OUT_NODATA)
                            .map_err(|e| ToolError::Execution(format!("write cell: {e}")))?;
                    }
                    continue;
                }

                let baseline = prm.method.baseline(&mut ref_vals);

                // Emit an anomaly for each input slice (over the full stack, not
                // just the reference range).
                for (si, &(ri, band)) in slices.iter().enumerate() {
                    let ras = &rasters[ri];
                    let v = ras.get(band, r as isize, c as isize);
                    let out_v = if v != ras.nodata && v.is_finite() {
                        prm.method.anomaly(v, &baseline)
                    } else {
                        OUT_NODATA
                    };
                    out.set(si as isize, r as isize, c as isize, out_v)
                        .map_err(|e| ToolError::Execution(format!("write cell: {e}")))?;
                }
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let out_path = write_or_store_output(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("slices".to_string(), json!(n_slices));
        outputs.insert("method".to_string(), json!(prm.method.label()));
        outputs.insert("reference_start".to_string(), json!(ref_start + 1));
        outputs.insert("reference_end".to_string(), json!(ref_end));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    DifferenceFromMean,
    PercentOfMean,
    ZScore,
    DifferenceFromMedian,
}

/// Per-pixel baseline statistics computed from the reference slices.
struct Baseline {
    center: f64,
    std: f64,
}

impl Method {
    fn label(&self) -> &'static str {
        match self {
            Method::DifferenceFromMean => "difference_from_mean",
            Method::PercentOfMean => "percent_of_mean",
            Method::ZScore => "z_score",
            Method::DifferenceFromMedian => "difference_from_median",
        }
    }

    fn is_median(&self) -> bool {
        matches!(self, Method::DifferenceFromMedian)
    }

    /// Computes the baseline from the valid reference values. `vals` may be
    /// reordered (median sorts in place).
    fn baseline(&self, vals: &mut [f64]) -> Baseline {
        if self.is_median() {
            Baseline {
                center: median(vals),
                std: 0.0,
            }
        } else {
            let n = vals.len() as f64;
            let mean = vals.iter().sum::<f64>() / n;
            let var = if vals.len() > 1 {
                vals.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / n
            } else {
                0.0
            };
            Baseline {
                center: mean,
                std: var.sqrt(),
            }
        }
    }

    /// Transforms one slice value into its anomaly against the baseline.
    fn anomaly(&self, v: f64, b: &Baseline) -> f64 {
        match self {
            Method::DifferenceFromMean | Method::DifferenceFromMedian => v - b.center,
            Method::PercentOfMean => {
                if b.center == 0.0 {
                    OUT_NODATA
                } else {
                    100.0 * v / b.center
                }
            }
            Method::ZScore => {
                if b.std == 0.0 {
                    0.0
                } else {
                    (v - b.center) / b.std
                }
            }
        }
    }
}

/// Median of a slice (sorts in place); lower-middle for even counts.
fn median(vals: &mut [f64]) -> f64 {
    vals.sort_by(f64::total_cmp);
    let n = vals.len();
    if n % 2 == 1 {
        vals[n / 2]
    } else {
        0.5 * (vals[n / 2 - 1] + vals[n / 2])
    }
}

struct Params {
    method: Method,
    reference_range: Option<(usize, usize)>,
    min_valid: usize,
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

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let method = match args.get("method").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("difference_from_mean") => Method::DifferenceFromMean,
        Some("percent_of_mean") => Method::PercentOfMean,
        Some("z_score") => Method::ZScore,
        Some("difference_from_median") => Method::DifferenceFromMedian,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'method' must be one of difference_from_mean|percent_of_mean|z_score|difference_from_median, got '{o}'"
            )))
        }
    };

    let reference_range = match args.get("reference_range").and_then(Value::as_str) {
        None => None,
        Some(s) if s.trim().is_empty() => None,
        Some(s) => Some(parse_range(s)?),
    };

    let min_valid = match args.get("min_valid") {
        None | Some(Value::Null) => 1,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(1).max(1) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'min_valid' must be an integer".into()))?
            .max(1),
        _ => {
            return Err(ToolError::Validation(
                "'min_valid' must be an integer".into(),
            ))
        }
    };

    Ok(Params {
        method,
        reference_range,
        min_valid,
    })
}

/// Parses a 1-based inclusive "start-end" range (e.g. "1-12").
fn parse_range(s: &str) -> Result<(usize, usize), ToolError> {
    let parts: Vec<&str> = s.split(['-', ':']).map(str::trim).collect();
    if parts.len() != 2 {
        return Err(ToolError::Validation(format!(
            "'reference_range' must be 'start-end', got '{s}'"
        )));
    }
    let start = parts[0]
        .parse::<usize>()
        .map_err(|_| ToolError::Validation(format!("bad reference_range start '{}'", parts[0])))?;
    let end = parts[1]
        .parse::<usize>()
        .map_err(|_| ToolError::Validation(format!("bad reference_range end '{}'", parts[1])))?;
    if start < 1 || end < start {
        return Err(ToolError::Validation(format!(
            "'reference_range' must have 1 <= start <= end, got '{s}'"
        )));
    }
    Ok((start, end))
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{memory_store, CrsInfo};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds an in-memory multiband raster from a per-band, row-major buffer.
    fn multiband(cols: usize, rows: usize, bands: usize, data: Vec<Vec<f64>>) -> String {
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
        for (b, band) in data.iter().enumerate() {
            for row in 0..rows {
                for col in 0..cols {
                    r.set(
                        b as isize,
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

    fn run(args: serde_json::Value) -> Raster {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = MultidimensionalAnomalyTool.run(&args, &ctx()).unwrap();
        load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    /// difference_from_mean: each output slice equals value minus the per-pixel
    /// temporal mean, and the anomalies sum to ~0 across slices.
    #[test]
    fn difference_from_mean_is_value_minus_mean() {
        // single pixel, 4 slices: values 2,4,6,8 -> mean 5.
        let bands = vec![vec![2.0], vec![4.0], vec![6.0], vec![8.0]];
        let path = multiband(1, 1, 4, bands);
        let r = run(json!({ "input": path, "method": "difference_from_mean" }));
        assert_eq!(r.bands, 4, "one output band per input slice");
        assert_eq!(r.get(0, 0, 0), -3.0);
        assert_eq!(r.get(1, 0, 0), -1.0);
        assert_eq!(r.get(2, 0, 0), 1.0);
        assert_eq!(r.get(3, 0, 0), 3.0);
        let sum: f64 = (0..4).map(|b| r.get(b, 0, 0)).sum();
        assert!(sum.abs() < 1e-9, "difference-from-mean anomalies sum to 0");
    }

    /// z_score: standardized deviations; a flat series yields all zeros.
    #[test]
    fn z_score_standardizes_and_flat_is_zero() {
        // values 2,4,6,8 -> mean 5, population std = sqrt(5) ~= 2.23607.
        let bands = vec![vec![2.0], vec![4.0], vec![6.0], vec![8.0]];
        let path = multiband(1, 1, 4, bands);
        let r = run(json!({ "input": path, "method": "z_score" }));
        let std = 5.0_f64.sqrt();
        assert!((r.get(0, 0, 0) - (-3.0 / std)).abs() < 1e-6);
        assert!((r.get(3, 0, 0) - (3.0 / std)).abs() < 1e-6);

        // constant series -> std 0 -> all z-scores 0.
        let flat = multiband(1, 1, 3, vec![vec![7.0], vec![7.0], vec![7.0]]);
        let rf = run(json!({ "input": flat, "method": "z_score" }));
        for b in 0..3 {
            assert_eq!(rf.get(b, 0, 0), 0.0);
        }
    }

    /// difference_from_median: baseline is the median, not the mean.
    #[test]
    fn difference_from_median_uses_median() {
        // values 1,2,100 -> median 2 (robust to the outlier).
        let bands = vec![vec![1.0], vec![2.0], vec![100.0]];
        let path = multiband(1, 1, 3, bands);
        let r = run(json!({ "input": path, "method": "difference_from_median" }));
        assert_eq!(r.get(0, 0, 0), -1.0);
        assert_eq!(r.get(1, 0, 0), 0.0);
        assert_eq!(r.get(2, 0, 0), 98.0);
    }

    /// percent_of_mean: value as a percentage of the temporal mean.
    #[test]
    fn percent_of_mean() {
        // values 5,15 -> mean 10 -> 50%, 150%.
        let bands = vec![vec![5.0], vec![15.0]];
        let path = multiband(1, 1, 2, bands);
        let r = run(json!({ "input": path, "method": "percent_of_mean" }));
        assert_eq!(r.get(0, 0, 0), 50.0);
        assert_eq!(r.get(1, 0, 0), 150.0);
    }

    /// reference_range restricts the baseline to a subset of slices, but
    /// anomalies are still emitted for every input slice.
    #[test]
    fn reference_range_restricts_baseline() {
        // 4 slices: 10,20,100,200. Baseline from slices 1-2 -> mean 15.
        let bands = vec![vec![10.0], vec![20.0], vec![100.0], vec![200.0]];
        let path = multiband(1, 1, 4, bands);
        let r = run(
            json!({ "input": path, "method": "difference_from_mean", "reference_range": "1-2" }),
        );
        assert_eq!(r.bands, 4, "still 4 output slices");
        assert_eq!(r.get(0, 0, 0), -5.0); // 10 - 15
        assert_eq!(r.get(1, 0, 0), 5.0); // 20 - 15
        assert_eq!(r.get(2, 0, 0), 85.0); // 100 - 15
        assert_eq!(r.get(3, 0, 0), 185.0); // 200 - 15
    }

    /// A pixel whose reference observations are all no-data becomes no-data
    /// across every output slice (pass-through of the missing-data condition).
    #[test]
    fn all_nodata_pixel_stays_nodata() {
        // 2 pixels, 2 slices. Pixel 0 valid, pixel 1 all no-data.
        let bands = vec![vec![4.0, -9999.0], vec![6.0, -9999.0]];
        let path = multiband(2, 1, 2, bands);
        let r = run(json!({ "input": path, "method": "difference_from_mean" }));
        assert_eq!(r.get(0, 0, 0), -1.0); // pixel 0, slice 0: 4 - 5
        assert_eq!(r.get(0, 0, 1), -9999.0); // pixel 1 no-data
        assert_eq!(r.get(1, 0, 1), -9999.0);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            MultidimensionalAnomalyTool.validate(&args)
        };
        assert!(bad(json!({})).is_err()); // missing input
        assert!(bad(json!({ "input": "a.tif", "method": "bogus" })).is_err());
        assert!(bad(json!({ "input": "a.tif", "reference_range": "3" })).is_err());
        assert!(bad(json!({ "input": "a.tif", "reference_range": "5-2" })).is_err());
        assert!(bad(json!({ "input": "a.tif", "method": "z_score" })).is_ok());
        assert!(bad(
            json!({ "input": "a.tif", "method": "percent_of_mean", "reference_range": "1-3" })
        )
        .is_ok());
    }
}
