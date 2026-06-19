//! DEM smoothing filters, ported from `lidar`'s `filtering.py`
//! (`MeanFilter`, `MedianFilter`, `GaussianFilter`).
//!
//! These mirror the `scipy.ndimage` routines the Python package uses: a box
//! `convolve` (mean), `median_filter`, and a separable `gaussian_filter`, all
//! with the SciPy default `reflect` boundary mode. As in the Python code the
//! filter runs over raw cell values; the no-data tag is preserved on the output.

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};

use crate::common::{
    band_to_vec, load_input_raster, parse_optional_output, raster_like_with_data,
    write_or_store_output,
};

/// Smooths a DEM with a mean, median, or Gaussian filter.
pub struct DemFilterTool;

impl Tool for DemFilterTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "dem_filter",
            display_name: "DEM Filter",
            summary: "Smooth a DEM with a mean, median, or Gaussian filter (lidar filtering.py).",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input DEM raster file path.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output raster path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "filter",
                    description: "Filter type: 'mean', 'median', or 'gaussian' (default 'mean').",
                    required: false,
                },
                ToolParamSpec {
                    name: "kernel_size",
                    description: "Window size in pixels for mean/median filters (default 3).",
                    required: false,
                },
                ToolParamSpec {
                    name: "sigma",
                    description: "Standard deviation for the Gaussian filter (default 1.0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to filter (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args.get("input").and_then(Value::as_str).is_none() {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        if let Some(f) = args.get("filter").and_then(Value::as_str) {
            match f.to_lowercase().as_str() {
                "mean" | "median" | "gaussian" => {}
                other => {
                    return Err(ToolError::Validation(format!(
                        "unknown filter '{other}' (expected mean, median, or gaussian)"
                    )))
                }
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args
            .get("input")
            .and_then(Value::as_str)
            .ok_or_else(|| ToolError::Validation("missing required parameter 'input'".to_string()))?;
        let output = parse_optional_output(args, "output")?;
        let filter = args
            .get("filter")
            .and_then(Value::as_str)
            .unwrap_or("mean")
            .to_lowercase();
        let kernel_size = args
            .get("kernel_size")
            .and_then(Value::as_u64)
            .unwrap_or(3)
            .max(1) as usize;
        let sigma = args.get("sigma").and_then(Value::as_f64).unwrap_or(1.0);
        let band_1based = args.get("band").and_then(Value::as_u64).unwrap_or(1).max(1);
        let band = (band_1based - 1) as isize;

        let raster = load_input_raster(input)?;
        if band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "band {band_1based} out of range (raster has {} band(s))",
                raster.bands
            )));
        }

        let rows = raster.rows;
        let cols = raster.cols;
        let nodata = raster.nodata;
        let src = band_to_vec(&raster, band);

        ctx.progress.info(&format!("applying {filter} filter"));
        let filtered = match filter.as_str() {
            "mean" => mean_filter(&src, rows, cols, kernel_size),
            "median" => median_filter(&src, rows, cols, kernel_size),
            "gaussian" => gaussian_filter(&src, rows, cols, sigma),
            _ => unreachable!("validate() restricts the filter type"),
        };
        ctx.progress.progress(1.0);

        let out_raster = raster_like_with_data(&raster, filtered, nodata, raster.data_type)?;
        let out_path = write_or_store_output(out_raster, output)?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("filter".to_string(), json!(filter));
        Ok(ToolRunResult { outputs })
    }
}

/// Reflects an index into `[0, n)` using SciPy's `reflect` mode (edge value
/// duplicated): `-1 -> 0`, `n -> n-1`.
fn reflect(mut i: isize, n: usize) -> usize {
    if n == 1 {
        return 0;
    }
    let n = n as isize;
    loop {
        if i < 0 {
            i = -1 - i;
        } else if i >= n {
            i = 2 * n - 1 - i;
        } else {
            break;
        }
    }
    i as usize
}

/// Box mean filter (`scipy.ndimage.convolve` with a uniform `k x k` kernel).
fn mean_filter(src: &[f64], rows: usize, cols: usize, k: usize) -> Vec<f64> {
    let mut out = vec![0.0_f64; rows * cols];
    let radius = (k / 2) as isize;
    let norm = 1.0 / (k * k) as f64;
    for r in 0..rows {
        for c in 0..cols {
            let mut acc = 0.0;
            for dr in -radius..=radius {
                let rr = reflect(r as isize + dr, rows);
                for dc in -radius..=radius {
                    let cc = reflect(c as isize + dc, cols);
                    acc += src[rr * cols + cc];
                }
            }
            out[r * cols + c] = acc * norm;
        }
    }
    out
}

/// Median filter (`scipy.ndimage.median_filter`, `k x k` window, reflect mode).
fn median_filter(src: &[f64], rows: usize, cols: usize, k: usize) -> Vec<f64> {
    let mut out = vec![0.0_f64; rows * cols];
    let radius = (k / 2) as isize;
    let mut window: Vec<f64> = Vec::with_capacity(k * k);
    for r in 0..rows {
        for c in 0..cols {
            window.clear();
            for dr in -radius..=radius {
                let rr = reflect(r as isize + dr, rows);
                for dc in -radius..=radius {
                    let cc = reflect(c as isize + dc, cols);
                    window.push(src[rr * cols + cc]);
                }
            }
            window.sort_by(|a, b| a.total_cmp(b));
            // SciPy returns the lower-middle element for even-sized windows.
            out[r * cols + c] = window[(window.len() - 1) / 2];
        }
    }
    out
}

/// Separable Gaussian filter (`scipy.ndimage.gaussian_filter`, `truncate=4.0`,
/// reflect mode).
fn gaussian_filter(src: &[f64], rows: usize, cols: usize, sigma: f64) -> Vec<f64> {
    if sigma <= 0.0 {
        return src.to_vec();
    }
    let radius = (4.0 * sigma + 0.5) as isize;
    let mut kernel = vec![0.0_f64; (2 * radius + 1) as usize];
    let mut sum = 0.0;
    for (j, w) in kernel.iter_mut().enumerate() {
        let x = j as isize - radius;
        let val = (-(x * x) as f64 / (2.0 * sigma * sigma)).exp();
        *w = val;
        sum += val;
    }
    for w in &mut kernel {
        *w /= sum;
    }

    // Horizontal pass.
    let mut tmp = vec![0.0_f64; rows * cols];
    for r in 0..rows {
        for c in 0..cols {
            let mut acc = 0.0;
            for (j, &w) in kernel.iter().enumerate() {
                let dc = j as isize - radius;
                let cc = reflect(c as isize + dc, cols);
                acc += w * src[r * cols + cc];
            }
            tmp[r * cols + c] = acc;
        }
    }
    // Vertical pass.
    let mut out = vec![0.0_f64; rows * cols];
    for r in 0..rows {
        for c in 0..cols {
            let mut acc = 0.0;
            for (j, &w) in kernel.iter().enumerate() {
                let dr = j as isize - radius;
                let rr = reflect(r as isize + dr, rows);
                acc += w * tmp[rr * cols + c];
            }
            out[r * cols + c] = acc;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mean_of_constant_is_constant() {
        let src = vec![5.0; 9];
        let out = mean_filter(&src, 3, 3, 3);
        for v in out {
            assert!((v - 5.0).abs() < 1e-9);
        }
    }

    #[test]
    fn median_removes_spike() {
        // A single high spike in the center is replaced by the neighborhood median.
        let mut src = vec![1.0; 9];
        src[4] = 100.0;
        let out = median_filter(&src, 3, 3, 3);
        assert_eq!(out[4], 1.0);
    }

    #[test]
    fn gaussian_preserves_constant() {
        let src = vec![7.0; 25];
        let out = gaussian_filter(&src, 5, 5, 1.0);
        for v in out {
            assert!((v - 7.0).abs() < 1e-6);
        }
    }

    #[test]
    fn reflect_indexing() {
        assert_eq!(reflect(-1, 5), 0);
        assert_eq!(reflect(-2, 5), 1);
        assert_eq!(reflect(5, 5), 4);
        assert_eq!(reflect(6, 5), 3);
    }
}
