//! GeoLibre tool: per-band statistics plus the cross-band covariance and
//! correlation matrices for a stack of co-registered rasters.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Band Collection Statistics* (Spatial
//! Analyst). `cell_statistics` (#411) reduces a raster stack **per pixel**; this
//! summarizes the stack **per band** and, crucially, computes the cross-band
//! covariance and correlation matrices — the standard precursor to PCA /
//! multivariate image classification and to diagnosing band redundancy. Nothing
//! in the bundled suite emits that matrix.
//!
//! Inputs are one multiband raster (each band is a layer) or a comma-separated
//! list of co-registered single/multi-band rasters. Per-band min/max/mean/std use
//! every valid cell of that band; the covariance and correlation matrices use only
//! cells that are valid in **all** bands (list-wise complete observations), so the
//! matrix is internally consistent.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;

use crate::common::{load_input_raster, parse_optional_output, write_text_output};

pub struct BandCollectionStatisticsTool;

impl Tool for BandCollectionStatisticsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "band_collection_statistics",
            display_name: "Band Collection Statistics",
            summary: "Per-band min/max/mean/std plus the cross-band covariance and correlation matrices for a raster stack (like ArcGIS Band Collection Statistics) — the multivariate summary (PCA / classification precursor) the bundled suite lacks; distinct from per-pixel cell_statistics.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "inputs",
                    description: "One multiband raster (each band is a layer) or a comma-separated list of co-registered rasters.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional CSV path for the statistics report. The matrices are always returned in the result.",
                    required: false,
                },
                ToolParamSpec {
                    name: "detail",
                    description: "'detailed' (default) includes covariance/correlation matrices; 'brief' returns per-band statistics only.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        parse_inputs(args)?;
        parse_detail(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let paths = parse_inputs(args)?;
        let output = parse_optional_output(args, "output")?;
        let detailed = parse_detail(args)?;

        let rasters: Vec<Raster> = paths
            .iter()
            .map(|p| load_input_raster(p))
            .collect::<Result<_, _>>()?;
        let (rows, cols) = (rasters[0].rows, rasters[0].cols);
        let base = &rasters[0];
        let aligned = |a: f64, b: f64| (a - b).abs() <= 1e-6 * a.abs().max(b.abs()).max(1.0);
        let mut layers: Vec<(usize, isize)> = Vec::new();
        for (i, r) in rasters.iter().enumerate() {
            if r.rows != rows || r.cols != cols {
                return Err(ToolError::Validation(format!(
                    "raster {i} is {}x{}, expected {rows}x{cols}",
                    r.rows, r.cols
                )));
            }
            if !aligned(r.x_min, base.x_min)
                || !aligned(r.y_min, base.y_min)
                || !aligned(r.cell_size_x, base.cell_size_x)
                || !aligned(r.cell_size_y, base.cell_size_y)
            {
                return Err(ToolError::Validation(format!(
                    "raster {i} is not co-registered with input 0 (origin/resolution differ)"
                )));
            }
            for b in 0..r.bands {
                layers.push((i, b as isize));
            }
        }
        let nb = layers.len();
        if nb < 2 {
            return Err(ToolError::Validation(format!(
                "need at least 2 bands to compute a covariance matrix, got {nb}"
            )));
        }

        ctx.progress
            .info(&format!("{nb} band(s) over {rows}x{cols}"));

        // Per-band running stats over each band's own valid cells.
        let mut count = vec![0u64; nb];
        let mut sum = vec![0.0f64; nb];
        let mut sumsq = vec![0.0f64; nb];
        let mut minv = vec![f64::INFINITY; nb];
        let mut maxv = vec![f64::NEG_INFINITY; nb];

        // Cross-band accumulators over list-wise complete cells (valid in all).
        let mut complete = 0u64;
        let mut csum = vec![0.0f64; nb];
        let mut cross = vec![0.0f64; nb * nb];

        let mut vals = vec![0.0f64; nb];
        for r in 0..rows {
            for c in 0..cols {
                let mut all_valid = true;
                for (k, &(ri, band)) in layers.iter().enumerate() {
                    let ras = &rasters[ri];
                    let v = ras.get(band, r as isize, c as isize);
                    if v != ras.nodata && v.is_finite() {
                        vals[k] = v;
                        count[k] += 1;
                        sum[k] += v;
                        sumsq[k] += v * v;
                        minv[k] = minv[k].min(v);
                        maxv[k] = maxv[k].max(v);
                    } else {
                        vals[k] = f64::NAN;
                        all_valid = false;
                    }
                }
                if all_valid {
                    complete += 1;
                    for i in 0..nb {
                        csum[i] += vals[i];
                        for j in 0..nb {
                            cross[i * nb + j] += vals[i] * vals[j];
                        }
                    }
                }
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let mut mean = vec![f64::NAN; nb];
        let mut std = vec![f64::NAN; nb];
        for k in 0..nb {
            if count[k] > 0 {
                let n = count[k] as f64;
                mean[k] = sum[k] / n;
                std[k] = (sumsq[k] / n - mean[k] * mean[k]).max(0.0).sqrt();
            } else {
                minv[k] = f64::NAN;
                maxv[k] = f64::NAN;
            }
        }

        // Population covariance / correlation over complete cells.
        let mut cov = vec![f64::NAN; nb * nb];
        let mut corr = vec![f64::NAN; nb * nb];
        let mut cmean = vec![f64::NAN; nb];
        if complete > 0 {
            let n = complete as f64;
            for i in 0..nb {
                cmean[i] = csum[i] / n;
            }
            for i in 0..nb {
                for j in 0..nb {
                    cov[i * nb + j] = cross[i * nb + j] / n - cmean[i] * cmean[j];
                }
            }
            for i in 0..nb {
                for j in 0..nb {
                    let denom = (cov[i * nb + i] * cov[j * nb + j]).sqrt();
                    corr[i * nb + j] = if denom > 0.0 {
                        (cov[i * nb + j] / denom).clamp(-1.0, 1.0)
                    } else {
                        f64::NAN
                    };
                }
            }
        }

        // Build outputs.
        let mut outputs = BTreeMap::new();
        outputs.insert("band_count".to_string(), json!(nb));
        outputs.insert("complete_cells".to_string(), json!(complete));
        let per_band: Vec<Value> = (0..nb)
            .map(|k| {
                json!({
                    "band": k + 1,
                    "count": count[k],
                    "min": finite_or_null(minv[k]),
                    "max": finite_or_null(maxv[k]),
                    "mean": finite_or_null(mean[k]),
                    "std": finite_or_null(std[k]),
                })
            })
            .collect();
        outputs.insert("bands".to_string(), json!(per_band));
        if detailed {
            outputs.insert("covariance".to_string(), matrix_json(&cov, nb));
            outputs.insert("correlation".to_string(), matrix_json(&corr, nb));
        }

        if let Some(path) = output {
            let mut csv = String::from("band,count,min,max,mean,std\n");
            for k in 0..nb {
                csv.push_str(&format!(
                    "{},{},{},{},{},{}\n",
                    k + 1,
                    count[k],
                    minv[k],
                    maxv[k],
                    mean[k],
                    std[k]
                ));
            }
            if detailed {
                csv.push_str("\ncovariance\n");
                csv.push_str(&matrix_csv(&cov, nb));
                csv.push_str("\ncorrelation\n");
                csv.push_str(&matrix_csv(&corr, nb));
            }
            write_text_output(&csv, path)?;
            outputs.insert("output".to_string(), json!(path));
        }

        Ok(ToolRunResult { outputs })
    }
}

fn finite_or_null(v: f64) -> Value {
    if v.is_finite() {
        json!(v)
    } else {
        Value::Null
    }
}

fn matrix_json(m: &[f64], n: usize) -> Value {
    let rows: Vec<Vec<Value>> = (0..n)
        .map(|i| (0..n).map(|j| finite_or_null(m[i * n + j])).collect())
        .collect();
    json!(rows)
}

fn matrix_csv(m: &[f64], n: usize) -> String {
    let mut s = String::new();
    for i in 0..n {
        let row: Vec<String> = (0..n).map(|j| m[i * n + j].to_string()).collect();
        s.push_str(&row.join(","));
        s.push('\n');
    }
    s
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

fn parse_detail(args: &ToolArgs) -> Result<bool, ToolError> {
    Ok(
        match args.get("detail").and_then(Value::as_str).map(str::trim) {
            None | Some("") | Some("detailed") => true,
            Some("brief") => false,
            Some(o) => {
                return Err(ToolError::Validation(format!(
                    "'detail' must be 'detailed' or 'brief', got '{o}'"
                )))
            }
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

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
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        BandCollectionStatisticsTool.run(&args, &ctx()).unwrap()
    }

    /// Per-band mean/min/max over a 2-band, 4-cell stack.
    #[test]
    fn per_band_stats() {
        let b0 = vec![1.0, 2.0, 3.0, 4.0];
        let b1 = vec![10.0, 20.0, 30.0, 40.0];
        let out = run(json!({ "inputs": multiband(2, 2, 2, vec![b0, b1]) }));
        let bands = out.outputs["bands"].as_array().unwrap();
        assert_eq!(bands[0]["mean"].as_f64().unwrap(), 2.5);
        assert_eq!(bands[0]["min"].as_f64().unwrap(), 1.0);
        assert_eq!(bands[1]["max"].as_f64().unwrap(), 40.0);
        assert_eq!(bands[1]["mean"].as_f64().unwrap(), 25.0);
    }

    /// Two perfectly linearly-related bands correlate at +1; the covariance
    /// diagonal equals each band's population variance.
    #[test]
    fn perfect_correlation() {
        // band1 = 10 * band0 -> corr = 1.
        let b0 = vec![1.0, 2.0, 3.0, 4.0];
        let b1 = vec![10.0, 20.0, 30.0, 40.0];
        let out = run(json!({ "inputs": multiband(2, 2, 2, vec![b0, b1]) }));
        let corr = out.outputs["correlation"].as_array().unwrap();
        let c01 = corr[0].as_array().unwrap()[1].as_f64().unwrap();
        assert!((c01 - 1.0).abs() < 1e-9, "corr should be 1, got {c01}");
        // var(band0) of 1,2,3,4 (pop) = 1.25.
        let cov = out.outputs["covariance"].as_array().unwrap();
        let v00 = cov[0].as_array().unwrap()[0].as_f64().unwrap();
        assert!(
            (v00 - 1.25).abs() < 1e-9,
            "var band0 should be 1.25, got {v00}"
        );
    }

    /// Anti-correlated bands correlate at -1.
    #[test]
    fn negative_correlation() {
        let b0 = vec![1.0, 2.0, 3.0, 4.0];
        let b1 = vec![4.0, 3.0, 2.0, 1.0];
        let out = run(json!({ "inputs": multiband(2, 2, 2, vec![b0, b1]) }));
        let corr = out.outputs["correlation"].as_array().unwrap();
        let c01 = corr[0].as_array().unwrap()[1].as_f64().unwrap();
        assert!((c01 + 1.0).abs() < 1e-9, "corr should be -1, got {c01}");
    }

    /// brief mode omits the matrices.
    #[test]
    fn brief_omits_matrices() {
        let b0 = vec![1.0, 2.0];
        let b1 = vec![3.0, 4.0];
        let out = run(json!({ "inputs": multiband(2, 1, 2, vec![b0, b1]), "detail": "brief" }));
        assert!(!out.outputs.contains_key("covariance"));
        assert!(out.outputs.contains_key("bands"));
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            BandCollectionStatisticsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "inputs": "a.tif", "detail": "verbose" })).is_err());
        // single-band input is rejected at run time (need >=2 bands).
        let single = multiband(2, 1, 1, vec![vec![1.0, 2.0]]);
        let args: ToolArgs = serde_json::from_value(json!({ "inputs": single })).unwrap();
        assert!(BandCollectionStatisticsTool.run(&args, &ctx()).is_err());
    }
}
