//! GeoLibre tool: per-pixel statistic across a stack of aligned rasters.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Cell Statistics* (Spatial Analyst /
//! Image Analyst). For each cell it reduces the values at that position across
//! all input rasters (and all their bands) into a single output value.
//!
//! The repo has neighbourhood and zonal reducers but no plain multi-raster
//! *local* reducer:
//! * `focal_statistics` reduces a moving window **within one** raster;
//! * `zonal_histogram` reduces per zone, not per cell;
//! * `find_argument_statistics` returns the argmax/argmin **slice index**, not
//!   the statistic value, and only that family.
//!
//! In the bundled whitebox suite `sum_overlay` / `weighted_sum` only sum; none
//! gives median / majority / minority / variety / percentile across a stack.
//!
//! Inputs are one multiband raster (each band is a layer) or a comma-separated
//! list of co-registered rasters (every band of every raster is a layer).
//! Statistics: `mean` (default), `majority`, `maximum`, `median`, `minimum`,
//! `minority`, `percentile`, `range`, `std` (population), `sum`, `variety`.
//! `ignore_nodata` (default true) skips no-data observations per cell; when
//! false a single no-data observation makes the whole cell no-data. Cells with
//! no valid observations are no-data.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::common::{load_input_raster, parse_optional_output, raster_like_with_data};

pub struct CellStatisticsTool;

impl Tool for CellStatisticsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "cell_statistics",
            display_name: "Cell Statistics",
            summary: "Per-pixel statistic across a stack of aligned rasters (like ArcGIS Cell Statistics): mean, majority, maximum, median, minimum, minority, percentile, range, std, sum, or variety — the local multi-raster reducer the bundled sum_overlay/weighted_sum (sum only) and the neighbourhood/zonal tools don't provide.",
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
                    description: "Output single-band raster of the reduced value. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "statistic",
                    description: "'mean' (default), 'majority', 'maximum', 'median', 'minimum', 'minority', 'percentile', 'range', 'std', 'sum', or 'variety'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "ignore_nodata",
                    description: "Skip no-data observations per cell (default true). When false, any no-data observation makes the cell no-data.",
                    required: false,
                },
                ToolParamSpec {
                    name: "percentile_value",
                    description: "Percentile in [0, 100] for statistic 'percentile' (default 50 = median).",
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
        // list of (raster_index, band) layers.
        let rasters: Vec<Raster> = paths
            .iter()
            .map(|p| load_input_raster(p))
            .collect::<Result<_, _>>()?;
        let (rows, cols) = (rasters[0].rows, rasters[0].cols);
        let mut layers: Vec<(usize, isize)> = Vec::new();
        for (i, r) in rasters.iter().enumerate() {
            if r.rows != rows || r.cols != cols {
                return Err(ToolError::Validation(format!(
                    "raster {i} is {}x{}, expected {rows}x{cols}",
                    r.rows, r.cols
                )));
            }
            for b in 0..r.bands {
                layers.push((i, b as isize));
            }
        }
        let n_layers = layers.len();
        if n_layers < 2 {
            return Err(ToolError::Validation(format!(
                "need at least 2 layers (bands across inputs), got {n_layers}"
            )));
        }

        ctx.progress.info(&format!(
            "{} input(s), {n_layers} layer(s), {}x{}, statistic {}",
            rasters.len(),
            rows,
            cols,
            prm.statistic.label()
        ));

        let nodata = -9999.0_f64;
        let mut out = vec![nodata; rows * cols];

        // Reusable per-pixel buffer of valid observations.
        let mut vals: Vec<f64> = Vec::with_capacity(n_layers);
        for r in 0..rows {
            for c in 0..cols {
                vals.clear();
                let mut had_nodata = false;
                for &(ri, band) in &layers {
                    let ras = &rasters[ri];
                    let v = ras.get(band, r as isize, c as isize);
                    if v != ras.nodata && v.is_finite() {
                        vals.push(v);
                    } else {
                        had_nodata = true;
                    }
                }
                if vals.is_empty() || (!prm.ignore_nodata && had_nodata) {
                    continue;
                }
                out[r * cols + c] = reduce(&mut vals, &prm);
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let out_r = raster_like_with_data(&rasters[0], out, nodata, DataType::F32)?;
        let out_path = crate::common::write_or_store_output(out_r, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("layers".to_string(), json!(n_layers));
        outputs.insert("statistic".to_string(), json!(prm.statistic.label()));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        Ok(ToolRunResult { outputs })
    }
}

/// Reduces a cell's valid observations to a single value. `vals` may be
/// reordered (median/percentile sort it in place).
fn reduce(vals: &mut [f64], prm: &Params) -> f64 {
    match prm.statistic {
        Statistic::Mean => vals.iter().sum::<f64>() / vals.len() as f64,
        Statistic::Sum => vals.iter().sum(),
        Statistic::Maximum => vals.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        Statistic::Minimum => vals.iter().copied().fold(f64::INFINITY, f64::min),
        Statistic::Range => {
            let mn = vals.iter().copied().fold(f64::INFINITY, f64::min);
            let mx = vals.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            mx - mn
        }
        Statistic::Std => {
            let n = vals.len() as f64;
            let mean = vals.iter().sum::<f64>() / n;
            let var = vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
            var.sqrt()
        }
        Statistic::Median => percentile(vals, 50.0),
        Statistic::Percentile => percentile(vals, prm.percentile_value),
        Statistic::Variety => variety(vals),
        Statistic::Majority => extreme_frequency(vals, true),
        Statistic::Minority => extreme_frequency(vals, false),
    }
}

/// Linear-interpolated percentile (matches numpy's default) after sorting.
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

/// Number of distinct values (exact bit-equality on the f64s).
fn variety(vals: &[f64]) -> f64 {
    let mut seen: Vec<u64> = vals.iter().map(|v| v.to_bits()).collect();
    seen.sort_unstable();
    seen.dedup();
    seen.len() as f64
}

/// Most (`want_max`) or least (`!want_max`) frequent value. Ties resolve to the
/// smaller value, matching ArcGIS's deterministic behaviour.
fn extreme_frequency(vals: &[f64], want_max: bool) -> f64 {
    let mut counts: BTreeMap<u64, (usize, f64)> = BTreeMap::new();
    for &v in vals {
        let e = counts.entry(v.to_bits()).or_insert((0, v));
        e.0 += 1;
    }
    let mut best_count = if want_max { 0usize } else { usize::MAX };
    let mut best_val = f64::NAN;
    // BTreeMap iterates by key ascending; since duplicate values share a key,
    // the first value hit at the extreme count is the smallest — the tie rule.
    let mut ordered: Vec<(f64, usize)> = counts.values().map(|(c, v)| (*v, *c)).collect();
    ordered.sort_by(|a, b| a.0.total_cmp(&b.0));
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
    Majority,
    Maximum,
    Median,
    Minimum,
    Minority,
    Percentile,
    Range,
    Std,
    Sum,
    Variety,
}

impl Statistic {
    fn label(&self) -> &'static str {
        match self {
            Statistic::Mean => "mean",
            Statistic::Majority => "majority",
            Statistic::Maximum => "maximum",
            Statistic::Median => "median",
            Statistic::Minimum => "minimum",
            Statistic::Minority => "minority",
            Statistic::Percentile => "percentile",
            Statistic::Range => "range",
            Statistic::Std => "std",
            Statistic::Sum => "sum",
            Statistic::Variety => "variety",
        }
    }
}

struct Params {
    statistic: Statistic,
    ignore_nodata: bool,
    percentile_value: f64,
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
    let statistic = match args.get("statistic").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("mean") => Statistic::Mean,
        Some("majority") => Statistic::Majority,
        Some("maximum") | Some("max") => Statistic::Maximum,
        Some("median") => Statistic::Median,
        Some("minimum") | Some("min") => Statistic::Minimum,
        Some("minority") => Statistic::Minority,
        Some("percentile") => Statistic::Percentile,
        Some("range") => Statistic::Range,
        Some("std") | Some("stddev") | Some("standard_deviation") => Statistic::Std,
        Some("sum") => Statistic::Sum,
        Some("variety") => Statistic::Variety,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'statistic' must be one of mean|majority|maximum|median|minimum|minority|percentile|range|std|sum|variety, got '{o}'"
            )))
        }
    };

    let ignore_nodata = parse_optional_bool(args, "ignore_nodata")?.unwrap_or(true);

    let percentile_value = match parse_optional_f64(args, "percentile_value")? {
        None => 50.0,
        Some(p) if (0.0..=100.0).contains(&p) => p,
        Some(_) => {
            return Err(ToolError::Validation(
                "'percentile_value' must be in [0, 100]".to_string(),
            ))
        }
    };

    Ok(Params {
        statistic,
        ignore_nodata,
        percentile_value,
    })
}

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
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
        Some(Value::Number(n)) => Ok(Some(n.as_f64().unwrap_or(0.0) != 0.0)),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
}

fn parse_optional_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(Some(n.as_f64().unwrap_or(f64::NAN))),
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
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> Raster {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CellStatisticsTool.run(&args, &ctx()).unwrap();
        load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    /// mean / sum / max / min / range over a 4-layer stack for one pixel.
    #[test]
    fn basic_reducers() {
        // 1x1 raster, values 2,4,6,8 across 4 bands.
        let bands = vec![vec![2.0], vec![4.0], vec![6.0], vec![8.0]];
        let path = multiband(1, 1, 4, bands);
        assert_eq!(
            run(json!({ "inputs": path.clone(), "statistic": "mean" })).get(0, 0, 0),
            5.0
        );
        assert_eq!(
            run(json!({ "inputs": path.clone(), "statistic": "sum" })).get(0, 0, 0),
            20.0
        );
        assert_eq!(
            run(json!({ "inputs": path.clone(), "statistic": "maximum" })).get(0, 0, 0),
            8.0
        );
        assert_eq!(
            run(json!({ "inputs": path.clone(), "statistic": "minimum" })).get(0, 0, 0),
            2.0
        );
        assert_eq!(
            run(json!({ "inputs": path, "statistic": "range" })).get(0, 0, 0),
            6.0
        );
    }

    /// median and percentile interpolate like numpy.
    #[test]
    fn median_and_percentile() {
        // values 10,20,30,40,50 -> median 30; 25th pct = 20.
        let bands = vec![vec![10.0], vec![20.0], vec![30.0], vec![40.0], vec![50.0]];
        let path = multiband(1, 1, 5, bands);
        assert_eq!(
            run(json!({ "inputs": path.clone(), "statistic": "median" })).get(0, 0, 0),
            30.0
        );
        assert_eq!(
            run(json!({ "inputs": path, "statistic": "percentile", "percentile_value": 25 }))
                .get(0, 0, 0),
            20.0
        );
    }

    /// std is the population standard deviation.
    #[test]
    fn population_std() {
        // values 2,4,4,4,5,5,7,9 -> mean 5, pop std 2.
        let bands: Vec<Vec<f64>> = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]
            .iter()
            .map(|&v| vec![v])
            .collect();
        let path = multiband(1, 1, 8, bands);
        let s = run(json!({ "inputs": path, "statistic": "std" })).get(0, 0, 0);
        assert!((s - 2.0).abs() < 1e-9, "pop std {s} != 2");
    }

    /// majority / minority / variety on a categorical stack.
    #[test]
    fn categorical_reducers() {
        // values 1,1,1,2,2,3 -> majority 1, minority 3, variety 3.
        let bands: Vec<Vec<f64>> = [1.0, 1.0, 1.0, 2.0, 2.0, 3.0]
            .iter()
            .map(|&v| vec![v])
            .collect();
        let path = multiband(1, 1, 6, bands);
        assert_eq!(
            run(json!({ "inputs": path.clone(), "statistic": "majority" })).get(0, 0, 0),
            1.0
        );
        assert_eq!(
            run(json!({ "inputs": path.clone(), "statistic": "minority" })).get(0, 0, 0),
            3.0
        );
        assert_eq!(
            run(json!({ "inputs": path, "statistic": "variety" })).get(0, 0, 0),
            3.0
        );
    }

    /// ignore_nodata=true skips nulls; =false makes the whole cell no-data.
    #[test]
    fn nodata_handling() {
        // values 4, nodata, 8 -> mean of 4,8 = 6 when ignoring.
        let bands = vec![vec![4.0], vec![-9999.0], vec![8.0]];
        let path = multiband(1, 1, 3, bands);
        let ignore =
            run(json!({ "inputs": path.clone(), "statistic": "mean", "ignore_nodata": true }));
        assert_eq!(ignore.get(0, 0, 0), 6.0);
        let strict = run(json!({ "inputs": path, "statistic": "mean", "ignore_nodata": false }));
        assert_eq!(strict.get(0, 0, 0), -9999.0, "one nodata -> cell nodata");
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CellStatisticsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err()); // missing inputs
        assert!(bad(json!({ "inputs": "a.tif", "statistic": "bogus" })).is_err());
        assert!(bad(
            json!({ "inputs": "a.tif", "statistic": "percentile", "percentile_value": 150 })
        )
        .is_err());
        assert!(bad(json!({ "inputs": "a.tif", "statistic": "mean" })).is_ok());
    }
}
