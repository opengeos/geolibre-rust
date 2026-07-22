//! GeoLibre tool: per-pixel *argument* statistics over a raster stack.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Find Argument Statistics* (Image
//! Analyst). Where the bundled stack tools return aggregated **values**
//! (`image_stack_profile` samples a series; cell-statistics compute value
//! aggregates), this tool answers *where* a condition occurs along the stack:
//! the slice index (or date) of the maximum or minimum, the position of the
//! median, the number of slices that satisfy a threshold, or the longest
//! consecutive run that does. "When does NDVI peak?", "how many weeks was soil
//! moisture below X?", and "longest dry spell?" are all argument statistics.
//!
//! Inputs are one multiband raster (each band is a time slice) or a
//! comma-separated list of co-registered rasters (every band of every raster is
//! a slice, in order). Statistics:
//! * **argmax** / **argmin** — slice of the maximum / minimum value.
//! * **median_position** — slice whose value is the ordered-middle (median).
//! * **duration** — count of slices satisfying `comparison threshold`.
//! * **longest_run** — longest consecutive run of satisfying slices.
//!
//! For `argmax`/`argmin`/`median_position` the output is the 0-based slice index,
//! or the matching `dates` value (e.g. day-of-year) when `dates` is supplied.
//! For `duration`/`longest_run` the output is a slice count. Pixels with fewer
//! than `min_valid` non-no-data observations become no-data.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::common::{load_input_raster, parse_optional_output, raster_like_with_data};

pub struct FindArgumentStatisticsTool;

impl Tool for FindArgumentStatisticsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "find_argument_statistics",
            display_name: "Find Argument Statistics",
            summary: "Per-pixel argument statistics over a raster stack (like ArcGIS Find Argument Statistics): the slice index (or date) of the max/min, the median position, the count of slices past a threshold (duration), or the longest consecutive run — the argument-position answers the bundled cell-statistics and image_stack_profile value aggregates can't give.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "inputs",
                    description: "One multiband raster (each band is a slice) or a comma-separated list of co-registered rasters, in time order.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output single-band raster of slice indices / dates / counts. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "statistic",
                    description: "'argmax' (default), 'argmin', 'median_position', 'duration', or 'longest_run'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "threshold",
                    description: "Numeric threshold for 'duration'/'longest_run' (required for those).",
                    required: false,
                },
                ToolParamSpec {
                    name: "comparison",
                    description: "Comparison for 'duration'/'longest_run': '>' (default), '>=', '<', or '<='.",
                    required: false,
                },
                ToolParamSpec {
                    name: "dates",
                    description: "Optional comma-separated numeric dates (e.g. day-of-year) per slice; argmax/argmin/median output the matching date instead of the slice index.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_valid",
                    description: "Minimum non-no-data observations per pixel (default 1).",
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
        if let Some(d) = &prm.dates {
            if d.len() != n_slices {
                return Err(ToolError::Validation(format!(
                    "{} dates for {n_slices} slices",
                    d.len()
                )));
            }
        }

        ctx.progress.info(&format!(
            "{} input(s), {n_slices} slice(s), {}x{}, statistic {}",
            rasters.len(),
            rows,
            cols,
            prm.statistic.label()
        ));

        let nodata = -9999.0_f64;
        let mut out = vec![nodata; rows * cols];

        // Reusable per-pixel buffers: (slice_index, value) for valid observations.
        let mut idxs: Vec<usize> = Vec::with_capacity(n_slices);
        let mut vals: Vec<f64> = Vec::with_capacity(n_slices);
        for r in 0..rows {
            for c in 0..cols {
                idxs.clear();
                vals.clear();
                for (si, &(ri, band)) in slices.iter().enumerate() {
                    let ras = &rasters[ri];
                    let v = ras.get(band, r as isize, c as isize);
                    if v != ras.nodata && v.is_finite() {
                        idxs.push(si);
                        vals.push(v);
                    }
                }
                if vals.len() < prm.min_valid {
                    continue;
                }
                let cell = r * cols + c;
                out[cell] = match prm.statistic {
                    Statistic::ArgMax => prm.emit(arg_extreme(&idxs, &vals, true)),
                    Statistic::ArgMin => prm.emit(arg_extreme(&idxs, &vals, false)),
                    Statistic::MedianPosition => prm.emit(median_position(&idxs, &vals)),
                    Statistic::Duration => count_matches(&vals, &prm) as f64,
                    Statistic::LongestRun => longest_run(&idxs, &vals, &prm) as f64,
                };
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let out_r = raster_like_with_data(&rasters[0], out, nodata, DataType::F32)?;
        let out_path = crate::common::write_or_store_output(out_r, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("slices".to_string(), json!(n_slices));
        outputs.insert("statistic".to_string(), json!(prm.statistic.label()));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        Ok(ToolRunResult { outputs })
    }
}

/// Slice index of the max (`want_max`) or min value; first occurrence on ties.
fn arg_extreme(idxs: &[usize], vals: &[f64], want_max: bool) -> usize {
    let mut best = 0usize;
    for i in 1..vals.len() {
        let better = if want_max {
            vals[i] > vals[best]
        } else {
            vals[i] < vals[best]
        };
        if better {
            best = i;
        }
    }
    idxs[best]
}

/// Slice index whose value is the ordered middle (median) of the series. For an
/// even count the lower of the two central values is chosen, so the result is
/// always an actual observed slice.
fn median_position(idxs: &[usize], vals: &[f64]) -> usize {
    let mut order: Vec<usize> = (0..vals.len()).collect();
    order.sort_by(|&a, &b| vals[a].total_cmp(&vals[b]));
    let mid = (vals.len() - 1) / 2;
    idxs[order[mid]]
}

/// Number of valid slices satisfying `comparison threshold`.
fn count_matches(vals: &[f64], prm: &Params) -> usize {
    vals.iter().filter(|&&v| prm.matches(v)).count()
}

/// Longest run of consecutive (by original slice order) satisfying slices.
fn longest_run(idxs: &[usize], vals: &[f64], prm: &Params) -> usize {
    let mut best = 0usize;
    let mut cur = 0usize;
    let mut prev: Option<usize> = None;
    for (k, &v) in vals.iter().enumerate() {
        let consecutive = prev.map(|p| idxs[k] == p + 1).unwrap_or(false);
        if prm.matches(v) {
            cur = if consecutive { cur + 1 } else { 1 };
            best = best.max(cur);
        } else {
            cur = 0;
        }
        prev = Some(idxs[k]);
    }
    best
}

// ── Parameters ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Statistic {
    ArgMax,
    ArgMin,
    MedianPosition,
    Duration,
    LongestRun,
}

impl Statistic {
    fn label(&self) -> &'static str {
        match self {
            Statistic::ArgMax => "argmax",
            Statistic::ArgMin => "argmin",
            Statistic::MedianPosition => "median_position",
            Statistic::Duration => "duration",
            Statistic::LongestRun => "longest_run",
        }
    }
    fn is_threshold_mode(&self) -> bool {
        matches!(self, Statistic::Duration | Statistic::LongestRun)
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Comparison {
    Gt,
    Ge,
    Lt,
    Le,
}

struct Params {
    statistic: Statistic,
    threshold: f64,
    comparison: Comparison,
    dates: Option<Vec<f64>>,
    min_valid: usize,
}

impl Params {
    /// Maps a chosen slice index to the emitted value: the date when `dates` is
    /// supplied, otherwise the slice index itself.
    fn emit(&self, slice: usize) -> f64 {
        match &self.dates {
            Some(d) => d[slice],
            None => slice as f64,
        }
    }
    fn matches(&self, v: f64) -> bool {
        match self.comparison {
            Comparison::Gt => v > self.threshold,
            Comparison::Ge => v >= self.threshold,
            Comparison::Lt => v < self.threshold,
            Comparison::Le => v <= self.threshold,
        }
    }
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
        None | Some("") | Some("argmax") => Statistic::ArgMax,
        Some("argmin") => Statistic::ArgMin,
        Some("median_position") => Statistic::MedianPosition,
        Some("duration") => Statistic::Duration,
        Some("longest_run") => Statistic::LongestRun,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'statistic' must be one of argmax|argmin|median_position|duration|longest_run, got '{o}'"
            )))
        }
    };

    let comparison = match args
        .get("comparison")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some(">") | Some("gt") => Comparison::Gt,
        Some(">=") | Some("ge") => Comparison::Ge,
        Some("<") | Some("lt") => Comparison::Lt,
        Some("<=") | Some("le") => Comparison::Le,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'comparison' must be one of >|>=|<|<=, got '{o}'"
            )))
        }
    };

    // threshold: accept a JSON number or a numeric string.
    let threshold_val = match args.get("threshold") {
        None | Some(Value::Null) => None,
        Some(Value::Number(n)) => Some(n.as_f64().unwrap_or(f64::NAN)),
        Some(Value::String(s)) if s.trim().is_empty() => None,
        Some(Value::String(s)) => Some(
            s.trim()
                .parse::<f64>()
                .map_err(|_| ToolError::Validation("'threshold' must be a number".into()))?,
        ),
        Some(_) => return Err(ToolError::Validation("'threshold' must be a number".into())),
    };
    if statistic.is_threshold_mode() && threshold_val.is_none() {
        return Err(ToolError::Validation(format!(
            "statistic '{}' requires a 'threshold'",
            statistic.label()
        )));
    }
    let threshold = threshold_val.unwrap_or(0.0);

    let dates = match args.get("dates").and_then(Value::as_str) {
        None => None,
        Some(s) if s.trim().is_empty() => None,
        Some(s) => Some(
            s.split(',')
                .map(str::trim)
                .filter(|x| !x.is_empty())
                .map(|x| {
                    x.parse::<f64>()
                        .map_err(|_| ToolError::Validation(format!("date '{x}' is not a number")))
                })
                .collect::<Result<Vec<_>, _>>()?,
        ),
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
        _ => 1,
    };

    Ok(Params {
        statistic,
        threshold,
        comparison,
        dates,
        min_valid,
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
        let out = FindArgumentStatisticsTool.run(&args, &ctx()).unwrap();
        load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    /// argmax over a multiband raster recovers the exact band of the maximum.
    #[test]
    fn argmax_recovers_peak_band() {
        // 2x1 raster, 4 bands. Pixel 0 peaks at band 2; pixel 1 peaks at band 0.
        let bands = vec![
            vec![1.0, 9.0], // band 0
            vec![2.0, 3.0], // band 1
            vec![8.0, 1.0], // band 2
            vec![4.0, 2.0], // band 3
        ];
        let path = multiband(2, 1, 4, bands);
        let r = run(json!({ "inputs": path, "statistic": "argmax" }));
        assert_eq!(r.get(0, 0, 0), 2.0, "pixel 0 max at band 2");
        assert_eq!(r.get(0, 0, 1), 0.0, "pixel 1 max at band 0");
    }

    /// argmin picks the minimum slice; dates map the slice to a date value.
    #[test]
    fn argmin_and_dates() {
        let bands = vec![vec![5.0], vec![1.0], vec![3.0]];
        let path = multiband(1, 1, 3, bands);
        let idx = run(json!({ "inputs": path.clone(), "statistic": "argmin" }));
        assert_eq!(idx.get(0, 0, 0), 1.0, "min at slice 1");
        // day-of-year dates -> output the date at the argmin slice.
        let doy = run(json!({ "inputs": path, "statistic": "argmin", "dates": "10,20,30" }));
        assert_eq!(doy.get(0, 0, 0), 20.0, "date of argmin slice");
    }

    /// duration counts matching slices; longest_run finds the longest streak.
    #[test]
    fn duration_and_longest_run() {
        // series per single pixel: 1,5,6,2,7,8  with threshold > 4  -> matches at
        // slices 1,2,4,5 => duration 4; longest consecutive run = 2 (slices 4,5
        // and 1,2 are both length 2).
        let bands = vec![
            vec![1.0],
            vec![5.0],
            vec![6.0],
            vec![2.0],
            vec![7.0],
            vec![8.0],
        ];
        let path = multiband(1, 1, 6, bands);
        let dur = run(
            json!({ "inputs": path.clone(), "statistic": "duration", "threshold": 4.0, "comparison": ">" }),
        );
        assert_eq!(dur.get(0, 0, 0), 4.0, "4 slices above 4");
        let run_r = run(
            json!({ "inputs": path, "statistic": "longest_run", "threshold": 4.0, "comparison": ">" }),
        );
        assert_eq!(run_r.get(0, 0, 0), 2.0, "longest streak above 4 is 2");
    }

    /// median_position returns the slice of the ordered-middle value.
    #[test]
    fn median_position_is_middle_slice() {
        // values 10,50,30,40,20 -> sorted middle is 30 at original slice 2.
        let bands = vec![vec![10.0], vec![50.0], vec![30.0], vec![40.0], vec![20.0]];
        let path = multiband(1, 1, 5, bands);
        let r = run(json!({ "inputs": path, "statistic": "median_position" }));
        assert_eq!(r.get(0, 0, 0), 2.0, "median value 30 is at slice 2");
    }

    /// Pixels with too few valid observations become no-data.
    #[test]
    fn insufficient_valid_is_nodata() {
        let bands = vec![vec![-9999.0], vec![-9999.0]];
        let path = multiband(1, 1, 2, bands);
        let r = run(json!({ "inputs": path, "statistic": "argmax", "min_valid": 1 }));
        assert_eq!(r.get(0, 0, 0), -9999.0, "no valid observations -> nodata");
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            FindArgumentStatisticsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err()); // missing inputs
        assert!(bad(json!({ "inputs": "a.tif", "statistic": "bogus" })).is_err());
        assert!(bad(json!({ "inputs": "a.tif", "statistic": "duration" })).is_err()); // no threshold
        assert!(bad(json!({ "inputs": "a.tif", "comparison": "!=" })).is_err());
        assert!(bad(json!({ "inputs": "a.tif", "statistic": "argmax" })).is_ok());
        assert!(
            bad(json!({ "inputs": "a.tif", "statistic": "duration", "threshold": 4.0 })).is_ok()
        );
    }
}
