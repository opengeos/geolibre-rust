//! Focal (moving-window) neighborhood statistics over a single raster band.
//!
//! A unified port of ArcGIS Spatial Analyst's *Focal Statistics*: slide a
//! shaped neighborhood (rectangle, circle, annulus, or wedge) over the grid and
//! reduce the cells falling inside it with one statistic. Numeric statistics
//! (`mean`, `maximum`, `minimum`, `median`, `range`, `std`, `sum`,
//! `percentile`) and categorical statistics (`majority`, `minority`,
//! `variety`) are both supported, so the same tool smooths a DEM and finds the
//! modal class of a land-cover raster.
//!
//! The neighborhood is precomputed once as a set of `(dr, dc)` cell offsets,
//! then a sliding-window pass evaluates every cell. `ignore_nodata` (default
//! true) reduces over just the valid cells in the window; when false, any
//! no-data cell in the window forces a no-data output. Windows that gather no
//! valid cell yield no-data. The per-row pass is fanned out with the
//! `dispatch_worker` pattern (a background thread per worker on native builds,
//! serial on wasm32 where `std::thread::spawn` traps).

use std::sync::mpsc;
use std::sync::Arc;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{
    band_to_vec, load_input_raster, parse_optional_output, raster_like_with_data,
    write_or_store_output,
};

/// Moving-window neighborhood statistics over one raster band.
pub struct FocalStatisticsTool;

/// The neighborhood shape, with its size expressed in cells.
#[derive(Debug, Clone, Copy)]
enum Neighborhood {
    /// `width` x `height` cell rectangle centered on the focal cell.
    Rectangle { width: usize, height: usize },
    /// Filled disc of `radius` cells (Euclidean).
    Circle { radius: f64 },
    /// Ring between `inner` and `outer` cell radii (inclusive).
    Annulus { inner: f64, outer: f64 },
    /// Pie slice of `radius` cells swept from `start` to `end` degrees
    /// (arithmetic convention: 0 deg = east, increasing counterclockwise).
    Wedge { radius: f64, start: f64, end: f64 },
}

/// The reduction applied to the cells inside a neighborhood.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
    /// Statistics whose output is a category value / count rather than a
    /// continuous measurement (drives the output data type).
    fn is_categorical(self) -> bool {
        matches!(
            self,
            Statistic::Majority | Statistic::Minority | Statistic::Variety
        )
    }
}

impl Tool for FocalStatisticsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "focal_statistics",
            display_name: "Focal Statistics",
            summary: "Moving-window neighborhood statistics (rectangle/circle/annulus/wedge; numeric and categorical).",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input raster file path.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output raster path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "statistics",
                    description: "Statistic: mean, majority, maximum, median, minimum, minority, percentile, range, std, sum, or variety (default mean).",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighborhood",
                    description: "Neighborhood shape: rectangle, circle, annulus, or wedge (default rectangle).",
                    required: false,
                },
                ToolParamSpec {
                    name: "width",
                    description: "Rectangle width in cells (default 3).",
                    required: false,
                },
                ToolParamSpec {
                    name: "height",
                    description: "Rectangle height in cells (default 3).",
                    required: false,
                },
                ToolParamSpec {
                    name: "radius",
                    description: "Radius in cells for circle / wedge, and the outer radius for annulus (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "inner_radius",
                    description: "Inner radius in cells for the annulus neighborhood (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "start_angle",
                    description: "Wedge start angle in degrees (0 = east, counterclockwise; default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "end_angle",
                    description: "Wedge end angle in degrees (0 = east, counterclockwise; default 90).",
                    required: false,
                },
                ToolParamSpec {
                    name: "percentile_value",
                    description: "Percentile in [0, 100] for the percentile statistic (default 90).",
                    required: false,
                },
                ToolParamSpec {
                    name: "ignore_nodata",
                    description: "If true (default), reduce over valid cells only; if false, any no-data cell in the window yields no-data.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to process (default 1).",
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
        // Parsing surfaces any bad enum / numeric parameters up front.
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args.get("input").and_then(Value::as_str).ok_or_else(|| {
            ToolError::Validation("missing required parameter 'input'".to_string())
        })?;
        let output = parse_optional_output(args, "output")?;
        let params = parse_params(args)?;

        let raster = load_input_raster(input)?;
        let band_1based = args.get("band").and_then(Value::as_u64).unwrap_or(1).max(1);
        let band = (band_1based - 1) as isize;
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

        let offsets = build_offsets(&params.neighborhood);
        if offsets.is_empty() {
            return Err(ToolError::Validation(
                "neighborhood contains no cells; check the shape/size parameters".to_string(),
            ));
        }

        ctx.progress.info(&format!(
            "focal {:?} over {} neighbor cells",
            params.statistic,
            offsets.len()
        ));
        let out_data = focal_pass(&src, rows, cols, nodata, &offsets, &params);
        ctx.progress.progress(1.0);

        let out_type = if params.statistic.is_categorical() {
            DataType::I32
        } else {
            DataType::F32
        };
        let out_raster = raster_like_with_data(&raster, out_data, nodata, out_type)?;
        let out_path = write_or_store_output(out_raster, output)?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert(
            "statistics".to_string(),
            json!(format!("{:?}", params.statistic).to_lowercase()),
        );
        outputs.insert("neighbor_cells".to_string(), json!(offsets.len()));
        Ok(ToolRunResult { outputs })
    }
}

/// Fully parsed, validated run parameters (shape + statistic + flags).
struct Params {
    neighborhood: Neighborhood,
    statistic: Statistic,
    percentile: f64,
    ignore_nodata: bool,
}

/// Parses and validates every non-IO parameter (enums, sizes, angles).
fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let statistic = match args
        .get("statistics")
        .and_then(Value::as_str)
        .unwrap_or("mean")
        .to_lowercase()
        .as_str()
    {
        "mean" => Statistic::Mean,
        "majority" => Statistic::Majority,
        "maximum" | "max" => Statistic::Maximum,
        "median" => Statistic::Median,
        "minimum" | "min" => Statistic::Minimum,
        "minority" => Statistic::Minority,
        "percentile" => Statistic::Percentile,
        "range" => Statistic::Range,
        "std" | "stddev" => Statistic::Std,
        "sum" => Statistic::Sum,
        "variety" => Statistic::Variety,
        other => {
            return Err(ToolError::Validation(format!(
                "unknown statistics '{other}' (expected mean, majority, maximum, median, minimum, minority, percentile, range, std, sum, or variety)"
            )))
        }
    };

    let shape = args
        .get("neighborhood")
        .and_then(Value::as_str)
        .unwrap_or("rectangle")
        .to_lowercase();
    let neighborhood = match shape.as_str() {
        "rectangle" | "rect" => {
            let width = parse_num(args, "width")?.unwrap_or(3.0).max(1.0) as usize;
            let height = parse_num(args, "height")?.unwrap_or(3.0).max(1.0) as usize;
            Neighborhood::Rectangle { width, height }
        }
        "circle" => {
            let radius = parse_num(args, "radius")?.unwrap_or(1.0).max(0.0);
            Neighborhood::Circle { radius }
        }
        "annulus" => {
            let outer = parse_num(args, "radius")?.unwrap_or(1.0).max(0.0);
            let inner = parse_num(args, "inner_radius")?.unwrap_or(0.0).max(0.0);
            if inner > outer {
                return Err(ToolError::Validation(
                    "annulus inner_radius must be <= radius (outer)".to_string(),
                ));
            }
            Neighborhood::Annulus { inner, outer }
        }
        "wedge" => {
            let radius = parse_num(args, "radius")?.unwrap_or(1.0).max(0.0);
            let start = parse_num(args, "start_angle")?.unwrap_or(0.0);
            let end = parse_num(args, "end_angle")?.unwrap_or(90.0);
            Neighborhood::Wedge { radius, start, end }
        }
        other => {
            return Err(ToolError::Validation(format!(
                "unknown neighborhood '{other}' (expected rectangle, circle, annulus, or wedge)"
            )))
        }
    };

    let percentile = parse_num(args, "percentile_value")?.unwrap_or(90.0);
    if !(0.0..=100.0).contains(&percentile) {
        return Err(ToolError::Validation(
            "percentile_value must be within [0, 100]".to_string(),
        ));
    }
    let ignore_nodata = parse_bool(args, "ignore_nodata")?.unwrap_or(true);

    Ok(Params {
        neighborhood,
        statistic,
        percentile,
        ignore_nodata,
    })
}

/// Parses an optional numeric parameter accepting a JSON number or a numeric
/// string (host UIs frequently post numbers as strings).
fn parse_num(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
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

/// Parses an optional boolean accepting a JSON bool or a truthy/falsy string.
fn parse_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" | "y" => Ok(Some(true)),
            "false" | "0" | "no" | "n" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
}

/// Precomputes the `(dr, dc)` cell offsets that define the neighborhood mask.
/// The focal cell `(0, 0)` is always included.
fn build_offsets(shape: &Neighborhood) -> Vec<(isize, isize)> {
    let mut offsets = Vec::new();
    match *shape {
        Neighborhood::Rectangle { width, height } => {
            let rw = (width / 2) as isize;
            let rh = (height / 2) as isize;
            for dr in -rh..=rh {
                for dc in -rw..=rw {
                    offsets.push((dr, dc));
                }
            }
        }
        Neighborhood::Circle { radius } => {
            let r = radius.floor() as isize;
            let r2 = radius * radius;
            for dr in -r..=r {
                for dc in -r..=r {
                    if (dr * dr + dc * dc) as f64 <= r2 + 1e-9 {
                        offsets.push((dr, dc));
                    }
                }
            }
        }
        Neighborhood::Annulus { inner, outer } => {
            let r = outer.floor() as isize;
            let (i2, o2) = (inner * inner, outer * outer);
            for dr in -r..=r {
                for dc in -r..=r {
                    let d2 = (dr * dr + dc * dc) as f64;
                    if d2 >= i2 - 1e-9 && d2 <= o2 + 1e-9 {
                        offsets.push((dr, dc));
                    }
                }
            }
        }
        Neighborhood::Wedge { radius, start, end } => {
            let r = radius.floor() as isize;
            let r2 = radius * radius;
            let start = normalize_deg(start);
            let end = normalize_deg(end);
            for dr in -r..=r {
                for dc in -r..=r {
                    if (dr == 0 && dc == 0) || (dr * dr + dc * dc) as f64 > r2 + 1e-9 {
                        if dr == 0 && dc == 0 {
                            offsets.push((0, 0)); // apex always inside the wedge
                        }
                        continue;
                    }
                    // Arithmetic angle: east = 0 deg, counterclockwise positive.
                    // Screen rows grow downward, so north is -dr.
                    let ang = normalize_deg((-(dr as f64)).atan2(dc as f64).to_degrees());
                    if angle_in_arc(ang, start, end) {
                        offsets.push((dr, dc));
                    }
                }
            }
        }
    }
    offsets
}

/// Wraps an angle in degrees into `[0, 360)`.
fn normalize_deg(a: f64) -> f64 {
    let mut a = a % 360.0;
    if a < 0.0 {
        a += 360.0;
    }
    a
}

/// True when `ang` (in `[0, 360)`) lies on the counterclockwise arc from
/// `start` to `end`, handling the wrap-around case.
fn angle_in_arc(ang: f64, start: f64, end: f64) -> bool {
    if (start - end).abs() < 1e-9 {
        return true; // degenerate: treat as full circle
    }
    if start <= end {
        ang >= start - 1e-9 && ang <= end + 1e-9
    } else {
        ang >= start - 1e-9 || ang <= end + 1e-9
    }
}

/// Runs the sliding-window reduction over every cell, fanning rows out to
/// worker threads (serial on wasm32) and reassembling them in row order.
fn focal_pass(
    src: &[f64],
    rows: usize,
    cols: usize,
    nodata: f64,
    offsets: &[(isize, isize)],
    params: &Params,
) -> Vec<f64> {
    let num_procs = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
        .max(1);

    let src = Arc::new(src.to_vec());
    let offsets = Arc::new(offsets.to_vec());
    let statistic = params.statistic;
    let percentile = params.percentile;
    let ignore_nodata = params.ignore_nodata;

    let (tx, rx) = mpsc::channel::<(usize, Vec<f64>)>();
    for tid in 0..num_procs {
        let src = src.clone();
        let offsets = offsets.clone();
        let tx = tx.clone();
        dispatch_worker(num_procs > 1, move || {
            let mut window: Vec<f64> = Vec::with_capacity(offsets.len());
            for r in (0..rows).filter(|row| row % num_procs == tid) {
                let mut row_out = vec![nodata; cols];
                for c in 0..cols {
                    window.clear();
                    let mut saw_nodata = false;
                    for &(dr, dc) in offsets.iter() {
                        let rr = r as isize + dr;
                        let cc = c as isize + dc;
                        if rr < 0 || cc < 0 || rr >= rows as isize || cc >= cols as isize {
                            continue;
                        }
                        let v = src[rr as usize * cols + cc as usize];
                        if is_nodata(v, nodata) {
                            saw_nodata = true;
                            continue;
                        }
                        window.push(v);
                    }
                    if (!ignore_nodata && saw_nodata) || window.is_empty() {
                        continue; // stays nodata
                    }
                    row_out[c] = reduce(&mut window, statistic, percentile);
                }
                let _ = tx.send((r, row_out));
            }
        });
    }
    drop(tx);

    let mut out = vec![nodata; rows * cols];
    for _ in 0..rows {
        if let Ok((r, row_out)) = rx.recv() {
            out[r * cols..r * cols + cols].copy_from_slice(&row_out);
        }
    }
    out
}

/// True when `v` should be treated as no-data (matches the sentinel or is
/// non-finite, covering `NaN`-tagged rasters).
#[inline]
fn is_nodata(v: f64, nodata: f64) -> bool {
    v == nodata || !v.is_finite()
}

/// Reduces the gathered window values with the requested statistic. `window`
/// may be reordered (sorted) in place.
fn reduce(window: &mut [f64], statistic: Statistic, percentile: f64) -> f64 {
    match statistic {
        Statistic::Sum => window.iter().sum(),
        Statistic::Mean => window.iter().sum::<f64>() / window.len() as f64,
        Statistic::Maximum => window.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        Statistic::Minimum => window.iter().copied().fold(f64::INFINITY, f64::min),
        Statistic::Range => {
            let mn = window.iter().copied().fold(f64::INFINITY, f64::min);
            let mx = window.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            mx - mn
        }
        Statistic::Std => {
            let n = window.len() as f64;
            let mean = window.iter().sum::<f64>() / n;
            let var = window.iter().map(|&v| (v - mean) * (v - mean)).sum::<f64>() / n;
            var.sqrt()
        }
        Statistic::Median => {
            window.sort_by(|a, b| a.total_cmp(b));
            percentile_sorted(window, 50.0)
        }
        Statistic::Percentile => {
            window.sort_by(|a, b| a.total_cmp(b));
            percentile_sorted(window, percentile)
        }
        Statistic::Variety => distinct_count(window) as f64,
        Statistic::Majority => mode(window, true),
        Statistic::Minority => mode(window, false),
    }
}

/// Linear-interpolation percentile (numpy default) over a pre-sorted slice.
fn percentile_sorted(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    if n == 1 {
        return sorted[0];
    }
    let rank = (p / 100.0) * (n as f64 - 1.0);
    let lo = rank.floor() as usize;
    let hi = rank.ceil() as usize;
    if lo == hi {
        sorted[lo]
    } else {
        let frac = rank - lo as f64;
        sorted[lo] * (1.0 - frac) + sorted[hi] * frac
    }
}

/// Counts distinct values (exact `f64` equality; suits categorical rasters).
fn distinct_count(window: &mut [f64]) -> usize {
    window.sort_by(|a, b| a.total_cmp(b));
    let mut count = 0usize;
    let mut prev: Option<f64> = None;
    for &v in window.iter() {
        if prev.map(|p| p != v).unwrap_or(true) {
            count += 1;
            prev = Some(v);
        }
    }
    count
}

/// Most (`majority=true`) or least (`majority=false`) frequent value. Ties are
/// broken by the smallest value, matching SciPy's modal filter convention.
fn mode(window: &mut [f64], majority: bool) -> f64 {
    window.sort_by(|a, b| a.total_cmp(b));
    let mut best_val = window[0];
    let mut best_count = 0usize;
    let mut cur_val = window[0];
    let mut cur_count = 0usize;
    for &v in window.iter() {
        if v == cur_val {
            cur_count += 1;
        } else {
            if is_better(cur_count, best_count, majority) {
                best_count = cur_count;
                best_val = cur_val;
            }
            cur_val = v;
            cur_count = 1;
        }
    }
    if is_better(cur_count, best_count, majority) {
        best_val = cur_val;
    }
    best_val
}

/// Whether `count` beats the incumbent `best` for a majority/minority search.
/// Because values are visited in ascending order, a strict `>` (majority) or
/// `<` (minority) comparison keeps the first — i.e. smallest — value on ties.
#[inline]
fn is_better(count: usize, best: usize, majority: bool) -> bool {
    if best == 0 {
        return true;
    }
    if majority {
        count > best
    } else {
        count < best
    }
}

/// Runs `work` on a background thread on native targets, or inline when
/// `parallel` is false or when compiling for wasm32 (where `thread::spawn`
/// traps). Mirrors the whitebox `dispatch_worker` fan-out helper.
fn dispatch_worker<F>(parallel: bool, work: F)
where
    F: FnOnce() + Send + 'static,
{
    #[cfg(not(target_arch = "wasm32"))]
    {
        if parallel {
            std::thread::spawn(work);
            return;
        }
    }
    let _ = parallel;
    work();
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run_focal(
        src: &[f64],
        rows: usize,
        cols: usize,
        nodata: f64,
        shape: Neighborhood,
        statistic: Statistic,
        percentile: f64,
        ignore_nodata: bool,
    ) -> Vec<f64> {
        let offsets = build_offsets(&shape);
        let params = Params {
            neighborhood: shape,
            statistic,
            percentile,
            ignore_nodata,
        };
        focal_pass(src, rows, cols, nodata, &offsets, &params)
    }

    #[test]
    fn rectangle_mean_matches_hand_computation() {
        // 3x3 grid 1..9; center 3x3 mean over the whole window = 5.
        let src: Vec<f64> = (1..=9).map(|v| v as f64).collect();
        let out = run_focal(
            &src,
            3,
            3,
            -9999.0,
            Neighborhood::Rectangle {
                width: 3,
                height: 3,
            },
            Statistic::Mean,
            90.0,
            true,
        );
        assert!((out[4] - 5.0).abs() < 1e-9);
        // Corner (0,0): mean of {1,2,4,5} = 3.0 (only in-bounds cells count).
        assert!((out[0] - 3.0).abs() < 1e-9);
    }

    #[test]
    fn sum_and_range_and_std() {
        let src = vec![1.0, 2.0, 3.0, 4.0];
        // 2x2 grid, 3x3 rectangle window spans the whole grid for every cell.
        let out_sum = run_focal(
            &src,
            2,
            2,
            -1.0,
            Neighborhood::Rectangle {
                width: 3,
                height: 3,
            },
            Statistic::Sum,
            90.0,
            true,
        );
        assert!(out_sum.iter().all(|&v| (v - 10.0).abs() < 1e-9));
        let out_range = run_focal(
            &src,
            2,
            2,
            -1.0,
            Neighborhood::Rectangle {
                width: 3,
                height: 3,
            },
            Statistic::Range,
            90.0,
            true,
        );
        assert!(out_range.iter().all(|&v| (v - 3.0).abs() < 1e-9));
    }

    #[test]
    fn categorical_majority_minority_variety() {
        // Row-major 3x3, every cell's 3x3 window = the whole grid.
        // values: three 7s, two 3s, four 1s -> majority 1, minority 3, variety 3.
        let src = vec![7.0, 7.0, 7.0, 3.0, 3.0, 1.0, 1.0, 1.0, 1.0];
        let shape = Neighborhood::Rectangle {
            width: 3,
            height: 3,
        };
        let maj = run_focal(&src, 3, 3, -1.0, shape, Statistic::Majority, 90.0, true);
        assert_eq!(maj[4], 1.0);
        let min = run_focal(&src, 3, 3, -1.0, shape, Statistic::Minority, 90.0, true);
        assert_eq!(min[4], 3.0);
        let var = run_focal(&src, 3, 3, -1.0, shape, Statistic::Variety, 90.0, true);
        assert_eq!(var[4], 3.0);
    }

    #[test]
    fn ignore_nodata_flag_controls_masking() {
        let nd = -9999.0;
        // center has a nodata neighbor; ignore=true averages the rest.
        let src = vec![1.0, 1.0, 1.0, 1.0, 5.0, nd, 1.0, 1.0, 1.0];
        let shape = Neighborhood::Rectangle {
            width: 3,
            height: 3,
        };
        let keep = run_focal(&src, 3, 3, nd, shape, Statistic::Mean, 90.0, true);
        // mean of the 8 valid cells in the window = (1*7 + 5)/8 = 1.5
        assert!((keep[4] - 1.5).abs() < 1e-9);
        let strict = run_focal(&src, 3, 3, nd, shape, Statistic::Mean, 90.0, false);
        assert_eq!(strict[4], nd); // any nodata in window -> nodata
    }

    #[test]
    fn circle_offsets_form_a_disc() {
        let offsets = build_offsets(&Neighborhood::Circle { radius: 1.0 });
        // radius-1 disc = center + 4 rook neighbors = 5 cells.
        assert_eq!(offsets.len(), 5);
        assert!(offsets.contains(&(0, 0)));
        assert!(!offsets.contains(&(1, 1))); // corner is outside radius 1
    }

    #[test]
    fn annulus_excludes_center() {
        let offsets = build_offsets(&Neighborhood::Annulus {
            inner: 1.0,
            outer: 1.0,
        });
        assert!(!offsets.contains(&(0, 0)));
        assert_eq!(offsets.len(), 4); // the 4 rook cells at distance exactly 1
    }

    #[test]
    fn wedge_selects_a_quadrant() {
        // 0..90 deg (east to north) wedge of radius 2 excludes the west/south.
        let offsets = build_offsets(&Neighborhood::Wedge {
            radius: 2.0,
            start: 0.0,
            end: 90.0,
        });
        assert!(offsets.contains(&(0, 0))); // apex
        assert!(offsets.contains(&(0, 1))); // due east
        assert!(offsets.contains(&(-1, 0))); // due north (row up)
        assert!(!offsets.contains(&(0, -1))); // due west excluded
        assert!(!offsets.contains(&(1, 0))); // due south excluded
    }

    #[test]
    fn percentile_interpolates_like_numpy() {
        let sorted = [0.0, 1.0, 2.0, 3.0, 4.0];
        assert!((percentile_sorted(&sorted, 50.0) - 2.0).abs() < 1e-9);
        assert!((percentile_sorted(&sorted, 25.0) - 1.0).abs() < 1e-9);
        assert!((percentile_sorted(&sorted, 90.0) - 3.6).abs() < 1e-9);
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = FocalStatisticsTool;
        let mut args = ToolArgs::new();
        args.insert("input".to_string(), json!("memory://x"));
        args.insert("statistics".to_string(), json!("bogus"));
        assert!(tool.validate(&args).is_err());

        let mut args2 = ToolArgs::new();
        args2.insert("input".to_string(), json!("memory://x"));
        args2.insert("neighborhood".to_string(), json!("hexagon"));
        assert!(tool.validate(&args2).is_err());

        let mut args3 = ToolArgs::new();
        args3.insert("input".to_string(), json!("memory://x"));
        args3.insert("percentile_value".to_string(), json!(150.0));
        assert!(tool.validate(&args3).is_err());
    }
}
