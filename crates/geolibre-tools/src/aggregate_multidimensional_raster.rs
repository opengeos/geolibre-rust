//! GeoLibre tool: reduce a raster cube along its dimension into bins.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Aggregate Multidimensional Raster*
//! (Image Analyst / Spatial Analyst).
//!
//! ## Why the catalog needs it
//!
//! Turning a dense time series into composites is the first step of nearly
//! every climate, phenology or monitoring workflow: 365 daily NDVI slices
//! become 12 monthly means, or a decade of scenes becomes one value per year.
//!
//! Neither registry can do it. `cell_statistics` (round 17) reduces a stack but
//! only ever to **one** layer — it has no concept of a dimension, so it cannot
//! produce a per-bin output. `multidimensional_anomaly` (round 16) is the
//! opposite operation: it rewrites every slice against a baseline and preserves
//! the slice count. `time_series_forecast` builds an H3 space-time cube from
//! *vector* points, not a raster stack.
//!
//! ## Binning
//!
//! Slices are grouped by their dimension coordinate:
//!
//! * `all` — one bin holding every slice (the `cell_statistics` case).
//! * `interval_value` — fixed-width bins of `interval_value` coordinate units.
//! * `interval_count` — that many equal-width bins spanning the range.
//! * `interval_ranges` — explicit `"start:end,start:end"` bins, half-open on
//!   the upper edge except the last, which includes it.
//!
//! With no `dimension_values` the coordinate is the 1-based slice index, so
//! `interval_value: 3` means "every three slices".
//!
//! **Scope:** ArcGIS also offers calendar keywords (`MONTHLY`,
//! `RECURRING_MONTHLY`, …). Those need real dates with calendar arithmetic, and
//! this catalog's cubes carry plain numeric coordinates, so they are not
//! offered rather than being approximated with 30.44-day bins that would
//! silently disagree with a calendar month. `interval_ranges` expresses any of
//! them exactly when the caller knows the dates.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::args_common::{bool_or, choice_or, f64_or, opt_f64, opt_usize};
use crate::common::{parse_optional_output, write_or_store_output};
use crate::cube::{load_cube, Cube};
use crate::raster_stack::raster_like_multiband;

pub struct AggregateMultidimensionalRasterTool;

impl Tool for AggregateMultidimensionalRasterTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "aggregate_multidimensional_raster",
            display_name: "Aggregate Multidimensional Raster",
            summary: "Bins the slices of a raster cube along its dimension and reduces each bin to one output band — daily scenes into monthly composites, a decade into annual means (ArcGIS Aggregate Multidimensional Raster). Round 17's cell_statistics reduces a stack only to a single layer and has no notion of a dimension, while round 16's multidimensional_anomaly performs the opposite operation and preserves the slice count.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "One multiband raster (each band is a slice) or a comma-separated list of co-registered rasters.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output multiband raster, one band per bin. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "dimension",
                    description: "Name of the dimension, used only in the report (default 'slice').",
                    required: false,
                },
                ToolParamSpec {
                    name: "dimension_values",
                    description: "Comma-separated coordinate of each slice, strictly increasing. Defaults to the 1-based slice index.",
                    required: false,
                },
                ToolParamSpec {
                    name: "aggregation_definition",
                    description: "'all' (default), 'interval_value', 'interval_count', or 'interval_ranges'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "interval_value",
                    description: "Bin width in dimension units, for 'interval_value'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "interval_count",
                    description: "Number of equal-width bins, for 'interval_count'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "interval_ranges",
                    description: "Explicit bins as 'start:end,start:end', for 'interval_ranges'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "aggregation_method",
                    description: "'mean' (default), 'majority', 'maximum', 'median', 'minimum', 'minority', 'percentile', 'range', 'std', 'sum', or 'variety'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "percentile_value",
                    description: "Percentile in [0, 100] for method 'percentile' (default 50).",
                    required: false,
                },
                ToolParamSpec {
                    name: "ignore_nodata",
                    description: "Skip no-data slices per cell (default true). When false, one no-data slice makes the whole bin no-data at that cell.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        crate::raster_stack::parse_input_paths(args, "input")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let prm = parse_params(args)?;
        let output = parse_optional_output(args, "output")?;
        let cube = load_cube(args, "input", Some("dimension_values"), Some("dimension"), 1)?;

        let bins = build_bins(&cube, &prm)?;
        if bins.is_empty() {
            return Err(ToolError::Execution(
                "the binning produced no non-empty bins".to_string(),
            ));
        }
        ctx.progress.info(&format!(
            "{} slice(s) over '{}' -> {} bin(s), method {}",
            cube.len(),
            cube.dimension,
            bins.len(),
            prm.method.label()
        ));

        let (rows, cols) = (cube.rows, cube.cols);
        let nodata = -9999.0_f64;
        let mut bands: Vec<Vec<f64>> = vec![vec![nodata; rows * cols]; bins.len()];
        let mut vals: Vec<f64> = Vec::new();

        for (bi, bin) in bins.iter().enumerate() {
            for r in 0..rows {
                for c in 0..cols {
                    vals.clear();
                    let mut had_nodata = false;
                    for &s in &bin.slices {
                        match cube.get(s, r, c) {
                            Some(v) => vals.push(v),
                            None => had_nodata = true,
                        }
                    }
                    if vals.is_empty() || (!prm.ignore_nodata && had_nodata) {
                        continue;
                    }
                    bands[bi][r * cols + c] = reduce(&mut vals, &prm);
                }
            }
            ctx.progress
                .progress((bi as f64 + 1.0) / bins.len() as f64);
        }

        let out_raster = raster_like_multiband(cube.template(), &bands, nodata, DataType::F32)?;
        let out_path = write_or_store_output(out_raster, output)?;

        let bin_report: Vec<Value> = bins
            .iter()
            .map(|b| {
                json!({
                    "start": b.start,
                    "end": b.end,
                    "slice_count": b.slices.len(),
                })
            })
            .collect();

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("dimension".to_string(), json!(cube.dimension));
        outputs.insert("input_slices".to_string(), json!(cube.len()));
        outputs.insert("bin_count".to_string(), json!(bins.len()));
        outputs.insert("bins".to_string(), Value::Array(bin_report));
        outputs.insert("aggregation_method".to_string(), json!(prm.method.label()));
        Ok(ToolRunResult { outputs })
    }
}

/// One output bin: its coordinate span and the slices falling inside it.
struct Bin {
    start: f64,
    end: f64,
    slices: Vec<usize>,
}

/// Groups the cube's slices into bins. Empty bins are dropped, so an output
/// band always has data behind it.
fn build_bins(cube: &Cube, prm: &Params) -> Result<Vec<Bin>, ToolError> {
    let n = cube.len();
    let coords: Vec<f64> = (0..n).map(|s| cube.coord(s)).collect();
    let (lo, hi) = (coords[0], coords[n - 1]);

    let mut edges: Vec<(f64, f64)> = Vec::new();
    match prm.definition {
        Definition::All => edges.push((lo, hi)),
        Definition::IntervalValue => {
            let w = prm.interval_value.expect("validated");
            let mut start = lo;
            // `<=` rather than `<` so the final partial bin is still emitted
            // when the span is not a whole multiple of the width.
            while start <= hi {
                edges.push((start, start + w));
                start += w;
            }
        }
        Definition::IntervalCount => {
            let k = prm.interval_count.expect("validated");
            // A degenerate span (one slice, or all coordinates equal) has no
            // width to divide, so it collapses to a single bin.
            let w = if hi > lo {
                (hi - lo) / k as f64
            } else {
                f64::INFINITY
            };
            if !w.is_finite() {
                edges.push((lo, hi));
            } else {
                for i in 0..k {
                    edges.push((lo + i as f64 * w, lo + (i + 1) as f64 * w));
                }
            }
        }
        Definition::IntervalRanges => {
            edges = prm.ranges.clone().expect("validated");
        }
    }

    let last = edges.len().saturating_sub(1);
    let mut bins = Vec::new();
    for (i, &(start, end)) in edges.iter().enumerate() {
        // Half-open [start, end) so a slice on a shared edge lands in exactly
        // one bin; the last bin closes so the final coordinate is not lost.
        let slices: Vec<usize> = (0..n)
            .filter(|&s| {
                let c = coords[s];
                c >= start && (c < end || (i == last && c <= end))
            })
            .collect();
        if !slices.is_empty() {
            bins.push(Bin { start, end, slices });
        }
    }
    Ok(bins)
}

/// Reduces a bin's values at one cell. `vals` may be reordered.
fn reduce(vals: &mut [f64], prm: &Params) -> f64 {
    match prm.method {
        Method::Mean => vals.iter().sum::<f64>() / vals.len() as f64,
        Method::Sum => vals.iter().sum(),
        Method::Maximum => vals.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        Method::Minimum => vals.iter().copied().fold(f64::INFINITY, f64::min),
        Method::Range => {
            let mn = vals.iter().copied().fold(f64::INFINITY, f64::min);
            let mx = vals.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            mx - mn
        }
        Method::Std => {
            let n = vals.len() as f64;
            let mean = vals.iter().sum::<f64>() / n;
            (vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n).sqrt()
        }
        Method::Median => percentile(vals, 50.0),
        Method::Percentile => percentile(vals, prm.percentile_value),
        Method::Variety => {
            let mut seen: Vec<u64> = vals.iter().map(|v| v.to_bits()).collect();
            seen.sort_unstable();
            seen.dedup();
            seen.len() as f64
        }
        Method::Majority => extreme_frequency(vals, true),
        Method::Minority => extreme_frequency(vals, false),
    }
}

/// Linear-interpolated percentile, matching `cell_statistics`.
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

/// Most or least frequent value; ties resolve to the smaller value, matching
/// `cell_statistics`.
fn extreme_frequency(vals: &[f64], want_max: bool) -> f64 {
    let mut counts: BTreeMap<u64, (usize, f64)> = BTreeMap::new();
    for &v in vals {
        let e = counts.entry(v.to_bits()).or_insert((0, v));
        e.0 += 1;
    }
    let mut ordered: Vec<(f64, usize)> = counts.values().map(|(c, v)| (*v, *c)).collect();
    ordered.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut best_count = if want_max { 0usize } else { usize::MAX };
    let mut best_val = f64::NAN;
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
enum Definition {
    All,
    IntervalValue,
    IntervalCount,
    IntervalRanges,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
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

impl Method {
    fn label(self) -> &'static str {
        match self {
            Method::Mean => "mean",
            Method::Majority => "majority",
            Method::Maximum => "maximum",
            Method::Median => "median",
            Method::Minimum => "minimum",
            Method::Minority => "minority",
            Method::Percentile => "percentile",
            Method::Range => "range",
            Method::Std => "std",
            Method::Sum => "sum",
            Method::Variety => "variety",
        }
    }
}

struct Params {
    definition: Definition,
    interval_value: Option<f64>,
    interval_count: Option<usize>,
    ranges: Option<Vec<(f64, f64)>>,
    method: Method,
    percentile_value: f64,
    ignore_nodata: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let definition = match choice_or(
        args,
        "aggregation_definition",
        &["all", "interval_value", "interval_count", "interval_ranges"],
        "all",
    )? {
        "interval_value" => Definition::IntervalValue,
        "interval_count" => Definition::IntervalCount,
        "interval_ranges" => Definition::IntervalRanges,
        _ => Definition::All,
    };

    let interval_value = opt_f64(args, "interval_value")?;
    let interval_count = opt_usize(args, "interval_count")?;
    let ranges = parse_ranges(args)?;

    // Each definition needs its own parameter; falling back to `all` when it is
    // missing would silently collapse the cube to one band.
    match definition {
        Definition::IntervalValue => match interval_value {
            Some(v) if v > 0.0 && v.is_finite() => {}
            _ => {
                return Err(ToolError::Validation(
                    "'interval_value' must be a positive number when aggregation_definition is \
                     'interval_value'"
                        .to_string(),
                ))
            }
        },
        Definition::IntervalCount => match interval_count {
            Some(k) if k >= 1 => {}
            _ => {
                return Err(ToolError::Validation(
                    "'interval_count' must be at least 1 when aggregation_definition is \
                     'interval_count'"
                        .to_string(),
                ))
            }
        },
        Definition::IntervalRanges => {
            if ranges.is_none() {
                return Err(ToolError::Validation(
                    "'interval_ranges' is required when aggregation_definition is \
                     'interval_ranges'"
                        .to_string(),
                ));
            }
        }
        Definition::All => {}
    }

    let method = match choice_or(
        args,
        "aggregation_method",
        &[
            "mean",
            "majority",
            "maximum",
            "median",
            "minimum",
            "minority",
            "percentile",
            "range",
            "std",
            "sum",
            "variety",
        ],
        "mean",
    )? {
        "majority" => Method::Majority,
        "maximum" => Method::Maximum,
        "median" => Method::Median,
        "minimum" => Method::Minimum,
        "minority" => Method::Minority,
        "percentile" => Method::Percentile,
        "range" => Method::Range,
        "std" => Method::Std,
        "sum" => Method::Sum,
        "variety" => Method::Variety,
        _ => Method::Mean,
    };

    let percentile_value = f64_or(args, "percentile_value", 50.0)?;
    if !(0.0..=100.0).contains(&percentile_value) {
        return Err(ToolError::Validation(format!(
            "'percentile_value' must be in [0, 100], got {percentile_value}"
        )));
    }
    let ignore_nodata = bool_or(args, "ignore_nodata", true)?;

    Ok(Params {
        definition,
        interval_value,
        interval_count,
        ranges,
        method,
        percentile_value,
        ignore_nodata,
    })
}

/// Parses `"start:end,start:end"` into ordered, non-degenerate spans.
fn parse_ranges(args: &ToolArgs) -> Result<Option<Vec<(f64, f64)>>, ToolError> {
    let Some(s) = args.get("interval_ranges").and_then(Value::as_str) else {
        return Ok(None);
    };
    if s.trim().is_empty() {
        return Ok(None);
    }
    let mut out = Vec::new();
    for part in s.split(',').map(str::trim).filter(|p| !p.is_empty()) {
        let (a, b) = part.split_once(':').ok_or_else(|| {
            ToolError::Validation(format!(
                "'interval_ranges' entry '{part}' must be 'start:end'"
            ))
        })?;
        let start: f64 = a.trim().parse().map_err(|_| {
            ToolError::Validation(format!("'interval_ranges' start '{a}' is not a number"))
        })?;
        let end: f64 = b.trim().parse().map_err(|_| {
            ToolError::Validation(format!("'interval_ranges' end '{b}' is not a number"))
        })?;
        if !(start.is_finite() && end.is_finite()) || end <= start {
            return Err(ToolError::Validation(format!(
                "'interval_ranges' entry '{part}' must have end greater than start"
            )));
        }
        out.push((start, end));
    }
    if out.is_empty() {
        return Ok(None);
    }
    Ok(Some(out))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::load_input_raster;
    use crate::cube::test_support::cube_raster;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::Raster;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn run(args: Value) -> (Raster, BTreeMap<String, Value>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = AggregateMultidimensionalRasterTool.run(&args, &ctx()).unwrap();
        let r = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (r, out.outputs)
    }

    /// A 1x1 cube of six slices, values 1..6.
    fn six_slices() -> String {
        cube_raster(
            1,
            1,
            &[
                vec![1.0],
                vec![2.0],
                vec![3.0],
                vec![4.0],
                vec![5.0],
                vec![6.0],
            ],
        )
    }

    /// The default collapses the whole cube to one band, matching
    /// `cell_statistics`.
    #[test]
    fn all_reduces_to_one_band() {
        let (out, outputs) = run(json!({ "input": six_slices() }));
        assert_eq!(out.bands, 1);
        assert_eq!(outputs["bin_count"].as_u64().unwrap(), 1);
        assert!((out.get(0, 0, 0) - 3.5).abs() < 1e-6, "mean of 1..6 is 3.5");
    }

    /// This is the operation the tool exists for: many slices into a few
    /// composites, one output band per bin.
    #[test]
    fn interval_value_bins_by_coordinate() {
        // Slice indices 1..6 as coordinates; width 2 gives [1,3) [3,5) [5,7].
        let (out, outputs) = run(json!({
            "input": six_slices(),
            "aggregation_definition": "interval_value", "interval_value": 2
        }));
        assert_eq!(out.bands, 3, "expected three bins");
        assert_eq!(outputs["bin_count"].as_u64().unwrap(), 3);
        assert!((out.get(0, 0, 0) - 1.5).abs() < 1e-6, "bin 1: mean(1,2)");
        assert!((out.get(1, 0, 0) - 3.5).abs() < 1e-6, "bin 2: mean(3,4)");
        assert!((out.get(2, 0, 0) - 5.5).abs() < 1e-6, "bin 3: mean(5,6)");
    }

    /// Real dimension coordinates (years) drive the binning, not slice order.
    #[test]
    fn dimension_values_drive_the_bins() {
        // Two slices in 2020, one in 2021, three in 2022.
        let (out, outputs) = run(json!({
            "input": six_slices(),
            "dimension": "year",
            "dimension_values": "2020.0, 2020.5, 2021.0, 2022.0, 2022.3, 2022.6",
            "aggregation_definition": "interval_value", "interval_value": 1.0
        }));
        assert_eq!(outputs["dimension"].as_str().unwrap(), "year");
        assert_eq!(out.bands, 3, "three calendar years");
        let bins = outputs["bins"].as_array().unwrap();
        assert_eq!(bins[0]["slice_count"].as_u64().unwrap(), 2);
        assert_eq!(bins[1]["slice_count"].as_u64().unwrap(), 1);
        assert_eq!(bins[2]["slice_count"].as_u64().unwrap(), 3);
        assert!((out.get(0, 0, 0) - 1.5).abs() < 1e-6);
        assert!((out.get(1, 0, 0) - 3.0).abs() < 1e-6);
        assert!((out.get(2, 0, 0) - 5.0).abs() < 1e-6, "mean(4,5,6)");
    }

    #[test]
    fn interval_count_splits_the_span_evenly() {
        let (out, _) = run(json!({
            "input": six_slices(),
            "aggregation_definition": "interval_count", "interval_count": 2
        }));
        assert_eq!(out.bands, 2);
        // Coordinates 1..6 split at 3.5 -> {1,2,3} and {4,5,6}.
        assert!((out.get(0, 0, 0) - 2.0).abs() < 1e-6);
        assert!((out.get(1, 0, 0) - 5.0).abs() < 1e-6);
    }

    #[test]
    fn interval_ranges_are_explicit() {
        let (out, outputs) = run(json!({
            "input": six_slices(),
            "aggregation_definition": "interval_ranges",
            "interval_ranges": "1:3, 4:7"
        }));
        assert_eq!(out.bands, 2);
        assert!((out.get(0, 0, 0) - 1.5).abs() < 1e-6, "[1,3) = 1,2");
        assert!((out.get(1, 0, 0) - 5.0).abs() < 1e-6, "[4,7] = 4,5,6");
        assert_eq!(outputs["bins"][0]["slice_count"].as_u64().unwrap(), 2);
    }

    /// Every reducer is available per bin, not just the mean.
    #[test]
    fn reducers_are_applied_per_bin() {
        let src = six_slices();
        let max = run(json!({
            "input": src.clone(), "aggregation_definition": "interval_value",
            "interval_value": 3, "aggregation_method": "maximum"
        }))
        .0;
        assert_eq!(max.get(0, 0, 0), 3.0);
        assert_eq!(max.get(1, 0, 0), 6.0);

        let sum = run(json!({
            "input": src.clone(), "aggregation_definition": "interval_value",
            "interval_value": 3, "aggregation_method": "sum"
        }))
        .0;
        assert_eq!(sum.get(0, 0, 0), 6.0);
        assert_eq!(sum.get(1, 0, 0), 15.0);

        let p25 = run(json!({
            "input": src, "aggregation_method": "percentile", "percentile_value": 25
        }))
        .0;
        // 25th percentile of 1..6 with linear interpolation is 2.25.
        assert!((p25.get(0, 0, 0) - 2.25).abs() < 1e-6);
    }

    /// No-data handling matches `cell_statistics`: skipped by default, fatal to
    /// the bin when `ignore_nodata` is false.
    #[test]
    fn nodata_handling() {
        let src = cube_raster(1, 1, &[vec![2.0], vec![-9999.0], vec![4.0]]);
        let ignore = run(json!({ "input": src.clone() })).0;
        assert!((ignore.get(0, 0, 0) - 3.0).abs() < 1e-6);
        let strict = run(json!({ "input": src, "ignore_nodata": false })).0;
        assert_eq!(strict.get(0, 0, 0), -9999.0);
    }

    /// Empty bins are dropped so every output band has data behind it.
    #[test]
    fn empty_bins_are_dropped() {
        let (out, outputs) = run(json!({
            "input": six_slices(),
            "aggregation_definition": "interval_ranges",
            "interval_ranges": "1:3, 100:200, 4:7"
        }));
        assert_eq!(out.bands, 2, "the 100:200 bin holds no slices");
        assert_eq!(outputs["bin_count"].as_u64().unwrap(), 2);
    }

    /// Cells vary independently across the grid.
    #[test]
    fn works_per_cell() {
        // 2 cells, 4 slices. Cell 0 counts up, cell 1 counts down.
        let src = cube_raster(
            2,
            1,
            &[
                vec![1.0, 40.0],
                vec![2.0, 30.0],
                vec![3.0, 20.0],
                vec![4.0, 10.0],
            ],
        );
        let (out, _) = run(json!({
            "input": src, "aggregation_definition": "interval_value", "interval_value": 2
        }));
        assert_eq!(out.bands, 2);
        assert!((out.get(0, 0, 0) - 1.5).abs() < 1e-6);
        assert!((out.get(0, 0, 1) - 35.0).abs() < 1e-6);
        assert!((out.get(1, 0, 0) - 3.5).abs() < 1e-6);
        assert!((out.get(1, 0, 1) - 15.0).abs() < 1e-6);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            AggregateMultidimensionalRasterTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        // Each definition needs its own parameter rather than silently
        // collapsing to 'all'.
        assert!(bad(json!({"input": "a.tif", "aggregation_definition": "interval_value"})).is_err());
        assert!(bad(json!({"input": "a.tif", "aggregation_definition": "interval_count"})).is_err());
        assert!(bad(json!({"input": "a.tif", "aggregation_definition": "interval_ranges"})).is_err());
        assert!(bad(
            json!({"input": "a.tif", "aggregation_definition": "interval_value", "interval_value": 0})
        )
        .is_err());
        assert!(bad(json!({"input": "a.tif", "interval_ranges": "5:3"})).is_err());
        assert!(bad(json!({"input": "a.tif", "interval_ranges": "bogus"})).is_err());
        assert!(bad(json!({"input": "a.tif", "aggregation_method": "mode"})).is_err());
        assert!(bad(json!({"input": "a.tif", "percentile_value": 101})).is_err());
        assert!(bad(json!({"input": "a.tif"})).is_ok());
    }
}
