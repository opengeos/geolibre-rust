//! Slice a continuous raster into N ordinal zones (ArcGIS `Slice`).
//!
//! A single-band continuous surface is classified into integer zones using one
//! of five methods:
//!   * **Equal interval** — the value range is split into `number_zones` bands of
//!     equal width.
//!   * **Equal area (quantile)** — breaks fall on the data quantiles so each zone
//!     holds (as nearly as possible) the same number of cells.
//!   * **Natural breaks (Jenks)** — Fisher–Jenks dynamic programming on a binned
//!     histogram, minimising within-class variance. O(k·bins²), deterministic,
//!     and cheap enough for WASM.
//!   * **Geometric interval** — class widths grow as a geometric series (log
//!     spacing over the range), good for skewed / exponential data.
//!   * **Standard deviation** — breaks at `mean ± m·(σ·class_interval_size)` for
//!     integer `m`, clipped to the data range; `number_zones` is data-driven.
//!
//! Output is an I32 zone raster; zone ids run from `base_output_zone` upward.

use std::collections::BTreeMap;

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

/// No-data value written to the output zone raster (well outside any zone id).
const OUT_NODATA: f64 = -9999.0;
/// Histogram resolution for the Fisher–Jenks dynamic program.
const JENKS_BINS: usize = 256;

/// Classification method for [`SliceRasterTool`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SliceType {
    EqualInterval,
    EqualArea,
    NaturalBreaks,
    GeometricInterval,
    StdDev,
}

impl SliceType {
    /// Parses a slice-type string, accepting ArcGIS-style and snake_case spellings.
    pub fn parse(s: &str) -> Result<SliceType, ToolError> {
        let key: String = s
            .trim()
            .to_ascii_lowercase()
            .chars()
            .map(|c| if c == ' ' || c == '-' { '_' } else { c })
            .collect();
        match key.as_str() {
            "equal_interval" | "equal" | "interval" => Ok(SliceType::EqualInterval),
            "equal_area" | "quantile" | "quantiles" => Ok(SliceType::EqualArea),
            "natural_breaks" | "natural_break" | "jenks" | "naturalbreaks" => {
                Ok(SliceType::NaturalBreaks)
            }
            "geometric_interval" | "geometric" | "geometrical_interval" => {
                Ok(SliceType::GeometricInterval)
            }
            "std_dev" | "stddev" | "standard_deviation" | "std_deviation" | "sd" => {
                Ok(SliceType::StdDev)
            }
            other => Err(ToolError::Validation(format!(
                "unknown slice_type '{other}' (expected one of: equal_interval, equal_area, \
                 natural_breaks, geometric_interval, std_dev)"
            ))),
        }
    }
}

/// Parsed, validated parameters.
struct Params {
    band: isize,
    number_zones: usize,
    slice_type: SliceType,
    base_output_zone: i64,
    class_interval_size: f64,
}

/// Parses an optional numeric parameter, accepting a JSON number or a numeric
/// string (host UIs post strings).
fn parse_optional_number(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s.trim().parse::<f64>().map(Some).map_err(|_| {
            ToolError::Validation(format!("parameter '{key}' must be a number, got '{s}'"))
        }),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a number"
        ))),
    }
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let band_1based = parse_optional_number(args, "band")?.unwrap_or(1.0).max(1.0) as usize;
    let number_zones = parse_optional_number(args, "number_zones")?.unwrap_or(10.0);
    if number_zones < 2.0 || number_zones.fract() != 0.0 {
        return Err(ToolError::Validation(
            "parameter 'number_zones' must be an integer >= 2".to_string(),
        ));
    }
    let slice_type = match args.get("slice_type").and_then(Value::as_str) {
        Some(s) => SliceType::parse(s)?,
        None => SliceType::EqualInterval,
    };
    let base_output_zone = parse_optional_number(args, "base_output_zone")?.unwrap_or(1.0) as i64;
    let class_interval_size = parse_optional_number(args, "class_interval_size")?.unwrap_or(1.0);
    if class_interval_size <= 0.0 {
        return Err(ToolError::Validation(
            "parameter 'class_interval_size' must be > 0".to_string(),
        ));
    }
    Ok(Params {
        band: (band_1based - 1) as isize,
        number_zones: number_zones as usize,
        slice_type,
        base_output_zone,
        class_interval_size,
    })
}

/// Interior break values (ascending) that separate the zones. `k` zones produce
/// at most `k-1` interior breaks; a value `v` lands in zone `i` when
/// `breaks[i-1] < v <= breaks[i]`.
pub fn compute_breaks(
    sorted_valid: &[f64],
    number_zones: usize,
    slice_type: SliceType,
    class_interval_size: f64,
) -> Vec<f64> {
    if sorted_valid.is_empty() {
        return Vec::new();
    }
    let min = sorted_valid[0];
    let max = sorted_valid[sorted_valid.len() - 1];
    if max <= min {
        return Vec::new();
    }
    match slice_type {
        SliceType::EqualInterval => {
            let step = (max - min) / number_zones as f64;
            (1..number_zones).map(|i| min + step * i as f64).collect()
        }
        SliceType::EqualArea => quantile_breaks(sorted_valid, number_zones),
        SliceType::NaturalBreaks => jenks_breaks(sorted_valid, number_zones, JENKS_BINS),
        SliceType::GeometricInterval => geometric_breaks(min, max, number_zones),
        SliceType::StdDev => {
            let (mean, std) = mean_std(sorted_valid);
            stddev_breaks(min, max, mean, std, class_interval_size)
        }
    }
}

/// Equal-area (quantile) breaks: the `i/k` quantiles of the sorted data.
fn quantile_breaks(sorted: &[f64], k: usize) -> Vec<f64> {
    let n = sorted.len();
    let mut breaks = Vec::with_capacity(k.saturating_sub(1));
    for i in 1..k {
        // Nearest-rank quantile position.
        let pos = (i as f64 / k as f64) * n as f64;
        let idx = (pos.ceil() as usize).clamp(1, n) - 1;
        breaks.push(sorted[idx]);
    }
    dedup_monotonic(breaks)
}

/// Geometric-interval breaks: interval widths grow as a geometric series. Data
/// are shifted to be strictly positive so the log spacing is well defined.
fn geometric_breaks(min: f64, max: f64, k: usize) -> Vec<f64> {
    let shift = if min <= 0.0 { 1.0 - min } else { 0.0 };
    let smin = min + shift;
    let smax = max + shift;
    let ratio = (smax / smin).powf(1.0 / k as f64);
    (1..k)
        .map(|i| smin * ratio.powi(i as i32) - shift)
        .collect()
}

/// Standard-deviation breaks at `mean ± m·(σ·iv)` for integer `m`, strictly
/// inside `(min, max)`. `number_zones` is ignored; the count follows the spread.
fn stddev_breaks(min: f64, max: f64, mean: f64, std: f64, iv: f64) -> Vec<f64> {
    if std <= 0.0 {
        return Vec::new();
    }
    let width = std * iv;
    let mut breaks = Vec::new();
    // m = 0 places a break at the mean; expand outward until outside the range.
    for m in 0.. {
        let v = mean + m as f64 * width;
        if v >= max {
            break;
        }
        if v > min {
            breaks.push(v);
        }
    }
    for m in 1.. {
        let v = mean - m as f64 * width;
        if v <= min {
            break;
        }
        if v < max {
            breaks.push(v);
        }
    }
    breaks.sort_by(|a, b| a.partial_cmp(b).unwrap());
    dedup_monotonic(breaks)
}

/// Mean and population standard deviation of the values.
fn mean_std(values: &[f64]) -> (f64, f64) {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let var = values.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / n;
    (mean, var.sqrt())
}

/// Fisher–Jenks natural breaks via dynamic programming over a binned histogram.
/// Returns up to `k-1` interior break values (bin edges) in ascending order.
#[allow(clippy::needless_range_loop)] // index-driven DP over the cost/split tables
fn jenks_breaks(sorted: &[f64], k: usize, nbins: usize) -> Vec<f64> {
    let min = sorted[0];
    let max = sorted[sorted.len() - 1];
    let nbins = nbins.max(k);
    let bin_width = (max - min) / nbins as f64;

    // Weighted histogram: weight[b] = number of cells in bin b.
    let mut weight = vec![0.0_f64; nbins];
    for &v in sorted {
        let mut b = ((v - min) / bin_width) as usize;
        if b >= nbins {
            b = nbins - 1;
        }
        weight[b] += 1.0;
    }
    let center = |b: usize| min + (b as f64 + 0.5) * bin_width;

    // Prefix sums over bins (1-based) for O(1) weighted sum-of-squared-deviations.
    let mut pw = vec![0.0_f64; nbins + 1];
    let mut pwx = vec![0.0_f64; nbins + 1];
    let mut pwx2 = vec![0.0_f64; nbins + 1];
    for b in 0..nbins {
        let w = weight[b];
        let x = center(b);
        pw[b + 1] = pw[b] + w;
        pwx[b + 1] = pwx[b] + w * x;
        pwx2[b + 1] = pwx2[b] + w * x * x;
    }
    // SSD of contiguous bins i..=j (1-based, inclusive).
    let ssd = |i: usize, j: usize| -> f64 {
        let w = pw[j] - pw[i - 1];
        if w <= 0.0 {
            return 0.0;
        }
        let wx = pwx[j] - pwx[i - 1];
        let wx2 = pwx2[j] - pwx2[i - 1];
        (wx2 - wx * wx / w).max(0.0)
    };

    let n = nbins;
    let k = k.min(n);
    // cost[i] for the current class-count; split[m][i] = optimal previous boundary.
    let mut prev = vec![f64::INFINITY; n + 1];
    for i in 1..=n {
        prev[i] = ssd(1, i);
    }
    let mut splits: Vec<Vec<usize>> = vec![vec![0; n + 1]; k + 1];
    for m in 2..=k {
        let mut cur = vec![f64::INFINITY; n + 1];
        for i in m..=n {
            let mut best = f64::INFINITY;
            let mut best_j = m - 1;
            for j in (m - 1)..=(i - 1) {
                let c = prev[j] + ssd(j + 1, i);
                if c < best {
                    best = c;
                    best_j = j;
                }
            }
            cur[i] = best;
            splits[m][i] = best_j;
        }
        prev = cur;
    }

    // Backtrack the interior boundaries (1-based bin indices) into break values.
    // A boundary after bin `j` splits weight[0..j] from weight[j..]; snapping the
    // break to the midpoint of the empty gap between the two populated bins it
    // separates keeps the boundary "natural" even when bins are empty (ties).
    let mut breaks = Vec::with_capacity(k.saturating_sub(1));
    let mut i = n;
    for m in (2..=k).rev() {
        let j = splits[m][i];
        let lo = (0..j).rev().find(|&b| weight[b] > 0.0);
        let hi = (j..n).find(|&b| weight[b] > 0.0);
        let brk = match (lo, hi) {
            (Some(lo), Some(hi)) => (center(lo) + center(hi)) / 2.0,
            _ => min + j as f64 * bin_width,
        };
        breaks.push(brk);
        i = j;
    }
    breaks.reverse();
    dedup_monotonic(breaks)
}

/// Drops non-increasing / duplicate breaks so the sequence is strictly ascending.
fn dedup_monotonic(breaks: Vec<f64>) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::with_capacity(breaks.len());
    for b in breaks {
        if out.last().map(|&last| b > last).unwrap_or(true) {
            out.push(b);
        }
    }
    out
}

/// Assigns a value to a zone id given ascending interior `breaks`.
fn zone_of(v: f64, breaks: &[f64], base: i64) -> i64 {
    let offset = breaks.partition_point(|&b| b < v);
    base + offset as i64
}

/// Classify a continuous raster band into ordinal zones (ArcGIS `Slice`).
pub struct SliceRasterTool;

impl Tool for SliceRasterTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "slice_raster",
            display_name: "Slice Raster",
            summary: "Classify a continuous raster into N ordinal zones (equal interval, equal area, natural breaks, geometric interval, or standard deviation).",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input continuous raster file path.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output zone raster path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to classify (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "number_zones",
                    description: "Number of output zones (default 10). Ignored for std_dev.",
                    required: false,
                },
                ToolParamSpec {
                    name: "slice_type",
                    description: "Classification method: equal_interval, equal_area, natural_breaks, geometric_interval, or std_dev (default equal_interval).",
                    required: false,
                },
                ToolParamSpec {
                    name: "base_output_zone",
                    description: "Zone id assigned to the lowest zone (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "class_interval_size",
                    description: "Std-dev interval size in σ units (default 1.0). Only used by std_dev.",
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
        if params.band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "band {} out of range (raster has {} band(s))",
                params.band + 1,
                raster.bands
            )));
        }
        let nodata = raster.nodata;
        let cells = band_to_vec(&raster, params.band);

        ctx.progress.info("gathering valid values");
        let mut sorted: Vec<f64> = cells
            .iter()
            .copied()
            .filter(|&v| v != nodata && v.is_finite())
            .collect();
        if sorted.is_empty() {
            return Err(ToolError::Execution(
                "raster band contains no valid (non-nodata) values".to_string(),
            ));
        }
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

        ctx.progress.info("computing class breaks");
        let breaks = compute_breaks(
            &sorted,
            params.number_zones,
            params.slice_type,
            params.class_interval_size,
        );

        ctx.progress.info("assigning zones");
        let mut zone_data = vec![OUT_NODATA; cells.len()];
        let n_zones = breaks.len() + 1;
        let mut counts = vec![0_u64; n_zones];
        let mut valid_cells = 0_u64;
        for (i, &v) in cells.iter().enumerate() {
            if v != nodata && v.is_finite() {
                let zone = zone_of(v, &breaks, params.base_output_zone);
                zone_data[i] = zone as f64;
                counts[(zone - params.base_output_zone) as usize] += 1;
                valid_cells += 1;
            }
        }
        ctx.progress.progress(0.9);

        let out_raster = raster_like_with_data(&raster, zone_data, OUT_NODATA, DataType::I32)?;
        let out_path = write_or_store_output(out_raster, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("zones".to_string(), json!(n_zones));
        outputs.insert("breaks".to_string(), json!(breaks));
        outputs.insert("zone_counts".to_string(), json!(counts));
        outputs.insert("valid_cells".to_string(), json!(valid_cells));
        outputs.insert("min".to_string(), json!(sorted[0]));
        outputs.insert("max".to_string(), json!(sorted[sorted.len() - 1]));
        ctx.progress.progress(1.0);
        Ok(ToolRunResult { outputs })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn raster_from(cols: usize, rows: usize, data: Vec<f64>, nodata: f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata,
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

    fn sorted(data: &[f64]) -> Vec<f64> {
        let mut v = data.to_vec();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v
    }

    fn is_ascending(breaks: &[f64]) -> bool {
        breaks.windows(2).all(|w| w[1] > w[0])
    }

    #[test]
    fn equal_interval_splits_range_evenly() {
        let data: Vec<f64> = (0..=100).map(|i| i as f64).collect();
        let breaks = compute_breaks(&sorted(&data), 4, SliceType::EqualInterval, 1.0);
        assert_eq!(breaks.len(), 3);
        assert!((breaks[0] - 25.0).abs() < 1e-9);
        assert!((breaks[1] - 50.0).abs() < 1e-9);
        assert!((breaks[2] - 75.0).abs() < 1e-9);
        assert!(is_ascending(&breaks));
    }

    #[test]
    fn equal_area_gives_balanced_zone_counts() {
        // 100 distinct values, 5 quantile zones -> ~20 cells each.
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let s = sorted(&data);
        let breaks = compute_breaks(&s, 5, SliceType::EqualArea, 1.0);
        assert!(is_ascending(&breaks));
        let mut counts = vec![0usize; breaks.len() + 1];
        for &v in &data {
            counts[zone_of(v, &breaks, 1) as usize - 1] += 1;
        }
        for c in &counts {
            assert!(
                (*c as i64 - 20).abs() <= 1,
                "unbalanced zone counts: {counts:?}"
            );
        }
    }

    #[test]
    fn natural_breaks_separates_clusters() {
        // Two tight clusters around 0 and 100 -> one break between them.
        let mut data = vec![0.5_f64; 50];
        data.extend(vec![99.5_f64; 50]);
        let breaks = compute_breaks(&sorted(&data), 2, SliceType::NaturalBreaks, 1.0);
        assert_eq!(breaks.len(), 1);
        assert!(
            breaks[0] > 5.0 && breaks[0] < 95.0,
            "break not between clusters: {breaks:?}"
        );
    }

    #[test]
    fn geometric_interval_widths_grow() {
        let breaks = compute_breaks(
            &sorted(&[1.0, 1000.0]),
            3,
            SliceType::GeometricInterval,
            1.0,
        );
        assert_eq!(breaks.len(), 2);
        assert!(is_ascending(&breaks));
        // Later interval is wider than the first (geometric growth).
        let w0 = breaks[0] - 1.0;
        let w1 = breaks[1] - breaks[0];
        assert!(w1 > w0, "geometric widths did not grow: {breaks:?}");
    }

    #[test]
    fn std_dev_breaks_around_mean() {
        // Symmetric spread: mean 50, breaks straddle the mean.
        let data: Vec<f64> = (0..=100).map(|i| i as f64).collect();
        let breaks = compute_breaks(&sorted(&data), 10, SliceType::StdDev, 1.0);
        assert!(is_ascending(&breaks));
        assert!(
            breaks.iter().any(|&b| (b - 50.0).abs() < 1e-6),
            "no break at the mean: {breaks:?}"
        );
    }

    #[test]
    fn zone_assignment_covers_full_range() {
        let data: Vec<f64> = (0..=100).map(|i| i as f64).collect();
        let breaks = compute_breaks(&sorted(&data), 4, SliceType::EqualInterval, 1.0);
        assert_eq!(zone_of(0.0, &breaks, 1), 1);
        assert_eq!(zone_of(100.0, &breaks, 1), 4);
    }

    #[test]
    fn end_to_end_zone_counts_sum_to_valid_cells() {
        let nodata = -9999.0;
        let mut data: Vec<f64> = (0..99).map(|i| i as f64).collect();
        data.push(nodata); // one nodata cell
        let path = raster_from(10, 10, data, nodata);
        let out = SliceRasterTool
            .run(
                &serde_json::from_value(json!({
                    "input": path,
                    "number_zones": 5,
                    "slice_type": "equal_area"
                }))
                .unwrap(),
                &ctx(),
            )
            .unwrap();
        let counts = out.outputs["zone_counts"].as_array().unwrap();
        let sum: u64 = counts.iter().map(|c| c.as_u64().unwrap()).sum();
        assert_eq!(sum, out.outputs["valid_cells"].as_u64().unwrap());
        assert_eq!(sum, 99);
        // Every zone is non-empty for a quantile classification.
        assert!(counts.iter().all(|c| c.as_u64().unwrap() > 0));
        // Read back the raster and confirm nodata is preserved.
        let r = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        assert_eq!(r.get(0, 9, 9), OUT_NODATA);
    }

    #[test]
    fn rejects_bad_parameters() {
        assert!(SliceRasterTool
            .validate(&serde_json::from_value(json!({})).unwrap())
            .is_err());
        // number_zones < 2
        assert!(SliceRasterTool
            .validate(&serde_json::from_value(json!({"input": "x", "number_zones": 1})).unwrap())
            .is_err());
        // unknown slice_type
        assert!(SliceRasterTool
            .validate(
                &serde_json::from_value(json!({"input": "x", "slice_type": "bogus"})).unwrap()
            )
            .is_err());
        // numeric string is accepted
        assert!(SliceRasterTool
            .validate(&serde_json::from_value(json!({"input": "x", "number_zones": "8"})).unwrap())
            .is_ok());
    }
}
