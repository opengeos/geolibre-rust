//! GeoLibre tool: non-overlapping block reduction at unchanged resolution.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Block Statistics* (Spatial Analyst).
//!
//! This fills the missing quadrant of the neighbourhood-reduction family:
//!
//! * `focal_statistics` reduces an **overlapping** moving window;
//! * `zonal_statistics` / `zonal_characterization` reduce by an explicit **zone
//!   layer**;
//! * `cell_statistics` reduces **across a stack** of aligned rasters;
//! * this reduces over **non-overlapping blocks**, writing each block's value
//!   back to every cell in the block at the input resolution.
//!
//! The bundled `aggregate_raster` is not a substitute: it *reduces resolution*
//! by the aggregation factor and offers only mean/sum/min/max/range. Keeping the
//! original grid is the whole point here — `input - block_statistics(input,
//! mean)` is the standard local-residual surface, and it cannot be expressed
//! with `aggregate_raster` without a resample round-trip that reintroduces
//! interpolation error. The categorical statistics (majority, minority, median,
//! variety) are likewise absent there.
//!
//! Neighbourhood shapes follow ArcGIS: `rectangle`, `circle`, `annulus` and
//! `wedge`. The `irregular`/`weight` kernel-file shapes are **not** supported —
//! the repo has no kernel-file format — and are rejected with a clear error
//! rather than silently falling back to a rectangle.

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};

/// Reduces non-overlapping blocks of a raster and broadcasts each result back
/// across its block.
pub struct BlockStatisticsTool;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Statistic {
    Mean,
    Majority,
    Maximum,
    Median,
    Minimum,
    Minority,
    Range,
    Std,
    Sum,
    Variety,
}

impl Statistic {
    fn parse(s: &str) -> Result<Self, ToolError> {
        match s.trim().to_ascii_lowercase().as_str() {
            "mean" => Ok(Self::Mean),
            "majority" => Ok(Self::Majority),
            "maximum" | "max" => Ok(Self::Maximum),
            "median" => Ok(Self::Median),
            "minimum" | "min" => Ok(Self::Minimum),
            "minority" => Ok(Self::Minority),
            "range" => Ok(Self::Range),
            "std" | "stddev" | "standard_deviation" => Ok(Self::Std),
            "sum" => Ok(Self::Sum),
            "variety" => Ok(Self::Variety),
            other => Err(ToolError::Validation(format!(
                "unknown statistic '{other}' (expected one of: mean, majority, maximum, median, minimum, minority, range, std, sum, variety)"
            ))),
        }
    }

    /// Applies the statistic to a block's valid values. `values` may be
    /// reordered. Returns `None` for an empty block.
    fn apply(self, values: &mut [f64]) -> Option<f64> {
        if values.is_empty() {
            return None;
        }
        let n = values.len() as f64;
        Some(match self {
            Self::Mean => values.iter().sum::<f64>() / n,
            Self::Sum => values.iter().sum::<f64>(),
            Self::Maximum => values.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            Self::Minimum => values.iter().copied().fold(f64::INFINITY, f64::min),
            Self::Range => {
                let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                let min = values.iter().copied().fold(f64::INFINITY, f64::min);
                max - min
            }
            Self::Std => {
                let mean = values.iter().sum::<f64>() / n;
                (values.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / n).sqrt()
            }
            Self::Median => {
                values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                let mid = values.len() / 2;
                if values.len() % 2 == 1 {
                    values[mid]
                } else {
                    (values[mid - 1] + values[mid]) / 2.0
                }
            }
            Self::Variety => distinct_counts(values).len() as f64,
            Self::Majority => extreme_by_count(values, true)?,
            Self::Minority => extreme_by_count(values, false)?,
        })
    }
}

/// Groups equal values, returning `(value, count)` sorted ascending by value so
/// ties break deterministically (lowest value wins) rather than by hash order.
fn distinct_counts(values: &mut [f64]) -> Vec<(f64, usize)> {
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut out: Vec<(f64, usize)> = Vec::new();
    for &v in values.iter() {
        match out.last_mut() {
            Some((last, count)) if *last == v => *count += 1,
            _ => out.push((v, 1)),
        }
    }
    out
}

/// Most (or least) frequent value. Ties resolve to the lowest value, which keeps
/// results reproducible across platforms and WASM.
fn extreme_by_count(values: &mut [f64], most: bool) -> Option<f64> {
    let counts = distinct_counts(values);
    let mut best: Option<(f64, usize)> = None;
    for (v, c) in counts {
        best = Some(match best {
            None => (v, c),
            Some((bv, bc)) => {
                let better = if most { c > bc } else { c < bc };
                if better {
                    (v, c)
                } else {
                    (bv, bc)
                }
            }
        });
    }
    best.map(|(v, _)| v)
}

#[derive(Clone, Copy, Debug)]
enum Neighborhood {
    Rectangle {
        width: usize,
        height: usize,
    },
    Circle {
        radius: f64,
    },
    Annulus {
        inner: f64,
        outer: f64,
    },
    Wedge {
        radius: f64,
        start_deg: f64,
        end_deg: f64,
    },
}

impl Neighborhood {
    /// Block extent in cells: the bounding box the block tiles the raster with.
    fn block_size(self) -> (usize, usize) {
        match self {
            Self::Rectangle { width, height } => (width, height),
            Self::Circle { radius } | Self::Wedge { radius, .. } => {
                let d = (radius * 2.0).floor() as usize + 1;
                (d, d)
            }
            Self::Annulus { outer, .. } => {
                let d = (outer * 2.0).floor() as usize + 1;
                (d, d)
            }
        }
    }

    /// Whether a cell at `(d_col, d_row)` from the block's centre participates.
    fn includes(self, d_col: f64, d_row: f64) -> bool {
        match self {
            Self::Rectangle { .. } => true,
            Self::Circle { radius } => d_col * d_col + d_row * d_row <= radius * radius,
            Self::Annulus { inner, outer } => {
                let d2 = d_col * d_col + d_row * d_row;
                d2 <= outer * outer && d2 >= inner * inner
            }
            Self::Wedge {
                radius,
                start_deg,
                end_deg,
            } => {
                let d2 = d_col * d_col + d_row * d_row;
                if d2 > radius * radius {
                    return false;
                }
                if d2 == 0.0 {
                    return true;
                }
                // Screen rows increase downward, so negate d_row to get a
                // conventional counter-clockwise-from-east map angle.
                let mut angle = (-d_row).atan2(d_col).to_degrees();
                if angle < 0.0 {
                    angle += 360.0;
                }
                let start = start_deg.rem_euclid(360.0);
                let mut end = end_deg.rem_euclid(360.0);
                if start == end {
                    return true;
                }
                if end < start {
                    // Wedge wraps through 0 degrees.
                    end += 360.0;
                    if angle < start {
                        angle += 360.0;
                    }
                }
                angle >= start && angle <= end
            }
        }
    }
}

impl Tool for BlockStatisticsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "block_statistics",
            display_name: "Block Statistics",
            summary: "Partitions a raster into non-overlapping blocks, reduces each by a statistic, and writes that value back to every cell in the block at the input resolution (ArcGIS Block Statistics). Complements the shipped focal_statistics (overlapping window), zonal_statistics (by zone layer) and cell_statistics (across a stack); unlike the bundled aggregate_raster it preserves cell size, so input - block mean gives a local residual surface, and it adds the categorical statistics majority/minority/median/variety.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output raster path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "statistic",
                    description: "One of mean (default), majority, maximum, median, minimum, minority, range, std, sum, variety.",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighborhood",
                    description: "Block shape: rectangle (default), circle, annulus, or wedge.",
                    required: false,
                },
                ToolParamSpec {
                    name: "size",
                    description: "Block dimensions in cells: 'width,height' for rectangle (default '3,3'), 'radius' for circle/wedge, or 'inner,outer' for annulus.",
                    required: false,
                },
                ToolParamSpec {
                    name: "start_angle",
                    description: "Wedge start angle in degrees counter-clockwise from east (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "end_angle",
                    description: "Wedge end angle in degrees counter-clockwise from east (default 90).",
                    required: false,
                },
                ToolParamSpec {
                    name: "ignore_nodata",
                    description: "If true (default), blocks are computed from valid cells only; if false, any no-data cell makes the whole block no-data.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read (default 1).",
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
        if params.band == 0 || params.band as usize > raster.bands {
            return Err(ToolError::Validation(format!(
                "band {} out of range (raster has {} band(s))",
                params.band, raster.bands
            )));
        }
        let band = (params.band - 1) as isize;

        let rows = raster.rows;
        let cols = raster.cols;
        let nodata = raster.nodata;
        let (block_w, block_h) = params.neighborhood.block_size();

        ctx.progress.info("reducing blocks");
        let mut data = vec![nodata; rows * cols];
        let mut blocks_written = 0_u64;
        let mut blocks_empty = 0_u64;

        let mut block_row = 0;
        while block_row < rows {
            let mut block_col = 0;
            while block_col < cols {
                // Centre of the block's bounding box, used by the shape tests.
                let centre_col = block_col as f64 + (block_w as f64 - 1.0) / 2.0;
                let centre_row = block_row as f64 + (block_h as f64 - 1.0) / 2.0;

                let mut values: Vec<f64> = Vec::new();
                let mut members: Vec<(usize, usize)> = Vec::new();
                let mut saw_nodata = false;

                for r in block_row..(block_row + block_h).min(rows) {
                    for c in block_col..(block_col + block_w).min(cols) {
                        if !params
                            .neighborhood
                            .includes(c as f64 - centre_col, r as f64 - centre_row)
                        {
                            continue;
                        }
                        members.push((r, c));
                        let v = raster.get(band, r as isize, c as isize);
                        if v == nodata || !v.is_finite() {
                            saw_nodata = true;
                        } else {
                            values.push(v);
                        }
                    }
                }

                let result = if saw_nodata && !params.ignore_nodata {
                    None
                } else {
                    params.statistic.apply(&mut values)
                };

                match result {
                    Some(v) => {
                        for (r, c) in &members {
                            data[r * cols + c] = v;
                        }
                        blocks_written += 1;
                    }
                    None => blocks_empty += 1,
                }

                block_col += block_w;
            }
            block_row += block_h;
            ctx.progress
                .progress((block_row.min(rows) as f64) / rows.max(1) as f64);
        }

        // Variety and majority/minority of integer classes stay integral, but
        // mean/std do not, so always widen to F32 for a faithful result.
        let out = raster_like_with_data(&raster, data, nodata, DataType::F32)?;
        let out_path = write_or_store_output(out, output)?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("blocks_written".to_string(), json!(blocks_written));
        outputs.insert("blocks_empty".to_string(), json!(blocks_empty));
        outputs.insert("block_width".to_string(), json!(block_w));
        outputs.insert("block_height".to_string(), json!(block_h));
        Ok(ToolRunResult { outputs })
    }
}

struct Params {
    statistic: Statistic,
    neighborhood: Neighborhood,
    ignore_nodata: bool,
    band: u64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let statistic = match opt_str(args, "statistic")? {
        Some(s) => Statistic::parse(s)?,
        None => Statistic::Mean,
    };

    let shape = opt_str(args, "neighborhood")?
        .unwrap_or("rectangle")
        .trim()
        .to_ascii_lowercase();
    let size = opt_str(args, "size")?.unwrap_or("3,3");
    let nums = parse_size_list(size)?;

    let neighborhood = match shape.as_str() {
        "rectangle" | "rect" => {
            let (w, h) = match nums.len() {
                0 => (3.0, 3.0),
                1 => (nums[0], nums[0]),
                _ => (nums[0], nums[1]),
            };
            if w < 1.0 || h < 1.0 {
                return Err(ToolError::Validation(
                    "rectangle 'size' must be at least 1 cell in each dimension".to_string(),
                ));
            }
            Neighborhood::Rectangle {
                width: w as usize,
                height: h as usize,
            }
        }
        "circle" => {
            let radius = nums.first().copied().unwrap_or(1.0);
            if radius <= 0.0 {
                return Err(ToolError::Validation(
                    "circle 'size' (radius) must be positive".to_string(),
                ));
            }
            Neighborhood::Circle { radius }
        }
        "annulus" => {
            if nums.len() < 2 {
                return Err(ToolError::Validation(
                    "annulus 'size' must be 'inner,outer'".to_string(),
                ));
            }
            let (inner, outer) = (nums[0], nums[1]);
            if inner < 0.0 || outer <= inner {
                return Err(ToolError::Validation(
                    "annulus 'size' must satisfy 0 <= inner < outer".to_string(),
                ));
            }
            Neighborhood::Annulus { inner, outer }
        }
        "wedge" => {
            let radius = nums.first().copied().unwrap_or(1.0);
            if radius <= 0.0 {
                return Err(ToolError::Validation(
                    "wedge 'size' (radius) must be positive".to_string(),
                ));
            }
            Neighborhood::Wedge {
                radius,
                start_deg: opt_f64(args, "start_angle")?.unwrap_or(0.0),
                end_deg: opt_f64(args, "end_angle")?.unwrap_or(90.0),
            }
        }
        "irregular" | "weight" => {
            return Err(ToolError::Validation(
                "neighborhood 'irregular'/'weight' needs a kernel file, which this tool does not support; use rectangle, circle, annulus or wedge".to_string(),
            ))
        }
        other => {
            return Err(ToolError::Validation(format!(
                "unknown neighborhood '{other}' (expected rectangle, circle, annulus or wedge)"
            )))
        }
    };

    Ok(Params {
        statistic,
        neighborhood,
        ignore_nodata: opt_bool(args, "ignore_nodata")?.unwrap_or(true),
        band: opt_u64(args, "band")?.unwrap_or(1),
    })
}

fn parse_size_list(s: &str) -> Result<Vec<f64>, ToolError> {
    let mut out = Vec::new();
    for part in s.split(',') {
        let t = part.trim();
        if t.is_empty() {
            continue;
        }
        let v = t.parse::<f64>().map_err(|_| {
            ToolError::Validation(format!("parameter 'size' has non-numeric component '{t}'"))
        })?;
        // NaN/inf must be rejected here: `NaN as usize` saturates to 0, which
        // would make `block_size()` return (0, 0) and leave the block loops in
        // `run` unable to advance — an infinite spin on caller-controlled input.
        if !v.is_finite() {
            return Err(ToolError::Validation(format!(
                "parameter 'size' component '{t}' must be finite"
            )));
        }
        out.push(v);
    }
    Ok(out)
}

fn opt_str<'a>(args: &'a ToolArgs, key: &str) -> Result<Option<&'a str>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => Ok(Some(s.as_str())),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a string when provided"
        ))),
    }
}

fn opt_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
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

fn opt_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => n.as_u64().map(Some).ok_or_else(|| {
            ToolError::Validation(format!("parameter '{key}' must be a positive integer"))
        }),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s.trim().parse::<u64>().map(Some).map_err(|_| {
            ToolError::Validation(format!("parameter '{key}' must be a positive integer"))
        }),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a positive integer"
        ))),
    }
}

fn opt_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
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
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
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

    fn raster(cols: usize, rows: usize, data: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
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

    fn run_with(path: String, extra: Value) -> Raster {
        let mut obj = serde_json::Map::new();
        obj.insert("input".to_string(), json!(path));
        if let Value::Object(m) = extra {
            for (k, v) in m {
                obj.insert(k, v);
            }
        }
        let args: ToolArgs = serde_json::from_value(Value::Object(obj)).unwrap();
        let out = BlockStatisticsTool.run(&args, &ctx()).unwrap();
        load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    /// The defining property: output keeps the input's dimensions, and every
    /// cell of a block carries that block's statistic.
    #[test]
    fn broadcasts_block_mean_at_input_resolution() {
        // 4x4, blocks of 2x2. Top-left block values 1,2,5,6 -> mean 3.5.
        let path = raster(
            4,
            4,
            &[
                1.0, 2.0, 3.0, 4.0, //
                5.0, 6.0, 7.0, 8.0, //
                9.0, 10.0, 11.0, 12.0, //
                13.0, 14.0, 15.0, 16.0,
            ],
        );
        let out = run_with(path, json!({ "size": "2,2", "statistic": "mean" }));
        assert_eq!(out.rows, 4, "resolution must be preserved");
        assert_eq!(out.cols, 4);
        for (r, c) in [(0, 0), (0, 1), (1, 0), (1, 1)] {
            assert_eq!(out.get(0, r, c), 3.5);
        }
        // Bottom-right block: 11,12,15,16 -> 13.5
        assert_eq!(out.get(0, 3, 3), 13.5);
    }

    /// Blocks are non-overlapping: neighbouring blocks get different values,
    /// which is what distinguishes this from a moving window.
    #[test]
    fn blocks_do_not_overlap() {
        let path = raster(4, 1, &[0.0, 0.0, 10.0, 10.0]);
        let out = run_with(path, json!({ "size": "2,1", "statistic": "mean" }));
        assert_eq!(out.get(0, 0, 0), 0.0);
        assert_eq!(out.get(0, 0, 1), 0.0);
        assert_eq!(out.get(0, 0, 2), 10.0);
        assert_eq!(out.get(0, 0, 3), 10.0);
    }

    /// Categorical statistics the bundled aggregate_raster cannot do.
    #[test]
    fn categorical_statistics() {
        let path = raster(2, 2, &[5.0, 5.0, 5.0, 7.0]);
        assert_eq!(
            run_with(
                path.clone(),
                json!({ "size": "2,2", "statistic": "majority" })
            )
            .get(0, 0, 0),
            5.0
        );
        assert_eq!(
            run_with(
                path.clone(),
                json!({ "size": "2,2", "statistic": "minority" })
            )
            .get(0, 0, 0),
            7.0
        );
        assert_eq!(
            run_with(
                path.clone(),
                json!({ "size": "2,2", "statistic": "variety" })
            )
            .get(0, 0, 0),
            2.0
        );
        assert_eq!(
            run_with(path, json!({ "size": "2,2", "statistic": "median" })).get(0, 0, 0),
            5.0
        );
    }

    /// Partial edge blocks are computed from whatever cells they contain.
    #[test]
    fn partial_edge_blocks_are_computed() {
        // 3 columns with a 2-wide block leaves a 1-wide remainder.
        let path = raster(3, 1, &[2.0, 4.0, 9.0]);
        let out = run_with(path, json!({ "size": "2,1", "statistic": "mean" }));
        assert_eq!(out.get(0, 0, 0), 3.0);
        assert_eq!(out.get(0, 0, 1), 3.0);
        assert_eq!(out.get(0, 0, 2), 9.0, "remainder block uses its own cells");
    }

    /// ignore_nodata=false poisons the whole block; the default ignores it.
    #[test]
    fn nodata_handling_follows_flag() {
        let path = raster(2, 2, &[4.0, 6.0, 8.0, -9999.0]);
        let ignored = run_with(path.clone(), json!({ "size": "2,2" }));
        assert_eq!(ignored.get(0, 0, 0), 6.0, "mean of 4,6,8");
        let strict = run_with(path, json!({ "size": "2,2", "ignore_nodata": false }));
        assert_eq!(strict.get(0, 0, 0), strict.nodata);
    }

    /// A circular neighbourhood excludes the block's corners.
    #[test]
    fn circle_neighborhood_excludes_corners() {
        // radius 1 over a 3x3 block: the 4 corners fall outside.
        let path = raster(3, 3, &[100.0, 1.0, 100.0, 1.0, 1.0, 1.0, 100.0, 1.0, 100.0]);
        let out = run_with(
            path,
            json!({ "neighborhood": "circle", "size": "1", "statistic": "maximum" }),
        );
        // Only the plus-shaped cells participate, all of which are 1.
        assert_eq!(out.get(0, 1, 1), 1.0);
        // Corners were never members, so they stay no-data.
        assert_eq!(out.get(0, 0, 0), out.nodata);
    }

    #[test]
    fn rejects_bad_parameters() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(BlockStatisticsTool.validate(&args).is_err());

        let path = raster(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        for bad in [
            json!({ "input": path.clone(), "statistic": "nonsense" }),
            json!({ "input": path.clone(), "neighborhood": "hexagon" }),
            json!({ "input": path.clone(), "neighborhood": "annulus", "size": "5" }),
            json!({ "input": path.clone(), "neighborhood": "irregular" }),
            json!({ "input": path.clone(), "size": "a,b" }),
            // NaN/inf would saturate to a zero block size and spin forever.
            json!({ "input": path.clone(), "size": "nan,nan" }),
            json!({ "input": path.clone(), "size": "inf,3" }),
        ] {
            let args: ToolArgs = serde_json::from_value(bad).unwrap();
            assert!(BlockStatisticsTool.validate(&args).is_err());
        }
    }
}
