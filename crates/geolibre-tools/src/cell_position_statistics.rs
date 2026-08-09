//! GeoLibre tool: per-cell *position* and rank statistics across a raster stack.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Highest Position*, *Lowest Position*,
//! *Popularity* and *Rank* (Spatial Analyst / Image Analyst Local toolset),
//! consolidated behind one `statistic` parameter — the same way
//! `cell_statistics` consolidates the value-based Local tools.
//!
//! ## Why this is not `cell_statistics`
//!
//! `cell_statistics` reduces a stack to a **value**: the maximum NDVI across
//! twelve monthly composites. It can never say **which month** that maximum
//! came from, and that is a different question with different uses — phenology,
//! best-month compositing, change timing, "which scenario dominates here".
//! The bundled `find_argument_statistics` works along a single multiband
//! raster's argument axis and has neither the popularity/rank semantics nor the
//! separate-selector-raster input ArcGIS defines.
//!
//! ## Semantics (matching ArcGIS)
//!
//! * `highest_position` / `lowest_position` — the **1-based** layer index
//!   holding the extreme value. Ties resolve to the lowest index.
//! * `popularity` — the value occurring exactly *n* times, where *n* comes from
//!   the `selector` raster or constant. If **no** value occurs *n* times, or if
//!   more than one value does, the cell is no-data. That tie rule is ArcGIS's,
//!   and it is why popularity is not just "majority with a count".
//! * `rank` — the *n*-th smallest value in the cell's stack, 1-based. A rank
//!   outside `1..=count` yields no-data.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::args_common::{bool_or, choice_or, opt_f64};
use crate::common::{load_input_raster, parse_optional_output};
use crate::raster_stack::{
    check_alignment_refs, load_stack, parse_band_policy, parse_input_paths, write_stack_result,
};

pub struct CellPositionStatisticsTool;

impl Tool for CellPositionStatisticsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "cell_position_statistics",
            display_name: "Cell Position Statistics",
            summary: "Per-cell position and rank statistics across a stack of aligned rasters (ArcGIS Highest Position, Lowest Position, Popularity and Rank): which layer holds the extreme value, the value occurring a given number of times, or the nth-smallest value. cell_statistics reduces a stack to a value and can never report which layer that value came from, and the bundled find_argument_statistics works along one multiband raster's argument axis without the popularity/rank selector semantics.",
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
                    description: "Output raster. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "statistic",
                    description: "'highest_position' (default), 'lowest_position', 'popularity', or 'rank'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "selector",
                    description: "Required for 'popularity' and 'rank': a co-registered raster path, or a constant number, giving the per-cell occurrence count (popularity) or 1-based rank to select.",
                    required: false,
                },
                ToolParamSpec {
                    name: "ignore_nodata",
                    description: "Skip no-data observations per cell (default true). When false, any no-data observation makes the cell no-data.",
                    required: false,
                },
                ToolParamSpec {
                    name: "process_as_multiband",
                    description: "'single_band' (default): every band of every input is one layer of a single stack. 'multi_band': band i of every input forms its own stack, giving a multiband result.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        let stat = parse_statistic(args)?;
        parse_input_paths(args, "inputs")?;
        parse_band_policy(args, "process_as_multiband")?;
        // Presence, not string-ness: `run` accepts a bare number through
        // `opt_f64`, and the parameter is documented as "a raster path, or a
        // constant number". Testing `as_str` alone rejected `{"selector": 3}`.
        let selector_given = match args.get("selector") {
            None | Some(Value::Null) => false,
            Some(Value::String(s)) => !s.trim().is_empty(),
            Some(_) => true,
        };
        if stat.needs_selector() && !selector_given {
            return Err(ToolError::Validation(format!(
                "statistic '{}' requires 'selector' (a raster path or a constant number)",
                stat.label()
            )));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let statistic = parse_statistic(args)?;
        let paths = parse_input_paths(args, "inputs")?;
        let policy = parse_band_policy(args, "process_as_multiband")?;
        let ignore_nodata = bool_or(args, "ignore_nodata", true)?;
        let output = parse_optional_output(args, "output")?;

        // Position statistics are meaningful with a single layer (the answer is
        // trivially layer 1), but the tool exists to compare layers, so hold to
        // the same 2-layer floor as cell_statistics.
        let stack = load_stack(&paths, policy, 2)?;
        let (rows, cols) = (stack.rows, stack.cols);
        let selector = parse_selector(args, statistic, &stack)?;

        ctx.progress.info(&format!(
            "{} input(s), {} layer(s), {rows}x{cols}, statistic {}, {}",
            paths.len(),
            stack.total_layers(),
            statistic.label(),
            policy.label()
        ));

        let nodata = -9999.0_f64;
        let n_groups = stack.group_count();
        let mut bands: Vec<Vec<f64>> = Vec::with_capacity(n_groups);
        let mut vals: Vec<f64> = Vec::new();
        let mut resolved = 0usize;

        for g in 0..n_groups {
            let mut out = vec![nodata; rows * cols];
            for r in 0..rows {
                for c in 0..cols {
                    let had_nodata = stack.cell_values(g, r, c, &mut vals);
                    if vals.is_empty() || (!ignore_nodata && had_nodata) {
                        continue;
                    }
                    // Positions are reported against the *layer* order, so the
                    // extreme search must run before any sort. `rank` and
                    // `popularity` sort/count a copy instead.
                    let v = match statistic {
                        Statistic::HighestPosition => Some(extreme_position(&vals, true)),
                        Statistic::LowestPosition => Some(extreme_position(&vals, false)),
                        Statistic::Popularity => {
                            selector_at(&selector, r, c).and_then(|n| popularity(&vals, n))
                        }
                        Statistic::Rank => {
                            selector_at(&selector, r, c).and_then(|n| rank(&mut vals, n))
                        }
                    };
                    if let Some(v) = v {
                        out[r * cols + c] = v;
                        resolved += 1;
                    }
                }
                ctx.progress
                    .progress((g as f64 + (r as f64 + 1.0) / rows as f64) / n_groups as f64);
            }
            bands.push(out);
        }

        // Positions are small integers; popularity/rank return input values, so
        // they must keep float precision.
        let dtype = if statistic.is_position() {
            DataType::I32
        } else {
            DataType::F32
        };
        let out_path = write_stack_result(stack.template(), bands, nodata, dtype, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("statistic".to_string(), json!(statistic.label()));
        outputs.insert("layers".to_string(), json!(stack.total_layers()));
        outputs.insert("bands".to_string(), json!(n_groups));
        outputs.insert("resolved_cells".to_string(), json!(resolved));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        Ok(ToolRunResult { outputs })
    }
}

// ── Reducers ────────────────────────────────────────────────────────────────

/// 1-based index of the maximum (`want_max`) or minimum value. Ties resolve to
/// the lowest index, matching ArcGIS.
fn extreme_position(vals: &[f64], want_max: bool) -> f64 {
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
    best as f64 + 1.0
}

/// The value occurring exactly `n` times. `None` when no value does, or when
/// more than one does (ArcGIS resolves that tie to NoData rather than picking).
fn popularity(vals: &[f64], n: f64) -> Option<f64> {
    if n < 1.0 || n.fract() != 0.0 || n > vals.len() as f64 {
        return None;
    }
    let want = n as usize;
    let mut counts: BTreeMap<u64, (usize, f64)> = BTreeMap::new();
    for &v in vals {
        let e = counts.entry(v.to_bits()).or_insert((0, v));
        e.0 += 1;
    }
    let mut hit = None;
    for (count, value) in counts.values() {
        if *count == want {
            if hit.is_some() {
                return None; // ambiguous → NoData
            }
            hit = Some(*value);
        }
    }
    hit
}

/// The `n`-th smallest value, 1-based. `vals` is sorted in place.
fn rank(vals: &mut [f64], n: f64) -> Option<f64> {
    if n < 1.0 || n.fract() != 0.0 || n > vals.len() as f64 {
        return None;
    }
    vals.sort_by(f64::total_cmp);
    Some(vals[n as usize - 1])
}

// ── Parameters ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Statistic {
    HighestPosition,
    LowestPosition,
    Popularity,
    Rank,
}

impl Statistic {
    fn label(self) -> &'static str {
        match self {
            Statistic::HighestPosition => "highest_position",
            Statistic::LowestPosition => "lowest_position",
            Statistic::Popularity => "popularity",
            Statistic::Rank => "rank",
        }
    }

    fn is_position(self) -> bool {
        matches!(self, Statistic::HighestPosition | Statistic::LowestPosition)
    }

    fn needs_selector(self) -> bool {
        matches!(self, Statistic::Popularity | Statistic::Rank)
    }
}

fn parse_statistic(args: &ToolArgs) -> Result<Statistic, ToolError> {
    let choice = choice_or(
        args,
        "statistic",
        &["highest_position", "lowest_position", "popularity", "rank"],
        "highest_position",
    )?;
    Ok(match choice {
        "highest_position" => Statistic::HighestPosition,
        "lowest_position" => Statistic::LowestPosition,
        "popularity" => Statistic::Popularity,
        _ => Statistic::Rank,
    })
}

/// The popularity/rank selector: a constant, a raster, or absent.
enum Selector {
    None,
    Constant(f64),
    Raster(Box<Raster>),
}

fn parse_selector(
    args: &ToolArgs,
    statistic: Statistic,
    stack: &crate::raster_stack::Stack,
) -> Result<Selector, ToolError> {
    if !statistic.needs_selector() {
        return Ok(Selector::None);
    }
    // A bare number is a constant; anything else is a path. Try the numeric
    // reading first so "3" never gets probed as a filename.
    if let Some(v) = opt_f64(args, "selector").ok().flatten() {
        return Ok(Selector::Constant(v));
    }
    let path = args
        .get("selector")
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| {
            ToolError::Validation(format!(
                "statistic '{}' requires 'selector'",
                statistic.label()
            ))
        })?;
    let raster = load_input_raster(path)?;
    // The selector is read cell-by-cell against the stack grid, so it has to be
    // on that grid; without this check a mismatched selector would silently
    // sample the wrong location.
    check_alignment_refs(&[stack.template(), &raster])?;
    Ok(Selector::Raster(Box::new(raster)))
}

fn selector_at(selector: &Selector, row: usize, col: usize) -> Option<f64> {
    match selector {
        Selector::None => None,
        Selector::Constant(v) => Some(*v),
        Selector::Raster(r) => {
            let v = r.get(0, row as isize, col as isize);
            (v != r.nodata && v.is_finite()).then_some(v)
        }
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

    /// Builds an in-memory raster from per-band, row-major buffers.
    fn multiband(rows: usize, cols: usize, data: Vec<Vec<f64>>) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: data.len(),
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

    /// Single-band convenience wrapper.
    fn raster(data: &[f64], rows: usize, cols: usize) -> String {
        multiband(rows, cols, vec![data.to_vec()])
    }

    fn run(args: Value) -> (Raster, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = CellPositionStatisticsTool.run(&args, &ctx()).unwrap();
        let path = res.outputs["output"].as_str().unwrap().to_string();
        (load_input_raster(&path).unwrap(), res)
    }

    #[test]
    fn highest_position_reports_the_winning_layer_not_its_value() {
        // Three 1x3 layers. Cell 0 peaks in layer 1, cell 1 in layer 2,
        // cell 2 in layer 3 — the whole point of the tool.
        let a = raster(&[9.0, 1.0, 1.0], 1, 3);
        let b = raster(&[1.0, 9.0, 1.0], 1, 3);
        let c = raster(&[1.0, 1.0, 9.0], 1, 3);
        let (out, _) = run(json!({
            "inputs": format!("{a},{b},{c}"),
            "statistic": "highest_position",
        }));
        assert_eq!(out.get(0, 0, 0), 1.0);
        assert_eq!(out.get(0, 0, 1), 2.0);
        assert_eq!(out.get(0, 0, 2), 3.0);
    }

    #[test]
    fn lowest_position_and_ties_resolve_to_the_lowest_index() {
        let a = raster(&[5.0, 2.0], 1, 2);
        let b = raster(&[5.0, 7.0], 1, 2);
        let (out, _) = run(json!({
            "inputs": format!("{a},{b}"),
            "statistic": "lowest_position",
        }));
        // Cell 0 is a tie (5 vs 5) → layer 1; cell 1 min is layer 1 (2 < 7).
        assert_eq!(out.get(0, 0, 0), 1.0);
        assert_eq!(out.get(0, 0, 1), 1.0);
    }

    #[test]
    fn highest_position_ties_also_resolve_to_the_lowest_index() {
        let a = raster(&[4.0], 1, 1);
        let b = raster(&[4.0], 1, 1);
        let (out, _) = run(json!({
            "inputs": format!("{a},{b}"),
            "statistic": "highest_position",
        }));
        assert_eq!(out.get(0, 0, 0), 1.0);
    }

    #[test]
    fn rank_returns_the_nth_smallest_value() {
        let a = raster(&[30.0], 1, 1);
        let b = raster(&[10.0], 1, 1);
        let c = raster(&[20.0], 1, 1);
        let inputs = format!("{a},{b},{c}");
        for (n, expect) in [(1.0, 10.0), (2.0, 20.0), (3.0, 30.0)] {
            let (out, _) = run(json!({
                "inputs": inputs, "statistic": "rank", "selector": n,
            }));
            assert_eq!(out.get(0, 0, 0), expect, "rank {n}");
        }
    }

    #[test]
    fn rank_outside_the_layer_count_is_nodata() {
        let a = raster(&[1.0], 1, 1);
        let b = raster(&[2.0], 1, 1);
        let (out, _) = run(json!({
            "inputs": format!("{a},{b}"), "statistic": "rank", "selector": 5,
        }));
        assert_eq!(out.get(0, 0, 0), out.nodata);
    }

    #[test]
    fn popularity_returns_the_value_with_exactly_n_occurrences() {
        // Layers: 7,7,7,3 → the value occurring 3 times is 7; occurring once, 3.
        let a = raster(&[7.0], 1, 1);
        let b = raster(&[7.0], 1, 1);
        let c = raster(&[7.0], 1, 1);
        let d = raster(&[3.0], 1, 1);
        let inputs = format!("{a},{b},{c},{d}");
        let (out3, _) = run(json!({
            "inputs": inputs, "statistic": "popularity", "selector": 3,
        }));
        assert_eq!(out3.get(0, 0, 0), 7.0);
        let (out1, _) = run(json!({
            "inputs": inputs, "statistic": "popularity", "selector": 1,
        }));
        assert_eq!(out1.get(0, 0, 0), 3.0);
    }

    #[test]
    fn popularity_is_nodata_when_two_values_share_the_count() {
        // 5,5,8,8 — both values occur twice, so "the" value is ambiguous.
        let a = raster(&[5.0], 1, 1);
        let b = raster(&[5.0], 1, 1);
        let c = raster(&[8.0], 1, 1);
        let d = raster(&[8.0], 1, 1);
        let (out, _) = run(json!({
            "inputs": format!("{a},{b},{c},{d}"),
            "statistic": "popularity",
            "selector": 2,
        }));
        assert_eq!(out.get(0, 0, 0), out.nodata);
    }

    #[test]
    fn popularity_is_nodata_when_no_value_has_that_count() {
        let a = raster(&[1.0], 1, 1);
        let b = raster(&[2.0], 1, 1);
        let (out, _) = run(json!({
            "inputs": format!("{a},{b}"), "statistic": "popularity", "selector": 2,
        }));
        assert_eq!(out.get(0, 0, 0), out.nodata);
    }

    #[test]
    fn selector_may_be_a_raster_varying_per_cell() {
        let a = raster(&[10.0, 10.0], 1, 2);
        let b = raster(&[20.0, 20.0], 1, 2);
        let c = raster(&[30.0, 30.0], 1, 2);
        let sel = raster(&[1.0, 3.0], 1, 2);
        let (out, _) = run(json!({
            "inputs": format!("{a},{b},{c}"),
            "statistic": "rank",
            "selector": sel,
        }));
        assert_eq!(out.get(0, 0, 0), 10.0);
        assert_eq!(out.get(0, 0, 1), 30.0);
    }

    #[test]
    fn nodata_layers_are_skipped_and_shift_the_reported_position() {
        // Layer 1 is no-data at the cell, so the surviving layers are 2 and 3
        // and the max among them is reported by its position within the
        // *valid* observations — the documented ignore_nodata behaviour.
        let a = raster(&[-9999.0], 1, 1);
        let b = raster(&[5.0], 1, 1);
        let c = raster(&[9.0], 1, 1);
        let (out, _) = run(json!({
            "inputs": format!("{a},{b},{c}"),
            "statistic": "highest_position",
        }));
        assert_eq!(out.get(0, 0, 0), 2.0);
    }

    #[test]
    fn ignore_nodata_false_makes_the_whole_cell_nodata() {
        let a = raster(&[-9999.0], 1, 1);
        let b = raster(&[5.0], 1, 1);
        let (out, _) = run(json!({
            "inputs": format!("{a},{b}"),
            "statistic": "highest_position",
            "ignore_nodata": false,
        }));
        assert_eq!(out.get(0, 0, 0), out.nodata);
    }

    #[test]
    fn multi_band_policy_produces_one_output_band_per_input_band() {
        // Two 2-band inputs; band 1 peaks in input 1, band 2 peaks in input 2.
        let a = multiband(1, 1, vec![vec![9.0], vec![1.0]]);
        let b = multiband(1, 1, vec![vec![1.0], vec![9.0]]);
        let (out, res) = run(json!({
            "inputs": format!("{a},{b}"),
            "statistic": "highest_position",
            "process_as_multiband": "multi_band",
        }));
        assert_eq!(res.outputs["bands"], json!(2));
        assert_eq!(out.get(0, 0, 0), 1.0);
        assert_eq!(out.get(1, 0, 0), 2.0);
    }

    #[test]
    fn a_numeric_selector_passes_validation() {
        // Regression: validate() tested Value::as_str, so a JSON number was
        // rejected even though run() accepts it through opt_f64.
        let a = raster(&[1.0], 1, 1);
        let b = raster(&[2.0], 1, 1);
        let args: ToolArgs = serde_json::from_value(json!({
            "inputs": format!("{a},{b}"), "statistic": "rank", "selector": 1,
        }))
        .unwrap();
        assert!(CellPositionStatisticsTool.validate(&args).is_ok());
    }

    #[test]
    fn rejects_missing_selector_and_bad_statistic() {
        let a = raster(&[1.0], 1, 1);
        let b = raster(&[2.0], 1, 1);
        let inputs = format!("{a},{b}");

        let args: ToolArgs =
            serde_json::from_value(json!({"inputs": inputs, "statistic": "rank"})).unwrap();
        assert!(CellPositionStatisticsTool.validate(&args).is_err());

        let args: ToolArgs =
            serde_json::from_value(json!({"inputs": inputs, "statistic": "nope"})).unwrap();
        assert!(CellPositionStatisticsTool.validate(&args).is_err());
    }
}
