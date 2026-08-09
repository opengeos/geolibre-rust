//! GeoLibre tool: per-cell count of how many rasters in a stack compare a given
//! way against a reference raster.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Equal To Frequency*, *Greater Than
//! Frequency* and *Less Than Frequency* (Spatial Analyst / Image Analyst Local
//! toolset), consolidated behind one `comparison` parameter and extended with
//! the three inclusive/negated operators that cost nothing to add.
//!
//! ## Why neither registry covers this
//!
//! This is a *comparative count against a reference surface*, and nothing in
//! the catalog performs one. `cell_statistics` reduces a stack but has no
//! reference input. `raster_calculator` can express a single comparison but not
//! an aggregation across an arbitrary-length stack. `conditional_evaluation` is
//! a one-condition selector, not a counter.
//!
//! The question it answers is the standard exceedance-frequency question: how
//! many of these 30 daily grids exceeded the threshold surface, how many of
//! these 40 model years came in below the baseline. That is a routine
//! climate/hydrology need with no current answer in the catalog.
//!
//! ## Floating-point equality
//!
//! Exact `==` on floating-point rasters is almost always the wrong test, so
//! `equal`/`not_equal` take an absolute `tolerance` (default 0, preserving
//! exact behaviour for integer-valued data). This is deliberate: without it the
//! `equal` mode would silently return zeros on any resampled float input.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::args_common::{bool_or, choice_or, f64_or, req_str};
use crate::common::{load_input_raster, parse_optional_output};
use crate::raster_stack::{
    check_alignment, load_stack, parse_band_policy, parse_input_paths, write_stack_result,
};

pub struct FrequencyComparisonTool;

impl Tool for FrequencyComparisonTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "frequency_comparison",
            display_name: "Frequency Comparison",
            summary: "Counts, per cell, how many rasters in a stack are equal to / greater than / less than the corresponding cell of a reference value raster (ArcGIS Equal To Frequency, Greater Than Frequency and Less Than Frequency, plus the inclusive and negated operators). Neither registry can aggregate a comparison across a stack: cell_statistics has no reference input, raster_calculator expresses one comparison but cannot count across an arbitrary-length stack, and conditional_evaluation is a single-condition selector.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "value_raster",
                    description: "Reference raster each stack layer is compared against, co-registered with the inputs.",
                    required: true,
                },
                ToolParamSpec {
                    name: "inputs",
                    description: "One multiband raster (each band is a layer) or a comma-separated list of co-registered rasters.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output count raster (0..n). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "comparison",
                    description: "'equal' (default), 'greater', 'less', 'greater_equal', 'less_equal', or 'not_equal'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Absolute tolerance for 'equal'/'not_equal' on floating-point data (default 0 = exact).",
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
        req_str(args, "value_raster")?;
        parse_input_paths(args, "inputs")?;
        parse_comparison(args)?;
        parse_band_policy(args, "process_as_multiband")?;
        let tol = f64_or(args, "tolerance", 0.0)?;
        if tol < 0.0 {
            return Err(ToolError::Validation(
                "'tolerance' must be >= 0".to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let value_path = req_str(args, "value_raster")?.to_string();
        let paths = parse_input_paths(args, "inputs")?;
        let comparison = parse_comparison(args)?;
        let tolerance = f64_or(args, "tolerance", 0.0)?;
        let ignore_nodata = bool_or(args, "ignore_nodata", true)?;
        let policy = parse_band_policy(args, "process_as_multiband")?;
        let output = parse_optional_output(args, "output")?;

        // One layer is a meaningful comparison here (unlike a stack reducer),
        // so the floor is 1 rather than cell_statistics' 2.
        let stack = load_stack(&paths, policy, 1)?;
        let value = load_input_raster(&value_path)?;
        // The reference is read cell-by-cell against the stack grid, so a
        // mismatched grid would silently compare the wrong locations.
        check_alignment(&[stack.template().clone(), value.clone()])?;

        let (rows, cols) = (stack.rows, stack.cols);
        ctx.progress.info(&format!(
            "{} input(s), {} layer(s), {rows}x{cols}, comparison {}, {}",
            paths.len(),
            stack.total_layers(),
            comparison.label(),
            policy.label()
        ));

        let nodata = -9999.0_f64;
        let n_groups = stack.group_count();
        let mut bands: Vec<Vec<f64>> = Vec::with_capacity(n_groups);
        let mut vals: Vec<f64> = Vec::new();

        for g in 0..n_groups {
            let mut out = vec![nodata; rows * cols];
            for r in 0..rows {
                for c in 0..cols {
                    let reference = value.get(0, r as isize, c as isize);
                    if reference == value.nodata || !reference.is_finite() {
                        continue;
                    }
                    let had_nodata = stack.cell_values(g, r, c, &mut vals);
                    if vals.is_empty() || (!ignore_nodata && had_nodata) {
                        continue;
                    }
                    let count = vals
                        .iter()
                        .filter(|v| comparison.holds(**v, reference, tolerance))
                        .count();
                    out[r * cols + c] = count as f64;
                }
                ctx.progress
                    .progress((g as f64 + (r as f64 + 1.0) / rows as f64) / n_groups as f64);
            }
            bands.push(out);
        }

        let out_path = write_stack_result(stack.template(), bands, nodata, DataType::I32, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("comparison".to_string(), json!(comparison.label()));
        outputs.insert("layers".to_string(), json!(stack.total_layers()));
        outputs.insert("bands".to_string(), json!(n_groups));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Comparison {
    Equal,
    NotEqual,
    Greater,
    GreaterEqual,
    Less,
    LessEqual,
}

impl Comparison {
    fn label(self) -> &'static str {
        match self {
            Comparison::Equal => "equal",
            Comparison::NotEqual => "not_equal",
            Comparison::Greater => "greater",
            Comparison::GreaterEqual => "greater_equal",
            Comparison::Less => "less",
            Comparison::LessEqual => "less_equal",
        }
    }

    /// Does `value` compare this way against `reference`?
    ///
    /// `tolerance` only affects the equality operators; using it to fuzz the
    /// ordering comparisons would make `greater` and `less_equal` overlap.
    fn holds(self, value: f64, reference: f64, tolerance: f64) -> bool {
        match self {
            Comparison::Equal => (value - reference).abs() <= tolerance,
            Comparison::NotEqual => (value - reference).abs() > tolerance,
            Comparison::Greater => value > reference,
            Comparison::GreaterEqual => value >= reference,
            Comparison::Less => value < reference,
            Comparison::LessEqual => value <= reference,
        }
    }
}

fn parse_comparison(args: &ToolArgs) -> Result<Comparison, ToolError> {
    let choice = choice_or(
        args,
        "comparison",
        &[
            "equal",
            "not_equal",
            "greater",
            "greater_equal",
            "less",
            "less_equal",
        ],
        "equal",
    )?;
    Ok(match choice {
        "equal" => Comparison::Equal,
        "not_equal" => Comparison::NotEqual,
        "greater" => Comparison::Greater,
        "greater_equal" => Comparison::GreaterEqual,
        "less" => Comparison::Less,
        _ => Comparison::LessEqual,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
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

    fn raster(data: &[f64], rows: usize, cols: usize) -> String {
        multiband(rows, cols, vec![data.to_vec()])
    }

    fn run(args: Value) -> Raster {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = FrequencyComparisonTool.run(&args, &ctx()).unwrap();
        load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    #[test]
    fn counts_layers_above_the_reference_surface() {
        // Reference 5. Stack 3, 6, 9 → one below, two above.
        let value = raster(&[5.0], 1, 1);
        let stack = multiband(1, 1, vec![vec![3.0], vec![6.0], vec![9.0]]);
        let out = run(json!({
            "value_raster": value, "inputs": stack, "comparison": "greater",
        }));
        assert_eq!(out.get(0, 0, 0), 2.0);
        let out = run(json!({
            "value_raster": value, "inputs": stack, "comparison": "less",
        }));
        assert_eq!(out.get(0, 0, 0), 1.0);
    }

    #[test]
    fn the_reference_varies_per_cell() {
        // Two cells with different thresholds — the whole reason the reference
        // is a raster and not a scalar.
        let value = raster(&[2.0, 8.0], 1, 2);
        let stack = multiband(1, 2, vec![vec![5.0, 5.0], vec![9.0, 9.0]]);
        let out = run(json!({
            "value_raster": value, "inputs": stack, "comparison": "greater",
        }));
        assert_eq!(out.get(0, 0, 0), 2.0); // 5 and 9 both exceed 2
        assert_eq!(out.get(0, 0, 1), 1.0); // only 9 exceeds 8
    }

    #[test]
    fn equal_is_exact_by_default_and_tolerant_on_request() {
        let value = raster(&[1.0], 1, 1);
        let stack = multiband(1, 1, vec![vec![1.0], vec![1.05], vec![2.0]]);
        let exact = run(json!({
            "value_raster": value, "inputs": stack, "comparison": "equal",
        }));
        assert_eq!(exact.get(0, 0, 0), 1.0);
        let tolerant = run(json!({
            "value_raster": value, "inputs": stack, "comparison": "equal", "tolerance": 0.1,
        }));
        assert_eq!(tolerant.get(0, 0, 0), 2.0);
    }

    #[test]
    fn inclusive_and_negated_operators_are_consistent() {
        // greater + less_equal must partition the stack, and equal + not_equal
        // likewise — a property that catches an operator wired to the wrong arm.
        let value = raster(&[5.0], 1, 1);
        let stack = multiband(1, 1, vec![vec![4.0], vec![5.0], vec![6.0]]);
        let gt = run(json!({"value_raster": value, "inputs": stack, "comparison": "greater"}));
        let le = run(json!({"value_raster": value, "inputs": stack, "comparison": "less_equal"}));
        assert_eq!(gt.get(0, 0, 0) + le.get(0, 0, 0), 3.0);
        let eq = run(json!({"value_raster": value, "inputs": stack, "comparison": "equal"}));
        let ne = run(json!({"value_raster": value, "inputs": stack, "comparison": "not_equal"}));
        assert_eq!(eq.get(0, 0, 0) + ne.get(0, 0, 0), 3.0);
        assert_eq!(eq.get(0, 0, 0), 1.0);
    }

    #[test]
    fn nodata_in_the_reference_makes_the_cell_nodata() {
        let value = raster(&[-9999.0, 5.0], 1, 2);
        let stack = multiband(1, 2, vec![vec![9.0, 9.0], vec![9.0, 9.0]]);
        let out = run(json!({
            "value_raster": value, "inputs": stack, "comparison": "greater",
        }));
        assert_eq!(out.get(0, 0, 0), out.nodata);
        assert_eq!(out.get(0, 0, 1), 2.0);
    }

    #[test]
    fn nodata_layers_are_skipped_unless_ignore_nodata_is_false() {
        let value = raster(&[5.0], 1, 1);
        let stack = multiband(1, 1, vec![vec![-9999.0], vec![9.0]]);
        let skipped = run(json!({
            "value_raster": value, "inputs": stack, "comparison": "greater",
        }));
        assert_eq!(skipped.get(0, 0, 0), 1.0);
        let strict = run(json!({
            "value_raster": value, "inputs": stack, "comparison": "greater",
            "ignore_nodata": false,
        }));
        assert_eq!(strict.get(0, 0, 0), strict.nodata);
    }

    #[test]
    fn multi_band_policy_counts_each_band_separately() {
        let value = multiband(1, 1, vec![vec![5.0], vec![5.0]]);
        let a = multiband(1, 1, vec![vec![9.0], vec![1.0]]);
        let b = multiband(1, 1, vec![vec![9.0], vec![1.0]]);
        let args: ToolArgs = serde_json::from_value(json!({
            "value_raster": value,
            "inputs": format!("{a},{b}"),
            "comparison": "greater",
            "process_as_multiband": "multi_band",
        }))
        .unwrap();
        let res = FrequencyComparisonTool.run(&args, &ctx()).unwrap();
        assert_eq!(res.outputs["bands"], json!(2));
        let out = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        assert_eq!(out.get(0, 0, 0), 2.0); // band 1: both 9 > 5
        assert_eq!(out.get(1, 0, 0), 0.0); // band 2: neither 1 > 5
    }

    #[test]
    fn rejects_bad_parameters() {
        let value = raster(&[1.0], 1, 1);
        let stack = raster(&[1.0], 1, 1);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            FrequencyComparisonTool.validate(&args).is_err()
        };
        assert!(bad(json!({"inputs": stack})));
        assert!(bad(
            json!({"value_raster": value, "inputs": stack, "comparison": "nope"})
        ));
        assert!(bad(
            json!({"value_raster": value, "inputs": stack, "tolerance": -1})
        ));
    }

    #[test]
    fn rejects_a_misaligned_reference_raster() {
        let value = raster(&[1.0, 2.0], 1, 2);
        let stack = multiband(1, 1, vec![vec![1.0], vec![2.0]]);
        let args: ToolArgs = serde_json::from_value(json!({
            "value_raster": value, "inputs": stack,
        }))
        .unwrap();
        assert!(FrequencyComparisonTool.run(&args, &ctx()).is_err());
    }
}
