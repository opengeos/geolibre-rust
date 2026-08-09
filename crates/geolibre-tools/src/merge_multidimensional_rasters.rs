//! GeoLibre tool: concatenate raster cubes along a shared dimension.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Merge Multidimensional Rasters*
//! (Multidimension), and the complement of `subset_multidimensional_raster`.
//!
//! ## Why this is the other half
//!
//! The round-18 cube tools all **reduce** a cube
//! (`aggregate_multidimensional_raster`) or **analyse** one
//! (`multidimensional_raster_correlation`, `multidimensional_principal_components`,
//! `dimensional_moving_statistics`). Nothing built a bigger one.
//!
//! Multi-year and multi-sensor archives arrive as separate files per period —
//! one cube per year, per tile, per acquisition campaign — and every analysis
//! tool wants the single joined cube. `composite_bands` stacks single-band
//! rasters into *bands*, which is a different axis and carries no dimension
//! coordinates.
//!
//! ## Why `cube::load_cube` is not enough on its own
//!
//! `load_cube` already concatenates several inputs, but it requires the
//! coordinates to be **strictly increasing** — by design, because binning and
//! lagged correlation assume dimension order. Merging is exactly the case where
//! that does not hold yet: two archives overlap, or arrive out of order. This
//! tool is the step that sorts and de-duplicates them into a cube `load_cube`
//! will then accept.
//!
//! ## Overlap resolution
//!
//! Duplicate coordinates are resolved per cell by `resolve_overlap`. Coordinates
//! are compared with a **relative tolerance**, not `==`: float coordinates read
//! from different files are not bit-identical, and exact comparison would leave
//! near-duplicate slices sitting next to each other and break the strictly
//! increasing invariant downstream.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::args_common::choice_or;
use crate::common::{load_input_raster, parse_optional_output, write_or_store_output};
use crate::raster_stack::{check_alignment_refs, raster_like_multiband};

const RESOLVE: [&str; 6] = ["first", "last", "mean", "min", "max", "error"];

pub struct MergeMultidimensionalRastersTool;

impl Tool for MergeMultidimensionalRastersTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "merge_multidimensional_rasters",
            display_name: "Merge Multidimensional Rasters",
            summary: "Joins several co-registered raster cubes into one, ordered by dimension coordinate, resolving duplicate coordinates per cell (ArcGIS Merge Multidimensional Rasters). aggregate_multidimensional_raster reduces a cube and composite_bands stacks single-band rasters into bands, but nothing joined cubes along their shared dimension.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "inputs",
                    description: "Comma- or semicolon-separated cube paths (multiband rasters). At least two.",
                    required: true,
                },
                ToolParamSpec {
                    name: "dimension_values",
                    description: "Optional per-slice coordinates covering every slice of every input, in input order. Need not be sorted or unique — that is what this tool resolves. Without them the inputs are concatenated in the order given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "dimension",
                    description: "Name of the dimension, used only in reporting (default 'slice').",
                    required: false,
                },
                ToolParamSpec {
                    name: "resolve_overlap",
                    description: "Duplicate-coordinate rule: 'first' (default), 'last', 'mean', 'min', 'max', or 'error'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Relative tolerance for treating two coordinates as equal (default 1e-9).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output cube raster, one band per merged slice, in dimension order. If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        let inputs = parse_inputs(args)?;
        if inputs.len() < 2 {
            return Err(ToolError::Validation(format!(
                "'inputs' must list at least 2 cubes to merge, got {}",
                inputs.len()
            )));
        }
        choice_or(args, "resolve_overlap", &RESOLVE, "first")?;
        if let Some(t) = args.get("tolerance").and_then(Value::as_f64) {
            if !t.is_finite() || t < 0.0 {
                return Err(ToolError::Validation(
                    "'tolerance' must be a non-negative number".to_string(),
                ));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let inputs = parse_inputs(args)?;
        if inputs.len() < 2 {
            return Err(ToolError::Validation(
                "'inputs' must list at least 2 cubes to merge".to_string(),
            ));
        }
        let resolve = choice_or(args, "resolve_overlap", &RESOLVE, "first")?;
        let tolerance = args
            .get("tolerance")
            .and_then(Value::as_f64)
            .unwrap_or(1e-9)
            .max(0.0);
        let dimension = args
            .get("dimension")
            .and_then(Value::as_str)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .unwrap_or("slice")
            .to_string();

        let rasters: Vec<Raster> = inputs
            .iter()
            .map(|p| load_input_raster(p))
            .collect::<Result<_, _>>()?;
        let refs: Vec<&Raster> = rasters.iter().collect();
        check_alignment_refs(&refs)?;

        // (raster index, band) for every slice, in input order.
        let mut slices: Vec<(usize, isize)> = Vec::new();
        for (i, r) in rasters.iter().enumerate() {
            for b in 0..r.bands {
                slices.push((i, b as isize));
            }
        }
        let total = slices.len();

        let coords = parse_coords(args, "dimension_values", total)?;
        let (rows, cols) = (rasters[0].rows, rasters[0].cols);
        let template = &rasters[0];
        let nodata = template.nodata;

        // Group slices by coordinate. Without coordinates there is nothing to
        // sort or de-duplicate by, so input order is the only defensible
        // ordering — and it is stated in the docs rather than implied.
        let groups: Vec<(f64, Vec<usize>)> = match &coords {
            None => slices
                .iter()
                .enumerate()
                .map(|(i, _)| ((i + 1) as f64, vec![i]))
                .collect(),
            Some(cs) => {
                let mut order: Vec<usize> = (0..total).collect();
                order.sort_by(|&a, &b| cs[a].total_cmp(&cs[b]));
                let mut out: Vec<(f64, Vec<usize>)> = Vec::new();
                for idx in order {
                    let c = cs[idx];
                    // Relative comparison: coordinates from different files are
                    // never bit-identical, and == would leave near-duplicates
                    // adjacent and break the strictly-increasing invariant the
                    // downstream cube loader enforces.
                    match out.last_mut() {
                        Some((prev, members)) if close(*prev, c, tolerance) => members.push(idx),
                        _ => out.push((c, vec![idx])),
                    }
                }
                out
            }
        };

        let collisions = groups.iter().filter(|(_, m)| m.len() > 1).count();
        if collisions > 0 && resolve == "error" {
            return Err(ToolError::Execution(format!(
                "{collisions} coordinate(s) appear in more than one input; set 'resolve_overlap' \
                 to first/last/mean/min/max to say how they should be combined"
            )));
        }
        ctx.progress.info(&format!(
            "{} input(s), {total} slice(s) -> {} along '{dimension}' ({collisions} overlap(s))",
            inputs.len(),
            groups.len()
        ));

        let value = |slice: usize, r: usize, c: usize| -> Option<f64> {
            let (ri, band) = slices[slice];
            let raster = &rasters[ri];
            let v = raster.get(band, r as isize, c as isize);
            (v != raster.nodata && v.is_finite()).then_some(v)
        };

        let mut bands: Vec<Vec<f64>> = Vec::with_capacity(groups.len());
        for (_, members) in &groups {
            let mut buf = vec![nodata; rows * cols];
            for r in 0..rows {
                for c in 0..cols {
                    let vals: Vec<f64> = members.iter().filter_map(|&s| value(s, r, c)).collect();
                    if vals.is_empty() {
                        continue;
                    }
                    buf[r * cols + c] = match resolve {
                        "last" => vals[vals.len() - 1],
                        "mean" => vals.iter().sum::<f64>() / vals.len() as f64,
                        "min" => vals.iter().cloned().fold(f64::INFINITY, f64::min),
                        "max" => vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
                        // "first" and "error" (which only gets here when there
                        // were no collisions) both take the earliest input.
                        _ => vals[0],
                    };
                }
            }
            bands.push(buf);
        }

        let out_coords: Vec<f64> = groups.iter().map(|(c, _)| *c).collect();
        let out = raster_like_multiband(template, &bands, nodata, DataType::F32)?;
        let out_path = write_or_store_output(out, parse_optional_output(args, "output")?)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_count".to_string(), json!(inputs.len()));
        outputs.insert("input_slices".to_string(), json!(total));
        outputs.insert("output_slices".to_string(), json!(groups.len()));
        outputs.insert("overlaps_resolved".to_string(), json!(collisions));
        outputs.insert("resolve_overlap".to_string(), json!(resolve));
        outputs.insert("dimension".to_string(), json!(dimension));
        outputs.insert("has_coordinates".to_string(), json!(coords.is_some()));
        if coords.is_some() {
            outputs.insert("dimension_values".to_string(), json!(out_coords));
        }
        Ok(ToolRunResult { outputs })
    }
}

/// Relative closeness, so the comparison behaves the same at year 2000 and at
/// depth 0.001.
fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * a.abs().max(b.abs()).max(1.0)
}

fn parse_inputs(args: &ToolArgs) -> Result<Vec<String>, ToolError> {
    match args.get("inputs") {
        Some(Value::String(s)) => Ok(s
            .split([',', ';'])
            .map(str::trim)
            .filter(|p| !p.is_empty())
            .map(str::to_string)
            .collect()),
        Some(Value::Array(items)) => items
            .iter()
            .map(|v| {
                v.as_str().map(str::to_string).ok_or_else(|| {
                    ToolError::Validation("every entry of 'inputs' must be a string".to_string())
                })
            })
            .collect(),
        Some(_) => Err(ToolError::Validation(
            "'inputs' must be a delimited string or an array of strings".to_string(),
        )),
        None => Err(ToolError::Validation(
            "missing required parameter 'inputs'".to_string(),
        )),
    }
}

/// Unlike `cube::parse_coords`, this deliberately accepts unsorted and
/// duplicated values — resolving exactly those is the tool's job.
fn parse_coords(
    args: &ToolArgs,
    key: &str,
    n_slices: usize,
) -> Result<Option<Vec<f64>>, ToolError> {
    let Some(s) = args.get(key).and_then(Value::as_str) else {
        return Ok(None);
    };
    if s.trim().is_empty() {
        return Ok(None);
    }
    let vals: Vec<f64> = s
        .split([',', ';'])
        .map(str::trim)
        .filter(|t| !t.is_empty())
        .map(|t| {
            t.parse::<f64>()
                .map_err(|_| ToolError::Validation(format!("'{key}' entry '{t}' is not a number")))
        })
        .collect::<Result<_, _>>()?;
    if vals.len() != n_slices {
        return Err(ToolError::Validation(format!(
            "'{key}' has {} value(s) but the inputs hold {n_slices} slice(s) in total; supply one \
             coordinate per slice across every input, in input order",
            vals.len()
        )));
    }
    if let Some(bad) = vals.iter().find(|v| !v.is_finite()) {
        return Err(ToolError::Validation(format!(
            "'{key}' contains a non-finite value ({bad})"
        )));
    }
    Ok(Some(vals))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cube::test_support::cube_raster;
    use wbcore::{AllowAllCapabilities, ProgressSink};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn run(args: Value) -> (Raster, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = MergeMultidimensionalRastersTool.run(&args, &ctx()).unwrap();
        let out = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        (out, res)
    }

    fn values(r: &Raster) -> Vec<f64> {
        (0..r.bands).map(|b| r.get(b as isize, 0, 0)).collect()
    }

    #[test]
    fn concatenates_two_cubes_in_input_order_without_coordinates() {
        let a = cube_raster(1, 1, &[vec![1.0], vec![2.0]]);
        let b = cube_raster(1, 1, &[vec![3.0]]);
        let (out, res) = run(json!({"inputs": format!("{a},{b}")}));
        assert_eq!(values(&out), vec![1.0, 2.0, 3.0]);
        assert_eq!(res.outputs["output_slices"], json!(3));
        assert_eq!(res.outputs["has_coordinates"], json!(false));
    }

    #[test]
    fn coordinates_put_the_slices_in_dimension_order() {
        // The later archive is listed FIRST; the merge must still sort by year.
        let late = cube_raster(1, 1, &[vec![30.0], vec![40.0]]);
        let early = cube_raster(1, 1, &[vec![10.0], vec![20.0]]);
        let (out, res) = run(json!({
            "inputs": format!("{late},{early}"),
            "dimension_values": "2002,2003,2000,2001",
        }));
        assert_eq!(values(&out), vec![10.0, 20.0, 30.0, 40.0]);
        assert_eq!(
            res.outputs["dimension_values"],
            json!([2000.0, 2001.0, 2002.0, 2003.0])
        );
    }

    #[test]
    fn overlapping_coordinates_collapse_to_one_slice() {
        let a = cube_raster(1, 1, &[vec![10.0], vec![20.0]]);
        let b = cube_raster(1, 1, &[vec![99.0], vec![30.0]]);
        let (out, res) = run(json!({
            "inputs": format!("{a},{b}"), "dimension_values": "2000,2001,2001,2002",
        }));
        assert_eq!(res.outputs["output_slices"], json!(3));
        assert_eq!(res.outputs["overlaps_resolved"], json!(1));
        // 'first' keeps the earlier input's 20.0, not 99.0.
        assert_eq!(values(&out), vec![10.0, 20.0, 30.0]);
    }

    #[test]
    fn resolve_last_prefers_the_later_input() {
        let a = cube_raster(1, 1, &[vec![20.0]]);
        let b = cube_raster(1, 1, &[vec![99.0]]);
        let (out, _) = run(json!({
            "inputs": format!("{a},{b}"), "dimension_values": "2001,2001",
            "resolve_overlap": "last",
        }));
        assert_eq!(values(&out), vec![99.0]);
    }

    #[test]
    fn resolve_mean_min_and_max_combine_the_duplicates() {
        let a = cube_raster(1, 1, &[vec![10.0]]);
        let b = cube_raster(1, 1, &[vec![30.0]]);
        let pair = format!("{a},{b}");
        for (rule, expect) in [("mean", 20.0), ("min", 10.0), ("max", 30.0)] {
            let (out, _) = run(json!({
                "inputs": pair, "dimension_values": "5,5", "resolve_overlap": rule,
            }));
            assert_eq!(values(&out), vec![expect], "rule {rule}");
        }
    }

    #[test]
    fn resolve_error_refuses_to_guess() {
        let a = cube_raster(1, 1, &[vec![10.0]]);
        let b = cube_raster(1, 1, &[vec![30.0]]);
        let args: ToolArgs = serde_json::from_value(json!({
            "inputs": format!("{a},{b}"), "dimension_values": "5,5",
            "resolve_overlap": "error",
        }))
        .unwrap();
        let err = MergeMultidimensionalRastersTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("resolve_overlap"), "{err}");
    }

    #[test]
    fn nearly_equal_coordinates_are_treated_as_one() {
        // The exact-equality trap: 2001.0 and 2001.0000000001 read from two
        // files must not survive as two adjacent slices, which would break the
        // strictly-increasing invariant cube::load_cube enforces downstream.
        let a = cube_raster(1, 1, &[vec![10.0]]);
        let b = cube_raster(1, 1, &[vec![30.0]]);
        let (out, res) = run(json!({
            "inputs": format!("{a},{b}"),
            "dimension_values": "2001.0,2001.0000000001",
        }));
        assert_eq!(res.outputs["output_slices"], json!(1));
        assert_eq!(values(&out), vec![10.0]);
    }

    #[test]
    fn a_tighter_tolerance_keeps_them_apart() {
        let a = cube_raster(1, 1, &[vec![10.0]]);
        let b = cube_raster(1, 1, &[vec![30.0]]);
        let (_, res) = run(json!({
            "inputs": format!("{a},{b}"),
            "dimension_values": "2001.0,2001.0000000001", "tolerance": 0.0,
        }));
        assert_eq!(res.outputs["output_slices"], json!(2));
    }

    #[test]
    fn merging_fills_a_gap_left_by_a_no_data_cell_in_one_input() {
        // The practical reason for 'mean'/'last': two overlapping archives
        // where one has a cloud gap.
        let a = cube_raster(2, 1, &[vec![10.0, -9999.0]]);
        let b = cube_raster(2, 1, &[vec![-9999.0, 20.0]]);
        let (out, _) = run(json!({
            "inputs": format!("{a},{b}"), "dimension_values": "5,5",
            "resolve_overlap": "mean",
        }));
        assert_eq!(out.get(0, 0, 0), 10.0);
        assert_eq!(out.get(0, 0, 1), 20.0, "the gap is filled from the other cube");
    }

    #[test]
    fn a_grid_mismatch_is_refused() {
        let a = cube_raster(1, 1, &[vec![1.0]]);
        let b = cube_raster(2, 2, &[vec![1.0, 2.0, 3.0, 4.0]]);
        let args: ToolArgs =
            serde_json::from_value(json!({"inputs": format!("{a},{b}")})).unwrap();
        assert!(MergeMultidimensionalRastersTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn the_merged_cube_is_accepted_by_the_downstream_cube_loader() {
        // The whole point: load_cube demands strictly increasing coordinates,
        // and the merge output plus its reported coordinates must satisfy it.
        let late = cube_raster(1, 1, &[vec![30.0]]);
        let early = cube_raster(1, 1, &[vec![10.0], vec![20.0]]);
        let (out, res) = run(json!({
            "inputs": format!("{late},{early}"),
            "dimension_values": "2002,2000,2001",
        }));
        let coords: Vec<String> = res.outputs["dimension_values"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap().to_string())
            .collect();
        let id = wbraster::memory_store::put_raster(out);
        let path = wbraster::memory_store::make_raster_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": path, "dimension_values": coords.join(","),
        }))
        .unwrap();
        let cube = crate::cube::load_cube(&args, "input", Some("dimension_values"), None, 1)
            .expect("merged cube must satisfy load_cube's ordering invariant");
        assert_eq!(cube.len(), 3);
        assert_eq!(cube.coord(0), 2000.0);
    }

    #[test]
    fn rejects_bad_parameters() {
        let a = cube_raster(1, 1, &[vec![1.0]]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            MergeMultidimensionalRastersTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        // A single cube is not a merge.
        assert!(bad(json!({"inputs": a})));
        assert!(bad(json!({"inputs": "a.tif,b.tif", "resolve_overlap": "average"})));
        assert!(bad(json!({"inputs": "a.tif,b.tif", "tolerance": -1})));
    }

    #[test]
    fn a_coordinate_count_mismatch_names_the_totals() {
        let a = cube_raster(1, 1, &[vec![1.0], vec![2.0]]);
        let b = cube_raster(1, 1, &[vec![3.0]]);
        let args: ToolArgs = serde_json::from_value(json!({
            "inputs": format!("{a},{b}"), "dimension_values": "1,2",
        }))
        .unwrap();
        let err = MergeMultidimensionalRastersTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("3 slice(s)"), "{err}");
    }
}
