//! GeoLibre tool: slice a raster cube down to a subset of its dimension.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Subset Multidimensional Raster*
//! (Multidimension), which also covers *Select By Dimension*.
//!
//! ## Why the cube tools needed this
//!
//! Round 18 added `cube.rs` and four cube consumers —
//! `aggregate_multidimensional_raster`, `dimensional_moving_statistics`,
//! `multidimensional_raster_correlation`, `multidimensional_principal_components`
//! — and every one of them takes a **whole** cube. There was no way to narrow one
//! first, so a study window, a season, or a single decade could not be analysed
//! without regenerating the input outside the catalog, and every experiment paid
//! for the full stack (which for the PCA and correlation tools is the expensive
//! part).
//!
//! `slice_raster` is unrelated despite the name: it reclassifies cell *values*
//! into slices. Nothing subsetted the dimension axis.
//!
//! ## Scope note: no variable filtering
//!
//! ArcGIS cubes carry named variables. This stack's cube model
//! (`cube.rs`, established in round 16 by `multidimensional_anomaly`) is an
//! ordered stack of co-registered slices with optional numeric coordinates and
//! no variable metadata — there is nowhere to read variable names from. Filtering
//! is therefore by dimension coordinate, index, and step. Where a caller needs
//! per-variable cubes, they are separate inputs in this model already.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::args_common::usize_or;
use crate::common::{parse_optional_output, write_or_store_output};
use crate::cube::load_cube;
use crate::raster_stack::raster_like_multiband;

pub struct SubsetMultidimensionalRasterTool;

impl Tool for SubsetMultidimensionalRasterTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "subset_multidimensional_raster",
            display_name: "Subset Multidimensional Raster",
            summary: "Narrows a multidimensional raster cube to a subset of its slices, by dimension-coordinate range, explicit indices, or a step (ArcGIS Subset Multidimensional Raster / Select By Dimension). The round-18 cube tools all consume a whole cube, so there was no way to analyse a study window without regenerating the input outside the catalog.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Multidimensional raster, or a comma-separated list of co-registered rasters whose bands concatenate into one cube.",
                    required: true,
                },
                ToolParamSpec {
                    name: "dimension_values",
                    description: "Optional strictly increasing per-slice coordinates (dates, depths, wavelengths), one per slice. Required by 'dimension_range'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "dimension",
                    description: "Name of the dimension, used only in reporting (default 'slice').",
                    required: false,
                },
                ToolParamSpec {
                    name: "dimension_range",
                    description: "Inclusive coordinate range as 'min,max'. Requires 'dimension_values'. Mutually exclusive with 'indices'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "indices",
                    description: "Explicit 0-based slice indices, e.g. '0,2,5'. Mutually exclusive with 'dimension_range'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "step",
                    description: "Keep every Nth slice after the other filters (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output cube raster, one band per retained slice. If omitted, stored in memory.",
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
        let has_range = str_param(args, "dimension_range").is_some();
        let has_indices = str_param(args, "indices").is_some();
        // Silently preferring one would mean something completely different
        // from what the caller asked for.
        if has_range && has_indices {
            return Err(ToolError::Validation(
                "'dimension_range' and 'indices' are mutually exclusive: one selects by \
                 coordinate, the other by position"
                    .to_string(),
            ));
        }
        if has_range && str_param(args, "dimension_values").is_none() {
            return Err(ToolError::Validation(
                "'dimension_range' selects by coordinate, so 'dimension_values' is required; \
                 without it there are no coordinates to compare against (use 'indices' to select \
                 by position instead)"
                    .to_string(),
            ));
        }
        if let Some(r) = str_param(args, "dimension_range") {
            parse_range(r)?;
        }
        if usize_or(args, "step", 1)? == 0 {
            return Err(ToolError::Validation(
                "'step' must be at least 1".to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let step = usize_or(args, "step", 1)?.max(1);
        let cube = load_cube(
            args,
            "input",
            Some("dimension_values"),
            Some("dimension"),
            1,
        )?;
        let n = cube.len();

        let mut keep: Vec<usize> = (0..n).collect();

        if let Some(r) = str_param(args, "dimension_range") {
            let (lo, hi) = parse_range(r)?;
            if cube.coords.is_none() {
                // Falling back to index semantics here would quietly answer a
                // different question than the one asked.
                return Err(ToolError::Validation(
                    "'dimension_range' needs 'dimension_values'; the cube has no coordinates"
                        .to_string(),
                ));
            }
            keep.retain(|&i| {
                let c = cube.coord(i);
                c >= lo && c <= hi
            });
        } else if let Some(spec) = str_param(args, "indices") {
            let requested = parse_indices(spec)?;
            for &i in &requested {
                if i >= n {
                    return Err(ToolError::Validation(format!(
                        "index {i} is out of range; the cube has {n} slice(s) (0-{})",
                        n - 1
                    )));
                }
            }
            keep = requested;
        }

        if step > 1 {
            keep = keep.into_iter().step_by(step).collect();
        }

        if keep.is_empty() {
            return Err(ToolError::Execution(
                "the subset selected no slices; widen 'dimension_range', check 'indices', or \
                 lower 'step'"
                    .to_string(),
            ));
        }

        let template = cube.template();
        let nodata = template.nodata;
        let (rows, cols) = (cube.rows, cube.cols);
        ctx.progress.info(&format!(
            "{n} slice(s) -> {} along '{}'",
            keep.len(),
            cube.dimension
        ));

        let mut bands: Vec<Vec<f64>> = Vec::with_capacity(keep.len());
        for &s in &keep {
            let mut buf = vec![nodata; rows * cols];
            for r in 0..rows {
                for c in 0..cols {
                    if let Some(v) = cube.get(s, r, c) {
                        buf[r * cols + c] = v;
                    }
                }
            }
            bands.push(buf);
        }

        let kept_coords: Vec<f64> = keep.iter().map(|&i| cube.coord(i)).collect();
        let out = raster_like_multiband(template, &bands, nodata, DataType::F32)?;
        let out_path = write_or_store_output(out, parse_optional_output(args, "output")?)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_slices".to_string(), json!(n));
        outputs.insert("output_slices".to_string(), json!(keep.len()));
        outputs.insert("dimension".to_string(), json!(cube.dimension));
        outputs.insert("indices".to_string(), json!(keep));
        outputs.insert("coord_min".to_string(), json!(kept_coords[0]));
        outputs.insert(
            "coord_max".to_string(),
            json!(kept_coords[kept_coords.len() - 1]),
        );
        outputs.insert(
            "has_coordinates".to_string(),
            json!(cube.coords.is_some()),
        );
        Ok(ToolRunResult { outputs })
    }
}

fn str_param<'a>(args: &'a ToolArgs, key: &str) -> Option<&'a str> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
}

fn parse_range(s: &str) -> Result<(f64, f64), ToolError> {
    let parts: Vec<f64> = s
        .split([',', ';'])
        .map(str::trim)
        .filter(|p| !p.is_empty())
        .map(|p| {
            p.parse::<f64>().map_err(|_| {
                ToolError::Validation(format!("'dimension_range' entry '{p}' is not a number"))
            })
        })
        .collect::<Result<_, _>>()?;
    if parts.len() != 2 {
        return Err(ToolError::Validation(
            "'dimension_range' must be 'min,max'".to_string(),
        ));
    }
    if !(parts[0] <= parts[1]) {
        return Err(ToolError::Validation(format!(
            "'dimension_range' min ({}) must not exceed max ({})",
            parts[0], parts[1]
        )));
    }
    Ok((parts[0], parts[1]))
}

fn parse_indices(s: &str) -> Result<Vec<usize>, ToolError> {
    let mut v: Vec<usize> = s
        .split([',', ';'])
        .map(str::trim)
        .filter(|p| !p.is_empty())
        .map(|p| {
            p.parse::<usize>()
                .map_err(|_| ToolError::Validation(format!("'indices' entry '{p}' is not a non-negative integer")))
        })
        .collect::<Result<_, _>>()?;
    if v.is_empty() {
        return Err(ToolError::Validation(
            "'indices' listed no slices".to_string(),
        ));
    }
    v.sort_unstable();
    v.dedup();
    Ok(v)
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

    /// A 1x1 cube of five slices valued 10, 20, 30, 40, 50.
    fn five() -> String {
        cube_raster(
            1,
            1,
            &[vec![10.0], vec![20.0], vec![30.0], vec![40.0], vec![50.0]],
        )
    }

    fn run(args: Value) -> (Raster, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = SubsetMultidimensionalRasterTool.run(&args, &ctx()).unwrap();
        let out = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        (out, res)
    }

    fn values(r: &Raster) -> Vec<f64> {
        (0..r.bands).map(|b| r.get(b as isize, 0, 0)).collect()
    }

    #[test]
    fn without_filters_the_cube_passes_through_unchanged() {
        let (out, res) = run(json!({"input": five()}));
        assert_eq!(out.bands, 5);
        assert_eq!(values(&out), vec![10.0, 20.0, 30.0, 40.0, 50.0]);
        assert_eq!(res.outputs["output_slices"], json!(5));
    }

    #[test]
    fn a_dimension_range_keeps_the_slices_inside_it_inclusively() {
        let (out, res) = run(json!({
            "input": five(), "dimension_values": "2000,2001,2002,2003,2004",
            "dimension_range": "2001,2003",
        }));
        assert_eq!(values(&out), vec![20.0, 30.0, 40.0]);
        assert_eq!(res.outputs["coord_min"], json!(2001.0));
        assert_eq!(res.outputs["coord_max"], json!(2003.0));
        assert_eq!(res.outputs["indices"], json!([1, 2, 3]));
    }

    #[test]
    fn a_range_matching_a_coordinate_exactly_includes_it() {
        // Inclusive bounds: an exclusive comparison would silently drop the
        // endpoint year a user explicitly asked for.
        let (out, _) = run(json!({
            "input": five(), "dimension_values": "2000,2001,2002,2003,2004",
            "dimension_range": "2004,2004",
        }));
        assert_eq!(values(&out), vec![50.0]);
    }

    #[test]
    fn explicit_indices_select_by_position() {
        let (out, _) = run(json!({"input": five(), "indices": "0,2,4"}));
        assert_eq!(values(&out), vec![10.0, 30.0, 50.0]);
    }

    #[test]
    fn indices_are_sorted_and_deduplicated() {
        let (out, res) = run(json!({"input": five(), "indices": "4,0,4,2"}));
        assert_eq!(values(&out), vec![10.0, 30.0, 50.0]);
        assert_eq!(res.outputs["indices"], json!([0, 2, 4]));
    }

    #[test]
    fn step_thins_the_retained_slices() {
        let (out, _) = run(json!({"input": five(), "step": 2}));
        assert_eq!(values(&out), vec![10.0, 30.0, 50.0]);
    }

    #[test]
    fn step_applies_after_the_range_not_before() {
        // Range first (slices 1..3), then every 2nd of those.
        let (out, _) = run(json!({
            "input": five(), "dimension_values": "2000,2001,2002,2003,2004",
            "dimension_range": "2001,2003", "step": 2,
        }));
        assert_eq!(values(&out), vec![20.0, 40.0]);
    }

    #[test]
    fn no_data_survives_the_subset() {
        let path = cube_raster(2, 1, &[vec![1.0, -9999.0], vec![3.0, 4.0]]);
        let (out, _) = run(json!({"input": path, "indices": "0"}));
        assert_eq!(out.get(0, 0, 0), 1.0);
        assert_eq!(out.get(0, 0, 1), out.nodata);
    }

    #[test]
    fn an_out_of_range_index_is_reported_not_clamped() {
        let args: ToolArgs =
            serde_json::from_value(json!({"input": five(), "indices": "0,9"})).unwrap();
        let err = SubsetMultidimensionalRasterTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("out of range"), "{err}");
    }

    #[test]
    fn a_range_selecting_nothing_errors_rather_than_writing_an_empty_cube() {
        let args: ToolArgs = serde_json::from_value(json!({
            "input": five(), "dimension_values": "2000,2001,2002,2003,2004",
            "dimension_range": "3000,3001",
        }))
        .unwrap();
        let err = SubsetMultidimensionalRasterTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("no slices"), "{err}");
    }

    #[test]
    fn a_range_without_coordinates_is_refused_rather_than_read_as_indices() {
        // Treating 2001..2003 as slice positions would answer a completely
        // different question without saying so.
        let args: ToolArgs = serde_json::from_value(json!({
            "input": five(), "dimension_range": "2001,2003",
        }))
        .unwrap();
        assert!(SubsetMultidimensionalRasterTool.validate(&args).is_err());
    }

    #[test]
    fn subsetting_then_aggregating_matches_aggregating_the_subset_directly() {
        // The point of the tool: the narrowed cube must behave exactly like a
        // cube that only ever held those slices.
        let full = five();
        let (narrowed, _) = run(json!({"input": full, "indices": "1,2,3"}));
        let direct = cube_raster(1, 1, &[vec![20.0], vec![30.0], vec![40.0]]);
        let narrowed_id = wbraster::memory_store::put_raster(narrowed);
        let narrowed_path = wbraster::memory_store::make_raster_memory_path(&narrowed_id);

        let agg = |p: &str| -> f64 {
            let args: ToolArgs = serde_json::from_value(json!({
                "input": p, "aggregation_method": "mean",
            }))
            .unwrap();
            let res = crate::aggregate_multidimensional_raster::AggregateMultidimensionalRasterTool
                .run(&args, &ctx())
                .unwrap();
            let r = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
            r.get(0, 0, 0)
        };
        assert!((agg(&narrowed_path) - agg(&direct)).abs() < 1e-9);
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SubsetMultidimensionalRasterTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"input": "a.tif", "step": 0})));
        assert!(bad(json!({"input": "a.tif", "dimension_range": "5"})));
        assert!(bad(json!({"input": "a.tif", "dimension_range": "9,1", "dimension_values": "1,2"})));
        // Mutually exclusive selectors.
        assert!(bad(json!({
            "input": "a.tif", "indices": "0", "dimension_range": "1,2",
            "dimension_values": "1,2",
        })));
    }
}
