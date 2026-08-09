//! GeoLibre tool: stack N single-band rasters into one multiband raster.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Composite Bands* (Data Management).
//!
//! ## Why the catalog needed this
//!
//! The registry is full of multiband **consumers** —
//! `band_collection_statistics`, `principal_component_analysis`,
//! `spectral_index`, `spectral_angle_mapper`, `linear_spectral_unmixing`, the
//! OBIA suite, the cube tools — and had no general way to *assemble* a
//! multiband raster. The only builder was the bundled `create_colour_composite`,
//! hard-wired to exactly three bands for RGB display (`split_colour_composite`
//! is its inverse).
//!
//! That mattered because satellite imagery arrives one file per band: a Landsat
//! or Sentinel delivery is a directory of single-band GeoTIFFs. Without this
//! tool none of the consumers above could be pointed at real data.
//!
//! ## No-data policy
//!
//! `any` (the default) writes no-data wherever *any* contributing band is
//! no-data, which is what the spectral tools want — a pixel missing one band is
//! not a usable spectrum. `all` is the permissive alternative, keeping a cell
//! whenever at least one band has a value.
//!
//! ## Output type
//!
//! Output defaults to the widest input data type rather than `F32`. A 32-bit
//! float stops representing integers exactly at 2^24, so an integer class-code
//! band composited into an `F32` stack would start rounding — the same trap
//! `observer_points` hit in round 18.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::args_common::{band_index, choice_or};
use crate::common::{load_input_raster, parse_optional_output, write_or_store_output};
use crate::raster_stack::raster_like_multiband;

const NODATA_POLICIES: [&str; 2] = ["any", "all"];

pub struct CompositeBandsTool;

impl Tool for CompositeBandsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "composite_bands",
            display_name: "Composite Bands",
            summary: "Stacks several single-band rasters into one multiband raster, in the order given (ArcGIS Composite Bands). The bundled create_colour_composite is fixed at three RGB bands, so there was no way to assemble the N-band input that band_collection_statistics, principal_component_analysis, spectral_index and the OBIA suite all require.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "inputs",
                    description: "Comma- or semicolon-separated raster paths, in output band order. At least one.",
                    required: true,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band taken from each input when an input is itself multiband (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "nodata_policy",
                    description: "'any' (default): a cell is no-data if any band is no-data. 'all': only if every band is.",
                    required: false,
                },
                ToolParamSpec {
                    name: "nodata",
                    description: "No-data value for the output (default: the first input's no-data).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output multiband raster. If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        let inputs = parse_inputs(args)?;
        if inputs.is_empty() {
            return Err(ToolError::Validation(
                "'inputs' must list at least one raster".to_string(),
            ));
        }
        band_index(args, "band")?;
        choice_or(args, "nodata_policy", &NODATA_POLICIES, "any")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let inputs = parse_inputs(args)?;
        if inputs.is_empty() {
            return Err(ToolError::Validation(
                "'inputs' must list at least one raster".to_string(),
            ));
        }
        let band = band_index(args, "band")?;
        let policy = choice_or(args, "nodata_policy", &NODATA_POLICIES, "any")?;
        let output = parse_optional_output(args, "output")?;

        let first = load_input_raster(&inputs[0])?;
        let rows = first.rows;
        let cols = first.cols;
        let out_nodata = args
            .get("nodata")
            .and_then(Value::as_f64)
            .unwrap_or(first.nodata);
        ctx.progress.info(&format!(
            "compositing {} raster(s) into {rows}x{cols}",
            inputs.len()
        ));

        let mut bands: Vec<Vec<f64>> = Vec::with_capacity(inputs.len());
        // Per-cell validity across bands, tracked as we go so the no-data
        // policy can be applied in one pass at the end.
        let mut valid_count = vec![0_u32; rows * cols];
        let mut widest = first.data_type;

        for (i, path) in inputs.iter().enumerate() {
            let raster = if i == 0 {
                first.clone()
            } else {
                load_input_raster(path)?
            };
            // Compare against the reference explicitly rather than via
            // check_alignment_refs, so the message can name the offending file.
            if raster.rows != rows || raster.cols != cols {
                return Err(ToolError::Validation(format!(
                    "input {i} ('{path}') is {}x{}, but input 0 is {rows}x{cols}; composite_bands \
                     does not resample — align the inputs first (warp_raster / reproject_raster)",
                    raster.rows, raster.cols
                )));
            }
            let aligned = |a: f64, b: f64| (a - b).abs() <= 1e-6 * a.abs().max(b.abs()).max(1.0);
            if !aligned(raster.x_min, first.x_min)
                || !aligned(raster.y_min, first.y_min)
                || !aligned(raster.cell_size_x, first.cell_size_x)
                || !aligned(raster.cell_size_y, first.cell_size_y)
            {
                return Err(ToolError::Validation(format!(
                    "input {i} ('{path}') has a different origin or cell size from input 0; \
                     composite_bands does not resample"
                )));
            }
            if raster.crs.epsg.is_some()
                && first.crs.epsg.is_some()
                && raster.crs.epsg != first.crs.epsg
            {
                return Err(ToolError::Validation(format!(
                    "input {i} ('{path}') is EPSG {:?} but input 0 is EPSG {:?}",
                    raster.crs.epsg, first.crs.epsg
                )));
            }
            if band as usize >= raster.bands {
                return Err(ToolError::Validation(format!(
                    "band {} out of range for input {i} ('{path}'), which has {} band(s)",
                    band + 1,
                    raster.bands
                )));
            }

            widest = wider(widest, raster.data_type);
            let mut buf = vec![out_nodata; rows * cols];
            for r in 0..rows {
                for c in 0..cols {
                    let v = raster.get(band, r as isize, c as isize);
                    if v == raster.nodata || !v.is_finite() {
                        continue;
                    }
                    buf[r * cols + c] = v;
                    valid_count[r * cols + c] += 1;
                }
            }
            bands.push(buf);
        }

        let n = bands.len() as u32;
        let mut masked = 0_u64;
        for cell in 0..rows * cols {
            let keep = match policy {
                "all" => valid_count[cell] > 0,
                _ => valid_count[cell] == n,
            };
            if !keep {
                masked += 1;
                for b in bands.iter_mut() {
                    b[cell] = out_nodata;
                }
            }
        }

        let out = raster_like_multiband(&first, &bands, out_nodata, widest)?;
        let out_path = write_or_store_output(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("band_count".to_string(), json!(n));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        outputs.insert("masked_cells".to_string(), json!(masked));
        outputs.insert("data_type".to_string(), json!(format!("{widest:?}")));
        Ok(ToolRunResult { outputs })
    }
}

/// Accepts `inputs` as a delimited string (the repo convention, see
/// `build_seamlines`) or as a JSON array, since both read naturally from a host
/// UI and rejecting one of them would be a gratuitous papercut.
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

/// Ranks the data types by representable range so a stack never narrows below
/// its widest contributor.
fn wider(a: DataType, b: DataType) -> DataType {
    let rank = |d: DataType| match d {
        DataType::U8 | DataType::I8 => 0,
        DataType::U16 | DataType::I16 => 1,
        DataType::U32 | DataType::I32 => 2,
        DataType::F32 => 3,
        _ => 4,
    };
    if rank(b) > rank(a) {
        b
    } else {
        a
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

    #[allow(clippy::too_many_arguments)]
    fn raster_at(
        rows: usize,
        cols: usize,
        data: &[f64],
        x_min: f64,
        cell: f64,
        dt: DataType,
        epsg: Option<u32>,
    ) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min,
            y_min: 0.0,
            cell_size: cell,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: dt,
            crs: CrsInfo {
                epsg,
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

    fn raster(rows: usize, cols: usize, data: &[f64]) -> String {
        raster_at(rows, cols, data, 0.0, 1.0, DataType::F64, Some(3857))
    }

    fn run(args: Value) -> (Raster, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = CompositeBandsTool.run(&args, &ctx()).unwrap();
        let out = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        (out, res)
    }

    #[test]
    fn stacks_inputs_in_the_order_given() {
        let a = raster(1, 2, &[1.0, 2.0]);
        let b = raster(1, 2, &[10.0, 20.0]);
        let c = raster(1, 2, &[100.0, 200.0]);
        let (out, res) = run(json!({"inputs": format!("{a},{b},{c}")}));
        assert_eq!(out.bands, 3);
        assert_eq!(res.outputs["band_count"], json!(3));
        assert_eq!(out.get(0, 0, 0), 1.0);
        assert_eq!(out.get(1, 0, 0), 10.0);
        assert_eq!(out.get(2, 0, 0), 100.0);
        assert_eq!(out.get(2, 0, 1), 200.0);
    }

    #[test]
    fn accepts_a_json_array_as_well_as_a_delimited_string() {
        let a = raster(1, 2, &[1.0, 2.0]);
        let b = raster(1, 2, &[3.0, 4.0]);
        let (out, _) = run(json!({"inputs": [a, b]}));
        assert_eq!(out.bands, 2);
        assert_eq!(out.get(1, 0, 1), 4.0);
    }

    #[test]
    fn a_single_input_is_a_valid_one_band_composite() {
        let a = raster(1, 2, &[5.0, 6.0]);
        let (out, _) = run(json!({"inputs": a}));
        assert_eq!(out.bands, 1);
        assert_eq!(out.get(0, 0, 0), 5.0);
    }

    #[test]
    fn nodata_policy_any_masks_the_whole_stack_at_that_cell() {
        // The spectral tools' requirement: a pixel missing one band is not a
        // usable spectrum, so every band must be masked there.
        let a = raster(1, 2, &[1.0, -9999.0]);
        let b = raster(1, 2, &[2.0, 7.0]);
        let (out, res) = run(json!({"inputs": format!("{a},{b}")}));
        assert_eq!(out.get(0, 0, 1), out.nodata);
        assert_eq!(out.get(1, 0, 1), out.nodata, "band 1 must be masked too");
        assert_eq!(out.get(0, 0, 0), 1.0, "the valid cell survives");
        assert_eq!(res.outputs["masked_cells"], json!(1));
    }

    #[test]
    fn nodata_policy_all_keeps_a_cell_with_one_valid_band() {
        let a = raster(1, 2, &[1.0, -9999.0]);
        let b = raster(1, 2, &[2.0, 7.0]);
        let (out, res) = run(json!({
            "inputs": format!("{a},{b}"), "nodata_policy": "all",
        }));
        assert_eq!(out.get(1, 0, 1), 7.0);
        assert_eq!(out.get(0, 0, 1), out.nodata, "the missing band stays no-data");
        assert_eq!(res.outputs["masked_cells"], json!(0));
    }

    #[test]
    fn a_grid_mismatch_names_the_offending_input() {
        let a = raster(1, 2, &[1.0, 2.0]);
        let b = raster(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        let args: ToolArgs =
            serde_json::from_value(json!({"inputs": format!("{a},{b}")})).unwrap();
        let err = CompositeBandsTool.run(&args, &ctx()).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("input 1"), "got: {msg}");
        assert!(msg.contains("does not resample"), "got: {msg}");
    }

    #[test]
    fn a_shifted_origin_is_rejected_even_at_the_same_size() {
        // Same rows/cols, different origin: silently stacking these would
        // misregister every band by a cell.
        let a = raster_at(1, 2, &[1.0, 2.0], 0.0, 1.0, DataType::F64, Some(3857));
        let b = raster_at(1, 2, &[3.0, 4.0], 5.0, 1.0, DataType::F64, Some(3857));
        let args: ToolArgs =
            serde_json::from_value(json!({"inputs": format!("{a},{b}")})).unwrap();
        let err = CompositeBandsTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("origin or cell size"), "{err}");
    }

    #[test]
    fn a_crs_mismatch_is_rejected() {
        let a = raster_at(1, 2, &[1.0, 2.0], 0.0, 1.0, DataType::F64, Some(3857));
        let b = raster_at(1, 2, &[3.0, 4.0], 0.0, 1.0, DataType::F64, Some(4326));
        let args: ToolArgs =
            serde_json::from_value(json!({"inputs": format!("{a},{b}")})).unwrap();
        assert!(CompositeBandsTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn the_output_type_does_not_narrow_below_the_widest_input() {
        // An I32 band carrying class codes above 2^24 would start rounding if
        // the stack were written as F32.
        let a = raster_at(1, 1, &[20_000_001.0], 0.0, 1.0, DataType::I32, Some(3857));
        let b = raster_at(1, 1, &[1.0], 0.0, 1.0, DataType::U8, Some(3857));
        let (out, _) = run(json!({"inputs": format!("{a},{b}")}));
        assert_eq!(out.get(0, 0, 0), 20_000_001.0);
    }

    #[test]
    fn band_selects_from_a_multiband_input() {
        let mut r = Raster::new(RasterConfig {
            cols: 1,
            rows: 1,
            bands: 2,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F64,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        r.set(0, 0, 0, 11.0).unwrap();
        r.set(1, 0, 0, 22.0).unwrap();
        let id = wbraster::memory_store::put_raster(r);
        let path = wbraster::memory_store::make_raster_memory_path(&id);
        let (out, _) = run(json!({"inputs": path, "band": 2}));
        assert_eq!(out.get(0, 0, 0), 22.0);
    }

    #[test]
    fn an_out_of_range_band_names_the_input() {
        let a = raster(1, 1, &[1.0]);
        let args: ToolArgs = serde_json::from_value(json!({"inputs": a, "band": 5})).unwrap();
        let err = CompositeBandsTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("band 5"), "{err}");
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CompositeBandsTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"inputs": ""})));
        assert!(bad(json!({"inputs": "a.tif", "band": 0})));
        assert!(bad(json!({"inputs": "a.tif", "nodata_policy": "some"})));
        assert!(bad(json!({"inputs": 42})));
    }
}
