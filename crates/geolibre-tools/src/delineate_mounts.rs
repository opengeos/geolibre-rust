//! Delineation of the nested hierarchy of elevated features (mounts), ported
//! from `lidar`'s `mounts.py`.
//!
//! A mount is an inverted depression: flip the DEM, extract sinks on the flipped
//! surface, then delineate nested depressions. The output id/level rasters
//! describe the mount hierarchy of the original DEM.

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{
    band_to_vec, load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
    write_text_output,
};
use crate::delineate_depressions::{delineate_core, depression_csv, depression_props_map};
use crate::extract_sinks::extract_sinks_core;
use crate::polygonize::{polygonize_to_geojson, PolygonizeParams};

/// Flips a DEM so that peaks become pits: `flipped = -dem + max + delta`.
/// No-data cells are preserved (the Python original mishandles them; we keep
/// them as no-data).
pub fn flip_dem(dem: &[f64], nodata: f64, delta: f64) -> Vec<f64> {
    let max_elev = dem
        .iter()
        .cloned()
        .filter(|&v| v != nodata)
        .fold(f64::NEG_INFINITY, f64::max);
    dem.iter()
        .map(|&v| if v == nodata { nodata } else { -v + max_elev + delta })
        .collect()
}

/// Delineates the nested hierarchy of mounts in a DEM (`lidar` `DelineateMounts`).
pub struct DelineateMountsTool;

impl Tool for DelineateMountsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "delineate_mounts",
            display_name: "Delineate Mounts",
            summary: "Delineate the nested hierarchy of elevated features (mounts) in a DEM (lidar mounts.py).",
            category: ToolCategory::Hydrology,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input DEM raster file path.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional mount-id raster path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "level_output",
                    description: "Optional output path for the mount-level raster.",
                    required: false,
                },
                ToolParamSpec {
                    name: "csv_output",
                    description: "Optional output path for the mount attribute CSV.",
                    required: false,
                },
                ToolParamSpec {
                    name: "vector_output",
                    description: "Optional output path for mount polygons as GeoJSON (attributes joined).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_size",
                    description: "Minimum number of pixels for a mount (default 10).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_height",
                    description: "Minimum mount height (accepted for compatibility; unused, matching lidar).",
                    required: false,
                },
                ToolParamSpec {
                    name: "interval",
                    description: "Slicing interval for the level-set method (default 0.3).",
                    required: false,
                },
                ToolParamSpec {
                    name: "delta",
                    description: "Base value added when flipping the DEM (default 100).",
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
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args
            .get("input")
            .and_then(Value::as_str)
            .ok_or_else(|| ToolError::Validation("missing required parameter 'input'".to_string()))?;
        let output = parse_optional_output(args, "output")?;
        let level_output = parse_optional_output(args, "level_output")?;
        let csv_output = parse_optional_output(args, "csv_output")?;
        let vector_output = parse_optional_output(args, "vector_output")?;
        let min_size = args.get("min_size").and_then(Value::as_u64).unwrap_or(10);
        let interval = args.get("interval").and_then(Value::as_f64).unwrap_or(0.3);
        let delta = args.get("delta").and_then(Value::as_f64).unwrap_or(100.0);

        let raster = load_input_raster(input)?;
        let rows = raster.rows;
        let cols = raster.cols;
        let nodata = raster.nodata;
        let resolution = raster.cell_size_x;
        let dem = band_to_vec(&raster, 0);

        ctx.progress.info("flipping DEM");
        let flipped = flip_dem(&dem, nodata, delta);

        ctx.progress.info("extracting sinks on flipped DEM");
        let sinks = extract_sinks_core(&flipped, rows, cols, nodata, min_size, 0.0);
        ctx.progress.progress(0.5);

        ctx.progress.info("delineating nested mounts");
        let result = delineate_core(&sinks.sink, rows, cols, resolution, min_size, interval);
        ctx.progress.progress(0.9);

        let mut outputs = std::collections::BTreeMap::new();

        if let Some(path) = vector_output {
            let pmap = depression_props_map(&result.depressions);
            let geojson = polygonize_to_geojson(&PolygonizeParams {
                labels: &result.obj_image,
                rows,
                cols,
                x_min: raster.x_min,
                y_max: raster.y_max(),
                cell_size_x: raster.cell_size_x,
                cell_size_y: raster.cell_size_y,
                epsg: raster.crs.epsg,
                props_by_id: &pmap,
            });
            write_text_output(&geojson, path)?;
            outputs.insert("vector".to_string(), json!(path));
        }

        let obj_raster = raster_like_with_data(&raster, result.obj_image, 0.0, DataType::I32)?;
        let obj_path = write_or_store_output(obj_raster, output)?;
        outputs.insert("output".to_string(), json!(obj_path));

        if let Some(path) = level_output {
            let r = raster_like_with_data(&raster, result.level_image, 0.0, DataType::I32)?;
            outputs.insert("level".to_string(), json!(write_or_store_output(r, Some(path))?));
        }
        if let Some(path) = csv_output {
            write_text_output(&depression_csv(&result.depressions), path)?;
            outputs.insert("csv".to_string(), json!(path));
        }

        outputs.insert("mounts".to_string(), json!(result.depressions.len()));
        ctx.progress.progress(1.0);
        Ok(ToolRunResult { outputs })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flip_inverts_relief() {
        let nodata = -9999.0;
        let dem = vec![1.0, 2.0, 3.0, nodata];
        let flipped = flip_dem(&dem, nodata, 100.0);
        // max is 3; flipped = -v + 3 + 100.
        assert_eq!(flipped[0], 102.0);
        assert_eq!(flipped[2], 100.0);
        assert_eq!(flipped[3], nodata);
        // Relief is inverted: the lowest cell becomes the highest.
        assert!(flipped[0] > flipped[2]);
    }

    #[test]
    fn detects_a_mount() {
        // A plateau with a central graded peak: flipping makes it a pit with
        // relief (a flat-topped peak would flip to a flat pit and produce no
        // slices, matching the Python behavior).
        let rows = 7;
        let cols = 7;
        let nodata = -9999.0;
        let mut dem = vec![1.0; rows * cols];
        for r in 1..6 {
            for c in 1..6 {
                let cheb = (r as i64 - 3).abs().max((c as i64 - 3).abs()) as f64;
                dem[r * cols + c] = 4.0 - cheb; // peak 4 at center, sloping to 2
            }
        }
        let result = {
            let flipped = flip_dem(&dem, nodata, 100.0);
            let sinks = extract_sinks_core(&flipped, rows, cols, nodata, 1, 0.0);
            delineate_core(&sinks.sink, rows, cols, 1.0, 1, 1.0)
        };
        assert!(!result.depressions.is_empty());
    }
}
