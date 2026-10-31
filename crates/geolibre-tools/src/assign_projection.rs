//! GeoLibre tools: assign a coordinate reference system to a raster, vector, or
//! LiDAR dataset **without** transforming its coordinates.
//!
//! The whitebox catalog advertises `assign_projection_{raster,vector,lidar}` but
//! the WASM binary never implemented them, so in GeoLibre's local (WASM) mode
//! they rendered a dead-end form with no parameters (GeoLibre#1355). These
//! implementations fill that gap. Unlike `reproject_*` (which warps the grid /
//! reprojects geometries), assign-projection only rewrites the stored CRS
//! metadata to the given EPSG code, for datasets whose coordinates are already
//! in that CRS but that carry a missing or wrong projection tag.

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata,
    ToolParamSpec, ToolRunResult,
};
use wbraster::CrsInfo;

use crate::common::{load_input_raster, parse_optional_output, write_or_store_output};
use crate::lidar_common::{load_input_cloud, write_or_store_cloud};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// The three `input`/`epsg`/`output` parameter specs shared by every
/// assign-projection tool. `kind` names the dataset in the descriptions.
fn assign_params(kind: &'static str) -> Vec<ToolParamSpec> {
    let input_desc: &'static str = match kind {
        "raster" => "Input raster file path (or in-memory handle).",
        "vector" => "Input vector file path, format auto-detected (or in-memory handle).",
        _ => "Input LiDAR file path (LAS/LAZ/COPC, or in-memory handle).",
    };
    let output_desc: &'static str = match kind {
        "raster" => {
            "Optional output raster path. If omitted, the result is stored in memory."
        }
        "vector" => {
            "Optional output vector path (driver from its extension). If omitted, stored in memory."
        }
        _ => "Optional output LiDAR path. If omitted, the result is stored in memory.",
    };
    vec![
        ToolParamSpec {
            name: "input",
            description: input_desc,
            required: true,
        },
        ToolParamSpec {
            name: "epsg",
            description: "EPSG code of the CRS to assign, e.g. 4326 (WGS84) or 3857 (Web Mercator). Coordinates are not transformed.",
            required: true,
        },
        ToolParamSpec {
            name: "output",
            description: output_desc,
            required: false,
        },
    ]
}

/// Validates that `input` is a non-empty string and `epsg` is a positive code.
fn validate_assign(args: &ToolArgs) -> Result<(), ToolError> {
    if args
        .get("input")
        .and_then(Value::as_str)
        .map(str::trim)
        .unwrap_or("")
        .is_empty()
    {
        return Err(ToolError::Validation(
            "missing required string parameter 'input'".to_string(),
        ));
    }
    parse_epsg(args)?;
    Ok(())
}

/// Reads the required `input` path parameter.
fn input_path(args: &ToolArgs) -> Result<&str, ToolError> {
    args.get("input")
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation("missing required parameter 'input'".to_string()))
}

/// Assigns a CRS to a raster's metadata without warping the grid.
pub struct AssignProjectionRasterTool;

impl Tool for AssignProjectionRasterTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "assign_projection_raster",
            display_name: "Assign Projection Raster",
            summary: "Assign a coordinate reference system to a raster without warping its cells.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: assign_params("raster"),
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        validate_assign(args)
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = input_path(args)?;
        let epsg = parse_epsg(args)?;
        let output = parse_optional_output(args, "output")?;

        let mut raster = load_input_raster(input)?;
        ctx.progress.info(&format!("assigning EPSG:{epsg}"));
        raster.crs = CrsInfo::from_epsg(epsg);
        let (rows, cols) = (raster.rows, raster.cols);
        let out_path = write_or_store_output(raster, output)?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("epsg".to_string(), json!(epsg));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        Ok(ToolRunResult { outputs })
    }
}

/// Assigns a CRS to a vector layer's metadata without reprojecting geometries.
pub struct AssignProjectionVectorTool;

impl Tool for AssignProjectionVectorTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "assign_projection_vector",
            display_name: "Assign Projection Vector",
            summary: "Assign a coordinate reference system to a vector layer without reprojecting its geometries.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: assign_params("vector"),
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        validate_assign(args)
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = input_path(args)?;
        let epsg = parse_epsg(args)?;
        let output = parse_optional_str(args, "output")?;

        let mut layer = load_input_layer(input)?;
        ctx.progress.info(&format!("assigning EPSG:{epsg}"));
        layer.assign_crs_epsg(epsg);
        let feature_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("epsg".to_string(), json!(epsg));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        Ok(ToolRunResult { outputs })
    }
}

/// Assigns a CRS to a LiDAR point cloud's metadata without reprojecting points.
pub struct AssignProjectionLidarTool;

impl Tool for AssignProjectionLidarTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "assign_projection_lidar",
            display_name: "Assign Projection Lidar",
            summary: "Assign a coordinate reference system to a LiDAR point cloud without reprojecting its points.",
            category: ToolCategory::Lidar,
            license_tier: LicenseTier::Open,
            params: assign_params("lidar"),
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        validate_assign(args)
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = input_path(args)?;
        let epsg = parse_epsg(args)?;
        let output = parse_optional_str(args, "output")?;

        let mut cloud = load_input_cloud(input)?;
        ctx.progress.info(&format!("assigning EPSG:{epsg}"));
        cloud.assign_crs_epsg(epsg);
        let point_count = cloud.point_count();
        let out_path = write_or_store_cloud(cloud, output)?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("epsg".to_string(), json!(epsg));
        outputs.insert("point_count".to_string(), json!(point_count));
        Ok(ToolRunResult { outputs })
    }
}

/// Parses the required positive `epsg` parameter (accepts a number or a numeric
/// string). Shared by every assign-projection tool.
fn parse_epsg(args: &ToolArgs) -> Result<u32, ToolError> {
    let raw = args
        .get("epsg")
        .ok_or_else(|| ToolError::Validation("missing required parameter 'epsg'".to_string()))?;
    let code = match raw {
        Value::Number(n) => n.as_u64(),
        Value::String(s) => s.trim().parse::<u64>().ok(),
        _ => None,
    };
    match code {
        Some(c) if c > 0 && c <= u32::MAX as u64 => Ok(c as u32),
        _ => Err(ToolError::Validation(
            "parameter 'epsg' must be a positive EPSG code".to_string(),
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn args(input: &str, epsg: Value) -> ToolArgs {
        serde_json::from_value(json!({ "input": input, "epsg": epsg })).unwrap()
    }

    #[test]
    fn assigns_crs_to_a_raster_without_warping() {
        use wbraster::{DataType, Raster, RasterConfig};
        let cfg = RasterConfig {
            cols: 2,
            rows: 1,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: -9999.0,
            data_type: DataType::F64,
            // Start mislabeled as 4326 to prove the tool overwrites it.
            crs: CrsInfo::from_epsg(4326),
            metadata: Vec::new(),
        };
        let raster = Raster::from_data(cfg, vec![1.0, 2.0]).unwrap();
        let id = wbraster::memory_store::put_raster(raster);
        let input = wbraster::memory_store::make_raster_memory_path(&id);

        let out = AssignProjectionRasterTool
            .run(&args(&input, json!(32610)), &ctx())
            .unwrap();
        assert_eq!(out.outputs["epsg"], json!(32610));

        let result = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        assert_eq!(result.crs.epsg, Some(32610));
        // The grid is untouched: same cells, same values.
        assert_eq!((result.rows, result.cols), (1, 2));
        assert_eq!(result.get(0, 0, 0), 1.0);
        assert_eq!(result.get(0, 0, 1), 2.0);
    }

    #[test]
    fn assigns_crs_to_a_vector_layer() {
        use wbvector::{memory_store, Layer};
        let mut layer = Layer::new("test");
        layer.set_crs_epsg(Some(4326));
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let out = AssignProjectionVectorTool
            .run(&args(&input, json!(3857)), &ctx())
            .unwrap();
        assert_eq!(out.outputs["epsg"], json!(3857));

        let result =
            crate::vector_common::load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        assert_eq!(result.crs_epsg(), Some(3857));
    }

    #[test]
    fn assigns_crs_to_a_point_cloud() {
        use wblidar::{memory_store, point::PointRecord, PointCloud};
        let cloud = PointCloud {
            points: vec![PointRecord::default(), PointRecord::default()],
            crs: None,
        };
        let id = memory_store::put_lidar(cloud);
        let input = memory_store::make_lidar_memory_path(&id);

        let out = AssignProjectionLidarTool
            .run(&args(&input, json!(26910)), &ctx())
            .unwrap();
        assert_eq!(out.outputs["epsg"], json!(26910));
        assert_eq!(out.outputs["point_count"], json!(2));

        let result =
            crate::lidar_common::load_input_cloud(out.outputs["output"].as_str().unwrap()).unwrap();
        assert_eq!(result.crs.and_then(|c| c.epsg), Some(26910));
    }

    #[test]
    fn rejects_missing_input_and_bad_epsg() {
        let missing: ToolArgs = serde_json::from_value(json!({ "epsg": 4326 })).unwrap();
        assert!(AssignProjectionVectorTool.validate(&missing).is_err());

        let bad_epsg = args("memory://vector/whatever", json!(0));
        assert!(AssignProjectionRasterTool.validate(&bad_epsg).is_err());

        let str_epsg = args("memory://vector/whatever", json!("32610"));
        // A numeric string is accepted by parse_epsg, so validation passes.
        assert!(AssignProjectionLidarTool.validate(&str_epsg).is_ok());
    }
}
