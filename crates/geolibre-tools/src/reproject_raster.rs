//! GeoLibre tool: reproject (warp) a raster to a target EPSG CRS.
//!
//! The whitebox suite ships `reproject_vector` but no raster warp; this fills
//! that gap by wrapping `wbraster`'s built-in reprojection. The source CRS is
//! read from the raster's own metadata, so only a destination EPSG is required.

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata,
    ToolParamSpec, ToolRunResult,
};
use wbraster::ResampleMethod;

use crate::common::{load_input_raster, parse_optional_output, write_or_store_output};

/// Warps a raster into a target EPSG coordinate reference system, resampling
/// the grid with the requested method.
pub struct ReprojectRasterTool;

impl Tool for ReprojectRasterTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "reproject_raster",
            display_name: "Reproject Raster",
            summary: "Reproject (warp) a raster into a target EPSG coordinate reference system.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input raster file path. Must carry a source CRS (EPSG or WKT).",
                    required: true,
                },
                ToolParamSpec {
                    name: "epsg",
                    description: "Destination EPSG code, e.g. 3857 (Web Mercator) or 4326 (WGS84).",
                    required: true,
                },
                ToolParamSpec {
                    name: "method",
                    description: "Resampling method: nearest (default), bilinear, cubic, lanczos, average, min, max, mode, median.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output raster path. If omitted, the result is stored in memory.",
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
        if parse_epsg(args).is_err() {
            return Err(ToolError::Validation(
                "missing or invalid required parameter 'epsg' (a positive EPSG code)".to_string(),
            ));
        }
        if let Some(m) = args.get("method").and_then(Value::as_str) {
            parse_resample(m)?;
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args
            .get("input")
            .and_then(Value::as_str)
            .ok_or_else(|| ToolError::Validation("missing required parameter 'input'".to_string()))?;
        let dst_epsg = parse_epsg(args)?;
        let method = match args.get("method").and_then(Value::as_str) {
            Some(m) => parse_resample(m)?,
            None => ResampleMethod::Nearest,
        };
        let output = parse_optional_output(args, "output")?;

        let raster = load_input_raster(input)?;
        let src_epsg = raster.crs.epsg;
        if src_epsg.is_none() && raster.crs.wkt.is_none() && raster.crs.proj4.is_none() {
            return Err(ToolError::Validation(
                "input raster has no source CRS (EPSG/WKT/PROJ); cannot reproject".to_string(),
            ));
        }

        ctx.progress
            .info(&format!("reprojecting to EPSG:{dst_epsg}"));
        let reprojected = raster
            .reproject_to_epsg(dst_epsg, method)
            .map_err(|e| ToolError::Execution(format!("reprojection failed: {e}")))?;

        let (rows, cols) = (reprojected.rows, reprojected.cols);
        let out_path = write_or_store_output(reprojected, output)?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("dst_epsg".to_string(), json!(dst_epsg));
        if let Some(src) = src_epsg {
            outputs.insert("src_epsg".to_string(), json!(src));
        }
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        Ok(ToolRunResult { outputs })
    }
}

/// Parses the required positive `epsg` parameter (accepts a number or a numeric
/// string).
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

/// Maps a resampling-method name to a `ResampleMethod`.
pub fn parse_resample(name: &str) -> Result<ResampleMethod, ToolError> {
    match name.trim().to_ascii_lowercase().as_str() {
        "nearest" | "nn" => Ok(ResampleMethod::Nearest),
        "bilinear" | "linear" => Ok(ResampleMethod::Bilinear),
        "cubic" | "bicubic" => Ok(ResampleMethod::Cubic),
        "lanczos" => Ok(ResampleMethod::Lanczos),
        "average" | "mean" => Ok(ResampleMethod::Average),
        "min" | "minimum" => Ok(ResampleMethod::Min),
        "max" | "maximum" => Ok(ResampleMethod::Max),
        "mode" => Ok(ResampleMethod::Mode),
        "median" => Ok(ResampleMethod::Median),
        other => Err(ToolError::Validation(format!(
            "unknown resampling method '{other}'"
        ))),
    }
}
