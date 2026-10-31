//! Example GeoLibre tool: min-max normalize a raster band to the range [0, 1].
//!
//! Use this as a template for new tools. The pattern is:
//!   1. `metadata()` declares the tool id, summary, and parameters.
//!   2. `validate()` cheaply checks the args before any heavy work.
//!   3. `run()` reads inputs (file path or `memory://` handle), does the work,
//!      and writes the output (to a file path, or stores it in memory).
//!
//! Input/output rasters are passed as file paths, which the WASI runner backs
//! with its in-memory `/work` filesystem. The default `manifest()` (derived from
//! `metadata()`) is sufficient, so it is not implemented here.

use std::path::Path;
use std::sync::Arc;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata,
    ToolParamSpec, ToolRunResult,
};
use wbraster::{memory_store, Raster, RasterFormat};

/// Linearly rescales one band of a raster so its minimum maps to 0 and its
/// maximum maps to 1. No-data cells are preserved.
pub struct RasterNormalizeTool;

impl Tool for RasterNormalizeTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "raster_normalize",
            display_name: "Raster Normalize",
            summary: "Min-max normalize a raster band to the range [0, 1].",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input raster file path.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output raster path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to normalize (default 1).",
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
        let output = parse_optional_output(args)?;
        let band_1based = args.get("band").and_then(Value::as_u64).unwrap_or(1).max(1);
        let band = (band_1based - 1) as isize;

        let mut raster = load_input_raster(input)?;
        if band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "band {band_1based} out of range (raster has {} band(s))",
                raster.bands
            )));
        }

        let nodata = raster.nodata;
        let rows = raster.rows as isize;
        let cols = raster.cols as isize;

        ctx.progress.info("computing band range");
        let (mut min, mut max) = (f64::INFINITY, f64::NEG_INFINITY);
        for row in 0..rows {
            for col in 0..cols {
                let v = raster.get(band, row, col);
                if v != nodata && v.is_finite() {
                    min = min.min(v);
                    max = max.max(v);
                }
            }
        }

        if !min.is_finite() || !max.is_finite() {
            return Err(ToolError::Execution(
                "raster band contains no valid (non-nodata) values".to_string(),
            ));
        }
        let range = max - min;

        ctx.progress.info("normalizing");
        for row in 0..rows {
            for col in 0..cols {
                let v = raster.get(band, row, col);
                if v != nodata && v.is_finite() {
                    let scaled = if range == 0.0 { 0.0 } else { (v - min) / range };
                    raster
                        .set(band, row, col, scaled)
                        .map_err(|e| ToolError::Execution(format!("failed writing cell: {e}")))?;
                }
            }
            ctx.progress.progress((row as f64 + 1.0) / rows as f64);
        }

        let out_path = write_or_store_output(raster, output)?;
        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("min".to_string(), json!(min));
        outputs.insert("max".to_string(), json!(max));
        Ok(ToolRunResult { outputs })
    }
}

/// Parses an optional `output` string parameter (absent / null / empty -> None).
fn parse_optional_output(args: &ToolArgs) -> Result<Option<&str>, ToolError> {
    match args.get("output") {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => Ok(Some(s.as_str())),
        Some(_) => Err(ToolError::Validation(
            "parameter 'output' must be a string when provided".to_string(),
        )),
    }
}

/// Loads a raster from a file path or an in-memory (`memory://`) handle.
fn load_input_raster(path: &str) -> Result<Raster, ToolError> {
    if memory_store::raster_is_memory_path(path) {
        let id = memory_store::raster_path_to_id(path)
            .ok_or_else(|| ToolError::Validation("malformed in-memory raster path".to_string()))?;
        let arc: Arc<Raster> = memory_store::get_raster_arc_by_id(id)
            .ok_or_else(|| ToolError::Validation(format!("unknown in-memory raster id '{id}'")))?;
        return Ok((*arc).clone());
    }

    Raster::read(path)
        .map_err(|e| ToolError::Execution(format!("failed reading input raster: {e}")))
}

/// Writes the raster to `output_path`, or stores it in memory and returns a
/// `memory://` handle when no path is given.
fn write_or_store_output(raster: Raster, output_path: Option<&str>) -> Result<String, ToolError> {
    match output_path {
        Some(output_path) => {
            if let Some(parent) = Path::new(output_path).parent() {
                if !parent.as_os_str().is_empty() {
                    std::fs::create_dir_all(parent).map_err(|e| {
                        ToolError::Execution(format!("failed creating output directory: {e}"))
                    })?;
                }
            }
            let fmt = RasterFormat::for_output_path(output_path)
                .map_err(|e| ToolError::Validation(format!("unsupported output path: {e}")))?;
            raster
                .write(output_path, fmt)
                .map_err(|e| ToolError::Execution(format!("failed writing output raster: {e}")))?;
            Ok(output_path.to_string())
        }
        None => {
            let id = memory_store::put_raster(raster);
            Ok(memory_store::make_raster_memory_path(&id))
        }
    }
}
