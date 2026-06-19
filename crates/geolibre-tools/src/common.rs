//! Shared helpers for GeoLibre tools: raster input/output over file paths or
//! in-memory (`memory://`) handles, and small numeric utilities.
//!
//! These mirror the private helpers in `raster_normalize.rs` but are shared so
//! the LiDAR/DEM tools (filters, sink extraction, depression delineation) can
//! reuse one implementation.

use std::path::Path;
use std::sync::Arc;

use serde_json::Value;
use wbcore::{ToolArgs, ToolError};
use wbraster::{memory_store, DataType, Raster, RasterFormat};

/// Parses an optional string `output` parameter (absent / null / empty -> None).
pub fn parse_optional_output<'a>(
    args: &'a ToolArgs,
    key: &str,
) -> Result<Option<&'a str>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => Ok(Some(s.as_str())),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a string when provided"
        ))),
    }
}

/// Loads a raster from a file path or an in-memory (`memory://`) handle.
pub fn load_input_raster(path: &str) -> Result<Raster, ToolError> {
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

/// Writes `raster` to `output_path`, or stores it in memory and returns a
/// `memory://` handle when no path is given.
pub fn write_or_store_output(
    raster: Raster,
    output_path: Option<&str>,
) -> Result<String, ToolError> {
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

/// Writes a plain-text artifact (e.g. a CSV attribute table) to a file path,
/// creating parent directories as needed. In-memory output is not supported for
/// text artifacts, so a path is required.
pub fn write_text_output(text: &str, output_path: &str) -> Result<(), ToolError> {
    if let Some(parent) = Path::new(output_path).parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|e| {
                ToolError::Execution(format!("failed creating output directory: {e}"))
            })?;
        }
    }
    std::fs::write(output_path, text)
        .map_err(|e| ToolError::Execution(format!("failed writing text output: {e}")))
}

/// Reads one band of a raster into a dense, row-major `f64` buffer. No-data
/// cells are returned verbatim (the raster's `nodata` value), so callers should
/// compare against `raster.nodata`.
pub fn band_to_vec(raster: &Raster, band: isize) -> Vec<f64> {
    let rows = raster.rows;
    let cols = raster.cols;
    let mut out = vec![0.0_f64; rows * cols];
    for row in 0..rows {
        for col in 0..cols {
            out[row * cols + col] = raster.get(band, row as isize, col as isize);
        }
    }
    out
}

/// Builds a new single-band raster that copies `template`'s geometry and CRS but
/// uses the supplied data buffer, no-data value, and data type.
pub fn raster_like_with_data(
    template: &Raster,
    data: Vec<f64>,
    nodata: f64,
    data_type: DataType,
) -> Result<Raster, ToolError> {
    let mut out = Raster::new_like(template);
    out.nodata = nodata;
    out.data_type = data_type;
    let rows = out.rows;
    let cols = out.cols;
    if data.len() != rows * cols {
        return Err(ToolError::Execution(format!(
            "output buffer length {} does not match {rows}x{cols}",
            data.len()
        )));
    }
    for row in 0..rows {
        for col in 0..cols {
            out.set(0, row as isize, col as isize, data[row * cols + col])
                .map_err(|e| ToolError::Execution(format!("failed writing cell: {e}")))?;
        }
    }
    Ok(out)
}
