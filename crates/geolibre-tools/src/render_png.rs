//! GeoLibre tool: render a single raster band to a PNG image through a colormap.
//!
//! The whitebox suite has RGB composite tools but nothing that renders one
//! continuous band through a colormap to a web-displayable PNG. No-data cells
//! become fully transparent, so the image overlays cleanly.

use std::path::Path;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata,
    ToolParamSpec, ToolRunResult,
};

use crate::common::load_input_raster;
use crate::render::{encode_png_rgba, normalize, Colormap};

/// Renders one band of a raster to an RGBA PNG using a named colormap, stretched
/// across the band's value range (or an explicit min/max).
pub struct RenderPngTool;

impl Tool for RenderPngTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "render_raster_png",
            display_name: "Render Raster to PNG",
            summary: "Render a single raster band to a PNG image through a colormap (no-data becomes transparent).",
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
                    description: "Output PNG file path (e.g. /work/preview.png).",
                    required: true,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to render (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "colormap",
                    description: "Colormap: viridis (default), magma, turbo, terrain, or grayscale.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min",
                    description: "Value mapped to the low end of the colormap (default: band minimum).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max",
                    description: "Value mapped to the high end of the colormap (default: band maximum).",
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
        if args.get("output").and_then(Value::as_str).is_none() {
            return Err(ToolError::Validation(
                "missing required string parameter 'output'".to_string(),
            ));
        }
        if let Some(c) = args.get("colormap").and_then(Value::as_str) {
            Colormap::parse(c)?;
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args
            .get("input")
            .and_then(Value::as_str)
            .ok_or_else(|| ToolError::Validation("missing required parameter 'input'".to_string()))?;
        let output = args
            .get("output")
            .and_then(Value::as_str)
            .ok_or_else(|| ToolError::Validation("missing required parameter 'output'".to_string()))?;
        let band_1based = args.get("band").and_then(Value::as_u64).unwrap_or(1).max(1);
        let band = (band_1based - 1) as isize;
        let colormap = match args.get("colormap").and_then(Value::as_str) {
            Some(c) => Colormap::parse(c)?,
            None => Colormap::Viridis,
        };

        let raster = load_input_raster(input)?;
        if band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "band {band_1based} out of range (raster has {} band(s))",
                raster.bands
            )));
        }

        // Resolve the stretch range from args or the band's statistics.
        let stats = raster
            .statistics_band(band)
            .map_err(|e| ToolError::Execution(format!("failed computing band statistics: {e}")))?;
        let min = args.get("min").and_then(Value::as_f64).unwrap_or(stats.min);
        let max = args.get("max").and_then(Value::as_f64).unwrap_or(stats.max);
        if !min.is_finite() || !max.is_finite() {
            return Err(ToolError::Execution(
                "raster band has no finite values to render".to_string(),
            ));
        }

        let rows = raster.rows as isize;
        let cols = raster.cols as isize;
        let nodata = raster.nodata;

        ctx.progress.info("rendering");
        let mut rgba = vec![0u8; (rows * cols * 4) as usize];
        for row in 0..rows {
            for col in 0..cols {
                let v = raster.get(band, row, col);
                let idx = ((row * cols + col) * 4) as usize;
                if v == nodata || !v.is_finite() {
                    continue; // leave fully transparent
                }
                let [r, g, b] = colormap.rgb(normalize(v, min, max));
                rgba[idx] = r;
                rgba[idx + 1] = g;
                rgba[idx + 2] = b;
                rgba[idx + 3] = 255;
            }
            ctx.progress.progress((row as f64 + 1.0) / rows as f64);
        }

        let png = encode_png_rgba(&rgba, cols as u32, rows as u32)?;
        write_bytes(output, &png)?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(output));
        outputs.insert("width".to_string(), json!(cols));
        outputs.insert("height".to_string(), json!(rows));
        outputs.insert("min".to_string(), json!(min));
        outputs.insert("max".to_string(), json!(max));
        Ok(ToolRunResult { outputs })
    }
}

/// Writes bytes to a path, creating parent directories as needed.
fn write_bytes(path: &str, bytes: &[u8]) -> Result<(), ToolError> {
    if let Some(parent) = Path::new(path).parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|e| {
                ToolError::Execution(format!("failed creating output directory: {e}"))
            })?;
        }
    }
    std::fs::write(path, bytes)
        .map_err(|e| ToolError::Execution(format!("failed writing output file: {e}")))
}
