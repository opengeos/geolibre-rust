//! GeoLibre tool: convert a vector dataset between formats.
//!
//! Reads any format `wbvector` understands (GeoJSON, Shapefile, FlatGeobuf,
//! GeoPackage, GeoParquet, GML, GPX, KML, ...) and writes it back out in the
//! format implied by the output path's extension. One tool instead of a matrix
//! of format-specific converters.

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata,
    ToolParamSpec, ToolRunResult,
};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Reads a vector dataset and writes it to another vector format.
pub struct VectorConvertTool;

impl Tool for VectorConvertTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "vector_convert",
            display_name: "Vector Convert",
            summary: "Convert a vector dataset to another format (the output extension picks the driver).",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector file path (format auto-detected) or in-memory handle.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path; the driver is taken from its extension (.geojson, .fgb, .shp, .gpkg, .parquet, ...). If omitted, the layer is stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args.get("input").and_then(Value::as_str).map(str::trim).unwrap_or("").is_empty() {
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
        let output = parse_optional_str(args, "output")?;

        ctx.progress.info("reading input vector");
        let layer = load_input_layer(input)?;
        let feature_count = layer.len();

        ctx.progress.info("writing output vector");
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        Ok(ToolRunResult { outputs })
    }
}
