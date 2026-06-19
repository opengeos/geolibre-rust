//! GeoLibre tools: read and write GeoParquet, the cloud-native columnar vector
//! format. The whitebox suite has no GeoParquet I/O tool (the capability lives
//! in `wbvector` behind a feature flag); these expose it.
//!
//! - `write_geoparquet`: read any supported vector format, write `.parquet`.
//! - `read_geoparquet`: read `.parquet`, write any supported vector format
//!   (or store in memory).

use parquet::basic::{Compression, ZstdLevel};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata,
    ToolParamSpec, ToolRunResult,
};
use wbvector::geoparquet::{write_with_options, GeoParquetWriteOptions};
use wbvector::Layer;

use crate::hilbert::hilbert_for_point;
use crate::vector_common::{ensure_parent_dir, load_input_layer, parse_optional_str, write_or_store_layer};

/// Converts any supported vector dataset to GeoParquet.
pub struct WriteGeoParquetTool;

impl Tool for WriteGeoParquetTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "write_geoparquet",
            display_name: "Write GeoParquet",
            summary: "Convert a vector dataset (GeoJSON, Shapefile, FlatGeobuf, GeoPackage, ...) to GeoParquet, Hilbert-sorted and ZSTD-compressed by default.",
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
                    description: "Output GeoParquet file path (e.g. /work/data.parquet).",
                    required: true,
                },
                ToolParamSpec {
                    name: "compression",
                    description: "Column compression: zstd (default), snappy, gzip, or uncompressed.",
                    required: false,
                },
                ToolParamSpec {
                    name: "hilbert_sort",
                    description: "Sort features along a Hilbert curve for spatial locality (default true).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "output")?;
        if let Some(c) = args.get("compression").and_then(Value::as_str) {
            parse_compression(c)?;
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = require_str(args, "output")?;
        let compression = match args.get("compression").and_then(Value::as_str) {
            Some(c) => parse_compression(c)?,
            None => Compression::ZSTD(ZstdLevel::default()),
        };
        let hilbert = args
            .get("hilbert_sort")
            .and_then(Value::as_bool)
            .unwrap_or(true);

        ctx.progress.info("reading input vector");
        let mut layer = load_input_layer(input)?;
        let feature_count = layer.len();

        if hilbert {
            ctx.progress.info("Hilbert-sorting features");
            hilbert_sort(&mut layer);
        }

        ctx.progress.info("writing GeoParquet");
        let options = GeoParquetWriteOptions::new().with_compression(compression);
        ensure_parent_dir(output)?;
        write_with_options(&layer, output, &options)
            .map_err(|e| ToolError::Execution(format!("failed writing GeoParquet: {e}")))?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(output));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("compression".to_string(), json!(compression_name(compression)));
        outputs.insert("hilbert_sorted".to_string(), json!(hilbert));
        Ok(ToolRunResult { outputs })
    }
}

/// Reorders a layer's features along a Hilbert curve computed from each
/// feature's bounding-box center over the dataset extent. Features without a
/// geometry keep their relative order at the end.
fn hilbert_sort(layer: &mut Layer) {
    // Dataset extent over all feature bounding boxes.
    let (mut min_x, mut min_y, mut max_x, mut max_y) =
        (f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for f in &layer.features {
        if let Some(b) = f.geometry.as_ref().and_then(|g| g.bbox()) {
            min_x = min_x.min(b.min_x);
            min_y = min_y.min(b.min_y);
            max_x = max_x.max(b.max_x);
            max_y = max_y.max(b.max_y);
        }
    }
    if !min_x.is_finite() {
        return; // nothing georeferenced to sort
    }

    // Key each feature; geometry-less features sort last (key = u64::MAX).
    let mut indexed: Vec<(u64, usize)> = layer
        .features
        .iter()
        .enumerate()
        .map(|(i, f)| {
            let key = match f.geometry.as_ref().and_then(|g| g.bbox()) {
                Some(b) => hilbert_for_point(
                    0.5 * (b.min_x + b.max_x),
                    0.5 * (b.min_y + b.max_y),
                    min_x,
                    min_y,
                    max_x,
                    max_y,
                ),
                None => u64::MAX,
            };
            (key, i)
        })
        .collect();
    indexed.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    let original = std::mem::take(&mut layer.features);
    let mut slots: Vec<Option<_>> = original.into_iter().map(Some).collect();
    layer.features = indexed
        .into_iter()
        .map(|(_, i)| slots[i].take().expect("each index used once"))
        .collect();
}

/// Parses a compression-codec name into a parquet `Compression`.
fn parse_compression(name: &str) -> Result<Compression, ToolError> {
    match name.trim().to_ascii_lowercase().as_str() {
        "zstd" => Ok(Compression::ZSTD(ZstdLevel::default())),
        "snappy" | "snap" => Ok(Compression::SNAPPY),
        "gzip" | "gz" => Ok(Compression::GZIP(Default::default())),
        "uncompressed" | "none" | "off" => Ok(Compression::UNCOMPRESSED),
        other => Err(ToolError::Validation(format!(
            "unknown compression '{other}' (expected zstd, snappy, gzip, or uncompressed)"
        ))),
    }
}

/// Short stable name for a compression codec, for the result JSON.
fn compression_name(c: Compression) -> &'static str {
    match c {
        Compression::ZSTD(_) => "zstd",
        Compression::SNAPPY => "snappy",
        Compression::GZIP(_) => "gzip",
        Compression::UNCOMPRESSED => "uncompressed",
        _ => "other",
    }
}

/// Reads GeoParquet and writes it to another vector format (or memory).
pub struct ReadGeoParquetTool;

impl Tool for ReadGeoParquetTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "read_geoparquet",
            display_name: "Read GeoParquet",
            summary: "Read a GeoParquet file and convert it to another vector format (or store it in memory).",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input GeoParquet file path (e.g. /work/data.parquet).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path; format is taken from its extension (.geojson, .fgb, .shp, .gpkg, ...). If omitted, the layer is stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;

        ctx.progress.info("reading GeoParquet");
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

/// Fetches a required, non-empty string parameter.
fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    match args.get(key).and_then(Value::as_str) {
        Some(s) if !s.trim().is_empty() => Ok(s),
        _ => Err(ToolError::Validation(format!(
            "missing required string parameter '{key}'"
        ))),
    }
}
