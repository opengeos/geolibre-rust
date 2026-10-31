//! GeoLibre tool: extract a bbox/zoom subset of a PMTiles archive.
//!
//! Drives the sans-IO [`geolibre_pmtiles::extract::Extractor`] with plain file
//! reads, so it works both natively (a downloaded planet build on disk) and in
//! the WASI runner over `/work`. Remote extraction over HTTP uses the same
//! engine through the `geolibre-wasm` `PmtilesExtractor` binding instead.

use std::fs::File;
use std::io::{Read, Seek, SeekFrom};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata,
    ToolParamSpec, ToolRunResult,
};

use geolibre_pmtiles::extract::{ExtractOptions, Extractor};
use geolibre_pmtiles::LonLatBounds;

/// Extracts the tiles intersecting a bbox across a zoom range into a new,
/// self-contained PMTiles archive.
pub struct PmtilesExtractTool;

impl Tool for PmtilesExtractTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "pmtiles_extract",
            display_name: "PMTiles Extract",
            summary: "Extract a bbox/zoom subset of a PMTiles archive into a new archive (e.g. an offline basemap from a Protomaps planet build).",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input .pmtiles file path.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output .pmtiles file path.",
                    required: true,
                },
                ToolParamSpec {
                    name: "bbox",
                    description: "WGS84 bounding box as 'min_lon,min_lat,max_lon,max_lat'.",
                    required: true,
                },
                ToolParamSpec {
                    name: "min_zoom",
                    description: "Lowest zoom level to include (default 0, so the map stays usable zoomed out).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_zoom",
                    description: "Highest zoom level to include (default: the source archive's max zoom).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_tiles",
                    description: "Abort if the selection addresses more tiles than this (default 2000000).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["input", "output"] {
            if args.get(key).and_then(Value::as_str).map(str::trim).unwrap_or("").is_empty() {
                return Err(ToolError::Validation(format!(
                    "missing required string parameter '{key}'"
                )));
            }
        }
        parse_bbox(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args.get("input").and_then(Value::as_str).ok_or_else(|| {
            ToolError::Validation("missing required parameter 'input'".to_string())
        })?;
        let output = args.get("output").and_then(Value::as_str).ok_or_else(|| {
            ToolError::Validation("missing required parameter 'output'".to_string())
        })?;
        let bbox = parse_bbox(args)?;
        let min_zoom = optional_u8(args, "min_zoom")?.unwrap_or(0);
        // 30 exceeds any real archive; the extractor clamps to the source.
        let max_zoom = optional_u8(args, "max_zoom")?.unwrap_or(30);

        let mut opts = ExtractOptions::new(bbox, min_zoom, max_zoom);
        if let Some(cap) = args.get("max_tiles").and_then(Value::as_u64) {
            opts.max_tiles = cap;
        }

        let mut file = File::open(input)
            .map_err(|e| ToolError::Execution(format!("cannot open {input}: {e}")))?;
        let file_len = file
            .metadata()
            .map_err(|e| ToolError::Execution(format!("cannot stat {input}: {e}")))?
            .len();

        let mut extractor = Extractor::new(opts)
            .map_err(|e| ToolError::Execution(e.to_string()))?;
        ctx.progress.info("reading header and directories");
        let mut reported_data = false;
        while !extractor.is_done() {
            let wants = extractor.wanted();
            if wants.is_empty() {
                return Err(ToolError::Execution(
                    "extractor stalled: nothing wanted but not done".to_string(),
                ));
            }
            for range in wants {
                if range.offset >= file_len {
                    return Err(ToolError::Execution(format!(
                        "archive truncated: needs bytes at {} but file is {} bytes",
                        range.offset, file_len
                    )));
                }
                let len = range.length.min(file_len - range.offset) as usize;
                let mut buf = vec![0u8; len];
                file.seek(SeekFrom::Start(range.offset))
                    .and_then(|_| file.read_exact(&mut buf))
                    .map_err(|e| ToolError::Execution(format!("read {input}: {e}")))?;
                extractor
                    .feed(range.offset, &buf)
                    .map_err(|e| ToolError::Execution(e.to_string()))?;
            }
            let p = extractor.progress();
            if p.phase == "data" && !reported_data {
                reported_data = true;
                ctx.progress.info(&format!(
                    "copying {} tiles ({} blobs, {} bytes)",
                    p.tiles_selected, p.blobs_total, p.data_bytes_total
                ));
            }
        }

        let progress = extractor.progress();
        let archive = extractor
            .finish()
            .map_err(|e| ToolError::Execution(e.to_string()))?;
        let archive_len = archive.len();
        std::fs::write(output, archive)
            .map_err(|e| ToolError::Execution(format!("write {output}: {e}")))?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(output));
        outputs.insert("tiles".to_string(), json!(progress.tiles_selected));
        outputs.insert("bytes".to_string(), json!(archive_len));
        Ok(ToolRunResult { outputs })
    }
}

/// Parses `bbox` = "min_lon,min_lat,max_lon,max_lat".
fn parse_bbox(args: &ToolArgs) -> Result<LonLatBounds, ToolError> {
    let raw = args
        .get("bbox")
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| {
            ToolError::Validation("missing required string parameter 'bbox'".to_string())
        })?;
    let parts: Vec<f64> = raw
        .split(',')
        .map(|p| p.trim().parse::<f64>())
        .collect::<Result<_, _>>()
        .map_err(|_| {
            ToolError::Validation(format!(
                "bbox must be 'min_lon,min_lat,max_lon,max_lat', got '{raw}'"
            ))
        })?;
    if parts.len() != 4 {
        return Err(ToolError::Validation(format!(
            "bbox must have 4 comma-separated numbers, got {}",
            parts.len()
        )));
    }
    Ok(LonLatBounds {
        min_lon: parts[0],
        min_lat: parts[1],
        max_lon: parts[2],
        max_lat: parts[3],
    })
}

fn optional_u8(args: &ToolArgs, key: &str) -> Result<Option<u8>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(v) => {
            let n = v
                .as_u64()
                .or_else(|| v.as_str().and_then(|s| s.trim().parse::<u64>().ok()))
                .filter(|&n| n <= 30)
                .ok_or_else(|| {
                    ToolError::Validation(format!("'{key}' must be an integer in 0..=30"))
                })?;
            Ok(Some(n as u8))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use geolibre_pmtiles::writer::{build_png, Tile};
    use wbcore::{AllowAllCapabilities, ProgressSink};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn args(v: Value) -> ToolArgs {
        serde_json::from_value(v).unwrap()
    }

    #[test]
    fn validate_catches_missing_and_malformed_params() {
        let tool = PmtilesExtractTool;
        assert!(tool.validate(&args(json!({}))).is_err());
        assert!(tool
            .validate(&args(json!({
                "input": "a.pmtiles", "output": "b.pmtiles", "bbox": "1,2,3"
            })))
            .is_err());
        assert!(tool
            .validate(&args(json!({
                "input": "a.pmtiles", "output": "b.pmtiles", "bbox": "5,5,60,55"
            })))
            .is_ok());
    }

    #[test]
    fn extracts_from_a_file_end_to_end() {
        let dir = std::env::temp_dir().join(format!("pmx-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let input = dir.join("src.pmtiles");
        let output = dir.join("out.pmtiles");

        let mut tiles = Vec::new();
        for z in 0..=3u8 {
            let n = 1u32 << z;
            for y in 0..n {
                for x in 0..n {
                    tiles.push(Tile { z, x, y, data: format!("t{z}/{x}/{y}").into_bytes() });
                }
            }
        }
        let world =
            LonLatBounds { min_lon: -180.0, min_lat: -85.0, max_lon: 180.0, max_lat: 85.0 };
        std::fs::write(&input, build_png(tiles, &world, 0, 3).unwrap()).unwrap();

        let tool = PmtilesExtractTool;
        let a = args(json!({
            "input": input.to_str().unwrap(),
            "output": output.to_str().unwrap(),
            "bbox": "5,5,60,55",
        }));
        let result = tool.run(&a, &ctx()).unwrap();
        assert!(result.outputs["tiles"].as_u64().unwrap() > 0);

        let out = std::fs::read(&output).unwrap();
        let h = geolibre_pmtiles::format::Header::parse(&out).unwrap();
        assert_eq!((h.min_zoom, h.max_zoom), (0, 3));
        std::fs::remove_dir_all(&dir).ok();
    }
}
