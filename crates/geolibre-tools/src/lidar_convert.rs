//! GeoLibre tool: convert a LiDAR point cloud between LAS, LAZ and COPC.
//!
//! The output format follows the output path's extension (`.las`, `.laz`, or
//! `.copc.laz`), so one tool covers every direction: LAS -> LAZ compression,
//! LAZ -> LAS decompression, and building a Cloud Optimized Point Cloud
//! (an octree-ordered LAZ that viewers can stream) from any of them. GeoLibre's
//! layer menu uses it to export a loaded point cloud in the format the user
//! picks.
//!
//! The cloud is read whole (a COPC writer needs every point to build its
//! octree), so on 32-bit targets it refuses inputs whose header declares more
//! points than fit in memory, before decoding anything.

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wblidar::memory_store;

use crate::lidar_common::{load_input_cloud, write_or_store_cloud};

/// Cap on input points on 32-bit targets (wasm32). A decoded point is 336
/// bytes and a wasm32 allocation tops out at 2 GB, so 5M points (1.7 GB) is the
/// most one `Vec<PointRecord>` can hold with room left for the writer.
#[cfg(target_pointer_width = "32")]
const MAX_POINTS: u64 = 5_000_000;
#[cfg(not(target_pointer_width = "32"))]
const MAX_POINTS: u64 = u64::MAX;

const OUTPUT_EXTENSIONS: [&str; 4] = [".las", ".laz", ".copc.laz", ".copc.las"];

pub struct LidarConvertTool;

impl Tool for LidarConvertTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "lidar_convert",
            display_name: "LiDAR Convert",
            summary: "Convert a LiDAR point cloud between LAS, LAZ and COPC. The output format follows the output file extension: .las (uncompressed), .laz (LASzip compressed) or .copc.laz (Cloud Optimized Point Cloud, streamable by web viewers).",
            category: ToolCategory::Lidar,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input LAS/LAZ/COPC point cloud.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output point cloud path; .las, .laz or .copc.laz picks the format. If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        input_path(args)?;
        output_path(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = input_path(args)?;
        let output = output_path(args)?;

        if !memory_store::lidar_is_memory_path(input) {
            // Header-only; refuses a cloud too large to decode before trying.
            if let Ok(count) = wblidar::frontend::read_point_count(input) {
                if count > MAX_POINTS {
                    return Err(ToolError::Execution(format!(
                        "{count} points is more than this build can convert in memory (at most {MAX_POINTS}); thin the cloud first, e.g. with lidar_grid_thin"
                    )));
                }
            }
        }

        let cloud = load_input_cloud(input)?;
        if cloud.points.is_empty() {
            return Err(ToolError::Execution(
                "input point cloud has no points".to_string(),
            ));
        }
        let points = cloud.points.len();
        if cloud.crs.is_none() {
            ctx.progress
                .info("input has no readable CRS metadata; the output will have none");
        }
        let out_path = write_or_store_cloud(cloud, output)?;
        ctx.progress.info(&format!("wrote {points} points"));

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("points".to_string(), json!(points));
        Ok(ToolRunResult { outputs })
    }
}

fn input_path(args: &ToolArgs) -> Result<&str, ToolError> {
    args.get("input")
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation("missing required parameter 'input'".to_string()))
}

/// The output path, if given, after checking its extension names a format.
fn output_path(args: &ToolArgs) -> Result<Option<&str>, ToolError> {
    let output = crate::common::parse_optional_output(args, "output")?;
    if let Some(path) = output {
        let lower = path.to_ascii_lowercase();
        if !OUTPUT_EXTENSIONS.iter().any(|ext| lower.ends_with(ext)) {
            return Err(ToolError::Validation(
                "'output' must end in .las, .laz or .copc.laz".to_string(),
            ));
        }
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wblidar::{PointCloud, PointRecord};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn cloud_path(n: usize) -> String {
        let mut cloud = PointCloud::default();
        cloud.points = (0..n)
            .map(|i| {
                let mut p = PointRecord::default();
                p.x = 500_000.0 + i as f64;
                p.y = 4_000_000.0 + (i % 7) as f64;
                p.z = 100.0 + (i % 13) as f64;
                p.classification = 2;
                p
            })
            .collect();
        memory_store::make_lidar_memory_path(&memory_store::put_lidar(cloud))
    }

    fn convert(input: &str, output: &str) -> usize {
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "output": output })).unwrap();
        let out = LidarConvertTool.run(&args, &ctx()).unwrap();
        out.outputs["points"].as_u64().unwrap() as usize
    }

    /// LAS -> LAZ -> COPC -> LAS keeps every point, and each file reads back.
    #[test]
    fn round_trips_through_every_format() {
        let dir = std::env::temp_dir().join(format!("lidar_convert_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let las = dir.join("a.las").to_string_lossy().to_string();
        let laz = dir.join("b.laz").to_string_lossy().to_string();
        let copc = dir.join("c.copc.laz").to_string_lossy().to_string();
        let back = dir.join("d.las").to_string_lossy().to_string();

        assert_eq!(convert(&cloud_path(500), &las), 500);
        assert_eq!(convert(&las, &laz), 500);
        assert_eq!(convert(&laz, &copc), 500);
        assert_eq!(convert(&copc, &back), 500);

        let copc_bytes = std::fs::read(&copc).unwrap();
        assert!(
            copc_bytes.windows(4).any(|w| w == b"copc"),
            "COPC output has no copc VLR"
        );
        let cloud = load_input_cloud(&back).unwrap();
        assert_eq!(cloud.points.len(), 500);
        assert!(cloud.points.iter().all(|p| p.classification == 2));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn rejects_an_unknown_output_extension() {
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": cloud_path(3), "output": "x.ply" })).unwrap();
        assert!(LidarConvertTool.validate(&args).is_err());
    }
}
