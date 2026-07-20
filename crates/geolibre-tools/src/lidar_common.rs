//! Shared LiDAR I/O for GeoLibre tools: load a `PointCloud` from a file path or
//! an in-memory (`memory://`) handle, and write it back to a path or memory.
//!
//! Mirrors `common.rs` (raster) and `vector_common.rs` (vector) but for
//! `wblidar::PointCloud`.

use std::sync::Arc;

use wbcore::ToolError;
use wblidar::{memory_store, PointCloud};

/// Loads a point cloud from a file path (LAS/LAZ/COPC, format auto-detected) or
/// an in-memory (`memory://`) handle.
pub fn load_input_cloud(path: &str) -> Result<PointCloud, ToolError> {
    if memory_store::lidar_is_memory_path(path) {
        let id = memory_store::lidar_path_to_id(path)
            .ok_or_else(|| ToolError::Validation("malformed in-memory lidar path".to_string()))?;
        let arc: Arc<PointCloud> = memory_store::get_lidar_arc_by_id(id)
            .ok_or_else(|| ToolError::Validation(format!("unknown in-memory lidar id '{id}'")))?;
        return Ok((*arc).clone());
    }
    PointCloud::read(path)
        .map_err(|e| ToolError::Execution(format!("failed reading input lidar: {e}")))
}

/// Writes `cloud` to `output_path` using the format implied by its extension, or
/// stores it in memory and returns a `memory://` handle when no path is given.
pub fn write_or_store_cloud(
    cloud: PointCloud,
    output_path: Option<&str>,
) -> Result<String, ToolError> {
    match output_path {
        Some(path) => {
            if let Some(parent) = std::path::Path::new(path).parent() {
                if !parent.as_os_str().is_empty() {
                    std::fs::create_dir_all(parent).map_err(|e| {
                        ToolError::Execution(format!("failed creating output directory: {e}"))
                    })?;
                }
            }
            cloud
                .write(path)
                .map_err(|e| ToolError::Execution(format!("failed writing output lidar: {e}")))?;
            Ok(path.to_string())
        }
        None => {
            let id = memory_store::put_lidar(cloud);
            Ok(memory_store::make_lidar_memory_path(&id))
        }
    }
}
