//! Shared vector I/O for GeoLibre tools: load a `Layer` from a file path or an
//! in-memory (`memory://`) handle, and write it back to a path or memory.
//!
//! Mirrors `common.rs` (the raster equivalent) but for `wbvector::Layer`.

use std::sync::Arc;

use serde_json::Value;
use wbcore::{ToolArgs, ToolError};
use wbvector::{memory_store, Coord, Geometry, Layer, VectorFormat};

/// Even-odd ray-casting test: is `(x, y)` inside the closed ring `coords`?
/// Points exactly on an edge are treated as inside (deterministic tie-break).
pub fn ring_contains(coords: &[Coord], x: f64, y: f64) -> bool {
    let n = coords.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = (coords[i].x, coords[i].y);
        let (xj, yj) = (coords[j].x, coords[j].y);
        if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi) {
            inside = !inside;
        }
        j = i;
    }
    inside
}

/// Tests whether `(x, y)` falls inside a polygon (inside its exterior and
/// outside every hole).
pub fn polygon_contains(exterior: &[Coord], interiors: &[Vec<Coord>], x: f64, y: f64) -> bool {
    if !ring_contains(exterior, x, y) {
        return false;
    }
    for hole in interiors {
        if ring_contains(hole, x, y) {
            return false;
        }
    }
    true
}

/// Tests whether `(x, y)` falls inside a `Polygon`/`MultiPolygon` geometry.
/// Non-areal geometries always return `false`.
pub fn geometry_contains_point(geom: &Geometry, x: f64, y: f64) -> bool {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            let holes: Vec<Vec<Coord>> = interiors.iter().map(|r| r.0.clone()).collect();
            polygon_contains(&exterior.0, &holes, x, y)
        }
        Geometry::MultiPolygon(polys) => polys.iter().any(|(ext, hs)| {
            let holes: Vec<Vec<Coord>> = hs.iter().map(|r| r.0.clone()).collect();
            polygon_contains(&ext.0, &holes, x, y)
        }),
        Geometry::GeometryCollection(gs) => gs.iter().any(|g| geometry_contains_point(g, x, y)),
        _ => false,
    }
}

/// Parses an optional string parameter (absent / null / empty -> None).
pub fn parse_optional_str<'a>(args: &'a ToolArgs, key: &str) -> Result<Option<&'a str>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => Ok(Some(s.as_str())),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a string when provided"
        ))),
    }
}

/// Loads a vector layer from a file path (format auto-detected) or an in-memory
/// (`memory://`) handle.
pub fn load_input_layer(path: &str) -> Result<Layer, ToolError> {
    if memory_store::vector_is_memory_path(path) {
        let id = memory_store::vector_path_to_id(path)
            .ok_or_else(|| ToolError::Validation("malformed in-memory vector path".to_string()))?;
        let arc: Arc<Layer> = memory_store::get_vector_arc_by_id(id)
            .ok_or_else(|| ToolError::Validation(format!("unknown in-memory vector id '{id}'")))?;
        return Ok((*arc).clone());
    }
    wbvector::read(path)
        .map_err(|e| ToolError::Execution(format!("failed reading input vector: {e}")))
}

/// Writes `layer` to `output_path` using the format implied by its extension,
/// or stores it in memory and returns a `memory://` handle when no path is given.
pub fn write_or_store_layer(layer: Layer, output_path: Option<&str>) -> Result<String, ToolError> {
    match output_path {
        Some(path) => {
            let fmt = VectorFormat::detect(path)
                .map_err(|e| ToolError::Validation(format!("unsupported output path: {e}")))?;
            ensure_parent_dir(path)?;
            wbvector::write(&layer, path, fmt)
                .map_err(|e| ToolError::Execution(format!("failed writing output vector: {e}")))?;
            Ok(path.to_string())
        }
        None => {
            let id = memory_store::put_vector(layer);
            Ok(memory_store::make_vector_memory_path(&id))
        }
    }
}

/// Creates the parent directory of `path` if needed.
pub fn ensure_parent_dir(path: &str) -> Result<(), ToolError> {
    if let Some(parent) = std::path::Path::new(path).parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|e| {
                ToolError::Execution(format!("failed creating output directory: {e}"))
            })?;
        }
    }
    Ok(())
}
