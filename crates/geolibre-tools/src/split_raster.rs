//! GeoLibre tool: split one raster into multiple standalone raster tiles.
//!
//! Pure-Rust counterpart of ArcGIS Data Management's *Split Raster*. Unlike
//! `raster_to_tiles` / `write_pmtiles` (which build web-tile pyramids), this
//! writes a set of independent raster files — an `N×M` grid, fixed-pixel tiles,
//! or one clip per polygon — a common preprocessing step for ML chipping,
//! parallel processing, and data delivery.
//!
//! Tile geometry (CRS, cell size, and the per-tile origin) is preserved so each
//! output georeferences correctly. An optional pixel `overlap` pads every tile.
//! In polygon mode, cells whose centers fall outside the clipping polygon are
//! set to no-data.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{Raster, RasterConfig, RasterFormat};

use crate::common::load_input_raster;
use crate::vector_common::{geometry_contains_point, load_input_layer, parse_optional_str};

pub struct SplitRasterTool;

impl Tool for SplitRasterTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "split_raster",
            display_name: "Split Raster",
            summary: "Split a raster into standalone tiles by tile count, fixed pixel size, or polygon features, with optional pixel overlap — like ArcGIS Split Raster. Each tile preserves CRS and georeferencing.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output_dir",
                    description: "Output directory for the tile files (created if needed).",
                    required: true,
                },
                ToolParamSpec {
                    name: "base_name",
                    description: "Output filename prefix (default 'tile').",
                    required: false,
                },
                ToolParamSpec {
                    name: "split_method",
                    description: "'count' (N×M grid, default), 'size' (fixed-pixel tiles), or 'polygons' (one clip per feature).",
                    required: false,
                },
                ToolParamSpec {
                    name: "num_x",
                    description: "Number of tile columns for 'count' (default 2).",
                    required: false,
                },
                ToolParamSpec {
                    name: "num_y",
                    description: "Number of tile rows for 'count' (default 2).",
                    required: false,
                },
                ToolParamSpec {
                    name: "tile_size_x",
                    description: "Tile width in pixels for 'size' (default 256).",
                    required: false,
                },
                ToolParamSpec {
                    name: "tile_size_y",
                    description: "Tile height in pixels for 'size' (default 256).",
                    required: false,
                },
                ToolParamSpec {
                    name: "polygons",
                    description: "Polygon vector layer for 'polygons' mode (one tile per feature).",
                    required: false,
                },
                ToolParamSpec {
                    name: "overlap",
                    description: "Pixel overlap padded onto each tile on all sides (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "format",
                    description: "Output file extension/driver: 'tif' (default) or 'png'.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "output_dir")?;
        let method = parse_method(args)?;
        if matches!(method, Method::Polygons) && parse_optional_str(args, "polygons")?.is_none() {
            return Err(ToolError::Validation(
                "split_method 'polygons' requires a 'polygons' layer".to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output_dir = require_str(args, "output_dir")?;
        let base_name = parse_optional_str(args, "base_name")?.unwrap_or("tile");
        let method = parse_method(args)?;
        let overlap = parse_usize(args, "overlap")?.unwrap_or(0) as isize;
        let ext = match parse_optional_str(args, "format")?.map(str::to_ascii_lowercase) {
            None => "tif".to_string(),
            Some(e) if e == "tif" || e == "tiff" || e == "png" => e,
            Some(o) => {
                return Err(ToolError::Validation(format!(
                    "'format' must be 'tif' or 'png', got '{o}'"
                )))
            }
        };

        let raster = load_input_raster(input)?;
        let rows = raster.rows as isize;
        let cols = raster.cols as isize;
        std::fs::create_dir_all(output_dir)
            .map_err(|e| ToolError::Execution(format!("failed creating output dir: {e}")))?;

        // Build the list of (name, col0, row0, col1, row1, optional clip geometry).
        let mut windows: Vec<Window> = Vec::new();
        match method {
            Method::Count => {
                let nx = parse_usize(args, "num_x")?.unwrap_or(2).max(1) as isize;
                let ny = parse_usize(args, "num_y")?.unwrap_or(2).max(1) as isize;
                for iy in 0..ny {
                    for ix in 0..nx {
                        let c0 = ix * cols / nx;
                        let c1 = (ix + 1) * cols / nx;
                        let r0 = iy * rows / ny;
                        let r1 = (iy + 1) * rows / ny;
                        if c1 > c0 && r1 > r0 {
                            windows.push(Window {
                                name: format!("{base_name}_r{iy}_c{ix}.{ext}"),
                                c0,
                                r0,
                                c1,
                                r1,
                                clip: None,
                            });
                        }
                    }
                }
            }
            Method::Size => {
                let tx = parse_usize(args, "tile_size_x")?.unwrap_or(256).max(1) as isize;
                let ty = parse_usize(args, "tile_size_y")?.unwrap_or(256).max(1) as isize;
                let mut iy = 0;
                let mut r0 = 0;
                while r0 < rows {
                    let r1 = (r0 + ty).min(rows);
                    let mut ix = 0;
                    let mut c0 = 0;
                    while c0 < cols {
                        let c1 = (c0 + tx).min(cols);
                        windows.push(Window {
                            name: format!("{base_name}_r{iy}_c{ix}.{ext}"),
                            c0,
                            r0,
                            c1,
                            r1,
                            clip: None,
                        });
                        c0 = c1;
                        ix += 1;
                    }
                    r0 = r1;
                    iy += 1;
                }
            }
            Method::Polygons => {
                let path = parse_optional_str(args, "polygons")?.unwrap();
                let layer = load_input_layer(path)?;
                for (i, f) in layer.features.iter().enumerate() {
                    let Some(geom) = f.geometry.as_ref() else {
                        continue;
                    };
                    let Some(bb) = geom.bbox() else { continue };
                    // Map world bbox to the cell window.
                    let (c0, c1) = col_range(&raster, bb.min_x, bb.max_x);
                    let (r0, r1) = row_range(&raster, bb.min_y, bb.max_y);
                    if c1 > c0 && r1 > r0 {
                        windows.push(Window {
                            name: format!("{base_name}_{i}.{ext}"),
                            c0,
                            r0,
                            c1,
                            r1,
                            clip: Some(geom.clone()),
                        });
                    }
                }
            }
        }

        if windows.is_empty() {
            return Err(ToolError::Execution(
                "no tiles produced (check split parameters / polygon overlap)".to_string(),
            ));
        }

        ctx.progress
            .info(&format!("writing {} tile(s)", windows.len()));
        let mut written = Vec::new();
        let n = windows.len();
        for (idx, w) in windows.iter().enumerate() {
            // Apply overlap padding, clamped to the raster.
            let c0 = (w.c0 - overlap).max(0);
            let r0 = (w.r0 - overlap).max(0);
            let c1 = (w.c1 + overlap).min(cols);
            let r1 = (w.r1 + overlap).min(rows);
            let path = format!("{}/{}", output_dir.trim_end_matches('/'), w.name);
            write_tile(&raster, c0, r0, c1, r1, w.clip.as_ref(), &path)?;
            written.push(path);
            ctx.progress.progress((idx as f64 + 1.0) / n as f64);
        }

        let mut outputs = BTreeMap::new();
        outputs.insert("output_dir".to_string(), json!(output_dir));
        outputs.insert("tile_count".to_string(), json!(written.len()));
        outputs.insert("tiles".to_string(), json!(written));
        Ok(ToolRunResult { outputs })
    }
}

struct Window {
    name: String,
    c0: isize,
    r0: isize,
    c1: isize,
    r1: isize,
    clip: Option<wbvector::Geometry>,
}

/// Writes the cell window `[c0,c1) × [r0,r1)` of `src` to `path`, preserving
/// CRS, cell size, and the tile's georeferenced origin. When `clip` is given,
/// cells whose centers fall outside it become no-data.
fn write_tile(
    src: &Raster,
    c0: isize,
    r0: isize,
    c1: isize,
    r1: isize,
    clip: Option<&wbvector::Geometry>,
    path: &str,
) -> Result<(), ToolError> {
    let tcols = (c1 - c0) as usize;
    let trows = (r1 - r0) as usize;
    let cell_x = src.cell_size_x.abs();
    let cell_y = src.cell_size_y.abs();
    // Tile south edge: source y_min + (source_rows - r1) rows of cell_y.
    let tile_y_min = src.y_min + (src.rows as isize - r1) as f64 * cell_y;
    let tile_x_min = src.x_min + c0 as f64 * cell_x;

    let mut tile = Raster::new(RasterConfig {
        cols: tcols,
        rows: trows,
        bands: src.bands,
        x_min: tile_x_min,
        y_min: tile_y_min,
        cell_size: cell_x,
        cell_size_y: Some(cell_y),
        nodata: src.nodata,
        data_type: src.data_type,
        crs: src.crs.clone(),
        metadata: src.metadata.clone(),
    });

    for band in 0..src.bands as isize {
        for tr in 0..trows {
            let sr = r0 + tr as isize;
            for tc in 0..tcols {
                let sc = c0 + tc as isize;
                let mut v = src.get(band, sr, sc);
                if let Some(g) = clip {
                    let (x, y) = cell_center(src, sr, sc);
                    if !geometry_contains_point(g, x, y) {
                        v = src.nodata;
                    }
                }
                tile.set(band, tr as isize, tc as isize, v)
                    .map_err(|e| ToolError::Execution(format!("failed writing tile cell: {e}")))?;
            }
        }
    }

    let fmt = RasterFormat::for_output_path(path)
        .map_err(|e| ToolError::Validation(format!("unsupported output path: {e}")))?;
    tile.write(path, fmt)
        .map_err(|e| ToolError::Execution(format!("failed writing tile {path}: {e}")))
}

/// Cell-center world coordinates for a source cell.
fn cell_center(r: &Raster, row: isize, col: isize) -> (f64, f64) {
    let x = r.x_min + (col as f64 + 0.5) * r.cell_size_x.abs();
    let y = r.y_min + (r.rows as f64 - 1.0 - row as f64 + 0.5) * r.cell_size_y.abs();
    (x, y)
}

/// Clamps a world x-range to source column indices `[c0, c1)`.
fn col_range(r: &Raster, min_x: f64, max_x: f64) -> (isize, isize) {
    let cx = r.cell_size_x.abs();
    let cols = r.cols as isize;
    let a = ((min_x - r.x_min) / cx).floor() as isize;
    let b = ((max_x - r.x_min) / cx).ceil() as isize;
    (a.max(0).min(cols), b.max(0).min(cols))
}

/// Clamps a world y-range to source row indices `[r0, r1)` (row 0 = top).
fn row_range(r: &Raster, min_y: f64, max_y: f64) -> (isize, isize) {
    let cy = r.cell_size_y.abs();
    let rows = r.rows as isize;
    // world y -> row: row = rows - (y - y_min)/cy   (top-down)
    let top = (r.rows as f64 - (max_y - r.y_min) / cy).floor() as isize;
    let bot = (r.rows as f64 - (min_y - r.y_min) / cy).ceil() as isize;
    (top.max(0).min(rows), bot.max(0).min(rows))
}

enum Method {
    Count,
    Size,
    Polygons,
}

fn parse_method(args: &ToolArgs) -> Result<Method, ToolError> {
    match args
        .get("split_method")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some("count") => Ok(Method::Count),
        Some("size") => Ok(Method::Size),
        Some("polygons") => Ok(Method::Polygons),
        Some(o) => Err(ToolError::Validation(format!(
            "'split_method' must be 'count', 'size', or 'polygons', got '{o}'"
        ))),
    }
}

fn parse_usize(args: &ToolArgs, key: &str) -> Result<Option<usize>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => n
            .as_i64()
            .filter(|v| *v >= 0)
            .map(|v| Some(v as usize))
            .ok_or_else(|| {
                ToolError::Validation(format!("'{key}' must be a non-negative integer"))
            }),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => {
            s.trim().parse::<usize>().map(Some).map_err(|_| {
                ToolError::Validation(format!("'{key}' must be a non-negative integer"))
            })
        }
        Some(_) => Err(ToolError::Validation(format!(
            "'{key}' must be a non-negative integer"
        ))),
    }
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{DataType, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn raster_path(rows: usize, cols: usize) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: Default::default(),
            metadata: Default::default(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, (row * cols + col) as f64)
                    .unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn tmp_dir(tag: &str) -> String {
        use std::sync::atomic::{AtomicU64, Ordering};
        static SEQ: AtomicU64 = AtomicU64::new(0);
        let n = SEQ.fetch_add(1, Ordering::Relaxed);
        let d = std::env::temp_dir().join(format!("split_raster_{tag}_{}_{n}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        d.to_str().unwrap().to_string()
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        SplitRasterTool.run(&args, &ctx()).unwrap()
    }

    #[test]
    fn count_mode_makes_grid_of_tiles() {
        let dir = tmp_dir("count");
        let out = run(json!({
            "input": raster_path(4, 4), "output_dir": dir,
            "split_method": "count", "num_x": 2, "num_y": 2
        }));
        assert_eq!(out.outputs["tile_count"], json!(4));
        // Each 2x2 tile georeferences and re-reads.
        let first = out.outputs["tiles"].as_array().unwrap()[0]
            .as_str()
            .unwrap();
        let t = load_input_raster(first).unwrap();
        assert_eq!(t.cols, 2);
        assert_eq!(t.rows, 2);
    }

    #[test]
    fn size_mode_covers_all_pixels() {
        let dir = tmp_dir("size");
        let out = run(json!({
            "input": raster_path(5, 5), "output_dir": dir,
            "split_method": "size", "tile_size_x": 2, "tile_size_y": 2
        }));
        // 5/2 -> 3 tiles per axis -> 9 tiles.
        assert_eq!(out.outputs["tile_count"], json!(9));
    }

    #[test]
    fn count_tile_preserves_origin() {
        let dir = tmp_dir("origin");
        let out = run(json!({
            "input": raster_path(4, 4), "output_dir": dir,
            "split_method": "count", "num_x": 2, "num_y": 1
        }));
        // The east tile's x_min must be shifted by 2 cells.
        let tiles = out.outputs["tiles"].as_array().unwrap();
        let east = load_input_raster(tiles[1].as_str().unwrap()).unwrap();
        assert!((east.x_min - 2.0).abs() < 1e-9, "east x_min {}", east.x_min);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SplitRasterTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "input": "x.tif" })).is_err(),
            "missing output_dir"
        );
        assert!(
            bad(json!({ "input": "x.tif", "output_dir": "/d", "split_method": "polygons" }))
                .is_err(),
            "polygons needs layer"
        );
        assert!(bad(json!({ "input": "x.tif", "output_dir": "/d" })).is_ok());
    }
}
