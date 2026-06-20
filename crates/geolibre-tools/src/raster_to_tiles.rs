//! GeoLibre tool: slice a raster into a Web Mercator XYZ tile pyramid of PNGs.
//!
//! The suite has no web-map tiling output. This reprojects the input to
//! EPSG:3857, then for each requested zoom level samples the standard slippy-map
//! tile grid into 256x256 RGBA PNGs written as `{output_dir}/{z}/{x}/{y}.png`.
//! Fully transparent (no-data) tiles are skipped so the pyramid stays sparse.

use std::path::Path;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata,
    ToolParamSpec, ToolRunResult,
};
use wbraster::{NodataPolicy, Raster, ResampleMethod};

use crate::common::load_input_raster;
use crate::render::{encode_png_rgba, normalize, Colormap};
use crate::reproject_raster::parse_resample;

/// EPSG:3857 (Web Mercator) world half-extent in meters: pi * 6378137.
const ORIGIN: f64 = 20_037_508.342_789_244;
const TILE_SIZE: usize = 256;
const WEB_MERCATOR_EPSG: u32 = 3857;
/// Safety cap so a too-wide zoom range cannot generate an unbounded pyramid.
const MAX_TILES: usize = 4096;

/// Renders a raster into a Web Mercator XYZ PNG tile pyramid.
pub struct RasterToTilesTool;

impl Tool for RasterToTilesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "raster_to_tiles",
            display_name: "Raster to XYZ Tiles",
            summary: "Slice a raster into a Web Mercator (EPSG:3857) XYZ PNG tile pyramid for web maps.",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input raster file path. Must carry a source CRS (EPSG or WKT).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output_dir",
                    description: "Output directory for the {z}/{x}/{y}.png tile tree (e.g. /work/tiles).",
                    required: true,
                },
                ToolParamSpec {
                    name: "min_zoom",
                    description: "Minimum zoom level (default: a single native zoom matching the raster resolution).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_zoom",
                    description: "Maximum zoom level (default: same as min_zoom).",
                    required: false,
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
                    name: "method",
                    description: "Resampling method for reprojection and sampling: bilinear (default), nearest, cubic.",
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
        if args.get("output_dir").and_then(Value::as_str).is_none() {
            return Err(ToolError::Validation(
                "missing required string parameter 'output_dir'".to_string(),
            ));
        }
        if let Some(c) = args.get("colormap").and_then(Value::as_str) {
            Colormap::parse(c)?;
        }
        if let Some(m) = args.get("method").and_then(Value::as_str) {
            parse_resample(m)?;
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args
            .get("input")
            .and_then(Value::as_str)
            .ok_or_else(|| ToolError::Validation("missing required parameter 'input'".to_string()))?;
        let output_dir = args
            .get("output_dir")
            .and_then(Value::as_str)
            .ok_or_else(|| {
                ToolError::Validation("missing required parameter 'output_dir'".to_string())
            })?
            .trim_end_matches('/');
        let band_1based = args.get("band").and_then(Value::as_u64).unwrap_or(1).max(1);
        let band = (band_1based - 1) as isize;
        let colormap = match args.get("colormap").and_then(Value::as_str) {
            Some(c) => Colormap::parse(c)?,
            None => Colormap::Viridis,
        };
        let method = match args.get("method").and_then(Value::as_str) {
            Some(m) => parse_resample(m)?,
            None => ResampleMethod::Bilinear,
        };

        let raster = load_input_raster(input)?;
        if band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "band {band_1based} out of range (raster has {} band(s))",
                raster.bands
            )));
        }

        // Work in Web Mercator: reproject unless the source is already 3857.
        ctx.progress.info("reprojecting to EPSG:3857");
        let merc: Raster = if raster.crs.epsg == Some(WEB_MERCATOR_EPSG) {
            raster
        } else {
            if raster.crs.epsg.is_none()
                && raster.crs.wkt.is_none()
                && raster.crs.proj4.is_none()
            {
                return Err(ToolError::Validation(
                    "input raster has no source CRS (EPSG/WKT/PROJ); cannot tile".to_string(),
                ));
            }
            raster
                .reproject_to_epsg(WEB_MERCATOR_EPSG, method)
                .map_err(|e| ToolError::Execution(format!("reprojection to 3857 failed: {e}")))?
        };

        let stats = merc
            .statistics_band(band)
            .map_err(|e| ToolError::Execution(format!("failed computing band statistics: {e}")))?;
        let min = args.get("min").and_then(Value::as_f64).unwrap_or(stats.min);
        let max = args.get("max").and_then(Value::as_f64).unwrap_or(stats.max);
        if !min.is_finite() || !max.is_finite() {
            return Err(ToolError::Execution(
                "raster band has no finite values to render".to_string(),
            ));
        }

        // Resolve the zoom range. Default to a single native zoom that matches
        // the reprojected pixel size.
        let native = native_zoom(merc.cell_size_x.abs());
        let min_zoom = args
            .get("min_zoom")
            .and_then(Value::as_u64)
            .map(|z| z as u32)
            .unwrap_or(native);
        let max_zoom = args
            .get("max_zoom")
            .and_then(Value::as_u64)
            .map(|z| z as u32)
            .unwrap_or(min_zoom)
            .max(min_zoom);
        if max_zoom > 24 {
            return Err(ToolError::Validation(
                "max_zoom must be <= 24".to_string(),
            ));
        }

        let ex = merc.extent();
        let mut written = 0usize;
        let mut total = 0usize;
        let mut zoom_summaries = Vec::new();

        for z in min_zoom..=max_zoom {
            let n = 1u64 << z; // tiles per side
            let span = (2.0 * ORIGIN) / n as f64; // tile width in meters
            let (tx_min, tx_max) = tile_range(ex.x_min, ex.x_max, -ORIGIN, span, n);
            // Y tile index increases southward (row 0 at the top / north).
            let (ty_min, ty_max) = tile_range(ORIGIN - ex.y_max, ORIGIN - ex.y_min, 0.0, span, n);

            let count_this_zoom =
                ((tx_max - tx_min + 1) as usize) * ((ty_max - ty_min + 1) as usize);
            total += count_this_zoom;
            if total > MAX_TILES {
                return Err(ToolError::Validation(format!(
                    "tile pyramid would exceed {MAX_TILES} tiles; narrow the zoom range (min_zoom/max_zoom)"
                )));
            }

            let mut z_written = 0usize;
            for tx in tx_min..=tx_max {
                for ty in ty_min..=ty_max {
                    let x0 = -ORIGIN + tx as f64 * span;
                    let y_top = ORIGIN - ty as f64 * span;
                    if let Some(png) =
                        render_tile(&merc, band, x0, y_top, span, method, colormap, min, max)?
                    {
                        let path = format!("{output_dir}/{z}/{tx}/{ty}.png");
                        write_bytes(&path, &png)?;
                        z_written += 1;
                    }
                }
            }
            written += z_written;
            zoom_summaries.push(json!({ "zoom": z, "tiles_written": z_written }));
            ctx.progress
                .progress((z - min_zoom + 1) as f64 / (max_zoom - min_zoom + 1) as f64);
        }

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output_dir".to_string(), json!(output_dir));
        outputs.insert("min_zoom".to_string(), json!(min_zoom));
        outputs.insert("max_zoom".to_string(), json!(max_zoom));
        outputs.insert("tiles_written".to_string(), json!(written));
        outputs.insert("zooms".to_string(), json!(zoom_summaries));
        Ok(ToolRunResult { outputs })
    }
}

/// Renders one 256x256 tile. Returns `None` if every pixel is no-data / outside
/// the raster, so the caller can skip writing an empty tile.
#[allow(clippy::too_many_arguments)]
fn render_tile(
    merc: &Raster,
    band: isize,
    x0: f64,
    y_top: f64,
    span: f64,
    method: ResampleMethod,
    colormap: Colormap,
    min: f64,
    max: f64,
) -> Result<Option<Vec<u8>>, ToolError> {
    let px_span = span / TILE_SIZE as f64;
    let mut rgba = vec![0u8; TILE_SIZE * TILE_SIZE * 4];
    let mut any = false;
    for py in 0..TILE_SIZE {
        let wy = y_top - (py as f64 + 0.5) * px_span;
        for px in 0..TILE_SIZE {
            let wx = x0 + (px as f64 + 0.5) * px_span;
            if let Some(v) = merc.sample_world(band, wx, wy, method, NodataPolicy::Strict) {
                let idx = (py * TILE_SIZE + px) * 4;
                let [r, g, b] = colormap.rgb(normalize(v, min, max));
                rgba[idx] = r;
                rgba[idx + 1] = g;
                rgba[idx + 2] = b;
                rgba[idx + 3] = 255;
                any = true;
            }
        }
    }
    if !any {
        return Ok(None);
    }
    Ok(Some(encode_png_rgba(&rgba, TILE_SIZE as u32, TILE_SIZE as u32)?))
}

/// Inclusive tile-index range covering `[lo, hi]` along an axis whose tile 0
/// starts at `axis_origin`, with `span`-meter tiles and `n` tiles total.
fn tile_range(lo: f64, hi: f64, axis_origin: f64, span: f64, n: u64) -> (u64, u64) {
    let max_idx = n - 1;
    let a = (((lo - axis_origin) / span).floor()).clamp(0.0, max_idx as f64) as u64;
    // Subtract a tiny epsilon so a coordinate exactly on a tile boundary does
    // not pull in the next (empty) tile.
    let b = (((hi - axis_origin) / span - 1e-9).floor()).clamp(0.0, max_idx as f64) as u64;
    (a.min(b), a.max(b))
}

/// Zoom level whose 256px tiles have a pixel size closest to `cell_size_m`.
fn native_zoom(cell_size_m: f64) -> u32 {
    if cell_size_m <= 0.0 {
        return 0;
    }
    // res(z) = (2*ORIGIN) / (256 * 2^z); solve for z so res ~= cell_size_m.
    let z = ((2.0 * ORIGIN / (TILE_SIZE as f64 * cell_size_m)).log2()).round();
    z.clamp(0.0, 24.0) as u32
}

/// Writes bytes to a path, creating parent directories as needed.
fn write_bytes(path: &str, bytes: &[u8]) -> Result<(), ToolError> {
    if let Some(parent) = Path::new(path).parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|e| {
                ToolError::Execution(format!("failed creating tile directory: {e}"))
            })?;
        }
    }
    std::fs::write(path, bytes)
        .map_err(|e| ToolError::Execution(format!("failed writing tile: {e}")))
}
