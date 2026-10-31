//! GeoLibre tool: render a raster into a single PMTiles archive.
//!
//! Like `raster_to_tiles`, but instead of a `{z}/{x}/{y}.png` directory tree it
//! packs the Web Mercator tile pyramid into one `.pmtiles` file (the modern
//! single-file web-map tile format), so a whole basemap layer is one download.

use std::f64::consts::PI;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata,
    ToolParamSpec, ToolRunResult,
};
use wbraster::{Raster, ResampleMethod};

use crate::common::{load_input_raster, write_bytes};
use crate::pmtiles::{self, LonLatBounds, Tile};
use crate::raster_to_tiles::{
    native_zoom, render_tile, tile_range, MAX_TILES, ORIGIN, TILE_SIZE, WEB_MERCATOR_EPSG,
};
use crate::render::Colormap;
use crate::reproject_raster::parse_resample;

/// Web Mercator (EPSG:3857) sphere radius, for the meters -> lon/lat inverse.
const R: f64 = 6_378_137.0;

/// Renders a raster into a single PMTiles archive.
pub struct WritePmTilesTool;

impl Tool for WritePmTilesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "write_pmtiles",
            display_name: "Raster to PMTiles",
            summary: "Render a raster into a single PMTiles archive (Web Mercator PNG tile pyramid).",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec { name: "input", description: "Input raster file path. Must carry a source CRS (EPSG or WKT).", required: true },
                ToolParamSpec { name: "output", description: "Output PMTiles file path (e.g. /work/dem.pmtiles).", required: true },
                ToolParamSpec { name: "min_zoom", description: "Minimum zoom level (default: a single native zoom matching the raster resolution).", required: false },
                ToolParamSpec { name: "max_zoom", description: "Maximum zoom level (default: same as min_zoom).", required: false },
                ToolParamSpec { name: "band", description: "1-based band to render (default 1).", required: false },
                ToolParamSpec { name: "colormap", description: "Colormap: viridis (default), magma, turbo, terrain, or grayscale.", required: false },
                ToolParamSpec { name: "method", description: "Resampling method: bilinear (default), nearest, cubic.", required: false },
                ToolParamSpec { name: "min", description: "Value mapped to the low end of the colormap (default: band minimum).", required: false },
                ToolParamSpec { name: "max", description: "Value mapped to the high end of the colormap (default: band maximum).", required: false },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["input", "output"] {
            if args.get(key).and_then(Value::as_str).is_none() {
                return Err(ToolError::Validation(format!(
                    "missing required string parameter '{key}'"
                )));
            }
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
        let input = require_str(args, "input")?;
        let output = require_str(args, "output")?;
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

        ctx.progress.info("reprojecting to EPSG:3857");
        let merc: Raster = if raster.crs.epsg == Some(WEB_MERCATOR_EPSG) {
            raster
        } else {
            if raster.crs.epsg.is_none() && raster.crs.wkt.is_none() && raster.crs.proj4.is_none() {
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

        let native = native_zoom(merc.cell_size_x.abs());
        let min_zoom = args.get("min_zoom").and_then(Value::as_u64).map(|z| z as u32).unwrap_or(native);
        let max_zoom = args
            .get("max_zoom")
            .and_then(Value::as_u64)
            .map(|z| z as u32)
            .unwrap_or(min_zoom)
            .max(min_zoom);
        if max_zoom > 24 {
            return Err(ToolError::Validation("max_zoom must be <= 24".to_string()));
        }

        let ex = merc.extent();
        let mut tiles: Vec<Tile> = Vec::new();
        let mut total = 0usize;
        for z in min_zoom..=max_zoom {
            let n = 1u64 << z;
            let span = (2.0 * ORIGIN) / n as f64;
            let (tx_min, tx_max) = tile_range(ex.x_min, ex.x_max, -ORIGIN, span, n);
            let (ty_min, ty_max) = tile_range(ORIGIN - ex.y_max, ORIGIN - ex.y_min, 0.0, span, n);
            total += ((tx_max - tx_min + 1) as usize) * ((ty_max - ty_min + 1) as usize);
            if total > MAX_TILES {
                return Err(ToolError::Validation(format!(
                    "tile pyramid would exceed {MAX_TILES} tiles; narrow the zoom range (min_zoom/max_zoom)"
                )));
            }
            for tx in tx_min..=tx_max {
                for ty in ty_min..=ty_max {
                    let x0 = -ORIGIN + tx as f64 * span;
                    let y_top = ORIGIN - ty as f64 * span;
                    if let Some(png) =
                        render_tile(&merc, band, x0, y_top, span, method, colormap, min, max)?
                    {
                        tiles.push(Tile { z: z as u8, x: tx as u32, y: ty as u32, data: png });
                    }
                }
            }
            ctx.progress.progress((z - min_zoom + 1) as f64 / (max_zoom - min_zoom + 1) as f64);
        }

        let tile_count = tiles.len();
        if tile_count == 0 {
            return Err(ToolError::Execution(
                "no tiles produced (raster band is entirely no-data?)".to_string(),
            ));
        }

        let bounds = LonLatBounds {
            min_lon: merc_to_lon(ex.x_min),
            min_lat: merc_to_lat(ex.y_min),
            max_lon: merc_to_lon(ex.x_max),
            max_lat: merc_to_lat(ex.y_max),
        };

        ctx.progress.info("packing PMTiles archive");
        let archive = pmtiles::build(tiles, &bounds, min_zoom as u8, max_zoom as u8)?;
        write_bytes(output, &archive)?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(output));
        outputs.insert("min_zoom".to_string(), json!(min_zoom));
        outputs.insert("max_zoom".to_string(), json!(max_zoom));
        outputs.insert("tiles".to_string(), json!(tile_count));
        outputs.insert("bytes".to_string(), json!(archive.len()));
        let _ = TILE_SIZE; // referenced for clarity that tiles are 256px
        Ok(ToolRunResult { outputs })
    }
}

/// Web Mercator meters -> longitude (degrees).
fn merc_to_lon(x: f64) -> f64 {
    (x / ORIGIN) * 180.0
}

/// Web Mercator meters -> latitude (degrees).
fn merc_to_lat(y: f64) -> f64 {
    (2.0 * (y / R).exp().atan() - PI / 2.0).to_degrees()
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}
