//! PMTiles writing, re-exported from the shared `geolibre-pmtiles` core.
//!
//! The v3 format code (header, directories, varints) plus the bbox extract
//! engine live in `crates/geolibre-pmtiles` so the wasm-bindgen browser
//! library can use them without pulling this crate's heavy raster/vector
//! dependencies. This module keeps the original `pmtiles::build` surface for
//! the tiling tools.

pub use geolibre_pmtiles::writer::Tile;
pub use geolibre_pmtiles::LonLatBounds;

use wbcore::ToolError;

/// Builds a complete PMTiles v3 archive from PNG tiles. See
/// [`geolibre_pmtiles::writer::build_png`].
pub fn build(
    tiles: Vec<Tile>,
    bounds: &LonLatBounds,
    min_zoom: u8,
    max_zoom: u8,
) -> Result<Vec<u8>, ToolError> {
    geolibre_pmtiles::writer::build_png(tiles, bounds, min_zoom, max_zoom)
        .map_err(|e| ToolError::Execution(e.to_string()))
}
