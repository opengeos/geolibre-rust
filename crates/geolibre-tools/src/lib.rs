//! New geospatial tools that extend the `whitebox_next_gen` suite.
//!
//! Each tool implements [`wbcore::Tool`] (the same trait whitebox's own tools
//! use), so they plug into the registry alongside `register_default_tools` and
//! are exposed through the WASI runner exactly like the built-in tools.
//!
//! Add a new tool by creating a module with a `Tool` impl and pushing it in
//! [`geolibre_tools`].

mod common;
mod delineate_depressions;
mod delineate_mounts;
mod dem_filter;
mod extract_sinks;
mod fill;
mod geoparquet_io;
mod hilbert;
mod polygonize;
mod pmtiles;
mod raster_normalize;
mod raster_to_tiles;
mod regions;
mod render;
mod render_png;
mod render_vector_png;
mod reproject_raster;
mod spectral_index;
mod vector_common;
mod vector_convert;
mod write_pmtiles;

use wbcore::Tool;

/// Returns every GeoLibre-authored tool as a boxed [`Tool`].
///
/// The binding layer (e.g. `geolibre-cli`) registers these into the same
/// registry as whitebox's built-in tools:
///
/// ```ignore
/// let mut registry = ToolRegistry::new();
/// register_default_tools(&mut registry);            // whitebox's ~733 tools
/// for tool in geolibre_tools::geolibre_tools() {     // plus GeoLibre's new ones
///     registry.register(tool);
/// }
/// ```
pub fn geolibre_tools() -> Vec<Box<dyn Tool>> {
    vec![
        Box::new(raster_normalize::RasterNormalizeTool),
        Box::new(dem_filter::DemFilterTool),
        Box::new(extract_sinks::ExtractSinksTool),
        Box::new(delineate_depressions::DelineateDepressionsTool),
        Box::new(delineate_mounts::DelineateMountsTool),
        Box::new(reproject_raster::ReprojectRasterTool),
        Box::new(render_png::RenderPngTool),
        Box::new(raster_to_tiles::RasterToTilesTool),
        Box::new(geoparquet_io::WriteGeoParquetTool),
        Box::new(geoparquet_io::ReadGeoParquetTool),
        Box::new(spectral_index::SpectralIndexTool),
        Box::new(vector_convert::VectorConvertTool),
        Box::new(render_vector_png::RenderVectorPngTool),
        Box::new(write_pmtiles::WritePmTilesTool),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exposes_at_least_one_tool() {
        let tools = geolibre_tools();
        assert!(!tools.is_empty());
        // Every tool must have a non-empty id.
        for tool in &tools {
            assert!(!tool.metadata().id.is_empty());
        }
    }
}
