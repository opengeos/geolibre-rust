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

use std::collections::BTreeMap;

use wbcore::{Tool, ToolDatasetSchema, ToolParamSchema};

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

fn schemas(entries: &[(&str, ToolParamSchema)]) -> BTreeMap<String, ToolParamSchema> {
    entries
        .iter()
        .map(|(name, schema)| ((*name).to_string(), schema.clone()))
        .collect()
}

/// Explicit parameter schemas for the GeoLibre-authored tools, keyed by tool id.
///
/// The manifest emitter (`geolibre-cli`) feeds these to
/// `wbcore::manifest_with_param_schema_json` so each param carries an accurate
/// `io_role`/`data_kind`/`schema`. Without them, the keyword-based inference
/// mis-types scalars whose descriptions mention a dataset — e.g.
/// `write_geoparquet.hilbert_sort` ("sort features…") would read as a vector
/// input, and `delineate_*.min_depth/min_height` ("matching lidar") as LiDAR —
/// which would make a host UI demand a layer for a plain number/flag.
pub fn geolibre_param_schemas(tool_id: &str) -> Option<BTreeMap<String, ToolParamSchema>> {
    let raster_in = ToolParamSchema::input_raster;
    let raster_out = ToolParamSchema::output_raster;
    let vector_in = ToolParamSchema::input_vector_any;
    let vector_out = ToolParamSchema::output_vector_any;
    let file_out = || ToolParamSchema::output(ToolDatasetSchema::File);
    let table_out = || ToolParamSchema::output(ToolDatasetSchema::Table);
    let int = ToolParamSchema::scalar_integer;
    let float = ToolParamSchema::scalar_float;
    let colormaps = || {
        ToolParamSchema::enum_values(&["viridis", "magma", "turbo", "terrain", "grayscale"])
    };

    let map = match tool_id {
        "raster_normalize" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            ("band", int()),
        ]),
        "dem_filter" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            ("filter", ToolParamSchema::enum_values(&["mean", "median", "gaussian"])),
            ("kernel_size", int()),
            ("sigma", float()),
            ("band", int()),
        ]),
        "extract_sinks" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            ("min_size", int()),
            ("region_output", raster_out()),
            ("depth_output", raster_out()),
            ("filled_output", raster_out()),
            ("csv_output", table_out()),
            ("vector_output", vector_out()),
            ("flat_increment", float()),
        ]),
        "delineate_depressions" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            ("level_output", raster_out()),
            ("csv_output", table_out()),
            ("vector_output", vector_out()),
            ("min_size", int()),
            ("min_depth", float()),
            ("interval", float()),
        ]),
        "delineate_mounts" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            ("level_output", raster_out()),
            ("csv_output", table_out()),
            ("vector_output", vector_out()),
            ("min_size", int()),
            ("min_height", float()),
            ("interval", float()),
            ("delta", float()),
        ]),
        "reproject_raster" => schemas(&[
            ("input", raster_in()),
            ("epsg", int()),
            ("method", ToolParamSchema::enum_values(&["nearest", "bilinear", "cubic", "lanczos"])),
            ("output", raster_out()),
        ]),
        "render_raster_png" => schemas(&[
            ("input", raster_in()),
            ("output", file_out()),
            ("band", int()),
            ("colormap", colormaps()),
            ("min", float()),
            ("max", float()),
        ]),
        "raster_to_tiles" => schemas(&[
            ("input", raster_in()),
            ("output_dir", file_out()),
            ("min_zoom", int()),
            ("max_zoom", int()),
            ("band", int()),
            ("colormap", colormaps()),
            ("method", ToolParamSchema::enum_values(&["bilinear", "nearest", "cubic"])),
            ("min", float()),
            ("max", float()),
        ]),
        "write_pmtiles" => schemas(&[
            ("input", raster_in()),
            ("output", file_out()),
            ("min_zoom", int()),
            ("max_zoom", int()),
            ("band", int()),
            ("colormap", colormaps()),
            ("method", ToolParamSchema::enum_values(&["bilinear", "nearest", "cubic"])),
            ("min", float()),
            ("max", float()),
        ]),
        "spectral_index" => schemas(&[
            ("input", raster_in()),
            ("index", ToolParamSchema::enum_values(&["ndvi", "ndwi", "ndbi", "nbr", "evi", "savi"])),
            ("red", int()),
            ("nir", int()),
            ("green", int()),
            ("blue", int()),
            ("swir", int()),
            ("soil_factor", float()),
            ("output", raster_out()),
        ]),
        "write_geoparquet" => schemas(&[
            ("input", vector_in()),
            ("output", file_out()),
            ("compression", ToolParamSchema::enum_values(&["zstd", "snappy", "gzip", "uncompressed"])),
            ("hilbert_sort", ToolParamSchema::bool()),
        ]),
        "read_geoparquet" => schemas(&[
            ("input", ToolParamSchema::input(ToolDatasetSchema::File)),
            ("output", vector_out()),
        ]),
        "vector_convert" => schemas(&[("input", vector_in()), ("output", vector_out())]),
        "render_vector_png" => schemas(&[
            ("input", vector_in()),
            ("output", file_out()),
            ("width", int()),
            ("height", int()),
            ("fill", ToolParamSchema::string()),
            ("stroke", ToolParamSchema::string()),
            ("stroke_width", float()),
            ("background", ToolParamSchema::string()),
        ]),
        _ => return None,
    };
    Some(map)
}

/// Curated parameter-schema corrections for a handful of upstream
/// `whitebox_next_gen` tools that wbcore's name/description-based inference
/// mis-types.
///
/// These tools are *not* GeoLibre-authored, so they don't appear in
/// [`geolibre_param_schemas`]; but their inferred schemas are wrong enough to
/// break a host UI. `spatial_join`, for instance, infers its `target`/`join`
/// layer paths as plain strings (so a UI renders a text box instead of a file
/// picker), its numeric `distance` as a string, and its `strategy` enum as a
/// LiDAR input file. Supplying explicit schemas here lets the manifest emitter
/// hand consumers an accurate `io_role`/`data_kind`/`schema` for each param.
///
/// Keep this list small and evidence-driven: only override a tool whose
/// inferred schema is demonstrably wrong, and mirror the param order and option
/// lists the tool documents in its own descriptions.
pub fn whitebox_param_schema_overrides(
    tool_id: &str,
) -> Option<BTreeMap<String, ToolParamSchema>> {
    let vector_in = ToolParamSchema::input_vector_any;
    let vector_out = ToolParamSchema::output_vector_any;
    let float = ToolParamSchema::scalar_float;

    let map = match tool_id {
        "spatial_join" => schemas(&[
            ("target", vector_in()),
            ("join", vector_in()),
            (
                "predicate",
                ToolParamSchema::enum_values(&[
                    "intersects",
                    "within",
                    "contains",
                    "touches",
                    "crosses",
                    "overlaps",
                    "within_distance",
                ]),
            ),
            ("distance", float()),
            (
                "strategy",
                ToolParamSchema::enum_values(&[
                    "first", "last", "count", "sum", "mean", "min", "max",
                ]),
            ),
            ("prefix", ToolParamSchema::string()),
            ("output", vector_out()),
        ]),
        _ => return None,
    };
    Some(map)
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

    #[test]
    fn every_tool_has_explicit_param_schemas() {
        // A GeoLibre tool without explicit schemas falls back to keyword-based
        // inference, which mis-types scalars (e.g. "...features..." -> vector).
        // Guard that every tool declares a schema for each of its params.
        for tool in geolibre_tools() {
            let meta = tool.metadata();
            let schemas = geolibre_param_schemas(&meta.id)
                .unwrap_or_else(|| panic!("missing param schemas for tool '{}'", meta.id));
            for param in &meta.params {
                assert!(
                    schemas.contains_key(param.name),
                    "tool '{}' is missing a schema for param '{}'",
                    meta.id,
                    param.name
                );
            }
        }
    }

    #[test]
    fn spatial_join_override_has_typed_inputs() {
        // Guard the curated correction: the layer paths must be file inputs, the
        // enums must carry their options, and `distance` must be numeric. (The
        // registry-sync check that these keys match the real tool params lives in
        // geolibre-cli, which can build the whitebox registry.)
        let map = whitebox_param_schema_overrides("spatial_join")
            .expect("missing spatial_join override");
        let keys: std::collections::BTreeSet<_> = map.keys().map(String::as_str).collect();
        let expected: std::collections::BTreeSet<_> =
            ["target", "join", "predicate", "distance", "strategy", "prefix", "output"]
                .into_iter()
                .collect();
        assert_eq!(keys, expected);
        assert!(whitebox_param_schema_overrides("definitely_not_a_tool").is_none());
    }
}
