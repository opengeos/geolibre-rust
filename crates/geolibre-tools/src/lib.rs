//! New geospatial tools that extend the `whitebox_next_gen` suite.
//!
//! Each tool implements [`wbcore::Tool`] (the same trait whitebox's own tools
//! use), so they plug into the registry alongside `register_default_tools` and
//! are exposed through the WASI runner exactly like the built-in tools.
//!
//! Add a new tool by creating a module with a `Tool` impl and pushing it in
//! [`geolibre_tools`].

mod aggregate_polygons;
mod apportion_polygon;
mod assign_projection;
mod build_balanced_zones;
mod cartogram;
mod central_feature;
mod collapse_dual_lines_to_centerline;
mod common;
mod corridor;
mod count_overlapping_features;
mod cut_fill;
mod delineate_built_up_areas;
mod delineate_depressions;
mod delineate_mounts;
mod dem_filter;
mod directional_distribution;
mod eliminate_polygons;
mod emerging_hot_spot_analysis;
mod expand_shrink;
mod extract_sinks;
mod fill;
mod generate_transects_along_lines;
mod geographically_weighted_regression;
mod geoparquet_io;
mod h3_polyfill;
mod interpolate_shape;
mod line_of_sight;
mod h3_to_vector;
mod hilbert;
mod incremental_spatial_autocorrelation;
mod lidar_common;
mod multiple_ring_buffer;
mod polygonize;
mod pmtiles;
mod pmtiles_extract;
mod polygon_neighbors;
mod raster_normalize;
mod raster_to_h3;
mod raster_to_tiles;
mod regions;
mod regularize_building_footprints;
mod render;
mod render_png;
mod render_vector_png;
mod reproject_raster;
mod ripleys_k;
mod simplify_shared_edges;
mod smooth_natural_features;
mod smooth_shared_edges;
mod spectral_index;
mod split_by_attributes;
mod subdivide_polygon;
mod tabulate_intersection;
mod thin_road_network;
mod vector_common;
mod vector_convert;
mod vector_to_h3;
mod vector_to_pmtiles;
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
        Box::new(assign_projection::AssignProjectionRasterTool),
        Box::new(assign_projection::AssignProjectionVectorTool),
        Box::new(assign_projection::AssignProjectionLidarTool),
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
        Box::new(regularize_building_footprints::RegularizeBuildingFootprintsTool),
        Box::new(smooth_natural_features::SmoothNaturalFeaturesTool),
        Box::new(eliminate_polygons::EliminatePolygonsTool),
        Box::new(simplify_shared_edges::SimplifySharedEdgesTool),
        Box::new(smooth_shared_edges::SmoothSharedEdgesTool),
        Box::new(emerging_hot_spot_analysis::EmergingHotSpotAnalysisTool),
        Box::new(line_of_sight::LineOfSightTool),
        Box::new(corridor::CorridorTool),
        Box::new(interpolate_shape::InterpolateShapeTool),
        Box::new(collapse_dual_lines_to_centerline::CollapseDualLinesToCenterlineTool),
        Box::new(count_overlapping_features::CountOverlappingFeaturesTool),
        Box::new(subdivide_polygon::SubdividePolygonTool),
        Box::new(generate_transects_along_lines::GenerateTransectsAlongLinesTool),
        Box::new(polygon_neighbors::PolygonNeighborsTool),
        Box::new(split_by_attributes::SplitByAttributesTool),
        Box::new(incremental_spatial_autocorrelation::IncrementalSpatialAutocorrelationTool),
        Box::new(apportion_polygon::ApportionPolygonTool),
        Box::new(central_feature::CentralFeatureTool),
        Box::new(expand_shrink::ExpandShrinkTool),
        Box::new(delineate_built_up_areas::DelineateBuiltUpAreasTool),
        Box::new(aggregate_polygons::AggregatePolygonsTool),
        Box::new(multiple_ring_buffer::MultipleRingBufferTool),
        Box::new(directional_distribution::DirectionalDistributionTool),
        Box::new(tabulate_intersection::TabulateIntersectionTool),
        Box::new(cut_fill::CutFillTool),
        Box::new(ripleys_k::RipleysKTool),
        Box::new(geographically_weighted_regression::GeographicallyWeightedRegressionTool),
        Box::new(build_balanced_zones::BuildBalancedZonesTool),
        Box::new(cartogram::CartogramTool),
        Box::new(thin_road_network::ThinRoadNetworkTool),
        Box::new(vector_to_h3::VectorToH3Tool),
        Box::new(h3_to_vector::H3ToVectorTool),
        Box::new(h3_polyfill::H3PolyfillTool),
        Box::new(raster_to_h3::RasterToH3Tool),
        Box::new(render_vector_png::RenderVectorPngTool),
        Box::new(write_pmtiles::WritePmTilesTool),
        Box::new(vector_to_pmtiles::VectorToPmTilesTool),
        Box::new(pmtiles_extract::PmtilesExtractTool),
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
    let lidar_in = ToolParamSchema::input_lidar;
    let lidar_out = || ToolParamSchema::output(ToolDatasetSchema::Lidar);
    let int = ToolParamSchema::scalar_integer;
    let float = ToolParamSchema::scalar_float;
    let colormaps = || {
        ToolParamSchema::enum_values(&["viridis", "magma", "turbo", "terrain", "grayscale"])
    };

    let map = match tool_id {
        "assign_projection_raster" => schemas(&[
            ("input", raster_in()),
            ("epsg", int()),
            ("output", raster_out()),
        ]),
        "assign_projection_vector" => schemas(&[
            ("input", vector_in()),
            ("epsg", int()),
            ("output", vector_out()),
        ]),
        "assign_projection_lidar" => schemas(&[
            ("input", lidar_in()),
            ("epsg", int()),
            ("output", lidar_out()),
        ]),
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
        "vector_to_pmtiles" => schemas(&[
            ("input", vector_in()),
            ("output", file_out()),
            ("min_zoom", int()),
            ("max_zoom", int()),
            ("layer_name", ToolParamSchema::string()),
            ("simplify", ToolParamSchema::bool()),
            ("drop_rate", float()),
        ]),
        "pmtiles_extract" => schemas(&[
            ("input", ToolParamSchema::input(ToolDatasetSchema::File)),
            ("output", file_out()),
            ("bbox", ToolParamSchema::string()),
            ("min_zoom", int()),
            ("max_zoom", int()),
            ("max_tiles", int()),
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
        "regularize_building_footprints" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("method", ToolParamSchema::enum_values(&[
                "right_angles",
                "right_angles_and_diagonals",
                "any_angle",
                "circle",
            ])),
            ("tolerance", float()),
            ("diagonal_penalty", float()),
            ("min_radius", float()),
            ("max_radius", float()),
        ]),
        "smooth_natural_features" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("segment_length", float()),
            ("iterations", int()),
            ("preserve_area", ToolParamSchema::bool()),
        ]),
        "eliminate_polygons" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("max_area", float()),
            ("where", ToolParamSchema::string()),
            ("exclude", ToolParamSchema::string()),
            ("strategy", ToolParamSchema::enum_values(&["longest_border", "largest_area"])),
            ("tolerance", float()),
        ]),
        "simplify_shared_edges" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("tolerance", float()),
            ("simplify_boundary", ToolParamSchema::bool()),
            ("snap_tolerance", float()),
        ]),
        "smooth_shared_edges" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("algorithm", ToolParamSchema::enum_values(&["paek", "bezier"])),
            ("tolerance", float()),
            ("smooth_boundary", ToolParamSchema::bool()),
            ("snap_tolerance", float()),
        ]),
        "emerging_hot_spot_analysis" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("time_field", ToolParamSchema::string()),
            ("time_step", ToolParamSchema::string()),
            ("value_field", ToolParamSchema::string()),
            ("resolution", int()),
            ("neighborhood", int()),
            ("time_window", int()),
        ]),
        "line_of_sight" => schemas(&[
            ("dem", raster_in()),
            ("observers", vector_in()),
            ("targets", vector_in()),
            ("output", vector_out()),
            ("observer_offset", float()),
            ("target_offset", float()),
            ("pair_field", ToolParamSchema::string()),
            ("band", int()),
        ]),
        "corridor" => schemas(&[
            ("cost1", raster_in()),
            ("cost2", raster_in()),
            ("cost", raster_in()),
            ("source1", raster_in()),
            ("source2", raster_in()),
            ("output", raster_out()),
            ("threshold", float()),
            ("percent", float()),
            ("band", int()),
        ]),
        "interpolate_shape" => schemas(&[
            ("input", vector_in()),
            ("surface", raster_in()),
            ("output", vector_out()),
            ("sample_distance", float()),
            ("method", ToolParamSchema::enum_values(&["bilinear", "nearest"])),
            ("attributes", ToolParamSchema::string()),
            ("band", int()),
        ]),
        "collapse_dual_lines_to_centerline" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("min_width", float()),
            ("max_width", float()),
            ("attribute", ToolParamSchema::string()),
            ("sample_distance", float()),
            ("min_overlap", float()),
        ]),
        "count_overlapping_features" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("min_count", int()),
            ("id_field", ToolParamSchema::string()),
            ("report_ids", table_out()),
        ]),
        "subdivide_polygon" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("method", ToolParamSchema::enum_values(&["equal_parts", "equal_areas"])),
            ("num_parts", int()),
            ("target_area", float()),
            ("angle", float()),
        ]),
        "generate_transects_along_lines" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("interval", float()),
            ("length", float()),
            ("offset", float()),
            ("include_ends", ToolParamSchema::bool()),
        ]),
        "polygon_neighbors" => schemas(&[
            ("input", vector_in()),
            ("output", table_out()),
            ("id_field", ToolParamSchema::string()),
            ("both_sides", ToolParamSchema::bool()),
            ("snap_tolerance", float()),
        ]),
        "split_by_attributes" => schemas(&[
            ("input", vector_in()),
            ("output_dir", file_out()),
            ("fields", ToolParamSchema::string()),
            ("format", ToolParamSchema::enum_values(&["geojson", "fgb", "parquet", "shp"])),
        ]),
        "incremental_spatial_autocorrelation" => schemas(&[
            ("input", vector_in()),
            ("field", ToolParamSchema::string()),
            ("output", table_out()),
            ("begin_distance", float()),
            ("increment", float()),
            ("num_bands", int()),
        ]),
        "apportion_polygon" => schemas(&[
            ("target", vector_in()),
            ("source", vector_in()),
            ("fields", ToolParamSchema::string()),
            ("output", vector_out()),
            ("method", ToolParamSchema::enum_values(&["area", "weight"])),
            ("weight_field", ToolParamSchema::string()),
            ("suffix", ToolParamSchema::string()),
        ]),
        "central_feature" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("statistic", ToolParamSchema::enum_values(&["central_feature", "linear_directional_mean"])),
            ("weight_field", ToolParamSchema::string()),
            ("case_field", ToolParamSchema::string()),
            ("distance", ToolParamSchema::enum_values(&["euclidean", "manhattan"])),
            ("orientation_only", ToolParamSchema::bool()),
        ]),
        "expand_shrink" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            ("classes", ToolParamSchema::string()),
            ("cells", int()),
            ("mode", ToolParamSchema::enum_values(&["expand", "shrink"])),
            ("band", int()),
        ]),
        "delineate_built_up_areas" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("grouping_distance", float()),
            ("min_building_count", int()),
            ("min_area", float()),
            ("simplify_tolerance", float()),
        ]),
        "aggregate_polygons" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("aggregation_distance", float()),
            ("min_area", float()),
            ("min_hole_size", float()),
            ("barrier", vector_in()),
        ]),
        "multiple_ring_buffer" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("distances", ToolParamSchema::string()),
            ("ring_type", ToolParamSchema::enum_values(&["rings", "disks"])),
            ("dissolve", ToolParamSchema::enum_values(&["none", "per_ring"])),
            ("distance_field", ToolParamSchema::string()),
        ]),
        "directional_distribution" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("statistic", ToolParamSchema::enum_values(&[
                "mean_center",
                "median_center",
                "central_feature",
                "standard_distance",
                "standard_deviational_ellipse",
            ])),
            ("weight_field", ToolParamSchema::string()),
            ("case_field", ToolParamSchema::string()),
            ("n_std", int()),
        ]),
        "tabulate_intersection" => schemas(&[
            ("input", vector_in()),
            ("class_features", vector_in()),
            ("output", vector_out()),
            ("class_field", ToolParamSchema::string()),
            ("sum_fields", ToolParamSchema::string()),
            ("zone_field", ToolParamSchema::string()),
        ]),
        "cut_fill" => schemas(&[
            ("input", raster_in()),
            ("after", raster_in()),
            ("plane", float()),
            ("output", raster_out()),
            ("band", int()),
            ("tolerance", float()),
            ("region_output", raster_out()),
            ("csv_output", table_out()),
        ]),
        "ripleys_k" => schemas(&[
            ("input", vector_in()),
            ("output", table_out()),
            ("distance_bands", int()),
            ("max_distance", float()),
            ("permutations", int()),
            ("weight_field", ToolParamSchema::string()),
            ("seed", int()),
        ]),
        "geographically_weighted_regression" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("y_field", ToolParamSchema::string()),
            ("x_fields", ToolParamSchema::string()),
            ("kernel", ToolParamSchema::enum_values(&["gaussian", "bisquare"])),
            ("bandwidth_type", ToolParamSchema::enum_values(&["adaptive", "fixed"])),
            ("bandwidth", float()),
        ]),
        "build_balanced_zones" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("zones", int()),
            ("criterion", ToolParamSchema::enum_values(&["homogeneity", "equal_count", "equal_sum"])),
            ("fields", ToolParamSchema::string()),
            ("contiguity", ToolParamSchema::enum_values(&["rook", "queen"])),
            ("tolerance", float()),
        ]),
        "cartogram" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("value_field", ToolParamSchema::string()),
            ("method", ToolParamSchema::enum_values(&["non_contiguous", "dorling"])),
            ("iterations", int()),
        ]),
        "thin_road_network" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("min_length", float()),
            ("hierarchy_field", ToolParamSchema::string()),
            ("visibility_field", ToolParamSchema::string()),
            ("keep_only", ToolParamSchema::bool()),
            ("snap_tolerance", float()),
        ]),
        "vector_to_h3" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("resolution", int()),
        ]),
        "h3_to_vector" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("field", ToolParamSchema::string()),
        ]),
        "h3_polyfill" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("resolution", int()),
        ]),
        "raster_to_h3" => schemas(&[
            ("input", raster_in()),
            ("output", vector_out()),
            ("resolution", int()),
            ("band", int()),
            ("aggregate", ToolParamSchema::enum_values(&["mean", "sum", "min", "max", "count", "median"])),
        ]),
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
}
