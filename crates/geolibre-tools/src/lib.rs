//! New geospatial tools that extend the `whitebox_next_gen` suite.
//!
//! Each tool implements [`wbcore::Tool`] (the same trait whitebox's own tools
//! use), so they plug into the registry alongside `register_default_tools` and
//! are exposed through the WASI runner exactly like the built-in tools.
//!
//! Add a new tool by creating a module with a `Tool` impl and pushing it in
//! [`geolibre_tools`].

mod aggregate_points;
mod aggregate_polygons;
mod apportion_polygon;
mod assign_projection;
mod boundary_clean;
mod build_balanced_zones;
mod calculate_motion_statistics;
mod cartogram;
mod central_feature;
mod collapse_dual_lines_to_centerline;
mod colocation_analysis;
mod common;
mod corridor;
mod count_overlapping_features;
mod create_spatially_balanced_points;
mod cut_fill;
mod delineate_built_up_areas;
mod delineate_depressions;
mod delineate_mounts;
mod dem_filter;
mod detect_feature_changes;
mod detect_image_anomalies;
mod directional_distribution;
mod eliminate_polygons;
mod emerging_hot_spot_analysis;
mod expand_shrink;
mod extract_sinks;
mod fill;
mod find_identical;
mod find_space_time_matches;
mod fuzzy_overlay;
mod generate_od_links;
mod generate_transects_along_lines;
mod geographically_weighted_regression;
mod geoparquet_io;
mod h3_polyfill;
mod hdbscan;
mod interpolate_shape;
mod line_of_sight;
mod h3_to_vector;
mod hilbert;
mod incremental_spatial_autocorrelation;
mod integrate;
mod lidar_common;
mod multiple_ring_buffer;
mod neighborhood_summary_statistics;
mod path_distance;
mod polygonize;
mod pmtiles;
mod pmtiles_extract;
mod polygon_neighbors;
mod raster_normalize;
mod raster_to_h3;
mod raster_to_tiles;
mod reconstruct_tracks;
mod regions;
mod regularize_building_footprints;
mod remove_overlap_multiple;
mod render;
mod render_png;
mod render_vector_png;
mod reproject_raster;
mod resolve_building_conflicts;
mod ripleys_k;
mod rubbersheet_features;
mod similarity_search;
mod simplify_shared_edges;
mod smooth_natural_features;
mod smooth_shared_edges;
mod snap_tracks;
mod solar_radiation;
mod spectral_index;
mod split_by_attributes;
mod storage_capacity;
mod subdivide_polygon;
mod tabulate_intersection;
mod thin_road_network;
mod time_series_clustering;
mod trace_proximity_events;
mod vector_common;
mod vector_convert;
mod vector_to_h3;
mod vector_to_pmtiles;
mod write_pmtiles;

mod sort_features;

mod calculate_composite_index;

mod calculate_rates;

mod color_polygons;

mod dice;

mod spatial_outlier_detection;

mod bivariate_spatial_association;

mod generate_trend_raster;

mod warp_raster;

mod weighted_voronoi;

mod pycnophylactic_interpolation;

mod cost_connectivity;

mod locate_regions;

mod edgematch_features;

mod landtrendr;

mod local_outlier_analysis;

mod collapse_hydro_polygon;

mod change_point_detection;

mod time_series_forecast;
mod resolve_road_conflicts;
mod presence_only_prediction;
mod topo_to_raster;
mod collapse_road_detail;
mod analyze_changes_ccdc;
mod space_time_kernel_density;
mod geotagged_photos_to_points;
mod darcy_flow;
mod time_series_cross_correlation;
mod generalized_linear_regression;
mod interpolate_with_barriers;
mod convert_coordinate_notation;
mod repair_geometry;
mod grid_index_features;
mod local_bivariate_relationships;
mod dimension_reduction;
mod feature_outline_masks;
mod line_density;
mod pairwise_comparison_weights;
mod kernel_density_ratio;
mod detect_incidents;
mod find_argument_statistics;
mod strip_map_index_features;
mod zonal_histogram;

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
        Box::new(reconstruct_tracks::ReconstructTracksTool),
        Box::new(solar_radiation::SolarRadiationTool),
        Box::new(hdbscan::HdbscanTool),
        Box::new(colocation_analysis::ColocationAnalysisTool),
        Box::new(similarity_search::SimilaritySearchTool),
        Box::new(detect_feature_changes::DetectFeatureChangesTool),
        Box::new(integrate::IntegrateTool),
        Box::new(rubbersheet_features::RubbersheetFeaturesTool),
        Box::new(snap_tracks::SnapTracksTool),
        Box::new(remove_overlap_multiple::RemoveOverlapMultipleTool),
        Box::new(fuzzy_overlay::FuzzyOverlayTool),
        Box::new(aggregate_points::AggregatePointsTool),
        Box::new(generate_od_links::GenerateOdLinksTool),
        Box::new(neighborhood_summary_statistics::NeighborhoodSummaryStatisticsTool),
        Box::new(storage_capacity::StorageCapacityTool),
        Box::new(find_space_time_matches::FindSpaceTimeMatchesTool),
        Box::new(create_spatially_balanced_points::CreateSpatiallyBalancedPointsTool),
        Box::new(find_identical::FindIdenticalTool),
        Box::new(path_distance::PathDistanceTool),
        Box::new(time_series_clustering::TimeSeriesClusteringTool),
        Box::new(trace_proximity_events::TraceProximityEventsTool),
        Box::new(detect_image_anomalies::DetectImageAnomaliesTool),
        Box::new(resolve_building_conflicts::ResolveBuildingConflictsTool),
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
        Box::new(boundary_clean::BoundaryCleanTool),
        Box::new(calculate_motion_statistics::CalculateMotionStatisticsTool),
        Box::new(sort_features::SortFeaturesTool),
        Box::new(calculate_composite_index::CalculateCompositeIndexTool),
        Box::new(calculate_rates::CalculateRatesTool),
        Box::new(color_polygons::ColorPolygonsTool),
        Box::new(dice::DiceTool),
        Box::new(spatial_outlier_detection::SpatialOutlierDetectionTool),
        Box::new(bivariate_spatial_association::BivariateSpatialAssociationTool),
        Box::new(generate_trend_raster::GenerateTrendRasterTool),
        Box::new(warp_raster::WarpRasterTool),
        Box::new(weighted_voronoi::WeightedVoronoiTool),
        Box::new(pycnophylactic_interpolation::PycnophylacticInterpolationTool),
        Box::new(cost_connectivity::CostConnectivityTool),
        Box::new(locate_regions::LocateRegionsTool),
        Box::new(edgematch_features::EdgematchFeaturesTool),
        Box::new(landtrendr::LandtrendrTool),
        Box::new(local_outlier_analysis::LocalOutlierAnalysisTool),
        Box::new(collapse_hydro_polygon::CollapseHydroPolygonTool),
        Box::new(change_point_detection::ChangePointDetectionTool),
        Box::new(time_series_forecast::TimeSeriesForecastTool),
        Box::new(resolve_road_conflicts::ResolveRoadConflictsTool),
        Box::new(presence_only_prediction::PresenceOnlyPredictionTool),
        Box::new(topo_to_raster::TopoToRasterTool),
        Box::new(collapse_road_detail::CollapseRoadDetailTool),
        Box::new(analyze_changes_ccdc::AnalyzeChangesCcdcTool),
        Box::new(space_time_kernel_density::SpaceTimeKernelDensityTool),
        Box::new(geotagged_photos_to_points::GeotaggedPhotosToPointsTool),
        Box::new(darcy_flow::DarcyFlowTool),
        Box::new(time_series_cross_correlation::TimeSeriesCrossCorrelationTool),
        Box::new(generalized_linear_regression::GeneralizedLinearRegressionTool),
        Box::new(interpolate_with_barriers::InterpolateWithBarriersTool),
        Box::new(convert_coordinate_notation::ConvertCoordinateNotationTool),
        Box::new(repair_geometry::RepairGeometryTool),
        Box::new(grid_index_features::GridIndexFeaturesTool),
        Box::new(local_bivariate_relationships::LocalBivariateRelationshipsTool),
        Box::new(dimension_reduction::DimensionReductionTool),
        Box::new(feature_outline_masks::FeatureOutlineMasksTool),
        Box::new(line_density::LineDensityTool),
        Box::new(pairwise_comparison_weights::PairwiseComparisonWeightsTool),
        Box::new(kernel_density_ratio::KernelDensityRatioTool),
        Box::new(detect_incidents::DetectIncidentsTool),
        Box::new(find_argument_statistics::FindArgumentStatisticsTool),
        Box::new(strip_map_index_features::StripMapIndexFeaturesTool),
        Box::new(zonal_histogram::ZonalHistogramTool),
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
        "boundary_clean" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            ("method", ToolParamSchema::enum_values(&["majority", "expand_shrink"])),
            ("neighbors", int()),
            ("threshold", ToolParamSchema::enum_values(&["majority", "half"])),
            ("iterations", int()),
            ("sort", ToolParamSchema::enum_values(&["descending", "ascending", "none"])),
            ("band", int()),
        ]),
        "expand_shrink" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            ("classes", ToolParamSchema::string()),
            ("cells", int()),
            ("mode", ToolParamSchema::enum_values(&["expand", "shrink"])),
            ("band", int()),
        ]),
        "resolve_building_conflicts" => schemas(&[
            ("buildings", vector_in()),
            ("barriers", vector_in()),
            ("output", vector_out()),
            ("barrier_width", float()),
            ("gap", float()),
            ("min_size", float()),
            ("hide", ToolParamSchema::bool()),
        ]),
        "detect_image_anomalies" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            ("mode", ToolParamSchema::enum_values(&["global", "local"])),
            ("window", int()),
            ("threshold", float()),
            ("mask_output", raster_out()),
        ]),
        "trace_proximity_events" => schemas(&[
            ("input", vector_in()),
            ("track_field", ToolParamSchema::string()),
            ("time_field", ToolParamSchema::string()),
            ("output", vector_out()),
            ("search_distance", float()),
            ("min_duration", ToolParamSchema::string()),
            ("entities", ToolParamSchema::string()),
            ("depth", int()),
        ]),
        "time_series_clustering" => schemas(&[
            ("input", vector_in()),
            ("time_field", ToolParamSchema::string()),
            ("output", vector_out()),
            ("num_clusters", int()),
            (
                "characteristic",
                ToolParamSchema::enum_values(&["value", "profile", "correlation"]),
            ),
            ("time_step", ToolParamSchema::string()),
            ("value_field", ToolParamSchema::string()),
            ("resolution", int()),
            ("seed", int()),
        ]),
        "path_distance" => schemas(&[
            ("source", raster_in()),
            ("output", raster_out()),
            ("cost", raster_in()),
            ("surface", raster_in()),
            (
                "vertical_factor",
                ToolParamSchema::enum_values(&[
                    "tobler",
                    "linear",
                    "sym_linear",
                    "inverse_linear",
                    "binary",
                ]),
            ),
            ("slope_factor", float()),
            ("zero_factor", float()),
            ("max_slope", float()),
            ("band", int()),
        ]),
        "find_identical" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("fields", ToolParamSchema::string()),
            ("compare_geometry", ToolParamSchema::bool()),
            ("xy_tolerance", float()),
            ("mode", ToolParamSchema::enum_values(&["report", "delete"])),
        ]),
        "create_spatially_balanced_points" => schemas(&[
            ("constraint", vector_in()),
            ("output", vector_out()),
            ("num_points", int()),
            ("probability", raster_in()),
            ("seed", int()),
        ]),
        "find_space_time_matches" => schemas(&[
            ("primary", vector_in()),
            ("secondary", vector_in()),
            ("time_field", ToolParamSchema::string()),
            ("output", vector_out()),
            ("search_distance", float()),
            ("time_window", ToolParamSchema::string()),
            (
                "temporal_relationship",
                ToolParamSchema::enum_values(&["either", "before", "after"]),
            ),
            ("primary_id_field", ToolParamSchema::string()),
            ("secondary_id_field", ToolParamSchema::string()),
        ]),
        "storage_capacity" => schemas(&[
            ("dem", raster_in()),
            ("output", file_out()),
            ("zones", vector_in()),
            ("zone_id_field", ToolParamSchema::string()),
            ("num_levels", int()),
            ("increment", float()),
            ("min_elevation", float()),
            ("max_elevation", float()),
            ("band", int()),
        ]),
        "neighborhood_summary_statistics" => schemas(&[
            ("input", vector_in()),
            ("fields", ToolParamSchema::string()),
            ("output", vector_out()),
            (
                "neighborhood",
                ToolParamSchema::enum_values(&["knn", "distance_band", "contiguity"]),
            ),
            ("neighbors", int()),
            ("distance", float()),
            ("weights", ToolParamSchema::enum_values(&["uniform", "inverse_distance"])),
        ]),
        "generate_od_links" => schemas(&[
            ("origins", vector_in()),
            ("destinations", vector_in()),
            ("output", vector_out()),
            ("num_nearest", int()),
            ("search_distance", float()),
            ("id_field", ToolParamSchema::string()),
            ("origin_id_field", ToolParamSchema::string()),
            ("dest_id_field", ToolParamSchema::string()),
        ]),
        "aggregate_points" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("aggregation_distance", float()),
            ("min_points", int()),
            ("method", ToolParamSchema::enum_values(&["convex_hull", "buffer"])),
            ("sum_fields", ToolParamSchema::string()),
        ]),
        "fuzzy_overlay" => schemas(&[
            ("input", raster_in()),
            ("inputs", ToolParamSchema::string()),
            ("output", raster_out()),
            (
                "function",
                ToolParamSchema::enum_values(&[
                    "linear", "gaussian", "small", "large", "ms_small", "ms_large",
                ]),
            ),
            (
                "overlay",
                ToolParamSchema::enum_values(&["and", "or", "product", "sum", "gamma"]),
            ),
            ("midpoint", float()),
            ("spread", float()),
            ("min", float()),
            ("max", float()),
            ("gamma", float()),
            ("band", int()),
        ]),
        "calculate_motion_statistics" => schemas(&[
            ("input", vector_in()),
            ("track_field", ToolParamSchema::string()),
            ("time_field", ToolParamSchema::string()),
            ("output", vector_out()),
            ("window", int()),
            ("idle_distance", float()),
            ("idle_duration", float()),
        ]),
        "sort_features" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("method", ToolParamSchema::enum_values(&["hilbert", "attribute"])),
            ("fields", ToolParamSchema::string()),
            ("index_field", ToolParamSchema::string()),
        ]),
        "calculate_composite_index" => schemas(&[
            ("input", vector_in()),
            ("fields", ToolParamSchema::string()),
            ("output", vector_out()),
            ("scaling", ToolParamSchema::enum_values(&["minmax", "zscore", "percentile", "none"])),
            ("weights", ToolParamSchema::string()),
            ("combine", ToolParamSchema::enum_values(&["mean", "sum", "geometric_mean"])),
            ("output_range", ToolParamSchema::enum_values(&["minmax", "zero_to_100", "zscore", "none"])),
        ]),
        "calculate_rates" => schemas(&[
            ("input", vector_in()),
            ("count_field", ToolParamSchema::string()),
            ("population_field", ToolParamSchema::string()),
            ("output", vector_out()),
            ("method", ToolParamSchema::enum_values(&["crude", "eb_global", "eb_spatial"])),
            ("per", float()),
            ("neighbors", int()),
        ]),
        "color_polygons" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("field", ToolParamSchema::string()),
            ("adjacency", ToolParamSchema::enum_values(&["edge", "edge_or_corner"])),
            ("snap_tolerance", float()),
        ]),
        "dice" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("vertex_limit", int()),
        ]),
        "spatial_outlier_detection" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("neighbors", int()),
            ("percent_outlier", float()),
            ("threshold", float()),
        ]),
        "bivariate_spatial_association" => schemas(&[
            ("input", vector_in()),
            ("x_field", ToolParamSchema::string()),
            ("y_field", ToolParamSchema::string()),
            ("output", vector_out()),
            ("neighbors", int()),
            ("permutations", int()),
            ("seed", int()),
        ]),
        "generate_trend_raster" => schemas(&[
            ("inputs", ToolParamSchema::string()),
            ("output", raster_out()),
            ("times", ToolParamSchema::string()),
            ("method", ToolParamSchema::enum_values(&["linear", "mann_kendall"])),
            ("intercept_output", raster_out()),
            ("significance_output", raster_out()),
            ("min_valid", int()),
            ("band", int()),
        ]),
        "warp_raster" => schemas(&[
            ("input", raster_in()),
            ("gcps", ToolParamSchema::string()),
            ("output", raster_out()),
            ("transform", ToolParamSchema::enum_values(&["poly1", "poly2", "poly3"])),
            ("resampling", ToolParamSchema::enum_values(&["nearest", "bilinear"])),
            ("cell_size", float()),
            ("epsg", int()),
            ("band", int()),
        ]),
        "weighted_voronoi" => schemas(&[
            ("input", vector_in()),
            ("output", raster_out()),
            ("weight_field", ToolParamSchema::string()),
            ("weight_type", ToolParamSchema::enum_values(&["multiplicative", "additive", "power"])),
            ("cell_size", float()),
            ("margin", float()),
            ("epsg", int()),
        ]),
        "pycnophylactic_interpolation" => schemas(&[
            ("input", vector_in()),
            ("count_field", ToolParamSchema::string()),
            ("output", raster_out()),
            ("cell_size", float()),
            ("iterations", int()),
            ("tolerance", float()),
            ("non_negative", ToolParamSchema::bool()),
        ]),
        "cost_connectivity" => schemas(&[
            ("sources", vector_in()),
            ("cost", raster_in()),
            ("output", vector_out()),
            ("connections", ToolParamSchema::enum_values(&["mst", "all_neighbors"])),
            ("id_field", ToolParamSchema::string()),
            ("band", int()),
        ]),
        "locate_regions" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            ("total_area", float()),
            ("num_regions", int()),
            ("shape", float()),
            ("min_distance", float()),
            ("band", int()),
        ]),
        "edgematch_features" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("tolerance", float()),
            ("method", ToolParamSchema::enum_values(&["midpoint", "move_endpoint"])),
            ("match_fields", ToolParamSchema::string()),
            ("links", vector_out()),
        ]),
        "landtrendr" => schemas(&[
            ("inputs", ToolParamSchema::string()),
            ("output", raster_out()),
            ("years", ToolParamSchema::string()),
            ("magnitude_output", raster_out()),
            ("duration_output", raster_out()),
            ("direction", ToolParamSchema::enum_values(&["loss", "gain"])),
            ("max_segments", int()),
            ("spike_threshold", float()),
            ("min_valid", int()),
            ("band", int()),
        ]),
        "local_outlier_analysis" => schemas(&[
            ("input", vector_in()),
            ("time_field", ToolParamSchema::string()),
            ("output", vector_out()),
            ("value_field", ToolParamSchema::string()),
            ("time_step", float()),
            ("resolution", int()),
            ("kring", int()),
            ("time_window", int()),
            ("permutations", int()),
            ("seed", int()),
        ]),
        "collapse_hydro_polygon" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("collapse_width", float()),
            ("sample_distance", float()),
            ("min_length", float()),
            ("retained", vector_out()),
        ]),
        "change_point_detection" => schemas(&[
            ("input", vector_in()),
            ("time_field", ToolParamSchema::string()),
            ("output", vector_out()),
            ("value_field", ToolParamSchema::string()),
            ("change_type", ToolParamSchema::enum_values(&["mean", "slope"])),
            ("method", ToolParamSchema::enum_values(&["auto", "defined"])),
            ("num_change_points", int()),
            ("sensitivity", float()),
            ("time_step", float()),
            ("resolution", int()),
        ]),
        "time_series_forecast" => schemas(&[
            ("input", vector_in()),
            ("time_field", ToolParamSchema::string()),
            ("output", vector_out()),
            ("value_field", ToolParamSchema::string()),
            ("steps", int()),
            ("model", ToolParamSchema::enum_values(&["auto", "exp_smoothing", "linear", "parabolic"])),
            ("holdout", int()),
            ("time_step", float()),
            ("resolution", int()),
        ]),
        "reconstruct_tracks" => schemas(&[
            ("input", vector_in()),
            ("track_field", ToolParamSchema::string()),
            ("time_field", ToolParamSchema::string()),
            ("output", vector_out()),
            ("time_gap", float()),
            ("distance_gap", float()),
            ("dwells", vector_out()),
            ("dwell_distance", float()),
            ("dwell_min_duration", float()),
        ]),
        "hdbscan" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("min_cluster_size", int()),
            ("min_samples", int()),
        ]),
        "colocation_analysis" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("category_field", ToolParamSchema::string()),
            ("category_a", ToolParamSchema::string()),
            ("category_b", ToolParamSchema::string()),
            ("neighbors", int()),
            ("weight", ToolParamSchema::enum_values(&["gaussian", "uniform"])),
            ("permutations", int()),
            ("seed", int()),
        ]),
        "similarity_search" => schemas(&[
            ("reference", vector_in()),
            ("candidates", vector_in()),
            ("fields", ToolParamSchema::string()),
            ("output", vector_out()),
            ("match_method", ToolParamSchema::enum_values(&["euclidean", "cosine"])),
            ("most_or_least", ToolParamSchema::enum_values(&["most", "least", "both"])),
            ("num_results", int()),
        ]),
        "detect_feature_changes" => schemas(&[
            ("update", vector_in()),
            ("base", vector_in()),
            ("output", vector_out()),
            ("search_distance", float()),
            ("spatial_tolerance", float()),
            ("compare_fields", ToolParamSchema::string()),
        ]),
        "integrate" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("tolerance", float()),
            ("snap_to_edges", ToolParamSchema::bool()),
        ]),
        "rubbersheet_features" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("links", vector_in()),
            ("target", vector_in()),
            ("search_distance", float()),
            ("method", ToolParamSchema::enum_values(&["linear", "idw"])),
            ("power", float()),
        ]),
        "remove_overlap_multiple" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("method", ToolParamSchema::enum_values(&["center_line", "thiessen"])),
            ("grid_resolution", int()),
        ]),
        "snap_tracks" => schemas(&[
            ("input", vector_in()),
            ("network", vector_in()),
            ("track_field", ToolParamSchema::string()),
            ("time_field", ToolParamSchema::string()),
            ("output", vector_out()),
            ("search_distance", float()),
            ("max_candidates", int()),
        ]),
        "solar_radiation" => schemas(&[
            ("dem", raster_in()),
            ("output", raster_out()),
            ("direct_output", raster_out()),
            ("diffuse_output", raster_out()),
            ("start_day", ToolParamSchema::string()),
            ("end_day", ToolParamSchema::string()),
            ("day_interval", int()),
            ("time_step", float()),
            ("latitude", float()),
            ("transmittivity", float()),
            ("diffuse_proportion", float()),
            ("horizon_distance", int()),
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
        "find_argument_statistics" => schemas(&[
            ("inputs", ToolParamSchema::string()),
            ("output", raster_out()),
            ("statistic", ToolParamSchema::enum_values(&["argmax", "argmin", "median_position", "duration", "longest_run"])),
            ("threshold", float()),
            ("comparison", ToolParamSchema::enum_values(&[">", ">=", "<", "<="])),
            ("dates", ToolParamSchema::string()),
            ("min_valid", int()),
        ]),
        "detect_incidents" => schemas(&[
            ("input", vector_in()),
            ("track_field", ToolParamSchema::string()),
            ("time_field", ToolParamSchema::string()),
            ("start_condition", ToolParamSchema::string()),
            ("end_condition", ToolParamSchema::string()),
            ("mode", ToolParamSchema::enum_values(&["points", "segments"])),
            ("output", vector_out()),
        ]),
        "kernel_density_ratio" => schemas(&[
            ("input", vector_in()),
            ("denominator", vector_in()),
            ("output", raster_out()),
            ("weight_field", ToolParamSchema::string()),
            ("denominator_weight_field", ToolParamSchema::string()),
            ("bandwidth", float()),
            ("cell_size", float()),
            ("log_ratio", ToolParamSchema::bool()),
            ("denominator_floor", float()),
            ("epsg", int()),
        ]),
        "pairwise_comparison_weights" => schemas(&[
            ("matrix", ToolParamSchema::string()),
            ("input", ToolParamSchema::input(ToolDatasetSchema::Table)),
            ("criteria", ToolParamSchema::string()),
            ("output", table_out()),
        ]),
        "line_density" => schemas(&[
            ("input", vector_in()),
            ("output", raster_out()),
            ("weight_field", ToolParamSchema::string()),
            ("search_radius", float()),
            ("cell_size", float()),
            ("area_units", ToolParamSchema::string()),
        ]),
        "feature_outline_masks" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("margin", float()),
            ("mask_kind", ToolParamSchema::enum_values(&["exact", "convex_hull", "box"])),
            ("masked_layer", vector_in()),
            ("id_field", ToolParamSchema::string()),
        ]),
        "dimension_reduction" => schemas(&[
            ("input", vector_in()),
            ("fields", ToolParamSchema::string()),
            ("output", vector_out()),
            ("table", table_out()),
            ("num_components", int()),
            ("min_variance", float()),
            ("standardize", ToolParamSchema::bool()),
        ]),
        "local_bivariate_relationships" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("field1", ToolParamSchema::string()),
            ("field2", ToolParamSchema::string()),
            ("neighbors", int()),
            ("permutations", int()),
            ("significance", float()),
            ("seed", int()),
        ]),
        "strip_map_index_features" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("page_length", float()),
            ("page_width", float()),
            ("overlap", float()),
            ("orientation", ToolParamSchema::enum_values(&["along_line", "horizontal", "vertical"])),
            ("start_page", int()),
        ]),
        "grid_index_features" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("mode", ToolParamSchema::enum_values(&["grid", "strip"])),
            ("x_min", float()),
            ("y_min", float()),
            ("x_max", float()),
            ("y_max", float()),
            ("tile_width", float()),
            ("tile_height", float()),
            ("page_size", ToolParamSchema::enum_values(&["a0", "a1", "a2", "a3", "a4", "letter", "legal", "tabloid"])),
            ("map_scale", float()),
            ("origin_x", float()),
            ("origin_y", float()),
            ("naming", ToolParamSchema::enum_values(&["alphanumeric", "sequential"])),
            ("intersect_only", ToolParamSchema::bool()),
            ("route", vector_in()),
            ("overlap", float()),
            ("epsg", int()),
        ]),
        "repair_geometry" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("check_only", ToolParamSchema::bool()),
        ]),
        "convert_coordinate_notation" => schemas(&[
            ("input", vector_in()),
            ("input_notation", ToolParamSchema::enum_values(&["DD", "DMS", "DDM", "UTM", "MGRS"])),
            ("output_notation", ToolParamSchema::enum_values(&["DD", "DMS", "DDM", "UTM", "MGRS"])),
            ("coord_field", ToolParamSchema::string()),
            ("output_field", ToolParamSchema::string()),
            ("precision", int()),
            ("update_geometry", ToolParamSchema::bool()),
            ("output", vector_out()),
        ]),
        "interpolate_with_barriers" => schemas(&[
            ("input", vector_in()),
            ("field", ToolParamSchema::string()),
            ("output", raster_out()),
            ("barriers", vector_in()),
            ("method", ToolParamSchema::enum_values(&["idw", "local_polynomial"])),
            ("power", float()),
            ("bandwidth", float()),
            ("radius", float()),
            ("cell_size", float()),
        ]),
        "generalized_linear_regression" => schemas(&[
            ("input", vector_in()),
            ("dependent_field", ToolParamSchema::string()),
            ("explanatory_fields", ToolParamSchema::string()),
            ("family", ToolParamSchema::enum_values(&["gaussian", "poisson", "logistic"])),
            ("output", vector_out()),
            ("report", table_out()),
        ]),
        "time_series_cross_correlation" => schemas(&[
            ("input", raster_in()),
            ("secondary", raster_in()),
            ("output", raster_out()),
            ("corr_output", raster_out()),
            ("corr0_output", raster_out()),
            ("pvalue_output", raster_out()),
            ("min_lag", int()),
            ("max_lag", int()),
            ("detrend", ToolParamSchema::bool()),
            ("deseasonalize", ToolParamSchema::bool()),
            ("season_length", int()),
            ("min_valid", int()),
            ("band", int()),
        ]),
        "darcy_flow" => schemas(&[
            ("input", raster_in()),
            ("transmissivity", raster_in()),
            ("porosity", raster_in()),
            ("output", raster_out()),
            ("direction", raster_out()),
            ("band", int()),
            ("seeds", vector_in()),
            ("streamlines", vector_out()),
            ("step", float()),
            ("max_steps", int()),
        ]),
        "geotagged_photos_to_points" => schemas(&[
            ("input", ToolParamSchema::string()),
            ("output", vector_out()),
            ("recursive", ToolParamSchema::bool()),
            ("only_geotagged", ToolParamSchema::bool()),
        ]),
        "space_time_kernel_density" => schemas(&[
            ("input", vector_in()),
            ("output", raster_out()),
            ("time_field", ToolParamSchema::string()),
            ("time_step", ToolParamSchema::string()),
            ("temporal_bandwidth", ToolParamSchema::string()),
            ("spatial_bandwidth", float()),
            ("cell_size", float()),
            ("weight_field", ToolParamSchema::string()),
            ("spatial_kernel", ToolParamSchema::enum_values(&["epanechnikov", "quartic"])),
            ("temporal_kernel", ToolParamSchema::enum_values(&["triangular", "epanechnikov"])),
            ("epsg", int()),
        ]),
        "analyze_changes_ccdc" => schemas(&[
            ("input", ToolParamSchema::string()),
            ("output", raster_out()),
            ("dates", ToolParamSchema::string()),
            ("period", float()),
            ("harmonic_order", int()),
            ("change_threshold", float()),
            ("min_consecutive", int()),
            ("min_observations", int()),
            ("band", int()),
            ("break_date_output", raster_out()),
            ("rmse_output", raster_out()),
            ("slope_output", raster_out()),
            ("amplitude_output", raster_out()),
        ]),
        "collapse_road_detail" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("collapse_distance", float()),
            ("road_class_field", ToolParamSchema::string()),
            ("snap_tolerance", float()),
        ]),
        "topo_to_raster" => schemas(&[
            ("contours", vector_in()),
            ("points", vector_in()),
            ("streams", vector_in()),
            ("output", raster_out()),
            ("elevation_field", ToolParamSchema::string()),
            ("cell_size", float()),
            ("tension", float()),
            ("iterations", int()),
            ("tolerance", float()),
            ("enforce_drainage", ToolParamSchema::bool()),
            ("stream_burn", float()),
        ]),
        "presence_only_prediction" => schemas(&[
            ("input", vector_in()),
            ("explanatory", raster_in()),
            ("output", raster_out()),
            ("report", file_out()),
            ("features", ToolParamSchema::string()),
            ("background", int()),
            ("regularization", float()),
            ("hinge_knots", int()),
            ("seed", int()),
        ]),
        "resolve_road_conflicts" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("symbol_width", float()),
            ("symbol_width_field", ToolParamSchema::string()),
            ("hierarchy_field", ToolParamSchema::string()),
            ("scale", float()),
            ("gap", float()),
            ("max_iter", int()),
            ("pin_endpoints", ToolParamSchema::bool()),
            ("links", vector_out()),
        ]),
        "zonal_histogram" => schemas(&[
            ("zones", raster_in()),
            ("value", raster_in()),
            ("output", table_out()),
            ("mode", ToolParamSchema::enum_values(&["classes", "bins"])),
            ("bins", int()),
            ("percent", ToolParamSchema::bool()),
            ("zone_band", int()),
            ("value_band", int()),
            ("long_output", table_out()),
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
