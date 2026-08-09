//! New geospatial tools that extend the `whitebox_next_gen` suite.
//!
//! Each tool implements [`wbcore::Tool`] (the same trait whitebox's own tools
//! use), so they plug into the registry alongside `register_default_tools` and
//! are exposed through the WASI runner exactly like the built-in tools.
//!
//! Add a new tool by creating a module with a `Tool` impl and pushing it in
//! [`geolibre_tools`].

mod add_surface_information;
mod adjust_3d_z;
mod adjust_stream_to_raster;
mod aggregate_points;
mod aggregate_polygons;
mod apportion_polygon;
mod assign_projection;
mod attribute_uncertainty;
mod block_statistics;
mod boundary_clean;
mod buffer_3d;
mod build_balanced_zones;
mod build_seamlines;
mod calculate_adjacent_fields;
mod calculate_distance_band;
mod calculate_grid_convergence_angle;
mod calculate_missing_z_values;
mod calculate_motion_statistics;
mod calculate_transformation_errors;
mod calculate_utm_zone;
mod cartogram;
mod causal_inference_analysis;
mod cell_records_to_sectors;
mod central_feature;
mod classification_accuracy_assessment;
mod collapse_dual_lines_to_centerline;
mod collect_events;
mod colocation_analysis;
mod combine;
mod common;
mod compare_spatial_weights;
mod construct_sight_lines;
mod corridor;
mod count_overlapping_features;
mod create_cartographic_partitions;
mod create_overpass;
mod create_underpass;
mod create_routes;
mod create_spatially_balanced_points;
mod cul_de_sac_masks;
mod cut_fill;
mod delineate_built_up_areas;
mod delineate_depressions;
mod delineate_mounts;
mod dem_filter;
mod dendrogram;
mod densify_sampling_network;
mod detect_feature_changes;
mod detect_image_anomalies;
mod diffusion_interpolation_with_barriers;
mod directional_distribution;
mod directional_trend;
mod dissolve_route_events;
mod eighty_twenty_analysis;
mod eliminate_polygon_part;
mod eliminate_polygons;
mod emerging_hot_spot_analysis;
mod empirical_bayesian_kriging;
mod enforce_river_monotonicity;
mod estimate_time_to_event;
mod evaluate_bin_sizes;
mod excel_to_table;
mod euclidean_direction;
mod expand_shrink;
mod exploratory_interpolation;
mod exploratory_regression;
mod extract_locations_from_text;
mod extract_sinks;
mod feature_vertices_to_points;
mod features_to_gpx;
mod features_to_gtfs;
mod fill;
mod fill_missing_values;
mod fill_spill_merge;
mod fill_spill_merge_core;
mod find_dwell_locations;
mod find_identical;
mod find_meeting_locations;
mod find_space_time_matches;
mod flip_line;
mod focal_flow;
mod focal_statistics;
mod forest_based_forecast;
mod compute_sar_indices;
mod convert_sar_units;
mod flatten_interferogram;
mod frequency_comparison;
mod fuzzy_overlay;
mod gaussian_geostatistical_simulations;
mod generate_breach_lines;
mod generate_near_table;
mod generate_od_links;
mod generate_points_along_3d_lines;
mod generate_spatial_weights_matrix;
mod generate_subset_polygons;
mod generate_transects_along_lines;
mod geographically_weighted_regression;
mod geoparquet_io;
mod getis_ord_general_g;
mod gpx_to_features;
mod graphic_buffer;
mod h3_polyfill;
mod h3_to_vector;
mod hdbscan;
mod hilbert;
mod idw_3d;
mod incremental_spatial_autocorrelation;
mod inside_3d;
mod integrate;
mod interpolate_from_spatiotemporal_points;
mod interpolate_shape;
mod kernel_interpolation_with_barriers;
mod kml_to_features;
mod lidar_common;
mod line_of_sight;
mod local_polynomial_interpolation;
mod locate_lines_along_routes;
mod matched_filter_target_detection;
mod maximum_likelihood_classification;
mod median_center;
mod merge_divided_roads;
mod mgwr;
mod minimum_bounding_volume;
mod multicriteria_overlay;
mod multilook;
mod multiple_ring_buffer;
mod near_3d;
mod neighborhood_summary_statistics;
mod non_maximum_suppression;
mod optics_clustering;
mod optimal_corridor_connections;
mod optimal_interpolation;
mod path_distance;
mod pmtiles;
mod pivot_table;
mod pmtiles_extract;
mod point_statistics;
mod points_to_line;
mod polygon_neighbors;
mod polygonize;
mod raster_normalize;
mod raster_to_h3;
mod raster_to_tiles;
mod reclassify_field;
mod reconstruct_tracks;
mod regions;
mod regularize_adjacent_building_footprint;
mod regularize_building_footprints;
mod remove_overlap_multiple;
mod render;
mod render_png;
mod render_vector_png;
mod reproject_raster;
mod resolve_building_conflicts;
mod ripleys_k;
mod rubbersheet_features;
mod sar_coherence;
mod similarity_search;
mod simplify_3d_line;
mod simplify_by_circular_arcs;
mod simplify_building;
mod simplify_shared_edges;
mod slice_raster;
mod smooth_natural_features;
mod smooth_shared_edges;
mod snap_tracks;
mod solar_radiation;
mod spatial_eigenvector_filtering;
mod spatially_constrained_multivariate_clustering;
mod spectral_index;
mod split_by_attributes;
mod split_by_features;
mod split_line_at_point;
mod stack_profile;
mod storage_capacity;
mod subdivide_polygon;
mod subset_features;
mod summarize_nearby;
mod summarize_categorical_raster;
mod summarize_percent_change;
mod summarize_within;
mod summary_statistics;
mod surface_volume;
mod tabulate_intersection;
mod thin_road_network;
mod time_series_clustering;
mod time_series_smoothing;
mod trace_proximity_events;
mod transform_features;
mod transform_fields;
mod transform_route_events;
mod trim_line;
mod cost_back_link;
mod difference_3d;
mod enclose_multipatch;
mod generate_network_swm;
mod intersect_3d;
mod intervisibility;
mod is_closed_3d;
mod mesh3d;
mod multipatch_footprint;
mod raster_stack;
mod simplify_by_tangent_segments;
mod union_3d;
mod unwrap_phase;
mod vector_common;
mod vector_convert;
mod vector_to_h3;
mod vector_to_pmtiles;
mod voxel_isosurface;
mod write_pmtiles;

mod sort_features;

mod calculate_central_meridian_and_parallels;
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

mod align_features;
mod analyze_changes_ccdc;
mod band_collection_statistics;
mod calculate_polygon_main_angle;
mod add_z_information;
mod args_common;
mod apply_radiometric_calibration;
mod cell_position_statistics;
mod cell_statistics;
mod collapse_road_detail;
mod compute_accuracy_for_object_detection;
mod contour_with_barriers;
mod convert_coordinate_notation;
mod create_spatial_sampling_locations;
mod darcy_flow;
mod detect_graphic_conflict;
mod detect_incidents;
mod dimension_reduction;
mod disperse_markers;
mod extract_scanned_features;
mod feature_outline_masks;
mod find_argument_statistics;
mod generalized_linear_regression;
mod geodetic_densify;
mod geotagged_photos_to_points;
mod grid_index_features;
mod gtfs_to_features;
mod hotspot_common;
mod identify_narrow_polygons;
mod interpolate_with_barriers;
mod intersecting_layers_masks;
mod kernel_density_ratio;
mod las_height_metrics;
mod line_density;
mod line_statistics;
mod local_bivariate_relationships;
mod merge_lines_by_pseudo_node;
mod multidimensional_anomaly;
mod multivariate_clustering;
mod optimized_hot_spot_analysis;
mod optimized_outlier_analysis;
mod pairwise_comparison_weights;
mod percentile_contours;
mod points_to_path;
mod porous_puff;
mod predict_using_trend_raster;
mod presence_only_prediction;
mod propagate_displacement;
mod repair_geometry;
mod resolve_road_conflicts;
mod space_time_kernel_density;
mod spatial_association_between_zones;
mod strip_map_index_features;
mod table_to_geometry;
mod time_series_cross_correlation;
mod time_series_forecast;
mod topo_to_raster;
mod zonal_characterization;
mod zonal_fill;
mod zonal_geometry;
mod zonal_histogram;

mod calculate_transit_service_frequency;
mod feature_to_line;
mod group_by_proximity;
mod hot_spot_analysis_comparison;
mod polygon_volume;
mod rescale_by_function;
mod split_raster;
mod surface_parameters;

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
        Box::new(build_seamlines::BuildSeamlinesTool),
        Box::new(merge_divided_roads::MergeDividedRoadsTool),
        Box::new(spatially_constrained_multivariate_clustering::SpatiallyConstrainedMultivariateClusteringTool),
        Box::new(optimal_corridor_connections::OptimalCorridorConnectionsTool),
        Box::new(adjust_stream_to_raster::AdjustStreamToRasterTool),
        Box::new(generate_breach_lines::GenerateBreachLinesTool),
        Box::new(idw_3d::Idw3dTool),
        Box::new(local_polynomial_interpolation::LocalPolynomialInterpolationTool),
        Box::new(kernel_interpolation_with_barriers::KernelInterpolationWithBarriersTool),
        Box::new(minimum_bounding_volume::MinimumBoundingVolumeTool),
        Box::new(voxel_isosurface::VoxelIsosurfaceTool),
        Box::new(inside_3d::Inside3dTool),
        Box::new(cost_back_link::CostBackLinkTool),
        Box::new(simplify_by_tangent_segments::SimplifyByTangentSegmentsTool),
        Box::new(difference_3d::Difference3dTool),
        Box::new(add_z_information::AddZInformationTool),
        Box::new(enclose_multipatch::EncloseMultipatchTool),
        Box::new(generate_network_swm::GenerateNetworkSwmTool),
        Box::new(intersect_3d::Intersect3dTool),
        Box::new(intervisibility::IntervisibilityTool),
        Box::new(is_closed_3d::IsClosed3dTool),
        Box::new(multipatch_footprint::MultipatchFootprintTool),
        Box::new(union_3d::Union3dTool),
        Box::new(unwrap_phase::UnwrapPhaseTool),
        Box::new(stack_profile::StackProfileTool),
        Box::new(generate_points_along_3d_lines::GeneratePointsAlong3dLinesTool),
        Box::new(split_by_features::SplitByFeaturesTool),
        Box::new(construct_sight_lines::ConstructSightLinesTool),
        Box::new(summarize_percent_change::SummarizePercentChangeTool),
        Box::new(summarize_categorical_raster::SummarizeCategoricalRasterTool),
        Box::new(dissolve_route_events::DissolveRouteEventsTool),
        Box::new(calculate_transformation_errors::CalculateTransformationErrorsTool),
        Box::new(subset_features::SubsetFeaturesTool),
        Box::new(pivot_table::PivotTableTool),
        Box::new(classification_accuracy_assessment::ClassificationAccuracyAssessmentTool),
        Box::new(diffusion_interpolation_with_barriers::DiffusionInterpolationWithBarriersTool),
        Box::new(compare_spatial_weights::CompareSpatialWeightsTool),
        Box::new(spatial_eigenvector_filtering::SpatialEigenvectorFilteringTool),
        Box::new(enforce_river_monotonicity::EnforceRiverMonotonicityTool),
        Box::new(calculate_grid_convergence_angle::CalculateGridConvergenceAngleTool),
        Box::new(maximum_likelihood_classification::MaximumLikelihoodClassificationTool),
        Box::new(focal_statistics::FocalStatisticsTool),
        Box::new(multicriteria_overlay::MulticriteriaOverlayTool),
        Box::new(multilook::MultilookTool),
        Box::new(surface_volume::SurfaceVolumeTool),
        Box::new(generate_spatial_weights_matrix::GenerateSpatialWeightsMatrixTool),
        Box::new(generate_subset_polygons::GenerateSubsetPolygonsTool),
        Box::new(calculate_distance_band::CalculateDistanceBandTool),
        Box::new(point_statistics::PointStatisticsTool),
        Box::new(kml_to_features::KmlToFeaturesTool),
        Box::new(graphic_buffer::GraphicBufferTool),
        Box::new(features_to_gpx::FeaturesToGpxTool),
        Box::new(reclassify_field::ReclassifyFieldTool),
        Box::new(points_to_line::PointsToLineTool),
        Box::new(collect_events::CollectEventsTool),
        Box::new(slice_raster::SliceRasterTool),
        Box::new(gpx_to_features::GpxToFeaturesTool),
        Box::new(summarize_within::SummarizeWithinTool),
        Box::new(summary_statistics::SummaryStatisticsTool),
        Box::new(median_center::MedianCenterTool),
        Box::new(flip_line::FlipLineTool),
        Box::new(focal_flow::FocalFlowTool),
        Box::new(cell_records_to_sectors::CellRecordsToSectorsTool),
        Box::new(estimate_time_to_event::EstimateTimeToEventTool),
        Box::new(transform_route_events::TransformRouteEventsTool),
        Box::new(eighty_twenty_analysis::EightyTwentyAnalysisTool),
        Box::new(extract_locations_from_text::ExtractLocationsFromTextTool),
        Box::new(create_cartographic_partitions::CreateCartographicPartitionsTool),
        Box::new(calculate_utm_zone::CalculateUtmZoneTool),
        Box::new(locate_lines_along_routes::LocateLinesAlongRoutesTool),
        Box::new(optimal_interpolation::OptimalInterpolationTool),
        Box::new(calculate_adjacent_fields::CalculateAdjacentFieldsTool),
        Box::new(dendrogram::DendrogramTool),
        Box::new(densify_sampling_network::DensifySamplingNetworkTool),
        Box::new(feature_vertices_to_points::FeatureVerticesToPointsTool),
        Box::new(create_overpass::CreateOverpassTool),
        Box::new(create_underpass::CreateUnderpassTool),
        Box::new(features_to_gtfs::FeaturesToGtfsTool),
        Box::new(adjust_3d_z::Adjust3dZTool),
        Box::new(attribute_uncertainty::AttributeUncertaintyTool),
        Box::new(split_line_at_point::SplitLineAtPointTool),
        Box::new(directional_trend::DirectionalTrendTool),
        Box::new(evaluate_bin_sizes::EvaluateBinSizesTool),
        Box::new(excel_to_table::ExcelToTableTool),
        Box::new(add_surface_information::AddSurfaceInformationTool),
        Box::new(create_routes::CreateRoutesTool),
        Box::new(exploratory_interpolation::ExploratoryInterpolationTool),
        Box::new(getis_ord_general_g::GetisOrdGeneralGTool),
        Box::new(matched_filter_target_detection::MatchedFilterTargetDetectionTool),
        Box::new(trim_line::TrimLineTool),
        Box::new(fill_missing_values::FillMissingValuesTool),
        Box::new(time_series_smoothing::TimeSeriesSmoothingTool),
        Box::new(combine::CombineTool),
        Box::new(assign_projection::AssignProjectionRasterTool),
        Box::new(assign_projection::AssignProjectionVectorTool),
        Box::new(assign_projection::AssignProjectionLidarTool),
        Box::new(raster_normalize::RasterNormalizeTool),
        Box::new(dem_filter::DemFilterTool),
        Box::new(extract_sinks::ExtractSinksTool),
        Box::new(fill_spill_merge::FillSpillMergeTool),
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
        Box::new(regularize_adjacent_building_footprint::RegularizeAdjacentBuildingFootprintTool),
        Box::new(smooth_natural_features::SmoothNaturalFeaturesTool),
        Box::new(eliminate_polygons::EliminatePolygonsTool),
        Box::new(eliminate_polygon_part::EliminatePolygonPartTool),
        Box::new(simplify_3d_line::Simplify3dLineTool),
        Box::new(simplify_building::SimplifyBuildingTool),
        Box::new(simplify_by_circular_arcs::SimplifyByCircularArcsTool),
        Box::new(simplify_shared_edges::SimplifySharedEdgesTool),
        Box::new(smooth_shared_edges::SmoothSharedEdgesTool),
        Box::new(emerging_hot_spot_analysis::EmergingHotSpotAnalysisTool),
        Box::new(line_of_sight::LineOfSightTool),
        Box::new(corridor::CorridorTool),
        Box::new(interpolate_from_spatiotemporal_points::InterpolateFromSpatiotemporalPointsTool),
        Box::new(interpolate_shape::InterpolateShapeTool),
        Box::new(collapse_dual_lines_to_centerline::CollapseDualLinesToCenterlineTool),
        Box::new(count_overlapping_features::CountOverlappingFeaturesTool),
        Box::new(non_maximum_suppression::NonMaximumSuppressionTool),
        Box::new(subdivide_polygon::SubdividePolygonTool),
        Box::new(generate_transects_along_lines::GenerateTransectsAlongLinesTool),
        Box::new(polygon_neighbors::PolygonNeighborsTool),
        Box::new(split_by_attributes::SplitByAttributesTool),
        Box::new(incremental_spatial_autocorrelation::IncrementalSpatialAutocorrelationTool),
        Box::new(apportion_polygon::ApportionPolygonTool),
        Box::new(central_feature::CentralFeatureTool),
        Box::new(expand_shrink::ExpandShrinkTool),
        Box::new(euclidean_direction::EuclideanDirectionTool),
        Box::new(reconstruct_tracks::ReconstructTracksTool),
        Box::new(solar_radiation::SolarRadiationTool),
        Box::new(hdbscan::HdbscanTool),
        Box::new(optics_clustering::OpticsClusteringTool),
        Box::new(colocation_analysis::ColocationAnalysisTool),
        Box::new(similarity_search::SimilaritySearchTool),
        Box::new(detect_feature_changes::DetectFeatureChangesTool),
        Box::new(integrate::IntegrateTool),
        Box::new(rubbersheet_features::RubbersheetFeaturesTool),
        Box::new(snap_tracks::SnapTracksTool),
        Box::new(remove_overlap_multiple::RemoveOverlapMultipleTool),
        Box::new(forest_based_forecast::ForestBasedForecastTool),
        Box::new(compute_sar_indices::ComputeSarIndicesTool),
        Box::new(convert_sar_units::ConvertSarUnitsTool),
        Box::new(flatten_interferogram::FlattenInterferogramTool),
        Box::new(frequency_comparison::FrequencyComparisonTool),
        Box::new(fuzzy_overlay::FuzzyOverlayTool),
        Box::new(aggregate_points::AggregatePointsTool),
        Box::new(generate_od_links::GenerateOdLinksTool),
        Box::new(generate_near_table::GenerateNearTableTool),
        Box::new(near_3d::Near3dTool),
        Box::new(neighborhood_summary_statistics::NeighborhoodSummaryStatisticsTool),
        Box::new(storage_capacity::StorageCapacityTool),
        Box::new(find_space_time_matches::FindSpaceTimeMatchesTool),
        Box::new(create_spatially_balanced_points::CreateSpatiallyBalancedPointsTool),
        Box::new(find_dwell_locations::FindDwellLocationsTool),
        Box::new(find_identical::FindIdenticalTool),
        Box::new(path_distance::PathDistanceTool),
        Box::new(time_series_clustering::TimeSeriesClusteringTool),
        Box::new(trace_proximity_events::TraceProximityEventsTool),
        Box::new(find_meeting_locations::FindMeetingLocationsTool),
        Box::new(detect_image_anomalies::DetectImageAnomaliesTool),
        Box::new(resolve_building_conflicts::ResolveBuildingConflictsTool),
        Box::new(delineate_built_up_areas::DelineateBuiltUpAreasTool),
        Box::new(aggregate_polygons::AggregatePolygonsTool),
        Box::new(multiple_ring_buffer::MultipleRingBufferTool),
        Box::new(directional_distribution::DirectionalDistributionTool),
        Box::new(tabulate_intersection::TabulateIntersectionTool),
        Box::new(summarize_nearby::SummarizeNearbyTool),
        Box::new(cul_de_sac_masks::CulDeSacMasksTool),
        Box::new(cut_fill::CutFillTool),
        Box::new(ripleys_k::RipleysKTool),
        Box::new(geographically_weighted_regression::GeographicallyWeightedRegressionTool),
        Box::new(mgwr::MgwrTool),
        Box::new(buffer_3d::Buffer3dTool),
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
        Box::new(block_statistics::BlockStatisticsTool),
        Box::new(calculate_missing_z_values::CalculateMissingZValuesTool),
        Box::new(calculate_motion_statistics::CalculateMotionStatisticsTool),
        Box::new(sort_features::SortFeaturesTool),
        Box::new(calculate_central_meridian_and_parallels::CalculateCentralMeridianAndParallelsTool),
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
        Box::new(porous_puff::PorousPuffTool),
        Box::new(predict_using_trend_raster::PredictUsingTrendRasterTool),
        Box::new(time_series_cross_correlation::TimeSeriesCrossCorrelationTool),
        Box::new(generalized_linear_regression::GeneralizedLinearRegressionTool),
        Box::new(interpolate_with_barriers::InterpolateWithBarriersTool),
        Box::new(convert_coordinate_notation::ConvertCoordinateNotationTool),
        Box::new(repair_geometry::RepairGeometryTool),
        Box::new(grid_index_features::GridIndexFeaturesTool),
        Box::new(local_bivariate_relationships::LocalBivariateRelationshipsTool),
        Box::new(dimension_reduction::DimensionReductionTool),
        Box::new(feature_outline_masks::FeatureOutlineMasksTool),
        Box::new(intersecting_layers_masks::IntersectingLayersMasksTool),
        Box::new(line_density::LineDensityTool),
        Box::new(pairwise_comparison_weights::PairwiseComparisonWeightsTool),
        Box::new(kernel_density_ratio::KernelDensityRatioTool),
        Box::new(detect_incidents::DetectIncidentsTool),
        Box::new(find_argument_statistics::FindArgumentStatisticsTool),
        Box::new(las_height_metrics::LasHeightMetricsTool),
        Box::new(apply_radiometric_calibration::ApplyRadiometricCalibrationTool),
        Box::new(cell_position_statistics::CellPositionStatisticsTool),
        Box::new(cell_statistics::CellStatisticsTool),
        Box::new(multidimensional_anomaly::MultidimensionalAnomalyTool),
        Box::new(propagate_displacement::PropagateDisplacementTool),
        Box::new(empirical_bayesian_kriging::EmpiricalBayesianKrigingTool),
        Box::new(gaussian_geostatistical_simulations::GaussianGeostatisticalSimulationsTool),
        Box::new(exploratory_regression::ExploratoryRegressionTool),
        Box::new(causal_inference_analysis::CausalInferenceAnalysisTool),
        Box::new(align_features::AlignFeaturesTool),
        Box::new(multivariate_clustering::MultivariateClusteringTool),
        Box::new(table_to_geometry::TableToGeometryTool),
        Box::new(transform_features::TransformFeaturesTool),
        Box::new(transform_fields::TransformFieldsTool),
        Box::new(detect_graphic_conflict::DetectGraphicConflictTool),
        Box::new(disperse_markers::DisperseMarkersTool),
        Box::new(geodetic_densify::GeodeticDensifyTool),
        Box::new(strip_map_index_features::StripMapIndexFeaturesTool),
        Box::new(zonal_histogram::ZonalHistogramTool),
        Box::new(points_to_path::PointsToPathTool),
        Box::new(extract_scanned_features::ExtractScannedFeaturesTool),
        Box::new(gtfs_to_features::GtfsToFeaturesTool),
        Box::new(create_spatial_sampling_locations::CreateSpatialSamplingLocationsTool),
        Box::new(compute_accuracy_for_object_detection::ComputeAccuracyForObjectDetectionTool),
        Box::new(contour_with_barriers::ContourWithBarriersTool),
        Box::new(percentile_contours::PercentileContoursTool),
        Box::new(spatial_association_between_zones::SpatialAssociationBetweenZonesTool),
        Box::new(merge_lines_by_pseudo_node::MergeLinesByPseudoNodeTool),
        Box::new(identify_narrow_polygons::IdentifyNarrowPolygonsTool),
        Box::new(zonal_geometry::ZonalGeometryTool),
        Box::new(zonal_characterization::ZonalCharacterizationTool),
        Box::new(zonal_fill::ZonalFillTool),
        Box::new(calculate_polygon_main_angle::CalculatePolygonMainAngleTool),
        Box::new(band_collection_statistics::BandCollectionStatisticsTool),
        Box::new(line_statistics::LineStatisticsTool),
        Box::new(optimized_hot_spot_analysis::OptimizedHotSpotAnalysisTool),
        Box::new(optimized_outlier_analysis::OptimizedOutlierAnalysisTool),
        Box::new(polygon_volume::PolygonVolumeTool),
        Box::new(hot_spot_analysis_comparison::HotSpotAnalysisComparisonTool),
        Box::new(group_by_proximity::GroupByProximityTool),
        Box::new(feature_to_line::FeatureToLineTool),
        Box::new(split_raster::SplitRasterTool),
        Box::new(surface_parameters::SurfaceParametersTool),
        Box::new(sar_coherence::SarCoherenceTool),
        Box::new(rescale_by_function::RescaleByFunctionTool),
        Box::new(calculate_transit_service_frequency::CalculateTransitServiceFrequencyTool),
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
    let json_out = || ToolParamSchema::output(ToolDatasetSchema::Json);
    let lidar_in = ToolParamSchema::input_lidar;
    let lidar_out = || ToolParamSchema::output(ToolDatasetSchema::Lidar);
    let int = ToolParamSchema::scalar_integer;
    let float = ToolParamSchema::scalar_float;
    let colormaps =
        || ToolParamSchema::enum_values(&["viridis", "magma", "turbo", "terrain", "grayscale"]);

    let map = match tool_id {
        "build_seamlines" => schemas(&[
            ("inputs", ToolParamSchema::string()),
            ("output", vector_out()),
            ("footprints", vector_in()),
            ("method", ToolParamSchema::enum_values(&["voronoi", "order", "radiometry", "edge_detection"])),
            ("sort_ascending", ToolParamSchema::bool()),
            ("cell_size", float()),
            ("band", int()),
            ("min_region_size", int()),
            ("blend_width", float()),
            ("blend_type", ToolParamSchema::enum_values(&["both", "inside", "outside"])),
        ]),
        "merge_divided_roads" => schemas(&[
            ("input", vector_in()),
            ("merge_field", ToolParamSchema::string()),
            ("merge_distance", float()),
            ("output", vector_out()),
            ("output_displacement_features", vector_out()),
            ("character_field", ToolParamSchema::string()),
            ("output_table", table_out()),
        ]),
        "spatially_constrained_multivariate_clustering" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("analysis_fields", ToolParamSchema::string()),
            ("number_of_clusters", int()),
            ("neighborhood", ToolParamSchema::enum_values(&["contiguity_edges", "contiguity_edges_corners", "knn"])),
            ("number_of_neighbors", int()),
            ("constraint", ToolParamSchema::enum_values(&["none", "feature_count", "attribute_value"])),
            ("constraint_field", ToolParamSchema::string()),
            ("min_constraint", float()),
            ("max_constraint", float()),
            ("scale_data", ToolParamSchema::bool()),
            ("output_table", table_out()),
        ]),
        "optimal_corridor_connections" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("cost_raster", raster_in()),
            ("barriers", vector_in()),
            ("corridor_width", float()),
            ("output_lines", vector_out()),
            ("neighbor_option", ToolParamSchema::enum_values(&["spanning_tree", "all_pairs"])),
            ("cell_size", float()),
        ]),
        "adjust_stream_to_raster" => schemas(&[
            ("input", vector_in()),
            ("dem", raster_in()),
            ("output", vector_out()),
            ("snap_distance", float()),
            ("channel_threshold", float()),
            ("remove_disconnected", ToolParamSchema::bool()),
            ("output_stream_raster", raster_out()),
            ("output_flow_direction", raster_out()),
            ("output_split_points", vector_out()),
        ]),
        "generate_breach_lines" => schemas(&[
            ("input", vector_in()),
            ("dem", raster_in()),
            ("output", vector_out()),
            ("connection_points", vector_in()),
            ("method", ToolParamSchema::enum_values(&["minimum_breaching_cost", "shortest_path", "minimum_elevation_change"])),
            ("max_length", float()),
        ]),
        "idw_3d" => schemas(&[
            ("input", vector_in()),
            ("value_field", ToolParamSchema::string()),
            ("output", raster_out()),
            ("z_field", ToolParamSchema::string()),
            ("power", float()),
            ("elev_inflation_factor", float()),
            ("x_spacing", float()),
            ("y_spacing", float()),
            ("z_spacing", float()),
            ("z_min", float()),
            ("z_max", float()),
            ("neighbors", int()),
            ("search_radius", float()),
            ("output_cv_features", vector_out()),
            ("epsg", int()),
        ]),
        "kernel_interpolation_with_barriers" => schemas(&[
            ("input", vector_in()),
            ("z_field", ToolParamSchema::string()),
            ("output", raster_out()),
            ("output_error", raster_out()),
            ("barriers", vector_in()),
            ("cell_size", float()),
            (
                "kernel",
                ToolParamSchema::enum_values(&[
                    "exponential",
                    "gaussian",
                    "quartic",
                    "epanechnikov",
                    "polynomial5",
                    "constant",
                ]),
            ),
            ("bandwidth", float()),
            ("power", int()),
            ("ridge", float()),
            ("max_neighbors", int()),
        ]),
        "local_polynomial_interpolation" => schemas(&[
            ("input", vector_in()),
            ("z_field", ToolParamSchema::string()),
            ("output", raster_out()),
            ("cell_size", float()),
            ("order", int()),
            ("kernel", ToolParamSchema::enum_values(&["exponential", "gaussian", "quartic", "epanechnikov", "fifth_order", "constant"])),
            ("bandwidth", float()),
            ("neighbors", int()),
            ("weight_field", ToolParamSchema::string()),
            ("condition_number", float()),
            ("output_type", ToolParamSchema::enum_values(&["prediction", "standard_error", "condition_number"])),
            ("epsg", int()),
        ]),
        "multipatch_footprint" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("group_field", ToolParamSchema::string()),
            ("simplify_tolerance", float()),
        ]),
        "add_z_information" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("properties", ToolParamSchema::string()),
            ("noise_filtering", float()),
        ]),
        "enclose_multipatch" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("skip_closed", ToolParamSchema::bool()),
            ("drop_failed", ToolParamSchema::bool()),
        ]),
        "intervisibility" => schemas(&[
            ("input", vector_in()),
            // Comma-separated layer list, as stack_profile.surfaces does.
            ("obstructions", ToolParamSchema::string()),
            ("output", vector_out()),
            ("visible_field", ToolParamSchema::string()),
            ("observer_offset", float()),
            ("target_offset", float()),
            ("visible_only", ToolParamSchema::bool()),
        ]),
        "simplify_by_tangent_segments" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("max_offset", float()),
            ("anchor_points", vector_in()),
            ("anchor_tolerance", float()),
            ("min_run", int()),
            ("preserve_endpoints", ToolParamSchema::bool()),
        ]),
        "cost_back_link" => schemas(&[
            // Dual-typed: a source raster or a point vector layer.
            ("source", ToolParamSchema::string()),
            ("cost", raster_in()),
            ("surface", raster_in()),
            ("output", raster_out()),
            ("out_distance", raster_out()),
            ("max_distance", float()),
        ]),
        "generate_network_swm" => schemas(&[
            ("input", vector_in()),
            ("network", vector_in()),
            ("output", table_out()),
            ("id_field", ToolParamSchema::string()),
            ("impedance_field", ToolParamSchema::string()),
            ("distance_cutoff", float()),
            ("max_neighbors", int()),
            (
                "conceptualization",
                ToolParamSchema::enum_values(&["fixed", "inverse"]),
            ),
            ("exponent", float()),
            ("row_standardization", ToolParamSchema::bool()),
            ("search_tolerance", float()),
        ]),
        "intersect_3d" => schemas(&[
            ("input", vector_in()),
            ("input2", vector_in()),
            ("output", table_out()),
            ("mode", ToolParamSchema::enum_values(&["table", "solid"])),
            ("resolution", int()),
        ]),
        "difference_3d" => schemas(&[
            ("input", vector_in()),
            ("subtract", vector_in()),
            ("output", table_out()),
            ("resolution", int()),
            ("keep_geometry", ToolParamSchema::bool()),
        ]),
        "is_closed_3d" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("closed_only", ToolParamSchema::bool()),
        ]),
        "union_3d" => schemas(&[
            ("input", vector_in()),
            ("output", table_out()),
            ("group_field", ToolParamSchema::string()),
            ("output_overlap_table", table_out()),
            ("resolution", int()),
        ]),
        "inside_3d" => schemas(&[
            ("target", vector_in()),
            ("container", vector_in()),
            ("output", table_out()),
            ("mode", ToolParamSchema::enum_values(&["simple", "complex"])),
            ("output_features", vector_out()),
        ]),
        "voxel_isosurface" => schemas(&[
            ("input", raster_in()),
            ("values", ToolParamSchema::string()),
            ("output", vector_out()),
            ("z_min", float()),
            ("z_spacing", float()),
            ("close_boundaries", ToolParamSchema::bool()),
            ("smooth", int()),
        ]),
        "minimum_bounding_volume" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("z_value", ToolParamSchema::string()),
            ("geometry_type", ToolParamSchema::enum_values(&["convex_hull", "sphere", "envelope"])),
            ("group_option", ToolParamSchema::enum_values(&["none", "all", "list"])),
            ("group_field", ToolParamSchema::string()),
            ("mbv_fields", ToolParamSchema::bool()),
        ]),
        "stack_profile" => schemas(&[
            ("input", vector_in()),
            ("surfaces", ToolParamSchema::string()),
            ("output", table_out()),
            ("sample_distance", float()),
            ("line_id_field", ToolParamSchema::string()),
            ("method", ToolParamSchema::enum_values(&["bilinear", "nearest"])),
        ]),
        "generate_points_along_3d_lines" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("method", ToolParamSchema::enum_values(&["distance", "percentage", "distance_field"])),
            ("distance", float()),
            ("percentage", float()),
            ("distance_field", ToolParamSchema::string()),
            ("include_end_points", ToolParamSchema::bool()),
            ("add_chainage", ToolParamSchema::bool()),
        ]),
        "split_by_features" => schemas(&[
            ("input", vector_in()),
            ("split_features", vector_in()),
            ("split_field", ToolParamSchema::string()),
            ("output_dir", file_out()),
            ("output_format", ToolParamSchema::enum_values(&["geojson", "shp", "gpkg", "parquet", "csv"])),
        ]),
        "construct_sight_lines" => schemas(&[
            ("observers", vector_in()),
            ("targets", vector_in()),
            ("output", vector_out()),
            ("observer_height_field", ToolParamSchema::string()),
            ("target_height_field", ToolParamSchema::string()),
            ("join_field", ToolParamSchema::string()),
            ("sample_distance", float()),
            ("output_direction", ToolParamSchema::bool()),
            ("distance_method", ToolParamSchema::enum_values(&["2d", "3d"])),
        ]),
        "summarize_categorical_raster" => schemas(&[
            ("input", raster_in()),
            ("output", table_out()),
            ("aoi", vector_in()),
            ("aoi_id_field", ToolParamSchema::string()),
            (
                "area_units",
                ToolParamSchema::enum_values(&["map_units", "hectares", "square_kilometers"]),
            ),
            ("include_nodata", ToolParamSchema::bool()),
            ("bands", ToolParamSchema::string()),
        ]),
        "summarize_percent_change" => schemas(&[
            ("input", vector_in()),
            ("current_features", vector_in()),
            ("previous_features", vector_in()),
            ("output", vector_out()),
            ("search_radius", float()),
        ]),
        "dissolve_route_events" => schemas(&[
            ("input", vector_in()),
            ("output", table_out()),
            ("route_id_field", ToolParamSchema::string()),
            ("from_measure_field", ToolParamSchema::string()),
            ("to_measure_field", ToolParamSchema::string()),
            ("dissolve_fields", ToolParamSchema::string()),
            ("mode", ToolParamSchema::enum_values(&["dissolve", "concatenate"])),
            ("tolerance", float()),
            ("separator", ToolParamSchema::string()),
        ]),
        "calculate_transformation_errors" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("method", ToolParamSchema::enum_values(&["affine", "similarity", "projective"])),
            ("keep_geometry", ToolParamSchema::bool()),
        ]),
        "subset_features" => schemas(&[
            ("input", vector_in()),
            ("output_training", vector_out()),
            ("output_test", vector_out()),
            ("size", float()),
            ("size_method", ToolParamSchema::enum_values(&["percentage", "absolute"])),
            ("seed", int()),
            ("group_field", ToolParamSchema::string()),
        ]),
        "pivot_table" => schemas(&[
            ("input", vector_in()),
            ("output", table_out()),
            ("fields", ToolParamSchema::string()),
            ("pivot_field", ToolParamSchema::string()),
            ("value_field", ToolParamSchema::string()),
            ("aggregate", ToolParamSchema::enum_values(&["first", "sum", "mean", "min", "max", "count"])),
        ]),
        "classification_accuracy_assessment" => schemas(&[
            ("points", vector_in()),
            ("class_field", ToolParamSchema::string()),
            ("input", raster_in()),
            ("classified_field", ToolParamSchema::string()),
            ("band", int()),
            ("output", vector_out()),
        ]),
        "diffusion_interpolation_with_barriers" => schemas(&[
            ("input", vector_in()),
            ("z_field", ToolParamSchema::string()),
            ("output", raster_out()),
            ("cell_size", float()),
            ("barriers", vector_in()),
            ("bandwidth", float()),
            ("number_iterations", int()),
            ("weight_field", ToolParamSchema::string()),
            ("epsg", int()),
        ]),
        "compare_spatial_weights" => schemas(&[
            ("input", vector_in()),
            ("output", table_out()),
            ("input_fields", ToolParamSchema::string()),
            ("id_field", ToolParamSchema::string()),
            ("methods", ToolParamSchema::string()),
            ("number_of_neighbors", int()),
            ("threshold_distance", float()),
        ]),
        "spatial_eigenvector_filtering" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            (
                "method",
                ToolParamSchema::enum_values(&[
                    "contiguity_edges",
                    "contiguity_edges_corners",
                    "knn",
                ]),
            ),
            ("number_of_neighbors", int()),
            ("min_autocorrelation", float()),
            ("max_components", int()),
        ]),
        "enforce_river_monotonicity" => schemas(&[
            ("input", vector_in()),
            ("surface", raster_in()),
            ("output", raster_out()),
            ("tolerance", float()),
            ("sample_distance", float()),
            ("band", int()),
        ]),
        "calculate_grid_convergence_angle" => schemas(&[
            ("input", vector_in()),
            ("angle_field", ToolParamSchema::string()),
            (
                "angle_type",
                ToolParamSchema::enum_values(&["geographic", "arithmetic"]),
            ),
            ("central_meridian", float()),
            ("output", vector_out()),
        ]),
        "maximum_likelihood_classification" => schemas(&[
            ("input", raster_in()),
            ("training", raster_in()),
            ("output", raster_out()),
            ("prob_output", raster_out()),
            (
                "a_priori",
                ToolParamSchema::enum_values(&["equal", "sample"]),
            ),
            ("reject_fraction", float()),
        ]),
        "focal_statistics" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            (
                "statistics",
                ToolParamSchema::enum_values(&[
                    "mean",
                    "majority",
                    "maximum",
                    "median",
                    "minimum",
                    "minority",
                    "percentile",
                    "range",
                    "std",
                    "sum",
                    "variety",
                ]),
            ),
            (
                "neighborhood",
                ToolParamSchema::enum_values(&["rectangle", "circle", "annulus", "wedge"]),
            ),
            ("width", int()),
            ("height", int()),
            ("radius", float()),
            ("inner_radius", float()),
            ("start_angle", float()),
            ("end_angle", float()),
            ("percentile_value", float()),
            ("ignore_nodata", ToolParamSchema::bool()),
            ("band", int()),
        ]),
        "multilook" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            ("complex", ToolParamSchema::bool()),
            (
                "input_domain",
                ToolParamSchema::enum_values(&["intensity", "amplitude"]),
            ),
            ("range_looks", int()),
            ("azimuth_looks", int()),
            ("auto_looks", ToolParamSchema::bool()),
            (
                "output_units",
                ToolParamSchema::enum_values(&["amplitude", "intensity", "db"]),
            ),
            ("statistic", ToolParamSchema::enum_values(&["mean", "median"])),
        ]),
        "multicriteria_overlay" => schemas(&[
            (
                "inputs",
                ToolParamSchema::input_multiple(ToolDatasetSchema::Raster),
            ),
            ("output", raster_out()),
            (
                "method",
                ToolParamSchema::enum_values(&[
                    "weighted_sum",
                    "weighted_geometric_mean",
                    "owa",
                    "topsis",
                ]),
            ),
            ("weights", ToolParamSchema::string()),
            ("order_weights", ToolParamSchema::string()),
            ("from_scale", float()),
            ("to_scale", float()),
        ]),
        "surface_volume" => schemas(&[
            ("input", raster_in()),
            ("reference_plane", float()),
            (
                "direction",
                ToolParamSchema::enum_values(&["above", "below", "both"]),
            ),
            ("band", int()),
            ("output", table_out()),
        ]),
        "generate_subset_polygons" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("output_points", vector_out()),
            ("min_points_per_subset", int()),
            ("max_points_per_subset", int()),
            (
                "coincident_points",
                ToolParamSchema::enum_values(&["single", "all"]),
            ),
            ("clip_to_hull", ToolParamSchema::bool()),
        ]),
        "generate_spatial_weights_matrix" => schemas(&[
            ("input", vector_in()),
            ("output", table_out()),
            ("id_field", ToolParamSchema::string()),
            (
                "method",
                ToolParamSchema::enum_values(&[
                    "knn",
                    "fixed_distance_band",
                    "inverse_distance",
                    "contiguity_edges",
                    "contiguity_edges_corners",
                    "delaunay",
                ]),
            ),
            ("number_of_neighbors", int()),
            ("threshold_distance", float()),
            ("exponent", float()),
            ("row_standardization", ToolParamSchema::bool()),
            ("snap_tolerance", float()),
        ]),
        "calculate_distance_band" => schemas(&[
            ("input", vector_in()),
            ("neighbors", int()),
            (
                "distance_method",
                ToolParamSchema::enum_values(&["euclidean", "manhattan"]),
            ),
            ("output", file_out()),
        ]),
        "point_statistics" => schemas(&[
            ("input", vector_in()),
            ("field", ToolParamSchema::string()),
            ("output", raster_out()),
            (
                "statistic",
                ToolParamSchema::enum_values(&[
                    "mean", "majority", "maximum", "median", "minimum", "minority", "range", "std",
                    "sum", "variety",
                ]),
            ),
            (
                "neighborhood",
                ToolParamSchema::enum_values(&["circle", "rectangle"]),
            ),
            ("radius", float()),
            ("cell_size", float()),
            ("epsg", int()),
        ]),
        "kml_to_features" => schemas(&[
            ("input", ToolParamSchema::input(ToolDatasetSchema::File)),
            ("points_output", vector_out()),
            ("lines_output", vector_out()),
            ("polygons_output", vector_out()),
        ]),
        "graphic_buffer" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("distance", float()),
            (
                "cap",
                ToolParamSchema::enum_values(&["round", "square", "butt"]),
            ),
            (
                "join",
                ToolParamSchema::enum_values(&["round", "miter", "bevel"]),
            ),
            ("miter_limit", float()),
            ("dissolve", ToolParamSchema::bool()),
        ]),
        "features_to_gpx" => schemas(&[
            ("input", vector_in()),
            ("output", file_out()),
            ("name_field", ToolParamSchema::string()),
            ("description_field", ToolParamSchema::string()),
            ("z_field", ToolParamSchema::string()),
            ("date_field", ToolParamSchema::string()),
        ]),
        "reclassify_field" => schemas(&[
            ("input", vector_in()),
            ("field", ToolParamSchema::string()),
            ("output", vector_out()),
            (
                "method",
                ToolParamSchema::enum_values(&[
                    "equal_interval",
                    "defined_interval",
                    "quantile",
                    "natural_breaks",
                    "std_dev",
                    "geometric_interval",
                ]),
            ),
            ("classes", int()),
            ("interval", float()),
            ("std_dev_interval", float()),
            ("class_field", ToolParamSchema::string()),
            ("break_field", ToolParamSchema::string()),
        ]),
        "points_to_line" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("line_field", ToolParamSchema::string()),
            ("sort_field", ToolParamSchema::string()),
            ("close_line", ToolParamSchema::bool()),
        ]),
        "collect_events" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("tolerance", float()),
        ]),
        "slice_raster" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            ("band", int()),
            ("number_zones", int()),
            (
                "slice_type",
                ToolParamSchema::enum_values(&[
                    "equal_interval",
                    "equal_area",
                    "natural_breaks",
                    "geometric_interval",
                    "std_dev",
                ]),
            ),
            ("base_output_zone", int()),
            ("class_interval_size", float()),
        ]),
        "gpx_to_features" => schemas(&[
            ("input", ToolParamSchema::input(ToolDatasetSchema::File)),
            ("points_output", vector_out()),
            ("lines_output", vector_out()),
        ]),
        "summary_statistics" => schemas(&[
            ("input", vector_in()),
            ("statistics", ToolParamSchema::string()),
            ("case_fields", ToolParamSchema::string()),
            ("output", table_out()),
        ]),
        "median_center" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("weight_field", ToolParamSchema::string()),
            ("case_field", ToolParamSchema::string()),
            ("attribute_fields", ToolParamSchema::string()),
        ]),
        "flip_line" => schemas(&[("input", vector_in()), ("output", vector_out())]),
        "focal_flow" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            ("threshold", float()),
            ("band", int()),
        ]),
        "cell_records_to_sectors" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            (
                "output_type",
                ToolParamSchema::enum_values(&["wedge", "line"]),
            ),
            ("azimuth_field", ToolParamSchema::string()),
            ("beamwidth_field", ToolParamSchema::string()),
            ("radius_field", ToolParamSchema::string()),
            ("azimuth", float()),
            ("beamwidth", float()),
            ("radius", float()),
            ("segments", int()),
        ]),
        "estimate_time_to_event" => schemas(&[
            ("input", vector_in()),
            ("age_field", ToolParamSchema::string()),
            ("event_field", ToolParamSchema::string()),
            ("stratify_field", ToolParamSchema::string()),
            ("output", vector_out()),
        ]),
        "transform_route_events" => schemas(&[
            ("source_events", vector_in()),
            ("source_routes", vector_in()),
            ("target_routes", vector_in()),
            ("output", vector_out()),
            ("route_id_field", ToolParamSchema::string()),
            ("measure_field", ToolParamSchema::string()),
            ("target_id_field", ToolParamSchema::string()),
            ("cluster_tolerance", float()),
        ]),
        "eighty_twenty_analysis" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("cluster_tolerance", float()),
            ("weight_field", ToolParamSchema::string()),
            ("threshold", float()),
        ]),
        "extract_locations_from_text" => schemas(&[
            ("input_text", ToolParamSchema::string()),
            ("notations", ToolParamSchema::string()),
            ("output", vector_out()),
        ]),
        "create_cartographic_partitions" => schemas(&[
            ("input", vector_in()),
            ("feature_count", int()),
            ("output", vector_out()),
        ]),
        "calculate_utm_zone" => schemas(&[
            ("input", vector_in()),
            ("zone_field", ToolParamSchema::string()),
            ("epsg_field", ToolParamSchema::string()),
            ("output", vector_out()),
        ]),

        "locate_lines_along_routes" => schemas(&[
            ("input_features", vector_in()),
            ("routes", vector_in()),
            ("route_id_field", ToolParamSchema::string()),
            ("tolerance", float()),
            ("output", vector_out()),
        ]),
        "optimal_interpolation" => schemas(&[
            ("background", raster_in()),
            ("observations", vector_in()),
            ("field", ToolParamSchema::string()),
            ("output", raster_out()),
            ("correlation_length", float()),
            ("background_error_variance", float()),
            ("obs_error_variance", float()),
            ("error_field", ToolParamSchema::string()),
            ("band", int()),
            ("max_obs", int()),
            ("analysis_error", raster_out()),
        ]),
        "calculate_adjacent_fields" => schemas(&[
            ("input", vector_in()),
            ("page_name_field", ToolParamSchema::string()),
            ("output", vector_out()),
            ("include_diagonal", ToolParamSchema::bool()),
            ("snap_tolerance", float()),
        ]),
        "dendrogram" => schemas(&[
            ("input", vector_in()),
            ("class_field", ToolParamSchema::string()),
            ("fields", ToolParamSchema::string()),
            ("output", table_out()),
            ("distance", ToolParamSchema::enum_values(&["variance", "mean_only"])),
            ("standardize", ToolParamSchema::bool()),
            ("output_text", file_out()),
        ]),
        "densify_sampling_network" => schemas(&[
            ("prediction_error", raster_in()),
            ("output", vector_out()),
            ("count", int()),
            ("inhibition_distance", float()),
            ("mask", vector_in()),
            ("band", int()),
        ]),
        "feature_vertices_to_points" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            (
                "point_location",
                ToolParamSchema::enum_values(&[
                    "ALL",
                    "START",
                    "END",
                    "BOTH_ENDS",
                    "MID",
                    "DANGLE",
                ]),
            ),
        ]),
        "create_underpass" => schemas(&[
            ("above", vector_in()),
            ("below", vector_in()),
            ("output", vector_out()),
            ("output_lines", vector_out()),
            ("margin_along", float()),
            ("margin_across", float()),
            ("min_angle", float()),
        ]),
        "create_overpass" => schemas(&[
            ("above", vector_in()),
            ("below", vector_in()),
            ("output", vector_out()),
            ("output_decoration", vector_out()),
            ("margin_along", float()),
            ("margin_across", float()),
            (
                "wing_type",
                ToolParamSchema::enum_values(&["none", "perpendicular", "parallel"]),
            ),
        ]),
        "features_to_gtfs" => schemas(&[
            ("stops_input", vector_in()),
            ("shapes_input", vector_in()),
            (
                "output_dir",
                ToolParamSchema::output(ToolDatasetSchema::File),
            ),
        ]),
        "adjust_3d_z" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("factor", float()),
            ("offset", float()),
            (
                "from_unit",
                ToolParamSchema::enum_values(&[
                    "meters",
                    "feet",
                    "us_feet",
                    "centimeters",
                    "millimeters",
                    "kilometers",
                    "miles",
                    "yards",
                    "inches",
                ]),
            ),
            (
                "to_unit",
                ToolParamSchema::enum_values(&[
                    "meters",
                    "feet",
                    "us_feet",
                    "centimeters",
                    "millimeters",
                    "kilometers",
                    "miles",
                    "yards",
                    "inches",
                ]),
            ),
        ]),
        "attribute_uncertainty" => schemas(&[
            ("input", vector_in()),
            ("value_field", ToolParamSchema::string()),
            ("error_field", ToolParamSchema::string()),
            ("lower_field", ToolParamSchema::string()),
            ("upper_field", ToolParamSchema::string()),
            (
                "distribution",
                ToolParamSchema::enum_values(&["NORMAL", "UNIFORM"]),
            ),
            ("iterations", int()),
            ("seed", int()),
            ("output", vector_out()),
        ]),
        "split_line_at_point" => schemas(&[
            ("input_lines", vector_in()),
            ("point_features", vector_in()),
            ("search_radius", float()),
            ("output", vector_out()),
        ]),
        "directional_trend" => schemas(&[
            ("input", vector_in()),
            ("field", ToolParamSchema::string()),
            ("azimuth", ToolParamSchema::string()),
            ("order", int()),
            ("output", vector_out()),
        ]),
        "evaluate_bin_sizes" => schemas(&[
            ("input", vector_in()),
            ("output", table_out()),
            ("bin_shape", ToolParamSchema::enum_values(&["hexagon", "square"])),
            ("sizes", ToolParamSchema::string()),
            ("steps", int()),
            ("analysis_field", ToolParamSchema::string()),
        ]),
        "excel_to_table" => schemas(&[
            ("input", ToolParamSchema::input(ToolDatasetSchema::File)),
            ("output", table_out()),
            ("sheet", ToolParamSchema::string()),
            ("cell_range", ToolParamSchema::string()),
            ("field_names_row", int()),
        ]),
        "add_surface_information" => schemas(&[
            ("input", vector_in()),
            ("surface", raster_in()),
            ("output", vector_out()),
            ("properties", ToolParamSchema::string()),
            ("sample_distance", float()),
            (
                "method",
                ToolParamSchema::enum_values(&["bilinear", "nearest"]),
            ),
            ("band", int()),
        ]),
        "create_routes" => schemas(&[
            ("input", vector_in()),
            ("route_id_field", ToolParamSchema::string()),
            ("output", vector_out()),
            (
                "measure_source",
                ToolParamSchema::enum_values(&["LENGTH", "ONE_FIELD", "TWO_FIELDS"]),
            ),
            ("from_measure_field", ToolParamSchema::string()),
            ("to_measure_field", ToolParamSchema::string()),
            ("measure_factor", float()),
            ("measure_offset", float()),
            ("ignore_gaps", ToolParamSchema::bool()),
            ("snap_tolerance", float()),
        ]),
        "exploratory_interpolation" => schemas(&[
            ("input", vector_in()),
            ("field", ToolParamSchema::string()),
            ("output", table_out()),
            ("methods", ToolParamSchema::string()),
            (
                "criterion",
                ToolParamSchema::enum_values(&["rmse", "mae", "me"]),
            ),
            ("power", float()),
            ("output_raster", raster_out()),
            ("cell_size", float()),
        ]),
        "getis_ord_general_g" => schemas(&[
            ("input", vector_in()),
            ("field", ToolParamSchema::string()),
            (
                "weights",
                ToolParamSchema::enum_values(&["distance_band", "k_nearest", "queen", "rook"]),
            ),
            ("distance_band", float()),
            ("k", int()),
            ("row_standardize", ToolParamSchema::bool()),
            ("output", table_out()),
        ]),
        "matched_filter_target_detection" => schemas(&[
            ("input", raster_in()),
            ("target_spectrum", ToolParamSchema::string()),
            ("output", raster_out()),
            ("method", ToolParamSchema::enum_values(&["cem", "ace"])),
            ("threshold", float()),
            ("mask_output", raster_out()),
        ]),
        "trim_line" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("dangle_length", float()),
            ("keep_short", ToolParamSchema::bool()),
            ("snap_tolerance", float()),
        ]),
        "fill_missing_values" => schemas(&[
            ("input", vector_in()),
            ("fill_field", ToolParamSchema::string()),
            ("output", vector_out()),
            (
                "estimator",
                ToolParamSchema::enum_values(&["mean", "median", "min", "max", "temporal_trend"]),
            ),
            (
                "neighbourhood",
                ToolParamSchema::enum_values(&["knn", "distance_band"]),
            ),
            ("k", int()),
            ("search_radius", float()),
            ("time_field", ToolParamSchema::string()),
            ("time_window", float()),
            ("flag_field", ToolParamSchema::string()),
        ]),
        "time_series_smoothing" => schemas(&[
            ("input", vector_in()),
            ("value_field", ToolParamSchema::string()),
            ("time_field", ToolParamSchema::string()),
            ("id_field", ToolParamSchema::string()),
            (
                "method",
                ToolParamSchema::enum_values(&["moving_average", "local_linear"]),
            ),
            ("window", int()),
            ("bandwidth", int()),
            (
                "alignment",
                ToolParamSchema::enum_values(&["backward", "centered", "forward"]),
            ),
            ("output_field", ToolParamSchema::string()),
            ("output", vector_out()),
        ]),
        "combine" => schemas(&[
            ("inputs", ToolParamSchema::string()),
            ("output", raster_out()),
            ("csv_output", table_out()),
            ("band", int()),
        ]),
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
            (
                "filter",
                ToolParamSchema::enum_values(&["mean", "median", "gaussian"]),
            ),
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
        "fill_spill_merge" => schemas(&[
            ("dem", raster_in()),
            ("water_level", float()),
            ("surface_water", raster_in()),
            ("ocean_level", float()),
            ("edge_outlet", ToolParamSchema::bool()),
            ("output", raster_out()),
            ("water_surface", raster_out()),
            ("flood_extent", raster_out()),
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
            (
                "method",
                ToolParamSchema::enum_values(&["nearest", "bilinear", "cubic", "lanczos"]),
            ),
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
            (
                "method",
                ToolParamSchema::enum_values(&["bilinear", "nearest", "cubic"]),
            ),
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
            (
                "method",
                ToolParamSchema::enum_values(&["bilinear", "nearest", "cubic"]),
            ),
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
            (
                "index",
                ToolParamSchema::enum_values(&["ndvi", "ndwi", "ndbi", "nbr", "evi", "savi"]),
            ),
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
            (
                "compression",
                ToolParamSchema::enum_values(&["zstd", "snappy", "gzip", "uncompressed"]),
            ),
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
            (
                "method",
                ToolParamSchema::enum_values(&[
                    "right_angles",
                    "right_angles_and_diagonals",
                    "any_angle",
                    "circle",
                ]),
            ),
            ("tolerance", float()),
            ("diagonal_penalty", float()),
            ("min_radius", float()),
            ("max_radius", float()),
        ]),
        "regularize_adjacent_building_footprint" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("group", ToolParamSchema::string()),
            (
                "method",
                ToolParamSchema::enum_values(&["right_angles", "right_angles_and_diagonals"]),
            ),
            ("tolerance", float()),
            ("precision", float()),
            ("adjacency_distance", float()),
        ]),
        "smooth_natural_features" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("segment_length", float()),
            ("iterations", int()),
            ("preserve_area", ToolParamSchema::bool()),
        ]),
        "transform_features" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("dx", float()),
            ("dy", float()),
            ("angle", float()),
            ("scale_x", float()),
            ("scale_y", float()),
            (
                "mirror_axis",
                ToolParamSchema::enum_values(&["NONE", "X", "Y"]),
            ),
            (
                "anchor",
                ToolParamSchema::enum_values(&["CENTROID", "ORIGIN", "XY"]),
            ),
            ("anchor_x", float()),
            ("anchor_y", float()),
        ]),
        "eliminate_polygons" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("max_area", float()),
            ("where", ToolParamSchema::string()),
            ("exclude", ToolParamSchema::string()),
            (
                "strategy",
                ToolParamSchema::enum_values(&["longest_border", "largest_area"]),
            ),
            ("tolerance", float()),
        ]),
        "eliminate_polygon_part" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            (
                "condition",
                ToolParamSchema::enum_values(&["AREA", "PERCENT"]),
            ),
            ("min_area", float()),
            ("percentage", float()),
            (
                "part_option",
                ToolParamSchema::enum_values(&["CONTAINED_ONLY", "ANY"]),
            ),
        ]),
        "simplify_3d_line" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("tolerance", float()),
            ("z_factor", float()),
        ]),
        "simplify_by_circular_arcs" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("tolerance", float()),
            ("mode", ToolParamSchema::enum_values(&["arcs", "tangent"])),
            ("min_arc_angle", float()),
            ("max_radius", float()),
            ("densify_output", ToolParamSchema::bool()),
            ("output_arcs", table_out()),
        ]),
        "simplify_building" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("tolerance", float()),
            ("minimum_area", float()),
            ("keep_collapsed_points", ToolParamSchema::bool()),
            ("corner_tolerance", float()),
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
            (
                "algorithm",
                ToolParamSchema::enum_values(&["paek", "bezier"]),
            ),
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
        "interpolate_from_spatiotemporal_points" => schemas(&[
            ("input", vector_in()),
            ("value_field", ToolParamSchema::string()),
            ("time_field", ToolParamSchema::string()),
            ("output", raster_out()),
            (
                "time_step",
                ToolParamSchema::enum_values(&[
                    "daily",
                    "weekly",
                    "monthly",
                    "quarterly",
                    "yearly",
                ]),
            ),
            ("cell_size", float()),
            (
                "method",
                ToolParamSchema::enum_values(&["idw", "nearest", "mean", "median"]),
            ),
            ("power", float()),
            ("neighbors", int()),
            ("min_points", int()),
        ]),
        "interpolate_shape" => schemas(&[
            ("input", vector_in()),
            ("surface", raster_in()),
            ("output", vector_out()),
            ("sample_distance", float()),
            (
                "method",
                ToolParamSchema::enum_values(&["bilinear", "nearest"]),
            ),
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
        "non_maximum_suppression" => schemas(&[
            ("input", vector_in()),
            ("confidence_score_field", ToolParamSchema::string()),
            ("output", vector_out()),
            ("max_overlap_ratio", float()),
            ("class_value_field", ToolParamSchema::string()),
        ]),
        "subdivide_polygon" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            (
                "method",
                ToolParamSchema::enum_values(&["equal_parts", "equal_areas"]),
            ),
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
            (
                "format",
                ToolParamSchema::enum_values(&["geojson", "fgb", "parquet", "shp"]),
            ),
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
            (
                "statistic",
                ToolParamSchema::enum_values(&["central_feature", "linear_directional_mean"]),
            ),
            ("weight_field", ToolParamSchema::string()),
            ("case_field", ToolParamSchema::string()),
            (
                "distance",
                ToolParamSchema::enum_values(&["euclidean", "manhattan"]),
            ),
            ("orientation_only", ToolParamSchema::bool()),
        ]),
        "block_statistics" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            (
                "statistic",
                ToolParamSchema::enum_values(&[
                    "mean", "majority", "maximum", "median", "minimum", "minority", "range", "std",
                    "sum", "variety",
                ]),
            ),
            (
                "neighborhood",
                ToolParamSchema::enum_values(&["rectangle", "circle", "annulus", "wedge"]),
            ),
            ("size", ToolParamSchema::string()),
            ("start_angle", float()),
            ("end_angle", float()),
            ("ignore_nodata", ToolParamSchema::bool()),
            ("band", int()),
        ]),
        "boundary_clean" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            (
                "method",
                ToolParamSchema::enum_values(&["majority", "expand_shrink"]),
            ),
            ("neighbors", int()),
            (
                "threshold",
                ToolParamSchema::enum_values(&["majority", "half"]),
            ),
            ("iterations", int()),
            (
                "sort",
                ToolParamSchema::enum_values(&["descending", "ascending", "none"]),
            ),
            ("band", int()),
        ]),
        "euclidean_direction" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            ("output_distance", raster_out()),
            ("output_back_direction", raster_out()),
            ("barriers", raster_in()),
            ("max_distance", float()),
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
        "find_meeting_locations" => schemas(&[
            ("input", vector_in()),
            ("track_field", ToolParamSchema::string()),
            ("time_field", ToolParamSchema::string()),
            ("search_distance", float()),
            ("min_meeting_duration", float()),
            ("max_meeting_duration", float()),
            ("min_participants", int()),
            ("time_step", float()),
            ("output", vector_out()),
            ("output_area", vector_out()),
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
        "find_dwell_locations" => schemas(&[
            ("input", vector_in()),
            ("track_field", ToolParamSchema::string()),
            ("time_field", ToolParamSchema::string()),
            ("output", vector_out()),
            ("distance_tolerance", float()),
            ("time_tolerance", float()),
            (
                "output_type",
                ToolParamSchema::enum_values(&[
                    "dwell_features",
                    "mean_centers",
                    "convex_hulls",
                    "all_features",
                ]),
            ),
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
        "near_3d" => schemas(&[
            ("input", vector_in()),
            ("near_features", vector_in()),
            ("output", vector_out()),
            ("search_radius", float()),
            ("location", ToolParamSchema::bool()),
            ("angle", ToolParamSchema::bool()),
            ("delta", ToolParamSchema::bool()),
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
            (
                "weights",
                ToolParamSchema::enum_values(&["uniform", "inverse_distance"]),
            ),
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
        "generate_near_table" => schemas(&[
            ("input", vector_in()),
            ("near_features", vector_in()),
            ("output", vector_out()),
            ("search_radius", float()),
            ("closest_count", int()),
            ("angle", ToolParamSchema::bool()),
            ("location", ToolParamSchema::bool()),
        ]),
        "aggregate_points" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("aggregation_distance", float()),
            ("min_points", int()),
            (
                "method",
                ToolParamSchema::enum_values(&["convex_hull", "buffer"]),
            ),
            ("sum_fields", ToolParamSchema::string()),
        ]),
        "forest_based_forecast" => schemas(&[
            ("input", vector_in()),
            ("location_field", ToolParamSchema::string()),
            ("time_field", ToolParamSchema::string()),
            ("value_field", ToolParamSchema::string()),
            ("output", table_out()),
            ("forecast_steps", int()),
            ("time_window", int()),
            ("validation_steps", int()),
            ("n_trees", int()),
            ("min_leaf_size", int()),
            ("max_depth", int()),
            ("seed", int()),
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
        "calculate_missing_z_values" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("placeholder", float()),
            ("method", ToolParamSchema::enum_values(&["linear", "nearest"])),
            ("extrapolate", ToolParamSchema::bool()),
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
            (
                "method",
                ToolParamSchema::enum_values(&["hilbert", "attribute"]),
            ),
            ("fields", ToolParamSchema::string()),
            ("index_field", ToolParamSchema::string()),
        ]),
        "calculate_central_meridian_and_parallels" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("field", ToolParamSchema::string()),
            ("standard_offset", float()),
        ]),
        "calculate_composite_index" => schemas(&[
            ("input", vector_in()),
            ("fields", ToolParamSchema::string()),
            ("output", vector_out()),
            (
                "scaling",
                ToolParamSchema::enum_values(&["minmax", "zscore", "percentile", "none"]),
            ),
            ("weights", ToolParamSchema::string()),
            (
                "combine",
                ToolParamSchema::enum_values(&["mean", "sum", "geometric_mean"]),
            ),
            (
                "output_range",
                ToolParamSchema::enum_values(&["minmax", "zero_to_100", "zscore", "none"]),
            ),
        ]),
        "calculate_rates" => schemas(&[
            ("input", vector_in()),
            ("count_field", ToolParamSchema::string()),
            ("population_field", ToolParamSchema::string()),
            ("output", vector_out()),
            (
                "method",
                ToolParamSchema::enum_values(&["crude", "eb_global", "eb_spatial"]),
            ),
            ("per", float()),
            ("neighbors", int()),
        ]),
        "color_polygons" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("field", ToolParamSchema::string()),
            (
                "adjacency",
                ToolParamSchema::enum_values(&["edge", "edge_or_corner"]),
            ),
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
            (
                "method",
                ToolParamSchema::enum_values(&["linear", "mann_kendall"]),
            ),
            ("intercept_output", raster_out()),
            ("significance_output", raster_out()),
            ("min_valid", int()),
            ("band", int()),
        ]),
        "warp_raster" => schemas(&[
            ("input", raster_in()),
            ("gcps", ToolParamSchema::string()),
            ("output", raster_out()),
            (
                "transform",
                ToolParamSchema::enum_values(&["poly1", "poly2", "poly3"]),
            ),
            (
                "resampling",
                ToolParamSchema::enum_values(&["nearest", "bilinear"]),
            ),
            ("cell_size", float()),
            ("epsg", int()),
            ("band", int()),
        ]),
        "weighted_voronoi" => schemas(&[
            ("input", vector_in()),
            ("output", raster_out()),
            ("weight_field", ToolParamSchema::string()),
            (
                "weight_type",
                ToolParamSchema::enum_values(&["multiplicative", "additive", "power"]),
            ),
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
            (
                "connections",
                ToolParamSchema::enum_values(&["mst", "all_neighbors"]),
            ),
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
            (
                "method",
                ToolParamSchema::enum_values(&["midpoint", "move_endpoint"]),
            ),
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
            (
                "change_type",
                ToolParamSchema::enum_values(&["mean", "slope"]),
            ),
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
            (
                "model",
                ToolParamSchema::enum_values(&["auto", "exp_smoothing", "linear", "parabolic"]),
            ),
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
        "optics_clustering" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("min_features_cluster", int()),
            ("search_distance", float()),
            ("cluster_sensitivity", float()),
        ]),
        "colocation_analysis" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("category_field", ToolParamSchema::string()),
            ("category_a", ToolParamSchema::string()),
            ("category_b", ToolParamSchema::string()),
            ("neighbors", int()),
            (
                "weight",
                ToolParamSchema::enum_values(&["gaussian", "uniform"]),
            ),
            ("permutations", int()),
            ("seed", int()),
        ]),
        "similarity_search" => schemas(&[
            ("reference", vector_in()),
            ("candidates", vector_in()),
            ("fields", ToolParamSchema::string()),
            ("output", vector_out()),
            (
                "match_method",
                ToolParamSchema::enum_values(&["euclidean", "cosine"]),
            ),
            (
                "most_or_least",
                ToolParamSchema::enum_values(&["most", "least", "both"]),
            ),
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
        "align_features" => schemas(&[
            ("input", vector_in()),
            ("target", vector_in()),
            ("output", vector_out()),
            ("search_distance", float()),
            ("match_field", ToolParamSchema::string()),
            ("target_match_field", ToolParamSchema::string()),
        ]),
        "remove_overlap_multiple" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            (
                "method",
                ToolParamSchema::enum_values(&["center_line", "thiessen"]),
            ),
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
            (
                "ring_type",
                ToolParamSchema::enum_values(&["rings", "disks"]),
            ),
            (
                "dissolve",
                ToolParamSchema::enum_values(&["none", "per_ring"]),
            ),
            ("distance_field", ToolParamSchema::string()),
        ]),
        "directional_distribution" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            (
                "statistic",
                ToolParamSchema::enum_values(&[
                    "mean_center",
                    "median_center",
                    "central_feature",
                    "standard_distance",
                    "standard_deviational_ellipse",
                ]),
            ),
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
        "summarize_within" => schemas(&[
            ("polygons", vector_in()),
            ("input", vector_in()),
            ("output", vector_out()),
            ("fields", ToolParamSchema::string()),
            ("keep_all", ToolParamSchema::bool()),
            ("shape_sum", ToolParamSchema::bool()),
            ("group_field", ToolParamSchema::string()),
        ]),
        "summarize_nearby" => schemas(&[
            ("input", vector_in()),
            ("summary_features", vector_in()),
            ("output", vector_out()),
            ("distances", ToolParamSchema::string()),
            ("sum_fields", ToolParamSchema::string()),
            ("id_field", ToolParamSchema::string()),
        ]),
        "cul_de_sac_masks" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("reference_scale", float()),
            ("symbol_width", float()),
            ("margin", float()),
            ("page_unit", ToolParamSchema::enum_values(&["points", "mm", "inches"])),
            ("map_units_per_meter", float()),
            ("tolerance", float()),
            ("attributes", ToolParamSchema::enum_values(&["ids_only", "all"])),
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
            (
                "kernel",
                ToolParamSchema::enum_values(&["gaussian", "bisquare"]),
            ),
            (
                "bandwidth_type",
                ToolParamSchema::enum_values(&["adaptive", "fixed"]),
            ),
            ("bandwidth", float()),
        ]),
        "mgwr" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("dependent_field", ToolParamSchema::string()),
            ("explanatory_fields", ToolParamSchema::string()),
            (
                "kernel",
                ToolParamSchema::enum_values(&["gaussian", "bisquare"]),
            ),
            (
                "bandwidth_type",
                ToolParamSchema::enum_values(&["adaptive", "fixed"]),
            ),
            ("tolerance", float()),
            ("max_iterations", int()),
        ]),
        "buffer_3d" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("distance", float()),
            ("distance_field", ToolParamSchema::string()),
            ("shape", ToolParamSchema::enum_values(&["round", "flat"])),
            ("quality", int()),
        ]),
        "build_balanced_zones" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("zones", int()),
            (
                "criterion",
                ToolParamSchema::enum_values(&["homogeneity", "equal_count", "equal_sum"]),
            ),
            ("fields", ToolParamSchema::string()),
            (
                "contiguity",
                ToolParamSchema::enum_values(&["rook", "queen"]),
            ),
            ("tolerance", float()),
        ]),
        "cartogram" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("value_field", ToolParamSchema::string()),
            (
                "method",
                ToolParamSchema::enum_values(&["non_contiguous", "dorling"]),
            ),
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
            (
                "aggregate",
                ToolParamSchema::enum_values(&["mean", "sum", "min", "max", "count", "median"]),
            ),
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
            (
                "statistic",
                ToolParamSchema::enum_values(&[
                    "argmax",
                    "argmin",
                    "median_position",
                    "duration",
                    "longest_run",
                ]),
            ),
            ("threshold", float()),
            (
                "comparison",
                ToolParamSchema::enum_values(&[">", ">=", "<", "<="]),
            ),
            ("dates", ToolParamSchema::string()),
            ("min_valid", int()),
        ]),
        "las_height_metrics" => schemas(&[
            ("input", lidar_in()),
            ("output", raster_out()),
            ("metrics", ToolParamSchema::string()),
            ("height_percentiles", ToolParamSchema::string()),
            ("min_height", float()),
            ("min_points", int()),
            ("cell_size", float()),
        ]),
        "unwrap_phase" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            (
                "method",
                ToolParamSchema::enum_values(&["least_squares_pcg"]),
            ),
            ("coherence", raster_in()),
            ("coherence_threshold", float()),
            ("max_iterations", int()),
            ("tolerance", float()),
            ("reference_row", int()),
            ("reference_col", int()),
            ("band", int()),
        ]),
        "flatten_interferogram" => schemas(&[
            ("input", raster_in()),
            ("dem", raster_in()),
            ("output", raster_out()),
            ("perpendicular_baseline", float()),
            ("wavelength", float()),
            // Dual-typed: a constant in degrees or a raster path.
            ("incidence_angle", ToolParamSchema::string()),
            ("slant_range", float()),
            ("reference_elevation", float()),
            ("out_topographic_phase", raster_out()),
            ("band", int()),
        ]),
        "convert_sar_units" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            (
                "conversion",
                ToolParamSchema::enum_values(&[
                    "linear_to_db",
                    "db_to_linear",
                    "amplitude_to_intensity",
                    "intensity_to_amplitude",
                    "complex_to_intensity",
                    "phase_to_displacement",
                ]),
            ),
            ("wavelength", float()),
            ("band", int()),
        ]),
        "compute_sar_indices" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            (
                "index",
                ToolParamSchema::enum_values(&["rvi", "rfdi", "csi", "dpsvi"]),
            ),
            ("polarization_bands", ToolParamSchema::string()),
            (
                "input_units",
                ToolParamSchema::enum_values(&["linear", "db"]),
            ),
        ]),
        "apply_radiometric_calibration" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            (
                "calibration_type",
                ToolParamSchema::enum_values(&["sigma0", "beta0", "gamma0"]),
            ),
            ("calibration_constant", float()),
            ("calibration_lut", raster_in()),
            // Dual-typed: a constant in degrees or a raster path.
            ("incidence_angle", ToolParamSchema::string()),
            (
                "input_units",
                ToolParamSchema::enum_values(&["dn", "amplitude", "intensity", "db"]),
            ),
            (
                "output_units",
                ToolParamSchema::enum_values(&["linear", "db"]),
            ),
            ("band", int()),
        ]),
        "cell_position_statistics" => schemas(&[
            (
                "inputs",
                ToolParamSchema::input_multiple(ToolDatasetSchema::Raster),
            ),
            ("output", raster_out()),
            (
                "statistic",
                ToolParamSchema::enum_values(&[
                    "highest_position",
                    "lowest_position",
                    "popularity",
                    "rank",
                ]),
            ),
            // Dual-typed: a raster path or a bare constant, so it stays a
            // string rather than claiming to be a raster input.
            ("selector", ToolParamSchema::string()),
            ("ignore_nodata", ToolParamSchema::bool()),
            (
                "process_as_multiband",
                ToolParamSchema::enum_values(&["single_band", "multi_band"]),
            ),
        ]),
        "frequency_comparison" => schemas(&[
            ("value_raster", raster_in()),
            (
                "inputs",
                ToolParamSchema::input_multiple(ToolDatasetSchema::Raster),
            ),
            ("output", raster_out()),
            (
                "comparison",
                ToolParamSchema::enum_values(&[
                    "equal",
                    "not_equal",
                    "greater",
                    "greater_equal",
                    "less",
                    "less_equal",
                ]),
            ),
            ("tolerance", float()),
            ("ignore_nodata", ToolParamSchema::bool()),
            (
                "process_as_multiband",
                ToolParamSchema::enum_values(&["single_band", "multi_band"]),
            ),
        ]),
        "cell_statistics" => schemas(&[
            (
                "inputs",
                ToolParamSchema::input_multiple(ToolDatasetSchema::Raster),
            ),
            ("output", raster_out()),
            (
                "statistic",
                ToolParamSchema::enum_values(&[
                    "mean",
                    "majority",
                    "maximum",
                    "median",
                    "minimum",
                    "minority",
                    "percentile",
                    "range",
                    "std",
                    "sum",
                    "variety",
                ]),
            ),
            ("ignore_nodata", ToolParamSchema::bool()),
            ("percentile_value", float()),
        ]),
        "multidimensional_anomaly" => schemas(&[
            ("input", ToolParamSchema::string()),
            ("output", raster_out()),
            (
                "method",
                ToolParamSchema::enum_values(&[
                    "difference_from_mean",
                    "percent_of_mean",
                    "z_score",
                    "difference_from_median",
                ]),
            ),
            ("reference_range", ToolParamSchema::string()),
            ("min_valid", int()),
        ]),
        "detect_incidents" => schemas(&[
            ("input", vector_in()),
            ("track_field", ToolParamSchema::string()),
            ("time_field", ToolParamSchema::string()),
            ("start_condition", ToolParamSchema::string()),
            ("end_condition", ToolParamSchema::string()),
            (
                "mode",
                ToolParamSchema::enum_values(&["points", "segments"]),
            ),
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
            (
                "mask_kind",
                ToolParamSchema::enum_values(&["exact", "convex_hull", "box"]),
            ),
            ("masked_layer", vector_in()),
            ("id_field", ToolParamSchema::string()),
        ]),
        "intersecting_layers_masks" => schemas(&[
            ("masked_layer", vector_in()),
            ("masking_layer", vector_in()),
            ("output", vector_out()),
            ("margin", float()),
            (
                "mask_kind",
                ToolParamSchema::enum_values(&["exact", "convex_hull", "box"]),
            ),
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
            (
                "orientation",
                ToolParamSchema::enum_values(&["along_line", "horizontal", "vertical"]),
            ),
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
            (
                "page_size",
                ToolParamSchema::enum_values(&[
                    "a0", "a1", "a2", "a3", "a4", "letter", "legal", "tabloid",
                ]),
            ),
            ("map_scale", float()),
            ("origin_x", float()),
            ("origin_y", float()),
            (
                "naming",
                ToolParamSchema::enum_values(&["alphanumeric", "sequential"]),
            ),
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
            (
                "input_notation",
                ToolParamSchema::enum_values(&["DD", "DMS", "DDM", "UTM", "MGRS"]),
            ),
            (
                "output_notation",
                ToolParamSchema::enum_values(&["DD", "DMS", "DDM", "UTM", "MGRS"]),
            ),
            ("coord_field", ToolParamSchema::string()),
            ("output_field", ToolParamSchema::string()),
            ("precision", int()),
            ("update_geometry", ToolParamSchema::bool()),
            ("output", vector_out()),
        ]),
        "geodetic_densify" => schemas(&[
            ("input", vector_in()),
            (
                "geodetic_type",
                ToolParamSchema::enum_values(&["geodesic", "rhumb"]),
            ),
            ("max_segment_length", float()),
            ("vertices_per_segment", int()),
            ("output", vector_out()),
        ]),
        "interpolate_with_barriers" => schemas(&[
            ("input", vector_in()),
            ("field", ToolParamSchema::string()),
            ("output", raster_out()),
            ("barriers", vector_in()),
            (
                "method",
                ToolParamSchema::enum_values(&["idw", "local_polynomial"]),
            ),
            ("power", float()),
            ("bandwidth", float()),
            ("radius", float()),
            ("cell_size", float()),
        ]),
        "generalized_linear_regression" => schemas(&[
            ("input", vector_in()),
            ("dependent_field", ToolParamSchema::string()),
            ("explanatory_fields", ToolParamSchema::string()),
            (
                "family",
                ToolParamSchema::enum_values(&["gaussian", "poisson", "logistic"]),
            ),
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
        "predict_using_trend_raster" => schemas(&[
            ("input", raster_in()),
            ("intercept", raster_in()),
            ("output", raster_out()),
            ("times", ToolParamSchema::string()),
            ("start", float()),
            ("end", float()),
            ("interval", float()),
            ("band", int()),
        ]),
        "porous_puff" => schemas(&[
            ("magnitude", raster_in()),
            ("direction", raster_in()),
            ("x", float()),
            ("y", float()),
            ("mass", float()),
            ("porosity", ToolParamSchema::string()),
            ("thickness", ToolParamSchema::string()),
            ("dispersivity_long", float()),
            ("dispersivity_trans", float()),
            ("retardation", float()),
            ("decay", float()),
            ("time", ToolParamSchema::string()),
            ("output", raster_out()),
            ("band", int()),
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
            (
                "spatial_kernel",
                ToolParamSchema::enum_values(&["epanechnikov", "quartic"]),
            ),
            (
                "temporal_kernel",
                ToolParamSchema::enum_values(&["triangular", "epanechnikov"]),
            ),
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
        "propagate_displacement" => schemas(&[
            ("input", vector_in()),
            ("links", vector_in()),
            ("output", vector_out()),
            (
                "adjustment_style",
                ToolParamSchema::enum_values(&["auto", "preserve_orientation", "solid"]),
            ),
            ("search_distance", float()),
        ]),
        "empirical_bayesian_kriging" => schemas(&[
            ("input", vector_in()),
            ("field", ToolParamSchema::string()),
            ("output", raster_out()),
            ("cell_size", float()),
            ("subset_size", int()),
            ("overlap", float()),
            ("simulations", int()),
            (
                "semivariogram",
                ToolParamSchema::enum_values(&["power", "linear", "exponential"]),
            ),
            (
                "transform",
                ToolParamSchema::enum_values(&["none", "log_empirical"]),
            ),
            ("error_output", raster_out()),
            ("seed", int()),
        ]),
        "gaussian_geostatistical_simulations" => schemas(&[
            ("input", vector_in()),
            ("value_field", ToolParamSchema::string()),
            ("output", raster_out()),
            ("num_realizations", int()),
            ("cell_size", float()),
            (
                "variogram_model",
                ToolParamSchema::enum_values(&["exponential", "spherical", "gaussian"]),
            ),
            ("nugget", float()),
            ("sill", float()),
            ("range", float()),
            ("max_neighbors", int()),
            ("seed", int()),
            ("output_mean", raster_out()),
            ("output_std", raster_out()),
        ]),
        "exploratory_regression" => schemas(&[
            ("input", vector_in()),
            ("dependent_field", ToolParamSchema::string()),
            ("explanatory_fields", ToolParamSchema::string()),
            ("min_vars", int()),
            ("max_vars", int()),
            ("max_coef_p", float()),
            ("min_adj_r2", float()),
            ("max_vif", float()),
            ("min_jb_p", float()),
            ("min_moran_p", float()),
            ("neighbors", int()),
            ("output", table_out()),
        ]),
        "causal_inference_analysis" => schemas(&[
            ("input", vector_in()),
            ("outcome_field", ToolParamSchema::string()),
            ("treatment_field", ToolParamSchema::string()),
            ("confounding_fields", ToolParamSchema::string()),
            (
                "method",
                ToolParamSchema::enum_values(&["ps_matching", "ipw", "regression_adjustment"]),
            ),
            ("add_spatial_confounders", ToolParamSchema::bool()),
            ("balance_threshold", float()),
            ("output", vector_out()),
            ("seed", int()),
        ]),
        "multivariate_clustering" => schemas(&[
            ("input", vector_in()),
            ("fields", ToolParamSchema::string()),
            ("output", vector_out()),
            ("num_clusters", int()),
            (
                "method",
                ToolParamSchema::enum_values(&["kmeans", "kmedoids"]),
            ),
            ("seed", int()),
        ]),
        "table_to_geometry" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            (
                "mode",
                ToolParamSchema::enum_values(&["xy_to_line", "bearing_distance", "ellipse"]),
            ),
            (
                "line_type",
                ToolParamSchema::enum_values(&["geodesic", "rhumb", "planar"]),
            ),
            ("vertex_spacing", float()),
            ("polygon_output", ToolParamSchema::bool()),
            ("start_x", ToolParamSchema::string()),
            ("start_y", ToolParamSchema::string()),
            ("end_x", ToolParamSchema::string()),
            ("end_y", ToolParamSchema::string()),
            ("x", ToolParamSchema::string()),
            ("y", ToolParamSchema::string()),
            ("bearing", ToolParamSchema::string()),
            ("distance", ToolParamSchema::string()),
            ("major", ToolParamSchema::string()),
            ("minor", ToolParamSchema::string()),
            ("azimuth", ToolParamSchema::string()),
        ]),
        "transform_fields" => schemas(&[
            ("input", vector_in()),
            ("fields", ToolParamSchema::string()),
            (
                "transform",
                ToolParamSchema::enum_values(&[
                    "zscore", "minmax", "robust", "log", "log1p", "sqrt", "boxcox", "inverse",
                    "bin", "onehot",
                ]),
            ),
            ("output", vector_out()),
            ("bins", int()),
            (
                "bin_method",
                ToolParamSchema::enum_values(&["equal_interval", "quantile", "std_dev"]),
            ),
            ("boxcox_lambda", float()),
            ("suffix", ToolParamSchema::string()),
            ("drop_input", ToolParamSchema::bool()),
        ]),
        "detect_graphic_conflict" => schemas(&[
            ("input", vector_in()),
            ("conflict", vector_in()),
            ("symbol_width", float()),
            ("conflict_symbol_width", float()),
            ("conflict_distance", float()),
            ("line_connection_allowance", float()),
            ("output", vector_out()),
        ]),
        "disperse_markers" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("min_spacing", float()),
            (
                "pattern",
                ToolParamSchema::enum_values(&["expanded", "ring", "cross", "square"]),
            ),
            ("seed", int()),
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
        "extract_scanned_features" => schemas(&[
            ("input", raster_in()),
            ("output", vector_out()),
            (
                "feature_type",
                ToolParamSchema::enum_values(&["lines", "polygons"]),
            ),
            ("foreground_value", float()),
            ("threshold", float()),
            ("noise_size", int()),
            ("hole_size", int()),
            ("gap_distance", float()),
            ("simplify_tolerance", float()),
            ("band", int()),
        ]),
        "gtfs_to_features" => schemas(&[
            ("input", ToolParamSchema::input(ToolDatasetSchema::File)),
            ("stops_output", vector_out()),
            ("shapes_output", vector_out()),
            ("frequency", ToolParamSchema::bool()),
            ("start_time", ToolParamSchema::string()),
            ("end_time", ToolParamSchema::string()),
        ]),
        "create_spatial_sampling_locations" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            (
                "method",
                ToolParamSchema::enum_values(&[
                    "simple_random",
                    "stratified",
                    "systematic",
                    "cluster",
                ]),
            ),
            ("num_samples", int()),
            ("strata_field", ToolParamSchema::string()),
            (
                "allocation",
                ToolParamSchema::enum_values(&["proportional", "equal", "population_field"]),
            ),
            ("population_field", ToolParamSchema::string()),
            (
                "bin_shape",
                ToolParamSchema::enum_values(&["square", "hexagon", "triangle"]),
            ),
            ("bin_size", float()),
            ("num_clusters", int()),
            ("min_distance", float()),
            ("seed", int()),
        ]),
        "compute_accuracy_for_object_detection" => schemas(&[
            ("detected", vector_in()),
            ("ground_truth", vector_in()),
            ("output", table_out()),
            ("detected_class_field", ToolParamSchema::string()),
            ("ground_truth_class_field", ToolParamSchema::string()),
            ("confidence_field", ToolParamSchema::string()),
            ("min_iou", float()),
        ]),
        "contour_with_barriers" => schemas(&[
            ("input", raster_in()),
            ("output", vector_out()),
            ("barriers", vector_in()),
            ("interval", float()),
            ("base", float()),
            ("levels", ToolParamSchema::string()),
            ("band", int()),
        ]),
        "percentile_contours" => schemas(&[
            ("input", raster_in()),
            ("output", vector_out()),
            ("percentiles", ToolParamSchema::string()),
            ("mode", ToolParamSchema::enum_values(&["value", "volume"])),
            ("ignore_negative", ToolParamSchema::bool()),
            ("smooth", ToolParamSchema::bool()),
            ("band", int()),
        ]),
        "spatial_association_between_zones" => schemas(&[
            ("zones1", raster_in()),
            ("zones2", raster_in()),
            ("output", table_out()),
            ("band1", int()),
            ("band2", int()),
        ]),
        "merge_lines_by_pseudo_node" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("dissolve_fields", ToolParamSchema::string()),
            ("snap_tolerance", float()),
        ]),
        "identify_narrow_polygons" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("width_tolerance", float()),
            ("min_narrow_area", float()),
            ("narrow_only", ToolParamSchema::bool()),
        ]),
        // note: `narrow_area`, `min_width`, `is_narrow` are output attributes, not params
        "points_to_path" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("order_field", ToolParamSchema::string()),
            ("group_field", ToolParamSchema::string()),
            ("natural_sort", ToolParamSchema::bool()),
            ("close_path", ToolParamSchema::bool()),
        ]),
        "zonal_geometry" => schemas(&[
            ("zones", raster_in()),
            ("output", raster_out()),
            (
                "measure",
                ToolParamSchema::enum_values(&[
                    "area",
                    "perimeter",
                    "thickness",
                    "centroid_x",
                    "centroid_y",
                    "major_axis",
                    "minor_axis",
                    "orientation",
                ]),
            ),
            ("as_table", ToolParamSchema::bool()),
        ]),
        "zonal_characterization" => schemas(&[
            ("zones", raster_in()),
            ("rasters", ToolParamSchema::string()),
            ("output", table_out()),
            ("percentile", float()),
            ("zone_band", int()),
            ("ignore_nodata", ToolParamSchema::bool()),
        ]),
        "zonal_fill" => schemas(&[
            ("zones", raster_in()),
            ("weight", raster_in()),
            ("output", raster_out()),
        ]),
        "calculate_polygon_main_angle" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("angle_field", ToolParamSchema::string()),
            (
                "convention",
                ToolParamSchema::enum_values(&["arithmetic", "geographic", "graphic"]),
            ),
        ]),
        "band_collection_statistics" => schemas(&[
            (
                "inputs",
                ToolParamSchema::input_multiple(ToolDatasetSchema::Raster),
            ),
            ("output", file_out()),
            (
                "detail",
                ToolParamSchema::enum_values(&["detailed", "brief"]),
            ),
        ]),
        "line_statistics" => schemas(&[
            ("input", vector_in()),
            ("field", ToolParamSchema::string()),
            ("output", raster_out()),
            (
                "statistic",
                ToolParamSchema::enum_values(&[
                    "mean", "majority", "maximum", "median", "minimum", "minority", "range",
                    "variety", "length",
                ]),
            ),
            ("search_radius", float()),
            ("cell_size", float()),
        ]),
        "optimized_hot_spot_analysis" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("analysis_field", ToolParamSchema::string()),
            (
                "aggregation",
                ToolParamSchema::enum_values(&["fishnet", "snap"]),
            ),
            ("cell_size", float()),
            ("distance_band", float()),
        ]),
        "optimized_outlier_analysis" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("analysis_field", ToolParamSchema::string()),
            (
                "aggregation",
                ToolParamSchema::enum_values(&["fishnet", "snap"]),
            ),
            ("cell_size", float()),
            ("distance_band", float()),
            ("fdr_alpha", float()),
        ]),
        "polygon_volume" => schemas(&[
            ("surface", raster_in()),
            ("input", vector_in()),
            ("height_field", ToolParamSchema::string()),
            (
                "direction",
                ToolParamSchema::enum_values(&["above", "below", "both"]),
            ),
            ("band", int()),
            ("output", vector_out()),
        ]),
        "hot_spot_analysis_comparison" => schemas(&[
            ("input1", vector_in()),
            ("input2", vector_in()),
            ("bin_field", ToolParamSchema::string()),
            ("significance", float()),
            ("match_field", ToolParamSchema::string()),
            ("tolerance", float()),
            ("permutations", int()),
            ("seed", int()),
            ("output", vector_out()),
        ]),
        "group_by_proximity" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            (
                "relationship",
                ToolParamSchema::enum_values(&["near", "intersects"]),
            ),
            ("spatial_near_distance", float()),
            ("attribute_field", ToolParamSchema::string()),
            ("group_field", ToolParamSchema::string()),
        ]),
        "feature_to_line" => schemas(&[
            ("input", vector_in()),
            ("output", vector_out()),
            ("cluster_tolerance", float()),
            ("attributes", ToolParamSchema::bool()),
        ]),
        "split_raster" => schemas(&[
            ("input", raster_in()),
            ("output_dir", file_out()),
            ("base_name", ToolParamSchema::string()),
            (
                "split_method",
                ToolParamSchema::enum_values(&["count", "size", "polygons"]),
            ),
            ("num_x", int()),
            ("num_y", int()),
            ("tile_size_x", int()),
            ("tile_size_y", int()),
            ("polygons", vector_in()),
            ("overlap", int()),
            ("format", ToolParamSchema::enum_values(&["tif", "png"])),
        ]),
        "sar_coherence" => schemas(&[
            ("reference", raster_in()),
            ("secondary", raster_in()),
            ("output", raster_out()),
            ("output_phase", raster_out()),
            ("window_size", ToolParamSchema::string()),
            ("bias_correction", ToolParamSchema::bool()),
        ]),
        "surface_parameters" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            (
                "parameter",
                ToolParamSchema::enum_values(&[
                    "slope",
                    "aspect",
                    "mean_curvature",
                    "profile_curvature",
                    "tangential_curvature",
                    "plan_curvature",
                    "gaussian_curvature",
                    "casorati_curvature",
                    "contour_geodesic_torsion",
                ]),
            ),
            (
                "fit",
                ToolParamSchema::enum_values(&["quadratic", "biquadratic"]),
            ),
            ("neighborhood_distance", float()),
            (
                "neighborhood_type",
                ToolParamSchema::enum_values(&["fixed", "adaptive"]),
            ),
            ("band", int()),
        ]),
        "rescale_by_function" => schemas(&[
            ("input", raster_in()),
            ("output", raster_out()),
            (
                "function",
                ToolParamSchema::enum_values(&[
                    "linear",
                    "inverse_linear",
                    "power",
                    "exponential",
                    "logarithmic",
                    "logistic",
                    "gaussian",
                    "near",
                    "small",
                    "large",
                    "symmetric_linear",
                ]),
            ),
            ("from_scale", float()),
            ("to_scale", float()),
            ("param1", float()),
            ("param2", float()),
            ("low_threshold", float()),
            ("high_threshold", float()),
            ("value_below", float()),
            ("value_above", float()),
            ("band", int()),
        ]),
        "calculate_transit_service_frequency" => schemas(&[
            ("input", ToolParamSchema::input(ToolDatasetSchema::File)),
            ("target", ToolParamSchema::enum_values(&["stops", "lines"])),
            ("date", ToolParamSchema::string()),
            ("start_time", ToolParamSchema::string()),
            ("duration_minutes", float()),
            (
                "count",
                ToolParamSchema::enum_values(&["departures", "arrivals"]),
            ),
            ("output", vector_out()),
        ]),
        // ── Whitebox tools with params mistyped by keyword inference ─────────────
        //
        // `wbcore::manifest_with_io_schema_json` infers a param's type from
        // keywords in its name and description. Five params across the first four
        // tools below are boolean flags whose descriptions happen to contain the
        // words "output" or "if true / writes" that trigger the inference to
        // classify them as output dataset paths instead of scalars. Providing
        // explicit schemas here short-circuits the inference for the whole tool.
        //
        // Additional corrections bundled with each tool:
        //   buffer_vector     – `mitre_limit` inferred as string, is a float
        //   individual_tree_segmentation – `veg_classes` inferred as text-file
        //                        input (is a plain string); `grid_cell_size` and
        //                        `grid_refine_iterations` inferred as raster
        //                        inputs (are scalars); `output_id_mode` inferred
        //                        as a file output (is an enum selector)
        //   lidar_tile        – `origin_x`/`origin_y` inferred as raster inputs
        //                        (are coordinate floats)
        //
        // Count-like params (`quadrant_segments`, `adaptive_neighbors`,
        // `adaptive_sector_count`, `grid_refine_iterations`, `max_iterations`,
        // `min_cluster_points`, `threads`, `seed`, `min_points_in_tile`) are
        // integers upstream — they are all read with `parse_usize_alias` — so
        // they get `int()` rather than `float()`.
        "buffer_vector" => schemas(&[
            ("input",             vector_in()),
            ("distance",          float()),
            ("quadrant_segments", int()),
            ("cap_style",         ToolParamSchema::enum_values(&["round", "flat", "square"])),
            ("join_style",        ToolParamSchema::enum_values(&["round", "bevel", "mitre"])),
            ("mitre_limit",       float()),      // was: string
            ("dissolve",          ToolParamSchema::bool()), // was: vector input
            ("output",            vector_out()),
        ]),
        "individual_tree_segmentation" => schemas(&[
            ("input",                  lidar_in()),
            ("only_use_veg",           ToolParamSchema::bool()),
            ("veg_classes",            ToolParamSchema::string()), // was: text dataset input
            ("min_height",             float()),
            ("max_height",             float()),
            ("bandwidth_min",          float()),
            ("bandwidth_max",          float()),
            ("adaptive_bandwidth",     ToolParamSchema::bool()),
            ("adaptive_neighbors",     int()),
            ("adaptive_sector_count",  int()),
            ("grid_acceleration",      ToolParamSchema::bool()),
            ("grid_cell_size",         float()),   // was: raster input
            ("grid_refine_exact",      ToolParamSchema::bool()),
            ("grid_refine_iterations", int()),     // was: raster input
            ("tile_size",              float()),
            ("tile_overlap",           float()),
            ("vertical_bandwidth",     float()),
            ("max_iterations",         int()),
            ("convergence_tol",        float()),
            ("min_cluster_points",     int()),
            ("mode_merge_dist",        float()),
            ("threads",                int()),
            ("simd",                   ToolParamSchema::bool()),
            ("output_id_mode",         ToolParamSchema::enum_values(&[ // was: file output
                "rgb", "user_data", "point_source_id",
                "rgb+user_data", "rgb+point_source_id",
            ])),
            ("output_sidecar_csv",     ToolParamSchema::bool()), // was: table output
            ("seed",                   int()),
            ("output",                 lidar_out()),
        ]),
        "las_to_shapefile" => schemas(&[
            ("input",             lidar_in()),
            ("output",            vector_out()),
            ("output_multipoint", ToolParamSchema::bool()), // was: vector output
        ]),
        "lidar_tile" => schemas(&[
            ("input",              lidar_in()),
            ("tile_width",         float()),
            ("tile_height",        float()),
            ("origin_x",           float()),   // was: raster input
            ("origin_y",           float()),   // was: raster input
            ("min_points_in_tile", int()),
            ("output_laz_format",  ToolParamSchema::bool()), // was: lidar output
            ("output_directory",   file_out()),
        ]),
        "lidar_tile_footprint" => schemas(&[
            ("input",        lidar_in()),
            ("output",       vector_out()),
            ("output_hulls", ToolParamSchema::bool()), // was: file output
        ]),
        // `download_osm_vector` reads no dataset input at all (it queries
        // Overpass for an extent), yet the inference gave it two. The
        // `split_output_by_geometry` flag became a *vector input*, so a host that
        // resolves inputs before running the tool asked for a layer to satisfy a
        // checkbox and could not run the tool at all; `provenance_output` writes a
        // sidecar but was inferred as an existing JSON input; and the EPSG codes,
        // timeout and counts all came out as floats.
        //
        // The types below are the ones the tool's own arg parsing reads
        // (`wbtools_oss/src/tools/gis/osm_download.rs`): `as_u64` for each int,
        // `as_f64` for the two floats, and `as_str` for `cache_dir`, which is a
        // directory path rather than a dataset. The two enum lists
        // (`OSM_FILTER_PRESETS`, and the `overpass_profile` validation) are
        // spelled out because an explicit table short-circuits the inference for
        // the whole tool, which had been parsing them out of the descriptions.
        "download_osm_vector" => schemas(&[
            ("west",                     float()),
            ("south",                    float()),
            ("east",                     float()),
            ("north",                    float()),
            ("input_extent_epsg",        int()),   // was: float
            ("filter_preset",            ToolParamSchema::enum_values(&[
                "all", "roads", "buildings", "water", "landuse", "trails",
                "parks", "rail", "amenities", "boundaries", "transit", "poi",
            ])),
            ("include_tags",             ToolParamSchema::string()),
            ("include_key_values",       ToolParamSchema::string()),
            ("filter_key",               ToolParamSchema::string()),
            ("filter_key_value",         ToolParamSchema::string()),
            ("include_points",           ToolParamSchema::bool()),
            ("include_lines",            ToolParamSchema::bool()),
            ("include_polygons",         ToolParamSchema::bool()),
            ("clip_to_extent",           ToolParamSchema::bool()),
            ("split_output_by_geometry", ToolParamSchema::bool()), // was: vector input
            ("output_epsg",              int()),   // was: file output
            ("overpass_profile",         ToolParamSchema::enum_values(&[
                "main", "kumi", "fr", "custom",
            ])),
            ("overpass_url",             ToolParamSchema::string()),
            ("timeout_seconds",          int()),   // was: float
            ("max_elements",             int()),   // was: float
            ("chunk_large_aoi",          ToolParamSchema::bool()),
            ("chunk_max_area_deg2",      float()),
            ("max_chunk_count",          int()),   // was: float
            ("chunk_parallel_requests",  int()),   // was: float
            ("cache_dir",                ToolParamSchema::string()), // was: json input
            ("cache_ttl_hours",          int()),   // was: float
            ("provenance_output",        json_out()), // was: json input
            ("output",                   vector_out()),
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
    fn corrected_whitebox_schemas_keep_their_exact_types() {
        // The five whitebox overrides above exist only to correct keyword
        // inference, so asserting that the param names are present is not
        // enough — a regression back to a dataset or float schema would still
        // pass. Pin the exact type of every param they declare.
        let float = ToolParamSchema::scalar_float();
        let int = ToolParamSchema::scalar_integer();
        let boolean = ToolParamSchema::bool();
        let lidar_in = ToolParamSchema::input_lidar();
        let lidar_out = ToolParamSchema::output(ToolDatasetSchema::Lidar);
        let vector_out = ToolParamSchema::output_vector_any();

        let expected: &[(&str, &[(&str, ToolParamSchema)])] = &[
            (
                "buffer_vector",
                &[
                    ("input", ToolParamSchema::input_vector_any()),
                    ("distance", float.clone()),
                    ("quadrant_segments", int.clone()),
                    (
                        "cap_style",
                        ToolParamSchema::enum_values(&["round", "flat", "square"]),
                    ),
                    (
                        "join_style",
                        ToolParamSchema::enum_values(&["round", "bevel", "mitre"]),
                    ),
                    ("mitre_limit", float.clone()),
                    ("dissolve", boolean.clone()),
                    ("output", vector_out.clone()),
                ],
            ),
            (
                "individual_tree_segmentation",
                &[
                    ("input", lidar_in.clone()),
                    ("only_use_veg", boolean.clone()),
                    ("veg_classes", ToolParamSchema::string()),
                    ("min_height", float.clone()),
                    ("max_height", float.clone()),
                    ("bandwidth_min", float.clone()),
                    ("bandwidth_max", float.clone()),
                    ("adaptive_bandwidth", boolean.clone()),
                    ("adaptive_neighbors", int.clone()),
                    ("adaptive_sector_count", int.clone()),
                    ("grid_acceleration", boolean.clone()),
                    ("grid_cell_size", float.clone()),
                    ("grid_refine_exact", boolean.clone()),
                    ("grid_refine_iterations", int.clone()),
                    ("tile_size", float.clone()),
                    ("tile_overlap", float.clone()),
                    ("vertical_bandwidth", float.clone()),
                    ("max_iterations", int.clone()),
                    ("convergence_tol", float.clone()),
                    ("min_cluster_points", int.clone()),
                    ("mode_merge_dist", float.clone()),
                    ("threads", int.clone()),
                    ("simd", boolean.clone()),
                    (
                        "output_id_mode",
                        ToolParamSchema::enum_values(&[
                            "rgb",
                            "user_data",
                            "point_source_id",
                            "rgb+user_data",
                            "rgb+point_source_id",
                        ]),
                    ),
                    ("output_sidecar_csv", boolean.clone()),
                    ("seed", int.clone()),
                    ("output", lidar_out.clone()),
                ],
            ),
            (
                "las_to_shapefile",
                &[
                    ("input", lidar_in.clone()),
                    ("output", vector_out.clone()),
                    ("output_multipoint", boolean.clone()),
                ],
            ),
            (
                "lidar_tile",
                &[
                    ("input", lidar_in.clone()),
                    ("tile_width", float.clone()),
                    ("tile_height", float.clone()),
                    ("origin_x", float.clone()),
                    ("origin_y", float.clone()),
                    ("min_points_in_tile", int.clone()),
                    ("output_laz_format", boolean.clone()),
                    (
                        "output_directory",
                        ToolParamSchema::output(ToolDatasetSchema::File),
                    ),
                ],
            ),
            (
                "lidar_tile_footprint",
                &[
                    ("input", lidar_in.clone()),
                    ("output", vector_out.clone()),
                    ("output_hulls", boolean.clone()),
                ],
            ),
            (
                "download_osm_vector",
                &[
                    ("west", float.clone()),
                    ("south", float.clone()),
                    ("east", float.clone()),
                    ("north", float.clone()),
                    ("input_extent_epsg", int.clone()),
                    (
                        "filter_preset",
                        ToolParamSchema::enum_values(&[
                            "all",
                            "roads",
                            "buildings",
                            "water",
                            "landuse",
                            "trails",
                            "parks",
                            "rail",
                            "amenities",
                            "boundaries",
                            "transit",
                            "poi",
                        ]),
                    ),
                    ("include_tags", ToolParamSchema::string()),
                    ("include_key_values", ToolParamSchema::string()),
                    ("filter_key", ToolParamSchema::string()),
                    ("filter_key_value", ToolParamSchema::string()),
                    ("include_points", boolean.clone()),
                    ("include_lines", boolean.clone()),
                    ("include_polygons", boolean.clone()),
                    ("clip_to_extent", boolean.clone()),
                    ("split_output_by_geometry", boolean.clone()),
                    ("output_epsg", int.clone()),
                    (
                        "overpass_profile",
                        ToolParamSchema::enum_values(&["main", "kumi", "fr", "custom"]),
                    ),
                    ("overpass_url", ToolParamSchema::string()),
                    ("timeout_seconds", int.clone()),
                    ("max_elements", int.clone()),
                    ("chunk_large_aoi", boolean.clone()),
                    ("chunk_max_area_deg2", float.clone()),
                    ("max_chunk_count", int.clone()),
                    ("chunk_parallel_requests", int.clone()),
                    ("cache_dir", ToolParamSchema::string()),
                    ("cache_ttl_hours", int.clone()),
                    (
                        "provenance_output",
                        ToolParamSchema::output(ToolDatasetSchema::Json),
                    ),
                    ("output", vector_out.clone()),
                ],
            ),
        ];

        for (tool_id, params) in expected {
            let schemas = geolibre_param_schemas(tool_id)
                .unwrap_or_else(|| panic!("missing param schemas for tool '{tool_id}'"));
            assert_eq!(
                schemas.len(),
                params.len(),
                "tool '{tool_id}' declares a different number of params than expected"
            );
            for (name, want) in *params {
                let got = schemas
                    .get(*name)
                    .unwrap_or_else(|| panic!("tool '{tool_id}' is missing param '{name}'"));
                assert_eq!(got, want, "tool '{tool_id}' param '{name}' has the wrong schema");
            }
        }
    }
}
