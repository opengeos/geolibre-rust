//! Sink (maximum depression extent) extraction, ported from `lidar`'s
//! `filling.py` `ExtractSinks`.
//!
//! Pipeline: fill depressions (Wang & Liu), subtract the original DEM to get the
//! fill depth, group the filled cells into components of more than `min_size`
//! pixels, then emit the sink elevations, the labeled regions, the fill depth,
//! the filled DEM, and a per-region attribute CSV.

use std::collections::HashMap;

use serde_json::{json, Map, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{
    band_to_vec, load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
    write_text_output,
};
use crate::fill::fill_depressions_wang_and_liu;
use crate::polygonize::{polygonize_to_geojson, PolygonizeParams};
use crate::regions::{region_group, region_props, RegionProps};

/// Arrays and attributes produced by sink extraction. All buffers are row-major
/// of length `rows * cols`.
pub struct SinkResult {
    /// DEM elevation inside sinks, `0` elsewhere (no-data value `0`).
    pub sink: Vec<f64>,
    /// 1-based region id per cell, `0` background.
    pub region: Vec<f64>,
    /// Fill depth (`filled - dem`) inside sinks, `0` elsewhere.
    pub depth: Vec<f64>,
    /// Depression-filled DEM (no-data preserved).
    pub filled: Vec<f64>,
    /// Per-region attributes (ascending region id).
    pub props: Vec<RegionProps>,
}

/// Core of `ExtractSinks`, independent of raster I/O so other tools (e.g.
/// `DelineateMounts`) can reuse it.
pub fn extract_sinks_core(
    dem: &[f64],
    rows: usize,
    cols: usize,
    nodata: f64,
    min_size: u64,
    flat_increment: f64,
) -> SinkResult {
    let filled = fill_depressions_wang_and_liu(dem, rows, cols, nodata, flat_increment);

    // Fill depth; 0 outside depressions and at no-data cells.
    let mut diff = vec![0.0_f64; rows * cols];
    for i in 0..rows * cols {
        if dem[i] != nodata && filled[i] != nodata {
            let d = filled[i] - dem[i];
            diff[i] = if d > 0.0 { d } else { 0.0 };
        }
    }

    // Group filled cells into depression regions larger than min_size.
    let (labels, count) = region_group(&diff, rows, cols, min_size, nodata);
    let props = region_props(&labels, count, dem, rows, cols);

    let mut sink = vec![0.0_f64; rows * cols];
    let mut region = vec![0.0_f64; rows * cols];
    let mut depth = vec![0.0_f64; rows * cols];
    for i in 0..rows * cols {
        if labels[i] != 0 && dem[i] != nodata {
            sink[i] = dem[i];
            region[i] = labels[i] as f64;
            depth[i] = diff[i];
        }
    }

    SinkResult {
        sink,
        region,
        depth,
        filled,
        props,
    }
}

/// Builds the per-region attribute CSV emitted by `ExtractSinks`.
pub fn sink_csv(props: &[RegionProps], resolution: f64) -> String {
    let mut csv = String::from(
        "region_id,count,area,volume,avg_depth,max_depth,min_elev,max_elev,perimeter,\
         major_axis,minor_axis,elongatedness,eccentricity,orientation,area_bbox_ratio\n",
    );
    let res2 = resolution * resolution;
    for p in props {
        let count = p.area as f64;
        let size = count * res2;
        let max_depth = p.max_intensity - p.min_intensity;
        let mean_depth = (p.max_intensity * count - p.sum_intensity) / count;
        let volume = mean_depth * count * res2;
        let perimeter = p.perimeter * resolution;
        let major_axis = p.major_axis_length * resolution;
        let mut minor_axis = p.minor_axis_length * resolution;
        if minor_axis == 0.0 {
            minor_axis = resolution;
        }
        let elongatedness = major_axis / minor_axis;
        // lidar divides by the literal 3.1415 (not PI); kept for output parity.
        #[allow(clippy::approx_constant)]
        let orientation = p.orientation / 3.1415 * 180.0;
        csv.push_str(&format!(
            "{},{},{:.2},{:.2},{:.2},{:.2},{:.2},{:.2},{:.2},{:.2},{:.2},{:.2},{:.2},{:.2},{:.2}\n",
            p.label,
            p.area,
            size,
            volume,
            mean_depth,
            max_depth,
            p.min_intensity,
            p.max_intensity,
            perimeter,
            major_axis,
            minor_axis,
            elongatedness,
            p.eccentricity,
            orientation,
            p.extent,
        ));
    }
    csv
}

/// Builds a per-region attribute table (keyed by region id) for joining onto
/// polygonized features. Mirrors the columns of [`sink_csv`].
#[allow(clippy::approx_constant)]
pub fn region_props_map(props: &[RegionProps], resolution: f64) -> HashMap<i64, Map<String, Value>> {
    let res2 = resolution * resolution;
    let mut map = HashMap::new();
    for p in props {
        let count = p.area as f64;
        let max_depth = p.max_intensity - p.min_intensity;
        let mean_depth = (p.max_intensity * count - p.sum_intensity) / count;
        let major_axis = p.major_axis_length * resolution;
        let mut minor_axis = p.minor_axis_length * resolution;
        if minor_axis == 0.0 {
            minor_axis = resolution;
        }
        let mut attrs = Map::new();
        attrs.insert("region_id".to_string(), json!(p.label));
        attrs.insert("count".to_string(), json!(p.area));
        attrs.insert("area".to_string(), json!(count * res2));
        attrs.insert("volume".to_string(), json!(mean_depth * count * res2));
        attrs.insert("avg_depth".to_string(), json!(mean_depth));
        attrs.insert("max_depth".to_string(), json!(max_depth));
        attrs.insert("min_elev".to_string(), json!(p.min_intensity));
        attrs.insert("max_elev".to_string(), json!(p.max_intensity));
        attrs.insert("perimeter".to_string(), json!(p.perimeter * resolution));
        attrs.insert("major_axis".to_string(), json!(major_axis));
        attrs.insert("minor_axis".to_string(), json!(minor_axis));
        attrs.insert("elongatedness".to_string(), json!(major_axis / minor_axis));
        attrs.insert("eccentricity".to_string(), json!(p.eccentricity));
        attrs.insert("orientation".to_string(), json!(p.orientation / 3.1415 * 180.0));
        attrs.insert("area_bbox_ratio".to_string(), json!(p.extent));
        map.insert(p.label as i64, attrs);
    }
    map
}

/// Extracts sinks from a DEM (`lidar` `ExtractSinks`).
pub struct ExtractSinksTool;

impl Tool for ExtractSinksTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "extract_sinks",
            display_name: "Extract Sinks",
            summary: "Extract sinks (maximum depression extent) from a DEM (lidar filling.py).",
            category: ToolCategory::Hydrology,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input DEM raster file path.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional sink raster path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_size",
                    description: "Minimum number of pixels for a sink (kept if strictly greater; default 10).",
                    required: false,
                },
                ToolParamSpec {
                    name: "region_output",
                    description: "Optional output path for the labeled-region raster.",
                    required: false,
                },
                ToolParamSpec {
                    name: "depth_output",
                    description: "Optional output path for the fill-depth raster.",
                    required: false,
                },
                ToolParamSpec {
                    name: "filled_output",
                    description: "Optional output path for the depression-filled DEM.",
                    required: false,
                },
                ToolParamSpec {
                    name: "csv_output",
                    description: "Optional output path for the per-region attribute CSV.",
                    required: false,
                },
                ToolParamSpec {
                    name: "vector_output",
                    description: "Optional output path for region polygons as GeoJSON (attributes joined).",
                    required: false,
                },
                ToolParamSpec {
                    name: "flat_increment",
                    description: "Flat increment added during filling to enforce drainage (default 0).",
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
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args
            .get("input")
            .and_then(Value::as_str)
            .ok_or_else(|| ToolError::Validation("missing required parameter 'input'".to_string()))?;
        let output = parse_optional_output(args, "output")?;
        let region_output = parse_optional_output(args, "region_output")?;
        let depth_output = parse_optional_output(args, "depth_output")?;
        let filled_output = parse_optional_output(args, "filled_output")?;
        let csv_output = parse_optional_output(args, "csv_output")?;
        let vector_output = parse_optional_output(args, "vector_output")?;
        let min_size = args.get("min_size").and_then(Value::as_u64).unwrap_or(10);
        let flat_increment = args.get("flat_increment").and_then(Value::as_f64).unwrap_or(0.0);

        let raster = load_input_raster(input)?;
        let rows = raster.rows;
        let cols = raster.cols;
        let nodata = raster.nodata;
        let resolution = raster.cell_size_x;
        let dem = band_to_vec(&raster, 0);

        ctx.progress.info("filling depressions");
        let result = extract_sinks_core(&dem, rows, cols, nodata, min_size, flat_increment);
        ctx.progress.progress(0.7);

        let mut outputs = std::collections::BTreeMap::new();

        let sink_raster = raster_like_with_data(&raster, result.sink, 0.0, raster.data_type)?;
        let sink_path = write_or_store_output(sink_raster, output)?;
        outputs.insert("output".to_string(), json!(sink_path));

        if let Some(path) = vector_output {
            let pmap = region_props_map(&result.props, resolution);
            let geojson = polygonize_to_geojson(&PolygonizeParams {
                labels: &result.region,
                rows,
                cols,
                x_min: raster.x_min,
                y_max: raster.y_max(),
                cell_size_x: raster.cell_size_x,
                cell_size_y: raster.cell_size_y,
                epsg: raster.crs.epsg,
                props_by_id: &pmap,
            });
            write_text_output(&geojson, path)?;
            outputs.insert("vector".to_string(), json!(path));
        }
        if let Some(path) = region_output {
            let r = raster_like_with_data(&raster, result.region, 0.0, DataType::I32)?;
            outputs.insert("region".to_string(), json!(write_or_store_output(r, Some(path))?));
        }
        if let Some(path) = depth_output {
            let r = raster_like_with_data(&raster, result.depth, 0.0, raster.data_type)?;
            outputs.insert("depth".to_string(), json!(write_or_store_output(r, Some(path))?));
        }
        if let Some(path) = filled_output {
            let r = raster_like_with_data(&raster, result.filled, nodata, raster.data_type)?;
            outputs.insert("filled".to_string(), json!(write_or_store_output(r, Some(path))?));
        }
        if let Some(path) = csv_output {
            write_text_output(&sink_csv(&result.props, resolution), path)?;
            outputs.insert("csv".to_string(), json!(path));
        }

        outputs.insert("regions".to_string(), json!(result.props.len()));
        ctx.progress.progress(1.0);
        Ok(ToolRunResult { outputs })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extracts_a_simple_sink() {
        // 5x5 plateau at 10 with a 3x3 pit (value 1) in the middle.
        let rows = 5;
        let cols = 5;
        let nodata = -9999.0;
        let mut dem = vec![10.0; rows * cols];
        for r in 1..4 {
            for c in 1..4 {
                dem[r * cols + c] = 1.0;
            }
        }
        let result = extract_sinks_core(&dem, rows, cols, nodata, 1, 0.0);
        // The 9-cell pit should be one region with sink elevations preserved.
        assert_eq!(result.props.len(), 1);
        assert_eq!(result.props[0].area, 9);
        assert_eq!(result.sink[2 * cols + 2], 1.0);
        assert_eq!(result.region[2 * cols + 2], 1.0);
        // Filled center is raised to the rim.
        assert_eq!(result.filled[2 * cols + 2], 10.0);
        assert!(result.depth[2 * cols + 2] > 0.0);
    }
}
