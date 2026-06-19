//! Nested-depression delineation via the level-set method, ported from
//! `lidar`'s `slicing.py` (`levelSet`, `updateLevel`, `obj_to_level`,
//! `DelineateDepressions`).
//!
//! Starting from a sink raster (elevations inside sinks, `0` elsewhere), each
//! connected sink region is sliced from its highest elevation downward at a
//! fixed `interval`. Every time a depression splits into two or more sub-basins
//! a new, deeper-level depression is recorded, building a nested hierarchy. The
//! tool emits a unique-id raster, a depression-level raster, and an attribute
//! CSV.

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
use crate::polygonize::{polygonize_to_geojson, PolygonizeParams};
use crate::regions::{region_group, region_props, RegionProps};

/// A single (possibly nested) depression.
#[derive(Debug, Clone)]
pub struct Depression {
    pub id: i64,
    pub level: i64,
    pub count: u64,
    pub size: f64,
    pub volume: f64,
    pub mean_depth: f64,
    pub max_depth: f64,
    pub min_elev: f64,
    pub max_elev: f64,
    pub in_nbr_id: Vec<i64>,
    pub region_id: u32,
    pub perimeter: f64,
    pub major_axis: f64,
    pub minor_axis: f64,
    pub elongatedness: f64,
    pub eccentricity: f64,
    pub orientation: f64,
    pub area_bbox_ratio: f64,
}

/// Builds a [`Depression`] from a region's properties at one slice level.
// lidar divides orientation by the literal 3.1415 (not PI); kept for parity.
#[allow(clippy::approx_constant)]
fn make_dep(
    obj: &RegionProps,
    unique_id: i64,
    level: i64,
    region_id: u32,
    resolution: f64,
) -> Depression {
    let count = obj.area as f64;
    let res2 = resolution * resolution;
    let max_depth = obj.max_intensity - obj.min_intensity;
    let mean_depth = (obj.max_intensity * count - obj.sum_intensity) / count;
    let mut minor_axis = obj.minor_axis_length * resolution;
    if minor_axis == 0.0 {
        minor_axis = resolution;
    }
    let major_axis = obj.major_axis_length * resolution;
    Depression {
        id: unique_id,
        level,
        count: obj.area,
        size: count * res2,
        volume: mean_depth * count * res2,
        mean_depth,
        max_depth,
        min_elev: obj.min_intensity,
        max_elev: obj.max_intensity,
        in_nbr_id: Vec::new(),
        region_id,
        perimeter: obj.perimeter * resolution,
        major_axis,
        minor_axis,
        elongatedness: major_axis / minor_axis,
        eccentricity: obj.eccentricity,
        orientation: obj.orientation / 3.1415 * 180.0,
        area_bbox_ratio: obj.extent,
    }
}

/// Writes `unique_id` into `level_img` at every cell of region `label`.
fn write_label(level_img: &mut [f64], labels: &[u32], label: u32, unique_id: i64) {
    for (i, &l) in labels.iter().enumerate() {
        if l == label {
            level_img[i] = unique_id as f64;
        }
    }
}

/// Level-set identification of nested depressions within one region.
///
/// `img` is the region's elevation subimage (row-major, `0` outside the region
/// mask). `interval` is negative (top-down). Returns the region-local object-id
/// image and the depressions found, with `id` values starting at `obj_uid + 1`.
#[allow(clippy::too_many_arguments)]
pub fn level_set(
    img: &[f64],
    rows: usize,
    cols: usize,
    region_id: u32,
    obj_uid: i64,
    min_size: u64,
    interval: f64,
    resolution: f64,
) -> (Vec<f64>, Vec<Depression>) {
    let n = rows * cols;
    let mut level_img = vec![0.0_f64; n];

    // Background (0) becomes a large sentinel; elevations are positive.
    let max_elev = img
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let big = sentinel_nodata(max_elev);
    let mut work: Vec<f64> = img.iter().map(|&v| if v == 0.0 { big } else { v }).collect();
    let min_elev = work.iter().cloned().fold(f64::INFINITY, f64::min);

    let mut dep_list: Vec<Depression> = Vec::new();
    let mut unique_id = obj_uid;

    // Degenerate 1-pixel-wide region: one depression over the elevation cells.
    // (The Python special-case here computes garbage stats over the bbox; we
    // compute sensible stats over the actual region instead.)
    if rows == 1 || cols == 1 {
        let mask: Vec<bool> = img.iter().map(|&v| v != 0.0).collect();
        let labels: Vec<u32> = mask.iter().map(|&m| if m { 1 } else { 0 }).collect();
        let props = region_props(&labels, 1, img, rows, cols);
        if let Some(obj) = props.first() {
            unique_id += 1;
            dep_list.push(make_dep(obj, unique_id, 1, region_id, resolution));
            write_label(&mut level_img, &labels, 1, unique_id);
        }
        return (level_img, dep_list);
    }

    let mut parent_ids: HashMap<i64, u32> = HashMap::new();
    let mut nbr_ids: HashMap<i64, Vec<usize>> = HashMap::new();

    let mut elev = max_elev;
    while elev > min_elev {
        for v in work.iter_mut() {
            if *v > elev {
                *v = 0.0;
            }
        }
        let (labels, count) = region_group(&work, rows, cols, min_size, big);
        if count == 0 {
            break;
        }
        let objects = region_props(&labels, count, &work, rows, cols);

        for (i, obj) in objects.iter().enumerate() {
            let (row, col) = obj.first_coord;
            if parent_ids.is_empty() {
                // First (maximum-extent) depression of the region.
                unique_id += 1;
                dep_list.push(make_dep(obj, unique_id, 1, region_id, resolution));
                parent_ids.insert(unique_id, 0);
                nbr_ids.insert(unique_id, Vec::new());
                write_label(&mut level_img, &labels, obj.label, unique_id);
            } else {
                let parent_id = level_img[row * cols + col] as i64;
                if let Some(cnt) = parent_ids.get_mut(&parent_id) {
                    *cnt += 1;
                    nbr_ids.entry(parent_id).or_default().push(i);
                }
                // parent_id not in parent_ids (already split & popped, or
                // background) is skipped, matching the non-crashing assumption.
            }
        }

        // Promote any parent that split into two or more children this slice.
        let keys: Vec<i64> = parent_ids.keys().cloned().collect();
        for key in keys {
            let split = parent_ids.get(&key).copied().unwrap_or(0) > 1;
            if split {
                let children = nbr_ids.get(&key).cloned().unwrap_or_default();
                for new_key in children {
                    let obj = &objects[new_key];
                    unique_id += 1;
                    dep_list.push(make_dep(obj, unique_id, 1, region_id, resolution));
                    let parent_idx = (key - 1 - obj_uid) as usize;
                    if let Some(parent) = dep_list.get_mut(parent_idx) {
                        parent.in_nbr_id.push(unique_id);
                    }
                    parent_ids.insert(unique_id, 0);
                    nbr_ids.insert(unique_id, Vec::new());
                    write_label(&mut level_img, &labels, obj.label, unique_id);
                }
                parent_ids.remove(&key);
            } else {
                parent_ids.insert(key, 0);
                nbr_ids.insert(key, Vec::new());
            }
        }

        elev += interval; // interval is negative
    }

    update_level(&mut dep_list, obj_uid);
    (level_img, dep_list)
}

/// Assigns hierarchy levels: a leaf is level 1; a parent is one more than its
/// deepest child. Mirrors `slicing.py` `updateLevel`.
fn update_level(dep_list: &mut [Depression], obj_uid: i64) {
    let levels_by_index = |dep_list: &[Depression], id: i64| -> i64 {
        let idx = (id - 1 - obj_uid) as usize;
        dep_list.get(idx).map(|d| d.level).unwrap_or(0)
    };
    for i in (0..dep_list.len()).rev() {
        if dep_list[i].in_nbr_id.is_empty() {
            dep_list[i].level = 1;
        } else {
            let mut max_child = 0;
            for &id in &dep_list[i].in_nbr_id.clone() {
                let l = levels_by_index(dep_list, id);
                if l > max_child {
                    max_child = l;
                }
            }
            dep_list[i].level = max_child + 1;
        }
    }
}

/// `lidar`'s `get_min_max_nodata` sentinel: `10^(floor(log10(max))+2) - 1`.
fn sentinel_nodata(max_elev: f64) -> f64 {
    if max_elev <= 0.0 {
        return 1.0e9;
    }
    let exp = max_elev.log10().floor() as i32 + 2;
    10.0_f64.powi(exp) - 1.0
}

/// Converts a region-local object-id image to a depression-level image using
/// the global depression list (`slicing.py` `obj_to_level`).
fn obj_to_level(obj_img: &[f64], global: &[Depression]) -> Vec<f64> {
    obj_img
        .iter()
        .map(|&v| {
            let id = v as i64;
            if id > 0 {
                global
                    .get((id - 1) as usize)
                    .map(|d| d.level as f64)
                    .unwrap_or(0.0)
            } else {
                0.0
            }
        })
        .collect()
}

/// Result of delineation: full-image id and level rasters plus the depression
/// list. Buffers are row-major of length `rows * cols`.
pub struct DelineateResult {
    pub obj_image: Vec<f64>,
    pub level_image: Vec<f64>,
    pub depressions: Vec<Depression>,
}

/// Core of `DelineateDepressions`, independent of raster I/O.
///
/// `sink` is the sink raster (positive elevations inside sinks, `0` elsewhere).
/// `interval` is the positive slicing interval (negated internally).
pub fn delineate_core(
    sink: &[f64],
    rows: usize,
    cols: usize,
    resolution: f64,
    min_size: u64,
    interval: f64,
) -> DelineateResult {
    let neg_interval = -interval.abs();

    // Elevation image with 0 background (matches lidar's mutated working array).
    let elev_img: Vec<f64> = sink.iter().map(|&v| if v > 0.0 { v } else { 0.0 }).collect();
    let (labels, count) = region_group(&elev_img, rows, cols, min_size, 0.0);
    let regions = region_props(&labels, count, &elev_img, rows, cols);

    let mut obj_image = vec![0.0_f64; rows * cols];
    let mut level_image = vec![0.0_f64; rows * cols];
    let mut global: Vec<Depression> = Vec::new();
    let mut obj_uid: i64 = 0;

    for region in &regions {
        let (min_row, min_col, max_row, max_col) = region.bbox;
        let br = max_row - min_row;
        let bc = max_col - min_col;

        // Region-local elevation subimage (0 outside the region mask).
        let mut img = vec![0.0_f64; br * bc];
        for r in 0..br {
            for c in 0..bc {
                let gi = (min_row + r) * cols + (min_col + c);
                if labels[gi] == region.label {
                    img[r * bc + c] = elev_img[gi];
                }
            }
        }

        let (out_obj, dep_list) = level_set(
            &img,
            br,
            bc,
            region.label,
            obj_uid,
            min_size,
            neg_interval,
            resolution,
        );

        obj_uid += dep_list.len() as i64;
        global.extend(dep_list);
        let level_obj = obj_to_level(&out_obj, &global);

        // Composite region buffers back into the full image.
        for r in 0..br {
            for c in 0..bc {
                let gi = (min_row + r) * cols + (min_col + c);
                let oid = out_obj[r * bc + c];
                if oid > 0.0 {
                    obj_image[gi] = oid;
                }
                let lvl = level_obj[r * bc + c];
                if lvl > 0.0 {
                    level_image[gi] = lvl;
                }
            }
        }
    }

    DelineateResult {
        obj_image,
        level_image,
        depressions: global,
    }
}

/// Builds the depression attribute CSV emitted by `DelineateDepressions`.
pub fn depression_csv(deps: &[Depression]) -> String {
    let mut csv = String::from(
        "id,level,count,area,volume,avg_depth,max_depth,min_elev,max_elev,children_id,\
         region_id,perimeter,major_axis,minor_axis,elongatedness,eccentricity,orientation,\
         area_bbox_ratio\n",
    );
    for d in deps {
        let children = if d.in_nbr_id.is_empty() {
            "[]".to_string()
        } else {
            let parts: Vec<String> = d.in_nbr_id.iter().map(|id| id.to_string()).collect();
            format!("[{}]", parts.join(":"))
        };
        csv.push_str(&format!(
            "{},{},{},{:.2},{:.2},{:.2},{:.2},{:.2},{:.2},{},{},{:.2},{:.2},{:.2},{:.2},{:.2},{:.2},{:.2}\n",
            d.id,
            d.level,
            d.count,
            d.size,
            d.volume,
            d.mean_depth,
            d.max_depth,
            d.min_elev,
            d.max_elev,
            children,
            d.region_id,
            d.perimeter,
            d.major_axis,
            d.minor_axis,
            d.elongatedness,
            d.eccentricity,
            d.orientation,
            d.area_bbox_ratio,
        ));
    }
    csv
}

/// Builds a depression attribute table keyed by depression id, for joining onto
/// polygonized features (mirrors the [`depression_csv`] columns).
pub fn depression_props_map(deps: &[Depression]) -> HashMap<i64, Map<String, Value>> {
    let mut map = HashMap::new();
    for d in deps {
        let children: Vec<i64> = d.in_nbr_id.clone();
        let mut attrs = Map::new();
        attrs.insert("id".to_string(), json!(d.id));
        attrs.insert("level".to_string(), json!(d.level));
        attrs.insert("count".to_string(), json!(d.count));
        attrs.insert("area".to_string(), json!(d.size));
        attrs.insert("volume".to_string(), json!(d.volume));
        attrs.insert("avg_depth".to_string(), json!(d.mean_depth));
        attrs.insert("max_depth".to_string(), json!(d.max_depth));
        attrs.insert("min_elev".to_string(), json!(d.min_elev));
        attrs.insert("max_elev".to_string(), json!(d.max_elev));
        attrs.insert("children_id".to_string(), json!(children));
        attrs.insert("region_id".to_string(), json!(d.region_id));
        attrs.insert("perimeter".to_string(), json!(d.perimeter));
        attrs.insert("major_axis".to_string(), json!(d.major_axis));
        attrs.insert("minor_axis".to_string(), json!(d.minor_axis));
        attrs.insert("elongatedness".to_string(), json!(d.elongatedness));
        attrs.insert("eccentricity".to_string(), json!(d.eccentricity));
        attrs.insert("orientation".to_string(), json!(d.orientation));
        attrs.insert("area_bbox_ratio".to_string(), json!(d.area_bbox_ratio));
        map.insert(d.id, attrs);
    }
    map
}

/// Delineates nested depressions from a sink raster (`lidar` `DelineateDepressions`).
pub struct DelineateDepressionsTool;

impl Tool for DelineateDepressionsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "delineate_depressions",
            display_name: "Delineate Depressions",
            summary: "Delineate nested depressions from a sink raster using the level-set method (lidar slicing.py).",
            category: ToolCategory::Hydrology,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input sink raster (e.g. the output of extract_sinks).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional depression-id raster path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "level_output",
                    description: "Optional output path for the depression-level raster.",
                    required: false,
                },
                ToolParamSpec {
                    name: "csv_output",
                    description: "Optional output path for the depression attribute CSV.",
                    required: false,
                },
                ToolParamSpec {
                    name: "vector_output",
                    description: "Optional output path for depression polygons as GeoJSON (attributes joined).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_size",
                    description: "Minimum number of pixels for a depression (default 10).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_depth",
                    description: "Minimum depression depth (accepted for compatibility; unused, matching lidar).",
                    required: false,
                },
                ToolParamSpec {
                    name: "interval",
                    description: "Slicing interval for the level-set method (default 0.3).",
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
        let level_output = parse_optional_output(args, "level_output")?;
        let csv_output = parse_optional_output(args, "csv_output")?;
        let vector_output = parse_optional_output(args, "vector_output")?;
        let min_size = args.get("min_size").and_then(Value::as_u64).unwrap_or(10);
        let interval = args.get("interval").and_then(Value::as_f64).unwrap_or(0.3);

        let raster = load_input_raster(input)?;
        let rows = raster.rows;
        let cols = raster.cols;
        let resolution = raster.cell_size_x;
        let sink = band_to_vec(&raster, 0);

        ctx.progress.info("slicing depressions (level-set)");
        let result = delineate_core(&sink, rows, cols, resolution, min_size, interval);
        ctx.progress.progress(0.8);

        let mut outputs = std::collections::BTreeMap::new();

        if let Some(path) = vector_output {
            let pmap = depression_props_map(&result.depressions);
            let geojson = polygonize_to_geojson(&PolygonizeParams {
                labels: &result.obj_image,
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

        let obj_raster = raster_like_with_data(&raster, result.obj_image, 0.0, DataType::I32)?;
        let obj_path = write_or_store_output(obj_raster, output)?;
        outputs.insert("output".to_string(), json!(obj_path));

        if let Some(path) = level_output {
            let r = raster_like_with_data(&raster, result.level_image, 0.0, DataType::I32)?;
            outputs.insert("level".to_string(), json!(write_or_store_output(r, Some(path))?));
        }
        if let Some(path) = csv_output {
            write_text_output(&depression_csv(&result.depressions), path)?;
            outputs.insert("csv".to_string(), json!(path));
        }

        outputs.insert("depressions".to_string(), json!(result.depressions.len()));
        ctx.progress.progress(1.0);
        Ok(ToolRunResult { outputs })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_basin_one_depression() {
        // A 6x6 sink that is one simple bowl: a single depression, level 1.
        let rows = 6;
        let cols = 6;
        let mut sink = vec![0.0_f64; rows * cols];
        for r in 1..5 {
            for c in 1..5 {
                // Bowl: deeper in the middle.
                let d = ((r as f64 - 2.5).abs()).max((c as f64 - 2.5).abs());
                sink[r * cols + c] = 5.0 + d; // 5..7
            }
        }
        let result = delineate_core(&sink, rows, cols, 1.0, 1, 0.3);
        assert!(!result.depressions.is_empty());
        // The maximum-extent depression covers the whole basin.
        assert_eq!(result.depressions[0].region_id, 1);
        assert!(result.depressions[0].count >= 1);
    }

    #[test]
    fn two_basins_split_into_nested_hierarchy() {
        // Two pits joined by a saddle: top-level depression splits into two,
        // giving a parent at level 2 with two level-1 children.
        let rows = 5;
        let cols = 9;
        let mut sink = vec![0.0_f64; rows * cols];
        // Fill a connected sink region; two deep wells at columns 2 and 6.
        for r in 1..4 {
            for c in 1..8 {
                sink[r * cols + c] = 10.0;
            }
        }
        for r in 1..4 {
            sink[r * cols + 2] = 6.0; // well A
            sink[r * cols + 6] = 6.0; // well B
        }
        let result = delineate_core(&sink, rows, cols, 1.0, 1, 1.0);
        let max_level = result.depressions.iter().map(|d| d.level).max().unwrap_or(0);
        assert!(max_level >= 2, "expected a nested hierarchy, got level {max_level}");
    }
}
