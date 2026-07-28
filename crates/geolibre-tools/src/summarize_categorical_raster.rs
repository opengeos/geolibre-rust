//! GeoLibre tool: per-slice class area/count table for a categorical raster.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Summarize Categorical Raster* (Image
//! Analyst).
//!
//! The nearest bundled tool is `raster_area`, which stops short in two ways: it
//! writes each class's total **back into the raster cells** rather than emitting
//! a table, and it is strictly single-band. The point of this tool is the time
//! series — given a multi-band/multidimensional classification stack (an annual
//! land-cover product, a monthly flood extent, a `landtrendr` or
//! `analyze_changes_ccdc` output), it produces a tidy table of class areas per
//! slice, which is the input to every land-cover-change chart.
//!
//! `cross_tabulation` compares two categorical rasters against each other and
//! `zonal_characterization` requires an explicit zone layer; neither gives
//! per-slice class totals for a single stack.
//!
//! Cell area is integrated **per row** for geographic (EPSG:4326) rasters, where
//! a degree of longitude shrinks toward the poles. Treating degrees as a
//! constant area is the classic bug here and is covered by a unit test.

use std::collections::{BTreeMap, BTreeSet};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Layer};

use crate::common::load_input_raster;
use crate::vector_common::{
    geometry_contains_point, load_input_layer, parse_optional_str, write_or_store_layer,
};

/// Mean Earth radius (IUGG), used for geographic cell-area integration.
const EARTH_RADIUS_M: f64 = 6_371_008.8;

/// Tabulates class counts and areas per slice of a categorical raster.
pub struct SummarizeCategoricalRasterTool;

impl Tool for SummarizeCategoricalRasterTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "summarize_categorical_raster",
            display_name: "Summarize Categorical Raster",
            summary: "Tabulates the cell count, area and share of every class in a categorical raster, one row per slice for multi-band/multidimensional stacks, optionally restricted to an area of interest (ArcGIS Summarize Categorical Raster). The bundled raster_area writes class totals back into raster cells rather than a table and is single-band only, so it cannot produce the per-date class-area series that land-cover-change analysis needs. Cell area is integrated per row for geographic rasters.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input categorical raster (single- or multi-band).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output table path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "aoi",
                    description: "Optional polygon layer restricting the summary to an area of interest.",
                    required: false,
                },
                ToolParamSpec {
                    name: "aoi_id_field",
                    description: "Field on 'aoi' identifying each area; when given, the summary is produced per AOI feature.",
                    required: false,
                },
                ToolParamSpec {
                    name: "area_units",
                    description: "Area units: map_units (default), hectares, or square_kilometers.",
                    required: false,
                },
                ToolParamSpec {
                    name: "include_nodata",
                    description: "If true, no-data is reported as its own class row (default false).",
                    required: false,
                },
                ToolParamSpec {
                    name: "bands",
                    description: "Comma-separated 1-based band list to summarize (default: all bands).",
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
        parse_area_units(args)?;
        parse_bands(args)?;
        parse_optional_bool(args, "include_nodata")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args.get("input").and_then(Value::as_str).ok_or_else(|| {
            ToolError::Validation("missing required parameter 'input'".to_string())
        })?;
        let output = parse_optional_str(args, "output")?;
        let aoi_path = parse_optional_str(args, "aoi")?;
        let aoi_id_field = parse_optional_str(args, "aoi_id_field")?;
        let (unit_label, unit_scale) = parse_area_units(args)?;
        let include_nodata = parse_optional_bool(args, "include_nodata")?.unwrap_or(false);
        let requested_bands = parse_bands(args)?;

        let raster = load_input_raster(input)?;
        let rows = raster.rows;
        let cols = raster.cols;
        let nodata = raster.nodata;

        let bands: Vec<usize> = match requested_bands {
            Some(list) => {
                for b in &list {
                    if *b == 0 || *b > raster.bands {
                        return Err(ToolError::Validation(format!(
                            "band {b} out of range (raster has {} band(s))",
                            raster.bands
                        )));
                    }
                }
                list
            }
            None => (1..=raster.bands).collect(),
        };

        // Per-row cell area in map units (or m^2 for geographic rasters).
        let row_area = row_areas(&raster);

        // AOI mask: None means "every cell participates". Otherwise each cell
        // carries EVERY AOI feature containing it, so overlaps are preserved.
        let aoi = match aoi_path {
            Some(path) => {
                let layer = load_input_layer(path)?;
                Some(build_aoi_mask(&raster, &layer, aoi_id_field)?)
            }
            None => None,
        };

        ctx.progress.info("tabulating classes");
        // (band, aoi_id, class) -> (cell_count, area)
        let mut acc: BTreeMap<(usize, usize, i64), (u64, f64)> = BTreeMap::new();
        // (band, aoi_id) -> total area, for the percentage column.
        let mut totals: BTreeMap<(usize, usize), (u64, f64)> = BTreeMap::new();

        for (bi, band_1based) in bands.iter().enumerate() {
            let band = (*band_1based - 1) as isize;
            for row in 0..rows {
                let area = row_area[row] * unit_scale;
                for col in 0..cols {
                    // A cell may belong to SEVERAL overlapping AOI features, so
                    // it is accumulated once per containing feature. A scalar
                    // "last writer wins" mask would make each per-feature
                    // summary depend on feature order and silently drop the
                    // overlap from all but one of them.
                    let slots: &[usize] = match &aoi {
                        None => &[0_usize],
                        Some((mask, _)) => {
                            let s = &mask[row * cols + col];
                            if s.is_empty() {
                                continue;
                            }
                            s.as_slice()
                        }
                    };
                    let v = raster.get(band, row as isize, col as isize);
                    let is_nodata = v == nodata || !v.is_finite();
                    if is_nodata && !include_nodata {
                        continue;
                    }
                    // Categorical rasters carry integral class codes; rounding
                    // keeps float-backed inputs (F32 land cover) from splitting
                    // one class across several near-equal keys.
                    let class = if is_nodata {
                        i64::MIN
                    } else {
                        v.round() as i64
                    };
                    for aoi_id in slots {
                        let e = acc
                            .entry((*band_1based, *aoi_id, class))
                            .or_insert((0, 0.0));
                        e.0 += 1;
                        e.1 += area;
                        let t = totals.entry((*band_1based, *aoi_id)).or_insert((0, 0.0));
                        t.0 += 1;
                        t.1 += area;
                    }
                }
            }
            ctx.progress
                .progress((bi as f64 + 1.0) / bands.len().max(1) as f64);
        }

        let aoi_labels = aoi.as_ref().map(|(_, labels)| labels.clone());

        let mut out = Layer::new("summarize_categorical_raster");
        out.add_field(FieldDef::new("band", FieldType::Integer));
        if aoi_labels.is_some() {
            out.add_field(FieldDef::new("aoi_id", FieldType::Text));
        }
        out.add_field(FieldDef::new("class", FieldType::Integer));
        out.add_field(FieldDef::new("cell_count", FieldType::Integer));
        out.add_field(FieldDef::new(
            format!("area_{unit_label}"),
            FieldType::Float,
        ));
        out.add_field(FieldDef::new("percent", FieldType::Float));

        let mut class_set: BTreeSet<i64> = BTreeSet::new();
        for ((band, aoi_id, class), (count, area)) in &acc {
            class_set.insert(*class);
            let total_area = totals.get(&(*band, *aoi_id)).map(|t| t.1).unwrap_or(0.0);
            let percent = if total_area > 0.0 {
                area / total_area * 100.0
            } else {
                0.0
            };
            let mut fields: Vec<(String, FieldValue)> =
                vec![("band".into(), FieldValue::Integer(*band as i64))];
            if let Some(labels) = &aoi_labels {
                let label = labels.get(*aoi_id).cloned().unwrap_or_default();
                fields.push(("aoi_id".into(), FieldValue::Text(label)));
            }
            // NoData rows are flagged with the sentinel class; surface them as
            // NULL rather than a nonsense integer.
            if *class == i64::MIN {
                fields.push(("class".into(), FieldValue::Null));
            } else {
                fields.push(("class".into(), FieldValue::Integer(*class)));
            }
            fields.push(("cell_count".into(), FieldValue::Integer(*count as i64)));
            fields.push((format!("area_{unit_label}"), FieldValue::Float(*area)));
            fields.push(("percent".into(), FieldValue::Float(percent)));

            let refs: Vec<(&str, FieldValue)> = fields
                .iter()
                .map(|(k, v)| (k.as_str(), v.clone()))
                .collect();
            out.add_feature(None, &refs)
                .map_err(|e| ToolError::Execution(format!("failed writing summary row: {e}")))?;
        }

        let row_count = acc.len();
        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("row_count".to_string(), json!(row_count));
        outputs.insert("band_count".to_string(), json!(bands.len()));
        outputs.insert("class_count".to_string(), json!(class_set.len()));
        Ok(ToolRunResult { outputs })
    }
}

/// Per-row cell area. Projected rasters get a constant `|cx * cy|`; geographic
/// (EPSG:4326) rasters are integrated per row in square metres, since a cell's
/// ground area shrinks toward the poles.
fn row_areas(raster: &wbraster::Raster) -> Vec<f64> {
    let rows = raster.rows;
    let cx = raster.cell_size_x.abs();
    let cy = raster.cell_size_y.abs();

    let geographic = raster.crs.epsg == Some(4326);
    if !geographic {
        return vec![cx * cy; rows];
    }

    let y_max = raster.y_min + rows as f64 * cy;
    let d_lon = cx.to_radians();
    let r2 = EARTH_RADIUS_M * EARTH_RADIUS_M;
    (0..rows)
        .map(|row| {
            let lat_top = (y_max - row as f64 * cy).clamp(-90.0, 90.0).to_radians();
            let lat_bottom = (y_max - (row as f64 + 1.0) * cy)
                .clamp(-90.0, 90.0)
                .to_radians();
            r2 * d_lon * (lat_top.sin() - lat_bottom.sin()).abs()
        })
        .collect()
}

/// Rasterises the AOI polygons onto the raster grid, returning per cell the list
/// of AOI features containing it (empty outside every AOI) and the label list.
///
/// A list rather than a scalar because AOI features may overlap: a cell in an
/// overlap belongs to each of them, and collapsing it to one would make the
/// per-feature summaries order-dependent.
fn build_aoi_mask(
    raster: &wbraster::Raster,
    layer: &Layer,
    id_field: Option<&str>,
) -> Result<(Vec<Vec<usize>>, Vec<String>), ToolError> {
    let rows = raster.rows;
    let cols = raster.cols;
    let cx = raster.cell_size_x.abs();
    let cy = raster.cell_size_y.abs();
    let y_max = raster.y_min + rows as f64 * cy;

    let mut labels: Vec<String> = Vec::new();
    let mut mask: Vec<Vec<usize>> = vec![Vec::new(); rows * cols];

    let id_idx = match id_field {
        Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
            ToolError::Validation(format!("aoi_id_field '{f}' not found on the AOI layer"))
        })?),
        None => None,
    };

    for feature in layer.iter() {
        let Some(geom) = feature.geometry.as_ref() else {
            continue;
        };
        // Without an id field every AOI feature collapses into one mask.
        let slot = if id_field.is_some() { labels.len() } else { 0 };
        let label = match id_idx {
            Some(i) => feature
                .attributes
                .get(i)
                .map(field_to_string)
                .unwrap_or_default(),
            None => "aoi".to_string(),
        };
        if id_field.is_some() || labels.is_empty() {
            labels.push(label);
        }

        // Only scan the window this feature can possibly cover: the ray cast
        // walks every ring, so testing the whole grid per feature is the
        // dominant cost on a multi-feature AOI.
        let Some((min_x, min_y, max_x, max_y)) = geometry_bounds(geom) else {
            continue;
        };
        let col_lo =
            (((min_x - raster.x_min) / cx).floor() as isize).clamp(0, cols as isize) as usize;
        let col_hi =
            (((max_x - raster.x_min) / cx).ceil() as isize + 1).clamp(0, cols as isize) as usize;
        let row_lo = (((y_max - max_y) / cy).floor() as isize).clamp(0, rows as isize) as usize;
        let row_hi = (((y_max - min_y) / cy).ceil() as isize + 1).clamp(0, rows as isize) as usize;

        for row in row_lo..row_hi {
            // Sample at the cell centre.
            let y = y_max - (row as f64 + 0.5) * cy;
            for col in col_lo..col_hi {
                let x = raster.x_min + (col as f64 + 0.5) * cx;
                if geometry_contains_point(geom, x, y) {
                    let cell = &mut mask[row * cols + col];
                    if !cell.contains(&slot) {
                        cell.push(slot);
                    }
                }
            }
        }
    }

    if labels.is_empty() {
        return Err(ToolError::Validation(
            "AOI layer contains no features".to_string(),
        ));
    }
    Ok((mask, labels))
}

/// Axis-aligned bounds of a geometry, or `None` when it holds no coordinates.
fn geometry_bounds(geom: &wbvector::Geometry) -> Option<(f64, f64, f64, f64)> {
    let mut b = (
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );
    let mut seen = false;
    let mut visit = |cs: &[wbvector::Coord], b: &mut (f64, f64, f64, f64), seen: &mut bool| {
        for c in cs {
            b.0 = b.0.min(c.x);
            b.1 = b.1.min(c.y);
            b.2 = b.2.max(c.x);
            b.3 = b.3.max(c.y);
            *seen = true;
        }
    };
    fn walk(
        g: &wbvector::Geometry,
        b: &mut (f64, f64, f64, f64),
        seen: &mut bool,
        visit: &mut impl FnMut(&[wbvector::Coord], &mut (f64, f64, f64, f64), &mut bool),
    ) {
        match g {
            wbvector::Geometry::Point(c) => visit(std::slice::from_ref(c), b, seen),
            wbvector::Geometry::MultiPoint(cs) | wbvector::Geometry::LineString(cs) => {
                visit(cs, b, seen)
            }
            wbvector::Geometry::MultiLineString(parts) => {
                for cs in parts {
                    visit(cs, b, seen);
                }
            }
            wbvector::Geometry::Polygon { exterior, .. } => visit(&exterior.0, b, seen),
            wbvector::Geometry::MultiPolygon(parts) => {
                for (ext, _) in parts {
                    visit(&ext.0, b, seen);
                }
            }
            wbvector::Geometry::GeometryCollection(gs) => {
                for g in gs {
                    walk(g, b, seen, visit);
                }
            }
        }
    }
    walk(geom, &mut b, &mut seen, &mut visit);
    if seen {
        Some(b)
    } else {
        None
    }
}

fn field_to_string(v: &FieldValue) -> String {
    match v {
        FieldValue::Text(s) => s.clone(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) => f.to_string(),
        FieldValue::Boolean(b) => b.to_string(),
        FieldValue::Null => String::new(),
        other => format!("{other:?}"),
    }
}

fn parse_area_units(args: &ToolArgs) -> Result<(&'static str, f64), ToolError> {
    match parse_optional_str(args, "area_units")? {
        None => Ok(("map_units", 1.0)),
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "map_units" | "map" => Ok(("map_units", 1.0)),
            "hectares" | "ha" => Ok(("hectares", 1.0 / 10_000.0)),
            "square_kilometers" | "sq_km" | "km2" => Ok(("square_kilometers", 1.0 / 1_000_000.0)),
            other => Err(ToolError::Validation(format!(
                "unknown area_units '{other}' (expected map_units, hectares or square_kilometers)"
            ))),
        },
    }
}

fn parse_bands(args: &ToolArgs) -> Result<Option<Vec<usize>>, ToolError> {
    let Some(s) = parse_optional_str(args, "bands")? else {
        return Ok(None);
    };
    let mut out = Vec::new();
    for part in s.split(',') {
        let t = part.trim();
        if t.is_empty() {
            continue;
        }
        out.push(t.parse::<usize>().map_err(|_| {
            ToolError::Validation(format!("parameter 'bands' has non-integer component '{t}'"))
        })?);
    }
    if out.is_empty() {
        return Ok(None);
    }
    // Deduplicate while preserving order: "1,1" would otherwise iterate band 1
    // twice and double every count and area for that band.
    let mut seen = std::collections::BTreeSet::new();
    out.retain(|b| seen.insert(*b));
    Ok(Some(out))
}

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
    use wbvector::{Coord, Geometry, Ring};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn raster_epsg(cols: usize, rows: usize, bands: usize, data: &[f64], epsg: u32) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(epsg),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for b in 0..bands {
            for row in 0..rows {
                for col in 0..cols {
                    let idx = b * rows * cols + row * cols + col;
                    r.set(b as isize, row as isize, col as isize, data[idx])
                        .unwrap();
                }
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn raster(cols: usize, rows: usize, bands: usize, data: &[f64]) -> String {
        raster_epsg(cols, rows, bands, data, 3857)
    }

    fn run_with(extra: Value) -> Layer {
        let args: ToolArgs = serde_json::from_value(extra).unwrap();
        let out = SummarizeCategoricalRasterTool.run(&args, &ctx()).unwrap();
        load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    fn get(l: &Layer, f: &wbvector::Feature, key: &str) -> Option<FieldValue> {
        l.schema.field_index(key).map(|i| f.attributes[i].clone())
    }

    fn get_int(l: &Layer, f: &wbvector::Feature, key: &str) -> Option<i64> {
        match get(l, f, key) {
            Some(FieldValue::Integer(i)) => Some(i),
            _ => None,
        }
    }

    fn get_float(l: &Layer, f: &wbvector::Feature, key: &str) -> Option<f64> {
        match get(l, f, key) {
            Some(FieldValue::Float(x)) => Some(x),
            _ => None,
        }
    }

    /// Counts and areas per class, and percentages summing to 100.
    #[test]
    fn tabulates_class_counts_and_areas() {
        // 2x2 of unit cells: three cells of class 1, one of class 2.
        let path = raster(2, 2, 1, &[1.0, 1.0, 1.0, 2.0]);
        let layer = run_with(json!({ "input": path }));
        assert_eq!(layer.iter().count(), 2);

        let c1 = layer
            .features
            .iter()
            .find(|f| get_int(&layer, f, "class") == Some(1))
            .unwrap();
        assert_eq!(get_int(&layer, c1, "cell_count"), Some(3));
        assert_eq!(get_float(&layer, c1, "area_map_units"), Some(3.0));
        assert!((get_float(&layer, c1, "percent").unwrap() - 75.0).abs() < 1e-9);

        let total: f64 = layer
            .features
            .iter()
            .map(|f| get_float(&layer, f, "percent").unwrap())
            .sum();
        assert!((total - 100.0).abs() < 1e-9, "percentages sum to 100");
    }

    /// The reason the tool exists: one row set per slice, so a stack yields a
    /// class-area time series.
    #[test]
    fn emits_one_row_set_per_band() {
        // Band 1: all class 1. Band 2: half class 1, half class 2.
        let data = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0];
        let path = raster(2, 2, 2, &data);
        let layer = run_with(json!({ "input": path }));

        let b1: Vec<_> = layer
            .features
            .iter()
            .filter(|f| get_int(&layer, f, "band") == Some(1))
            .collect();
        let b2: Vec<_> = layer
            .features
            .iter()
            .filter(|f| get_int(&layer, f, "band") == Some(2))
            .collect();
        assert_eq!(b1.len(), 1, "band 1 has a single class");
        assert_eq!(b2.len(), 2, "band 2 has two classes");
        assert_eq!(get_int(&layer, b2[0], "cell_count"), Some(2));
    }

    /// Area units convert.
    #[test]
    fn area_units_convert() {
        // 100x100 unit cells = 10,000 map units = 1 hectare.
        let path = raster(100, 100, 1, &vec![1.0; 10_000]);
        let layer = run_with(json!({ "input": path, "area_units": "hectares" }));
        let f = layer.iter().next().unwrap();
        assert!((get_float(&layer, f, "area_hectares").unwrap() - 1.0).abs() < 1e-9);
    }

    /// No-data is excluded by default and reported as a NULL class when asked.
    #[test]
    fn nodata_handling() {
        let path = raster(2, 2, 1, &[1.0, 1.0, 1.0, -9999.0]);
        let default = run_with(json!({ "input": path.clone() }));
        assert_eq!(default.iter().count(), 1);
        assert_eq!(
            get_int(&default, default.iter().next().unwrap(), "cell_count"),
            Some(3)
        );

        let with_nd = run_with(json!({ "input": path, "include_nodata": true }));
        assert_eq!(with_nd.iter().count(), 2);
        assert!(
            with_nd
                .iter()
                .any(|f| matches!(get(&with_nd, f, "class"), Some(FieldValue::Null))),
            "no-data row is reported with a NULL class"
        );
    }

    /// Geographic rasters integrate cell area per row: a cell near the pole
    /// covers far less ground than one at the equator. Treating degrees as a
    /// constant area would make these equal.
    #[test]
    fn geographic_cell_area_shrinks_toward_the_pole() {
        // Two rows of one cell each, spanning 88..90 degrees north.
        let mut r = Raster::new(RasterConfig {
            cols: 1,
            rows: 2,
            bands: 1,
            x_min: 0.0,
            y_min: 88.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(4326),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        // Row 0 (89..90 N) is class 1; row 1 (88..89 N) is class 2.
        r.set(0, 0, 0, 1.0).unwrap();
        r.set(0, 1, 0, 2.0).unwrap();
        let id = wbraster::memory_store::put_raster(r);
        let path = wbraster::memory_store::make_raster_memory_path(&id);

        let layer = run_with(json!({ "input": path }));
        let a1 = layer
            .features
            .iter()
            .find(|f| get_int(&layer, f, "class") == Some(1))
            .and_then(|f| get_float(&layer, f, "area_map_units"))
            .unwrap();
        let a2 = layer
            .features
            .iter()
            .find(|f| get_int(&layer, f, "class") == Some(2))
            .and_then(|f| get_float(&layer, f, "area_map_units"))
            .unwrap();
        assert!(
            a1 < a2 * 0.5,
            "the polar cell (89-90N) must be far smaller than the 88-89N cell, got {a1} vs {a2}"
        );
    }

    /// Repeating a band must not double-count it.
    #[test]
    fn duplicate_band_entries_are_deduplicated() {
        let path = raster(2, 2, 1, &[1.0, 1.0, 1.0, 1.0]);
        let layer = run_with(json!({ "input": path, "bands": "1,1" }));
        assert_eq!(layer.iter().count(), 1);
        assert_eq!(
            get_int(&layer, layer.iter().next().unwrap(), "cell_count"),
            Some(4),
            "band 1 listed twice must still count 4 cells, not 8"
        );
    }

    /// An AOI restricts the summary to the cells it covers.
    #[test]
    fn aoi_restricts_the_summary() {
        // 4x1 row, classes 1,1,2,2. AOI covers only the left half.
        let path = raster(4, 1, 1, &[1.0, 1.0, 2.0, 2.0]);
        let mut aoi = Layer::new("aoi");
        aoi.add_field(FieldDef::new("name", FieldType::Text));
        aoi.add_feature(
            Some(Geometry::Polygon {
                exterior: Ring(vec![
                    Coord::xy(0.0, 0.0),
                    Coord::xy(2.0, 0.0),
                    Coord::xy(2.0, 1.0),
                    Coord::xy(0.0, 1.0),
                    Coord::xy(0.0, 0.0),
                ]),
                interiors: vec![],
            }),
            &[("name", FieldValue::Text("west".into()))],
        )
        .unwrap();
        let aid = wbvector::memory_store::put_vector(aoi);
        let aoi_path = wbvector::memory_store::make_vector_memory_path(&aid);

        let layer = run_with(json!({ "input": path, "aoi": aoi_path }));
        assert_eq!(layer.iter().count(), 1, "only class 1 falls inside the AOI");
        assert_eq!(
            get_int(&layer, layer.iter().next().unwrap(), "class"),
            Some(1)
        );
        assert_eq!(
            get_int(&layer, layer.iter().next().unwrap(), "cell_count"),
            Some(2)
        );
    }

    /// Overlapping AOI features each get a complete summary: a cell in the
    /// overlap counts for BOTH, rather than only whichever feature happened to
    /// be rasterized last.
    #[test]
    fn overlapping_aoi_features_both_count_the_shared_cells() {
        // 4x1 row, all class 1. Two AOIs overlapping on the middle two cells.
        let path = raster(4, 1, 1, &[1.0, 1.0, 1.0, 1.0]);
        let mut aoi = Layer::new("aoi");
        aoi.add_field(FieldDef::new("name", FieldType::Text));
        let mut add = |x0: f64, x1: f64, name: &str| {
            aoi.add_feature(
                Some(Geometry::Polygon {
                    exterior: Ring(vec![
                        Coord::xy(x0, 0.0),
                        Coord::xy(x1, 0.0),
                        Coord::xy(x1, 1.0),
                        Coord::xy(x0, 1.0),
                        Coord::xy(x0, 0.0),
                    ]),
                    interiors: vec![],
                }),
                &[("name", FieldValue::Text(name.into()))],
            )
            .unwrap();
        };
        add(0.0, 3.0, "west");
        add(1.0, 4.0, "east");
        let aid = wbvector::memory_store::put_vector(aoi);
        let aoi_path = wbvector::memory_store::make_vector_memory_path(&aid);

        let layer = run_with(json!({
            "input": path, "aoi": aoi_path, "aoi_id_field": "name"
        }));

        let ai = layer.schema.field_index("aoi_id").unwrap();
        let mut counts = std::collections::BTreeMap::new();
        for f in layer.iter() {
            if let FieldValue::Text(name) = &f.attributes[ai] {
                counts.insert(name.clone(), get_int(&layer, f, "cell_count").unwrap());
            }
        }
        // Each AOI spans 3 cells; the middle two are shared, so both must
        // report 3 rather than one of them losing the overlap.
        assert_eq!(counts.get("west").copied(), Some(3), "counts {counts:?}");
        assert_eq!(counts.get("east").copied(), Some(3), "counts {counts:?}");
    }

    #[test]
    fn rejects_bad_parameters() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(SummarizeCategoricalRasterTool.validate(&args).is_err());

        let path = raster(2, 2, 1, &[1.0; 4]);
        for bad in [
            json!({ "input": path.clone(), "area_units": "furlongs" }),
            json!({ "input": path.clone(), "bands": "x" }),
        ] {
            let args: ToolArgs = serde_json::from_value(bad).unwrap();
            assert!(SummarizeCategoricalRasterTool.validate(&args).is_err());
        }

        // Out-of-range band is caught at run time, when the raster is known.
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": path, "bands": "9" })).unwrap();
        assert!(SummarizeCategoricalRasterTool.run(&args, &ctx()).is_err());
    }
}
