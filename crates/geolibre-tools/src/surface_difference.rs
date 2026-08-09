//! GeoLibre tool: polygons of where one surface lies above, below, or level
//! with another.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Surface Difference* (3D Analyst).
//!
//! ## Why the catalog needs it
//!
//! "Where did the glacier thin, by how much, and over what footprint?" is a
//! question about *regions*, not cells. The deliverable is a polygon layer you
//! can symbolise, join, sort by volume and hand to someone — the extent of each
//! area of change, carrying its own volume and depth statistics.
//!
//! The catalog measures the same quantities but never delineates them:
//!
//! * `cut_fill` produces the signed change raster, the scene totals, and a
//!   per-region label **raster** plus a CSV — no geometry, so a region cannot
//!   be overlaid, clipped, or joined to anything;
//! * `surface_volume` and `polygon_volume` measure against a reference plane or
//!   within polygons you already have, rather than discovering the polygons;
//! * `percentile_contours` traces one threshold at a time from a single raster
//!   and knows nothing about a second surface.
//!
//! ## Method
//!
//! The two co-registered surfaces are differenced, each cell is classified as
//! `above` / `below` / `coincident` against `tolerance`, and the resulting
//! three-class raster is traced with the shared `polygonize` machinery. Every
//! polygon carries its class, area, volume, and the minimum, maximum and mean
//! difference inside it.
//!
//! The reference may be a second raster or a constant elevation, so a surface
//! can be compared against a flat datum without first materialising one.

use std::collections::BTreeMap;
use std::collections::HashMap;

use serde_json::{json, Map, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};
use wbvector::{Feature, FieldDef, FieldType, FieldValue, GeometryType, Layer};

use crate::args_common::{band_index, f64_or, opt_f64, req_str};
use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};
use crate::geojson_geom::geometry_from_json;
use crate::polygonize::{polygonize_to_geojson, PolygonizeParams};
use crate::raster_stack::check_alignment_refs;
use crate::vector_common::write_or_store_layer;

pub struct SurfaceDifferenceTool;

impl Tool for SurfaceDifferenceTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "surface_difference",
            display_name: "Surface Difference",
            summary: "Delineates where one surface lies above, below or level with another as polygons, each carrying its class, area, volume and difference statistics (ArcGIS Surface Difference). cut_fill measures the same quantities but emits only a region-label raster and a CSV, so its regions have no geometry to overlay, clip or join; surface_volume and polygon_volume measure against a plane or within polygons you already have rather than discovering them.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Surface raster to classify.",
                    required: true,
                },
                ToolParamSpec {
                    name: "reference",
                    description: "Reference surface: a co-registered raster path, or a constant elevation to compare against a flat datum.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon layer of the classified regions. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_raster",
                    description: "Output signed difference raster (input minus reference). Always produced; stored in memory when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Differences within this of zero count as coincident (default 0.0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_area",
                    description: "Discard regions smaller than this area in map units squared.",
                    required: false,
                },
                ToolParamSpec {
                    name: "include_coincident",
                    description: "Emit polygons for the unchanged class too (default false). Off by default because on most scenes it is one enormous background polygon.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "Band of the input surface (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "reference_band",
                    description: "Band of the reference surface (default 0).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        // Dual-typed: only the constant form can be checked without the grids.
        if reference_path(args)?.is_none() {
            reference_constant(args)?;
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input_path = req_str(args, "input")?.to_string();
        let prm = parse_params(args)?;
        let band = band_index(args, "band")?;
        let ref_band = band_index(args, "reference_band")?;
        let output = parse_optional_output(args, "output")?;
        let out_raster = parse_optional_output(args, "output_raster")?;

        let surface = load_input_raster(&input_path)?;
        let (rows, cols) = (surface.rows, surface.cols);

        // Reference: a raster or a flat datum.
        let (reference, ref_label) = match reference_path(args)? {
            Some(p) => {
                let r = load_input_raster(p)?;
                check_alignment_refs(&[&surface, &r])?;
                (Some(r), "raster".to_string())
            }
            None => (None, format!("constant {}", reference_constant(args)?)),
        };
        let ref_const = if reference.is_none() {
            reference_constant(args)?
        } else {
            0.0
        };

        // Signed difference, and the class of each cell.
        let nodata = -9999.0_f64;
        let mut diff = vec![nodata; rows * cols];
        let mut labels = vec![0.0f64; rows * cols];
        let (mut n_above, mut n_below, mut n_same) = (0usize, 0usize, 0usize);

        for r in 0..rows {
            for c in 0..cols {
                let i = r * cols + c;
                let a = surface.get(band, r as isize, c as isize);
                if a == surface.nodata || !a.is_finite() {
                    continue;
                }
                let b = match &reference {
                    None => ref_const,
                    Some(rr) => {
                        let v = rr.get(ref_band, r as isize, c as isize);
                        if v == rr.nodata || !v.is_finite() {
                            continue;
                        }
                        v
                    }
                };
                let d = a - b;
                diff[i] = d;
                // Labels feed `polygonize`, which treats 0 as background, so
                // the classes start at 1.
                labels[i] = if d > prm.tolerance {
                    n_above += 1;
                    Class::Above.label_value()
                } else if d < -prm.tolerance {
                    n_below += 1;
                    Class::Below.label_value()
                } else {
                    n_same += 1;
                    Class::Coincident.label_value()
                };
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        if !prm.include_coincident {
            for v in labels.iter_mut() {
                if *v == Class::Coincident.label_value() {
                    *v = 0.0;
                }
            }
        }

        ctx.progress.info(&format!(
            "{rows}x{cols} against {ref_label}: {n_above} above, {n_below} below, {n_same} level"
        ));

        let diff_path = write_or_store_output(
            raster_like_with_data(&surface, diff.clone(), nodata, DataType::F32)?,
            out_raster,
        )?;

        // Find the connected components first and give each a unique label, so
        // `polygonize`'s feature id maps back to exactly one component. Tracing
        // the three-class raster directly and then trying to re-identify each
        // ring's component afterwards would be guesswork.
        let components = components_of(&labels, &diff, rows, cols, nodata);
        let mut unique = vec![0.0f64; rows * cols];
        for (i, comp) in components.iter().enumerate() {
            for &cell in &comp.cells {
                unique[cell] = (i + 1) as f64;
            }
        }

        let cell_area = surface.cell_size_x * surface.cell_size_y;
        let props: HashMap<i64, Map<String, Value>> = HashMap::new();
        let geojson = polygonize_to_geojson(&PolygonizeParams {
            labels: &unique,
            rows,
            cols,
            x_min: surface.x_min,
            y_max: surface.y_min + rows as f64 * surface.cell_size_y,
            cell_size_x: surface.cell_size_x,
            cell_size_y: surface.cell_size_y,
            epsg: surface.crs.epsg,
            props_by_id: &props,
        });

        let mut layer = Layer::new("surface_difference");
        layer.geom_type = Some(GeometryType::Polygon);
        if let Some(e) = surface.crs.epsg {
            layer = layer.with_crs_epsg(e);
        }
        layer.add_field(FieldDef::new("id", FieldType::Integer));
        layer.add_field(FieldDef::new("class", FieldType::Text));
        layer.add_field(FieldDef::new("cell_count", FieldType::Integer));
        layer.add_field(FieldDef::new("area", FieldType::Float));
        layer.add_field(FieldDef::new("volume", FieldType::Float));
        layer.add_field(FieldDef::new("min_difference", FieldType::Float));
        layer.add_field(FieldDef::new("max_difference", FieldType::Float));
        layer.add_field(FieldDef::new("mean_difference", FieldType::Float));

        let parsed: Value = serde_json::from_str(&geojson).map_err(|e| {
            ToolError::Execution(format!("polygonize produced invalid GeoJSON: {e}"))
        })?;
        let feats = parsed
            .get("features")
            .and_then(Value::as_array)
            .map(Vec::as_slice)
            .unwrap_or_default();

        let mut fid = 0u64;
        let mut total_volume = 0.0;
        for f in feats {
            let Some(geom) = f.get("geometry").and_then(geometry_from_json) else {
                continue;
            };
            let id = f
                .get("properties")
                .and_then(|p| p.get("id"))
                .and_then(Value::as_i64)
                .unwrap_or(0);
            if id <= 0 || id as usize > components.len() {
                continue;
            }
            let comp = &components[id as usize - 1];
            let class = comp.class;
            let st = &comp.stats;
            let area = st.count as f64 * cell_area;
            if let Some(min) = prm.min_area {
                if area < min {
                    continue;
                }
            }
            let volume = st.abs_sum * cell_area;
            total_volume += volume;

            let mut feat = Feature::with_geometry(fid, geom, layer.schema.len());
            feat.set_by_index(0, FieldValue::Integer(fid as i64));
            feat.set_by_index(1, FieldValue::Text(class.name().to_string()));
            feat.set_by_index(2, FieldValue::Integer(st.count as i64));
            feat.set_by_index(3, FieldValue::Float(area));
            feat.set_by_index(4, FieldValue::Float(volume));
            feat.set_by_index(5, FieldValue::Float(st.min));
            feat.set_by_index(6, FieldValue::Float(st.max));
            feat.set_by_index(7, FieldValue::Float(st.sum / st.count as f64));
            layer.push(feat);
            fid += 1;
        }

        let feature_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("output_raster".to_string(), json!(diff_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("above_cells".to_string(), json!(n_above));
        outputs.insert("below_cells".to_string(), json!(n_below));
        outputs.insert("coincident_cells".to_string(), json!(n_same));
        outputs.insert("above_volume".to_string(), json!(class_volume(&diff, nodata, cell_area, true)));
        outputs.insert("below_volume".to_string(), json!(class_volume(&diff, nodata, cell_area, false)));
        outputs.insert("total_polygon_volume".to_string(), json!(total_volume));
        outputs.insert("reference".to_string(), json!(ref_label));
        Ok(ToolRunResult { outputs })
    }
}

/// Cut or fill volume over the whole scene.
fn class_volume(diff: &[f64], nodata: f64, cell_area: f64, above: bool) -> f64 {
    diff.iter()
        .filter(|&&d| d != nodata && d.is_finite() && ((d > 0.0) == above) && d != 0.0)
        .map(|d| d.abs() * cell_area)
        .sum::<f64>()
}

/// Difference statistics for one connected component.
struct Stats {
    count: usize,
    sum: f64,
    abs_sum: f64,
    min: f64,
    max: f64,
}

/// One connected run of same-class cells, with its difference statistics.
struct Component {
    class: Class,
    cells: Vec<usize>,
    stats: Stats,
}

/// Finds the connected components of the class raster, using the same
/// 4-connectivity rule `polygonize` applies, so relabelling them uniquely makes
/// the traced rings correspond one-to-one.
fn components_of(
    labels: &[f64],
    diff: &[f64],
    rows: usize,
    cols: usize,
    nodata: f64,
) -> Vec<Component> {
    let mut seen = vec![false; rows * cols];
    let mut out: Vec<Component> = Vec::new();
    let mut stack: Vec<usize> = Vec::new();

    for start in 0..rows * cols {
        if seen[start] || labels[start] == 0.0 {
            continue;
        }
        let code = labels[start];
        seen[start] = true;
        stack.clear();
        stack.push(start);

        let mut cells: Vec<usize> = Vec::new();
        while let Some(i) = stack.pop() {
            cells.push(i);
            let (r, c) = (i / cols, i % cols);
            let mut push = |rr: usize, cc: usize, st: &mut Vec<usize>, sn: &mut Vec<bool>| {
                let j = rr * cols + cc;
                if !sn[j] && labels[j] == code {
                    sn[j] = true;
                    st.push(j);
                }
            };
            if r > 0 {
                push(r - 1, c, &mut stack, &mut seen);
            }
            if r + 1 < rows {
                push(r + 1, c, &mut stack, &mut seen);
            }
            if c > 0 {
                push(r, c - 1, &mut stack, &mut seen);
            }
            if c + 1 < cols {
                push(r, c + 1, &mut stack, &mut seen);
            }
        }

        let mut st = Stats {
            count: 0,
            sum: 0.0,
            abs_sum: 0.0,
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
        };
        for &i in &cells {
            let d = diff[i];
            if d == nodata || !d.is_finite() {
                continue;
            }
            st.count += 1;
            st.sum += d;
            st.abs_sum += d.abs();
            st.min = st.min.min(d);
            st.max = st.max.max(d);
        }
        // A component whose cells are all no-data in the difference raster has
        // nothing to report; dropping it keeps the labelling contiguous with
        // what actually gets traced.
        let Some(class) = Class::from_label(code) else {
            continue;
        };
        if st.count == 0 {
            continue;
        }
        out.push(Component {
            class,
            cells,
            stats: st,
        });
    }
    out
}

/// Which side of the reference a cell falls on.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Class {
    Below,
    Coincident,
    Above,
}

impl Class {
    fn label_value(self) -> f64 {
        match self {
            Class::Below => 1.0,
            Class::Coincident => 2.0,
            Class::Above => 3.0,
        }
    }

    fn from_label(v: f64) -> Option<Self> {
        match v as i64 {
            1 => Some(Class::Below),
            2 => Some(Class::Coincident),
            3 => Some(Class::Above),
            _ => None,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Class::Below => "below",
            Class::Coincident => "coincident",
            Class::Above => "above",
        }
    }
}

/// The `reference` parameter as a raster path, or `None` when it is a constant.
fn reference_path(args: &ToolArgs) -> Result<Option<&str>, ToolError> {
    let raw = args.get("reference");
    if raw.is_none() || matches!(raw, Some(Value::Null)) {
        return Err(ToolError::Validation(
            "missing required parameter 'reference'".to_string(),
        ));
    }
    Ok(raw
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty() && s.parse::<f64>().is_err()))
}

/// The `reference` parameter as a constant elevation.
fn reference_constant(args: &ToolArgs) -> Result<f64, ToolError> {
    let v = opt_f64(args, "reference")?.ok_or_else(|| {
        ToolError::Validation("missing required parameter 'reference'".to_string())
    })?;
    if !v.is_finite() {
        return Err(ToolError::Validation(
            "'reference' must be a raster path or a finite elevation".to_string(),
        ));
    }
    Ok(v)
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    tolerance: f64,
    min_area: Option<f64>,
    include_coincident: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let tolerance = f64_or(args, "tolerance", 0.0)?;
    if !tolerance.is_finite() || tolerance < 0.0 {
        return Err(ToolError::Validation(format!(
            "'tolerance' must be non-negative, got {tolerance}"
        )));
    }
    let min_area = crate::args_common::opt_positive_f64(args, "min_area")?;
    let include_coincident = crate::args_common::bool_or(args, "include_coincident", false)?;
    Ok(Params {
        tolerance,
        min_area,
        include_coincident,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vector_common::load_input_layer;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn raster_of(cols: usize, rows: usize, vals: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 10.0,
            cell_size_y: Some(10.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(32610),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, vals[row * cols + col])
                    .unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, BTreeMap<String, Value>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SurfaceDifferenceTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (layer, out.outputs)
    }

    fn field(layer: &Layer, f: &Feature, name: &str) -> FieldValue {
        f.attributes[layer.schema.field_index(name).unwrap()].clone()
    }

    /// One raised block and one excavated block become two polygons with the
    /// right class, area and volume — the deliverable `cut_fill` cannot make.
    #[test]
    fn delineates_a_mound_and_a_pit() {
        let (rows, cols) = (10, 10);
        let before = vec![100.0; rows * cols];
        let mut after = before.clone();
        // A 3x3 mound 2 m high, and a 2x2 pit 5 m deep, well apart.
        for r in 1..4 {
            for c in 1..4 {
                after[r * cols + c] = 102.0;
            }
        }
        for r in 6..8 {
            for c in 6..8 {
                after[r * cols + c] = 95.0;
            }
        }
        let (layer, outputs) = run(json!({
            "input": raster_of(cols, rows, &after),
            "reference": raster_of(cols, rows, &before)
        }));

        assert_eq!(layer.len(), 2, "expected one mound and one pit polygon");
        let mut saw_above = false;
        let mut saw_below = false;
        for f in layer.iter() {
            let FieldValue::Text(class) = field(&layer, f, "class") else {
                panic!("class must be text")
            };
            let FieldValue::Float(area) = field(&layer, f, "area") else {
                panic!()
            };
            let FieldValue::Float(volume) = field(&layer, f, "volume") else {
                panic!()
            };
            match class.as_str() {
                "above" => {
                    saw_above = true;
                    // 9 cells of 100 m^2, 2 m high.
                    assert!((area - 900.0).abs() < 1e-6, "mound area {area}");
                    assert!((volume - 1800.0).abs() < 1e-3, "mound volume {volume}");
                }
                "below" => {
                    saw_below = true;
                    // 4 cells of 100 m^2, 5 m deep.
                    assert!((area - 400.0).abs() < 1e-6, "pit area {area}");
                    assert!((volume - 2000.0).abs() < 1e-3, "pit volume {volume}");
                }
                other => panic!("unexpected class {other}"),
            }
        }
        assert!(saw_above && saw_below);
        assert!((outputs["above_volume"].as_f64().unwrap() - 1800.0).abs() < 1e-3);
        assert!((outputs["below_volume"].as_f64().unwrap() - 2000.0).abs() < 1e-3);
    }

    /// Separate regions of the same class stay separate polygons, each with its
    /// own statistics — the whole point of delineating rather than totalling.
    #[test]
    fn separate_regions_are_separate_polygons() {
        let (rows, cols) = (10, 10);
        let before = vec![50.0; rows * cols];
        let mut after = before.clone();
        // Two disconnected mounds of different heights.
        after[1 * cols + 1] = 51.0;
        for r in 6..8 {
            for c in 6..8 {
                after[r * cols + c] = 54.0;
            }
        }
        let (layer, _) = run(json!({
            "input": raster_of(cols, rows, &after),
            "reference": raster_of(cols, rows, &before)
        }));
        assert_eq!(layer.len(), 2, "two disconnected mounds -> two polygons");

        let mut volumes: Vec<f64> = layer
            .iter()
            .map(|f| match field(&layer, f, "volume") {
                FieldValue::Float(v) => v,
                _ => panic!(),
            })
            .collect();
        volumes.sort_by(f64::total_cmp);
        // 1 cell x 1 m = 100 m^3; 4 cells x 4 m = 1600 m^3.
        assert!((volumes[0] - 100.0).abs() < 1e-3, "{volumes:?}");
        assert!((volumes[1] - 1600.0).abs() < 1e-3, "{volumes:?}");
    }

    /// The tolerance decides what counts as unchanged.
    #[test]
    fn tolerance_absorbs_small_differences() {
        let (rows, cols) = (6, 6);
        let before = vec![10.0; rows * cols];
        let mut after = before.clone();
        for r in 1..3 {
            for c in 1..3 {
                after[r * cols + c] = 10.3; // a 0.3 m difference
            }
        }
        let src = raster_of(cols, rows, &after);
        let refr = raster_of(cols, rows, &before);
        assert_eq!(
            run(json!({"input": src, "reference": refr})).0.len(),
            1,
            "with no tolerance the small change is a region"
        );
        let (loose, _) = run(json!({
            "input": raster_of(cols, rows, &after),
            "reference": raster_of(cols, rows, &before),
            "tolerance": 0.5
        }));
        assert_eq!(loose.len(), 0, "a 0.5 m tolerance should absorb it");
    }

    /// A constant reference compares a surface against a flat datum without
    /// materialising one.
    #[test]
    fn constant_reference_is_a_flat_datum() {
        let (rows, cols) = (6, 6);
        let mut z = vec![5.0; rows * cols];
        for r in 2..4 {
            for c in 2..4 {
                z[r * cols + c] = 15.0;
            }
        }
        let (layer, outputs) = run(json!({
            "input": raster_of(cols, rows, &z), "reference": 10.0
        }));
        assert!(outputs["reference"].as_str().unwrap().contains("constant"));
        // Everything except the block is below 10; the block is above.
        let classes: Vec<String> = layer
            .iter()
            .map(|f| match field(&layer, f, "class") {
                FieldValue::Text(t) => t,
                _ => panic!(),
            })
            .collect();
        assert!(classes.contains(&"above".to_string()));
        assert!(classes.contains(&"below".to_string()));
    }

    /// The coincident class is suppressed by default and available on request.
    #[test]
    fn coincident_polygons_are_opt_in() {
        let (rows, cols) = (6, 6);
        let before = vec![1.0; rows * cols];
        let mut after = before.clone();
        after[0] = 2.0;
        let (default, _) = run(json!({
            "input": raster_of(cols, rows, &after),
            "reference": raster_of(cols, rows, &before)
        }));
        assert_eq!(default.len(), 1, "only the changed cell by default");

        let (with_same, _) = run(json!({
            "input": raster_of(cols, rows, &after),
            "reference": raster_of(cols, rows, &before),
            "include_coincident": true
        }));
        assert_eq!(with_same.len(), 2, "the unchanged background is added");
    }

    /// The difference raster is emitted even without a path.
    #[test]
    fn difference_raster_is_always_emitted() {
        let (rows, cols) = (3, 3);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": raster_of(cols, rows, &[7.0; 9]),
            "reference": raster_of(cols, rows, &[4.0; 9])
        }))
        .unwrap();
        let out = SurfaceDifferenceTool.run(&args, &ctx()).unwrap();
        let diff = load_input_raster(out.outputs["output_raster"].as_str().unwrap()).unwrap();
        assert_eq!(diff.get(0, 1, 1), 3.0);
    }

    /// The area filter drops small regions.
    #[test]
    fn min_area_filters_small_regions() {
        let (rows, cols) = (10, 10);
        let before = vec![0.0; rows * cols];
        let mut after = before.clone();
        after[11] = 1.0; // one cell = 100 m^2
        for r in 6..9 {
            for c in 6..9 {
                after[r * cols + c] = 1.0; // nine cells = 900 m^2
            }
        }
        let (all, _) = run(json!({
            "input": raster_of(cols, rows, &after),
            "reference": raster_of(cols, rows, &before)
        }));
        assert_eq!(all.len(), 2);
        let (big, _) = run(json!({
            "input": raster_of(cols, rows, &after),
            "reference": raster_of(cols, rows, &before),
            "min_area": 500.0
        }));
        assert_eq!(big.len(), 1, "min_area should drop the single cell");
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SurfaceDifferenceTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({"input": "a.tif"})).is_err());
        assert!(bad(json!({"input": "a.tif", "reference": -1.0, "tolerance": -1})).is_err());
        assert!(bad(json!({"input": "a.tif", "reference": "b.tif"})).is_ok());
        assert!(bad(json!({"input": "a.tif", "reference": 12.5})).is_ok());
    }
}
