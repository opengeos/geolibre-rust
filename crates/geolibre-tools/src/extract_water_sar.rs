//! GeoLibre tool: delineate open water from SAR backscatter.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Extract Water* (Image Analyst).
//!
//! ## Why the catalog needs it
//!
//! Calm water is a specular reflector: it bounces the radar pulse away from the
//! sensor and comes back almost black, typically 10–20 dB below any land
//! surface. That makes SAR the definitive flood-mapping instrument — it is the
//! only sensor that sees the ground through the storm cloud that caused the
//! flood, day or night.
//!
//! Neither registry can do this from SAR. `spectral_index`'s NDWI and MNDWI are
//! *reflectance* indices needing green/NIR/SWIR bands, which a radar does not
//! have; `depth_to_water` is a terrain index that models where water *should*
//! be from a DEM rather than observing where it is; and thresholding by hand
//! with `raster_calculator` fails across scenes because the water/land split
//! moves with wind, incidence angle and calibration.
//!
//! ## Method
//!
//! Backscatter is converted to decibels — where the water and land modes are
//! close to Gaussian and roughly equally spread, which is what makes a
//! histogram split work at all — and the threshold is chosen by **Otsu's
//! method**, maximising between-class variance. That is a per-scene decision
//! with no hand-tuned constant, so the same call works on a calm lake and a
//! wind-roughened estuary.
//!
//! An optional DEM supplies two physical plausibility filters that remove the
//! classic false positives:
//!
//! * **radar shadow** on steep terrain is as dark as water but cannot hold any,
//!   so cells steeper than `max_slope` are rejected;
//! * **height above the local minimum** rejects dark patches perched far above
//!   the surrounding drainage, such as smooth asphalt and dry sand.
//!
//! Surviving cells are grouped into 4-connected regions, filtered by area, and
//! traced to polygons. Both the polygons and the binary mask are written.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};
use wbvector::{Feature, FieldDef, FieldType, FieldValue, GeometryType, Layer};

use crate::args_common::{band_index, f64_or, opt_f64, opt_positive_f64, req_str, usize_or};
use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};
use crate::raster_stack::check_alignment_refs;
use crate::sar_common::{
    connected_regions, otsu_threshold, power_to_db, rasterize_mask, regions_to_geometries,
    MaskSide, SarUnits,
};
use crate::solar_radiation::slope_aspect;
use crate::vector_common::{load_input_layer, write_or_store_layer};

pub struct ExtractWaterSarTool;

impl Tool for ExtractWaterSarTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "extract_water_sar",
            display_name: "Extract Water (SAR)",
            summary: "Delineates open water and flood extent from SAR backscatter with a per-scene Otsu threshold on decibels, plus optional DEM slope and height-above-minimum filters that reject radar shadow and dry dark surfaces (ArcGIS Extract Water). Neither registry can do this from radar: spectral_index's NDWI/MNDWI need green/NIR/SWIR reflectance bands a SAR does not have, depth_to_water models where water should be from terrain rather than observing it, and a hand-set raster_calculator threshold does not transfer between scenes because the water/land split moves with wind and incidence angle.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "SAR backscatter raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output water polygon layer. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_raster",
                    description: "Output binary water mask (1 water, 0 land). Always produced; stored in memory when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "dem",
                    description: "Optional DEM, co-registered with the input, enabling the slope and height filters.",
                    required: false,
                },
                ToolParamSpec {
                    name: "threshold_db",
                    description: "Fixed threshold in dB; cells at or below it are water. Overrides the automatic Otsu threshold.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_area",
                    description: "Discard water bodies smaller than this area in map units squared.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_slope",
                    description: "With a DEM, reject water on terrain steeper than this many degrees (default 5). Standing water cannot occupy a steep slope, but radar shadow there is just as dark.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_height_above_min",
                    description: "With a DEM, reject candidate cells more than this many z units above the lowest DEM cell in their own connected region.",
                    required: false,
                },
                ToolParamSpec {
                    name: "mask_features",
                    description: "Polygon layer restricting the analysis area.",
                    required: false,
                },
                ToolParamSpec {
                    name: "mask_type",
                    description: "'land_polygon' (default; analyse outside the polygons) or 'water_polygon' (analyse inside them).",
                    required: false,
                },
                ToolParamSpec {
                    name: "histogram_bins",
                    description: "Bins used for the Otsu histogram (default 256).",
                    required: false,
                },
                ToolParamSpec {
                    name: "input_units",
                    description: "'intensity' (default), 'dn', 'amplitude', or 'db'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "Band to threshold (default 0).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input_path = req_str(args, "input")?.to_string();
        let prm = parse_params(args)?;
        let band = band_index(args, "band")?;
        let output = parse_optional_output(args, "output")?;
        let out_raster = parse_optional_output(args, "output_raster")?;

        let raster = load_input_raster(&input_path)?;
        let (rows, cols) = (raster.rows, raster.cols);

        // Decibels: the water and land modes are close to Gaussian there, which
        // is the assumption Otsu's between-class variance criterion rests on.
        // Splitting a linear-power histogram instead puts the threshold far too
        // close to zero, because the land mode has a long right tail.
        let mut db = vec![f64::NAN; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let v = raster.get(band, r as isize, c as isize);
                if v == raster.nodata || !v.is_finite() {
                    continue;
                }
                if let Some(p) = prm.input_units.to_power(v) {
                    if let Some(d) = power_to_db(p) {
                        db[r * cols + c] = d;
                    }
                }
            }
        }

        let analyse = match args.get("mask_features").and_then(Value::as_str) {
            Some(p) if !p.trim().is_empty() => {
                let layer = load_input_layer(p.trim())?;
                rasterize_mask(&raster, &layer, prm.mask_side)
            }
            _ => vec![true; rows * cols],
        };

        // Optional terrain filters.
        let dem = match args.get("dem").and_then(Value::as_str) {
            Some(p) if !p.trim().is_empty() => {
                let d = load_input_raster(p.trim())?;
                check_alignment_refs(&[&raster, &d])?;
                Some(d)
            }
            _ => None,
        };
        let (elev, slope_deg) = match &dem {
            None => (None, None),
            Some(d) => {
                let (csx, csy) = (raster.cell_size_x, raster.cell_size_y);
                let mut z = vec![f64::NAN; rows * cols];
                for r in 0..rows {
                    for c in 0..cols {
                        let v = d.get(0, r as isize, c as isize);
                        if v != d.nodata && v.is_finite() {
                            z[r * cols + c] = v;
                        }
                    }
                }
                // Horn's operator assumes one cell size; use the mean when the
                // grid is close to square, which the alignment check ensures it
                // shares with the SAR raster.
                let (s, _) = slope_aspect(&z, rows, cols, (csx + csy) / 2.0);
                let deg: Vec<f64> = s.iter().map(|v| v.to_degrees()).collect();
                (Some(z), Some(deg))
            }
        };

        // Threshold: fixed if given, else Otsu over the analysis area only —
        // including masked-out land would move the split.
        let sample: Vec<f64> = (0..rows * cols)
            .filter(|&i| analyse[i] && db[i].is_finite())
            .map(|i| db[i])
            .collect();
        if sample.is_empty() {
            return Err(ToolError::Execution(
                "no valid cells in the analysis area".to_string(),
            ));
        }
        let (threshold, source) = match prm.fixed_threshold {
            Some(t) => (t, "fixed"),
            None => match otsu_threshold(&sample, prm.bins) {
                Some(t) => (t, "otsu"),
                None => {
                    return Err(ToolError::Execution(
                        "backscatter has no spread to threshold; the scene is uniform, so no \
                         water/land split exists. Supply 'threshold_db' to force one."
                            .to_string(),
                    ))
                }
            },
        };
        ctx.progress
            .info(&format!("{source} threshold {threshold:.2} dB"));

        // Candidate water: dark, inside the analysis area, not too steep.
        let mut flag = vec![false; rows * cols];
        let mut rejected_slope = 0usize;
        for i in 0..rows * cols {
            if !analyse[i] || !db[i].is_finite() || db[i] > threshold {
                continue;
            }
            if let Some(sd) = &slope_deg {
                if sd[i].is_finite() && sd[i] > prm.max_slope {
                    rejected_slope += 1;
                    continue;
                }
            }
            flag[i] = true;
        }

        let min_cells = match prm.min_area {
            None => 1,
            Some(a) => ((a / (raster.cell_size_x * raster.cell_size_y)).ceil() as usize).max(1),
        };
        let mut regions = connected_regions(&flag, &db, rows, cols, min_cells);

        // Height-above-minimum, applied per region: a genuine water body is
        // level, so cells far above its own lowest point are not part of it.
        let mut rejected_height = 0usize;
        if let (Some(z), Some(limit)) = (&elev, prm.max_height_above_min) {
            for reg in regions.iter_mut() {
                let base = reg
                    .cells
                    .iter()
                    .filter_map(|&c| z[c].is_finite().then_some(z[c]))
                    .fold(f64::INFINITY, f64::min);
                if !base.is_finite() {
                    continue;
                }
                let before = reg.cells.len();
                let mut sum = 0.0;
                reg.cells.retain(|&c| {
                    let keep = !z[c].is_finite() || z[c] - base <= limit;
                    if keep {
                        sum += db[c];
                    }
                    keep
                });
                reg.value_sum = sum;
                rejected_height += before - reg.cells.len();
            }
            regions.retain(|r| r.cells.len() >= min_cells);
        }

        ctx.progress.info(&format!(
            "{} water region(s); rejected {rejected_slope} steep and {rejected_height} raised cell(s)",
            regions.len()
        ));

        // Mask raster reflects the regions that actually survived, so it and
        // the polygons cannot disagree.
        let nodata = -9999.0_f64;
        let mut mask = vec![0.0f64; rows * cols];
        for i in 0..rows * cols {
            if !analyse[i] || !db[i].is_finite() {
                mask[i] = nodata;
            }
        }
        for reg in &regions {
            for &c in &reg.cells {
                mask[c] = 1.0;
            }
        }
        let mask_path = write_or_store_output(
            raster_like_with_data(&raster, mask, nodata, DataType::F32)?,
            out_raster,
        )?;

        // Polygons.
        let mut layer = Layer::new("water");
        layer.geom_type = Some(GeometryType::Polygon);
        if let Some(e) = raster.crs.epsg {
            layer = layer.with_crs_epsg(e);
        }
        layer.add_field(FieldDef::new("id", FieldType::Integer));
        layer.add_field(FieldDef::new("cell_count", FieldType::Integer));
        layer.add_field(FieldDef::new("area", FieldType::Float));
        layer.add_field(FieldDef::new("mean_db", FieldType::Float));
        layer.add_field(FieldDef::new("threshold_db", FieldType::Float));

        let cell_area = raster.cell_size_x * raster.cell_size_y;
        let geoms = regions_to_geometries(&raster, &regions, rows, cols)?;
        let mut fid = 0u64;
        for (idx, geom) in geoms {
            let reg = &regions[idx];
            let mut f = Feature::with_geometry(fid, geom, layer.schema.len());
            f.set_by_index(0, FieldValue::Integer(fid as i64));
            f.set_by_index(1, FieldValue::Integer(reg.cells.len() as i64));
            f.set_by_index(2, FieldValue::Float(reg.cells.len() as f64 * cell_area));
            f.set_by_index(3, FieldValue::Float(reg.mean()));
            f.set_by_index(4, FieldValue::Float(threshold));
            layer.push(f);
            fid += 1;
        }
        let total_area = regions
            .iter()
            .map(|r| r.cells.len() as f64 * cell_area)
            .sum::<f64>();
        let count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("output_raster".to_string(), json!(mask_path));
        outputs.insert("threshold_db".to_string(), json!(threshold));
        outputs.insert("threshold_source".to_string(), json!(source));
        outputs.insert("water_body_count".to_string(), json!(count));
        outputs.insert("water_area".to_string(), json!(total_area));
        outputs.insert("rejected_steep_cells".to_string(), json!(rejected_slope));
        outputs.insert("rejected_raised_cells".to_string(), json!(rejected_height));
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    fixed_threshold: Option<f64>,
    min_area: Option<f64>,
    max_slope: f64,
    max_height_above_min: Option<f64>,
    mask_side: MaskSide,
    bins: usize,
    input_units: SarUnits,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let fixed_threshold = opt_f64(args, "threshold_db")?;
    if let Some(t) = fixed_threshold {
        if !t.is_finite() {
            return Err(ToolError::Validation(
                "'threshold_db' must be a finite number of decibels".to_string(),
            ));
        }
    }
    let min_area = opt_positive_f64(args, "min_area")?;
    let max_slope = f64_or(args, "max_slope", 5.0)?;
    if !(0.0..=90.0).contains(&max_slope) {
        return Err(ToolError::Validation(format!(
            "'max_slope' must be in [0, 90] degrees, got {max_slope}"
        )));
    }
    let max_height_above_min = match opt_f64(args, "max_height_above_min")? {
        None => None,
        Some(v) if v >= 0.0 && v.is_finite() => Some(v),
        Some(v) => {
            return Err(ToolError::Validation(format!(
                "'max_height_above_min' must be non-negative, got {v}"
            )))
        }
    };
    let mask_side = MaskSide::parse(args.get("mask_type").and_then(Value::as_str).unwrap_or(""))?;
    let bins = usize_or(args, "histogram_bins", 256)?;
    if bins < 2 {
        return Err(ToolError::Validation(
            "'histogram_bins' must be at least 2".to_string(),
        ));
    }
    let input_units = SarUnits::parse(args.get("input_units").and_then(Value::as_str).unwrap_or(""))?;

    Ok(Params {
        fixed_threshold,
        min_area,
        max_slope,
        max_height_above_min,
        mask_side,
        bins,
        input_units,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
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

    fn jitter(i: usize) -> f64 {
        let mut x = (i as u64).wrapping_mul(6364136223846793005).wrapping_add(11);
        x ^= x >> 33;
        x = x.wrapping_mul(0xff51afd7ed558ccd);
        x ^= x >> 33;
        (x % 2003) as f64 / 2003.0
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
        let out = ExtractWaterSarTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (layer, out.outputs)
    }

    /// A lake in a land scene: Otsu must find the split without being told
    /// where it is, and the mapped area must match the true lake area.
    #[test]
    fn otsu_finds_a_lake() {
        let (rows, cols) = (40, 40);
        // Land at about -8 dB (0.16 linear), water at about -20 dB (0.01).
        let mut v: Vec<f64> = (0..rows * cols)
            .map(|i| 0.16 * (0.8 + 0.4 * jitter(i)))
            .collect();
        // A 10x10 lake = 100 cells of 100 m^2 = 10 000 m^2.
        for r in 10..20 {
            for c in 10..20 {
                v[r * cols + c] = 0.01 * (0.8 + 0.4 * jitter(r * cols + c));
            }
        }
        let (layer, outputs) = run(json!({ "input": raster_of(cols, rows, &v) }));
        assert_eq!(outputs["threshold_source"].as_str().unwrap(), "otsu");
        assert_eq!(layer.len(), 1, "expected exactly one lake");
        let area = outputs["water_area"].as_f64().unwrap();
        assert!(
            (area - 10_000.0).abs() < 1.0,
            "mapped area {area} should be the true 10000 m^2"
        );
        let t = outputs["threshold_db"].as_f64().unwrap();
        assert!(
            (-20.0..-8.0).contains(&t),
            "threshold {t} did not land between the water and land modes"
        );
    }

    /// The mask raster and the polygons must describe the same water — a caller
    /// that uses one and then the other must not see a different answer.
    #[test]
    fn mask_raster_agrees_with_polygons() {
        let (rows, cols) = (30, 30);
        let mut v = vec![0.16; rows * cols];
        for r in 5..12 {
            for c in 6..14 {
                v[r * cols + c] = 0.01;
            }
        }
        let args: ToolArgs =
            serde_json::from_value(json!({"input": raster_of(cols, rows, &v)})).unwrap();
        let out = ExtractWaterSarTool.run(&args, &ctx()).unwrap();
        let mask = load_input_raster(out.outputs["output_raster"].as_str().unwrap()).unwrap();
        let flagged = (0..rows)
            .flat_map(|r| (0..cols).map(move |c| (r, c)))
            .filter(|&(r, c)| mask.get(0, r as isize, c as isize) == 1.0)
            .count();
        assert_eq!(flagged, 7 * 8, "mask should flag exactly the 56 water cells");
        let area = out.outputs["water_area"].as_f64().unwrap();
        assert!((area - flagged as f64 * 100.0).abs() < 1e-6);
    }

    /// Radar shadow on a steep slope is as dark as water; the DEM filter is
    /// what tells them apart, and this is the tool's main false-positive guard.
    #[test]
    fn dem_slope_filter_rejects_radar_shadow() {
        let (rows, cols) = (30, 40);
        let mut v = vec![0.16; rows * cols];
        let mut z = vec![100.0; rows * cols];
        // A genuine lake on flat ground.
        for r in 5..12 {
            for c in 4..12 {
                v[r * cols + c] = 0.01;
            }
        }
        // A dark patch on a very steep hillside — radar shadow, not water.
        for r in 18..25 {
            for c in 24..32 {
                v[r * cols + c] = 0.01;
            }
        }
        for r in 0..rows {
            for c in 20..cols {
                z[r * cols + c] = 100.0 + 30.0 * (c - 20) as f64; // ~72 degrees
            }
        }
        let src = raster_of(cols, rows, &v);
        let dem = raster_of(cols, rows, &z);

        let (no_dem, _) = run(json!({ "input": src.clone() }));
        assert_eq!(no_dem.len(), 2, "without a DEM both dark patches map as water");

        let (with_dem, outputs) = run(json!({
            "input": src, "dem": dem, "max_slope": 5.0
        }));
        assert_eq!(
            with_dem.len(),
            1,
            "the slope filter should leave only the real lake"
        );
        assert!(outputs["rejected_steep_cells"].as_u64().unwrap() > 0);
    }

    /// A fixed threshold overrides Otsu.
    #[test]
    fn fixed_threshold_overrides_otsu() {
        let (rows, cols) = (20, 20);
        let mut v = vec![0.16; rows * cols];
        for r in 4..10 {
            for c in 4..10 {
                v[r * cols + c] = 0.01;
            }
        }
        // -30 dB is below even the water mode, so nothing qualifies.
        let (layer, outputs) = run(json!({
            "input": raster_of(cols, rows, &v), "threshold_db": -30.0
        }));
        assert_eq!(outputs["threshold_source"].as_str().unwrap(), "fixed");
        assert_eq!(layer.len(), 0, "an impossible threshold must find no water");
    }

    /// The area filter drops small ponds.
    #[test]
    fn min_area_filters_small_bodies() {
        let (rows, cols) = (30, 30);
        let mut v = vec![0.16; rows * cols];
        // Big lake: 6x6 = 36 cells = 3600 m^2.
        for r in 4..10 {
            for c in 4..10 {
                v[r * cols + c] = 0.01;
            }
        }
        // Tiny pond: 2x2 = 4 cells = 400 m^2.
        for r in 20..22 {
            for c in 20..22 {
                v[r * cols + c] = 0.01;
            }
        }
        let src = raster_of(cols, rows, &v);
        assert_eq!(run(json!({"input": src.clone()})).0.len(), 2);
        let (filtered, _) = run(json!({"input": src, "min_area": 1000.0}));
        assert_eq!(filtered.len(), 1, "min_area should drop the 400 m^2 pond");
    }

    /// A uniform scene has no water/land split; say so rather than inventing
    /// a threshold.
    #[test]
    fn uniform_scene_is_an_error() {
        let args: ToolArgs =
            serde_json::from_value(json!({"input": raster_of(8, 8, &[0.16; 64])})).unwrap();
        let err = ExtractWaterSarTool.run(&args, &ctx()).unwrap_err();
        assert!(
            format!("{err:?}").contains("uniform"),
            "expected a uniform-scene error, got {err:?}"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ExtractWaterSarTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({"input": "a.tif", "max_slope": 120})).is_err());
        assert!(bad(json!({"input": "a.tif", "min_area": -5})).is_err());
        assert!(bad(json!({"input": "a.tif", "max_height_above_min": -1})).is_err());
        assert!(bad(json!({"input": "a.tif", "histogram_bins": 1})).is_err());
        assert!(bad(json!({"input": "a.tif", "mask_type": "cloud"})).is_err());
        assert!(bad(json!({"input": "a.tif", "input_units": "watts"})).is_err());
        assert!(bad(json!({"input": "a.tif", "threshold_db": -15})).is_ok());
    }
}
