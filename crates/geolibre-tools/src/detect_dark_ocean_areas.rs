//! GeoLibre tool: find anomalously dark patches on water in SAR imagery.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Detect Dark Ocean Areas* (Image
//! Analyst).
//!
//! ## Why the catalog needs it
//!
//! A film of oil damps the short capillary waves that scatter radar back to the
//! sensor, so a slick appears as a dark patch 3–10 dB below the surrounding
//! sea. This is the operational basis of satellite oil-spill monitoring, and
//! the same signature marks natural seeps, algal films and low-wind zones.
//!
//! It is the exact complement of `detect_bright_ocean_objects` — the same
//! scene, the anomaly in the other direction — and neither registry has either.
//! It is also **not** the same problem as `extract_water_sar`: that one splits
//! a bimodal water/land histogram, whereas here everything in view is already
//! water and the slick is a modest departure from a single unimodal
//! distribution. Running an Otsu split on an all-water scene simply bisects the
//! sea clutter and reports half the ocean as a slick, which is why this tool
//! defaults to a statistical departure from the background instead.
//!
//! ## Method
//!
//! The background sea level and spread are estimated over the analysis area in
//! decibels, and a cell is flagged when it sits `threshold` standard deviations
//! *below* the mean — or, optionally, below an Otsu split or a fixed dB value.
//! Flagged cells are grouped into 4-connected regions and filtered by area,
//! since a real slick is large and contiguous while speckle is not.
//!
//! Each region is reported with its **contrast**: how many dB below the
//! background it sits. That is the number an analyst uses to separate a genuine
//! mineral-oil slick (typically more than 5 dB) from a look-alike such as a
//! wind-shadow or a rain cell.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;
use wbvector::{Feature, FieldDef, FieldType, FieldValue, GeometryType, Layer};

use crate::args_common::{band_index, choice_or, f64_or, opt_f64, opt_positive_f64, req_str, usize_or};
use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};
use crate::sar_common::{
    connected_regions, otsu_threshold, power_to_db, rasterize_mask, regions_to_geometries,
    MaskSide, SarUnits,
};
use crate::vector_common::{load_input_layer, write_or_store_layer};

pub struct DetectDarkOceanAreasTool;

impl Tool for DetectDarkOceanAreasTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "detect_dark_ocean_areas",
            display_name: "Detect Dark Ocean Areas",
            summary: "Finds anomalously dark patches on water in SAR imagery — oil slicks, natural seeps, algal films, low-wind zones — as regions sitting a set number of standard deviations below the local sea background, reported with their dB contrast (ArcGIS Detect Dark Ocean Areas). The complement of detect_bright_ocean_objects, and distinct from extract_water_sar: an all-water scene is unimodal, so an Otsu split would simply bisect the sea clutter and report half the ocean as a slick.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "SAR backscatter raster over water.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon layer of dark areas. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_raster",
                    description: "Output binary mask (1 dark, 0 background). Always produced; stored in memory when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'statistical' (default; a set number of standard deviations below the background mean), 'otsu', or 'fixed'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "threshold",
                    description: "For 'statistical', how many standard deviations below the mean counts as dark (default 2.0). For 'fixed', the dB level at or below which a cell is dark.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_area",
                    description: "Discard dark areas smaller than this area in map units squared. Slicks are large and contiguous; speckle is not.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_contrast_db",
                    description: "Discard regions sitting less than this many dB below the background. A mineral-oil slick is typically more than 5 dB down.",
                    required: false,
                },
                ToolParamSpec {
                    name: "mask_features",
                    description: "Polygon layer restricting the analysis to water. Land is bright and would inflate the background estimate.",
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
                    description: "1-based band to search (default 1).",
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

        let sample: Vec<f64> = (0..rows * cols)
            .filter(|&i| analyse[i] && db[i].is_finite())
            .map(|i| db[i])
            .collect();
        if sample.is_empty() {
            return Err(ToolError::Execution(
                "no valid cells in the analysis area".to_string(),
            ));
        }
        let n = sample.len() as f64;
        let background_mean = sample.iter().sum::<f64>() / n;
        let background_sd =
            (sample.iter().map(|v| (v - background_mean).powi(2)).sum::<f64>() / n).sqrt();

        let threshold = match prm.method {
            Method::Statistical => {
                // An exact `<= 0.0` test is not enough: summing squared
                // deviations of identical values still leaves a few ulps of
                // rounding (about 3.6e-15 dB on a uniform scene), so the guard
                // would pass and the threshold would land inside the noise —
                // splitting the scene arbitrarily. Compare against a relative
                // epsilon instead.
                if background_sd <= 1e-9 * background_mean.abs().max(1.0) {
                    return Err(ToolError::Execution(
                        "backscatter is perfectly uniform, so nothing can be anomalously dark. \
                         Use method 'fixed' to force a level."
                            .to_string(),
                    ));
                }
                background_mean - prm.threshold * background_sd
            }
            Method::Otsu => otsu_threshold(&sample, prm.bins).ok_or_else(|| {
                ToolError::Execution(
                    "backscatter has no spread to threshold; the scene is uniform.".to_string(),
                )
            })?,
            Method::Fixed => prm.threshold,
        };
        ctx.progress.info(&format!(
            "background {background_mean:.2} +/- {background_sd:.2} dB, {} threshold {threshold:.2} dB",
            prm.method.label()
        ));

        let mut flag = vec![false; rows * cols];
        for i in 0..rows * cols {
            if analyse[i] && db[i].is_finite() && db[i] <= threshold {
                flag[i] = true;
            }
        }

        let cell_area = raster.cell_size_x * raster.cell_size_y;
        let min_cells = match prm.min_area {
            None => 1,
            Some(a) => ((a / cell_area).ceil() as usize).max(1),
        };
        let mut regions = connected_regions(&flag, &db, rows, cols, min_cells);

        // Contrast filter: a large dark region that is only marginally below
        // the background is a wind shadow, not a slick.
        if let Some(min_contrast) = prm.min_contrast_db {
            regions.retain(|r| background_mean - r.mean() >= min_contrast);
        }
        ctx.progress
            .info(&format!("{} dark area(s)", regions.len()));

        // The mask reflects the regions that survived every filter, so it
        // cannot disagree with the polygons.
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

        let mut layer = Layer::new("dark_ocean_areas");
        layer.geom_type = Some(GeometryType::Polygon);
        if let Some(e) = raster.crs.epsg {
            layer = layer.with_crs_epsg(e);
        }
        layer.add_field(FieldDef::new("id", FieldType::Integer));
        layer.add_field(FieldDef::new("cell_count", FieldType::Integer));
        layer.add_field(FieldDef::new("area", FieldType::Float));
        layer.add_field(FieldDef::new("mean_db", FieldType::Float));
        layer.add_field(FieldDef::new("contrast_db", FieldType::Float));
        layer.add_field(FieldDef::new("background_db", FieldType::Float));

        let geoms = regions_to_geometries(&raster, &regions, rows, cols)?;
        for (fid, (idx, geom)) in geoms.into_iter().enumerate() {
            let fid = fid as u64;
            let reg = &regions[idx];
            let mut f = Feature::with_geometry(fid, geom, layer.schema.len());
            f.set_by_index(0, FieldValue::Integer(fid as i64));
            f.set_by_index(1, FieldValue::Integer(reg.cells.len() as i64));
            f.set_by_index(2, FieldValue::Float(reg.cells.len() as f64 * cell_area));
            f.set_by_index(3, FieldValue::Float(reg.mean()));
            f.set_by_index(4, FieldValue::Float(background_mean - reg.mean()));
            f.set_by_index(5, FieldValue::Float(background_mean));
            layer.push(f);
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
        outputs.insert("dark_area_count".to_string(), json!(count));
        outputs.insert("dark_area".to_string(), json!(total_area));
        outputs.insert("threshold_db".to_string(), json!(threshold));
        outputs.insert("background_mean_db".to_string(), json!(background_mean));
        outputs.insert("background_sd_db".to_string(), json!(background_sd));
        outputs.insert("method".to_string(), json!(prm.method.label()));
        Ok(ToolRunResult { outputs })
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Method {
    Statistical,
    Otsu,
    Fixed,
}

impl Method {
    fn label(self) -> &'static str {
        match self {
            Method::Statistical => "statistical",
            Method::Otsu => "otsu",
            Method::Fixed => "fixed",
        }
    }
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    method: Method,
    threshold: f64,
    min_area: Option<f64>,
    min_contrast_db: Option<f64>,
    mask_side: MaskSide,
    bins: usize,
    input_units: SarUnits,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let method = match choice_or(
        args,
        "method",
        &["statistical", "otsu", "fixed"],
        "statistical",
    )? {
        "otsu" => Method::Otsu,
        "fixed" => Method::Fixed,
        _ => Method::Statistical,
    };
    // The default only makes sense for the statistical form; 'fixed' reads the
    // same parameter as an absolute dB level, so it must be supplied.
    let threshold = match method {
        Method::Fixed => opt_f64(args, "threshold")?.ok_or_else(|| {
            ToolError::Validation(
                "method 'fixed' needs 'threshold' as the dB level at or below which a cell is dark"
                    .to_string(),
            )
        })?,
        _ => f64_or(args, "threshold", 2.0)?,
    };
    if !threshold.is_finite() {
        return Err(ToolError::Validation(
            "'threshold' must be a finite number".to_string(),
        ));
    }
    if method == Method::Statistical && threshold <= 0.0 {
        return Err(ToolError::Validation(format!(
            "method 'statistical' needs a positive 'threshold' in standard deviations, got {threshold}"
        )));
    }

    let min_area = opt_positive_f64(args, "min_area")?;
    let min_contrast_db = match opt_f64(args, "min_contrast_db")? {
        None => None,
        Some(v) if v >= 0.0 && v.is_finite() => Some(v),
        Some(v) => {
            return Err(ToolError::Validation(format!(
                "'min_contrast_db' must be non-negative, got {v}"
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
        method,
        threshold,
        min_area,
        min_contrast_db,
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
    use wbraster::Raster;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn jitter(i: usize) -> f64 {
        let mut x = (i as u64).wrapping_mul(6364136223846793005).wrapping_add(23);
        x ^= x >> 33;
        x = x.wrapping_mul(0xff51afd7ed558ccd);
        x ^= x >> 33;
        (x % 4001) as f64 / 4001.0
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

    /// Sea at about -12 dB with a slick 8 dB below it.
    fn slick_scene(rows: usize, cols: usize) -> Vec<f64> {
        let sea = 0.063_f64; // about -12 dB
        let mut v: Vec<f64> = (0..rows * cols)
            .map(|i| sea * (0.85 + 0.3 * jitter(i)))
            .collect();
        for r in 8..18 {
            for c in 10..24 {
                v[r * cols + c] = sea * 0.16; // about -8 dB contrast
            }
        }
        v
    }

    fn run(args: Value) -> (Layer, BTreeMap<String, Value>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = DetectDarkOceanAreasTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (layer, out.outputs)
    }

    /// The slick is found and its contrast is reported in dB.
    #[test]
    fn finds_a_slick_and_reports_contrast() {
        let (rows, cols) = (40, 40);
        let (layer, outputs) = run(json!({
            "input": raster_of(cols, rows, &slick_scene(rows, cols)),
            "min_area": 2000.0
        }));
        assert_eq!(layer.len(), 1, "expected exactly one slick");
        let ci = layer.schema.field_index("contrast_db").unwrap();
        let FieldValue::Float(contrast) = layer.iter().next().unwrap().attributes[ci] else {
            panic!("contrast must be a float")
        };
        // The slick is 0.16x the sea power; against a background dominated by
        // open sea that is roughly 7-8 dB of contrast.
        assert!(
            (5.0..10.0).contains(&contrast),
            "implausible slick contrast {contrast} dB"
        );
        let area = outputs["dark_area"].as_f64().unwrap();
        assert!(
            (area - 14_000.0).abs() < 2000.0,
            "mapped slick area {area} should be near the true 14000 m^2"
        );
    }

    /// The reason this tool does not default to Otsu: on an all-water scene
    /// there is no second mode, so an Otsu split bisects the clutter and calls
    /// a large fraction of open sea a slick. The statistical default must not.
    #[test]
    fn statistical_default_beats_otsu_on_a_clean_sea() {
        let (rows, cols) = (40, 40);
        // No slick at all — just sea clutter.
        let sea = 0.063_f64;
        let v: Vec<f64> = (0..rows * cols)
            .map(|i| sea * (0.85 + 0.3 * jitter(i)))
            .collect();
        let src = raster_of(cols, rows, &v);
        let total = (rows * cols) as f64 * 100.0;

        let (_, stat) = run(json!({"input": src.clone(), "min_area": 2000.0}));
        let stat_area = stat["dark_area"].as_f64().unwrap();
        assert!(
            stat_area < 0.05 * total,
            "statistical method flagged {stat_area} of {total} on a clean sea"
        );

        let (_, otsu) = run(json!({"input": src, "method": "otsu", "min_area": 2000.0}));
        let otsu_area = otsu["dark_area"].as_f64().unwrap();
        assert!(
            otsu_area > 5.0 * stat_area.max(1.0),
            "otsu should over-report on a unimodal scene ({otsu_area} vs {stat_area}); \
             if it no longer does, the module doc's rationale needs revisiting"
        );
    }

    /// The contrast filter separates a real slick from a shallow wind shadow.
    #[test]
    fn contrast_filter_rejects_look_alikes() {
        let (rows, cols) = (40, 50);
        let sea = 0.063_f64;
        let mut v: Vec<f64> = (0..rows * cols)
            .map(|i| sea * (0.9 + 0.2 * jitter(i)))
            .collect();
        // A genuine slick, ~8 dB down.
        for r in 6..16 {
            for c in 5..17 {
                v[r * cols + c] = sea * 0.16;
            }
        }
        // A wind shadow, only ~3.5 dB down — dark enough to be detected, too
        // shallow to be mineral oil.
        for r in 24..34 {
            for c in 30..42 {
                v[r * cols + c] = sea * 0.45;
            }
        }
        let src = raster_of(cols, rows, &v);

        // The slick itself inflates the background spread, so a 1-sigma cut
        // would already miss the shallow patch; 0.8 catches both.
        let (both, _) = run(json!({
            "input": src.clone(), "threshold": 0.8, "min_area": 2000.0
        }));
        assert_eq!(both.len(), 2, "both dark patches should be found first");

        let (strong, _) = run(json!({
            "input": src, "threshold": 0.8, "min_area": 2000.0, "min_contrast_db": 5.0
        }));
        assert_eq!(
            strong.len(),
            1,
            "the contrast filter should keep only the real slick"
        );
    }

    /// The mask raster and the polygons describe the same areas.
    #[test]
    fn mask_raster_agrees_with_polygons() {
        let (rows, cols) = (40, 40);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": raster_of(cols, rows, &slick_scene(rows, cols)), "min_area": 2000.0
        }))
        .unwrap();
        let out = DetectDarkOceanAreasTool.run(&args, &ctx()).unwrap();
        let mask = load_input_raster(out.outputs["output_raster"].as_str().unwrap()).unwrap();
        let flagged = (0..rows)
            .flat_map(|r| (0..cols).map(move |c| (r, c)))
            .filter(|&(r, c)| mask.get(0, r as isize, c as isize) == 1.0)
            .count();
        let area = out.outputs["dark_area"].as_f64().unwrap();
        assert!((area - flagged as f64 * 100.0).abs() < 1e-6);
    }

    /// A fixed dB level is used verbatim.
    #[test]
    fn fixed_method_uses_the_given_level() {
        let (rows, cols) = (30, 30);
        let (_, outputs) = run(json!({
            "input": raster_of(cols, rows, &slick_scene(rows, cols)),
            "method": "fixed", "threshold": -17.0
        }));
        assert_eq!(outputs["method"].as_str().unwrap(), "fixed");
        assert!((outputs["threshold_db"].as_f64().unwrap() + 17.0).abs() < 1e-9);
    }

    /// A perfectly uniform sea has no anomaly; say so rather than flag half of
    /// it.
    #[test]
    fn uniform_sea_is_an_error() {
        let args: ToolArgs =
            serde_json::from_value(json!({"input": raster_of(8, 8, &[0.06; 64])})).unwrap();
        let err = DetectDarkOceanAreasTool.run(&args, &ctx()).unwrap_err();
        assert!(
            format!("{err:?}").contains("uniform"),
            "expected a uniform-scene error, got {err:?}"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            DetectDarkOceanAreasTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({"input": "a.tif", "method": "magic"})).is_err());
        // 'fixed' has no sensible default level.
        assert!(bad(json!({"input": "a.tif", "method": "fixed"})).is_err());
        assert!(bad(json!({"input": "a.tif", "threshold": -1.0})).is_err());
        assert!(bad(json!({"input": "a.tif", "min_contrast_db": -2})).is_err());
        assert!(bad(json!({"input": "a.tif", "min_area": -1})).is_err());
        assert!(bad(json!({"input": "a.tif"})).is_ok());
        assert!(bad(json!({"input": "a.tif", "method": "fixed", "threshold": -20})).is_ok());
    }
}
