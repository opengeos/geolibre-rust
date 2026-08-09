//! GeoLibre tool: terrain-flattened SAR backscatter from a DEM.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Apply Radiometric Terrain
//! Flattening* (Image Analyst).
//!
//! ## Why the catalog needs it
//!
//! Round 17's `apply_radiometric_calibration` normalises digital numbers to
//! sigma nought or gamma nought using the **ellipsoid** incidence angle — a
//! flat-earth assumption. Over any real terrain that assumption fails badly:
//! slopes tilted toward the sensor are compressed into fewer pixels and look
//! bright, slopes tilted away are stretched and look dark, and the resulting
//! pattern tracks topography rather than land cover. On mountainous scenes the
//! effect routinely exceeds 6 dB, which is larger than the difference between
//! forest and bare soil.
//!
//! Terrain flattening replaces the flat-earth reference area with the true
//! local illuminated area derived from a DEM. Neither registry has anything
//! for it: grepping both for `terrain_flat`, `rtc` and `radiometric` finds only
//! round 17's ellipsoid-based calibration.
//!
//! ## Method
//!
//! Slope and aspect come from Horn's 3x3 gradient (shared with
//! `solar_radiation`). The local incidence angle between the terrain normal and
//! the direction to the sensor is
//!
//! ```text
//! cos(theta_loc) = cos(theta) * cos(slope)
//!                + sin(theta) * sin(slope) * cos(look_azimuth - uphill_aspect)
//! ```
//!
//! and the flattened backscatter follows from replacing `theta` with
//! `theta_loc` in the usual normalisation:
//!
//! ```text
//! sigma0_flat = beta0 * sin(theta_loc)
//! gamma0_flat = beta0 * tan(theta_loc)
//! ```
//!
//! The scattering-area ratio `sin(theta) / sin(theta_loc)` — how much larger
//! the true illuminated area is than the flat-earth one — is written as a
//! secondary output, along with a layover/shadow mask (layover where
//! `theta_loc <= 0`, radar shadow where `theta_loc >= 90 degrees`; in both
//! cases no valid backscatter can be recovered and the cell is no-data).
//!
//! ## Deliberate scope
//!
//! This is the **local-incidence-angle** approximation, not the full
//! Small (2011) integration of DEM facet areas in radar geometry. The exact
//! method needs the sensor's orbit state vectors and a range-Doppler transform
//! between map and radar geometry; neither is available from a map-projected
//! DEM plus a raster, which is all this catalog can be handed. The
//! approximation is standard, is what most map-projected implementations use,
//! and removes the great majority of the topographic signal.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::args_common::{band_index, choice_or, f64_or, opt_f64, req_str};
use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};
use crate::raster_stack::check_alignment_refs;
use crate::sar_common::{power_to_db, SarUnits};
use crate::solar_radiation::slope_aspect;

pub struct RadiometricTerrainFlatteningTool;

impl Tool for RadiometricTerrainFlatteningTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "radiometric_terrain_flattening",
            display_name: "Radiometric Terrain Flattening",
            summary: "Removes the topographic brightness signal from calibrated SAR backscatter by replacing the flat-earth reference area with the DEM-derived local illuminated area, and flags layover and radar shadow (ArcGIS Apply Radiometric Terrain Flattening). Round 17's apply_radiometric_calibration normalises with the ellipsoid incidence angle, so on sloping ground its sigma0/gamma0 tracks topography rather than land cover — routinely by more than 6 dB, which exceeds the forest/bare-soil difference it is meant to measure.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Calibrated SAR backscatter raster (beta nought by default; see 'input_calibration').",
                    required: true,
                },
                ToolParamSpec {
                    name: "dem",
                    description: "Digital elevation model, co-registered with the input raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output terrain-flattened backscatter raster. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_scattering_area",
                    description: "Output ratio of the true illuminated area to the flat-earth area. Always produced; stored in memory when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_geometry_mask",
                    description: "Output mask: 0 normal, 1 layover, 2 radar shadow. Always produced; stored in memory when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "incidence_angle",
                    description: "Ellipsoid incidence angle in degrees: a constant, or the path to a per-cell raster. Required.",
                    required: true,
                },
                ToolParamSpec {
                    name: "look_azimuth",
                    description: "Compass azimuth in degrees from the target toward the sensor (default 270, i.e. a right-looking ascending pass). A right-looking descending Sentinel-1 pass is about 283.",
                    required: false,
                },
                ToolParamSpec {
                    name: "input_calibration",
                    description: "What the input already is: 'beta0' (default), 'sigma0_ellipsoid', or 'gamma0_ellipsoid'. The latter two are converted back to beta nought using the ellipsoid angle before flattening.",
                    required: false,
                },
                ToolParamSpec {
                    name: "calibration_type",
                    description: "Output convention: 'gamma0' (default) or 'sigma0'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "input_units",
                    description: "'intensity' (default), 'dn', 'amplitude', or 'db'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_units",
                    description: "'linear' (default) or 'db'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "z_factor",
                    description: "Multiplier converting DEM z units into the raster's XY units (default 1.0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "Band of the input raster to flatten (default 0).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        req_str(args, "dem")?;
        // Dual-typed: validate the constant form only, since the raster form
        // can only be checked once the grids are known.
        if incidence_path(args)?.is_none() {
            incidence_constant(args)?;
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input_path = req_str(args, "input")?.to_string();
        let dem_path = req_str(args, "dem")?.to_string();
        let prm = parse_params(args)?;
        let band = band_index(args, "band")?;
        let output = parse_optional_output(args, "output")?;
        let out_area = parse_optional_output(args, "output_scattering_area")?;
        let out_mask = parse_optional_output(args, "output_geometry_mask")?;

        let raster = load_input_raster(&input_path)?;
        let dem = load_input_raster(&dem_path)?;
        check_alignment_refs(&[&raster, &dem])?;
        let (rows, cols) = (raster.rows, raster.cols);

        // Per-cell ellipsoid incidence angle, from a constant or a raster.
        let theta = incidence_field(args, &raster, rows, cols)?;

        // Terrain derivatives. Horn's operator wants a single cell size; a
        // strongly anisotropic grid would bias the gradient, so require near-
        // square cells rather than silently averaging them.
        let (csx, csy) = (raster.cell_size_x, raster.cell_size_y);
        if (csx - csy).abs() > 1e-6 * csx.abs().max(csy.abs()).max(1.0) {
            return Err(ToolError::Validation(format!(
                "terrain flattening needs square cells; got {csx} x {csy}. Resample the DEM first."
            )));
        }
        let mut z = vec![f64::NAN; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let v = dem.get(0, r as isize, c as isize);
                if v != dem.nodata && v.is_finite() {
                    z[r * cols + c] = v * prm.z_factor;
                }
            }
        }
        let (slope, aspect_down) = slope_aspect(&z, rows, cols, csx);

        ctx.progress.info(&format!(
            "{rows}x{cols}, {} -> {}, look azimuth {} deg",
            prm.input_calibration.label(),
            prm.output_calibration.label(),
            prm.look_azimuth.to_degrees()
        ));

        let nodata = -9999.0_f64;
        let mut out = vec![nodata; rows * cols];
        let mut area = vec![nodata; rows * cols];
        let mut mask = vec![nodata; rows * cols];
        let (mut n_layover, mut n_shadow, mut n_valid) = (0usize, 0usize, 0usize);

        for r in 0..rows {
            for c in 0..cols {
                let i = r * cols + c;
                let raw = raster.get(band, r as isize, c as isize);
                if raw == raster.nodata || !raw.is_finite() || !z[i].is_finite() {
                    continue;
                }
                let Some(power) = prm.input_units.to_power(raw) else {
                    continue;
                };
                let th = theta[i];
                if !th.is_finite() || th <= 0.0 || th >= PI / 2.0 {
                    continue;
                }

                // Undo whatever ellipsoid normalisation the input already has,
                // so flattening always starts from beta nought.
                let beta0 = match prm.input_calibration {
                    Calibration::Beta0 => power,
                    Calibration::Sigma0 => power / th.sin(),
                    Calibration::Gamma0 => power / th.tan(),
                };

                // Local incidence angle, projected into the range plane.
                //
                // Only the range direction is foreshortened — azimuth sampling
                // is unaffected by terrain — so the quantity that scales the
                // pixel's ground area is the slope component along the look
                // direction, not the full 3-D angle between the look vector and
                // the terrain normal. Using the 3-D angle would also make
                // layover undetectable: an unsigned angle is never negative, so
                // the "tilted past the beam" case could not be expressed.
                //
                // `slope_range` is positive where the ground descends toward
                // the sensor (the surface faces it) and negative where it
                // descends away.
                let slope_range =
                    (slope[i].tan() * (aspect_down[i] - prm.look_azimuth).cos()).atan();
                let theta_loc = th - slope_range;

                // Layover: the slope leans past the look direction, so several
                // ground cells fold into one range bin and no per-cell
                // backscatter exists. Shadow: the slope faces away steeply
                // enough that the sensor never illuminates it.
                if theta_loc <= 0.0 {
                    mask[i] = 1.0;
                    n_layover += 1;
                    continue;
                }
                if theta_loc >= PI / 2.0 {
                    mask[i] = 2.0;
                    n_shadow += 1;
                    continue;
                }
                mask[i] = 0.0;

                let sin_loc = theta_loc.sin();
                let flattened = match prm.output_calibration {
                    Calibration::Sigma0 => beta0 * sin_loc,
                    // gamma0 = sigma0 / cos(theta_loc)
                    Calibration::Gamma0 => beta0 * theta_loc.tan(),
                    Calibration::Beta0 => beta0,
                };
                area[i] = th.sin() / sin_loc;

                out[i] = if prm.db_output {
                    match power_to_db(flattened) {
                        Some(v) => v,
                        None => continue,
                    }
                } else {
                    flattened
                };
                n_valid += 1;
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        ctx.progress.info(&format!(
            "{n_valid} valid cell(s), {n_layover} layover, {n_shadow} shadow"
        ));

        let out_path = write_or_store_output(
            raster_like_with_data(&raster, out, nodata, DataType::F32)?,
            output,
        )?;
        // Secondary outputs are always emitted: gating on a supplied path would
        // silently drop them for in-memory callers.
        let area_path = write_or_store_output(
            raster_like_with_data(&raster, area, nodata, DataType::F32)?,
            out_area,
        )?;
        let mask_path = write_or_store_output(
            raster_like_with_data(&raster, mask, nodata, DataType::F32)?,
            out_mask,
        )?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("output_scattering_area".to_string(), json!(area_path));
        outputs.insert("output_geometry_mask".to_string(), json!(mask_path));
        outputs.insert(
            "calibration_type".to_string(),
            json!(prm.output_calibration.label()),
        );
        outputs.insert(
            "output_units".to_string(),
            json!(if prm.db_output { "db" } else { "linear" }),
        );
        outputs.insert("valid_cells".to_string(), json!(n_valid));
        outputs.insert("layover_cells".to_string(), json!(n_layover));
        outputs.insert("shadow_cells".to_string(), json!(n_shadow));
        Ok(ToolRunResult { outputs })
    }
}

/// Per-cell ellipsoid incidence angle in radians, from a constant or a raster.
///
/// The parameter is dual-typed, so a plain number is tried first and only a
/// value that is not a number is treated as a path. Parsing it unconditionally
/// as a float is how round 17's `flatten_interferogram` made its documented
/// raster form unusable.
fn incidence_field(
    args: &ToolArgs,
    template: &Raster,
    rows: usize,
    cols: usize,
) -> Result<Vec<f64>, ToolError> {
    // A path is a string that does not parse as a number; anything else — a
    // JSON number or a numeric string — is the constant form. Reading the
    // scalar unconditionally is how round 17's `flatten_interferogram` shipped
    // a documented raster parameter that could never be used.
    let Some(path) = incidence_path(args)? else {
        let deg = incidence_constant(args)?;
        return Ok(vec![deg.to_radians(); rows * cols]);
    };
    let inc = load_input_raster(path)?;
    check_alignment_refs(&[template, &inc])?;
    let mut out = vec![f64::NAN; rows * cols];
    for r in 0..rows {
        for c in 0..cols {
            let v = inc.get(0, r as isize, c as isize);
            if v != inc.nodata && v.is_finite() {
                out[r * cols + c] = v.to_radians();
            }
        }
    }
    Ok(out)
}

/// The `incidence_angle` parameter as a raster path, or `None` when it is a
/// constant.
///
/// Reads the raw JSON rather than going through `parse_optional_str`, which
/// rejects a non-string outright — the parameter is legitimately either a
/// number or a path.
fn incidence_path(args: &ToolArgs) -> Result<Option<&str>, ToolError> {
    Ok(args
        .get("incidence_angle")
        .and_then(serde_json::Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty() && s.parse::<f64>().is_err()))
}

/// The `incidence_angle` parameter as a constant in degrees.
fn incidence_constant(args: &ToolArgs) -> Result<f64, ToolError> {
    let deg = opt_f64(args, "incidence_angle")?.ok_or_else(|| {
        ToolError::Validation("missing required parameter 'incidence_angle'".to_string())
    })?;
    if !(0.0..90.0).contains(&deg) {
        return Err(ToolError::Validation(format!(
            "'incidence_angle' must be in [0, 90) degrees, got {deg}"
        )));
    }
    Ok(deg)
}

/// Which area normalisation a backscatter raster carries.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Calibration {
    Beta0,
    Sigma0,
    Gamma0,
}

impl Calibration {
    fn label(self) -> &'static str {
        match self {
            Calibration::Beta0 => "beta0",
            Calibration::Sigma0 => "sigma0",
            Calibration::Gamma0 => "gamma0",
        }
    }
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    input_calibration: Calibration,
    output_calibration: Calibration,
    input_units: SarUnits,
    db_output: bool,
    look_azimuth: f64,
    z_factor: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let input_calibration = match choice_or(
        args,
        "input_calibration",
        &["beta0", "sigma0_ellipsoid", "gamma0_ellipsoid"],
        "beta0",
    )? {
        "sigma0_ellipsoid" => Calibration::Sigma0,
        "gamma0_ellipsoid" => Calibration::Gamma0,
        _ => Calibration::Beta0,
    };
    let output_calibration = match choice_or(args, "calibration_type", &["gamma0", "sigma0"], "gamma0")?
    {
        "sigma0" => Calibration::Sigma0,
        _ => Calibration::Gamma0,
    };
    let input_units = SarUnits::parse(
        args.get("input_units")
            .and_then(serde_json::Value::as_str)
            .unwrap_or(""),
    )?;
    let db_output = choice_or(args, "output_units", &["linear", "db"], "linear")? == "db";

    let look_azimuth = f64_or(args, "look_azimuth", 270.0)?;
    if !look_azimuth.is_finite() {
        return Err(ToolError::Validation(
            "'look_azimuth' must be a finite number of degrees".to_string(),
        ));
    }

    let z_factor = match opt_f64(args, "z_factor")? {
        None => 1.0,
        Some(z) if z > 0.0 && z.is_finite() => z,
        Some(z) => {
            return Err(ToolError::Validation(format!(
                "'z_factor' must be positive, got {z}"
            )))
        }
    };

    Ok(Params {
        input_calibration,
        output_calibration,
        input_units,
        db_output,
        look_azimuth: look_azimuth.to_radians(),
        z_factor,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
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

    fn run(args: Value) -> (Raster, BTreeMap<String, Value>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = RadiometricTerrainFlatteningTool.run(&args, &ctx()).unwrap();
        let r = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (r, out.outputs)
    }

    /// Over flat ground the local incidence angle equals the ellipsoid angle,
    /// so flattening must reproduce the plain gamma nought exactly — this is
    /// the analytic anchor the whole tool is checked against.
    #[test]
    fn flat_terrain_reproduces_ellipsoid_gamma0() {
        let (rows, cols) = (5, 5);
        let beta0 = 0.25;
        let sar = raster_of(cols, rows, &[beta0; 25]);
        let dem = raster_of(cols, rows, &[100.0; 25]); // perfectly flat
        let (out, outputs) = run(json!({
            "input": sar, "dem": dem, "incidence_angle": 35.0
        }));
        let want = beta0 * 35.0_f64.to_radians().tan();
        let got = out.get(0, 2, 2);
        assert!(
            (got - want).abs() < 1e-6 * want,
            "flat ground: {got} != ellipsoid gamma0 {want}"
        );
        // No topography means no layover and no shadow.
        assert_eq!(outputs["layover_cells"].as_u64().unwrap(), 0);
        assert_eq!(outputs["shadow_cells"].as_u64().unwrap(), 0);
        // ...and a scattering-area ratio of exactly 1.
        let area =
            load_input_raster(outputs["output_scattering_area"].as_str().unwrap()).unwrap();
        assert!((area.get(0, 2, 2) - 1.0).abs() < 1e-6);
    }

    /// The whole point of the tool: a uniform-beta0 scene over a hillside must
    /// come out *flatter* than the ellipsoid product, because the topographic
    /// brightness pattern is what gets removed.
    ///
    /// Slopes facing the sensor and slopes facing away have different local
    /// incidence angles, so the ellipsoid gamma0 is constant while the true
    /// scattering areas differ — a scene of constant *radar brightness* must
    /// therefore produce a varying flattened backscatter that mirrors the
    /// terrain, and the scattering-area ratio must exceed 1 on far-facing
    /// slopes and fall below 1 on near-facing ones.
    #[test]
    fn slope_facing_the_sensor_differs_from_slope_facing_away() {
        let (rows, cols) = (5, 9);
        // A symmetric ridge running north-south: west flank rises east, east
        // flank falls. Look azimuth 270 (sensor due west).
        let mut z = vec![0.0; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let d = (c as i32 - 4).abs() as f64;
                z[r * cols + c] = 100.0 - 5.0 * d;
            }
        }
        let sar = raster_of(cols, rows, &[0.25; 45]);
        let dem = raster_of(cols, rows, &z);
        let (_, outputs) = run(json!({
            "input": sar, "dem": dem, "incidence_angle": 35.0, "look_azimuth": 270.0
        }));
        let area =
            load_input_raster(outputs["output_scattering_area"].as_str().unwrap()).unwrap();

        // Column 2 is on the west flank, which descends toward the sensor in
        // the west and therefore *faces* it; column 6 is on the east flank,
        // facing away. A slope facing the sensor is foreshortened — many metres
        // of ground squeeze into few range bins — so each pixel covers more
        // ground than the flat-earth assumption allows and the ratio exceeds 1.
        // A slope facing away is stretched, so its ratio falls below 1.
        let west = area.get(0, 2, 2);
        let east = area.get(0, 2, 6);
        assert!(
            west > 1.001,
            "foreshortened near slope should enlarge the illuminated area, got {west}"
        );
        assert!(
            east < 0.999,
            "stretched far slope should shrink the illuminated area, got {east}"
        );

        // And the correction must act in the opposite sense to the brightness
        // error it removes: the over-bright near slope is pulled down.
        let (out, _) = run(json!({
            "input": raster_of(cols, rows, &[0.25; 45]),
            "dem": raster_of(cols, rows, &z),
            "incidence_angle": 35.0, "look_azimuth": 270.0
        }));
        let ellipsoid = 0.25 * 35.0_f64.to_radians().tan();
        assert!(
            out.get(0, 2, 2) < ellipsoid,
            "near slope should be darkened relative to the ellipsoid product"
        );
        assert!(
            out.get(0, 2, 6) > ellipsoid,
            "far slope should be brightened relative to the ellipsoid product"
        );
    }

    /// A slope steeper than the incidence angle and tilted into the look
    /// direction is layover: no per-cell backscatter exists there.
    #[test]
    fn detects_layover_on_steep_near_slopes() {
        let (rows, cols) = (3, 7);
        // Ground descending toward the west, i.e. facing the sensor at azimuth
        // 270, at 45 degrees — well past a 20 degree incidence angle, so the
        // surface has tilted through the beam and folds in range.
        let mut z = vec![0.0; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                z[r * cols + c] = 10.0 * c as f64;
            }
        }
        let sar = raster_of(cols, rows, &[0.25; 21]);
        let dem = raster_of(cols, rows, &z);
        let (_, outputs) = run(json!({
            "input": sar, "dem": dem, "incidence_angle": 20.0, "look_azimuth": 270.0
        }));
        assert!(
            outputs["layover_cells"].as_u64().unwrap() > 0,
            "a 45-degree slope into a 20-degree look must produce layover"
        );
        let mask =
            load_input_raster(outputs["output_geometry_mask"].as_str().unwrap()).unwrap();
        assert_eq!(mask.get(0, 1, 3), 1.0, "interior cell should be flagged 1");
    }

    /// Input already normalised on the ellipsoid is converted back to beta
    /// nought first, so all three input conventions agree on flat ground.
    #[test]
    fn input_calibration_round_trips() {
        let (rows, cols) = (3, 3);
        let theta = 40.0_f64.to_radians();
        let beta0 = 0.4;
        let dem = raster_of(cols, rows, &[50.0; 9]);
        let from_beta = run(json!({
            "input": raster_of(cols, rows, &[beta0; 9]), "dem": dem.clone(),
            "incidence_angle": 40.0
        }))
        .0
        .get(0, 1, 1);
        let from_sigma = run(json!({
            "input": raster_of(cols, rows, &[beta0 * theta.sin(); 9]), "dem": dem.clone(),
            "incidence_angle": 40.0, "input_calibration": "sigma0_ellipsoid"
        }))
        .0
        .get(0, 1, 1);
        let from_gamma = run(json!({
            "input": raster_of(cols, rows, &[beta0 * theta.tan(); 9]), "dem": dem,
            "incidence_angle": 40.0, "input_calibration": "gamma0_ellipsoid"
        }))
        .0
        .get(0, 1, 1);
        assert!(
            (from_beta - from_sigma).abs() < 1e-5 * from_beta
                && (from_beta - from_gamma).abs() < 1e-5 * from_beta,
            "conventions disagree: beta {from_beta}, sigma {from_sigma}, gamma {from_gamma}"
        );
    }

    /// dB output is the decibel form of the linear answer.
    #[test]
    fn db_output_matches_linear() {
        let (rows, cols) = (3, 3);
        let dem = raster_of(cols, rows, &[0.0; 9]);
        let lin = run(json!({
            "input": raster_of(cols, rows, &[0.3; 9]), "dem": dem.clone(),
            "incidence_angle": 30.0
        }))
        .0
        .get(0, 1, 1);
        let db = run(json!({
            "input": raster_of(cols, rows, &[0.3; 9]), "dem": dem,
            "incidence_angle": 30.0, "output_units": "db"
        }))
        .0
        .get(0, 1, 1);
        assert!(
            (db - 10.0 * lin.log10()).abs() < 1e-4,
            "dB {db} does not match linear {lin}"
        );
    }

    /// A per-cell incidence-angle raster must actually be usable — a documented
    /// parameter with no test is how round 17 shipped a broken one.
    #[test]
    fn accepts_a_raster_incidence_angle() {
        let (rows, cols) = (3, 3);
        // 30 degrees on the left column, 50 on the right.
        let inc: Vec<f64> = (0..9).map(|i| if i % 3 == 0 { 30.0 } else { 50.0 }).collect();
        let sar = raster_of(cols, rows, &[0.25; 9]);
        let dem = raster_of(cols, rows, &[0.0; 9]);
        let (out, _) = run(json!({
            "input": sar, "dem": dem, "incidence_angle": raster_of(cols, rows, &inc)
        }));
        let left = out.get(0, 1, 0);
        let right = out.get(0, 1, 2);
        assert!(
            (left - 0.25 * 30.0_f64.to_radians().tan()).abs() < 1e-6,
            "left cell used the wrong angle: {left}"
        );
        assert!(
            (right - 0.25 * 50.0_f64.to_radians().tan()).abs() < 1e-6,
            "right cell used the wrong angle: {right}"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            RadiometricTerrainFlatteningTool.validate(&args)
        };
        assert!(bad(json!({"dem": "d.tif", "incidence_angle": 30})).is_err());
        assert!(bad(json!({"input": "a.tif", "incidence_angle": 30})).is_err());
        assert!(bad(json!({"input": "a.tif", "dem": "d.tif"})).is_err());
        let ok_base = json!({"input": "a.tif", "dem": "d.tif", "incidence_angle": 30});
        assert!(bad(ok_base.clone()).is_ok());
        let mut m = ok_base.as_object().unwrap().clone();
        m.insert("calibration_type".into(), json!("beta0"));
        assert!(bad(Value::Object(m.clone())).is_err());
        m.insert("calibration_type".into(), json!("sigma0"));
        m.insert("z_factor".into(), json!(-1.0));
        assert!(bad(Value::Object(m.clone())).is_err());
        m.insert("z_factor".into(), json!(1.0));
        m.insert("input_units".into(), json!("watts"));
        assert!(bad(Value::Object(m)).is_err());
    }
}
