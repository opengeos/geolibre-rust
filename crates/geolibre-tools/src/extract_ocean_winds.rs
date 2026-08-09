//! GeoLibre tool: ocean surface wind speed from calibrated SAR backscatter.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Extract Ocean Winds* (Image Analyst),
//! implemented with the CMOD5.N geophysical model function.
//!
//! ## The next link in the SAR chain
//!
//! Rounds 17-18 built the radiometric chain — `apply_radiometric_calibration`,
//! `convert_sar_units`, `radiometric_terrain_flattening`, `multilook`, the
//! speckle filters, `detect_dark_ocean_areas`, `detect_bright_ocean_objects`,
//! `extract_water_sar`. Every one of those stays in the **radiometric** domain:
//! they measure and classify backscatter. None converts backscatter into a
//! **geophysical** quantity.
//!
//! Wind retrieval is the flagship SAR ocean product, and it is what makes the
//! dark patches `detect_dark_ocean_areas` finds interpretable at all: low wind
//! and an oil slick look identical in raw backscatter, and only a wind field
//! tells them apart.
//!
//! ## CMOD5.N is a closed-form function, not a model file
//!
//! CMOD5.N (Hersbach 2007) expresses
//! `sigma0 = B0 * (1 + B1 cos(phi) + B2 cos(2 phi))^1.6`, where `B0`, `B1` and
//! `B2` are published polynomials in wind speed and incidence angle with 28
//! fixed coefficients. No training, no model file, no external service — which
//! is what makes it implementable here at all.
//!
//! ## Why wind direction is an input
//!
//! The GMF has **two** unknowns (speed and the wind direction relative to the
//! radar look) and a single-polarization scene supplies **one** observation per
//! cell. Direction is therefore a required prior, taken from a weather model or
//! supplied as a constant; it cannot be retrieved here. Direction retrieval
//! needs wind-streak image analysis or a scatterometer, both separate problems.
//!
//! ## Inversion, and why it is not a plain bisection
//!
//! CMOD5.N is monotonic in wind speed over the **operational** range (up to
//! about 25 m/s at every geometry), but it is *not* globally monotonic: at low
//! incidence angles the upwind curve saturates and turns over above roughly
//! 30 m/s, so sigma0 at 50 m/s can be *lower* than at the peak. A bisection over
//! the whole `[0.2, max_wind]` bracket would therefore converge on the wrong
//! side of that peak, and its endpoint test would reject observations it should
//! have solved.
//!
//! The inversion instead scans forward in 0.5 m/s steps and bisects inside the
//! **first rising interval** that brackets the observation, returning the lowest
//! consistent wind speed — the standard convention for C-band retrieval, where
//! the high branch is unreliable anyway. Bisection within that interval is
//! deliberate over Newton: no derivative, and it cannot diverge.
//!
//! A cell that no rising interval brackets is *not solvable* — written as
//! no-data and counted, never clamped, because clamping would fabricate a
//! plausible-looking wind field out of a calibration error.
//!
//! ## Scope, stated plainly
//!
//! * **VV polarization** is the core deliverable. `polarization = "hh"` applies
//!   a Thompson-style polarization ratio, which is an approximation and is
//!   documented as such rather than presented as equivalent.
//! * Output is the **10 m equivalent neutral wind** CMOD5.N is defined against,
//!   not a true 10 m wind under non-neutral stratification.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::args_common::{choice_or, opt_positive_f64, req_str};
use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};
use crate::raster_stack::check_alignment_refs;
use crate::sar_common::{rasterize_mask, MaskSide, SarUnits};
use crate::vector_common::load_input_layer;

const POLARIZATIONS: [&str; 2] = ["vv", "hh"];

/// Lower end of the inversion bracket. CMOD is not defined at zero wind and the
/// GMF flattens below roughly 0.5 m/s.
const MIN_WIND: f64 = 0.2;

/// Upper limit accepted for `max_wind`. It sets the per-cell scan length, and
/// no real retrieval goes near it (the strongest recorded surface winds are
/// about 30 m/s below it).
const MAX_WIND_CEILING: f64 = 100.0;

pub struct ExtractOceanWindsTool;

impl Tool for ExtractOceanWindsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "extract_ocean_winds",
            display_name: "Extract Ocean Winds",
            summary: "Retrieves 10 m equivalent-neutral ocean wind speed from calibrated SAR backscatter by inverting the CMOD5.N geophysical model function (ArcGIS Extract Ocean Winds). The whole shipped SAR chain stays in the radiometric domain and nothing converted backscatter into a geophysical quantity, which is also what separates a low-wind patch from an oil slick.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Calibrated backscatter (sigma0) raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "units",
                    description: "Units of 'input': 'dn', 'amplitude', 'intensity', or 'db' (default 'db').",
                    required: false,
                },
                ToolParamSpec {
                    name: "incidence_angle",
                    description: "Per-cell incidence angle in degrees: a raster path, or a constant number.",
                    required: true,
                },
                ToolParamSpec {
                    name: "wind_direction",
                    description: "Meteorological wind direction in degrees (the direction the wind blows FROM), as a raster path or a constant. A required prior: CMOD has two unknowns and one observation.",
                    required: true,
                },
                ToolParamSpec {
                    name: "look_direction",
                    description: "Radar look azimuth in degrees, as a raster path or a constant. Used to convert wind direction into the relative angle CMOD takes.",
                    required: true,
                },
                ToolParamSpec {
                    name: "polarization",
                    description: "'vv' (default, native CMOD5.N) or 'hh' (applies an approximate polarization ratio).",
                    required: false,
                },
                ToolParamSpec {
                    name: "polarization_ratio_alpha",
                    description: "Alpha in the HH-to-VV polarization ratio (default 0.6). Only used when polarization is 'hh'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "water_mask",
                    description: "Optional polygon layer or raster restricting output to water; land returns are meaningless here.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_wind",
                    description: "Upper end of the inversion bracket, m/s (default 50).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output wind-speed raster in m/s. If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        for k in ["incidence_angle", "wind_direction", "look_direction"] {
            if args.get(k).is_none() || matches!(args.get(k), Some(Value::Null)) {
                return Err(ToolError::Validation(format!(
                    "missing required parameter '{k}' (a raster path or a constant number)"
                )));
            }
        }
        if let Some(u) = args.get("units").and_then(Value::as_str) {
            SarUnits::parse(u)?;
        }
        match args.get("water_mask") {
            None | Some(Value::Null) | Some(Value::String(_)) => {}
            Some(other) => {
                return Err(ToolError::Validation(format!(
                    "'water_mask' must be a raster or vector path; got {other}"
                )))
            }
        }
        choice_or(args, "polarization", &POLARIZATIONS, "vv")?;
        if let Some(m) = opt_positive_f64(args, "max_wind")? {
            // The scan length is (max_wind - MIN_WIND) / SCAN_STEP evaluations
            // PER CELL, so an unbounded value does not fail — it stops
            // responding. 100 m/s is above any real retrieval.
            if m <= MIN_WIND || m > MAX_WIND_CEILING {
                return Err(ToolError::Validation(format!(
                    "'max_wind' must be between {MIN_WIND} and {MAX_WIND_CEILING} m/s, got {m}"
                )));
            }
        }
        if let Some(a) = opt_positive_f64(args, "polarization_ratio_alpha")? {
            if !a.is_finite() {
                return Err(ToolError::Validation(
                    "'polarization_ratio_alpha' must be finite".to_string(),
                ));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?.to_string();
        let units = match args.get("units").and_then(Value::as_str) {
            Some(u) => SarUnits::parse(u)?,
            None => SarUnits::parse("db")?,
        };
        let polarization = choice_or(args, "polarization", &POLARIZATIONS, "vv")?;
        let alpha = opt_positive_f64(args, "polarization_ratio_alpha")?.unwrap_or(0.6);
        let max_wind = opt_positive_f64(args, "max_wind")?.unwrap_or(50.0);

        let sigma = load_input_raster(&input)?;
        let rows = sigma.rows;
        let cols = sigma.cols;

        let incidence = Field::load(args, "incidence_angle", &sigma)?;
        let wind_dir = Field::load(args, "wind_direction", &sigma)?;
        let look_dir = Field::load(args, "look_direction", &sigma)?;

        // Water mask: a polygon layer via rasterize_mask, or any non-zero cell
        // of a raster.
        let mask: Option<Vec<bool>> = match args.get("water_mask").and_then(Value::as_str) {
            None => None,
            Some(spec) if spec.trim().is_empty() => None,
            Some(spec) => Some(match load_input_raster(spec) {
                Ok(r) => {
                    check_alignment_refs(&[&sigma, &r])?;
                    (0..rows * cols)
                        .map(|i| {
                            let v = r.get(0, (i / cols) as isize, (i % cols) as isize);
                            v != r.nodata && v.is_finite() && v != 0.0
                        })
                        .collect()
                }
                Err(_) => {
                    // The parameter is a WATER mask, so its polygons delineate
                    // water and the cells to analyse are the ones inside them.
                    // (MaskSide has no "inside" variant; passing one made this
                    // whole branch fail at runtime.)
                    let layer = load_input_layer(spec)?;
                    // rasterize_mask compares cell-centre coordinates directly
                    // against the polygon coordinates, so a mask in a different
                    // CRS masks the wrong cells (or none) without any error.
                    match (sigma.crs.epsg, layer.crs_epsg()) {
                        (Some(img), Some(mask_epsg)) if img != mask_epsg => {
                            return Err(ToolError::Validation(format!(
                                "'water_mask' is EPSG:{mask_epsg} but 'input' is EPSG:{img}; \
                                 reproject the mask to the image CRS first"
                            )))
                        }
                        (Some(_), Some(_)) => {}
                        // One or both CRSs are undeclared, which is legitimate
                        // (a GeoJSON without a `crs` member is CRS84 by spec),
                        // so requiring the metadata would reject valid input.
                        // The harm a mismatch actually causes is masking the
                        // wrong cells, and the detectable form of that is a
                        // mask whose extent does not meet the raster's at all.
                        _ => {
                            let r_x1 = sigma.x_min + sigma.cols as f64 * sigma.cell_size_x;
                            let r_y1 = sigma.y_min + sigma.rows as f64 * sigma.cell_size_y;
                            if let Some((mx0, my0, mx1, my1)) = layer_bounds(&layer) {
                                if mx1 < sigma.x_min
                                    || mx0 > r_x1
                                    || my1 < sigma.y_min
                                    || my0 > r_y1
                                {
                                    return Err(ToolError::Validation(format!(
                                        "'water_mask' spans ({mx0:.3}, {my0:.3})-({mx1:.3}, \
                                         {my1:.3}) but the image spans ({:.3}, {:.3})-({r_x1:.3}, \
                                         {r_y1:.3}); they do not overlap, which usually means the \
                                         two are in different coordinate systems. Declare or \
                                         reproject the CRS of both.",
                                        sigma.x_min, sigma.y_min
                                    )));
                                }
                            }
                        }
                    }
                    rasterize_mask(&sigma, &layer, MaskSide::WaterPolygon)
                }
            }),
        };

        ctx.progress.info(&format!(
            "{rows}x{cols}, CMOD5.N inversion, {polarization} polarization"
        ));

        let nodata = -9999.0_f64;
        let mut wind = vec![nodata; rows * cols];
        let mut valid = 0_u64;
        let mut unsolvable = 0_u64;
        let mut masked = 0_u64;
        let (mut sum, mut lo_v, mut hi_v) = (0.0_f64, f64::MAX, f64::MIN);

        for row in 0..rows {
            for col in 0..cols {
                let idx = row * cols + col;
                if mask.as_ref().is_some_and(|m| !m[idx]) {
                    masked += 1;
                    continue;
                }
                let raw = sigma.get(0, row as isize, col as isize);
                if raw == sigma.nodata || !raw.is_finite() {
                    continue;
                }
                // `to_power` converts every supported unit to linear power,
                // which is what CMOD is defined on.
                let Some(mut s0) = units.to_power(raw) else {
                    continue;
                };
                let (Some(theta), Some(wdir), Some(ldir)) = (
                    incidence.get(&sigma, row, col),
                    wind_dir.get(&sigma, row, col),
                    look_dir.get(&sigma, row, col),
                ) else {
                    continue;
                };
                // CMOD is fitted over roughly 16-66 degrees; outside that the
                // polynomials are extrapolation, not measurement.
                if !(15.0..=70.0).contains(&theta) {
                    continue;
                }
                if polarization == "hh" {
                    // Thompson-style ratio: sigma0_VV = sigma0_HH * PR.
                    s0 *= polarization_ratio(theta, alpha);
                }
                // Relative wind direction: the angle between where the wind
                // comes from and where the radar looks.
                let phi = normalize_deg(wdir - ldir);

                match invert(s0, phi, theta, max_wind) {
                    Some(v) => {
                        wind[idx] = v;
                        valid += 1;
                        sum += v;
                        lo_v = lo_v.min(v);
                        hi_v = hi_v.max(v);
                    }
                    None => unsolvable += 1,
                }
            }
            ctx.progress
                .progress((row as f64 + 1.0) / rows.max(1) as f64);
        }

        let out = raster_like_with_data(&sigma, wind, nodata, DataType::F32)?;
        let out_path = write_or_store_output(out, parse_optional_output(args, "output")?)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("valid_cells".to_string(), json!(valid));
        outputs.insert("unsolvable_cells".to_string(), json!(unsolvable));
        outputs.insert("masked_cells".to_string(), json!(masked));
        outputs.insert("polarization".to_string(), json!(polarization));
        outputs.insert(
            "mean_wind_speed".to_string(),
            json!(if valid > 0 { sum / valid as f64 } else { 0.0 }),
        );
        outputs.insert(
            "min_wind_speed".to_string(),
            json!(if valid > 0 { lo_v } else { 0.0 }),
        );
        outputs.insert(
            "max_wind_speed".to_string(),
            json!(if valid > 0 { hi_v } else { 0.0 }),
        );
        Ok(ToolRunResult { outputs })
    }
}

/// A parameter that may be a raster or a single constant.
///
/// `vector_common::parse_optional_str` errors on a non-string value, so it
/// cannot express this dual type; the raw JSON is read instead (the round-18
/// lesson). A string that parses as a number is treated as a constant, so
/// `"35"` and `35` behave the same.
enum Field {
    Constant(f64),
    Raster(Box<Raster>),
}

impl Field {
    fn load(args: &ToolArgs, key: &str, template: &Raster) -> Result<Self, ToolError> {
        match args.get(key) {
            Some(Value::Number(n)) => n
                .as_f64()
                .filter(|v| v.is_finite())
                .map(Field::Constant)
                .ok_or_else(|| ToolError::Validation(format!("'{key}' must be a finite number"))),
            Some(Value::String(s)) => {
                let t = s.trim();
                if t.is_empty() {
                    return Err(ToolError::Validation(format!("'{key}' is empty")));
                }
                if let Ok(v) = t.parse::<f64>() {
                    if v.is_finite() {
                        return Ok(Field::Constant(v));
                    }
                }
                let r = load_input_raster(t)?;
                check_alignment_refs(&[template, &r])?;
                Ok(Field::Raster(Box::new(r)))
            }
            _ => Err(ToolError::Validation(format!(
                "'{key}' must be a raster path or a number"
            ))),
        }
    }

    fn get(&self, _template: &Raster, row: usize, col: usize) -> Option<f64> {
        match self {
            Field::Constant(v) => Some(*v),
            Field::Raster(r) => {
                let v = r.get(0, row as isize, col as isize);
                (v != r.nodata && v.is_finite()).then_some(v)
            }
        }
    }
}

/// Axis-aligned extent of a layer's geometry, or `None` when it has none.
fn layer_bounds(layer: &wbvector::Layer) -> Option<(f64, f64, f64, f64)> {
    let mut b = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
    let mut seen = false;
    for f in layer.iter() {
        if let Some(g) = f.geometry.as_ref() {
            let (x0, y0, x1, y1) = crate::sar_common::envelope(g);
            if x0.is_finite() && y0.is_finite() {
                seen = true;
                b.0 = b.0.min(x0);
                b.1 = b.1.min(y0);
                b.2 = b.2.max(x1);
                b.3 = b.3.max(y1);
            }
        }
    }
    seen.then_some(b)
}

/// Wraps an angle difference into [0, 360).
fn normalize_deg(d: f64) -> f64 {
    let m = d % 360.0;
    if m < 0.0 {
        m + 360.0
    } else {
        m
    }
}

/// HH-to-VV polarization ratio (Thompson form):
/// `PR = (1 + 2 tan^2 theta)^2 / (1 + alpha tan^2 theta)^2`.
///
/// An approximation, not an equivalence — see the module doc.
fn polarization_ratio(theta_deg: f64, alpha: f64) -> f64 {
    let t = theta_deg.to_radians().tan();
    let t2 = t * t;
    let num = 1.0 + 2.0 * t2;
    let den = 1.0 + alpha * t2;
    (num * num) / (den * den)
}

/// CMOD5.N coefficients (Hersbach 2007), 1-indexed as published so the formulas
/// below can be read against the paper directly. Index 0 is unused padding.
const C: [f64; 29] = [
    0.0, -0.6878, -0.7957, 0.3380, -0.1728, 0.0000, 0.0040, 0.1103, 0.0159, 6.7329, 2.7713,
    -2.2885, 0.4971, -0.7250, 0.0450, 0.0066, 0.3222, 0.0120, 22.7000, 2.0813, 3.0000, 8.3659,
    -3.3428, 1.3236, 6.2437, 2.3893, 0.3249, 4.1590, 1.6930,
];

/// The CMOD5.N forward model: linear sigma0 for a wind speed `v` (m/s), a
/// relative wind direction `phi` (degrees) and an incidence angle `theta`
/// (degrees).
fn cmod5n_forward(v: f64, phi_deg: f64, theta_deg: f64) -> f64 {
    const ZPOW: f64 = 1.6;
    const THETM: f64 = 40.0;
    const THETHR: f64 = 25.0;

    let y0 = C[19];
    let pn = C[20];
    let a = y0 - (y0 - 1.0) / pn;
    let b = 1.0 / (pn * (y0 - 1.0).powf(pn - 1.0));

    let fi = phi_deg.to_radians();
    let csfi = fi.cos();
    let cs2fi = 2.0 * csfi * csfi - 1.0;

    let x = (theta_deg - THETM) / THETHR;
    let xx = x * x;

    // B0: wind speed and incidence angle.
    let a0 = C[1] + C[2] * x + C[3] * xx + C[4] * x * xx;
    let a1 = C[5] + C[6] * x;
    let a2 = C[7] + C[8] * x;
    let gam = C[9] + C[10] * x + C[11] * xx;
    let s0 = C[12] + C[13] * x;

    let s = a2 * v;
    let s_vec = if s < s0 { s0 } else { s };
    let mut a3 = 1.0 / (1.0 + (-s_vec).exp());
    if s < s0 {
        // The published low-wind continuation, which keeps B0 smooth as the
        // sigmoid argument falls below s0.
        a3 *= (s / s0).powf(s0 * (1.0 - a3));
    }
    let b0 = a3.powf(gam) * 10.0_f64.powf(a0 + a1 * v);

    // B1: the upwind/downwind asymmetry.
    let mut b1 = C[15] * v * (0.5 + x - (4.0 * (x + C[16] + C[17] * v)).tanh());
    b1 = C[14] * (1.0 + x) - b1;
    b1 /= (0.34 * (v - C[18])).exp() + 1.0;

    // B2: the upwind/crosswind modulation.
    let v0 = C[21] + C[22] * x + C[23] * xx;
    let d1 = C[24] + C[25] * x + C[26] * xx;
    let d2 = C[27] + C[28] * x;
    let mut v2 = v / v0 + 1.0;
    if v2 < y0 {
        v2 = a + b * (v2 - 1.0).powf(pn);
    }
    let b2 = (-d1 + d2 * v2) * (-v2).exp();

    // The directional term can in principle go negative at extreme B1/B2,
    // where a fractional power is NaN. Flooring at zero keeps the forward
    // model total, and the inversion then simply finds no bracket there.
    let base = 1.0 + b1 * csfi + b2 * cs2fi;
    if base <= 0.0 {
        return 0.0;
    }
    b0 * base.powf(ZPOW)
}

/// Coarse scan step, m/s, used to find the sub-interval holding the solution.
const SCAN_STEP: f64 = 0.5;

/// Inverts the GMF for wind speed.
///
/// **CMOD5.N is not globally monotonic in wind speed.** At low incidence angles
/// the upwind curve saturates and turns over above roughly 30 m/s, so a plain
/// bisection over the whole bracket can converge on the wrong side of the peak
/// (or reject a perfectly good observation because the endpoint value is below
/// the peak). This scans forward in `SCAN_STEP` increments and bisects inside
/// the **first** interval that brackets the observation on a rising segment,
/// which returns the lowest consistent wind speed — the standard convention for
/// C-band retrieval, where the high branch is unreliable anyway.
///
/// Returns `None` when no rising interval brackets the observation: the cell is
/// genuinely unsolvable, and clamping it to an endpoint would manufacture a
/// plausible-looking wind out of a calibration error.
fn invert(sigma0: f64, phi_deg: f64, theta_deg: f64, max_wind: f64) -> Option<f64> {
    if !sigma0.is_finite() || sigma0 <= 0.0 {
        return None;
    }
    let f = |v: f64| cmod5n_forward(v, phi_deg, theta_deg);

    let mut prev_v = MIN_WIND;
    let mut prev_f = f(prev_v);
    if sigma0 < prev_f {
        return None; // below what even the slowest wind produces
    }
    // Step count is computed up front rather than accumulating `v += step`: a
    // float loop that stops advancing is the round-18 non-termination trap.
    let steps = (((max_wind - MIN_WIND) / SCAN_STEP).ceil() as usize).max(1);
    for i in 1..=steps {
        let v = (MIN_WIND + i as f64 * SCAN_STEP).min(max_wind);
        let fv = f(v);
        if fv > prev_f && prev_f <= sigma0 && sigma0 <= fv {
            // Rising segment that brackets the observation: bisect inside it.
            let (mut lo, mut hi) = (prev_v, v);
            for _ in 0..40 {
                if hi - lo < 1e-3 {
                    break;
                }
                let mid = 0.5 * (lo + hi);
                if f(mid) < sigma0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            return Some(0.5 * (lo + hi));
        }
        prev_v = v;
        prev_f = fv;
    }
    None
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

    fn raster(rows: usize, cols: usize, data: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F64,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, data[row * cols + col])
                    .unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: Value) -> (Raster, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = ExtractOceanWindsTool.run(&args, &ctx()).unwrap();
        let out = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        (out, res)
    }

    #[test]
    fn the_forward_model_increases_over_the_operational_wind_range() {
        // Monotonic up to 25 m/s at every geometry, which is the range C-band
        // retrieval is actually usable over.
        for theta in [20.0, 35.0, 50.0] {
            for phi in [0.0, 90.0, 180.0] {
                let mut prev = f64::NEG_INFINITY;
                for i in 1..=50 {
                    let v = i as f64 * 0.5;
                    let s = cmod5n_forward(v, phi, theta);
                    assert!(
                        s > prev,
                        "not monotonic at v={v}, phi={phi}, theta={theta}: {s} <= {prev}"
                    );
                    prev = s;
                }
            }
        }
    }

    #[test]
    fn the_forward_model_saturates_and_turns_over_at_low_incidence() {
        // Documents the real GMF behaviour the inversion has to cope with:
        // upwind at 20 degrees, sigma0 peaks near 30 m/s and then FALLS. A
        // plain bisection over the whole bracket would land on the wrong side
        // of that peak, which is why invert() scans for the first rising
        // interval instead.
        let peak = (1..=100)
            .map(|i| {
                let v = i as f64 * 0.5;
                (v, cmod5n_forward(v, 0.0, 20.0))
            })
            .fold((0.0_f64, f64::NEG_INFINITY), |acc, x| {
                if x.1 > acc.1 {
                    x
                } else {
                    acc
                }
            });
        assert!(
            peak.0 < 45.0,
            "expected a turnover below the bracket end, peak at {} m/s",
            peak.0
        );
        assert!(
            cmod5n_forward(50.0, 0.0, 20.0) < peak.1,
            "the curve should fall away past its peak"
        );
    }

    #[test]
    fn inversion_picks_the_lowest_root_past_the_turnover() {
        // At a geometry with a turnover, an observation matched by two wind
        // speeds must resolve to the lower one rather than an arbitrary side.
        let truth = 8.0;
        let s0 = cmod5n_forward(truth, 0.0, 20.0);
        let got = invert(s0, 0.0, 20.0, 50.0).expect("should invert");
        assert!((got - truth).abs() < 0.05, "expected {truth}, got {got}");
    }

    #[test]
    fn the_forward_model_is_finite_across_the_domain() {
        for theta in [16.0, 25.0, 40.0, 55.0, 66.0] {
            for phi in [0.0, 45.0, 90.0, 135.0, 180.0, 270.0] {
                for v in [0.2, 1.0, 10.0, 25.0, 50.0] {
                    let s = cmod5n_forward(v, phi, theta);
                    assert!(
                        s.is_finite() && s >= 0.0,
                        "v={v} phi={phi} theta={theta}: {s}"
                    );
                }
            }
        }
    }

    #[test]
    fn upwind_backscatter_exceeds_crosswind() {
        // The B1/B2 directional terms exist to produce exactly this asymmetry;
        // if they were dropped the retrieval would be direction-blind.
        let up = cmod5n_forward(10.0, 0.0, 35.0);
        let cross = cmod5n_forward(10.0, 90.0, 35.0);
        assert!(up > cross, "upwind {up} should exceed crosswind {cross}");
    }

    #[test]
    fn inversion_recovers_the_wind_speed_the_forward_model_was_given() {
        // The round-trip that makes the tool trustworthy.
        for theta in [25.0, 35.0, 45.0] {
            for phi in [0.0, 60.0, 120.0, 180.0] {
                for truth in [3.0, 7.5, 12.0, 20.0] {
                    let s0 = cmod5n_forward(truth, phi, theta);
                    let got = invert(s0, phi, theta, 50.0).expect("should invert");
                    assert!(
                        (got - truth).abs() < 0.02,
                        "theta={theta} phi={phi}: expected {truth}, got {got}"
                    );
                }
            }
        }
    }

    #[test]
    fn an_out_of_bracket_backscatter_is_unsolvable_not_clamped() {
        // A sigma0 far above anything 50 m/s can produce must not silently
        // return 50 — that would turn a calibration error into a hurricane.
        assert!(invert(1e6, 0.0, 35.0, 50.0).is_none());
        assert!(invert(1e-12, 0.0, 35.0, 50.0).is_none());
        assert!(invert(-1.0, 0.0, 35.0, 50.0).is_none());
    }

    #[test]
    fn a_uniform_scene_retrieves_a_uniform_wind_field() {
        let truth = 9.0;
        let s0 = cmod5n_forward(truth, 45.0, 35.0);
        let db = 10.0 * s0.log10();
        let (out, res) = run(json!({
            "input": raster(2, 2, &[db; 4]), "units": "db",
            "incidence_angle": 35.0, "wind_direction": 45.0, "look_direction": 0.0,
        }));
        assert_eq!(res.outputs["valid_cells"], json!(4));
        for r in 0..2 {
            for c in 0..2 {
                assert!(
                    (out.get(0, r, c) - truth).abs() < 0.05,
                    "{}",
                    out.get(0, r, c)
                );
            }
        }
    }

    #[test]
    fn intensity_units_give_the_same_answer_as_db() {
        let truth = 12.0;
        let s0 = cmod5n_forward(truth, 30.0, 40.0);
        let common = json!({
            "incidence_angle": 40.0, "wind_direction": 30.0, "look_direction": 0.0,
        });
        let mut a = common.clone();
        a["input"] = json!(raster(1, 1, &[10.0 * s0.log10()]));
        a["units"] = json!("db");
        let mut b = common;
        b["input"] = json!(raster(1, 1, &[s0]));
        b["units"] = json!("intensity");
        let (ra, _) = run(a);
        let (rb, _) = run(b);
        assert!((ra.get(0, 0, 0) - rb.get(0, 0, 0)).abs() < 1e-3);
    }

    #[test]
    fn relative_direction_is_wind_minus_look() {
        // Same physical geometry expressed two ways must give one answer.
        let s0 = cmod5n_forward(8.0, 30.0, 35.0);
        let db = 10.0 * s0.log10();
        let img = raster(1, 1, &[db]);
        let (a, _) = run(json!({
            "input": img, "units": "db", "incidence_angle": 35.0,
            "wind_direction": 30.0, "look_direction": 0.0,
        }));
        let (b, _) = run(json!({
            "input": img, "units": "db", "incidence_angle": 35.0,
            "wind_direction": 120.0, "look_direction": 90.0,
        }));
        assert!((a.get(0, 0, 0) - b.get(0, 0, 0)).abs() < 1e-6);
    }

    #[test]
    fn a_negative_relative_angle_wraps_rather_than_failing() {
        let s0 = cmod5n_forward(8.0, 350.0, 35.0);
        let db = 10.0 * s0.log10();
        let (out, res) = run(json!({
            "input": raster(1, 1, &[db]), "units": "db", "incidence_angle": 35.0,
            "wind_direction": 10.0, "look_direction": 20.0,
        }));
        assert_eq!(res.outputs["valid_cells"], json!(1));
        assert!((out.get(0, 0, 0) - 8.0).abs() < 0.05);
    }

    #[test]
    fn per_cell_incidence_and_direction_rasters_are_honoured() {
        // Two cells at different incidence angles carrying the sigma0 that
        // 6 and 14 m/s produce there; both must come back correctly.
        let s_a = cmod5n_forward(6.0, 0.0, 25.0);
        let s_b = cmod5n_forward(14.0, 0.0, 45.0);
        let (out, _) = run(json!({
            "input": raster(1, 2, &[10.0 * s_a.log10(), 10.0 * s_b.log10()]),
            "units": "db",
            "incidence_angle": raster(1, 2, &[25.0, 45.0]),
            "wind_direction": raster(1, 2, &[0.0, 0.0]),
            "look_direction": 0.0,
        }));
        assert!(
            (out.get(0, 0, 0) - 6.0).abs() < 0.05,
            "{}",
            out.get(0, 0, 0)
        );
        assert!(
            (out.get(0, 0, 1) - 14.0).abs() < 0.05,
            "{}",
            out.get(0, 0, 1)
        );
    }

    #[test]
    fn incidence_angles_outside_the_fitted_range_are_skipped() {
        // CMOD is fitted over roughly 16-66 degrees; beyond that the
        // polynomials are extrapolation, not measurement.
        let s0 = cmod5n_forward(10.0, 0.0, 35.0);
        let (_, res) = run(json!({
            "input": raster(1, 1, &[10.0 * s0.log10()]), "units": "db",
            "incidence_angle": 5.0, "wind_direction": 0.0, "look_direction": 0.0,
        }));
        assert_eq!(res.outputs["valid_cells"], json!(0));
    }

    #[test]
    fn a_water_mask_excludes_land_cells() {
        let s0 = cmod5n_forward(10.0, 0.0, 35.0);
        let db = 10.0 * s0.log10();
        let (out, res) = run(json!({
            "input": raster(1, 2, &[db, db]), "units": "db",
            "incidence_angle": 35.0, "wind_direction": 0.0, "look_direction": 0.0,
            "water_mask": raster(1, 2, &[1.0, 0.0]),
        }));
        assert_eq!(res.outputs["masked_cells"], json!(1));
        assert_eq!(res.outputs["valid_cells"], json!(1));
        assert_eq!(out.get(0, 0, 1), out.nodata);
    }

    #[test]
    fn a_polygon_water_mask_excludes_land_cells() {
        // The vector branch takes a different path from the raster one: the
        // raster load fails, a layer is loaded, and rasterize_mask builds the
        // mask. Only the raster branch was covered.
        let s0 = cmod5n_forward(10.0, 0.0, 35.0);
        let db = 10.0 * s0.log10();
        let mut l = wbvector::Layer::new("water")
            .with_geom_type(wbvector::GeometryType::Polygon)
            .with_crs_epsg(3857);
        // Covers only the western half of a 1x2 grid (cells at x = 0.5, 1.5).
        l.add_feature(
            Some(wbvector::Geometry::polygon(
                vec![
                    wbvector::Coord::xy(0.0, 0.0),
                    wbvector::Coord::xy(1.0, 0.0),
                    wbvector::Coord::xy(1.0, 1.0),
                    wbvector::Coord::xy(0.0, 1.0),
                    wbvector::Coord::xy(0.0, 0.0),
                ],
                Vec::new(),
            )),
            &[],
        )
        .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let mask = wbvector::memory_store::make_vector_memory_path(&id);
        let (out, res) = run(json!({
            "input": raster(1, 2, &[db, db]), "units": "db",
            "incidence_angle": 35.0, "wind_direction": 0.0, "look_direction": 0.0,
            "water_mask": mask,
        }));
        assert_eq!(res.outputs["masked_cells"], json!(1));
        assert_eq!(res.outputs["valid_cells"], json!(1));
        assert_eq!(out.get(0, 0, 1), out.nodata, "the land cell must be masked");
    }

    #[test]
    fn a_water_mask_in_a_different_crs_is_refused() {
        // rasterize_mask compares raw coordinates, so a WGS84 mask over a
        // Web Mercator scene silently masks the wrong cells.
        let mut l = wbvector::Layer::new("water")
            .with_geom_type(wbvector::GeometryType::Polygon)
            .with_crs_epsg(4326);
        l.add_feature(
            Some(wbvector::Geometry::polygon(
                vec![
                    wbvector::Coord::xy(0.0, 0.0),
                    wbvector::Coord::xy(1.0, 0.0),
                    wbvector::Coord::xy(1.0, 1.0),
                    wbvector::Coord::xy(0.0, 0.0),
                ],
                Vec::new(),
            )),
            &[],
        )
        .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let mask = wbvector::memory_store::make_vector_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": raster(1, 2, &[-20.0, -20.0]), "units": "db",
            "incidence_angle": 35.0, "wind_direction": 0.0, "look_direction": 0.0,
            "water_mask": mask,
        }))
        .unwrap();
        let err = ExtractOceanWindsTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("EPSG:4326"), "{err}");
    }

    #[test]
    fn an_undeclared_mask_that_does_not_meet_the_image_is_refused() {
        // With no CRS on either side the EPSG check cannot fire, so a
        // disjoint extent is the detectable signature of a mismatch.
        let mut l = wbvector::Layer::new("water").with_geom_type(wbvector::GeometryType::Polygon);
        l.add_feature(
            Some(wbvector::Geometry::polygon(
                vec![
                    wbvector::Coord::xy(500_000.0, 4_000_000.0),
                    wbvector::Coord::xy(500_100.0, 4_000_000.0),
                    wbvector::Coord::xy(500_100.0, 4_000_100.0),
                    wbvector::Coord::xy(500_000.0, 4_000_000.0),
                ],
                Vec::new(),
            )),
            &[],
        )
        .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let mask = wbvector::memory_store::make_vector_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": raster(1, 2, &[-20.0, -20.0]), "units": "db",
            "incidence_angle": 35.0, "wind_direction": 0.0, "look_direction": 0.0,
            "water_mask": mask,
        }))
        .unwrap();
        let err = ExtractOceanWindsTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("do not overlap"), "{err}");
    }

    #[test]
    fn nodata_cells_stay_nodata() {
        let (out, res) = run(json!({
            "input": raster(1, 2, &[-9999.0, -20.0]), "units": "db",
            "incidence_angle": 35.0, "wind_direction": 0.0, "look_direction": 0.0,
        }));
        assert_eq!(out.get(0, 0, 0), out.nodata);
        assert_eq!(res.outputs["valid_cells"], json!(1));
    }

    #[test]
    fn hh_polarization_reports_a_different_wind_than_vv() {
        // The ratio is an approximation, but it must actually be applied —
        // silently ignoring it would mislabel HH scenes as native VV.
        let s0 = cmod5n_forward(10.0, 0.0, 35.0);
        let db = 10.0 * s0.log10();
        let img = raster(1, 1, &[db]);
        let base = json!({
            "input": img, "units": "db", "incidence_angle": 35.0,
            "wind_direction": 0.0, "look_direction": 0.0,
        });
        let (vv, _) = run(base.clone());
        let mut hh_args = base;
        hh_args["polarization"] = json!("hh");
        let (hh, res) = run(hh_args);
        assert_eq!(res.outputs["polarization"], json!("hh"));
        assert!(
            (vv.get(0, 0, 0) - hh.get(0, 0, 0)).abs() > 0.5,
            "the polarization ratio was not applied"
        );
    }

    #[test]
    fn the_polarization_ratio_is_above_one_and_grows_with_incidence() {
        assert!(polarization_ratio(20.0, 0.6) > 1.0);
        assert!(polarization_ratio(45.0, 0.6) > polarization_ratio(20.0, 0.6));
    }

    #[test]
    fn unsolvable_cells_are_counted_and_left_nodata() {
        let (out, res) = run(json!({
            // +40 dB is far beyond anything 50 m/s produces.
            "input": raster(1, 1, &[40.0]), "units": "db",
            "incidence_angle": 35.0, "wind_direction": 0.0, "look_direction": 0.0,
        }));
        assert_eq!(res.outputs["unsolvable_cells"], json!(1));
        assert_eq!(res.outputs["valid_cells"], json!(0));
        assert_eq!(out.get(0, 0, 0), out.nodata);
    }

    #[test]
    fn a_constant_supplied_as_a_string_behaves_like_a_number() {
        let s0 = cmod5n_forward(9.0, 0.0, 35.0);
        let (out, _) = run(json!({
            "input": raster(1, 1, &[10.0 * s0.log10()]), "units": "db",
            "incidence_angle": "35", "wind_direction": "0", "look_direction": "0",
        }));
        assert!((out.get(0, 0, 0) - 9.0).abs() < 0.05);
    }

    #[test]
    fn rejects_bad_parameters() {
        let img = raster(1, 1, &[-20.0]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ExtractOceanWindsTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"input": img.clone()})));
        // Direction is a required prior, not an optional refinement.
        assert!(bad(json!({"input": img.clone(), "incidence_angle": 35.0})));
        assert!(bad(json!({
            "input": img.clone(), "incidence_angle": 35.0, "wind_direction": 0.0,
        })));
        let full = json!({
            "input": img, "incidence_angle": 35.0, "wind_direction": 0.0,
            "look_direction": 0.0,
        });
        let with = |k: &str, v: Value| {
            let mut m = full.clone();
            m[k] = v;
            m
        };
        assert!(bad(with("units", json!("watts"))));
        assert!(bad(with("polarization", json!("hv"))));
        assert!(bad(with("max_wind", json!(0.1))));
        assert!(bad(with("max_wind", json!(-5))));
        // Unbounded above, the per-cell scan would run ~2e9 evaluations and
        // the tool would stop responding rather than fail.
        assert!(bad(with("max_wind", json!(1e9))));
        // A non-string water_mask must not be silently dropped.
        assert!(bad(with("water_mask", json!(1))));
    }
}
