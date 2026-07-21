//! GeoLibre tool: incoming solar radiation over a DEM.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Raster Solar Radiation* (Spatial
//! Analyst) — a flagship capability (solar siting, agriculture, ecology) absent
//! from the bundled suite, whose building blocks (`horizon_angle`,
//! `time_in_daylight`, `shadow_image`) exist but are never integrated into
//! energy units.
//!
//! For every cell the tool computes total incoming shortwave radiation (Wh/m²)
//! over a date range by integrating over sun positions:
//!
//! 1. **Slope & aspect** from the DEM (Horn's 3×3 gradient).
//! 2. **Horizon angles** in 16 azimuth sectors (a capped terrain ray-scan) give
//!    self-/cast-shadowing and a **sky-view factor** for the diffuse term.
//! 3. For each sampled day and time step with the sun up: the direct-normal
//!    irradiance `I0·E0·τ^m` (air mass `m = 1/sin(alt)`) projected onto the slope
//!    via the incidence angle and zeroed where the terrain horizon blocks the
//!    sun, plus an isotropic diffuse term weighted by the sky-view factor.
//! 4. Instantaneous irradiance × the time step (and × the day sampling interval)
//!    accumulates to the period total; `direct_output` / `diffuse_output` split
//!    the components.
//!
//! Everything is deterministic — dates are explicit parameters, no `Date::now`.
//! Use a **projected DEM in metres** (slope/horizon need real distances);
//! `latitude` is taken from the parameter, or derived for EPSG:4326 / 3857.
//! Cost is roughly cells × sampled-days × steps × horizon rays, so this suits
//! moderate DEMs (or a coarse `time_step` / `day_interval`).

use std::collections::BTreeMap;
use std::f64::consts::PI;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};

/// Number of azimuth sectors for horizon / sky-view computation.
const N_AZIMUTH: usize = 16;
/// Solar constant (W/m²).
const SOLAR_CONSTANT: f64 = 1367.0;

pub struct SolarRadiationTool;

impl Tool for SolarRadiationTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "solar_radiation",
            display_name: "Solar Radiation",
            summary: "Incoming solar radiation (Wh/m²) over a DEM for a date range: direct + diffuse, accounting for slope/aspect, horizon shading, and sun position through time — like ArcGIS Raster Solar Radiation.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "dem",
                    description: "Input elevation raster (projected, metres — slope and horizon need real distances).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output total-radiation raster (Wh/m²). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "direct_output",
                    description: "Optional output raster for the direct-radiation component (Wh/m²).",
                    required: false,
                },
                ToolParamSpec {
                    name: "diffuse_output",
                    description: "Optional output raster for the diffuse-radiation component (Wh/m²).",
                    required: false,
                },
                ToolParamSpec {
                    name: "start_day",
                    description: "First day of year (1-365) or an ISO date (YYYY-MM-DD). Default 1.",
                    required: false,
                },
                ToolParamSpec {
                    name: "end_day",
                    description: "Last day of year (1-365) or an ISO date. Default 365.",
                    required: false,
                },
                ToolParamSpec {
                    name: "day_interval",
                    description: "Sample every this many days across the range (each sample represents that many days). Default: span/12 (about 12 samples).",
                    required: false,
                },
                ToolParamSpec {
                    name: "time_step",
                    description: "Hours between sun-position samples through each day. Default 0.5.",
                    required: false,
                },
                ToolParamSpec {
                    name: "latitude",
                    description: "Latitude in degrees. Default: derived from the DEM centre for EPSG:4326 / 3857.",
                    required: false,
                },
                ToolParamSpec {
                    name: "transmittivity",
                    description: "Atmospheric transmittivity (0-1) for a unit air mass. Default 0.6.",
                    required: false,
                },
                ToolParamSpec {
                    name: "diffuse_proportion",
                    description: "Fraction of extraterrestrial horizontal irradiance treated as diffuse sky radiation (0-1). Default 0.3.",
                    required: false,
                },
                ToolParamSpec {
                    name: "horizon_distance",
                    description: "Maximum horizon ray length in cells (bounds shadow computation). Default 100.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based DEM band (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args
            .get("dem")
            .and_then(Value::as_str)
            .map(str::trim)
            .unwrap_or("")
            .is_empty()
        {
            return Err(ToolError::Validation(
                "missing required string parameter 'dem'".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let dem_path = args
            .get("dem")
            .and_then(Value::as_str)
            .filter(|s| !s.trim().is_empty())
            .ok_or_else(|| ToolError::Validation("missing required parameter 'dem'".to_string()))?;
        let output = parse_optional_output(args, "output")?;
        let direct_output = parse_optional_output(args, "direct_output")?;
        let diffuse_output = parse_optional_output(args, "diffuse_output")?;
        let prm = parse_params(args)?;

        let dem = load_input_raster(dem_path)?;
        if prm.band as usize >= dem.bands {
            return Err(ToolError::Validation(format!(
                "band {} out of range (raster has {} band(s))",
                prm.band + 1,
                dem.bands
            )));
        }
        let rows = dem.rows;
        let cols = dem.cols;
        let nodata = dem.nodata;
        let cell = dem.cell_size_x.min(dem.cell_size_y).max(f64::MIN_POSITIVE);

        let latitude = resolve_latitude(&dem, prm.latitude)?.to_radians();

        // Elevation grid.
        let mut z = vec![f64::NAN; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let v = dem.get(prm.band, r as isize, c as isize);
                z[r * cols + c] = if v == nodata || v.is_nan() {
                    f64::NAN
                } else {
                    v
                };
            }
        }

        ctx.progress.info("computing slope, aspect, and horizons");
        let (slope, aspect) = slope_aspect(&z, rows, cols, cell);
        // Horizon angle (radians) per cell per azimuth sector, and sky-view.
        let horizons = compute_horizons(&z, rows, cols, cell, prm.horizon_distance);
        let svf = sky_view(&horizons, rows, cols, &slope);

        // Day / time sampling schedule.
        let (d0, d1) = (prm.start_day, prm.end_day);
        let span = (d1 - d0 + 1).max(1);
        let day_interval = prm.day_interval.unwrap_or((span / 12).max(1));
        let mut sample_days = Vec::new();
        let mut d = d0;
        while d <= d1 {
            sample_days.push(d);
            d += day_interval;
        }

        ctx.progress.info(&format!(
            "integrating {} sample day(s) x ~{:.0} step(s)/day",
            sample_days.len(),
            24.0 / prm.time_step
        ));

        let n = rows * cols;
        let mut total = vec![0.0f64; n];
        let mut direct = vec![0.0f64; n];
        let mut diffuse = vec![0.0f64; n];

        for &day in &sample_days {
            let decl = declination(day);
            let e0 = 1.0 + 0.033 * (2.0 * PI * day as f64 / 365.0).cos();
            let i_ext = SOLAR_CONSTANT * e0;
            let mut hour = 0.0;
            while hour < 24.0 {
                let (alt, az) = sun_position(latitude, decl, hour);
                if alt > 0.017 {
                    // > ~1°
                    let sin_alt = alt.sin();
                    let m = (1.0 / sin_alt).min(38.0);
                    let dni = i_ext * prm.transmittivity.powf(m);
                    let dhi = prm.diffuse_proportion * i_ext * sin_alt; // diffuse horizontal
                                                                        // Nearest azimuth sector for the shadow test.
                    let sector =
                        (((az / (2.0 * PI)) * N_AZIMUTH as f64).round() as usize) % N_AZIMUTH;
                    for idx in 0..n {
                        if z[idx].is_nan() {
                            continue;
                        }
                        let shadowed = alt < horizons[idx * N_AZIMUTH + sector];
                        let dir = if shadowed {
                            0.0
                        } else {
                            let cos_i = slope[idx].cos() * sin_alt
                                + slope[idx].sin() * alt.cos() * (az - aspect[idx]).cos();
                            dni * cos_i.max(0.0)
                        };
                        let dif = dhi * svf[idx];
                        let e = (dir + dif) * prm.time_step * day_interval as f64;
                        direct[idx] += dir * prm.time_step * day_interval as f64;
                        diffuse[idx] += dif * prm.time_step * day_interval as f64;
                        total[idx] += e;
                    }
                }
                hour += prm.time_step;
            }
        }

        // Mask nodata cells back out.
        for idx in 0..n {
            if z[idx].is_nan() {
                total[idx] = nodata;
                direct[idx] = nodata;
                diffuse[idx] = nodata;
            }
        }

        let out_raster = raster_like_with_data(&dem, total, nodata, DataType::F32)?;
        let out_path = write_or_store_output(out_raster, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        if let Some(p) = direct_output {
            let r = raster_like_with_data(&dem, direct, nodata, DataType::F32)?;
            outputs.insert(
                "direct_output".to_string(),
                json!(write_or_store_output(r, Some(p))?),
            );
        }
        if let Some(p) = diffuse_output {
            let r = raster_like_with_data(&dem, diffuse, nodata, DataType::F32)?;
            outputs.insert(
                "diffuse_output".to_string(),
                json!(write_or_store_output(r, Some(p))?),
            );
        }
        outputs.insert("sample_days".to_string(), json!(sample_days.len()));
        outputs.insert("latitude".to_string(), json!(latitude.to_degrees()));
        Ok(ToolRunResult { outputs })
    }
}

// ── Solar geometry ───────────────────────────────────────────────────────────

/// Solar declination (radians) for a day of year.
fn declination(day: i64) -> f64 {
    0.409_28 * (2.0 * PI * (284.0 + day as f64) / 365.0).sin()
}

/// Sun altitude and azimuth (radians; azimuth from north, clockwise) at a solar
/// hour, for a latitude and declination.
fn sun_position(lat: f64, decl: f64, hour: f64) -> (f64, f64) {
    let h = (hour - 12.0) * 15.0 * PI / 180.0; // hour angle
    let sin_alt = lat.sin() * decl.sin() + lat.cos() * decl.cos() * h.cos();
    let alt = sin_alt.clamp(-1.0, 1.0).asin();
    let cos_alt = alt.cos().max(1e-9);
    let cos_az = (decl.sin() * lat.cos() - decl.cos() * lat.sin() * h.cos()) / cos_alt;
    let sin_az = -decl.cos() * h.sin() / cos_alt;
    let mut az = sin_az.atan2(cos_az);
    if az < 0.0 {
        az += 2.0 * PI;
    }
    (alt, az)
}

// ── Terrain derivatives ──────────────────────────────────────────────────────

/// Slope and aspect (radians; aspect from north, clockwise) per cell via Horn's
/// 3×3 gradient. Flat cells get aspect 0.
fn slope_aspect(z: &[f64], rows: usize, cols: usize, cell: f64) -> (Vec<f64>, Vec<f64>) {
    let mut slope = vec![0.0; rows * cols];
    let mut aspect = vec![0.0; rows * cols];
    let at = |r: isize, c: isize| -> Option<f64> {
        if r < 0 || c < 0 || r >= rows as isize || c >= cols as isize {
            return None;
        }
        let v = z[r as usize * cols + c as usize];
        if v.is_nan() {
            None
        } else {
            Some(v)
        }
    };
    for r in 0..rows as isize {
        for c in 0..cols as isize {
            let idx = r as usize * cols + c as usize;
            let center = match at(r, c) {
                Some(v) => v,
                None => continue,
            };
            // Fill missing neighbours with the centre value.
            let g = |dr, dc| at(r + dr, c + dc).unwrap_or(center);
            let (nw, n, ne) = (g(-1, -1), g(-1, 0), g(-1, 1));
            let (w, e) = (g(0, -1), g(0, 1));
            let (sw, s, se) = (g(1, -1), g(1, 0), g(1, 1));
            // dz/dEast and dz/dNorth (row increases south, so north = r-1).
            let dz_de = ((ne + 2.0 * e + se) - (nw + 2.0 * w + sw)) / (8.0 * cell);
            let dz_dn = ((nw + 2.0 * n + ne) - (sw + 2.0 * s + se)) / (8.0 * cell);
            let grad = (dz_de * dz_de + dz_dn * dz_dn).sqrt();
            slope[idx] = grad.atan();
            if grad > 1e-12 {
                // Aspect = downhill direction (−gradient), from north clockwise.
                let mut a = (-dz_de).atan2(-dz_dn);
                if a < 0.0 {
                    a += 2.0 * PI;
                }
                aspect[idx] = a;
            }
        }
    }
    (slope, aspect)
}

/// Horizon angle (radians above horizontal) per cell for each of `N_AZIMUTH`
/// sectors, by scanning a ray up to `max_cells` cells.
fn compute_horizons(z: &[f64], rows: usize, cols: usize, cell: f64, max_cells: usize) -> Vec<f64> {
    let mut horizons = vec![0.0f64; rows * cols * N_AZIMUTH];
    for s in 0..N_AZIMUTH {
        let az = 2.0 * PI * s as f64 / N_AZIMUTH as f64; // from north, clockwise
                                                         // Step direction in grid: east = +col, north = -row.
        let de = az.sin(); // east component
        let dn = az.cos(); // north component
        let (dc, dr) = (de, -dn); // column step, row step
        for r in 0..rows {
            for c in 0..cols {
                let idx = r * cols + c;
                let z0 = z[idx];
                if z0.is_nan() {
                    continue;
                }
                let mut max_ang = 0.0f64;
                for k in 1..=max_cells {
                    let rr = r as f64 + dr * k as f64;
                    let cc = c as f64 + dc * k as f64;
                    let ri = rr.round() as isize;
                    let ci = cc.round() as isize;
                    if ri < 0 || ci < 0 || ri >= rows as isize || ci >= cols as isize {
                        break;
                    }
                    let zt = z[ri as usize * cols + ci as usize];
                    if zt.is_nan() {
                        continue;
                    }
                    let dist = (k as f64) * cell;
                    let ang = ((zt - z0) / dist).atan();
                    if ang > max_ang {
                        max_ang = ang;
                    }
                }
                horizons[idx * N_AZIMUTH + s] = max_ang;
            }
        }
    }
    horizons
}

/// Sky-view factor per cell: the sky fraction above the horizon, further reduced
/// by the surface tilt.
fn sky_view(horizons: &[f64], rows: usize, cols: usize, slope: &[f64]) -> Vec<f64> {
    let n = rows * cols;
    let mut svf = vec![0.0f64; n];
    for idx in 0..n {
        let mut acc = 0.0;
        for s in 0..N_AZIMUTH {
            let h = horizons[idx * N_AZIMUTH + s];
            // Fraction of the sky dome above horizon angle h in this sector.
            acc += h.cos().powi(2);
        }
        let terrain = acc / N_AZIMUTH as f64;
        // Tilted surface sees (1+cos β)/2 of the isotropic sky.
        svf[idx] = terrain * (1.0 + slope[idx].cos()) / 2.0;
    }
    svf
}

// ── Latitude ─────────────────────────────────────────────────────────────────

fn resolve_latitude(dem: &Raster, param: Option<f64>) -> Result<f64, ToolError> {
    if let Some(lat) = param {
        return Ok(lat);
    }
    let cy = dem.y_min + dem.rows as f64 * dem.cell_size_y * 0.5;
    match dem.crs.epsg {
        Some(4326) => Ok(cy),
        Some(3857) => {
            // Inverse spherical Mercator to latitude.
            let r = 6_378_137.0;
            Ok((cy / r).sinh().atan().to_degrees())
        }
        _ => Err(ToolError::Validation(
            "cannot derive latitude from this CRS; pass 'latitude' in degrees".to_string(),
        )),
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    start_day: i64,
    end_day: i64,
    day_interval: Option<i64>,
    time_step: f64,
    latitude: Option<f64>,
    transmittivity: f64,
    diffuse_proportion: f64,
    horizon_distance: usize,
    band: isize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let start_day = parse_day(args, "start_day")?.unwrap_or(1);
    let end_day = parse_day(args, "end_day")?.unwrap_or(365);
    if !(1..=366).contains(&start_day) || !(1..=366).contains(&end_day) || start_day > end_day {
        return Err(ToolError::Validation(
            "'start_day'/'end_day' must be day-of-year 1-365 with start <= end".to_string(),
        ));
    }
    let day_interval = match opt_f64(args, "day_interval")? {
        None => None,
        Some(v) if v >= 1.0 => Some(v as i64),
        Some(_) => {
            return Err(ToolError::Validation(
                "'day_interval' must be >= 1".to_string(),
            ))
        }
    };
    let time_step = opt_f64(args, "time_step")?.unwrap_or(0.5);
    if !(time_step > 0.0 && time_step <= 12.0) {
        return Err(ToolError::Validation(
            "'time_step' must be in (0, 12] hours".to_string(),
        ));
    }
    let latitude = opt_f64(args, "latitude")?;
    if let Some(l) = latitude {
        if !(-90.0..=90.0).contains(&l) {
            return Err(ToolError::Validation(
                "'latitude' must be between -90 and 90".to_string(),
            ));
        }
    }
    let transmittivity = opt_f64(args, "transmittivity")?.unwrap_or(0.6);
    if !(0.0..=1.0).contains(&transmittivity) {
        return Err(ToolError::Validation(
            "'transmittivity' must be between 0 and 1".to_string(),
        ));
    }
    let diffuse_proportion = opt_f64(args, "diffuse_proportion")?.unwrap_or(0.3);
    if !(0.0..=1.0).contains(&diffuse_proportion) {
        return Err(ToolError::Validation(
            "'diffuse_proportion' must be between 0 and 1".to_string(),
        ));
    }
    let horizon_distance = opt_f64(args, "horizon_distance")?.unwrap_or(100.0).max(1.0) as usize;
    let band_1based = opt_f64(args, "band")?.map(|v| v as i64).unwrap_or(1).max(1);
    Ok(Params {
        start_day,
        end_day,
        day_interval,
        time_step,
        latitude,
        transmittivity,
        diffuse_proportion,
        horizon_distance,
        band: (band_1based - 1) as isize,
    })
}

/// Parses a day parameter: a day-of-year integer or an ISO date (its DOY).
fn parse_day(args: &ToolArgs, key: &str) -> Result<Option<i64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_i64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => {
            let s = s.trim();
            if let Ok(v) = s.parse::<i64>() {
                return Ok(Some(v));
            }
            // ISO YYYY-MM-DD -> day of year.
            let y: i64 = s.get(0..4).and_then(|v| v.parse().ok()).ok_or_else(|| {
                ToolError::Validation(format!("'{key}' must be a day-of-year or ISO date"))
            })?;
            let m: i64 = s
                .get(5..7)
                .and_then(|v| v.parse().ok())
                .ok_or_else(|| ToolError::Validation(format!("'{key}' has a bad month")))?;
            let d: i64 = s
                .get(8..10)
                .and_then(|v| v.parse().ok())
                .ok_or_else(|| ToolError::Validation(format!("'{key}' has a bad day")))?;
            Ok(Some(day_of_year(y, m, d)))
        }
        Some(_) => Err(ToolError::Validation(format!(
            "'{key}' must be a day-of-year or ISO date"
        ))),
    }
}

fn day_of_year(y: i64, m: i64, d: i64) -> i64 {
    let leap = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
    let days = [
        31,
        if leap { 29 } else { 28 },
        31,
        30,
        31,
        30,
        31,
        31,
        30,
        31,
        30,
        31,
    ];
    let mut doy = d;
    let months = (m - 1).clamp(0, 11) as usize;
    doy += days.iter().take(months).sum::<i64>();
    doy
}

fn opt_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be a number"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a number"
        ))),
    }
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

    /// A tilted DEM: elevation increases with a chosen gradient so the whole
    /// surface faces a known aspect.
    fn tilted_dem(cols: usize, rows: usize, dz_dcol: f64, dz_drow: f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 30.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(
                    0,
                    row as isize,
                    col as isize,
                    1000.0 + dz_dcol * col as f64 + dz_drow * row as f64,
                )
                .unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SolarRadiationTool.run(&args, &ctx()).unwrap();
        let r = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, r)
    }

    fn mean(r: &Raster) -> f64 {
        let mut s = 0.0;
        let mut n = 0;
        for row in 0..r.rows {
            for col in 0..r.cols {
                let v = r.get(0, row as isize, col as isize);
                if v != r.nodata && !v.is_nan() {
                    s += v;
                    n += 1;
                }
            }
        }
        s / n as f64
    }

    /// In the northern hemisphere, a south-facing slope receives much more solar
    /// energy than a north-facing one.
    #[test]
    fn south_slope_beats_north_slope_north_hemisphere() {
        // South-facing: elevation increases toward the north (row 0). Row 0 is
        // north, so higher north -> slope faces south. dz/drow < 0 means value
        // decreases with row (south), i.e. north is higher -> faces south.
        // Steep slopes on the winter solstice, when a low sun makes aspect
        // dominate (a north-facing slope gets almost no direct sun).
        let south = tilted_dem(20, 20, 0.0, -30.0); // north higher -> south-facing
        let north = tilted_dem(20, 20, 0.0, 30.0); // south higher -> north-facing
        let common = json!({
            "output": null, "start_day": 355, "end_day": 355, "time_step": 1.0,
            "latitude": 45.0, "horizon_distance": 5,
        });
        let (_o1, rs) = run({
            let mut m = common.clone();
            m["dem"] = json!(south);
            m
        });
        let (_o2, rn) = run({
            let mut m = common.clone();
            m["dem"] = json!(north);
            m
        });
        let (es, en) = (mean(&rs), mean(&rn));
        assert!(
            es > en * 1.5,
            "winter south slope {es} should far beat north slope {en}"
        );
    }

    /// A flat surface at mid-latitude on the summer solstice gets a plausible
    /// daily total (a few kWh/m²).
    #[test]
    fn flat_surface_plausible_magnitude() {
        let flat = tilted_dem(10, 10, 0.0, 0.0);
        let (_o, r) = run(json!({
            "dem": flat, "start_day": 172, "end_day": 172, "time_step": 0.5,
            "latitude": 40.0, "horizon_distance": 3,
        }));
        let e = mean(&r); // Wh/m² for one day
        assert!(
            (2000.0..12000.0).contains(&e),
            "flat summer daily total {e} Wh/m² is outside the plausible 2-12 kWh range"
        );
    }

    #[test]
    fn sun_is_up_at_noon_summer() {
        let (alt, _az) = sun_position(45.0_f64.to_radians(), declination(172), 12.0);
        // Noon summer solstice at 45N: altitude ~ 68°.
        assert!(
            (alt.to_degrees() - 68.0).abs() < 3.0,
            "noon altitude {} off",
            alt.to_degrees()
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SolarRadiationTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "dem": "d.tif", "transmittivity": 2.0 })).is_err());
        assert!(bad(json!({ "dem": "d.tif", "start_day": 300, "end_day": 100 })).is_err());
        assert!(bad(json!({ "dem": "d.tif", "latitude": 45.0 })).is_ok());
    }
}
