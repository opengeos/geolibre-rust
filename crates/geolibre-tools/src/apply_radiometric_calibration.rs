//! GeoLibre tool: calibrate SAR digital numbers to physical backscatter.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Apply Radiometric Calibration*
//! (Image Analyst).
//!
//! ## Why the catalog needs it
//!
//! Raw SAR digital numbers are sensor- and scene-specific. Two acquisitions of
//! the same field, from the same satellite, can differ by an order of magnitude
//! in DN with no change on the ground. Every published backscatter threshold
//! (flood mapping, biomass, soil moisture) is stated in sigma nought or gamma
//! nought, so applying one to uncalibrated DNs is meaningless.
//!
//! This is the *first* step of a real SAR workflow, and it is the missing front
//! end of the chain the catalog already has: `multilook`, `sar_coherence`, the
//! Lee speckle filters, and `compute_sar_indices` all consume whatever values
//! they are handed. `rescale_value_range` is a linear stretch with no
//! radiometric meaning and cannot substitute.
//!
//! ## The three conventions
//!
//! ```text
//! beta0  = DN^2 / A^2          (radar brightness, geometry-independent)
//! sigma0 = beta0 * sin(theta)  (per unit ground area)
//! gamma0 = beta0 * tan(theta)  (per unit area perpendicular to look direction)
//! ```
//!
//! where `A` is the sensor calibration constant (scalar or a per-range LUT) and
//! `theta` the local incidence angle. Sigma nought and gamma nought therefore
//! *require* an incidence angle; asking for them without one is an error rather
//! than a silent fallback to beta nought, because the difference is a factor of
//! two or more over a typical swath.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::args_common::{band_index, choice_or, opt_f64, req_str};
use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};
use crate::raster_stack::check_alignment;

pub struct ApplyRadiometricCalibrationTool;

impl Tool for ApplyRadiometricCalibrationTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "apply_radiometric_calibration",
            display_name: "Apply Radiometric Calibration",
            summary: "Converts raw SAR digital numbers into calibrated backscatter — beta nought, sigma nought or gamma nought — so acquisitions become physically comparable and published thresholds apply (ArcGIS Apply Radiometric Calibration). Nothing in either registry does this: rescale_value_range is a linear stretch with no radiometric meaning, and multilook, sar_coherence, the Lee speckle filters and compute_sar_indices all consume whatever values they are handed.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "SAR raster in digital numbers, amplitude, intensity or dB (see 'input_units').",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output calibrated backscatter raster. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "calibration_type",
                    description: "'sigma0' (default), 'beta0', or 'gamma0'. sigma0 and gamma0 require 'incidence_angle'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "calibration_constant",
                    description: "Sensor calibration constant A (default 1.0). Ignored when 'calibration_lut' is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "calibration_lut",
                    description: "Optional co-registered raster of per-cell calibration constants A, overriding 'calibration_constant'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "incidence_angle",
                    description: "Local incidence angle in degrees: a constant, or a co-registered raster path. Required for 'sigma0' and 'gamma0'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "input_units",
                    description: "'dn' (default), 'amplitude', 'intensity', or 'db'. dn and amplitude are squared to power; db is delinearised first.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_units",
                    description: "'linear' (default) or 'db'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to calibrate (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        let kind = parse_calibration(args)?;
        choice_or(
            args,
            "input_units",
            &["dn", "amplitude", "intensity", "db"],
            "dn",
        )?;
        choice_or(args, "output_units", &["linear", "db"], "linear")?;
        band_index(args, "band")?;
        if let Some(a) = opt_f64(args, "calibration_constant")? {
            if a <= 0.0 {
                return Err(ToolError::Validation(
                    "'calibration_constant' must be > 0".to_string(),
                ));
            }
        }
        // Refuse rather than silently degrading to beta0: over a Sentinel-1
        // swath the incidence-angle factor spans more than 2x.
        if kind.needs_incidence() && args.get("incidence_angle").is_none() {
            return Err(ToolError::Validation(format!(
                "'{}' requires 'incidence_angle' (a constant in degrees or a raster path)",
                kind.label()
            )));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?.to_string();
        let kind = parse_calibration(args)?;
        let input_units = choice_or(
            args,
            "input_units",
            &["dn", "amplitude", "intensity", "db"],
            "dn",
        )?;
        let db_output = choice_or(args, "output_units", &["linear", "db"], "linear")? == "db";
        let band = band_index(args, "band")?;
        let output = parse_optional_output(args, "output")?;

        let raster = load_input_raster(&input)?;
        let (rows, cols) = (raster.rows, raster.cols);
        if band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "'band' {} is out of range; '{input}' has {} band(s)",
                band + 1,
                raster.bands
            )));
        }

        let cal = parse_cell_source(args, "calibration_lut", "calibration_constant", &raster)?
            .unwrap_or(CellSource::Constant(1.0));
        let incidence = parse_cell_source(args, "incidence_angle", "incidence_angle", &raster)?;
        if kind.needs_incidence() && incidence.is_none() {
            return Err(ToolError::Validation(format!(
                "'{}' requires 'incidence_angle'",
                kind.label()
            )));
        }

        ctx.progress.info(&format!(
            "{rows}x{cols}, {} from {input_units}, output {}",
            kind.label(),
            if db_output { "dB" } else { "linear" }
        ));

        let nodata = -9999.0_f64;
        let mut out = vec![nodata; rows * cols];
        let mut dropped = 0usize;

        for r in 0..rows {
            for c in 0..cols {
                let raw = raster.get(band, r as isize, c as isize);
                if raw == raster.nodata || !raw.is_finite() {
                    continue;
                }
                // Everything downstream is a power, so normalise first.
                let Some(power) = to_power(raw, input_units) else {
                    dropped += 1;
                    continue;
                };
                let Some(a) = cal.at(r, c) else {
                    continue;
                };
                if a <= 0.0 {
                    dropped += 1;
                    continue;
                }
                let beta0 = power / (a * a);

                let value = match kind {
                    Calibration::Beta0 => Some(beta0),
                    Calibration::Sigma0 | Calibration::Gamma0 => {
                        let Some(theta_deg) = incidence.as_ref().and_then(|s| s.at(r, c)) else {
                            continue;
                        };
                        // Outside (0, 90) the trigonometric factors are
                        // meaningless (tan blows up at 90 exactly).
                        if theta_deg <= 0.0 || theta_deg >= 90.0 {
                            dropped += 1;
                            continue;
                        }
                        let theta = theta_deg.to_radians();
                        Some(match kind {
                            Calibration::Sigma0 => beta0 * theta.sin(),
                            _ => beta0 * theta.tan(),
                        })
                    }
                };

                let Some(v) = value else { continue };
                if db_output {
                    if v > 0.0 {
                        out[r * cols + c] = 10.0 * v.log10();
                    } else {
                        dropped += 1;
                    }
                } else {
                    out[r * cols + c] = v;
                }
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let out_r = raster_like_with_data(&raster, out, nodata, DataType::F32)?;
        let out_path = write_or_store_output(out_r, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("calibration_type".to_string(), json!(kind.label()));
        outputs.insert("input_units".to_string(), json!(input_units));
        outputs.insert(
            "output_units".to_string(),
            json!(if db_output { "db" } else { "linear" }),
        );
        outputs.insert("out_of_domain_cells".to_string(), json!(dropped));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        Ok(ToolRunResult { outputs })
    }
}

/// Normalises an input value to linear power.
fn to_power(v: f64, units: &str) -> Option<f64> {
    match units {
        // DN and amplitude are both field quantities: power is their square.
        "dn" | "amplitude" => Some(v * v),
        "intensity" => (v >= 0.0).then_some(v),
        _ => Some(10.0_f64.powf(v / 10.0)),
    }
}

// ── Parameters ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
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

    fn needs_incidence(self) -> bool {
        !matches!(self, Calibration::Beta0)
    }
}

fn parse_calibration(args: &ToolArgs) -> Result<Calibration, ToolError> {
    Ok(
        match choice_or(
            args,
            "calibration_type",
            &["sigma0", "beta0", "gamma0"],
            "sigma0",
        )? {
            "beta0" => Calibration::Beta0,
            "gamma0" => Calibration::Gamma0,
            _ => Calibration::Sigma0,
        },
    )
}

/// A per-cell quantity supplied either as a constant or as a raster.
enum CellSource {
    Constant(f64),
    Raster(Box<Raster>),
}

impl CellSource {
    fn at(&self, row: usize, col: usize) -> Option<f64> {
        match self {
            CellSource::Constant(v) => Some(*v),
            CellSource::Raster(r) => {
                let v = r.get(0, row as isize, col as isize);
                (v != r.nodata && v.is_finite()).then_some(v)
            }
        }
    }
}

/// Resolves a dual-typed parameter into a per-cell source.
///
/// `raster_key` and `constant_key` may name the same parameter (as they do for
/// `incidence_angle`, where `39` and a raster path are both legal), so a value
/// only counts as a path when it is a **non-numeric** string. That keeps the
/// resolution order-independent, which matters because when the keys *differ*
/// (`calibration_lut` vs `calibration_constant`) the raster must win — a caller
/// supplying a LUT alongside a leftover constant expects the LUT.
fn parse_cell_source(
    args: &ToolArgs,
    raster_key: &str,
    constant_key: &str,
    template: &Raster,
) -> Result<Option<CellSource>, ToolError> {
    let path = args
        .get(raster_key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty() && s.parse::<f64>().is_err());
    if let Some(path) = path {
        let raster = load_input_raster(path)?;
        // Read cell-by-cell against the input grid, so a mismatched grid would
        // silently calibrate against the wrong locations.
        check_alignment(&[template.clone(), raster.clone()])?;
        return Ok(Some(CellSource::Raster(Box::new(raster))));
    }
    Ok(opt_f64(args, constant_key)?.map(CellSource::Constant))
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

    /// Output rasters are F32, so value comparisons use a **relative**
    /// tolerance rather than an absolute one.
    fn close(actual: f64, expect: f64) -> bool {
        (actual - expect).abs() <= 1e-6 * expect.abs().max(1.0)
    }

    fn raster(data: &[f64]) -> String {
        let cols = data.len();
        let mut r = Raster::new(RasterConfig {
            cols,
            rows: 1,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
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
        for (c, v) in data.iter().enumerate() {
            r.set(0, 0, c as isize, *v).unwrap();
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: Value) -> (Raster, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = ApplyRadiometricCalibrationTool.run(&args, &ctx()).unwrap();
        let out = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        (out, res)
    }

    #[test]
    fn beta0_squares_the_dn_and_divides_by_the_constant_squared() {
        // DN 10, A = 2 gives 100 / 4 = 25.
        let (out, _) = run(json!({
            "input": raster(&[10.0]),
            "calibration_type": "beta0",
            "calibration_constant": 2.0,
        }));
        assert!(close(out.get(0, 0, 0), 25.0));
    }

    #[test]
    fn sigma0_and_gamma0_apply_the_incidence_factors() {
        let src = raster(&[10.0]);
        let beta0 = 100.0;
        let theta = 30.0_f64;
        let (sig, _) = run(json!({
            "input": src, "calibration_type": "sigma0", "incidence_angle": theta,
        }));
        assert!(close(sig.get(0, 0, 0), beta0 * theta.to_radians().sin()));
        let (gam, _) = run(json!({
            "input": src, "calibration_type": "gamma0", "incidence_angle": theta,
        }));
        assert!(close(gam.get(0, 0, 0), beta0 * theta.to_radians().tan()));
        // The conventions must actually differ, or the checks are vacuous.
        assert!((sig.get(0, 0, 0) - gam.get(0, 0, 0)).abs() > 1.0);
    }

    #[test]
    fn sigma0_without_an_incidence_angle_is_rejected_not_silently_beta0() {
        let args: ToolArgs = serde_json::from_value(json!({
            "input": raster(&[10.0]), "calibration_type": "sigma0",
        }))
        .unwrap();
        assert!(ApplyRadiometricCalibrationTool.validate(&args).is_err());
    }

    #[test]
    fn incidence_angle_may_vary_across_the_swath() {
        // The reason a constant is not enough: the factor spans the swath.
        let (out, _) = run(json!({
            "input": raster(&[10.0, 10.0]),
            "calibration_type": "sigma0",
            "incidence_angle": raster(&[20.0, 45.0]),
        }));
        assert!(close(out.get(0, 0, 0), 100.0 * 20.0_f64.to_radians().sin()));
        assert!(close(out.get(0, 0, 1), 100.0 * 45.0_f64.to_radians().sin()));
    }

    #[test]
    fn a_calibration_lut_overrides_the_constant() {
        let (out, _) = run(json!({
            "input": raster(&[10.0, 10.0]),
            "calibration_type": "beta0",
            "calibration_constant": 2.0,
            "calibration_lut": raster(&[1.0, 5.0]),
        }));
        assert!(close(out.get(0, 0, 0), 100.0));
        assert!(close(out.get(0, 0, 1), 4.0));
    }

    #[test]
    fn input_units_change_the_power_normalisation() {
        // Intensity 100, amplitude 10 and 20 dB all describe the same power.
        let (from_amp, _) = run(json!({
            "input": raster(&[10.0]), "calibration_type": "beta0", "input_units": "amplitude",
        }));
        let (from_int, _) = run(json!({
            "input": raster(&[100.0]), "calibration_type": "beta0", "input_units": "intensity",
        }));
        assert!(close(from_amp.get(0, 0, 0), from_int.get(0, 0, 0)));
        let (from_db, _) = run(json!({
            "input": raster(&[20.0]), "calibration_type": "beta0", "input_units": "db",
        }));
        assert!(close(from_db.get(0, 0, 0), from_int.get(0, 0, 0)));
    }

    #[test]
    fn db_output_is_the_log_of_the_linear_answer() {
        let src = raster(&[10.0]);
        let (lin, _) = run(json!({"input": src, "calibration_type": "beta0"}));
        let (db, _) = run(json!({
            "input": src, "calibration_type": "beta0", "output_units": "db",
        }));
        assert!(close(db.get(0, 0, 0), 10.0 * lin.get(0, 0, 0).log10()));
    }

    #[test]
    fn zero_power_is_nodata_in_db_output_not_negative_infinity() {
        let (out, res) = run(json!({
            "input": raster(&[0.0]), "calibration_type": "beta0", "output_units": "db",
        }));
        assert_eq!(out.get(0, 0, 0), out.nodata);
        assert_eq!(res.outputs["out_of_domain_cells"], json!(1));
    }

    #[test]
    fn out_of_range_incidence_angles_are_dropped() {
        let (out, res) = run(json!({
            "input": raster(&[10.0, 10.0]),
            "calibration_type": "gamma0",
            "incidence_angle": raster(&[0.0, 90.0]),
        }));
        assert_eq!(out.get(0, 0, 0), out.nodata);
        assert_eq!(out.get(0, 0, 1), out.nodata);
        assert_eq!(res.outputs["out_of_domain_cells"], json!(2));
    }

    #[test]
    fn nodata_passes_through() {
        let (out, _) = run(json!({
            "input": raster(&[-9999.0, 10.0]), "calibration_type": "beta0",
        }));
        assert_eq!(out.get(0, 0, 0), out.nodata);
        assert!(close(out.get(0, 0, 1), 100.0));
    }

    #[test]
    fn rejects_bad_parameters() {
        let src = raster(&[1.0]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ApplyRadiometricCalibrationTool.validate(&args).is_err()
        };
        assert!(bad(json!({"input": src, "calibration_type": "nope"})));
        assert!(bad(
            json!({"input": src, "calibration_type": "beta0", "calibration_constant": 0})
        ));
        assert!(bad(
            json!({"input": src, "calibration_type": "beta0", "input_units": "watts"})
        ));
        assert!(bad(
            json!({"input": src, "calibration_type": "beta0", "band": 0})
        ));
    }
}
