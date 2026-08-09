//! GeoLibre tool: subtract the sensor's thermal noise floor from SAR
//! backscatter.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Remove Thermal Noise* (Image
//! Analyst).
//!
//! ## Why the catalog needs it
//!
//! Every SAR receiver contributes additive thermal noise. Over bright targets
//! it is irrelevant; over the dark surfaces that SAR is most often used to
//! study — calm water, flooded ground, dry sand, radar shadow — it is a large
//! fraction of the measured power, and in Sentinel-1 IW it varies by several dB
//! *across the swath*, leaving the characteristic bright/dark banding at
//! sub-swath boundaries.
//!
//! That banding breaks exactly the tools the catalog already ships: an Otsu
//! threshold for water extraction latches onto a sub-swath edge instead of the
//! shoreline, and `compute_sar_indices` computes cross-polarisation ratios in
//! which the numerator is mostly noise. Neither registry has anything for it —
//! `lee_filter` and its relatives suppress *multiplicative* speckle, which is a
//! different quantity requiring a different (and non-commuting) correction.
//!
//! ## Noise sources
//!
//! The real noise-equivalent sigma nought lives in the product's own metadata
//! (Sentinel-1 ships it as a per-swath LUT in `noise-*.xml`), which this
//! catalog is never handed — tools here receive a raster, not a SAFE archive.
//! Three explicit forms are accepted instead, in order of precedence:
//!
//! * `noise_raster` — a full per-cell noise-power raster, the exact case;
//! * `noise_profile` — a comma-separated across-range profile, linearly
//!   interpolated across the columns, which is the shape a Sentinel-1 LUT
//!   actually has;
//! * `noise_constant` — a single noise-equivalent sigma nought.
//!
//! ## Clamping
//!
//! Subtraction in the power domain can drive a cell to zero or below wherever
//! the signal is at or under the noise floor. Those cells are clamped to
//! `floor_fraction` of the local noise power rather than to zero, so a
//! subsequent dB conversion stays finite instead of producing `-inf`; the count
//! of clamped cells is reported, because a large fraction means the scene is
//! noise-limited and the result should not be trusted.

use std::collections::BTreeMap;

use serde_json::{json, Value};
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

pub struct RemoveThermalNoiseTool;

impl Tool for RemoveThermalNoiseTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "remove_thermal_noise",
            display_name: "Remove Thermal Noise",
            summary: "Subtracts the receiver's additive thermal noise floor from SAR backscatter, removing the across-swath banding that dominates dark surfaces such as calm water, flooded ground and radar shadow (ArcGIS Remove Thermal Noise). Nothing in either registry does this: the bundled lee/frost/kuan/refined-lee filters suppress multiplicative speckle, a different quantity. Without it an Otsu water threshold latches onto a sub-swath boundary rather than a shoreline, and cross-polarisation ratios are computed from a numerator that is largely noise.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Calibrated SAR backscatter raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output denoised raster. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_snr",
                    description: "Output signal-to-noise ratio in dB. Always produced; stored in memory when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "noise_raster",
                    description: "Per-cell noise power raster, co-registered with the input. Takes precedence over the other two forms.",
                    required: false,
                },
                ToolParamSpec {
                    name: "noise_profile",
                    description: "Comma-separated across-range noise powers, linearly interpolated across the raster's columns. At least two values.",
                    required: false,
                },
                ToolParamSpec {
                    name: "noise_constant",
                    description: "Single noise-equivalent sigma nought applied to every cell.",
                    required: false,
                },
                ToolParamSpec {
                    name: "noise_units",
                    description: "Units the noise values are in: 'intensity' (default), 'dn', 'amplitude', or 'db'.",
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
                    name: "floor_fraction",
                    description: "Cells driven to or below zero are clamped to this fraction of the local noise power (default 0.01). Keeps a dB conversion finite.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "Band to denoise (default 0).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        let prm = parse_params(args)?;
        // One noise source must be given; silently defaulting to "no noise"
        // would make the tool a no-op that looks like it worked.
        if args.get("noise_raster").and_then(Value::as_str).is_none()
            && prm.profile.is_none()
            && prm.constant.is_none()
        {
            return Err(ToolError::Validation(
                "one of 'noise_raster', 'noise_profile' or 'noise_constant' is required"
                    .to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        self.validate(args)?;
        let input_path = req_str(args, "input")?.to_string();
        let prm = parse_params(args)?;
        let band = band_index(args, "band")?;
        let output = parse_optional_output(args, "output")?;
        let out_snr = parse_optional_output(args, "output_snr")?;

        let raster = load_input_raster(&input_path)?;
        let (rows, cols) = (raster.rows, raster.cols);

        let noise = noise_field(args, &prm, &raster, rows, cols)?;

        ctx.progress.info(&format!(
            "{rows}x{cols}, noise from {}, input {}",
            noise.1,
            prm.input_units.label()
        ));

        let nodata = -9999.0_f64;
        let mut out = vec![nodata; rows * cols];
        let mut snr = vec![nodata; rows * cols];
        let (mut clamped, mut valid) = (0usize, 0usize);

        for r in 0..rows {
            for c in 0..cols {
                let i = r * cols + c;
                let raw = raster.get(band, r as isize, c as isize);
                if raw == raster.nodata || !raw.is_finite() {
                    continue;
                }
                let Some(power) = prm.input_units.to_power(raw) else {
                    continue;
                };
                let n = noise.0[i];
                if !n.is_finite() || n < 0.0 {
                    continue;
                }

                let mut denoised = power - n;
                if denoised <= 0.0 {
                    // The measurement is at or below the noise floor: no signal
                    // was recovered. Clamp rather than emit zero so a dB
                    // conversion downstream stays finite.
                    denoised = (n * prm.floor_fraction).max(f64::MIN_POSITIVE);
                    clamped += 1;
                }
                valid += 1;

                out[i] = if prm.db_output {
                    match power_to_db(denoised) {
                        Some(v) => v,
                        None => continue,
                    }
                } else {
                    denoised
                };
                if n > 0.0 {
                    if let Some(v) = power_to_db(denoised / n) {
                        snr[i] = v;
                    }
                }
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let clamped_fraction = if valid == 0 {
            0.0
        } else {
            clamped as f64 / valid as f64
        };
        ctx.progress.info(&format!(
            "{valid} valid cell(s), {clamped} clamped at the noise floor ({:.1}%)",
            100.0 * clamped_fraction
        ));

        let out_path = write_or_store_output(
            raster_like_with_data(&raster, out, nodata, DataType::F32)?,
            output,
        )?;
        let snr_path = write_or_store_output(
            raster_like_with_data(&raster, snr, nodata, DataType::F32)?,
            out_snr,
        )?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("output_snr".to_string(), json!(snr_path));
        outputs.insert("noise_source".to_string(), json!(noise.1));
        outputs.insert("valid_cells".to_string(), json!(valid));
        outputs.insert("clamped_cells".to_string(), json!(clamped));
        outputs.insert("clamped_fraction".to_string(), json!(clamped_fraction));
        outputs.insert(
            "output_units".to_string(),
            json!(if prm.db_output { "db" } else { "linear" }),
        );
        Ok(ToolRunResult { outputs })
    }
}

/// Builds the per-cell noise power field, and names the source it came from.
fn noise_field(
    args: &ToolArgs,
    prm: &Params,
    template: &Raster,
    rows: usize,
    cols: usize,
) -> Result<(Vec<f64>, &'static str), ToolError> {
    if let Some(path) = args.get("noise_raster").and_then(Value::as_str) {
        let path = path.trim();
        if !path.is_empty() {
            let nr = load_input_raster(path)?;
            check_alignment_refs(&[template, &nr])?;
            let mut out = vec![f64::NAN; rows * cols];
            for r in 0..rows {
                for c in 0..cols {
                    let v = nr.get(0, r as isize, c as isize);
                    if v != nr.nodata && v.is_finite() {
                        // A noise raster is in the same representation as the
                        // profile/constant forms, so it goes through the same
                        // unit conversion.
                        if let Some(p) = prm.noise_units.to_power(v) {
                            out[r * cols + c] = p;
                        }
                    }
                }
            }
            return Ok((out, "raster"));
        }
    }

    if let Some(profile) = &prm.profile {
        let mut out = vec![0.0f64; rows * cols];
        for c in 0..cols {
            // Map column centres onto the profile's own sample positions, so a
            // profile of any length stretches across the full swath.
            let t = if cols == 1 {
                0.0
            } else {
                c as f64 / (cols - 1) as f64
            };
            let pos = t * (profile.len() - 1) as f64;
            let lo = pos.floor() as usize;
            let hi = (lo + 1).min(profile.len() - 1);
            let frac = pos - lo as f64;
            let v = profile[lo] * (1.0 - frac) + profile[hi] * frac;
            for r in 0..rows {
                out[r * cols + c] = v;
            }
        }
        return Ok((out, "profile"));
    }

    let k = prm.constant.expect("validate guarantees a noise source");
    Ok((vec![k; rows * cols], "constant"))
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    input_units: SarUnits,
    noise_units: SarUnits,
    db_output: bool,
    floor_fraction: f64,
    /// Noise profile already converted to power.
    profile: Option<Vec<f64>>,
    /// Noise constant already converted to power.
    constant: Option<f64>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let input_units = SarUnits::parse(args.get("input_units").and_then(Value::as_str).unwrap_or(""))?;
    let noise_units = SarUnits::parse(args.get("noise_units").and_then(Value::as_str).unwrap_or(""))?;
    let db_output = choice_or(args, "output_units", &["linear", "db"], "linear")? == "db";

    let floor_fraction = f64_or(args, "floor_fraction", 0.01)?;
    if !(0.0..1.0).contains(&floor_fraction) {
        return Err(ToolError::Validation(format!(
            "'floor_fraction' must be in [0, 1), got {floor_fraction}"
        )));
    }

    let profile = match args.get("noise_profile").and_then(Value::as_str) {
        None => None,
        Some(s) if s.trim().is_empty() => None,
        Some(s) => {
            let vals: Result<Vec<f64>, _> = s
                .split(',')
                .map(str::trim)
                .filter(|t| !t.is_empty())
                .map(|t| {
                    t.parse::<f64>().map_err(|_| {
                        ToolError::Validation(format!(
                            "'noise_profile' entry '{t}' is not a number"
                        ))
                    })
                })
                .collect();
            let vals = vals?;
            if vals.len() < 2 {
                return Err(ToolError::Validation(
                    "'noise_profile' needs at least two values".to_string(),
                ));
            }
            let mut powers = Vec::with_capacity(vals.len());
            for v in vals {
                powers.push(noise_units.to_power(v).ok_or_else(|| {
                    ToolError::Validation(format!("'noise_profile' value {v} is not a valid power"))
                })?);
            }
            Some(powers)
        }
    };

    let constant = match opt_f64(args, "noise_constant")? {
        None => None,
        Some(v) => Some(noise_units.to_power(v).ok_or_else(|| {
            ToolError::Validation(format!("'noise_constant' {v} is not a valid power"))
        })?),
    };

    Ok(Params {
        input_units,
        noise_units,
        db_output,
        floor_fraction,
        profile,
        constant,
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

    fn raster_of(cols: usize, rows: usize, vals: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
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
        let out = RemoveThermalNoiseTool.run(&args, &ctx()).unwrap();
        let r = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (r, out.outputs)
    }

    /// A constant noise floor is subtracted exactly in the power domain.
    #[test]
    fn subtracts_a_constant_floor() {
        let src = raster_of(3, 1, &[0.5, 0.2, 0.05]);
        let (out, outputs) = run(json!({
            "input": src, "noise_constant": 0.02
        }));
        assert!((out.get(0, 0, 0) - 0.48).abs() < 1e-6);
        assert!((out.get(0, 0, 1) - 0.18).abs() < 1e-6);
        assert!((out.get(0, 0, 2) - 0.03).abs() < 1e-6);
        assert_eq!(outputs["clamped_cells"].as_u64().unwrap(), 0);
        assert_eq!(outputs["noise_source"].as_str().unwrap(), "constant");
    }

    /// This is the reason the tool exists: an across-range noise ramp makes a
    /// uniform scene look like it has a brightness gradient, and removing it
    /// must flatten the scene back out.
    #[test]
    fn removes_across_swath_banding() {
        // True backscatter is a constant 0.30 everywhere. The receiver adds a
        // noise floor rising from 0.02 on the near edge to 0.10 on the far edge,
        // so the *measured* power ramps from 0.32 to 0.40 — a 1 dB tilt across
        // the swath that has nothing to do with the ground.
        let cols = 9;
        let measured: Vec<f64> = (0..cols)
            .map(|c| 0.30 + 0.02 + 0.08 * c as f64 / (cols - 1) as f64)
            .collect();
        let before_spread = measured.last().unwrap() - measured.first().unwrap();
        assert!(before_spread > 0.07, "fixture has no banding to remove");

        let (out, _) = run(json!({
            "input": raster_of(cols, 1, &measured),
            "noise_profile": "0.02,0.10"
        }));
        for c in 0..cols {
            let v = out.get(0, 0, c as isize);
            assert!(
                (v - 0.30).abs() < 1e-6,
                "column {c} still banded: {v} != 0.30"
            );
        }
    }

    /// The profile stretches across the columns however many samples it has.
    #[test]
    fn profile_interpolates_across_columns() {
        let cols = 5;
        // Measured = signal 1.0 + noise ramp 0.0 -> 0.4.
        let measured: Vec<f64> = (0..cols)
            .map(|c| 1.0 + 0.4 * c as f64 / (cols - 1) as f64)
            .collect();
        let (out, outputs) = run(json!({
            "input": raster_of(cols, 1, &measured),
            "noise_profile": "0.0, 0.2, 0.4"
        }));
        assert_eq!(outputs["noise_source"].as_str().unwrap(), "profile");
        for c in 0..cols {
            assert!(
                (out.get(0, 0, c as isize) - 1.0).abs() < 1e-6,
                "column {c} not flattened"
            );
        }
    }

    /// Cells at or below the noise floor are clamped, not zeroed, so dB stays
    /// finite — and the caller is told how many.
    #[test]
    fn clamps_cells_below_the_floor() {
        let src = raster_of(3, 1, &[0.5, 0.02, 0.001]);
        let (out, outputs) = run(json!({
            "input": src, "noise_constant": 0.02, "output_units": "db"
        }));
        assert_eq!(outputs["clamped_cells"].as_u64().unwrap(), 2);
        for c in 0..3 {
            let v = out.get(0, 0, c);
            assert!(v.is_finite() && v != -9999.0, "cell {c} lost to -inf: {v}");
        }
        // The clamped floor is 1% of the noise power = 2e-4 -> about -37 dB.
        assert!((out.get(0, 0, 2) - 10.0 * (0.0002f64).log10()).abs() < 1e-4);
    }

    /// A noise raster overrides the scalar forms and is read per cell.
    #[test]
    fn noise_raster_takes_precedence() {
        let src = raster_of(2, 1, &[1.0, 1.0]);
        let noise = raster_of(2, 1, &[0.1, 0.5]);
        let (out, outputs) = run(json!({
            "input": src, "noise_raster": noise, "noise_constant": 0.9
        }));
        assert_eq!(outputs["noise_source"].as_str().unwrap(), "raster");
        assert!((out.get(0, 0, 0) - 0.9).abs() < 1e-6);
        assert!((out.get(0, 0, 1) - 0.5).abs() < 1e-6);
    }

    /// Noise given in dB is converted before subtraction.
    #[test]
    fn noise_units_are_honoured() {
        // -10 dB is 0.1 in linear power.
        let src = raster_of(1, 1, &[0.5]);
        let (out, _) = run(json!({
            "input": src, "noise_constant": -10.0, "noise_units": "db"
        }));
        assert!(
            (out.get(0, 0, 0) - 0.4).abs() < 1e-6,
            "expected 0.5 - 0.1, got {}",
            out.get(0, 0, 0)
        );
    }

    /// SNR is emitted even when no path is supplied.
    #[test]
    fn emits_snr_without_a_path() {
        let src = raster_of(1, 1, &[1.1]);
        let args: ToolArgs =
            serde_json::from_value(json!({"input": src, "noise_constant": 0.1})).unwrap();
        let out = RemoveThermalNoiseTool.run(&args, &ctx()).unwrap();
        let snr =
            load_input_raster(out.outputs["output_snr"].as_str().unwrap()).unwrap();
        // signal 1.0 over noise 0.1 = 10 dB.
        assert!((snr.get(0, 0, 0) - 10.0).abs() < 1e-4);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            RemoveThermalNoiseTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        // No noise source at all must be an error, not a silent no-op.
        assert!(bad(json!({"input": "a.tif"})).is_err());
        assert!(bad(json!({"input": "a.tif", "noise_profile": "0.1"})).is_err());
        assert!(bad(json!({"input": "a.tif", "noise_profile": "0.1,x"})).is_err());
        assert!(bad(json!({"input": "a.tif", "noise_constant": 0.1, "floor_fraction": 1.0})).is_err());
        assert!(bad(json!({"input": "a.tif", "noise_constant": 0.1, "noise_units": "watts"})).is_err());
        assert!(bad(json!({"input": "a.tif", "noise_constant": 0.1})).is_ok());
        assert!(bad(json!({"input": "a.tif", "noise_profile": "0.1,0.2"})).is_ok());
    }
}
