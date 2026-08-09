//! GeoLibre tool: convert a SAR raster between amplitude, intensity, linear
//! power, decibel and displacement representations.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Convert SAR Units* (Image Analyst).
//!
//! ## Why this is not `raster_calculator`
//!
//! The arithmetic is one line; the *semantics* are the tool. SAR products
//! arrive in mutually incompatible units and every threshold, filter and index
//! in the catalog silently assumes one of them. Running `lee_filter` on dB data
//! instead of intensity is simply wrong — speckle is multiplicative in
//! intensity and additive in dB, which is the entire premise those filters are
//! built on. Likewise the ratio-based SAR indices are meaningless on dB input.
//!
//! `raster_calculator` will happily evaluate `10 * log10(x)` and hand back
//! `-inf` for every zero-valued cell, which then poisons any downstream
//! statistic. Here non-positive input maps to **no-data** instead, which is the
//! guard that makes the conversion safe to chain.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::args_common::{band_index, opt_choice, opt_positive_f64, req_str};
use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};

/// Sentinel-1 C-band wavelength in metres, the default for phase→displacement.
const C_BAND_WAVELENGTH_M: f64 = 0.05546576;

pub struct ConvertSarUnitsTool;

impl Tool for ConvertSarUnitsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "convert_sar_units",
            display_name: "Convert SAR Units",
            summary: "Converts a SAR raster between amplitude, intensity, linear power, decibel, complex-to-intensity and unwrapped-phase-to-displacement representations (ArcGIS Convert SAR Units). raster_calculator can evaluate the arithmetic but carries none of the semantics or guards — it returns -inf for zero-valued cells on a log conversion, and nothing in the catalog records which unit a SAR raster is in, so speckle filters and ratio indices are routinely applied to the wrong representation.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input SAR raster. For 'complex_to_intensity' this is a two-band I/Q raster (the multilook / sar_coherence convention).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output raster path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "conversion",
                    description: "'linear_to_db', 'db_to_linear', 'amplitude_to_intensity', 'intensity_to_amplitude', 'complex_to_intensity', or 'phase_to_displacement'.",
                    required: true,
                },
                ToolParamSpec {
                    name: "wavelength",
                    description: "Radar wavelength in metres for 'phase_to_displacement' (default 0.05546576, Sentinel-1 C-band).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read (default 1). Ignored for 'complex_to_intensity', which always reads bands 1 and 2 as I and Q.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        parse_conversion(args)?;
        opt_positive_f64(args, "wavelength")?;
        band_index(args, "band")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?.to_string();
        let conversion = parse_conversion(args)?;
        let wavelength = opt_positive_f64(args, "wavelength")?.unwrap_or(C_BAND_WAVELENGTH_M);
        let band = band_index(args, "band")?;
        let output = parse_optional_output(args, "output")?;

        let raster = load_input_raster(&input)?;
        let (rows, cols) = (raster.rows, raster.cols);

        if conversion == Conversion::ComplexToIntensity && raster.bands < 2 {
            return Err(ToolError::Validation(format!(
                "'complex_to_intensity' needs a two-band I/Q raster; '{input}' has {} band(s)",
                raster.bands
            )));
        }
        if band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "'band' {} is out of range; '{input}' has {} band(s)",
                band + 1,
                raster.bands
            )));
        }

        ctx.progress.info(&format!(
            "{rows}x{cols}, conversion {}",
            conversion.label()
        ));

        let nodata = -9999.0_f64;
        let mut out = vec![nodata; rows * cols];
        let mut dropped = 0usize;

        for r in 0..rows {
            for c in 0..cols {
                let idx = r * cols + c;
                let v = match conversion {
                    Conversion::ComplexToIntensity => {
                        let i = raster.get(0, r as isize, c as isize);
                        let q = raster.get(1, r as isize, c as isize);
                        if !valid(i, raster.nodata) || !valid(q, raster.nodata) {
                            continue;
                        }
                        Some(i * i + q * q)
                    }
                    _ => {
                        let v = raster.get(band, r as isize, c as isize);
                        if !valid(v, raster.nodata) {
                            continue;
                        }
                        conversion.apply(v, wavelength)
                    }
                };
                match v {
                    Some(v) => out[idx] = v,
                    // Domain violations (log of a non-positive value, sqrt of a
                    // negative intensity) become no-data rather than -inf/NaN.
                    None => dropped += 1,
                }
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let out_r = raster_like_with_data(&raster, out, nodata, DataType::F32)?;
        let out_path = write_or_store_output(out_r, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("conversion".to_string(), json!(conversion.label()));
        outputs.insert("out_of_domain_cells".to_string(), json!(dropped));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        if conversion == Conversion::PhaseToDisplacement {
            outputs.insert("wavelength".to_string(), json!(wavelength));
        }
        Ok(ToolRunResult { outputs })
    }
}

fn valid(v: f64, nodata: f64) -> bool {
    v != nodata && v.is_finite()
}

// ── Parameters ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Conversion {
    LinearToDb,
    DbToLinear,
    AmplitudeToIntensity,
    IntensityToAmplitude,
    ComplexToIntensity,
    PhaseToDisplacement,
}

impl Conversion {
    fn label(self) -> &'static str {
        match self {
            Conversion::LinearToDb => "linear_to_db",
            Conversion::DbToLinear => "db_to_linear",
            Conversion::AmplitudeToIntensity => "amplitude_to_intensity",
            Conversion::IntensityToAmplitude => "intensity_to_amplitude",
            Conversion::ComplexToIntensity => "complex_to_intensity",
            Conversion::PhaseToDisplacement => "phase_to_displacement",
        }
    }

    /// Applies the conversion, returning `None` for out-of-domain input.
    fn apply(self, v: f64, wavelength: f64) -> Option<f64> {
        match self {
            // log10 of a non-positive power is undefined; -inf would silently
            // destroy every downstream statistic, so drop the cell instead.
            Conversion::LinearToDb => (v > 0.0).then(|| 10.0 * v.log10()),
            Conversion::DbToLinear => Some(10.0_f64.powf(v / 10.0)),
            Conversion::AmplitudeToIntensity => Some(v * v),
            // Intensity is a power and cannot be negative.
            Conversion::IntensityToAmplitude => (v >= 0.0).then(|| v.sqrt()),
            Conversion::ComplexToIntensity => None, // handled by the caller
            // Line-of-sight displacement from unwrapped phase. The sign
            // convention is ArcGIS's: positive phase = motion away from the
            // sensor = negative displacement.
            Conversion::PhaseToDisplacement => {
                Some(-v * wavelength / (4.0 * std::f64::consts::PI))
            }
        }
    }
}

fn parse_conversion(args: &ToolArgs) -> Result<Conversion, ToolError> {
    // Required: an unset conversion is an error, not an arbitrary default.
    let raw = opt_choice(args, "conversion").ok_or_else(|| {
        ToolError::Validation("missing required parameter 'conversion'".to_string())
    })?;
    Ok(match raw.as_str() {
        "linear_to_db" => Conversion::LinearToDb,
        "db_to_linear" => Conversion::DbToLinear,
        "amplitude_to_intensity" => Conversion::AmplitudeToIntensity,
        "intensity_to_amplitude" => Conversion::IntensityToAmplitude,
        "complex_to_intensity" => Conversion::ComplexToIntensity,
        "phase_to_displacement" => Conversion::PhaseToDisplacement,
        other => {
            return Err(ToolError::Validation(format!(
                "'conversion' must be one of linear_to_db|db_to_linear|amplitude_to_intensity|\
                 intensity_to_amplitude|complex_to_intensity|phase_to_displacement, got '{other}'"
            )))
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Output rasters are F32, so value comparisons use a **relative**
    /// tolerance. An absolute 1e-9 check would fail on any value that is not
    /// exactly representable in f32 (lambda/2, 0.8, ...) even when the
    /// arithmetic is exactly right.
    fn close(actual: f64, expect: f64) -> bool {
        (actual - expect).abs() <= 1e-6 * expect.abs().max(1.0)
    }

    fn multiband(rows: usize, cols: usize, data: Vec<Vec<f64>>) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: data.len(),
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
        for (b, band) in data.iter().enumerate() {
            for row in 0..rows {
                for col in 0..cols {
                    r.set(
                        b as isize,
                        row as isize,
                        col as isize,
                        band[row * cols + col],
                    )
                    .unwrap();
                }
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn raster(data: &[f64]) -> String {
        multiband(1, data.len(), vec![data.to_vec()])
    }

    fn run(args: Value) -> (Raster, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = ConvertSarUnitsTool.run(&args, &ctx()).unwrap();
        let out = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        (out, res)
    }

    #[test]
    fn linear_and_db_round_trip() {
        let (db, _) = run(json!({"input": raster(&[1.0, 10.0, 100.0]), "conversion": "linear_to_db"}));
        assert!(close(db.get(0, 0, 0), 0.0));
        assert!(close(db.get(0, 0, 1), 10.0));
        assert!(close(db.get(0, 0, 2), 20.0));

        let (lin, _) = run(json!({"input": raster(&[0.0, 10.0, 20.0]), "conversion": "db_to_linear"}));
        assert!(close(lin.get(0, 0, 0), 1.0));
        assert!(close(lin.get(0, 0, 1), 10.0));
        assert!(close(lin.get(0, 0, 2), 100.0));
    }

    #[test]
    fn non_positive_power_becomes_nodata_not_negative_infinity() {
        // The whole reason this tool exists rather than raster_calculator.
        let (out, res) = run(json!({"input": raster(&[0.0, -3.0, 4.0]), "conversion": "linear_to_db"}));
        assert_eq!(out.get(0, 0, 0), out.nodata);
        assert_eq!(out.get(0, 0, 1), out.nodata);
        assert!(out.get(0, 0, 2).is_finite() && out.get(0, 0, 2) != out.nodata);
        assert_eq!(res.outputs["out_of_domain_cells"], json!(2));
    }

    #[test]
    fn amplitude_and_intensity_round_trip() {
        let (inten, _) = run(json!({"input": raster(&[2.0, 3.0]), "conversion": "amplitude_to_intensity"}));
        assert!(close(inten.get(0, 0, 0), 4.0));
        assert!(close(inten.get(0, 0, 1), 9.0));

        let (amp, _) = run(json!({"input": raster(&[4.0, 9.0]), "conversion": "intensity_to_amplitude"}));
        assert!(close(amp.get(0, 0, 0), 2.0));
        assert!(close(amp.get(0, 0, 1), 3.0));
    }

    #[test]
    fn negative_intensity_is_rejected_rather_than_producing_nan() {
        let (out, res) = run(json!({"input": raster(&[-1.0]), "conversion": "intensity_to_amplitude"}));
        assert_eq!(out.get(0, 0, 0), out.nodata);
        assert_eq!(res.outputs["out_of_domain_cells"], json!(1));
    }

    #[test]
    fn complex_to_intensity_uses_both_bands() {
        // I=3, Q=4 gives intensity 25 — not 3, and not 7.
        let src = multiband(1, 1, vec![vec![3.0], vec![4.0]]);
        let (out, _) = run(json!({"input": src, "conversion": "complex_to_intensity"}));
        assert!(close(out.get(0, 0, 0), 25.0));
    }

    #[test]
    fn complex_to_intensity_rejects_single_band_input() {
        let args: ToolArgs = serde_json::from_value(
            json!({"input": raster(&[1.0]), "conversion": "complex_to_intensity"}),
        )
        .unwrap();
        assert!(ConvertSarUnitsTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn phase_to_displacement_uses_the_quarter_wavelength_relation() {
        // A full 2-pi fringe is half a wavelength of line-of-sight motion, and
        // the sign is negative for positive phase.
        let lambda = 0.05546576;
        let (out, _) = run(json!({
            "input": raster(&[2.0 * std::f64::consts::PI]),
            "conversion": "phase_to_displacement",
            "wavelength": lambda,
        }));
        assert!(close(out.get(0, 0, 0), -lambda / 2.0));
    }

    #[test]
    fn nodata_survives_every_conversion() {
        for conv in [
            "linear_to_db",
            "db_to_linear",
            "amplitude_to_intensity",
            "intensity_to_amplitude",
            "phase_to_displacement",
        ] {
            let (out, _) = run(json!({"input": raster(&[-9999.0, 5.0]), "conversion": conv}));
            assert_eq!(out.get(0, 0, 0), out.nodata, "conversion {conv}");
            assert_ne!(out.get(0, 0, 1), out.nodata, "conversion {conv}");
        }
    }

    #[test]
    fn rejects_missing_or_unknown_conversion() {
        let src = raster(&[1.0]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ConvertSarUnitsTool.validate(&args).is_err()
        };
        assert!(bad(json!({"input": src})));
        assert!(bad(json!({"input": src, "conversion": "nope"})));
        assert!(bad(
            json!({"input": src, "conversion": "linear_to_db", "wavelength": 0})
        ));
    }
}
