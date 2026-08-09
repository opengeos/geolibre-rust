//! GeoLibre tool: remove topographic phase from a wrapped interferogram.
//!
//! Pure-Rust counterpart of the flattening half of ArcGIS Pro's *Generate
//! Interferogram* (Image Analyst) — its `TOPO` mode and
//! `out_flattened_interferogram` output.
//!
//! ## Why only the flattening half
//!
//! `sar_coherence` (round 16) **already produces the interferogram**: its
//! `output_phase` parameter emits "the interferometric phase raster in radians,
//! wrapped to (-pi, pi]" from the complex conjugate product, which is exactly
//! Generate Interferogram's core. What it cannot do — and what nothing else in
//! either registry can do — is remove the topographic component.
//!
//! That component dominates. Over any relief the fringes you see are mostly
//! terrain, not the deformation signal, so an unflattened interferogram is
//! close to unusable for deformation work and `unwrap_phase` would be
//! unwrapping mostly-terrain fringes.
//!
//! `raster_calculator` can subtract two rasters but cannot *simulate* the
//! topographic phase, which is the actual content of this tool.
//!
//! ## The relation
//!
//! ```text
//! phi_topo = -(4*pi / lambda) * (B_perp * (h - h_ref)) / (R * sin(theta))
//! phi_flat = wrap(phi_input - phi_topo)
//! ```
//!
//! `h_ref` defaults to the DEM mean so the output is centred rather than
//! carrying an arbitrary global offset that the caller would have to remove.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::args_common::{f64_or, opt_f64, opt_positive_f64, req_str};
use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};
use crate::raster_stack::check_alignment_refs;
use crate::vector_common::parse_optional_str;

/// Sentinel-1 C-band wavelength in metres.
const C_BAND_WAVELENGTH_M: f64 = 0.05546576;

pub struct FlattenInterferogramTool;

impl Tool for FlattenInterferogramTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "flatten_interferogram",
            display_name: "Flatten Interferogram",
            summary: "Removes the topographic phase component from a wrapped interferogram using a DEM and the interferometric baseline, leaving the deformation and atmosphere signal (the TOPO flattening mode of ArcGIS Generate Interferogram). sar_coherence already emits the interferogram itself via output_phase, but nothing in either registry can simulate and subtract the topographic phase — and over any relief that component dominates, so an unflattened interferogram is unusable for deformation work.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Wrapped interferometric phase raster in radians (as sar_coherence's output_phase emits), or a two-band I/Q interferogram.",
                    required: true,
                },
                ToolParamSpec {
                    name: "dem",
                    description: "Co-registered elevation raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Flattened wrapped-phase raster, re-wrapped to (-pi, pi]. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "perpendicular_baseline",
                    description: "Perpendicular baseline B_perp in metres. Required: it sets the height-to-phase sensitivity, and there is no defensible default.",
                    required: true,
                },
                ToolParamSpec {
                    name: "wavelength",
                    description: "Radar wavelength in metres (default 0.05546576, Sentinel-1 C-band).",
                    required: false,
                },
                ToolParamSpec {
                    name: "incidence_angle",
                    description: "Incidence angle in degrees: a constant (default 39.0, typical Sentinel-1 IW mid-swath) or a co-registered raster path.",
                    required: false,
                },
                ToolParamSpec {
                    name: "slant_range",
                    description: "Sensor-to-ground slant range R in metres (default 800000).",
                    required: false,
                },
                ToolParamSpec {
                    name: "reference_elevation",
                    description: "Elevation treated as zero topographic phase. Default: the DEM mean, so the output carries no arbitrary global offset.",
                    required: false,
                },
                ToolParamSpec {
                    name: "out_topographic_phase",
                    description: "Simulated topographic phase raster, WRAPPED to (-pi, pi] like the flattened output — not the absolute simulated phase. Always produced; falls back to an in-memory handle when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band holding the phase for single-band input (default 1). Ignored for two-band I/Q input.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        req_str(args, "dem")?;
        opt_f64(args, "perpendicular_baseline")?.ok_or_else(|| {
            ToolError::Validation(
                "missing required parameter 'perpendicular_baseline' (metres)".to_string(),
            )
        })?;
        opt_positive_f64(args, "wavelength")?;
        opt_positive_f64(args, "slant_range")?;
        crate::args_common::band_index(args, "band")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?.to_string();
        let dem_path = req_str(args, "dem")?.to_string();
        let output = parse_optional_output(args, "output")?;
        let topo_out = parse_optional_output(args, "out_topographic_phase")?;
        let b_perp = opt_f64(args, "perpendicular_baseline")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'perpendicular_baseline'".to_string())
        })?;
        let wavelength = opt_positive_f64(args, "wavelength")?.unwrap_or(C_BAND_WAVELENGTH_M);
        let slant_range = opt_positive_f64(args, "slant_range")?.unwrap_or(800_000.0);
        let band = crate::args_common::band_index(args, "band")?;

        let phase_raster = load_input_raster(&input)?;
        let dem = load_input_raster(&dem_path)?;
        check_alignment_refs(&[&phase_raster, &dem])?;

        // A two-band input is an I/Q interferogram; recover its argument. The
        // auto-detection applies ONLY when the caller did not name a band —
        // otherwise a multi-band real-phase stack with band = 2 would silently
        // return atan2(band1, band0) instead, which looks plausible and is
        // wrong. An explicit band on exactly two bands is ambiguous, so it is
        // rejected rather than guessed.
        let band_given = crate::args_common::opt_usize(args, "band")?.is_some();
        let complex_input = phase_raster.bands >= 2 && !band_given;
        let (rows, cols) = (phase_raster.rows, phase_raster.cols);
        if band_given && phase_raster.bands >= 2 {
            ctx.progress
                .info("'band' was supplied, so the input is read as real phase rather than I/Q");
        }
        if !complex_input && band as usize >= phase_raster.bands {
            return Err(ToolError::Validation(format!(
                "'band' {} is out of range; '{input}' has {} band(s)",
                band + 1,
                phase_raster.bands
            )));
        }

        // Incidence angle: constant or raster, resolved the same way
        // apply_radiometric_calibration does.
        let incidence_path =
            parse_optional_str(args, "incidence_angle")?.filter(|s| s.parse::<f64>().is_err());
        let incidence_raster = match incidence_path {
            Some(p) => {
                let r = load_input_raster(p)?;
                check_alignment_refs(&[&phase_raster, &r])?;
                Some(r)
            }
            None => None,
        };
        // Only read the scalar form when no raster was supplied: `f64_or` would
        // otherwise re-parse the raster path as a number and fail the run.
        let incidence_const = if incidence_raster.is_some() {
            f64::NAN // unused; every cell reads the raster
        } else {
            f64_or(args, "incidence_angle", 39.0)?
        };

        // Reference elevation: the DEM mean unless told otherwise, so the
        // flattened field is centred rather than offset by an arbitrary datum.
        let reference = match opt_f64(args, "reference_elevation")? {
            Some(v) => v,
            None => dem_mean(&dem)
                .ok_or_else(|| ToolError::Execution("DEM holds no valid elevations".to_string()))?,
        };

        ctx.progress.info(&format!(
            "{rows}x{cols}, B_perp {b_perp} m, lambda {wavelength} m, h_ref {reference:.2} m"
        ));

        let nodata = -9999.0_f64;
        let mut flat = vec![nodata; rows * cols];
        let mut topo = vec![nodata; rows * cols];
        let mut valid = 0_u64;

        for r in 0..rows {
            for c in 0..cols {
                let idx = r * cols + c;
                let Some(phi) = read_phase(&phase_raster, complex_input, band, r, c) else {
                    continue;
                };
                let h = dem.get(0, r as isize, c as isize);
                if h == dem.nodata || !h.is_finite() {
                    continue;
                }
                let theta_deg = match &incidence_raster {
                    Some(ir) => {
                        let v = ir.get(0, r as isize, c as isize);
                        if v == ir.nodata || !v.is_finite() {
                            continue;
                        }
                        v
                    }
                    None => incidence_const,
                };
                if theta_deg <= 0.0 || theta_deg >= 90.0 {
                    continue;
                }
                let phi_topo = -(4.0 * std::f64::consts::PI / wavelength)
                    * (b_perp * (h - reference))
                    / (slant_range * theta_deg.to_radians().sin());
                topo[idx] = wrap(phi_topo);
                flat[idx] = wrap(phi - phi_topo);
                valid += 1;
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let flat_raster = raster_like_with_data(&phase_raster, flat, nodata, DataType::F32)?;
        let flat_path = write_or_store_output(flat_raster, output)?;
        // Repo convention (create_overpass): emit the secondary output
        // unconditionally so a caller with no scratch path still gets it.
        let topo_raster = raster_like_with_data(&phase_raster, topo, nodata, DataType::F32)?;
        let topo_path = write_or_store_output(topo_raster, topo_out)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(flat_path));
        outputs.insert("out_topographic_phase".to_string(), json!(topo_path));
        outputs.insert("reference_elevation".to_string(), json!(reference));
        outputs.insert("valid_cells".to_string(), json!(valid));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        Ok(ToolRunResult { outputs })
    }
}

/// Reads the input phase, from either a real phase band or an I/Q pair.
fn read_phase(r: &Raster, complex_input: bool, band: isize, row: usize, col: usize) -> Option<f64> {
    if complex_input {
        let i = r.get(0, row as isize, col as isize);
        let q = r.get(1, row as isize, col as isize);
        if i == r.nodata || q == r.nodata || !i.is_finite() || !q.is_finite() {
            return None;
        }
        Some(q.atan2(i))
    } else {
        let v = r.get(band, row as isize, col as isize);
        (v != r.nodata && v.is_finite()).then_some(v)
    }
}

/// Wraps a phase into (-pi, pi], the convention `sar_coherence` already uses.
pub(crate) fn wrap(phi: f64) -> f64 {
    use std::f64::consts::PI;
    let mut p = (phi + PI).rem_euclid(2.0 * PI) - PI;
    // rem_euclid maps exactly -pi to -pi; the convention is the open end there.
    if p <= -PI {
        p += 2.0 * PI;
    }
    p
}

fn dem_mean(dem: &Raster) -> Option<f64> {
    let mut sum = 0.0;
    let mut n = 0u64;
    for r in 0..dem.rows {
        for c in 0..dem.cols {
            let v = dem.get(0, r as isize, c as isize);
            if v != dem.nodata && v.is_finite() {
                sum += v;
                n += 1;
            }
        }
    }
    (n > 0).then(|| sum / n as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use std::f64::consts::PI;
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

    fn raster(rows: usize, cols: usize, bands: Vec<Vec<f64>>) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: bands.len(),
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
        for (b, band) in bands.iter().enumerate() {
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

    fn run(args: Value) -> (Raster, Raster, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = FlattenInterferogramTool.run(&args, &ctx()).unwrap();
        let flat = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        let topo =
            load_input_raster(res.outputs["out_topographic_phase"].as_str().unwrap()).unwrap();
        (flat, topo, res)
    }

    /// Forward-simulates the topographic phase the tool should remove.
    fn simulate(h: f64, h_ref: f64, b_perp: f64, lambda: f64, theta_deg: f64, r: f64) -> f64 {
        -(4.0 * PI / lambda) * (b_perp * (h - h_ref)) / (r * theta_deg.to_radians().sin())
    }

    const LAMBDA: f64 = 0.05546576;
    const THETA: f64 = 39.0;
    const RANGE: f64 = 800_000.0;
    const BPERP: f64 = 150.0;

    #[test]
    fn a_pure_topographic_signal_flattens_to_zero() {
        // Build a DEM ramp, forward-simulate its wrapped topographic phase,
        // and check the tool recovers a flat field.
        let heights: Vec<f64> = (0..40).map(|i| i as f64 * 50.0).collect();
        let h_ref = heights.iter().sum::<f64>() / heights.len() as f64;
        let phase: Vec<f64> = heights
            .iter()
            .map(|h| wrap(simulate(*h, h_ref, BPERP, LAMBDA, THETA, RANGE)))
            .collect();

        let (flat, _, _) = run(json!({
            "input": raster(1, 40, vec![phase]),
            "dem": raster(1, 40, vec![heights]),
            "perpendicular_baseline": BPERP,
        }));
        for c in 0..40 {
            assert!(
                flat.get(0, 0, c).abs() < 1e-4,
                "cell {c} left residual {}",
                flat.get(0, 0, c)
            );
        }
    }

    #[test]
    fn a_superimposed_deformation_signal_survives_flattening() {
        // The point of the tool: terrain goes, signal stays.
        let heights: Vec<f64> = (0..40).map(|i| i as f64 * 50.0).collect();
        let h_ref = heights.iter().sum::<f64>() / heights.len() as f64;
        let deformation = 0.7_f64;
        let phase: Vec<f64> = heights
            .iter()
            .map(|h| wrap(simulate(*h, h_ref, BPERP, LAMBDA, THETA, RANGE) + deformation))
            .collect();

        let (flat, _, _) = run(json!({
            "input": raster(1, 40, vec![phase]),
            "dem": raster(1, 40, vec![heights]),
            "perpendicular_baseline": BPERP,
        }));
        for c in 0..40 {
            assert!(
                (flat.get(0, 0, c) - deformation).abs() < 1e-4,
                "cell {c} gave {} not {deformation}",
                flat.get(0, 0, c)
            );
        }
    }

    #[test]
    fn flat_terrain_leaves_the_input_phase_untouched() {
        // Constant DEM means zero topographic phase everywhere.
        let phase = vec![0.3, -1.2, 2.9];
        let (flat, topo, _) = run(json!({
            "input": raster(1, 3, vec![phase.clone()]),
            "dem": raster(1, 3, vec![vec![100.0; 3]]),
            "perpendicular_baseline": BPERP,
        }));
        for (c, want) in phase.iter().enumerate() {
            assert!((flat.get(0, 0, c as isize) - want).abs() < 1e-5);
            assert!(topo.get(0, 0, c as isize).abs() < 1e-9);
        }
    }

    #[test]
    fn a_zero_baseline_removes_nothing() {
        // B_perp = 0 means no height sensitivity at all.
        let phase = vec![0.5, -0.5];
        let (flat, _, _) = run(json!({
            "input": raster(1, 2, vec![phase.clone()]),
            "dem": raster(1, 2, vec![vec![0.0, 1000.0]]),
            "perpendicular_baseline": 0.0,
        }));
        for (c, want) in phase.iter().enumerate() {
            assert!((flat.get(0, 0, c as isize) - want).abs() < 1e-5);
        }
    }

    #[test]
    fn the_output_stays_wrapped_into_the_principal_interval() {
        let heights: Vec<f64> = (0..50).map(|i| i as f64 * 200.0).collect();
        let (flat, topo, _) = run(json!({
            "input": raster(1, 50, vec![vec![0.0; 50]]),
            "dem": raster(1, 50, vec![heights]),
            "perpendicular_baseline": 300.0,
        }));
        for c in 0..50 {
            for v in [flat.get(0, 0, c), topo.get(0, 0, c)] {
                assert!(
                    v > -PI - 1e-12 && v <= PI + 1e-12,
                    "phase {v} escaped (-pi, pi]"
                );
            }
        }
    }

    #[test]
    fn complex_iq_input_is_accepted_and_its_argument_used() {
        // I/Q encoding of phase 0.6 must behave like the real phase 0.6.
        let phi = 0.6_f64;
        let (from_complex, _, _) = run(json!({
            "input": raster(1, 1, vec![vec![phi.cos()], vec![phi.sin()]]),
            "dem": raster(1, 1, vec![vec![100.0]]),
            "perpendicular_baseline": BPERP,
        }));
        assert!((from_complex.get(0, 0, 0) - phi).abs() < 1e-5);
    }

    #[test]
    fn nodata_in_the_dem_leaves_the_cell_unresolved() {
        let (flat, _, res) = run(json!({
            "input": raster(1, 2, vec![vec![0.1, 0.2]]),
            "dem": raster(1, 2, vec![vec![-9999.0, 100.0]]),
            "perpendicular_baseline": BPERP,
        }));
        assert_eq!(flat.get(0, 0, 0), flat.nodata);
        assert_eq!(res.outputs["valid_cells"], json!(1));
    }

    #[test]
    fn the_topographic_phase_output_exists_without_a_path() {
        let args: ToolArgs = serde_json::from_value(json!({
            "input": raster(1, 1, vec![vec![0.0]]),
            "dem": raster(1, 1, vec![vec![10.0]]),
            "perpendicular_baseline": BPERP,
        }))
        .unwrap();
        let res = FlattenInterferogramTool.run(&args, &ctx()).unwrap();
        let p = res.outputs["out_topographic_phase"].as_str().unwrap();
        assert!(load_input_raster(p).is_ok());
    }

    #[test]
    fn a_raster_incidence_angle_is_accepted() {
        // Regression: `f64_or` ran unconditionally on `incidence_angle`, so a
        // raster path was re-parsed as a number and the whole run failed —
        // the documented raster form was entirely unusable.
        let heights = vec![0.0, 100.0, 200.0];
        let (flat, _, _) = run(json!({
            "input": raster(1, 3, vec![vec![0.0, 0.0, 0.0]]),
            "dem": raster(1, 3, vec![heights]),
            "perpendicular_baseline": BPERP,
            "incidence_angle": raster(1, 3, vec![vec![30.0, 39.0, 45.0]]),
        }));
        for c in 0..3 {
            assert_ne!(flat.get(0, 0, c), flat.nodata, "cell {c} unresolved");
        }
    }

    #[test]
    fn an_explicit_band_is_honoured_on_a_multiband_input() {
        // Regression: bands >= 2 was treated as I/Q unconditionally, so a
        // real-phase stack silently returned atan2 of two unrelated bands.
        let (flat, _, _) = run(json!({
            "input": raster(1, 1, vec![vec![0.0], vec![0.25]]),
            "dem": raster(1, 1, vec![vec![100.0]]),
            "perpendicular_baseline": BPERP,
            "band": 2,
        }));
        assert!(
            (flat.get(0, 0, 0) - 0.25).abs() < 1e-5,
            "band 2 was not read as real phase: {}",
            flat.get(0, 0, 0)
        );
    }

    #[test]
    fn rejects_bad_parameters() {
        let r = raster(1, 1, vec![vec![0.0]]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            FlattenInterferogramTool.validate(&args).is_err()
        };
        assert!(bad(json!({"dem": r})));
        assert!(bad(json!({"input": r, "dem": r})));
        assert!(bad(
            json!({"input": r, "dem": r, "perpendicular_baseline": 100, "wavelength": 0})
        ));
    }
}
