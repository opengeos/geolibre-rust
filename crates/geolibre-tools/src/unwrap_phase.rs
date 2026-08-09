//! GeoLibre tool: resolve 2-pi ambiguities in wrapped interferometric phase.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Unwrap Phase* (Image Analyst).
//!
//! ## The gap
//!
//! Interferometric phase is intrinsically wrapped into `(-pi, pi]`. Every
//! downstream use — deformation, DEM generation, displacement conversion via
//! `convert_sar_units` — needs the continuous absolute field. Nothing in either
//! registry unwraps: `raster_calculator` cannot express a global solve, and no
//! bundled filter (`lee_filter`, `enhanced_lee_filter`, `refined_lee_filter`,
//! `non_local_means_filter`) is a substitute, because unwrapping is not
//! smoothing.
//!
//! With `sar_coherence` (phase), `flatten_interferogram` (topography removal)
//! and this tool, the catalog has a minimum viable InSAR chain.
//!
//! ## Method: weighted least-squares by PCG
//!
//! Wrapped differences between 4-neighbours are the *true* gradient wherever
//! the field is smooth enough, because wrapping the difference undoes the
//! wrapping of the operands. The unwrapped field is then whichever surface
//! best matches that gradient in the least-squares sense — the solution of a
//! discrete Poisson equation, solved with preconditioned conjugate gradient.
//!
//! Pure `Vec<f64>` arithmetic; no linear-algebra crate, satisfying the repo's
//! no-heavy-dependencies constraint. Deterministic and iteration-capped, so
//! WASM-safe.
//!
//! ## Scope, deliberately
//!
//! Least-squares unwrapping is *global and smooth*: it never leaves residual
//! 2-pi jumps, but it also spreads the error from any genuine discontinuity
//! (a shear fault, a layover shadow) across the scene rather than localising
//! it. Branch-cut and minimum-cost-flow unwrappers trade differently; this
//! ships the method ArcGIS names ("Least squares PCG") and says so rather than
//! implying a residue-exact result.
//!
//! The solution is fixed only up to an additive constant — that is inherent to
//! unwrapping, not a limitation here — so a reference pixel is pinned to zero
//! and reported.
//!
//! One consequence of that: only the connected component containing the
//! reference pixel is pinned. If masking splits the valid pixels into several
//! components, each of the others keeps its own arbitrary constant, so values
//! are comparable *within* a component but not *across* components.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::args_common::{choice_or, f64_or, opt_usize, req_str, usize_or};
use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};
use crate::flatten_interferogram::wrap;
use crate::raster_stack::check_alignment_refs;
use crate::vector_common::parse_optional_str;

pub struct UnwrapPhaseTool;

impl Tool for UnwrapPhaseTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "unwrap_phase",
            display_name: "Unwrap Phase",
            summary: "Converts wrapped interferometric phase in (-pi, pi] into a continuous absolute field by weighted least squares, solved with preconditioned conjugate gradient, with optional coherence weighting and masking (ArcGIS Unwrap Phase). Nothing in either registry unwraps phase — raster_calculator cannot express the global solve and the bundled speckle filters smooth rather than unwrap — so sar_coherence's phase output and flatten_interferogram currently have no continuous endpoint.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Wrapped-phase raster in radians, or a two-band I/Q interferogram.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Unwrapped continuous phase raster. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'least_squares_pcg' (default), matching the method ArcGIS supports.",
                    required: false,
                },
                ToolParamSpec {
                    name: "coherence",
                    description: "Optional co-registered coherence raster used as per-edge solve weights, so low-coherence pixels pull the solution less.",
                    required: false,
                },
                ToolParamSpec {
                    name: "coherence_threshold",
                    description: "Mask out pixels whose coherence falls below this (default 0.3). Ignored when no coherence raster is given. Masking that splits the valid pixels into disconnected components leaves each component other than the reference pixel's with its own arbitrary additive constant.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_iterations",
                    description: "PCG iteration cap (default 500).",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Relative residual at which the solve stops (default 1e-8).",
                    required: false,
                },
                ToolParamSpec {
                    name: "reference_row",
                    description: "Row of the pixel pinned to zero, fixing the solution's additive constant. Default: the first valid pixel.",
                    required: false,
                },
                ToolParamSpec {
                    name: "reference_col",
                    description: "Column of the pinned pixel. Default: the first valid pixel.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band holding the phase for single-band input (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        choice_or(args, "method", &["least_squares_pcg"], "least_squares_pcg")?;
        let iters = usize_or(args, "max_iterations", 500)?;
        if iters == 0 {
            return Err(ToolError::Validation(
                "'max_iterations' must be at least 1".to_string(),
            ));
        }
        let tol = f64_or(args, "tolerance", 1e-8)?;
        if tol <= 0.0 {
            return Err(ToolError::Validation("'tolerance' must be > 0".to_string()));
        }
        crate::args_common::band_index(args, "band")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?.to_string();
        let output = parse_optional_output(args, "output")?;
        let max_iter = usize_or(args, "max_iterations", 500)?;
        let tolerance = f64_or(args, "tolerance", 1e-8)?;
        let threshold = f64_or(args, "coherence_threshold", 0.3)?;
        let band = crate::args_common::band_index(args, "band")?;

        let raster = load_input_raster(&input)?;
        let (rows, cols) = (raster.rows, raster.cols);
        // I/Q auto-detection applies only when the caller did not name a band.
        // Without that guard a three-band phase stack, or two unrelated phase
        // images, would be read as interleaved I and Q and the tool would
        // return atan2 of two unrelated phases rather than failing.
        let band_given = opt_usize(args, "band")?.is_some();
        let complex_input = raster.bands == 2 && !band_given;
        if raster.bands > 2 && !band_given {
            // Three or more bands is not an I/Q pair, and guessing bands 1 and
            // 2 would return plausible but wrong phase.
            return Err(ToolError::Validation(format!(
                "'input' has {} bands; supply 'band' to select the phase band, or pass a \
                 two-band I/Q raster",
                raster.bands
            )));
        }
        if !complex_input && band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "'band' {} is out of range; '{input}' has {} band(s)",
                band + 1,
                raster.bands
            )));
        }

        let coherence = match parse_optional_str(args, "coherence")? {
            Some(p) => {
                let c = load_input_raster(p)?;
                check_alignment_refs(&[&raster, &c])?;
                Some(c)
            }
            None => None,
        };

        // Valid mask: finite phase, and coherent enough when coherence is given.
        let mut phase = vec![f64::NAN; rows * cols];
        let mut weight = vec![0.0_f64; rows * cols];
        let mut valid = 0_u64;
        for r in 0..rows {
            for c in 0..cols {
                let idx = r * cols + c;
                let Some(p) = read_phase(&raster, complex_input, band, r, c) else {
                    continue;
                };
                let w = match &coherence {
                    None => 1.0,
                    Some(cr) => {
                        let g = cr.get(0, r as isize, c as isize);
                        if g == cr.nodata || !g.is_finite() || g < threshold {
                            continue;
                        }
                        g.clamp(0.0, 1.0)
                    }
                };
                phase[idx] = p;
                weight[idx] = w;
                valid += 1;
            }
        }
        if valid == 0 {
            return Err(ToolError::Execution(
                "no pixels survived the validity and coherence mask".to_string(),
            ));
        }

        let reference = pick_reference(args, &weight, rows, cols)?;
        ctx.progress.info(&format!(
            "{rows}x{cols}, {valid} valid pixel(s), reference ({}, {})",
            reference / cols,
            reference % cols
        ));

        let (solution, iterations, residual) = solve(
            &phase, &weight, rows, cols, reference, max_iter, tolerance, ctx,
        );

        let nodata = -9999.0_f64;
        let mut out = vec![nodata; rows * cols];
        for i in 0..rows * cols {
            if weight[i] > 0.0 {
                if !solution[i].is_finite() {
                    // A breakdown that still produced non-finite values must
                    // not be published as a phase field.
                    return Err(ToolError::Execution(
                        "the least-squares solve did not produce a finite field; try a larger \
                         'tolerance', fewer masked pixels, or a different reference pixel"
                            .to_string(),
                    ));
                }
                out[i] = solution[i];
            }
        }

        let out_raster = raster_like_with_data(&raster, out, nodata, DataType::F32)?;
        let out_path = write_or_store_output(out_raster, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("valid_cells".to_string(), json!(valid));
        outputs.insert("iterations".to_string(), json!(iterations));
        outputs.insert("residual".to_string(), json!(residual));
        outputs.insert("reference_row".to_string(), json!(reference / cols));
        outputs.insert("reference_col".to_string(), json!(reference % cols));
        outputs.insert("converged".to_string(), json!(residual <= tolerance));
        Ok(ToolRunResult { outputs })
    }
}

/// Weighted least-squares unwrap by preconditioned conjugate gradient.
///
/// Solves `A x = b`, where `A` is the weighted graph Laplacian over valid
/// 4-neighbour edges and `b` the divergence of the **wrapped** phase
/// differences. Wrapping each difference is the whole trick: over a smooth
/// field it recovers the true gradient regardless of how the operands wrapped.
#[allow(clippy::too_many_arguments)]
fn solve(
    phase: &[f64],
    weight: &[f64],
    rows: usize,
    cols: usize,
    reference: usize,
    max_iter: usize,
    tolerance: f64,
    ctx: &ToolContext,
) -> (Vec<f64>, usize, f64) {
    let n = rows * cols;
    // Edge weights: the harmonic-friendly product of the two endpoint weights,
    // zero whenever either endpoint is masked.
    let edge = |a: usize, b: usize| -> f64 {
        if weight[a] > 0.0 && weight[b] > 0.0 {
            weight[a] * weight[b]
        } else {
            0.0
        }
    };
    let neighbours = |i: usize| -> [(Option<usize>, f64); 4] {
        let (r, c) = (i / cols, i % cols);
        let right = (c + 1 < cols).then(|| i + 1);
        let left = (c > 0).then(|| i - 1);
        let down = (r + 1 < rows).then(|| i + cols);
        let up = (r > 0).then(|| i - cols);
        [
            (right, right.map_or(0.0, |j| edge(i, j))),
            (left, left.map_or(0.0, |j| edge(i, j))),
            (down, down.map_or(0.0, |j| edge(i, j))),
            (up, up.map_or(0.0, |j| edge(i, j))),
        ]
    };

    // Right-hand side: divergence of the wrapped gradient.
    let mut b = vec![0.0_f64; n];
    for i in 0..n {
        if weight[i] <= 0.0 {
            continue;
        }
        let mut acc = 0.0;
        for (j, w) in neighbours(i) {
            let Some(j) = j else { continue };
            if w <= 0.0 {
                continue;
            }
            acc += w * wrap(phase[j] - phase[i]);
        }
        b[i] = acc;
    }

    // Applies A. The reference pixel's row is replaced by the identity, which
    // pins the otherwise-singular system's null space (the additive constant).
    let apply = |x: &[f64], out: &mut Vec<f64>| {
        for i in 0..n {
            if weight[i] <= 0.0 {
                out[i] = 0.0;
                continue;
            }
            if i == reference {
                out[i] = x[i];
                continue;
            }
            let mut acc = 0.0;
            for (j, w) in neighbours(i) {
                let Some(j) = j else { continue };
                if w <= 0.0 {
                    continue;
                }
                acc += w * (x[j] - x[i]);
            }
            out[i] = acc;
        }
    };

    // Diagonal of A, for Jacobi preconditioning.
    let mut diag = vec![1.0_f64; n];
    for i in 0..n {
        if weight[i] <= 0.0 || i == reference {
            continue;
        }
        let d: f64 = neighbours(i)
            .iter()
            .map(|(j, w)| if j.is_some() { *w } else { 0.0 })
            .sum();
        diag[i] = if d > 0.0 { -d } else { 1.0 };
    }

    let mut rhs = b.clone();
    rhs[reference] = 0.0; // pinned to zero

    let mut x = vec![0.0_f64; n];
    let mut ax = vec![0.0_f64; n];
    apply(&x, &mut ax);
    let mut r: Vec<f64> = (0..n).map(|i| rhs[i] - ax[i]).collect();
    let mut z: Vec<f64> = (0..n).map(|i| r[i] / diag[i]).collect();
    let mut p = z.clone();
    let mut rz: f64 = (0..n).map(|i| r[i] * z[i]).sum();
    let b_norm = rhs.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-300);

    let mut ap = vec![0.0_f64; n];
    let mut iterations = 0usize;
    let mut residual = (r.iter().map(|v| v * v).sum::<f64>().sqrt()) / b_norm;

    while iterations < max_iter && residual > tolerance {
        apply(&p, &mut ap);
        let pap: f64 = (0..n).map(|i| p[i] * ap[i]).sum();
        if pap.abs() < 1e-300 {
            break; // exhausted the Krylov space
        }
        let alpha = rz / pap;
        if !alpha.is_finite() {
            break;
        }
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        for i in 0..n {
            z[i] = r[i] / diag[i];
        }
        let rz_next: f64 = (0..n).map(|i| r[i] * z[i]).sum();
        // The system is indefinite (masked rows are zero, the pinned row is
        // +1 identity, the rest are a negated Laplacian), so the Jacobi
        // preconditioner mixes signs and `rz` is not guaranteed positive. A
        // zero here would make beta inf/NaN, poison x, and — because
        // `NaN > tolerance` is false — exit the loop quietly and write NaN
        // into the output raster.
        if rz.abs() < 1e-300 || !rz_next.is_finite() {
            break;
        }
        let beta = rz_next / rz;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        rz = rz_next;
        iterations += 1;
        residual = (r.iter().map(|v| v * v).sum::<f64>().sqrt()) / b_norm;
        if iterations.is_multiple_of(32) {
            ctx.progress.progress(iterations as f64 / max_iter as f64);
        }
    }

    // Restore the pinned pixel's original phase as the datum, so the returned
    // field agrees with the input there rather than being offset by it.
    let shift = phase[reference] - x[reference];
    for v in x.iter_mut() {
        *v += shift;
    }
    (x, iterations, residual)
}

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

fn pick_reference(
    args: &ToolArgs,
    weight: &[f64],
    rows: usize,
    cols: usize,
) -> Result<usize, ToolError> {
    match (
        opt_usize(args, "reference_row")?,
        opt_usize(args, "reference_col")?,
    ) {
        (Some(r), Some(c)) => {
            if r >= rows || c >= cols {
                return Err(ToolError::Validation(format!(
                    "reference pixel ({r}, {c}) is outside the {rows}x{cols} grid"
                )));
            }
            if weight[r * cols + c] <= 0.0 {
                return Err(ToolError::Validation(format!(
                    "reference pixel ({r}, {c}) is masked out; pin a valid pixel instead"
                )));
            }
            Ok(r * cols + c)
        }
        (None, None) => weight
            .iter()
            .position(|w| *w > 0.0)
            .ok_or_else(|| ToolError::Execution("no valid pixel to pin".to_string())),
        _ => Err(ToolError::Validation(
            "supply both 'reference_row' and 'reference_col', or neither".to_string(),
        )),
    }
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

    fn run(args: Value) -> (Raster, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = UnwrapPhaseTool.run(&args, &ctx()).unwrap();
        let out = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        (out, res)
    }

    #[test]
    fn a_wrapped_linear_ramp_is_recovered_up_to_the_pinned_datum() {
        // The canonical check: a ramp spanning many fringes, wrapped, must
        // come back linear.
        let rows = 8;
        let cols = 24;
        let truth: Vec<f64> = (0..rows * cols)
            .map(|i| 0.9 * (i % cols) as f64 + 0.2 * (i / cols) as f64)
            .collect();
        let wrapped: Vec<f64> = truth.iter().map(|v| wrap(*v)).collect();
        // The input really is wrapped — otherwise the test is vacuous.
        assert!(wrapped.iter().zip(&truth).any(|(w, t)| (w - t).abs() > 1.0));

        let (out, res) = run(json!({
            "input": raster(rows, cols, vec![wrapped]),
            "max_iterations": 4000,
            "tolerance": 1e-12,
        }));
        assert_eq!(res.outputs["reference_row"], json!(0));

        // Compare against truth after removing the shared additive constant.
        let offset = out.get(0, 0, 0) - truth[0];
        for r in 0..rows {
            for c in 0..cols {
                let got = out.get(0, r as isize, c as isize) - offset;
                let want = truth[r * cols + c];
                assert!(
                    (got - want).abs() < 1e-3,
                    "({r},{c}) gave {got}, expected {want}"
                );
            }
        }
    }

    #[test]
    fn a_constant_field_unwraps_to_itself() {
        let (out, _) = run(json!({
            "input": raster(4, 4, vec![vec![0.75; 16]]),
        }));
        for r in 0..4 {
            for c in 0..4 {
                assert!((out.get(0, r, c) - 0.75).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn the_result_is_continuous_where_the_input_jumped_by_two_pi() {
        // The defining property: no neighbouring pair may differ by ~2*pi.
        let cols = 20;
        let truth: Vec<f64> = (0..cols).map(|i| 1.1 * i as f64).collect();
        let wrapped: Vec<f64> = truth.iter().map(|v| wrap(*v)).collect();
        // The wrapped input does contain 2-pi jumps.
        assert!(wrapped.windows(2).any(|w| (w[1] - w[0]).abs() > 4.0));

        let (out, _) = run(json!({
            "input": raster(1, cols, vec![wrapped]),
            "max_iterations": 4000,
            "tolerance": 1e-12,
        }));
        for c in 1..cols {
            let d = out.get(0, 0, c as isize) - out.get(0, 0, c as isize - 1);
            assert!(d.abs() < PI, "jump of {d} survived at column {c}");
        }
    }

    #[test]
    fn low_coherence_pixels_are_masked_out() {
        let coh = vec![0.9, 0.1, 0.9, 0.9];
        let (out, res) = run(json!({
            "input": raster(1, 4, vec![vec![0.1, 0.2, 0.3, 0.4]]),
            "coherence": raster(1, 4, vec![coh]),
            "coherence_threshold": 0.5,
        }));
        assert_eq!(out.get(0, 0, 1), out.nodata);
        assert_eq!(res.outputs["valid_cells"], json!(3));
    }

    #[test]
    fn an_explicit_reference_pixel_sets_the_datum() {
        let cols = 10;
        let wrapped: Vec<f64> = (0..cols).map(|i| wrap(0.8 * i as f64)).collect();
        let (out, res) = run(json!({
            "input": raster(1, cols, vec![wrapped.clone()]),
            "reference_row": 0,
            "reference_col": 5,
            "max_iterations": 2000,
        }));
        assert_eq!(res.outputs["reference_col"], json!(5));
        // The pinned pixel keeps its input phase exactly.
        assert!((out.get(0, 0, 5) - wrapped[5]).abs() < 1e-4);
    }

    #[test]
    fn complex_iq_input_is_accepted() {
        let phi = 0.4_f64;
        let (out, _) = run(json!({
            "input": raster(1, 2, vec![
                vec![phi.cos(), phi.cos()],
                vec![phi.sin(), phi.sin()],
            ]),
        }));
        assert!((out.get(0, 0, 0) - phi).abs() < 1e-5);
    }

    #[test]
    fn a_masked_reference_pixel_is_rejected_with_a_useful_message() {
        let args: ToolArgs = serde_json::from_value(json!({
            "input": raster(1, 3, vec![vec![0.1, 0.2, 0.3]]),
            "coherence": raster(1, 3, vec![vec![0.9, 0.0, 0.9]]),
            "coherence_threshold": 0.5,
            "reference_row": 0,
            "reference_col": 1,
        }))
        .unwrap();
        let err = UnwrapPhaseTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err:?}").contains("masked"), "got {err:?}");
    }

    #[test]
    fn a_mask_that_splits_the_grid_leaves_each_component_self_consistent() {
        // Only the reference pixel's component is pinned; the others keep an
        // arbitrary constant. Values must still be continuous WITHIN each side.
        let cols = 7;
        let truth: Vec<f64> = (0..cols).map(|i| 1.1 * i as f64).collect();
        let wrapped: Vec<f64> = truth.iter().map(|v| wrap(*v)).collect();
        // Column 3 is masked out, splitting the row in two.
        let coh: Vec<f64> = (0..cols).map(|i| if i == 3 { 0.0 } else { 0.9 }).collect();
        let (out, res) = run(json!({
            "input": raster(1, cols, vec![wrapped]),
            "coherence": raster(1, cols, vec![coh]),
            "coherence_threshold": 0.5,
            "max_iterations": 2000,
        }));
        assert_eq!(res.outputs["valid_cells"], json!(6));
        assert_eq!(out.get(0, 0, 3), out.nodata);
        // Within each component the gradient is recovered.
        for (a, b) in [(0, 1), (1, 2), (4, 5), (5, 6)] {
            let d = out.get(0, 0, b) - out.get(0, 0, a);
            assert!((d - 1.1).abs() < 1e-2, "gap {a}->{b} gave {d}");
        }
    }

    #[test]
    fn a_three_band_input_without_a_band_is_rejected() {
        // Three bands is not an I/Q pair; guessing bands 1 and 2 would return
        // plausible but wrong phase.
        let args: ToolArgs = serde_json::from_value(json!({
            "input": raster(1, 1, vec![vec![0.1], vec![0.2], vec![0.3]]),
        }))
        .unwrap();
        assert!(UnwrapPhaseTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn rejects_bad_parameters() {
        let r = raster(1, 2, vec![vec![0.0, 0.1]]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            UnwrapPhaseTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"input": r, "method": "branch_cut"})));
        assert!(bad(json!({"input": r, "max_iterations": 0})));
        assert!(bad(json!({"input": r, "tolerance": 0})));
    }
}
