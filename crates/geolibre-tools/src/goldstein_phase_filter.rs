//! GeoLibre tool: adaptive spectral filtering of interferometric SAR phase.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Apply Complex Data Filter* (Image
//! Analyst), which offers the Goldstein phase filter.
//!
//! ## Why the catalog needs it
//!
//! Round 17 added `unwrap_phase`, but phase unwrapping is only as good as the
//! phase it is handed. A raw interferogram is dense with noise-induced
//! *residues* — 2x2 loops whose wrapped phase differences do not sum to zero —
//! and every residue is a branch cut the unwrapper has to route around. Past a
//! few thousand residues the unwrapped surface is unusable.
//!
//! Goldstein filtering is the standard fix, and nothing in either registry can
//! stand in for it. The whole bundled `*_filter` family (`lee_filter`,
//! `frost_filter`, `kuan_filter`, `refined_lee_filter`, `wiener_filter`,
//! `gaussian_filter`) operates on real-valued **intensity**. Smoothing wrapped
//! phase directly is worse than doing nothing: averaging values that live on a
//! circle drags `+3.1` and `-3.1` radians — which are 0.08 radians apart — to a
//! meaningless 0.
//!
//! ## The algorithm (Goldstein & Werner 1998)
//!
//! The interferogram is cut into overlapping square patches. Each patch is
//! transformed with a 2-D FFT, its spectrum is smoothed to estimate where the
//! signal sits, and every bin is scaled by that smoothed magnitude raised to
//! `alpha`:
//!
//! ```text
//! Z'(u,v) = Z(u,v) * |S{ Z(u,v) }| ^ alpha
//! ```
//!
//! Because a locally-planar fringe pattern concentrates into a few bins while
//! noise spreads across all of them, this sharpens the fringe and suppresses
//! the rest. `alpha = 0` is a no-op; `alpha = 1` is maximally aggressive.
//! Patches are recombined under a triangular (Bartlett) weight so the overlap
//! blends without seams.
//!
//! Filtering is done on the **complex** signal, so the circular-mean problem
//! above never arises.
//!
//! ## Inputs and outputs
//!
//! Input is either a two-band I/Q raster or a single-band wrapped-phase raster
//! (in radians, as `sar_coherence`'s `output_phase` emits), in which case unit
//! amplitude is assumed. Both a filtered two-band complex raster and the
//! filtered wrapped phase are always written. The residue count before and
//! after is reported, which is the quantity that decides whether unwrapping
//! will succeed.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::args_common::{f64_or, usize_or};
use crate::common::{load_input_raster, parse_optional_output, write_or_store_output};
use crate::fft2::{fft2, next_pow2, Cpx};
use crate::raster_stack::raster_like_multiband;

pub struct GoldsteinPhaseFilterTool;

impl Tool for GoldsteinPhaseFilterTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "goldstein_phase_filter",
            display_name: "Goldstein Phase Filter",
            summary: "Adaptive spectral (Goldstein) filtering of interferometric SAR phase, the step that makes phase unwrapping viable (ArcGIS Apply Complex Data Filter). Round 17's unwrap_phase fails on a raw interferogram because noise residues force branch cuts; the entire bundled lee/frost/kuan/wiener/gaussian filter family works on real-valued intensity and cannot be used here, since averaging wrapped phase across the +/-pi discontinuity destroys it. Filters the complex signal patch-by-patch in the frequency domain and reports the residue count before and after.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Interferogram: a two-band I/Q raster, or a single-band wrapped-phase raster in radians (unit amplitude is assumed).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output two-band (I, Q) filtered complex raster. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_phase",
                    description: "Output single-band filtered wrapped phase in radians, in (-pi, pi]. Always produced; stored in memory when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "alpha",
                    description: "Filter strength in [0, 1] (default 0.5). 0 leaves the phase untouched; 1 is maximally aggressive.",
                    required: false,
                },
                ToolParamSpec {
                    name: "outer_window_size",
                    description: "FFT patch size in cells (default 32). Rounded up to a power of two.",
                    required: false,
                },
                ToolParamSpec {
                    name: "inner_window_size",
                    description: "Step between patches in cells (default = outer/4). Smaller means more overlap and smoother blending; must not exceed the outer window.",
                    required: false,
                },
                ToolParamSpec {
                    name: "spectrum_smoothing",
                    description: "Width in bins of the moving-average applied to the patch spectrum before exponentiation (default 3, must be odd).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        crate::args_common::req_str(args, "input")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input_path = crate::args_common::req_str(args, "input")?.to_string();
        let prm = parse_params(args)?;
        let output = parse_optional_output(args, "output")?;
        let output_phase = parse_optional_output(args, "output_phase")?;

        let raster = load_input_raster(&input_path)?;
        let (rows, cols) = (raster.rows, raster.cols);
        if rows == 0 || cols == 0 {
            return Err(ToolError::Validation("input raster is empty".to_string()));
        }
        let (mut re, mut im) = read_complex(&raster);

        ctx.progress.info(&format!(
            "{rows}x{cols}, {} input, alpha {}, patch {} step {}",
            if raster.bands >= 2 {
                "complex I/Q"
            } else {
                "wrapped phase"
            },
            prm.alpha,
            prm.patch,
            prm.step
        ));

        let residues_before = count_residues(&re, &im, rows, cols);

        filter_patches(&mut re, &mut im, rows, cols, &prm, ctx);

        let residues_after = count_residues(&re, &im, rows, cols);
        ctx.progress.info(&format!(
            "residues {residues_before} -> {residues_after}"
        ));

        // Complex output: band 0 = I, band 1 = Q.
        let nodata = -9999.0_f64;
        let mut bands: Vec<Vec<f64>> = vec![vec![nodata; rows * cols]; 2];
        let mut phase = vec![nodata; rows * cols];
        for i in 0..rows * cols {
            if re[i].is_finite() && im[i].is_finite() {
                bands[0][i] = re[i];
                bands[1][i] = im[i];
                phase[i] = im[i].atan2(re[i]);
            }
        }

        let complex_raster = raster_like_multiband(&raster, &bands, nodata, DataType::F32)?;
        let out_path = write_or_store_output(complex_raster, output)?;

        let phase_raster =
            crate::common::raster_like_with_data(&raster, phase, nodata, DataType::F32)?;
        // Emitted unconditionally: gating on a supplied path would silently
        // drop the secondary output for in-memory callers.
        let phase_path = write_or_store_output(phase_raster, output_phase)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("output_phase".to_string(), json!(phase_path));
        outputs.insert("alpha".to_string(), json!(prm.alpha));
        outputs.insert("outer_window_size".to_string(), json!(prm.patch));
        outputs.insert("inner_window_size".to_string(), json!(prm.step));
        outputs.insert("residues_before".to_string(), json!(residues_before));
        outputs.insert("residues_after".to_string(), json!(residues_after));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        Ok(ToolRunResult { outputs })
    }
}

/// Reads the input as complex samples, mapping no-data to NaN.
///
/// A two-band raster is I/Q; a single-band raster is wrapped phase, promoted to
/// unit amplitude. Rasters with more than two bands are ambiguous, so only the
/// first two are used — matching `sar_coherence`'s I/Q convention.
fn read_complex(r: &Raster) -> (Vec<f64>, Vec<f64>) {
    let (rows, cols) = (r.rows, r.cols);
    let nd = r.nodata;
    let mut re = vec![f64::NAN; rows * cols];
    let mut im = vec![f64::NAN; rows * cols];
    let complex = r.bands >= 2;
    for row in 0..rows {
        for col in 0..cols {
            let i = row * cols + col;
            let a = r.get(0, row as isize, col as isize);
            if a == nd || !a.is_finite() {
                continue;
            }
            if complex {
                let b = r.get(1, row as isize, col as isize);
                if b == nd || !b.is_finite() {
                    continue;
                }
                re[i] = a;
                im[i] = b;
            } else {
                re[i] = a.cos();
                im[i] = a.sin();
            }
        }
    }
    (re, im)
}

/// Runs the patch-wise spectral filter in place.
fn filter_patches(
    re: &mut [f64],
    im: &mut [f64],
    rows: usize,
    cols: usize,
    prm: &Params,
    ctx: &ToolContext,
) {
    let n = prm.patch;
    let step = prm.step;
    // Accumulators for the weighted overlap-add recombination.
    let mut acc_re = vec![0.0f64; rows * cols];
    let mut acc_im = vec![0.0f64; rows * cols];
    let mut acc_w = vec![0.0f64; rows * cols];

    let weight = bartlett(n);
    let mut buf: Vec<Cpx> = vec![(0.0, 0.0); n * n];
    let mut mag = vec![0.0f64; n * n];

    // Patch origins. The last patch in each direction is pulled back inside the
    // raster so edge cells still get full coverage.
    let origins = |extent: usize| -> Vec<usize> {
        let mut v = Vec::new();
        if extent <= n {
            v.push(0);
            return v;
        }
        let mut o = 0usize;
        while o + n < extent {
            v.push(o);
            o += step;
        }
        v.push(extent - n);
        v
    };
    let row_origins = origins(rows);
    let col_origins = origins(cols);
    let total = row_origins.len() * col_origins.len();
    let mut done = 0usize;

    for &r0 in &row_origins {
        for &c0 in &col_origins {
            // Load the patch; cells outside the raster or flagged no-data
            // contribute zero, which the weight accumulator accounts for.
            for pr in 0..n {
                for pc in 0..n {
                    let (gr, gc) = (r0 + pr, c0 + pc);
                    let v = if gr < rows && gc < cols {
                        let i = gr * cols + gc;
                        if re[i].is_finite() && im[i].is_finite() {
                            (re[i], im[i])
                        } else {
                            (0.0, 0.0)
                        }
                    } else {
                        (0.0, 0.0)
                    };
                    buf[pr * n + pc] = v;
                }
            }

            fft2(&mut buf, n, n, false);

            // Smoothed spectral magnitude, then the Goldstein scaling.
            for (k, c) in buf.iter().enumerate() {
                mag[k] = (c.0 * c.0 + c.1 * c.1).sqrt();
            }
            let smoothed = smooth_wrapped(&mag, n, prm.smoothing);
            // Normalising by the peak keeps the gain scale-free, so `alpha`
            // means the same thing whatever the interferogram's amplitude is.
            let peak = smoothed.iter().copied().fold(0.0f64, f64::max);
            if peak > 0.0 {
                for (k, c) in buf.iter_mut().enumerate() {
                    let g = (smoothed[k] / peak).powf(prm.alpha);
                    c.0 *= g;
                    c.1 *= g;
                }
            }

            fft2(&mut buf, n, n, true);

            for pr in 0..n {
                let gr = r0 + pr;
                if gr >= rows {
                    continue;
                }
                for pc in 0..n {
                    let gc = c0 + pc;
                    if gc >= cols {
                        continue;
                    }
                    let w = weight[pr] * weight[pc];
                    let i = gr * cols + gc;
                    let v = buf[pr * n + pc];
                    acc_re[i] += w * v.0;
                    acc_im[i] += w * v.1;
                    acc_w[i] += w;
                }
            }

            done += 1;
            ctx.progress.progress(done as f64 / total as f64);
        }
    }

    for i in 0..rows * cols {
        if !re[i].is_finite() || !im[i].is_finite() {
            continue;
        }
        if acc_w[i] > 0.0 {
            re[i] = acc_re[i] / acc_w[i];
            im[i] = acc_im[i] / acc_w[i];
        }
    }
}

/// Triangular (Bartlett) window of length `n`, strictly positive so every cell
/// in a patch contributes.
fn bartlett(n: usize) -> Vec<f64> {
    if n <= 1 {
        return vec![1.0; n.max(1)];
    }
    let half = (n - 1) as f64 / 2.0;
    (0..n)
        .map(|i| 1.0 - ((i as f64 - half) / (half + 1.0)).abs())
        .collect()
}

/// Moving average of a square spectrum with wraparound at the Nyquist edges
/// (the spectrum is periodic, so wrapping is the correct boundary rule).
fn smooth_wrapped(mag: &[f64], n: usize, width: usize) -> Vec<f64> {
    if width <= 1 {
        return mag.to_vec();
    }
    let half = (width / 2) as isize;
    let ni = n as isize;
    let mut out = vec![0.0f64; n * n];
    for r in 0..ni {
        for c in 0..ni {
            let mut sum = 0.0;
            let mut count = 0.0;
            for dr in -half..=half {
                for dc in -half..=half {
                    let rr = (r + dr).rem_euclid(ni) as usize;
                    let cc = (c + dc).rem_euclid(ni) as usize;
                    sum += mag[rr * n + cc];
                    count += 1.0;
                }
            }
            out[r as usize * n + c as usize] = sum / count;
        }
    }
    out
}

/// Counts phase residues: 2x2 loops whose wrapped phase differences do not sum
/// to zero. This is the standard interferogram quality metric — every residue
/// is an obstacle the unwrapper must route a branch cut around.
fn count_residues(re: &[f64], im: &[f64], rows: usize, cols: usize) -> usize {
    if rows < 2 || cols < 2 {
        return 0;
    }
    let phase = |i: usize| im[i].atan2(re[i]);
    let mut residues = 0usize;
    for r in 0..rows - 1 {
        for c in 0..cols - 1 {
            let idx = [r * cols + c, r * cols + c + 1, (r + 1) * cols + c + 1, (r + 1) * cols + c];
            if idx
                .iter()
                .any(|&i| !re[i].is_finite() || !im[i].is_finite())
            {
                continue;
            }
            let mut sum = 0.0;
            for k in 0..4 {
                sum += wrap(phase(idx[(k + 1) % 4]) - phase(idx[k]));
            }
            // A closed loop must sum to 0 or +/-2pi; anything but 0 is a residue.
            if (sum / (2.0 * std::f64::consts::PI)).abs() > 0.5 {
                residues += 1;
            }
        }
    }
    residues
}

/// Wraps an angle into (-pi, pi].
fn wrap(a: f64) -> f64 {
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut x = a % two_pi;
    if x > std::f64::consts::PI {
        x -= two_pi;
    } else if x <= -std::f64::consts::PI {
        x += two_pi;
    }
    x
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    alpha: f64,
    patch: usize,
    step: usize,
    smoothing: usize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let alpha = f64_or(args, "alpha", 0.5)?;
    if !(0.0..=1.0).contains(&alpha) {
        return Err(ToolError::Validation(format!(
            "'alpha' must be in [0, 1], got {alpha}"
        )));
    }

    let requested = usize_or(args, "outer_window_size", 32)?;
    if requested < 4 {
        return Err(ToolError::Validation(format!(
            "'outer_window_size' must be at least 4, got {requested}"
        )));
    }
    let patch = next_pow2(requested).ok_or_else(|| {
        ToolError::Validation(format!(
            "'outer_window_size' {requested} is too large for any FFT size"
        ))
    })?;

    let step = match crate::args_common::opt_usize(args, "inner_window_size")? {
        None => (patch / 4).max(1),
        Some(0) => {
            return Err(ToolError::Validation(
                "'inner_window_size' must be at least 1".to_string(),
            ))
        }
        Some(s) if s > patch => {
            return Err(ToolError::Validation(format!(
                "'inner_window_size' ({s}) must not exceed the outer window ({patch})"
            )))
        }
        Some(s) => s,
    };

    let smoothing = usize_or(args, "spectrum_smoothing", 3)?;
    if smoothing == 0 || smoothing % 2 == 0 {
        return Err(ToolError::Validation(format!(
            "'spectrum_smoothing' must be a positive odd number, got {smoothing}"
        )));
    }
    // `smooth_wrapped` costs patch^2 * width^2 per patch, so an unbounded width
    // blocks the call — 999 on a 32-cell patch is about 1e9 operations for
    // every patch in the raster. A width beyond the patch is meaningless
    // anyway: the moving average just wraps over the whole spectrum repeatedly.
    if smoothing > patch {
        return Err(ToolError::Validation(format!(
            "'spectrum_smoothing' ({smoothing}) must not exceed the FFT patch size \
             ({patch}); a wider average simply wraps over the whole spectrum"
        )));
    }

    Ok(Params {
        alpha,
        patch,
        step,
        smoothing,
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

    /// Deterministic value noise in [-1, 1] — no RNG, so WASM behaves like
    /// native and the test is reproducible.
    fn noise(i: usize, salt: u64) -> f64 {
        let mut x = (i as u64).wrapping_mul(6364136223846793005).wrapping_add(salt);
        x ^= x >> 33;
        x = x.wrapping_mul(0xff51afd7ed558ccd);
        x ^= x >> 33;
        ((x % 20001) as f64 / 10000.0) - 1.0
    }

    fn phase_raster(cols: usize, rows: usize, vals: &[f64]) -> String {
        make_raster(cols, rows, 1, vals)
    }

    fn make_raster(cols: usize, rows: usize, bands: usize, vals: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands,
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
        for b in 0..bands {
            for row in 0..rows {
                for col in 0..cols {
                    r.set(
                        b as isize,
                        row as isize,
                        col as isize,
                        vals[b * rows * cols + row * cols + col],
                    )
                    .unwrap();
                }
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: Value) -> (Raster, BTreeMap<String, Value>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = GoldsteinPhaseFilterTool.run(&args, &ctx()).unwrap();
        let phase = load_input_raster(out.outputs["output_phase"].as_str().unwrap()).unwrap();
        (phase, out.outputs)
    }

    /// A noisy fringe ramp: the filter must cut the residue count sharply.
    /// This is the property the tool exists for — unwrapping succeeds or fails
    /// on residue count.
    #[test]
    fn reduces_residues_on_a_noisy_ramp() {
        let (rows, cols) = (64, 64);
        let mut vals = vec![0.0; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let i = r * cols + c;
                // ~4 fringes across the scene, decorrelated by additive
                // *complex* noise — which is how interferometric noise actually
                // arises. Perturbing the phase directly cannot produce residues
                // unless the perturbation exceeds pi, whereas complex noise of
                // comparable power to the signal randomises the phase outright.
                let clean = 8.0 * std::f64::consts::PI * c as f64 / cols as f64;
                let amp = 1.2;
                let re = clean.cos() + amp * noise(i, 7);
                let im = clean.sin() + amp * noise(i, 991);
                vals[i] = im.atan2(re);
            }
        }
        let src = phase_raster(cols, rows, &vals);
        let (_, outputs) = run(json!({ "input": src, "alpha": 0.8 }));
        let before = outputs["residues_before"].as_u64().unwrap();
        let after = outputs["residues_after"].as_u64().unwrap();
        assert!(before > 50, "fixture is not noisy enough: {before} residues");
        assert!(
            after * 2 < before,
            "filter must at least halve residues, got {before} -> {after}"
        );
    }

    /// alpha = 0 is the identity: the spectral gain is 1 everywhere, so the
    /// overlap-add must reconstruct the input phase.
    #[test]
    fn alpha_zero_preserves_phase() {
        let (rows, cols) = (32, 32);
        let mut vals = vec![0.0; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                vals[r * cols + c] = wrap(0.3 * r as f64 + 0.11 * c as f64);
            }
        }
        let src = phase_raster(cols, rows, &vals);
        let (out, _) = run(json!({ "input": src, "alpha": 0.0, "outer_window_size": 16 }));
        for r in 0..rows {
            for c in 0..cols {
                let got = out.get(0, r as isize, c as isize);
                let want = vals[r * cols + c];
                // Compare on the circle: +pi and -pi are the same phase.
                let d = wrap(got - want).abs();
                assert!(d < 1e-4, "cell ({r},{c}): {got} != {want} (delta {d})");
            }
        }
    }

    /// A clean fringe pattern must survive filtering — the tool must suppress
    /// noise without eating the signal it is there to preserve.
    #[test]
    fn clean_fringes_survive() {
        let (rows, cols) = (48, 48);
        let mut vals = vec![0.0; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                vals[r * cols + c] = wrap(6.0 * std::f64::consts::PI * c as f64 / cols as f64);
            }
        }
        let src = phase_raster(cols, rows, &vals);
        let (out, outputs) = run(json!({ "input": src, "alpha": 1.0, "outer_window_size": 16 }));
        assert_eq!(outputs["residues_before"].as_u64().unwrap(), 0);
        assert_eq!(outputs["residues_after"].as_u64().unwrap(), 0);
        // Away from the patch edges the fringe gradient must be intact.
        let mut worst: f64 = 0.0;
        for r in 8..rows - 8 {
            for c in 8..cols - 8 {
                let got = out.get(0, r as isize, c as isize);
                worst = worst.max(wrap(got - vals[r * cols + c]).abs());
            }
        }
        assert!(worst < 0.2, "clean fringe distorted by {worst} rad");
    }

    /// Two-band input is read as I/Q rather than as phase.
    #[test]
    fn accepts_complex_iq_input() {
        let (rows, cols) = (16, 16);
        let mut vals = vec![0.0; 2 * rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let p = 0.2 * c as f64;
                vals[r * cols + c] = 5.0 * p.cos();
                vals[rows * cols + r * cols + c] = 5.0 * p.sin();
            }
        }
        let src = make_raster(cols, rows, 2, &vals);
        let (out, outputs) = run(json!({ "input": src, "alpha": 0.0, "outer_window_size": 8 }));
        // Phase, not amplitude, must come through.
        let got = out.get(0, 8, 5);
        assert!(
            wrap(got - 0.2 * 5.0).abs() < 1e-3,
            "I/Q phase not recovered: {got}"
        );
        let complex = load_input_raster(outputs["output"].as_str().unwrap()).unwrap();
        assert_eq!(complex.bands, 2, "complex output must have I and Q bands");
    }

    #[test]
    fn window_is_rounded_to_a_power_of_two() {
        let src = phase_raster(8, 8, &[0.0; 64]);
        let (_, outputs) = run(json!({ "input": src, "outer_window_size": 20 }));
        assert_eq!(outputs["outer_window_size"].as_u64().unwrap(), 32);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            GoldsteinPhaseFilterTool.validate(&args)
        };
        assert!(bad(json!({})).is_err()); // missing input
        assert!(bad(json!({"input": "a.tif", "alpha": 1.5})).is_err());
        assert!(bad(json!({"input": "a.tif", "alpha": -0.1})).is_err());
        assert!(bad(json!({"input": "a.tif", "outer_window_size": 2})).is_err());
        assert!(bad(json!({"input": "a.tif", "spectrum_smoothing": 4})).is_err());
        // Unbounded smoothing would block the call for minutes per patch.
        assert!(
            bad(json!({"input": "a.tif", "outer_window_size": 32, "spectrum_smoothing": 999}))
                .is_err()
        );
        // inner window larger than the (rounded) outer window
        assert!(
            bad(json!({"input": "a.tif", "outer_window_size": 16, "inner_window_size": 64}))
                .is_err()
        );
        assert!(bad(json!({"input": "a.tif", "alpha": 0.7})).is_ok());
    }
}
