//! GeoLibre tool: automatic image-to-image co-registration.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Apply Coregistration* (Image
//! Analyst).
//!
//! ## Why the catalog needs it
//!
//! `sar_coherence`'s own documentation says its inputs "must already be
//! co-registered" — and nothing in either registry could produce that.
//! `flatten_interferogram` and `goldstein_phase_filter` inherit the same
//! precondition: a sub-pixel misregistration between two SAR acquisitions
//! destroys coherence outright, because the interferometric phase is the
//! difference of two signals that must describe the *same* ground cell.
//!
//! The nearest existing tools do not close the gap:
//! * `georeference_raster_from_control_points` needs control points a human has
//!   already picked;
//! * `image_correlation` reports a single whole-image similarity score and
//!   applies no transform;
//! * `warp_raster` applies a transform someone else has to supply.
//!
//! This tool measures the misregistration and applies it.
//!
//! ## Method
//!
//! A regular grid of tiles is cut from the reference. Each tile is matched
//! against the secondary by **normalised cross-correlation** over an integer
//! search window, then refined to sub-pixel accuracy by fitting a parabola
//! through the correlation peak and its two neighbours in each axis. The
//! resulting tie points are fitted with a least-squares transform (translation,
//! affine, or second-order polynomial) under iterative outlier trimming, and
//! the secondary is resampled onto the reference grid.
//!
//! Outlier rejection is deterministic trimming (refit, drop residuals beyond
//! `outlier_sigma` standard deviations, repeat) rather than RANSAC, so the
//! result does not depend on a random seed — the WASM builds have no RNG the
//! native build shares.
//!
//! **Every band** of the secondary is resampled with the fitted model, so a
//! two-band I/Q raster stays a valid complex image (I and Q are interpolated
//! separately, which is correct for the sub-pixel shifts co-registration
//! removes).

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Feature, Geometry, GeometryType, Layer};

use crate::args_common::{choice_or, f64_or, req_str, usize_or};
use crate::common::{load_input_raster, parse_optional_output};
use crate::raster_stack::raster_like_multiband;
use crate::vector_common::write_or_store_layer;

pub struct CoregisterRastersTool;

impl Tool for CoregisterRastersTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "coregister_rasters",
            display_name: "Coregister Rasters",
            summary: "Automatically aligns a secondary raster to a reference by normalised cross-correlation of a tile grid, sub-pixel peak refinement, a least-squares transform fit with outlier trimming, and resampling (ArcGIS Apply Coregistration). This is the missing precondition for the catalog's whole InSAR chain: sar_coherence, flatten_interferogram and goldstein_phase_filter all require co-registered input and nothing produced it. georeference_raster_from_control_points needs hand-picked control points, image_correlation only scores similarity, and warp_raster needs a transform supplied to it.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "reference",
                    description: "Reference raster. The output is written on this raster's grid.",
                    required: true,
                },
                ToolParamSpec {
                    name: "secondary",
                    description: "Raster to align to the reference. All of its bands are resampled.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output resampled secondary raster on the reference grid. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tie_points",
                    description: "Output point layer of accepted tie points with their measured shift and correlation. Always produced; stored in memory when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band used for matching (default 1). A two-band input is matched on its complex magnitude unless this is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "transform",
                    description: "'translation', 'affine' (default), or 'polynomial2'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tile_size",
                    description: "Side of each matching tile in cells (default 32).",
                    required: false,
                },
                ToolParamSpec {
                    name: "grid_size",
                    description: "Number of tiles per axis (default 8, so 64 tiles).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_shift",
                    description: "Half-width of the integer search window in cells (default 8).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_correlation",
                    description: "Reject tie points whose peak normalised correlation is below this (default 0.3).",
                    required: false,
                },
                ToolParamSpec {
                    name: "outlier_sigma",
                    description: "Trim tie points whose fit residual exceeds this many standard deviations (default 3.0). Set to 0 to disable trimming.",
                    required: false,
                },
                ToolParamSpec {
                    name: "resample",
                    description: "'bilinear' (default) or 'nearest'.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "reference")?;
        req_str(args, "secondary")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let ref_path = req_str(args, "reference")?.to_string();
        let sec_path = req_str(args, "secondary")?.to_string();
        let prm = parse_params(args)?;
        let output = parse_optional_output(args, "output")?;
        let tie_out = parse_optional_output(args, "tie_points")?;

        let reference = load_input_raster(&ref_path)?;
        let secondary = load_input_raster(&sec_path)?;
        if reference.crs.epsg.is_some()
            && secondary.crs.epsg.is_some()
            && reference.crs.epsg != secondary.crs.epsg
        {
            return Err(ToolError::Validation(format!(
                "reference CRS (EPSG {:?}) differs from secondary (EPSG {:?}); reproject first",
                reference.crs.epsg, secondary.crs.epsg
            )));
        }

        let (rows, cols) = (reference.rows, reference.cols);
        let ref_grey = match_band(&reference, prm.band);
        let sec_grey = match_band(&secondary, prm.band);

        // Nominal reference-pixel -> secondary-pixel mapping from the two
        // geotransforms. The matcher then measures the residual misregistration
        // on top of this, so rasters with different origins still work.
        let nominal = NominalMap::new(&reference, &secondary);

        ctx.progress.info(&format!(
            "reference {rows}x{cols}, secondary {}x{}, {} tiles of {}, search +/-{}",
            secondary.rows,
            secondary.cols,
            prm.grid_size * prm.grid_size,
            prm.tile_size,
            prm.max_shift
        ));

        let ties = measure_tie_points(
            &ref_grey,
            rows,
            cols,
            &sec_grey,
            secondary.rows,
            secondary.cols,
            &nominal,
            &prm,
            ctx,
        );
        let measured = ties.len();
        if measured < prm.transform.min_points() {
            return Err(ToolError::Execution(format!(
                "only {measured} tie point(s) passed min_correlation {}; a {} fit needs {}. \
                 Try a larger tile_size, a larger max_shift, or a lower min_correlation.",
                prm.min_correlation,
                prm.transform.label(),
                prm.transform.min_points()
            )));
        }

        let fit = fit_transform(&ties, prm.transform, prm.outlier_sigma)?;
        ctx.progress.info(&format!(
            "{} fit on {}/{} tie points, RMSE {:.4} px",
            prm.transform.label(),
            fit.used.len(),
            measured,
            fit.rmse
        ));

        // Resample every band of the secondary onto the reference grid.
        let nodata = -9999.0_f64;
        let mut bands: Vec<Vec<f64>> = Vec::with_capacity(secondary.bands);
        for b in 0..secondary.bands {
            let src = read_band(&secondary, b as isize);
            let mut dst = vec![nodata; rows * cols];
            for r in 0..rows {
                for c in 0..cols {
                    let (sx, sy) = fit.apply(c as f64 + 0.5, r as f64 + 0.5);
                    if let Some(v) = sample(
                        &src,
                        secondary.rows,
                        secondary.cols,
                        sx,
                        sy,
                        prm.bilinear,
                    ) {
                        dst[r * cols + c] = v;
                    }
                }
            }
            bands.push(dst);
        }
        let warped = raster_like_multiband(&reference, &bands, nodata, DataType::F32)?;
        let out_path = crate::common::write_or_store_output(warped, output)?;

        let tie_layer = tie_point_layer(&reference, &ties, &fit);
        let tie_path = write_or_store_layer(tie_layer, tie_out)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("tie_points".to_string(), json!(tie_path));
        outputs.insert("transform".to_string(), json!(prm.transform.label()));
        outputs.insert("tie_points_measured".to_string(), json!(measured));
        outputs.insert("tie_points_used".to_string(), json!(fit.used.len()));
        outputs.insert("rmse_pixels".to_string(), json!(fit.rmse));
        outputs.insert("mean_shift_col".to_string(), json!(fit.mean_shift.0));
        outputs.insert("mean_shift_row".to_string(), json!(fit.mean_shift.1));
        outputs.insert("bands".to_string(), json!(secondary.bands));
        Ok(ToolRunResult { outputs })
    }
}

/// The nominal reference-pixel -> secondary-pixel map implied by the two
/// geotransforms (before any measured misregistration).
struct NominalMap {
    /// secondary_col = (ref_x - sec_x_min) / sec_cell_x
    sec_x_min: f64,
    sec_y_max: f64,
    sec_cell_x: f64,
    sec_cell_y: f64,
    ref_x_min: f64,
    ref_y_max: f64,
    ref_cell_x: f64,
    ref_cell_y: f64,
}

impl NominalMap {
    fn new(reference: &Raster, secondary: &Raster) -> Self {
        NominalMap {
            sec_x_min: secondary.x_min,
            sec_y_max: secondary.y_min + secondary.rows as f64 * secondary.cell_size_y,
            sec_cell_x: secondary.cell_size_x,
            sec_cell_y: secondary.cell_size_y,
            ref_x_min: reference.x_min,
            ref_y_max: reference.y_min + reference.rows as f64 * reference.cell_size_y,
            ref_cell_x: reference.cell_size_x,
            ref_cell_y: reference.cell_size_y,
        }
    }

    /// Maps a reference pixel centre to the nominal secondary pixel coordinate.
    fn map(&self, ref_col: f64, ref_row: f64) -> (f64, f64) {
        let x = self.ref_x_min + ref_col * self.ref_cell_x;
        let y = self.ref_y_max - ref_row * self.ref_cell_y;
        (
            (x - self.sec_x_min) / self.sec_cell_x,
            (self.sec_y_max - y) / self.sec_cell_y,
        )
    }
}

/// One accepted match: reference pixel centre, the secondary pixel it maps to,
/// and the correlation peak value.
struct TiePoint {
    ref_col: f64,
    ref_row: f64,
    sec_col: f64,
    sec_row: f64,
    correlation: f64,
}

/// The band used for matching, with no-data as NaN.
///
/// A two-band raster with no explicit band request is treated as complex I/Q
/// and matched on its magnitude — matching on I alone would correlate the
/// carrier rather than the scene.
fn match_band(r: &Raster, band: Option<isize>) -> Vec<f64> {
    let (rows, cols) = (r.rows, r.cols);
    let nd = r.nodata;
    let mut out = vec![f64::NAN; rows * cols];
    let complex = r.bands == 2 && band.is_none();
    let b = band.unwrap_or(0);
    for row in 0..rows {
        for col in 0..cols {
            let i = row * cols + col;
            let a = r.get(b, row as isize, col as isize);
            if a == nd || !a.is_finite() {
                continue;
            }
            if complex {
                let q = r.get(1, row as isize, col as isize);
                if q == nd || !q.is_finite() {
                    continue;
                }
                out[i] = (a * a + q * q).sqrt();
            } else {
                out[i] = a;
            }
        }
    }
    out
}

fn read_band(r: &Raster, band: isize) -> Vec<f64> {
    let (rows, cols) = (r.rows, r.cols);
    let nd = r.nodata;
    let mut out = vec![f64::NAN; rows * cols];
    for row in 0..rows {
        for col in 0..cols {
            let v = r.get(band, row as isize, col as isize);
            if v != nd && v.is_finite() {
                out[row * cols + col] = v;
            }
        }
    }
    out
}

/// Matches a grid of reference tiles against the secondary.
#[allow(clippy::too_many_arguments)]
fn measure_tie_points(
    refi: &[f64],
    ref_rows: usize,
    ref_cols: usize,
    sec: &[f64],
    sec_rows: usize,
    sec_cols: usize,
    nominal: &NominalMap,
    prm: &Params,
    ctx: &ToolContext,
) -> Vec<TiePoint> {
    let t = prm.tile_size;
    let g = prm.grid_size;
    let mut ties = Vec::new();
    let mut done = 0usize;

    for gi in 0..g {
        for gj in 0..g {
            done += 1;
            ctx.progress.progress(done as f64 / (g * g) as f64);

            // Tile origins span a g x g grid, inset by `max_shift` so every
            // tile can search its full window in both directions. Without the
            // inset the edge tiles can only search inward, their peak clamps to
            // zero shift, and the mean of the fit is dragged towards the
            // origin in proportion to how many tiles sit on the border.
            if ref_rows < t || ref_cols < t {
                continue;
            }
            let (lo_r, hi_r) = inset_span(ref_rows, t, prm.max_shift);
            let (lo_c, hi_c) = inset_span(ref_cols, t, prm.max_shift);
            let r0 = lo_r + if g == 1 { (hi_r - lo_r) / 2 } else { gi * (hi_r - lo_r) / (g - 1) };
            let c0 = lo_c + if g == 1 { (hi_c - lo_c) / 2 } else { gj * (hi_c - lo_c) / (g - 1) };

            // Reference tile, mean-centred; a flat tile carries no information
            // and would make the normalisation singular.
            let mut tile = Vec::with_capacity(t * t);
            let mut ok = true;
            for rr in 0..t {
                for cc in 0..t {
                    let v = refi[(r0 + rr) * ref_cols + c0 + cc];
                    if !v.is_finite() {
                        ok = false;
                        break;
                    }
                    tile.push(v);
                }
                if !ok {
                    break;
                }
            }
            if !ok {
                continue;
            }
            let (tile_mean, tile_norm) = mean_and_norm(&tile);
            if tile_norm <= 0.0 {
                continue;
            }

            // Nominal position of the tile's top-left corner in the secondary.
            let (nsc, nsr) = nominal.map(c0 as f64, r0 as f64);
            let base_c = nsc.round() as isize;
            let base_r = nsr.round() as isize;

            // Integer NCC search.
            let m = prm.max_shift as isize;
            let mut best = (f64::NEG_INFINITY, 0isize, 0isize);
            let mut surface = vec![f64::NEG_INFINITY; ((2 * m + 1) * (2 * m + 1)) as usize];
            for dr in -m..=m {
                for dc in -m..=m {
                    let sr = base_r + dr;
                    let sc = base_c + dc;
                    if sr < 0 || sc < 0 {
                        continue;
                    }
                    let (sr, sc) = (sr as usize, sc as usize);
                    if sr + t > sec_rows || sc + t > sec_cols {
                        continue;
                    }
                    let Some(score) =
                        ncc(&tile, tile_mean, tile_norm, sec, sec_rows, sec_cols, sr, sc, t)
                    else {
                        continue;
                    };
                    surface[((dr + m) * (2 * m + 1) + (dc + m)) as usize] = score;
                    if score > best.0 {
                        best = (score, dr, dc);
                    }
                }
            }
            if best.0 < prm.min_correlation {
                continue;
            }
            // A peak sitting on the edge of the search window is not a peak:
            // the true maximum may lie outside, and the parabolic refinement
            // has no neighbour on that side. Accepting these biases the fit
            // towards zero shift, so they are dropped instead.
            if best.1.abs() == m || best.2.abs() == m {
                continue;
            }

            // Sub-pixel refinement: parabola through the peak and its two
            // neighbours on each axis. Skipped on the search-window edge, where
            // one neighbour is missing.
            let at = |dr: isize, dc: isize| -> f64 {
                if dr < -m || dr > m || dc < -m || dc > m {
                    return f64::NEG_INFINITY;
                }
                surface[((dr + m) * (2 * m + 1) + (dc + m)) as usize]
            };
            let sub_c = parabolic(at(best.1, best.2 - 1), best.0, at(best.1, best.2 + 1));
            let sub_r = parabolic(at(best.1 - 1, best.2), best.0, at(best.1 + 1, best.2));

            ties.push(TiePoint {
                // Tile centre, in reference pixel coordinates.
                ref_col: c0 as f64 + t as f64 / 2.0,
                ref_row: r0 as f64 + t as f64 / 2.0,
                sec_col: (base_c + best.2) as f64 + sub_c + t as f64 / 2.0,
                sec_row: (base_r + best.1) as f64 + sub_r + t as f64 / 2.0,
                correlation: best.0,
            });
        }
    }
    ties
}

/// The inclusive range of tile origins along one axis, inset by `margin` so a
/// tile's search window stays inside the raster.
///
/// Falls back to the un-inset range when the raster is too small to afford the
/// margin; the edge-peak rejection then discards whatever could not be searched
/// properly.
fn inset_span(extent: usize, tile: usize, margin: usize) -> (usize, usize) {
    let full_hi = extent - tile;
    if full_hi < 2 * margin {
        return (0, full_hi);
    }
    (margin, full_hi - margin)
}

/// Mean and L2 norm of the mean-centred values.
fn mean_and_norm(v: &[f64]) -> (f64, f64) {
    let n = v.len() as f64;
    let mean = v.iter().sum::<f64>() / n;
    let norm = v.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>().sqrt();
    (mean, norm)
}

/// Normalised cross-correlation of a mean-centred tile against a window of the
/// secondary. Returns `None` if the window has no-data or is flat.
#[allow(clippy::too_many_arguments)]
fn ncc(
    tile: &[f64],
    tile_mean: f64,
    tile_norm: f64,
    sec: &[f64],
    _sec_rows: usize,
    sec_cols: usize,
    sr: usize,
    sc: usize,
    t: usize,
) -> Option<f64> {
    let mut sum = 0.0;
    let mut sum_sq = 0.0;
    for rr in 0..t {
        for cc in 0..t {
            let v = sec[(sr + rr) * sec_cols + sc + cc];
            if !v.is_finite() {
                return None;
            }
            sum += v;
            sum_sq += v * v;
        }
    }
    let n = (t * t) as f64;
    let mean = sum / n;
    let var = sum_sq - n * mean * mean;
    if var <= 0.0 {
        return None;
    }
    let norm = var.sqrt();

    let mut dot = 0.0;
    for rr in 0..t {
        for cc in 0..t {
            let a = tile[rr * t + cc] - tile_mean;
            let b = sec[(sr + rr) * sec_cols + sc + cc] - mean;
            dot += a * b;
        }
    }
    Some(dot / (tile_norm * norm))
}

/// Sub-pixel offset of a parabola through three samples, clamped to +/-1 cell.
fn parabolic(left: f64, peak: f64, right: f64) -> f64 {
    if !left.is_finite() || !right.is_finite() {
        return 0.0;
    }
    let denom = left - 2.0 * peak + right;
    if denom.abs() < 1e-12 {
        return 0.0;
    }
    (0.5 * (left - right) / denom).clamp(-1.0, 1.0)
}

/// Which spatial model maps reference pixels to secondary pixels.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum TransformKind {
    Translation,
    Affine,
    Polynomial2,
}

impl TransformKind {
    fn label(self) -> &'static str {
        match self {
            TransformKind::Translation => "translation",
            TransformKind::Affine => "affine",
            TransformKind::Polynomial2 => "polynomial2",
        }
    }

    /// Basis terms evaluated at a reference pixel coordinate.
    fn basis(self, x: f64, y: f64) -> Vec<f64> {
        match self {
            TransformKind::Translation => vec![1.0],
            TransformKind::Affine => vec![1.0, x, y],
            TransformKind::Polynomial2 => vec![1.0, x, y, x * x, x * y, y * y],
        }
    }

    fn n_terms(self) -> usize {
        self.basis(0.0, 0.0).len()
    }

    /// Fewest tie points that determine the model.
    fn min_points(self) -> usize {
        self.n_terms()
    }
}

/// A fitted misregistration model plus its quality.
struct Fit {
    kind: TransformKind,
    /// Coefficients for secondary column and row respectively.
    cx: Vec<f64>,
    cy: Vec<f64>,
    /// Indices of the tie points that survived trimming.
    used: Vec<usize>,
    rmse: f64,
    mean_shift: (f64, f64),
}

impl Fit {
    /// Maps a reference pixel coordinate to a secondary pixel coordinate.
    fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        let b = self.kind.basis(x, y);
        let mut sx = 0.0;
        let mut sy = 0.0;
        for (i, t) in b.iter().enumerate() {
            sx += self.cx[i] * t;
            sy += self.cy[i] * t;
        }
        // A pure translation still needs the identity part of the mapping: the
        // fit models only the constant offset, so the reference position must
        // be carried through.
        if self.kind == TransformKind::Translation {
            (x + sx, y + sy)
        } else {
            (sx, sy)
        }
    }
}

/// Least-squares fit with iterative outlier trimming.
fn fit_transform(
    ties: &[TiePoint],
    kind: TransformKind,
    outlier_sigma: f64,
) -> Result<Fit, ToolError> {
    let mut active: Vec<usize> = (0..ties.len()).collect();
    let mut cx;
    let mut cy;
    let mut rmse;

    // Three trimming passes: enough to shed gross blunders without eroding a
    // legitimately noisy set.
    let mut pass = 0;
    loop {
        let n_terms = kind.n_terms();
        let mut rows_a: Vec<Vec<f64>> = Vec::with_capacity(active.len());
        let mut bx: Vec<f64> = Vec::with_capacity(active.len());
        let mut by: Vec<f64> = Vec::with_capacity(active.len());
        for &i in &active {
            let t = &ties[i];
            rows_a.push(kind.basis(t.ref_col, t.ref_row));
            // Translation models the *offset*; the others model the absolute
            // secondary coordinate.
            if kind == TransformKind::Translation {
                bx.push(t.sec_col - t.ref_col);
                by.push(t.sec_row - t.ref_row);
            } else {
                bx.push(t.sec_col);
                by.push(t.sec_row);
            }
        }
        cx = solve_least_squares(&rows_a, &bx, n_terms)?;
        cy = solve_least_squares(&rows_a, &by, n_terms)?;

        // Residuals in pixels.
        let mut resid: Vec<f64> = Vec::with_capacity(active.len());
        for (k, &i) in active.iter().enumerate() {
            let px: f64 = rows_a[k].iter().zip(&cx).map(|(a, b)| a * b).sum();
            let py: f64 = rows_a[k].iter().zip(&cy).map(|(a, b)| a * b).sum();
            let (ex, ey) = if kind == TransformKind::Translation {
                (bx[k] - px, by[k] - py)
            } else {
                (ties[i].sec_col - px, ties[i].sec_row - py)
            };
            resid.push((ex * ex + ey * ey).sqrt());
        }
        rmse = (resid.iter().map(|r| r * r).sum::<f64>() / resid.len() as f64).sqrt();

        pass += 1;
        if outlier_sigma <= 0.0 || pass >= 3 || rmse <= 0.0 {
            break;
        }
        let cutoff = outlier_sigma * rmse;
        let keep: Vec<usize> = active
            .iter()
            .zip(&resid)
            .filter(|(_, &r)| r <= cutoff)
            .map(|(&i, _)| i)
            .collect();
        // Never trim below what the model needs, and stop when nothing moved.
        if keep.len() < kind.min_points() || keep.len() == active.len() {
            break;
        }
        active = keep;
    }

    let mean_shift = {
        let n = active.len() as f64;
        let sx: f64 = active.iter().map(|&i| ties[i].sec_col - ties[i].ref_col).sum();
        let sy: f64 = active.iter().map(|&i| ties[i].sec_row - ties[i].ref_row).sum();
        (sx / n, sy / n)
    };

    Ok(Fit {
        kind,
        cx,
        cy,
        used: active,
        rmse,
        mean_shift,
    })
}

/// Solves `A x = b` in the least-squares sense via the normal equations with
/// Gaussian elimination and partial pivoting.
fn solve_least_squares(a: &[Vec<f64>], b: &[f64], n: usize) -> Result<Vec<f64>, ToolError> {
    // Normal equations: (A^T A) x = A^T b.
    let mut ata = vec![vec![0.0f64; n]; n];
    let mut atb = vec![0.0f64; n];
    for (row, &bi) in a.iter().zip(b) {
        for i in 0..n {
            atb[i] += row[i] * bi;
            for j in 0..n {
                ata[i][j] += row[i] * row[j];
            }
        }
    }
    // Tiny Tikhonov term: a degenerate tie-point layout (all tiles on one line)
    // makes A^T A singular, and a slightly-biased answer beats an error.
    for (i, r) in ata.iter_mut().enumerate() {
        r[i] += 1e-9;
    }

    for col in 0..n {
        let pivot = (col..n)
            .max_by(|&r1, &r2| ata[r1][col].abs().total_cmp(&ata[r2][col].abs()))
            .unwrap();
        if ata[pivot][col].abs() < 1e-12 {
            return Err(ToolError::Execution(
                "tie points are degenerate; the transform is not determined".to_string(),
            ));
        }
        ata.swap(col, pivot);
        atb.swap(col, pivot);
        for r in col + 1..n {
            let f = ata[r][col] / ata[col][col];
            if f == 0.0 {
                continue;
            }
            for c in col..n {
                ata[r][c] -= f * ata[col][c];
            }
            atb[r] -= f * atb[col];
        }
    }
    let mut x = vec![0.0f64; n];
    for i in (0..n).rev() {
        let mut s = atb[i];
        for j in i + 1..n {
            s -= ata[i][j] * x[j];
        }
        x[i] = s / ata[i][i];
    }
    Ok(x)
}

/// Samples a band at a fractional pixel coordinate.
fn sample(src: &[f64], rows: usize, cols: usize, x: f64, y: f64, bilinear: bool) -> Option<f64> {
    if !x.is_finite() || !y.is_finite() {
        return None;
    }
    if !bilinear {
        let c = x.floor() as isize;
        let r = y.floor() as isize;
        if r < 0 || c < 0 || r as usize >= rows || c as usize >= cols {
            return None;
        }
        let v = src[r as usize * cols + c as usize];
        return v.is_finite().then_some(v);
    }
    // Pixel centres sit at (col + 0.5); shift so the interpolation weights are
    // measured between centres rather than between corners.
    let fx = x - 0.5;
    let fy = y - 0.5;
    let c0 = fx.floor();
    let r0 = fy.floor();
    let tx = fx - c0;
    let ty = fy - r0;
    let (c0, r0) = (c0 as isize, r0 as isize);
    let mut acc = 0.0;
    let mut wsum = 0.0;
    for (dr, wy) in [(0isize, 1.0 - ty), (1, ty)] {
        for (dc, wx) in [(0isize, 1.0 - tx), (1, tx)] {
            let r = r0 + dr;
            let c = c0 + dc;
            if r < 0 || c < 0 || r as usize >= rows || c as usize >= cols {
                continue;
            }
            let v = src[r as usize * cols + c as usize];
            if !v.is_finite() {
                continue;
            }
            let w = wx * wy;
            acc += w * v;
            wsum += w;
        }
    }
    // Require most of the kernel to be present, so edge cells do not get a
    // value extrapolated from a single corner.
    (wsum > 0.5).then(|| acc / wsum)
}

/// Builds the tie-point output layer in the reference's CRS.
fn tie_point_layer(reference: &Raster, ties: &[TiePoint], fit: &Fit) -> Layer {
    let mut layer = Layer::new("coregister_tie_points");
    layer.geom_type = Some(GeometryType::Point);
    if let Some(e) = reference.crs.epsg {
        layer = layer.with_crs_epsg(e);
    }
    layer.add_field(FieldDef::new("id", FieldType::Integer));
    layer.add_field(FieldDef::new("ref_col", FieldType::Float));
    layer.add_field(FieldDef::new("ref_row", FieldType::Float));
    layer.add_field(FieldDef::new("shift_col", FieldType::Float));
    layer.add_field(FieldDef::new("shift_row", FieldType::Float));
    layer.add_field(FieldDef::new("correlation", FieldType::Float));
    layer.add_field(FieldDef::new("used", FieldType::Boolean));
    layer.add_field(FieldDef::new("residual", FieldType::Float));

    let y_max = reference.y_min + reference.rows as f64 * reference.cell_size_y;
    for (i, t) in ties.iter().enumerate() {
        let x = reference.x_min + t.ref_col * reference.cell_size_x;
        let y = y_max - t.ref_row * reference.cell_size_y;
        let (px, py) = fit.apply(t.ref_col, t.ref_row);
        let residual = ((t.sec_col - px).powi(2) + (t.sec_row - py).powi(2)).sqrt();

        let mut f = Feature::with_geometry(
            i as u64,
            Geometry::Point(Coord::xy(x, y)),
            layer.schema.len(),
        );
        f.set_by_index(0, FieldValue::Integer(i as i64));
        f.set_by_index(1, FieldValue::Float(t.ref_col));
        f.set_by_index(2, FieldValue::Float(t.ref_row));
        f.set_by_index(3, FieldValue::Float(t.sec_col - t.ref_col));
        f.set_by_index(4, FieldValue::Float(t.sec_row - t.ref_row));
        f.set_by_index(5, FieldValue::Float(t.correlation));
        f.set_by_index(6, FieldValue::Boolean(fit.used.contains(&i)));
        f.set_by_index(7, FieldValue::Float(residual));
        layer.push(f);
    }
    layer
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    band: Option<isize>,
    transform: TransformKind,
    tile_size: usize,
    grid_size: usize,
    max_shift: usize,
    min_correlation: f64,
    outlier_sigma: f64,
    bilinear: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let band = match crate::args_common::opt_usize(args, "band")? {
        None => None,
        Some(b) => Some(b as isize),
    };
    let transform = match choice_or(
        args,
        "transform",
        &["translation", "affine", "polynomial2"],
        "affine",
    )? {
        "translation" => TransformKind::Translation,
        "polynomial2" => TransformKind::Polynomial2,
        _ => TransformKind::Affine,
    };
    let tile_size = usize_or(args, "tile_size", 32)?;
    if tile_size < 4 {
        return Err(ToolError::Validation(format!(
            "'tile_size' must be at least 4, got {tile_size}"
        )));
    }
    let grid_size = usize_or(args, "grid_size", 8)?;
    if grid_size == 0 {
        return Err(ToolError::Validation(
            "'grid_size' must be at least 1".to_string(),
        ));
    }
    let max_shift = usize_or(args, "max_shift", 8)?;
    if max_shift == 0 {
        return Err(ToolError::Validation(
            "'max_shift' must be at least 1".to_string(),
        ));
    }
    let min_correlation = f64_or(args, "min_correlation", 0.3)?;
    if !(-1.0..=1.0).contains(&min_correlation) {
        return Err(ToolError::Validation(format!(
            "'min_correlation' must be in [-1, 1], got {min_correlation}"
        )));
    }
    let outlier_sigma = f64_or(args, "outlier_sigma", 3.0)?;
    if outlier_sigma < 0.0 {
        return Err(ToolError::Validation(
            "'outlier_sigma' must not be negative (0 disables trimming)".to_string(),
        ));
    }
    let bilinear = choice_or(args, "resample", &["bilinear", "nearest"], "bilinear")? == "bilinear";

    Ok(Params {
        band,
        transform,
        tile_size,
        grid_size,
        max_shift,
        min_correlation,
        outlier_sigma,
        bilinear,
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

    /// A smooth, texture-rich scene — correlation needs structure to lock onto.
    fn scene(x: f64, y: f64) -> f64 {
        50.0 + 20.0 * (x / 7.0).sin() * (y / 9.0).cos()
            + 10.0 * (x / 23.0 + y / 17.0).sin()
            + 5.0 * (x / 3.5).cos()
    }

    fn raster_from(cols: usize, rows: usize, bands: usize, f: impl Fn(usize, usize, usize) -> f64) -> String {
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
                    r.set(b as isize, row as isize, col as isize, f(b, row, col))
                        .unwrap();
                }
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: Value) -> (Raster, BTreeMap<String, Value>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CoregisterRastersTool.run(&args, &ctx()).unwrap();
        let raster = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (raster, out.outputs)
    }

    /// A known integer shift must be measured to well under a pixel and undone.
    #[test]
    fn recovers_a_known_integer_shift() {
        let (rows, cols) = (96, 96);
        let reference = raster_from(cols, rows, 1, |_, r, c| scene(c as f64, r as f64));
        // Secondary is the scene shifted by (+3 cols, -2 rows): the value at
        // secondary (r, c) is the reference value at (r + 2, c - 3).
        let secondary =
            raster_from(cols, rows, 1, |_, r, c| scene(c as f64 - 3.0, r as f64 + 2.0));

        let (out, outputs) = run(json!({
            "reference": reference, "secondary": secondary,
            "tile_size": 24, "grid_size": 5, "max_shift": 6
        }));

        let dc = outputs["mean_shift_col"].as_f64().unwrap();
        let dr = outputs["mean_shift_row"].as_f64().unwrap();
        assert!(
            (dc - 3.0).abs() < 0.3,
            "column shift {dc} should be about +3"
        );
        assert!((dr + 2.0).abs() < 0.3, "row shift {dr} should be about -2");
        assert!(
            outputs["rmse_pixels"].as_f64().unwrap() < 0.5,
            "fit RMSE too high: {}",
            outputs["rmse_pixels"]
        );

        // The warped secondary must now agree with the reference in the
        // interior (edges lose data that shifted in from outside).
        let mut worst: f64 = 0.0;
        for r in 12..rows - 12 {
            for c in 12..cols - 12 {
                let got = out.get(0, r as isize, c as isize);
                worst = worst.max((got - scene(c as f64, r as f64)).abs());
            }
        }
        assert!(worst < 1.5, "aligned image still differs by {worst}");
    }

    /// Sub-pixel shifts are what actually matter for InSAR, and the parabolic
    /// refinement must resolve them.
    #[test]
    fn recovers_a_subpixel_shift() {
        let (rows, cols) = (80, 80);
        let reference = raster_from(cols, rows, 1, |_, r, c| scene(c as f64, r as f64));
        let secondary =
            raster_from(cols, rows, 1, |_, r, c| scene(c as f64 - 1.4, r as f64 - 0.6));

        let (_, outputs) = run(json!({
            "reference": reference, "secondary": secondary,
            "tile_size": 24, "grid_size": 4, "max_shift": 5
        }));
        let dc = outputs["mean_shift_col"].as_f64().unwrap();
        let dr = outputs["mean_shift_row"].as_f64().unwrap();
        assert!(
            (dc - 1.4).abs() < 0.35,
            "sub-pixel column shift {dc} should be about +1.4"
        );
        assert!(
            (dr - 0.6).abs() < 0.35,
            "sub-pixel row shift {dr} should be about +0.6"
        );
    }

    /// Already-aligned input must be left alone rather than nudged.
    #[test]
    fn identity_when_already_aligned() {
        let (rows, cols) = (72, 72);
        let reference = raster_from(cols, rows, 1, |_, r, c| scene(c as f64, r as f64));
        let secondary = raster_from(cols, rows, 1, |_, r, c| scene(c as f64, r as f64));
        let (_, outputs) = run(json!({
            "reference": reference, "secondary": secondary,
            "tile_size": 24, "grid_size": 4, "max_shift": 4
        }));
        assert!(outputs["mean_shift_col"].as_f64().unwrap().abs() < 0.15);
        assert!(outputs["mean_shift_row"].as_f64().unwrap().abs() < 0.15);
    }

    /// Every band of a complex secondary is warped, not just the matched one.
    #[test]
    fn warps_all_bands_of_a_complex_raster() {
        let (rows, cols) = (72, 72);
        let reference = raster_from(cols, rows, 2, |b, r, c| {
            let p = scene(c as f64, r as f64);
            if b == 0 {
                p
            } else {
                p * 0.5
            }
        });
        let secondary = raster_from(cols, rows, 2, |b, r, c| {
            let p = scene(c as f64 - 2.0, r as f64);
            if b == 0 {
                p
            } else {
                p * 0.5
            }
        });
        let (out, outputs) = run(json!({
            "reference": reference, "secondary": secondary,
            "tile_size": 24, "grid_size": 4, "max_shift": 5
        }));
        assert_eq!(out.bands, 2, "both bands must survive");
        assert_eq!(outputs["bands"].as_u64().unwrap(), 2);
        // Band 1 is band 0 halved, and must be aligned identically.
        let b0 = out.get(0, 36, 36);
        let b1 = out.get(1, 36, 36);
        assert!(
            (b1 - b0 * 0.5).abs() < 1e-3,
            "band 1 ({b1}) is not band 0 ({b0}) halved — bands warped inconsistently"
        );
    }

    /// The tie-point layer records the measurements and which ones were used.
    #[test]
    fn emits_tie_points_even_without_a_path() {
        let (rows, cols) = (64, 64);
        let reference = raster_from(cols, rows, 1, |_, r, c| scene(c as f64, r as f64));
        let secondary = raster_from(cols, rows, 1, |_, r, c| scene(c as f64 - 2.0, r as f64));
        let args: ToolArgs = serde_json::from_value(json!({
            "reference": reference, "secondary": secondary,
            "tile_size": 20, "grid_size": 3, "max_shift": 4
        }))
        .unwrap();
        let out = CoregisterRastersTool.run(&args, &ctx()).unwrap();
        let path = out.outputs["tie_points"].as_str().unwrap();
        assert!(!path.is_empty(), "tie points must be produced without a path");
        let layer = crate::vector_common::load_input_layer(path).unwrap();
        assert_eq!(layer.len(), 9, "3x3 grid should yield 9 tie points");
        let used_idx = layer.schema.field_index("used").unwrap();
        let used = layer
            .iter()
            .filter(|f| matches!(f.attributes[used_idx], FieldValue::Boolean(true)))
            .count();
        assert!(used >= 3, "expected most tie points to be used, got {used}");
    }

    /// A featureless scene gives no reliable matches and must fail loudly
    /// rather than return a bogus transform.
    #[test]
    fn errors_when_no_tie_points_pass() {
        let (rows, cols) = (64, 64);
        let flat = raster_from(cols, rows, 1, |_, _, _| 7.0);
        let flat2 = raster_from(cols, rows, 1, |_, _, _| 7.0);
        let args: ToolArgs = serde_json::from_value(json!({
            "reference": flat, "secondary": flat2, "tile_size": 16, "grid_size": 3
        }))
        .unwrap();
        let err = CoregisterRastersTool.run(&args, &ctx()).unwrap_err();
        assert!(
            format!("{err:?}").contains("tie point"),
            "expected a tie-point failure, got {err:?}"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CoregisterRastersTool.validate(&args)
        };
        assert!(bad(json!({"secondary": "b.tif"})).is_err());
        assert!(bad(json!({"reference": "a.tif"})).is_err());
        let base = |extra: Value| {
            let mut m = serde_json::Map::new();
            m.insert("reference".into(), json!("a.tif"));
            m.insert("secondary".into(), json!("b.tif"));
            if let Value::Object(o) = extra {
                m.extend(o);
            }
            Value::Object(m)
        };
        assert!(bad(base(json!({"transform": "spline"}))).is_err());
        assert!(bad(base(json!({"tile_size": 2}))).is_err());
        assert!(bad(base(json!({"max_shift": 0}))).is_err());
        assert!(bad(base(json!({"min_correlation": 2.0}))).is_err());
        assert!(bad(base(json!({"outlier_sigma": -1.0}))).is_err());
        assert!(bad(base(json!({"resample": "cubic"}))).is_err());
        assert!(bad(base(json!({"transform": "polynomial2"}))).is_ok());
    }
}
