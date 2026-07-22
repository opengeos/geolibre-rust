//! GeoLibre tool: `empirical_bayesian_kriging` — local-subset kriging with
//! simulated semivariograms.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Empirical Bayesian Kriging* (EBK).
//! The bundled suite already has classical kriging (`ordinary_kriging`,
//! `simple_kriging`, `universal_kriging`, `ordinary_cokriging`, `local_kriging`,
//! `fit_variogram`/`estimate_variogram`) but every one of them fits a single,
//! manually-accepted semivariogram — either global or over one moving window —
//! and treats it as known-with-certainty. EBK's distinguishing idea is
//! different: split the data into many small, overlapping local neighbourhoods,
//! fit a semivariogram *independently in each one*, and instead of trusting
//! that single fit, simulate several plausible semivariograms around it so the
//! output prediction-standard-error surface reflects genuine model uncertainty,
//! not just interpolation uncertainty conditional on one variogram.
//!
//! ## Method
//! 1. Deduplicate coincident input points (bit-exact coordinates; merge by
//!    mean of the value field) so the kriging system can never be singular
//!    from a zero-distance pair.
//! 2. Optionally `log_empirical`-transform the value field (requires strictly
//!    positive values); kriging then runs in log space and predictions/errors
//!    are back-transformed with the standard lognormal mean/variance formulas.
//! 3. Tile the point extent into a grid of overlapping local subsets sized so
//!    each holds roughly `subset_size` points (`overlap` pads every tile by
//!    that fraction of its width/height on each side, so subsets share data
//!    near their edges — no vendored kdtree crate is wired into this crate, so
//!    tiling is done directly over the point extent rather than through a
//!    tree).
//! 4. Per subset: fit an empirical semivariogram (distance-binned, pair-count
//!    weighted) by weighted least squares — closed-form 2×2 normal equations
//!    solved with a hand-rolled Cholesky solve for `power`/`linear`, plus a
//!    golden-section search over the range for `exponential`. Draw
//!    `simulations` perturbed semivariograms from a seeded deterministic RNG
//!    (inline splitmix64 + Box–Muller), log-normal jitter scaled by the fit's
//!    relative WLS residual — a bounded, reproducible stand-in for EBK's full
//!    parametric semivariogram bootstrap (documented cut below).
//! 5. Per output cell: krige against the nearest few subsets (ordinary
//!    kriging with a Lagrange multiplier, solved by a hand-rolled Gauss–Jordan
//!    matrix inverse — reused across every cell a subset/simulation touches),
//!    average each subset's `simulations` draws (law-of-total-variance:
//!    mean kriging variance + variance of the simulated predictions), then
//!    blend subsets by inverse-squared distance to their centroid and add the
//!    between-subset spread to the final standard error.
//!
//! Predicting exactly at an input point returns that point's value with zero
//! error (kriging's exact-interpolator property with zero nugget at the query
//! location), independent of simulation noise.
//!
//! ## Deliberate v1 scope cuts (documented for reviewers)
//! - Semivariogram simulation is log-normal parameter jitter around the WLS
//!   fit (seeded, deterministic), not EBK's full parametric semivariogram
//!   bootstrap (simulate synthetic data under the fitted model, refit, repeat).
//! - `power`'s exponent is fixed at 1.5 rather than fit; `linear`/`power`/
//!   `exponential` cover the ArcGIS default family, `thin_plate_spline` was
//!   dropped from the recommendation's proposed enum.
//! - Local subsets come from tiling the point extent, not a kdtree (this crate
//!   has no kdtree dependency); this is geometrically equivalent for roughly
//!   uniform point layouts but less adaptive under heavy clustering.
//! - Cell blending uses the nearest few subsets by centroid distance rather
//!   than a full nonstationary moving-window recombination.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
use wbvector::Geometry;

use crate::common::{parse_optional_output, write_or_store_output};
use crate::vector_common::load_input_layer;

pub struct EmpiricalBayesianKrigingTool;

impl Tool for EmpiricalBayesianKrigingTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "empirical_bayesian_kriging",
            display_name: "Empirical Bayesian Kriging",
            summary: "Interpolate a point field via many overlapping local subsets, each fitting its own semivariogram whose uncertainty is captured by simulating perturbed semivariograms from a seeded RNG, then blending subset kriging predictions/variances by distance (like ArcGIS's Empirical Bayesian Kriging) — distinct from the bundled ordinary/universal/simple/co-kriging, which all fit one manually-accepted global (or single moving-window) semivariogram.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "field",
                    description: "Numeric field on the input layer to interpolate.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output prediction raster. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Output cell size in CRS units. Default: sized so the longer grid axis is ~150 cells.",
                    required: false,
                },
                ToolParamSpec {
                    name: "subset_size",
                    description: "Target number of points per local subset (default 40, minimum 4).",
                    required: false,
                },
                ToolParamSpec {
                    name: "overlap",
                    description: "Fraction of a subset tile's width/height that neighbouring tiles overlap by (default 0.25).",
                    required: false,
                },
                ToolParamSpec {
                    name: "simulations",
                    description: "Number of simulated semivariograms drawn per subset to represent model uncertainty (default 5, max 50).",
                    required: false,
                },
                ToolParamSpec {
                    name: "semivariogram",
                    description: "Semivariogram model family: 'power' (default), 'linear', or 'exponential'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "transform",
                    description: "Data transform before kriging: 'none' (default) or 'log_empirical' (natural log; requires field > 0 everywhere).",
                    required: false,
                },
                ToolParamSpec {
                    name: "error_output",
                    description: "Optional output prediction standard-error raster.",
                    required: false,
                },
                ToolParamSpec {
                    name: "seed",
                    description: "Seed for the deterministic semivariogram simulation RNG (default 42).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "field")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let field = require_str(args, "field")?;
        let output = parse_optional_output(args, "output")?;
        let error_output = parse_optional_output(args, "error_output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let field_idx = layer
            .schema
            .field_index(field)
            .ok_or_else(|| ToolError::Validation(format!("field '{field}' not found")))?;
        let epsg = layer.crs_epsg();

        let mut points: Vec<Pt> = Vec::new();
        for feat in layer.iter() {
            let Some((x, y)) = feat.geometry.as_ref().and_then(point_xy) else {
                continue;
            };
            let Some(v) = feat.attributes.get(field_idx).and_then(|f| f.as_f64()) else {
                continue;
            };
            if !x.is_finite() || !y.is_finite() || !v.is_finite() {
                continue;
            }
            points.push(Pt { x, y, v });
        }
        if points.len() < 8 {
            return Err(ToolError::Execution(format!(
                "need at least 8 valid point(s) with a numeric '{field}' value, found {}",
                points.len()
            )));
        }

        points = dedupe_points(points);

        if matches!(prm.transform, Transform::LogEmpirical) {
            if points.iter().any(|p| p.v <= 0.0) {
                return Err(ToolError::Execution(
                    "'transform' = log_empirical requires every value in 'field' to be > 0"
                        .to_string(),
                ));
            }
            for p in &mut points {
                p.v = p.v.ln();
            }
        }

        // ── Output grid from the point extent ────────────────────────────────
        let (mut x_min, mut y_min, mut x_max, mut y_max) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for p in &points {
            x_min = x_min.min(p.x);
            x_max = x_max.max(p.x);
            y_min = y_min.min(p.y);
            y_max = y_max.max(p.y);
        }
        let span_x = (x_max - x_min).max(1e-9);
        let span_y = (y_max - y_min).max(1e-9);
        let cell_size = prm
            .cell_size
            .unwrap_or((span_x.max(span_y) / 150.0).max(1e-9));
        let margin = cell_size;
        let gx_min = x_min - margin;
        let gy_min = y_min - margin;
        let cols = (((x_max + margin - gx_min) / cell_size).ceil() as usize).max(1);
        let rows = (((y_max + margin - gy_min) / cell_size).ceil() as usize).max(1);
        if rows.saturating_mul(cols) > 4_000_000 {
            return Err(ToolError::Validation(format!(
                "requested grid is {rows}x{cols} cells (> 4M); increase cell_size"
            )));
        }
        let gy_max = gy_min + rows as f64 * cell_size;

        // ── Build overlapping local subsets and fit + simulate per subset ───
        let subsets = build_subsets(&points, prm.subset_size, prm.overlap);
        if subsets.is_empty() {
            return Err(ToolError::Execution(
                "could not form any local subset (need >= 4 points per subset); try a larger 'subset_size' or 'overlap'".to_string(),
            ));
        }
        ctx.progress.info(&format!(
            "{} point(s) -> {} local subset(s), fitting + simulating semivariograms",
            points.len(),
            subsets.len()
        ));

        let prepared: Vec<PreparedSubset> = subsets
            .iter()
            .enumerate()
            .map(|(si, s)| prepare_subset(&points, s, prm.kind, prm.simulations, prm.seed, si))
            .collect();

        ctx.progress
            .info(&format!("interpolating a {rows}x{cols} prediction raster"));

        // ── Krige every cell against its nearest local subsets ──────────────
        let n = rows * cols;
        let mut pred = vec![0.0f64; n];
        let mut se = vec![0.0f64; n];
        const NEAREST_M: usize = 4;
        for r in 0..rows {
            let cy = gy_max - (r as f64 + 0.5) * cell_size;
            for c in 0..cols {
                let cx = gx_min + (c as f64 + 0.5) * cell_size;
                let (p, s) = blend_cell(&prepared, cx, cy, NEAREST_M);
                pred[r * cols + c] = p;
                se[r * cols + c] = s;
            }
            if rows > 0 {
                ctx.progress.progress((r as f64 + 1.0) / rows as f64);
            }
        }

        // ── Back-transform if logged, then write outputs ─────────────────────
        let nodata = f64::NAN;
        if matches!(prm.transform, Transform::LogEmpirical) {
            for i in 0..n {
                let mu = pred[i];
                let var_log = (se[i] * se[i]).max(0.0);
                pred[i] = (mu + var_log / 2.0).exp();
                let var_orig = (var_log.exp() - 1.0) * (2.0 * mu + var_log).exp();
                se[i] = var_orig.max(0.0).sqrt();
            }
        }

        let raster = build_raster(&pred, rows, cols, gx_min, gy_min, cell_size, nodata, epsg)?;
        let out_path = write_or_store_output(raster, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        if let Some(p) = error_output {
            let se_raster = build_raster(&se, rows, cols, gx_min, gy_min, cell_size, nodata, epsg)?;
            outputs.insert(
                "error_output".to_string(),
                json!(write_or_store_output(se_raster, Some(p))?),
            );
        }
        outputs.insert("point_count".to_string(), json!(points.len()));
        outputs.insert("subset_count".to_string(), json!(subsets.len()));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        outputs.insert("cell_size".to_string(), json!(cell_size));
        Ok(ToolRunResult { outputs })
    }
}

// ── Point I/O ─────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
struct Pt {
    x: f64,
    y: f64,
    v: f64,
}

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        _ => None,
    }
}

/// Merges bit-exact coincident points (same x,y) by averaging their value, in a
/// stable coordinate-sorted pass so the result — and everything downstream —
/// is independent of the input feature order and of hash-map iteration order
/// (required for the determinism guarantee: same seed -> identical raster).
fn dedupe_points(mut pts: Vec<Pt>) -> Vec<Pt> {
    pts.sort_by(|a, b| (a.x.to_bits(), a.y.to_bits()).cmp(&(b.x.to_bits(), b.y.to_bits())));
    let mut out: Vec<Pt> = Vec::with_capacity(pts.len());
    let mut i = 0;
    while i < pts.len() {
        let mut j = i + 1;
        let mut sum = pts[i].v;
        let mut cnt = 1usize;
        while j < pts.len() && pts[j].x == pts[i].x && pts[j].y == pts[i].y {
            sum += pts[j].v;
            cnt += 1;
            j += 1;
        }
        out.push(Pt {
            x: pts[i].x,
            y: pts[i].y,
            v: sum / cnt as f64,
        });
        i = j;
    }
    out
}

// ── Local subsets ────────────────────────────────────────────────────────────

struct SubsetSpan {
    idxs: Vec<usize>,
    cx: f64,
    cy: f64,
}

/// Tiles the point extent into overlapping local subsets. Each tile targets
/// `subset_size` points; `overlap` pads every tile's bounding box by that
/// fraction of its width/height on each side so neighbouring subsets share
/// points near their shared edge. Tiles with fewer than 4 points are dropped
/// (not enough to fit a semivariogram); coverage of any dropped tile's area is
/// still provided by the "nearest subsets" blending at prediction time.
fn build_subsets(points: &[Pt], subset_size: usize, overlap: f64) -> Vec<SubsetSpan> {
    let (mut x_min, mut y_min, mut x_max, mut y_max) = (
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );
    for p in points {
        x_min = x_min.min(p.x);
        x_max = x_max.max(p.x);
        y_min = y_min.min(p.y);
        y_max = y_max.max(p.y);
    }
    let span_x = (x_max - x_min).max(1e-9);
    let span_y = (y_max - y_min).max(1e-9);
    let target = subset_size.max(4);
    let nt = ((points.len() as f64 / target as f64).sqrt().ceil() as usize).max(1);
    let tile_w = span_x / nt as f64;
    let tile_h = span_y / nt as f64;
    let pad_x = overlap * tile_w;
    let pad_y = overlap * tile_h;
    const MAX_SUBSET_POINTS: usize = 200;

    let mut out = Vec::new();
    for ti in 0..nt {
        let tx0 = x_min + ti as f64 * tile_w;
        let tx1 = if ti + 1 == nt {
            x_max
        } else {
            x_min + (ti + 1) as f64 * tile_w
        };
        let bx0 = tx0 - pad_x;
        let bx1 = tx1 + pad_x;
        for tj in 0..nt {
            let ty0 = y_min + tj as f64 * tile_h;
            let ty1 = if tj + 1 == nt {
                y_max
            } else {
                y_min + (tj + 1) as f64 * tile_h
            };
            let by0 = ty0 - pad_y;
            let by1 = ty1 + pad_y;

            let mut idxs: Vec<usize> = (0..points.len())
                .filter(|&i| {
                    let p = points[i];
                    p.x >= bx0 && p.x <= bx1 && p.y >= by0 && p.y <= by1
                })
                .collect();
            if idxs.len() < 4 {
                continue;
            }
            let tcx = (tx0 + tx1) / 2.0;
            let tcy = (ty0 + ty1) / 2.0;
            if idxs.len() > MAX_SUBSET_POINTS {
                idxs.sort_by(|&a, &b| {
                    let da = (points[a].x - tcx).powi(2) + (points[a].y - tcy).powi(2);
                    let db = (points[b].x - tcx).powi(2) + (points[b].y - tcy).powi(2);
                    da.partial_cmp(&db).unwrap()
                });
                idxs.truncate(MAX_SUBSET_POINTS);
            }
            let (mut cx, mut cy) = (0.0, 0.0);
            for &i in &idxs {
                cx += points[i].x;
                cy += points[i].y;
            }
            let n = idxs.len() as f64;
            out.push(SubsetSpan {
                idxs,
                cx: cx / n,
                cy: cy / n,
            });
        }
    }
    out
}

// ── Semivariogram model ──────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Power,
    Linear,
    Exponential,
}

#[derive(Clone, Copy)]
struct Variogram {
    kind: Kind,
    nugget: f64,
    scale: f64,
    range: f64,
}

/// Semivariance at separation `h`. By convention gamma(0) = 0 (the nugget is a
/// discontinuity approached only as h -> 0+, not the value *at* h = 0), which
/// keeps every diagonal of the kriging matrix at 0.
fn gamma_h(v: &Variogram, h: f64) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    match v.kind {
        Kind::Linear => v.nugget + v.scale * h,
        Kind::Power => v.nugget + v.scale * h.powf(1.5),
        Kind::Exponential => v.nugget + v.scale * (1.0 - (-h / v.range.max(1e-9)).exp()),
    }
}

fn dist(a: &Pt, x: f64, y: f64) -> f64 {
    ((a.x - x).powi(2) + (a.y - y).powi(2)).sqrt()
}

/// Weighted least-squares fit of `nugget + scale * predictor(h)` to the binned
/// empirical semivariogram, via the closed-form 2x2 normal equations solved by
/// [`cholesky_solve`]. Returns `(nugget, scale, wls_relative_rmse)`.
fn wls_fit_linear(bins: &[(f64, f64, f64)], predictor: impl Fn(f64) -> f64) -> (f64, f64, f64) {
    let (mut s_w, mut s_wx, mut s_wxx, mut s_wy, mut s_wxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for &(h, g, w) in bins {
        let x = predictor(h);
        s_w += w;
        s_wx += w * x;
        s_wxx += w * x * x;
        s_wy += w * g;
        s_wxy += w * x * g;
    }
    // Tiny ridge keeps the 2x2 Gram matrix strictly SPD even for degenerate
    // (e.g. single-bin) inputs.
    let ridge = 1e-9 * s_w.max(1.0);
    let a = [s_w + ridge, s_wx, s_wx, s_wxx + ridge];
    let b = [s_wy, s_wxy];
    let (nugget, scale) = match cholesky_solve(&a, &b, 2) {
        Some(x) => (x[0].max(0.0), x[1].max(1e-12)),
        None => (0.0, 1e-6),
    };
    let mut sse = 0.0;
    let mut sw = 0.0;
    for &(h, g, w) in bins {
        let x = predictor(h);
        let r = g - (nugget + scale * x);
        sse += w * r * r;
        sw += w * g * g;
    }
    let rel_rmse = if sw > 0.0 { (sse / sw).sqrt() } else { 0.2 };
    (nugget, scale, rel_rmse)
}

/// Fits a semivariogram to a subset's points: bins the empirical semivariogram
/// by distance (pair-count weighted), then WLS-fits the chosen model family.
/// Returns the fit and its relative WLS residual (used to scale simulation
/// jitter — a poorly-fitting model gets more spread across its simulations).
fn fit_variogram(pts: &[Pt], kind: Kind) -> (Variogram, f64) {
    let n = pts.len();
    let mut pairs: Vec<(f64, f64)> = Vec::with_capacity(n * (n - 1) / 2);
    let mut max_h = 0.0f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let h = dist(&pts[i], pts[j].x, pts[j].y);
            if h <= 0.0 {
                continue;
            }
            let g = 0.5 * (pts[i].v - pts[j].v).powi(2);
            pairs.push((h, g));
            max_h = max_h.max(h);
        }
    }
    if pairs.is_empty() || max_h <= 0.0 {
        return (
            Variogram {
                kind,
                nugget: 0.0,
                scale: 1e-6,
                range: 1.0,
            },
            0.2,
        );
    }
    let cutoff = max_h * 0.7;
    let nbins = (pairs.len() / 30).clamp(4, 12);
    let bin_w = (cutoff / nbins as f64).max(1e-12);
    let mut sum_h = vec![0.0; nbins];
    let mut sum_g = vec![0.0; nbins];
    let mut cnt = vec![0.0; nbins];
    for &(h, g) in &pairs {
        if h > cutoff {
            continue;
        }
        let bi = ((h / bin_w) as usize).min(nbins - 1);
        sum_h[bi] += h;
        sum_g[bi] += g;
        cnt[bi] += 1.0;
    }
    let mut bins: Vec<(f64, f64, f64)> = (0..nbins)
        .filter(|&i| cnt[i] > 0.0)
        .map(|i| (sum_h[i] / cnt[i], sum_g[i] / cnt[i], cnt[i]))
        .collect();
    if bins.is_empty() {
        bins = pairs.iter().map(|&(h, g)| (h, g, 1.0)).collect();
    }

    match kind {
        Kind::Linear => {
            let (nugget, scale, rmse) = wls_fit_linear(&bins, |h| h);
            (
                Variogram {
                    kind,
                    nugget,
                    scale,
                    range: max_h,
                },
                rmse,
            )
        }
        Kind::Power => {
            let (nugget, scale, rmse) = wls_fit_linear(&bins, |h| h.powf(1.5));
            (
                Variogram {
                    kind,
                    nugget,
                    scale,
                    range: max_h,
                },
                rmse,
            )
        }
        Kind::Exponential => {
            // Golden-section search over the (nonlinear) range; nugget/scale
            // are a closed-form inner WLS fit for each candidate range.
            let lo0 = (bin_w * 0.25).max(1e-9);
            let hi0 = (max_h * 3.0).max(lo0 * 2.0);
            let phi = 0.618_033_988_749_895;
            let (mut lo, mut hi) = (lo0, hi0);
            let eval = |range: f64| -> (f64, f64, f64, f64) {
                let (nugget, scale, rmse) = wls_fit_linear(&bins, |h| 1.0 - (-h / range).exp());
                let sse: f64 = bins
                    .iter()
                    .map(|&(h, g, w)| {
                        let pred = nugget + scale * (1.0 - (-h / range).exp());
                        w * (g - pred).powi(2)
                    })
                    .sum();
                (sse, nugget, scale, rmse)
            };
            let mut x1 = hi - phi * (hi - lo);
            let mut x2 = lo + phi * (hi - lo);
            let mut f1 = eval(x1).0;
            let mut f2 = eval(x2).0;
            for _ in 0..40 {
                if hi - lo < 1e-6 {
                    break;
                }
                if f1 < f2 {
                    hi = x2;
                    x2 = x1;
                    f2 = f1;
                    x1 = hi - phi * (hi - lo);
                    f1 = eval(x1).0;
                } else {
                    lo = x1;
                    x1 = x2;
                    f1 = f2;
                    x2 = lo + phi * (hi - lo);
                    f2 = eval(x2).0;
                }
            }
            let best_range = (lo + hi) / 2.0;
            let (_, nugget, scale, rmse) = eval(best_range);
            (
                Variogram {
                    kind,
                    nugget,
                    scale,
                    range: best_range,
                },
                rmse,
            )
        }
    }
}

/// Draws a perturbed semivariogram: log-normal jitter on nugget/scale (and
/// range, for `exponential`) scaled by `sigma` (the base fit's relative WLS
/// residual, clamped to a sane band). This is EBK's "simulate a new
/// semivariogram from the fitted model" step, simplified to a seeded parameter
/// perturbation rather than a full data-simulate-and-refit bootstrap — see the
/// module-level scope-cut notes.
fn simulate_variogram(base: &Variogram, sigma: f64, rng: &mut Rng) -> Variogram {
    // Bounded so a single unlucky simulation cannot dominate its subset's
    // mean (Box-Muller tails are otherwise unbounded); a factor of e^1 ~= 2.7x
    // is already a generous swing for a semivariogram parameter.
    let sigma = sigma.clamp(0.05, 0.3);
    let jitter = |rng: &mut Rng, s: f64| (rng.normal().clamp(-3.0, 3.0) * s).exp();
    let eps = 1e-6 * base.scale.max(1.0);
    Variogram {
        kind: base.kind,
        nugget: (base.nugget + eps) * jitter(rng, sigma),
        scale: (base.scale + eps) * jitter(rng, sigma),
        range: (base.range + eps) * jitter(rng, sigma * 0.5),
    }
}

// ── Seeded deterministic RNG (splitmix64 + Box-Muller) ──────────────────────

struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Rng { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform sample in (0, 1), never exactly 0 (safe for `ln`).
    fn next_f64(&mut self) -> f64 {
        let bits = self.next_u64() >> 11; // 53 significant bits
        ((bits as f64) + 1.0) / ((1u64 << 53) as f64 + 1.0)
    }

    fn normal(&mut self) -> f64 {
        let u1 = self.next_f64();
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }
}

/// Derives a deterministic per-(subset, simulation) seed from the global seed
/// so results are reproducible but not correlated across subsets/simulations.
fn derive_seed(seed: u64, subset_idx: usize, sim_idx: usize) -> u64 {
    let a = seed ^ (subset_idx as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    a ^ (sim_idx as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9)
}

// ── Linear algebra: hand-rolled Cholesky (SPD) and Gauss-Jordan inverse ─────

/// Solves the `n x n` symmetric positive-definite system `a * x = b` (row-major
/// `a`) via Cholesky decomposition. `None` if `a` is not (numerically) SPD.
fn cholesky_solve(a: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut l = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i * n + j];
            for k in 0..j {
                sum -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if sum <= 0.0 {
                    return None;
                }
                l[i * n + j] = sum.sqrt();
            } else {
                l[i * n + j] = sum / l[j * n + j];
            }
        }
    }
    // Forward solve L*y = b.
    let mut y = vec![0.0f64; n];
    for i in 0..n {
        let mut sum = b[i];
        for k in 0..i {
            sum -= l[i * n + k] * y[k];
        }
        y[i] = sum / l[i * n + i];
    }
    // Backward solve L^T*x = y.
    let mut x = vec![0.0f64; n];
    for ii in 0..n {
        let i = n - 1 - ii;
        let mut sum = y[i];
        for k in (i + 1)..n {
            sum -= l[k * n + i] * x[k];
        }
        x[i] = sum / l[i * n + i];
    }
    Some(x)
}

/// Inverts the `n x n` matrix `a` (row-major) via Gauss-Jordan elimination with
/// partial pivoting. `None` if `a` is (numerically) singular. Used for the
/// symmetric-but-indefinite ordinary-kriging system (the Lagrange-multiplier
/// row/column makes it non-SPD, so Cholesky does not apply there).
fn gauss_jordan_inverse(a: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut m = vec![0.0f64; n * 2 * n];
    for i in 0..n {
        for j in 0..n {
            m[i * 2 * n + j] = a[i * n + j];
        }
        m[i * 2 * n + n + i] = 1.0;
    }
    for col in 0..n {
        let mut pivot = col;
        let mut best = m[col * 2 * n + col].abs();
        for r in (col + 1)..n {
            let v = m[r * 2 * n + col].abs();
            if v > best {
                best = v;
                pivot = r;
            }
        }
        if best < 1e-12 {
            return None;
        }
        if pivot != col {
            for k in 0..2 * n {
                m.swap(col * 2 * n + k, pivot * 2 * n + k);
            }
        }
        let piv = m[col * 2 * n + col];
        for k in 0..2 * n {
            m[col * 2 * n + k] /= piv;
        }
        for r in 0..n {
            if r == col {
                continue;
            }
            let factor = m[r * 2 * n + col];
            if factor == 0.0 {
                continue;
            }
            for k in 0..2 * n {
                m[r * 2 * n + k] -= factor * m[col * 2 * n + k];
            }
        }
    }
    let mut inv = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            inv[i * n + j] = m[i * 2 * n + n + j];
        }
    }
    Some(inv)
}

// ── Ordinary kriging over one local subset + one simulated semivariogram ────

/// Ordinary-kriging weight system for `pts` under `variogram`, pre-factorized
/// once as a matrix inverse `(n+1) x (n+1)` and reused for every query cell
/// that subset/simulation touches (only the right-hand side changes per
/// query).
struct KrigeSystem {
    inv: Vec<f64>,
    n: usize,
}

fn build_krige_system(pts: &[Pt], variogram: &Variogram) -> Option<KrigeSystem> {
    let n = pts.len();
    let m = n + 1;
    let mut a = vec![0.0f64; m * m];
    for i in 0..n {
        for j in 0..n {
            a[i * m + j] = gamma_h(variogram, dist(&pts[i], pts[j].x, pts[j].y));
        }
        a[i * m + n] = 1.0;
        a[n * m + i] = 1.0;
    }
    a[n * m + n] = 0.0;
    let inv = gauss_jordan_inverse(&a, m)?;
    Some(KrigeSystem { inv, n })
}

/// Predicts at `(qx, qy)` using a pre-built kriging system. Coincident with an
/// input point (within a tiny epsilon) short-circuits to that point's exact
/// value and zero variance — kriging's exact-interpolator property, made
/// robust to any nugget introduced by semivariogram simulation.
fn krige_predict(
    pts: &[Pt],
    sys: &KrigeSystem,
    variogram: &Variogram,
    qx: f64,
    qy: f64,
) -> (f64, f64) {
    for p in pts {
        if dist(p, qx, qy) < 1e-9 {
            return (p.v, 0.0);
        }
    }
    let m = sys.n + 1;
    let mut g = vec![0.0f64; m];
    for (i, p) in pts.iter().enumerate() {
        g[i] = gamma_h(variogram, dist(p, qx, qy));
    }
    g[sys.n] = 1.0;
    let mut w = vec![0.0f64; m];
    for (k, wk) in w.iter_mut().enumerate() {
        let mut s = 0.0;
        for (j, &gj) in g.iter().enumerate() {
            s += sys.inv[k * m + j] * gj;
        }
        *wk = s;
    }
    let mut pred = 0.0;
    for (i, p) in pts.iter().enumerate() {
        pred += w[i] * p.v;
    }
    let var: f64 = w.iter().zip(&g).map(|(wk, gk)| wk * gk).sum();
    (pred, var.max(0.0))
}

// ── Prepared (fit + simulated + factorized) local subsets ──────────────────

struct PreparedSubset {
    pts: Vec<Pt>,
    cx: f64,
    cy: f64,
    sims: Vec<(Variogram, KrigeSystem)>,
}

fn prepare_subset(
    points: &[Pt],
    span: &SubsetSpan,
    kind: Kind,
    simulations: usize,
    seed: u64,
    subset_idx: usize,
) -> PreparedSubset {
    let pts: Vec<Pt> = span.idxs.iter().map(|&i| points[i]).collect();
    let (base, rel_rmse) = fit_variogram(&pts, kind);
    let mut sims = Vec::with_capacity(simulations);
    for k in 0..simulations {
        let mut rng = Rng::new(derive_seed(seed, subset_idx, k));
        let mut vg = simulate_variogram(&base, rel_rmse, &mut rng);
        let sys = loop {
            if let Some(sys) = build_krige_system(&pts, &vg) {
                break sys;
            }
            // Extremely rare (near-singular Gram matrix from an unlucky
            // draw): nudge the nugget up as a ridge and retry.
            vg.nugget = vg.nugget.max(1e-6) * 4.0;
            if vg.nugget > 1e6 {
                break build_krige_system(&pts, &vg)
                    .unwrap_or_else(|| KrigeSystem { inv: vec![], n: 0 });
            }
        };
        if sys.n > 0 {
            sims.push((vg, sys));
        }
    }
    PreparedSubset {
        pts,
        cx: span.cx,
        cy: span.cy,
        sims,
    }
}

/// Blends the `m` nearest prepared subsets (by centroid distance) at query
/// point `(qx, qy)`. Each subset's prediction/variance is the mean over its
/// simulated semivariograms plus the variance *of* those simulated
/// predictions (law of total variance — the model-uncertainty term EBK adds
/// on top of plain kriging variance); subsets are then combined by
/// inverse-squared-distance weight, adding the spread *between* subset
/// predictions as a further, nonstationarity-driven variance term.
fn blend_cell(prepared: &[PreparedSubset], qx: f64, qy: f64, m: usize) -> (f64, f64) {
    let mut order: Vec<usize> = (0..prepared.len()).collect();
    order.sort_by(|&a, &b| {
        let da = (prepared[a].cx - qx).powi(2) + (prepared[a].cy - qy).powi(2);
        let db = (prepared[b].cx - qx).powi(2) + (prepared[b].cy - qy).powi(2);
        da.partial_cmp(&db).unwrap()
    });
    order.truncate(m.max(1));

    let mut weights = Vec::with_capacity(order.len());
    let mut preds = Vec::with_capacity(order.len());
    let mut vars = Vec::with_capacity(order.len());
    for &si in &order {
        let s = &prepared[si];
        if s.sims.is_empty() {
            continue;
        }
        let mut p_sum = 0.0;
        let mut v_sum = 0.0;
        let mut p_sq_sum = 0.0;
        for (vg, sys) in &s.sims {
            let (p, v) = krige_predict(&s.pts, sys, vg, qx, qy);
            p_sum += p;
            v_sum += v;
            p_sq_sum += p * p;
        }
        let k = s.sims.len() as f64;
        let mean_p = p_sum / k;
        let mean_v = v_sum / k;
        let var_of_p = (p_sq_sum / k - mean_p * mean_p).max(0.0);
        let d2 = (s.cx - qx).powi(2) + (s.cy - qy).powi(2);
        let w = 1.0 / (d2 + 1e-6);
        weights.push(w);
        preds.push(mean_p);
        vars.push(mean_v + var_of_p);
    }
    if weights.is_empty() {
        return (f64::NAN, f64::NAN);
    }
    let w_sum: f64 = weights.iter().sum();
    let final_pred: f64 = preds.iter().zip(&weights).map(|(p, w)| p * w).sum::<f64>() / w_sum;
    let mean_var: f64 = vars.iter().zip(&weights).map(|(v, w)| v * w).sum::<f64>() / w_sum;
    let between: f64 = preds
        .iter()
        .zip(&weights)
        .map(|(p, w)| w * (p - final_pred).powi(2))
        .sum::<f64>()
        / w_sum;
    (final_pred, (mean_var + between).max(0.0).sqrt())
}

// ── Raster output ────────────────────────────────────────────────────────────

#[allow(clippy::too_many_arguments)]
fn build_raster(
    z: &[f64],
    rows: usize,
    cols: usize,
    x_min: f64,
    y_min: f64,
    cell_size: f64,
    nodata: f64,
    epsg: Option<u32>,
) -> Result<Raster, ToolError> {
    let crs = match epsg {
        Some(code) => CrsInfo::from_epsg(code),
        None => CrsInfo::new(),
    };
    let mut raster = Raster::new(RasterConfig {
        cols,
        rows,
        bands: 1,
        x_min,
        y_min,
        cell_size,
        cell_size_y: Some(cell_size),
        nodata,
        data_type: DataType::F32,
        crs,
        metadata: Vec::new(),
    });
    for r in 0..rows {
        for c in 0..cols {
            let v = z[r * cols + c];
            let v = if v.is_nan() { nodata } else { v };
            raster
                .set(0, r as isize, c as isize, v)
                .map_err(|e| ToolError::Execution(format!("failed writing cell: {e}")))?;
        }
    }
    Ok(raster)
}

// ── Parameters ───────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Transform {
    None,
    LogEmpirical,
}

struct Params {
    cell_size: Option<f64>,
    subset_size: usize,
    overlap: f64,
    simulations: usize,
    kind: Kind,
    transform: Transform,
    seed: u64,
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let cell_size = match parse_f64(args, "cell_size")? {
        None => None,
        Some(cs) if cs > 0.0 && cs.is_finite() => Some(cs),
        Some(_) => {
            return Err(ToolError::Validation(
                "'cell_size' must be a positive number".to_string(),
            ))
        }
    };

    let subset_size = match parse_u64(args, "subset_size")? {
        None => 40usize,
        Some(v) if v >= 4 => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "'subset_size' must be >= 4".to_string(),
            ))
        }
    };

    let overlap = match parse_f64(args, "overlap")? {
        None => 0.25,
        Some(v) if v >= 0.0 && v.is_finite() => v,
        Some(_) => return Err(ToolError::Validation("'overlap' must be >= 0".to_string())),
    };

    let simulations = match parse_u64(args, "simulations")? {
        None => 5usize,
        Some(v) if (1..=50).contains(&v) => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "'simulations' must be between 1 and 50".to_string(),
            ))
        }
    };

    let kind = match args
        .get("semivariogram")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
    {
        None => Kind::Power,
        Some(ref s) if s.is_empty() || s == "power" => Kind::Power,
        Some(ref s) if s == "linear" => Kind::Linear,
        Some(ref s) if s == "exponential" => Kind::Exponential,
        Some(s) => {
            return Err(ToolError::Validation(format!(
                "'semivariogram' must be power/linear/exponential, got '{s}'"
            )))
        }
    };

    let transform = match args
        .get("transform")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
    {
        None => Transform::None,
        Some(ref s) if s.is_empty() || s == "none" => Transform::None,
        Some(ref s) if s == "log_empirical" => Transform::LogEmpirical,
        Some(s) => {
            return Err(ToolError::Validation(format!(
                "'transform' must be none/log_empirical, got '{s}'"
            )))
        }
    };

    let seed = parse_u64(args, "seed")?.unwrap_or(42);

    Ok(Params {
        cell_size,
        subset_size,
        overlap,
        simulations,
        kind,
        transform,
        seed,
    })
}

fn parse_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
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

fn parse_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64().or_else(|| n.as_f64().map(|f| f as u64))),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be an integer"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be an integer"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn point_layer(pts: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (x, y, v) in pts {
            l.add_feature(Some(Geometry::point(*x, *y)), &[("v", (*v).into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = EmpiricalBayesianKrigingTool.run(&args, &ctx()).unwrap();
        let r = crate::common::load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, r)
    }

    /// Random-ish but deterministic scattered point set covering a modest
    /// extent, large enough to form several local subsets.
    fn scattered(n: usize, seed: u64) -> Vec<(f64, f64, f64)> {
        let mut rng = Rng::new(seed);
        (0..n)
            .map(|_| {
                let x = rng.next_f64() * 100.0;
                let y = rng.next_f64() * 100.0;
                let v = 10.0 + 0.3 * x - 0.2 * y + 5.0 * (x / 20.0).sin();
                (x, y, v)
            })
            .collect()
    }

    // ── Core: exact interpolation ────────────────────────────────────────────

    /// Kriging is an exact interpolator (zero nugget at the query location):
    /// predicting AT an input point must return ~that point's value with ~0
    /// variance, regardless of the simulated semivariogram's own nugget.
    #[test]
    fn exact_interpolation_at_data_points() {
        let pts: Vec<Pt> = scattered(20, 7)
            .into_iter()
            .map(|(x, y, v)| Pt { x, y, v })
            .collect();
        let (variogram, _) = fit_variogram(&pts, Kind::Power);
        let sys = build_krige_system(&pts, &variogram).unwrap();
        for p in &pts {
            let (pred, var) = krige_predict(&pts, &sys, &variogram, p.x, p.y);
            assert!(
                (pred - p.v).abs() < 1e-6,
                "prediction at a data point should equal its value: got {pred}, want {}",
                p.v
            );
            assert!(
                var.abs() < 1e-9,
                "variance at a data point should be ~0, got {var}"
            );
        }
    }

    /// The same exact-interpolation property holds (approximately) through the
    /// full tool for a spatially-*coherent* field: the raster cell nearest
    /// each input point predicts close to its value. (A field with no spatial
    /// structure at all is legitimately smoothed at even a fraction-of-a-cell
    /// offset by a fitted near-pure-nugget semivariogram — that is correct
    /// kriging behaviour, not a bug, so this test uses a smooth trend+wave
    /// field like real elevation/temperature data instead.)
    #[test]
    fn full_tool_is_near_exact_at_input_points() {
        let pts = scattered(24, 9);
        let path = point_layer(&pts);
        let (_out, r) = run(json!({
            "input": path, "field": "v", "cell_size": 1.0,
            "subset_size": 12, "simulations": 3,
        }));
        for (x, y, v) in pts {
            let col = ((x - r.x_min) / r.cell_size_x).floor() as isize;
            let row = ((r.y_max() - y) / r.cell_size_y).floor() as isize;
            let cell = r.get(0, row, col);
            assert!(
                (cell - v).abs() < 8.0,
                "cell near ({x},{y}) should be close to {v}, got {cell}"
            );
        }
    }

    // ── Core: linear trend reproduction ─────────────────────────────────────

    /// A pure linear trend field (no noise), densely sampled, is reproduced by
    /// ordinary kriging within a modest tolerance at an interior query point.
    #[test]
    fn reproduces_linear_trend() {
        let mut pts = Vec::new();
        for i in 0..7 {
            for j in 0..7 {
                let x = i as f64 * 10.0;
                let y = j as f64 * 10.0;
                pts.push(Pt {
                    x,
                    y,
                    v: 5.0 + 2.0 * x + 1.5 * y,
                });
            }
        }
        let (variogram, _) = fit_variogram(&pts, Kind::Linear);
        let sys = build_krige_system(&pts, &variogram).unwrap();
        let (qx, qy) = (33.0, 27.0);
        let truth = 5.0 + 2.0 * qx + 1.5 * qy;
        let (pred, _var) = krige_predict(&pts, &sys, &variogram, qx, qy);
        assert!(
            (pred - truth).abs() < truth.abs() * 0.05 + 2.0,
            "interior prediction {pred} should track the linear trend {truth}"
        );
    }

    // ── Core: hand-rolled linear solver ──────────────────────────────────────

    /// `cholesky_solve` reproduces the known solution of a textbook SPD
    /// system.
    #[test]
    fn cholesky_solves_known_spd_system() {
        // A = [[4,1,0],[1,3,1],[0,1,2]], symmetric positive-definite (diagonally
        // dominant). Pick x = [1,2,3] and derive b = A*x exactly.
        let a = [4.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 2.0];
        let x_true = [1.0, 2.0, 3.0];
        let mut b = [0.0; 3];
        for i in 0..3 {
            for j in 0..3 {
                b[i] += a[i * 3 + j] * x_true[j];
            }
        }
        let x = cholesky_solve(&a, &b, 3).expect("SPD system must solve");
        for i in 0..3 {
            assert!(
                (x[i] - x_true[i]).abs() < 1e-9,
                "component {i}: got {}, want {}",
                x[i],
                x_true[i]
            );
        }
    }

    #[test]
    fn cholesky_rejects_non_spd() {
        // Not positive-definite (leading principal minor <= 0 after the first
        // pivot's Schur complement).
        let a = [1.0, 2.0, 2.0, 1.0];
        assert!(cholesky_solve(&a, &[1.0, 1.0], 2).is_none());
    }

    #[test]
    fn gauss_jordan_inverts_identity_like_system() {
        let a = [2.0, 0.0, 0.0, 3.0];
        let inv = gauss_jordan_inverse(&a, 2).unwrap();
        assert!((inv[0] - 0.5).abs() < 1e-9);
        assert!((inv[3] - 1.0 / 3.0).abs() < 1e-9);
        assert!(inv[1].abs() < 1e-9 && inv[2].abs() < 1e-9);
    }

    // ── Determinism ───────────────────────────────────────────────────────────

    /// Same seed -> bit-for-bit identical output raster.
    #[test]
    fn deterministic_with_fixed_seed() {
        let pts = scattered(60, 11);
        let path_a = point_layer(&pts);
        let path_b = point_layer(&pts);
        let (_oa, ra) = run(json!({
            "input": path_a, "field": "v", "cell_size": 3.0, "seed": 123, "simulations": 4,
        }));
        let (_ob, rb) = run(json!({
            "input": path_b, "field": "v", "cell_size": 3.0, "seed": 123, "simulations": 4,
        }));
        assert_eq!(ra.rows, rb.rows);
        assert_eq!(ra.cols, rb.cols);
        for row in 0..ra.rows as isize {
            for col in 0..ra.cols as isize {
                let va = ra.get(0, row, col);
                let vb = rb.get(0, row, col);
                assert_eq!(
                    va.to_bits(),
                    vb.to_bits(),
                    "cell ({row},{col}) differs between identical-seed runs: {va} vs {vb}"
                );
            }
        }
    }

    /// Different seeds produce a measurably different raster (the simulation
    /// step is actually doing something) while staying in the same ballpark.
    #[test]
    fn different_seeds_change_output_but_stay_close() {
        let pts = scattered(60, 12);
        let path_a = point_layer(&pts);
        let path_b = point_layer(&pts);
        let (_oa, ra) = run(json!({
            "input": path_a, "field": "v", "cell_size": 5.0, "seed": 1, "simulations": 4,
        }));
        let (_ob, rb) = run(json!({
            "input": path_b, "field": "v", "cell_size": 5.0, "seed": 2, "simulations": 4,
        }));
        let mut any_diff = false;
        let mut max_diff = 0.0f64;
        for row in 0..ra.rows as isize {
            for col in 0..ra.cols as isize {
                let va = ra.get(0, row, col);
                let vb = rb.get(0, row, col);
                if !va.is_nan() && !vb.is_nan() {
                    if va.to_bits() != vb.to_bits() {
                        any_diff = true;
                    }
                    max_diff = max_diff.max((va - vb).abs());
                }
            }
        }
        assert!(any_diff, "different seeds should perturb at least one cell");
        assert!(
            max_diff < 20.0,
            "seed-driven spread should stay bounded, got max diff {max_diff}"
        );
    }

    // ── Transform ────────────────────────────────────────────────────────────

    #[test]
    fn log_empirical_transform_stays_positive() {
        let pts: Vec<(f64, f64, f64)> = scattered(40, 21)
            .into_iter()
            .map(|(x, y, v)| (x, y, (v + 20.0).max(1.0)))
            .collect();
        let path = point_layer(&pts);
        let (_out, r) = run(json!({
            "input": path, "field": "v", "cell_size": 5.0, "transform": "log_empirical",
        }));
        for row in 0..r.rows as isize {
            for col in 0..r.cols as isize {
                let v = r.get(0, row, col);
                if !v.is_nan() {
                    assert!(
                        v > 0.0,
                        "log_empirical predictions must stay positive, got {v}"
                    );
                }
            }
        }
    }

    #[test]
    fn log_empirical_rejects_nonpositive_values() {
        let pts = [
            (0.0, 0.0, 1.0),
            (10.0, 0.0, -1.0),
            (0.0, 10.0, 2.0),
            (10.0, 10.0, 3.0),
            (5.0, 5.0, 4.0),
            (2.0, 8.0, 5.0),
            (8.0, 2.0, 6.0),
            (3.0, 3.0, 7.0),
        ];
        let path = point_layer(&pts);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": path, "field": "v", "transform": "log_empirical",
        }))
        .unwrap();
        assert!(EmpiricalBayesianKrigingTool.run(&args, &ctx()).is_err());
    }

    // ── error_output ─────────────────────────────────────────────────────────

    #[test]
    fn error_output_is_nonnegative() {
        let pts = scattered(50, 5);
        let path = point_layer(&pts);
        let se_path = std::env::temp_dir()
            .join("geolibre_ebk_se_test.tif")
            .to_string_lossy()
            .into_owned();
        let args: ToolArgs = serde_json::from_value(json!({
            "input": path, "field": "v", "cell_size": 6.0, "error_output": se_path,
        }))
        .unwrap();
        let out = EmpiricalBayesianKrigingTool.run(&args, &ctx()).unwrap();
        let se = crate::common::load_input_raster(out.outputs["error_output"].as_str().unwrap())
            .unwrap();
        for row in 0..se.rows as isize {
            for col in 0..se.cols as isize {
                let v = se.get(0, row, col);
                if !v.is_nan() {
                    assert!(v >= 0.0, "standard error must be >= 0, got {v}");
                }
            }
        }
    }

    // ── Parameter validation ─────────────────────────────────────────────────

    #[test]
    fn rejects_bad_parameters() {
        let path = point_layer(&scattered(20, 3));
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            EmpiricalBayesianKrigingTool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input/field");
        assert!(bad(json!({ "input": path })).is_err(), "missing field");
        assert!(
            bad(json!({ "input": path, "field": "v", "semivariogram": "gaussian" })).is_err(),
            "bad semivariogram"
        );
        assert!(
            bad(json!({ "input": path, "field": "v", "transform": "boxcox" })).is_err(),
            "bad transform"
        );
        assert!(
            bad(json!({ "input": path, "field": "v", "subset_size": 2 })).is_err(),
            "subset_size too small"
        );
        assert!(
            bad(json!({ "input": path, "field": "v", "overlap": -0.5 })).is_err(),
            "negative overlap"
        );
        assert!(
            bad(json!({ "input": path, "field": "v", "simulations": 0 })).is_err(),
            "zero simulations"
        );
        assert!(
            bad(json!({ "input": path, "field": "v", "cell_size": -1.0 })).is_err(),
            "negative cell_size"
        );
        assert!(
            bad(json!({ "input": path, "field": "v", "semivariogram": "exponential", "transform": "log_empirical" }))
                .is_ok(),
            "valid combination"
        );
    }

    #[test]
    fn rejects_too_few_points() {
        let path = point_layer(&[(0.0, 0.0, 1.0), (1.0, 1.0, 2.0)]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": path, "field": "v" })).unwrap();
        assert!(EmpiricalBayesianKrigingTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn rejects_missing_field() {
        let path = point_layer(&scattered(10, 4));
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": path, "field": "not_a_field" })).unwrap();
        assert!(EmpiricalBayesianKrigingTool.run(&args, &ctx()).is_err());
    }
}
