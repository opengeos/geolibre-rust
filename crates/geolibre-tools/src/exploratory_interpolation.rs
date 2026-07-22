//! GeoLibre tool: exploratory interpolation (cross-validate & rank methods).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Exploratory Interpolation*
//! (Geostatistical Analyst). The bundled suite ships several interpolators
//! (`idw_interpolation`, the kriging/variogram family, `thin_plate_spline`,
//! `natural_neighbour_interpolation`, `tin_interpolation`, plus GeoLibre's own
//! `interpolate_with_barriers` and `empirical_bayesian_kriging`), but nothing
//! that runs *several* methods, leave-one-out cross-validates each, and ranks
//! them. Users have to hand-run each interpolator and compare by eye.
//!
//! This tool orchestrates a family of deterministic, WASM-safe interpolators:
//!
//! * `idw` — inverse-distance weighting, `w = 1 / d^power`;
//! * `nearest` — nearest-neighbour (Thiessen) assignment;
//! * `trend1` — global first-order polynomial trend surface (planar least
//!   squares); and
//! * `trend2` — global second-order polynomial trend surface.
//!
//! For every requested method it runs **leave-one-out cross-validation
//! (LOOCV)**: each sample is held out in turn and predicted from the remaining
//! samples, yielding, per method, the classic geostatistical validation
//! metrics — mean error `ME` (prediction bias), mean absolute error `MAE`,
//! root-mean-square error `RMSE`, the error range, and the Pearson correlation
//! between predicted and observed values. Methods are then ranked by the
//! chosen `criterion` (`rmse` ascending, `mae` ascending, or `me` by smallest
//! absolute bias). The winner is reported, and when `output_raster` is given
//! the tool fits that method on *all* samples and writes its prediction grid.
//!
//! `RMSSE` (root-mean-square standardized error) is intentionally out of scope
//! for v1: it needs a per-prediction kriging standard error, which these
//! deterministic interpolators do not produce.
//!
//! The primary output is a ranked comparison table — CSV when `output` ends in
//! `.csv`, otherwise a geometry-less vector table — one row per method.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
use wbvector::{Feature, FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::common::{write_or_store_output, write_text_output};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

const OUT_NODATA: f64 = -9999.0;
/// Hard cap on winner-raster grid dimensions to keep a single run tractable.
const MAX_DIM: usize = 4000;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Idw,
    Nearest,
    Trend1,
    Trend2,
}

impl Method {
    fn id(self) -> &'static str {
        match self {
            Method::Idw => "idw",
            Method::Nearest => "nearest",
            Method::Trend1 => "trend1",
            Method::Trend2 => "trend2",
        }
    }

    fn parse(s: &str) -> Option<Method> {
        match s.trim().to_ascii_lowercase().as_str() {
            "idw" | "inverse_distance" => Some(Method::Idw),
            "nearest" | "nearest_neighbour" | "nearest_neighbor" | "thiessen" => {
                Some(Method::Nearest)
            }
            "trend" | "trend1" | "polynomial1" | "linear" => Some(Method::Trend1),
            "trend2" | "polynomial2" | "quadratic" => Some(Method::Trend2),
            _ => None,
        }
    }

    /// Number of polynomial terms for the trend methods (else 0).
    fn trend_terms(self) -> usize {
        match self {
            Method::Trend1 => 3, // 1, x, y
            Method::Trend2 => 6, // 1, x, y, x^2, xy, y^2
            _ => 0,
        }
    }
}

const ALL_METHODS: [Method; 4] = [Method::Idw, Method::Nearest, Method::Trend1, Method::Trend2];

#[derive(Clone, Copy, PartialEq, Eq)]
enum Criterion {
    Rmse,
    Mae,
    Me,
}

impl Criterion {
    fn id(self) -> &'static str {
        match self {
            Criterion::Rmse => "rmse",
            Criterion::Mae => "mae",
            Criterion::Me => "me",
        }
    }
}

pub struct ExploratoryInterpolationTool;

impl Tool for ExploratoryInterpolationTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "exploratory_interpolation",
            display_name: "Exploratory Interpolation",
            summary: "Run several interpolation methods (IDW, nearest-neighbour, first/second-order trend surfaces), leave-one-out cross-validate each, and rank them by prediction accuracy (RMSE/MAE) or bias (ME) — like ArcGIS Exploratory Interpolation. Reports a ranked per-method comparison table and can emit the best method's prediction raster. Subsumes a standalone cross-validation request.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer of measurements.",
                    required: true,
                },
                ToolParamSpec {
                    name: "field",
                    description: "Numeric field on the points holding the value to interpolate.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output ranked comparison table — a CSV (extension .csv) or a geometry-less vector table. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "methods",
                    description: "Comma/pipe-separated subset of methods to compare: 'idw', 'nearest', 'trend1', 'trend2'. Default: all four.",
                    required: false,
                },
                ToolParamSpec {
                    name: "criterion",
                    description: "Ranking criterion: 'rmse' (default, lowest error), 'mae' (lowest absolute error), or 'me' (least bias, smallest |mean error|).",
                    required: false,
                },
                ToolParamSpec {
                    name: "power",
                    description: "IDW distance-decay exponent (default 2).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_raster",
                    description: "Optional output GeoTIFF: the best-ranked method fit on all samples, evaluated over a grid spanning the sample extent.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Winner-raster cell size in CRS units (degrees for a geographic CRS). Default: sized so the longer axis spans ~256 cells.",
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
        let field = require_str(args, "field")?.to_string();
        let output = parse_optional_str(args, "output")?;
        let output_raster = crate::common::parse_optional_output(args, "output_raster")?;
        let prm = parse_params(args)?;

        // ── Load samples ─────────────────────────────────────────────────────
        let layer = load_input_layer(input)?;
        let field_idx = layer
            .schema
            .field_index(&field)
            .ok_or_else(|| ToolError::Validation(format!("field '{field}' not found")))?;

        let mut sx: Vec<f64> = Vec::new();
        let mut sy: Vec<f64> = Vec::new();
        let mut sv: Vec<f64> = Vec::new();
        for feat in layer.iter() {
            let Some(geom) = feat.geometry.as_ref() else {
                continue;
            };
            let Some((x, y)) = point_xy(geom) else {
                continue;
            };
            let Some(v) = feat.attributes.get(field_idx).and_then(|f| f.as_f64()) else {
                continue;
            };
            if !v.is_finite() {
                continue;
            }
            sx.push(x);
            sy.push(y);
            sv.push(v);
        }
        let n = sx.len();
        if n < 3 {
            return Err(ToolError::Execution(format!(
                "need at least 3 point features with a finite '{field}' value, got {n}"
            )));
        }

        // Global coordinate normalization keeps the trend design matrix
        // well-conditioned; the same transform is reused for every fold.
        let (cx, cy, scale) = normalize_params(&sx, &sy);

        ctx.progress.info(&format!(
            "{n} samples, {} method(s), LOOCV, rank by {}",
            prm.methods.len(),
            prm.criterion.id()
        ));

        // ── Leave-one-out cross-validate every requested method ──────────────
        let mut results: Vec<MethodResult> = Vec::with_capacity(prm.methods.len());
        for (mi, &m) in prm.methods.iter().enumerate() {
            let mut se = 0.0_f64; // sum of squared error
            let mut sae = 0.0_f64; // sum of absolute error
            let mut serr = 0.0_f64; // sum of signed error
            let mut emin = f64::INFINITY;
            let mut emax = f64::NEG_INFINITY;
            // For the predicted-vs-observed Pearson correlation.
            let (mut sp, mut so, mut spp, mut soo, mut spo) = (0.0, 0.0, 0.0, 0.0, 0.0);
            let mut valid = 0usize;

            for i in 0..n {
                let Some(pred) = predict(
                    m, i, &sx, &sy, &sv, prm.power, cx, cy, scale, /*loo=*/ true,
                ) else {
                    continue;
                };
                let obs = sv[i];
                let err = pred - obs;
                se += err * err;
                sae += err.abs();
                serr += err;
                emin = emin.min(err);
                emax = emax.max(err);
                sp += pred;
                so += obs;
                spp += pred * pred;
                soo += obs * obs;
                spo += pred * obs;
                valid += 1;
            }
            if valid == 0 {
                return Err(ToolError::Execution(format!(
                    "method '{}' produced no cross-validation predictions (too few samples?)",
                    m.id()
                )));
            }
            let nf = valid as f64;
            let rmse = (se / nf).sqrt();
            let mae = sae / nf;
            let me = serr / nf;
            // Pearson r between predictions and observations.
            let cov = spo / nf - (sp / nf) * (so / nf);
            let vp = (spp / nf - (sp / nf).powi(2)).max(0.0);
            let vo = (soo / nf - (so / nf).powi(2)).max(0.0);
            let pearson = if vp > 0.0 && vo > 0.0 {
                cov / (vp.sqrt() * vo.sqrt())
            } else {
                f64::NAN
            };
            results.push(MethodResult {
                method: m,
                n: valid,
                me,
                mae,
                rmse,
                err_min: emin,
                err_max: emax,
                pearson,
            });
            ctx.progress
                .progress((mi as f64 + 1.0) / prm.methods.len() as f64);
        }

        // ── Rank by the chosen criterion (ascending goodness) ────────────────
        let key = |r: &MethodResult| -> f64 {
            match prm.criterion {
                Criterion::Rmse => r.rmse,
                Criterion::Mae => r.mae,
                Criterion::Me => r.me.abs(),
            }
        };
        // Stable order: primary criterion, then RMSE, then method id for
        // determinism.
        let mut order: Vec<usize> = (0..results.len()).collect();
        order.sort_by(|&a, &b| {
            key(&results[a])
                .partial_cmp(&key(&results[b]))
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(
                    results[a]
                        .rmse
                        .partial_cmp(&results[b].rmse)
                        .unwrap_or(std::cmp::Ordering::Equal),
                )
                .then(results[a].method.id().cmp(results[b].method.id()))
        });
        let mut rank_of = vec![0usize; results.len()];
        for (rank, &idx) in order.iter().enumerate() {
            rank_of[idx] = rank + 1;
        }
        let best = order[0];
        let best_method = results[best].method;

        // ── Comparison table (one row per method, in rank order) ─────────────
        let mut table = Layer::new("exploratory_interpolation");
        for (name, ty) in [
            ("rank", FieldType::Integer),
            ("method", FieldType::Text),
            ("n", FieldType::Integer),
            ("me", FieldType::Float),
            ("mae", FieldType::Float),
            ("rmse", FieldType::Float),
            ("err_min", FieldType::Float),
            ("err_max", FieldType::Float),
            ("pearson_r", FieldType::Float),
            ("best", FieldType::Boolean),
        ] {
            table.add_field(FieldDef::new(name, ty));
        }
        let mut csv = String::from("rank,method,n,me,mae,rmse,err_min,err_max,pearson_r,best\n");
        for (fid, &idx) in order.iter().enumerate() {
            let r = &results[idx];
            let is_best = idx == best;
            table.push(Feature {
                fid: fid as u64,
                geometry: None,
                attributes: vec![
                    FieldValue::Integer(rank_of[idx] as i64),
                    FieldValue::Text(r.method.id().to_string()),
                    FieldValue::Integer(r.n as i64),
                    FieldValue::Float(r.me),
                    FieldValue::Float(r.mae),
                    FieldValue::Float(r.rmse),
                    FieldValue::Float(r.err_min),
                    FieldValue::Float(r.err_max),
                    FieldValue::Float(r.pearson),
                    FieldValue::Boolean(is_best),
                ],
            });
            csv.push_str(&format!(
                "{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{}\n",
                rank_of[idx],
                r.method.id(),
                r.n,
                r.me,
                r.mae,
                r.rmse,
                r.err_min,
                r.err_max,
                r.pearson,
                is_best
            ));
        }

        let out_path = match output {
            Some(p) if p.to_ascii_lowercase().ends_with(".csv") => {
                write_text_output(&csv, p)?;
                p.to_string()
            }
            Some(p) => write_or_store_layer(table, Some(p))?,
            None => write_or_store_layer(table, None)?,
        };

        // ── Optional: winner prediction raster (fit on all samples) ──────────
        let mut raster_path: Option<String> = None;
        let mut grid_meta: Option<(usize, usize, f64)> = None;
        if output_raster.is_some() {
            let (rp, rows, cols, cell) = build_winner_raster(
                best_method,
                &sx,
                &sy,
                &sv,
                prm.power,
                cx,
                cy,
                scale,
                &layer,
                prm.cell_size,
                output_raster,
            )?;
            raster_path = Some(rp);
            grid_meta = Some((rows, cols, cell));
        }

        // ── Result payload ───────────────────────────────────────────────────
        let table_json: Vec<Value> = order
            .iter()
            .map(|&idx| {
                let r = &results[idx];
                json!({
                    "rank": rank_of[idx],
                    "method": r.method.id(),
                    "n": r.n,
                    "me": r.me,
                    "mae": r.mae,
                    "rmse": r.rmse,
                    "err_min": r.err_min,
                    "err_max": r.err_max,
                    "pearson_r": r.pearson,
                    "best": idx == best,
                })
            })
            .collect();

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("observations".to_string(), json!(n));
        outputs.insert("criterion".to_string(), json!(prm.criterion.id()));
        outputs.insert("methods_evaluated".to_string(), json!(results.len()));
        outputs.insert("best_method".to_string(), json!(best_method.id()));
        outputs.insert("best_rmse".to_string(), json!(results[best].rmse));
        outputs.insert("best_mae".to_string(), json!(results[best].mae));
        outputs.insert("best_me".to_string(), json!(results[best].me));
        outputs.insert("table".to_string(), json!(table_json));
        if let Some(rp) = raster_path {
            outputs.insert("output_raster".to_string(), json!(rp));
            if let Some((rows, cols, cell)) = grid_meta {
                outputs.insert("raster_rows".to_string(), json!(rows));
                outputs.insert("raster_cols".to_string(), json!(cols));
                outputs.insert("raster_cell_size".to_string(), json!(cell));
            }
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── Cross-validation prediction ─────────────────────────────────────────────

struct MethodResult {
    method: Method,
    n: usize,
    me: f64,
    mae: f64,
    rmse: f64,
    err_min: f64,
    err_max: f64,
    pearson: f64,
}

/// Predicts the value at sample `target` from the other samples.
///
/// When `loo` is true the target sample is excluded (leave-one-out); when
/// false it is included (used when evaluating the fitted surface on a grid,
/// where `target` is a synthetic query index equal to `n`). Returns `None`
/// only if a trend fit is under-determined.
#[allow(clippy::too_many_arguments)]
fn predict(
    m: Method,
    target: usize,
    sx: &[f64],
    sy: &[f64],
    sv: &[f64],
    power: f64,
    cx: f64,
    cy: f64,
    scale: f64,
    loo: bool,
) -> Option<f64> {
    let qx = sx[target];
    let qy = sy[target];
    match m {
        Method::Idw => {
            let mut num = 0.0;
            let mut den = 0.0;
            for j in 0..sx.len() {
                if loo && j == target {
                    continue;
                }
                let dx = sx[j] - qx;
                let dy = sy[j] - qy;
                let d2 = dx * dx + dy * dy;
                if d2 <= 0.0 {
                    // Coincident sample: exact hit.
                    return Some(sv[j]);
                }
                let w = 1.0 / d2.powf(power * 0.5);
                num += w * sv[j];
                den += w;
            }
            if den > 0.0 {
                Some(num / den)
            } else {
                None
            }
        }
        Method::Nearest => {
            let mut best_d = f64::INFINITY;
            let mut best_v = None;
            for j in 0..sx.len() {
                if loo && j == target {
                    continue;
                }
                let dx = sx[j] - qx;
                let dy = sy[j] - qy;
                let d2 = dx * dx + dy * dy;
                if d2 < best_d {
                    best_d = d2;
                    best_v = Some(sv[j]);
                }
            }
            best_v
        }
        Method::Trend1 | Method::Trend2 => {
            let terms = m.trend_terms();
            let coeffs = fit_trend(m, target, sx, sy, sv, cx, cy, scale, loo)?;
            let basis = trend_basis(terms, (qx - cx) / scale, (qy - cy) / scale);
            let mut z = 0.0;
            for (c, b) in coeffs.iter().zip(basis.iter()) {
                z += c * b;
            }
            Some(z)
        }
    }
}

/// The polynomial basis vector for a normalized coordinate.
fn trend_basis(terms: usize, x: f64, y: f64) -> Vec<f64> {
    if terms == 3 {
        vec![1.0, x, y]
    } else {
        vec![1.0, x, y, x * x, x * y, y * y]
    }
}

/// Fits a trend surface (ordinary least squares) via the normal equations
/// `(XᵀX) c = Xᵀz`, optionally excluding `target`. Returns `None` when
/// under-determined or singular.
#[allow(clippy::too_many_arguments)]
fn fit_trend(
    m: Method,
    target: usize,
    sx: &[f64],
    sy: &[f64],
    sv: &[f64],
    cx: f64,
    cy: f64,
    scale: f64,
    loo: bool,
) -> Option<Vec<f64>> {
    let terms = m.trend_terms();
    let mut ata = vec![vec![0.0_f64; terms]; terms];
    let mut atb = vec![0.0_f64; terms];
    let mut used = 0usize;
    for j in 0..sx.len() {
        if loo && j == target {
            continue;
        }
        let b = trend_basis(terms, (sx[j] - cx) / scale, (sy[j] - cy) / scale);
        for r in 0..terms {
            for c in 0..terms {
                ata[r][c] += b[r] * b[c];
            }
            atb[r] += b[r] * sv[j];
        }
        used += 1;
    }
    if used < terms {
        return None;
    }
    solve_linear(ata, atb)
}

/// Solves a small dense linear system by Gaussian elimination with partial
/// pivoting. Returns `None` if the matrix is singular.
#[allow(clippy::needless_range_loop)]
fn solve_linear(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        // Partial pivot.
        let mut piv = col;
        let mut best = a[col][col].abs();
        for r in (col + 1)..n {
            if a[r][col].abs() > best {
                best = a[r][col].abs();
                piv = r;
            }
        }
        if best < 1e-12 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        let pivot = a[col][col];
        for r in (col + 1)..n {
            let f = a[r][col] / pivot;
            if f == 0.0 {
                continue;
            }
            for c in col..n {
                a[r][c] -= f * a[col][c];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for c in (i + 1)..n {
            s -= a[i][c] * x[c];
        }
        x[i] = s / a[i][i];
    }
    Some(x)
}

/// Centre/scale so normalized coordinates are ~O(1) (conditioning for the
/// trend solve). Scale is half the larger coordinate span (>= 1e-9).
fn normalize_params(sx: &[f64], sy: &[f64]) -> (f64, f64, f64) {
    let n = sx.len() as f64;
    let cx = sx.iter().sum::<f64>() / n;
    let cy = sy.iter().sum::<f64>() / n;
    let (mut xmin, mut xmax, mut ymin, mut ymax) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for i in 0..sx.len() {
        xmin = xmin.min(sx[i]);
        xmax = xmax.max(sx[i]);
        ymin = ymin.min(sy[i]);
        ymax = ymax.max(sy[i]);
    }
    let span = (xmax - xmin).max(ymax - ymin);
    let scale = if span.is_finite() && span > 0.0 {
        0.5 * span
    } else {
        1.0
    };
    (cx, cy, scale.max(1e-9))
}

// ── Winner raster ────────────────────────────────────────────────────────────

/// Fits `m` on all samples and evaluates it over a grid spanning the sample
/// bounding box. Returns `(path, rows, cols, cell)`.
#[allow(clippy::too_many_arguments)]
fn build_winner_raster(
    m: Method,
    sx: &[f64],
    sy: &[f64],
    sv: &[f64],
    power: f64,
    cx: f64,
    cy: f64,
    scale: f64,
    layer: &Layer,
    cell_size: Option<f64>,
    output: Option<&str>,
) -> Result<(String, usize, usize, f64), ToolError> {
    let (mut x_min, mut y_min, mut x_max, mut y_max) = (
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );
    for i in 0..sx.len() {
        x_min = x_min.min(sx[i]);
        y_min = y_min.min(sy[i]);
        x_max = x_max.max(sx[i]);
        y_max = y_max.max(sy[i]);
    }
    let width = (x_max - x_min).max(0.0);
    let height = (y_max - y_min).max(0.0);
    let span = width.max(height);
    if span <= 0.0 {
        return Err(ToolError::Execution(
            "all sample points are coincident; cannot build a winner raster".to_string(),
        ));
    }
    let cell = match cell_size {
        Some(c) => c,
        None => span / 256.0,
    };
    if !(cell.is_finite() && cell > 0.0) {
        return Err(ToolError::Validation(
            "'cell_size' must be a positive number".to_string(),
        ));
    }
    let ox = x_min - cell;
    let oy_top = y_max + cell;
    let cols = (((width + 2.0 * cell) / cell).ceil() as usize).max(1);
    let rows = (((height + 2.0 * cell) / cell).ceil() as usize).max(1);
    if cols > MAX_DIM || rows > MAX_DIM {
        return Err(ToolError::Validation(format!(
            "winner grid {rows}x{cols} exceeds the {MAX_DIM} cap; increase 'cell_size'"
        )));
    }

    // Pre-fit the trend once for the whole grid (idw/nearest need no fit).
    let terms = m.trend_terms();
    let trend_coeffs = if terms > 0 {
        // loo=false so `target` is never excluded; the dummy 0 is unused.
        Some(
            fit_trend(m, 0, sx, sy, sv, cx, cy, scale, false).ok_or_else(|| {
                ToolError::Execution(format!("trend fit for '{}' is singular", m.id()))
            })?,
        )
    } else {
        None
    };

    let mut data = vec![OUT_NODATA; rows * cols];
    // For idw/nearest, append a synthetic query slot so `predict` can reuse the
    // sample arrays with loo=false.
    let mut qx = sx.to_vec();
    let mut qy = sy.to_vec();
    qx.push(0.0);
    qy.push(0.0);
    let qidx = sx.len();
    for r in 0..rows {
        let wy = oy_top - (r as f64 + 0.5) * cell;
        for c in 0..cols {
            let wx = ox + (c as f64 + 0.5) * cell;
            let v = match m {
                Method::Trend1 | Method::Trend2 => {
                    let basis = trend_basis(terms, (wx - cx) / scale, (wy - cy) / scale);
                    let coeffs = trend_coeffs.as_ref().unwrap();
                    let mut z = 0.0;
                    for (cf, b) in coeffs.iter().zip(basis.iter()) {
                        z += cf * b;
                    }
                    Some(z)
                }
                _ => {
                    qx[qidx] = wx;
                    qy[qidx] = wy;
                    predict(m, qidx, &qx, &qy, sv, power, cx, cy, scale, false)
                }
            };
            if let Some(z) = v {
                if z.is_finite() {
                    data[r * cols + c] = z;
                }
            }
        }
    }

    let mut out = Raster::new(RasterConfig {
        cols,
        rows,
        bands: 1,
        x_min: ox,
        y_min: oy_top - rows as f64 * cell,
        cell_size: cell,
        cell_size_y: Some(cell),
        nodata: OUT_NODATA,
        data_type: DataType::F32,
        crs: CrsInfo {
            epsg: layer.crs_epsg(),
            wkt: None,
            proj4: None,
        },
        metadata: Vec::new(),
    });
    for r in 0..rows {
        for c in 0..cols {
            out.set(0, r as isize, c as isize, data[r * cols + c])
                .map_err(|e| ToolError::Execution(format!("failed writing cell: {e}")))?;
        }
    }
    let path = write_or_store_output(out, output)?;
    Ok((path, rows, cols, cell))
}

// ── Geometry / parameters ────────────────────────────────────────────────────

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        _ => None,
    }
}

struct Params {
    methods: Vec<Method>,
    criterion: Criterion,
    power: f64,
    cell_size: Option<f64>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let methods = match args.get("methods").and_then(Value::as_str).map(str::trim) {
        None | Some("") => ALL_METHODS.to_vec(),
        Some(list) => {
            let mut out: Vec<Method> = Vec::new();
            for tok in list.split(['|', ',', ';']) {
                let tok = tok.trim();
                if tok.is_empty() {
                    continue;
                }
                let m = Method::parse(tok).ok_or_else(|| {
                    ToolError::Validation(format!(
                        "unknown method '{tok}'; expected any of idw, nearest, trend1, trend2"
                    ))
                })?;
                if !out.contains(&m) {
                    out.push(m);
                }
            }
            if out.is_empty() {
                return Err(ToolError::Validation(
                    "'methods' listed no valid method".to_string(),
                ));
            }
            out
        }
    };
    let criterion = match args.get("criterion").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("rmse") => Criterion::Rmse,
        Some("mae") => Criterion::Mae,
        Some("me") => Criterion::Me,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'criterion' must be 'rmse', 'mae', or 'me', got '{o}'"
            )))
        }
    };
    let power = opt_pos(args, "power")?.unwrap_or(2.0);
    let cell_size = opt_pos(args, "cell_size")?;
    Ok(Params {
        methods,
        criterion,
        power,
        cell_size,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
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

fn opt_pos(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match opt_f64(args, key)? {
        Some(v) if v > 0.0 && v.is_finite() => Ok(Some(v)),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a positive number"
        ))),
        None => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn point_layer(rows: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (x, y, v) in rows {
            l.add_feature(Some(Geometry::point(*x, *y)), &[("v", (*v).into())])
                .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        ExploratoryInterpolationTool.run(&args, &ctx()).unwrap()
    }

    /// A perfectly planar field z = 3 + 2x - y must be recovered exactly by the
    /// first-order trend surface (LOOCV RMSE ~ 0), so trend1 ranks first and its
    /// RMSE beats IDW and nearest on the same data.
    #[test]
    fn planar_field_ranks_trend_first() {
        let mut rows = Vec::new();
        for gx in 0..6 {
            for gy in 0..6 {
                let x = gx as f64 * 10.0;
                let y = gy as f64 * 10.0;
                rows.push((x, y, 3.0 + 2.0 * x - y));
            }
        }
        let pts = point_layer(&rows);
        let out = run(json!({ "input": pts, "field": "v", "criterion": "rmse" }));
        assert_eq!(out.outputs["best_method"], json!("trend1"));
        // trend1's LOOCV RMSE on an exactly planar field is numerically ~0.
        let table = out.outputs["table"].as_array().unwrap();
        let trend1 = table.iter().find(|r| r["method"] == "trend1").unwrap();
        assert!(
            trend1["rmse"].as_f64().unwrap() < 1e-6,
            "planar field -> trend1 RMSE ~0, got {}",
            trend1["rmse"]
        );
        // Its Pearson r with observations is ~1.
        assert!(trend1["pearson_r"].as_f64().unwrap() > 0.999);
        // Rank 1 row is trend1.
        assert_eq!(table[0]["method"], json!("trend1"));
        assert_eq!(table[0]["rank"], json!(1));
    }

    /// A method subset restricts the compared methods; the LOOCV count equals
    /// the sample total for every method.
    #[test]
    fn method_subset_and_counts() {
        let mut rows = Vec::new();
        for i in 0..12 {
            let x = i as f64;
            rows.push((x, 0.5 * x, x * x)); // curved field
        }
        let pts = point_layer(&rows);
        let out = run(json!({
            "input": pts, "field": "v", "methods": "idw,nearest"
        }));
        assert_eq!(out.outputs["methods_evaluated"], json!(2));
        let table = out.outputs["table"].as_array().unwrap();
        assert_eq!(table.len(), 2);
        for r in table {
            assert_eq!(r["n"], json!(12));
            assert!(r["rmse"].as_f64().unwrap() >= 0.0);
        }
        // Only the two requested methods appear.
        let names: Vec<&str> = table
            .iter()
            .map(|r| r["method"].as_str().unwrap())
            .collect();
        assert!(names.contains(&"idw"));
        assert!(names.contains(&"nearest"));
        assert!(!names.contains(&"trend1"));
    }

    /// The `me` criterion ranks by least |bias|, and the winner raster is
    /// written when requested. On a quadratic field the second-order trend
    /// should win and its surface should match the truth near a sample.
    #[test]
    fn me_criterion_and_winner_raster() {
        let mut rows = Vec::new();
        for gx in 0..6 {
            for gy in 0..6 {
                let x = gx as f64;
                let y = gy as f64;
                rows.push((x, y, 5.0 + x * x + y * y + x * y));
            }
        }
        let pts = point_layer(&rows);
        let tif =
            std::env::temp_dir().join(format!("exp_interp_winner_{}.tif", std::process::id()));
        let tif = tif.to_str().unwrap().to_string();
        let out = run(json!({
            "input": pts, "field": "v", "criterion": "me",
            "output_raster": tif,
        }));
        // trend2 recovers a quadratic field -> least bias.
        assert_eq!(out.outputs["best_method"], json!("trend2"));
        assert!(out.outputs.contains_key("output_raster"));
        let rpath = out.outputs["output_raster"].as_str().unwrap();
        let r = crate::common::load_input_raster(rpath).unwrap();
        // Sample the winner raster near (2,3) and compare to the true surface.
        let col = ((2.0 - r.x_min) / r.cell_size_x).floor() as isize;
        let y_top = r.y_min + r.rows as f64 * r.cell_size_y;
        let row = ((y_top - 3.0) / r.cell_size_y).floor() as isize;
        let v = r.get(0, row, col);
        let truth = 5.0 + 2.0 * 2.0 + 3.0 * 3.0 + 2.0 * 3.0;
        assert!(
            (v - truth).abs() < 1.0,
            "winner raster near (2,3) ~= {truth}, got {v}"
        );
        let _ = std::fs::remove_file(&tif);
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ExploratoryInterpolationTool.validate(&args)
        };
        assert!(bad(json!({})).is_err()); // no input
        assert!(bad(json!({ "input": "a.geojson" })).is_err()); // no field
        assert!(bad(json!({ "input": "a.geojson", "field": "v", "methods": "spline" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "field": "v", "criterion": "aic" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "field": "v", "power": -1 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "field": "v" })).is_ok());
        assert!(
            bad(json!({ "input": "a.geojson", "field": "v", "methods": "idw|trend2" })).is_ok()
        );
    }
}
