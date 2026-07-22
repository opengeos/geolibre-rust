//! GeoLibre tool: multiscale geographically weighted regression (MGWR).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Multiscale Geographically Weighted
//! Regression* (Spatial Statistics). Standard GWR (`geographically_weighted_regression`)
//! fits a local linear regression at every feature, but forces every term
//! (intercept and every explanatory variable) to share one bandwidth — one
//! spatial scale for every process. In reality some relationships are highly
//! local (a short-range amenity effect) while others are near-constant across
//! the whole study area (a regional or global driver). MGWR lets **each term
//! find its own optimal bandwidth**, via back-fitting (Fotheringham, Yang &
//! Kalogirou 2017):
//!
//! 1. **Warm start** — fit ordinary single-bandwidth GWR (reusing the same
//!    kernel-weighted least-squares core and kdtree-free O(n) neighbour sweep
//!    as `geographically_weighted_regression`) to get initial per-location
//!    coefficients and an initial shared bandwidth.
//! 2. **Back-fit** — repeat until convergence: for every term `k` in turn,
//!    compute the *partial residual* (the response with every other term's
//!    current local contribution removed), then re-estimate `k`'s local
//!    coefficient by a **univariate** weighted regression of the partial
//!    residual on term `k`'s column, searching for `k`'s own bandwidth by
//!    **golden-section search minimizing AICc**. Convergence is judged by the
//!    relative change in residual sum of squares (RSS) between sweeps
//!    (Fotheringham et al.'s "score of change").
//!
//! Terms whose relationship barely varies across space settle on a **large**
//! bandwidth (approaching the whole study area — "regional"/"global" scale);
//! terms with a strongly local relationship settle on a **small** one
//! ("local" scale). Each location's fit (the inner `i in 0..n` sweep inside
//! `univariate_fit`/`full_fit`) is independent of every other location's, so
//! the per-feature work is embarrassingly parallel; it runs sequentially here
//! to stay dependency-free and portable to WASM (no thread pool / rayon).
//!
//! Output: per-feature local coefficients (one field per term), `predicted`,
//! `residual`, `local_r2`, and `condition_number` (a local multicollinearity
//! diagnostic, both computed with the warm-start reference bandwidth's
//! neighbourhood). The tool report carries each term's optimal bandwidth with
//! a local/regional/global interpretation, backfitting convergence status,
//! effective parameters (sum of each term's hat-matrix trace), and global
//! AICc / R² — plus the warm-start single-bandwidth AICc for direct
//! comparison (MGWR should never do worse).
//!
//! Linear algebra (Gaussian elimination with partial pivoting) and the kernel
//! machinery are hand-rolled and self-contained here — no new crate, and no
//! reach into the sibling `geographically_weighted_regression` module.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Upper bound on the adaptive-bandwidth neighbour count searched during
/// back-fitting, to keep the O(n) per-location neighbour sweep bounded on
/// large inputs. Never applied silently — a capped run logs it once.
const MAX_ADAPTIVE_NEIGHBORS: usize = 400;

pub struct MgwrTool;

impl Tool for MgwrTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "mgwr",
            display_name: "Multiscale Geographically Weighted Regression",
            summary: "GWR where each explanatory term gets its own bandwidth via AICc back-fitting (MGWR): per-feature coefficients, local R2, residuals, condition number, plus each term's optimal bandwidth with a local/regional/global interpretation and global AICc / R2 diagnostics.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input feature layer (points, or other geometries via their representative point).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "dependent_field",
                    description: "Dependent (response) numeric field.",
                    required: true,
                },
                ToolParamSpec {
                    name: "explanatory_fields",
                    description: "Comma-separated explanatory numeric field(s), each fit with its own bandwidth.",
                    required: true,
                },
                ToolParamSpec {
                    name: "kernel",
                    description: "Distance-decay kernel: 'gaussian' (default) or 'bisquare'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "bandwidth_type",
                    description: "'adaptive' (default; bandwidth = k nearest neighbours) or 'fixed' (bandwidth = a distance). Applies to every term's own bandwidth search.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Back-fitting convergence tolerance: stop when the relative change in RSS between sweeps drops below this (default 0.001).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_iterations",
                    description: "Maximum back-fitting sweeps over all terms (default 30).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["input", "dependent_field", "explanatory_fields"] {
            if args
                .get(key)
                .and_then(Value::as_str)
                .map(str::trim)
                .unwrap_or("")
                .is_empty()
            {
                return Err(ToolError::Validation(format!(
                    "missing required string parameter '{key}'"
                )));
            }
        }
        parse_params(args)?;
        Ok(())
    }

    // Range loops index parallel matrices/vectors (design rows, per-term
    // coefficient columns); the index form is clearer than zipping several
    // slices here.
    #[allow(clippy::needless_range_loop)]
    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let y_field = require_str(args, "dependent_field")?;
        let x_fields: Vec<String> = require_str(args, "explanatory_fields")?
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(str::to_string)
            .collect();
        if x_fields.is_empty() {
            return Err(ToolError::Validation(
                "'explanatory_fields' must list at least one explanatory field".to_string(),
            ));
        }
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let schema = layer.schema.clone();

        // Collect locations and the design matrix (intercept + explanatory).
        let mut locs: Vec<(f64, f64)> = Vec::new();
        let mut xs: Vec<Vec<f64>> = Vec::new();
        let mut ys: Vec<f64> = Vec::new();
        let mut idx_map: Vec<usize> = Vec::new();
        for (fi, feature) in layer.features.iter().enumerate() {
            let Some((px, py)) = feature.geometry.as_ref().and_then(rep_point) else {
                continue;
            };
            let y = feature
                .get(&schema, y_field)
                .ok()
                .and_then(FieldValue::as_f64);
            let mut row = Vec::with_capacity(x_fields.len() + 1);
            row.push(1.0); // intercept term
            let mut ok = y.is_some();
            for xf in &x_fields {
                match feature.get(&schema, xf).ok().and_then(FieldValue::as_f64) {
                    Some(v) if v.is_finite() => row.push(v),
                    _ => ok = false,
                }
            }
            if let (true, Some(y)) = (ok && y.unwrap().is_finite(), y) {
                locs.push((px, py));
                xs.push(row);
                ys.push(y);
                idx_map.push(fi);
            }
        }
        let n = ys.len();
        let p = x_fields.len() + 1;
        if n <= p + 2 {
            return Err(ToolError::Execution(format!(
                "need more than {} observations with valid y and x values, found {n}",
                p + 2
            )));
        }

        // Mean-centre every explanatory column (not the intercept). A
        // nonzero-mean regressor is structurally collinear with the
        // intercept under any local window — the intercept can always
        // partly compensate for a mis-estimated local slope by shifting —
        // a well-documented source of local multicollinearity in GWR-family
        // models. Left uncentred, that collinearity can destabilise
        // back-fitting badly: terms "collude" to explain each other's
        // leftover local bias with an ever-smaller neighbourhood instead of
        // settling at their true scale. Centring is a pure reparametrisation
        // (it changes nothing about the fitted values), so the per-location
        // slope coefficients back-fitting produces are already on the
        // original scale; only the intercept needs translating back at the
        // end.
        let col_mean: Vec<f64> = (0..p)
            .map(|c| {
                if c == 0 {
                    0.0
                } else {
                    xs.iter().map(|row| row[c]).sum::<f64>() / n as f64
                }
            })
            .collect();
        for row in &mut xs {
            for c in 1..p {
                row[c] -= col_mean[c];
            }
        }

        let data = Data { locs, xs, ys, p };
        let diag = bbox_diag(&data.locs);

        let neighbor_cap = if n > MAX_ADAPTIVE_NEIGHBORS {
            ctx.progress.info(&format!(
                "adaptive neighbour count capped at {MAX_ADAPTIVE_NEIGHBORS} (of {n} observations) to bound back-fitting cost"
            ));
            MAX_ADAPTIVE_NEIGHBORS
        } else {
            n
        };

        // ── Warm start: single shared bandwidth (same idea as GWR) — used
        // only as the AICc baseline MGWR is compared against, and as the
        // reference neighbourhood for the local R2 / condition-number
        // diagnostics. ──
        ctx.progress
            .info("MGWR: fitting single-bandwidth warm start");
        let (init_bw, init_fit) =
            optimize_full_bandwidth(&data, prm.kernel, prm.adaptive, neighbor_cap).ok_or_else(
                || {
                    ToolError::Execution(
                        "initial single-bandwidth GWR fit failed (try a different bandwidth_type)"
                            .to_string(),
                    )
                },
            )?;
        ctx.progress.info(&format!(
            "warm start bandwidth {init_bw:.4}, AICc {:.2}",
            init_fit.aicc
        ));

        // ── Back-fitting: each term finds its own bandwidth. ──
        // Back-fitting starts from the *global* OLS coefficients (constant
        // across every location), not the single-bandwidth GWR surface:
        // starting from an already spatially-varying fit would feed every
        // term's very first partial residual a locally noisy signal, and
        // the AICc search would then favour the smallest allowed
        // neighbourhood for *every* term regardless of its true scale. The
        // flat OLS start (the standard MGWR initialisation, Fotheringham et
        // al. 2017) only lets a term shrink its bandwidth once the other
        // terms' contributions are themselves reasonably stable.
        let ols_beta = ols_fit(&data).ok_or_else(|| {
            ToolError::Execution("initial OLS fit failed (design matrix is singular)".to_string())
        })?;
        let mut beta = vec![ols_beta; n]; // n x p, per-location coefficients
        let mut bandwidths = vec![init_bw; p];
        let mut term_trace = vec![init_fit.trace / p as f64; p];
        // Search over the same neighbour-count range the single-bandwidth
        // warm start used (`p + 2` up to the capped sample size): giving
        // MGWR's per-term search a *smaller* range than the shared-bandwidth
        // baseline would artificially handicap it in the AICc comparison
        // below (it could never reproduce the baseline's own best fit, let
        // alone improve on it by letting a term deviate from it).
        let (lo_adaptive, hi_adaptive) = ((p + 2) as f64, neighbor_cap.max(p + 2) as f64);
        let (lo_fixed, hi_fixed) = (
            (diag / (n as f64).sqrt() * 0.1).max(f64::MIN_POSITIVE),
            diag.max(f64::MIN_POSITIVE),
        );

        let mut prev_rss: f64 = (0..n)
            .map(|i| {
                let pred: f64 = (0..p).map(|c| data.xs[i][c] * beta[i][c]).sum();
                (data.ys[i] - pred).powi(2)
            })
            .sum();
        let mut converged = false;
        let mut iterations_done = 0usize;

        for iter in 0..prm.max_iterations {
            for k in 0..p {
                let partial: Vec<f64> = (0..n)
                    .map(|i| {
                        let mut r = data.ys[i];
                        for j in 0..p {
                            if j != k {
                                r -= data.xs[i][j] * beta[i][j];
                            }
                        }
                        r
                    })
                    .collect();
                let (lo, hi) = if prm.adaptive {
                    (lo_adaptive, hi_adaptive)
                } else {
                    (lo_fixed, hi_fixed)
                };
                // Score each candidate bandwidth for term k by the AICc of
                // the *whole* additive model — this term's candidate fit
                // plus every other term held at its current estimate — not
                // just this term's own isolated univariate fit. Scoring in
                // isolation would let every term chase the smallest
                // bandwidth that best explains its own partial residual
                // (which, early in back-fitting, still carries the other
                // terms' not-yet-converged local noise), collapsing every
                // term onto a near-interpolating neighbourhood. Charging
                // each candidate against the model's *cumulative* effective
                // parameters (every other term's current trace + this
                // one's) is what lets an already-flexible model push a
                // genuinely global term like a constant coefficient back
                // out to a large bandwidth.
                let other_trace: f64 = (0..p).filter(|&j| j != k).map(|j| term_trace[j]).sum();
                let coarse_steps = if prm.adaptive {
                    18.min(hi_adaptive as usize - lo_adaptive as usize).max(1)
                } else {
                    18
                };
                let found = search_bandwidth(lo, hi, coarse_steps, |b| {
                    let fit = univariate_fit(
                        &data,
                        k,
                        &partial,
                        b,
                        prm.kernel,
                        prm.adaptive,
                        neighbor_cap,
                    )?;
                    let rss_full: f64 = (0..n)
                        .map(|i| (partial[i] - data.xs[i][k] * fit.coef[i]).powi(2))
                        .sum();
                    let aicc = model_aicc(n, rss_full, other_trace + fit.trace);
                    aicc.is_finite().then_some(aicc)
                });
                let Some((best_b, _)) = found else {
                    continue; // keep this term's previous bandwidth/coefficients
                };
                if let Some(fit_k) = univariate_fit(
                    &data,
                    k,
                    &partial,
                    best_b,
                    prm.kernel,
                    prm.adaptive,
                    neighbor_cap,
                ) {
                    bandwidths[k] = if prm.adaptive {
                        best_b.round().clamp(lo_adaptive, hi_adaptive)
                    } else {
                        best_b
                    };
                    for i in 0..n {
                        beta[i][k] = fit_k.coef[i];
                    }
                    term_trace[k] = fit_k.trace;
                }
            }
            let rss_now: f64 = (0..n)
                .map(|i| {
                    let pred: f64 = (0..p).map(|c| data.xs[i][c] * beta[i][c]).sum();
                    (data.ys[i] - pred).powi(2)
                })
                .sum();
            iterations_done = iter + 1;
            let soc = (prev_rss - rss_now).abs() / rss_now.max(1e-12);
            ctx.progress.info(&format!(
                "back-fit sweep {iterations_done}: rss {rss_now:.6}, score of change {soc:.6}"
            ));
            prev_rss = rss_now;
            if soc < prm.tolerance {
                converged = true;
                break;
            }
        }
        if !converged {
            ctx.progress.info(&format!(
                "back-fitting stopped at max_iterations ({}) without reaching tolerance {:.6}",
                prm.max_iterations, prm.tolerance
            ));
        }

        // ── Final diagnostics from the converged per-term coefficients. ──
        let mut predicted = vec![0.0; n];
        let mut residual = vec![0.0; n];
        for i in 0..n {
            let pred: f64 = (0..p).map(|c| data.xs[i][c] * beta[i][c]).sum();
            predicted[i] = pred;
            residual[i] = data.ys[i] - pred;
        }
        let rss: f64 = residual.iter().map(|r| r * r).sum();
        let ybar = data.ys.iter().sum::<f64>() / n as f64;
        let tss: f64 = data.ys.iter().map(|y| (y - ybar).powi(2)).sum();
        let r2 = if tss > 0.0 { 1.0 - rss / tss } else { 0.0 };
        let trace_s: f64 = term_trace.iter().sum();
        let nf = n as f64;
        let aicc = model_aicc(n, rss, trace_s);
        let adj_r2 = if nf - trace_s - 1.0 > 0.0 {
            1.0 - (1.0 - r2) * (nf - 1.0) / (nf - trace_s - 1.0)
        } else {
            r2
        };

        // Local R2 and condition number: neighbourhood defined by the
        // warm-start reference bandwidth, evaluated against the final
        // back-fit residuals (not re-fit per term).
        let (local_r2, condition_number) = local_diagnostics(
            &data,
            &residual,
            init_bw,
            prm.kernel,
            prm.adaptive,
            neighbor_cap,
        );

        // ── Write output fields. ──
        let mut field_names: Vec<String> = vec![
            "predicted".into(),
            "residual".into(),
            "local_r2".into(),
            "condition_number".into(),
        ];
        field_names.push("b_intercept".into());
        for xf in &x_fields {
            field_names.push(format!("b_{xf}"));
        }
        for name in &field_names {
            layer.add_field(FieldDef::new(name.clone(), FieldType::Float));
        }
        let mut row_for_feature: Vec<Option<usize>> = vec![None; layer.features.len()];
        for (row, &fi) in idx_map.iter().enumerate() {
            row_for_feature[fi] = Some(row);
        }
        let extra = field_names.len();
        for (fi, feature) in layer.features.iter_mut().enumerate() {
            match row_for_feature[fi] {
                Some(r) => {
                    feature.attributes.push(FieldValue::Float(predicted[r]));
                    feature.attributes.push(FieldValue::Float(residual[r]));
                    feature.attributes.push(FieldValue::Float(local_r2[r]));
                    feature
                        .attributes
                        .push(FieldValue::Float(condition_number[r]));
                    // Translate the intercept back from the centred fit to
                    // the original explanatory-variable scale; the slope
                    // coefficients (c >= 1) are unaffected by centring.
                    let intercept_true =
                        beta[r][0] - (1..p).map(|c| beta[r][c] * col_mean[c]).sum::<f64>();
                    feature.attributes.push(FieldValue::Float(intercept_true));
                    for c in 1..p {
                        feature.attributes.push(FieldValue::Float(beta[r][c]));
                    }
                }
                None => {
                    for _ in 0..extra {
                        feature.attributes.push(FieldValue::Null);
                    }
                }
            }
        }

        ctx.progress.info(&format!(
            "MGWR: R2 {r2:.4}, adjusted R2 {adj_r2:.4}, AICc {aicc:.2} (single-bandwidth warm start AICc {:.2}), effective params {trace_s:.2}",
            init_fit.aicc
        ));

        let out_path = write_or_store_layer(layer, output)?;

        let mut term_names: Vec<String> = vec!["intercept".to_string()];
        term_names.extend(x_fields.iter().cloned());
        let variable_bandwidths: Vec<Value> = term_names
            .iter()
            .zip(bandwidths.iter())
            .map(|(name, &b)| {
                json!({
                    "field": name,
                    "bandwidth": b,
                    "scale": interpret_scale(b, prm.adaptive, n, diag),
                })
            })
            .collect();

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("observations".to_string(), json!(n));
        outputs.insert("terms".to_string(), json!(p));
        outputs.insert("kernel".to_string(), json!(prm.kernel.as_str()));
        outputs.insert(
            "bandwidth_type".to_string(),
            json!(if prm.adaptive { "adaptive" } else { "fixed" }),
        );
        outputs.insert(
            "variable_bandwidths".to_string(),
            json!(variable_bandwidths),
        );
        outputs.insert("iterations".to_string(), json!(iterations_done));
        outputs.insert("converged".to_string(), json!(converged));
        outputs.insert("tolerance".to_string(), json!(prm.tolerance));
        outputs.insert("max_iterations".to_string(), json!(prm.max_iterations));
        outputs.insert("effective_params".to_string(), json!(trace_s));
        outputs.insert("r2".to_string(), json!(r2));
        outputs.insert("adjusted_r2".to_string(), json!(adj_r2));
        outputs.insert("aicc".to_string(), json!(aicc));
        outputs.insert("residual_ss".to_string(), json!(rss));
        outputs.insert("single_bandwidth".to_string(), json!(init_bw));
        outputs.insert("single_bandwidth_aicc".to_string(), json!(init_fit.aicc));
        Ok(ToolRunResult { outputs })
    }
}

// ── MGWR core ──────────────────────────────────────────────────────────────

struct Data {
    locs: Vec<(f64, f64)>,
    xs: Vec<Vec<f64>>, // n rows, each length p (intercept first)
    ys: Vec<f64>,
    p: usize,
}

/// Plain (unweighted) global least squares: `β = (XᵀX)⁻¹ Xᵀy`. Used only to
/// give back-fitting a flat, spatially-constant starting surface.
fn ols_fit(data: &Data) -> Option<Vec<f64>> {
    let n = data.ys.len();
    let p = data.p;
    let mut a = vec![vec![0.0; p]; p];
    let mut rhs = vec![0.0; p];
    for i in 0..n {
        let xi = &data.xs[i];
        let yi = data.ys[i];
        for r in 0..p {
            for c in 0..p {
                a[r][c] += xi[r] * xi[c];
            }
            rhs[r] += xi[r] * yi;
        }
    }
    solve(&a, &rhs)
}

/// Result of a full multivariate local fit at one bandwidth (the warm start).
struct FullFit {
    trace: f64,
    aicc: f64,
}

/// Fits single-bandwidth GWR (all terms sharing `bandwidth`) — used only to
/// warm-start the back-fitting coefficients and as the AICc baseline that
/// MGWR is compared against.
// Range loops index parallel matrices/vectors (design rows, distances); the
// index form is clearer than zipping several slices here.
#[allow(clippy::needless_range_loop)]
fn full_fit(
    data: &Data,
    bandwidth: f64,
    kernel: Kernel,
    adaptive: bool,
    neighbor_cap: usize,
) -> Option<FullFit> {
    let n = data.ys.len();
    let p = data.p;
    let mut residual = vec![0.0; n];
    let mut trace = 0.0;
    let mut solved = 0usize;

    for i in 0..n {
        let (xi, yi) = data.locs[i];
        let d: Vec<f64> = (0..n)
            .map(|j| {
                let (xj, yj) = data.locs[j];
                ((xi - xj).powi(2) + (yi - yj).powi(2)).sqrt()
            })
            .collect();
        let b = if adaptive {
            let k = (bandwidth.round() as usize).clamp(p + 1, neighbor_cap.min(n));
            let mut sorted = d.clone();
            sorted.sort_by(f64::total_cmp);
            sorted[k - 1].max(f64::MIN_POSITIVE)
        } else {
            bandwidth
        };

        let mut a = vec![vec![0.0; p]; p];
        let mut rhs = vec![0.0; p];
        for j in 0..n {
            let w = kernel.weight(d[j], b);
            if w <= 0.0 {
                continue;
            }
            let xj = &data.xs[j];
            let yj = data.ys[j];
            for r in 0..p {
                let wr = w * xj[r];
                for c in 0..p {
                    a[r][c] += wr * xj[c];
                }
                rhs[r] += wr * yj;
            }
        }
        let Some(beta) = solve(&a, &rhs) else {
            continue;
        };
        let Some(z) = solve(&a, &data.xs[i]) else {
            continue;
        };
        let s_ii: f64 = (0..p).map(|c| data.xs[i][c] * z[c]).sum();
        let yhat_i: f64 = (0..p).map(|c| data.xs[i][c] * beta[c]).sum();
        residual[i] = data.ys[i] - yhat_i;
        trace += s_ii;
        solved += 1;
    }

    if solved == 0 {
        return None;
    }
    let rss: f64 = residual.iter().map(|r| r * r).sum();
    let aicc = model_aicc(n, rss, trace);
    Some(FullFit { trace, aicc })
}

/// AICc-minimizing shared bandwidth over a coarse candidate grid — the same
/// approach `geographically_weighted_regression` uses when no bandwidth is
/// supplied, reused here purely to warm-start back-fitting.
fn optimize_full_bandwidth(
    data: &Data,
    kernel: Kernel,
    adaptive: bool,
    neighbor_cap: usize,
) -> Option<(f64, FullFit)> {
    let n = data.ys.len();
    let p = data.p;
    let candidates: Vec<f64> = if adaptive {
        let lo = p + 2;
        let hi = neighbor_cap.max(lo);
        let steps = 18.min(hi - lo).max(1);
        (0..=steps)
            .map(|s| (lo + s * (hi - lo) / steps) as f64)
            .collect()
    } else {
        let diag = bbox_diag(&data.locs);
        let lo = (diag / (n as f64).sqrt() * 0.5).max(f64::MIN_POSITIVE);
        let hi = diag.max(lo * 2.0);
        (0..=18).map(|s| lo + (hi - lo) * s as f64 / 18.0).collect()
    };

    let mut best: Option<(f64, FullFit)> = None;
    for &b in &candidates {
        if let Some(fit) = full_fit(data, b, kernel, adaptive, neighbor_cap) {
            if fit.aicc.is_finite() && best.as_ref().is_none_or(|(_, bf)| fit.aicc < bf.aicc) {
                best = Some((b, fit));
            }
        }
    }
    best
}

/// Result of a univariate (single-term) weighted local fit, used inside
/// back-fitting for one term against its current partial residual.
struct UniFit {
    coef: Vec<f64>,
    trace: f64,
}

/// The corrected AIC (Fotheringham et al.) of a fit summarized by its
/// residual sum of squares and effective number of parameters (hat-matrix
/// trace) over `n` observations.
fn model_aicc(n: usize, rss: f64, trace: f64) -> f64 {
    let nf = n as f64;
    let sigma = (rss / nf).sqrt().max(f64::MIN_POSITIVE);
    let denom = nf - 2.0 - trace;
    if denom > 0.0 {
        2.0 * nf * sigma.ln() + nf * (2.0 * std::f64::consts::PI).ln() + nf * (nf + trace) / denom
    } else {
        f64::INFINITY
    }
}

/// Fits term `col`'s local coefficient against `target` (the partial
/// residual with every other term's contribution removed) at bandwidth `b`:
/// a weighted regression through the origin, `β_i = Σ w x target / Σ w x²`.
#[allow(clippy::needless_range_loop)]
fn univariate_fit(
    data: &Data,
    col: usize,
    target: &[f64],
    b: f64,
    kernel: Kernel,
    adaptive: bool,
    neighbor_cap: usize,
) -> Option<UniFit> {
    let n = data.ys.len();
    let mut coef = vec![0.0; n];
    let mut trace = 0.0;
    let mut solved = 0usize;
    for i in 0..n {
        let (xi, yi) = data.locs[i];
        let d: Vec<f64> = (0..n)
            .map(|j| {
                let (xj, yj) = data.locs[j];
                ((xi - xj).powi(2) + (yi - yj).powi(2)).sqrt()
            })
            .collect();
        let bw = if adaptive {
            let k = (b.round() as usize).clamp(2, neighbor_cap.min(n));
            let mut sorted = d.clone();
            sorted.sort_by(f64::total_cmp);
            sorted[k - 1].max(f64::MIN_POSITIVE)
        } else {
            b
        };
        let mut sxx = 0.0;
        let mut sxy = 0.0;
        for j in 0..n {
            let w = kernel.weight(d[j], bw);
            if w <= 0.0 {
                continue;
            }
            let xkj = data.xs[j][col];
            sxx += w * xkj * xkj;
            sxy += w * xkj * target[j];
        }
        if sxx <= 1e-12 {
            continue;
        }
        let beta_i = sxy / sxx;
        coef[i] = beta_i;
        let xki = data.xs[i][col];
        let w_ii = kernel.weight(0.0, bw);
        trace += w_ii * xki * xki / sxx;
        solved += 1;
    }
    if solved == 0 {
        return None;
    }
    Some(UniFit { coef, trace })
}

/// Local R² and a local design-matrix condition number, both defined over
/// the neighbourhood implied by the warm-start reference bandwidth: R² uses
/// the final MGWR residuals (not a per-term re-fit); the condition number is
/// the 1-norm estimate `cond(A) = ||A||₁ · ||A⁻¹||₁` of the local weighted
/// normal-equations matrix `A = XᵀW_iX`, a standard multicollinearity
/// diagnostic for the explanatory variables near each feature.
#[allow(clippy::needless_range_loop)]
fn local_diagnostics(
    data: &Data,
    residual: &[f64],
    reference_bandwidth: f64,
    kernel: Kernel,
    adaptive: bool,
    neighbor_cap: usize,
) -> (Vec<f64>, Vec<f64>) {
    let n = data.ys.len();
    let p = data.p;
    let mut local_r2 = vec![0.0; n];
    let mut condition_number = vec![f64::NAN; n];
    for i in 0..n {
        let (xi, yi) = data.locs[i];
        let d: Vec<f64> = (0..n)
            .map(|j| {
                let (xj, yj) = data.locs[j];
                ((xi - xj).powi(2) + (yi - yj).powi(2)).sqrt()
            })
            .collect();
        let bw = if adaptive {
            let k = (reference_bandwidth.round() as usize).clamp(p + 1, neighbor_cap.min(n));
            let mut sorted = d.clone();
            sorted.sort_by(f64::total_cmp);
            sorted[k - 1].max(f64::MIN_POSITIVE)
        } else {
            reference_bandwidth
        };
        let (mut sw, mut swy, mut swy2, mut rss_w) = (0.0, 0.0, 0.0, 0.0);
        let mut a = vec![vec![0.0; p]; p];
        for j in 0..n {
            let w = kernel.weight(d[j], bw);
            if w <= 0.0 {
                continue;
            }
            sw += w;
            swy += w * data.ys[j];
            swy2 += w * data.ys[j] * data.ys[j];
            rss_w += w * residual[j] * residual[j];
            let xj = &data.xs[j];
            for r in 0..p {
                for c in 0..p {
                    a[r][c] += w * xj[r] * xj[c];
                }
            }
        }
        if sw > 0.0 {
            let ybar_w = swy / sw;
            let tss_w = swy2 - sw * ybar_w * ybar_w;
            local_r2[i] = if tss_w > 0.0 {
                (1.0 - rss_w / tss_w).clamp(0.0, 1.0)
            } else {
                0.0
            };
        }
        condition_number[i] = condition_number_1norm(&a);
    }
    (local_r2, condition_number)
}

/// 1-norm condition number estimate `||A||₁ · ||A⁻¹||₁` of a small dense
/// matrix (`A⁻¹` built column-by-column by solving `A x = e_k`).
fn condition_number_1norm(a: &[Vec<f64>]) -> f64 {
    let p = a.len();
    let mut ainv_cols: Vec<Vec<f64>> = Vec::with_capacity(p);
    for k in 0..p {
        let mut e = vec![0.0; p];
        e[k] = 1.0;
        match solve(a, &e) {
            Some(col) => ainv_cols.push(col),
            None => return f64::INFINITY,
        }
    }
    let norm1_a = (0..p)
        .map(|j| (0..p).map(|i| a[i][j].abs()).sum::<f64>())
        .fold(0.0_f64, f64::max);
    let norm1_ainv = (0..p)
        .map(|j| (0..p).map(|i| ainv_cols[j][i].abs()).sum::<f64>())
        .fold(0.0_f64, f64::max);
    norm1_a * norm1_ainv
}

/// Deterministic golden-section search for the bandwidth in `[lo, hi]`
/// minimizing `f` (here, a term's AICc). Returns `None` if `f` never returns
/// a finite value. `max_iter` bounds the number of `f` evaluations.
fn golden_section_search(
    mut lo: f64,
    mut hi: f64,
    max_iter: usize,
    stop_width: f64,
    mut f: impl FnMut(f64) -> Option<f64>,
) -> Option<(f64, f64)> {
    // Not `hi <= lo`: that would treat a NaN bound as "proceed anyway"
    // instead of falling back to the single-point evaluation below.
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    let degenerate = !(hi > lo);
    if degenerate {
        let v = f(lo)?;
        return Some((lo, v));
    }
    const GR: f64 = 0.618_033_988_749_895; // (sqrt(5)-1)/2
    let mut c = hi - GR * (hi - lo);
    let mut d = lo + GR * (hi - lo);
    let mut fc = f(c);
    let mut fd = f(d);
    let mut best: Option<(f64, f64)> = None;
    for eval in [(c, fc), (d, fd)] {
        if let (x, Some(v)) = eval {
            if best.is_none_or(|(_, bv)| v < bv) {
                best = Some((x, v));
            }
        }
    }
    for _ in 0..max_iter {
        if (hi - lo).abs() < stop_width {
            break;
        }
        let go_left = match (fc, fd) {
            (Some(vc), Some(vd)) => vc < vd,
            (Some(_), None) => true,
            (None, Some(_)) => false,
            (None, None) => break,
        };
        if go_left {
            hi = d;
            d = c;
            fd = fc;
            c = hi - GR * (hi - lo);
            fc = f(c);
            if let Some(v) = fc {
                if best.is_none_or(|(_, bv)| v < bv) {
                    best = Some((c, v));
                }
            }
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = lo + GR * (hi - lo);
            fd = f(d);
            if let Some(v) = fd {
                if best.is_none_or(|(_, bv)| v < bv) {
                    best = Some((d, v));
                }
            }
        }
    }
    best
}

/// Bandwidth search used for each term's back-fitting step: a coarse,
/// evenly-spaced scan of `[lo, hi]` followed by a golden-section refinement
/// bracketed around the best grid point. AICc surfaces here are often
/// nearly flat over wide stretches (a term whose true coefficient barely
/// varies gives almost the same fit at every large-enough bandwidth) with
/// only floating-point-scale differences between neighbouring candidates —
/// exactly the case golden-section search alone assumes cannot happen
/// (strict unimodality). A blind golden-section run can drift into one of
/// those plateaus and settle for a tie no better than the true optimum
/// while never sampling far enough away to find it; scanning the whole
/// range first guarantees the refinement starts from the right
/// neighbourhood.
fn search_bandwidth(
    lo: f64,
    hi: f64,
    coarse_steps: usize,
    mut f: impl FnMut(f64) -> Option<f64>,
) -> Option<(f64, f64)> {
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    let degenerate = !(hi > lo);
    if degenerate {
        return f(lo).map(|v| (lo, v));
    }
    let steps = coarse_steps.max(1);
    let mut candidates = Vec::with_capacity(steps + 1);
    let mut best: Option<(f64, f64)> = None;
    for s in 0..=steps {
        let b = lo + (hi - lo) * s as f64 / steps as f64;
        if let Some(v) = f(b) {
            if best.is_none_or(|(_, bv)| v < bv) {
                best = Some((b, v));
            }
        }
        candidates.push(b);
    }
    let (best_b, best_v) = best?;
    let idx = candidates
        .iter()
        .position(|&b| (b - best_b).abs() < f64::EPSILON)
        .unwrap_or(0);
    let bracket_lo = candidates[idx.saturating_sub(1)];
    let bracket_hi = candidates[(idx + 1).min(candidates.len() - 1)];
    let tol = ((bracket_hi - bracket_lo) * 1e-3).max(f64::MIN_POSITIVE);
    if let Some((rb, rv)) = golden_section_search(bracket_lo, bracket_hi, 12, tol, &mut f) {
        if rv < best_v {
            return Some((rb, rv));
        }
    }
    Some((best_b, best_v))
}

/// Classifies a term's optimal bandwidth relative to the whole study area:
/// near-global bandwidths mean the relationship barely varies across space,
/// small ones mean it is highly local. Thresholds follow the informal
/// convention used in the MGWR literature (Fotheringham et al. 2017).
fn interpret_scale(bandwidth: f64, adaptive: bool, n: usize, diag: f64) -> &'static str {
    let ratio = if adaptive {
        bandwidth / n.max(1) as f64
    } else {
        bandwidth / diag.max(f64::MIN_POSITIVE)
    };
    if ratio < 0.25 {
        "local"
    } else if ratio < 0.75 {
        "regional"
    } else {
        "global"
    }
}

fn bbox_diag(locs: &[(f64, f64)]) -> f64 {
    let (mut minx, mut miny, mut maxx, mut maxy) = (
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );
    for &(x, y) in locs {
        minx = minx.min(x);
        miny = miny.min(y);
        maxx = maxx.max(x);
        maxy = maxy.max(y);
    }
    ((maxx - minx).powi(2) + (maxy - miny).powi(2)).sqrt()
}

/// Solves the small dense system `a x = b` (Gaussian elimination with partial
/// pivoting). Returns `None` if `a` is singular.
#[allow(clippy::needless_range_loop)]
fn solve(a_in: &[Vec<f64>], b_in: &[f64]) -> Option<Vec<f64>> {
    let n = b_in.len();
    let mut a: Vec<Vec<f64>> = a_in.to_vec();
    let mut b = b_in.to_vec();
    for col in 0..n {
        let mut piv = col;
        let mut best = a[col][col].abs();
        for r in col + 1..n {
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
        let d = a[col][col];
        for r in col + 1..n {
            let f = a[r][col] / d;
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
        for c in i + 1..n {
            s -= a[i][c] * x[c];
        }
        x[i] = s / a[i][i];
    }
    Some(x)
}

// ── Kernel ───────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kernel {
    Gaussian,
    Bisquare,
}

impl Kernel {
    fn as_str(self) -> &'static str {
        match self {
            Self::Gaussian => "gaussian",
            Self::Bisquare => "bisquare",
        }
    }
    fn weight(self, d: f64, b: f64) -> f64 {
        if b <= 0.0 {
            return 0.0;
        }
        let t = d / b;
        match self {
            Self::Gaussian => (-0.5 * t * t).exp(),
            Self::Bisquare => {
                if t < 1.0 {
                    let u = 1.0 - t * t;
                    u * u
                } else {
                    0.0
                }
            }
        }
    }
}

// ── Helpers / parameters ────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .map(str::trim)
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

fn rep_point(geom: &Geometry) -> Option<(f64, f64)> {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0usize;
    collect(geom, &mut |x, y| {
        sx += x;
        sy += y;
        n += 1;
    });
    (n > 0).then(|| (sx / n as f64, sy / n as f64))
}

fn collect(geom: &Geometry, f: &mut impl FnMut(f64, f64)) {
    match geom {
        Geometry::Point(c) => f(c.x, c.y),
        Geometry::MultiPoint(cs) | Geometry::LineString(cs) => cs.iter().for_each(|c| f(c.x, c.y)),
        Geometry::MultiLineString(ls) => ls.iter().flatten().for_each(|c| f(c.x, c.y)),
        Geometry::Polygon { exterior, .. } => exterior.coords().iter().for_each(|c| f(c.x, c.y)),
        Geometry::MultiPolygon(parts) => parts
            .iter()
            .flat_map(|(e, _)| e.coords())
            .for_each(|c| f(c.x, c.y)),
        Geometry::GeometryCollection(gs) => gs.iter().for_each(|g| collect(g, f)),
    }
}

struct Params {
    kernel: Kernel,
    adaptive: bool,
    tolerance: f64,
    max_iterations: usize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let kernel = match parse_optional_str(args, "kernel")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("gaussian") => Kernel::Gaussian,
        Some("bisquare") => Kernel::Bisquare,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown kernel '{other}' (expected gaussian or bisquare)"
            )))
        }
    };
    let adaptive = match parse_optional_str(args, "bandwidth_type")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("adaptive") => true,
        Some("fixed") => false,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown bandwidth_type '{other}' (expected adaptive or fixed)"
            )))
        }
    };
    let tolerance = parse_optional_f64(args, "tolerance")?.unwrap_or(1e-3);
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'tolerance' must be a positive number".to_string(),
        ));
    }
    let max_iterations = parse_optional_u64(args, "max_iterations")?.unwrap_or(30);
    if max_iterations == 0 {
        return Err(ToolError::Validation(
            "parameter 'max_iterations' must be at least 1".to_string(),
        ));
    }
    Ok(Params {
        kernel,
        adaptive,
        tolerance,
        max_iterations: max_iterations as usize,
    })
}

fn parse_optional_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
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

fn parse_optional_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
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
    use wbvector::{memory_store, FieldDef, FieldType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn layer_with(rows: &[(f64, f64, f64, f64, f64)]) -> String {
        // (x, y, x1, x2, yval)
        let mut layer = Layer::new("pts");
        layer.add_field(FieldDef::new("x1", FieldType::Float));
        layer.add_field(FieldDef::new("x2", FieldType::Float));
        layer.add_field(FieldDef::new("yv", FieldType::Float));
        for &(px, py, x1, x2, yv) in rows {
            layer
                .add_feature(
                    Some(Geometry::point(px, py)),
                    &[
                        ("x1", FieldValue::Float(x1)),
                        ("x2", FieldValue::Float(x2)),
                        ("yv", FieldValue::Float(yv)),
                    ],
                )
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        MgwrTool.run(&args, &ctx()).unwrap()
    }

    /// x1's coefficient varies smoothly across space (a location-dependent
    /// gradient in `px`); x2's coefficient is constant everywhere. MGWR
    /// should recover a *larger* (more global) bandwidth for the constant
    /// term x2 than for the spatially-varying term x1, back-fitting should
    /// converge, and MGWR's AICc should be no worse than the single-shared-
    /// bandwidth warm start's AICc on the same data.
    #[test]
    fn recovers_distinct_per_variable_bandwidths_and_beats_single_bandwidth_aicc() {
        let mut rows = Vec::new();
        let mut idx = 0i32;
        for gx in 0..12i32 {
            for gy in 0..8i32 {
                let px = gx as f64;
                let py = gy as f64;
                let fi = idx as f64;
                // x1 and x2 vary from point to point but are deliberately
                // *not* spatially smooth: high-frequency sinusoids in the
                // point index look "noise-like" from any small spatial
                // window (consecutive/nearby points get very different
                // values), so a small neighbourhood cannot find a spurious
                // *local* correlation between x1 or x2 and whatever smooth
                // spatial trend is actually attributable to a different
                // term. A lower-frequency (smooth) explanatory value would
                // let a term whose true coefficient is constant "borrow" a
                // small neighbourhood's positional trend and fake a
                // spatially varying fit, defeating the point of this test.
                // Both have a realistic nonzero mean, like real attribute
                // data (income, distance, population, ...) — the tool
                // mean-centres explanatory columns internally before
                // back-fitting for exactly this reason (a nonzero-mean
                // regressor is otherwise structurally collinear with the
                // intercept under any local window), so this doesn't need
                // to be pre-centred here.
                let x1 = 2.0 + 1.5 * (7.3 * fi + 0.3).sin();
                let x2 = 3.0 + 1.2 * (5.9 * fi + 1.1).cos();
                let beta1 = 1.0 + 0.2 * px; // smooth spatial gradient in px
                let beta2 = 4.0; // constant everywhere
                let intercept = 8.0;
                // A deterministic (non-random) perturbation so the fit is
                // never exactly noise-free: with a perfectly exact linear
                // relationship the RSS is ~0 at every bandwidth, which makes
                // AICc's ln(sigma) term numerically degenerate. Real data
                // always carries some unexplained variation.
                let noise = 0.5 * (9.7 * fi + 2.0).sin() + 0.3 * (11.3 * fi).cos();
                idx += 1;
                let y = intercept + beta1 * x1 + beta2 * x2 + noise;
                rows.push((px, py, x1, x2, y));
            }
        }
        let out = run(json!({
            "input": layer_with(&rows),
            "dependent_field": "yv",
            "explanatory_fields": "x1,x2",
            "bandwidth_type": "adaptive",
        }));

        assert!(
            out.outputs["converged"].as_bool().unwrap(),
            "back-fitting should converge on this well-behaved data: {out:?}",
        );
        assert!(out.outputs["iterations"].as_u64().unwrap() >= 1);

        let vb = out.outputs["variable_bandwidths"].as_array().unwrap();
        let bw = |field: &str| -> f64 {
            vb.iter()
                .find(|v| v["field"] == field)
                .unwrap_or_else(|| panic!("no bandwidth entry for {field}"))["bandwidth"]
                .as_f64()
                .unwrap()
        };
        let b_x1 = bw("x1");
        let b_x2 = bw("x2");
        assert!(
            b_x1 < b_x2,
            "spatially-varying x1 should get a smaller bandwidth than constant x2: b_x1={b_x1}, b_x2={b_x2}"
        );

        let mgwr_aicc = out.outputs["aicc"].as_f64().unwrap();
        let single_aicc = out.outputs["single_bandwidth_aicc"].as_f64().unwrap();
        assert!(mgwr_aicc.is_finite() && single_aicc.is_finite());
        assert!(
            mgwr_aicc <= single_aicc + 1e-6,
            "MGWR AICc ({mgwr_aicc}) should be <= single-bandwidth GWR AICc ({single_aicc})"
        );
    }

    #[test]
    fn deterministic_across_runs() {
        let mut rows = Vec::new();
        for i in 0..40 {
            let px = (i % 8) as f64;
            let py = (i / 8) as f64;
            let x1 = (i as f64 * 0.37).sin() * 3.0 + px;
            let x2 = py + (i as f64 * 0.11).cos();
            rows.push((px, py, x1, x2, 4.0 + 1.5 * x1 + 0.7 * x2));
        }
        let args = json!({
            "input": layer_with(&rows),
            "dependent_field": "yv",
            "explanatory_fields": "x1,x2",
            "max_iterations": 15,
        });
        let out1 = run(args.clone());
        let out2 = run(args);
        assert_eq!(out1.outputs["aicc"], out2.outputs["aicc"]);
        assert_eq!(
            out1.outputs["variable_bandwidths"],
            out2.outputs["variable_bandwidths"]
        );
        assert_eq!(out1.outputs["r2"], out2.outputs["r2"]);
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = MgwrTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "input": "p.geojson", "dependent_field": "yv" })).is_err(),
            "no explanatory_fields"
        );
        assert!(bad(json!({
            "input": "p.geojson", "dependent_field": "yv", "explanatory_fields": "x1", "kernel": "tri"
        }))
        .is_err());
        assert!(bad(json!({
            "input": "p.geojson", "dependent_field": "yv", "explanatory_fields": "x1", "bandwidth_type": "diag"
        }))
        .is_err());
        assert!(bad(json!({
            "input": "p.geojson", "dependent_field": "yv", "explanatory_fields": "x1", "tolerance": 0
        }))
        .is_err());
        assert!(bad(json!({
            "input": "p.geojson", "dependent_field": "yv", "explanatory_fields": "x1", "max_iterations": 0
        }))
        .is_err());
        assert!(bad(json!({
            "input": "p.geojson", "dependent_field": "yv", "explanatory_fields": "x1"
        }))
        .is_ok());
    }
}
