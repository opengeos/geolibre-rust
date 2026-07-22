//! GeoLibre tool: presence-only (MaxEnt-style) species/suitability prediction.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Presence-only Prediction (MaxEnt)*
//! (Spatial Statistics). Species-distribution / habitat-suitability modelling
//! from **presence points + explanatory rasters** is a genuine first for the
//! browser stack — the ecology community otherwise relies on the Java MaxEnt jar
//! or R's `dismo`/`maxnet`. It extends the suitability pipeline beyond the
//! un-fitted `fuzzy_overlay` / `calculate_composite_index` into a *fitted* model.
//!
//! ## Method
//! MaxEnt is, up to the link function, equivalent to a logistic regression that
//! contrasts **presence** locations against a large **background** sample of the
//! landscape (Fithian & Hastie 2013; Renner et al. 2015). We:
//!   1. sample a seeded, reproducible background set of cells over the raster
//!      stack (reservoir sampling with a splitmix64 PRNG — no `Date::now`/OS RNG,
//!      so WASM runs are deterministic);
//!   2. sample every explanatory raster at each presence + background location
//!      (bilinear);
//!   3. standardize each covariate to the background distribution, then apply a
//!      **basis expansion** (linear / quadratic / hinge feature classes);
//!   4. fit an **L1-regularized (lasso) logistic regression** of presence-vs-
//!      background by proximal-gradient descent (ISTA + soft-thresholding) — no
//!      linear-algebra crate, deterministic, and the L1 penalty performs feature
//!      selection just like MaxEnt's regularization;
//!   5. evaluate the fitted model over the raster stack to produce a
//!      **probability-of-presence surface**, plus a coefficient / variable-
//!      importance report.
//!
//! ## Link function
//! The output surface is the **logistic** probability of the fitted presence-vs-
//! background classifier (relative habitat suitability in `[0,1]`, centred near
//! 0.5 because presence and background classes are weighted equally). This is the
//! "logistic output" convention; MaxEnt's newer cloglog output is a monotone
//! re-scaling that preserves ranking (hence AUC) — see the PR for the scope note.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
    write_text_output,
};
use crate::vector_common::load_input_layer;
use wbvector::Geometry;

const NODATA: f64 = -9999.0;

pub struct PresenceOnlyPredictionTool;

impl Tool for PresenceOnlyPredictionTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "presence_only_prediction",
            display_name: "Presence-only Prediction",
            summary: "MaxEnt-style species/habitat-suitability modelling: fit an L1-regularized logistic model contrasting presence points against a seeded background sample of explanatory rasters (linear/quadratic/hinge basis expansion), then predict a probability-of-presence surface plus a variable-importance report — like ArcGIS Presence-only Prediction (MaxEnt).",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input presence point vector layer (observed occurrences).",
                    required: true,
                },
                ToolParamSpec {
                    name: "explanatory",
                    description: "Comma-separated list of explanatory raster paths (the covariate stack). Sampled at presence + background locations; the first raster defines the output grid.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output probability-of-presence raster path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "report",
                    description: "Optional output path for a CSV coefficient / variable-importance table.",
                    required: false,
                },
                ToolParamSpec {
                    name: "features",
                    description: "Comma-separated basis-expansion feature classes: any of 'linear', 'quadratic', 'hinge'. Default 'linear,quadratic'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "background",
                    description: "Number of background samples drawn from the raster stack. Default 1000 (capped at the number of valid cells).",
                    required: false,
                },
                ToolParamSpec {
                    name: "regularization",
                    description: "L1 (lasso) regularization strength on the standardized features. Larger = simpler model. Default 0.01.",
                    required: false,
                },
                ToolParamSpec {
                    name: "hinge_knots",
                    description: "Number of interior knots per covariate for hinge features (only used when 'hinge' is enabled). Default 10.",
                    required: false,
                },
                ToolParamSpec {
                    name: "seed",
                    description: "Random seed for the deterministic background sample (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_output(args, "output")?;
        let report = parse_optional_output(args, "report")?;
        let prm = parse_params(args)?;

        // ── Load inputs ──────────────────────────────────────────────────────
        let layer = load_input_layer(input)?;
        let presence_xy: Vec<(f64, f64)> = layer
            .iter()
            .filter_map(|f| f.geometry.as_ref().and_then(point_xy))
            .collect();
        if presence_xy.is_empty() {
            return Err(ToolError::Validation(
                "input layer has no point geometries".into(),
            ));
        }

        let rasters: Vec<Raster> = prm
            .explanatory
            .iter()
            .map(|p| load_input_raster(p))
            .collect::<Result<_, _>>()?;
        let ncov = rasters.len();
        let reference = &rasters[0];
        ctx.progress.info(&format!(
            "{} presence point(s), {} explanatory raster(s)",
            presence_xy.len(),
            ncov
        ));

        // ── Background sample (seeded reservoir over valid cells) ────────────
        let mut rng = prm.seed ^ 0x2545_F491_4F6C_DD1D;
        let mut reservoir: Vec<(f64, f64)> = Vec::with_capacity(prm.background);
        let mut seen: u64 = 0;
        for row in 0..reference.rows as isize {
            for col in 0..reference.cols as isize {
                let x = reference.col_center_x(col);
                let y = reference.row_center_y(row);
                if sample_stack(&rasters, x, y).is_none() {
                    continue; // some covariate has no data here
                }
                seen += 1;
                if reservoir.len() < prm.background {
                    reservoir.push((x, y));
                } else {
                    let j = next_u64(&mut rng) % seen;
                    if (j as usize) < prm.background {
                        reservoir[j as usize] = (x, y);
                    }
                }
            }
        }
        if reservoir.len() < 2 {
            return Err(ToolError::Execution(
                "fewer than 2 valid background cells in the raster stack".into(),
            ));
        }

        // ── Raw covariate matrices at presence + background ─────────────────
        let mut pres_raw: Vec<Vec<f64>> = Vec::new();
        let mut skipped = 0usize;
        for &(x, y) in &presence_xy {
            match sample_stack(&rasters, x, y) {
                Some(v) => pres_raw.push(v),
                None => skipped += 1,
            }
        }
        if pres_raw.len() < 2 {
            return Err(ToolError::Execution(
                "fewer than 2 presence points fall on valid data in every explanatory raster"
                    .into(),
            ));
        }
        let bg_raw: Vec<Vec<f64>> = reservoir
            .iter()
            .map(|&(x, y)| sample_stack(&rasters, x, y).expect("reservoir cells are valid"))
            .collect();

        // ── Standardize covariates to the background distribution ───────────
        let (cov_mean, cov_std) = standardize_params(&bg_raw, ncov);

        // ── Build the basis-expansion feature specs ─────────────────────────
        let specs = build_specs(&bg_raw, &cov_mean, &cov_std, &prm);
        let nfeat = specs.len();

        // ── Design matrix (presence rows then background rows) ──────────────
        let mut x_rows: Vec<Vec<f64>> = Vec::with_capacity(pres_raw.len() + bg_raw.len());
        for raw in pres_raw.iter().chain(bg_raw.iter()) {
            x_rows.push(expand(raw, &cov_mean, &cov_std, &specs));
        }
        let n_p = pres_raw.len();
        let n_b = bg_raw.len();

        // Column-standardize the features so the L1 penalty is scale-fair and
        // the coefficients are directly comparable for variable importance.
        let (feat_mean, feat_std) = standardize_params(&x_rows, nfeat);
        for row in &mut x_rows {
            for j in 0..nfeat {
                row[j] = (row[j] - feat_mean[j]) / feat_std[j];
            }
        }

        // Labels + class-balanced weights (presence mass == background mass).
        let mut y = vec![1.0f64; n_p];
        y.extend(std::iter::repeat_n(0.0, n_b));
        let mut w = vec![0.5 / n_p as f64; n_p];
        w.extend(std::iter::repeat_n(0.5 / n_b as f64, n_b));

        // ── Fit L1-regularized logistic regression (ISTA) ───────────────────
        let fit = fit_lasso_logistic(&x_rows, &y, &w, prm.lambda);
        ctx.progress.info(&format!(
            "fitted {nfeat} feature(s) in {} iteration(s) (converged={})",
            fit.iterations, fit.converged
        ));

        // ── Training-set diagnostics ────────────────────────────────────────
        let probs: Vec<f64> = x_rows.iter().map(|r| sigmoid(fit.eta(r))).collect();
        let mean_p_pres = probs[..n_p].iter().sum::<f64>() / n_p as f64;
        let mean_p_bg = probs[n_p..].iter().sum::<f64>() / n_b as f64;
        let train_auc = auc(&probs[..n_p], &probs[n_p..]);

        // ── Variable importance (sum of |standardized coefficient| per cov) ──
        let mut importance = vec![0.0f64; ncov];
        for (j, s) in specs.iter().enumerate() {
            importance[s.cov] += fit.beta[j].abs();
        }
        let imp_total: f64 = importance.iter().sum::<f64>().max(f64::EPSILON);
        let imp_pct: Vec<f64> = importance.iter().map(|v| 100.0 * v / imp_total).collect();

        // ── Predict the probability surface over the reference grid ─────────
        let rows = reference.rows;
        let cols = reference.cols;
        let mut surface = vec![NODATA; rows * cols];
        for row in 0..rows {
            for col in 0..cols {
                let x = reference.col_center_x(col as isize);
                let y = reference.row_center_y(row as isize);
                if let Some(raw) = sample_stack(&rasters, x, y) {
                    let feat = expand(&raw, &cov_mean, &cov_std, &specs);
                    let mut std_feat = vec![0.0f64; nfeat];
                    for j in 0..nfeat {
                        std_feat[j] = (feat[j] - feat_mean[j]) / feat_std[j];
                    }
                    surface[row * cols + col] = sigmoid(fit.eta(&std_feat));
                }
            }
        }
        let out_raster = raster_like_with_data(reference, surface, NODATA, DataType::F32)?;
        let out_path = write_or_store_output(out_raster, output)?;

        // ── Report ──────────────────────────────────────────────────────────
        let cov_names: Vec<String> = prm.explanatory.iter().map(|p| basename(p)).collect();
        let report_path = match report {
            Some(path) => {
                write_text_output(&build_report(&fit, &specs, &cov_names, &imp_pct), path)?;
                Some(path.to_string())
            }
            None => None,
        };

        // ── Result ──────────────────────────────────────────────────────────
        let mut imp_map = serde_json::Map::new();
        for (name, pct) in cov_names.iter().zip(imp_pct.iter()) {
            imp_map.insert(name.clone(), json!((pct * 100.0).round() / 100.0));
        }

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        if let Some(p) = report_path {
            outputs.insert("report".to_string(), json!(p));
        }
        outputs.insert("n_presence".to_string(), json!(n_p));
        outputs.insert("n_presence_skipped".to_string(), json!(skipped));
        outputs.insert("n_background".to_string(), json!(n_b));
        outputs.insert("n_covariates".to_string(), json!(ncov));
        outputs.insert("n_features".to_string(), json!(nfeat));
        outputs.insert("iterations".to_string(), json!(fit.iterations));
        outputs.insert("converged".to_string(), json!(fit.converged));
        outputs.insert(
            "mean_prob_presence".to_string(),
            json!((mean_p_pres * 1e6).round() / 1e6),
        );
        outputs.insert(
            "mean_prob_background".to_string(),
            json!((mean_p_bg * 1e6).round() / 1e6),
        );
        outputs.insert(
            "training_auc".to_string(),
            json!((train_auc * 1e6).round() / 1e6),
        );
        outputs.insert("importance".to_string(), Value::Object(imp_map));
        Ok(ToolRunResult { outputs })
    }
}

// ── Raster stack sampling ──────────────────────────────────────────────────────

/// Samples every raster in the stack at world coordinates `(x, y)` (bilinear).
/// Returns `None` if any raster has no data there.
fn sample_stack(rasters: &[Raster], x: f64, y: f64) -> Option<Vec<f64>> {
    let mut out = Vec::with_capacity(rasters.len());
    for r in rasters {
        out.push(sample_bilinear(r, x, y)?);
    }
    Some(out)
}

fn cell_value(r: &Raster, row: isize, col: isize) -> Option<f64> {
    if row < 0 || col < 0 || row >= r.rows as isize || col >= r.cols as isize {
        return None;
    }
    let v = r.get(0, row, col);
    if v == r.nodata || v.is_nan() {
        None
    } else {
        Some(v)
    }
}

fn sample_bilinear(r: &Raster, x: f64, y: f64) -> Option<f64> {
    let fx = (x - r.x_min) / r.cell_size_x - 0.5;
    let fy = (r.y_max() - y) / r.cell_size_y - 0.5;
    let col0 = fx.floor() as isize;
    let row0 = fy.floor() as isize;
    let tx = fx - col0 as f64;
    let ty = fy - row0 as f64;

    let v00 = cell_value(r, row0, col0);
    let v01 = cell_value(r, row0, col0 + 1);
    let v10 = cell_value(r, row0 + 1, col0);
    let v11 = cell_value(r, row0 + 1, col0 + 1);

    if let (Some(a), Some(b), Some(c), Some(d)) = (v00, v01, v10, v11) {
        let top = a * (1.0 - tx) + b * tx;
        let bot = c * (1.0 - tx) + d * tx;
        return Some(top * (1.0 - ty) + bot * ty);
    }
    // Fall back to the nearest valid corner by bilinear weight.
    [
        (v00, (1.0 - tx) * (1.0 - ty)),
        (v01, tx * (1.0 - ty)),
        (v10, (1.0 - tx) * ty),
        (v11, tx * ty),
    ]
    .into_iter()
    .filter_map(|(v, wt)| v.map(|v| (v, wt)))
    .max_by(|a, b| a.1.total_cmp(&b.1))
    .map(|(v, _)| v)
}

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        _ => None,
    }
}

// ── Standardization & basis expansion ──────────────────────────────────────────

/// Column-wise mean and (population) standard deviation over a set of rows.
/// A zero/degenerate std is replaced by 1 so the column becomes constant-0.
fn standardize_params(rows: &[Vec<f64>], ncol: usize) -> (Vec<f64>, Vec<f64>) {
    let n = rows.len().max(1) as f64;
    let mut mean = vec![0.0f64; ncol];
    for r in rows {
        for j in 0..ncol {
            mean[j] += r[j];
        }
    }
    for m in &mut mean {
        *m /= n;
    }
    let mut var = vec![0.0f64; ncol];
    for r in rows {
        for j in 0..ncol {
            let d = r[j] - mean[j];
            var[j] += d * d;
        }
    }
    let std: Vec<f64> = var
        .iter()
        .map(|v| {
            let s = (v / n).sqrt();
            if s.is_finite() && s > 1e-12 {
                s
            } else {
                1.0
            }
        })
        .collect();
    (mean, std)
}

#[derive(Clone)]
enum Kind {
    Linear,
    Quadratic,
    HingeFwd { knot: f64, hi: f64 },
    HingeRev { knot: f64, lo: f64 },
}

#[derive(Clone)]
struct FeatSpec {
    cov: usize,
    kind: Kind,
    label: String,
}

/// Standardized covariate value `z = (raw - mean) / std`.
#[inline]
fn zval(raw: f64, mean: f64, std: f64) -> f64 {
    (raw - mean) / std
}

/// Builds the feature specs for every covariate given the requested classes.
/// Hinge knots are placed at interior fractions of the background's standardized
/// covariate range.
fn build_specs(
    bg_raw: &[Vec<f64>],
    cov_mean: &[f64],
    cov_std: &[f64],
    prm: &Params,
) -> Vec<FeatSpec> {
    let ncov = cov_mean.len();
    let mut specs = Vec::new();
    for c in 0..ncov {
        if prm.features.linear {
            specs.push(FeatSpec {
                cov: c,
                kind: Kind::Linear,
                label: "linear".into(),
            });
        }
        if prm.features.quadratic {
            specs.push(FeatSpec {
                cov: c,
                kind: Kind::Quadratic,
                label: "quadratic".into(),
            });
        }
        if prm.features.hinge {
            // Range of standardized covariate c over the background.
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for r in bg_raw {
                let z = zval(r[c], cov_mean[c], cov_std[c]);
                lo = lo.min(z);
                hi = hi.max(z);
            }
            if hi > lo {
                let k = prm.hinge_knots.max(1);
                for i in 1..=k {
                    let knot = lo + (hi - lo) * (i as f64) / (k as f64 + 1.0);
                    specs.push(FeatSpec {
                        cov: c,
                        kind: Kind::HingeFwd { knot, hi },
                        label: format!("hinge_fwd_{i}"),
                    });
                    specs.push(FeatSpec {
                        cov: c,
                        kind: Kind::HingeRev { knot, lo },
                        label: format!("hinge_rev_{i}"),
                    });
                }
            }
        }
    }
    specs
}

/// Expands a raw covariate vector into the (un-standardized) feature vector.
fn expand(raw: &[f64], cov_mean: &[f64], cov_std: &[f64], specs: &[FeatSpec]) -> Vec<f64> {
    specs
        .iter()
        .map(|s| {
            let z = zval(raw[s.cov], cov_mean[s.cov], cov_std[s.cov]);
            match s.kind {
                Kind::Linear => z,
                Kind::Quadratic => z * z,
                Kind::HingeFwd { knot, hi } => ((z - knot) / (hi - knot)).clamp(0.0, 1.0),
                Kind::HingeRev { knot, lo } => ((knot - z) / (knot - lo)).clamp(0.0, 1.0),
            }
        })
        .collect()
}

// ── L1-regularized logistic regression (proximal gradient / ISTA) ──────────────

struct Fit {
    intercept: f64,
    beta: Vec<f64>,
    iterations: usize,
    converged: bool,
}

impl Fit {
    #[inline]
    fn eta(&self, x: &[f64]) -> f64 {
        self.intercept + x.iter().zip(&self.beta).map(|(a, b)| a * b).sum::<f64>()
    }
}

#[inline]
fn sigmoid(z: f64) -> f64 {
    1.0 / (1.0 + (-z.clamp(-40.0, 40.0)).exp())
}

#[inline]
fn soft_threshold(z: f64, t: f64) -> f64 {
    if z > t {
        z - t
    } else if z < -t {
        z + t
    } else {
        0.0
    }
}

/// Fits `min_β Σ w_i loss(y_i, β·x_i) + λ ||β||_1` (intercept unpenalized) by
/// proximal-gradient descent with a fixed Lipschitz step. Deterministic.
fn fit_lasso_logistic(x: &[Vec<f64>], y: &[f64], w: &[f64], lambda: f64) -> Fit {
    let n = x.len();
    let p = if n > 0 { x[0].len() } else { 0 };
    let mut intercept = 0.0f64;
    let mut beta = vec![0.0f64; p];

    // Lipschitz constant of the smooth-part gradient: L = 0.25 * Σ w_i (1+‖x_i‖²).
    let mut lip = 0.0f64;
    for (row, &wi) in x.iter().zip(w) {
        lip += wi * (1.0 + row.iter().map(|v| v * v).sum::<f64>());
    }
    lip *= 0.25;
    let step = if lip > 0.0 { 1.0 / lip } else { 1.0 };

    let max_iter = 5000usize;
    let tol = 1e-7;
    let mut iterations = 0;
    let mut converged = false;
    for _ in 0..max_iter {
        iterations += 1;
        // Gradient of the weighted logistic loss.
        let mut g0 = 0.0f64;
        let mut g = vec![0.0f64; p];
        for i in 0..n {
            let eta = intercept + x[i].iter().zip(&beta).map(|(a, b)| a * b).sum::<f64>();
            let r = w[i] * (sigmoid(eta) - y[i]);
            g0 += r;
            for j in 0..p {
                g[j] += r * x[i][j];
            }
        }
        intercept -= step * g0;
        let mut max_delta = step * g0.abs();
        for j in 0..p {
            let old = beta[j];
            beta[j] = soft_threshold(beta[j] - step * g[j], step * lambda);
            max_delta = max_delta.max((beta[j] - old).abs());
        }
        if max_delta < tol {
            converged = true;
            break;
        }
    }

    Fit {
        intercept,
        beta,
        iterations,
        converged,
    }
}

/// Mann-Whitney AUC of positive vs negative scores (0.5 for ties).
fn auc(pos: &[f64], neg: &[f64]) -> f64 {
    if pos.is_empty() || neg.is_empty() {
        return 0.5;
    }
    let mut wins = 0.0f64;
    for &p in pos {
        for &q in neg {
            if p > q {
                wins += 1.0;
            } else if p == q {
                wins += 0.5;
            }
        }
    }
    wins / (pos.len() as f64 * neg.len() as f64)
}

// ── Report ─────────────────────────────────────────────────────────────────────

fn build_report(fit: &Fit, specs: &[FeatSpec], cov_names: &[String], imp_pct: &[f64]) -> String {
    let mut s = String::from("covariate,feature,coefficient\n");
    s.push_str(&format!("intercept,intercept,{:.6}\n", fit.intercept));
    for (j, spec) in specs.iter().enumerate() {
        s.push_str(&format!(
            "{},{},{:.6}\n",
            cov_names[spec.cov], spec.label, fit.beta[j]
        ));
    }
    s.push_str("\n# variable importance (sum |standardized coefficient|, percent)\n");
    s.push_str("covariate,importance_percent\n");
    for (name, pct) in cov_names.iter().zip(imp_pct.iter()) {
        s.push_str(&format!("{name},{pct:.2}\n"));
    }
    s
}

fn basename(path: &str) -> String {
    std::path::Path::new(path)
        .file_stem()
        .and_then(|s| s.to_str())
        .map(|s| s.to_string())
        .unwrap_or_else(|| path.to_string())
}

// ── Splitmix64 PRNG (matches colocation_analysis / spatially_balanced) ─────────

fn next_u64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

// ── Parameters ──────────────────────────────────────────────────────────────────

struct FeatureSet {
    linear: bool,
    quadratic: bool,
    hinge: bool,
}

struct Params {
    explanatory: Vec<String>,
    features: FeatureSet,
    background: usize,
    lambda: f64,
    hinge_knots: usize,
    seed: u64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let raw = require_str(args, "explanatory")?;
    let explanatory: Vec<String> = raw
        .split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();
    if explanatory.is_empty() {
        return Err(ToolError::Validation(
            "'explanatory' must list at least one raster path".into(),
        ));
    }

    let features = parse_features(args)?;
    let background = opt_usize(args, "background")?.unwrap_or(1000);
    if background < 2 {
        return Err(ToolError::Validation(
            "'background' must be at least 2".into(),
        ));
    }
    let lambda = opt_f64(args, "regularization")?.unwrap_or(0.01);
    if lambda < 0.0 || !lambda.is_finite() {
        return Err(ToolError::Validation(
            "'regularization' must be a non-negative number".into(),
        ));
    }
    let hinge_knots = opt_usize(args, "hinge_knots")?.unwrap_or(10).max(1);
    let seed = opt_usize(args, "seed")?.unwrap_or(1) as u64;

    Ok(Params {
        explanatory,
        features,
        background,
        lambda,
        hinge_knots,
        seed,
    })
}

fn parse_features(args: &ToolArgs) -> Result<FeatureSet, ToolError> {
    let spec = match args.get("features") {
        None | Some(Value::Null) => "linear,quadratic".to_string(),
        Some(Value::String(s)) if s.trim().is_empty() => "linear,quadratic".to_string(),
        Some(Value::String(s)) => s.clone(),
        Some(_) => {
            return Err(ToolError::Validation(
                "'features' must be a comma-separated string".into(),
            ))
        }
    };
    let mut set = FeatureSet {
        linear: false,
        quadratic: false,
        hinge: false,
    };
    for tok in spec.split(',').map(|t| t.trim().to_lowercase()) {
        match tok.as_str() {
            "" => {}
            "linear" => set.linear = true,
            "quadratic" => set.quadratic = true,
            "hinge" => set.hinge = true,
            other => {
                return Err(ToolError::Validation(format!(
                    "unknown feature class '{other}' (expected linear, quadratic, or hinge)"
                )))
            }
        }
    }
    if !(set.linear || set.quadratic || set.hinge) {
        return Err(ToolError::Validation(
            "'features' must enable at least one of linear, quadratic, hinge".into(),
        ));
    }
    Ok(set)
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

fn opt_usize(args: &ToolArgs, key: &str) -> Result<Option<usize>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => n.as_u64().map(|v| Some(v as usize)).ok_or_else(|| {
            ToolError::Validation(format!("parameter '{key}' must be a whole number"))
        }),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s.trim().parse::<usize>().map(Some).map_err(|_| {
            ToolError::Validation(format!("parameter '{key}' must be a whole number"))
        }),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a whole number"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{memory_store, RasterConfig};
    use wbvector::{FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// A single-band raster from a row-major closure `f(row,col) -> value`.
    fn make_raster(rows: usize, cols: usize, f: impl Fn(usize, usize) -> f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: Default::default(),
            metadata: Default::default(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, f(row, col)).unwrap();
            }
        }
        let id = memory_store::put_raster(r);
        memory_store::make_raster_memory_path(&id)
    }

    /// Presence points at the given world coordinates.
    fn make_points(pts: &[(f64, f64)]) -> String {
        let mut l = Layer::new("pres").with_geom_type(GeometryType::Point);
        l.add_field(FieldDef::new("id", FieldType::Integer));
        for (i, &(x, y)) in pts.iter().enumerate() {
            l.add_feature(Some(Geometry::point(x, y)), &[("id", (i as i64).into())])
                .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        PresenceOnlyPredictionTool.run(&args, &ctx()).unwrap()
    }

    /// A grid whose "elevation" rises to the north (row 0). Presence points are
    /// planted in the high-elevation band, so the model must recover the signal:
    /// AUC ≫ 0.5 and mean presence probability ≫ background.
    fn planted() -> (String, Vec<(f64, f64)>) {
        let rows = 40usize;
        let cols = 40usize;
        // elevation = (rows-row): high in the north, low in the south.
        let elev = make_raster(rows, cols, |row, _| (rows - row) as f64);
        // Presence in the top 6 rows (high elevation), spread across columns.
        let mut pts = Vec::new();
        for row in 0..6 {
            for col in (0..cols).step_by(2) {
                pts.push((col as f64 + 0.5, (rows - 1 - row) as f64 + 0.5));
            }
        }
        (elev, pts)
    }

    #[test]
    fn recovers_planted_elevation_signal() {
        let (elev, pts) = planted();
        let out = run(json!({
            "input": make_points(&pts),
            "explanatory": elev,
            "features": "linear,quadratic",
            "background": 400,
            "seed": 7,
        }));
        let auc = out.outputs["training_auc"].as_f64().unwrap();
        let mp = out.outputs["mean_prob_presence"].as_f64().unwrap();
        let mb = out.outputs["mean_prob_background"].as_f64().unwrap();
        assert!(auc > 0.85, "planted signal AUC should be high, got {auc}");
        assert!(
            mp > mb + 0.2,
            "presence prob {mp} should far exceed background {mb}"
        );
    }

    #[test]
    fn deterministic_same_seed() {
        let (elev, pts) = planted();
        let mk = |seed: u64| {
            let out = run(json!({
                "input": make_points(&pts), "explanatory": elev.clone(),
                "background": 300, "seed": seed,
            }));
            let path = out.outputs["output"].as_str().unwrap().to_string();
            let r = load_input_raster(&path).unwrap();
            crate::common::band_to_vec(&r, 0)
        };
        assert_eq!(mk(11), mk(11), "same seed must give identical surface");
        assert_ne!(mk(11), mk(12), "different seed should differ");
    }

    #[test]
    fn noise_covariate_gets_low_importance() {
        // Covariate 0 = elevation signal; covariate 1 = deterministic checkerboard
        // noise unrelated to presence. Importance should concentrate on the signal.
        let (elev, pts) = planted();
        let noise = make_raster(40, 40, |row, col| ((row * 7 + col * 13) % 11) as f64);
        let out = run(json!({
            "input": make_points(&pts),
            "explanatory": format!("{elev},{noise}"),
            "features": "linear,quadratic",
            "background": 400,
            "seed": 3,
        }));
        let imp = out.outputs["importance"].as_object().unwrap();
        let vals: Vec<f64> = imp.values().map(|v| v.as_f64().unwrap()).collect();
        let hi = vals.iter().cloned().fold(f64::MIN, f64::max);
        let lo = vals.iter().cloned().fold(f64::MAX, f64::min);
        assert!(
            hi > 70.0 && lo < 30.0,
            "signal should dominate importance, got {vals:?}"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            PresenceOnlyPredictionTool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input+explanatory");
        assert!(
            bad(json!({ "input": "a.geojson" })).is_err(),
            "missing explanatory"
        );
        assert!(
            bad(json!({ "input": "a.geojson", "explanatory": "r.tif", "features": "bogus" }))
                .is_err(),
            "bad feature class"
        );
        assert!(
            bad(json!({ "input": "a.geojson", "explanatory": "r.tif", "background": 1 })).is_err(),
            "background too small"
        );
        assert!(
            bad(json!({ "input": "a.geojson", "explanatory": "r.tif" })).is_ok(),
            "minimal valid args"
        );
    }
}
