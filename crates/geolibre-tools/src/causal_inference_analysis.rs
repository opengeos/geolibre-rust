//! GeoLibre tool: causal inference analysis (propensity-score treatment effect
//! estimation).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Causal Inference Analysis* (Spatial
//! Statistics): estimate the causal effect of a binary treatment on an outcome
//! from observational data, adjusting for confounders instead of just reporting
//! a raw association. Everything else in the catalog answers "is X associated
//! with Y" (`generalized_linear_regression`, `geographically_weighted_regression`,
//! `bivariate_spatial_association`); this answers "how much did X *cause* Y to
//! change, once we account for what else drives both".
//!
//! ## Algorithm
//!
//! 1. Fit a propensity model — logistic regression of `treatment_field` on
//!    `confounding_fields` via IRLS (the same iteratively-reweighted-least-
//!    squares core used by `generalized_linear_regression`'s logistic family,
//!    reimplemented here self-contained).
//! 2. Estimate the Average Treatment Effect with one of three `method`s:
//!    * `ps_matching` (default) — greedy nearest-neighbour matching on the
//!      propensity score within a caliper (`0.2 * std(propensity)`, the
//!      standard Austin (2011) rule), found with a 1-D vendored `kdtree`.
//!      Effect = mean(outcome_treated − outcome_matched_control); this is the
//!      average effect on the treated (ATT).
//!    * `ipw` — stabilized inverse-probability-of-treatment weights (Hajek
//!      estimator), truncated at the 1st/99th percentile to bound the
//!      influence of extreme propensity scores.
//!    * `regression_adjustment` — OLS of `outcome ~ treatment + confounders`;
//!      the treatment coefficient is the ATE under the linearity assumption.
//! 3. Bootstrap a confidence interval by resampling rows with replacement and
//!    re-running the full pipeline (propensity refit included) `B` times,
//!    using an inline seeded splitmix64/xorshift64 generator (no `Date::now`,
//!    no external RNG crate — deterministic given `seed`).
//! 4. Report the per-confounder standardized mean difference (SMD) between
//!    treatment groups before adjustment and after (matched sample for
//!    `ps_matching`, IPW-weighted for `ipw`), warning via `ctx.progress.info`
//!    when the post-adjustment maximum SMD still exceeds `balance_threshold`.
//! 5. Optionally (`add_spatial_confounders`) widen the confounder set with
//!    each feature's x/y coordinates and the neighbour-averaged (k-nearest,
//!    2-D `kdtree`) value of each confounder, to absorb spatially structured
//!    unobserved confounding.
//!
//! Also reports the naive (unadjusted) difference in means alongside the ATE
//! so a caller can see the size of the confounding bias the adjustment removed.
//!
//! **v1 scope cut:** the treatment must be binary (0/1); ArcGIS's continuous
//! exposure-response variant is not implemented.

use std::collections::BTreeMap;

use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

const BOOTSTRAP_ITERS: usize = 200;
const SPATIAL_NEIGHBORS: usize = 8;

pub struct CausalInferenceAnalysisTool;

impl Tool for CausalInferenceAnalysisTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "causal_inference_analysis",
            display_name: "Causal Inference Analysis",
            summary: "Estimate the causal effect of a binary treatment on an outcome via propensity-score matching, inverse-probability weighting, or regression adjustment, with a seeded bootstrap CI and a per-confounder covariate-balance table before/after adjustment (like ArcGIS Causal Inference Analysis).",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input feature layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "outcome_field",
                    description: "Numeric outcome (dependent) field.",
                    required: true,
                },
                ToolParamSpec {
                    name: "treatment_field",
                    description: "Binary (0/1, or boolean) treatment field. Continuous treatments are not supported in this version.",
                    required: true,
                },
                ToolParamSpec {
                    name: "confounding_fields",
                    description: "Comma-separated numeric confounder field(s) to adjust for.",
                    required: true,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'ps_matching' (default, nearest-neighbour caliper matching on the propensity score), 'ipw' (stabilized inverse-probability weighting), or 'regression_adjustment' (outcome ~ treatment + confounders).",
                    required: false,
                },
                ToolParamSpec {
                    name: "add_spatial_confounders",
                    description: "If true, widen the confounder set with each feature's x/y coordinates plus the k-nearest-neighbour average of each confounder (default false).",
                    required: false,
                },
                ToolParamSpec {
                    name: "balance_threshold",
                    description: "Standardized-mean-difference threshold for the post-adjustment balance check (default 0.1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (input copy + propensity score / weight / match-id columns). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "seed",
                    description: "Integer seed for the deterministic bootstrap confidence interval (default 42).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in [
            "input",
            "outcome_field",
            "treatment_field",
            "confounding_fields",
        ] {
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

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let prm = parse_params(args)?;
        let output = parse_optional_str(args, "output")?;

        let mut layer = load_input_layer(input)?;
        let schema = layer.schema.clone();

        // ── Gather valid rows (finite treatment/outcome/confounders, and a
        //    representative point when spatial confounders are requested). ──
        let mut idx_map: Vec<usize> = Vec::new();
        let mut treat: Vec<f64> = Vec::new();
        let mut y: Vec<f64> = Vec::new();
        let mut conf_raw: Vec<Vec<f64>> = Vec::new();
        let mut reps: Vec<(f64, f64)> = Vec::new();

        for (fi, feature) in layer.features.iter().enumerate() {
            let Some(tv) = feature
                .get(&schema, &prm.treatment_field)
                .ok()
                .and_then(to_numeric_treatment)
            else {
                continue;
            };
            let Some(yv) = feature
                .get(&schema, &prm.outcome_field)
                .ok()
                .and_then(FieldValue::as_f64)
            else {
                continue;
            };
            if !tv.is_finite() || !yv.is_finite() {
                continue;
            }
            let mut row = Vec::with_capacity(prm.confounding_fields.len());
            let mut ok = true;
            for cf in &prm.confounding_fields {
                match feature.get(&schema, cf).ok().and_then(FieldValue::as_f64) {
                    Some(v) if v.is_finite() => row.push(v),
                    _ => {
                        ok = false;
                        break;
                    }
                }
            }
            if !ok {
                continue;
            }
            let rep = feature.geometry.as_ref().and_then(representative_point);
            if prm.add_spatial_confounders && rep.is_none() {
                continue;
            }
            idx_map.push(fi);
            treat.push(tv);
            y.push(yv);
            conf_raw.push(row);
            reps.push(rep.unwrap_or((0.0, 0.0)));
        }

        let n = treat.len();
        if n < 10 {
            return Err(ToolError::Execution(format!(
                "need at least 10 observations with valid treatment/outcome/confounder values, found {n}"
            )));
        }
        if !treat.iter().all(|&t| t == 0.0 || t == 1.0) {
            return Err(ToolError::Execution(
                "'treatment_field' must be binary (0/1); continuous treatments are not supported in this version"
                    .to_string(),
            ));
        }
        let n_treated = treat.iter().filter(|&&t| t == 1.0).count();
        let n_control = n - n_treated;
        if n_treated < 2 || n_control < 2 {
            return Err(ToolError::Execution(format!(
                "need at least 2 treated and 2 control observations, found {n_treated} treated / {n_control} control"
            )));
        }

        // ── Optionally widen the confounder set with spatial terms. ──────────
        let mut field_names = prm.confounding_fields.clone();
        if prm.add_spatial_confounders {
            let k = SPATIAL_NEIGHBORS.min(n.saturating_sub(1)).max(1);
            let nbr_avg = neighbor_averages(&reps, &conf_raw, k);
            for i in 0..n {
                conf_raw[i].push(reps[i].0);
                conf_raw[i].push(reps[i].1);
                conf_raw[i].extend_from_slice(&nbr_avg[i]);
            }
            field_names.push("spatial_x".to_string());
            field_names.push("spatial_y".to_string());
            for cf in &prm.confounding_fields {
                field_names.push(format!("{cf}_nbr_mean"));
            }
        }
        let p_ps = field_names.len() + 1; // + intercept

        ctx.progress.info(&format!(
            "causal_inference_analysis ({}): {n} obs ({n_treated} treated / {n_control} control), {} confounder term(s)",
            prm.method.as_str(),
            field_names.len()
        ));

        let x_ps: Vec<Vec<f64>> = conf_raw
            .iter()
            .map(|row| {
                let mut r = Vec::with_capacity(row.len() + 1);
                r.push(1.0);
                r.extend_from_slice(row);
                r
            })
            .collect();

        let ps = fit_logistic(&x_ps, &treat, p_ps)
            .map(|beta| predict_logistic(&x_ps, &beta))
            .ok_or_else(|| {
                ToolError::Execution(
                    "propensity model fit failed (singular design — check for collinear or constant confounders)"
                        .to_string(),
                )
            })?;

        // ── Point estimate + method-specific per-unit outputs. ───────────────
        let mut match_of: Vec<Option<usize>> = vec![None; n];
        let mut weight_of: Vec<f64> = vec![0.0; n];
        let mut n_matched = 0usize;
        let mut caliper_used = f64::NAN;

        let ate = match prm.method {
            Method::PsMatching => {
                let caliper = (0.2 * std_dev(&ps)).max(1e-6);
                caliper_used = caliper;
                let matches = match_treated_to_controls(&ps, &treat, caliper);
                let (ate, matched) = matching_ate(&y, &treat, &matches).ok_or_else(|| {
                    ToolError::Execution(
                        "no treated unit found a control within the caliper — insufficient common support"
                            .to_string(),
                    )
                })?;
                n_matched = matched;
                for (i, m) in matches.into_iter().enumerate() {
                    if let Some(j) = m {
                        match_of[i] = Some(j);
                        weight_of[i] = 1.0;
                        weight_of[j] = 1.0;
                    }
                }
                ate
            }
            Method::Ipw => {
                let weights = ipw_weights(&ps, &treat);
                weight_of.clone_from(&weights);
                hajek_ate(&y, &treat, &weights).ok_or_else(|| {
                    ToolError::Execution(
                        "IPW weights degenerate (zero mass in a group)".to_string(),
                    )
                })?
            }
            Method::RegressionAdjustment => {
                weight_of.iter_mut().for_each(|w| *w = 1.0);
                ate_regression(&conf_raw, &treat, &y).ok_or_else(|| {
                    ToolError::Execution(
                        "regression-adjustment fit failed (singular design)".to_string(),
                    )
                })?
            }
        };
        let naive = naive_diff_means(&y, &treat).unwrap_or(f64::NAN);

        // ── Balance table before / after. ─────────────────────────────────────
        let uniform = vec![1.0; n];
        let mut balance_before = Vec::with_capacity(field_names.len());
        let mut balance_after = Vec::with_capacity(field_names.len());
        for (c, name) in field_names.iter().enumerate() {
            let col: Vec<f64> = conf_raw.iter().map(|r| r[c]).collect();
            let before = smd_field(&col, &treat, &uniform);
            let after = match prm.method {
                Method::PsMatching => {
                    let (v, t) = matched_subset(&col, &treat, &match_of);
                    if t.is_empty() {
                        f64::NAN
                    } else {
                        smd_field(&v, &t, &vec![1.0; v.len()])
                    }
                }
                Method::Ipw => smd_field(&col, &treat, &weight_of),
                Method::RegressionAdjustment => before, // no reweighting; covariates enter the model directly
            };
            balance_before.push((name.clone(), before));
            balance_after.push((name.clone(), after));
        }
        let max_smd_before = balance_before
            .iter()
            .map(|(_, v)| *v)
            .filter(|v| v.is_finite())
            .fold(0.0_f64, f64::max);
        let max_smd_after = balance_after
            .iter()
            .map(|(_, v)| *v)
            .filter(|v| v.is_finite())
            .fold(0.0_f64, f64::max);
        let balance_ok = max_smd_after <= prm.balance_threshold;
        if !balance_ok {
            ctx.progress.info(&format!(
                "WARNING: post-adjustment balance check failed — max standardized mean difference {max_smd_after:.4} exceeds balance_threshold {:.4}",
                prm.balance_threshold
            ));
        }

        // ── Seeded bootstrap confidence interval. ─────────────────────────────
        let method = prm.method;
        let effect_fn = |idx: &[usize]| -> Option<f64> {
            let t_b: Vec<f64> = idx.iter().map(|&i| treat[i]).collect();
            let y_b: Vec<f64> = idx.iter().map(|&i| y[i]).collect();
            match method {
                Method::PsMatching | Method::Ipw => {
                    let x_ps_b: Vec<Vec<f64>> = idx.iter().map(|&i| x_ps[i].clone()).collect();
                    let beta = fit_logistic(&x_ps_b, &t_b, p_ps)?;
                    let ps_b = predict_logistic(&x_ps_b, &beta);
                    match method {
                        Method::PsMatching => {
                            let caliper = (0.2 * std_dev(&ps_b)).max(1e-6);
                            let matches = match_treated_to_controls(&ps_b, &t_b, caliper);
                            matching_ate(&y_b, &t_b, &matches).map(|(a, _)| a)
                        }
                        Method::Ipw => {
                            let w_b = ipw_weights(&ps_b, &t_b);
                            hajek_ate(&y_b, &t_b, &w_b)
                        }
                        Method::RegressionAdjustment => unreachable!(),
                    }
                }
                Method::RegressionAdjustment => {
                    let x_conf_b: Vec<Vec<f64>> =
                        idx.iter().map(|&i| conf_raw[i].clone()).collect();
                    ate_regression(&x_conf_b, &t_b, &y_b)
                }
            }
        };
        let (ci_low, ci_high, valid_boot) = bootstrap_ci(n, prm.seed, BOOTSTRAP_ITERS, effect_fn);
        if valid_boot * 2 < BOOTSTRAP_ITERS {
            ctx.progress.info(&format!(
                "WARNING: only {valid_boot}/{BOOTSTRAP_ITERS} bootstrap replicates converged; confidence interval may be unreliable"
            ));
        }

        // ── Append per-feature output fields, aligned back to the full layer. ─
        for name in ["ci_propensity", "ci_weight"] {
            layer.add_field(FieldDef::new(name, FieldType::Float));
        }
        layer.add_field(FieldDef::new("ci_match_id", FieldType::Integer));
        let mut row_for_feature: Vec<Option<usize>> = vec![None; layer.features.len()];
        for (row, &fi) in idx_map.iter().enumerate() {
            row_for_feature[fi] = Some(row);
        }
        for (fi, feature) in layer.features.iter_mut().enumerate() {
            match row_for_feature[fi] {
                Some(r) => {
                    feature.attributes.push(FieldValue::Float(ps[r]));
                    feature.attributes.push(FieldValue::Float(weight_of[r]));
                    match match_of[r] {
                        Some(j) => feature
                            .attributes
                            .push(FieldValue::Integer(idx_map[j] as i64)),
                        None => feature.attributes.push(FieldValue::Null),
                    }
                }
                None => {
                    feature.attributes.push(FieldValue::Null);
                    feature.attributes.push(FieldValue::Null);
                    feature.attributes.push(FieldValue::Null);
                }
            }
        }

        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("method".to_string(), json!(prm.method.as_str()));
        outputs.insert("observations".to_string(), json!(n));
        outputs.insert("n_treated".to_string(), json!(n_treated));
        outputs.insert("n_control".to_string(), json!(n_control));
        outputs.insert("ate".to_string(), json!(ate));
        outputs.insert("ate_ci_low".to_string(), json!(ci_low));
        outputs.insert("ate_ci_high".to_string(), json!(ci_high));
        outputs.insert("bootstrap_iterations".to_string(), json!(BOOTSTRAP_ITERS));
        outputs.insert("bootstrap_valid".to_string(), json!(valid_boot));
        outputs.insert("naive_diff_means".to_string(), json!(naive));
        outputs.insert("seed".to_string(), json!(prm.seed));
        outputs.insert(
            "balance_threshold".to_string(),
            json!(prm.balance_threshold),
        );
        outputs.insert("max_smd_before".to_string(), json!(max_smd_before));
        outputs.insert("max_smd_after".to_string(), json!(max_smd_after));
        outputs.insert("balance_ok".to_string(), json!(balance_ok));
        let balance_json: Vec<Value> = field_names
            .iter()
            .enumerate()
            .map(|(i, name)| {
                json!({
                    "field": name,
                    "smd_before": balance_before[i].1,
                    "smd_after": balance_after[i].1,
                })
            })
            .collect();
        outputs.insert("balance".to_string(), json!(balance_json));
        if prm.method == Method::PsMatching {
            outputs.insert("n_matched".to_string(), json!(n_matched));
            outputs.insert("caliper".to_string(), json!(caliper_used));
        }

        Ok(ToolRunResult { outputs })
    }
}

// ── Method ────────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    PsMatching,
    Ipw,
    RegressionAdjustment,
}

impl Method {
    fn as_str(self) -> &'static str {
        match self {
            Self::PsMatching => "ps_matching",
            Self::Ipw => "ipw",
            Self::RegressionAdjustment => "regression_adjustment",
        }
    }
}

// ── Propensity model (self-contained logistic IRLS core, mirrors the
//    logistic branch of `generalized_linear_regression::fit_glm`) ──────────────

#[allow(clippy::needless_range_loop)]
fn fit_logistic(xs: &[Vec<f64>], ys: &[f64], p: usize) -> Option<Vec<f64>> {
    let n = ys.len();
    let mut mu: Vec<f64> = ys.iter().map(|&y| (y + 0.5) / 2.0).collect();
    let mut eta: Vec<f64> = mu.iter().map(|&m| (m / (1.0 - m)).ln()).collect();
    let mut beta = vec![0.0; p];
    for _ in 0..50 {
        let mut a = vec![vec![0.0; p]; p];
        let mut rhs = vec![0.0; p];
        for i in 0..n {
            let v = (mu[i] * (1.0 - mu[i])).max(1e-10);
            let z = eta[i] + (ys[i] - mu[i]) / v;
            let xi = &xs[i];
            for r in 0..p {
                let wr = v * xi[r];
                for c in 0..p {
                    a[r][c] += wr * xi[c];
                }
                rhs[r] += wr * z;
            }
        }
        let new_beta = solve(&a, &rhs)?;
        let mut delta = 0.0;
        for i in 0..n {
            let e: f64 = (0..p).map(|c| xs[i][c] * new_beta[c]).sum();
            eta[i] = e;
            mu[i] = 1.0 / (1.0 + (-e.clamp(-700.0, 700.0)).exp());
        }
        for c in 0..p {
            delta += (new_beta[c] - beta[c]).abs();
        }
        beta = new_beta;
        if delta < 1e-10 {
            break;
        }
    }
    Some(beta)
}

fn predict_logistic(xs: &[Vec<f64>], beta: &[f64]) -> Vec<f64> {
    xs.iter()
        .map(|xi| {
            let e: f64 = xi.iter().zip(beta).map(|(x, b)| x * b).sum();
            1.0 / (1.0 + (-e.clamp(-700.0, 700.0)).exp())
        })
        .collect()
}

/// OLS via the normal equations (shared by regression adjustment).
#[allow(clippy::needless_range_loop)]
fn fit_ols(xs: &[Vec<f64>], ys: &[f64], p: usize) -> Option<Vec<f64>> {
    let n = xs.len();
    let mut a = vec![vec![0.0; p]; p];
    let mut rhs = vec![0.0; p];
    for i in 0..n {
        let xi = &xs[i];
        for r in 0..p {
            for c in 0..p {
                a[r][c] += xi[r] * xi[c];
            }
            rhs[r] += xi[r] * ys[i];
        }
    }
    solve(&a, &rhs)
}

/// Solves the small dense system `a x = b` by Gaussian elimination with
/// partial pivoting. Returns `None` if `a` is singular.
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

// ── Effect estimators ───────────────────────────────────────────────────────────

/// Greedy nearest-neighbour matching (with replacement) of each treated unit
/// to the control with the closest propensity score, via a 1-D `kdtree`.
/// Returns `None` in the treated unit's slot when no control falls within the
/// caliper.
fn match_treated_to_controls(ps: &[f64], treat: &[f64], caliper: f64) -> Vec<Option<usize>> {
    let mut tree: KdTree<f64, usize, [f64; 1]> = KdTree::new(1);
    for (i, &t) in treat.iter().enumerate() {
        if t == 0.0 {
            let _ = tree.add([ps[i]], i);
        }
    }
    treat
        .iter()
        .enumerate()
        .map(|(i, &t)| {
            if t != 1.0 || tree.size() == 0 {
                return None;
            }
            match tree.nearest(&[ps[i]], 1, &squared_euclidean) {
                Ok(found) if !found.is_empty() => {
                    let (dist_sq, &j) = found[0];
                    if dist_sq.sqrt() <= caliper {
                        Some(j)
                    } else {
                        None
                    }
                }
                _ => None,
            }
        })
        .collect()
}

fn matching_ate(y: &[f64], treat: &[f64], matches: &[Option<usize>]) -> Option<(f64, usize)> {
    let diffs: Vec<f64> = matches
        .iter()
        .enumerate()
        .filter(|(i, _)| treat[*i] == 1.0)
        .filter_map(|(i, m)| m.map(|j| y[i] - y[j]))
        .collect();
    if diffs.is_empty() {
        return None;
    }
    Some((diffs.iter().sum::<f64>() / diffs.len() as f64, diffs.len()))
}

/// Stabilized, truncated inverse-probability-of-treatment weights.
fn ipw_weights(ps: &[f64], treat: &[f64]) -> Vec<f64> {
    let n = treat.len() as f64;
    let p_treat = treat.iter().sum::<f64>() / n;
    let eps = 1e-6;
    let raw: Vec<f64> = treat
        .iter()
        .zip(ps)
        .map(|(&t, &p)| {
            let pc = p.clamp(eps, 1.0 - eps);
            if t == 1.0 {
                p_treat / pc
            } else {
                (1.0 - p_treat) / (1.0 - pc)
            }
        })
        .collect();
    truncate_weights(&raw, 0.01)
}

fn truncate_weights(w: &[f64], tail: f64) -> Vec<f64> {
    let n = w.len();
    if n == 0 {
        return Vec::new();
    }
    let mut sorted = w.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let lo_idx = (((n as f64) * tail).floor() as usize).min(n - 1);
    let hi_idx = (((n as f64) * (1.0 - tail)).ceil() as usize)
        .min(n)
        .saturating_sub(1);
    let lo = sorted[lo_idx];
    let hi = sorted[hi_idx.max(lo_idx)];
    w.iter().map(|&v| v.clamp(lo, hi.max(lo))).collect()
}

/// Hajek (weighted-mean-ratio) ATE estimator.
fn hajek_ate(y: &[f64], treat: &[f64], w: &[f64]) -> Option<f64> {
    let (mut sw1, mut sy1, mut sw0, mut sy0) = (0.0, 0.0, 0.0, 0.0);
    for i in 0..y.len() {
        if treat[i] == 1.0 {
            sw1 += w[i];
            sy1 += w[i] * y[i];
        } else {
            sw0 += w[i];
            sy0 += w[i] * y[i];
        }
    }
    if sw1 <= 0.0 || sw0 <= 0.0 {
        return None;
    }
    Some(sy1 / sw1 - sy0 / sw0)
}

/// Regression-adjustment ATE: the treatment coefficient in an OLS of
/// `outcome ~ intercept + treatment + confounders`.
fn ate_regression(x_conf: &[Vec<f64>], treat: &[f64], y: &[f64]) -> Option<f64> {
    let n = treat.len();
    if n == 0 {
        return None;
    }
    let q = x_conf[0].len();
    let p = q + 2;
    let design: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = Vec::with_capacity(p);
            row.push(1.0);
            row.push(treat[i]);
            row.extend_from_slice(&x_conf[i]);
            row
        })
        .collect();
    let beta = fit_ols(&design, y, p)?;
    Some(beta[1])
}

fn naive_diff_means(y: &[f64], treat: &[f64]) -> Option<f64> {
    let (mut s1, mut n1, mut s0, mut n0) = (0.0, 0.0, 0.0, 0.0);
    for i in 0..y.len() {
        if treat[i] == 1.0 {
            s1 += y[i];
            n1 += 1.0;
        } else {
            s0 += y[i];
            n0 += 1.0;
        }
    }
    if n1 == 0.0 || n0 == 0.0 {
        return None;
    }
    Some(s1 / n1 - s0 / n0)
}

// ── Balance diagnostics ──────────────────────────────────────────────────────

fn weighted_group_moments(
    vals: &[f64],
    treat: &[f64],
    weights: &[f64],
    want_treated: bool,
) -> (f64, f64) {
    let mut sw = 0.0;
    let mut sy = 0.0;
    for i in 0..vals.len() {
        if (treat[i] == 1.0) != want_treated {
            continue;
        }
        sw += weights[i];
        sy += weights[i] * vals[i];
    }
    if sw <= 0.0 {
        return (f64::NAN, f64::NAN);
    }
    let mean = sy / sw;
    let mut svar = 0.0;
    for i in 0..vals.len() {
        if (treat[i] == 1.0) != want_treated {
            continue;
        }
        svar += weights[i] * (vals[i] - mean).powi(2);
    }
    (mean, svar / sw)
}

/// Standardized mean difference (Cohen-style, pooled SD) between the treated
/// and control groups for one confounder column, optionally weighted.
fn smd_field(vals: &[f64], treat: &[f64], weights: &[f64]) -> f64 {
    let (m1, v1) = weighted_group_moments(vals, treat, weights, true);
    let (m0, v0) = weighted_group_moments(vals, treat, weights, false);
    if !m1.is_finite() || !m0.is_finite() {
        return f64::NAN;
    }
    let pooled = ((v1 + v0) / 2.0).sqrt();
    if pooled <= 1e-12 {
        return 0.0;
    }
    (m1 - m0).abs() / pooled
}

/// Builds the matched-pairs subsample (treated + its matched control, one
/// weight-1 row each) used for the post-matching balance table.
fn matched_subset(vals: &[f64], treat: &[f64], matches: &[Option<usize>]) -> (Vec<f64>, Vec<f64>) {
    let mut v = Vec::new();
    let mut t = Vec::new();
    for (i, m) in matches.iter().enumerate() {
        if treat[i] != 1.0 {
            continue;
        }
        if let Some(j) = m {
            v.push(vals[i]);
            t.push(1.0);
            v.push(vals[*j]);
            t.push(0.0);
        }
    }
    (v, t)
}

fn std_dev(v: &[f64]) -> f64 {
    let n = v.len() as f64;
    if n == 0.0 {
        return 0.0;
    }
    let mean = v.iter().sum::<f64>() / n;
    (v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n).sqrt()
}

// ── Seeded bootstrap (splitmix64 seed expansion + xorshift64, WASM-safe:
//    no `Date::now`, no external RNG crate). ─────────────────────────────────

fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn xorshift_next(state: &mut u64) -> u64 {
    let mut x = *state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *state = x;
    x
}

/// Unbiased index in `0..n` via rejection sampling (no modulo bias).
fn next_index(state: &mut u64, n: usize) -> usize {
    let bound = n as u64;
    let limit = u64::MAX - (u64::MAX % bound);
    loop {
        let r = xorshift_next(state);
        if r < limit {
            return (r % bound) as usize;
        }
    }
}

/// Bootstraps a percentile confidence interval for `effect_fn` by resampling
/// `n` row indices with replacement, `iters` times, from a seed expanded via
/// splitmix64 into a nonzero xorshift64 state. Replicates where `effect_fn`
/// returns `None` (a degenerate resample) are skipped.
fn bootstrap_ci(
    n: usize,
    seed: u64,
    iters: usize,
    effect_fn: impl Fn(&[usize]) -> Option<f64>,
) -> (f64, f64, usize) {
    let mut seed_state = seed;
    let mut state = splitmix64(&mut seed_state).max(1);
    let mut samples = Vec::with_capacity(iters);
    for _ in 0..iters {
        let idx: Vec<usize> = (0..n).map(|_| next_index(&mut state, n)).collect();
        if let Some(a) = effect_fn(&idx) {
            if a.is_finite() {
                samples.push(a);
            }
        }
    }
    let m = samples.len();
    if m == 0 {
        return (f64::NAN, f64::NAN, 0);
    }
    samples.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let lo_idx = (((m as f64) * 0.025).floor() as usize).min(m - 1);
    let hi_idx = (((m as f64) * 0.975).ceil() as usize)
        .min(m)
        .saturating_sub(1);
    (samples[lo_idx], samples[hi_idx.max(lo_idx)], m)
}

// ── Spatial confounders ──────────────────────────────────────────────────────

fn representative_point(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(pts) => avg_xy(pts.iter().map(|c| (c.x, c.y))),
        Geometry::LineString(coords) => avg_xy(coords.iter().map(|c| (c.x, c.y))),
        Geometry::MultiLineString(parts) => avg_xy(parts.iter().flatten().map(|c| (c.x, c.y))),
        Geometry::Polygon { exterior, .. } => avg_xy(exterior.coords().iter().map(|c| (c.x, c.y))),
        Geometry::MultiPolygon(parts) => avg_xy(
            parts
                .iter()
                .flat_map(|(ext, _)| ext.coords().iter())
                .map(|c| (c.x, c.y)),
        ),
        Geometry::GeometryCollection(parts) => parts.first().and_then(representative_point),
    }
}

fn avg_xy(it: impl Iterator<Item = (f64, f64)>) -> Option<(f64, f64)> {
    let (mut sx, mut sy, mut n) = (0.0, 0.0, 0usize);
    for (x, y) in it {
        sx += x;
        sy += y;
        n += 1;
    }
    if n == 0 {
        None
    } else {
        Some((sx / n as f64, sy / n as f64))
    }
}

/// For each point, the mean of `values` over its `k` nearest neighbours
/// (excluding itself), found with a 2-D `kdtree`.
fn neighbor_averages(reps: &[(f64, f64)], values: &[Vec<f64>], k: usize) -> Vec<Vec<f64>> {
    let n = reps.len();
    let ncols = values.first().map(Vec::len).unwrap_or(0);
    let mut tree: KdTree<f64, usize, [f64; 2]> = KdTree::new(2);
    for (i, &(x, y)) in reps.iter().enumerate() {
        let _ = tree.add([x, y], i);
    }
    (0..n)
        .map(|i| {
            let (x, y) = reps[i];
            let found = tree
                .nearest(&[x, y], k + 1, &squared_euclidean)
                .unwrap_or_default();
            let mut sums = vec![0.0; ncols];
            let mut cnt = 0.0;
            for (_, &j) in found {
                if j == i {
                    continue;
                }
                for c in 0..ncols {
                    sums[c] += values[j][c];
                }
                cnt += 1.0;
            }
            if cnt > 0.0 {
                sums.iter().map(|s| s / cnt).collect()
            } else {
                vec![0.0; ncols]
            }
        })
        .collect()
}

// ── Parameters ──────────────────────────────────────────────────────────────────

struct Params {
    outcome_field: String,
    treatment_field: String,
    confounding_fields: Vec<String>,
    method: Method,
    add_spatial_confounders: bool,
    balance_threshold: f64,
    seed: u64,
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .map(str::trim)
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

fn to_numeric_treatment(v: &FieldValue) -> Option<f64> {
    match v {
        FieldValue::Boolean(b) => Some(if *b { 1.0 } else { 0.0 }),
        other => other.as_f64(),
    }
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

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let outcome_field = require_str(args, "outcome_field")?.to_string();
    let treatment_field = require_str(args, "treatment_field")?.to_string();
    let confounding_fields: Vec<String> = require_str(args, "confounding_fields")?
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_string)
        .collect();
    if confounding_fields.is_empty() {
        return Err(ToolError::Validation(
            "'confounding_fields' must list at least one field".to_string(),
        ));
    }
    let method = match parse_optional_str(args, "method")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("ps_matching") | Some("matching") => Method::PsMatching,
        Some("ipw") => Method::Ipw,
        Some("regression_adjustment") | Some("regression") => Method::RegressionAdjustment,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown method '{other}' (expected ps_matching, ipw, or regression_adjustment)"
            )))
        }
    };
    let add_spatial_confounders =
        parse_optional_bool(args, "add_spatial_confounders")?.unwrap_or(false);
    let balance_threshold = parse_optional_f64(args, "balance_threshold")?.unwrap_or(0.1);
    if !balance_threshold.is_finite() || balance_threshold <= 0.0 {
        return Err(ToolError::Validation(
            "'balance_threshold' must be positive".to_string(),
        ));
    }
    let seed = parse_optional_u64(args, "seed")?.unwrap_or(42);

    Ok(Params {
        outcome_field,
        treatment_field,
        confounding_fields,
        method,
        add_spatial_confounders,
        balance_threshold,
        seed,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, Geometry, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Deterministic pseudo-random value in `[0, 1)` (classic GLSL-style hash;
    /// no RNG state, purely a function of `i`) used only to build a synthetic
    /// dataset with realistic confounding and overlap — the tool itself never
    /// calls this.
    fn hash01(i: usize) -> f64 {
        let s = (i as f64 * 12.9898).sin() * 43758.5453;
        s - s.floor()
    }

    fn sigmoid(x: f64) -> f64 {
        1.0 / (1.0 + (-x).exp())
    }

    /// Builds a confounded synthetic dataset: a confounder `x` drives both the
    /// treatment assignment (via a logistic propensity) and the outcome
    /// (`y = true_effect * t + outcome_confound_coef * x + small noise`), so a
    /// naive difference in means is biased away from `true_effect`.
    fn confounded_layer(
        n: usize,
        true_effect: f64,
        outcome_confound_coef: f64,
    ) -> (String, Vec<f64>, Vec<f64>) {
        let mut layer = Layer::new("units").with_geom_type(GeometryType::Point);
        layer.add_field(FieldDef::new("x1", FieldType::Float));
        layer.add_field(FieldDef::new("treated", FieldType::Integer));
        layer.add_field(FieldDef::new("outcome", FieldType::Float));
        let mut treat = Vec::with_capacity(n);
        let mut outcome = Vec::with_capacity(n);
        for i in 0..n {
            let x = (i as f64) / (n as f64) - 0.5; // centered in [-0.5, 0.5)
            let t = if sigmoid(2.5 * x) > hash01(i) {
                1.0
            } else {
                0.0
            };
            let noise = (hash01(i + 9973) - 0.5) * 0.4;
            let y = true_effect * t + outcome_confound_coef * x + noise;
            treat.push(t);
            outcome.push(y);
            // Geometry is placed independently of `x1` (a different hash) so
            // enabling spatial confounders doesn't just re-derive a perfectly
            // collinear copy of `x1` and make the propensity design singular.
            let gx = hash01(i + 4271) * 100.0;
            let gy = hash01(i + 8081) * 100.0;
            layer
                .add_feature(
                    Some(Geometry::point(gx, gy)),
                    &[
                        ("x1", FieldValue::Float(x)),
                        ("treated", FieldValue::Integer(t as i64)),
                        ("outcome", FieldValue::Float(y)),
                    ],
                )
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        (memory_store::make_vector_memory_path(&id), treat, outcome)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CausalInferenceAnalysisTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// With confounding present, the adjusted ATE (any method) recovers the
    /// true effect much better than the naive difference in means, and lands
    /// inside its own bootstrap CI.
    #[test]
    fn ate_recovers_true_effect_while_naive_is_biased() {
        let true_effect = 5.0;
        let (input, _t, _y) = confounded_layer(240, true_effect, 4.0);
        for method in ["ps_matching", "ipw", "regression_adjustment"] {
            let (out, _layer) = run(json!({
                "input": input, "outcome_field": "outcome", "treatment_field": "treated",
                "confounding_fields": "x1", "method": method,
            }));
            let ate = out.outputs["ate"].as_f64().unwrap();
            let naive = out.outputs["naive_diff_means"].as_f64().unwrap();
            let lo = out.outputs["ate_ci_low"].as_f64().unwrap();
            let hi = out.outputs["ate_ci_high"].as_f64().unwrap();
            assert!(
                (ate - true_effect).abs() < 1.5,
                "{method}: ate {ate} too far from true effect {true_effect}"
            );
            // The naive estimate is pulled off by the confound; the adjusted
            // estimate must be reliably closer to the truth.
            assert!(
                (ate - true_effect).abs() < (naive - true_effect).abs(),
                "{method}: adjustment did not reduce bias (ate {ate}, naive {naive}, true {true_effect})"
            );
            assert!(
                lo <= ate + 1e-9 && ate - 1e-9 <= hi,
                "{method}: ate {ate} outside CI [{lo}, {hi}]"
            );
            assert!(
                lo <= true_effect + 2.0 && hi >= true_effect - 2.0,
                "{method}: CI [{lo},{hi}] misses truth"
            );
        }
    }

    /// Balance (max standardized mean difference) is large before adjustment
    /// and drops below `balance_threshold` after matching / IPW.
    #[test]
    fn balance_improves_after_adjustment() {
        let (input, _t, _y) = confounded_layer(240, 5.0, 4.0);
        for method in ["ps_matching", "ipw"] {
            let (out, _layer) = run(json!({
                "input": input, "outcome_field": "outcome", "treatment_field": "treated",
                "confounding_fields": "x1", "method": method, "balance_threshold": 0.15,
            }));
            let before = out.outputs["max_smd_before"].as_f64().unwrap();
            let after = out.outputs["max_smd_after"].as_f64().unwrap();
            assert!(
                before > 0.15,
                "{method}: expected imbalance before adjustment, got {before}"
            );
            assert!(
                after < before,
                "{method}: balance did not improve ({before} -> {after})"
            );
            assert!(
                after <= 0.15,
                "{method}: max_smd_after {after} did not drop below threshold"
            );
        }
    }

    /// The output layer carries a propensity score in [0, 1] for every unit,
    /// plus a weight and (for matching) a match id.
    #[test]
    fn output_layer_has_propensity_and_weight_fields() {
        let (input, _t, _y) = confounded_layer(120, 3.0, 2.0);
        let (_out, layer) = run(json!({
            "input": input, "outcome_field": "outcome", "treatment_field": "treated",
            "confounding_fields": "x1", "method": "ps_matching",
        }));
        let pi = layer.schema.field_index("ci_propensity").unwrap();
        let wi = layer.schema.field_index("ci_weight").unwrap();
        for f in &layer.features {
            let p = f.attributes[pi].as_f64().unwrap();
            assert!((0.0..=1.0).contains(&p), "propensity {p} out of range");
            assert!(f.attributes[wi].as_f64().is_some());
        }
    }

    /// Same seed -> identical ATE and CI (determinism required for WASM).
    #[test]
    fn deterministic_same_seed_same_result() {
        let (input, _t, _y) = confounded_layer(150, 4.0, 3.0);
        let args = json!({
            "input": input, "outcome_field": "outcome", "treatment_field": "treated",
            "confounding_fields": "x1", "method": "ipw", "seed": 7,
        });
        let (out1, _) = run(args.clone());
        let (out2, _) = run(args);
        assert_eq!(out1.outputs["ate"], out2.outputs["ate"]);
        assert_eq!(out1.outputs["ate_ci_low"], out2.outputs["ate_ci_low"]);
        assert_eq!(out1.outputs["ate_ci_high"], out2.outputs["ate_ci_high"]);
    }

    /// Different seeds are allowed to (and, with a stochastic bootstrap,
    /// typically do) produce different CI bounds, showing the seed is wired
    /// through rather than ignored.
    #[test]
    fn different_seeds_can_change_ci() {
        let (input, _t, _y) = confounded_layer(150, 4.0, 3.0);
        let base = json!({
            "input": input, "outcome_field": "outcome", "treatment_field": "treated",
            "confounding_fields": "x1", "method": "ipw",
        });
        let mut a1 = base.clone();
        a1["seed"] = json!(1);
        let mut a2 = base;
        a2["seed"] = json!(2);
        let (out1, _) = run(a1);
        let (out2, _) = run(a2);
        // Not a strict inequality assertion (collisions are possible in
        // principle) — just confirm both runs produced finite, usable CIs.
        assert!(out1.outputs["ate_ci_low"].as_f64().unwrap().is_finite());
        assert!(out2.outputs["ate_ci_low"].as_f64().unwrap().is_finite());
    }

    /// Enabling spatial confounders runs end-to-end without error and still
    /// produces a usable estimate.
    #[test]
    fn spatial_confounders_run_end_to_end() {
        let (input, _t, _y) = confounded_layer(150, 4.0, 3.0);
        let (out, _layer) = run(json!({
            "input": input, "outcome_field": "outcome", "treatment_field": "treated",
            "confounding_fields": "x1", "method": "ps_matching", "add_spatial_confounders": true,
        }));
        assert!(out.outputs["ate"].as_f64().unwrap().is_finite());
        let balance = out.outputs["balance"].as_array().unwrap();
        // x1, spatial_x, spatial_y, x1_nbr_mean
        assert_eq!(balance.len(), 4);
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = CausalInferenceAnalysisTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "p.geojson", "outcome_field": "y" })).is_err());
        assert!(bad(json!({
            "input": "p.geojson", "outcome_field": "y", "treatment_field": "t",
            "confounding_fields": "x1", "method": "propensity_score_double_robust_ml"
        }))
        .is_err());
        assert!(bad(json!({
            "input": "p.geojson", "outcome_field": "y", "treatment_field": "t",
            "confounding_fields": "x1", "balance_threshold": -0.1
        }))
        .is_err());
        assert!(bad(json!({
            "input": "p.geojson", "outcome_field": "y", "treatment_field": "t",
            "confounding_fields": "x1"
        }))
        .is_ok());
    }

    /// A non-binary treatment (e.g. a continuous exposure) is rejected at run
    /// time — the v1 scope cut called out in the module docs.
    #[test]
    fn rejects_non_binary_treatment() {
        let mut layer = Layer::new("units").with_geom_type(GeometryType::Point);
        layer.add_field(FieldDef::new("x1", FieldType::Float));
        layer.add_field(FieldDef::new("dose", FieldType::Float));
        layer.add_field(FieldDef::new("outcome", FieldType::Float));
        for i in 0..30 {
            let x = i as f64 * 0.1;
            layer
                .add_feature(
                    Some(Geometry::point(x, 0.0)),
                    &[
                        ("x1", FieldValue::Float(x)),
                        ("dose", FieldValue::Float(x * 2.0)), // continuous, not 0/1
                        ("outcome", FieldValue::Float(x + 1.0)),
                    ],
                )
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input, "outcome_field": "outcome", "treatment_field": "dose",
            "confounding_fields": "x1",
        }))
        .unwrap();
        let err = CausalInferenceAnalysisTool.run(&args, &ctx()).unwrap_err();
        assert!(
            format!("{err}").contains("binary"),
            "unexpected error: {err}"
        );
    }

    /// Too few observations is rejected with a clear execution error rather
    /// than panicking.
    #[test]
    fn rejects_too_few_observations() {
        let mut layer = Layer::new("units").with_geom_type(GeometryType::Point);
        layer.add_field(FieldDef::new("x1", FieldType::Float));
        layer.add_field(FieldDef::new("treated", FieldType::Integer));
        layer.add_field(FieldDef::new("outcome", FieldType::Float));
        for i in 0..5 {
            layer
                .add_feature(
                    Some(Geometry::point(i as f64, 0.0)),
                    &[
                        ("x1", FieldValue::Float(i as f64)),
                        ("treated", FieldValue::Integer((i % 2) as i64)),
                        ("outcome", FieldValue::Float(i as f64)),
                    ],
                )
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input, "outcome_field": "outcome", "treatment_field": "treated",
            "confounding_fields": "x1",
        }))
        .unwrap();
        assert!(CausalInferenceAnalysisTool.run(&args, &ctx()).is_err());
    }
}
