//! GeoLibre tool: generalized linear regression (global OLS / Poisson / logistic).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Generalized Linear Regression* (Spatial
//! Statistics): fit ONE global model relating a dependent field to explanatory
//! fields, with full model diagnostics. The catalog already ships a *local*
//! model (`geographically_weighted_regression`); this is its global companion.
//! whitebox's bundled regressions are raster image classifiers (random forest,
//! SVM, KNN over rasters), not attribute regression with inferential statistics.
//!
//! Three model families are supported via iteratively reweighted least squares
//! (IRLS) on the small `n × p` design matrix (intercept + explanatory terms):
//!
//! - `gaussian` — ordinary least squares (identity link); the Gaussian fit is an
//!   exact one-step solve of the normal equations `β = (XᵀX)⁻¹ Xᵀy`.
//! - `poisson` — log link for counts; working weights `μ`, response
//!   `η + (y−μ)/μ`.
//! - `logistic` — logit link for a 0/1 dependent; working weights `μ(1−μ)`.
//!
//! The normal-equation solve reuses the same Gaussian-elimination solver used by
//! `geographically_weighted_regression`.
//!
//! Each output feature keeps its attributes and gains `glr_estimated` (the fitted
//! value μ), `glr_residual` (raw `y − μ`) and `glr_std_resid` (the Pearson
//! residual `(y − μ)/√Var(μ)`). Features missing a value get nulls.
//!
//! Diagnostics (in the run outputs, and optionally a CSV `report`): per-term
//! coefficient, standard error, t/z statistic, two-sided probability, and the
//! variance-inflation factor (VIF) that flags multicollinearity; plus a
//! model block — AICc, R²/adjusted R² (Gaussian) or deviance & McFadden pseudo-R²
//! (Poisson/logistic), and the studentized Koenker (Breusch–Pagan) statistic for
//! heteroskedasticity. Residual spatial autocorrelation can then be checked by
//! handing `glr_residual` to the bundled Moran's I (left to the caller for v1).

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue};

use crate::common::write_text_output;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct GeneralizedLinearRegressionTool;

impl Tool for GeneralizedLinearRegressionTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "generalized_linear_regression",
            display_name: "Generalized Linear Regression",
            summary: "Global regression (gaussian/OLS, poisson, or logistic) of a dependent field on explanatory fields via IRLS, with fitted + residual output fields and full diagnostics: per-term coefficient, standard error, t/z, probability and VIF, plus AICc, R²/deviance and the studentized Koenker (Breusch–Pagan) heteroskedasticity test — like ArcGIS Generalized Linear Regression.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input feature layer (points, or other geometries — geometry is passed through unchanged).",
                    required: true,
                },
                ToolParamSpec {
                    name: "dependent_field",
                    description: "Dependent (response) numeric field. For 'logistic' it must be 0/1.",
                    required: true,
                },
                ToolParamSpec {
                    name: "explanatory_fields",
                    description: "Comma-separated explanatory numeric field(s).",
                    required: true,
                },
                ToolParamSpec {
                    name: "family",
                    description: "Model family: 'gaussian' (OLS, default), 'poisson' (log link, counts) or 'logistic' (logit link, 0/1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "report",
                    description: "Optional CSV path for the coefficient + model diagnostics report.",
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
        parse_family(args)?;
        let x_fields = parse_x_fields(args)?;
        if x_fields.is_empty() {
            return Err(ToolError::Validation(
                "'explanatory_fields' must list at least one field".to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let y_field = require_str(args, "dependent_field")?.to_string();
        let x_fields = parse_x_fields(args)?;
        if x_fields.is_empty() {
            return Err(ToolError::Validation(
                "'explanatory_fields' must list at least one field".to_string(),
            ));
        }
        let family = parse_family(args)?;
        let output = parse_optional_str(args, "output")?;
        let report = parse_optional_str(args, "report")?.map(str::to_string);

        let mut layer = load_input_layer(input)?;
        let schema = layer.schema.clone();

        // Build the design matrix (intercept first) from features with all values.
        let mut xs: Vec<Vec<f64>> = Vec::new();
        let mut ys: Vec<f64> = Vec::new();
        let mut idx_map: Vec<usize> = Vec::new();
        for (fi, feature) in layer.features.iter().enumerate() {
            let y = feature
                .get(&schema, &y_field)
                .ok()
                .and_then(FieldValue::as_f64);
            let mut row = Vec::with_capacity(x_fields.len() + 1);
            row.push(1.0);
            let mut ok = y.map(f64::is_finite).unwrap_or(false);
            for xf in &x_fields {
                match feature.get(&schema, xf).ok().and_then(FieldValue::as_f64) {
                    Some(v) if v.is_finite() => row.push(v),
                    _ => ok = false,
                }
            }
            if let (true, Some(y)) = (ok, y) {
                xs.push(row);
                ys.push(y);
                idx_map.push(fi);
            }
        }
        let n = ys.len();
        let p = x_fields.len() + 1;
        if n <= p + 1 {
            return Err(ToolError::Execution(format!(
                "need more than {} observations with valid values, found {n}",
                p + 1
            )));
        }
        if family == Family::Logistic && !ys.iter().all(|&y| y == 0.0 || y == 1.0) {
            return Err(ToolError::Execution(
                "'logistic' family requires a 0/1 dependent field".to_string(),
            ));
        }
        if family == Family::Poisson && ys.iter().any(|&y| y < 0.0) {
            return Err(ToolError::Execution(
                "'poisson' family requires a non-negative dependent field".to_string(),
            ));
        }

        ctx.progress
            .info(&format!("GLR ({}): {n} obs, {p} term(s)", family.as_str()));

        let fit = fit_glm(&xs, &ys, p, family).ok_or_else(|| {
            ToolError::Execution(
                "GLR fit failed (singular design matrix — check for collinear or constant fields)"
                    .to_string(),
            )
        })?;
        let diag = diagnostics(&xs, &ys, p, family, &fit, &x_fields);

        // Append per-feature output fields, keeping alignment with input features.
        for name in ["glr_estimated", "glr_residual", "glr_std_resid"] {
            layer.add_field(FieldDef::new(name, FieldType::Float));
        }
        let mut row_for_feature: Vec<Option<usize>> = vec![None; layer.features.len()];
        for (row, &fi) in idx_map.iter().enumerate() {
            row_for_feature[fi] = Some(row);
        }
        for (fi, feature) in layer.features.iter_mut().enumerate() {
            match row_for_feature[fi] {
                Some(r) => {
                    feature.attributes.push(FieldValue::Float(fit.mu[r]));
                    feature.attributes.push(FieldValue::Float(fit.residual[r]));
                    feature
                        .attributes
                        .push(FieldValue::Float(fit.pearson_resid[r]));
                }
                None => {
                    for _ in 0..3 {
                        feature.attributes.push(FieldValue::Null);
                    }
                }
            }
        }

        ctx.progress.info(&format!(
            "AICc {:.2}, {}",
            diag.aicc,
            match family {
                Family::Gaussian => format!("R2 {:.4}, adj R2 {:.4}", diag.r2, diag.adj_r2),
                _ => format!(
                    "deviance {:.2}, pseudo-R2 {:.4}",
                    diag.deviance, diag.pseudo_r2
                ),
            }
        ));

        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("family".to_string(), json!(family.as_str()));
        outputs.insert("observations".to_string(), json!(n));
        outputs.insert("terms".to_string(), json!(p));
        outputs.insert("iterations".to_string(), json!(fit.iterations));
        outputs.insert("aicc".to_string(), json!(diag.aicc));
        outputs.insert("aic".to_string(), json!(diag.aic));
        outputs.insert("log_likelihood".to_string(), json!(diag.log_likelihood));
        outputs.insert("deviance".to_string(), json!(diag.deviance));
        outputs.insert("null_deviance".to_string(), json!(diag.null_deviance));
        outputs.insert("residual_ss".to_string(), json!(diag.rss));
        if family == Family::Gaussian {
            outputs.insert("r2".to_string(), json!(diag.r2));
            outputs.insert("adjusted_r2".to_string(), json!(diag.adj_r2));
        } else {
            outputs.insert("pseudo_r2".to_string(), json!(diag.pseudo_r2));
        }
        outputs.insert("koenker_bp".to_string(), json!(diag.koenker_stat));
        outputs.insert("koenker_df".to_string(), json!(diag.koenker_df));
        outputs.insert("koenker_p".to_string(), json!(diag.koenker_p));
        outputs.insert(
            "max_vif".to_string(),
            json!(diag
                .terms
                .iter()
                .map(|t| t.vif)
                .filter(|v| v.is_finite())
                .fold(0.0_f64, f64::max)),
        );
        // Per-term diagnostics as a compact array the UI can table.
        let terms_json: Vec<Value> = diag
            .terms
            .iter()
            .map(|t| {
                json!({
                    "term": t.name,
                    "coefficient": t.coef,
                    "std_error": t.std_error,
                    "statistic": t.statistic,
                    "probability": t.probability,
                    "vif": t.vif,
                })
            })
            .collect();
        outputs.insert("coefficients".to_string(), json!(terms_json));

        if let Some(path) = report {
            write_text_output(&report_csv(&diag, family, n, p), &path)?;
            outputs.insert("report".to_string(), json!(path));
        }

        Ok(ToolRunResult { outputs })
    }
}

// ── Model family ───────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Family {
    Gaussian,
    Poisson,
    Logistic,
}

impl Family {
    fn as_str(self) -> &'static str {
        match self {
            Self::Gaussian => "gaussian",
            Self::Poisson => "poisson",
            Self::Logistic => "logistic",
        }
    }

    /// Inverse link: linear predictor η → mean μ.
    fn inv_link(self, eta: f64) -> f64 {
        match self {
            Self::Gaussian => eta,
            Self::Poisson => eta.clamp(-700.0, 700.0).exp(),
            Self::Logistic => 1.0 / (1.0 + (-eta.clamp(-700.0, 700.0)).exp()),
        }
    }

    /// Variance function V(μ).
    fn variance(self, mu: f64) -> f64 {
        match self {
            Self::Gaussian => 1.0,
            Self::Poisson => mu.max(1e-10),
            Self::Logistic => (mu * (1.0 - mu)).max(1e-10),
        }
    }
}

// ── IRLS fit ────────────────────────────────────────────────────────────────────

struct GlmFit {
    beta: Vec<f64>,
    mu: Vec<f64>,
    residual: Vec<f64>,      // raw y - μ
    pearson_resid: Vec<f64>, // (y - μ)/√V(μ)
    xtwx_inv: Vec<Vec<f64>>, // (XᵀWX)⁻¹ at convergence
    iterations: usize,
}

/// Fits a GLM by iteratively reweighted least squares. For the Gaussian family
/// this is a single exact solve of the normal equations. Returns `None` on a
/// singular design.
#[allow(clippy::needless_range_loop)]
fn fit_glm(xs: &[Vec<f64>], ys: &[f64], p: usize, family: Family) -> Option<GlmFit> {
    let n = ys.len();
    // Initial mean (kept strictly inside the link domain).
    let mut mu: Vec<f64> = ys
        .iter()
        .map(|&y| match family {
            Family::Gaussian => y,
            Family::Poisson => (y + 0.1).max(1e-3),
            Family::Logistic => (y + 0.5) / 2.0,
        })
        .collect();
    let mut eta: Vec<f64> = mu
        .iter()
        .map(|&m| match family {
            Family::Gaussian => m,
            Family::Poisson => m.ln(),
            Family::Logistic => (m / (1.0 - m)).ln(),
        })
        .collect();

    let mut beta = vec![0.0; p];
    let mut xtwx_inv = vec![vec![0.0; p]; p];
    let max_iter = if family == Family::Gaussian { 1 } else { 50 };
    let mut iterations = 0;

    for it in 0..max_iter {
        iterations = it + 1;
        // Working weights w and working response z (canonical links).
        // dμ/dη = V(μ) for canonical links, so w = V(μ), z = η + (y-μ)/V(μ).
        let mut a = vec![vec![0.0; p]; p];
        let mut rhs = vec![0.0; p];
        for i in 0..n {
            let (w, z) = match family {
                Family::Gaussian => (1.0, ys[i]),
                _ => {
                    let v = family.variance(mu[i]);
                    (v, eta[i] + (ys[i] - mu[i]) / v)
                }
            };
            let xi = &xs[i];
            for r in 0..p {
                let wr = w * xi[r];
                for c in 0..p {
                    a[r][c] += wr * xi[c];
                }
                rhs[r] += wr * z;
            }
        }
        let new_beta = solve(&a, &rhs)?;
        // Refresh η, μ.
        for i in 0..n {
            let e: f64 = (0..p).map(|c| xs[i][c] * new_beta[c]).sum();
            eta[i] = e;
            mu[i] = family.inv_link(e);
        }
        let delta: f64 = (0..p).map(|c| (new_beta[c] - beta[c]).abs()).sum();
        beta = new_beta;
        // Cache the inverse of the final weighted cross-product for covariance.
        xtwx_inv = invert(&a)?;
        if delta < 1e-10 {
            break;
        }
    }

    let residual: Vec<f64> = (0..n).map(|i| ys[i] - mu[i]).collect();
    let pearson_resid: Vec<f64> = (0..n)
        .map(|i| residual[i] / family.variance(mu[i]).sqrt())
        .collect();

    Some(GlmFit {
        beta,
        mu,
        residual,
        pearson_resid,
        xtwx_inv,
        iterations,
    })
}

// ── Diagnostics ──────────────────────────────────────────────────────────────────

struct TermDiag {
    name: String,
    coef: f64,
    std_error: f64,
    statistic: f64,
    probability: f64,
    vif: f64,
}

struct Diagnostics {
    terms: Vec<TermDiag>,
    aic: f64,
    aicc: f64,
    log_likelihood: f64,
    deviance: f64,
    null_deviance: f64,
    rss: f64,
    r2: f64,
    adj_r2: f64,
    pseudo_r2: f64,
    koenker_stat: f64,
    koenker_df: usize,
    koenker_p: f64,
}

#[allow(clippy::needless_range_loop)]
fn diagnostics(
    xs: &[Vec<f64>],
    ys: &[f64],
    p: usize,
    family: Family,
    fit: &GlmFit,
    x_fields: &[String],
) -> Diagnostics {
    let n = ys.len();
    let nf = n as f64;
    let rss: f64 = fit.residual.iter().map(|r| r * r).sum();

    // Dispersion: σ̂² for Gaussian (RSS/(n-p)); fixed 1 for Poisson/logistic.
    let dispersion = match family {
        Family::Gaussian => rss / (nf - p as f64).max(1.0),
        _ => 1.0,
    };

    // Standard errors from the (XᵀWX)⁻¹ diagonal scaled by dispersion.
    let mut terms = Vec::with_capacity(p);
    let term_names: Vec<String> = std::iter::once("Intercept".to_string())
        .chain(x_fields.iter().cloned())
        .collect();
    let use_t = family == Family::Gaussian;
    let resid_df = (nf - p as f64).max(1.0);
    for c in 0..p {
        let var = fit.xtwx_inv[c][c] * dispersion;
        let se = if var > 0.0 { var.sqrt() } else { f64::NAN };
        let stat = fit.beta[c] / se;
        let prob = if stat.is_finite() {
            if use_t {
                two_sided_t(stat.abs(), resid_df)
            } else {
                2.0 * (1.0 - normal_cdf(stat.abs()))
            }
        } else {
            f64::NAN
        };
        // VIF (auxiliary OLS of each explanatory on the others).
        let vif = if c == 0 {
            f64::NAN // intercept has no VIF
        } else {
            vif_for(xs, p, c)
        };
        terms.push(TermDiag {
            name: term_names[c].clone(),
            coef: fit.beta[c],
            std_error: se,
            statistic: stat,
            probability: prob,
            vif,
        });
    }

    // Log-likelihood, deviance and information criteria.
    let (log_likelihood, deviance, null_deviance, k) = match family {
        Family::Gaussian => {
            let ll = -0.5 * nf * ((2.0 * std::f64::consts::PI).ln() + (rss / nf).ln() + 1.0);
            let ybar = ys.iter().sum::<f64>() / nf;
            let tss: f64 = ys.iter().map(|y| (y - ybar).powi(2)).sum();
            // "Deviance" for Gaussian is the RSS; null deviance is the TSS.
            (ll, rss, tss, p + 1) // + variance parameter
        }
        Family::Poisson => {
            let ll: f64 = (0..n)
                .map(|i| {
                    let mu = fit.mu[i].max(1e-12);
                    ys[i] * mu.ln() - mu - ln_factorial(ys[i])
                })
                .sum();
            let dev: f64 = 2.0
                * (0..n)
                    .map(|i| {
                        let y = ys[i];
                        let mu = fit.mu[i].max(1e-12);
                        let t = if y > 0.0 { y * (y / mu).ln() } else { 0.0 };
                        t - (y - mu)
                    })
                    .sum::<f64>();
            let ybar = ys.iter().sum::<f64>() / nf;
            let null_dev: f64 = 2.0
                * (0..n)
                    .map(|i| {
                        let y = ys[i];
                        let t = if y > 0.0 { y * (y / ybar).ln() } else { 0.0 };
                        t - (y - ybar)
                    })
                    .sum::<f64>();
            (ll, dev, null_dev, p)
        }
        Family::Logistic => {
            let ll: f64 = (0..n)
                .map(|i| {
                    let mu = fit.mu[i].clamp(1e-12, 1.0 - 1e-12);
                    ys[i] * mu.ln() + (1.0 - ys[i]) * (1.0 - mu).ln()
                })
                .sum();
            let dev: f64 = -2.0
                * (0..n)
                    .map(|i| {
                        let mu = fit.mu[i].clamp(1e-12, 1.0 - 1e-12);
                        ys[i] * mu.ln() + (1.0 - ys[i]) * (1.0 - mu).ln()
                    })
                    .sum::<f64>();
            let ybar = (ys.iter().sum::<f64>() / nf).clamp(1e-12, 1.0 - 1e-12);
            let null_dev: f64 = -2.0
                * (0..n)
                    .map(|i| ys[i] * ybar.ln() + (1.0 - ys[i]) * (1.0 - ybar).ln())
                    .sum::<f64>();
            (ll, dev, null_dev, p)
        }
    };
    let kf = k as f64;
    let aic = -2.0 * log_likelihood + 2.0 * kf;
    let aicc = if nf - kf - 1.0 > 0.0 {
        aic + 2.0 * kf * (kf + 1.0) / (nf - kf - 1.0)
    } else {
        f64::INFINITY
    };

    let (r2, adj_r2, pseudo_r2) = match family {
        Family::Gaussian => {
            let r2 = if null_deviance > 0.0 {
                1.0 - rss / null_deviance
            } else {
                0.0
            };
            let adj = if nf - p as f64 > 0.0 {
                1.0 - (1.0 - r2) * (nf - 1.0) / (nf - p as f64)
            } else {
                r2
            };
            (r2, adj, f64::NAN)
        }
        _ => {
            let pseudo = if null_deviance > 0.0 {
                1.0 - deviance / null_deviance
            } else {
                0.0
            };
            (f64::NAN, f64::NAN, pseudo)
        }
    };

    // Studentized Koenker (Breusch–Pagan) test: regress squared residuals on the
    // explanatory design; statistic = n·R²_aux ~ χ² with (p-1) df.
    let (koenker_stat, koenker_df, koenker_p) = koenker_bp(xs, &fit.residual, p);

    Diagnostics {
        terms,
        aic,
        aicc,
        log_likelihood,
        deviance,
        null_deviance,
        rss,
        r2,
        adj_r2,
        pseudo_r2,
        koenker_stat,
        koenker_df,
        koenker_p,
    }
}

/// Variance inflation factor for explanatory column `c` (1-based within design,
/// i.e. c ∈ 1..p): 1/(1−R²) from an OLS of that column on the remaining
/// explanatory columns plus an intercept.
#[allow(clippy::needless_range_loop)]
fn vif_for(xs: &[Vec<f64>], p: usize, c: usize) -> f64 {
    if p <= 2 {
        return 1.0; // single explanatory term — nothing to be collinear with
    }
    let n = xs.len();
    // Predictors: intercept + all explanatory columns except c.
    let cols: Vec<usize> = (0..p).filter(|&j| j != c && j != 0).collect();
    let q = cols.len() + 1; // + intercept
    let mut a = vec![vec![0.0; q]; q];
    let mut rhs = vec![0.0; q];
    let mut sy = 0.0;
    let mut syy = 0.0;
    for i in 0..n {
        let mut row = Vec::with_capacity(q);
        row.push(1.0);
        for &j in &cols {
            row.push(xs[i][j]);
        }
        let target = xs[i][c];
        for r in 0..q {
            for cc in 0..q {
                a[r][cc] += row[r] * row[cc];
            }
            rhs[r] += row[r] * target;
        }
        sy += target;
        syy += target * target;
    }
    let Some(beta) = solve(&a, &rhs) else {
        return f64::NAN;
    };
    // R² of the auxiliary regression.
    let mut rss = 0.0;
    for i in 0..n {
        let mut yhat = beta[0];
        for (k, &j) in cols.iter().enumerate() {
            yhat += beta[k + 1] * xs[i][j];
        }
        rss += (xs[i][c] - yhat).powi(2);
    }
    let tss = syy - sy * sy / n as f64;
    if tss <= 0.0 {
        return f64::NAN;
    }
    let r2 = 1.0 - rss / tss;
    if r2 >= 1.0 {
        f64::INFINITY
    } else {
        1.0 / (1.0 - r2)
    }
}

/// Studentized Koenker (Breusch–Pagan) heteroskedasticity test.
#[allow(clippy::needless_range_loop)]
fn koenker_bp(xs: &[Vec<f64>], residual: &[f64], p: usize) -> (f64, usize, f64) {
    let n = xs.len();
    if n <= p || p < 2 {
        return (f64::NAN, 0, f64::NAN);
    }
    // Auxiliary OLS of e² on the full design (intercept + explanatory).
    let g: Vec<f64> = residual.iter().map(|r| r * r).collect();
    let mut a = vec![vec![0.0; p]; p];
    let mut rhs = vec![0.0; p];
    for i in 0..n {
        let xi = &xs[i];
        for r in 0..p {
            for c in 0..p {
                a[r][c] += xi[r] * xi[c];
            }
            rhs[r] += xi[r] * g[i];
        }
    }
    let Some(beta) = solve(&a, &rhs) else {
        return (f64::NAN, p - 1, f64::NAN);
    };
    let gbar = g.iter().sum::<f64>() / n as f64;
    let tss: f64 = g.iter().map(|v| (v - gbar).powi(2)).sum();
    let rss: f64 = (0..n)
        .map(|i| {
            let yhat: f64 = (0..p).map(|c| xs[i][c] * beta[c]).sum();
            (g[i] - yhat).powi(2)
        })
        .sum();
    if tss <= 0.0 {
        return (0.0, p - 1, 1.0);
    }
    let r2 = 1.0 - rss / tss;
    let stat = n as f64 * r2;
    let df = p - 1;
    let pval = chi2_sf(stat, df as f64);
    (stat, df, pval)
}

// ── CSV report ────────────────────────────────────────────────────────────────

fn report_csv(diag: &Diagnostics, family: Family, n: usize, p: usize) -> String {
    let mut s = String::new();
    s.push_str("section,name,value\n");
    s.push_str(&format!("model,family,{}\n", family.as_str()));
    s.push_str(&format!("model,observations,{n}\n"));
    s.push_str(&format!("model,terms,{p}\n"));
    s.push_str(&format!("model,aic,{:.6}\n", diag.aic));
    s.push_str(&format!("model,aicc,{:.6}\n", diag.aicc));
    s.push_str(&format!(
        "model,log_likelihood,{:.6}\n",
        diag.log_likelihood
    ));
    s.push_str(&format!("model,deviance,{:.6}\n", diag.deviance));
    s.push_str(&format!("model,null_deviance,{:.6}\n", diag.null_deviance));
    if family == Family::Gaussian {
        s.push_str(&format!("model,r2,{:.6}\n", diag.r2));
        s.push_str(&format!("model,adjusted_r2,{:.6}\n", diag.adj_r2));
        s.push_str(&format!("model,residual_ss,{:.6}\n", diag.rss));
    } else {
        s.push_str(&format!("model,pseudo_r2,{:.6}\n", diag.pseudo_r2));
    }
    s.push_str(&format!("koenker,statistic,{:.6}\n", diag.koenker_stat));
    s.push_str(&format!("koenker,df,{}\n", diag.koenker_df));
    s.push_str(&format!("koenker,probability,{:.6}\n", diag.koenker_p));
    // Coefficient table.
    s.push_str("term,coefficient,std_error,statistic,probability,vif\n");
    for t in &diag.terms {
        s.push_str(&format!(
            "{},{:.8},{:.8},{:.6},{:.6},{:.6}\n",
            t.name, t.coef, t.std_error, t.statistic, t.probability, t.vif
        ));
    }
    s
}

// ── Linear algebra (shared with GWR's solver) ───────────────────────────────────

/// Solves the small dense system `a x = b` by Gaussian elimination with partial
/// pivoting. Returns `None` if `a` is singular.
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

/// Inverts a small dense matrix by solving against each identity column.
#[allow(clippy::needless_range_loop)]
fn invert(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    let mut inv = vec![vec![0.0; n]; n];
    for c in 0..n {
        let mut e = vec![0.0; n];
        e[c] = 1.0;
        let col = solve(a, &e)?;
        for r in 0..n {
            inv[r][c] = col[r];
        }
    }
    Some(inv)
}

// ── Special functions (dependency-free) ─────────────────────────────────────────

/// Standard normal CDF via the error function.
fn normal_cdf(x: f64) -> f64 {
    0.5 * (1.0 + erf(x / std::f64::consts::SQRT_2))
}

/// Error function (Abramowitz & Stegun 7.1.26, |error| < 1.5e-7).
fn erf(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.327_591_1 * x);
    let y = 1.0
        - (((((1.061_405_429 * t - 1.453_152_027) * t) + 1.421_413_741) * t - 0.284_496_736) * t
            + 0.254_829_592)
            * t
            * (-x * x).exp();
    sign * y
}

/// Two-sided p-value for a Student-t statistic with `df` degrees of freedom,
/// via the regularized incomplete beta function.
fn two_sided_t(t: f64, df: f64) -> f64 {
    let x = df / (df + t * t);
    betai(df / 2.0, 0.5, x)
}

/// Upper-tail probability of the χ² distribution: P(X > x) for `df` degrees of
/// freedom, via the regularized upper incomplete gamma function.
fn chi2_sf(x: f64, df: f64) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    gammq(df / 2.0, x / 2.0)
}

/// Regularized incomplete beta I_x(a, b) (Numerical Recipes).
fn betai(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let bt = (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        bt * betacf(a, b, x) / a
    } else {
        1.0 - bt * betacf(b, a, 1.0 - x) / b
    }
}

fn betacf(a: f64, b: f64, x: f64) -> f64 {
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < 1e-30 {
        d = 1e-30;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..200 {
        let m = m as f64;
        let m2 = 2.0 * m;
        let mut aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < 1e-30 {
            d = 1e-30;
        }
        c = 1.0 + aa / c;
        if c.abs() < 1e-30 {
            c = 1e-30;
        }
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < 1e-30 {
            d = 1e-30;
        }
        c = 1.0 + aa / c;
        if c.abs() < 1e-30 {
            c = 1e-30;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-12 {
            break;
        }
    }
    h
}

/// Regularized upper incomplete gamma Q(a, x) = 1 − P(a, x) (Numerical Recipes).
fn gammq(a: f64, x: f64) -> f64 {
    if x < 0.0 || a <= 0.0 {
        return f64::NAN;
    }
    if x < a + 1.0 {
        1.0 - gser(a, x)
    } else {
        gcf(a, x)
    }
}

fn gser(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    let gln = ln_gamma(a);
    let mut ap = a;
    let mut sum = 1.0 / a;
    let mut del = sum;
    for _ in 0..500 {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * 1e-14 {
            break;
        }
    }
    sum * (-x + a * x.ln() - gln).exp()
}

fn gcf(a: f64, x: f64) -> f64 {
    let gln = ln_gamma(a);
    let mut b = x + 1.0 - a;
    let mut c = 1e30;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..500 {
        let an = -(i as f64) * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < 1e-30 {
            d = 1e-30;
        }
        c = b + an / c;
        if c.abs() < 1e-30 {
            c = 1e-30;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-14 {
            break;
        }
    }
    (-x + a * x.ln() - gln).exp() * h
}

/// Natural log of the gamma function (Lanczos approximation).
fn ln_gamma(x: f64) -> f64 {
    const G: [f64; 6] = [
        76.180_091_729_471_46,
        -86.505_320_329_416_77,
        24.014_098_240_830_91,
        -1.231_739_572_450_155,
        0.001_208_650_973_866_179,
        -0.000_005_395_239_384_953,
    ];
    let mut xx = x;
    let mut tmp = x + 5.5;
    tmp -= (x + 0.5) * tmp.ln();
    let mut ser = 1.000_000_000_190_015;
    for g in G {
        xx += 1.0;
        ser += g / xx;
    }
    -tmp + (2.506_628_274_631_000_5 * ser / x).ln()
}

/// ln(y!) for a non-negative (possibly non-integer) count via ln Γ(y+1).
fn ln_factorial(y: f64) -> f64 {
    ln_gamma(y + 1.0)
}

// ── Parameters ──────────────────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .map(str::trim)
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

fn parse_x_fields(args: &ToolArgs) -> Result<Vec<String>, ToolError> {
    Ok(require_str(args, "explanatory_fields")?
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_string)
        .collect())
}

fn parse_family(args: &ToolArgs) -> Result<Family, ToolError> {
    match parse_optional_str(args, "family")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("gaussian") | Some("ols") | Some("normal") => Ok(Family::Gaussian),
        Some("poisson") => Ok(Family::Poisson),
        Some("logistic") | Some("binomial") => Ok(Family::Logistic),
        Some(other) => Err(ToolError::Validation(format!(
            "unknown family '{other}' (expected gaussian, poisson or logistic)"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Geometry, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a point layer: each row is (x1, x2, y). x2 is Null when NaN.
    fn layer_with(rows: &[(f64, f64, f64)]) -> String {
        let mut layer = Layer::new("pts").with_geom_type(GeometryType::Point);
        layer.add_field(FieldDef::new("x1", FieldType::Float));
        layer.add_field(FieldDef::new("x2", FieldType::Float));
        layer.add_field(FieldDef::new("yv", FieldType::Float));
        for (i, &(x1, x2, yv)) in rows.iter().enumerate() {
            layer
                .add_feature(
                    Some(Geometry::point(i as f64, 0.0)),
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

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = GeneralizedLinearRegressionTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn coef(out: &ToolRunResult, name: &str) -> f64 {
        out.outputs["coefficients"]
            .as_array()
            .unwrap()
            .iter()
            .find(|t| t["term"] == name)
            .unwrap()["coefficient"]
            .as_f64()
            .unwrap()
    }

    /// Gaussian/OLS exactly recovers a known linear relationship y = 3 + 2·x1.
    #[test]
    fn gaussian_recovers_linear_relationship() {
        let mut rows = Vec::new();
        for i in 0..40 {
            let x1 = (i as f64 * 0.37).sin() * 5.0 + i as f64 * 0.1;
            rows.push((x1, 0.0, 3.0 + 2.0 * x1));
        }
        let (out, layer) = run(json!({
            "input": layer_with(&rows), "dependent_field": "yv",
            "explanatory_fields": "x1", "family": "gaussian",
        }));
        assert!(out.outputs["r2"].as_f64().unwrap() > 0.9999);
        assert!((coef(&out, "Intercept") - 3.0).abs() < 1e-6);
        assert!((coef(&out, "x1") - 2.0).abs() < 1e-6);
        // Output fields present and residuals ~0.
        let ri = layer.schema.field_index("glr_residual").unwrap();
        for f in &layer.features {
            assert!(f.attributes[ri].as_f64().unwrap().abs() < 1e-6);
        }
    }

    /// Logistic recovers the sign/direction of a separating relationship and
    /// classifies the training data well.
    #[test]
    fn logistic_fits_binary_outcome() {
        // y = 1 when x1 large. Add slight noise-free separation with overlap.
        let mut rows = Vec::new();
        for i in 0..60 {
            let x1 = i as f64 * 0.2 - 6.0;
            let y = if x1 + (i % 3) as f64 * 0.1 > 0.0 {
                1.0
            } else {
                0.0
            };
            rows.push((x1, 0.0, y));
        }
        let (out, layer) = run(json!({
            "input": layer_with(&rows), "dependent_field": "yv",
            "explanatory_fields": "x1", "family": "logistic",
        }));
        // Positive slope on x1 (higher x1 → higher probability of 1).
        assert!(coef(&out, "x1") > 0.0, "slope {}", coef(&out, "x1"));
        assert!(out.outputs["pseudo_r2"].as_f64().unwrap() > 0.5);
        // Fitted probabilities classify the training data with high accuracy.
        let ei = layer.schema.field_index("glr_estimated").unwrap();
        let yi = layer.schema.field_index("yv").unwrap();
        let correct = layer
            .features
            .iter()
            .filter(|f| {
                let p = f.attributes[ei].as_f64().unwrap();
                let y = f.attributes[yi].as_f64().unwrap();
                (p >= 0.5) == (y >= 0.5)
            })
            .count();
        assert!(correct as f64 / layer.features.len() as f64 > 0.9);
    }

    /// Poisson recovers a log-linear count relationship.
    #[test]
    fn poisson_recovers_log_linear_rate() {
        // log(μ) = 0.5 + 0.3·x1; counts are the rounded means (deterministic).
        let mut rows = Vec::new();
        for i in 0..50 {
            let x1 = i as f64 * 0.1;
            let mu = (0.5 + 0.3 * x1).exp();
            rows.push((x1, 0.0, mu.round()));
        }
        let (out, _l) = run(json!({
            "input": layer_with(&rows), "dependent_field": "yv",
            "explanatory_fields": "x1", "family": "poisson",
        }));
        // Slope should be near 0.3 (rounding introduces small error).
        assert!(
            (coef(&out, "x1") - 0.3).abs() < 0.1,
            "slope {}",
            coef(&out, "x1")
        );
        assert!(out.outputs["aicc"].as_f64().unwrap().is_finite());
    }

    /// VIF flags two nearly collinear explanatory fields.
    #[test]
    fn vif_flags_multicollinearity() {
        let mut rows = Vec::new();
        for i in 0..40 {
            let x1 = (i as f64 * 0.31).sin() * 3.0 + i as f64 * 0.05;
            let x2 = x1 * 2.0 + 0.001 * (i as f64 * 0.7).cos(); // ~perfectly collinear
            rows.push((x1, x2, 1.0 + x1 + 0.5 * x2));
        }
        let (out, _l) = run(json!({
            "input": layer_with(&rows), "dependent_field": "yv",
            "explanatory_fields": "x1,x2", "family": "gaussian",
        }));
        assert!(
            out.outputs["max_vif"].as_f64().unwrap() > 10.0,
            "max_vif {}",
            out.outputs["max_vif"]
        );
    }

    /// Features missing an explanatory value are passed through with null outputs.
    #[test]
    fn passes_through_missing_rows_as_null() {
        let mut rows: Vec<(f64, f64, f64)> = (0..30)
            .map(|i| {
                let x1 = i as f64 * 0.2;
                (x1, 0.0, 2.0 + 1.5 * x1)
            })
            .collect();
        // One row with a NaN explanatory becomes Null in the layer.
        rows.push((f64::NAN, 0.0, 99.0));
        let (out, layer) = run(json!({
            "input": layer_with(&rows), "dependent_field": "yv",
            "explanatory_fields": "x1", "family": "gaussian",
        }));
        assert_eq!(out.outputs["observations"].as_u64().unwrap(), 30);
        let ei = layer.schema.field_index("glr_estimated").unwrap();
        let last = layer.features.last().unwrap();
        assert!(last.attributes[ei].is_null());
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = GeneralizedLinearRegressionTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "p.geojson", "dependent_field": "yv" })).is_err());
        assert!(bad(json!({
            "input": "p.geojson", "dependent_field": "yv", "explanatory_fields": "x1",
            "family": "weibull"
        }))
        .is_err());
        assert!(bad(json!({
            "input": "p.geojson", "dependent_field": "yv", "explanatory_fields": "x1"
        }))
        .is_ok());
    }

    #[test]
    fn logistic_rejects_non_binary_dependent() {
        let rows: Vec<(f64, f64, f64)> = (0..30).map(|i| (i as f64 * 0.1, 0.0, i as f64)).collect();
        let args: ToolArgs = serde_json::from_value(json!({
            "input": layer_with(&rows), "dependent_field": "yv",
            "explanatory_fields": "x1", "family": "logistic",
        }))
        .unwrap();
        assert!(GeneralizedLinearRegressionTool.run(&args, &ctx()).is_err());
    }
}
