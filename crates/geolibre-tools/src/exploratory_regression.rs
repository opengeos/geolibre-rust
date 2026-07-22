//! GeoLibre tool: exploratory regression (search over candidate OLS models).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Exploratory Regression* (Spatial
//! Statistics). `generalized_linear_regression` fits ONE model the caller
//! specifies; finding a properly-specified model still means manually trying
//! variable subsets by hand. This tool automates that search: it enumerates
//! every `k`-sized combination of candidate explanatory fields (for `k` in
//! `[min_vars, max_vars]`), fits an OLS model for each, and screens it against
//! five diagnostics — adjusted R², coefficient significance, multicollinearity
//! (VIF), residual normality (Jarque–Bera) and residual spatial
//! autocorrelation (Moran's I) — reporting every evaluated model plus the
//! ranked subset that passes all five.
//!
//! The OLS core (normal-equation solve, standard errors/p-values via the
//! Student-t distribution, and the VIF auxiliary regression) reuses the exact
//! linear-algebra and diagnostic routines from `generalized_linear_regression`
//! (`solve`, `invert`, `vif_for`, `two_sided_t`, `chi2_sf`, `normal_cdf`) —
//! the Gaussian branch of that tool's IRLS fit is already an exact one-step
//! normal-equation solve, so this file only adds the small `Xᵀy`/`XᵀX`
//! assembly around those shared primitives. Residual spatial autocorrelation
//! reuses `incremental_spatial_autocorrelation`'s representative-point
//! extraction (`representative_xy`) and its Esri randomization-variance
//! formula for Moran's I, generalized from that tool's symmetric
//! fixed-distance weights to the asymmetric k-nearest-neighbour weights this
//! tool needs (general `S1`/`S2` from each feature's in-degree + out-degree,
//! which specializes back to their `S1=2·S0`, `S2=4·Σk²` shortcut for
//! symmetric binary weights).
//!
//! Because Moran's I needs an O(n²) neighbour search, it is computed only for
//! models that already pass the four cheap (O(n) or O(n·p²)) screens —
//! exactly as the issue asks. The candidate-combination search is capped by
//! an internal budget (`MAX_MODELS`); if the full `Σ C(n, k)` search would
//! exceed it, evaluation stops deterministically (ascending model size, then
//! lexicographic field order) and the exact number of *un*evaluated
//! combinations is logged and returned — never silently dropped.
//!
//! Output is a report table — one row per evaluated model with all
//! diagnostics and a `passed` flag (CSV when `output` ends in `.csv`,
//! otherwise a geometry-less vector table) — plus a ranked `passing_models`
//! list (by adjusted R², ties broken by lower VIF) and a per-candidate-field
//! `variable_summary`: how often each field was significant, with which sign.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Feature, FieldDef, FieldType, FieldValue, Layer};

use crate::common::write_text_output;
use crate::generalized_linear_regression::{
    chi2_sf, invert, normal_cdf, solve, two_sided_t, vif_for,
};
use crate::incremental_spatial_autocorrelation::representative_xy;
use crate::vector_common::{load_input_layer, write_or_store_layer};

/// Safety cap on the number of `Xᵀy`/VIF/JB model evaluations. Enumeration is
/// deterministic (ascending model size, then lexicographic field order), so
/// hitting the cap always drops the same, reported tail of combinations.
const MAX_MODELS: usize = 20_000;

pub struct ExploratoryRegressionTool;

impl Tool for ExploratoryRegressionTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "exploratory_regression",
            display_name: "Exploratory Regression",
            summary: "Enumerate OLS models over every combination of candidate explanatory fields and screen each on adjusted R2, coefficient significance, VIF multicollinearity, Jarque-Bera residual normality and residual Moran's I spatial autocorrelation, reporting every model plus the ranked passing set and a per-variable significance summary — like ArcGIS Exploratory Regression.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input feature layer with geometry (used for the residual Moran's I neighbourhood).",
                    required: true,
                },
                ToolParamSpec {
                    name: "dependent_field",
                    description: "Dependent (response) numeric field.",
                    required: true,
                },
                ToolParamSpec {
                    name: "explanatory_fields",
                    description: "Comma-separated candidate explanatory numeric fields to search over.",
                    required: true,
                },
                ToolParamSpec {
                    name: "min_vars",
                    description: "Minimum explanatory variables per candidate model (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_vars",
                    description: "Maximum explanatory variables per candidate model (default: min(candidate count, 6)).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_coef_p",
                    description: "Maximum allowed p-value for every explanatory coefficient (default 0.05).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_adj_r2",
                    description: "Minimum adjusted R-squared (default 0.5).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_vif",
                    description: "Maximum allowed variance inflation factor per term (default 7.5).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_jb_p",
                    description: "Minimum Jarque-Bera p-value for residual normality (default 0.1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_moran_p",
                    description: "Minimum p-value for residual spatial autocorrelation (Moran's I); above this the residuals are not significantly clustered (default 0.1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighbors",
                    description: "Number of nearest neighbours defining the spatial weights for the residual Moran's I test (default 8).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output report table path — a CSV (extension .csv) or a geometry-less vector table.",
                    required: true,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["input", "dependent_field", "explanatory_fields", "output"] {
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
        let y_field = require_str(args, "dependent_field")?;
        let candidates = parse_candidate_fields(args)?;
        if candidates.is_empty() {
            return Err(ToolError::Validation(
                "'explanatory_fields' must list at least one candidate field".to_string(),
            ));
        }
        if candidates.iter().any(|f| f == y_field) {
            return Err(ToolError::Validation(
                "'dependent_field' must not also appear in 'explanatory_fields'".to_string(),
            ));
        }
        parse_params(args, candidates.len())?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let y_field = require_str(args, "dependent_field")?.to_string();
        let candidates = parse_candidate_fields(args)?;
        let output = require_str(args, "output")?.to_string();
        let prm = parse_params(args, candidates.len())?;

        let layer = load_input_layer(input)?;
        let schema = layer.schema.clone();
        let y_idx = schema.field_index(&y_field).ok_or_else(|| {
            ToolError::Validation(format!("dependent field '{y_field}' not found"))
        })?;
        let mut cand_idx = Vec::with_capacity(candidates.len());
        for f in &candidates {
            let idx = schema.field_index(f).ok_or_else(|| {
                ToolError::Validation(format!("explanatory field '{f}' not found"))
            })?;
            cand_idx.push(idx);
        }

        // Complete-case sample: every candidate model is compared on the same
        // rows, so a row needs the dependent field, EVERY candidate field, and
        // a geometry (for the residual Moran's I neighbourhood).
        let mut ys: Vec<f64> = Vec::new();
        let mut xs_all: Vec<Vec<f64>> = Vec::new();
        let mut coords: Vec<(f64, f64)> = Vec::new();
        for feature in layer.iter() {
            let Some(y) = feature
                .attributes
                .get(y_idx)
                .and_then(FieldValue::as_f64)
                .filter(|v| v.is_finite())
            else {
                continue;
            };
            let mut row = Vec::with_capacity(candidates.len());
            let mut ok = true;
            for &idx in &cand_idx {
                match feature.attributes.get(idx).and_then(FieldValue::as_f64) {
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
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let Some(xy) = representative_xy(geom) else {
                continue;
            };
            ys.push(y);
            xs_all.push(row);
            coords.push(xy);
        }
        let n = ys.len();
        let min_needed = prm.max_vars + 3;
        if n <= min_needed {
            return Err(ToolError::Execution(format!(
                "need more than {min_needed} complete observations (dependent + all candidate fields + geometry), found {n}"
            )));
        }
        if prm.neighbors >= n {
            return Err(ToolError::Execution(format!(
                "'neighbors' ({}) must be less than the number of complete observations ({n})",
                prm.neighbors
            )));
        }

        let n_candidates = candidates.len();
        let total_possible: u128 = (prm.min_vars..=prm.max_vars)
            .map(|k| count_combinations(n_candidates, k))
            .sum();

        ctx.progress.info(&format!(
            "exploratory regression: {n} obs, {n_candidates} candidate field(s), {total_possible} candidate model(s) in [{}, {}] vars",
            prm.min_vars, prm.max_vars
        ));

        let mut results: Vec<ModelResult> = Vec::new();
        let mut budget_hit = false;
        'outer: for k in prm.min_vars..=prm.max_vars {
            let mut combo: Vec<usize> = (0..k).collect();
            loop {
                if results.len() >= MAX_MODELS {
                    budget_hit = true;
                    break 'outer;
                }
                results.push(evaluate_model(&combo, &xs_all, &ys, &coords, n, &prm));
                if !next_combination(&mut combo, n_candidates) {
                    break;
                }
            }
        }
        let dropped = if budget_hit {
            total_possible.saturating_sub(results.len() as u128)
        } else {
            0
        };
        if dropped > 0 {
            ctx.progress.info(&format!(
                "combination budget ({MAX_MODELS}) reached: {dropped} candidate model(s) were NOT evaluated (narrow min_vars/max_vars or the candidate field list to cover them)"
            ));
        }

        // ── Report table (one row per evaluated model) ──────────────────────
        let mut table = Layer::new("exploratory_regression");
        for (name, ty) in [
            ("model_id", FieldType::Integer),
            ("variables", FieldType::Text),
            ("n_vars", FieldType::Integer),
            ("observations", FieldType::Integer),
            ("adj_r2", FieldType::Float),
            ("r2", FieldType::Float),
            ("max_coef_p", FieldType::Float),
            ("max_vif", FieldType::Float),
            ("jb_stat", FieldType::Float),
            ("jb_p", FieldType::Float),
            ("moran_i", FieldType::Float),
            ("moran_z", FieldType::Float),
            ("moran_p", FieldType::Float),
            ("passed", FieldType::Boolean),
        ] {
            table.add_field(FieldDef::new(name, ty));
        }

        let mut csv = String::from(
            "model_id,variables,n_vars,observations,adj_r2,r2,max_coef_p,max_vif,jb_stat,jb_p,moran_i,moran_z,moran_p,passed\n",
        );
        let mut passed_count = 0usize;
        for (id, r) in results.iter().enumerate() {
            let vars_str = r
                .combo
                .iter()
                .map(|&i| candidates[i].as_str())
                .collect::<Vec<_>>()
                .join("|");
            if r.passed {
                passed_count += 1;
            }
            table.push(Feature {
                fid: id as u64,
                geometry: None,
                attributes: vec![
                    FieldValue::Integer(id as i64),
                    FieldValue::Text(vars_str.clone()),
                    FieldValue::Integer(r.combo.len() as i64),
                    FieldValue::Integer(r.n as i64),
                    FieldValue::Float(r.adj_r2),
                    FieldValue::Float(r.r2),
                    FieldValue::Float(r.max_coef_p),
                    FieldValue::Float(r.max_vif),
                    FieldValue::Float(r.jb_stat),
                    FieldValue::Float(r.jb_p),
                    opt_float(r.moran_i),
                    opt_float(r.moran_z),
                    opt_float(r.moran_p),
                    FieldValue::Boolean(r.passed),
                ],
            });
            csv.push_str(&format!(
                "{id},{vars_str},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{},{},{},{}\n",
                r.combo.len(),
                r.n,
                r.adj_r2,
                r.r2,
                r.max_coef_p,
                r.max_vif,
                r.jb_stat,
                r.jb_p,
                fmt_opt(r.moran_i),
                fmt_opt(r.moran_z),
                fmt_opt(r.moran_p),
                r.passed
            ));
        }

        let out_path = if output.to_ascii_lowercase().ends_with(".csv") {
            write_text_output(&csv, &output)?;
            output.clone()
        } else {
            write_or_store_layer(table, Some(output.as_str()))?
        };

        // ── Ranked passing models (by adjusted R2 desc, ties by lower VIF) ──
        let mut passing_idx: Vec<usize> =
            (0..results.len()).filter(|&i| results[i].passed).collect();
        passing_idx.sort_by(|&a, &b| {
            results[b]
                .adj_r2
                .total_cmp(&results[a].adj_r2)
                .then(results[a].max_vif.total_cmp(&results[b].max_vif))
        });
        let passing_models: Vec<Value> = passing_idx
            .iter()
            .enumerate()
            .map(|(rank, &i)| {
                let r = &results[i];
                json!({
                    "rank": rank + 1,
                    "variables": r.combo.iter().map(|&ci| candidates[ci].clone()).collect::<Vec<_>>(),
                    "adj_r2": r.adj_r2,
                    "r2": r.r2,
                    "max_coef_p": r.max_coef_p,
                    "max_vif": r.max_vif,
                    "jb_p": r.jb_p,
                    "moran_p": r.moran_p,
                })
            })
            .collect();

        // ── Per-variable summary ─────────────────────────────────────────────
        // (times_in_model, times_significant, positive_count, negative_count)
        let mut var_stats: Vec<(usize, usize, usize, usize)> = vec![(0, 0, 0, 0); n_candidates];
        let mut var_in_passing = vec![0usize; n_candidates];
        for r in &results {
            if r.passed {
                for &ci in &r.combo {
                    var_in_passing[ci] += 1;
                }
            }
            for t in &r.terms {
                let e = &mut var_stats[t.field_idx];
                e.0 += 1;
                if t.p.is_finite() && t.p <= prm.max_coef_p {
                    e.1 += 1;
                    if t.coef > 0.0 {
                        e.2 += 1;
                    } else if t.coef < 0.0 {
                        e.3 += 1;
                    }
                }
            }
        }
        let variable_summary: Vec<Value> = candidates
            .iter()
            .enumerate()
            .map(|(i, name)| {
                let (times, sig, pos, neg) = var_stats[i];
                let pct_sig = if times > 0 {
                    sig as f64 / times as f64
                } else {
                    0.0
                };
                let dominant_sign = if sig == 0 {
                    "none"
                } else if pos > 0 && neg == 0 {
                    "positive"
                } else if neg > 0 && pos == 0 {
                    "negative"
                } else {
                    "mixed"
                };
                json!({
                    "field": name,
                    "times_evaluated": times,
                    "times_significant": sig,
                    "pct_significant": pct_sig,
                    "positive_count": pos,
                    "negative_count": neg,
                    "dominant_sign": dominant_sign,
                    "times_in_passing_model": var_in_passing[i],
                })
            })
            .collect();

        ctx.progress.info(&format!(
            "{} model(s) evaluated, {passed_count} passed all thresholds",
            results.len()
        ));

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("observations".to_string(), json!(n));
        outputs.insert("candidate_fields".to_string(), json!(n_candidates));
        outputs.insert(
            "models_possible".to_string(),
            json!(total_possible.to_string()),
        );
        outputs.insert("models_evaluated".to_string(), json!(results.len()));
        outputs.insert("models_dropped".to_string(), json!(dropped.to_string()));
        outputs.insert("models_passed".to_string(), json!(passed_count));
        outputs.insert("passing_models".to_string(), json!(passing_models));
        outputs.insert("variable_summary".to_string(), json!(variable_summary));
        // Keep the inlined table small; the full report always lives at `output`.
        if results.len() <= 1000 {
            let table_json: Vec<Value> = results
                .iter()
                .enumerate()
                .map(|(id, r)| {
                    json!({
                        "model_id": id,
                        "variables": r.combo.iter().map(|&ci| candidates[ci].clone()).collect::<Vec<_>>(),
                        "adj_r2": r.adj_r2,
                        "max_coef_p": r.max_coef_p,
                        "max_vif": r.max_vif,
                        "jb_p": r.jb_p,
                        "moran_p": r.moran_p,
                        "passed": r.passed,
                    })
                })
                .collect();
            outputs.insert("table".to_string(), json!(table_json));
        }

        Ok(ToolRunResult { outputs })
    }
}

// ── Model evaluation ─────────────────────────────────────────────────────────

struct TermResult {
    field_idx: usize,
    coef: f64,
    p: f64,
}

struct ModelResult {
    combo: Vec<usize>,
    n: usize,
    adj_r2: f64,
    r2: f64,
    max_coef_p: f64,
    max_vif: f64,
    jb_stat: f64,
    jb_p: f64,
    moran_i: Option<f64>,
    moran_z: Option<f64>,
    moran_p: Option<f64>,
    passed: bool,
    terms: Vec<TermResult>,
}

/// Fits and screens the OLS model over `combo` (candidate-field indices).
/// Always returns a row (never silently dropped) — a singular or degenerate
/// design comes back with NaN diagnostics and `passed = false`.
#[allow(clippy::too_many_arguments)]
fn evaluate_model(
    combo: &[usize],
    xs_all: &[Vec<f64>],
    ys: &[f64],
    coords: &[(f64, f64)],
    n: usize,
    prm: &Params,
) -> ModelResult {
    let p = combo.len() + 1;
    let failed = |terms: Vec<TermResult>| ModelResult {
        combo: combo.to_vec(),
        n,
        adj_r2: f64::NAN,
        r2: f64::NAN,
        max_coef_p: f64::NAN,
        max_vif: f64::NAN,
        jb_stat: f64::NAN,
        jb_p: f64::NAN,
        moran_i: None,
        moran_z: None,
        moran_p: None,
        passed: false,
        terms,
    };

    let mut xs = vec![vec![0.0; p]; n];
    for i in 0..n {
        xs[i][0] = 1.0;
        for (j, &ci) in combo.iter().enumerate() {
            xs[i][j + 1] = xs_all[i][ci];
        }
    }
    let Some(fit) = ols_fit(&xs, ys, p) else {
        return failed(Vec::new());
    };
    let rss: f64 = fit.residual.iter().map(|r| r * r).sum();
    let ybar = ys.iter().sum::<f64>() / n as f64;
    let tss: f64 = ys.iter().map(|y| (y - ybar).powi(2)).sum();
    if tss <= 0.0 {
        return failed(Vec::new());
    }
    let r2 = 1.0 - rss / tss;
    let dof = (n as f64 - p as f64).max(1.0);
    let adj_r2 = 1.0 - (1.0 - r2) * (n as f64 - 1.0) / dof;
    let dispersion = rss / dof;

    let mut max_coef_p = 0.0_f64;
    let mut max_vif = 1.0_f64;
    let mut terms = Vec::with_capacity(combo.len());
    for c in 1..p {
        let var = fit.xtx_inv[c][c] * dispersion;
        let se = if var > 0.0 { var.sqrt() } else { f64::NAN };
        let t = fit.beta[c] / se;
        let pval = if t.is_finite() {
            two_sided_t(t.abs(), dof)
        } else {
            f64::NAN
        };
        let vif = vif_for(&xs, p, c);
        max_coef_p = max_coef_p.max(if pval.is_finite() {
            pval
        } else {
            f64::INFINITY
        });
        max_vif = max_vif.max(if vif.is_finite() { vif } else { f64::INFINITY });
        terms.push(TermResult {
            field_idx: combo[c - 1],
            coef: fit.beta[c],
            p: pval,
        });
    }

    let (jb_stat, jb_p) = jarque_bera(&fit.residual);

    let cheap_pass = adj_r2.is_finite()
        && adj_r2 >= prm.min_adj_r2
        && max_coef_p.is_finite()
        && max_coef_p <= prm.max_coef_p
        && max_vif.is_finite()
        && max_vif <= prm.max_vif
        && jb_p.is_finite()
        && jb_p >= prm.min_jb_p;

    let (moran_i, moran_z, moran_p) = if cheap_pass {
        let m = moran_i_knn(coords, &fit.residual, prm.neighbors);
        (Some(m.i), Some(m.z), Some(m.p))
    } else {
        (None, None, None)
    };
    let passed = cheap_pass
        && moran_p
            .map(|p| p.is_finite() && p >= prm.min_moran_p)
            .unwrap_or(false);

    ModelResult {
        combo: combo.to_vec(),
        n,
        adj_r2,
        r2,
        max_coef_p,
        max_vif,
        jb_stat,
        jb_p,
        moran_i,
        moran_z,
        moran_p,
        passed,
        terms,
    }
}

// ── OLS core (reuses `solve`/`invert` from generalized_linear_regression) ───

struct OlsFit {
    beta: Vec<f64>,
    residual: Vec<f64>,
    xtx_inv: Vec<Vec<f64>>,
}

/// Exact one-step normal-equation OLS solve — the same computation as the
/// Gaussian branch of `generalized_linear_regression::fit_glm`, assembled
/// directly here since that function isn't itself reusable (it also runs
/// IRLS for the other GLM families we don't need).
#[allow(clippy::needless_range_loop)]
fn ols_fit(xs: &[Vec<f64>], ys: &[f64], p: usize) -> Option<OlsFit> {
    let n = xs.len();
    let mut xtx = vec![vec![0.0; p]; p];
    let mut xty = vec![0.0; p];
    for i in 0..n {
        let xi = &xs[i];
        for r in 0..p {
            for c in 0..p {
                xtx[r][c] += xi[r] * xi[c];
            }
            xty[r] += xi[r] * ys[i];
        }
    }
    let beta = solve(&xtx, &xty)?;
    let xtx_inv = invert(&xtx)?;
    let residual: Vec<f64> = (0..n).map(|i| ys[i] - dot(&xs[i], &beta)).collect();
    Some(OlsFit {
        beta,
        residual,
        xtx_inv,
    })
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Jarque-Bera normality statistic and its chi-squared(2) p-value.
fn jarque_bera(residual: &[f64]) -> (f64, f64) {
    let n = residual.len() as f64;
    if n < 8.0 {
        return (f64::NAN, f64::NAN);
    }
    let mean = residual.iter().sum::<f64>() / n;
    let dev: Vec<f64> = residual.iter().map(|r| r - mean).collect();
    let m2 = dev.iter().map(|d| d * d).sum::<f64>() / n;
    if m2 <= 0.0 {
        return (f64::NAN, f64::NAN);
    }
    let m3 = dev.iter().map(|d| d.powi(3)).sum::<f64>() / n;
    let m4 = dev.iter().map(|d| d.powi(4)).sum::<f64>() / n;
    let skew = m3 / m2.powf(1.5);
    let kurt = m4 / (m2 * m2);
    let jb = (n / 6.0) * (skew * skew + (kurt - 3.0).powi(2) / 4.0);
    (jb, chi2_sf(jb, 2.0))
}

// ── Residual Moran's I over k-nearest-neighbour weights ─────────────────────

struct MoranResult {
    i: f64,
    z: f64,
    p: f64,
}

/// Global Moran's I of `residual` over a binary k-nearest-neighbour weight
/// matrix (asymmetric in general — `i` neighbouring `j` doesn't imply the
/// reverse), with Esri's randomization-variance z-score. This is the same
/// `S0/S1/S2` + kurtosis formula `incremental_spatial_autocorrelation` uses
/// for symmetric fixed-distance weights, generalized to asymmetric weights
/// via each feature's row-sum (out-degree) + column-sum (in-degree):
/// `S1 = 1/2 Σ(w_ij + w_ji)²`, `S2 = Σ(w_i· + w_·i)²`. It specializes back to
/// their `S1 = 2·S0`, `S2 = 4·Σkᵢ²` shortcut when the weights are symmetric
/// binary (their fixed-distance case).
#[allow(clippy::needless_range_loop)]
fn moran_i_knn(coords: &[(f64, f64)], residual: &[f64], k: usize) -> MoranResult {
    let n = residual.len();
    let nan = MoranResult {
        i: f64::NAN,
        z: f64::NAN,
        p: f64::NAN,
    };
    if k == 0 || n <= k + 1 {
        return nan;
    }
    let nf = n as f64;
    let mean = residual.iter().sum::<f64>() / nf;
    let dev: Vec<f64> = residual.iter().map(|r| r - mean).collect();
    let m2: f64 = dev.iter().map(|d| d * d).sum();
    let m4: f64 = dev.iter().map(|d| d.powi(4)).sum();
    if m2 <= 0.0 {
        return nan;
    }

    // Binary kNN weights: w[i][j] = 1 when j is one of i's k nearest.
    let mut w = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        let mut dists: Vec<(f64, usize)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| {
                let dx = coords[i].0 - coords[j].0;
                let dy = coords[i].1 - coords[j].1;
                ((dx * dx + dy * dy).sqrt(), j)
            })
            .collect();
        dists.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        for &(_, j) in dists.iter().take(k) {
            w[i][j] = 1.0;
        }
    }

    let s0 = (n * k) as f64;
    let mut cross = 0.0;
    let mut row_sum = vec![0.0_f64; n];
    let mut col_sum = vec![0.0_f64; n];
    for i in 0..n {
        for j in 0..n {
            let wij = w[i][j];
            if wij != 0.0 {
                cross += wij * dev[i] * dev[j];
                row_sum[i] += wij;
                col_sum[j] += wij;
            }
        }
    }
    let i_val = (nf / s0) * (cross / m2);

    let mut s1 = 0.0;
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let s = w[i][j] + w[j][i];
            if s != 0.0 {
                s1 += 0.5 * s * s;
            }
        }
    }
    let s2: f64 = (0..n).map(|i| (row_sum[i] + col_sum[i]).powi(2)).sum();

    let kurt = nf * m4 / (m2 * m2);
    let e_i = -1.0 / (nf - 1.0);
    let n2 = nf * nf;
    let a = nf * ((n2 - 3.0 * nf + 3.0) * s1 - nf * s2 + 3.0 * s0 * s0);
    let b = kurt * ((n2 - nf) * s1 - 2.0 * nf * s2 + 6.0 * s0 * s0);
    let denom = (nf - 1.0) * (nf - 2.0) * (nf - 3.0) * s0 * s0;
    let var = if denom != 0.0 {
        (a - b) / denom - e_i * e_i
    } else {
        f64::NAN
    };
    let (z, p) = if var > 0.0 {
        let z = (i_val - e_i) / var.sqrt();
        (z, 2.0 * (1.0 - normal_cdf(z.abs())))
    } else {
        (f64::NAN, f64::NAN)
    };
    MoranResult { i: i_val, z, p }
}

// ── Combinatorics ────────────────────────────────────────────────────────────

/// Advances `comb` (a strictly increasing 0-based index tuple of size `k`,
/// chosen from `0..n`) to the next combination in lexicographic order.
/// Returns `false` once `comb` was the last one.
fn next_combination(comb: &mut [usize], n: usize) -> bool {
    let k = comb.len();
    for i in (0..k).rev() {
        if comb[i] < n - k + i {
            comb[i] += 1;
            for j in i + 1..k {
                comb[j] = comb[j - 1] + 1;
            }
            return true;
        }
    }
    false
}

/// `C(n, k)` computed iteratively in `u128` (exact at every step since each
/// partial product is itself a binomial coefficient).
fn count_combinations(n: usize, k: usize) -> u128 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut result: u128 = 1;
    for i in 0..k {
        result = result * (n - i) as u128 / (i + 1) as u128;
    }
    result
}

// ── Field-value helpers ──────────────────────────────────────────────────────

fn opt_float(v: Option<f64>) -> FieldValue {
    match v {
        Some(x) if x.is_finite() => FieldValue::Float(x),
        _ => FieldValue::Null,
    }
}

fn fmt_opt(v: Option<f64>) -> String {
    match v {
        Some(x) if x.is_finite() => format!("{x:.6}"),
        _ => String::new(),
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    min_vars: usize,
    max_vars: usize,
    max_coef_p: f64,
    min_adj_r2: f64,
    max_vif: f64,
    min_jb_p: f64,
    min_moran_p: f64,
    neighbors: usize,
}

fn parse_params(args: &ToolArgs, n_candidates: usize) -> Result<Params, ToolError> {
    let min_vars = match parse_optional_u64(args, "min_vars")? {
        None => 1,
        Some(v) if v >= 1 => v as usize,
        Some(_) => return Err(ToolError::Validation("'min_vars' must be >= 1".to_string())),
    };
    let default_max = n_candidates.min(6).max(min_vars.min(n_candidates));
    let max_vars = match parse_optional_u64(args, "max_vars")? {
        None => default_max,
        Some(v) if v >= 1 => v as usize,
        Some(_) => return Err(ToolError::Validation("'max_vars' must be >= 1".to_string())),
    };
    if min_vars > max_vars {
        return Err(ToolError::Validation(
            "'min_vars' must be <= 'max_vars'".to_string(),
        ));
    }
    if max_vars > n_candidates {
        return Err(ToolError::Validation(format!(
            "'max_vars' ({max_vars}) cannot exceed the number of candidate explanatory fields ({n_candidates})"
        )));
    }

    let max_coef_p = parse_optional_f64(args, "max_coef_p")?.unwrap_or(0.05);
    if max_coef_p <= 0.0 || max_coef_p > 1.0 {
        return Err(ToolError::Validation(
            "'max_coef_p' must be in (0, 1]".to_string(),
        ));
    }
    let min_adj_r2 = parse_optional_f64(args, "min_adj_r2")?.unwrap_or(0.5);
    if min_adj_r2 > 1.0 || min_adj_r2.is_nan() {
        return Err(ToolError::Validation(
            "'min_adj_r2' must be <= 1".to_string(),
        ));
    }
    let max_vif = parse_optional_f64(args, "max_vif")?.unwrap_or(7.5);
    if max_vif < 1.0 || max_vif.is_nan() {
        return Err(ToolError::Validation("'max_vif' must be >= 1".to_string()));
    }
    let min_jb_p = parse_optional_f64(args, "min_jb_p")?.unwrap_or(0.1);
    if !(0.0..=1.0).contains(&min_jb_p) {
        return Err(ToolError::Validation(
            "'min_jb_p' must be in [0, 1]".to_string(),
        ));
    }
    let min_moran_p = parse_optional_f64(args, "min_moran_p")?.unwrap_or(0.1);
    if !(0.0..=1.0).contains(&min_moran_p) {
        return Err(ToolError::Validation(
            "'min_moran_p' must be in [0, 1]".to_string(),
        ));
    }
    let neighbors = match parse_optional_u64(args, "neighbors")? {
        None => 8,
        Some(v) if v >= 1 => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "'neighbors' must be >= 1".to_string(),
            ))
        }
    };

    Ok(Params {
        min_vars,
        max_vars,
        max_coef_p,
        min_adj_r2,
        max_vif,
        min_jb_p,
        min_moran_p,
        neighbors,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn parse_candidate_fields(args: &ToolArgs) -> Result<Vec<String>, ToolError> {
    Ok(require_str(args, "explanatory_fields")?
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_string)
        .collect())
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
    use wbvector::{memory_store, FieldValue as FV, Geometry, GeometryType, Layer as WbLayer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// A layer with a scrambled point scatter (decorrelates spatial adjacency
    /// from row index) and the given numeric fields per row.
    fn layer_with(fields: &[&str], rows: &[Vec<f64>]) -> String {
        let mut layer = WbLayer::new("pts").with_geom_type(GeometryType::Point);
        for f in fields {
            layer.add_field(FieldDef::new(*f, FieldType::Float));
        }
        for (i, row) in rows.iter().enumerate() {
            let n = rows.len();
            // Two coprime-to-n multipliers scatter points so geometric
            // neighbours are not the same as index neighbours.
            let px = ((i * 13) % n) as f64;
            let py = ((i * 29) % n) as f64;
            let attrs: Vec<(&str, FV)> = fields
                .iter()
                .zip(row.iter())
                .map(|(&f, &v)| (f, FV::Float(v)))
                .collect();
            layer
                .add_feature(Some(Geometry::point(px, py)), &attrs)
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    // The report ("output") is a CSV; parse it into `field -> value` rows
    // keyed by header name (the `variables` column never contains a comma —
    // fields are joined with `|` — so a plain split is safe).
    fn parse_csv(text: &str) -> Vec<BTreeMap<String, String>> {
        let mut lines = text.lines();
        let header: Vec<String> = lines
            .next()
            .unwrap()
            .split(',')
            .map(str::to_string)
            .collect();
        lines
            .filter(|l| !l.is_empty())
            .map(|l| {
                header
                    .iter()
                    .cloned()
                    .zip(l.split(',').map(str::to_string))
                    .collect()
            })
            .collect()
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Vec<BTreeMap<String, String>>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ExploratoryRegressionTool.run(&args, &ctx()).unwrap();
        let path = out.outputs["output"].as_str().unwrap();
        let rows = parse_csv(&std::fs::read_to_string(path).unwrap());
        (out, rows)
    }

    fn tmp_output(name: &str) -> String {
        std::env::temp_dir()
            .join(format!(
                "exploratory_regression_{name}_{}.csv",
                std::process::id()
            ))
            .to_str()
            .unwrap()
            .to_string()
    }

    /// Deterministic xorshift64* generator — no wall-clock/unseeded RNG, per
    /// repo policy. Only used to synthesize near-normal test noise (a raw
    /// sinusoid's marginal distribution is far from normal and trips the
    /// Jarque-Bera screen even for a tiny, well-fit perturbation).
    fn xorshift(seed: &mut u64) -> f64 {
        *seed ^= *seed << 13;
        *seed ^= *seed >> 7;
        *seed ^= *seed << 17;
        (*seed >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Approximately N(0,1) via the sum of 12 uniforms (classic Irwin-Hall /
    /// CLT construction), still fully deterministic.
    fn pseudo_normal(seed: &mut u64) -> f64 {
        (0..12).map(|_| xorshift(seed)).sum::<f64>() - 6.0
    }

    fn find_model<'a>(passing: &'a [Value], vars: &[&str]) -> Option<&'a Value> {
        passing.iter().find(|m| {
            let got: Vec<&str> = m["variables"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_str().unwrap())
                .collect();
            got == vars
        })
    }

    /// y = 2*x1 + small noise; x2 is unrelated. The {x1} model should pass
    /// every threshold and its adjusted R2 should beat the {x2}-only model.
    #[test]
    fn x1_model_passes_and_beats_x2() {
        let n = 60;
        let mut seed = 0x9E3779B97F4A7C15_u64;
        let mut rows = Vec::with_capacity(n);
        for i in 0..n {
            let x1 = (i as f64 * 0.37).sin() * 5.0 + i as f64 * 0.1;
            let x2 = (i as f64 * 0.19).cos() * 3.0;
            let noise = 0.3 * pseudo_normal(&mut seed);
            let y = 2.0 * x1 + noise;
            rows.push(vec![x1, x2, y]);
        }
        let input = layer_with(&["x1", "x2", "yv"], &rows);
        let (out, table) = run(json!({
            "input": input, "dependent_field": "yv", "explanatory_fields": "x1,x2",
            "min_vars": 1, "max_vars": 2, "neighbors": 6, "output": tmp_output("x1_beats_x2"),
        }));

        // No silent truncation on a tiny search: every combination evaluated.
        assert_eq!(out.outputs["models_possible"], json!("3"));
        assert_eq!(out.outputs["models_evaluated"], json!(3));
        assert_eq!(out.outputs["models_dropped"], json!("0"));
        assert_eq!(table.len(), 3);

        let passing = out.outputs["passing_models"].as_array().unwrap();
        let x1_only = find_model(passing, &["x1"]).expect("{x1} model should pass");
        let x1_adj_r2 = x1_only["adj_r2"].as_f64().unwrap();
        assert!(x1_adj_r2 > 0.9, "adj_r2 {x1_adj_r2}");
        assert!(x1_only["max_coef_p"].as_f64().unwrap() <= 0.05);

        // Whether or not {x2} alone passes, x1 must explain far more variance.
        let x2_row = table.iter().find(|r| r["variables"] == "x2").unwrap();
        let x1_row = table.iter().find(|r| r["variables"] == "x1").unwrap();
        let x2_adj_r2: f64 = x2_row["adj_r2"].parse().unwrap();
        let x1_adj_r2_csv: f64 = x1_row["adj_r2"].parse().unwrap();
        assert!(
            (x1_adj_r2_csv - x1_adj_r2).abs() < 1e-6,
            "CSV and JSON adj_r2 for {{x1}} should agree"
        );
        assert!(
            x1_adj_r2 > x2_adj_r2,
            "x1 adj_r2 {x1_adj_r2} should beat x2 adj_r2 {x2_adj_r2}"
        );

        // Variable summary: x1 significant every time it appears, x2 less so.
        let summary = out.outputs["variable_summary"].as_array().unwrap();
        let x1_sum = summary.iter().find(|v| v["field"] == "x1").unwrap();
        assert_eq!(x1_sum["times_evaluated"], json!(2));
        assert!(x1_sum["times_significant"].as_u64().unwrap() >= 1);
        assert_eq!(x1_sum["dominant_sign"], json!("positive"));
    }

    /// Two near-perfectly collinear explanatory fields trip the VIF screen
    /// when combined, even though each alone (or the pair) predicts y well.
    #[test]
    fn collinear_pair_trips_vif() {
        let n = 40;
        let mut rows = Vec::with_capacity(n);
        for i in 0..n {
            let xa = (i as f64 * 0.31).sin() * 3.0 + i as f64 * 0.05;
            let xb = xa * 2.0 + 0.001 * (i as f64 * 0.7).cos(); // ~ perfectly collinear
            let y = 1.0 + xa + 0.5 * xb;
            rows.push(vec![xa, xb, y]);
        }
        let input = layer_with(&["xa", "xb", "yv"], &rows);
        let (out, table) = run(json!({
            "input": input, "dependent_field": "yv", "explanatory_fields": "xa,xb",
            "min_vars": 1, "max_vars": 2, "neighbors": 5, "output": tmp_output("vif"),
        }));
        assert_eq!(out.outputs["models_evaluated"], json!(3));

        let combo_row = table.iter().find(|r| r["variables"] == "xa|xb").unwrap();
        let max_vif: f64 = combo_row["max_vif"].parse().unwrap();
        assert!(
            max_vif > 7.5,
            "max_vif {max_vif} should exceed the default threshold"
        );
        assert_eq!(combo_row["passed"], "false");

        // It also shows up in the passing_models list as absent.
        let passing = out.outputs["passing_models"].as_array().unwrap();
        assert!(find_model(passing, &["xa", "xb"]).is_none());
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = ExploratoryRegressionTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "p.geojson", "dependent_field": "yv" })).is_err());
        // dependent_field also listed as a candidate.
        assert!(bad(json!({
            "input": "p.geojson", "dependent_field": "yv", "explanatory_fields": "yv,x1",
            "output": "out.csv"
        }))
        .is_err());
        // min_vars > max_vars.
        assert!(bad(json!({
            "input": "p.geojson", "dependent_field": "yv", "explanatory_fields": "x1,x2",
            "min_vars": 2, "max_vars": 1, "output": "out.csv"
        }))
        .is_err());
        // max_vars exceeds candidate count.
        assert!(bad(json!({
            "input": "p.geojson", "dependent_field": "yv", "explanatory_fields": "x1",
            "max_vars": 2, "output": "out.csv"
        }))
        .is_err());
        // missing required output.
        assert!(bad(json!({
            "input": "p.geojson", "dependent_field": "yv", "explanatory_fields": "x1,x2"
        }))
        .is_err());
        assert!(bad(json!({
            "input": "p.geojson", "dependent_field": "yv", "explanatory_fields": "x1,x2",
            "output": "out.csv"
        }))
        .is_ok());
    }
}
