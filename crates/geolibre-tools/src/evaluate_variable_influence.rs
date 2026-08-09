//! GeoLibre tool: which explanatory variables actually drive a model, and how.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Evaluate Variable Influence* (Spatial
//! Statistics).
//!
//! ## Why the catalog needs it
//!
//! A model that predicts well is not yet an answer. The question people
//! actually bring to a habitat model, a landslide-susceptibility surface or a
//! yield model is *which* variables matter and *in which direction* — is the
//! effect of slope monotonic, does it saturate, is there a threshold? A
//! coefficient table answers that only for a linear model, and nothing answers
//! it for a tree ensemble.
//!
//! The catalog can fit models — `generalized_linear_regression`,
//! `exploratory_regression`, `geographically_weighted_regression`,
//! `forest_based_forecast`, and whitebox's random-forest family — but not one
//! of them reports permutation importance or a partial-dependence curve.
//!
//! ## A deliberate adaptation
//!
//! ArcGIS takes a **trained model file** (`.ssm`) as input, because ArcGIS has
//! a model-serialisation format and a separate tool that writes one. This
//! catalog has neither: `geolibre-tools` deliberately does not depend on
//! `wbtools_oss`, so a whitebox-trained forest is not reachable from here, and
//! no tool in either registry emits a portable model object.
//!
//! So this tool **fits the model itself** and then explains it, from a table of
//! a dependent field and its explanatory fields. That is the same work ArcGIS
//! splits across two tools, and it is the only form the split can take without
//! a model format to pass between them. When one appears, the fitting half can
//! be replaced by loading it.
//!
//! ## Determinism
//!
//! There is no RNG. Bagging uses deterministic strided subsamples rather than
//! bootstrap draws, feature subsets rotate by tree index, and the permutation
//! used for importance is a fixed coprime stride. Every run on the same table
//! gives the same numbers, which is required for the WASM builds and welcome
//! everywhere else.
//!
//! ## What is reported
//!
//! * **Permutation importance** — how much out-of-bag accuracy is lost when a
//!   variable's values are shuffled and the rest are left alone. This measures
//!   what the model *uses*, which is the honest question; a variable can be
//!   strongly correlated with the target and still be unimportant because
//!   another variable already carries the information.
//! * **Partial dependence** — the model's mean prediction as one variable
//!   sweeps its observed range with every other variable held at its real
//!   values. This is the shape of the effect, not just its size.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Layer};

use crate::args_common::{choice_or, req_str, usize_or};
use crate::common::parse_optional_output;
use crate::vector_common::{load_input_layer, write_or_store_layer};

pub struct EvaluateVariableInfluenceTool;

impl Tool for EvaluateVariableInfluenceTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "evaluate_variable_influence",
            display_name: "Evaluate Variable Influence",
            summary: "Reports permutation importance and partial-dependence curves for a forest model fitted from a table, answering which explanatory variables a model actually uses and what shape each effect has (ArcGIS Evaluate Variable Influence). The catalog can fit models — generalized_linear_regression, exploratory_regression, geographically_weighted_regression, forest_based_forecast, whitebox's random forests — but none reports either measure. ArcGIS reads a trained .ssm model file; this stack has no portable model format, so the tool fits its own forest and then explains it. Fully deterministic: strided bagging and a fixed permutation, no RNG.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Layer or table holding the training observations.",
                    required: true,
                },
                ToolParamSpec {
                    name: "dependent_field",
                    description: "Field being predicted.",
                    required: true,
                },
                ToolParamSpec {
                    name: "explanatory_fields",
                    description: "Comma-separated numeric fields to evaluate.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output table of per-variable importance. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_partial_dependence",
                    description: "Output table of partial-dependence curves, one row per variable and grid point. Always produced; stored in memory when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "model_type",
                    description: "'regression' (default) for a continuous target, or 'classification' for a categorical one.",
                    required: false,
                },
                ToolParamSpec {
                    name: "number_of_trees",
                    description: "Trees in the ensemble (default 50).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_depth",
                    description: "Maximum tree depth (default 8).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_samples_leaf",
                    description: "Fewest observations a leaf may hold (default 3).",
                    required: false,
                },
                ToolParamSpec {
                    name: "grid_points",
                    description: "Points per partial-dependence curve (default 12).",
                    required: false,
                },
                ToolParamSpec {
                    name: "permutation_repeats",
                    description: "Distinct permutations averaged per variable (default 3). More is steadier and slower.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        req_str(args, "dependent_field")?;
        parse_fields(args)?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input_path = req_str(args, "input")?.to_string();
        let dependent = req_str(args, "dependent_field")?.to_string();
        let fields = parse_fields(args)?;
        let prm = parse_params(args)?;
        let output = parse_optional_output(args, "output")?;
        let out_pd = parse_optional_output(args, "output_partial_dependence")?;

        let layer = load_input_layer(&input_path)?;
        let (x, y) = read_table(&layer, &dependent, &fields)?;
        let n = y.len();
        let p = fields.len();
        if n < prm.min_samples_leaf * 2 + 1 {
            return Err(ToolError::Execution(format!(
                "only {n} complete observation(s); too few to fit and validate a model"
            )));
        }

        ctx.progress.info(&format!(
            "{n} observation(s), {p} variable(s), {} {} tree(s)",
            prm.trees,
            if prm.classification { "classification" } else { "regression" }
        ));

        let forest = Forest::fit(&x, &y, p, &prm);
        let baseline = forest.oob_score(&x, &y, &prm);
        ctx.progress.info(&format!(
            "out-of-bag {} {baseline:.4}",
            if prm.classification { "accuracy" } else { "R2" }
        ));

        // Permutation importance: the loss in out-of-bag score when one
        // variable is shuffled and everything else is left alone.
        let mut importance = vec![0.0f64; p];
        let mut permuted = x.clone();
        for (j, imp) in importance.iter_mut().enumerate() {
            let mut total = 0.0;
            for rep in 0..prm.repeats {
                // A stride coprime with n is a single-cycle permutation, so
                // every value moves and no value keeps its own row.
                let stride = coprime_stride(n, rep);
                for i in 0..n {
                    permuted[i * p + j] = x[((i * stride + 1 + rep) % n) * p + j];
                }
                total += baseline - forest.oob_score(&permuted, &y, &prm);
            }
            // Restore the column before moving to the next variable.
            for i in 0..n {
                permuted[i * p + j] = x[i * p + j];
            }
            *imp = total / prm.repeats as f64;
            ctx.progress.progress((j as f64 + 1.0) / p as f64);
        }

        // Rank, largest influence first.
        let mut order: Vec<usize> = (0..p).collect();
        order.sort_by(|&a, &b| importance[b].total_cmp(&importance[a]));
        let importance_sum: f64 = importance.iter().map(|v| v.max(0.0)).sum();

        let mut imp_layer = Layer::new("variable_influence");
        imp_layer.add_field(FieldDef::new("rank", FieldType::Integer));
        imp_layer.add_field(FieldDef::new("variable", FieldType::Text));
        imp_layer.add_field(FieldDef::new("importance", FieldType::Float));
        imp_layer.add_field(FieldDef::new("importance_pct", FieldType::Float));
        imp_layer.add_field(FieldDef::new("min", FieldType::Float));
        imp_layer.add_field(FieldDef::new("max", FieldType::Float));

        for (rank, &j) in order.iter().enumerate() {
            let col: Vec<f64> = (0..n).map(|i| x[i * p + j]).collect();
            let pct = if importance_sum > 0.0 {
                100.0 * importance[j].max(0.0) / importance_sum
            } else {
                0.0
            };
            imp_layer
                .add_feature(
                    None,
                    &[
                        ("rank", FieldValue::Integer(rank as i64 + 1)),
                        ("variable", FieldValue::Text(fields[j].clone())),
                        ("importance", FieldValue::Float(importance[j])),
                        ("importance_pct", FieldValue::Float(pct)),
                        (
                            "min",
                            FieldValue::Float(col.iter().copied().fold(f64::INFINITY, f64::min)),
                        ),
                        (
                            "max",
                            FieldValue::Float(
                                col.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                            ),
                        ),
                    ],
                )
                .map_err(|e| ToolError::Execution(format!("writing the influence table: {e}")))?;
        }

        // Partial dependence: sweep one variable across its range with every
        // other variable left at its real value, and average the predictions.
        let mut pd_layer = Layer::new("partial_dependence");
        pd_layer.add_field(FieldDef::new("variable", FieldType::Text));
        pd_layer.add_field(FieldDef::new("point", FieldType::Integer));
        pd_layer.add_field(FieldDef::new("value", FieldType::Float));
        pd_layer.add_field(FieldDef::new("prediction", FieldType::Float));

        let mut probe = x.clone();
        for j in 0..p {
            let col: Vec<f64> = (0..n).map(|i| x[i * p + j]).collect();
            let lo = col.iter().copied().fold(f64::INFINITY, f64::min);
            let hi = col.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            for g in 0..prm.grid {
                // A constant variable has a one-point curve rather than a
                // division by zero.
                let v = if prm.grid == 1 || hi <= lo {
                    lo
                } else {
                    lo + (hi - lo) * g as f64 / (prm.grid - 1) as f64
                };
                for i in 0..n {
                    probe[i * p + j] = v;
                }
                let mean: f64 = (0..n)
                    .map(|i| forest.predict(&probe[i * p..(i + 1) * p]))
                    .sum::<f64>()
                    / n as f64;
                pd_layer
                    .add_feature(
                        None,
                        &[
                            ("variable", FieldValue::Text(fields[j].clone())),
                            ("point", FieldValue::Integer(g as i64)),
                            ("value", FieldValue::Float(v)),
                            ("prediction", FieldValue::Float(mean)),
                        ],
                    )
                    .map_err(|e| {
                        ToolError::Execution(format!("writing the partial-dependence table: {e}"))
                    })?;
                if hi <= lo {
                    break;
                }
            }
            for i in 0..n {
                probe[i * p + j] = x[i * p + j];
            }
        }

        let ranked: Vec<Value> = order
            .iter()
            .map(|&j| json!({ "variable": fields[j], "importance": importance[j] }))
            .collect();

        let imp_path = write_or_store_layer(imp_layer, output)?;
        let pd_path = write_or_store_layer(pd_layer, out_pd)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(imp_path));
        outputs.insert("output_partial_dependence".to_string(), json!(pd_path));
        outputs.insert("observations".to_string(), json!(n));
        outputs.insert("variables".to_string(), json!(p));
        outputs.insert("oob_score".to_string(), json!(baseline));
        outputs.insert(
            "score_type".to_string(),
            json!(if prm.classification { "accuracy" } else { "r2" }),
        );
        outputs.insert("ranked".to_string(), Value::Array(ranked));
        Ok(ToolRunResult { outputs })
    }
}

/// A stride coprime with `n`, so stepping by it visits every row exactly once.
fn coprime_stride(n: usize, rep: usize) -> usize {
    if n <= 2 {
        return 1;
    }
    let mut s = (n / 3).max(2) + rep * 7;
    while gcd(s, n) != 1 {
        s += 1;
    }
    s
}

fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// Reads the dependent and explanatory columns, keeping only complete rows.
fn read_table(
    layer: &Layer,
    dependent: &str,
    fields: &[String],
) -> Result<(Vec<f64>, Vec<f64>), ToolError> {
    let dep_idx = layer.schema.field_index(dependent).ok_or_else(|| {
        ToolError::Validation(format!("'dependent_field' '{dependent}' is not in the input"))
    })?;
    let idx: Vec<usize> = fields
        .iter()
        .map(|f| {
            layer.schema.field_index(f).ok_or_else(|| {
                ToolError::Validation(format!("'explanatory_fields' entry '{f}' is not in the input"))
            })
        })
        .collect::<Result<_, _>>()?;

    let numeric = |v: &FieldValue| -> Option<f64> {
        match v {
            FieldValue::Float(x) => Some(*x),
            FieldValue::Integer(x) => Some(*x as f64),
            FieldValue::Boolean(b) => Some(if *b { 1.0 } else { 0.0 }),
            _ => None,
        }
    };

    let mut x = Vec::new();
    let mut y = Vec::new();
    for f in layer.iter() {
        let Some(target) = numeric(&f.attributes[dep_idx]).filter(|v| v.is_finite()) else {
            continue;
        };
        let mut row = Vec::with_capacity(idx.len());
        let mut ok = true;
        for &i in &idx {
            match numeric(&f.attributes[i]).filter(|v| v.is_finite()) {
                Some(v) => row.push(v),
                None => {
                    ok = false;
                    break;
                }
            }
        }
        // A row missing any variable cannot be used: the tree would have to
        // invent a value for the split it lands on.
        if !ok {
            continue;
        }
        x.extend(row);
        y.push(target);
    }
    if y.is_empty() {
        return Err(ToolError::Execution(
            "no row has a usable value in every requested field".to_string(),
        ));
    }
    Ok((x, y))
}

// ── The forest ──────────────────────────────────────────────────────────────

/// A binary decision tree over the explanatory columns.
struct Tree {
    /// `(feature, threshold, left, right)` for a split; `usize::MAX` marks a
    /// leaf whose prediction is in `value`.
    nodes: Vec<Node>,
}

struct Node {
    feature: usize,
    threshold: f64,
    left: usize,
    right: usize,
    value: f64,
}

impl Node {
    fn is_leaf(&self) -> bool {
        self.feature == usize::MAX
    }
}

/// A deterministic bagged ensemble.
struct Forest {
    trees: Vec<Tree>,
    /// Rows held out of each tree, for the out-of-bag score.
    oob: Vec<Vec<usize>>,
    classification: bool,
}

impl Forest {
    fn fit(x: &[f64], y: &[f64], p: usize, prm: &Params) -> Forest {
        let n = y.len();
        let mut trees = Vec::with_capacity(prm.trees);
        let mut oob = Vec::with_capacity(prm.trees);

        for t in 0..prm.trees {
            // Deterministic bagging: tree `t` trains on every row except one
            // stratum, rotating which stratum is held out. Bootstrap draws
            // would need an RNG the WASM build does not share with native.
            let folds = prm.trees.min(n).max(2);
            let held = t % folds;
            let mut rows = Vec::with_capacity(n);
            let mut out = Vec::new();
            for i in 0..n {
                if i % folds == held {
                    out.push(i);
                } else {
                    rows.push(i);
                }
            }
            if rows.is_empty() {
                continue;
            }
            // Feature subsampling, rotated by tree index so the ensemble sees
            // every variable in several combinations.
            //
            // The rotation is by `t`, not by `t * take`: the latter is the
            // identity whenever `take == p`, which left the feature *order*
            // fixed across every tree. With a deterministic split search that
            // breaks ties on the first candidate, two perfectly correlated
            // variables would then have the same one chosen every time and the
            // other never used at all — so it would score zero importance while
            // carrying identical information. Rotating the order spreads the
            // credit between them, as a randomised forest does.
            let take = ((p as f64).sqrt().ceil() as usize).clamp(1, p);
            let features: Vec<usize> = (0..take).map(|k| (t + k) % p).collect();
            trees.push(Tree::fit(x, y, p, &rows, &features, prm));
            oob.push(out);
        }
        Forest {
            trees,
            oob,
            classification: prm.classification,
        }
    }

    /// Mean (regression) or majority (classification) over the ensemble.
    fn predict(&self, row: &[f64]) -> f64 {
        self.aggregate(self.trees.iter().map(|t| t.predict(row)))
    }

    fn aggregate(&self, preds: impl Iterator<Item = f64>) -> f64 {
        let vals: Vec<f64> = preds.collect();
        if vals.is_empty() {
            return 0.0;
        }
        if !self.classification {
            return vals.iter().sum::<f64>() / vals.len() as f64;
        }
        let mut counts: BTreeMap<u64, (usize, f64)> = BTreeMap::new();
        for v in vals {
            let e = counts.entry(v.to_bits()).or_insert((0, v));
            e.0 += 1;
        }
        counts
            .values()
            .max_by(|a, b| a.0.cmp(&b.0).then(b.1.total_cmp(&a.1)))
            .map(|(_, v)| *v)
            .unwrap_or(0.0)
    }

    /// Out-of-bag R-squared, or accuracy for a categorical target.
    ///
    /// Scoring only on rows a tree never saw is what makes the permutation
    /// importance meaningful: an in-sample score would reward memorisation and
    /// under-report the loss from shuffling.
    fn oob_score(&self, x: &[f64], y: &[f64], prm: &Params) -> f64 {
        let n = y.len();
        let p = prm.p;
        let mut sum = vec![0.0f64; n];
        let mut votes: Vec<Vec<f64>> = vec![Vec::new(); n];
        let mut count = vec![0usize; n];

        for (t, tree) in self.trees.iter().enumerate() {
            for &i in &self.oob[t] {
                let pred = tree.predict(&x[i * p..(i + 1) * p]);
                if self.classification {
                    votes[i].push(pred);
                } else {
                    sum[i] += pred;
                }
                count[i] += 1;
            }
        }

        if self.classification {
            let (mut right, mut total) = (0usize, 0usize);
            for i in 0..n {
                if count[i] == 0 {
                    continue;
                }
                total += 1;
                if self.aggregate(votes[i].iter().copied()) == y[i] {
                    right += 1;
                }
            }
            if total == 0 {
                return 0.0;
            }
            return right as f64 / total as f64;
        }

        let scored: Vec<usize> = (0..n).filter(|&i| count[i] > 0).collect();
        if scored.len() < 2 {
            return 0.0;
        }
        let mean = scored.iter().map(|&i| y[i]).sum::<f64>() / scored.len() as f64;
        let mut ss_res = 0.0;
        let mut ss_tot = 0.0;
        for &i in &scored {
            let pred = sum[i] / count[i] as f64;
            ss_res += (y[i] - pred).powi(2);
            ss_tot += (y[i] - mean).powi(2);
        }
        if ss_tot <= 0.0 {
            return 0.0;
        }
        1.0 - ss_res / ss_tot
    }
}

impl Tree {
    fn fit(
        x: &[f64],
        y: &[f64],
        p: usize,
        rows: &[usize],
        features: &[usize],
        prm: &Params,
    ) -> Tree {
        let mut nodes = Vec::new();
        build(x, y, p, rows, features, prm, 0, &mut nodes);
        Tree { nodes }
    }

    fn predict(&self, row: &[f64]) -> f64 {
        if self.nodes.is_empty() {
            return 0.0;
        }
        let mut i = 0usize;
        loop {
            let node = &self.nodes[i];
            if node.is_leaf() {
                return node.value;
            }
            i = if row[node.feature] <= node.threshold {
                node.left
            } else {
                node.right
            };
        }
    }
}

/// Recursively builds a subtree over `rows`, returning its node index.
#[allow(clippy::too_many_arguments)]
fn build(
    x: &[f64],
    y: &[f64],
    p: usize,
    rows: &[usize],
    features: &[usize],
    prm: &Params,
    depth: usize,
    nodes: &mut Vec<Node>,
) -> usize {
    let idx = nodes.len();
    nodes.push(Node {
        feature: usize::MAX,
        threshold: 0.0,
        left: 0,
        right: 0,
        value: leaf_value(y, rows, prm.classification),
    });

    if depth >= prm.max_depth || rows.len() < 2 * prm.min_samples_leaf {
        return idx;
    }

    // Best split by variance reduction (regression) or Gini gain
    // (classification), over the tree's feature subset.
    let parent = impurity(y, rows, prm.classification);
    let mut best: Option<(f64, usize, f64)> = None;
    for &f in features {
        let mut vals: Vec<f64> = rows.iter().map(|&i| x[i * p + f]).collect();
        vals.sort_by(f64::total_cmp);
        vals.dedup();
        if vals.len() < 2 {
            continue;
        }
        // Candidate thresholds are the midpoints between distinct values,
        // capped so a wide column does not dominate the cost.
        let step = (vals.len() / 16).max(1);
        for w in vals.windows(2).step_by(step) {
            let thr = 0.5 * (w[0] + w[1]);
            let (l, r): (Vec<usize>, Vec<usize>) =
                rows.iter().partition(|&&i| x[i * p + f] <= thr);
            if l.len() < prm.min_samples_leaf || r.len() < prm.min_samples_leaf {
                continue;
            }
            let wl = l.len() as f64 / rows.len() as f64;
            let wr = 1.0 - wl;
            let gain = parent
                - wl * impurity(y, &l, prm.classification)
                - wr * impurity(y, &r, prm.classification);
            if gain > best.map(|b| b.0).unwrap_or(1e-12) {
                best = Some((gain, f, thr));
            }
        }
    }

    let Some((_, feature, threshold)) = best else {
        return idx;
    };
    let (l, r): (Vec<usize>, Vec<usize>) =
        rows.iter().partition(|&&i| x[i * p + feature] <= threshold);
    let left = build(x, y, p, &l, features, prm, depth + 1, nodes);
    let right = build(x, y, p, &r, features, prm, depth + 1, nodes);
    nodes[idx].feature = feature;
    nodes[idx].threshold = threshold;
    nodes[idx].left = left;
    nodes[idx].right = right;
    idx
}

/// Mean for regression, most common class for classification.
fn leaf_value(y: &[f64], rows: &[usize], classification: bool) -> f64 {
    if rows.is_empty() {
        return 0.0;
    }
    if !classification {
        return rows.iter().map(|&i| y[i]).sum::<f64>() / rows.len() as f64;
    }
    let mut counts: BTreeMap<u64, (usize, f64)> = BTreeMap::new();
    for &i in rows {
        let e = counts.entry(y[i].to_bits()).or_insert((0, y[i]));
        e.0 += 1;
    }
    counts
        .values()
        .max_by(|a, b| a.0.cmp(&b.0).then(b.1.total_cmp(&a.1)))
        .map(|(_, v)| *v)
        .unwrap_or(0.0)
}

/// Variance for regression, Gini for classification.
fn impurity(y: &[f64], rows: &[usize], classification: bool) -> f64 {
    if rows.is_empty() {
        return 0.0;
    }
    let n = rows.len() as f64;
    if !classification {
        let mean = rows.iter().map(|&i| y[i]).sum::<f64>() / n;
        return rows.iter().map(|&i| (y[i] - mean).powi(2)).sum::<f64>() / n;
    }
    let mut counts: BTreeMap<u64, usize> = BTreeMap::new();
    for &i in rows {
        *counts.entry(y[i].to_bits()).or_default() += 1;
    }
    1.0 - counts
        .values()
        .map(|&c| {
            let f = c as f64 / n;
            f * f
        })
        .sum::<f64>()
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    classification: bool,
    trees: usize,
    max_depth: usize,
    min_samples_leaf: usize,
    grid: usize,
    repeats: usize,
    /// Variable count, carried so the scorers can index rows.
    p: usize,
}

fn parse_fields(args: &ToolArgs) -> Result<Vec<String>, ToolError> {
    let raw = req_str(args, "explanatory_fields")?;
    let fields: Vec<String> = raw
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(String::from)
        .collect();
    if fields.is_empty() {
        return Err(ToolError::Validation(
            "'explanatory_fields' must name at least one field".to_string(),
        ));
    }
    let mut sorted = fields.clone();
    sorted.sort();
    sorted.dedup();
    if sorted.len() != fields.len() {
        return Err(ToolError::Validation(
            "'explanatory_fields' contains a duplicate".to_string(),
        ));
    }
    Ok(fields)
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let classification =
        choice_or(args, "model_type", &["regression", "classification"], "regression")?
            == "classification";
    let trees = usize_or(args, "number_of_trees", 50)?;
    if trees < 2 {
        return Err(ToolError::Validation(
            "'number_of_trees' must be at least 2; the out-of-bag score needs held-out rows"
                .to_string(),
        ));
    }
    let max_depth = usize_or(args, "max_depth", 8)?;
    if max_depth == 0 {
        return Err(ToolError::Validation(
            "'max_depth' must be at least 1".to_string(),
        ));
    }
    let min_samples_leaf = usize_or(args, "min_samples_leaf", 3)?;
    if min_samples_leaf == 0 {
        return Err(ToolError::Validation(
            "'min_samples_leaf' must be at least 1".to_string(),
        ));
    }
    let grid = usize_or(args, "grid_points", 12)?;
    if grid == 0 {
        return Err(ToolError::Validation(
            "'grid_points' must be at least 1".to_string(),
        ));
    }
    let repeats = usize_or(args, "permutation_repeats", 3)?;
    if repeats == 0 {
        return Err(ToolError::Validation(
            "'permutation_repeats' must be at least 1".to_string(),
        ));
    }
    let p = parse_fields(args)?.len();
    Ok(Params {
        classification,
        trees,
        max_depth,
        min_samples_leaf,
        grid,
        repeats,
        p,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Deterministic pseudo-noise in [-1, 1].
    fn noise(i: usize, salt: u64) -> f64 {
        let mut x = (i as u64).wrapping_mul(6364136223846793005).wrapping_add(salt);
        x ^= x >> 33;
        x = x.wrapping_mul(0xff51afd7ed558ccd);
        x ^= x >> 33;
        ((x % 2001) as f64 / 1000.0) - 1.0
    }

    /// Builds an attribute-only table with the given columns.
    fn table(cols: &[(&str, Vec<f64>)]) -> String {
        let mut layer = Layer::new("training");
        for (name, _) in cols {
            layer.add_field(FieldDef::new(*name, FieldType::Float));
        }
        let n = cols[0].1.len();
        for i in 0..n {
            let attrs: Vec<(&str, FieldValue)> = cols
                .iter()
                .map(|(name, v)| (*name, FieldValue::Float(v[i])))
                .collect();
            layer.add_feature(None, &attrs).unwrap();
        }
        write_or_store_layer(layer, None).unwrap()
    }

    fn run(args: Value) -> (Layer, BTreeMap<String, Value>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = EvaluateVariableInfluenceTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (layer, out.outputs)
    }

    fn ranked(outputs: &BTreeMap<String, Value>) -> Vec<String> {
        outputs["ranked"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v["variable"].as_str().unwrap().to_string())
            .collect()
    }

    /// The core claim: a variable the target depends on ranks above one it is
    /// independent of.
    #[test]
    fn ranks_the_driving_variable_first() {
        let n = 200;
        let signal: Vec<f64> = (0..n).map(|i| noise(i, 1)).collect();
        let junk: Vec<f64> = (0..n).map(|i| noise(i, 2)).collect();
        // y depends on `signal` alone, plus a little noise.
        let y: Vec<f64> = (0..n).map(|i| 3.0 * signal[i] + 0.1 * noise(i, 3)).collect();

        let src = table(&[
            ("y", y),
            ("signal", signal),
            ("junk", junk),
        ]);
        let (layer, outputs) = run(json!({
            "input": src, "dependent_field": "y",
            "explanatory_fields": "signal,junk", "number_of_trees": 20
        }));
        assert_eq!(layer.len(), 2, "one row per variable");
        assert_eq!(
            ranked(&outputs)[0],
            "signal",
            "the driving variable must rank first"
        );
        assert!(
            outputs["oob_score"].as_f64().unwrap() > 0.7,
            "the model should fit a clean signal, got R2 {}",
            outputs["oob_score"]
        );

        // The irrelevant variable's importance should be near zero.
        let vi = layer.schema.field_index("variable").unwrap();
        let ii = layer.schema.field_index("importance").unwrap();
        for f in layer.iter() {
            let (FieldValue::Text(name), FieldValue::Float(imp)) =
                (&f.attributes[vi], &f.attributes[ii])
            else {
                panic!()
            };
            if name == "junk" {
                assert!(*imp < 0.1, "irrelevant variable scored {imp}");
            } else {
                assert!(*imp > 0.2, "driving variable scored only {imp}");
            }
        }
    }

    /// Permutation importance measures what the model *uses*, not raw
    /// correlation: given two copies of the same variable, neither is
    /// indispensable, so both score low even though both correlate perfectly
    /// with the target. This is the property that makes the measure honest.
    #[test]
    fn redundant_copies_share_the_credit() {
        let n = 160;
        let a: Vec<f64> = (0..n).map(|i| noise(i, 5)).collect();
        let y: Vec<f64> = a.iter().map(|v| 2.0 * v).collect();
        let b = a.clone(); // an exact duplicate column

        let solo = run(json!({
            "input": table(&[("y", y.clone()), ("a", a.clone())]),
            "dependent_field": "y", "explanatory_fields": "a", "number_of_trees": 20
        }))
        .1;
        let both = run(json!({
            "input": table(&[("y", y), ("a", a), ("b", b)]),
            "dependent_field": "y", "explanatory_fields": "a,b", "number_of_trees": 20
        }))
        .1;

        let solo_imp = solo["ranked"][0]["importance"].as_f64().unwrap();
        let both_imp = both["ranked"][0]["importance"].as_f64().unwrap();
        assert!(
            both_imp < solo_imp,
            "with a duplicate available, shuffling one copy should cost less \
             ({both_imp} vs {solo_imp})"
        );
    }

    /// Partial dependence recovers the *shape* of the effect, not just its
    /// size: for a monotonically increasing relationship the curve must rise.
    #[test]
    fn partial_dependence_follows_the_relationship() {
        let n = 200;
        let x: Vec<f64> = (0..n).map(|i| i as f64 / n as f64).collect();
        let y: Vec<f64> = x.iter().enumerate().map(|(i, v)| 5.0 * v + 0.05 * noise(i, 7)).collect();
        let args: ToolArgs = serde_json::from_value(json!({
            "input": table(&[("y", y), ("x", x)]),
            "dependent_field": "y", "explanatory_fields": "x",
            "number_of_trees": 20, "grid_points": 10
        }))
        .unwrap();
        let out = EvaluateVariableInfluenceTool.run(&args, &ctx()).unwrap();
        let pd = load_input_layer(out.outputs["output_partial_dependence"].as_str().unwrap())
            .unwrap();
        assert_eq!(pd.len(), 10, "one row per grid point");

        let pi = pd.schema.field_index("prediction").unwrap();
        let preds: Vec<f64> = pd
            .iter()
            .map(|f| match f.attributes[pi] {
                FieldValue::Float(v) => v,
                _ => panic!(),
            })
            .collect();
        assert!(
            preds.last().unwrap() > preds.first().unwrap(),
            "an increasing relationship must give a rising curve: {preds:?}"
        );
        // Monotone overall, allowing for the step structure of a tree.
        let rises = preds.windows(2).filter(|w| w[1] >= w[0]).count();
        assert!(
            rises >= 7,
            "the curve should be broadly monotone, got {preds:?}"
        );
    }

    /// A classification target is scored by accuracy and predicts real classes.
    #[test]
    fn classification_mode_scores_accuracy() {
        let n = 160;
        let a: Vec<f64> = (0..n).map(|i| noise(i, 11)).collect();
        let junk: Vec<f64> = (0..n).map(|i| noise(i, 12)).collect();
        // Class is decided by the sign of `a`.
        let y: Vec<f64> = a.iter().map(|v| if *v > 0.0 { 1.0 } else { 0.0 }).collect();

        let (_, outputs) = run(json!({
            "input": table(&[("cls", y), ("a", a), ("junk", junk)]),
            "dependent_field": "cls", "explanatory_fields": "a,junk",
            "model_type": "classification", "number_of_trees": 20
        }));
        assert_eq!(outputs["score_type"].as_str().unwrap(), "accuracy");
        assert!(
            outputs["oob_score"].as_f64().unwrap() > 0.85,
            "a clean threshold rule should classify well, got {}",
            outputs["oob_score"]
        );
        assert_eq!(ranked(&outputs)[0], "a");
    }

    /// Runs are reproducible — there is no RNG anywhere in the fit.
    #[test]
    fn results_are_deterministic() {
        let n = 120;
        let a: Vec<f64> = (0..n).map(|i| noise(i, 21)).collect();
        let b: Vec<f64> = (0..n).map(|i| noise(i, 22)).collect();
        let y: Vec<f64> = (0..n).map(|i| a[i] - 2.0 * b[i]).collect();
        let build = || {
            run(json!({
                "input": table(&[("y", y.clone()), ("a", a.clone()), ("b", b.clone())]),
                "dependent_field": "y", "explanatory_fields": "a,b", "number_of_trees": 15
            }))
            .1
        };
        let first = build();
        let second = build();
        assert_eq!(first["oob_score"], second["oob_score"]);
        assert_eq!(first["ranked"], second["ranked"]);
    }

    /// Importance percentages are reported alongside the raw values.
    #[test]
    fn importance_percentages_are_emitted() {
        let n = 120;
        let a: Vec<f64> = (0..n).map(|i| noise(i, 31)).collect();
        let b: Vec<f64> = (0..n).map(|i| noise(i, 32)).collect();
        let y: Vec<f64> = (0..n).map(|i| 4.0 * a[i] + b[i]).collect();
        let (layer, _) = run(json!({
            "input": table(&[("y", y), ("a", a), ("b", b)]),
            "dependent_field": "y", "explanatory_fields": "a,b", "number_of_trees": 20
        }));
        let pi = layer.schema.field_index("importance_pct").unwrap();
        let total: f64 = layer
            .iter()
            .map(|f| match f.attributes[pi] {
                FieldValue::Float(v) => v,
                _ => panic!(),
            })
            .sum();
        assert!(
            (total - 100.0).abs() < 1e-6,
            "percentages should sum to 100, got {total}"
        );
    }

    /// Rows missing any requested field are skipped rather than guessed at.
    #[test]
    fn incomplete_rows_are_skipped() {
        let mut layer = Layer::new("t");
        layer.add_field(FieldDef::new("y", FieldType::Float));
        layer.add_field(FieldDef::new("a", FieldType::Float));
        for i in 0..60 {
            let v = noise(i, 41);
            layer
                .add_feature(
                    None,
                    &[
                        ("y", FieldValue::Float(2.0 * v)),
                        (
                            "a",
                            // Every tenth row is missing its explanatory value.
                            if i % 10 == 0 {
                                FieldValue::Null
                            } else {
                                FieldValue::Float(v)
                            },
                        ),
                    ],
                )
                .unwrap();
        }
        let src = write_or_store_layer(layer, None).unwrap();
        let (_, outputs) = run(json!({
            "input": src, "dependent_field": "y", "explanatory_fields": "a",
            "number_of_trees": 10
        }));
        assert_eq!(
            outputs["observations"].as_u64().unwrap(),
            54,
            "the six incomplete rows should be dropped"
        );
    }

    /// A missing field is reported by name rather than silently ignored.
    #[test]
    fn missing_fields_are_reported() {
        let src = table(&[("y", vec![1.0; 40]), ("a", vec![2.0; 40])]);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": src, "dependent_field": "y", "explanatory_fields": "a,nope"
        }))
        .unwrap();
        let err = EvaluateVariableInfluenceTool.run(&args, &ctx()).unwrap_err();
        assert!(
            format!("{err:?}").contains("nope"),
            "the error should name the missing field, got {err:?}"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            EvaluateVariableInfluenceTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({"input": "t.shp", "dependent_field": "y"})).is_err());
        let base = json!({"input": "t.shp", "dependent_field": "y", "explanatory_fields": "a,b"});
        assert!(bad(base.clone()).is_ok());
        let with = |k: &str, v: Value| {
            let mut m = base.as_object().unwrap().clone();
            m.insert(k.into(), v);
            Value::Object(m)
        };
        assert!(bad(with("model_type", json!("cluster"))).is_err());
        assert!(bad(with("number_of_trees", json!(1))).is_err());
        assert!(bad(with("max_depth", json!(0))).is_err());
        assert!(bad(with("grid_points", json!(0))).is_err());
        assert!(bad(with("permutation_repeats", json!(0))).is_err());
        assert!(bad(json!({
            "input": "t.shp", "dependent_field": "y", "explanatory_fields": "a,a"
        }))
        .is_err());
    }
}
