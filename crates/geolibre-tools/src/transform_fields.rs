//! GeoLibre tool: standardize, transform, and encode attribute fields.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Standardize Field* / *Transform
//! Field* / *Encode Field* (Data Management) tools. The repo's statistics
//! suite (`calculate_composite_index`, `dimension_reduction`,
//! `generalized_linear_regression`, …) each re-implement ad-hoc scaling
//! internally; this tool provides a single, named, documentable set of
//! attribute-prep transforms that can be chained ahead of them.
//!
//! Every transform is pure attribute math over the layer schema/features in
//! two passes: fit statistics over the non-null values of a field, then apply
//! the transform elementwise (null-aware — nulls and out-of-domain values stay
//! null rather than panicking or poisoning the fit). The same transform is
//! applied to every field named in `fields`; new fields are appended as
//! `<field><suffix>` (or, for `onehot`, one `<field><suffix><category>` field
//! per distinct category). The fitted parameters (mean/sd, median/IQR,
//! Box-Cox lambda, bin breaks, one-hot categories) are returned in the run
//! outputs so the transform is documentable and, where relevant, invertible.

use std::collections::{BTreeMap, HashMap, HashSet};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Feature, FieldDef, FieldType, FieldValue, Layer, Schema};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Cap on distinct categories materialized by the `onehot` transform; beyond
/// this the long tail is grouped into a single `<field><suffix>other` column.
const MAX_ONEHOT_CATEGORIES: usize = 50;

pub struct TransformFieldsTool;

impl Tool for TransformFieldsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "transform_fields",
            display_name: "Transform Fields",
            summary: "Standardize, transform, and encode attribute fields (like ArcGIS Standardize/Transform/Encode Field): z-score, min-max, robust, log/log1p/sqrt/boxcox/inverse, equal-interval/quantile/std-dev binning, and one-hot encoding.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (or table) with the fields to transform.",
                    required: true,
                },
                ToolParamSpec {
                    name: "fields",
                    description: "Comma-separated field names to transform; the same transform is applied to each.",
                    required: true,
                },
                ToolParamSpec {
                    name: "transform",
                    description: "Transform to apply: zscore, minmax, robust, log, log1p, sqrt, boxcox, inverse, bin, onehot.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector layer with the transformed fields added. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "bins",
                    description: "Number of classes for 'bin' (default 5).",
                    required: false,
                },
                ToolParamSpec {
                    name: "bin_method",
                    description: "Classification method for 'bin': equal_interval (default), quantile, or std_dev.",
                    required: false,
                },
                ToolParamSpec {
                    name: "boxcox_lambda",
                    description: "Fixed lambda for 'boxcox'; if omitted, fit by maximum likelihood (golden-section search).",
                    required: false,
                },
                ToolParamSpec {
                    name: "suffix",
                    description: "Output field name suffix (default depends on transform, e.g. '_z' for zscore).",
                    required: false,
                },
                ToolParamSpec {
                    name: "drop_input",
                    description: "If true, remove the original field(s) from the output after transforming (default false).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "fields")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let n = layer.len();

        let mut field_indices = Vec::with_capacity(prm.fields.len());
        for f in &prm.fields {
            let idx = layer
                .schema
                .field_index(f)
                .ok_or_else(|| ToolError::Validation(format!("field '{f}' not found")))?;
            field_indices.push(idx);
        }

        ctx.progress.info(&format!(
            "transform_fields: '{}' over {} field(s), {} feature(s)",
            prm.transform.as_str(),
            prm.fields.len(),
            n
        ));

        let suffix = prm
            .suffix
            .clone()
            .unwrap_or_else(|| prm.transform.default_suffix().to_string());

        let mut out_columns: Vec<OutputColumn> = Vec::new();
        let mut field_stats = serde_json::Map::new();

        for (field_name, &idx) in prm.fields.iter().zip(&field_indices) {
            match prm.transform {
                Transform::OneHot => {
                    let (cols, stats) = build_onehot_columns(&layer, idx, field_name, &suffix, ctx);
                    out_columns.extend(cols);
                    field_stats.insert(field_name.clone(), stats);
                }
                Transform::Bin => {
                    let raw = read_numeric_column(&layer, idx, n);
                    let (col, stats) =
                        build_bin_column(&raw, field_name, &suffix, prm.bins, prm.bin_method);
                    out_columns.push(col);
                    field_stats.insert(field_name.clone(), stats);
                }
                _ => {
                    let raw = read_numeric_column(&layer, idx, n);
                    let (col, stats) = build_scalar_column(
                        &raw,
                        field_name,
                        &suffix,
                        prm.transform,
                        prm.boxcox_lambda,
                    );
                    out_columns.push(col);
                    field_stats.insert(field_name.clone(), stats);
                }
            }
        }

        for col in &out_columns {
            layer.add_field(FieldDef::new(col.name.clone(), col.field_type));
        }
        for (fi, feat) in layer.features.iter_mut().enumerate() {
            for col in &out_columns {
                feat.attributes.push(col.values[fi].clone());
            }
        }

        if prm.drop_input {
            layer = drop_fields(layer, &prm.fields);
        }

        let output_field_names: Vec<String> = out_columns.iter().map(|c| c.name.clone()).collect();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("transform".to_string(), json!(prm.transform.as_str()));
        outputs.insert("fields".to_string(), json!(prm.fields));
        outputs.insert("output_fields".to_string(), json!(output_field_names));
        outputs.insert("field_stats".to_string(), Value::Object(field_stats));
        Ok(ToolRunResult { outputs })
    }
}

// ── Output column model ─────────────────────────────────────────────────────

struct OutputColumn {
    name: String,
    field_type: FieldType,
    values: Vec<FieldValue>,
}

fn read_numeric_column(layer: &Layer, idx: usize, n: usize) -> Vec<Option<f64>> {
    (0..n)
        .map(|fi| {
            layer.features[fi]
                .attributes
                .get(idx)
                .and_then(FieldValue::as_f64)
        })
        .collect()
}

/// Rebuilds `layer` with the named fields removed (the vector data model has
/// no in-place field-removal, so we copy the retained columns into a fresh
/// schema/feature set).
fn drop_fields(layer: Layer, drop_names: &[String]) -> Layer {
    let keep_idx: Vec<usize> = layer
        .schema
        .fields()
        .iter()
        .enumerate()
        .filter(|(_, f)| !drop_names.iter().any(|d| d == &f.name))
        .map(|(i, _)| i)
        .collect();
    let mut new_schema = Schema::new();
    for &i in &keep_idx {
        new_schema.add_field(layer.schema.fields()[i].clone());
    }
    let features: Vec<Feature> = layer
        .features
        .iter()
        .map(|f| Feature {
            fid: f.fid,
            geometry: f.geometry.clone(),
            attributes: keep_idx.iter().map(|&i| f.attributes[i].clone()).collect(),
        })
        .collect();
    Layer {
        name: layer.name,
        geom_type: layer.geom_type,
        crs: layer.crs,
        schema: new_schema,
        features,
        extent: layer.extent,
    }
}

// ── Scalar transforms (zscore/minmax/robust/log/log1p/sqrt/boxcox/inverse) ──

fn build_scalar_column(
    raw: &[Option<f64>],
    field_name: &str,
    suffix: &str,
    transform: Transform,
    fixed_lambda: Option<f64>,
) -> (OutputColumn, Value) {
    let vals: Vec<f64> = raw.iter().filter_map(|v| *v).collect();
    let name = format!("{field_name}{suffix}");

    let (values, stats): (Vec<Option<f64>>, Value) = match transform {
        Transform::ZScore => {
            let (mean, sd) = mean_sd(&vals);
            let values = raw
                .iter()
                .map(|v| v.map(|x| if sd > 0.0 { (x - mean) / sd } else { 0.0 }))
                .collect();
            (
                values,
                json!({ "transform": "zscore", "mean": mean, "sd": sd, "n": vals.len() }),
            )
        }
        Transform::MinMax => {
            let (lo, hi) = min_max(&vals);
            let range = hi - lo;
            let values = raw
                .iter()
                .map(|v| v.map(|x| if range > 0.0 { (x - lo) / range } else { 0.0 }))
                .collect();
            (
                values,
                json!({ "transform": "minmax", "min": lo, "max": hi }),
            )
        }
        Transform::Robust => {
            let mut sorted = vals.clone();
            sorted.sort_by(|a, b| a.total_cmp(b));
            let median = percentile(&sorted, 50.0);
            let q1 = percentile(&sorted, 25.0);
            let q3 = percentile(&sorted, 75.0);
            let iqr = q3 - q1;
            let values = raw
                .iter()
                .map(|v| v.map(|x| if iqr > 0.0 { (x - median) / iqr } else { 0.0 }))
                .collect();
            (
                values,
                json!({ "transform": "robust", "median": median, "q1": q1, "q3": q3, "iqr": iqr }),
            )
        }
        Transform::Log => {
            let mut n_guarded = 0usize;
            let values = raw
                .iter()
                .map(|v| match v {
                    Some(x) if *x > 0.0 => Some(x.ln()),
                    Some(_) => {
                        n_guarded += 1;
                        None
                    }
                    None => None,
                })
                .collect();
            (
                values,
                json!({ "transform": "log", "n_guarded_nonpositive": n_guarded }),
            )
        }
        Transform::Log1p => {
            let mut n_guarded = 0usize;
            let values = raw
                .iter()
                .map(|v| match v {
                    Some(x) if *x > -1.0 => Some((1.0 + x).ln()),
                    Some(_) => {
                        n_guarded += 1;
                        None
                    }
                    None => None,
                })
                .collect();
            (
                values,
                json!({ "transform": "log1p", "n_guarded_le_neg1": n_guarded }),
            )
        }
        Transform::Sqrt => {
            let mut n_guarded = 0usize;
            let values = raw
                .iter()
                .map(|v| match v {
                    Some(x) if *x >= 0.0 => Some(x.sqrt()),
                    Some(_) => {
                        n_guarded += 1;
                        None
                    }
                    None => None,
                })
                .collect();
            (
                values,
                json!({ "transform": "sqrt", "n_guarded_negative": n_guarded }),
            )
        }
        Transform::Inverse => {
            let mut n_guarded = 0usize;
            let values = raw
                .iter()
                .map(|v| match v {
                    Some(x) if *x != 0.0 => Some(1.0 / x),
                    Some(_) => {
                        n_guarded += 1;
                        None
                    }
                    None => None,
                })
                .collect();
            (
                values,
                json!({ "transform": "inverse", "n_guarded_zero": n_guarded }),
            )
        }
        Transform::BoxCox => {
            let positive: Vec<f64> = vals.iter().copied().filter(|x| *x > 0.0).collect();
            let n_guarded = vals.len() - positive.len();
            let lambda = match fixed_lambda {
                Some(l) => l,
                None if positive.len() >= 2 => fit_boxcox_lambda(&positive),
                None => 1.0,
            };
            let values = raw
                .iter()
                .map(|v| match v {
                    Some(x) if *x > 0.0 => Some(boxcox_transform(*x, lambda)),
                    _ => None,
                })
                .collect();
            (
                values,
                json!({ "transform": "boxcox", "lambda": lambda, "lambda_fitted": fixed_lambda.is_none(), "n_guarded_nonpositive": n_guarded }),
            )
        }
        Transform::Bin | Transform::OneHot => unreachable!("handled separately"),
    };

    let field_values = values
        .into_iter()
        .map(|v| v.map(FieldValue::Float).unwrap_or(FieldValue::Null))
        .collect();
    (
        OutputColumn {
            name,
            field_type: FieldType::Float,
            values: field_values,
        },
        stats,
    )
}

// ── Binning ──────────────────────────────────────────────────────────────────

fn build_bin_column(
    raw: &[Option<f64>],
    field_name: &str,
    suffix: &str,
    k: u64,
    method: BinMethod,
) -> (OutputColumn, Value) {
    let vals: Vec<f64> = raw.iter().filter_map(|v| *v).collect();
    let name = format!("{field_name}{suffix}");
    let k = k.max(2) as usize;

    let breaks: Vec<f64> = if vals.is_empty() {
        Vec::new()
    } else {
        match method {
            BinMethod::EqualInterval => {
                let (lo, hi) = min_max(&vals);
                let step = (hi - lo) / k as f64;
                if step > 0.0 {
                    (1..k).map(|i| lo + step * i as f64).collect()
                } else {
                    Vec::new()
                }
            }
            BinMethod::Quantile => {
                let mut sorted = vals.clone();
                sorted.sort_by(|a, b| a.total_cmp(b));
                (1..k)
                    .map(|i| percentile(&sorted, 100.0 * i as f64 / k as f64))
                    .collect()
            }
            BinMethod::StdDev => {
                // Simplified standard-deviation classification: `k` equal-width
                // bins of one sample standard deviation each, centered on the
                // mean (breaks at mean + sd*(i - k/2)). This is not ArcGIS's
                // exact standard-deviation-interval algorithm (which snaps to
                // a fixed interval size like 1 or 1/2 sd and can vary the bin
                // count), but it is deterministic and keeps the mean centered.
                let (mean, sd) = mean_sd(&vals);
                if sd > 0.0 {
                    let half = k as f64 / 2.0;
                    (1..k).map(|i| mean + sd * (i as f64 - half)).collect()
                } else {
                    Vec::new()
                }
            }
        }
    };

    let values: Vec<FieldValue> = raw
        .iter()
        .map(|v| match v {
            Some(x) => FieldValue::Integer(assign_class(*x, &breaks)),
            None => FieldValue::Null,
        })
        .collect();

    let stats = json!({
        "transform": "bin",
        "bin_method": method.as_str(),
        "bins": k,
        "breaks": breaks,
    });

    (
        OutputColumn {
            name,
            field_type: FieldType::Integer,
            values,
        },
        stats,
    )
}

fn assign_class(x: f64, breaks: &[f64]) -> i64 {
    let mut c = 1i64;
    for &b in breaks {
        if x > b {
            c += 1;
        } else {
            break;
        }
    }
    c
}

// ── One-hot encoding ─────────────────────────────────────────────────────────

fn category_key(v: &FieldValue) -> Option<String> {
    match v {
        FieldValue::Null => None,
        FieldValue::Text(s) | FieldValue::Date(s) | FieldValue::DateTime(s) => Some(s.clone()),
        FieldValue::Integer(i) => Some(i.to_string()),
        FieldValue::Float(f) => Some(format!("{f}")),
        FieldValue::Boolean(b) => Some(b.to_string()),
        FieldValue::Blob(_) => None,
    }
}

fn sanitize_category(s: &str) -> String {
    let mut out: String = s
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() {
                c.to_ascii_lowercase()
            } else {
                '_'
            }
        })
        .collect();
    out.truncate(32);
    if out.is_empty() {
        out.push_str("na");
    }
    out
}

fn build_onehot_columns(
    layer: &Layer,
    idx: usize,
    field_name: &str,
    suffix: &str,
    ctx: &ToolContext,
) -> (Vec<OutputColumn>, Value) {
    let n = layer.len();
    let mut counts: HashMap<String, usize> = HashMap::new();
    let mut raw_keys: Vec<Option<String>> = Vec::with_capacity(n);
    for feat in layer.features.iter() {
        let key = feat.attributes.get(idx).and_then(category_key);
        if let Some(k) = &key {
            *counts.entry(k.clone()).or_insert(0) += 1;
        }
        raw_keys.push(key);
    }

    let mut cats: Vec<(String, usize)> = counts.into_iter().collect();
    // Descending frequency, alphabetical tie-break, for a deterministic order.
    cats.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    let capped = cats.len() > MAX_ONEHOT_CATEGORIES;
    let kept: Vec<String> = cats
        .iter()
        .take(MAX_ONEHOT_CATEGORIES)
        .map(|(k, _)| k.clone())
        .collect();
    if capped {
        ctx.progress.info(&format!(
            "transform_fields: field '{field_name}' has {} distinct categories; capped to {MAX_ONEHOT_CATEGORIES}, remainder grouped into '{field_name}{suffix}other'",
            cats.len()
        ));
    }
    let kept_set: HashSet<&str> = kept.iter().map(|s| s.as_str()).collect();

    let mut cols = Vec::with_capacity(kept.len() + usize::from(capped));
    for cat in &kept {
        let col_name = format!("{field_name}{suffix}{}", sanitize_category(cat));
        let values: Vec<FieldValue> = raw_keys
            .iter()
            .map(|k| match k {
                Some(kk) if kk == cat => FieldValue::Integer(1),
                Some(_) => FieldValue::Integer(0),
                None => FieldValue::Null,
            })
            .collect();
        cols.push(OutputColumn {
            name: col_name,
            field_type: FieldType::Integer,
            values,
        });
    }
    if capped {
        let col_name = format!("{field_name}{suffix}other");
        let values: Vec<FieldValue> = raw_keys
            .iter()
            .map(|k| match k {
                Some(kk) if !kept_set.contains(kk.as_str()) => FieldValue::Integer(1),
                Some(_) => FieldValue::Integer(0),
                None => FieldValue::Null,
            })
            .collect();
        cols.push(OutputColumn {
            name: col_name,
            field_type: FieldType::Integer,
            values,
        });
    }

    let stats = json!({
        "transform": "onehot",
        "categories": kept,
        "n_categories_total": cats.len(),
        "capped": capped,
    });
    (cols, stats)
}

// ── Box-Cox fitting (golden-section search on the profile log-likelihood) ───

fn boxcox_transform(x: f64, lambda: f64) -> f64 {
    if lambda.abs() < 1e-8 {
        x.ln()
    } else {
        (x.powf(lambda) - 1.0) / lambda
    }
}

fn boxcox_loglik(vals: &[f64], lambda: f64) -> f64 {
    let n = vals.len() as f64;
    if n <= 0.0 {
        return f64::NEG_INFINITY;
    }
    let transformed: Vec<f64> = vals.iter().map(|&x| boxcox_transform(x, lambda)).collect();
    let mean = transformed.iter().sum::<f64>() / n;
    let var = transformed.iter().map(|t| (t - mean).powi(2)).sum::<f64>() / n;
    if var <= 0.0 {
        return f64::NEG_INFINITY;
    }
    let sum_log_x: f64 = vals.iter().map(|x| x.ln()).sum();
    -0.5 * n * var.ln() + (lambda - 1.0) * sum_log_x
}

fn fit_boxcox_lambda(positive_vals: &[f64]) -> f64 {
    let f = |lambda: f64| boxcox_loglik(positive_vals, lambda);
    golden_section_max(-5.0, 5.0, &f, 60)
}

/// Golden-section search for the argmax of a unimodal `f` over `[a, b]`.
fn golden_section_max<F: Fn(f64) -> f64>(mut a: f64, mut b: f64, f: &F, iters: usize) -> f64 {
    let gr = (5f64.sqrt() - 1.0) / 2.0; // ~0.618
    let mut c = b - gr * (b - a);
    let mut d = a + gr * (b - a);
    let mut fc = f(c);
    let mut fd = f(d);
    for _ in 0..iters {
        if fc > fd {
            b = d;
            d = c;
            fd = fc;
            c = b - gr * (b - a);
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + gr * (b - a);
            fd = f(d);
        }
    }
    (a + b) / 2.0
}

// ── Basic statistics helpers ─────────────────────────────────────────────────

fn mean_sd(vals: &[f64]) -> (f64, f64) {
    let n = vals.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    let mean = vals.iter().sum::<f64>() / n as f64;
    let var = if n > 1 {
        vals.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0)
    } else {
        0.0
    };
    (mean, var.sqrt())
}

fn min_max(vals: &[f64]) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &v in vals {
        if v.is_finite() {
            lo = lo.min(v);
            hi = hi.max(v);
        }
    }
    if !lo.is_finite() {
        (0.0, 0.0)
    } else {
        (lo, hi)
    }
}

/// Linear-interpolation percentile (numpy `'linear'` convention) over an
/// already-sorted slice.
fn percentile(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let n = sorted.len();
    if n == 1 {
        return sorted[0];
    }
    let rank = p / 100.0 * (n as f64 - 1.0);
    let lo = rank.floor() as usize;
    let hi = rank.ceil() as usize;
    if lo == hi {
        sorted[lo]
    } else {
        let frac = rank - lo as f64;
        sorted[lo] + (sorted[hi] - sorted[lo]) * frac
    }
}

// ── Parameters ───────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Transform {
    ZScore,
    MinMax,
    Robust,
    Log,
    Log1p,
    Sqrt,
    BoxCox,
    Inverse,
    Bin,
    OneHot,
}

impl Transform {
    fn as_str(self) -> &'static str {
        match self {
            Transform::ZScore => "zscore",
            Transform::MinMax => "minmax",
            Transform::Robust => "robust",
            Transform::Log => "log",
            Transform::Log1p => "log1p",
            Transform::Sqrt => "sqrt",
            Transform::BoxCox => "boxcox",
            Transform::Inverse => "inverse",
            Transform::Bin => "bin",
            Transform::OneHot => "onehot",
        }
    }

    fn default_suffix(self) -> &'static str {
        match self {
            Transform::ZScore => "_z",
            Transform::MinMax => "_minmax",
            Transform::Robust => "_robust",
            Transform::Log => "_log",
            Transform::Log1p => "_log1p",
            Transform::Sqrt => "_sqrt",
            Transform::BoxCox => "_boxcox",
            Transform::Inverse => "_inv",
            Transform::Bin => "_bin",
            Transform::OneHot => "_",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum BinMethod {
    EqualInterval,
    Quantile,
    StdDev,
}

impl BinMethod {
    fn as_str(self) -> &'static str {
        match self {
            BinMethod::EqualInterval => "equal_interval",
            BinMethod::Quantile => "quantile",
            BinMethod::StdDev => "std_dev",
        }
    }
}

struct Params {
    fields: Vec<String>,
    transform: Transform,
    bins: u64,
    bin_method: BinMethod,
    boxcox_lambda: Option<f64>,
    suffix: Option<String>,
    drop_input: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let fields_str = require_str(args, "fields")?;
    let fields: Vec<String> = fields_str
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_string)
        .collect();
    if fields.is_empty() {
        return Err(ToolError::Validation(
            "'fields' must name at least one field".to_string(),
        ));
    }

    let transform_str = require_str(args, "transform")?;
    let transform = match transform_str {
        "zscore" => Transform::ZScore,
        "minmax" => Transform::MinMax,
        "robust" => Transform::Robust,
        "log" => Transform::Log,
        "log1p" => Transform::Log1p,
        "sqrt" => Transform::Sqrt,
        "boxcox" => Transform::BoxCox,
        "inverse" => Transform::Inverse,
        "bin" => Transform::Bin,
        "onehot" => Transform::OneHot,
        o => {
            return Err(ToolError::Validation(format!(
                "'transform' must be one of zscore/minmax/robust/log/log1p/sqrt/boxcox/inverse/bin/onehot, got '{o}'"
            )))
        }
    };

    let bins = parse_optional_u64(args, "bins")?.unwrap_or(5);
    if transform == Transform::Bin && bins < 2 {
        return Err(ToolError::Validation(
            "'bins' must be at least 2".to_string(),
        ));
    }

    let bin_method = match args
        .get("bin_method")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some("equal_interval") => BinMethod::EqualInterval,
        Some("quantile") => BinMethod::Quantile,
        Some("std_dev") => BinMethod::StdDev,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'bin_method' must be equal_interval/quantile/std_dev, got '{o}'"
            )))
        }
    };

    let boxcox_lambda = parse_optional_f64(args, "boxcox_lambda")?;
    let suffix = parse_optional_str(args, "suffix")?.map(str::to_string);
    let drop_input = parse_optional_bool(args, "drop_input")?.unwrap_or(false);

    Ok(Params {
        fields,
        transform,
        bins,
        bin_method,
        boxcox_lambda,
        suffix,
        drop_input,
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

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Geometry, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn numeric_layer(vals: &[f64]) -> String {
        let mut l = Layer::new("t")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (i, v) in vals.iter().enumerate() {
            l.add_feature(Some(Geometry::point(i as f64, 0.0)), &[("v", (*v).into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn categorical_layer(cats: &[&str]) -> String {
        let mut l = Layer::new("t")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("cat", FieldType::Text));
        for (i, c) in cats.iter().enumerate() {
            l.add_feature(
                Some(Geometry::point(i as f64, 0.0)),
                &[("cat", (*c).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> Layer {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = TransformFieldsTool.run(&args, &ctx()).unwrap();
        load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    fn col_f64(l: &Layer, name: &str) -> Vec<Option<f64>> {
        let i = l.schema.field_index(name).unwrap();
        l.features
            .iter()
            .map(|f| f.attributes[i].as_f64())
            .collect()
    }

    /// zscore: recomputing mean/sd on the output field independently should
    /// give mean ~0, sd ~1.
    #[test]
    fn zscore_output_has_mean_zero_sd_one() {
        let vals = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let l = run(json!({
            "input": numeric_layer(&vals),
            "fields": "v",
            "transform": "zscore",
        }));
        let z: Vec<f64> = col_f64(&l, "v_z").into_iter().map(|v| v.unwrap()).collect();
        let (mean, sd) = mean_sd(&z);
        assert!(mean.abs() < 1e-9, "mean should be ~0, got {mean}");
        assert!((sd - 1.0).abs() < 1e-9, "sd should be ~1, got {sd}");
    }

    /// minmax output stays within [0, 1] and hits both endpoints.
    #[test]
    fn minmax_output_in_unit_range() {
        let vals = [3.0, 1.0, 8.0, 5.0, 2.0];
        let l = run(json!({
            "input": numeric_layer(&vals),
            "fields": "v",
            "transform": "minmax",
        }));
        let m: Vec<f64> = col_f64(&l, "v_minmax")
            .into_iter()
            .map(|v| v.unwrap())
            .collect();
        for x in &m {
            assert!((0.0..=1.0).contains(x), "value {x} out of [0,1]");
        }
        assert!((m.iter().cloned().fold(f64::INFINITY, f64::min) - 0.0).abs() < 1e-9);
        assert!((m.iter().cloned().fold(f64::NEG_INFINITY, f64::max) - 1.0).abs() < 1e-9);
    }

    /// quantile binning should produce classes of roughly equal size.
    #[test]
    fn quantile_bin_produces_near_equal_counts() {
        let vals: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let l = run(json!({
            "input": numeric_layer(&vals),
            "fields": "v",
            "transform": "bin",
            "bins": 5,
            "bin_method": "quantile",
        }));
        let idx = l.schema.field_index("v_bin").unwrap();
        let mut counts = [0i64; 6]; // classes are 1-based
        for f in &l.features {
            let c = f.attributes[idx].as_i64().unwrap();
            counts[c as usize] += 1;
        }
        for (c, &count) in counts.iter().enumerate().skip(1) {
            assert!(
                (count - 20).abs() <= 1,
                "class {c} has {count} features, expected ~20"
            );
        }
    }

    /// onehot columns for a fully-populated categorical field sum to 1 per row.
    #[test]
    fn onehot_columns_sum_to_one_per_row() {
        let cats = ["a", "b", "a", "c", "b", "a"];
        let l = run(json!({
            "input": categorical_layer(&cats),
            "fields": "cat",
            "transform": "onehot",
        }));
        let onehot_indices: Vec<usize> = l
            .schema
            .fields()
            .iter()
            .enumerate()
            .filter(|(_, f)| f.name.starts_with("cat_"))
            .map(|(i, _)| i)
            .collect();
        assert_eq!(onehot_indices.len(), 3, "expects one column per category");
        for f in &l.features {
            let sum: i64 = onehot_indices
                .iter()
                .map(|&i| f.attributes[i].as_i64().unwrap())
                .sum();
            assert_eq!(sum, 1, "one-hot row should sum to 1");
        }
    }

    /// log of non-positive values is guarded to null, not a panic/NaN.
    #[test]
    fn log_of_nonpositive_is_null() {
        let vals = [1.0, -2.0, 0.0, 10.0];
        let l = run(json!({
            "input": numeric_layer(&vals),
            "fields": "v",
            "transform": "log",
        }));
        let log_vals = col_f64(&l, "v_log");
        assert!(log_vals[0].is_some());
        assert!(log_vals[1].is_none(), "log(-2) should be null");
        assert!(log_vals[2].is_none(), "log(0) should be null");
        assert!(log_vals[3].is_some());
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            TransformFieldsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err()); // no fields/transform
        assert!(bad(json!({ "input": "a.geojson", "fields": "v" })).is_err()); // no transform
        assert!(
            bad(json!({ "input": "a.geojson", "fields": "v", "transform": "cube_root" })).is_err()
        );
        assert!(bad(json!({
            "input": "a.geojson", "fields": "v", "transform": "bin", "bin_method": "kmeans"
        }))
        .is_err());
        assert!(bad(json!({
            "input": "a.geojson", "fields": "v", "transform": "boxcox", "boxcox_lambda": "abc"
        }))
        .is_err());
        assert!(bad(json!({ "input": "a.geojson", "fields": "v", "transform": "zscore" })).is_ok());
    }

    #[test]
    fn rejects_unknown_field() {
        let vals = [1.0, 2.0, 3.0];
        let args: ToolArgs = serde_json::from_value(json!({
            "input": numeric_layer(&vals),
            "fields": "does_not_exist",
            "transform": "zscore",
        }))
        .unwrap();
        assert!(TransformFieldsTool.run(&args, &ctx()).is_err());
    }
}
