//! GeoLibre tool: combine several numeric attributes into a single composite
//! index.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Calculate Composite Index* (Spatial
//! Statistics) — the standard methodology for vulnerability / deprivation / SDG
//! indices. Nothing in the bundled whitebox suite builds an attribute index
//! (`weighted_overlay` / `weighted_sum` are raster-only); this is the vector-side
//! sibling of `fuzzy_overlay` and feeds straight into the hot-spot and
//! regionalization tools already in the suite.
//!
//! Each input variable is first scaled (min-max, z-score, or percentile), with
//! an optional per-variable `reverse` so that "higher = worse" variables count
//! the right way. The scaled variables are combined (weighted mean, sum, or
//! geometric mean) and the result is optionally rescaled (min-max, 0–100, or
//! z-score). Outputs: the composite `index`, its `index_rank` and
//! `index_pctl`, and one `<field>_scaled` column per input variable.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CalculateCompositeIndexTool;

impl Tool for CalculateCompositeIndexTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "calculate_composite_index",
            display_name: "Calculate Composite Index",
            summary: "Combine several numeric attributes into a composite index (like ArcGIS Calculate Composite Index): per-variable scaling (min-max/z-score/percentile) with optional reversal, weighted combination (mean/sum/geometric mean), and output rescaling — the vector-side sibling of fuzzy_overlay.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (or table) with numeric variable fields.",
                    required: true,
                },
                ToolParamSpec {
                    name: "fields",
                    description: "Comma-separated variable fields; append ':reverse' to a field whose high values should lower the index (e.g. 'income:reverse,unemployment').",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector layer with the index fields added. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "scaling",
                    description: "Per-variable scaling: 'minmax' (default), 'zscore', 'percentile', or 'none'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "weights",
                    description: "Comma-separated weights matching 'fields' (default: equal weights).",
                    required: false,
                },
                ToolParamSpec {
                    name: "combine",
                    description: "How to combine scaled variables: 'mean' (default), 'sum', or 'geometric_mean'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_range",
                    description: "Rescale the combined index: 'minmax' (0–1, default), 'zero_to_100', 'zscore', or 'none'.",
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
        let n = layer.features.len();

        // Resolve variable field indices.
        let mut var_idx = Vec::with_capacity(prm.vars.len());
        for v in &prm.vars {
            let idx = layer
                .schema
                .field_index(&v.name)
                .ok_or_else(|| ToolError::Validation(format!("field '{}' not found", v.name)))?;
            var_idx.push(idx);
        }
        let weights = match &prm.weights {
            Some(w) if w.len() != prm.vars.len() => {
                return Err(ToolError::Validation(format!(
                    "{} weights given for {} fields",
                    w.len(),
                    prm.vars.len()
                )))
            }
            Some(w) => w.clone(),
            None => vec![1.0; prm.vars.len()],
        };

        // Read raw values per variable.
        let mut raw: Vec<Vec<Option<f64>>> = vec![vec![None; n]; prm.vars.len()];
        for (fi, feat) in layer.features.iter().enumerate() {
            for (vi, &idx) in var_idx.iter().enumerate() {
                raw[vi][fi] = feat.attributes.get(idx).and_then(|v| v.as_f64());
            }
        }

        ctx.progress.info(&format!(
            "composite index over {} variable(s)",
            prm.vars.len()
        ));

        // Scale each variable (missing values -> 0.5 for bounded scalings, 0 for z-score).
        let mut scaled: Vec<Vec<f64>> = Vec::with_capacity(prm.vars.len());
        for (vi, col) in raw.iter().enumerate() {
            let mut s = scale_column(col, prm.scaling);
            if prm.vars[vi].reverse {
                reverse_scaled(&mut s, prm.scaling);
            }
            scaled.push(s);
        }

        // Combine per feature.
        let wsum: f64 = weights.iter().sum();
        let mut combined = vec![0.0f64; n];
        for fi in 0..n {
            combined[fi] = match prm.combine {
                Combine::Mean => {
                    if wsum > 0.0 {
                        (0..prm.vars.len())
                            .map(|vi| weights[vi] * scaled[vi][fi])
                            .sum::<f64>()
                            / wsum
                    } else {
                        0.0
                    }
                }
                Combine::Sum => (0..prm.vars.len())
                    .map(|vi| weights[vi] * scaled[vi][fi])
                    .sum::<f64>(),
                Combine::GeometricMean => {
                    // Weighted geometric mean; clamp to >0 to keep it defined.
                    if wsum > 0.0 {
                        let logsum: f64 = (0..prm.vars.len())
                            .map(|vi| weights[vi] * scaled[vi][fi].max(1e-9).ln())
                            .sum();
                        (logsum / wsum).exp()
                    } else {
                        0.0
                    }
                }
            };
        }

        // Rescale the combined index.
        let index = rescale_output(&combined, prm.output_range);

        // Ranks (1 = highest index) and percentiles.
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| index[b].total_cmp(&index[a]));
        let mut rank = vec![0i64; n];
        let mut pctl = vec![0.0f64; n];
        for (r, &fi) in order.iter().enumerate() {
            rank[fi] = (r + 1) as i64;
            pctl[fi] = if n > 1 {
                100.0 * (1.0 - r as f64 / (n as f64 - 1.0))
            } else {
                100.0
            };
        }

        // Append fields.
        layer.add_field(FieldDef::new("index", FieldType::Float));
        layer.add_field(FieldDef::new("index_rank", FieldType::Integer));
        layer.add_field(FieldDef::new("index_pctl", FieldType::Float));
        for v in &prm.vars {
            layer.add_field(FieldDef::new(
                format!("{}_scaled", v.name),
                FieldType::Float,
            ));
        }
        for fi in 0..n {
            let feat = &mut layer.features[fi];
            feat.attributes.push(FieldValue::Float(index[fi]));
            feat.attributes.push(FieldValue::Integer(rank[fi]));
            feat.attributes.push(FieldValue::Float(pctl[fi]));
            for s in &scaled {
                feat.attributes.push(FieldValue::Float(s[fi]));
            }
        }

        let (imin, imax) = min_max(&index);
        let out_path = write_or_store_layer(layer, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("index_min".to_string(), json!(imin));
        outputs.insert("index_max".to_string(), json!(imax));
        Ok(ToolRunResult { outputs })
    }
}

// ── Scaling ──────────────────────────────────────────────────────────────────

fn scale_column(col: &[Option<f64>], scaling: Scaling) -> Vec<f64> {
    let vals: Vec<f64> = col.iter().filter_map(|v| *v).collect();
    match scaling {
        Scaling::None => col.iter().map(|v| v.unwrap_or(0.0)).collect(),
        Scaling::MinMax => {
            let (lo, hi) = min_max(&vals);
            let range = hi - lo;
            col.iter()
                .map(|v| match v {
                    Some(x) if range > 0.0 => (x - lo) / range,
                    Some(_) => 0.5,
                    None => 0.5,
                })
                .collect()
        }
        Scaling::ZScore => {
            let mean = if vals.is_empty() {
                0.0
            } else {
                vals.iter().sum::<f64>() / vals.len() as f64
            };
            let var = if vals.len() > 1 {
                vals.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (vals.len() as f64 - 1.0)
            } else {
                0.0
            };
            let sd = var.sqrt();
            col.iter()
                .map(|v| match v {
                    Some(x) if sd > 0.0 => (x - mean) / sd,
                    _ => 0.0,
                })
                .collect()
        }
        Scaling::Percentile => {
            // Rank each value to [0,1] by its position in the sorted set.
            let mut sorted = vals.clone();
            sorted.sort_by(|a, b| a.total_cmp(b));
            let m = sorted.len();
            col.iter()
                .map(|v| match v {
                    Some(x) if m > 1 => {
                        // fraction of values <= x
                        let cnt = sorted.iter().filter(|&&s| s <= *x).count();
                        (cnt as f64 - 1.0) / (m as f64 - 1.0)
                    }
                    Some(_) => 0.5,
                    None => 0.5,
                })
                .collect()
        }
    }
}

/// Reverses a scaled column so high inputs become low contributions.
fn reverse_scaled(s: &mut [f64], scaling: Scaling) {
    match scaling {
        Scaling::MinMax | Scaling::Percentile => {
            for v in s.iter_mut() {
                *v = 1.0 - *v;
            }
        }
        Scaling::ZScore | Scaling::None => {
            for v in s.iter_mut() {
                *v = -*v;
            }
        }
    }
}

fn rescale_output(combined: &[f64], range: OutRange) -> Vec<f64> {
    match range {
        OutRange::None => combined.to_vec(),
        OutRange::MinMax | OutRange::ZeroTo100 => {
            let (lo, hi) = min_max(combined);
            let span = hi - lo;
            let scale = if matches!(range, OutRange::ZeroTo100) {
                100.0
            } else {
                1.0
            };
            combined
                .iter()
                .map(|x| {
                    if span > 0.0 {
                        scale * (x - lo) / span
                    } else {
                        0.0
                    }
                })
                .collect()
        }
        OutRange::ZScore => {
            let n = combined.len();
            let mean = if n == 0 {
                0.0
            } else {
                combined.iter().sum::<f64>() / n as f64
            };
            let var = if n > 1 {
                combined.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0)
            } else {
                0.0
            };
            let sd = var.sqrt();
            combined
                .iter()
                .map(|x| if sd > 0.0 { (x - mean) / sd } else { 0.0 })
                .collect()
        }
    }
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

// ── Parameters ───────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Scaling {
    MinMax,
    ZScore,
    Percentile,
    None,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Combine {
    Mean,
    Sum,
    GeometricMean,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum OutRange {
    MinMax,
    ZeroTo100,
    ZScore,
    None,
}

struct VarSpec {
    name: String,
    reverse: bool,
}

struct Params {
    vars: Vec<VarSpec>,
    scaling: Scaling,
    weights: Option<Vec<f64>>,
    combine: Combine,
    output_range: OutRange,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let fields = require_str(args, "fields")?;
    let vars: Vec<VarSpec> = fields
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|spec| {
            let mut it = spec.splitn(2, ':');
            let name = it.next().unwrap().trim().to_string();
            let reverse = matches!(it.next().map(|m| m.trim().to_ascii_lowercase()), Some(m) if m == "reverse" || m == "rev");
            VarSpec { name, reverse }
        })
        .collect();
    if vars.is_empty() {
        return Err(ToolError::Validation(
            "'fields' must name at least one variable".to_string(),
        ));
    }
    let scaling = match args.get("scaling").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("minmax") => Scaling::MinMax,
        Some("zscore") => Scaling::ZScore,
        Some("percentile") => Scaling::Percentile,
        Some("none") => Scaling::None,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'scaling' must be minmax/zscore/percentile/none, got '{o}'"
            )))
        }
    };
    let weights = match args.get("weights").and_then(Value::as_str) {
        None => None,
        Some(s) if s.trim().is_empty() => None,
        Some(s) => Some(
            s.split(',')
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(|w| {
                    w.parse::<f64>()
                        .map_err(|_| ToolError::Validation(format!("weight '{w}' is not a number")))
                })
                .collect::<Result<Vec<_>, _>>()?,
        ),
    };
    let combine = match args.get("combine").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("mean") => Combine::Mean,
        Some("sum") => Combine::Sum,
        Some("geometric_mean") => Combine::GeometricMean,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'combine' must be mean/sum/geometric_mean, got '{o}'"
            )))
        }
    };
    let output_range = match args
        .get("output_range")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some("minmax") => OutRange::MinMax,
        Some("zero_to_100") | Some("0to100") => OutRange::ZeroTo100,
        Some("zscore") => OutRange::ZScore,
        Some("none") => OutRange::None,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'output_range' must be minmax/zero_to_100/zscore/none, got '{o}'"
            )))
        }
    };
    Ok(Params {
        vars,
        scaling,
        weights,
        combine,
        output_range,
    })
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
    use wbvector::{memory_store, Geometry, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn layer(rows: &[(f64, f64)]) -> String {
        let mut l = Layer::new("z")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("a", FieldType::Float));
        l.add_field(FieldDef::new("b", FieldType::Float));
        for (i, (a, b)) in rows.iter().enumerate() {
            l.add_feature(
                Some(Geometry::point(i as f64, 0.0)),
                &[("a", (*a).into()), ("b", (*b).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> Layer {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CalculateCompositeIndexTool.run(&args, &ctx()).unwrap();
        load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    fn col(l: &Layer, name: &str) -> Vec<f64> {
        let i = l.schema.field_index(name).unwrap();
        l.features
            .iter()
            .map(|f| f.attributes[i].as_f64().unwrap())
            .collect()
    }

    /// Min-max scaling + equal-weight mean produces a monotone index; highest
    /// combined inputs get rank 1.
    #[test]
    fn builds_monotone_index() {
        // a and b both increase together, so the index should increase too.
        let rows = [(0.0, 0.0), (5.0, 5.0), (10.0, 10.0)];
        let l = run(json!({ "input": layer(&rows), "fields": "a,b" }));
        let idx = col(&l, "index");
        assert!(idx[0] < idx[1] && idx[1] < idx[2], "index should increase");
        assert!((idx[0] - 0.0).abs() < 1e-9 && (idx[2] - 1.0).abs() < 1e-9);
        let rank = l.schema.field_index("index_rank").unwrap();
        assert_eq!(
            l.features[2].attributes[rank].as_i64(),
            Some(1),
            "highest = rank 1"
        );
    }

    /// A reversed variable inverts its contribution: high 'a' should lower the
    /// index.
    #[test]
    fn reverse_inverts_contribution() {
        let rows = [(0.0, 0.0), (10.0, 0.0)]; // only 'a' varies
        let l = run(json!({ "input": layer(&rows), "fields": "a:reverse,b" }));
        let idx = col(&l, "index");
        assert!(idx[0] > idx[1], "reversed high 'a' should reduce the index");
    }

    /// Weights bias the combination toward the weighted variable.
    #[test]
    fn weights_bias_combination() {
        let rows = [(10.0, 0.0), (0.0, 10.0)];
        // Weight 'a' heavily: row 0 (high a) should outrank row 1.
        let l = run(json!({ "input": layer(&rows), "fields": "a,b", "weights": "9,1" }));
        let idx = col(&l, "index");
        assert!(idx[0] > idx[1]);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CalculateCompositeIndexTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err()); // no fields
        assert!(bad(json!({ "input": "a.geojson", "fields": "a", "scaling": "log" })).is_err());
        assert!(
            bad(json!({ "input": "a.geojson", "fields": "a,b", "combine": "median" })).is_err()
        );
        assert!(bad(json!({ "input": "a.geojson", "fields": "a,b" })).is_ok());
    }
}
