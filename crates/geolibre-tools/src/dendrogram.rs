//! GeoLibre tool: hierarchical class-separability tree from cluster signatures.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Dendrogram* (Spatial Analyst). The
//! shipped `multivariate_clustering` and the bundled `k_means_clustering` /
//! `modified_k_means_clustering` all emit class assignments with no way to
//! judge whether the resulting classes are actually distinct. A user who asks
//! for 8 clusters gets 8 clusters, with nothing to reveal that three of them
//! are near-duplicates and the real structure is 6.
//!
//! This is that diagnostic. Classes are characterized by their mean vector in
//! the chosen feature space (and, under `distance = variance`, their pooled
//! within-class spread), then agglomeratively merged by average linkage. The
//! merge *distance* is what carries the information: classes that join very
//! early are the ones a clustering run failed to separate.
//!
//! Class counts here are small (tens, not thousands), so the naive O(n^3)
//! agglomeration is entirely adequate and avoids pulling in a clustering crate.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Distance between two class signatures.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Distance {
    /// Euclidean distance between class means only.
    MeanOnly,
    /// Mean distance normalized by the pooled within-class spread, so a pair of
    /// tight classes reads as better separated than a pair of diffuse ones at
    /// the same mean separation.
    Variance,
}

pub struct DendrogramTool;

impl Tool for DendrogramTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "dendrogram",
            display_name: "Dendrogram",
            summary: "Build a hierarchical merge tree over class signatures to show how separable the classes are, like ArcGIS Dendrogram.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Clustered feature layer carrying a class field and the numeric fields defining the feature space.",
                    required: true,
                },
                ToolParamSpec {
                    name: "class_field",
                    description: "Field holding the class / cluster identifier.",
                    required: true,
                },
                ToolParamSpec {
                    name: "fields",
                    description: "Comma-separated numeric fields defining the feature space.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output merge table path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distance",
                    description: "Distance between signatures: 'variance' (default; normalized by within-class spread) or 'mean_only'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "standardize",
                    description: "Z-standardize each field before computing distances so large-range fields do not dominate (default true).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_text",
                    description: "Optional path for a plain-text indented rendering of the tree.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "class_field")?;
        let fields = parse_fields(args)?;
        if fields.is_empty() {
            return Err(ToolError::Validation(
                "'fields' must name at least one numeric field".to_string(),
            ));
        }
        parse_distance(args)?;
        parse_optional_bool(args, "standardize")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let class_field = require_str(args, "class_field")?;
        let fields = parse_fields(args)?;
        let output = parse_optional_str(args, "output")?;
        let distance = parse_distance(args)?;
        let standardize = parse_optional_bool(args, "standardize")?.unwrap_or(true);
        let text_out = parse_optional_str(args, "output_text")?.map(String::from);

        let layer = load_input_layer(input)?;
        if layer.schema.field_index(class_field).is_none() {
            return Err(ToolError::Validation(format!(
                "class_field '{class_field}' not found on the input layer"
            )));
        }
        for f in &fields {
            if layer.schema.field_index(f).is_none() {
                return Err(ToolError::Validation(format!(
                    "field '{f}' not found on the input layer"
                )));
            }
        }

        // Collect per-class value vectors.
        let mut by_class: BTreeMap<String, Vec<Vec<f64>>> = BTreeMap::new();
        for feat in layer.features.iter() {
            let class = match feat.get(&layer.schema, class_field) {
                Ok(v) => field_string(v),
                Err(_) => continue,
            };
            if class.is_empty() {
                continue;
            }
            let mut row = Vec::with_capacity(fields.len());
            let mut ok = true;
            for f in &fields {
                match feat.get(&layer.schema, f).ok().and_then(|v| v.as_f64()) {
                    Some(x) if x.is_finite() => row.push(x),
                    // A row with a missing predictor cannot be placed in the
                    // feature space; drop it rather than imputing silently.
                    _ => {
                        ok = false;
                        break;
                    }
                }
            }
            if ok {
                by_class.entry(class).or_default().push(row);
            }
        }

        if by_class.len() < 2 {
            return Err(ToolError::Validation(format!(
                "need at least 2 classes to build a dendrogram, found {}",
                by_class.len()
            )));
        }

        // Optional z-standardization across the whole dataset, so a field
        // measured in millions does not swamp one measured in fractions.
        let dim = fields.len();
        let (centers, scales) = if standardize {
            let all: Vec<&Vec<f64>> = by_class.values().flatten().collect();
            let n = all.len() as f64;
            let mut mean = vec![0.0; dim];
            for r in &all {
                for d in 0..dim {
                    mean[d] += r[d];
                }
            }
            for m in mean.iter_mut() {
                *m /= n;
            }
            let mut sd = vec![0.0; dim];
            for r in &all {
                for d in 0..dim {
                    sd[d] += (r[d] - mean[d]).powi(2);
                }
            }
            for s in sd.iter_mut() {
                *s = (*s / n).sqrt();
                if *s <= f64::EPSILON {
                    *s = 1.0; // constant field: leave it unscaled rather than dividing by 0
                }
            }
            (mean, sd)
        } else {
            (vec![0.0; dim], vec![1.0; dim])
        };

        // Build one signature per class.
        let mut sigs: Vec<Sig> = Vec::with_capacity(by_class.len());
        for (name, rows) in &by_class {
            let n = rows.len() as f64;
            let mut mean = vec![0.0; dim];
            for r in rows {
                for d in 0..dim {
                    mean[d] += (r[d] - centers[d]) / scales[d];
                }
            }
            for m in mean.iter_mut() {
                *m /= n;
            }
            let mut var = 0.0;
            for r in rows {
                for d in 0..dim {
                    let z = (r[d] - centers[d]) / scales[d];
                    var += (z - mean[d]).powi(2);
                }
            }
            var /= (rows.len() * dim).max(1) as f64;
            sigs.push(Sig {
                labels: vec![name.clone()],
                mean,
                spread: var.sqrt(),
                count: rows.len(),
            });
        }

        ctx.progress.info(&format!(
            "agglomerating {} class signature(s) over {dim} field(s)",
            sigs.len()
        ));

        // Agglomerative merge, average linkage.
        let mut out = Layer::new("dendrogram");
        out.add_field(FieldDef::new("step", FieldType::Integer));
        out.add_field(FieldDef::new("merged_a", FieldType::Text));
        out.add_field(FieldDef::new("merged_b", FieldType::Text));
        out.add_field(FieldDef::new("distance", FieldType::Float));
        out.add_field(FieldDef::new("size", FieldType::Integer));

        let mut steps: Vec<(String, String, f64)> = Vec::new();
        let mut step = 0usize;
        let mut first_merge = f64::NAN;

        while sigs.len() > 1 {
            let mut best: Option<(f64, usize, usize)> = None;
            for i in 0..sigs.len() {
                for j in (i + 1)..sigs.len() {
                    let d = signature_distance(&sigs[i], &sigs[j], distance);
                    if best.is_none_or(|(bd, _, _)| d < bd) {
                        best = Some((d, i, j));
                    }
                }
            }
            let (d, i, j) = best.expect("at least one pair while len > 1");
            if step == 0 {
                first_merge = d;
            }
            let a_label = sigs[i].labels.join("+");
            let b_label = sigs[j].labels.join("+");

            // Merge j into i (weighted mean), then drop j.
            let (ni, nj) = (sigs[i].count as f64, sigs[j].count as f64);
            let total = ni + nj;
            for d2 in 0..dim {
                sigs[i].mean[d2] = (sigs[i].mean[d2] * ni + sigs[j].mean[d2] * nj) / total;
            }
            sigs[i].spread = (sigs[i].spread * ni + sigs[j].spread * nj) / total;
            sigs[i].count += sigs[j].count;
            let mut jl = sigs[j].labels.clone();
            sigs[i].labels.append(&mut jl);
            let size = sigs[i].count;
            sigs.remove(j);

            step += 1;
            out.add_feature(
                None,
                &[
                    ("step", FieldValue::Integer(step as i64)),
                    ("merged_a", FieldValue::Text(a_label.clone())),
                    ("merged_b", FieldValue::Text(b_label.clone())),
                    ("distance", FieldValue::Float(d)),
                    ("size", FieldValue::Integer(size as i64)),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed writing merge row: {e}")))?;
            steps.push((a_label, b_label, d));
        }

        if let Some(path) = &text_out {
            crate::common::write_text_output(&render_text(&steps), path)?;
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("class_count".to_string(), json!(by_class.len()));
        outputs.insert("merge_count".to_string(), json!(step));
        // The smallest merge distance is the headline diagnostic: a near-zero
        // value means two classes were not actually separated.
        outputs.insert("min_merge_distance".to_string(), json!(first_merge));
        if let Some(p) = text_out {
            outputs.insert("output_text".to_string(), json!(p));
        }
        Ok(ToolRunResult { outputs })
    }
}

struct Sig {
    labels: Vec<String>,
    mean: Vec<f64>,
    spread: f64,
    count: usize,
}

fn signature_distance(a: &Sig, b: &Sig, mode: Distance) -> f64 {
    let euclid: f64 = a
        .mean
        .iter()
        .zip(b.mean.iter())
        .map(|(x, y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt();
    match mode {
        Distance::MeanOnly => euclid,
        Distance::Variance => {
            // Pooled spread; guard the zero case so identical tight classes do
            // not produce an infinite separation.
            let pooled = ((a.spread.powi(2) + b.spread.powi(2)) / 2.0).sqrt();
            if pooled <= f64::EPSILON {
                euclid
            } else {
                euclid / pooled
            }
        }
    }
}

/// Indented text rendering of the merge sequence, deepest merge last.
fn render_text(steps: &[(String, String, f64)]) -> String {
    let mut s = String::from("step  distance  merged\n");
    for (i, (a, b, d)) in steps.iter().enumerate() {
        s.push_str(&format!("{:>4}  {:>8.4}  {} + {}\n", i + 1, d, a, b));
    }
    s
}

fn field_string(v: &FieldValue) -> String {
    match v {
        FieldValue::Text(s) | FieldValue::Date(s) | FieldValue::DateTime(s) => s.clone(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(x) => x.to_string(),
        FieldValue::Boolean(b) => b.to_string(),
        FieldValue::Null | FieldValue::Blob(_) => String::new(),
    }
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    match args.get(key).and_then(Value::as_str) {
        Some(s) if !s.trim().is_empty() => Ok(s),
        _ => Err(ToolError::Validation(format!(
            "missing required string parameter '{key}'"
        ))),
    }
}

fn parse_fields(args: &ToolArgs) -> Result<Vec<String>, ToolError> {
    let raw = require_str(args, "fields")?;
    Ok(raw
        .split(',')
        .map(str::trim)
        .filter(|t| !t.is_empty())
        .map(String::from)
        .collect())
}

fn parse_distance(args: &ToolArgs) -> Result<Distance, ToolError> {
    match args.get("distance").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("variance") => Ok(Distance::Variance),
        Some("mean_only") => Ok(Distance::MeanOnly),
        Some(o) => Err(ToolError::Validation(format!(
            "'distance' must be variance or mean_only, got '{o}'"
        ))),
    }
}

fn parse_optional_bool(args: &ToolArgs, k: &str) -> Result<Option<bool>, ToolError> {
    match args.get(k) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{k}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{k}' must be a boolean"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Coord, Geometry, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Points tagged (class, v1, v2).
    fn layer_of(items: Vec<(&str, f64, f64)>) -> String {
        let mut l = Layer::new("c")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("cls", FieldType::Text));
        l.add_field(FieldDef::new("v1", FieldType::Float));
        l.add_field(FieldDef::new("v2", FieldType::Float));
        for (i, (c, a, b)) in items.into_iter().enumerate() {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(i as f64, 0.0))),
                &[
                    ("cls", FieldValue::Text(c.to_string())),
                    ("v1", FieldValue::Float(a)),
                    ("v2", FieldValue::Float(b)),
                ],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = DendrogramTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn merge_at(layer: &Layer, step: usize) -> (String, String, f64) {
        let si = layer.schema.field_index("step").unwrap();
        let ai = layer.schema.field_index("merged_a").unwrap();
        let bi = layer.schema.field_index("merged_b").unwrap();
        let di = layer.schema.field_index("distance").unwrap();
        for f in layer.features.iter() {
            if f.attributes[si].as_f64() == Some(step as f64) {
                return (
                    field_string(&f.attributes[ai]),
                    field_string(&f.attributes[bi]),
                    f.attributes[di].as_f64().unwrap(),
                );
            }
        }
        panic!("no step {step}");
    }

    /// THE diagnostic: two near-identical classes must merge first, and at a
    /// far smaller distance than the genuinely separate one.
    #[test]
    fn near_duplicate_classes_merge_first() {
        let input = layer_of(vec![
            ("a", 0.0, 0.0),
            ("a", 0.1, 0.1),
            ("b", 0.05, 0.05), // essentially the same cloud as a
            ("b", 0.15, 0.15),
            ("c", 100.0, 100.0), // far away
            ("c", 100.1, 100.1),
        ]);
        let (out, layer) = run(json!({
            "input": input, "class_field": "cls", "fields": "v1,v2"
        }));
        let (a, b, d1) = merge_at(&layer, 1);
        let mut merged = [a.as_str(), b.as_str()];
        merged.sort();
        assert_eq!(merged, ["a", "b"], "a and b are the duplicate pair");

        let (_, _, d2) = merge_at(&layer, 2);
        assert!(d1 < d2, "duplicate pair must merge closer than the outlier");
        assert_eq!(out.outputs["merge_count"], json!(2));
        assert_eq!(out.outputs["class_count"], json!(3));
    }

    /// min_merge_distance is the headline number and matches step 1.
    #[test]
    fn reports_min_merge_distance() {
        let input = layer_of(vec![("a", 0.0, 0.0), ("b", 1.0, 0.0), ("c", 50.0, 0.0)]);
        let (out, layer) = run(json!({
            "input": input, "class_field": "cls", "fields": "v1,v2",
            "distance": "mean_only", "standardize": false
        }));
        let (_, _, d1) = merge_at(&layer, 1);
        assert!((out.outputs["min_merge_distance"].as_f64().unwrap() - d1).abs() < 1e-9);
        assert!((d1 - 1.0).abs() < 1e-9, "a and b are 1 apart");
    }

    /// Standardization stops a large-range field from dominating the metric.
    ///
    /// Unstandardized, `v2`'s thousand-unit gap swamps everything and the merge
    /// distance carries that raw magnitude. Standardized, each field
    /// contributes on a unit-variance scale, so the distance collapses to
    /// order 1 and both fields actually participate.
    ///
    /// Note this does NOT assert the merge *order* flips: with only three
    /// classes, rescaling every dimension to unit variance makes the pairwise
    /// distances near-equal, so order is a knife-edge here and would make a
    /// brittle assertion. The scale change is the robust property.
    #[test]
    fn standardize_balances_field_ranges() {
        let input = layer_of(vec![
            ("a", 0.0, 0.0),
            ("b", 1000.0, 0.0),
            ("c", 1010.0, 1000.0),
        ]);

        let (raw, raw_layer) = run(json!({
            "input": input, "class_field": "cls", "fields": "v1,v2",
            "distance": "mean_only", "standardize": false
        }));
        // b and c sit 10 apart in v1 but 1000 apart in v2; unstandardized the
        // raw magnitudes decide and the distance is on v2's scale.
        let raw_d = raw.outputs["min_merge_distance"].as_f64().unwrap();
        assert!(
            raw_d > 100.0,
            "unstandardized distance should carry the raw field magnitude, got {raw_d}"
        );
        assert_eq!(raw_layer.len(), 2, "3 classes produce 2 merges");

        let (std, _) = run(json!({
            "input": input, "class_field": "cls", "fields": "v1,v2",
            "distance": "mean_only", "standardize": true
        }));
        let std_d = std.outputs["min_merge_distance"].as_f64().unwrap();
        assert!(
            std_d < 5.0,
            "standardized distance should be order 1, got {std_d}"
        );
    }

    /// Merged labels accumulate so the tree is readable.
    #[test]
    fn labels_accumulate_across_merges() {
        let input = layer_of(vec![("a", 0.0, 0.0), ("b", 0.1, 0.0), ("c", 50.0, 0.0)]);
        let (_, layer) = run(json!({
            "input": input, "class_field": "cls", "fields": "v1,v2"
        }));
        let (a, b, _) = merge_at(&layer, 2);
        // One side of the final merge is the combined a+b group.
        assert!(
            a.contains('+') || b.contains('+'),
            "expected a combined label, got '{a}' and '{b}'"
        );
    }

    /// Rows with a missing predictor are dropped, not imputed.
    #[test]
    fn rows_with_missing_values_are_dropped() {
        let mut l = Layer::new("c")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("cls", FieldType::Text));
        l.add_field(FieldDef::new("v1", FieldType::Float));
        for (c, v) in [("a", Some(0.0)), ("a", None), ("b", Some(10.0))] {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(0.0, 0.0))),
                &[
                    ("cls", FieldValue::Text(c.to_string())),
                    ("v1", v.map(FieldValue::Float).unwrap_or(FieldValue::Null)),
                ],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        let (out, _) = run(json!({
            "input": memory_store::make_vector_memory_path(&id),
            "class_field": "cls", "fields": "v1"
        }));
        assert_eq!(out.outputs["class_count"], json!(2));
    }

    #[test]
    fn rejects_bad_parameters() {
        let p = layer_of(vec![("a", 0.0, 0.0), ("b", 1.0, 1.0)]);
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            DendrogramTool.validate(&args).is_err()
        };
        assert!(bad(json!({ "class_field": "cls", "fields": "v1" })));
        assert!(bad(json!({ "input": p, "fields": "v1" })));
        assert!(bad(json!({ "input": p, "class_field": "cls" })));
        assert!(bad(
            json!({ "input": p, "class_field": "cls", "fields": "v1", "distance": "ward" })
        ));

        // A single class cannot form a tree.
        let one = layer_of(vec![("a", 0.0, 0.0), ("a", 1.0, 1.0)]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": one, "class_field": "cls", "fields": "v1" }))
                .unwrap();
        assert!(matches!(
            DendrogramTool.run(&args, &ctx()).unwrap_err(),
            ToolError::Validation(_)
        ));
    }
}
