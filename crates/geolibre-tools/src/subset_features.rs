//! GeoLibre tool: random train/test split of a feature layer.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Subset Features* (Data Management /
//! Geostatistical Analyst). The repo has grown a real modelling and validation
//! family — `classification_accuracy_assessment`,
//! `compute_accuracy_for_object_detection`, `presence_only_prediction`,
//! `forest_based_forecast`, `maximum_likelihood_classification`,
//! `generalized_linear_regression` — and every one of them assumes the caller
//! already holds out data. Nothing produced that split: the bundled
//! `random_sample` works on rasters and `lidar_classify_subset` is lidar-only.
//!
//! The shuffle is a partial Fisher-Yates driven by a seeded splitmix64 stream,
//! never `rand`/`Date::now`, so the same seed reproduces the same split in the
//! browser and natively. `group_field` stratifies, preserving each group's
//! share of the training set — the difference between an honest validation and
//! one that accidentally holds out an entire class.

use std::collections::BTreeMap;
use std::collections::HashMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Default seed. Fixed (not time-derived) so a plain run is reproducible.
const DEFAULT_SEED: u64 = 0x5EED_0DEF_A017_C0DE;

pub struct SubsetFeaturesTool;

impl Tool for SubsetFeaturesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "subset_features",
            display_name: "Subset Features",
            summary: "Randomly split a feature layer into a training subset and its complement (test subset), by percentage or absolute count, using a seeded deterministic RNG. Optionally stratified by a group field. Like ArcGIS Subset Features.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (any geometry type).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output_training",
                    description: "Output path for the training subset. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_test",
                    description: "Optional output path for the complement (test) subset. If omitted, the test subset is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "size",
                    description: "Training subset size: a percentage (0-100) or an absolute feature count, per 'size_method' (default 50).",
                    required: false,
                },
                ToolParamSpec {
                    name: "size_method",
                    description: "'percentage' (default) or 'absolute'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "seed",
                    description: "RNG seed for a reproducible split (default fixed constant).",
                    required: false,
                },
                ToolParamSpec {
                    name: "group_field",
                    description: "Optional field to stratify on: each group contributes its own proportional share of the training subset.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        let size = parse_optional_f64(args, "size")?.unwrap_or(50.0);
        let method = parse_method(args)?;
        if size < 0.0 {
            return Err(ToolError::Validation(
                "'size' must be non-negative".to_string(),
            ));
        }
        if method == SizeMethod::Percentage && size > 100.0 {
            return Err(ToolError::Validation(
                "'size' must be within 0-100 when 'size_method' is 'percentage'".to_string(),
            ));
        }
        parse_optional_u64(args, "seed")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let out_train = parse_optional_str(args, "output_training")?;
        let out_test = parse_optional_str(args, "output_test")?;
        let size = parse_optional_f64(args, "size")?.unwrap_or(50.0);
        let method = parse_method(args)?;
        let seed = parse_optional_u64(args, "seed")?.unwrap_or(DEFAULT_SEED);
        let group_field = parse_optional_str(args, "group_field")?;

        let layer = load_input_layer(input)?;
        let n = layer.features.len();
        if n == 0 {
            return Err(ToolError::Execution("input has no features".to_string()));
        }
        let gidx = match group_field {
            Some(f) => Some(
                layer
                    .schema
                    .field_index(f)
                    .ok_or_else(|| ToolError::Validation(format!("group_field '{f}' not found")))?,
            ),
            None => None,
        };

        // Bucket feature indices: one bucket overall, or one per group value.
        let mut buckets: Vec<Vec<usize>> = Vec::new();
        match gidx {
            None => buckets.push((0..n).collect()),
            Some(g) => {
                let mut pos: HashMap<String, usize> = HashMap::new();
                for (i, f) in layer.features.iter().enumerate() {
                    let key = key_of(f.attributes.get(g));
                    let b = *pos.entry(key).or_insert_with(|| {
                        buckets.push(Vec::new());
                        buckets.len() - 1
                    });
                    buckets[b].push(i);
                }
            }
        }

        // Partial Fisher-Yates per bucket, taking that bucket's share.
        let mut is_training = vec![false; n];
        let mut rng = Rng::new(seed);
        for bucket in &mut buckets {
            let bn = bucket.len();
            let take = match method {
                SizeMethod::Percentage => {
                    ((size / 100.0) * bn as f64).round().clamp(0.0, bn as f64) as usize
                }
                SizeMethod::Absolute => {
                    // Absolute counts are per-group when stratifying, scaled by
                    // the group's share so the total still lands on `size`.
                    let want = if gidx.is_some() {
                        (size * bn as f64 / n as f64).round()
                    } else {
                        size.round()
                    };
                    want.clamp(0.0, bn as f64) as usize
                }
            };
            for k in 0..take {
                let j = k + (rng.next_u64() as usize) % (bn - k);
                bucket.swap(k, j);
                is_training[bucket[k]] = true;
            }
        }

        let n_train = is_training.iter().filter(|&&t| t).count();
        ctx.progress.info(&format!(
            "split {n} feature(s) into {n_train} training / {} test",
            n - n_train
        ));

        let training = build_subset(&layer, &is_training, true, "training")?;
        let test = build_subset(&layer, &is_training, false, "test")?;
        let train_path = write_or_store_layer(training, out_train)?;
        let test_path = write_or_store_layer(test, out_test)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output_training".to_string(), json!(train_path));
        outputs.insert("output_test".to_string(), json!(test_path));
        // 'output' mirrors the training set so generic callers that read a
        // single 'output' key still get the primary result.
        outputs.insert("output".to_string(), json!(train_path));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("training_count".to_string(), json!(n_train));
        outputs.insert("test_count".to_string(), json!(n - n_train));
        outputs.insert("group_count".to_string(), json!(buckets.len()));
        Ok(ToolRunResult { outputs })
    }
}

/// Copies the features whose `is_training` flag equals `want` into a new layer,
/// preserving all attributes and appending `ORIG_FID`.
fn build_subset(
    layer: &Layer,
    is_training: &[bool],
    want: bool,
    name: &str,
) -> Result<Layer, ToolError> {
    let mut out = Layer::new(name);
    if let Some(gt) = layer.geom_type {
        out = out.with_geom_type(gt);
    }
    if let Some(epsg) = layer.crs_epsg() {
        out = out.with_crs_epsg(epsg);
    }
    for fd in layer.schema.fields() {
        out.add_field(fd.clone());
    }
    out.add_field(FieldDef::new("ORIG_FID", FieldType::Integer));

    let names: Vec<String> = layer
        .schema
        .fields()
        .iter()
        .map(|f| f.name.clone())
        .collect();
    for (i, feat) in layer.features.iter().enumerate() {
        if is_training[i] != want {
            continue;
        }
        let mut attrs: Vec<(&str, FieldValue)> = names
            .iter()
            .enumerate()
            .map(|(fi, nm)| {
                (
                    nm.as_str(),
                    feat.attributes.get(fi).cloned().unwrap_or(FieldValue::Null),
                )
            })
            .collect();
        attrs.push(("ORIG_FID", FieldValue::Integer(i as i64)));
        out.add_feature(feat.geometry.clone(), &attrs)
            .map_err(|e| ToolError::Execution(format!("failed adding feature: {e}")))?;
    }
    Ok(out)
}

fn key_of(v: Option<&FieldValue>) -> String {
    match v {
        None | Some(FieldValue::Null) => "NULL".to_string(),
        Some(FieldValue::Integer(i)) => i.to_string(),
        Some(FieldValue::Float(f)) => format!("{f}"),
        Some(FieldValue::Text(s)) => s.clone(),
        Some(FieldValue::Boolean(b)) => b.to_string(),
        Some(FieldValue::Date(s)) | Some(FieldValue::DateTime(s)) => s.clone(),
        Some(FieldValue::Blob(b)) => format!("blob[{}]", b.len()),
    }
}

// ── Deterministic RNG (splitmix64) ──────────────────────────────────────────

struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

// ── Params ──────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum SizeMethod {
    Percentage,
    Absolute,
}

fn parse_method(args: &ToolArgs) -> Result<SizeMethod, ToolError> {
    match args
        .get("size_method")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        None | Some("") | Some("percentage") | Some("percentage_of_input") => {
            Ok(SizeMethod::Percentage)
        }
        Some("absolute") | Some("absolute_value") => Ok(SizeMethod::Absolute),
        Some(o) => Err(ToolError::Validation(format!(
            "'size_method' must be 'percentage' or 'absolute', got '{o}'"
        ))),
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
    Ok(parse_optional_f64(args, key)?.map(|v| v.abs() as u64))
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
    use wbvector::{Geometry, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn points(n: usize, groups: &[&str]) -> String {
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        let grouped = !groups.is_empty();
        if grouped {
            l.add_field(FieldDef::new("grp", FieldType::Text));
        }
        for i in 0..n {
            let a: Vec<(&str, FieldValue)> = if grouped {
                vec![("grp", FieldValue::Text(groups[i % groups.len()].to_string()))]
            } else {
                vec![]
            };
            l.add_feature(Some(Geometry::point(i as f64, 0.0)), &a)
                .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SubsetFeaturesTool.run(&args, &ctx()).unwrap();
        let tr = load_input_layer(out.outputs["output_training"].as_str().unwrap()).unwrap();
        let te = load_input_layer(out.outputs["output_test"].as_str().unwrap()).unwrap();
        (out, tr, te)
    }

    #[test]
    fn percentage_split_partitions_without_overlap() {
        let input = points(100, &[]);
        let (out, tr, te) = run(json!({ "input": input, "size": 70 }));
        assert_eq!(out.outputs["training_count"], json!(70));
        assert_eq!(out.outputs["test_count"], json!(30));
        // The two subsets must be a true partition: disjoint and complete.
        let idx = tr.schema.field_index("ORIG_FID").unwrap();
        let mut all: Vec<i64> = tr
            .iter()
            .chain(te.iter())
            .map(|f| f.attributes[idx].as_i64().unwrap())
            .collect();
        all.sort_unstable();
        all.dedup();
        assert_eq!(all.len(), 100);
    }

    #[test]
    fn same_seed_reproduces_the_split() {
        let input = points(50, &[]);
        let ids = |l: &Layer| -> Vec<i64> {
            let i = l.schema.field_index("ORIG_FID").unwrap();
            l.iter().map(|f| f.attributes[i].as_i64().unwrap()).collect()
        };
        let (_o, a, _t) = run(json!({ "input": input.clone(), "size": 40, "seed": 7 }));
        let (_o, b, _t) = run(json!({ "input": input.clone(), "size": 40, "seed": 7 }));
        let (_o, c, _t) = run(json!({ "input": input, "size": 40, "seed": 8 }));
        assert_eq!(ids(&a), ids(&b));
        assert_ne!(ids(&a), ids(&c));
    }

    #[test]
    fn absolute_size_takes_exactly_that_many() {
        let input = points(30, &[]);
        let (out, _tr, _te) = run(json!({
            "input": input, "size": 12, "size_method": "absolute"
        }));
        assert_eq!(out.outputs["training_count"], json!(12));
        assert_eq!(out.outputs["test_count"], json!(18));
    }

    #[test]
    fn stratified_split_preserves_group_shares() {
        // 60 features across 3 groups (20 each); 50% must take 10 from each,
        // not 30 from whichever group the shuffle happened to favour.
        let input = points(60, &["a", "b", "c"]);
        let (out, tr, _te) = run(json!({
            "input": input, "size": 50, "group_field": "grp"
        }));
        assert_eq!(out.outputs["group_count"], json!(3));
        assert_eq!(out.outputs["training_count"], json!(30));
        let g = tr.schema.field_index("grp").unwrap();
        let mut counts: BTreeMap<String, usize> = BTreeMap::new();
        for f in tr.iter() {
            *counts
                .entry(f.attributes[g].as_str().unwrap().to_string())
                .or_default() += 1;
        }
        assert_eq!(counts.values().copied().collect::<Vec<_>>(), vec![10, 10, 10]);
    }

    #[test]
    fn oversized_absolute_request_is_clamped() {
        let input = points(10, &[]);
        let (out, _tr, te) = run(json!({
            "input": input, "size": 999, "size_method": "absolute"
        }));
        assert_eq!(out.outputs["training_count"], json!(10));
        assert_eq!(te.features.len(), 0);
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SubsetFeaturesTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "p.shp", "size": 150 })).is_err());
        assert!(bad(json!({ "input": "p.shp", "size": -1 })).is_err());
        assert!(bad(json!({ "input": "p.shp", "size_method": "bogus" })).is_err());
        assert!(bad(json!({ "input": "p.shp", "size": 150, "size_method": "absolute" })).is_ok());
    }
}
