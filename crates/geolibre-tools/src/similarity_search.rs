//! GeoLibre tool: rank candidate features by attribute similarity.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Similarity Search* (Spatial
//! Statistics). Nothing bundled does attribute-space similarity ranking — the
//! classifiers predict classes; this ranks. "Find the census tracts most like
//! this one" for site selection and market analysis; a natural companion to
//! `geographically_weighted_regression` / `build_balanced_zones`.
//!
//! Each chosen numeric field is **z-standardized** over the combined
//! reference+candidate distribution, so fields with different units contribute
//! equally. The reference profile is the mean of the reference features'
//! standardized vectors. Every candidate is scored against it — Euclidean
//! distance (smaller = more similar) or cosine similarity (larger = more
//! similar) — and ranked. `most_or_least` keeps the most-similar, least-similar,
//! or both ends; `num_results` caps how many. Output copies the candidates with
//! a `sim_rank`, `sim_score`, `match` label, and per-field standardized
//! differences (`<field>_zd`).

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Feature, FieldDef, FieldType, FieldValue, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct SimilaritySearchTool;

impl Tool for SimilaritySearchTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "similarity_search",
            display_name: "Similarity Search",
            summary: "Rank candidate features by attribute similarity to one or more reference features across z-standardized numeric fields (Euclidean or cosine), most/least/both — like ArcGIS Similarity Search.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "reference",
                    description: "Reference vector layer (one or more features defining the profile to match).",
                    required: true,
                },
                ToolParamSpec {
                    name: "candidates",
                    description: "Candidate vector layer to rank by similarity to the reference.",
                    required: true,
                },
                ToolParamSpec {
                    name: "fields",
                    description: "Comma-separated numeric field(s) to compare (must exist on both layers).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (ranked candidates). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "match_method",
                    description: "'euclidean' (z-standardized distance; default) or 'cosine' (angular similarity).",
                    required: false,
                },
                ToolParamSpec {
                    name: "most_or_least",
                    description: "Keep the 'most' similar (default), 'least' similar, or 'both' ends.",
                    required: false,
                },
                ToolParamSpec {
                    name: "num_results",
                    description: "Number of results to keep per end (0 = all candidates). Default 10.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "reference")?;
        require_str(args, "candidates")?;
        require_str(args, "fields")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let reference = require_str(args, "reference")?;
        let candidates = require_str(args, "candidates")?;
        let fields_arg = require_str(args, "fields")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let ref_layer = load_input_layer(reference)?;
        let cand_layer = load_input_layer(candidates)?;

        let field_names: Vec<String> = fields_arg
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect();
        if field_names.is_empty() {
            return Err(ToolError::Validation(
                "'fields' must name at least one field".to_string(),
            ));
        }
        let ref_idx = field_indices(&ref_layer, &field_names, "reference")?;
        let cand_idx = field_indices(&cand_layer, &field_names, "candidates")?;
        let nf = field_names.len();

        // Raw field vectors for reference and candidate features.
        let ref_vecs = collect_vectors(&ref_layer, &ref_idx);
        let cand_vecs = collect_vectors(&cand_layer, &cand_idx);
        if ref_vecs.is_empty() {
            return Err(ToolError::Execution(
                "no reference features with all fields present".to_string(),
            ));
        }
        if cand_vecs.is_empty() {
            return Err(ToolError::Execution(
                "no candidate features with all fields present".to_string(),
            ));
        }

        // Standardize each field over the combined distribution.
        let mut mean = vec![0.0; nf];
        let mut std = vec![1.0; nf];
        let combined: Vec<&Vec<f64>> = ref_vecs
            .iter()
            .map(|(_, v)| v)
            .chain(cand_vecs.iter().map(|(_, v)| v))
            .collect();
        let m = combined.len() as f64;
        for f in 0..nf {
            let mu = combined.iter().map(|v| v[f]).sum::<f64>() / m;
            let var = combined.iter().map(|v| (v[f] - mu).powi(2)).sum::<f64>() / m;
            mean[f] = mu;
            std[f] = var.sqrt().max(1e-12);
        }
        let z = |v: &[f64]| -> Vec<f64> { (0..nf).map(|f| (v[f] - mean[f]) / std[f]).collect() };

        // Reference profile = mean of the reference z-vectors.
        let mut profile = vec![0.0; nf];
        for (_, v) in &ref_vecs {
            let zv = z(v);
            for f in 0..nf {
                profile[f] += zv[f];
            }
        }
        for p in profile.iter_mut() {
            *p /= ref_vecs.len() as f64;
        }

        ctx.progress.info(&format!(
            "ranking {} candidate(s) against a {}-feature reference profile ({} field(s))",
            cand_vecs.len(),
            ref_vecs.len(),
            nf
        ));

        // Score each candidate.
        struct Scored {
            feat: usize,
            score: f64, // distance (euclidean) or -cosine so smaller = more similar
            zd: Vec<f64>,
        }
        let mut scored: Vec<Scored> = cand_vecs
            .iter()
            .map(|(fi, v)| {
                let zv = z(v);
                let zd: Vec<f64> = (0..nf).map(|f| zv[f] - profile[f]).collect();
                let score = match prm.method {
                    Method::Euclidean => zd.iter().map(|d| d * d).sum::<f64>().sqrt(),
                    Method::Cosine => -cosine(&zv, &profile),
                };
                Scored {
                    feat: *fi,
                    score,
                    zd,
                }
            })
            .collect();
        // Sort most-similar first (smallest score).
        scored.sort_by(|a, b| a.score.total_cmp(&b.score));

        // Select which candidates to emit and their rank/label.
        let total = scored.len();
        let cap = if prm.num_results == 0 {
            total
        } else {
            prm.num_results.min(total)
        };
        // (index into scored, rank, label)
        let mut chosen: Vec<(usize, i64, &'static str)> = Vec::new();
        if matches!(prm.select, Select::Most | Select::Both) {
            for (r, s) in scored.iter().enumerate().take(cap) {
                chosen.push((r, (r + 1) as i64, "most"));
                let _ = s;
            }
        }
        if matches!(prm.select, Select::Least | Select::Both) {
            for (r, _s) in scored.iter().rev().enumerate().take(cap) {
                let idx = total - 1 - r;
                chosen.push((idx, (r + 1) as i64, "least"));
            }
        }

        // Build the output layer: candidate schema + similarity columns.
        let mut out = Layer::new("similarity");
        out.geom_type = cand_layer.geom_type;
        out.schema = cand_layer.schema.clone();
        if let Some(epsg) = cand_layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("sim_rank", FieldType::Integer));
        out.add_field(FieldDef::new("sim_score", FieldType::Float));
        out.add_field(FieldDef::new("match", FieldType::Text));
        for name in &field_names {
            out.add_field(FieldDef::new(format!("{name}_zd"), FieldType::Float));
        }

        for (si, rank, label) in &chosen {
            let s = &scored[*si];
            let src = &cand_layer.features[s.feat];
            let mut attrs = src.attributes.clone();
            attrs.push(FieldValue::Integer(*rank));
            let report_score = match prm.method {
                Method::Euclidean => s.score,
                Method::Cosine => -s.score, // report cosine similarity itself
            };
            attrs.push(FieldValue::Float(report_score));
            attrs.push(FieldValue::Text((*label).to_string()));
            for d in &s.zd {
                attrs.push(FieldValue::Float(*d));
            }
            out.push(Feature {
                fid: 0,
                geometry: src.geometry.clone(),
                attributes: attrs,
            });
        }

        let feature_count = out.len();
        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("candidate_count".to_string(), json!(cand_vecs.len()));
        outputs.insert("result_count".to_string(), json!(feature_count));
        outputs.insert("field_count".to_string(), json!(nf));
        Ok(ToolRunResult { outputs })
    }
}

// ── Helpers ──────────────────────────────────────────────────────────────────

fn cosine(a: &[f64], b: &[f64]) -> f64 {
    let dot: f64 = a.iter().zip(b).map(|(x, y)| x * y).sum();
    let na: f64 = a.iter().map(|x| x * x).sum::<f64>().sqrt();
    let nb: f64 = b.iter().map(|x| x * x).sum::<f64>().sqrt();
    if na < 1e-12 || nb < 1e-12 {
        0.0
    } else {
        dot / (na * nb)
    }
}

fn field_indices(layer: &Layer, names: &[String], which: &str) -> Result<Vec<usize>, ToolError> {
    names
        .iter()
        .map(|n| {
            layer.schema.field_index(n).ok_or_else(|| {
                ToolError::Validation(format!("field '{n}' not found on {which} layer"))
            })
        })
        .collect()
}

/// Collects (feature_index, field values) for features where every field parses.
fn collect_vectors(layer: &Layer, idx: &[usize]) -> Vec<(usize, Vec<f64>)> {
    let mut out = Vec::new();
    for (fi, feature) in layer.features.iter().enumerate() {
        let mut v = Vec::with_capacity(idx.len());
        let mut ok = true;
        for &i in idx {
            match feature.attributes.get(i).and_then(FieldValue::as_f64) {
                Some(x) if x.is_finite() => v.push(x),
                _ => {
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            out.push((fi, v));
        }
    }
    out
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy)]
enum Method {
    Euclidean,
    Cosine,
}

#[derive(Clone, Copy)]
enum Select {
    Most,
    Least,
    Both,
}

struct Params {
    method: Method,
    select: Select,
    num_results: usize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let method = match parse_optional_str(args, "match_method")? {
        None => Method::Euclidean,
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "euclidean" => Method::Euclidean,
            "cosine" => Method::Cosine,
            other => {
                return Err(ToolError::Validation(format!(
                    "'match_method' must be 'euclidean' or 'cosine', got '{other}'"
                )))
            }
        },
    };
    let select = match parse_optional_str(args, "most_or_least")? {
        None => Select::Most,
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "most" => Select::Most,
            "least" => Select::Least,
            "both" => Select::Both,
            other => {
                return Err(ToolError::Validation(format!(
                    "'most_or_least' must be 'most', 'least', or 'both', got '{other}'"
                )))
            }
        },
    };
    let num_results = match args.get("num_results") {
        None | Some(Value::Null) => 10,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(10) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 10,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'num_results' must be an integer".into()))?,
        Some(_) => {
            return Err(ToolError::Validation(
                "'num_results' must be an integer".to_string(),
            ))
        }
    };
    Ok(Params {
        method,
        select,
        num_results,
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
    use wbvector::{memory_store, Geometry, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn layer(rows: &[(&str, f64, f64)]) -> String {
        let mut l = Layer::new("x")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.add_field(FieldDef::new("a", FieldType::Float));
        l.add_field(FieldDef::new("b", FieldType::Float));
        for (i, (n, a, b)) in rows.iter().enumerate() {
            l.add_feature(
                Some(Geometry::point(i as f64, 0.0)),
                &[
                    ("name", (*n).into()),
                    ("a", (*a).into()),
                    ("b", (*b).into()),
                ],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SimilaritySearchTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn names_in_order(layer: &Layer) -> Vec<String> {
        let ni = layer.schema.field_index("name").unwrap();
        let ri = layer.schema.field_index("sim_rank").unwrap();
        let mut v: Vec<(i64, String)> = layer
            .iter()
            .map(|f| {
                (
                    f.attributes[ri].as_i64().unwrap(),
                    f.attributes[ni].as_str().unwrap().to_string(),
                )
            })
            .collect();
        v.sort_by_key(|(r, _)| *r);
        v.into_iter().map(|(_, n)| n).collect()
    }

    /// The candidate matching the reference ranks first; the opposite ranks last.
    #[test]
    fn ranks_most_similar_first() {
        let reference = layer(&[("ref", 10.0, 10.0)]);
        let candidates = layer(&[
            ("twin", 10.0, 10.1), // almost identical
            ("near", 9.0, 11.0),  // close
            ("far", -50.0, 80.0), // very different
        ]);
        let (_o, out) = run(json!({
            "reference": reference, "candidates": candidates, "fields": "a,b",
            "num_results": 0,
        }));
        let order = names_in_order(&out);
        assert_eq!(
            order.first().map(String::as_str),
            Some("twin"),
            "twin should rank most similar"
        );
        assert_eq!(
            order.last().map(String::as_str),
            Some("far"),
            "far should rank least similar"
        );
    }

    /// num_results caps the output.
    #[test]
    fn num_results_caps_output() {
        let reference = layer(&[("ref", 0.0, 0.0)]);
        let candidates = layer(&[
            ("a", 1.0, 1.0),
            ("b", 2.0, 2.0),
            ("c", 3.0, 3.0),
            ("d", 4.0, 4.0),
        ]);
        let (out, _l) = run(json!({
            "reference": reference, "candidates": candidates, "fields": "a,b", "num_results": 2,
        }));
        assert_eq!(out.outputs["result_count"], json!(2));
    }

    /// 'both' returns the top and bottom ends, labelled.
    #[test]
    fn both_ends() {
        let reference = layer(&[("ref", 0.0, 0.0)]);
        let candidates = layer(&[("close", 0.1, 0.1), ("mid", 5.0, 5.0), ("far", 20.0, 20.0)]);
        let (_o, out) = run(json!({
            "reference": reference, "candidates": candidates, "fields": "a,b",
            "most_or_least": "both", "num_results": 1,
        }));
        let mi = out.schema.field_index("match").unwrap();
        let ni = out.schema.field_index("name").unwrap();
        let mut most = None;
        let mut least = None;
        for f in out.iter() {
            let label = f.attributes[mi].as_str().unwrap();
            let name = f.attributes[ni].as_str().unwrap().to_string();
            if label == "most" {
                most = Some(name);
            } else {
                least = Some(name);
            }
        }
        assert_eq!(most.as_deref(), Some("close"));
        assert_eq!(least.as_deref(), Some("far"));
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SimilaritySearchTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "reference": "r.geojson", "candidates": "c.geojson" })).is_err());
        assert!(bad(json!({ "reference": "r.geojson", "candidates": "c.geojson", "fields": "a", "match_method": "manhattan" })).is_err());
        assert!(
            bad(json!({ "reference": "r.geojson", "candidates": "c.geojson", "fields": "a" }))
                .is_ok()
        );
    }
}
