//! GeoLibre tool: statistically compare two hot-spot result surfaces.
//!
//! Pure-Rust counterpart of ArcGIS Spatial Statistics' *Hot Spot Analysis
//! Comparison*. The repo already computes hot spots (`optimized_hot_spot_analysis`,
//! `emerging_hot_spot_analysis`, bundled `getis_ord_gi_star`); this tool takes
//! two such result layers, matches their locations, and classifies where the
//! hot/cold pattern changed between them.
//!
//! Each input carries a per-feature significance value in `bin_field` (positive =
//! hot, negative = cold, |value| below `significance` = not significant).
//! Features are joined by an id (`match_field`) or, failing that, by nearest
//! centroid within `tolerance`. Every matched location is labeled with a
//! transition category (e.g. `hot_to_hot`, `hot_to_cold`, `none_to_hot`) and an
//! agreement flag. A seeded permutation test reports whether the two patterns
//! agree more than expected by chance.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct HotSpotAnalysisComparisonTool;

impl Tool for HotSpotAnalysisComparisonTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "hot_spot_analysis_comparison",
            display_name: "Hot Spot Analysis Comparison",
            summary: "Compare two Getis-Ord hot-spot result layers: match their locations, classify each transition (hot_to_hot, hot_to_cold, none_to_hot, ...), and run a seeded permutation test of overall agreement — like ArcGIS Hot Spot Analysis Comparison.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input1",
                    description: "First hot-spot result layer (points or polygons).",
                    required: true,
                },
                ToolParamSpec {
                    name: "input2",
                    description: "Second hot-spot result layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "bin_field",
                    description: "Significance field present on both inputs (positive = hot, negative = cold).",
                    required: true,
                },
                ToolParamSpec {
                    name: "significance",
                    description: "Absolute-value threshold below which a location is 'none' (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "match_field",
                    description: "Optional id field present on both inputs used to join locations; otherwise nearest centroid.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Maximum centroid distance for nearest-location matching (map units). Unlimited if omitted.",
                    required: false,
                },
                ToolParamSpec {
                    name: "permutations",
                    description: "Number of permutations for the global agreement test (default 499).",
                    required: false,
                },
                ToolParamSpec {
                    name: "seed",
                    description: "Seed for the deterministic permutation RNG (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input1")?;
        require_str(args, "input2")?;
        require_str(args, "bin_field")?;
        parse_optional_f64(args, "significance")?;
        parse_optional_f64(args, "tolerance")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input1 = require_str(args, "input1")?;
        let input2 = require_str(args, "input2")?;
        let bin_field = require_str(args, "bin_field")?;
        let significance = parse_optional_f64(args, "significance")?
            .unwrap_or(0.0)
            .abs();
        let match_field = parse_optional_str(args, "match_field")?;
        let tolerance = parse_optional_f64(args, "tolerance")?;
        let permutations = parse_optional_usize(args, "permutations")?.unwrap_or(499);
        let seed = parse_optional_u64(args, "seed")?.unwrap_or(1);
        let output = parse_optional_str(args, "output")?;

        let l1 = load_input_layer(input1)?;
        let l2 = load_input_layer(input2)?;
        let b1 = l1.schema.field_index(bin_field).ok_or_else(|| {
            ToolError::Validation(format!("bin_field '{bin_field}' not found in input1"))
        })?;
        let b2 = l2.schema.field_index(bin_field).ok_or_else(|| {
            ToolError::Validation(format!("bin_field '{bin_field}' not found in input2"))
        })?;

        // Feature centroids + bin values.
        let f1: Vec<(Coord, f64, Option<Geometry>, Option<String>)> =
            features(&l1, b1, match_field);
        let f2: Vec<(Coord, f64, Option<Geometry>, Option<String>)> =
            features(&l2, b2, match_field);

        // Match each input1 feature to an input2 feature.
        let mut matched: Vec<Match> = Vec::new();
        for (c1, v1, g1, k1) in &f1 {
            let mut best: Option<(usize, f64)> = None;
            for (idx, (c2, _v2, _g2, k2)) in f2.iter().enumerate() {
                let ok = match (match_field, k1, k2) {
                    (Some(_), Some(a), Some(b)) => a == b,
                    (Some(_), _, _) => false,
                    (None, _, _) => true,
                };
                if !ok {
                    continue;
                }
                let d = (c1.x - c2.x).hypot(c1.y - c2.y);
                if best.map(|(_, bd)| d < bd).unwrap_or(true) {
                    best = Some((idx, d));
                }
                if match_field.is_some() {
                    break; // first key match is the join
                }
            }
            if let Some((idx, d)) = best {
                if tolerance.map(|t| d <= t).unwrap_or(true) {
                    matched.push(Match {
                        geom: g1.clone(),
                        cat1: classify(*v1, significance),
                        cat2: classify(f2[idx].1, significance),
                        bin1: *v1,
                        bin2: f2[idx].1,
                    });
                }
            }
        }

        if matched.is_empty() {
            return Err(ToolError::Execution(
                "no locations matched between the two inputs".to_string(),
            ));
        }

        let observed = matched.iter().filter(|m| m.cat1 == m.cat2).count();
        let agreement_fraction = observed as f64 / matched.len() as f64;

        // Permutation test: shuffle cat2 labels, count agreement >= observed.
        let cats2: Vec<Category> = matched.iter().map(|m| m.cat2).collect();
        let cats1: Vec<Category> = matched.iter().map(|m| m.cat1).collect();
        let mut rng = SplitMix64::new(seed);
        let mut ge = 0usize;
        for _ in 0..permutations {
            let mut perm = cats2.clone();
            shuffle(&mut perm, &mut rng);
            let agree = cats1
                .iter()
                .zip(perm.iter())
                .filter(|(a, b)| a == b)
                .count();
            if agree >= observed {
                ge += 1;
            }
        }
        let p_value = (ge as f64 + 1.0) / (permutations as f64 + 1.0);

        ctx.progress.info(&format!(
            "{} matched, agreement {:.3}, p={:.4}",
            matched.len(),
            agreement_fraction,
            p_value
        ));

        // Output layer.
        let mut out = Layer::new("hot_spot_comparison");
        if let Some(gt) = l1.geom_type {
            out = out.with_geom_type(gt);
        } else {
            out = out.with_geom_type(GeometryType::Point);
        }
        if let Some(epsg) = l1.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("cat1", FieldType::Text));
        out.add_field(FieldDef::new("cat2", FieldType::Text));
        out.add_field(FieldDef::new("category", FieldType::Text));
        out.add_field(FieldDef::new("agreement", FieldType::Integer));
        out.add_field(FieldDef::new("bin1", FieldType::Float));
        out.add_field(FieldDef::new("bin2", FieldType::Float));

        let mut transitions: BTreeMap<String, u64> = BTreeMap::new();
        for m in &matched {
            let category = format!("{}_to_{}", m.cat1.label(), m.cat2.label());
            *transitions.entry(category.clone()).or_insert(0) += 1;
            out.add_feature(
                m.geom.clone(),
                &[
                    ("cat1", FieldValue::Text(m.cat1.label().to_string())),
                    ("cat2", FieldValue::Text(m.cat2.label().to_string())),
                    ("category", FieldValue::Text(category)),
                    (
                        "agreement",
                        FieldValue::Integer(if m.cat1 == m.cat2 { 1 } else { 0 }),
                    ),
                    ("bin1", FieldValue::Float(m.bin1)),
                    ("bin2", FieldValue::Float(m.bin2)),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed adding feature: {e}")))?;
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("n_matched".to_string(), json!(matched.len()));
        outputs.insert("agreement_fraction".to_string(), json!(agreement_fraction));
        outputs.insert("p_value".to_string(), json!(p_value));
        outputs.insert(
            "transitions".to_string(),
            json!(transitions
                .into_iter()
                .map(|(k, v)| json!({ "category": k, "count": v }))
                .collect::<Vec<_>>()),
        );
        Ok(ToolRunResult { outputs })
    }
}

struct Match {
    geom: Option<Geometry>,
    cat1: Category,
    cat2: Category,
    bin1: f64,
    bin2: f64,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Category {
    Hot,
    Cold,
    None,
}

impl Category {
    fn label(self) -> &'static str {
        match self {
            Category::Hot => "hot",
            Category::Cold => "cold",
            Category::None => "none",
        }
    }
}

fn classify(v: f64, significance: f64) -> Category {
    if !v.is_finite() || v.abs() <= significance {
        Category::None
    } else if v > 0.0 {
        Category::Hot
    } else {
        Category::Cold
    }
}

/// Extracts (centroid, bin value, geometry, optional match key) per feature.
fn features(
    layer: &Layer,
    bin_idx: usize,
    match_field: Option<&str>,
) -> Vec<(Coord, f64, Option<Geometry>, Option<String>)> {
    let kidx = match_field.and_then(|f| layer.schema.field_index(f));
    layer
        .features
        .iter()
        .map(|f| {
            let c = f
                .geometry
                .as_ref()
                .map(centroid)
                .unwrap_or(Coord::xy(f64::NAN, f64::NAN));
            let v = f
                .attributes
                .get(bin_idx)
                .and_then(FieldValue::as_f64)
                .unwrap_or(f64::NAN);
            let key = kidx.and_then(|k| f.attributes.get(k)).map(|fv| match fv {
                FieldValue::Text(s) | FieldValue::Date(s) | FieldValue::DateTime(s) => s.clone(),
                other => other.to_string(),
            });
            (c, v, f.geometry.clone(), key)
        })
        .collect()
}

fn centroid(geom: &Geometry) -> Coord {
    let cs = geom.all_coords();
    if cs.is_empty() {
        return Coord::xy(f64::NAN, f64::NAN);
    }
    let (mut sx, mut sy) = (0.0, 0.0);
    for c in &cs {
        sx += c.x;
        sy += c.y;
    }
    let n = cs.len() as f64;
    Coord::xy(sx / n, sy / n)
}

// ── Deterministic RNG ───────────────────────────────────────────────────────

struct SplitMix64(u64);

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        SplitMix64(seed.wrapping_add(0x9E3779B97F4A7C15))
    }
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

fn shuffle<T>(v: &mut [T], rng: &mut SplitMix64) {
    for i in (1..v.len()).rev() {
        let j = rng.below(i + 1);
        v.swap(i, j);
    }
}

// ── Params ──────────────────────────────────────────────────────────────────

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

fn parse_optional_usize(args: &ToolArgs, key: &str) -> Result<Option<usize>, ToolError> {
    Ok(parse_optional_f64(args, key)?.and_then(|v| {
        if v.fract() == 0.0 && v >= 0.0 {
            Some(v as usize)
        } else {
            None
        }
    }))
}

fn parse_optional_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    Ok(parse_optional_f64(args, key)?.and_then(|v| {
        if v.fract() == 0.0 && v >= 0.0 {
            Some(v as u64)
        } else {
            None
        }
    }))
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

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Points with a `gi` significance field at matching coordinates.
    fn layer(pts: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("h")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("gi", FieldType::Float));
        for (x, y, g) in pts {
            l.add_feature(
                Some(Geometry::point(*x, *y)),
                &[("gi", FieldValue::Float(*g))],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = HotSpotAnalysisComparisonTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    #[test]
    fn identical_patterns_agree_fully() {
        let a = layer(&[(0.0, 0.0, 2.0), (1.0, 0.0, -2.0), (2.0, 0.0, 0.0)]);
        let out = run(json!({ "input1": a.clone(), "input2": a, "bin_field": "gi" })).0;
        assert_eq!(out.outputs["n_matched"], json!(3));
        assert!((out.outputs["agreement_fraction"].as_f64().unwrap() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn transitions_are_labeled() {
        // hot->cold at (0,0); cold->cold at (1,0).
        let a = layer(&[(0.0, 0.0, 3.0), (1.0, 0.0, -3.0)]);
        let b = layer(&[(0.0, 0.0, -3.0), (1.0, 0.0, -3.0)]);
        let (out, l) = run(json!({ "input1": a, "input2": b, "bin_field": "gi" }));
        let cats: Vec<String> = l
            .iter()
            .map(|f| f.get(&l.schema, "category").unwrap().to_string())
            .collect();
        assert!(cats.contains(&"hot_to_cold".to_string()));
        assert!(cats.contains(&"cold_to_cold".to_string()));
        assert!(out.outputs["agreement_fraction"].as_f64().unwrap() < 1.0);
    }

    #[test]
    fn significance_threshold_demotes_to_none() {
        // |gi|=1 below threshold 1.5 -> none.
        let a = layer(&[(0.0, 0.0, 1.0)]);
        let b = layer(&[(0.0, 0.0, 1.0)]);
        let (_o, l) =
            run(json!({ "input1": a, "input2": b, "bin_field": "gi", "significance": 1.5 }));
        assert_eq!(
            l.features[0].get(&l.schema, "cat1").unwrap(),
            &FieldValue::Text("none".into())
        );
    }

    #[test]
    fn deterministic_p_value() {
        let a = layer(&[(0.0, 0.0, 2.0), (1.0, 0.0, -2.0), (2.0, 0.0, 2.0)]);
        let b = layer(&[(0.0, 0.0, 2.0), (1.0, 0.0, 2.0), (2.0, 0.0, -2.0)]);
        let p1 =
            run(json!({ "input1": a.clone(), "input2": b.clone(), "bin_field": "gi", "seed": 7 }))
                .0
                .outputs["p_value"]
                .as_f64()
                .unwrap();
        let p2 = run(json!({ "input1": a, "input2": b, "bin_field": "gi", "seed": 7 }))
            .0
            .outputs["p_value"]
            .as_f64()
            .unwrap();
        assert_eq!(p1, p2);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            HotSpotAnalysisComparisonTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "input1": "a", "input2": "b" })).is_err(),
            "needs bin_field"
        );
        assert!(bad(json!({ "input1": "a", "input2": "b", "bin_field": "gi" })).is_ok());
    }
}
