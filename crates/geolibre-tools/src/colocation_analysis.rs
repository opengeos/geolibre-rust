//! GeoLibre tool: local colocation quotient between two point categories.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Colocation Analysis* (Spatial
//! Statistics). The bundled point-pattern statistics (`ripleys_k`,
//! `nearest_neighbour_index`, `quadrat_count_test`) each treat one population at
//! a time; nothing measures the **asymmetric** association between two
//! categories — e.g. whether coffee shops (A) disproportionately sit among
//! transit stops (B).
//!
//! For each category-A point the tool computes the **local colocation quotient**
//! (Leslie & Kronenfeld 2011): the (kernel-weighted) fraction of its `neighbors`
//! nearest neighbours that are category B, divided by B's global share of the
//! other points. `CLQ > 1` means A is drawn toward B (colocation); `< 1` means A
//! avoids B. Significance comes from a **conditional permutation test**: with the
//! neighbour structure fixed, the category labels are reshuffled (a seeded
//! splitmix64 RNG, reproducible in WASM) and the local CLQ recomputed, giving a
//! one-sided p-value and a class (`colocated` / `isolated` / `none`).
//!
//! Output copies the input points and adds `local_clq`, `p_value`, and
//! `coloc_type`; non-A points get an empty class. A projected CRS is recommended
//! (distances and kernel bandwidth in its units).

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct ColocationAnalysisTool;

impl Tool for ColocationAnalysisTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "colocation_analysis",
            display_name: "Colocation Analysis",
            summary: "Local colocation quotient (Leslie & Kronenfeld) measuring whether category-A points are disproportionately found among category-B neighbours, with a seeded permutation-test p-value — like ArcGIS Colocation Analysis.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer with a category field.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (a copy of the input with colocation fields). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "category_field",
                    description: "Field holding each point's category.",
                    required: true,
                },
                ToolParamSpec {
                    name: "category_a",
                    description: "Focal category (A): the points whose colocation with B is measured.",
                    required: true,
                },
                ToolParamSpec {
                    name: "category_b",
                    description: "Neighbour category (B): the points A is tested against.",
                    required: true,
                },
                ToolParamSpec {
                    name: "neighbors",
                    description: "Number of nearest neighbours k (default 8).",
                    required: false,
                },
                ToolParamSpec {
                    name: "weight",
                    description: "Neighbour weighting: 'gaussian' (distance-decay, default) or 'uniform'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "permutations",
                    description: "Number of label permutations for the significance test (default 99; 0 disables).",
                    required: false,
                },
                ToolParamSpec {
                    name: "seed",
                    description: "Seed for the permutation RNG (default 1), for reproducible results.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "category_field")?;
        require_str(args, "category_a")?;
        require_str(args, "category_b")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let cat_idx = layer
            .schema
            .field_index(&prm.category_field)
            .ok_or_else(|| {
                ToolError::Validation(format!("category_field '{}' not found", prm.category_field))
            })?;

        // Collect points, categories, and the mapping back to features.
        let mut pts: Vec<(f64, f64)> = Vec::new();
        let mut is_b: Vec<bool> = Vec::new();
        let mut is_a: Vec<bool> = Vec::new();
        let mut feat_of: Vec<usize> = Vec::new();
        for (fi, feature) in layer.features.iter().enumerate() {
            let Some((x, y)) = feature.geometry.as_ref().and_then(point_xy) else {
                continue;
            };
            let cat = feature
                .attributes
                .get(cat_idx)
                .map(value_string)
                .unwrap_or_default();
            pts.push((x, y));
            is_a.push(cat == prm.category_a);
            is_b.push(cat == prm.category_b);
            feat_of.push(fi);
        }
        let n = pts.len();
        let n_b = is_b.iter().filter(|&&b| b).count();
        let n_a = is_a.iter().filter(|&&a| a).count();
        if n_a == 0 || n_b == 0 || n < 3 {
            return Err(ToolError::Execution(format!(
                "need points of both categories (A={n_a}, B={n_b}) and >= 3 points total"
            )));
        }
        let k = prm.neighbors.min(n - 1);
        let same_cat = prm.category_a == prm.category_b;
        // Expected B share among a point's potential neighbours.
        let expected = if same_cat {
            (n_b as f64 - 1.0) / (n as f64 - 1.0)
        } else {
            n_b as f64 / (n as f64 - 1.0)
        };

        ctx.progress.info(&format!(
            "{n_a} A-point(s), {n_b} B-point(s); k={k}, {} permutation(s)",
            prm.permutations
        ));

        // Precompute each A-point's k nearest neighbours and weights.
        struct Focal {
            point: usize,      // index in pts
            neigh: Vec<usize>, // neighbour indices
            w: Vec<f64>,       // neighbour weights
            wsum: f64,
        }
        let mut focals: Vec<Focal> = Vec::new();
        for i in 0..n {
            if !is_a[i] {
                continue;
            }
            let mut ds: Vec<(f64, usize)> = (0..n)
                .filter(|&j| j != i)
                .map(|j| (dist(pts[i], pts[j]), j))
                .collect();
            ds.sort_by(|a, b| a.0.total_cmp(&b.0));
            ds.truncate(k);
            let bandwidth = ds.last().map(|&(d, _)| d).unwrap_or(1.0).max(1e-9);
            let (neigh, w): (Vec<usize>, Vec<f64>) = ds
                .iter()
                .map(|&(d, j)| {
                    let weight = match prm.weight {
                        Weight::Uniform => 1.0,
                        Weight::Gaussian => (-0.5 * (d / bandwidth).powi(2)).exp(),
                    };
                    (j, weight)
                })
                .unzip();
            let wsum = w.iter().sum::<f64>().max(1e-12);
            focals.push(Focal {
                point: i,
                neigh,
                w,
                wsum,
            });
        }

        // Observed local CLQ.
        let clq = |f: &Focal, labels: &[bool]| -> f64 {
            let bw: f64 = f
                .neigh
                .iter()
                .zip(&f.w)
                .filter(|(&j, _)| labels[j])
                .map(|(_, &w)| w)
                .sum();
            (bw / f.wsum) / expected.max(1e-12)
        };
        let observed: Vec<f64> = focals.iter().map(|f| clq(f, &is_b)).collect();

        // Conditional permutation test: reshuffle labels, recompute.
        let mut ge = vec![1usize; focals.len()]; // count sim >= obs (incl. observed)
        let mut le = vec![1usize; focals.len()];
        if prm.permutations > 0 {
            let mut labels = is_b.clone();
            let mut rng = prm.seed;
            for _ in 0..prm.permutations {
                fisher_yates(&mut labels, &mut rng);
                for (fi, f) in focals.iter().enumerate() {
                    let sim = clq(f, &labels);
                    if sim >= observed[fi] {
                        ge[fi] += 1;
                    }
                    if sim <= observed[fi] {
                        le[fi] += 1;
                    }
                }
            }
        }

        // Write results onto a copy of the input.
        layer.add_field(FieldDef::new("local_clq", FieldType::Float));
        layer.add_field(FieldDef::new("p_value", FieldType::Float));
        layer.add_field(FieldDef::new("coloc_type", FieldType::Text));
        let denom = (prm.permutations + 1) as f64;
        let mut per_feat: Vec<(f64, f64, &str)> = vec![(f64::NAN, f64::NAN, ""); layer.len()];
        let mut colocated = 0usize;
        let mut isolated = 0usize;
        for (fi, f) in focals.iter().enumerate() {
            let obs = observed[fi];
            let (p, class) = if prm.permutations == 0 {
                (
                    f64::NAN,
                    if obs > 1.0 {
                        "colocated"
                    } else if obs < 1.0 {
                        "isolated"
                    } else {
                        "none"
                    },
                )
            } else {
                let p_high = ge[fi] as f64 / denom;
                let p_low = le[fi] as f64 / denom;
                if obs > 1.0 && p_high <= 0.05 {
                    colocated += 1;
                    (p_high, "colocated")
                } else if obs < 1.0 && p_low <= 0.05 {
                    isolated += 1;
                    (p_low, "isolated")
                } else {
                    (p_high.min(p_low), "none")
                }
            };
            per_feat[feat_of[f.point]] = (obs, p, class);
        }
        for (fi, feature) in layer.features.iter_mut().enumerate() {
            let (clq, p, class) = per_feat[fi];
            feature.attributes.push(FieldValue::Float(clq));
            feature.attributes.push(FieldValue::Float(p));
            feature.attributes.push(FieldValue::Text(class.to_string()));
        }

        let global_clq = observed.iter().sum::<f64>() / observed.len() as f64;
        ctx.progress.info(&format!(
            "global CLQ {global_clq:.3}; {colocated} colocated, {isolated} isolated A-point(s)"
        ));

        let feature_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("focal_count".to_string(), json!(focals.len()));
        outputs.insert("global_clq".to_string(), json!(global_clq));
        outputs.insert("colocated_count".to_string(), json!(colocated));
        outputs.insert("isolated_count".to_string(), json!(isolated));
        Ok(ToolRunResult { outputs })
    }
}

// ── Helpers ──────────────────────────────────────────────────────────────────

fn dist(a: (f64, f64), b: (f64, f64)) -> f64 {
    (a.0 - b.0).hypot(a.1 - b.1)
}

/// In-place Fisher-Yates shuffle with a seeded splitmix64 RNG.
fn fisher_yates<T>(v: &mut [T], rng: &mut u64) {
    let n = v.len();
    for i in (1..n).rev() {
        let j = (next_u64(rng) % (i as u64 + 1)) as usize;
        v.swap(i, j);
    }
}

fn next_u64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => {
            let n = cs.len() as f64;
            Some((
                cs.iter().map(|c| c.x).sum::<f64>() / n,
                cs.iter().map(|c| c.y).sum::<f64>() / n,
            ))
        }
        _ => None,
    }
}

fn value_string(fv: &FieldValue) -> String {
    if let Some(i) = fv.as_i64() {
        i.to_string()
    } else if let Some(f) = fv.as_f64() {
        format!("{f}")
    } else {
        fv.as_str().unwrap_or("").to_string()
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy)]
enum Weight {
    Uniform,
    Gaussian,
}

struct Params {
    category_field: String,
    category_a: String,
    category_b: String,
    neighbors: usize,
    weight: Weight,
    permutations: usize,
    seed: u64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let category_field = require_str(args, "category_field")?.to_string();
    let category_a = require_str(args, "category_a")?.to_string();
    let category_b = require_str(args, "category_b")?.to_string();
    let neighbors = match parse_opt_u64(args, "neighbors")? {
        None => 8,
        Some(v) if v >= 1 => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "'neighbors' must be >= 1".to_string(),
            ))
        }
    };
    let weight = match parse_optional_str(args, "weight")? {
        None => Weight::Gaussian,
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "gaussian" => Weight::Gaussian,
            "uniform" => Weight::Uniform,
            other => {
                return Err(ToolError::Validation(format!(
                    "'weight' must be 'gaussian' or 'uniform', got '{other}'"
                )))
            }
        },
    };
    let permutations = parse_opt_u64(args, "permutations")?.unwrap_or(99) as usize;
    let seed = parse_opt_u64(args, "seed")?.unwrap_or(1);
    Ok(Params {
        category_field,
        category_a,
        category_b,
        neighbors,
        weight,
        permutations,
        seed,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn parse_opt_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64()),
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
    use wbvector::{memory_store, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn layer_of(pts: &[(f64, f64, &str)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("cat", FieldType::Text));
        for &(x, y, c) in pts {
            l.add_feature(Some(Geometry::point(x, y)), &[("cat", c.into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        ColocationAnalysisTool.run(&args, &ctx()).unwrap()
    }

    /// A points placed right next to B points -> global CLQ well above 1.
    #[test]
    fn colocated_categories_give_high_clq() {
        let mut pts = Vec::new();
        // Each A is tightly surrounded by 4 B points, sites spread far apart,
        // against a large C background so B is not the whole map.
        for i in 0..10 {
            let x = i as f64 * 200.0;
            pts.push((x, 0.0, "A"));
            for j in 0..4 {
                pts.push((x + 0.3 * (j as f64 - 1.5), 0.3, "B"));
            }
        }
        for i in 0..80 {
            pts.push(((i % 10) as f64 * 200.0 + 50.0, 50.0 + i as f64, "C"));
        }
        let out = run(json!({
            "input": layer_of(&pts), "category_field": "cat",
            "category_a": "A", "category_b": "B", "neighbors": 4, "seed": 42,
        }));
        let g = out.outputs["global_clq"].as_f64().unwrap();
        assert!(g > 1.5, "colocated A/B should give CLQ >> 1, got {g}");
        assert!(
            out.outputs["colocated_count"].as_u64().unwrap() >= 3,
            "A points surrounded by B should be significantly colocated"
        );
    }

    /// A points kept away from B -> global CLQ below 1 (isolation).
    #[test]
    fn segregated_categories_give_low_clq() {
        let mut pts = Vec::new();
        // A cluster on the left, B cluster on the right.
        for i in 0..15 {
            pts.push((i as f64, 0.0, "A"));
            pts.push((1000.0 + i as f64, 0.0, "B"));
        }
        let out = run(json!({
            "input": layer_of(&pts), "category_field": "cat",
            "category_a": "A", "category_b": "B", "neighbors": 3, "seed": 1,
        }));
        let g = out.outputs["global_clq"].as_f64().unwrap();
        assert!(g < 0.5, "segregated A/B should give CLQ < 1, got {g}");
    }

    /// Deterministic: same seed -> identical global CLQ / counts.
    #[test]
    fn deterministic_with_seed() {
        let mut pts = Vec::new();
        for i in 0..20 {
            pts.push((
                i as f64 * 3.0,
                (i % 3) as f64,
                if i % 2 == 0 { "A" } else { "B" },
            ));
        }
        let args = json!({
            "input": layer_of(&pts), "category_field": "cat",
            "category_a": "A", "category_b": "B", "neighbors": 4, "seed": 7,
        });
        let a = run(args.clone());
        let b = run(args);
        assert_eq!(a.outputs["global_clq"], b.outputs["global_clq"]);
        assert_eq!(a.outputs["colocated_count"], b.outputs["colocated_count"]);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ColocationAnalysisTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "input": "a.geojson", "category_field": "c", "category_a": "A" })).is_err()
        );
        assert!(bad(json!({ "input": "a.geojson", "category_field": "c", "category_a": "A", "category_b": "B" })).is_ok());
    }
}
