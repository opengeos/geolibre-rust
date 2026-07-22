//! GeoLibre tool: seeded Monte-Carlo error propagation over an uncertain
//! attribute (ArcGIS Pro *Attribute Uncertainty*, Spatial Statistics).
//!
//! Each feature carries an estimate of some quantity (`value_field`) that is
//! only known up to an uncertainty. The uncertainty is described either by a
//! symmetric error field (`error_field` — a standard error / margin) or by an
//! explicit lower/upper bound pair (`lower_field` / `upper_field`). This tool
//! draws `iterations` seeded realizations of the value from that per-feature
//! error distribution (`NORMAL` or `UNIFORM`), then summarizes the realized
//! values with per-feature statistics: mean, standard deviation, the 5th /
//! 50th / 95th percentiles, and the coefficient of variation (a scale-free
//! stability score). The output is the input layer with those statistics
//! appended as new fields, so downstream analysis can carry the propagated
//! uncertainty forward.
//!
//! All randomness is a single seeded splitmix64 stream (no `Date::now` / `rand`
//! crate), so the output is bit-for-bit reproducible in native builds and in
//! WASM. The stream is advanced per feature in schema order, and each feature's
//! draws are independent, so the result is independent of iteration count only
//! in the limit — for a fixed seed and fixed `iterations` it is exactly
//! reproducible.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// z-score of a two-sided 95% confidence interval, used to convert an
/// `[lower, upper]` bound pair into a normal standard deviation.
const Z95: f64 = 1.959_963_984_540_054;

/// Seeded Monte-Carlo attribute error propagation.
pub struct AttributeUncertaintyTool;

impl Tool for AttributeUncertaintyTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "attribute_uncertainty",
            display_name: "Attribute Uncertainty",
            summary: "Seeded Monte-Carlo error propagation over an uncertain attribute (like ArcGIS Attribute Uncertainty): draw NORMAL/UNIFORM realizations from a per-feature error field or bound pair and append mean/std/p5/p50/p95/CV statistics.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (or table) with the uncertain attribute.",
                    required: true,
                },
                ToolParamSpec {
                    name: "value_field",
                    description: "Numeric field holding the estimate (distribution mean/center) for each feature.",
                    required: true,
                },
                ToolParamSpec {
                    name: "error_field",
                    description: "Numeric field holding the symmetric uncertainty: a standard error for NORMAL, or a half-width for UNIFORM. Provide this OR both lower_field/upper_field.",
                    required: false,
                },
                ToolParamSpec {
                    name: "lower_field",
                    description: "Numeric field with the lower bound of the value. Use together with upper_field instead of error_field.",
                    required: false,
                },
                ToolParamSpec {
                    name: "upper_field",
                    description: "Numeric field with the upper bound of the value. Use together with lower_field instead of error_field.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distribution",
                    description: "Error distribution: NORMAL (default) or UNIFORM.",
                    required: false,
                },
                ToolParamSpec {
                    name: "iterations",
                    description: "Number of Monte-Carlo realizations per feature (default 1000).",
                    required: false,
                },
                ToolParamSpec {
                    name: "seed",
                    description: "Seed for the deterministic RNG (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector layer with the appended statistics. If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "value_field")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let n = layer.len();

        let value_idx = field_index(&layer, &prm.value_field)?;
        let (spread_source, lo_idx, hi_idx, err_idx) = match &prm.spread {
            SpreadSpec::Error(f) => (SpreadKind::Error, None, None, Some(field_index(&layer, f)?)),
            SpreadSpec::Bounds { lower, upper } => (
                SpreadKind::Bounds,
                Some(field_index(&layer, lower)?),
                Some(field_index(&layer, upper)?),
                None,
            ),
        };

        ctx.progress.info(&format!(
            "attribute_uncertainty: {} feature(s), {} {} draw(s) each, seed {}",
            n,
            prm.iterations,
            prm.distribution.as_str(),
            prm.seed
        ));

        // Per-feature statistics columns (null where the value is missing).
        let mut mean = vec![FieldValue::Null; n];
        let mut std = vec![FieldValue::Null; n];
        let mut p5 = vec![FieldValue::Null; n];
        let mut p50 = vec![FieldValue::Null; n];
        let mut p95 = vec![FieldValue::Null; n];
        let mut cv = vec![FieldValue::Null; n];

        let mut draws: Vec<f64> = Vec::with_capacity(prm.iterations);
        let mut n_evaluated = 0usize;
        let mut n_skipped = 0usize;

        for fi in 0..n {
            // Re-seed per feature from a mix of the base seed and the feature
            // index so the stream (and thus each feature's result) is stable
            // regardless of how many draws earlier features consumed.
            let mut rng = Rng::new(prm.seed ^ splitmix(fi as u64 + 1));

            let center = read_f64(&layer, fi, value_idx);
            let sigma = match spread_source {
                SpreadKind::Error => read_f64(&layer, fi, err_idx.unwrap()),
                SpreadKind::Bounds => {
                    match (
                        read_f64(&layer, fi, lo_idx.unwrap()),
                        read_f64(&layer, fi, hi_idx.unwrap()),
                    ) {
                        (Some(lo), Some(hi)) => Some(hi - lo), // full width
                        _ => None,
                    }
                }
            };

            // center from value_field; for bounds, midpoint is used when the
            // value field itself is missing.
            let (mu, spread) = match (center, sigma) {
                (Some(c), Some(s)) => (c, s.abs()),
                (None, Some(s)) if spread_source == SpreadKind::Bounds => {
                    // Fall back to the midpoint of the bounds as the center.
                    let lo = read_f64(&layer, fi, lo_idx.unwrap());
                    let hi = read_f64(&layer, fi, hi_idx.unwrap());
                    match (lo, hi) {
                        (Some(lo), Some(hi)) => ((lo + hi) / 2.0, s.abs()),
                        _ => {
                            n_skipped += 1;
                            continue;
                        }
                    }
                }
                _ => {
                    n_skipped += 1;
                    continue;
                }
            };

            draws.clear();
            for _ in 0..prm.iterations {
                let v = match prm.distribution {
                    Distribution::Normal => {
                        // spread is a standard error for NORMAL; for bounds it
                        // is the full CI width, converted to sd via z95.
                        let sd = match spread_source {
                            SpreadKind::Error => spread,
                            SpreadKind::Bounds => spread / (2.0 * Z95),
                        };
                        mu + sd * rng.standard_normal()
                    }
                    Distribution::Uniform => {
                        // spread is a half-width for the error field, or the
                        // full width for bounds.
                        let half = match spread_source {
                            SpreadKind::Error => spread,
                            SpreadKind::Bounds => spread / 2.0,
                        };
                        mu + (2.0 * rng.f64() - 1.0) * half
                    }
                };
                draws.push(v);
            }

            let (m, s) = mean_sd(&draws);
            draws.sort_by(|a, b| a.total_cmp(b));
            mean[fi] = FieldValue::Float(m);
            std[fi] = FieldValue::Float(s);
            p5[fi] = FieldValue::Float(percentile(&draws, 5.0));
            p50[fi] = FieldValue::Float(percentile(&draws, 50.0));
            p95[fi] = FieldValue::Float(percentile(&draws, 95.0));
            cv[fi] = FieldValue::Float(if m.abs() > 0.0 { s / m.abs() } else { 0.0 });
            n_evaluated += 1;
        }

        // Append the statistics columns.
        let base = &prm.value_field;
        let columns: [(&str, Vec<FieldValue>); 6] = [
            ("mc_mean", mean),
            ("mc_std", std),
            ("mc_p5", p5),
            ("mc_p50", p50),
            ("mc_p95", p95),
            ("mc_cv", cv),
        ];
        let mut out_names = Vec::with_capacity(columns.len());
        for (suffix, _) in &columns {
            let name = format!("{base}_{suffix}");
            layer.add_field(FieldDef::new(name.clone(), FieldType::Float));
            out_names.push(name);
        }
        for (fi, feat) in layer.features.iter_mut().enumerate() {
            for (_, values) in &columns {
                feat.attributes.push(values[fi].clone());
            }
        }

        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("evaluated".to_string(), json!(n_evaluated));
        outputs.insert("skipped".to_string(), json!(n_skipped));
        outputs.insert("iterations".to_string(), json!(prm.iterations));
        outputs.insert("distribution".to_string(), json!(prm.distribution.as_str()));
        outputs.insert("seed".to_string(), json!(prm.seed));
        outputs.insert("output_fields".to_string(), json!(out_names));
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Distribution {
    Normal,
    Uniform,
}

impl Distribution {
    fn as_str(self) -> &'static str {
        match self {
            Distribution::Normal => "NORMAL",
            Distribution::Uniform => "UNIFORM",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum SpreadKind {
    Error,
    Bounds,
}

enum SpreadSpec {
    Error(String),
    Bounds { lower: String, upper: String },
}

struct Params {
    value_field: String,
    spread: SpreadSpec,
    distribution: Distribution,
    iterations: usize,
    seed: u64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let value_field = require_str(args, "value_field")?.to_string();

    let error_field = parse_optional_str(args, "error_field")?.map(str::to_string);
    let lower_field = parse_optional_str(args, "lower_field")?.map(str::to_string);
    let upper_field = parse_optional_str(args, "upper_field")?.map(str::to_string);

    let spread = match (error_field, lower_field, upper_field) {
        (Some(e), None, None) => SpreadSpec::Error(e),
        (None, Some(l), Some(u)) => SpreadSpec::Bounds { lower: l, upper: u },
        (Some(_), Some(_), _) | (Some(_), _, Some(_)) => {
            return Err(ToolError::Validation(
                "provide either 'error_field' OR both 'lower_field'/'upper_field', not both"
                    .to_string(),
            ));
        }
        (None, Some(_), None) | (None, None, Some(_)) => {
            return Err(ToolError::Validation(
                "'lower_field' and 'upper_field' must be provided together".to_string(),
            ));
        }
        (None, None, None) => {
            return Err(ToolError::Validation(
                "specify the uncertainty via 'error_field' or 'lower_field'+'upper_field'"
                    .to_string(),
            ));
        }
    };

    let distribution = match parse_optional_str(args, "distribution")? {
        None => Distribution::Normal,
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "normal" => Distribution::Normal,
            "uniform" => Distribution::Uniform,
            other => {
                return Err(ToolError::Validation(format!(
                    "unknown distribution '{other}' (expected NORMAL or UNIFORM)"
                )))
            }
        },
    };

    let iterations = parse_optional_u64(args, "iterations")?.unwrap_or(1000);
    if iterations < 2 {
        return Err(ToolError::Validation(
            "'iterations' must be at least 2".to_string(),
        ));
    }
    let iterations = iterations as usize;

    let seed = parse_optional_u64(args, "seed")?.unwrap_or(1);

    Ok(Params {
        value_field,
        spread,
        distribution,
        iterations,
        seed,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn parse_optional_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64().or_else(|| n.as_f64().map(|f| f.max(0.0) as u64))),
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

fn field_index(layer: &Layer, name: &str) -> Result<usize, ToolError> {
    layer
        .schema
        .field_index(name)
        .ok_or_else(|| ToolError::Validation(format!("field '{name}' not found")))
}

fn read_f64(layer: &Layer, fi: usize, idx: usize) -> Option<f64> {
    layer.features[fi]
        .attributes
        .get(idx)
        .and_then(FieldValue::as_f64)
        .filter(|v| v.is_finite())
}

// ── Statistics ────────────────────────────────────────────────────────────────

fn mean_sd(xs: &[f64]) -> (f64, f64) {
    let n = xs.len() as f64;
    if n == 0.0 {
        return (0.0, 0.0);
    }
    let mean = xs.iter().sum::<f64>() / n;
    let var = xs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
    (mean, var.max(0.0).sqrt())
}

/// Linear-interpolation percentile over a slice sorted ascending.
fn percentile(sorted: &[f64], pct: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    if sorted.len() == 1 {
        return sorted[0];
    }
    let rank = (pct / 100.0) * (sorted.len() - 1) as f64;
    let lo = rank.floor() as usize;
    let hi = rank.ceil() as usize;
    if lo == hi {
        return sorted[lo];
    }
    let frac = rank - lo as f64;
    sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

// ── Deterministic RNG (splitmix64) ──────────────────────────────────────────

fn splitmix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

struct Rng {
    state: u64,
    /// Cached second value of a Box-Muller pair.
    spare: Option<f64>,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Self {
            state: seed,
            spare: None,
        }
    }
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// Uniform in [0, 1).
    fn f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    /// One draw from the standard normal via Box-Muller (with a cached spare).
    fn standard_normal(&mut self) -> f64 {
        if let Some(v) = self.spare.take() {
            return v;
        }
        // Guard u1 away from 0 so ln() is finite.
        let mut u1 = self.f64();
        while u1 <= f64::MIN_POSITIVE {
            u1 = self.f64();
        }
        let u2 = self.f64();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = std::f64::consts::TAU * u2;
        self.spare = Some(r * theta.sin());
        r * theta.cos()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink, ToolContext};
    use wbvector::{memory_store, Coord, Feature, Geometry, GeometryType, Layer, Schema};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a point layer with a value field and an error field.
    fn build_layer(rows: &[(f64, f64)]) -> Layer {
        let mut schema = Schema::new();
        schema.add_field(FieldDef::new("val", FieldType::Float));
        schema.add_field(FieldDef::new("err", FieldType::Float));
        let mut features = Vec::new();
        for (i, (v, e)) in rows.iter().enumerate() {
            let mut f = Feature::new();
            f.fid = i as u64;
            f.geometry = Some(Geometry::Point(Coord::xy(i as f64, 0.0)));
            f.attributes = vec![FieldValue::Float(*v), FieldValue::Float(*e)];
            features.push(f);
        }
        Layer {
            name: "t".to_string(),
            geom_type: Some(GeometryType::Point),
            crs: None,
            schema,
            features,
            extent: None,
        }
    }

    /// Converts a JSON object into `ToolArgs` (a `BTreeMap`).
    fn to_args(v: Value) -> BTreeMap<String, Value> {
        v.as_object()
            .unwrap()
            .iter()
            .map(|(k, val)| (k.clone(), val.clone()))
            .collect()
    }

    fn run(layer: Layer, args: Value) -> Layer {
        let id = memory_store::put_vector(layer);
        let path = memory_store::make_vector_memory_path(&id);
        let mut map = to_args(args);
        map.insert("input".to_string(), json!(path));
        let res = AttributeUncertaintyTool.run(&map, &ctx()).unwrap();
        let out = res.outputs["output"].as_str().unwrap().to_string();
        let oid = memory_store::vector_path_to_id(&out).unwrap();
        (*memory_store::get_vector_arc_by_id(oid).unwrap()).clone()
    }

    fn col(layer: &Layer, name: &str) -> Vec<Option<f64>> {
        let idx = layer.schema.field_index(name).unwrap();
        layer
            .features
            .iter()
            .map(|f| f.attributes.get(idx).and_then(FieldValue::as_f64))
            .collect()
    }

    #[test]
    fn attribute_uncertainty_recovers_mean_and_std_normal() {
        // Large iteration count -> MC mean ~ value, MC std ~ error.
        let layer = build_layer(&[(100.0, 10.0)]);
        let out = run(
            layer,
            json!({ "value_field": "val", "error_field": "err",
                    "distribution": "NORMAL", "iterations": 20000, "seed": 7 }),
        );
        let mean = col(&out, "val_mc_mean")[0].unwrap();
        let std = col(&out, "val_mc_std")[0].unwrap();
        assert!((mean - 100.0).abs() < 0.5, "mean {mean} not near 100");
        assert!((std - 10.0).abs() < 0.5, "std {std} not near 10");
        // p50 near the mean, p5/p95 roughly +/- 1.645 sd.
        let p5 = col(&out, "val_mc_p5")[0].unwrap();
        let p95 = col(&out, "val_mc_p95")[0].unwrap();
        assert!((p5 - (100.0 - 1.645 * 10.0)).abs() < 1.0);
        assert!((p95 - (100.0 + 1.645 * 10.0)).abs() < 1.0);
    }

    #[test]
    fn attribute_uncertainty_zero_error_is_degenerate() {
        let layer = build_layer(&[(42.0, 0.0)]);
        let out = run(
            layer,
            json!({ "value_field": "val", "error_field": "err", "iterations": 1000, "seed": 3 }),
        );
        assert!((col(&out, "val_mc_mean")[0].unwrap() - 42.0).abs() < 1e-9);
        assert!(col(&out, "val_mc_std")[0].unwrap().abs() < 1e-9);
        assert!(col(&out, "val_mc_cv")[0].unwrap().abs() < 1e-9);
    }

    #[test]
    fn attribute_uncertainty_deterministic_by_seed() {
        let a = run(
            build_layer(&[(5.0, 2.0)]),
            json!({ "value_field": "val", "error_field": "err", "iterations": 500, "seed": 11 }),
        );
        let b = run(
            build_layer(&[(5.0, 2.0)]),
            json!({ "value_field": "val", "error_field": "err", "iterations": 500, "seed": 11 }),
        );
        assert_eq!(col(&a, "val_mc_mean"), col(&b, "val_mc_mean"));
        assert_eq!(col(&a, "val_mc_std"), col(&b, "val_mc_std"));
        assert_eq!(col(&a, "val_mc_p95"), col(&b, "val_mc_p95"));
    }

    #[test]
    fn attribute_uncertainty_uniform_bounds() {
        // Uniform in [lower, upper] -> mean ~ midpoint, all draws within bounds.
        let mut schema = Schema::new();
        schema.add_field(FieldDef::new("val", FieldType::Float));
        schema.add_field(FieldDef::new("lo", FieldType::Float));
        schema.add_field(FieldDef::new("hi", FieldType::Float));
        let mut f = Feature::new();
        f.geometry = Some(Geometry::Point(Coord::xy(0.0, 0.0)));
        f.attributes = vec![
            FieldValue::Float(10.0),
            FieldValue::Float(0.0),
            FieldValue::Float(20.0),
        ];
        let layer = Layer {
            name: "t".to_string(),
            geom_type: Some(GeometryType::Point),
            crs: None,
            schema,
            features: vec![f],
            extent: None,
        };
        let out = run(
            layer,
            json!({ "value_field": "val", "lower_field": "lo", "upper_field": "hi",
                    "distribution": "UNIFORM", "iterations": 20000, "seed": 5 }),
        );
        let mean = col(&out, "val_mc_mean")[0].unwrap();
        assert!(
            (mean - 10.0).abs() < 0.3,
            "mean {mean} not near midpoint 10"
        );
        let p5 = col(&out, "val_mc_p5")[0].unwrap();
        let p95 = col(&out, "val_mc_p95")[0].unwrap();
        assert!(p5 >= 0.0 && p95 <= 20.0, "percentiles outside bounds");
        // Uniform std ~ width / sqrt(12) = 20/3.4641 ~ 5.77.
        let std = col(&out, "val_mc_std")[0].unwrap();
        assert!((std - 5.7735).abs() < 0.2, "std {std} not near 5.77");
    }

    #[test]
    fn attribute_uncertainty_passes_through_missing_values() {
        let mut layer = build_layer(&[(100.0, 10.0), (0.0, 0.0)]);
        // Make the second feature's value null.
        layer.features[1].attributes[0] = FieldValue::Null;
        let out = run(
            layer,
            json!({ "value_field": "val", "error_field": "err", "iterations": 100, "seed": 1 }),
        );
        // First feature computed, second stays null.
        assert!(col(&out, "val_mc_mean")[0].is_some());
        assert!(col(&out, "val_mc_mean")[1].is_none());
    }

    #[test]
    fn attribute_uncertainty_rejects_bad_parameters() {
        // Missing value_field.
        assert!(AttributeUncertaintyTool
            .validate(&to_args(json!({ "input": "x", "error_field": "err" })))
            .is_err());
        // No spread specification at all.
        assert!(AttributeUncertaintyTool
            .validate(&to_args(json!({ "input": "x", "value_field": "val" })))
            .is_err());
        // Both error and bounds.
        assert!(AttributeUncertaintyTool
            .validate(&to_args(json!({ "input": "x", "value_field": "val",
                               "error_field": "err", "lower_field": "lo", "upper_field": "hi" })))
            .is_err());
        // Only one bound.
        assert!(AttributeUncertaintyTool
            .validate(&to_args(
                json!({ "input": "x", "value_field": "val", "lower_field": "lo" })
            ))
            .is_err());
        // Bad distribution.
        assert!(AttributeUncertaintyTool
            .validate(&to_args(json!({ "input": "x", "value_field": "val",
                               "error_field": "err", "distribution": "poisson" })))
            .is_err());
    }
}
