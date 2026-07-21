//! GeoLibre tool: compute smoothed rates from count and population fields.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Calculate Rates* (Spatial
//! Statistics). Raw rates from small-population areas are unstable and dominate
//! downstream hot-spot analysis; empirical-Bayes smoothing shrinks each area's
//! rate toward a reference mean by an amount that depends on its population.
//! The bundled suite has the hot-spot statistics (`getis_ord_gi_star`,
//! `local_morans_i_lisa`) but nothing to stabilise the input rates — this is the
//! standard preparation step before them.
//!
//! Three methods:
//! * **crude** — `count / population * per`.
//! * **eb_global** — Marshall global empirical Bayes: shrink each crude rate
//!   toward the global mean rate, more for smaller populations.
//! * **eb_spatial** — the same shrinkage but toward a *local* reference rate
//!   computed over each area's k-nearest neighbours (by centroid), so smoothing
//!   respects spatial structure.
//!
//! Output fields (original attributes preserved): `crude_rate`, `smooth_rate`
//! (the EB estimate, equal to the crude rate for `method=crude`), and `rate_se`
//! (Poisson standard error of the crude rate).

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CalculateRatesTool;

impl Tool for CalculateRatesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "calculate_rates",
            display_name: "Calculate Rates",
            summary: "Compute smoothed rates from count and population fields (like ArcGIS Calculate Rates): crude rate, Marshall global empirical Bayes, or spatial empirical Bayes over k-nearest neighbours — the standard stabilization step before hot-spot analysis of disease/crime rates.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon or point layer with count and population fields.",
                    required: true,
                },
                ToolParamSpec {
                    name: "count_field",
                    description: "Field holding the event count (numerator).",
                    required: true,
                },
                ToolParamSpec {
                    name: "population_field",
                    description: "Field holding the population at risk (denominator).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output layer with rate fields added. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'crude', 'eb_global' (default), or 'eb_spatial'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "per",
                    description: "Rate multiplier, e.g. 100000 for per-100k (default 100000).",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighbors",
                    description: "eb_spatial only: number of nearest neighbours for the local reference rate (default 8).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "count_field")?;
        require_str(args, "population_field")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let n = layer.features.len();
        let count_idx = layer.schema.field_index(&prm.count_field).ok_or_else(|| {
            ToolError::Validation(format!("count_field '{}' not found", prm.count_field))
        })?;
        let pop_idx = layer
            .schema
            .field_index(&prm.population_field)
            .ok_or_else(|| {
                ToolError::Validation(format!(
                    "population_field '{}' not found",
                    prm.population_field
                ))
            })?;

        // Read counts, populations, and centroids.
        let mut o = vec![0.0f64; n];
        let mut p = vec![0.0f64; n];
        let mut cx = vec![0.0f64; n];
        let mut cy = vec![0.0f64; n];
        for (i, feat) in layer.features.iter().enumerate() {
            o[i] = feat
                .attributes
                .get(count_idx)
                .and_then(|v| v.as_f64())
                .unwrap_or(0.0);
            p[i] = feat
                .attributes
                .get(pop_idx)
                .and_then(|v| v.as_f64())
                .unwrap_or(0.0);
            if let Some(bb) = feat.geometry.as_ref().and_then(|g| g.bbox()) {
                cx[i] = (bb.min_x + bb.max_x) / 2.0;
                cy[i] = (bb.min_y + bb.max_y) / 2.0;
            }
        }

        // Crude rate (per-unit); guard against zero population.
        let crude: Vec<f64> = (0..n)
            .map(|i| if p[i] > 0.0 { o[i] / p[i] } else { 0.0 })
            .collect();

        ctx.progress
            .info(&format!("{} area(s), method {}", n, prm.method.label()));

        let smooth: Vec<f64> = match prm.method {
            Method::Crude => crude.clone(),
            Method::EbGlobal => {
                let idx: Vec<usize> = (0..n).collect();
                (0..n)
                    .map(|i| eb_estimate(i, &idx, &o, &p, &crude))
                    .collect()
            }
            Method::EbSpatial => (0..n)
                .map(|i| {
                    let mut nbrs = knn(i, &cx, &cy, prm.neighbors);
                    nbrs.push(i); // include self in the local set
                    eb_estimate(i, &nbrs, &o, &p, &crude)
                })
                .collect(),
        };

        // Poisson SE of the crude rate: sqrt(O)/P, scaled by `per`.
        let se: Vec<f64> = (0..n)
            .map(|i| if p[i] > 0.0 { o[i].sqrt() / p[i] } else { 0.0 })
            .collect();

        layer.add_field(FieldDef::new("crude_rate", FieldType::Float));
        layer.add_field(FieldDef::new("smooth_rate", FieldType::Float));
        layer.add_field(FieldDef::new("rate_se", FieldType::Float));
        for i in 0..n {
            let f = &mut layer.features[i];
            f.attributes.push(FieldValue::Float(crude[i] * prm.per));
            f.attributes.push(FieldValue::Float(smooth[i] * prm.per));
            f.attributes.push(FieldValue::Float(se[i] * prm.per));
        }

        let out_path = write_or_store_layer(layer, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        Ok(ToolRunResult { outputs })
    }
}

/// Marshall empirical-Bayes estimate for area `i`, using the reference set
/// `refset` (indices) to estimate the prior mean and between-area variance.
fn eb_estimate(i: usize, refset: &[usize], o: &[f64], p: &[f64], crude: &[f64]) -> f64 {
    let sum_o: f64 = refset.iter().map(|&j| o[j]).sum();
    let sum_p: f64 = refset.iter().map(|&j| p[j]).sum();
    if sum_p <= 0.0 || p[i] <= 0.0 {
        return crude[i];
    }
    let m = sum_o / sum_p; // reference mean rate
    let p_bar = sum_p / refset.len() as f64;
    // Population-weighted variance of the crude rates about m.
    let s2 = refset
        .iter()
        .map(|&j| p[j] * (crude[j] - m).powi(2))
        .sum::<f64>()
        / sum_p;
    let a = (s2 - m / p_bar).max(0.0); // between-area variance (non-negative)
    let c = a / (a + m / p[i]); // shrinkage weight in [0, 1]
    m + c * (crude[i] - m)
}

/// Indices of the `k` nearest neighbours of `i` by centroid (brute force).
fn knn(i: usize, cx: &[f64], cy: &[f64], k: usize) -> Vec<usize> {
    let mut d: Vec<(f64, usize)> = (0..cx.len())
        .filter(|&j| j != i)
        .map(|j| ((cx[j] - cx[i]).hypot(cy[j] - cy[i]), j))
        .collect();
    d.sort_by(|a, b| a.0.total_cmp(&b.0));
    d.into_iter().take(k).map(|(_, j)| j).collect()
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Crude,
    EbGlobal,
    EbSpatial,
}

impl Method {
    fn label(&self) -> &'static str {
        match self {
            Method::Crude => "crude",
            Method::EbGlobal => "eb_global",
            Method::EbSpatial => "eb_spatial",
        }
    }
}

struct Params {
    count_field: String,
    population_field: String,
    method: Method,
    per: f64,
    neighbors: usize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let count_field = require_str(args, "count_field")?.to_string();
    let population_field = require_str(args, "population_field")?.to_string();
    let method = match args.get("method").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("eb_global") => Method::EbGlobal,
        Some("crude") => Method::Crude,
        Some("eb_spatial") => Method::EbSpatial,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'method' must be crude/eb_global/eb_spatial, got '{o}'"
            )))
        }
    };
    let per = opt_f64(args, "per")?.unwrap_or(100_000.0);
    if per <= 0.0 {
        return Err(ToolError::Validation("'per' must be positive".to_string()));
    }
    let neighbors = opt_usize(args, "neighbors")?.unwrap_or(8).max(1);
    Ok(Params {
        count_field,
        population_field,
        method,
        per,
        neighbors,
    })
}

fn opt_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
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

fn opt_usize(args: &ToolArgs, key: &str) -> Result<Option<usize>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64().map(|v| v as usize)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
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
    use wbvector::{memory_store, Geometry, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn layer(rows: &[(f64, f64, f64, f64)]) -> String {
        let mut l = Layer::new("z")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("cnt", FieldType::Float));
        l.add_field(FieldDef::new("pop", FieldType::Float));
        for (x, y, c, p) in rows {
            l.add_feature(
                Some(Geometry::point(*x, *y)),
                &[("cnt", (*c).into()), ("pop", (*p).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> Layer {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CalculateRatesTool.run(&args, &ctx()).unwrap();
        load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    fn col(l: &Layer, name: &str) -> Vec<f64> {
        let i = l.schema.field_index(name).unwrap();
        l.features
            .iter()
            .map(|f| f.attributes[i].as_f64().unwrap())
            .collect()
    }

    /// Crude rate = count / population * per.
    #[test]
    fn crude_rate_is_exact() {
        let rows = [(0.0, 0.0, 5.0, 1000.0), (1.0, 0.0, 20.0, 2000.0)];
        let l = run(json!({
            "input": layer(&rows), "count_field": "cnt", "population_field": "pop",
            "method": "crude", "per": 100000.0,
        }));
        let cr = col(&l, "crude_rate");
        assert!((cr[0] - 500.0).abs() < 1e-6); // 5/1000*100k
        assert!((cr[1] - 1000.0).abs() < 1e-6); // 20/2000*100k
                                                // For crude method, smooth_rate == crude_rate.
        assert_eq!(col(&l, "smooth_rate"), cr);
    }

    /// Global EB shrinks a tiny-population area's extreme rate toward the mean
    /// more than a large-population area's rate.
    #[test]
    fn eb_global_shrinks_small_populations_more() {
        // One tiny area with a wildly high rate, many stable areas near the mean.
        let mut rows = vec![(0.0, 0.0, 5.0, 10.0)]; // rate 0.5 (huge), pop 10
        for i in 1..30 {
            rows.push((i as f64, 0.0, 10.0, 1000.0)); // rate 0.01, pop 1000
        }
        let l = run(json!({
            "input": layer(&rows), "count_field": "cnt", "population_field": "pop",
            "method": "eb_global", "per": 1.0,
        }));
        let crude = col(&l, "crude_rate");
        let smooth = col(&l, "smooth_rate");
        // The tiny area's smoothed rate is pulled far below its crude 0.5...
        assert!(
            smooth[0] < crude[0] * 0.8,
            "tiny-pop extreme rate must shrink toward mean"
        );
        // ...while a large-pop stable area barely moves.
        assert!(
            (smooth[5] - crude[5]).abs() < crude[5] * 0.2,
            "stable area barely shrinks"
        );
    }

    /// Spatial EB runs and keeps rates finite/non-negative.
    #[test]
    fn eb_spatial_runs() {
        let rows: Vec<(f64, f64, f64, f64)> = (0..20)
            .map(|i| {
                (
                    (i % 5) as f64,
                    (i / 5) as f64,
                    (i % 3 + 1) as f64,
                    100.0 + 10.0 * i as f64,
                )
            })
            .collect();
        let l = run(json!({
            "input": layer(&rows), "count_field": "cnt", "population_field": "pop",
            "method": "eb_spatial", "neighbors": 4, "per": 1000.0,
        }));
        for v in col(&l, "smooth_rate") {
            assert!(v.is_finite() && v >= 0.0);
        }
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CalculateRatesTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "count_field": "c" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "count_field": "c", "population_field": "p", "method": "x" })).is_err());
        assert!(bad(
            json!({ "input": "a.geojson", "count_field": "c", "population_field": "p", "per": -1 })
        )
        .is_err());
        assert!(
            bad(json!({ "input": "a.geojson", "count_field": "c", "population_field": "p" }))
                .is_ok()
        );
    }
}
