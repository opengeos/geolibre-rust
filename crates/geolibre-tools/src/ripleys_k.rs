//! GeoLibre tool: Ripley's K multi-distance point-pattern analysis.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Multi-Distance Spatial Cluster
//! Analysis (Ripley's K)* (Spatial Statistics): characterize whether a point
//! pattern is clustered or dispersed **across a range of distances**, which the
//! bundled single-scale statistics (`nearest_neighbour_index`,
//! `quadrat_count_test`) cannot.
//!
//! For each distance `d`, Ripley's K is
//!
//! ```text
//! K(d) = A * Σ_{i≠j} w_i w_j · 1[dist(i,j) ≤ d] / (Σw)² − Σw²)
//! ```
//!
//! with `A` the study-area (the bounding rectangle of the points). The
//! variance-stabilized L transform `L(d) = √(K(d)/π)` is compared to its CSR
//! expectation `d`: `L(d) − d > 0` means more points fall within `d` than
//! expected under complete spatial randomness (**clustering**); `< 0` means
//! **dispersion**.
//!
//! Significance comes from a Monte-Carlo envelope: `permutations` sets of `n`
//! uniform-random points are drawn in the same study area (a deterministic,
//! seeded splitmix64 RNG, so results are reproducible in WASM), and the min/max
//! `L(d)` across simulations give the lower/upper envelope. Observed `L(d)`
//! outside the envelope is significant clustering (above) or dispersion (below)
//! at that scale. Because the same uncorrected estimator is applied to both the
//! observed and the simulated patterns, edge effects cancel in the comparison —
//! so no explicit edge correction is applied.
//!
//! Output is a table (`distance, observed_k, expected_k, observed_l,
//! expected_l`, plus `lower_l, upper_l` when `permutations > 0`), written to a
//! CSV when `output` is given and always returned in the tool result. The pair
//! counting is O(n²) per realization, so this suits moderate point counts.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldValue, Geometry};

use crate::common::write_text_output;
use crate::vector_common::{load_input_layer, parse_optional_str};

pub struct RipleysKTool;

impl Tool for RipleysKTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "ripleys_k",
            display_name: "Ripley's K",
            summary: "Multi-distance point-pattern analysis (Ripley's K / L function) with Monte-Carlo complete-spatial-randomness envelopes, to detect clustering or dispersion across a range of distances.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer (other geometries use their vertex-mean representative point).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional CSV path for the K/L table. The table is also returned in the tool result.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distance_bands",
                    description: "Number of distance bands to evaluate (default 10).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_distance",
                    description: "Largest distance to evaluate, in CRS units. Default: a quarter of the shorter side of the point bounding box.",
                    required: false,
                },
                ToolParamSpec {
                    name: "permutations",
                    description: "Number of CSR simulations for the significance envelope (default 99; 0 disables the envelope).",
                    required: false,
                },
                ToolParamSpec {
                    name: "weight_field",
                    description: "Optional numeric field weighting each point.",
                    required: false,
                },
                ToolParamSpec {
                    name: "seed",
                    description: "Seed for the CSR RNG (default 1), for reproducible envelopes.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args
            .get("input")
            .and_then(Value::as_str)
            .map(str::trim)
            .unwrap_or("")
            .is_empty()
        {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args
            .get("input")
            .and_then(Value::as_str)
            .filter(|s| !s.trim().is_empty())
            .ok_or_else(|| {
                ToolError::Validation("missing required parameter 'input'".to_string())
            })?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let schema = &layer.schema;
        let mut pts: Vec<(f64, f64, f64)> = Vec::new(); // x, y, weight
        for feature in &layer.features {
            let Some((x, y)) = feature.geometry.as_ref().and_then(rep_point) else {
                continue;
            };
            let w = match &prm.weight_field {
                Some(f) => match feature.get(schema, f).ok().and_then(FieldValue::as_f64) {
                    Some(v) if v.is_finite() && v > 0.0 => v,
                    _ => continue,
                },
                None => 1.0,
            };
            pts.push((x, y, w));
        }
        let n = pts.len();
        if n < 2 {
            return Err(ToolError::Execution(
                "need at least 2 points for Ripley's K".to_string(),
            ));
        }

        let (minx, miny, maxx, maxy) = bbox(&pts);
        let (width, height) = (maxx - minx, maxy - miny);
        let area = width * height;
        if area <= 0.0 {
            return Err(ToolError::Execution(
                "points are collinear or coincident (study area has zero size)".to_string(),
            ));
        }
        let max_distance = match prm.max_distance {
            Some(d) => d,
            None => 0.25 * width.min(height),
        };
        if !(max_distance > 0.0 && max_distance.is_finite()) {
            return Err(ToolError::Validation(
                "could not determine a positive 'max_distance'; pass it explicitly".to_string(),
            ));
        }
        let step = max_distance / prm.distance_bands as f64;

        ctx.progress.info(&format!(
            "{n} point(s), {} band(s) to {max_distance:.3}, {} permutation(s)",
            prm.distance_bands, prm.permutations
        ));

        // Weight normalization: Σw² over ordered i≠j pairs = (Σw)² − Σw².
        let sum_w: f64 = pts.iter().map(|p| p.2).sum();
        let sum_w2: f64 = pts.iter().map(|p| p.2 * p.2).sum();
        let denom = sum_w * sum_w - sum_w2;
        let norm = if denom > 0.0 { area / denom } else { 0.0 };

        let observed_k = k_function(&pts, step, prm.distance_bands, norm);

        // CSR envelope.
        let mut lower = vec![f64::INFINITY; prm.distance_bands];
        let mut upper = vec![f64::NEG_INFINITY; prm.distance_bands];
        if prm.permutations > 0 {
            let weights: Vec<f64> = pts.iter().map(|p| p.2).collect();
            let mut rng = prm.seed;
            for _ in 0..prm.permutations {
                let sim: Vec<(f64, f64, f64)> = weights
                    .iter()
                    .map(|&w| {
                        (
                            minx + next_f64(&mut rng) * width,
                            miny + next_f64(&mut rng) * height,
                            w,
                        )
                    })
                    .collect();
                let k = k_function(&sim, step, prm.distance_bands, norm);
                for b in 0..prm.distance_bands {
                    let l = l_transform(k[b]);
                    lower[b] = lower[b].min(l);
                    upper[b] = upper[b].max(l);
                }
            }
        }

        // Assemble the table.
        let mut distances = Vec::with_capacity(prm.distance_bands);
        let mut obs_l = Vec::with_capacity(prm.distance_bands);
        let mut csv = String::from("distance,observed_k,expected_k,observed_l,expected_l");
        if prm.permutations > 0 {
            csv.push_str(",lower_l,upper_l");
        }
        csv.push('\n');
        let mut significant = 0usize;
        for b in 0..prm.distance_bands {
            let d = (b + 1) as f64 * step;
            let ok = observed_k[b];
            let ol = l_transform(ok);
            distances.push(d);
            obs_l.push(ol);
            csv.push_str(&format!(
                "{d:.6},{ok:.6},{:.6},{ol:.6},{d:.6}",
                std::f64::consts::PI * d * d
            ));
            if prm.permutations > 0 {
                csv.push_str(&format!(",{:.6},{:.6}", lower[b], upper[b]));
                if ol > upper[b] || ol < lower[b] {
                    significant += 1;
                }
            }
            csv.push('\n');
        }

        let out_path = match output {
            Some(path) => {
                write_text_output(&csv, path)?;
                Some(path.to_string())
            }
            None => None,
        };

        let mut outputs = BTreeMap::new();
        if let Some(p) = &out_path {
            outputs.insert("output".to_string(), json!(p));
        }
        outputs.insert("point_count".to_string(), json!(n));
        outputs.insert("study_area".to_string(), json!(area));
        outputs.insert("max_distance".to_string(), json!(max_distance));
        outputs.insert("distance_bands".to_string(), json!(prm.distance_bands));
        outputs.insert("permutations".to_string(), json!(prm.permutations));
        outputs.insert("distances".to_string(), json!(distances));
        outputs.insert("observed_l".to_string(), json!(obs_l));
        outputs.insert("expected_l".to_string(), json!(distances));
        if prm.permutations > 0 {
            outputs.insert("lower_l".to_string(), json!(lower));
            outputs.insert("upper_l".to_string(), json!(upper));
            outputs.insert("significant_bands".to_string(), json!(significant));
        }
        outputs.insert("table_csv".to_string(), json!(csv));
        Ok(ToolRunResult { outputs })
    }
}

// ── K function ────────────────────────────────────────────────────────────────

/// Cumulative Ripley's K per distance band. `norm = A / ((Σw)² − Σw²)`.
fn k_function(pts: &[(f64, f64, f64)], step: f64, bands: usize, norm: f64) -> Vec<f64> {
    let mut binned = vec![0.0f64; bands];
    let n = pts.len();
    for i in 0..n {
        let (xi, yi, wi) = pts[i];
        for &(xj, yj, wj) in pts.iter().skip(i + 1) {
            let d = ((xi - xj).powi(2) + (yi - yj).powi(2)).sqrt();
            let b = (d / step).floor() as isize;
            if b >= 0 && (b as usize) < bands {
                // Ordered pairs i≠j contribute both directions -> 2·w_i·w_j.
                binned[b as usize] += 2.0 * wi * wj;
            }
        }
    }
    // Prefix-sum to cumulative counts, then scale.
    let mut cum = 0.0;
    binned
        .iter()
        .map(|&c| {
            cum += c;
            cum * norm
        })
        .collect()
}

fn l_transform(k: f64) -> f64 {
    (k.max(0.0) / std::f64::consts::PI).sqrt()
}

fn bbox(pts: &[(f64, f64, f64)]) -> (f64, f64, f64, f64) {
    let mut minx = f64::INFINITY;
    let mut miny = f64::INFINITY;
    let mut maxx = f64::NEG_INFINITY;
    let mut maxy = f64::NEG_INFINITY;
    for &(x, y, _) in pts {
        minx = minx.min(x);
        miny = miny.min(y);
        maxx = maxx.max(x);
        maxy = maxy.max(y);
    }
    (minx, miny, maxx, maxy)
}

/// Representative point of a geometry: the coordinate for a point, otherwise the
/// mean of its vertices.
fn rep_point(geom: &Geometry) -> Option<(f64, f64)> {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0usize;
    collect(geom, &mut |x, y| {
        sx += x;
        sy += y;
        n += 1;
    });
    (n > 0).then(|| (sx / n as f64, sy / n as f64))
}

fn collect(geom: &Geometry, f: &mut impl FnMut(f64, f64)) {
    match geom {
        Geometry::Point(c) => f(c.x, c.y),
        Geometry::MultiPoint(cs) | Geometry::LineString(cs) => cs.iter().for_each(|c| f(c.x, c.y)),
        Geometry::MultiLineString(ls) => ls.iter().flatten().for_each(|c| f(c.x, c.y)),
        Geometry::Polygon { exterior, .. } => exterior.coords().iter().for_each(|c| f(c.x, c.y)),
        Geometry::MultiPolygon(parts) => parts
            .iter()
            .flat_map(|(e, _)| e.coords())
            .for_each(|c| f(c.x, c.y)),
        Geometry::GeometryCollection(gs) => gs.iter().for_each(|g| collect(g, f)),
    }
}

// ── Deterministic RNG (splitmix64) ────────────────────────────────────────────

fn next_u64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn next_f64(state: &mut u64) -> f64 {
    (next_u64(state) >> 11) as f64 / (1u64 << 53) as f64
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    distance_bands: usize,
    max_distance: Option<f64>,
    permutations: usize,
    weight_field: Option<String>,
    seed: u64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let distance_bands = match parse_optional_f64(args, "distance_bands")? {
        None => 10,
        Some(v) if v.fract() == 0.0 && (1.0..=1000.0).contains(&v) => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "parameter 'distance_bands' must be an integer between 1 and 1000".to_string(),
            ))
        }
    };
    let max_distance = parse_optional_f64(args, "max_distance")?;
    if let Some(d) = max_distance {
        if !(d > 0.0 && d.is_finite()) {
            return Err(ToolError::Validation(
                "parameter 'max_distance' must be a positive number".to_string(),
            ));
        }
    }
    let permutations = match parse_optional_f64(args, "permutations")? {
        None => 99,
        Some(v) if v.fract() == 0.0 && (0.0..=100_000.0).contains(&v) => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "parameter 'permutations' must be a non-negative integer".to_string(),
            ))
        }
    };
    let weight_field = parse_optional_str(args, "weight_field")?.map(str::to_string);
    let seed = match parse_optional_f64(args, "seed")? {
        None => 1,
        Some(v) if v.fract() == 0.0 && v >= 0.0 && v.is_finite() => v as u64,
        Some(_) => {
            return Err(ToolError::Validation(
                "parameter 'seed' must be a non-negative integer".to_string(),
            ))
        }
    };
    Ok(Params {
        distance_bands,
        max_distance,
        permutations,
        weight_field,
        seed,
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

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn points_layer(pts: &[(f64, f64)]) -> String {
        let mut layer = Layer::new("pts");
        for &(x, y) in pts {
            layer.add_feature(Some(Geometry::point(x, y)), &[]).unwrap();
        }
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        RipleysKTool.run(&args, &ctx()).unwrap()
    }

    fn floats(out: &ToolRunResult, key: &str) -> Vec<f64> {
        out.outputs[key]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap())
            .collect()
    }

    /// A regular grid is dispersed: at short distances observed L is below the
    /// CSR expectation (fewer close pairs than random).
    #[test]
    fn regular_grid_is_dispersed_at_short_range() {
        let mut pts = Vec::new();
        for i in 0..8 {
            for j in 0..8 {
                pts.push((i as f64 * 10.0, j as f64 * 10.0));
            }
        }
        let out = run(
            json!({ "input": points_layer(&pts), "max_distance": 8.0, "distance_bands": 4, "permutations": 0 }),
        );
        let d = floats(&out, "distances");
        let l = floats(&out, "observed_l");
        // At distances below the 10-unit spacing, no pairs -> L = 0 < d.
        assert!(
            l[0] < d[0],
            "grid should be dispersed: L({})={} !< {}",
            d[0],
            l[0],
            d[0]
        );
    }

    /// Two tight clusters are strongly clustered: at a distance spanning a
    /// cluster, observed L greatly exceeds the CSR expectation.
    #[test]
    fn clusters_are_detected_as_clustering() {
        let mut pts = Vec::new();
        // Cluster A around (0,0), cluster B around (100,100); 20 pts each.
        let mut s = 42u64;
        for _ in 0..20 {
            pts.push((next_f64(&mut s) * 2.0, next_f64(&mut s) * 2.0));
            pts.push((
                100.0 + next_f64(&mut s) * 2.0,
                100.0 + next_f64(&mut s) * 2.0,
            ));
        }
        let out = run(
            json!({ "input": points_layer(&pts), "max_distance": 20.0, "distance_bands": 10, "permutations": 0 }),
        );
        let d = floats(&out, "distances");
        let l = floats(&out, "observed_l");
        // At ~5 units (within a cluster of radius ~2) many pairs are close ->
        // strong clustering: observed L well above expected d.
        let b = 2; // distance ~ 6
        assert!(
            l[b] > d[b] * 2.0,
            "expected clustering: L({})={} vs d={}",
            d[b],
            l[b],
            d[b]
        );
    }

    /// The CSR envelope is deterministic: the same seed reproduces it exactly.
    #[test]
    fn envelope_is_reproducible_with_seed() {
        let pts: Vec<(f64, f64)> = (0..30)
            .map(|i| ((i * 7 % 11) as f64 * 9.0, (i * 13 % 11) as f64 * 9.0))
            .collect();
        let a = run(
            json!({ "input": points_layer(&pts), "permutations": 20, "seed": 7, "max_distance": 30.0 }),
        );
        let b = run(
            json!({ "input": points_layer(&pts), "permutations": 20, "seed": 7, "max_distance": 30.0 }),
        );
        assert_eq!(floats(&a, "upper_l"), floats(&b, "upper_l"));
        assert_eq!(floats(&a, "lower_l"), floats(&b, "lower_l"));
        // A different seed generally gives a different envelope.
        let c = run(
            json!({ "input": points_layer(&pts), "permutations": 20, "seed": 99, "max_distance": 30.0 }),
        );
        assert_ne!(floats(&a, "upper_l"), floats(&c, "upper_l"));
    }

    /// Envelope brackets the observed L for a roughly-random pattern, and the
    /// clustered pattern breaks above the envelope.
    #[test]
    fn envelope_flags_significant_clustering() {
        let mut pts = Vec::new();
        let mut s = 5u64;
        for _ in 0..40 {
            pts.push((next_f64(&mut s), next_f64(&mut s))); // one tight blob in [0,1]^2
        }
        // spread the study area so the blob is a strong cluster
        pts.push((50.0, 50.0));
        let out = run(
            json!({ "input": points_layer(&pts), "permutations": 50, "distance_bands": 8, "seed": 3 }),
        );
        assert!(out.outputs["significant_bands"].as_u64().unwrap() >= 1);
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = RipleysKTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "p.geojson", "distance_bands": 0 })).is_err());
        assert!(bad(json!({ "input": "p.geojson", "distance_bands": 2.5 })).is_err());
        assert!(bad(json!({ "input": "p.geojson", "max_distance": -1 })).is_err());
        assert!(bad(json!({ "input": "p.geojson", "permutations": -5 })).is_err());
        assert!(bad(json!({ "input": "p.geojson" })).is_ok());
    }
}
