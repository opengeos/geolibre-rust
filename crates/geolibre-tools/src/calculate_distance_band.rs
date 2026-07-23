//! GeoLibre tool: distance band from a neighbour count.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Calculate Distance Band from Neighbor
//! Count* (Spatial Statistics). It reports the minimum, average, and maximum
//! distance at which every feature has at least `neighbors` neighbours — the
//! diagnostic used to choose a distance-band threshold before running any
//! distance-based spatial statistic (`getis_ord_general_g`, incremental
//! Moran's I, distance-band spatial weights, …).
//!
//! For each feature the distance to its `neighbors`-th nearest neighbour is
//! computed with a kd-tree; the three summary statistics are the min / mean /
//! max of those per-feature distances. Using the **maximum** as the band
//! guarantees no feature is left with fewer than `neighbors` neighbours.
//!
//! `distance_method` is `euclidean` (default) or `manhattan`. An optional
//! `output` CSV writes the per-feature k-th-neighbour distance. Use a projected
//! CRS so distances are in its linear units.

use std::collections::BTreeMap;

use kdtree::KdTree;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::Geometry;

use crate::common::write_text_output;
use crate::vector_common::{load_input_layer, parse_optional_str};

pub struct CalculateDistanceBandTool;

impl Tool for CalculateDistanceBandTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "calculate_distance_band",
            display_name: "Calculate Distance Band",
            summary: "Report the minimum, average, and maximum distance at which every feature has at least N neighbours (like ArcGIS Calculate Distance Band from Neighbor Count) — the threshold diagnostic that guarantees no feature is left with zero neighbours before distance-based spatial stats. The bundled generate_spatial_weights_matrix builds the matrix but never reports this band.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point/feature vector layer (projected CRS; distances in its units).",
                    required: true,
                },
                ToolParamSpec {
                    name: "neighbors",
                    description: "Neighbour count k (>= 1). Default 1.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distance_method",
                    description: "'euclidean' (default) or 'manhattan'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional CSV path for a feature_id,kth_neighbor_distance table.",
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
        let prm = parse_params(args)?;
        let output = parse_optional_str(args, "output")?;

        let layer = load_input_layer(input)?;

        // Representative point per feature; remember which feature it came from.
        // Non-point/multipoint geometries have no representative location and are
        // skipped; the count is surfaced so the band isn't silently computed over
        // a subset (which would mislead threshold selection).
        let mut pts: Vec<[f64; 2]> = Vec::new();
        let mut feat_of: Vec<usize> = Vec::new();
        let mut skipped_non_point = 0usize;
        for (fi, feature) in layer.features.iter().enumerate() {
            if let Some((x, y)) = feature.geometry.as_ref().and_then(point_xy) {
                pts.push([x, y]);
                feat_of.push(fi);
            } else {
                skipped_non_point += 1;
            }
        }
        if skipped_non_point > 0 {
            ctx.progress.info(&format!(
                "skipped {skipped_non_point} non-point feature(s) with no representative location"
            ));
        }
        let n = pts.len();
        if n <= prm.neighbors {
            return Err(ToolError::Execution(format!(
                "need more than 'neighbors' ({}) point features, found {n}",
                prm.neighbors
            )));
        }

        ctx.progress.info(&format!(
            "{n} point(s), {}-nearest, {} distance",
            prm.neighbors,
            prm.method.label()
        ));

        // kd-tree over the representative points.
        let mut tree: KdTree<f64, usize, [f64; 2]> = KdTree::new(2);
        for (i, p) in pts.iter().enumerate() {
            tree.add(*p, i)
                .map_err(|e| ToolError::Execution(format!("kd-tree insert failed: {e:?}")))?;
        }

        // Distance to the k-th nearest neighbour (excluding self) per feature.
        // Query for neighbors+1 then drop the feature's own index, so a *different*
        // feature at distance 0 (coincident point) still counts as a neighbour.
        let mut kth: Vec<f64> = Vec::with_capacity(n);
        for (i, p) in pts.iter().enumerate() {
            let found = match prm.method {
                Distance::Euclidean => {
                    tree.nearest(p, prm.neighbors + 1, &kdtree::distance::squared_euclidean)
                }
                Distance::Manhattan => tree.nearest(p, prm.neighbors + 1, &manhattan),
            }
            .map_err(|e| ToolError::Execution(format!("kd-tree query failed: {e:?}")))?;
            let kth_d = found
                .into_iter()
                .filter(|(_, &j)| j != i)
                .take(prm.neighbors)
                .last()
                .map(|(d, _)| match prm.method {
                    Distance::Euclidean => d.sqrt(),
                    Distance::Manhattan => d,
                })
                .ok_or_else(|| {
                    ToolError::Execution(format!(
                        "feature {i} has no {}-th neighbour",
                        prm.neighbors
                    ))
                })?;
            kth.push(kth_d);
        }

        let (mut mn, mut mx, mut sum) = (f64::INFINITY, f64::NEG_INFINITY, 0.0f64);
        for &d in &kth {
            mn = mn.min(d);
            mx = mx.max(d);
            sum += d;
        }
        let avg = sum / n as f64;

        ctx.progress.info(&format!(
            "distance band: min {mn:.4}, avg {avg:.4}, max {mx:.4}"
        ));

        let mut out_path: Option<String> = None;
        if let Some(path) = output {
            let mut csv = String::from("feature_id,kth_neighbor_distance\n");
            for (i, &d) in kth.iter().enumerate() {
                csv.push_str(&format!("{},{d}\n", feat_of[i]));
            }
            write_text_output(&csv, path)?;
            out_path = Some(path.to_string());
        }

        let mut outputs = BTreeMap::new();
        outputs.insert("neighbors".to_string(), json!(prm.neighbors));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("min_distance".to_string(), json!(mn));
        outputs.insert("avg_distance".to_string(), json!(avg));
        outputs.insert("max_distance".to_string(), json!(mx));
        outputs.insert("distance_method".to_string(), json!(prm.method.label()));
        outputs.insert("skipped_non_point".to_string(), json!(skipped_non_point));
        if let Some(p) = out_path {
            outputs.insert("output".to_string(), json!(p));
        }
        Ok(ToolRunResult { outputs })
    }
}

/// L1 (Manhattan) distance for the kd-tree query.
fn manhattan(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).sum()
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

// ── Parameters ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Distance {
    Euclidean,
    Manhattan,
}

impl Distance {
    fn label(&self) -> &'static str {
        match self {
            Distance::Euclidean => "euclidean",
            Distance::Manhattan => "manhattan",
        }
    }
}

struct Params {
    neighbors: usize,
    method: Distance,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    // Accept only integers >= 1: reject non-integers (3.5), negatives, and 0
    // rather than silently coercing them to 1.
    let invalid_neighbors =
        || ToolError::Validation("'neighbors' must be an integer >= 1".to_string());
    let neighbors = match args.get("neighbors") {
        None | Some(Value::Null) => 1,
        Some(Value::Number(n)) => {
            // Integer JSON numbers, or whole-valued floats like 3.0.
            let k = if let Some(u) = n.as_u64() {
                u
            } else {
                let f = n.as_f64().ok_or_else(invalid_neighbors)?;
                if f.fract() != 0.0 || f < 1.0 {
                    return Err(invalid_neighbors());
                }
                f as u64
            };
            if k < 1 {
                return Err(invalid_neighbors());
            }
            k as usize
        }
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => {
            let k = s.trim().parse::<usize>().map_err(|_| invalid_neighbors())?;
            if k < 1 {
                return Err(invalid_neighbors());
            }
            k
        }
        Some(_) => return Err(invalid_neighbors()),
    };
    let method = match args
        .get("distance_method")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some("euclidean") => Distance::Euclidean,
        Some("manhattan") => Distance::Manhattan,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'distance_method' must be 'euclidean' or 'manhattan', got '{o}'"
            )))
        }
    };
    Ok(Params { neighbors, method })
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

    fn layer_of(pts: &[(f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        for &(x, y) in pts {
            l.add_feature(Some(Geometry::point(x, y)), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        CalculateDistanceBandTool.run(&args, &ctx()).unwrap()
    }

    /// A regular unit-spaced 1-D chain: 1st-neighbour distance is 1 for every
    /// point, so min = avg = max = 1.
    #[test]
    fn unit_chain_first_neighbor() {
        let pts: Vec<(f64, f64)> = (0..6).map(|i| (i as f64, 0.0)).collect();
        let out = run(json!({ "input": layer_of(&pts), "neighbors": 1 }));
        assert!((out.outputs["min_distance"].as_f64().unwrap() - 1.0).abs() < 1e-9);
        assert!((out.outputs["max_distance"].as_f64().unwrap() - 1.0).abs() < 1e-9);
        assert!((out.outputs["avg_distance"].as_f64().unwrap() - 1.0).abs() < 1e-9);
    }

    /// 2nd-nearest neighbour on the same chain: interior points have their 2nd
    /// neighbour at distance 2, endpoints at distance 2 as well (0->2). So max=2.
    #[test]
    fn second_neighbor_band() {
        let pts: Vec<(f64, f64)> = (0..6).map(|i| (i as f64, 0.0)).collect();
        let out = run(json!({ "input": layer_of(&pts), "neighbors": 2 }));
        // Endpoint 0: neighbours at 1 and 2 -> 2nd = 2. Interior: 1 and 1(other side)? no,
        // neighbours of point 2 are {1,3} at dist 1, then {0,4} at dist 2 -> 2nd nearest = 1.
        // So min = 1 (interior), max = 2 (endpoints).
        assert!((out.outputs["min_distance"].as_f64().unwrap() - 1.0).abs() < 1e-9);
        assert!((out.outputs["max_distance"].as_f64().unwrap() - 2.0).abs() < 1e-9);
    }

    /// Manhattan distance on a diagonal pair: L1 of (0,0)->(3,4) is 7.
    #[test]
    fn manhattan_metric() {
        let pts = vec![(0.0, 0.0), (3.0, 4.0)];
        let out =
            run(json!({ "input": layer_of(&pts), "neighbors": 1, "distance_method": "manhattan" }));
        assert!((out.outputs["max_distance"].as_f64().unwrap() - 7.0).abs() < 1e-9);
        // Euclidean would be 5.
        let oute = run(json!({ "input": layer_of(&pts), "neighbors": 1 }));
        assert!((oute.outputs["max_distance"].as_f64().unwrap() - 5.0).abs() < 1e-9);
    }

    /// Duplicate (coincident) points: a second point at distance 0 is a real
    /// neighbour, so the 1st-neighbour distance for a duplicated location is 0.
    #[test]
    fn coincident_points_zero_distance() {
        let pts = vec![(0.0, 0.0), (0.0, 0.0), (10.0, 0.0)];
        let out = run(json!({ "input": layer_of(&pts), "neighbors": 1 }));
        assert!((out.outputs["min_distance"].as_f64().unwrap()).abs() < 1e-9);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CalculateDistanceBandTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "distance_method": "chebyshev" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "neighbors": 3 })).is_ok());
        // Strict neighbors: reject non-integers, 0, and negatives; accept 3.0.
        assert!(bad(json!({ "input": "a.geojson", "neighbors": 3.5 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "neighbors": 0 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "neighbors": -2 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "neighbors": 3.0 })).is_ok());
    }
}
