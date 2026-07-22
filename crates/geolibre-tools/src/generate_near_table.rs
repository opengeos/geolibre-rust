//! GeoLibre tool: k-nearest / within-radius proximity table.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Generate Near Table* (Analysis).
//! The bundled whitebox `near` writes a single `NEAR_FID`/`NEAR_DIST` per input
//! feature — one nearest neighbour, distance only, no k-nearest, no near-angle,
//! no near-XY. This tool emits a standalone proximity table: for every input
//! feature, the `closest_count` nearest near-layer features (or every near
//! feature within `search_radius`), each row carrying
//!
//! * `IN_FID`    — input feature index (0-based),
//! * `NEAR_FID`  — near feature index (0-based),
//! * `NEAR_DIST` — planar distance between representative points,
//! * `NEAR_RANK` — 1 = nearest, 2 = second nearest, …,
//! * `NEAR_ANGLE`— (when `angle=true`) arithmetic bearing in degrees, east = 0,
//!   counter-clockwise, in `(-180, 180]`,
//! * `NEAR_X` / `NEAR_Y` — (when `location=true`) the near feature's location.
//!
//! Each output feature is a `Point` at the matched near location, so the table
//! renders directly with `render_vector_png`. Distances are planar (straight
//! line) between representative points (centroid for lines/polygons); geodesic
//! distance and true edge-to-edge distance are out of scope for v1.
//!
//! Neighbour search uses the vendored `kdtree` (same crate as
//! `neighborhood_summary_statistics`). When `input` and `near_features` resolve
//! to the same path, a feature never matches itself.

use std::collections::BTreeMap;

use geo::{Centroid, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct GenerateNearTableTool;

impl Tool for GenerateNearTableTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "generate_near_table",
            display_name: "Generate Near Table",
            summary: "For each input feature, list the k nearest (or all within a search radius) near-layer features with distance, rank, optional bearing angle and near-XY — a multi-neighbour proximity table, like ArcGIS Generate Near Table.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (points, lines, or polygons) whose proximity is measured.",
                    required: true,
                },
                ToolParamSpec {
                    name: "near_features",
                    description: "Near vector layer to search against (points, lines, or polygons).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output point vector path (one Point per near match). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "search_radius",
                    description: "Only keep near features within this distance (map units). Required when closest_count = 0 (all within radius).",
                    required: false,
                },
                ToolParamSpec {
                    name: "closest_count",
                    description: "Number of nearest near features per input feature (k, default 1). Use 0 to return every near feature within search_radius.",
                    required: false,
                },
                ToolParamSpec {
                    name: "angle",
                    description: "Add NEAR_ANGLE (arithmetic bearing in degrees, east = 0, CCW). Default false.",
                    required: false,
                },
                ToolParamSpec {
                    name: "location",
                    description: "Add NEAR_X / NEAR_Y columns for the near location. Default false.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["input", "near_features"] {
            if args
                .get(key)
                .and_then(Value::as_str)
                .map(str::trim)
                .unwrap_or("")
                .is_empty()
            {
                return Err(ToolError::Validation(format!(
                    "missing required string parameter '{key}'"
                )));
            }
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input_path = args.get("input").and_then(Value::as_str).unwrap();
        let near_path = args.get("near_features").and_then(Value::as_str).unwrap();
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let input = load_input_layer(input_path)?;
        let near = load_input_layer(near_path)?;

        // Representative point per feature (None -> no geometry, skipped).
        let in_reps: Vec<Option<(f64, f64)>> = input
            .features
            .iter()
            .map(|f| f.geometry.as_ref().and_then(representative_point))
            .collect();
        let near_reps: Vec<Option<(f64, f64)>> = near
            .features
            .iter()
            .map(|f| f.geometry.as_ref().and_then(representative_point))
            .collect();

        if near_reps.iter().all(Option::is_none) {
            return Err(ToolError::Execution(
                "near_features has no usable geometry".to_string(),
            ));
        }

        // Same-layer self-exclusion when both paths resolve to the same source.
        let self_join = input_path == near_path;

        // Build the kd-tree over near representative points.
        let mut tree: KdTree<f64, usize, [f64; 2]> = KdTree::new(2);
        for (j, r) in near_reps.iter().enumerate() {
            if let Some((x, y)) = r {
                tree.add([*x, *y], j)
                    .map_err(|e| ToolError::Execution(format!("kd-tree insert failed: {e:?}")))?;
            }
        }

        ctx.progress.info(&format!(
            "matching {} input features against {} near features",
            in_reps.len(),
            near_reps.len()
        ));

        // ── Output schema: table columns on a Point layer ────────────────────────
        let mut out = Layer::new("near_table").with_geom_type(GeometryType::Point);
        if let Some(epsg) = input.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("IN_FID", FieldType::Integer));
        out.add_field(FieldDef::new("NEAR_FID", FieldType::Integer));
        out.add_field(FieldDef::new("NEAR_DIST", FieldType::Float));
        out.add_field(FieldDef::new("NEAR_RANK", FieldType::Integer));
        if prm.angle {
            out.add_field(FieldDef::new("NEAR_ANGLE", FieldType::Float));
        }
        if prm.location {
            out.add_field(FieldDef::new("NEAR_X", FieldType::Float));
            out.add_field(FieldDef::new("NEAR_Y", FieldType::Float));
        }

        let radius2 = prm.search_radius.map(|r| r * r);
        let mut row_count = 0usize;
        let mut matched_inputs = 0usize;
        let mut total_dist = 0.0;

        for (i, rep) in in_reps.iter().enumerate() {
            let Some((x, y)) = rep.map(|(x, y)| (x, y)) else {
                continue;
            };

            // Candidate near features (index, squared distance), sorted ascending.
            let cands: Vec<(f64, usize)> = if prm.closest_count == 0 {
                // "all within radius" mode (radius guaranteed present by parse).
                let r2 = radius2.unwrap();
                tree.within(&[x, y], r2, &squared_euclidean)
                    .map_err(|e| ToolError::Execution(format!("kd-tree query failed: {e:?}")))?
                    .into_iter()
                    .map(|(d2, &j)| (d2, j))
                    .filter(|&(_, j)| !(self_join && j == i))
                    .collect()
            } else {
                // k nearest, then radius filter. Query a couple extra to survive
                // dropping the self-match in a self-join.
                let want = prm.closest_count + if self_join { 1 } else { 0 };
                tree.nearest(&[x, y], want, &squared_euclidean)
                    .map_err(|e| ToolError::Execution(format!("kd-tree query failed: {e:?}")))?
                    .into_iter()
                    .map(|(d2, &j)| (d2, j))
                    .filter(|&(_, j)| !(self_join && j == i))
                    .filter(|&(d2, _)| radius2.is_none_or(|r2| d2 <= r2))
                    .take(prm.closest_count)
                    .collect()
            };

            if !cands.is_empty() {
                matched_inputs += 1;
            }

            for (rank, (d2, j)) in cands.into_iter().enumerate() {
                let dist = d2.sqrt();
                let (nx, ny) = near_reps[j].unwrap();
                let mut attrs: Vec<(&str, FieldValue)> = vec![
                    ("IN_FID", FieldValue::Integer(i as i64)),
                    ("NEAR_FID", FieldValue::Integer(j as i64)),
                    ("NEAR_DIST", FieldValue::Float(dist)),
                    ("NEAR_RANK", FieldValue::Integer(rank as i64 + 1)),
                ];
                if prm.angle {
                    let ang = (ny - y).atan2(nx - x).to_degrees();
                    attrs.push(("NEAR_ANGLE", FieldValue::Float(ang)));
                }
                if prm.location {
                    attrs.push(("NEAR_X", FieldValue::Float(nx)));
                    attrs.push(("NEAR_Y", FieldValue::Float(ny)));
                }
                out.add_feature(Some(Geometry::Point(Coord::xy(nx, ny))), &attrs)
                    .map_err(|e| ToolError::Execution(format!("failed writing row: {e}")))?;
                row_count += 1;
                total_dist += dist;
            }
        }

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_count".to_string(), json!(input.features.len()));
        outputs.insert("near_count".to_string(), json!(near.features.len()));
        outputs.insert("row_count".to_string(), json!(row_count));
        outputs.insert("matched_input_count".to_string(), json!(matched_inputs));
        outputs.insert(
            "mean_distance".to_string(),
            json!(if row_count > 0 {
                total_dist / row_count as f64
            } else {
                0.0
            }),
        );
        Ok(ToolRunResult { outputs })
    }
}

// ── Representative points ─────────────────────────────────────────────────────

fn representative_point(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => {
            let (sx, sy) = cs
                .iter()
                .fold((0.0, 0.0), |(ax, ay), c| (ax + c.x, ay + c.y));
            let k = cs.len() as f64;
            Some((sx / k, sy / k))
        }
        Geometry::LineString(cs) if !cs.is_empty() => {
            let ls = LineString::new(cs.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect());
            ls.centroid().map(|p| (p.x(), p.y()))
        }
        Geometry::MultiLineString(parts) => {
            let mls = geo::MultiLineString(
                parts
                    .iter()
                    .map(|cs| {
                        LineString::new(cs.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect())
                    })
                    .collect(),
            );
            mls.centroid().map(|p| (p.x(), p.y()))
        }
        Geometry::Polygon { .. } | Geometry::MultiPolygon(_) => to_multipolygon(geom)
            .and_then(|mp| mp.centroid())
            .map(|p| (p.x(), p.y())),
        _ => None,
    }
}

fn to_multipolygon(geom: &Geometry) -> Option<MultiPolygon> {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => Some(MultiPolygon(vec![rings_to_polygon(exterior, interiors)])),
        Geometry::MultiPolygon(parts) => Some(MultiPolygon(
            parts.iter().map(|(e, i)| rings_to_polygon(e, i)).collect(),
        )),
        _ => None,
    }
}

fn rings_to_polygon(exterior: &Ring, interiors: &[Ring]) -> Polygon {
    Polygon::new(
        ring_to_linestring(exterior),
        interiors.iter().map(ring_to_linestring).collect(),
    )
}

fn ring_to_linestring(ring: &Ring) -> LineString {
    LineString::new(
        ring.coords()
            .iter()
            .map(|c| GeoCoord { x: c.x, y: c.y })
            .collect(),
    )
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    search_radius: Option<f64>,
    closest_count: usize,
    angle: bool,
    location: bool,
}

fn parse_optional_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => {
            Ok(Some(s.trim().parse::<f64>().map_err(|_| {
                ToolError::Validation(format!("'{key}' must be a number"))
            })?))
        }
        Some(_) => Err(ToolError::Validation(format!("'{key}' must be a number"))),
    }
}

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<bool, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(false),
        Some(Value::Bool(b)) => Ok(*b),
        Some(Value::String(s)) => match s.trim().to_lowercase().as_str() {
            "" => Ok(false),
            "true" | "1" | "yes" => Ok(true),
            "false" | "0" | "no" => Ok(false),
            _ => Err(ToolError::Validation(format!("'{key}' must be a boolean"))),
        },
        Some(Value::Number(n)) => Ok(n.as_f64().map(|v| v != 0.0).unwrap_or(false)),
        Some(_) => Err(ToolError::Validation(format!("'{key}' must be a boolean"))),
    }
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let search_radius = parse_optional_f64(args, "search_radius")?;
    if let Some(r) = search_radius {
        if r.is_nan() || r <= 0.0 {
            return Err(ToolError::Validation(
                "'search_radius' must be positive".to_string(),
            ));
        }
    }
    let closest_count = match args.get("closest_count") {
        None | Some(Value::Null) => 1,
        Some(Value::Number(n)) => n.as_u64().ok_or_else(|| {
            ToolError::Validation("'closest_count' must be a non-negative integer".into())
        })? as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'closest_count' must be an integer".into()))?,
        Some(_) => {
            return Err(ToolError::Validation(
                "'closest_count' must be a number".to_string(),
            ))
        }
    };
    if closest_count == 0 && search_radius.is_none() {
        return Err(ToolError::Validation(
            "'search_radius' is required when 'closest_count' = 0 (all within radius)".to_string(),
        ));
    }
    Ok(Params {
        search_radius,
        closest_count,
        angle: parse_optional_bool(args, "angle")?,
        location: parse_optional_bool(args, "location")?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::memory_store;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn point_layer(pts: &[(f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("id", FieldType::Integer));
        for (i, (x, y)) in pts.iter().enumerate() {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(*x, *y))),
                &[("id", (i as i64).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = GenerateNearTableTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn col(layer: &Layer, name: &str) -> usize {
        layer.schema.field_index(name).unwrap()
    }

    /// Default k=1: each input gets its single nearest near feature, with a
    /// correct 3-4-5 distance.
    #[test]
    fn nearest_only_default() {
        let input = point_layer(&[(0.0, 0.0)]);
        let near = point_layer(&[(3.0, 4.0), (100.0, 0.0)]);
        let (out, layer) = run(json!({ "input": input, "near_features": near }));
        assert_eq!(out.outputs["row_count"], json!(1));
        let f = layer.iter().next().unwrap();
        assert_eq!(f.attributes[col(&layer, "NEAR_FID")].as_i64().unwrap(), 0);
        assert!((f.attributes[col(&layer, "NEAR_DIST")].as_f64().unwrap() - 5.0).abs() < 1e-9);
        assert_eq!(f.attributes[col(&layer, "NEAR_RANK")].as_i64().unwrap(), 1);
    }

    /// closest_count = k returns the k nearest ordered by ascending rank.
    #[test]
    fn k_nearest_ranked() {
        let input = point_layer(&[(0.0, 0.0)]);
        let near = point_layer(&[(1.0, 0.0), (2.0, 0.0), (50.0, 0.0)]);
        let (out, layer) =
            run(json!({ "input": input, "near_features": near, "closest_count": 2 }));
        assert_eq!(out.outputs["row_count"], json!(2));
        let ranks: Vec<i64> = layer
            .iter()
            .map(|f| f.attributes[col(&layer, "NEAR_RANK")].as_i64().unwrap())
            .collect();
        assert_eq!(ranks, vec![1, 2]);
        // Rank 1 is the closest (x=1), rank 2 is next (x=2).
        let fids: Vec<i64> = layer
            .iter()
            .map(|f| f.attributes[col(&layer, "NEAR_FID")].as_i64().unwrap())
            .collect();
        assert_eq!(fids, vec![0, 1]);
    }

    /// search_radius with closest_count = 0 keeps every near feature within the
    /// band and drops the rest.
    #[test]
    fn all_within_radius() {
        let input = point_layer(&[(0.0, 0.0)]);
        let near = point_layer(&[(3.0, 0.0), (8.0, 0.0), (50.0, 0.0)]);
        let (out, _l) = run(json!({
            "input": input, "near_features": near, "closest_count": 0, "search_radius": 10.0
        }));
        // 3 and 8 within 10; 50 excluded.
        assert_eq!(out.outputs["row_count"], json!(2));
    }

    /// search_radius also caps the k-nearest mode; an input with no near feature
    /// in range produces no rows.
    #[test]
    fn radius_filters_k_nearest_and_unmatched_pass_through() {
        let input = point_layer(&[(0.0, 0.0), (1000.0, 1000.0)]);
        let near = point_layer(&[(2.0, 0.0)]);
        let (out, layer) = run(json!({
            "input": input, "near_features": near, "closest_count": 5, "search_radius": 10.0
        }));
        // Only the first input has a near feature within 10 units.
        assert_eq!(out.outputs["row_count"], json!(1));
        assert_eq!(out.outputs["matched_input_count"], json!(1));
        let f = layer.iter().next().unwrap();
        assert_eq!(f.attributes[col(&layer, "IN_FID")].as_i64().unwrap(), 0);
    }

    /// angle=true adds NEAR_ANGLE with the arithmetic bearing (east = 0, CCW).
    #[test]
    fn near_angle_bearing() {
        let input = point_layer(&[(0.0, 0.0)]);
        let near = point_layer(&[(0.0, 5.0)]); // due north -> 90 degrees.
        let (_o, layer) = run(json!({
            "input": input, "near_features": near, "angle": true, "location": true
        }));
        let f = layer.iter().next().unwrap();
        assert!((f.attributes[col(&layer, "NEAR_ANGLE")].as_f64().unwrap() - 90.0).abs() < 1e-9);
        assert!((f.attributes[col(&layer, "NEAR_X")].as_f64().unwrap() - 0.0).abs() < 1e-9);
        assert!((f.attributes[col(&layer, "NEAR_Y")].as_f64().unwrap() - 5.0).abs() < 1e-9);
    }

    /// A self-join (input path == near path) never matches a feature to itself.
    #[test]
    fn self_join_excludes_self() {
        let layer_path = point_layer(&[(0.0, 0.0), (10.0, 0.0), (20.0, 0.0)]);
        let (out, layer) = run(json!({
            "input": layer_path.clone(), "near_features": layer_path, "closest_count": 1
        }));
        assert_eq!(out.outputs["row_count"], json!(3));
        for f in layer.iter() {
            let inf = f.attributes[col(&layer, "IN_FID")].as_i64().unwrap();
            let nf = f.attributes[col(&layer, "NEAR_FID")].as_i64().unwrap();
            assert_ne!(inf, nf, "no feature matches itself");
        }
    }

    #[test]
    fn rejects_missing_near_features() {
        let input = point_layer(&[(0.0, 0.0)]);
        let args: ToolArgs = serde_json::from_value(json!({ "input": input })).unwrap();
        assert!(GenerateNearTableTool.validate(&args).is_err());
    }

    #[test]
    fn rejects_all_within_radius_without_radius() {
        let input = point_layer(&[(0.0, 0.0)]);
        let near = point_layer(&[(1.0, 0.0)]);
        let args: ToolArgs = serde_json::from_value(
            json!({ "input": input, "near_features": near, "closest_count": 0 }),
        )
        .unwrap();
        assert!(GenerateNearTableTool.validate(&args).is_err());
    }
}
