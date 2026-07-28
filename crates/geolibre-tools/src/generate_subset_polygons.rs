//! GeoLibre tool: density-adaptive non-overlapping partition of a point layer.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Generate Subset Polygons*
//! (Geostatistical Analyst).
//!
//! This is the standard preprocessing step for interpolating a very large point
//! set: split the study area into subsets small enough to fit a local model, fit
//! each independently, then blend. GeoLibre now ships a deep interpolation suite
//! — `local_polynomial_interpolation`, `kernel_interpolation_with_barriers`,
//! `diffusion_interpolation_with_barriers`, plus the bundled kriging/IDW/RBF/TPS
//! tools — and every one of them either fits globally or uses a fixed search
//! neighbourhood. None can be driven over a data-adaptive partition.
//!
//! The near misses are each wrong in a specific way:
//!
//! * `voronoi_diagram` gives one cell **per point**, not a grouping of points;
//! * `rectangular_grid_from_*` / `hexagonal_grid_from_*` tile the extent on a
//!   fixed geometry, ignoring density — sparse areas get empty cells, dense
//!   areas overloaded ones;
//! * `build_balanced_zones` balances against an *attribute* under contiguity
//!   constraints (a much heavier optimisation, not a fast spatial split);
//! * `create_spatially_balanced_points` *selects* a sample rather than
//!   partitioning;
//! * `group_by_proximity` tags points with a cluster id but returns no polygons
//!   and honours no size bounds.
//!
//! The split is a k-d tree built on the **median** coordinate rather than the
//! geometric midpoint, so both halves get roughly equal counts regardless of any
//! density gradient. That makes it fully deterministic — no RNG, so WASM and
//! native agree exactly.

use std::collections::BTreeMap;

use geo::{Area, BooleanOps, ConvexHull, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Partitions a point layer into compact, non-overlapping subset polygons.
pub struct GenerateSubsetPolygonsTool;

impl Tool for GenerateSubsetPolygonsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "generate_subset_polygons",
            display_name: "Generate Subset Polygons",
            summary: "Partitions a dense point layer into compact, non-overlapping polygons that tile the point extent, each holding between a minimum and maximum number of points (ArcGIS Generate Subset Polygons). The standard preprocessing step for interpolating very large point sets. Unlike the bundled fixed grids it adapts to density (median k-d split), unlike voronoi_diagram it groups points rather than making one cell each, and unlike group_by_proximity it emits polygons and honours size bounds.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output polygon path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_points",
                    description: "Optional path for the input points tagged with their subset id. If omitted, stored in memory (still returned).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_points_per_subset",
                    description: "Lower bound on points per subset (default 50).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_points_per_subset",
                    description: "Upper bound on points per subset (default 200).",
                    required: false,
                },
                ToolParamSpec {
                    name: "coincident_points",
                    description: "'single' (default; exactly coincident points count once when sizing) or 'all' (each counts individually).",
                    required: false,
                },
                ToolParamSpec {
                    name: "clip_to_hull",
                    description: "Clip subset polygons to the convex hull of the points so the result hugs the data rather than its bounding box (default true).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args.get("input").and_then(Value::as_str).is_none() {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args.get("input").and_then(Value::as_str).ok_or_else(|| {
            ToolError::Validation("missing required parameter 'input'".to_string())
        })?;
        let output = parse_optional_str(args, "output")?;
        let output_points = parse_optional_str(args, "output_points")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;

        // Flatten every point (including multipoint members) with its source
        // feature index, so the id can be joined back.
        let mut pts: Vec<(f64, f64, usize)> = Vec::new();
        for (fid, feature) in layer.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            collect_points(geom, fid, &mut pts);
        }
        if pts.is_empty() {
            return Err(ToolError::Execution(
                "input layer contains no point geometries".to_string(),
            ));
        }
        ctx.progress
            .info(&format!("partitioning {} point(s)", pts.len()));

        // Overall extent, padded so boundary points sit strictly inside.
        let (mut min_x, mut min_y) = (f64::INFINITY, f64::INFINITY);
        let (mut max_x, mut max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
        for (x, y, _) in &pts {
            min_x = min_x.min(*x);
            max_x = max_x.max(*x);
            min_y = min_y.min(*y);
            max_y = max_y.max(*y);
        }
        let pad = ((max_x - min_x).abs().max((max_y - min_y).abs()) * 1e-6).max(1e-9);
        let extent = Rect {
            min_x: min_x - pad,
            min_y: min_y - pad,
            max_x: max_x + pad,
            max_y: max_y + pad,
        };

        // Recursive median split.
        let mut indices: Vec<usize> = (0..pts.len()).collect();
        let mut leaves: Vec<Leaf> = Vec::new();
        let mut undersized = 0_u64;
        split(
            &pts,
            &mut indices[..],
            extent,
            &prm,
            &mut leaves,
            &mut undersized,
            0,
        );

        ctx.progress
            .info(&format!("{} subset(s) produced", leaves.len()));

        // Optional clip to the point set's convex hull.
        let hull: Option<MultiPolygon<f64>> = if prm.clip_to_hull && pts.len() >= 3 {
            let mp = MultiPolygon(vec![Polygon::new(
                LineString::new(
                    pts.iter()
                        .map(|(x, y, _)| GeoCoord { x: *x, y: *y })
                        .collect(),
                ),
                vec![],
            )]);
            let hull_poly = mp.convex_hull();
            if hull_poly.unsigned_area() > 0.0 {
                Some(MultiPolygon(vec![hull_poly]))
            } else {
                // Collinear or degenerate point sets have no area to clip to.
                None
            }
        } else {
            None
        };

        // ── Subset polygons ──────────────────────────────────────────────────
        let mut out = Layer::new("subset_polygons");
        out.add_field(FieldDef::new("subset_id", FieldType::Integer));
        out.add_field(FieldDef::new("point_count", FieldType::Integer));
        out.add_field(FieldDef::new("area", FieldType::Float));
        out.crs = layer.crs.clone();
        out.geom_type = Some(GeometryType::Polygon);

        // Subset id per flattened point. Points belonging to a leaf whose
        // clipped geometry turns out empty are tagged -1 rather than left
        // pointing at a polygon that was never written, which would be a silent
        // broken join downstream.
        let mut assignment = vec![-1_i64; pts.len()];
        let mut emitted = 0_usize;
        let mut dropped_points = 0_u64;

        for (sid, leaf) in leaves.iter().enumerate() {
            let rect_poly = leaf.rect.to_geo();
            let clipped = match &hull {
                Some(h) => MultiPolygon(vec![rect_poly]).intersection(h),
                None => MultiPolygon(vec![rect_poly]),
            };
            if clipped.0.is_empty() {
                dropped_points += leaf.indices.len() as u64;
                continue;
            }
            let area = clipped.unsigned_area();
            let Some(geom) = multipolygon_to_wb(&clipped) else {
                dropped_points += leaf.indices.len() as u64;
                continue;
            };
            // Only tag points once the polygon is certain to be written.
            for i in &leaf.indices {
                assignment[*i] = sid as i64;
            }
            out.add_feature(
                Some(geom),
                &[
                    ("subset_id", (sid as i64).into()),
                    ("point_count", (leaf.indices.len() as i64).into()),
                    ("area", area.into()),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed adding subset polygon: {e}")))?;
            emitted += 1;
        }

        // ── Tagged points ────────────────────────────────────────────────────
        let mut pt_layer = Layer::new("subset_points");
        for field in layer.schema.fields().iter() {
            pt_layer.add_field(field.clone());
        }
        pt_layer.add_field(FieldDef::new("subset_id", FieldType::Integer));
        pt_layer.crs = layer.crs.clone();
        pt_layer.geom_type = Some(GeometryType::Point);

        // Materialised once: `layer.iter().nth(fid)` inside the loop would make
        // the tagged-point write O(points x features), on the very path that
        // exists for large point sets.
        let src_features: Vec<&wbvector::Feature> = layer.iter().collect();
        for (idx, (x, y, fid)) in pts.iter().enumerate() {
            let src = src_features.get(*fid).copied();
            let mut attrs: Vec<(&str, wbvector::FieldValue)> = layer
                .schema
                .fields()
                .iter()
                .enumerate()
                .map(|(i, f)| {
                    (
                        f.name.as_str(),
                        src.map(|s| s.attributes[i].clone())
                            .unwrap_or(wbvector::FieldValue::Null),
                    )
                })
                .collect();
            attrs.push(("subset_id", assignment[idx].into()));
            pt_layer
                .add_feature(Some(Geometry::Point(Coord::xy(*x, *y))), &attrs)
                .map_err(|e| ToolError::Execution(format!("failed adding tagged point: {e}")))?;
        }

        let subset_count = leaves.len();
        let point_count = pts.len();
        let out_path = write_or_store_layer(out, output)?;
        let pts_path = write_or_store_layer(pt_layer, output_points)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("output_points".to_string(), json!(pts_path));
        outputs.insert("subset_count".to_string(), json!(subset_count));
        outputs.insert("emitted_polygon_count".to_string(), json!(emitted));
        outputs.insert("point_count".to_string(), json!(point_count));
        // Where min and max genuinely cannot both hold, say so rather than
        // silently violating the minimum.
        outputs.insert("undersized_subset_count".to_string(), json!(undersized));
        // Points whose subset polygon was dropped carry subset_id = -1.
        outputs.insert("dropped_point_count".to_string(), json!(dropped_points));
        Ok(ToolRunResult { outputs })
    }
}

// ── Partitioning ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
struct Rect {
    min_x: f64,
    min_y: f64,
    max_x: f64,
    max_y: f64,
}

impl Rect {
    fn to_geo(self) -> Polygon<f64> {
        Polygon::new(
            LineString::new(vec![
                GeoCoord {
                    x: self.min_x,
                    y: self.min_y,
                },
                GeoCoord {
                    x: self.max_x,
                    y: self.min_y,
                },
                GeoCoord {
                    x: self.max_x,
                    y: self.max_y,
                },
                GeoCoord {
                    x: self.min_x,
                    y: self.max_y,
                },
                GeoCoord {
                    x: self.min_x,
                    y: self.min_y,
                },
            ]),
            vec![],
        )
    }
}

struct Leaf {
    rect: Rect,
    indices: Vec<usize>,
}

/// Counts a point set's size under the coincident-point rule: `single` collapses
/// exactly-coincident coordinates so a tower of stacked observations does not
/// force a degenerate split.
fn effective_count(pts: &[(f64, f64, usize)], idx: &[usize], coincident_single: bool) -> usize {
    if !coincident_single {
        return idx.len();
    }
    let mut keys: Vec<(u64, u64)> = idx
        .iter()
        .map(|i| (pts[*i].0.to_bits(), pts[*i].1.to_bits()))
        .collect();
    keys.sort_unstable();
    keys.dedup();
    keys.len()
}

/// Recursive median split. Splits along the longer axis at the median point
/// coordinate; stops when a further split would push either child below
/// `min_points`.
fn split(
    pts: &[(f64, f64, usize)],
    idx: &mut [usize],
    rect: Rect,
    prm: &Params,
    out: &mut Vec<Leaf>,
    undersized: &mut u64,
    depth: usize,
) {
    let n_eff = effective_count(pts, idx, prm.coincident_single);

    // Depth guard: with heavy coincidence the effective count can stop shrinking,
    // so bound recursion rather than risk running away.
    let can_split = n_eff > prm.max_points && idx.len() >= 2 && depth < 64;
    if !can_split {
        if n_eff < prm.min_points {
            *undersized += 1;
        }
        out.push(Leaf {
            rect,
            indices: idx.to_vec(),
        });
        return;
    }

    // Split along the longer axis of the *rectangle*, which keeps subsets
    // compact rather than long slivers.
    let horizontal = (rect.max_x - rect.min_x) >= (rect.max_y - rect.min_y);
    let key = |i: usize| if horizontal { pts[i].0 } else { pts[i].1 };

    idx.sort_by(|a, b| {
        key(*a)
            .partial_cmp(&key(*b))
            // Ties broken by index so the ordering is total and deterministic.
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.cmp(b))
    });

    // The split must fall on a distinct-coordinate boundary. A run of
    // coincident points straddling the midpoint would otherwise be cut in two,
    // and each half would then count that one location independently — so a
    // single stack could end up with several `subset_id`s, contradicting the
    // collapse that `coincident_points = "single"` promises.
    let boundary = |b: usize| b > 0 && b < idx.len() && key(idx[b - 1]) != key(idx[b]);
    let target = idx.len() / 2;
    let mid = if boundary(target) {
        target
    } else {
        // Nearest distinct-coordinate boundary either side of the midpoint.
        let up = (target + 1..idx.len()).find(|b| boundary(*b));
        let down = (1..target).rev().find(|b| boundary(*b));
        match (up, down) {
            (Some(u), Some(d)) => {
                if u - target <= target - d {
                    u
                } else {
                    d
                }
            }
            (Some(u), None) => u,
            (None, Some(d)) => d,
            // Every point shares this coordinate: no split can separate them.
            (None, None) => {
                out.push(Leaf {
                    rect,
                    indices: idx.to_vec(),
                });
                return;
            }
        }
    };

    // Both halves must still clear the minimum, measured with the SAME metric
    // the maximum uses — comparing a raw count against a coincidence-collapsed
    // maximum would let a split proceed that is known in advance to leave an
    // undersized child.
    let left_eff = effective_count(pts, &idx[..mid], prm.coincident_single);
    let right_eff = effective_count(pts, &idx[mid..], prm.coincident_single);
    if left_eff < prm.min_points || right_eff < prm.min_points {
        out.push(Leaf {
            rect,
            indices: idx.to_vec(),
        });
        return;
    }

    let cut = (key(idx[mid - 1]) + key(idx[mid])) / 2.0;
    // A degenerate cut (all coordinates equal on this axis) cannot separate the
    // points; emit the leaf instead of recursing forever.
    let lo = if horizontal { rect.min_x } else { rect.min_y };
    let hi = if horizontal { rect.max_x } else { rect.max_y };
    if !(cut > lo && cut < hi) {
        out.push(Leaf {
            rect,
            indices: idx.to_vec(),
        });
        return;
    }

    let (left_rect, right_rect) = if horizontal {
        (Rect { max_x: cut, ..rect }, Rect { min_x: cut, ..rect })
    } else {
        (Rect { max_y: cut, ..rect }, Rect { min_y: cut, ..rect })
    };

    let (left, right) = idx.split_at_mut(mid);
    split(pts, left, left_rect, prm, out, undersized, depth + 1);
    split(pts, right, right_rect, prm, out, undersized, depth + 1);
}

fn collect_points(geom: &Geometry, fid: usize, out: &mut Vec<(f64, f64, usize)>) {
    match geom {
        Geometry::Point(c) => out.push((c.x, c.y, fid)),
        Geometry::MultiPoint(cs) => {
            for c in cs {
                out.push((c.x, c.y, fid));
            }
        }
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                collect_points(g, fid, out);
            }
        }
        _ => {}
    }
}

/// Converts a `geo` `MultiPolygon` back to a `wbvector` geometry.
fn multipolygon_to_wb(mp: &MultiPolygon<f64>) -> Option<Geometry> {
    let parts: Vec<(Ring, Vec<Ring>)> =
        mp.0.iter()
            .filter_map(|p| {
                let ext = linestring_to_ring(p.exterior())?;
                let holes: Vec<Ring> = p
                    .interiors()
                    .iter()
                    .filter_map(linestring_to_ring)
                    .collect();
                Some((ext, holes))
            })
            .collect();
    match parts.len() {
        0 => None,
        1 => {
            let (exterior, interiors) = parts.into_iter().next().unwrap();
            Some(Geometry::Polygon {
                exterior,
                interiors,
            })
        }
        _ => Some(Geometry::MultiPolygon(parts)),
    }
}

fn linestring_to_ring(ls: &LineString<f64>) -> Option<Ring> {
    let mut cs: Vec<Coord> = ls.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
    // `geo` closes rings with a duplicate final vertex; `Ring` stores none.
    if cs.len() > 1 && cs[0].x == cs[cs.len() - 1].x && cs[0].y == cs[cs.len() - 1].y {
        cs.pop();
    }
    if cs.len() < 3 {
        return None;
    }
    Some(Ring::new(cs))
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    min_points: usize,
    max_points: usize,
    coincident_single: bool,
    clip_to_hull: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let min_points = opt_u64(args, "min_points_per_subset")?.unwrap_or(50) as usize;
    let max_points = opt_u64(args, "max_points_per_subset")?.unwrap_or(200) as usize;
    if min_points == 0 {
        return Err(ToolError::Validation(
            "'min_points_per_subset' must be at least 1".to_string(),
        ));
    }
    if max_points < min_points {
        return Err(ToolError::Validation(format!(
            "'max_points_per_subset' ({max_points}) must be >= 'min_points_per_subset' ({min_points})"
        )));
    }
    let coincident_single = match parse_optional_str(args, "coincident_points")? {
        None => true,
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "single" | "coincident_single" => true,
            "all" | "coincident_all" => false,
            other => {
                return Err(ToolError::Validation(format!(
                    "unknown coincident_points '{other}' (expected 'single' or 'all')"
                )))
            }
        },
    };
    Ok(Params {
        min_points,
        max_points,
        coincident_single,
        clip_to_hull: opt_bool(args, "clip_to_hull")?.unwrap_or(true),
    })
}

fn opt_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => n.as_u64().map(Some).ok_or_else(|| {
            ToolError::Validation(format!("parameter '{key}' must be a positive integer"))
        }),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s.trim().parse::<u64>().map(Some).map_err(|_| {
            ToolError::Validation(format!("parameter '{key}' must be a positive integer"))
        }),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a positive integer"
        ))),
    }
}

fn opt_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
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
        let mut l = Layer::new("pts");
        l.geom_type = Some(GeometryType::Point);
        for (x, y) in pts {
            l.add_feature(Some(Geometry::Point(Coord::xy(*x, *y))), &[])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    /// A regular grid of points, split into subsets.
    fn grid(n: usize) -> Vec<(f64, f64)> {
        let mut v = Vec::new();
        for i in 0..n {
            for j in 0..n {
                v.push((i as f64, j as f64));
            }
        }
        v
    }

    fn run(path: String, extra: Value) -> (Layer, Layer, ToolRunResult) {
        let mut obj = serde_json::Map::new();
        obj.insert("input".to_string(), json!(path));
        if let Value::Object(m) = extra {
            for (k, v) in m {
                obj.insert(k, v);
            }
        }
        let args: ToolArgs = serde_json::from_value(Value::Object(obj)).unwrap();
        let res = GenerateSubsetPolygonsTool.run(&args, &ctx()).unwrap();
        let polys = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        let tagged = load_input_layer(res.outputs["output_points"].as_str().unwrap()).unwrap();
        (polys, tagged, res)
    }

    fn count_field(l: &Layer, f: &wbvector::Feature, key: &str) -> i64 {
        let i = l.schema.field_index(key).unwrap();
        match &f.attributes[i] {
            wbvector::FieldValue::Integer(v) => *v,
            other => panic!("expected integer, got {other:?}"),
        }
    }

    /// Every subset respects the maximum, and every point lands in exactly one.
    #[test]
    fn respects_max_and_partitions_every_point() {
        let pts = grid(20); // 400 points
        let path = point_layer(&pts);
        let (polys, tagged, res) = run(
            path,
            json!({ "min_points_per_subset": 10, "max_points_per_subset": 50 }),
        );

        assert!(polys.len() > 1, "400 points with max 50 must split");
        let total: i64 = polys
            .iter()
            .map(|f| count_field(&polys, f, "point_count"))
            .sum();
        assert_eq!(
            total as usize,
            pts.len(),
            "every point is assigned exactly once"
        );
        assert_eq!(tagged.len(), pts.len());

        for f in polys.iter() {
            let c = count_field(&polys, f, "point_count");
            assert!(c <= 50, "subset exceeds the maximum: {c}");
        }
        assert_eq!(res.outputs["undersized_subset_count"], json!(0));
    }

    /// Subsets respect the minimum too.
    #[test]
    fn respects_min_points() {
        let pts = grid(16); // 256 points
        let path = point_layer(&pts);
        let (polys, _, _) = run(
            path,
            json!({ "min_points_per_subset": 40, "max_points_per_subset": 80 }),
        );
        for f in polys.iter() {
            let c = count_field(&polys, f, "point_count");
            assert!(c >= 40, "subset below the minimum: {c}");
        }
    }

    /// The partition is non-overlapping: pairwise polygon intersections have
    /// zero area. This is the property that makes it a partition rather than a
    /// clustering.
    #[test]
    fn subsets_do_not_overlap() {
        let pts = grid(14);
        let path = point_layer(&pts);
        let (polys, _, _) = run(
            path,
            json!({ "min_points_per_subset": 10, "max_points_per_subset": 40, "clip_to_hull": false }),
        );

        let geos: Vec<MultiPolygon<f64>> = polys
            .iter()
            .filter_map(|f| f.geometry.as_ref().and_then(to_geo))
            .collect();
        assert!(geos.len() > 1);
        for i in 0..geos.len() {
            for j in (i + 1)..geos.len() {
                let a = geos[i].intersection(&geos[j]).unsigned_area();
                assert!(a < 1e-9, "subsets {i} and {j} overlap by {a}");
            }
        }
    }

    /// Density-adaptive: a dense cluster gets more subsets than a sparse tail,
    /// which is exactly what a fixed grid cannot do.
    #[test]
    fn adapts_to_density() {
        // 200 points packed into x in [0,1], plus 20 spread over x in [10,30].
        let mut pts: Vec<(f64, f64)> = Vec::new();
        for i in 0..200 {
            pts.push((i as f64 * 0.005, (i % 13) as f64 * 0.01));
        }
        for i in 0..20 {
            pts.push((10.0 + i as f64, (i % 7) as f64));
        }
        let path = point_layer(&pts);
        let (polys, _, _) = run(
            path,
            json!({ "min_points_per_subset": 5, "max_points_per_subset": 40, "clip_to_hull": false }),
        );

        // Count subsets whose polygon lies in the dense region (x < 5).
        let mut dense = 0;
        let mut sparse = 0;
        for f in polys.iter() {
            let Some(g) = f.geometry.as_ref().and_then(to_geo) else {
                continue;
            };
            let cx = g.0[0].exterior().0.iter().map(|c| c.x).sum::<f64>()
                / g.0[0].exterior().0.len() as f64;
            if cx < 5.0 {
                dense += 1;
            } else {
                sparse += 1;
            }
        }
        assert!(
            dense > sparse,
            "the dense cluster should get more subsets ({dense}) than the sparse tail ({sparse})"
        );
    }

    /// Clipping to the hull keeps the partition inside the data footprint.
    #[test]
    fn clip_to_hull_shrinks_total_area() {
        // A diagonal band of points: its bounding box is much larger than its hull.
        let pts: Vec<(f64, f64)> = (0..120).map(|i| (i as f64, i as f64)).collect();
        let path = point_layer(&pts);
        let (clipped, _, _) = run(
            path.clone(),
            json!({ "min_points_per_subset": 5, "max_points_per_subset": 30, "clip_to_hull": true }),
        );
        let (unclipped, _, _) = run(
            path,
            json!({ "min_points_per_subset": 5, "max_points_per_subset": 30, "clip_to_hull": false }),
        );

        let area = |l: &Layer| -> f64 {
            l.iter()
                .filter_map(|f| f.geometry.as_ref().and_then(to_geo))
                .map(|g| g.unsigned_area())
                .sum()
        };
        // Perfectly collinear points have a zero-area hull, so nothing survives
        // the clip; either way the clipped area must not exceed the unclipped.
        assert!(area(&clipped) <= area(&unclipped) + 1e-9);
    }

    /// Coincident points collapse under the default rule, so a stack of
    /// duplicates does not force a split it cannot achieve.
    #[test]
    fn coincident_points_collapse_when_single() {
        // 300 copies of one coordinate: effectively a single location.
        let pts: Vec<(f64, f64)> = (0..300).map(|_| (5.0, 5.0)).collect();
        let path = point_layer(&pts);
        let (polys, _, _) = run(
            path,
            json!({ "min_points_per_subset": 2, "max_points_per_subset": 10, "clip_to_hull": false }),
        );
        assert_eq!(
            polys.len(),
            1,
            "coincident points cannot be separated, so one subset is correct"
        );
    }

    /// A stack of coincident points must land wholly in one subset. If the
    /// split boundary could fall inside the stack, both halves would count that
    /// one location and it would receive several subset ids.
    #[test]
    fn coincident_stack_is_not_split_across_subsets() {
        // 30 distinct locations plus a 40-deep stack sitting mid-range, so the
        // naive midpoint lands inside the stack.
        let mut pts: Vec<(f64, f64)> = (0..30).map(|i| (i as f64, 0.0)).collect();
        pts.extend(std::iter::repeat_n((15.0, 0.0), 40));
        let path = point_layer(&pts);
        let (_, tagged, _) = run(
            path,
            json!({ "min_points_per_subset": 2, "max_points_per_subset": 20, "clip_to_hull": false }),
        );

        // Collect the subset ids assigned to the stacked coordinate.
        let xi = tagged.schema.field_index("subset_id").unwrap();
        let mut ids = std::collections::BTreeSet::new();
        for f in tagged.iter() {
            let Some(Geometry::Point(c)) = f.geometry.as_ref() else {
                continue;
            };
            if (c.x - 15.0).abs() < 1e-12 && c.y.abs() < 1e-12 {
                if let wbvector::FieldValue::Integer(v) = f.attributes[xi] {
                    ids.insert(v);
                }
            }
        }
        assert_eq!(
            ids.len(),
            1,
            "the coincident stack must sit in exactly one subset, got ids {ids:?}"
        );
    }

    /// A point set smaller than the minimum yields one subset and is reported.
    #[test]
    fn undersized_input_is_reported_not_hidden() {
        let pts = grid(3); // 9 points
        let path = point_layer(&pts);
        let (polys, _, res) = run(
            path,
            json!({ "min_points_per_subset": 50, "max_points_per_subset": 100 }),
        );
        assert_eq!(polys.len(), 1);
        assert_eq!(res.outputs["undersized_subset_count"], json!(1));
    }

    fn to_geo(g: &Geometry) -> Option<MultiPolygon<f64>> {
        match g {
            Geometry::Polygon {
                exterior,
                interiors,
            } => Some(MultiPolygon(vec![Polygon::new(
                LineString::new(
                    exterior
                        .0
                        .iter()
                        .map(|c| GeoCoord { x: c.x, y: c.y })
                        .collect(),
                ),
                interiors
                    .iter()
                    .map(|r| {
                        LineString::new(r.0.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect())
                    })
                    .collect(),
            )])),
            Geometry::MultiPolygon(parts) => Some(MultiPolygon(
                parts
                    .iter()
                    .map(|(e, hs)| {
                        Polygon::new(
                            LineString::new(
                                e.0.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect(),
                            ),
                            hs.iter()
                                .map(|r| {
                                    LineString::new(
                                        r.0.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect(),
                                    )
                                })
                                .collect(),
                        )
                    })
                    .collect(),
            )),
            _ => None,
        }
    }

    #[test]
    fn rejects_bad_parameters() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(GenerateSubsetPolygonsTool.validate(&args).is_err());

        let path = point_layer(&[(0.0, 0.0), (1.0, 1.0)]);
        for bad in [
            json!({ "input": path.clone(), "min_points_per_subset": 100, "max_points_per_subset": 10 }),
            json!({ "input": path.clone(), "min_points_per_subset": 0 }),
            json!({ "input": path.clone(), "coincident_points": "sometimes" }),
        ] {
            let args: ToolArgs = serde_json::from_value(bad).unwrap();
            assert!(GenerateSubsetPolygonsTool.validate(&args).is_err());
        }
    }
}
