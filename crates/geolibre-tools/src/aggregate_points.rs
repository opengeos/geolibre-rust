//! GeoLibre tool: aggregate clusters of points into polygons.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Aggregate Points* (Cartography).
//! `vector_hex_binning` imposes an arbitrary grid and `concave_hull` produces
//! one hull for the whole layer; neither yields a polygon per *cluster* at a
//! chosen distance. This is the point analogue of the shipped
//! `aggregate_polygons` and `delineate_built_up_areas`, reusing the same
//! `geo` `BooleanOps` buffer-union machinery.
//!
//! Points within `aggregation_distance` of each other are grouped by a
//! grid-hashed union-find (single-link clustering); each cluster of at least
//! `min_points` points becomes one polygon:
//!
//! * `convex_hull` (default) — the convex hull of the cluster's points (falls
//!   back to a buffered union when the hull is degenerate: fewer than 3 points
//!   or collinear);
//! * `buffer` — the union of `aggregation_distance / 2` buffers around the
//!   points, eroded back to tighten (a filled footprint of the cluster).
//!
//! Each output polygon carries `point_count` and, for any `sum_fields`, the
//! per-cluster sum of that numeric field.

use std::collections::{BTreeMap, HashMap};

use geo::{Area, Buffer, ConvexHull, LineString, MultiPoint, MultiPolygon, Point, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

const SLIVER_AREA_EPS: f64 = 1e-9;

#[derive(Clone, Copy, PartialEq)]
enum Method {
    ConvexHull,
    Buffer,
}

pub struct AggregatePointsTool;

impl Tool for AggregatePointsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "aggregate_points",
            display_name: "Aggregate Points",
            summary: "Group points that fall within an aggregation distance of each other into cluster polygons (convex hull or buffered union), each carrying a point count and optional per-cluster field sums, like ArcGIS Aggregate Points.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer (Point or MultiPoint).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "aggregation_distance",
                    description: "Maximum distance between points for them to join the same cluster (map units).",
                    required: true,
                },
                ToolParamSpec {
                    name: "min_points",
                    description: "Minimum points a cluster needs to be output (default 2; set 1 to keep singletons).",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'convex_hull' (default) or 'buffer' (buffered union eroded back).",
                    required: false,
                },
                ToolParamSpec {
                    name: "sum_fields",
                    description: "Comma-separated numeric fields to sum per cluster (output columns 'sum_<field>').",
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

        // Resolve sum-field indices up front.
        let sum_idx: Vec<(String, usize)> = prm
            .sum_fields
            .iter()
            .map(|f| {
                layer
                    .schema
                    .field_index(f)
                    .map(|i| (f.clone(), i))
                    .ok_or_else(|| ToolError::Validation(format!("sum_field '{f}' not found")))
            })
            .collect::<Result<_, _>>()?;

        // Collect points (each vertex of a Point/MultiPoint) with feature index.
        let mut pts: Vec<(f64, f64, usize)> = Vec::new();
        for (fidx, feature) in layer.features.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            match geom {
                wbvector::Geometry::Point(c) => pts.push((c.x, c.y, fidx)),
                wbvector::Geometry::MultiPoint(cs) => {
                    for c in cs {
                        pts.push((c.x, c.y, fidx));
                    }
                }
                _ => {}
            }
        }
        if pts.is_empty() {
            return Err(ToolError::Execution(
                "no point features in input".to_string(),
            ));
        }

        ctx.progress
            .info(&format!("clustering {} point(s)", pts.len()));

        // ── Grid-hashed single-link clustering (union-find) ──────────────────────
        let d = prm.aggregation_distance;
        let cell = d.max(1e-9);
        let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
        for (i, (x, y, _)) in pts.iter().enumerate() {
            grid.entry(((x / cell).floor() as i64, (y / cell).floor() as i64))
                .or_default()
                .push(i);
        }
        let mut uf = UnionFind::new(pts.len());
        let d2 = d * d;
        for (i, (x, y, _)) in pts.iter().enumerate() {
            let (gx, gy) = ((x / cell).floor() as i64, (y / cell).floor() as i64);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    if let Some(bucket) = grid.get(&(gx + dx, gy + dy)) {
                        for &j in bucket {
                            if j <= i {
                                continue;
                            }
                            let (xj, yj, _) = pts[j];
                            if (x - xj).powi(2) + (y - yj).powi(2) <= d2 {
                                uf.union(i, j);
                            }
                        }
                    }
                }
            }
        }

        // Gather clusters by root.
        let mut clusters: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..pts.len() {
            clusters.entry(uf.find(i)).or_default().push(i);
        }

        // ── Build the output layer ───────────────────────────────────────────────
        let mut out = Layer::new("aggregated_points").with_geom_type(GeometryType::MultiPolygon);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("cluster_id", FieldType::Integer));
        out.add_field(FieldDef::new("point_count", FieldType::Integer));
        for (name, _) in &sum_idx {
            out.add_field(FieldDef::new(format!("sum_{name}"), FieldType::Float));
        }

        let radius = d / 2.0;
        let mut cluster_id = 0i64;
        let mut kept = 0usize;
        let mut clustered_points = 0usize;
        // Sort roots for deterministic output ordering.
        let mut roots: Vec<usize> = clusters.keys().cloned().collect();
        roots.sort_unstable();
        for root in roots {
            let members = &clusters[&root];
            if members.len() < prm.min_points {
                continue;
            }
            let coords: Vec<Point> = members
                .iter()
                .map(|&i| Point::new(pts[i].0, pts[i].1))
                .collect();
            let mp = build_polygon(prm.method, &coords, radius);
            if mp.unsigned_area() <= SLIVER_AREA_EPS {
                continue;
            }

            let mut attrs: Vec<(String, FieldValue)> = vec![
                ("cluster_id".to_string(), cluster_id.into()),
                ("point_count".to_string(), (members.len() as i64).into()),
            ];
            // Per-cluster field sums.
            for (name, idx) in &sum_idx {
                let s: f64 = members
                    .iter()
                    .map(|&i| {
                        layer.features[pts[i].2].attributes[*idx]
                            .as_f64()
                            .unwrap_or(0.0)
                    })
                    .sum();
                attrs.push((format!("sum_{name}"), FieldValue::Float(s)));
            }
            let attr_refs: Vec<(&str, FieldValue)> =
                attrs.iter().map(|(k, v)| (k.as_str(), v.clone())).collect();

            out.add_feature(Some(multipolygon_to_geometry(&mp)), &attr_refs)
                .map_err(|e| ToolError::Execution(format!("failed writing cluster: {e}")))?;
            cluster_id += 1;
            kept += 1;
            clustered_points += members.len();
        }

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_points".to_string(), json!(pts.len()));
        outputs.insert("cluster_count".to_string(), json!(kept));
        outputs.insert("clustered_points".to_string(), json!(clustered_points));
        Ok(ToolRunResult { outputs })
    }
}

/// Build a cluster polygon: convex hull (with buffered-union fallback for
/// degenerate hulls) or a buffered union eroded back.
fn build_polygon(method: Method, coords: &[Point], radius: f64) -> MultiPolygon {
    match method {
        Method::ConvexHull => {
            if coords.len() >= 3 {
                let hull = MultiPoint(coords.to_vec()).convex_hull();
                if hull.unsigned_area() > SLIVER_AREA_EPS {
                    return MultiPolygon(vec![hull]);
                }
            }
            // 1–2 points or collinear: buffered union so the polygon is valid.
            buffered_union(coords, radius.max(1e-6))
        }
        Method::Buffer => {
            let unioned = buffered_union(coords, radius.max(1e-6));
            // Erode back by half the radius to tighten the footprint without
            // dropping small clusters entirely.
            let eroded = unioned.buffer(-radius * 0.5);
            if eroded.unsigned_area() > SLIVER_AREA_EPS {
                eroded
            } else {
                unioned
            }
        }
    }
}

fn buffered_union(coords: &[Point], radius: f64) -> MultiPolygon {
    MultiPoint(coords.to_vec()).buffer(radius)
}

// ── Union-Find ─────────────────────────────────────────────────────────────────

struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn new(n: usize) -> UnionFind {
        UnionFind {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }
    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        match self.rank[ra].cmp(&self.rank[rb]) {
            std::cmp::Ordering::Less => self.parent[ra] = rb,
            std::cmp::Ordering::Greater => self.parent[rb] = ra,
            std::cmp::Ordering::Equal => {
                self.parent[rb] = ra;
                self.rank[ra] += 1;
            }
        }
    }
}

// ── geo -> wbvector conversion ───────────────────────────────────────────────

fn multipolygon_to_geometry(mp: &MultiPolygon) -> wbvector::Geometry {
    wbvector::Geometry::MultiPolygon(mp.0.iter().map(polygon_to_rings).collect())
}

fn polygon_to_rings(poly: &Polygon) -> (Ring, Vec<Ring>) {
    (
        linestring_to_ring(poly.exterior()),
        poly.interiors().iter().map(linestring_to_ring).collect(),
    )
}

fn linestring_to_ring(ls: &LineString) -> Ring {
    let mut coords: Vec<Coord> = ls.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
    if coords.len() >= 2 && coords.first() == coords.last() {
        coords.pop();
    }
    Ring::new(coords)
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    aggregation_distance: f64,
    min_points: usize,
    method: Method,
    sum_fields: Vec<String>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let aggregation_distance = match args.get("aggregation_distance") {
        Some(Value::Number(n)) => n.as_f64().unwrap_or(0.0),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'aggregation_distance' must be a number".into()))?,
        _ => {
            return Err(ToolError::Validation(
                "missing required numeric parameter 'aggregation_distance'".to_string(),
            ))
        }
    };
    if aggregation_distance.is_nan() || aggregation_distance <= 0.0 {
        return Err(ToolError::Validation(
            "'aggregation_distance' must be positive".to_string(),
        ));
    }
    let min_points = match args.get("min_points") {
        None | Some(Value::Null) => 2,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(2).max(1) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 2,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'min_points' must be an integer".into()))?
            .max(1),
        Some(_) => {
            return Err(ToolError::Validation(
                "'min_points' must be a number".to_string(),
            ))
        }
    };
    let method = match parse_optional_str(args, "method")?.map(|s| s.trim().to_lowercase()) {
        None => Method::ConvexHull,
        Some(s) if s.is_empty() || s == "convex_hull" || s == "convexhull" => Method::ConvexHull,
        Some(s) if s == "buffer" => Method::Buffer,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "'method' must be 'convex_hull' or 'buffer', got '{other}'"
            )))
        }
    };
    let sum_fields = parse_optional_str(args, "sum_fields")?
        .map(|s| {
            s.split(',')
                .map(str::trim)
                .filter(|f| !f.is_empty())
                .map(str::to_string)
                .collect()
        })
        .unwrap_or_default();
    Ok(Params {
        aggregation_distance,
        min_points,
        method,
        sum_fields,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Geometry};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn layer_of(pts: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("val", FieldType::Float));
        for (x, y, v) in pts {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(*x, *y))),
                &[("val", (*v).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = AggregatePointsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// Two well-separated tight clusters -> two polygons, each counting 4 points.
    #[test]
    fn two_clusters_from_two_blobs() {
        let mut pts = Vec::new();
        for (cx, cy) in [(0.0, 0.0), (1000.0, 1000.0)] {
            for (dx, dy) in [(0.0, 0.0), (5.0, 0.0), (0.0, 5.0), (5.0, 5.0)] {
                pts.push((cx + dx, cy + dy, 1.0));
            }
        }
        let input = layer_of(&pts);
        let (out, layer) = run(json!({ "input": input, "aggregation_distance": 20.0 }));
        assert_eq!(out.outputs["cluster_count"], json!(2));
        let cidx = layer.schema.field_index("point_count").unwrap();
        for f in layer.iter() {
            assert_eq!(f.attributes[cidx].as_i64().unwrap(), 4);
        }
    }

    /// A distance that spans the gap merges everything into one cluster.
    #[test]
    fn large_distance_merges_all() {
        let pts = [(0.0, 0.0, 1.0), (1000.0, 1000.0, 1.0), (0.0, 5.0, 1.0)];
        let input = layer_of(&pts);
        let (out, _l) = run(json!({ "input": input, "aggregation_distance": 2000.0 }));
        assert_eq!(out.outputs["cluster_count"], json!(1));
    }

    /// min_points drops clusters below the threshold (an isolated singleton).
    #[test]
    fn min_points_drops_small_clusters() {
        let pts = [
            (0.0, 0.0, 1.0),
            (5.0, 0.0, 1.0),
            (5.0, 5.0, 1.0),
            (5000.0, 5000.0, 1.0), // lone outlier
        ];
        let input = layer_of(&pts);
        let (out, _l) =
            run(json!({ "input": input, "aggregation_distance": 20.0, "min_points": 2 }));
        assert_eq!(out.outputs["cluster_count"], json!(1), "outlier dropped");
        assert_eq!(out.outputs["clustered_points"], json!(3));
    }

    /// sum_fields aggregates a numeric field per cluster.
    #[test]
    fn sum_fields_aggregates() {
        let pts = [(0.0, 0.0, 10.0), (5.0, 0.0, 20.0), (0.0, 5.0, 30.0)];
        let input = layer_of(&pts);
        let (_out, layer) = run(json!({
            "input": input, "aggregation_distance": 20.0, "sum_fields": "val"
        }));
        let sidx = layer.schema.field_index("sum_val").unwrap();
        let f = layer.iter().next().unwrap();
        assert!((f.attributes[sidx].as_f64().unwrap() - 60.0).abs() < 1e-6);
    }

    /// The buffer method also yields one polygon per cluster.
    #[test]
    fn buffer_method_runs() {
        let pts = [
            (0.0, 0.0, 1.0),
            (5.0, 0.0, 1.0),
            (0.0, 5.0, 1.0),
            (5.0, 5.0, 1.0),
        ];
        let input = layer_of(&pts);
        let (out, _l) = run(json!({
            "input": input, "aggregation_distance": 20.0, "method": "buffer"
        }));
        assert_eq!(out.outputs["cluster_count"], json!(1));
    }

    #[test]
    fn rejects_missing_distance() {
        let input = layer_of(&[(0.0, 0.0, 1.0)]);
        let args: ToolArgs = serde_json::from_value(json!({ "input": input })).unwrap();
        assert!(AggregatePointsTool.validate(&args).is_err());
    }
}
