//! GeoLibre tool: adaptive, density-driven partition polygons for tiling large jobs.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Create Cartographic Partitions*
//! (Cartography). Large generalization / masking / conflict-detection jobs cannot
//! be run on a whole dataset at once; they are tiled into partitions that are
//! processed one at a time. A uniform fishnet (`grid_index_features`) is a poor
//! fit when features cluster — dense tiles blow the memory budget while sparse
//! tiles waste passes. This tool instead produces *adaptive* partitions: each one
//! holds roughly the same number of features (denser where the data clusters),
//! and together they tile the input extent with no gaps or overlaps so the tiling
//! is coverage-safe.
//!
//! **Algorithm.** Each feature is reduced to a representative point (its centroid,
//! falling back to its bounding-box centre). Starting from one cell covering the
//! whole extent, the tool recursively bisects any cell that still holds more than
//! `feature_count` points. Each split cuts the cell's *longer* axis at the
//! *median* of the points inside it (a kd-tree / median-split), so every split
//! divides the points roughly in half regardless of how they are distributed —
//! that is what makes the partitions adaptive to density. Cutting at the median
//! (rather than the geometric midpoint) means a tight cluster is carved into many
//! small cells while empty space stays as one large cell. Splitting stops when a
//! cell is small enough or its points are coincident (unsplittable). The leaf
//! rectangles are emitted as the partition polygons.
//!
//! Because every split divides a rectangle into two rectangles along a line, the
//! leaves exactly tile the root extent: their union is the extent, and no two
//! overlap. Each feature is assigned (by its representative point) to exactly one
//! partition, so the partition feature counts sum to the number of placed
//! features. Every partition carries `partition_id` (0-based) and `feature_count`
//! (the number of features whose representative point falls inside it).
//!
//! The process is fully deterministic: the same input and `feature_count` always
//! yield the same partitions.

use std::collections::BTreeMap;

use geo::{Centroid, Coord as GeoCoord, LineString, MultiLineString, MultiPoint, Point, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CreateCartographicPartitionsTool;

impl Tool for CreateCartographicPartitionsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "create_cartographic_partitions",
            display_name: "Create Cartographic Partitions",
            summary: "Produce adaptive, non-overlapping partition polygons that each hold roughly a target feature count (denser where data clusters) so large generalization / masking jobs can be tiled coverage-safely — a pure-Rust counterpart of ArcGIS Create Cartographic Partitions.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (any geometry). Each feature is reduced to a representative point and assigned to one partition.",
                    required: true,
                },
                ToolParamSpec {
                    name: "feature_count",
                    description: "Target maximum number of features per partition. A cell is split while it still holds more than this many features. Default 1000.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon vector of partitions (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        let input = parse_optional_str(args, "input")?;
        if input.is_none() {
            return Err(ToolError::Validation(
                "parameter 'input' is required".to_string(),
            ));
        }
        parse_feature_count(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = parse_optional_str(args, "input")?
            .ok_or_else(|| ToolError::Validation("parameter 'input' is required".to_string()))?;
        let output = parse_optional_str(args, "output")?;
        let feature_count = parse_feature_count(args)?;

        let layer = load_input_layer(input)?;
        let epsg = layer.crs_epsg();

        // Reduce each feature to a representative point.
        let pts: Vec<(f64, f64)> = layer
            .iter()
            .filter_map(|f| f.geometry.as_ref().and_then(representative_point))
            .collect();
        if pts.is_empty() {
            return Err(ToolError::Validation(
                "input layer has no features with usable geometry".to_string(),
            ));
        }

        // Root extent covering every representative point (padded if degenerate).
        let extent = points_extent(&pts);

        // Median-split recursion into leaf cells.
        let leaves = partition(extent, pts.clone(), feature_count);

        // Emit one polygon per leaf cell.
        let mut out = new_partition_layer(epsg);
        let mut max_features = 0usize;
        let mut min_features = usize::MAX;
        for (pid, leaf) in leaves.iter().enumerate() {
            max_features = max_features.max(leaf.count);
            min_features = min_features.min(leaf.count);
            let r = &leaf.rect;
            out.push(Feature {
                fid: 0,
                geometry: Some(Geometry::polygon(
                    rect_coords(r.x0, r.y0, r.x1, r.y1),
                    vec![],
                )),
                attributes: vec![
                    FieldValue::Integer(pid as i64),
                    FieldValue::Integer(leaf.count as i64),
                ],
            });
        }
        if min_features == usize::MAX {
            min_features = 0;
        }

        let partition_count = leaves.len();
        let placed: usize = leaves.iter().map(|l| l.count).sum();
        let out_path = write_or_store_layer(out, output)?;

        ctx.progress.info(&format!(
            "{partition_count} partition(s), {placed} feature(s) placed (target {feature_count}/partition)"
        ));

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("partition_count".to_string(), json!(partition_count));
        outputs.insert("feature_count_target".to_string(), json!(feature_count));
        outputs.insert("features_placed".to_string(), json!(placed));
        outputs.insert("max_partition_features".to_string(), json!(max_features));
        outputs.insert("min_partition_features".to_string(), json!(min_features));
        Ok(ToolRunResult { outputs })
    }
}

// ── Partitioning ─────────────────────────────────────────────────────────────

#[derive(Clone, Copy)]
struct Rect {
    x0: f64,
    y0: f64,
    x1: f64,
    y1: f64,
}

struct Leaf {
    rect: Rect,
    count: usize,
}

/// Recursively median-splits `pts` within `root` until every leaf holds at most
/// `feature_count` points (or its points are coincident and cannot be split).
fn partition(root: Rect, pts: Vec<(f64, f64)>, feature_count: usize) -> Vec<Leaf> {
    let mut leaves = Vec::new();
    // Explicit stack of (cell rectangle, points inside it) to avoid deep recursion.
    let mut stack: Vec<(Rect, Vec<(f64, f64)>)> = vec![(root, pts)];

    while let Some((rect, cell_pts)) = stack.pop() {
        if cell_pts.len() <= feature_count {
            leaves.push(Leaf {
                rect,
                count: cell_pts.len(),
            });
            continue;
        }

        // Prefer splitting the longer side of the cell; fall back to the other
        // axis if the points cannot be separated along the preferred one.
        let wider_is_x = (rect.x1 - rect.x0) >= (rect.y1 - rect.y0);
        let axes = if wider_is_x {
            [Axis::X, Axis::Y]
        } else {
            [Axis::Y, Axis::X]
        };

        let mut split_done = false;
        for axis in axes {
            if let Some(sv) = choose_split(&cell_pts, axis) {
                let (mut left_pts, mut right_pts) = (Vec::new(), Vec::new());
                for &p in &cell_pts {
                    let c = match axis {
                        Axis::X => p.0,
                        Axis::Y => p.1,
                    };
                    if c <= sv {
                        left_pts.push(p);
                    } else {
                        right_pts.push(p);
                    }
                }
                let (left_rect, right_rect) = split_rect(rect, axis, sv);
                stack.push((left_rect, left_pts));
                stack.push((right_rect, right_pts));
                split_done = true;
                break;
            }
        }

        if !split_done {
            // Points are coincident on both axes: emit as an (oversized) leaf.
            leaves.push(Leaf {
                rect,
                count: cell_pts.len(),
            });
        }
    }

    // Deterministic order: sort leaves by lower-left corner.
    leaves.sort_by(|a, b| {
        a.rect
            .y0
            .partial_cmp(&b.rect.y0)
            .unwrap()
            .then(a.rect.x0.partial_cmp(&b.rect.x0).unwrap())
    });
    leaves
}

#[derive(Clone, Copy)]
enum Axis {
    X,
    Y,
}

/// Picks a split value along `axis` that separates the points into two non-empty
/// groups as evenly as possible, returning `None` when every point shares the
/// same coordinate on that axis (no clean split exists).
///
/// The value is placed strictly between two distinct, adjacent sorted
/// coordinates nearest the median, so assigning by `coord <= value` reproduces
/// the geometric split exactly (no point lands on the wrong side of the cut).
fn choose_split(pts: &[(f64, f64)], axis: Axis) -> Option<f64> {
    let mut coords: Vec<f64> = pts
        .iter()
        .map(|p| match axis {
            Axis::X => p.0,
            Axis::Y => p.1,
        })
        .collect();
    coords.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = coords.len();
    let mid = n / 2;

    // Search outward from the median index for the boundary between two distinct
    // adjacent coordinates that is closest to the middle (best balance).
    for off in 0..n {
        for &i in &[mid + off, mid.wrapping_sub(off)] {
            if i >= 1 && i < n && coords[i - 1] < coords[i] {
                return Some((coords[i - 1] + coords[i]) * 0.5);
            }
        }
    }
    None
}

/// Divides `rect` into two rectangles along `axis` at coordinate `sv`.
fn split_rect(rect: Rect, axis: Axis, sv: f64) -> (Rect, Rect) {
    match axis {
        Axis::X => (Rect { x1: sv, ..rect }, Rect { x0: sv, ..rect }),
        Axis::Y => (Rect { y1: sv, ..rect }, Rect { y0: sv, ..rect }),
    }
}

/// Bounding box of the representative points, padded so it always has area.
fn points_extent(pts: &[(f64, f64)]) -> Rect {
    let (mut x0, mut y0) = (f64::INFINITY, f64::INFINITY);
    let (mut x1, mut y1) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for &(x, y) in pts {
        x0 = x0.min(x);
        y0 = y0.min(y);
        x1 = x1.max(x);
        y1 = y1.max(y);
    }
    // Pad degenerate (zero-width/height) extents so leaf polygons stay valid.
    let px = if (x1 - x0).abs() > 0.0 {
        (x1 - x0) * 1e-6
    } else {
        0.5
    };
    let py = if (y1 - y0).abs() > 0.0 {
        (y1 - y0) * 1e-6
    } else {
        0.5
    };
    Rect {
        x0: x0 - px,
        y0: y0 - py,
        x1: x1 + px,
        y1: y1 + py,
    }
}

// ── Geometry helpers ─────────────────────────────────────────────────────────

/// A feature's representative point: its geometric centroid, falling back to the
/// centre of its bounding box when the centroid is undefined (e.g. degenerate).
fn representative_point(g: &Geometry) -> Option<(f64, f64)> {
    if let Some(gg) = to_geo_geometry(g) {
        if let Some(c) = gg.centroid() {
            if c.x().is_finite() && c.y().is_finite() {
                return Some((c.x(), c.y()));
            }
        }
    }
    let bb = g.bbox()?;
    Some(((bb.min_x + bb.max_x) * 0.5, (bb.min_y + bb.max_y) * 0.5))
}

/// Closed, counter-clockwise corner ring for an axis-aligned rectangle.
fn rect_coords(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Coord> {
    vec![
        Coord::xy(x0, y0),
        Coord::xy(x1, y0),
        Coord::xy(x1, y1),
        Coord::xy(x0, y1),
    ]
}

fn new_partition_layer(epsg: Option<u32>) -> Layer {
    let mut l = Layer::new("cartographic_partitions").with_geom_type(GeometryType::Polygon);
    if let Some(e) = epsg {
        l = l.with_crs_epsg(e);
    }
    l.add_field(FieldDef::new("partition_id", FieldType::Integer));
    l.add_field(FieldDef::new("feature_count", FieldType::Integer));
    l
}

/// Converts a `wbvector` geometry to a `geo` geometry for centroid computation.
fn to_geo_geometry(g: &Geometry) -> Option<geo::Geometry<f64>> {
    let gc = |c: &Coord| GeoCoord { x: c.x, y: c.y };
    Some(match g {
        Geometry::Point(c) => geo::Geometry::Point(Point::new(c.x, c.y)),
        Geometry::MultiPoint(cs) => geo::Geometry::MultiPoint(MultiPoint(
            cs.iter().map(|c| Point::new(c.x, c.y)).collect(),
        )),
        Geometry::LineString(cs) => {
            geo::Geometry::LineString(LineString::new(cs.iter().map(gc).collect()))
        }
        Geometry::MultiLineString(parts) => geo::Geometry::MultiLineString(MultiLineString(
            parts
                .iter()
                .map(|p| LineString::new(p.iter().map(gc).collect()))
                .collect(),
        )),
        Geometry::Polygon {
            exterior,
            interiors,
        } => geo::Geometry::Polygon(rings_to_geo(exterior, interiors)),
        Geometry::MultiPolygon(parts) => geo::Geometry::MultiPolygon(geo::MultiPolygon(
            parts.iter().map(|(e, h)| rings_to_geo(e, h)).collect(),
        )),
        Geometry::GeometryCollection(_) => return None,
    })
}

fn rings_to_geo(exterior: &wbvector::Ring, interiors: &[wbvector::Ring]) -> Polygon<f64> {
    let ring_ls = |r: &wbvector::Ring| {
        LineString::new(
            r.coords()
                .iter()
                .map(|c| GeoCoord { x: c.x, y: c.y })
                .collect(),
        )
    };
    Polygon::new(ring_ls(exterior), interiors.iter().map(ring_ls).collect())
}

// ── Parameters ───────────────────────────────────────────────────────────────

/// Parses `feature_count` (JSON number or numeric string), defaulting to 1000.
fn parse_feature_count(args: &ToolArgs) -> Result<usize, ToolError> {
    let v = match args.get("feature_count") {
        None | Some(Value::Null) => return Ok(1000),
        Some(Value::Number(n)) => n.as_f64(),
        Some(Value::String(s)) if s.trim().is_empty() => return Ok(1000),
        Some(Value::String(s)) => s.trim().parse::<f64>().ok(),
        Some(_) => None,
    };
    match v {
        Some(x) if x.is_finite() && x >= 1.0 => Ok(x as usize),
        _ => Err(ToolError::Validation(
            "parameter 'feature_count' must be an integer >= 1".to_string(),
        )),
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

    /// Builds an in-memory point layer and returns its `memory://` path.
    fn point_layer(coords: &[(f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("id", FieldType::Integer));
        for (i, &(x, y)) in coords.iter().enumerate() {
            l.add_feature(Some(Geometry::point(x, y)), &[("id", (i as i64).into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CreateCartographicPartitionsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn poly_area(g: &Geometry) -> f64 {
        if let Geometry::Polygon { exterior, .. } = g {
            let c = exterior.coords();
            let n = c.len();
            let mut a = 0.0;
            for i in 0..n {
                let j = (i + 1) % n;
                a += c[i].x * c[j].y - c[j].x * c[i].y;
            }
            (a * 0.5).abs()
        } else {
            0.0
        }
    }

    fn rect_of(f: &Feature) -> (f64, f64, f64, f64) {
        let bb = f.geometry.as_ref().unwrap().bbox().unwrap();
        (bb.min_x, bb.min_y, bb.max_x, bb.max_y)
    }

    /// Two disjoint rectangles overlap iff they overlap on both axes (open test).
    fn rects_overlap(a: (f64, f64, f64, f64), b: (f64, f64, f64, f64)) -> bool {
        let eps = 1e-9;
        a.0 < b.2 - eps && b.0 < a.2 - eps && a.1 < b.3 - eps && b.1 < a.3 - eps
    }

    /// Core property: every partition holds at most `feature_count` features, the
    /// counts sum to the input size, and partitions tile the extent without gaps
    /// or overlaps (their areas sum to the root extent area).
    #[test]
    fn partitions_respect_target_and_tile_extent() {
        // 400 points on a 20×20 lattice.
        let mut coords = Vec::new();
        for i in 0..20 {
            for j in 0..20 {
                coords.push((i as f64, j as f64));
            }
        }
        let path = point_layer(&coords);
        let (out, layer) = run(json!({ "input": path, "feature_count": 50 }));

        let cnt_i = layer.schema.field_index("feature_count").unwrap();
        // Every partition is within the target.
        for f in &layer.features {
            assert!(f.attributes[cnt_i].as_i64().unwrap() <= 50);
        }
        // Counts sum to the number of input features.
        let summed: i64 = layer
            .features
            .iter()
            .map(|f| f.attributes[cnt_i].as_i64().unwrap())
            .sum();
        assert_eq!(summed, 400);
        assert_eq!(out.outputs["features_placed"], json!(400));

        // Partitions tile the extent: no two overlap, and their areas sum to the
        // root extent area.
        let rects: Vec<_> = layer.features.iter().map(rect_of).collect();
        for a in 0..rects.len() {
            for b in (a + 1)..rects.len() {
                assert!(!rects_overlap(rects[a], rects[b]), "partitions overlap");
            }
        }
        let total_area: f64 = layer
            .features
            .iter()
            .filter_map(|f| f.geometry.as_ref())
            .map(poly_area)
            .sum();
        let ext_area = {
            let x0 = rects.iter().map(|r| r.0).fold(f64::INFINITY, f64::min);
            let y0 = rects.iter().map(|r| r.1).fold(f64::INFINITY, f64::min);
            let x1 = rects.iter().map(|r| r.2).fold(f64::NEG_INFINITY, f64::max);
            let y1 = rects.iter().map(|r| r.3).fold(f64::NEG_INFINITY, f64::max);
            (x1 - x0) * (y1 - y0)
        };
        assert!(
            (total_area - ext_area).abs() < 1e-6,
            "summed partition area {total_area} must equal extent area {ext_area}"
        );
    }

    /// Denser regions get smaller partitions: a tight cluster yields more, smaller
    /// cells than an equally-populous but sparse region.
    #[test]
    fn adapts_to_density() {
        let mut coords = Vec::new();
        // Dense cluster of 100 points in a tiny [0,1]² box.
        for i in 0..10 {
            for j in 0..10 {
                coords.push((i as f64 * 0.1, j as f64 * 0.1));
            }
        }
        // Sparse 100 points spread across a wide [100,200]² box.
        for i in 0..10 {
            for j in 0..10 {
                coords.push((100.0 + i as f64 * 10.0, 100.0 + j as f64 * 10.0));
            }
        }
        let path = point_layer(&coords);
        let (_out, layer) = run(json!({ "input": path, "feature_count": 25 }));
        // 200 points / 25 target -> at least 8 partitions.
        assert!(layer.features.len() >= 8);
        // The cluster region (small coords) is covered by partitions far smaller
        // than those over the sparse region.
        let cluster_area: f64 = layer
            .features
            .iter()
            .filter(|f| rect_of(f).2 < 50.0)
            .filter_map(|f| f.geometry.as_ref())
            .map(poly_area)
            .sum();
        let sparse_area: f64 = layer
            .features
            .iter()
            .filter(|f| rect_of(f).0 > 50.0)
            .filter_map(|f| f.geometry.as_ref())
            .map(poly_area)
            .sum();
        assert!(
            cluster_area < sparse_area,
            "cluster partitions ({cluster_area}) should be smaller than sparse ({sparse_area})"
        );
    }

    /// A dataset already under the target yields a single partition covering all.
    #[test]
    fn small_input_single_partition() {
        let path = point_layer(&[(0.0, 0.0), (1.0, 1.0), (2.0, 0.5)]);
        let (out, layer) = run(json!({ "input": path, "feature_count": 1000 }));
        assert_eq!(out.outputs["partition_count"], json!(1));
        assert_eq!(layer.features.len(), 1);
        let cnt_i = layer.schema.field_index("feature_count").unwrap();
        assert_eq!(layer.features[0].attributes[cnt_i].as_i64(), Some(3));
    }

    /// Coincident points cannot be split: they stay in one (oversized) partition
    /// instead of looping forever.
    #[test]
    fn coincident_points_dont_loop() {
        let path = point_layer(&[(5.0, 5.0), (5.0, 5.0), (5.0, 5.0), (5.0, 5.0)]);
        let (out, layer) = run(json!({ "input": path, "feature_count": 1 }));
        assert_eq!(out.outputs["partition_count"], json!(1));
        let cnt_i = layer.schema.field_index("feature_count").unwrap();
        assert_eq!(layer.features[0].attributes[cnt_i].as_i64(), Some(4));
    }

    /// Polygon input is reduced to centroids and assigned coverage-safely.
    #[test]
    fn polygon_input_uses_centroids() {
        let mut l = Layer::new("polys")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("id", FieldType::Integer));
        for i in 0..8 {
            let x = i as f64 * 10.0;
            l.add_feature(
                Some(Geometry::polygon(rect_coords(x, 0.0, x + 5.0, 5.0), vec![])),
                &[("id", (i as i64).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        let path = memory_store::make_vector_memory_path(&id);
        let (out, _layer) = run(json!({ "input": path, "feature_count": 3 }));
        assert_eq!(out.outputs["features_placed"], json!(8));
        assert!(out.outputs["partition_count"].as_i64().unwrap() >= 3);
    }

    /// The output inherits the input CRS.
    #[test]
    fn preserves_crs() {
        let path = point_layer(&[(0.0, 0.0), (1.0, 1.0)]);
        let (_out, layer) = run(json!({ "input": path }));
        assert_eq!(layer.crs_epsg(), Some(3857));
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CreateCartographicPartitionsTool.validate(&args)
        };
        // Missing required input.
        assert!(bad(json!({ "feature_count": 100 })).is_err());
        // feature_count below 1.
        assert!(bad(json!({ "input": "x.geojson", "feature_count": 0 })).is_err());
        // Non-numeric feature_count.
        assert!(bad(json!({ "input": "x.geojson", "feature_count": "lots" })).is_err());
        // A valid call validates.
        assert!(bad(json!({ "input": "x.geojson", "feature_count": 500 })).is_ok());
    }
}
