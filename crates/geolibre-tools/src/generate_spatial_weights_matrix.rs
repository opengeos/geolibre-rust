//! GeoLibre tool: build a reusable spatial-weights table.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Generate Spatial Weights Matrix*
//! (Spatial Statistics). Every geolibre spatial-stats tool (`global_morans_i`,
//! `getis_ord_gi_star`, `local_morans_i_lisa`, …) rebuilds a neighbour
//! structure internally and throws it away. This tool emits that structure as
//! an inspectable, shareable long-format table — one row per neighbouring pair:
//!
//! * `origin_id`   — the focal feature's id,
//! * `neighbor_id` — a neighbour's id,
//! * `weight`      — the spatial weight linking them.
//!
//! The neighbour set is built from one of six conceptualizations:
//!
//! * `knn` — the k nearest features (kd-tree), weight 1,
//! * `fixed_distance_band` — every feature within `threshold_distance`, weight 1,
//! * `inverse_distance` — every feature within `threshold_distance` (or all
//!   others when unset), weight `1 / d^exponent`,
//! * `contiguity_edges` — polygons sharing a border segment (rook),
//! * `contiguity_edges_corners` — polygons sharing a border **or** a corner
//!   vertex (queen),
//! * `delaunay` — features adjacent in the Delaunay triangulation of their
//!   representative points.
//!
//! With `row_standardization`, each origin's weights are divided by their row
//! sum so every origin's weights total 1. Point-based methods use a
//! representative point per feature (centroid for lines/polygons); contiguity
//! requires polygon geometry. Output is a geometry-less attribute table (or a
//! CSV when the path ends in `.csv`), not ESRI's proprietary `.swm`.
//! Deterministic.

use std::collections::{BTreeMap, HashMap, HashSet};

use geo::{Centroid, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Feature, FieldDef, FieldType, FieldValue, Geometry, Layer, Ring};

use crate::common::write_text_output;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct GenerateSpatialWeightsMatrixTool;

impl Tool for GenerateSpatialWeightsMatrixTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "generate_spatial_weights_matrix",
            display_name: "Generate Spatial Weights Matrix",
            summary: "Emit a reusable long-format neighbour/weights table (origin_id, neighbor_id, weight) from a chosen conceptualization — KNN, fixed distance band, inverse distance, contiguity (edges / edges+corners), or Delaunay — with optional row standardization, like ArcGIS Generate Spatial Weights Matrix.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (points/lines/polygons; contiguity requires polygons).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output table path — a CSV (extension .csv) or a geometry-less vector table. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "id_field",
                    description: "Field holding each feature's unique id. Default: the 0-based feature index.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "Conceptualization of spatial relationships: 'knn', 'fixed_distance_band', 'inverse_distance', 'contiguity_edges', 'contiguity_edges_corners', or 'delaunay'. Default 'knn'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "number_of_neighbors",
                    description: "Number of nearest neighbours k (method = knn). Default 8.",
                    required: false,
                },
                ToolParamSpec {
                    name: "threshold_distance",
                    description: "Distance band cutoff in map units. Required for fixed_distance_band; optional cutoff for inverse_distance (unset = all other features).",
                    required: false,
                },
                ToolParamSpec {
                    name: "exponent",
                    description: "Distance decay exponent p for inverse_distance (weight = 1/d^p). Default 1.",
                    required: false,
                },
                ToolParamSpec {
                    name: "row_standardization",
                    description: "Divide each origin's weights by their row sum so every origin's weights total 1. Default false.",
                    required: false,
                },
                ToolParamSpec {
                    name: "snap_tolerance",
                    description: "Quantize vertices onto a grid of this size (CRS units) before matching borders for contiguity. Default 0 (exact coordinates).",
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
        let n = layer.features.len();

        // Resolve the per-feature id strings.
        let id_idx = match &prm.id_field {
            Some(f) => Some(
                layer
                    .schema
                    .field_index(f)
                    .ok_or_else(|| ToolError::Validation(format!("id_field '{f}' not found")))?,
            ),
            None => None,
        };
        let ids: Vec<String> = (0..n)
            .map(|i| match id_idx {
                Some(fi) => field_key(&layer.features[i].attributes[fi]),
                None => i.to_string(),
            })
            .collect();

        // adjacency[i] = list of (neighbour index, raw weight).
        let adjacency: Vec<Vec<(usize, f64)>> = match prm.method {
            Method::Knn => neighbors_knn(&layer, prm.number_of_neighbors)?,
            Method::FixedDistanceBand => {
                neighbors_distance_band(&layer, prm.threshold_distance.unwrap(), false, 1.0)?
            }
            Method::InverseDistance => {
                neighbors_inverse_distance(&layer, prm.threshold_distance, prm.exponent)?
            }
            Method::ContiguityEdges => neighbors_contiguity(&layer, false, prm.snap_tolerance)?,
            Method::ContiguityEdgesCorners => {
                neighbors_contiguity(&layer, true, prm.snap_tolerance)?
            }
            Method::Delaunay => neighbors_delaunay(&layer)?,
        };

        // Optional row standardization: each origin's weights sum to 1.
        let mut adjacency = adjacency;
        if prm.row_standardization {
            for row in adjacency.iter_mut() {
                let sum: f64 = row.iter().map(|(_, w)| *w).sum();
                if sum > 0.0 {
                    for (_, w) in row.iter_mut() {
                        *w /= sum;
                    }
                }
            }
        }

        let pair_count: usize = adjacency.iter().map(Vec::len).sum();
        let isolates = adjacency.iter().filter(|r| r.is_empty()).count();
        ctx.progress.info(&format!(
            "{n} feature(s); {pair_count} neighbour link(s); {isolates} isolate(s)"
        ));

        // ── Emit the long-format weights table ────────────────────────────────
        let mut table = Layer::new("spatial_weights");
        table.add_field(FieldDef::new("origin_id", FieldType::Text));
        table.add_field(FieldDef::new("neighbor_id", FieldType::Text));
        table.add_field(FieldDef::new("weight", FieldType::Float));

        let mut csv = String::from("origin_id,neighbor_id,weight\n");
        let mut rows = 0usize;
        for (i, row) in adjacency.iter().enumerate() {
            // Deterministic ordering: by neighbour index.
            let mut sorted = row.clone();
            sorted.sort_by(|a, b| a.0.cmp(&b.0));
            for (j, w) in sorted {
                table.push(Feature {
                    fid: 0,
                    geometry: None,
                    attributes: vec![
                        FieldValue::Text(ids[i].clone()),
                        FieldValue::Text(ids[j].clone()),
                        FieldValue::Float(w),
                    ],
                });
                csv.push_str(&format!("{},{},{}\n", ids[i], ids[j], w));
                rows += 1;
            }
        }

        let out_path = match output {
            Some(path) if path.to_ascii_lowercase().ends_with(".csv") => {
                write_text_output(&csv, path)?;
                path.to_string()
            }
            other => write_or_store_layer(table, other)?,
        };

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("link_count".to_string(), json!(pair_count));
        outputs.insert("isolate_count".to_string(), json!(isolates));
        outputs.insert("row_count".to_string(), json!(rows));
        outputs.insert(
            "mean_neighbors".to_string(),
            json!(if n > 0 {
                pair_count as f64 / n as f64
            } else {
                0.0
            }),
        );
        Ok(ToolRunResult { outputs })
    }
}

// ── Neighbour builders ────────────────────────────────────────────────────────

/// Representative point per feature (`None` when no usable geometry).
fn representative_points(layer: &Layer) -> Vec<Option<(f64, f64)>> {
    layer
        .features
        .iter()
        .map(|f| f.geometry.as_ref().and_then(representative_point))
        .collect()
}

/// Build a kd-tree over the representative points that exist.
fn build_tree(reps: &[Option<(f64, f64)>]) -> Result<KdTree<f64, usize, [f64; 2]>, ToolError> {
    let mut tree: KdTree<f64, usize, [f64; 2]> = KdTree::new(2);
    for (i, r) in reps.iter().enumerate() {
        if let Some((x, y)) = r {
            tree.add([*x, *y], i)
                .map_err(|e| ToolError::Execution(format!("kd-tree insert failed: {e:?}")))?;
        }
    }
    Ok(tree)
}

fn neighbors_knn(layer: &Layer, k: usize) -> Result<Vec<Vec<(usize, f64)>>, ToolError> {
    let reps = representative_points(layer);
    if reps.iter().all(Option::is_none) {
        return Err(ToolError::Execution(
            "input has no usable point geometry for knn".to_string(),
        ));
    }
    let tree = build_tree(&reps)?;
    let n = reps.len();
    let mut out = vec![Vec::new(); n];
    for (i, rep) in reps.iter().enumerate() {
        let Some((x, y)) = rep else { continue };
        // Ask for k+1 so a feature can drop its own self-match.
        let found = tree
            .nearest(&[*x, *y], k + 1, &squared_euclidean)
            .map_err(|e| ToolError::Execution(format!("kd-tree query failed: {e:?}")))?;
        for (_d2, &j) in found.into_iter().filter(|(_, &j)| j != i).take(k) {
            out[i].push((j, 1.0));
        }
    }
    Ok(out)
}

fn neighbors_distance_band(
    layer: &Layer,
    threshold: f64,
    inverse: bool,
    exponent: f64,
) -> Result<Vec<Vec<(usize, f64)>>, ToolError> {
    let reps = representative_points(layer);
    if reps.iter().all(Option::is_none) {
        return Err(ToolError::Execution(
            "input has no usable point geometry".to_string(),
        ));
    }
    let tree = build_tree(&reps)?;
    let r2 = threshold * threshold;
    let n = reps.len();
    let mut out = vec![Vec::new(); n];
    for (i, rep) in reps.iter().enumerate() {
        let Some((x, y)) = rep else { continue };
        let found = tree
            .within(&[*x, *y], r2, &squared_euclidean)
            .map_err(|e| ToolError::Execution(format!("kd-tree query failed: {e:?}")))?;
        for (d2, &j) in found {
            if j == i {
                continue;
            }
            let w = if inverse {
                let d = d2.sqrt();
                if d > 0.0 {
                    1.0 / d.powf(exponent)
                } else {
                    continue;
                }
            } else {
                1.0
            };
            out[i].push((j, w));
        }
    }
    Ok(out)
}

fn neighbors_inverse_distance(
    layer: &Layer,
    threshold: Option<f64>,
    exponent: f64,
) -> Result<Vec<Vec<(usize, f64)>>, ToolError> {
    if let Some(t) = threshold {
        return neighbors_distance_band(layer, t, true, exponent);
    }
    // No cutoff: every other feature contributes 1/d^p (dense; O(n^2)).
    let reps = representative_points(layer);
    if reps.iter().all(Option::is_none) {
        return Err(ToolError::Execution(
            "input has no usable point geometry".to_string(),
        ));
    }
    let n = reps.len();
    let mut out = vec![Vec::new(); n];
    for i in 0..n {
        let Some((xi, yi)) = reps[i] else { continue };
        for (j, rep) in reps.iter().enumerate() {
            if j == i {
                continue;
            }
            let Some((xj, yj)) = rep else { continue };
            let d = (xi - xj).hypot(yi - yj);
            if d > 0.0 {
                out[i].push((j, 1.0 / d.powf(exponent)));
            }
        }
    }
    Ok(out)
}

fn neighbors_delaunay(layer: &Layer) -> Result<Vec<Vec<(usize, f64)>>, ToolError> {
    let reps = representative_points(layer);
    // Collect the indices that have a representative point, keeping a mapping.
    let present: Vec<usize> = reps
        .iter()
        .enumerate()
        .filter_map(|(i, r)| r.map(|_| i))
        .collect();
    let pts: Vec<(f64, f64)> = present.iter().map(|&i| reps[i].unwrap()).collect();
    let n = reps.len();
    let mut sets: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    for t in delaunay(&pts) {
        for &(a, b) in &[(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let (ga, gb) = (present[a], present[b]);
            sets[ga].insert(gb);
            sets[gb].insert(ga);
        }
    }
    Ok(sets
        .into_iter()
        .map(|s| s.into_iter().map(|j| (j, 1.0)).collect())
        .collect())
}

/// Polygon contiguity via shared edges (rook) and, when `include_corners`,
/// shared vertices too (queen). Symmetric by construction.
fn neighbors_contiguity(
    layer: &Layer,
    include_corners: bool,
    snap: f64,
) -> Result<Vec<Vec<(usize, f64)>>, ToolError> {
    let n = layer.features.len();
    let mut edge_feats: HashMap<(Key, Key), HashSet<usize>> = HashMap::new();
    let mut vert_feats: HashMap<Key, HashSet<usize>> = HashMap::new();
    let mut any_poly = false;

    for (fidx, feature) in layer.features.iter().enumerate() {
        let Some(geom) = feature.geometry.as_ref() else {
            continue;
        };
        let rings = polygon_rings(geom, snap);
        if rings.is_empty() {
            continue;
        }
        any_poly = true;
        for ring in &rings {
            let m = ring.len();
            for i in 0..m {
                let a = ring[i];
                let b = ring[(i + 1) % m];
                vert_feats.entry(key(a)).or_default().insert(fidx);
                if key(a) == key(b) {
                    continue;
                }
                edge_feats.entry(edge_key(a, b)).or_default().insert(fidx);
            }
        }
    }
    if !any_poly {
        return Err(ToolError::Execution(
            "contiguity requires polygon geometry, but none was found".to_string(),
        ));
    }

    let mut sets: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    let link = |feats: &HashSet<usize>, sets: &mut Vec<HashSet<usize>>| {
        if feats.len() < 2 {
            return;
        }
        let list: Vec<usize> = feats.iter().copied().collect();
        for i in 0..list.len() {
            for j in (i + 1)..list.len() {
                sets[list[i]].insert(list[j]);
                sets[list[j]].insert(list[i]);
            }
        }
    };
    for feats in edge_feats.values() {
        link(feats, &mut sets);
    }
    if include_corners {
        for feats in vert_feats.values() {
            link(feats, &mut sets);
        }
    }

    Ok(sets
        .into_iter()
        .map(|s| s.into_iter().map(|j| (j, 1.0)).collect())
        .collect())
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

// ── Contiguity keys and rings ─────────────────────────────────────────────────

#[derive(Clone, Copy)]
struct P {
    x: f64,
    y: f64,
}

type Key = (u64, u64);

fn key(p: P) -> Key {
    (p.x.to_bits(), p.y.to_bits())
}

fn edge_key(a: P, b: P) -> (Key, Key) {
    let (ka, kb) = (key(a), key(b));
    if ka <= kb {
        (ka, kb)
    } else {
        (kb, ka)
    }
}

fn canonical(x: f64, y: f64, snap: f64) -> P {
    if snap > 0.0 {
        P {
            x: (x / snap).round() * snap,
            y: (y / snap).round() * snap,
        }
    } else {
        P { x, y }
    }
}

/// All rings (exterior + interiors) of a polygon geometry as canonical vertex
/// chains without the closing duplicate.
fn polygon_rings(geom: &Geometry, snap: f64) -> Vec<Vec<P>> {
    let ring_pts = |ring: &Ring| -> Vec<P> {
        let mut pts: Vec<P> = Vec::with_capacity(ring.len());
        for c in ring.coords() {
            let p = canonical(c.x, c.y, snap);
            if pts.last().is_none_or(|l| key(*l) != key(p)) {
                pts.push(p);
            }
        }
        while pts.len() >= 2 && key(pts[0]) == key(*pts.last().unwrap()) {
            pts.pop();
        }
        pts
    };
    let mut out = Vec::new();
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            out.push(ring_pts(exterior));
            out.extend(interiors.iter().map(&ring_pts));
        }
        Geometry::MultiPolygon(parts) => {
            for (ext, holes) in parts {
                out.push(ring_pts(ext));
                out.extend(holes.iter().map(&ring_pts));
            }
        }
        _ => {}
    }
    out.retain(|r| r.len() >= 3);
    out
}

fn field_key(fv: &FieldValue) -> String {
    if let Some(i) = fv.as_i64() {
        i.to_string()
    } else if let Some(f) = fv.as_f64() {
        format!("{f}")
    } else {
        fv.as_str().unwrap_or("").to_string()
    }
}

// ── Delaunay triangulation (Bowyer–Watson) ────────────────────────────────────

/// Returns triangles as index triples into `pts` (super-triangle removed).
fn delaunay(pts: &[(f64, f64)]) -> Vec<[usize; 3]> {
    let n = pts.len();
    if n < 3 {
        return Vec::new();
    }
    let (mut minx, mut miny, mut maxx, mut maxy) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
    for &(x, y) in pts {
        minx = minx.min(x);
        miny = miny.min(y);
        maxx = maxx.max(x);
        maxy = maxy.max(y);
    }
    let dmax = (maxx - minx).max(maxy - miny).max(1.0) * 20.0;
    let (mx, my) = ((minx + maxx) / 2.0, (miny + maxy) / 2.0);
    let mut verts: Vec<(f64, f64)> = pts.to_vec();
    let s0 = verts.len();
    verts.push((mx - dmax, my - dmax));
    verts.push((mx + dmax, my - dmax));
    verts.push((mx, my + dmax));
    let mut tris: Vec<[usize; 3]> = vec![[s0, s0 + 1, s0 + 2]];

    for i in 0..n {
        let p = verts[i];
        let mut bad: Vec<usize> = Vec::new();
        for (ti, t) in tris.iter().enumerate() {
            if in_circumcircle(p, verts[t[0]], verts[t[1]], verts[t[2]]) {
                bad.push(ti);
            }
        }
        let mut boundary: Vec<(usize, usize)> = Vec::new();
        for &ti in &bad {
            let t = tris[ti];
            for e in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                let shared = bad.iter().any(|&oj| {
                    oj != ti && {
                        let o = tris[oj];
                        let edges = [(o[0], o[1]), (o[1], o[2]), (o[2], o[0])];
                        edges
                            .iter()
                            .any(|&(u, v)| (u == e.0 && v == e.1) || (u == e.1 && v == e.0))
                    }
                });
                if !shared {
                    boundary.push(e);
                }
            }
        }
        bad.sort_unstable_by(|a, b| b.cmp(a));
        for &ti in &bad {
            tris.swap_remove(ti);
        }
        for (u, v) in boundary {
            tris.push([u, v, i]);
        }
    }
    tris.retain(|t| t.iter().all(|&v| v < s0));
    tris
}

fn in_circumcircle(p: (f64, f64), a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> bool {
    let ax = a.0 - p.0;
    let ay = a.1 - p.1;
    let bx = b.0 - p.0;
    let by = b.1 - p.1;
    let cx = c.0 - p.0;
    let cy = c.1 - p.1;
    let det = (ax * ax + ay * ay) * (bx * cy - cx * by) - (bx * bx + by * by) * (ax * cy - cx * ay)
        + (cx * cx + cy * cy) * (ax * by - bx * ay);
    let area2 = (b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1);
    if area2 > 0.0 {
        det > 0.0
    } else {
        det < 0.0
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Knn,
    FixedDistanceBand,
    InverseDistance,
    ContiguityEdges,
    ContiguityEdgesCorners,
    Delaunay,
}

struct Params {
    id_field: Option<String>,
    method: Method,
    number_of_neighbors: usize,
    threshold_distance: Option<f64>,
    exponent: f64,
    row_standardization: bool,
    snap_tolerance: f64,
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
        Some(Value::Number(n)) => Ok(n.as_f64().map(|v| v != 0.0).unwrap_or(false)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(false),
            "true" | "1" | "yes" => Ok(true),
            "false" | "0" | "no" => Ok(false),
            _ => Err(ToolError::Validation(format!("'{key}' must be a boolean"))),
        },
        Some(_) => Err(ToolError::Validation(format!("'{key}' must be a boolean"))),
    }
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let id_field = parse_optional_str(args, "id_field")?.map(str::to_string);

    let method = match parse_optional_str(args, "method")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("knn") | Some("k_nearest_neighbors") => Method::Knn,
        Some("fixed_distance_band") | Some("distance_band") => Method::FixedDistanceBand,
        Some("inverse_distance") | Some("idw") => Method::InverseDistance,
        Some("contiguity_edges") | Some("edges") | Some("rook") => Method::ContiguityEdges,
        Some("contiguity_edges_corners") | Some("edges_corners") | Some("queen") => {
            Method::ContiguityEdgesCorners
        }
        Some("delaunay") | Some("triangulation") => Method::Delaunay,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "'method' must be one of knn, fixed_distance_band, inverse_distance, contiguity_edges, contiguity_edges_corners, delaunay (got '{other}')"
            )))
        }
    };

    let number_of_neighbors = match args.get("number_of_neighbors") {
        None | Some(Value::Null) => 8,
        Some(Value::Number(n)) => n.as_u64().ok_or_else(|| {
            ToolError::Validation("'number_of_neighbors' must be a positive integer".into())
        })? as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 8,
        Some(Value::String(s)) => s.trim().parse::<usize>().map_err(|_| {
            ToolError::Validation("'number_of_neighbors' must be an integer".into())
        })?,
        Some(_) => {
            return Err(ToolError::Validation(
                "'number_of_neighbors' must be a number".into(),
            ))
        }
    };
    if method == Method::Knn && number_of_neighbors == 0 {
        return Err(ToolError::Validation(
            "'number_of_neighbors' must be >= 1 for knn".to_string(),
        ));
    }

    let threshold_distance = parse_optional_f64(args, "threshold_distance")?;
    if let Some(t) = threshold_distance {
        if !(t.is_finite() && t > 0.0) {
            return Err(ToolError::Validation(
                "'threshold_distance' must be positive".to_string(),
            ));
        }
    }
    if method == Method::FixedDistanceBand && threshold_distance.is_none() {
        return Err(ToolError::Validation(
            "'threshold_distance' is required for fixed_distance_band".to_string(),
        ));
    }

    let exponent = parse_optional_f64(args, "exponent")?.unwrap_or(1.0);
    if !(exponent.is_finite() && exponent > 0.0) {
        return Err(ToolError::Validation(
            "'exponent' must be a positive number".to_string(),
        ));
    }

    let snap_tolerance = parse_optional_f64(args, "snap_tolerance")?.unwrap_or(0.0);
    if !(snap_tolerance >= 0.0 && snap_tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "'snap_tolerance' must be a non-negative number".to_string(),
        ));
    }

    Ok(Params {
        id_field,
        method,
        number_of_neighbors,
        threshold_distance,
        exponent,
        row_standardization: parse_optional_bool(args, "row_standardization")?,
        snap_tolerance,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Coord, GeometryType};

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

    fn square(x0: f64, y0: f64, s: f64) -> Geometry {
        Geometry::polygon(
            vec![
                Coord::xy(x0, y0),
                Coord::xy(x0 + s, y0),
                Coord::xy(x0 + s, y0 + s),
                Coord::xy(x0, y0 + s),
            ],
            vec![],
        )
    }

    fn poly_layer(named: &[(&str, Geometry)]) -> String {
        let mut l = Layer::new("polys")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        for (n, g) in named {
            l.add_feature(Some(g.clone()), &[("name", (*n).into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = GenerateSpatialWeightsMatrixTool.run(&args, &ctx()).unwrap();
        let table = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, table)
    }

    fn rows(table: &Layer) -> Vec<(String, String, f64)> {
        let oi = table.schema.field_index("origin_id").unwrap();
        let ni = table.schema.field_index("neighbor_id").unwrap();
        let wi = table.schema.field_index("weight").unwrap();
        table
            .iter()
            .map(|f| {
                (
                    f.attributes[oi].as_str().unwrap().to_string(),
                    f.attributes[ni].as_str().unwrap().to_string(),
                    f.attributes[wi].as_f64().unwrap(),
                )
            })
            .collect()
    }

    /// KNN gives exactly k neighbours per origin.
    #[test]
    fn knn_gives_exactly_k_neighbors() {
        // 5 collinear points; k=2.
        let input = point_layer(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0)]);
        let (out, table) = run(json!({
            "input": input, "id_field": "id", "method": "knn", "number_of_neighbors": 2
        }));
        assert_eq!(out.outputs["row_count"], json!(10)); // 5 origins * 2
        let mut per_origin: BTreeMap<String, usize> = BTreeMap::new();
        for (o, _n, _w) in rows(&table) {
            *per_origin.entry(o).or_default() += 1;
        }
        assert!(per_origin.values().all(|&c| c == 2), "each origin has k=2");
    }

    /// Row standardization makes each origin's weights sum to 1.
    #[test]
    fn row_standardized_weights_sum_to_one() {
        let input = point_layer(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (5.0, 0.0)]);
        let (_o, table) = run(json!({
            "input": input, "id_field": "id", "method": "knn",
            "number_of_neighbors": 2, "row_standardization": true
        }));
        let mut sums: BTreeMap<String, f64> = BTreeMap::new();
        for (o, _n, w) in rows(&table) {
            *sums.entry(o).or_default() += w;
        }
        for (o, s) in sums {
            assert!((s - 1.0).abs() < 1e-9, "origin {o} sums to {s}, not 1");
        }
    }

    /// Contiguity (queen) is symmetric: i~j implies j~i.
    #[test]
    fn contiguity_is_symmetric() {
        // 2x2 grid of squares; centre-touching produces edge + corner links.
        let input = poly_layer(&[
            ("A", square(0.0, 0.0, 10.0)),
            ("B", square(10.0, 0.0, 10.0)),
            ("C", square(0.0, 10.0, 10.0)),
            ("D", square(10.0, 10.0, 10.0)),
        ]);
        let (_o, table) = run(json!({
            "input": input, "id_field": "name", "method": "contiguity_edges_corners"
        }));
        let set: HashSet<(String, String)> =
            rows(&table).into_iter().map(|(o, n, _)| (o, n)).collect();
        for (o, n) in &set {
            assert!(
                set.contains(&(n.clone(), o.clone())),
                "missing reverse link {n}->{o}"
            );
        }
        // In a 2x2 grid every square touches all three others (queen).
        assert_eq!(set.len(), 12, "4 squares * 3 neighbours each");
    }

    /// Edge-only (rook) contiguity excludes diagonal (corner-only) neighbours.
    #[test]
    fn rook_excludes_diagonal() {
        let input = poly_layer(&[
            ("A", square(0.0, 0.0, 10.0)),
            ("B", square(10.0, 0.0, 10.0)),
            ("C", square(0.0, 10.0, 10.0)),
            ("D", square(10.0, 10.0, 10.0)),
        ]);
        let (_o, table) = run(json!({
            "input": input, "id_field": "name", "method": "contiguity_edges"
        }));
        let set: HashSet<(String, String)> =
            rows(&table).into_iter().map(|(o, n, _)| (o, n)).collect();
        // A(0,0) and D(10,10) touch only at a corner -> not rook neighbours.
        assert!(!set.contains(&("A".into(), "D".into())));
        assert!(!set.contains(&("B".into(), "C".into())));
        // A-B share an edge.
        assert!(set.contains(&("A".into(), "B".into())));
        assert_eq!(set.len(), 8, "4 edges, both directions");
    }

    /// Fixed distance band links every feature within the threshold, weight 1.
    #[test]
    fn distance_band_within_threshold() {
        let input = point_layer(&[(0.0, 0.0), (3.0, 0.0), (100.0, 0.0)]);
        let (out, table) = run(json!({
            "input": input, "id_field": "id",
            "method": "fixed_distance_band", "threshold_distance": 5.0
        }));
        // 0<->1 within 5; 2 isolated.
        assert_eq!(out.outputs["link_count"], json!(2));
        assert_eq!(out.outputs["isolate_count"], json!(1));
        assert!(rows(&table)
            .iter()
            .all(|(_, _, w)| (*w - 1.0).abs() < 1e-12));
    }

    /// Inverse-distance weights decay with distance (1/d).
    #[test]
    fn inverse_distance_decays() {
        let input = point_layer(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]);
        let (_o, table) = run(json!({
            "input": input, "id_field": "id", "method": "inverse_distance"
        }));
        // From origin 0: neighbour at d=1 weight 1.0, neighbour at d=2 weight 0.5.
        let r: Vec<_> = rows(&table)
            .into_iter()
            .filter(|(o, _, _)| o == "0")
            .collect();
        let w1 = r.iter().find(|(_, n, _)| n == "1").unwrap().2;
        let w2 = r.iter().find(|(_, n, _)| n == "2").unwrap().2;
        assert!((w1 - 1.0).abs() < 1e-9);
        assert!((w2 - 0.5).abs() < 1e-9);
    }

    /// Delaunay neighbours are symmetric and non-empty for a simple point set.
    #[test]
    fn delaunay_symmetric() {
        let input = point_layer(&[(0.0, 0.0), (10.0, 0.0), (5.0, 8.0), (5.0, 3.0)]);
        let (out, table) = run(json!({
            "input": input, "id_field": "id", "method": "delaunay"
        }));
        assert!(out.outputs["link_count"].as_u64().unwrap() >= 6);
        let set: HashSet<(String, String)> =
            rows(&table).into_iter().map(|(o, n, _)| (o, n)).collect();
        for (o, n) in &set {
            assert!(set.contains(&(n.clone(), o.clone())), "asymmetric {o}->{n}");
        }
    }

    #[test]
    fn rejects_missing_input() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(GenerateSpatialWeightsMatrixTool.validate(&args).is_err());
    }

    #[test]
    fn rejects_distance_band_without_threshold() {
        let input = point_layer(&[(0.0, 0.0)]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "method": "fixed_distance_band" }))
                .unwrap();
        assert!(GenerateSpatialWeightsMatrixTool.validate(&args).is_err());
    }

    #[test]
    fn rejects_bad_method() {
        let input = point_layer(&[(0.0, 0.0)]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "method": "nonsense" })).unwrap();
        assert!(GenerateSpatialWeightsMatrixTool.validate(&args).is_err());
    }
}
