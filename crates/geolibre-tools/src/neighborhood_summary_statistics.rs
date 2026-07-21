//! GeoLibre tool: per-feature summary statistics of neighbouring attribute values.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Neighborhood Summary Statistics*
//! (Spatial Statistics). Every weights-based workflow — spatial-regression
//! preparation, smoothing, anomaly screening — starts from neighbour summaries
//! (the classic "spatial lag" columns), yet the bundled suite computes only
//! global/local indices (`global_morans_i`, `getis_ord_gi_star`) and never
//! exports the neighbour statistics themselves.
//!
//! For every feature this computes mean / median / std / min / max / sum of the
//! chosen numeric `fields` over that feature's neighbours, defined by one of:
//!
//! * `knn` — the k nearest features by representative-point distance;
//! * `distance_band` — every feature within a distance;
//! * `contiguity` — polygons sharing an edge (the shared-edge hashing from
//!   `polygon_neighbors`).
//!
//! With `weights = inverse_distance` the mean and std are distance-weighted
//! (sum/min/max/median stay unweighted). The focal feature is excluded; each
//! output row copies the input attributes and adds `<field>_nbr_<stat>` columns
//! plus `nbr_count`.

use std::collections::{BTreeMap, HashMap, HashSet};

use geo::{Centroid, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

const STATS: [&str; 6] = ["mean", "median", "std", "min", "max", "sum"];

#[derive(Clone, Copy, PartialEq)]
enum Neighborhood {
    Knn,
    DistanceBand,
    Contiguity,
}

#[derive(Clone, Copy, PartialEq)]
enum Weights {
    Uniform,
    InverseDistance,
}

pub struct NeighborhoodSummaryStatisticsTool;

impl Tool for NeighborhoodSummaryStatisticsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "neighborhood_summary_statistics",
            display_name: "Neighborhood Summary Statistics",
            summary: "For each feature, summarize chosen numeric fields over its neighbours (k nearest, within a distance, or shared-edge contiguity), adding <field>_nbr_mean/median/std/min/max/sum and a neighbour count, like ArcGIS Neighborhood Summary Statistics.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (points, lines, or polygons).",
                    required: true,
                },
                ToolParamSpec {
                    name: "fields",
                    description: "Comma-separated numeric fields to summarize over neighbours.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (input copy + neighbour-statistic columns). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighborhood",
                    description: "'knn' (default), 'distance_band', or 'contiguity' (polygons sharing an edge).",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighbors",
                    description: "k for knn (default 8).",
                    required: false,
                },
                ToolParamSpec {
                    name: "distance",
                    description: "Band distance for distance_band (map units).",
                    required: false,
                },
                ToolParamSpec {
                    name: "weights",
                    description: "'uniform' (default) or 'inverse_distance' (weights the mean and std).",
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
        if args
            .get("fields")
            .and_then(Value::as_str)
            .map(str::trim)
            .unwrap_or("")
            .is_empty()
        {
            return Err(ToolError::Validation(
                "missing required string parameter 'fields'".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args.get("input").and_then(Value::as_str).unwrap();
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let n = layer.features.len();
        if n == 0 {
            return Err(ToolError::Execution("input has no features".to_string()));
        }

        // Resolve field indices.
        let field_idx: Vec<(String, usize)> = prm
            .fields
            .iter()
            .map(|f| {
                layer
                    .schema
                    .field_index(f)
                    .map(|i| (f.clone(), i))
                    .ok_or_else(|| ToolError::Validation(format!("field '{f}' not found")))
            })
            .collect::<Result<_, _>>()?;

        // Representative point per feature (None if no geometry).
        let reps: Vec<Option<(f64, f64)>> = layer
            .features
            .iter()
            .map(|f| f.geometry.as_ref().and_then(representative_point))
            .collect();

        ctx.progress.info(&format!("building {} neighbourhoods", n));

        // Neighbour lists (indices) per feature.
        let neighbors: Vec<Vec<usize>> = match prm.neighborhood {
            Neighborhood::Contiguity => contiguity_neighbors(&layer),
            Neighborhood::Knn => knn_neighbors(&reps, prm.neighbors),
            Neighborhood::DistanceBand => {
                let d = prm.distance.ok_or_else(|| {
                    ToolError::Validation(
                        "'distance' is required for neighborhood='distance_band'".to_string(),
                    )
                })?;
                distance_band_neighbors(&reps, d)
            }
        };

        // ── Build output layer: input schema + stat columns + nbr_count ──────────
        let mut out = wbvector::Layer::new("neighborhood_stats");
        if let Some(gt) = layer.geom_type {
            out = out.with_geom_type(gt);
        }
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for field in layer.schema.fields() {
            out.add_field(field.clone());
        }
        for (name, _) in &field_idx {
            for stat in STATS {
                out.add_field(FieldDef::new(
                    format!("{name}_nbr_{stat}"),
                    FieldType::Float,
                ));
            }
        }
        out.add_field(FieldDef::new("nbr_count", FieldType::Integer));

        for i in 0..n {
            let mut attrs = layer.features[i].attributes.clone();
            let nbrs = &neighbors[i];

            for (_, fidx) in &field_idx {
                // Collect neighbour values and weights.
                let mut vals: Vec<f64> = Vec::with_capacity(nbrs.len());
                let mut wts: Vec<f64> = Vec::with_capacity(nbrs.len());
                for &j in nbrs {
                    let v = layer.features[j].attributes[*fidx].as_f64();
                    let Some(v) = v else { continue };
                    if !v.is_finite() {
                        continue;
                    }
                    let w = match prm.weights {
                        Weights::Uniform => 1.0,
                        Weights::InverseDistance => inverse_distance_weight(reps[i], reps[j]),
                    };
                    vals.push(v);
                    wts.push(w);
                }
                let s = summarize(&vals, &wts);
                attrs.push(FieldValue::Float(s.mean));
                attrs.push(FieldValue::Float(s.median));
                attrs.push(FieldValue::Float(s.std));
                attrs.push(FieldValue::Float(s.min));
                attrs.push(FieldValue::Float(s.max));
                attrs.push(FieldValue::Float(s.sum));
            }
            attrs.push(FieldValue::Integer(nbrs.len() as i64));

            out.push(wbvector::Feature {
                fid: 0,
                geometry: layer.features[i].geometry.clone(),
                attributes: attrs,
            });
        }

        let total_nbrs: usize = neighbors.iter().map(Vec::len).sum();
        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert(
            "mean_neighbors".to_string(),
            json!(if n > 0 {
                total_nbrs as f64 / n as f64
            } else {
                0.0
            }),
        );
        Ok(ToolRunResult { outputs })
    }
}

// ── Neighbour construction ────────────────────────────────────────────────────

fn knn_neighbors(reps: &[Option<(f64, f64)>], k: usize) -> Vec<Vec<usize>> {
    let n = reps.len();
    let mut out = vec![Vec::new(); n];
    for i in 0..n {
        let Some(pi) = reps[i] else { continue };
        let mut ds: Vec<(f64, usize)> = (0..n)
            .filter(|&j| j != i && reps[j].is_some())
            .map(|j| {
                let pj = reps[j].unwrap();
                ((pi.0 - pj.0).powi(2) + (pi.1 - pj.1).powi(2), j)
            })
            .collect();
        ds.sort_by(|a, b| a.0.total_cmp(&b.0));
        ds.truncate(k);
        out[i] = ds.into_iter().map(|(_, j)| j).collect();
    }
    out
}

fn distance_band_neighbors(reps: &[Option<(f64, f64)>], d: f64) -> Vec<Vec<usize>> {
    let n = reps.len();
    let d2 = d * d;
    // Grid-hash for scalability.
    let cell = d.max(1e-9);
    let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    for (i, r) in reps.iter().enumerate() {
        if let Some((x, y)) = r {
            grid.entry(((x / cell).floor() as i64, (y / cell).floor() as i64))
                .or_default()
                .push(i);
        }
    }
    let mut out = vec![Vec::new(); n];
    for i in 0..n {
        let Some((x, y)) = reps[i] else { continue };
        let (gx, gy) = ((x / cell).floor() as i64, (y / cell).floor() as i64);
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(bucket) = grid.get(&(gx + dx, gy + dy)) {
                    for &j in bucket {
                        if j == i {
                            continue;
                        }
                        let (xj, yj) = reps[j].unwrap();
                        if (x - xj).powi(2) + (y - yj).powi(2) <= d2 {
                            out[i].push(j);
                        }
                    }
                }
            }
        }
    }
    out
}

type Key = (u64, u64);

fn contiguity_neighbors(layer: &wbvector::Layer) -> Vec<Vec<usize>> {
    let n = layer.features.len();
    // edge -> set of features touching it.
    let mut edge_feats: HashMap<(Key, Key), HashSet<usize>> = HashMap::new();
    for (fidx, feature) in layer.features.iter().enumerate() {
        let Some(geom) = feature.geometry.as_ref() else {
            continue;
        };
        for ring in polygon_rings(geom) {
            let m = ring.len();
            for i in 0..m {
                let a = ring[i];
                let b = ring[(i + 1) % m];
                if a == b {
                    continue;
                }
                edge_feats.entry(edge_key(a, b)).or_default().insert(fidx);
            }
        }
    }
    let mut sets: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    for feats in edge_feats.values() {
        if feats.len() < 2 {
            continue;
        }
        let list: Vec<usize> = feats.iter().copied().collect();
        for a in 0..list.len() {
            for b in 0..list.len() {
                if a != b {
                    sets[list[a]].insert(list[b]);
                }
            }
        }
    }
    sets.into_iter().map(|s| s.into_iter().collect()).collect()
}

// ── Statistics ────────────────────────────────────────────────────────────────

struct Summary {
    mean: f64,
    median: f64,
    std: f64,
    min: f64,
    max: f64,
    sum: f64,
}

fn summarize(vals: &[f64], wts: &[f64]) -> Summary {
    if vals.is_empty() {
        let nan = f64::NAN;
        return Summary {
            mean: nan,
            median: nan,
            std: nan,
            min: nan,
            max: nan,
            sum: 0.0,
        };
    }
    let sum: f64 = vals.iter().sum();
    let wsum: f64 = wts.iter().sum();
    let mean = if wsum > 0.0 {
        vals.iter().zip(wts).map(|(v, w)| v * w).sum::<f64>() / wsum
    } else {
        sum / vals.len() as f64
    };
    // Weighted population variance around the weighted mean.
    let var = if wsum > 0.0 {
        vals.iter()
            .zip(wts)
            .map(|(v, w)| w * (v - mean).powi(2))
            .sum::<f64>()
            / wsum
    } else {
        0.0
    };
    let std = var.max(0.0).sqrt();
    let min = vals.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut sorted = vals.to_vec();
    sorted.sort_by(f64::total_cmp);
    let median = if sorted.len() % 2 == 1 {
        sorted[sorted.len() / 2]
    } else {
        let m = sorted.len() / 2;
        (sorted[m - 1] + sorted[m]) / 2.0
    };
    Summary {
        mean,
        median,
        std,
        min,
        max,
        sum,
    }
}

fn inverse_distance_weight(a: Option<(f64, f64)>, b: Option<(f64, f64)>) -> f64 {
    match (a, b) {
        (Some(a), Some(b)) => {
            let d = ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt();
            if d <= 1e-9 {
                1e9
            } else {
                1.0 / d
            }
        }
        _ => 1.0,
    }
}

// ── Representative points & polygon rings ─────────────────────────────────────

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

fn key_of(x: f64, y: f64) -> Key {
    (x.to_bits(), y.to_bits())
}

fn edge_key(a: (f64, f64), b: (f64, f64)) -> (Key, Key) {
    let (ka, kb) = (key_of(a.0, a.1), key_of(b.0, b.1));
    if ka <= kb {
        (ka, kb)
    } else {
        (kb, ka)
    }
}

/// Polygon rings as (x,y) vertex chains without the closing duplicate.
fn polygon_rings(geom: &Geometry) -> Vec<Vec<(f64, f64)>> {
    let ring_pts = |ring: &Ring| -> Vec<(f64, f64)> {
        let mut pts: Vec<(f64, f64)> = Vec::new();
        for c in ring.coords() {
            let p = (c.x, c.y);
            if pts
                .last()
                .is_none_or(|l| key_of(l.0, l.1) != key_of(p.0, p.1))
            {
                pts.push(p);
            }
        }
        while pts.len() >= 2
            && key_of(pts[0].0, pts[0].1) == key_of(pts.last().unwrap().0, pts.last().unwrap().1)
        {
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

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    fields: Vec<String>,
    neighborhood: Neighborhood,
    neighbors: usize,
    distance: Option<f64>,
    weights: Weights,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let fields: Vec<String> = parse_optional_str(args, "fields")?
        .map(|s| {
            s.split(',')
                .map(str::trim)
                .filter(|f| !f.is_empty())
                .map(str::to_string)
                .collect()
        })
        .unwrap_or_default();
    if fields.is_empty() {
        return Err(ToolError::Validation(
            "'fields' must list at least one field".to_string(),
        ));
    }
    let neighborhood =
        match parse_optional_str(args, "neighborhood")?.map(|s| s.trim().to_lowercase()) {
            None => Neighborhood::Knn,
            Some(s) if s.is_empty() || s == "knn" => Neighborhood::Knn,
            Some(s) if s == "distance_band" || s == "distance" => Neighborhood::DistanceBand,
            Some(s) if s == "contiguity" => Neighborhood::Contiguity,
            Some(other) => {
                return Err(ToolError::Validation(format!(
                    "'neighborhood' must be knn|distance_band|contiguity, got '{other}'"
                )))
            }
        };
    let neighbors = match args.get("neighbors") {
        None | Some(Value::Null) => 8,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(8).max(1) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 8,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'neighbors' must be an integer".into()))?
            .max(1),
        Some(_) => return Err(ToolError::Validation("'neighbors' must be a number".into())),
    };
    let distance = match args.get("distance") {
        None | Some(Value::Null) => None,
        Some(Value::Number(n)) => n.as_f64(),
        Some(Value::String(s)) if s.trim().is_empty() => None,
        Some(Value::String(s)) => Some(
            s.trim()
                .parse::<f64>()
                .map_err(|_| ToolError::Validation("'distance' must be a number".into()))?,
        ),
        Some(_) => return Err(ToolError::Validation("'distance' must be a number".into())),
    };
    if let Some(d) = distance {
        if d.is_nan() || d <= 0.0 {
            return Err(ToolError::Validation("'distance' must be positive".into()));
        }
    }
    let weights = match parse_optional_str(args, "weights")?.map(|s| s.trim().to_lowercase()) {
        None => Weights::Uniform,
        Some(s) if s.is_empty() || s == "uniform" => Weights::Uniform,
        Some(s) if s == "inverse_distance" || s == "idw" => Weights::InverseDistance,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "'weights' must be uniform|inverse_distance, got '{other}'"
            )))
        }
    };
    Ok(Params {
        fields,
        neighborhood,
        neighbors,
        distance,
        weights,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Coord, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn point_layer(pts: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (x, y, v) in pts {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(*x, *y))),
                &[("v", (*v).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn square(x0: f64, y0: f64, s: f64, v: f64) -> (Geometry, f64) {
        (
            Geometry::polygon(
                vec![
                    Coord::xy(x0, y0),
                    Coord::xy(x0 + s, y0),
                    Coord::xy(x0 + s, y0 + s),
                    Coord::xy(x0, y0 + s),
                ],
                vec![],
            ),
            v,
        )
    }

    fn poly_layer(sqs: Vec<(Geometry, f64)>) -> String {
        let mut l = Layer::new("polys")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (g, v) in sqs {
            l.add_feature(Some(g), &[("v", v.into())]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> Layer {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = NeighborhoodSummaryStatisticsTool
            .run(&args, &ctx())
            .unwrap();
        load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    /// knn=2: the middle of three collinear points averages its two neighbours.
    #[test]
    fn knn_mean_of_two_neighbors() {
        let input = point_layer(&[(0.0, 0.0, 10.0), (1.0, 0.0, 20.0), (2.0, 0.0, 30.0)]);
        let layer =
            run(json!({ "input": input, "fields": "v", "neighborhood": "knn", "neighbors": 2 }));
        let mi = layer.schema.field_index("v_nbr_mean").unwrap();
        let ci = layer.schema.field_index("nbr_count").unwrap();
        // Feature 1 (x=1) neighbours are x=0 (10) and x=2 (30) -> mean 20.
        let f = layer.iter().nth(1).unwrap();
        assert!((f.attributes[mi].as_f64().unwrap() - 20.0).abs() < 1e-9);
        assert_eq!(f.attributes[ci].as_i64().unwrap(), 2);
    }

    /// distance_band includes exactly the points within the band.
    #[test]
    fn distance_band_counts_and_sum() {
        // x=0 has neighbours within 1.5: x=1 (v=20) only. x=10 is isolated.
        let input = point_layer(&[(0.0, 0.0, 10.0), (1.0, 0.0, 20.0), (10.0, 0.0, 99.0)]);
        let layer = run(json!({
            "input": input, "fields": "v", "neighborhood": "distance_band", "distance": 1.5
        }));
        let si = layer.schema.field_index("v_nbr_sum").unwrap();
        let ci = layer.schema.field_index("nbr_count").unwrap();
        let f0 = layer.iter().next().unwrap();
        assert_eq!(f0.attributes[ci].as_i64().unwrap(), 1);
        assert!((f0.attributes[si].as_f64().unwrap() - 20.0).abs() < 1e-9);
        let f2 = layer.iter().nth(2).unwrap();
        assert_eq!(f2.attributes[ci].as_i64().unwrap(), 0, "isolated point");
    }

    /// contiguity: a 3-in-a-row of squares — the middle touches two neighbours.
    #[test]
    fn contiguity_shared_edges() {
        let input = poly_layer(vec![
            square(0.0, 0.0, 10.0, 5.0),
            square(10.0, 0.0, 10.0, 15.0),
            square(20.0, 0.0, 10.0, 25.0),
        ]);
        let layer = run(json!({ "input": input, "fields": "v", "neighborhood": "contiguity" }));
        let mi = layer.schema.field_index("v_nbr_mean").unwrap();
        let ci = layer.schema.field_index("nbr_count").unwrap();
        // Middle square touches squares 0 (5) and 2 (25) -> mean 15, count 2.
        let mid = layer.iter().nth(1).unwrap();
        assert_eq!(mid.attributes[ci].as_i64().unwrap(), 2);
        assert!((mid.attributes[mi].as_f64().unwrap() - 15.0).abs() < 1e-9);
        // End squares touch one neighbour.
        assert_eq!(
            layer.iter().next().unwrap().attributes[ci]
                .as_i64()
                .unwrap(),
            1
        );
    }

    /// std/min/max/median columns are present and correct.
    #[test]
    fn all_stats_present() {
        let input = point_layer(&[
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 10.0),
            (2.0, 0.0, 20.0),
            (3.0, 0.0, 30.0),
        ]);
        let layer =
            run(json!({ "input": input, "fields": "v", "neighborhood": "knn", "neighbors": 3 }));
        for stat in ["mean", "median", "std", "min", "max", "sum"] {
            assert!(layer.schema.field_index(&format!("v_nbr_{stat}")).is_some());
        }
        // Feature 0 neighbours (knn=3): v=10,20,30 -> min 10, max 30, median 20, sum 60.
        let f = layer.iter().next().unwrap();
        let get = |s: &str| {
            f.attributes[layer.schema.field_index(&format!("v_nbr_{s}")).unwrap()]
                .as_f64()
                .unwrap()
        };
        assert!((get("min") - 10.0).abs() < 1e-9);
        assert!((get("max") - 30.0).abs() < 1e-9);
        assert!((get("median") - 20.0).abs() < 1e-9);
        assert!((get("sum") - 60.0).abs() < 1e-9);
    }

    #[test]
    fn rejects_missing_fields() {
        let input = point_layer(&[(0.0, 0.0, 1.0)]);
        let args: ToolArgs = serde_json::from_value(json!({ "input": input })).unwrap();
        assert!(NeighborhoodSummaryStatisticsTool.validate(&args).is_err());
    }
}
