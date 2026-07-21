//! GeoLibre tool: polygon adjacency (contiguity) table with shared borders.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Polygon Neighbors* (Analysis). The
//! bundled suite tests topology rules but never exports adjacency as data.
//! Contiguity tables feed spatial-weights construction, redistricting QA, and
//! region-merge heuristics — and the shared-edge machinery already exists inside
//! `eliminate_polygons` / `simplify_shared_edges`; this exposes it as a table.
//!
//! Every polygon boundary is decomposed into undirected edges (endpoint-keyed,
//! optionally snapped to a grid so near-coincident borders match). For each edge
//! shared by two features their pair accumulates that edge's length; every
//! vertex shared by two features contributes to their node count. The result is
//! one row per neighbouring pair — `src_id, nbr_id, length, node_count` — where
//! `length > 0` marks edge (rook) neighbours and `length == 0, node_count > 0`
//! marks node-only (corner) neighbours. `both_sides` emits each pair once or in
//! both directions (ArcGIS-compatible). Output is a geometry-less attribute
//! table (or a CSV when the path ends in `.csv`).

use std::collections::{BTreeMap, HashMap, HashSet};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Feature, FieldDef, FieldType, FieldValue, Geometry, Layer, Ring};

use crate::common::write_text_output;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct PolygonNeighborsTool;

impl Tool for PolygonNeighborsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "polygon_neighbors",
            display_name: "Polygon Neighbors",
            summary: "Build a polygon adjacency (contiguity) table: neighbouring pairs with the length of their shared border and a count of point-touches, like ArcGIS Polygon Neighbors.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon vector layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output table path — a CSV (extension .csv) or a geometry-less vector table. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "id_field",
                    description: "Field identifying each polygon in the table. Default: the feature index.",
                    required: false,
                },
                ToolParamSpec {
                    name: "both_sides",
                    description: "Emit each neighbour pair twice (src→nbr and nbr→src). Default false (once per pair).",
                    required: false,
                },
                ToolParamSpec {
                    name: "snap_tolerance",
                    description: "Quantize vertices onto a grid of this size (CRS units) before matching borders. Default 0 (exact coordinates).",
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
        let id_idx = match &prm.id_field {
            Some(f) => Some(
                layer
                    .schema
                    .field_index(f)
                    .ok_or_else(|| ToolError::Validation(format!("id_field '{f}' not found")))?,
            ),
            None => None,
        };

        // ── Build edge → features and vertex → features maps ──────────────────
        let mut edge_feats: HashMap<(Key, Key), HashSet<usize>> = HashMap::new();
        let mut edge_len: HashMap<(Key, Key), f64> = HashMap::new();
        let mut vert_feats: HashMap<Key, HashSet<usize>> = HashMap::new();
        let mut ids: Vec<String> = Vec::new();
        let mut poly_count = 0usize;

        for (fidx, feature) in layer.features.iter().enumerate() {
            let id = match id_idx {
                Some(i) => field_key(&feature.attributes[i]),
                None => fidx.to_string(),
            };
            ids.push(id);
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let rings = polygon_rings(geom, prm.snap_tolerance);
            if rings.is_empty() {
                continue;
            }
            poly_count += 1;
            for ring in &rings {
                let n = ring.len();
                for i in 0..n {
                    let a = ring[i];
                    let b = ring[(i + 1) % n];
                    vert_feats.entry(key(a)).or_default().insert(fidx);
                    if key(a) == key(b) {
                        continue;
                    }
                    let e = edge_key(a, b);
                    edge_feats.entry(e).or_default().insert(fidx);
                    edge_len.entry(e).or_insert_with(|| dist(a, b));
                }
            }
        }

        // ── Accumulate per-pair shared length and shared nodes ────────────────
        let mut pairs: BTreeMap<(usize, usize), (f64, HashSet<Key>)> = BTreeMap::new();
        for (e, feats) in &edge_feats {
            if feats.len() < 2 {
                continue;
            }
            let len = edge_len[e];
            let list: Vec<usize> = feats.iter().copied().collect();
            for i in 0..list.len() {
                for j in (i + 1)..list.len() {
                    let (a, b) = order(list[i], list[j]);
                    pairs.entry((a, b)).or_default().0 += len;
                }
            }
        }
        for (v, feats) in &vert_feats {
            if feats.len() < 2 {
                continue;
            }
            let list: Vec<usize> = feats.iter().copied().collect();
            for i in 0..list.len() {
                for j in (i + 1)..list.len() {
                    let (a, b) = order(list[i], list[j]);
                    pairs.entry((a, b)).or_default().1.insert(*v);
                }
            }
        }

        ctx.progress.info(&format!(
            "{poly_count} polygon(s); {} neighbour pair(s)",
            pairs.len()
        ));

        // ── Emit the adjacency table ──────────────────────────────────────────
        let mut table = Layer::new("polygon_neighbors");
        table.add_field(FieldDef::new("src_id", FieldType::Text));
        table.add_field(FieldDef::new("nbr_id", FieldType::Text));
        table.add_field(FieldDef::new("length", FieldType::Float));
        table.add_field(FieldDef::new("node_count", FieldType::Integer));

        let mut csv = String::from("src_id,nbr_id,length,node_count\n");
        let mut rows = 0usize;
        let mut edge_neighbors = 0usize;
        let push_row =
            |src: &str, nbr: &str, len: f64, nodes: i64, table: &mut Layer, csv: &mut String| {
                table.push(Feature {
                    fid: 0,
                    geometry: None,
                    attributes: vec![
                        FieldValue::Text(src.to_string()),
                        FieldValue::Text(nbr.to_string()),
                        FieldValue::Float(len),
                        FieldValue::Integer(nodes),
                    ],
                });
                csv.push_str(&format!("{src},{nbr},{len},{nodes}\n"));
            };

        for ((a, b), (len, nodeset)) in &pairs {
            let nodes = nodeset.len() as i64;
            if *len > 0.0 {
                edge_neighbors += 1;
            }
            push_row(&ids[*a], &ids[*b], *len, nodes, &mut table, &mut csv);
            rows += 1;
            if prm.both_sides {
                push_row(&ids[*b], &ids[*a], *len, nodes, &mut table, &mut csv);
                rows += 1;
            }
        }

        // ── Write (CSV if the path ends in .csv, else a table layer) ──────────
        let out_path = match output {
            Some(path) if path.to_ascii_lowercase().ends_with(".csv") => {
                write_text_output(&csv, path)?;
                path.to_string()
            }
            other => write_or_store_layer(table, other)?,
        };

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("polygon_count".to_string(), json!(poly_count));
        outputs.insert("pair_count".to_string(), json!(pairs.len()));
        outputs.insert("edge_neighbor_pairs".to_string(), json!(edge_neighbors));
        outputs.insert("row_count".to_string(), json!(rows));
        Ok(ToolRunResult { outputs })
    }
}

// ── Keys and geometry ────────────────────────────────────────────────────────

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

fn order(a: usize, b: usize) -> (usize, usize) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

fn dist(a: P, b: P) -> f64 {
    (a.x - b.x).hypot(a.y - b.y)
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

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    id_field: Option<String>,
    both_sides: bool,
    snap_tolerance: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let id_field = parse_optional_str(args, "id_field")?.map(str::to_string);
    let both_sides = parse_optional_bool(args, "both_sides")?.unwrap_or(false);
    let snap_tolerance = match args.get("snap_tolerance") {
        None | Some(Value::Null) => 0.0,
        Some(Value::Number(n)) => n.as_f64().unwrap_or(0.0),
        Some(Value::String(s)) if s.trim().is_empty() => 0.0,
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'snap_tolerance' must be a number".into()))?,
        Some(_) => {
            return Err(ToolError::Validation(
                "'snap_tolerance' must be a number".into(),
            ))
        }
    };
    if !(snap_tolerance >= 0.0 && snap_tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "'snap_tolerance' must be a non-negative number".to_string(),
        ));
    }
    Ok(Params {
        id_field,
        both_sides,
        snap_tolerance,
    })
}

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
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
    use wbvector::{memory_store, Coord, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
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

    fn layer_of(named: &[(&str, Geometry)]) -> String {
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
        let out = PolygonNeighborsTool.run(&args, &ctx()).unwrap();
        let table = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, table)
    }

    fn rows(table: &Layer) -> Vec<(String, String, f64, i64)> {
        let si = table.schema.field_index("src_id").unwrap();
        let ni = table.schema.field_index("nbr_id").unwrap();
        let li = table.schema.field_index("length").unwrap();
        let ci = table.schema.field_index("node_count").unwrap();
        table
            .iter()
            .map(|f| {
                (
                    f.attributes[si].as_str().unwrap().to_string(),
                    f.attributes[ni].as_str().unwrap().to_string(),
                    f.attributes[li].as_f64().unwrap(),
                    f.attributes[ci].as_i64().unwrap(),
                )
            })
            .collect()
    }

    /// Two squares sharing a full edge -> one pair with that border length.
    #[test]
    fn edge_neighbors_share_border_length() {
        // A [0,10]^2, B [10,20]x[0,10]: shared vertical edge x=10, length 10.
        let input = layer_of(&[
            ("A", square(0.0, 0.0, 10.0)),
            ("B", square(10.0, 0.0, 10.0)),
        ]);
        let (out, table) = run(json!({ "input": input, "id_field": "name" }));
        assert_eq!(out.outputs["pair_count"], json!(1));
        let r = rows(&table);
        assert_eq!(r.len(), 1);
        assert_eq!((r[0].0.as_str(), r[0].1.as_str()), ("A", "B"));
        assert!(
            (r[0].2 - 10.0).abs() < 1e-9,
            "shared length {} != 10",
            r[0].2
        );
    }

    /// both_sides emits each pair in both directions.
    #[test]
    fn both_sides_doubles_rows() {
        let input = layer_of(&[
            ("A", square(0.0, 0.0, 10.0)),
            ("B", square(10.0, 0.0, 10.0)),
        ]);
        let (out, table) = run(json!({ "input": input, "id_field": "name", "both_sides": true }));
        assert_eq!(out.outputs["row_count"], json!(2));
        let mut pairs: Vec<(String, String)> =
            rows(&table).into_iter().map(|r| (r.0, r.1)).collect();
        pairs.sort();
        assert_eq!(
            pairs,
            vec![("A".into(), "B".into()), ("B".into(), "A".into())]
        );
    }

    /// Squares touching only at a corner are node-only neighbours (length 0).
    #[test]
    fn corner_touch_is_node_neighbor() {
        // A [0,10]^2 and C [10,20]x[10,20] touch only at (10,10).
        let input = layer_of(&[
            ("A", square(0.0, 0.0, 10.0)),
            ("C", square(10.0, 10.0, 10.0)),
        ]);
        let (_o, table) = run(json!({ "input": input, "id_field": "name" }));
        let r = rows(&table);
        assert_eq!(r.len(), 1, "one node-only pair");
        assert_eq!(r[0].2, 0.0, "corner touch has zero shared length");
        assert!(r[0].3 >= 1, "node_count should be >= 1");
    }

    /// Disjoint squares have no neighbours.
    #[test]
    fn disjoint_squares_have_no_neighbors() {
        let input = layer_of(&[
            ("A", square(0.0, 0.0, 5.0)),
            ("B", square(100.0, 100.0, 5.0)),
        ]);
        let (out, _t) = run(json!({ "input": input, "id_field": "name" }));
        assert_eq!(out.outputs["pair_count"], json!(0));
    }

    #[test]
    fn rejects_missing_input() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(PolygonNeighborsTool.validate(&args).is_err());
    }
}
