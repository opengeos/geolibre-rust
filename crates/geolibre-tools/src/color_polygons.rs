//! GeoLibre tool: assign an adjacency-safe color index to polygons.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Calculate Color Theorem Field*
//! (Cartography). Given a polygon coverage, it assigns each polygon a small
//! integer so that no two adjacent polygons share a value — instant
//! choropleth-safe styling for the repo's `render_vector_png` / PMTiles outputs
//! (parcels, admin units, `build_balanced_zones` output). Nothing in the bundled
//! whitebox suite assigns map colors; the shared-edge adjacency it needs already
//! exists inside `polygon_neighbors`.
//!
//! Adjacency is built by decomposing every polygon boundary into undirected
//! edges (optionally snapped to a grid): two features sharing an edge are
//! `edge` (rook) neighbours; with `adjacency=edge_or_corner` a shared vertex
//! also makes them neighbours (queen). Coloring uses the DSATUR heuristic
//! (colour the most-saturated uncoloured polygon next), which stays within 4–6
//! colours on planar subdivision maps. The colour index (1-based) is written to
//! the `field` attribute.

use std::collections::{BTreeSet, HashMap, HashSet};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct ColorPolygonsTool;

impl Tool for ColorPolygonsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "color_polygons",
            display_name: "Color Polygons",
            summary: "Assign each polygon a small integer colour index so no two adjacent polygons share a value (like ArcGIS Calculate Color Theorem Field) — adjacency-safe choropleth styling via DSATUR graph colouring over shared-edge (or shared-corner) contiguity.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon layer with the colour field added. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "field",
                    description: "Name of the colour index field to write (default 'color_id').",
                    required: false,
                },
                ToolParamSpec {
                    name: "adjacency",
                    description: "'edge' (rook, shared border; default) or 'edge_or_corner' (queen, shared border or vertex).",
                    required: false,
                },
                ToolParamSpec {
                    name: "snap_tolerance",
                    description: "Grid size to snap vertices to before matching borders (default 0 = exact).",
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

        let mut layer = load_input_layer(input)?;
        let n = layer.features.len();

        // Build adjacency: edges (and optionally vertices) shared by two features.
        let mut edge_feats: HashMap<(Key, Key), Vec<usize>> = HashMap::new();
        let mut vert_feats: HashMap<Key, Vec<usize>> = HashMap::new();
        for (fidx, feat) in layer.features.iter().enumerate() {
            let Some(geom) = feat.geometry.as_ref() else {
                continue;
            };
            for ring in polygon_rings(geom, prm.snap_tolerance) {
                let m = ring.len();
                for i in 0..m {
                    let a = ring[i];
                    let b = ring[(i + 1) % m];
                    push_unique(edge_feats.entry(edge_key(a, b)).or_default(), fidx);
                    if prm.queen {
                        push_unique(vert_feats.entry(key(a)).or_default(), fidx);
                    }
                }
            }
        }

        let mut adj: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];
        for feats in edge_feats.values() {
            add_clique(&mut adj, feats);
        }
        if prm.queen {
            for feats in vert_feats.values() {
                add_clique(&mut adj, feats);
            }
        }

        ctx.progress.info(&format!(
            "colouring {n} polygon(s) over the adjacency graph"
        ));

        let colors = dsatur(&adj);
        let num_colors = colors.iter().copied().max().map(|c| c + 1).unwrap_or(0);

        // Verify the coloring is proper (defensive: report conflicts if any).
        let mut conflicts = 0usize;
        for a in 0..n {
            for &b in &adj[a] {
                if b > a && colors[a] == colors[b] {
                    conflicts += 1;
                }
            }
        }

        layer.add_field(FieldDef::new(prm.field.clone(), FieldType::Integer));
        for (i, feat) in layer.features.iter_mut().enumerate() {
            feat.attributes
                .push(FieldValue::Integer(colors[i] as i64 + 1));
        }

        let out_path = write_or_store_layer(layer, output)?;
        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("num_colors".to_string(), json!(num_colors));
        outputs.insert("conflicts".to_string(), json!(conflicts));
        Ok(ToolRunResult { outputs })
    }
}

/// DSATUR greedy graph coloring. Returns a 0-based colour per vertex.
fn dsatur(adj: &[BTreeSet<usize>]) -> Vec<usize> {
    let n = adj.len();
    let mut color = vec![usize::MAX; n];
    let degree: Vec<usize> = adj.iter().map(|a| a.len()).collect();
    // Saturation = number of distinct colours among coloured neighbours.
    let mut sat_colors: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    for _ in 0..n {
        // Pick the uncoloured vertex with max saturation, then max degree, then
        // smallest index (deterministic).
        let mut best: Option<usize> = None;
        for v in 0..n {
            if color[v] != usize::MAX {
                continue;
            }
            best = Some(match best {
                None => v,
                Some(b) => {
                    let key_v = (sat_colors[v].len(), degree[v]);
                    let key_b = (sat_colors[b].len(), degree[b]);
                    if key_v > key_b {
                        v
                    } else {
                        b
                    }
                }
            });
        }
        let Some(v) = best else { break };
        // Smallest colour not used by a neighbour.
        let used: HashSet<usize> = adj[v]
            .iter()
            .filter_map(|&u| (color[u] != usize::MAX).then_some(color[u]))
            .collect();
        let mut c = 0;
        while used.contains(&c) {
            c += 1;
        }
        color[v] = c;
        for &u in &adj[v] {
            sat_colors[u].insert(c);
        }
    }
    color
}

fn add_clique(adj: &mut [BTreeSet<usize>], feats: &[usize]) {
    for i in 0..feats.len() {
        for j in (i + 1)..feats.len() {
            let (a, b) = (feats[i], feats[j]);
            if a != b {
                adj[a].insert(b);
                adj[b].insert(a);
            }
        }
    }
}

fn push_unique(v: &mut Vec<usize>, x: usize) {
    if !v.contains(&x) {
        v.push(x);
    }
}

// ── Geometry helpers (shared-edge adjacency, mirrors polygon_neighbors) ─────────

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

/// All rings of a polygon geometry as canonical vertex chains (no closing dup).
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

// ── Parameters ───────────────────────────────────────────────────────────────

struct Params {
    field: String,
    queen: bool,
    snap_tolerance: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let field = parse_optional_str(args, "field")?
        .unwrap_or("color_id")
        .to_string();
    let queen = match args.get("adjacency").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("edge") => false,
        Some("edge_or_corner") => true,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'adjacency' must be 'edge' or 'edge_or_corner', got '{o}'"
            )))
        }
    };
    let snap_tolerance = match args.get("snap_tolerance") {
        None | Some(Value::Null) => 0.0,
        Some(Value::Number(n)) => n.as_f64().unwrap_or(0.0),
        Some(Value::String(s)) if s.trim().is_empty() => 0.0,
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'snap_tolerance' must be a number".into()))?,
        _ => {
            return Err(ToolError::Validation(
                "'snap_tolerance' must be a number".into(),
            ))
        }
    };
    if snap_tolerance < 0.0 {
        return Err(ToolError::Validation(
            "'snap_tolerance' must be >= 0".to_string(),
        ));
    }
    Ok(Params {
        field,
        queen,
        snap_tolerance,
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

    fn square(x0: f64, y0: f64, s: f64) -> Geometry {
        let ring = Ring::new(vec![
            Coord::xy(x0, y0),
            Coord::xy(x0 + s, y0),
            Coord::xy(x0 + s, y0 + s),
            Coord::xy(x0, y0 + s),
            Coord::xy(x0, y0),
        ]);
        Geometry::Polygon {
            exterior: ring,
            interiors: vec![],
        }
    }

    fn grid_layer(cols: usize, rows: usize) -> String {
        let mut l = Layer::new("grid")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("gid", FieldType::Integer));
        let mut gid = 0i64;
        for r in 0..rows {
            for c in 0..cols {
                l.add_feature(
                    Some(square(c as f64, r as f64, 1.0)),
                    &[("gid", gid.into())],
                )
                .unwrap();
                gid += 1;
            }
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ColorPolygonsTool.run(&args, &ctx()).unwrap();
        let l = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, l)
    }

    /// A rook-adjacency grid is 2-colourable (checkerboard); no adjacent pair
    /// shares a colour.
    #[test]
    fn grid_is_two_colorable_rook() {
        let (out, l) = run(json!({ "input": grid_layer(5, 5), "adjacency": "edge" }));
        assert_eq!(
            out.outputs["conflicts"],
            json!(0),
            "coloring must be proper"
        );
        assert_eq!(
            out.outputs["num_colors"],
            json!(2),
            "rook grid needs 2 colors"
        );
        let cf = l.schema.field_index("color_id").unwrap();
        // Neighbours (i, i+1) horizontally must differ.
        let colors: Vec<i64> = l
            .features
            .iter()
            .map(|f| f.attributes[cf].as_i64().unwrap())
            .collect();
        for r in 0..5 {
            for c in 0..4 {
                assert_ne!(colors[r * 5 + c], colors[r * 5 + c + 1]);
            }
        }
    }

    /// Queen adjacency (diagonal touches count) needs more than 2 colours on a
    /// grid, and is still a proper colouring.
    #[test]
    fn grid_queen_is_proper() {
        let (out, _l) = run(json!({ "input": grid_layer(4, 4), "adjacency": "edge_or_corner" }));
        assert_eq!(out.outputs["conflicts"], json!(0));
        assert!(
            out.outputs["num_colors"].as_i64().unwrap() >= 4,
            "queen grid needs >= 4 colors"
        );
    }

    /// Custom field name is honoured.
    #[test]
    fn custom_field_name() {
        let (_out, l) = run(json!({ "input": grid_layer(3, 3), "field": "hue" }));
        assert!(l.schema.field_index("hue").is_some());
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ColorPolygonsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "adjacency": "touch" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "snap_tolerance": -1 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "adjacency": "edge_or_corner" })).is_ok());
    }
}
