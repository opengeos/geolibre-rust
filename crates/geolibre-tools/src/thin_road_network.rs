//! GeoLibre tool: road-network generalization (thinning) for small scales.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Thin Road Network* (Cartography):
//! reduce the density of a road network for display at smaller scales while
//! **preserving connectivity and hierarchy**. OSM road extracts (e.g. from the
//! bundled `download_osm_vector`) are far too dense to draw at small scale; the
//! suite has network *analysis* (`build_network_topology`, shortest paths) but
//! no network *generalization*.
//!
//! The tool is non-destructive: every input road keeps its geometry and gains a
//! visibility flag (`visibility_field`, 1 = kept, 0 = thinned out), exactly like
//! ArcGIS. Optionally `keep_only` filters the output to the visible roads.
//!
//! The algorithm builds the road graph (junction nodes from snapped segment
//! endpoints, one edge per road feature) and greedily hides the least important
//! short roads while never breaking connectivity:
//!
//! 1. Candidates are roads shorter than `min_length` (the small-scale
//!    invisibility threshold).
//! 2. They are considered least-important first — highest `hierarchy_field`
//!    class number (e.g. residential before motorway), then shortest.
//! 3. A candidate is hidden only if it is **not a bridge** in the current
//!    visible graph (its removal leaves its endpoints connected), so the visible
//!    network keeps exactly the connected components it started with. Important
//!    and load-bearing roads survive because they become bridges as their
//!    neighbours are thinned.
//!
//! Non-line features pass through as visible. Bridge detection is a reachability
//! check per candidate — O(E²) worst case, suited to moderate networks.
//!
//! Scope for v1: merging divided roads (dual carriageways) into a single
//! centreline — the separate ArcGIS *Merge Divided Roads* tool — is deferred;
//! this tool covers density thinning.

use std::collections::{BTreeMap, HashMap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct ThinRoadNetworkTool;

impl Tool for ThinRoadNetworkTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "thin_road_network",
            display_name: "Thin Road Network",
            summary: "Generalize a road network for small-scale display: hide short, low-hierarchy roads while preserving connectivity, flagging each road visible/thinned (non-destructive) like ArcGIS Thin Road Network.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input road (line) vector layer, format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_length",
                    description: "Roads shorter than this (CRS units) are candidates to thin out. Longer roads are always kept.",
                    required: true,
                },
                ToolParamSpec {
                    name: "hierarchy_field",
                    description: "Optional numeric road-class field; lower value = more important (kept longer). Candidates are thinned highest-class-number first.",
                    required: false,
                },
                ToolParamSpec {
                    name: "visibility_field",
                    description: "Name of the output visibility flag field (1 = kept, 0 = thinned). Default 'visible'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "keep_only",
                    description: "If true, output only the visible roads instead of flagging all of them. Default false (non-destructive).",
                    required: false,
                },
                ToolParamSpec {
                    name: "snap_tolerance",
                    description: "Distance within which road endpoints are treated as the same junction. Default 1e-6.",
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
        let schema = layer.schema.clone();

        // Build one edge per line feature.
        let mut node_id: HashMap<(i64, i64), usize> = HashMap::new();
        let mut edges: Vec<EdgeMeta> = Vec::new();
        let inv_tol = 1.0 / prm.snap_tolerance;
        let key = |x: f64, y: f64, map: &mut HashMap<(i64, i64), usize>| {
            let k = ((x * inv_tol).round() as i64, (y * inv_tol).round() as i64);
            let next = map.len();
            *map.entry(k).or_insert(next)
        };
        for (fi, feature) in layer.features.iter().enumerate() {
            let Some((a, b, len)) = feature.geometry.as_ref().and_then(line_endpoints) else {
                continue; // non-line: passes through as visible
            };
            let na = key(a.0, a.1, &mut node_id);
            let nb = key(b.0, b.1, &mut node_id);
            let hierarchy = prm
                .hierarchy_field
                .as_ref()
                .and_then(|f| feature.get(&schema, f).ok().and_then(FieldValue::as_f64))
                .unwrap_or(0.0);
            edges.push(EdgeMeta {
                fi,
                a: na,
                b: nb,
                len,
                hierarchy,
            });
        }
        let n_nodes = node_id.len();
        let n_roads = edges.len();

        // Visible adjacency: node -> list of (edge index, other node).
        let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n_nodes];
        for (ei, e) in edges.iter().enumerate() {
            adj[e.a].push((ei, e.b));
            adj[e.b].push((ei, e.a));
        }
        let mut visible = vec![true; edges.len()];

        // Candidate order: shorter than min_length, least important first
        // (highest hierarchy class number, then shortest).
        let mut candidates: Vec<usize> = (0..edges.len())
            .filter(|&ei| edges[ei].len < prm.min_length)
            .collect();
        candidates.sort_by(|&i, &j| {
            edges[j]
                .hierarchy
                .total_cmp(&edges[i].hierarchy)
                .then(edges[i].len.total_cmp(&edges[j].len))
        });

        ctx.progress.info(&format!(
            "{n_roads} road(s), {n_nodes} junction(s); {} candidate(s) under min_length",
            candidates.len()
        ));

        // Hide a candidate when doing so neither fragments the core network nor
        // erases a component's last road. `adj` holds only still-visible edges,
        // so `adj[x].len()` is x's current degree (edge `ei` still counted).
        let mut thinned = 0usize;
        for ei in candidates {
            let (a, b) = (edges[ei].a, edges[ei].b);
            let (da, db) = (adj[a].len(), adj[b].len());
            let hide = if a == b {
                true // self-loop: removing it changes no connectivity
            } else if da <= 1 && db <= 1 {
                false // the only road of its component: keep a representative
            } else if da <= 1 || db <= 1 {
                true // dead-end spur: trimming a leaf off the network
            } else {
                // Both endpoints stay in the network: hide only if not a bridge.
                still_connected(&adj, &visible, a, b, ei)
            };
            if hide {
                visible[ei] = false;
                remove_edge(&mut adj, ei, a, b);
                thinned += 1;
            }
        }

        // Map edge visibility back to features (non-line features stay visible).
        let mut feature_visible = vec![true; layer.features.len()];
        for (ei, e) in edges.iter().enumerate() {
            if !visible[ei] {
                feature_visible[e.fi] = false;
            }
        }

        let kept = feature_visible.iter().filter(|&&v| v).count();
        ctx.progress.info(&format!(
            "thinned {thinned} road(s), {kept} feature(s) kept visible"
        ));

        // Write the visibility flag, or filter to visible.
        if prm.keep_only {
            let keep: Vec<bool> = feature_visible.clone();
            let mut idx = 0;
            layer.features.retain(|_| {
                let k = keep[idx];
                idx += 1;
                k
            });
        } else {
            layer.add_field(FieldDef::new(&prm.visibility_field, FieldType::Integer));
            for (fi, feature) in layer.features.iter_mut().enumerate() {
                feature
                    .attributes
                    .push(FieldValue::Integer(feature_visible[fi] as i64));
            }
        }

        let feature_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("road_count".to_string(), json!(n_roads));
        outputs.insert("junction_count".to_string(), json!(n_nodes));
        outputs.insert("thinned_count".to_string(), json!(thinned));
        outputs.insert("visible_count".to_string(), json!(kept));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        Ok(ToolRunResult { outputs })
    }
}

struct EdgeMeta {
    fi: usize,
    a: usize,
    b: usize,
    len: f64,
    hierarchy: f64,
}

/// Removes edge `ei` from both endpoints' adjacency lists.
fn remove_edge(adj: &mut [Vec<(usize, usize)>], ei: usize, a: usize, b: usize) {
    adj[a].retain(|&(e, _)| e != ei);
    adj[b].retain(|&(e, _)| e != ei);
}

/// True if `a` can still reach `b` over the visible graph without edge `skip`.
fn still_connected(
    adj: &[Vec<(usize, usize)>],
    visible: &[bool],
    a: usize,
    b: usize,
    skip: usize,
) -> bool {
    let mut seen = vec![false; adj.len()];
    let mut stack = vec![a];
    seen[a] = true;
    while let Some(u) = stack.pop() {
        if u == b {
            return true;
        }
        for &(ei, v) in &adj[u] {
            if ei == skip || !visible[ei] || seen[v] {
                continue;
            }
            seen[v] = true;
            stack.push(v);
        }
    }
    false
}

/// (start point, end point, total length) of a line geometry.
type Endpoints = ((f64, f64), (f64, f64), f64);

/// Overall endpoints (snapped) and total length of a line geometry.
fn line_endpoints(geom: &Geometry) -> Option<Endpoints> {
    let parts: Vec<&Vec<Coord>> = match geom {
        Geometry::LineString(cs) => vec![cs],
        Geometry::MultiLineString(ls) => ls.iter().collect(),
        _ => return None,
    };
    let mut first: Option<(f64, f64)> = None;
    let mut last: (f64, f64) = (0.0, 0.0);
    let mut len = 0.0;
    for coords in parts {
        if coords.len() < 2 {
            continue;
        }
        if first.is_none() {
            first = Some((coords[0].x, coords[0].y));
        }
        last = {
            let c = coords.last().unwrap();
            (c.x, c.y)
        };
        for w in coords.windows(2) {
            len += ((w[0].x - w[1].x).powi(2) + (w[0].y - w[1].y).powi(2)).sqrt();
        }
    }
    first.map(|f| (f, last, len))
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    min_length: f64,
    hierarchy_field: Option<String>,
    visibility_field: String,
    keep_only: bool,
    snap_tolerance: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let min_length = parse_optional_f64(args, "min_length")?.ok_or_else(|| {
        ToolError::Validation("missing required parameter 'min_length'".to_string())
    })?;
    if !(min_length > 0.0 && min_length.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'min_length' must be a positive number".to_string(),
        ));
    }
    let hierarchy_field = parse_optional_str(args, "hierarchy_field")?.map(str::to_string);
    let visibility_field = parse_optional_str(args, "visibility_field")?
        .map(str::to_string)
        .unwrap_or_else(|| "visible".to_string());
    let keep_only = parse_optional_bool(args, "keep_only")?.unwrap_or(false);
    let snap_tolerance = parse_optional_f64(args, "snap_tolerance")?.unwrap_or(1e-6);
    if !(snap_tolerance > 0.0 && snap_tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'snap_tolerance' must be a positive number".to_string(),
        ));
    }
    Ok(Params {
        min_length,
        hierarchy_field,
        visibility_field,
        keep_only,
        snap_tolerance,
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

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "true" | "yes" | "1" => Ok(Some(true)),
            "false" | "no" | "0" => Ok(Some(false)),
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
    use wbvector::{memory_store, FieldDef, FieldType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn line(a: (f64, f64), b: (f64, f64)) -> Geometry {
        Geometry::line_string(vec![Coord::xy(a.0, a.1), Coord::xy(b.0, b.1)])
    }

    fn run(build: impl FnOnce(&mut Layer), args: serde_json::Value) -> (ToolRunResult, Layer) {
        let mut layer = Layer::new("roads");
        build(&mut layer);
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);
        let mut v = args;
        v["input"] = json!(input);
        let args: ToolArgs = serde_json::from_value(v).unwrap();
        let out = ThinRoadNetworkTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn vis(layer: &Layer, i: usize) -> i64 {
        match layer.features[i].get(&layer.schema, "visible").unwrap() {
            FieldValue::Integer(v) => *v,
            other => panic!("visible should be integer, got {other:?}"),
        }
    }

    #[test]
    fn keeps_connectivity_thinning_a_redundant_parallel_road() {
        // A--B connected by two parallel roads: a long one (len 10) and a short
        // detour via C (two 3-unit segments). With min_length 5 the short
        // segments are candidates; one can be hidden but connectivity via the
        // long road stays, so the network never splits.
        let (out, layer) = run(
            |l| {
                l.add_feature(Some(line((0.0, 0.0), (10.0, 0.0))), &[])
                    .unwrap(); // long A-B
                l.add_feature(Some(line((0.0, 0.0), (5.0, 3.0))), &[])
                    .unwrap(); // A-C (len ~5.83)
                l.add_feature(Some(line((5.0, 3.0), (10.0, 0.0))), &[])
                    .unwrap(); // C-B
            },
            json!({ "min_length": 7.0 }),
        );
        // The two short detour roads are candidates; at least one is thinned.
        assert!(out.outputs["thinned_count"].as_u64().unwrap() >= 1);
        // The long A-B road is longer than min_length -> always visible.
        assert_eq!(vis(&layer, 0), 1);
    }

    #[test]
    fn keeps_the_last_road_of_a_component() {
        // A lone short road (its own component, two dead-ends) is a representative
        // and must never be erased, however short.
        let (out, layer) = run(
            |l| {
                l.add_feature(Some(line((0.0, 0.0), (1.0, 0.0))), &[])
                    .unwrap();
            },
            json!({ "min_length": 100.0 }),
        );
        assert_eq!(out.outputs["thinned_count"], json!(0));
        assert_eq!(vis(&layer, 0), 1);
    }

    #[test]
    fn prunes_a_short_dead_end_spur() {
        // A long backbone A-B (kept) with a short dead-end spur B-C. The spur's
        // free end C is a leaf, so trimming it does not fragment the network.
        let (out, layer) = run(
            |l| {
                l.add_feature(Some(line((0.0, 0.0), (100.0, 0.0))), &[])
                    .unwrap(); // A-B backbone
                l.add_feature(Some(line((100.0, 0.0), (110.0, 0.0))), &[])
                    .unwrap(); // B-C spur (len 10)
            },
            json!({ "min_length": 50.0 }),
        );
        assert_eq!(out.outputs["thinned_count"], json!(1));
        assert_eq!(vis(&layer, 0), 1, "backbone kept");
        assert_eq!(vis(&layer, 1), 0, "dead-end spur pruned");
    }

    #[test]
    fn hierarchy_thins_low_class_before_high() {
        // A and B joined by a motorway (class 1) and a residential (class 3),
        // both short and parallel. The residential is thinned first; the
        // motorway then becomes a bridge and is kept.
        let (_, layer) = run(
            |l| {
                l.add_field(FieldDef::new("class", FieldType::Integer));
                l.add_feature(
                    Some(line((0.0, 0.0), (4.0, 0.0))),
                    &[("class", FieldValue::Integer(1))],
                )
                .unwrap(); // motorway
                l.add_feature(
                    Some(line((0.0, 0.0), (4.0, 0.0))),
                    &[("class", FieldValue::Integer(3))],
                )
                .unwrap(); // residential
            },
            json!({ "min_length": 10.0, "hierarchy_field": "class" }),
        );
        assert_eq!(vis(&layer, 0), 1, "motorway should be kept");
        assert_eq!(vis(&layer, 1), 0, "residential should be thinned");
    }

    #[test]
    fn keep_only_filters_output() {
        let (out, layer) = run(
            |l| {
                l.add_feature(Some(line((0.0, 0.0), (10.0, 0.0))), &[])
                    .unwrap();
                l.add_feature(Some(line((0.0, 0.0), (5.0, 3.0))), &[])
                    .unwrap();
                l.add_feature(Some(line((5.0, 3.0), (10.0, 0.0))), &[])
                    .unwrap();
            },
            json!({ "min_length": 7.0, "keep_only": true }),
        );
        let visible = out.outputs["visible_count"].as_u64().unwrap() as usize;
        assert_eq!(layer.len(), visible, "keep_only should drop thinned roads");
        assert!(layer.len() < 3);
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = ThinRoadNetworkTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "input": "r.geojson" })).is_err(),
            "missing min_length"
        );
        assert!(bad(json!({ "input": "r.geojson", "min_length": 0 })).is_err());
        assert!(bad(json!({ "input": "r.geojson", "min_length": -5 })).is_err());
        assert!(bad(json!({ "input": "r.geojson", "min_length": 100 })).is_ok());
    }
}
