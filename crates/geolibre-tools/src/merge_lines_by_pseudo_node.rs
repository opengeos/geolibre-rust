//! GeoLibre tool: dissolve line chains that meet at pseudo-nodes.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Merge Lines By Pseudo Node*
//! (Topographic Production; equivalent to *Unsplit Line* in Data Management /
//! Analysis): collapse runs of line segments joined end-to-end at *pseudo-nodes*
//! — nodes where exactly two line ends meet (degree 2) — into single continuous
//! polylines, while preserving true intersections (degree ≥ 3). The standard
//! topology cleanup after raster stream/road vectorization or tiled imports,
//! where a single real road arrives chopped into many small segments.
//!
//! The bundled `merge_line_segments` stitches coincident endpoints purely
//! geometrically; it does not stop at real intersections or honour attribute
//! breaks. This tool is pseudo-node-aware and attribute-respecting:
//!
//! - Chains are only extended *through* degree-2 nodes, so a T- or X-junction
//!   (degree ≥ 3) always ends a chain and stays a node.
//! - With `dissolve_fields`, two segments merge only when they agree on every
//!   listed field, so a highway does not fuse with the residential street it
//!   happens to touch end-to-end.
//!
//! Each merged feature keeps the attributes of its first (lowest original id)
//! segment and gains a `merged_count` of how many segments were fused.
//! Non-line features pass through untouched; `MultiLineString` inputs are
//! exploded into their parts before merging.

use std::collections::{BTreeMap, HashMap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct MergeLinesByPseudoNodeTool;

impl Tool for MergeLinesByPseudoNodeTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "merge_lines_by_pseudo_node",
            display_name: "Merge Lines By Pseudo Node",
            summary: "Dissolve chains of line segments that meet at degree-2 pseudo-nodes into single polylines, preserving true intersections and (optionally) attribute breaks.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input line vector file path, format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "dissolve_fields",
                    description: "Comma-separated attribute field names; two segments merge only where they agree on all of them. If omitted, segments merge across every degree-2 node regardless of attributes.",
                    required: false,
                },
                ToolParamSpec {
                    name: "snap_tolerance",
                    description: "Endpoints within this distance (CRS units) are treated as the same node. Default 1e-6.",
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
        let schema = layer.schema.clone();

        // Resolve dissolve-field indices up front.
        let mut dissolve_idx = Vec::new();
        for name in &prm.dissolve_fields {
            let idx = schema.field_index(name).ok_or_else(|| {
                ToolError::Validation(format!("dissolve field '{name}' not found in input schema"))
            })?;
            dissolve_idx.push(idx);
        }

        // Explode features into individual polylines, remembering their source
        // feature (for attributes). Non-line features are collected for
        // pass-through.
        let mut lines: Vec<LineRef> = Vec::new();
        let mut passthrough: Vec<Feature> = Vec::new();
        for (fidx, feature) in layer.features.iter().enumerate() {
            match &feature.geometry {
                Some(Geometry::LineString(cs)) => {
                    if cs.len() >= 2 {
                        lines.push(LineRef {
                            coords: cs.clone(),
                            feat: fidx,
                        });
                    }
                }
                Some(Geometry::MultiLineString(parts)) => {
                    for cs in parts {
                        if cs.len() >= 2 {
                            lines.push(LineRef {
                                coords: cs.clone(),
                                feat: fidx,
                            });
                        }
                    }
                }
                _ => passthrough.push(feature.clone()),
            }
        }
        let input_line_count = lines.len();

        // Assign a node id to every distinct endpoint (snapped) and record which
        // line-ends touch it.
        let mut nodes = NodeIndex::new(prm.snap_tolerance);
        let mut ends: Vec<(usize, usize)> = Vec::with_capacity(lines.len()); // (start_node, end_node)
        for l in &lines {
            let a = nodes.id_of(&l.coords[0]);
            let b = nodes.id_of(l.coords.last().unwrap());
            ends.push((a, b));
        }
        // Node -> incident line-ends: (line_index, which_end 0=start 1=end).
        let mut incidence: HashMap<usize, Vec<(usize, u8)>> = HashMap::new();
        for (li, &(a, b)) in ends.iter().enumerate() {
            incidence.entry(a).or_default().push((li, 0));
            incidence.entry(b).or_default().push((li, 1));
        }

        ctx.progress.info(&format!(
            "{input_line_count} line segment(s), {} node(s)",
            incidence.len()
        ));

        // Walk chains: start each unused line, extend through mergeable
        // pseudo-nodes on both ends, and emit the concatenated polyline.
        let dissolve_key = |li: usize| -> Vec<String> {
            let f = &layer.features[lines[li].feat];
            dissolve_idx
                .iter()
                .map(|&i| {
                    f.attributes
                        .get(i)
                        .map(field_value_string)
                        .unwrap_or_default()
                })
                .collect()
        };

        let mut used = vec![false; lines.len()];
        let mut out_lines: Vec<(Vec<Coord>, usize, usize)> = Vec::new(); // (coords, source_feat, merged_count)
        for seed in 0..lines.len() {
            if used[seed] {
                continue;
            }
            used[seed] = true;
            let mut coords = lines[seed].coords.clone();
            let (mut head, mut tail) = ends[seed];
            let (mut head_line, mut tail_line) = (seed, seed);
            let seed_key = dissolve_key(seed);
            let mut members = vec![seed];

            // Extend forward from the tail node (the current tail segment is
            // `tail_line`, so `next_link` skips it and finds the continuation).
            loop {
                match next_link(tail, tail_line, &incidence, &used, &dissolve_key, &seed_key) {
                    Some((nl, other_end)) => {
                        used[nl] = true;
                        members.push(nl);
                        // Append nl's coords oriented to continue from `tail`.
                        let (na, nb) = ends[nl];
                        if other_end == 0 {
                            // nl starts at `tail`: append its tail-ward coords.
                            append_skip_first(&mut coords, &lines[nl].coords, false);
                            tail = nb;
                        } else {
                            append_skip_first(&mut coords, &lines[nl].coords, true);
                            tail = na;
                        }
                        tail_line = nl;
                        if tail == head {
                            break; // closed loop
                        }
                    }
                    None => break,
                }
            }
            // Extend backward from the head node (prepend).
            loop {
                if tail == head {
                    break;
                }
                match next_link(head, head_line, &incidence, &used, &dissolve_key, &seed_key) {
                    Some((nl, other_end)) => {
                        used[nl] = true;
                        members.push(nl);
                        let (na, nb) = ends[nl];
                        // Prepend nl so its far end becomes the new head.
                        if other_end == 1 {
                            // nl ends at `head`: prepend its coords as-is.
                            prepend_skip_last(&mut coords, &lines[nl].coords, false);
                            head = na;
                        } else {
                            prepend_skip_last(&mut coords, &lines[nl].coords, true);
                            head = nb;
                        }
                        head_line = nl;
                    }
                    None => break,
                }
            }

            out_lines.push((coords, lines[seed].feat, members.len()));
        }

        let merged_chains = out_lines.iter().filter(|(_, _, n)| *n > 1).count();

        // Build the output layer: merged lines first (with merged_count), then
        // pass-through non-line features.
        let mut out_schema = schema.clone();
        if out_schema.field_index("merged_count").is_none() {
            out_schema.add_field(FieldDef::new("merged_count", FieldType::Integer));
        }
        let mc_idx = out_schema.field_index("merged_count").unwrap();

        let mut out_features: Vec<Feature> =
            Vec::with_capacity(out_lines.len() + passthrough.len());
        for (coords, feat, count) in out_lines {
            let mut f = layer.features[feat].clone();
            f.geometry = Some(Geometry::LineString(coords));
            f.attributes.resize(out_schema.len(), FieldValue::Null);
            f.set_by_index(mc_idx, FieldValue::Integer(count as i64));
            f.fid = out_features.len() as u64;
            out_features.push(f);
        }
        for mut f in passthrough {
            f.attributes.resize(out_schema.len(), FieldValue::Null);
            f.set_by_index(mc_idx, FieldValue::Integer(1));
            f.fid = out_features.len() as u64;
            out_features.push(f);
        }

        ctx.progress.info(&format!(
            "merged {input_line_count} segment(s) into {} line(s) ({merged_chains} chain(s) dissolved)",
            out_features.len() - passthrough_len(&layer)
        ));

        let mut out_layer = wbvector::Layer::new(layer.name);
        out_layer.schema = out_schema;
        out_layer.crs = layer.crs;
        out_layer.geom_type = Some(GeometryType::LineString);
        out_layer.features = out_features;

        let feature_count = out_layer.len();
        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_line_count".to_string(), json!(input_line_count));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("merged_chains".to_string(), json!(merged_chains));
        Ok(ToolRunResult { outputs })
    }
}

/// Counts the pass-through (non-line) features of a layer (for the log line).
fn passthrough_len(layer: &wbvector::Layer) -> usize {
    layer
        .features
        .iter()
        .filter(|f| {
            !matches!(
                f.geometry,
                Some(Geometry::LineString(_)) | Some(Geometry::MultiLineString(_))
            )
        })
        .count()
}

/// A single polyline together with its source feature index.
struct LineRef {
    coords: Vec<Coord>,
    feat: usize,
}

/// Finds the one other unused line that continues a chain through node `n`.
///
/// Returns `Some((line_index, which_end))` when `n` is a mergeable pseudo-node:
/// exactly two distinct line-ends touch it, one of them the current chain, the
/// other an unused line whose dissolve key matches. `which_end` is that line's
/// end (0 or 1) incident to `n`. Returns `None` at real junctions (degree ≠ 2),
/// dead ends, attribute breaks, or when the continuation is already used.
fn next_link(
    n: usize,
    current: usize,
    incidence: &HashMap<usize, Vec<(usize, u8)>>,
    used: &[bool],
    dissolve_key: &impl Fn(usize) -> Vec<String>,
    seed_key: &[String],
) -> Option<(usize, u8)> {
    let inc = incidence.get(&n)?;
    if inc.len() != 2 {
        return None; // dead end or true junction
    }
    // The continuation is the incident end that is not the current line.
    let cand = inc.iter().find(|&&(li, _)| li != current)?;
    let (nl, end) = *cand;
    if nl == current || used[nl] {
        return None;
    }
    if dissolve_key(nl) != seed_key {
        return None; // attribute break
    }
    Some((nl, end))
}

/// Appends `src` onto `dst`, skipping the shared first vertex; `reversed`
/// consumes `src` back-to-front (so a line joined at its far end reads forward).
fn append_skip_first(dst: &mut Vec<Coord>, src: &[Coord], reversed: bool) {
    if reversed {
        for c in src.iter().rev().skip(1) {
            dst.push(c.clone());
        }
    } else {
        for c in src.iter().skip(1) {
            dst.push(c.clone());
        }
    }
}

/// Prepends `src` before `dst`, skipping the shared last vertex; `reversed`
/// consumes `src` back-to-front.
fn prepend_skip_last(dst: &mut Vec<Coord>, src: &[Coord], reversed: bool) {
    let mut head: Vec<Coord> = if reversed {
        src.iter().rev().take(src.len() - 1).cloned().collect()
    } else {
        src[..src.len() - 1].to_vec()
    };
    head.append(dst);
    *dst = head;
}

// ── Node index (endpoint snapping) ──────────────────────────────────────────

/// Maps snapped endpoint coordinates to integer node ids. Coordinates are
/// quantized onto a grid of cell `tol`, and each new grid cell (and its already
/// seen 8-neighbours, so points straddling a cell edge still merge) yields one
/// node id.
struct NodeIndex {
    tol: f64,
    map: HashMap<(i64, i64), usize>,
    next: usize,
}

impl NodeIndex {
    fn new(tol: f64) -> Self {
        Self {
            tol: tol.max(f64::MIN_POSITIVE),
            map: HashMap::new(),
            next: 0,
        }
    }

    fn cell(&self, c: &Coord) -> (i64, i64) {
        (
            (c.x / self.tol).round() as i64,
            (c.y / self.tol).round() as i64,
        )
    }

    fn id_of(&mut self, c: &Coord) -> usize {
        let (cx, cy) = self.cell(c);
        // Reuse a node from this cell or any of the 8 neighbours (handles points
        // that snap to adjacent grid cells).
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(&id) = self.map.get(&(cx + dx, cy + dy)) {
                    self.map.insert((cx, cy), id);
                    return id;
                }
            }
        }
        let id = self.next;
        self.next += 1;
        self.map.insert((cx, cy), id);
        id
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    dissolve_fields: Vec<String>,
    snap_tolerance: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let dissolve_fields = parse_optional_str(args, "dissolve_fields")?
        .map(|s| {
            s.split(',')
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(str::to_string)
                .collect::<Vec<_>>()
        })
        .unwrap_or_default();
    let snap_tolerance = parse_optional_f64(args, "snap_tolerance")?.unwrap_or(1e-6);
    if !(snap_tolerance > 0.0 && snap_tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'snap_tolerance' must be a positive number".to_string(),
        ));
    }
    Ok(Params {
        dissolve_fields,
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

fn field_value_string(v: &FieldValue) -> String {
    match v {
        FieldValue::Text(s) | FieldValue::Date(s) | FieldValue::DateTime(s) => s.clone(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) => f.to_string(),
        FieldValue::Boolean(b) => b.to_string(),
        FieldValue::Null | FieldValue::Blob(_) => String::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn line(pts: &[(f64, f64)]) -> Geometry {
        Geometry::LineString(pts.iter().map(|&(x, y)| Coord::xy(x, y)).collect())
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = MergeLinesByPseudoNodeTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn len_of(g: &Geometry) -> usize {
        match g {
            Geometry::LineString(cs) => cs.len(),
            _ => 0,
        }
    }

    /// Three segments in a straight chain (two pseudo-nodes) merge into one.
    #[test]
    fn merges_a_straight_chain() {
        let mut layer = Layer::new("roads");
        layer
            .add_feature(Some(line(&[(0.0, 0.0), (1.0, 0.0)])), &[])
            .unwrap();
        layer
            .add_feature(Some(line(&[(1.0, 0.0), (2.0, 0.0)])), &[])
            .unwrap();
        layer
            .add_feature(Some(line(&[(2.0, 0.0), (3.0, 0.0)])), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input }));
        assert_eq!(out.outputs["feature_count"], json!(1));
        assert_eq!(out.outputs["merged_chains"], json!(1));
        // 3 segments -> 4 vertices in one line.
        assert_eq!(len_of(layer.features[0].geometry.as_ref().unwrap()), 4);
        assert_eq!(
            layer.features[0]
                .get(&layer.schema, "merged_count")
                .unwrap()
                .as_i64(),
            Some(3)
        );
    }

    /// A T-junction (degree-3 node) is preserved: the three arms stay separate.
    #[test]
    fn preserves_true_intersection() {
        let mut layer = Layer::new("roads");
        // Two collinear arms and one branching arm all meet at (1,0).
        layer
            .add_feature(Some(line(&[(0.0, 0.0), (1.0, 0.0)])), &[])
            .unwrap();
        layer
            .add_feature(Some(line(&[(1.0, 0.0), (2.0, 0.0)])), &[])
            .unwrap();
        layer
            .add_feature(Some(line(&[(1.0, 0.0), (1.0, 1.0)])), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, _) = run_tool(json!({ "input": input }));
        // Degree-3 node -> nothing merges.
        assert_eq!(out.outputs["feature_count"], json!(3));
        assert_eq!(out.outputs["merged_chains"], json!(0));
    }

    /// `dissolve_fields` stops a merge where an attribute differs.
    #[test]
    fn dissolve_fields_break_a_chain() {
        let mut layer = Layer::new("roads");
        layer.add_field(FieldDef::new("cls", FieldType::Text));
        layer
            .add_feature(
                Some(line(&[(0.0, 0.0), (1.0, 0.0)])),
                &[("cls", "a".into())],
            )
            .unwrap();
        // pseudo-node at (1,0) but different class -> no merge.
        layer
            .add_feature(
                Some(line(&[(1.0, 0.0), (2.0, 0.0)])),
                &[("cls", "b".into())],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (with, _) = run_tool(json!({ "input": input, "dissolve_fields": "cls" }));
        assert_eq!(with.outputs["feature_count"], json!(2));
        // Without the field constraint they merge.
        let id2 = {
            let mut layer = Layer::new("roads");
            layer.add_field(FieldDef::new("cls", FieldType::Text));
            layer
                .add_feature(
                    Some(line(&[(0.0, 0.0), (1.0, 0.0)])),
                    &[("cls", "a".into())],
                )
                .unwrap();
            layer
                .add_feature(
                    Some(line(&[(1.0, 0.0), (2.0, 0.0)])),
                    &[("cls", "b".into())],
                )
                .unwrap();
            memory_store::put_vector(layer)
        };
        let (without, _) =
            run_tool(json!({ "input": memory_store::make_vector_memory_path(&id2) }));
        assert_eq!(without.outputs["feature_count"], json!(1));
    }

    /// A closed loop of pseudo-nodes collapses to a single closed line.
    #[test]
    fn closes_a_loop() {
        let mut layer = Layer::new("roads");
        layer
            .add_feature(Some(line(&[(0.0, 0.0), (2.0, 0.0)])), &[])
            .unwrap();
        layer
            .add_feature(Some(line(&[(2.0, 0.0), (2.0, 2.0)])), &[])
            .unwrap();
        layer
            .add_feature(Some(line(&[(2.0, 2.0), (0.0, 0.0)])), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input }));
        assert_eq!(out.outputs["feature_count"], json!(1));
        let g = layer.features[0].geometry.as_ref().unwrap();
        if let Geometry::LineString(cs) = g {
            assert_eq!(cs.first(), cs.last(), "loop should be closed");
        } else {
            panic!("expected LineString");
        }
    }

    /// Non-line features pass through untouched.
    #[test]
    fn passes_points_through() {
        let mut layer = Layer::new("mixed");
        layer
            .add_feature(Some(Geometry::point(5.0, 5.0)), &[])
            .unwrap();
        layer
            .add_feature(Some(line(&[(0.0, 0.0), (1.0, 0.0)])), &[])
            .unwrap();
        layer
            .add_feature(Some(line(&[(1.0, 0.0), (2.0, 0.0)])), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input }));
        // one merged line + one point.
        assert_eq!(out.outputs["feature_count"], json!(2));
        assert!(layer
            .features
            .iter()
            .any(|f| matches!(f.geometry, Some(Geometry::Point(_)))));
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = MergeLinesByPseudoNodeTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "x.geojson", "snap_tolerance": 0 })).is_err());
        assert!(bad(json!({ "input": "x.geojson" })).is_ok());
        assert!(bad(json!({ "input": "x.geojson", "snap_tolerance": "0.5" })).is_ok());
    }
}
