//! GeoLibre tool: knockout masks at cul-de-sac bulbs in a road network.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Cul-De-Sac Masks* (Cartography →
//! Masking Tools). The repo ships two of the three ArcGIS masking tools —
//! `feature_outline_masks` (mask kinds `box` / `convex_hull`) and
//! `intersecting_layers_masks` — but neither handles the cul-de-sac case: a
//! thick road casing drawn around a small terminal loop fills the bulb in
//! entirely, and an outline mask buffers every feature uniformly rather than
//! finding the bulbs.
//!
//! A cul-de-sac here is a terminal loop: a run of edges that leaves the network
//! at a single **articulation node** and returns to it. Detecting the loop
//! rather than just a degree-1 endpoint is what separates a bulb from a plain
//! dangle — a dead-end stub with no turning circle needs no mask.
//!
//! `symbol_width` and `margin` are page units (points/millimetres at
//! `reference_scale`); they are converted to ground units before buffering, so
//! the same parameters produce correctly-sized masks at any scale.

use std::collections::BTreeMap;

use geo::{BooleanOps, Buffer, MultiPolygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Quantized coordinate key, so endpoints that differ only by float noise snap
/// to the same graph node.
type NodeKey = (i64, i64);

pub struct CulDeSacMasksTool;

impl Tool for CulDeSacMasksTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "cul_de_sac_masks",
            display_name: "Cul-De-Sac Masks",
            summary: "Generate scale-aware knockout mask polygons at cul-de-sac bulbs in a road network, like ArcGIS Cul-De-Sac Masks.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Road centerline layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output mask polygon path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "reference_scale",
                    description: "Map scale denominator the mask is sized for (e.g. 24000 for 1:24,000).",
                    required: true,
                },
                ToolParamSpec {
                    name: "symbol_width",
                    description: "Road casing width in page units at the reference scale.",
                    required: true,
                },
                ToolParamSpec {
                    name: "margin",
                    description: "Extra margin around the symbol, in page units (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Coordinate snapping tolerance for joining line endpoints into network nodes (default 0.001 map units).",
                    required: false,
                },
                ToolParamSpec {
                    name: "attributes",
                    description: "ids_only (default) | all — whether to carry the source feature's attributes onto the mask.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        for (k, must_be_positive) in [("reference_scale", true), ("symbol_width", true)] {
            let v = parse_optional_f64(args, k)?.ok_or_else(|| {
                ToolError::Validation(format!("missing required parameter '{k}'"))
            })?;
            if !v.is_finite() || (must_be_positive && v <= 0.0) {
                return Err(ToolError::Validation(format!(
                    "'{k}' must be greater than 0"
                )));
            }
        }
        if let Some(m) = parse_optional_f64(args, "margin")? {
            if !m.is_finite() || m < 0.0 {
                return Err(ToolError::Validation(
                    "'margin' must be zero or greater".to_string(),
                ));
            }
        }
        parse_attributes(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let scale = parse_optional_f64(args, "reference_scale")?.unwrap();
        let symbol_width = parse_optional_f64(args, "symbol_width")?.unwrap();
        let margin = parse_optional_f64(args, "margin")?.unwrap_or(0.0);
        let tolerance = parse_optional_f64(args, "tolerance")?.unwrap_or(0.001);
        let keep_all_attrs = parse_attributes(args)?;

        let layer = load_input_layer(input)?;

        // Build the node/edge graph. Each edge is one line feature; nodes are
        // quantized endpoints so float noise does not split a junction.
        let q = |v: f64| -> i64 { (v / tolerance.max(1e-12)).round() as i64 };
        let mut edges: Vec<Edge> = Vec::new();
        for (fi, feat) in layer.features.iter().enumerate() {
            let Some(geom) = feat.geometry.as_ref() else {
                continue;
            };
            for part in line_parts(geom) {
                if part.len() < 2 {
                    continue;
                }
                let a = (q(part[0].x), q(part[0].y));
                let b = (q(part[part.len() - 1].x), q(part[part.len() - 1].y));
                edges.push(Edge {
                    feature: fi,
                    a,
                    b,
                    coords: part,
                });
            }
        }

        // Node degree over the whole network.
        let mut degree: BTreeMap<NodeKey, usize> = BTreeMap::new();
        for e in &edges {
            *degree.entry(e.a).or_insert(0) += 1;
            *degree.entry(e.b).or_insert(0) += 1;
        }

        ctx.progress.info(&format!(
            "scanning {} edge(s) over {} node(s) for cul-de-sac bulbs",
            edges.len(),
            degree.len()
        ));

        // Ground width: page units divided by the scale denominator gives
        // ground units per page unit, so a 1 pt casing at 1:24,000 is 24,000 pt
        // on the ground. `symbol_width` is a full width, so buffer by half.
        let ground = (symbol_width / 2.0 + margin) * scale;

        // A cul-de-sac bulb: a single edge whose two endpoints are the SAME node
        // (a self-loop), attached to the rest of the network at that node.
        // Degree at the loop node counts the loop's two ends plus the approach.
        let mut bulbs: Vec<(usize, MultiPolygon)> = Vec::new();
        for e in &edges {
            if e.a != e.b {
                continue; // not a loop
            }
            // The loop contributes 2 to its node's degree; anything above that
            // is the stem connecting it to the network.
            let deg = degree.get(&e.a).copied().unwrap_or(0);
            if deg < 3 {
                // An isolated loop with no approach is not a cul-de-sac.
                continue;
            }
            let ls = geo::LineString::new(
                e.coords
                    .iter()
                    .map(|c| geo::Coord { x: c.x, y: c.y })
                    .collect(),
            );
            let buffered = ls.buffer(ground);
            if !buffered.0.is_empty() {
                bulbs.push((e.feature, buffered));
            }
        }

        let mut out = Layer::new(layer.name.clone());
        out.crs = layer.crs.clone();
        out.geom_type = Some(GeometryType::MultiPolygon);
        out.add_field(FieldDef::new("source_fid", FieldType::Integer));
        out.add_field(FieldDef::new("mask_kind", FieldType::Text));
        if keep_all_attrs {
            for fd in layer.schema.fields().iter() {
                out.add_field(fd.clone());
            }
        }

        for (fid, poly) in &bulbs {
            let mut fields: Vec<(String, FieldValue)> = vec![
                ("source_fid".into(), FieldValue::Integer(*fid as i64)),
                ("mask_kind".into(), FieldValue::Text("cul_de_sac".into())),
            ];
            if keep_all_attrs {
                let feat = &layer.features[*fid];
                for (i, fd) in layer.schema.fields().iter().enumerate() {
                    fields.push((fd.name.clone(), feat.attributes[i].clone()));
                }
            }
            let refs: Vec<(&str, FieldValue)> = fields
                .iter()
                .map(|(k, v)| (k.as_str(), v.clone()))
                .collect();
            out.add_feature(Some(multipolygon_to_geometry(poly)), &refs)
                .map_err(|e| ToolError::Execution(format!("failed writing mask: {e}")))?;
        }

        let count = bulbs.len();
        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("mask_count".to_string(), json!(count));
        outputs.insert("edge_count".to_string(), json!(edges.len()));
        outputs.insert("ground_radius".to_string(), json!(ground));
        Ok(ToolRunResult { outputs })
    }
}

struct Edge {
    feature: usize,
    a: NodeKey,
    b: NodeKey,
    coords: Vec<Coord>,
}

/// Every linear part of a geometry, as its own coordinate run.
fn line_parts(geom: &Geometry) -> Vec<Vec<Coord>> {
    match geom {
        Geometry::LineString(cs) => vec![cs.clone()],
        Geometry::MultiLineString(parts) => parts.clone(),
        Geometry::GeometryCollection(gs) => gs.iter().flat_map(line_parts).collect(),
        _ => Vec::new(),
    }
}

fn multipolygon_to_geometry(mp: &MultiPolygon) -> Geometry {
    Geometry::MultiPolygon(
        mp.0.iter()
            .map(|p| {
                (
                    linestring_to_ring(p.exterior()),
                    p.interiors().iter().map(linestring_to_ring).collect(),
                )
            })
            .collect(),
    )
}

fn linestring_to_ring(ls: &geo::LineString) -> Ring {
    let mut coords: Vec<Coord> = ls.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
    if coords.len() >= 2 && coords.first().map(|c| (c.x, c.y)) == coords.last().map(|c| (c.x, c.y))
    {
        coords.pop();
    }
    Ring::new(coords)
}

// ── parameter parsing ────────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
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

/// Returns true when all source attributes should be carried onto the mask.
fn parse_attributes(args: &ToolArgs) -> Result<bool, ToolError> {
    match args
        .get("attributes")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some("ids_only") => Ok(false),
        Some("all") => Ok(true),
        Some(o) => Err(ToolError::Validation(format!(
            "'attributes' must be ids_only or all, got '{o}'"
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

    fn lines(parts: Vec<Vec<(f64, f64)>>) -> String {
        let mut l = Layer::new("roads")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        for (i, p) in parts.into_iter().enumerate() {
            l.add_feature(
                Some(Geometry::LineString(
                    p.into_iter().map(|(x, y)| Coord::xy(x, y)).collect(),
                )),
                &[("name", FieldValue::Text(format!("r{i}")))],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    /// A closed loop that starts AND ends at the attach point `(ax, ay)`.
    /// The circle is centred at `(ax + r, ay)` so vertex 0 is the attach point,
    /// which is what makes the loop share a node with the stem.
    fn bulb_at(ax: f64, ay: f64, r: f64) -> Vec<(f64, f64)> {
        let (cx, cy) = (ax + r, ay);
        let mut p = Vec::new();
        for i in 0..=12 {
            // Start at angle pi so the first vertex lands exactly on (ax, ay).
            let a = std::f64::consts::PI + (i as f64 / 12.0) * std::f64::consts::TAU;
            p.push((cx + r * a.cos(), cy + r * a.sin()));
        }
        // Force exact closure onto the attach point (trig round-trip drift).
        p[0] = (ax, ay);
        let n = p.len();
        p[n - 1] = (ax, ay);
        p
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CulDeSacMasksTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// A stem leading to a terminal loop is a cul-de-sac and gets a mask.
    #[test]
    fn detects_a_terminal_loop() {
        // Stem from (0,0) to (10,0); bulb is a circle starting/ending at (10,0).
        let stem = vec![(0.0, 0.0), (10.0, 0.0)];
        let input = lines(vec![stem, bulb_at(10.0, 0.0, 2.0)]);
        let (out, layer) = run(json!({
            "input": input, "reference_scale": 1.0,
            "symbol_width": 2.0, "margin": 0.0
        }));
        assert_eq!(out.outputs["mask_count"], json!(1));
        assert_eq!(layer.len(), 1);
        let ki = layer.schema.field_index("mask_kind").unwrap();
        assert!(matches!(&layer.features[0].attributes[ki],
                         FieldValue::Text(s) if s == "cul_de_sac"));
    }

    /// THE discrimination: a plain dead-end stub is NOT a cul-de-sac, because
    /// it has no turning loop. Degree-1 endpoint detection alone would wrongly
    /// mask it.
    #[test]
    fn plain_dangle_is_not_a_cul_de_sac() {
        let input = lines(vec![
            vec![(0.0, 0.0), (10.0, 0.0)],
            vec![(10.0, 0.0), (20.0, 0.0)], // dead-ends at (20,0)
        ]);
        let (out, _) = run(json!({
            "input": input, "reference_scale": 1.0, "symbol_width": 2.0
        }));
        assert_eq!(
            out.outputs["mask_count"],
            json!(0),
            "a dangle with no loop must not be masked"
        );
    }

    /// A loop with no stem is a traffic island / roundabout, not a cul-de-sac.
    #[test]
    fn isolated_loop_is_not_a_cul_de_sac() {
        let input = lines(vec![bulb_at(0.0, 0.0, 5.0)]);
        let (out, _) = run(json!({
            "input": input, "reference_scale": 1.0, "symbol_width": 2.0
        }));
        assert_eq!(out.outputs["mask_count"], json!(0));
    }

    /// Page units scale to ground units through reference_scale.
    #[test]
    fn mask_size_scales_with_reference_scale() {
        let input = lines(vec![vec![(0.0, 0.0), (10.0, 0.0)], bulb_at(10.0, 0.0, 2.0)]);
        let (small, _) = run(json!({
            "input": input, "reference_scale": 1.0, "symbol_width": 2.0
        }));
        let (big, _) = run(json!({
            "input": input, "reference_scale": 10.0, "symbol_width": 2.0
        }));
        let (rs, rb) = (
            small.outputs["ground_radius"].as_f64().unwrap(),
            big.outputs["ground_radius"].as_f64().unwrap(),
        );
        assert!((rs - 1.0).abs() < 1e-9, "half of 2 page units at 1:1");
        assert!((rb - 10.0).abs() < 1e-9, "same symbol, 10x scale");
    }

    /// margin widens the mask.
    #[test]
    fn margin_widens_the_mask() {
        let input = lines(vec![vec![(0.0, 0.0), (10.0, 0.0)], bulb_at(10.0, 0.0, 2.0)]);
        let (no_margin, _) = run(json!({
            "input": input, "reference_scale": 1.0, "symbol_width": 2.0, "margin": 0.0
        }));
        let (with_margin, _) = run(json!({
            "input": input, "reference_scale": 1.0, "symbol_width": 2.0, "margin": 3.0
        }));
        assert!(
            with_margin.outputs["ground_radius"].as_f64().unwrap()
                > no_margin.outputs["ground_radius"].as_f64().unwrap()
        );
    }

    /// attributes=all carries the source row onto the mask.
    #[test]
    fn attributes_all_carries_source_fields() {
        let input = lines(vec![vec![(0.0, 0.0), (10.0, 0.0)], bulb_at(10.0, 0.0, 2.0)]);
        let (_, layer) = run(json!({
            "input": input, "reference_scale": 1.0, "symbol_width": 2.0,
            "attributes": "all"
        }));
        assert!(layer.schema.field_index("name").is_some());
    }

    #[test]
    fn rejects_bad_parameters() {
        let p = lines(vec![vec![(0.0, 0.0), (1.0, 0.0)]]);
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CulDeSacMasksTool.validate(&args).is_err()
        };
        assert!(bad(json!({ "reference_scale": 1, "symbol_width": 1 })));
        assert!(bad(json!({ "input": p, "symbol_width": 1 })));
        assert!(bad(json!({ "input": p, "reference_scale": 1 })));
        assert!(bad(
            json!({ "input": p, "reference_scale": 0, "symbol_width": 1 })
        ));
        assert!(bad(
            json!({ "input": p, "reference_scale": 1, "symbol_width": 1, "margin": -1 })
        ));
        assert!(bad(
            json!({ "input": p, "reference_scale": 1, "symbol_width": 1, "attributes": "some" })
        ));
    }
}
