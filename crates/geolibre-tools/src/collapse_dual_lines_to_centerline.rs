//! GeoLibre tool: collapse paired dual carriageways into single centerlines.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Collapse Dual Lines To Centerline*
//! (Cartography; related to *Merge Divided Roads*). Completes the road-
//! generalization arc started by `thin_road_network`: the bundled
//! `river_centerlines` derives a centerline from polygons or a raster, not from
//! **paired line features**, and nothing bundled collapses dual-line geometry.
//! Essential for generalizing OSM dual carriageways (each direction a separate
//! way) for small-scale mapping.
//!
//! The algorithm:
//!
//! 1. Densify every line at `sample_distance` so parallelism can be measured
//!    vertex-by-vertex.
//! 2. For each candidate pair of lines (bounding boxes within `max_width`, and
//!    an optional `attribute` value matching), measure the *directed overlap* in
//!    both directions — the fraction of one line's vertices whose perpendicular
//!    distance to the other line lands in `[min_width, max_width]`. A pair whose
//!    smaller overlap clears `min_overlap` is a dual carriageway.
//! 3. Greedily accept mutually-best pairs (each line collapses at most once),
//!    largest overlap first.
//! 4. The centerline is the ordered midpoints between each vertex of one line
//!    and its nearest point on the other, over the overlapping span.
//! 5. Reconnect: snap each new centerline endpoint to the nearest surviving line
//!    endpoint within `max_width`, so the network stays routable.
//!
//! Unpaired lines pass through unchanged. Distances are in the layer's CRS units
//! (use a projected CRS). Non-line geometries pass through.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldValue, Geometry, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CollapseDualLinesToCenterlineTool;

impl Tool for CollapseDualLinesToCenterlineTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "collapse_dual_lines_to_centerline",
            display_name: "Collapse Dual Lines To Centerline",
            summary: "Detect pairs of roughly parallel lines within a width tolerance (dual carriageways, road casings) and replace each pair with a single centerline, reconnecting junctions — like ArcGIS Collapse Dual Lines To Centerline.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input line vector layer (projected CRS; widths in its units).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output line vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_width",
                    description: "Minimum separation between paired lines, in CRS units (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_width",
                    description: "Maximum separation between paired lines, in CRS units. Required (the collapse width).",
                    required: true,
                },
                ToolParamSpec {
                    name: "attribute",
                    description: "Optional field that must match for two lines to pair (e.g. road name or class).",
                    required: false,
                },
                ToolParamSpec {
                    name: "sample_distance",
                    description: "Densification interval for parallelism testing, in CRS units. Default: max_width / 2.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_overlap",
                    description: "Minimum fraction of a line that must run parallel to its partner to pair (0-1, default 0.5).",
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
        let attr_idx = match &prm.attribute {
            Some(a) => Some(layer.schema.field_index(a).ok_or_else(|| {
                ToolError::Validation(format!("attribute field '{a}' not found"))
            })?),
            None => None,
        };

        // Extract line features (index, densified polyline, attribute key).
        let mut lines: Vec<LineRec> = Vec::new();
        for (fidx, feature) in layer.features.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            for chain in line_chains(geom) {
                if chain.len() < 2 {
                    continue;
                }
                let dense = densify(&chain, prm.sample_distance);
                let bbox = bbox_of(&dense);
                let key = attr_idx.map(|i| field_key(&feature.attributes[i]));
                lines.push(LineRec {
                    feature: fidx,
                    raw: chain,
                    dense,
                    bbox,
                    key,
                });
            }
        }

        ctx.progress
            .info(&format!("{} line(s) to test for dual pairs", lines.len()));

        // ── Candidate pairs (bbox-prefiltered) with directed overlap ──────────
        let mut candidates: Vec<(f64, usize, usize)> = Vec::new();
        for i in 0..lines.len() {
            for j in (i + 1)..lines.len() {
                if lines[i].key != lines[j].key {
                    continue;
                }
                if !bbox_within(&lines[i].bbox, &lines[j].bbox, prm.max_width) {
                    continue;
                }
                let ov_ij = directed_overlap(&lines[i].dense, &lines[j].raw, &prm);
                if ov_ij < prm.min_overlap {
                    continue;
                }
                let ov_ji = directed_overlap(&lines[j].dense, &lines[i].raw, &prm);
                let ov = ov_ij.min(ov_ji);
                if ov >= prm.min_overlap {
                    candidates.push((ov, i, j));
                }
            }
        }
        // Greedy: strongest overlap first, each line used once.
        candidates.sort_by(|a, b| b.0.total_cmp(&a.0));
        let mut consumed = vec![false; lines.len()];
        let mut centerlines: Vec<Vec<P>> = Vec::new();
        let mut pair_features: Vec<usize> = Vec::new();
        for (_ov, i, j) in &candidates {
            if consumed[*i] || consumed[*j] {
                continue;
            }
            let cl = centerline(&lines[*i].dense, &lines[*j].raw, &prm);
            if cl.len() >= 2 {
                consumed[*i] = true;
                consumed[*j] = true;
                centerlines.push(cl);
                pair_features.push(lines[*i].feature);
            }
        }

        // ── Assemble output: centerlines + pass-through (unpaired) lines ──────
        // Collect surviving endpoints for reconnection snapping.
        let mut survivor_endpoints: Vec<P> = Vec::new();
        for (li, rec) in lines.iter().enumerate() {
            if !consumed[li] {
                survivor_endpoints.push(rec.raw[0]);
                survivor_endpoints.push(rec.raw[rec.raw.len() - 1]);
            }
        }
        // Snap each centerline end to the nearest survivor endpoint within width.
        for cl in centerlines.iter_mut() {
            snap_end(cl, 0, &survivor_endpoints, prm.max_width);
            let last = cl.len() - 1;
            snap_end(cl, last, &survivor_endpoints, prm.max_width);
        }

        let mut out = Layer::new("centerlines").with_geom_type(wbvector::GeometryType::LineString);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(wbvector::FieldDef::new("source", wbvector::FieldType::Text));

        let mut centerline_count = 0usize;
        for cl in &centerlines {
            let coords: Vec<Coord> = cl.iter().map(|p| Coord::xy(p.x, p.y)).collect();
            out.add_feature(
                Some(Geometry::line_string(coords)),
                &[("source", "centerline".into())],
            )
            .map_err(|e| ToolError::Execution(format!("failed writing centerline: {e}")))?;
            centerline_count += 1;
        }
        let mut passthrough = 0usize;
        for (li, rec) in lines.iter().enumerate() {
            if consumed[li] {
                continue;
            }
            let coords: Vec<Coord> = rec.raw.iter().map(|p| Coord::xy(p.x, p.y)).collect();
            out.add_feature(
                Some(Geometry::line_string(coords)),
                &[("source", "kept".into())],
            )
            .map_err(|e| ToolError::Execution(format!("failed writing line: {e}")))?;
            passthrough += 1;
        }

        ctx.progress.info(&format!(
            "{centerline_count} centerline(s) from {} collapsed line(s); {passthrough} kept",
            centerline_count * 2
        ));

        let feature_count = out.len();
        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_lines".to_string(), json!(lines.len()));
        outputs.insert("pairs_collapsed".to_string(), json!(centerline_count));
        outputs.insert("centerline_count".to_string(), json!(centerline_count));
        outputs.insert("passthrough_count".to_string(), json!(passthrough));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        Ok(ToolRunResult { outputs })
    }
}

// ── Line pairing ─────────────────────────────────────────────────────────────

struct LineRec {
    feature: usize,
    raw: Vec<P>,
    dense: Vec<P>,
    bbox: [f64; 4], // xmin, ymin, xmax, ymax
    key: Option<String>,
}

/// Fraction of `a`'s vertices whose nearest point on polyline `b` lies within
/// `[min_width, max_width]` — a directed measure of how much of `a` runs
/// parallel to `b` at carriageway spacing.
fn directed_overlap(a: &[P], b: &[P], prm: &Params) -> f64 {
    if a.is_empty() || b.len() < 2 {
        return 0.0;
    }
    let mut within = 0usize;
    for &p in a {
        let (d, _) = nearest_on_polyline(p, b);
        if d >= prm.min_width && d <= prm.max_width {
            within += 1;
        }
    }
    within as f64 / a.len() as f64
}

/// The centerline of a dual pair: for each vertex of `a` whose nearest point on
/// `b` is within the width band, the midpoint between them, in `a`'s order.
fn centerline(a: &[P], b: &[P], prm: &Params) -> Vec<P> {
    let mut out: Vec<P> = Vec::new();
    for &p in a {
        let (d, q) = nearest_on_polyline(p, b);
        if d >= prm.min_width && d <= prm.max_width {
            let mid = P {
                x: (p.x + q.x) * 0.5,
                y: (p.y + q.y) * 0.5,
            };
            if out.last().is_none_or(|l| dist(*l, mid) > 1e-9) {
                out.push(mid);
            }
        }
    }
    out
}

/// Nearest point on polyline `b` to `p`, and its distance.
fn nearest_on_polyline(p: P, b: &[P]) -> (f64, P) {
    let mut best_d = f64::INFINITY;
    let mut best = b[0];
    for w in b.windows(2) {
        let q = project_point_seg(p, w[0], w[1]);
        let d = dist(p, q);
        if d < best_d {
            best_d = d;
            best = q;
        }
    }
    (best_d, best)
}

fn project_point_seg(p: P, a: P, b: P) -> P {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let len2 = dx * dx + dy * dy;
    if len2 <= 0.0 {
        return a;
    }
    let t = (((p.x - a.x) * dx + (p.y - a.y) * dy) / len2).clamp(0.0, 1.0);
    P {
        x: a.x + t * dx,
        y: a.y + t * dy,
    }
}

/// Moves endpoint `idx` of `cl` to the nearest survivor endpoint within `width`.
fn snap_end(cl: &mut [P], idx: usize, endpoints: &[P], width: f64) {
    let p = cl[idx];
    let mut best_d = width;
    let mut best: Option<P> = None;
    for &e in endpoints {
        let d = dist(p, e);
        if d < best_d {
            best_d = d;
            best = Some(e);
        }
    }
    if let Some(e) = best {
        cl[idx] = e;
    }
}

// ── Geometry helpers ─────────────────────────────────────────────────────────

#[derive(Clone, Copy)]
struct P {
    x: f64,
    y: f64,
}

fn dist(a: P, b: P) -> f64 {
    (a.x - b.x).hypot(a.y - b.y)
}

fn line_chains(geom: &Geometry) -> Vec<Vec<P>> {
    match geom {
        Geometry::LineString(cs) => vec![to_pts(cs)],
        Geometry::MultiLineString(lines) => lines.iter().map(|l| to_pts(l)).collect(),
        _ => Vec::new(),
    }
}

fn to_pts(cs: &[Coord]) -> Vec<P> {
    let mut out: Vec<P> = Vec::with_capacity(cs.len());
    for c in cs {
        let p = P { x: c.x, y: c.y };
        if out.last().is_none_or(|l| dist(*l, p) > 1e-12) {
            out.push(p);
        }
    }
    out
}

fn densify(pts: &[P], step: f64) -> Vec<P> {
    if pts.len() < 2 || step <= 0.0 {
        return pts.to_vec();
    }
    let mut out = Vec::with_capacity(pts.len() * 2);
    for w in pts.windows(2) {
        let (a, b) = (w[0], w[1]);
        out.push(a);
        let d = dist(a, b);
        let pieces = (d / step).ceil().max(1.0) as usize;
        for j in 1..pieces {
            let t = j as f64 / pieces as f64;
            out.push(P {
                x: a.x + (b.x - a.x) * t,
                y: a.y + (b.y - a.y) * t,
            });
        }
    }
    out.push(pts[pts.len() - 1]);
    out
}

fn bbox_of(pts: &[P]) -> [f64; 4] {
    let mut b = [
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    ];
    for p in pts {
        b[0] = b[0].min(p.x);
        b[1] = b[1].min(p.y);
        b[2] = b[2].max(p.x);
        b[3] = b[3].max(p.y);
    }
    b
}

fn bbox_within(a: &[f64; 4], b: &[f64; 4], pad: f64) -> bool {
    a[0] - pad <= b[2] && b[0] - pad <= a[2] && a[1] - pad <= b[3] && b[1] - pad <= a[3]
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
    min_width: f64,
    max_width: f64,
    attribute: Option<String>,
    sample_distance: f64,
    min_overlap: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let max_width = parse_optional_f64(args, "max_width")?.ok_or_else(|| {
        ToolError::Validation("required parameter 'max_width' is missing".to_string())
    })?;
    if !(max_width > 0.0 && max_width.is_finite()) {
        return Err(ToolError::Validation(
            "'max_width' must be a positive number".to_string(),
        ));
    }
    let min_width = parse_optional_f64(args, "min_width")?.unwrap_or(0.0);
    if !(min_width >= 0.0 && min_width < max_width) {
        return Err(ToolError::Validation(
            "'min_width' must be >= 0 and < 'max_width'".to_string(),
        ));
    }
    let attribute = parse_optional_str(args, "attribute")?.map(str::to_string);
    let sample_distance = parse_optional_f64(args, "sample_distance")?.unwrap_or(max_width / 2.0);
    if !(sample_distance > 0.0 && sample_distance.is_finite()) {
        return Err(ToolError::Validation(
            "'sample_distance' must be a positive number".to_string(),
        ));
    }
    let min_overlap = parse_optional_f64(args, "min_overlap")?.unwrap_or(0.5);
    if !(0.0..=1.0).contains(&min_overlap) {
        return Err(ToolError::Validation(
            "'min_overlap' must be between 0 and 1".to_string(),
        ));
    }
    Ok(Params {
        min_width,
        max_width,
        attribute,
        sample_distance,
        min_overlap,
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

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn line(coords: &[(f64, f64)]) -> Geometry {
        Geometry::line_string(coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect())
    }

    fn layer_of(geoms: Vec<Geometry>) -> String {
        let mut l = Layer::new("lines")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(26918);
        for g in geoms {
            l.add_feature(Some(g), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CollapseDualLinesToCenterlineTool
            .run(&args, &ctx())
            .unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// Two parallel carriageways 10 apart collapse to a centerline between them.
    #[test]
    fn collapses_parallel_pair() {
        // Both run east along y=0 and y=10 from x=0..100.
        let north = line(&[(0.0, 10.0), (100.0, 10.0)]);
        let south = line(&[(0.0, 0.0), (100.0, 0.0)]);
        let input = layer_of(vec![north, south]);
        let (out, layer) = run(json!({ "input": input, "max_width": 15.0, "min_overlap": 0.5 }));
        assert_eq!(out.outputs["pairs_collapsed"], json!(1));
        assert_eq!(out.outputs["passthrough_count"], json!(0));
        // The centerline should sit at y=5.
        let cl = match layer.features[0].geometry.as_ref().unwrap() {
            Geometry::LineString(cs) => cs.clone(),
            other => panic!("expected line, got {other:?}"),
        };
        for c in &cl {
            assert!((c.y - 5.0).abs() < 1e-6, "centerline not at y=5: {}", c.y);
        }
    }

    /// Lines too far apart (beyond max_width) are not paired.
    #[test]
    fn far_lines_are_not_paired() {
        let north = line(&[(0.0, 50.0), (100.0, 50.0)]);
        let south = line(&[(0.0, 0.0), (100.0, 0.0)]);
        let input = layer_of(vec![north, south]);
        let (out, _l) = run(json!({ "input": input, "max_width": 15.0 }));
        assert_eq!(out.outputs["pairs_collapsed"], json!(0));
        assert_eq!(out.outputs["passthrough_count"], json!(2));
    }

    /// Crossing (non-parallel) lines are not paired.
    #[test]
    fn crossing_lines_are_not_paired() {
        let a = line(&[(0.0, 0.0), (100.0, 100.0)]);
        let b = line(&[(0.0, 100.0), (100.0, 0.0)]);
        let input = layer_of(vec![a, b]);
        let (out, _l) = run(json!({ "input": input, "max_width": 15.0 }));
        assert_eq!(out.outputs["pairs_collapsed"], json!(0));
    }

    /// The attribute filter prevents pairing lines with different values.
    #[test]
    fn attribute_must_match_to_pair() {
        let mut l = Layer::new("lines")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(26918);
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.add_feature(
            Some(line(&[(0.0, 10.0), (100.0, 10.0)])),
            &[("name", "A".into())],
        )
        .unwrap();
        l.add_feature(
            Some(line(&[(0.0, 0.0), (100.0, 0.0)])),
            &[("name", "B".into())],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);
        let (out, _l) = run(json!({ "input": input, "max_width": 15.0, "attribute": "name" }));
        assert_eq!(
            out.outputs["pairs_collapsed"],
            json!(0),
            "different names must not pair"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CollapseDualLinesToCenterlineTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err()); // no max_width
        assert!(bad(json!({ "input": "a.geojson", "max_width": 0 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "max_width": 10, "min_width": 20 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "max_width": 10 })).is_ok());
    }
}
