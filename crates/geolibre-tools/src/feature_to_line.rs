//! GeoLibre tool: planarize line / polygon-boundary features into noded,
//! non-overlapping line segments split at every mutual and self intersection.
//!
//! Pure-Rust counterpart of ArcGIS Data Management's *Feature To Line*. The
//! bundled `polygons_to_lines` only cracks polygon rings and `split_vector_lines`
//! splits by a separate break layer; neither nodes an arbitrary set of lines at
//! their crossings — the standard precursor to topology building, routing
//! graphs, and clean overlay. `merge_lines_by_pseudo_node` is the inverse.
//!
//! Every input line and polygon boundary is decomposed into segments; each
//! segment is split at all points where another segment crosses it or touches
//! its interior (T-junctions and collinear-overlap endpoints included). The
//! result is a set of noded two-vertex line features, each tagged with the fid
//! of the source feature it came from.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// One directed segment with the source feature index it came from.
struct Seg {
    a: Coord,
    b: Coord,
    src: usize,
}

pub struct FeatureToLineTool;

impl Tool for FeatureToLineTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "feature_to_line",
            display_name: "Feature To Line",
            summary: "Planarize line/polygon-boundary features into noded, non-overlapping line segments split at every mutual and self intersection, like ArcGIS Feature To Line — the topology-building precursor the bundled polygons_to_lines and split_vector_lines don't provide.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input line or polygon vector layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output line vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cluster_tolerance",
                    description: "Snap split points to this map-unit grid to merge near-coincident nodes (default 1e-9).",
                    required: false,
                },
                ToolParamSpec {
                    name: "attributes",
                    description: "When true (default), tag each output segment with 'src_fid', the source feature id.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_optional_f64(args, "cluster_tolerance")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let tol = parse_optional_f64(args, "cluster_tolerance")?
            .unwrap_or(1e-9)
            .max(0.0);
        let attributes = parse_optional_bool(args, "attributes")?.unwrap_or(true);

        let layer = load_input_layer(input)?;

        // Collect all input segments with their source feature index.
        let mut segs: Vec<Seg> = Vec::new();
        for (fi, f) in layer.features.iter().enumerate() {
            if let Some(g) = f.geometry.as_ref() {
                collect_segments(g, fi, &mut segs);
            }
        }
        if segs.is_empty() {
            return Err(ToolError::Execution(
                "input has no line or polygon-boundary geometry".to_string(),
            ));
        }

        // Per-segment split parameters (t in (0,1)); 0/1 always implied.
        let mut splits: Vec<Vec<f64>> = vec![Vec::new(); segs.len()];

        // 1-D sweep on x to prune the pairwise intersection search.
        let mut order: Vec<usize> = (0..segs.len()).collect();
        order.sort_by(|&i, &j| {
            seg_min_x(&segs[i])
                .partial_cmp(&seg_min_x(&segs[j]))
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        ctx.progress
            .info(&format!("noding {} segment(s)", segs.len()));
        for (oi, &i) in order.iter().enumerate() {
            let i_maxx = seg_max_x(&segs[i]) + tol;
            for &j in order.iter().skip(oi + 1) {
                if seg_min_x(&segs[j]) - tol > i_maxx {
                    break; // no further segment can overlap in x
                }
                if !y_overlap(&segs[i], &segs[j], tol) {
                    continue;
                }
                node_pair(&segs, i, j, tol, &mut splits);
            }
        }

        // Emit split sub-segments.
        let mut out = Layer::new("feature_to_line").with_geom_type(GeometryType::LineString);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        if attributes {
            out.add_field(FieldDef::new("src_fid", FieldType::Integer));
        }

        let mut n_out = 0usize;
        for (si, seg) in segs.iter().enumerate() {
            let mut ts: Vec<f64> = splits[si]
                .iter()
                .copied()
                .filter(|t| *t > 1e-12 && *t < 1.0 - 1e-12)
                .collect();
            ts.push(0.0);
            ts.push(1.0);
            ts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            ts.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
            for w in ts.windows(2) {
                let p0 = lerp(&seg.a, &seg.b, w[0]);
                let p1 = lerp(&seg.a, &seg.b, w[1]);
                if (p0.x - p1.x).abs() < 1e-15 && (p0.y - p1.y).abs() < 1e-15 {
                    continue; // zero-length
                }
                let attrs: Vec<(&str, FieldValue)> = if attributes {
                    vec![("src_fid", FieldValue::Integer(seg.src as i64))]
                } else {
                    vec![]
                };
                out.add_feature(Some(Geometry::LineString(vec![p0, p1])), &attrs)
                    .map_err(|e| ToolError::Execution(format!("failed adding line: {e}")))?;
                n_out += 1;
            }
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_segments".to_string(), json!(segs.len()));
        outputs.insert("output_lines".to_string(), json!(n_out));
        Ok(ToolRunResult { outputs })
    }
}

/// Records the split parameter(s) where segments `i` and `j` meet, on whichever
/// segment's interior the meeting point falls.
fn node_pair(segs: &[Seg], i: usize, j: usize, tol: f64, splits: &mut [Vec<f64>]) {
    let (a, b) = (&segs[i], &segs[j]);
    let denom = (a.b.x - a.a.x) * (b.b.y - b.a.y) - (a.b.y - a.a.y) * (b.b.x - b.a.x);
    if denom.abs() > 1e-12 {
        // Proper (non-parallel) intersection.
        let t = ((b.a.x - a.a.x) * (b.b.y - b.a.y) - (b.a.y - a.a.y) * (b.b.x - b.a.x)) / denom;
        let u = ((b.a.x - a.a.x) * (a.b.y - a.a.y) - (b.a.y - a.a.y) * (a.b.x - a.a.x)) / denom;
        if (-1e-9..=1.0 + 1e-9).contains(&t) && (-1e-9..=1.0 + 1e-9).contains(&u) {
            splits[i].push(t.clamp(0.0, 1.0));
            splits[j].push(u.clamp(0.0, 1.0));
        }
    } else {
        // Parallel/collinear: project each segment's endpoints onto the other.
        add_endpoint_splits(a, b, i, tol, splits);
        add_endpoint_splits(b, a, j, tol, splits);
    }
}

/// For the collinear case, project `other`'s endpoints onto `seg` and record
/// interior split params where they lie on the segment.
fn add_endpoint_splits(seg: &Seg, other: &Seg, seg_idx: usize, tol: f64, splits: &mut [Vec<f64>]) {
    for p in [&other.a, &other.b] {
        if let Some(t) = project_on(seg, p, tol) {
            splits[seg_idx].push(t);
        }
    }
}

/// Returns the parameter `t` of point `p`'s projection on `seg` when `p` lies on
/// the segment (within `tol`) and strictly interior.
fn project_on(seg: &Seg, p: &Coord, tol: f64) -> Option<f64> {
    let dx = seg.b.x - seg.a.x;
    let dy = seg.b.y - seg.a.y;
    let len2 = dx * dx + dy * dy;
    if len2 <= f64::EPSILON {
        return None;
    }
    let t = ((p.x - seg.a.x) * dx + (p.y - seg.a.y) * dy) / len2;
    if !(1e-9..=1.0 - 1e-9).contains(&t) {
        return None;
    }
    let px = seg.a.x + t * dx;
    let py = seg.a.y + t * dy;
    if (p.x - px).hypot(p.y - py) <= tol.max(1e-9) {
        Some(t)
    } else {
        None
    }
}

fn collect_segments(geom: &Geometry, src: usize, out: &mut Vec<Seg>) {
    fn push(coords: &[Coord], src: usize, out: &mut Vec<Seg>) {
        for w in coords.windows(2) {
            if w[0] != w[1] {
                out.push(Seg {
                    a: w[0].clone(),
                    b: w[1].clone(),
                    src,
                });
            }
        }
    }
    match geom {
        Geometry::LineString(cs) => push(cs, src, out),
        Geometry::MultiLineString(ls) => ls.iter().for_each(|l| push(l, src, out)),
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            push(&exterior.0, src, out);
            interiors.iter().for_each(|r| push(&r.0, src, out));
        }
        Geometry::MultiPolygon(ps) => {
            for (e, hs) in ps {
                push(&e.0, src, out);
                hs.iter().for_each(|r| push(&r.0, src, out));
            }
        }
        Geometry::GeometryCollection(gs) => gs.iter().for_each(|g| collect_segments(g, src, out)),
        _ => {}
    }
}

fn lerp(a: &Coord, b: &Coord, t: f64) -> Coord {
    Coord::xy(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t)
}

fn seg_min_x(s: &Seg) -> f64 {
    s.a.x.min(s.b.x)
}
fn seg_max_x(s: &Seg) -> f64 {
    s.a.x.max(s.b.x)
}
fn y_overlap(a: &Seg, b: &Seg, tol: f64) -> bool {
    let (a0, a1) = (a.a.y.min(a.b.y), a.a.y.max(a.b.y));
    let (b0, b1) = (b.a.y.min(b.b.y), b.a.y.max(b.b.y));
    a0 - tol <= b1 && b0 - tol <= a1
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
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!("'{key}' must be a boolean"))),
        },
        Some(Value::Number(n)) => Ok(Some(n.as_f64().unwrap_or(0.0) != 0.0)),
        Some(_) => Err(ToolError::Validation(format!("'{key}' must be a boolean"))),
    }
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn lines(ls: &[Vec<(f64, f64)>]) -> String {
        let mut l = Layer::new("l")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        for coords in ls {
            let cs: Vec<Coord> = coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
            l.add_feature(Some(Geometry::LineString(cs)), &[]).unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = FeatureToLineTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// A plus-shape (two crossing lines) nodes into 4 segments.
    #[test]
    fn crossing_lines_split_into_four() {
        // Horizontal (-1,0)-(1,0) and vertical (0,-1)-(0,1) cross at origin.
        let input = lines(&[vec![(-1.0, 0.0), (1.0, 0.0)], vec![(0.0, -1.0), (0.0, 1.0)]]);
        let (out, _l) = run(json!({ "input": input }));
        assert_eq!(out.outputs["output_lines"], json!(4));
    }

    /// A T-junction: the vertical touches the horizontal's interior -> 3 pieces.
    #[test]
    fn t_junction_splits_the_crossed_line() {
        let input = lines(&[vec![(-1.0, 0.0), (1.0, 0.0)], vec![(0.0, 0.0), (0.0, 1.0)]]);
        let (out, _l) = run(json!({ "input": input }));
        // Horizontal splits into 2, vertical stays 1 -> 3 total.
        assert_eq!(out.outputs["output_lines"], json!(3));
    }

    /// Non-intersecting distant lines pass through unchanged.
    #[test]
    fn disjoint_lines_unchanged() {
        let input = lines(&[vec![(0.0, 0.0), (1.0, 0.0)], vec![(0.0, 10.0), (1.0, 10.0)]]);
        let (out, _l) = run(json!({ "input": input }));
        assert_eq!(out.outputs["output_lines"], json!(2));
    }

    /// A multi-vertex line yields one segment per span.
    #[test]
    fn polyline_becomes_per_span_segments() {
        let input = lines(&[vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]]);
        let (out, _l) = run(json!({ "input": input }));
        assert_eq!(out.outputs["output_lines"], json!(2));
    }

    #[test]
    fn rejects_missing_input() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            FeatureToLineTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "l.shp" })).is_ok());
    }
}
