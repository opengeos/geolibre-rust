//! GeoLibre tool: split polylines at snapped point locations.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Split Line At Point* (Data
//! Management). The bundled Whitebox suite can split lines at their mutual
//! intersections (`split_with_lines`) or by a maximum segment length
//! (`split_vector_lines`), but none split a polyline at arbitrary point
//! locations — the standard operation for inserting nodes at gauges, junctions,
//! access points, or sampling stations.
//!
//! For each input line and each point that falls within `search_radius` of the
//! line, the point is projected onto the nearest segment; where the
//! perpendicular projection lands strictly inside the line (not on an endpoint)
//! a split node is inserted at the projected location. Each resulting
//! sub-segment is emitted as its own `LineString` feature carrying the parent
//! line's attributes plus a `src_fid` (source feature index) and `seg_index`
//! (0-based position of the sub-segment along the parent). Lines with no point
//! within the radius pass through unchanged as a single feature.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct SplitLineAtPointTool;

impl Tool for SplitLineAtPointTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "split_line_at_point",
            display_name: "Split Line At Point",
            summary: "Split polyline features at the locations of point features that fall within a search radius, inserting a node at each snapped point and emitting one feature per resulting sub-segment (parent attributes preserved) — like ArcGIS Split Line At Point.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input_lines",
                    description: "Input polyline vector layer to be split.",
                    required: true,
                },
                ToolParamSpec {
                    name: "point_features",
                    description: "Point (or multipoint) vector layer whose locations define the split nodes.",
                    required: true,
                },
                ToolParamSpec {
                    name: "search_radius",
                    description: "Maximum distance, in the lines' CRS units, between a point and a line for the point to split that line. Required, positive.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output line vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["input_lines", "point_features"] {
            if args
                .get(key)
                .and_then(Value::as_str)
                .map(str::trim)
                .unwrap_or("")
                .is_empty()
            {
                return Err(ToolError::Validation(format!(
                    "missing required string parameter '{key}'"
                )));
            }
        }
        parse_radius(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input_lines = require_str(args, "input_lines")?;
        let point_features = require_str(args, "point_features")?;
        let output = parse_optional_str(args, "output")?;
        let radius = parse_radius(args)?;

        let lines = load_input_layer(input_lines)?;
        let points_layer = load_input_layer(point_features)?;

        // Collect all split points from the point layer.
        let mut pts: Vec<P> = Vec::new();
        for f in lines_points(&points_layer) {
            pts.push(f);
        }
        ctx.progress
            .info(&format!("loaded {} split point(s)", pts.len()));

        // Build the output layer: preserve the parent schema, then append the
        // two provenance fields (with collision-safe names).
        let mut out = Layer::new("split_lines").with_geom_type(GeometryType::LineString);
        if let Some(epsg) = lines.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for def in lines.schema.fields() {
            out.add_field(def.clone());
        }
        let src_fid_name = unique_field_name(&out, "src_fid");
        let seg_index_name = unique_field_name_with(&out, "seg_index", &src_fid_name);
        out.add_field(FieldDef::new(&src_fid_name, FieldType::Integer));
        out.add_field(FieldDef::new(&seg_index_name, FieldType::Integer));

        let mut input_line_count = 0usize;
        let mut split_count = 0usize; // number of split nodes actually inserted
        for (fidx, feature) in lines.features.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let chains = line_chains(geom);
            if chains.is_empty() {
                continue;
            }
            input_line_count += 1;
            let mut seg_index = 0i64;
            for chain in chains {
                if chain.len() < 2 {
                    continue;
                }
                let cuts = split_positions(&chain, &pts, radius);
                split_count += cuts.len();
                for sub in split_chain(&chain, &cuts) {
                    if sub.len() < 2 {
                        continue;
                    }
                    let cs: Vec<Coord> = sub.iter().map(|p| Coord::xy(p.x, p.y)).collect();
                    let mut attrs = feature.attributes.clone();
                    // Pad/truncate defensively so length matches the parent field
                    // count regardless of source layer quirks.
                    attrs.resize(lines.schema.len(), wbvector::FieldValue::Null);
                    attrs.push(wbvector::FieldValue::Integer(fidx as i64));
                    attrs.push(wbvector::FieldValue::Integer(seg_index));
                    out.push(Feature {
                        fid: out.len() as u64,
                        geometry: Some(Geometry::line_string(cs)),
                        attributes: attrs,
                    });
                    seg_index += 1;
                }
            }
        }

        let output_feature_count = out.len();
        ctx.progress.info(&format!(
            "split {input_line_count} line(s) at {split_count} node(s) into {output_feature_count} feature(s)"
        ));

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_line_count".to_string(), json!(input_line_count));
        outputs.insert("point_count".to_string(), json!(pts.len()));
        outputs.insert("split_count".to_string(), json!(split_count));
        outputs.insert(
            "output_feature_count".to_string(),
            json!(output_feature_count),
        );
        Ok(ToolRunResult { outputs })
    }
}

// ── Core geometry ─────────────────────────────────────────────────────────────

const EPS: f64 = 1e-9;

#[derive(Clone, Copy, Debug)]
struct P {
    x: f64,
    y: f64,
}

fn dist(a: P, b: P) -> f64 {
    (a.x - b.x).hypot(a.y - b.y)
}

fn interp(a: P, b: P, t: f64) -> P {
    P {
        x: a.x + (b.x - a.x) * t,
        y: a.y + (b.y - a.y) * t,
    }
}

/// Arc-length positions along `chain` at which points in `pts` (within
/// `radius`) should split the line. Positions are strictly interior
/// (0 < s < total), sorted ascending and de-duplicated.
fn split_positions(chain: &[P], pts: &[P], radius: f64) -> Vec<f64> {
    // Cumulative arc-length at each vertex.
    let mut cum = vec![0.0f64; chain.len()];
    for i in 1..chain.len() {
        cum[i] = cum[i - 1] + dist(chain[i - 1], chain[i]);
    }
    let total = *cum.last().unwrap();
    if total <= EPS {
        return Vec::new();
    }

    let mut positions: Vec<f64> = Vec::new();
    for &pt in pts {
        // Project the point onto every segment; keep the nearest.
        let mut best_dist = f64::INFINITY;
        let mut best_arc = 0.0f64;
        for i in 1..chain.len() {
            let a = chain[i - 1];
            let b = chain[i];
            let seg_len = dist(a, b);
            if seg_len <= EPS {
                continue;
            }
            let dx = b.x - a.x;
            let dy = b.y - a.y;
            let t = (((pt.x - a.x) * dx + (pt.y - a.y) * dy) / (seg_len * seg_len)).clamp(0.0, 1.0);
            let proj = interp(a, b, t);
            let d = dist(pt, proj);
            if d < best_dist {
                best_dist = d;
                best_arc = cum[i - 1] + t * seg_len;
            }
        }
        if best_dist <= radius && best_arc > EPS && best_arc < total - EPS {
            positions.push(best_arc);
        }
    }

    positions.sort_by(f64::total_cmp);
    positions.dedup_by(|a, b| (*a - *b).abs() <= EPS.max(total * 1e-9));
    positions
}

/// Splits a polyline into sub-chains at the given (sorted, interior) arc-length
/// positions, inserting the interpolated node at each cut. A cut that lands on
/// an existing vertex simply splits there without duplicating the vertex.
fn split_chain(chain: &[P], cuts: &[f64]) -> Vec<Vec<P>> {
    if cuts.is_empty() {
        return vec![chain.to_vec()];
    }
    let mut cum = vec![0.0f64; chain.len()];
    for i in 1..chain.len() {
        cum[i] = cum[i - 1] + dist(chain[i - 1], chain[i]);
    }

    let mut result: Vec<Vec<P>> = Vec::new();
    let mut current: Vec<P> = vec![chain[0]];
    let mut ci = 0usize;

    for i in 1..chain.len() {
        let (a0, a1) = (cum[i - 1], cum[i]);
        let seg_len = (a1 - a0).max(EPS);
        while ci < cuts.len() && cuts[ci] < a1 - EPS {
            let c = cuts[ci];
            let t = ((c - a0) / seg_len).clamp(0.0, 1.0);
            if t <= EPS {
                // Cut on the start vertex of this segment: split there only if
                // the current sub-chain has actual extent.
                if current.len() > 1 {
                    let last = *current.last().unwrap();
                    result.push(std::mem::replace(&mut current, vec![last]));
                }
            } else {
                let p = interp(chain[i - 1], chain[i], t);
                current.push(p);
                result.push(std::mem::replace(&mut current, vec![p]));
            }
            ci += 1;
        }
        current.push(chain[i]);
    }
    result.push(current);
    result
}

// ── Layer helpers ─────────────────────────────────────────────────────────────

/// Extracts every vertex of point/multipoint features as split candidates.
fn lines_points(layer: &Layer) -> Vec<P> {
    let mut out = Vec::new();
    for f in layer.features.iter() {
        let Some(geom) = f.geometry.as_ref() else {
            continue;
        };
        match geom {
            Geometry::Point(c) => out.push(P { x: c.x, y: c.y }),
            Geometry::MultiPoint(cs) => {
                for c in cs {
                    out.push(P { x: c.x, y: c.y });
                }
            }
            // Tolerate stray non-point geometries by using their vertices.
            other => {
                for c in other.all_coords() {
                    out.push(P { x: c.x, y: c.y });
                }
            }
        }
    }
    out
}

/// Splits a line geometry into one vertex-list per part, dropping consecutive
/// duplicate vertices.
fn line_chains(geom: &Geometry) -> Vec<Vec<P>> {
    let to_pts = |cs: &[Coord]| -> Vec<P> {
        let mut out: Vec<P> = Vec::with_capacity(cs.len());
        for c in cs {
            let p = P { x: c.x, y: c.y };
            if out.last().is_none_or(|l| dist(*l, p) > 1e-12) {
                out.push(p);
            }
        }
        out
    };
    match geom {
        Geometry::LineString(cs) => vec![to_pts(cs)],
        Geometry::MultiLineString(lines) => lines.iter().map(|l| to_pts(l)).collect(),
        _ => Vec::new(),
    }
}

fn unique_field_name(layer: &Layer, base: &str) -> String {
    if layer.schema.field_index(base).is_none() {
        return base.to_string();
    }
    let mut i = 1;
    loop {
        let cand = format!("{base}_{i}");
        if layer.schema.field_index(&cand).is_none() {
            return cand;
        }
        i += 1;
    }
}

fn unique_field_name_with(layer: &Layer, base: &str, taken: &str) -> String {
    let mut cand = unique_field_name(layer, base);
    while cand == taken {
        cand = format!("{cand}_");
    }
    cand
}

// ── Parameters ────────────────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

fn parse_radius(args: &ToolArgs) -> Result<f64, ToolError> {
    let radius = parse_optional_f64(args, "search_radius")?.ok_or_else(|| {
        ToolError::Validation("required parameter 'search_radius' is missing".into())
    })?;
    if !(radius > 0.0 && radius.is_finite()) {
        return Err(ToolError::Validation(
            "'search_radius' must be a positive number".to_string(),
        ));
    }
    Ok(radius)
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
    use wbvector::{memory_store, FieldValue};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn line_layer(coords: &[(f64, f64)]) -> String {
        let mut l = Layer::new("lines")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        let cs = coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
        l.add_feature(Some(Geometry::line_string(cs)), &[("name", "road".into())])
            .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn point_layer(coords: &[(f64, f64)]) -> String {
        let mut l = Layer::new("points")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        for &(x, y) in coords {
            l.add_feature(Some(Geometry::point(x, y)), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SplitLineAtPointTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn total_length(layer: &Layer) -> f64 {
        let mut sum = 0.0;
        for f in layer.iter() {
            if let Some(Geometry::LineString(cs)) = f.geometry.as_ref() {
                for w in cs.windows(2) {
                    sum += (w[0].x - w[1].x).hypot(w[0].y - w[1].y);
                }
            }
        }
        sum
    }

    /// A point on a horizontal line splits it into two sub-segments whose
    /// combined length equals the original (length conservation), and the split
    /// node lands at the projected x.
    #[test]
    fn splits_at_snapped_point() {
        let line = line_layer(&[(0.0, 0.0), (100.0, 0.0)]);
        let points = point_layer(&[(40.0, 5.0)]); // 5 units off the line
        let (out, layer) = run(json!({
            "input_lines": line, "point_features": points, "search_radius": 10.0,
        }));
        assert_eq!(out.outputs["split_count"], json!(1));
        assert_eq!(out.outputs["output_feature_count"], json!(2));
        // Length conserved.
        assert!((total_length(&layer) - 100.0).abs() < 1e-6);
        // The shared node sits at x=40 (perpendicular projection).
        let mut endpoints: Vec<f64> = Vec::new();
        for f in layer.iter() {
            if let Some(Geometry::LineString(cs)) = f.geometry.as_ref() {
                endpoints.push(cs[0].x);
                endpoints.push(cs.last().unwrap().x);
            }
        }
        assert!(endpoints.iter().any(|&x| (x - 40.0).abs() < 1e-6));
    }

    /// Parent attributes are carried to every sub-segment.
    #[test]
    fn preserves_parent_attributes() {
        let line = line_layer(&[(0.0, 0.0), (100.0, 0.0)]);
        let points = point_layer(&[(30.0, 0.0), (60.0, 0.0)]);
        let (_o, layer) = run(json!({
            "input_lines": line, "point_features": points, "search_radius": 1.0,
        }));
        assert_eq!(layer.len(), 3); // two cuts -> three segments
        let ni = layer.schema.field_index("name").unwrap();
        for f in layer.iter() {
            assert_eq!(f.attributes[ni], FieldValue::Text("road".into()));
        }
        // seg_index runs 0,1,2.
        let si = layer.schema.field_index("seg_index").unwrap();
        let mut idxs: Vec<i64> = layer
            .iter()
            .map(|f| f.attributes[si].as_i64().unwrap())
            .collect();
        idxs.sort();
        assert_eq!(idxs, vec![0, 1, 2]);
    }

    /// A point beyond the search radius leaves the line untouched (pass-through).
    #[test]
    fn non_matching_point_passes_through() {
        let line = line_layer(&[(0.0, 0.0), (100.0, 0.0)]);
        let points = point_layer(&[(50.0, 50.0)]); // 50 units away
        let (out, layer) = run(json!({
            "input_lines": line, "point_features": points, "search_radius": 10.0,
        }));
        assert_eq!(out.outputs["split_count"], json!(0));
        assert_eq!(layer.len(), 1);
        assert!((total_length(&layer) - 100.0).abs() < 1e-6);
    }

    /// A point projecting onto a line endpoint does not create a zero-length
    /// segment (no split at endpoints).
    #[test]
    fn no_split_at_endpoints() {
        let line = line_layer(&[(0.0, 0.0), (100.0, 0.0)]);
        let points = point_layer(&[(0.0, 2.0), (100.0, 2.0)]);
        let (out, layer) = run(json!({
            "input_lines": line, "point_features": points, "search_radius": 10.0,
        }));
        assert_eq!(out.outputs["split_count"], json!(0));
        assert_eq!(layer.len(), 1);
    }

    /// Two points near the same location collapse to a single split.
    #[test]
    fn coincident_points_dedupe() {
        let line = line_layer(&[(0.0, 0.0), (100.0, 0.0)]);
        let points = point_layer(&[(50.0, 1.0), (50.0, -1.0)]);
        let (out, _layer) = run(json!({
            "input_lines": line, "point_features": points, "search_radius": 5.0,
        }));
        assert_eq!(out.outputs["output_feature_count"], json!(2));
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SplitLineAtPointTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input_lines": "a.geojson" })).is_err()); // no points
        assert!(bad(json!({ "input_lines": "a.geojson", "point_features": "p.geojson" })).is_err()); // no radius
        assert!(bad(
            json!({ "input_lines": "a.geojson", "point_features": "p.geojson", "search_radius": 0 })
        )
        .is_err());
        assert!(bad(
            json!({ "input_lines": "a.geojson", "point_features": "p.geojson", "search_radius": 5 })
        )
        .is_ok());
        // String-encoded radius accepted (host UIs post strings).
        assert!(bad(json!({
            "input_lines": "a.geojson", "point_features": "p.geojson", "search_radius": "5"
        }))
        .is_ok());
    }
}
