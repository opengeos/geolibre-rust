//! GeoLibre tool: extract selected feature vertices as a point layer.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Feature Vertices To Points* (Data
//! Management). The bundled Whitebox `extract_nodes` only covers the ALL-vertices
//! case; the selective locations (start / end / both-ends / midpoint / dangle)
//! are not extractable as a point layer anywhere else, and dangle nodes are only
//! *reported* inside topology validation rather than materialised as points.
//!
//! Each input feature is decomposed into "parts" (one per LineString, per polygon
//! ring, or per point); the requested `point_location` selector then emits points:
//!
//! * `ALL`       — every vertex of every part.
//! * `START`     — the first vertex of each part.
//! * `END`       — the last vertex of each part.
//! * `BOTH_ENDS` — the first and last vertex of each part.
//! * `MID`       — the geometric midpoint (by arc length) of each part's boundary;
//!   this need not fall on a vertex.
//! * `DANGLE`    — line endpoints that are *not* shared with any other line
//!   endpoint in the dataset (node degree 1), computed by endpoint-degree
//!   counting — the same idea topology validation uses to flag dangles.
//!
//! Every output point carries a parent-feature id (`orig_fid`) so the points can
//! be joined back to their source feature. All original attributes are copied
//! through, mirroring ArcGIS behaviour.

use std::collections::{BTreeMap, HashMap, HashSet};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct FeatureVerticesToPointsTool;

impl Tool for FeatureVerticesToPointsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "feature_vertices_to_points",
            display_name: "Feature Vertices To Points",
            summary: "Create a point layer from selected feature vertices — ALL, START, END, BOTH_ENDS, MID (arc-length midpoint), or DANGLE (unshared line endpoints) — with a parent-FID field, like ArcGIS Feature Vertices To Points.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input line, polygon, or point vector layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output point vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "point_location",
                    description: "Which vertices to extract: ALL | START | END | BOTH_ENDS | MID | DANGLE. Default ALL.",
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
        parse_location(args)?;
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
        let location = parse_location(args)?;

        let layer = load_input_layer(input)?;

        // Output schema: copy every input field, then a (uniquely-named) parent id.
        let mut out = Layer::new("vertices").with_geom_type(GeometryType::Point);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for def in layer.schema.fields() {
            out.add_field(def.clone());
        }
        let fid_field = unique_field_name("orig_fid", layer.schema.fields());
        out.add_field(FieldDef::new(fid_field.as_str(), FieldType::Integer));

        // DANGLE needs a global node-degree map built over every open line part.
        let dangles: Option<HashSet<NodeKey>> = if location == PointLocation::Dangle {
            Some(dangle_nodes(&layer))
        } else {
            None
        };

        let mut point_count = 0usize;
        for feature in layer.features.iter() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let selected = select_points(geom, location, dangles.as_ref());
            for c in selected {
                let mut attrs = feature.attributes.clone();
                // Guard against a schema/feature length mismatch from odd inputs.
                attrs.resize(layer.schema.len(), FieldValue::Null);
                attrs.push(FieldValue::Integer(feature.fid as i64));
                out.push(wbvector::Feature {
                    fid: point_count as u64,
                    geometry: Some(Geometry::Point(c)),
                    attributes: attrs,
                });
                point_count += 1;
            }
        }

        ctx.progress.info(&format!(
            "extracted {point_count} point(s) using location '{}'",
            location.as_str()
        ));

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("point_count".to_string(), json!(point_count));
        outputs.insert(
            "input_feature_count".to_string(),
            json!(layer.features.len()),
        );
        outputs.insert("point_location".to_string(), json!(location.as_str()));
        Ok(ToolRunResult { outputs })
    }
}

/// Selects the output coordinates for one geometry under `location`.
fn select_points(
    geom: &Geometry,
    location: PointLocation,
    dangles: Option<&HashSet<NodeKey>>,
) -> Vec<Coord> {
    let mut out = Vec::new();
    for part in parts(geom) {
        let cs = &part.coords;
        if cs.is_empty() {
            continue;
        }
        match location {
            PointLocation::All => out.extend(cs.iter().cloned()),
            PointLocation::Start => out.push(cs[0].clone()),
            PointLocation::End => out.push(cs[cs.len() - 1].clone()),
            PointLocation::BothEnds => {
                out.push(cs[0].clone());
                if cs.len() > 1 {
                    out.push(cs[cs.len() - 1].clone());
                }
            }
            PointLocation::Mid => {
                if let Some(mid) = arc_length_midpoint(cs, part.closed) {
                    out.push(mid);
                }
            }
            PointLocation::Dangle => {
                // Only open line parts with real endpoints can dangle.
                if part.closed || cs.len() < 2 {
                    continue;
                }
                let set = dangles.expect("dangle set present in DANGLE mode");
                for c in [&cs[0], &cs[cs.len() - 1]] {
                    if set.contains(&node_key(c)) {
                        out.push(c.clone());
                    }
                }
            }
        }
    }
    out
}

/// A part is a single coordinate chain plus whether it is a closed ring.
struct Part {
    coords: Vec<Coord>,
    closed: bool,
}

/// Decomposes a geometry into coordinate chains. Polygon rings drop their closing
/// duplicate vertex so each ring is a sequence of distinct vertices; consecutive
/// duplicate coordinates in line chains are collapsed.
fn parts(geom: &Geometry) -> Vec<Part> {
    let line = |cs: &[Coord]| Part {
        coords: dedup_consecutive(cs),
        closed: false,
    };
    let ring = |cs: &[Coord]| Part {
        coords: open_ring(cs),
        closed: true,
    };
    match geom {
        Geometry::Point(c) => vec![Part {
            coords: vec![c.clone()],
            closed: false,
        }],
        Geometry::MultiPoint(cs) => cs
            .iter()
            .map(|c| Part {
                coords: vec![c.clone()],
                closed: false,
            })
            .collect(),
        Geometry::LineString(cs) => vec![line(cs)],
        Geometry::MultiLineString(ls) => ls.iter().map(|l| line(l)).collect(),
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            let mut v = vec![ring(exterior.coords())];
            v.extend(interiors.iter().map(|r| ring(r.coords())));
            v
        }
        Geometry::MultiPolygon(polys) => {
            let mut v = Vec::new();
            for (ext, holes) in polys {
                v.push(ring(ext.coords()));
                v.extend(holes.iter().map(|r| ring(r.coords())));
            }
            v
        }
        Geometry::GeometryCollection(geoms) => geoms.iter().flat_map(parts).collect(),
    }
}

/// Collapses runs of coincident consecutive coordinates.
fn dedup_consecutive(cs: &[Coord]) -> Vec<Coord> {
    let mut out: Vec<Coord> = Vec::with_capacity(cs.len());
    for c in cs {
        if out.last().is_none_or(|l| !same_xy(l, c)) {
            out.push(c.clone());
        }
    }
    out
}

/// Returns a ring's distinct vertices (drops the closing duplicate).
fn open_ring(cs: &[Coord]) -> Vec<Coord> {
    let mut v = dedup_consecutive(cs);
    if v.len() >= 2 && same_xy(&v[0], &v[v.len() - 1]) {
        v.pop();
    }
    v
}

/// Geometric midpoint of a chain by arc length. For a closed ring the traversed
/// boundary includes the closing segment back to the first vertex.
fn arc_length_midpoint(cs: &[Coord], closed: bool) -> Option<Coord> {
    if cs.is_empty() {
        return None;
    }
    if cs.len() == 1 {
        return Some(cs[0].clone());
    }
    // Build the walk (append the closing vertex for a ring).
    let mut walk: Vec<&Coord> = cs.iter().collect();
    if closed {
        walk.push(&cs[0]);
    }
    let mut total = 0.0f64;
    for w in walk.windows(2) {
        total += dist(w[0], w[1]);
    }
    if total <= 0.0 {
        return Some(cs[0].clone());
    }
    let target = total / 2.0;
    let mut acc = 0.0f64;
    for w in walk.windows(2) {
        let seg = dist(w[0], w[1]);
        if acc + seg >= target {
            let t = if seg > 0.0 { (target - acc) / seg } else { 0.0 };
            return Some(Coord::xy(
                w[0].x + (w[1].x - w[0].x) * t,
                w[0].y + (w[1].y - w[0].y) * t,
            ));
        }
        acc += seg;
    }
    Some(cs[cs.len() - 1].clone())
}

/// Builds the set of node keys that are dangles (endpoint-degree exactly 1) over
/// every open line part in the layer.
fn dangle_nodes(layer: &Layer) -> HashSet<NodeKey> {
    let mut degree: HashMap<NodeKey, u32> = HashMap::new();
    for feature in layer.features.iter() {
        let Some(geom) = feature.geometry.as_ref() else {
            continue;
        };
        for part in parts(geom) {
            if part.closed || part.coords.len() < 2 {
                continue;
            }
            let cs = &part.coords;
            *degree.entry(node_key(&cs[0])).or_insert(0) += 1;
            *degree.entry(node_key(&cs[cs.len() - 1])).or_insert(0) += 1;
        }
    }
    degree
        .into_iter()
        .filter(|&(_, d)| d == 1)
        .map(|(k, _)| k)
        .collect()
}

// ── Geometry helpers ─────────────────────────────────────────────────────────

/// Quantised coordinate key used for endpoint coincidence (dangle) tests.
type NodeKey = (i64, i64);

/// Snap scale for node coincidence: ~1e-6 units (0.1 m at UTM scale, ~0.1 m at
/// degree scale). Cleanly shared endpoints are exactly equal and match at any
/// scale; this tolerates minor floating-point noise.
const NODE_SCALE: f64 = 1.0e6;

fn node_key(c: &Coord) -> NodeKey {
    (
        (c.x * NODE_SCALE).round() as i64,
        (c.y * NODE_SCALE).round() as i64,
    )
}

fn same_xy(a: &Coord, b: &Coord) -> bool {
    (a.x - b.x).abs() <= f64::EPSILON && (a.y - b.y).abs() <= f64::EPSILON
}

fn dist(a: &Coord, b: &Coord) -> f64 {
    (a.x - b.x).hypot(a.y - b.y)
}

/// Picks a field name not already present in `existing` by appending underscores.
fn unique_field_name(base: &str, existing: &[FieldDef]) -> String {
    let taken: HashSet<&str> = existing.iter().map(|d| d.name.as_str()).collect();
    let mut name = base.to_string();
    while taken.contains(name.as_str()) {
        name.push('_');
    }
    name
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum PointLocation {
    All,
    Start,
    End,
    BothEnds,
    Mid,
    Dangle,
}

impl PointLocation {
    fn as_str(self) -> &'static str {
        match self {
            Self::All => "ALL",
            Self::Start => "START",
            Self::End => "END",
            Self::BothEnds => "BOTH_ENDS",
            Self::Mid => "MID",
            Self::Dangle => "DANGLE",
        }
    }
}

fn parse_location(args: &ToolArgs) -> Result<PointLocation, ToolError> {
    let raw = match parse_optional_str(args, "point_location")? {
        None => return Ok(PointLocation::All),
        Some(s) => s,
    };
    match raw.trim().to_ascii_uppercase().as_str() {
        "ALL" => Ok(PointLocation::All),
        "START" => Ok(PointLocation::Start),
        "END" => Ok(PointLocation::End),
        "BOTH_ENDS" | "BOTH" => Ok(PointLocation::BothEnds),
        "MID" | "MIDPOINT" => Ok(PointLocation::Mid),
        "DANGLE" | "DANGLES" => Ok(PointLocation::Dangle),
        other => Err(ToolError::Validation(format!(
            "parameter 'point_location' must be one of ALL, START, END, BOTH_ENDS, MID, DANGLE (got '{other}')"
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

    fn line_layer(lines: &[Vec<(f64, f64)>]) -> String {
        let mut l = Layer::new("lines")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        for (i, coords) in lines.iter().enumerate() {
            let cs = coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
            l.add_feature(
                Some(Geometry::line_string(cs)),
                &[("name", FieldValue::Text(format!("l{i}")))],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn poly_layer(exterior: &[(f64, f64)]) -> String {
        let mut l = Layer::new("polys")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        let cs = exterior.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
        l.add_feature(Some(Geometry::polygon(cs, vec![])), &[])
            .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = FeatureVerticesToPointsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn pts(layer: &Layer) -> Vec<(f64, f64)> {
        layer
            .features
            .iter()
            .map(|f| match f.geometry.as_ref().unwrap() {
                Geometry::Point(c) => (c.x, c.y),
                other => panic!("expected point, got {other:?}"),
            })
            .collect()
    }

    /// ALL emits every vertex and copies the parent id + attributes.
    #[test]
    fn all_emits_every_vertex() {
        let input = line_layer(&[vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]]);
        let (out, layer) = run(json!({ "input": input, "point_location": "ALL" }));
        assert_eq!(out.outputs["point_count"], json!(3));
        assert_eq!(pts(&layer), vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]);
        // orig_fid + copied "name" attribute present.
        let fi = layer.schema.field_index("orig_fid").unwrap();
        let ni = layer.schema.field_index("name").unwrap();
        assert_eq!(layer.features[0].attributes[fi], FieldValue::Integer(0));
        assert_eq!(
            layer.features[0].attributes[ni],
            FieldValue::Text("l0".into())
        );
    }

    /// START/END/BOTH_ENDS pick the right endpoints (per part).
    #[test]
    fn endpoints_selection() {
        let input = line_layer(&[vec![(0.0, 0.0), (5.0, 0.0), (10.0, 0.0)]]);
        let (_o, s) = run(json!({ "input": &input, "point_location": "START" }));
        assert_eq!(pts(&s), vec![(0.0, 0.0)]);
        let (_o, e) = run(json!({ "input": &input, "point_location": "END" }));
        assert_eq!(pts(&e), vec![(10.0, 0.0)]);
        let (_o, b) = run(json!({ "input": &input, "point_location": "BOTH_ENDS" }));
        assert_eq!(pts(&b), vec![(0.0, 0.0), (10.0, 0.0)]);
    }

    /// MID is the arc-length midpoint, not necessarily a vertex.
    #[test]
    fn mid_is_arclength_midpoint() {
        // Total length 10; midpoint at x=5 which is NOT a vertex.
        let input = line_layer(&[vec![(0.0, 0.0), (2.0, 0.0), (10.0, 0.0)]]);
        let (out, layer) = run(json!({ "input": input, "point_location": "MID" }));
        assert_eq!(out.outputs["point_count"], json!(1));
        let (x, y) = pts(&layer)[0];
        assert!((x - 5.0).abs() < 1e-9 && y.abs() < 1e-9, "mid at ({x},{y})");
    }

    /// Polygon ALL drops the closing duplicate vertex; MID walks the full ring.
    #[test]
    fn polygon_ring_vertices() {
        // Closed square (5 coords incl. closing dup) -> 4 distinct vertices.
        let input = poly_layer(&[
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 10.0),
            (0.0, 10.0),
            (0.0, 0.0),
        ]);
        let (out, _l) = run(json!({ "input": input, "point_location": "ALL" }));
        assert_eq!(out.outputs["point_count"], json!(4));
    }

    /// DANGLE keeps only endpoints not shared with another line's endpoint.
    #[test]
    fn dangle_keeps_unshared_endpoints() {
        // Two lines meeting at (5,0): that node is shared (degree 2, not a dangle).
        // Free ends at (0,0) and (10,0) are dangles.
        let input = line_layer(&[vec![(0.0, 0.0), (5.0, 0.0)], vec![(5.0, 0.0), (10.0, 0.0)]]);
        let (out, layer) = run(json!({ "input": input, "point_location": "DANGLE" }));
        assert_eq!(out.outputs["point_count"], json!(2));
        let mut got = pts(&layer);
        got.sort_by(|a, b| a.0.total_cmp(&b.0));
        assert_eq!(got, vec![(0.0, 0.0), (10.0, 0.0)]);
    }

    /// Non-line geometries pass through DANGLE producing no points.
    #[test]
    fn dangle_ignores_polygons() {
        let input = poly_layer(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 0.0)]);
        let (out, _l) = run(json!({ "input": input, "point_location": "DANGLE" }));
        assert_eq!(out.outputs["point_count"], json!(0));
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            FeatureVerticesToPointsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err()); // no input
        assert!(bad(json!({ "input": "a.geojson", "point_location": "NOPE" })).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_ok()); // default ALL
        assert!(bad(json!({ "input": "a.geojson", "point_location": "dangle" })).is_ok());
    }
}
