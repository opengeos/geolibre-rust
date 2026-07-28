//! GeoLibre tool: true 3D nearest-neighbour proximity.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Near 3D* (3D Analyst). The bundled
//! `near` and the shipped `generate_near_table` both measure proximity in plan
//! view, which produces confidently wrong answers on 3D data: two utility lines
//! crossing at the same map location but 8 m apart vertically are reported as
//! 0 m apart. That inverts the result for exactly the clearance and
//! conflict-detection questions 3D proximity is used to answer.
//!
//! Distances here are straight-line 3D Euclidean and are refined exactly
//! against every candidate segment. Endpoint-only spatial pruning is unsafe
//! for long segments whose interior crosses the search radius.
//!
//! Vertices with no Z are treated as `z = 0` so a 2D near-layer degrades to a
//! planar answer rather than failing.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// A 3D point.
type P3 = [f64; 3];

pub struct Near3dTool;

impl Tool for Near3dTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "near_3d",
            display_name: "Near 3D",
            summary: "Find the nearest feature in true 3D distance, reporting distance, the nearest location, angles and per-axis deltas, like ArcGIS Near 3D.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input feature layer (3D vertices; missing Z is treated as 0).",
                    required: true,
                },
                ToolParamSpec {
                    name: "near_features",
                    description: "Comma-separated list of layers to search against. Omit to self-join the input layer.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "search_radius",
                    description: "Only consider near features within this 3D distance (map units). Unlimited when omitted.",
                    required: false,
                },
                ToolParamSpec {
                    name: "location",
                    description: "Append the nearest point's coordinates as near_x / near_y / near_z (default false).",
                    required: false,
                },
                ToolParamSpec {
                    name: "angle",
                    description: "Append the horizontal bearing (near_bearing, degrees from north) and vertical angle (near_vert_angle, degrees from horizontal) to the nearest point (default false).",
                    required: false,
                },
                ToolParamSpec {
                    name: "delta",
                    description: "Append per-axis offsets near_dx / near_dy / near_dz (default false).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        if let Some(r) = parse_optional_f64(args, "search_radius")? {
            if !r.is_finite() || r <= 0.0 {
                return Err(ToolError::Validation(
                    "'search_radius' must be greater than 0".to_string(),
                ));
            }
        }
        for k in ["location", "angle", "delta"] {
            parse_optional_bool(args, k)?;
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let radius = parse_optional_f64(args, "search_radius")?;
        let want_location = parse_optional_bool(args, "location")?.unwrap_or(false);
        let want_angle = parse_optional_bool(args, "angle")?.unwrap_or(false);
        let want_delta = parse_optional_bool(args, "delta")?.unwrap_or(false);

        let layer = load_input_layer(input)?;

        // Near layers: explicit list, or the input itself (self-join).
        let near_paths: Vec<String> = match parse_optional_str(args, "near_features")? {
            Some(s) => s
                .split(',')
                .map(str::trim)
                .filter(|t| !t.is_empty())
                .map(String::from)
                .collect(),
            None => Vec::new(),
        };
        let self_join = near_paths.is_empty();

        // Flatten every near feature into segments, tagged with (layer, feature).
        let mut near_segs: Vec<Seg> = Vec::new();
        if self_join {
            collect_segments(&layer, 0, &mut near_segs);
        } else {
            for (li, p) in near_paths.iter().enumerate() {
                let nl = load_input_layer(p)?;
                collect_segments(&nl, li, &mut near_segs);
            }
        }

        ctx.progress.info(&format!(
            "matching {} feature(s) against {} near segment(s) in 3D",
            layer.len(),
            near_segs.len()
        ));

        let mut out = Layer::new(layer.name.clone());
        out.crs = layer.crs.clone();
        out.geom_type = layer.geom_type;
        for fd in layer.schema.fields().iter() {
            out.add_field(fd.clone());
        }
        out.add_field(FieldDef::new("near_fid", FieldType::Integer));
        out.add_field(FieldDef::new("near_layer", FieldType::Integer));
        out.add_field(FieldDef::new("near_dist", FieldType::Float));
        if want_location {
            for f in ["near_x", "near_y", "near_z"] {
                out.add_field(FieldDef::new(f, FieldType::Float));
            }
        }
        if want_angle {
            for f in ["near_bearing", "near_vert_angle"] {
                out.add_field(FieldDef::new(f, FieldType::Float));
            }
        }
        if want_delta {
            for f in ["near_dx", "near_dy", "near_dz"] {
                out.add_field(FieldDef::new(f, FieldType::Float));
            }
        }

        let radius2 = radius.map(|r| r * r);
        let mut matched = 0usize;

        for (fi, feat) in layer.features.iter().enumerate() {
            let mut src_pts: Vec<P3> = Vec::new();
            if let Some(g) = feat.geometry.as_ref() {
                input_vertices_3d(g, &mut src_pts);
            }
            // `src` is the input vertex that ends up closest; it anchors the
            // reported deltas and angles.
            let mut src: Option<P3> = None;

            // Minimise over every input vertex and near segment. Endpoint-only
            // kd-tree candidates are not sufficient: a long segment can cross
            // the search sphere while both endpoints lie far outside it.
            let mut best: Option<(f64, P3, usize, usize)> = None;
            for p in &src_pts {
                for s in &near_segs {
                    if self_join && s.feature == fi {
                        continue;
                    }
                    let (d2, q) = point_seg_dist2(*p, s.a, s.b);
                    if radius2.is_some_and(|r2| d2 > r2) {
                        continue;
                    }
                    if best.is_none_or(|(bd, _, _, _)| d2 < bd) {
                        best = Some((d2, q, s.feature, s.layer));
                        src = Some(*p);
                    }
                }
            }
            let best = best.map(|(d2, q, f, l)| (d2.sqrt(), q, f, l));

            let mut fields: Vec<(String, FieldValue)> = layer
                .schema
                .fields()
                .iter()
                .enumerate()
                .map(|(i, fd)| (fd.name.clone(), feat.attributes[i].clone()))
                .collect();

            match (src, best) {
                (Some(p), Some((dist, q, nfid, nlayer))) => {
                    matched += 1;
                    fields.push(("near_fid".into(), FieldValue::Integer(nfid as i64)));
                    fields.push(("near_layer".into(), FieldValue::Integer(nlayer as i64)));
                    fields.push(("near_dist".into(), FieldValue::Float(dist)));
                    if want_location {
                        fields.push(("near_x".into(), FieldValue::Float(q[0])));
                        fields.push(("near_y".into(), FieldValue::Float(q[1])));
                        fields.push(("near_z".into(), FieldValue::Float(q[2])));
                    }
                    if want_angle {
                        let (dx, dy, dz) = (q[0] - p[0], q[1] - p[1], q[2] - p[2]);
                        // Bearing: degrees clockwise from north.
                        let bearing = (dx.atan2(dy).to_degrees() + 360.0) % 360.0;
                        let horiz = (dx * dx + dy * dy).sqrt();
                        let vert = dz.atan2(horiz).to_degrees();
                        fields.push(("near_bearing".into(), FieldValue::Float(bearing)));
                        fields.push(("near_vert_angle".into(), FieldValue::Float(vert)));
                    }
                    if want_delta {
                        fields.push(("near_dx".into(), FieldValue::Float(q[0] - p[0])));
                        fields.push(("near_dy".into(), FieldValue::Float(q[1] - p[1])));
                        fields.push(("near_dz".into(), FieldValue::Float(q[2] - p[2])));
                    }
                }
                _ => {
                    // No match inside the radius (or unusable geometry): emit -1 /
                    // nulls so the row survives and the miss is explicit.
                    fields.push(("near_fid".into(), FieldValue::Integer(-1)));
                    fields.push(("near_layer".into(), FieldValue::Integer(-1)));
                    fields.push(("near_dist".into(), FieldValue::Null));
                    if want_location {
                        for f in ["near_x", "near_y", "near_z"] {
                            fields.push((f.into(), FieldValue::Null));
                        }
                    }
                    if want_angle {
                        for f in ["near_bearing", "near_vert_angle"] {
                            fields.push((f.into(), FieldValue::Null));
                        }
                    }
                    if want_delta {
                        for f in ["near_dx", "near_dy", "near_dz"] {
                            fields.push((f.into(), FieldValue::Null));
                        }
                    }
                }
            }

            let refs: Vec<(&str, FieldValue)> = fields
                .iter()
                .map(|(k, v)| (k.as_str(), v.clone()))
                .collect();
            out.add_feature(feat.geometry.clone(), &refs)
                .map_err(|e| ToolError::Execution(format!("failed writing feature: {e}")))?;
            ctx.progress
                .progress((fi as f64 + 1.0) / layer.len().max(1) as f64);
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(layer.len()));
        outputs.insert("features_matched".to_string(), json!(matched));
        outputs.insert("near_segments".to_string(), json!(near_segs.len()));
        Ok(ToolRunResult { outputs })
    }
}

/// One near-feature segment (a degenerate segment `a == b` represents a point).
struct Seg {
    a: P3,
    b: P3,
    feature: usize,
    layer: usize,
}

fn z_of(c: &Coord) -> f64 {
    c.z.unwrap_or(0.0)
}

fn collect_segments(layer: &Layer, layer_idx: usize, out: &mut Vec<Seg>) {
    for (fi, f) in layer.features.iter().enumerate() {
        if let Some(g) = f.geometry.as_ref() {
            push_geom_segments(g, fi, layer_idx, out);
        }
    }
}

fn push_geom_segments(geom: &Geometry, fi: usize, li: usize, out: &mut Vec<Seg>) {
    let mut run = |cs: &[Coord]| {
        if cs.is_empty() {
            return;
        }
        if cs.len() == 1 {
            let p = [cs[0].x, cs[0].y, z_of(&cs[0])];
            out.push(Seg {
                a: p,
                b: p,
                feature: fi,
                layer: li,
            });
            return;
        }
        for w in cs.windows(2) {
            out.push(Seg {
                a: [w[0].x, w[0].y, z_of(&w[0])],
                b: [w[1].x, w[1].y, z_of(&w[1])],
                feature: fi,
                layer: li,
            });
        }
    };
    match geom {
        Geometry::Point(c) => {
            let p = [c.x, c.y, z_of(c)];
            out.push(Seg {
                a: p,
                b: p,
                feature: fi,
                layer: li,
            });
        }
        Geometry::MultiPoint(cs) => {
            for c in cs {
                let p = [c.x, c.y, z_of(c)];
                out.push(Seg {
                    a: p,
                    b: p,
                    feature: fi,
                    layer: li,
                });
            }
        }
        Geometry::LineString(cs) => run(cs),
        Geometry::MultiLineString(parts) => parts.iter().for_each(|cs| run(cs)),
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            run_ring(exterior.coords(), &mut run);
            for r in interiors {
                run_ring(r.coords(), &mut run);
            }
        }
        Geometry::MultiPolygon(parts) => {
            for (e, hs) in parts {
                run_ring(e.coords(), &mut run);
                for r in hs {
                    run_ring(r.coords(), &mut run);
                }
            }
        }
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                push_geom_segments(g, fi, li, out);
            }
        }
    }
}

/// Rings are stored without the closing duplicate, so close them explicitly or
/// the final edge (last vertex back to first) would be missing from the index.
fn run_ring(coords: &[Coord], run: &mut impl FnMut(&[Coord])) {
    if coords.len() < 2 {
        run(coords);
        return;
    }
    let mut closed: Vec<Coord> = coords.to_vec();
    closed.push(coords[0].clone());
    run(&closed);
}

/// Collects every 3D vertex of an input feature.
///
/// The nearest distance is minimised over **all** of these, not over a single
/// representative point. Taking only the first vertex would measure from
/// wherever the geometry happens to start, which for a line or ring is
/// arbitrary — and would give the wrong answer for this module's own
/// motivating case (two utility lines crossing 8 m apart vertically, where the
/// crossing is rarely at vertex 0).
fn input_vertices_3d(geom: &Geometry, out: &mut Vec<P3>) {
    match geom {
        Geometry::Point(c) => out.push([c.x, c.y, z_of(c)]),
        Geometry::MultiPoint(cs) | Geometry::LineString(cs) => {
            out.extend(cs.iter().map(|c| [c.x, c.y, z_of(c)]))
        }
        Geometry::MultiLineString(parts) => {
            out.extend(parts.iter().flatten().map(|c| [c.x, c.y, z_of(c)]))
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            out.extend(exterior.coords().iter().map(|c| [c.x, c.y, z_of(c)]));
            for r in interiors {
                out.extend(r.coords().iter().map(|c| [c.x, c.y, z_of(c)]));
            }
        }
        Geometry::MultiPolygon(parts) => {
            for (e, hs) in parts {
                out.extend(e.coords().iter().map(|c| [c.x, c.y, z_of(c)]));
                for r in hs {
                    out.extend(r.coords().iter().map(|c| [c.x, c.y, z_of(c)]));
                }
            }
        }
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                input_vertices_3d(g, out);
            }
        }
    }
}

/// Squared 3D distance from `p` to segment `a`-`b`, plus the closest point.
/// A degenerate segment collapses to point-to-point.
fn point_seg_dist2(p: P3, a: P3, b: P3) -> (f64, P3) {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ap = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let len2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
    let t = if len2 <= f64::EPSILON {
        0.0
    } else {
        ((ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]) / len2).clamp(0.0, 1.0)
    };
    let q = [a[0] + ab[0] * t, a[1] + ab[1] * t, a[2] + ab[2] * t];
    let d = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2], q)
}

// ── parameter parsing ────────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    match args.get(key).and_then(Value::as_str) {
        Some(s) if !s.trim().is_empty() => Ok(s),
        _ => Err(ToolError::Validation(format!(
            "missing required string parameter '{key}'"
        ))),
    }
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

fn parse_optional_bool(args: &ToolArgs, k: &str) -> Result<Option<bool>, ToolError> {
    match args.get(k) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{k}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{k}' must be a boolean"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn points(items: Vec<(f64, f64, f64)>) -> String {
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        for (x, y, z) in items {
            l.add_feature(Some(Geometry::Point(Coord::xyz(x, y, z))), &[])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn lines(items: Vec<Vec<Coord>>) -> String {
        let mut l = Layer::new("l")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        for cs in items {
            l.add_feature(Some(Geometry::LineString(cs)), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = Near3dTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn num(layer: &Layer, row: usize, field: &str) -> f64 {
        let i = layer.schema.field_index(field).unwrap();
        layer.features[row].attributes[i].as_f64().unwrap()
    }

    /// THE regression: features coincident in plan view but separated in Z must
    /// report the vertical separation, not 0.
    #[test]
    fn vertical_separation_is_measured() {
        let input = points(vec![(0.0, 0.0, 0.0)]);
        let near = points(vec![(0.0, 0.0, 8.0)]);
        let (_, layer) = run(json!({ "input": input, "near_features": near }));
        assert!(
            (num(&layer, 0, "near_dist") - 8.0).abs() < 1e-9,
            "2D proximity would report 0; got {}",
            num(&layer, 0, "near_dist")
        );
    }

    /// The nearest point on a long segment lies between vertices, so a
    /// vertex-only answer would be wrong.
    #[test]
    fn refines_to_segment_interior() {
        // Segment from (0,0,0) to (100,0,0); query point above its midpoint.
        let near = lines(vec![vec![
            Coord::xyz(0.0, 0.0, 0.0),
            Coord::xyz(100.0, 0.0, 0.0),
        ]]);
        let input = points(vec![(50.0, 0.0, 3.0)]);
        let (_, layer) = run(json!({
            "input": input, "near_features": near, "location": true
        }));
        assert!(
            (num(&layer, 0, "near_dist") - 3.0).abs() < 1e-9,
            "nearest vertex is 50 units away; exact answer is 3"
        );
        assert!((num(&layer, 0, "near_x") - 50.0).abs() < 1e-9);
    }

    /// Picks the genuinely closest of several candidates in 3D.
    #[test]
    fn picks_closest_in_3d() {
        let input = points(vec![(0.0, 0.0, 0.0)]);
        // A is closer in plan (2 away) but 20 up; B is 5 away in plan, same Z.
        let near = points(vec![(2.0, 0.0, 20.0), (5.0, 0.0, 0.0)]);
        let (_, layer) = run(json!({ "input": input, "near_features": near }));
        assert_eq!(num(&layer, 0, "near_fid"), 1.0, "should choose B");
        assert!((num(&layer, 0, "near_dist") - 5.0).abs() < 1e-9);
    }

    /// The module doc's motivating case: two utility lines crossing 8 m apart
    /// vertically, where the crossing is NOT at either line's first vertex.
    /// Measuring from a single representative vertex would report the distance
    /// to an arbitrary line end instead of the clearance at the crossing.
    #[test]
    fn measures_clearance_at_the_crossing_not_the_first_vertex() {
        // Input runs west->east at z=0, crossing x=50 at its midpoint.
        let input = lines(vec![vec![
            Coord::xyz(0.0, 0.0, 0.0),
            Coord::xyz(50.0, 0.0, 0.0),
            Coord::xyz(100.0, 0.0, 0.0),
        ]]);
        // Near line runs south->north at x=50, 8 m above.
        let near = lines(vec![vec![
            Coord::xyz(50.0, -100.0, 8.0),
            Coord::xyz(50.0, 100.0, 8.0),
        ]]);
        let (_, layer) = run(json!({ "input": input, "near_features": near }));
        let d = num(&layer, 0, "near_dist");
        assert!(
            (d - 8.0).abs() < 1e-9,
            "expected the 8 m vertical clearance at the crossing, got {d} \
             (first-vertex reduction would report ~50)"
        );
    }

    /// search_radius excludes far features and the miss is explicit.
    #[test]
    fn search_radius_excludes() {
        let input = points(vec![(0.0, 0.0, 0.0)]);
        let near = points(vec![(0.0, 0.0, 100.0)]);
        let (out, layer) = run(json!({
            "input": input, "near_features": near, "search_radius": 10
        }));
        assert_eq!(num(&layer, 0, "near_fid"), -1.0);
        assert_eq!(out.outputs["features_matched"], json!(0));
    }

    /// A self-join must not match a feature to itself.
    #[test]
    fn self_join_skips_self() {
        let input = points(vec![(0.0, 0.0, 0.0), (0.0, 0.0, 4.0)]);
        let (out, layer) = run(json!({ "input": input }));
        assert_eq!(out.outputs["features_matched"], json!(2));
        assert_eq!(num(&layer, 0, "near_fid"), 1.0);
        assert_eq!(num(&layer, 1, "near_fid"), 0.0);
        assert!((num(&layer, 0, "near_dist") - 4.0).abs() < 1e-9);
    }

    /// Deltas and angles describe the offset to the nearest point.
    #[test]
    fn reports_deltas_and_angles() {
        let input = points(vec![(0.0, 0.0, 0.0)]);
        // Due north, 10 out and 10 up -> bearing 0, vertical angle 45.
        let near = points(vec![(0.0, 10.0, 10.0)]);
        let (_, layer) = run(json!({
            "input": input, "near_features": near, "delta": true, "angle": true
        }));
        assert!((num(&layer, 0, "near_dy") - 10.0).abs() < 1e-9);
        assert!((num(&layer, 0, "near_dz") - 10.0).abs() < 1e-9);
        assert!(num(&layer, 0, "near_bearing").abs() < 1e-6);
        assert!((num(&layer, 0, "near_vert_angle") - 45.0).abs() < 1e-6);
    }

    /// A 2D near-layer (no Z) degrades to a planar answer instead of failing.
    #[test]
    fn handles_2d_near_layer() {
        let mut l = Layer::new("p2")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_feature(Some(Geometry::Point(Coord::xy(3.0, 4.0))), &[])
            .unwrap();
        let id = memory_store::put_vector(l);
        let near = memory_store::make_vector_memory_path(&id);
        let input = points(vec![(0.0, 0.0, 0.0)]);
        let (_, layer) = run(json!({ "input": input, "near_features": near }));
        assert!((num(&layer, 0, "near_dist") - 5.0).abs() < 1e-9);
    }

    #[test]
    fn rejects_bad_parameters() {
        let p = points(vec![(0.0, 0.0, 0.0)]);
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            Near3dTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({ "input": p, "search_radius": 0 })));
        assert!(bad(json!({ "input": p, "search_radius": -5 })));
        assert!(bad(json!({ "input": p, "angle": "maybe" })));
    }

    /// The distance primitive itself, including the degenerate case.
    #[test]
    fn point_segment_distance() {
        let (d2, q) = point_seg_dist2([5.0, 0.0, 3.0], [0.0, 0.0, 0.0], [10.0, 0.0, 0.0]);
        assert!((d2.sqrt() - 3.0).abs() < 1e-9);
        assert!((q[0] - 5.0).abs() < 1e-9);
        // Beyond the end of the segment clamps to the endpoint.
        let (d2, q) = point_seg_dist2([20.0, 0.0, 0.0], [0.0, 0.0, 0.0], [10.0, 0.0, 0.0]);
        assert!((d2.sqrt() - 10.0).abs() < 1e-9);
        assert!((q[0] - 10.0).abs() < 1e-9);
        // Degenerate segment -> point-to-point, finite.
        let (d2, _) = point_seg_dist2([3.0, 4.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert!((d2.sqrt() - 5.0).abs() < 1e-9);
    }
}
