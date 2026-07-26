//! GeoLibre tool: volumetric buffers around 3D points and lines.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Buffer 3D* (3D Analyst). The bundled
//! `buffer_vector` produces a 2D polygon regardless of input Z, which is the
//! wrong answer for any genuinely volumetric question — clearance around a
//! power line, a sensor detection envelope, a proximity zone around a flight
//! path. In plan view the buffer of a line at 100 m and one at 500 m are
//! identical.
//!
//! # Output representation
//!
//! `wbvector` is 2.5D: it has no mesh or solid type, but `Coord` carries an
//! optional Z and polygon rings may hold 3D vertices. The buffer is therefore
//! emitted as a **triangulated surface**: a `MultiPolygon` whose parts are
//! individual 3-vertex triangles with full XYZ coordinates. That keeps the
//! result inside the existing type system, round-trips through GeoJSON, and is
//! directly consumable by a renderer — at the cost of being a boundary
//! representation rather than a true solid.
//!
//! Two deliberate limitations follow, and they are stated in the output rather
//! than hidden: overlapping buffers are **not** unioned (there is no 3D boolean
//! op available), and `triangle_count` is reported so a caller can see the
//! tessellation cost before rendering.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

type P3 = [f64; 3];

/// End-cap style for line buffers.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Shape {
    /// Hemispherical caps (a capsule) — the true constant-distance envelope.
    Round,
    /// Flat caps (a cylinder) — cheaper, and what you want for a corridor with
    /// a defined start and end.
    Flat,
}

pub struct Buffer3dTool;

impl Tool for Buffer3dTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "buffer_3d",
            display_name: "Buffer 3D",
            summary: "Buffer 3D points and lines into volumetric geometry (spheres and capsules) emitted as a triangulated surface, like ArcGIS Buffer 3D.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input 3D point or polyline layer (missing Z is treated as 0).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distance",
                    description: "Buffer radius in map units.",
                    required: true,
                },
                ToolParamSpec {
                    name: "distance_field",
                    description: "Optional per-feature radius field, overriding 'distance'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "shape",
                    description: "round (default; hemispherical caps) or flat (cylinder with flat caps).",
                    required: false,
                },
                ToolParamSpec {
                    name: "quality",
                    description: "Tessellation resolution: segments around the circumference (default 12, min 4, max 64).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        let d = parse_optional_f64(args, "distance")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'distance'".to_string())
        })?;
        if !d.is_finite() || d <= 0.0 {
            return Err(ToolError::Validation(
                "'distance' must be greater than 0".to_string(),
            ));
        }
        if let Some(q) = parse_optional_f64(args, "quality")? {
            if !q.is_finite() || !(4.0..=64.0).contains(&q) {
                return Err(ToolError::Validation(
                    "'quality' must be between 4 and 64".to_string(),
                ));
            }
        }
        parse_shape(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let default_dist = parse_optional_f64(args, "distance")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'distance'".to_string())
        })?;
        let dist_field = parse_optional_str(args, "distance_field")?.map(String::from);
        let shape = parse_shape(args)?;
        let quality = parse_optional_f64(args, "quality")?.unwrap_or(12.0) as usize;

        let layer = load_input_layer(input)?;
        if let Some(f) = &dist_field {
            if layer.schema.field_index(f).is_none() {
                return Err(ToolError::Validation(format!(
                    "distance_field '{f}' not found on the input layer"
                )));
            }
        }

        let mut out = Layer::new(layer.name.clone());
        out.crs = layer.crs.clone();
        out.geom_type = Some(GeometryType::MultiPolygon);
        for fd in layer.schema.fields().iter() {
            out.add_field(fd.clone());
        }
        out.add_field(FieldDef::new("buffer_radius", FieldType::Float));
        out.add_field(FieldDef::new("triangle_count", FieldType::Integer));

        ctx.progress
            .info(&format!("buffering {} feature(s) in 3D", layer.len()));

        let mut total_triangles = 0usize;
        let mut buffered = 0usize;

        for (fi, feat) in layer.features.iter().enumerate() {
            let Some(geom) = feat.geometry.as_ref() else {
                continue;
            };
            let radius = match &dist_field {
                Some(f) => feat
                    .get(&layer.schema, f)
                    .ok()
                    .and_then(|v| v.as_f64())
                    .filter(|d| d.is_finite() && *d > 0.0)
                    .unwrap_or(default_dist),
                None => default_dist,
            };

            let tris = buffer_geometry(geom, radius, shape, quality);
            if tris.is_empty() {
                continue;
            }
            total_triangles += tris.len();
            buffered += 1;

            let mut fields: Vec<(String, FieldValue)> = layer
                .schema
                .fields()
                .iter()
                .enumerate()
                .map(|(i, fd)| (fd.name.clone(), feat.attributes[i].clone()))
                .collect();
            fields.push(("buffer_radius".into(), FieldValue::Float(radius)));
            fields.push((
                "triangle_count".into(),
                FieldValue::Integer(tris.len() as i64),
            ));
            let refs: Vec<(&str, FieldValue)> = fields
                .iter()
                .map(|(k, v)| (k.as_str(), v.clone()))
                .collect();
            out.add_feature(Some(triangles_to_geometry(&tris)), &refs)
                .map_err(|e| ToolError::Execution(format!("failed writing buffer: {e}")))?;
            ctx.progress
                .progress((fi as f64 + 1.0) / layer.len().max(1) as f64);
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(buffered));
        outputs.insert("triangle_count".to_string(), json!(total_triangles));
        // Made explicit so a caller does not mistake this for a solid-modelling
        // result: overlapping buffers are emitted separately, not merged.
        outputs.insert("unioned".to_string(), json!(false));
        outputs.insert(
            "representation".to_string(),
            json!("triangulated_surface_3d"),
        );
        Ok(ToolRunResult { outputs })
    }
}

/// One triangle of the output surface.
type Tri = [P3; 3];

fn z_of(c: &Coord) -> f64 {
    c.z.unwrap_or(0.0)
}

fn buffer_geometry(geom: &Geometry, r: f64, shape: Shape, q: usize) -> Vec<Tri> {
    let mut tris = Vec::new();
    match geom {
        Geometry::Point(c) => sphere(&mut tris, [c.x, c.y, z_of(c)], r, q),
        Geometry::MultiPoint(cs) => {
            for c in cs {
                sphere(&mut tris, [c.x, c.y, z_of(c)], r, q);
            }
        }
        Geometry::LineString(cs) => line_buffer(&mut tris, cs, r, shape, q),
        Geometry::MultiLineString(parts) => {
            for cs in parts {
                line_buffer(&mut tris, cs, r, shape, q);
            }
        }
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                tris.extend(buffer_geometry(g, r, shape, q));
            }
        }
        // Polygons have no defined 3D buffer here; they pass through unbuffered.
        _ => {}
    }
    tris
}

/// UV sphere at `c` with radius `r`, `q` segments around and `q/2` rings up.
fn sphere(out: &mut Vec<Tri>, c: P3, r: f64, q: usize) {
    let rings = (q / 2).max(2);
    let mut grid: Vec<Vec<P3>> = Vec::with_capacity(rings + 1);
    for i in 0..=rings {
        let phi = std::f64::consts::PI * i as f64 / rings as f64; // 0..pi
        let mut row = Vec::with_capacity(q + 1);
        for j in 0..=q {
            let theta = std::f64::consts::TAU * j as f64 / q as f64;
            row.push([
                c[0] + r * phi.sin() * theta.cos(),
                c[1] + r * phi.sin() * theta.sin(),
                c[2] + r * phi.cos(),
            ]);
        }
        grid.push(row);
    }
    for i in 0..rings {
        for j in 0..q {
            let (a, b, cc, d) = (
                grid[i][j],
                grid[i][j + 1],
                grid[i + 1][j + 1],
                grid[i + 1][j],
            );
            // The polar rows collapse to a point, so emit one triangle there
            // rather than a degenerate quad.
            if i == 0 {
                out.push([a, cc, d]);
            } else if i == rings - 1 {
                out.push([a, b, cc]);
            } else {
                out.push([a, b, cc]);
                out.push([a, cc, d]);
            }
        }
    }
}

/// Capsule (round) or cylinder (flat) along each segment of a polyline.
fn line_buffer(out: &mut Vec<Tri>, coords: &[Coord], r: f64, shape: Shape, q: usize) {
    if coords.len() < 2 {
        if let Some(c) = coords.first() {
            sphere(out, [c.x, c.y, z_of(c)], r, q);
        }
        return;
    }
    for w in coords.windows(2) {
        let a = [w[0].x, w[0].y, z_of(&w[0])];
        let b = [w[1].x, w[1].y, z_of(&w[1])];
        cylinder(out, a, b, r, q, shape == Shape::Flat);
    }
    if shape == Shape::Round {
        // A sphere at every vertex covers the joints between consecutive
        // segments, which avoids needing a 3D mitre computation.
        for c in coords {
            sphere(out, [c.x, c.y, z_of(c)], r, q);
        }
    }
}

/// Open cylinder from `a` to `b`, optionally capped with flat discs.
fn cylinder(out: &mut Vec<Tri>, a: P3, b: P3, r: f64, q: usize, capped: bool) {
    let axis = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let len = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    if len <= f64::EPSILON {
        sphere(out, a, r, q);
        return;
    }
    let n = [axis[0] / len, axis[1] / len, axis[2] / len];
    let (u, v) = orthonormal_basis(n);

    let ring = |p: P3| -> Vec<P3> {
        (0..=q)
            .map(|j| {
                let t = std::f64::consts::TAU * j as f64 / q as f64;
                [
                    p[0] + r * (u[0] * t.cos() + v[0] * t.sin()),
                    p[1] + r * (u[1] * t.cos() + v[1] * t.sin()),
                    p[2] + r * (u[2] * t.cos() + v[2] * t.sin()),
                ]
            })
            .collect()
    };
    let ra = ring(a);
    let rb = ring(b);
    for j in 0..q {
        out.push([ra[j], ra[j + 1], rb[j + 1]]);
        out.push([ra[j], rb[j + 1], rb[j]]);
    }
    if capped {
        for j in 0..q {
            out.push([a, ra[j], ra[j + 1]]);
            out.push([b, rb[j + 1], rb[j]]);
        }
    }
}

/// Two unit vectors perpendicular to `n` and to each other. Picking the seed
/// axis by smallest component avoids the degenerate case where the seed is
/// parallel to `n` (which would give a zero-length cross product).
fn orthonormal_basis(n: P3) -> (P3, P3) {
    let seed = if n[0].abs() <= n[1].abs() && n[0].abs() <= n[2].abs() {
        [1.0, 0.0, 0.0]
    } else if n[1].abs() <= n[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let u = normalize(cross(seed, n));
    let v = normalize(cross(n, u));
    (u, v)
}

fn cross(a: P3, b: P3) -> P3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalize(a: P3) -> P3 {
    let l = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    if l <= f64::EPSILON {
        [1.0, 0.0, 0.0]
    } else {
        [a[0] / l, a[1] / l, a[2] / l]
    }
}

/// Each triangle becomes one part of a MultiPolygon, with 3D vertices.
fn triangles_to_geometry(tris: &[Tri]) -> Geometry {
    Geometry::MultiPolygon(
        tris.iter()
            .map(|t| {
                (
                    Ring::new(
                        t.iter()
                            .map(|p| Coord::xyz(p[0], p[1], p[2]))
                            .collect::<Vec<_>>(),
                    ),
                    Vec::new(),
                )
            })
            .collect(),
    )
}

// ── parameter parsing ────────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

fn parse_shape(args: &ToolArgs) -> Result<Shape, ToolError> {
    match args.get("shape").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("round") => Ok(Shape::Round),
        Some("flat") => Ok(Shape::Flat),
        Some(o) => Err(ToolError::Validation(format!(
            "'shape' must be round or flat, got '{o}'"
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

    fn point_layer(items: Vec<(f64, f64, f64)>) -> String {
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("r", FieldType::Float));
        for (x, y, z) in items {
            l.add_feature(
                Some(Geometry::Point(Coord::xyz(x, y, z))),
                &[("r", FieldValue::Float(7.0))],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn line_layer(cs: Vec<Coord>) -> String {
        let mut l = Layer::new("l")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_feature(Some(Geometry::LineString(cs)), &[]).unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = Buffer3dTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// Every emitted vertex of a point buffer sits on the sphere of the given
    /// radius — the defining property of the geometry.
    #[test]
    fn sphere_vertices_lie_on_the_radius() {
        let (_, layer) = run(json!({
            "input": point_layer(vec![(10.0, 20.0, 30.0)]),
            "distance": 5.0, "quality": 8
        }));
        let Geometry::MultiPolygon(parts) = layer.features[0].geometry.as_ref().unwrap() else {
            panic!("expected MultiPolygon");
        };
        assert!(!parts.is_empty());
        for (ring, _) in parts {
            for c in ring.coords() {
                let d =
                    ((c.x - 10.0).powi(2) + (c.y - 20.0).powi(2) + (c.z.unwrap() - 30.0).powi(2))
                        .sqrt();
                assert!(
                    (d - 5.0).abs() < 1e-6,
                    "vertex at distance {d} from the centre, expected 5"
                );
            }
        }
    }

    /// THE regression versus a 2D buffer: Z is carried, so the envelope is
    /// centred on the feature's elevation rather than flattened to the ground.
    #[test]
    fn buffer_is_centred_on_feature_elevation() {
        let (_, layer) = run(json!({
            "input": point_layer(vec![(0.0, 0.0, 500.0)]),
            "distance": 10.0, "quality": 8
        }));
        let Geometry::MultiPolygon(parts) = layer.features[0].geometry.as_ref().unwrap() else {
            panic!("expected MultiPolygon");
        };
        let zs: Vec<f64> = parts
            .iter()
            .flat_map(|(r, _)| r.coords().iter().map(|c| c.z.unwrap()))
            .collect();
        let (lo, hi) = (
            zs.iter().cloned().fold(f64::MAX, f64::min),
            zs.iter().cloned().fold(f64::MIN, f64::max),
        );
        assert!((lo - 490.0).abs() < 1e-6, "bottom at {lo}, expected 490");
        assert!((hi - 510.0).abs() < 1e-6, "top at {hi}, expected 510");
    }

    /// A cylinder's side vertices are exactly `distance` from the axis, which
    /// is what makes a line buffer a constant-clearance envelope.
    #[test]
    fn cylinder_sides_are_at_the_radius() {
        let mut tris = Vec::new();
        cylinder(&mut tris, [0.0, 0.0, 0.0], [100.0, 0.0, 0.0], 3.0, 8, false);
        assert!(!tris.is_empty());
        for t in &tris {
            for p in t {
                // Distance from the x-axis.
                let d = (p[1] * p[1] + p[2] * p[2]).sqrt();
                assert!((d - 3.0).abs() < 1e-9, "side vertex at {d}, expected 3");
            }
        }
    }

    /// A vertical line must buffer correctly — the degenerate case for a naive
    /// basis that seeds on the Z axis.
    #[test]
    fn vertical_line_does_not_degenerate() {
        let (out, _) = run(json!({
            "input": line_layer(vec![
                Coord::xyz(0.0, 0.0, 0.0),
                Coord::xyz(0.0, 0.0, 100.0),
            ]),
            "distance": 4.0, "quality": 8
        }));
        assert!(out.outputs["triangle_count"].as_f64().unwrap() > 0.0);

        // And the basis itself is well-formed for a Z-aligned axis.
        let (u, v) = orthonormal_basis([0.0, 0.0, 1.0]);
        let dot = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
        assert!(dot.abs() < 1e-12, "basis vectors must be perpendicular");
        for w in [u, v] {
            let l = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
            assert!((l - 1.0).abs() < 1e-12, "basis vectors must be unit length");
        }
    }

    /// round adds joint spheres; flat does not.
    #[test]
    fn shape_controls_cap_geometry() {
        let cs = vec![
            Coord::xyz(0.0, 0.0, 0.0),
            Coord::xyz(10.0, 0.0, 0.0),
            Coord::xyz(20.0, 10.0, 0.0),
        ];
        let (round, _) = run(json!({
            "input": line_layer(cs.clone()), "distance": 2.0,
            "quality": 8, "shape": "round"
        }));
        let (flat, _) = run(json!({
            "input": line_layer(cs), "distance": 2.0,
            "quality": 8, "shape": "flat"
        }));
        assert!(
            round.outputs["triangle_count"].as_f64().unwrap()
                > flat.outputs["triangle_count"].as_f64().unwrap(),
            "round caps add joint spheres"
        );
    }

    /// distance_field overrides the default radius per feature.
    #[test]
    fn distance_field_overrides_default() {
        let (_, layer) = run(json!({
            "input": point_layer(vec![(0.0, 0.0, 0.0)]),
            "distance": 1.0, "distance_field": "r", "quality": 8
        }));
        let ri = layer.schema.field_index("buffer_radius").unwrap();
        assert_eq!(layer.features[0].attributes[ri].as_f64(), Some(7.0));
    }

    /// quality controls tessellation density.
    #[test]
    fn quality_controls_triangle_count() {
        let input = point_layer(vec![(0.0, 0.0, 0.0)]);
        let (coarse, _) = run(json!({ "input": input, "distance": 5.0, "quality": 4 }));
        let (fine, _) = run(json!({ "input": input, "distance": 5.0, "quality": 32 }));
        assert!(
            fine.outputs["triangle_count"].as_f64().unwrap()
                > coarse.outputs["triangle_count"].as_f64().unwrap()
        );
    }

    /// The no-union limitation is reported, not hidden.
    #[test]
    fn reports_representation_and_union_limitation() {
        let (out, _) = run(json!({
            "input": point_layer(vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]),
            "distance": 10.0, "quality": 6
        }));
        // Two overlapping spheres, emitted as two separate features.
        assert_eq!(out.outputs["feature_count"], json!(2));
        assert_eq!(out.outputs["unioned"], json!(false));
        assert_eq!(
            out.outputs["representation"],
            json!("triangulated_surface_3d")
        );
    }

    #[test]
    fn rejects_bad_parameters() {
        let p = point_layer(vec![(0.0, 0.0, 0.0)]);
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            Buffer3dTool.validate(&args).is_err()
        };
        assert!(bad(json!({ "input": p })));
        assert!(bad(json!({ "input": p, "distance": 0 })));
        assert!(bad(json!({ "input": p, "distance": -1 })));
        assert!(bad(json!({ "input": p, "distance": 1, "quality": 2 })));
        assert!(bad(json!({ "input": p, "distance": 1, "quality": 100 })));
        assert!(bad(json!({ "input": p, "distance": 1, "shape": "square" })));
        assert!(bad(json!({ "distance": 1 })));
    }
}
