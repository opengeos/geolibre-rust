//! GeoLibre tool: Z-aware Douglas-Peucker simplification of 3D polylines.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Simplify 3D Line* (3D Analyst). The
//! repo carries 3D data through `interpolate_shape`, `add_surface_information`,
//! `adjust_3d_z` and `las_height_metrics`, but has no Z-aware geometry
//! *operator*: every vector operation available works in plan view and discards
//! Z. For simplification that is not cosmetic. The bundled `simplify_features`
//! measures perpendicular distance in 2D, so a line draped over a ridge can have
//! its whole vertical profile flattened while the planimetric tolerance is
//! respected — the crest vertex is planimetrically redundant and gets dropped.
//!
//! This tool splits on the perpendicular distance from a vertex to the 3D chord
//! joining the segment endpoints, computed as
//!
//! ```text
//!   d = |(p - a) x (b - a)| / |b - a|
//! ```
//!
//! Retained vertices keep their Z **exactly** — nothing is re-interpolated — and
//! endpoints are always kept. `z_factor` scales Z before the distance test, for
//! data whose Z unit differs from its XY unit (e.g. metres of elevation on a
//! degree-based CRS, or feet on metres).
//!
//! Vertices with no Z are treated as 2D for the distance test (the Z term drops
//! out), so a mixed layer degrades to plain Douglas-Peucker rather than failing.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Geometry, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct Simplify3dLineTool;

impl Tool for Simplify3dLineTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "simplify_3d_line",
            display_name: "Simplify 3D Line",
            summary: "Reduce the vertex count of 3D polylines using a Douglas-Peucker split that measures deviation in three dimensions, preserving the vertical profile, like ArcGIS Simplify 3D Line.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polyline layer (3D vertices; 2D vertices are handled as a planar special case).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Maximum 3D deviation in map units. Vertices closer than this to the chord are removed.",
                    required: true,
                },
                ToolParamSpec {
                    name: "z_factor",
                    description: "Multiplier applied to Z before the distance test, for data whose Z unit differs from XY (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        let tol = parse_optional_f64(args, "tolerance")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'tolerance'".to_string())
        })?;
        if !tol.is_finite() || tol <= 0.0 {
            return Err(ToolError::Validation(
                "'tolerance' must be greater than 0".to_string(),
            ));
        }
        if let Some(zf) = parse_optional_f64(args, "z_factor")? {
            if !zf.is_finite() || zf < 0.0 {
                return Err(ToolError::Validation(
                    "'z_factor' must be a non-negative number".to_string(),
                ));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let tolerance = parse_optional_f64(args, "tolerance")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'tolerance'".to_string())
        })?;
        let z_factor = parse_optional_f64(args, "z_factor")?.unwrap_or(1.0);

        let layer = load_input_layer(input)?;
        let mut out = Layer::new(layer.name.clone());
        out.crs = layer.crs.clone();
        out.geom_type = layer.geom_type;
        for fd in layer.schema.fields().iter() {
            out.add_field(fd.clone());
        }

        ctx.progress
            .info(&format!("simplifying {} feature(s)", layer.len()));

        let mut vertices_before = 0usize;
        let mut vertices_after = 0usize;
        let mut simplified = 0usize;

        for (fi, feat) in layer.features.iter().enumerate() {
            let geom = match feat.geometry.as_ref() {
                Some(g) => {
                    let (new_geom, before, after, touched) =
                        simplify_geometry(g, tolerance, z_factor);
                    vertices_before += before;
                    vertices_after += after;
                    if touched {
                        simplified += 1;
                    }
                    Some(new_geom)
                }
                None => None,
            };
            let fields: Vec<(&str, wbvector::FieldValue)> = layer
                .schema
                .fields()
                .iter()
                .enumerate()
                .map(|(i, fd)| (fd.name.as_str(), feat.attributes[i].clone()))
                .collect();
            out.add_feature(geom, &fields)
                .map_err(|e| ToolError::Execution(format!("failed writing feature: {e}")))?;
            ctx.progress
                .progress((fi as f64 + 1.0) / layer.len().max(1) as f64);
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(layer.len()));
        outputs.insert("features_simplified".to_string(), json!(simplified));
        outputs.insert("vertices_before".to_string(), json!(vertices_before));
        outputs.insert("vertices_after".to_string(), json!(vertices_after));
        outputs.insert(
            "vertices_removed".to_string(),
            json!(vertices_before.saturating_sub(vertices_after)),
        );
        Ok(ToolRunResult { outputs })
    }
}

/// Simplifies every linear part of a geometry. Returns
/// `(geometry, vertices_before, vertices_after, changed)`. Non-linear geometry
/// passes through untouched.
fn simplify_geometry(
    geom: &Geometry,
    tolerance: f64,
    z_factor: f64,
) -> (Geometry, usize, usize, bool) {
    match geom {
        Geometry::LineString(cs) => {
            let before = cs.len();
            let simplified = douglas_peucker_3d(cs, tolerance, z_factor);
            let after = simplified.len();
            (
                Geometry::LineString(simplified),
                before,
                after,
                after < before,
            )
        }
        Geometry::MultiLineString(parts) => {
            let mut before = 0;
            let mut after = 0;
            let mut out = Vec::with_capacity(parts.len());
            for cs in parts {
                before += cs.len();
                let s = douglas_peucker_3d(cs, tolerance, z_factor);
                after += s.len();
                out.push(s);
            }
            (
                Geometry::MultiLineString(out),
                before,
                after,
                after < before,
            )
        }
        other => (other.clone(), 0, 0, false),
    }
}

/// Recursive Douglas-Peucker using 3D point-to-chord distance.
fn douglas_peucker_3d(coords: &[Coord], tolerance: f64, z_factor: f64) -> Vec<Coord> {
    if coords.len() <= 2 {
        return coords.to_vec();
    }
    let mut keep = vec![false; coords.len()];
    keep[0] = true;
    keep[coords.len() - 1] = true;
    dp_recurse(coords, 0, coords.len() - 1, tolerance, z_factor, &mut keep);
    coords
        .iter()
        .zip(keep.iter())
        .filter_map(|(c, k)| if *k { Some(c.clone()) } else { None })
        .collect()
}

fn dp_recurse(
    coords: &[Coord],
    first: usize,
    last: usize,
    tolerance: f64,
    z_factor: f64,
    keep: &mut [bool],
) {
    let mut stack = vec![(first, last)];
    while let Some((first, last)) = stack.pop() {
        if last <= first + 1 {
            continue;
        }
        let mut max_d = -1.0_f64;
        let mut max_i = first;
        for i in (first + 1)..last {
            let d = point_chord_distance_3d(&coords[i], &coords[first], &coords[last], z_factor);
            if d > max_d {
                max_d = d;
                max_i = i;
            }
        }
        if max_d > tolerance {
            keep[max_i] = true;
            stack.push((max_i, last));
            stack.push((first, max_i));
        }
    }
}

/// Perpendicular distance from `p` to the 3D segment chord `a`-`b`, via the
/// cross-product formula. A degenerate (zero-length) chord falls back to plain
/// point-to-point distance so coincident endpoints cannot produce NaN.
fn point_chord_distance_3d(p: &Coord, a: &Coord, b: &Coord, z_factor: f64) -> f64 {
    let z = |c: &Coord| c.z.unwrap_or(0.0) * z_factor;
    let (px, py, pz) = (p.x, p.y, z(p));
    let (ax, ay, az) = (a.x, a.y, z(a));
    let (bx, by, bz) = (b.x, b.y, z(b));

    let (ux, uy, uz) = (bx - ax, by - ay, bz - az);
    let (vx, vy, vz) = (px - ax, py - ay, pz - az);
    let chord2 = ux * ux + uy * uy + uz * uz;
    if chord2 <= f64::EPSILON {
        return (vx * vx + vy * vy + vz * vz).sqrt();
    }
    // |v x u| / |u|
    let cx = vy * uz - vz * uy;
    let cy = vz * ux - vx * uz;
    let cz = vx * uy - vy * ux;
    ((cx * cx + cy * cy + cz * cz) / chord2).sqrt()
}

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

    fn line_layer(coords: Vec<Coord>) -> String {
        let mut l = Layer::new("lines")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_feature(Some(Geometry::LineString(coords)), &[])
            .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = Simplify3dLineTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn coords_of(layer: &Layer, i: usize) -> Vec<Coord> {
        match layer.features[i].geometry.as_ref().unwrap() {
            Geometry::LineString(cs) => cs.clone(),
            other => panic!("expected LineString, got {other:?}"),
        }
    }

    /// THE regression this tool exists for: a vertex that is planimetrically
    /// redundant but vertically significant must survive.
    #[test]
    fn keeps_vertically_significant_vertex() {
        // Straight in plan (y=0, x = 0,5,10) but the middle vertex is 50 units up.
        let cs = vec![
            Coord::xyz(0.0, 0.0, 0.0),
            Coord::xyz(5.0, 0.0, 50.0),
            Coord::xyz(10.0, 0.0, 0.0),
        ];
        let (_, layer) = run(json!({ "input": line_layer(cs), "tolerance": 1.0 }));
        let out = coords_of(&layer, 0);
        assert_eq!(out.len(), 3, "the ridge crest must not be dropped");
        assert_eq!(out[1].z, Some(50.0), "retained Z must be exact");
    }

    /// A vertex that is collinear in 3D is genuinely redundant and goes.
    #[test]
    fn removes_collinear_3d_vertex() {
        let cs = vec![
            Coord::xyz(0.0, 0.0, 0.0),
            Coord::xyz(5.0, 0.0, 5.0),
            Coord::xyz(10.0, 0.0, 10.0),
        ];
        let (out, layer) = run(json!({ "input": line_layer(cs), "tolerance": 1.0 }));
        assert_eq!(coords_of(&layer, 0).len(), 2);
        assert_eq!(out.outputs["vertices_removed"], json!(1));
    }

    /// z_factor scales the vertical term: a small Z bump becomes significant.
    #[test]
    fn z_factor_amplifies_vertical_deviation() {
        let cs = vec![
            Coord::xyz(0.0, 0.0, 0.0),
            Coord::xyz(5.0, 0.0, 0.5),
            Coord::xyz(10.0, 0.0, 0.0),
        ];
        // Unscaled the 0.5 bump is under tolerance -> dropped.
        let (_, plain) = run(json!({ "input": line_layer(cs.clone()), "tolerance": 1.0 }));
        assert_eq!(coords_of(&plain, 0).len(), 2);
        // Scaled 10x it clears tolerance -> kept.
        let (_, scaled) =
            run(json!({ "input": line_layer(cs), "tolerance": 1.0, "z_factor": 10.0 }));
        assert_eq!(coords_of(&scaled, 0).len(), 3);
    }

    /// Endpoints are never removed, even below tolerance.
    #[test]
    fn endpoints_always_retained() {
        let cs = vec![
            Coord::xyz(0.0, 0.0, 0.0),
            Coord::xyz(0.001, 0.0, 0.0),
            Coord::xyz(0.002, 0.0, 0.0),
        ];
        let (_, layer) = run(json!({ "input": line_layer(cs), "tolerance": 100.0 }));
        let out = coords_of(&layer, 0);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].x, 0.0);
        assert_eq!(out[1].x, 0.002);
    }

    /// 2D vertices (no Z) degrade to plain Douglas-Peucker instead of failing.
    #[test]
    fn handles_2d_vertices() {
        let cs = vec![
            Coord::xy(0.0, 0.0),
            Coord::xy(5.0, 0.1),
            Coord::xy(10.0, 0.0),
        ];
        let (_, layer) = run(json!({ "input": line_layer(cs), "tolerance": 1.0 }));
        assert_eq!(coords_of(&layer, 0).len(), 2);
    }

    /// Non-linear geometry passes through untouched.
    #[test]
    fn passes_through_non_line_geometry() {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_feature(Some(Geometry::Point(Coord::xyz(1.0, 2.0, 3.0))), &[])
            .unwrap();
        let id = memory_store::put_vector(l);
        let path = memory_store::make_vector_memory_path(&id);
        let (_, layer) = run(json!({ "input": path, "tolerance": 1.0 }));
        assert!(matches!(
            layer.features[0].geometry.as_ref().unwrap(),
            Geometry::Point(_)
        ));
    }

    #[test]
    fn rejects_bad_parameters() {
        let p = line_layer(vec![Coord::xyz(0.0, 0.0, 0.0), Coord::xyz(1.0, 1.0, 1.0)]);
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            Simplify3dLineTool.validate(&args).is_err()
        };
        assert!(bad(json!({ "input": p })), "tolerance is required");
        assert!(bad(json!({ "input": p, "tolerance": 0 })));
        assert!(bad(json!({ "input": p, "tolerance": -1 })));
        assert!(bad(json!({ "input": p, "tolerance": 1, "z_factor": -2 })));
        assert!(bad(json!({ "tolerance": 1 })));
    }

    /// The distance formula itself, including the degenerate-chord fallback.
    #[test]
    fn distance_formula() {
        let a = Coord::xyz(0.0, 0.0, 0.0);
        let b = Coord::xyz(10.0, 0.0, 0.0);
        // 3 units above the midpoint of a horizontal chord.
        let p = Coord::xyz(5.0, 0.0, 3.0);
        assert!((point_chord_distance_3d(&p, &a, &b, 1.0) - 3.0).abs() < 1e-9);
        // Degenerate chord -> point-to-point distance, not NaN.
        let d = point_chord_distance_3d(&p, &a, &a, 1.0);
        assert!(d.is_finite());
        assert!((d - (25.0_f64 + 9.0).sqrt()).abs() < 1e-9);
    }
}
