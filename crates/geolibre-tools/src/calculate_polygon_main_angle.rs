//! GeoLibre tool: dominant orientation angle of each polygon.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Calculate Polygon Main Angle*
//! (Cartography). It writes each polygon's dominant axis direction to an attribute
//! field, the standard input for rotating marker symbols and labels to follow
//! feature orientation (buildings, parcels, agricultural fields).
//!
//! No per-polygon orientation tool exists in either registry: the bundled
//! `patch_orientation` operates on raster patches, not vector polygons. This pairs
//! with the repo's footprint pipeline (`regularize_building_footprints`,
//! `regularize_adjacent_building_footprint` #410).
//!
//! The main angle is the direction of the longer side of the polygon's
//! minimum-area bounding rectangle, found by rotating calipers over the convex
//! hull. It is reported in one of three conventions:
//! * `arithmetic` (default) — counterclockwise from the positive x-axis;
//! * `geographic` — clockwise from north;
//! * `graphic` — clockwise from the positive x-axis (screen space).
//!
//! A line direction is defined modulo 180°, so every value is in `[0, 180)`.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CalculatePolygonMainAngleTool;

impl Tool for CalculatePolygonMainAngleTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "calculate_polygon_main_angle",
            display_name: "Calculate Polygon Main Angle",
            summary: "Write each polygon's dominant orientation (the long-side direction of its minimum-area bounding rectangle, via rotating calipers) to an angle field — like ArcGIS Calculate Polygon Main Angle. The per-polygon vector orientation the bundled raster-only patch_orientation lacks; feeds symbol/label rotation.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon layer with the angle field added. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "angle_field",
                    description: "Name of the angle field to write (default 'MainAngle').",
                    required: false,
                },
                ToolParamSpec {
                    name: "convention",
                    description: "'arithmetic' (CCW from east; default), 'geographic' (CW from north), or 'graphic' (CW from east).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_convention(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let field = parse_optional_str(args, "angle_field")?
            .unwrap_or("MainAngle")
            .to_string();
        let convention = parse_convention(args)?;

        let mut layer = load_input_layer(input)?;
        let n = layer.features.len();

        let mut angles = Vec::with_capacity(n);
        let mut degenerate = 0usize;
        for feat in &layer.features {
            let angle = feat
                .geometry
                .as_ref()
                .and_then(exterior_points)
                .and_then(|pts| main_angle_arithmetic(&pts));
            match angle {
                Some(a) => angles.push(convention.convert(a)),
                None => {
                    degenerate += 1;
                    angles.push(0.0);
                }
            }
        }

        ctx.progress
            .info(&format!("computed main angle for {n} polygon(s)"));

        layer.add_field(FieldDef::new(field.clone(), FieldType::Float));
        for (i, feat) in layer.features.iter_mut().enumerate() {
            feat.attributes.push(FieldValue::Float(angles[i]));
        }

        let out_path = write_or_store_layer(layer, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("angle_field".to_string(), json!(field));
        outputs.insert("degenerate".to_string(), json!(degenerate));
        Ok(ToolRunResult { outputs })
    }
}

/// All exterior-ring vertices of a (multi)polygon; other geometry -> None.
fn exterior_points(geom: &Geometry) -> Option<Vec<(f64, f64)>> {
    let ring_pts =
        |ring: &Ring| -> Vec<(f64, f64)> { ring.coords().iter().map(|c| (c.x, c.y)).collect() };
    let mut pts: Vec<(f64, f64)> = Vec::new();
    match geom {
        Geometry::Polygon { exterior, .. } => pts.extend(ring_pts(exterior)),
        Geometry::MultiPolygon(parts) => {
            for (ext, _holes) in parts {
                pts.extend(ring_pts(ext));
            }
        }
        _ => return None,
    }
    if pts.len() >= 3 {
        Some(pts)
    } else {
        None
    }
}

/// Main-axis angle (arithmetic, degrees in [0,180)) of the minimum-area bounding
/// rectangle of `pts`. Returns None for a degenerate (collinear/point) set.
fn main_angle_arithmetic(pts: &[(f64, f64)]) -> Option<f64> {
    let hull = convex_hull(pts);
    if hull.len() < 3 {
        return None;
    }
    let m = hull.len();
    let mut best_area = f64::INFINITY;
    let mut best_angle = 0.0;
    for i in 0..m {
        let (ax, ay) = hull[i];
        let (bx, by) = hull[(i + 1) % m];
        let (ex, ey) = (bx - ax, by - ay);
        let len = (ex * ex + ey * ey).sqrt();
        if len <= 0.0 {
            continue;
        }
        // Unit edge direction (ux,uy) and its perpendicular (-uy,ux).
        let (ux, uy) = (ex / len, ey / len);
        let (mut min_u, mut max_u, mut min_v, mut max_v) = (
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        );
        for &(px, py) in &hull {
            let u = px * ux + py * uy;
            let v = -px * uy + py * ux;
            min_u = min_u.min(u);
            max_u = max_u.max(u);
            min_v = min_v.min(v);
            max_v = max_v.max(v);
        }
        let w = max_u - min_u;
        let h = max_v - min_v;
        let area = w * h;
        if area < best_area {
            best_area = area;
            // Longer side direction: the edge axis if width >= height, else perp.
            best_angle = if w >= h {
                uy.atan2(ux).to_degrees()
            } else {
                ux.atan2(-uy).to_degrees()
            };
        }
    }
    // Normalize direction modulo 180 into [0,180).
    let mut a = best_angle % 180.0;
    if a < 0.0 {
        a += 180.0;
    }
    Some(a)
}

/// Andrew's monotone-chain convex hull (counterclockwise, no repeated endpoint).
fn convex_hull(pts: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut p: Vec<(f64, f64)> = pts.to_vec();
    p.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    p.dedup();
    let n = p.len();
    if n < 3 {
        return p;
    }
    let cross = |o: (f64, f64), a: (f64, f64), b: (f64, f64)| {
        (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
    };
    let mut hull: Vec<(f64, f64)> = Vec::with_capacity(2 * n);
    // Lower hull.
    for &pt in &p {
        while hull.len() >= 2 && cross(hull[hull.len() - 2], hull[hull.len() - 1], pt) <= 0.0 {
            hull.pop();
        }
        hull.push(pt);
    }
    // Upper hull.
    let lower = hull.len() + 1;
    for &pt in p.iter().rev() {
        while hull.len() >= lower && cross(hull[hull.len() - 2], hull[hull.len() - 1], pt) <= 0.0 {
            hull.pop();
        }
        hull.push(pt);
    }
    hull.pop();
    hull
}

// ── Conventions ───────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Convention {
    Arithmetic,
    Geographic,
    Graphic,
}

impl Convention {
    /// Converts an arithmetic angle (CCW from east, [0,180)) into this
    /// convention, keeping the result in [0,180).
    fn convert(self, arithmetic: f64) -> f64 {
        let a = match self {
            Convention::Arithmetic => arithmetic,
            Convention::Geographic => 90.0 - arithmetic,
            Convention::Graphic => -arithmetic,
        };
        let mut v = a % 180.0;
        if v < 0.0 {
            v += 180.0;
        }
        v
    }
}

fn parse_convention(args: &ToolArgs) -> Result<Convention, ToolError> {
    Ok(
        match args
            .get("convention")
            .and_then(Value::as_str)
            .map(str::trim)
        {
            None | Some("") | Some("arithmetic") => Convention::Arithmetic,
            Some("geographic") => Convention::Geographic,
            Some("graphic") => Convention::Graphic,
            Some(o) => {
                return Err(ToolError::Validation(format!(
                    "'convention' must be arithmetic|geographic|graphic, got '{o}'"
                )))
            }
        },
    )
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
    use wbvector::{memory_store, Coord, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// A layer with one rectangle rotated by `deg` degrees about the origin.
    fn rect_layer(w: f64, h: f64, deg: f64) -> String {
        let rad = deg.to_radians();
        let (cs, sn) = (rad.cos(), rad.sin());
        let rot = |x: f64, y: f64| (x * cs - y * sn, x * sn + y * cs);
        let corners = [(0.0, 0.0), (w, 0.0), (w, h), (0.0, h)];
        let mut ring: Vec<Coord> = corners
            .iter()
            .map(|&(x, y)| {
                let (rx, ry) = rot(x, y);
                Coord::xy(rx, ry)
            })
            .collect();
        ring.push(ring[0].clone());
        let mut l = Layer::new("r")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("id", FieldType::Integer));
        l.add_feature(
            Some(Geometry::Polygon {
                exterior: Ring::new(ring),
                interiors: vec![],
            }),
            &[("id", 1i64.into())],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> Layer {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CalculatePolygonMainAngleTool.run(&args, &ctx()).unwrap();
        load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    fn angle_of(l: &Layer, field: &str) -> f64 {
        let idx = l.schema.field_index(field).unwrap();
        l.features[0].attributes[idx].as_f64().unwrap()
    }

    /// A wide horizontal rectangle has arithmetic main angle ~0.
    #[test]
    fn horizontal_rect_is_zero() {
        let l = run(json!({ "input": rect_layer(10.0, 2.0, 0.0) }));
        let a = angle_of(&l, "MainAngle");
        assert!(!(1.0..=179.0).contains(&a), "expected ~0, got {a}");
    }

    /// A rectangle rotated 30° has arithmetic main angle ~30.
    #[test]
    fn rotated_rect_recovers_angle() {
        let l = run(json!({ "input": rect_layer(10.0, 3.0, 30.0) }));
        let a = angle_of(&l, "MainAngle");
        assert!((a - 30.0).abs() < 2.0, "expected ~30, got {a}");
    }

    /// Geographic convention of a horizontal (east-pointing) feature is ~90.
    #[test]
    fn geographic_of_horizontal_is_ninety() {
        let l = run(json!({ "input": rect_layer(10.0, 2.0, 0.0), "convention": "geographic" }));
        let a = angle_of(&l, "MainAngle");
        assert!((a - 90.0).abs() < 1.0, "expected ~90, got {a}");
    }

    /// Custom field name is honoured.
    #[test]
    fn custom_field() {
        let l = run(json!({ "input": rect_layer(5.0, 1.0, 45.0), "angle_field": "rot" }));
        assert!(l.schema.field_index("rot").is_some());
        let a = angle_of(&l, "rot");
        assert!((a - 45.0).abs() < 2.0, "expected ~45, got {a}");
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CalculatePolygonMainAngleTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "convention": "polar" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "convention": "geographic" })).is_ok());
    }
}
