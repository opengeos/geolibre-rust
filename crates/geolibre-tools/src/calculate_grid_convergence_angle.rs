//! GeoLibre tool: per-feature grid convergence angle (grid north vs true north).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Calculate Grid Convergence Angle*
//! (Cartography). For each input feature it computes the centroid, takes its
//! geographic longitude/latitude and a central meridian, and stamps an attribute
//! field with the angle between grid north and true north — the value needed to
//! rotate north arrows, labels, and directional symbols correctly on a projected
//! map.
//!
//! The convergence is the standard transverse-Mercator/UTM approximation:
//!
//! ```text
//! γ = atan( tan(λ − λ0) · sin(φ) )
//! ```
//!
//! with λ the feature longitude, λ0 the central meridian, and φ the latitude
//! (all in radians; the result is returned in degrees). Two output conventions
//! are offered: *geographic* (degrees clockwise from north, positive east — the
//! raw γ) and *arithmetic* (degrees counter-clockwise from east, `90 − γ`).
//!
//! Note: the bundled `convergence_index` is an unrelated terrain flow-convergence
//! raster metric, not grid convergence.
//!
//! Scope for v1: geometry coordinates are treated as decimal-degree lon/lat
//! (EPSG:4326), matching how `calculate_utm_zone` reads geometry. Reprojecting a
//! projected input to geographic first would need PROJ, which the WASM-first
//! stack excludes; if the layer carries a non-geographic CRS the coordinate
//! values are still used as lon/lat and this limitation is surfaced in progress.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CalculateGridConvergenceAngleTool;

impl Tool for CalculateGridConvergenceAngleTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "calculate_grid_convergence_angle",
            display_name: "Calculate Grid Convergence Angle",
            summary: "Compute, per feature, the grid convergence angle between grid north and true north from the centroid lon/lat and a central meridian (default the feature's own UTM central meridian), writing it to a field for correct symbol/label rotation on projected maps — like ArcGIS Calculate Grid Convergence Angle. Geographic (clockwise from north) or arithmetic (counter-clockwise from east) convention. Pure math, no PROJ; unrelated to the bundled convergence_index terrain metric. Geometry is read as lon/lat degrees.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (any geometry type) in geographic coordinates (lon/lat degrees, EPSG:4326).",
                    required: true,
                },
                ToolParamSpec {
                    name: "angle_field",
                    description: "Name of the float field to hold the convergence angle in degrees. Default 'GRID_CONV'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "angle_type",
                    description: "Angle convention: 'geographic' (degrees clockwise from north, positive east) or 'arithmetic' (degrees counter-clockwise from east, 90 − geographic). Default 'geographic'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "central_meridian",
                    description: "Central meridian longitude in degrees. If omitted, each feature uses its own UTM-zone central meridian derived from its longitude.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;

        // v1 treats coordinates as lon/lat degrees. Warn (rather than reject) on a
        // non-geographic CRS so the caller knows the assumption in play.
        if let Some(epsg) = layer.crs_epsg() {
            if !is_geographic_epsg(epsg) {
                ctx.progress.info(&format!(
                    "input CRS EPSG:{epsg} is not geographic; coordinates are treated as lon/lat degrees (v1 supports geographic input only)"
                ));
            }
        }

        // Copy the schema and append the output field (reusing an existing field of
        // the same name so a re-run overwrites in place rather than duplicating).
        let mut out = Layer::new("grid_convergence");
        out.geom_type = layer.geom_type;
        out.schema = layer.schema.clone();
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        let angle_idx = ensure_field(&mut out, &prm.angle_field, FieldType::Float);
        let width = out.schema.fields().len();

        let mut annotated = 0usize;
        let mut skipped = 0usize;
        for feature in layer.iter() {
            let mut attrs = feature.attributes.clone();
            attrs.resize(width, FieldValue::Null);

            match feature.geometry.as_ref().and_then(centroid) {
                Some((lon, lat)) => {
                    let cm = prm
                        .central_meridian
                        .unwrap_or_else(|| utm_central_meridian(lon));
                    let angle = convergence_angle(lon, lat, cm, prm.angle_type);
                    attrs[angle_idx] = FieldValue::Float(angle);
                    annotated += 1;
                }
                None => {
                    attrs[angle_idx] = FieldValue::Null;
                    skipped += 1;
                }
            }

            out.push(Feature {
                fid: 0,
                geometry: feature.geometry.clone(),
                attributes: attrs,
            });
        }

        ctx.progress.info(&format!(
            "{annotated} feature(s) annotated, {skipped} skipped"
        ));

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(annotated + skipped));
        outputs.insert("annotated_count".to_string(), json!(annotated));
        outputs.insert("angle_field".to_string(), json!(prm.angle_field));
        outputs.insert("angle_type".to_string(), json!(prm.angle_type.as_str()));
        Ok(ToolRunResult { outputs })
    }
}

// ── Convergence math ──────────────────────────────────────────────────────────

/// Grid convergence angle in degrees for a feature at (`lon`, `lat`) relative to
/// central meridian `cm` (all degrees), in the requested convention.
fn convergence_angle(lon: f64, lat: f64, cm: f64, angle_type: AngleType) -> f64 {
    let dlon = (lon - cm).to_radians();
    let phi = lat.to_radians();
    // Geographic γ: degrees clockwise from north, positive east of the meridian.
    let geographic = (dlon.tan() * phi.sin()).atan().to_degrees();
    match angle_type {
        AngleType::Geographic => geographic,
        // Arithmetic: degrees counter-clockwise from east.
        AngleType::Arithmetic => 90.0 - geographic,
    }
}

/// UTM-zone central meridian (degrees) for a longitude: the mid-meridian of the
/// 6°-wide zone containing `lon`.
fn utm_central_meridian(lon: f64) -> f64 {
    let lon = ((lon + 180.0).rem_euclid(360.0)) - 180.0;
    ((lon + 180.0) / 6.0).floor() * 6.0 + 3.0 - 180.0
}

/// True for geographic (angular, lon/lat) CRSes we can accept directly.
fn is_geographic_epsg(epsg: u32) -> bool {
    matches!(epsg, 4326 | 4269 | 4267 | 4258 | 4283 | 4322 | 4324 | 4030)
}

// ── Centroid ──────────────────────────────────────────────────────────────────

/// Centroid of a geometry as (lon, lat). Polygons use the area-weighted (shoelace)
/// centroid of their exterior rings; lines and points fall back to the mean of
/// their vertices. Returns None for empty geometries.
fn centroid(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::Polygon { exterior, .. } => {
            ring_centroid(exterior.coords()).or_else(|| vertex_mean(geom))
        }
        Geometry::MultiPolygon(parts) => {
            let mut sx = 0.0;
            let mut sy = 0.0;
            let mut sa = 0.0;
            for (ext, _holes) in parts {
                if let Some((cx, cy, a)) = ring_signed_centroid(ext.coords()) {
                    sx += cx * a;
                    sy += cy * a;
                    sa += a;
                }
            }
            if sa.abs() > f64::EPSILON {
                Some((sx / sa, sy / sa))
            } else {
                vertex_mean(geom)
            }
        }
        _ => vertex_mean(geom),
    }
}

/// Area-weighted centroid of a single closed ring, or None if degenerate.
fn ring_centroid(coords: &[Coord]) -> Option<(f64, f64)> {
    ring_signed_centroid(coords).map(|(cx, cy, _a)| (cx, cy))
}

/// Signed shoelace centroid of a ring: returns (cx, cy, signed_area). Edges wrap
/// from the last vertex to the first, so the ring may be stored either closed
/// (first == last) or open — the closing edge is included either way.
fn ring_signed_centroid(coords: &[Coord]) -> Option<(f64, f64, f64)> {
    let n = coords.len();
    if n < 3 {
        return None;
    }
    let mut area2 = 0.0;
    let mut cx = 0.0;
    let mut cy = 0.0;
    for i in 0..n {
        let p = &coords[i];
        let q = &coords[(i + 1) % n];
        let cross = p.x * q.y - q.x * p.y;
        area2 += cross;
        cx += (p.x + q.x) * cross;
        cy += (p.y + q.y) * cross;
    }
    if area2.abs() < 1e-15 {
        return None;
    }
    let a = area2 / 2.0;
    Some((cx / (6.0 * a), cy / (6.0 * a), a))
}

/// Mean of every vertex in a geometry (used for points, lines, and degenerate
/// polygons where the shoelace area vanishes).
fn vertex_mean(geom: &Geometry) -> Option<(f64, f64)> {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0u64;
    accumulate(geom, &mut sx, &mut sy, &mut n);
    (n > 0).then(|| (sx / n as f64, sy / n as f64))
}

fn accumulate(geom: &Geometry, sx: &mut f64, sy: &mut f64, n: &mut u64) {
    let mut add = |c: &Coord| {
        *sx += c.x;
        *sy += c.y;
        *n += 1;
    };
    match geom {
        Geometry::Point(c) => add(c),
        Geometry::LineString(cs) | Geometry::MultiPoint(cs) => cs.iter().for_each(add),
        Geometry::MultiLineString(lines) => lines.iter().flatten().for_each(add),
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            exterior.coords().iter().for_each(&mut add);
            interiors
                .iter()
                .for_each(|r| r.coords().iter().for_each(&mut add));
        }
        Geometry::MultiPolygon(polys) => {
            for (ext, holes) in polys {
                ext.coords().iter().for_each(&mut add);
                holes
                    .iter()
                    .for_each(|r| r.coords().iter().for_each(&mut add));
            }
        }
        Geometry::GeometryCollection(geoms) => {
            for g in geoms {
                accumulate(g, sx, sy, n);
            }
        }
    }
}

// ── Fields ────────────────────────────────────────────────────────────────────

/// Returns the index of a field named `name`, adding it with `ty` if absent.
fn ensure_field(layer: &mut Layer, name: &str, ty: FieldType) -> usize {
    if let Some(idx) = layer.schema.field_index(name) {
        return idx;
    }
    layer.add_field(FieldDef::new(name, ty));
    layer.schema.field_index(name).expect("field just added")
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum AngleType {
    Geographic,
    Arithmetic,
}

impl AngleType {
    fn as_str(self) -> &'static str {
        match self {
            AngleType::Geographic => "geographic",
            AngleType::Arithmetic => "arithmetic",
        }
    }
}

struct Params {
    angle_field: String,
    angle_type: AngleType,
    central_meridian: Option<f64>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let angle_field = parse_optional_str(args, "angle_field")?
        .map(str::to_string)
        .unwrap_or_else(|| "GRID_CONV".to_string());

    let angle_type = match parse_optional_str(args, "angle_type")? {
        None => AngleType::Geographic,
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "geographic" => AngleType::Geographic,
            "arithmetic" => AngleType::Arithmetic,
            other => {
                return Err(ToolError::Validation(format!(
                    "angle_type must be 'geographic' or 'arithmetic', got '{other}'"
                )))
            }
        },
    };

    let central_meridian = match args.get("central_meridian") {
        None | Some(Value::Null) => None,
        Some(v) => Some(v.as_f64().ok_or_else(|| {
            ToolError::Validation("central_meridian must be a number when provided".into())
        })?),
    };

    Ok(Params {
        angle_field,
        angle_type,
        central_meridian,
    })
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
    use wbvector::{memory_store, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Point layer (lon/lat, EPSG:4326) with a "name" field.
    fn point_layer(rows: &[(&str, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("name", FieldType::Text));
        for (name, lon, lat) in rows {
            l.add_feature(
                Some(Geometry::point(*lon, *lat)),
                &[("name", (*name).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CalculateGridConvergenceAngleTool
            .run(&args, &ctx())
            .unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn field_f64(layer: &Layer, feat: usize, name: &str) -> f64 {
        let idx = layer.schema.field_index(name).unwrap();
        layer.features[feat].attributes[idx].as_f64().unwrap()
    }

    /// At the equator convergence is zero for any longitude offset.
    #[test]
    fn zero_at_equator() {
        let input = point_layer(&[("a", 5.0, 0.0), ("b", -30.0, 0.0)]);
        let (out, layer) = run(json!({ "input": input, "central_meridian": 0.0 }));
        assert_eq!(out.outputs["annotated_count"], json!(2));
        assert!(field_f64(&layer, 0, "GRID_CONV").abs() < 1e-9);
        assert!(field_f64(&layer, 1, "GRID_CONV").abs() < 1e-9);
    }

    /// Explicit deterministic value: lon=1, lat=60, cm=0 -> atan(tan(1°)·sin(60°)).
    #[test]
    fn known_value_and_sign_symmetry() {
        let input = point_layer(&[("east", 1.0, 60.0), ("west", -1.0, 60.0)]);
        let (_out, layer) = run(json!({ "input": input, "central_meridian": 0.0 }));
        let expected = ((1.0_f64.to_radians().tan()) * (60.0_f64.to_radians().sin()))
            .atan()
            .to_degrees();
        assert!((expected - 0.866_025).abs() < 1e-4, "sanity: {expected}");
        let east = field_f64(&layer, 0, "GRID_CONV");
        let west = field_f64(&layer, 1, "GRID_CONV");
        // East of the meridian is positive; west is the negative mirror.
        assert!((east - expected).abs() < 1e-6);
        assert!((west + expected).abs() < 1e-6);
        assert!(east > 0.0 && west < 0.0);
    }

    /// Magnitude grows with latitude for the same longitude offset.
    #[test]
    fn magnitude_increases_with_latitude() {
        let input = point_layer(&[("lo", 1.0, 20.0), ("hi", 1.0, 70.0)]);
        let (_out, layer) = run(json!({ "input": input, "central_meridian": 0.0 }));
        let lo = field_f64(&layer, 0, "GRID_CONV").abs();
        let hi = field_f64(&layer, 1, "GRID_CONV").abs();
        assert!(hi > lo);
    }

    /// Arithmetic convention is 90 − geographic.
    #[test]
    fn arithmetic_is_complement() {
        let input = point_layer(&[("p", 1.0, 60.0)]);
        let (_g, geo) = run(json!({ "input": input.clone(), "central_meridian": 0.0 }));
        let (_a, ari) = run(json!({
            "input": input, "central_meridian": 0.0, "angle_type": "arithmetic",
        }));
        let g = field_f64(&geo, 0, "GRID_CONV");
        let a = field_f64(&ari, 0, "GRID_CONV");
        assert!((a - (90.0 - g)).abs() < 1e-9);
    }

    /// Custom field name and the derived per-feature central meridian.
    #[test]
    fn custom_field_and_derived_meridian() {
        // lon=4 is in UTM zone 31 (cm=3); offset = 1°, lat=45.
        let input = point_layer(&[("p", 4.0, 45.0)]);
        let (out, layer) = run(json!({ "input": input, "angle_field": "conv" }));
        let expected = ((1.0_f64.to_radians().tan()) * (45.0_f64.to_radians().sin()))
            .atan()
            .to_degrees();
        assert!((field_f64(&layer, 0, "conv") - expected).abs() < 1e-6);
        assert_eq!(out.outputs["angle_field"], json!("conv"));
    }

    /// A feature with no geometry passes through with a null angle.
    #[test]
    fn passes_through_null_geometry() {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.push(Feature {
            fid: 0,
            geometry: None,
            attributes: vec![FieldValue::Text("empty".into())],
        });
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);
        let (out, layer) = run(json!({ "input": input }));
        assert_eq!(out.outputs["annotated_count"], json!(0));
        let idx = layer.schema.field_index("GRID_CONV").unwrap();
        assert!(matches!(
            layer.features[0].attributes[idx],
            FieldValue::Null
        ));
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CalculateGridConvergenceAngleTool.validate(&args)
        };
        // Missing input.
        assert!(bad(json!({})).is_err());
        // Bad angle_type.
        assert!(bad(json!({ "input": "a.geojson", "angle_type": "polar" })).is_err());
        // Valid.
        assert!(bad(json!({ "input": "a.geojson" })).is_ok());
        assert!(bad(json!({ "input": "a.geojson", "angle_type": "arithmetic" })).is_ok());
    }
}
