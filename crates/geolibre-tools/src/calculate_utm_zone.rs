//! GeoLibre tool: calculate the optimal UTM zone / EPSG code per feature.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Calculate UTM Zone* (Cartography). For
//! each input feature it computes the centroid, derives the UTM zone from the
//! centroid longitude, and stamps two attribute fields: the zone number (1–60)
//! and the WGS84 UTM EPSG code (326## for the northern hemisphere, 327## for the
//! southern). This lets a map series pick the least-distorted projection per
//! sheet without PROJ.
//!
//! The zone is `floor((lon + 180) / 6) + 1`, clamped to 1–60, with the two
//! standard exceptions that ArcGIS also applies:
//!
//! * **Norway** — zone 32 is widened west over zone 31 for latitudes 56°–64°N,
//!   longitudes 3°–12°E.
//! * **Svalbard** — for 72°–84°N the odd zones 31/33/35/37 are widened and the
//!   even zones 32/34/36 suppressed.
//!
//! Scope for v1: geometry coordinates are treated as decimal-degree lon/lat
//! (EPSG:4326), matching how `convert_coordinate_notation` reads point geometry.
//! Reprojecting a projected input to geographic first would need PROJ, which the
//! WASM-first stack excludes; features are validated to carry a geographic (or
//! unspecified) CRS and rejected otherwise.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CalculateUtmZoneTool;

impl Tool for CalculateUtmZoneTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "calculate_utm_zone",
            display_name: "Calculate UTM Zone",
            summary: "Compute the optimal UTM zone and WGS84 EPSG code for each feature from its centroid longitude/latitude, adding a zone (1–60) and an EPSG (326##/327##) field — for adaptive-projection map series, like ArcGIS Calculate UTM Zone. Pure math with the Norway/Svalbard exceptions; no PROJ. Geometry is read as lon/lat degrees.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (any geometry type) in geographic coordinates (lon/lat degrees, EPSG:4326).",
                    required: true,
                },
                ToolParamSpec {
                    name: "zone_field",
                    description: "Name of the integer field to hold the UTM zone number (1–60). Default 'utm_zone'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "epsg_field",
                    description: "Name of the integer field to hold the WGS84 UTM EPSG code (326## north, 327## south). Default 'utm_epsg'.",
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

        // Guard against a projected CRS: the centroid math assumes lon/lat degrees.
        // EPSG:4326 (and the legacy 4269 NAD83) are geographic; treat an
        // unspecified CRS as geographic too (memory layers, plain GeoJSON).
        if let Some(epsg) = layer.crs_epsg() {
            if !is_geographic_epsg(epsg) {
                return Err(ToolError::Validation(format!(
                    "input CRS EPSG:{epsg} is not geographic; calculate_utm_zone needs lon/lat degrees (EPSG:4326). Reproject first."
                )));
            }
        }

        // Copy the schema and append the two output fields (skipping names that
        // already exist so a re-run overwrites in place rather than duplicating).
        let mut out = Layer::new("utm_zone");
        out.geom_type = layer.geom_type;
        out.schema = layer.schema.clone();
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        let zone_idx = ensure_field(&mut out, &prm.zone_field, FieldType::Integer);
        let epsg_idx = ensure_field(&mut out, &prm.epsg_field, FieldType::Integer);
        let width = out.schema.fields().len();

        let mut stamped = 0usize;
        let mut skipped = 0usize;
        for feature in layer.iter() {
            let mut attrs = feature.attributes.clone();
            attrs.resize(width, FieldValue::Null);

            match feature.geometry.as_ref().and_then(centroid) {
                Some((lon, lat)) => {
                    let zone = utm_zone(lon, lat);
                    let epsg = utm_epsg(zone, lat >= 0.0);
                    attrs[zone_idx] = FieldValue::Integer(zone as i64);
                    attrs[epsg_idx] = FieldValue::Integer(epsg as i64);
                    stamped += 1;
                }
                None => {
                    attrs[zone_idx] = FieldValue::Null;
                    attrs[epsg_idx] = FieldValue::Null;
                    skipped += 1;
                }
            }

            out.push(Feature {
                fid: 0,
                geometry: feature.geometry.clone(),
                attributes: attrs,
            });
        }

        ctx.progress
            .info(&format!("{stamped} feature(s) stamped, {skipped} skipped"));

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("stamped".to_string(), json!(stamped));
        outputs.insert("skipped".to_string(), json!(skipped));
        outputs.insert("zone_field".to_string(), json!(prm.zone_field));
        outputs.insert("epsg_field".to_string(), json!(prm.epsg_field));
        Ok(ToolRunResult { outputs })
    }
}

// ── UTM zone math ─────────────────────────────────────────────────────────────

/// UTM zone (1–60) for a centroid at (`lon`, `lat`) in degrees, applying the
/// standard Norway and Svalbard exceptions.
fn utm_zone(lon: f64, lat: f64) -> i32 {
    // Wrap longitude into [-180, 180) so features at exactly 180° map to zone 60.
    let lon = ((lon + 180.0).rem_euclid(360.0)) - 180.0;
    let mut zone = ((lon + 180.0) / 6.0).floor() as i32 + 1;
    zone = zone.clamp(1, 60);

    // Norway: zone 32 is widened west across the 6°E meridian for 56°–64°N.
    if (56.0..64.0).contains(&lat) && (3.0..12.0).contains(&lon) {
        zone = 32;
    }

    // Svalbard: odd zones widened, even zones suppressed, for 72°–84°N.
    if (72.0..84.0).contains(&lat) {
        zone = if (0.0..9.0).contains(&lon) {
            31
        } else if (9.0..21.0).contains(&lon) {
            33
        } else if (21.0..33.0).contains(&lon) {
            35
        } else if (33.0..42.0).contains(&lon) {
            37
        } else {
            zone
        };
    }

    zone
}

/// WGS84 UTM EPSG code: 326## for the northern hemisphere, 327## for the southern.
fn utm_epsg(zone: i32, is_north: bool) -> i32 {
    if is_north {
        32600 + zone
    } else {
        32700 + zone
    }
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
/// around from the last vertex to the first, so the ring may be stored either
/// closed (first == last) or open — the closing edge is included either way (it
/// contributes zero when the ring is already closed).
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

struct Params {
    zone_field: String,
    epsg_field: String,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let zone_field = parse_optional_str(args, "zone_field")?
        .map(str::to_string)
        .unwrap_or_else(|| "utm_zone".to_string());
    let epsg_field = parse_optional_str(args, "epsg_field")?
        .map(str::to_string)
        .unwrap_or_else(|| "utm_epsg".to_string());
    if zone_field == epsg_field {
        return Err(ToolError::Validation(
            "zone_field and epsg_field must be different field names".into(),
        ));
    }
    Ok(Params {
        zone_field,
        epsg_field,
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
        let out = CalculateUtmZoneTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn field_i64(layer: &Layer, feat: usize, name: &str) -> i64 {
        let idx = layer.schema.field_index(name).unwrap();
        layer.features[feat].attributes[idx].as_i64().unwrap()
    }

    /// Zone formula and EPSG for a handful of known locations.
    #[test]
    fn known_zones_and_epsg() {
        // White House 18N, Eiffel Tower 31N, Sydney Opera House 56S.
        let input = point_layer(&[
            ("wh", -77.0365, 38.8977),
            ("eiffel", 2.2945, 48.8584),
            ("sydney", 151.2153, -33.8568),
        ]);
        let (out, layer) = run(json!({ "input": input }));
        assert_eq!(out.outputs["stamped"], json!(3));
        assert_eq!(field_i64(&layer, 0, "utm_zone"), 18);
        assert_eq!(field_i64(&layer, 0, "utm_epsg"), 32618);
        assert_eq!(field_i64(&layer, 1, "utm_zone"), 31);
        assert_eq!(field_i64(&layer, 1, "utm_epsg"), 32631);
        assert_eq!(field_i64(&layer, 2, "utm_zone"), 56);
        assert_eq!(field_i64(&layer, 2, "utm_epsg"), 32756); // southern -> 327##
    }

    /// Norway widens zone 32 west of 6°E for 56°–64°N.
    #[test]
    fn norway_exception() {
        // Bergen ~ (5.32, 60.39) would be zone 31 by the plain formula, but is 32.
        assert_eq!(utm_zone(5.32, 60.39), 32);
        // Just south of the band, the plain formula holds.
        assert_eq!(utm_zone(5.32, 55.0), 31);
    }

    /// Svalbard suppresses even zones; 20°E at 78°N is zone 33, not 34.
    #[test]
    fn svalbard_exception() {
        assert_eq!(utm_zone(20.0, 78.0), 33);
        assert_eq!(utm_zone(8.0, 78.0), 31);
        assert_eq!(utm_zone(25.0, 78.0), 35);
        assert_eq!(utm_zone(40.0, 78.0), 37);
    }

    /// A polygon is stamped from its area centroid, not a stray vertex.
    #[test]
    fn polygon_uses_area_centroid() {
        let mut l = Layer::new("poly")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("name", FieldType::Text));
        let ring = Geometry::polygon(
            vec![
                Coord::xy(2.0, 40.0),
                Coord::xy(4.0, 40.0),
                Coord::xy(4.0, 42.0),
                Coord::xy(2.0, 42.0),
                Coord::xy(2.0, 40.0),
            ],
            vec![],
        );
        l.add_feature(Some(ring), &[("name", "sq".into())]).unwrap();
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);
        let (out, layer) = run(json!({ "input": input }));
        assert_eq!(out.outputs["stamped"], json!(1));
        // Centroid lon = 3 -> zone floor((3+180)/6)+1 = 31.
        assert_eq!(field_i64(&layer, 0, "utm_zone"), 31);
        assert_eq!(field_i64(&layer, 0, "utm_epsg"), 32631);
    }

    /// An unclosed ring (first != last, as file readers store them) must still
    /// yield the correct area centroid — the shoelace loop wraps the closing edge.
    #[test]
    fn polygon_ring_stored_unclosed() {
        let mut l = Layer::new("poly")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("name", FieldType::Text));
        // Square centred at lon 15, lat 5 (zone 33), stored WITHOUT a closing vertex.
        let ring = Geometry::polygon(
            vec![
                Coord::xy(14.0, 4.0),
                Coord::xy(16.0, 4.0),
                Coord::xy(16.0, 6.0),
                Coord::xy(14.0, 6.0),
            ],
            vec![],
        );
        l.add_feature(Some(ring), &[("name", "sq".into())]).unwrap();
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);
        let (_out, layer) = run(json!({ "input": input }));
        // Centroid lon = 15 -> zone floor((15+180)/6)+1 = 33.
        assert_eq!(field_i64(&layer, 0, "utm_zone"), 33);
        assert_eq!(field_i64(&layer, 0, "utm_epsg"), 32633);
    }

    /// A feature with no geometry passes through with null zone/EPSG.
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
        assert_eq!(out.outputs["stamped"], json!(0));
        assert_eq!(out.outputs["skipped"], json!(1));
        let idx = layer.schema.field_index("utm_zone").unwrap();
        assert!(matches!(
            layer.features[0].attributes[idx],
            FieldValue::Null
        ));
    }

    /// Custom field names are honoured.
    #[test]
    fn custom_field_names() {
        let input = point_layer(&[("wh", -77.0365, 38.8977)]);
        let (_out, layer) = run(json!({
            "input": input, "zone_field": "zone", "epsg_field": "epsg",
        }));
        assert_eq!(field_i64(&layer, 0, "zone"), 18);
        assert_eq!(field_i64(&layer, 0, "epsg"), 32618);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CalculateUtmZoneTool.validate(&args)
        };
        // Missing input.
        assert!(bad(json!({})).is_err());
        // Same field name for both outputs.
        assert!(bad(json!({
            "input": "a.geojson", "zone_field": "z", "epsg_field": "z",
        }))
        .is_err());
        // Valid.
        assert!(bad(json!({ "input": "a.geojson" })).is_ok());
    }

    /// A projected input CRS is rejected (centroid math needs lon/lat).
    #[test]
    fn rejects_projected_crs() {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.add_feature(Some(Geometry::point(0.0, 0.0)), &[("name", "a".into())])
            .unwrap();
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({ "input": input })).unwrap();
        assert!(CalculateUtmZoneTool.run(&args, &ctx()).is_err());
    }
}
