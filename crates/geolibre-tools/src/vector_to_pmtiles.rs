//! GeoLibre tool: pack a vector layer into a single PMTiles archive.
//!
//! The vector counterpart to `write_pmtiles`: instead of a PNG raster pyramid
//! it builds a Mapbox Vector Tile (MVT) pyramid, so a whole styleable overlay
//! is one download. Tiling (clip, simplify, quantize, MVT encode) is delegated
//! to `freestiler-core`; this module adapts `wbvector::Layer` to that engine's
//! input types and hands the resulting tiles to the shared PMTiles writer.

use std::collections::BTreeMap;

use freestiler_core::engine::{compute_all_bounds, generate_tiles, SilentReporter, TileConfig};
use freestiler_core::pmtiles_writer::TileFormat;
use freestiler_core::tiler::{
    Feature as MvtFeature, Geometry as MvtGeometry, LayerData, PropertyValue,
};
use geo::{Coord as GeoCoord, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata,
    ToolParamSpec, ToolRunResult,
};
use wbvector::{Coord, FieldType, FieldValue, Geometry, Layer, Ring};

use crate::common::write_bytes;
use crate::pmtiles::{self, LonLatBounds, Tile};
use crate::vector_common::{load_input_layer, parse_optional_str};

/// The tiling engine works in lon/lat, so anything else is reprojected first.
const WGS84_EPSG: u32 = 4326;

/// Beyond this, tile counts explode while MVT detail is already at the 4096
/// quantization grid's limit.
const MAX_ZOOM: u64 = 18;

/// Web Mercator's latitude limit.
const MERCATOR_LAT_LIMIT: f64 = 85.051_128_779_806_59;

/// Packs a vector layer into a single PMTiles archive of MVT tiles.
pub struct VectorToPmTilesTool;

impl Tool for VectorToPmTilesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "vector_to_pmtiles",
            display_name: "Vector to PMTiles",
            summary: "Pack a vector layer into a single PMTiles archive (Mapbox Vector Tile pyramid).",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec { name: "input", description: "Input vector file path (GeoJSON, Shapefile, GeoPackage, FlatGeobuf, GeoParquet, ...).", required: true },
                ToolParamSpec { name: "output", description: "Output PMTiles file path (e.g. /work/roads.pmtiles).", required: true },
                ToolParamSpec { name: "min_zoom", description: "Minimum zoom level (default 0).", required: false },
                ToolParamSpec { name: "max_zoom", description: "Maximum zoom level (default 14, maximum 18).", required: false },
                ToolParamSpec { name: "layer_name", description: "Name of the layer inside the tiles, used when styling (default: the input layer's name).", required: false },
                ToolParamSpec { name: "simplify", description: "Simplify geometries per zoom level (default true).", required: false },
                ToolParamSpec { name: "drop_rate", description: "Rate at which features are dropped at low zooms, as in tippecanoe (default: no dropping).", required: false },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["input", "output"] {
            if args.get(key).and_then(Value::as_str).is_none() {
                return Err(ToolError::Validation(format!(
                    "missing required string parameter '{key}'"
                )));
            }
        }
        let (min_zoom, max_zoom) = zoom_range(args)?;
        if min_zoom > max_zoom {
            return Err(ToolError::Validation(format!(
                "min_zoom ({min_zoom}) must not exceed max_zoom ({max_zoom})"
            )));
        }
        if let Some(rate) = args.get("drop_rate") {
            let rate = rate
                .as_f64()
                .ok_or_else(|| ToolError::Validation("drop_rate must be a number".to_string()))?;
            if !(rate.is_finite() && rate > 0.0) {
                return Err(ToolError::Validation(
                    "drop_rate must be a finite number greater than 0".to_string(),
                ));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = require_str(args, "output")?;
        let (min_zoom, max_zoom) = zoom_range(args)?;
        let simplify = args.get("simplify").and_then(Value::as_bool).unwrap_or(true);
        let drop_rate = args.get("drop_rate").and_then(Value::as_f64);

        let layer = load_input_layer(input)?;
        if layer.features.is_empty() {
            return Err(ToolError::Validation(
                "input layer has no features to tile".to_string(),
            ));
        }

        let layer = to_wgs84(layer, ctx)?;
        let name = match parse_optional_str(args, "layer_name")? {
            Some(n) => n.to_string(),
            None if !layer.name.trim().is_empty() => layer.name.clone(),
            None => "layer".to_string(),
        };

        let field_names: Vec<String> = layer
            .schema
            .fields()
            .iter()
            .map(|f| f.name.clone())
            .collect();
        let field_types: Vec<String> = layer
            .schema
            .fields()
            .iter()
            .map(|f| tilejson_field_type(&f.field_type).to_string())
            .collect();

        ctx.progress.info("converting geometries");
        let features = convert_features(&layer);
        if features.is_empty() {
            return Err(ToolError::Execution(
                "no tileable geometries found (all features were empty or unsupported)".to_string(),
            ));
        }
        let feature_count = features.len();

        let data = LayerData {
            name: name.clone(),
            features,
            prop_names: field_names.clone(),
            prop_types: field_types
                .iter()
                .map(|t| t.to_lowercase())
                .collect(),
            min_zoom: min_zoom as u8,
            max_zoom: max_zoom as u8,
        };
        let layers = [data];

        let (west, south, east, north) = compute_all_bounds(&layers);
        if ![west, south, east, north].iter().all(|v| v.is_finite()) {
            return Err(ToolError::Execution(
                "input geometries have no finite extent".to_string(),
            ));
        }

        let config = TileConfig {
            tile_format: TileFormat::Mvt,
            min_zoom: min_zoom as u8,
            max_zoom: max_zoom as u8,
            base_zoom: None,
            simplification: simplify,
            drop_rate,
            cluster_distance: None,
            cluster_maxzoom: None,
            coalesce: false,
        };

        ctx.progress.info("generating vector tiles");
        let generated = generate_tiles(&layers, &config, &SilentReporter)
            .map_err(|e| ToolError::Execution(format!("vector tiling failed: {e}")))?;
        if generated.is_empty() {
            return Err(ToolError::Execution(
                "tiling produced no tiles for the requested zoom range".to_string(),
            ));
        }
        let tile_count = generated.len();

        let tiles: Vec<Tile> = generated
            .into_iter()
            .map(|(c, data)| Tile { z: c.z, x: c.x, y: c.y, data })
            .collect();
        let bounds = clamp_bounds(west, south, east, north);

        ctx.progress.info("packing PMTiles archive");
        let metadata = tilejson(&name, &field_names, &field_types, &bounds, min_zoom, max_zoom);
        let archive = pmtiles::build_vector(
            tiles,
            &bounds,
            min_zoom as u8,
            max_zoom as u8,
            metadata.to_string().as_bytes(),
        )?;
        write_bytes(output, &archive)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(output));
        outputs.insert("layer_name".to_string(), json!(name));
        outputs.insert("min_zoom".to_string(), json!(min_zoom));
        outputs.insert("max_zoom".to_string(), json!(max_zoom));
        outputs.insert("features".to_string(), json!(feature_count));
        outputs.insert("tiles".to_string(), json!(tile_count));
        outputs.insert("bytes".to_string(), json!(archive.len()));
        Ok(ToolRunResult { outputs })
    }
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

/// Clamps a raw data extent to what a PMTiles header and TileJSON can express.
/// Source data legitimately carries out-of-range longitudes — datasets that
/// unwrap across the antimeridian rather than splitting at it are common, e.g.
/// Alaska reaching past -180 — and the tiler clips such geometry at the
/// projection edge. The tiled area therefore never exceeds these limits even
/// when the input extent does, so reporting the raw extent would advertise
/// coverage that no tile in the archive actually provides.
fn clamp_bounds(west: f64, south: f64, east: f64, north: f64) -> LonLatBounds {
    LonLatBounds {
        min_lon: west.clamp(-180.0, 180.0),
        min_lat: south.clamp(-MERCATOR_LAT_LIMIT, MERCATOR_LAT_LIMIT),
        max_lon: east.clamp(-180.0, 180.0),
        max_lat: north.clamp(-MERCATOR_LAT_LIMIT, MERCATOR_LAT_LIMIT),
    }
}

fn zoom_range(args: &ToolArgs) -> Result<(u64, u64), ToolError> {
    let min_zoom = args.get("min_zoom").and_then(Value::as_u64).unwrap_or(0);
    let max_zoom = args.get("max_zoom").and_then(Value::as_u64).unwrap_or(14);
    if max_zoom > MAX_ZOOM {
        return Err(ToolError::Validation(format!(
            "max_zoom must be <= {MAX_ZOOM}"
        )));
    }
    Ok((min_zoom, max_zoom))
}

/// Reprojects `layer` to WGS84 unless it is already there. A layer with no
/// declared CRS is taken at face value as lon/lat, which is what the formats
/// that omit it (notably GeoJSON, whose spec fixes CRS84) mean by it.
fn to_wgs84(layer: Layer, ctx: &ToolContext) -> Result<Layer, ToolError> {
    match layer.crs.as_ref().and_then(|c| c.epsg) {
        Some(WGS84_EPSG) => Ok(layer),
        Some(_) => {
            ctx.progress.info("reprojecting to EPSG:4326");
            wbvector::reproject::layer_to_epsg(&layer, WGS84_EPSG)
                .map_err(|e| ToolError::Execution(format!("reprojection to 4326 failed: {e}")))
        }
        None => {
            ctx.progress
                .info("input has no declared CRS; assuming EPSG:4326 (lon/lat)");
            Ok(layer)
        }
    }
}

/// Flattens `layer`'s features into the tiling engine's representation.
/// Geometry collections contribute one feature per member, each carrying a copy
/// of the parent's attributes; empty and null geometries are dropped, since the
/// MVT encoder has nothing to emit for them.
fn convert_features(layer: &Layer) -> Vec<MvtFeature> {
    let mut out = Vec::with_capacity(layer.features.len());
    for feature in &layer.features {
        let Some(geometry) = feature.geometry.as_ref() else {
            continue;
        };
        let mut geoms = Vec::new();
        flatten_geometry(geometry, &mut geoms);
        let properties: Vec<PropertyValue> =
            feature.attributes.iter().map(convert_value).collect();
        for geom in geoms {
            out.push(MvtFeature {
                id: Some(feature.fid),
                geometry: geom,
                properties: properties.clone(),
            });
        }
    }
    out
}

fn flatten_geometry(geom: &Geometry, out: &mut Vec<MvtGeometry>) {
    match geom {
        Geometry::Point(c) => out.push(MvtGeometry::Point(Point::new(c.x, c.y))),
        Geometry::MultiPoint(cs) if !cs.is_empty() => out.push(MvtGeometry::MultiPoint(
            MultiPoint(cs.iter().map(|c| Point::new(c.x, c.y)).collect()),
        )),
        Geometry::LineString(cs) if cs.len() >= 2 => {
            out.push(MvtGeometry::LineString(to_line_string(cs)))
        }
        Geometry::MultiLineString(lines) => {
            let lines: Vec<LineString<f64>> = lines
                .iter()
                .filter(|cs| cs.len() >= 2)
                .map(|cs| to_line_string(cs))
                .collect();
            if !lines.is_empty() {
                out.push(MvtGeometry::MultiLineString(MultiLineString(lines)));
            }
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            if let Some(poly) = to_polygon(exterior, interiors) {
                out.push(MvtGeometry::Polygon(poly));
            }
        }
        Geometry::MultiPolygon(polys) => {
            let polys: Vec<Polygon<f64>> = polys
                .iter()
                .filter_map(|(exterior, interiors)| to_polygon(exterior, interiors))
                .collect();
            if !polys.is_empty() {
                out.push(MvtGeometry::MultiPolygon(MultiPolygon(polys)));
            }
        }
        Geometry::GeometryCollection(members) => {
            for member in members {
                flatten_geometry(member, out);
            }
        }
        _ => {}
    }
}

fn to_line_string(coords: &[Coord]) -> LineString<f64> {
    LineString(
        coords
            .iter()
            .map(|c| GeoCoord { x: c.x, y: c.y })
            .collect(),
    )
}

/// A ring needs 4 positions to bound any area (3 distinct plus the repeated
/// closing vertex); anything shorter is degenerate and yields no polygon.
fn to_polygon(exterior: &Ring, interiors: &[Ring]) -> Option<Polygon<f64>> {
    if exterior.0.len() < 4 {
        return None;
    }
    Some(Polygon::new(
        to_line_string(&exterior.0),
        interiors
            .iter()
            .filter(|r| r.0.len() >= 4)
            .map(|r| to_line_string(&r.0))
            .collect(),
    ))
}

fn convert_value(value: &FieldValue) -> PropertyValue {
    match value {
        FieldValue::Integer(v) => PropertyValue::Int(*v),
        FieldValue::Float(v) => PropertyValue::Double(*v),
        FieldValue::Text(v) => PropertyValue::String(v.clone()),
        FieldValue::Boolean(v) => PropertyValue::Bool(*v),
        FieldValue::Date(v) | FieldValue::DateTime(v) => PropertyValue::String(v.clone()),
        // MVT has no binary property type, and inlining blobs into every tile
        // that touches the feature would bloat the archive regardless.
        FieldValue::Blob(_) | FieldValue::Null => PropertyValue::Null,
    }
}

/// TileJSON's `fields` vocabulary, which is narrower than the source schema's.
fn tilejson_field_type(field_type: &FieldType) -> &'static str {
    match field_type {
        FieldType::Integer | FieldType::Float => "Number",
        FieldType::Boolean => "Boolean",
        _ => "String",
    }
}

/// Builds the TileJSON metadata block. The `vector_layers` array is what lets a
/// renderer discover the layer id and its fields; MapLibre will not style a
/// vector source without it.
fn tilejson(
    name: &str,
    field_names: &[String],
    field_types: &[String],
    bounds: &LonLatBounds,
    min_zoom: u64,
    max_zoom: u64,
) -> Value {
    let fields: serde_json::Map<String, Value> = field_names
        .iter()
        .zip(field_types)
        .map(|(n, t)| (n.clone(), json!(t)))
        .collect();
    json!({
        "name": name,
        "type": "overlay",
        "format": "pbf",
        "minzoom": min_zoom,
        "maxzoom": max_zoom,
        "bounds": [bounds.min_lon, bounds.min_lat, bounds.max_lon, bounds.max_lat],
        "vector_layers": [{
            "id": name,
            "fields": fields,
            "minzoom": min_zoom,
            "maxzoom": max_zoom,
        }],
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring(coords: &[(f64, f64)]) -> Ring {
        Ring(
            coords
                .iter()
                .map(|(x, y)| Coord { x: *x, y: *y, z: None, m: None })
                .collect(),
        )
    }

    #[test]
    fn flattens_geometry_collections_into_members() {
        let gc = Geometry::GeometryCollection(vec![
            Geometry::Point(Coord { x: 1.0, y: 2.0, z: None, m: None }),
            Geometry::GeometryCollection(vec![Geometry::Point(Coord {
                x: 3.0,
                y: 4.0,
                z: None,
                m: None,
            })]),
        ]);
        let mut out = Vec::new();
        flatten_geometry(&gc, &mut out);
        assert_eq!(out.len(), 2, "nested collections should flatten recursively");
    }

    #[test]
    fn drops_degenerate_geometries() {
        let mut out = Vec::new();
        // A single-vertex line and a 3-position ring cannot be rendered.
        flatten_geometry(
            &Geometry::LineString(vec![Coord { x: 0.0, y: 0.0, z: None, m: None }]),
            &mut out,
        );
        flatten_geometry(
            &Geometry::Polygon {
                exterior: ring(&[(0.0, 0.0), (1.0, 0.0), (0.0, 0.0)]),
                interiors: vec![],
            },
            &mut out,
        );
        assert!(out.is_empty(), "degenerate geometries should be dropped");
    }

    #[test]
    fn keeps_polygon_but_drops_degenerate_hole() {
        let mut out = Vec::new();
        flatten_geometry(
            &Geometry::Polygon {
                exterior: ring(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)]),
                interiors: vec![ring(&[(0.5, 0.5), (1.0, 0.5), (0.5, 0.5)])],
            },
            &mut out,
        );
        match out.as_slice() {
            [MvtGeometry::Polygon(p)] => assert!(p.interiors().is_empty()),
            _ => panic!("expected a single polygon with its degenerate hole dropped"),
        }
    }

    #[test]
    fn tilejson_declares_vector_layers_for_styling() {
        let meta = tilejson(
            "roads",
            &["name".to_string(), "lanes".to_string()],
            &["String".to_string(), "Number".to_string()],
            &LonLatBounds { min_lon: -1.0, min_lat: -2.0, max_lon: 3.0, max_lat: 4.0 },
            0,
            14,
        );
        let layer = &meta["vector_layers"][0];
        assert_eq!(layer["id"], "roads");
        assert_eq!(layer["fields"]["lanes"], "Number");
        assert_eq!(meta["format"], "pbf");
    }

    #[test]
    fn bounds_are_clamped_to_the_tiled_area() {
        // The US states dataset unwraps Alaska past the antimeridian to
        // -188.9; the tiler clips there, so the header must not claim it.
        let b = clamp_bounds(-188.905, 17.93, -65.627, 71.35);
        assert_eq!(b.min_lon, -180.0);
        assert_eq!(b.max_lon, -65.627);
        assert_eq!(b.min_lat, 17.93);

        let polar = clamp_bounds(-10.0, -89.0, 10.0, 89.0);
        assert!(polar.min_lat > -85.06 && polar.max_lat < 85.06);
    }

    #[test]
    fn in_range_bounds_pass_through_untouched() {
        let b = clamp_bounds(-1.5, -2.5, 3.5, 4.5);
        assert_eq!((b.min_lon, b.min_lat, b.max_lon, b.max_lat), (-1.5, -2.5, 3.5, 4.5));
    }

    #[test]
    fn blobs_become_null_to_keep_attributes_aligned() {
        let values = [
            FieldValue::Text("a".to_string()),
            FieldValue::Blob(vec![1, 2, 3]),
            FieldValue::Integer(7),
        ];
        let converted: Vec<PropertyValue> = values.iter().map(convert_value).collect();
        assert_eq!(converted.len(), 3, "positions must line up with prop_names");
        assert!(matches!(converted[1], PropertyValue::Null));
    }
}
