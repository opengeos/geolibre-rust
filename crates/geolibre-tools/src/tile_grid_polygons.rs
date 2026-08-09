//! GeoLibre tool: web-map tile boundaries as a polygon layer.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Map Server Cache Tiling Scheme To
//! Polygons* (Cartography), which also covers *Create Vector Tile Index* and
//! *Generate Tile Cache Tiling Scheme*.
//!
//! ## Why the catalog needed this
//!
//! Three separate ArcGIS tools produce tile-boundary polygons and the registry
//! had none — no `tiling`, `xyz` or `quadkey` id anywhere in ~1,100 tools.
//!
//! That is conspicuous given what this crate is. GeoLibre reasons about tile
//! pyramids constantly — `raster_to_tiles`, `select_tiles_by_polygon`,
//! `write_pmtiles`, `vector_to_pmtiles`, `pmtiles_extract` — but could not
//! *draw* one, which made tiling problems hard to debug: no way to see which
//! tiles a dataset covers, check an extent against tile boundaries, or build a
//! tile index to drive a partitioned job.
//!
//! The bundled `rectangular_grid_from_*` / `hexagonal_grid_from_*` pair and
//! `grid_index_features` build arbitrary grids in the input CRS. A tile grid is
//! a different object: fixed to the WebMercator quadtree, aligned to the global
//! origin, and carrying `z`/`x`/`y` identity. An arbitrary grid cannot be made
//! to coincide with one.
//!
//! ## Two traps this has to avoid
//!
//! **Latitude clamping.** WebMercator is undefined at the poles; `tan(pi/2)`
//! diverges. Latitudes are clamped to +/-85.0511287798066 (the square-world
//! limit) before projecting, or an unclamped input silently yields garbage tile
//! indices instead of an error.
//!
//! **Tile-count overflow.** A world grid at z=22 is 2^44 tiles. The count is
//! accumulated in `u64` and checked against `max_tiles` *before* anything is
//! allocated — computing it in `usize` would overflow on 32-bit wasm and the
//! wrapped value would sail past the guard.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::args_common::{choice_or, opt_usize, usize_or};
use crate::common::load_input_raster;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// The latitude where WebMercator's y equals its x half-extent, giving the
/// square world the tile scheme assumes.
const MAX_LAT: f64 = 85.051_128_779_806_59;
/// Half the WebMercator world extent, in metres.
const HALF_EXTENT: f64 = 20_037_508.342_789_244;
/// Highest zoom accepted. 2^30 tiles per axis already exceeds any real use and
/// keeps `2^z` far inside f64's exact-integer range.
const MAX_ZOOM: usize = 30;

const SCHEMES: [&str; 2] = ["xyz", "tms"];
const OUT_CRS: [&str; 2] = ["epsg:3857", "epsg:4326"];

pub struct TileGridPolygonsTool;

impl Tool for TileGridPolygonsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "tile_grid_polygons",
            display_name: "Tile Grid Polygons",
            summary: "Emits web-map tile boundaries as polygons for a zoom level (or range) and extent, each carrying z/x/y and its quadkey (ArcGIS Map Server Cache Tiling Scheme To Polygons, Create Vector Tile Index, Generate Tile Cache Tiling Scheme). The bundled rectangular/hexagonal grid tools build arbitrary grids in the input CRS and cannot coincide with the WebMercator quadtree the PMTiles and raster_to_tiles tools use.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "zoom",
                    description: "Zoom level. Use with 'max_zoom' to emit a pyramid of levels.",
                    required: true,
                },
                ToolParamSpec {
                    name: "max_zoom",
                    description: "Optional highest zoom; when given, every level from 'zoom' to 'max_zoom' is emitted.",
                    required: false,
                },
                ToolParamSpec {
                    name: "extent",
                    description: "Optional raster or vector layer whose bounds clip the grid. Without it, the whole world at that zoom (subject to 'max_tiles').",
                    required: false,
                },
                ToolParamSpec {
                    name: "bbox",
                    description: "Optional explicit extent as 'min_lon,min_lat,max_lon,max_lat' in degrees. Takes precedence over 'extent'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "scheme",
                    description: "'xyz' (default, Google/OSM, y counts down from the north) or 'tms' (y counts up from the south).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_crs",
                    description: "'EPSG:3857' (default, square tiles) or 'EPSG:4326' (lon/lat quadrilaterals).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_tiles",
                    description: "Hard cap on the number of tiles (default 1000000); exceeding it is an error naming the computed count.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon layer. If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        let zoom = opt_usize(args, "zoom")?.ok_or_else(|| {
            ToolError::Validation("missing required integer parameter 'zoom'".to_string())
        })?;
        if zoom > MAX_ZOOM {
            return Err(ToolError::Validation(format!(
                "'zoom' must be between 0 and {MAX_ZOOM}, got {zoom}"
            )));
        }
        if let Some(mz) = opt_usize(args, "max_zoom")? {
            if mz > MAX_ZOOM {
                return Err(ToolError::Validation(format!(
                    "'max_zoom' must be between 0 and {MAX_ZOOM}, got {mz}"
                )));
            }
            if mz < zoom {
                return Err(ToolError::Validation(format!(
                    "'max_zoom' ({mz}) must be at least 'zoom' ({zoom})"
                )));
            }
        }
        choice_or(args, "scheme", &SCHEMES, "xyz")?;
        choice_or(args, "output_crs", &OUT_CRS, "epsg:3857")?;
        if let Some(b) = parse_optional_str(args, "bbox")? {
            parse_bbox(b)?;
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let zoom = opt_usize(args, "zoom")?
            .ok_or_else(|| ToolError::Validation("missing required parameter 'zoom'".into()))?;
        let max_zoom = opt_usize(args, "max_zoom")?.unwrap_or(zoom);
        let scheme = choice_or(args, "scheme", &SCHEMES, "xyz")?;
        let out_crs = choice_or(args, "output_crs", &OUT_CRS, "epsg:3857")?;
        let max_tiles = usize_or(args, "max_tiles", 1_000_000)? as u64;
        let output = parse_optional_str(args, "output")?;

        // bbox wins over extent: an explicit box is a deliberate override.
        let bounds = match parse_optional_str(args, "bbox")? {
            Some(b) => parse_bbox(b)?,
            None => match parse_optional_str(args, "extent")? {
                Some(spec) => extent_bounds(spec)?,
                None => (-180.0, -MAX_LAT, 180.0, MAX_LAT),
            },
        };
        let (min_lon, min_lat, max_lon, max_lat) = bounds;

        // Count in u64 BEFORE allocating. usize would wrap on 32-bit wasm at
        // high zooms and the wrapped value would pass the cap check.
        let mut ranges = Vec::new();
        let mut total: u64 = 0;
        for z in zoom..=max_zoom {
            let n = 1_u64 << z;
            let x0 = lon_to_tile_x(min_lon, z).min(n - 1);
            let x1 = lon_to_tile_x(max_lon, z).min(n - 1);
            // y is inverted relative to latitude, so max_lat gives the low y.
            let y0 = lat_to_tile_y(max_lat, z).min(n - 1);
            let y1 = lat_to_tile_y(min_lat, z).min(n - 1);
            let (x0, x1) = (x0.min(x1), x0.max(x1));
            let (y0, y1) = (y0.min(y1), y0.max(y1));
            let count = (x1 - x0 + 1).saturating_mul(y1 - y0 + 1);
            total = total.saturating_add(count);
            if total > max_tiles {
                return Err(ToolError::Validation(format!(
                    "the requested grid needs at least {total} tiles, over the 'max_tiles' cap of \
                     {max_tiles}; narrow the extent, lower the zoom, or raise 'max_tiles'"
                )));
            }
            ranges.push((z, x0, x1, y0, y1));
        }
        ctx.progress
            .info(&format!("zoom {zoom}-{max_zoom}: {total} tile(s)"));

        let epsg: u32 = if out_crs == "epsg:4326" { 4326 } else { 3857 };
        let mut layer = Layer::new("tile_grid")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(epsg);
        layer.add_field(FieldDef::new("z", FieldType::Integer));
        layer.add_field(FieldDef::new("x", FieldType::Integer));
        layer.add_field(FieldDef::new("y", FieldType::Integer));
        layer.add_field(FieldDef::new("quadkey", FieldType::Text));
        layer.add_field(FieldDef::new("tile_id", FieldType::Text));

        for (z, x0, x1, y0, y1) in ranges {
            let n = 1_u64 << z;
            for x in x0..=x1 {
                for y in y0..=y1 {
                    // The quadkey is always defined on XYZ order; converting a
                    // TMS y back first keeps it meaningful under both schemes.
                    let quadkey = quadkey(x, y, z);
                    // Corners come from the XYZ y, then the label is flipped
                    // for TMS. Deriving geometry from a flipped y would mirror
                    // every tile vertically.
                    let (west, north) = tile_nw(x, y, z, epsg);
                    let (east, south) = tile_nw(x + 1, y + 1, z, epsg);
                    let out_y = if scheme == "tms" { n - 1 - y } else { y };

                    layer.push(Feature {
                        fid: 0,
                        geometry: Some(Geometry::polygon(
                            vec![
                                Coord::xy(west, south),
                                Coord::xy(east, south),
                                Coord::xy(east, north),
                                Coord::xy(west, north),
                                Coord::xy(west, south),
                            ],
                            Vec::new(),
                        )),
                        attributes: vec![
                            FieldValue::Integer(z as i64),
                            FieldValue::Integer(x as i64),
                            FieldValue::Integer(out_y as i64),
                            FieldValue::Text(quadkey),
                            FieldValue::Text(format!("{z}/{x}/{out_y}")),
                        ],
                    });
                }
            }
        }

        let tile_count = layer.features.len();
        let out_path = write_or_store_layer(layer, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("tile_count".to_string(), json!(tile_count));
        outputs.insert("min_zoom".to_string(), json!(zoom));
        outputs.insert("max_zoom".to_string(), json!(max_zoom));
        outputs.insert("scheme".to_string(), json!(scheme));
        outputs.insert("output_crs".to_string(), json!(epsg));
        outputs.insert(
            "bbox".to_string(),
            json!([min_lon, min_lat, max_lon, max_lat]),
        );
        Ok(ToolRunResult { outputs })
    }
}

fn lon_to_tile_x(lon: f64, z: usize) -> u64 {
    let n = (1_u64 << z) as f64;
    let lon = lon.clamp(-180.0, 180.0);
    let x = (lon + 180.0) / 360.0 * n;
    // A longitude of exactly +180 lands on tile n, which does not exist.
    (x.floor().max(0.0) as u64).min((1_u64 << z) - 1)
}

fn lat_to_tile_y(lat: f64, z: usize) -> u64 {
    let n = (1_u64 << z) as f64;
    // Clamping first is what keeps tan() away from its pole.
    let lat = lat.clamp(-MAX_LAT, MAX_LAT).to_radians();
    let y = (1.0 - (lat.tan() + 1.0 / lat.cos()).ln() / std::f64::consts::PI) / 2.0 * n;
    (y.floor().max(0.0) as u64).min((1_u64 << z) - 1)
}

/// North-west corner of tile (x, y) at zoom z, in the requested CRS.
///
/// Accepts x/y one past the last tile so a caller can ask for the south-east
/// corner as the next tile's north-west corner.
fn tile_nw(x: u64, y: u64, z: usize, epsg: u32) -> (f64, f64) {
    let n = (1_u64 << z) as f64;
    if epsg == 4326 {
        let lon = x as f64 / n * 360.0 - 180.0;
        let t = std::f64::consts::PI * (1.0 - 2.0 * y as f64 / n);
        let lat = t.sinh().atan().to_degrees();
        (lon, lat)
    } else {
        // WebMercator tile bounds are exact linear divisions of the extent, so
        // no trigonometry is involved and no precision is lost.
        let span = 2.0 * HALF_EXTENT / n;
        (-HALF_EXTENT + x as f64 * span, HALF_EXTENT - y as f64 * span)
    }
}

/// Bing-style quadkey: interleave the x/y bits into base-4 digits, most
/// significant zoom level first.
fn quadkey(x: u64, y: u64, z: usize) -> String {
    let mut s = String::with_capacity(z);
    for i in (1..=z).rev() {
        let mask = 1_u64 << (i - 1);
        let mut digit = b'0';
        if x & mask != 0 {
            digit += 1;
        }
        if y & mask != 0 {
            digit += 2;
        }
        s.push(digit as char);
    }
    s
}

fn parse_bbox(s: &str) -> Result<(f64, f64, f64, f64), ToolError> {
    let parts: Vec<f64> = s
        .split([',', ';'])
        .map(str::trim)
        .filter(|p| !p.is_empty())
        .map(|p| {
            p.parse::<f64>().map_err(|_| {
                ToolError::Validation(format!("'bbox' entry '{p}' is not a number"))
            })
        })
        .collect::<Result<_, _>>()?;
    if parts.len() != 4 {
        return Err(ToolError::Validation(
            "'bbox' must be 'min_lon,min_lat,max_lon,max_lat'".to_string(),
        ));
    }
    let (a, b, c, d) = (parts[0], parts[1], parts[2], parts[3]);
    if a >= c || b >= d {
        return Err(ToolError::Validation(format!(
            "'bbox' is empty or inverted: min ({a}, {b}) must be strictly below max ({c}, {d})"
        )));
    }
    Ok((a, b, c, d))
}

/// Bounds of an extent layer, converted to lon/lat degrees.
///
/// A layer declaring EPSG:3857 is converted; one declaring 4326 (or nothing,
/// with degree-range coordinates) is used as-is. Anything else is refused —
/// guessing at an unknown projection is how a tile grid ends up in the wrong
/// hemisphere.
fn extent_bounds(spec: &str) -> Result<(f64, f64, f64, f64), ToolError> {
    let (min_x, min_y, max_x, max_y, epsg) = if let Ok(r) = load_input_raster(spec) {
        (
            r.x_min,
            r.y_min,
            r.x_min + r.cols as f64 * r.cell_size_x,
            r.y_min + r.rows as f64 * r.cell_size_y,
            r.crs.epsg,
        )
    } else {
        let layer = load_input_layer(spec)?;
        let mut b = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
        let mut seen = false;
        for f in layer.iter() {
            if let Some(g) = f.geometry.as_ref() {
                for c in geometry_coords(g) {
                    seen = true;
                    b.0 = b.0.min(c.0);
                    b.1 = b.1.min(c.1);
                    b.2 = b.2.max(c.0);
                    b.3 = b.3.max(c.1);
                }
            }
        }
        if !seen {
            return Err(ToolError::Execution(
                "'extent' layer has no geometry to take bounds from".to_string(),
            ));
        }
        (b.0, b.1, b.2, b.3, layer.crs_epsg())
    };

    match epsg {
        Some(3857) => Ok((
            mercator_x_to_lon(min_x),
            mercator_y_to_lat(min_y),
            mercator_x_to_lon(max_x),
            mercator_y_to_lat(max_y),
        )),
        Some(4326) | None => {
            if !(-180.0..=180.0).contains(&min_x) || !(-90.0..=90.0).contains(&min_y) {
                return Err(ToolError::Validation(format!(
                    "'extent' has coordinates ({min_x:.3}, {min_y:.3}) outside the lon/lat range \
                     but declares no usable CRS; reproject it to EPSG:4326 or EPSG:3857 first"
                )));
            }
            Ok((min_x, min_y, max_x, max_y))
        }
        Some(other) => Err(ToolError::Validation(format!(
            "'extent' is EPSG:{other}; tile grids are defined on WebMercator, so supply an extent \
             in EPSG:4326 or EPSG:3857 (or pass 'bbox' in degrees)"
        ))),
    }
}

fn mercator_x_to_lon(x: f64) -> f64 {
    (x / HALF_EXTENT * 180.0).clamp(-180.0, 180.0)
}

fn mercator_y_to_lat(y: f64) -> f64 {
    let t = y / HALF_EXTENT * std::f64::consts::PI;
    t.sinh().atan().to_degrees().clamp(-MAX_LAT, MAX_LAT)
}

fn geometry_coords(g: &Geometry) -> Vec<(f64, f64)> {
    match g {
        Geometry::Point(p) => vec![(p.x, p.y)],
        Geometry::MultiPoint(ps) => ps.iter().map(|p| (p.x, p.y)).collect(),
        Geometry::LineString(cs) => cs.iter().map(|c| (c.x, c.y)).collect(),
        Geometry::MultiLineString(ls) => {
            ls.iter().flat_map(|l| l.iter().map(|c| (c.x, c.y))).collect()
        }
        Geometry::Polygon { exterior, .. } => exterior.0.iter().map(|c| (c.x, c.y)).collect(),
        Geometry::MultiPolygon(ps) => ps
            .iter()
            .flat_map(|(e, _)| e.0.iter().map(|c| (c.x, c.y)))
            .collect(),
        Geometry::GeometryCollection(gs) => gs.iter().flat_map(geometry_coords).collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use wbcore::{AllowAllCapabilities, ProgressSink};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = TileGridPolygonsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    fn attr(layer: &Layer, i: usize, name: &str) -> FieldValue {
        let k = layer.schema.field_index(name).unwrap();
        layer.features[i].attributes[k].clone()
    }

    #[test]
    fn zoom_zero_is_a_single_world_tile() {
        let (layer, res) = run(json!({"zoom": 0}));
        assert_eq!(layer.features.len(), 1);
        assert_eq!(res.outputs["tile_count"], json!(1));
        assert_eq!(attr(&layer, 0, "tile_id"), FieldValue::Text("0/0/0".into()));
    }

    #[test]
    fn zoom_two_covers_the_world_in_sixteen_tiles() {
        let (layer, _) = run(json!({"zoom": 2}));
        assert_eq!(layer.features.len(), 16);
    }

    #[test]
    fn the_zoom_zero_tile_spans_the_whole_mercator_extent() {
        let (layer, _) = run(json!({"zoom": 0}));
        let Some(Geometry::Polygon { exterior, .. }) = layer.features[0].geometry.as_ref() else {
            panic!("expected a polygon");
        };
        let xs: Vec<f64> = exterior.0.iter().map(|c| c.x).collect();
        let ys: Vec<f64> = exterior.0.iter().map(|c| c.y).collect();
        let min_x = xs.iter().cloned().fold(f64::MAX, f64::min);
        let max_y = ys.iter().cloned().fold(f64::MIN, f64::max);
        assert!((min_x + HALF_EXTENT).abs() < 1e-6, "got {min_x}");
        assert!((max_y - HALF_EXTENT).abs() < 1e-6, "got {max_y}");
    }

    #[test]
    fn epsg_4326_output_spans_the_full_longitude_range() {
        let (layer, res) = run(json!({"zoom": 0, "output_crs": "EPSG:4326"}));
        assert_eq!(res.outputs["output_crs"], json!(4326));
        let Some(Geometry::Polygon { exterior, .. }) = layer.features[0].geometry.as_ref() else {
            panic!()
        };
        let xs: Vec<f64> = exterior.0.iter().map(|c| c.x).collect();
        assert!((xs.iter().cloned().fold(f64::MAX, f64::min) + 180.0).abs() < 1e-9);
        assert!((xs.iter().cloned().fold(f64::MIN, f64::max) - 180.0).abs() < 1e-9);
        // Latitude tops out at the WebMercator limit, not at 90.
        let ys: Vec<f64> = exterior.0.iter().map(|c| c.y).collect();
        let top = ys.iter().cloned().fold(f64::MIN, f64::max);
        assert!((top - MAX_LAT).abs() < 1e-6, "got {top}");
    }

    #[test]
    fn quadkeys_match_the_bing_reference_values() {
        // The canonical z=1 quadrant labelling: NW=0, NE=1, SW=2, SE=3.
        assert_eq!(quadkey(0, 0, 1), "0");
        assert_eq!(quadkey(1, 0, 1), "1");
        assert_eq!(quadkey(0, 1, 1), "2");
        assert_eq!(quadkey(1, 1, 1), "3");
        // A published multi-level example.
        assert_eq!(quadkey(3, 5, 3), "213");
        assert_eq!(quadkey(0, 0, 0), "");
    }

    #[test]
    fn tms_flips_y_but_not_the_geometry() {
        // The trap: deriving corners from a flipped y would mirror every tile.
        // Same tile under both schemes must occupy the same ground.
        let (xyz, _) = run(json!({"zoom": 1, "scheme": "xyz"}));
        let (tms, _) = run(json!({"zoom": 1, "scheme": "tms"}));
        let find = |l: &Layer, x: i64, y: i64| -> Vec<(f64, f64)> {
            let xi = l.schema.field_index("x").unwrap();
            let yi = l.schema.field_index("y").unwrap();
            for f in l.iter() {
                if f.attributes[xi] == FieldValue::Integer(x)
                    && f.attributes[yi] == FieldValue::Integer(y)
                {
                    let Some(Geometry::Polygon { exterior, .. }) = f.geometry.as_ref() else {
                        panic!()
                    };
                    return exterior.0.iter().map(|c| (c.x, c.y)).collect();
                }
            }
            panic!("tile {x},{y} not found");
        };
        // XYZ (0,0) is the north-west tile; in TMS the same ground is (0,1).
        assert_eq!(find(&xyz, 0, 0), find(&tms, 0, 1));
    }

    #[test]
    fn a_bbox_narrows_the_grid_to_the_covering_tiles() {
        // A small box near 0,0 at z=2 straddles no more than a couple of tiles.
        let (layer, _) = run(json!({"zoom": 2, "bbox": "1,1,2,2"}));
        assert_eq!(layer.features.len(), 1);
        // East of the prime meridian and north of the equator at z=2 is x=2,y=1.
        assert_eq!(attr(&layer, 0, "x"), FieldValue::Integer(2));
        assert_eq!(attr(&layer, 0, "y"), FieldValue::Integer(1));
    }

    #[test]
    fn a_zoom_range_emits_every_level() {
        let (layer, res) = run(json!({"zoom": 0, "max_zoom": 2}));
        // 1 + 4 + 16
        assert_eq!(layer.features.len(), 21);
        assert_eq!(res.outputs["min_zoom"], json!(0));
        assert_eq!(res.outputs["max_zoom"], json!(2));
    }

    #[test]
    fn an_oversized_request_is_refused_before_allocating() {
        // z=20 world-wide is ~1.1e12 tiles. This must fail fast with a count,
        // not attempt the allocation.
        let args: ToolArgs = serde_json::from_value(json!({"zoom": 20})).unwrap();
        let err = TileGridPolygonsTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("max_tiles"), "{err}");
    }

    #[test]
    fn a_high_zoom_tile_count_does_not_overflow_the_guard() {
        // z=30 world-wide is 2^60 tiles, which wraps a 32-bit usize to a small
        // number. The u64 accumulation is what keeps the cap meaningful.
        let args: ToolArgs =
            serde_json::from_value(json!({"zoom": 30, "max_tiles": 1000000})).unwrap();
        assert!(TileGridPolygonsTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn a_polar_latitude_is_clamped_rather_than_producing_garbage() {
        // tan(pi/2) diverges; without clamping the tile index is nonsense.
        let (layer, _) = run(json!({"zoom": 1, "bbox": "-10,-89.9,10,89.9"}));
        let yi = layer.schema.field_index("y").unwrap();
        for f in layer.iter() {
            let FieldValue::Integer(y) = f.attributes[yi] else {
                panic!()
            };
            assert!((0..2).contains(&y), "y {y} out of range at z=1");
        }
    }

    #[test]
    fn an_extent_layer_in_3857_is_converted_before_tiling() {
        // Roughly a degree box around the origin, expressed in metres.
        let mut l = Layer::new("e")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_feature(Some(Geometry::point(111_319.5, 111_325.1)), &[])
            .unwrap();
        l.add_feature(Some(Geometry::point(222_639.0, 222_684.2)), &[])
            .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let path = wbvector::memory_store::make_vector_memory_path(&id);
        let (_, res) = run(json!({"zoom": 4, "extent": path}));
        let bbox = res.outputs["bbox"].as_array().unwrap();
        // 111319.5 m easting is 1 degree of longitude.
        assert!((bbox[0].as_f64().unwrap() - 1.0).abs() < 1e-3, "{bbox:?}");
    }

    #[test]
    fn an_extent_in_an_unsupported_crs_is_refused() {
        let mut l = Layer::new("e")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(32610);
        l.add_feature(Some(Geometry::point(550000.0, 4180000.0)), &[])
            .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let path = wbvector::memory_store::make_vector_memory_path(&id);
        let args: ToolArgs =
            serde_json::from_value(json!({"zoom": 4, "extent": path})).unwrap();
        let err = TileGridPolygonsTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("EPSG:32610"), "{err}");
    }

    #[test]
    fn tiles_are_contiguous_with_no_gap_between_neighbours() {
        // A gap or overlap here would mean raster_to_tiles and this tool
        // disagree about where a tile boundary is.
        let (layer, _) = run(json!({"zoom": 2, "bbox": "-179,-60,179,60"}));
        let xi = layer.schema.field_index("x").unwrap();
        let yi = layer.schema.field_index("y").unwrap();
        let bounds = |x: i64, y: i64| -> (f64, f64, f64, f64) {
            for f in layer.iter() {
                if f.attributes[xi] == FieldValue::Integer(x)
                    && f.attributes[yi] == FieldValue::Integer(y)
                {
                    let Some(Geometry::Polygon { exterior, .. }) = f.geometry.as_ref() else {
                        panic!()
                    };
                    let xs: Vec<f64> = exterior.0.iter().map(|c| c.x).collect();
                    let ys: Vec<f64> = exterior.0.iter().map(|c| c.y).collect();
                    return (
                        xs.iter().cloned().fold(f64::MAX, f64::min),
                        ys.iter().cloned().fold(f64::MAX, f64::min),
                        xs.iter().cloned().fold(f64::MIN, f64::max),
                        ys.iter().cloned().fold(f64::MIN, f64::max),
                    );
                }
            }
            panic!("missing tile {x},{y}");
        };
        let a = bounds(1, 1);
        let b = bounds(2, 1);
        assert!((a.2 - b.0).abs() < 1e-6, "east of (1,1) != west of (2,1)");
        let c = bounds(1, 2);
        assert!((a.1 - c.3).abs() < 1e-6, "south of (1,1) != north of (1,2)");
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            TileGridPolygonsTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"zoom": 40})));
        assert!(bad(json!({"zoom": 5, "max_zoom": 3})));
        assert!(bad(json!({"zoom": 2, "scheme": "wmts"})));
        assert!(bad(json!({"zoom": 2, "output_crs": "EPSG:4269"})));
        assert!(bad(json!({"zoom": 2, "bbox": "1,2,3"})));
        // Inverted box: silently swapping it would hide a caller's bug.
        assert!(bad(json!({"zoom": 2, "bbox": "10,10,0,0"})));
    }
}
