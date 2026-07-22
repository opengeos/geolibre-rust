//! GeoLibre tool: build a map-book index grid with page naming.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Grid Index Features* (Cartography),
//! with a `mode=strip` variant that mirrors *Strip Map Index Features*. The
//! bundled `rectangular_grid_from_*` tools only draw a bare fishnet — no page
//! semantics (scale-derived sizing, page naming, or a data-intersection filter),
//! and there is no equivalent at all for the rotated strip-map pages that follow
//! a route. Map-book / atlas indexing is a common cartographic-output need.
//!
//! **Grid mode** (default). The tool tiles a rectangular extent into page-sized
//! rectangles and emits one polygon per page carrying `page_name`, `row`, `col`,
//! and `page_number`:
//!
//! - The extent comes from explicit `x_min`/`y_min`/`x_max`/`y_max`, or (when
//!   those are absent) from the bounding box of the `input` layer.
//! - Tile size in CRS units comes from explicit `tile_width`/`tile_height`, or
//!   is derived from a paper `page_size` preset (A0..A4, letter, legal, tabloid)
//!   scaled by the map-scale denominator `map_scale` (ground size = paper metres
//!   × scale; the CRS is assumed to be metric).
//! - The grid is aligned to an `origin` (defaults to the extent's lower-left
//!   corner). Columns run west→east, rows are numbered north→south so the
//!   top-left page is row 1 / column 1.
//! - `naming` = `alphanumeric` (default) names pages by column letter + row
//!   number (`A1`, `B3`, ... spreadsheet-style past `Z`), or `sequential` names
//!   them by their reading-order page number.
//! - `intersect_only` drops pages that do not intersect the `input` features
//!   (fast bbox prefilter, then a precise `geo` intersection test).
//!
//! **Strip mode** (`mode=strip`). Walks a `route` polyline (defaulting to the
//! first line in `input`) and places overlapping rectangular pages centred on
//! the route and rotated to the local bearing, so a long linear feature (a
//! highway, a pipeline, a river) is covered by a sequence of oriented map pages.
//! Page length runs along the route (`tile_height`), width across it
//! (`tile_width`), and consecutive pages overlap by the `overlap` fraction.

use std::collections::BTreeMap;

use geo::{Coord as GeoCoord, Intersects, LineString, MultiLineString, MultiPoint, Point, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct GridIndexFeaturesTool;

impl Tool for GridIndexFeaturesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "grid_index_features",
            display_name: "Grid Index Features",
            summary: "Build a map-series index grid of page-sized rectangles with page names (A1, B3, ... or sequential), scale-derived sizing, and an optional data-intersection filter — like ArcGIS Grid Index Features, plus a strip mode that follows a route with rotated pages like ArcGIS Strip Map Index Features.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer. Supplies the extent when no explicit extent is given, the route for strip mode, and the features tested by intersect_only.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon vector of index pages (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "mode",
                    description: "grid (default): tile a rectangular extent. strip: follow a route line with rotated overlapping pages.",
                    required: false,
                },
                ToolParamSpec {
                    name: "x_min",
                    description: "Explicit extent minimum X (CRS units). Give all four of x_min/y_min/x_max/y_max to override the input extent.",
                    required: false,
                },
                ToolParamSpec {
                    name: "y_min",
                    description: "Explicit extent minimum Y (CRS units).",
                    required: false,
                },
                ToolParamSpec {
                    name: "x_max",
                    description: "Explicit extent maximum X (CRS units).",
                    required: false,
                },
                ToolParamSpec {
                    name: "y_max",
                    description: "Explicit extent maximum Y (CRS units).",
                    required: false,
                },
                ToolParamSpec {
                    name: "tile_width",
                    description: "Page width in CRS units. Give both tile_width and tile_height, or use page_size + map_scale.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tile_height",
                    description: "Page height in CRS units (page length along the route in strip mode).",
                    required: false,
                },
                ToolParamSpec {
                    name: "page_size",
                    description: "Paper page preset (a0..a4, letter, legal, tabloid) whose metric dimensions are multiplied by map_scale to size each page. Requires map_scale.",
                    required: false,
                },
                ToolParamSpec {
                    name: "map_scale",
                    description: "Map-scale denominator (e.g. 24000 for 1:24,000). With page_size, tile size = paper metres × map_scale (metric CRS assumed).",
                    required: false,
                },
                ToolParamSpec {
                    name: "origin_x",
                    description: "X of the grid origin the tile boundaries align to. Default: extent minimum X.",
                    required: false,
                },
                ToolParamSpec {
                    name: "origin_y",
                    description: "Y of the grid origin the tile boundaries align to. Default: extent minimum Y.",
                    required: false,
                },
                ToolParamSpec {
                    name: "naming",
                    description: "Page naming scheme: alphanumeric (default; column letter + row number, e.g. A1) or sequential (reading-order number).",
                    required: false,
                },
                ToolParamSpec {
                    name: "intersect_only",
                    description: "When true, drop pages that do not intersect any input feature. Requires input. Default false.",
                    required: false,
                },
                ToolParamSpec {
                    name: "route",
                    description: "Strip mode: line vector layer to follow. Defaults to the first line feature in input.",
                    required: false,
                },
                ToolParamSpec {
                    name: "overlap",
                    description: "Strip mode: fraction [0,1) by which consecutive pages overlap along the route. Default 0.05.",
                    required: false,
                },
                ToolParamSpec {
                    name: "epsg",
                    description: "Output EPSG code when the extent is given explicitly and there is no input layer to inherit a CRS from.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        // Load the input layer once when we have a path (extent / route / filter).
        let input_layer = match prm.input.as_deref() {
            Some(p) => Some(load_input_layer(p)?),
            None => None,
        };

        // Resolve the output CRS: input's CRS, else the explicit epsg.
        let epsg = input_layer.as_ref().and_then(|l| l.crs_epsg()).or(prm.epsg);

        let (layer, page_count, mode_label) = match prm.mode {
            Mode::Grid => build_grid(&prm, input_layer.as_ref(), epsg, ctx)?,
            Mode::Strip => build_strip(&prm, input_layer.as_ref(), epsg, ctx)?,
        };

        let out_path = write_or_store_layer(layer, output)?;

        ctx.progress
            .info(&format!("{page_count} index page(s) ({mode_label})"));

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("page_count".to_string(), json!(page_count));
        outputs.insert("mode".to_string(), json!(mode_label));
        Ok(ToolRunResult { outputs })
    }
}

// ── Grid mode ────────────────────────────────────────────────────────────────

/// Builds the rectangular page grid, returning the layer, page count and label.
fn build_grid(
    prm: &Params,
    input: Option<&Layer>,
    epsg: Option<u32>,
    ctx: &ToolContext,
) -> Result<(Layer, usize, &'static str), ToolError> {
    let extent = resolve_extent(prm, input)?;
    let (tw, th) = resolve_tile_size(prm)?;

    let origin_x = prm.origin_x.unwrap_or(extent.x_min);
    let origin_y = prm.origin_y.unwrap_or(extent.y_min);

    // Column/row index ranges (aligned to the origin) that cover the extent.
    const EPS: f64 = 1e-9;
    let i_start = ((extent.x_min - origin_x) / tw + EPS).floor() as i64;
    let i_end = ((extent.x_max - origin_x) / tw - EPS).ceil() as i64 - 1;
    let j_start = ((extent.y_min - origin_y) / th + EPS).floor() as i64;
    let j_end = ((extent.y_max - origin_y) / th - EPS).ceil() as i64 - 1;
    let n_cols = (i_end - i_start + 1).max(1);
    let n_rows = (j_end - j_start + 1).max(1);

    // Optional intersect_only: pre-convert input features to `geo` geometries.
    let filter = if prm.intersect_only {
        let layer = input.ok_or_else(|| {
            ToolError::Validation("intersect_only requires an input layer".to_string())
        })?;
        Some(collect_geo_features(layer))
    } else {
        None
    };

    let mut layer = new_page_layer(epsg);
    let mut page_number = 0usize;
    // Iterate reading order: rows north→south, columns west→east.
    for r in 0..n_rows {
        // Row r (1-based from the top) corresponds to grid j index j_end - r.
        let j = j_end - r;
        let y0 = origin_y + j as f64 * th;
        let y1 = y0 + th;
        let row_no = r + 1;
        for c in 0..n_cols {
            let i = i_start + c;
            let x0 = origin_x + i as f64 * tw;
            let x1 = x0 + tw;
            let col_no = c + 1;

            if let Some(feats) = filter.as_ref() {
                if !rect_intersects_any(x0, y0, x1, y1, feats) {
                    continue;
                }
            }

            page_number += 1;
            let page_name = match prm.naming {
                Naming::Alphanumeric => format!("{}{}", column_letters(col_no as usize), row_no),
                Naming::Sequential => page_number.to_string(),
            };
            let coords = rect_coords(x0, y0, x1, y1);
            layer.push(Feature {
                fid: 0,
                geometry: Some(Geometry::polygon(coords, vec![])),
                attributes: vec![
                    FieldValue::Text(page_name),
                    FieldValue::Integer(row_no),
                    FieldValue::Integer(col_no),
                    FieldValue::Integer(page_number as i64),
                ],
            });
        }
    }

    ctx.progress
        .info(&format!("grid {n_cols}×{n_rows}, tile {tw:.3}×{th:.3}"));
    Ok((layer, page_number, "grid"))
}

// ── Strip mode ───────────────────────────────────────────────────────────────

/// Builds strip-map pages: rotated rectangles walked along a route line.
fn build_strip(
    prm: &Params,
    input: Option<&Layer>,
    epsg: Option<u32>,
    _ctx: &ToolContext,
) -> Result<(Layer, usize, &'static str), ToolError> {
    let (tw, th) = resolve_tile_size(prm)?;

    // Route: an explicit route layer, else the input layer.
    let route_layer = match prm.route.as_deref() {
        Some(p) => Some(load_input_layer(p)?),
        None => input.cloned(),
    };
    let route_layer = route_layer.ok_or_else(|| {
        ToolError::Validation("strip mode requires a route (or input) line layer".to_string())
    })?;
    let route = first_polyline(&route_layer).ok_or_else(|| {
        ToolError::Validation("strip mode found no line geometry in the route layer".to_string())
    })?;

    let epsg = route_layer.crs_epsg().or(epsg);
    let mut layer = new_page_layer(epsg);

    // Cumulative arc-length parameterization of the route.
    let total = polyline_length(&route);
    if total <= 0.0 {
        return Err(ToolError::Validation(
            "strip mode route has zero length".to_string(),
        ));
    }
    // Advance one page length along the route per step, minus the overlap.
    let step = (th * (1.0 - prm.overlap)).max(th * 1e-3);
    let mut page_number = 0usize;
    let mut s = 0.0;
    loop {
        let center_s = (s + th * 0.5).min(total);
        let (cx, cy) = point_at_distance(&route, center_s);
        let (dx, dy) = tangent_at_distance(&route, center_s);
        // Rotated rectangle: `th` along the route direction, `tw` across it.
        let (px, py) = (-dy, dx); // unit perpendicular
        let hl = th * 0.5;
        let hw = tw * 0.5;
        let corners = [
            (cx - dx * hl - px * hw, cy - dy * hl - py * hw),
            (cx + dx * hl - px * hw, cy + dy * hl - py * hw),
            (cx + dx * hl + px * hw, cy + dy * hl + py * hw),
            (cx - dx * hl + px * hw, cy - dy * hl + py * hw),
        ];
        let coords: Vec<Coord> = corners.iter().map(|&(x, y)| Coord::xy(x, y)).collect();

        page_number += 1;
        let page_name = match prm.naming {
            Naming::Alphanumeric => format!("{}1", column_letters(page_number)),
            Naming::Sequential => page_number.to_string(),
        };
        layer.push(Feature {
            fid: 0,
            geometry: Some(Geometry::polygon(coords, vec![])),
            attributes: vec![
                FieldValue::Text(page_name),
                FieldValue::Integer(1),
                FieldValue::Integer(page_number as i64),
                FieldValue::Integer(page_number as i64),
            ],
        });

        if s + th >= total {
            break;
        }
        s += step;
    }

    Ok((layer, page_number, "strip"))
}

// ── Output layer schema ──────────────────────────────────────────────────────

fn new_page_layer(epsg: Option<u32>) -> Layer {
    let mut l = Layer::new("grid_index").with_geom_type(GeometryType::Polygon);
    if let Some(e) = epsg {
        l = l.with_crs_epsg(e);
    }
    l.add_field(FieldDef::new("page_name", FieldType::Text));
    l.add_field(FieldDef::new("row", FieldType::Integer));
    l.add_field(FieldDef::new("col", FieldType::Integer));
    l.add_field(FieldDef::new("page_number", FieldType::Integer));
    l
}

/// Closed, counter-clockwise corner ring for an axis-aligned rectangle.
fn rect_coords(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Coord> {
    vec![
        Coord::xy(x0, y0),
        Coord::xy(x1, y0),
        Coord::xy(x1, y1),
        Coord::xy(x0, y1),
    ]
}

/// Spreadsheet-style column letters: 1→A, 26→Z, 27→AA, 28→AB, ...
fn column_letters(mut n: usize) -> String {
    let mut out = Vec::new();
    while n > 0 {
        let rem = (n - 1) % 26;
        out.push((b'A' + rem as u8) as char);
        n = (n - 1) / 26;
    }
    out.iter().rev().collect()
}

// ── Extent / tile size resolution ────────────────────────────────────────────

struct Extent {
    x_min: f64,
    y_min: f64,
    x_max: f64,
    y_max: f64,
}

fn resolve_extent(prm: &Params, input: Option<&Layer>) -> Result<Extent, ToolError> {
    if let (Some(x_min), Some(y_min), Some(x_max), Some(y_max)) =
        (prm.x_min, prm.y_min, prm.x_max, prm.y_max)
    {
        if x_max <= x_min || y_max <= y_min {
            return Err(ToolError::Validation(
                "explicit extent must have x_max > x_min and y_max > y_min".to_string(),
            ));
        }
        return Ok(Extent {
            x_min,
            y_min,
            x_max,
            y_max,
        });
    }
    let layer = input.ok_or_else(|| {
        ToolError::Validation(
            "grid mode needs an input layer or all four of x_min/y_min/x_max/y_max".to_string(),
        )
    })?;
    layer_extent(layer)
}

/// Union of the bounding boxes of every feature in the layer.
fn layer_extent(layer: &Layer) -> Result<Extent, ToolError> {
    let (mut x_min, mut y_min) = (f64::INFINITY, f64::INFINITY);
    let (mut x_max, mut y_max) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for f in layer.iter() {
        if let Some(bb) = f.geometry.as_ref().and_then(|g| g.bbox()) {
            x_min = x_min.min(bb.min_x);
            y_min = y_min.min(bb.min_y);
            x_max = x_max.max(bb.max_x);
            y_max = y_max.max(bb.max_y);
        }
    }
    if !(x_min.is_finite() && x_max > x_min && y_max > y_min) {
        return Err(ToolError::Validation(
            "could not derive a valid extent from the input layer".to_string(),
        ));
    }
    Ok(Extent {
        x_min,
        y_min,
        x_max,
        y_max,
    })
}

/// Metric paper dimensions (portrait width, height) for the page presets.
fn page_preset(name: &str) -> Option<(f64, f64)> {
    Some(match name.trim().to_ascii_lowercase().as_str() {
        "a0" => (0.841, 1.189),
        "a1" => (0.594, 0.841),
        "a2" => (0.420, 0.594),
        "a3" => (0.297, 0.420),
        "a4" => (0.210, 0.297),
        "letter" => (0.2159, 0.2794),
        "legal" => (0.2159, 0.3556),
        "tabloid" => (0.2794, 0.4318),
        _ => return None,
    })
}

fn resolve_tile_size(prm: &Params) -> Result<(f64, f64), ToolError> {
    if let (Some(w), Some(h)) = (prm.tile_width, prm.tile_height) {
        return Ok((w, h));
    }
    if let (Some(page), Some(scale)) = (prm.page_size.as_deref(), prm.map_scale) {
        let (pw, ph) = page_preset(page).ok_or_else(|| {
            ToolError::Validation(format!(
                "unknown page_size '{page}' (use a0..a4, letter, legal, tabloid)"
            ))
        })?;
        return Ok((pw * scale, ph * scale));
    }
    Err(ToolError::Validation(
        "tile size needs tile_width+tile_height, or page_size+map_scale".to_string(),
    ))
}

// ── intersect_only: geo geometry conversion + test ───────────────────────────

/// A cached input feature: its bbox (for prefiltering) and a `geo` geometry.
struct GeoFeat {
    bb: [f64; 4], // x_min, y_min, x_max, y_max
    geom: geo::Geometry<f64>,
}

fn collect_geo_features(layer: &Layer) -> Vec<GeoFeat> {
    layer
        .iter()
        .filter_map(|f| {
            let g = f.geometry.as_ref()?;
            let bb = g.bbox()?;
            let geom = to_geo_geometry(g)?;
            Some(GeoFeat {
                bb: [bb.min_x, bb.min_y, bb.max_x, bb.max_y],
                geom,
            })
        })
        .collect()
}

/// True when the axis-aligned rectangle intersects any input feature.
fn rect_intersects_any(x0: f64, y0: f64, x1: f64, y1: f64, feats: &[GeoFeat]) -> bool {
    let rect = geo::Geometry::Polygon(Polygon::new(
        LineString::new(vec![
            GeoCoord { x: x0, y: y0 },
            GeoCoord { x: x1, y: y0 },
            GeoCoord { x: x1, y: y1 },
            GeoCoord { x: x0, y: y1 },
            GeoCoord { x: x0, y: y0 },
        ]),
        vec![],
    ));
    for f in feats {
        // Fast reject on bbox overlap.
        if f.bb[0] > x1 || f.bb[2] < x0 || f.bb[1] > y1 || f.bb[3] < y0 {
            continue;
        }
        if rect.intersects(&f.geom) {
            return true;
        }
    }
    false
}

/// Converts a `wbvector` geometry to a `geo` geometry (polygonal/linear/point).
fn to_geo_geometry(g: &Geometry) -> Option<geo::Geometry<f64>> {
    let gc = |c: &Coord| GeoCoord { x: c.x, y: c.y };
    Some(match g {
        Geometry::Point(c) => geo::Geometry::Point(Point::new(c.x, c.y)),
        Geometry::MultiPoint(cs) => geo::Geometry::MultiPoint(MultiPoint(
            cs.iter().map(|c| Point::new(c.x, c.y)).collect(),
        )),
        Geometry::LineString(cs) => {
            geo::Geometry::LineString(LineString::new(cs.iter().map(gc).collect()))
        }
        Geometry::MultiLineString(parts) => geo::Geometry::MultiLineString(MultiLineString(
            parts
                .iter()
                .map(|p| LineString::new(p.iter().map(gc).collect()))
                .collect(),
        )),
        Geometry::Polygon {
            exterior,
            interiors,
        } => geo::Geometry::Polygon(rings_to_geo(exterior, interiors)),
        Geometry::MultiPolygon(parts) => geo::Geometry::MultiPolygon(geo::MultiPolygon(
            parts.iter().map(|(e, h)| rings_to_geo(e, h)).collect(),
        )),
        Geometry::GeometryCollection(_) => return None,
    })
}

fn rings_to_geo(exterior: &wbvector::Ring, interiors: &[wbvector::Ring]) -> Polygon<f64> {
    let ring_ls = |r: &wbvector::Ring| {
        LineString::new(
            r.coords()
                .iter()
                .map(|c| GeoCoord { x: c.x, y: c.y })
                .collect(),
        )
    };
    Polygon::new(ring_ls(exterior), interiors.iter().map(ring_ls).collect())
}

// ── Route geometry (strip mode) ──────────────────────────────────────────────

/// First polyline in the layer as a flat vertex list.
fn first_polyline(layer: &Layer) -> Option<Vec<(f64, f64)>> {
    for f in layer.iter() {
        match f.geometry.as_ref() {
            Some(Geometry::LineString(cs)) if cs.len() >= 2 => {
                return Some(cs.iter().map(|c| (c.x, c.y)).collect());
            }
            Some(Geometry::MultiLineString(parts)) => {
                if let Some(p) = parts.iter().find(|p| p.len() >= 2) {
                    return Some(p.iter().map(|c| (c.x, c.y)).collect());
                }
            }
            _ => {}
        }
    }
    None
}

fn polyline_length(pts: &[(f64, f64)]) -> f64 {
    pts.windows(2)
        .map(|w| (w[1].0 - w[0].0).hypot(w[1].1 - w[0].1))
        .sum()
}

/// Point at arc-length `s` along the polyline (clamped to its ends).
fn point_at_distance(pts: &[(f64, f64)], s: f64) -> (f64, f64) {
    if s <= 0.0 {
        return pts[0];
    }
    let mut acc = 0.0;
    for w in pts.windows(2) {
        let seg = (w[1].0 - w[0].0).hypot(w[1].1 - w[0].1);
        if acc + seg >= s {
            let t = if seg > 0.0 { (s - acc) / seg } else { 0.0 };
            return (
                w[0].0 + (w[1].0 - w[0].0) * t,
                w[0].1 + (w[1].1 - w[0].1) * t,
            );
        }
        acc += seg;
    }
    *pts.last().unwrap()
}

/// Unit tangent (route direction) at arc-length `s`.
fn tangent_at_distance(pts: &[(f64, f64)], s: f64) -> (f64, f64) {
    let mut acc = 0.0;
    for w in pts.windows(2) {
        let seg = (w[1].0 - w[0].0).hypot(w[1].1 - w[0].1);
        if acc + seg >= s || acc + seg >= polyline_length(pts) {
            let (dx, dy) = (w[1].0 - w[0].0, w[1].1 - w[0].1);
            let len = dx.hypot(dy);
            if len > 0.0 {
                return (dx / len, dy / len);
            }
        }
        acc += seg;
    }
    (1.0, 0.0)
}

// ── Parameters ───────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq)]
enum Mode {
    Grid,
    Strip,
}

#[derive(Clone, Copy, PartialEq)]
enum Naming {
    Alphanumeric,
    Sequential,
}

struct Params {
    input: Option<String>,
    mode: Mode,
    x_min: Option<f64>,
    y_min: Option<f64>,
    x_max: Option<f64>,
    y_max: Option<f64>,
    tile_width: Option<f64>,
    tile_height: Option<f64>,
    page_size: Option<String>,
    map_scale: Option<f64>,
    origin_x: Option<f64>,
    origin_y: Option<f64>,
    naming: Naming,
    intersect_only: bool,
    route: Option<String>,
    overlap: f64,
    epsg: Option<u32>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let input = parse_optional_str(args, "input")?.map(str::to_string);
    let route = parse_optional_str(args, "route")?.map(str::to_string);
    let page_size = parse_optional_str(args, "page_size")?.map(str::to_string);

    let mode = match parse_optional_str(args, "mode")?.map(|s| s.trim().to_ascii_lowercase()) {
        None => Mode::Grid,
        Some(s) if s == "grid" => Mode::Grid,
        Some(s) if s == "strip" => Mode::Strip,
        Some(s) => {
            return Err(ToolError::Validation(format!(
                "mode must be 'grid' or 'strip', got '{s}'"
            )))
        }
    };

    let naming = match parse_optional_str(args, "naming")?.map(|s| s.trim().to_ascii_lowercase()) {
        None => Naming::Alphanumeric,
        Some(s) if s == "alphanumeric" => Naming::Alphanumeric,
        Some(s) if s == "sequential" => Naming::Sequential,
        Some(s) => {
            return Err(ToolError::Validation(format!(
                "naming must be 'alphanumeric' or 'sequential', got '{s}'"
            )))
        }
    };

    let tile_width = opt_pos(args, "tile_width")?;
    let tile_height = opt_pos(args, "tile_height")?;
    let map_scale = opt_pos(args, "map_scale")?;
    let overlap = match opt_f64(args, "overlap")? {
        None => 0.05,
        Some(v) if (0.0..1.0).contains(&v) => v,
        Some(_) => {
            return Err(ToolError::Validation(
                "overlap must be in [0, 1)".to_string(),
            ))
        }
    };

    // Validate the page_size preset name eagerly when provided.
    if let Some(page) = page_size.as_deref() {
        if page_preset(page).is_none() {
            return Err(ToolError::Validation(format!(
                "unknown page_size '{page}' (use a0..a4, letter, legal, tabloid)"
            )));
        }
    }

    let intersect_only = opt_bool(args, "intersect_only")?.unwrap_or(false);
    let epsg = opt_epsg(args, "epsg")?;

    Ok(Params {
        input,
        mode,
        x_min: opt_f64(args, "x_min")?,
        y_min: opt_f64(args, "y_min")?,
        x_max: opt_f64(args, "x_max")?,
        y_max: opt_f64(args, "y_max")?,
        tile_width,
        tile_height,
        page_size,
        map_scale,
        origin_x: opt_f64(args, "origin_x")?,
        origin_y: opt_f64(args, "origin_y")?,
        naming,
        intersect_only,
        route,
        overlap,
        epsg,
    })
}

fn opt_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
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

fn opt_pos(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match opt_f64(args, key)? {
        Some(v) if v > 0.0 && v.is_finite() => Ok(Some(v)),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a positive number"
        ))),
        None => Ok(None),
    }
}

fn opt_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
}

fn opt_epsg(args: &ToolArgs, key: &str) -> Result<Option<u32>, ToolError> {
    match opt_f64(args, key)? {
        Some(v) if v > 0.0 && v.fract() == 0.0 => Ok(Some(v as u32)),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a positive integer EPSG code"
        ))),
        None => Ok(None),
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

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = GridIndexFeaturesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn poly_area(g: &Geometry) -> f64 {
        // Shoelace area of a simple polygon exterior ring.
        if let Geometry::Polygon { exterior, .. } = g {
            let c = exterior.coords();
            let n = c.len();
            let mut a = 0.0;
            for i in 0..n {
                let j = (i + 1) % n;
                a += c[i].x * c[j].y - c[j].x * c[i].y;
            }
            (a * 0.5).abs()
        } else {
            0.0
        }
    }

    /// A full grid tiles the extent exactly: no gaps/overlaps, and the summed
    /// page area equals the extent area when the tile size divides it.
    #[test]
    fn full_grid_tiles_extent_exactly() {
        let (out, layer) = run(json!({
            "x_min": 0.0, "y_min": 0.0, "x_max": 100.0, "y_max": 60.0,
            "tile_width": 25.0, "tile_height": 20.0, "epsg": 3857,
        }));
        // 100/25 = 4 cols, 60/20 = 3 rows -> 12 pages.
        assert_eq!(out.outputs["page_count"], json!(12));
        assert_eq!(layer.features.len(), 12);
        let total: f64 = layer
            .features
            .iter()
            .filter_map(|f| f.geometry.as_ref())
            .map(poly_area)
            .sum();
        assert!(
            (total - 100.0 * 60.0).abs() < 1e-6,
            "summed page area {total} should equal the extent area 6000"
        );
    }

    /// Alphanumeric naming: top-left page is A1, and row numbering runs top→down.
    #[test]
    fn alphanumeric_naming_top_left_is_a1() {
        let (_out, layer) = run(json!({
            "x_min": 0.0, "y_min": 0.0, "x_max": 100.0, "y_max": 60.0,
            "tile_width": 25.0, "tile_height": 20.0,
        }));
        let name_i = layer.schema.field_index("page_name").unwrap();
        let row_i = layer.schema.field_index("row").unwrap();
        let col_i = layer.schema.field_index("col").unwrap();
        // First feature is reading order (row 1, col 1) = top-left.
        let f0 = &layer.features[0];
        assert_eq!(f0.attributes[name_i].as_str(), Some("A1"));
        assert_eq!(f0.attributes[row_i].as_i64(), Some(1));
        assert_eq!(f0.attributes[col_i].as_i64(), Some(1));
        // The top row (row 1) must sit at the top of the extent (max Y).
        let ymax = |f: &Feature| f.geometry.as_ref().unwrap().bbox().unwrap().max_y;
        assert!((ymax(f0) - 60.0).abs() < 1e-9);
        // Last column of the top row is "D1" (4 columns).
        let d1 = layer
            .features
            .iter()
            .find(|f| f.attributes[name_i].as_str() == Some("D1"))
            .unwrap();
        assert_eq!(d1.attributes[col_i].as_i64(), Some(4));
    }

    /// intersect_only drops pages that miss the input; every kept page hits it.
    #[test]
    fn intersect_only_keeps_only_hit_pages() {
        // A single point near the lower-left of a 4×3 grid.
        let mut pts = Layer::new("p")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        pts.add_field(FieldDef::new("id", FieldType::Integer));
        pts.add_feature(Some(Geometry::point(10.0, 10.0)), &[("id", 1i64.into())])
            .unwrap();
        let id = memory_store::put_vector(pts);
        let path = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run(json!({
            "input": path, "intersect_only": true,
            "x_min": 0.0, "y_min": 0.0, "x_max": 100.0, "y_max": 60.0,
            "tile_width": 25.0, "tile_height": 20.0,
        }));
        // Only the page containing (10,10) survives.
        assert_eq!(out.outputs["page_count"], json!(1));
        assert_eq!(layer.features.len(), 1);
        let bb = layer.features[0].geometry.as_ref().unwrap().bbox().unwrap();
        assert!(bb.min_x <= 10.0 && 10.0 <= bb.max_x && bb.min_y <= 10.0 && 10.0 <= bb.max_y);
    }

    /// Sequential naming numbers pages 1..N in reading order.
    #[test]
    fn sequential_naming_numbers_pages() {
        let (_out, layer) = run(json!({
            "x_min": 0.0, "y_min": 0.0, "x_max": 50.0, "y_max": 40.0,
            "tile_width": 25.0, "tile_height": 20.0, "naming": "sequential",
        }));
        let name_i = layer.schema.field_index("page_name").unwrap();
        let pn_i = layer.schema.field_index("page_number").unwrap();
        for (k, f) in layer.features.iter().enumerate() {
            assert_eq!(f.attributes[pn_i].as_i64(), Some(k as i64 + 1));
            assert_eq!(
                f.attributes[name_i].as_str(),
                Some((k + 1).to_string().as_str())
            );
        }
    }

    /// Strip mode walks a straight route and covers its full length.
    #[test]
    fn strip_mode_covers_route() {
        let mut line = Layer::new("r")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        line.add_field(FieldDef::new("id", FieldType::Integer));
        line.add_feature(
            Some(Geometry::line_string(vec![
                Coord::xy(0.0, 0.0),
                Coord::xy(1000.0, 0.0),
            ])),
            &[("id", 1i64.into())],
        )
        .unwrap();
        let id = memory_store::put_vector(line);
        let path = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run(json!({
            "mode": "strip", "route": path,
            "tile_width": 200.0, "tile_height": 300.0, "overlap": 0.0,
        }));
        // 1000 / 300 -> 4 pages (0-300, 300-600, 600-900, 900-1000).
        assert_eq!(out.outputs["page_count"], json!(4));
        assert!(layer.features.len() >= 4);
        // Axis-aligned route -> each page is 200 wide (Y) × 300 long (X).
        let bb = layer.features[0].geometry.as_ref().unwrap().bbox().unwrap();
        assert!((bb.max_x - bb.min_x - 300.0).abs() < 1e-6);
        assert!((bb.max_y - bb.min_y - 200.0).abs() < 1e-6);
    }

    /// Column-letter helper rolls over past Z like spreadsheet columns.
    #[test]
    fn column_letters_roll_over() {
        assert_eq!(column_letters(1), "A");
        assert_eq!(column_letters(26), "Z");
        assert_eq!(column_letters(27), "AA");
        assert_eq!(column_letters(28), "AB");
    }

    /// page_size + map_scale derives a metric tile size.
    #[test]
    fn page_size_scale_derives_tile() {
        // A4 portrait = 0.210 × 0.297 m; at 1:10000 -> 2100 × 2970 m.
        let (_out, layer) = run(json!({
            "x_min": 0.0, "y_min": 0.0, "x_max": 4200.0, "y_max": 2970.0,
            "page_size": "a4", "map_scale": 10000.0,
        }));
        // 4200/2100 = 2 cols, 2970/2970 = 1 row -> 2 pages.
        assert_eq!(layer.features.len(), 2);
        let bb = layer.features[0].geometry.as_ref().unwrap().bbox().unwrap();
        assert!((bb.max_x - bb.min_x - 2100.0).abs() < 1e-3);
        assert!((bb.max_y - bb.min_y - 2970.0).abs() < 1e-3);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            GridIndexFeaturesTool.validate(&args)
        };
        // Unknown mode.
        assert!(bad(json!({ "mode": "spiral" })).is_err());
        // Unknown page preset.
        assert!(bad(json!({ "page_size": "a9", "map_scale": 1000.0 })).is_err());
        // Negative tile size.
        assert!(bad(json!({ "tile_width": -5.0 })).is_err());
        // Overlap out of range.
        assert!(bad(json!({ "overlap": 1.5 })).is_err());
        // A valid explicit-extent + tile grid validates.
        assert!(bad(json!({
            "x_min": 0.0, "y_min": 0.0, "x_max": 10.0, "y_max": 10.0,
            "tile_width": 5.0, "tile_height": 5.0
        }))
        .is_ok());
    }

    /// Missing both extent and input fails at run time (grid mode).
    #[test]
    fn grid_without_extent_or_input_errors() {
        let args: ToolArgs = serde_json::from_value(json!({
            "tile_width": 5.0, "tile_height": 5.0
        }))
        .unwrap();
        assert!(GridIndexFeaturesTool.run(&args, &ctx()).is_err());
    }
}
