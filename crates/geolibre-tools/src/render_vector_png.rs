//! GeoLibre tool: draw a vector layer to a PNG map image.
//!
//! The vector analog of `render_raster_png`: rasterize points, lines, and
//! polygons (with holes) to an RGBA PNG so a layer can be previewed inline
//! (e.g. in a notebook) without a mapping library. Single fill/stroke colors;
//! the image aspect ratio follows the data unless a height is given.

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata,
    ToolParamSpec, ToolRunResult,
};
use wbvector::{Coord, Geometry, Ring};

use crate::common::write_bytes;
use crate::render::encode_png_rgba;
use crate::vector_common::load_input_layer;

type Rgba = [u8; 4];

/// Renders a vector layer to an RGBA PNG.
pub struct RenderVectorPngTool;

impl Tool for RenderVectorPngTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "render_vector_png",
            display_name: "Render Vector to PNG",
            summary: "Draw a vector layer (points/lines/polygons) to a PNG map image.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec { name: "input", description: "Input vector file path or in-memory handle.", required: true },
                ToolParamSpec { name: "output", description: "Output PNG file path (e.g. /work/map.png).", required: true },
                ToolParamSpec { name: "width", description: "Image width in pixels (default 800).", required: false },
                ToolParamSpec { name: "height", description: "Image height in pixels (default: from the data aspect ratio).", required: false },
                ToolParamSpec { name: "fill", description: "Polygon/point fill color as #rrggbb or #rrggbbaa (default semi-transparent blue).", required: false },
                ToolParamSpec { name: "stroke", description: "Outline/line color as #rrggbb or #rrggbbaa (default dark blue).", required: false },
                ToolParamSpec { name: "stroke_width", description: "Outline/line width in pixels (default 1).", required: false },
                ToolParamSpec { name: "background", description: "Background color as #rrggbb[aa] or 'transparent' (default transparent).", required: false },
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
        for key in ["fill", "stroke", "background"] {
            if let Some(c) = args.get(key).and_then(Value::as_str) {
                parse_color(c)?;
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = require_str(args, "output")?;
        let fill = color_or(args, "fill", [51, 136, 255, 150])?;
        let stroke = color_or(args, "stroke", [20, 60, 120, 255])?;
        let background = color_or(args, "background", [0, 0, 0, 0])?;
        let stroke_width = args.get("stroke_width").and_then(Value::as_u64).unwrap_or(1).max(1) as i64;
        let width = args.get("width").and_then(Value::as_u64).unwrap_or(800).clamp(16, 8192) as i64;
        let pad: f64 = 8.0;

        let layer = load_input_layer(input)?;
        let extent = layer_extent(&layer)
            .ok_or_else(|| ToolError::Execution("layer has no drawable geometry".to_string()))?;
        let (min_x, min_y, max_x, max_y) = extent;

        // Height follows the data aspect ratio unless overridden.
        let span_x = (max_x - min_x).max(f64::MIN_POSITIVE);
        let span_y = (max_y - min_y).max(f64::MIN_POSITIVE);
        let height = match args.get("height").and_then(Value::as_u64) {
            Some(h) => (h as i64).clamp(16, 8192),
            None => {
                let inner_w = (width as f64 - 2.0 * pad).max(1.0);
                let h = inner_w * span_y / span_x + 2.0 * pad;
                (h.round() as i64).clamp(16, 8192)
            }
        };

        let inner_w = width as f64 - 2.0 * pad;
        let inner_h = height as f64 - 2.0 * pad;
        // World -> pixel. Y is flipped (north up). Degenerate spans center the data.
        let to_px = |c: &Coord| -> (f64, f64) {
            let fx = if max_x > min_x { (c.x - min_x) / span_x } else { 0.5 };
            let fy = if max_y > min_y { (max_y - c.y) / span_y } else { 0.5 };
            (pad + fx * inner_w, pad + fy * inner_h)
        };

        ctx.progress.info("rasterizing features");
        let mut canvas = Canvas::new(width, height, background);
        for feature in &layer.features {
            if let Some(geom) = &feature.geometry {
                draw_geometry(&mut canvas, geom, &to_px, fill, stroke, stroke_width);
            }
        }

        let png = encode_png_rgba(&canvas.rgba, width as u32, height as u32)?;
        write_bytes(output, &png)?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(output));
        outputs.insert("width".to_string(), json!(width));
        outputs.insert("height".to_string(), json!(height));
        outputs.insert("feature_count".to_string(), json!(layer.features.len()));
        Ok(ToolRunResult { outputs })
    }
}

/// A simple RGBA8 raster canvas with straight-alpha "source over" compositing.
struct Canvas {
    w: i64,
    h: i64,
    rgba: Vec<u8>,
}

impl Canvas {
    fn new(w: i64, h: i64, background: Rgba) -> Self {
        let mut rgba = vec![0u8; (w * h * 4) as usize];
        if background[3] != 0 {
            for px in rgba.chunks_exact_mut(4) {
                px.copy_from_slice(&background);
            }
        }
        Self { w, h, rgba }
    }

    fn blend(&mut self, x: i64, y: i64, c: Rgba) {
        if x < 0 || y < 0 || x >= self.w || y >= self.h || c[3] == 0 {
            return;
        }
        let i = ((y * self.w + x) * 4) as usize;
        let sa = c[3] as f64 / 255.0;
        let da = self.rgba[i + 3] as f64 / 255.0;
        let out_a = sa + da * (1.0 - sa);
        for ch in 0..3 {
            let s = c[ch] as f64;
            let d = self.rgba[i + ch] as f64;
            let v = if out_a > 0.0 {
                (s * sa + d * da * (1.0 - sa)) / out_a
            } else {
                0.0
            };
            self.rgba[i + ch] = v.round().clamp(0.0, 255.0) as u8;
        }
        self.rgba[i + 3] = (out_a * 255.0).round().clamp(0.0, 255.0) as u8;
    }

    /// Fills a square of side `2*r+1` centered at `(x, y)` (used to thicken
    /// lines and draw points).
    fn dot(&mut self, x: i64, y: i64, r: i64, c: Rgba) {
        for dy in -r..=r {
            for dx in -r..=r {
                self.blend(x + dx, y + dy, c);
            }
        }
    }
}

/// Recursively draws a geometry: polygons are filled (even-odd, holes honored)
/// then outlined; lines are stroked; points are dots.
fn draw_geometry<F: Fn(&Coord) -> (f64, f64)>(
    canvas: &mut Canvas,
    geom: &Geometry,
    to_px: &F,
    fill: Rgba,
    stroke: Rgba,
    sw: i64,
) {
    match geom {
        Geometry::Point(c) => {
            let (x, y) = to_px(c);
            canvas.dot(x.round() as i64, y.round() as i64, sw.max(2), fill);
        }
        Geometry::MultiPoint(cs) => {
            for c in cs {
                let (x, y) = to_px(c);
                canvas.dot(x.round() as i64, y.round() as i64, sw.max(2), fill);
            }
        }
        Geometry::LineString(cs) => draw_line_string(canvas, cs, to_px, stroke, sw),
        Geometry::MultiLineString(lines) => {
            for cs in lines {
                draw_line_string(canvas, cs, to_px, stroke, sw);
            }
        }
        Geometry::Polygon { exterior, interiors } => {
            draw_polygon(canvas, exterior, interiors, to_px, fill, stroke, sw)
        }
        Geometry::MultiPolygon(polys) => {
            for (exterior, interiors) in polys {
                draw_polygon(canvas, exterior, interiors, to_px, fill, stroke, sw);
            }
        }
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                draw_geometry(canvas, g, to_px, fill, stroke, sw);
            }
        }
    }
}

fn draw_line_string<F: Fn(&Coord) -> (f64, f64)>(
    canvas: &mut Canvas,
    coords: &[Coord],
    to_px: &F,
    stroke: Rgba,
    sw: i64,
) {
    let pts: Vec<(f64, f64)> = coords.iter().map(to_px).collect();
    for seg in pts.windows(2) {
        draw_segment(canvas, seg[0], seg[1], stroke, sw);
    }
}

fn draw_polygon<F: Fn(&Coord) -> (f64, f64)>(
    canvas: &mut Canvas,
    exterior: &Ring,
    interiors: &[Ring],
    to_px: &F,
    fill: Rgba,
    stroke: Rgba,
    sw: i64,
) {
    // Even-odd scanline fill over the exterior plus all interior rings; an
    // even number of crossings inside a hole flips back to "outside", so holes
    // are honored automatically.
    let rings: Vec<Vec<(f64, f64)>> = std::iter::once(exterior)
        .chain(interiors)
        .map(|r| r.coords().iter().map(to_px).collect())
        .collect();
    if fill[3] != 0 {
        fill_rings(canvas, &rings, fill);
    }
    for ring in &rings {
        for seg in ring.windows(2) {
            draw_segment(canvas, seg[0], seg[1], stroke, sw);
        }
        // Close the ring if the data did not repeat the first vertex.
        if let (Some(&first), Some(&last)) = (ring.first(), ring.last()) {
            if first != last {
                draw_segment(canvas, last, first, stroke, sw);
            }
        }
    }
}

/// Even-odd polygon fill across all rings.
fn fill_rings(canvas: &mut Canvas, rings: &[Vec<(f64, f64)>], fill: Rgba) {
    let mut y_min = f64::INFINITY;
    let mut y_max = f64::NEG_INFINITY;
    for ring in rings {
        for &(_, y) in ring {
            y_min = y_min.min(y);
            y_max = y_max.max(y);
        }
    }
    if !y_min.is_finite() {
        return;
    }
    let y0 = (y_min.floor() as i64).max(0);
    let y1 = (y_max.ceil() as i64).min(canvas.h - 1);
    for y in y0..=y1 {
        let yc = y as f64 + 0.5;
        let mut xs: Vec<f64> = Vec::new();
        for ring in rings {
            let len = ring.len();
            if len < 2 {
                continue;
            }
            // Iterate every edge including the closing one (`(i + 1) % len`) so a
            // ring that does not repeat its first vertex still fills correctly;
            // a ring that does repeat it just yields a degenerate (no-crossing)
            // closing edge.
            for i in 0..len {
                let (x0, ya) = ring[i];
                let (x1, yb) = ring[(i + 1) % len];
                // Half-open edge test avoids double-counting shared vertices.
                if (ya <= yc && yb > yc) || (yb <= yc && ya > yc) {
                    let t = (yc - ya) / (yb - ya);
                    xs.push(x0 + t * (x1 - x0));
                }
            }
        }
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mut i = 0;
        while i + 1 < xs.len() {
            let xa = xs[i].round() as i64;
            let xb = xs[i + 1].round() as i64;
            for x in xa..xb {
                canvas.blend(x, y, fill);
            }
            i += 2;
        }
    }
}

/// Bresenham line between two pixel points, thickened to `sw` pixels.
fn draw_segment(canvas: &mut Canvas, a: (f64, f64), b: (f64, f64), c: Rgba, sw: i64) {
    let r = (sw - 1).max(0) / 2;
    let (mut x0, mut y0) = (a.0.round() as i64, a.1.round() as i64);
    let (x1, y1) = (b.0.round() as i64, b.1.round() as i64);
    let dx = (x1 - x0).abs();
    let dy = -(y1 - y0).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let sy = if y0 < y1 { 1 } else { -1 };
    let mut err = dx + dy;
    loop {
        if r == 0 {
            canvas.blend(x0, y0, c);
        } else {
            canvas.dot(x0, y0, r, c);
        }
        if x0 == x1 && y0 == y1 {
            break;
        }
        let e2 = 2 * err;
        if e2 >= dy {
            err += dy;
            x0 += sx;
        }
        if e2 <= dx {
            err += dx;
            y0 += sy;
        }
    }
}

/// Bounding box over every feature geometry: `(min_x, min_y, max_x, max_y)`.
fn layer_extent(layer: &wbvector::Layer) -> Option<(f64, f64, f64, f64)> {
    let (mut min_x, mut min_y, mut max_x, mut max_y) =
        (f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut any = false;
    for f in &layer.features {
        if let Some(b) = f.geometry.as_ref().and_then(|g| g.bbox()) {
            min_x = min_x.min(b.min_x);
            min_y = min_y.min(b.min_y);
            max_x = max_x.max(b.max_x);
            max_y = max_y.max(b.max_y);
            any = true;
        }
    }
    any.then_some((min_x, min_y, max_x, max_y))
}

fn color_or(args: &ToolArgs, key: &str, default: Rgba) -> Result<Rgba, ToolError> {
    match args.get(key).and_then(Value::as_str) {
        Some(c) => parse_color(c),
        None => Ok(default),
    }
}

/// Parses `#rrggbb`, `#rrggbbaa`, or `transparent`/`none` into RGBA.
fn parse_color(s: &str) -> Result<Rgba, ToolError> {
    let s = s.trim();
    if s.eq_ignore_ascii_case("transparent") || s.eq_ignore_ascii_case("none") {
        return Ok([0, 0, 0, 0]);
    }
    let hex = s.strip_prefix('#').unwrap_or(s);
    // Guard before byte-slicing: a multi-byte char could otherwise make a 6/8
    // byte string slice across a char boundary and panic.
    if !hex.is_ascii() {
        return Err(color_err(s));
    }
    let byte = |i: usize| u8::from_str_radix(&hex[i..i + 2], 16).ok();
    match hex.len() {
        6 => match (byte(0), byte(2), byte(4)) {
            (Some(r), Some(g), Some(b)) => Ok([r, g, b, 255]),
            _ => Err(color_err(s)),
        },
        8 => match (byte(0), byte(2), byte(4), byte(6)) {
            (Some(r), Some(g), Some(b), Some(a)) => Ok([r, g, b, a]),
            _ => Err(color_err(s)),
        },
        _ => Err(color_err(s)),
    }
}

fn color_err(s: &str) -> ToolError {
    ToolError::Validation(format!(
        "invalid color '{s}' (expected #rrggbb, #rrggbbaa, or transparent)"
    ))
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}
