//! GeoLibre tool: strip-map atlas pages that follow a route.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Strip Map Index Features*
//! (Cartography). `grid_index_features` tiles a rectangular extent into a
//! regular fishnet; this tool instead walks each input polyline, accumulating
//! arc length, and emits a rotated rectangular page for every
//! `page_length * (1 - overlap)` units travelled. Each page's long axis is
//! aligned to the local chord bearing of the route (or snapped to the nearest
//! horizontal/vertical axis), so a highway, river, pipeline, or trail is
//! covered end-to-end by a numbered, overlapping sequence of map-book pages —
//! the standard "strip map" atlas layout for long, narrow corridors that a
//! regular grid index would slice at arbitrary angles.
//!
//! Algorithm: for each input line (each part of a `MultiLineString` walked
//! independently), an entry station starts at 0 and advances by
//! `page_length * (1 - overlap)` until it passes the line's end. Each page
//! spans stations `[entry, entry + page_length]` (the last page's exit is
//! clamped to the line's total length); the chord between the entry and exit
//! points on the route gives the page's centre and bearing. Rectangle corners
//! are then built `page_length` long (along the chord) by `page_width` wide
//! (across it) and rotated to that bearing — or, for `horizontal`/`vertical`
//! orientation, snapped to whichever axis (east-west or north-south) is
//! nearest the chord bearing. Pages are numbered sequentially from
//! `start_page` in route-following (left-to-right) reading order, carrying
//! the source feature's index as `line_id` and the applied rotation in
//! degrees.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct StripMapIndexFeaturesTool;

impl Tool for StripMapIndexFeaturesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "strip_map_index_features",
            display_name: "Strip Map Index Features",
            summary: "Walk each input line and emit a numbered sequence of rotated, overlapping page rectangles oriented to the local route bearing, for strip-map atlases along a highway, river, pipeline, or trail — like ArcGIS Strip Map Index Features.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input line vector layer to index.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon vector path of index pages (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "page_length",
                    description: "Along-route length of each page, in CRS units. Required.",
                    required: true,
                },
                ToolParamSpec {
                    name: "page_width",
                    description: "Cross-route width of each page, in CRS units. Required.",
                    required: true,
                },
                ToolParamSpec {
                    name: "overlap",
                    description: "Fraction [0, 1) of page_length that consecutive pages overlap along the route. Default 0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "orientation",
                    description: "along_line (default): rotate each page to its local chord bearing. horizontal/vertical: snap rotation to the nearest east-west/north-south axis.",
                    required: false,
                },
                ToolParamSpec {
                    name: "start_page",
                    description: "First page number to assign. Default 1.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args
            .get("input")
            .and_then(Value::as_str)
            .map(str::trim)
            .unwrap_or("")
            .is_empty()
        {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args
            .get("input")
            .and_then(Value::as_str)
            .filter(|s| !s.trim().is_empty())
            .ok_or_else(|| {
                ToolError::Validation("missing required parameter 'input'".to_string())
            })?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;

        let mut out = Layer::new("strip_map_index").with_geom_type(GeometryType::Polygon);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("page_number", FieldType::Integer));
        out.add_field(FieldDef::new("line_id", FieldType::Integer));
        out.add_field(FieldDef::new("rotation", FieldType::Float));

        let mut page_number = prm.start_page;
        let mut max_deviation_deg = 0.0f64;
        for (fidx, feature) in layer.features.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            for chain in line_chains(geom) {
                page_number += emit_pages(
                    &mut out,
                    &chain,
                    fidx as i64,
                    page_number,
                    &prm,
                    &mut max_deviation_deg,
                )?;
            }
        }
        let page_count = (page_number - prm.start_page) as usize;

        ctx.progress.info(&format!(
            "generated {page_count} strip-map page(s) starting at {}",
            prm.start_page
        ));

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("page_count".to_string(), json!(page_count));
        Ok(ToolRunResult { outputs })
    }
}

/// Emits pages along one polyline, returning the number produced. Updates
/// `max_deviation_deg` with the largest snap deviation seen (0 for
/// `along_line` orientation, which never deviates from the chord bearing).
fn emit_pages(
    out: &mut Layer,
    pts: &[P],
    line_id: i64,
    first_page_number: i64,
    prm: &Params,
    max_deviation_deg: &mut f64,
) -> Result<i64, ToolError> {
    if pts.len() < 2 {
        return Ok(0);
    }
    let mut cum = vec![0.0f64; pts.len()];
    for i in 1..pts.len() {
        cum[i] = cum[i - 1] + dist(pts[i - 1], pts[i]);
    }
    let total = *cum.last().unwrap();
    if total <= 0.0 {
        return Ok(0);
    }

    let step = prm.page_length * (1.0 - prm.overlap);
    let mut entry = 0.0f64;
    let mut count = 0i64;
    while entry < total - 1e-9 {
        let exit = (entry + prm.page_length).min(total);
        let entry_pt = point_at(pts, &cum, entry);
        let exit_pt = point_at(pts, &cum, exit);

        // Chord bearing (compass: 0 = north/+y, clockwise); degenerate chords
        // (zero-length last page) fall back to the local tangent.
        let (chord_dx, chord_dy) = (exit_pt.x - entry_pt.x, exit_pt.y - entry_pt.y);
        let chord_len = chord_dx.hypot(chord_dy);
        let (dir_x, dir_y) = if chord_len > 1e-12 {
            (chord_dx / chord_len, chord_dy / chord_len)
        } else {
            point_and_tangent(pts, &cum, entry)
        };
        let bearing = bearing_deg(dir_x, dir_y);

        let rotation = match prm.orientation {
            Orientation::AlongLine => bearing,
            Orientation::Horizontal => nearest_axis(bearing, &[90.0, 270.0]),
            Orientation::Vertical => nearest_axis(bearing, &[0.0, 180.0]),
        };
        let deviation = angle_diff(bearing, rotation);
        if deviation > *max_deviation_deg {
            *max_deviation_deg = deviation;
        }

        let rot_rad = rotation.to_radians();
        // Direction unit vector for the (possibly snapped) rotation, using the
        // same compass convention: 0deg = +y, 90deg = +x.
        let (rdx, rdy) = (rot_rad.sin(), rot_rad.cos());
        let (px, py) = (rdy, -rdx); // perpendicular (rotate -90deg)

        let cx = (entry_pt.x + exit_pt.x) * 0.5;
        let cy = (entry_pt.y + exit_pt.y) * 0.5;
        let hl = prm.page_length * 0.5;
        let hw = prm.page_width * 0.5;
        // Counter-clockwise ring: back-left, back-right, front-right, front-left.
        let corners = [
            (cx - rdx * hl - px * hw, cy - rdy * hl - py * hw),
            (cx - rdx * hl + px * hw, cy - rdy * hl + py * hw),
            (cx + rdx * hl + px * hw, cy + rdy * hl + py * hw),
            (cx + rdx * hl - px * hw, cy + rdy * hl - py * hw),
        ];
        let coords: Vec<Coord> = corners.iter().map(|&(x, y)| Coord::xy(x, y)).collect();

        out.push(Feature {
            fid: 0,
            geometry: Some(Geometry::polygon(coords, vec![])),
            attributes: vec![
                FieldValue::Integer(first_page_number + count),
                FieldValue::Integer(line_id),
                FieldValue::Float(rotation),
            ],
        });
        count += 1;
        entry += step;
    }
    Ok(count)
}

/// Point at arc-length `s` along the polyline (clamped to its ends).
fn point_at(pts: &[P], cum: &[f64], s: f64) -> P {
    point_and_tangent_impl(pts, cum, s).0
}

/// Local unit tangent at arc-length `s` (used only as a degenerate-chord
/// fallback, so it need not carry the point).
fn point_and_tangent(pts: &[P], cum: &[f64], s: f64) -> (f64, f64) {
    point_and_tangent_impl(pts, cum, s).1
}

fn point_and_tangent_impl(pts: &[P], cum: &[f64], s: f64) -> (P, (f64, f64)) {
    let total = *cum.last().unwrap_or(&0.0);
    let s = s.clamp(0.0, total);
    let mut i = 0;
    while i + 1 < pts.len() && cum[i + 1] < s {
        i += 1;
    }
    let (a, b) = (pts[i], pts[(i + 1).min(pts.len() - 1)]);
    let seg_len = dist(a, b);
    if seg_len <= 0.0 {
        return (a, (1.0, 0.0));
    }
    let t = ((s - cum[i]) / seg_len).clamp(0.0, 1.0);
    let p = P {
        x: a.x + (b.x - a.x) * t,
        y: a.y + (b.y - a.y) * t,
    };
    let tangent = ((b.x - a.x) / seg_len, (b.y - a.y) / seg_len);
    (p, tangent)
}

/// Compass bearing (0-360, 0 = north/+y, clockwise) of a direction vector.
fn bearing_deg(dx: f64, dy: f64) -> f64 {
    let deg = dx.atan2(dy).to_degrees();
    (deg + 360.0) % 360.0
}

/// Minimal absolute difference between two compass bearings, in [0, 180].
fn angle_diff(a: f64, b: f64) -> f64 {
    let d = (a - b).rem_euclid(360.0);
    d.min(360.0 - d)
}

/// Chooses whichever candidate axis bearing is closest (circularly) to `bearing`.
fn nearest_axis(bearing: f64, axes: &[f64]) -> f64 {
    axes.iter()
        .copied()
        .min_by(|&a, &b| {
            angle_diff(bearing, a)
                .partial_cmp(&angle_diff(bearing, b))
                .unwrap()
        })
        .unwrap_or(bearing)
}

// ── Geometry helpers ─────────────────────────────────────────────────────────

#[derive(Clone, Copy)]
struct P {
    x: f64,
    y: f64,
}

fn dist(a: P, b: P) -> f64 {
    (a.x - b.x).hypot(a.y - b.y)
}

fn line_chains(geom: &Geometry) -> Vec<Vec<P>> {
    let to_pts = |cs: &[Coord]| -> Vec<P> {
        let mut out: Vec<P> = Vec::with_capacity(cs.len());
        for c in cs {
            let p = P { x: c.x, y: c.y };
            if out.last().is_none_or(|l| dist(*l, p) > 1e-12) {
                out.push(p);
            }
        }
        out
    };
    match geom {
        Geometry::LineString(cs) => vec![to_pts(cs)],
        Geometry::MultiLineString(lines) => lines.iter().map(|l| to_pts(l)).collect(),
        _ => Vec::new(),
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq)]
enum Orientation {
    AlongLine,
    Horizontal,
    Vertical,
}

struct Params {
    page_length: f64,
    page_width: f64,
    overlap: f64,
    orientation: Orientation,
    start_page: i64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let page_length = parse_optional_f64(args, "page_length")?.ok_or_else(|| {
        ToolError::Validation("required parameter 'page_length' is missing".into())
    })?;
    if !(page_length > 0.0 && page_length.is_finite()) {
        return Err(ToolError::Validation(
            "'page_length' must be a positive number".to_string(),
        ));
    }
    let page_width = parse_optional_f64(args, "page_width")?.ok_or_else(|| {
        ToolError::Validation("required parameter 'page_width' is missing".into())
    })?;
    if !(page_width > 0.0 && page_width.is_finite()) {
        return Err(ToolError::Validation(
            "'page_width' must be a positive number".to_string(),
        ));
    }
    let overlap = match parse_optional_f64(args, "overlap")? {
        None => 0.0,
        Some(v) if (0.0..1.0).contains(&v) => v,
        Some(_) => {
            return Err(ToolError::Validation(
                "'overlap' must be in [0, 1)".to_string(),
            ))
        }
    };
    let orientation =
        match parse_optional_str(args, "orientation")?.map(|s| s.trim().to_ascii_lowercase()) {
            None => Orientation::AlongLine,
            Some(s) if s == "along_line" => Orientation::AlongLine,
            Some(s) if s == "horizontal" => Orientation::Horizontal,
            Some(s) if s == "vertical" => Orientation::Vertical,
            Some(s) => {
                return Err(ToolError::Validation(format!(
                    "'orientation' must be 'along_line', 'horizontal', or 'vertical', got '{s}'"
                )))
            }
        };
    let start_page = parse_optional_i64(args, "start_page")?.unwrap_or(1);

    Ok(Params {
        page_length,
        page_width,
        overlap,
        orientation,
        start_page,
    })
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

fn parse_optional_i64(args: &ToolArgs, key: &str) -> Result<Option<i64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_i64().or_else(|| n.as_f64().map(|f| f as i64))),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<i64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be an integer"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be an integer"
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

    fn line_layer(coords: &[(f64, f64)]) -> String {
        let mut l = Layer::new("lines")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        let cs = coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
        l.add_feature(Some(Geometry::line_string(cs)), &[]).unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = StripMapIndexFeaturesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn poly_area(g: &Geometry) -> f64 {
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

    /// A straight horizontal line of length L with page_length p yields
    /// ceil(L / (p*(1-overlap))) pages, each an axis-aligned rectangle of the
    /// right width/length and rotation 90deg (east, along the line).
    #[test]
    fn straight_line_page_count_and_shape() {
        let input = line_layer(&[(0.0, 0.0), (1000.0, 0.0)]);
        let (out, layer) = run(json!({
            "input": input, "page_length": 300.0, "page_width": 100.0, "overlap": 0.0,
        }));
        // 1000 / 300 -> ceil = 4 pages.
        assert_eq!(out.outputs["page_count"], json!(4));
        assert_eq!(layer.features.len(), 4);
        for i in 0..3 {
            // Full interior pages: area = page_length * page_width.
            let a = poly_area(layer.features[i].geometry.as_ref().unwrap());
            assert!(
                (a - 300.0 * 100.0).abs() < 1e-6,
                "page {i} area {a} unexpected"
            );
        }
        let rot_i = layer.schema.field_index("rotation").unwrap();
        for f in &layer.features {
            let rot = f.attributes[rot_i].as_f64().unwrap();
            assert!(
                (rot - 90.0).abs() < 1e-6,
                "expected east bearing, got {rot}"
            );
        }
    }

    /// Overlap increases page count and the step between successive entry
    /// stations shrinks accordingly (pages advance by page_length*(1-overlap)).
    #[test]
    fn overlap_increases_page_count() {
        let input = line_layer(&[(0.0, 0.0), (1000.0, 0.0)]);
        let (out, _layer) = run(json!({
            "input": input, "page_length": 300.0, "page_width": 100.0, "overlap": 0.5,
        }));
        // step = 150 -> ceil(1000/150) = 7 pages.
        assert_eq!(out.outputs["page_count"], json!(7));
    }

    /// A diagonal 45deg line rotates each page to the chord bearing.
    #[test]
    fn diagonal_line_rotates_pages() {
        let input = line_layer(&[(0.0, 0.0), (1000.0, 1000.0)]);
        let (_out, layer) = run(json!({
            "input": input, "page_length": 200.0, "page_width": 50.0,
        }));
        let rot_i = layer.schema.field_index("rotation").unwrap();
        for f in &layer.features {
            let rot = f.attributes[rot_i].as_f64().unwrap();
            assert!(
                (rot - 45.0).abs() < 1e-6,
                "expected 45deg bearing, got {rot}"
            );
        }
    }

    /// horizontal/vertical orientation snaps rotation to the nearest axis.
    #[test]
    fn orientation_snaps_to_axis() {
        let input = line_layer(&[(0.0, 0.0), (1000.0, 1000.0)]); // 45deg chord
        let (_out, h_layer) = run(json!({
            "input": input, "page_length": 200.0, "page_width": 50.0, "orientation": "horizontal",
        }));
        let rot_i = h_layer.schema.field_index("rotation").unwrap();
        for f in &h_layer.features {
            let rot = f.attributes[rot_i].as_f64().unwrap();
            assert!(
                (rot - 90.0).abs() < 1e-6,
                "horizontal orientation should snap to 90deg, got {rot}"
            );
        }

        let (_out, v_layer) = run(json!({
            "input": input, "page_length": 200.0, "page_width": 50.0, "orientation": "vertical",
        }));
        for f in &v_layer.features {
            let rot = f.attributes[rot_i].as_f64().unwrap();
            assert!(
                (rot - 0.0).abs() < 1e-6,
                "vertical orientation should snap to 0deg, got {rot}"
            );
        }
    }

    /// Page numbers run sequentially from start_page, and line_id matches the
    /// source feature index.
    #[test]
    fn page_numbering_starts_at_start_page() {
        let input = line_layer(&[(0.0, 0.0), (600.0, 0.0)]);
        let (_out, layer) = run(json!({
            "input": input, "page_length": 300.0, "page_width": 100.0, "start_page": 5,
        }));
        let pn_i = layer.schema.field_index("page_number").unwrap();
        let lid_i = layer.schema.field_index("line_id").unwrap();
        let mut nums: Vec<i64> = layer
            .features
            .iter()
            .map(|f| f.attributes[pn_i].as_i64().unwrap())
            .collect();
        nums.sort_unstable();
        assert_eq!(nums, vec![5, 6]);
        for f in &layer.features {
            assert_eq!(f.attributes[lid_i].as_i64(), Some(0));
        }
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            StripMapIndexFeaturesTool.validate(&args)
        };
        assert!(bad(json!({})).is_err()); // missing input
        assert!(bad(json!({ "input": "a.geojson" })).is_err()); // missing page_length/page_width
        assert!(
            bad(json!({ "input": "a.geojson", "page_length": 0.0, "page_width": 10.0 })).is_err()
        );
        assert!(
            bad(json!({ "input": "a.geojson", "page_length": 10.0, "page_width": -1.0 })).is_err()
        );
        assert!(bad(
            json!({ "input": "a.geojson", "page_length": 10.0, "page_width": 5.0, "overlap": 1.0 })
        )
        .is_err());
        assert!(bad(json!({ "input": "a.geojson", "page_length": 10.0, "page_width": 5.0, "orientation": "diagonal" })).is_err());
        assert!(
            bad(json!({ "input": "a.geojson", "page_length": 10.0, "page_width": 5.0 })).is_ok()
        );
    }
}
