//! GeoLibre tool: subdivide polygons into equal-area (or target-area) parts.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Subdivide Polygon* (Data Management).
//! Nothing bundled splits a polygon into equal-area pieces. Used for parcel
//! pre-division, sampling-frame construction, and splitting oversized polygons
//! for parallel processing or tiling — a natural companion to the repo's tiling
//! identity.
//!
//! Each input polygon is cut into parallel strips at a chosen `angle`:
//!
//! 1. Rotate the polygon by `-angle` so the cuts become vertical.
//! 2. Walk left→right; for each target cumulative area, **binary-search** the
//!    cut position `x` where the area of the polygon left of `x` equals the
//!    target (the area-left-of-`x` function is monotonic, so bisection converges
//!    fast). Each strip is the polygon clipped between successive cuts via `geo`
//!    `BooleanOps` intersection with a half-plane rectangle.
//! 3. Rotate the strips back.
//!
//! `method = equal_parts` makes `num_parts` strips of area `total/num_parts`;
//! `method = equal_areas` makes strips of `target_area` each with the remainder
//! in the last strip. Attributes are copied to every strip, which also carries a
//! `part` index and its `part_area`. Only the parallel-strip subdivision is
//! implemented (ArcGIS's stacked-blocks layout is not). Distances/areas are in
//! the layer's CRS units — use a projected (equal-area) CRS for true equal areas.

use std::collections::BTreeMap;

use geo::{Area, BooleanOps, BoundingRect, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Binary-search iterations for each cut (2^-50 of the width is plenty).
const BISECT_ITERS: usize = 50;
/// Upper bound on parts per polygon, to bound runtime for tiny target areas.
const MAX_PARTS: usize = 10_000;

pub struct SubdividePolygonTool;

impl Tool for SubdividePolygonTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "subdivide_polygon",
            display_name: "Subdivide Polygon",
            summary: "Divide each polygon into equal-area parts (a set number of parts, or a target area each) using straight parallel cuts at a given angle, like ArcGIS Subdivide Polygon.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon vector layer (use a projected/equal-area CRS for true equal areas).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'equal_parts' (num_parts strips of equal area, default) or 'equal_areas' (strips of target_area, remainder in the last).",
                    required: false,
                },
                ToolParamSpec {
                    name: "num_parts",
                    description: "Number of equal-area parts (for method=equal_parts). Default 2.",
                    required: false,
                },
                ToolParamSpec {
                    name: "target_area",
                    description: "Area of each part in CRS units (required for method=equal_areas).",
                    required: false,
                },
                ToolParamSpec {
                    name: "angle",
                    description: "Cut orientation in degrees (0 = vertical cuts / horizontal strips). Default 0.",
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

        let mut out = Layer::new("subdivided").with_geom_type(GeometryType::MultiPolygon);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for field in layer.schema.fields() {
            out.add_field(field.clone());
        }
        out.add_field(FieldDef::new("part", FieldType::Integer));
        out.add_field(FieldDef::new("part_area", FieldType::Float));

        let (theta_cos, theta_sin) = {
            let r = prm.angle.to_radians();
            (r.cos(), r.sin())
        };

        let mut total_parts = 0usize;
        let mut subdivided = 0usize;
        let mut passthrough = 0usize;
        for feature in layer.features.iter() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let Some(mp) = to_multipolygon(geom) else {
                // Non-polygon: pass through with null part.
                emit(&mut out, feature, geom.clone(), -1, f64::NAN)?;
                passthrough += 1;
                continue;
            };
            let strips = subdivide(&mp, &prm, theta_cos, theta_sin);
            if strips.len() <= 1 {
                emit(&mut out, feature, geom.clone(), 0, mp.unsigned_area())?;
                passthrough += 1;
                continue;
            }
            for (k, strip) in strips.iter().enumerate() {
                let area = strip.unsigned_area();
                emit(
                    &mut out,
                    feature,
                    multipolygon_to_geometry(strip),
                    k as i64,
                    area,
                )?;
                total_parts += 1;
            }
            subdivided += 1;
        }

        ctx.progress.info(&format!(
            "{subdivided} polygon(s) subdivided into {total_parts} part(s); {passthrough} passed through"
        ));

        let feature_count = out.len();
        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("subdivided_count".to_string(), json!(subdivided));
        outputs.insert("part_count".to_string(), json!(total_parts));
        outputs.insert("passthrough_count".to_string(), json!(passthrough));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        Ok(ToolRunResult { outputs })
    }
}

/// Emits one output feature copying `src`'s attributes plus part/part_area.
fn emit(
    out: &mut Layer,
    src: &wbvector::Feature,
    geom: Geometry,
    part: i64,
    part_area: f64,
) -> Result<(), ToolError> {
    let mut attrs = src.attributes.clone();
    attrs.push(part.into());
    attrs.push(wbvector::FieldValue::Float(part_area));
    out.push(wbvector::Feature {
        fid: 0,
        geometry: Some(geom),
        attributes: attrs,
    });
    Ok(())
}

// ── Subdivision ──────────────────────────────────────────────────────────────

/// Cuts a polygon into parallel equal-area (or target-area) strips.
fn subdivide(mp: &MultiPolygon, prm: &Params, c: f64, s: f64) -> Vec<MultiPolygon> {
    let total = mp.unsigned_area();
    if total <= 0.0 {
        return vec![mp.clone()];
    }
    // Determine the target area per strip and the strip count.
    let (num, target) = match prm.method {
        Method::EqualParts => (prm.num_parts, total / prm.num_parts as f64),
        Method::EqualAreas => {
            let t = prm.target_area.unwrap_or(total);
            let n = (total / t).ceil() as usize;
            (n.clamp(1, MAX_PARTS), t)
        }
    };
    if num <= 1 {
        return vec![mp.clone()];
    }

    // Rotate into the cut frame (cuts become vertical lines x = const).
    let rot = rotate_mp(mp, c, s);
    let Some(rect) = rot.bounding_rect() else {
        return vec![mp.clone()];
    };
    let (min_x, max_x) = (rect.min().x, rect.max().x);
    let (min_y, max_y) = (rect.min().y, rect.max().y);
    let pad = (max_x - min_x + max_y - min_y).abs() + 1.0;

    // Cumulative-area targets at each internal cut.
    let mut cuts = vec![min_x - pad];
    for k in 1..num {
        let goal = (k as f64) * target;
        if goal >= total {
            break;
        }
        let x = bisect_area(&rot, min_x, max_x, min_y - pad, max_y + pad, pad, goal);
        cuts.push(x);
    }
    cuts.push(max_x + pad);

    // Clip each strip and rotate back (inverse rotation: cos, -sin).
    let mut strips = Vec::new();
    for w in cuts.windows(2) {
        let clip = axis_rect(w[0], w[1], min_y - pad, max_y + pad);
        let strip = rot.intersection(&clip);
        if strip.unsigned_area() > 0.0 {
            strips.push(rotate_mp(&strip, c, -s));
        }
    }
    if strips.is_empty() {
        vec![mp.clone()]
    } else {
        strips
    }
}

/// Binary-searches the vertical cut `x` in `[min_x, max_x]` at which the area of
/// `mp` left of `x` equals `goal`.
fn bisect_area(
    mp: &MultiPolygon,
    min_x: f64,
    max_x: f64,
    lo_y: f64,
    hi_y: f64,
    pad: f64,
    goal: f64,
) -> f64 {
    let mut lo = min_x;
    let mut hi = max_x;
    for _ in 0..BISECT_ITERS {
        let mid = 0.5 * (lo + hi);
        let left = axis_rect(min_x - pad, mid, lo_y, hi_y);
        let a = mp.intersection(&left).unsigned_area();
        if a < goal {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// An axis-aligned rectangle polygon spanning `[x0, x1] × [y0, y1]`.
fn axis_rect(x0: f64, x1: f64, y0: f64, y1: f64) -> MultiPolygon {
    let ext = LineString::new(vec![
        GeoCoord { x: x0, y: y0 },
        GeoCoord { x: x1, y: y0 },
        GeoCoord { x: x1, y: y1 },
        GeoCoord { x: x0, y: y1 },
        GeoCoord { x: x0, y: y0 },
    ]);
    MultiPolygon(vec![Polygon::new(ext, vec![])])
}

/// Rotates every coordinate of `mp` about the origin by the angle whose cosine
/// is `c` and sine is `s`.
fn rotate_mp(mp: &MultiPolygon, c: f64, s: f64) -> MultiPolygon {
    let rot_ls = |ls: &LineString| {
        LineString::new(
            ls.0.iter()
                .map(|p| GeoCoord {
                    x: p.x * c - p.y * s,
                    y: p.x * s + p.y * c,
                })
                .collect(),
        )
    };
    MultiPolygon(
        mp.0.iter()
            .map(|poly| {
                Polygon::new(
                    rot_ls(poly.exterior()),
                    poly.interiors().iter().map(rot_ls).collect(),
                )
            })
            .collect(),
    )
}

// ── geo <-> wbvector conversion ──────────────────────────────────────────────

fn to_multipolygon(geom: &Geometry) -> Option<MultiPolygon> {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => Some(MultiPolygon(vec![rings_to_polygon(exterior, interiors)])),
        Geometry::MultiPolygon(parts) => Some(MultiPolygon(
            parts.iter().map(|(e, i)| rings_to_polygon(e, i)).collect(),
        )),
        _ => None,
    }
}

fn rings_to_polygon(exterior: &wbvector::Ring, interiors: &[wbvector::Ring]) -> Polygon {
    Polygon::new(
        ring_to_linestring(exterior),
        interiors.iter().map(ring_to_linestring).collect(),
    )
}

fn ring_to_linestring(ring: &wbvector::Ring) -> LineString {
    LineString::new(
        ring.coords()
            .iter()
            .map(|c| GeoCoord { x: c.x, y: c.y })
            .collect(),
    )
}

fn multipolygon_to_geometry(mp: &MultiPolygon) -> Geometry {
    Geometry::MultiPolygon(mp.0.iter().map(polygon_to_rings).collect())
}

fn polygon_to_rings(poly: &Polygon) -> (wbvector::Ring, Vec<wbvector::Ring>) {
    (
        linestring_to_ring(poly.exterior()),
        poly.interiors().iter().map(linestring_to_ring).collect(),
    )
}

fn linestring_to_ring(ls: &LineString) -> wbvector::Ring {
    let mut coords: Vec<Coord> = ls.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
    if coords.len() >= 2 && coords.first() == coords.last() {
        coords.pop();
    }
    wbvector::Ring::new(coords)
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    EqualParts,
    EqualAreas,
}

struct Params {
    method: Method,
    num_parts: usize,
    target_area: Option<f64>,
    angle: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let method = match parse_optional_str(args, "method")? {
        None => Method::EqualParts,
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "equal_parts" | "number_of_equal_parts" => Method::EqualParts,
            "equal_areas" | "equal_area" => Method::EqualAreas,
            other => {
                return Err(ToolError::Validation(format!(
                    "'method' must be 'equal_parts' or 'equal_areas', got '{other}'"
                )))
            }
        },
    };
    let num_parts = match parse_optional_f64(args, "num_parts")? {
        None => 2,
        Some(v) if v.fract() == 0.0 && v >= 2.0 && v <= MAX_PARTS as f64 => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(format!(
                "'num_parts' must be an integer between 2 and {MAX_PARTS}"
            )))
        }
    };
    let target_area = parse_optional_f64(args, "target_area")?;
    if let Some(t) = target_area {
        if !(t > 0.0 && t.is_finite()) {
            return Err(ToolError::Validation(
                "'target_area' must be a positive number".to_string(),
            ));
        }
    }
    if method == Method::EqualAreas && target_area.is_none() {
        return Err(ToolError::Validation(
            "method=equal_areas requires 'target_area'".to_string(),
        ));
    }
    let angle = parse_optional_f64(args, "angle")?.unwrap_or(0.0);
    if !angle.is_finite() {
        return Err(ToolError::Validation(
            "'angle' must be a finite number".to_string(),
        ));
    }
    Ok(Params {
        method,
        num_parts,
        target_area,
        angle,
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

    fn rect(w: f64, h: f64) -> Geometry {
        Geometry::polygon(
            vec![
                Coord::xy(0.0, 0.0),
                Coord::xy(w, 0.0),
                Coord::xy(w, h),
                Coord::xy(0.0, h),
            ],
            vec![],
        )
    }

    fn layer_of(g: Geometry) -> String {
        let mut l = Layer::new("polys")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.add_feature(Some(g), &[("name", "p".into())]).unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SubdividePolygonTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn part_areas(layer: &Layer) -> Vec<f64> {
        let idx = layer.schema.field_index("part_area").unwrap();
        layer
            .iter()
            .map(|f| f.attributes[idx].as_f64().unwrap())
            .collect()
    }

    /// A 100×10 rectangle into 4 equal parts -> four 250-area strips.
    #[test]
    fn four_equal_parts() {
        let input = layer_of(rect(100.0, 10.0));
        let (out, layer) = run(json!({ "input": input, "num_parts": 4 }));
        assert_eq!(out.outputs["part_count"], json!(4));
        let areas = part_areas(&layer);
        assert_eq!(areas.len(), 4);
        for a in &areas {
            assert!((a - 250.0).abs() < 1e-3, "part area {a} != 250");
        }
        // Total area conserved.
        let total: f64 = areas.iter().sum();
        assert!((total - 1000.0).abs() < 1e-3, "total {total} != 1000");
    }

    /// Cutting at 90° gives the same equal areas (strips run the other way).
    #[test]
    fn angle_rotates_cut_direction() {
        let input = layer_of(rect(100.0, 40.0));
        let (out, layer) = run(json!({ "input": input, "num_parts": 5, "angle": 90.0 }));
        assert_eq!(out.outputs["part_count"], json!(5));
        for a in part_areas(&layer) {
            assert!((a - 800.0).abs() < 1e-2, "part area {a} != 800");
        }
    }

    /// equal_areas with a target splits into ceil(total/target) parts.
    #[test]
    fn equal_areas_target() {
        let input = layer_of(rect(100.0, 10.0)); // area 1000
        let (out, layer) = run(json!({
            "input": input, "method": "equal_areas", "target_area": 300.0,
        }));
        // ceil(1000/300) = 4 parts: 300,300,300,100.
        assert_eq!(out.outputs["part_count"], json!(4));
        let mut areas = part_areas(&layer);
        areas.sort_by(f64::total_cmp);
        assert!(
            (areas[0] - 100.0).abs() < 1e-2,
            "remainder {} != 100",
            areas[0]
        );
        for a in &areas[1..] {
            assert!((a - 300.0).abs() < 1e-2, "part {a} != 300");
        }
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SubdividePolygonTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "num_parts": 1 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "method": "equal_areas" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "method": "spiral" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "num_parts": 3 })).is_ok());
    }
}
