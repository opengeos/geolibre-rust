//! GeoLibre tool: build knockout mask polygons (and optional decoration wing
//! ticks) where an "above" line crosses a "below" line, so the lower line can be
//! masked at the crossing to read as a bridge/overpass.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Create Overpass* / *Create Underpass*
//! (Cartography). Given two line layers — the features that should appear to
//! pass *over* (`above`) and those masked *under* (`below`) — the tool finds
//! every place an above-line segment crosses a below-line segment and emits, at
//! each crossing, an oriented rectangular mask centered on the crossing point:
//!
//! - the rectangle is aligned with the local direction of the *above* line at
//!   the crossing, so the knockout follows the road/rail that stays on top;
//! - it extends `margin_along` in each direction along that line (total length
//!   `2 * margin_along`) and `margin_across` to each side (total width
//!   `2 * margin_across`), in CRS units.
//!
//! The mask is a plain rectangle (area exactly `4 * margin_along * margin_across`),
//! which a renderer draws in the background color on top of the below layer to
//! erase it under the overpass. Optionally, decorative *wing ticks* are emitted
//! as a companion line layer, controlled by `wing_type`:
//!
//! - `none`          — no decoration lines.
//! - `perpendicular` — a tick across each end of the mask (the classic
//!   bridge-parapet ticks), perpendicular to the above line.
//! - `parallel`      — two casing lines running the length of the mask, offset
//!   `margin_across` to each side of the above line.
//!
//! Only proper crossings count (an above segment and a below segment whose
//! interiors intersect); shared endpoints and merely touching lines are ignored,
//! matching the cartographic intent (a bridge is where lines genuinely cross).
//! Non-line geometries in either input are skipped.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CreateOverpassTool;

impl Tool for CreateOverpassTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "create_overpass",
            display_name: "Create Overpass",
            summary: "Build knockout mask polygons (and optional wing-tick decoration lines) where an 'above' line crosses a 'below' line, so the lower line is masked at the crossing.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "above",
                    description: "Line vector of features that pass OVER the crossing (the mask is oriented along these). Format auto-detected, or in-memory handle.",
                    required: true,
                },
                ToolParamSpec {
                    name: "below",
                    description: "Line vector of features masked UNDER the crossing (the lines to be knocked out). Format auto-detected, or in-memory handle.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output path for the mask polygons (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_decoration",
                    description: "Optional output path for the wing-tick decoration lines. If omitted, stored in memory (still returned as 'output_decoration').",
                    required: false,
                },
                ToolParamSpec {
                    name: "margin_along",
                    description: "Half-length of each mask along the above line, in CRS units (mask length is 2x this). Default 1.0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "margin_across",
                    description: "Half-width of each mask across the above line, in CRS units (mask width is 2x this). Default 1.0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "wing_type",
                    description: "Decoration style: 'none', 'perpendicular' (default; a tick across each mask end) or 'parallel' (casing lines along each side).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["above", "below"] {
            if args
                .get(key)
                .and_then(Value::as_str)
                .map(str::trim)
                .unwrap_or("")
                .is_empty()
            {
                return Err(ToolError::Validation(format!(
                    "missing required string parameter '{key}'"
                )));
            }
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let above_path = required_str(args, "above")?;
        let below_path = required_str(args, "below")?;
        let output = parse_optional_str(args, "output")?;
        let output_decoration = parse_optional_str(args, "output_decoration")?;
        let prm = parse_params(args)?;

        let above = load_input_layer(above_path)?;
        let below = load_input_layer(below_path)?;
        let crs = above.crs.clone();

        // Flatten each layer to (feature_index, segment) pairs.
        let above_segs = collect_segments(&above);
        let below_segs = collect_segments(&below);
        ctx.progress.info(&format!(
            "{} above segment(s), {} below segment(s)",
            above_segs.len(),
            below_segs.len()
        ));

        // Every proper crossing of an above segment with a below segment.
        let mut crossings: Vec<Crossing> = Vec::new();
        for a in &above_segs {
            for b in &below_segs {
                if !seg_bbox_overlap(a, b) {
                    continue;
                }
                if let Some(point) = segment_intersection(a.p1, a.p2, b.p1, b.p2) {
                    // Local direction of the ABOVE line at the crossing.
                    let (dx, dy) = unit(a.p2.x - a.p1.x, a.p2.y - a.p1.y);
                    crossings.push(Crossing {
                        point,
                        ux: dx,
                        uy: dy,
                        above_fid: a.fid,
                        below_fid: b.fid,
                    });
                }
            }
        }
        ctx.progress
            .info(&format!("{} crossing(s) found", crossings.len()));

        // Build the mask polygon layer.
        let mut mask_layer = Layer::new("overpass_mask");
        mask_layer.add_field(FieldDef::new("above_fid", FieldType::Integer));
        mask_layer.add_field(FieldDef::new("below_fid", FieldType::Integer));
        mask_layer.add_field(FieldDef::new("cross_x", FieldType::Float));
        mask_layer.add_field(FieldDef::new("cross_y", FieldType::Float));
        mask_layer.crs = crs.clone();
        mask_layer.geom_type = Some(GeometryType::Polygon);

        // Build the decoration line layer.
        let mut deco_layer = Layer::new("overpass_decoration");
        deco_layer.add_field(FieldDef::new("above_fid", FieldType::Integer));
        deco_layer.add_field(FieldDef::new("below_fid", FieldType::Integer));
        deco_layer.crs = crs;
        deco_layer.geom_type = Some(GeometryType::MultiLineString);

        for c in &crossings {
            let (ring, wings) = build_mask(c, prm.margin_along, prm.margin_across, prm.wing_type);
            mask_layer
                .add_feature(
                    Some(Geometry::Polygon {
                        exterior: wbvector::Ring::new(ring),
                        interiors: Vec::new(),
                    }),
                    &[
                        ("above_fid", (c.above_fid as i64).into()),
                        ("below_fid", (c.below_fid as i64).into()),
                        ("cross_x", c.point.x.into()),
                        ("cross_y", c.point.y.into()),
                    ],
                )
                .map_err(|e| ToolError::Execution(format!("failed adding mask feature: {e}")))?;

            if !wings.is_empty() {
                deco_layer
                    .add_feature(
                        Some(Geometry::MultiLineString(wings)),
                        &[
                            ("above_fid", (c.above_fid as i64).into()),
                            ("below_fid", (c.below_fid as i64).into()),
                        ],
                    )
                    .map_err(|e| {
                        ToolError::Execution(format!("failed adding decoration feature: {e}"))
                    })?;
            }
        }

        let crossing_count = crossings.len();
        let decoration_count = deco_layer.len();
        let mask_path = write_or_store_layer(mask_layer, output)?;
        let deco_path = write_or_store_layer(deco_layer, output_decoration)?;

        ctx.progress.info(&format!(
            "wrote {} mask polygon(s) and {} decoration feature(s)",
            crossing_count, decoration_count
        ));

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(mask_path));
        outputs.insert("output_decoration".to_string(), json!(deco_path));
        outputs.insert("crossing_count".to_string(), json!(crossing_count));
        outputs.insert("decoration_count".to_string(), json!(decoration_count));
        outputs.insert("wing_type".to_string(), json!(prm.wing_type.as_str()));
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum WingType {
    None,
    Perpendicular,
    Parallel,
}

impl WingType {
    fn as_str(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Perpendicular => "perpendicular",
            Self::Parallel => "parallel",
        }
    }
}

struct Params {
    margin_along: f64,
    margin_across: f64,
    wing_type: WingType,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let margin_along = parse_optional_f64(args, "margin_along")?.unwrap_or(1.0);
    let margin_across = parse_optional_f64(args, "margin_across")?.unwrap_or(1.0);
    for (name, v) in [
        ("margin_along", margin_along),
        ("margin_across", margin_across),
    ] {
        if !(v > 0.0 && v.is_finite()) {
            return Err(ToolError::Validation(format!(
                "parameter '{name}' must be a positive number"
            )));
        }
    }
    let wing_type = match parse_optional_str(args, "wing_type")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("perpendicular") => WingType::Perpendicular,
        Some("none") => WingType::None,
        Some("parallel") => WingType::Parallel,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown wing_type '{other}' (expected none, perpendicular or parallel)"
            )))
        }
    };
    Ok(Params {
        margin_along,
        margin_across,
        wing_type,
    })
}

/// Parses an optional numeric parameter, accepting a JSON number or a numeric
/// string (host UIs often post form values as strings).
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

fn required_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

// ── Geometry helpers ──────────────────────────────────────────────────────────

#[derive(Clone, Copy)]
struct Pt {
    x: f64,
    y: f64,
}

/// A single line segment tagged with its source feature index.
struct Seg {
    p1: Pt,
    p2: Pt,
    fid: usize,
}

/// Flattens every line feature of a layer into its constituent segments.
fn collect_segments(layer: &Layer) -> Vec<Seg> {
    let mut segs = Vec::new();
    for (fid, feature) in layer.features.iter().enumerate() {
        let Some(geom) = feature.geometry.as_ref() else {
            continue;
        };
        match geom {
            Geometry::LineString(cs) => push_line(cs, fid, &mut segs),
            Geometry::MultiLineString(lines) => {
                for cs in lines {
                    push_line(cs, fid, &mut segs);
                }
            }
            _ => {} // non-line geometries are ignored
        }
    }
    segs
}

fn push_line(cs: &[Coord], fid: usize, segs: &mut Vec<Seg>) {
    for w in cs.windows(2) {
        segs.push(Seg {
            p1: Pt {
                x: w[0].x,
                y: w[0].y,
            },
            p2: Pt {
                x: w[1].x,
                y: w[1].y,
            },
            fid,
        });
    }
}

/// Cheap axis-aligned bbox rejection before the exact intersection test.
fn seg_bbox_overlap(a: &Seg, b: &Seg) -> bool {
    let (a_minx, a_maxx) = min_max(a.p1.x, a.p2.x);
    let (a_miny, a_maxy) = min_max(a.p1.y, a.p2.y);
    let (b_minx, b_maxx) = min_max(b.p1.x, b.p2.x);
    let (b_miny, b_maxy) = min_max(b.p1.y, b.p2.y);
    a_minx <= b_maxx && a_maxx >= b_minx && a_miny <= b_maxy && a_maxy >= b_miny
}

fn min_max(a: f64, b: f64) -> (f64, f64) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

/// Proper intersection of two segments. Returns the crossing point only when the
/// interiors cross (parameters strictly inside `(0, 1)` on both segments), so a
/// shared endpoint or a mere touch does not count as an overpass crossing.
fn segment_intersection(p1: Pt, p2: Pt, q1: Pt, q2: Pt) -> Option<Pt> {
    let r = (p2.x - p1.x, p2.y - p1.y);
    let s = (q2.x - q1.x, q2.y - q1.y);
    let denom = r.0 * s.1 - r.1 * s.0;
    if denom.abs() < 1e-12 {
        return None; // parallel or degenerate
    }
    let qp = (q1.x - p1.x, q1.y - p1.y);
    let t = (qp.0 * s.1 - qp.1 * s.0) / denom;
    let u = (qp.0 * r.1 - qp.1 * r.0) / denom;
    // Strictly interior on both segments (an endpoint touch is not a crossing).
    if t > 1e-9 && t < 1.0 - 1e-9 && u > 1e-9 && u < 1.0 - 1e-9 {
        Some(Pt {
            x: p1.x + t * r.0,
            y: p1.y + t * r.1,
        })
    } else {
        None
    }
}

/// Normalizes a vector; falls back to +x for a degenerate zero vector.
fn unit(dx: f64, dy: f64) -> (f64, f64) {
    let len = dx.hypot(dy);
    if len < 1e-12 {
        (1.0, 0.0)
    } else {
        (dx / len, dy / len)
    }
}

struct Crossing {
    point: Pt,
    ux: f64,
    uy: f64,
    above_fid: usize,
    below_fid: usize,
}

/// Builds the oriented rectangular mask ring (CCW, unclosed — `Ring` closes it)
/// and any decoration wing lines for one crossing.
fn build_mask(
    c: &Crossing,
    along: f64,
    across: f64,
    wing: WingType,
) -> (Vec<Coord>, Vec<Vec<Coord>>) {
    let (ux, uy) = (c.ux, c.uy);
    let (vx, vy) = (-uy, ux); // perpendicular (left of direction)
    let p = c.point;
    // corner = P + su*along*u + sv*across*v
    let corner = |su: f64, sv: f64| {
        Coord::xy(
            p.x + su * along * ux + sv * across * vx,
            p.y + su * along * uy + sv * across * vy,
        )
    };
    // CCW: back-right, front-right, front-left, back-left.
    let ring = vec![
        corner(-1.0, -1.0),
        corner(1.0, -1.0),
        corner(1.0, 1.0),
        corner(-1.0, 1.0),
    ];

    let wings = match wing {
        WingType::None => Vec::new(),
        WingType::Perpendicular => vec![
            // tick across each end of the mask
            vec![corner(1.0, -1.0), corner(1.0, 1.0)],
            vec![corner(-1.0, -1.0), corner(-1.0, 1.0)],
        ],
        WingType::Parallel => vec![
            // casing line down each side
            vec![corner(-1.0, -1.0), corner(1.0, -1.0)],
            vec![corner(-1.0, 1.0), corner(1.0, 1.0)],
        ],
    };
    (ring, wings)
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn line_layer(name: &str, lines: &[Vec<(f64, f64)>]) -> String {
        let mut layer = Layer::new(name);
        for l in lines {
            let coords: Vec<Coord> = l.iter().map(|(x, y)| Coord::xy(*x, *y)).collect();
            layer
                .add_feature(Some(Geometry::LineString(coords)), &[])
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CreateOverpassTool.run(&args, &ctx()).unwrap();
        let mask = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        let deco = load_input_layer(out.outputs["output_decoration"].as_str().unwrap()).unwrap();
        (out, mask, deco)
    }

    fn poly_area(g: &Geometry) -> f64 {
        // Shoelace over the exterior ring.
        let Geometry::Polygon { exterior, .. } = g else {
            return 0.0;
        };
        let cs = &exterior.0;
        let n = cs.len();
        let mut a = 0.0;
        for i in 0..n {
            let j = (i + 1) % n;
            a += cs[i].x * cs[j].y - cs[j].x * cs[i].y;
        }
        a.abs() / 2.0
    }

    /// A horizontal above line crossing a vertical below line yields one mask,
    /// centered on the crossing, with area exactly 4*along*across.
    #[test]
    fn single_crossing_mask_area_and_center() {
        let above = line_layer("above", &[vec![(-10.0, 0.0), (10.0, 0.0)]]);
        let below = line_layer("below", &[vec![(0.0, -10.0), (0.0, 10.0)]]);
        let (out, mask, _) = run_tool(json!({
            "above": above, "below": below,
            "margin_along": 3.0, "margin_across": 2.0,
        }));
        assert_eq!(out.outputs["crossing_count"], json!(1));
        assert_eq!(mask.len(), 1);
        let g = mask.features[0].geometry.as_ref().unwrap();
        assert!((poly_area(g) - 4.0 * 3.0 * 2.0).abs() < 1e-9);
        // Mask spans [-3,3] x [-2,2] around origin: bbox check.
        let cs = if let Geometry::Polygon { exterior, .. } = g {
            &exterior.0
        } else {
            panic!()
        };
        let minx = cs.iter().map(|c| c.x).fold(f64::INFINITY, f64::min);
        let maxx = cs.iter().map(|c| c.x).fold(f64::NEG_INFINITY, f64::max);
        let miny = cs.iter().map(|c| c.y).fold(f64::INFINITY, f64::min);
        let maxy = cs.iter().map(|c| c.y).fold(f64::NEG_INFINITY, f64::max);
        assert!((minx + 3.0).abs() < 1e-9 && (maxx - 3.0).abs() < 1e-9);
        assert!((miny + 2.0).abs() < 1e-9 && (maxy - 2.0).abs() < 1e-9);
    }

    /// The mask orients along the above line: a 45-degree above line rotates the
    /// rectangle so its long axis follows that direction (checked via extent).
    #[test]
    fn mask_orients_along_above_line() {
        let above = line_layer("above", &[vec![(-10.0, -10.0), (10.0, 10.0)]]);
        let below = line_layer("below", &[vec![(-10.0, 10.0), (10.0, -10.0)]]);
        let (_, mask, _) = run_tool(json!({
            "above": above, "below": below,
            "margin_along": 4.0, "margin_across": 1.0,
        }));
        let g = mask.features[0].geometry.as_ref().unwrap();
        // Area invariant holds regardless of rotation.
        assert!((poly_area(g) - 4.0 * 4.0 * 1.0).abs() < 1e-9);
        // Long axis is diagonal, so extent > the axis-aligned 8x2 would give.
        let cs = if let Geometry::Polygon { exterior, .. } = g {
            &exterior.0
        } else {
            panic!()
        };
        let maxx = cs.iter().map(|c| c.x).fold(f64::NEG_INFINITY, f64::max);
        // Corner reaches ~ 4/sqrt2 + 1/sqrt2 in x for a 45-deg orientation.
        assert!(maxx > 3.0, "expected diagonal orientation, got maxx={maxx}");
    }

    /// Two below lines crossing one above line produce two masks.
    #[test]
    fn multiple_crossings() {
        let above = line_layer("above", &[vec![(-20.0, 0.0), (20.0, 0.0)]]);
        let below = line_layer(
            "below",
            &[
                vec![(-5.0, -5.0), (-5.0, 5.0)],
                vec![(5.0, -5.0), (5.0, 5.0)],
            ],
        );
        let (out, mask, _) = run_tool(json!({ "above": above, "below": below }));
        assert_eq!(out.outputs["crossing_count"], json!(2));
        assert_eq!(mask.len(), 2);
    }

    /// Wing types: perpendicular -> 2 ticks; parallel -> 2 casing lines; none -> 0.
    #[test]
    fn wing_types_control_decoration() {
        let above = line_layer("above", &[vec![(-10.0, 0.0), (10.0, 0.0)]]);
        let below = line_layer("below", &[vec![(0.0, -10.0), (0.0, 10.0)]]);

        let (out_p, _, deco_p) = run_tool(json!({
            "above": above.clone(), "below": below.clone(), "wing_type": "perpendicular"
        }));
        assert_eq!(deco_p.len(), 1);
        if let Some(Geometry::MultiLineString(ls)) = deco_p.features[0].geometry.as_ref() {
            assert_eq!(ls.len(), 2);
        } else {
            panic!("expected multilinestring");
        }
        assert_eq!(out_p.outputs["decoration_count"], json!(1));

        let (_, _, deco_par) = run_tool(json!({
            "above": above.clone(), "below": below.clone(), "wing_type": "parallel"
        }));
        assert_eq!(deco_par.len(), 1);

        let (out_n, _, deco_n) = run_tool(json!({
            "above": above, "below": below, "wing_type": "none"
        }));
        assert_eq!(deco_n.len(), 0);
        assert_eq!(out_n.outputs["decoration_count"], json!(0));
    }

    /// Lines that only touch at an endpoint (T-junction), or are parallel and
    /// non-crossing, produce no mask.
    #[test]
    fn touching_or_parallel_produce_no_mask() {
        // T-junction: below line ends exactly on the above line.
        let above = line_layer("above", &[vec![(-10.0, 0.0), (10.0, 0.0)]]);
        let below = line_layer("below", &[vec![(0.0, 0.0), (0.0, 10.0)]]);
        let (out, mask, _) = run_tool(json!({ "above": above, "below": below }));
        assert_eq!(out.outputs["crossing_count"], json!(0));
        assert_eq!(mask.len(), 0);

        // Parallel, offset lines never cross.
        let a2 = line_layer("above", &[vec![(-10.0, 0.0), (10.0, 0.0)]]);
        let b2 = line_layer("below", &[vec![(-10.0, 5.0), (10.0, 5.0)]]);
        let (out2, mask2, _) = run_tool(json!({ "above": a2, "below": b2 }));
        assert_eq!(out2.outputs["crossing_count"], json!(0));
        assert_eq!(mask2.len(), 0);
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = CreateOverpassTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing inputs must fail");
        assert!(
            bad(json!({ "above": "a.geojson" })).is_err(),
            "missing below"
        );
        assert!(
            bad(json!({ "above": "a.geojson", "below": "b.geojson", "margin_along": 0 })).is_err()
        );
        assert!(
            bad(json!({ "above": "a.geojson", "below": "b.geojson", "margin_across": -1 }))
                .is_err()
        );
        assert!(
            bad(json!({ "above": "a.geojson", "below": "b.geojson", "wing_type": "bogus" }))
                .is_err()
        );
        assert!(bad(json!({ "above": "a.geojson", "below": "b.geojson" })).is_ok());
        assert!(
            bad(json!({ "above": "a.geojson", "below": "b.geojson", "margin_along": "2.5" }))
                .is_ok(),
            "numeric strings ok"
        );
    }
}
