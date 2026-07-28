//! GeoLibre tool: underpass masks and trimmed lines at line crossings.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Create Underpass* (Cartography).
//!
//! GeoLibre already ships `create_overpass`, which handles the mirror case:
//! masking the *lower* feature so the *upper* one reads as continuous. Underpass
//! is its sibling and was simply absent.
//!
//! The two are not interchangeable. `create_overpass` produces a mask that hides
//! the bottom line under the top line's casing. An underpass instead needs the
//! **decorative gap**: the lower line is cut back on both sides of the crossing
//! so the break itself communicates that the feature dives under. A road network
//! styled with only `create_overpass` renders every grade separation
//! identically; this is the tool that distinguishes a tunnel portal from a
//! bridge.
//!
//! Two outputs are produced:
//!
//! * `output` — the mask polygons, oriented along the *above* line, one per
//!   resolved crossing, carrying both participants' feature ids;
//! * `output_lines` — the `below` layer with the underpass gaps actually cut
//!   out, ready to draw. This is the step that has no counterpart in
//!   `create_overpass`.
//!
//! Scope note: ArcGIS's `where_clause` attribute filter is not implemented;
//! filter the `above` layer upstream instead.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Builds underpass masks and cuts the corresponding gaps into the lower lines.
pub struct CreateUnderpassTool;

impl Tool for CreateUnderpassTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "create_underpass",
            display_name: "Create Underpass",
            summary: "At line crossings, builds mask polygons oriented along the upper feature and cuts a matching gap out of the lower feature so it reads as passing beneath (ArcGIS Create Underpass). Sibling of the shipped create_overpass, which masks the lower line so the upper one reads continuous; an underpass instead needs the decorative break, which is what distinguishes a tunnel portal from a bridge. Emits both the masks and the trimmed lower lines.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "above",
                    description: "Line vector of features that pass OVER the crossing (masks are oriented along these). Format auto-detected, or in-memory handle.",
                    required: true,
                },
                ToolParamSpec {
                    name: "below",
                    description: "Line vector of features that pass UNDER the crossing (the lines that get the gap). Format auto-detected, or in-memory handle.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output path for the underpass mask polygons. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_lines",
                    description: "Optional output path for the below layer with underpass gaps cut out. If omitted, stored in memory (still returned as 'output_lines').",
                    required: false,
                },
                ToolParamSpec {
                    name: "margin_along",
                    description: "Half-length of each mask along the above line, in CRS units (the gap is 2x this). Default 1.0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "margin_across",
                    description: "Half-width of each mask across the above line, in CRS units. Default 1.0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_angle",
                    description: "Skip crossings where the two lines meet at less than this angle in degrees, since a near-parallel mask is unstable. Default 15.",
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
        let output_lines = parse_optional_str(args, "output_lines")?;
        let prm = parse_params(args)?;

        let above = load_input_layer(above_path)?;
        let below = load_input_layer(below_path)?;
        let crs = above.crs.clone();

        let above_segs = collect_segments(&above);
        let below_segs = collect_segments(&below);
        ctx.progress.info(&format!(
            "{} above segment(s), {} below segment(s)",
            above_segs.len(),
            below_segs.len()
        ));

        // Every proper crossing of an above segment with a below segment.
        let min_cos = prm.min_angle.to_radians().cos();
        let mut masks: Vec<Mask> = Vec::new();
        let mut skipped_shallow = 0_u64;
        for a in &above_segs {
            for b in &below_segs {
                if !seg_bbox_overlap(a, b) {
                    continue;
                }
                let Some(point) = segment_intersection(a.p1, a.p2, b.p1, b.p2) else {
                    continue;
                };
                let (ux, uy) = unit(a.p2.x - a.p1.x, a.p2.y - a.p1.y);
                let (bx, by) = unit(b.p2.x - b.p1.x, b.p2.y - b.p1.y);
                // Near-parallel crossings give an unstable mask orientation, so
                // skip them rather than emit a degenerate rectangle.
                if (ux * bx + uy * by).abs() > min_cos {
                    skipped_shallow += 1;
                    continue;
                }
                masks.push(Mask {
                    cx: point.x,
                    cy: point.y,
                    ux,
                    uy,
                    along: prm.margin_along,
                    across: prm.margin_across,
                    above_fid: a.fid,
                    below_fid: b.fid,
                });
            }
        }
        ctx.progress
            .info(&format!("{} crossing(s) masked", masks.len()));

        // ── Mask polygons ────────────────────────────────────────────────────
        let mut mask_layer = Layer::new("underpass_mask");
        mask_layer.add_field(FieldDef::new("above_fid", FieldType::Integer));
        mask_layer.add_field(FieldDef::new("below_fid", FieldType::Integer));
        mask_layer.add_field(FieldDef::new("cross_x", FieldType::Float));
        mask_layer.add_field(FieldDef::new("cross_y", FieldType::Float));
        mask_layer.crs = crs.clone();
        mask_layer.geom_type = Some(GeometryType::Polygon);

        for m in &masks {
            mask_layer
                .add_feature(
                    Some(Geometry::Polygon {
                        exterior: Ring::new(m.ring()),
                        interiors: Vec::new(),
                    }),
                    &[
                        ("above_fid", (m.above_fid as i64).into()),
                        ("below_fid", (m.below_fid as i64).into()),
                        ("cross_x", m.cx.into()),
                        ("cross_y", m.cy.into()),
                    ],
                )
                .map_err(|e| ToolError::Execution(format!("failed adding mask feature: {e}")))?;
        }

        // ── Trimmed below lines ──────────────────────────────────────────────
        ctx.progress.info("cutting underpass gaps");
        let mut lines_layer = Layer::new("underpass_lines");
        for field in below.schema.fields().iter() {
            lines_layer.add_field(field.clone());
        }
        lines_layer.add_field(FieldDef::new("gap_count", FieldType::Integer));
        lines_layer.crs = crs;
        lines_layer.geom_type = Some(GeometryType::MultiLineString);

        let mut total_gaps = 0_u64;
        let mut trimmed_features = 0_u64;

        for (fid, feature) in below.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            // Only the masks generated by this feature's own crossings apply.
            let applicable: Vec<&Mask> = masks.iter().filter(|m| m.below_fid == fid).collect();

            let parts: Vec<Vec<Coord>> = match geom {
                Geometry::LineString(cs) => cut_line(cs, &applicable),
                Geometry::MultiLineString(lines) => lines
                    .iter()
                    .flat_map(|cs| cut_line(cs, &applicable))
                    .collect(),
                // Non-line geometry passes through untouched.
                _ => {
                    let mut attrs: Vec<(&str, wbvector::FieldValue)> = below
                        .schema
                        .fields()
                        .iter()
                        .enumerate()
                        .map(|(i, f)| (f.name.as_str(), feature.attributes[i].clone()))
                        .collect();
                    attrs.push(("gap_count", 0_i64.into()));
                    lines_layer
                        .add_feature(Some(geom.clone()), &attrs)
                        .map_err(|e| {
                            ToolError::Execution(format!("failed adding pass-through feature: {e}"))
                        })?;
                    continue;
                }
            };

            if parts.is_empty() {
                // Entirely consumed by masks — nothing left to draw.
                continue;
            }
            // A line split into N parts has N-1 gaps cut into it.
            let gaps = parts.len().saturating_sub(1) as i64;
            if gaps > 0 {
                total_gaps += gaps as u64;
                trimmed_features += 1;
            }

            let mut attrs: Vec<(&str, wbvector::FieldValue)> = below
                .schema
                .fields()
                .iter()
                .enumerate()
                .map(|(i, f)| (f.name.as_str(), feature.attributes[i].clone()))
                .collect();
            attrs.push(("gap_count", gaps.into()));

            lines_layer
                .add_feature(Some(Geometry::MultiLineString(parts)), &attrs)
                .map_err(|e| {
                    ToolError::Execution(format!("failed adding trimmed line feature: {e}"))
                })?;
        }

        let mask_count = masks.len();
        let mask_path = write_or_store_layer(mask_layer, output)?;
        let lines_path = write_or_store_layer(lines_layer, output_lines)?;

        ctx.progress.info(&format!(
            "wrote {mask_count} mask polygon(s) and cut {total_gaps} gap(s)"
        ));

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(mask_path));
        outputs.insert("output_lines".to_string(), json!(lines_path));
        outputs.insert("crossing_count".to_string(), json!(mask_count));
        outputs.insert("gap_count".to_string(), json!(total_gaps));
        outputs.insert("trimmed_feature_count".to_string(), json!(trimmed_features));
        outputs.insert("skipped_shallow_count".to_string(), json!(skipped_shallow));
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    margin_along: f64,
    margin_across: f64,
    min_angle: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let margin_along = parse_optional_f64(args, "margin_along")?.unwrap_or(1.0);
    let margin_across = parse_optional_f64(args, "margin_across")?.unwrap_or(1.0);
    let min_angle = parse_optional_f64(args, "min_angle")?.unwrap_or(15.0);
    if margin_along <= 0.0 || margin_across <= 0.0 {
        return Err(ToolError::Validation(
            "'margin_along' and 'margin_across' must be positive".to_string(),
        ));
    }
    if !(0.0..90.0).contains(&min_angle) {
        return Err(ToolError::Validation(
            "'min_angle' must be in [0, 90) degrees".to_string(),
        ));
    }
    Ok(Params {
        margin_along,
        margin_across,
        min_angle,
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

fn required_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

// ── Geometry ──────────────────────────────────────────────────────────────────

#[derive(Clone, Copy)]
struct Pt {
    x: f64,
    y: f64,
}

struct Seg {
    p1: Pt,
    p2: Pt,
    fid: usize,
}

/// An oriented rectangular mask centred on a crossing.
struct Mask {
    cx: f64,
    cy: f64,
    ux: f64,
    uy: f64,
    along: f64,
    across: f64,
    above_fid: usize,
    below_fid: usize,
}

impl Mask {
    /// CCW ring (unclosed — `Ring::new` closes it).
    fn ring(&self) -> Vec<Coord> {
        let (vx, vy) = (-self.uy, self.ux);
        let corner = |su: f64, sv: f64| {
            Coord::xy(
                self.cx + su * self.along * self.ux + sv * self.across * vx,
                self.cy + su * self.along * self.uy + sv * self.across * vy,
            )
        };
        vec![
            corner(-1.0, -1.0),
            corner(1.0, -1.0),
            corner(1.0, 1.0),
            corner(-1.0, 1.0),
        ]
    }

    /// Projects a point into the mask's local (along, across) frame.
    fn local(&self, x: f64, y: f64) -> (f64, f64) {
        let (dx, dy) = (x - self.cx, y - self.cy);
        let (vx, vy) = (-self.uy, self.ux);
        (dx * self.ux + dy * self.uy, dx * vx + dy * vy)
    }
}

fn collect_segments(layer: &Layer) -> Vec<Seg> {
    let mut segs = Vec::new();
    for (fid, feature) in layer.iter().enumerate() {
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
            _ => {}
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

/// Proper interior intersection of two segments (a shared endpoint does not
/// count as a crossing).
fn segment_intersection(p1: Pt, p2: Pt, q1: Pt, q2: Pt) -> Option<Pt> {
    let r = (p2.x - p1.x, p2.y - p1.y);
    let s = (q2.x - q1.x, q2.y - q1.y);
    let denom = r.0 * s.1 - r.1 * s.0;
    if denom.abs() < 1e-12 {
        return None;
    }
    let qp = (q1.x - p1.x, q1.y - p1.y);
    let t = (qp.0 * s.1 - qp.1 * s.0) / denom;
    let u = (qp.0 * r.1 - qp.1 * r.0) / denom;
    if t > 1e-9 && t < 1.0 - 1e-9 && u > 1e-9 && u < 1.0 - 1e-9 {
        Some(Pt {
            x: p1.x + t * r.0,
            y: p1.y + t * r.1,
        })
    } else {
        None
    }
}

fn unit(dx: f64, dy: f64) -> (f64, f64) {
    let len = dx.hypot(dy);
    if len < 1e-12 {
        (1.0, 0.0)
    } else {
        (dx / len, dy / len)
    }
}

/// Liang-Barsky clip of one segment against a mask rectangle, in the mask's
/// local frame. Returns the `[t0, t1]` sub-interval of the segment that lies
/// **inside** the rectangle, or `None` when it misses entirely.
fn segment_inside_interval(a: (f64, f64), b: (f64, f64), mask: &Mask) -> Option<(f64, f64)> {
    let (mut t0, mut t1) = (0.0_f64, 1.0_f64);
    let d = (b.0 - a.0, b.1 - a.1);
    // Four half-planes: along >= -L, along <= L, across >= -W, across <= W.
    let checks = [
        (-d.0, a.0 + mask.along),
        (d.0, mask.along - a.0),
        (-d.1, a.1 + mask.across),
        (d.1, mask.across - a.1),
    ];
    for (p, q) in checks {
        if p.abs() < 1e-15 {
            // Parallel to this boundary: reject only if already outside.
            if q < 0.0 {
                return None;
            }
            continue;
        }
        let r = q / p;
        if p < 0.0 {
            if r > t1 {
                return None;
            }
            if r > t0 {
                t0 = r;
            }
        } else {
            if r < t0 {
                return None;
            }
            if r < t1 {
                t1 = r;
            }
        }
    }
    if t1 <= t0 {
        None
    } else {
        Some((t0, t1))
    }
}

/// Cuts every mask rectangle out of a polyline, returning the surviving spans.
fn cut_line(coords: &[Coord], masks: &[&Mask]) -> Vec<Vec<Coord>> {
    if coords.len() < 2 {
        return Vec::new();
    }
    if masks.is_empty() {
        return vec![coords.to_vec()];
    }

    let mut parts: Vec<Vec<Coord>> = Vec::new();
    let mut current: Vec<Coord> = Vec::new();

    for w in coords.windows(2) {
        let (px, py) = (w[0].x, w[0].y);
        let (qx, qy) = (w[1].x, w[1].y);

        // Collect the covered sub-intervals of this segment.
        let mut covered: Vec<(f64, f64)> = Vec::new();
        for m in masks {
            let a = m.local(px, py);
            let b = m.local(qx, qy);
            if let Some(iv) = segment_inside_interval(a, b, m) {
                covered.push(iv);
            }
        }
        covered.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
        // Merge overlapping intervals so adjacent crossings cut one clean gap.
        let mut merged: Vec<(f64, f64)> = Vec::new();
        for iv in covered {
            match merged.last_mut() {
                Some(last) if iv.0 <= last.1 => last.1 = last.1.max(iv.1),
                _ => merged.push(iv),
            }
        }

        let lerp = |t: f64| Coord::xy(px + t * (qx - px), py + t * (qy - py));

        if current.is_empty() {
            current.push(Coord::xy(px, py));
        }
        let mut cursor = 0.0_f64;
        for (s, e) in merged {
            if s > cursor {
                // Keep the span before the gap, then close the part.
                current.push(lerp(s));
            }
            if !current.is_empty() && current.len() >= 2 {
                parts.push(std::mem::take(&mut current));
            } else {
                current.clear();
            }
            cursor = e.max(cursor);
            // Resume after the gap.
            if cursor < 1.0 {
                current.push(lerp(cursor));
            }
        }
        if cursor < 1.0 {
            current.push(Coord::xy(qx, qy));
        }
    }

    if current.len() >= 2 {
        parts.push(current);
    }
    // Drop degenerate zero-length remnants left at mask boundaries.
    parts.retain(|p| p.len() >= 2 && polyline_length(p) > 1e-12);
    parts
}

fn polyline_length(cs: &[Coord]) -> f64 {
    cs.windows(2)
        .map(|w| (w[1].x - w[0].x).hypot(w[1].y - w[0].y))
        .sum()
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

    fn line_layer(name: &str, lines: &[Vec<(f64, f64)>]) -> String {
        let mut l = Layer::new(name);
        l.geom_type = Some(GeometryType::LineString);
        for pts in lines {
            let cs: Vec<Coord> = pts.iter().map(|(x, y)| Coord::xy(*x, *y)).collect();
            l.add_feature(Some(Geometry::LineString(cs)), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(above: String, below: String, extra: Value) -> (Layer, Layer, ToolRunResult) {
        let mut obj = serde_json::Map::new();
        obj.insert("above".to_string(), json!(above));
        obj.insert("below".to_string(), json!(below));
        if let Value::Object(m) = extra {
            for (k, v) in m {
                obj.insert(k, v);
            }
        }
        let args: ToolArgs = serde_json::from_value(Value::Object(obj)).unwrap();
        let res = CreateUnderpassTool.run(&args, &ctx()).unwrap();
        let masks = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        let lines = load_input_layer(res.outputs["output_lines"].as_str().unwrap()).unwrap();
        (masks, lines, res)
    }

    /// A single perpendicular crossing yields one mask and one gap, splitting
    /// the below line into two parts. This is the behaviour that distinguishes
    /// underpass from overpass.
    #[test]
    fn crossing_cuts_a_gap_in_the_lower_line() {
        let above = line_layer("above", &[vec![(5.0, 0.0), (5.0, 10.0)]]);
        let below = line_layer("below", &[vec![(0.0, 5.0), (10.0, 5.0)]]);
        let (masks, lines, res) = run(
            above,
            below,
            json!({ "margin_along": 1.0, "margin_across": 1.0 }),
        );

        assert_eq!(masks.len(), 1, "one crossing -> one mask");
        assert_eq!(res.outputs["gap_count"], json!(1));
        assert_eq!(lines.len(), 1);

        let geom = lines.iter().next().unwrap().geometry.as_ref().unwrap();
        let Geometry::MultiLineString(parts) = geom else {
            panic!("expected a MultiLineString, got {geom:?}");
        };
        assert_eq!(parts.len(), 2, "the gap splits the line in two");

        // The gap spans the mask's across-width (1.0 either side of x=5).
        let left_end = parts[0].last().unwrap().x;
        let right_start = parts[1].first().unwrap().x;
        assert!(
            (left_end - 4.0).abs() < 1e-9,
            "left part ends at x=4, got {left_end}"
        );
        assert!(
            (right_start - 6.0).abs() < 1e-9,
            "right part resumes at x=6, got {right_start}"
        );
    }

    /// Total remaining length equals the original minus the gap width, so the
    /// cut removes exactly what the mask covers and nothing more.
    #[test]
    fn removed_length_matches_the_mask_width() {
        let above = line_layer("above", &[vec![(5.0, 0.0), (5.0, 10.0)]]);
        let below = line_layer("below", &[vec![(0.0, 5.0), (10.0, 5.0)]]);
        let (_, lines, _) = run(above, below, json!({ "margin_across": 1.5 }));

        let geom = lines.iter().next().unwrap().geometry.as_ref().unwrap();
        let Geometry::MultiLineString(parts) = geom else {
            panic!("expected MultiLineString");
        };
        let remaining: f64 = parts.iter().map(|p| polyline_length(p)).sum();
        // Original length 10; the mask is 1.5 either side of the crossing => 3.
        assert!(
            (remaining - 7.0).abs() < 1e-9,
            "expected 7.0 of line to survive, got {remaining}"
        );
    }

    /// Two crossings on one line cut two gaps, leaving three parts.
    #[test]
    fn multiple_crossings_cut_multiple_gaps() {
        let above = line_layer(
            "above",
            &[vec![(3.0, 0.0), (3.0, 10.0)], vec![(7.0, 0.0), (7.0, 10.0)]],
        );
        let below = line_layer("below", &[vec![(0.0, 5.0), (10.0, 5.0)]]);
        let (masks, lines, res) = run(above, below, json!({ "margin_across": 0.5 }));

        assert_eq!(masks.len(), 2);
        assert_eq!(res.outputs["gap_count"], json!(2));
        let geom = lines.iter().next().unwrap().geometry.as_ref().unwrap();
        let Geometry::MultiLineString(parts) = geom else {
            panic!("expected MultiLineString");
        };
        assert_eq!(parts.len(), 3);
    }

    /// Lines that do not cross are passed through whole, with no gap.
    #[test]
    fn non_crossing_lines_are_untouched() {
        let above = line_layer("above", &[vec![(0.0, 0.0), (10.0, 0.0)]]);
        let below = line_layer("below", &[vec![(0.0, 5.0), (10.0, 5.0)]]);
        let (masks, lines, res) = run(above, below, json!({}));

        assert_eq!(masks.len(), 0);
        assert_eq!(res.outputs["gap_count"], json!(0));
        let geom = lines.iter().next().unwrap().geometry.as_ref().unwrap();
        let Geometry::MultiLineString(parts) = geom else {
            panic!("expected MultiLineString");
        };
        assert_eq!(parts.len(), 1);
        assert!((polyline_length(&parts[0]) - 10.0).abs() < 1e-9);
    }

    /// Near-parallel crossings are skipped: the mask orientation would be
    /// unstable and the resulting gap meaningless.
    #[test]
    fn shallow_crossings_are_skipped() {
        // Two nearly-parallel lines that still technically cross.
        let above = line_layer("above", &[vec![(0.0, 0.0), (10.0, 0.2)]]);
        let below = line_layer("below", &[vec![(0.0, 0.2), (10.0, 0.0)]]);
        let (masks, _, res) = run(above, below, json!({ "min_angle": 15.0 }));
        assert_eq!(masks.len(), 0, "shallow crossing must not produce a mask");
        assert!(res.outputs["skipped_shallow_count"].as_u64().unwrap() >= 1);
    }

    /// The mask rectangle is oriented along the ABOVE line, so its long axis
    /// follows the upper feature regardless of the lower one's direction.
    #[test]
    fn mask_is_oriented_along_the_above_line() {
        // Above runs east-west, so the mask's along-axis is horizontal.
        let above = line_layer("above", &[vec![(0.0, 5.0), (10.0, 5.0)]]);
        let below = line_layer("below", &[vec![(5.0, 0.0), (5.0, 10.0)]]);
        let (masks, _, _) = run(
            above,
            below,
            json!({ "margin_along": 3.0, "margin_across": 0.5 }),
        );
        let geom = masks.iter().next().unwrap().geometry.as_ref().unwrap();
        let Geometry::Polygon { exterior, .. } = geom else {
            panic!("expected Polygon");
        };
        let xs: Vec<f64> = exterior.0.iter().map(|c| c.x).collect();
        let ys: Vec<f64> = exterior.0.iter().map(|c| c.y).collect();
        let width = xs.iter().cloned().fold(f64::MIN, f64::max)
            - xs.iter().cloned().fold(f64::MAX, f64::min);
        let height = ys.iter().cloned().fold(f64::MIN, f64::max)
            - ys.iter().cloned().fold(f64::MAX, f64::min);
        assert!(
            (width - 6.0).abs() < 1e-9,
            "along-axis spans 2*3, got {width}"
        );
        assert!(
            (height - 1.0).abs() < 1e-9,
            "across-axis spans 2*0.5, got {height}"
        );
    }

    #[test]
    fn rejects_bad_parameters() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(CreateUnderpassTool.validate(&args).is_err());

        let a = line_layer("a", &[vec![(0.0, 0.0), (1.0, 1.0)]]);
        let b = line_layer("b", &[vec![(0.0, 1.0), (1.0, 0.0)]]);
        for bad in [
            json!({ "above": a.clone(), "below": b.clone(), "margin_along": -1.0 }),
            json!({ "above": a.clone(), "below": b.clone(), "margin_across": 0.0 }),
            json!({ "above": a.clone(), "below": b.clone(), "min_angle": 95.0 }),
            json!({ "above": a.clone(), "below": b.clone(), "margin_along": "wide" }),
        ] {
            let args: ToolArgs = serde_json::from_value(bad).unwrap();
            assert!(CreateUnderpassTool.validate(&args).is_err());
        }
    }
}
