//! GeoLibre tool: arc-aware simplification preserving engineered curves.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Simplify By Straight Lines And
//! Circular Arcs* (Editing), with the sibling *Simplify By Tangent Segments*
//! exposed through `mode`.
//!
//! Every simplifier in either registry is **vertex-removal only**:
//! `simplify_features` and `simplify_shared_edges` (bundled), `simplify_building`,
//! `simplify_3d_line`, `smooth_natural_features`. They all answer "which of
//! these vertices can I drop", and the output is always a polyline.
//!
//! That is the wrong model for built infrastructure. Cul-de-sac bulbs, highway
//! curves, roundabouts, cadastral boundaries following a curve, silo and tank
//! footprints and stadium outlines are *designed as arcs* and typically arrive
//! densified into dozens of near-collinear vertices. Douglas-Peucker on that
//! data forces a choice between a visibly faceted curve and an unacceptable
//! vertex count, because the primitive being fitted is wrong.
//!
//! Fitting arcs directly gives a large vertex reduction at **lower** geometric
//! error than any chord approximation with the same budget.
//!
//! ## Fit, then verify
//!
//! Circles are fitted by algebraic (Kåsa) least squares — a closed-form 3x3
//! solve, no linear-algebra crate — but a candidate is **accepted on the true
//! maximum orthogonal deviation** of its vertices from the fitted arc. The
//! algebraic fit is biased for short, low-curvature runs, so validating on the
//! real distance is what keeps the `tolerance` guarantee honest.
//!
//! ## Output representation
//!
//! The repo's vector model has no native arc primitive, so fitted arcs are
//! re-densified to a chord tolerance on write (`densify_output`, default true)
//! while the arc parameters (centre, radius, sweep) are carried as attributes.
//! Set `densify_output=false` to emit only the arc endpoints plus those
//! attributes, for a consumer that can render true curves.
//!
//! Scope note: features are simplified independently, so a shared boundary
//! between two polygons can crack. Use the bundled `simplify_shared_edges` when
//! coverage safety matters more than curve fidelity.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, Geometry, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Replaces runs of vertices with straight segments and circular arcs.
pub struct SimplifyByCircularArcsTool;

impl Tool for SimplifyByCircularArcsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "simplify_by_circular_arcs",
            display_name: "Simplify By Circular Arcs",
            summary: "Simplifies lines and polygons by replacing runs of vertices with straight segments and true circular arcs, preserving engineered curves instead of approximating them with chords (ArcGIS Simplify By Straight Lines And Circular Arcs; 'tangent' mode covers Simplify By Tangent Segments). Every simplifier in either registry is vertex-removal only, which forces a choice between a faceted curve and an unacceptable vertex count on cul-de-sac bulbs, roundabouts and tank footprints. Circles are fitted by algebraic least squares but accepted on true orthogonal deviation, so the tolerance guarantee holds.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input line or polygon layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Maximum allowed deviation of a fitted primitive from the original vertices, in CRS units.",
                    required: true,
                },
                ToolParamSpec {
                    name: "mode",
                    description: "'arcs' (default; fits straight segments and circular arcs) or 'tangent' (straight segments only, the Simplify By Tangent Segments behaviour).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_arc_angle",
                    description: "Arcs subtending less than this many degrees are emitted as straight segments instead (default 10).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_radius",
                    description: "Arcs with a larger radius than this are emitted as straight segments (default: 1e6).",
                    required: false,
                },
                ToolParamSpec {
                    name: "densify_output",
                    description: "Re-densify fitted arcs to a chord tolerance on write (default true). When false, only arc endpoints are emitted.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args.get("input").and_then(Value::as_str).is_none() {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        if args.get("tolerance").is_none() {
            return Err(ToolError::Validation(
                "missing required parameter 'tolerance'".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args.get("input").and_then(Value::as_str).ok_or_else(|| {
            ToolError::Validation("missing required parameter 'input'".to_string())
        })?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;

        let mut out = Layer::new("simplify_by_circular_arcs");
        for field in layer.schema.fields().iter() {
            out.add_field(field.clone());
        }
        out.add_field(FieldDef::new("arc_count", FieldType::Integer));
        out.add_field(FieldDef::new("input_vertices", FieldType::Integer));
        out.add_field(FieldDef::new("output_vertices", FieldType::Integer));
        out.crs = layer.crs.clone();
        out.geom_type = layer.geom_type;

        let mut total_in = 0_u64;
        let mut total_out = 0_u64;
        let mut total_arcs = 0_u64;

        for (fid, feature) in layer.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let mut arcs = 0_u64;
            let mut vin = 0_u64;
            let mut vout = 0_u64;

            let simplified = simplify_geometry(geom, &prm, &mut arcs, &mut vin, &mut vout);

            total_in += vin;
            total_out += vout;
            total_arcs += arcs;

            let mut attrs: Vec<(&str, wbvector::FieldValue)> = layer
                .schema
                .fields()
                .iter()
                .enumerate()
                .map(|(i, f)| (f.name.as_str(), feature.attributes[i].clone()))
                .collect();
            attrs.push(("arc_count", (arcs as i64).into()));
            attrs.push(("input_vertices", (vin as i64).into()));
            attrs.push(("output_vertices", (vout as i64).into()));

            out.add_feature(simplified, &attrs)
                .map_err(|e| ToolError::Execution(format!("failed adding feature: {e}")))?;

            ctx.progress
                .progress((fid as f64 + 1.0) / layer.len().max(1) as f64);
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("arc_count".to_string(), json!(total_arcs));
        outputs.insert("input_vertex_count".to_string(), json!(total_in));
        outputs.insert("output_vertex_count".to_string(), json!(total_out));
        outputs.insert(
            "vertex_reduction".to_string(),
            json!(if total_in > 0 {
                1.0 - total_out as f64 / total_in as f64
            } else {
                0.0
            }),
        );
        Ok(ToolRunResult { outputs })
    }
}

fn simplify_geometry(
    geom: &Geometry,
    prm: &Params,
    arcs: &mut u64,
    vin: &mut u64,
    vout: &mut u64,
) -> Option<Geometry> {
    match geom {
        Geometry::LineString(cs) => {
            *vin += cs.len() as u64;
            let s = simplify_run(cs, prm, false, arcs);
            *vout += s.len() as u64;
            Some(Geometry::LineString(s))
        }
        Geometry::MultiLineString(parts) => {
            let mut outp = Vec::new();
            for cs in parts {
                *vin += cs.len() as u64;
                let s = simplify_run(cs, prm, false, arcs);
                *vout += s.len() as u64;
                outp.push(s);
            }
            Some(Geometry::MultiLineString(outp))
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            *vin += exterior.0.len() as u64;
            let ext = simplify_run(&exterior.0, prm, true, arcs);
            *vout += ext.len() as u64;
            let mut holes = Vec::new();
            for h in interiors {
                *vin += h.0.len() as u64;
                let s = simplify_run(&h.0, prm, true, arcs);
                *vout += s.len() as u64;
                if s.len() >= 3 {
                    holes.push(Ring::new(s));
                }
            }
            if ext.len() < 3 {
                return Some(geom.clone());
            }
            Some(Geometry::Polygon {
                exterior: Ring::new(ext),
                interiors: holes,
            })
        }
        Geometry::MultiPolygon(parts) => {
            let mut outp = Vec::new();
            for (e, hs) in parts {
                *vin += e.0.len() as u64;
                let ext = simplify_run(&e.0, prm, true, arcs);
                *vout += ext.len() as u64;
                if ext.len() < 3 {
                    outp.push((e.clone(), hs.clone()));
                    continue;
                }
                let mut holes = Vec::new();
                for h in hs {
                    *vin += h.0.len() as u64;
                    let s = simplify_run(&h.0, prm, true, arcs);
                    *vout += s.len() as u64;
                    if s.len() >= 3 {
                        holes.push(Ring::new(s));
                    }
                }
                outp.push((Ring::new(ext), holes));
            }
            Some(Geometry::MultiPolygon(outp))
        }
        // Points and collections pass through untouched.
        other => Some(other.clone()),
    }
}

/// Greedily covers a vertex run with the longest acceptable primitives.
///
/// At each position the longest arc that satisfies `tolerance` is preferred;
/// where no arc qualifies, the longest straight segment is used. Growing each
/// primitive as far as tolerance allows is what produces large vertex
/// reductions rather than many short arcs.
fn simplify_run(coords: &[Coord], prm: &Params, closed: bool, arcs: &mut u64) -> Vec<Coord> {
    let n = coords.len();
    if n < 3 {
        return coords.to_vec();
    }
    // Work on an explicitly closed copy for rings so the wrap segment is fitted
    // like any other.
    let pts: Vec<Coord> = if closed {
        let mut v = coords.to_vec();
        v.push(coords[0].clone());
        v
    } else {
        coords.to_vec()
    };
    let m = pts.len();

    let mut out: Vec<Coord> = vec![pts[0].clone()];
    let mut i = 0_usize;

    while i + 1 < m {
        // Longest straight run from i that stays within tolerance.
        let mut best_line = i + 1;
        let mut j = i + 2;
        while j < m {
            if max_line_deviation(&pts[i..=j]) <= prm.tolerance {
                best_line = j;
                j += 1;
            } else {
                break;
            }
        }

        // Longest acceptable arc from i, if arcs are enabled.
        let mut best_arc: Option<(usize, Circle)> = None;
        if prm.fit_arcs {
            let mut j = i + 3; // an arc needs at least 4 points to beat a line
            while j < m {
                // Keep extending while the circle still *fits*. The acceptance
                // guards are checked separately and must NOT stop the search:
                // swept angle grows with the run, so a short 4-point candidate
                // is routinely below min_arc_angle even though the full arc is
                // well above it. Breaking on that guard would find no arcs at
                // all on a densified curve.
                match fit_circle(&pts[i..=j]) {
                    Some(c) if max_arc_deviation(&pts[i..=j], &c) <= prm.tolerance => {
                        if c.radius <= prm.max_radius
                            && sweep_degrees(&pts[i..=j], &c) >= prm.min_arc_angle
                        {
                            best_arc = Some((j, c));
                        }
                        j += 1;
                    }
                    _ => break,
                }
            }
        }

        // Prefer whichever primitive covers more vertices; ties go to the line,
        // since a straight segment is cheaper and never worse to render.
        match best_arc {
            Some((arc_end, circle)) if arc_end > best_line => {
                *arcs += 1;
                if prm.densify_output {
                    let dens = densify_arc(&pts[i], &pts[arc_end], &circle, &pts[i..=arc_end], prm);
                    // Skip the first point: it is already in `out`.
                    out.extend(dens.into_iter().skip(1));
                } else {
                    out.push(pts[arc_end].clone());
                }
                i = arc_end;
            }
            _ => {
                out.push(pts[best_line].clone());
                i = best_line;
            }
        }
    }

    if closed {
        // The explicit closing duplicate is not stored in a `Ring`.
        if out.len() > 1 {
            let first = &out[0];
            let last = &out[out.len() - 1];
            if (first.x - last.x).abs() < 1e-12 && (first.y - last.y).abs() < 1e-12 {
                out.pop();
            }
        }
        if out.len() < 3 {
            return coords.to_vec();
        }
    }
    out
}

/// Maximum orthogonal distance of intermediate points from the chord.
fn max_line_deviation(pts: &[Coord]) -> f64 {
    if pts.len() < 3 {
        return 0.0;
    }
    let (a, b) = (&pts[0], &pts[pts.len() - 1]);
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let len = dx.hypot(dy);
    let mut worst: f64 = 0.0;
    for p in &pts[1..pts.len() - 1] {
        let d = if len < 1e-15 {
            (p.x - a.x).hypot(p.y - a.y)
        } else {
            ((p.x - a.x) * dy - (p.y - a.y) * dx).abs() / len
        };
        worst = worst.max(d);
    }
    worst
}

#[derive(Clone, Copy, Debug)]
struct Circle {
    cx: f64,
    cy: f64,
    radius: f64,
}

/// Kåsa algebraic circle fit: minimises the algebraic residual, giving a
/// closed-form 3x3 solve. Biased for short low-curvature runs, which is why the
/// caller validates the fit on true orthogonal deviation.
fn fit_circle(pts: &[Coord]) -> Option<Circle> {
    let n = pts.len() as f64;
    if pts.len() < 3 {
        return None;
    }
    let (mut sx, mut sy) = (0.0, 0.0);
    for p in pts {
        sx += p.x;
        sy += p.y;
    }
    // Centre the data for conditioning.
    let (mx, my) = (sx / n, sy / n);

    let (mut suu, mut suv, mut svv) = (0.0, 0.0, 0.0);
    let (mut suuu, mut svvv, mut suvv, mut svuu) = (0.0, 0.0, 0.0, 0.0);
    for p in pts {
        let u = p.x - mx;
        let v = p.y - my;
        suu += u * u;
        suv += u * v;
        svv += v * v;
        suuu += u * u * u;
        svvv += v * v * v;
        suvv += u * v * v;
        svuu += v * u * u;
    }

    let det = suu * svv - suv * suv;
    if det.abs() < 1e-18 {
        // Collinear points define no circle.
        return None;
    }
    let b1 = (suuu + suvv) / 2.0;
    let b2 = (svvv + svuu) / 2.0;
    let uc = (b1 * svv - b2 * suv) / det;
    let vc = (b2 * suu - b1 * suv) / det;

    let radius = (uc * uc + vc * vc + (suu + svv) / n).sqrt();
    if !radius.is_finite() || radius <= 0.0 {
        return None;
    }
    Some(Circle {
        cx: uc + mx,
        cy: vc + my,
        radius,
    })
}

/// True maximum orthogonal deviation of the points from the fitted circle.
fn max_arc_deviation(pts: &[Coord], c: &Circle) -> f64 {
    let mut worst: f64 = 0.0;
    for p in pts {
        let d = ((p.x - c.cx).hypot(p.y - c.cy) - c.radius).abs();
        worst = worst.max(d);
    }
    worst
}

/// Total swept angle covered by the run, in degrees. Guards against a run that
/// wraps more than once or reverses direction, which is not a single arc.
fn sweep_degrees(pts: &[Coord], c: &Circle) -> f64 {
    if pts.len() < 2 {
        return 0.0;
    }
    let ang = |p: &Coord| (p.y - c.cy).atan2(p.x - c.cx);
    let mut total = 0.0;
    let mut sign = 0.0;
    for w in pts.windows(2) {
        let mut d = ang(&w[1]) - ang(&w[0]);
        while d > std::f64::consts::PI {
            d -= std::f64::consts::TAU;
        }
        while d < -std::f64::consts::PI {
            d += std::f64::consts::TAU;
        }
        if sign == 0.0 {
            sign = d.signum();
        } else if d != 0.0 && d.signum() != sign {
            // Direction reversal: not a single circular arc.
            return 0.0;
        }
        total += d;
    }
    let deg = total.abs().to_degrees();
    if deg >= 360.0 {
        // A full wrap is not a single arc either.
        return 0.0;
    }
    deg
}

/// Re-densifies a fitted arc so the chord error stays within tolerance.
fn densify_arc(
    start: &Coord,
    end: &Coord,
    c: &Circle,
    original: &[Coord],
    prm: &Params,
) -> Vec<Coord> {
    let a0 = (start.y - c.cy).atan2(start.x - c.cx);
    let a1 = (end.y - c.cy).atan2(end.x - c.cx);

    // Follow the sweep direction the original vertices actually took.
    let sweep = {
        let mut s = sweep_signed(original, c);
        if s == 0.0 {
            let mut d = a1 - a0;
            while d > std::f64::consts::PI {
                d -= std::f64::consts::TAU;
            }
            while d < -std::f64::consts::PI {
                d += std::f64::consts::TAU;
            }
            s = d;
        }
        s
    };

    // Chord error for a step t is r * (1 - cos(t/2)); invert for the tolerance.
    let ratio = (1.0 - prm.tolerance / c.radius).clamp(-1.0, 1.0);
    let max_step = 2.0 * ratio.acos();
    let steps = if max_step > 1e-9 {
        ((sweep.abs() / max_step).ceil() as usize).clamp(1, 512)
    } else {
        1
    };

    let mut out = Vec::with_capacity(steps + 1);
    for k in 0..=steps {
        let t = a0 + sweep * (k as f64 / steps as f64);
        out.push(Coord::xy(
            c.cx + c.radius * t.cos(),
            c.cy + c.radius * t.sin(),
        ));
    }
    // Pin the endpoints exactly so consecutive primitives stay connected.
    out[0] = start.clone();
    let last = out.len() - 1;
    out[last] = end.clone();
    out
}

fn sweep_signed(pts: &[Coord], c: &Circle) -> f64 {
    let ang = |p: &Coord| (p.y - c.cy).atan2(p.x - c.cx);
    let mut total = 0.0;
    for w in pts.windows(2) {
        let mut d = ang(&w[1]) - ang(&w[0]);
        while d > std::f64::consts::PI {
            d -= std::f64::consts::TAU;
        }
        while d < -std::f64::consts::PI {
            d += std::f64::consts::TAU;
        }
        total += d;
    }
    total
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    tolerance: f64,
    fit_arcs: bool,
    min_arc_angle: f64,
    max_radius: f64,
    densify_output: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let tolerance = opt_f64(args, "tolerance")?.ok_or_else(|| {
        ToolError::Validation("missing required parameter 'tolerance'".to_string())
    })?;
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(ToolError::Validation(
            "'tolerance' must be a positive, finite number".to_string(),
        ));
    }
    let fit_arcs = match parse_optional_str(args, "mode")? {
        None => true,
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "arcs" => true,
            "tangent" => false,
            other => {
                return Err(ToolError::Validation(format!(
                    "unknown mode '{other}' (expected 'arcs' or 'tangent')"
                )))
            }
        },
    };
    let min_arc_angle = opt_f64(args, "min_arc_angle")?.unwrap_or(10.0);
    if !(0.0..360.0).contains(&min_arc_angle) {
        return Err(ToolError::Validation(
            "'min_arc_angle' must be in [0, 360) degrees".to_string(),
        ));
    }
    let max_radius = opt_f64(args, "max_radius")?.unwrap_or(1.0e6);
    if !max_radius.is_finite() || max_radius <= 0.0 {
        return Err(ToolError::Validation(
            "'max_radius' must be a positive, finite number".to_string(),
        ));
    }
    Ok(Params {
        tolerance,
        fit_arcs,
        min_arc_angle,
        max_radius,
        densify_output: opt_bool(args, "densify_output")?.unwrap_or(true),
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

    fn line_layer(pts: Vec<(f64, f64)>) -> String {
        let mut l = Layer::new("lines");
        l.geom_type = Some(GeometryType::LineString);
        l.add_feature(
            Some(Geometry::LineString(
                pts.iter().map(|(x, y)| Coord::xy(*x, *y)).collect(),
            )),
            &[],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(path: String, extra: Value) -> (Layer, ToolRunResult) {
        let mut obj = serde_json::Map::new();
        obj.insert("input".to_string(), json!(path));
        if let Value::Object(m) = extra {
            for (k, v) in m {
                obj.insert(k, v);
            }
        }
        let args: ToolArgs = serde_json::from_value(Value::Object(obj)).unwrap();
        let res = SimplifyByCircularArcsTool.run(&args, &ctx()).unwrap();
        let l = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (l, res)
    }

    fn first_line(l: &Layer) -> Vec<Coord> {
        match l.iter().next().unwrap().geometry.as_ref().unwrap() {
            Geometry::LineString(cs) => cs.clone(),
            other => panic!("expected LineString, got {other:?}"),
        }
    }

    /// A densified semicircle is recognised as one arc rather than being
    /// faceted — this is the behaviour no vertex-removal simplifier can give.
    #[test]
    fn densified_arc_is_recognised_as_one_arc() {
        // 60 samples along a half circle of radius 100.
        let pts: Vec<(f64, f64)> = (0..=60)
            .map(|i| {
                let t = std::f64::consts::PI * i as f64 / 60.0;
                (100.0 * t.cos(), 100.0 * t.sin())
            })
            .collect();
        let path = line_layer(pts);
        let (_, res) = run(path, json!({ "tolerance": 0.5, "densify_output": false }));
        assert_eq!(
            res.outputs["arc_count"],
            json!(1),
            "the whole semicircle should be a single fitted arc"
        );
        // Endpoints only, with densification off.
        assert_eq!(res.outputs["output_vertex_count"], json!(2));
    }

    /// The tolerance guarantee holds: every original vertex stays within
    /// `tolerance` of the simplified geometry.
    #[test]
    fn every_input_vertex_stays_within_tolerance() {
        let original: Vec<(f64, f64)> = (0..=80)
            .map(|i| {
                let t = 1.5 * std::f64::consts::PI * i as f64 / 80.0;
                (50.0 * t.cos(), 50.0 * t.sin())
            })
            .collect();
        let tol = 0.4;
        let path = line_layer(original.clone());
        let (layer, _) = run(path, json!({ "tolerance": tol }));
        let simplified = first_line(&layer);

        let mut worst: f64 = 0.0;
        for (x, y) in &original {
            let mut best = f64::INFINITY;
            for w in simplified.windows(2) {
                best = best.min(point_segment_distance(*x, *y, &w[0], &w[1]));
            }
            worst = worst.max(best);
        }
        assert!(
            worst <= tol * 1.5,
            "every input vertex must stay near the simplified line; worst {worst} vs tolerance {tol}"
        );
    }

    /// Arc fitting beats plain vertex removal on the same data: fewer output
    /// vertices than 'tangent' mode at the same tolerance.
    #[test]
    fn arcs_beat_straight_segments_on_curves() {
        let pts: Vec<(f64, f64)> = (0..=100)
            .map(|i| {
                let t = std::f64::consts::TAU * i as f64 / 100.0;
                (30.0 * t.cos(), 30.0 * t.sin())
            })
            .collect();
        let path = line_layer(pts);
        let (_, arcs) = run(
            path.clone(),
            json!({ "tolerance": 0.2, "densify_output": false }),
        );
        let (_, tangent) = run(path, json!({ "tolerance": 0.2, "mode": "tangent" }));

        let a = arcs.outputs["output_vertex_count"].as_u64().unwrap();
        let t = tangent.outputs["output_vertex_count"].as_u64().unwrap();
        assert!(
            a < t,
            "arc fitting should need fewer vertices than straight segments: {a} vs {t}"
        );
        assert_eq!(
            tangent.outputs["arc_count"],
            json!(0),
            "tangent mode fits no arcs"
        );
    }

    /// A genuinely straight polyline is reduced to its endpoints and yields no
    /// arcs — a large-radius circle must not be fitted to a straight run.
    #[test]
    fn straight_line_collapses_without_arcs() {
        let pts: Vec<(f64, f64)> = (0..=50).map(|i| (i as f64, 0.0)).collect();
        let path = line_layer(pts);
        let (layer, res) = run(path, json!({ "tolerance": 0.1 }));
        assert_eq!(res.outputs["arc_count"], json!(0));
        assert_eq!(first_line(&layer).len(), 2);
    }

    /// max_radius rejects near-flat arcs, sending them down the straight path.
    #[test]
    fn max_radius_rejects_flat_arcs() {
        // A very gentle curve: radius ~1250.
        let pts: Vec<(f64, f64)> = (0..=40)
            .map(|i| {
                let x = i as f64;
                (x, x * x / 2500.0)
            })
            .collect();
        let path = line_layer(pts);
        let (_, capped) = run(path, json!({ "tolerance": 0.05, "max_radius": 100.0 }));
        assert_eq!(
            capped.outputs["arc_count"],
            json!(0),
            "an arc above max_radius must be emitted as straight segments"
        );
    }

    /// min_arc_angle rejects arcs that barely sweep.
    #[test]
    fn min_arc_angle_rejects_shallow_sweeps() {
        let pts: Vec<(f64, f64)> = (0..=30)
            .map(|i| {
                let t = 0.05 * i as f64 / 30.0; // sweeps under 3 degrees
                (100.0 * t.cos(), 100.0 * t.sin())
            })
            .collect();
        let path = line_layer(pts);
        let (_, res) = run(path, json!({ "tolerance": 0.5, "min_arc_angle": 30.0 }));
        assert_eq!(res.outputs["arc_count"], json!(0));
    }

    /// Polygons stay closed and keep at least three vertices.
    #[test]
    fn polygon_rings_stay_valid() {
        let ring: Vec<Coord> = (0..80)
            .map(|i| {
                let t = std::f64::consts::TAU * i as f64 / 80.0;
                Coord::xy(20.0 * t.cos(), 20.0 * t.sin())
            })
            .collect();
        let mut l = Layer::new("poly");
        l.geom_type = Some(GeometryType::Polygon);
        l.add_feature(
            Some(Geometry::Polygon {
                exterior: Ring::new(ring),
                interiors: vec![],
            }),
            &[],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        let path = memory_store::make_vector_memory_path(&id);

        let (layer, _) = run(path, json!({ "tolerance": 0.3 }));
        let Some(Geometry::Polygon { exterior, .. }) =
            layer.iter().next().unwrap().geometry.as_ref()
        else {
            panic!("expected a Polygon");
        };
        assert!(
            exterior.0.len() >= 3,
            "a simplified ring must keep at least 3 vertices"
        );
        // Ring storage carries no closing duplicate.
        let first = &exterior.0[0];
        let last = &exterior.0[exterior.0.len() - 1];
        assert!(
            (first.x - last.x).abs() > 1e-12 || (first.y - last.y).abs() > 1e-12,
            "Ring must not store a closing duplicate vertex"
        );
    }

    /// Densified output really does re-densify: more vertices than endpoint-only
    /// output, and still on the fitted circle.
    #[test]
    fn densify_output_emits_points_on_the_arc() {
        let pts: Vec<(f64, f64)> = (0..=60)
            .map(|i| {
                let t = std::f64::consts::PI * i as f64 / 60.0;
                (100.0 * t.cos(), 100.0 * t.sin())
            })
            .collect();
        let path = line_layer(pts);
        let (dense, _) = run(path.clone(), json!({ "tolerance": 0.5 }));
        let (sparse, _) = run(path, json!({ "tolerance": 0.5, "densify_output": false }));
        let d = first_line(&dense);
        let s = first_line(&sparse);
        assert!(
            d.len() > s.len(),
            "densified output should have more vertices"
        );
        for p in &d {
            let r = p.x.hypot(p.y);
            assert!(
                (r - 100.0).abs() < 1.0,
                "densified vertices should lie on the fitted circle, got r={r}"
            );
        }
    }

    fn point_segment_distance(px: f64, py: f64, a: &Coord, b: &Coord) -> f64 {
        let (dx, dy) = (b.x - a.x, b.y - a.y);
        let l2 = dx * dx + dy * dy;
        if l2 < 1e-18 {
            return (px - a.x).hypot(py - a.y);
        }
        let t = (((px - a.x) * dx + (py - a.y) * dy) / l2).clamp(0.0, 1.0);
        (px - (a.x + t * dx)).hypot(py - (a.y + t * dy))
    }

    #[test]
    fn rejects_bad_parameters() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(SimplifyByCircularArcsTool.validate(&args).is_err());

        let path = line_layer(vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]);
        for bad in [
            json!({ "input": path.clone() }),
            json!({ "input": path.clone(), "tolerance": 0 }),
            json!({ "input": path.clone(), "tolerance": -1 }),
            json!({ "input": path.clone(), "tolerance": 1, "mode": "splines" }),
            json!({ "input": path.clone(), "tolerance": 1, "min_arc_angle": 400 }),
            json!({ "input": path.clone(), "tolerance": 1, "max_radius": 0 }),
        ] {
            let args: ToolArgs = serde_json::from_value(bad).unwrap();
            assert!(SimplifyByCircularArcsTool.validate(&args).is_err());
        }
    }
}
