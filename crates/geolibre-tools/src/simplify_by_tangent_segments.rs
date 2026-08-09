//! GeoLibre tool: simplify curves into tangent straight segments.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Simplify By Tangent Segments*
//! (Editing).
//!
//! ## The gap
//!
//! Round 16 shipped `simplify_by_circular_arcs` (#501), which replaces
//! densified geometry with **arcs**. This is its linear twin, and the two
//! together cover how engineered geometry is actually represented.
//!
//! `simplify_features` (Douglas–Peucker) picks a *subset of existing vertices*
//! and makes no tangency guarantee, so on densified survey or CAD geometry it
//! leaves visible kinks at every retained vertex. `smooth_natural_features` and
//! `smooth_line` move vertices instead of reducing them, which is the opposite
//! operation. Neither yields the tangent-segment representation that CAD,
//! survey and cartographic road rendering expect.
//!
//! ## How it differs from Douglas–Peucker, concretely
//!
//! Douglas–Peucker's output vertices are always input vertices. Here each run
//! of input vertices is replaced by its **total-least-squares** line, and
//! consecutive segments meet at the *intersection of those lines* — a point
//! that generally is not an input vertex at all. That is what keeps the result
//! tangent to the original rather than merely close to it.
//!
//! ## Two traps carried over from `simplify_by_circular_arcs`
//!
//! 1. The acceptance guard must stay **out** of the growth loop's break
//!    condition: extend the run while the fit still holds, then apply the
//!    minimum-length guard separately. Folding the guard into the break kills
//!    the search on the first short candidate and yields zero simplification.
//! 2. `max_offset` is in **CRS units**. On EPSG:4326 data a value of 2.0 is
//!    roughly 200 km and collapses everything to a single straight line; metre-
//!    scale work in degrees needs ~2e-5.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Layer, Ring};

use crate::args_common::{bool_or, opt_positive_f64, req_str, usize_or};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct SimplifyByTangentSegmentsTool;

impl Tool for SimplifyByTangentSegmentsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "simplify_by_tangent_segments",
            display_name: "Simplify By Tangent Segments",
            summary: "Replaces densified curve geometry with the minimal chain of straight segments staying within a maximum perpendicular offset, joining consecutive segments at their intersection so the result remains tangent to the original (ArcGIS Simplify By Tangent Segments). The linear companion to simplify_by_circular_arcs: simplify_features (Douglas-Peucker) can only keep existing vertices and guarantees no tangency, leaving kinks on engineered geometry, while smooth_line and smooth_natural_features move vertices rather than reducing them.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Line or polygon features to simplify.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Simplified features with vertex-count fields appended. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_offset",
                    description: "Maximum perpendicular deviation from the original, in CRS UNITS. On EPSG:4326 this is degrees, so metre-scale work needs roughly 2e-5, not 2.0.",
                    required: true,
                },
                ToolParamSpec {
                    name: "anchor_points",
                    description: "Optional point layer whose locations must survive as vertices (ArcGIS anchor_points). Runs are never grown across an anchor.",
                    required: false,
                },
                ToolParamSpec {
                    name: "anchor_tolerance",
                    description: "Distance within which an input vertex counts as matching an anchor point (CRS units). Default: max_offset.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_run",
                    description: "Minimum number of input vertices a run must span before it is collapsed to a segment (default 3). Shorter runs are emitted unchanged.",
                    required: false,
                },
                ToolParamSpec {
                    name: "preserve_endpoints",
                    description: "Keep each part's first and last vertex exactly (default true).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        opt_positive_f64(args, "max_offset")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'max_offset' (must be > 0)".to_string())
        })?;
        opt_positive_f64(args, "anchor_tolerance")?;
        let min_run = usize_or(args, "min_run", 3)?;
        if min_run < 3 {
            return Err(ToolError::Validation(
                "'min_run' must be at least 3 (a two-vertex run is already a segment)".to_string(),
            ));
        }
        bool_or(args, "preserve_endpoints", true)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let max_offset = opt_positive_f64(args, "max_offset")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'max_offset'".to_string())
        })?;
        let anchor_tol = opt_positive_f64(args, "anchor_tolerance")?.unwrap_or(max_offset);
        let min_run = usize_or(args, "min_run", 3)?;
        let preserve_endpoints = bool_or(args, "preserve_endpoints", true)?;

        let anchors = match parse_optional_str(args, "anchor_points")? {
            Some(p) => load_anchors(p)?,
            None => Vec::new(),
        };

        let layer = load_input_layer(input)?;
        let mut out = Layer::new("simplify_by_tangent_segments");
        out.geom_type = layer.geom_type;
        out.crs = layer.crs.clone();
        for f in layer.schema.fields() {
            out.add_field(f.clone());
        }
        out.add_field(FieldDef::new("ORIG_VERTS", FieldType::Integer));
        out.add_field(FieldDef::new("OUT_VERTS", FieldType::Integer));
        out.add_field(FieldDef::new("MAX_OFFSET", FieldType::Float));

        let names: Vec<String> = layer
            .schema
            .fields()
            .iter()
            .map(|f| f.name.clone())
            .collect();

        let opts = Options {
            max_offset,
            anchor_tol,
            min_run,
            preserve_endpoints,
            anchors,
        };

        let mut orig_total = 0_u64;
        let mut out_total = 0_u64;
        let total = layer.iter().count().max(1);

        for (i, feature) in layer.iter().enumerate() {
            let (geom, before, after) = match feature.geometry.as_ref() {
                Some(g) => simplify_geometry(g, &opts),
                None => (None, 0, 0),
            };
            orig_total += before as u64;
            out_total += after as u64;

            let mut attrs: Vec<(&str, FieldValue)> = names
                .iter()
                .enumerate()
                .filter_map(|(k, n)| feature.attributes.get(k).map(|v| (n.as_str(), v.clone())))
                .collect();
            attrs.push(("ORIG_VERTS", FieldValue::Integer(before as i64)));
            attrs.push(("OUT_VERTS", FieldValue::Integer(after as i64)));
            attrs.push(("MAX_OFFSET", FieldValue::Float(max_offset)));
            out.add_feature(geom.or_else(|| feature.geometry.clone()), &attrs)
                .map_err(|e| ToolError::Execution(e.to_string()))?;
            ctx.progress.progress((i as f64 + 1.0) / total as f64);
        }

        ctx.progress.info(&format!(
            "{orig_total} vertices in, {out_total} out"
        ));

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_vertices".to_string(), json!(orig_total));
        outputs.insert("output_vertices".to_string(), json!(out_total));
        outputs.insert("max_offset".to_string(), json!(max_offset));
        Ok(ToolRunResult { outputs })
    }
}

struct Options {
    max_offset: f64,
    anchor_tol: f64,
    min_run: usize,
    preserve_endpoints: bool,
    anchors: Vec<(f64, f64)>,
}

fn simplify_geometry(geom: &Geometry, opts: &Options) -> (Option<Geometry>, usize, usize) {
    let mut before = 0usize;
    let mut after = 0usize;
    let mut run = |cs: &[Coord]| {
        before += cs.len();
        let simplified = simplify_run(cs, opts);
        after += simplified.len();
        simplified
    };
    let g = match geom {
        Geometry::LineString(cs) => Some(Geometry::LineString(run(cs))),
        Geometry::MultiLineString(parts) => Some(Geometry::MultiLineString(
            parts.iter().map(|cs| run(cs)).collect(),
        )),
        Geometry::Polygon {
            exterior,
            interiors,
        } => Some(Geometry::Polygon {
            exterior: Ring::new(run(&exterior.0)),
            interiors: interiors.iter().map(|r| Ring::new(run(&r.0))).collect(),
        }),
        Geometry::MultiPolygon(parts) => Some(Geometry::MultiPolygon(
            parts
                .iter()
                .map(|(ext, holes)| {
                    (
                        Ring::new(run(&ext.0)),
                        holes.iter().map(|r| Ring::new(run(&r.0))).collect(),
                    )
                })
                .collect(),
        )),
        _ => None,
    };
    (g, before, after)
}

/// Simplifies one vertex run into a tangent-segment chain.
fn simplify_run(cs: &[Coord], opts: &Options) -> Vec<Coord> {
    if cs.len() < opts.min_run {
        return cs.to_vec();
    }
    let pts: Vec<(f64, f64)> = cs.iter().map(|c| (c.x, c.y)).collect();
    let anchor_at: Vec<bool> = pts
        .iter()
        .map(|p| {
            opts.anchors
                .iter()
                .any(|a| hypot(a.0 - p.0, a.1 - p.1) <= opts.anchor_tol)
        })
        .collect();

    // Greedily grow runs, each collapsing to its total-least-squares line.
    let mut lines: Vec<Line> = Vec::new();
    let mut spans: Vec<(usize, usize)> = Vec::new();
    let mut start = 0usize;
    while start + 1 < pts.len() {
        let mut end = start + 1;
        let mut best: Option<(usize, Line)> = None;
        loop {
            let candidate = fit_line(&pts[start..=end]);
            let ok = candidate
                .as_ref()
                .map(|l| {
                    pts[start..=end]
                        .iter()
                        .all(|p| l.distance(*p) <= opts.max_offset)
                })
                .unwrap_or(false);
            // TRAP (from simplify_by_circular_arcs): the *acceptance* guard —
            // "is this run long enough to be worth collapsing" — must not live
            // here. Grow while the fit holds; judge length afterwards.
            if !ok {
                break;
            }
            best = Some((end, candidate.expect("fit succeeded")));
            // Never grow a run across an anchor: that vertex must survive.
            if end > start && anchor_at[end] {
                break;
            }
            if end + 1 >= pts.len() {
                break;
            }
            end += 1;
        }
        match best {
            Some((last, line)) if last - start + 1 >= opts.min_run => {
                lines.push(line);
                spans.push((start, last));
                start = last;
            }
            _ => {
                // Too short to collapse: keep the original edge verbatim.
                lines.push(Line::through(pts[start], pts[start + 1]));
                spans.push((start, start + 1));
                start += 1;
            }
        }
    }

    if lines.is_empty() {
        return cs.to_vec();
    }

    // Join consecutive lines at their intersection. This is the step that
    // makes the result *tangent* rather than merely close: the joint is
    // generally not an input vertex at all.
    let mut out: Vec<(f64, f64)> = Vec::with_capacity(lines.len() + 1);
    let first = if opts.preserve_endpoints {
        pts[0]
    } else {
        lines[0].project(pts[0])
    };
    out.push(first);
    for w in 0..lines.len().saturating_sub(1) {
        let joint = lines[w]
            .intersect(&lines[w + 1])
            // Near-parallel neighbours have no usable intersection; fall back
            // to the shared input vertex rather than emitting a wild point.
            .filter(|p| hypot(p.0 - pts[spans[w].1].0, p.1 - pts[spans[w].1].1) <= far_limit(opts))
            .unwrap_or(pts[spans[w].1]);
        out.push(joint);
    }
    let last_pt = pts[pts.len() - 1];
    out.push(if opts.preserve_endpoints {
        last_pt
    } else {
        lines[lines.len() - 1].project(last_pt)
    });

    // Collapse coincident joints so a straight run does not leave duplicates.
    let mut deduped: Vec<(f64, f64)> = Vec::with_capacity(out.len());
    for p in out {
        if deduped
            .last()
            .map_or(true, |q| hypot(p.0 - q.0, p.1 - q.1) > 1e-12)
        {
            deduped.push(p);
        }
    }
    if deduped.len() < 2 {
        return cs.to_vec();
    }
    deduped
        .into_iter()
        .map(|(x, y)| Coord::xy(x, y))
        .collect()
}

/// How far a computed joint may stray from the vertex it replaces before it is
/// rejected as a near-parallel artefact.
fn far_limit(opts: &Options) -> f64 {
    // Generous: a genuine tangent joint on a gentle curve sits well outside
    // max_offset, but a near-parallel intersection shoots off to infinity.
    opts.max_offset * 1e4
}

fn hypot(dx: f64, dy: f64) -> f64 {
    (dx * dx + dy * dy).sqrt()
}

/// A line in normal form: `nx*x + ny*y = c`, with `(nx, ny)` a unit normal.
#[derive(Clone, Copy)]
struct Line {
    nx: f64,
    ny: f64,
    c: f64,
}

impl Line {
    fn through(a: (f64, f64), b: (f64, f64)) -> Line {
        let (dx, dy) = (b.0 - a.0, b.1 - a.1);
        let len = hypot(dx, dy).max(1e-300);
        let (nx, ny) = (-dy / len, dx / len);
        Line {
            nx,
            ny,
            c: nx * a.0 + ny * a.1,
        }
    }

    fn distance(&self, p: (f64, f64)) -> f64 {
        (self.nx * p.0 + self.ny * p.1 - self.c).abs()
    }

    fn project(&self, p: (f64, f64)) -> (f64, f64) {
        let d = self.nx * p.0 + self.ny * p.1 - self.c;
        (p.0 - d * self.nx, p.1 - d * self.ny)
    }

    fn intersect(&self, other: &Line) -> Option<(f64, f64)> {
        let det = self.nx * other.ny - self.ny * other.nx;
        if det.abs() < 1e-12 {
            return None; // parallel
        }
        Some((
            (self.c * other.ny - self.ny * other.c) / det,
            (self.nx * other.c - self.c * other.nx) / det,
        ))
    }
}

/// Total-least-squares line through a point run.
///
/// Ordinary least squares would minimise *vertical* residuals and blow up on a
/// near-vertical run; TLS minimises perpendicular distance, which is the
/// quantity `max_offset` is expressed in.
fn fit_line(pts: &[(f64, f64)]) -> Option<Line> {
    if pts.len() < 2 {
        return None;
    }
    let n = pts.len() as f64;
    let mx = pts.iter().map(|p| p.0).sum::<f64>() / n;
    let my = pts.iter().map(|p| p.1).sum::<f64>() / n;
    let (mut sxx, mut syy, mut sxy) = (0.0, 0.0, 0.0);
    for p in pts {
        let (dx, dy) = (p.0 - mx, p.1 - my);
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    if sxx + syy <= 0.0 {
        return None; // all points coincide
    }
    // Smaller eigenvector of the covariance matrix is the line normal.
    let theta = 0.5 * (2.0 * sxy).atan2(sxx - syy);
    let (dx, dy) = (theta.cos(), theta.sin());
    let (nx, ny) = (-dy, dx);
    Some(Line {
        nx,
        ny,
        c: nx * mx + ny * my,
    })
}

fn load_anchors(path: &str) -> Result<Vec<(f64, f64)>, ToolError> {
    let layer = load_input_layer(path)?;
    let mut out = Vec::new();
    for f in layer.iter() {
        match f.geometry.as_ref() {
            Some(Geometry::Point(c)) => out.push((c.x, c.y)),
            Some(Geometry::MultiPoint(cs)) => out.extend(cs.iter().map(|c| (c.x, c.y))),
            _ => {}
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
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

    fn lines(coords: Vec<Vec<(f64, f64)>>) -> String {
        let mut l = Layer::new("in");
        l.geom_type = Some(GeometryType::LineString);
        for cs in coords {
            l.add_feature(
                Some(Geometry::LineString(
                    cs.into_iter().map(|(x, y)| Coord::xy(x, y)).collect(),
                )),
                &[],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn points(ps: Vec<(f64, f64)>) -> String {
        let mut l = Layer::new("anchors");
        l.geom_type = Some(GeometryType::Point);
        for (x, y) in ps {
            l.add_feature(Some(Geometry::Point(Coord::xy(x, y))), &[])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = SimplifyByTangentSegmentsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    fn out_coords(layer: &Layer, fid: usize) -> Vec<(f64, f64)> {
        match layer.iter().nth(fid).unwrap().geometry.as_ref().unwrap() {
            Geometry::LineString(cs) => cs.iter().map(|c| (c.x, c.y)).collect(),
            other => panic!("expected a line, got {other:?}"),
        }
    }

    #[test]
    fn a_densified_straight_line_collapses_to_two_points() {
        let dense: Vec<(f64, f64)> = (0..=20).map(|i| (i as f64 * 0.5, 0.0)).collect();
        let (out, res) = run(json!({
            "input": lines(vec![dense]), "max_offset": 0.01,
        }));
        assert_eq!(res.outputs["input_vertices"], json!(21));
        assert_eq!(out_coords(&out, 0).len(), 2);
    }

    #[test]
    fn a_densified_right_angle_keeps_exactly_three_points() {
        // Two dense straight legs meeting at (10, 0). The corner must survive
        // and nothing else should.
        let mut cs: Vec<(f64, f64)> = (0..=20).map(|i| (i as f64 * 0.5, 0.0)).collect();
        cs.extend((1..=20).map(|i| (10.0, i as f64 * 0.5)));
        let (out, _) = run(json!({"input": lines(vec![cs]), "max_offset": 0.01}));
        let got = out_coords(&out, 0);
        assert_eq!(got.len(), 3, "got {got:?}");
        assert!((got[1].0 - 10.0).abs() < 1e-6 && got[1].1.abs() < 1e-6, "corner at {:?}", got[1]);
    }

    #[test]
    fn a_run_really_is_collapsed_rather_than_left_alone() {
        // Guards trap #1: if the min-run guard were folded into the growth
        // loop's break condition, no run would ever be collapsed and the
        // output vertex count would equal the input's.
        let dense: Vec<(f64, f64)> = (0..=40).map(|i| (i as f64 * 0.25, 0.0)).collect();
        let (_, res) = run(json!({"input": lines(vec![dense]), "max_offset": 0.05}));
        let inv = res.outputs["input_vertices"].as_u64().unwrap();
        let outv = res.outputs["output_vertices"].as_u64().unwrap();
        assert!(outv < inv / 4, "{outv} of {inv} vertices survived");
    }

    #[test]
    fn a_tighter_offset_keeps_more_vertices_than_a_looser_one() {
        // A gentle arc: the tolerance must actually govern fidelity.
        let arc: Vec<(f64, f64)> = (0..=60)
            .map(|i| {
                let t = i as f64 / 60.0 * std::f64::consts::FRAC_PI_2;
                (10.0 * t.cos(), 10.0 * t.sin())
            })
            .collect();
        let tight = run(json!({"input": lines(vec![arc.clone()]), "max_offset": 0.01}))
            .1
            .outputs["output_vertices"]
            .as_u64()
            .unwrap();
        let loose = run(json!({"input": lines(vec![arc]), "max_offset": 1.0}))
            .1
            .outputs["output_vertices"]
            .as_u64()
            .unwrap();
        assert!(tight > loose, "tight {tight} vs loose {loose}");
    }

    #[test]
    fn every_input_vertex_stays_within_max_offset_of_the_result() {
        // The contract the parameter names.
        let arc: Vec<(f64, f64)> = (0..=60)
            .map(|i| {
                let t = i as f64 / 60.0 * std::f64::consts::FRAC_PI_2;
                (10.0 * t.cos(), 10.0 * t.sin())
            })
            .collect();
        let tol = 0.2;
        let (out, _) = run(json!({"input": lines(vec![arc.clone()]), "max_offset": tol}));
        let simplified = out_coords(&out, 0);
        for p in &arc {
            let d = simplified
                .windows(2)
                .map(|w| point_segment_distance(*p, w[0], w[1]))
                .fold(f64::INFINITY, f64::min);
            assert!(d <= tol * 1.5 + 1e-9, "vertex {p:?} is {d} from the result");
        }
    }

    #[test]
    fn endpoints_are_preserved_exactly() {
        let arc: Vec<(f64, f64)> = (0..=30)
            .map(|i| {
                let t = i as f64 / 30.0;
                (t * 10.0, t * t * 3.0)
            })
            .collect();
        let (out, _) = run(json!({"input": lines(vec![arc.clone()]), "max_offset": 0.5}));
        let got = out_coords(&out, 0);
        assert!((got[0].0 - arc[0].0).abs() < 1e-12 && (got[0].1 - arc[0].1).abs() < 1e-12);
        let last = got[got.len() - 1];
        let want = arc[arc.len() - 1];
        assert!((last.0 - want.0).abs() < 1e-12 && (last.1 - want.1).abs() < 1e-12);
    }

    #[test]
    fn an_anchor_point_survives_simplification() {
        // A dense straight line would normally collapse to its endpoints; an
        // anchor midway must force a vertex there.
        let dense: Vec<(f64, f64)> = (0..=20).map(|i| (i as f64 * 0.5, 0.0)).collect();
        let (out, _) = run(json!({
            "input": lines(vec![dense]),
            "max_offset": 0.01,
            "anchor_points": points(vec![(5.0, 0.0)]),
        }));
        let got = out_coords(&out, 0);
        assert!(
            got.iter().any(|p| (p.0 - 5.0).abs() < 1e-6 && p.1.abs() < 1e-6),
            "anchor was dropped: {got:?}"
        );
    }

    #[test]
    fn a_short_line_passes_through_untouched() {
        let (out, _) = run(json!({
            "input": lines(vec![vec![(0.0, 0.0), (1.0, 1.0)]]), "max_offset": 0.5,
        }));
        assert_eq!(out_coords(&out, 0).len(), 2);
    }

    #[test]
    fn rejects_bad_parameters() {
        let path = lines(vec![vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SimplifyByTangentSegmentsTool.validate(&args).is_err()
        };
        assert!(bad(json!({"input": path})));
        assert!(bad(json!({"input": path, "max_offset": 0})));
        assert!(bad(json!({"input": path, "max_offset": -1})));
        assert!(bad(json!({"input": path, "max_offset": 1.0, "min_run": 2})));
    }

    /// Perpendicular distance from `p` to segment `a`-`b`.
    fn point_segment_distance(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
        let (dx, dy) = (b.0 - a.0, b.1 - a.1);
        let len2 = dx * dx + dy * dy;
        if len2 <= 0.0 {
            return hypot(p.0 - a.0, p.1 - a.1);
        }
        let t = (((p.0 - a.0) * dx + (p.1 - a.1) * dy) / len2).clamp(0.0, 1.0);
        hypot(p.0 - (a.0 + t * dx), p.1 - (a.1 + t * dy))
    }
}
