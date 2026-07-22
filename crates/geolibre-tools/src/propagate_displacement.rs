//! GeoLibre tool: propagate displacement-cartography shifts onto nearby features.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Propagate Displacement* (Cartography).
//! `resolve_road_conflicts` moves roads to clear symbol-width conflicts and can
//! emit its per-vertex shift as a `links` layer; `resolve_building_conflicts`
//! displaces buildings that collide with barriers. Neither one *propagates* that
//! shift onto features that were left behind — a building beside a road that just
//! moved, a hydro line, a parcel edge — which shears the map. This tool closes
//! that gap: it interpolates the sparse `links` vectors into a smooth
//! displacement field with inverse-distance weighting (the same IDW core behind
//! `rubbersheet_features`, without a border) and applies it to `input`.
//!
//! `links` is a two-vertex line layer where each line's first vertex is the
//! original location and the last vertex is the displaced location — exactly the
//! schema `resolve_road_conflicts`'s optional `links` output writes, and the same
//! `links` layer `rubbersheet_features` already consumes.
//!
//! `adjustment_style` controls how a feature absorbs the sampled field:
//! * `solid` — rigid: translate the whole feature by the field sampled once at
//!   its representative point.
//! * `preserve_orientation` — rigid: translate the whole feature by the field
//!   *averaged* over every vertex of its footprint (steadier for large or
//!   elongated shapes than a single sample), still a pure translation so right
//!   angles and orientation are exactly preserved.
//! * `auto` (default) — rigid features (points, polygons — buildings, parcels)
//!   use the `preserve_orientation` averaged translation; flexible line features
//!   (hydro, boundaries) are bent: each vertex is displaced by the field sampled
//!   there, then lightly smoothed along arc length so the line doesn't kink where
//!   the field varies quickly between vertices.
//!
//! `search_distance` caps a link's influence: links farther than this from a
//! sample point are excluded, so far-away features are left in place instead of
//! being dragged by a distant, unrelated shift (mirrors "nothing to reference"
//! outside a link's zone of relevance rather than falling back to a global mean).
//!
//! v1 scope: `GeometryCollection` features are handled as a single rigid unit
//! (not split into rigid/flexible parts) since ArcGIS features are rarely mixed
//! geometry collections in practice.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Geometry, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// IDW power for the displacement-field interpolation (matches
/// `rubbersheet_features`'s default).
const IDW_POWER: f64 = 2.0;

pub struct PropagateDisplacementTool;

impl Tool for PropagateDisplacementTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "propagate_displacement",
            display_name: "Propagate Displacement",
            summary: "Propagate the displacement introduced by resolve_road_conflicts (or any two-point displacement-links layer) onto nearby features by interpolating the sparse links into a smooth IDW field: rigid features (buildings, parcels) get a single averaged translation, flexible lines get a per-vertex, arc-length-smoothed bend, like ArcGIS Propagate Displacement.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer to adjust (any geometry type).",
                    required: true,
                },
                ToolParamSpec {
                    name: "links",
                    description: "Displacement links: a 2-vertex line layer per link, first vertex = original location, last vertex = displaced location (the schema resolve_road_conflicts's 'links' output and rubbersheet_features's 'links' input both use).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "adjustment_style",
                    description: "'auto' (default: rigid translate for points/polygons, per-vertex bend for lines), 'preserve_orientation' (rigid translate by the field averaged over the footprint), or 'solid' (rigid translate by the field sampled once at the feature's representative point).",
                    required: false,
                },
                ToolParamSpec {
                    name: "search_distance",
                    description: "Maximum distance (CRS units) a link may influence a sample point; farther links are ignored. Default: unlimited (every link contributes everywhere).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "links")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let links_path = require_str(args, "links")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let link_layer = load_input_layer(links_path)?;
        let links = links_from_layer(&link_layer);
        if links.is_empty() {
            return Err(ToolError::Execution(
                "no displacement links (nothing to propagate)".to_string(),
            ));
        }

        ctx.progress.info(&format!(
            "{} displacement link(s); adjustment_style={}",
            links.len(),
            prm.style.name()
        ));

        let field = DisplacementField {
            links,
            power: IDW_POWER,
            search_distance: prm.search_distance,
        };

        let mut rigid_count = 0usize;
        let mut flexible_count = 0usize;
        let mut total_disp = 0.0f64;
        let mut n_touched = 0u64;

        for feature in layer.features.iter_mut() {
            let Some(geom) = feature.geometry.take() else {
                continue;
            };
            let rigid = match prm.style {
                Style::Solid | Style::PreserveOrientation => true,
                Style::Auto => !is_flexible(&geom),
            };
            let new_geom = if rigid {
                rigid_count += 1;
                let (dx, dy) = match prm.style {
                    Style::Solid => {
                        let (cx, cy) = representative_xy(&geom).unwrap_or((0.0, 0.0));
                        field.sample(cx, cy)
                    }
                    Style::PreserveOrientation | Style::Auto => {
                        averaged_displacement(&geom, &field)
                    }
                };
                translate_geometry(&geom, dx, dy, &mut total_disp, &mut n_touched)
            } else {
                flexible_count += 1;
                flexible_displace(&geom, &field, &mut total_disp, &mut n_touched)
            };
            feature.geometry = Some(new_geom);
        }
        layer.extent = None;

        let mean_disp = if n_touched > 0 {
            total_disp / n_touched as f64
        } else {
            0.0
        };
        ctx.progress.info(&format!(
            "{rigid_count} rigid, {flexible_count} flexible feature(s); mean displacement {mean_disp:.4}"
        ));

        let feature_count = layer.len();
        let link_count = field.links.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("link_count".to_string(), json!(link_count));
        outputs.insert("rigid_count".to_string(), json!(rigid_count));
        outputs.insert("flexible_count".to_string(), json!(flexible_count));
        outputs.insert("mean_displacement".to_string(), json!(mean_disp));
        Ok(ToolRunResult { outputs })
    }
}

// ── Links ────────────────────────────────────────────────────────────────────

/// A displacement link: a source point and its displacement vector.
#[derive(Clone, Copy)]
struct Link {
    fx: f64,
    fy: f64,
    dx: f64,
    dy: f64,
}

/// Read links from a line layer: each line's first vertex is the original
/// location, last vertex the displaced location. Matches
/// `resolve_road_conflicts`'s `links` output and `rubbersheet_features`'s
/// `links` input.
fn links_from_layer(layer: &Layer) -> Vec<Link> {
    let mut out = Vec::new();
    for feature in layer.features.iter() {
        match feature.geometry.as_ref() {
            Some(Geometry::LineString(cs)) if cs.len() >= 2 => {
                let a = &cs[0];
                let b = &cs[cs.len() - 1];
                out.push(Link {
                    fx: a.x,
                    fy: a.y,
                    dx: b.x - a.x,
                    dy: b.y - a.y,
                });
            }
            Some(Geometry::MultiLineString(lines)) => {
                for cs in lines {
                    if cs.len() >= 2 {
                        let a = &cs[0];
                        let b = &cs[cs.len() - 1];
                        out.push(Link {
                            fx: a.x,
                            fy: a.y,
                            dx: b.x - a.x,
                            dy: b.y - a.y,
                        });
                    }
                }
            }
            _ => {}
        }
    }
    out
}

// ── Displacement field (IDW) ────────────────────────────────────────────────

struct DisplacementField {
    links: Vec<Link>,
    power: f64,
    search_distance: Option<f64>,
}

impl DisplacementField {
    /// Displacement (dx, dy) to apply at world (x, y): inverse-distance-weighted
    /// blend of every link within `search_distance` (or all links, unbounded).
    /// No links in range -> no displacement (the feature stays put).
    fn sample(&self, x: f64, y: f64) -> (f64, f64) {
        let mut sx = 0.0;
        let mut sy = 0.0;
        let mut sw = 0.0;
        for l in &self.links {
            let d = (x - l.fx).hypot(y - l.fy);
            if let Some(sd) = self.search_distance {
                if d > sd {
                    continue;
                }
            }
            if d < 1e-9 {
                return (l.dx, l.dy);
            }
            let w = 1.0 / d.powf(self.power);
            sx += w * l.dx;
            sy += w * l.dy;
            sw += w;
        }
        if sw > 0.0 {
            (sx / sw, sy / sw)
        } else {
            (0.0, 0.0)
        }
    }
}

// ── Rigid (translate) handling ──────────────────────────────────────────────

/// Field displacement averaged over every vertex of `geom`'s footprint.
fn averaged_displacement(geom: &Geometry, field: &DisplacementField) -> (f64, f64) {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0u64;
    walk_coords(geom, &mut |x, y| {
        let (dx, dy) = field.sample(x, y);
        sx += dx;
        sy += dy;
        n += 1;
    });
    if n > 0 {
        (sx / n as f64, sy / n as f64)
    } else {
        (0.0, 0.0)
    }
}

/// Translate every vertex of `geom` by the fixed vector `(dx, dy)`, recording
/// each vertex's displacement magnitude into the running mean.
fn translate_geometry(
    geom: &Geometry,
    dx: f64,
    dy: f64,
    total_disp: &mut f64,
    n_touched: &mut u64,
) -> Geometry {
    let mag = (dx * dx + dy * dy).sqrt();
    map_geometry(geom, &mut |x, y| {
        *total_disp += mag;
        *n_touched += 1;
        (x + dx, y + dy)
    })
}

// ── Flexible (per-vertex bend) handling ─────────────────────────────────────

/// True for the geometry types treated as flexible (bent per vertex) under
/// `adjustment_style = auto`; everything else (points, polygons, and — as a
/// simplification — geometry collections) is rigid.
fn is_flexible(geom: &Geometry) -> bool {
    matches!(geom, Geometry::LineString(_) | Geometry::MultiLineString(_))
}

fn flexible_displace(
    geom: &Geometry,
    field: &DisplacementField,
    total_disp: &mut f64,
    n_touched: &mut u64,
) -> Geometry {
    match geom {
        Geometry::LineString(cs) => {
            Geometry::LineString(displace_chain(cs, field, total_disp, n_touched))
        }
        Geometry::MultiLineString(lines) => Geometry::MultiLineString(
            lines
                .iter()
                .map(|cs| displace_chain(cs, field, total_disp, n_touched))
                .collect(),
        ),
        other => other.clone(),
    }
}

/// Displace each vertex of a line chain by the field sampled there, then
/// lightly smooth the displacement along arc length so the line doesn't kink
/// between widely-varying samples. Endpoints are recorded raw (no neighbour to
/// average with); interior vertices blend with their neighbours, inverse-
/// distance weighted by segment length, so unevenly digitized lines don't bias
/// the smoothing toward the densely-vertexed side.
fn displace_chain(
    cs: &[Coord],
    field: &DisplacementField,
    total_disp: &mut f64,
    n_touched: &mut u64,
) -> Vec<Coord> {
    if cs.is_empty() {
        return Vec::new();
    }
    let raw: Vec<(f64, f64)> = cs.iter().map(|c| field.sample(c.x, c.y)).collect();
    let smoothed = smooth_along_arclength(&raw, cs);
    cs.iter()
        .zip(smoothed.iter())
        .map(|(c, &(dx, dy))| {
            *total_disp += (dx * dx + dy * dy).sqrt();
            *n_touched += 1;
            Coord::xy(c.x + dx, c.y + dy)
        })
        .collect()
}

fn smooth_along_arclength(raw: &[(f64, f64)], coords: &[Coord]) -> Vec<(f64, f64)> {
    let n = raw.len();
    if n < 3 {
        return raw.to_vec();
    }
    let mut out = raw.to_vec();
    for k in 1..n - 1 {
        let d_prev = (coords[k].x - coords[k - 1].x)
            .hypot(coords[k].y - coords[k - 1].y)
            .max(1e-9);
        let d_next = (coords[k + 1].x - coords[k].x)
            .hypot(coords[k + 1].y - coords[k].y)
            .max(1e-9);
        let w_prev = 1.0 / d_prev;
        let w_next = 1.0 / d_next;
        let w_self = w_prev + w_next;
        let total = 2.0 * w_self;
        out[k].0 = (w_self * raw[k].0 + w_prev * raw[k - 1].0 + w_next * raw[k + 1].0) / total;
        out[k].1 = (w_self * raw[k].1 + w_prev * raw[k - 1].1 + w_next * raw[k + 1].1) / total;
    }
    out
}

// ── Geometry helpers ─────────────────────────────────────────────────────────

fn map_geometry(geom: &Geometry, f: &mut impl FnMut(f64, f64) -> (f64, f64)) -> Geometry {
    let map_chain = |cs: &[Coord], f: &mut dyn FnMut(f64, f64) -> (f64, f64)| -> Vec<Coord> {
        cs.iter()
            .map(|c| {
                let (x, y) = f(c.x, c.y);
                Coord::xy(x, y)
            })
            .collect()
    };
    match geom {
        Geometry::Point(c) => {
            let (x, y) = f(c.x, c.y);
            Geometry::Point(Coord::xy(x, y))
        }
        Geometry::MultiPoint(cs) => Geometry::MultiPoint(map_chain(cs, f)),
        Geometry::LineString(cs) => Geometry::LineString(map_chain(cs, f)),
        Geometry::MultiLineString(lines) => {
            Geometry::MultiLineString(lines.iter().map(|l| map_chain(l, f)).collect())
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => Geometry::Polygon {
            exterior: Ring::new(map_chain(exterior.coords(), f)),
            interiors: interiors
                .iter()
                .map(|r| Ring::new(map_chain(r.coords(), f)))
                .collect(),
        },
        Geometry::MultiPolygon(parts) => Geometry::MultiPolygon(
            parts
                .iter()
                .map(|(e, holes)| {
                    (
                        Ring::new(map_chain(e.coords(), f)),
                        holes
                            .iter()
                            .map(|r| Ring::new(map_chain(r.coords(), f)))
                            .collect(),
                    )
                })
                .collect(),
        ),
        Geometry::GeometryCollection(gs) => {
            Geometry::GeometryCollection(gs.iter().map(|g| map_geometry(g, f)).collect())
        }
    }
}

fn walk_coords(geom: &Geometry, add: &mut impl FnMut(f64, f64)) {
    match geom {
        Geometry::Point(c) => add(c.x, c.y),
        Geometry::LineString(cs) | Geometry::MultiPoint(cs) => {
            cs.iter().for_each(|c| add(c.x, c.y))
        }
        Geometry::MultiLineString(ls) => ls.iter().flatten().for_each(|c| add(c.x, c.y)),
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            exterior.coords().iter().for_each(|c| add(c.x, c.y));
            interiors
                .iter()
                .for_each(|r| r.coords().iter().for_each(|c| add(c.x, c.y)));
        }
        Geometry::MultiPolygon(ps) => {
            for (e, h) in ps {
                e.coords().iter().for_each(|c| add(c.x, c.y));
                h.iter()
                    .for_each(|r| r.coords().iter().for_each(|c| add(c.x, c.y)));
            }
        }
        Geometry::GeometryCollection(gs) => gs.iter().for_each(|g| walk_coords(g, add)),
    }
}

fn representative_xy(geom: &Geometry) -> Option<(f64, f64)> {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0u64;
    walk_coords(geom, &mut |x, y| {
        sx += x;
        sy += y;
        n += 1;
    });
    (n > 0).then(|| (sx / n as f64, sy / n as f64))
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Style {
    Auto,
    PreserveOrientation,
    Solid,
}

impl Style {
    fn name(self) -> &'static str {
        match self {
            Style::Auto => "auto",
            Style::PreserveOrientation => "preserve_orientation",
            Style::Solid => "solid",
        }
    }
}

struct Params {
    style: Style,
    search_distance: Option<f64>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let style = match parse_optional_str(args, "adjustment_style")? {
        None => Style::Auto,
        Some(s) => match s.trim().to_ascii_lowercase().replace(['-', ' '], "_").as_str() {
            "" | "auto" => Style::Auto,
            "preserve_orientation" => Style::PreserveOrientation,
            "solid" => Style::Solid,
            other => {
                return Err(ToolError::Validation(format!(
                    "'adjustment_style' must be 'auto', 'preserve_orientation', or 'solid', got '{other}'"
                )))
            }
        },
    };
    let search_distance = parse_optional_f64(args, "search_distance")?;
    if let Some(sd) = search_distance {
        if !(sd.is_finite() && sd > 0.0) {
            return Err(ToolError::Validation(
                "'search_distance' must be a positive number".to_string(),
            ));
        }
    }
    Ok(Params {
        style,
        search_distance,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
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
    use wbvector::{memory_store, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn square(x0: f64, y0: f64, s: f64) -> Geometry {
        Geometry::polygon(
            vec![
                Coord::xy(x0, y0),
                Coord::xy(x0 + s, y0),
                Coord::xy(x0 + s, y0 + s),
                Coord::xy(x0, y0 + s),
            ],
            vec![],
        )
    }

    fn building_layer(sqs: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("b")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        for &(x, y, s) in sqs {
            l.add_feature(Some(square(x, y, s)), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn line_layer(lines: &[Vec<(f64, f64)>]) -> String {
        let mut l = Layer::new("ln")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        for ln in lines {
            l.add_feature(
                Some(Geometry::LineString(
                    ln.iter().map(|(x, y)| Coord::xy(*x, *y)).collect(),
                )),
                &[],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    type LinkPair = ((f64, f64), (f64, f64));

    fn link_layer(links: &[LinkPair]) -> String {
        let mut l = Layer::new("links")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        for &(a, b) in links {
            l.add_feature(
                Some(Geometry::line_string(vec![
                    Coord::xy(a.0, a.1),
                    Coord::xy(b.0, b.1),
                ])),
                &[],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = PropagateDisplacementTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn poly_pts(g: &Geometry) -> Vec<(f64, f64)> {
        match g {
            Geometry::Polygon { exterior, .. } => {
                exterior.coords().iter().map(|c| (c.x, c.y)).collect()
            }
            other => panic!("expected polygon, got {other:?}"),
        }
    }

    fn line_pts(g: &Geometry) -> Vec<(f64, f64)> {
        match g {
            Geometry::LineString(cs) => cs.iter().map(|c| (c.x, c.y)).collect(),
            other => panic!("expected linestring, got {other:?}"),
        }
    }

    /// A building near a single displacement link translates by ~the link vector
    /// (a single link's IDW field is constant everywhere).
    #[test]
    fn building_near_link_translates_by_link_vector() {
        let b = building_layer(&[(0.0, 0.0, 10.0)]);
        let links = link_layer(&[((5.0, 5.0), (15.0, 12.0))]); // dx=10, dy=7
        let (out, layer) = run(json!({ "input": b, "links": links }));
        assert_eq!(out.outputs["rigid_count"], json!(1));
        let pts = poly_pts(layer.features[0].geometry.as_ref().unwrap());
        // First vertex (0,0) -> should land near (10,7).
        let (x, y) = pts[0];
        assert!(
            (x - 10.0).abs() < 1e-6 && (y - 7.0).abs() < 1e-6,
            "vertex -> ({x},{y}), want (10,7)"
        );
    }

    /// With a search_distance, a building far outside every link's influence is
    /// left unmoved.
    #[test]
    fn far_building_unmoved_with_search_distance() {
        let b = building_layer(&[(10_000.0, 10_000.0, 10.0)]);
        let links = link_layer(&[((0.0, 0.0), (100.0, 100.0))]);
        let (out, layer) = run(json!({
            "input": b, "links": links, "search_distance": 50.0
        }));
        assert_eq!(out.outputs["mean_displacement"], json!(0.0));
        let pts = poly_pts(layer.features[0].geometry.as_ref().unwrap());
        assert_eq!(pts[0], (10_000.0, 10_000.0));
    }

    /// preserve_orientation keeps a rectangle axis-aligned and congruent (pure
    /// translation) even under a non-uniform (2-link) field.
    #[test]
    fn preserve_orientation_keeps_rectangle_congruent() {
        // A 20 x 8 rectangle, with two links giving a non-uniform field so a
        // naive per-vertex displacement would distort it.
        let mut l = Layer::new("r")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_feature(
            Some(Geometry::polygon(
                vec![
                    Coord::xy(0.0, 0.0),
                    Coord::xy(20.0, 0.0),
                    Coord::xy(20.0, 8.0),
                    Coord::xy(0.0, 8.0),
                ],
                vec![],
            )),
            &[],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);
        let links = link_layer(&[((-50.0, 4.0), (-40.0, 4.0)), ((70.0, 4.0), (70.0, 40.0))]);
        let (out, layer) = run(json!({
            "input": input, "links": links, "adjustment_style": "preserve_orientation"
        }));
        assert_eq!(out.outputs["rigid_count"], json!(1));
        let pts = poly_pts(layer.features[0].geometry.as_ref().unwrap());
        // Congruent: same edge lengths as the original 20x8 rectangle, and still
        // axis-aligned (bottom edge 0-1 shares y, right edge 1-2 shares x).
        let width = (pts[1].0 - pts[0].0).abs();
        let height = (pts[2].1 - pts[1].1).abs();
        assert!((width - 20.0).abs() < 1e-6, "width changed: {width}");
        assert!((height - 8.0).abs() < 1e-6, "height changed: {height}");
        assert!((pts[0].1 - pts[1].1).abs() < 1e-9, "not axis-aligned");
        assert!((pts[1].0 - pts[2].0).abs() < 1e-9, "not axis-aligned");
    }

    /// Under auto, a line bends per-vertex (non-uniform shift) rather than
    /// translating rigidly, when the field varies along its length.
    #[test]
    fn auto_bends_line_per_vertex() {
        let ln = line_layer(&[vec![(0.0, 0.0), (50.0, 0.0), (100.0, 0.0)]]);
        // Two links with very different y-displacement at each end of the line.
        let links = link_layer(&[((0.0, 0.0), (0.0, 20.0)), ((100.0, 0.0), (100.0, 0.0))]);
        let (out, layer) = run(json!({ "input": ln, "links": links, "adjustment_style": "auto" }));
        assert_eq!(out.outputs["flexible_count"], json!(1));
        let pts = line_pts(layer.features[0].geometry.as_ref().unwrap());
        let dy0 = pts[0].1;
        let dy1 = pts[1].1;
        let dy2 = pts[2].1;
        // First vertex shifted close to +20, last close to 0, middle in between
        // and distinct from both -> a bend, not a uniform rigid shift.
        assert!(dy0 > 15.0, "start not shifted enough: {dy0}");
        assert!(dy2.abs() < 1.0, "end shifted too much: {dy2}");
        assert!(
            dy1 > dy2 + 1.0 && dy1 < dy0 - 1.0,
            "middle vertex not between start/end: {dy1} (start {dy0}, end {dy2})"
        );
    }

    /// solid and preserve_orientation both force a rigid translation even for a
    /// line feature (unlike auto, which would bend it).
    #[test]
    fn solid_forces_rigid_translation_on_lines() {
        let ln = line_layer(&[vec![(0.0, 0.0), (50.0, 0.0), (100.0, 0.0)]]);
        let links = link_layer(&[((0.0, 0.0), (0.0, 20.0)), ((100.0, 0.0), (100.0, 0.0))]);
        let (out, layer) = run(json!({ "input": ln, "links": links, "adjustment_style": "solid" }));
        assert_eq!(out.outputs["rigid_count"], json!(1));
        let pts = line_pts(layer.features[0].geometry.as_ref().unwrap());
        let dy0 = pts[0].1;
        let dy1 = pts[1].1;
        let dy2 = pts[2].1;
        assert!(
            (dy0 - dy1).abs() < 1e-9 && (dy1 - dy2).abs() < 1e-9,
            "solid should shift every vertex identically: {dy0}, {dy1}, {dy2}"
        );
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            PropagateDisplacementTool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing everything");
        assert!(
            bad(json!({ "input": "a.geojson" })).is_err(),
            "missing links"
        );
        assert!(
            bad(json!({ "links": "l.geojson" })).is_err(),
            "missing input"
        );
        assert!(
            bad(json!({ "input": "a.geojson", "links": "l.geojson", "adjustment_style": "bogus" }))
                .is_err(),
            "bad adjustment_style"
        );
        assert!(
            bad(json!({ "input": "a.geojson", "links": "l.geojson", "search_distance": -5.0 }))
                .is_err(),
            "negative search_distance"
        );
        assert!(
            bad(json!({ "input": "a.geojson", "links": "l.geojson" })).is_ok(),
            "minimal valid params"
        );
    }
}
