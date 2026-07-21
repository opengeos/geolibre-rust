//! GeoLibre tool: displace buildings that conflict with symbolized roads.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Resolve Building Conflicts*
//! (Cartography). It completes the generalization family alongside the shipped
//! `regularize_building_footprints`, `delineate_built_up_areas`,
//! `thin_road_network`, and `collapse_dual_lines_to_centerline` — no bundled
//! tool performs displacement cartography.
//!
//! For small-scale mapping, building footprints that would graphically collide
//! with road symbology (the line widened to `barrier_width`, plus a `gap`) are
//! pushed away from the nearest barrier; building–building overlaps are then
//! relaxed apart. A building that still cannot be placed is shrunk about its
//! centroid down to `min_size`, and if it still conflicts it is hidden. A few
//! relaxation passes interleave barrier displacement and mutual spacing.
//!
//! Each output building carries a `status` of `unchanged` / `displaced` /
//! `shrunk` / `hidden`. The algorithm is deterministic (fixed feature order).

use std::collections::BTreeMap;

use geo::{
    Area, BooleanOps, Buffer, Coord as GeoCoord, Intersects, LineString, MultiLineString,
    MultiPolygon, Polygon,
};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

const MAX_PASSES: usize = 40;
const SLIVER: f64 = 1e-9;

pub struct ResolveBuildingConflictsTool;

impl Tool for ResolveBuildingConflictsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "resolve_building_conflicts",
            display_name: "Resolve Building Conflicts",
            summary: "Displace, shrink, or hide building footprints that graphically conflict with symbolized road barriers (and each other) for small-scale mapping, tagging each with a status, like ArcGIS Resolve Building Conflicts.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "buildings",
                    description: "Input building polygon layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "barriers",
                    description: "Barrier line layer (e.g. roads) buildings must clear.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output building layer with a 'status' field. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "barrier_width",
                    description: "Symbolized barrier width in map units (buffered by half this each side). Default 0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "gap",
                    description: "Minimum spacing to keep between a building and a barrier (map units). Default 0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_size",
                    description: "Minimum building area (map units²) before a conflicting building is hidden rather than shrunk. Default 0 (never hide).",
                    required: false,
                },
                ToolParamSpec {
                    name: "hide",
                    description: "Whether to drop hidden buildings from the output (default false: keep them with status='hidden').",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["buildings", "barriers"] {
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
        let buildings_path = args.get("buildings").and_then(Value::as_str).unwrap();
        let barriers_path = args.get("barriers").and_then(Value::as_str).unwrap();
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let blayer = load_input_layer(buildings_path)?;
        let barlayer = load_input_layer(barriers_path)?;

        // Forbidden zone = barrier lines buffered by (barrier_width/2 + gap).
        let mut lines: Vec<LineString> = Vec::new();
        for f in barlayer.iter() {
            if let Some(geom) = f.geometry.as_ref() {
                collect_lines(geom, &mut lines);
            }
        }
        let clearance = prm.barrier_width * 0.5 + prm.gap;
        let forbidden: MultiPolygon = if lines.is_empty() || clearance <= 0.0 {
            MultiPolygon(Vec::new())
        } else {
            MultiLineString(lines.clone()).buffer(clearance)
        };
        // Barrier segments for displacement directions.
        let segs: Vec<[f64; 4]> = lines
            .iter()
            .flat_map(|ls| ls.0.windows(2).map(|w| [w[0].x, w[0].y, w[1].x, w[1].y]))
            .collect();

        // Load buildings as mutable geometries.
        let mut builds: Vec<Build> = Vec::new();
        for (fidx, f) in blayer.features.iter().enumerate() {
            let Some(geom) = f.geometry.as_ref() else {
                continue;
            };
            let Some(poly) = to_polygon(geom) else {
                continue;
            };
            if poly.unsigned_area() <= 0.0 {
                continue;
            }
            builds.push(Build {
                src: fidx,
                poly,
                status: Status::Unchanged,
            });
        }
        if builds.is_empty() {
            return Err(ToolError::Execution(
                "no building polygons in input".to_string(),
            ));
        }

        ctx.progress.info(&format!(
            "resolving conflicts for {} building(s)",
            builds.len()
        ));

        // ── Relaxation passes: barrier displacement + mutual spacing ─────────────
        for _ in 0..MAX_PASSES {
            let mut moved = false;
            // Push each conflicting building away from the nearest barrier.
            for build in builds.iter_mut() {
                if build.status == Status::Hidden {
                    continue;
                }
                if !forbidden.0.is_empty() && forbidden.intersects(&build.poly) {
                    let (cx, cy) = centroid(&build.poly);
                    if let Some((nx, ny)) = push_dir(cx, cy, &segs) {
                        let step = building_radius(&build.poly).max(clearance) * 0.5;
                        translate(&mut build.poly, nx * step, ny * step);
                        if build.status == Status::Unchanged {
                            build.status = Status::Displaced;
                        }
                        moved = true;
                    }
                }
            }
            // Relax building-building overlaps.
            for i in 0..builds.len() {
                for j in (i + 1)..builds.len() {
                    if builds[i].status == Status::Hidden || builds[j].status == Status::Hidden {
                        continue;
                    }
                    if builds[i].poly.intersects(&builds[j].poly)
                        && builds[i].poly.intersection(&builds[j].poly).unsigned_area() > SLIVER
                    {
                        let (cix, ciy) = centroid(&builds[i].poly);
                        let (cjx, cjy) = centroid(&builds[j].poly);
                        let (mut dx, mut dy) = (cjx - cix, cjy - ciy);
                        let len = (dx * dx + dy * dy).sqrt();
                        if len < 1e-9 {
                            dx = 1.0;
                            dy = 0.0;
                        } else {
                            dx /= len;
                            dy /= len;
                        }
                        let step = (building_radius(&builds[i].poly)
                            + building_radius(&builds[j].poly))
                            * 0.25;
                        translate(&mut builds[i].poly, -dx * step, -dy * step);
                        translate(&mut builds[j].poly, dx * step, dy * step);
                        if builds[i].status == Status::Unchanged {
                            builds[i].status = Status::Displaced;
                        }
                        if builds[j].status == Status::Unchanged {
                            builds[j].status = Status::Displaced;
                        }
                        moved = true;
                    }
                }
            }
            if !moved {
                break;
            }
        }

        // ── Last resort: shrink, then hide, buildings still on a barrier ─────────
        for b in builds.iter_mut() {
            if b.status == Status::Hidden || forbidden.0.is_empty() {
                continue;
            }
            let mut guard = 0;
            while forbidden.intersects(&b.poly) && guard < 20 {
                if b.poly.unsigned_area() <= prm.min_size {
                    b.status = Status::Hidden;
                    break;
                }
                shrink(&mut b.poly, 0.8);
                b.status = Status::Shrunk;
                guard += 1;
            }
            if guard >= 20 && forbidden.intersects(&b.poly) {
                b.status = Status::Hidden;
            }
        }

        // ── Build output ─────────────────────────────────────────────────────────
        let mut out = Layer::new("resolved_buildings").with_geom_type(GeometryType::Polygon);
        if let Some(epsg) = blayer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for field in blayer.schema.fields() {
            out.add_field(field.clone());
        }
        out.add_field(FieldDef::new("status", FieldType::Text));

        let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
        let mut kept = 0usize;
        for b in &builds {
            let status = b.status.name();
            *counts.entry(status).or_default() += 1;
            if b.status == Status::Hidden && prm.hide {
                continue;
            }
            let mut attrs = blayer.features[b.src].attributes.clone();
            attrs.push(FieldValue::Text(status.to_string()));
            out.push(wbvector::Feature {
                fid: 0,
                geometry: Some(polygon_to_geometry(&b.poly)),
                attributes: attrs,
            });
            kept += 1;
        }

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("building_count".to_string(), json!(builds.len()));
        outputs.insert("output_count".to_string(), json!(kept));
        outputs.insert(
            "displaced".to_string(),
            json!(counts.get("displaced").copied().unwrap_or(0)),
        );
        outputs.insert(
            "shrunk".to_string(),
            json!(counts.get("shrunk").copied().unwrap_or(0)),
        );
        outputs.insert(
            "hidden".to_string(),
            json!(counts.get("hidden").copied().unwrap_or(0)),
        );
        Ok(ToolRunResult { outputs })
    }
}

// ── Building model ─────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq)]
enum Status {
    Unchanged,
    Displaced,
    Shrunk,
    Hidden,
}
impl Status {
    fn name(self) -> &'static str {
        match self {
            Status::Unchanged => "unchanged",
            Status::Displaced => "displaced",
            Status::Shrunk => "shrunk",
            Status::Hidden => "hidden",
        }
    }
}

struct Build {
    src: usize,
    poly: Polygon,
    status: Status,
}

fn centroid(p: &Polygon) -> (f64, f64) {
    // Vertex centroid of the exterior (stable, cheap; good enough for a push).
    let pts = &p.exterior().0;
    let n = pts.len().max(1) as f64;
    let (sx, sy) = pts
        .iter()
        .fold((0.0, 0.0), |(ax, ay), c| (ax + c.x, ay + c.y));
    (sx / n, sy / n)
}

fn building_radius(p: &Polygon) -> f64 {
    let (cx, cy) = centroid(p);
    p.exterior()
        .0
        .iter()
        .map(|c| ((c.x - cx).powi(2) + (c.y - cy).powi(2)).sqrt())
        .fold(0.0, f64::max)
        .max(1e-6)
}

fn translate(p: &mut Polygon, dx: f64, dy: f64) {
    let ext: Vec<GeoCoord> = p
        .exterior()
        .0
        .iter()
        .map(|c| GeoCoord {
            x: c.x + dx,
            y: c.y + dy,
        })
        .collect();
    let holes: Vec<LineString> = p
        .interiors()
        .iter()
        .map(|r| {
            LineString::new(
                r.0.iter()
                    .map(|c| GeoCoord {
                        x: c.x + dx,
                        y: c.y + dy,
                    })
                    .collect(),
            )
        })
        .collect();
    *p = Polygon::new(LineString::new(ext), holes);
}

fn shrink(p: &mut Polygon, factor: f64) {
    let (cx, cy) = centroid(p);
    let scale = |c: &GeoCoord| GeoCoord {
        x: cx + (c.x - cx) * factor,
        y: cy + (c.y - cy) * factor,
    };
    let ext = LineString::new(p.exterior().0.iter().map(scale).collect());
    let holes: Vec<LineString> = p
        .interiors()
        .iter()
        .map(|r| LineString::new(r.0.iter().map(scale).collect()))
        .collect();
    *p = Polygon::new(ext, holes);
}

/// Unit vector pushing a point away from its nearest barrier segment.
fn push_dir(x: f64, y: f64, segs: &[[f64; 4]]) -> Option<(f64, f64)> {
    let mut best = f64::INFINITY;
    let mut np = (0.0, 0.0);
    for s in segs {
        let (px, py) = nearest_on_seg(x, y, s[0], s[1], s[2], s[3]);
        let d = (x - px).powi(2) + (y - py).powi(2);
        if d < best {
            best = d;
            np = (px, py);
        }
    }
    if !best.is_finite() {
        return None;
    }
    let (mut dx, mut dy) = (x - np.0, y - np.1);
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-9 {
        // On the barrier: push perpendicular to the nearest segment.
        return Some((0.0, 1.0));
    }
    dx /= len;
    dy /= len;
    Some((dx, dy))
}

fn nearest_on_seg(px: f64, py: f64, ax: f64, ay: f64, bx: f64, by: f64) -> (f64, f64) {
    let dx = bx - ax;
    let dy = by - ay;
    let len2 = dx * dx + dy * dy;
    let t = if len2 <= 0.0 {
        0.0
    } else {
        (((px - ax) * dx + (py - ay) * dy) / len2).clamp(0.0, 1.0)
    };
    (ax + t * dx, ay + t * dy)
}

// ── geo <-> wbvector ───────────────────────────────────────────────────────────

fn to_polygon(geom: &Geometry) -> Option<Polygon> {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => Some(Polygon::new(
            ring_to_ls(exterior),
            interiors.iter().map(ring_to_ls).collect(),
        )),
        Geometry::MultiPolygon(parts) if !parts.is_empty() => {
            // Use the largest part as the representative footprint.
            let mut best: Option<Polygon> = None;
            let mut best_a = 0.0;
            for (e, i) in parts {
                let poly = Polygon::new(ring_to_ls(e), i.iter().map(ring_to_ls).collect());
                let a = poly.unsigned_area();
                if a > best_a {
                    best_a = a;
                    best = Some(poly);
                }
            }
            best
        }
        _ => None,
    }
}

fn ring_to_ls(r: &Ring) -> LineString {
    LineString::new(
        r.coords()
            .iter()
            .map(|c| GeoCoord { x: c.x, y: c.y })
            .collect(),
    )
}

fn collect_lines(geom: &Geometry, out: &mut Vec<LineString>) {
    match geom {
        Geometry::LineString(cs) => out.push(LineString::new(
            cs.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect(),
        )),
        Geometry::MultiLineString(parts) => {
            for cs in parts {
                out.push(LineString::new(
                    cs.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect(),
                ));
            }
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            out.push(ring_to_ls(exterior));
            for h in interiors {
                out.push(ring_to_ls(h));
            }
        }
        Geometry::MultiPolygon(parts) => {
            for (e, i) in parts {
                out.push(ring_to_ls(e));
                for h in i {
                    out.push(ring_to_ls(h));
                }
            }
        }
        _ => {}
    }
}

fn polygon_to_geometry(p: &Polygon) -> Geometry {
    let ext: Vec<Coord> = ls_to_coords(p.exterior());
    let holes: Vec<Ring> = p
        .interiors()
        .iter()
        .map(|r| Ring::new(ls_to_coords(r)))
        .collect();
    Geometry::Polygon {
        exterior: Ring::new(ext),
        interiors: holes,
    }
}

fn ls_to_coords(ls: &LineString) -> Vec<Coord> {
    let mut cs: Vec<Coord> = ls.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
    if cs.len() >= 2 && cs.first() == cs.last() {
        cs.pop();
    }
    cs
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    barrier_width: f64,
    gap: f64,
    min_size: f64,
    hide: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let barrier_width = f64_param(args, "barrier_width")?.unwrap_or(0.0).max(0.0);
    let gap = f64_param(args, "gap")?.unwrap_or(0.0).max(0.0);
    let min_size = f64_param(args, "min_size")?.unwrap_or(0.0).max(0.0);
    let hide = match args.get("hide") {
        None | Some(Value::Null) => false,
        Some(Value::Bool(b)) => *b,
        Some(Value::String(s)) => matches!(s.trim().to_lowercase().as_str(), "true" | "1" | "yes"),
        Some(_) => false,
    };
    if barrier_width == 0.0 && gap == 0.0 {
        return Err(ToolError::Validation(
            "set 'barrier_width' and/or 'gap' to define the clearance around barriers".to_string(),
        ));
    }
    Ok(Params {
        barrier_width,
        gap,
        min_size,
        hide,
    })
}

fn f64_param(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("'{key}' must be a number"))),
        Some(_) => Err(ToolError::Validation(format!("'{key}' must be a number"))),
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
        l.add_field(FieldDef::new("id", FieldType::Integer));
        for (n, (x, y, s)) in sqs.iter().enumerate() {
            l.add_feature(Some(square(*x, *y, *s)), &[("id", (n as i64).into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn line_layer(lines: &[Vec<(f64, f64)>]) -> String {
        let mut l = Layer::new("r")
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

    fn forbidden_of(barriers: &str, clearance: f64) -> MultiPolygon {
        let bar = load_input_layer(barriers).unwrap();
        let mut lines = Vec::new();
        for f in bar.iter() {
            if let Some(g) = f.geometry.as_ref() {
                collect_lines(g, &mut lines);
            }
        }
        MultiLineString(lines).buffer(clearance)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ResolveBuildingConflictsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// A building sitting on a road is displaced clear of the barrier buffer.
    #[test]
    fn displaces_conflicting_building() {
        // Road along y=0; building straddling it.
        let b = building_layer(&[(-5.0, -5.0, 10.0)]);
        let r = line_layer(&[vec![(-100.0, 0.0), (100.0, 0.0)]]);
        let (out, layer) = run(json!({
            "buildings": b, "barriers": r, "barrier_width": 4.0, "gap": 2.0
        }));
        assert_eq!(out.outputs["displaced"], json!(1));
        // The moved building must not intersect the forbidden zone.
        let forbidden = forbidden_of(&r, 4.0);
        let f = layer.iter().next().unwrap();
        let poly = to_polygon(f.geometry.as_ref().unwrap()).unwrap();
        assert!(!forbidden.intersects(&poly), "building still on barrier");
        let si = layer.schema.field_index("status").unwrap();
        assert_eq!(f.attributes[si].as_str(), Some("displaced"));
    }

    /// A building far from the road is left unchanged.
    #[test]
    fn leaves_clear_building_unchanged() {
        let b = building_layer(&[(0.0, 100.0, 10.0)]);
        let r = line_layer(&[vec![(-100.0, 0.0), (100.0, 0.0)]]);
        let (out, layer) = run(json!({
            "buildings": b, "barriers": r, "barrier_width": 4.0, "gap": 2.0
        }));
        assert_eq!(out.outputs["displaced"], json!(0));
        let si = layer.schema.field_index("status").unwrap();
        assert_eq!(
            layer.iter().next().unwrap().attributes[si].as_str(),
            Some("unchanged")
        );
    }

    /// After resolving, no output building intersects the barrier buffer.
    #[test]
    fn resolves_a_row_of_buildings() {
        // Several buildings along a road.
        let mut sqs = Vec::new();
        for i in 0..6 {
            sqs.push((i as f64 * 12.0 - 30.0, -4.0, 8.0));
        }
        let b = building_layer(&sqs);
        let r = line_layer(&[vec![(-100.0, 0.0), (100.0, 0.0)]]);
        let (_out, layer) = run(json!({
            "buildings": b, "barriers": r, "barrier_width": 3.0, "gap": 1.5
        }));
        let forbidden = forbidden_of(&r, 3.0);
        for f in layer.iter() {
            let poly = to_polygon(f.geometry.as_ref().unwrap()).unwrap();
            assert!(
                !forbidden.intersects(&poly),
                "a building remained on the barrier"
            );
        }
    }

    /// A building trapped against a barrier with a tiny min_size is shrunk/hidden.
    #[test]
    fn shrinks_or_hides_trapped_building() {
        // Building boxed by two near-parallel roads with no room to move.
        let b = building_layer(&[(-3.0, -3.0, 6.0)]);
        let r = line_layer(&[
            vec![(-100.0, 5.0), (100.0, 5.0)],
            vec![(-100.0, -5.0), (100.0, -5.0)],
        ]);
        let (out, _l) = run(json!({
            "buildings": b, "barriers": r, "barrier_width": 6.0, "gap": 2.0, "min_size": 1.0
        }));
        let shrunk = out.outputs["shrunk"].as_i64().unwrap();
        let hidden = out.outputs["hidden"].as_i64().unwrap();
        assert!(
            shrunk + hidden >= 1,
            "trapped building was shrunk or hidden"
        );
    }

    #[test]
    fn rejects_missing_clearance() {
        let b = building_layer(&[(0.0, 0.0, 10.0)]);
        let r = line_layer(&[vec![(0.0, 0.0), (10.0, 0.0)]]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "buildings": b, "barriers": r })).unwrap();
        assert!(ResolveBuildingConflictsTool.validate(&args).is_err());
    }
}
