//! GeoLibre tool: resolve overlaps among polygons into a clean planar partition.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Remove Overlap (Multiple)* (Analysis).
//! `count_overlapping_features` *finds and counts* overlap regions but nothing
//! reallocates them; trade areas, service territories and buffered zones
//! routinely need mutual exclusivity. This flattens the inputs into disjoint
//! regions (the same incremental `BooleanOps` overlay used by
//! `count_overlapping_features`), then reassigns every multi-contributor region
//! to exactly one input so the outputs tile the union with no gaps and no
//! overlaps, conserving total area.
//!
//! Two division rules:
//!
//! * `center_line` (default) — each point of an overlap goes to the contributor
//!   it lies *deepest* inside (largest distance to that contributor's boundary),
//!   so the dividing line is the equal-depth locus between the two originals,
//!   ArcGIS's centre-line behaviour.
//! * `thiessen` — each point goes to the contributor whose generator centroid is
//!   nearest, giving straight Thiessen boundaries.
//!
//! Overlap regions are split on a grid (`grid_resolution` cells across the
//! longer side); winner cells are merged into row rectangles, unioned, clipped
//! to the region, and any residual sliver is folded into the largest winner so
//! area is conserved exactly. Each input feature yields at most one output
//! feature carrying its original attributes.

use std::collections::BTreeMap;

use geo::{
    Area, BooleanOps, BoundingRect, Centroid, Coord as GeoCoord, LineString, MultiPolygon, Point,
    Polygon,
};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Areas below this are treated as numerical slivers, not real geometry.
const SLIVER_AREA_EPS: f64 = 1e-9;
/// Default grid cells across the longer side of an overlap region.
const DEFAULT_GRID: usize = 40;
/// Clamp on grid cells per axis, to bound runtime.
const MAX_GRID: usize = 200;

#[derive(Clone, Copy, PartialEq)]
enum Method {
    CenterLine,
    Thiessen,
}

pub struct RemoveOverlapMultipleTool;

impl Tool for RemoveOverlapMultipleTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "remove_overlap_multiple",
            display_name: "Remove Overlap (Multiple)",
            summary: "Reallocate every overlap among a set of polygons to exactly one feature, producing a gap-free, overlap-free partition of their union (centre-line or Thiessen division), like ArcGIS Remove Overlap (Multiple).",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon vector layer (features may overlap).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'center_line' (assign each overlap point to the feature it lies deepest inside, default) or 'thiessen' (nearest generator centroid).",
                    required: false,
                },
                ToolParamSpec {
                    name: "grid_resolution",
                    description: "Cells across the longer side of each overlap region when dividing (default 40; higher = smoother dividing lines, slower).",
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

        // Collect input polygons (with the source feature index for attributes).
        let mut inputs: Vec<(usize, MultiPolygon)> = Vec::new();
        for (fidx, feature) in layer.features.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let Some(mp) = to_multipolygon(geom) else {
                continue;
            };
            if mp.unsigned_area() <= 0.0 {
                continue;
            }
            inputs.push((fidx, mp));
        }
        if inputs.is_empty() {
            return Err(ToolError::Execution(
                "no polygon features in input".to_string(),
            ));
        }

        // Precompute each contributor's boundary segments and centroid, keyed by
        // its position in `inputs`.
        let centroids: Vec<Point> = inputs
            .iter()
            .map(|(_, mp)| mp.centroid().unwrap_or(Point::new(0.0, 0.0)))
            .collect();
        let bboxes: Vec<[f64; 4]> = inputs.iter().map(|(_, mp)| bbox(mp)).collect();

        ctx.progress
            .info(&format!("overlaying {} polygon(s)", inputs.len()));

        // ── Incremental disjoint-region overlay (contributor = index in inputs) ──
        let mut regions: Vec<Region> = Vec::new();
        for (ci, (_, poly)) in inputs.iter().enumerate() {
            let p_bbox = bboxes[ci];
            let mut remaining = poly.clone();
            let mut next: Vec<Region> = Vec::with_capacity(regions.len() + 1);
            for r in regions.drain(..) {
                if !bbox_overlap(&r.bbox, &p_bbox) {
                    next.push(r);
                    continue;
                }
                let inter = r.geom.intersection(poly);
                if inter.unsigned_area() > SLIVER_AREA_EPS {
                    let mut ids = r.ids.clone();
                    ids.push(ci);
                    next.push(Region::new(inter, ids));
                    let diff = r.geom.difference(poly);
                    if diff.unsigned_area() > SLIVER_AREA_EPS {
                        next.push(Region::new(diff, r.ids));
                    }
                    remaining = remaining.difference(&r.geom);
                } else {
                    next.push(r);
                }
            }
            if remaining.unsigned_area() > SLIVER_AREA_EPS {
                next.push(Region::new(remaining, vec![ci]));
            }
            regions = next;
        }

        let overlap_regions = regions.iter().filter(|r| r.ids.len() > 1).count();
        ctx.progress.info(&format!(
            "{} disjoint region(s), {overlap_regions} shared; reallocating",
            regions.len()
        ));

        // ── Reallocate each region to a single contributor ───────────────────────
        // Accumulate winning pieces per contributor, then union at the end.
        let mut pieces: Vec<Vec<MultiPolygon>> = vec![Vec::new(); inputs.len()];
        for region in &regions {
            if region.ids.len() == 1 {
                pieces[region.ids[0]].push(region.geom.clone());
                continue;
            }
            for (winner, piece) in split_region(region, &inputs, &centroids, prm.method, prm.grid) {
                if piece.unsigned_area() > SLIVER_AREA_EPS {
                    pieces[winner].push(piece);
                }
            }
        }

        // ── Build output: one feature per input that kept any area ───────────────
        let mut out = Layer::new("no_overlap").with_geom_type(GeometryType::MultiPolygon);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for field in layer.schema.fields() {
            out.add_field(field.clone());
        }

        let mut kept = 0usize;
        for (ci, (src_idx, _)) in inputs.iter().enumerate() {
            let parts = &pieces[ci];
            if parts.is_empty() {
                continue;
            }
            let mut merged = parts[0].clone();
            for p in &parts[1..] {
                merged = merged.union(p);
            }
            if merged.unsigned_area() <= SLIVER_AREA_EPS {
                continue;
            }
            out.push(wbvector::Feature {
                fid: 0,
                geometry: Some(multipolygon_to_geometry(&merged)),
                attributes: layer.features[*src_idx].attributes.clone(),
            });
            kept += 1;
        }

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_count".to_string(), json!(inputs.len()));
        outputs.insert("output_count".to_string(), json!(kept));
        outputs.insert("overlap_regions".to_string(), json!(overlap_regions));
        Ok(ToolRunResult { outputs })
    }
}

// ── Region split ────────────────────────────────────────────────────────────

/// Divide a multi-contributor `region` among its contributors, returning one
/// clipped `MultiPolygon` per winner. Grid the region, assign each cell to a
/// winner, merge same-winner cells along rows into rectangles, union and clip
/// each winner to the region, then fold any residual sliver into the largest
/// winner so the pieces exactly retile the region.
fn split_region(
    region: &Region,
    inputs: &[(usize, MultiPolygon)],
    centroids: &[Point],
    method: Method,
    grid: usize,
) -> Vec<(usize, MultiPolygon)> {
    let contributors = &region.ids;
    let [minx, miny, maxx, maxy] = region.bbox;
    let w = maxx - minx;
    let h = maxy - miny;
    if !(w > 0.0 && h > 0.0) {
        return vec![(contributors[0], region.geom.clone())];
    }
    let (nx, ny) = if w >= h {
        let nx = grid.clamp(1, MAX_GRID);
        (nx, ((nx as f64 * h / w).ceil() as usize).clamp(1, MAX_GRID))
    } else {
        let ny = grid.clamp(1, MAX_GRID);
        (((ny as f64 * w / h).ceil() as usize).clamp(1, MAX_GRID), ny)
    };
    let cw = w / nx as f64;
    let ch = h / ny as f64;

    // Per-winner rectangles from same-winner runs within each row.
    let mut rects: Vec<Vec<Polygon>> = vec![Vec::new(); contributors.len()];
    for j in 0..ny {
        let y0 = miny + j as f64 * ch;
        let cy = y0 + ch * 0.5;
        let mut run_winner: Option<usize> = None;
        let mut run_x0 = minx;
        for i in 0..nx {
            let cx = minx + (i as f64 + 0.5) * cw;
            let local = if point_in_mp(cx, cy, &region.geom) {
                Some(pick_winner(cx, cy, contributors, inputs, centroids, method))
            } else {
                None
            };
            if local != run_winner {
                if let Some(k) = run_winner {
                    rects[k].push(rect(run_x0, y0, minx + i as f64 * cw, y0 + ch));
                }
                run_winner = local;
                run_x0 = minx + i as f64 * cw;
            }
        }
        if let Some(k) = run_winner {
            rects[k].push(rect(run_x0, y0, maxx, y0 + ch));
        }
    }

    // Union each winner's rectangles, then clip to the region.
    let mut out: Vec<(usize, MultiPolygon)> = Vec::new();
    let mut covered = MultiPolygon(Vec::new());
    for (k, rs) in rects.iter().enumerate() {
        if rs.is_empty() {
            continue;
        }
        let mut mp = MultiPolygon(vec![rs[0].clone()]);
        for r in &rs[1..] {
            mp = mp.union(&MultiPolygon(vec![r.clone()]));
        }
        let clipped = mp.intersection(&region.geom);
        if clipped.unsigned_area() > SLIVER_AREA_EPS {
            covered = covered.union(&clipped);
            out.push((contributors[k], clipped));
        }
    }

    if out.is_empty() {
        // Degenerate (region thinner than one cell): give it to contributor 0.
        return vec![(contributors[0], region.geom.clone())];
    }

    // Fold any residual (cells whose centres fell outside the region, boundary
    // stair-steps) into the largest winner so area is conserved exactly.
    let residual = region.geom.difference(&covered);
    if residual.unsigned_area() > SLIVER_AREA_EPS {
        let big = out
            .iter()
            .enumerate()
            .max_by(|a, b| {
                a.1 .1
                    .unsigned_area()
                    .partial_cmp(&b.1 .1.unsigned_area())
                    .unwrap()
            })
            .map(|(idx, _)| idx)
            .unwrap();
        out[big].1 = out[big].1.union(&residual);
    }
    out
}

/// Choose which contributor a point belongs to.
fn pick_winner(
    x: f64,
    y: f64,
    contributors: &[usize],
    inputs: &[(usize, MultiPolygon)],
    centroids: &[Point],
    method: Method,
) -> usize {
    let mut best = 0usize;
    match method {
        // Deepest containment: largest distance to the contributor's boundary.
        // Depth ties (a point equidistant from both boundaries, e.g. near an
        // edge shared by both originals) are broken by nearest generator
        // centroid so the divide stays a straight centre line, not biased to
        // whichever contributor was overlaid first.
        Method::CenterLine => {
            let mut best_depth = f64::NEG_INFINITY;
            let mut best_cdist = f64::INFINITY;
            for (k, &ci) in contributors.iter().enumerate() {
                let d = dist_to_boundary(x, y, &inputs[ci].1);
                let c = centroids[ci];
                let cd = (x - c.x()).powi(2) + (y - c.y()).powi(2);
                let deeper = d > best_depth + 1e-9;
                let tie_closer = (d - best_depth).abs() <= 1e-9 && cd < best_cdist;
                if deeper || tie_closer {
                    best_depth = d;
                    best_cdist = cd;
                    best = k;
                }
            }
        }
        // Nearest generator centroid.
        Method::Thiessen => {
            let mut best_d = f64::INFINITY;
            for (k, &ci) in contributors.iter().enumerate() {
                let c = centroids[ci];
                let d = (x - c.x()).powi(2) + (y - c.y()).powi(2);
                if d < best_d {
                    best_d = d;
                    best = k;
                }
            }
        }
    }
    best
}

/// Minimum distance from a point to the boundary of a MultiPolygon.
fn dist_to_boundary(x: f64, y: f64, mp: &MultiPolygon) -> f64 {
    let mut best = f64::INFINITY;
    for poly in &mp.0 {
        for ring in std::iter::once(poly.exterior()).chain(poly.interiors()) {
            let pts = &ring.0;
            for w in pts.windows(2) {
                let d = dist_point_seg(x, y, w[0].x, w[0].y, w[1].x, w[1].y);
                if d < best {
                    best = d;
                }
            }
        }
    }
    best
}

fn dist_point_seg(px: f64, py: f64, ax: f64, ay: f64, bx: f64, by: f64) -> f64 {
    let dx = bx - ax;
    let dy = by - ay;
    let len2 = dx * dx + dy * dy;
    let t = if len2 <= 0.0 {
        0.0
    } else {
        (((px - ax) * dx + (py - ay) * dy) / len2).clamp(0.0, 1.0)
    };
    let cx = ax + t * dx;
    let cy = ay + t * dy;
    ((px - cx).powi(2) + (py - cy).powi(2)).sqrt()
}

/// Even-odd point-in-MultiPolygon test (exterior minus holes).
fn point_in_mp(x: f64, y: f64, mp: &MultiPolygon) -> bool {
    for poly in &mp.0 {
        if point_in_ring(x, y, poly.exterior())
            && !poly.interiors().iter().any(|h| point_in_ring(x, y, h))
        {
            return true;
        }
    }
    false
}

fn point_in_ring(x: f64, y: f64, ring: &LineString) -> bool {
    let pts = &ring.0;
    let mut inside = false;
    let n = pts.len();
    if n < 3 {
        return false;
    }
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = (pts[i].x, pts[i].y);
        let (xj, yj) = (pts[j].x, pts[j].y);
        if (yi > y) != (yj > y) {
            let xcross = (xj - xi) * (y - yi) / (yj - yi) + xi;
            if x < xcross {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
    Polygon::new(
        LineString::new(vec![
            GeoCoord { x: x0, y: y0 },
            GeoCoord { x: x1, y: y0 },
            GeoCoord { x: x1, y: y1 },
            GeoCoord { x: x0, y: y1 },
            GeoCoord { x: x0, y: y0 },
        ]),
        vec![],
    )
}

// ── Regions ──────────────────────────────────────────────────────────────────

struct Region {
    geom: MultiPolygon,
    ids: Vec<usize>,
    bbox: [f64; 4],
}

impl Region {
    fn new(geom: MultiPolygon, ids: Vec<usize>) -> Region {
        let bbox = bbox(&geom);
        Region { geom, ids, bbox }
    }
}

fn bbox(mp: &MultiPolygon) -> [f64; 4] {
    match mp.bounding_rect() {
        Some(r) => [r.min().x, r.min().y, r.max().x, r.max().y],
        None => [
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        ],
    }
}

fn bbox_overlap(a: &[f64; 4], b: &[f64; 4]) -> bool {
    a[0] <= b[2] && b[0] <= a[2] && a[1] <= b[3] && b[1] <= a[3]
}

// ── geo <-> wbvector conversion ───────────────────────────────────────────────

fn to_multipolygon(geom: &wbvector::Geometry) -> Option<MultiPolygon> {
    match geom {
        wbvector::Geometry::Polygon {
            exterior,
            interiors,
        } => Some(MultiPolygon(vec![rings_to_polygon(exterior, interiors)])),
        wbvector::Geometry::MultiPolygon(parts) => Some(MultiPolygon(
            parts.iter().map(|(e, i)| rings_to_polygon(e, i)).collect(),
        )),
        _ => None,
    }
}

fn rings_to_polygon(exterior: &Ring, interiors: &[Ring]) -> Polygon {
    Polygon::new(
        ring_to_linestring(exterior),
        interiors.iter().map(ring_to_linestring).collect(),
    )
}

fn ring_to_linestring(ring: &Ring) -> LineString {
    LineString::new(
        ring.coords()
            .iter()
            .map(|c| GeoCoord { x: c.x, y: c.y })
            .collect(),
    )
}

fn multipolygon_to_geometry(mp: &MultiPolygon) -> wbvector::Geometry {
    wbvector::Geometry::MultiPolygon(mp.0.iter().map(polygon_to_rings).collect())
}

fn polygon_to_rings(poly: &Polygon) -> (Ring, Vec<Ring>) {
    (
        linestring_to_ring(poly.exterior()),
        poly.interiors().iter().map(linestring_to_ring).collect(),
    )
}

fn linestring_to_ring(ls: &LineString) -> Ring {
    let mut coords: Vec<Coord> = ls.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
    if coords.len() >= 2 && coords.first() == coords.last() {
        coords.pop();
    }
    Ring::new(coords)
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    method: Method,
    grid: usize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let method = match parse_optional_str(args, "method")?.map(|s| s.trim().to_lowercase()) {
        None => Method::CenterLine,
        Some(s)
            if s.is_empty() || s == "center_line" || s == "centerline" || s == "centre_line" =>
        {
            Method::CenterLine
        }
        Some(s) if s == "thiessen" || s == "voronoi" => Method::Thiessen,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "'method' must be 'center_line' or 'thiessen', got '{other}'"
            )))
        }
    };
    let grid = match args.get("grid_resolution") {
        None | Some(Value::Null) => DEFAULT_GRID,
        Some(Value::Number(n)) => (n.as_u64().unwrap_or(DEFAULT_GRID as u64) as usize).max(1),
        Some(Value::String(s)) if s.trim().is_empty() => DEFAULT_GRID,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'grid_resolution' must be an integer".to_string()))?
            .max(1),
        Some(_) => {
            return Err(ToolError::Validation(
                "'grid_resolution' must be a number".to_string(),
            ))
        }
    };
    Ok(Params { method, grid })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, Geometry};

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

    fn layer_of(named: Vec<(&str, Geometry)>) -> String {
        let mut l = Layer::new("polys")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        for (name, g) in named {
            l.add_feature(Some(g), &[("name", name.into())]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = RemoveOverlapMultipleTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn geom_area(geom: &Geometry) -> f64 {
        to_multipolygon(geom)
            .map(|m| m.unsigned_area())
            .unwrap_or(0.0)
    }

    /// Total output area equals the union area; outputs do not overlap.
    #[test]
    fn conserves_union_area_and_removes_overlap() {
        let a = square(0.0, 0.0, 10.0);
        let b = square(5.0, 0.0, 10.0);
        let union_area = {
            let ma = to_multipolygon(&a).unwrap();
            let mb = to_multipolygon(&b).unwrap();
            ma.union(&mb).unsigned_area()
        };
        let input = layer_of(vec![("A", a), ("B", b)]);
        let (out, layer) = run(json!({ "input": input }));
        assert_eq!(out.outputs["output_count"], json!(2));
        assert_eq!(out.outputs["overlap_regions"], json!(1));

        let mps: Vec<MultiPolygon> = layer
            .iter()
            .map(|f| to_multipolygon(f.geometry.as_ref().unwrap()).unwrap())
            .collect();
        let total: f64 = mps.iter().map(|m| m.unsigned_area()).sum();
        assert!(
            (total - union_area).abs() < union_area * 1e-3,
            "total {total} != union {union_area}"
        );
        // Pairwise overlap must be ~zero.
        let overlap = mps[0].intersection(&mps[1]).unsigned_area();
        assert!(overlap < union_area * 1e-3, "residual overlap {overlap}");
    }

    /// The two symmetric squares split their overlap roughly in half.
    #[test]
    fn center_line_splits_evenly() {
        let a = square(0.0, 0.0, 10.0);
        let b = square(5.0, 0.0, 10.0); // overlap x in [5,10], area 50
        let input = layer_of(vec![("A", a), ("B", b)]);
        let (_out, layer) = run(json!({ "input": input, "grid_resolution": 60 }));
        let mut areas = BTreeMap::new();
        let nidx = layer.schema.field_index("name").unwrap();
        for f in layer.iter() {
            let name = f.attributes[nidx].as_str().unwrap().to_string();
            areas.insert(name, geom_area(f.geometry.as_ref().unwrap()));
        }
        // Each original owns 50 exclusive + ~25 of the 50 overlap = ~75.
        let a = areas["A"];
        let b = areas["B"];
        assert!((a - 75.0).abs() < 5.0, "A area {a} not ~75");
        assert!((b - 75.0).abs() < 5.0, "B area {b} not ~75");
    }

    /// Attributes survive to the output.
    #[test]
    fn preserves_attributes() {
        let a = square(0.0, 0.0, 10.0);
        let b = square(5.0, 0.0, 10.0);
        let input = layer_of(vec![("east", a), ("west", b)]);
        let (_out, layer) = run(json!({ "input": input }));
        let nidx = layer.schema.field_index("name").unwrap();
        let names: Vec<String> = layer
            .iter()
            .map(|f| f.attributes[nidx].as_str().unwrap().to_string())
            .collect();
        assert!(names.contains(&"east".to_string()));
        assert!(names.contains(&"west".to_string()));
    }

    /// Disjoint inputs pass through unchanged (area preserved, still two features).
    #[test]
    fn disjoint_inputs_unchanged() {
        let a = square(0.0, 0.0, 5.0);
        let b = square(100.0, 100.0, 5.0);
        let input = layer_of(vec![("A", a), ("B", b)]);
        let (out, layer) = run(json!({ "input": input }));
        assert_eq!(out.outputs["overlap_regions"], json!(0));
        assert_eq!(layer.iter().count(), 2);
        let total: f64 = layer
            .iter()
            .map(|f| geom_area(f.geometry.as_ref().unwrap()))
            .sum();
        assert!((total - 50.0).abs() < 1e-6, "area {total} != 50");
    }

    #[test]
    fn thiessen_method_runs() {
        let a = square(0.0, 0.0, 10.0);
        let b = square(5.0, 0.0, 10.0);
        let input = layer_of(vec![("A", a), ("B", b)]);
        let (out, _l) = run(json!({ "input": input, "method": "thiessen" }));
        assert_eq!(out.outputs["output_count"], json!(2));
    }

    #[test]
    fn rejects_missing_input() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(RemoveOverlapMultipleTool.validate(&args).is_err());
    }

    #[test]
    fn rejects_bad_method() {
        let a = square(0.0, 0.0, 10.0);
        let input = layer_of(vec![("A", a)]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "method": "bogus" })).unwrap();
        assert!(RemoveOverlapMultipleTool.validate(&args).is_err());
    }
}
