//! GeoLibre tool: fixed-width least-cost corridor polygons connecting regions.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Optimal Corridor Connections*
//! (Spatial Analyst). Two adjacent tools exist and neither produces this
//! output:
//!
//!   * GeoLibre's `cost_connectivity` (ArcGIS *Optimal Region Connections*)
//!     returns least-cost **lines** between regions — a zero-width network.
//!   * GeoLibre's `corridor` returns an accumulated-cost **surface** for a
//!     single source/destination pair, leaving the user to threshold and
//!     vectorise it by hand, once per pair.
//!
//! The deliverable in habitat-connectivity and greenway planning is neither: it
//! is a set of corridor **polygons** of a specified width spanning the whole
//! region network at once. That is what this produces, on top of machinery the
//! repo already has — a Dijkstra cost-distance sweep per region, a minimum
//! spanning tree over the region graph to decide which pairs actually need a
//! corridor, and `geo`'s `BooleanOps` to union the swept corridor into clean,
//! non-overlapping polygons.

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::collections::BTreeMap;

use geo::{Area, BooleanOps, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::common::load_input_raster;
use crate::vector_common::{
    geometry_contains_point, load_input_layer, parse_optional_str, write_or_store_layer,
};

/// Segments per quarter turn when rounding a corridor end/joint.
const ARC_STEPS: usize = 8;

pub struct OptimalCorridorConnectionsTool;

impl Tool for OptimalCorridorConnectionsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "optimal_corridor_connections",
            display_name: "Optimal Corridor Connections",
            summary: "Compute the network of fixed-width least-cost corridor polygons that optimally connects every input region across a cost surface, plus the corridor centerlines. Unlike cost_connectivity (zero-width lines) or corridor (a single-pair cost surface), this emits ready-to-use corridor polygons for the whole region network. Like ArcGIS Optimal Corridor Connections.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Region polygon features to connect.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output corridor polygon path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cost_raster",
                    description: "Optional cost surface. Without it, cost is uniform (Euclidean).",
                    required: false,
                },
                ToolParamSpec {
                    name: "barriers",
                    description: "Optional polygon features treated as impassable.",
                    required: false,
                },
                ToolParamSpec {
                    name: "corridor_width",
                    description: "Corridor width in map units (default: 10 cells).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_lines",
                    description: "Optional path for the corridor centerlines.",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighbor_option",
                    description: "'spanning_tree' (default) keeps only the corridors needed to connect every region; 'all_pairs' emits a corridor for every region pair.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Analysis cell size when no cost raster is supplied (default: region extent / 200).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_neighbor_option(args)?;
        for key in ["corridor_width", "cell_size"] {
            if let Some(v) = parse_optional_f64(args, key)? {
                if v <= 0.0 {
                    return Err(ToolError::Validation(format!("'{key}' must be positive")));
                }
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let all_pairs = parse_neighbor_option(args)? == NeighborOption::AllPairs;

        let regions = load_input_layer(input)?;
        let region_geoms: Vec<&Geometry> =
            regions.iter().filter_map(|f| f.geometry.as_ref()).collect();
        if region_geoms.len() < 2 {
            return Err(ToolError::Execution(
                "at least 2 region features are required to build a corridor".to_string(),
            ));
        }

        // Barriers are loaded up front: with no cost raster the analysis extent
        // is derived from the inputs, and a barrier is exactly the thing a
        // corridor has to detour *around* — so the grid must be wide enough to
        // hold that detour, not just the regions.
        let barrier_geoms: Vec<Geometry> = match parse_optional_str(args, "barriers")? {
            Some(p) => load_input_layer(p)?
                .iter()
                .filter_map(|f| f.geometry.clone())
                .collect(),
            None => Vec::new(),
        };

        // Analysis grid: the cost raster if given, else one derived from the
        // region (and barrier) extent.
        let grid = match parse_optional_str(args, "cost_raster")? {
            Some(p) => load_input_raster(p)?,
            None => {
                let (mut min_x, mut min_y) = (f64::INFINITY, f64::INFINITY);
                let (mut max_x, mut max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
                for b in region_geoms
                    .iter()
                    .map(|g| g.bbox())
                    .chain(barrier_geoms.iter().map(|g| g.bbox()))
                    .flatten()
                {
                    min_x = min_x.min(b.min_x);
                    min_y = min_y.min(b.min_y);
                    max_x = max_x.max(b.max_x);
                    max_y = max_y.max(b.max_y);
                }
                if !min_x.is_finite() {
                    return Err(ToolError::Execution(
                        "regions have no usable extent".to_string(),
                    ));
                }
                // Pad so corridors are not clipped at the very edge.
                let cell = parse_optional_f64(args, "cell_size")?
                    .unwrap_or_else(|| ((max_x - min_x).max(max_y - min_y) / 200.0).max(1e-9));
                let pad = cell * 5.0;
                let (min_x, min_y) = (min_x - pad, min_y - pad);
                let (max_x, max_y) = (max_x + pad, max_y + pad);
                let cols = (((max_x - min_x) / cell).ceil() as usize).max(2);
                let rows = (((max_y - min_y) / cell).ceil() as usize).max(2);
                let mut crs = CrsInfo::default();
                if let Some(e) = regions.crs_epsg() {
                    crs = CrsInfo::from_epsg(e);
                }
                let mut r = Raster::new(RasterConfig {
                    cols,
                    rows,
                    bands: 1,
                    x_min: min_x,
                    y_min: max_y - rows as f64 * cell,
                    cell_size: cell,
                    cell_size_y: Some(cell),
                    nodata: -9999.0,
                    data_type: DataType::F32,
                    crs,
                    metadata: Default::default(),
                });
                for row in 0..rows {
                    for col in 0..cols {
                        r.set(0, row as isize, col as isize, 1.0).map_err(|e| {
                            ToolError::Execution(format!("failed building analysis grid: {e}"))
                        })?;
                    }
                }
                r
            }
        };
        let rows = grid.rows;
        let cols = grid.cols;
        let width = parse_optional_f64(args, "corridor_width")?
            .unwrap_or_else(|| 10.0 * grid.cell_size_x.max(grid.cell_size_y));

        // Impassable cells: barrier polygons plus cost NoData.
        let mut blocked = vec![false; rows * cols];
        for row in 0..rows {
            for col in 0..cols {
                let v = grid.get(0, row as isize, col as isize);
                if v == grid.nodata || !v.is_finite() || v < 0.0 {
                    blocked[row * cols + col] = true;
                }
            }
        }
        for g in &barrier_geoms {
            mark_cells(&grid, g, &mut blocked);
        }

        // Region seed cells.
        let mut seeds: Vec<Vec<usize>> = Vec::with_capacity(region_geoms.len());
        for g in &region_geoms {
            let mut mask = vec![false; rows * cols];
            mark_cells(&grid, g, &mut mask);
            let cells: Vec<usize> = mask
                .iter()
                .enumerate()
                .filter(|&(i, &m)| m && !blocked[i])
                .map(|(i, _)| i)
                .collect();
            seeds.push(cells);
        }
        if seeds.iter().filter(|s| !s.is_empty()).count() < 2 {
            return Err(ToolError::Execution(
                "fewer than 2 regions have passable cells on the analysis grid".to_string(),
            ));
        }
        ctx.progress
            .info(&format!("cost-distance from {} region(s)", seeds.len()));

        // Cost-distance surface + back-links per region.
        let mut fields: Vec<(Vec<f64>, Vec<usize>)> = Vec::with_capacity(seeds.len());
        for (i, s) in seeds.iter().enumerate() {
            fields.push(cost_distance(&grid, &blocked, s));
            ctx.progress.progress((i as f64 + 1.0) / seeds.len() as f64);
        }

        // Pair costs: the cheapest cell where the two cost fields meet.
        let n = seeds.len();
        let mut pair_best: Vec<Vec<Option<(f64, usize)>>> = vec![vec![None; n]; n];
        for a in 0..n {
            for b in (a + 1)..n {
                if seeds[a].is_empty() || seeds[b].is_empty() {
                    continue;
                }
                let mut best: Option<(f64, usize)> = None;
                for i in 0..(rows * cols) {
                    let (da, db) = (fields[a].0[i], fields[b].0[i]);
                    if !da.is_finite() || !db.is_finite() {
                        continue;
                    }
                    let total = da + db;
                    if best.is_none_or(|(bc, _)| total < bc) {
                        best = Some((total, i));
                    }
                }
                pair_best[a][b] = best;
                pair_best[b][a] = best;
            }
        }

        // Which pairs get a corridor: a spanning tree by default, so the output
        // is the minimum network that connects everything rather than O(n^2)
        // overlapping corridors.
        let edges = if all_pairs {
            let mut e = Vec::new();
            for (a, row) in pair_best.iter().enumerate() {
                for (b, entry) in row.iter().enumerate().skip(a + 1) {
                    if let Some((c, cell)) = entry {
                        e.push((a, b, *c, *cell));
                    }
                }
            }
            e
        } else {
            minimum_spanning_tree(&pair_best, n)
        };
        if edges.is_empty() {
            return Err(ToolError::Execution(
                "no region pair could be connected across the cost surface; \
                 barriers may fully separate them"
                    .to_string(),
            ));
        }

        let mut out = Layer::new("optimal_corridors").with_geom_type(GeometryType::Polygon);
        let mut lines = Layer::new("corridor_centerlines").with_geom_type(GeometryType::LineString);
        if let Some(epsg) = regions.crs_epsg() {
            out = out.with_crs_epsg(epsg);
            lines = lines.with_crs_epsg(epsg);
        }
        for l in [&mut out, &mut lines] {
            l.add_field(FieldDef::new("FROM_REGION", FieldType::Integer));
            l.add_field(FieldDef::new("TO_REGION", FieldType::Integer));
            l.add_field(FieldDef::new("ACCUM_COST", FieldType::Float));
            l.add_field(FieldDef::new("LENGTH", FieldType::Float));
        }
        out.add_field(FieldDef::new("AREA", FieldType::Float));

        let mut total_area = 0.0;
        let mut total_length = 0.0;
        let mut emitted = 0usize;

        for (a, b, cost, meet) in edges {
            // Trace the meeting cell back to both regions and join the halves.
            let mut path_a = trace_back(&fields[a].1, meet);
            path_a.reverse();
            let path_b = trace_back(&fields[b].1, meet);
            let mut cells = path_a;
            cells.extend(path_b.into_iter().skip(1));
            if cells.len() < 2 {
                continue;
            }
            let pts: Vec<(f64, f64)> = cells
                .iter()
                .map(|&i| pixel_center(&grid, i / cols, i % cols))
                .collect();
            let length: f64 = pts
                .windows(2)
                .map(|w| (w[1].0 - w[0].0).hypot(w[1].1 - w[0].1))
                .sum();

            let mp = thicken(&pts, width / 2.0);
            let area = mp.unsigned_area();
            total_area += area;
            total_length += length;
            emitted += 1;

            let attrs = [
                ("FROM_REGION", FieldValue::Integer(a as i64)),
                ("TO_REGION", FieldValue::Integer(b as i64)),
                ("ACCUM_COST", FieldValue::Float(cost)),
                ("LENGTH", FieldValue::Float(length)),
            ];
            let mut poly_attrs = attrs.to_vec();
            poly_attrs.push(("AREA", FieldValue::Float(area)));
            out.add_feature(Some(multipolygon_to_geometry(&mp)), &poly_attrs)
                .map_err(|e| ToolError::Execution(format!("failed adding corridor: {e}")))?;
            lines
                .add_feature(
                    Some(Geometry::line_string(
                        pts.iter().map(|&(x, y)| Coord::xy(x, y)).collect(),
                    )),
                    &attrs,
                )
                .map_err(|e| ToolError::Execution(format!("failed adding centerline: {e}")))?;
        }
        if emitted == 0 {
            return Err(ToolError::Execution(
                "no corridor could be traced between the regions".to_string(),
            ));
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("region_count".to_string(), json!(n));
        outputs.insert("corridor_count".to_string(), json!(emitted));
        outputs.insert("corridor_width".to_string(), json!(width));
        outputs.insert("total_corridor_area".to_string(), json!(total_area));
        outputs.insert("total_centerline_length".to_string(), json!(total_length));
        if matches!(args.get("output_lines"), Some(v) if !v.is_null()) {
            let p = parse_optional_str(args, "output_lines")?;
            outputs.insert("output_lines".to_string(), json!(write_or_store_layer(lines, p)?));
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── Cost distance ───────────────────────────────────────────────────────────

/// Multi-source Dijkstra over the 8-neighbourhood; returns the accumulated
/// cost per cell and the back-link index used to trace a route home.
fn cost_distance(grid: &Raster, blocked: &[bool], seeds: &[usize]) -> (Vec<f64>, Vec<usize>) {
    let rows = grid.rows;
    let cols = grid.cols;
    let n = rows * cols;
    let diag = grid.cell_size_x.hypot(grid.cell_size_y);
    let straight = grid.cell_size_x.min(grid.cell_size_y);

    let mut dist = vec![f64::INFINITY; n];
    let mut back = vec![usize::MAX; n];
    let mut done = vec![false; n];
    let mut heap: BinaryHeap<Node> = BinaryHeap::new();
    for &s in seeds {
        dist[s] = 0.0;
        heap.push(Node { cost: 0.0, idx: s });
    }

    while let Some(Node { cost, idx }) = heap.pop() {
        if done[idx] {
            continue;
        }
        done[idx] = true;
        let (r, c) = (idx / cols, idx % cols);
        let here = grid.get(0, r as isize, c as isize).max(0.0);
        for (dr, dc) in [
            (-1_isize, -1_isize), (-1, 0), (-1, 1),
            (0, -1),                       (0, 1),
            (1, -1),              (1, 0),  (1, 1),
        ] {
            let (nr, nc) = (r as isize + dr, c as isize + dc);
            if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                continue;
            }
            let nidx = nr as usize * cols + nc as usize;
            if done[nidx] || blocked[nidx] {
                continue;
            }
            let there = grid.get(0, nr, nc).max(0.0);
            let step = if dr != 0 && dc != 0 { diag } else { straight };
            // Mean cost across the step, the standard cost-distance accrual.
            let nd = cost + (here + there) / 2.0 * step;
            if nd < dist[nidx] {
                dist[nidx] = nd;
                back[nidx] = idx;
                heap.push(Node {
                    cost: nd,
                    idx: nidx,
                });
            }
        }
    }
    (dist, back)
}

fn trace_back(back: &[usize], from: usize) -> Vec<usize> {
    let mut path = vec![from];
    let mut cur = from;
    while back[cur] != usize::MAX {
        cur = back[cur];
        path.push(cur);
    }
    path
}

/// Minimum spanning tree over the region graph (Prim), returning the retained
/// edges as `(a, b, cost, meeting_cell)`.
fn minimum_spanning_tree(
    pair: &[Vec<Option<(f64, usize)>>],
    n: usize,
) -> Vec<(usize, usize, f64, usize)> {
    let mut in_tree = vec![false; n];
    let mut edges = Vec::new();
    in_tree[0] = true;
    for _ in 1..n {
        let mut best: Option<(usize, usize, f64, usize)> = None;
        for a in 0..n {
            if !in_tree[a] {
                continue;
            }
            for b in 0..n {
                if in_tree[b] {
                    continue;
                }
                if let Some((c, cell)) = pair[a][b] {
                    if best.is_none_or(|(_, _, bc, _)| c < bc) {
                        best = Some((a, b, c, cell));
                    }
                }
            }
        }
        match best {
            // Disconnected graph: stop rather than inventing an edge.
            None => break,
            Some(e) => {
                in_tree[e.1] = true;
                edges.push(e);
            }
        }
    }
    edges
}

// ── Corridor geometry ───────────────────────────────────────────────────────

/// Thickens a centerline into a corridor polygon of half-width `hw`, by
/// unioning a capsule per segment. `BooleanOps` merges them into one clean
/// outline, which is exactly what a round-joined buffer is.
fn thicken(pts: &[(f64, f64)], hw: f64) -> MultiPolygon {
    let mut acc = MultiPolygon(Vec::new());
    for w in pts.windows(2) {
        let (a, b) = (w[0], w[1]);
        let (dx, dy) = (b.0 - a.0, b.1 - a.1);
        let len = dx.hypot(dy);
        if len <= 0.0 {
            continue;
        }
        let (nx, ny) = (-dy / len * hw, dx / len * hw);
        let rect = Polygon::new(
            LineString::new(vec![
                GeoCoord { x: a.0 + nx, y: a.1 + ny },
                GeoCoord { x: b.0 + nx, y: b.1 + ny },
                GeoCoord { x: b.0 - nx, y: b.1 - ny },
                GeoCoord { x: a.0 - nx, y: a.1 - ny },
                GeoCoord { x: a.0 + nx, y: a.1 + ny },
            ]),
            vec![],
        );
        acc = acc.union(&MultiPolygon(vec![rect]));
    }
    // Round the joints and the two ends so the corridor has no notches.
    for &(x, y) in pts {
        acc = acc.union(&MultiPolygon(vec![circle(x, y, hw)]));
    }
    acc
}

fn circle(cx: f64, cy: f64, r: f64) -> Polygon {
    let steps = ARC_STEPS * 4;
    let mut coords: Vec<GeoCoord> = (0..steps)
        .map(|i| {
            let t = std::f64::consts::TAU * i as f64 / steps as f64;
            GeoCoord {
                x: cx + r * t.cos(),
                y: cy + r * t.sin(),
            }
        })
        .collect();
    coords.push(coords[0]);
    Polygon::new(LineString::new(coords), vec![])
}

/// Marks every grid cell whose centre falls inside `geom`.
fn mark_cells(grid: &Raster, geom: &Geometry, mask: &mut [bool]) {
    let Some(bb) = geom.bbox() else { return };
    let cols = grid.cols as isize;
    let rows = grid.rows as isize;
    let c0 = (((bb.min_x - grid.x_min) / grid.cell_size_x).floor() as isize).max(0);
    let c1 = (((bb.max_x - grid.x_min) / grid.cell_size_x).ceil() as isize).min(cols - 1);
    let r0 = (((grid.y_max() - bb.max_y) / grid.cell_size_y).floor() as isize).max(0);
    let r1 = (((grid.y_max() - bb.min_y) / grid.cell_size_y).ceil() as isize).min(rows - 1);
    for r in r0..=r1 {
        for c in c0..=c1 {
            let (x, y) = pixel_center(grid, r as usize, c as usize);
            if geometry_contains_point(geom, x, y) {
                mask[(r * cols + c) as usize] = true;
            }
        }
    }
}

fn pixel_center(r: &Raster, row: usize, col: usize) -> (f64, f64) {
    (
        r.x_min + (col as f64 + 0.5) * r.cell_size_x,
        r.y_max() - (row as f64 + 0.5) * r.cell_size_y,
    )
}

fn multipolygon_to_geometry(mp: &MultiPolygon) -> Geometry {
    if mp.0.len() == 1 {
        let (exterior, interiors) = polygon_to_rings(&mp.0[0]);
        Geometry::Polygon {
            exterior,
            interiors,
        }
    } else {
        Geometry::MultiPolygon(mp.0.iter().map(polygon_to_rings).collect())
    }
}

fn polygon_to_rings(poly: &Polygon) -> (Ring, Vec<Ring>) {
    (
        linestring_to_ring(poly.exterior()),
        poly.interiors().iter().map(linestring_to_ring).collect(),
    )
}

fn linestring_to_ring(ls: &LineString) -> Ring {
    let mut coords: Vec<Coord> = ls.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
    if coords.len() >= 2 && coords.first().map(|c| (c.x, c.y)) == coords.last().map(|c| (c.x, c.y))
    {
        coords.pop();
    }
    Ring::new(coords)
}

struct Node {
    cost: f64,
    idx: usize,
}
impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.cost == other.cost && self.idx == other.idx
    }
}
impl Eq for Node {}
impl Ord for Node {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .cost
            .total_cmp(&self.cost)
            .then_with(|| other.idx.cmp(&self.idx))
    }
}
impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// ── Params ──────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum NeighborOption {
    SpanningTree,
    AllPairs,
}

fn parse_neighbor_option(args: &ToolArgs) -> Result<NeighborOption, ToolError> {
    match args
        .get("neighbor_option")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        None | Some("") | Some("spanning_tree") => Ok(NeighborOption::SpanningTree),
        Some("all_pairs") => Ok(NeighborOption::AllPairs),
        Some(o) => Err(ToolError::Validation(format!(
            "'neighbor_option' must be 'spanning_tree' or 'all_pairs', got '{o}'"
        ))),
    }
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

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn square(cx: f64, cy: f64, h: f64) -> Geometry {
        Geometry::polygon(
            vec![
                Coord::xy(cx - h, cy - h),
                Coord::xy(cx + h, cy - h),
                Coord::xy(cx + h, cy + h),
                Coord::xy(cx - h, cy + h),
                Coord::xy(cx - h, cy - h),
            ],
            vec![],
        )
    }

    /// Regions at the given centres, each a 6-unit square.
    fn regions(centres: &[(f64, f64)]) -> String {
        let mut l = Layer::new("r")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        for (i, (x, y)) in centres.iter().enumerate() {
            l.add_feature(
                Some(square(*x, *y, 3.0)),
                &[("name", FieldValue::Text(format!("r{i}")))],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = OptimalCorridorConnectionsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    #[test]
    fn three_regions_yield_a_two_corridor_spanning_tree() {
        // A spanning tree over n regions has exactly n-1 edges — the point of
        // the default mode, versus the 3 corridors all_pairs would give.
        let (out, layer) = run(json!({
            "input": regions(&[(10.0, 50.0), (50.0, 50.0), (90.0, 50.0)]),
            "cell_size": 2.0, "corridor_width": 6.0
        }));
        assert_eq!(out.outputs["region_count"], json!(3));
        assert_eq!(out.outputs["corridor_count"], json!(2));
        assert_eq!(layer.features.len(), 2);
    }

    #[test]
    fn all_pairs_emits_every_region_pair() {
        let (out, _l) = run(json!({
            "input": regions(&[(10.0, 50.0), (50.0, 50.0), (90.0, 50.0)]),
            "cell_size": 2.0, "corridor_width": 6.0, "neighbor_option": "all_pairs"
        }));
        assert_eq!(out.outputs["corridor_count"], json!(3));
    }

    #[test]
    fn corridor_area_scales_with_the_requested_width() {
        // Doubling the width should roughly double the swept area — the
        // property that distinguishes this from a zero-width connectivity line.
        let input = regions(&[(10.0, 50.0), (90.0, 50.0)]);
        let area = |w: f64| -> f64 {
            let (out, _l) = run(json!({
                "input": input.clone(), "cell_size": 2.0, "corridor_width": w
            }));
            out.outputs["total_corridor_area"].as_f64().unwrap()
        };
        let (a1, a2) = (area(4.0), area(8.0));
        let ratio = a2 / a1;
        assert!(ratio > 1.7 && ratio < 2.3, "area ratio {ratio} not ~2");
    }

    #[test]
    fn corridors_are_polygons_with_positive_area() {
        let (out, layer) = run(json!({
            "input": regions(&[(10.0, 50.0), (90.0, 50.0)]),
            "cell_size": 2.0, "corridor_width": 6.0
        }));
        assert!(out.outputs["total_corridor_area"].as_f64().unwrap() > 0.0);
        let a = layer.schema.field_index("AREA").unwrap();
        for f in layer.iter() {
            assert!(f.attributes[a].as_f64().unwrap() > 0.0);
            assert!(matches!(
                f.geometry.as_ref().unwrap(),
                Geometry::Polygon { .. } | Geometry::MultiPolygon(_)
            ));
        }
    }

    #[test]
    fn a_barrier_lengthens_the_corridor() {
        // A wall between two regions must force a detour, not be ignored.
        let input = regions(&[(10.0, 50.0), (90.0, 50.0)]);
        let mut b = Layer::new("b")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        // Vertical wall at x=50 spanning y 20..80, leaving gaps top and bottom.
        b.add_feature(
            Some(Geometry::polygon(
                vec![
                    Coord::xy(48.0, 20.0),
                    Coord::xy(52.0, 20.0),
                    Coord::xy(52.0, 80.0),
                    Coord::xy(48.0, 80.0),
                    Coord::xy(48.0, 20.0),
                ],
                vec![],
            )),
            &[],
        )
        .unwrap();
        let id = wbvector::memory_store::put_vector(b);
        let barriers = wbvector::memory_store::make_vector_memory_path(&id);

        let len = |bar: Option<&str>| -> f64 {
            let mut a = json!({
                "input": input.clone(), "cell_size": 2.0, "corridor_width": 4.0
            });
            if let Some(p) = bar {
                a["barriers"] = json!(p);
            }
            let (out, _l) = run(a);
            out.outputs["total_centerline_length"].as_f64().unwrap()
        };
        let (open, walled) = (len(None), len(Some(&barriers)));
        assert!(
            walled > open * 1.2,
            "barrier did not force a detour: {open} -> {walled}"
        );
    }

    #[test]
    fn a_cost_surface_steers_the_corridor() {
        // Cheap band along y=20, expensive elsewhere: the corridor should
        // detour into the cheap band rather than run straight.
        let mut r = Raster::new(RasterConfig {
            cols: 50,
            rows: 50,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 2.0,
            cell_size_y: Some(2.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo::from_epsg(3857),
            metadata: Default::default(),
        });
        for row in 0..50 {
            for col in 0..50 {
                let y = 99.0 - row as f64 * 2.0;
                let cost = if (y - 20.0).abs() < 5.0 { 1.0 } else { 50.0 };
                r.set(0, row as isize, col as isize, cost).unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        let cost = wbraster::memory_store::make_raster_memory_path(&id);

        let (_o, layer) = run(json!({
            "input": regions(&[(10.0, 50.0), (90.0, 50.0)]),
            "cost_raster": cost, "corridor_width": 4.0, "output_lines": ""
        }));
        // The corridor polygon must reach down into the cheap band.
        let mut min_y = f64::INFINITY;
        for f in layer.iter() {
            for c in f.geometry.as_ref().unwrap().all_coords() {
                min_y = min_y.min(c.y);
            }
        }
        assert!(min_y < 30.0, "corridor stayed at y >= {min_y}, ignoring cost");
    }

    #[test]
    fn fewer_than_two_regions_is_rejected() {
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": regions(&[(10.0, 10.0)]) })).unwrap();
        assert!(OptimalCorridorConnectionsTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            OptimalCorridorConnectionsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "r.shp", "corridor_width": 0 })).is_err());
        assert!(bad(json!({ "input": "r.shp", "neighbor_option": "star" })).is_err());
        assert!(bad(json!({ "input": "r.shp" })).is_ok());
    }
}
