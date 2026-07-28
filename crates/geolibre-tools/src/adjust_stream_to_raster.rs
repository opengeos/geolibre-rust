//! GeoLibre tool: snap stream centerlines onto a DEM's true flow path.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Adjust Stream To Raster* (Spatial
//! Analyst). The catalog only conflated hydrography in one direction — it bent
//! the **DEM** to match the streams. `enforce_river_monotonicity` forces
//! elevations to decrease downstream along given lines, and the bundled
//! `burn_streams`, `burn_streams_at_roads`, `raise_walls` and
//! `topological_breach_burn` all carve or fence the raster.
//!
//! The inverse had no coverage anywhere. It is the case you hit whenever the
//! DEM is the *more* accurate dataset — a lidar-derived surface paired with
//! coarse or dated hydrography — where burning the streams in would gouge false
//! channels across the true valley floor. Adjusting the vector instead
//! preserves the surface and fixes the line.
//!
//! Method: derive D8 pointers and flow accumulation from the DEM, snap each
//! line's nodes to the **nearest** channelised cell within `snap_distance`
//! (nearest, not highest-accumulation — accumulation grows monotonically
//! downstream, so maximising it would drag every node toward the downstream
//! edge of its search window), then trace the pointer chain between
//! consecutive nodes. Where the pointer trace
//! fails to arrive (a tributary junction, a mismatch too large for the snap
//! window) the tool falls back to a least-cost route that prefers
//! high-accumulation cells, and records a split point there so the failure is
//! visible rather than silently smoothed over.
//!
//! **Scope note:** D8 assigns exactly one downstream neighbour per cell, so a
//! traced path never forks. Braided channels therefore cannot produce true
//! divergences here; `output_split_points` marks trace failures instead, which
//! is the honest signal this method can give.

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::common::{load_input_raster, raster_like_with_data, write_or_store_output};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// D8 offsets, indexed by pointer code 0..7.
const D8: [(isize, isize); 8] = [
    (-1, -1), (-1, 0), (-1, 1),
    (0, 1),   (1, 1),  (1, 0),
    (1, -1),  (0, -1),
];

pub struct AdjustStreamToRasterTool;

impl Tool for AdjustStreamToRasterTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "adjust_stream_to_raster",
            display_name: "Adjust Stream To Raster",
            summary: "Move stream polyline vertices onto the flow path implied by a DEM so vector hydrography and the elevation surface agree, by snapping nodes to high-accumulation cells and tracing the D8 pointer chain between them. The inverse of stream burning. Like ArcGIS Adjust Stream To Raster.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Stream polyline features to adjust.",
                    required: true,
                },
                ToolParamSpec {
                    name: "dem",
                    description: "Surface raster defining the true flow path.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output adjusted stream path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "snap_distance",
                    description: "Search radius in map units for snapping nodes to a high-accumulation cell (default: 5 cells).",
                    required: false,
                },
                ToolParamSpec {
                    name: "channel_threshold",
                    description: "Flow-accumulation value (in cells) at or above which a cell counts as channelised for snapping (default: 1% of the DEM's maximum accumulation).",
                    required: false,
                },
                ToolParamSpec {
                    name: "remove_disconnected",
                    description: "Drop lines that cannot be routed at all (default true). When false they are passed through unchanged and flagged.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_stream_raster",
                    description: "Optional raster marking the adjusted stream cells.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_flow_direction",
                    description: "Optional D8 pointer raster derived from the DEM.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_split_points",
                    description: "Optional points where the pointer trace failed and a least-cost fallback was used.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "dem")?;
        if let Some(d) = parse_optional_f64(args, "snap_distance")? {
            if d < 0.0 {
                return Err(ToolError::Validation(
                    "'snap_distance' must be non-negative".to_string(),
                ));
            }
        }
        if let Some(t) = parse_optional_f64(args, "channel_threshold")? {
            if t < 0.0 {
                return Err(ToolError::Validation(
                    "'channel_threshold' must be non-negative".to_string(),
                ));
            }
        }
        parse_optional_bool(args, "remove_disconnected")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let dem_path = require_str(args, "dem")?;
        let output = parse_optional_str(args, "output")?;
        let remove_disconnected = parse_optional_bool(args, "remove_disconnected")?.unwrap_or(true);

        let layer = load_input_layer(input)?;
        if layer.features.is_empty() {
            return Err(ToolError::Execution("input has no features".to_string()));
        }
        let dem = load_input_raster(dem_path)?;
        let rows = dem.rows;
        let cols = dem.cols;
        let snap_distance = parse_optional_f64(args, "snap_distance")?
            .unwrap_or_else(|| 5.0 * dem.cell_size_x.max(dem.cell_size_y));

        ctx.progress.info("deriving D8 pointers and accumulation");
        let (pointer, accum) = d8_pointer_and_accumulation(&dem);
        let max_accum = accum.iter().cloned().fold(0.0_f64, f64::max);
        let channel_threshold = parse_optional_f64(args, "channel_threshold")?
            .unwrap_or_else(|| (max_accum * 0.01).max(2.0));

        let mut out = Layer::new("adjusted_streams").with_geom_type(GeometryType::LineString);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for fd in layer.schema.fields() {
            out.add_field(fd.clone());
        }
        out.add_field(FieldDef::new("ORIG_FID", FieldType::Integer));
        out.add_field(FieldDef::new("SNAP_FROM", FieldType::Float));
        out.add_field(FieldDef::new("SNAP_TO", FieldType::Float));
        out.add_field(FieldDef::new("ROUTED", FieldType::Boolean));
        let names: Vec<String> = layer
            .schema
            .fields()
            .iter()
            .map(|f| f.name.clone())
            .collect();

        let mut splits = Layer::new("split_points").with_geom_type(GeometryType::Point);
        if let Some(epsg) = layer.crs_epsg() {
            splits = splits.with_crs_epsg(epsg);
        }
        splits.add_field(FieldDef::new("ORIG_FID", FieldType::Integer));
        splits.add_field(FieldDef::new("REASON", FieldType::Text));

        let mut stream_cells = vec![0.0_f64; rows * cols];
        let (mut adjusted, mut dropped, mut fallbacks) = (0usize, 0usize, 0usize);
        let mut total_shift = 0.0_f64;
        let mut shift_samples = 0usize;

        for (fid, feat) in layer.iter().enumerate() {
            let Some(g) = &feat.geometry else { continue };
            let paths = line_paths(g);
            if paths.is_empty() {
                continue;
            }
            for path in paths {
                // Snap every vertex that can be snapped; consecutive snapped
                // nodes bracket the segments we then trace.
                let mut nodes: Vec<(usize, usize, f64)> = Vec::new();
                for c in &path {
                    if let Some((r, cc, moved)) =
                        snap_to_flow(&dem, &accum, c.x, c.y, snap_distance, channel_threshold)
                    {
                        if nodes.last().map(|&(pr, pc, _)| (pr, pc)) != Some((r, cc)) {
                            nodes.push((r, cc, moved));
                        }
                    }
                }
                if nodes.len() < 2 {
                    dropped += 1;
                    if !remove_disconnected {
                        emit(&mut out, &names, feat, fid, g.clone(), 0.0, 0.0, false)?;
                    }
                    continue;
                }

                let mut cells: Vec<(usize, usize)> = vec![(nodes[0].0, nodes[0].1)];
                let mut routed = true;
                for w in nodes.windows(2) {
                    let (a, b) = ((w[0].0, w[0].1), (w[1].0, w[1].1));
                    match trace_d8(&pointer, rows, cols, a, b) {
                        Some(seg) => cells.extend(seg.into_iter().skip(1)),
                        None => {
                            // The pointer chain missed: fall back to a
                            // least-cost route that prefers channelised cells.
                            fallbacks += 1;
                            let (x, y) = pixel_center(&dem, a.0, a.1);
                            splits
                                .add_feature(
                                    Some(Geometry::point(x, y)),
                                    &[
                                        ("ORIG_FID", FieldValue::Integer(fid as i64)),
                                        (
                                            "REASON",
                                            FieldValue::Text("d8_trace_failed".to_string()),
                                        ),
                                    ],
                                )
                                .map_err(|e| {
                                    ToolError::Execution(format!("failed adding split point: {e}"))
                                })?;
                            match least_cost(&dem, &accum, a, b) {
                                Some(seg) => cells.extend(seg.into_iter().skip(1)),
                                None => {
                                    routed = false;
                                    break;
                                }
                            }
                        }
                    }
                }
                if !routed || cells.len() < 2 {
                    dropped += 1;
                    if !remove_disconnected {
                        emit(&mut out, &names, feat, fid, g.clone(), 0.0, 0.0, false)?;
                    }
                    continue;
                }

                for &(r, c) in &cells {
                    stream_cells[r * cols + c] = 1.0;
                }
                let coords: Vec<Coord> = cells
                    .iter()
                    .map(|&(r, c)| {
                        let (x, y) = pixel_center(&dem, r, c);
                        Coord::xyz(x, y, dem.get(0, r as isize, c as isize))
                    })
                    .collect();
                let (snap_from, snap_to) = (nodes[0].2, nodes[nodes.len() - 1].2);
                total_shift += snap_from + snap_to;
                shift_samples += 2;
                adjusted += 1;
                emit(
                    &mut out,
                    &names,
                    feat,
                    fid,
                    Geometry::line_string(coords),
                    snap_from,
                    snap_to,
                    true,
                )?;
            }
            ctx.progress
                .progress((fid as f64 + 1.0) / layer.features.len() as f64);
        }

        if adjusted == 0 {
            return Err(ToolError::Execution(
                "no stream line could be routed onto the DEM's flow path; \
                 check the CRS match and increase 'snap_distance'"
                    .to_string(),
            ));
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("adjusted_count".to_string(), json!(adjusted));
        outputs.insert("dropped_count".to_string(), json!(dropped));
        outputs.insert("fallback_count".to_string(), json!(fallbacks));
        outputs.insert("snap_distance".to_string(), json!(snap_distance));
        outputs.insert("channel_threshold".to_string(), json!(channel_threshold));
        if shift_samples > 0 {
            outputs.insert(
                "mean_node_shift".to_string(),
                json!(total_shift / shift_samples as f64),
            );
        }

        // Presence of the key requests the artifact; an empty value means
        // "produce it, but hand it back in memory rather than to a file".
        let requested = |k: &str| matches!(args.get(k), Some(v) if !v.is_null());
        if requested("output_stream_raster") {
            let p = parse_optional_str(args, "output_stream_raster")?;
            let r = raster_like_with_data(&dem, stream_cells, 0.0, DataType::F32)?;
            outputs.insert(
                "output_stream_raster".to_string(),
                json!(write_or_store_output(r, p)?),
            );
        }
        if requested("output_flow_direction") {
            let p = parse_optional_str(args, "output_flow_direction")?;
            let data: Vec<f64> = pointer
                .iter()
                .map(|&d| if d < 0 { -1.0 } else { d as f64 })
                .collect();
            let r = raster_like_with_data(&dem, data, -1.0, DataType::F32)?;
            outputs.insert(
                "output_flow_direction".to_string(),
                json!(write_or_store_output(r, p)?),
            );
        }
        if requested("output_split_points") {
            let p = parse_optional_str(args, "output_split_points")?;
            outputs.insert(
                "output_split_points".to_string(),
                json!(write_or_store_layer(splits, p)?),
            );
        }
        Ok(ToolRunResult { outputs })
    }
}

#[allow(clippy::too_many_arguments)]
fn emit(
    out: &mut Layer,
    names: &[String],
    feat: &wbvector::Feature,
    fid: usize,
    geom: Geometry,
    snap_from: f64,
    snap_to: f64,
    routed: bool,
) -> Result<(), ToolError> {
    let mut attrs: Vec<(&str, FieldValue)> = names
        .iter()
        .enumerate()
        .map(|(i, nm)| {
            (
                nm.as_str(),
                feat.attributes.get(i).cloned().unwrap_or(FieldValue::Null),
            )
        })
        .collect();
    attrs.push(("ORIG_FID", FieldValue::Integer(fid as i64)));
    attrs.push(("SNAP_FROM", FieldValue::Float(snap_from)));
    attrs.push(("SNAP_TO", FieldValue::Float(snap_to)));
    attrs.push(("ROUTED", FieldValue::Boolean(routed)));
    out.add_feature(Some(geom), &attrs)
        .map_err(|e| ToolError::Execution(format!("failed adding feature: {e}")))?;
    Ok(())
}

// ── Flow derivation ─────────────────────────────────────────────────────────

/// D8 pointer (index into `D8`, or -1 for none) plus flow accumulation in cells.
fn d8_pointer_and_accumulation(dem: &Raster) -> (Vec<i8>, Vec<f64>) {
    let rows = dem.rows;
    let cols = dem.cols;
    let nodata = dem.nodata;
    let mut pointer = vec![-1_i8; rows * cols];
    let mut order: Vec<(f64, usize)> = Vec::with_capacity(rows * cols);

    for r in 0..rows {
        for c in 0..cols {
            let z = dem.get(0, r as isize, c as isize);
            if z == nodata || !z.is_finite() {
                continue;
            }
            order.push((z, r * cols + c));
            // Steepest descent, normalised by the step length so diagonals
            // are not unfairly favoured.
            let mut best = (0.0_f64, -1_i8);
            for (k, (dr, dc)) in D8.iter().enumerate() {
                let (nr, nc) = (r as isize + dr, c as isize + dc);
                if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                    continue;
                }
                let nz = dem.get(0, nr, nc);
                if nz == nodata || !nz.is_finite() {
                    continue;
                }
                let run = if *dr != 0 && *dc != 0 {
                    dem.cell_size_x.hypot(dem.cell_size_y)
                } else {
                    dem.cell_size_x.min(dem.cell_size_y)
                };
                let slope = (z - nz) / run;
                if slope > best.0 {
                    best = (slope, k as i8);
                }
            }
            pointer[r * cols + c] = best.1;
        }
    }

    // Accumulate from high to low so every cell is settled before it drains.
    order.sort_by(|a, b| b.0.total_cmp(&a.0));
    let mut accum = vec![1.0_f64; rows * cols];
    for &(_, idx) in &order {
        let p = pointer[idx];
        if p < 0 {
            continue;
        }
        let (r, c) = (idx / cols, idx % cols);
        let (dr, dc) = D8[p as usize];
        let (nr, nc) = (r as isize + dr, c as isize + dc);
        if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
            continue;
        }
        let nidx = nr as usize * cols + nc as usize;
        accum[nidx] += accum[idx];
    }
    (pointer, accum)
}

/// Snaps a map coordinate to the **nearest** channelised cell within `radius`,
/// returning the cell and how far it moved.
///
/// Nearest-above-threshold, not highest-accumulation: accumulation grows
/// monotonically downstream, so maximising it over the window would drag every
/// node toward the window's downstream edge and shorten the line. The offset
/// this tool exists to fix is *across* the channel, so the along-channel
/// position must be preserved.
fn snap_to_flow(
    dem: &Raster,
    accum: &[f64],
    x: f64,
    y: f64,
    radius: f64,
    channel_threshold: f64,
) -> Option<(usize, usize, f64)> {
    let (c0, r0) = dem.world_to_pixel(x, y)?;
    let cols = dem.cols;
    let rows = dem.rows;
    let win_c = (radius / dem.cell_size_x).ceil() as isize;
    let win_r = (radius / dem.cell_size_y).ceil() as isize;
    // Nearest channel cell, and — as a fallback for windows containing no
    // channel at all — the most channelised cell in the window.
    let mut nearest: Option<(usize, usize, f64, f64)> = None; // r, c, accum, dist
    let mut strongest: Option<(usize, usize, f64, f64)> = None;
    for r in (r0 - win_r)..=(r0 + win_r) {
        for c in (c0 - win_c)..=(c0 + win_c) {
            if r < 0 || c < 0 || r >= rows as isize || c >= cols as isize {
                continue;
            }
            let (cx, cy) = pixel_center(dem, r as usize, c as usize);
            let d = (cx - x).hypot(cy - y);
            if d > radius {
                continue;
            }
            let a = accum[r as usize * cols + c as usize];
            if strongest.is_none_or(|(_, _, ba, bd)| a > ba || (a == ba && d < bd)) {
                strongest = Some((r as usize, c as usize, a, d));
            }
            if a >= channel_threshold
                && nearest.is_none_or(|(_, _, ba, bd)| d < bd || (d == bd && a > ba))
            {
                nearest = Some((r as usize, c as usize, a, d));
            }
        }
    }
    nearest.or(strongest).map(|(r, c, _, d)| (r, c, d))
}

/// Follows the D8 chain from `a`, returning the cell path if it reaches `b`.
fn trace_d8(
    pointer: &[i8],
    rows: usize,
    cols: usize,
    a: (usize, usize),
    b: (usize, usize),
) -> Option<Vec<(usize, usize)>> {
    let budget = rows * cols;
    let mut path = vec![a];
    let (mut r, mut c) = a;
    for _ in 0..budget {
        if (r, c) == b {
            return Some(path);
        }
        let p = pointer[r * cols + c];
        if p < 0 {
            return None;
        }
        let (dr, dc) = D8[p as usize];
        let (nr, nc) = (r as isize + dr, c as isize + dc);
        if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
            return None;
        }
        r = nr as usize;
        c = nc as usize;
        path.push((r, c));
    }
    None
}

/// Least-cost route preferring high-accumulation (channelised) cells.
///
/// The search is confined to the bounding box of `a` and `b` inflated by a
/// margin, and the buffers are sized to that window. This runs once per failed
/// D8 trace, so allocating and sweeping the whole raster each time would
/// multiply the cost of a mismatched network by the number of fallbacks.
fn least_cost(
    dem: &Raster,
    accum: &[f64],
    a: (usize, usize),
    b: (usize, usize),
) -> Option<Vec<(usize, usize)>> {
    let rows = dem.rows;
    let cols = dem.cols;
    let nodata = dem.nodata;
    let diag = dem.cell_size_x.hypot(dem.cell_size_y);

    // Window: the a/b bbox, inflated by its own size (min 8 cells) so the route
    // has room to bend around obstacles without covering the whole DEM.
    let (r_lo, r_hi) = (a.0.min(b.0), a.0.max(b.0));
    let (c_lo, c_hi) = (a.1.min(b.1), a.1.max(b.1));
    let pad_r = (r_hi - r_lo).max(8);
    let pad_c = (c_hi - c_lo).max(8);
    let wr0 = r_lo.saturating_sub(pad_r);
    let wc0 = c_lo.saturating_sub(pad_c);
    let wr1 = (r_hi + pad_r).min(rows - 1);
    let wc1 = (c_hi + pad_c).min(cols - 1);
    let (wrows, wcols) = (wr1 - wr0 + 1, wc1 - wc0 + 1);
    let n = wrows * wcols;
    // Window-local <-> raster index.
    let local = |r: usize, c: usize| (r - wr0) * wcols + (c - wc0);
    let global = |i: usize| (wr0 + i / wcols, wc0 + i % wcols);

    let mut dist = vec![f64::INFINITY; n];
    let mut back = vec![usize::MAX; n];
    let mut done = vec![false; n];
    let mut heap: BinaryHeap<Node> = BinaryHeap::new();
    let start = local(a.0, a.1);
    let goal = local(b.0, b.1);
    dist[start] = 0.0;
    heap.push(Node {
        cost: 0.0,
        idx: start,
    });

    while let Some(Node { cost, idx }) = heap.pop() {
        if done[idx] {
            continue;
        }
        done[idx] = true;
        if idx == goal {
            let mut path = vec![global(idx)];
            let mut cur = idx;
            while back[cur] != usize::MAX {
                cur = back[cur];
                path.push(global(cur));
            }
            path.reverse();
            return Some(path);
        }
        let (r, c) = global(idx);
        for (dr, dc) in D8 {
            let (nr, nc) = (r as isize + dr, c as isize + dc);
            if nr < wr0 as isize || nc < wc0 as isize || nr > wr1 as isize || nc > wc1 as isize {
                continue;
            }
            let (nr, nc) = (nr as usize, nc as usize);
            let nidx = local(nr, nc);
            if done[nidx] {
                continue;
            }
            let nz = dem.get(0, nr as isize, nc as isize);
            if nz == nodata || !nz.is_finite() {
                continue;
            }
            // Per-axis step length, so anisotropic cells are costed correctly.
            let step = match (dr != 0, dc != 0) {
                (true, true) => diag,
                (true, false) => dem.cell_size_y,
                (false, true) => dem.cell_size_x,
                (false, false) => continue,
            };
            // Cheap where flow concentrates, expensive on hillslopes.
            let nd = cost + step / (1.0 + accum[nr * cols + nc]);
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
    None
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

fn pixel_center(r: &Raster, row: usize, col: usize) -> (f64, f64) {
    (
        r.x_min + (col as f64 + 0.5) * r.cell_size_x,
        r.y_max() - (row as f64 + 0.5) * r.cell_size_y,
    )
}

fn line_paths(g: &Geometry) -> Vec<Vec<Coord>> {
    match g {
        Geometry::LineString(cs) if cs.len() >= 2 => vec![cs.clone()],
        Geometry::MultiLineString(ls) => ls.iter().filter(|l| l.len() >= 2).cloned().collect(),
        Geometry::GeometryCollection(gs) => gs.iter().flat_map(line_paths).collect(),
        _ => Vec::new(),
    }
}

// ── Params ──────────────────────────────────────────────────────────────────

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
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
    use wbraster::RasterConfig;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// 40x40 DEM over (0,0)-(40,40), cell 1.
    fn dem(f: impl Fn(f64, f64) -> f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols: 40,
            rows: 40,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: Default::default(),
            metadata: Default::default(),
        });
        for row in 0..40 {
            for col in 0..40 {
                let x = 0.5 + col as f64;
                let y = 39.5 - row as f64;
                r.set(0, row as isize, col as isize, f(x, y)).unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    /// A V-shaped valley whose thalweg runs along y = 20, falling to the east.
    fn valley() -> String {
        dem(|x, y| 100.0 - x * 0.5 + (y - 20.0).abs() * 2.0)
    }

    /// A straight line offset from the true thalweg by `dy`.
    fn stream(dy: f64) -> String {
        let mut l = Layer::new("s")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.add_feature(
            Some(Geometry::line_string(vec![
                Coord::xy(5.5, 20.5 + dy),
                Coord::xy(20.5, 20.5 + dy),
                Coord::xy(34.5, 20.5 + dy),
            ])),
            &[("name", FieldValue::Text("creek".into()))],
        )
        .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = AdjustStreamToRasterTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// Mean absolute offset of a line's vertices from the true thalweg y = 20.
    fn mean_offset(layer: &Layer) -> f64 {
        let mut sum = 0.0;
        let mut n = 0;
        for f in layer.iter() {
            for c in f.geometry.as_ref().unwrap().all_coords() {
                sum += (c.y - 20.0).abs();
                n += 1;
            }
        }
        sum / n as f64
    }

    #[test]
    fn an_offset_stream_is_pulled_onto_the_thalweg() {
        // The line sits 3 units off the valley floor; adjusting must move it on.
        let (out, layer) = run(json!({
            "input": stream(3.0), "dem": valley(), "snap_distance": 6.0
        }));
        assert_eq!(out.outputs["adjusted_count"], json!(1));
        let off = mean_offset(&layer);
        assert!(off < 1.5, "adjusted line still {off} off the thalweg");
    }

    #[test]
    fn a_stream_already_on_the_thalweg_barely_moves() {
        let (out, layer) = run(json!({
            "input": stream(0.0), "dem": valley(), "snap_distance": 6.0
        }));
        assert!(mean_offset(&layer) < 1.0);
        assert!(out.outputs["mean_node_shift"].as_f64().unwrap() < 2.0);
    }

    #[test]
    fn attributes_survive_the_adjustment() {
        let (_o, layer) = run(json!({
            "input": stream(3.0), "dem": valley(), "snap_distance": 6.0
        }));
        let n = layer.schema.field_index("name").unwrap();
        assert_eq!(layer.features[0].attributes[n].as_str(), Some("creek"));
        let routed = layer.schema.field_index("ROUTED").unwrap();
        assert_eq!(layer.features[0].attributes[routed].as_bool(), Some(true));
    }

    #[test]
    fn output_vertices_carry_dem_elevations() {
        let (_o, layer) = run(json!({
            "input": stream(2.0), "dem": valley(), "snap_distance": 6.0
        }));
        let coords = layer.features[0].geometry.as_ref().unwrap().all_coords();
        assert!(coords.iter().all(|c| c.z.is_some()));
        // Elevation must fall downstream along the adjusted line.
        let zs: Vec<f64> = coords.iter().map(|c| c.z.unwrap()).collect();
        assert!(zs.first().unwrap() > zs.last().unwrap());
    }

    #[test]
    fn optional_rasters_are_written_on_request() {
        let args: ToolArgs = serde_json::from_value(json!({
            "input": stream(2.0), "dem": valley(), "snap_distance": 6.0,
            "output_stream_raster": "", "output_flow_direction": ""
        }))
        .unwrap();
        let out = AdjustStreamToRasterTool.run(&args, &ctx()).unwrap();
        let sr = crate::common::load_input_raster(
            out.outputs["output_stream_raster"].as_str().unwrap(),
        )
        .unwrap();
        let marked = (0..sr.rows)
            .flat_map(|r| (0..sr.cols).map(move |c| (r, c)))
            .filter(|&(r, c)| sr.get(0, r as isize, c as isize) > 0.0)
            .count();
        assert!(marked > 0, "no stream cells marked");
        assert!(out.outputs.contains_key("output_flow_direction"));
    }

    #[test]
    fn split_points_layer_is_emitted_when_requested() {
        let args: ToolArgs = serde_json::from_value(json!({
            "input": stream(3.0), "dem": valley(), "snap_distance": 6.0,
            "output_split_points": ""
        }))
        .unwrap();
        let out = AdjustStreamToRasterTool.run(&args, &ctx()).unwrap();
        // May legitimately be empty; the contract is that the layer exists.
        let sp = load_input_layer(out.outputs["output_split_points"].as_str().unwrap()).unwrap();
        assert!(sp.schema.field_index("REASON").is_some());
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            AdjustStreamToRasterTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "s.shp" })).is_err());
        assert!(bad(json!({ "input": "s.shp", "dem": "z.tif", "snap_distance": -1 })).is_err());
        assert!(bad(json!({ "input": "s.shp", "dem": "z.tif" })).is_ok());
    }
}
