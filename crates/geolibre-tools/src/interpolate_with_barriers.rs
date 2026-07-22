//! GeoLibre tool: barrier-aware point interpolation via cost (geodesic) distance.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Kernel Interpolation With Barriers*
//! and *Spline With Barriers* (Geostatistical / Spatial Analyst). The entire
//! bundled interpolation suite — `idw_interpolation`, the kriging + variogram
//! tools, `natural_neighbour_interpolation`, `thin_plate_spline`,
//! `tin_interpolation` — is barrier-blind: a sample on one side of an estuary,
//! ridge, fault, or wall contaminates estimates on the other side because they
//! all weight by straight-line Euclidean distance.
//!
//! This tool respects **absolute barriers**. Barriers (polylines or polygon
//! boundaries) are rasterized into an impassable mask over the interpolation
//! grid; influence then travels around them, not across. Concretely, the
//! straight-line sample→cell distance is replaced by a **cost / geodesic
//! distance** computed with a multi-source-style 8-connected Dijkstra over the
//! free-space grid (the same least-cost engine as `path_distance` /
//! `cost_connectivity` / `corridor`). On those geodesic distances the tool then
//! applies either:
//!
//! * `idw` (default) — inverse-distance weighting, `w = 1 / d^power`; or
//! * `local_polynomial` — a first-order (planar) weighted least-squares fit with
//!   a Gaussian kernel of the geodesic distance, `w = exp(-(d/bandwidth)^2)`.
//!
//! A cell not reachable from any sample within `radius` (geodesic) becomes
//! no-data — so barrier-enclosed pockets with no interior samples are left
//! empty rather than bleeding values through the wall. Distances are metres for
//! a geographic CRS (degrees are scaled at the grid's mean latitude) and CRS
//! units otherwise. Output is an F32 raster.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
use wbvector::{Coord, Geometry};

use crate::common::write_or_store_output;
use crate::vector_common::{load_input_layer, parse_optional_str};

const OUT_NODATA: f64 = -9999.0;
/// Hard cap on grid dimensions to keep a single run tractable.
const MAX_DIM: usize = 4000;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Idw,
    LocalPolynomial,
}

pub struct InterpolateWithBarriersTool;

impl Tool for InterpolateWithBarriersTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "interpolate_with_barriers",
            display_name: "Interpolate With Barriers",
            summary: "Interpolate point measurements to a raster while respecting absolute barriers (shorelines, ridges, faults, walls): influence travels around barriers via cost/geodesic distance, not across them — IDW or a first-order local-polynomial kernel, like ArcGIS Kernel Interpolation With Barriers / Spline With Barriers. The bundled idw/kriging/natural-neighbour/spline tools are all barrier-blind.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer of measurements.",
                    required: true,
                },
                ToolParamSpec {
                    name: "field",
                    description: "Numeric field on the points holding the value to interpolate.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output interpolated raster (GeoTIFF). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "barriers",
                    description: "Optional barrier vector layer (polylines or polygons). Rasterized as impassable cells so influence must route around them. Omit for a barrier-free control.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'idw' (default) or 'local_polynomial' (first-order weighted planar fit with a Gaussian geodesic-distance kernel).",
                    required: false,
                },
                ToolParamSpec {
                    name: "power",
                    description: "IDW distance-decay exponent (default 2).",
                    required: false,
                },
                ToolParamSpec {
                    name: "bandwidth",
                    description: "Gaussian kernel bandwidth for local_polynomial (distance units). Default: the search radius, or the grid extent otherwise.",
                    required: false,
                },
                ToolParamSpec {
                    name: "radius",
                    description: "Search radius as a geodesic (cost) distance: samples beyond it do not influence a cell, and cells unreachable within it become no-data. Metres for a geographic CRS, CRS units otherwise. Default: unbounded.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Output cell size in CRS units (degrees for a geographic CRS). Default: sized so the longer axis spans ~256 cells.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "field")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let field = require_str(args, "field")?.to_string();
        let output = parse_optional_output(args, "output")?;
        let barriers_path = parse_optional_str(args, "barriers")?;
        let prm = parse_params(args)?;

        // ── Load samples ─────────────────────────────────────────────────────
        let layer = load_input_layer(input)?;
        let field_idx = layer
            .schema
            .field_index(&field)
            .ok_or_else(|| ToolError::Validation(format!("field '{field}' not found")))?;

        let mut sx: Vec<f64> = Vec::new();
        let mut sy: Vec<f64> = Vec::new();
        let mut sv: Vec<f64> = Vec::new();
        for feat in layer.iter() {
            let Some(geom) = feat.geometry.as_ref() else {
                continue;
            };
            let Some((x, y)) = point_xy(geom) else {
                continue;
            };
            let Some(v) = feat.attributes.get(field_idx).and_then(|f| f.as_f64()) else {
                continue;
            };
            if !v.is_finite() {
                continue;
            }
            sx.push(x);
            sy.push(y);
            sv.push(v);
        }
        if sx.is_empty() {
            return Err(ToolError::Execution(format!(
                "no point features with a finite '{field}' value"
            )));
        }
        let n_samples = sx.len();

        // ── Build the interpolation grid over the sample bounding box ────────
        let (mut x_min, mut y_min, mut x_max, mut y_max) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for i in 0..n_samples {
            x_min = x_min.min(sx[i]);
            y_min = y_min.min(sy[i]);
            x_max = x_max.max(sx[i]);
            y_max = y_max.max(sy[i]);
        }
        let width = (x_max - x_min).max(0.0);
        let height = (y_max - y_min).max(0.0);
        let span = width.max(height);
        if span <= 0.0 {
            return Err(ToolError::Execution(
                "all sample points are coincident; cannot build a grid".to_string(),
            ));
        }
        let cell = match prm.cell_size {
            Some(c) => c,
            None => span / 256.0,
        };
        if !(cell.is_finite() && cell > 0.0) {
            return Err(ToolError::Validation(
                "'cell_size' must be a positive number".to_string(),
            ));
        }
        // One-cell margin so edge samples sit inside the grid.
        let ox = x_min - cell;
        let oy_top = y_max + cell; // world y at the top edge (row 0)
        let cols = (((width + 2.0 * cell) / cell).ceil() as usize).max(1);
        let rows = (((height + 2.0 * cell) / cell).ceil() as usize).max(1);
        if cols > MAX_DIM || rows > MAX_DIM {
            return Err(ToolError::Validation(format!(
                "grid {rows}x{cols} exceeds the {MAX_DIM} cap; increase 'cell_size'"
            )));
        }
        let n_cells = rows * cols;

        // ── Distance metric (metres for geographic, CRS units otherwise) ─────
        let geographic = layer.crs_epsg().map(|e| e == 4326).unwrap_or(true);
        let (step_x, step_y) = if geographic {
            let mean_lat = 0.5 * (y_min + y_max);
            (
                cell * 111_320.0 * mean_lat.to_radians().cos().abs(),
                cell * 110_540.0,
            )
        } else {
            (cell, cell)
        };
        let step_diag = step_x.hypot(step_y);
        let d_min = 0.5 * step_x.min(step_y); // IDW clamp: avoid a divide-by-zero singularity

        // ── Rasterize barriers into an impassable mask ───────────────────────
        let mut passable = vec![true; n_cells];
        let mut barrier_cells = 0usize;
        if let Some(bpath) = barriers_path {
            let blayer = load_input_layer(bpath)?;
            for feat in blayer.iter() {
                if let Some(geom) = feat.geometry.as_ref() {
                    rasterize_barrier(geom, ox, oy_top, cell, rows, cols, &mut passable);
                }
            }
            barrier_cells = passable.iter().filter(|&&p| !p).count();
        }

        // Snap samples to cells; force source cells passable so a sample landing
        // on a barrier still seeds the field.
        let mut source_cell: Vec<Option<usize>> = Vec::with_capacity(n_samples);
        for i in 0..n_samples {
            let cidx = world_to_cell(sx[i], sy[i], ox, oy_top, cell, rows, cols);
            if let Some(c) = cidx {
                passable[c] = true;
            }
            source_cell.push(cidx);
        }

        ctx.progress.info(&format!(
            "{n_samples} sample(s) -> {rows}x{cols} grid, {barrier_cells} barrier cell(s), method {}",
            prm.method_label()
        ));

        // ── Accumulators ─────────────────────────────────────────────────────
        // IDW: numerator/denominator. local_polynomial: nine weighted moments of
        // the planar system z = c0 + c1*X + c2*Y (X,Y centred & scaled for
        // conditioning), plus a per-cell contributing-sample count.
        let is_poly = prm.method == Method::LocalPolynomial;
        let mut num = vec![0.0_f64; n_cells];
        let mut den = vec![0.0_f64; n_cells];
        // local_polynomial only: nine weighted moments plus the value range of
        // each cell's contributors, so the planar fit can be clamped and never
        // extrapolate wildly beyond nearby data.
        let (mut moments, mut lo, mut hi) = if is_poly {
            (
                vec![[0.0_f64; 9]; n_cells],
                vec![f64::INFINITY; n_cells],
                vec![f64::NEG_INFINITY; n_cells],
            )
        } else {
            (Vec::new(), Vec::new(), Vec::new())
        };
        let mut count = vec![0u32; n_cells];

        let bandwidth = prm
            .bandwidth
            .unwrap_or_else(|| prm.radius.unwrap_or_else(|| span.max(cell) * step_x / cell));
        let cx0 = ox + 0.5 * (cols as f64) * cell; // grid centre (world coords)
        let cy0 = oy_top - 0.5 * (rows as f64) * cell;
        let scale = cell; // keep X,Y ~ O(cols)

        let neigh = neighbours(step_x, step_y, step_diag);
        let radius = prm.radius.unwrap_or(f64::INFINITY);

        // ── Per-source bounded Dijkstra (geodesic distance around barriers) ──
        let mut dist = vec![f64::INFINITY; n_cells];
        let mut settled = vec![false; n_cells];
        let mut heap: BinaryHeap<Node> = BinaryHeap::new();
        let mut touched: Vec<usize> = Vec::new();

        for (s, &src_opt) in source_cell.iter().enumerate() {
            let Some(src) = src_opt else { continue };
            heap.clear();
            touched.clear();
            dist[src] = 0.0;
            heap.push(Node {
                cost: 0.0,
                idx: src,
            });
            touched.push(src);

            while let Some(Node { cost: acc, idx }) = heap.pop() {
                if settled[idx] {
                    continue;
                }
                settled[idx] = true;

                // Deposit this sample's contribution to the settled cell.
                deposit(
                    idx,
                    acc,
                    s,
                    &sx,
                    &sy,
                    &sv,
                    is_poly,
                    prm.power,
                    bandwidth,
                    d_min,
                    cx0,
                    cy0,
                    scale,
                    &mut num,
                    &mut den,
                    &mut moments,
                    &mut lo,
                    &mut hi,
                    &mut count,
                );

                let r = (idx / cols) as isize;
                let c = (idx % cols) as isize;
                for &Neighbour {
                    dr,
                    dc,
                    cost,
                    ortho_a,
                    ortho_b,
                } in &neigh
                {
                    let nr = r + dr;
                    let nc = c + dc;
                    if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                        continue;
                    }
                    let nidx = nr as usize * cols + nc as usize;
                    if !passable[nidx] {
                        continue;
                    }
                    // Anti-leak: forbid a diagonal step that would slip through
                    // the corner of a thin barrier.
                    if let (Some(a), Some(b)) = (ortho_a, ortho_b) {
                        let ai = (r + a.0) as usize * cols + (c + a.1) as usize;
                        let bi = (r + b.0) as usize * cols + (c + b.1) as usize;
                        if !passable[ai] || !passable[bi] {
                            continue;
                        }
                    }
                    let nd = acc + cost;
                    if nd > radius {
                        continue;
                    }
                    if nd < dist[nidx] {
                        dist[nidx] = nd;
                        heap.push(Node {
                            cost: nd,
                            idx: nidx,
                        });
                        touched.push(nidx);
                    }
                }
            }

            // Reset only the cells this source touched.
            for &t in &touched {
                dist[t] = f64::INFINITY;
                settled[t] = false;
            }
            if s % 64 == 0 {
                ctx.progress.progress((s as f64 + 1.0) / n_samples as f64);
            }
        }

        // ── Finalize each cell ───────────────────────────────────────────────
        let cxc: Vec<f64> = (0..cols).map(|c| ox + (c as f64 + 0.5) * cell).collect();
        let cyc: Vec<f64> = (0..rows)
            .map(|r| oy_top - (r as f64 + 0.5) * cell)
            .collect();
        let mut data = vec![OUT_NODATA; n_cells];
        let mut filled = 0usize;
        for (r, &cy) in cyc.iter().enumerate() {
            for (c, &cx) in cxc.iter().enumerate() {
                let idx = r * cols + c;
                if count[idx] == 0 {
                    continue;
                }
                let value = if is_poly {
                    let m = &moments[idx];
                    let xc = (cx - cx0) / scale;
                    let yc = (cy - cy0) / scale;
                    finalize_poly(m, count[idx], xc, yc, lo[idx], hi[idx])
                } else if den[idx] > 0.0 {
                    num[idx] / den[idx]
                } else {
                    continue;
                };
                if value.is_finite() {
                    data[idx] = value;
                    filled += 1;
                }
            }
        }
        if filled == 0 {
            return Err(ToolError::Execution(
                "no cell was reachable from any sample within the search radius".to_string(),
            ));
        }

        // ── Emit raster ──────────────────────────────────────────────────────
        let mut out = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: ox,
            y_min: oy_top - rows as f64 * cell,
            cell_size: cell,
            cell_size_y: Some(cell),
            nodata: OUT_NODATA,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: layer.crs_epsg(),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for r in 0..rows {
            for c in 0..cols {
                out.set(0, r as isize, c as isize, data[r * cols + c])
                    .map_err(|e| ToolError::Execution(format!("failed writing cell: {e}")))?;
            }
        }
        let out_path = write_or_store_output(out, output)?;

        let (mut vmin, mut vmax) = (f64::INFINITY, f64::NEG_INFINITY);
        for &v in &data {
            if v != OUT_NODATA {
                vmin = vmin.min(v);
                vmax = vmax.max(v);
            }
        }

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("samples".to_string(), json!(n_samples));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        outputs.insert("cell_size".to_string(), json!(cell));
        outputs.insert("barrier_cells".to_string(), json!(barrier_cells));
        outputs.insert("filled_cells".to_string(), json!(filled));
        outputs.insert("nodata_cells".to_string(), json!(n_cells - filled));
        outputs.insert("value_min".to_string(), json!(vmin));
        outputs.insert("value_max".to_string(), json!(vmax));
        Ok(ToolRunResult { outputs })
    }
}

// ── Contribution / finalize ─────────────────────────────────────────────────

#[allow(clippy::too_many_arguments)]
fn deposit(
    idx: usize,
    d: f64,
    s: usize,
    sx: &[f64],
    sy: &[f64],
    sv: &[f64],
    is_poly: bool,
    power: f64,
    bandwidth: f64,
    d_min: f64,
    cx0: f64,
    cy0: f64,
    scale: f64,
    num: &mut [f64],
    den: &mut [f64],
    moments: &mut [[f64; 9]],
    lo: &mut [f64],
    hi: &mut [f64],
    count: &mut [u32],
) {
    count[idx] += 1;
    if is_poly {
        let z = sv[s];
        lo[idx] = lo[idx].min(z);
        hi[idx] = hi[idx].max(z);
        let w = (-(d / bandwidth).powi(2)).exp();
        if w <= 0.0 {
            return;
        }
        let x = (sx[s] - cx0) / scale;
        let y = (sy[s] - cy0) / scale;
        let m = &mut moments[idx];
        m[0] += w;
        m[1] += w * x;
        m[2] += w * y;
        m[3] += w * x * x;
        m[4] += w * x * y;
        m[5] += w * y * y;
        m[6] += w * z;
        m[7] += w * x * z;
        m[8] += w * y * z;
    } else {
        let de = d.max(d_min);
        let w = 1.0 / de.powf(power);
        num[idx] += w * sv[s];
        den[idx] += w;
    }
}

/// Solves the weighted planar least-squares system and evaluates it at the cell
/// centre, clamped to the value range of the contributing samples so the fit
/// never extrapolates wildly. Falls back to the weighted mean when
/// under-determined or singular.
fn finalize_poly(m: &[f64; 9], count: u32, xc: f64, yc: f64, lo: f64, hi: f64) -> f64 {
    let mean = if m[0] > 0.0 { m[6] / m[0] } else { OUT_NODATA };
    if count < 3 {
        return mean;
    }
    // Symmetric 3x3 normal-equation matrix A c = b.
    let a = [[m[0], m[1], m[2]], [m[1], m[3], m[4]], [m[2], m[4], m[5]]];
    let b = [m[6], m[7], m[8]];
    match solve3(a, b) {
        Some(c) => (c[0] + c[1] * xc + c[2] * yc).clamp(lo, hi),
        None => mean,
    }
}

/// Solves a 3x3 linear system by cofactor expansion; `None` if near-singular.
fn solve3(a: [[f64; 3]; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if det.abs() < 1e-12 {
        return None;
    }
    let mut out = [0.0; 3];
    for i in 0..3 {
        let mut ai = a;
        for r in 0..3 {
            ai[r][i] = b[r];
        }
        let d = ai[0][0] * (ai[1][1] * ai[2][2] - ai[1][2] * ai[2][1])
            - ai[0][1] * (ai[1][0] * ai[2][2] - ai[1][2] * ai[2][0])
            + ai[0][2] * (ai[1][0] * ai[2][1] - ai[1][1] * ai[2][0]);
        out[i] = d / det;
    }
    Some(out)
}

// ── Grid / geometry helpers ─────────────────────────────────────────────────

struct Neighbour {
    dr: isize,
    dc: isize,
    cost: f64,
    ortho_a: Option<(isize, isize)>,
    ortho_b: Option<(isize, isize)>,
}

fn neighbours(step_x: f64, step_y: f64, step_diag: f64) -> [Neighbour; 8] {
    [
        Neighbour {
            dr: -1,
            dc: 0,
            cost: step_y,
            ortho_a: None,
            ortho_b: None,
        },
        Neighbour {
            dr: 1,
            dc: 0,
            cost: step_y,
            ortho_a: None,
            ortho_b: None,
        },
        Neighbour {
            dr: 0,
            dc: -1,
            cost: step_x,
            ortho_a: None,
            ortho_b: None,
        },
        Neighbour {
            dr: 0,
            dc: 1,
            cost: step_x,
            ortho_a: None,
            ortho_b: None,
        },
        Neighbour {
            dr: -1,
            dc: -1,
            cost: step_diag,
            ortho_a: Some((-1, 0)),
            ortho_b: Some((0, -1)),
        },
        Neighbour {
            dr: -1,
            dc: 1,
            cost: step_diag,
            ortho_a: Some((-1, 0)),
            ortho_b: Some((0, 1)),
        },
        Neighbour {
            dr: 1,
            dc: -1,
            cost: step_diag,
            ortho_a: Some((1, 0)),
            ortho_b: Some((0, -1)),
        },
        Neighbour {
            dr: 1,
            dc: 1,
            cost: step_diag,
            ortho_a: Some((1, 0)),
            ortho_b: Some((0, 1)),
        },
    ]
}

fn world_to_cell(
    x: f64,
    y: f64,
    ox: f64,
    oy_top: f64,
    cell: f64,
    rows: usize,
    cols: usize,
) -> Option<usize> {
    let col = ((x - ox) / cell).floor();
    let row = ((oy_top - y) / cell).floor();
    if col < 0.0 || row < 0.0 || col >= cols as f64 || row >= rows as f64 {
        return None;
    }
    Some(row as usize * cols + col as usize)
}

/// Marks every cell touched by a barrier geometry impassable. Lines mark cells
/// along each segment; polygons mark their ring boundaries (the shoreline/wall),
/// leaving interiors traversable.
fn rasterize_barrier(
    geom: &Geometry,
    ox: f64,
    oy_top: f64,
    cell: f64,
    rows: usize,
    cols: usize,
    passable: &mut [bool],
) {
    let mut mark_seg = |a: &Coord, b: &Coord, passable: &mut [bool]| {
        // Step along the segment in half-cell increments.
        let dx = b.x - a.x;
        let dy = b.y - a.y;
        let seg = (dx.hypot(dy) / (0.5 * cell)).ceil().max(1.0) as usize;
        for k in 0..=seg {
            let t = k as f64 / seg as f64;
            let px = a.x + t * dx;
            let py = a.y + t * dy;
            if let Some(idx) = world_to_cell(px, py, ox, oy_top, cell, rows, cols) {
                passable[idx] = false;
            }
        }
    };
    let mark_ring = |coords: &[Coord],
                     passable: &mut [bool],
                     mark_seg: &mut dyn FnMut(&Coord, &Coord, &mut [bool])| {
        for w in coords.windows(2) {
            mark_seg(&w[0], &w[1], passable);
        }
    };
    match geom {
        Geometry::LineString(cs) => mark_ring(cs, passable, &mut mark_seg),
        Geometry::MultiLineString(lines) => {
            for l in lines {
                mark_ring(l, passable, &mut mark_seg);
            }
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            mark_ring(exterior.coords(), passable, &mut mark_seg);
            for ring in interiors {
                mark_ring(ring.coords(), passable, &mut mark_seg);
            }
        }
        Geometry::MultiPolygon(polys) => {
            for (ext, ints) in polys {
                mark_ring(ext.coords(), passable, &mut mark_seg);
                for ring in ints {
                    mark_ring(ring.coords(), passable, &mut mark_seg);
                }
            }
        }
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                rasterize_barrier(g, ox, oy_top, cell, rows, cols, passable);
            }
        }
        _ => {}
    }
}

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        _ => None,
    }
}

// ── Dijkstra node (min-heap via reversed Ord) ───────────────────────────────

#[derive(Clone, Copy)]
struct Node {
    cost: f64,
    idx: usize,
}
impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.cost == other.cost
    }
}
impl Eq for Node {}
impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Node {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .cost
            .partial_cmp(&self.cost)
            .unwrap_or(Ordering::Equal)
            .then(self.idx.cmp(&other.idx))
    }
}

// ── Parameters ───────────────────────────────────────────────────────────────

struct Params {
    method: Method,
    power: f64,
    bandwidth: Option<f64>,
    radius: Option<f64>,
    cell_size: Option<f64>,
}

impl Params {
    fn method_label(&self) -> &'static str {
        match self.method {
            Method::Idw => "idw",
            Method::LocalPolynomial => "local_polynomial",
        }
    }
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let method = match args.get("method").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("idw") => Method::Idw,
        Some("local_polynomial") => Method::LocalPolynomial,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'method' must be 'idw' or 'local_polynomial', got '{o}'"
            )))
        }
    };
    let power = opt_pos(args, "power")?.unwrap_or(2.0);
    let bandwidth = opt_pos(args, "bandwidth")?;
    let radius = opt_pos(args, "radius")?;
    let cell_size = opt_pos(args, "cell_size")?;
    Ok(Params {
        method,
        power,
        bandwidth,
        radius,
        cell_size,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn parse_optional_output<'a>(args: &'a ToolArgs, key: &str) -> Result<Option<&'a str>, ToolError> {
    crate::common::parse_optional_output(args, key)
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

fn opt_pos(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match opt_f64(args, key)? {
        Some(v) if v > 0.0 && v.is_finite() => Ok(Some(v)),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a positive number"
        ))),
        None => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{FieldDef, FieldType, GeometryType, Layer, Ring};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a projected point layer with (x, y, value) rows.
    fn point_layer(rows: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (x, y, v) in rows {
            l.add_feature(Some(Geometry::point(*x, *y)), &[("v", (*v).into())])
                .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    /// A single vertical barrier line at the given x from y0 to y1.
    fn wall_layer(x: f64, y0: f64, y1: f64) -> String {
        let mut l = Layer::new("wall")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_feature(
            Some(Geometry::line_string(vec![
                Coord::xy(x, y0),
                Coord::xy(x, y1),
            ])),
            &[],
        )
        .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = InterpolateWithBarriersTool.run(&args, &ctx()).unwrap();
        let r = crate::common::load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, r)
    }

    /// Samples the interpolated value nearest a world coordinate.
    fn sample_at(r: &Raster, x: f64, y: f64) -> f64 {
        let col = ((x - r.x_min) / r.cell_size_x).floor() as isize;
        let y_top = r.y_min + r.rows as f64 * r.cell_size_y;
        let row = ((y_top - y) / r.cell_size_y).floor() as isize;
        r.get(0, row, col)
    }

    /// Core property: with a high-value cluster on the left of a near-total
    /// vertical barrier and low values on the right, a probe cell just right of
    /// the barrier is estimated MUCH lower when the barrier is respected than
    /// when it is ignored.
    #[test]
    fn barrier_lowers_opposite_side_estimate() {
        // Left cluster (high), right cluster (low); barrier at x=50, y in [5,95].
        let mut rows = Vec::new();
        for &y in &[20.0, 40.0, 60.0, 80.0] {
            rows.push((20.0, y, 100.0)); // left, high
            rows.push((80.0, y, 0.0)); // right, low
        }
        let pts = point_layer(&rows);

        // Barrier-free control.
        let (_o0, free) = run(json!({
            "input": pts, "field": "v", "method": "idw", "power": 2.0, "cell_size": 2.0
        }));
        // Barrier-aware. Probe cell (52, 50) sits just right of the wall.
        let (o1, blocked) = run(json!({
            "input": pts, "field": "v", "barriers": wall_layer(50.0, 5.0, 95.0),
            "method": "idw", "power": 2.0, "cell_size": 2.0
        }));

        let free_v = sample_at(&free, 52.0, 50.0);
        let blocked_v = sample_at(&blocked, 52.0, 50.0);
        assert!(
            blocked_v < free_v - 10.0,
            "barrier-aware ({blocked_v:.2}) should be well below barrier-free ({free_v:.2})"
        );
        // Right-side, gap-free probe should be dominated by the low cluster.
        assert!(
            blocked_v < 30.0,
            "probe just past the wall should read low, got {blocked_v:.2}"
        );
        assert!(o1.outputs["barrier_cells"].as_u64().unwrap() > 0);
    }

    /// A sample sitting exactly on a cell reproduces its value there (IDW clamp).
    #[test]
    fn exact_hit_recovers_value() {
        let pts = point_layer(&[(10.0, 10.0, 42.0), (90.0, 90.0, 7.0)]);
        let (_o, r) = run(json!({ "input": pts, "field": "v", "cell_size": 2.0 }));
        let v = sample_at(&r, 10.0, 10.0);
        assert!(
            (v - 42.0).abs() < 1.0,
            "near a sample the value ~= 42, got {v}"
        );
    }

    /// A polygon barrier enclosing a pocket with no interior samples leaves that
    /// pocket as no-data (values cannot bleed through the ring).
    #[test]
    fn enclosed_pocket_is_nodata() {
        // Samples form a frame around a central 30..70 box; a polygon ring seals
        // the box, and the search radius is too short to route around it.
        let mut rows = Vec::new();
        for &x in &[10.0, 90.0] {
            for &y in &[10.0, 50.0, 90.0] {
                rows.push((x, y, 50.0));
            }
        }
        let pts = point_layer(&rows);
        let mut l = Layer::new("poly")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        let ring = Ring::new(vec![
            Coord::xy(30.0, 30.0),
            Coord::xy(70.0, 30.0),
            Coord::xy(70.0, 70.0),
            Coord::xy(30.0, 70.0),
            Coord::xy(30.0, 30.0),
        ]);
        l.push(wbvector::Feature {
            fid: 0,
            geometry: Some(Geometry::Polygon {
                exterior: ring,
                interiors: vec![],
            }),
            attributes: vec![],
        });
        let id = wbvector::memory_store::put_vector(l);
        let barriers = wbvector::memory_store::make_vector_memory_path(&id);

        let (o, r) = run(json!({
            "input": pts, "field": "v", "barriers": barriers,
            "cell_size": 2.0, "radius": 5.0
        }));
        assert!(o.outputs["nodata_cells"].as_u64().unwrap() > 0);
        // Centre of the sealed box, with no interior sample and radius 5, is
        // unreachable -> no-data.
        assert_eq!(sample_at(&r, 50.0, 50.0), OUT_NODATA);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            InterpolateWithBarriersTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err()); // no field
        assert!(bad(json!({ "input": "a.geojson", "field": "v", "method": "kriging" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "field": "v", "power": -1 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "field": "v" })).is_ok());
    }
}
