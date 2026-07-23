//! GeoLibre tool: heat-equation (diffusion) interpolation of scattered point
//! values that flows *around* no-flux barriers.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Diffusion Interpolation With Barriers*
//! (Geostatistical Analyst). Scattered samples are rasterized onto an output
//! grid as fixed (Dirichlet) cells; every other cell relaxes toward the average
//! of its valid neighbours by iterating the Laplace / heat equation, so the
//! interpolated field is the steady state of heat diffusing from the samples.
//!
//! Barriers (polylines or polygons) are rasterized into a `blocked` mask: a
//! blocked cell neither diffuses nor is written out, and it is treated as a
//! **no-flux (Neumann) boundary** by its neighbours — a neighbour that is
//! blocked or off-grid is simply excluded from the local average and the
//! remaining neighbours are renormalised (reflecting the cell's own value). The
//! grid edges are no-flux for the same reason. Influence therefore spreads by
//! *diffusion distance* around obstacles, not through them.
//!
//! How this differs from the sibling tools already in the suite:
//!
//! * `interpolate_with_barriers` routes influence along the **shortest
//!   non-barrier path** (a cost/geodesic-distance IDW / local-polynomial
//!   kernel). Diffusion instead solves the heat equation, so influence decays
//!   with the (physically correct) diffusion distance and a barrier *slows*
//!   cross-flow rather than simply lengthening a path.
//! * the bundled `anisotropic_diffusion_filter` is a raster **smoother** applied
//!   to an already-dense image; this tool is a **sparse-point interpolator** that
//!   conditions the field on scattered samples.
//!
//! For a geographic (EPSG:4326) input the geometry is projected to a local
//! equirectangular metre frame centred on the extent (so `cell_size` and the
//! diffusion stencil are isotropic in true metres); the output raster's
//! georeferencing is converted back to degrees with distinct x/y cell sizes so
//! it still overlays the input. For a projected input everything is in the CRS's
//! native units. Output is an F32 raster; the diffusion is deterministic (no
//! RNG, no time).

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
use wbvector::{Coord, Geometry};

use crate::common::{parse_optional_output, write_or_store_output};
use crate::vector_common::{load_input_layer, parse_optional_str};

/// Mean Earth radius (metres) for the local equirectangular projection.
const EARTH_R: f64 = 6_371_000.0;
const OUT_NODATA: f64 = -9999.0;
/// Hard cap on grid dimensions to keep a single run tractable.
const MAX_DIM: usize = 4000;

pub struct DiffusionInterpolationWithBarriersTool;

impl Tool for DiffusionInterpolationWithBarriersTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "diffusion_interpolation_with_barriers",
            display_name: "Diffusion Interpolation With Barriers",
            summary: "Interpolate scattered point measurements to a raster by simulating heat diffusion that flows AROUND absolute (no-flux) barriers rather than through them (like ArcGIS Diffusion Interpolation With Barriers). Samples are pinned as fixed cells; every other cell relaxes toward the mean of its non-barrier neighbours over `number_iterations` Jacobi sweeps of the Laplace/heat equation, so influence decays with diffusion distance. Distinct from the shortest-non-barrier-path kernel of interpolate_with_barriers and from the raster smoother anisotropic_diffusion_filter.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer of sample measurements (Point / MultiPoint).",
                    required: true,
                },
                ToolParamSpec {
                    name: "z_field",
                    description: "Numeric field on the points holding the value to interpolate.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output interpolated raster (GeoTIFF). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Output cell size / resolution in map units (metres for a geographic CRS, CRS units otherwise). Required.",
                    required: true,
                },
                ToolParamSpec {
                    name: "barriers",
                    description: "Optional barrier vector layer (polylines or polygons) that blocks diffusion. Polylines block the cells their segments cross; polygons block the cells whose centre falls inside. Omit for a barrier-free control.",
                    required: false,
                },
                ToolParamSpec {
                    name: "bandwidth",
                    description: "Per-iteration relaxation weight lambda in (0, 0.25] (default 0.2): each free cell moves this fraction of the way toward its neighbour average per sweep.",
                    required: false,
                },
                ToolParamSpec {
                    name: "number_iterations",
                    description: "Number of Jacobi relaxation sweeps of the heat equation (default 100). More iterations spread influence farther.",
                    required: false,
                },
                ToolParamSpec {
                    name: "weight_field",
                    description: "Optional numeric field giving a per-point weight (default 1), used only to average coincident samples that fall in the same output cell.",
                    required: false,
                },
                ToolParamSpec {
                    name: "epsg",
                    description: "Override the input CRS EPSG code (e.g. when the layer is unlabeled). 4326 triggers the local metre projection.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "z_field")?;
        // cell_size is required and must be positive.
        match opt_f64(args, "cell_size")? {
            Some(v) if v > 0.0 && v.is_finite() => {}
            Some(_) => {
                return Err(ToolError::Validation(
                    "'cell_size' must be a positive number".to_string(),
                ))
            }
            None => {
                return Err(ToolError::Validation(
                    "missing required parameter 'cell_size'".to_string(),
                ))
            }
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let z_field = require_str(args, "z_field")?.to_string();
        let output = parse_optional_output(args, "output")?;
        let barriers_path = parse_optional_str(args, "barriers")?;
        let prm = parse_params(args)?;
        let cell = match opt_f64(args, "cell_size")? {
            Some(v) if v > 0.0 && v.is_finite() => v,
            _ => {
                return Err(ToolError::Validation(
                    "'cell_size' must be a positive number".to_string(),
                ))
            }
        };

        // ── Load samples (x, y, value, weight) in native coordinates ─────────
        let layer = load_input_layer(input)?;
        let zidx = layer
            .schema
            .field_index(&z_field)
            .ok_or_else(|| ToolError::Validation(format!("field '{z_field}' not found")))?;
        let widx =
            match prm.weight_field.as_deref() {
                Some(wf) => Some(layer.schema.field_index(wf).ok_or_else(|| {
                    ToolError::Validation(format!("weight_field '{wf}' not found"))
                })?),
                None => None,
            };

        // raw: (x, y, value, weight) in native coordinates.
        let mut raw: Vec<(f64, f64, f64, f64)> = Vec::new();
        for feat in layer.iter() {
            let Some(geom) = feat.geometry.as_ref() else {
                continue;
            };
            let Some(v) = feat.attributes.get(zidx).and_then(|a| a.as_f64()) else {
                continue;
            };
            if !v.is_finite() {
                continue;
            }
            let w = match widx {
                Some(i) => feat
                    .attributes
                    .get(i)
                    .and_then(|a| a.as_f64())
                    .filter(|w| w.is_finite() && *w > 0.0)
                    .unwrap_or(1.0),
                None => 1.0,
            };
            match geom {
                Geometry::Point(c) => raw.push((c.x, c.y, v, w)),
                Geometry::MultiPoint(cs) => {
                    for c in cs {
                        raw.push((c.x, c.y, v, w));
                    }
                }
                _ => {}
            }
        }
        if raw.is_empty() {
            return Err(ToolError::Execution(format!(
                "input contains no point features with a finite '{z_field}' value"
            )));
        }
        let sample_count = raw.len();

        // ── Native bounding box of the samples ───────────────────────────────
        let (mut nxmin, mut nymin, mut nxmax, mut nymax) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for (x, y, _, _) in &raw {
            nxmin = nxmin.min(*x);
            nxmax = nxmax.max(*x);
            nymin = nymin.min(*y);
            nymax = nymax.max(*y);
        }

        // ── Local metre projection for a geographic CRS; identity otherwise ──
        let epsg = prm.epsg.or_else(|| layer.crs_epsg());
        let geographic = epsg == Some(4326);
        let (kx, ky) = if geographic {
            let lat0 = 0.5 * (nymin + nymax);
            let ky = EARTH_R * std::f64::consts::PI / 180.0;
            let kx = ky * lat0.to_radians().cos().max(1e-9);
            (kx, ky)
        } else {
            (1.0, 1.0)
        };
        let (lon0, lat0) = (nxmin, nymin); // projection origin (native)
        let fwd_x = |x: f64| (x - lon0) * kx;
        let fwd_y = |y: f64| (y - lat0) * ky;

        // Working-frame (metre / native) sample points.
        let pts: Vec<(f64, f64, f64, f64)> = raw
            .iter()
            .map(|(x, y, v, w)| (fwd_x(*x), fwd_y(*y), *v, *w))
            .collect();

        // ── Grid extent = sample bbox padded by ~2 cells each side ───────────
        let ext_w = (nxmax - nxmin) * kx;
        let ext_h = (nymax - nymin) * ky;
        let pad = 2.0 * cell;
        let gxmin = -pad;
        let gymin = -pad;
        let gxmax = ext_w + pad;
        let gymax_raw = ext_h + pad;
        let cols = (((gxmax - gxmin) / cell).ceil() as usize).max(1);
        let rows = (((gymax_raw - gymin) / cell).ceil() as usize).max(1);
        let gymax = gymin + rows as f64 * cell; // snap top edge to whole cells
        if cols > MAX_DIM || rows > MAX_DIM {
            return Err(ToolError::Validation(format!(
                "grid {rows}x{cols} exceeds the {MAX_DIM} cap; increase 'cell_size'"
            )));
        }
        let n_cells = rows * cols;

        // ── Rasterize samples into fixed (Dirichlet) cells ───────────────────
        // Accumulate a weighted mean of the values landing in each cell.
        let mut wsum = vec![0.0_f64; n_cells];
        let mut vsum = vec![0.0_f64; n_cells];
        let mut fixed = vec![false; n_cells];
        let mut global_wsum = 0.0_f64;
        let mut global_vsum = 0.0_f64;
        for &(px, py, v, w) in &pts {
            global_wsum += w;
            global_vsum += w * v;
            let col = ((px - gxmin) / cell).floor();
            let row = ((gymax - py) / cell).floor();
            if col < 0.0 || row < 0.0 || col >= cols as f64 || row >= rows as f64 {
                continue;
            }
            let idx = row as usize * cols + col as usize;
            wsum[idx] += w;
            vsum[idx] += w * v;
            fixed[idx] = true;
        }
        let global_mean = if global_wsum > 0.0 {
            global_vsum / global_wsum
        } else {
            0.0
        };

        // ── Rasterize barriers into a blocked mask ───────────────────────────
        let mut blocked = vec![false; n_cells];
        if let Some(bpath) = barriers_path {
            let blayer = load_input_layer(bpath)?;
            for feat in blayer.iter() {
                if let Some(geom) = feat.geometry.as_ref() {
                    rasterize_barrier(
                        geom,
                        &fwd_x,
                        &fwd_y,
                        gxmin,
                        gymax,
                        cell,
                        rows,
                        cols,
                        &mut blocked,
                    );
                }
            }
        }
        // A sample landing on a barrier still seeds the field: unblock fixed cells.
        for idx in 0..n_cells {
            if fixed[idx] {
                blocked[idx] = false;
            }
        }
        let barrier_cells = blocked.iter().filter(|&&b| b).count();

        // ── Initialise the field ─────────────────────────────────────────────
        // Fixed cells = weighted sample mean; free non-blocked cells = global mean.
        let mut cur = vec![global_mean; n_cells];
        for idx in 0..n_cells {
            if fixed[idx] && wsum[idx] > 0.0 {
                cur[idx] = vsum[idx] / wsum[idx];
            }
        }

        ctx.progress.info(&format!(
            "{sample_count} sample(s) -> {rows}x{cols} grid, {barrier_cells} barrier cell(s), \
             {} iteration(s), lambda={:.3}",
            prm.number_iterations, prm.bandwidth
        ));

        // ── Iterate Jacobi relaxations of the heat equation ──────────────────
        let mut next = cur.clone();
        let lambda = prm.bandwidth;
        for it in 0..prm.number_iterations {
            for r in 0..rows {
                for c in 0..cols {
                    let idx = r * cols + c;
                    if blocked[idx] || fixed[idx] {
                        next[idx] = cur[idx];
                        continue;
                    }
                    // Average of in-grid, non-blocked 4-neighbours (no-flux
                    // Neumann boundary at barriers and grid edges: invalid
                    // neighbours are excluded and the rest renormalised).
                    let mut sum = 0.0;
                    let mut n = 0usize;
                    if r > 0 {
                        let ni = idx - cols;
                        if !blocked[ni] {
                            sum += cur[ni];
                            n += 1;
                        }
                    }
                    if r + 1 < rows {
                        let ni = idx + cols;
                        if !blocked[ni] {
                            sum += cur[ni];
                            n += 1;
                        }
                    }
                    if c > 0 {
                        let ni = idx - 1;
                        if !blocked[ni] {
                            sum += cur[ni];
                            n += 1;
                        }
                    }
                    if c + 1 < cols {
                        let ni = idx + 1;
                        if !blocked[ni] {
                            sum += cur[ni];
                            n += 1;
                        }
                    }
                    if n == 0 {
                        next[idx] = cur[idx];
                    } else {
                        let avg = sum / n as f64;
                        next[idx] = cur[idx] + lambda * (avg - cur[idx]);
                    }
                }
            }
            std::mem::swap(&mut cur, &mut next);
            if it % 16 == 0 {
                ctx.progress
                    .progress((it as f64 + 1.0) / prm.number_iterations as f64);
            }
        }

        // ── Emit raster: non-blocked cells get their value, blocked -> nodata ─
        let mut data = vec![OUT_NODATA; n_cells];
        let (mut vmin, mut vmax) = (f64::INFINITY, f64::NEG_INFINITY);
        let mut filled = 0usize;
        for idx in 0..n_cells {
            if blocked[idx] {
                continue;
            }
            let v = cur[idx];
            if v.is_finite() {
                data[idx] = v;
                vmin = vmin.min(v);
                vmax = vmax.max(v);
                filled += 1;
            }
        }

        // ── Georeference back to native units ────────────────────────────────
        let (out_cell_x, out_cell_y) = if geographic {
            (cell / kx, cell / ky)
        } else {
            (cell, cell)
        };
        let out_xmin = lon0 + gxmin / kx;
        let out_ymin = lat0 + gymin / ky;

        let mut out = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: out_xmin,
            y_min: out_ymin,
            cell_size: out_cell_x,
            cell_size_y: Some(out_cell_y),
            nodata: OUT_NODATA,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg,
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

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("sample_count".to_string(), json!(sample_count));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        outputs.insert("iterations".to_string(), json!(prm.number_iterations));
        outputs.insert("barrier_cells".to_string(), json!(barrier_cells));
        outputs.insert("filled_cells".to_string(), json!(filled));
        outputs.insert(
            "value_min".to_string(),
            json!(if vmin.is_finite() { vmin } else { 0.0 }),
        );
        outputs.insert(
            "value_max".to_string(),
            json!(if vmax.is_finite() { vmax } else { 0.0 }),
        );
        Ok(ToolRunResult { outputs })
    }
}

// ── Barrier rasterization ────────────────────────────────────────────────────

/// Marks cells blocked by a barrier geometry. Polylines block the cells their
/// segments cross (sampled at half-cell steps); polygons block the cells whose
/// centre falls inside the polygon (exterior minus holes).
#[allow(clippy::too_many_arguments)]
fn rasterize_barrier(
    geom: &Geometry,
    fwd_x: &impl Fn(f64) -> f64,
    fwd_y: &impl Fn(f64) -> f64,
    gxmin: f64,
    gymax: f64,
    cell: f64,
    rows: usize,
    cols: usize,
    blocked: &mut [bool],
) {
    let cell_of = |px: f64, py: f64| -> Option<usize> {
        let col = ((px - gxmin) / cell).floor();
        let row = ((gymax - py) / cell).floor();
        if col < 0.0 || row < 0.0 || col >= cols as f64 || row >= rows as f64 {
            return None;
        }
        Some(row as usize * cols + col as usize)
    };
    let mark_seg = |a: &Coord, b: &Coord, blocked: &mut [bool]| {
        let (ax, ay) = (fwd_x(a.x), fwd_y(a.y));
        let (bx, by) = (fwd_x(b.x), fwd_y(b.y));
        let dx = bx - ax;
        let dy = by - ay;
        let seg = (dx.hypot(dy) / (0.5 * cell)).ceil().max(1.0) as usize;
        for k in 0..=seg {
            let t = k as f64 / seg as f64;
            if let Some(idx) = cell_of(ax + t * dx, ay + t * dy) {
                blocked[idx] = true;
            }
        }
    };
    let mark_line = |coords: &[Coord], blocked: &mut [bool]| {
        for w in coords.windows(2) {
            mark_seg(&w[0], &w[1], blocked);
        }
    };
    // Projected rings for a polygon fill test.
    let fill_polygon =
        |exterior: &[(f64, f64)], holes: &[Vec<(f64, f64)>], blocked: &mut [bool]| {
            if exterior.len() < 3 {
                return;
            }
            let (mut xmn, mut ymn, mut xmx, mut ymx) = (
                f64::INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
            );
            for &(x, y) in exterior {
                xmn = xmn.min(x);
                xmx = xmx.max(x);
                ymn = ymn.min(y);
                ymx = ymx.max(y);
            }
            let c0 = (((xmn - gxmin) / cell).floor() as isize).max(0) as usize;
            let c1 = (((xmx - gxmin) / cell).ceil() as isize).min(cols as isize) as usize;
            let r0 = (((gymax - ymx) / cell).floor() as isize).max(0) as usize;
            let r1 = (((gymax - ymn) / cell).ceil() as isize).min(rows as isize) as usize;
            for r in r0..r1 {
                let cy = gymax - (r as f64 + 0.5) * cell;
                for c in c0..c1 {
                    let cx = gxmin + (c as f64 + 0.5) * cell;
                    if point_in_ring(cx, cy, exterior)
                        && !holes.iter().any(|h| point_in_ring(cx, cy, h))
                    {
                        blocked[r * cols + c] = true;
                    }
                }
            }
        };
    let project_ring = |coords: &[Coord]| -> Vec<(f64, f64)> {
        coords.iter().map(|c| (fwd_x(c.x), fwd_y(c.y))).collect()
    };
    match geom {
        Geometry::LineString(cs) => mark_line(cs, blocked),
        Geometry::MultiLineString(lines) => {
            for l in lines {
                mark_line(l, blocked);
            }
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            let ext = project_ring(exterior.coords());
            let holes: Vec<Vec<(f64, f64)>> =
                interiors.iter().map(|r| project_ring(r.coords())).collect();
            fill_polygon(&ext, &holes, blocked);
        }
        Geometry::MultiPolygon(polys) => {
            for (ext_ring, int_rings) in polys {
                let ext = project_ring(ext_ring.coords());
                let holes: Vec<Vec<(f64, f64)>> =
                    int_rings.iter().map(|r| project_ring(r.coords())).collect();
                fill_polygon(&ext, &holes, blocked);
            }
        }
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                rasterize_barrier(g, fwd_x, fwd_y, gxmin, gymax, cell, rows, cols, blocked);
            }
        }
        _ => {}
    }
}

/// Ray-casting point-in-polygon test against a single ring of (x, y) vertices.
fn point_in_ring(px: f64, py: f64, ring: &[(f64, f64)]) -> bool {
    let n = ring.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = ring[i];
        let (xj, yj) = ring[j];
        if ((yi > py) != (yj > py)) && (px < (xj - xi) * (py - yi) / (yj - yi + f64::EPSILON) + xi)
        {
            inside = !inside;
        }
        j = i;
    }
    inside
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    bandwidth: f64,
    number_iterations: usize,
    weight_field: Option<String>,
    epsg: Option<u32>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let bandwidth = match opt_f64(args, "bandwidth")? {
        None => 0.2,
        Some(v) if v > 0.0 && v <= 0.25 && v.is_finite() => v,
        Some(_) => {
            return Err(ToolError::Validation(
                "'bandwidth' must be in the interval (0, 0.25]".to_string(),
            ))
        }
    };
    let number_iterations = match opt_f64(args, "number_iterations")? {
        None => 100,
        Some(v) if v >= 1.0 && v.is_finite() => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "'number_iterations' must be a positive integer".to_string(),
            ))
        }
    };
    let weight_field = match args
        .get("weight_field")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") => None,
        Some(s) => Some(s.to_string()),
    };
    let epsg = match opt_f64(args, "epsg")? {
        Some(v) if v > 0.0 => Some(v as u32),
        _ => None,
    };
    Ok(Params {
        bandwidth,
        number_iterations,
        weight_field,
        epsg,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
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

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a projected point layer from (x, y, value) triples.
    fn point_layer(pts: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (x, y, v) in pts {
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
        let out = DiffusionInterpolationWithBarriersTool
            .run(&args, &ctx())
            .unwrap();
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

    /// Two samples, no barrier: the field obeys the maximum principle — every
    /// filled cell lies within [min, max] — and a cell between them is strictly
    /// intermediate.
    #[test]
    fn maximum_principle_and_monotonic() {
        let pts = point_layer(&[(10.0, 50.0, 0.0), (90.0, 50.0, 100.0)]);
        let (o, r) = run(json!({
            "input": pts, "z_field": "v", "cell_size": 4.0, "number_iterations": 300
        }));
        // Every filled cell within [0, 100].
        for row in 0..r.rows as isize {
            for col in 0..r.cols as isize {
                let v = r.get(0, row, col);
                if v != r.nodata {
                    assert!(
                        (-1e-6..=100.0 + 1e-6).contains(&v),
                        "cell value {v} outside [0,100]"
                    );
                }
            }
        }
        // Midpoint cell is strictly between the two samples.
        let mid = sample_at(&r, 50.0, 50.0);
        assert!(
            mid > 1.0 && mid < 99.0,
            "mid cell should be intermediate, got {mid}"
        );
        // Closer to the high sample reads higher than closer to the low sample.
        let near_hi = sample_at(&r, 78.0, 50.0);
        let near_lo = sample_at(&r, 22.0, 50.0);
        assert!(
            near_hi > mid && mid > near_lo,
            "field should be monotonic left->right"
        );
        assert_eq!(o.outputs["barrier_cells"].as_u64().unwrap(), 0);
    }

    /// A constant sample field diffuses to that constant everywhere.
    #[test]
    fn constant_field() {
        let c = 7.5;
        let pts = point_layer(&[
            (10.0, 10.0, c),
            (90.0, 10.0, c),
            (10.0, 90.0, c),
            (90.0, 90.0, c),
            (50.0, 50.0, c),
        ]);
        let (_o, r) = run(json!({
            "input": pts, "z_field": "v", "cell_size": 5.0, "number_iterations": 50
        }));
        for row in 0..r.rows as isize {
            for col in 0..r.cols as isize {
                let v = r.get(0, row, col);
                if v != r.nodata {
                    assert!((v - c).abs() < 1e-6, "constant field cell {v} != {c}");
                }
            }
        }
    }

    /// A wall spanning the full grid height between two samples partitions the
    /// grid: cells just left of the wall stay near the left sample and cells
    /// just right stay near the right sample (the barrier blocks cross-flow).
    #[test]
    fn barrier_blocks_cross_diffusion() {
        let pts = point_layer(&[(25.0, 50.0, 100.0), (75.0, 50.0, 0.0)]);
        // Wall at x=50 spanning well beyond the padded grid height.
        let barriers = wall_layer(50.0, -50.0, 150.0);
        let (o, r) = run(json!({
            "input": pts, "z_field": "v", "barriers": barriers,
            "cell_size": 4.0, "number_iterations": 400
        }));
        assert!(o.outputs["barrier_cells"].as_u64().unwrap() > 0);
        let left = sample_at(&r, 46.0, 50.0);
        let right = sample_at(&r, 54.0, 50.0);
        assert!(
            left > right + 50.0,
            "barrier should keep left ({left:.2}) high and right ({right:.2}) low"
        );
        // Contrast: without the wall the two probes are close together.
        let (_o2, r2) = run(json!({
            "input": pts, "z_field": "v", "cell_size": 4.0, "number_iterations": 400
        }));
        let left_free = sample_at(&r2, 46.0, 50.0);
        let right_free = sample_at(&r2, 54.0, 50.0);
        assert!(
            (left_free - right_free).abs() < (left - right),
            "barrier should widen the left/right gap vs the barrier-free control"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            DiffusionInterpolationWithBarriersTool.validate(&args)
        };
        assert!(bad(json!({})).is_err()); // missing input
        assert!(bad(json!({ "input": "a.geojson" })).is_err()); // missing z_field
        assert!(bad(json!({ "input": "a.geojson", "z_field": "v" })).is_err()); // missing cell_size
        assert!(bad(json!({ "input": "a.geojson", "z_field": "v", "cell_size": 0.0 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "z_field": "v", "cell_size": -1.0 })).is_err());
        assert!(bad(
            json!({ "input": "a.geojson", "z_field": "v", "cell_size": 1.0, "bandwidth": 0.0 })
        )
        .is_err());
        assert!(bad(
            json!({ "input": "a.geojson", "z_field": "v", "cell_size": 1.0, "bandwidth": 0.3 })
        )
        .is_err());
        assert!(bad(json!({ "input": "a.geojson", "z_field": "v", "cell_size": 1.0 })).is_ok());
        assert!(bad(
            json!({ "input": "a.geojson", "z_field": "v", "cell_size": 1.0, "bandwidth": 0.25 })
        )
        .is_ok());
    }
}
