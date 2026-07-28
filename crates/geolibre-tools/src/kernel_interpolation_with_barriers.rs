//! GeoLibre tool: kernel interpolation over geodesic distances around barriers.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Kernel Interpolation With Barriers*
//! (Geostatistical Analyst).
//!
//! GeoLibre ships `diffusion_interpolation_with_barriers`, the *other*
//! barrier-aware interpolator, and the two solve the problem differently:
//!
//! * **Diffusion** integrates a heat equation to steady state, treating barriers
//!   as no-flux boundaries. Smooth and physically motivated, but the result is
//!   only implicitly a function of the observations.
//! * **Kernel** (this tool) measures geodesic distance *around* barriers and
//!   then applies an explicit local kernel — a weighted polynomial fit, like
//!   `local_polynomial_interpolation` but with straight-line distance replaced by
//!   around-the-barrier distance.
//!
//! They produce different surfaces; the kernel variant is the one that is robust
//! to sparse data and gives a predictable, inspectable weighting. ArcGIS ships
//! both for that reason.
//!
//! Nothing else in either registry measures interpolation distance around
//! barriers: `idw_interpolation`, `radial_basis_function_interpolation`,
//! `thin_plate_spline`, `natural_neighbour_interpolation`, the kriging family
//! and `local_polynomial_interpolation` all use Euclidean distance. All of them
//! will happily interpolate a shoreline salinity value straight across a
//! headland, or a groundwater level across a fault.
//!
//! Geodesic distances come from one Dijkstra per observation over the
//! 8-connected grid with barrier cells excluded, bounded by the kernel
//! bandwidth. Fully deterministic — no RNG.

use std::collections::{BTreeMap, BinaryHeap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
use wbvector::{Coord, Geometry, Layer};

use crate::common::write_or_store_output;
use crate::vector_common::{load_input_layer, parse_optional_str};

/// Interpolates point observations with a local kernel weighted by
/// barrier-respecting geodesic distance.
pub struct KernelInterpolationWithBarriersTool;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kernel {
    Exponential,
    Gaussian,
    Quartic,
    Epanechnikov,
    Polynomial5,
    Constant,
}

impl Kernel {
    /// Weight for a normalised distance `t = d / bandwidth` in `[0, 1]`.
    fn weight(self, t: f64) -> f64 {
        let t = t.clamp(0.0, 1.0);
        match self {
            Self::Constant => 1.0,
            Self::Exponential => (-3.0 * t).exp(),
            Self::Gaussian => (-3.0 * t * t).exp(),
            Self::Quartic => {
                let u = 1.0 - t * t;
                u * u
            }
            Self::Epanechnikov => 1.0 - t * t,
            // ArcGIS's "fifth-order polynomial" kernel.
            Self::Polynomial5 => {
                let u = 1.0 - t * t;
                u * u * u
            }
        }
    }

    fn parse(s: &str) -> Result<Self, ToolError> {
        match s.trim().to_ascii_lowercase().as_str() {
            "exponential" => Ok(Self::Exponential),
            "gaussian" => Ok(Self::Gaussian),
            "quartic" => Ok(Self::Quartic),
            "epanechnikov" => Ok(Self::Epanechnikov),
            "polynomial5" | "fifth_order" => Ok(Self::Polynomial5),
            "constant" => Ok(Self::Constant),
            other => Err(ToolError::Validation(format!(
                "unknown kernel '{other}' (expected exponential, gaussian, quartic, epanechnikov, polynomial5 or constant)"
            ))),
        }
    }
}

impl Tool for KernelInterpolationWithBarriersTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "kernel_interpolation_with_barriers",
            display_name: "Kernel Interpolation With Barriers",
            summary: "Interpolates point observations with a moving-window kernel-weighted polynomial fit in which distance is measured as the shortest path AROUND absolute barriers rather than straight through them (ArcGIS Kernel Interpolation With Barriers). Complements the shipped diffusion_interpolation_with_barriers, which solves a heat equation instead of applying an explicit kernel. Every other interpolator in either registry uses Euclidean distance, so all of them interpolate straight across headlands and faults.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point layer carrying the observations.",
                    required: true,
                },
                ToolParamSpec {
                    name: "z_field",
                    description: "Numeric field to interpolate.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output raster path for predicted values. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_error",
                    description: "Optional path for the local weighted residual RMS raster — the weighted spread of neighbour residuals about the local fit. This is a goodness-of-fit measure, NOT an ArcGIS-style prediction standard error (it accounts for neither effective degrees of freedom nor leverage). If omitted, stored in memory (still returned).",
                    required: false,
                },
                ToolParamSpec {
                    name: "barriers",
                    description: "Optional polyline layer of absolute barriers. With none supplied the tool degenerates to plain kernel interpolation.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Output cell size in CRS units. Defaults to roughly 1/100 of the extent's larger side.",
                    required: false,
                },
                ToolParamSpec {
                    name: "kernel",
                    description: "exponential, gaussian, quartic, epanechnikov, polynomial5 (default) or constant.",
                    required: false,
                },
                ToolParamSpec {
                    name: "bandwidth",
                    description: "Kernel bandwidth in CRS units. Defaults to a data-driven value from mean nearest-neighbour spacing.",
                    required: false,
                },
                ToolParamSpec {
                    name: "power",
                    description: "Local polynomial order: 0 (weighted mean) or 1 (default, planar fit).",
                    required: false,
                },
                ToolParamSpec {
                    name: "ridge",
                    description: "Ridge parameter stabilising the local least-squares fit where a neighbourhood is degenerate (default 1e-6).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_neighbors",
                    description: "Cap on the number of observations contributing to each local fit (default 32).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["input", "z_field"] {
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
        let input = required_str(args, "input")?;
        let z_field = required_str(args, "z_field")?;
        let output = parse_optional_str(args, "output")?;
        let output_error = parse_optional_str(args, "output_error")?;
        let barriers_path = parse_optional_str(args, "barriers")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let zi = layer.schema.field_index(z_field).ok_or_else(|| {
            ToolError::Validation(format!("z_field '{z_field}' not found on the input layer"))
        })?;

        // Observations: (x, y, z).
        let mut obs: Vec<(f64, f64, f64)> = Vec::new();
        for feature in layer.iter() {
            let Some(Geometry::Point(c)) = feature.geometry.as_ref() else {
                continue;
            };
            let Some(z) = numeric(&feature.attributes[zi]) else {
                continue;
            };
            obs.push((c.x, c.y, z));
        }
        if obs.is_empty() {
            return Err(ToolError::Execution(format!(
                "no point features with a numeric '{z_field}' value"
            )));
        }
        ctx.progress.info(&format!("{} observation(s)", obs.len()));

        // ── Output grid ──────────────────────────────────────────────────────
        let (mut min_x, mut min_y) = (f64::INFINITY, f64::INFINITY);
        let (mut max_x, mut max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
        for (x, y, _) in &obs {
            min_x = min_x.min(*x);
            max_x = max_x.max(*x);
            min_y = min_y.min(*y);
            max_y = max_y.max(*y);
        }
        let span_x = (max_x - min_x).abs();
        let span_y = (max_y - min_y).abs();
        let span = span_x.max(span_y);
        let cell = match prm.cell_size {
            Some(c) => c,
            None => {
                if span > 0.0 {
                    span / 100.0
                } else {
                    1.0
                }
            }
        };
        if !cell.is_finite() || cell <= 0.0 {
            return Err(ToolError::Validation(
                "'cell_size' must be a positive, finite number".to_string(),
            ));
        }
        let cols = ((span_x / cell).ceil() as usize + 1).max(1);
        let rows = ((span_y / cell).ceil() as usize + 1).max(1);
        if rows.saturating_mul(cols) > 40_000_000 {
            return Err(ToolError::Validation(format!(
                "requested grid is {rows}x{cols}; increase 'cell_size'"
            )));
        }
        // One distance grid is allocated per observation, so the real cost is
        // the product. Bounding only the grid turns a large point layer into an
        // OOM instead of a validation error.
        const MAX_DISTANCE_CELLS: usize = 200_000_000;
        let budget = rows.saturating_mul(cols).saturating_mul(obs.len());
        if budget > MAX_DISTANCE_CELLS {
            return Err(ToolError::Validation(format!(
                "{} observation(s) over a {rows}x{cols} grid needs {budget} distance cells (limit {MAX_DISTANCE_CELLS}); increase 'cell_size' or reduce the observation count",
                obs.len()
            )));
        }
        let y_max = min_y + rows as f64 * cell;

        // Bandwidth: default from mean nearest-neighbour spacing, which adapts
        // to how dense the observations actually are.
        let bandwidth = match prm.bandwidth {
            Some(b) => b,
            None => default_bandwidth(&obs),
        };
        if !bandwidth.is_finite() || bandwidth <= 0.0 {
            return Err(ToolError::Validation(
                "'bandwidth' must be a positive, finite number".to_string(),
            ));
        }
        ctx.progress
            .info(&format!("grid {rows}x{cols}, bandwidth {bandwidth:.4}"));

        // ── Barrier mask ─────────────────────────────────────────────────────
        let blocked = match barriers_path {
            Some(path) => {
                let bl = load_input_layer(path)?;
                rasterize_barriers(&bl, min_x, y_max, cell, rows, cols)
            }
            None => vec![false; rows * cols],
        };

        // ── Geodesic distance from each observation ──────────────────────────
        // One Dijkstra per observation, bounded by the bandwidth so the search
        // stays local rather than sweeping the whole grid.
        ctx.progress.info("computing geodesic distances");
        let mut dists: Vec<Vec<f64>> = Vec::with_capacity(obs.len());
        for (oi, (ox, oy, _)) in obs.iter().enumerate() {
            let col = (((ox - min_x) / cell).floor() as isize).clamp(0, cols as isize - 1) as usize;
            let row = (((y_max - oy) / cell).floor() as isize).clamp(0, rows as isize - 1) as usize;
            dists.push(dijkstra_from(
                row, col, rows, cols, cell, &blocked, bandwidth,
            ));
            ctx.progress
                .progress((oi as f64 + 1.0) / obs.len() as f64 * 0.5);
        }

        // ── Local weighted polynomial fit per cell ───────────────────────────
        ctx.progress.info("fitting local kernels");
        let nodata = -9999.0_f64;
        let mut pred = vec![nodata; rows * cols];
        let mut err = vec![nodata; rows * cols];
        let mut filled = 0_u64;
        let mut unreachable = 0_u64;

        for row in 0..rows {
            for col in 0..cols {
                let idx = row * cols + col;
                if blocked[idx] {
                    continue;
                }
                // Gather reachable observations within the bandwidth.
                let mut near: Vec<(usize, f64)> = Vec::new();
                for (oi, d) in dists.iter().enumerate() {
                    let dv = d[idx];
                    if dv.is_finite() && dv <= bandwidth {
                        near.push((oi, dv));
                    }
                }
                if near.is_empty() {
                    // Enclosed by barriers with no observation inside: no
                    // honest value exists, so leave it no-data rather than
                    // silently extrapolating.
                    unreachable += 1;
                    continue;
                }
                near.sort_by(|a, b| {
                    a.1.partial_cmp(&b.1)
                        .unwrap_or(std::cmp::Ordering::Equal)
                        .then_with(|| a.0.cmp(&b.0))
                });
                near.truncate(prm.max_neighbors);

                let cx = min_x + (col as f64 + 0.5) * cell;
                let cy = y_max - (row as f64 + 0.5) * cell;

                if let Some((value, stderr)) = local_fit(
                    &obs, &near, cx, cy, bandwidth, prm.kernel, prm.power, prm.ridge,
                ) {
                    pred[idx] = value;
                    err[idx] = stderr;
                    filled += 1;
                }
            }
            ctx.progress
                .progress(0.5 + (row as f64 + 1.0) / rows as f64 * 0.5);
        }

        let make = |data: Vec<f64>| -> Result<Raster, ToolError> {
            let mut r = Raster::new(RasterConfig {
                cols,
                rows,
                bands: 1,
                x_min: min_x,
                y_min: min_y,
                cell_size: cell,
                cell_size_y: Some(cell),
                nodata,
                data_type: DataType::F32,
                // Carry both the EPSG and the WKT: a layer defined only by WKT
                // would otherwise produce an unreferenced raster.
                crs: CrsInfo {
                    epsg: layer.crs_epsg(),
                    wkt: layer.crs_wkt().map(str::to_string),
                    proj4: None,
                },
                metadata: Vec::new(),
            });
            for row in 0..rows {
                for col in 0..cols {
                    r.set(0, row as isize, col as isize, data[row * cols + col])
                        .map_err(|e| ToolError::Execution(format!("failed writing cell: {e}")))?;
                }
            }
            Ok(r)
        };

        let out_path = write_or_store_output(make(pred)?, output)?;
        let err_path = write_or_store_output(make(err)?, output_error)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("output_error".to_string(), json!(err_path));
        outputs.insert(
            "output_error_semantics".to_string(),
            json!("local weighted residual RMS (not a prediction standard error)"),
        );
        outputs.insert("observation_count".to_string(), json!(obs.len()));
        outputs.insert("filled_cells".to_string(), json!(filled));
        outputs.insert("unreachable_cells".to_string(), json!(unreachable));
        outputs.insert("bandwidth".to_string(), json!(bandwidth));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        Ok(ToolRunResult { outputs })
    }
}

/// Weighted least-squares fit of a local polynomial, returning the value at
/// `(cx, cy)` and the weighted RMS of the neighbour residuals about that fit.
///
/// Deliberately *not* a prediction standard error: no effective-dof or leverage
/// correction is applied, and for `power == 0` the residuals are taken about the
/// same weighted mean being reported, so the value goes to zero on a locally
/// constant neighbourhood. Named and documented as residual RMS so it is not
/// mistaken for ArcGIS's uncertainty surface.
#[allow(clippy::too_many_arguments)]
fn local_fit(
    obs: &[(f64, f64, f64)],
    near: &[(usize, f64)],
    cx: f64,
    cy: f64,
    bandwidth: f64,
    kernel: Kernel,
    power: usize,
    ridge: f64,
) -> Option<(f64, f64)> {
    let weights: Vec<f64> = near
        .iter()
        .map(|(_, d)| kernel.weight(d / bandwidth))
        .collect();
    let wsum: f64 = weights.iter().sum();
    if wsum <= 0.0 {
        return None;
    }

    if power == 0 || near.len() < 3 {
        // Order 0 collapses to a weighted mean. Also the fallback when there
        // are too few neighbours to define a plane.
        let mean: f64 = near
            .iter()
            .zip(&weights)
            .map(|((oi, _), w)| obs[*oi].2 * w)
            .sum::<f64>()
            / wsum;
        let var: f64 = near
            .iter()
            .zip(&weights)
            .map(|((oi, _), w)| {
                let r = obs[*oi].2 - mean;
                w * r * r
            })
            .sum::<f64>()
            / wsum;
        return Some((mean, var.max(0.0).sqrt()));
    }

    // Order 1: fit z = a + b*(x - cx) + c*(y - cy) by weighted least squares.
    // Normal equations are 3x3, small enough to solve directly — no linear
    // algebra crate needed, which matters for the no-heavy-deps constraint.
    let mut ata = [[0.0_f64; 3]; 3];
    let mut atb = [0.0_f64; 3];
    for ((oi, _), w) in near.iter().zip(&weights) {
        let (x, y, z) = obs[*oi];
        let basis = [1.0, x - cx, y - cy];
        for i in 0..3 {
            for j in 0..3 {
                ata[i][j] += w * basis[i] * basis[j];
            }
            atb[i] += w * basis[i] * z;
        }
    }
    // Ridge on the diagonal keeps a degenerate (e.g. collinear) neighbourhood
    // solvable instead of producing a wild fit.
    for (i, row) in ata.iter_mut().enumerate() {
        row[i] += ridge * wsum;
    }

    let coef = solve3(ata, atb)?;
    // The fit is centred on the cell, so the value there is just the constant.
    let value = coef[0];

    let var: f64 = near
        .iter()
        .zip(&weights)
        .map(|((oi, _), w)| {
            let (x, y, z) = obs[*oi];
            let fitted = coef[0] + coef[1] * (x - cx) + coef[2] * (y - cy);
            let r = z - fitted;
            w * r * r
        })
        .sum::<f64>()
        / wsum;

    Some((value, var.max(0.0).sqrt()))
}

/// Gaussian elimination with partial pivoting on a 3x3 system.
fn solve3(mut a: [[f64; 3]; 3], mut b: [f64; 3]) -> Option<[f64; 3]> {
    for col in 0..3 {
        let mut piv = col;
        for r in (col + 1)..3 {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-14 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        for r in (col + 1)..3 {
            let f = a[r][col] / a[col][col];
            let pivot_row = a[col];
            for (c, v) in a[r].iter_mut().enumerate().skip(col) {
                *v -= f * pivot_row[c];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = [0.0_f64; 3];
    for i in (0..3).rev() {
        let mut s = b[i];
        for j in (i + 1)..3 {
            s -= a[i][j] * x[j];
        }
        x[i] = s / a[i][i];
    }
    if x.iter().all(|v| v.is_finite()) {
        Some(x)
    } else {
        None
    }
}

/// Dijkstra over the 8-connected grid from one seed cell, skipping blocked
/// cells and stopping past `max_dist`.
fn dijkstra_from(
    seed_row: usize,
    seed_col: usize,
    rows: usize,
    cols: usize,
    cell: f64,
    blocked: &[bool],
    max_dist: f64,
) -> Vec<f64> {
    #[derive(PartialEq)]
    struct Node(f64, usize);
    impl Eq for Node {}
    impl Ord for Node {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering {
            other
                .0
                .partial_cmp(&self.0)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| other.1.cmp(&self.1))
        }
    }
    impl PartialOrd for Node {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }

    let mut dist = vec![f64::INFINITY; rows * cols];
    let seed = seed_row * cols + seed_col;
    // A seed that landed on a barrier still anchors its own cell, otherwise the
    // observation would be silently dropped.
    dist[seed] = 0.0;
    let mut heap = BinaryHeap::new();
    heap.push(Node(0.0, seed));

    let diag = cell * std::f64::consts::SQRT_2;
    let steps: [(isize, isize, f64); 8] = [
        (-1, 0, cell),
        (1, 0, cell),
        (0, -1, cell),
        (0, 1, cell),
        (-1, -1, diag),
        (-1, 1, diag),
        (1, -1, diag),
        (1, 1, diag),
    ];

    while let Some(Node(d, i)) = heap.pop() {
        if d > dist[i] || d > max_dist {
            continue;
        }
        let row = i / cols;
        let col = i % cols;
        for (dr, dc, cost) in steps {
            let nr = row as isize + dr;
            let nc = col as isize + dc;
            if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                continue;
            }
            let j = nr as usize * cols + nc as usize;
            if blocked[j] {
                continue;
            }
            // No corner cutting: a diagonal step passes between two orthogonal
            // cells, and a rasterized diagonal barrier is exactly such a
            // staircase. Without this the 8-connected expansion crosses the wall
            // at its corners even though the burn left no cell gaps.
            if dr != 0 && dc != 0 {
                let a = (row as isize + dr) as usize * cols + col;
                let b = row * cols + (col as isize + dc) as usize;
                if blocked[a] && blocked[b] {
                    continue;
                }
            }
            let nd = d + cost;
            if nd < dist[j] && nd <= max_dist {
                dist[j] = nd;
                heap.push(Node(nd, j));
            }
        }
    }
    dist
}

/// Burns barrier polylines onto the grid with a DDA walk along each segment.
fn rasterize_barriers(
    layer: &Layer,
    min_x: f64,
    y_max: f64,
    cell: f64,
    rows: usize,
    cols: usize,
) -> Vec<bool> {
    let mut blocked = vec![false; rows * cols];
    let mut burn = |x: f64, y: f64| {
        let col = ((x - min_x) / cell).floor();
        let row = ((y_max - y) / cell).floor();
        if col >= 0.0 && row >= 0.0 && (col as usize) < cols && (row as usize) < rows {
            blocked[row as usize * cols + col as usize] = true;
        }
    };

    let mut walk = |cs: &[Coord]| {
        for w in cs.windows(2) {
            let (x0, y0, x1, y1) = (w[0].x, w[0].y, w[1].x, w[1].y);
            let len = (x1 - x0).hypot(y1 - y0);
            // Half-cell steps guarantee no gaps a diagonal path could slip
            // through, which would defeat the barrier entirely.
            let steps = ((len / (cell * 0.5)).ceil() as usize).max(1);
            for s in 0..=steps {
                let t = s as f64 / steps as f64;
                burn(x0 + t * (x1 - x0), y0 + t * (y1 - y0));
            }
        }
    };

    // Every ring of every part is burned, including polygon holes, MultiPolygon
    // members and nested collections — otherwise those barrier inputs are
    // silently invisible while the shipped interpolate_with_barriers honours
    // them.
    fn walk_geometry(geom: &Geometry, walk: &mut impl FnMut(&[Coord])) {
        match geom {
            Geometry::LineString(cs) => walk(cs),
            Geometry::MultiLineString(ls) => {
                for cs in ls {
                    walk(cs);
                }
            }
            Geometry::Polygon {
                exterior,
                interiors,
            } => {
                walk(&exterior.0);
                for hole in interiors {
                    walk(&hole.0);
                }
            }
            Geometry::MultiPolygon(parts) => {
                for (ext, holes) in parts {
                    walk(&ext.0);
                    for hole in holes {
                        walk(&hole.0);
                    }
                }
            }
            Geometry::GeometryCollection(gs) => {
                for g in gs {
                    walk_geometry(g, walk);
                }
            }
            // Points bound nothing, so they cannot act as a barrier.
            Geometry::Point(_) | Geometry::MultiPoint(_) => {}
        }
    }

    for feature in layer.iter() {
        if let Some(geom) = feature.geometry.as_ref() {
            walk_geometry(geom, &mut walk);
        }
    }
    blocked
}

/// Mean nearest-neighbour spacing, scaled up so a default neighbourhood holds
/// several observations rather than just the closest one.
///
/// The scan is all-pairs, so it is evaluated over at most `SAMPLE` probe points
/// (deterministically strided, not sampled at random). A default does not need
/// the exact mean, and an uncapped O(n^2) pass makes the tool appear to hang on
/// a large layer before any grid work starts.
fn default_bandwidth(obs: &[(f64, f64, f64)]) -> f64 {
    if obs.len() < 2 {
        return 1.0;
    }
    const SAMPLE: usize = 512;
    let stride = obs.len().div_ceil(SAMPLE).max(1);
    let mut total = 0.0;
    let mut probes = 0usize;
    for (i, (xi, yi, _)) in obs.iter().enumerate() {
        if i % stride != 0 {
            continue;
        }
        probes += 1;
        let mut best = f64::INFINITY;
        for (j, (xj, yj, _)) in obs.iter().enumerate() {
            if i == j {
                continue;
            }
            let d = (xi - xj).hypot(yi - yj);
            if d < best {
                best = d;
            }
        }
        if best.is_finite() {
            total += best;
        }
    }
    let mean = total / probes.max(1) as f64;
    if mean > 0.0 {
        mean * 4.0
    } else {
        1.0
    }
}

fn numeric(v: &wbvector::FieldValue) -> Option<f64> {
    match v {
        wbvector::FieldValue::Integer(i) => Some(*i as f64),
        wbvector::FieldValue::Float(f) => Some(*f),
        wbvector::FieldValue::Text(s) => s.trim().parse::<f64>().ok(),
        _ => None,
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    cell_size: Option<f64>,
    kernel: Kernel,
    bandwidth: Option<f64>,
    power: usize,
    ridge: f64,
    max_neighbors: usize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let kernel = match parse_optional_str(args, "kernel")? {
        Some(s) => Kernel::parse(s)?,
        None => Kernel::Polynomial5,
    };
    let power = opt_f64(args, "power")?.unwrap_or(1.0);
    if power != 0.0 && power != 1.0 {
        return Err(ToolError::Validation("'power' must be 0 or 1".to_string()));
    }
    let ridge = opt_f64(args, "ridge")?.unwrap_or(1e-6);
    // NaN passes every `<` comparison, so test finiteness explicitly: a NaN
    // ridge poisons the normal equations and `NaN as usize` truncates
    // max_neighbors to 0, emptying every neighbourhood.
    if !ridge.is_finite() || ridge < 0.0 {
        return Err(ToolError::Validation(
            "'ridge' must be non-negative".to_string(),
        ));
    }
    let max_neighbors = opt_f64(args, "max_neighbors")?.unwrap_or(32.0);
    if !max_neighbors.is_finite() || max_neighbors < 1.0 {
        return Err(ToolError::Validation(
            "'max_neighbors' must be at least 1".to_string(),
        ));
    }
    let cell_size = opt_f64(args, "cell_size")?;
    if let Some(c) = cell_size {
        if c <= 0.0 {
            return Err(ToolError::Validation(
                "'cell_size' must be positive".to_string(),
            ));
        }
    }
    let bandwidth = opt_f64(args, "bandwidth")?;
    if let Some(b) = bandwidth {
        if b <= 0.0 {
            return Err(ToolError::Validation(
                "'bandwidth' must be positive".to_string(),
            ));
        }
    }
    Ok(Params {
        cell_size,
        kernel,
        bandwidth,
        power: power as usize,
        ridge,
        max_neighbors: max_neighbors as usize,
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

fn required_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::load_input_raster;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, FieldValue, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn points(pts: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("obs");
        l.geom_type = Some(GeometryType::Point);
        l.add_field(FieldDef::new("z", FieldType::Float));
        for (x, y, z) in pts {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(*x, *y))),
                &[("z", FieldValue::Float(*z))],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn barrier(lines: &[Vec<(f64, f64)>]) -> String {
        let mut l = Layer::new("barriers");
        l.geom_type = Some(GeometryType::LineString);
        for pts in lines {
            let cs: Vec<Coord> = pts.iter().map(|(x, y)| Coord::xy(*x, *y)).collect();
            l.add_feature(Some(Geometry::LineString(cs)), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(extra: Value) -> (Raster, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(extra).unwrap();
        let res = KernelInterpolationWithBarriersTool
            .run(&args, &ctx())
            .unwrap();
        let r = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        (r, res)
    }

    /// With no barriers the surface honours the data: a constant field
    /// interpolates to that constant everywhere.
    #[test]
    fn constant_field_interpolates_to_the_constant() {
        let pts = points(&[
            (0.0, 0.0, 7.0),
            (10.0, 0.0, 7.0),
            (0.0, 10.0, 7.0),
            (10.0, 10.0, 7.0),
            (5.0, 5.0, 7.0),
        ]);
        let (out, _) = run(json!({
            "input": pts, "z_field": "z", "cell_size": 1.0, "bandwidth": 20.0
        }));
        let v = out.get(0, 5, 5);
        // The ridge term shrinks the fit toward zero very slightly, so the
        // tolerance is scaled to the default ridge (1e-6) rather than to
        // machine epsilon.
        assert!(
            (v - 7.0).abs() < 1e-3,
            "constant field must interpolate to ~7, got {v}"
        );
    }

    /// A first-order fit reproduces a planar trend, which a weighted mean
    /// cannot. This is what `power = 1` buys.
    #[test]
    fn planar_trend_is_reproduced_by_first_order_fit() {
        // z = x, sampled on a grid.
        let mut pts = Vec::new();
        for i in 0..=10 {
            for j in 0..=10 {
                pts.push((i as f64, j as f64, i as f64));
            }
        }
        let path = points(&pts);
        let (out, _) = run(json!({
            "input": path, "z_field": "z", "cell_size": 1.0,
            "bandwidth": 3.0, "power": 1, "kernel": "gaussian"
        }));
        // Sample near the middle of the grid; the true value at x=5 is 5.
        let cols = out.cols;
        let rows = out.rows;
        let mid_row = rows / 2;
        let mid_col = cols / 2;
        let v = out.get(0, mid_row as isize, mid_col as isize);
        let x = out.x_min + (mid_col as f64 + 0.5) * out.cell_size_x;
        assert!(
            (v - x).abs() < 0.5,
            "planar trend should give ~{x} at that cell, got {v}"
        );
    }

    /// The whole point of the tool: a barrier changes the result. Two
    /// observations either side of a wall must not blend across it the way an
    /// unblocked run does.
    #[test]
    fn barrier_changes_the_interpolated_surface() {
        // Low values on the left, high on the right, wall down the middle.
        let pts = points(&[
            (0.0, 5.0, 0.0),
            (1.0, 2.0, 0.0),
            (1.0, 8.0, 0.0),
            (10.0, 5.0, 100.0),
            (9.0, 2.0, 100.0),
            (9.0, 8.0, 100.0),
        ]);
        let wall = barrier(&[vec![(5.0, -1.0), (5.0, 11.0)]]);

        let (plain, _) = run(json!({
            "input": pts.clone(), "z_field": "z", "cell_size": 0.5, "bandwidth": 30.0
        }));
        let (walled, res) = run(json!({
            "input": pts, "z_field": "z", "cell_size": 0.5,
            "bandwidth": 30.0, "barriers": wall
        }));

        // Sample a cell just left of the wall.
        let col = ((4.0 - plain.x_min) / plain.cell_size_x) as isize;
        let row = (plain.rows / 2) as isize;
        let a = plain.get(0, row, col);
        let b = walled.get(0, row, col);
        assert!(a.is_finite() && b.is_finite());
        assert!(
            b < a - 1.0,
            "with a barrier the left side should stay closer to its own low values: unblocked {a}, walled {b}"
        );
        assert_eq!(res.outputs["observation_count"], json!(6));
    }

    /// A 45-degree barrier must be impermeable. A rasterized diagonal wall is a
    /// staircase, and an 8-connected search that ignores corners crosses it
    /// between the blocked cells even though the burn leaves no cell gaps — the
    /// axis-aligned tests cannot catch that.
    #[test]
    fn diagonal_barrier_is_not_permeable() {
        // Observations either side of a wall running corner to corner.
        let pts = points(&[(1.0, 9.0, 0.0), (9.0, 1.0, 100.0)]);
        let wall = barrier(&[vec![(-1.0, -1.0), (11.0, 11.0)]]);

        let (plain, _) = run(json!({
            "input": pts.clone(), "z_field": "z", "cell_size": 0.5, "bandwidth": 60.0
        }));
        let (walled, _) = run(json!({
            "input": pts, "z_field": "z", "cell_size": 0.5,
            "bandwidth": 60.0, "barriers": wall
        }));

        // Sample well inside the upper-left half, which holds only the 0-valued
        // observation. Without the wall the 100 leaks across; with it, it cannot.
        let col = ((2.0 - plain.x_min) / plain.cell_size_x) as isize;
        let row = ((plain.y_min + plain.rows as f64 * plain.cell_size_y - 8.0) / plain.cell_size_y)
            as isize;
        let a = plain.get(0, row, col);
        let b = walled.get(0, row, col);
        assert!(a.is_finite() && b.is_finite());
        assert!(
            a > 1.0,
            "without the wall the 100-valued observation should reach this cell, got {a}"
        );
        // A merely-reduced value would still pass if the 100 leaked across with
        // a longer detour, so require the far observation to be shut out
        // entirely: only the 0-valued side remains reachable.
        assert!(
            b.abs() < 1e-6,
            "the diagonal wall must be impermeable, so this cell sees only the 0-valued observation; got {b}"
        );
    }

    /// A polygon barrier is honoured, not silently ignored — only line barriers
    /// used to be walked, so a MultiPolygon barrier was invisible.
    #[test]
    fn polygon_barriers_are_rasterized() {
        let pts = points(&[(0.0, 5.0, 0.0), (10.0, 5.0, 100.0)]);
        let mut l = Layer::new("barrier_poly");
        l.geom_type = Some(GeometryType::Polygon);
        l.add_feature(
            Some(Geometry::MultiPolygon(vec![(
                wbvector::Ring::new(vec![
                    Coord::xy(4.5, -5.0),
                    Coord::xy(5.5, -5.0),
                    Coord::xy(5.5, 15.0),
                    Coord::xy(4.5, 15.0),
                ]),
                vec![],
            )])),
            &[],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        let poly_barrier = memory_store::make_vector_memory_path(&id);

        let (plain, _) = run(json!({
            "input": pts.clone(), "z_field": "z", "cell_size": 0.5, "bandwidth": 40.0
        }));
        let (walled, _) = run(json!({
            "input": pts, "z_field": "z", "cell_size": 0.5,
            "bandwidth": 40.0, "barriers": poly_barrier
        }));
        let col = ((4.0 - plain.x_min) / plain.cell_size_x) as isize;
        let row = (plain.rows / 2) as isize;
        assert!(
            walled.get(0, row, col) < plain.get(0, row, col) - 1.0,
            "a polygon barrier must block influence just as a line barrier does"
        );
    }

    /// A region fully enclosed by barriers with no observation inside gets
    /// no-data rather than a silently extrapolated value.
    #[test]
    fn fully_enclosed_region_is_nodata() {
        // Observations outside; a closed box in the middle with none inside.
        let pts = points(&[
            (0.0, 0.0, 1.0),
            (10.0, 0.0, 1.0),
            (0.0, 10.0, 1.0),
            (10.0, 10.0, 1.0),
        ]);
        let box_wall = barrier(&[vec![
            (4.0, 4.0),
            (6.0, 4.0),
            (6.0, 6.0),
            (4.0, 6.0),
            (4.0, 4.0),
        ]]);
        let (out, res) = run(json!({
            "input": pts, "z_field": "z", "cell_size": 0.25,
            "bandwidth": 50.0, "barriers": box_wall
        }));
        // The centre of the box is unreachable from every observation.
        let col = ((5.0 - out.x_min) / out.cell_size_x) as isize;
        let row =
            ((out.y_min + out.rows as f64 * out.cell_size_y - 5.0) / out.cell_size_y) as isize;
        assert_eq!(
            out.get(0, row, col),
            out.nodata,
            "an enclosed cell with no observation inside must be no-data"
        );
        assert!(res.outputs["unreachable_cells"].as_u64().unwrap() > 0);
    }

    /// The error surface is produced and is non-negative where defined.
    #[test]
    fn error_surface_is_emitted() {
        let pts = points(&[
            (0.0, 0.0, 1.0),
            (4.0, 0.0, 5.0),
            (0.0, 4.0, 3.0),
            (4.0, 4.0, 9.0),
        ]);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": pts, "z_field": "z", "cell_size": 1.0, "bandwidth": 10.0
        }))
        .unwrap();
        let res = KernelInterpolationWithBarriersTool
            .run(&args, &ctx())
            .unwrap();
        let err = load_input_raster(res.outputs["output_error"].as_str().unwrap()).unwrap();
        let mut seen = false;
        for r in 0..err.rows {
            for c in 0..err.cols {
                let v = err.get(0, r as isize, c as isize);
                if v != err.nodata {
                    assert!(v >= 0.0, "standard error must be non-negative, got {v}");
                    seen = true;
                }
            }
        }
        assert!(seen, "expected some defined error values");
    }

    /// Every kernel is accepted and produces a finite surface.
    #[test]
    fn all_kernels_run() {
        let pts = points(&[
            (0.0, 0.0, 1.0),
            (4.0, 0.0, 5.0),
            (0.0, 4.0, 3.0),
            (4.0, 4.0, 9.0),
        ]);
        for k in [
            "exponential",
            "gaussian",
            "quartic",
            "epanechnikov",
            "polynomial5",
            "constant",
        ] {
            let (out, _) = run(json!({
                "input": pts.clone(), "z_field": "z", "cell_size": 1.0,
                "bandwidth": 10.0, "kernel": k
            }));
            assert!(
                out.get(0, 2, 2).is_finite(),
                "kernel '{k}' produced a non-finite value"
            );
        }
    }

    #[test]
    fn rejects_bad_parameters() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(KernelInterpolationWithBarriersTool.validate(&args).is_err());

        let pts = points(&[(0.0, 0.0, 1.0), (1.0, 1.0, 2.0)]);
        for bad in [
            json!({ "input": pts.clone(), "z_field": "z", "kernel": "bessel" }),
            json!({ "input": pts.clone(), "z_field": "z", "power": 3 }),
            json!({ "input": pts.clone(), "z_field": "z", "cell_size": -1 }),
            json!({ "input": pts.clone(), "z_field": "z", "bandwidth": 0 }),
            json!({ "input": pts.clone(), "z_field": "z", "ridge": -1 }),
            // NaN slips past every `<` test, so it needs its own guard.
            json!({ "input": pts.clone(), "z_field": "z", "ridge": "nan" }),
            json!({ "input": pts.clone(), "z_field": "z", "max_neighbors": "nan" }),
        ] {
            let args: ToolArgs = serde_json::from_value(bad).unwrap();
            assert!(KernelInterpolationWithBarriersTool.validate(&args).is_err());
        }

        // A missing z_field is caught at run time, when the schema is known.
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": pts, "z_field": "nope" })).unwrap();
        assert!(KernelInterpolationWithBarriersTool
            .run(&args, &ctx())
            .is_err());
    }
}
