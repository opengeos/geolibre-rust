//! GeoLibre tool: `topo_to_raster` — hydrologically aware DEM interpolation
//! from contour lines, spot-height points, and (optionally) stream lines.
//!
//! Pure-Rust, dependency-free counterpart of ArcGIS Pro's *Topo To Raster*
//! (the ANUDEM method of Hutchinson). None of the ~791 bundled tools nor the
//! generic interpolators (IDW / TIN / spline / kriging) are drainage aware:
//! interpolating contours with a generic method terraces the surface and leaves
//! closed depressions between adjacent contour rings. `topo_to_raster` instead
//! fits a smooth **minimum-curvature (thin-plate) surface** that honours the
//! input elevations as fixed constraints, then optionally enforces drainage by
//! removing spurious interior sinks.
//!
//! ## Method (single-resolution finite-difference relaxation)
//! 1. Every contour segment is sampled onto the output grid and every spot
//!    height is dropped onto its cell; those cells become fixed Dirichlet
//!    constraints holding the known elevation.
//! 2. The remaining cells are solved by successive over-relaxation (SOR) of a
//!    discrete surface-fitting operator that blends the biharmonic (thin-plate,
//!    `∇⁴z = 0` → smooth, no terracing) and the Laplacian (membrane,
//!    `∇²z = 0`) via a `tension` weight. Iteration stops on a max-change
//!    tolerance or an iteration cap.
//! 3. With `enforce_drainage` (default on) the converged surface is passed
//!    through Wang & Liu priority-flood filling to remove closed interior
//!    depressions, so the DEM drains monotonically to its edges. Optional
//!    `streams` are burned a small amount into the surface first so flow is
//!    routed along them.
//!
//! ## Deliberate v1 scope cuts (documented for reviewers)
//! - Single-resolution SOR, not the coarse-to-fine multigrid of full ANUDEM.
//!   Convergence is accelerated with the membrane (tension) term rather than a
//!   grid hierarchy.
//! - Drainage enforcement fills every closed depression rather than ANUDEM's
//!   selective, morphology-preserving sink handling. Stream carving is a simple
//!   pre-fill burn, not full monotonic-downhill enforcement along each line.
//! - Boundary polygons and a curvature/roughness auto-estimator are not
//!   implemented; the analysis extent is the union of the input data.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{memory_store, CrsInfo, DataType, Raster, RasterConfig, RasterFormat};
use wbvector::{Coord, Geometry, Layer};

use crate::fill::fill_depressions_wang_and_liu;
use crate::vector_common::{load_input_layer, parse_optional_str};

pub struct TopoToRasterTool;

impl Tool for TopoToRasterTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "topo_to_raster",
            display_name: "Topo To Raster",
            summary: "Interpolate a hydrologically correct DEM from contour lines, spot-height points, and optional stream lines using a minimum-curvature (thin-plate) relaxation with drainage enforcement — like ArcGIS Topo To Raster (ANUDEM). Unlike the bundled IDW/TIN/spline interpolators it does not terrace contours and it removes spurious closed depressions.",
            category: ToolCategory::Hydrology,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "contours",
                    description: "Input contour line vector layer (elevation read from elevation_field). At least one of contours/points is required.",
                    required: false,
                },
                ToolParamSpec {
                    name: "points",
                    description: "Input spot-height point vector layer (elevation read from elevation_field).",
                    required: false,
                },
                ToolParamSpec {
                    name: "streams",
                    description: "Optional stream line vector layer used for drainage enforcement (burned into the surface before sink removal).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output DEM raster path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "elevation_field",
                    description: "Name of the numeric elevation field on the contour and point layers. Default: auto-detect (HEIGHT/ELEV/ELEVATION/Z, else the first numeric field).",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Output cell size in CRS units. Default: sized so the longer grid axis is ~400 cells.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tension",
                    description: "Surface stiffness in [0,1]: 0 = pure thin-plate (smoothest), 1 = taut membrane. Default 0.35.",
                    required: false,
                },
                ToolParamSpec {
                    name: "iterations",
                    description: "Maximum SOR iterations (default 2000). Iteration stops early once the max cell change falls below tolerance.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Convergence tolerance: stop when the largest cell change in an iteration is below this (CRS elevation units). Default 0.01.",
                    required: false,
                },
                ToolParamSpec {
                    name: "enforce_drainage",
                    description: "Remove spurious interior sinks after interpolation so the DEM drains to its edges (default true).",
                    required: false,
                },
                ToolParamSpec {
                    name: "stream_burn",
                    description: "Depth (elevation units) to lower stream cells before sink removal, to route flow along them. Default 0 (auto = 1% of relief when streams are supplied).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        let prm = parse_params(args)?;
        if prm.contours.is_none() && prm.points.is_none() {
            return Err(ToolError::Validation(
                "at least one of 'contours' or 'points' must be provided".to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let prm = parse_params(args)?;
        let output = parse_optional_str(args, "output")?;

        if prm.contours.is_none() && prm.points.is_none() {
            return Err(ToolError::Validation(
                "at least one of 'contours' or 'points' must be provided".to_string(),
            ));
        }

        // ── Collect elevation constraints from contours and points ───────────
        let mut samples: Vec<Sample> = Vec::new();
        let mut epsg: Option<u32> = None;
        let mut n_contour_features = 0usize;
        let mut n_point_features = 0usize;

        if let Some(path) = &prm.contours {
            let layer = load_input_layer(path)?;
            epsg = epsg.or_else(|| layer.crs_epsg());
            let field = resolve_elevation_field(&layer, prm.elevation_field.as_deref())?;
            for feature in layer.iter() {
                let Some(z) = feature.attributes.get(field).and_then(|v| v.as_f64()) else {
                    continue;
                };
                let Some(geom) = feature.geometry.as_ref() else {
                    continue;
                };
                let before = samples.len();
                collect_line_samples(geom, z, &mut samples);
                if samples.len() > before {
                    n_contour_features += 1;
                }
            }
        }

        if let Some(path) = &prm.points {
            let layer = load_input_layer(path)?;
            epsg = epsg.or_else(|| layer.crs_epsg());
            let field = resolve_elevation_field(&layer, prm.elevation_field.as_deref())?;
            for feature in layer.iter() {
                let Some(z) = feature.attributes.get(field).and_then(|v| v.as_f64()) else {
                    continue;
                };
                let Some(geom) = feature.geometry.as_ref() else {
                    continue;
                };
                collect_point_samples(geom, z, &mut samples);
                n_point_features += 1;
            }
        }

        if samples.len() < 2 {
            return Err(ToolError::Execution(
                "need at least two valid elevation samples to interpolate a surface".to_string(),
            ));
        }

        // ── Establish the output grid from the data extent ───────────────────
        let (mut min_x, mut min_y, mut max_x, mut max_y) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for s in &samples {
            min_x = min_x.min(s.x);
            min_y = min_y.min(s.y);
            max_x = max_x.max(s.x);
            max_y = max_y.max(s.y);
        }
        let width = (max_x - min_x).max(f64::EPSILON);
        let height = (max_y - min_y).max(f64::EPSILON);

        let cell_size = match prm.cell_size {
            Some(cs) => cs,
            None => (width.max(height) / 400.0).max(f64::MIN_POSITIVE),
        };
        // One-cell margin so edge contours sit inside the grid.
        let margin = cell_size;
        let x_min = min_x - margin;
        let y_min = min_y - margin;
        let cols = (((max_x + margin - x_min) / cell_size).ceil() as usize).max(2);
        let rows = (((max_y + margin - y_min) / cell_size).ceil() as usize).max(2);

        if rows.saturating_mul(cols) > 40_000_000 {
            return Err(ToolError::Validation(format!(
                "requested grid is {rows}x{cols} cells (> 40M); increase cell_size"
            )));
        }
        let y_max = y_min + rows as f64 * cell_size;

        ctx.progress.info(&format!(
            "interpolating a {rows}x{cols} DEM from {} sample(s)",
            samples.len()
        ));

        // ── Rasterize constraints (mean elevation per fixed cell) ────────────
        let n = rows * cols;
        let mut sum = vec![0.0_f64; n];
        let mut cnt = vec![0u32; n];
        let mut mean_all = 0.0_f64;
        for s in &samples {
            let col = ((s.x - x_min) / cell_size).floor() as isize;
            // Row 0 is the north (top) edge: invert y.
            let row = ((y_max - s.y) / cell_size).floor() as isize;
            if row < 0 || col < 0 || row as usize >= rows || col as usize >= cols {
                continue;
            }
            let i = row as usize * cols + col as usize;
            sum[i] += s.z;
            cnt[i] += 1;
            mean_all += s.z;
        }
        mean_all /= samples.len() as f64;

        let mut fixed = vec![false; n];
        let mut z = vec![mean_all; n];
        let mut n_fixed = 0usize;
        for i in 0..n {
            if cnt[i] > 0 {
                z[i] = sum[i] / cnt[i] as f64;
                fixed[i] = true;
                n_fixed += 1;
            }
        }
        if n_fixed == 0 {
            return Err(ToolError::Execution(
                "no elevation sample fell inside the output grid".to_string(),
            ));
        }

        // ── Solve the surface by SOR relaxation ──────────────────────────────
        let (iters_run, final_delta) = relax_surface(
            &mut z,
            &fixed,
            rows,
            cols,
            prm.tension,
            prm.iterations,
            prm.tolerance,
            ctx,
        );
        ctx.progress.info(&format!(
            "relaxation stopped after {iters_run} iteration(s), max change {final_delta:.4}"
        ));

        // ── Optional stream burn + drainage enforcement ─────────────────────
        let mut sinks_removed = 0usize;
        if prm.enforce_drainage {
            let (mut z_min, mut z_max) = (f64::INFINITY, f64::NEG_INFINITY);
            for &v in &z {
                z_min = z_min.min(v);
                z_max = z_max.max(v);
            }
            let relief = (z_max - z_min).max(1e-9);

            if let Some(path) = &prm.streams {
                let burn = if prm.stream_burn > 0.0 {
                    prm.stream_burn
                } else {
                    0.01 * relief
                };
                let layer = load_input_layer(path)?;
                let mut stream_samples: Vec<Sample> = Vec::new();
                for feature in layer.iter() {
                    if let Some(geom) = feature.geometry.as_ref() {
                        collect_line_samples(geom, 0.0, &mut stream_samples);
                    }
                }
                for s in &stream_samples {
                    let col = ((s.x - x_min) / cell_size).floor() as isize;
                    let row = ((y_max - s.y) / cell_size).floor() as isize;
                    if row < 0 || col < 0 || row as usize >= rows || col as usize >= cols {
                        continue;
                    }
                    z[row as usize * cols + col as usize] -= burn;
                }
                ctx.progress.info(&format!(
                    "burned {} stream cell(s) by {burn:.3}",
                    stream_samples.len()
                ));
            }

            let nodata = f64::NAN;
            let filled = fill_depressions_wang_and_liu(&z, rows, cols, nodata, 0.0);
            for i in 0..n {
                if filled[i] > z[i] + 1e-6 {
                    sinks_removed += 1;
                }
                z[i] = filled[i];
            }
            ctx.progress.info(&format!(
                "drainage enforcement raised {sinks_removed} sink cell(s)"
            ));
        }

        // ── Assemble and write the output raster ────────────────────────────
        let (mut z_min, mut z_max) = (f64::INFINITY, f64::NEG_INFINITY);
        for &v in &z {
            z_min = z_min.min(v);
            z_max = z_max.max(v);
        }
        let nodata = -9999.0_f64;
        let raster = build_raster(&z, rows, cols, x_min, y_min, cell_size, nodata, epsg)?;
        let out_path = write_or_store_output(raster, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        outputs.insert("cell_size".to_string(), json!(cell_size));
        outputs.insert("fixed_cells".to_string(), json!(n_fixed));
        outputs.insert("contour_features".to_string(), json!(n_contour_features));
        outputs.insert("point_features".to_string(), json!(n_point_features));
        outputs.insert("iterations".to_string(), json!(iters_run));
        outputs.insert("max_change".to_string(), json!(final_delta));
        outputs.insert("sinks_removed".to_string(), json!(sinks_removed));
        outputs.insert("z_min".to_string(), json!(z_min));
        outputs.insert("z_max".to_string(), json!(z_max));
        Ok(ToolRunResult { outputs })
    }
}

// ── Surface relaxation ───────────────────────────────────────────────────────

/// Successive over-relaxation of a blended biharmonic (thin-plate) / Laplacian
/// (membrane) surface operator. Fixed cells hold their constraint value.
///
/// Returns `(iterations_run, final_max_change)`.
#[allow(clippy::too_many_arguments)]
fn relax_surface(
    z: &mut [f64],
    fixed: &[bool],
    rows: usize,
    cols: usize,
    tension: f64,
    max_iters: usize,
    tolerance: f64,
    ctx: &ToolContext,
) -> (usize, f64) {
    const OMEGA: f64 = 1.5; // over-relaxation factor
    let t = tension.clamp(0.0, 1.0);
    let idx = |r: usize, c: usize| r * cols + c;
    let mut last_delta = 0.0;

    for iter in 0..max_iters {
        let mut max_delta = 0.0_f64;
        for r in 0..rows {
            for c in 0..cols {
                let i = idx(r, c);
                if fixed[i] {
                    continue;
                }
                let membrane = membrane_target(z, r, c, rows, cols, idx);
                // Blend in the biharmonic (thin-plate) target only in the deep
                // interior where the 13-point stencil is fully available.
                let target = if t >= 1.0 || r < 2 || c < 2 || r + 2 >= rows || c + 2 >= cols {
                    membrane
                } else {
                    let plate = biharmonic_target(z, r, c, idx);
                    t * membrane + (1.0 - t) * plate
                };
                let old = z[i];
                let new = old + OMEGA * (target - old);
                z[i] = new;
                let d = (new - old).abs();
                if d > max_delta {
                    max_delta = d;
                }
            }
        }
        last_delta = max_delta;
        if iter % 64 == 0 {
            ctx.progress
                .progress((iter as f64 + 1.0) / max_iters as f64);
        }
        if max_delta < tolerance {
            return (iter + 1, max_delta);
        }
    }
    (max_iters, last_delta)
}

/// Membrane (Laplacian) target: average of the four orthogonal neighbours,
/// with edges reflected (Neumann boundary) by substituting the centre value.
#[inline]
fn membrane_target(
    z: &[f64],
    r: usize,
    c: usize,
    rows: usize,
    cols: usize,
    idx: impl Fn(usize, usize) -> usize,
) -> f64 {
    let center = z[idx(r, c)];
    let n = if r > 0 { z[idx(r - 1, c)] } else { center };
    let s = if r + 1 < rows {
        z[idx(r + 1, c)]
    } else {
        center
    };
    let w = if c > 0 { z[idx(r, c - 1)] } else { center };
    let e = if c + 1 < cols {
        z[idx(r, c + 1)]
    } else {
        center
    };
    0.25 * (n + s + w + e)
}

/// Biharmonic (thin-plate) target solving `∇⁴z = 0` for the centre cell using
/// the 13-point stencil. Only valid when `2 <= r < rows-2` and `2 <= c < cols-2`.
#[inline]
fn biharmonic_target(z: &[f64], r: usize, c: usize, idx: impl Fn(usize, usize) -> usize) -> f64 {
    let o = z[idx(r - 1, c)] + z[idx(r + 1, c)] + z[idx(r, c - 1)] + z[idx(r, c + 1)];
    let d =
        z[idx(r - 1, c - 1)] + z[idx(r - 1, c + 1)] + z[idx(r + 1, c - 1)] + z[idx(r + 1, c + 1)];
    let f = z[idx(r - 2, c)] + z[idx(r + 2, c)] + z[idx(r, c - 2)] + z[idx(r, c + 2)];
    (8.0 * o - 2.0 * d - f) / 20.0
}

// ── Geometry sampling ────────────────────────────────────────────────────────

struct Sample {
    x: f64,
    y: f64,
    z: f64,
}

/// Densely samples a line geometry (Line/MultiLine) at ~1 unit per segment step
/// (the caller's grid resolution is finer than any real gap once rasterized) so
/// every crossed cell picks up the contour elevation. Steps are in coordinate
/// units; a fixed 24-sub-step per segment plus the endpoints keeps short and
/// long segments both represented without needing the cell size here.
fn collect_line_samples(geom: &Geometry, z: f64, out: &mut Vec<Sample>) {
    let mut push_line = |coords: &[Coord]| {
        for w in coords.windows(2) {
            let (a, b) = (&w[0], &w[1]);
            out.push(Sample { x: a.x, y: a.y, z });
            let dx = b.x - a.x;
            let dy = b.y - a.y;
            let steps = 24usize;
            for k in 1..steps {
                let t = k as f64 / steps as f64;
                out.push(Sample {
                    x: a.x + dx * t,
                    y: a.y + dy * t,
                    z,
                });
            }
        }
        if let Some(last) = coords.last() {
            out.push(Sample {
                x: last.x,
                y: last.y,
                z,
            });
        }
    };
    match geom {
        Geometry::LineString(cs) => push_line(cs),
        Geometry::MultiLineString(lines) => {
            for l in lines {
                push_line(l);
            }
        }
        Geometry::Point(c) => out.push(Sample { x: c.x, y: c.y, z }),
        Geometry::MultiPoint(cs) => {
            for c in cs {
                out.push(Sample { x: c.x, y: c.y, z });
            }
        }
        _ => {}
    }
}

fn collect_point_samples(geom: &Geometry, z: f64, out: &mut Vec<Sample>) {
    match geom {
        Geometry::Point(c) => out.push(Sample { x: c.x, y: c.y, z }),
        Geometry::MultiPoint(cs) => {
            for c in cs {
                out.push(Sample { x: c.x, y: c.y, z });
            }
        }
        _ => {}
    }
}

// ── Field resolution ─────────────────────────────────────────────────────────

/// Resolves the elevation field index: the caller-named field if given,
/// otherwise a common elevation name, otherwise the first numeric field.
fn resolve_elevation_field(layer: &Layer, name: Option<&str>) -> Result<usize, ToolError> {
    if let Some(name) = name {
        return layer
            .schema
            .field_index(name)
            .ok_or_else(|| ToolError::Validation(format!("elevation_field '{name}' not found")));
    }
    for candidate in [
        "HEIGHT",
        "ELEV",
        "ELEVATION",
        "Z",
        "height",
        "elev",
        "elevation",
        "z",
    ] {
        if let Some(i) = layer.schema.field_index(candidate) {
            return Ok(i);
        }
    }
    use wbvector::FieldType;
    for (i, f) in layer.schema.fields().iter().enumerate() {
        if matches!(f.field_type, FieldType::Integer | FieldType::Float) {
            return Ok(i);
        }
    }
    Err(ToolError::Validation(
        "no elevation field found; specify 'elevation_field'".to_string(),
    ))
}

// ── Raster output ────────────────────────────────────────────────────────────

#[allow(clippy::too_many_arguments)]
fn build_raster(
    z: &[f64],
    rows: usize,
    cols: usize,
    x_min: f64,
    y_min: f64,
    cell_size: f64,
    nodata: f64,
    epsg: Option<u32>,
) -> Result<Raster, ToolError> {
    let crs = match epsg {
        Some(code) => CrsInfo::from_epsg(code),
        None => CrsInfo::new(),
    };
    let mut raster = Raster::new(RasterConfig {
        cols,
        rows,
        bands: 1,
        x_min,
        y_min,
        cell_size,
        cell_size_y: Some(cell_size),
        nodata,
        data_type: DataType::F32,
        crs,
        metadata: Vec::new(),
    });
    for r in 0..rows {
        for c in 0..cols {
            raster
                .set(0, r as isize, c as isize, z[r * cols + c])
                .map_err(|e| ToolError::Execution(format!("failed writing cell: {e}")))?;
        }
    }
    Ok(raster)
}

fn write_or_store_output(raster: Raster, output_path: Option<&str>) -> Result<String, ToolError> {
    match output_path {
        Some(output_path) => {
            if let Some(parent) = std::path::Path::new(output_path).parent() {
                if !parent.as_os_str().is_empty() {
                    std::fs::create_dir_all(parent).map_err(|e| {
                        ToolError::Execution(format!("failed creating output directory: {e}"))
                    })?;
                }
            }
            let fmt = RasterFormat::for_output_path(output_path)
                .map_err(|e| ToolError::Validation(format!("unsupported output path: {e}")))?;
            raster
                .write(output_path, fmt)
                .map_err(|e| ToolError::Execution(format!("failed writing output raster: {e}")))?;
            Ok(output_path.to_string())
        }
        None => {
            let id = memory_store::put_raster(raster);
            Ok(memory_store::make_raster_memory_path(&id))
        }
    }
}

// ── Parameters ───────────────────────────────────────────────────────────────

struct Params {
    contours: Option<String>,
    points: Option<String>,
    streams: Option<String>,
    elevation_field: Option<String>,
    cell_size: Option<f64>,
    tension: f64,
    iterations: usize,
    tolerance: f64,
    enforce_drainage: bool,
    stream_burn: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let contours = parse_optional_str(args, "contours")?.map(str::to_string);
    let points = parse_optional_str(args, "points")?.map(str::to_string);
    let streams = parse_optional_str(args, "streams")?.map(str::to_string);
    let elevation_field = parse_optional_str(args, "elevation_field")?.map(str::to_string);
    let cell_size = opt_pos(args, "cell_size")?;
    let tension = opt_f64(args, "tension")?.unwrap_or(0.35);
    if !(0.0..=1.0).contains(&tension) {
        return Err(ToolError::Validation(
            "parameter 'tension' must be in [0, 1]".to_string(),
        ));
    }
    let iterations = opt_f64(args, "iterations")?
        .map(|v| v.max(1.0) as usize)
        .unwrap_or(2000);
    let tolerance = opt_pos(args, "tolerance")?.unwrap_or(0.01);
    let enforce_drainage = opt_bool(args, "enforce_drainage")?.unwrap_or(true);
    let stream_burn = opt_f64(args, "stream_burn")?.unwrap_or(0.0).max(0.0);
    Ok(Params {
        contours,
        points,
        streams,
        elevation_field,
        cell_size,
        tension,
        iterations,
        tolerance,
        enforce_drainage,
        stream_burn,
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

fn opt_pos(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match opt_f64(args, key)? {
        Some(v) if v > 0.0 && v.is_finite() => Ok(Some(v)),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a positive number"
        ))),
        None => Ok(None),
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
    use wbvector::{FieldDef, FieldType, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a contour line layer: each row is (z, [(x,y), ...]).
    fn contour_layer(lines: &[(f64, Vec<(f64, f64)>)]) -> String {
        let mut l = Layer::new("contours")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("HEIGHT", FieldType::Float));
        for (z, pts) in lines {
            let coords: Vec<Coord> = pts.iter().map(|(x, y)| Coord::xy(*x, *y)).collect();
            l.add_feature(
                Some(Geometry::line_string(coords)),
                &[("HEIGHT", (*z).into())],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn point_layer(pts: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("Z", FieldType::Float));
        for (x, y, z) in pts {
            l.add_feature(Some(Geometry::point(*x, *y)), &[("Z", (*z).into())])
                .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = TopoToRasterTool.run(&args, &ctx()).unwrap();
        let raster =
            crate::common::load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, raster)
    }

    /// A planar surface z = 2x is harmonic and biharmonic, so between contour
    /// constraints the interior must be reconstructed as a plane — no terracing.
    /// We verify the along-x gradient is a constant ~2 across the interior
    /// (rather than the flat steps a naive contour rasterization would give).
    #[test]
    fn reconstructs_a_plane_between_contours() {
        // Vertical contour lines at x = 0,10,...,100 with z = 2*x.
        let mut lines = Vec::new();
        for xi in 0..=10 {
            let x = xi as f64 * 10.0;
            lines.push((2.0 * x, vec![(x, 0.0), (x, 100.0)]));
        }
        let path = contour_layer(&lines);
        let (_out, raster) = run(json!({
            "contours": path,
            "cell_size": 5.0,
            "tension": 0.0,
            "iterations": 6000,
            "tolerance": 1e-6,
            "enforce_drainage": false,
        }));
        // Read a mid-row along x; the finite-difference gradient must be ~2
        // everywhere the stencil is fully interior (no flat "terrace" steps).
        let row = raster.rows as isize / 2;
        let cell = raster.cell_size_x;
        let mut prev: Option<f64> = None;
        let mut checked = 0;
        for col in 4..(raster.cols as isize - 4) {
            let v = raster.get(0, row, col);
            if let Some(p) = prev {
                let grad = (v - p) / cell;
                assert!(
                    (grad - 2.0).abs() < 0.25,
                    "interior gradient {grad} should be ~2 (plane, no terracing) at col {col}"
                );
                checked += 1;
            }
            prev = Some(v);
        }
        assert!(checked > 5, "expected several interior gradient checks");
    }

    /// Points-only input interpolates a bowl and the centre lies between the
    /// low centre point and the high rim.
    #[test]
    fn interpolates_from_points() {
        let pts = [
            (0.0, 0.0, 100.0),
            (100.0, 0.0, 100.0),
            (0.0, 100.0, 100.0),
            (100.0, 100.0, 100.0),
            (50.0, 50.0, 0.0), // low centre
        ];
        let path = point_layer(&pts);
        let (_out, raster) = run(json!({
            "points": path,
            "cell_size": 5.0,
            "enforce_drainage": false,
        }));
        let col = ((50.0 - raster.x_min) / raster.cell_size_x).floor() as isize;
        let row = ((raster.y_max() - 50.0) / raster.cell_size_y).floor() as isize;
        let center = raster.get(0, row, col);
        let corner = raster.get(0, 1, 1);
        assert!(center < corner, "bowl centre {center} below rim {corner}");
        assert!(center >= -1.0, "centre should not undershoot far below 0");
    }

    /// Drainage enforcement removes an interior sink: after enforcement no
    /// interior cell is strictly lower than all of its 4 neighbours.
    #[test]
    fn drainage_removes_interior_sinks() {
        // A ring of high contours around a single deep low point creates a pit.
        // Outer square contour at z = 50.
        let lines = vec![(
            50.0,
            vec![
                (0.0, 0.0),
                (100.0, 0.0),
                (100.0, 100.0),
                (0.0, 100.0),
                (0.0, 0.0),
            ],
        )];
        let cpath = contour_layer(&lines);
        let ppath = point_layer(&[(50.0, 50.0, 0.0)]); // deep central pit
        let (out, raster) = run(json!({
            "contours": cpath,
            "points": ppath,
            "cell_size": 5.0,
            "enforce_drainage": true,
        }));
        assert!(
            out.outputs["sinks_removed"].as_u64().unwrap() > 0,
            "the central pit should be filled"
        );
        // Verify no strict interior sink remains.
        let rows = raster.rows as isize;
        let cols = raster.cols as isize;
        let mut strict_sinks = 0;
        for r in 1..rows - 1 {
            for c in 1..cols - 1 {
                let v = raster.get(0, r, c);
                let mut is_sink = true;
                for (dr, dc) in [(-1, 0), (1, 0), (0, -1), (0, 1)] {
                    if raster.get(0, r + dr, c + dc) <= v {
                        is_sink = false;
                        break;
                    }
                }
                if is_sink {
                    strict_sinks += 1;
                }
            }
        }
        assert_eq!(
            strict_sinks, 0,
            "no interior sink should remain after enforcement"
        );
    }

    #[test]
    fn rejects_no_input() {
        let args: ToolArgs = serde_json::from_value(json!({ "cell_size": 5.0 })).unwrap();
        assert!(TopoToRasterTool.validate(&args).is_err());
        assert!(TopoToRasterTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn rejects_bad_tension() {
        let path = point_layer(&[(0.0, 0.0, 1.0), (10.0, 10.0, 2.0)]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "points": path, "tension": 5.0 })).unwrap();
        assert!(TopoToRasterTool.validate(&args).is_err());
    }
}
