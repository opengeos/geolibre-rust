//! GeoLibre tool: analytic 2-D Gaussian-puff contaminant transport over a
//! groundwater velocity field.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Porous Puff* (Spatial Analyst —
//! Groundwater). `darcy_flow` (#213) produces a Darcy velocity field but stops
//! at "where does the water go" — this tool answers "where does an
//! instantaneous contaminant release end up, and at what concentration".
//!
//! Given the Darcy velocity **magnitude** and **direction** rasters produced by
//! `darcy_flow`, a release point, and aquifer properties (porosity, thickness,
//! dispersivity, retardation, decay), the tool:
//!
//! 1. Reconstructs the east/north velocity components from magnitude+direction
//!    and marches the puff centre from the release point along the field,
//!    reusing `darcy_flow`'s bilinear-sampled, RK2-style step integrator —
//!    except stepping is driven by the *seepage* velocity
//!    `v = magnitude / porosity / retardation` and each step accumulates both
//!    elapsed time and travelled path length `L`, stopping exactly at the
//!    requested elapsed time `t` (extrapolating in a straight line with the
//!    last sampled velocity once the field can no longer be sampled, e.g. the
//!    plume has left the input grid, or once `max_steps` is exhausted).
//! 2. Evaluates the closed-form instantaneous 2-D Gaussian puff solution in
//!    coordinates aligned with the puff's final direction of travel:
//!
//!    ```text
//!    C(x,y,t) = M·exp(-decay·t) / (2π·σL·σT·n(x,y)·b(x,y)·R)
//!               · exp( -x'² / (2σL²) - y'² / (2σT²) )
//!    ```
//!
//!    where `x'`/`y'` are the along-flow / cross-flow offsets of the cell from
//!    the puff centre, `σL² = 2·αL·L`, `σT² = 2·αT·L` (`αL`/`αT` the
//!    longitudinal/transverse dispersivities and `L` the travelled path
//!    length), and `n(x,y)`/`b(x,y)` are the (optionally spatially varying)
//!    porosity and thickness. This is the standard instantaneous point-source
//!    solution to the 2-D depth-averaged retarded advection-dispersion
//!    equation with first-order decay; it reduces to the textbook closed form
//!    for a uniform velocity field, where `L = (magnitude/porosity/retardation)·t`
//!    and the puff centre sits at `release + L·(unit flow direction)`.
//!
//! A comma-separated `time` list produces one raster per elapsed time: the
//! last time is written to `output`, and every time (including the last) is
//! reported in the `time_series` output.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{load_input_raster, parse_optional_output, raster_like_with_data};

/// Output no-data sentinel (concentration is never negative, so this is
/// unambiguous).
const OUT_NODATA: f64 = -9999.0;

pub struct PorousPuffTool;

impl Tool for PorousPuffTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "porous_puff",
            display_name: "Porous Puff",
            summary: "Analytic advection-dispersion of an instantaneous contaminant slug through groundwater (like ArcGIS Porous Puff): march the puff centre from a release point along a Darcy velocity field (magnitude+direction, as produced by darcy_flow), then evaluate the closed-form 2-D Gaussian puff concentration in path-aligned coordinates, scaled by porosity, thickness, retardation, and first-order decay.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "magnitude",
                    description: "Darcy velocity magnitude raster, as produced by darcy_flow's 'output'. Cell size/extent are taken from this raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "direction",
                    description: "Darcy flow-direction raster (azimuth degrees from north, clockwise, down-gradient), as produced by darcy_flow's 'direction'. Must be co-registered with 'magnitude'.",
                    required: true,
                },
                ToolParamSpec {
                    name: "x",
                    description: "Release point X coordinate (map units, same CRS as the velocity rasters).",
                    required: true,
                },
                ToolParamSpec {
                    name: "y",
                    description: "Release point Y coordinate (map units).",
                    required: true,
                },
                ToolParamSpec {
                    name: "mass",
                    description: "Released contaminant mass (must be positive).",
                    required: true,
                },
                ToolParamSpec {
                    name: "porosity",
                    description: "Effective porosity: either a numeric constant (0, 1] or a raster path co-registered with 'magnitude'.",
                    required: true,
                },
                ToolParamSpec {
                    name: "thickness",
                    description: "Aquifer thickness: either a positive numeric constant or a raster path co-registered with 'magnitude'.",
                    required: true,
                },
                ToolParamSpec {
                    name: "dispersivity_long",
                    description: "Longitudinal dispersivity (length units; must be positive).",
                    required: true,
                },
                ToolParamSpec {
                    name: "dispersivity_trans",
                    description: "Transverse dispersivity (length units; must be positive).",
                    required: true,
                },
                ToolParamSpec {
                    name: "retardation",
                    description: "Retardation factor (default 1, must be >= 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "decay",
                    description: "First-order decay constant (default 0, must be >= 0; same time units as 'time').",
                    required: false,
                },
                ToolParamSpec {
                    name: "time",
                    description: "Comma-separated list of one or more positive elapsed times. The last time is written to 'output'; all times (including the last) are reported in the 'time_series' output.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output concentration raster path for the last time in 'time'. If omitted, stored in memory.",
                    required: true,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read from the magnitude/direction rasters (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "step",
                    description: "Puff-centre marching step length in CRS units (default: the magnitude raster's cell size).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_steps",
                    description: "Maximum marching steps per requested time (default 1000).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "magnitude")?;
        require_str(args, "direction")?;
        require_str(args, "output")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let mag_path = require_str(args, "magnitude")?;
        let dir_path = require_str(args, "direction")?;
        let out_path = parse_optional_output(args, "output")?;
        let prm = parse_params(args)?;

        let mag_r = load_input_raster(mag_path)?;
        let dir_r = load_input_raster(dir_path)?;
        if prm.band < 0 || prm.band as usize >= mag_r.bands {
            return Err(ToolError::Validation(format!(
                "band {} out of range (magnitude raster has {} band(s))",
                prm.band + 1,
                mag_r.bands
            )));
        }
        if dir_r.rows != mag_r.rows || dir_r.cols != mag_r.cols {
            return Err(ToolError::Execution(format!(
                "direction raster ({}x{}) is not co-registered with the magnitude raster ({}x{})",
                dir_r.rows, dir_r.cols, mag_r.rows, mag_r.cols
            )));
        }

        let rows = mag_r.rows;
        let cols = mag_r.cols;
        let cx = mag_r.cell_size_x.abs().max(f64::MIN_POSITIVE);
        let cy = mag_r.cell_size_y.abs().max(f64::MIN_POSITIVE);
        let geom = GridGeom {
            rows,
            cols,
            x_min: mag_r.x_min,
            y_max: mag_r.y_min + rows as f64 * cy,
            cx,
            cy,
        };

        // East/north velocity components reconstructed from magnitude+direction
        // (NaN marks invalid/no-data cells).
        let mnd = mag_r.nodata;
        let dnd = dir_r.nodata;
        let mut vel_e = vec![f64::NAN; rows * cols];
        let mut vel_n = vec![f64::NAN; rows * cols];
        for row in 0..rows {
            for col in 0..cols {
                let m = mag_r.get(prm.band, row as isize, col as isize);
                let d = dir_r.get(prm.band, row as isize, col as isize);
                if m == mnd || d == dnd || !m.is_finite() || !d.is_finite() || m < 0.0 {
                    continue;
                }
                let az = d.to_radians();
                let idx = row * cols + col;
                vel_e[idx] = m * az.sin();
                vel_n[idx] = m * az.cos();
            }
        }

        // Porosity / thickness grids: a constant fills every cell; a raster must
        // be co-registered with the magnitude raster.
        let porosity = prm.porosity.to_grid(&mag_r, "porosity")?;
        let thickness = prm.thickness.to_grid(&mag_r, "thickness")?;

        let field = VelocityField {
            geom: &geom,
            ve: &vel_e,
            vn: &vel_n,
        };

        ctx.progress
            .info("marching the puff centre along the velocity field");

        let step = prm.step.unwrap_or_else(|| cx.min(cy));
        let mut results = Vec::with_capacity(prm.times.len());
        for (i, &t) in prm.times.iter().enumerate() {
            let march = march_puff(
                &field,
                &porosity,
                &geom,
                prm.x0,
                prm.y0,
                t,
                prm.retardation,
                step,
                prm.max_steps,
            );
            let sigma_l = (2.0 * prm.dispersivity_long * march.length).sqrt();
            let sigma_t = (2.0 * prm.dispersivity_trans * march.length).sqrt();

            let mut data = vec![OUT_NODATA; rows * cols];
            let decay_factor = (-prm.decay * t).exp();
            let mut peak = 0.0f64;
            for row in 0..rows {
                for col in 0..cols {
                    let idx = row * cols + col;
                    let n_here = porosity[idx];
                    let b_here = thickness[idx];
                    if !n_here.is_finite() || n_here <= 0.0 || !b_here.is_finite() || b_here <= 0.0
                    {
                        continue;
                    }
                    let (cxp, cyp) = geom.cell_center(row, col);
                    let dx = cxp - march.cx;
                    let dy = cyp - march.cy;
                    // Project onto the (longitudinal, transverse) axes defined by
                    // the puff's final unit flow direction.
                    let xl = dx * march.ux + dy * march.uy;
                    let yt = -dx * march.uy + dy * march.ux;

                    let c = if sigma_l < 1e-9 || sigma_t < 1e-9 {
                        // Degenerate (near-zero travel/dispersion): all mass sits
                        // in the cell containing the puff centre.
                        if dx.hypot(dy) <= 0.5 * cx.max(cy) {
                            prm.mass * decay_factor / (prm.retardation * n_here * b_here * cx * cy)
                        } else {
                            0.0
                        }
                    } else {
                        let amp = prm.mass * decay_factor
                            / (2.0 * PI * sigma_l * sigma_t * n_here * b_here * prm.retardation);
                        amp * (-(xl * xl) / (2.0 * sigma_l * sigma_l)
                            - (yt * yt) / (2.0 * sigma_t * sigma_t))
                            .exp()
                    };
                    data[idx] = c;
                    peak = peak.max(c);
                }
            }
            ctx.progress
                .progress((i as f64 + 1.0) / prm.times.len() as f64);

            let raster = raster_like_with_data(&mag_r, data, OUT_NODATA, DataType::F32)?;
            results.push((t, march, sigma_l, sigma_t, peak, raster));
        }

        let mut outputs = BTreeMap::new();
        let mut time_series = Vec::with_capacity(results.len());
        let last_idx = results.len() - 1;
        let mut primary_out_path = String::new();
        for (i, (t, march, sigma_l, sigma_t, peak, raster)) in results.into_iter().enumerate() {
            let path = if i == last_idx {
                let p = crate::common::write_or_store_output(raster, out_path)?;
                primary_out_path = p.clone();
                p
            } else {
                let derived = out_path.map(|p| derive_time_path(p, t));
                crate::common::write_or_store_output(raster, derived.as_deref())?
            };
            time_series.push(json!({
                "time": t,
                "output": path,
                "center_x": march.cx,
                "center_y": march.cy,
                "path_length": march.length,
                "sigma_long": sigma_l,
                "sigma_trans": sigma_t,
                "peak_concentration": peak,
            }));
        }

        outputs.insert("output".to_string(), json!(primary_out_path));
        outputs.insert("time_series".to_string(), json!(time_series));

        Ok(ToolRunResult { outputs })
    }
}

/// Derives a sibling output path for a non-primary time by inserting
/// `_t<time>` before the file extension (e.g. `plume.tif` -> `plume_t5.tif`).
fn derive_time_path(base: &str, t: f64) -> String {
    let suffix = if t.fract().abs() < 1e-9 {
        format!("_t{}", t as i64)
    } else {
        format!("_t{t}")
    };
    match base.rfind('.') {
        Some(pos) if pos > base.rfind('/').unwrap_or(0) => {
            format!("{}{suffix}{}", &base[..pos], &base[pos..])
        }
        _ => format!("{base}{suffix}"),
    }
}

// ── Grid geometry + bilinear sampling ────────────────────────────────────────

/// Shared coordinate geometry for a raster grid: cell size, extent, and the
/// row-major layout. `y_max` is the raster's north edge (`y_min + rows*cy`);
/// row 0 is the north row, matching this repo's raster convention (see
/// `darcy_flow.rs`).
struct GridGeom {
    rows: usize,
    cols: usize,
    x_min: f64,
    y_max: f64,
    cx: f64,
    cy: f64,
}

impl GridGeom {
    fn cell_center(&self, row: usize, col: usize) -> (f64, f64) {
        let x = self.x_min + (col as f64 + 0.5) * self.cx;
        let y = self.y_max - (row as f64 + 0.5) * self.cy;
        (x, y)
    }

    /// Fractional (row, col) for map coordinate `(x, y)`, or `None` outside the
    /// cell-centre grid.
    fn frac(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let col_f = (x - self.x_min) / self.cx - 0.5;
        let row_f = (self.y_max - y) / self.cy - 0.5;
        if col_f < 0.0
            || row_f < 0.0
            || col_f > (self.cols - 1) as f64
            || row_f > (self.rows - 1) as f64
        {
            return None;
        }
        Some((row_f, col_f))
    }

    /// Bilinearly samples `grid` (row-major, `rows*cols`) at `(x, y)`. Returns
    /// `None` outside the grid or when any contributing cell is `NaN`.
    fn bilinear(&self, grid: &[f64], x: f64, y: f64) -> Option<f64> {
        let (row_f, col_f) = self.frac(x, y)?;
        let c0 = col_f.floor() as usize;
        let r0 = row_f.floor() as usize;
        let c1 = (c0 + 1).min(self.cols - 1);
        let r1 = (r0 + 1).min(self.rows - 1);
        let fx = col_f - c0 as f64;
        let fy = row_f - r0 as f64;
        let mut acc = 0.0;
        for (r, c, w) in [
            (r0, c0, (1.0 - fx) * (1.0 - fy)),
            (r0, c1, fx * (1.0 - fy)),
            (r1, c0, (1.0 - fx) * fy),
            (r1, c1, fx * fy),
        ] {
            let v = grid[r * self.cols + c];
            if v.is_nan() {
                return None;
            }
            acc += w * v;
        }
        Some(acc)
    }
}

/// A continuous Darcy-velocity field backed by east/north component grids
/// (`NaN` marks invalid cells), sampled bilinearly over cell centres.
struct VelocityField<'a> {
    geom: &'a GridGeom,
    ve: &'a [f64],
    vn: &'a [f64],
}

impl VelocityField<'_> {
    fn sample(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let e = self.geom.bilinear(self.ve, x, y)?;
        let n = self.geom.bilinear(self.vn, x, y)?;
        Some((e, n))
    }
}

// ── Puff-centre marching ─────────────────────────────────────────────────────

/// Result of marching the puff centre from the release point for elapsed time
/// `t`: final position, unit flow direction there, and the total travelled
/// path length (used for the dispersion spread).
struct March {
    cx: f64,
    cy: f64,
    ux: f64,
    uy: f64,
    length: f64,
}

/// Marches the puff centre from `(sx, sy)` for elapsed time `t`, stepping a
/// fixed spatial `step` along the local retarded seepage velocity
/// `magnitude/porosity/retardation`, accumulating travelled distance. Once the
/// field can no longer be sampled (e.g. the plume has left the grid), the last
/// known direction/speed is frozen and the march continues in a straight line
/// so the requested elapsed time is always reached. `max_steps` is a safety
/// cap; once exhausted the remaining time is also extrapolated in a straight
/// line at the last known speed.
#[allow(clippy::too_many_arguments)]
fn march_puff(
    field: &VelocityField,
    porosity: &[f64],
    geom: &GridGeom,
    sx: f64,
    sy: f64,
    t: f64,
    retardation: f64,
    step: f64,
    max_steps: usize,
) -> March {
    let mut x = sx;
    let mut y = sy;
    let mut length = 0.0f64;
    let mut remaining = t;
    // Default direction (north) if the field is never sampleable at all —
    // only affects the (zero-length, degenerate) output.
    let mut last_dir = (0.0f64, 1.0f64);
    let mut last_speed = 0.0f64;
    let mut have_dir = false;

    for _ in 0..max_steps {
        if remaining <= 0.0 {
            break;
        }
        let sampled = field.sample(x, y).and_then(|(ve, vn)| {
            let mag = ve.hypot(vn);
            if mag < 1e-12 {
                return None;
            }
            let n_here = geom.bilinear(porosity, x, y).unwrap_or(f64::NAN);
            if !n_here.is_finite() || n_here <= 0.0 {
                return None;
            }
            let vc = mag / n_here / retardation;
            Some(((ve / mag, vn / mag), vc))
        });
        let (dir, speed) = match sampled {
            Some((d, s)) => {
                last_dir = d;
                last_speed = s;
                have_dir = true;
                (d, s)
            }
            None if have_dir => (last_dir, last_speed),
            None => break, // never sampleable: puff does not move.
        };
        if speed < 1e-12 {
            break;
        }

        let dt_full = step / speed;
        if dt_full >= remaining {
            let d = speed * remaining;
            x += d * dir.0;
            y += d * dir.1;
            length += d;
            remaining = 0.0;
            break;
        }
        x += step * dir.0;
        y += step * dir.1;
        length += step;
        remaining -= dt_full;
    }

    // max_steps exhausted before reaching t: extrapolate the rest in a
    // straight line at the last known (or default) speed/direction.
    if remaining > 0.0 && have_dir && last_speed > 1e-12 {
        let d = last_speed * remaining;
        x += d * last_dir.0;
        y += d * last_dir.1;
        length += d;
    }

    March {
        cx: x,
        cy: y,
        ux: last_dir.0,
        uy: last_dir.1,
        length,
    }
}

// ── Porosity/thickness: raster or constant ───────────────────────────────────

enum ScalarSource {
    Const(f64),
    Path(String),
}

impl ScalarSource {
    /// Expands to a dense `rows*cols` grid matching `template`'s dimensions. A
    /// constant fills every cell; a raster is loaded and must be co-registered
    /// with `template`; its no-data cells become `NaN`.
    fn to_grid(&self, template: &wbraster::Raster, label: &str) -> Result<Vec<f64>, ToolError> {
        let rows = template.rows;
        let cols = template.cols;
        match self {
            ScalarSource::Const(v) => Ok(vec![*v; rows * cols]),
            ScalarSource::Path(p) => {
                let r = load_input_raster(p)?;
                if r.rows != rows || r.cols != cols {
                    return Err(ToolError::Execution(format!(
                        "{label} raster ({}x{}) is not co-registered with the magnitude raster ({rows}x{cols})",
                        r.rows, r.cols
                    )));
                }
                let nd = r.nodata;
                let mut g = vec![f64::NAN; rows * cols];
                for row in 0..rows {
                    for col in 0..cols {
                        let v = r.get(0, row as isize, col as isize);
                        if v != nd && v.is_finite() {
                            g[row * cols + col] = v;
                        }
                    }
                }
                Ok(g)
            }
        }
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    x0: f64,
    y0: f64,
    mass: f64,
    porosity: ScalarSource,
    thickness: ScalarSource,
    dispersivity_long: f64,
    dispersivity_trans: f64,
    retardation: f64,
    decay: f64,
    times: Vec<f64>,
    band: isize,
    step: Option<f64>,
    max_steps: usize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let x0 = require_f64(args, "x")?;
    let y0 = require_f64(args, "y")?;
    let mass = require_pos_f64(args, "mass")?;
    let porosity = require_scalar_source(args, "porosity")?;
    if let ScalarSource::Const(v) = &porosity {
        if *v <= 0.0 {
            return Err(ToolError::Validation(
                "parameter 'porosity' constant must be positive".into(),
            ));
        }
    }
    let thickness = require_scalar_source(args, "thickness")?;
    if let ScalarSource::Const(v) = &thickness {
        if *v <= 0.0 {
            return Err(ToolError::Validation(
                "parameter 'thickness' constant must be positive".into(),
            ));
        }
    }
    let dispersivity_long = require_pos_f64(args, "dispersivity_long")?;
    let dispersivity_trans = require_pos_f64(args, "dispersivity_trans")?;
    let retardation = match opt_f64(args, "retardation")? {
        Some(v) if v >= 1.0 && v.is_finite() => v,
        Some(_) => {
            return Err(ToolError::Validation(
                "parameter 'retardation' must be >= 1".into(),
            ))
        }
        None => 1.0,
    };
    let decay = match opt_f64(args, "decay")? {
        Some(v) if v >= 0.0 && v.is_finite() => v,
        Some(_) => {
            return Err(ToolError::Validation(
                "parameter 'decay' must be >= 0".into(),
            ))
        }
        None => 0.0,
    };
    let times = require_times(args, "time")?;

    let band_1based = match args.get("band") {
        None | Some(Value::Null) => 1,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(1).max(1) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'band' must be an integer".into()))?
            .max(1),
        _ => return Err(ToolError::Validation("'band' must be an integer".into())),
    };
    let step = opt_pos(args, "step")?;
    let max_steps = match args.get("max_steps") {
        None | Some(Value::Null) => 1000,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(1000).max(1) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 1000,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'max_steps' must be a positive integer".into()))?
            .max(1),
        _ => {
            return Err(ToolError::Validation(
                "'max_steps' must be a positive integer".into(),
            ))
        }
    };

    Ok(Params {
        x0,
        y0,
        mass,
        porosity,
        thickness,
        dispersivity_long,
        dispersivity_trans,
        retardation,
        decay,
        times,
        band: (band_1based - 1) as isize,
        step,
        max_steps,
    })
}

fn require_scalar_source(args: &ToolArgs, key: &str) -> Result<ScalarSource, ToolError> {
    match args.get(key) {
        Some(Value::Number(n)) => n
            .as_f64()
            .filter(|v| v.is_finite())
            .map(ScalarSource::Const)
            .ok_or_else(|| ToolError::Validation(format!("parameter '{key}' must be a number"))),
        Some(Value::String(s)) => {
            let s = s.trim();
            if s.is_empty() {
                return Err(ToolError::Validation(format!(
                    "missing required parameter '{key}'"
                )));
            }
            match s.parse::<f64>() {
                Ok(v) if v.is_finite() => Ok(ScalarSource::Const(v)),
                _ => Ok(ScalarSource::Path(s.to_string())),
            }
        }
        _ => Err(ToolError::Validation(format!(
            "missing required parameter '{key}'"
        ))),
    }
}

fn require_times(args: &ToolArgs, key: &str) -> Result<Vec<f64>, ToolError> {
    let raw = args
        .get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))?;
    let mut times = Vec::new();
    for tok in raw.split(',') {
        let tok = tok.trim();
        if tok.is_empty() {
            continue;
        }
        let v: f64 = tok
            .parse()
            .map_err(|_| ToolError::Validation(format!("invalid time value '{tok}'")))?;
        if !(v > 0.0 && v.is_finite()) {
            return Err(ToolError::Validation(format!(
                "time value '{tok}' must be a positive number"
            )));
        }
        times.push(v);
    }
    if times.is_empty() {
        return Err(ToolError::Validation(format!(
            "parameter '{key}' must list at least one positive time"
        )));
    }
    Ok(times)
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

fn require_f64(args: &ToolArgs, key: &str) -> Result<f64, ToolError> {
    opt_f64(args, key)?
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

fn require_pos_f64(args: &ToolArgs, key: &str) -> Result<f64, ToolError> {
    match opt_f64(args, key)? {
        Some(v) if v > 0.0 && v.is_finite() => Ok(v),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a positive number"
        ))),
        None => Err(ToolError::Validation(format!(
            "missing required parameter '{key}'"
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
    use wbraster::{memory_store, CrsInfo, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a single-band raster (EPSG:3857) from row-major data.
    fn raster_from(
        cols: usize,
        rows: usize,
        cell_size: f64,
        x_min: f64,
        y_min: f64,
        data: Vec<f64>,
    ) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min,
            y_min,
            cell_size,
            cell_size_y: Some(cell_size),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, data[row * cols + col])
                    .unwrap();
            }
        }
        let id = memory_store::put_raster(r);
        memory_store::make_raster_memory_path(&id)
    }

    fn constant(cols: usize, rows: usize, cell_size: f64, v: f64) -> String {
        raster_from(cols, rows, cell_size, 0.0, 0.0, vec![v; cols * rows])
    }

    fn get(r: &Raster, row: usize, col: usize) -> f64 {
        r.get(0, row as isize, col as isize)
    }

    /// Uniform eastward flow: the puff centre should advect exactly
    /// `(magnitude/porosity/retardation)*t` east of the release point, and the
    /// concentration there should match the analytic Gaussian peak.
    #[test]
    fn uniform_field_advects_and_peaks() {
        let (cols, rows) = (81, 41);
        let cell = 1.0;
        let mag = 2.0;
        let n0 = 0.5;
        let mass = 1000.0;
        let a_l = 2.0;
        let a_t = 0.4;
        let t = 5.0;

        let magnitude = constant(cols, rows, cell, mag);
        let direction = constant(cols, rows, cell, 90.0); // azimuth 90 = due east
        let porosity = constant(cols, rows, cell, n0);
        let thickness = constant(cols, rows, cell, 1.0);

        // Release near the west edge so the plume stays inside the grid.
        let sx = 5.0;
        let sy = rows as f64 / 2.0;

        let args: ToolArgs = serde_json::from_value(json!({
            "magnitude": magnitude, "direction": direction,
            "x": sx, "y": sy, "mass": mass,
            "porosity": porosity, "thickness": thickness,
            "dispersivity_long": a_l, "dispersivity_trans": a_t,
            "time": t.to_string(),
        }))
        .unwrap();
        let out = PorousPuffTool.run(&args, &ctx()).unwrap();
        let raster = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();

        let v = mag / n0; // retardation = 1
        let expected_dx = v * t;
        let center_x = sx + expected_dx;
        let center_y = sy;

        // Analytic peak amplitude at the puff centre.
        let length = expected_dx; // uniform field: path length == straight-line distance
        let sigma_l = (2.0 * a_l * length).sqrt();
        let sigma_t = (2.0 * a_t * length).sqrt();
        let expected_peak =
            mass / (2.0 * PI * sigma_l * sigma_t * n0 * thickness_const(&thickness));

        // Sample the raster cell nearest the analytic centre.
        let col = ((center_x - 0.0) / cell - 0.5).round().max(0.0) as usize;
        let row = (rows as f64 - (center_y - 0.0) / cell - 0.5)
            .round()
            .max(0.0) as usize;
        let sampled = get(&raster, row.min(rows - 1), col.min(cols - 1));

        assert!(
            (sampled - expected_peak).abs() / expected_peak < 0.1,
            "sampled {sampled} vs analytic peak {expected_peak}"
        );

        // The reported centre/path length match the closed form too.
        let ts = out.outputs["time_series"].as_array().unwrap();
        assert_eq!(ts.len(), 1);
        let center_x_out = ts[0]["center_x"].as_f64().unwrap();
        let center_y_out = ts[0]["center_y"].as_f64().unwrap();
        assert!((center_x_out - center_x).abs() < 1e-6);
        assert!((center_y_out - center_y).abs() < 1e-6);
        let length_out = ts[0]["path_length"].as_f64().unwrap();
        assert!((length_out - length).abs() < 1e-6);
    }

    fn thickness_const(path: &str) -> f64 {
        let r = load_input_raster(path).unwrap();
        get(&r, 0, 0)
    }

    /// The discretized mass integral (sum of C·n·b·cell_area) approximately
    /// conserves the released mass when decay is zero and retardation is 1.
    #[test]
    fn mass_is_approximately_conserved() {
        let (cols, rows) = (81, 61);
        let cell = 1.0;
        let mag = 1.5;
        let n0 = 0.3;
        let b0 = 2.0;
        let mass = 500.0;
        let a_l = 1.5;
        let a_t = 0.5;
        let t = 8.0;

        let magnitude = constant(cols, rows, cell, mag);
        let direction = constant(cols, rows, cell, 90.0);
        let porosity = constant(cols, rows, cell, n0);
        let thickness = constant(cols, rows, cell, b0);

        let sx = 10.0;
        let sy = rows as f64 / 2.0;

        let args: ToolArgs = serde_json::from_value(json!({
            "magnitude": magnitude, "direction": direction,
            "x": sx, "y": sy, "mass": mass,
            "porosity": porosity, "thickness": thickness,
            "dispersivity_long": a_l, "dispersivity_trans": a_t,
            "time": t.to_string(),
        }))
        .unwrap();
        let out = PorousPuffTool.run(&args, &ctx()).unwrap();
        let raster = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();

        let mut total = 0.0;
        for row in 0..rows {
            for col in 0..cols {
                let c = get(&raster, row, col);
                if c > OUT_NODATA + 1.0 {
                    total += c * n0 * b0 * cell * cell;
                }
            }
        }
        let ratio = total / mass;
        assert!(
            (ratio - 1.0).abs() < 0.05,
            "mass ratio {ratio} should be close to 1.0 (grid extent may be clipping the plume)"
        );
    }

    /// Positive decay reduces the peak concentration by exactly `exp(-decay*t)`
    /// relative to zero decay, all else equal.
    #[test]
    fn decay_reduces_peak_by_exp_factor() {
        let (cols, rows) = (41, 31);
        let cell = 1.0;
        let magnitude = constant(cols, rows, cell, 1.0);
        let direction = constant(cols, rows, cell, 90.0);
        let porosity = constant(cols, rows, cell, 0.3);
        let thickness = constant(cols, rows, cell, 1.0);
        let t: f64 = 6.0;
        let decay: f64 = 0.05;

        let base = json!({
            "magnitude": magnitude, "direction": direction,
            "x": 5.0, "y": (rows as f64) / 2.0, "mass": 100.0,
            "porosity": porosity, "thickness": thickness,
            "dispersivity_long": 1.0, "dispersivity_trans": 0.3,
            "time": t.to_string(),
        });

        let mut zero_decay = base.clone();
        zero_decay["decay"] = json!(0.0);
        let args: ToolArgs = serde_json::from_value(zero_decay).unwrap();
        let out0 = PorousPuffTool.run(&args, &ctx()).unwrap();
        let peak0 = out0.outputs["time_series"][0]["peak_concentration"]
            .as_f64()
            .unwrap();

        let mut with_decay = base;
        with_decay["decay"] = json!(decay);
        let args: ToolArgs = serde_json::from_value(with_decay).unwrap();
        let out1 = PorousPuffTool.run(&args, &ctx()).unwrap();
        let peak1 = out1.outputs["time_series"][0]["peak_concentration"]
            .as_f64()
            .unwrap();

        let expected = peak0 * (-decay * t).exp();
        assert!(
            (peak1 - expected).abs() / expected < 1e-6,
            "decayed peak {peak1} != expected {expected}"
        );
    }

    #[test]
    fn rejects_bad_parameters() {
        let (cols, rows) = (10, 10);
        let magnitude = constant(cols, rows, 1.0, 1.0);
        let direction = constant(cols, rows, 1.0, 90.0);

        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            PorousPuffTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "magnitude": magnitude, "direction": direction })).is_err()); // missing rest
        assert!(bad(json!({
            "magnitude": magnitude, "direction": direction,
            "x": 1.0, "y": 1.0, "mass": -5.0,
            "porosity": 0.3, "thickness": 1.0,
            "dispersivity_long": 1.0, "dispersivity_trans": 0.3,
            "time": "5", "output": "out.tif",
        }))
        .is_err()); // negative mass
        assert!(bad(json!({
            "magnitude": magnitude, "direction": direction,
            "x": 1.0, "y": 1.0, "mass": 5.0,
            "porosity": 0.0, "thickness": 1.0,
            "dispersivity_long": 1.0, "dispersivity_trans": 0.3,
            "time": "5", "output": "out.tif",
        }))
        .is_err()); // non-positive porosity constant
        assert!(bad(json!({
            "magnitude": magnitude, "direction": direction,
            "x": 1.0, "y": 1.0, "mass": 5.0,
            "porosity": 0.3, "thickness": 1.0,
            "dispersivity_long": -1.0, "dispersivity_trans": 0.3,
            "time": "5", "output": "out.tif",
        }))
        .is_err()); // negative dispersivity
        assert!(bad(json!({
            "magnitude": magnitude, "direction": direction,
            "x": 1.0, "y": 1.0, "mass": 5.0,
            "porosity": 0.3, "thickness": 1.0,
            "dispersivity_long": 1.0, "dispersivity_trans": 0.3,
            "time": "0", "output": "out.tif",
        }))
        .is_err()); // non-positive time
        assert!(bad(json!({
            "magnitude": magnitude, "direction": direction,
            "x": 1.0, "y": 1.0, "mass": 5.0,
            "porosity": 0.3, "thickness": 1.0,
            "dispersivity_long": 1.0, "dispersivity_trans": 0.3,
            "time": "5", "output": "out.tif",
        }))
        .is_ok());
    }
}
