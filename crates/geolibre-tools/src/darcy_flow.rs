//! GeoLibre tool: steady-state groundwater Darcy velocity + particle tracking.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Darcy Flow* and *Particle Track*
//! (Spatial Analyst — Groundwater). The ~791 bundled whitebox IDs contain zero
//! groundwater/hydrogeology tools, so a Darcy velocity field opens a niche no
//! WASM tool covers while reusing the repo's DEM/raster machinery.
//!
//! Given a groundwater **head** raster `h`, a **transmissivity** raster `T`, and
//! an effective-**porosity** raster `n` (all co-registered), the tool applies
//! Darcy's law cell by cell:
//!
//! ```text
//!   v = -(T / n) · ∇h
//! ```
//!
//! The head gradient `∇h = (∂h/∂East, ∂h/∂North)` is a 3×3 Horn finite-difference
//! gradient (exact for a planar surface, matching the repo's slope/aspect
//! convention). Two rasters are produced: a **magnitude** raster (the Darcy
//! velocity `|v| = (T/n)·|∇h|`, the `output`) and a **direction** raster
//! (azimuth in degrees from north, clockwise, pointing down-gradient — the way
//! water flows). No-data / non-positive-porosity cells stay no-data.
//!
//! With a `seeds` point layer the tool also runs **advective particle tracking**:
//! from each seed it integrates the velocity field with a midpoint (RK2) scheme,
//! stepping a fixed `step` distance along the local flow direction for up to
//! `max_steps` steps until the path leaves the grid or stalls, emitting one
//! streamline polyline per seed (the `streamlines` output).

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::common::{load_input_raster, parse_optional_output, raster_like_with_data};
use crate::vector_common::{load_input_layer, write_or_store_layer};

/// Output no-data sentinel for the magnitude and direction rasters (magnitude is
/// non-negative and direction is in [0, 360), so a negative value is unambiguous).
const OUT_NODATA: f64 = -9999.0;

pub struct DarcyFlowTool;

impl Tool for DarcyFlowTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "darcy_flow",
            display_name: "Darcy Flow",
            summary: "Steady-state groundwater Darcy velocity field from head, transmissivity, and porosity rasters — a flow-magnitude (Darcy velocity) raster and a down-gradient flow-direction raster (v = -(T/n)·∇head), with optional advective particle tracking (RK2) from seed points into streamline polylines — like ArcGIS Darcy Flow / Particle Track.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Groundwater head (water-table elevation) raster. Cell size is taken from this raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "transmissivity",
                    description: "Transmissivity raster (co-registered to the head raster).",
                    required: true,
                },
                ToolParamSpec {
                    name: "porosity",
                    description: "Effective-porosity raster (co-registered to the head raster; values must be > 0).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output Darcy-velocity magnitude raster path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "direction",
                    description: "Optional output flow-direction raster path (azimuth degrees, from north clockwise, down-gradient). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read from each input raster (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "seeds",
                    description: "Optional point vector of particle seeds. When given, streamlines are traced from each seed.",
                    required: false,
                },
                ToolParamSpec {
                    name: "streamlines",
                    description: "Optional output line vector path for traced particle streamlines (requires 'seeds').",
                    required: false,
                },
                ToolParamSpec {
                    name: "step",
                    description: "Particle-tracking step length in CRS units (default: the head cell size).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_steps",
                    description: "Maximum integration steps per particle (default 1000).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "transmissivity")?;
        require_str(args, "porosity")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let head_path = require_str(args, "input")?;
        let trans_path = require_str(args, "transmissivity")?;
        let poro_path = require_str(args, "porosity")?;
        let out_mag = parse_optional_output(args, "output")?;
        let out_dir = parse_optional_output(args, "direction")?;
        let seeds_path = parse_optional_output(args, "seeds")?;
        let streamlines_path = parse_optional_output(args, "streamlines")?;
        let prm = parse_params(args)?;

        let head = load_input_raster(head_path)?;
        let trans = load_input_raster(trans_path)?;
        let poro = load_input_raster(poro_path)?;

        if prm.band < 0 || prm.band as usize >= head.bands {
            return Err(ToolError::Validation(format!(
                "band {} out of range (head raster has {} band(s))",
                prm.band + 1,
                head.bands
            )));
        }
        for (name, r) in [("transmissivity", &trans), ("porosity", &poro)] {
            if r.rows != head.rows || r.cols != head.cols {
                return Err(ToolError::Execution(format!(
                    "{name} raster ({}x{}) is not co-registered with the head raster ({}x{})",
                    r.rows, r.cols, head.rows, head.cols
                )));
            }
        }

        let rows = head.rows;
        let cols = head.cols;
        let cx = head.cell_size_x.abs().max(f64::MIN_POSITIVE);
        let cy = head.cell_size_y.abs().max(f64::MIN_POSITIVE);
        let hnd = head.nodata;
        let tnd = trans.nodata;
        let pnd = poro.nodata;

        // Read the three bands into flat grids (NaN marks invalid cells).
        let read = |r: &wbraster::Raster, nodata: f64| -> Vec<f64> {
            let mut g = vec![f64::NAN; rows * cols];
            for row in 0..rows {
                for col in 0..cols {
                    let v = r.get(prm.band, row as isize, col as isize);
                    if v != nodata && v.is_finite() {
                        g[row * cols + col] = v;
                    }
                }
            }
            g
        };
        let h = read(&head, hnd);
        let t = read(&trans, tnd);
        let p = read(&poro, pnd);

        ctx.progress
            .info("computing head gradient and Darcy velocity");

        // Horn 3×3 finite-difference gradient of the head surface. Missing
        // neighbours fall back to the centre value (exact for a planar surface
        // on the interior). Row increases south, so north = r-1.
        let at = |g: &[f64], r: isize, c: isize| -> Option<f64> {
            if r < 0 || c < 0 || r >= rows as isize || c >= cols as isize {
                return None;
            }
            let v = g[r as usize * cols + c as usize];
            if v.is_nan() {
                None
            } else {
                Some(v)
            }
        };

        let mut magnitude = vec![OUT_NODATA; rows * cols];
        let mut direction = vec![OUT_NODATA; rows * cols];
        // East / north velocity components (map units), NaN where invalid.
        let mut vel_e = vec![f64::NAN; rows * cols];
        let mut vel_n = vec![f64::NAN; rows * cols];

        let mut computed = 0usize;
        let (mut vmin, mut vmax, mut vsum) = (f64::INFINITY, f64::NEG_INFINITY, 0.0f64);

        for r in 0..rows as isize {
            for c in 0..cols as isize {
                let idx = r as usize * cols + c as usize;
                let head_c = match at(&h, r, c) {
                    Some(v) => v,
                    None => continue,
                };
                // Transmissivity / porosity must be valid and porosity > 0.
                let (Some(tc), Some(pc)) = (at(&t, r, c), at(&p, r, c)) else {
                    continue;
                };
                if pc <= 0.0 {
                    continue;
                }

                let g = |dr, dc| at(&h, r + dr, c + dc).unwrap_or(head_c);
                let (nw, n, ne) = (g(-1, -1), g(-1, 0), g(-1, 1));
                let (w, e) = (g(0, -1), g(0, 1));
                let (sw, s, se) = (g(1, -1), g(1, 0), g(1, 1));
                let dh_de = ((ne + 2.0 * e + se) - (nw + 2.0 * w + sw)) / (8.0 * cx);
                let dh_dn = ((nw + 2.0 * n + ne) - (sw + 2.0 * s + se)) / (8.0 * cy);

                // Darcy velocity v = -(T/n)·∇h.
                let k = tc / pc;
                let ve = -k * dh_de;
                let vn = -k * dh_dn;
                let mag = ve.hypot(vn);

                magnitude[idx] = mag;
                vel_e[idx] = ve;
                vel_n[idx] = vn;
                // Azimuth from north, clockwise, in the direction of flow.
                direction[idx] = if mag > 1e-12 {
                    let mut a = ve.atan2(vn).to_degrees();
                    if a < 0.0 {
                        a += 360.0;
                    }
                    a
                } else {
                    0.0
                };

                computed += 1;
                vmin = vmin.min(mag);
                vmax = vmax.max(mag);
                vsum += mag;
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        if computed == 0 {
            return Err(ToolError::Execution(
                "no valid cells (check no-data and that porosity > 0)".to_string(),
            ));
        }

        let mag_raster = raster_like_with_data(&head, magnitude, OUT_NODATA, DataType::F32)?;
        let dir_raster = raster_like_with_data(&head, direction, OUT_NODATA, DataType::F32)?;
        let mag_out = crate::common::write_or_store_output(mag_raster, out_mag)?;
        let dir_out = crate::common::write_or_store_output(dir_raster, out_dir)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(mag_out));
        outputs.insert("direction".to_string(), json!(dir_out));
        outputs.insert("cells_computed".to_string(), json!(computed));
        outputs.insert("min_velocity".to_string(), json!(vmin));
        outputs.insert("max_velocity".to_string(), json!(vmax));
        outputs.insert("mean_velocity".to_string(), json!(vsum / computed as f64));

        // ── Optional particle tracking ──────────────────────────────────────
        if let Some(seeds_path) = seeds_path {
            let seeds = load_input_layer(seeds_path)?;
            let field = VelocityField {
                ve: &vel_e,
                vn: &vel_n,
                rows,
                cols,
                x_min: head.x_min,
                y_max: head.y_min + rows as f64 * cy,
                cx,
                cy,
            };
            let step = prm.step.unwrap_or_else(|| cx.min(cy));

            let mut lines = Layer::new("streamlines").with_geom_type(GeometryType::LineString);
            if let Some(epsg) = head.crs.epsg {
                lines = lines.with_crs_epsg(epsg);
            }
            lines.add_field(FieldDef::new("seed_id", FieldType::Integer));
            lines.add_field(FieldDef::new("n_points", FieldType::Integer));
            lines.add_field(FieldDef::new("length", FieldType::Float));

            let mut streamline_count = 0usize;
            for (sid, feat) in seeds.iter().enumerate() {
                let Some((sx, sy)) = feat.geometry.as_ref().and_then(seed_xy) else {
                    continue;
                };
                let path = trace_particle(&field, sx, sy, step, prm.max_steps);
                if path.len() < 2 {
                    continue;
                }
                let length = path
                    .windows(2)
                    .map(|w| (w[1].0 - w[0].0).hypot(w[1].1 - w[0].1))
                    .sum::<f64>();
                let coords: Vec<Coord> = path.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
                lines.push(Feature {
                    fid: 0,
                    geometry: Some(Geometry::line_string(coords)),
                    attributes: vec![
                        FieldValue::Integer(sid as i64),
                        FieldValue::Integer(path.len() as i64),
                        FieldValue::Float(length),
                    ],
                });
                streamline_count += 1;
            }

            ctx.progress
                .info(&format!("{streamline_count} streamline(s) traced"));
            let lines_out = write_or_store_layer(lines, streamlines_path)?;
            outputs.insert("streamlines".to_string(), json!(lines_out));
            outputs.insert("streamline_count".to_string(), json!(streamline_count));
        }

        Ok(ToolRunResult { outputs })
    }
}

// ── Particle tracking ──────────────────────────────────────────────────────────

/// A continuous velocity field backed by two flat grids of east/north components
/// (NaN marks invalid cells). Sampling is bilinear over cell centres.
struct VelocityField<'a> {
    ve: &'a [f64],
    vn: &'a [f64],
    rows: usize,
    cols: usize,
    x_min: f64,
    y_max: f64,
    cx: f64,
    cy: f64,
}

impl VelocityField<'_> {
    /// Bilinearly samples the velocity at map coordinate `(x, y)`. Returns `None`
    /// outside the cell-centre grid or when any contributing cell is invalid.
    fn sample(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let col_f = (x - self.x_min) / self.cx - 0.5;
        let row_f = (self.y_max - y) / self.cy - 0.5;
        if col_f < 0.0
            || row_f < 0.0
            || col_f > (self.cols - 1) as f64
            || row_f > (self.rows - 1) as f64
        {
            return None;
        }
        let c0 = col_f.floor() as usize;
        let r0 = row_f.floor() as usize;
        let c1 = (c0 + 1).min(self.cols - 1);
        let r1 = (r0 + 1).min(self.rows - 1);
        let fx = col_f - c0 as f64;
        let fy = row_f - r0 as f64;

        let mut ve = 0.0;
        let mut vn = 0.0;
        for (r, c, w) in [
            (r0, c0, (1.0 - fx) * (1.0 - fy)),
            (r0, c1, fx * (1.0 - fy)),
            (r1, c0, (1.0 - fx) * fy),
            (r1, c1, fx * fy),
        ] {
            let e = self.ve[r * self.cols + c];
            let n = self.vn[r * self.cols + c];
            if e.is_nan() || n.is_nan() {
                return None;
            }
            ve += w * e;
            vn += w * n;
        }
        Some((ve, vn))
    }

    /// The unit flow direction at `(x, y)`, or `None` if unsampleable or stalled.
    fn unit_dir(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let (ve, vn) = self.sample(x, y)?;
        let mag = ve.hypot(vn);
        if mag < 1e-12 {
            None
        } else {
            Some((ve / mag, vn / mag))
        }
    }
}

/// Traces one particle from `(sx, sy)` by advancing a fixed `step` distance along
/// the local flow direction with a midpoint (RK2) scheme, up to `max_steps`.
fn trace_particle(
    field: &VelocityField,
    sx: f64,
    sy: f64,
    step: f64,
    max_steps: usize,
) -> Vec<(f64, f64)> {
    let mut path = vec![(sx, sy)];
    let (mut x, mut y) = (sx, sy);
    for _ in 0..max_steps {
        let Some((dx1, dy1)) = field.unit_dir(x, y) else {
            break;
        };
        // Midpoint direction.
        let (mx, my) = (x + 0.5 * step * dx1, y + 0.5 * step * dy1);
        let (dx, dy) = field.unit_dir(mx, my).unwrap_or((dx1, dy1));
        let (nx, ny) = (x + step * dx, y + step * dy);
        path.push((nx, ny));
        // Stop if the point left the sampleable field.
        if field.sample(nx, ny).is_none() {
            break;
        }
        x = nx;
        y = ny;
    }
    path
}

fn seed_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        _ => None,
    }
}

// ── Parameters ──────────────────────────────────────────────────────────────────

struct Params {
    band: isize,
    step: Option<f64>,
    max_steps: usize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
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
        band: (band_1based - 1) as isize,
        step,
        max_steps,
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

    /// Builds a single-band raster (cell size 1, EPSG:3857) from row-major data.
    fn raster_from(cols: usize, rows: usize, data: Vec<f64>) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
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

    fn constant(cols: usize, rows: usize, v: f64) -> String {
        raster_from(cols, rows, vec![v; cols * rows])
    }

    fn get(r: &Raster, row: usize, col: usize) -> f64 {
        r.get(0, row as isize, col as isize)
    }

    /// On a planar head surface with a known linear gradient, the Darcy velocity
    /// magnitude equals the analytic (T/n)·slope and the direction points
    /// down-gradient — the core correctness property.
    #[test]
    fn planar_head_matches_analytic_velocity() {
        let (cols, rows) = (6, 5);
        // Head increases eastward: h = 3·col (slope 3 per cell, cell size 1).
        let slope = 3.0;
        let mut head = vec![0.0; cols * rows];
        for row in 0..rows {
            for col in 0..cols {
                head[row * cols + col] = slope * col as f64;
            }
        }
        let head = raster_from(cols, rows, head);
        let trans = constant(cols, rows, 2.0); // T = 2
        let poro = constant(cols, rows, 0.5); // n = 0.5  ->  T/n = 4

        let args: ToolArgs = serde_json::from_value(json!({
            "input": head, "transmissivity": trans, "porosity": poro,
        }))
        .unwrap();
        let out = DarcyFlowTool.run(&args, &ctx()).unwrap();
        let mag = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        let dir = load_input_raster(out.outputs["direction"].as_str().unwrap()).unwrap();

        let expected = 4.0 * slope; // (T/n)·slope = 12
                                    // Check an interior cell (Horn is exact there for a plane).
        assert!(
            (get(&mag, 2, 3) - expected).abs() < 1e-6,
            "magnitude {} != analytic {expected}",
            get(&mag, 2, 3)
        );
        // Head rises eastward -> water flows west -> azimuth 270°.
        assert!(
            (get(&dir, 2, 3) - 270.0).abs() < 1e-6,
            "direction {} != 270",
            get(&dir, 2, 3)
        );
    }

    /// A flat head surface produces zero Darcy velocity everywhere.
    #[test]
    fn flat_head_gives_zero_velocity() {
        let (cols, rows) = (5, 5);
        let head = constant(cols, rows, 42.0);
        let trans = constant(cols, rows, 3.0);
        let poro = constant(cols, rows, 0.3);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": head, "transmissivity": trans, "porosity": poro,
        }))
        .unwrap();
        let out = DarcyFlowTool.run(&args, &ctx()).unwrap();
        let mag = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        for row in 0..rows {
            for col in 0..cols {
                assert!(get(&mag, row, col).abs() < 1e-9);
            }
        }
        assert_eq!(out.outputs["max_velocity"], json!(0.0));
    }

    /// Particle tracking traces a streamline that flows down-gradient (westward,
    /// since head rises eastward).
    #[test]
    fn particle_tracking_flows_down_gradient() {
        let (cols, rows) = (10, 3);
        let mut head = vec![0.0; cols * rows];
        for row in 0..rows {
            for col in 0..cols {
                head[row * cols + col] = 2.0 * col as f64; // rises eastward
            }
        }
        let head = raster_from(cols, rows, head);
        let trans = constant(cols, rows, 1.0);
        let poro = constant(cols, rows, 0.25);

        // Seed near the east edge; the particle should march west.
        let mut seeds = Layer::new("seeds")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        seeds
            .add_feature(Some(Geometry::point(8.5, 1.5)), &[])
            .unwrap();
        let seeds_path = {
            let id = wbvector::memory_store::put_vector(seeds);
            wbvector::memory_store::make_vector_memory_path(&id)
        };

        let args: ToolArgs = serde_json::from_value(json!({
            "input": head, "transmissivity": trans, "porosity": poro,
            "seeds": seeds_path, "step": 0.5, "max_steps": 50,
        }))
        .unwrap();
        let out = DarcyFlowTool.run(&args, &ctx()).unwrap();
        assert_eq!(out.outputs["streamline_count"], json!(1));
        let lines = load_input_layer(out.outputs["streamlines"].as_str().unwrap()).unwrap();
        let geom = lines.features[0].geometry.as_ref().unwrap();
        let Geometry::LineString(cs) = geom else {
            panic!("expected a line string");
        };
        assert!(cs.len() >= 2);
        // The streamline advances westward (x decreases).
        assert!(
            cs.last().unwrap().x < cs.first().unwrap().x - 1.0,
            "streamline should flow west (down-gradient)"
        );
    }

    /// Rasters that are not co-registered are rejected.
    #[test]
    fn rejects_mismatched_dimensions() {
        let head = constant(5, 5, 1.0);
        let trans = constant(4, 5, 1.0); // wrong width
        let poro = constant(5, 5, 0.3);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": head, "transmissivity": trans, "porosity": poro,
        }))
        .unwrap();
        assert!(DarcyFlowTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            DarcyFlowTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "h.tif" })).is_err()); // missing T, n
        assert!(bad(json!({ "input": "h.tif", "transmissivity": "t.tif" })).is_err());
        assert!(bad(
            json!({ "input": "h.tif", "transmissivity": "t.tif", "porosity": "p.tif", "step": -1 })
        )
        .is_err());
        assert!(
            bad(json!({ "input": "h.tif", "transmissivity": "t.tif", "porosity": "p.tif" }))
                .is_ok()
        );
    }
}
