//! GeoLibre tool: volumetric integration of a raster surface relative to a
//! horizontal reference plane — a pure-Rust counterpart of ArcGIS 3D Analyst's
//! *Surface Volume* tool.
//!
//! For a single elevation surface and a reference-plane elevation it reports,
//! for the portion of the surface ABOVE and/or BELOW that plane:
//!   * **2D area** — the planimetric footprint of the qualifying cells
//!     (`count * cellArea`);
//!   * **3D surface area** — the draped area of those cells, each cell's footprint
//!     scaled by `sqrt(1 + (dz/dx)^2 + (dz/dy)^2)` (the tilt of the local surface);
//!   * **volume** — the prism between the surface and the plane,
//!     `sum(|z - plane| * cellArea)` over the qualifying cells.
//!
//! This is the earthworks cut/fill and reservoir-capacity primitive that neither
//! `cut_fill` (which differences two surfaces) nor `add_surface_information`
//! (per-feature Z statistics) provides: integration of one surface against a
//! horizontal datum.
//!
//! A cell qualifies for ABOVE when `z >= plane` and for BELOW when `z <= plane`,
//! so a plane placed at the surface minimum yields ABOVE coverage over every
//! valid cell (volume `sum((z - min) * cellArea)`), and a plane at the maximum
//! yields BELOW coverage over every valid cell. No-data cells are ignored.
//!
//! The gradient used for the 3D area is a central difference over the immediate
//! orthogonal neighbours (one-sided at the raster edge, zero where a neighbour is
//! no-data), matching the standard finite-difference slope. The result is a small
//! results table (one row per requested direction), returned in the run result
//! and, when `output` is given, written as CSV.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;

use crate::common::{load_input_raster, write_text_output};
use crate::vector_common::parse_optional_str;

/// Which side(s) of the reference plane to integrate.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Direction {
    Above,
    Below,
    Both,
}

impl Direction {
    fn parse(s: &str) -> Result<Direction, ToolError> {
        match s.trim().to_ascii_lowercase().as_str() {
            "above" => Ok(Direction::Above),
            "below" => Ok(Direction::Below),
            "both" => Ok(Direction::Both),
            other => Err(ToolError::Validation(format!(
                "parameter 'direction' must be 'above', 'below', or 'both' (got '{other}')"
            ))),
        }
    }

    /// The individual directions this selection expands to, in table order.
    fn parts(self) -> &'static [Direction] {
        match self {
            Direction::Above => &[Direction::Above],
            Direction::Below => &[Direction::Below],
            Direction::Both => &[Direction::Above, Direction::Below],
        }
    }

    fn label(self) -> &'static str {
        match self {
            Direction::Above => "Above",
            Direction::Below => "Below",
            Direction::Both => "Both",
        }
    }
}

/// Integrated area/volume totals for one direction.
struct Totals {
    n_cells: u64,
    area_2d: f64,
    area_3d: f64,
    volume: f64,
}

pub struct SurfaceVolumeTool;

impl Tool for SurfaceVolumeTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "surface_volume",
            display_name: "Surface Volume",
            summary: "Integrate a raster surface against a horizontal reference plane, reporting 2D area, 3D surface area, and volume for the portion above and/or below the plane (like ArcGIS's Surface Volume). Earthworks cut/fill and reservoir-capacity analysis.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input elevation surface raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "reference_plane",
                    description: "Reference-plane elevation to integrate against (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "direction",
                    description: "Which side of the plane to integrate: 'above', 'below', or 'both' (default 'above').",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output CSV path for the results table. If omitted, results are returned in the run result only.",
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

        let raster = load_input_raster(input)?;
        let band_1based = prm.band;
        if (band_1based as usize) > raster.bands {
            return Err(ToolError::Validation(format!(
                "band {band_1based} out of range (raster has {} band(s))",
                raster.bands
            )));
        }
        let band = (band_1based - 1) as isize;

        let cell_area = raster.cell_size_x.abs() * raster.cell_size_y.abs();
        if !cell_area.is_finite() || cell_area <= 0.0 {
            return Err(ToolError::Execution(
                "raster has a non-positive cell size; cannot compute areas".to_string(),
            ));
        }

        ctx.progress
            .info(&format!("integrating surface at plane z = {}", prm.plane));

        let parts = prm.direction.parts();
        let mut above = Totals::zero();
        let mut below = Totals::zero();
        let want_above = parts.contains(&Direction::Above);
        let want_below = parts.contains(&Direction::Below);

        let rows = raster.rows as isize;
        let cols = raster.cols as isize;
        let nodata = raster.nodata;
        let dx = raster.cell_size_x.abs();
        let dy = raster.cell_size_y.abs();

        for row in 0..rows {
            for col in 0..cols {
                let z = raster.get(band, row, col);
                if z == nodata || !z.is_finite() {
                    continue;
                }
                // 3D-area scale from the local finite-difference gradient.
                let factor = surface_factor(&raster, band, row, col, dx, dy, nodata);
                let cell_3d = cell_area * factor;

                if want_above && z >= prm.plane {
                    above.n_cells += 1;
                    above.area_2d += cell_area;
                    above.area_3d += cell_3d;
                    above.volume += (z - prm.plane) * cell_area;
                }
                if want_below && z <= prm.plane {
                    below.n_cells += 1;
                    below.area_2d += cell_area;
                    below.area_3d += cell_3d;
                    below.volume += (prm.plane - z) * cell_area;
                }
            }
            ctx.progress.progress((row as f64 + 1.0) / rows as f64);
        }

        // Build the results table (one row per requested direction).
        let mut csv = String::from("direction,area_2d,area_3d,volume,n_cells\n");
        let mut outputs = BTreeMap::new();
        for &part in parts {
            let t = match part {
                Direction::Above => &above,
                Direction::Below => &below,
                Direction::Both => unreachable!("parts() never yields Both"),
            };
            csv.push_str(&format!(
                "{},{},{},{},{}\n",
                part.label(),
                t.area_2d,
                t.area_3d,
                t.volume,
                t.n_cells
            ));
            let key = part.label().to_ascii_lowercase();
            outputs.insert(format!("{key}_area_2d"), json!(t.area_2d));
            outputs.insert(format!("{key}_area_3d"), json!(t.area_3d));
            outputs.insert(format!("{key}_volume"), json!(t.volume));
            outputs.insert(format!("{key}_n_cells"), json!(t.n_cells));
        }

        if let Some(path) = output {
            write_text_output(&csv, path)?;
            outputs.insert("output".to_string(), json!(path));
        }
        outputs.insert("reference_plane".to_string(), json!(prm.plane));
        outputs.insert("direction".to_string(), json!(prm.direction.label()));
        outputs.insert("cell_area".to_string(), json!(cell_area));

        Ok(ToolRunResult { outputs })
    }
}

impl Totals {
    fn zero() -> Totals {
        Totals {
            n_cells: 0,
            area_2d: 0.0,
            area_3d: 0.0,
            volume: 0.0,
        }
    }
}

/// Per-cell 3D-area scale factor `sqrt(1 + (dz/dx)^2 + (dz/dy)^2)` from a central
/// finite-difference gradient. Falls back to a one-sided difference at the raster
/// edge and treats no-data neighbours as absent (0 slope on that axis).
fn surface_factor(
    raster: &Raster,
    band: isize,
    row: isize,
    col: isize,
    dx: f64,
    dy: f64,
    nodata: f64,
) -> f64 {
    let center = raster.get(band, row, col);
    let left = sample(raster, band, row, col - 1, nodata);
    let right = sample(raster, band, row, col + 1, nodata);
    let up = sample(raster, band, row - 1, col, nodata);
    let down = sample(raster, band, row + 1, col, nodata);

    let gx = one_dim_gradient(left, right, center, dx).unwrap_or(0.0);
    let gy = one_dim_gradient(up, down, center, dy).unwrap_or(0.0);
    (1.0 + gx * gx + gy * gy).sqrt()
}

/// Reads a neighbour cell, returning `None` for out-of-bounds or no-data.
fn sample(raster: &Raster, band: isize, row: isize, col: isize, nodata: f64) -> Option<f64> {
    if row < 0 || col < 0 || row >= raster.rows as isize || col >= raster.cols as isize {
        return None;
    }
    let v = raster.get(band, row, col);
    if v == nodata || !v.is_finite() {
        None
    } else {
        Some(v)
    }
}

/// Gradient across one axis from the two opposite neighbours `a` (lower index)
/// and `b` (higher index) at spacing `step`, falling back to a one-sided slope
/// against `center` when only one neighbour is valid. Returns `None` when the
/// slope cannot be estimated (both neighbours absent).
fn one_dim_gradient(a: Option<f64>, b: Option<f64>, center: f64, step: f64) -> Option<f64> {
    match (a, b) {
        (Some(a), Some(b)) => Some((b - a) / (2.0 * step)),
        (Some(a), None) => Some((center - a) / step),
        (None, Some(b)) => Some((b - center) / step),
        (None, None) => None,
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    plane: f64,
    direction: Direction,
    band: u64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let plane = parse_optional_f64(args, "reference_plane")?.unwrap_or(0.0);
    if !plane.is_finite() {
        return Err(ToolError::Validation(
            "parameter 'reference_plane' must be finite".to_string(),
        ));
    }

    let direction = match parse_optional_str(args, "direction")? {
        None => Direction::Above,
        Some(s) => Direction::parse(s)?,
    };

    let band = match parse_optional_f64(args, "band")? {
        None => 1,
        Some(v) if v.fract() == 0.0 && v >= 1.0 => v as u64,
        Some(_) => {
            return Err(ToolError::Validation(
                "parameter 'band' must be a positive integer".to_string(),
            ))
        }
    };

    Ok(Params {
        plane,
        direction,
        band,
    })
}

/// Parses an optional numeric parameter, accepting a JSON number or a numeric
/// string (host UIs post scalars as strings).
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

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{memory_store, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a single-band raster from row-major elevations with unit cells.
    fn raster_path(rows: usize, cols: usize, cell: f64, vals: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: cell,
            cell_size_y: Some(cell),
            nodata: -9999.0,
            data_type: wbraster::DataType::F32,
            crs: Default::default(),
            metadata: Default::default(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, vals[row * cols + col])
                    .unwrap();
            }
        }
        let id = memory_store::put_raster(r);
        memory_store::make_raster_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        SurfaceVolumeTool.run(&args, &ctx()).unwrap()
    }

    #[test]
    fn flat_surface_volume_and_area_are_exact() {
        // 4x4 flat surface at z = 5, unit cells: 16 cells, area 16.
        let input = raster_path(4, 4, 1.0, &[5.0; 16]);
        let out = run(json!({ "input": input, "reference_plane": 0.0, "direction": "above" }));
        assert_eq!(out.outputs["above_n_cells"], json!(16));
        assert!((out.outputs["above_area_2d"].as_f64().unwrap() - 16.0).abs() < 1e-9);
        // Flat surface: 3D area == 2D area.
        assert!((out.outputs["above_area_3d"].as_f64().unwrap() - 16.0).abs() < 1e-9);
        // Volume = height 5 * area 16 = 80.
        assert!((out.outputs["above_volume"].as_f64().unwrap() - 80.0).abs() < 1e-9);
    }

    #[test]
    fn nonunit_cell_size_scales_area_and_volume() {
        // 2x2 flat at z = 3, cell size 10 -> cell area 100, total area 400.
        let input = raster_path(2, 2, 10.0, &[3.0; 4]);
        let out = run(json!({ "input": input, "reference_plane": 1.0, "direction": "above" }));
        assert!((out.outputs["above_area_2d"].as_f64().unwrap() - 400.0).abs() < 1e-6);
        // Volume = (3 - 1) * 400 = 800.
        assert!((out.outputs["above_volume"].as_f64().unwrap() - 800.0).abs() < 1e-6);
        assert!((out.outputs["cell_area"].as_f64().unwrap() - 100.0).abs() < 1e-9);
    }

    #[test]
    fn plane_at_min_covers_all_cells_above() {
        // Ramp 1..=9; plane at the min (1) -> every cell qualifies for above.
        let vals: Vec<f64> = (1..=9).map(|i| i as f64).collect();
        let input = raster_path(3, 3, 1.0, &vals);
        let out = run(json!({ "input": input, "reference_plane": 1.0, "direction": "above" }));
        assert_eq!(out.outputs["above_n_cells"], json!(9));
        // Volume = sum(z - 1) = sum(0,1,2,...,8) = 36.
        assert!((out.outputs["above_volume"].as_f64().unwrap() - 36.0).abs() < 1e-9);
    }

    #[test]
    fn plane_at_max_covers_all_cells_below() {
        let vals: Vec<f64> = (1..=9).map(|i| i as f64).collect();
        let input = raster_path(3, 3, 1.0, &vals);
        let out = run(json!({ "input": input, "reference_plane": 9.0, "direction": "below" }));
        assert_eq!(out.outputs["below_n_cells"], json!(9));
        // Volume = sum(9 - z) = sum(8,7,...,0) = 36.
        assert!((out.outputs["below_volume"].as_f64().unwrap() - 36.0).abs() < 1e-9);
    }

    #[test]
    fn both_directions_partition_around_the_plane() {
        // Values 0,10 twice; plane 5: above gets the 10s, below gets the 0s.
        let input = raster_path(2, 2, 1.0, &[0.0, 10.0, 10.0, 0.0]);
        let out = run(json!({ "input": input, "reference_plane": 5.0, "direction": "both" }));
        assert_eq!(out.outputs["above_n_cells"], json!(2));
        assert_eq!(out.outputs["below_n_cells"], json!(2));
        // Above volume = 2 * (10 - 5) = 10; below = 2 * (5 - 0) = 10.
        assert!((out.outputs["above_volume"].as_f64().unwrap() - 10.0).abs() < 1e-9);
        assert!((out.outputs["below_volume"].as_f64().unwrap() - 10.0).abs() < 1e-9);
    }

    #[test]
    fn nodata_cells_are_ignored() {
        // One no-data cell; only 3 valid cells contribute.
        let input = raster_path(2, 2, 1.0, &[4.0, 4.0, 4.0, -9999.0]);
        let out = run(json!({ "input": input, "reference_plane": 0.0, "direction": "above" }));
        assert_eq!(out.outputs["above_n_cells"], json!(3));
        assert!((out.outputs["above_volume"].as_f64().unwrap() - 12.0).abs() < 1e-9);
    }

    #[test]
    fn tilted_surface_3d_area_exceeds_2d_area() {
        // A constant east-west slope of 1 per cell: 3D area = 2D * sqrt(2).
        let vals: Vec<f64> = (0..16).map(|i| (i % 4) as f64).collect();
        let input = raster_path(4, 4, 1.0, &vals);
        let out = run(json!({ "input": input, "reference_plane": -1.0, "direction": "above" }));
        let a2 = out.outputs["above_area_2d"].as_f64().unwrap();
        let a3 = out.outputs["above_area_3d"].as_f64().unwrap();
        assert!(a3 > a2, "3D area {a3} should exceed 2D area {a2}");
        // Interior cells have gradient 1 -> factor sqrt(2); overall > 2D.
        assert!(
            a3 / a2 > 1.2,
            "expected a meaningful tilt, got ratio {}",
            a3 / a2
        );
    }

    #[test]
    fn csv_output_is_written() {
        let input = raster_path(2, 2, 1.0, &[1.0, 2.0, 3.0, 4.0]);
        let dir = std::env::temp_dir();
        let path = dir.join("surface_volume_test.csv");
        let path_str = path.to_str().unwrap();
        let out = run(json!({
            "input": input, "reference_plane": 0.0, "direction": "both", "output": path_str
        }));
        assert_eq!(out.outputs["output"], json!(path_str));
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.starts_with("direction,area_2d,area_3d,volume,n_cells\n"));
        assert!(text.contains("Above,"));
        assert!(text.contains("Below,"));
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = SurfaceVolumeTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input");
        assert!(
            bad(json!({ "input": "x.tif", "direction": "sideways" })).is_err(),
            "bad direction"
        );
        assert!(
            bad(json!({ "input": "x.tif", "reference_plane": "abc" })).is_err(),
            "non-numeric plane"
        );
        assert!(
            bad(json!({ "input": "x.tif", "band": 0 })).is_err(),
            "band must be >= 1"
        );
        assert!(
            bad(json!({ "input": "x.tif", "reference_plane": 100, "direction": "both" })).is_ok()
        );
    }
}
