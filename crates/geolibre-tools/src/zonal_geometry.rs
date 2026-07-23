//! GeoLibre tool: per-zone geometric measures from a categorical zone raster.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Zonal Geometry* / *Zonal Geometry As
//! Table* (Spatial Analyst). The bundled whitebox suite ships `zonal_statistics`
//! (per-zone summaries of a **value** raster) and the repo added `zonal_histogram`
//! (#248, per-zone value **distribution**); neither characterizes the **shape** of
//! each zone. This is the missing third leg: it measures the geometry of every
//! zone — area, perimeter, thickness, centroid, and standard-deviational ellipse
//! axes + orientation.
//!
//! Following ArcGIS, every distinct integer value in the zone raster is one zone
//! (zones need not be contiguous). For each zone the tool accumulates:
//! * `area` — cell count × cell area;
//! * `perimeter` — total length of zone/background boundary edges;
//! * `thickness` — 2 × the largest inscribed radius (max distance-to-boundary), approximated with a chamfer distance transform;
//! * `centroid_x` / `centroid_y` — mean of the zone's cell-center coordinates;
//! * `major_axis` / `minor_axis` — 4σ diameters of the best-fit (standard-deviational) ellipse from the coordinate covariance;
//! * `orientation` — ellipse major-axis azimuth in degrees, arithmetic ([0,180)).
//!
//! Raster mode writes the chosen `measure` back to every cell of its zone; table
//! mode (`as_table=true`) returns one row per zone with all measures (and writes a
//! CSV when `output` is a path).

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_text_output,
};

pub struct ZonalGeometryTool;

impl Tool for ZonalGeometryTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "zonal_geometry",
            display_name: "Zonal Geometry",
            summary: "Per-zone geometric measures from a categorical zone raster (like ArcGIS Zonal Geometry): area, perimeter, thickness, centroid, and standard-deviational ellipse major/minor axis + orientation — the zone-shape counterpart the bundled value-based zonal_statistics and zonal_histogram don't provide.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "zones",
                    description: "Categorical zone raster (each distinct integer value is one zone).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output raster (raster mode) or CSV path (as_table mode). If omitted, stored in memory / returned in the result only.",
                    required: false,
                },
                ToolParamSpec {
                    name: "measure",
                    description: "Raster-mode measure: 'area' (default), 'perimeter', 'thickness', 'centroid_x', 'centroid_y', 'major_axis', 'minor_axis', or 'orientation'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "as_table",
                    description: "When true, emit one row per zone with every measure instead of a per-cell raster (default false).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "zones")?;
        parse_measure(args)?;
        parse_bool(args, "as_table")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let zones_path = require_str(args, "zones")?;
        let output = parse_optional_output(args, "output")?;
        let measure = parse_measure(args)?;
        let as_table = parse_bool(args, "as_table")?.unwrap_or(false);

        let zr = load_input_raster(zones_path)?;
        let (rows, cols) = (zr.rows, zr.cols);
        let cell_x = zr.cell_size_x.abs();
        let cell_y = zr.cell_size_y.abs();
        let cell_area = cell_x * cell_y;
        // Chamfer transform assumes a roughly square cell for the diagonal step.
        let cell_avg = 0.5 * (cell_x + cell_y);

        // Read zones into a dense i64 label buffer (nodata -> None).
        let mut zone: Vec<Option<i64>> = vec![None; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let v = zr.get(0, r as isize, c as isize);
                if v != zr.nodata && v.is_finite() {
                    zone[r * cols + c] = Some(v.round() as i64);
                }
            }
        }

        // Per-zone accumulators.
        struct Acc {
            n: u64,
            perim_edges_x: u64, // vertical boundary edges (contribute cell_y each)
            perim_edges_y: u64, // horizontal boundary edges (contribute cell_x each)
            sx: f64,
            sy: f64,
            sxx: f64,
            syy: f64,
            sxy: f64,
            max_dt: f64,
        }
        let mut acc: BTreeMap<i64, Acc> = BTreeMap::new();

        // Distance-to-boundary transform (chamfer) shared across all zones: seed
        // every boundary cell at half a cell, propagate inward. A cell's value is
        // its distance to the nearest zone boundary, so per-zone max is the
        // largest inscribed radius.
        let big = f64::INFINITY;
        let mut dt = vec![big; rows * cols];
        let neighbour = |r: usize, c: usize, dr: isize, dc: isize| -> Option<Option<i64>> {
            let nr = r as isize + dr;
            let nc = c as isize + dc;
            if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                return None; // off-grid
            }
            Some(zone[nr as usize * cols + nc as usize])
        };
        for r in 0..rows {
            for c in 0..cols {
                let Some(z) = zone[r * cols + c] else {
                    continue;
                };
                // Boundary if any 4-neighbour is a different zone / nodata / edge.
                let mut boundary = false;
                for (dr, dc) in [(-1, 0), (1, 0), (0, -1), (0, 1)] {
                    match neighbour(r, c, dr, dc) {
                        None => boundary = true,
                        Some(nz) if nz != Some(z) => boundary = true,
                        _ => {}
                    }
                }
                if boundary {
                    dt[r * cols + c] = 0.5 * cell_avg;
                }
            }
        }
        chamfer(&mut dt, &zone, rows, cols, cell_avg);

        for r in 0..rows {
            for c in 0..cols {
                let Some(z) = zone[r * cols + c] else {
                    continue;
                };
                let e = acc.entry(z).or_insert(Acc {
                    n: 0,
                    perim_edges_x: 0,
                    perim_edges_y: 0,
                    sx: 0.0,
                    sy: 0.0,
                    sxx: 0.0,
                    syy: 0.0,
                    sxy: 0.0,
                    max_dt: 0.0,
                });
                e.n += 1;
                // Cell-center coordinates (raster y increases downward -> flip).
                let cx = zr.x_min + (c as f64 + 0.5) * cell_x;
                let cy = zr.y_min + (rows as f64 - 1.0 - r as f64 + 0.5) * cell_y;
                e.sx += cx;
                e.sy += cy;
                e.sxx += cx * cx;
                e.syy += cy * cy;
                e.sxy += cx * cy;
                let d = dt[r * cols + c];
                if d.is_finite() && d > e.max_dt {
                    e.max_dt = d;
                }
                // Boundary edges: count faces touching a different zone / edge.
                if neighbour(r, c, 0, -1)
                    .map(|nz| nz != Some(z))
                    .unwrap_or(true)
                {
                    e.perim_edges_x += 1;
                }
                if neighbour(r, c, 0, 1)
                    .map(|nz| nz != Some(z))
                    .unwrap_or(true)
                {
                    e.perim_edges_x += 1;
                }
                if neighbour(r, c, -1, 0)
                    .map(|nz| nz != Some(z))
                    .unwrap_or(true)
                {
                    e.perim_edges_y += 1;
                }
                if neighbour(r, c, 1, 0)
                    .map(|nz| nz != Some(z))
                    .unwrap_or(true)
                {
                    e.perim_edges_y += 1;
                }
            }
        }

        if acc.is_empty() {
            return Err(ToolError::Execution(
                "zone raster has no valid cells".to_string(),
            ));
        }

        ctx.progress
            .info(&format!("{} zone(s) over {rows}x{cols}", acc.len()));

        // Finalize measures per zone.
        struct ZoneStats {
            zone: i64,
            area: f64,
            perimeter: f64,
            thickness: f64,
            centroid_x: f64,
            centroid_y: f64,
            major_axis: f64,
            minor_axis: f64,
            orientation: f64,
        }
        let mut stats: Vec<ZoneStats> = Vec::with_capacity(acc.len());
        for (z, a) in &acc {
            let n = a.n as f64;
            let area = n * cell_area;
            let perimeter = a.perim_edges_x as f64 * cell_y + a.perim_edges_y as f64 * cell_x;
            let thickness = 2.0 * a.max_dt;
            let mx = a.sx / n;
            let my = a.sy / n;
            let cxx = (a.sxx / n - mx * mx).max(0.0);
            let cyy = (a.syy / n - my * my).max(0.0);
            let cxy = a.sxy / n - mx * my;
            let half = 0.5 * (cxx + cyy);
            let disc = (0.5 * (cxx - cyy)).hypot(cxy);
            let l1 = (half + disc).max(0.0);
            let l2 = (half - disc).max(0.0);
            // 4σ diameters of the best-fit ellipse.
            let major_axis = 4.0 * l1.sqrt();
            let minor_axis = 4.0 * l2.sqrt();
            // Major-axis azimuth (arithmetic, degrees in [0,180)).
            let mut orientation = if cxy == 0.0 && cxx >= cyy {
                0.0
            } else if cxy == 0.0 {
                90.0
            } else {
                0.5 * (2.0 * cxy).atan2(cxx - cyy).to_degrees()
            };
            while orientation < 0.0 {
                orientation += 180.0;
            }
            while orientation >= 180.0 {
                orientation -= 180.0;
            }
            stats.push(ZoneStats {
                zone: *z,
                area,
                perimeter,
                thickness,
                centroid_x: mx,
                centroid_y: my,
                major_axis,
                minor_axis,
                orientation,
            });
        }

        let mut outputs = BTreeMap::new();
        outputs.insert("zone_count".to_string(), json!(stats.len()));

        if as_table {
            let mut csv = String::from(
                "zone,area,perimeter,thickness,centroid_x,centroid_y,major_axis,minor_axis,orientation\n",
            );
            let mut zones_json = Vec::new();
            for s in &stats {
                csv.push_str(&format!(
                    "{},{},{},{},{},{},{},{},{}\n",
                    s.zone,
                    s.area,
                    s.perimeter,
                    s.thickness,
                    s.centroid_x,
                    s.centroid_y,
                    s.major_axis,
                    s.minor_axis,
                    s.orientation
                ));
                zones_json.push(json!({
                    "zone": s.zone,
                    "area": s.area,
                    "perimeter": s.perimeter,
                    "thickness": s.thickness,
                    "centroid_x": s.centroid_x,
                    "centroid_y": s.centroid_y,
                    "major_axis": s.major_axis,
                    "minor_axis": s.minor_axis,
                    "orientation": s.orientation,
                }));
            }
            if let Some(path) = output {
                write_text_output(&csv, path)?;
                outputs.insert("output".to_string(), json!(path));
            }
            outputs.insert("zones".to_string(), json!(zones_json));
            outputs.insert("measure".to_string(), json!("table"));
        } else {
            let lut: BTreeMap<i64, f64> = stats
                .iter()
                .map(|s| {
                    (
                        s.zone,
                        measure.value(
                            s.area,
                            s.perimeter,
                            s.thickness,
                            s.centroid_x,
                            s.centroid_y,
                            s.major_axis,
                            s.minor_axis,
                            s.orientation,
                        ),
                    )
                })
                .collect();
            let nodata = -9999.0_f64;
            let mut data = vec![nodata; rows * cols];
            for (i, z) in zone.iter().enumerate() {
                if let Some(z) = z {
                    if let Some(v) = lut.get(z) {
                        data[i] = *v;
                    }
                }
            }
            let out_r = raster_like_with_data(&zr, data, nodata, DataType::F32)?;
            let out_path = crate::common::write_or_store_output(out_r, output)?;
            outputs.insert("output".to_string(), json!(out_path));
            outputs.insert("measure".to_string(), json!(measure.label()));
        }
        Ok(ToolRunResult { outputs })
    }
}

/// Two-pass chamfer distance transform over valid (in-zone) cells. Seeds are the
/// cells already set finite; nodata cells are barriers left at +inf.
fn chamfer(dt: &mut [f64], zone: &[Option<i64>], rows: usize, cols: usize, cell: f64) {
    let ortho = cell;
    let diag = cell * std::f64::consts::SQRT_2;
    let idx = |r: usize, c: usize| r * cols + c;
    // Forward pass (top-left to bottom-right).
    for r in 0..rows {
        for c in 0..cols {
            if zone[idx(r, c)].is_none() {
                continue;
            }
            let mut best = dt[idx(r, c)];
            let mut relax = |nr: isize, nc: isize, w: f64| {
                if nr >= 0 && nc >= 0 && (nr as usize) < rows && (nc as usize) < cols {
                    let ni = idx(nr as usize, nc as usize);
                    if zone[ni].is_some() {
                        best = best.min(dt[ni] + w);
                    }
                }
            };
            let (ri, ci) = (r as isize, c as isize);
            relax(ri - 1, ci, ortho);
            relax(ri, ci - 1, ortho);
            relax(ri - 1, ci - 1, diag);
            relax(ri - 1, ci + 1, diag);
            dt[idx(r, c)] = best;
        }
    }
    // Backward pass (bottom-right to top-left).
    for r in (0..rows).rev() {
        for c in (0..cols).rev() {
            if zone[idx(r, c)].is_none() {
                continue;
            }
            let mut best = dt[idx(r, c)];
            let mut relax = |nr: isize, nc: isize, w: f64| {
                if nr >= 0 && nc >= 0 && (nr as usize) < rows && (nc as usize) < cols {
                    let ni = idx(nr as usize, nc as usize);
                    if zone[ni].is_some() {
                        best = best.min(dt[ni] + w);
                    }
                }
            };
            let (ri, ci) = (r as isize, c as isize);
            relax(ri + 1, ci, ortho);
            relax(ri, ci + 1, ortho);
            relax(ri + 1, ci + 1, diag);
            relax(ri + 1, ci - 1, diag);
            dt[idx(r, c)] = best;
        }
    }
}

// ── Measure selection ─────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Measure {
    Area,
    Perimeter,
    Thickness,
    CentroidX,
    CentroidY,
    MajorAxis,
    MinorAxis,
    Orientation,
}

impl Measure {
    fn label(&self) -> &'static str {
        match self {
            Measure::Area => "area",
            Measure::Perimeter => "perimeter",
            Measure::Thickness => "thickness",
            Measure::CentroidX => "centroid_x",
            Measure::CentroidY => "centroid_y",
            Measure::MajorAxis => "major_axis",
            Measure::MinorAxis => "minor_axis",
            Measure::Orientation => "orientation",
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn value(
        &self,
        area: f64,
        perimeter: f64,
        thickness: f64,
        cx: f64,
        cy: f64,
        major: f64,
        minor: f64,
        orientation: f64,
    ) -> f64 {
        match self {
            Measure::Area => area,
            Measure::Perimeter => perimeter,
            Measure::Thickness => thickness,
            Measure::CentroidX => cx,
            Measure::CentroidY => cy,
            Measure::MajorAxis => major,
            Measure::MinorAxis => minor,
            Measure::Orientation => orientation,
        }
    }
}

fn parse_measure(args: &ToolArgs) -> Result<Measure, ToolError> {
    Ok(match args.get("measure").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("area") => Measure::Area,
        Some("perimeter") => Measure::Perimeter,
        Some("thickness") => Measure::Thickness,
        Some("centroid_x") => Measure::CentroidX,
        Some("centroid_y") => Measure::CentroidY,
        Some("major_axis") => Measure::MajorAxis,
        Some("minor_axis") => Measure::MinorAxis,
        Some("orientation") => Measure::Orientation,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'measure' must be one of area|perimeter|thickness|centroid_x|centroid_y|major_axis|minor_axis|orientation, got '{o}'"
            )))
        }
    })
}

fn parse_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!("'{key}' must be a boolean"))),
        },
        Some(Value::Number(n)) => Ok(Some(n.as_f64().unwrap_or(0.0) != 0.0)),
        Some(_) => Err(ToolError::Validation(format!("'{key}' must be a boolean"))),
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
    use wbraster::{CrsInfo, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a single-band integer-valued raster from a row-major buffer.
    fn zone_raster(cols: usize, rows: usize, data: &[f64]) -> String {
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
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        ZonalGeometryTool.run(&args, &ctx()).unwrap()
    }

    /// Area = cell count; a 2x3 block of zone 1 has area 6, the rest zone 2.
    #[test]
    fn area_counts_cells() {
        // 3x3: top-left 2x2 is zone 1 (4 cells), rest zone 2 (5 cells).
        let data = [1.0, 1.0, 2.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0];
        let out = run(json!({ "zones": zone_raster(3, 3, &data), "as_table": true }));
        let zones = out.outputs["zones"].as_array().unwrap();
        let z1 = zones.iter().find(|z| z["zone"] == json!(1)).unwrap();
        let z2 = zones.iter().find(|z| z["zone"] == json!(2)).unwrap();
        assert_eq!(z1["area"].as_f64().unwrap(), 4.0);
        assert_eq!(z2["area"].as_f64().unwrap(), 5.0);
    }

    /// Perimeter of a solid 3x3 single-zone raster is its outer boundary = 12.
    #[test]
    fn perimeter_of_full_block() {
        let data = [1.0; 9];
        let out = run(json!({ "zones": zone_raster(3, 3, &data), "as_table": true }));
        let z = &out.outputs["zones"].as_array().unwrap()[0];
        assert_eq!(z["perimeter"].as_f64().unwrap(), 12.0);
    }

    /// An elongated horizontal zone has major>minor axis and ~0° orientation.
    #[test]
    fn orientation_of_horizontal_bar() {
        // 1x5 bar of zone 1 embedded in a 3x5 raster of zone 2.
        let mut data = vec![2.0; 15];
        for c in 0..5 {
            data[5 + c] = 1.0; // middle row
        }
        let out = run(json!({ "zones": zone_raster(5, 3, &data), "as_table": true }));
        let zones = out.outputs["zones"].as_array().unwrap();
        let z1 = zones.iter().find(|z| z["zone"] == json!(1)).unwrap();
        assert!(z1["major_axis"].as_f64().unwrap() > z1["minor_axis"].as_f64().unwrap());
        let ori = z1["orientation"].as_f64().unwrap();
        assert!(
            !(5.0..=175.0).contains(&ori),
            "horizontal bar orientation ~0, got {ori}"
        );
    }

    /// Raster mode paints each zone's cells with the chosen measure.
    #[test]
    fn raster_mode_paints_area() {
        let data = [1.0, 1.0, 2.0, 2.0];
        let out = run(json!({ "zones": zone_raster(2, 2, &data), "measure": "area" }));
        let r = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        // Each zone has 2 cells -> area 2 everywhere.
        for row in 0..2 {
            for col in 0..2 {
                assert_eq!(r.get(0, row, col), 2.0);
            }
        }
    }

    /// Thickness grows with zone width: a 5x5 block is thicker than a 1x5 bar.
    #[test]
    fn thickness_tracks_width() {
        let bar = {
            let mut d = vec![2.0; 15];
            for c in 0..5 {
                d[5 + c] = 1.0;
            }
            run(json!({ "zones": zone_raster(5, 3, &d), "as_table": true }))
        };
        let block = run(json!({ "zones": zone_raster(5, 5, &[1.0; 25]), "as_table": true }));
        let bar_t = bar.outputs["zones"].as_array().unwrap()[0]["thickness"]
            .as_f64()
            .unwrap();
        // block: single zone row
        let block_t = block.outputs["zones"].as_array().unwrap()[0]["thickness"]
            .as_f64()
            .unwrap();
        assert!(
            block_t > bar_t,
            "5x5 block ({block_t}) thicker than 1x5 bar ({bar_t})"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ZonalGeometryTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "zones": "z.tif", "measure": "bogus" })).is_err());
        assert!(bad(json!({ "zones": "z.tif", "measure": "thickness" })).is_ok());
    }
}
