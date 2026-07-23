//! GeoLibre tool: enforce a monotonically non-increasing elevation profile
//! downstream along mapped river polylines.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Enforce River Monotonicity*
//! (3D Analyst). Raw DEMs contain small rises along river channels (sensor
//! noise, bridges, overhanging vegetation) that break the physical assumption
//! that water flows downhill. Hydraulic / flood models need channel profiles
//! that never rise going downstream.
//!
//! The bundled `fill_depressions` / `fill_burn` remove *closed* pits but do
//! **not** enforce a monotone profile along a known river vector. This tool
//! reads a river polyline layer plus a DEM, and for each line:
//!
//! 1. Samples DEM Z (bilinear) at the two endpoints and orients the line so we
//!    traverse from the higher-Z endpoint to the lower-Z endpoint — downstream.
//! 2. Densifies the oriented line to `sample_distance` (≈ one point per cell).
//! 3. Marches downstream carrying `running_min`, the current enforced elevation
//!    ceiling. Each accepted step drops the ceiling by at least `tolerance`
//!    (default 0 = flat-allowed, non-increasing) and follows the terrain down
//!    wherever it naturally descends faster. Each visited cell is *carved* to
//!    `min(existing, ceiling)` — never raised.
//!
//! Only cells under the river are ever touched; the rest of the DEM is copied
//! through unchanged. The result is a corrected DEM whose sampled profile is
//! monotonically non-increasing (strictly decreasing if `tolerance` > 0) from
//! headwater to outlet along every river line. The vector and raster must
//! share a CRS.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;
use wbvector::{Coord, Geometry};

use crate::common::{load_input_raster, write_or_store_output};
use crate::vector_common::{load_input_layer, parse_optional_str};

pub struct EnforceRiverMonotonicityTool;

impl Tool for EnforceRiverMonotonicityTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "enforce_river_monotonicity",
            display_name: "Enforce River Monotonicity",
            summary: "Carve a DEM so elevations decrease monotonically downstream along mapped river polylines — a hydro-enforcement step for hydraulic and flood modeling, like ArcGIS Enforce River Monotonicity. Unlike fill_depressions/fill_burn (which remove closed pits), this enforces a monotone profile along a known river vector.",
            category: ToolCategory::Hydrology,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input river polyline layer (LineString / MultiLineString) sharing the DEM's CRS.",
                    required: true,
                },
                ToolParamSpec {
                    name: "surface",
                    description: "Elevation raster (DEM) to correct.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output corrected DEM path (format from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Minimum elevation drop enforced per densified step, in Z units (default 0 = non-increasing). Positive values force a strictly decreasing profile.",
                    required: false,
                },
                ToolParamSpec {
                    name: "sample_distance",
                    description: "Densification interval in CRS units along each river line. Default: the raster cell size (≈ one sample per cell).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based DEM band to correct (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "surface")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let surface = require_str(args, "surface")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        ctx.progress.info("reading DEM and river lines");
        let mut dem = load_input_raster(surface)?;
        let layer = load_input_layer(input)?;

        let step = prm
            .sample_distance
            .unwrap_or_else(|| dem.cell_size_x.min(dem.cell_size_y))
            .max(f64::MIN_POSITIVE);

        // Collect every line part (LineString + each part of a MultiLineString).
        let mut parts: Vec<Vec<Coord>> = Vec::new();
        for feature in &layer.features {
            match feature.geometry.as_ref() {
                Some(Geometry::LineString(cs)) => {
                    if cs.len() >= 2 {
                        parts.push(cs.clone());
                    }
                }
                Some(Geometry::MultiLineString(lines)) => {
                    for l in lines {
                        if l.len() >= 2 {
                            parts.push(l.clone());
                        }
                    }
                }
                _ => {}
            }
        }

        ctx.progress.info(&format!(
            "enforcing monotonicity on {} river line(s)",
            parts.len()
        ));

        let mut rivers_processed = 0usize;
        let mut cells_lowered = 0usize;
        let mut max_drop = 0.0f64;
        let mut total_drop = 0.0f64;

        for coords in &parts {
            if enforce_one_line(
                &mut dem,
                coords,
                prm.band,
                step,
                prm.tolerance,
                &mut cells_lowered,
                &mut max_drop,
                &mut total_drop,
            ) {
                rivers_processed += 1;
            }
        }

        ctx.progress.info(&format!(
            "{rivers_processed} river(s) processed, {cells_lowered} cell(s) lowered"
        ));

        let out_path = write_or_store_output(dem, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("rivers_processed".to_string(), json!(rivers_processed));
        outputs.insert("cells_lowered".to_string(), json!(cells_lowered));
        outputs.insert("max_drop".to_string(), json!(max_drop));
        outputs.insert("total_drop".to_string(), json!(total_drop));
        Ok(ToolRunResult { outputs })
    }
}

/// Carves one river line into the DEM. Returns true if at least one valid
/// sample was marched along the line (so it counts as "processed").
#[allow(clippy::too_many_arguments)]
fn enforce_one_line(
    dem: &mut Raster,
    coords: &[Coord],
    band: isize,
    step: f64,
    tolerance: f64,
    cells_lowered: &mut usize,
    max_drop: &mut f64,
    total_drop: &mut f64,
) -> bool {
    // Orient headwater→outlet: traverse from the higher-Z endpoint to the
    // lower-Z endpoint so a line drawn either way yields the same correction.
    let z_first = endpoint_z(dem, band, coords.first());
    let z_last = endpoint_z(dem, band, coords.last());
    let reversed = matches!((z_first, z_last), (Some(a), Some(b)) if a < b);

    let dense = densify(coords, step);
    let iter: Box<dyn Iterator<Item = &Coord>> = if reversed {
        Box::new(dense.iter().rev())
    } else {
        Box::new(dense.iter())
    };

    let mut running_min: Option<f64> = None;
    let mut marched = false;

    for p in iter {
        // Skip points whose bilinear sample is nodata / off-raster; this does
        // not break the running_min chain — we simply move on.
        let Some(z) = sample_bilinear(dem, band, p.x, p.y) else {
            continue;
        };
        marched = true;

        let ceiling = match running_min {
            None => z,
            Some(rm) => (rm - tolerance).min(z),
        };
        running_min = Some(ceiling);

        // Carve the DEM cell under this point: never raise it.
        let Some((col, row)) = dem.world_to_pixel(p.x, p.y) else {
            continue;
        };
        let existing = dem.get(band, row, col);
        if existing == dem.nodata || existing.is_nan() {
            continue;
        }
        if ceiling < existing {
            let drop = existing - ceiling;
            let _ = dem.set(band, row, col, ceiling);
            *cells_lowered += 1;
            *total_drop += drop;
            if drop > *max_drop {
                *max_drop = drop;
            }
        }
    }

    marched
}

fn endpoint_z(dem: &Raster, band: isize, c: Option<&Coord>) -> Option<f64> {
    c.and_then(|c| sample_bilinear(dem, band, c.x, c.y))
}

/// Densifies a coordinate chain so no segment is longer than `max_len`.
fn densify(coords: &[Coord], max_len: f64) -> Vec<Coord> {
    let n = coords.len();
    if n < 2 || max_len <= 0.0 {
        return coords.to_vec();
    }
    let mut out = Vec::with_capacity(n * 2);
    for i in 0..n - 1 {
        let (ax, ay) = (coords[i].x, coords[i].y);
        let (bx, by) = (coords[i + 1].x, coords[i + 1].y);
        out.push(Coord::xy(ax, ay));
        let d = (bx - ax).hypot(by - ay);
        let pieces = (d / max_len).ceil().max(1.0) as usize;
        for j in 1..pieces {
            let t = j as f64 / pieces as f64;
            out.push(Coord::xy(ax + (bx - ax) * t, ay + (by - ay) * t));
        }
    }
    out.push(Coord::xy(coords[n - 1].x, coords[n - 1].y));
    out
}

// ── Surface sampling (bilinear) ───────────────────────────────────────────────

fn cell_value(dem: &Raster, band: isize, row: isize, col: isize) -> Option<f64> {
    if row < 0 || col < 0 || row >= dem.rows as isize || col >= dem.cols as isize {
        return None;
    }
    let v = dem.get(band, row, col);
    if v == dem.nodata || v.is_nan() {
        None
    } else {
        Some(v)
    }
}

fn sample_bilinear(dem: &Raster, band: isize, x: f64, y: f64) -> Option<f64> {
    let fx = (x - dem.x_min) / dem.cell_size_x - 0.5;
    let fy = (dem.y_max() - y) / dem.cell_size_y - 0.5;
    let col0 = fx.floor() as isize;
    let row0 = fy.floor() as isize;
    let tx = fx - col0 as f64;
    let ty = fy - row0 as f64;

    let v00 = cell_value(dem, band, row0, col0);
    let v01 = cell_value(dem, band, row0, col0 + 1);
    let v10 = cell_value(dem, band, row0 + 1, col0);
    let v11 = cell_value(dem, band, row0 + 1, col0 + 1);

    if let (Some(a), Some(b), Some(c), Some(d)) = (v00, v01, v10, v11) {
        let top = a * (1.0 - tx) + b * tx;
        let bot = c * (1.0 - tx) + d * tx;
        return Some(top * (1.0 - ty) + bot * ty);
    }
    // Near a nodata/edge: fall back to the nearest valid corner by weight.
    let candidates = [
        (v00, (1.0 - tx) * (1.0 - ty)),
        (v01, tx * (1.0 - ty)),
        (v10, (1.0 - tx) * ty),
        (v11, tx * ty),
    ];
    candidates
        .iter()
        .filter_map(|(v, w)| v.map(|v| (v, *w)))
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(v, _)| v)
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    tolerance: f64,
    sample_distance: Option<f64>,
    band: isize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let tolerance = parse_optional_f64(args, "tolerance")?.unwrap_or(0.0);
    if !(tolerance >= 0.0 && tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "'tolerance' must be a finite number >= 0".to_string(),
        ));
    }

    let sample_distance = parse_optional_f64(args, "sample_distance")?;
    if let Some(v) = sample_distance {
        if !(v > 0.0 && v.is_finite()) {
            return Err(ToolError::Validation(
                "'sample_distance' must be a positive number".to_string(),
            ));
        }
    }

    let band_1based = parse_optional_f64(args, "band")?
        .map(|v| v as i64)
        .unwrap_or(1);
    if band_1based < 1 {
        return Err(ToolError::Validation("'band' must be >= 1".to_string()));
    }

    Ok(Params {
        tolerance,
        sample_distance,
        band: (band_1based - 1) as isize,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
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

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
    use wbvector::{memory_store, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// A DEM descending west→east (z = 100 - 10*col) with a raised "bump" cell
    /// at (bump_row, bump_col) that breaks the downstream profile.
    fn bump_dem(cols: usize, rows: usize, bump_row: usize, bump_col: usize, bump_z: f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
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
                r.set(0, row as isize, col as isize, 100.0 - 10.0 * col as f64)
                    .unwrap();
            }
        }
        r.set(0, bump_row as isize, bump_col as isize, bump_z)
            .unwrap();
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn line_layer(coords: &[(f64, f64)]) -> String {
        let mut l = Layer::new("rivers")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        let cs = coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
        l.add_feature(Some(Geometry::line_string(cs)), &[]).unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (Raster, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = EnforceRiverMonotonicityTool.run(&args, &ctx()).unwrap();
        let dem = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (dem, out)
    }

    #[test]
    fn carves_bump_and_profile_is_monotone() {
        // 10x5 DEM, bump at row 2 col 5 raised to 200. River across row 2.
        let dem_path = bump_dem(10, 5, 2, 5, 200.0);
        let input = line_layer(&[(0.5, 2.5), (9.5, 2.5)]);
        let (dem, out) = run(json!({ "input": input, "surface": dem_path }));

        // The bump cell was lowered.
        assert!(
            out.outputs["cells_lowered"].as_u64().unwrap() >= 1,
            "expected at least one lowered cell"
        );
        assert!(out.outputs["max_drop"].as_f64().unwrap() > 100.0);
        assert_eq!(out.outputs["rivers_processed"].as_u64().unwrap(), 1);

        // Sampled profile along the river (row 2, cols 0..9) is non-increasing.
        let mut prev = f64::INFINITY;
        for col in 0..10 {
            let z = dem.get(0, 2, col);
            assert!(z <= prev + 1e-9, "profile rose at col {col}: {z} > {prev}");
            prev = z;
        }
        // The bump cell (was 200) is now no higher than its upstream neighbor.
        assert!(dem.get(0, 2, 5) <= dem.get(0, 2, 4) + 1e-9);
    }

    #[test]
    fn cells_off_the_river_are_unchanged() {
        let dem_path = bump_dem(10, 5, 2, 5, 200.0);
        let input = line_layer(&[(0.5, 2.5), (9.5, 2.5)]);
        let (dem, _) = run(json!({ "input": input, "surface": dem_path }));
        // Rows other than the river row keep the original descending values.
        for row in [0, 1, 3, 4] {
            for col in 0..10 {
                let expected = 100.0 - 10.0 * col as f64;
                assert!(
                    (dem.get(0, row, col) - expected).abs() < 1e-6,
                    "off-river cell ({row},{col}) changed"
                );
            }
        }
    }

    #[test]
    fn orientation_is_endpoint_z_independent() {
        // Same river drawn upstream→downstream and reversed → identical DEM.
        let dem_a = bump_dem(10, 5, 2, 5, 200.0);
        let dem_b = bump_dem(10, 5, 2, 5, 200.0);
        let fwd = line_layer(&[(0.5, 2.5), (9.5, 2.5)]);
        let rev = line_layer(&[(9.5, 2.5), (0.5, 2.5)]);
        let (a, _) = run(json!({ "input": fwd, "surface": dem_a }));
        let (b, _) = run(json!({ "input": rev, "surface": dem_b }));
        for row in 0..5 {
            for col in 0..10 {
                assert!(
                    (a.get(0, row, col) - b.get(0, row, col)).abs() < 1e-9,
                    "orientation changed cell ({row},{col})"
                );
            }
        }
    }

    #[test]
    fn positive_tolerance_forces_strict_decrease() {
        let dem_path = bump_dem(10, 5, 2, 5, 200.0);
        let input = line_layer(&[(0.5, 2.5), (9.5, 2.5)]);
        let (dem, _) = run(json!({ "input": input, "surface": dem_path, "tolerance": 1.0 }));
        let mut prev = f64::INFINITY;
        for col in 0..10 {
            let z = dem.get(0, 2, col);
            assert!(z < prev + 1e-9, "profile not decreasing at col {col}");
            prev = z;
        }
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            EnforceRiverMonotonicityTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "surface": "d.tif" })).is_err());
        assert!(bad(json!({ "input": "r.geojson" })).is_err());
        assert!(bad(json!({ "input": "r.geojson", "surface": "d.tif", "tolerance": -1 })).is_err());
        assert!(
            bad(json!({ "input": "r.geojson", "surface": "d.tif", "sample_distance": 0 })).is_err()
        );
        assert!(bad(json!({ "input": "r.geojson", "surface": "d.tif" })).is_ok());
    }
}
