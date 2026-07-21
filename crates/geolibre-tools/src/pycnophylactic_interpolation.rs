//! GeoLibre tool: Tobler's pycnophylactic (mass-preserving) areal interpolation.
//!
//! Pure-Rust counterpart of ArcGIS Pro's areal-interpolation family
//! (Geostatistical Analyst). The GeoLibre `apportion_polygon` transfers
//! attributes by area weighting, which assumes uniform density within each
//! source zone. Pycnophylactic interpolation (Tobler 1979) removes that artifact
//! by producing a *smooth* population surface under the volume-preserving
//! constraint that each zone's cells still sum to its original count. Nothing in
//! the bundled suite does this (the kriging family interpolates point samples,
//! not areal aggregates).
//!
//! The zones are rasterised, each cell seeded with `count / zone_cell_count`, and
//! then iterated: a 3×3 mean smoothing pass followed by a per-zone additive
//! correction that restores each zone's exact total (with non-negativity
//! redistribution), until the maximum change falls below `tolerance` or
//! `iterations` is reached. The output density raster preserves each zone's mass.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
use wbvector::{Geometry, Ring};

use crate::common::parse_optional_output;
use crate::vector_common::load_input_layer;

pub struct PycnophylacticInterpolationTool;

impl Tool for PycnophylacticInterpolationTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "pycnophylactic_interpolation",
            display_name: "Pycnophylactic Interpolation",
            summary: "Tobler's mass-preserving areal interpolation: turn zone-aggregated counts into a smooth density surface whose per-zone sums exactly match the input (like ArcGIS Areal Interpolation) — the smooth alternative to the uniform-density assumption of apportion_polygon.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon layer with a count field (extensive variable).",
                    required: true,
                },
                ToolParamSpec {
                    name: "count_field",
                    description: "Field holding each zone's total count (population, etc.).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output density raster (per-cell value; zone sums equal the input). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Output cell size in CRS units (default: extent / 200).",
                    required: false,
                },
                ToolParamSpec {
                    name: "iterations",
                    description: "Maximum smoothing iterations (default 100).",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Stop when the max per-cell change falls below this fraction of the mean density (default 0.001).",
                    required: false,
                },
                ToolParamSpec {
                    name: "non_negative",
                    description: "Clamp negatives to 0 and redistribute to preserve mass (default true).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "count_field")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_output(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let cidx = layer.schema.field_index(&prm.count_field).ok_or_else(|| {
            ToolError::Validation(format!("count_field '{}' not found", prm.count_field))
        })?;

        // Zones: polygon rings + count + bbox.
        let mut zones: Vec<Zone> = Vec::new();
        let (mut xmin, mut ymin, mut xmax, mut ymax) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for feat in &layer.features {
            let Some(geom) = feat.geometry.as_ref() else {
                continue;
            };
            let polys = polygons(geom);
            if polys.is_empty() {
                continue;
            }
            let count = feat
                .attributes
                .get(cidx)
                .and_then(|v| v.as_f64())
                .unwrap_or(0.0);
            let mut zb = (
                f64::INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
            );
            for (ext, _) in &polys {
                for c in ext {
                    zb.0 = zb.0.min(c.0);
                    zb.1 = zb.1.min(c.1);
                    zb.2 = zb.2.max(c.0);
                    zb.3 = zb.3.max(c.1);
                }
            }
            xmin = xmin.min(zb.0);
            ymin = ymin.min(zb.1);
            xmax = xmax.max(zb.2);
            ymax = ymax.max(zb.3);
            zones.push(Zone {
                polys,
                count,
                bbox: zb,
            });
        }
        if zones.is_empty() {
            return Err(ToolError::Execution(
                "no polygon zones in input".to_string(),
            ));
        }

        let cell = prm
            .cell_size
            .unwrap_or((((xmax - xmin).max(ymax - ymin)) / 200.0).max(1e-9));
        let cols = (((xmax - xmin) / cell).ceil() as usize).max(1);
        let rows = (((ymax - ymin) / cell).ceil() as usize).max(1);
        let n = rows * cols;

        // Rasterise: assign each cell to the first zone whose polygon contains it.
        let mut zone_of = vec![usize::MAX; n];
        for r in 0..rows {
            let cy = ymax - (r as f64 + 0.5) * cell;
            for c in 0..cols {
                let cx = xmin + (c as f64 + 0.5) * cell;
                for (zi, z) in zones.iter().enumerate() {
                    if cx < z.bbox.0 || cx > z.bbox.2 || cy < z.bbox.1 || cy > z.bbox.3 {
                        continue;
                    }
                    if z.contains(cx, cy) {
                        zone_of[r * cols + c] = zi;
                        break;
                    }
                }
            }
        }
        // Cells per zone.
        let mut zn = vec![0usize; zones.len()];
        for &z in &zone_of {
            if z != usize::MAX {
                zn[z] += 1;
            }
        }

        // Seed uniform density.
        let mut grid = vec![0.0f64; n];
        for i in 0..n {
            let z = zone_of[i];
            if z != usize::MAX && zn[z] > 0 {
                grid[i] = zones[z].count / zn[z] as f64;
            }
        }

        ctx.progress.info(&format!(
            "{} zone(s), {rows}x{cols}; smoothing",
            zones.len()
        ));

        let total_mass: f64 = zones.iter().map(|z| z.count).sum();
        let mean_density = if n > 0 { total_mass / n as f64 } else { 0.0 };
        let stop_delta = (prm.tolerance * mean_density.abs()).max(1e-12);

        for it in 0..prm.iterations {
            // 3x3 mean smoothing over valid cells.
            let mut sm = grid.clone();
            for r in 0..rows {
                for c in 0..cols {
                    let idx = r * cols + c;
                    if zone_of[idx] == usize::MAX {
                        continue;
                    }
                    let mut sum = 0.0;
                    let mut cnt = 0.0;
                    for dr in -1i32..=1 {
                        for dc in -1i32..=1 {
                            let nr = r as i32 + dr;
                            let nc = c as i32 + dc;
                            if nr < 0 || nc < 0 || nr >= rows as i32 || nc >= cols as i32 {
                                continue;
                            }
                            let nidx = nr as usize * cols + nc as usize;
                            if zone_of[nidx] != usize::MAX {
                                sum += grid[nidx];
                                cnt += 1.0;
                            }
                        }
                    }
                    if cnt > 0.0 {
                        sm[idx] = sum / cnt;
                    }
                }
            }

            // Per-zone additive correction to restore exact mass.
            let mut zsum = vec![0.0f64; zones.len()];
            for i in 0..n {
                let z = zone_of[i];
                if z != usize::MAX {
                    zsum[z] += sm[i];
                }
            }
            let mut max_delta = 0.0f64;
            for i in 0..n {
                let z = zone_of[i];
                if z == usize::MAX || zn[z] == 0 {
                    continue;
                }
                let correction = (zones[z].count - zsum[z]) / zn[z] as f64;
                let newv = sm[i] + correction;
                max_delta = max_delta.max((newv - grid[i]).abs());
                sm[i] = newv;
            }

            // Non-negativity: clamp and rescale positives to preserve zone mass.
            if prm.non_negative {
                enforce_non_negative(&mut sm, &zone_of, &zones, zn.len());
            }

            grid = sm;
            if max_delta < stop_delta {
                ctx.progress
                    .info(&format!("converged after {} iteration(s)", it + 1));
                break;
            }
        }

        // Final per-zone mass (for the report).
        let mut final_sum = vec![0.0f64; zones.len()];
        for i in 0..n {
            let z = zone_of[i];
            if z != usize::MAX {
                final_sum[z] += grid[i];
            }
        }
        let max_mass_err = zones
            .iter()
            .zip(&final_sum)
            .map(|(z, &s)| (z.count - s).abs())
            .fold(0.0, f64::max);

        // Output raster (nodata outside all zones).
        let nodata = -9999.0;
        let out_data: Vec<f64> = (0..n)
            .map(|i| {
                if zone_of[i] == usize::MAX {
                    nodata
                } else {
                    grid[i]
                }
            })
            .collect();
        let mut out = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: xmin,
            y_min: ymin,
            cell_size: cell,
            cell_size_y: Some(cell),
            nodata,
            data_type: DataType::F32,
            crs: match layer.crs_epsg() {
                Some(e) => CrsInfo {
                    epsg: Some(e),
                    wkt: None,
                    proj4: None,
                },
                None => CrsInfo {
                    epsg: None,
                    wkt: None,
                    proj4: None,
                },
            },
            metadata: Vec::new(),
        });
        for r in 0..rows {
            for c in 0..cols {
                out.set(0, r as isize, c as isize, out_data[r * cols + c])
                    .map_err(|e| ToolError::Execution(format!("write failed: {e}")))?;
            }
        }

        let out_path = crate::common::write_or_store_output(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("zone_count".to_string(), json!(zones.len()));
        outputs.insert("total_mass".to_string(), json!(total_mass));
        outputs.insert("max_zone_mass_error".to_string(), json!(max_mass_err));
        Ok(ToolRunResult { outputs })
    }
}

/// Clamps negatives to 0 and rescales each zone's positive cells so the zone
/// mass is preserved.
fn enforce_non_negative(grid: &mut [f64], zone_of: &[usize], zones: &[Zone], _nzones: usize) {
    // Per zone: target mass, current positive mass.
    let mut target = vec![0.0f64; zones.len()];
    for (zi, z) in zones.iter().enumerate() {
        target[zi] = z.count;
    }
    let mut pos_sum = vec![0.0f64; zones.len()];
    for i in 0..grid.len() {
        let z = zone_of[i];
        if z == usize::MAX {
            continue;
        }
        if grid[i] < 0.0 {
            grid[i] = 0.0;
        } else {
            pos_sum[z] += grid[i];
        }
    }
    for i in 0..grid.len() {
        let z = zone_of[i];
        if z == usize::MAX || grid[i] <= 0.0 {
            continue;
        }
        if pos_sum[z] > 0.0 {
            grid[i] *= target[z] / pos_sum[z];
        }
    }
}

type PolyRings = (Vec<(f64, f64)>, Vec<Vec<(f64, f64)>>);

struct Zone {
    polys: Vec<PolyRings>,
    count: f64,
    bbox: (f64, f64, f64, f64),
}

impl Zone {
    /// Point-in-polygon (even-odd) across all parts, honouring holes.
    fn contains(&self, x: f64, y: f64) -> bool {
        for (ext, holes) in &self.polys {
            if point_in_ring(ext, x, y) && !holes.iter().any(|h| point_in_ring(h, x, y)) {
                return true;
            }
        }
        false
    }
}

fn point_in_ring(ring: &[(f64, f64)], x: f64, y: f64) -> bool {
    let n = ring.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = ring[i];
        let (xj, yj) = ring[j];
        if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi) {
            inside = !inside;
        }
        j = i;
    }
    inside
}

fn polygons(geom: &Geometry) -> Vec<PolyRings> {
    let ring_pts =
        |r: &Ring| -> Vec<(f64, f64)> { r.coords().iter().map(|c| (c.x, c.y)).collect() };
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            vec![(
                ring_pts(exterior),
                interiors.iter().map(&ring_pts).collect(),
            )]
        }
        Geometry::MultiPolygon(parts) => parts
            .iter()
            .map(|(e, h)| (ring_pts(e), h.iter().map(&ring_pts).collect()))
            .collect(),
        _ => Vec::new(),
    }
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

struct Params {
    count_field: String,
    cell_size: Option<f64>,
    iterations: usize,
    tolerance: f64,
    non_negative: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let count_field = require_str(args, "count_field")?.to_string();
    let cell_size = match args.get("cell_size") {
        None | Some(Value::Null) => None,
        Some(Value::Number(n)) => n.as_f64().filter(|v| *v > 0.0),
        Some(Value::String(s)) if s.trim().is_empty() => None,
        Some(Value::String(s)) => Some(
            s.trim()
                .parse::<f64>()
                .map_err(|_| ToolError::Validation("'cell_size' must be a number".into()))?,
        ),
        _ => None,
    };
    let iterations = match args.get("iterations") {
        None | Some(Value::Null) => 100,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(100).max(1) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 100,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'iterations' must be an integer".into()))?
            .max(1),
        _ => 100,
    };
    let tolerance = match args.get("tolerance") {
        None | Some(Value::Null) => 0.001,
        Some(Value::Number(n)) => n.as_f64().unwrap_or(0.001).max(0.0),
        Some(Value::String(s)) if s.trim().is_empty() => 0.001,
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'tolerance' must be a number".into()))?
            .max(0.0),
        _ => 0.001,
    };
    let non_negative = match args.get("non_negative") {
        None | Some(Value::Null) => true,
        Some(Value::Bool(b)) => *b,
        Some(Value::String(s)) => {
            !matches!(s.trim().to_ascii_lowercase().as_str(), "false" | "0" | "no")
        }
        _ => true,
    };
    Ok(Params {
        count_field,
        cell_size,
        iterations,
        tolerance,
        non_negative,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Coord, FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Geometry {
        Geometry::Polygon {
            exterior: Ring::new(vec![
                Coord::xy(x0, y0),
                Coord::xy(x1, y0),
                Coord::xy(x1, y1),
                Coord::xy(x0, y1),
            ]),
            interiors: vec![],
        }
    }

    fn two_zone_layer() -> String {
        let mut l = Layer::new("z")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("pop", FieldType::Float));
        l.add_feature(Some(rect(0.0, 0.0, 10.0, 10.0)), &[("pop", 1000.0.into())])
            .unwrap();
        l.add_feature(Some(rect(10.0, 0.0, 20.0, 10.0)), &[("pop", 100.0.into())])
            .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = PycnophylacticInterpolationTool.run(&args, &ctx()).unwrap();
        let r = crate::common::load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, r)
    }

    /// Mass is preserved: each zone's cell values still sum to its input count.
    #[test]
    fn preserves_zone_mass() {
        let (out, _r) = run(json!({
            "input": two_zone_layer(), "count_field": "pop", "cell_size": 1.0, "iterations": 50,
        }));
        let err = out.outputs["max_zone_mass_error"].as_f64().unwrap();
        assert!(err < 1.0, "zone mass must be preserved, max error {err}");
        assert!((out.outputs["total_mass"].as_f64().unwrap() - 1100.0).abs() < 1e-6);
    }

    /// The surface is smooth: density varies gradually near the shared border
    /// instead of jumping (the whole point vs. uniform apportionment).
    #[test]
    fn produces_a_smooth_surface() {
        let (_out, r) = run(json!({
            "input": two_zone_layer(), "count_field": "pop", "cell_size": 1.0, "iterations": 80,
        }));
        // Sample a row across the shared border at x=10. The high-density zone
        // near the border should be pulled DOWN below its uniform value (10),
        // and the low zone pulled UP above its uniform (1) -> smoothing.
        let mid_row = (r.rows / 2) as isize;
        let read = |cx: f64| {
            let col = ((cx - r.x_min) / r.cell_size_x).floor() as isize;
            r.get(0, mid_row, col)
        };
        let left_interior = read(2.0); // deep in high zone
        let left_border = read(9.0); // high zone near border
        assert!(
            left_border < left_interior,
            "high-density zone should thin toward the low-density neighbour ({left_border} vs {left_interior})"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            PycnophylacticInterpolationTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "count_field": "pop" })).is_ok());
    }
}
