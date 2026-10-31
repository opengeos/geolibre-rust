//! GeoLibre tool: fill each zone with the minimum weight value on its boundary.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Zonal Fill* (Spatial Analyst). For each
//! zone in a categorical raster it finds the smallest value of a co-registered
//! weight raster among that zone's boundary cells, then writes that minimum across
//! every cell of the zone. It is a building block for constrained-surface work
//! (flooding each zone up to the lowest point on its rim), complementing the
//! repo's depression tooling (`fill_spill_merge` #291, `storage_capacity`).
//!
//! No equivalent exists in the bundled whitebox suite: `fill` fills DEM pits by
//! flow, not per-zone boundary minima, and `zonal_statistics` summarizes values
//! but does not paint a boundary-min back into the zone.
//!
//! A boundary cell of a zone is one whose 4-neighbourhood contains a different
//! zone, a no-data cell, or the grid edge. Zones follow ArcGIS semantics: each
//! distinct integer value in the zone raster is one zone (need not be contiguous).

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{load_input_raster, parse_optional_output, raster_like_with_data};

pub struct ZonalFillTool;

impl Tool for ZonalFillTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "zonal_fill",
            display_name: "Zonal Fill",
            summary: "Fill each zone of a categorical raster with the minimum value of a weight raster found along that zone's boundary (like ArcGIS Zonal Fill) — a constrained-surface building block absent from the bundled suite (fill is flow-based, zonal_statistics only summarizes).",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "zones",
                    description: "Categorical zone raster (each distinct integer value is one zone).",
                    required: true,
                },
                ToolParamSpec {
                    name: "weight",
                    description: "Co-registered weight raster whose boundary minimum fills each zone.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output raster. If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "zones")?;
        require_str(args, "weight")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let zones_path = require_str(args, "zones")?;
        let weight_path = require_str(args, "weight")?;
        let output = parse_optional_output(args, "output")?;

        let zr = load_input_raster(zones_path)?;
        let wr = load_input_raster(weight_path)?;
        let (rows, cols) = (zr.rows, zr.cols);
        if wr.rows != rows || wr.cols != cols {
            return Err(ToolError::Validation(format!(
                "weight raster is {}x{}, expected {rows}x{cols} to match zones",
                wr.rows, wr.cols
            )));
        }
        let aligned = |a: f64, b: f64| (a - b).abs() <= 1e-6 * a.abs().max(b.abs()).max(1.0);
        if !aligned(wr.x_min, zr.x_min)
            || !aligned(wr.y_min, zr.y_min)
            || !aligned(wr.cell_size_x, zr.cell_size_x)
            || !aligned(wr.cell_size_y, zr.cell_size_y)
        {
            return Err(ToolError::Validation(
                "weight raster is not co-registered with the zone raster (origin/resolution differ)"
                    .to_string(),
            ));
        }

        // Read zone labels.
        let zone_at = |r: usize, c: usize| -> Option<i64> {
            let v = zr.get(0, r as isize, c as isize);
            if v != zr.nodata && v.is_finite() {
                Some(v.round() as i64)
            } else {
                None
            }
        };

        // Per-zone boundary minimum of the weight raster.
        let mut zmin: BTreeMap<i64, f64> = BTreeMap::new();
        for r in 0..rows {
            for c in 0..cols {
                let Some(z) = zone_at(r, c) else { continue };
                let mut boundary = false;
                for (dr, dc) in [(-1isize, 0isize), (1, 0), (0, -1), (0, 1)] {
                    let nr = r as isize + dr;
                    let nc = c as isize + dc;
                    let off_grid =
                        nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize;
                    if off_grid || zone_at(nr as usize, nc as usize) != Some(z) {
                        boundary = true;
                    }
                }
                if !boundary {
                    continue;
                }
                let w = wr.get(0, r as isize, c as isize);
                if w == wr.nodata || !w.is_finite() {
                    continue;
                }
                let e = zmin.entry(z).or_insert(f64::INFINITY);
                if w < *e {
                    *e = w;
                }
            }
        }

        if zmin.is_empty() {
            return Err(ToolError::Execution(
                "no zone had a valid weight value on its boundary".to_string(),
            ));
        }
        ctx.progress.info(&format!("filled {} zone(s)", zmin.len()));

        let nodata = -9999.0_f64;
        let mut data = vec![nodata; rows * cols];
        let mut zones_without_boundary = 0usize;
        for r in 0..rows {
            for c in 0..cols {
                let Some(z) = zone_at(r, c) else { continue };
                match zmin.get(&z) {
                    Some(v) if v.is_finite() => data[r * cols + c] = *v,
                    _ => zones_without_boundary += 1,
                }
            }
        }

        let out_r = raster_like_with_data(&zr, data, nodata, DataType::F32)?;
        let out_path = crate::common::write_or_store_output(out_r, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("zone_count".to_string(), json!(zmin.len()));
        outputs.insert(
            "cells_without_valid_boundary".to_string(),
            json!(zones_without_boundary),
        );
        Ok(ToolRunResult { outputs })
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

    fn raster(cols: usize, rows: usize, data: &[f64]) -> String {
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

    fn run(z: String, w: String) -> Raster {
        let args: ToolArgs = serde_json::from_value(json!({ "zones": z, "weight": w })).unwrap();
        let out = ZonalFillTool.run(&args, &ctx()).unwrap();
        load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    /// Every cell of a zone gets the minimum weight on that zone's boundary. In a
    /// single-zone 3x3 raster all cells are boundary cells, so fill = global min.
    #[test]
    fn single_zone_fills_global_min() {
        let z = raster(3, 3, &[1.0; 9]);
        let w = raster(3, 3, &[5.0, 3.0, 9.0, 7.0, 1.0, 8.0, 6.0, 4.0, 2.0]);
        let out = run(z, w);
        // The center (value 1) is interior, not boundary; boundary min = 2.
        for row in 0..3 {
            for col in 0..3 {
                assert_eq!(out.get(0, row, col), 2.0, "cell ({row},{col})");
            }
        }
    }

    /// Two side-by-side zones each fill to their own boundary minimum.
    #[test]
    fn two_zones_independent() {
        // 2x2: left column zone 1, right column zone 2.
        let z = raster(2, 2, &[1.0, 2.0, 1.0, 2.0]);
        let w = raster(2, 2, &[10.0, 20.0, 30.0, 40.0]);
        let out = run(z, w);
        // Zone 1 boundary min = 10, zone 2 = 20 (all cells are boundary here).
        assert_eq!(out.get(0, 0, 0), 10.0);
        assert_eq!(out.get(0, 1, 0), 10.0);
        assert_eq!(out.get(0, 0, 1), 20.0);
        assert_eq!(out.get(0, 1, 1), 20.0);
    }

    #[test]
    fn rejects_misregistered() {
        let z = raster(2, 2, &[1.0; 4]);
        let w = raster(3, 3, &[1.0; 9]);
        let args: ToolArgs = serde_json::from_value(json!({ "zones": z, "weight": w })).unwrap();
        assert!(ZonalFillTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ZonalFillTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "zones": "z.tif" })).is_err());
        assert!(bad(json!({ "zones": "z.tif", "weight": "w.tif" })).is_ok());
    }
}
