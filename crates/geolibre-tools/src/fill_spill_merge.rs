//! GeoLibre tool: Fill-Spill-Merge — volume-aware depression filling.
//!
//! Distributes a finite quantity of surface water across a DEM to produce
//! realistic standing-water (lake / inundation) extents. Unlike the depression
//! filling and breaching tools bundled from whitebox (which assume unlimited
//! water and raise every pit to its spill point), this tool is *water-volume
//! aware*: water flows downhill into pits; each depression fills only as far as
//! its available water allows; excess **spills** into the neighbouring
//! depression, and adjacent flooded depressions **merge** into larger lakes.
//!
//! Reimplements RichDEM's `FillSpillMerge` (Barnes, Callaghan & Wickert, 2020,
//! *Earth Surface Dynamics* 8, 431–445) from the paper — no RichDEM source is
//! copied (RichDEM is GPL-3; this crate is MIT). The numerics live in
//! [`crate::fill_spill_merge_core`]; this file is the raster-I/O wrapper.
//!
//! Water is supplied either as a uniform depth (`water_level`) applied to every
//! cell, or as a per-cell depth raster (`surface_water`). The grid edge and
//! NoData cells act as free outlets; a `ocean_level` may additionally flood
//! coastal areas at or below that elevation. Outputs are the standing water
//! depth (`output`), the hydraulic surface `dem + wtd` (`water_surface`), and an
//! optional flood-extent mask (`flood_extent`).

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{
    band_to_vec, load_input_raster, parse_optional_output, raster_like_with_data,
    write_or_store_output,
};
use crate::fill_spill_merge_core::fill_spill_merge;

/// Fill-Spill-Merge volume-aware depression filling (RichDEM `FillSpillMerge`).
pub struct FillSpillMergeTool;

impl Tool for FillSpillMergeTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "fill_spill_merge",
            display_name: "Fill-Spill-Merge",
            summary: "Route a finite volume of surface water across a DEM to produce realistic lake/inundation extents (partial depression filling with spill and merge), reimplementing RichDEM's Fill-Spill-Merge.",
            category: ToolCategory::Hydrology,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "dem",
                    description: "Input DEM raster (elevations).",
                    required: true,
                },
                ToolParamSpec {
                    name: "water_level",
                    description: "Uniform surface water depth applied to every cell. Provide this or 'surface_water'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "surface_water",
                    description: "Per-cell surface water depth raster (aligned to the DEM). Overrides 'water_level' where both are given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "ocean_level",
                    description: "Optional sea level: cells at or below this elevation that connect to the grid edge become ocean (coastal inundation baseline).",
                    required: false,
                },
                ToolParamSpec {
                    name: "edge_outlet",
                    description: "Treat every grid-border cell as a free outlet (default true). Water reaching the border leaves the DEM.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional standing-water-depth raster path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "water_surface",
                    description: "Optional output path for the hydraulic surface (DEM + water depth).",
                    required: false,
                },
                ToolParamSpec {
                    name: "flood_extent",
                    description: "Optional output path for a flood-extent mask (1 where standing water is present, else 0).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args.get("dem").and_then(Value::as_str).is_none() {
            return Err(ToolError::Validation(
                "missing required string parameter 'dem'".to_string(),
            ));
        }
        let has_level = parse_optional_f64(args, "water_level")?.is_some();
        let has_raster = parse_optional_output(args, "surface_water")?.is_some();
        if !has_level && !has_raster {
            return Err(ToolError::Validation(
                "provide surface water via 'water_level' (uniform depth) or 'surface_water' (raster)"
                    .to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let dem_path = args
            .get("dem")
            .and_then(Value::as_str)
            .ok_or_else(|| ToolError::Validation("missing required parameter 'dem'".to_string()))?;
        let output = parse_optional_output(args, "output")?;
        let water_surface_out = parse_optional_output(args, "water_surface")?;
        let flood_extent_out = parse_optional_output(args, "flood_extent")?;
        let surface_water_path = parse_optional_output(args, "surface_water")?;
        let water_level = parse_optional_f64(args, "water_level")?;
        let ocean_level = parse_optional_f64(args, "ocean_level")?;
        let edge_outlet = parse_optional_bool(args, "edge_outlet")?.unwrap_or(true);

        let raster = load_input_raster(dem_path)?;
        let rows = raster.rows;
        let cols = raster.cols;
        let nodata = raster.nodata;
        let dem = band_to_vec(&raster, 0);

        // Build the initial surface-water buffer: a per-cell raster if given,
        // otherwise a uniform depth.
        let initial_wtd = if let Some(sw_path) = surface_water_path {
            let sw = load_input_raster(sw_path)?;
            if sw.rows != rows || sw.cols != cols {
                return Err(ToolError::Validation(
                    "'surface_water' raster must match the DEM dimensions".to_string(),
                ));
            }
            let mut buf = band_to_vec(&sw, 0);
            let sw_nodata = sw.nodata;
            for v in buf.iter_mut() {
                if *v == sw_nodata || *v < 0.0 {
                    *v = 0.0;
                }
            }
            buf
        } else {
            let level = water_level.unwrap_or(0.0).max(0.0);
            vec![level; rows * cols]
        };

        ctx.progress.info("building depression hierarchy");
        let result = fill_spill_merge(
            &dem,
            rows,
            cols,
            nodata,
            &initial_wtd,
            ocean_level,
            edge_outlet,
        )
        .map_err(ToolError::Execution)?;
        ctx.progress.progress(0.8);

        let mut outputs = std::collections::BTreeMap::new();

        // Standing water depth (NoData preserved at DEM NoData cells).
        let mut wtd = result.wtd.clone();
        for i in 0..rows * cols {
            if dem[i] == nodata {
                wtd[i] = nodata;
            }
        }
        let wtd_raster = raster_like_with_data(&raster, wtd, nodata, DataType::F32)?;
        let wtd_path = write_or_store_output(wtd_raster, output)?;
        outputs.insert("output".to_string(), json!(wtd_path));

        if let Some(path) = water_surface_out {
            let mut surf = vec![nodata; rows * cols];
            for i in 0..rows * cols {
                if dem[i] != nodata {
                    surf[i] = dem[i] + result.wtd[i];
                }
            }
            let r = raster_like_with_data(&raster, surf, nodata, DataType::F32)?;
            outputs.insert(
                "water_surface".to_string(),
                json!(write_or_store_output(r, Some(path))?),
            );
        }

        if let Some(path) = flood_extent_out {
            let mut mask = vec![0.0f64; rows * cols];
            for i in 0..rows * cols {
                if dem[i] != nodata && result.wtd[i] > 0.0 {
                    mask[i] = 1.0;
                }
            }
            let r = raster_like_with_data(&raster, mask, 0.0, DataType::U8)?;
            outputs.insert(
                "flood_extent".to_string(),
                json!(write_or_store_output(r, Some(path))?),
            );
        }

        // Summary metrics useful to a caller/UI.
        let flooded_cells = result.wtd.iter().filter(|&&v| v > 0.0).count();
        let standing_volume: f64 = result
            .wtd
            .iter()
            .enumerate()
            .filter(|(i, _)| dem[*i] != nodata)
            .map(|(_, &v)| v)
            .sum();
        outputs.insert("flooded_cells".to_string(), json!(flooded_cells));
        outputs.insert("standing_water_volume".to_string(), json!(standing_volume));
        outputs.insert(
            "ocean_outflow_volume".to_string(),
            json!(result.ocean_volume),
        );
        outputs.insert(
            "depression_count".to_string(),
            json!(result.depression_count),
        );

        ctx.progress.progress(1.0);
        Ok(ToolRunResult { outputs })
    }
}

/// Parses an optional numeric parameter that may arrive as a JSON number or a
/// numeric string (host UIs post strings).
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

/// Parses an optional boolean parameter that may arrive as a JSON bool or a
/// string ("true"/"false").
fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
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
    use wbraster::{memory_store, CrsInfo, DataType, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn make_dem() -> Raster {
        // 5x5 bowl with a deep centre; border acts as outlet.
        #[rustfmt::skip]
        let data = vec![
            10.0,10.0,10.0,10.0,10.0,
            10.0, 5.0, 5.0, 5.0,10.0,
            10.0, 5.0, 0.0, 5.0,10.0,
            10.0, 5.0, 5.0, 5.0,10.0,
            10.0,10.0,10.0,10.0,10.0,
        ];
        let cfg = RasterConfig {
            cols: 5,
            rows: 5,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo::from_epsg(3857),
            metadata: Vec::new(),
        };
        Raster::from_data(cfg, data).unwrap()
    }

    fn run_tool(args: ToolArgs) -> ToolRunResult {
        let tool = FillSpillMergeTool;
        tool.validate(&args).expect("validate");
        tool.run(&args, &ctx()).expect("run")
    }

    #[test]
    fn uniform_water_pools_in_depression() {
        let id = memory_store::put_raster(make_dem());
        let dem_path = memory_store::make_raster_memory_path(&id);
        let mut args = ToolArgs::new();
        args.insert("dem".to_string(), json!(dem_path));
        args.insert("water_level".to_string(), json!(0.3));
        let out = run_tool(args);
        let flooded = out.outputs.get("flooded_cells").unwrap().as_u64().unwrap();
        assert!(flooded > 0, "water should pool in the depression");
        let standing = out
            .outputs
            .get("standing_water_volume")
            .unwrap()
            .as_f64()
            .unwrap();
        assert!(standing > 0.0);
    }

    #[test]
    fn excess_water_flows_to_ocean() {
        let id = memory_store::put_raster(make_dem());
        let dem_path = memory_store::make_raster_memory_path(&id);
        let mut args = ToolArgs::new();
        args.insert("dem".to_string(), json!(dem_path));
        args.insert("water_level".to_string(), json!(100.0));
        let out = run_tool(args);
        let ocean = out
            .outputs
            .get("ocean_outflow_volume")
            .unwrap()
            .as_f64()
            .unwrap();
        assert!(ocean > 0.0, "excess water should drain to the ocean");
    }

    #[test]
    fn rejects_missing_water() {
        let id = memory_store::put_raster(make_dem());
        let dem_path = memory_store::make_raster_memory_path(&id);
        let mut args = ToolArgs::new();
        args.insert("dem".to_string(), json!(dem_path));
        let tool = FillSpillMergeTool;
        assert!(tool.validate(&args).is_err());
    }

    #[test]
    fn rejects_missing_dem() {
        let tool = FillSpillMergeTool;
        let mut args = ToolArgs::new();
        args.insert("water_level".to_string(), json!(1.0));
        assert!(tool.validate(&args).is_err());
    }
}
