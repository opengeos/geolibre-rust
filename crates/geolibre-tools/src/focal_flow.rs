//! GeoLibre tool: 8-neighbour inflow bitmask from a surface raster.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Focal Flow* (Spatial Analyst).
//!
//! The bundled hydrology suite is all *downslope* routing — `d8_pointer`,
//! `d8_flow_accum`, `dinf_flow_accum`, `fd8_flow_accum`, `mdinf_flow_accum`,
//! `rho8_flow_accum` — every one of which answers "where does this cell send
//! its water". Focal Flow answers the inverse, local question: *which of my
//! eight neighbours drain into me*. That is not a D8 pointer read backwards,
//! because it is not limited to a single steepest direction: **every** neighbour
//! above the threshold is recorded, so convergent flow into pits, sinks and
//! saddles is captured where D8 collapses it to one arbitrary path.
//!
//! Each output cell is an 8-bit mask. Starting at the east neighbour and
//! proceeding **counter-clockwise**, the neighbours take weights 1, 2, 4, 8, 16,
//! 32, 64, 128, summed into one value. `0` means no neighbour flows in (a local
//! high); `255` means all eight do (a pit). The encoding matches ArcGIS so
//! results are directly comparable.

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};

/// Neighbour offsets and their bit weights, counter-clockwise from east.
///
/// `(d_row, d_col, weight)` with row increasing **downward** (north-up raster),
/// so "north" is `d_row = -1`. Counter-clockwise from east in map space is
/// therefore E, NE, N, NW, W, SW, S, SE.
const NEIGHBORS: [(isize, isize, u16); 8] = [
    (0, 1, 1),   // E
    (-1, 1, 2),  // NE
    (-1, 0, 4),  // N
    (-1, -1, 8), // NW
    (0, -1, 16), // W
    (1, -1, 32), // SW
    (1, 0, 64),  // S
    (1, 1, 128), // SE
];

/// Encodes, per cell, which of its eight neighbours flow into it.
pub struct FocalFlowTool;

impl Tool for FocalFlowTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "focal_flow",
            display_name: "Focal Flow",
            summary: "Encodes which of a cell's eight neighbours flow into it as an 8-bit mask (ArcGIS Focal Flow). Unlike the bundled D8/D-infinity accumulation tools, which route each cell downslope along one direction, this records every inflowing neighbour — capturing convergence into pits, sinks and saddles. Bit weights run counter-clockwise from east (1, 2, 4, 8, 16, 32, 64, 128); 0 is a local high, 255 a pit.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input surface raster (elevation or any continuous surface).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output raster path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "threshold",
                    description: "A neighbour counts as flowing in only if it exceeds the centre cell by more than this (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args.get("input").and_then(Value::as_str).is_none() {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        parse_threshold(args)?;
        parse_band(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args.get("input").and_then(Value::as_str).ok_or_else(|| {
            ToolError::Validation("missing required parameter 'input'".to_string())
        })?;
        let output = parse_optional_output(args, "output")?;
        let threshold = parse_threshold(args)?.unwrap_or(0.0);
        let band_1based = parse_band(args)?.unwrap_or(1);

        let raster = load_input_raster(input)?;
        if band_1based == 0 || band_1based as usize > raster.bands {
            return Err(ToolError::Validation(format!(
                "band {band_1based} out of range (raster has {} band(s))",
                raster.bands
            )));
        }
        let band = (band_1based - 1) as isize;

        let rows = raster.rows;
        let cols = raster.cols;
        let nodata = raster.nodata;
        // The mask is 0..=255, so -1 is a safe out-of-range no-data marker.
        let out_nodata = -1.0_f64;

        ctx.progress.info("computing focal flow");
        let mut data = vec![out_nodata; rows * cols];
        let mut pits = 0_u64;
        let mut peaks = 0_u64;
        let mut valid = 0_u64;

        for row in 0..rows {
            for col in 0..cols {
                let centre = raster.get(band, row as isize, col as isize);
                if centre == nodata || !centre.is_finite() {
                    continue;
                }
                let mut mask = 0_u16;
                for (d_row, d_col, weight) in NEIGHBORS {
                    let n_row = row as isize + d_row;
                    let n_col = col as isize + d_col;
                    // Off-grid neighbours contribute no bit, matching ArcGIS.
                    if n_row < 0 || n_col < 0 || n_row >= rows as isize || n_col >= cols as isize {
                        continue;
                    }
                    let neighbour = raster.get(band, n_row, n_col);
                    if neighbour == nodata || !neighbour.is_finite() {
                        continue;
                    }
                    if neighbour - centre > threshold {
                        mask |= weight;
                    }
                }
                data[row * cols + col] = mask as f64;
                valid += 1;
                if mask == 255 {
                    pits += 1;
                } else if mask == 0 {
                    peaks += 1;
                }
            }
            ctx.progress
                .progress((row as f64 + 1.0) / rows.max(1) as f64);
        }

        let out = raster_like_with_data(&raster, data, out_nodata, DataType::I16)?;
        let out_path = write_or_store_output(out, output)?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("valid_cells".to_string(), json!(valid));
        outputs.insert("pit_cells".to_string(), json!(pits));
        outputs.insert("peak_cells".to_string(), json!(peaks));
        Ok(ToolRunResult { outputs })
    }
}

fn parse_threshold(args: &ToolArgs) -> Result<Option<f64>, ToolError> {
    match args.get("threshold") {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s.trim().parse::<f64>().map(Some).map_err(|_| {
            ToolError::Validation("parameter 'threshold' must be a number".to_string())
        }),
        Some(_) => Err(ToolError::Validation(
            "parameter 'threshold' must be a number".to_string(),
        )),
    }
}

fn parse_band(args: &ToolArgs) -> Result<Option<u64>, ToolError> {
    match args.get("band") {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => n.as_u64().map(Some).ok_or_else(|| {
            ToolError::Validation("parameter 'band' must be a positive integer".to_string())
        }),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s.trim().parse::<u64>().map(Some).map_err(|_| {
            ToolError::Validation("parameter 'band' must be a positive integer".to_string())
        }),
        Some(_) => Err(ToolError::Validation(
            "parameter 'band' must be a positive integer".to_string(),
        )),
    }
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

    fn run_with(path: String, extra: Value) -> Raster {
        let mut obj = serde_json::Map::new();
        obj.insert("input".to_string(), json!(path));
        if let Value::Object(m) = extra {
            for (k, v) in m {
                obj.insert(k, v);
            }
        }
        let args: ToolArgs = serde_json::from_value(Value::Object(obj)).unwrap();
        let out = FocalFlowTool.run(&args, &ctx()).unwrap();
        load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    /// A cell lower than all eight neighbours is a pit: every bit set.
    #[test]
    fn pit_centre_has_all_bits_set() {
        let path = raster(3, 3, &[10.0, 10.0, 10.0, 10.0, 1.0, 10.0, 10.0, 10.0, 10.0]);
        let out = run_with(path, json!({}));
        assert_eq!(out.get(0, 1, 1), 255.0);
    }

    /// A cell higher than all eight neighbours is a local high: no bits set.
    #[test]
    fn peak_centre_has_no_bits_set() {
        let path = raster(3, 3, &[1.0, 1.0, 1.0, 1.0, 10.0, 1.0, 1.0, 1.0, 1.0]);
        let out = run_with(path, json!({}));
        assert_eq!(out.get(0, 1, 1), 0.0);
    }

    /// On a plane tilting up toward the west, only the three western neighbours
    /// (NW=8, W=16, SW=32) are above the centre => 8 + 16 + 32 = 56.
    #[test]
    fn planar_surface_sets_upslope_bits_only() {
        // Column 0 is highest, column 2 lowest; rows identical.
        let path = raster(3, 3, &[3.0, 2.0, 1.0, 3.0, 2.0, 1.0, 3.0, 2.0, 1.0]);
        let out = run_with(path, json!({}));
        assert_eq!(out.get(0, 1, 1), 56.0);
    }

    /// The threshold suppresses shallow gradients: a 1-unit rise does not count
    /// once the threshold reaches 1.
    #[test]
    fn threshold_suppresses_shallow_gradient() {
        let path = raster(3, 3, &[10.0, 10.0, 10.0, 10.0, 9.0, 10.0, 10.0, 10.0, 10.0]);
        let strict = run_with(path.clone(), json!({ "threshold": 1.0 }));
        assert_eq!(
            strict.get(0, 1, 1),
            0.0,
            "1-unit rise must not exceed threshold 1"
        );
        let loose = run_with(path, json!({ "threshold": 0.5 }));
        assert_eq!(loose.get(0, 1, 1), 255.0);
    }

    /// No-data neighbours contribute no bit, and a no-data centre stays no-data.
    #[test]
    fn nodata_is_excluded_and_propagated() {
        let path = raster(
            3,
            3,
            &[10.0, 10.0, 10.0, 10.0, 1.0, -9999.0, 10.0, 10.0, 10.0],
        );
        let out = run_with(path, json!({}));
        // East neighbour (weight 1) is no-data => 255 - 1 = 254.
        assert_eq!(out.get(0, 1, 1), 254.0);
        // The no-data cell itself is written as the output no-data marker.
        assert_eq!(out.get(0, 1, 2), out.nodata);
    }

    /// Edge cells simply have fewer neighbours; off-grid ones contribute nothing.
    #[test]
    fn edge_cells_ignore_offgrid_neighbours() {
        let path = raster(3, 3, &[1.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]);
        let out = run_with(path, json!({}));
        // Corner (0,0) has only E (1), S (64) and SE (128) in range => 193.
        assert_eq!(out.get(0, 0, 0), 193.0);
    }

    #[test]
    fn rejects_bad_parameters() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(FocalFlowTool.validate(&args).is_err());

        let path = raster(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": path.clone(), "threshold": "not-a-number" }))
                .unwrap();
        assert!(FocalFlowTool.validate(&args).is_err());

        // An out-of-range band is caught at run time, before the `band - 1`
        // subtraction that would otherwise underflow.
        for bad_band in [0, 2] {
            let args: ToolArgs =
                serde_json::from_value(json!({ "input": path.clone(), "band": bad_band })).unwrap();
            assert!(
                FocalFlowTool.run(&args, &ctx()).is_err(),
                "band {bad_band} on a single-band raster must be rejected"
            );
        }
    }
}
