//! GeoLibre tool: grow or shrink selected classes of a categorical raster.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Expand* and *Shrink* (Spatial
//! Analyst). The bundled morphology works on binary/greyscale rasters
//! (`buffer_raster`, opening/closing) and `nibble` fills masked areas — but none
//! grow or shrink **specific classes** of a categorical raster while leaving the
//! others intact. This strengthens the raster→vector cleanup pipeline
//! (expand/shrink → `polygonize` → `regularize_building_footprints` /
//! `smooth_natural_features`).
//!
//! `expand` dilates the selected classes by `cells` cells: each iteration, a
//! non-selected cell that touches a selected cell (8-connected) adopts the most
//! frequent selected class among its neighbours. `shrink` erodes them: each
//! iteration, a selected cell that touches a non-selected cell adopts the most
//! frequent non-selected class among its neighbours. No-data cells are barriers
//! (never grown into, never used as a shrink target). Ties go to the smaller
//! class value. All other cells keep their value.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{load_input_raster, parse_optional_output, raster_like_with_data};

pub struct ExpandShrinkTool;

impl Tool for ExpandShrinkTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "expand_shrink",
            display_name: "Expand / Shrink",
            summary: "Grow (expand) or shrink selected classes of a categorical raster by a number of cells, leaving other classes intact — the standard cleanup for classified land cover before polygonization, like ArcGIS Expand / Shrink.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input categorical (integer-valued) raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output raster path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "classes",
                    description: "Comma-separated class value(s) to expand or shrink.",
                    required: true,
                },
                ToolParamSpec {
                    name: "cells",
                    description: "Number of cells to expand or shrink by (>= 1).",
                    required: true,
                },
                ToolParamSpec {
                    name: "mode",
                    description: "'expand' (grow the selected classes; default) or 'shrink' (erode them).",
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
        let output = parse_optional_output(args, "output")?;
        let prm = parse_params(args)?;

        let raster = load_input_raster(input)?;
        if prm.band < 0 || prm.band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "band {} out of range (raster has {} band(s))",
                prm.band + 1,
                raster.bands
            )));
        }
        let rows = raster.rows;
        let cols = raster.cols;
        let nodata = raster.nodata;

        // Read the band into a flat grid.
        let mut grid = vec![0.0f64; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                grid[r * cols + c] = raster.get(prm.band, r as isize, c as isize);
            }
        }

        let is_selected = |v: f64| prm.classes.contains(&v);
        let is_valid = |v: f64| v != nodata && !v.is_nan();

        ctx.progress.info(&format!(
            "{} class(es) by {} cell(s) ({})",
            prm.classes.len(),
            prm.cells,
            match prm.mode {
                Mode::Expand => "expand",
                Mode::Shrink => "shrink",
            }
        ));

        let mut changed = 0usize;
        for _ in 0..prm.cells {
            let snapshot = grid.clone();
            let mut updates: Vec<(usize, f64)> = Vec::new();
            for r in 0..rows {
                for c in 0..cols {
                    let idx = r * cols + c;
                    let v = snapshot[idx];
                    if !is_valid(v) {
                        continue;
                    }
                    let sel = is_selected(v);
                    // Expand acts on non-selected boundary cells; shrink on
                    // selected boundary cells.
                    match prm.mode {
                        Mode::Expand if sel => continue,
                        Mode::Shrink if !sel => continue,
                        _ => {}
                    }
                    // Tally neighbour classes of the opposite kind.
                    let mut tally: BTreeMap<u64, (f64, usize)> = BTreeMap::new();
                    for (dr, dc) in NEIGH8 {
                        let nr = r as isize + dr;
                        let nc = c as isize + dc;
                        if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                            continue;
                        }
                        let nv = snapshot[nr as usize * cols + nc as usize];
                        if !is_valid(nv) {
                            continue;
                        }
                        let wanted = match prm.mode {
                            Mode::Expand => is_selected(nv),  // a selected neighbour
                            Mode::Shrink => !is_selected(nv), // a non-selected neighbour
                        };
                        if wanted {
                            let e = tally.entry(nv.to_bits()).or_insert((nv, 0));
                            e.1 += 1;
                        }
                    }
                    if tally.is_empty() {
                        continue;
                    }
                    // Pick the most frequent; ties -> smaller class value.
                    let best = tally
                        .values()
                        .max_by(|a, b| a.1.cmp(&b.1).then(b.0.total_cmp(&a.0)))
                        .map(|(v, _)| *v)
                        .unwrap();
                    updates.push((idx, best));
                }
            }
            if updates.is_empty() {
                break;
            }
            for (idx, v) in updates {
                grid[idx] = v;
                changed += 1;
            }
        }

        let out_raster = raster_like_with_data(&raster, grid, nodata, DataType::F32)?;
        let out_path = write_or_store(out_raster, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("cells_changed".to_string(), json!(changed));
        Ok(ToolRunResult { outputs })
    }
}

const NEIGH8: [(isize, isize); 8] = [
    (-1, -1),
    (-1, 0),
    (-1, 1),
    (0, -1),
    (0, 1),
    (1, -1),
    (1, 0),
    (1, 1),
];

fn write_or_store(raster: wbraster::Raster, output: Option<&str>) -> Result<String, ToolError> {
    crate::common::write_or_store_output(raster, output)
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    Expand,
    Shrink,
}

struct Params {
    classes: Vec<f64>,
    cells: usize,
    mode: Mode,
    band: isize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    // Accept a comma-separated string ("6" or "6,7") or a bare number (6).
    let classes_str = match args.get("classes") {
        Some(Value::String(s)) if !s.trim().is_empty() => s.trim().to_string(),
        Some(Value::Number(n)) => n.to_string(),
        _ => {
            return Err(ToolError::Validation(
                "required parameter 'classes' is missing".into(),
            ))
        }
    };
    let classes: Vec<f64> = classes_str
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| {
            s.parse::<f64>()
                .map_err(|_| ToolError::Validation(format!("class '{s}' is not a number")))
        })
        .collect::<Result<_, _>>()?;
    if classes.is_empty() {
        return Err(ToolError::Validation(
            "'classes' must name at least one class value".to_string(),
        ));
    }
    let cells = match args.get("cells") {
        Some(Value::Number(n)) => n.as_u64().unwrap_or(0) as usize,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'cells' must be a positive integer".into()))?,
        _ => {
            return Err(ToolError::Validation(
                "required parameter 'cells' is missing".to_string(),
            ))
        }
    };
    if cells < 1 {
        return Err(ToolError::Validation("'cells' must be >= 1".to_string()));
    }
    let mode = match args.get("mode").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("expand") => Mode::Expand,
        Some("shrink") => Mode::Shrink,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "'mode' must be 'expand' or 'shrink', got '{other}'"
            )))
        }
    };
    let band_1based = match args.get("band") {
        None | Some(Value::Null) => 1,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(1).max(1) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'band' must be an integer".into()))?
            .max(1),
        _ => 1,
    };
    Ok(Params {
        classes,
        cells,
        mode,
        band: (band_1based - 1) as isize,
    })
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

    fn raster_from(cols: usize, rows: usize, data: Vec<f64>) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -1.0,
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

    fn run(args: serde_json::Value) -> Raster {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ExpandShrinkTool.run(&args, &ctx()).unwrap();
        load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    fn count(r: &Raster, v: f64) -> usize {
        let mut n = 0;
        for row in 0..r.rows {
            for col in 0..r.cols {
                if r.get(0, row as isize, col as isize) == v {
                    n += 1;
                }
            }
        }
        n
    }

    /// A single class-2 cell in a field of 0s expands to a 3x3 block after one
    /// cell of expansion.
    #[test]
    fn expand_grows_a_single_cell() {
        let mut data = vec![0.0; 25]; // 5x5
        data[12] = 2.0; // center
        let input = raster_from(5, 5, data);
        let r = run(json!({ "input": input, "classes": "2", "cells": 1, "mode": "expand" }));
        // 8-connected expansion by 1 -> the center plus its 8 neighbours = 9.
        assert_eq!(count(&r, 2.0), 9, "expand by 1 should give a 3x3 block");
    }

    /// Shrinking a 3x3 block of class 2 by one cell leaves only the center.
    #[test]
    fn shrink_erodes_a_block() {
        let mut data = vec![0.0; 25];
        for (r, c) in [
            (1, 1),
            (1, 2),
            (1, 3),
            (2, 1),
            (2, 2),
            (2, 3),
            (3, 1),
            (3, 2),
            (3, 3),
        ] {
            data[r * 5 + c] = 2.0;
        }
        let input = raster_from(5, 5, data);
        let r = run(json!({ "input": input, "classes": "2", "cells": 1, "mode": "shrink" }));
        // Only the fully-interior center cell survives.
        assert_eq!(
            count(&r, 2.0),
            1,
            "shrink by 1 should leave only the center"
        );
    }

    /// Expansion never shrinks the class; shrink never grows it.
    #[test]
    fn monotone_area_change() {
        let mut data = vec![0.0; 100]; // 10x10
        for r in 3..7 {
            for c in 3..7 {
                data[r * 10 + c] = 5.0;
            }
        }
        let before = 16;
        let input = raster_from(10, 10, data);
        let exp =
            run(json!({ "input": input.clone(), "classes": "5", "cells": 2, "mode": "expand" }));
        let shr = run(json!({ "input": input, "classes": "5", "cells": 1, "mode": "shrink" }));
        assert!(count(&exp, 5.0) > before, "expand should grow the class");
        assert!(count(&shr, 5.0) < before, "shrink should reduce the class");
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ExpandShrinkTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.tif", "classes": "2" })).is_err()); // no cells
        assert!(bad(json!({ "input": "a.tif", "classes": "2", "cells": 0 })).is_err());
        assert!(bad(json!({ "input": "a.tif", "classes": "x", "cells": 1 })).is_err());
        assert!(
            bad(json!({ "input": "a.tif", "classes": "2", "cells": 1, "mode": "grow" })).is_err()
        );
        assert!(bad(json!({ "input": "a.tif", "classes": "2,3", "cells": 2 })).is_ok());
    }
}
