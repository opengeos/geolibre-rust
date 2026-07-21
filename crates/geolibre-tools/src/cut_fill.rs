//! GeoLibre tool: volumetric change between two DEM surfaces (cut / fill).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Cut Fill* (Spatial Analyst) and
//! *Surface Volume* (3D Analyst): measure how much material was moved between a
//! `before` and an `after` surface — earthworks, erosion / deposition, stockpile
//! measurement, lidar change detection — none of which the bundled raster suite
//! covers. A single surface can instead be compared against a reference `plane`
//! (surface-volume mode).
//!
//! The per-cell elevation change is `Δz = after − before` (two-surface mode) or
//! `Δz = surface − plane` (plane mode). Cells with `|Δz| ≤ tolerance` are
//! *unchanged*; `Δz > 0` is **fill** (material added, surface raised) and
//! `Δz < 0` is **cut** (material removed). Volumes accumulate `|Δz| × cell area`:
//!
//! - `fill_volume` — total material added,
//! - `cut_volume` — total material removed (reported positive),
//! - `net_volume` — `fill − cut` (positive = net gain).
//!
//! The primary output raster holds the signed `Δz` (nodata where either input
//! is nodata). Optionally, contiguous cut and fill patches are labelled into a
//! `region_output` raster (4-connectivity, cut and fill never share a region)
//! and a `csv_output` table lists each region's type, cell count, and volume —
//! matching the ArcGIS Cut Fill region output.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{
    band_to_vec, load_input_raster, parse_optional_output, raster_like_with_data,
    write_or_store_output, write_text_output,
};

/// No-data sentinel for the labelled region raster (region ids are >= 1).
const REGION_NODATA: f64 = -1.0;

pub struct CutFillTool;

impl Tool for CutFillTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "cut_fill",
            display_name: "Cut Fill",
            summary: "Volumetric change between two DEM surfaces (before/after) or a surface and a reference plane: signed elevation-change raster plus cut, fill, and net volumes, with optional per-region volume labelling.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "The 'before' (or the single surface for plane mode) raster file path.",
                    required: true,
                },
                ToolParamSpec {
                    name: "after",
                    description: "The 'after' raster (two-surface mode). Must share the grid of 'input'. Provide this or 'plane'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "plane",
                    description: "Reference plane elevation (surface-volume mode): change is measured as input minus this value. Provide this or 'after'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output raster path for the signed elevation change (Δz). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read from the input raster(s). Default 1.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Cells with |Δz| at or below this value (CRS z-units) count as unchanged. Default 0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "region_output",
                    description: "Optional raster path for contiguous cut/fill region labels (region ids >= 1; unchanged cells are nodata).",
                    required: false,
                },
                ToolParamSpec {
                    name: "csv_output",
                    description: "Optional CSV path for the per-region volume table (region_id, type, cell_count, volume).",
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
        let has_after = parse_optional_output(args, "after")?.is_some();
        let has_plane = parse_optional_f64(args, "plane")?.is_some();
        if has_after == has_plane {
            return Err(ToolError::Validation(
                "provide exactly one of 'after' (two-surface mode) or 'plane' (surface-volume mode)"
                    .to_string(),
            ));
        }
        parse_band(args)?;
        if let Some(t) = parse_optional_f64(args, "tolerance")? {
            if t < 0.0 || !t.is_finite() {
                return Err(ToolError::Validation(
                    "parameter 'tolerance' must be a non-negative number".to_string(),
                ));
            }
        }
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
        let after_path = parse_optional_output(args, "after")?;
        let plane = parse_optional_f64(args, "plane")?;
        let output = parse_optional_output(args, "output")?;
        let region_output = parse_optional_output(args, "region_output")?;
        let csv_output = parse_optional_output(args, "csv_output")?;
        let band = parse_band(args)?;
        let tolerance = parse_optional_f64(args, "tolerance")?.unwrap_or(0.0);
        if tolerance < 0.0 || !tolerance.is_finite() {
            return Err(ToolError::Validation(
                "parameter 'tolerance' must be a non-negative number".to_string(),
            ));
        }

        let before = load_input_raster(input)?;
        let band_idx = band as isize;
        if band >= before.bands {
            return Err(ToolError::Validation(format!(
                "band {} out of range (input has {} band(s))",
                band + 1,
                before.bands
            )));
        }
        let (rows, cols) = (before.rows, before.cols);
        let cell_area = before.cell_size_x * before.cell_size_y;
        let before_nd = before.nodata;
        let before_data = band_to_vec(&before, band_idx);

        // Build Δz per cell (NaN marks nodata / invalid).
        let mut diff = vec![f64::NAN; rows * cols];
        match (after_path, plane) {
            (Some(after_path), None) => {
                let after = load_input_raster(after_path)?;
                if after.rows != rows || after.cols != cols {
                    return Err(ToolError::Validation(format!(
                        "'after' grid {}x{} does not match 'input' {}x{}",
                        after.rows, after.cols, rows, cols
                    )));
                }
                if band >= after.bands {
                    return Err(ToolError::Validation(format!(
                        "band {} out of range ('after' has {} band(s))",
                        band + 1,
                        after.bands
                    )));
                }
                let after_nd = after.nodata;
                let after_data = band_to_vec(&after, band_idx);
                for i in 0..rows * cols {
                    let (b, a) = (before_data[i], after_data[i]);
                    if b != before_nd && a != after_nd && b.is_finite() && a.is_finite() {
                        diff[i] = a - b;
                    }
                }
            }
            (None, Some(plane)) => {
                for i in 0..rows * cols {
                    let b = before_data[i];
                    if b != before_nd && b.is_finite() {
                        diff[i] = b - plane;
                    }
                }
            }
            // validate() guarantees exactly one branch.
            _ => unreachable!("validate ensures exactly one of after/plane"),
        }

        // Classify and accumulate volumes.
        let (mut cut_volume, mut fill_volume) = (0.0, 0.0);
        let (mut cut_cells, mut fill_cells, mut unchanged_cells) = (0usize, 0usize, 0usize);
        let mut cls = vec![0i8; rows * cols]; // -1 cut, +1 fill, 0 unchanged/nodata
        for i in 0..rows * cols {
            let d = diff[i];
            if !d.is_finite() {
                continue;
            }
            if d > tolerance {
                fill_volume += d * cell_area;
                fill_cells += 1;
                cls[i] = 1;
            } else if d < -tolerance {
                cut_volume += -d * cell_area;
                cut_cells += 1;
                cls[i] = -1;
            } else {
                unchanged_cells += 1;
            }
        }
        let net_volume = fill_volume - cut_volume;

        ctx.progress.info(&format!(
            "cut {cut_volume:.3}, fill {fill_volume:.3}, net {net_volume:.3} (cell area {cell_area})"
        ));

        // Primary output: signed Δz raster (NaN -> nodata).
        let diff_data: Vec<f64> = diff
            .iter()
            .map(|&d| if d.is_finite() { d } else { before_nd })
            .collect();
        let diff_raster = raster_like_with_data(&before, diff_data, before_nd, DataType::F32)?;
        let out_path = write_or_store_output(diff_raster, output)?;

        // Optional region labelling (needed for region_output and/or csv_output).
        let mut region_count = 0usize;
        if region_output.is_some() || csv_output.is_some() {
            let (labels, regions) = label_regions(&cls, rows, cols, &diff, cell_area);
            region_count = regions.len();
            if let Some(path) = region_output {
                let region_data: Vec<f64> = labels
                    .iter()
                    .map(|&l| if l > 0 { l as f64 } else { REGION_NODATA })
                    .collect();
                let region_raster =
                    raster_like_with_data(&before, region_data, REGION_NODATA, DataType::F32)?;
                write_or_store_output(region_raster, Some(path))?;
            }
            if let Some(path) = csv_output {
                let mut csv = String::from("region_id,type,cell_count,volume\n");
                for r in &regions {
                    csv.push_str(&format!(
                        "{},{},{},{:.6}\n",
                        r.id,
                        if r.is_fill { "fill" } else { "cut" },
                        r.cells,
                        r.volume
                    ));
                }
                write_text_output(&csv, path)?;
            }
        }

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("cut_volume".to_string(), json!(cut_volume));
        outputs.insert("fill_volume".to_string(), json!(fill_volume));
        outputs.insert("net_volume".to_string(), json!(net_volume));
        outputs.insert("cut_cells".to_string(), json!(cut_cells));
        outputs.insert("fill_cells".to_string(), json!(fill_cells));
        outputs.insert("unchanged_cells".to_string(), json!(unchanged_cells));
        outputs.insert("cell_area".to_string(), json!(cell_area));
        if region_output.is_some() || csv_output.is_some() {
            outputs.insert("region_count".to_string(), json!(region_count));
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── Region labelling ──────────────────────────────────────────────────────────

struct Region {
    id: u32,
    is_fill: bool,
    cells: usize,
    volume: f64,
}

/// Labels contiguous (4-connected) cut and fill patches; cut and fill never
/// share a region. Returns per-cell labels (0 = none) and per-region summaries.
fn label_regions(
    cls: &[i8],
    rows: usize,
    cols: usize,
    diff: &[f64],
    cell_area: f64,
) -> (Vec<u32>, Vec<Region>) {
    let mut labels = vec![0u32; rows * cols];
    let mut regions: Vec<Region> = Vec::new();
    let mut stack: Vec<usize> = Vec::new();
    for start in 0..rows * cols {
        if cls[start] == 0 || labels[start] != 0 {
            continue;
        }
        let sign = cls[start];
        let id = regions.len() as u32 + 1;
        let mut cells = 0usize;
        let mut volume = 0.0;
        stack.push(start);
        labels[start] = id;
        while let Some(i) = stack.pop() {
            cells += 1;
            volume += diff[i].abs() * cell_area;
            let (r, c) = (i / cols, i % cols);
            let push = |nr: isize, nc: isize, stack: &mut Vec<usize>, labels: &mut Vec<u32>| {
                if nr < 0 || nc < 0 || nr as usize >= rows || nc as usize >= cols {
                    return;
                }
                let j = nr as usize * cols + nc as usize;
                if cls[j] == sign && labels[j] == 0 {
                    labels[j] = id;
                    stack.push(j);
                }
            };
            push(r as isize - 1, c as isize, &mut stack, &mut labels);
            push(r as isize + 1, c as isize, &mut stack, &mut labels);
            push(r as isize, c as isize - 1, &mut stack, &mut labels);
            push(r as isize, c as isize + 1, &mut stack, &mut labels);
        }
        regions.push(Region {
            id,
            is_fill: sign > 0,
            cells,
            volume,
        });
    }
    (labels, regions)
}

// ── Parameters ────────────────────────────────────────────────────────────────

fn parse_band(args: &ToolArgs) -> Result<usize, ToolError> {
    match parse_optional_f64(args, "band")? {
        None => Ok(0),
        Some(v) if v.fract() == 0.0 && v >= 1.0 && v.is_finite() => Ok(v as usize - 1),
        Some(_) => Err(ToolError::Validation(
            "parameter 'band' must be a positive integer".to_string(),
        )),
    }
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
    use wbraster::{memory_store, CrsInfo, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a rows x cols raster (top-down, row-major) with 1 m cells and the
    /// given nodata value.
    fn raster(rows: usize, cols: usize, data: Vec<f64>, nodata: f64) -> Raster {
        let cfg = RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata,
            data_type: DataType::F64,
            crs: CrsInfo::from_epsg(32610),
            metadata: Vec::new(),
        };
        Raster::from_data(cfg, data).unwrap()
    }

    fn path(r: Raster) -> String {
        let id = memory_store::put_raster(r);
        memory_store::make_raster_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        CutFillTool.run(&args, &ctx()).unwrap()
    }

    #[test]
    fn two_surface_cut_and_fill_volumes() {
        // 2x2, 1 m cells. before all 10. after: two cells +2 (fill), one -3 (cut),
        // one unchanged. fill_volume = 2*2*1 = 4; cut_volume = 3*1 = 3; net = 1.
        let before = path(raster(2, 2, vec![10.0, 10.0, 10.0, 10.0], -9999.0));
        let after = path(raster(2, 2, vec![12.0, 12.0, 7.0, 10.0], -9999.0));
        let out = run(json!({ "input": before, "after": after }));
        assert!((out.outputs["fill_volume"].as_f64().unwrap() - 4.0).abs() < 1e-9);
        assert!((out.outputs["cut_volume"].as_f64().unwrap() - 3.0).abs() < 1e-9);
        assert!((out.outputs["net_volume"].as_f64().unwrap() - 1.0).abs() < 1e-9);
        assert_eq!(out.outputs["fill_cells"], json!(2));
        assert_eq!(out.outputs["cut_cells"], json!(1));
        assert_eq!(out.outputs["unchanged_cells"], json!(1));
    }

    #[test]
    fn nodata_cells_are_excluded() {
        let before = path(raster(2, 2, vec![10.0, 10.0, -9999.0, 10.0], -9999.0));
        let after = path(raster(2, 2, vec![12.0, 10.0, 5.0, -9999.0], -9999.0));
        let out = run(json!({ "input": before, "after": after }));
        // Only the top row is valid in both: (10->12 fill 2), (10->10 unchanged).
        assert!((out.outputs["fill_volume"].as_f64().unwrap() - 2.0).abs() < 1e-9);
        assert_eq!(out.outputs["cut_cells"], json!(0));
        assert_eq!(out.outputs["fill_cells"], json!(1));
    }

    #[test]
    fn surface_volume_against_a_plane() {
        // Surface values vs plane 5: cells 8,6 are above (+3,+1 fill), 2 below
        // (-3 cut), 5 equal (unchanged). fill = 4, cut = 3, net = 1.
        let surf = path(raster(2, 2, vec![8.0, 6.0, 2.0, 5.0], -9999.0));
        let out = run(json!({ "input": surf, "plane": 5.0 }));
        assert!((out.outputs["fill_volume"].as_f64().unwrap() - 4.0).abs() < 1e-9);
        assert!((out.outputs["cut_volume"].as_f64().unwrap() - 3.0).abs() < 1e-9);
        assert_eq!(out.outputs["unchanged_cells"], json!(1));
    }

    #[test]
    fn tolerance_treats_small_changes_as_unchanged() {
        let before = path(raster(1, 3, vec![10.0, 10.0, 10.0], -9999.0));
        let after = path(raster(1, 3, vec![10.4, 12.0, 8.0], -9999.0));
        // tolerance 0.5: the +0.4 cell is unchanged; +2 fill, -2 cut remain.
        let out = run(json!({ "input": before, "after": after, "tolerance": 0.5 }));
        assert_eq!(out.outputs["unchanged_cells"], json!(1));
        assert_eq!(out.outputs["fill_cells"], json!(1));
        assert_eq!(out.outputs["cut_cells"], json!(1));
    }

    #[test]
    fn labels_contiguous_regions() {
        // A row: fill, fill, unchanged, cut. -> region 1 (fill, 2 cells),
        // region 2 (cut, 1 cell).
        let before = path(raster(1, 4, vec![10.0, 10.0, 10.0, 10.0], -9999.0));
        let after = path(raster(1, 4, vec![12.0, 13.0, 10.0, 7.0], -9999.0));
        let csv = format!(
            "{}/cf_regions_{}.csv",
            std::env::temp_dir().display(),
            std::process::id()
        );
        let out = run(json!({ "input": before, "after": after, "csv_output": csv }));
        assert_eq!(out.outputs["region_count"], json!(2));
        let text = std::fs::read_to_string(&csv).unwrap();
        let _ = std::fs::remove_file(&csv);
        // fill region volume = (2+3)*1 = 5; cut region volume = 3.
        assert!(text.contains(",fill,2,5.000000"), "csv was:\n{text}");
        assert!(text.contains(",cut,1,3.000000"), "csv was:\n{text}");
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = CutFillTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input");
        assert!(
            bad(json!({ "input": "a.tif" })).is_err(),
            "need after or plane"
        );
        assert!(
            bad(json!({ "input": "a.tif", "after": "b.tif", "plane": 5.0 })).is_err(),
            "not both"
        );
        assert!(bad(json!({ "input": "a.tif", "plane": 5.0 })).is_ok());
        assert!(bad(json!({ "input": "a.tif", "after": "b.tif", "tolerance": -1.0 })).is_err());
        assert!(bad(json!({ "input": "a.tif", "after": "b.tif", "band": 0 })).is_err());
    }
}
