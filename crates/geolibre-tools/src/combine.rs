//! GeoLibre tool: assign a unique zone ID to each unique combination of values
//! across two or more co-registered categorical rasters.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Combine* (Spatial Analyst). The
//! bundled whitebox raster suite only offers `cross_tabulation` (a two-raster
//! contingency table) — it cannot build the *unique-condition units* raster
//! used in suitability modelling and hydrologic response units (HRUs), where
//! every distinct tuple of input class values receives a single new integer id.
//!
//! Given `inputs` (a comma-separated list of ≥2 aligned integer/categorical
//! rasters), the tool walks every cell, forms the tuple of input values, and
//! maps each *distinct* tuple to a dense id starting at 1 (matching ArcGIS,
//! whose `VALUE` column is 1-based). Any cell that is nodata in *any* input is
//! nodata in the output. The primary output is a signed-integer raster of these
//! ids; an optional `csv_output` writes the value-attribute table (VAT):
//! `value, count, <input_1>, <input_2>, …` — one row per unique combination.
//!
//! The 3-file pattern (this file + `lib.rs` registration + README row) matches
//! the other GeoLibre tools; raster I/O reuses `crate::common`.

use std::collections::{BTreeMap, HashMap};

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

/// No-data sentinel for the integer id output (ids are >= 1).
const ID_NODATA: f64 = -1.0;

/// Assign a unique integer id to every unique combination of values across a set
/// of co-registered categorical rasters (ArcGIS Combine).
pub struct CombineTool;

impl Tool for CombineTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "combine",
            display_name: "Combine",
            summary: "Assign a unique zone id to each unique combination of values across two or more co-registered categorical rasters, plus a value-attribute table (id, count, one column per input).",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "inputs",
                    description: "Comma-separated list of >=2 co-registered categorical/integer rasters to combine.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output raster path for the unique-condition id grid. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "csv_output",
                    description: "Optional CSV path for the value-attribute table (value, count, one column per input raster).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read from every input raster. Default 1.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        let paths = parse_inputs(args)?;
        if paths.len() < 2 {
            return Err(ToolError::Validation(
                "parameter 'inputs' must list at least 2 comma-separated rasters".to_string(),
            ));
        }
        parse_band(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let paths = parse_inputs(args)?;
        if paths.len() < 2 {
            return Err(ToolError::Validation(
                "parameter 'inputs' must list at least 2 comma-separated rasters".to_string(),
            ));
        }
        let output = parse_optional_output(args, "output")?;
        let csv_output = parse_optional_output(args, "csv_output")?;
        let band = parse_band(args)?;
        let band_idx = band as isize;

        // Load and align all inputs against the first raster's grid.
        let rasters: Vec<_> = paths
            .iter()
            .map(|p| load_input_raster(p))
            .collect::<Result<Vec<_>, _>>()?;
        let (rows, cols) = (rasters[0].rows, rasters[0].cols);
        for (i, r) in rasters.iter().enumerate() {
            if r.rows != rows || r.cols != cols {
                return Err(ToolError::Validation(format!(
                    "input {} is {}x{}, expected {rows}x{cols} — all inputs must be co-registered",
                    i, r.rows, r.cols
                )));
            }
            if band >= r.bands {
                return Err(ToolError::Validation(format!(
                    "band {} out of range (input {} has {} band(s))",
                    band + 1,
                    i,
                    r.bands
                )));
            }
        }

        let nodata: Vec<f64> = rasters.iter().map(|r| r.nodata).collect();
        let bands: Vec<Vec<f64>> = rasters.iter().map(|r| band_to_vec(r, band_idx)).collect();

        ctx.progress
            .info(&format!("combining {} categorical rasters", rasters.len()));

        // Map each unique value-tuple to a dense 1-based id; count cells per id
        // and keep insertion order so the VAT is stable.
        let mut lookup: HashMap<Vec<i64>, u32> = HashMap::new();
        let mut combos: Vec<Combo> = Vec::new();
        let mut ids = vec![ID_NODATA; rows * cols];
        let mut valid_cells = 0usize;

        for i in 0..rows * cols {
            let mut tuple: Vec<i64> = Vec::with_capacity(rasters.len());
            let mut ok = true;
            for (b, &nd) in bands.iter().zip(&nodata) {
                let v = b[i];
                if v == nd || !v.is_finite() {
                    ok = false;
                    break;
                }
                tuple.push(v.round() as i64);
            }
            if !ok {
                continue;
            }
            valid_cells += 1;
            let id = match lookup.get(&tuple) {
                Some(&id) => {
                    combos[(id - 1) as usize].count += 1;
                    id
                }
                None => {
                    let id = combos.len() as u32 + 1;
                    lookup.insert(tuple.clone(), id);
                    combos.push(Combo {
                        id,
                        values: tuple,
                        count: 1,
                    });
                    id
                }
            };
            ids[i] = id as f64;
        }

        if combos.is_empty() {
            return Err(ToolError::Execution(
                "no cell is valid across all inputs — nothing to combine".to_string(),
            ));
        }

        ctx.progress.info(&format!(
            "{} unique combinations over {valid_cells} valid cells",
            combos.len()
        ));

        // Primary output: integer id raster (nodata where any input is nodata).
        let id_raster = raster_like_with_data(&rasters[0], ids, ID_NODATA, DataType::I32)?;
        let out_path = write_or_store_output(id_raster, output)?;

        // Optional value-attribute table.
        if let Some(path) = csv_output {
            let labels = column_labels(&paths);
            let mut csv = String::from("value,count");
            for label in &labels {
                csv.push(',');
                csv.push_str(label);
            }
            csv.push('\n');
            for c in &combos {
                csv.push_str(&format!("{},{}", c.id, c.count));
                for v in &c.values {
                    csv.push_str(&format!(",{v}"));
                }
                csv.push('\n');
            }
            write_text_output(&csv, path)?;
        }

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("combinations".to_string(), json!(combos.len()));
        outputs.insert("valid_cells".to_string(), json!(valid_cells));
        outputs.insert("input_count".to_string(), json!(rasters.len()));
        Ok(ToolRunResult { outputs })
    }
}

/// One unique combination of input values.
struct Combo {
    id: u32,
    values: Vec<i64>,
    count: u64,
}

/// Parses the comma-separated `inputs` list into trimmed, non-empty paths.
fn parse_inputs(args: &ToolArgs) -> Result<Vec<String>, ToolError> {
    let list = args
        .get("inputs")
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| {
            ToolError::Validation("missing required string parameter 'inputs'".to_string())
        })?;
    Ok(list
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_string)
        .collect())
}

/// 1-based `band` param -> 0-based index (default band 1). Rejects non-positive
/// and non-integer values.
fn parse_band(args: &ToolArgs) -> Result<usize, ToolError> {
    match args.get("band") {
        None | Some(Value::Null) => Ok(0),
        Some(Value::Number(n)) => band_from_f64(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(0),
        Some(Value::String(s)) => match s.trim().parse::<f64>() {
            Ok(v) => band_from_f64(Some(v)),
            Err(_) => Err(ToolError::Validation(
                "parameter 'band' must be a positive integer".to_string(),
            )),
        },
        Some(_) => Err(ToolError::Validation(
            "parameter 'band' must be a positive integer".to_string(),
        )),
    }
}

fn band_from_f64(v: Option<f64>) -> Result<usize, ToolError> {
    match v {
        Some(v) if v.is_finite() && v.fract() == 0.0 && v >= 1.0 => Ok(v as usize - 1),
        _ => Err(ToolError::Validation(
            "parameter 'band' must be a positive integer".to_string(),
        )),
    }
}

/// Derives a stable, unique column header for each input from its file basename
/// (directory and extension stripped). Falls back to `in{i}` for memory handles
/// or empty stems, and de-duplicates collisions with a numeric suffix.
fn column_labels(paths: &[String]) -> Vec<String> {
    let mut labels = Vec::with_capacity(paths.len());
    let mut seen: HashMap<String, u32> = HashMap::new();
    for (i, p) in paths.iter().enumerate() {
        let stem = std::path::Path::new(p)
            .file_stem()
            .and_then(|s| s.to_str())
            .map(sanitize)
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| format!("in{}", i + 1));
        let label = match seen.get_mut(&stem) {
            Some(n) => {
                *n += 1;
                format!("{stem}_{n}")
            }
            None => {
                seen.insert(stem.clone(), 0);
                stem
            }
        };
        labels.push(label);
    }
    labels
}

/// Keeps CSV-safe header characters, mapping anything else to `_`.
fn sanitize(s: &str) -> String {
    s.chars()
        .map(|c| {
            if c.is_alphanumeric() || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect()
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

    fn read_raster(p: &str) -> Raster {
        let id = memory_store::raster_path_to_id(p).unwrap();
        (*memory_store::get_raster_arc_by_id(id).unwrap()).clone()
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        CombineTool.run(&args, &ctx()).unwrap()
    }

    #[test]
    fn unique_combinations_get_dense_ids() {
        // a: [1,1,2,2]  b: [1,2,1,2]  -> tuples (1,1)(1,2)(2,1)(2,2) all distinct
        let a = path(raster(2, 2, vec![1.0, 1.0, 2.0, 2.0], -9999.0));
        let b = path(raster(2, 2, vec![1.0, 2.0, 1.0, 2.0], -9999.0));
        let out = run(json!({ "inputs": format!("{a},{b}") }));
        assert_eq!(out.outputs["combinations"], json!(4));
        assert_eq!(out.outputs["valid_cells"], json!(4));
        // Output ids are the four distinct labels 1..=4 (dense, no gaps).
        let r = read_raster(out.outputs["output"].as_str().unwrap());
        let mut vals: Vec<i64> = (0..4)
            .map(|i| r.get(0, (i / 2) as isize, (i % 2) as isize).round() as i64)
            .collect();
        vals.sort();
        assert_eq!(vals, vec![1, 2, 3, 4]);
    }

    #[test]
    fn identical_tuples_share_one_id_and_count() {
        // Two cells share tuple (5,7); one cell has (5,8). -> 2 combos.
        let a = path(raster(1, 3, vec![5.0, 5.0, 5.0], -1.0));
        let b = path(raster(1, 3, vec![7.0, 7.0, 8.0], -1.0));
        let csv = format!(
            "{}/combine_vat_{}.csv",
            std::env::temp_dir().display(),
            std::process::id()
        );
        let out = run(json!({ "inputs": format!("{a},{b}"), "csv_output": csv }));
        assert_eq!(out.outputs["combinations"], json!(2));
        let r = read_raster(out.outputs["output"].as_str().unwrap());
        // First two cells share an id; third differs.
        assert_eq!(r.get(0, 0, 0).round() as i64, r.get(0, 0, 1).round() as i64);
        assert_ne!(r.get(0, 0, 0).round() as i64, r.get(0, 0, 2).round() as i64);
        let text = std::fs::read_to_string(&csv).unwrap();
        let _ = std::fs::remove_file(&csv);
        // VAT: id 1 has count 2 with values 5,7; id 2 count 1 with 5,8.
        assert!(text.contains("value,count,"), "header was:\n{text}");
        assert!(text.contains("1,2,5,7"), "vat was:\n{text}");
        assert!(text.contains("2,1,5,8"), "vat was:\n{text}");
    }

    #[test]
    fn counts_conserve_valid_cells() {
        // Cell counts across the VAT must sum to valid_cells.
        let a = path(raster(2, 3, vec![1.0, 1.0, 2.0, 2.0, 3.0, 3.0], -9999.0));
        let b = path(raster(2, 3, vec![9.0, 9.0, 9.0, 8.0, 8.0, 8.0], -9999.0));
        let out = run(json!({ "inputs": format!("{a},{b}") }));
        assert_eq!(out.outputs["valid_cells"], json!(6));
        // (1,9)(2,9)(2,8)(3,8) -> 4 combos.
        assert_eq!(out.outputs["combinations"], json!(4));
    }

    #[test]
    fn nodata_in_any_input_is_nodata_out() {
        let a = path(raster(1, 3, vec![1.0, -9999.0, 3.0], -9999.0));
        let b = path(raster(1, 3, vec![4.0, 5.0, -9999.0], -9999.0));
        let out = run(json!({ "inputs": format!("{a},{b}") }));
        // Only cell 0 is valid in both.
        assert_eq!(out.outputs["valid_cells"], json!(1));
        assert_eq!(out.outputs["combinations"], json!(1));
        let r = read_raster(out.outputs["output"].as_str().unwrap());
        assert_eq!(r.get(0, 0, 0).round() as i64, 1);
        assert_eq!(r.get(0, 0, 1), r.nodata);
        assert_eq!(r.get(0, 0, 2), r.nodata);
    }

    #[test]
    fn combines_three_rasters() {
        let a = path(raster(1, 2, vec![1.0, 1.0], -1.0));
        let b = path(raster(1, 2, vec![2.0, 2.0], -1.0));
        let c = path(raster(1, 2, vec![3.0, 4.0], -1.0));
        let out = run(json!({ "inputs": format!("{a},{b},{c}") }));
        assert_eq!(out.outputs["input_count"], json!(3));
        // (1,2,3) and (1,2,4) -> 2 combos.
        assert_eq!(out.outputs["combinations"], json!(2));
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = CombineTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing inputs");
        assert!(bad(json!({ "inputs": "" })).is_err(), "empty inputs");
        assert!(
            bad(json!({ "inputs": "only_one.tif" })).is_err(),
            "need >=2"
        );
        assert!(
            bad(json!({ "inputs": "a.tif,b.tif", "band": 0 })).is_err(),
            "band must be positive"
        );
        assert!(bad(json!({ "inputs": "a.tif,b.tif" })).is_ok());
    }

    #[test]
    fn mismatched_grids_are_rejected() {
        let a = path(raster(2, 2, vec![1.0, 2.0, 3.0, 4.0], -1.0));
        let b = path(raster(1, 2, vec![1.0, 2.0], -1.0));
        let args: ToolArgs =
            serde_json::from_value(json!({ "inputs": format!("{a},{b}") })).unwrap();
        assert!(CombineTool.run(&args, &ctx()).is_err());
    }
}
