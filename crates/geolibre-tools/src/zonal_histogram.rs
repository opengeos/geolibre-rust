//! GeoLibre tool: per-zone raster value distribution table (zonal histogram).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Zonal Histogram* (Spatial Analyst):
//! tabulate the full **distribution** of a value raster's cells within each
//! zone, as a cross-tab table of zones x classes/bins. This is distinct from
//! the bundled `zonal_statistics`, which only reports summary moments
//! (mean/min/max/stddev) per zone and throws away the shape of the
//! distribution, and from `cross_tabulation`, which needs *two* categorical
//! rasters and reports zone x zone rather than zone x a value histogram.
//! "What is the land-cover composition of each district" or "how are DEM
//! elevations distributed within each watershed" has no direct answer
//! elsewhere in the suite.
//!
//! `zones` and `value` must share the same grid (rows/cols/cell size). In
//! `classes` mode (the default, for integer/categorical value rasters) each
//! distinct rounded value is its own class. In `bins` mode (for continuous
//! value rasters) the valid value range is split into `bins` equal-width
//! buckets computed from the data in a single pass.
//!
//! The primary `output` is a wide CSV: one row per zone, one column per
//! class/bin, holding cell counts (or, with `percent`, each class's percentage
//! of that zone's valid cells). An optional `long_output` CSV lists the same
//! data as `zone,class,count,percent` rows, which is easier to plot or join.
//!
//! Scope for v1: `zones` must be a raster (the zone-raster path required by
//! the issue); rasterizing a polygon zone layer on the fly is not yet
//! implemented.

use std::collections::{BTreeMap, BTreeSet};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};

use crate::common::{band_to_vec, load_input_raster, write_text_output};

pub struct ZonalHistogramTool;

impl Tool for ZonalHistogramTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "zonal_histogram",
            display_name: "Zonal Histogram",
            summary: "Tabulate the distribution of a value raster's cells within each zone: class counts (categorical) or binned frequencies (continuous), as a zones x classes cross-tab table.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "zones",
                    description: "Zone raster file path (integer/categorical zone ids).",
                    required: true,
                },
                ToolParamSpec {
                    name: "value",
                    description: "Value raster file path; must share the grid of 'zones'.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output CSV path: one row per zone, one column per class/bin.",
                    required: true,
                },
                ToolParamSpec {
                    name: "mode",
                    description: "'classes' (default) treats each rounded value as its own class, for integer/categorical value rasters. 'bins' buckets continuous values into equal-width bins.",
                    required: false,
                },
                ToolParamSpec {
                    name: "bins",
                    description: "Number of equal-width bins for 'bins' mode. Default 10.",
                    required: false,
                },
                ToolParamSpec {
                    name: "percent",
                    description: "Emit each class's percentage of the zone's valid cells instead of raw cell counts. Default false.",
                    required: false,
                },
                ToolParamSpec {
                    name: "zone_band",
                    description: "1-based band to read from 'zones'. Default 1.",
                    required: false,
                },
                ToolParamSpec {
                    name: "value_band",
                    description: "1-based band to read from 'value'. Default 1.",
                    required: false,
                },
                ToolParamSpec {
                    name: "long_output",
                    description: "Optional CSV path for the long-format (zone,class,count,percent) variant.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["zones", "value", "output"] {
            if args
                .get(key)
                .and_then(Value::as_str)
                .map(str::trim)
                .unwrap_or("")
                .is_empty()
            {
                return Err(ToolError::Validation(format!(
                    "missing required string parameter '{key}'"
                )));
            }
        }
        parse_mode(args)?;
        if let Some(b) = parse_optional_u64(args, "bins")? {
            if b < 2 {
                return Err(ToolError::Validation(
                    "parameter 'bins' must be at least 2".to_string(),
                ));
            }
        }
        parse_band(args, "zone_band")?;
        parse_band(args, "value_band")?;
        parse_optional_bool(args, "percent")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let zones_path = require_str(args, "zones")?;
        let value_path = require_str(args, "value")?;
        let output = require_str(args, "output")?;
        let long_output = parse_optional_str(args, "long_output")?;
        let mode = parse_mode(args)?;
        let bins = parse_optional_u64(args, "bins")?.unwrap_or(10) as usize;
        let percent = parse_optional_bool(args, "percent")?.unwrap_or(false);
        let zone_band = parse_band(args, "zone_band")?;
        let value_band = parse_band(args, "value_band")?;

        let zones = load_input_raster(zones_path)?;
        let value = load_input_raster(value_path)?;
        if zones.rows != value.rows || zones.cols != value.cols {
            return Err(ToolError::Validation(format!(
                "'zones' grid {}x{} does not match 'value' grid {}x{}",
                zones.rows, zones.cols, value.rows, value.cols
            )));
        }
        if zone_band >= zones.bands {
            return Err(ToolError::Validation(format!(
                "zone_band {} out of range ('zones' has {} band(s))",
                zone_band + 1,
                zones.bands
            )));
        }
        if value_band >= value.bands {
            return Err(ToolError::Validation(format!(
                "value_band {} out of range ('value' has {} band(s))",
                value_band + 1,
                value.bands
            )));
        }

        let (rows, cols) = (zones.rows, zones.cols);
        let zone_nd = zones.nodata;
        let value_nd = value.nodata;
        let zone_data = band_to_vec(&zones, zone_band as isize);
        let value_data = band_to_vec(&value, value_band as isize);

        // A cell is counted when both zone and value are valid (not the
        // raster's nodata sentinel, and finite).
        let is_valid =
            |z: f64, v: f64| z != zone_nd && v != value_nd && z.is_finite() && v.is_finite();

        // For 'bins' mode, a first pass finds the value range actually being
        // tabulated (cells whose zone is also valid), so bin edges reflect the
        // data that lands in the table.
        let (bin_min, bin_max) = if mode == Mode::Bins {
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for i in 0..rows * cols {
                let (z, v) = (zone_data[i], value_data[i]);
                if is_valid(z, v) {
                    lo = lo.min(v);
                    hi = hi.max(v);
                }
            }
            if !lo.is_finite() || !hi.is_finite() {
                return Err(ToolError::Execution(
                    "no valid (non-nodata) overlapping cells between 'zones' and 'value'"
                        .to_string(),
                ));
            }
            (lo, hi)
        } else {
            (0.0, 0.0)
        };
        let bin_width = if bin_max > bin_min {
            (bin_max - bin_min) / bins as f64
        } else {
            0.0
        };

        ctx.progress.info(&format!(
            "tabulating {rows}x{cols} grid in {} mode",
            mode.as_str()
        ));

        // Single pass: accumulate zone x class counts and per-zone totals.
        let mut counts: BTreeMap<(i64, i64), u64> = BTreeMap::new();
        let mut zone_totals: BTreeMap<i64, u64> = BTreeMap::new();
        for i in 0..rows * cols {
            let (z, v) = (zone_data[i], value_data[i]);
            if !is_valid(z, v) {
                continue;
            }
            let zone_key = z.round() as i64;
            let class_key = match mode {
                Mode::Classes => v.round() as i64,
                Mode::Bins => {
                    if bin_width > 0.0 {
                        let idx = ((v - bin_min) / bin_width).floor() as i64;
                        idx.clamp(0, bins as i64 - 1)
                    } else {
                        0
                    }
                }
            };
            *counts.entry((zone_key, class_key)).or_insert(0) += 1;
            *zone_totals.entry(zone_key).or_insert(0) += 1;
        }

        if zone_totals.is_empty() {
            return Err(ToolError::Execution(
                "no valid (non-nodata) overlapping cells between 'zones' and 'value'".to_string(),
            ));
        }

        let zone_keys: Vec<i64> = zone_totals.keys().copied().collect();
        let class_keys: BTreeSet<i64> = counts.keys().map(|&(_, c)| c).collect();
        let class_keys: Vec<i64> = class_keys.into_iter().collect();
        let class_labels: Vec<String> = class_keys
            .iter()
            .map(|&c| class_label(mode, c, bin_min, bin_width, bins))
            .collect();

        // Wide CSV: one row per zone, one column per class/bin.
        let mut wide = String::from("zone");
        for label in &class_labels {
            wide.push(',');
            wide.push_str(label);
        }
        wide.push('\n');
        for &zone in &zone_keys {
            let total = *zone_totals.get(&zone).unwrap_or(&0) as f64;
            wide.push_str(&zone.to_string());
            for &class in &class_keys {
                let count = counts.get(&(zone, class)).copied().unwrap_or(0);
                wide.push(',');
                if percent {
                    let pct = if total > 0.0 {
                        count as f64 / total * 100.0
                    } else {
                        0.0
                    };
                    wide.push_str(&format!("{pct:.4}"));
                } else {
                    wide.push_str(&count.to_string());
                }
            }
            wide.push('\n');
        }
        write_text_output(&wide, output)?;

        if let Some(long_path) = long_output {
            let mut long = String::from("zone,class,count,percent\n");
            for &zone in &zone_keys {
                let total = *zone_totals.get(&zone).unwrap_or(&0) as f64;
                for (idx, &class) in class_keys.iter().enumerate() {
                    let count = counts.get(&(zone, class)).copied().unwrap_or(0);
                    if count == 0 {
                        continue;
                    }
                    let pct = if total > 0.0 {
                        count as f64 / total * 100.0
                    } else {
                        0.0
                    };
                    long.push_str(&format!("{zone},{},{count},{pct:.4}\n", class_labels[idx]));
                }
            }
            write_text_output(&long, long_path)?;
        }

        let total_valid_cells: u64 = zone_totals.values().sum();
        ctx.progress.info(&format!(
            "{} zone(s) x {} class(es), {total_valid_cells} valid cell(s)",
            zone_keys.len(),
            class_keys.len()
        ));

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(output));
        if let Some(long_path) = long_output {
            outputs.insert("long_output".to_string(), json!(long_path));
        }
        outputs.insert("zone_count".to_string(), json!(zone_keys.len()));
        outputs.insert("class_count".to_string(), json!(class_keys.len()));
        outputs.insert("total_valid_cells".to_string(), json!(total_valid_cells));
        outputs.insert("mode".to_string(), json!(mode.as_str()));
        Ok(ToolRunResult { outputs })
    }
}

// ── Mode ──────────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    Classes,
    Bins,
}

impl Mode {
    fn as_str(self) -> &'static str {
        match self {
            Self::Classes => "classes",
            Self::Bins => "bins",
        }
    }
}

/// Formats a class key as a CSV-safe column label: the class value itself in
/// `classes` mode, or a `lo_hi` half-open range in `bins` mode.
fn class_label(mode: Mode, class: i64, bin_min: f64, bin_width: f64, bins: usize) -> String {
    match mode {
        Mode::Classes => class.to_string(),
        Mode::Bins => {
            if bin_width > 0.0 {
                let lo = bin_min + class as f64 * bin_width;
                let hi = if class as usize + 1 >= bins {
                    bin_min + bins as f64 * bin_width
                } else {
                    bin_min + (class + 1) as f64 * bin_width
                };
                format!("{lo:.4}_{hi:.4}")
            } else {
                // Degenerate range (single distinct value): one bin covering it.
                format!("{bin_min:.4}_{bin_min:.4}")
            }
        }
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

fn parse_optional_str<'a>(args: &'a ToolArgs, key: &str) -> Result<Option<&'a str>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => Ok(Some(s.as_str())),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a string when provided"
        ))),
    }
}

fn parse_mode(args: &ToolArgs) -> Result<Mode, ToolError> {
    match parse_optional_str(args, "mode")? {
        None => Ok(Mode::Classes),
        Some(s) => match s.to_ascii_lowercase().as_str() {
            "classes" => Ok(Mode::Classes),
            "bins" => Ok(Mode::Bins),
            other => Err(ToolError::Validation(format!(
                "parameter 'mode' must be 'classes' or 'bins', got '{other}'"
            ))),
        },
    }
}

fn parse_band(args: &ToolArgs, key: &str) -> Result<usize, ToolError> {
    match parse_optional_f64(args, key)? {
        None => Ok(0),
        Some(v) if v.fract() == 0.0 && v >= 1.0 && v.is_finite() => Ok(v as usize - 1),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a positive integer"
        ))),
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

fn parse_optional_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64().or_else(|| n.as_f64().map(|f| f as u64))),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be an integer"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be an integer"
        ))),
    }
}

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

    /// Builds a rows x cols raster (row-major) with 1 m cells.
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

    fn tmp_csv(tag: &str) -> String {
        format!(
            "{}/zh_{tag}_{}.csv",
            std::env::temp_dir().display(),
            std::process::id()
        )
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        ZonalHistogramTool.run(&args, &ctx()).unwrap()
    }

    #[test]
    fn two_zone_two_class_exact_counts() {
        // 2x2 zone raster: zone 1 top row, zone 2 bottom row.
        // value raster: top row both class 10; bottom row one class 10, one 20.
        let zones = path(raster(2, 2, vec![1.0, 1.0, 2.0, 2.0], -9999.0));
        let value = path(raster(2, 2, vec![10.0, 10.0, 10.0, 20.0], -9999.0));
        let out_csv = tmp_csv("classes");
        let out = run(json!({ "zones": zones, "value": value, "output": out_csv }));
        assert_eq!(out.outputs["zone_count"], json!(2));
        assert_eq!(out.outputs["class_count"], json!(2));
        assert_eq!(out.outputs["total_valid_cells"], json!(4));

        let text = std::fs::read_to_string(&out_csv).unwrap();
        let _ = std::fs::remove_file(&out_csv);
        let mut lines = text.lines();
        assert_eq!(lines.next().unwrap(), "zone,10,20");
        assert_eq!(lines.next().unwrap(), "1,2,0");
        assert_eq!(lines.next().unwrap(), "2,1,1");
    }

    #[test]
    fn bins_mode_buckets_a_continuous_ramp() {
        // Single zone, values 0..9 (10 cells), 5 equal-width bins over [0,9]:
        // width 1.8, so each bin gets exactly 2 cells.
        let zones = path(raster(1, 10, vec![1.0; 10], -9999.0));
        let ramp: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let value = path(raster(1, 10, ramp, -9999.0));
        let out_csv = tmp_csv("bins");
        let out = run(json!({
            "zones": zones, "value": value, "output": out_csv,
            "mode": "bins", "bins": 5
        }));
        assert_eq!(out.outputs["mode"], json!("bins"));
        assert_eq!(out.outputs["class_count"], json!(5));

        let text = std::fs::read_to_string(&out_csv).unwrap();
        let _ = std::fs::remove_file(&out_csv);
        let mut lines = text.lines();
        let _header = lines.next().unwrap();
        let row = lines.next().unwrap();
        let counts: Vec<u64> = row.split(',').skip(1).map(|s| s.parse().unwrap()).collect();
        assert_eq!(counts, vec![2, 2, 2, 2, 2]);
        assert_eq!(counts.iter().sum::<u64>(), 10);
    }

    #[test]
    fn nodata_cells_are_skipped_from_either_raster() {
        // zone nodata at index 0, value nodata at index 3 -> only 2 valid cells.
        let zones = path(raster(1, 4, vec![-9999.0, 1.0, 1.0, 1.0], -9999.0));
        let value = path(raster(1, 4, vec![5.0, 5.0, 6.0, -1.0], -1.0));
        let out_csv = tmp_csv("nodata");
        let out = run(json!({ "zones": zones, "value": value, "output": out_csv }));
        assert_eq!(out.outputs["total_valid_cells"], json!(2));
        let text = std::fs::read_to_string(&out_csv).unwrap();
        let _ = std::fs::remove_file(&out_csv);
        assert!(text.contains("zone,5,6"), "csv was:\n{text}");
        assert!(text.contains("1,1,1"), "csv was:\n{text}");
    }

    #[test]
    fn percent_output_sums_to_one_hundred_per_zone() {
        let zones = path(raster(1, 4, vec![1.0, 1.0, 1.0, 1.0], -9999.0));
        let value = path(raster(1, 4, vec![10.0, 10.0, 10.0, 20.0], -9999.0));
        let out_csv = tmp_csv("percent");
        let out = run(json!({
            "zones": zones, "value": value, "output": out_csv, "percent": true
        }));
        let _ = out;
        let text = std::fs::read_to_string(&out_csv).unwrap();
        let _ = std::fs::remove_file(&out_csv);
        let row = text.lines().nth(1).unwrap();
        let pcts: Vec<f64> = row.split(',').skip(1).map(|s| s.parse().unwrap()).collect();
        let sum: f64 = pcts.iter().sum();
        assert!((sum - 100.0).abs() < 1e-6, "percentages: {pcts:?}");
        assert!((pcts[0] - 75.0).abs() < 1e-6);
        assert!((pcts[1] - 25.0).abs() < 1e-6);
    }

    #[test]
    fn long_output_matches_wide_output() {
        let zones = path(raster(2, 2, vec![1.0, 1.0, 2.0, 2.0], -9999.0));
        let value = path(raster(2, 2, vec![10.0, 10.0, 10.0, 20.0], -9999.0));
        let out_csv = tmp_csv("wide_long");
        let long_csv = tmp_csv("long");
        let out = run(json!({
            "zones": zones, "value": value, "output": out_csv, "long_output": long_csv
        }));
        assert_eq!(out.outputs["long_output"], json!(long_csv));
        let long_text = std::fs::read_to_string(&long_csv).unwrap();
        let _ = std::fs::remove_file(&out_csv);
        let _ = std::fs::remove_file(&long_csv);
        assert!(
            long_text.contains("1,10,2,100.0000"),
            "long csv:\n{long_text}"
        );
        assert!(
            long_text.contains("2,10,1,50.0000"),
            "long csv:\n{long_text}"
        );
        assert!(
            long_text.contains("2,20,1,50.0000"),
            "long csv:\n{long_text}"
        );
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = ZonalHistogramTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing everything");
        assert!(
            bad(json!({ "zones": "z.tif", "value": "v.tif" })).is_err(),
            "missing output"
        );
        assert!(
            bad(json!({ "zones": "z.tif", "value": "v.tif", "output": "o.csv", "mode": "bogus" }))
                .is_err(),
            "bad mode"
        );
        assert!(
            bad(json!({ "zones": "z.tif", "value": "v.tif", "output": "o.csv", "mode": "bins", "bins": 1 }))
                .is_err(),
            "bins too small"
        );
        assert!(
            bad(json!({ "zones": "z.tif", "value": "v.tif", "output": "o.csv", "zone_band": 0 }))
                .is_err(),
            "bad band"
        );
        assert!(bad(json!({ "zones": "z.tif", "value": "v.tif", "output": "o.csv" })).is_ok());
        assert!(bad(json!({
            "zones": "z.tif", "value": "v.tif", "output": "o.csv",
            "mode": "bins", "bins": 8, "percent": true
        }))
        .is_ok());
    }
}
