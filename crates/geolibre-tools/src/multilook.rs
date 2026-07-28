//! GeoLibre tool: incoherent range/azimuth averaging of SAR imagery.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Multilook* (Image Analyst).
//!
//! SAR is the largest domain gap in the catalog: neither registry ships any
//! SAR/InSAR tool. That is not for lack of adjacent capability — the bundled
//! suite has the speckle *filters* (`lee_filter`, `enhanced_lee_filter`,
//! `refined_lee_filter`) — but those operate on already-detected intensity.
//! Multilooking happens **upstream** of them, at the point where complex SLC
//! data becomes intensity, and nothing in either registry does that step. So the
//! speckle filters currently have no in-catalog path to the data they were
//! designed for.
//!
//! ## Complex-input convention
//!
//! The repo's raster model is real-valued, so single-look complex data arrives
//! as a **two-band I/Q pair** (band 1 = in-phase, band 2 = quadrature). This
//! convention is established here and inherited by the rest of the SAR chain
//! (see `sar_coherence`).
//!
//! ## The one correctness trap
//!
//! Detection (`I² + Q²`) must happen **before** averaging. Averaging complex
//! samples is *coherent* summation and gives a completely different result;
//! incoherent averaging of detected intensity is what actually reduces speckle.
//! `coherent_vs_incoherent_averaging_differ` pins this down.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster, RasterConfig};

use crate::common::{load_input_raster, parse_optional_output, write_or_store_output};

/// Reduces SAR speckle by incoherently averaging adjacent looks.
pub struct MultilookTool;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Units {
    Amplitude,
    Intensity,
    Db,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Stat {
    Mean,
    Median,
}

impl Tool for MultilookTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "multilook",
            display_name: "Multilook",
            summary: "Reduces SAR speckle by incoherently averaging adjacent looks in range and azimuth, producing square-pixel intensity imagery (ArcGIS Multilook). Neither registry ships any SAR tool; the bundled lee_filter family operates on already-detected intensity, so this supplies the upstream step that turns complex single-look data into the detected imagery those filters expect. Complex input is read as a two-band I/Q pair; detection happens before averaging, and the equivalent number of looks is reported for downstream speckle filtering.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input SAR raster: a two-band I/Q pair for complex (SLC) data, or a single band of already-detected data (see 'input_domain').",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output raster path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "input_domain",
                    description: "For non-complex input, what band 1 already holds: 'intensity' (default) or 'amplitude'. Amplitude is squared to intensity before averaging, since incoherent averaging is only correct in the intensity domain.",
                    required: false,
                },
                ToolParamSpec {
                    name: "complex",
                    description: "Treat the input as a complex I/Q pair (default: true when the raster has 2+ bands, false otherwise).",
                    required: false,
                },
                ToolParamSpec {
                    name: "range_looks",
                    description: "Number of looks in range (x). Default 1, or derived when auto_looks is set.",
                    required: false,
                },
                ToolParamSpec {
                    name: "azimuth_looks",
                    description: "Number of looks in azimuth (y). Default 1, or derived when auto_looks is set.",
                    required: false,
                },
                ToolParamSpec {
                    name: "auto_looks",
                    description: "If true (default) and looks are not given explicitly, derive them from the pixel-spacing ratio so output ground pixels are approximately square.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_units",
                    description: "amplitude, intensity (default), or db.",
                    required: false,
                },
                ToolParamSpec {
                    name: "statistic",
                    description: "mean (default) or median; median is more robust where strong scatterers dominate a look window.",
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
        parse_units(args)?;
        parse_stat(args)?;
        parse_input_domain(args)?;
        for k in ["range_looks", "azimuth_looks"] {
            if let Some(v) = opt_u64(args, k)? {
                if v == 0 {
                    return Err(ToolError::Validation(format!(
                        "parameter '{k}' must be at least 1"
                    )));
                }
            }
        }
        opt_bool(args, "complex")?;
        opt_bool(args, "auto_looks")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args.get("input").and_then(Value::as_str).ok_or_else(|| {
            ToolError::Validation("missing required parameter 'input'".to_string())
        })?;
        let output = parse_optional_output(args, "output")?;
        let units = parse_units(args)?;
        let stat = parse_stat(args)?;
        let amplitude_input = parse_input_domain(args)?;

        let raster = load_input_raster(input)?;
        let complex = opt_bool(args, "complex")?.unwrap_or(raster.bands >= 2);
        if complex && raster.bands < 2 {
            return Err(ToolError::Validation(
                "complex input needs two bands (I and Q); got 1".to_string(),
            ));
        }

        let rows = raster.rows;
        let cols = raster.cols;
        let nodata = raster.nodata;
        let cx = raster.cell_size_x.abs().max(f64::MIN_POSITIVE);
        let cy = raster.cell_size_y.abs().max(f64::MIN_POSITIVE);

        // Look counts: explicit values win; otherwise square the ground pixel.
        let auto = opt_bool(args, "auto_looks")?.unwrap_or(true);
        let (range_looks, azimuth_looks) = match (
            opt_u64(args, "range_looks")?,
            opt_u64(args, "azimuth_looks")?,
        ) {
            (Some(r), Some(a)) => (r as usize, a as usize),
            (r, a) if auto => {
                // Choose integer looks whose ratio best squares the pixel:
                // looks_x * cx ≈ looks_y * cy.
                let (dr, da) = square_pixel_looks(cx, cy);
                (
                    r.map(|v| v as usize).unwrap_or(dr),
                    a.map(|v| v as usize).unwrap_or(da),
                )
            }
            (r, a) => (
                r.map(|v| v as usize).unwrap_or(1),
                a.map(|v| v as usize).unwrap_or(1),
            ),
        };
        let range_looks = range_looks.max(1);
        let azimuth_looks = azimuth_looks.max(1);

        ctx.progress.info(&format!(
            "multilooking {range_looks} x {azimuth_looks} look(s)"
        ));

        // ── Detection: intensity BEFORE any averaging ────────────────────────
        let mut intensity = vec![f64::NAN; rows * cols];
        for row in 0..rows {
            for col in 0..cols {
                let i_val = raster.get(0, row as isize, col as isize);
                if i_val == nodata || !i_val.is_finite() {
                    continue;
                }
                let v = if complex {
                    let q_val = raster.get(1, row as isize, col as isize);
                    if q_val == nodata || !q_val.is_finite() {
                        continue;
                    }
                    i_val * i_val + q_val * q_val
                } else if amplitude_input {
                    // Amplitude must be squared first: averaging amplitudes is
                    // the wrong domain, and the dB conversion below assumes
                    // intensity (10*log10(I), not 20*log10(A)).
                    i_val * i_val
                } else {
                    i_val
                };
                intensity[row * cols + col] = v;
            }
        }

        // ── Incoherent averaging over non-overlapping look windows ───────────
        let out_rows = rows.div_ceil(azimuth_looks);
        let out_cols = cols.div_ceil(range_looks);
        let out_nodata = -9999.0_f64;
        let mut data = vec![out_nodata; out_rows * out_cols];
        let mut valid_windows = 0_u64;

        for orow in 0..out_rows {
            for ocol in 0..out_cols {
                let mut vals: Vec<f64> = Vec::with_capacity(range_looks * azimuth_looks);
                for r in (orow * azimuth_looks)..((orow + 1) * azimuth_looks).min(rows) {
                    for c in (ocol * range_looks)..((ocol + 1) * range_looks).min(cols) {
                        let v = intensity[r * cols + c];
                        if v.is_finite() {
                            vals.push(v);
                        }
                    }
                }
                if vals.is_empty() {
                    continue;
                }
                valid_windows += 1;
                let mean_intensity = match stat {
                    Stat::Mean => vals.iter().sum::<f64>() / vals.len() as f64,
                    Stat::Median => {
                        vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                        let mid = vals.len() / 2;
                        if vals.len() % 2 == 1 {
                            vals[mid]
                        } else {
                            (vals[mid - 1] + vals[mid]) / 2.0
                        }
                    }
                };
                data[orow * out_cols + ocol] = match units {
                    Units::Intensity => mean_intensity,
                    Units::Amplitude => mean_intensity.max(0.0).sqrt(),
                    // Guard log10 against zero/negative: writing -inf would
                    // poison every downstream statistic.
                    Units::Db => {
                        if mean_intensity > 0.0 {
                            10.0 * mean_intensity.log10()
                        } else {
                            out_nodata
                        }
                    }
                };
            }
            ctx.progress
                .progress((orow as f64 + 1.0) / out_rows.max(1) as f64);
        }

        // Output cell size scales by the look counts; georeferencing is kept.
        // div_ceil pads the southern edge when the look window does not divide
        // the row count, so drop the origin by the remainder — otherwise y_max
        // would drift north of the input extent.
        let y_pad = (out_rows * azimuth_looks - rows) as f64 * cy;
        let mut out = Raster::new(RasterConfig {
            cols: out_cols,
            rows: out_rows,
            bands: 1,
            x_min: raster.x_min,
            y_min: raster.y_min - y_pad,
            cell_size: cx * range_looks as f64,
            cell_size_y: Some(cy * azimuth_looks as f64),
            nodata: out_nodata,
            data_type: DataType::F32,
            crs: raster.crs.clone(),
            metadata: raster.metadata.clone(),
        });
        for row in 0..out_rows {
            for col in 0..out_cols {
                out.set(0, row as isize, col as isize, data[row * out_cols + col])
                    .map_err(|e| ToolError::Execution(format!("failed writing cell: {e}")))?;
            }
        }

        let out_path = write_or_store_output(out, output)?;

        // The equivalent number of looks is what the bundled speckle filters
        // take as a parameter, so report it explicitly.
        let enl = (range_looks * azimuth_looks) as f64;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("range_looks".to_string(), json!(range_looks));
        outputs.insert("azimuth_looks".to_string(), json!(azimuth_looks));
        outputs.insert("equivalent_number_of_looks".to_string(), json!(enl));
        outputs.insert("output_rows".to_string(), json!(out_rows));
        outputs.insert("output_cols".to_string(), json!(out_cols));
        outputs.insert("valid_windows".to_string(), json!(valid_windows));
        Ok(ToolRunResult { outputs })
    }
}

/// Picks small integer look counts whose ratio best squares the ground pixel:
/// `looks_x * cell_x ≈ looks_y * cell_y`.
fn square_pixel_looks(cell_x: f64, cell_y: f64) -> (usize, usize) {
    if !cell_x.is_finite() || !cell_y.is_finite() || cell_x <= 0.0 || cell_y <= 0.0 {
        return (1, 1);
    }
    let target = cell_y / cell_x; // looks_x / looks_y
    let mut best = (1_usize, 1_usize);
    let mut best_err = f64::INFINITY;
    for ly in 1..=8_usize {
        for lx in 1..=8_usize {
            let err = ((lx as f64 / ly as f64) - target).abs();
            // Prefer the smallest look counts among equally good ratios.
            if err < best_err - 1e-12 {
                best_err = err;
                best = (lx, ly);
            }
        }
    }
    best
}

fn parse_units(args: &ToolArgs) -> Result<Units, ToolError> {
    match opt_str(args, "output_units")? {
        None => Ok(Units::Intensity),
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "intensity" => Ok(Units::Intensity),
            "amplitude" => Ok(Units::Amplitude),
            "db" => Ok(Units::Db),
            other => Err(ToolError::Validation(format!(
                "unknown output_units '{other}' (expected amplitude, intensity or db)"
            ))),
        },
    }
}

/// Returns `true` when a non-complex input holds amplitude rather than intensity.
fn parse_input_domain(args: &ToolArgs) -> Result<bool, ToolError> {
    match opt_str(args, "input_domain")? {
        None => Ok(false),
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "intensity" => Ok(false),
            "amplitude" => Ok(true),
            other => Err(ToolError::Validation(format!(
                "unknown input_domain '{other}' (expected intensity or amplitude)"
            ))),
        },
    }
}

fn parse_stat(args: &ToolArgs) -> Result<Stat, ToolError> {
    match opt_str(args, "statistic")? {
        None => Ok(Stat::Mean),
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "mean" => Ok(Stat::Mean),
            "median" => Ok(Stat::Median),
            other => Err(ToolError::Validation(format!(
                "unknown statistic '{other}' (expected mean or median)"
            ))),
        },
    }
}

fn opt_str<'a>(args: &'a ToolArgs, key: &str) -> Result<Option<&'a str>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => Ok(Some(s.as_str())),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a string when provided"
        ))),
    }
}

fn opt_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => n.as_u64().map(Some).ok_or_else(|| {
            ToolError::Validation(format!("parameter '{key}' must be a positive integer"))
        }),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s.trim().parse::<u64>().map(Some).map_err(|_| {
            ToolError::Validation(format!("parameter '{key}' must be a positive integer"))
        }),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a positive integer"
        ))),
    }
}

fn opt_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
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
    use wbraster::CrsInfo;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// `data` is band-major: all of band 0, then all of band 1.
    fn raster_cells(
        cols: usize,
        rows: usize,
        bands: usize,
        data: &[f64],
        cell_x: f64,
        cell_y: f64,
    ) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: cell_x,
            cell_size_y: Some(cell_y),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for b in 0..bands {
            for row in 0..rows {
                for col in 0..cols {
                    r.set(
                        b as isize,
                        row as isize,
                        col as isize,
                        data[b * rows * cols + row * cols + col],
                    )
                    .unwrap();
                }
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn raster(cols: usize, rows: usize, bands: usize, data: &[f64]) -> String {
        raster_cells(cols, rows, bands, data, 1.0, 1.0)
    }

    fn run_with(extra: Value) -> (Raster, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(extra).unwrap();
        let res = MultilookTool.run(&args, &ctx()).unwrap();
        let r = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        (r, res)
    }

    /// Detection squares I and Q before averaging.
    #[test]
    fn detects_complex_input_to_intensity() {
        // Single pixel, I=3, Q=4 => intensity 25, amplitude 5.
        let path = raster(1, 1, 2, &[3.0, 4.0]);
        let (out, _) = run_with(json!({ "input": path.clone() }));
        assert!((out.get(0, 0, 0) - 25.0).abs() < 1e-6);

        let (amp, _) = run_with(json!({ "input": path, "output_units": "amplitude" }));
        assert!((amp.get(0, 0, 0) - 5.0).abs() < 1e-6);
    }

    /// The correctness trap: detection must precede averaging. Two samples with
    /// opposite phase cancel under *coherent* summation but not under the
    /// incoherent averaging this tool performs.
    #[test]
    fn coherent_vs_incoherent_averaging_differ() {
        // Two pixels: (I,Q) = (1,0) and (-1,0). Coherent mean is 0 => intensity
        // 0. Incoherent mean of intensities is (1 + 1)/2 = 1.
        let path = raster(2, 1, 2, &[1.0, -1.0, 0.0, 0.0]);
        let (out, res) = run_with(json!({ "input": path, "range_looks": 2, "azimuth_looks": 1 }));
        assert_eq!(res.outputs["output_cols"], json!(1));
        assert!(
            (out.get(0, 0, 0) - 1.0).abs() < 1e-6,
            "incoherent averaging must give 1.0, not 0.0 (got {})",
            out.get(0, 0, 0)
        );
    }

    /// Output is decimated by the look counts and the cell size scales to match,
    /// so the raster still covers the same ground extent.
    #[test]
    fn decimates_and_rescales_geotransform() {
        let path = raster(4, 4, 1, &[4.0; 16]);
        let (out, res) = run_with(json!({ "input": path, "range_looks": 2, "azimuth_looks": 2 }));
        assert_eq!(out.rows, 2);
        assert_eq!(out.cols, 2);
        assert!((out.cell_size_x - 2.0).abs() < 1e-9);
        assert!((out.cell_size_y - 2.0).abs() < 1e-9);
        assert_eq!(res.outputs["equivalent_number_of_looks"], json!(4.0));
    }

    /// Averaging reduces speckle: the variance of the multilooked image is
    /// lower than the single-look input's. This is the whole purpose.
    #[test]
    fn averaging_reduces_variance() {
        // Alternating high/low intensity, single band already detected.
        let data: Vec<f64> = (0..64)
            .map(|i| if i % 2 == 0 { 10.0 } else { 0.0 })
            .collect();
        let path = raster(8, 8, 1, &data);
        let (out, _) = run_with(
            json!({ "input": path, "complex": false, "range_looks": 2, "azimuth_looks": 2 }),
        );

        let input_var = variance(&data);
        let mut vals = Vec::new();
        for r in 0..out.rows {
            for c in 0..out.cols {
                vals.push(out.get(0, r as isize, c as isize));
            }
        }
        let out_var = variance(&vals);
        assert!(
            out_var < input_var,
            "multilooking must reduce variance: {out_var} !< {input_var}"
        );
    }

    fn variance(v: &[f64]) -> f64 {
        let m = v.iter().sum::<f64>() / v.len() as f64;
        v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / v.len() as f64
    }

    /// The median statistic is exercised, including the even-count average.
    #[test]
    fn median_statistic_differs_from_mean() {
        // One 3-cell window with an outlier: mean 34.33, median 2.
        let path = raster(3, 1, 1, &[1.0, 2.0, 100.0]);
        let (mean, _) = run_with(
            json!({ "input": path.clone(), "complex": false, "range_looks": 3, "azimuth_looks": 1 }),
        );
        let (median, _) = run_with(
            json!({ "input": path, "complex": false, "range_looks": 3, "azimuth_looks": 1, "statistic": "median" }),
        );
        assert!(
            (mean.get(0, 0, 0) - 103.0 / 3.0).abs() < 1e-4,
            "got {}",
            mean.get(0, 0, 0)
        );
        assert!(
            (median.get(0, 0, 0) - 2.0).abs() < 1e-6,
            "got {}",
            median.get(0, 0, 0)
        );

        // Even count takes the average of the middle pair: (2 + 100) / 2 = 51.
        let even = raster(2, 1, 1, &[2.0, 100.0]);
        let (m, _) = run_with(
            json!({ "input": even, "complex": false, "range_looks": 2, "azimuth_looks": 1, "statistic": "median" }),
        );
        assert!(
            (m.get(0, 0, 0) - 51.0).abs() < 1e-6,
            "got {}",
            m.get(0, 0, 0)
        );
    }

    /// Amplitude input is squared to intensity before averaging.
    #[test]
    fn amplitude_input_is_squared_before_averaging() {
        // Amplitudes 3 and 4 -> intensities 9 and 16 -> mean 12.5.
        let path = raster(2, 1, 1, &[3.0, 4.0]);
        let (out, _) = run_with(json!({
            "input": path, "complex": false, "input_domain": "amplitude",
            "range_looks": 2, "azimuth_looks": 1
        }));
        assert!(
            (out.get(0, 0, 0) - 12.5).abs() < 1e-5,
            "amplitudes must be squared before averaging, got {}",
            out.get(0, 0, 0)
        );
    }

    /// auto_looks squares the ground pixel from the spacing ratio.
    #[test]
    fn auto_looks_squares_the_ground_pixel() {
        // Range spacing 2, azimuth spacing 8 => needs 4 range looks per azimuth
        // look so 4*2 == 1*8.
        let path = raster_cells(8, 8, 1, &vec![1.0; 64], 2.0, 8.0);
        let (out, res) = run_with(json!({ "input": path, "complex": false }));
        assert_eq!(res.outputs["range_looks"], json!(4));
        assert_eq!(res.outputs["azimuth_looks"], json!(1));
        assert!(
            (out.cell_size_x - out.cell_size_y).abs() < 1e-9,
            "output pixels should be square: {} vs {}",
            out.cell_size_x,
            out.cell_size_y
        );
    }

    /// dB of a zero-intensity window is written as no-data, not -inf.
    #[test]
    fn db_units_guard_against_zero() {
        let path = raster(1, 1, 1, &[0.0]);
        let (out, _) = run_with(json!({ "input": path, "complex": false, "output_units": "db" }));
        let v = out.get(0, 0, 0);
        assert!(v.is_finite(), "must not emit -inf, got {v}");
        assert_eq!(v, out.nodata);
    }

    #[test]
    fn rejects_bad_parameters() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(MultilookTool.validate(&args).is_err());

        let path = raster(2, 2, 1, &[1.0; 4]);
        for bad in [
            json!({ "input": path.clone(), "output_units": "watts" }),
            json!({ "input": path.clone(), "statistic": "mode" }),
            json!({ "input": path.clone(), "range_looks": 0 }),
            json!({ "input": path.clone(), "input_domain": "power" }),
        ] {
            let args: ToolArgs = serde_json::from_value(bad).unwrap();
            assert!(MultilookTool.validate(&args).is_err());
        }

        // Asking for complex on a single-band raster is caught at run time.
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": path, "complex": true })).unwrap();
        assert!(MultilookTool.run(&args, &ctx()).is_err());
    }
}
