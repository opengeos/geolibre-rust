//! GeoLibre tool: LandTrendr temporal segmentation of a yearly image series.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Analyze Changes Using LandTrendr*
//! (Image Analyst) — the reference method for forest-disturbance and land-change
//! mapping from annual Landsat/Sentinel composites. Nothing bundled comes close:
//! `change_vector_analysis` is two-date, and the GeoLibre `generate_trend_raster`
//! fits a single global trend. Builds on `spectral_index` (feed it NBR/NDVI
//! stacks) and the streaming raster-stack machinery from `generate_trend_raster`.
//!
//! Each pixel's yearly series is despiked, then segmented: vertices are inserted
//! greedily at the point of maximum deviation from the current piecewise-linear
//! fit (up to `max_segments` segments), and the segment with the greatest
//! disturbance (a drop in a vegetation index, `direction=loss`, or a rise for
//! `gain`) is reported. Outputs: the disturbance **year** (primary), and optional
//! **magnitude** and **duration** rasters. Pixels with fewer than `min_valid`
//! observations become no-data.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::common::{load_input_raster, parse_optional_output, raster_like_with_data};

pub struct LandtrendrTool;

impl Tool for LandtrendrTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "landtrendr",
            display_name: "LandTrendr",
            summary: "LandTrendr temporal segmentation of a yearly image series (like ArcGIS Analyze Changes Using LandTrendr): despike, greedy vertex-based piecewise-linear fitting, and extraction of the greatest disturbance (year, magnitude, duration) per pixel. The per-pixel change-history segmentation the two-date change_vector_analysis and single-trend generate_trend_raster can't do.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "inputs",
                    description: "Comma-separated yearly raster paths in time order (e.g. spectral_index NBR outputs; >= 4).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output raster: year of the greatest disturbance per pixel. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "years",
                    description: "Comma-separated years matching 'inputs' (default 0,1,2,...).",
                    required: false,
                },
                ToolParamSpec {
                    name: "magnitude_output",
                    description: "Optional raster of the disturbance magnitude (index change).",
                    required: false,
                },
                ToolParamSpec {
                    name: "duration_output",
                    description: "Optional raster of the disturbance duration (years).",
                    required: false,
                },
                ToolParamSpec {
                    name: "direction",
                    description: "Disturbance direction in the index: 'loss' (a drop; default) or 'gain' (a rise).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_segments",
                    description: "Maximum number of segments (default 6).",
                    required: false,
                },
                ToolParamSpec {
                    name: "spike_threshold",
                    description: "Despike strength 0..1 (fraction of the series range; default 0.75; 1 disables).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_valid",
                    description: "Minimum valid observations per pixel (default 4).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read from each input (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        let paths = parse_inputs(args)?;
        if paths.len() < 4 {
            return Err(ToolError::Validation(
                "'inputs' needs at least 4 yearly rasters".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let paths = parse_inputs(args)?;
        let prm = parse_params(args)?;
        let output = parse_optional_output(args, "output")?;
        let mag_out = parse_optional_output(args, "magnitude_output")?;
        let dur_out = parse_optional_output(args, "duration_output")?;

        let years = match &prm.years {
            Some(y) if y.len() != paths.len() => {
                return Err(ToolError::Validation(format!(
                    "{} years for {} rasters",
                    y.len(),
                    paths.len()
                )))
            }
            Some(y) => y.clone(),
            None => (0..paths.len()).map(|i| i as f64).collect(),
        };

        let rasters: Vec<Raster> = paths
            .iter()
            .map(|p| load_input_raster(p))
            .collect::<Result<_, _>>()?;
        let (rows, cols) = (rasters[0].rows, rasters[0].cols);
        let band = prm.band;
        for (i, r) in rasters.iter().enumerate() {
            if r.rows != rows || r.cols != cols {
                return Err(ToolError::Validation(format!("raster {i} size mismatch")));
            }
            if band < 0 || band as usize >= r.bands {
                return Err(ToolError::Validation(format!(
                    "band {} out of range for raster {i}",
                    band + 1
                )));
            }
        }

        ctx.progress.info(&format!(
            "LandTrendr over {} yearly raster(s), {rows}x{cols}",
            rasters.len()
        ));

        let nodata = rasters[0].nodata;
        let mut year_r = vec![nodata; rows * cols];
        let mut mag_r = vec![nodata; rows * cols];
        let mut dur_r = vec![nodata; rows * cols];

        let mut ty: Vec<f64> = Vec::with_capacity(rasters.len());
        let mut vv: Vec<f64> = Vec::with_capacity(rasters.len());
        for r in 0..rows {
            for c in 0..cols {
                ty.clear();
                vv.clear();
                for (k, ras) in rasters.iter().enumerate() {
                    let v = ras.get(band, r as isize, c as isize);
                    if v != ras.nodata && v.is_finite() {
                        ty.push(years[k]);
                        vv.push(v);
                    }
                }
                if vv.len() < prm.min_valid {
                    continue;
                }
                despike(&ty, &mut vv, prm.spike_threshold);
                let verts = fit_vertices(&ty, &vv, prm.max_segments);
                if let Some(dist) = greatest_disturbance(&ty, &vv, &verts, prm.gain) {
                    let idx = r * cols + c;
                    year_r[idx] = dist.year;
                    mag_r[idx] = dist.magnitude;
                    dur_r[idx] = dist.duration;
                }
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let out_r = raster_like_with_data(&rasters[0], year_r, nodata, DataType::F32)?;
        let out_path = crate::common::write_or_store_output(out_r, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        if let Some(p) = mag_out {
            let r = raster_like_with_data(&rasters[0], mag_r, nodata, DataType::F32)?;
            outputs.insert(
                "magnitude_output".to_string(),
                json!(crate::common::write_or_store_output(r, Some(p))?),
            );
        }
        if let Some(p) = dur_out {
            let r = raster_like_with_data(&rasters[0], dur_r, nodata, DataType::F32)?;
            outputs.insert(
                "duration_output".to_string(),
                json!(crate::common::write_or_store_output(r, Some(p))?),
            );
        }
        Ok(ToolRunResult { outputs })
    }
}

/// Removes single-year spikes: an interior point far from the line between its
/// neighbours (by more than `threshold × range`) is replaced by that line.
fn despike(t: &[f64], v: &mut [f64], threshold: f64) {
    if threshold >= 1.0 || v.len() < 3 {
        return;
    }
    let (lo, hi) = min_max(v);
    let range = (hi - lo).max(1e-12);
    for i in 1..v.len() - 1 {
        let interp = lerp(t[i - 1], v[i - 1], t[i + 1], v[i + 1], t[i]);
        if (v[i] - interp).abs() > threshold * range {
            v[i] = interp;
        }
    }
}

/// Greedy vertex selection: start with the endpoints, repeatedly insert the
/// point of maximum deviation from the current piecewise-linear fit until
/// `max_segments` segments are reached. Returns sorted vertex indices.
fn fit_vertices(t: &[f64], v: &[f64], max_segments: usize) -> Vec<usize> {
    let n = v.len();
    let mut verts = vec![0usize, n - 1];
    while verts.len() <= max_segments {
        // For each segment (consecutive vertex pair), find the interior point of
        // max deviation from the segment line.
        let mut best_dev = 0.0;
        let mut best_pt = None;
        for w in verts.windows(2) {
            let (a, b) = (w[0], w[1]);
            for i in (a + 1)..b {
                let fit = lerp(t[a], v[a], t[b], v[b], t[i]);
                let dev = (v[i] - fit).abs();
                if dev > best_dev {
                    best_dev = dev;
                    best_pt = Some(i);
                }
            }
        }
        match best_pt {
            Some(p) if best_dev > 0.0 => {
                verts.push(p);
                verts.sort_unstable();
            }
            _ => break,
        }
    }
    verts
}

struct Disturbance {
    year: f64,
    magnitude: f64,
    duration: f64,
}

/// The segment (between consecutive vertices) with the greatest disturbance: for
/// `gain=false` (loss) the largest drop in value, for `gain=true` the largest
/// rise. Magnitude is the absolute value change; year is the segment start.
fn greatest_disturbance(t: &[f64], v: &[f64], verts: &[usize], gain: bool) -> Option<Disturbance> {
    let mut best: Option<Disturbance> = None;
    let mut best_change = 0.0f64;
    for w in verts.windows(2) {
        let (a, b) = (w[0], w[1]);
        let change = v[b] - v[a]; // positive = increase
        let disturbance = if gain { change } else { -change }; // positive = the sought direction
        if disturbance > best_change {
            best_change = disturbance;
            best = Some(Disturbance {
                year: t[a],
                magnitude: change.abs(),
                duration: t[b] - t[a],
            });
        }
    }
    best
}

fn lerp(x0: f64, y0: f64, x1: f64, y1: f64, x: f64) -> f64 {
    if x1 == x0 {
        (y0 + y1) / 2.0
    } else {
        y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    }
}

fn min_max(v: &[f64]) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &x in v {
        lo = lo.min(x);
        hi = hi.max(x);
    }
    (lo, hi)
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    years: Option<Vec<f64>>,
    gain: bool,
    max_segments: usize,
    spike_threshold: f64,
    min_valid: usize,
    band: isize,
}

fn parse_inputs(args: &ToolArgs) -> Result<Vec<String>, ToolError> {
    let s = args
        .get("inputs")
        .and_then(Value::as_str)
        .ok_or_else(|| ToolError::Validation("missing required parameter 'inputs'".to_string()))?;
    let paths: Vec<String> = s
        .split(',')
        .map(str::trim)
        .filter(|p| !p.is_empty())
        .map(String::from)
        .collect();
    if paths.is_empty() {
        return Err(ToolError::Validation("'inputs' is empty".to_string()));
    }
    Ok(paths)
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let years = match args.get("years").and_then(Value::as_str) {
        None => None,
        Some(s) if s.trim().is_empty() => None,
        Some(s) => Some(
            s.split(',')
                .map(str::trim)
                .filter(|x| !x.is_empty())
                .map(|x| {
                    x.parse::<f64>()
                        .map_err(|_| ToolError::Validation(format!("year '{x}' is not a number")))
                })
                .collect::<Result<Vec<_>, _>>()?,
        ),
    };
    let gain = match args.get("direction").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("loss") => false,
        Some("gain") => true,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'direction' must be 'loss' or 'gain', got '{o}'"
            )))
        }
    };
    let max_segments = match args.get("max_segments") {
        None | Some(Value::Null) => 6,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(6).max(1) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 6,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'max_segments' must be an integer".into()))?
            .max(1),
        _ => 6,
    };
    let spike_threshold = match args.get("spike_threshold") {
        None | Some(Value::Null) => 0.75,
        Some(Value::Number(n)) => n.as_f64().unwrap_or(0.75).clamp(0.0, 1.0),
        Some(Value::String(s)) if s.trim().is_empty() => 0.75,
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'spike_threshold' must be a number".into()))?
            .clamp(0.0, 1.0),
        _ => 0.75,
    };
    let min_valid = match args.get("min_valid") {
        None | Some(Value::Null) => 4,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(4).max(3) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 4,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'min_valid' must be an integer".into()))?
            .max(3),
        _ => 4,
    };
    let band_1based = args.get("band").and_then(Value::as_u64).unwrap_or(1).max(1);
    Ok(Params {
        years,
        gain,
        max_segments,
        spike_threshold,
        min_valid,
        band: (band_1based - 1) as isize,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, RasterConfig};

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

    /// A stable-then-drop-then-recover NBR series pinpoints the disturbance year.
    #[test]
    fn detects_disturbance_year() {
        // Years 2000..2007. NBR high (0.8) then a crash at 2003 to 0.1, slow recovery.
        let vals = [0.8, 0.8, 0.8, 0.1, 0.25, 0.4, 0.55, 0.7];
        let years = "2000,2001,2002,2003,2004,2005,2006,2007";
        let rasters: Vec<String> = vals.iter().map(|&v| raster_from(1, 1, vec![v])).collect();
        let args: ToolArgs = serde_json::from_value(json!({
            "inputs": rasters.join(","), "years": years, "direction": "loss",
        }))
        .unwrap();
        let out = LandtrendrTool.run(&args, &ctx()).unwrap();
        let yr = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        assert_eq!(
            yr.get(0, 0, 0),
            2002.0,
            "the disturbance segment starts at 2002 (drop into 2003)"
        );
    }

    /// The disturbance function itself picks the biggest drop.
    #[test]
    fn greatest_disturbance_picks_biggest_drop() {
        let t: Vec<f64> = (0..6).map(|i| i as f64).collect();
        let v = vec![1.0, 0.9, 0.85, 0.2, 0.3, 0.35]; // big drop 0.85 -> 0.2
        let verts = fit_vertices(&t, &v, 6);
        let d = greatest_disturbance(&t, &v, &verts, false).unwrap();
        assert!(
            d.magnitude > 0.5,
            "should catch the ~0.65 drop, got {}",
            d.magnitude
        );
    }

    /// Despiking removes a one-year spike.
    #[test]
    fn despike_removes_spike() {
        let t: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let mut v = vec![0.8, 0.8, 0.1, 0.8, 0.8]; // single-year dropout
        despike(&t, &mut v, 0.5);
        assert!(
            (v[2] - 0.8).abs() < 1e-9,
            "the spike should be interpolated back to ~0.8"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            LandtrendrTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "inputs": "a.tif,b.tif,c.tif" })).is_err()); // < 4
        assert!(bad(json!({ "inputs": "a.tif,b.tif,c.tif,d.tif", "direction": "up" })).is_err());
        assert!(bad(json!({ "inputs": "a.tif,b.tif,c.tif,d.tif" })).is_ok());
    }
}
