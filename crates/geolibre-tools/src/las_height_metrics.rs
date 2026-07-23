//! GeoLibre tool: canopy height-metric raster from lidar heights.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *LAS Height Metrics* (3D Analyst). For
//! forestry / canopy-structure work it rasterizes, per cell, a suite of
//! distributional statistics of the lidar point **heights**: `mean`, `std`,
//! `skewness`, `kurtosis` (excess), median absolute deviation (`mad`), and any
//! number of height `percentiles`.
//!
//! The bundled lidar suite covers pieces but not this grid: `height_above_ground`
//! normalizes heights and `lidar_point_stats` gives basic per-cell
//! count/elevation stats, but neither emits the skewness/kurtosis/MAD +
//! multi-percentile canopy-metrics stack in one pass.
//!
//! Input heights should be **above-ground** (run `height_above_ground` /
//! `normalize_lidar` first); the tool uses each point's Z as its height. Points
//! below `min_height` (ground/noise) are excluded; cells with fewer than
//! `min_points` are no-data. Output is a multi-band raster, one band per
//! requested metric then per percentile (ascending); the band order is returned
//! in `bands`. `cell_size` and heights are in the cloud's CRS units.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};

use crate::common::{parse_optional_output, write_or_store_output};
use crate::lidar_common::load_input_cloud;

const NODATA: f64 = -9999.0;

pub struct LasHeightMetricsTool;

impl Tool for LasHeightMetricsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "las_height_metrics",
            display_name: "LAS Height Metrics",
            summary: "Rasterize per-cell canopy height metrics from lidar heights (like ArcGIS LAS Height Metrics): mean, std, skewness, kurtosis, MAD, and height percentiles as a multi-band grid — the forestry moment/percentile suite the bundled height_above_ground and lidar_point_stats don't emit in one pass.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input LAS/LAZ point cloud of above-ground heights (normalize first with height_above_ground).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output multi-band raster path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "metrics",
                    description: "Comma-separated subset of mean,std,skewness,kurtosis,mad (default 'mean,std,skewness,kurtosis,mad').",
                    required: false,
                },
                ToolParamSpec {
                    name: "height_percentiles",
                    description: "Comma-separated percentiles in [0,100], e.g. '50,95,99' (default '50,95,99'; empty for none).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_height",
                    description: "Exclude points with height below this (ground/noise cutoff). Default 0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_points",
                    description: "Cells with fewer valid points become no-data. Default 4.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Output cell size in CRS units. Default 10.",
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

        let cloud = load_input_cloud(input)?;
        let n_pts = cloud.points.len();
        if n_pts == 0 {
            return Err(ToolError::Execution("point cloud is empty".to_string()));
        }

        // Bounds over included points (height >= min_height).
        let (mut min_x, mut min_y) = (f64::INFINITY, f64::INFINITY);
        let (mut max_x, mut max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
        let mut n_incl = 0usize;
        for p in &cloud.points {
            if p.z < prm.min_height {
                continue;
            }
            min_x = min_x.min(p.x);
            min_y = min_y.min(p.y);
            max_x = max_x.max(p.x);
            max_y = max_y.max(p.y);
            n_incl += 1;
        }
        if n_incl == 0 {
            return Err(ToolError::Execution(format!(
                "no points at or above min_height ({})",
                prm.min_height
            )));
        }

        let cell = prm.cell_size;
        let cols = (((max_x - min_x) / cell).ceil() as usize).max(1);
        let rows = (((max_y - min_y) / cell).ceil() as usize).max(1);
        let n_cells = rows * cols;

        ctx.progress.info(&format!(
            "{n_incl}/{n_pts} point(s) above {:.2}; {rows}x{cols} grid @ {cell}",
            prm.min_height
        ));

        // Bin heights per cell.
        let mut bins: Vec<Vec<f64>> = vec![Vec::new(); n_cells];
        for p in &cloud.points {
            if p.z < prm.min_height {
                continue;
            }
            let mut col = ((p.x - min_x) / cell).floor() as isize;
            let mut row_b = ((p.y - min_y) / cell).floor() as isize;
            col = col.clamp(0, cols as isize - 1);
            row_b = row_b.clamp(0, rows as isize - 1);
            let row = rows as isize - 1 - row_b; // row 0 = north (top)
            bins[(row as usize) * cols + col as usize].push(p.z);
        }

        // Compute each requested band.
        let band_names = prm.band_names();
        let n_bands = band_names.len();
        let mut band_data: Vec<Vec<f64>> = vec![vec![NODATA; n_cells]; n_bands];
        let mut filled = 0usize;
        for (cell_i, heights) in bins.iter_mut().enumerate() {
            if heights.len() < prm.min_points {
                continue;
            }
            filled += 1;
            let stats = CellStats::new(heights);
            let mut b = 0;
            for m in &prm.metrics {
                band_data[b][cell_i] = stats.metric(*m);
                b += 1;
            }
            for &p in &prm.percentiles {
                band_data[b][cell_i] = stats.percentile(p);
                b += 1;
            }
        }

        ctx.progress.info(&format!(
            "{filled}/{n_cells} cell(s) filled, {n_bands} band(s)"
        ));

        // Build the multi-band raster.
        let epsg = cloud.crs.as_ref().and_then(|c| c.epsg);
        let mut raster = Raster::new(RasterConfig {
            cols,
            rows,
            bands: n_bands,
            x_min: min_x,
            y_min: min_y,
            cell_size: cell,
            cell_size_y: Some(cell),
            nodata: NODATA,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg,
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for (bi, data) in band_data.iter().enumerate() {
            for r in 0..rows {
                for c in 0..cols {
                    raster
                        .set(bi as isize, r as isize, c as isize, data[r * cols + c])
                        .map_err(|e| ToolError::Execution(format!("failed writing cell: {e}")))?;
                }
            }
        }

        let out_path = write_or_store_output(raster, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("bands".to_string(), json!(band_names));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        outputs.insert("cells_filled".to_string(), json!(filled));
        outputs.insert("points_used".to_string(), json!(n_incl));
        Ok(ToolRunResult { outputs })
    }
}

// ── Per-cell statistics ───────────────────────────────────────────────────────

struct CellStats {
    n: f64,
    mean: f64,
    std: f64,
    m3: f64,
    m4: f64,
    sorted: Vec<f64>,
}

impl CellStats {
    /// Consumes the cell's heights (sorts in place for percentiles/MAD).
    fn new(heights: &mut [f64]) -> CellStats {
        let n = heights.len() as f64;
        let mean = heights.iter().sum::<f64>() / n;
        let (mut s2, mut s3, mut s4) = (0.0, 0.0, 0.0);
        for &h in heights.iter() {
            let d = h - mean;
            let d2 = d * d;
            s2 += d2;
            s3 += d2 * d;
            s4 += d2 * d2;
        }
        let var = s2 / n;
        heights.sort_by(f64::total_cmp);
        CellStats {
            n,
            mean,
            std: var.sqrt(),
            m3: s3 / n,
            m4: s4 / n,
            sorted: heights.to_vec(),
        }
    }

    fn metric(&self, m: Metric) -> f64 {
        match m {
            Metric::Mean => self.mean,
            Metric::Std => self.std,
            Metric::Skewness => {
                if self.std <= 0.0 {
                    0.0
                } else {
                    self.m3 / self.std.powi(3)
                }
            }
            Metric::Kurtosis => {
                if self.std <= 0.0 {
                    0.0
                } else {
                    self.m4 / self.std.powi(4) - 3.0 // excess kurtosis
                }
            }
            Metric::Mad => {
                let med = percentile_sorted(&self.sorted, 50.0);
                let mut dev: Vec<f64> = self.sorted.iter().map(|h| (h - med).abs()).collect();
                dev.sort_by(f64::total_cmp);
                percentile_sorted(&dev, 50.0)
            }
        }
    }

    fn percentile(&self, p: f64) -> f64 {
        let _ = self.n; // n retained for clarity/debugging
        percentile_sorted(&self.sorted, p)
    }
}

/// Linear-interpolated percentile of an already-sorted slice (numpy default).
fn percentile_sorted(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return NODATA;
    }
    if n == 1 {
        return sorted[0];
    }
    let rank = (p / 100.0) * (n - 1) as f64;
    let lo = rank.floor() as usize;
    let hi = rank.ceil() as usize;
    if lo == hi {
        sorted[lo]
    } else {
        let frac = rank - lo as f64;
        sorted[lo] * (1.0 - frac) + sorted[hi] * frac
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Metric {
    Mean,
    Std,
    Skewness,
    Kurtosis,
    Mad,
}

impl Metric {
    fn name(&self) -> &'static str {
        match self {
            Metric::Mean => "mean",
            Metric::Std => "std",
            Metric::Skewness => "skewness",
            Metric::Kurtosis => "kurtosis",
            Metric::Mad => "mad",
        }
    }
    fn parse(s: &str) -> Option<Metric> {
        match s.trim().to_ascii_lowercase().as_str() {
            "mean" => Some(Metric::Mean),
            "std" | "stddev" | "standard_deviation" => Some(Metric::Std),
            "skewness" | "skew" => Some(Metric::Skewness),
            "kurtosis" | "kurt" => Some(Metric::Kurtosis),
            "mad" | "median_absolute_deviation" => Some(Metric::Mad),
            _ => None,
        }
    }
}

struct Params {
    metrics: Vec<Metric>,
    percentiles: Vec<f64>,
    min_height: f64,
    min_points: usize,
    cell_size: f64,
}

impl Params {
    fn band_names(&self) -> Vec<String> {
        let mut names: Vec<String> = self.metrics.iter().map(|m| m.name().to_string()).collect();
        for &p in &self.percentiles {
            // p50 / p95 / p99; strip trailing .0 for whole percentiles.
            if (p.fract()).abs() < 1e-9 {
                names.push(format!("p{}", p as i64));
            } else {
                names.push(format!("p{p}"));
            }
        }
        names
    }
}

fn parse_optional_f64(args: &ToolArgs, key: &str, default: f64) -> Result<f64, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(default),
        Some(Value::Number(n)) => Ok(n.as_f64().unwrap_or(default)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(default),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be a number"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a number"
        ))),
    }
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let metrics = match args.get("metrics").and_then(Value::as_str) {
        None => vec![
            Metric::Mean,
            Metric::Std,
            Metric::Skewness,
            Metric::Kurtosis,
            Metric::Mad,
        ],
        Some(s) if s.trim().is_empty() => vec![
            Metric::Mean,
            Metric::Std,
            Metric::Skewness,
            Metric::Kurtosis,
            Metric::Mad,
        ],
        Some(s) => {
            let mut v = Vec::new();
            for tok in s.split(',').map(str::trim).filter(|t| !t.is_empty()) {
                v.push(Metric::parse(tok).ok_or_else(|| {
                    ToolError::Validation(format!(
                        "unknown metric '{tok}' (expected mean,std,skewness,kurtosis,mad)"
                    ))
                })?);
            }
            v
        }
    };

    let percentiles = match args.get("height_percentiles").and_then(Value::as_str) {
        None => vec![50.0, 95.0, 99.0],
        Some(s) if s.trim().is_empty() => Vec::new(),
        Some(s) => {
            let mut v = Vec::new();
            for tok in s.split(',').map(str::trim).filter(|t| !t.is_empty()) {
                let p: f64 = tok.parse().map_err(|_| {
                    ToolError::Validation(format!("percentile '{tok}' is not a number"))
                })?;
                if !(0.0..=100.0).contains(&p) {
                    return Err(ToolError::Validation(format!(
                        "percentile {p} out of range [0, 100]"
                    )));
                }
                v.push(p);
            }
            v
        }
    };

    if metrics.is_empty() && percentiles.is_empty() {
        return Err(ToolError::Validation(
            "no metrics or percentiles requested".to_string(),
        ));
    }

    let min_height = parse_optional_f64(args, "min_height", 0.0)?;
    let min_points = match args.get("min_points") {
        None | Some(Value::Null) => 4,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(4).max(1) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 4,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'min_points' must be an integer".into()))?
            .max(1),
        Some(_) => {
            return Err(ToolError::Validation(
                "'min_points' must be a number".into(),
            ))
        }
    };
    let cell_size = parse_optional_f64(args, "cell_size", 10.0)?;
    if !(cell_size > 0.0 && cell_size.is_finite()) {
        return Err(ToolError::Validation(
            "'cell_size' must be a positive number".to_string(),
        ));
    }

    Ok(Params {
        metrics,
        percentiles,
        min_height,
        min_points,
        cell_size,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wblidar::{memory_store, PointCloud, PointRecord};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn cloud_of(pts: &[(f64, f64, f64)]) -> String {
        let mut cloud = PointCloud::default();
        cloud.crs = Some(wblidar::Crs {
            epsg: Some(32610),
            wkt: None,
        });
        for &(x, y, z) in pts {
            let mut p = PointRecord::default();
            p.x = x;
            p.y = y;
            p.z = z;
            cloud.points.push(p);
        }
        let id = memory_store::put_lidar(cloud);
        memory_store::make_lidar_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = LasHeightMetricsTool.run(&args, &ctx()).unwrap();
        let r = crate::common::load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, r)
    }

    /// One cell full of known heights -> mean/std/percentile bands are exact.
    #[test]
    fn single_cell_moments() {
        // All points inside one 100x100 cell; heights 2,4,4,4,5,5,7,9.
        let hs = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let pts: Vec<(f64, f64, f64)> = hs.iter().map(|&z| (5.0, 5.0, z)).collect();
        let (out, r) = run(json!({
            "input": cloud_of(&pts), "cell_size": 100.0,
            "metrics": "mean,std", "height_percentiles": "50", "min_points": 1
        }));
        let names: Vec<String> = serde_json::from_value(out.outputs["bands"].clone()).unwrap();
        assert_eq!(names, vec!["mean", "std", "p50"]);
        // mean 5, pop std 2, median 4.5.
        assert!((r.get(0, 0, 0) - 5.0).abs() < 1e-5);
        assert!((r.get(1, 0, 0) - 2.0).abs() < 1e-5);
        assert!((r.get(2, 0, 0) - 4.5).abs() < 1e-5);
    }

    /// Cells below min_points are no-data.
    #[test]
    fn min_points_gate() {
        // 2 points in one cell, min_points 4 -> that cell is nodata.
        let (_o, r) = run(json!({
            "input": cloud_of(&[(1.0, 1.0, 3.0), (1.0, 1.0, 5.0)]),
            "cell_size": 100.0, "metrics": "mean", "height_percentiles": "", "min_points": 4
        }));
        assert_eq!(r.get(0, 0, 0), NODATA);
    }

    /// min_height excludes ground points from the statistics.
    #[test]
    fn min_height_excludes_ground() {
        // heights 0,0,0 (ground) + 10,10,10,10 (canopy); min_height 1 keeps canopy.
        let mut pts = vec![(2.0, 2.0, 0.0); 3];
        pts.extend(vec![(2.0, 2.0, 10.0); 4]);
        let (_o, r) = run(json!({
            "input": cloud_of(&pts), "cell_size": 100.0,
            "metrics": "mean", "height_percentiles": "", "min_height": 1.0, "min_points": 1
        }));
        assert!(
            (r.get(0, 0, 0) - 10.0).abs() < 1e-5,
            "mean should be canopy-only 10"
        );
    }

    /// Skewness sign is recovered (right-skewed heights -> positive skew).
    #[test]
    fn skewness_sign() {
        // right-skewed: many low, few high.
        let mut hs = vec![1.0; 10];
        hs.push(20.0);
        let pts: Vec<(f64, f64, f64)> = hs.iter().map(|&z| (3.0, 3.0, z)).collect();
        let (_o, r) = run(json!({
            "input": cloud_of(&pts), "cell_size": 100.0,
            "metrics": "skewness", "height_percentiles": "", "min_points": 1
        }));
        assert!(
            r.get(0, 0, 0) > 0.5,
            "right-skewed heights -> positive skewness"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            LasHeightMetricsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.laz", "metrics": "bogus" })).is_err());
        assert!(bad(json!({ "input": "a.laz", "height_percentiles": "150" })).is_err());
        assert!(bad(json!({ "input": "a.laz", "cell_size": 0 })).is_err());
        assert!(bad(json!({ "input": "a.laz" })).is_ok());
    }
}
