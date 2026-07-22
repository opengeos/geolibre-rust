//! GeoLibre tool: relative-risk ratio of two kernel density surfaces.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Calculate Kernel Density Ratio*
//! (Spatial Analyst). The bundled `heat_map` produces a single-layer density,
//! but raw case density is misleading wherever the at-risk population varies.
//! Nothing in the ~791-tool suite produces a density *ratio* with shared
//! bandwidth handling and stability controls — the standard relative-risk
//! surface in epidemiology (cases vs. population) and crime analysis
//! (incidents vs. exposure).
//!
//! Both point sets are evaluated with the SAME quartic (biweight) kernel and a
//! shared `bandwidth`, then divided cell-by-cell:
//! `ratio = numerator_density / denominator_density`. A `denominator_floor`
//! guards the division: cells whose denominator density falls below the floor
//! (or is zero) are written as no-data so no `inf`/`NaN` leaks into the output.
//! With `log_ratio` the result is natural-log transformed so over- and
//! under-representation are symmetric around 0.
//!
//! The output raster spans the union bounding box of both layers, padded by the
//! bandwidth so kernels near the edge are not clipped. Distances are haversine
//! metres for a geographic CRS (EPSG:4326 or untagged) and planar CRS units
//! otherwise, matching the rest of the movement/analysis family.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
use wbvector::{Geometry, Layer};

use crate::common::{parse_optional_output, write_or_store_output};
use crate::vector_common::{load_input_layer, parse_optional_str};

/// No-data sentinel for the ratio raster (finite, so `is_finite()` stays true).
const NODATA: f64 = -9999.0;

pub struct KernelDensityRatioTool;

impl Tool for KernelDensityRatioTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "kernel_density_ratio",
            display_name: "Kernel Density Ratio",
            summary: "Relative-risk ratio of two kernel density surfaces: evaluate a numerator point layer (cases) and a denominator point layer (population at risk) with the same quartic kernel and shared bandwidth, then divide cell-by-cell — like ArcGIS Calculate Kernel Density Ratio. A denominator floor writes no-data where the population is too thin (no inf/NaN leaks), and an optional natural-log transform makes over/under-representation symmetric. The density-ratio the single-layer heat_map cannot express.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Numerator point layer (e.g. cases / incidents).",
                    required: true,
                },
                ToolParamSpec {
                    name: "denominator",
                    description: "Denominator point layer (e.g. population / exposure at risk).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output ratio raster path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "weight_field",
                    description: "Numeric field weighting numerator points (default: each point counts 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "denominator_weight_field",
                    description: "Numeric field weighting denominator points (default: each point counts 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "bandwidth",
                    description: "Shared kernel bandwidth / search radius (metres for a geographic CRS, CRS units otherwise). Default: Silverman rule from the combined point spread.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Output cell size in CRS units (default: larger extent side / 200).",
                    required: false,
                },
                ToolParamSpec {
                    name: "log_ratio",
                    description: "If true, output the natural log of the ratio (symmetric around 0). Default false.",
                    required: false,
                },
                ToolParamSpec {
                    name: "denominator_floor",
                    description: "Cells whose denominator density is below this value (or zero) become no-data, avoiding division blow-ups. Default 0 (only exactly-zero denominator is masked).",
                    required: false,
                },
                ToolParamSpec {
                    name: "epsg",
                    description: "EPSG to tag the output (default: from the numerator layer).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "denominator")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let num_path = require_str(args, "input")?;
        let den_path = require_str(args, "denominator")?;
        let output = parse_optional_output(args, "output")?;
        let prm = parse_params(args)?;

        let num_layer = load_input_layer(num_path)?;
        let den_layer = load_input_layer(den_path)?;

        // Geographic (haversine metres) when the numerator CRS is lon/lat.
        let geographic = num_layer.crs_epsg().map(|e| e == 4326).unwrap_or(true);

        let num_pts = collect_points(&num_layer, prm.weight_field.as_deref())?;
        let den_pts = collect_points(&den_layer, prm.den_weight_field.as_deref())?;
        if num_pts.is_empty() {
            return Err(ToolError::Execution(
                "numerator layer has no point features".to_string(),
            ));
        }
        if den_pts.is_empty() {
            return Err(ToolError::Execution(
                "denominator layer has no point features".to_string(),
            ));
        }

        // Bandwidth: explicit, or a Silverman rule-of-thumb on the combined set.
        let bandwidth = match prm.bandwidth {
            Some(b) => b,
            None => default_bandwidth(&num_pts, &den_pts, geographic),
        };
        if !(bandwidth > 0.0 && bandwidth.is_finite()) {
            return Err(ToolError::Execution(
                "computed bandwidth is not positive; pass an explicit 'bandwidth'".to_string(),
            ));
        }

        // Union bbox of both layers, padded by the bandwidth (in CRS units).
        let (mut xmin, mut ymin, mut xmax, mut ymax) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for p in num_pts.iter().chain(den_pts.iter()) {
            xmin = xmin.min(p.x);
            xmax = xmax.max(p.x);
            ymin = ymin.min(p.y);
            ymax = ymax.max(p.y);
        }
        let mid_lat = (ymin + ymax) * 0.5;
        // Bandwidth converted to CRS units for padding / cell-range culling.
        let (pad_x, pad_y) = bandwidth_crs_units(bandwidth, geographic, mid_lat);
        xmin -= pad_x;
        xmax += pad_x;
        ymin -= pad_y;
        ymax += pad_y;
        let ext_w = (xmax - xmin).max(1e-12);
        let ext_h = (ymax - ymin).max(1e-12);

        let cell = prm
            .cell_size
            .unwrap_or((ext_w.max(ext_h) / 200.0).max(1e-12));
        let cols = ((ext_w / cell).ceil() as usize).max(1);
        let rows = ((ext_h / cell).ceil() as usize).max(1);

        ctx.progress.info(&format!(
            "{} numerator + {} denominator point(s) -> {rows}x{cols} ratio raster (bandwidth {:.4})",
            num_pts.len(),
            den_pts.len(),
            bandwidth
        ));

        // Splat each point's quartic kernel onto the two density grids.
        let mut num_dens = vec![0.0f64; rows * cols];
        let mut den_dens = vec![0.0f64; rows * cols];
        splat(
            &num_pts,
            &mut num_dens,
            rows,
            cols,
            xmin,
            ymax,
            cell,
            bandwidth,
            geographic,
        );
        splat(
            &den_pts,
            &mut den_dens,
            rows,
            cols,
            xmin,
            ymax,
            cell,
            bandwidth,
            geographic,
        );

        // Divide with the stability floor.
        let floor = prm.denominator_floor;
        let mut data = vec![NODATA; rows * cols];
        let mut valid = 0usize;
        let mut masked = 0usize;
        let mut min_v = f64::INFINITY;
        let mut max_v = f64::NEG_INFINITY;
        let mut sum_v = 0.0f64;
        for i in 0..rows * cols {
            let d = den_dens[i];
            if d <= 0.0 || d < floor {
                masked += 1;
                continue;
            }
            let mut r = num_dens[i] / d;
            if prm.log_ratio {
                // ln(0) guard: numerator can legitimately be 0 (no cases here).
                r = if r > 0.0 { r.ln() } else { NODATA };
                if r == NODATA {
                    masked += 1;
                    continue;
                }
            }
            if !r.is_finite() {
                masked += 1;
                continue;
            }
            data[i] = r;
            valid += 1;
            min_v = min_v.min(r);
            max_v = max_v.max(r);
            sum_v += r;
        }
        if valid == 0 {
            return Err(ToolError::Execution(
                "no valid cells: the denominator density is below the floor everywhere".to_string(),
            ));
        }
        let mean_v = sum_v / valid as f64;

        let out = build_raster(
            data,
            rows,
            cols,
            xmin,
            ymin,
            cell,
            prm.epsg.or_else(|| num_layer.crs_epsg()),
        )?;
        let out_path = write_or_store_output(out, output)?;

        ctx.progress.info(&format!(
            "{valid} valid cell(s), {masked} masked (below floor); ratio min {min_v:.4} mean {mean_v:.4} max {max_v:.4}"
        ));

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        outputs.insert("bandwidth".to_string(), json!(bandwidth));
        outputs.insert("valid_cells".to_string(), json!(valid));
        outputs.insert("masked_cells".to_string(), json!(masked));
        outputs.insert("ratio_min".to_string(), json!(min_v));
        outputs.insert("ratio_mean".to_string(), json!(mean_v));
        outputs.insert("ratio_max".to_string(), json!(max_v));
        Ok(ToolRunResult { outputs })
    }
}

// ── Point collection ─────────────────────────────────────────────────────────

struct Pt {
    x: f64,
    y: f64,
    w: f64,
}

/// Reads point features and their (optional) weights from a layer.
fn collect_points(layer: &Layer, weight_field: Option<&str>) -> Result<Vec<Pt>, ToolError> {
    let widx = match weight_field {
        Some(f) => Some(
            layer
                .schema
                .field_index(f)
                .ok_or_else(|| ToolError::Validation(format!("weight field '{f}' not found")))?,
        ),
        None => None,
    };
    let mut pts = Vec::new();
    for feat in layer.iter() {
        let Some((x, y)) = feat.geometry.as_ref().and_then(point_xy) else {
            continue;
        };
        let w = match widx {
            Some(i) => feat
                .attributes
                .get(i)
                .and_then(|v| v.as_f64())
                .unwrap_or(1.0),
            None => 1.0,
        };
        if w > 0.0 && w.is_finite() {
            pts.push(Pt { x, y, w });
        }
    }
    Ok(pts)
}

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        _ => None,
    }
}

// ── Kernel density ───────────────────────────────────────────────────────────

/// Accumulates each point's quartic (biweight) kernel onto the density grid.
/// The output grid is row-major with row 0 at `y_top` (north) descending.
#[allow(clippy::too_many_arguments)]
fn splat(
    pts: &[Pt],
    grid: &mut [f64],
    rows: usize,
    cols: usize,
    x_min: f64,
    y_top: f64,
    cell: f64,
    h: f64,
    geographic: bool,
) {
    // 2D quartic normalization: K(u) = (3/pi) (1-u^2)^2 for u<=1, over h^2.
    let norm = 3.0 / (std::f64::consts::PI * h * h);
    for p in pts {
        // Cell-index window that could fall within the bandwidth of this point.
        let mid_lat = p.y;
        let (pad_x, pad_y) = bandwidth_crs_units(h, geographic, mid_lat);
        let c_lo = (((p.x - pad_x - x_min) / cell).floor() as isize).max(0) as usize;
        let c_hi =
            (((p.x + pad_x - x_min) / cell).ceil() as isize).clamp(0, cols as isize) as usize;
        // Row 0 is the north edge; larger row -> smaller y.
        let r_lo = (((y_top - (p.y + pad_y)) / cell).floor() as isize).max(0) as usize;
        let r_hi =
            (((y_top - (p.y - pad_y)) / cell).ceil() as isize).clamp(0, rows as isize) as usize;
        for r in r_lo..r_hi {
            let cy = y_top - (r as f64 + 0.5) * cell;
            for c in c_lo..c_hi {
                let cx = x_min + (c as f64 + 0.5) * cell;
                let d = distance(cx, cy, p.x, p.y, geographic);
                let u = d / h;
                if u < 1.0 {
                    let k = (1.0 - u * u).powi(2);
                    grid[r * cols + c] += p.w * norm * k;
                }
            }
        }
    }
}

/// Silverman rule-of-thumb bandwidth from the combined point spread: the
/// standard distance about the mean center times `n^(-1/6)` (2-D optimal).
fn default_bandwidth(a: &[Pt], b: &[Pt], geographic: bool) -> f64 {
    let n = a.len() + b.len();
    if n == 0 {
        return 1.0;
    }
    let (mut mx, mut my) = (0.0, 0.0);
    for p in a.iter().chain(b.iter()) {
        mx += p.x;
        my += p.y;
    }
    mx /= n as f64;
    my /= n as f64;
    let mut sq = 0.0;
    for p in a.iter().chain(b.iter()) {
        let d = distance(mx, my, p.x, p.y, geographic);
        sq += d * d;
    }
    let sd = (sq / n as f64).sqrt();
    let bw = sd * (n as f64).powf(-1.0 / 6.0);
    if bw > 0.0 && bw.is_finite() {
        bw
    } else {
        1.0
    }
}

/// Converts a bandwidth in distance units into (x, y) padding in CRS units.
/// For a geographic CRS the bandwidth is metres, so it is turned back into
/// degrees (latitude-corrected in x); for a projected CRS it is already in
/// CRS units.
fn bandwidth_crs_units(h: f64, geographic: bool, lat: f64) -> (f64, f64) {
    if geographic {
        const M_PER_DEG: f64 = 111_320.0;
        let pad_y = h / M_PER_DEG;
        let coslat = lat.to_radians().cos().abs().max(1e-6);
        let pad_x = h / (M_PER_DEG * coslat);
        (pad_x, pad_y)
    } else {
        (h, h)
    }
}

fn distance(x0: f64, y0: f64, x1: f64, y1: f64, geographic: bool) -> f64 {
    if geographic {
        haversine(y0, x0, y1, x1)
    } else {
        (x1 - x0).hypot(y1 - y0)
    }
}

fn haversine(lat0: f64, lon0: f64, lat1: f64, lon1: f64) -> f64 {
    const R: f64 = 6_371_000.0;
    let (p0, p1) = (lat0.to_radians(), lat1.to_radians());
    let dphi = (lat1 - lat0).to_radians();
    let dlmb = (lon1 - lon0).to_radians();
    let a = (dphi / 2.0).sin().powi(2) + p0.cos() * p1.cos() * (dlmb / 2.0).sin().powi(2);
    2.0 * R * a.sqrt().asin()
}

// ── Raster construction ──────────────────────────────────────────────────────

fn build_raster(
    data: Vec<f64>,
    rows: usize,
    cols: usize,
    x_min: f64,
    y_min: f64,
    cell: f64,
    epsg: Option<u32>,
) -> Result<Raster, ToolError> {
    let mut out = Raster::new(RasterConfig {
        cols,
        rows,
        bands: 1,
        x_min,
        y_min,
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
    for r in 0..rows {
        for c in 0..cols {
            out.set(0, r as isize, c as isize, data[r * cols + c])
                .map_err(|e| ToolError::Execution(format!("write failed: {e}")))?;
        }
    }
    Ok(out)
}

// ── Parameters ───────────────────────────────────────────────────────────────

struct Params {
    weight_field: Option<String>,
    den_weight_field: Option<String>,
    bandwidth: Option<f64>,
    cell_size: Option<f64>,
    log_ratio: bool,
    denominator_floor: f64,
    epsg: Option<u32>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let weight_field = parse_optional_str(args, "weight_field")?.map(String::from);
    let den_weight_field = parse_optional_str(args, "denominator_weight_field")?.map(String::from);
    let bandwidth = opt_pos(args, "bandwidth")?;
    let cell_size = opt_pos(args, "cell_size")?;
    let log_ratio = opt_bool(args, "log_ratio")?;
    let denominator_floor = opt_f64(args, "denominator_floor")?.unwrap_or(0.0).max(0.0);
    let epsg = opt_u32(args, "epsg")?;
    Ok(Params {
        weight_field,
        den_weight_field,
        bandwidth,
        cell_size,
        log_ratio,
        denominator_floor,
        epsg,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn opt_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
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

fn opt_pos(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match opt_f64(args, key)? {
        Some(v) if v > 0.0 && v.is_finite() => Ok(Some(v)),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a positive number"
        ))),
        None => Ok(None),
    }
}

fn opt_u32(args: &ToolArgs, key: &str) -> Result<Option<u32>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64().map(|v| v as u32)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<u32>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be an integer"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be an integer"
        ))),
    }
}

fn opt_bool(args: &ToolArgs, key: &str) -> Result<bool, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(false),
        Some(Value::Bool(b)) => Ok(*b),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(false),
            "true" | "1" | "yes" => Ok(true),
            "false" | "0" | "no" => Ok(false),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(Value::Number(n)) => Ok(n.as_f64().map(|v| v != 0.0).unwrap_or(false)),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a projected (EPSG:3857) point layer from (x, y) pairs.
    fn pt_layer(pts: &[(f64, f64)]) -> String {
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("w", FieldType::Float));
        for (x, y) in pts {
            l.add_feature(Some(Geometry::point(*x, *y)), &[("w", (1.0f64).into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = KernelDensityRatioTool.run(&args, &ctx()).unwrap();
        let r = crate::common::load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, r)
    }

    /// Sample the raster value nearest a world coordinate (skips no-data).
    fn value_at(r: &Raster, x: f64, y: f64) -> Option<f64> {
        let col = ((x - r.x_min) / r.cell_size_x).floor() as isize;
        // Row 0 is the north edge.
        let y_top = r.y_min + r.rows as f64 * r.cell_size_y;
        let row = ((y_top - y) / r.cell_size_y).floor() as isize;
        if row < 0 || col < 0 || row >= r.rows as isize || col >= r.cols as isize {
            return None;
        }
        let v = r.get(0, row, col);
        if v == r.nodata {
            None
        } else {
            Some(v)
        }
    }

    /// Identical numerator and denominator layers -> ratio is 1 everywhere valid.
    #[test]
    fn identical_layers_give_unit_ratio() {
        let pts = [(0.0, 0.0), (100.0, 0.0), (0.0, 100.0), (100.0, 100.0)];
        let layer = pt_layer(&pts);
        let (out, r) = run(json!({
            "input": layer, "denominator": layer,
            "bandwidth": 60.0, "cell_size": 10.0,
        }));
        assert!(out.outputs["valid_cells"].as_u64().unwrap() > 0);
        let mut checked = 0;
        for row in 0..r.rows {
            for col in 0..r.cols {
                let v = r.get(0, row as isize, col as isize);
                if v != r.nodata {
                    assert!(
                        (v - 1.0).abs() < 1e-9,
                        "identical layers -> ratio 1, got {v}"
                    );
                    checked += 1;
                }
            }
        }
        assert!(checked > 10);
    }

    /// Where the numerator concentrates, the ratio exceeds 1; elsewhere it stays
    /// near 1. Numerator = denominator plus an extra cluster at one corner.
    #[test]
    fn concentration_raises_ratio() {
        let den = [(0.0, 0.0), (100.0, 0.0), (0.0, 100.0), (100.0, 100.0)];
        let num = [
            (0.0, 0.0),
            (100.0, 0.0),
            (0.0, 100.0),
            (100.0, 100.0),
            // extra cases piled onto the (0,0) corner
            (2.0, 2.0),
            (-2.0, 3.0),
            (3.0, -1.0),
        ];
        let (_o, r) = run(json!({
            "input": pt_layer(&num), "denominator": pt_layer(&den),
            "bandwidth": 40.0, "cell_size": 5.0,
        }));
        let hot = value_at(&r, 0.0, 0.0).expect("cell near cluster");
        let cold = value_at(&r, 100.0, 100.0).expect("cell far from cluster");
        assert!(
            hot > 1.5 && hot > cold * 1.3,
            "ratio should spike near the extra cluster: hot={hot}, cold={cold}"
        );
    }

    /// A denominator floor masks thin-population cells and NO inf/NaN leaks:
    /// the denominator points sit only in the left half, so right-half cells
    /// fall below the floor and become no-data.
    #[test]
    fn floor_masks_and_no_nonfinite() {
        let den = [(0.0, 0.0), (0.0, 50.0), (0.0, 100.0)]; // left edge only
        let num = [(0.0, 0.0), (200.0, 50.0)]; // one case far to the right
        let (out, r) = run(json!({
            "input": pt_layer(&num), "denominator": pt_layer(&den),
            "bandwidth": 30.0, "cell_size": 10.0, "denominator_floor": 1e-9,
        }));
        assert!(out.outputs["masked_cells"].as_u64().unwrap() > 0);
        let mut nodata = 0;
        for row in 0..r.rows {
            for col in 0..r.cols {
                let v = r.get(0, row as isize, col as isize);
                assert!(v.is_finite(), "no inf/NaN in output, got {v}");
                if v == r.nodata {
                    nodata += 1;
                }
            }
        }
        assert!(nodata > 0, "far-from-denominator cells should be masked");
        // The far-right case has no denominator support -> masked, not inf.
        assert!(value_at(&r, 200.0, 50.0).is_none());
    }

    /// log_ratio of identical layers is ~0 (ln 1).
    #[test]
    fn log_ratio_of_identical_is_zero() {
        let pts = [(0.0, 0.0), (100.0, 0.0), (0.0, 100.0), (100.0, 100.0)];
        let layer = pt_layer(&pts);
        let (_o, r) = run(json!({
            "input": layer, "denominator": layer,
            "bandwidth": 60.0, "cell_size": 10.0, "log_ratio": true,
        }));
        for row in 0..r.rows {
            for col in 0..r.cols {
                let v = r.get(0, row as isize, col as isize);
                if v != r.nodata {
                    assert!(v.abs() < 1e-9, "ln(1) should be 0, got {v}");
                }
            }
        }
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            KernelDensityRatioTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err()); // missing denominator
        assert!(bad(
            json!({ "input": "a.geojson", "denominator": "b.geojson", "bandwidth": -5.0 })
        )
        .is_err());
        assert!(bad(json!({ "input": "a.geojson", "denominator": "b.geojson" })).is_ok());
    }
}
