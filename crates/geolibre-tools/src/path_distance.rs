//! GeoLibre tool: cost distance with surface distance and vertical factors.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Path Distance* (Spatial Analyst).
//! The bundled `cost_distance` is planar and slope-blind; `corridor` already
//! ships an internal 8-connected Dijkstra, and this extends that same engine
//! with elevation-aware step costs so each move pays the true 3-D **surface
//! distance** (from a DEM) times a **vertical factor** (a slope-dependent
//! up/downhill penalty) times the friction cost — realistic travel-time /
//! energy surfaces for hiking, wildlife, and access modeling.
//!
//! Step cost between adjacent cells =
//! `surface_distance · vertical_factor(slope) · mean(friction)`, where
//! `surface_distance = √(planar_d² + Δz²)` and the vertical factor is one of:
//!
//! * `tobler` (default) — Tobler's hiking function, `exp(3.5·|S + 0.05|)` with
//!   `S = Δz / planar_d` (relative walking-time cost; anisotropic up vs down);
//! * `linear` — `max(zero_factor, 1 + slope_factor·S)` (uphill costs more);
//! * `sym_linear` — `1 + slope_factor·|S|` (both directions cost more);
//! * `inverse_linear` — `max(zero_factor, 1 − slope_factor·S)` (downhill costs more);
//! * `binary` — 1 within `max_slope` degrees, impassable beyond.
//!
//! Output is the accumulated least-cost surface (F32, no-data where
//! unreachable). Source cells start at 0.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{load_input_raster, parse_optional_output, raster_like_with_data};

const OUT_NODATA: f64 = -1.0;

#[derive(Clone, Copy, PartialEq)]
enum VerticalFactor {
    Tobler,
    Linear,
    SymLinear,
    InverseLinear,
    Binary,
}

pub struct PathDistanceTool;

impl Tool for PathDistanceTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "path_distance",
            display_name: "Path Distance",
            summary: "Accumulated least-cost distance where each step pays the true 3-D surface distance (from a DEM) times a slope-dependent vertical factor (Tobler hiking, linear, binary) times friction, like ArcGIS Path Distance.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "source",
                    description: "Source raster; cells with positive, non-nodata values are origins (accumulated cost 0).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output accumulated-cost raster. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cost",
                    description: "Optional friction/cost-of-passage raster (default uniform 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "surface",
                    description: "Optional elevation (DEM) raster; enables 3-D surface distance and slope-based vertical factors.",
                    required: false,
                },
                ToolParamSpec {
                    name: "vertical_factor",
                    description: "'tobler' (default), 'linear', 'sym_linear', 'inverse_linear', or 'binary'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "slope_factor",
                    description: "Slope penalty coefficient for the linear factors (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "zero_factor",
                    description: "Minimum (downhill) vertical factor floor for linear/inverse_linear (default 0.1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_slope",
                    description: "For 'binary': slopes steeper than this many degrees are impassable (default 30).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band for all rasters (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args
            .get("source")
            .and_then(Value::as_str)
            .map(str::trim)
            .unwrap_or("")
            .is_empty()
        {
            return Err(ToolError::Validation(
                "missing required string parameter 'source'".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let source_path = args.get("source").and_then(Value::as_str).unwrap();
        let output = parse_optional_output(args, "output")?;
        let prm = parse_params(args)?;

        let source = load_input_raster(source_path)?;
        let rows = source.rows;
        let cols = source.cols;
        let n = rows * cols;
        let band = prm.band;

        let cost = match &prm.cost_path {
            Some(p) => Some(load_input_raster(p)?),
            None => None,
        };
        let surface = match &prm.surface_path {
            Some(p) => Some(load_input_raster(p)?),
            None => None,
        };
        for (name, r) in [("cost", &cost), ("surface", &surface)] {
            if let Some(r) = r {
                if r.rows != rows || r.cols != cols {
                    return Err(ToolError::Validation(format!(
                        "'{name}' raster is {}x{}, expected {rows}x{cols}",
                        r.rows, r.cols
                    )));
                }
            }
        }

        // Materialize per-cell friction, elevation, and source flags.
        let s_nodata = source.nodata;
        let mut is_source = vec![false; n];
        let mut friction = vec![1.0_f64; n];
        let mut elev = vec![f64::NAN; n];
        let mut passable = vec![true; n];
        for row in 0..rows as isize {
            for col in 0..cols as isize {
                let i = row as usize * cols + col as usize;
                let sv = source.get(band, row, col);
                if sv != s_nodata && sv.is_finite() && sv > 0.0 {
                    is_source[i] = true;
                }
                if let Some(c) = &cost {
                    let f = c.get(band, row, col);
                    if f == c.nodata || !f.is_finite() {
                        passable[i] = false;
                    } else {
                        friction[i] = f.max(0.0);
                    }
                }
                if let Some(s) = &surface {
                    let z = s.get(band, row, col);
                    if z == s.nodata || !z.is_finite() {
                        passable[i] = false;
                    } else {
                        elev[i] = z;
                    }
                }
            }
        }

        ctx.progress.info("accumulating path distance from source");

        let cx = source.cell_size_x;
        let cy = source.cell_size_y;
        let d_orth = cx.min(cy);
        let d_diag = (cx * cx + cy * cy).sqrt();
        let neigh: [(isize, isize, f64); 8] = [
            (-1, 0, d_orth),
            (1, 0, d_orth),
            (0, -1, d_orth),
            (0, 1, d_orth),
            (-1, -1, d_diag),
            (-1, 1, d_diag),
            (1, -1, d_diag),
            (1, 1, d_diag),
        ];

        let mut dist = vec![f64::INFINITY; n];
        let mut done = vec![false; n];
        let mut heap: BinaryHeap<Node> = BinaryHeap::new();
        let mut source_count = 0usize;
        for i in 0..n {
            if is_source[i] && passable[i] {
                dist[i] = 0.0;
                heap.push(Node { cost: 0.0, idx: i });
                source_count += 1;
            }
        }
        if source_count == 0 {
            return Err(ToolError::Execution(
                "no valid source cells (positive and passable)".to_string(),
            ));
        }

        let use_surface = surface.is_some();
        while let Some(Node { cost: acc, idx }) = heap.pop() {
            if done[idx] {
                continue;
            }
            done[idx] = true;
            if !passable[idx] {
                continue;
            }
            let r = (idx / cols) as isize;
            let c = (idx % cols) as isize;
            for (dr, dc, planar) in neigh {
                let nr = r + dr;
                let nc = c + dc;
                if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                    continue;
                }
                let nidx = nr as usize * cols + nc as usize;
                if !passable[nidx] {
                    continue;
                }
                let dz = if use_surface {
                    elev[nidx] - elev[idx]
                } else {
                    0.0
                };
                let vf = vertical_factor(prm.vertical_factor, dz, planar, &prm);
                if !vf.is_finite() {
                    continue; // impassable slope
                }
                let surface_d = if use_surface {
                    (planar * planar + dz * dz).sqrt()
                } else {
                    planar
                };
                let mean_friction = 0.5 * (friction[idx] + friction[nidx]);
                let step = surface_d * vf * mean_friction;
                let nd = acc + step;
                if nd < dist[nidx] {
                    dist[nidx] = nd;
                    heap.push(Node {
                        cost: nd,
                        idx: nidx,
                    });
                }
            }
        }

        let data: Vec<f64> = dist
            .iter()
            .map(|&d| if d.is_finite() { d } else { OUT_NODATA })
            .collect();
        let out = raster_like_with_data(&source, data, OUT_NODATA, DataType::F32)?;
        let out_path = crate::common::write_or_store_output(out, output)?;

        let reachable = dist.iter().filter(|d| d.is_finite()).count();
        let max_cost = dist
            .iter()
            .cloned()
            .filter(|d| d.is_finite())
            .fold(0.0, f64::max);
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("source_cells".to_string(), json!(source_count));
        outputs.insert("reachable_cells".to_string(), json!(reachable));
        outputs.insert("max_cost".to_string(), json!(max_cost));
        Ok(ToolRunResult { outputs })
    }
}

/// The vertical (slope) factor for a move with vertical change `dz` over planar
/// distance `planar`. `S = dz / planar` is rise-over-run (tan of the slope).
fn vertical_factor(kind: VerticalFactor, dz: f64, planar: f64, prm: &Params) -> f64 {
    let s = if planar > 0.0 { dz / planar } else { 0.0 };
    match kind {
        // Tobler's hiking function: relative walking time ∝ exp(3.5·|S + 0.05|).
        VerticalFactor::Tobler => (3.5 * (s + 0.05).abs()).exp(),
        VerticalFactor::Linear => (1.0 + prm.slope_factor * s).max(prm.zero_factor),
        VerticalFactor::SymLinear => 1.0 + prm.slope_factor * s.abs(),
        VerticalFactor::InverseLinear => (1.0 - prm.slope_factor * s).max(prm.zero_factor),
        VerticalFactor::Binary => {
            let deg = s.atan().to_degrees().abs();
            if deg <= prm.max_slope {
                1.0
            } else {
                f64::INFINITY
            }
        }
    }
}

// ── Dijkstra node (min-heap via reversed Ord, as in corridor) ─────────────────

#[derive(Clone, Copy)]
struct Node {
    cost: f64,
    idx: usize,
}
impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.cost == other.cost
    }
}
impl Eq for Node {}
impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Node {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .cost
            .partial_cmp(&self.cost)
            .unwrap_or(Ordering::Equal)
            .then(self.idx.cmp(&other.idx))
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    cost_path: Option<String>,
    surface_path: Option<String>,
    vertical_factor: VerticalFactor,
    slope_factor: f64,
    zero_factor: f64,
    max_slope: f64,
    band: isize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let cost_path = str_param(args, "cost")?;
    let surface_path = str_param(args, "surface")?;
    let vertical_factor = match args
        .get("vertical_factor")
        .and_then(Value::as_str)
        .map(|s| s.trim().to_lowercase())
    {
        None => VerticalFactor::Tobler,
        Some(s) if s.is_empty() || s == "tobler" => VerticalFactor::Tobler,
        Some(s) if s == "linear" => VerticalFactor::Linear,
        Some(s) if s == "sym_linear" => VerticalFactor::SymLinear,
        Some(s) if s == "inverse_linear" => VerticalFactor::InverseLinear,
        Some(s) if s == "binary" => VerticalFactor::Binary,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "'vertical_factor' must be tobler|linear|sym_linear|inverse_linear|binary, got '{other}'"
            )))
        }
    };
    Ok(Params {
        cost_path,
        surface_path,
        vertical_factor,
        slope_factor: f64_param(args, "slope_factor")?.unwrap_or(1.0),
        zero_factor: f64_param(args, "zero_factor")?.unwrap_or(0.1),
        max_slope: f64_param(args, "max_slope")?.unwrap_or(30.0),
        band: args.get("band").and_then(Value::as_u64).unwrap_or(1).max(1) as isize - 1,
    })
}

fn str_param(args: &ToolArgs, key: &str) -> Result<Option<String>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => Ok(Some(s.clone())),
        Some(_) => Err(ToolError::Validation(format!("'{key}' must be a string"))),
    }
}

fn f64_param(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("'{key}' must be a number"))),
        Some(_) => Err(ToolError::Validation(format!("'{key}' must be a number"))),
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

    fn raster_of(rows: usize, cols: usize, vals: &[f64], nodata: f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata,
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
                r.set(0, row as isize, col as isize, vals[row * cols + col])
                    .unwrap();
            }
        }
        let id = memory_store::put_raster(r);
        memory_store::make_raster_memory_path(&id)
    }

    fn read(path: &str) -> (Vec<f64>, usize, usize) {
        let r = load_input_raster(path).unwrap();
        let mut v = Vec::new();
        for row in 0..r.rows as isize {
            for col in 0..r.cols as isize {
                v.push(r.get(0, row, col));
            }
        }
        (v, r.rows, r.cols)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Vec<f64>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = PathDistanceTool.run(&args, &ctx()).unwrap();
        let (v, _, _) = read(out.outputs["output"].as_str().unwrap());
        (out, v)
    }

    /// Flat surface, uniform friction, no DEM: planar cost distance from the
    /// left-column source grows by 1 per column (Tobler flat factor cancels in
    /// ratios but here we use no surface -> vf uses dz=0).
    #[test]
    fn planar_from_left_edge() {
        // 1x4 row, source at col 0. No surface -> planar distance, tobler flat.
        let src = raster_of(1, 4, &[1.0, 0.0, 0.0, 0.0], -1.0);
        let (out, v) = run(json!({ "source": src }));
        assert_eq!(out.outputs["source_cells"], json!(1));
        // With tobler and dz=0: vf = exp(3.5*0.05) constant k; cost[c] = c*1*k.
        let k = (3.5 * 0.05_f64).exp();
        for (c, &val) in v.iter().enumerate() {
            assert!(
                (val - c as f64 * k).abs() < 1e-6,
                "col {c}: {val} vs {}",
                c as f64 * k
            );
        }
    }

    /// A vertical cliff (impassable via binary max_slope) blocks propagation.
    #[test]
    fn binary_blocks_steep_slope() {
        // 1x3: elevations 0, 100, 0 with 1 m cells -> slope 89.4° between cells.
        let src = raster_of(1, 3, &[1.0, 0.0, 0.0], -1.0);
        let dem = raster_of(1, 3, &[0.0, 100.0, 0.0], -9999.0);
        let (out, v) = run(json!({
            "source": src, "surface": dem, "vertical_factor": "binary", "max_slope": 30.0
        }));
        // col 1 and 2 unreachable (cliff) -> nodata.
        assert_eq!(v[1], OUT_NODATA);
        assert_eq!(v[2], OUT_NODATA);
        assert_eq!(out.outputs["reachable_cells"], json!(1));
    }

    /// Surface distance exceeds planar when there is relief.
    #[test]
    fn surface_distance_adds_vertical() {
        // 1x2: flat friction, elevations 0 and 3, cells 4 m apart (set cell_size).
        // Use sym_linear with slope_factor 0 so vf=1 -> cost = surface distance.
        let src = raster_of(1, 2, &[1.0, 0.0], -1.0);
        let dem = raster_of(1, 2, &[0.0, 3.0], -9999.0);
        let (_out, v) = run(json!({
            "source": src, "surface": dem, "vertical_factor": "sym_linear", "slope_factor": 0.0
        }));
        // planar 1, dz 3 -> surface sqrt(1+9)=sqrt(10).
        assert!((v[1] - 10.0_f64.sqrt()).abs() < 1e-6, "{}", v[1]);
    }

    /// Tobler uphill costs more than downhill for the same |slope|.
    #[test]
    fn tobler_is_anisotropic() {
        // uphill move: dz +1 over planar 1 -> S=1
        let up = vertical_factor(VerticalFactor::Tobler, 1.0, 1.0, &dummy());
        // downhill move: dz -1 over planar 1 -> S=-1
        let down = vertical_factor(VerticalFactor::Tobler, -1.0, 1.0, &dummy());
        assert!(
            up > down,
            "uphill {up} should cost more than downhill {down}"
        );
    }

    fn dummy() -> Params {
        Params {
            cost_path: None,
            surface_path: None,
            vertical_factor: VerticalFactor::Tobler,
            slope_factor: 1.0,
            zero_factor: 0.1,
            max_slope: 30.0,
            band: 0,
        }
    }

    #[test]
    fn rejects_missing_source() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(PathDistanceTool.validate(&args).is_err());
    }
}
