//! GeoLibre tool: least-cost corridor between two accumulated-cost surfaces.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Corridor* (Spatial Analyst). The
//! bundled cost suite (`cost_distance`, `cost_allocation`, `cost_pathway`)
//! produces a single optimal path; wildlife-corridor and route-planning
//! workflows need the *corridor band* — the swath of near-optimal routes — which
//! no bundled tool computes.
//!
//! The corridor surface is the cell-wise sum of two accumulated-cost rasters,
//! one accumulated from source A and one from source B: each cell's value is the
//! total cost of the cheapest path A→cell→B that passes through it. The global
//! minimum of that surface is the least-cost-path cost; cells at (or near) the
//! minimum form the corridor. Thresholding — an absolute cost or a percentage
//! above the minimum — turns the surface into the near-optimal band.
//!
//! Two ways to supply the inputs:
//!
//! * **Direct** — `cost1` and `cost2`, two accumulated-cost rasters (e.g. from
//!   the bundled `cost_distance`).
//! * **Convenience** — a single friction raster `cost` plus two source rasters
//!   `source1`/`source2`; the tool accumulates cost from each source itself
//!   (8-connected Dijkstra, move cost = distance × mean friction of the two
//!   cells) and sums the results.
//!
//! Without a threshold the output is the summed corridor surface; with
//! `threshold` or `percent` it is that surface masked to the near-optimal band
//! (no-data outside). No-data propagates: a cell is no-data if either input is.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::common::{
    band_to_vec, load_input_raster, parse_optional_output, raster_like_with_data,
    write_or_store_output,
};

pub struct CorridorTool;

impl Tool for CorridorTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "corridor",
            display_name: "Corridor",
            summary: "Least-cost corridor between two accumulated-cost surfaces (ArcGIS Corridor): sum two cost-distance rasters into a corridor surface whose minimum is the least-cost path, then optionally threshold to the near-optimal band. Direct (two cost rasters) or convenience mode (one friction raster + two source rasters).",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "cost1",
                    description: "First accumulated-cost raster (from source A). Direct mode.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cost2",
                    description: "Second accumulated-cost raster (from source B). Direct mode.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cost",
                    description: "Friction/cost-of-passage raster. Convenience mode (with source1/source2).",
                    required: false,
                },
                ToolParamSpec {
                    name: "source1",
                    description: "First source raster (positive cells are sources). Convenience mode.",
                    required: false,
                },
                ToolParamSpec {
                    name: "source2",
                    description: "Second source raster (positive cells are sources). Convenience mode.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output corridor raster path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "threshold",
                    description: "Keep only cells whose corridor cost is at or below this absolute value (masks the rest to no-data).",
                    required: false,
                },
                ToolParamSpec {
                    name: "percent",
                    description: "Keep only cells within this percent above the minimum corridor cost (e.g. 5 = within 5%). Ignored if 'threshold' is set.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read from the input raster(s). Default 1.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let prm = parse_params(args)?;
        let output = parse_optional_output(args, "output")?;

        // ── Obtain the two accumulated-cost surfaces ──────────────────────────
        let (acc1, acc2, template) = match &prm.mode {
            Mode::Direct { cost1, cost2 } => {
                ctx.progress.info("reading accumulated-cost rasters");
                let r1 = load_input_raster(cost1)?;
                let r2 = load_input_raster(cost2)?;
                check_same_grid(&r1, &r2)?;
                let a1 = to_option_vec(&r1, prm.band);
                let a2 = to_option_vec(&r2, prm.band);
                (a1, a2, r1)
            }
            Mode::Convenience {
                cost,
                source1,
                source2,
            } => {
                ctx.progress.info("reading friction and source rasters");
                let friction = load_input_raster(cost)?;
                let s1 = load_input_raster(source1)?;
                let s2 = load_input_raster(source2)?;
                check_same_grid(&friction, &s1)?;
                check_same_grid(&friction, &s2)?;
                ctx.progress
                    .info("accumulating cost distance from each source");
                let a1 = accumulate_cost(&friction, &s1, prm.band);
                let a2 = accumulate_cost(&friction, &s2, prm.band);
                (a1, a2, friction)
            }
        };

        // ── Sum into the corridor surface (no-data if either is missing) ──────
        let n = acc1.len();
        let mut corridor = vec![None; n];
        for i in 0..n {
            if let (Some(a), Some(b)) = (acc1[i], acc2[i]) {
                corridor[i] = Some(a + b);
            }
        }
        let min_val = corridor
            .iter()
            .filter_map(|v| *v)
            .fold(f64::INFINITY, f64::min);
        if !min_val.is_finite() {
            return Err(ToolError::Execution(
                "corridor surface is empty (the two cost surfaces never overlap)".to_string(),
            ));
        }

        // ── Optional thresholding to the near-optimal band ────────────────────
        let cutoff = match (prm.threshold, prm.percent) {
            (Some(t), _) => Some(t),
            (None, Some(p)) => Some(min_val * (1.0 + p / 100.0)),
            (None, None) => None,
        };
        let nodata = template.nodata;
        let mut data = vec![nodata; n];
        let mut band_cells = 0usize;
        let mut valid_cells = 0usize;
        for i in 0..n {
            let Some(v) = corridor[i] else { continue };
            valid_cells += 1;
            match cutoff {
                Some(c) => {
                    if v <= c {
                        data[i] = v;
                        band_cells += 1;
                    }
                }
                None => {
                    data[i] = v;
                    band_cells += 1;
                }
            }
        }

        ctx.progress.info(&format!(
            "corridor min cost {min_val:.3}; {band_cells}/{valid_cells} cell(s) in the output band"
        ));

        let out_raster = raster_like_with_data(&template, data, nodata, DataType::F32)?;
        let out_path = write_or_store_output(out_raster, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("min_cost".to_string(), json!(min_val));
        outputs.insert("band_cells".to_string(), json!(band_cells));
        outputs.insert("valid_cells".to_string(), json!(valid_cells));
        if let Some(c) = cutoff {
            outputs.insert("cutoff".to_string(), json!(c));
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── Accumulated cost (convenience mode) ──────────────────────────────────────

/// Node in the Dijkstra frontier; ordered so the `BinaryHeap` (a max-heap) pops
/// the smallest accumulated cost first.
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
        // Reversed so the min-cost node is "greatest" for the max-heap.
        other
            .cost
            .partial_cmp(&self.cost)
            .unwrap_or(Ordering::Equal)
            .then(self.idx.cmp(&other.idx))
    }
}

/// Accumulated least-cost distance from every source cell over a friction
/// raster: 8-connected Dijkstra where the step cost is the geometric distance
/// times the mean friction of the two cells. Source cells (positive, non-nodata
/// in `source`) start at 0. Returns `None` for cells unreachable or on nodata
/// friction.
fn accumulate_cost(friction: &Raster, source: &Raster, band: isize) -> Vec<Option<f64>> {
    let rows = friction.rows;
    let cols = friction.cols;
    let n = rows * cols;
    let fr = band_to_vec_opt(friction, band);
    let src = band_to_vec(source, band);
    let s_nodata = source.nodata;

    let cx = friction.cell_size_x;
    let cy = friction.cell_size_y;
    let d_orth = cx.min(cy);
    let d_diag = (cx * cx + cy * cy).sqrt();

    let mut dist = vec![f64::INFINITY; n];
    let mut done = vec![false; n];
    let mut heap: BinaryHeap<Node> = BinaryHeap::new();
    for i in 0..n {
        let s = src[i];
        if s != s_nodata && s.is_finite() && s > 0.0 && fr[i].is_some() {
            dist[i] = 0.0;
            heap.push(Node { cost: 0.0, idx: i });
        }
    }

    // 8-neighbour offsets with their base geometric distance.
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

    while let Some(Node { cost, idx }) = heap.pop() {
        if done[idx] {
            continue;
        }
        done[idx] = true;
        let r = (idx / cols) as isize;
        let c = (idx % cols) as isize;
        let Some(f_here) = fr[idx] else { continue };
        for (dr, dc, base) in neigh {
            let nr = r + dr;
            let nc = c + dc;
            if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                continue;
            }
            let nidx = nr as usize * cols + nc as usize;
            let Some(f_n) = fr[nidx] else { continue };
            let step = base * 0.5 * (f_here + f_n);
            let nd = cost + step;
            if nd < dist[nidx] {
                dist[nidx] = nd;
                heap.push(Node {
                    cost: nd,
                    idx: nidx,
                });
            }
        }
    }

    dist.into_iter()
        .map(|d| if d.is_finite() { Some(d) } else { None })
        .collect()
}

// ── Raster helpers ───────────────────────────────────────────────────────────

/// Reads a band into `Option<f64>` cells (no-data / non-finite -> `None`).
fn to_option_vec(raster: &Raster, band: isize) -> Vec<Option<f64>> {
    band_to_vec_opt(raster, band)
}

fn band_to_vec_opt(raster: &Raster, band: isize) -> Vec<Option<f64>> {
    let raw = band_to_vec_band(raster, band);
    raw.into_iter()
        .map(|v| {
            if v == raster.nodata || !v.is_finite() {
                None
            } else {
                Some(v)
            }
        })
        .collect()
}

/// Like `common::band_to_vec` but for an explicit band index.
fn band_to_vec_band(raster: &Raster, band: isize) -> Vec<f64> {
    let rows = raster.rows;
    let cols = raster.cols;
    let mut out = vec![0.0_f64; rows * cols];
    for row in 0..rows {
        for col in 0..cols {
            out[row * cols + col] = raster.get(band, row as isize, col as isize);
        }
    }
    out
}

fn check_same_grid(a: &Raster, b: &Raster) -> Result<(), ToolError> {
    if a.rows != b.rows || a.cols != b.cols {
        return Err(ToolError::Validation(format!(
            "raster dimensions differ: {}x{} vs {}x{}",
            a.rows, a.cols, b.rows, b.cols
        )));
    }
    Ok(())
}

// ── Parameters ────────────────────────────────────────────────────────────────

enum Mode {
    Direct {
        cost1: String,
        cost2: String,
    },
    Convenience {
        cost: String,
        source1: String,
        source2: String,
    },
}

struct Params {
    mode: Mode,
    threshold: Option<f64>,
    percent: Option<f64>,
    band: isize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let cost1 = opt_str(args, "cost1");
    let cost2 = opt_str(args, "cost2");
    let cost = opt_str(args, "cost");
    let source1 = opt_str(args, "source1");
    let source2 = opt_str(args, "source2");

    let mode = if let (Some(cost1), Some(cost2)) = (&cost1, &cost2) {
        Mode::Direct {
            cost1: cost1.clone(),
            cost2: cost2.clone(),
        }
    } else if let (Some(cost), Some(source1), Some(source2)) = (&cost, &source1, &source2) {
        Mode::Convenience {
            cost: cost.clone(),
            source1: source1.clone(),
            source2: source2.clone(),
        }
    } else {
        return Err(ToolError::Validation(
            "provide either 'cost1'+'cost2' (direct) or 'cost'+'source1'+'source2' (convenience)"
                .to_string(),
        ));
    };

    let threshold = opt_f64(args, "threshold")?;
    let percent = opt_f64(args, "percent")?;
    if let Some(p) = percent {
        if p < 0.0 {
            return Err(ToolError::Validation(
                "'percent' must be non-negative".to_string(),
            ));
        }
    }
    let band_1based = opt_f64(args, "band")?.map(|v| v as i64).unwrap_or(1);
    if band_1based < 1 {
        return Err(ToolError::Validation("'band' must be >= 1".to_string()));
    }
    Ok(Params {
        mode,
        threshold,
        percent,
        band: (band_1based - 1) as isize,
    })
}

fn opt_str(args: &ToolArgs, key: &str) -> Option<String> {
    match args.get(key) {
        Some(Value::String(s)) if !s.trim().is_empty() => Some(s.trim().to_string()),
        _ => None,
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn raster_from(cols: usize, rows: usize, data: Vec<f64>, nodata: f64) -> String {
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
                r.set(0, row as isize, col as isize, data[row * cols + col])
                    .unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CorridorTool.run(&args, &ctx()).unwrap();
        let r = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, r)
    }

    /// Direct mode: corridor = cost1 + cost2, min at the meeting point.
    #[test]
    fn sums_two_cost_surfaces() {
        // A 1D corridor: distance from left (0..4) and from right (4..0).
        let cost1 = raster_from(5, 1, vec![0.0, 1.0, 2.0, 3.0, 4.0], -1.0);
        let cost2 = raster_from(5, 1, vec![4.0, 3.0, 2.0, 1.0, 0.0], -1.0);
        let (out, r) = run(json!({ "cost1": cost1, "cost2": cost2 }));
        // Every cell sums to 4 (the least-cost path length); min = 4.
        assert_eq!(out.outputs["min_cost"], json!(4.0));
        for c in 0..5 {
            assert_eq!(r.get(0, 0, c), 4.0);
        }
    }

    /// Percent threshold masks cells above (1+p%) * min.
    #[test]
    fn percent_threshold_masks_band() {
        // cost1+cost2: a valley at the centre, higher at the ends.
        let cost1 = raster_from(5, 1, vec![0.0, 1.0, 2.0, 3.0, 4.0], -1.0);
        let cost2 = raster_from(5, 1, vec![8.0, 5.0, 2.0, 5.0, 8.0], -1.0);
        // sums: 8, 6, 4, 8, 12 -> min 4
        let (out, r) = run(json!({ "cost1": cost1, "cost2": cost2, "percent": 60.0 }));
        assert_eq!(out.outputs["min_cost"], json!(4.0));
        // cutoff = 4 * 1.6 = 6.4 -> keep cells with sum <= 6.4: indices 1 (6) and 2 (4).
        assert_eq!(out.outputs["band_cells"], json!(2));
        assert_eq!(r.get(0, 0, 2), 4.0);
        assert_eq!(r.get(0, 0, 1), 6.0);
        assert_eq!(r.get(0, 0, 0), r.nodata);
    }

    /// No-data propagates: a cell missing in either input is missing in the sum.
    #[test]
    fn nodata_propagates() {
        let cost1 = raster_from(3, 1, vec![0.0, -1.0, 2.0], -1.0);
        let cost2 = raster_from(3, 1, vec![2.0, 1.0, 0.0], -1.0);
        let (out, r) = run(json!({ "cost1": cost1, "cost2": cost2 }));
        assert_eq!(r.get(0, 0, 1), r.nodata, "nodata should propagate");
        assert_eq!(out.outputs["valid_cells"], json!(2));
    }

    /// Convenience mode: uniform friction gives a distance-based corridor whose
    /// minimum lies on the straight line between the two sources.
    #[test]
    fn convenience_mode_accumulates_cost() {
        // 5x1 uniform friction 1; source1 at left, source2 at right.
        let cost = raster_from(5, 1, vec![1.0; 5], -1.0);
        let s1 = raster_from(5, 1, vec![1.0, 0.0, 0.0, 0.0, 0.0], 0.0);
        let s2 = raster_from(5, 1, vec![0.0, 0.0, 0.0, 0.0, 1.0], 0.0);
        let (out, r) = run(json!({ "cost": cost, "source1": s1, "source2": s2 }));
        // Distance from each end summed is constant (= 4) along the whole line.
        let min = out.outputs["min_cost"].as_f64().unwrap();
        assert!((min - 4.0).abs() < 1e-9, "min cost {min} should be 4");
        for c in 0..5 {
            assert!((r.get(0, 0, c) - 4.0).abs() < 1e-9);
        }
    }

    #[test]
    fn rejects_incomplete_inputs() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CorridorTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "cost1": "a.tif" })).is_err());
        assert!(bad(json!({ "cost": "f.tif", "source1": "s.tif" })).is_err());
        assert!(bad(json!({ "cost1": "a.tif", "cost2": "b.tif" })).is_ok());
    }
}
