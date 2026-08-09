//! GeoLibre tool: the backlink direction raster of an accumulative cost surface.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Cost Back Link* and *Path Back Link*
//! (Spatial Analyst Distance toolset).
//!
//! ## Why the intermediate matters
//!
//! `cost_distance` produces accumulated cost and `cost_pathway`,
//! `cost_allocation` and `path_distance` consume it, but the **backlink raster
//! itself is never exposed**. That raster is the reusable half of the solve:
//! with it, retracing an arbitrary number of destinations costs one grid walk
//! each against a *single* accumulation pass. Without it, every destination
//! forces a fresh cost solve, making a many-destination job O(destinations)
//! full solves instead of one.
//!
//! It is also a documented input to ArcGIS's own path tools, so its absence
//! blocks porting existing workflows.
//!
//! ## Encoding
//!
//! ArcGIS's convention, preserved exactly so the output is interchangeable:
//! `0` marks a source cell, and `1`–`8` give the direction of the **next cell
//! on the way back to the source**, starting at `1` = east and proceeding
//! counter-clockwise. Unreachable cells are no-data.
//!
//! ## Surface mode
//!
//! Supplying `surface` enables the Path Back Link variant: horizontal distance
//! between cell centres is replaced by true surface distance
//! `sqrt(dxy^2 + dz^2)`, so a route climbing a slope pays for the climb. This
//! is the one piece of Path Distance's factor machinery that needs no extra
//! rasters, and it is what makes the tool useful on terrain.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap};

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};
use wbvector::Geometry;

use crate::args_common::{opt_positive_f64, req_str};
use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};
use crate::raster_stack::check_alignment_refs;
use crate::vector_common::{load_input_layer, parse_optional_str};

/// Neighbour offsets in ArcGIS backlink order: 1 = east, counter-clockwise.
/// Index `k` here encodes as `k + 1`.
const NEIGHBOURS: [(isize, isize); 8] = [
    (0, 1),   // 1 east
    (-1, 1),  // 2 north-east
    (-1, 0),  // 3 north
    (-1, -1), // 4 north-west
    (0, -1),  // 5 west
    (1, -1),  // 6 south-west
    (1, 0),   // 7 south
    (1, 1),   // 8 south-east
];

pub struct CostBackLinkTool;

impl Tool for CostBackLinkTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "cost_back_link",
            display_name: "Cost Back Link",
            summary: "Emits the 0-8 neighbour-direction raster encoding, for every cell, the next step along the least-cost route back to the nearest source, plus the accumulated-cost surface (ArcGIS Cost Back Link and Path Back Link). cost_distance, cost_pathway and cost_allocation all compute this internally and discard it, so retracing N destinations currently costs N full cost solves instead of one accumulation pass plus N cheap walks.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "source",
                    description: "Source cells: a raster (any non-zero, non-no-data cell is a source) or a point/multipoint vector layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "cost",
                    description: "Cost-per-unit-distance raster. Optional when 'surface' is given (then cost defaults to 1 everywhere, giving pure surface distance).",
                    required: false,
                },
                ToolParamSpec {
                    name: "surface",
                    description: "Optional elevation raster enabling the Path Back Link variant: horizontal spacing is replaced by true surface distance sqrt(dxy^2 + dz^2).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output backlink raster (0 = source, 1-8 = direction toward the source, 1 = east counter-clockwise). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "out_distance",
                    description: "Accumulated-cost raster. Always produced; falls back to an in-memory handle when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_distance",
                    description: "Optional accumulation cutoff: cells costing more than this are left unreachable.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "source")?;
        opt_positive_f64(args, "max_distance")?;
        if args.get("cost").is_none() && args.get("surface").is_none() {
            return Err(ToolError::Validation(
                "supply 'cost', 'surface', or both: with neither there is nothing to accumulate"
                    .to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let source_spec = req_str(args, "source")?.to_string();
        let cost_path = parse_optional_str(args, "cost")?.map(str::to_string);
        let surface_path = parse_optional_str(args, "surface")?.map(str::to_string);
        let output = parse_optional_output(args, "output")?;
        let out_distance = parse_optional_output(args, "out_distance")?;
        let max_distance = opt_positive_f64(args, "max_distance")?.unwrap_or(f64::INFINITY);

        // The grid comes from whichever cost-like raster was supplied.
        let cost = match &cost_path {
            Some(p) => Some(load_input_raster(p)?),
            None => None,
        };
        let surface = match &surface_path {
            Some(p) => Some(load_input_raster(p)?),
            None => None,
        };
        let template = cost
            .as_ref()
            .or(surface.as_ref())
            .ok_or_else(|| ToolError::Validation("supply 'cost' or 'surface'".to_string()))?;
        if let (Some(c), Some(s)) = (&cost, &surface) {
            check_alignment_refs(&[c, s])?;
        }

        let rows = template.rows;
        let cols = template.cols;
        let (dx, dy) = (template.cell_size_x, template.cell_size_y);

        let sources = load_sources(&source_spec, template)?;
        if sources.is_empty() {
            return Err(ToolError::Execution(
                "'source' selected no cells on the cost grid".to_string(),
            ));
        }
        ctx.progress
            .info(&format!("{rows}x{cols}, {} source cell(s)", sources.len()));

        let nodata = -9999.0_f64;
        let mut accum = vec![f64::INFINITY; rows * cols];
        let mut backlink = vec![nodata; rows * cols];
        let mut heap: BinaryHeap<Step> = BinaryHeap::new();

        for &idx in &sources {
            accum[idx] = 0.0;
            backlink[idx] = 0.0; // ArcGIS marks source cells 0.
            heap.push(Step { cost: 0.0, idx });
        }

        let cell_cost = |idx: usize| -> Option<f64> {
            match &cost {
                None => Some(1.0),
                Some(c) => {
                    let v = c.get(0, (idx / cols) as isize, (idx % cols) as isize);
                    (v != c.nodata && v.is_finite() && v >= 0.0).then_some(v)
                }
            }
        };
        let elev = |idx: usize| -> Option<f64> {
            let s = surface.as_ref()?;
            let v = s.get(0, (idx / cols) as isize, (idx % cols) as isize);
            (v != s.nodata && v.is_finite()).then_some(v)
        };

        let mut settled = 0_u64;
        while let Some(Step { cost: c, idx }) = heap.pop() {
            // Stale heap entry from a later relaxation of the same cell.
            if c > accum[idx] {
                continue;
            }
            settled += 1;
            let (r, col) = ((idx / cols) as isize, (idx % cols) as isize);
            let from_cost = match cell_cost(idx) {
                Some(v) => v,
                None => continue,
            };

            for (k, (dr, dc)) in NEIGHBOURS.iter().enumerate() {
                let (nr, nc) = (r + dr, col + dc);
                if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                    continue;
                }
                let nidx = nr as usize * cols + nc as usize;
                let Some(to_cost) = cell_cost(nidx) else {
                    continue;
                };

                let mut span = ((*dr as f64 * dy).powi(2) + (*dc as f64 * dx).powi(2)).sqrt();
                if surface.is_some() {
                    // Path Back Link: pay for the climb, not just the plan
                    // distance. Missing elevation leaves the pair unusable.
                    let (Some(za), Some(zb)) = (elev(idx), elev(nidx)) else {
                        continue;
                    };
                    span = (span * span + (zb - za).powi(2)).sqrt();
                }
                // The standard trapezoidal accumulation: the mean of the two
                // cells' cost rates over the distance between their centres.
                let step = span * (from_cost + to_cost) / 2.0;
                let candidate = c + step;
                if candidate > max_distance || candidate >= accum[nidx] {
                    continue;
                }
                accum[nidx] = candidate;
                // The neighbour steps back toward `idx`, which is the OPPOSITE
                // of the direction we just travelled. NEIGHBOURS is arranged so
                // that opposite directions are 4 apart.
                backlink[nidx] = (((k + 4) % 8) + 1) as f64;
                heap.push(Step {
                    cost: candidate,
                    idx: nidx,
                });
            }
        }

        let mut reachable = 0_u64;
        let mut distance = vec![nodata; rows * cols];
        for i in 0..rows * cols {
            if accum[i].is_finite() {
                distance[i] = accum[i];
                reachable += 1;
            } else {
                backlink[i] = nodata;
            }
        }

        let link_raster = raster_like_with_data(template, backlink, nodata, DataType::I32)?;
        let link_path = write_or_store_output(link_raster, output)?;
        // Repo convention (create_overpass): emit the secondary output
        // unconditionally so a caller with no scratch path still gets it back.
        let dist_raster = raster_like_with_data(template, distance, nodata, DataType::F32)?;
        let dist_path = write_or_store_output(dist_raster, out_distance)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(link_path));
        outputs.insert("out_distance".to_string(), json!(dist_path));
        outputs.insert("source_cells".to_string(), json!(sources.len()));
        outputs.insert("reachable_cells".to_string(), json!(reachable));
        outputs.insert("settled_cells".to_string(), json!(settled));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        Ok(ToolRunResult { outputs })
    }
}

/// Min-heap entry. `BinaryHeap` is a max-heap, so the ordering is reversed.
struct Step {
    cost: f64,
    idx: usize,
}

impl PartialEq for Step {
    fn eq(&self, other: &Self) -> bool {
        self.cost == other.cost && self.idx == other.idx
    }
}
impl Eq for Step {}
impl Ord for Step {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .cost
            .total_cmp(&self.cost)
            .then_with(|| other.idx.cmp(&self.idx))
    }
}
impl PartialOrd for Step {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Resolves the source parameter to cell indices on the cost grid.
///
/// A raster source uses any non-zero, non-no-data cell (the ArcGIS rule); a
/// vector source snaps each point to the cell containing it.
fn load_sources(spec: &str, template: &Raster) -> Result<Vec<usize>, ToolError> {
    let cols = template.cols;
    let rows = template.rows;

    // Probe as a raster first. A raster that loads but fails alignment is a
    // real error and must surface as one — falling through to the vector path
    // would replace it with a misleading "not a vector layer" message.
    if let Ok(raster) = load_input_raster(spec) {
        check_alignment_refs(&[template, &raster])?;
        let mut out = Vec::new();
        for r in 0..rows {
            for c in 0..cols {
                let v = raster.get(0, r as isize, c as isize);
                if v != raster.nodata && v.is_finite() && v != 0.0 {
                    out.push(r * cols + c);
                }
            }
        }
        return Ok(out);
    }

    let layer = load_input_layer(spec)?;
    let mut out = Vec::new();
    let push = |x: f64, y: f64, out: &mut Vec<usize>| {
        // y is measured up from y_min, while rows count down from the top.
        let c = ((x - template.x_min) / template.cell_size_x).floor();
        let r = (rows as f64 - 1.0) - ((y - template.y_min) / template.cell_size_y).floor();
        if c >= 0.0 && r >= 0.0 && (c as usize) < cols && (r as usize) < rows {
            out.push(r as usize * cols + c as usize);
        }
    };
    for feature in layer.iter() {
        match feature.geometry.as_ref() {
            Some(Geometry::Point(p)) => push(p.x, p.y, &mut out),
            Some(Geometry::MultiPoint(ps)) => {
                for p in ps {
                    push(p.x, p.y, &mut out);
                }
            }
            _ => {}
        }
    }
    out.sort_unstable();
    out.dedup();
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
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

    fn raster(rows: usize, cols: usize, data: &[f64]) -> String {
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

    fn run(args: Value) -> (Raster, Raster, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = CostBackLinkTool.run(&args, &ctx()).unwrap();
        let link = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        let dist = load_input_raster(res.outputs["out_distance"].as_str().unwrap()).unwrap();
        (link, dist, res)
    }

    #[test]
    fn the_source_cell_is_marked_zero_and_costs_nothing() {
        let src = raster(1, 3, &[1.0, 0.0, 0.0]);
        let cost = raster(1, 3, &[1.0, 1.0, 1.0]);
        let (link, dist, _) = run(json!({"source": src, "cost": cost}));
        assert_eq!(link.get(0, 0, 0), 0.0);
        assert_eq!(dist.get(0, 0, 0), 0.0);
    }

    #[test]
    fn backlinks_point_back_toward_the_source() {
        // Source at the west end of a 1x4 row; every other cell must step
        // WEST (code 5) to get home.
        let src = raster(1, 4, &[1.0, 0.0, 0.0, 0.0]);
        let cost = raster(1, 4, &[1.0, 1.0, 1.0, 1.0]);
        let (link, _, _) = run(json!({"source": src, "cost": cost}));
        for c in 1..4 {
            assert_eq!(link.get(0, 0, c), 5.0, "cell {c} should point west");
        }
    }

    #[test]
    fn a_source_at_the_east_end_flips_every_backlink_to_east() {
        // The mirror case. If the opposite-direction arithmetic were wrong,
        // one of these two tests would pass and the other fail.
        let src = raster(1, 4, &[0.0, 0.0, 0.0, 1.0]);
        let cost = raster(1, 4, &[1.0, 1.0, 1.0, 1.0]);
        let (link, _, _) = run(json!({"source": src, "cost": cost}));
        for c in 0..3 {
            assert_eq!(link.get(0, 0, c), 1.0, "cell {c} should point east");
        }
    }

    #[test]
    fn accumulated_cost_uses_the_mean_of_the_two_cell_rates() {
        // Uniform cost 1 over unit cells: each step costs 1.
        let src = raster(1, 4, &[1.0, 0.0, 0.0, 0.0]);
        let cost = raster(1, 4, &[1.0, 1.0, 1.0, 1.0]);
        let (_, dist, _) = run(json!({"source": src, "cost": cost}));
        for c in 0..4 {
            assert!(
                (dist.get(0, 0, c) - c as f64).abs() < 1e-5,
                "cell {c} cost {}",
                dist.get(0, 0, c)
            );
        }
    }

    #[test]
    fn a_high_cost_barrier_reroutes_the_backlink() {
        // 3x3, source at (0,0). The middle column is very expensive, so the
        // cheapest way to (1,2) is around rather than straight through.
        let src = raster(3, 3, &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let cost = raster(3, 3, &[1.0, 1000.0, 1.0, 1.0, 1000.0, 1.0, 1.0, 1.0, 1.0]);
        let (_, dist, _) = run(json!({"source": src, "cost": cost}));
        // Going around the bottom must beat crossing the barrier.
        assert!(
            dist.get(0, 1, 2) < 500.0,
            "route did not avoid the barrier: {}",
            dist.get(0, 1, 2)
        );
    }

    #[test]
    fn max_distance_leaves_far_cells_unreachable() {
        let src = raster(1, 6, &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let cost = raster(1, 6, &[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
        let (link, dist, res) = run(json!({
            "source": src, "cost": cost, "max_distance": 2.5,
        }));
        assert!(dist.get(0, 0, 2).is_finite() && dist.get(0, 0, 2) != dist.nodata);
        assert_eq!(dist.get(0, 0, 5), dist.nodata);
        assert_eq!(link.get(0, 0, 5), link.nodata);
        assert_eq!(res.outputs["reachable_cells"], json!(3));
    }

    #[test]
    fn surface_mode_charges_for_the_climb() {
        // Flat vs a 1-unit rise per cell. With cost 1 everywhere the surface
        // run must be strictly longer: sqrt(2) per step rather than 1.
        let src = raster(1, 3, &[1.0, 0.0, 0.0]);
        let cost = raster(1, 3, &[1.0, 1.0, 1.0]);
        let flat = raster(1, 3, &[0.0, 0.0, 0.0]);
        let ramp = raster(1, 3, &[0.0, 1.0, 2.0]);
        let (_, d_flat, _) = run(json!({"source": src, "cost": cost, "surface": flat}));
        let (_, d_ramp, _) = run(json!({"source": src, "cost": cost, "surface": ramp}));
        assert!((d_flat.get(0, 0, 2) - 2.0).abs() < 1e-5);
        let expect = 2.0 * 2.0_f64.sqrt();
        assert!(
            (d_ramp.get(0, 0, 2) - expect).abs() < 1e-4,
            "expected {expect}, got {}",
            d_ramp.get(0, 0, 2)
        );
    }

    #[test]
    fn out_distance_is_produced_even_without_a_path() {
        // The round-16 lesson: a secondary output gated on a supplied path
        // silently vanishes for callers with no scratch directory.
        let src = raster(1, 2, &[1.0, 0.0]);
        let cost = raster(1, 2, &[1.0, 1.0]);
        let args: ToolArgs = serde_json::from_value(json!({"source": src, "cost": cost})).unwrap();
        let res = CostBackLinkTool.run(&args, &ctx()).unwrap();
        let p = res.outputs["out_distance"].as_str().unwrap();
        assert!(!p.is_empty());
        assert!(load_input_raster(p).is_ok());
    }

    #[test]
    fn diagonal_steps_use_the_diagonal_codes() {
        // Covers codes 2/4/6/8, which the axis-aligned tests never exercise.
        // Source at the north-west corner of a 2x2 grid: the opposite corner
        // steps back north-west, which is code 4.
        let src = raster(2, 2, &[1.0, 0.0, 0.0, 0.0]);
        let cost = raster(2, 2, &[1.0, 1.0, 1.0, 1.0]);
        let (link, _, _) = run(json!({"source": src, "cost": cost}));
        assert_eq!(link.get(0, 1, 1), 4.0, "expected a north-west backlink");
    }

    #[test]
    fn rejects_bad_parameters() {
        let src = raster(1, 2, &[1.0, 0.0]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CostBackLinkTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        // Neither cost nor surface: nothing to accumulate.
        assert!(bad(json!({"source": src})));
        assert!(bad(json!({"source": src, "cost": src, "max_distance": -1})));
    }
}
