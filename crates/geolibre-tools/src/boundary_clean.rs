//! GeoLibre tool: smooth and generalize the boundaries between zones of a
//! categorical raster.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Boundary Clean* and *Majority Filter*
//! (Spatial Analyst). The bundled whitebox suite has `clump`, `nibble`, and the
//! GeoLibre `expand_shrink`, but nothing that smooths ragged zone boundaries or
//! removes single-cell speckle from a classified raster — the standard cleanup
//! step before `raster_to_vector_polygons` →
//! `regularize_building_footprints` / `smooth_natural_features`.
//!
//! Two modes:
//!
//! * **majority** — each cell is replaced by the most frequent value among its
//!   4- or 8-connected neighbours, provided that value reaches a threshold
//!   (`majority` = a strict majority of the neighbours, `half` = at least half).
//!   Removes isolated speckle while leaving coherent regions intact. Runs
//!   `iterations` times.
//! * **expand_shrink** (Boundary Clean) — smooths zone boundaries by an
//!   expansion pass (higher-priority zones grow one cell into their neighbours)
//!   followed by a shrink pass (they retreat one cell). Zone priority is set by
//!   total cell count: `descending` lets larger zones dominate, `ascending`
//!   favours smaller zones, `none` breaks ties by class value only. The
//!   expand-then-shrink round-trip removes jagged one-cell protrusions and
//!   inlets while roughly conserving each zone's area.
//!
//! No-data cells are barriers: never overwritten, never used as replacement
//! values. Ties are broken toward the smaller class value so results are
//! deterministic.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{load_input_raster, parse_optional_output, raster_like_with_data};

pub struct BoundaryCleanTool;

impl Tool for BoundaryCleanTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "boundary_clean",
            display_name: "Boundary Clean",
            summary: "Smooth categorical raster zone boundaries and remove speckle via majority filtering or expand/shrink boundary cleaning, the standard cleanup before polygonization — like ArcGIS Boundary Clean / Majority Filter.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input categorical (integer-valued) raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output raster path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'majority' (neighbour-majority speckle removal; default) or 'expand_shrink' (Boundary Clean smoothing).",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighbors",
                    description: "Neighbourhood connectivity: 4 (rook) or 8 (queen; default).",
                    required: false,
                },
                ToolParamSpec {
                    name: "threshold",
                    description: "majority mode only: 'majority' (strict majority of neighbours; default) or 'half' (at least half).",
                    required: false,
                },
                ToolParamSpec {
                    name: "iterations",
                    description: "Number of passes to run (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "sort",
                    description: "expand_shrink mode only: zone priority by cell count — 'descending' (larger zones win; default), 'ascending' (smaller zones win), or 'none'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read (default 1).",
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

        let raster = load_input_raster(input)?;
        if prm.band < 0 || prm.band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "band {} out of range (raster has {} band(s))",
                prm.band + 1,
                raster.bands
            )));
        }
        let rows = raster.rows;
        let cols = raster.cols;
        let nodata = raster.nodata;

        let mut grid = vec![0.0f64; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                grid[r * cols + c] = raster.get(prm.band, r as isize, c as isize);
            }
        }
        let is_valid = |v: f64| v != nodata && !v.is_nan();

        ctx.progress.info(&format!(
            "boundary clean: {} mode, {} pass(es)",
            match prm.method {
                Method::Majority => "majority",
                Method::ExpandShrink => "expand_shrink",
            },
            prm.iterations
        ));

        let neigh: &[(isize, isize)] = if prm.neighbors == 4 { &NEIGH4 } else { &NEIGH8 };
        let mut changed = 0usize;

        match prm.method {
            Method::Majority => {
                for _ in 0..prm.iterations {
                    let snapshot = grid.clone();
                    let mut n = 0usize;
                    for r in 0..rows {
                        for c in 0..cols {
                            let idx = r * cols + c;
                            let v = snapshot[idx];
                            if !is_valid(v) {
                                continue;
                            }
                            if let Some(nv) = majority_replacement(
                                &snapshot,
                                rows,
                                cols,
                                r,
                                c,
                                v,
                                prm.min_count,
                                neigh,
                                &is_valid,
                            ) {
                                grid[idx] = nv;
                                n += 1;
                            }
                        }
                    }
                    changed += n;
                    if n == 0 {
                        break;
                    }
                }
            }
            Method::ExpandShrink => {
                // Zone priority from the original cell counts.
                let priority = zone_priority(&grid, prm.sort, &is_valid);
                for _ in 0..prm.iterations {
                    let before = grid.clone();
                    // Expansion: a cell is captured by the highest-priority
                    // neighbouring zone that outranks its own zone.
                    let n1 = morph_pass(&mut grid, rows, cols, neigh, &is_valid, &priority, true);
                    // Shrink: a cell retreats to the lowest-priority neighbouring
                    // zone that its own zone outranks.
                    let n2 = morph_pass(&mut grid, rows, cols, neigh, &is_valid, &priority, false);
                    changed += n1 + n2;
                    if grid == before {
                        break;
                    }
                }
            }
        }

        let out_raster = raster_like_with_data(&raster, grid, nodata, DataType::F32)?;
        let out_path = crate::common::write_or_store_output(out_raster, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("cells_changed".to_string(), json!(changed));
        Ok(ToolRunResult { outputs })
    }
}

const NEIGH4: [(isize, isize); 4] = [(-1, 0), (0, -1), (0, 1), (1, 0)];
const NEIGH8: [(isize, isize); 8] = [
    (-1, -1),
    (-1, 0),
    (-1, 1),
    (0, -1),
    (0, 1),
    (1, -1),
    (1, 0),
    (1, 1),
];

/// Returns the majority neighbour value for a cell if it reaches `min_count` and
/// differs from the current value, else `None`.
#[allow(clippy::too_many_arguments)]
fn majority_replacement(
    snap: &[f64],
    rows: usize,
    cols: usize,
    r: usize,
    c: usize,
    center: f64,
    min_count: usize,
    neigh: &[(isize, isize)],
    is_valid: &dyn Fn(f64) -> bool,
) -> Option<f64> {
    let mut tally: BTreeMap<u64, (f64, usize)> = BTreeMap::new();
    for (dr, dc) in neigh {
        let nr = r as isize + dr;
        let nc = c as isize + dc;
        if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
            continue;
        }
        let nv = snap[nr as usize * cols + nc as usize];
        if !is_valid(nv) {
            continue;
        }
        let e = tally.entry(nv.to_bits()).or_insert((nv, 0));
        e.1 += 1;
    }
    // Most frequent neighbour value; ties -> smaller class value.
    let (best_v, best_n) = tally
        .values()
        .copied()
        .max_by(|a, b| a.1.cmp(&b.1).then(b.0.total_cmp(&a.0)))?;
    if best_n >= min_count && best_v != center {
        Some(best_v)
    } else {
        None
    }
}

/// One morphology pass over the grid. When `expand` is true, each cell adopts the
/// highest-priority neighbouring zone that ranks *above* its own zone; when false
/// (shrink) it adopts the lowest-priority neighbour that ranks *below* it.
/// Updates are simultaneous (computed from a snapshot). Returns the change count.
fn morph_pass(
    grid: &mut [f64],
    rows: usize,
    cols: usize,
    neigh: &[(isize, isize)],
    is_valid: &dyn Fn(f64) -> bool,
    priority: &BTreeMap<u64, f64>,
    expand: bool,
) -> usize {
    let snap = grid.to_vec();
    let rank = |v: f64| -> (f64, f64) {
        // Compare by (priority, then smaller value wins on ties).
        (*priority.get(&v.to_bits()).unwrap_or(&0.0), -v)
    };
    let mut changed = 0usize;
    for r in 0..rows {
        for c in 0..cols {
            let idx = r * cols + c;
            let v = snap[idx];
            if !is_valid(v) {
                continue;
            }
            let self_rank = rank(v);
            let mut best: Option<(f64, (f64, f64))> = None;
            for (dr, dc) in neigh {
                let nr = r as isize + dr;
                let nc = c as isize + dc;
                if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                    continue;
                }
                let nv = snap[nr as usize * cols + nc as usize];
                if !is_valid(nv) || nv == v {
                    continue;
                }
                let nr_rank = rank(nv);
                let qualifies = if expand {
                    nr_rank > self_rank
                } else {
                    nr_rank < self_rank
                };
                if !qualifies {
                    continue;
                }
                let better = match best {
                    None => true,
                    Some((_, br)) => {
                        if expand {
                            nr_rank > br
                        } else {
                            nr_rank < br
                        }
                    }
                };
                if better {
                    best = Some((nv, nr_rank));
                }
            }
            if let Some((nv, _)) = best {
                grid[idx] = nv;
                changed += 1;
            }
        }
    }
    changed
}

/// Priority weight per zone value, derived from cell counts and the sort mode.
fn zone_priority(grid: &[f64], sort: Sort, is_valid: &dyn Fn(f64) -> bool) -> BTreeMap<u64, f64> {
    let mut counts: BTreeMap<u64, (f64, usize)> = BTreeMap::new();
    for &v in grid {
        if is_valid(v) {
            let e = counts.entry(v.to_bits()).or_insert((v, 0));
            e.1 += 1;
        }
    }
    counts
        .into_iter()
        .map(|(bits, (_, n))| {
            let p = match sort {
                Sort::Descending => n as f64,
                Sort::Ascending => -(n as f64),
                Sort::None => 0.0,
            };
            (bits, p)
        })
        .collect()
}

// ── Parameters ─────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Majority,
    ExpandShrink,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Sort {
    Descending,
    Ascending,
    None,
}

struct Params {
    method: Method,
    neighbors: u8,
    /// Minimum neighbour count for a majority replacement.
    min_count: usize,
    iterations: usize,
    sort: Sort,
    band: isize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let method = match args.get("method").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("majority") => Method::Majority,
        Some("expand_shrink") => Method::ExpandShrink,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "'method' must be 'majority' or 'expand_shrink', got '{other}'"
            )))
        }
    };
    let neighbors = match args.get("neighbors") {
        None | Some(Value::Null) => 8u8,
        Some(Value::Number(n)) => match n.as_u64() {
            Some(4) => 4,
            Some(8) => 8,
            _ => return Err(ToolError::Validation("'neighbors' must be 4 or 8".into())),
        },
        Some(Value::String(s)) if s.trim().is_empty() => 8,
        Some(Value::String(s)) => match s.trim() {
            "4" => 4,
            "8" => 8,
            _ => return Err(ToolError::Validation("'neighbors' must be 4 or 8".into())),
        },
        _ => return Err(ToolError::Validation("'neighbors' must be 4 or 8".into())),
    };
    let threshold_half = match args.get("threshold").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("majority") => false,
        Some("half") => true,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "'threshold' must be 'majority' or 'half', got '{other}'"
            )))
        }
    };
    // majority: strictly more than half of the neighbours; half: at least half.
    let min_count = if threshold_half {
        (neighbors as usize).div_ceil(2)
    } else {
        neighbors as usize / 2 + 1
    };
    let iterations = match args.get("iterations") {
        None | Some(Value::Null) => 1usize,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(1).max(1) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'iterations' must be a positive integer".into()))?
            .max(1),
        _ => {
            return Err(ToolError::Validation(
                "'iterations' must be an integer".into(),
            ))
        }
    };
    let sort = match args.get("sort").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("descending") => Sort::Descending,
        Some("ascending") => Sort::Ascending,
        Some("none") => Sort::None,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "'sort' must be 'descending', 'ascending', or 'none', got '{other}'"
            )))
        }
    };
    let band_1based = match args.get("band") {
        None | Some(Value::Null) => 1,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(1).max(1) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'band' must be an integer".into()))?
            .max(1),
        _ => 1,
    };
    Ok(Params {
        method,
        neighbors,
        min_count,
        iterations,
        sort,
        band: (band_1based - 1) as isize,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, Raster, RasterConfig};

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
            nodata: -1.0,
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

    fn run(args: serde_json::Value) -> Raster {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = BoundaryCleanTool.run(&args, &ctx()).unwrap();
        load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    fn count(r: &Raster, v: f64) -> usize {
        let mut n = 0;
        for row in 0..r.rows {
            for col in 0..r.cols {
                if r.get(0, row as isize, col as isize) == v {
                    n += 1;
                }
            }
        }
        n
    }

    /// A single stray cell of class 9 in a field of 0s is majority-filtered away.
    #[test]
    fn majority_removes_speckle() {
        let mut data = vec![0.0; 25]; // 5x5
        data[12] = 9.0; // isolated center speckle
        let input = raster_from(5, 5, data);
        let r = run(json!({ "input": input, "method": "majority" }));
        assert_eq!(count(&r, 9.0), 0, "isolated speckle should be removed");
        assert_eq!(
            count(&r, 0.0),
            25,
            "everything should become the background"
        );
    }

    /// A coherent 2x2 block is NOT eroded by a majority filter (its cells each
    /// still have enough same-class neighbours to hold, and the background can't
    /// reach a strict majority against a block corner under 8-connectivity...
    /// actually corners flip, so assert the block's interior-ish survival via a
    /// larger block).
    #[test]
    fn majority_preserves_solid_region() {
        let mut data = vec![0.0; 49]; // 7x7
        for r in 1..6 {
            for c in 1..6 {
                data[r * 7 + c] = 3.0; // solid 5x5 block of class 3
            }
        }
        let before = 25; // solid 5x5 block
        let input = raster_from(7, 7, data);
        let r = run(json!({ "input": input, "method": "majority", "neighbors": 8 }));
        // The solid interior must survive; only fringe corners may flip.
        assert!(
            count(&r, 3.0) as f64 >= before as f64 * 0.7,
            "a solid region must largely survive majority filtering"
        );
    }

    /// Expand/shrink boundary clean removes a one-cell protrusion sticking out of
    /// a large zone into a small one, and keeps the categorical palette.
    #[test]
    fn expand_shrink_smooths_protrusion() {
        // 7x7: left half class 1 (big), right half class 2, with a one-cell
        // spike of class 1 poking into the class-2 area.
        let mut data = vec![2.0; 49];
        for r in 0..7 {
            for c in 0..3 {
                data[r * 7 + c] = 1.0;
            }
        }
        data[3 * 7 + 4] = 1.0; // spike into class 2
        let input = raster_from(7, 7, data);
        let r = run(json!({ "input": input, "method": "expand_shrink", "sort": "descending" }));
        // Output values stay within the original class set.
        for row in 0..r.rows {
            for col in 0..r.cols {
                let v = r.get(0, row as isize, col as isize);
                assert!(v == 1.0 || v == 2.0, "no new classes introduced");
            }
        }
        // The spike should be reabsorbed by the surrounding class 2.
        assert_eq!(
            r.get(0, 3, 4),
            2.0,
            "the one-cell class-1 spike should be smoothed away"
        );
    }

    /// No-data cells are preserved and never overwritten.
    #[test]
    fn nodata_is_preserved() {
        let mut data = vec![0.0; 25];
        data[12] = -1.0; // nodata in the middle
        data[6] = 7.0;
        let input = raster_from(5, 5, data);
        let r = run(json!({ "input": input, "method": "majority" }));
        assert_eq!(r.get(0, 2, 2), -1.0, "nodata must be preserved");
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            BoundaryCleanTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.tif", "method": "clean" })).is_err());
        assert!(bad(json!({ "input": "a.tif", "neighbors": 6 })).is_err());
        assert!(bad(json!({ "input": "a.tif", "threshold": "most" })).is_err());
        assert!(bad(json!({ "input": "a.tif", "sort": "biggest" })).is_err());
        assert!(
            bad(json!({ "input": "a.tif", "method": "expand_shrink", "sort": "ascending" }))
                .is_ok()
        );
    }
}
