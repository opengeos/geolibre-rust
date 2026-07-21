//! GeoLibre tool: find the best contiguous region(s) of a target area from a
//! suitability surface.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Locate Regions* (Spatial Analyst). The
//! stack can build suitability surfaces (`fuzzy_overlay`, the bundled
//! `weighted_overlay`/`weighted_sum`) but has nothing to answer "give me the
//! best contiguous N hectares" — thresholding a suitability raster yields
//! fragmented, arbitrarily shaped blobs. `clump` labels *existing* regions; it
//! doesn't grow optimal ones.
//!
//! Regions are grown by best-first accretion from the highest-suitability seeds.
//! Each candidate cell's score blends its suitability with a compactness term
//! (`shape` weights a penalty on distance from the growing region's seed), so
//! higher `shape` yields rounder regions. A region grows until it reaches its
//! target cell count; the next region seeds from the best remaining cell outside
//! a `min_distance` buffer of the ones already chosen. Output is a raster of
//! 1-based region ids (no-data elsewhere), with per-region area and mean
//! suitability reported.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{load_input_raster, parse_optional_output, raster_like_with_data};

pub struct LocateRegionsTool;

impl Tool for LocateRegionsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "locate_regions",
            display_name: "Locate Regions",
            summary: "Find the best contiguous region(s) of a target area from a suitability raster (like ArcGIS Locate Regions): best-first region growing from suitability peaks with a shape/compactness control and inter-region spacing — the siting step that turns a fuzzy_overlay/weighted_overlay surface into actual regions, which clump and thresholding can't.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Suitability raster (higher = more suitable).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output raster of 1-based region ids (no-data elsewhere). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "total_area",
                    description: "Total target area across all regions, in CRS area units (default: 5% of the valid area).",
                    required: false,
                },
                ToolParamSpec {
                    name: "num_regions",
                    description: "Number of regions to locate (default 1). Each targets total_area / num_regions.",
                    required: false,
                },
                ToolParamSpec {
                    name: "shape",
                    description: "Compactness weight 0..1: 0 = pure suitability, 1 = strongly favour round regions (default 0.3).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_distance",
                    description: "Minimum distance between regions, in CRS units (default 0).",
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
                "band {} out of range",
                prm.band + 1
            )));
        }
        let rows = raster.rows;
        let cols = raster.cols;
        let nodata = raster.nodata;
        let cell_area = raster.cell_size_x * raster.cell_size_y;
        let cell_len = (raster.cell_size_x + raster.cell_size_y) / 2.0;

        let n = rows * cols;
        let mut suit = vec![f64::NAN; n];
        let mut valid_cells = 0usize;
        for r in 0..rows {
            for c in 0..cols {
                let v = raster.get(prm.band, r as isize, c as isize);
                if v != nodata && v.is_finite() {
                    suit[r * cols + c] = v;
                    valid_cells += 1;
                }
            }
        }
        if valid_cells == 0 {
            return Err(ToolError::Execution(
                "raster has no valid cells".to_string(),
            ));
        }

        // Target cells per region.
        let total_area = prm
            .total_area
            .unwrap_or(0.05 * valid_cells as f64 * cell_area);
        let per_region_cells = ((total_area / prm.num_regions as f64) / cell_area)
            .round()
            .max(1.0) as usize;
        let expected_radius =
            ((per_region_cells as f64 / std::f64::consts::PI).sqrt() * cell_len).max(cell_len);

        ctx.progress.info(&format!(
            "locating {} region(s) of ~{} cell(s) each",
            prm.num_regions, per_region_cells
        ));

        let mut labels = vec![0u32; n]; // 0 = unassigned
        let mut available = suit.iter().map(|s| !s.is_nan()).collect::<Vec<bool>>();
        let range = suit_range(&suit); // hoisted out of the hot loop

        let mut regions_info: Vec<(usize, f64)> = Vec::new(); // (cells, mean suitability)
        for region_id in 1..=prm.num_regions as u32 {
            // Seed: highest-suitability available cell.
            let Some(seed) = (0..n)
                .filter(|&i| available[i])
                .max_by(|&a, &b| suit[a].total_cmp(&suit[b]))
            else {
                break;
            };
            let sr = (seed / cols) as f64;
            let sc = (seed % cols) as f64;

            // Best-first growth.
            let mut heap: BinaryHeap<Cand> = BinaryHeap::new();
            heap.push(Cand {
                score: suit[seed],
                cell: seed,
            });
            let mut in_region = vec![false; n];
            let mut count = 0usize;
            let mut suit_sum = 0.0;
            while let Some(Cand { cell, .. }) = heap.pop() {
                if in_region[cell] || !available[cell] {
                    continue;
                }
                in_region[cell] = true;
                labels[cell] = region_id;
                count += 1;
                suit_sum += suit[cell];
                if count >= per_region_cells {
                    break;
                }
                let r = (cell / cols) as isize;
                let c = (cell % cols) as isize;
                for (dr, dc) in NEIGH4 {
                    let nr = r + dr;
                    let nc = c + dc;
                    if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                        continue;
                    }
                    let nidx = nr as usize * cols + nc as usize;
                    if !available[nidx] || in_region[nidx] || suit[nidx].is_nan() {
                        continue;
                    }
                    // Score = suitability minus a shape penalty on distance from the seed.
                    let dist = (((nr as f64 - sr) * raster.cell_size_y)
                        .hypot((nc as f64 - sc) * raster.cell_size_x))
                        / expected_radius;
                    let score = suit[nidx] - prm.shape * dist * range;
                    heap.push(Cand { score, cell: nidx });
                }
            }
            if count == 0 {
                break;
            }
            regions_info.push((count, suit_sum / count as f64));

            // Remove the region and its min_distance buffer from availability.
            let buffer_cells = (prm.min_distance / cell_len).ceil() as isize;
            for i in 0..n {
                if in_region[i] {
                    available[i] = false;
                }
            }
            if buffer_cells > 0 {
                let region_cells: Vec<usize> = (0..n).filter(|&i| in_region[i]).collect();
                for &rc in &region_cells {
                    let r = (rc / cols) as isize;
                    let c = (rc % cols) as isize;
                    for dr in -buffer_cells..=buffer_cells {
                        for dc in -buffer_cells..=buffer_cells {
                            let nr = r + dr;
                            let nc = c + dc;
                            if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                                continue;
                            }
                            if (dr as f64).hypot(dc as f64) <= buffer_cells as f64 {
                                available[nr as usize * cols + nc as usize] = false;
                            }
                        }
                    }
                }
            }
        }

        // Output raster: region id or nodata.
        let out_nodata = 0.0;
        let data: Vec<f64> = labels.iter().map(|&l| l as f64).collect();
        let out = raster_like_with_data(&raster, data, out_nodata, DataType::F32)?;
        let out_path = crate::common::write_or_store_output(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("regions_located".to_string(), json!(regions_info.len()));
        outputs.insert(
            "region_areas".to_string(),
            json!(regions_info
                .iter()
                .map(|(c, _)| *c as f64 * cell_area)
                .collect::<Vec<_>>()),
        );
        outputs.insert(
            "region_mean_suitability".to_string(),
            json!(regions_info.iter().map(|(_, m)| *m).collect::<Vec<_>>()),
        );
        Ok(ToolRunResult { outputs })
    }
}

const NEIGH4: [(isize, isize); 4] = [(-1, 0), (0, -1), (0, 1), (1, 0)];

/// Range (max-min) of valid suitability, used to scale the shape penalty.
fn suit_range(suit: &[f64]) -> f64 {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &s in suit {
        if s.is_finite() {
            lo = lo.min(s);
            hi = hi.max(s);
        }
    }
    if hi > lo {
        hi - lo
    } else {
        1.0
    }
}

/// Best-first candidate (max-heap by score).
struct Cand {
    score: f64,
    cell: usize,
}
impl PartialEq for Cand {
    fn eq(&self, o: &Self) -> bool {
        self.score == o.score
    }
}
impl Eq for Cand {}
impl PartialOrd for Cand {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Cand {
    fn cmp(&self, o: &Self) -> Ordering {
        self.score.total_cmp(&o.score).then(o.cell.cmp(&self.cell))
    }
}

struct Params {
    total_area: Option<f64>,
    num_regions: usize,
    shape: f64,
    min_distance: f64,
    band: isize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let total_area = opt_pos_f64(args, "total_area")?;
    let num_regions = match args.get("num_regions") {
        None | Some(Value::Null) => 1,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(1).max(1) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'num_regions' must be an integer".into()))?
            .max(1),
        _ => 1,
    };
    let shape = match args.get("shape") {
        None | Some(Value::Null) => 0.3,
        Some(Value::Number(n)) => n.as_f64().unwrap_or(0.3).clamp(0.0, 1.0),
        Some(Value::String(s)) if s.trim().is_empty() => 0.3,
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'shape' must be a number".into()))?
            .clamp(0.0, 1.0),
        _ => 0.3,
    };
    let min_distance = opt_pos_f64(args, "min_distance")?.unwrap_or(0.0);
    let band_1based = args.get("band").and_then(Value::as_u64).unwrap_or(1).max(1);
    Ok(Params {
        total_area,
        num_regions,
        shape,
        min_distance,
        band: (band_1based - 1) as isize,
    })
}

fn opt_pos_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64().filter(|v| *v >= 0.0)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be a number"))),
        _ => Ok(None),
    }
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

    fn run(args: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = LocateRegionsTool.run(&args, &ctx()).unwrap();
        let r = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, r)
    }

    /// The located region sits on the high-suitability bump and hits the target
    /// area.
    #[test]
    fn locates_on_the_suitability_peak() {
        // 20x20 low suitability with a high plateau in the top-left.
        let mut data = vec![1.0; 400];
        for r in 0..8 {
            for c in 0..8 {
                data[r * 20 + c] = 10.0;
            }
        }
        let input = raster_from(20, 20, data);
        let (out, r) = run(json!({
            "input": input, "total_area": 25.0, "num_regions": 1, "shape": 0.2,
        }));
        assert_eq!(out.outputs["regions_located"], json!(1));
        // The region should land in the high-suitability quadrant.
        let mut in_plateau = 0;
        let mut total = 0;
        for row in 0..r.rows {
            for col in 0..r.cols {
                if r.get(0, row as isize, col as isize) == 1.0 {
                    total += 1;
                    if row < 8 && col < 8 {
                        in_plateau += 1;
                    }
                }
            }
        }
        assert!(
            (20..=30).contains(&total),
            "area near target 25, got {total}"
        );
        assert!(
            in_plateau as f64 / total as f64 > 0.8,
            "region should sit on the high-suitability plateau"
        );
    }

    /// Two regions are disjoint and both are located.
    #[test]
    fn locates_multiple_disjoint_regions() {
        // Two separated high bumps.
        let mut data = vec![1.0; 400];
        for r in 2..6 {
            for c in 2..6 {
                data[r * 20 + c] = 9.0;
            }
        }
        for r in 14..18 {
            for c in 14..18 {
                data[r * 20 + c] = 9.0;
            }
        }
        let input = raster_from(20, 20, data);
        let (out, r) = run(json!({
            "input": input, "total_area": 24.0, "num_regions": 2, "min_distance": 3.0,
        }));
        assert_eq!(out.outputs["regions_located"], json!(2));
        // Labels 1 and 2 both present, no overlap.
        let mut seen = [0, 0];
        for row in 0..r.rows {
            for col in 0..r.cols {
                match r.get(0, row as isize, col as isize) as i32 {
                    1 => seen[0] += 1,
                    2 => seen[1] += 1,
                    _ => {}
                }
            }
        }
        assert!(seen[0] > 0 && seen[1] > 0, "both regions must be located");
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            LocateRegionsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.tif" })).is_ok());
    }
}
