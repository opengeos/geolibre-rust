//! GeoLibre tool: trace a cost backlink raster into least-cost polylines.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Optimal Path As Line* (Spatial
//! Analyst Distance toolset), which also covers *Cost Path As Polyline*.
//!
//! ## Why this is the missing half
//!
//! The bundled `cost_pathway` marks the least-cost route as a **raster** of
//! path cells. That is fine for display and useless for everything else: you
//! cannot measure it, join to it, snap to it, or serve it as vector tiles.
//!
//! `cost_back_link` (GeoLibre) already publishes the reusable half of the solve
//! — the 0-8 direction raster telling every cell which neighbour steps toward
//! the source. Nothing consumed it into geometry. `cost_connectivity` does emit
//! polylines, but only for the N-site least-cost *network* (MST / all-neighbour
//! pairs); the ordinary source-to-destination path had no vector form at all.
//!
//! Turning raster results into clean vector output is this crate's core
//! business (`polygonize`, `raster_to_vector_*`, `contour_list`), so the
//! cost-distance family being stuck as pixels was a real hole.
//!
//! ## Backlink encoding
//!
//! The ArcGIS convention that `cost_back_link` writes: `0` marks a source cell
//! and `1`-`8` give the direction of the **next cell on the way back to the
//! source**, `1` = east, proceeding counter-clockwise. Unreachable cells are
//! no-data. This module decodes exactly that; it does not accept the D8 flow
//! encoding (1,2,4,8,...), which is a different, non-overlapping convention.
//!
//! ## Path types
//!
//! * `each_cell` — one path per destination cell (the ArcGIS "each cell" rule).
//! * `each_zone` — one path per destination zone, starting from that zone's
//!   cheapest cell. Requires `accumulation` to pick the cheapest.
//! * `best_single` — the single cheapest path over all destinations.
//!
//! ## The cycle guard
//!
//! A backlink raster produced by `cost_back_link` is acyclic by construction —
//! it is a Dijkstra predecessor tree. A hand-edited or externally-produced one
//! need not be, and a naive walk would spin forever. Every trace therefore
//! carries a visited set and fails loudly rather than hanging.

use std::collections::{BTreeMap, HashSet};

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::args_common::{choice_or, req_str};
use crate::common::load_input_raster;
use crate::raster_stack::check_alignment_refs;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Neighbour offsets in ArcGIS backlink order: code `k + 1`, 1 = east,
/// counter-clockwise. Must stay identical to `cost_back_link::NEIGHBOURS`.
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

const PATH_TYPES: [&str; 3] = ["each_cell", "each_zone", "best_single"];

pub struct OptimalPathAsLineTool;

impl Tool for OptimalPathAsLineTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "optimal_path_as_line",
            display_name: "Optimal Path As Line",
            summary: "Traces a cost backlink raster from each destination back to its source and emits the least-cost routes as polylines (ArcGIS Optimal Path As Line / Cost Path As Polyline). cost_pathway only marks the route as raster cells and cost_connectivity only covers the N-site network case, so the ordinary source-to-destination path had no vector form.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "destination",
                    description: "Destinations to trace from: a raster (any non-zero, non-no-data cell) or a point/multipoint vector layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "backlink",
                    description: "Backlink direction raster from cost_back_link (0 = source, 1-8 = direction toward the source, 1 = east counter-clockwise).",
                    required: true,
                },
                ToolParamSpec {
                    name: "accumulation",
                    description: "Optional accumulated-cost raster (cost_back_link's out_distance). Supplies each path's total cost and is required by path_type 'each_zone' and 'best_single'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "path_type",
                    description: "One of 'each_cell' (default, one path per destination cell), 'each_zone' (one path per destination zone, from its cheapest cell), 'best_single' (the single cheapest path).",
                    required: false,
                },
                ToolParamSpec {
                    name: "zone_field",
                    description: "Attribute holding the zone id for a vector destination under path_type 'each_zone'. Raster destinations use the cell value as the zone.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polyline layer. If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "destination")?;
        req_str(args, "backlink")?;
        let path_type = choice_or(args, "path_type", &PATH_TYPES, "each_cell")?;
        if matches!(path_type, "each_zone" | "best_single") && args.get("accumulation").is_none() {
            return Err(ToolError::Validation(format!(
                "path_type '{path_type}' ranks destinations by cost, so 'accumulation' is required"
            )));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let dest_spec = req_str(args, "destination")?.to_string();
        let backlink_path = req_str(args, "backlink")?.to_string();
        let accum_path = parse_optional_str(args, "accumulation")?.map(str::to_string);
        let path_type = choice_or(args, "path_type", &PATH_TYPES, "each_cell")?;
        let zone_field = parse_optional_str(args, "zone_field")?.map(str::to_string);
        let output = parse_optional_str(args, "output")?;

        let backlink = load_input_raster(&backlink_path)?;
        let accum = match &accum_path {
            Some(p) => {
                let a = load_input_raster(p)?;
                check_alignment_refs(&[&backlink, &a])?;
                Some(a)
            }
            None => None,
        };

        let rows = backlink.rows;
        let cols = backlink.cols;
        let starts = load_destinations(&dest_spec, &backlink, zone_field.as_deref())?;
        if starts.is_empty() {
            return Err(ToolError::Execution(
                "'destination' selected no cells on the backlink grid".to_string(),
            ));
        }
        ctx.progress.info(&format!(
            "{rows}x{cols}, {} destination cell(s), path_type={path_type}",
            starts.len()
        ));

        let cost_at = |idx: usize| -> Option<f64> {
            let a = accum.as_ref()?;
            let v = a.get(0, (idx / cols) as isize, (idx % cols) as isize);
            (v != a.nodata && v.is_finite()).then_some(v)
        };

        // Reduce the destination set according to path_type before tracing:
        // 'each_zone' keeps the cheapest cell per zone, 'best_single' the
        // cheapest overall. Doing this first means we never trace a path we
        // are about to throw away.
        let selected: Vec<Dest> = match path_type {
            "each_cell" => starts,
            "each_zone" => {
                let mut best: BTreeMap<i64, (f64, Dest)> = BTreeMap::new();
                for d in starts {
                    let Some(c) = cost_at(d.idx) else { continue };
                    match best.get(&d.zone) {
                        Some((bc, _)) if *bc <= c => {}
                        _ => {
                            best.insert(d.zone, (c, d));
                        }
                    }
                }
                best.into_values().map(|(_, d)| d).collect()
            }
            _ => {
                let mut best: Option<(f64, Dest)> = None;
                for d in starts {
                    let Some(c) = cost_at(d.idx) else { continue };
                    if best.as_ref().is_none_or(|(bc, _)| c < *bc) {
                        best = Some((c, d));
                    }
                }
                best.into_iter().map(|(_, d)| d).collect()
            }
        };
        if selected.is_empty() {
            return Err(ToolError::Execution(
                "no destination cell had a valid accumulated cost; check that 'accumulation' \
                 matches the backlink raster"
                    .to_string(),
            ));
        }

        let mut out = Layer::new("optimal_paths").with_geom_type(GeometryType::LineString);
        if let Some(e) = backlink.crs.epsg {
            out = out.with_crs_epsg(e);
        }
        out.add_field(FieldDef::new("DestID", FieldType::Integer));
        out.add_field(FieldDef::new("Zone", FieldType::Integer));
        out.add_field(FieldDef::new("SrcRow", FieldType::Integer));
        out.add_field(FieldDef::new("SrcCol", FieldType::Integer));
        out.add_field(FieldDef::new("PathCost", FieldType::Float));
        out.add_field(FieldDef::new("PathCells", FieldType::Integer));

        // Cell centres. Rows count down from the top, so y is measured from the
        // raster's top edge rather than from y_min.
        let y_max = backlink.y_min + rows as f64 * backlink.cell_size_y;
        let cell_xy = |idx: usize| -> Coord {
            let r = idx / cols;
            let c = idx % cols;
            Coord::xy(
                backlink.x_min + (c as f64 + 0.5) * backlink.cell_size_x,
                y_max - (r as f64 + 0.5) * backlink.cell_size_y,
            )
        };

        let mut unreachable = 0_u64;
        let mut traced = 0_u64;
        for (dest_id, d) in selected.iter().enumerate() {
            let cells = trace(d.idx, &backlink, rows, cols)?;
            let Some(cells) = cells else {
                unreachable += 1;
                continue;
            };
            // A destination sitting on a source cell is a zero-length path, not
            // an error, but it has no line geometry to emit.
            if cells.len() < 2 {
                continue;
            }
            let source = *cells.last().expect("non-empty trace");
            let coords: Vec<Coord> = cells.iter().map(|&c| cell_xy(c)).collect();
            traced += 1;
            out.push(Feature {
                fid: 0,
                geometry: Some(Geometry::line_string(coords)),
                attributes: vec![
                    FieldValue::Integer(dest_id as i64),
                    FieldValue::Integer(d.zone),
                    FieldValue::Integer((source / cols) as i64),
                    FieldValue::Integer((source % cols) as i64),
                    match cost_at(d.idx) {
                        Some(c) => FieldValue::Float(c),
                        None => FieldValue::Null,
                    },
                    FieldValue::Integer(cells.len() as i64),
                ],
            });
        }

        let path_count = out.features.len();
        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("path_count".to_string(), json!(path_count));
        outputs.insert("destination_cells".to_string(), json!(selected.len()));
        outputs.insert("traced".to_string(), json!(traced));
        outputs.insert("unreachable".to_string(), json!(unreachable));
        Ok(ToolRunResult { outputs })
    }
}

/// A destination cell plus the zone it belongs to.
struct Dest {
    idx: usize,
    zone: i64,
}

/// Walks the backlink raster from `start` to a source cell (code 0).
///
/// Returns `Ok(None)` when the destination is unreachable (no-data backlink),
/// and an error when the raster contains a cycle or steps off the grid — both
/// signal a malformed backlink, and hanging or silently truncating would be
/// worse than saying so.
fn trace(
    start: usize,
    backlink: &Raster,
    rows: usize,
    cols: usize,
) -> Result<Option<Vec<usize>>, ToolError> {
    let code_at =
        |idx: usize| -> f64 { backlink.get(0, (idx / cols) as isize, (idx % cols) as isize) };
    let first = code_at(start);
    if first == backlink.nodata || !first.is_finite() {
        return Ok(None);
    }

    let mut path = vec![start];
    let mut seen: HashSet<usize> = HashSet::new();
    seen.insert(start);
    let mut cur = start;
    loop {
        let code = code_at(cur);
        if code == backlink.nodata || !code.is_finite() {
            return Err(ToolError::Execution(format!(
                "backlink raster has no-data at cell ({}, {}) mid-path; the accumulation and \
                 backlink rasters do not agree",
                cur / cols,
                cur % cols
            )));
        }
        let code = code.round() as i64;
        if code == 0 {
            return Ok(Some(path)); // reached a source
        }
        if !(1..=8).contains(&code) {
            return Err(ToolError::Execution(format!(
                "backlink value {code} at cell ({}, {}) is outside the 0-8 ArcGIS encoding; the \
                 D8 flow-direction encoding (1,2,4,8,...) is not accepted here",
                cur / cols,
                cur % cols
            )));
        }
        let (dr, dc) = NEIGHBOURS[(code - 1) as usize];
        let nr = (cur / cols) as isize + dr;
        let nc = (cur % cols) as isize + dc;
        if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
            return Err(ToolError::Execution(format!(
                "backlink at cell ({}, {}) points off the grid edge",
                cur / cols,
                cur % cols
            )));
        }
        let next = nr as usize * cols + nc as usize;
        if !seen.insert(next) {
            return Err(ToolError::Execution(format!(
                "backlink raster contains a cycle at cell ({nr}, {nc}); it is not a valid \
                 predecessor tree"
            )));
        }
        path.push(next);
        cur = next;
    }
}

/// Resolves the destination parameter to cell indices, each tagged with a zone.
///
/// Mirrors `cost_back_link::load_sources`: raster destinations use any non-zero,
/// non-no-data cell (with the cell value as the zone), vector destinations snap
/// each point to its containing cell.
fn load_destinations(
    spec: &str,
    template: &Raster,
    zone_field: Option<&str>,
) -> Result<Vec<Dest>, ToolError> {
    let rows = template.rows;
    let cols = template.cols;

    // Probe as a raster first, exactly as cost_back_link does: a raster that
    // loads but fails alignment is a real error and must surface as one rather
    // than degrading into a confusing "not a vector layer" message.
    if let Ok(raster) = load_input_raster(spec) {
        check_alignment_refs(&[template, &raster])?;
        let mut out = Vec::new();
        for r in 0..rows {
            for c in 0..cols {
                let v = raster.get(0, r as isize, c as isize);
                if v != raster.nodata && v.is_finite() && v != 0.0 {
                    out.push(Dest {
                        idx: r * cols + c,
                        zone: v.round() as i64,
                    });
                }
            }
        }
        return Ok(out);
    }

    let layer = load_input_layer(spec)?;
    let zone_idx = match zone_field {
        Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
            ToolError::Validation(format!(
                "zone_field '{f}' not found in the destination layer"
            ))
        })?),
        None => None,
    };
    let mut out = Vec::new();
    for (fid, feature) in layer.iter().enumerate() {
        let zone = match zone_idx {
            Some(i) => feature
                .attributes
                .get(i)
                .and_then(field_as_i64)
                .unwrap_or(fid as i64),
            None => fid as i64,
        };
        let push = |x: f64, y: f64, out: &mut Vec<Dest>| {
            let c = ((x - template.x_min) / template.cell_size_x).floor();
            let r = (rows as f64 - 1.0) - ((y - template.y_min) / template.cell_size_y).floor();
            if c >= 0.0 && r >= 0.0 && (c as usize) < cols && (r as usize) < rows {
                out.push(Dest {
                    idx: r as usize * cols + c as usize,
                    zone,
                });
            }
        };
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
    Ok(out)
}

fn field_as_i64(v: &FieldValue) -> Option<i64> {
    match v {
        FieldValue::Integer(i) => Some(*i),
        FieldValue::Float(f) if f.is_finite() => Some(f.round() as i64),
        FieldValue::Text(s) => s.trim().parse::<i64>().ok(),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, RasterConfig};

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
            data_type: DataType::F64,
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

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = OptimalPathAsLineTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    /// Source at the west end of a 1x4 row; every other cell points west (5).
    fn west_backlink() -> String {
        raster(1, 4, &[0.0, 5.0, 5.0, 5.0])
    }

    #[test]
    fn traces_a_straight_run_back_to_the_source() {
        let dest = raster(1, 4, &[0.0, 0.0, 0.0, 1.0]);
        let (layer, res) = run(json!({"destination": dest, "backlink": west_backlink()}));
        assert_eq!(res.outputs["path_count"], json!(1));
        let f = &layer.features[0];
        let Some(Geometry::LineString(coords)) = f.geometry.as_ref() else {
            panic!("expected a LineString");
        };
        // Four cell centres, destination first, source last.
        assert_eq!(coords.len(), 4);
        assert!(
            (coords[0].x - 3.5).abs() < 1e-9,
            "starts at the destination"
        );
        assert!((coords[3].x - 0.5).abs() < 1e-9, "ends at the source");
    }

    #[test]
    fn the_mirror_case_traces_east() {
        // If the direction decoding were off by the 4-apart opposite rule, one
        // of this pair would pass and the other would not.
        let backlink = raster(1, 4, &[1.0, 1.0, 1.0, 0.0]);
        let dest = raster(1, 4, &[1.0, 0.0, 0.0, 0.0]);
        let (layer, _) = run(json!({"destination": dest, "backlink": backlink}));
        let Some(Geometry::LineString(coords)) = layer.features[0].geometry.as_ref() else {
            panic!("expected a LineString");
        };
        assert!((coords[0].x - 0.5).abs() < 1e-9);
        assert!((coords[coords.len() - 1].x - 3.5).abs() < 1e-9);
    }

    #[test]
    fn a_diagonal_backlink_walks_diagonally() {
        // 2x2, source at the north-west corner, opposite corner steps NW (4).
        let backlink = raster(2, 2, &[0.0, 5.0, 3.0, 4.0]);
        let dest = raster(2, 2, &[0.0, 0.0, 0.0, 1.0]);
        let (layer, _) = run(json!({"destination": dest, "backlink": backlink}));
        let Some(Geometry::LineString(coords)) = layer.features[0].geometry.as_ref() else {
            panic!("expected a LineString");
        };
        assert_eq!(coords.len(), 2, "one diagonal hop to the source");
        assert!((coords[1].x - 0.5).abs() < 1e-9 && (coords[1].y - 1.5).abs() < 1e-9);
    }

    #[test]
    fn path_cost_comes_from_the_accumulation_raster() {
        let dest = raster(1, 4, &[0.0, 0.0, 0.0, 1.0]);
        let accum = raster(1, 4, &[0.0, 1.0, 2.0, 3.0]);
        let (layer, _) = run(json!({
            "destination": dest, "backlink": west_backlink(), "accumulation": accum,
        }));
        let idx = layer.schema.field_index("PathCost").unwrap();
        assert_eq!(layer.features[0].attributes[idx], FieldValue::Float(3.0));
    }

    #[test]
    fn each_cell_emits_one_path_per_destination() {
        let dest = raster(1, 4, &[0.0, 1.0, 1.0, 1.0]);
        let (layer, res) = run(json!({"destination": dest, "backlink": west_backlink()}));
        assert_eq!(res.outputs["path_count"], json!(3));
        assert_eq!(layer.features.len(), 3);
    }

    #[test]
    fn each_zone_keeps_only_the_cheapest_cell_per_zone() {
        // Two destination cells both in zone 7; only the cheaper one survives.
        let dest = raster(1, 4, &[0.0, 7.0, 0.0, 7.0]);
        let accum = raster(1, 4, &[0.0, 1.0, 2.0, 3.0]);
        let (layer, res) = run(json!({
            "destination": dest, "backlink": west_backlink(),
            "accumulation": accum, "path_type": "each_zone",
        }));
        assert_eq!(res.outputs["path_count"], json!(1));
        let idx = layer.schema.field_index("PathCost").unwrap();
        assert_eq!(
            layer.features[0].attributes[idx],
            FieldValue::Float(1.0),
            "must pick the cheaper of the two zone-7 cells"
        );
    }

    #[test]
    fn best_single_keeps_exactly_one_path_across_all_zones() {
        let dest = raster(1, 4, &[0.0, 1.0, 2.0, 3.0]);
        let accum = raster(1, 4, &[0.0, 5.0, 2.0, 9.0]);
        let (layer, res) = run(json!({
            "destination": dest, "backlink": west_backlink(),
            "accumulation": accum, "path_type": "best_single",
        }));
        assert_eq!(res.outputs["path_count"], json!(1));
        let idx = layer.schema.field_index("PathCost").unwrap();
        assert_eq!(layer.features[0].attributes[idx], FieldValue::Float(2.0));
    }

    #[test]
    fn a_destination_on_the_source_yields_no_line() {
        // Zero-length path: counted as a destination, but there is no geometry
        // to emit and that is not an error.
        let dest = raster(1, 4, &[1.0, 0.0, 0.0, 0.0]);
        let (layer, res) = run(json!({"destination": dest, "backlink": west_backlink()}));
        assert_eq!(res.outputs["path_count"], json!(0));
        assert!(layer.features.is_empty());
    }

    #[test]
    fn unreachable_destinations_are_counted_not_fatal() {
        // Cell 3 is no-data in the backlink: never reached by the cost solve.
        let backlink = raster(1, 4, &[0.0, 5.0, 5.0, -9999.0]);
        let dest = raster(1, 4, &[0.0, 0.0, 1.0, 1.0]);
        let (_, res) = run(json!({"destination": dest, "backlink": backlink}));
        assert_eq!(res.outputs["unreachable"], json!(1));
        assert_eq!(res.outputs["path_count"], json!(1));
    }

    #[test]
    fn a_cyclic_backlink_errors_instead_of_hanging() {
        // Two cells pointing at each other and no source anywhere. A naive walk
        // would loop forever; this must fail fast.
        let backlink = raster(1, 2, &[1.0, 5.0]);
        let dest = raster(1, 2, &[0.0, 1.0]);
        let args: ToolArgs =
            serde_json::from_value(json!({"destination": dest, "backlink": backlink})).unwrap();
        let err = OptimalPathAsLineTool.run(&args, &ctx()).unwrap_err();
        assert!(
            format!("{err}").contains("cycle"),
            "expected a cycle error, got: {err}"
        );
    }

    #[test]
    fn a_d8_flow_encoding_is_rejected_rather_than_silently_misread() {
        // 16 is a valid D8 flow code but outside the 0-8 backlink encoding.
        // Reading it as a backlink would produce confident nonsense.
        let backlink = raster(1, 2, &[0.0, 16.0]);
        let dest = raster(1, 2, &[0.0, 1.0]);
        let args: ToolArgs =
            serde_json::from_value(json!({"destination": dest, "backlink": backlink})).unwrap();
        let err = OptimalPathAsLineTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("0-8"), "got: {err}");
    }

    #[test]
    fn a_backlink_pointing_off_the_grid_errors() {
        // Cell 0 points west out of the raster.
        let backlink = raster(1, 2, &[5.0, 5.0]);
        let dest = raster(1, 2, &[0.0, 1.0]);
        let args: ToolArgs =
            serde_json::from_value(json!({"destination": dest, "backlink": backlink})).unwrap();
        let err = OptimalPathAsLineTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("off the grid"), "got: {err}");
    }

    #[test]
    fn rejects_bad_parameters() {
        let r = raster(1, 2, &[0.0, 5.0]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            OptimalPathAsLineTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"destination": r.clone()})));
        assert!(bad(
            json!({"destination": r.clone(), "backlink": r.clone(), "path_type": "nope"})
        ));
        // each_zone and best_single rank by cost, so accumulation is required.
        assert!(bad(
            json!({"destination": r.clone(), "backlink": r.clone(), "path_type": "each_zone"})
        ));
        assert!(bad(
            json!({"destination": r.clone(), "backlink": r, "path_type": "best_single"})
        ));
    }
}
