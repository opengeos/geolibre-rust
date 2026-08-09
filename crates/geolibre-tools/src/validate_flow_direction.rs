//! GeoLibre tool: find the cells where a flow-direction raster is broken.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Validate Flow Direction* (Spatial
//! Analyst).
//!
//! ## Why the catalog needs it
//!
//! A flow-direction raster is the input to nearly every hydrology tool the
//! catalog ships — accumulation, watersheds, stream extraction, distance to
//! stream. When one is malformed the downstream tools do not error: they return
//! a plausible-looking answer computed from a broken graph, or they hang. A
//! circular flow path is the worst case, because a routing loop is an infinite
//! loop.
//!
//! Such rasters are common in practice: they arrive from another package with a
//! different direction encoding, get resampled with bilinear interpolation
//! (which produces codes like 12 that mean nothing), get clipped so flow now
//! runs off an edge into nothing, or get hand-edited.
//!
//! The hydrology suite is large — `d8_pointer`, `d8_flow_accumulation`,
//! `fd8_flow_accumulation`, `watershed`, `basins`, `trace_downslope_flowpaths`
//! and more — and every one of them *consumes* a pointer raster. Not one checks
//! that the raster it was handed is sane.
//!
//! ## What is checked
//!
//! * **invalid_code** — a value that is not a power-of-two D8 direction.
//! * **undefined** — a zero direction away from the raster edge, i.e. a sink
//!   with no outlet, which stalls routing.
//! * **flows_off_edge** — flow leaving the raster. Normal on a true boundary,
//!   so it is reported as informational and can be switched off.
//! * **flows_to_nodata** — flow into a no-data cell, which silently truncates a
//!   flow path mid-basin.
//! * **circular_flow** — the cell lies on a cycle. Found by three-colour
//!   traversal so each cell is visited once and the check itself cannot loop.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::args_common::{band_index, bool_or, choice_or, req_str};
use crate::common::{load_input_raster, parse_optional_output};
use crate::vector_common::write_or_store_layer;

/// D8 codes in the ESRI/whitebox convention, and the (row, col) step each
/// means. Row increases southward.
const D8: [(f64, isize, isize); 8] = [
    (1.0, 0, 1),    // east
    (2.0, 1, 1),    // south-east
    (4.0, 1, 0),    // south
    (8.0, 1, -1),   // south-west
    (16.0, 0, -1),  // west
    (32.0, -1, -1), // north-west
    (64.0, -1, 0),  // north
    (128.0, -1, 1), // north-east
];

pub struct ValidateFlowDirectionTool;

impl Tool for ValidateFlowDirectionTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "validate_flow_direction",
            display_name: "Validate Flow Direction",
            summary: "Reports the cells that make a flow-direction raster unusable — codes that are not D8 directions, interior cells with no outlet, flow into no-data, flow off the edge, and cells on a circular flow path — as a point layer of problems (ArcGIS Validate Flow Direction). The bundled hydrology suite consumes pointer rasters everywhere (d8_flow_accumulation, watershed, basins, trace_downslope_flowpaths) but nothing validates one, so a raster broken by resampling, reclassing or clipping yields a plausible wrong answer rather than an error, and a routing loop can hang.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Flow-direction (pointer) raster to check.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output point layer, one point per problem cell. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "flow_direction_type",
                    description: "'d8' (default) for single power-of-two codes, or 'mfd' where a value may be the sum of several directions.",
                    required: false,
                },
                ToolParamSpec {
                    name: "report_edge_outflow",
                    description: "Report cells whose flow leaves the raster (default false). Normal on a true boundary, a real defect on a clipped interior.",
                    required: false,
                },
                ToolParamSpec {
                    name: "check_circular",
                    description: "Detect circular flow paths (default true).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band holding the directions (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input_path = req_str(args, "input")?.to_string();
        let prm = parse_params(args)?;
        let band = band_index(args, "band")?;
        let output = parse_optional_output(args, "output")?;

        let raster = load_input_raster(&input_path)?;
        let (rows, cols) = (raster.rows, raster.cols);

        // Directions, with no-data as None.
        let mut dir = vec![None::<f64>; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let v = raster.get(band, r as isize, c as isize);
                if v != raster.nodata && v.is_finite() {
                    dir[r * cols + c] = Some(v);
                }
            }
        }

        // Problem per cell; a cell reports only its most serious problem, so a
        // single defect is not counted several times.
        let mut problem: Vec<Option<Problem>> = vec![None; rows * cols];
        // Downstream neighbour, for the cycle walk. Only set for cells whose
        // direction is a single valid step landing on a valid cell.
        let mut next = vec![None::<usize>; rows * cols];

        for r in 0..rows {
            for c in 0..cols {
                let i = r * cols + c;
                let Some(code) = dir[i] else {
                    continue;
                };
                if code == 0.0 {
                    // A zero on the border is an outlet; inside it is a stall.
                    if r > 0 && c > 0 && r + 1 < rows && c + 1 < cols {
                        problem[i] = Some(Problem::Undefined);
                    }
                    continue;
                }
                let steps = decode(code, prm.mfd);
                if steps.is_empty() {
                    problem[i] = Some(Problem::InvalidCode);
                    continue;
                }

                // Under MFD a value is the sum of several directions, so every
                // component is followed; the single-successor cycle walk only
                // applies when exactly one direction is set.
                let mut off_edge = false;
                let mut to_nodata = false;
                let mut single: Option<usize> = None;
                for &(dr, dc) in &steps {
                    let nr = r as isize + dr;
                    let nc = c as isize + dc;
                    if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                        off_edge = true;
                        continue;
                    }
                    let j = nr as usize * cols + nc as usize;
                    if dir[j].is_none() {
                        to_nodata = true;
                        continue;
                    }
                    if steps.len() == 1 {
                        single = Some(j);
                    }
                }
                next[i] = single;

                // Flow into no-data truncates a path mid-basin, which is worse
                // than leaving the raster at its edge.
                if to_nodata {
                    problem[i] = Some(Problem::FlowsToNodata);
                } else if off_edge && prm.report_edge_outflow {
                    problem[i] = Some(Problem::FlowsOffEdge);
                }
            }
        }

        let mut cycle_cells = 0usize;
        if prm.check_circular {
            cycle_cells = mark_cycles(&next, &mut problem, rows * cols);
        }

        ctx.progress.info(&format!(
            "{rows}x{cols}, {} problem cell(s){}",
            problem.iter().filter(|p| p.is_some()).count(),
            if cycle_cells > 0 {
                format!(", {cycle_cells} on circular paths")
            } else {
                String::new()
            }
        ));

        // Output points at the centres of the problem cells.
        let mut layer = Layer::new("flow_direction_problems");
        layer.geom_type = Some(GeometryType::Point);
        if let Some(e) = raster.crs.epsg {
            layer = layer.with_crs_epsg(e);
        }
        layer.add_field(FieldDef::new("id", FieldType::Integer));
        layer.add_field(FieldDef::new("row", FieldType::Integer));
        layer.add_field(FieldDef::new("col", FieldType::Integer));
        layer.add_field(FieldDef::new("value", FieldType::Float));
        layer.add_field(FieldDef::new("problem", FieldType::Text));

        let (csx, csy) = (raster.cell_size_x, raster.cell_size_y);
        let y_max = raster.y_min + rows as f64 * csy;
        let mut counts: BTreeMap<&'static str, usize> = BTreeMap::new();
        let mut fid = 0u64;

        for r in 0..rows {
            for c in 0..cols {
                let i = r * cols + c;
                let Some(p) = problem[i] else {
                    continue;
                };
                *counts.entry(p.label()).or_default() += 1;
                let x = raster.x_min + (c as f64 + 0.5) * csx;
                let y = y_max - (r as f64 + 0.5) * csy;
                let mut f = Feature::with_geometry(
                    fid,
                    Geometry::Point(Coord::xy(x, y)),
                    layer.schema.len(),
                );
                f.set_by_index(0, FieldValue::Integer(fid as i64));
                f.set_by_index(1, FieldValue::Integer(r as i64));
                f.set_by_index(2, FieldValue::Integer(c as i64));
                f.set_by_index(
                    3,
                    match dir[i] {
                        Some(v) => FieldValue::Float(v),
                        None => FieldValue::Null,
                    },
                );
                f.set_by_index(4, FieldValue::Text(p.label().to_string()));
                layer.push(f);
                fid += 1;
            }
        }

        let problem_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("problem_count".to_string(), json!(problem_count));
        outputs.insert("valid".to_string(), json!(problem_count == 0));
        outputs.insert("counts_by_problem".to_string(), json!(counts));
        outputs.insert(
            "flow_direction_type".to_string(),
            json!(if prm.mfd { "mfd" } else { "d8" }),
        );
        Ok(ToolRunResult { outputs })
    }
}

/// Decodes a direction value into its (row, col) steps.
///
/// Under D8 the value must be exactly one of the eight codes. Under MFD it may
/// be their sum, so every set bit contributes a step; a value with bits outside
/// the eight directions is still invalid.
fn decode(code: f64, mfd: bool) -> Vec<(isize, isize)> {
    // Directions are integers; a fractional value cannot be one.
    if code < 0.0 || code.fract() != 0.0 || code > 255.0 {
        return Vec::new();
    }
    if !mfd {
        return D8
            .iter()
            .find(|(v, _, _)| *v == code)
            .map(|&(_, dr, dc)| vec![(dr, dc)])
            .unwrap_or_default();
    }
    let bits = code as u32;
    let mut out = Vec::new();
    for &(v, dr, dc) in &D8 {
        if bits & (v as u32) != 0 {
            out.push((dr, dc));
        }
    }
    out
}

/// Marks every cell lying on a cycle of the single-successor graph.
///
/// Three-colour traversal: white unvisited, grey on the current path, black
/// finished. Meeting a grey cell closes a cycle, and every cell from that point
/// back around the path is on it. Each cell is coloured once, so the detector
/// itself terminates even on the pathological input it exists to find.
fn mark_cycles(next: &[Option<usize>], problem: &mut [Option<Problem>], n: usize) -> usize {
    #[derive(Clone, Copy, PartialEq)]
    enum Colour {
        White,
        Grey,
        Black,
    }
    let mut colour = vec![Colour::White; n];
    let mut marked = 0usize;
    let mut path: Vec<usize> = Vec::new();

    for start in 0..n {
        if colour[start] != Colour::White {
            continue;
        }
        path.clear();
        let mut cur = start;
        loop {
            match colour[cur] {
                Colour::White => {
                    colour[cur] = Colour::Grey;
                    path.push(cur);
                    match next[cur] {
                        Some(j) => cur = j,
                        None => break,
                    }
                }
                Colour::Grey => {
                    // Closed a loop: everything from `cur` onward in the
                    // current path is on the cycle.
                    if let Some(pos) = path.iter().position(|&x| x == cur) {
                        for &cell in &path[pos..] {
                            // Do not overwrite a more specific diagnosis.
                            if problem[cell].is_none() {
                                problem[cell] = Some(Problem::CircularFlow);
                                marked += 1;
                            }
                        }
                    }
                    break;
                }
                Colour::Black => break,
            }
        }
        for &cell in &path {
            colour[cell] = Colour::Black;
        }
    }
    marked
}

/// What is wrong with a cell, most serious first.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Problem {
    InvalidCode,
    Undefined,
    FlowsToNodata,
    FlowsOffEdge,
    CircularFlow,
}

impl Problem {
    fn label(self) -> &'static str {
        match self {
            Problem::InvalidCode => "invalid_code",
            Problem::Undefined => "undefined",
            Problem::FlowsToNodata => "flows_to_nodata",
            Problem::FlowsOffEdge => "flows_off_edge",
            Problem::CircularFlow => "circular_flow",
        }
    }
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    mfd: bool,
    report_edge_outflow: bool,
    check_circular: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let mfd = choice_or(args, "flow_direction_type", &["d8", "mfd"], "d8")? == "mfd";
    Ok(Params {
        mfd,
        report_edge_outflow: bool_or(args, "report_edge_outflow", false)?,
        check_circular: bool_or(args, "check_circular", true)?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vector_common::load_input_layer;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
    use serde_json::Value;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn raster_of(cols: usize, rows: usize, vals: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 10.0,
            cell_size_y: Some(10.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(32610),
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
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, BTreeMap<String, Value>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ValidateFlowDirectionTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (layer, out.outputs)
    }

    fn problems(layer: &Layer) -> Vec<(i64, i64, String)> {
        let ri = layer.schema.field_index("row").unwrap();
        let ci = layer.schema.field_index("col").unwrap();
        let pi = layer.schema.field_index("problem").unwrap();
        let mut out: Vec<(i64, i64, String)> = layer
            .iter()
            .map(|f| {
                let (FieldValue::Integer(r), FieldValue::Integer(c), FieldValue::Text(p)) =
                    (&f.attributes[ri], &f.attributes[ci], &f.attributes[pi])
                else {
                    panic!("unexpected attribute types")
                };
                (*r, *c, p.clone())
            })
            .collect();
        out.sort();
        out
    }

    /// A well-formed pointer raster reports nothing.
    #[test]
    fn a_valid_raster_is_clean() {
        // Everything flows east; the last column leaves the raster, which is
        // normal boundary behaviour and off by default.
        let (rows, cols) = (3, 4);
        let (layer, outputs) = run(json!({
            "input": raster_of(cols, rows, &[1.0; 12])
        }));
        assert_eq!(layer.len(), 0, "a clean raster should report nothing");
        assert!(outputs["valid"].as_bool().unwrap());
    }

    /// A resampled or reclassed raster picks up codes that are not directions;
    /// these are the single most common defect and must be caught.
    #[test]
    fn invalid_codes_are_reported() {
        let (rows, cols) = (3, 3);
        let mut v = vec![1.0; 9];
        v[4] = 12.0; // not a power of two: what bilinear resampling produces
        v[7] = 3.5; // fractional
        let (layer, outputs) = run(json!({ "input": raster_of(cols, rows, &v) }));
        assert!(!outputs["valid"].as_bool().unwrap());
        let found = problems(&layer);
        assert!(found.contains(&(1, 1, "invalid_code".to_string())), "{found:?}");
        assert!(found.contains(&(2, 1, "invalid_code".to_string())), "{found:?}");
    }

    /// A circular flow path is the failure that hangs downstream routing, so it
    /// must be found — and the detector itself must terminate.
    #[test]
    fn circular_flow_is_detected() {
        // A 2x2 loop in the middle of a 4x4: (1,1)->E->(1,2)->S->(2,2)->W->
        // (2,1)->N->(1,1).
        let (rows, cols) = (4, 4);
        let mut v = vec![1.0; rows * cols];
        let at = |r: usize, c: usize| r * cols + c;
        v[at(1, 1)] = 1.0; // east
        v[at(1, 2)] = 4.0; // south
        v[at(2, 2)] = 16.0; // west
        v[at(2, 1)] = 64.0; // north
        let (layer, outputs) = run(json!({ "input": raster_of(cols, rows, &v) }));
        let found = problems(&layer);
        let cycle: Vec<&(i64, i64, String)> = found
            .iter()
            .filter(|(_, _, p)| p == "circular_flow")
            .collect();
        assert_eq!(cycle.len(), 4, "all four loop cells should be flagged: {found:?}");
        assert!(!outputs["valid"].as_bool().unwrap());
        assert_eq!(
            outputs["counts_by_problem"]["circular_flow"].as_u64().unwrap(),
            4
        );
    }

    /// A two-cell ping-pong is the shortest possible cycle.
    #[test]
    fn two_cell_loop_is_a_cycle() {
        let (rows, cols) = (3, 4);
        let mut v = vec![1.0; rows * cols];
        v[cols + 1] = 1.0; // east
        v[cols + 2] = 16.0; // west, straight back
        let (layer, _) = run(json!({ "input": raster_of(cols, rows, &v) }));
        let found = problems(&layer);
        assert_eq!(
            found
                .iter()
                .filter(|(_, _, p)| p == "circular_flow")
                .count(),
            2,
            "{found:?}"
        );
    }

    /// Flow into no-data truncates a path mid-basin.
    #[test]
    fn flow_into_nodata_is_reported() {
        let (rows, cols) = (3, 3);
        let mut v = vec![1.0; 9];
        v[5] = -9999.0; // the cell east of the centre is missing
        let (layer, _) = run(json!({ "input": raster_of(cols, rows, &v) }));
        let found = problems(&layer);
        assert!(
            found.contains(&(1, 1, "flows_to_nodata".to_string())),
            "{found:?}"
        );
    }

    /// An interior cell with no direction stalls routing; a border one is a
    /// legitimate outlet.
    #[test]
    fn interior_zero_is_undefined_but_border_zero_is_not() {
        let (rows, cols) = (3, 3);
        let mut v = vec![1.0; 9];
        v[4] = 0.0; // interior
        let (layer, _) = run(json!({ "input": raster_of(cols, rows, &v) }));
        assert!(problems(&layer).contains(&(1, 1, "undefined".to_string())));

        let mut border = vec![1.0; 9];
        border[0] = 0.0; // corner
        let (clean, _) = run(json!({ "input": raster_of(cols, rows, &border) }));
        assert!(
            !problems(&clean)
                .iter()
                .any(|(_, _, p)| p == "undefined"),
            "a border outlet is not a defect"
        );
    }

    /// Edge outflow is opt-in, because on a true boundary it is expected.
    #[test]
    fn edge_outflow_is_opt_in() {
        let (rows, cols) = (2, 2);
        let v = vec![1.0; 4]; // everything flows east, off the raster
        let (quiet, _) = run(json!({ "input": raster_of(cols, rows, &v) }));
        assert_eq!(quiet.len(), 0);

        let (loud, outputs) = run(json!({
            "input": raster_of(cols, rows, &v), "report_edge_outflow": true
        }));
        assert_eq!(loud.len(), 2, "the two right-hand cells flow off the edge");
        assert_eq!(
            outputs["counts_by_problem"]["flows_off_edge"].as_u64().unwrap(),
            2
        );
    }

    /// Under MFD a value may be the sum of several directions, so a code that
    /// is invalid for D8 can be perfectly valid here.
    #[test]
    fn mfd_accepts_summed_directions() {
        let (rows, cols) = (3, 3);
        let mut v = vec![1.0; 9];
        v[4] = 5.0; // east + south, legal under MFD
        let (d8, _) = run(json!({ "input": raster_of(cols, rows, &v) }));
        assert!(
            problems(&d8).contains(&(1, 1, "invalid_code".to_string())),
            "5 is not a D8 code"
        );

        let (mfd, outputs) = run(json!({
            "input": raster_of(cols, rows, &v), "flow_direction_type": "mfd"
        }));
        assert!(
            !problems(&mfd)
                .iter()
                .any(|(r, c, p)| *r == 1 && *c == 1 && p == "invalid_code"),
            "5 = east|south is valid under MFD"
        );
        assert_eq!(outputs["flow_direction_type"].as_str().unwrap(), "mfd");
    }

    /// The cycle check can be switched off.
    #[test]
    fn circular_check_can_be_disabled() {
        let (rows, cols) = (3, 4);
        let mut v = vec![1.0; rows * cols];
        v[cols + 1] = 1.0;
        v[cols + 2] = 16.0;
        let (off, _) = run(json!({
            "input": raster_of(cols, rows, &v), "check_circular": false
        }));
        assert!(!problems(&off).iter().any(|(_, _, p)| p == "circular_flow"));
    }

    /// A long chain that terminates is not a cycle — the traversal must not
    /// mistake a shared downstream path for a loop.
    #[test]
    fn converging_paths_are_not_cycles() {
        // Two rows both flowing east into a common column, then off the edge.
        let (rows, cols) = (2, 6);
        let mut v = vec![1.0; rows * cols];
        v[cols] = 128.0; // (1,0) flows north-east into row 0
        let (layer, outputs) = run(json!({ "input": raster_of(cols, rows, &v) }));
        assert!(
            !problems(&layer).iter().any(|(_, _, p)| p == "circular_flow"),
            "converging flow is not circular: {:?}",
            problems(&layer)
        );
        assert!(outputs["valid"].as_bool().unwrap());
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ValidateFlowDirectionTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({"input": "a.tif", "flow_direction_type": "dinf"})).is_err());
        assert!(bad(json!({"input": "a.tif", "flow_direction_type": "mfd"})).is_ok());
    }
}
