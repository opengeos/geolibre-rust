//! GeoLibre tool: vertical cross-section panels along a line through a stack
//! of surfaces.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Fence Diagram* (3D Analyst).
//!
//! ## The gap
//!
//! `profile`, `long_profile` and `stack_profile` sample surfaces along a line,
//! but they return 2D profile **lines** or tables. None of them produces the 3D
//! panel geometry between consecutive surfaces, which is the entire point: a
//! fence diagram shows the *units* — the stratigraphic layers, the atmospheric
//! slabs — not just their boundaries.
//!
//! `voxel_isosurface` (round 16) extracts a level surface from a volume, which
//! answers a different question.
//!
//! ## Panels, not lines
//!
//! Each panel spans one consecutive surface pair along the trace, emitted as a
//! triangle strip in the `buffer_3d` convention so it renders and composes like
//! any other 3D geometry in the catalog. Panels are open ribbons rather than
//! closed solids — they are section faces, and claiming a volume for them would
//! be wrong — so `is_closed_3d` will correctly report them as open.
//!
//! Surfaces are ordered **top to bottom** as given. Where a lower surface rises
//! above the one above it (a crossing, common in interpolated geology) the
//! panel is emitted anyway with `INVERTED` set, rather than being silently
//! dropped or flipped: the crossing is usually the interesting part.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::args_common::{band_index, opt_f64, opt_positive_f64, req_str};
use crate::common::load_input_raster;
use crate::mesh3d::{tri_area, triangles_to_geometry};
use crate::raster_stack::check_alignment;
use crate::surface_solid::{default_spacing, densify, sample_bilinear, strip};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct FenceDiagramTool;

impl Tool for FenceDiagramTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "fence_diagram",
            display_name: "Fence Diagram",
            summary: "Generates 3D vertical section panels where a line trace crosses a stack of surfaces — the standard way to visualise subsurface stratigraphy or atmospheric layering (ArcGIS Fence Diagram). profile, long_profile and stack_profile sample surfaces along a line but return 2D profile lines or tables, with no geometry for the space BETWEEN consecutive surfaces, and voxel_isosurface extracts a level surface from a volume rather than a section.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Line features defining the section trace.",
                    required: true,
                },
                ToolParamSpec {
                    name: "surfaces",
                    description: "Comma-separated list of co-registered surface rasters, ordered TOP to BOTTOM.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "3D section panels, one per trace per surface pair. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "sample_distance",
                    description: "Along-trace sampling spacing in CRS units. Default: the finest surface cell size.",
                    required: false,
                },
                ToolParamSpec {
                    name: "floor_height",
                    description: "Optional constant lower bound: a synthetic bottom surface at this elevation, closing the deepest panel.",
                    required: false,
                },
                ToolParamSpec {
                    name: "ceiling_height",
                    description: "Optional constant upper bound: a synthetic top surface at this elevation.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read from each surface raster (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        let surfaces = req_str(args, "surfaces")?;
        let n = surfaces.split(',').filter(|s| !s.trim().is_empty()).count();
        let bounds = [
            opt_f64(args, "floor_height")?.is_some(),
            opt_f64(args, "ceiling_height")?.is_some(),
        ]
        .iter()
        .filter(|b| **b)
        .count();
        // One surface plus one constant bound still defines a panel; a lone
        // surface on its own defines nothing to draw between.
        if n + bounds < 2 {
            return Err(ToolError::Validation(
                "a fence needs at least two bounding levels: supply two surfaces, or one surface \
                 plus 'floor_height' or 'ceiling_height'"
                    .to_string(),
            ));
        }
        opt_positive_f64(args, "sample_distance")?;
        band_index(args, "band")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let band = band_index(args, "band")?;
        let ceiling = opt_f64(args, "ceiling_height")?;
        let floor = opt_f64(args, "floor_height")?;

        let paths: Vec<&str> = req_str(args, "surfaces")?
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .collect();
        let rasters: Vec<Raster> = paths
            .iter()
            .map(|p| load_input_raster(p))
            .collect::<Result<_, _>>()?;
        if rasters.len() > 1 {
            check_alignment(&rasters)?;
        }
        let spacing = match opt_positive_f64(args, "sample_distance")? {
            Some(v) => v,
            None => default_spacing(&rasters.iter().collect::<Vec<_>>())?,
        };

        // Levels, top to bottom: an optional constant ceiling, the supplied
        // surfaces in order, then an optional constant floor.
        let mut levels: Vec<Level> = Vec::new();
        if let Some(z) = ceiling {
            levels.push(Level::Constant(z));
        }
        for (i, _) in rasters.iter().enumerate() {
            levels.push(Level::Surface(i));
        }
        if let Some(z) = floor {
            levels.push(Level::Constant(z));
        }

        let lines = load_input_layer(input)?;
        let mut out = Layer::new("fence_diagram");
        out.geom_type = Some(GeometryType::MultiPolygon);
        out.crs = lines.crs.clone();
        out.add_field(FieldDef::new("SRC_FID", FieldType::Integer));
        // LEVEL indices, not surface indices: a ceiling_height shifts every
        // surface by one, so these count bounding levels top to bottom.
        out.add_field(FieldDef::new("LEVEL_TOP", FieldType::Integer));
        out.add_field(FieldDef::new("LEVEL_BOTTOM", FieldType::Integer));
        out.add_field(FieldDef::new("PANEL_AREA", FieldType::Float));
        out.add_field(FieldDef::new("STATIONS", FieldType::Integer));
        out.add_field(FieldDef::new("INVERTED", FieldType::Boolean));

        let mut panels = 0_u64;
        let mut inverted_panels = 0_u64;
        let total = lines.iter().count().max(1);

        for (fid, feature) in lines.iter().enumerate() {
            for coords in line_parts(feature.geometry.as_ref()) {
                let stations = densify(&coords, spacing);
                if stations.len() < 2 {
                    continue;
                }

                // Sample every level at every station once, reusing the columns
                // across all panel pairs.
                let columns: Vec<Vec<Option<f64>>> = levels
                    .iter()
                    .map(|lv| {
                        stations
                            .iter()
                            .map(|(x, y)| match lv {
                                Level::Constant(z) => Some(*z),
                                Level::Surface(i) => sample_bilinear(&rasters[*i], band, *x, *y),
                            })
                            .collect()
                    })
                    .collect();

                for pair in 0..levels.len().saturating_sub(1) {
                    // Split into CONTIGUOUS runs of stations where both levels are
                    // defined. Simply skipping undefined stations would bridge an
                    // interior hole in a surface, drawing a panel across ground the
                    // data does not cover.
                    let mut runs: Vec<(Vec<[f64; 3]>, Vec<[f64; 3]>, bool)> = Vec::new();
                    let mut upper: Vec<[f64; 3]> = Vec::new();
                    let mut lower: Vec<[f64; 3]> = Vec::new();
                    let mut inverted = false;
                    for (s_i, (x, y)) in stations.iter().enumerate() {
                        match (columns[pair][s_i], columns[pair + 1][s_i]) {
                            (Some(zu), Some(zl)) => {
                                if zl > zu {
                                    inverted = true;
                                }
                                upper.push([*x, *y, zu]);
                                lower.push([*x, *y, zl]);
                            }
                            _ => {
                                if upper.len() >= 2 {
                                    runs.push((
                                        std::mem::take(&mut lower),
                                        std::mem::take(&mut upper),
                                        inverted,
                                    ));
                                } else {
                                    upper.clear();
                                    lower.clear();
                                }
                                inverted = false;
                            }
                        }
                    }
                    if upper.len() >= 2 {
                        runs.push((lower, upper, inverted));
                    }

                    for (lower, upper, inverted) in runs {
                        let tris = strip(&lower, &upper, false);
                        if tris.is_empty() {
                            continue;
                        }
                        let area: f64 = tris.iter().map(tri_area).sum();
                        if area <= 0.0 {
                            // Two coincident levels give a degenerate ribbon with
                            // no section to show.
                            continue;
                        }

                        out.add_feature(
                            Some(triangles_to_geometry(&tris)),
                            &[
                                ("SRC_FID", FieldValue::Integer(fid as i64)),
                                ("LEVEL_TOP", FieldValue::Integer(pair as i64)),
                                ("LEVEL_BOTTOM", FieldValue::Integer(pair as i64 + 1)),
                                ("PANEL_AREA", FieldValue::Float(area)),
                                ("STATIONS", FieldValue::Integer(upper.len() as i64)),
                                ("INVERTED", FieldValue::Boolean(inverted)),
                            ],
                        )
                        .map_err(|e| ToolError::Execution(e.to_string()))?;
                        panels += 1;
                        if inverted {
                            inverted_panels += 1;
                        }
                    }
                }
            }
            ctx.progress.progress((fid as f64 + 1.0) / total as f64);
        }

        if panels == 0 {
            return Err(ToolError::Execution(
                "no section panel could be built: check that the traces fall within the surfaces"
                    .to_string(),
            ));
        }
        ctx.progress.info(&format!(
            "{panels} panel(s) across {} level(s), {inverted_panels} inverted",
            levels.len()
        ));

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("panel_count".to_string(), json!(panels));
        outputs.insert("level_count".to_string(), json!(levels.len()));
        outputs.insert("inverted_panels".to_string(), json!(inverted_panels));
        outputs.insert("sample_distance".to_string(), json!(spacing));
        Ok(ToolRunResult { outputs })
    }
}

/// One bounding level of the fence: a sampled surface or a constant elevation.
enum Level {
    Surface(usize),
    Constant(f64),
}

/// Every part of a line feature, so a MultiLineString trace is sectioned in
/// full rather than only along its first part.
fn line_parts(geom: Option<&Geometry>) -> Vec<Vec<Coord>> {
    match geom {
        Some(Geometry::LineString(cs)) => vec![cs.clone()],
        Some(Geometry::MultiLineString(parts)) => parts.clone(),
        _ => Vec::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, RasterConfig};
    use wbvector::memory_store;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// 10x10 raster of 10-unit cells covering (0,0)-(100,100).
    fn surface(f: impl Fn(usize, usize) -> f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols: 10,
            rows: 10,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 10.0,
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
        for row in 0..10 {
            for col in 0..10 {
                r.set(0, row as isize, col as isize, f(row, col)).unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn flat(v: f64) -> String {
        surface(move |_, _| v)
    }

    fn trace(pts: Vec<(f64, f64)>) -> String {
        let mut l = Layer::new("trace");
        l.geom_type = Some(GeometryType::LineString);
        l.add_feature(
            Some(Geometry::LineString(
                pts.into_iter().map(|(x, y)| Coord::xy(x, y)).collect(),
            )),
            &[],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = FenceDiagramTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    fn num(layer: &Layer, fid: usize, name: &str) -> f64 {
        let i = layer.schema.field_index(name).unwrap();
        match &layer.iter().nth(fid).unwrap().attributes[i] {
            FieldValue::Float(v) => *v,
            FieldValue::Integer(v) => *v as f64,
            other => panic!("expected a number, got {other:?}"),
        }
    }

    #[test]
    fn a_panel_area_is_trace_length_times_layer_thickness() {
        // An 80-unit trace between flat surfaces 20 apart: area 1600.
        let (out, res) = run(json!({
            "input": trace(vec![(10.0, 50.0), (90.0, 50.0)]),
            "surfaces": format!("{},{}", flat(30.0), flat(10.0)),
            "sample_distance": 5.0,
        }));
        assert_eq!(res.outputs["panel_count"], json!(1));
        let a = num(&out, 0, "PANEL_AREA");
        assert!((a - 1600.0).abs() < 1.0, "expected 1600, got {a}");
    }

    #[test]
    fn three_surfaces_give_two_stacked_panels() {
        let (out, res) = run(json!({
            "input": trace(vec![(10.0, 50.0), (90.0, 50.0)]),
            "surfaces": format!("{},{},{}", flat(30.0), flat(20.0), flat(10.0)),
            "sample_distance": 10.0,
        }));
        assert_eq!(res.outputs["panel_count"], json!(2));
        assert_eq!(num(&out, 0, "LEVEL_TOP"), 0.0);
        assert_eq!(num(&out, 0, "LEVEL_BOTTOM"), 1.0);
        assert_eq!(num(&out, 1, "LEVEL_TOP"), 1.0);
        // Each 10-unit-thick layer over 80 units is 800.
        assert!((num(&out, 0, "PANEL_AREA") - 800.0).abs() < 1.0);
        assert!((num(&out, 1, "PANEL_AREA") - 800.0).abs() < 1.0);
    }

    #[test]
    fn the_panel_is_a_3d_ribbon_carrying_both_surface_elevations() {
        let (out, _) = run(json!({
            "input": trace(vec![(10.0, 50.0), (90.0, 50.0)]),
            "surfaces": format!("{},{}", flat(30.0), flat(10.0)),
            "sample_distance": 20.0,
        }));
        let geom = out.iter().next().unwrap().geometry.clone().unwrap();
        let tris = crate::inside_3d::collect_triangles(&geom);
        let zs: Vec<f64> = tris.iter().flatten().map(|v| v[2]).collect();
        assert!(zs.iter().any(|z| (z - 30.0).abs() < 1e-6), "no top Z");
        assert!(zs.iter().any(|z| (z - 10.0).abs() < 1e-6), "no bottom Z");
    }

    #[test]
    fn a_floor_height_closes_the_deepest_panel() {
        // One surface plus a floor is a valid fence.
        let (_, res) = run(json!({
            "input": trace(vec![(10.0, 50.0), (90.0, 50.0)]),
            "surfaces": flat(30.0),
            "floor_height": 0.0,
            "sample_distance": 10.0,
        }));
        assert_eq!(res.outputs["panel_count"], json!(1));
        assert_eq!(res.outputs["level_count"], json!(2));
    }

    #[test]
    fn a_ceiling_adds_a_panel_above_the_top_surface() {
        let (_, res) = run(json!({
            "input": trace(vec![(10.0, 50.0), (90.0, 50.0)]),
            "surfaces": flat(30.0),
            "ceiling_height": 50.0,
            "floor_height": 0.0,
            "sample_distance": 10.0,
        }));
        assert_eq!(res.outputs["level_count"], json!(3));
        assert_eq!(res.outputs["panel_count"], json!(2));
    }

    #[test]
    fn crossing_surfaces_are_flagged_rather_than_dropped() {
        // The lower surface rises above the upper one across the trace — a
        // real feature of interpolated geology, not an error to hide.
        let upper = surface(|_, col| col as f64 * 5.0);
        let lower = surface(|_, col| 40.0 - col as f64 * 3.0);
        let (out, res) = run(json!({
            "input": trace(vec![(5.0, 50.0), (95.0, 50.0)]),
            "surfaces": format!("{upper},{lower}"),
            "sample_distance": 5.0,
        }));
        assert_eq!(res.outputs["panel_count"], json!(1));
        assert_eq!(res.outputs["inverted_panels"], json!(1));
        let i = out.schema.field_index("INVERTED").unwrap();
        assert_eq!(
            out.iter().next().unwrap().attributes[i],
            FieldValue::Boolean(true)
        );
    }

    #[test]
    fn coincident_surfaces_produce_no_panel() {
        let args: ToolArgs = serde_json::from_value(json!({
            "input": trace(vec![(10.0, 50.0), (90.0, 50.0)]),
            "surfaces": format!("{},{}", flat(20.0), flat(20.0)),
        }))
        .unwrap();
        // A zero-thickness layer has no section to draw.
        assert!(FenceDiagramTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn a_finer_sample_distance_does_not_change_the_area_of_a_flat_layer() {
        // Sampling density must not bias the measurement.
        let mk = |s: f64| {
            let (out, _) = run(json!({
                "input": trace(vec![(10.0, 50.0), (90.0, 50.0)]),
                "surfaces": format!("{},{}", flat(30.0), flat(10.0)),
                "sample_distance": s,
            }));
            num(&out, 0, "PANEL_AREA")
        };
        assert!((mk(2.0) - mk(40.0)).abs() < 1.0);
    }

    #[test]
    fn rejects_bad_parameters() {
        let t = trace(vec![(10.0, 50.0), (90.0, 50.0)]);
        let s = flat(10.0);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            FenceDiagramTool.validate(&args).is_err()
        };
        assert!(bad(json!({"surfaces": s})));
        // One surface and no constant bound is not a fence.
        assert!(bad(json!({"input": t, "surfaces": s})));
        assert!(bad(
            json!({"input": t, "surfaces": format!("{s},{s}"), "sample_distance": 0})
        ));
    }
}
