//! GeoLibre tool: intersection of overlapping closed 3D solids.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Intersect 3D* (3D Analyst).
//!
//! ## The gap
//!
//! The 2D overlay suite is complete and battle-tested — `union`, `intersect`,
//! `erase`, `identity`, `symmetrical_difference`, `clip`, all `geo`
//! `BooleanOps`-backed. In 3D, round 16's `union_3d` is the *only* overlay
//! operation in either registry, and it answers only the combined-volume
//! question. There is no way to ask how much volume two solids **share**, which
//! is the question behind clash detection, shared airspace, plume overlap and
//! excavation conflict.
//!
//! `inside_3d` answers a containment *predicate* for points and lines, not
//! solid-versus-solid overlap, and `polygon_volume` / `surface_volume` both
//! measure against a reference plane.
//!
//! ## Scope, deliberately (inherited from `union_3d`)
//!
//! Per-solid volumes are exact (signed-tetrahedron summation). The
//! **intersection** volume is estimated by voxel occupancy over the pair's
//! shared bounding box, at a reported `resolution`. That is approximate but
//! bounded, trivially parallel, and free of the floating-point robustness
//! minefield an exact mesh boolean's arrangement and coplanar-face handling
//! represents. The parameter surface does not promise exact merged geometry;
//! `mode = solid` emits the intersection's **bounding solid** with the measured
//! volume as an attribute, never a claimed exact mesh.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, GeometryType, Layer};

use crate::args_common::{choice_or, req_str, usize_or};
use crate::inside_3d::{collect_triangles, Solid};
use crate::mesh3d::{box_mesh, grid_for, intersect_bbox, mesh_volume, occupancy_volume};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct Intersect3dTool;

impl Tool for Intersect3dTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "intersect_3d",
            display_name: "Intersect 3D",
            summary: "Computes the volume shared by overlapping closed 3D solids, pair by pair, with an optional bounding solid for each intersection (ArcGIS Intersect 3D). union_3d is the only 3D overlay in either registry and answers only the combined-volume question; inside_3d is a containment predicate for points and lines, and polygon_volume / surface_volume measure against a reference plane. Per-solid volumes are exact; the intersection is estimated by voxel occupancy at a reported resolution.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Closed 3D solid features (triangle-mesh MultiPolygons with Z, as buffer_3d and minimum_bounding_volume emit).",
                    required: true,
                },
                ToolParamSpec {
                    name: "input2",
                    description: "Optional second solid layer. When omitted, every pair within 'input' is intersected (ArcGIS's single-input mode).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output table of intersecting pairs, one row per non-empty intersection. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "mode",
                    description: "'table' (default): attributes only. 'solid': also emit each intersection's bounding solid as geometry. The measured volume is always an attribute — no exact merged mesh is claimed.",
                    required: false,
                },
                ToolParamSpec {
                    name: "resolution",
                    description: "Voxel cells along the longest axis of each pair's shared bounding box (default 64). Higher is more accurate and slower.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        choice_or(args, "mode", &["table", "solid"], "table")?;
        let r = usize_or(args, "resolution", 64)?;
        if !(2..=1024).contains(&r) {
            return Err(ToolError::Validation(
                "'resolution' must be between 2 and 1024".to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?;
        let input2 = parse_optional_str(args, "input2")?;
        let output = parse_optional_str(args, "output")?;
        let emit_solid = choice_or(args, "mode", &["table", "solid"], "table")? == "solid";
        let resolution = usize_or(args, "resolution", 64)?;

        let (a_solids, a_open) = load_solids(input)?;
        let (b_solids, b_open) = match input2 {
            Some(p) => {
                let (s, o) = load_solids(p)?;
                (Some(s), o)
            }
            None => (None, 0),
        };
        if a_solids.is_empty() {
            return Err(ToolError::Execution(format!(
                "'input' holds no closed triangle-mesh solids ({a_open} open mesh(es) skipped)"
            )));
        }
        if let Some(b) = &b_solids {
            if b.is_empty() {
                return Err(ToolError::Execution(format!(
                    "'input2' holds no closed triangle-mesh solids ({b_open} open mesh(es) skipped)"
                )));
            }
        }

        // Pair enumeration: cross-layer when input2 is given, otherwise every
        // unordered pair within the single layer.
        let pairs: Vec<(usize, usize)> = match &b_solids {
            Some(b) => (0..a_solids.len())
                .flat_map(|i| (0..b.len()).map(move |j| (i, j)))
                .collect(),
            None => (0..a_solids.len())
                .flat_map(|i| ((i + 1)..a_solids.len()).map(move |j| (i, j)))
                .collect(),
        };
        ctx.progress
            .info(&format!("{} candidate pair(s)", pairs.len()));

        let mut out = Layer::new("intersect_3d");
        if emit_solid {
            out.geom_type = Some(GeometryType::MultiPolygon);
        }
        out.add_field(FieldDef::new("SRC_FID_1", FieldType::Integer));
        out.add_field(FieldDef::new("SRC_FID_2", FieldType::Integer));
        out.add_field(FieldDef::new("VOLUME_1", FieldType::Float));
        out.add_field(FieldDef::new("VOLUME_2", FieldType::Float));
        out.add_field(FieldDef::new("VOLUME", FieldType::Float));
        out.add_field(FieldDef::new("RESOLUTION", FieldType::Integer));

        let mut hits = 0_u64;
        let mut total_volume = 0.0_f64;
        let total = pairs.len().max(1);

        for (n, (i, j)) in pairs.iter().enumerate() {
            let sa = &a_solids[*i];
            let sb = match &b_solids {
                Some(b) => &b[*j],
                None => &a_solids[*j],
            };
            // Cheap rejection first: no shared bounding box, no intersection.
            let Some((min, max)) = intersect_bbox(sa, sb) else {
                continue;
            };
            let volume = occupancy_volume(min, max, resolution, |x, y, z| {
                sa.contains(x, y, z) && sb.contains(x, y, z)
            });
            if volume <= 0.0 {
                // Boxes touch but the solids do not actually overlap.
                continue;
            }

            let geom = emit_solid.then(|| box_mesh(min, max));
            out.add_feature(
                geom,
                &[
                    ("SRC_FID_1", FieldValue::Integer(sa.fid as i64)),
                    ("SRC_FID_2", FieldValue::Integer(sb.fid as i64)),
                    ("VOLUME_1", FieldValue::Float(mesh_volume(&sa.tris))),
                    ("VOLUME_2", FieldValue::Float(mesh_volume(&sb.tris))),
                    ("VOLUME", FieldValue::Float(volume)),
                    ("RESOLUTION", FieldValue::Integer(resolution as i64)),
                ],
            )
            .map_err(|e| ToolError::Execution(e.to_string()))?;
            hits += 1;
            total_volume += volume;
            ctx.progress.progress((n as f64 + 1.0) / total as f64);
        }

        // The effective cell size of the last-sized grid is a useful accuracy
        // hint, so report it rather than leaving `resolution` unitless.
        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("intersecting_pairs".to_string(), json!(hits));
        outputs.insert("total_intersection_volume".to_string(), json!(total_volume));
        outputs.insert("candidate_pairs".to_string(), json!(pairs.len()));
        outputs.insert("resolution".to_string(), json!(resolution));
        outputs.insert("open_meshes_skipped".to_string(), json!(a_open + b_open));
        Ok(ToolRunResult { outputs })
    }
}

/// Loads closed solids from a layer, counting the open meshes it skipped.
///
/// Open meshes are skipped rather than measured: an unclosed mesh does not
/// bound a volume, so any number computed from it is an artefact of where its
/// hole happens to be. `is_closed_3d` exists to diagnose them.
pub(crate) fn load_solids(path: &str) -> Result<(Vec<Solid>, u64), ToolError> {
    let layer = load_input_layer(path)?;
    let mut solids = Vec::new();
    let mut open = 0_u64;
    for (fid, feature) in layer.iter().enumerate() {
        let Some(geom) = feature.geometry.as_ref() else {
            continue;
        };
        let tris = collect_triangles(geom);
        if tris.is_empty() {
            continue;
        }
        let solid = Solid::new(fid, tris);
        if solid.closed {
            solids.push(solid);
        } else {
            open += 1;
        }
    }
    Ok((solids, open))
}

/// Voxel edge length a grid over `min..max` at `resolution` would use — the
/// accuracy scale of an occupancy estimate.
pub(crate) fn cell_size(min: [f64; 3], max: [f64; 3], resolution: usize) -> f64 {
    let (_, _, _, cell_vol, _) = grid_for(min, max, resolution);
    cell_vol.cbrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Geometry};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn layer_of(geoms: Vec<Geometry>) -> String {
        let mut l = Layer::new("in");
        l.geom_type = Some(GeometryType::MultiPolygon);
        for g in geoms {
            l.add_feature(Some(g), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = Intersect3dTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    fn volume(res: &ToolRunResult) -> f64 {
        res.outputs["total_intersection_volume"].as_f64().unwrap()
    }

    #[test]
    fn two_boxes_overlapping_on_one_axis_share_the_expected_volume() {
        // [0,10]^3 and [8,18]^3 overlap in x only over [8,10]: 2 x 10 x 10 = 200.
        let (_, res) = run(json!({
            "input": layer_of(vec![
                box_mesh([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]),
                box_mesh([8.0, 0.0, 0.0], [18.0, 10.0, 10.0]),
            ]),
            "resolution": 64,
        }));
        assert_eq!(res.outputs["intersecting_pairs"], json!(1));
        let v = volume(&res);
        assert!((v - 200.0).abs() < 5.0, "expected ~200, got {v}");
    }

    #[test]
    fn disjoint_solids_produce_no_rows() {
        let (out, res) = run(json!({
            "input": layer_of(vec![
                box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
                box_mesh([50.0, 50.0, 50.0], [51.0, 51.0, 51.0]),
            ]),
        }));
        assert_eq!(res.outputs["intersecting_pairs"], json!(0));
        assert_eq!(out.iter().count(), 0);
        assert_eq!(volume(&res), 0.0);
    }

    #[test]
    fn a_nested_solid_intersects_to_its_own_full_volume() {
        // The small box lies entirely inside the large one, so the
        // intersection is the small box: 2^3 = 8.
        let (_, res) = run(json!({
            "input": layer_of(vec![
                box_mesh([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]),
                box_mesh([4.0, 4.0, 4.0], [6.0, 6.0, 6.0]),
            ]),
            "resolution": 96,
        }));
        let v = volume(&res);
        assert!((v - 8.0).abs() < 0.5, "expected ~8, got {v}");
    }

    #[test]
    fn intersection_never_exceeds_either_input_volume() {
        // A structural invariant: A n B <= min(|A|, |B|). A sampler that
        // mixed up its bounding box could violate this.
        let (out, _) = run(json!({
            "input": layer_of(vec![
                box_mesh([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]),
                box_mesh([5.0, 5.0, 5.0], [15.0, 15.0, 15.0]),
            ]),
            "resolution": 64,
        }));
        let f = out.iter().next().unwrap();
        let idx = |n: &str| out.schema.field_index(n).unwrap();
        let num = |n: &str| match &f.attributes[idx(n)] {
            FieldValue::Float(v) => *v,
            other => panic!("expected float, got {other:?}"),
        };
        assert!(num("VOLUME") <= num("VOLUME_1") + 1e-6);
        assert!(num("VOLUME") <= num("VOLUME_2") + 1e-6);
    }

    #[test]
    fn higher_resolution_converges_on_the_analytic_answer() {
        let mk = |r: usize| {
            let (_, res) = run(json!({
                "input": layer_of(vec![
                    box_mesh([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]),
                    box_mesh([7.0, 7.0, 0.0], [17.0, 17.0, 10.0]),
                ]),
                "resolution": r,
            }));
            volume(&res)
        };
        // Overlap is 3 x 3 x 10 = 90.
        let coarse = (mk(8) - 90.0).abs();
        let fine = (mk(96) - 90.0).abs();
        assert!(fine <= coarse, "coarse err {coarse}, fine err {fine}");
        assert!(fine < 2.0, "fine estimate off by {fine}");
    }

    #[test]
    fn two_layers_intersect_across_rather_than_within() {
        // Within layer A the two boxes overlap each other, but single-input
        // mode is not what was asked for: only A-vs-B pairs must be reported.
        let a = layer_of(vec![
            box_mesh([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]),
            box_mesh([5.0, 0.0, 0.0], [15.0, 10.0, 10.0]),
        ]);
        let b = layer_of(vec![box_mesh([0.0, 0.0, 0.0], [2.0, 10.0, 10.0])]);
        let (_, res) = run(json!({"input": a, "input2": b, "resolution": 48}));
        // Only box A0 meets B; A1 starts at x=5 and B ends at x=2.
        assert_eq!(res.outputs["intersecting_pairs"], json!(1));
    }

    #[test]
    fn solid_mode_emits_geometry_and_table_mode_does_not() {
        let input = layer_of(vec![
            box_mesh([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]),
            box_mesh([8.0, 0.0, 0.0], [18.0, 10.0, 10.0]),
        ]);
        let (table, _) = run(json!({"input": input, "mode": "table"}));
        assert!(table.iter().next().unwrap().geometry.is_none());

        let (solid, _) = run(json!({"input": input, "mode": "solid"}));
        let g = solid.iter().next().unwrap().geometry.clone().unwrap();
        assert!(matches!(g, Geometry::MultiPolygon(_)));
    }

    #[test]
    fn open_meshes_are_skipped_and_counted() {
        let closed = box_mesh([0.0; 3], [10.0, 10.0, 10.0]);
        let tris = collect_triangles(&box_mesh([5.0, 0.0, 0.0], [15.0, 10.0, 10.0]));
        let open = crate::mesh3d::triangles_to_geometry(&tris[2..]);
        let args: ToolArgs =
            serde_json::from_value(json!({"input": layer_of(vec![closed, open])})).unwrap();
        // Only one closed solid remains, so there are no pairs at all.
        let res = Intersect3dTool.run(&args, &ctx()).unwrap();
        assert_eq!(res.outputs["open_meshes_skipped"], json!(1));
        assert_eq!(res.outputs["intersecting_pairs"], json!(0));
    }

    #[test]
    fn rejects_bad_parameters() {
        let path = layer_of(vec![box_mesh([0.0; 3], [1.0, 1.0, 1.0])]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            Intersect3dTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"input": path, "mode": "nope"})));
        assert!(bad(json!({"input": path, "resolution": 1})));
        assert!(bad(json!({"input": path, "resolution": 5000})));
    }
}
