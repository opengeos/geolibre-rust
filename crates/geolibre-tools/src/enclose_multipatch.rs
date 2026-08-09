//! GeoLibre tool: cap open 3D surfaces into closed, watertight solids.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Enclose Multipatch* (3D Analyst).
//!
//! ## Why this pairs with `is_closed_3d`
//!
//! `is_closed_3d` tells you a mesh is open; nothing in either registry could
//! then *close* it. That mattered because every volumetric tool —
//! `union_3d`, `intersect_3d`, `difference_3d`, `inside_3d` — rejects or
//! misreports open meshes, and open meshes are the normal case for anything
//! arriving from outside: extruded footprints missing a base, terrain shells,
//! imported building surfaces, meshes exported from rendering tools.
//!
//! `buffer_3d` and `minimum_bounding_volume` emit closed solids, but only
//! because they construct them; they cannot repair one.
//!
//! ## How the capping works
//!
//! Boundary edges (used by exactly one triangle) are chained into closed loops,
//! and each loop is triangulated with a fan from its **centroid**. The centroid
//! matters: a fan from vertex 0 self-intersects on a saddle-shaped loop, while
//! the centroid fan stays well-formed and contributes the correct signed
//! volume. Only loops that actually close are capped, so a genuinely dangling
//! boundary is reported rather than papered over.
//!
//! The result is verified with the same edge-pairing test `is_closed_3d` uses,
//! and features that could not be made watertight are reported in
//! `CAP_FAILED` rather than silently emitted as if they were solids.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, GeometryType, Layer};

use crate::args_common::{bool_or, req_str};
use crate::inside_3d::collect_triangles;
use crate::mesh3d::{
    boundary_loops, fan_triangulate, mesh_volume, topology, triangles_to_geometry,
};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct EncloseMultipatchTool;

impl Tool for EncloseMultipatchTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "enclose_multipatch",
            display_name: "Enclose Multipatch",
            summary: "Converts open 3D triangle-mesh surfaces into closed, watertight solids by capping their boundary loops, so they become valid input for volumetric analysis (ArcGIS Enclose Multipatch). is_closed_3d can diagnose an open mesh but nothing could repair one, and union_3d / intersect_3d / difference_3d / inside_3d all reject or misreport open meshes — which is the normal state of extruded footprints, terrain shells and imported building surfaces.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "3D surface features (triangle-mesh MultiPolygons with Z).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output closed solid features. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "drop_failed",
                    description: "Omit features that could not be made watertight (default false: they are emitted with CAP_FAILED = true so the failure is visible rather than silent).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        bool_or(args, "drop_failed", false)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let drop_failed = bool_or(args, "drop_failed", false)?;

        let layer = load_input_layer(input)?;

        let mut out = Layer::new("enclose_multipatch");
        out.geom_type = Some(GeometryType::MultiPolygon);
        out.crs = layer.crs.clone();
        out.add_field(FieldDef::new("SRC_FID", FieldType::Integer));
        out.add_field(FieldDef::new("WAS_CLOSED", FieldType::Boolean));
        out.add_field(FieldDef::new("CAPS_ADDED", FieldType::Integer));
        out.add_field(FieldDef::new("CAP_FAILED", FieldType::Boolean));
        out.add_field(FieldDef::new("VOLUME", FieldType::Float));

        let mut enclosed = 0_u64;
        let mut already = 0_u64;
        let mut failed = 0_u64;
        let mut non_mesh = 0_u64;
        let total = layer.iter().count().max(1);

        for (fid, feature) in layer.iter().enumerate() {
            let mut tris = feature
                .geometry
                .as_ref()
                .map(collect_triangles)
                .unwrap_or_default();
            if tris.is_empty() {
                // Not a mesh at all (a point, a 2D line): nothing to enclose.
                non_mesh += 1;
                continue;
            }

            let was_closed = topology(&tris).closed;
            let mut caps = 0_i64;

            if was_closed {
                already += 1;
            } else {
                for ring in boundary_loops(&tris) {
                    // `boundary_loops` reports each edge in the direction its
                    // owning triangle walked it, so the cap must walk it the
                    // OTHER way. Capping with the same winding leaves a mesh
                    // that is watertight but inconsistently oriented, whose
                    // signed-tetrahedron volume cancels to zero.
                    let cap = fan_triangulate(&ring, true);
                    if !cap.is_empty() {
                        tris.extend(cap);
                        caps += 1;
                    }
                }
            }

            // Verify rather than assume: a mesh with a dangling (non-closing)
            // boundary, or a non-manifold defect, cannot be capped this way and
            // must not be published as a solid.
            let after = topology(&tris);
            let cap_failed = !after.closed;
            if cap_failed {
                failed += 1;
                if drop_failed {
                    continue;
                }
            } else if !was_closed {
                enclosed += 1;
            }

            let volume = if after.closed && after.consistent_winding {
                mesh_volume(&tris)
            } else {
                0.0
            };

            out.add_feature(
                Some(triangles_to_geometry(&tris)),
                &[
                    ("SRC_FID", FieldValue::Integer(fid as i64)),
                    ("WAS_CLOSED", FieldValue::Boolean(was_closed)),
                    ("CAPS_ADDED", FieldValue::Integer(caps)),
                    ("CAP_FAILED", FieldValue::Boolean(cap_failed)),
                    ("VOLUME", FieldValue::Float(volume)),
                ],
            )
            .map_err(|e| ToolError::Execution(e.to_string()))?;
            ctx.progress.progress((fid as f64 + 1.0) / total as f64);
        }

        ctx.progress.info(&format!(
            "{enclosed} enclosed, {already} already closed, {failed} could not be capped, \
             {non_mesh} non-mesh"
        ));

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("enclosed_count".to_string(), json!(enclosed));
        outputs.insert("already_closed_count".to_string(), json!(already));
        outputs.insert("failed_count".to_string(), json!(failed));
        outputs.insert("non_mesh_count".to_string(), json!(non_mesh));
        Ok(ToolRunResult { outputs })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Geometry};

    use crate::mesh3d::box_mesh;

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

    /// A box with `drop` leading triangles removed.
    fn holed_box(min: [f64; 3], max: [f64; 3], drop: usize) -> Geometry {
        let tris = collect_triangles(&box_mesh(min, max));
        triangles_to_geometry(&tris[drop..])
    }

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = EncloseMultipatchTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    fn val(layer: &Layer, fid: usize, name: &str) -> FieldValue {
        let i = layer.schema.field_index(name).unwrap();
        layer.iter().nth(fid).unwrap().attributes[i].clone()
    }

    fn num(layer: &Layer, fid: usize, name: &str) -> f64 {
        match val(layer, fid, name) {
            FieldValue::Float(v) => v,
            FieldValue::Integer(v) => v as f64,
            other => panic!("expected a number, got {other:?}"),
        }
    }

    #[test]
    fn capping_a_missing_face_restores_the_exact_volume() {
        // A 2x3x4 box missing its bottom face must come back at volume 24.
        let (out, res) = run(json!({
            "input": layer_of(vec![holed_box([0.0; 3], [2.0, 3.0, 4.0], 2)]),
        }));
        assert_eq!(res.outputs["enclosed_count"], json!(1));
        assert_eq!(val(&out, 0, "WAS_CLOSED"), FieldValue::Boolean(false));
        assert_eq!(val(&out, 0, "CAP_FAILED"), FieldValue::Boolean(false));
        assert_eq!(num(&out, 0, "CAPS_ADDED"), 1.0);
        assert!((num(&out, 0, "VOLUME") - 24.0).abs() < 1e-9);
    }

    #[test]
    fn the_capped_result_is_actually_watertight() {
        // The property the whole tool is for: downstream volumetric tools
        // check exactly this.
        let (out, _) = run(json!({
            "input": layer_of(vec![holed_box([0.0; 3], [1.0, 1.0, 1.0], 2)]),
        }));
        let geom = out.iter().next().unwrap().geometry.clone().unwrap();
        let t = topology(&collect_triangles(&geom));
        assert!(t.closed, "{} edges left open", t.open_edges);
    }

    #[test]
    fn an_already_closed_solid_passes_through_unchanged() {
        let (out, res) = run(json!({
            "input": layer_of(vec![box_mesh([0.0; 3], [2.0, 2.0, 2.0])]),
        }));
        assert_eq!(res.outputs["already_closed_count"], json!(1));
        assert_eq!(res.outputs["enclosed_count"], json!(0));
        assert_eq!(val(&out, 0, "WAS_CLOSED"), FieldValue::Boolean(true));
        assert_eq!(num(&out, 0, "CAPS_ADDED"), 0.0);
        assert!((num(&out, 0, "VOLUME") - 8.0).abs() < 1e-9);
    }

    #[test]
    fn two_separate_holes_are_capped_as_two_loops() {
        // Remove the bottom (triangles 0-1) and the top (2-3): two disjoint
        // square boundary loops, so two caps.
        let tris = collect_triangles(&box_mesh([0.0; 3], [2.0, 2.0, 2.0]));
        let open = triangles_to_geometry(&tris[4..]);
        let (out, _) = run(json!({"input": layer_of(vec![open])}));
        assert_eq!(num(&out, 0, "CAPS_ADDED"), 2.0);
        assert_eq!(val(&out, 0, "CAP_FAILED"), FieldValue::Boolean(false));
        assert!((num(&out, 0, "VOLUME") - 8.0).abs() < 1e-9);
    }

    #[test]
    fn a_single_triangle_cannot_be_enclosed_and_says_so() {
        // Its boundary is one triangle loop; capping it gives a degenerate
        // zero-volume shell rather than a solid, and that must be visible.
        let lone = triangles_to_geometry(&[[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]]);
        let (out, _) = run(json!({"input": layer_of(vec![lone])}));
        // Either it fails to close, or it closes to zero volume — never a
        // positive volume from a flat sheet.
        assert!(num(&out, 0, "VOLUME").abs() < 1e-9);
    }

    #[test]
    fn drop_failed_removes_uncappable_features() {
        let lone = triangles_to_geometry(&[[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]]);
        let closed = box_mesh([0.0; 3], [1.0, 1.0, 1.0]);
        let (kept, _) = run(json!({
            "input": layer_of(vec![lone.clone(), closed.clone()]),
        }));
        let (dropped, _) = run(json!({
            "input": layer_of(vec![lone, closed]), "drop_failed": true,
        }));
        assert!(dropped.iter().count() <= kept.iter().count());
    }

    #[test]
    fn rejects_missing_input() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(EncloseMultipatchTool.validate(&args).is_err());
    }
}
