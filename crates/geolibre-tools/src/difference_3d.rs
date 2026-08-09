//! GeoLibre tool: subtract 3D solids from other 3D solids.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Difference 3D* (3D Analyst).
//!
//! ## The gap
//!
//! Same root cause as `intersect_3d`: `union_3d` was the only 3D overlay in
//! either registry, and `erase` is polygon-only 2D. Without a 3D difference
//! there is no way to answer "how much of this solid survives once the cut
//! volume is removed" — excavation remaining after a cut body, plume volume
//! outside an exclusion envelope, buildable volume left after setback solids.
//!
//! ## Scope, deliberately (inherited from `union_3d`)
//!
//! Each minuend's own volume is **exact** (signed-tetrahedron summation). The
//! *removed* volume is estimated by voxel occupancy over the minuend's bounding
//! box at a reported `resolution`, and the remainder is the exact volume minus
//! that estimate. No exact merged mesh is produced or claimed.
//!
//! ## Why the union of subtrahends, not a running subtraction
//!
//! Subtracting each subtrahend's volume in turn double-counts wherever two
//! subtrahends overlap each other, which would report a *negative* remainder
//! for heavily overlapping cut bodies. Occupancy against "inside the minuend
//! and inside **any** subtrahend" is overlap-correct by construction — the same
//! reason `union_3d` samples rather than sums.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Layer};

use crate::args_common::{req_str, usize_or};
use crate::intersect_3d::load_solids;
use crate::mesh3d::{bbox_overlap, mesh_volume, occupancy_volume};
use crate::vector_common::{parse_optional_str, write_or_store_layer};

pub struct Difference3dTool;

impl Tool for Difference3dTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "difference_3d",
            display_name: "Difference 3D",
            summary: "Subtracts a set of 3D solids from another set, reporting the volume removed from and remaining in each minuend feature (ArcGIS Difference 3D). union_3d was the only 3D overlay in either registry and erase is polygon-only 2D, so there was no way to compute excavated volume remaining after a cut body or plume volume outside an exclusion envelope. Minuend volumes are exact; the removed volume is estimated by voxel occupancy at a reported resolution.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Minuend closed 3D solids (triangle-mesh MultiPolygons with Z).",
                    required: true,
                },
                ToolParamSpec {
                    name: "subtract",
                    description: "Subtrahend closed 3D solids to remove from each minuend.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output table, one row per minuend. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "resolution",
                    description: "Voxel cells along the longest axis of each minuend's bounding box (default 64). Higher is more accurate and slower.",
                    required: false,
                },
                ToolParamSpec {
                    name: "keep_geometry",
                    description: "Carry each minuend's original geometry onto its output row (default false: an attribute-only table). The geometry is the UNCUT minuend — no exact merged mesh is produced.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        req_str(args, "subtract")?;
        crate::args_common::opt_bool(args, "keep_geometry")?;
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
        let subtract = req_str(args, "subtract")?;
        let output = parse_optional_str(args, "output")?;
        let resolution = usize_or(args, "resolution", 64)?;
        let keep_geometry = crate::args_common::bool_or(args, "keep_geometry", false)?;

        let (minuends, min_open) = load_solids(input)?;
        let (subtrahends, sub_open) = load_solids(subtract)?;
        if minuends.is_empty() {
            return Err(ToolError::Execution(format!(
                "'input' holds no closed triangle-mesh solids ({min_open} open mesh(es) skipped)"
            )));
        }

        // An empty subtrahend set is not an error: nothing is removed, and the
        // remainder equals the input. Reporting that is more useful than
        // failing, since it is the natural result of a filtered cut layer.
        ctx.progress.info(&format!(
            "{} minuend(s), {} subtrahend(s)",
            minuends.len(),
            subtrahends.len()
        ));

        let mut out = Layer::new("difference_3d");
        if keep_geometry {
            out.geom_type = Some(wbvector::GeometryType::MultiPolygon);
        }
        out.add_field(FieldDef::new("SRC_FID", FieldType::Integer));
        out.add_field(FieldDef::new("VOLUME_IN", FieldType::Float));
        out.add_field(FieldDef::new("VOLUME_REMOVED", FieldType::Float));
        out.add_field(FieldDef::new("VOLUME_OUT", FieldType::Float));
        out.add_field(FieldDef::new("CUTTERS", FieldType::Integer));
        out.add_field(FieldDef::new("RESOLUTION", FieldType::Integer));

        let mut total_removed = 0.0_f64;
        let mut fully_consumed = 0_u64;
        let total = minuends.len().max(1);

        for (n, m) in minuends.iter().enumerate() {
            let volume_in = mesh_volume(&m.tris);
            // Only subtrahends whose boxes meet this minuend can remove
            // anything, and restricting the set also shrinks the inner test.
            let cutters: Vec<&crate::inside_3d::Solid> = subtrahends
                .iter()
                .filter(|s| bbox_overlap(m, s))
                .collect();

            let removed = if cutters.is_empty() {
                0.0
            } else {
                occupancy_volume(m.min, m.max, resolution, |x, y, z| {
                    m.contains(x, y, z) && cutters.iter().any(|s| s.contains(x, y, z))
                })
            };
            // The estimate can marginally exceed the exact volume at coarse
            // resolutions; clamping keeps VOLUME_OUT physically meaningful
            // instead of letting a sampling artefact go negative.
            let removed = removed.min(volume_in);
            let remaining = (volume_in - removed).max(0.0);

            if remaining <= 0.0 && volume_in > 0.0 {
                fully_consumed += 1;
            }
            total_removed += removed;

            let geom = keep_geometry.then(|| crate::mesh3d::triangles_to_geometry(&m.tris));
            out.add_feature(
                geom,
                &[
                    ("SRC_FID", FieldValue::Integer(m.fid as i64)),
                    ("VOLUME_IN", FieldValue::Float(volume_in)),
                    ("VOLUME_REMOVED", FieldValue::Float(removed)),
                    ("VOLUME_OUT", FieldValue::Float(remaining)),
                    ("CUTTERS", FieldValue::Integer(cutters.len() as i64)),
                    ("RESOLUTION", FieldValue::Integer(resolution as i64)),
                ],
            )
            .map_err(|e| ToolError::Execution(e.to_string()))?;
            ctx.progress.progress((n as f64 + 1.0) / total as f64);
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("minuend_count".to_string(), json!(minuends.len()));
        outputs.insert("subtrahend_count".to_string(), json!(subtrahends.len()));
        outputs.insert("total_volume_removed".to_string(), json!(total_removed));
        outputs.insert("fully_consumed".to_string(), json!(fully_consumed));
        outputs.insert("resolution".to_string(), json!(resolution));
        outputs.insert(
            "open_meshes_skipped".to_string(),
            json!(min_open + sub_open),
        );
        Ok(ToolRunResult { outputs })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Geometry, GeometryType};

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

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = Difference3dTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer_local(res.outputs["output"].as_str().unwrap());
        (layer, res)
    }

    fn load_input_layer_local(p: &str) -> Layer {
        crate::vector_common::load_input_layer(p).unwrap()
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
    fn a_half_overlapping_cutter_removes_half_the_volume() {
        // Minuend [0,10]^3 = 1000. Cutter covers x in [5,15]: removes 500.
        let (out, _) = run(json!({
            "input": layer_of(vec![box_mesh([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])]),
            "subtract": layer_of(vec![box_mesh([5.0, 0.0, 0.0], [15.0, 10.0, 10.0])]),
            "resolution": 64,
        }));
        assert!((num(&out, 0, "VOLUME_IN") - 1000.0).abs() < 1e-6);
        let removed = num(&out, 0, "VOLUME_REMOVED");
        assert!((removed - 500.0).abs() < 15.0, "removed {removed}");
        // The two halves must still add up.
        assert!(
            (num(&out, 0, "VOLUME_OUT") + removed - 1000.0).abs() < 1e-6,
            "in/out/removed are inconsistent"
        );
    }

    #[test]
    fn a_disjoint_cutter_removes_nothing() {
        let (out, res) = run(json!({
            "input": layer_of(vec![box_mesh([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])]),
            "subtract": layer_of(vec![box_mesh([50.0, 50.0, 50.0], [60.0, 60.0, 60.0])]),
        }));
        assert_eq!(num(&out, 0, "VOLUME_REMOVED"), 0.0);
        assert_eq!(num(&out, 0, "CUTTERS"), 0.0);
        assert!((num(&out, 0, "VOLUME_OUT") - 1000.0).abs() < 1e-6);
        assert_eq!(res.outputs["fully_consumed"], json!(0));
    }

    #[test]
    fn an_enclosing_cutter_consumes_the_minuend_entirely() {
        let (out, res) = run(json!({
            "input": layer_of(vec![box_mesh([4.0, 4.0, 4.0], [6.0, 6.0, 6.0])]),
            "subtract": layer_of(vec![box_mesh([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])]),
            "resolution": 64,
        }));
        assert!(num(&out, 0, "VOLUME_OUT") < 1e-6);
        assert_eq!(res.outputs["fully_consumed"], json!(1));
    }

    #[test]
    fn overlapping_cutters_are_not_double_counted() {
        // Two cutters each covering x in [5,15] and [6,16] overlap heavily.
        // Subtracting them in turn would remove more than the minuend holds
        // and drive the remainder negative; occupancy must not.
        let (out, _) = run(json!({
            "input": layer_of(vec![box_mesh([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])]),
            "subtract": layer_of(vec![
                box_mesh([5.0, 0.0, 0.0], [15.0, 10.0, 10.0]),
                box_mesh([6.0, 0.0, 0.0], [16.0, 10.0, 10.0]),
            ]),
            "resolution": 64,
        }));
        // The union of the cutters covers x in [5,10] inside the minuend: 500.
        let removed = num(&out, 0, "VOLUME_REMOVED");
        assert!((removed - 500.0).abs() < 20.0, "removed {removed}");
        assert!(num(&out, 0, "VOLUME_OUT") >= 0.0);
    }

    #[test]
    fn removed_volume_never_exceeds_the_minuend() {
        // The clamp that keeps VOLUME_OUT physically meaningful at coarse
        // resolutions.
        let (out, _) = run(json!({
            "input": layer_of(vec![box_mesh([0.0, 0.0, 0.0], [3.0, 3.0, 3.0])]),
            "subtract": layer_of(vec![box_mesh([-5.0, -5.0, -5.0], [8.0, 8.0, 8.0])]),
            "resolution": 4,
        }));
        assert!(num(&out, 0, "VOLUME_REMOVED") <= num(&out, 0, "VOLUME_IN") + 1e-9);
        assert!(num(&out, 0, "VOLUME_OUT") >= 0.0);
    }

    #[test]
    fn every_minuend_gets_a_row_even_when_untouched() {
        let (out, res) = run(json!({
            "input": layer_of(vec![
                box_mesh([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]),
                box_mesh([90.0, 90.0, 90.0], [92.0, 92.0, 92.0]),
            ]),
            "subtract": layer_of(vec![box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])]),
        }));
        assert_eq!(out.iter().count(), 2);
        assert_eq!(res.outputs["minuend_count"], json!(2));
        assert!(num(&out, 1, "VOLUME_REMOVED") == 0.0);
    }

    #[test]
    fn higher_resolution_converges() {
        let mk = |r: usize| {
            let (out, _) = run(json!({
                "input": layer_of(vec![box_mesh([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])]),
                "subtract": layer_of(vec![box_mesh([7.0, 0.0, 0.0], [17.0, 10.0, 10.0])]),
                "resolution": r,
            }));
            (num(&out, 0, "VOLUME_REMOVED") - 300.0).abs()
        };
        let coarse = mk(6);
        let fine = mk(96);
        assert!(fine <= coarse, "coarse err {coarse}, fine err {fine}");
        assert!(fine < 10.0, "fine estimate off by {fine}");
    }

    #[test]
    fn keep_geometry_carries_the_uncut_minuend() {
        let (out, _) = run(json!({
            "input": layer_of(vec![box_mesh([0.0, 0.0, 0.0], [2.0, 2.0, 2.0])]),
            "subtract": layer_of(vec![box_mesh([1.0, 0.0, 0.0], [3.0, 2.0, 2.0])]),
            "keep_geometry": true,
        }));
        assert!(out.iter().next().unwrap().geometry.is_some());
    }

    #[test]
    fn rejects_bad_parameters() {
        let path = layer_of(vec![box_mesh([0.0; 3], [1.0, 1.0, 1.0])]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            Difference3dTool.validate(&args).is_err()
        };
        assert!(bad(json!({"input": path})));
        assert!(bad(json!({"subtract": path})));
        assert!(bad(json!({"input": path, "subtract": path, "resolution": 0})));
    }
}
