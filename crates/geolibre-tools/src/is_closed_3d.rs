//! GeoLibre tool: test whether 3D features are watertight closed solids.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Is Closed 3D* (3D Analyst).
//!
//! ## Why this is a tool and not an internal check
//!
//! Every volumetric tool in the catalog — `union_3d`, `inside_3d`,
//! `minimum_bounding_volume`, `buffer_3d`, and round 17's `intersect_3d` /
//! `difference_3d` — assumes its input bounds a volume. Feed one an open
//! surface and you get a plausible-looking number that means nothing: the
//! signed-tetrahedron sum of an unclosed mesh is an artefact of where the hole
//! happens to be, and ray-cast containment flips parity arbitrarily.
//!
//! `union_3d` already skips open meshes and reports a count, but there was no
//! way to ask *which* features are open, or *how* open they are, before running
//! anything. This promotes the invariant `voxel_isosurface`'s own test asserts
//! internally into a user-facing precondition check.
//!
//! ## What "closed" means here
//!
//! A triangle mesh bounds a volume exactly when every undirected edge is shared
//! by exactly two triangles. `OPEN_EDGES` counts edges used once (holes) and
//! `NONMANIFOLD_EDGES` counts edges used three or more times (self-touching or
//! duplicated geometry) — two genuinely different defects that a single boolean
//! would conflate. `CONSISTENT_WINDING` is reported separately because a mesh
//! can be perfectly closed yet have a face flipped, which leaves the volume
//! wrong while every edge still pairs up.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Layer};

use crate::args_common::req_str;
use crate::inside_3d::collect_triangles;
use crate::mesh3d::{mesh_volume, topology};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct IsClosed3dTool;

impl Tool for IsClosed3dTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "is_closed_3d",
            display_name: "Is Closed 3D",
            summary: "Flags which 3D multipatch features are watertight closed solids and which are open surfaces, with diagnostic counts of unpaired boundary edges, non-manifold edges and winding consistency (ArcGIS Is Closed 3D). Every volumetric tool in the catalog — union_3d, inside_3d, minimum_bounding_volume, buffer_3d — assumes closed input and returns a meaningless number for an open mesh, and nothing lets you check that precondition first.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "3D features (triangle-mesh MultiPolygons with Z, as buffer_3d, minimum_bounding_volume and voxel_isosurface emit).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path: the input features with the diagnostic fields appended. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "closed_only",
                    description: "Emit only the features that are closed solids (default false: every feature is emitted with its verdict).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        crate::args_common::opt_bool(args, "closed_only")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let closed_only = crate::args_common::bool_or(args, "closed_only", false)?;

        let layer = load_input_layer(input)?;

        let mut out = Layer::new("is_closed_3d");
        out.geom_type = layer.geom_type;
        out.crs = layer.crs.clone();
        // Carry the input schema through so the verdict lands on the caller's
        // own attributes rather than replacing them.
        for f in layer.schema.fields() {
            out.add_field(f.clone());
        }
        out.add_field(FieldDef::new("IS_CLOSED", FieldType::Boolean));
        out.add_field(FieldDef::new("OPEN_EDGES", FieldType::Integer));
        out.add_field(FieldDef::new("NONMANIFOLD_EDGES", FieldType::Integer));
        out.add_field(FieldDef::new("CONSISTENT_WINDING", FieldType::Boolean));
        out.add_field(FieldDef::new("TRI_COUNT", FieldType::Integer));
        out.add_field(FieldDef::new("VOLUME", FieldType::Float));

        let field_names: Vec<String> = layer
            .schema
            .fields()
            .iter()
            .map(|f| f.name.clone())
            .collect();

        let mut closed_count = 0_u64;
        let mut open_count = 0_u64;
        let mut skipped = 0_u64;
        let total = layer.iter().count().max(1);

        for (i, feature) in layer.iter().enumerate() {
            let tris = feature
                .geometry
                .as_ref()
                .map(collect_triangles)
                .unwrap_or_default();
            if tris.is_empty() {
                // Not a mesh at all (a point, a 2D line): there is nothing to
                // judge, so it is neither closed nor a defect.
                skipped += 1;
                if !closed_only {
                    let mut attrs = carry(&field_names, feature);
                    attrs.push(("IS_CLOSED", FieldValue::Boolean(false)));
                    attrs.push(("OPEN_EDGES", FieldValue::Integer(0)));
                    attrs.push(("NONMANIFOLD_EDGES", FieldValue::Integer(0)));
                    attrs.push(("CONSISTENT_WINDING", FieldValue::Boolean(false)));
                    attrs.push(("TRI_COUNT", FieldValue::Integer(0)));
                    attrs.push(("VOLUME", FieldValue::Float(0.0)));
                    out.add_feature(feature.geometry.clone(), &attrs)
                        .map_err(exec)?;
                }
                continue;
            }

            let t = topology(&tris);
            if t.closed {
                closed_count += 1;
            } else {
                open_count += 1;
            }
            if closed_only && !t.closed {
                continue;
            }

            // Only a closed, consistently wound mesh has a meaningful volume;
            // reporting one for an open mesh is exactly the mistake this tool
            // exists to prevent.
            let volume = if t.closed && t.consistent_winding {
                mesh_volume(&tris)
            } else {
                0.0
            };

            let mut attrs = carry(&field_names, feature);
            attrs.push(("IS_CLOSED", FieldValue::Boolean(t.closed)));
            attrs.push(("OPEN_EDGES", FieldValue::Integer(t.open_edges as i64)));
            attrs.push((
                "NONMANIFOLD_EDGES",
                FieldValue::Integer(t.nonmanifold_edges as i64),
            ));
            attrs.push((
                "CONSISTENT_WINDING",
                FieldValue::Boolean(t.consistent_winding),
            ));
            attrs.push(("TRI_COUNT", FieldValue::Integer(tris.len() as i64)));
            attrs.push(("VOLUME", FieldValue::Float(volume)));
            out.add_feature(feature.geometry.clone(), &attrs)
                .map_err(exec)?;

            ctx.progress.progress((i as f64 + 1.0) / total as f64);
        }

        ctx.progress.info(&format!(
            "{closed_count} closed, {open_count} open, {skipped} non-mesh"
        ));
        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("closed_count".to_string(), json!(closed_count));
        outputs.insert("open_count".to_string(), json!(open_count));
        outputs.insert("non_mesh_count".to_string(), json!(skipped));
        Ok(ToolRunResult { outputs })
    }
}

/// Copies the input feature's attributes into name/value pairs for the output.
fn carry<'a>(names: &'a [String], feature: &wbvector::Feature) -> Vec<(&'a str, FieldValue)> {
    names
        .iter()
        .enumerate()
        .filter_map(|(i, n)| feature.attributes.get(i).map(|v| (n.as_str(), v.clone())))
        .collect()
}

fn exec(e: impl std::fmt::Display) -> ToolError {
    ToolError::Execution(e.to_string())
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

    /// A box with its bottom face removed — the canonical open mesh.
    fn open_box() -> Geometry {
        let tris = crate::inside_3d::collect_triangles(&box_mesh([0.0; 3], [1.0, 1.0, 1.0]));
        crate::mesh3d::triangles_to_geometry(&tris[2..])
    }

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = IsClosed3dTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    fn field(layer: &Layer, fid: usize, name: &str) -> FieldValue {
        let i = layer.schema.field_index(name).unwrap();
        layer.iter().nth(fid).unwrap().attributes[i].clone()
    }

    #[test]
    fn a_closed_box_is_reported_closed_with_its_exact_volume() {
        let (out, res) = run(json!({
            "input": layer_of(vec![box_mesh([0.0; 3], [2.0, 3.0, 4.0])]),
        }));
        assert_eq!(field(&out, 0, "IS_CLOSED"), FieldValue::Boolean(true));
        assert_eq!(field(&out, 0, "OPEN_EDGES"), FieldValue::Integer(0));
        assert_eq!(
            field(&out, 0, "CONSISTENT_WINDING"),
            FieldValue::Boolean(true)
        );
        match field(&out, 0, "VOLUME") {
            FieldValue::Float(v) => assert!((v - 24.0).abs() < 1e-9),
            other => panic!("expected a float volume, got {other:?}"),
        }
        assert_eq!(res.outputs["closed_count"], json!(1));
        assert_eq!(res.outputs["open_count"], json!(0));
    }

    #[test]
    fn an_open_box_is_flagged_with_its_boundary_edge_count() {
        let (out, res) = run(json!({"input": layer_of(vec![open_box()])}));
        assert_eq!(field(&out, 0, "IS_CLOSED"), FieldValue::Boolean(false));
        // A square hole has four unpaired edges.
        assert_eq!(field(&out, 0, "OPEN_EDGES"), FieldValue::Integer(4));
        assert_eq!(res.outputs["open_count"], json!(1));
    }

    #[test]
    fn an_open_mesh_reports_no_volume_rather_than_a_misleading_number() {
        // The precise failure this tool exists to prevent: the signed-tet sum
        // of an unclosed mesh is an artefact, so it must not be published.
        let (out, _) = run(json!({"input": layer_of(vec![open_box()])}));
        assert_eq!(field(&out, 0, "VOLUME"), FieldValue::Float(0.0));
    }

    #[test]
    fn closed_only_filters_the_output() {
        let (out, res) = run(json!({
            "input": layer_of(vec![box_mesh([0.0; 3], [1.0, 1.0, 1.0]), open_box()]),
            "closed_only": true,
        }));
        assert_eq!(out.iter().count(), 1);
        // Counts still reflect everything inspected, not just what was emitted.
        assert_eq!(res.outputs["closed_count"], json!(1));
        assert_eq!(res.outputs["open_count"], json!(1));
    }

    #[test]
    fn non_mesh_features_are_counted_separately_not_called_open() {
        let (_, res) = run(json!({
            "input": layer_of(vec![Geometry::Point(wbvector::Coord::xyz(0.0, 0.0, 0.0))]),
        }));
        assert_eq!(res.outputs["non_mesh_count"], json!(1));
        assert_eq!(res.outputs["open_count"], json!(0));
        assert_eq!(res.outputs["closed_count"], json!(0));
    }

    #[test]
    fn input_attributes_are_carried_through() {
        let mut l = Layer::new("in");
        l.geom_type = Some(GeometryType::MultiPolygon);
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.add_feature(
            Some(box_mesh([0.0; 3], [1.0, 1.0, 1.0])),
            &[("name", FieldValue::Text("tower".into()))],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        let path = memory_store::make_vector_memory_path(&id);

        let (out, _) = run(json!({"input": path}));
        assert_eq!(field(&out, 0, "name"), FieldValue::Text("tower".into()));
        assert_eq!(field(&out, 0, "IS_CLOSED"), FieldValue::Boolean(true));
    }

    #[test]
    fn rejects_missing_input() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(IsClosed3dTool.validate(&args).is_err());
    }
}
