//! GeoLibre tool: reverse the vertex order (start↔end) of polyline features.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Flip Line* (Editing): reverse the
//! digitized direction of line features by reversing their vertex order. This is
//! routine direction normalization for digitizing cleanup and hydrologic /
//! transport network editing, where a handful of arcs were captured pointing the
//! wrong way (e.g. a stream that must flow downhill, or a one-way road digitized
//! against its travel direction).
//!
//! `flip_image` in the bundled suite only reflects rasters; there is no vector
//! line-direction reversal in either catalog. This tool fills that gap:
//!
//! - Each `LineString`'s coordinate vector is reversed in place, so the former
//!   end vertex becomes the start and vice versa. The geometry (the set of
//!   points and the path traced) is identical — only its orientation flips.
//! - For `MultiLineString`, every part is reversed *and* the part order is
//!   reversed, so the whole multipart geometry reads back-to-front.
//! - Non-line features (points, polygons) pass through untouched.
//!
//! Attributes and schema are preserved exactly; only geometry orientation
//! changes. Reported `flipped_count` is how many line features were reversed.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Geometry};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct FlipLineTool;

impl Tool for FlipLineTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "flip_line",
            display_name: "Flip Line",
            summary: "Reverse the vertex order (start↔end) of polyline features to flip their digitized direction.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input line vector file path, format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
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
        let output = parse_optional_str(args, "output")?;

        let mut layer = load_input_layer(input)?;
        let feature_count = layer.features.len();

        let mut flipped = 0usize;
        for feature in &mut layer.features {
            match &mut feature.geometry {
                Some(Geometry::LineString(coords)) => {
                    if flip_coords(coords) {
                        flipped += 1;
                    }
                }
                Some(Geometry::MultiLineString(parts)) => {
                    let mut changed = false;
                    for part in parts.iter_mut() {
                        changed |= flip_coords(part);
                    }
                    // Reverse the part order too, so the multipart geometry as a
                    // whole reads end-to-start.
                    if parts.len() > 1 {
                        parts.reverse();
                        changed = true;
                    }
                    if changed {
                        flipped += 1;
                    }
                }
                _ => {}
            }
        }

        ctx.progress
            .info(&format!("flipped {flipped} of {feature_count} feature(s)"));

        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("flipped_count".to_string(), json!(flipped));
        Ok(ToolRunResult { outputs })
    }
}

/// Reverses a coordinate vector in place. Returns `true` when the order actually
/// changed (a line with 0 or 1 vertices is left untouched).
fn flip_coords(coords: &mut [Coord]) -> bool {
    if coords.len() < 2 {
        return false;
    }
    coords.reverse();
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn line(pts: &[(f64, f64)]) -> Geometry {
        Geometry::LineString(pts.iter().map(|&(x, y)| Coord::xy(x, y)).collect())
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = FlipLineTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn coords_of(g: &Geometry) -> Vec<(f64, f64)> {
        match g {
            Geometry::LineString(cs) => cs.iter().map(|c| (c.x, c.y)).collect(),
            _ => panic!("expected LineString"),
        }
    }

    /// The core property: vertex order is reversed, coordinate set unchanged.
    #[test]
    fn reverses_vertex_order() {
        let mut layer = Layer::new("streams");
        layer
            .add_feature(Some(line(&[(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)])), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input }));
        assert_eq!(out.outputs["flipped_count"], json!(1));
        assert_eq!(
            coords_of(layer.features[0].geometry.as_ref().unwrap()),
            vec![(2.0, 0.0), (1.0, 1.0), (0.0, 0.0)]
        );
    }

    /// Flipping twice is the identity.
    #[test]
    fn double_flip_is_identity() {
        let mut layer = Layer::new("streams");
        layer
            .add_feature(Some(line(&[(0.0, 0.0), (3.0, 4.0), (5.0, 1.0)])), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (_, once) = run_tool(json!({ "input": input }));
        let id2 = memory_store::put_vector(once);
        let (_, twice) = run_tool(json!({
            "input": memory_store::make_vector_memory_path(&id2)
        }));
        assert_eq!(
            coords_of(twice.features[0].geometry.as_ref().unwrap()),
            vec![(0.0, 0.0), (3.0, 4.0), (5.0, 1.0)]
        );
    }

    /// Attributes are preserved on the flipped feature.
    #[test]
    fn preserves_attributes() {
        let mut layer = Layer::new("streams");
        layer.add_field(FieldDef::new("name", FieldType::Text));
        layer
            .add_feature(
                Some(line(&[(0.0, 0.0), (1.0, 0.0)])),
                &[("name", "creek".into())],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (_, layer) = run_tool(json!({ "input": input }));
        assert_eq!(
            layer.features[0]
                .get(&layer.schema, "name")
                .unwrap()
                .as_str(),
            Some("creek")
        );
    }

    /// MultiLineString: each part reversed and the part order reversed.
    #[test]
    fn flips_multilinestring() {
        let mut layer = Layer::new("streams");
        let g = Geometry::MultiLineString(vec![
            vec![Coord::xy(0.0, 0.0), Coord::xy(1.0, 0.0)],
            vec![Coord::xy(2.0, 0.0), Coord::xy(3.0, 0.0)],
        ]);
        layer.add_feature(Some(g), &[]).unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input }));
        assert_eq!(out.outputs["flipped_count"], json!(1));
        match layer.features[0].geometry.as_ref().unwrap() {
            Geometry::MultiLineString(parts) => {
                assert_eq!(parts[0][0], Coord::xy(3.0, 0.0));
                assert_eq!(parts[0][1], Coord::xy(2.0, 0.0));
                assert_eq!(parts[1][0], Coord::xy(1.0, 0.0));
                assert_eq!(parts[1][1], Coord::xy(0.0, 0.0));
            }
            _ => panic!("expected MultiLineString"),
        }
    }

    /// Non-line features pass through untouched and are not counted as flipped.
    #[test]
    fn passes_points_through() {
        let mut layer = Layer::new("mixed");
        layer
            .add_feature(Some(Geometry::point(5.0, 5.0)), &[])
            .unwrap();
        layer
            .add_feature(Some(line(&[(0.0, 0.0), (1.0, 0.0)])), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input }));
        assert_eq!(out.outputs["feature_count"], json!(2));
        assert_eq!(out.outputs["flipped_count"], json!(1));
        assert!(matches!(
            layer.features[0].geometry,
            Some(Geometry::Point(_))
        ));
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = FlipLineTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "" })).is_err());
        assert!(bad(json!({ "input": "x.geojson" })).is_ok());
    }
}
