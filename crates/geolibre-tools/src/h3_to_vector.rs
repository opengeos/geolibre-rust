//! GeoLibre tool: turn H3 cell IDs into cell-boundary polygons.
//!
//! The inverse of [`crate::vector_to_h3`]: given a layer whose features carry an
//! H3 cell id in a text field (for example a CSV/table of ids, or the output of
//! `vector_to_h3`), emit one polygon per feature tracing that cell's boundary.
//! All of the input's other attributes are copied through unchanged, so a
//! `count`/value column rides along for rendering or joins.
//!
//! H3 work is done by [`h3o`](https://docs.rs/h3o), a pure-Rust reimplementation
//! of H3 with no C dependencies, so it compiles to the same wasm targets as the
//! rest of the suite. See issue #22.

use std::collections::BTreeMap;
use std::str::FromStr;

use h3o::CellIndex;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Default name of the field holding the H3 cell id.
const DEFAULT_FIELD: &str = "h3";

/// Builds H3 cell-boundary polygons from a column of H3 cell ids.
pub struct H3ToVectorTool;

impl Tool for H3ToVectorTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "h3_to_vector",
            display_name: "H3 Cells to Polygons",
            summary: "Turn a column of H3 cell ids into cell-boundary polygons, copying other attributes through.",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector file path (or in-memory handle) with a field of H3 cell ids.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path; the driver is taken from its extension (.geojson, .fgb, .parquet, ...). If omitted, the layer is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "field",
                    description: "Name of the text field holding the H3 cell id; default 'h3'.",
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
        let input = args.get("input").and_then(Value::as_str).ok_or_else(|| {
            ToolError::Validation("missing required parameter 'input'".to_string())
        })?;
        let output = parse_optional_str(args, "output")?;
        let field = parse_optional_str(args, "field")?.unwrap_or(DEFAULT_FIELD);

        ctx.progress.info("reading input vector");
        let layer = load_input_layer(input)?;
        let field_idx = layer.schema.field_index(field).ok_or_else(|| {
            ToolError::Validation(format!("input has no field '{field}' to read H3 ids from"))
        })?;
        let feature_count = layer.len();

        // Copy the input schema so every attribute rides along with the cell.
        let mut out = Layer::new("h3_cells")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        for def in layer.schema.fields() {
            out.add_field(def.clone());
        }
        // Guarantee an `h3` text column exists for the canonical id, even when the
        // source field is named something else.
        let h3_out = if layer.schema.field_index(DEFAULT_FIELD).is_some() {
            DEFAULT_FIELD
        } else {
            out.add_field(FieldDef::new(DEFAULT_FIELD, FieldType::Text));
            DEFAULT_FIELD
        };

        ctx.progress.info("building H3 cell polygons");
        let mut skipped = 0u64;
        for feature in layer.iter() {
            let Some(id) = feature.attributes.get(field_idx).and_then(FieldValue::as_str) else {
                skipped += 1;
                continue;
            };
            let Ok(cell) = CellIndex::from_str(id.trim()) else {
                skipped += 1; // not a valid H3 id
                continue;
            };
            // Carry every source attribute through, then set the canonical id.
            let mut attrs: Vec<(&str, FieldValue)> = layer
                .schema
                .fields()
                .iter()
                .enumerate()
                .map(|(i, def)| (def.name.as_str(), feature.attributes[i].clone()))
                .collect();
            attrs.push((h3_out, FieldValue::Text(cell.to_string())));

            out.add_feature(
                Some(Geometry::polygon(cell_polygon_ring(cell), Vec::new())),
                &attrs,
            )
            .map_err(|e| ToolError::Execution(format!("failed building H3 feature: {e}")))?;
        }

        ctx.progress.info("writing output vector");
        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("skipped".to_string(), json!(skipped));
        Ok(ToolRunResult { outputs })
    }
}

/// Builds the exterior ring (lon/lat coordinates) of an H3 cell's boundary,
/// without the closing duplicate vertex (matching `Ring`).
fn cell_polygon_ring(cell: CellIndex) -> Vec<Coord> {
    cell.boundary()
        .iter()
        .map(|ll| Coord::xy(ll.lng(), ll.lat()))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// A small layer of H3 ids (San Francisco and Paris at res 8) plus a value.
    fn ids_layer() -> Layer {
        let mut l = Layer::new("ids").with_geom_type(GeometryType::Point);
        l.add_field(FieldDef::new("h3", FieldType::Text));
        l.add_field(FieldDef::new("count", FieldType::Integer));
        let sf = h3o::LatLng::new(37.7749, -122.4194)
            .unwrap()
            .to_cell(h3o::Resolution::Eight);
        let paris = h3o::LatLng::new(48.8566, 2.3522)
            .unwrap()
            .to_cell(h3o::Resolution::Eight);
        l.add_feature(None, &[("h3", FieldValue::Text(sf.to_string())), ("count", FieldValue::Integer(2))])
            .unwrap();
        l.add_feature(None, &[("h3", FieldValue::Text(paris.to_string())), ("count", FieldValue::Integer(1))])
            .unwrap();
        l
    }

    #[test]
    fn builds_cell_polygons_and_copies_attributes() {
        let id = wbvector::memory_store::put_vector(ids_layer());
        let input = wbvector::memory_store::make_vector_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({ "input": input })).unwrap();

        let result = H3ToVectorTool.run(&args, &ctx()).unwrap();
        assert_eq!(result.outputs["feature_count"], json!(2));
        assert_eq!(result.outputs["skipped"], json!(0));

        let out = load_input_layer(result.outputs["output"].as_str().unwrap()).unwrap();
        assert_eq!(out.len(), 2);
        assert_eq!(out.geom_type, Some(GeometryType::Polygon));
        // The `count` attribute is carried through unchanged.
        let cidx = out.schema.field_index("count").unwrap();
        let total: i64 = out.iter().map(|f| f.attributes[cidx].as_i64().unwrap()).sum();
        assert_eq!(total, 3);
        // Each cell boundary is a hexagon/pentagon: >= 5 vertices.
        for f in out.iter() {
            match &f.geometry {
                Some(Geometry::Polygon { exterior, .. }) => assert!(exterior.len() >= 5),
                _ => panic!("expected polygon geometry"),
            }
        }
    }

    #[test]
    fn skips_invalid_ids_and_errors_on_missing_field() {
        // Missing field -> validation error.
        let mut l = Layer::new("x").with_geom_type(GeometryType::Point);
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.add_feature(None, &[("name", FieldValue::Text("nope".into()))]).unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let input = wbvector::memory_store::make_vector_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({ "input": input })).unwrap();
        assert!(H3ToVectorTool.run(&args, &ctx()).is_err());

        // Present field but unparsable value -> skipped, not an error.
        let mut l2 = Layer::new("y").with_geom_type(GeometryType::Point);
        l2.add_field(FieldDef::new("h3", FieldType::Text));
        l2.add_feature(None, &[("h3", FieldValue::Text("not-an-h3-id".into()))]).unwrap();
        let id2 = wbvector::memory_store::put_vector(l2);
        let input2 = wbvector::memory_store::make_vector_memory_path(&id2);
        let args2: ToolArgs = serde_json::from_value(json!({ "input": input2 })).unwrap();
        let r = H3ToVectorTool.run(&args2, &ctx()).unwrap();
        assert_eq!(r.outputs["skipped"], json!(1));
        assert_eq!(load_input_layer(r.outputs["output"].as_str().unwrap()).unwrap().len(), 0);
    }
}
