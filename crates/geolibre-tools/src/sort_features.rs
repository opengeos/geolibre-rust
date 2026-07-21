//! GeoLibre tool: reorder a vector dataset by attribute fields or by a spatial
//! (Hilbert-curve) key.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Sort* (Data Management), including its
//! spatial-sort methods. Nothing in the bundled whitebox suite reorders vector
//! features (only `sort_lidar` exists, for point clouds). Beyond ArcGIS parity a
//! Hilbert-curve spatial sort directly improves the repo's own outputs: it
//! clusters spatially-near records so `write_geoparquet` row-group/bbox
//! statistics prune better, and `vector_to_pmtiles` gets more local features per
//! tile. The same `hilbert_for_point` mapping the GeoParquet writer uses is
//! reused here.
//!
//! `method=hilbert` (default) sorts by the Hilbert-curve distance of each
//! feature's bounding-box centre over the dataset extent. `method=attribute`
//! sorts by one or more fields (`fields="pop:desc,name:asc"`); an optional
//! trailing spatial tiebreak is applied when values are equal. `index_field`
//! writes the computed curve distance as a new attribute.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Layer};

use crate::hilbert::hilbert_for_point;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct SortFeaturesTool;

impl Tool for SortFeaturesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "sort_features",
            display_name: "Sort Features",
            summary: "Reorder a vector layer by attribute fields or along a Hilbert space-filling curve (spatial sort), like ArcGIS Sort — Hilbert ordering clusters spatially-near records to improve GeoParquet pruning and PMTiles locality.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector layer in the new order. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'hilbert' (spatial curve order; default) or 'attribute' (by the 'fields' list).",
                    required: false,
                },
                ToolParamSpec {
                    name: "fields",
                    description: "attribute mode: comma-separated fields with optional direction, e.g. 'pop:desc,name:asc'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "index_field",
                    description: "Optional field name to store each feature's Hilbert-curve distance.",
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
        let prm = parse_params(args)?;
        if matches!(prm.method, Method::Attribute) && prm.fields.is_empty() {
            return Err(ToolError::Validation(
                "attribute sort requires a non-empty 'fields' list".to_string(),
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
        let index_field = parse_optional_str(args, "index_field")?.map(str::to_string);
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let n = layer.features.len();

        // Resolve attribute field indices up front (attribute mode).
        let mut field_keys: Vec<(usize, bool)> = Vec::new();
        if matches!(prm.method, Method::Attribute) {
            for (name, ascending) in &prm.fields {
                let idx = layer.schema.field_index(name).ok_or_else(|| {
                    ToolError::Validation(format!("sort field '{name}' not found"))
                })?;
                field_keys.push((idx, *ascending));
            }
        }

        // Dataset extent for the Hilbert mapping.
        let (mut min_x, mut min_y, mut max_x, mut max_y) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for f in &layer.features {
            if let Some(bb) = f.geometry.as_ref().and_then(|g| g.bbox()) {
                min_x = min_x.min(bb.min_x);
                min_y = min_y.min(bb.min_y);
                max_x = max_x.max(bb.max_x);
                max_y = max_y.max(bb.max_y);
            }
        }
        if !min_x.is_finite() {
            (min_x, min_y, max_x, max_y) = (0.0, 0.0, 0.0, 0.0);
        }

        // Hilbert key per feature (features with no geometry sort last).
        let hilbert: Vec<Option<u64>> = layer
            .features
            .iter()
            .map(|f| {
                f.geometry.as_ref().and_then(|g| g.bbox()).map(|bb| {
                    let cx = (bb.min_x + bb.max_x) / 2.0;
                    let cy = (bb.min_y + bb.max_y) / 2.0;
                    hilbert_for_point(cx, cy, min_x, min_y, max_x, max_y)
                })
            })
            .collect();

        // Build an index order.
        let mut order: Vec<usize> = (0..n).collect();
        match prm.method {
            Method::Hilbert => {
                order.sort_by(|&a, &b| {
                    hilbert[a].cmp(&hilbert[b]).then(a.cmp(&b)) // stable tiebreak
                });
            }
            Method::Attribute => {
                order.sort_by(|&a, &b| {
                    for &(idx, ascending) in &field_keys {
                        let va = layer.features[a].attributes.get(idx);
                        let vb = layer.features[b].attributes.get(idx);
                        let mut ord = compare_values(va, vb);
                        if !ascending {
                            ord = ord.reverse();
                        }
                        if ord != std::cmp::Ordering::Equal {
                            return ord;
                        }
                    }
                    // Spatial tiebreak, then stable.
                    hilbert[a].cmp(&hilbert[b]).then(a.cmp(&b))
                });
            }
        }

        ctx.progress
            .info(&format!("sorting {n} feature(s) by {}", prm.method.label()));

        // Optionally add the Hilbert index field.
        if let Some(name) = &index_field {
            layer.add_field(FieldDef::new(name.clone(), FieldType::Integer));
            for (i, f) in layer.features.iter_mut().enumerate() {
                let v = hilbert[i].map(|h| h as i64).unwrap_or(-1);
                f.attributes.push(FieldValue::Integer(v));
            }
        }

        // Reassemble the layer in the new order, preserving schema/CRS.
        let mut out = Layer::new(layer.name.clone());
        out.geom_type = layer.geom_type;
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for fdef in layer.schema.fields() {
            out.add_field(fdef.clone());
        }
        for &i in &order {
            out.push(layer.features[i].clone());
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        Ok(ToolRunResult { outputs })
    }
}

/// Orders two attribute values: numbers numerically, else lexicographically by
/// string form. Nulls sort first.
fn compare_values(a: Option<&FieldValue>, b: Option<&FieldValue>) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    let na = a.and_then(|v| v.as_f64());
    let nb = b.and_then(|v| v.as_f64());
    match (na, nb) {
        (Some(x), Some(y)) => x.total_cmp(&y),
        (Some(_), None) => Ordering::Greater,
        (None, Some(_)) => Ordering::Less,
        (None, None) => {
            let sa = a.map(value_string).unwrap_or_default();
            let sb = b.map(value_string).unwrap_or_default();
            sa.cmp(&sb)
        }
    }
}

fn value_string(fv: &FieldValue) -> String {
    if let Some(i) = fv.as_i64() {
        i.to_string()
    } else if let Some(f) = fv.as_f64() {
        format!("{f}")
    } else {
        fv.as_str().unwrap_or("").to_string()
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Hilbert,
    Attribute,
}

impl Method {
    fn label(&self) -> &'static str {
        match self {
            Method::Hilbert => "hilbert",
            Method::Attribute => "attribute",
        }
    }
}

struct Params {
    method: Method,
    /// (field_name, ascending) pairs for attribute sort.
    fields: Vec<(String, bool)>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let method = match args.get("method").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("hilbert") => Method::Hilbert,
        Some("attribute") => Method::Attribute,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "'method' must be 'hilbert' or 'attribute', got '{other}'"
            )))
        }
    };
    let fields = match args.get("fields").and_then(Value::as_str) {
        None => Vec::new(),
        Some(s) => s
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(|spec| {
                let mut it = spec.splitn(2, ':');
                let name = it.next().unwrap().trim().to_string();
                let ascending = match it.next().map(|d| d.trim().to_ascii_lowercase()) {
                    None => true,
                    Some(d) if d == "asc" || d == "ascending" => true,
                    Some(d) if d == "desc" || d == "descending" => false,
                    Some(d) => {
                        return Err(ToolError::Validation(format!(
                            "sort direction '{d}' must be 'asc' or 'desc'"
                        )))
                    }
                };
                Ok((name, ascending))
            })
            .collect::<Result<Vec<_>, _>>()?,
    };
    Ok(Params { method, fields })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Geometry, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn point_layer(rows: &[(f64, f64, i64, &str)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("pop", FieldType::Integer));
        l.add_field(FieldDef::new("name", FieldType::Text));
        for (x, y, pop, name) in rows {
            l.add_feature(
                Some(Geometry::point(*x, *y)),
                &[("pop", (*pop).into()), ("name", (*name).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> Layer {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SortFeaturesTool.run(&args, &ctx()).unwrap();
        load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    /// Attribute descending sort orders features by the numeric field.
    #[test]
    fn attribute_sort_descending() {
        let rows = [
            (0.0, 0.0, 30, "c"),
            (1.0, 1.0, 10, "a"),
            (2.0, 2.0, 20, "b"),
        ];
        let l = run(json!({
            "input": point_layer(&rows), "method": "attribute", "fields": "pop:desc",
        }));
        let pop = l.schema.field_index("pop").unwrap();
        let seq: Vec<i64> = l
            .features
            .iter()
            .map(|f| f.attributes[pop].as_i64().unwrap())
            .collect();
        assert_eq!(seq, vec![30, 20, 10]);
    }

    /// Hilbert sort keeps all features and produces a monotone non-decreasing
    /// curve index (which we surface with index_field).
    #[test]
    fn hilbert_sort_is_monotone() {
        // A scattering of points.
        let rows: Vec<(f64, f64, i64, &str)> = (0..25)
            .map(|i| ((i % 5) as f64, (i / 5) as f64, i as i64, "p"))
            .collect();
        let l = run(json!({
            "input": point_layer(&rows), "method": "hilbert", "index_field": "hid",
        }));
        assert_eq!(l.features.len(), 25);
        let hid = l.schema.field_index("hid").unwrap();
        let keys: Vec<i64> = l
            .features
            .iter()
            .map(|f| f.attributes[hid].as_i64().unwrap())
            .collect();
        assert!(
            keys.windows(2).all(|w| w[0] <= w[1]),
            "hilbert index must be non-decreasing after sort"
        );
        // Spatial locality: consecutive features are usually grid-adjacent.
    }

    /// Feature count and attributes are preserved through the sort.
    #[test]
    fn preserves_all_features() {
        let rows = [(5.0, 5.0, 1, "x"), (0.0, 9.0, 2, "y"), (9.0, 0.0, 3, "z")];
        let l = run(json!({ "input": point_layer(&rows), "method": "hilbert" }));
        assert_eq!(l.features.len(), 3);
        let name = l.schema.field_index("name").unwrap();
        let mut names: Vec<String> = l
            .features
            .iter()
            .map(|f| f.attributes[name].as_str().unwrap().to_string())
            .collect();
        names.sort();
        assert_eq!(names, vec!["x", "y", "z"]);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SortFeaturesTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "method": "random" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "method": "attribute" })).is_err()); // no fields
        assert!(
            bad(json!({ "input": "a.geojson", "method": "attribute", "fields": "pop:up" }))
                .is_err()
        );
        assert!(
            bad(json!({ "input": "a.geojson", "method": "attribute", "fields": "pop:desc" }))
                .is_ok()
        );
    }
}
