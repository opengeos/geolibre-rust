//! GeoLibre tool: reshape a wide attribute table to long form.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Transpose Fields* (Data Management),
//! and the exact inverse of GeoLibre's `pivot_table`.
//!
//! ## Why both directions are needed
//!
//! `pivot_table` goes **long to wide**, and its module doc explains why:
//! `multivariate_clustering`, `dimension_reduction`,
//! `spatially_constrained_multivariate_clustering`, `similarity_search` and
//! `calculate_composite_index` all want one row per feature and one column per
//! variable.
//!
//! The *time-series* half of the catalog wants the opposite — one row per
//! (entity, time) observation: `time_series_clustering`, `time_series_forecast`,
//! `time_series_smoothing`, `time_series_cross_correlation`,
//! `change_point_detection`, `summarize_percent_change`,
//! `emerging_hot_spot_analysis`, `estimate_time_to_event`.
//!
//! Real tables arrive wide (`POP_2000`, `POP_2010`, `POP_2020` as separate
//! columns) and there was no way to feed them to any of those. `summary_statistics`
//! aggregates groups but does not reshape.
//!
//! ## Output value typing
//!
//! The transposed columns are scanned **before** any row is written. If they are
//! uniformly integer the output value column is integer, uniformly numeric it is
//! float, otherwise text. Deciding per row instead would produce a column whose
//! type depends on which row happened to be written first, and coercing a text
//! column to float silently would drop data.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Feature, FieldDef, FieldType, FieldValue, Layer};

use crate::args_common::bool_or;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Ceiling on emitted rows. Every row carries a geometry clone, so the output
/// is resident in full before the writer runs.
const MAX_OUTPUT_ROWS: u64 = 20_000_000;

pub struct TransposeFieldsTool;

impl Tool for TransposeFieldsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "transpose_fields",
            display_name: "Transpose Fields",
            summary: "Reshapes a wide attribute table to long form: each listed value column becomes its own row, tagged with the column it came from (ArcGIS Transpose Fields). pivot_table only goes long-to-wide, so tables arriving with one column per period could not be fed to the time-series tools that need one row per observation.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer or table.",
                    required: true,
                },
                ToolParamSpec {
                    name: "transpose_fields",
                    description: "Comma- or semicolon-separated value columns to unpivot.",
                    required: true,
                },
                ToolParamSpec {
                    name: "value_field",
                    description: "Name of the output column holding the cell value (default 'VALUE').",
                    required: false,
                },
                ToolParamSpec {
                    name: "transposed_field",
                    description: "Name of the output column holding the source column name (default 'FIELD').",
                    required: false,
                },
                ToolParamSpec {
                    name: "field_labels",
                    description: "Optional labels written into 'transposed_field' instead of the raw column names, e.g. '2000,2010,2020'. Must have the same length as 'transpose_fields'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "retain_fields",
                    description: "Attributes copied onto every output row. Defaults to all non-transposed fields.",
                    required: false,
                },
                ToolParamSpec {
                    name: "drop_nulls",
                    description: "Skip output rows whose value is null (default false).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output layer. If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if parse_optional_str(args, "input")?.is_none() {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        let fields = parse_list(args, "transpose_fields")?;
        if fields.is_empty() {
            return Err(ToolError::Validation(
                "'transpose_fields' must name at least one column".to_string(),
            ));
        }
        let labels = parse_list(args, "field_labels")?;
        if !labels.is_empty() && labels.len() != fields.len() {
            return Err(ToolError::Validation(format!(
                "'field_labels' has {} entries but 'transpose_fields' has {}; they must match \
                 one-to-one",
                labels.len(),
                fields.len()
            )));
        }
        bool_or(args, "drop_nulls", false)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = parse_optional_str(args, "input")?
            .ok_or_else(|| ToolError::Validation("missing required parameter 'input'".into()))?;
        let fields = parse_list(args, "transpose_fields")?;
        let labels = parse_list(args, "field_labels")?;
        let value_name = parse_optional_str(args, "value_field")?
            .unwrap_or("VALUE")
            .to_string();
        let field_name = parse_optional_str(args, "transposed_field")?
            .unwrap_or("FIELD")
            .to_string();
        let retain_spec = parse_list(args, "retain_fields")?;
        let drop_nulls = bool_or(args, "drop_nulls", false)?;
        let output = parse_optional_str(args, "output")?;

        let layer = load_input_layer(input)?;

        // features x columns rows, each cloning the geometry. A 100-column
        // table over 100k polygons is 10M features held in memory before the
        // writer runs, so fail with a clear message rather than exhausting it.
        let projected = (layer.features.len() as u64).saturating_mul(fields.len() as u64);
        if projected > MAX_OUTPUT_ROWS {
            return Err(ToolError::Validation(format!(
                "this would emit {projected} rows ({} feature(s) x {} column(s)), over the \
                 {MAX_OUTPUT_ROWS}-row limit; transpose fewer columns or split the input",
                layer.features.len(),
                fields.len()
            )));
        }

        // Resolve the transposed columns up front so a typo fails before any
        // work rather than producing a silently empty column.
        let mut t_idx = Vec::with_capacity(fields.len());
        for f in &fields {
            let i = layer.schema.field_index(f).ok_or_else(|| {
                ToolError::Validation(format!(
                    "transpose field '{f}' not found in the input layer"
                ))
            })?;
            t_idx.push(i);
        }

        // Retained columns: an explicit list, or everything not transposed.
        let retain: Vec<(String, usize, FieldType)> = if retain_spec.is_empty() {
            layer
                .schema
                .fields()
                .iter()
                .enumerate()
                .filter(|(i, _)| !t_idx.contains(i))
                .map(|(i, f)| (f.name.clone(), i, f.field_type))
                .collect()
        } else {
            let mut v = Vec::with_capacity(retain_spec.len());
            for name in &retain_spec {
                let i = layer.schema.field_index(name).ok_or_else(|| {
                    ToolError::Validation(format!(
                        "retain field '{name}' not found in the input layer"
                    ))
                })?;
                v.push((name.clone(), i, layer.schema.fields()[i].field_type));
            }
            v
        };

        // A retained column colliding with an appended one would silently
        // misalign the output: wbvector keeps the first schema entry for a
        // duplicated name, so the second FieldDef is dropped while attributes
        // are still written positionally.
        let src_fid = "SRC_FID".to_string();
        let appended = [&field_name, &value_name, &src_fid];
        for name in appended {
            if retain.iter().any(|(n, _, _)| n == name) {
                return Err(ToolError::Validation(format!(
                    "the retained fields already include '{name}', which this tool appends; \
                     rename it or narrow 'retain_fields'"
                )));
            }
        }
        // The appended names must also be distinct from EACH OTHER. Setting
        // value_field == transposed_field, or either to SRC_FID, calls
        // add_field twice with one name; wbvector keeps the first entry, so
        // every attribute after the duplicate is read from the wrong column.
        for (i, a) in appended.iter().enumerate() {
            for b in appended.iter().skip(i + 1) {
                if a == b {
                    return Err(ToolError::Validation(format!(
                        "'{a}' is used for more than one output column; 'value_field', \
                         'transposed_field' and the SRC_FID column must all differ"
                    )));
                }
            }
        }

        // Decide the value column's type from ALL transposed columns at once.
        let value_type = unified_value_type(&layer, &t_idx);
        ctx.progress.info(&format!(
            "{} feature(s) x {} column(s) -> long form, value type {value_type:?}",
            layer.features.len(),
            fields.len()
        ));

        let mut out = Layer::new(format!("{}_long", layer.name));
        if let Some(gt) = layer.geom_type {
            out = out.with_geom_type(gt);
        }
        if let Some(e) = layer.crs_epsg() {
            out = out.with_crs_epsg(e);
        }
        for (name, _, ft) in &retain {
            out.add_field(FieldDef::new(name.as_str(), *ft));
        }
        out.add_field(FieldDef::new(field_name.as_str(), FieldType::Text));
        out.add_field(FieldDef::new(value_name.as_str(), value_type));
        out.add_field(FieldDef::new("SRC_FID", FieldType::Integer));

        let mut dropped = 0_u64;
        for (src_fid, feature) in layer.iter().enumerate() {
            for (k, &fi) in t_idx.iter().enumerate() {
                let raw = feature
                    .attributes
                    .get(fi)
                    .cloned()
                    .unwrap_or(FieldValue::Null);
                if drop_nulls && matches!(raw, FieldValue::Null) {
                    dropped += 1;
                    continue;
                }
                let label = labels.get(k).cloned().unwrap_or_else(|| fields[k].clone());
                let mut attrs: Vec<FieldValue> = retain
                    .iter()
                    .map(|(_, i, _)| {
                        feature
                            .attributes
                            .get(*i)
                            .cloned()
                            .unwrap_or(FieldValue::Null)
                    })
                    .collect();
                attrs.push(FieldValue::Text(label));
                attrs.push(coerce(raw, value_type));
                attrs.push(FieldValue::Integer(src_fid as i64));
                out.push(Feature {
                    fid: 0,
                    geometry: feature.geometry.clone(),
                    attributes: attrs,
                });
            }
        }

        let row_count = out.features.len();
        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("row_count".to_string(), json!(row_count));
        outputs.insert("source_features".to_string(), json!(layer.features.len()));
        outputs.insert("transposed_fields".to_string(), json!(fields.len()));
        outputs.insert("dropped_nulls".to_string(), json!(dropped));
        outputs.insert("value_type".to_string(), json!(format!("{value_type:?}")));
        Ok(ToolRunResult { outputs })
    }
}

/// Picks one type that can hold every transposed column's values.
///
/// Integer only when every source column is integer; float when they are all
/// numeric; text otherwise. Scanning declared schema types is not enough on its
/// own for text columns that happen to hold numbers, so the actual values are
/// checked too — a text column of numerals is still usable as a float column,
/// and forcing it to text would break the numeric time-series consumers.
fn unified_value_type(layer: &Layer, t_idx: &[usize]) -> FieldType {
    let mut all_integer = true;
    let mut all_numeric = true;
    for &i in t_idx {
        match layer.schema.fields()[i].field_type {
            FieldType::Integer => {}
            FieldType::Float => all_integer = false,
            _ => {
                // Declared non-numeric: fall back to inspecting the values.
                all_integer = false;
                for f in layer.iter() {
                    match f.attributes.get(i) {
                        None | Some(FieldValue::Null) => {}
                        Some(v) => {
                            if as_f64(v).is_none() {
                                all_numeric = false;
                                break;
                            }
                        }
                    }
                }
            }
        }
        if !all_numeric {
            break;
        }
    }
    if all_integer {
        FieldType::Integer
    } else if all_numeric {
        FieldType::Float
    } else {
        FieldType::Text
    }
}

fn coerce(v: FieldValue, target: FieldType) -> FieldValue {
    match (target, &v) {
        (_, FieldValue::Null) => FieldValue::Null,
        (FieldType::Integer, _) => {
            as_f64(&v).map_or(FieldValue::Null, |f| FieldValue::Integer(f.round() as i64))
        }
        (FieldType::Float, _) => as_f64(&v).map_or(FieldValue::Null, FieldValue::Float),
        (FieldType::Text, FieldValue::Text(s)) => FieldValue::Text(s.clone()),
        (FieldType::Text, _) => FieldValue::Text(display(&v)),
        _ => v,
    }
}

fn as_f64(v: &FieldValue) -> Option<f64> {
    match v {
        FieldValue::Integer(i) => Some(*i as f64),
        FieldValue::Float(f) => Some(*f),
        FieldValue::Boolean(b) => Some(if *b { 1.0 } else { 0.0 }),
        FieldValue::Text(s) => s.trim().parse::<f64>().ok(),
        _ => None,
    }
}

fn display(v: &FieldValue) -> String {
    match v {
        FieldValue::Null => String::new(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) => format!("{f}"),
        FieldValue::Text(s) => s.clone(),
        FieldValue::Boolean(b) => b.to_string(),
        FieldValue::Date(s) | FieldValue::DateTime(s) => s.clone(),
        FieldValue::Blob(b) => format!("blob[{}]", b.len()),
    }
}

/// Reads a list parameter given either as a delimited string or a JSON array.
fn parse_list(args: &ToolArgs, key: &str) -> Result<Vec<String>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(Vec::new()),
        Some(Value::String(s)) => Ok(s
            .split([',', ';'])
            .map(str::trim)
            .filter(|p| !p.is_empty())
            .map(str::to_string)
            .collect()),
        Some(Value::Array(items)) => items
            .iter()
            .map(|v| {
                v.as_str().map(str::to_string).ok_or_else(|| {
                    ToolError::Validation(format!("every entry of '{key}' must be a string"))
                })
            })
            .collect(),
        Some(_) => Err(ToolError::Validation(format!(
            "'{key}' must be a delimited string or an array of strings"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{Geometry, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Two sites, population in three year columns.
    fn wide() -> String {
        let mut l = Layer::new("sites")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("site", FieldType::Text));
        l.add_field(FieldDef::new("POP_2000", FieldType::Integer));
        l.add_field(FieldDef::new("POP_2010", FieldType::Integer));
        l.add_field(FieldDef::new("POP_2020", FieldType::Integer));
        for (i, (name, a, b, c)) in [("north", 10, 20, 30), ("south", 40, 50, 60)]
            .into_iter()
            .enumerate()
        {
            l.add_feature(
                Some(Geometry::point(i as f64, 0.0)),
                &[
                    ("site", name.into()),
                    ("POP_2000", (a as i64).into()),
                    ("POP_2010", (b as i64).into()),
                    ("POP_2020", (c as i64).into()),
                ],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = TransposeFieldsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    fn col(layer: &Layer, name: &str) -> Vec<FieldValue> {
        let i = layer.schema.field_index(name).unwrap();
        layer.iter().map(|f| f.attributes[i].clone()).collect()
    }

    #[test]
    fn produces_one_row_per_feature_and_column() {
        let (layer, res) = run(json!({
            "input": wide(), "transpose_fields": "POP_2000,POP_2010,POP_2020",
        }));
        assert_eq!(layer.features.len(), 6);
        assert_eq!(res.outputs["row_count"], json!(6));
        assert_eq!(res.outputs["source_features"], json!(2));
    }

    #[test]
    fn the_field_column_records_the_source_column_name() {
        let (layer, _) = run(json!({
            "input": wide(), "transpose_fields": "POP_2000,POP_2010,POP_2020",
        }));
        assert_eq!(
            col(&layer, "FIELD")[..3],
            [
                FieldValue::Text("POP_2000".into()),
                FieldValue::Text("POP_2010".into()),
                FieldValue::Text("POP_2020".into()),
            ]
        );
        assert_eq!(
            col(&layer, "VALUE")[..3],
            [
                FieldValue::Integer(10),
                FieldValue::Integer(20),
                FieldValue::Integer(30),
            ]
        );
    }

    #[test]
    fn field_labels_replace_the_raw_column_names() {
        let (layer, _) = run(json!({
            "input": wide(), "transpose_fields": "POP_2000,POP_2010,POP_2020",
            "field_labels": "2000,2010,2020",
        }));
        assert_eq!(col(&layer, "FIELD")[0], FieldValue::Text("2000".into()));
    }

    #[test]
    fn retained_attributes_and_geometry_are_copied_to_every_row() {
        let (layer, _) = run(json!({
            "input": wide(), "transpose_fields": "POP_2000,POP_2010,POP_2020",
        }));
        let sites = col(&layer, "site");
        assert_eq!(sites[0], FieldValue::Text("north".into()));
        assert_eq!(sites[2], FieldValue::Text("north".into()));
        assert_eq!(sites[3], FieldValue::Text("south".into()));
        assert!(layer.iter().all(|f| f.geometry.is_some()));
    }

    #[test]
    fn src_fid_points_back_at_the_source_feature() {
        let (layer, _) = run(json!({
            "input": wide(), "transpose_fields": "POP_2000,POP_2020",
        }));
        assert_eq!(
            col(&layer, "SRC_FID"),
            [
                FieldValue::Integer(0),
                FieldValue::Integer(0),
                FieldValue::Integer(1),
                FieldValue::Integer(1),
            ]
        );
    }

    #[test]
    fn custom_output_column_names_are_honoured() {
        let (layer, _) = run(json!({
            "input": wide(), "transpose_fields": "POP_2000,POP_2010",
            "value_field": "pop", "transposed_field": "year",
        }));
        assert!(layer.schema.field_index("pop").is_some());
        assert!(layer.schema.field_index("year").is_some());
        assert!(layer.schema.field_index("VALUE").is_none());
    }

    #[test]
    fn retain_fields_narrows_the_carried_attributes() {
        let (layer, _) = run(json!({
            "input": wide(), "transpose_fields": "POP_2000,POP_2010",
            "retain_fields": "site",
        }));
        assert!(layer.schema.field_index("site").is_some());
        // POP_2020 was neither transposed nor retained, so it must not appear.
        assert!(layer.schema.field_index("POP_2020").is_none());
    }

    #[test]
    fn a_mixed_numeric_set_widens_to_float_rather_than_truncating() {
        // One integer column and one float column: the shared value column has
        // to be float, or 1.5 would be written as 2 (or 1).
        let mut l = Layer::new("m");
        l.add_field(FieldDef::new("a", FieldType::Integer));
        l.add_field(FieldDef::new("b", FieldType::Float));
        l.add_feature(None, &[("a", 3i64.into()), ("b", 1.5f64.into())])
            .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let path = wbvector::memory_store::make_vector_memory_path(&id);
        let (layer, res) = run(json!({"input": path, "transpose_fields": "a,b"}));
        assert_eq!(res.outputs["value_type"], json!("Float"));
        assert_eq!(
            col(&layer, "VALUE"),
            [FieldValue::Float(3.0), FieldValue::Float(1.5)]
        );
    }

    #[test]
    fn a_genuinely_textual_column_forces_a_text_value_column() {
        // Coercing this to float would silently null out "high"/"low".
        let mut l = Layer::new("m");
        l.add_field(FieldDef::new("a", FieldType::Integer));
        l.add_field(FieldDef::new("b", FieldType::Text));
        l.add_feature(None, &[("a", 3i64.into()), ("b", "high".into())])
            .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let path = wbvector::memory_store::make_vector_memory_path(&id);
        let (layer, res) = run(json!({"input": path, "transpose_fields": "a,b"}));
        assert_eq!(res.outputs["value_type"], json!("Text"));
        assert_eq!(
            col(&layer, "VALUE"),
            [
                FieldValue::Text("3".into()),
                FieldValue::Text("high".into())
            ]
        );
    }

    #[test]
    fn a_text_column_holding_numerals_still_reads_as_numeric() {
        let mut l = Layer::new("m");
        l.add_field(FieldDef::new("a", FieldType::Text));
        l.add_feature(None, &[("a", "2.5".into())]).unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let path = wbvector::memory_store::make_vector_memory_path(&id);
        let (layer, res) = run(json!({"input": path, "transpose_fields": "a"}));
        assert_eq!(res.outputs["value_type"], json!("Float"));
        assert_eq!(col(&layer, "VALUE"), [FieldValue::Float(2.5)]);
    }

    #[test]
    fn drop_nulls_removes_empty_observations_and_counts_them() {
        let mut l = Layer::new("m");
        l.add_field(FieldDef::new("a", FieldType::Integer));
        l.add_field(FieldDef::new("b", FieldType::Integer));
        l.push(Feature {
            fid: 0,
            geometry: None,
            attributes: vec![FieldValue::Integer(1), FieldValue::Null],
        });
        let id = wbvector::memory_store::put_vector(l);
        let path = wbvector::memory_store::make_vector_memory_path(&id);
        let (layer, res) = run(json!({
            "input": path, "transpose_fields": "a,b", "drop_nulls": true,
        }));
        assert_eq!(layer.features.len(), 1);
        assert_eq!(res.outputs["dropped_nulls"], json!(1));
    }

    #[test]
    fn round_trips_back_through_pivot_table() {
        // The two tools are inverses; if the long form this emits cannot be
        // pivoted back to the original wide values, one of them is wrong.
        let (long, _) = run(json!({
            "input": wide(), "transpose_fields": "POP_2000,POP_2010,POP_2020",
        }));
        let id = wbvector::memory_store::put_vector(long);
        let long_path = wbvector::memory_store::make_vector_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": long_path,
            "fields": "site",
            "pivot_field": "FIELD",
            "value_field": "VALUE",
        }))
        .unwrap();
        let res = crate::pivot_table::PivotTableTool
            .run(&args, &ctx())
            .unwrap();
        let wide_again = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        assert_eq!(wide_again.features.len(), 2);
        let i = wide_again.schema.field_index("POP_2010").unwrap();
        let vals: Vec<f64> = wide_again
            .iter()
            .map(|f| as_f64(&f.attributes[i]).unwrap())
            .collect();
        assert_eq!(vals, vec![20.0, 50.0]);
    }

    #[test]
    fn a_retained_field_colliding_with_an_appended_one_is_rejected() {
        // wbvector keeps the first schema entry for a duplicate name, so this
        // would silently misalign every output attribute.
        let mut l = Layer::new("m");
        l.add_field(FieldDef::new("VALUE", FieldType::Text));
        l.add_field(FieldDef::new("a", FieldType::Integer));
        l.add_feature(None, &[("VALUE", "x".into()), ("a", 1i64.into())])
            .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let path = wbvector::memory_store::make_vector_memory_path(&id);
        let args: ToolArgs =
            serde_json::from_value(json!({"input": path, "transpose_fields": "a"})).unwrap();
        let err = TransposeFieldsTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("VALUE"), "{err}");
    }

    #[test]
    fn appended_columns_colliding_with_each_other_are_rejected() {
        // add_field would be called twice with one name; wbvector keeps the
        // first entry, so every later attribute reads from the wrong column.
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            TransposeFieldsTool.run(&args, &ctx()).is_err()
        };
        assert!(bad(json!({
            "input": wide(), "transpose_fields": "POP_2000,POP_2010",
            "value_field": "V", "transposed_field": "V",
        })));
        assert!(bad(json!({
            "input": wide(), "transpose_fields": "POP_2000,POP_2010",
            "value_field": "SRC_FID",
        })));
        assert!(bad(json!({
            "input": wide(), "transpose_fields": "POP_2000,POP_2010",
            "transposed_field": "SRC_FID",
        })));
    }

    #[test]
    fn an_unknown_column_name_fails_loudly() {
        let args: ToolArgs =
            serde_json::from_value(json!({"input": wide(), "transpose_fields": "NOPE"})).unwrap();
        let err = TransposeFieldsTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("NOPE"), "{err}");
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            TransposeFieldsTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"input": "x.shp"})));
        assert!(bad(json!({"input": "x.shp", "transpose_fields": ""})));
        // Mismatched label count would zip short and mislabel every row.
        assert!(bad(json!({
            "input": "x.shp", "transpose_fields": "a,b,c", "field_labels": "1,2",
        })));
    }
}
