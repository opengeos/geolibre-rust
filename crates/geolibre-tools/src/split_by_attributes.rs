//! GeoLibre tool: split a layer into one output file per attribute value.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Split By Attributes* (Analysis). No
//! bundled tool partitions a dataset by attribute into separate outputs — a
//! constant preprocessing chore (per-county files, per-class layers, per-year
//! exports) that otherwise needs external scripting.
//!
//! Features are grouped by the combined value of one or more `fields`; each
//! group is written to `output_dir/<sanitized_value>.<format>`, preserving the
//! input schema, geometry type, and CRS. A `split_summary.csv` lists every group
//! with its value, feature count, and filename. Filenames are sanitized to
//! alphanumerics and underscores; multi-field keys join with `__`.

use std::collections::BTreeMap;
use std::path::Path;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldValue, Layer, VectorFormat};

use crate::common::write_text_output;
use crate::vector_common::load_input_layer;

pub struct SplitByAttributesTool;

impl Tool for SplitByAttributesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "split_by_attributes",
            display_name: "Split By Attributes",
            summary: "Split a vector layer into one output file per unique value (or combination) of the given field(s), like ArcGIS Split By Attributes.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output_dir",
                    description: "Output directory; one file per group is written here (created if needed).",
                    required: true,
                },
                ToolParamSpec {
                    name: "fields",
                    description: "Comma-separated field name(s) to split on (a multi-field key joins values with '__').",
                    required: true,
                },
                ToolParamSpec {
                    name: "format",
                    description: "Output format extension: geojson (default), fgb, parquet, shp.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "output_dir")?;
        require_str(args, "fields")?;
        parse_format(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output_dir = require_str(args, "output_dir")?;
        let fields_arg = require_str(args, "fields")?;
        let (ext, format) = parse_format(args)?;

        let layer = load_input_layer(input)?;

        // Resolve the split field indices.
        let field_names: Vec<&str> = fields_arg
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .collect();
        if field_names.is_empty() {
            return Err(ToolError::Validation(
                "'fields' must name at least one field".to_string(),
            ));
        }
        let mut field_idx = Vec::new();
        for name in &field_names {
            let i = layer
                .schema
                .field_index(name)
                .ok_or_else(|| ToolError::Validation(format!("field '{name}' not found")))?;
            field_idx.push(i);
        }

        // Group feature indices by combined key.
        let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for (fidx, feature) in layer.features.iter().enumerate() {
            let key = field_idx
                .iter()
                .map(|&i| value_string(&feature.attributes[i]))
                .collect::<Vec<_>>()
                .join("__");
            groups.entry(key).or_default().push(fidx);
        }

        ctx.progress.info(&format!(
            "{} feature(s) -> {} group(s) by {}",
            layer.len(),
            groups.len(),
            field_names.join(", ")
        ));

        std::fs::create_dir_all(output_dir)
            .map_err(|e| ToolError::Execution(format!("failed creating output directory: {e}")))?;

        // Write one file per group, preserving schema/geometry/CRS.
        let mut summary = String::from("value,feature_count,filename\n");
        let mut used_names: std::collections::HashSet<String> = std::collections::HashSet::new();
        let mut files: Vec<String> = Vec::new();
        for (key, indices) in &groups {
            let mut sanitized = sanitize(key);
            // Ensure uniqueness after sanitization.
            let base = sanitized.clone();
            let mut n = 1;
            while !used_names.insert(sanitized.clone()) {
                n += 1;
                sanitized = format!("{base}_{n}");
            }
            let filename = format!("{sanitized}.{ext}");
            let path = Path::new(output_dir).join(&filename);
            let path_str = path.to_string_lossy().to_string();

            let mut group_layer = Layer::new(&sanitized);
            group_layer.schema = layer.schema.clone();
            group_layer.geom_type = layer.geom_type;
            if let Some(epsg) = layer.crs_epsg() {
                group_layer = group_layer.with_crs_epsg(epsg);
            }
            for &fi in indices {
                group_layer.push(layer.features[fi].clone());
            }
            wbvector::write(&group_layer, &path_str, format)
                .map_err(|e| ToolError::Execution(format!("failed writing group '{key}': {e}")))?;

            summary.push_str(&format!("{key},{},{filename}\n", indices.len()));
            files.push(path_str);
        }

        let summary_path = Path::new(output_dir).join("split_summary.csv");
        write_text_output(&summary, &summary_path.to_string_lossy())?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output_dir".to_string(), json!(output_dir));
        outputs.insert("group_count".to_string(), json!(groups.len()));
        outputs.insert("feature_count".to_string(), json!(layer.len()));
        outputs.insert(
            "summary".to_string(),
            json!(summary_path.to_string_lossy().to_string()),
        );
        outputs.insert("files".to_string(), json!(files));
        Ok(ToolRunResult { outputs })
    }
}

/// A stable string for a field value used as the group key.
fn value_string(fv: &FieldValue) -> String {
    match fv {
        FieldValue::Null => "NULL".to_string(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) => format!("{f}"),
        FieldValue::Text(s) => s.clone(),
        FieldValue::Boolean(b) => b.to_string(),
        other => other.as_f64().map(|v| v.to_string()).unwrap_or_default(),
    }
}

/// Sanitizes a key to a safe filename stem (alphanumerics and underscores).
fn sanitize(key: &str) -> String {
    let mut out = String::with_capacity(key.len());
    let mut last_us = false;
    for ch in key.chars() {
        if ch.is_ascii_alphanumeric() || ch == '-' {
            out.push(ch);
            last_us = false;
        } else if !last_us {
            out.push('_');
            last_us = true;
        }
    }
    let trimmed = out.trim_matches('_').to_string();
    if trimmed.is_empty() {
        "value".to_string()
    } else {
        trimmed
    }
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

/// Parses the `format` parameter to a `(extension, VectorFormat)`.
fn parse_format(args: &ToolArgs) -> Result<(&'static str, VectorFormat), ToolError> {
    let f = match args.get("format") {
        None | Some(Value::Null) => "geojson",
        Some(Value::String(s)) if s.trim().is_empty() => "geojson",
        Some(Value::String(s)) => s.trim(),
        Some(_) => {
            return Err(ToolError::Validation(
                "'format' must be a string".to_string(),
            ))
        }
    };
    match f.to_ascii_lowercase().as_str() {
        "geojson" | "json" => Ok(("geojson", VectorFormat::GeoJson)),
        "fgb" | "flatgeobuf" => Ok(("fgb", VectorFormat::FlatGeobuf)),
        "parquet" | "geoparquet" => Ok(("parquet", VectorFormat::GeoParquet)),
        "shp" | "shapefile" => Ok(("shp", VectorFormat::Shapefile)),
        other => Err(ToolError::Validation(format!(
            "unsupported 'format' '{other}' (use geojson, fgb, parquet, shp)"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, Geometry, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn input_layer() -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("cls", FieldType::Text));
        l.add_field(FieldDef::new("yr", FieldType::Integer));
        let rows = [
            ("water", 2020),
            ("water", 2021),
            ("urban", 2020),
            ("forest", 2020),
            ("urban", 2020),
        ];
        for (i, (c, y)) in rows.iter().enumerate() {
            l.add_feature(
                Some(Geometry::point(i as f64, 0.0)),
                &[("cls", (*c).into()), ("yr", (*y as i64).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn tmp_dir(name: &str) -> String {
        let d = std::env::temp_dir().join(format!("geolibre_split_{name}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        d.to_string_lossy().to_string()
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        SplitByAttributesTool.run(&args, &ctx()).unwrap()
    }

    #[test]
    fn splits_by_single_field() {
        let dir = tmp_dir("single");
        let out = run(json!({
            "input": input_layer(), "output_dir": dir, "fields": "cls",
        }));
        // 3 classes: water, urban, forest.
        assert_eq!(out.outputs["group_count"], json!(3));
        // Each group file loads back with the right feature count.
        let water = load_input_layer(&format!("{dir}/water.geojson")).unwrap();
        assert_eq!(water.len(), 2);
        let forest = load_input_layer(&format!("{dir}/forest.geojson")).unwrap();
        assert_eq!(forest.len(), 1);
        // Summary exists.
        assert!(std::path::Path::new(&format!("{dir}/split_summary.csv")).exists());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn splits_by_multiple_fields() {
        let dir = tmp_dir("multi");
        let out = run(json!({
            "input": input_layer(), "output_dir": dir, "fields": "cls,yr",
        }));
        // (water,2020),(water,2021),(urban,2020),(forest,2020) = 4 groups.
        // The "__" key separator sanitizes to a single "_" in the filename.
        assert_eq!(out.outputs["group_count"], json!(4));
        let urban = load_input_layer(&format!("{dir}/urban_2020.geojson")).unwrap();
        assert_eq!(urban.len(), 2);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn sanitize_makes_safe_names() {
        assert_eq!(sanitize("New York / NY"), "New_York_NY");
        assert_eq!(sanitize("a..b"), "a_b");
        assert_eq!(sanitize("///"), "value");
    }

    #[test]
    fn rejects_missing_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SplitByAttributesTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "output_dir": "/tmp/x" })).is_err());
        assert!(bad(
            json!({ "input": "a.geojson", "output_dir": "/tmp/x", "fields": "c", "format": "kml" })
        )
        .is_err());
        assert!(
            bad(json!({ "input": "a.geojson", "output_dir": "/tmp/x", "fields": "c" })).is_ok()
        );
    }
}
