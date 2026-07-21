//! GeoLibre tool: detect and remove duplicate vector features.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Find Identical* / *Delete Identical*
//! (Data Management). The bundled `remove_duplicates` is LiDAR-point-cloud only
//! and `eliminate_coincident_points` is points-only; nothing deduplicates
//! arbitrary vector features, a routine chore after merges and appends.
//!
//! Every feature gets a canonical key built from its geometry (all vertices,
//! optionally snapped to `xy_tolerance`) and/or the values of chosen `fields`.
//! Features with equal keys form a duplicate group. In `report` mode the output
//! copies the input and adds `dup_group` (a stable group id) and `dup_seq` (0
//! for the first feature of a group, 1.. for the rest); in `delete` mode only
//! the first feature of each group is kept.

use std::collections::{BTreeMap, HashMap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

#[derive(Clone, Copy, PartialEq)]
enum Mode {
    Report,
    Delete,
}

pub struct FindIdenticalTool;

impl Tool for FindIdenticalTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "find_identical",
            display_name: "Find Identical",
            summary: "Find features that are identical on geometry and/or chosen fields, grouping duplicates; 'report' adds dup_group/dup_seq columns, 'delete' keeps only the first of each group, like ArcGIS Find Identical / Delete Identical.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer to deduplicate.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "fields",
                    description: "Comma-separated fields to compare. Omit to compare geometry only; combine with 'compare_geometry'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "compare_geometry",
                    description: "Whether geometry is part of the identity test (default true; set false to compare fields only).",
                    required: false,
                },
                ToolParamSpec {
                    name: "xy_tolerance",
                    description: "Snap tolerance for geometry comparison (vertices are rounded to this grid before hashing). Default 0 (exact).",
                    required: false,
                },
                ToolParamSpec {
                    name: "mode",
                    description: "'report' (default; adds dup_group/dup_seq) or 'delete' (keeps the first of each group).",
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
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args.get("input").and_then(Value::as_str).unwrap();
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;

        let field_idx: Vec<usize> = prm
            .fields
            .iter()
            .map(|f| {
                layer
                    .schema
                    .field_index(f)
                    .ok_or_else(|| ToolError::Validation(format!("field '{f}' not found")))
            })
            .collect::<Result<_, _>>()?;

        if !prm.compare_geometry && field_idx.is_empty() {
            return Err(ToolError::Validation(
                "nothing to compare: enable compare_geometry or provide 'fields'".to_string(),
            ));
        }

        ctx.progress
            .info(&format!("hashing {} feature(s)", layer.features.len()));

        // Assign each feature to a group by canonical key, in first-seen order.
        let mut key_to_group: HashMap<String, usize> = HashMap::new();
        let mut group_of: Vec<usize> = Vec::with_capacity(layer.features.len());
        let mut seq_of: Vec<usize> = Vec::with_capacity(layer.features.len());
        let mut group_sizes: Vec<usize> = Vec::new();
        for feature in &layer.features {
            let mut key = String::new();
            if prm.compare_geometry {
                geometry_key(feature.geometry.as_ref(), prm.xy_tolerance, &mut key);
                key.push('|');
            }
            for &i in &field_idx {
                key.push_str(&field_key(&feature.attributes[i]));
                key.push('\u{1}');
            }
            let g = *key_to_group.entry(key).or_insert_with(|| {
                group_sizes.push(0);
                group_sizes.len() - 1
            });
            let seq = group_sizes[g];
            group_sizes[g] += 1;
            group_of.push(g);
            seq_of.push(seq);
        }

        let group_count = group_sizes.len();
        let duplicate_count = layer.features.len() - group_count;

        // Build output.
        let mut out = Layer::new("find_identical");
        if let Some(gt) = layer.geom_type {
            out = out.with_geom_type(gt);
        }
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for field in layer.schema.fields() {
            out.add_field(field.clone());
        }
        if prm.mode == Mode::Report {
            out.add_field(FieldDef::new("dup_group", FieldType::Integer));
            out.add_field(FieldDef::new("dup_seq", FieldType::Integer));
        }

        let mut kept = 0usize;
        for (fidx, feature) in layer.features.iter().enumerate() {
            if prm.mode == Mode::Delete && seq_of[fidx] != 0 {
                continue; // keep only the first of each group
            }
            let mut attrs = feature.attributes.clone();
            if prm.mode == Mode::Report {
                attrs.push(FieldValue::Integer(group_of[fidx] as i64));
                attrs.push(FieldValue::Integer(seq_of[fidx] as i64));
            }
            out.push(wbvector::Feature {
                fid: 0,
                geometry: feature.geometry.clone(),
                attributes: attrs,
            });
            kept += 1;
        }

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_count".to_string(), json!(layer.features.len()));
        outputs.insert("group_count".to_string(), json!(group_count));
        outputs.insert("duplicate_count".to_string(), json!(duplicate_count));
        outputs.insert("output_count".to_string(), json!(kept));
        Ok(ToolRunResult { outputs })
    }
}

// ── Canonical keys ────────────────────────────────────────────────────────────

fn snap(v: f64, tol: f64) -> i64 {
    if tol > 0.0 {
        (v / tol).round() as i64
    } else {
        v.to_bits() as i64
    }
}

/// Append a canonical, order-preserving encoding of the geometry to `key`.
fn geometry_key(geom: Option<&Geometry>, tol: f64, key: &mut String) {
    let Some(geom) = geom else {
        key.push_str("NULL");
        return;
    };
    let coord = |c: &Coord, key: &mut String| {
        key.push_str(&snap(c.x, tol).to_string());
        key.push(',');
        key.push_str(&snap(c.y, tol).to_string());
        key.push(';');
    };
    let ring = |r: &Ring, key: &mut String| {
        key.push('(');
        for c in r.coords() {
            coord(c, key);
        }
        key.push(')');
    };
    match geom {
        Geometry::Point(c) => {
            key.push('P');
            coord(c, key);
        }
        Geometry::MultiPoint(cs) => {
            key.push_str("MP");
            for c in cs {
                coord(c, key);
            }
        }
        Geometry::LineString(cs) => {
            key.push('L');
            for c in cs {
                coord(c, key);
            }
        }
        Geometry::MultiLineString(parts) => {
            key.push_str("ML");
            for cs in parts {
                key.push('(');
                for c in cs {
                    coord(c, key);
                }
                key.push(')');
            }
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            key.push_str("PG");
            ring(exterior, key);
            for h in interiors {
                ring(h, key);
            }
        }
        Geometry::MultiPolygon(parts) => {
            key.push_str("MPG");
            for (ext, holes) in parts {
                key.push('[');
                ring(ext, key);
                for h in holes {
                    ring(h, key);
                }
                key.push(']');
            }
        }
        Geometry::GeometryCollection(items) => {
            key.push_str("GC");
            for g in items {
                geometry_key(Some(g), tol, key);
                key.push('/');
            }
        }
    }
}

fn field_key(fv: &FieldValue) -> String {
    if let Some(i) = fv.as_i64() {
        format!("i{i}")
    } else if let Some(f) = fv.as_f64() {
        format!("f{f}")
    } else {
        format!("s{}", fv.as_str().unwrap_or(""))
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    fields: Vec<String>,
    compare_geometry: bool,
    xy_tolerance: f64,
    mode: Mode,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let fields: Vec<String> = parse_optional_str(args, "fields")?
        .map(|s| {
            s.split(',')
                .map(str::trim)
                .filter(|f| !f.is_empty())
                .map(str::to_string)
                .collect()
        })
        .unwrap_or_default();
    let compare_geometry = match args.get("compare_geometry") {
        None | Some(Value::Null) => true,
        Some(Value::Bool(b)) => *b,
        Some(Value::String(s)) => !matches!(s.trim().to_lowercase().as_str(), "false" | "0" | "no"),
        Some(_) => true,
    };
    let xy_tolerance = match args.get("xy_tolerance") {
        None | Some(Value::Null) => 0.0,
        Some(Value::Number(n)) => n.as_f64().unwrap_or(0.0).max(0.0),
        Some(Value::String(s)) if s.trim().is_empty() => 0.0,
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'xy_tolerance' must be a number".into()))?
            .max(0.0),
        Some(_) => {
            return Err(ToolError::Validation(
                "'xy_tolerance' must be a number".into(),
            ))
        }
    };
    let mode = match parse_optional_str(args, "mode")?.map(|s| s.trim().to_lowercase()) {
        None => Mode::Report,
        Some(s) if s.is_empty() || s == "report" => Mode::Report,
        Some(s) if s == "delete" => Mode::Delete,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "'mode' must be 'report' or 'delete', got '{other}'"
            )))
        }
    };
    Ok(Params {
        fields,
        compare_geometry,
        xy_tolerance,
        mode,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn layer_of(pts: &[(f64, f64, &str)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        for (x, y, n) in pts {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(*x, *y))),
                &[("name", (*n).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = FindIdenticalTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// Geometry-only: two coincident points are one group; report flags them.
    #[test]
    fn geometry_duplicates_grouped() {
        let input = layer_of(&[(0.0, 0.0, "a"), (0.0, 0.0, "b"), (5.0, 5.0, "c")]);
        let (out, layer) = run(json!({ "input": input }));
        assert_eq!(out.outputs["input_count"], json!(3));
        assert_eq!(out.outputs["group_count"], json!(2));
        assert_eq!(out.outputs["duplicate_count"], json!(1));
        let gi = layer.schema.field_index("dup_group").unwrap();
        let si = layer.schema.field_index("dup_seq").unwrap();
        let groups: Vec<i64> = layer
            .iter()
            .map(|f| f.attributes[gi].as_i64().unwrap())
            .collect();
        assert_eq!(groups[0], groups[1], "coincident points share a group");
        assert_ne!(groups[0], groups[2]);
        // dup_seq is 0 then 1 for the coincident pair.
        assert_eq!(
            layer.iter().nth(1).unwrap().attributes[si]
                .as_i64()
                .unwrap(),
            1
        );
    }

    /// delete mode keeps one feature per group.
    #[test]
    fn delete_keeps_first_per_group() {
        let input = layer_of(&[(0.0, 0.0, "a"), (0.0, 0.0, "b"), (5.0, 5.0, "c")]);
        let (out, layer) = run(json!({ "input": input, "mode": "delete" }));
        assert_eq!(out.outputs["output_count"], json!(2));
        assert_eq!(layer.iter().count(), 2);
        // No dup_* columns in delete mode.
        assert!(layer.schema.field_index("dup_group").is_none());
    }

    /// Field-only comparison ignores geometry.
    #[test]
    fn field_only_comparison() {
        // Distinct geometries but two share name "x".
        let input = layer_of(&[(0.0, 0.0, "x"), (9.0, 9.0, "x"), (5.0, 5.0, "y")]);
        let (out, _l) = run(json!({
            "input": input, "fields": "name", "compare_geometry": false
        }));
        assert_eq!(out.outputs["group_count"], json!(2), "two distinct names");
        assert_eq!(out.outputs["duplicate_count"], json!(1));
    }

    /// Geometry + fields: same location but different name are NOT identical.
    #[test]
    fn geometry_and_fields_combined() {
        let input = layer_of(&[(0.0, 0.0, "a"), (0.0, 0.0, "b")]);
        let (out, _l) = run(json!({ "input": input, "fields": "name" }));
        assert_eq!(
            out.outputs["group_count"],
            json!(2),
            "same xy, different name -> distinct"
        );
    }

    /// xy_tolerance snaps near-coincident points into one group.
    #[test]
    fn xy_tolerance_snaps() {
        let input = layer_of(&[(0.0, 0.0, "a"), (0.4, 0.3, "b")]);
        // exact -> 2 groups
        let (exact, _l) = run(json!({ "input": input }));
        assert_eq!(exact.outputs["group_count"], json!(2));
        // tolerance 1.0 -> both snap to (0,0) -> 1 group
        let (snapped, _l) = run(json!({ "input": input, "xy_tolerance": 1.0 }));
        assert_eq!(snapped.outputs["group_count"], json!(1));
    }

    #[test]
    fn rejects_missing_input() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(FindIdenticalTool.validate(&args).is_err());
    }

    #[test]
    fn rejects_nothing_to_compare() {
        let input = layer_of(&[(0.0, 0.0, "a")]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "compare_geometry": false })).unwrap();
        // No fields + no geometry comparison -> run errors.
        assert!(FindIdenticalTool.run(&args, &ctx()).is_err());
    }
}
