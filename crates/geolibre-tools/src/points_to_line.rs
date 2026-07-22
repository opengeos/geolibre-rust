//! GeoLibre tool: build polylines from ordered point features.
//!
//! Pure-Rust counterpart of ArcGIS's *Points To Line* (Data Management): take a
//! point layer and connect the points into polylines — one line per distinct
//! `line_field` value, with the vertices of each line ordered by `sort_field`,
//! and an optional `close_line` flag that returns each line to its first vertex.
//!
//! This is the ArcGIS-idiomatic sibling of the repo's `points_to_path` (which
//! follows the QGIS *Points to Path* naming/semantics with `group_field`,
//! `order_field`, `natural_sort`, `close_path`). The two share the connect-points
//! idea but differ in parameters and output contract:
//!
//! - Parameters use the ArcGIS names `line_field` / `sort_field` / `close_line`.
//! - The output line **carries the `line_field` value under its own field name**
//!   (as ArcGIS does), rather than a generic `group` column, plus a `point_count`.
//! - The whitebox suite has no equivalent: `points_along_lines` is the inverse
//!   (line -> points) and `create_routes` merges pre-existing *line* features.
//!
//! With no `line_field`, all points form a single line. With no `sort_field`,
//! the input feature order is used. `MultiPoint` features contribute each of
//! their points; non-point geometry is ignored. Output keeps the input CRS.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct PointsToLineTool;

/// Fallback output field name when the input has no `line_field`.
const DEFAULT_LINE_FIELD: &str = "line_id";

impl Tool for PointsToLineTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "points_to_line",
            display_name: "Points To Line",
            summary: "Build polyline features from points: one line per distinct line_field value, vertices ordered by sort_field, with an optional close-line flag.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector file path (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output line vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "line_field",
                    description: "Attribute whose distinct values partition the points into separate lines (one polyline per value). If omitted, all points form a single line.",
                    required: false,
                },
                ToolParamSpec {
                    name: "sort_field",
                    description: "Attribute that sets the vertex order within each line. Numeric values sort numerically, text lexicographically. If omitted, input feature order is used.",
                    required: false,
                },
                ToolParamSpec {
                    name: "close_line",
                    description: "Append the first vertex to the end of each line, closing it into a ring. Default false.",
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
        parse_optional_bool(args, "close_line")?;
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
        let line_field = parse_optional_str(args, "line_field")?;
        let sort_field = parse_optional_str(args, "sort_field")?;
        let close_line = parse_optional_bool(args, "close_line")?.unwrap_or(false);

        let layer = load_input_layer(input)?;
        let schema = &layer.schema;
        let line_idx = match line_field {
            Some(name) => Some(schema.field_index(name).ok_or_else(|| {
                ToolError::Validation(format!("line_field '{name}' not found in input schema"))
            })?),
            None => None,
        };
        let sort_idx = match sort_field {
            Some(name) => Some(schema.field_index(name).ok_or_else(|| {
                ToolError::Validation(format!("sort_field '{name}' not found in input schema"))
            })?),
            None => None,
        };

        // Collect points as (line key, sort key, feature index, sub index, coord).
        let mut pts: Vec<PointRef> = Vec::new();
        for (fidx, feature) in layer.features.iter().enumerate() {
            let line = match line_idx {
                Some(li) => feature
                    .attributes
                    .get(li)
                    .map(field_value_string)
                    .unwrap_or_default(),
                None => String::new(),
            };
            let sort = match sort_idx {
                Some(si) => SortKey::from_value(feature.attributes.get(si)),
                None => SortKey::Num(fidx as f64), // input order
            };
            for (sub, c) in point_coords(feature.geometry.as_ref())
                .into_iter()
                .enumerate()
            {
                pts.push(PointRef {
                    line: line.clone(),
                    sort: sort.clone(),
                    fidx,
                    sub,
                    coord: c,
                });
            }
        }
        let input_points = pts.len();

        // Partition into lines (BTreeMap keeps output deterministic by key).
        let mut lines: BTreeMap<String, Vec<PointRef>> = BTreeMap::new();
        for p in pts {
            lines.entry(p.line.clone()).or_default().push(p);
        }

        ctx.progress.info(&format!(
            "{input_points} point(s) -> {} line(s)",
            lines.len()
        ));

        // Output carries the line_field value under its own name (ArcGIS-style),
        // or a default field name when no line_field was supplied.
        let out_field = line_field.unwrap_or(DEFAULT_LINE_FIELD);
        let mut out = Layer::new("lines");
        out.geom_type = Some(GeometryType::LineString);
        out.crs = layer.crs.clone();
        out.add_field(FieldDef::new(out_field, FieldType::Text));
        out.add_field(FieldDef::new("point_count", FieldType::Integer));

        let mut skipped = 0usize;
        for (line_val, mut members) in lines {
            // Stable sort by (sort key, original feature index, sub index).
            members.sort_by(|a, b| {
                a.sort
                    .cmp(&b.sort)
                    .then(a.fidx.cmp(&b.fidx))
                    .then(a.sub.cmp(&b.sub))
            });
            if members.len() < 2 {
                skipped += 1;
                continue; // a line needs at least two vertices
            }
            let mut coords: Vec<Coord> = members.iter().map(|m| m.coord.clone()).collect();
            if close_line {
                coords.push(coords[0].clone());
            }
            let n = coords.len() as i64;
            out.add_feature(
                Some(Geometry::LineString(coords)),
                &[
                    (out_field, FieldValue::Text(line_val)),
                    ("point_count", FieldValue::Integer(n)),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed adding line: {e}")))?;
        }

        if skipped > 0 {
            ctx.progress.info(&format!(
                "{skipped} group(s) had <2 points and were skipped"
            ));
        }

        let feature_count = out.len();
        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_points".to_string(), json!(input_points));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("skipped_groups".to_string(), json!(skipped));
        Ok(ToolRunResult { outputs })
    }
}

/// One point ready to be ordered into a line.
struct PointRef {
    line: String,
    sort: SortKey,
    fidx: usize,
    sub: usize,
    coord: Coord,
}

/// A sortable vertex-order value: numeric values compare by magnitude and sort
/// before text values (which compare lexicographically).
#[derive(Clone)]
enum SortKey {
    Num(f64),
    Text(String),
}

impl SortKey {
    fn from_value(v: Option<&FieldValue>) -> Self {
        match v {
            Some(fv) => match fv.as_f64() {
                Some(n) => SortKey::Num(n),
                None => SortKey::Text(field_value_string(fv)),
            },
            None => SortKey::Text(String::new()),
        }
    }

    fn cmp(&self, other: &SortKey) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        match (self, other) {
            (SortKey::Num(a), SortKey::Num(b)) => a.partial_cmp(b).unwrap_or(Ordering::Equal),
            (SortKey::Num(_), SortKey::Text(_)) => Ordering::Less,
            (SortKey::Text(_), SortKey::Num(_)) => Ordering::Greater,
            (SortKey::Text(a), SortKey::Text(b)) => a.cmp(b),
        }
    }
}

/// Extracts point coordinates from a geometry: a `Point` yields one, a
/// `MultiPoint` yields each; other geometry types yield none.
fn point_coords(geom: Option<&Geometry>) -> Vec<Coord> {
    match geom {
        Some(Geometry::Point(c)) => vec![c.clone()],
        Some(Geometry::MultiPoint(cs)) => cs.clone(),
        _ => vec![],
    }
}

fn field_value_string(v: &FieldValue) -> String {
    match v {
        FieldValue::Text(s) | FieldValue::Date(s) | FieldValue::DateTime(s) => s.clone(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) => f.to_string(),
        FieldValue::Boolean(b) => b.to_string(),
        FieldValue::Null | FieldValue::Blob(_) => String::new(),
    }
}

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn run(input: &str, args: serde_json::Value) -> (ToolRunResult, Layer) {
        let mut m = args.as_object().unwrap().clone();
        m.insert("input".to_string(), json!(input));
        let args: ToolArgs = serde_json::from_value(Value::Object(m)).unwrap();
        let out = PointsToLineTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn line_coords(layer: &Layer, fid: usize) -> Vec<(f64, f64)> {
        match &layer.features[fid].geometry {
            Some(Geometry::LineString(cs)) => cs.iter().map(|c| (c.x, c.y)).collect(),
            _ => vec![],
        }
    }

    /// One line per distinct line_field value; vertex count == group size.
    #[test]
    fn one_line_per_group() {
        let mut layer = Layer::new("pts").with_geom_type(GeometryType::Point);
        layer.add_field(FieldDef::new("route", FieldType::Text));
        layer.add_field(FieldDef::new("seq", FieldType::Integer));
        for (x, y, r, s) in [
            (0.0, 0.0, "A", 1i64),
            (1.0, 0.0, "A", 2),
            (0.0, 5.0, "B", 1),
            (1.0, 5.0, "B", 2),
            (2.0, 5.0, "B", 3),
        ] {
            layer
                .add_feature(
                    Some(Geometry::point(x, y)),
                    &[("route", r.into()), ("seq", s.into())],
                )
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run(
            &input,
            json!({ "line_field": "route", "sort_field": "seq" }),
        );
        assert_eq!(out.outputs["feature_count"], json!(2));
        let pc = |fid: usize| {
            layer.features[fid]
                .get(&layer.schema, "point_count")
                .unwrap()
                .as_i64()
                .unwrap()
        };
        // Group A -> 2 vertices, group B -> 3 (BTreeMap orders keys A, B).
        assert_eq!(pc(0), 2);
        assert_eq!(pc(1), 3);
        // The line_field value is carried under its own field name.
        assert_eq!(
            layer.features[0].get(&layer.schema, "route").unwrap(),
            &FieldValue::Text("A".into())
        );
    }

    /// Vertices are ordered by sort_field, not by input file order.
    #[test]
    fn orders_vertices_by_sort_field() {
        let mut layer = Layer::new("pts").with_geom_type(GeometryType::Point);
        layer.add_field(FieldDef::new("seq", FieldType::Integer));
        // Added out of order: seq 3, 1, 2.
        for (x, s) in [(3.0, 3i64), (1.0, 1), (2.0, 2)] {
            layer
                .add_feature(Some(Geometry::point(x, 0.0)), &[("seq", s.into())])
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run(&input, json!({ "sort_field": "seq" }));
        assert_eq!(out.outputs["feature_count"], json!(1));
        // x should increase 1,2,3 after ordering by seq.
        assert_eq!(
            line_coords(&layer, 0),
            vec![(1.0, 0.0), (2.0, 0.0), (3.0, 0.0)]
        );
    }

    /// `close_line` appends the first vertex, adding one to the count.
    #[test]
    fn close_line_closes_the_ring() {
        let mut layer = Layer::new("pts").with_geom_type(GeometryType::Point);
        for (x, y) in [(0.0, 0.0), (2.0, 0.0), (1.0, 2.0)] {
            layer.add_feature(Some(Geometry::point(x, y)), &[]).unwrap();
        }
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (_, layer) = run(&input, json!({ "close_line": true }));
        let cs = line_coords(&layer, 0);
        assert_eq!(cs.len(), 4, "3 points + closing vertex");
        assert_eq!(cs.first(), cs.last(), "line is closed");
        assert_eq!(
            layer.features[0]
                .get(&layer.schema, "point_count")
                .unwrap()
                .as_i64(),
            Some(4)
        );
    }

    /// Without a line_field, all points connect into a single line in input order.
    #[test]
    fn single_line_in_input_order() {
        let mut layer = Layer::new("pts").with_geom_type(GeometryType::Point);
        for (x, y) in [(5.0, 0.0), (3.0, 0.0), (9.0, 0.0)] {
            layer.add_feature(Some(Geometry::point(x, y)), &[]).unwrap();
        }
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run(&input, json!({}));
        assert_eq!(out.outputs["feature_count"], json!(1));
        assert_eq!(
            line_coords(&layer, 0),
            vec![(5.0, 0.0), (3.0, 0.0), (9.0, 0.0)]
        );
    }

    /// Groups with a single point cannot form a line and are skipped.
    #[test]
    fn skips_single_point_groups() {
        let mut layer = Layer::new("pts").with_geom_type(GeometryType::Point);
        layer.add_field(FieldDef::new("route", FieldType::Text));
        for (x, y, r) in [(0.0, 0.0, "A"), (1.0, 0.0, "A"), (9.0, 9.0, "solo")] {
            layer
                .add_feature(Some(Geometry::point(x, y)), &[("route", r.into())])
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, _) = run(&input, json!({ "line_field": "route" }));
        assert_eq!(out.outputs["feature_count"], json!(1));
        assert_eq!(out.outputs["skipped_groups"], json!(1));
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = PointsToLineTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input");
        assert!(bad(json!({ "input": "x.geojson" })).is_ok());
        assert!(bad(json!({ "input": "x.geojson", "close_line": "maybe" })).is_err());
    }

    /// A missing line/sort field is a clear error.
    #[test]
    fn missing_field_errors() {
        let mut layer = Layer::new("pts").with_geom_type(GeometryType::Point);
        layer
            .add_feature(Some(Geometry::point(0.0, 0.0)), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "sort_field": "nope" })).unwrap();
        assert!(PointsToLineTool.run(&args, &ctx()).is_err());
    }
}
