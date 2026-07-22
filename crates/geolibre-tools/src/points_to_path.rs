//! GeoLibre tool: connect points into path (polyline) features.
//!
//! Pure-Rust counterpart of QGIS's *Points to Path* (Vector creation) — also the
//! analysis staple behind ArcGIS's *Points To Line* (Data Management): order a
//! point layer by a field and join the points into polylines, one line per group.
//!
//! The bundled suite has no general points-to-line converter: `points_along_lines`
//! goes the other way (line → points), `csv_points_to_vector` only reads points,
//! and the repo's `reconstruct_tracks` is time-specific (it *requires* a timestamp
//! field and adds gap-splitting and dwell detection). This tool is the plain,
//! general operation:
//!
//! - `order_field` sets the connection sequence within each path (numeric fields
//!   sort numerically, text lexicographically, or by `natural_sort` so `a9 < a10`).
//!   Omit it to keep the input feature order.
//! - `group_field` partitions the points into separate paths — one line per
//!   distinct value. Omit it to build a single path through every point.
//! - `close_path` appends the first vertex to close each path into a ring.
//!
//! Each output line carries its `group` value, the `begin`/`end` order values,
//! and the vertex `num_points`. Non-point features are ignored; a `MultiPoint`
//! contributes each of its points. Output keeps the input CRS.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct PointsToPathTool;

impl Tool for PointsToPathTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "points_to_path",
            display_name: "Points To Path",
            summary: "Connect points into polyline features, ordered by a field and (optionally) grouped into one path per value, with optional natural sort and path closing.",
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
                    name: "order_field",
                    description: "Attribute that sets the point connection order within each path. Numeric values sort numerically, text lexicographically. If omitted, input feature order is used.",
                    required: false,
                },
                ToolParamSpec {
                    name: "group_field",
                    description: "Attribute that partitions points into separate paths (one line per distinct value). If omitted, all points form a single path.",
                    required: false,
                },
                ToolParamSpec {
                    name: "natural_sort",
                    description: "Sort text order values naturally so 'a9' comes before 'a10'. Default false.",
                    required: false,
                },
                ToolParamSpec {
                    name: "close_path",
                    description: "Append the first vertex to the end of each path, closing it into a ring. Default false.",
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
        parse_optional_bool(args, "natural_sort")?;
        parse_optional_bool(args, "close_path")?;
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
        let order_field = parse_optional_str(args, "order_field")?;
        let group_field = parse_optional_str(args, "group_field")?;
        let natural_sort = parse_optional_bool(args, "natural_sort")?.unwrap_or(false);
        let close_path = parse_optional_bool(args, "close_path")?.unwrap_or(false);

        let layer = load_input_layer(input)?;
        let schema = &layer.schema;
        let order_idx = match order_field {
            Some(name) => Some(schema.field_index(name).ok_or_else(|| {
                ToolError::Validation(format!("order_field '{name}' not found in input schema"))
            })?),
            None => None,
        };
        let group_idx = match group_field {
            Some(name) => Some(schema.field_index(name).ok_or_else(|| {
                ToolError::Validation(format!("group_field '{name}' not found in input schema"))
            })?),
            None => None,
        };

        // Collect points as (group key, order key, feature index, sub index, coord).
        let mut pts: Vec<PointRef> = Vec::new();
        for (fidx, feature) in layer.features.iter().enumerate() {
            let group = match group_idx {
                Some(gi) => feature
                    .attributes
                    .get(gi)
                    .map(field_value_string)
                    .unwrap_or_default(),
                None => String::new(),
            };
            let order = match order_idx {
                Some(oi) => OrderKey::from_value(feature.attributes.get(oi)),
                None => OrderKey::Num(fidx as f64), // input order
            };
            for (sub, c) in point_coords(feature.geometry.as_ref())
                .into_iter()
                .enumerate()
            {
                pts.push(PointRef {
                    group: group.clone(),
                    order: order.clone(),
                    fidx,
                    sub,
                    coord: c,
                });
            }
        }
        let input_points = pts.len();

        // Partition into groups (BTreeMap keeps output deterministic by group).
        let mut groups: BTreeMap<String, Vec<PointRef>> = BTreeMap::new();
        for p in pts {
            groups.entry(p.group.clone()).or_default().push(p);
        }

        ctx.progress.info(&format!(
            "{input_points} point(s) -> {} path(s)",
            groups.len()
        ));

        let mut out = Layer::new("paths");
        out.geom_type = Some(GeometryType::LineString);
        out.crs = layer.crs.clone();
        out.add_field(FieldDef::new("group", FieldType::Text));
        out.add_field(FieldDef::new("begin", FieldType::Text));
        out.add_field(FieldDef::new("end", FieldType::Text));
        out.add_field(FieldDef::new("num_points", FieldType::Integer));

        let mut skipped = 0usize;
        for (group, mut members) in groups {
            // Stable sort by (order key, original feature index, sub index).
            members.sort_by(|a, b| {
                a.order
                    .cmp(&b.order, natural_sort)
                    .then(a.fidx.cmp(&b.fidx))
                    .then(a.sub.cmp(&b.sub))
            });
            if members.len() < 2 {
                skipped += 1;
                continue; // a path needs at least two vertices
            }
            let mut coords: Vec<Coord> = members.iter().map(|m| m.coord.clone()).collect();
            let begin = members.first().unwrap().order.label();
            let end = members.last().unwrap().order.label();
            if close_path {
                coords.push(coords[0].clone());
            }
            let n = coords.len() as i64;
            out.add_feature(
                Some(Geometry::LineString(coords)),
                &[
                    ("group", FieldValue::Text(group)),
                    ("begin", FieldValue::Text(begin)),
                    ("end", FieldValue::Text(end)),
                    ("num_points", FieldValue::Integer(n)),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed adding path: {e}")))?;
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

/// One point ready to be ordered into a path.
struct PointRef {
    group: String,
    order: OrderKey,
    fidx: usize,
    sub: usize,
    coord: Coord,
}

/// A sortable order value: numeric values compare by magnitude and sort before
/// text values (which compare lexicographically or naturally).
#[derive(Clone)]
enum OrderKey {
    Num(f64),
    Text(String),
}

impl OrderKey {
    fn from_value(v: Option<&FieldValue>) -> Self {
        match v {
            Some(fv) => match fv.as_f64() {
                Some(n) => OrderKey::Num(n),
                None => OrderKey::Text(field_value_string(fv)),
            },
            None => OrderKey::Text(String::new()),
        }
    }

    /// The value as text (for the begin/end attributes).
    fn label(&self) -> String {
        match self {
            OrderKey::Num(n) => {
                if n.fract() == 0.0 && n.abs() < 1e15 {
                    format!("{}", *n as i64)
                } else {
                    format!("{n}")
                }
            }
            OrderKey::Text(s) => s.clone(),
        }
    }

    fn cmp(&self, other: &OrderKey, natural: bool) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        match (self, other) {
            (OrderKey::Num(a), OrderKey::Num(b)) => a.partial_cmp(b).unwrap_or(Ordering::Equal),
            (OrderKey::Num(_), OrderKey::Text(_)) => Ordering::Less,
            (OrderKey::Text(_), OrderKey::Num(_)) => Ordering::Greater,
            (OrderKey::Text(a), OrderKey::Text(b)) => {
                if natural {
                    natural_cmp(a, b)
                } else {
                    a.cmp(b)
                }
            }
        }
    }
}

/// Natural comparison: splits each string into numeric and non-numeric runs and
/// compares them chunk by chunk, so "a9" < "a10".
fn natural_cmp(a: &str, b: &str) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    let (mut ai, mut bi) = (a.chars().peekable(), b.chars().peekable());
    loop {
        match (ai.peek().copied(), bi.peek().copied()) {
            (None, None) => return Ordering::Equal,
            (None, Some(_)) => return Ordering::Less,
            (Some(_), None) => return Ordering::Greater,
            (Some(ca), Some(cb)) => {
                if ca.is_ascii_digit() && cb.is_ascii_digit() {
                    let na = take_number(&mut ai);
                    let nb = take_number(&mut bi);
                    match na.cmp(&nb) {
                        Ordering::Equal => continue,
                        o => return o,
                    }
                } else {
                    match ca.cmp(&cb) {
                        Ordering::Equal => {
                            ai.next();
                            bi.next();
                        }
                        o => return o,
                    }
                }
            }
        }
    }
}

fn take_number(it: &mut std::iter::Peekable<std::str::Chars>) -> u128 {
    let mut n: u128 = 0;
    while let Some(&c) = it.peek() {
        if let Some(d) = c.to_digit(10) {
            n = n.saturating_mul(10).saturating_add(d as u128);
            it.next();
        } else {
            break;
        }
    }
    n
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
        let out = PointsToPathTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn line_coords(layer: &Layer, fid: usize) -> Vec<(f64, f64)> {
        match &layer.features[fid].geometry {
            Some(Geometry::LineString(cs)) => cs.iter().map(|c| (c.x, c.y)).collect(),
            _ => vec![],
        }
    }

    /// Points ordered by a numeric field connect in that order, not file order.
    #[test]
    fn orders_by_numeric_field() {
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

        let (out, layer) = run(&input, json!({ "order_field": "seq" }));
        assert_eq!(out.outputs["feature_count"], json!(1));
        // x should increase 1,2,3 after ordering by seq.
        assert_eq!(
            line_coords(&layer, 0),
            vec![(1.0, 0.0), (2.0, 0.0), (3.0, 0.0)]
        );
        assert_eq!(
            layer.features[0]
                .get(&layer.schema, "num_points")
                .unwrap()
                .as_i64(),
            Some(3)
        );
    }

    /// A group field partitions the points into one path per value.
    #[test]
    fn groups_into_separate_paths() {
        let mut layer = Layer::new("pts").with_geom_type(GeometryType::Point);
        layer.add_field(FieldDef::new("trip", FieldType::Text));
        layer.add_field(FieldDef::new("seq", FieldType::Integer));
        for (x, y, t, s) in [
            (0.0, 0.0, "A", 1i64),
            (1.0, 0.0, "A", 2),
            (0.0, 5.0, "B", 1),
            (1.0, 5.0, "B", 2),
            (2.0, 5.0, "B", 3),
        ] {
            layer
                .add_feature(
                    Some(Geometry::point(x, y)),
                    &[("trip", t.into()), ("seq", s.into())],
                )
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run(
            &input,
            json!({ "order_field": "seq", "group_field": "trip" }),
        );
        assert_eq!(out.outputs["feature_count"], json!(2));
        let np = |fid: usize| {
            layer.features[fid]
                .get(&layer.schema, "num_points")
                .unwrap()
                .as_i64()
                .unwrap()
        };
        // Group A -> 2 vertices, group B -> 3 (BTreeMap orders groups A, B).
        assert_eq!(np(0), 2);
        assert_eq!(np(1), 3);
        assert_eq!(
            layer.features[0].get(&layer.schema, "group").unwrap(),
            &FieldValue::Text("A".into())
        );
    }

    /// `close_path` appends the first vertex to close the ring.
    #[test]
    fn close_path_closes_the_ring() {
        let mut layer = Layer::new("pts").with_geom_type(GeometryType::Point);
        for (x, y) in [(0.0, 0.0), (2.0, 0.0), (1.0, 2.0)] {
            layer.add_feature(Some(Geometry::point(x, y)), &[]).unwrap();
        }
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (_, layer) = run(&input, json!({ "close_path": true }));
        let cs = line_coords(&layer, 0);
        assert_eq!(cs.len(), 4, "3 points + closing vertex");
        assert_eq!(cs.first(), cs.last(), "path is closed");
    }

    /// Natural sort orders text order values with embedded numbers correctly.
    #[test]
    fn natural_sort_handles_embedded_numbers() {
        let mut layer = Layer::new("pts").with_geom_type(GeometryType::Point);
        layer.add_field(FieldDef::new("label", FieldType::Text));
        // Lexicographic order would be a1, a10, a2; natural is a1, a2, a10.
        for (x, l) in [(0.0, "a1"), (2.0, "a10"), (1.0, "a2")] {
            layer
                .add_feature(Some(Geometry::point(x, 0.0)), &[("label", l.into())])
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (_, layer) = run(
            &input,
            json!({ "order_field": "label", "natural_sort": true }),
        );
        // Natural order -> x increases 0,1,2 (a1,a2,a10).
        assert_eq!(
            line_coords(&layer, 0),
            vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]
        );
    }

    /// Without an order field, points connect in input feature order.
    #[test]
    fn defaults_to_input_order() {
        let mut layer = Layer::new("pts").with_geom_type(GeometryType::Point);
        for (x, y) in [(5.0, 0.0), (3.0, 0.0), (9.0, 0.0)] {
            layer.add_feature(Some(Geometry::point(x, y)), &[]).unwrap();
        }
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (_, layer) = run(&input, json!({}));
        assert_eq!(
            line_coords(&layer, 0),
            vec![(5.0, 0.0), (3.0, 0.0), (9.0, 0.0)]
        );
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = PointsToPathTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input");
        assert!(bad(json!({ "input": "x.geojson" })).is_ok());
        assert!(bad(json!({ "input": "x.geojson", "close_path": "maybe" })).is_err());
    }

    /// A missing order/group field is a clear error.
    #[test]
    fn missing_field_errors() {
        let mut layer = Layer::new("pts").with_geom_type(GeometryType::Point);
        layer
            .add_feature(Some(Geometry::point(0.0, 0.0)), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "order_field": "nope" })).unwrap();
        assert!(PointsToPathTool.run(&args, &ctx()).is_err());
    }
}
