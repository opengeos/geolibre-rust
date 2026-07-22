//! GeoLibre tool: build m-enabled route features from line features.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Create Routes* (Linear Referencing):
//! dissolve input line features by a route-identifier field, stitch each route's
//! segments end-to-end into a single ordered polyline, and assign linear
//! *measures* (m-values) along it. The result is the entry point of the
//! linear-referencing pipeline — the m-enabled routes that event-location tools
//! (locate points/lines along routes, route calibration) consume.
//!
//! The bundled suite has no route builder: it can *locate* events on routes and
//! *recalibrate* existing measures, but nothing turns a plain road/rail/river
//! network into measured routes in the first place. This tool fills that gap.
//!
//! - `route_id_field` groups the input lines into routes (one output polyline per
//!   distinct value). Its value and type are copied onto the output feature.
//! - Each route's segments are ordered by a greedy nearest-endpoint walk that
//!   starts at a free end, so a road arriving as many small segments becomes one
//!   continuous, monotonically-measured line.
//! - `measure_source` decides how measures are assigned:
//!     * `LENGTH` — measures run 0..geometric-length along the route.
//!     * `ONE_FIELD` — the route's start measure is read from `from_measure_field`;
//!       the end measure is start + geometric length.
//!     * `TWO_FIELDS` — start/end measures are read from `from_measure_field` /
//!       `to_measure_field` and interpolated by cumulative length.
//! - `measure_factor` scales the natural (length-based) measure and
//!   `measure_offset` shifts every measure, matching the ArcGIS parameters.
//! - `ignore_gaps` (LENGTH/ONE_FIELD) omits the straight jumps between disjoint
//!   segments from the accumulated measure, so gaps do not inflate distances.
//!
//! Per-vertex m-values are written into the output geometry (OGC m-enabled
//! coordinates) and each route also carries `from_m`, `to_m`, `length`
//! (accumulated measure length, before the factor) and `n_segments`. Non-line
//! features are ignored; `MultiLineString` inputs are exploded into their parts.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CreateRoutesTool;

/// How measures are derived along a route.
#[derive(Clone, Copy, PartialEq)]
enum MeasureSource {
    Length,
    OneField,
    TwoFields,
}

impl Tool for CreateRoutesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "create_routes",
            display_name: "Create Routes",
            summary: "Dissolve line features by a route-identifier field, order their segments into a single polyline per route, and assign linear measures (m-values) from geometric length or from/to measure fields.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input line vector file path, format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "route_id_field",
                    description: "Attribute that identifies each route; input lines are grouped by it, one output route per distinct value.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output route vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "measure_source",
                    description: "How measures are assigned: LENGTH (0..geometric length), ONE_FIELD (start from a field, end = start + length), or TWO_FIELDS (start/end from fields, interpolated). Default LENGTH.",
                    required: false,
                },
                ToolParamSpec {
                    name: "from_measure_field",
                    description: "Numeric field giving each route's start measure (used by ONE_FIELD and TWO_FIELDS). Read from the route's first segment.",
                    required: false,
                },
                ToolParamSpec {
                    name: "to_measure_field",
                    description: "Numeric field giving each route's end measure (used by TWO_FIELDS). Read from the route's first segment.",
                    required: false,
                },
                ToolParamSpec {
                    name: "measure_factor",
                    description: "Multiplier applied to the natural (length-based) measure, e.g. 0.001 to turn metres into kilometres. Default 1.0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "measure_offset",
                    description: "Value added to every measure along the route (shifts the whole route's m-values). Default 0.0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "ignore_gaps",
                    description: "For LENGTH/ONE_FIELD, exclude the straight jumps between disjoint segments from the accumulated measure so gaps do not inflate distances. Default false.",
                    required: false,
                },
                ToolParamSpec {
                    name: "snap_tolerance",
                    description: "Endpoints within this distance (CRS units) are treated as coincident when stitching segments. Default 1e-6.",
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
        if args
            .get("route_id_field")
            .and_then(Value::as_str)
            .map(str::trim)
            .unwrap_or("")
            .is_empty()
        {
            return Err(ToolError::Validation(
                "missing required string parameter 'route_id_field'".to_string(),
            ));
        }
        parse_params(args)?;
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
        let route_id_field = args
            .get("route_id_field")
            .and_then(Value::as_str)
            .filter(|s| !s.trim().is_empty())
            .ok_or_else(|| {
                ToolError::Validation("missing required parameter 'route_id_field'".to_string())
            })?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let schema = &layer.schema;

        let rid_idx = schema.field_index(route_id_field).ok_or_else(|| {
            ToolError::Validation(format!(
                "route_id_field '{route_id_field}' not found in input schema"
            ))
        })?;
        let rid_type = schema.field(route_id_field).map(|f| f.field_type);

        // Resolve measure-field indices for the source that needs them.
        let from_idx = resolve_field(
            schema,
            prm.from_measure_field.as_deref(),
            "from_measure_field",
        )?;
        let to_idx = resolve_field(schema, prm.to_measure_field.as_deref(), "to_measure_field")?;
        match prm.measure_source {
            MeasureSource::OneField => {
                if from_idx.is_none() {
                    return Err(ToolError::Validation(
                        "measure_source ONE_FIELD requires 'from_measure_field'".to_string(),
                    ));
                }
            }
            MeasureSource::TwoFields => {
                if from_idx.is_none() || to_idx.is_none() {
                    return Err(ToolError::Validation(
                        "measure_source TWO_FIELDS requires both 'from_measure_field' and 'to_measure_field'".to_string(),
                    ));
                }
            }
            MeasureSource::Length => {}
        }

        // Group segments by route id, preserving first-seen order for determinism.
        // Each group keeps its segments' coords and the source feature index of the
        // first segment (whose attributes seed the route).
        let mut order: Vec<String> = Vec::new();
        let mut groups: BTreeMap<String, RouteGroup> = BTreeMap::new();
        let mut input_lines = 0usize;
        for (fidx, feature) in layer.features.iter().enumerate() {
            let parts: Vec<Vec<Coord>> = match &feature.geometry {
                Some(Geometry::LineString(cs)) if cs.len() >= 2 => vec![cs.clone()],
                Some(Geometry::MultiLineString(ps)) => {
                    ps.iter().filter(|cs| cs.len() >= 2).cloned().collect()
                }
                _ => continue,
            };
            if parts.is_empty() {
                continue;
            }
            let key = feature
                .attributes
                .get(rid_idx)
                .map(field_value_string)
                .unwrap_or_default();
            let entry = groups.entry(key.clone()).or_insert_with(|| {
                order.push(key.clone());
                RouteGroup {
                    segments: Vec::new(),
                    first_feat: fidx,
                }
            });
            for p in parts {
                input_lines += 1;
                entry.segments.push(p);
            }
        }

        ctx.progress.info(&format!(
            "{input_lines} line segment(s) -> {} route(s)",
            groups.len()
        ));

        // Build the output layer.
        let mut out = Layer::new(layer.name.clone());
        out.geom_type = Some(GeometryType::LineString);
        out.crs = layer.crs.clone();
        out.add_field(FieldDef::new(
            route_id_field,
            rid_type.unwrap_or(FieldType::Text),
        ));
        out.add_field(FieldDef::new("from_m", FieldType::Float));
        out.add_field(FieldDef::new("to_m", FieldType::Float));
        out.add_field(FieldDef::new("length", FieldType::Float));
        out.add_field(FieldDef::new("n_segments", FieldType::Integer));

        for key in &order {
            let group = &groups[key];
            let n_segments = group.segments.len();

            // Order and stitch the route's segments into one polyline, tracking
            // which edges are gaps (jumps between disjoint segments).
            let (mut coords, gap_edge) = stitch(&group.segments, prm.snap_tolerance);
            if coords.len() < 2 {
                continue;
            }

            // Cumulative measure length per vertex (gaps optionally excluded).
            let mut cum = vec![0.0f64; coords.len()];
            for i in 1..coords.len() {
                let seg = dist(&coords[i - 1], &coords[i]);
                let add = if gap_edge[i] && prm.ignore_gaps {
                    0.0
                } else {
                    seg
                };
                cum[i] = cum[i - 1] + add;
            }
            let length = *cum.last().unwrap();

            // Read source measures for the route (from its first segment's feature).
            let feat = &layer.features[group.first_feat];
            let field_val = |idx: Option<usize>| -> f64 {
                idx.and_then(|i| feat.attributes.get(i))
                    .and_then(FieldValue::as_f64)
                    .unwrap_or(0.0)
            };

            // Resolve start/end measures and per-vertex m-values.
            let (from_m, to_m) = match prm.measure_source {
                MeasureSource::Length => {
                    let f = prm.measure_offset;
                    let t = prm.measure_offset + prm.measure_factor * length;
                    (f, t)
                }
                MeasureSource::OneField => {
                    let start = field_val(from_idx);
                    let f = start + prm.measure_offset;
                    let t = f + prm.measure_factor * length;
                    (f, t)
                }
                MeasureSource::TwoFields => {
                    let f = field_val(from_idx) * prm.measure_factor + prm.measure_offset;
                    let t = field_val(to_idx) * prm.measure_factor + prm.measure_offset;
                    (f, t)
                }
            };

            for (i, c) in coords.iter_mut().enumerate() {
                let m = match prm.measure_source {
                    MeasureSource::TwoFields => {
                        if length > 0.0 {
                            from_m + (to_m - from_m) * (cum[i] / length)
                        } else {
                            from_m
                        }
                    }
                    _ => from_m + prm.measure_factor * cum[i],
                };
                c.m = Some(m);
            }

            let rid_value = feat
                .attributes
                .get(rid_idx)
                .cloned()
                .unwrap_or(FieldValue::Null);

            out.add_feature(
                Some(Geometry::LineString(coords)),
                &[
                    (route_id_field, rid_value),
                    ("from_m", FieldValue::Float(from_m)),
                    ("to_m", FieldValue::Float(to_m)),
                    ("length", FieldValue::Float(length)),
                    ("n_segments", FieldValue::Integer(n_segments as i64)),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed adding route: {e}")))?;
        }

        let route_count = out.len();
        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_lines".to_string(), json!(input_lines));
        outputs.insert("route_count".to_string(), json!(route_count));
        Ok(ToolRunResult { outputs })
    }
}

/// Segments and seed feature for one route.
struct RouteGroup {
    segments: Vec<Vec<Coord>>,
    first_feat: usize,
}

/// Parameters after parsing/validation.
struct Params {
    measure_source: MeasureSource,
    from_measure_field: Option<String>,
    to_measure_field: Option<String>,
    measure_factor: f64,
    measure_offset: f64,
    ignore_gaps: bool,
    snap_tolerance: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let measure_source = match parse_optional_str(args, "measure_source")? {
        None => MeasureSource::Length,
        Some(s) => match s.trim().to_ascii_uppercase().as_str() {
            "LENGTH" => MeasureSource::Length,
            "ONE_FIELD" | "ONEFIELD" => MeasureSource::OneField,
            "TWO_FIELDS" | "TWOFIELDS" => MeasureSource::TwoFields,
            other => {
                return Err(ToolError::Validation(format!(
                    "parameter 'measure_source' must be LENGTH, ONE_FIELD or TWO_FIELDS (got '{other}')"
                )))
            }
        },
    };
    let from_measure_field = parse_optional_str(args, "from_measure_field")?.map(str::to_string);
    let to_measure_field = parse_optional_str(args, "to_measure_field")?.map(str::to_string);
    let measure_factor = parse_optional_f64(args, "measure_factor")?.unwrap_or(1.0);
    if !measure_factor.is_finite() {
        return Err(ToolError::Validation(
            "parameter 'measure_factor' must be a finite number".to_string(),
        ));
    }
    let measure_offset = parse_optional_f64(args, "measure_offset")?.unwrap_or(0.0);
    if !measure_offset.is_finite() {
        return Err(ToolError::Validation(
            "parameter 'measure_offset' must be a finite number".to_string(),
        ));
    }
    let ignore_gaps = parse_optional_bool(args, "ignore_gaps")?.unwrap_or(false);
    let snap_tolerance = parse_optional_f64(args, "snap_tolerance")?.unwrap_or(1e-6);
    if !(snap_tolerance >= 0.0 && snap_tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'snap_tolerance' must be a non-negative number".to_string(),
        ));
    }
    Ok(Params {
        measure_source,
        from_measure_field,
        to_measure_field,
        measure_factor,
        measure_offset,
        ignore_gaps,
        snap_tolerance,
    })
}

/// Resolves an optional field name to its schema index, erroring if named but
/// absent.
fn resolve_field(
    schema: &wbvector::Schema,
    name: Option<&str>,
    label: &str,
) -> Result<Option<usize>, ToolError> {
    match name {
        None => Ok(None),
        Some(n) => schema.field_index(n).map(Some).ok_or_else(|| {
            ToolError::Validation(format!("{label} '{n}' not found in input schema"))
        }),
    }
}

/// Orders a route's segments into one polyline by a greedy nearest-endpoint walk,
/// returning the stitched coordinates and, for each vertex, whether the edge that
/// enters it is a gap (a jump longer than `snap_tol` between disjoint segments).
///
/// The walk starts at the segment endpoint that is farthest from every other
/// segment's endpoints (a natural free end), then repeatedly appends the unused
/// segment whose nearer endpoint is closest to the current tail, orienting it so
/// its measures increase along the route.
fn stitch(segments: &[Vec<Coord>], snap_tol: f64) -> (Vec<Coord>, Vec<bool>) {
    if segments.is_empty() {
        return (Vec::new(), Vec::new());
    }
    if segments.len() == 1 {
        let coords = segments[0].clone();
        let gaps = vec![false; coords.len()];
        return (coords, gaps);
    }

    let mut used = vec![false; segments.len()];

    // Choose a starting segment/orientation: the endpoint whose distance to the
    // nearest *other* segment endpoint is largest is the most free end.
    let mut start_seg = 0usize;
    let mut start_reversed = false;
    let mut best_free = -1.0f64;
    for (i, seg) in segments.iter().enumerate() {
        for (which, ep) in [(false, &seg[0]), (true, seg.last().unwrap())] {
            let mut nearest = f64::INFINITY;
            for (j, other) in segments.iter().enumerate() {
                if j == i {
                    continue;
                }
                for oep in [&other[0], other.last().unwrap()] {
                    nearest = nearest.min(dist(ep, oep));
                }
            }
            if nearest > best_free {
                best_free = nearest;
                start_seg = i;
                start_reversed = which;
            }
        }
    }

    let mut coords: Vec<Coord> = if start_reversed {
        segments[start_seg].iter().rev().cloned().collect()
    } else {
        segments[start_seg].clone()
    };
    let mut gap_edge = vec![false; coords.len()];
    used[start_seg] = true;

    // Greedily append the closest remaining segment.
    for _ in 1..segments.len() {
        let tail = coords.last().unwrap().clone();
        let mut best: Option<(usize, bool, f64)> = None; // (seg, reversed, dist)
        for (j, seg) in segments.iter().enumerate() {
            if used[j] {
                continue;
            }
            let d_start = dist(&tail, &seg[0]);
            let d_end = dist(&tail, seg.last().unwrap());
            let (rev, d) = if d_start <= d_end {
                (false, d_start)
            } else {
                (true, d_end)
            };
            if best.map(|(_, _, bd)| d < bd).unwrap_or(true) {
                best = Some((j, rev, d));
            }
        }
        let (j, rev, d) = best.unwrap();
        used[j] = true;
        let seg: Vec<Coord> = if rev {
            segments[j].iter().rev().cloned().collect()
        } else {
            segments[j].clone()
        };
        if d <= snap_tol {
            // Contiguous: drop the duplicate connecting vertex.
            for c in seg.into_iter().skip(1) {
                coords.push(c);
                gap_edge.push(false);
            }
        } else {
            // Gap: keep the connecting vertex; the edge into it is a gap.
            let mut first = true;
            for c in seg {
                coords.push(c);
                gap_edge.push(first);
                first = false;
            }
        }
    }
    (coords, gap_edge)
}

fn dist(a: &Coord, b: &Coord) -> f64 {
    let dx = a.x - b.x;
    let dy = a.y - b.y;
    (dx * dx + dy * dy).sqrt()
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

fn parse_optional_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be a number"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a number"
        ))),
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

    fn line(pts: &[(f64, f64)]) -> Geometry {
        Geometry::LineString(pts.iter().map(|&(x, y)| Coord::xy(x, y)).collect())
    }

    fn run(input: &str, args: serde_json::Value) -> (ToolRunResult, Layer) {
        let mut m = args.as_object().unwrap().clone();
        m.insert("input".to_string(), json!(input));
        let args: ToolArgs = serde_json::from_value(Value::Object(m)).unwrap();
        let out = CreateRoutesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn getf(layer: &Layer, fid: usize, name: &str) -> f64 {
        layer.features[fid]
            .get(&layer.schema, name)
            .unwrap()
            .as_f64()
            .unwrap()
    }

    fn ms(layer: &Layer, fid: usize) -> Vec<f64> {
        match &layer.features[fid].geometry {
            Some(Geometry::LineString(cs)) => cs.iter().map(|c| c.m.unwrap()).collect(),
            _ => vec![],
        }
    }

    /// Two collinear segments of one route stitch into a single measured line
    /// whose end measure equals its geometric length (LENGTH source).
    #[test]
    fn length_measures_match_geometry() {
        let mut layer = Layer::new("roads");
        layer.add_field(FieldDef::new("rid", FieldType::Text));
        layer
            .add_feature(
                Some(line(&[(0.0, 0.0), (3.0, 0.0)])),
                &[("rid", "A".into())],
            )
            .unwrap();
        layer
            .add_feature(
                Some(line(&[(3.0, 0.0), (3.0, 4.0)])),
                &[("rid", "A".into())],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run(&input, json!({ "route_id_field": "rid" }));
        assert_eq!(out.outputs["route_count"], json!(1));
        assert_eq!(out.outputs["input_lines"], json!(2));
        // 3 + 4 = 7 total length.
        assert!((getf(&layer, 0, "length") - 7.0).abs() < 1e-9);
        assert!((getf(&layer, 0, "from_m") - 0.0).abs() < 1e-9);
        assert!((getf(&layer, 0, "to_m") - 7.0).abs() < 1e-9);
        assert_eq!(getf(&layer, 0, "n_segments"), 2.0);
        // Per-vertex m-values increase 0, 3, 7.
        let m = ms(&layer, 0);
        assert_eq!(m.first().copied(), Some(0.0));
        assert!((m.last().unwrap() - 7.0).abs() < 1e-9);
        for w in m.windows(2) {
            assert!(w[1] >= w[0], "measures must be monotonic");
        }
    }

    /// Distinct route ids produce one route each; total length is conserved.
    #[test]
    fn groups_by_route_id() {
        let mut layer = Layer::new("roads");
        layer.add_field(FieldDef::new("rid", FieldType::Text));
        for (pts, r) in [
            (vec![(0.0, 0.0), (2.0, 0.0)], "A"),
            (vec![(0.0, 5.0), (2.0, 5.0)], "B"),
            (vec![(2.0, 5.0), (5.0, 5.0)], "B"),
        ] {
            layer
                .add_feature(Some(line(&pts)), &[("rid", r.into())])
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run(&input, json!({ "route_id_field": "rid" }));
        assert_eq!(out.outputs["route_count"], json!(2));
        // A: length 2, B: length 5. Total input length 7 preserved.
        let total: f64 = (0..2).map(|i| getf(&layer, i, "length")).sum();
        assert!((total - 7.0).abs() < 1e-9);
    }

    /// measure_factor and measure_offset scale/shift the measures.
    #[test]
    fn factor_and_offset_apply() {
        let mut layer = Layer::new("roads");
        layer.add_field(FieldDef::new("rid", FieldType::Text));
        layer
            .add_feature(
                Some(line(&[(0.0, 0.0), (100.0, 0.0)])),
                &[("rid", "A".into())],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (_, layer) = run(
            &input,
            json!({ "route_id_field": "rid", "measure_factor": 0.001, "measure_offset": 10.0 }),
        );
        // 100 m -> 0.1 km, offset 10 -> from 10, to 10.1.
        assert!((getf(&layer, 0, "from_m") - 10.0).abs() < 1e-9);
        assert!((getf(&layer, 0, "to_m") - 10.1).abs() < 1e-9);
    }

    /// TWO_FIELDS reads start/end measures from fields and interpolates.
    #[test]
    fn two_fields_interpolate() {
        let mut layer = Layer::new("roads");
        layer.add_field(FieldDef::new("rid", FieldType::Text));
        layer.add_field(FieldDef::new("m0", FieldType::Float));
        layer.add_field(FieldDef::new("m1", FieldType::Float));
        layer
            .add_feature(
                Some(line(&[(0.0, 0.0), (2.0, 0.0), (4.0, 0.0)])),
                &[
                    ("rid", "A".into()),
                    ("m0", 100.0.into()),
                    ("m1", 200.0.into()),
                ],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (_, layer) = run(
            &input,
            json!({
                "route_id_field": "rid",
                "measure_source": "TWO_FIELDS",
                "from_measure_field": "m0",
                "to_measure_field": "m1",
            }),
        );
        assert!((getf(&layer, 0, "from_m") - 100.0).abs() < 1e-9);
        assert!((getf(&layer, 0, "to_m") - 200.0).abs() < 1e-9);
        // Midpoint vertex at half length -> measure 150.
        let m = ms(&layer, 0);
        assert!((m[1] - 150.0).abs() < 1e-9, "midpoint measure interpolated");
    }

    /// ignore_gaps omits the jump between disjoint segments from the measure.
    #[test]
    fn ignore_gaps_excludes_jump() {
        let mut layer = Layer::new("roads");
        layer.add_field(FieldDef::new("rid", FieldType::Text));
        // Two 2-unit segments separated by a 10-unit gap.
        layer
            .add_feature(
                Some(line(&[(0.0, 0.0), (2.0, 0.0)])),
                &[("rid", "A".into())],
            )
            .unwrap();
        layer
            .add_feature(
                Some(line(&[(12.0, 0.0), (14.0, 0.0)])),
                &[("rid", "A".into())],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (_, with) = run(
            &input,
            json!({ "route_id_field": "rid", "ignore_gaps": true }),
        );
        // Only the two 2-unit segments count: length 4, gap of 10 excluded.
        assert!((getf(&with, 0, "length") - 4.0).abs() < 1e-9);

        let (_, without) = run(&input, json!({ "route_id_field": "rid" }));
        // Gap of 10 included -> length 14.
        assert!((getf(&without, 0, "length") - 14.0).abs() < 1e-9);
    }

    /// Non-line features are ignored; a route needs at least one line.
    #[test]
    fn ignores_non_lines() {
        let mut layer = Layer::new("mixed");
        layer.add_field(FieldDef::new("rid", FieldType::Text));
        layer
            .add_feature(Some(Geometry::point(0.0, 0.0)), &[("rid", "A".into())])
            .unwrap();
        layer
            .add_feature(
                Some(line(&[(0.0, 0.0), (1.0, 0.0)])),
                &[("rid", "B".into())],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, _) = run(&input, json!({ "route_id_field": "rid" }));
        // Only the line-bearing route B is built.
        assert_eq!(out.outputs["route_count"], json!(1));
        assert_eq!(out.outputs["input_lines"], json!(1));
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = CreateRoutesTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input");
        assert!(
            bad(json!({ "input": "x.geojson" })).is_err(),
            "missing route_id_field"
        );
        assert!(bad(json!({ "input": "x.geojson", "route_id_field": "rid" })).is_ok());
        assert!(
            bad(
                json!({ "input": "x.geojson", "route_id_field": "rid", "measure_source": "BOGUS" })
            )
            .is_err(),
            "bad measure_source"
        );
        assert!(
            bad(json!({ "input": "x.geojson", "route_id_field": "rid", "measure_factor": "nan" }))
                .is_err(),
            "non-finite factor"
        );
        assert!(
            bad(json!({ "input": "x.geojson", "route_id_field": "rid", "snap_tolerance": "0.5" }))
                .is_ok(),
            "numeric string tolerance accepted"
        );
    }

    /// A named-but-missing measure field is a clear error.
    #[test]
    fn missing_measure_field_errors() {
        let mut layer = Layer::new("roads");
        layer.add_field(FieldDef::new("rid", FieldType::Text));
        layer
            .add_feature(
                Some(line(&[(0.0, 0.0), (1.0, 0.0)])),
                &[("rid", "A".into())],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input,
            "route_id_field": "rid",
            "measure_source": "ONE_FIELD",
            "from_measure_field": "nope",
        }))
        .unwrap();
        assert!(CreateRoutesTool.run(&args, &ctx()).is_err());
    }
}
