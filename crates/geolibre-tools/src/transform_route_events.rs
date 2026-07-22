//! GeoLibre tool: transform linear-referencing events from a source route
//! system onto a target route system.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Transform Route Events* (Linear
//! Referencing). The bundle already has overlay/dissolve/split machinery for
//! events on a single route system, but nothing re-references events measured
//! against one route system onto a *different* one — a routine maintenance task
//! when a linear network is re-versioned (new geometry, new route identifiers).
//! This is distinct from route recalibration, which keeps the same routes.
//!
//! The transform is a geometric composition of two classic LR operations:
//!
//!   1. **Event -> XY (source side).** Each point event carries a route
//!      identifier (`route_id_field`) and a measure (`measure_field`). The
//!      matching source route is looked up by identifier and walked to the
//!      event's measure to recover an XY location. Measures are modelled as
//!      cumulative planar distance from the route's first vertex (the ArcGIS
//!      *Create Routes* `LENGTH` convention); stored M values are not required.
//!
//!   2. **XY -> event (target side).** The recovered XY is relocated onto the
//!      nearest target route within `cluster_tolerance`. The perpendicular
//!      offset (`transfer_dist`) and the along-route distance of the snap point
//!      become the event's new target identifier and target measure.
//!
//! Events whose source route is missing, whose measure is invalid, or whose XY
//! has no target route within the cluster tolerance are dropped and counted;
//! the remainder are emitted as points located on the target system carrying
//! their original attributes plus the transformed reference.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Re-references linear events from a source route system onto a target system.
pub struct TransformRouteEventsTool;

impl Tool for TransformRouteEventsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "transform_route_events",
            display_name: "Transform Route Events",
            summary: "Re-reference point events measured against a source route system onto a target route system: recover each event's XY from its source route and measure, then relocate it onto the nearest target route within a cluster tolerance, emitting the target route id and measure — like ArcGIS Transform Route Events.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "source_events",
                    description: "Input point event layer. Each feature carries a route identifier (route_id_field) and a measure (measure_field).",
                    required: true,
                },
                ToolParamSpec {
                    name: "source_routes",
                    description: "Source route system: line layer whose features are keyed by route_id_field. Measures are distance from each route's first vertex.",
                    required: true,
                },
                ToolParamSpec {
                    name: "target_routes",
                    description: "Target route system: line layer onto which events are relocated. Route identifiers come from target_id_field.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output point vector path (driver from extension). If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "route_id_field",
                    description: "Field naming the route identifier in source_events and source_routes (default 'route_id'). Absent field -> the source route feature index is used.",
                    required: false,
                },
                ToolParamSpec {
                    name: "measure_field",
                    description: "Numeric field in source_events holding the source-system measure (default 'measure').",
                    required: false,
                },
                ToolParamSpec {
                    name: "target_id_field",
                    description: "Field naming the route identifier in target_routes (default: same as route_id_field). Absent field -> the target route feature index is used.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cluster_tolerance",
                    description: "Maximum XY distance (CRS units) from the recovered event location to a target route for the event to transfer. 0 or omitted means no limit (always snap to the nearest target route).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["source_events", "source_routes", "target_routes"] {
            if args.get(key).and_then(Value::as_str).is_none() {
                return Err(ToolError::Validation(format!(
                    "missing required string parameter '{key}'"
                )));
            }
        }
        let tol = parse_optional_f64(args, "cluster_tolerance")?.unwrap_or(0.0);
        if tol < 0.0 {
            return Err(ToolError::Validation(
                "parameter 'cluster_tolerance' must be non-negative".to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let source_events = require_str(args, "source_events")?;
        let source_routes = require_str(args, "source_routes")?;
        let target_routes = require_str(args, "target_routes")?;
        let output = parse_optional_str(args, "output")?;
        let route_id_field = parse_optional_str(args, "route_id_field")?.unwrap_or("route_id");
        let measure_field = parse_optional_str(args, "measure_field")?.unwrap_or("measure");
        let target_id_field =
            parse_optional_str(args, "target_id_field")?.unwrap_or(route_id_field);
        let cluster_tolerance = parse_optional_f64(args, "cluster_tolerance")?.unwrap_or(0.0);

        let events = load_input_layer(source_events)?;
        let src_routes = load_input_layer(source_routes)?;
        let tgt_routes = load_input_layer(target_routes)?;

        // Build the source route index (id -> polyline with cumulative measures).
        ctx.progress.info("indexing source routes");
        let src_index = build_route_index(&src_routes, route_id_field);
        if src_index.is_empty() {
            return Err(ToolError::Execution(
                "source_routes contains no usable line geometry".to_string(),
            ));
        }

        // Build the target route list (id + polyline with cumulative measures).
        ctx.progress.info("indexing target routes");
        let tgt_index: Vec<Route> = build_route_index(&tgt_routes, target_id_field)
            .into_values()
            .collect();
        if tgt_index.is_empty() {
            return Err(ToolError::Execution(
                "target_routes contains no usable line geometry".to_string(),
            ));
        }

        // Prepare output schema: original event fields + transform results.
        let mut out = Layer::new("transform_route_events").with_geom_type(GeometryType::Point);
        if let Some(epsg) = tgt_routes.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for fd in events.schema.fields() {
            out.add_field(fd.clone());
        }
        out.add_field(FieldDef::new("t_route_id", FieldType::Text));
        out.add_field(FieldDef::new("t_measure", FieldType::Float));
        out.add_field(FieldDef::new("transfer_dist", FieldType::Float));

        let id_idx = events.schema.field_index(route_id_field);
        let m_idx = events.schema.field_index(measure_field);

        let mut transferred = 0usize;
        let mut unmatched_source = 0usize;
        let mut invalid_measure = 0usize;
        let mut unmatched_target = 0usize;
        let mut max_transfer_dist = 0.0f64;

        let n = events.features.len().max(1);
        for (i, feat) in events.features.iter().enumerate() {
            // Resolve the event's source route identifier.
            let route_key = id_idx
                .and_then(|idx| feat.attributes.get(idx))
                .map(field_value_key)
                .unwrap_or_default();
            let Some(route) = src_index.get(&route_key) else {
                unmatched_source += 1;
                continue;
            };

            // Resolve the event measure.
            let Some(measure) = m_idx
                .and_then(|idx| feat.attributes.get(idx))
                .and_then(FieldValue::as_f64)
                .filter(|m| m.is_finite())
            else {
                invalid_measure += 1;
                continue;
            };

            // Recover the XY location at `measure` along the source route.
            let Some((x, y)) = route.interpolate(measure) else {
                invalid_measure += 1;
                continue;
            };

            // Relocate onto the nearest target route within the cluster tolerance.
            let Some(snap) = nearest_on_routes(&tgt_index, x, y, cluster_tolerance) else {
                unmatched_target += 1;
                continue;
            };

            let mut attrs = feat.attributes.clone();
            // Pad to the (possibly wider) output schema before appending results.
            attrs.resize(events.schema.len(), FieldValue::Null);
            attrs.push(FieldValue::Text(snap.route_id));
            attrs.push(FieldValue::Float(snap.measure));
            attrs.push(FieldValue::Float(snap.dist));
            out.push(wbvector::Feature {
                fid: transferred as u64,
                geometry: Some(Geometry::point(snap.x, snap.y)),
                attributes: attrs,
            });
            transferred += 1;
            max_transfer_dist = max_transfer_dist.max(snap.dist);

            ctx.progress.progress((i as f64 + 1.0) / n as f64);
        }

        ctx.progress
            .info(&format!("transferred {transferred} event(s)"));

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("transferred".to_string(), json!(transferred));
        outputs.insert("unmatched_source".to_string(), json!(unmatched_source));
        outputs.insert("invalid_measure".to_string(), json!(invalid_measure));
        outputs.insert("unmatched_target".to_string(), json!(unmatched_target));
        outputs.insert(
            "max_transfer_dist".to_string(),
            json!(if transferred > 0 {
                max_transfer_dist
            } else {
                0.0
            }),
        );
        Ok(ToolRunResult { outputs })
    }
}

// ── Route model ───────────────────────────────────────────────────────────────

/// A route: an ordered vertex list with the cumulative planar distance
/// (measure) at each vertex. `measure[0] == 0`, `measure.last() == length`.
struct Route {
    id: String,
    pts: Vec<(f64, f64)>,
    cum: Vec<f64>,
}

impl Route {
    fn length(&self) -> f64 {
        self.cum.last().copied().unwrap_or(0.0)
    }

    /// XY at along-route distance `m`, clamped to `[0, length]`.
    fn interpolate(&self, m: f64) -> Option<(f64, f64)> {
        if self.pts.len() < 2 {
            return None;
        }
        let total = self.length();
        if total <= 0.0 {
            return Some(self.pts[0]);
        }
        let m = m.clamp(0.0, total);
        let mut i = 0;
        while i + 1 < self.pts.len() && self.cum[i + 1] < m {
            i += 1;
        }
        let (ax, ay) = self.pts[i];
        let (bx, by) = self.pts[(i + 1).min(self.pts.len() - 1)];
        let seg = self.cum[(i + 1).min(self.cum.len() - 1)] - self.cum[i];
        if seg <= 0.0 {
            return Some((ax, ay));
        }
        let t = ((m - self.cum[i]) / seg).clamp(0.0, 1.0);
        Some((ax + (bx - ax) * t, ay + (by - ay) * t))
    }
}

/// Result of relocating an XY onto a target route.
struct Snap {
    route_id: String,
    measure: f64,
    x: f64,
    y: f64,
    dist: f64,
}

/// Builds an id -> [`Route`] map. Line features sharing an identifier are
/// concatenated in input order; multi-part lines are concatenated part by part.
/// When the id field is absent, the route feature index is used as the id.
fn build_route_index(layer: &Layer, id_field: &str) -> BTreeMap<String, Route> {
    let id_idx = layer.schema.field_index(id_field);
    // Preserve input order of parts per id.
    let mut parts: Vec<(String, Vec<(f64, f64)>)> = Vec::new();
    let mut order: BTreeMap<String, usize> = BTreeMap::new();

    for (fi, feat) in layer.features.iter().enumerate() {
        let key = id_idx
            .and_then(|idx| feat.attributes.get(idx))
            .map(field_value_key)
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| fi.to_string());
        let Some(geom) = feat.geometry.as_ref() else {
            continue;
        };
        for line in geom_lines(geom) {
            if line.len() < 2 {
                continue;
            }
            if let Some(&pos) = order.get(&key) {
                parts[pos].1.extend(line);
            } else {
                order.insert(key.clone(), parts.len());
                parts.push((key.clone(), line));
            }
        }
    }

    let mut index = BTreeMap::new();
    for (id, pts) in parts {
        // Drop consecutive duplicate vertices so zero-length segments don't
        // create ambiguous measures.
        let mut clean: Vec<(f64, f64)> = Vec::with_capacity(pts.len());
        for p in pts {
            if clean.last().is_none_or(|&q| q != p) {
                clean.push(p);
            }
        }
        if clean.len() < 2 {
            continue;
        }
        let mut cum = vec![0.0f64; clean.len()];
        for i in 1..clean.len() {
            cum[i] = cum[i - 1] + dist(clean[i - 1], clean[i]);
        }
        index.insert(
            id.clone(),
            Route {
                id,
                pts: clean,
                cum,
            },
        );
    }
    index
}

/// Finds the nearest point across all target routes to `(x, y)`. When
/// `tolerance > 0`, candidates farther than the tolerance are rejected.
fn nearest_on_routes(routes: &[Route], x: f64, y: f64, tolerance: f64) -> Option<Snap> {
    let mut best: Option<Snap> = None;
    for route in routes {
        let Some((px, py, measure, d)) = nearest_on_route(route, x, y) else {
            continue;
        };
        if tolerance > 0.0 && d > tolerance {
            continue;
        }
        if best.as_ref().is_none_or(|b| d < b.dist) {
            best = Some(Snap {
                route_id: route.id.clone(),
                measure,
                x: px,
                y: py,
                dist: d,
            });
        }
    }
    best
}

/// Nearest point on a single route to `(x, y)`. Returns
/// `(snap_x, snap_y, along_route_measure, distance)`.
fn nearest_on_route(route: &Route, x: f64, y: f64) -> Option<(f64, f64, f64, f64)> {
    if route.pts.len() < 2 {
        let (px, py) = *route.pts.first()?;
        return Some((px, py, 0.0, dist((px, py), (x, y))));
    }
    let mut best: Option<(f64, f64, f64, f64)> = None;
    for i in 0..route.pts.len() - 1 {
        let (ax, ay) = route.pts[i];
        let (bx, by) = route.pts[i + 1];
        let (dx, dy) = (bx - ax, by - ay);
        let seg2 = dx * dx + dy * dy;
        let t = if seg2 <= 0.0 {
            0.0
        } else {
            (((x - ax) * dx + (y - ay) * dy) / seg2).clamp(0.0, 1.0)
        };
        let (px, py) = (ax + dx * t, ay + dy * t);
        let d = dist((px, py), (x, y));
        let measure = route.cum[i] + t * seg2.sqrt();
        if best.is_none_or(|b| d < b.3) {
            best = Some((px, py, measure, d));
        }
    }
    best
}

// ── Geometry helpers ──────────────────────────────────────────────────────────

/// Extracts every polyline part from a geometry as a vertex list.
fn geom_lines(geom: &Geometry) -> Vec<Vec<(f64, f64)>> {
    match geom {
        Geometry::LineString(cs) => vec![cs.iter().map(|c| (c.x, c.y)).collect()],
        Geometry::MultiLineString(parts) => parts
            .iter()
            .map(|cs| cs.iter().map(|c| (c.x, c.y)).collect())
            .collect(),
        _ => Vec::new(),
    }
}

fn dist(a: (f64, f64), b: (f64, f64)) -> f64 {
    let dx = a.0 - b.0;
    let dy = a.1 - b.1;
    (dx * dx + dy * dy).sqrt()
}

/// Stable string key for a route identifier field value. Integers and integral
/// floats collapse to the same key so `5` and `5.0` match.
fn field_value_key(v: &FieldValue) -> String {
    match v {
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) if f.fract() == 0.0 && f.is_finite() => format!("{}", *f as i64),
        FieldValue::Float(f) => format!("{f}"),
        FieldValue::Text(s) | FieldValue::Date(s) | FieldValue::DateTime(s) => s.clone(),
        FieldValue::Boolean(b) => b.to_string(),
        FieldValue::Null => String::new(),
        FieldValue::Blob(_) => String::new(),
    }
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

/// Parses an optional numeric parameter, accepting a JSON number or a numeric
/// string (host UIs often post form values as strings).
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

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;
    use wbcore::{AllowAllCapabilities, ProgressSink, ToolContext};
    use wbvector::{memory_store, Coord, FieldDef, FieldType, Geometry, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn to_args(v: serde_json::Value) -> ToolArgs {
        serde_json::from_value(v).expect("args object")
    }

    /// A straight route along +x from (x0,0) to (x0+len,0) with the given id.
    fn straight_route(id: &str, x0: f64, len: f64) -> Layer {
        let mut l = Layer::new("routes").with_geom_type(GeometryType::LineString);
        l.add_field(FieldDef::new("route_id", FieldType::Text));
        l.add_feature(
            Some(Geometry::line_string(vec![
                Coord::xy(x0, 0.0),
                Coord::xy(x0 + len, 0.0),
            ])),
            &[("route_id", id.into())],
        )
        .unwrap();
        l
    }

    fn events(pairs: &[(&str, f64)]) -> Layer {
        let mut l = Layer::new("events").with_geom_type(GeometryType::Point);
        l.add_field(FieldDef::new("route_id", FieldType::Text));
        l.add_field(FieldDef::new("measure", FieldType::Float));
        l.add_field(FieldDef::new("label", FieldType::Text));
        for (i, (rid, m)) in pairs.iter().enumerate() {
            l.add_feature(
                Some(Geometry::point(0.0, 0.0)),
                &[
                    ("route_id", (*rid).into()),
                    ("measure", (*m).into()),
                    ("label", format!("e{i}").into()),
                ],
            )
            .unwrap();
        }
        l
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        TransformRouteEventsTool
            .run(&to_args(args), &ctx())
            .expect("run should succeed")
    }

    /// Identity transform: source == target routes. Measures and locations must
    /// be preserved (transfer_dist ~ 0, t_measure ~ source measure).
    #[test]
    fn identity_transform_preserves_measures() {
        let ev = memory_store::make_vector_memory_path(&memory_store::put_vector(events(&[
            ("A", 3.0),
            ("A", 7.0),
        ])));
        let routes = memory_store::make_vector_memory_path(&memory_store::put_vector(
            straight_route("A", 0.0, 10.0),
        ));
        let src = memory_store::make_vector_memory_path(&memory_store::put_vector(straight_route(
            "A", 0.0, 10.0,
        )));

        let res = run(json!({
            "source_events": ev,
            "source_routes": routes,
            "target_routes": src,
        }));
        assert_eq!(res.outputs["transferred"], json!(2));

        let out_path = res.outputs["output"].as_str().unwrap();
        let layer = load_input_layer(out_path).unwrap();
        assert_eq!(layer.len(), 2);
        let tm = layer.schema.field_index("t_measure").unwrap();
        let td = layer.schema.field_index("transfer_dist").unwrap();
        let m0 = layer[0].attributes[tm].as_f64().unwrap();
        assert!((m0 - 3.0).abs() < 1e-9, "measure preserved, got {m0}");
        let d0 = layer[0].attributes[td].as_f64().unwrap();
        assert!(d0 < 1e-9, "no offset for identity, got {d0}");
        // Event geometry recovered at (3, 0).
        if let Some(Geometry::Point(c)) = &layer[0].geometry {
            assert!((c.x - 3.0).abs() < 1e-9 && c.y.abs() < 1e-9);
        } else {
            panic!("expected point geometry");
        }
    }

    /// A parallel target route shifted in +y: measures preserved along the
    /// same-shape route, offset equals the shift.
    #[test]
    fn parallel_shift_reports_offset() {
        let ev =
            memory_store::make_vector_memory_path(&memory_store::put_vector(events(&[("A", 4.0)])));
        let src = memory_store::make_vector_memory_path(&memory_store::put_vector(straight_route(
            "A", 0.0, 10.0,
        )));
        // Target route parallel to source but 2.0 above it.
        let mut tgt_layer = Layer::new("routes").with_geom_type(GeometryType::LineString);
        tgt_layer.add_field(FieldDef::new("route_id", FieldType::Text));
        tgt_layer
            .add_feature(
                Some(Geometry::line_string(vec![
                    Coord::xy(0.0, 2.0),
                    Coord::xy(10.0, 2.0),
                ])),
                &[("route_id", "T1".into())],
            )
            .unwrap();
        let tgt = memory_store::make_vector_memory_path(&memory_store::put_vector(tgt_layer));

        let res = run(json!({
            "source_events": ev,
            "source_routes": src,
            "target_routes": tgt,
        }));
        assert_eq!(res.outputs["transferred"], json!(1));
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        let tm = layer.schema.field_index("t_measure").unwrap();
        let td = layer.schema.field_index("transfer_dist").unwrap();
        let ti = layer.schema.field_index("t_route_id").unwrap();
        assert!((layer[0].attributes[tm].as_f64().unwrap() - 4.0).abs() < 1e-9);
        assert!((layer[0].attributes[td].as_f64().unwrap() - 2.0).abs() < 1e-9);
        assert_eq!(layer[0].attributes[ti].as_str().unwrap(), "T1");
    }

    /// Events whose source route id has no match are dropped and counted.
    #[test]
    fn non_matching_source_is_passed_over() {
        let ev = memory_store::make_vector_memory_path(&memory_store::put_vector(events(&[
            ("A", 5.0),
            ("GHOST", 5.0),
        ])));
        let src = memory_store::make_vector_memory_path(&memory_store::put_vector(straight_route(
            "A", 0.0, 10.0,
        )));
        let tgt = memory_store::make_vector_memory_path(&memory_store::put_vector(straight_route(
            "A", 0.0, 10.0,
        )));
        let res = run(json!({
            "source_events": ev,
            "source_routes": src,
            "target_routes": tgt,
        }));
        assert_eq!(res.outputs["transferred"], json!(1));
        assert_eq!(res.outputs["unmatched_source"], json!(1));
    }

    /// A target route beyond the cluster tolerance rejects the event.
    #[test]
    fn cluster_tolerance_rejects_far_target() {
        let ev =
            memory_store::make_vector_memory_path(&memory_store::put_vector(events(&[("A", 5.0)])));
        let src = memory_store::make_vector_memory_path(&memory_store::put_vector(straight_route(
            "A", 0.0, 10.0,
        )));
        // Target 100 units away in y.
        let tgt = memory_store::make_vector_memory_path(&memory_store::put_vector({
            let mut l = Layer::new("r").with_geom_type(GeometryType::LineString);
            l.add_field(FieldDef::new("route_id", FieldType::Text));
            l.add_feature(
                Some(Geometry::line_string(vec![
                    Coord::xy(0.0, 100.0),
                    Coord::xy(10.0, 100.0),
                ])),
                &[("route_id", "T".into())],
            )
            .unwrap();
            l
        }));
        let res = run(json!({
            "source_events": ev,
            "source_routes": src,
            "target_routes": tgt,
            "cluster_tolerance": 5.0,
        }));
        assert_eq!(res.outputs["transferred"], json!(0));
        assert_eq!(res.outputs["unmatched_target"], json!(1));
    }

    #[test]
    fn rejects_bad_parameters() {
        // Missing required inputs.
        assert!(TransformRouteEventsTool
            .validate(&to_args(json!({ "source_events": "a" })))
            .is_err());
        // Negative tolerance.
        assert!(TransformRouteEventsTool
            .validate(&to_args(json!({
                "source_events": "a",
                "source_routes": "b",
                "target_routes": "c",
                "cluster_tolerance": -1.0
            })))
            .is_err());
        // Valid.
        assert!(TransformRouteEventsTool
            .validate(&to_args(json!({
                "source_events": "a",
                "source_routes": "b",
                "target_routes": "c"
            })))
            .is_ok());
    }
}
