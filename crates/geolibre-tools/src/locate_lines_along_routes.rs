//! GeoLibre tool: locate line/polygon features along routes as from/to-measure
//! line events.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Locate Features Along Routes*
//! (Linear Referencing) for the **line/polygon-event** case. The bundled
//! `locate_points_along_routes` handles points only (single-measure point
//! events); this tool overlays line or polygon features on measured routes and
//! emits from-measure/to-measure line events (`RID`, `FMEAS`, `TMEAS`).
//!
//! A route is a polyline whose measure at any point is the cumulative
//! arc-length from its start. For each input feature the tool walks the
//! feature's geometry (line vertices, or polygon boundary rings), densifies it
//! so no gap exceeds the search `tolerance`, and projects every sample point
//! onto each route. Samples whose perpendicular distance to a route is within
//! `tolerance` contribute their measure; the located interval for that
//! feature/route pair is `[min measure, max measure]`. The event geometry is
//! the sub-portion of the route cut between those two measures, so its length
//! equals `TMEAS - FMEAS` exactly (a measure-conservation invariant).
//!
//! Scope for v1: one contiguous interval per (feature, route) pair. A single
//! feature that touches a route in several disjoint stretches yields the
//! measure span covering them, not one event per stretch.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Overlays line/polygon features on measured routes and emits from/to-measure
/// line events.
pub struct LocateLinesAlongRoutesTool;

impl Tool for LocateLinesAlongRoutesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "locate_lines_along_routes",
            display_name: "Locate Lines Along Routes",
            summary: "Overlay line or polygon features on measured routes to produce from/to-measure line events (RID, FMEAS, TMEAS), like ArcGIS Locate Features Along Routes (line/polygon events).",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input_features",
                    description: "Input line or polygon vector layer to locate along the routes.",
                    required: true,
                },
                ToolParamSpec {
                    name: "routes",
                    description: "Route (line) vector layer providing the linear-referencing measures.",
                    required: true,
                },
                ToolParamSpec {
                    name: "route_id_field",
                    description: "Field on the routes layer holding each route's identifier (RID). If omitted, the route's 0-based feature index is used.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Search radius in CRS units. A feature vertex within this distance of a route is treated as coincident. Required (> 0).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output line-event vector path (driver from extension). If omitted, the result is stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input_features")?;
        require_str(args, "routes")?;
        parse_tolerance(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input_features")?;
        let routes_path = require_str(args, "routes")?;
        let route_id_field = parse_optional_str(args, "route_id_field")?;
        let tolerance = parse_tolerance(args)?;
        let output = parse_optional_str(args, "output")?;

        let features = load_input_layer(input)?;
        let routes_layer = load_input_layer(routes_path)?;

        // Build measured routes.
        let rid_type = match route_id_field {
            Some(name) => {
                let idx = routes_layer.schema.field_index(name).ok_or_else(|| {
                    ToolError::Validation(format!(
                        "route_id_field '{name}' not found on the routes layer"
                    ))
                })?;
                routes_layer.schema.fields()[idx].field_type
            }
            None => FieldType::Integer,
        };
        let rid_field_idx = route_id_field.and_then(|name| routes_layer.schema.field_index(name));

        let routes: Vec<Route> = routes_layer
            .features
            .iter()
            .enumerate()
            .filter_map(|(i, f)| {
                let pts = concat_line(f.geometry.as_ref()?);
                if pts.len() < 2 {
                    return None;
                }
                let rid = match rid_field_idx {
                    Some(idx) => f.attributes.get(idx).cloned().unwrap_or(FieldValue::Null),
                    None => FieldValue::Integer(i as i64),
                };
                Some(Route::new(pts, rid))
            })
            .collect();

        if routes.is_empty() {
            return Err(ToolError::Execution(
                "routes layer contains no usable polylines".to_string(),
            ));
        }

        // Prepare output layer.
        let mut out = Layer::new("route_events").with_geom_type(GeometryType::LineString);
        if let Some(epsg) = routes_layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("FID", FieldType::Integer));
        out.add_field(FieldDef::new("RID", rid_type));
        out.add_field(FieldDef::new("FMEAS", FieldType::Float));
        out.add_field(FieldDef::new("TMEAS", FieldType::Float));

        let mut event_count = 0usize;
        for (fidx, feature) in features.features.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let samples = sample_feature(geom, tolerance);
            if samples.is_empty() {
                continue;
            }
            for route in &routes {
                // Collect measures of samples within tolerance of this route.
                let mut m_min = f64::INFINITY;
                let mut m_max = f64::NEG_INFINITY;
                for s in &samples {
                    let (dist, meas) = route.project(s);
                    if dist <= tolerance {
                        m_min = m_min.min(meas);
                        m_max = m_max.max(meas);
                    }
                }
                if m_min.is_finite() && m_max.is_finite() {
                    let sub = route.cut(m_min, m_max);
                    if sub.len() < 2 {
                        continue;
                    }
                    out.add_feature(
                        Some(Geometry::line_string(sub)),
                        &[
                            ("FID", FieldValue::Integer(fidx as i64)),
                            ("RID", route.rid.clone()),
                            ("FMEAS", FieldValue::Float(m_min)),
                            ("TMEAS", FieldValue::Float(m_max)),
                        ],
                    )
                    .map_err(|e| {
                        ToolError::Execution(format!("failed writing event feature: {e}"))
                    })?;
                    event_count += 1;
                }
            }
        }

        ctx.progress
            .info(&format!("located {event_count} line event(s)"));

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("event_count".to_string(), json!(event_count));
        outputs.insert("route_count".to_string(), json!(routes.len()));
        Ok(ToolRunResult { outputs })
    }
}

// ── Route (measured polyline) ─────────────────────────────────────────────────

struct Route {
    pts: Vec<Coord>,
    /// Cumulative measure at each vertex (measure = arc-length from start).
    cum: Vec<f64>,
    rid: FieldValue,
}

impl Route {
    fn new(pts: Vec<Coord>, rid: FieldValue) -> Self {
        let mut cum = vec![0.0f64; pts.len()];
        for i in 1..pts.len() {
            cum[i] = cum[i - 1] + dist(&pts[i - 1], &pts[i]);
        }
        Route { pts, cum, rid }
    }

    fn total(&self) -> f64 {
        *self.cum.last().unwrap_or(&0.0)
    }

    /// Nearest point on the route to `q`; returns (perpendicular distance,
    /// measure at that nearest point).
    fn project(&self, q: &Coord) -> (f64, f64) {
        let mut best_d = f64::INFINITY;
        let mut best_m = 0.0;
        for i in 0..self.pts.len() - 1 {
            let (a, b) = (&self.pts[i], &self.pts[i + 1]);
            let seg = dist(a, b);
            let (px, py, t) = if seg <= 0.0 {
                (a.x, a.y, 0.0)
            } else {
                let t = (((q.x - a.x) * (b.x - a.x) + (q.y - a.y) * (b.y - a.y)) / (seg * seg))
                    .clamp(0.0, 1.0);
                (a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t, t)
            };
            let d = (q.x - px).hypot(q.y - py);
            if d < best_d {
                best_d = d;
                best_m = self.cum[i] + t * seg;
            }
        }
        (best_d, best_m)
    }

    /// Coordinates of the route sub-portion between measures `m0` and `m1`.
    fn cut(&self, m0: f64, m1: f64) -> Vec<Coord> {
        let total = self.total();
        let lo = m0.clamp(0.0, total);
        let hi = m1.clamp(0.0, total);
        if hi <= lo {
            return Vec::new();
        }
        let mut out = vec![self.point_at(lo)];
        for i in 0..self.pts.len() {
            if self.cum[i] > lo && self.cum[i] < hi {
                out.push(self.pts[i].clone());
            }
        }
        out.push(self.point_at(hi));
        // Drop consecutive duplicates.
        let mut dedup: Vec<Coord> = Vec::with_capacity(out.len());
        for c in out {
            if dedup.last().is_none_or(|l| dist(l, &c) > 1e-12) {
                dedup.push(c);
            }
        }
        dedup
    }

    /// Interpolated point at measure `m` along the route.
    fn point_at(&self, m: f64) -> Coord {
        let total = self.total();
        let m = m.clamp(0.0, total);
        let mut i = 0;
        while i + 1 < self.pts.len() && self.cum[i + 1] < m {
            i += 1;
        }
        let (a, b) = (&self.pts[i], &self.pts[(i + 1).min(self.pts.len() - 1)]);
        let seg = dist(a, b);
        if seg <= 0.0 {
            return a.clone();
        }
        let t = ((m - self.cum[i]) / seg).clamp(0.0, 1.0);
        Coord::xy(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t)
    }
}

// ── Geometry helpers ──────────────────────────────────────────────────────────

fn dist(a: &Coord, b: &Coord) -> f64 {
    (a.x - b.x).hypot(a.y - b.y)
}

/// Flattens a line geometry into a single ordered vertex chain, concatenating
/// multi-part lines and dropping duplicate joints. Non-line geometry yields an
/// empty chain.
fn concat_line(geom: &Geometry) -> Vec<Coord> {
    let mut out: Vec<Coord> = Vec::new();
    let mut push = |c: &Coord| {
        if out.last().is_none_or(|l| dist(l, c) > 1e-12) {
            out.push(c.clone());
        }
    };
    match geom {
        Geometry::LineString(cs) => cs.iter().for_each(&mut push),
        Geometry::MultiLineString(lines) => {
            for l in lines {
                l.iter().for_each(&mut push);
            }
        }
        _ => {}
    }
    out
}

/// Sample points along an input feature's geometry, densified so no gap between
/// consecutive samples exceeds `step` (the search tolerance). For polygons the
/// boundary rings (exterior + holes) are walked. Points and multipoints yield no
/// samples (the point-event case is handled by `locate_points_along_routes`).
fn sample_feature(geom: &Geometry, step: f64) -> Vec<Coord> {
    let mut out: Vec<Coord> = Vec::new();
    let mut chains: Vec<Vec<Coord>> = Vec::new();
    match geom {
        Geometry::LineString(cs) => chains.push(cs.clone()),
        Geometry::MultiLineString(lines) => {
            for l in lines {
                chains.push(l.clone());
            }
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            chains.push(exterior.0.clone());
            for r in interiors {
                chains.push(r.0.clone());
            }
        }
        Geometry::MultiPolygon(polys) => {
            for (ext, ints) in polys {
                chains.push(ext.0.clone());
                for r in ints {
                    chains.push(r.0.clone());
                }
            }
        }
        _ => {}
    }
    let step = if step.is_finite() && step > 0.0 {
        step
    } else {
        f64::INFINITY
    };
    for chain in &chains {
        for w in chain.windows(2) {
            let (a, b) = (&w[0], &w[1]);
            out.push(a.clone());
            let seg = dist(a, b);
            if seg > step {
                let n = (seg / step).ceil() as usize;
                for k in 1..n {
                    let t = k as f64 / n as f64;
                    out.push(Coord::xy(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t));
                }
            }
        }
        if let Some(last) = chain.last() {
            out.push(last.clone());
        }
    }
    out
}

// ── Parameter parsing ─────────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn parse_tolerance(args: &ToolArgs) -> Result<f64, ToolError> {
    let v = match args.get("tolerance") {
        None | Some(Value::Null) => {
            return Err(ToolError::Validation(
                "required parameter 'tolerance' is missing".to_string(),
            ))
        }
        Some(Value::Number(n)) => n.as_f64(),
        Some(Value::String(s)) if s.trim().is_empty() => {
            return Err(ToolError::Validation(
                "required parameter 'tolerance' is missing".to_string(),
            ))
        }
        Some(Value::String(s)) => s.trim().parse::<f64>().ok(),
        Some(_) => {
            return Err(ToolError::Validation(
                "parameter 'tolerance' must be a number".to_string(),
            ))
        }
    };
    match v {
        Some(t) if t.is_finite() && t > 0.0 => Ok(t),
        _ => Err(ToolError::Validation(
            "'tolerance' must be a positive number".to_string(),
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::memory_store;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Single route from (0,0)-(100,0) with an integer RID field.
    fn route_layer() -> String {
        let mut l = Layer::new("routes")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("rid", FieldType::Integer));
        l.add_feature(
            Some(Geometry::line_string(vec![
                Coord::xy(0.0, 0.0),
                Coord::xy(100.0, 0.0),
            ])),
            &[("rid", FieldValue::Integer(7))],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn line_features(lines: &[Vec<(f64, f64)>]) -> String {
        let mut l = Layer::new("feats")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        for coords in lines {
            let cs = coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
            l.add_feature(Some(Geometry::line_string(cs)), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = LocateLinesAlongRoutesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn field(layer: &Layer, feat: usize, name: &str) -> FieldValue {
        let idx = layer.schema.field_index(name).unwrap();
        layer.features[feat].attributes[idx].clone()
    }

    fn line_len(layer: &Layer, feat: usize) -> f64 {
        match layer.features[feat].geometry.as_ref().unwrap() {
            Geometry::LineString(cs) => cs
                .windows(2)
                .map(|w| ((w[0].x - w[1].x).powi(2) + (w[0].y - w[1].y).powi(2)).sqrt())
                .sum(),
            other => panic!("expected line, got {other:?}"),
        }
    }

    /// A line lying just beside the route from x=20 to x=70 yields one event with
    /// FMEAS=20, TMEAS=70 and geometry length == TMEAS-FMEAS (measure conserved).
    #[test]
    fn locates_line_measures() {
        let routes = route_layer();
        let feats = line_features(&[vec![(20.0, 0.5), (70.0, 0.5)]]);
        let (out, layer) = run(json!({
            "input_features": feats, "routes": routes,
            "route_id_field": "rid", "tolerance": 2.0,
        }));
        assert_eq!(out.outputs["event_count"], json!(1));
        assert_eq!(field(&layer, 0, "RID"), FieldValue::Integer(7));
        let fm = field(&layer, 0, "FMEAS").as_f64().unwrap();
        let tm = field(&layer, 0, "TMEAS").as_f64().unwrap();
        assert!((fm - 20.0).abs() < 1e-6, "FMEAS {fm}");
        assert!((tm - 70.0).abs() < 1e-6, "TMEAS {tm}");
        assert!(
            (line_len(&layer, 0) - (tm - fm)).abs() < 1e-6,
            "measure not conserved by geometry length"
        );
    }

    /// A polygon straddling the route produces an event spanning its footprint.
    #[test]
    fn locates_polygon_footprint() {
        let routes = route_layer();
        let mut l = Layer::new("poly")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        // Box from x[30,60], y[-1,1] centred on the route.
        l.add_feature(
            Some(Geometry::polygon(
                vec![
                    Coord::xy(30.0, -1.0),
                    Coord::xy(60.0, -1.0),
                    Coord::xy(60.0, 1.0),
                    Coord::xy(30.0, 1.0),
                    Coord::xy(30.0, -1.0),
                ],
                vec![],
            )),
            &[],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        let feats = memory_store::make_vector_memory_path(&id);
        let (out, layer) = run(json!({
            "input_features": feats, "routes": routes,
            "route_id_field": "rid", "tolerance": 0.5,
        }));
        assert_eq!(out.outputs["event_count"], json!(1));
        let fm = field(&layer, 0, "FMEAS").as_f64().unwrap();
        let tm = field(&layer, 0, "TMEAS").as_f64().unwrap();
        assert!((fm - 30.0).abs() < 1e-6, "FMEAS {fm}");
        assert!((tm - 60.0).abs() < 1e-6, "TMEAS {tm}");
    }

    /// A feature far from every route contributes no event (pass-through).
    #[test]
    fn distant_feature_no_event() {
        let routes = route_layer();
        let feats = line_features(&[
            vec![(20.0, 0.5), (70.0, 0.5)],     // near -> 1 event
            vec![(0.0, 500.0), (100.0, 500.0)], // far -> 0 events
        ]);
        let (out, layer) = run(json!({
            "input_features": feats, "routes": routes, "tolerance": 2.0,
        }));
        assert_eq!(out.outputs["event_count"], json!(1));
        // With no route_id_field, RID falls back to the route's feature index.
        assert_eq!(field(&layer, 0, "RID"), FieldValue::Integer(0));
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            LocateLinesAlongRoutesTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input_features": "a.geojson" })).is_err()); // no routes
        assert!(bad(json!({ "input_features": "a.geojson", "routes": "r.geojson" })).is_err()); // no tolerance
        assert!(bad(
            json!({ "input_features": "a.geojson", "routes": "r.geojson", "tolerance": 0 })
        )
        .is_err()); // tolerance must be > 0
        assert!(bad(
            json!({ "input_features": "a.geojson", "routes": "r.geojson", "tolerance": 5 })
        )
        .is_ok());
    }
}
