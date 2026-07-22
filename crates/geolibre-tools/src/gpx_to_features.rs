//! GeoLibre tool: convert a GPX file into GIS features.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *GPX To Features* (Conversion): parse a
//! `.gpx` file into waypoint **points** (from `<wpt>`) and/or track/route
//! **polylines** (from each `<trkseg>`'s `<trkpt>` vertices and each `<rte>`'s
//! `<rtept>` vertices), carrying the `ele` (elevation), `time`, and `name`
//! attributes GPS receivers record.
//!
//! GPX is the single most common consumer / field-GPS exchange format, and
//! nothing in the repo or the bundled whitebox-wasm suite reads it. The work is
//! pure XML parsing plus geometry assembly — squarely in-stack — and the
//! GeoJSON/PMTiles outputs drop straight into the web map. The XML is read with
//! the pure-Rust [`quick-xml`](https://crates.io/crates/quick-xml) pull parser,
//! so it stays inside the GDAL/GEOS/PROJ-free, WASM-first stack.
//!
//! GPX coordinates are always WGS84, so output is EPSG:4326.

use std::collections::BTreeMap;

use quick_xml::events::{BytesStart, Event};
use quick_xml::Reader;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{parse_optional_str, write_or_store_layer};

pub struct GpxToFeaturesTool;

impl Tool for GpxToFeaturesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "gpx_to_features",
            display_name: "GPX To Features",
            summary: "Parse a GPX file into waypoint points and/or track/route polylines, carrying ele/time/name.",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Path to a .gpx file (GPS Exchange Format XML).",
                    required: true,
                },
                ToolParamSpec {
                    name: "points_output",
                    description: "Optional output point vector path for waypoints (driver from its extension). If neither output is given, waypoints are stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "lines_output",
                    description: "Optional output line vector path for tracks and routes (driver from its extension).",
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
        let path = args
            .get("input")
            .and_then(Value::as_str)
            .filter(|s| !s.trim().is_empty())
            .ok_or_else(|| {
                ToolError::Validation("missing required parameter 'input'".to_string())
            })?;
        let points_output = parse_optional_str(args, "points_output")?;
        let lines_output = parse_optional_str(args, "lines_output")?;

        let text = std::fs::read_to_string(path)
            .map_err(|e| ToolError::Execution(format!("failed reading {path}: {e}")))?;
        let parsed = parse_gpx(&text)?;

        // ── Waypoint points layer ───────────────────────────────────────────
        let mut points = Layer::new("gpx_waypoints")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        points.add_field(FieldDef::new("name", FieldType::Text));
        points.add_field(FieldDef::new("ele", FieldType::Float));
        points.add_field(FieldDef::new("time", FieldType::Text));
        for wpt in &parsed.waypoints {
            points
                .add_feature(
                    Some(Geometry::point(wpt.lon, wpt.lat)),
                    &[
                        (
                            "name",
                            FieldValue::Text(wpt.name.clone().unwrap_or_default()),
                        ),
                        (
                            "ele",
                            wpt.ele.map(FieldValue::Float).unwrap_or(FieldValue::Null),
                        ),
                        (
                            "time",
                            FieldValue::Text(wpt.time.clone().unwrap_or_default()),
                        ),
                    ],
                )
                .map_err(|e| ToolError::Execution(format!("failed adding waypoint: {e}")))?;
        }
        let waypoint_count = points.len();

        // ── Track / route line layer ────────────────────────────────────────
        let mut lines = Layer::new("gpx_lines")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(4326);
        lines.add_field(FieldDef::new("name", FieldType::Text));
        lines.add_field(FieldDef::new("type", FieldType::Text));
        lines.add_field(FieldDef::new("n_points", FieldType::Integer));
        let mut track_vertex_count = 0usize;
        for line in &parsed.lines {
            if line.coords.len() < 2 {
                continue;
            }
            if line.kind == LineKind::Track {
                track_vertex_count += line.coords.len();
            }
            let coords: Vec<Coord> = line.coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
            let n = coords.len() as i64;
            lines
                .add_feature(
                    Some(Geometry::LineString(coords)),
                    &[
                        (
                            "name",
                            FieldValue::Text(line.name.clone().unwrap_or_default()),
                        ),
                        ("type", FieldValue::Text(line.kind.label().to_string())),
                        ("n_points", FieldValue::Integer(n)),
                    ],
                )
                .map_err(|e| ToolError::Execution(format!("failed adding line: {e}")))?;
        }
        let line_count = lines.len();
        ctx.progress.info(&format!(
            "parsed {waypoint_count} waypoint(s), {line_count} track/route line(s)"
        ));

        // Write lines only when a path is given; write points unless the caller
        // asked only for lines. With no output path a layer is stored in memory
        // and its handle returned.
        let mut lines_path: Option<String> = None;
        if lines_output.is_some() {
            lines_path = Some(write_or_store_layer(lines, lines_output)?);
        }
        let write_points = points_output.is_some() || lines_output.is_none();
        let points_path = if write_points {
            Some(write_or_store_layer(points, points_output)?)
        } else {
            None
        };

        let mut outputs = BTreeMap::new();
        if let Some(p) = &points_path {
            outputs.insert("output".to_string(), json!(p));
            outputs.insert("points_output".to_string(), json!(p));
        }
        if let Some(p) = &lines_path {
            outputs.insert("lines_output".to_string(), json!(p));
        }
        outputs.insert("waypoint_count".to_string(), json!(waypoint_count));
        outputs.insert("point_count".to_string(), json!(waypoint_count));
        outputs.insert("line_count".to_string(), json!(line_count));
        outputs.insert("track_vertex_count".to_string(), json!(track_vertex_count));
        Ok(ToolRunResult { outputs })
    }
}

/// A parsed GPX waypoint (`<wpt>`) with its optional recorded attributes.
struct Waypoint {
    lat: f64,
    lon: f64,
    ele: Option<f64>,
    time: Option<String>,
    name: Option<String>,
}

/// A point element (`<wpt>`/`<trkpt>`/`<rtept>`) under construction while its
/// `lat`/`lon` attributes and `<ele>`/`<time>`/`<name>` children are collected.
#[derive(Default)]
struct PendingPoint {
    lat: Option<f64>,
    lon: Option<f64>,
    ele: Option<f64>,
    time: Option<String>,
    name: Option<String>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum LineKind {
    Track,
    Route,
}

impl LineKind {
    fn label(self) -> &'static str {
        match self {
            LineKind::Track => "track",
            LineKind::Route => "route",
        }
    }
}

/// A parsed polyline: one per `<trkseg>` (track) or per `<rte>` (route).
struct GpxLine {
    kind: LineKind,
    name: Option<String>,
    coords: Vec<(f64, f64)>, // (lon, lat)
}

/// The waypoints and lines extracted from a GPX document.
struct ParsedGpx {
    waypoints: Vec<Waypoint>,
    lines: Vec<GpxLine>,
}

/// Reads `lat`/`lon` attributes off a point element (`wpt`/`trkpt`/`rtept`).
fn read_lat_lon(e: &BytesStart) -> (Option<f64>, Option<f64>) {
    let (mut lat, mut lon) = (None, None);
    for attr in e.attributes().flatten() {
        let key = attr.key.local_name();
        // lat/lon are plain decimal numbers, so no XML unescaping is needed.
        let val = std::str::from_utf8(attr.value.as_ref())
            .ok()
            .and_then(|v| v.trim().parse::<f64>().ok());
        match key.as_ref() {
            b"lat" => lat = val,
            b"lon" => lon = val,
            _ => {}
        }
    }
    (lat, lon)
}

/// Streaming pull-parse of a GPX document into waypoints and track/route lines.
///
/// Tracks a small amount of context so the ambiguous `<name>` element is routed
/// to the enclosing point, track, or route, and so `<ele>`/`<time>` land on the
/// point currently being built.
fn parse_gpx(text: &str) -> Result<ParsedGpx, ToolError> {
    let mut reader = Reader::from_str(text);
    reader.config_mut().trim_text(true);

    let mut waypoints: Vec<Waypoint> = Vec::new();
    let mut lines: Vec<GpxLine> = Vec::new();

    // Point currently being built (from wpt/trkpt/rtept).
    let mut cur_pt: Option<PendingPoint> = None;
    // Line being assembled (track segment or route) and its name.
    let mut cur_line: Option<Vec<(f64, f64)>> = None;
    let mut track_name: Option<String> = None;
    let mut route_name: Option<String> = None;
    let mut in_trk = false;
    let mut in_rte = false;
    // Child text element currently being collected ("ele"/"time"/"name").
    let mut collecting: Option<Vec<u8>> = None;
    let mut text_buf = String::new();

    // Handles a Start (or Empty, via the caller) element open.
    macro_rules! on_open {
        ($e:expr) => {{
            let name = $e.local_name().as_ref().to_vec();
            match name.as_slice() {
                b"wpt" | b"trkpt" | b"rtept" => {
                    let (lat, lon) = read_lat_lon(&$e);
                    cur_pt = Some(PendingPoint {
                        lat,
                        lon,
                        ..Default::default()
                    });
                }
                b"trk" => {
                    in_trk = true;
                    track_name = None;
                }
                b"trkseg" => {
                    cur_line = Some(Vec::new());
                }
                b"rte" => {
                    in_rte = true;
                    route_name = None;
                    cur_line = Some(Vec::new());
                }
                b"ele" | b"time" | b"name" => {
                    collecting = Some(name);
                    text_buf.clear();
                }
                _ => {}
            }
        }};
    }

    // Handles a matching End element close.
    macro_rules! on_close {
        ($name:expr) => {{
            match $name.as_slice() {
                b"ele" => {
                    if let Some(p) = cur_pt.as_mut() {
                        p.ele = text_buf.trim().parse::<f64>().ok();
                    }
                    collecting = None;
                }
                b"time" => {
                    if let Some(p) = cur_pt.as_mut() {
                        let t = text_buf.trim();
                        if !t.is_empty() {
                            p.time = Some(t.to_string());
                        }
                    }
                    collecting = None;
                }
                b"name" => {
                    let t = text_buf.trim().to_string();
                    if let Some(p) = cur_pt.as_mut() {
                        if !t.is_empty() {
                            p.name = Some(t);
                        }
                    } else if in_rte {
                        if !t.is_empty() {
                            route_name = Some(t);
                        }
                    } else if in_trk {
                        if !t.is_empty() {
                            track_name = Some(t);
                        }
                    }
                    collecting = None;
                }
                b"wpt" => {
                    if let Some(PendingPoint {
                        lat: Some(lat),
                        lon: Some(lon),
                        ele,
                        time,
                        name,
                    }) = cur_pt.take()
                    {
                        waypoints.push(Waypoint {
                            lat,
                            lon,
                            ele,
                            time,
                            name,
                        });
                    }
                }
                b"trkpt" | b"rtept" => {
                    if let Some(PendingPoint {
                        lat: Some(lat),
                        lon: Some(lon),
                        ..
                    }) = cur_pt.take()
                    {
                        if let Some(line) = cur_line.as_mut() {
                            line.push((lon, lat));
                        }
                    }
                }
                b"trkseg" => {
                    if let Some(coords) = cur_line.take() {
                        lines.push(GpxLine {
                            kind: LineKind::Track,
                            name: track_name.clone(),
                            coords,
                        });
                    }
                }
                b"trk" => {
                    in_trk = false;
                    track_name = None;
                }
                b"rte" => {
                    if let Some(coords) = cur_line.take() {
                        lines.push(GpxLine {
                            kind: LineKind::Route,
                            name: route_name.clone(),
                            coords,
                        });
                    }
                    in_rte = false;
                    route_name = None;
                }
                _ => {}
            }
        }};
    }

    loop {
        match reader.read_event() {
            Ok(Event::Start(e)) => on_open!(e),
            Ok(Event::Empty(e)) => {
                let name = e.local_name().as_ref().to_vec();
                on_open!(e);
                on_close!(name);
            }
            Ok(Event::Text(e)) => {
                if collecting.is_some() {
                    let t = e
                        .decode()
                        .map_err(|err| ToolError::Execution(format!("GPX text decode: {err}")))?;
                    text_buf.push_str(&t);
                }
            }
            Ok(Event::End(e)) => {
                let name = e.local_name().as_ref().to_vec();
                on_close!(name);
            }
            Ok(Event::Eof) => break,
            Err(err) => {
                return Err(ToolError::Execution(format!(
                    "failed parsing GPX at byte {}: {err}",
                    reader.buffer_position()
                )))
            }
            _ => {}
        }
    }

    Ok(ParsedGpx { waypoints, lines })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::Layer as VLayer;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    const SAMPLE: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<gpx version="1.1" creator="test" xmlns="http://www.topografix.com/GPX/1/1">
  <wpt lat="35.9606" lon="-83.9207">
    <ele>280.5</ele>
    <time>2026-07-22T14:00:00Z</time>
    <name>Trailhead</name>
  </wpt>
  <wpt lat="35.9650" lon="-83.9100">
    <ele>310.0</ele>
    <name>Overlook</name>
  </wpt>
  <trk>
    <name>Morning Hike</name>
    <trkseg>
      <trkpt lat="35.9606" lon="-83.9207"><ele>280.5</ele></trkpt>
      <trkpt lat="35.9620" lon="-83.9180"><ele>290.0</ele></trkpt>
      <trkpt lat="35.9650" lon="-83.9100"><ele>310.0</ele></trkpt>
    </trkseg>
  </trk>
  <rte>
    <name>Planned Route</name>
    <rtept lat="35.9606" lon="-83.9207"/>
    <rtept lat="35.9700" lon="-83.9000"/>
  </rte>
</gpx>"#;

    fn write_gpx(body: &str) -> std::path::PathBuf {
        use std::sync::atomic::{AtomicU64, Ordering};
        static SEQ: AtomicU64 = AtomicU64::new(0);
        let n = SEQ.fetch_add(1, Ordering::Relaxed);
        let p = std::env::temp_dir().join(format!("gpx_test_{}_{n}.gpx", std::process::id()));
        std::fs::write(&p, body).unwrap();
        p
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        GpxToFeaturesTool.run(&args, &ctx()).unwrap()
    }

    fn load(path: &str) -> VLayer {
        crate::vector_common::load_input_layer(path).unwrap()
    }

    #[test]
    fn waypoint_count_equals_wpt_count() {
        let p = write_gpx(SAMPLE);
        let out = run(json!({ "input": p.to_str().unwrap() }));
        // Two <wpt> in the fixture.
        assert_eq!(out.outputs["waypoint_count"], json!(2));
        assert_eq!(out.outputs["point_count"], json!(2));
        let layer = load(out.outputs["output"].as_str().unwrap());
        assert_eq!(layer.len(), 2);
        let f = &layer.features[0];
        assert_eq!(
            f.get(&layer.schema, "name").unwrap(),
            &FieldValue::Text("Trailhead".into())
        );
        assert_eq!(
            f.get(&layer.schema, "ele").unwrap(),
            &FieldValue::Float(280.5)
        );
        assert_eq!(
            f.get(&layer.schema, "time").unwrap(),
            &FieldValue::Text("2026-07-22T14:00:00Z".into())
        );
        assert!(matches!(f.geometry, Some(Geometry::Point(_))));
    }

    #[test]
    fn track_vertex_count_equals_trkpt_count() {
        let p = write_gpx(SAMPLE);
        let lines_path = write_gpx("").with_extension("lines.geojson");
        let out = run(json!({
            "input": p.to_str().unwrap(),
            "lines_output": lines_path.to_str().unwrap(),
        }));
        // One track (3 trkpt) + one route (2 rtept) = 2 lines.
        assert_eq!(out.outputs["line_count"], json!(2));
        // Track vertices only: 3 trkpt.
        assert_eq!(out.outputs["track_vertex_count"], json!(3));
        let layer = load(out.outputs["lines_output"].as_str().unwrap());
        assert_eq!(layer.len(), 2);
        // Total vertices across both lines: 3 + 2 = 5.
        let total: usize = layer
            .features
            .iter()
            .map(|f| match &f.geometry {
                Some(Geometry::LineString(cs)) => cs.len(),
                _ => 0,
            })
            .sum();
        assert_eq!(total, 5);
        // The track carries its name and type.
        let trk = layer
            .features
            .iter()
            .find(|f| f.get(&layer.schema, "type").unwrap() == &FieldValue::Text("track".into()))
            .unwrap();
        assert_eq!(
            trk.get(&layer.schema, "name").unwrap(),
            &FieldValue::Text("Morning Hike".into())
        );
        assert_eq!(
            trk.get(&layer.schema, "n_points").unwrap(),
            &FieldValue::Integer(3)
        );
    }

    #[test]
    fn handles_self_closing_waypoint() {
        // Waypoint with no children, self-closed.
        let body = r#"<gpx version="1.1"><wpt lat="1.0" lon="2.0"/></gpx>"#;
        let p = write_gpx(body);
        let out = run(json!({ "input": p.to_str().unwrap() }));
        assert_eq!(out.outputs["waypoint_count"], json!(1));
        let layer = load(out.outputs["output"].as_str().unwrap());
        if let Some(Geometry::Point(c)) = &layer.features[0].geometry {
            assert!((c.x - 2.0).abs() < 1e-9 && (c.y - 1.0).abs() < 1e-9);
        } else {
            panic!("expected point");
        }
    }

    #[test]
    fn rejects_missing_input() {
        let tool = GpxToFeaturesTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "/some/file.gpx" })).is_ok());
    }
}
