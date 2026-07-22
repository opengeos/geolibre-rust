//! GeoLibre tool: convert a GTFS transit feed into GIS features.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *GTFS Stops To Features* /
//! *GTFS Shapes To Features* (Conversion) plus *Calculate Transit Service
//! Frequency* (Public Transit): parse a standard GTFS feed into stop **points**
//! (from `stops.txt`) and route **lines** (from `shapes.txt`, attributed via
//! `trips.txt`/`routes.txt`), optionally with per-stop service-frequency counts
//! from `stop_times.txt`.
//!
//! GTFS is the de-facto open standard for public-transit schedules, and neither
//! the repo nor the bundled whitebox-wasm suite touches transit data. The work
//! is pure CSV parsing plus geometry assembly — squarely in-stack — and the
//! GeoJSON/PMTiles outputs drop straight into the web map.
//!
//! `input` is a path to an **extracted** GTFS feed directory (a folder of the
//! `*.txt` files). GTFS coordinates are always WGS84, so output is EPSG:4326.
//! v1 reads an unzipped feed; zip support is future work.

use std::collections::{BTreeMap, HashMap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{parse_optional_str, write_or_store_layer};

pub struct GtfsToFeaturesTool;

impl Tool for GtfsToFeaturesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "gtfs_to_features",
            display_name: "GTFS To Features",
            summary: "Parse a GTFS transit feed into stop points and route lines, optionally with per-stop service-frequency counts.",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Path to an extracted GTFS feed directory (folder of *.txt files: stops.txt, shapes.txt, trips.txt, ...).",
                    required: true,
                },
                ToolParamSpec {
                    name: "stops_output",
                    description: "Optional output point vector path for stops (driver from its extension). If neither output is given, stops are stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "shapes_output",
                    description: "Optional output line vector path for route shapes (driver from its extension).",
                    required: false,
                },
                ToolParamSpec {
                    name: "frequency",
                    description: "When true, count each stop's scheduled departures from stop_times.txt into an n_departures field. Default false.",
                    required: false,
                },
                ToolParamSpec {
                    name: "start_time",
                    description: "Optional HH:MM:SS lower bound for the frequency window (GTFS times may exceed 24:00:00).",
                    required: false,
                },
                ToolParamSpec {
                    name: "end_time",
                    description: "Optional HH:MM:SS upper bound for the frequency window.",
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
        let dir = args
            .get("input")
            .and_then(Value::as_str)
            .filter(|s| !s.trim().is_empty())
            .ok_or_else(|| {
                ToolError::Validation("missing required parameter 'input'".to_string())
            })?;
        let stops_output = parse_optional_str(args, "stops_output")?;
        let shapes_output = parse_optional_str(args, "shapes_output")?;
        let frequency = parse_optional_bool(args, "frequency")?.unwrap_or(false);
        let win_start = parse_optional_str(args, "start_time")?
            .map(parse_gtfs_time)
            .transpose()?;
        let win_end = parse_optional_str(args, "end_time")?
            .map(parse_gtfs_time)
            .transpose()?;

        // Optional per-stop departure counts.
        let mut departures: HashMap<String, u64> = HashMap::new();
        if frequency {
            let table = read_table(dir, "stop_times.txt")?;
            let sid = table.col("stop_id")?;
            let dep = table.col("departure_time")?;
            for row in &table.rows {
                let stop = row.get(sid).cloned().unwrap_or_default();
                if let Some(t) = row.get(dep).and_then(|s| parse_gtfs_time(s).ok()) {
                    if win_start.map(|w| t >= w).unwrap_or(true)
                        && win_end.map(|w| t <= w).unwrap_or(true)
                    {
                        *departures.entry(stop).or_insert(0) += 1;
                    }
                }
            }
        }

        // ── Stops layer ─────────────────────────────────────────────────────
        let stops_tbl = read_table(dir, "stops.txt")?;
        let (i_id, i_name, i_lat, i_lon) = (
            stops_tbl.col("stop_id")?,
            stops_tbl.col("stop_name").unwrap_or(usize::MAX),
            stops_tbl.col("stop_lat")?,
            stops_tbl.col("stop_lon")?,
        );
        let i_loc = stops_tbl.col("location_type").ok();

        let mut stops = Layer::new("gtfs_stops")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        stops.add_field(FieldDef::new("stop_id", FieldType::Text));
        stops.add_field(FieldDef::new("stop_name", FieldType::Text));
        stops.add_field(FieldDef::new("location_type", FieldType::Integer));
        if frequency {
            stops.add_field(FieldDef::new("n_departures", FieldType::Integer));
        }
        for row in &stops_tbl.rows {
            let (Some(lat), Some(lon)) = (
                row.get(i_lat).and_then(|s| s.trim().parse::<f64>().ok()),
                row.get(i_lon).and_then(|s| s.trim().parse::<f64>().ok()),
            ) else {
                continue;
            };
            let id = row.get(i_id).cloned().unwrap_or_default();
            let name = if i_name != usize::MAX {
                row.get(i_name).cloned().unwrap_or_default()
            } else {
                String::new()
            };
            let mut attrs = vec![
                ("stop_id", FieldValue::Text(id.clone())),
                ("stop_name", FieldValue::Text(name)),
                (
                    "location_type",
                    i_loc
                        .and_then(|li| row.get(li))
                        .and_then(|s| s.trim().parse::<i64>().ok())
                        .map(FieldValue::Integer)
                        .unwrap_or(FieldValue::Integer(0)),
                ),
            ];
            if frequency {
                attrs.push((
                    "n_departures",
                    FieldValue::Integer(*departures.get(&id).unwrap_or(&0) as i64),
                ));
            }
            stops
                .add_feature(Some(Geometry::point(lon, lat)), &attrs)
                .map_err(|e| ToolError::Execution(format!("failed adding stop: {e}")))?;
        }
        let stop_count = stops.len();
        ctx.progress.info(&format!("parsed {stop_count} stop(s)"));

        // ── Shapes layer (optional; only if shapes.txt exists) ──────────────
        let mut shape_count = 0usize;
        let mut shapes_path: Option<String> = None;
        if shapes_output.is_some() {
            if let Ok(shapes_tbl) = read_table(dir, "shapes.txt") {
                let shapes = build_shapes(dir, &shapes_tbl)?;
                shape_count = shapes.len();
                let path = write_or_store_layer(shapes, shapes_output)?;
                shapes_path = Some(path);
            } else {
                ctx.progress
                    .info("shapes.txt not found; skipping route shapes");
            }
        }

        // Write stops unless the caller asked only for shapes. When no output
        // path is given the layer is stored in memory and its handle returned.
        let write_stops = stops_output.is_some() || shapes_output.is_none();
        let stops_path = if write_stops {
            Some(write_or_store_layer(stops, stops_output)?)
        } else {
            None
        };

        let mut outputs = BTreeMap::new();
        if let Some(p) = &stops_path {
            outputs.insert("output".to_string(), json!(p));
            outputs.insert("stops_output".to_string(), json!(p));
        }
        if let Some(p) = &shapes_path {
            outputs.insert("shapes_output".to_string(), json!(p));
        }
        outputs.insert("stop_count".to_string(), json!(stop_count));
        outputs.insert("shape_count".to_string(), json!(shape_count));
        if frequency {
            outputs.insert(
                "total_departures".to_string(),
                json!(departures.values().sum::<u64>()),
            );
        }
        Ok(ToolRunResult { outputs })
    }
}

/// Builds route-shape polylines from `shapes.txt`, attributed with the route
/// each shape belongs to (via `trips.txt` + `routes.txt` when present).
fn build_shapes(dir: &str, shapes_tbl: &Table) -> Result<Layer, ToolError> {
    let s_id = shapes_tbl.col("shape_id")?;
    let s_lat = shapes_tbl.col("shape_pt_lat")?;
    let s_lon = shapes_tbl.col("shape_pt_lon")?;
    let s_seq = shapes_tbl.col("shape_pt_sequence")?;

    // shape_id -> Vec<(sequence, lon, lat)>
    let mut by_shape: BTreeMap<String, Vec<(i64, f64, f64)>> = BTreeMap::new();
    for row in &shapes_tbl.rows {
        let (Some(lat), Some(lon), Some(seq)) = (
            row.get(s_lat).and_then(|s| s.trim().parse::<f64>().ok()),
            row.get(s_lon).and_then(|s| s.trim().parse::<f64>().ok()),
            row.get(s_seq).and_then(|s| s.trim().parse::<i64>().ok()),
        ) else {
            continue;
        };
        let id = row.get(s_id).cloned().unwrap_or_default();
        by_shape.entry(id).or_default().push((seq, lon, lat));
    }

    // shape_id -> route_id (first trip using the shape).
    let mut shape_route: HashMap<String, String> = HashMap::new();
    if let Ok(trips) = read_table(dir, "trips.txt") {
        if let (Ok(t_route), Ok(t_shape)) = (trips.col("route_id"), trips.col("shape_id")) {
            for row in &trips.rows {
                if let (Some(sh), Some(rt)) = (row.get(t_shape), row.get(t_route)) {
                    if !sh.is_empty() {
                        shape_route.entry(sh.clone()).or_insert_with(|| rt.clone());
                    }
                }
            }
        }
    }
    // route_id -> display name.
    let mut route_name: HashMap<String, String> = HashMap::new();
    if let Ok(routes) = read_table(dir, "routes.txt") {
        if let Ok(r_id) = routes.col("route_id") {
            let short = routes.col("route_short_name").ok();
            let long = routes.col("route_long_name").ok();
            for row in &routes.rows {
                let id = row.get(r_id).cloned().unwrap_or_default();
                let name = short
                    .and_then(|c| row.get(c))
                    .filter(|s| !s.is_empty())
                    .or_else(|| long.and_then(|c| row.get(c)))
                    .cloned()
                    .unwrap_or_default();
                route_name.insert(id, name);
            }
        }
    }

    let mut layer = Layer::new("gtfs_shapes")
        .with_geom_type(GeometryType::LineString)
        .with_crs_epsg(4326);
    layer.add_field(FieldDef::new("shape_id", FieldType::Text));
    layer.add_field(FieldDef::new("route_id", FieldType::Text));
    layer.add_field(FieldDef::new("route_name", FieldType::Text));
    for (id, mut pts) in by_shape {
        if pts.len() < 2 {
            continue;
        }
        pts.sort_by_key(|&(seq, _, _)| seq);
        let coords: Vec<Coord> = pts
            .iter()
            .map(|&(_, lon, lat)| Coord::xy(lon, lat))
            .collect();
        let route = shape_route.get(&id).cloned().unwrap_or_default();
        let rname = route_name.get(&route).cloned().unwrap_or_default();
        layer
            .add_feature(
                Some(Geometry::LineString(coords)),
                &[
                    ("shape_id", FieldValue::Text(id)),
                    ("route_id", FieldValue::Text(route)),
                    ("route_name", FieldValue::Text(rname)),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed adding shape: {e}")))?;
    }
    Ok(layer)
}

// ── Minimal CSV table ────────────────────────────────────────────────────────

/// A parsed CSV file: header column names + rows of string cells.
struct Table {
    header: HashMap<String, usize>,
    rows: Vec<Vec<String>>,
}

impl Table {
    /// Column index for `name`, or a validation error if absent.
    fn col(&self, name: &str) -> Result<usize, ToolError> {
        self.header
            .get(name)
            .copied()
            .ok_or_else(|| ToolError::Validation(format!("GTFS column '{name}' not found")))
    }
}

/// Reads and parses one GTFS `*.txt` CSV file from the feed directory.
fn read_table(dir: &str, file: &str) -> Result<Table, ToolError> {
    let path = std::path::Path::new(dir).join(file);
    let text = std::fs::read_to_string(&path)
        .map_err(|e| ToolError::Execution(format!("failed reading {}: {e}", path.display())))?;
    // Strip a UTF-8 BOM if present.
    let text = text.strip_prefix('\u{feff}').unwrap_or(&text);
    let mut lines = text.lines();
    let header_line = lines
        .next()
        .ok_or_else(|| ToolError::Execution(format!("{file} is empty")))?;
    let header: HashMap<String, usize> = split_csv(header_line)
        .into_iter()
        .enumerate()
        .map(|(i, name)| (name.trim().to_string(), i))
        .collect();
    let rows: Vec<Vec<String>> = lines
        .filter(|l| !l.trim().is_empty())
        .map(split_csv)
        .collect();
    Ok(Table { header, rows })
}

/// Splits one CSV line, honouring double-quoted fields (with `""` escapes) and
/// a trailing carriage return.
fn split_csv(line: &str) -> Vec<String> {
    let line = line.strip_suffix('\r').unwrap_or(line);
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut in_quotes = false;
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        if in_quotes {
            if c == '"' {
                if chars.peek() == Some(&'"') {
                    cur.push('"');
                    chars.next();
                } else {
                    in_quotes = false;
                }
            } else {
                cur.push(c);
            }
        } else {
            match c {
                '"' => in_quotes = true,
                ',' => out.push(std::mem::take(&mut cur)),
                _ => cur.push(c),
            }
        }
    }
    out.push(cur);
    out
}

/// Parses a GTFS `HH:MM:SS` time to seconds since midnight (hours may exceed 24).
fn parse_gtfs_time(s: &str) -> Result<i64, ToolError> {
    let parts: Vec<&str> = s.trim().split(':').collect();
    if parts.len() != 3 {
        return Err(ToolError::Validation(format!(
            "time '{s}' must be HH:MM:SS"
        )));
    }
    let h: i64 = parts[0].parse().map_err(|_| bad_time(s))?;
    let m: i64 = parts[1].parse().map_err(|_| bad_time(s))?;
    let sec: i64 = parts[2].parse().map_err(|_| bad_time(s))?;
    Ok(h * 3600 + m * 60 + sec)
}

fn bad_time(s: &str) -> ToolError {
    ToolError::Validation(format!("time '{s}' must be HH:MM:SS"))
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
    use wbvector::Layer as VLayer;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Writes a tiny GTFS feed to a unique temp dir and returns its path.
    ///
    /// The directory is keyed by process id *and* a per-call atomic counter, so
    /// concurrently-running tests never share (and clobber) the same feed.
    fn write_feed() -> std::path::PathBuf {
        use std::sync::atomic::{AtomicU64, Ordering};
        static SEQ: AtomicU64 = AtomicU64::new(0);
        let n = SEQ.fetch_add(1, Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!("gtfs_test_{}_{n}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(
            dir.join("stops.txt"),
            "stop_id,stop_name,stop_lat,stop_lon,location_type\n\
             A,\"Alpha, Main St\",40.0,-74.0,0\n\
             B,Beta,40.1,-74.1,0\n",
        )
        .unwrap();
        std::fs::write(
            dir.join("shapes.txt"),
            "shape_id,shape_pt_lat,shape_pt_lon,shape_pt_sequence\n\
             S1,40.0,-74.0,1\n\
             S1,40.05,-74.05,3\n\
             S1,40.1,-74.1,2\n",
        )
        .unwrap();
        std::fs::write(
            dir.join("trips.txt"),
            "route_id,trip_id,shape_id\nR1,T1,S1\n",
        )
        .unwrap();
        std::fs::write(
            dir.join("routes.txt"),
            "route_id,route_short_name,route_long_name,route_type\nR1,1,Red Line,1\n",
        )
        .unwrap();
        std::fs::write(
            dir.join("stop_times.txt"),
            "trip_id,departure_time,stop_id,stop_sequence\n\
             T1,06:00:00,A,1\n\
             T1,06:10:00,B,2\n\
             T1,07:00:00,A,1\n",
        )
        .unwrap();
        dir
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        GtfsToFeaturesTool.run(&args, &ctx()).unwrap()
    }

    fn load(path: &str) -> VLayer {
        crate::vector_common::load_input_layer(path).unwrap()
    }

    #[test]
    fn parses_stops_as_points() {
        let dir = write_feed();
        let out = run(json!({ "input": dir.to_str().unwrap() }));
        assert_eq!(out.outputs["stop_count"], json!(2));
        let layer = load(out.outputs["output"].as_str().unwrap());
        assert_eq!(layer.len(), 2);
        // The quoted name with an embedded comma parses as one field.
        let a = &layer.features[0];
        assert_eq!(
            a.get(&layer.schema, "stop_name").unwrap(),
            &FieldValue::Text("Alpha, Main St".into())
        );
        assert!(matches!(a.geometry, Some(Geometry::Point(_))));
    }

    #[test]
    fn builds_shapes_ordered_by_sequence() {
        let dir = write_feed();
        let shapes_path = dir.join("out_shapes.geojson");
        let out = run(json!({
            "input": dir.to_str().unwrap(),
            "shapes_output": shapes_path.to_str().unwrap(),
        }));
        assert_eq!(out.outputs["shape_count"], json!(1));
        let layer = load(out.outputs["shapes_output"].as_str().unwrap());
        let f = &layer.features[0];
        // Route joined from trips + routes.
        assert_eq!(
            f.get(&layer.schema, "route_name").unwrap(),
            &FieldValue::Text("1".into())
        );
        // Vertices are ordered by shape_pt_sequence (1,2,3), not file order.
        if let Some(Geometry::LineString(cs)) = &f.geometry {
            assert_eq!(cs.len(), 3);
            assert!(
                (cs[1].y - 40.1).abs() < 1e-9,
                "seq 2 should come before seq 3"
            );
        } else {
            panic!("expected LineString");
        }
    }

    #[test]
    fn frequency_counts_departures() {
        let dir = write_feed();
        let out = run(json!({ "input": dir.to_str().unwrap(), "frequency": true }));
        let layer = load(out.outputs["output"].as_str().unwrap());
        let dep = |name: &str| {
            layer
                .features
                .iter()
                .find(|f| {
                    f.get(&layer.schema, "stop_id").unwrap() == &FieldValue::Text(name.into())
                })
                .unwrap()
                .get(&layer.schema, "n_departures")
                .unwrap()
                .as_i64()
                .unwrap()
        };
        // Stop A has two departures (06:00, 07:00); B has one.
        assert_eq!(dep("A"), 2);
        assert_eq!(dep("B"), 1);
        assert_eq!(out.outputs["total_departures"], json!(3));
    }

    #[test]
    fn frequency_window_filters() {
        let dir = write_feed();
        let out = run(json!({
            "input": dir.to_str().unwrap(),
            "frequency": true,
            "start_time": "06:30:00",
            "end_time": "08:00:00",
        }));
        // Only A's 07:00 departure falls in the window.
        assert_eq!(out.outputs["total_departures"], json!(1));
    }

    #[test]
    fn rejects_missing_input() {
        let tool = GtfsToFeaturesTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "/some/dir" })).is_ok());
    }
}
