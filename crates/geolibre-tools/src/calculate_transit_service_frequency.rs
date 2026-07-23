//! GeoLibre tool: windowed transit service-frequency from a GTFS feed.
//!
//! Pure-Rust counterpart of ArcGIS Public Transit's *Calculate Transit Service
//! Frequency*. It extends the repo's GTFS suite (`gtfs_to_features`,
//! `features_to_gtfs`): for a chosen service **date** (resolved through
//! `calendar.txt` / `calendar_dates.txt`) and a **time window**, it counts the
//! scheduled transit trips serving each stop (or route line) and reports both
//! the raw count and trips-per-hour — the basis for transit accessibility maps.
//!
//! Unlike `gtfs_to_features` (which counts departures over an optional raw
//! time-of-day window with no calendar logic), this resolves which services
//! actually run on a given date, counts arrivals or departures, and can
//! aggregate to route polylines.
//!
//! `input` is a path to an extracted GTFS feed directory. Output is EPSG:4326.

use std::collections::{BTreeMap, HashMap, HashSet};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{parse_optional_str, write_or_store_layer};

pub struct CalculateTransitServiceFrequencyTool;

impl Tool for CalculateTransitServiceFrequencyTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "calculate_transit_service_frequency",
            display_name: "Calculate Transit Service Frequency",
            summary: "Windowed GTFS transit service frequency: for a service date (via calendar/calendar_dates) and time window, count trips serving each stop or route line and report trips-per-hour, like ArcGIS Calculate Transit Service Frequency.",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Path to an extracted GTFS feed directory.",
                    required: true,
                },
                ToolParamSpec {
                    name: "target",
                    description: "'stops' (point per stop, default) or 'lines' (polyline per route).",
                    required: false,
                },
                ToolParamSpec {
                    name: "date",
                    description: "Service date YYYYMMDD, resolved via calendar.txt/calendar_dates.txt. If omitted, all services count.",
                    required: false,
                },
                ToolParamSpec {
                    name: "start_time",
                    description: "Window start HH:MM or HH:MM:SS (GTFS times may exceed 24:00). If omitted, the whole day counts.",
                    required: false,
                },
                ToolParamSpec {
                    name: "duration_minutes",
                    description: "Window length in minutes (used with start_time; also sets the per-hour denominator).",
                    required: false,
                },
                ToolParamSpec {
                    name: "count",
                    description: "'departures' (default) or 'arrivals' — which stop-time field to count.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_target(args)?;
        parse_count(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let dir = require_str(args, "input")?;
        let target = parse_target(args)?;
        let use_arrivals = matches!(parse_count(args)?, CountField::Arrivals);
        let date = parse_optional_str(args, "date")?.map(str::to_string);
        let win_start = parse_optional_str(args, "start_time")?
            .map(parse_gtfs_time)
            .transpose()?;
        let duration = parse_optional_f64(args, "duration_minutes")?;
        let win_end = match (win_start, duration) {
            (Some(s), Some(d)) => Some(s + (d * 60.0) as i64),
            _ => None,
        };
        let hours = duration.map(|d| (d / 60.0).max(1e-9)).unwrap_or(1.0);
        let time_col = if use_arrivals {
            "arrival_time"
        } else {
            "departure_time"
        };

        // trips.txt -> trip -> (service_id, route_id, shape_id)
        let trips_tbl = read_table(dir, "trips.txt")?;
        let t_trip = trips_tbl.col("trip_id")?;
        let t_service = trips_tbl.col("service_id").ok();
        let t_route = trips_tbl.col("route_id").ok();
        let t_shape = trips_tbl.col("shape_id").ok();
        let mut trip_service: HashMap<String, String> = HashMap::new();
        let mut trip_route: HashMap<String, String> = HashMap::new();
        for row in &trips_tbl.rows {
            let id = row.get(t_trip).cloned().unwrap_or_default();
            if let Some(si) = t_service {
                trip_service.insert(id.clone(), row.get(si).cloned().unwrap_or_default());
            }
            if let Some(ri) = t_route {
                trip_route.insert(id.clone(), row.get(ri).cloned().unwrap_or_default());
            }
        }

        // Active services on `date` (all services if no date given).
        let active: Option<HashSet<String>> = match &date {
            Some(d) => Some(active_services(dir, d)?),
            None => None,
        };
        let is_active = |trip: &str| -> bool {
            match &active {
                None => true,
                Some(set) => trip_service
                    .get(trip)
                    .map(|s| set.contains(s))
                    .unwrap_or(false),
            }
        };
        let in_window = |t: i64| -> bool {
            win_start.map(|w| t >= w).unwrap_or(true) && win_end.map(|w| t < w).unwrap_or(true)
        };

        // Count qualifying stop-time events.
        let st = read_table(dir, "stop_times.txt")?;
        let s_trip = st.col("trip_id")?;
        let s_stop = st.col("stop_id")?;
        let s_time = st.col(time_col)?;

        let mut stop_counts: HashMap<String, u64> = HashMap::new();
        // For lines: earliest qualifying trip time per trip.
        let mut qualifying_trips: HashSet<String> = HashSet::new();
        for row in &st.rows {
            let trip = row.get(s_trip).cloned().unwrap_or_default();
            if !is_active(&trip) {
                continue;
            }
            let Some(t) = row.get(s_time).and_then(|s| parse_gtfs_time(s).ok()) else {
                continue;
            };
            if !in_window(t) {
                continue;
            }
            let stop = row.get(s_stop).cloned().unwrap_or_default();
            *stop_counts.entry(stop).or_insert(0) += 1;
            qualifying_trips.insert(trip);
        }

        let (layer, feature_count) = match target {
            Target::Stops => build_stops(dir, &stop_counts, hours)?,
            Target::Lines => build_lines(
                dir,
                &qualifying_trips,
                &trip_route,
                &trips_tbl,
                t_trip,
                t_shape,
                hours,
            )?,
        };

        ctx.progress
            .info(&format!("{feature_count} feature(s) for target {target:?}"));

        let total: u64 = stop_counts.values().sum();
        let out_path = write_or_store_layer(layer, parse_optional_str(args, "output")?)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("total_events".to_string(), json!(total));
        outputs.insert(
            "qualifying_trips".to_string(),
            json!(qualifying_trips.len()),
        );
        Ok(ToolRunResult { outputs })
    }
}

/// Builds the per-stop point layer with counts and trips-per-hour.
fn build_stops(
    dir: &str,
    counts: &HashMap<String, u64>,
    hours: f64,
) -> Result<(Layer, usize), ToolError> {
    let stops = read_table(dir, "stops.txt")?;
    let i_id = stops.col("stop_id")?;
    let i_name = stops.col("stop_name").unwrap_or(usize::MAX);
    let i_lat = stops.col("stop_lat")?;
    let i_lon = stops.col("stop_lon")?;

    let mut layer = Layer::new("transit_frequency_stops")
        .with_geom_type(GeometryType::Point)
        .with_crs_epsg(4326);
    layer.add_field(FieldDef::new("stop_id", FieldType::Text));
    layer.add_field(FieldDef::new("stop_name", FieldType::Text));
    layer.add_field(FieldDef::new("n_trips", FieldType::Integer));
    layer.add_field(FieldDef::new("trips_per_hour", FieldType::Float));

    for row in &stops.rows {
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
        let n = *counts.get(&id).unwrap_or(&0);
        layer
            .add_feature(
                Some(Geometry::point(lon, lat)),
                &[
                    ("stop_id", FieldValue::Text(id)),
                    ("stop_name", FieldValue::Text(name)),
                    ("n_trips", FieldValue::Integer(n as i64)),
                    ("trips_per_hour", FieldValue::Float(n as f64 / hours)),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed adding stop: {e}")))?;
    }
    let n = layer.len();
    Ok((layer, n))
}

/// Builds a per-route polyline layer whose frequency is the number of qualifying
/// trips on the route within the window.
fn build_lines(
    dir: &str,
    qualifying_trips: &HashSet<String>,
    trip_route: &HashMap<String, String>,
    trips_tbl: &Table,
    t_trip: usize,
    t_shape: Option<usize>,
    hours: f64,
) -> Result<(Layer, usize), ToolError> {
    // route -> qualifying trip count.
    let mut route_trips: HashMap<String, u64> = HashMap::new();
    for trip in qualifying_trips {
        if let Some(route) = trip_route.get(trip) {
            *route_trips.entry(route.clone()).or_insert(0) += 1;
        }
    }

    // route -> a representative shape_id (first trip's shape).
    let mut route_shape: HashMap<String, String> = HashMap::new();
    if let Some(si) = t_shape {
        for row in &trips_tbl.rows {
            let trip = row.get(t_trip).cloned().unwrap_or_default();
            if let (Some(route), Some(shape)) = (trip_route.get(&trip), row.get(si)) {
                if !shape.is_empty() {
                    route_shape
                        .entry(route.clone())
                        .or_insert_with(|| shape.clone());
                }
            }
        }
    }

    // shape_id -> ordered coords.
    let shapes: HashMap<String, Vec<Coord>> = match read_table(dir, "shapes.txt") {
        Ok(tbl) => shape_geometries(&tbl)?,
        Err(_) => HashMap::new(),
    };
    let route_name = route_names(dir);

    let mut layer = Layer::new("transit_frequency_lines")
        .with_geom_type(GeometryType::LineString)
        .with_crs_epsg(4326);
    layer.add_field(FieldDef::new("route_id", FieldType::Text));
    layer.add_field(FieldDef::new("route_name", FieldType::Text));
    layer.add_field(FieldDef::new("n_trips", FieldType::Integer));
    layer.add_field(FieldDef::new("trips_per_hour", FieldType::Float));

    let mut routes: Vec<(&String, &u64)> = route_trips.iter().collect();
    routes.sort_by(|a, b| a.0.cmp(b.0));
    for (route, n) in routes {
        let coords = route_shape
            .get(route)
            .and_then(|sid| shapes.get(sid))
            .cloned();
        let Some(coords) = coords else { continue };
        if coords.len() < 2 {
            continue;
        }
        let name = route_name.get(route).cloned().unwrap_or_default();
        layer
            .add_feature(
                Some(Geometry::LineString(coords)),
                &[
                    ("route_id", FieldValue::Text(route.clone())),
                    ("route_name", FieldValue::Text(name)),
                    ("n_trips", FieldValue::Integer(*n as i64)),
                    ("trips_per_hour", FieldValue::Float(*n as f64 / hours)),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed adding route: {e}")))?;
    }
    let n = layer.len();
    Ok((layer, n))
}

/// shape_id -> coords ordered by shape_pt_sequence.
fn shape_geometries(tbl: &Table) -> Result<HashMap<String, Vec<Coord>>, ToolError> {
    let s_id = tbl.col("shape_id")?;
    let s_lat = tbl.col("shape_pt_lat")?;
    let s_lon = tbl.col("shape_pt_lon")?;
    let s_seq = tbl.col("shape_pt_sequence")?;
    let mut by_shape: HashMap<String, Vec<(i64, f64, f64)>> = HashMap::new();
    for row in &tbl.rows {
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
    Ok(by_shape
        .into_iter()
        .map(|(id, mut pts)| {
            pts.sort_by_key(|&(seq, _, _)| seq);
            (
                id,
                pts.into_iter()
                    .map(|(_, lon, lat)| Coord::xy(lon, lat))
                    .collect(),
            )
        })
        .collect())
}

fn route_names(dir: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
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
                map.insert(id, name);
            }
        }
    }
    map
}

/// Resolves the set of service_ids active on `date` (YYYYMMDD) from
/// `calendar.txt` (weekday + date-range) and `calendar_dates.txt` (exceptions).
fn active_services(dir: &str, date: &str) -> Result<HashSet<String>, ToolError> {
    let ymd: i64 = date
        .trim()
        .parse()
        .map_err(|_| ToolError::Validation(format!("date '{date}' must be YYYYMMDD")))?;
    let weekday = weekday_of(date)?; // 0=Mon .. 6=Sun
    let mut active: HashSet<String> = HashSet::new();

    if let Ok(cal) = read_table(dir, "calendar.txt") {
        let c_service = cal.col("service_id")?;
        let day_cols = [
            "monday",
            "tuesday",
            "wednesday",
            "thursday",
            "friday",
            "saturday",
            "sunday",
        ];
        let dcol = cal.col(day_cols[weekday]).ok();
        let start = cal.col("start_date").ok();
        let end = cal.col("end_date").ok();
        for row in &cal.rows {
            let runs = dcol
                .and_then(|c| row.get(c))
                .map(|s| s.trim() == "1")
                .unwrap_or(false);
            let in_range = {
                let s_ok = start
                    .and_then(|c| row.get(c))
                    .and_then(|s| s.trim().parse::<i64>().ok())
                    .map(|s| ymd >= s)
                    .unwrap_or(true);
                let e_ok = end
                    .and_then(|c| row.get(c))
                    .and_then(|s| s.trim().parse::<i64>().ok())
                    .map(|e| ymd <= e)
                    .unwrap_or(true);
                s_ok && e_ok
            };
            if runs && in_range {
                active.insert(row.get(c_service).cloned().unwrap_or_default());
            }
        }
    }

    if let Ok(cd) = read_table(dir, "calendar_dates.txt") {
        if let (Ok(c_service), Ok(c_date), Ok(c_ex)) = (
            cd.col("service_id"),
            cd.col("date"),
            cd.col("exception_type"),
        ) {
            for row in &cd.rows {
                if row
                    .get(c_date)
                    .map(|d| d.trim() == date.trim())
                    .unwrap_or(false)
                {
                    let service = row.get(c_service).cloned().unwrap_or_default();
                    match row.get(c_ex).map(|s| s.trim()) {
                        Some("1") => {
                            active.insert(service);
                        }
                        Some("2") => {
                            active.remove(&service);
                        }
                        _ => {}
                    }
                }
            }
        }
    }
    Ok(active)
}

/// Day of week for a YYYYMMDD string, 0 = Monday .. 6 = Sunday (Sakamoto).
fn weekday_of(date: &str) -> Result<usize, ToolError> {
    if date.trim().len() != 8 {
        return Err(ToolError::Validation(format!(
            "date '{date}' must be YYYYMMDD"
        )));
    }
    let bad = || ToolError::Validation(format!("date '{date}' must be YYYYMMDD"));
    let y: i64 = date[0..4].parse().map_err(|_| bad())?;
    let m: i64 = date[4..6].parse().map_err(|_| bad())?;
    let d: i64 = date[6..8].parse().map_err(|_| bad())?;
    if !(1..=12).contains(&m) {
        return Err(bad());
    }
    // Sakamoto's algorithm (0 = Sunday).
    let t = [0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4];
    let mut yy = y;
    if m < 3 {
        yy -= 1;
    }
    let dow_sun0 = (yy + yy / 4 - yy / 100 + yy / 400 + t[(m - 1) as usize] + d) % 7;
    // Convert Sunday=0 -> Monday=0.
    Ok(((dow_sun0 + 6) % 7) as usize)
}

// ── Minimal CSV table (shared shape with gtfs_to_features) ───────────────────

struct Table {
    header: HashMap<String, usize>,
    rows: Vec<Vec<String>>,
}

impl Table {
    fn col(&self, name: &str) -> Result<usize, ToolError> {
        self.header
            .get(name)
            .copied()
            .ok_or_else(|| ToolError::Validation(format!("GTFS column '{name}' not found")))
    }
}

fn read_table(dir: &str, file: &str) -> Result<Table, ToolError> {
    let path = std::path::Path::new(dir).join(file);
    let text = std::fs::read_to_string(&path)
        .map_err(|e| ToolError::Execution(format!("failed reading {}: {e}", path.display())))?;
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

fn parse_gtfs_time(s: &str) -> Result<i64, ToolError> {
    let parts: Vec<&str> = s.trim().split(':').collect();
    if parts.len() < 2 || parts.len() > 3 {
        return Err(ToolError::Validation(format!(
            "time '{s}' must be HH:MM[:SS]"
        )));
    }
    let h: i64 = parts[0].parse().map_err(|_| bad_time(s))?;
    let m: i64 = parts[1].parse().map_err(|_| bad_time(s))?;
    let sec: i64 = parts.get(2).map(|p| p.parse().unwrap_or(0)).unwrap_or(0);
    Ok(h * 3600 + m * 60 + sec)
}

fn bad_time(s: &str) -> ToolError {
    ToolError::Validation(format!("time '{s}' must be HH:MM[:SS]"))
}

// ── Params ──────────────────────────────────────────────────────────────────

#[derive(Debug)]
enum Target {
    Stops,
    Lines,
}

enum CountField {
    Departures,
    Arrivals,
}

fn parse_target(args: &ToolArgs) -> Result<Target, ToolError> {
    match args.get("target").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("stops") => Ok(Target::Stops),
        Some("lines") => Ok(Target::Lines),
        Some(o) => Err(ToolError::Validation(format!(
            "'target' must be 'stops' or 'lines', got '{o}'"
        ))),
    }
}

fn parse_count(args: &ToolArgs) -> Result<CountField, ToolError> {
    match args.get("count").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("departures") => Ok(CountField::Departures),
        Some("arrivals") => Ok(CountField::Arrivals),
        Some(o) => Err(ToolError::Validation(format!(
            "'count' must be 'departures' or 'arrivals', got '{o}'"
        ))),
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

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn write_feed() -> std::path::PathBuf {
        use std::sync::atomic::{AtomicU64, Ordering};
        static SEQ: AtomicU64 = AtomicU64::new(0);
        let n = SEQ.fetch_add(1, Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!("transit_freq_{}_{n}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(
            dir.join("stops.txt"),
            "stop_id,stop_name,stop_lat,stop_lon\nA,Alpha,40.0,-74.0\nB,Beta,40.1,-74.1\n",
        )
        .unwrap();
        // 20200106 is a Monday.
        std::fs::write(
            dir.join("calendar.txt"),
            "service_id,monday,tuesday,wednesday,thursday,friday,saturday,sunday,start_date,end_date\n\
             WK,1,1,1,1,1,0,0,20200101,20201231\n\
             WE,0,0,0,0,0,1,1,20200101,20201231\n",
        )
        .unwrap();
        std::fs::write(
            dir.join("trips.txt"),
            "route_id,service_id,trip_id,shape_id\n\
             R1,WK,T1,S1\nR1,WK,T2,S1\nR1,WE,T3,S1\n",
        )
        .unwrap();
        std::fs::write(
            dir.join("routes.txt"),
            "route_id,route_short_name,route_long_name\nR1,1,Red\n",
        )
        .unwrap();
        std::fs::write(
            dir.join("shapes.txt"),
            "shape_id,shape_pt_lat,shape_pt_lon,shape_pt_sequence\n\
             S1,40.0,-74.0,1\nS1,40.1,-74.1,2\n",
        )
        .unwrap();
        std::fs::write(
            dir.join("stop_times.txt"),
            "trip_id,arrival_time,departure_time,stop_id,stop_sequence\n\
             T1,06:00:00,06:00:00,A,1\nT1,06:10:00,06:10:00,B,2\n\
             T2,07:00:00,07:00:00,A,1\nT2,07:10:00,07:10:00,B,2\n\
             T3,09:00:00,09:00:00,A,1\nT3,09:10:00,09:10:00,B,2\n",
        )
        .unwrap();
        dir
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CalculateTransitServiceFrequencyTool
            .run(&args, &ctx())
            .unwrap();
        let layer = crate::vector_common::load_input_layer(out.outputs["output"].as_str().unwrap())
            .unwrap();
        (out, layer)
    }

    #[test]
    fn weekday_resolution_matches_calendar() {
        // 20200106 is a Monday -> only WK services (T1, T2) count; WE (T3) excluded.
        let dir = write_feed();
        let (out, layer) = run(json!({
            "input": dir.to_str().unwrap(), "target": "stops", "date": "20200106"
        }));
        // Stop A served by T1 and T2 only -> 2.
        let a = layer
            .features
            .iter()
            .find(|f| f.get(&layer.schema, "stop_id").unwrap() == &FieldValue::Text("A".into()))
            .unwrap();
        assert_eq!(
            a.get(&layer.schema, "n_trips").unwrap().as_i64().unwrap(),
            2
        );
        assert_eq!(out.outputs["qualifying_trips"], json!(2));
    }

    #[test]
    fn time_window_and_per_hour() {
        let dir = write_feed();
        // Window 06:00 for 60 min -> only T1's 06:00 departure at A.
        let (_o, layer) = run(json!({
            "input": dir.to_str().unwrap(), "date": "20200106",
            "start_time": "06:00", "duration_minutes": 60
        }));
        let a = layer
            .features
            .iter()
            .find(|f| f.get(&layer.schema, "stop_id").unwrap() == &FieldValue::Text("A".into()))
            .unwrap();
        assert_eq!(
            a.get(&layer.schema, "n_trips").unwrap().as_i64().unwrap(),
            1
        );
        // 1 trip / 1 hour = 1.0 per hour.
        assert!(
            (a.get(&layer.schema, "trips_per_hour")
                .unwrap()
                .as_f64()
                .unwrap()
                - 1.0)
                .abs()
                < 1e-9
        );
    }

    #[test]
    fn lines_target_aggregates_to_route() {
        let dir = write_feed();
        // No date -> all 3 trips count on route R1.
        let (out, layer) = run(json!({ "input": dir.to_str().unwrap(), "target": "lines" }));
        assert_eq!(out.outputs["feature_count"], json!(1));
        let f = &layer.features[0];
        assert_eq!(
            f.get(&layer.schema, "n_trips").unwrap().as_i64().unwrap(),
            3
        );
        assert!(matches!(f.geometry, Some(Geometry::LineString(_))));
    }

    #[test]
    fn saturday_selects_weekend_service() {
        let dir = write_feed();
        // 20200104 is a Saturday -> WE service (T3) only.
        let (out, _l) = run(json!({
            "input": dir.to_str().unwrap(), "date": "20200104"
        }));
        assert_eq!(out.outputs["qualifying_trips"], json!(1));
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CalculateTransitServiceFrequencyTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "/d", "target": "bogus" })).is_err());
        assert!(bad(json!({ "input": "/d", "count": "middle" })).is_err());
        assert!(bad(json!({ "input": "/d" })).is_ok());
    }
}
