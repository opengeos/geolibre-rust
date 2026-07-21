//! GeoLibre tool: annotate timestamped track points with motion statistics.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Calculate Motion Statistics*
//! (GeoAnalytics Desktop). The GeoLibre movement suite can reconstruct, snap,
//! and contact-trace tracks (`reconstruct_tracks`, `snap_tracks`,
//! `trace_proximity_events`, `find_space_time_matches`), but nothing annotates
//! the individual points with speed, acceleration, bearing, distance travelled,
//! or idle time — the most common track derivative. There is no equivalent in
//! the bundled whitebox suite.
//!
//! Points are grouped by `track_field`, sorted by `time_field`, and each point
//! is annotated (original attributes preserved) with, relative to the previous
//! point in its track: `seq` (0-based index), `dist` (segment length),
//! `cum_dist` (cumulative), `dt` (seconds), `elapsed` (seconds since track
//! start), `speed`, `avg_speed` (over a trailing `window` of points), `accel`
//! (change in speed / dt), `bearing` (degrees clockwise from north), and `idle`
//! (1 when the mover travelled no more than `idle_distance` over the last
//! `idle_duration` seconds). Distances are haversine metres for a geographic CRS
//! and CRS units otherwise; `time_field` is numeric seconds or ISO-8601.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CalculateMotionStatisticsTool;

impl Tool for CalculateMotionStatisticsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "calculate_motion_statistics",
            display_name: "Calculate Motion Statistics",
            summary: "Annotate timestamped track points with motion statistics — speed, acceleration, bearing, segment and cumulative distance, elapsed time and idle flag — grouped by track and ordered by time, like ArcGIS Calculate Motion Statistics.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer of timestamped observations.",
                    required: true,
                },
                ToolParamSpec {
                    name: "track_field",
                    description: "Field identifying each track/mover.",
                    required: true,
                },
                ToolParamSpec {
                    name: "time_field",
                    description: "Field holding each point's time: numeric seconds or an ISO-8601 timestamp.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output point vector (original attributes + motion fields). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "window",
                    description: "Trailing number of points for the smoothed 'avg_speed' (default 1 = instantaneous).",
                    required: false,
                },
                ToolParamSpec {
                    name: "idle_distance",
                    description: "Idle radius: a point is idle when the mover travelled no more than this over the last idle_duration seconds. Default 0 (idle detection off).",
                    required: false,
                },
                ToolParamSpec {
                    name: "idle_duration",
                    description: "Look-back window (seconds) for idle detection. Default 60.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "track_field")?;
        require_str(args, "time_field")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let track_idx = layer.schema.field_index(&prm.track_field).ok_or_else(|| {
            ToolError::Validation(format!("track_field '{}' not found", prm.track_field))
        })?;
        let time_idx = layer.schema.field_index(&prm.time_field).ok_or_else(|| {
            ToolError::Validation(format!("time_field '{}' not found", prm.time_field))
        })?;
        let geographic = layer.crs_epsg().map(|e| e == 4326).unwrap_or(true);
        let n_in_fields = layer.schema.fields().len();

        // Gather feature indices per track (keeping a handle to the source row).
        let mut tracks: BTreeMap<String, Vec<Obs>> = BTreeMap::new();
        let mut skipped = 0usize;
        for (fid, feature) in layer.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                skipped += 1;
                continue;
            };
            let Some((x, y)) = point_xy(geom) else {
                skipped += 1;
                continue;
            };
            let Some(t) = feature.attributes.get(time_idx).and_then(parse_time_value) else {
                skipped += 1;
                continue;
            };
            let id = feature
                .attributes
                .get(track_idx)
                .map(value_string)
                .unwrap_or_default();
            tracks.entry(id).or_default().push(Obs { fid, x, y, t });
        }

        ctx.progress
            .info(&format!("{} track(s) to annotate", tracks.len()));

        // Output layer: original fields copied through, then motion fields.
        let mut out = Layer::new("motion").with_geom_type(GeometryType::Point);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for f in layer.schema.fields() {
            out.add_field(f.clone());
        }
        for (name, ty) in MOTION_FIELDS {
            out.add_field(FieldDef::new(*name, *ty));
        }

        let mut annotated = 0usize;
        for (_id, mut obs) in tracks {
            obs.sort_by(|a, b| a.t.total_cmp(&b.t));
            let mut cum_dist = 0.0f64;
            let mut prev_speed = 0.0f64;
            let t0 = obs.first().map(|o| o.t).unwrap_or(0.0);
            for i in 0..obs.len() {
                let o = obs[i];
                let (dist, dt, speed, accel, bearing) = if i == 0 {
                    (0.0, 0.0, 0.0, 0.0, 0.0)
                } else {
                    let p = obs[i - 1];
                    let d = distance(p.x, p.y, o.x, o.y, geographic);
                    let dt = o.t - p.t;
                    let speed = if dt > 0.0 { d / dt } else { 0.0 };
                    let accel = if dt > 0.0 {
                        (speed - prev_speed) / dt
                    } else {
                        0.0
                    };
                    let bearing = bearing_deg(p.x, p.y, o.x, o.y, geographic);
                    (d, dt, speed, accel, bearing)
                };
                cum_dist += dist;
                prev_speed = speed;

                // Trailing-window average speed.
                let avg_speed = window_speed(&obs, i, prm.window, geographic);
                // Idle detection over the look-back window.
                let idle = if prm.idle_distance > 0.0 {
                    is_idle(&obs, i, prm.idle_distance, prm.idle_duration, geographic)
                } else {
                    false
                };

                // Copy the original attributes, then append motion fields.
                let src = &layer.features[o.fid];
                let mut attrs: Vec<FieldValue> =
                    Vec::with_capacity(n_in_fields + MOTION_FIELDS.len());
                for k in 0..n_in_fields {
                    attrs.push(src.attributes.get(k).cloned().unwrap_or(FieldValue::Null));
                }
                attrs.push(FieldValue::Integer(i as i64));
                attrs.push(FieldValue::Float(dist));
                attrs.push(FieldValue::Float(cum_dist));
                attrs.push(FieldValue::Float(dt));
                attrs.push(FieldValue::Float(o.t - t0));
                attrs.push(FieldValue::Float(speed));
                attrs.push(FieldValue::Float(avg_speed));
                attrs.push(FieldValue::Float(accel));
                attrs.push(FieldValue::Float(bearing));
                attrs.push(FieldValue::Integer(idle as i64));

                out.push(Feature {
                    fid: 0,
                    geometry: Some(Geometry::point(o.x, o.y)),
                    attributes: attrs,
                });
                annotated += 1;
            }
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("annotated".to_string(), json!(annotated));
        outputs.insert("skipped".to_string(), json!(skipped));
        Ok(ToolRunResult { outputs })
    }
}

const MOTION_FIELDS: &[(&str, FieldType)] = &[
    ("seq", FieldType::Integer),
    ("dist", FieldType::Float),
    ("cum_dist", FieldType::Float),
    ("dt", FieldType::Float),
    ("elapsed", FieldType::Float),
    ("speed", FieldType::Float),
    ("avg_speed", FieldType::Float),
    ("accel", FieldType::Float),
    ("bearing", FieldType::Float),
    ("idle", FieldType::Integer),
];

#[derive(Clone, Copy)]
struct Obs {
    fid: usize,
    x: f64,
    y: f64,
    t: f64,
}

/// Average speed over the trailing `window` segments ending at point `i`:
/// total distance / total elapsed time across those points.
fn window_speed(obs: &[Obs], i: usize, window: usize, geographic: bool) -> f64 {
    if i == 0 {
        return 0.0;
    }
    let start = i.saturating_sub(window);
    let mut dist = 0.0;
    for k in start..i {
        dist += distance(obs[k].x, obs[k].y, obs[k + 1].x, obs[k + 1].y, geographic);
    }
    let dt = obs[i].t - obs[start].t;
    if dt > 0.0 {
        dist / dt
    } else {
        0.0
    }
}

/// True when, looking back `idle_duration` seconds, the mover stayed within
/// `idle_distance` of point `i` (i.e. barely moved).
fn is_idle(
    obs: &[Obs],
    i: usize,
    idle_distance: f64,
    idle_duration: f64,
    geographic: bool,
) -> bool {
    let now = obs[i];
    let mut k = i;
    let mut spanned = false;
    while k > 0 && now.t - obs[k - 1].t <= idle_duration {
        k -= 1;
        spanned = true;
        if distance(now.x, now.y, obs[k].x, obs[k].y, geographic) > idle_distance {
            return false;
        }
    }
    // Need at least one earlier point within the look-back window to judge.
    spanned
}

/// Initial bearing (degrees clockwise from north) from point 0 to point 1.
fn bearing_deg(x0: f64, y0: f64, x1: f64, y1: f64, geographic: bool) -> f64 {
    let deg = if geographic {
        let (lat0, lat1) = (y0.to_radians(), y1.to_radians());
        let dlon = (x1 - x0).to_radians();
        let yy = dlon.sin() * lat1.cos();
        let xx = lat0.cos() * lat1.sin() - lat0.sin() * lat1.cos() * dlon.cos();
        yy.atan2(xx).to_degrees()
    } else {
        // Clockwise from north (+y): atan2(dx, dy).
        (x1 - x0).atan2(y1 - y0).to_degrees()
    };
    (deg + 360.0) % 360.0
}

fn distance(x0: f64, y0: f64, x1: f64, y1: f64, geographic: bool) -> f64 {
    if geographic {
        haversine(y0, x0, y1, x1)
    } else {
        (x1 - x0).hypot(y1 - y0)
    }
}

fn haversine(lat0: f64, lon0: f64, lat1: f64, lon1: f64) -> f64 {
    const R: f64 = 6_371_000.0;
    let (p0, p1) = (lat0.to_radians(), lat1.to_radians());
    let dphi = (lat1 - lat0).to_radians();
    let dlmb = (lon1 - lon0).to_radians();
    let a = (dphi / 2.0).sin().powi(2) + p0.cos() * p1.cos() * (dlmb / 2.0).sin().powi(2);
    2.0 * R * a.sqrt().asin()
}

fn parse_time_value(fv: &FieldValue) -> Option<f64> {
    if let Some(n) = fv.as_f64() {
        return Some(n);
    }
    fv.as_str().and_then(parse_iso8601_seconds)
}

fn parse_iso8601_seconds(s: &str) -> Option<f64> {
    let s = s.trim();
    if s.len() < 10 {
        return None;
    }
    let b = s.as_bytes();
    let year: i64 = s.get(0..4)?.parse().ok()?;
    if b[4] != b'-' {
        return None;
    }
    let month: i64 = s.get(5..7)?.parse().ok()?;
    let day: i64 = s.get(8..10)?.parse().ok()?;
    let (mut hh, mut mm, mut ss) = (0i64, 0i64, 0i64);
    if s.len() >= 19 && (b[10] == b'T' || b[10] == b' ') {
        hh = s.get(11..13)?.parse().ok()?;
        mm = s.get(14..16)?.parse().ok()?;
        ss = s.get(17..19)?.parse().ok()?;
    }
    Some((days_from_civil(year, month, day) * 86400 + hh * 3600 + mm * 60 + ss) as f64)
}

fn days_from_civil(y: i64, m: i64, d: i64) -> i64 {
    let y = if m <= 2 { y - 1 } else { y };
    let era = if y >= 0 { y } else { y - 399 } / 400;
    let yoe = y - era * 400;
    let doy = (153 * (if m > 2 { m - 3 } else { m + 9 }) + 2) / 5 + d - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146097 + doe - 719468
}

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        _ => None,
    }
}

fn value_string(fv: &FieldValue) -> String {
    if let Some(i) = fv.as_i64() {
        i.to_string()
    } else if let Some(f) = fv.as_f64() {
        format!("{f}")
    } else {
        fv.as_str().unwrap_or("").to_string()
    }
}

struct Params {
    track_field: String,
    time_field: String,
    window: usize,
    idle_distance: f64,
    idle_duration: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let track_field = require_str(args, "track_field")?.to_string();
    let time_field = require_str(args, "time_field")?.to_string();
    let window = opt_usize(args, "window")?.unwrap_or(1).max(1);
    let idle_distance = opt_f64(args, "idle_distance")?.unwrap_or(0.0).max(0.0);
    let idle_duration = opt_f64(args, "idle_duration")?.unwrap_or(60.0).max(0.0);
    Ok(Params {
        track_field,
        time_field,
        window,
        idle_distance,
        idle_duration,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn opt_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
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

fn opt_usize(args: &ToolArgs, key: &str) -> Result<Option<usize>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64().map(|v| v as usize)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be an integer"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be an integer"
        ))),
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

    fn track_layer(rows: &[(&str, f64, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("id", FieldType::Text));
        l.add_field(FieldDef::new("t", FieldType::Float));
        for (id, t, x, y) in rows {
            l.add_feature(
                Some(Geometry::point(*x, *y)),
                &[("id", (*id).into()), ("t", (*t).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> Layer {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CalculateMotionStatisticsTool.run(&args, &ctx()).unwrap();
        load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    fn get(layer: &Layer, feat: usize, name: &str) -> f64 {
        let idx = layer.schema.field_index(name).unwrap();
        layer.features[feat].attributes[idx].as_f64().unwrap()
    }

    /// Speed, distance and cumulative distance are computed correctly along a
    /// straight constant-time track.
    #[test]
    fn computes_speed_and_distance() {
        let rows = [
            ("a", 0.0, 0.0, 0.0),
            ("a", 10.0, 30.0, 40.0), // +50 units in 10 s -> speed 5
            ("a", 20.0, 60.0, 80.0), // +50 units in 10 s -> speed 5
        ];
        let l = run(json!({ "input": track_layer(&rows), "track_field": "id", "time_field": "t" }));
        assert_eq!(l.features.len(), 3);
        // Row order is time-sorted within the track.
        assert!(
            (get(&l, 0, "speed") - 0.0).abs() < 1e-9,
            "first point has no speed"
        );
        assert!((get(&l, 1, "dist") - 50.0).abs() < 1e-6);
        assert!((get(&l, 1, "speed") - 5.0).abs() < 1e-6);
        assert!((get(&l, 2, "cum_dist") - 100.0).abs() < 1e-6);
        assert!((get(&l, 2, "elapsed") - 20.0).abs() < 1e-6);
    }

    /// Acceleration is the change in speed over dt.
    #[test]
    fn computes_acceleration() {
        let rows = [
            ("a", 0.0, 0.0, 0.0),
            ("a", 1.0, 1.0, 0.0), // speed 1
            ("a", 2.0, 4.0, 0.0), // speed 3 -> accel (3-1)/1 = 2
        ];
        let l = run(json!({ "input": track_layer(&rows), "track_field": "id", "time_field": "t" }));
        assert!((get(&l, 2, "accel") - 2.0).abs() < 1e-6);
    }

    /// Bearing due east in a projected CRS is 90 degrees.
    #[test]
    fn computes_bearing() {
        let rows = [("a", 0.0, 0.0, 0.0), ("a", 1.0, 10.0, 0.0)];
        let l = run(json!({ "input": track_layer(&rows), "track_field": "id", "time_field": "t" }));
        assert!((get(&l, 1, "bearing") - 90.0).abs() < 1e-6, "east = 90 deg");
    }

    /// A stationary mover is flagged idle once the look-back window is covered.
    #[test]
    fn flags_idle() {
        let rows = [
            ("a", 0.0, 0.0, 0.0),
            ("a", 30.0, 0.5, 0.5),
            ("a", 60.0, 0.5, 0.0), // barely moved over 60 s
        ];
        let l = run(json!({
            "input": track_layer(&rows), "track_field": "id", "time_field": "t",
            "idle_distance": 5.0, "idle_duration": 60.0,
        }));
        let idle_idx = l.schema.field_index("idle").unwrap();
        assert_eq!(
            l.features[2].attributes[idle_idx].as_i64(),
            Some(1),
            "stationary mover should be idle"
        );
    }

    /// Original attributes are preserved alongside the new motion fields.
    #[test]
    fn preserves_input_attributes() {
        let rows = [("a", 0.0, 0.0, 0.0), ("a", 1.0, 1.0, 0.0)];
        let l = run(json!({ "input": track_layer(&rows), "track_field": "id", "time_field": "t" }));
        assert!(l.schema.field_index("id").is_some());
        assert!(l.schema.field_index("t").is_some());
        assert!(l.schema.field_index("speed").is_some());
        assert_eq!(
            l.features[0].attributes[l.schema.field_index("id").unwrap()].as_str(),
            Some("a")
        );
    }

    #[test]
    fn rejects_missing_required() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CalculateMotionStatisticsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "track_field": "id" })).is_err());
        assert!(
            bad(json!({ "input": "a.geojson", "track_field": "id", "time_field": "t" })).is_ok()
        );
    }
}
