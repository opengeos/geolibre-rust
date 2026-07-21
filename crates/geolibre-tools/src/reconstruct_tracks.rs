//! GeoLibre tool: reconstruct movement tracks from timestamped points.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Reconstruct Tracks* and *Find Dwell
//! Locations* (GeoAnalytics). There are zero trajectory/movement tools among the
//! ~791 bundled IDs, yet movement data (GPS, AIS, wildlife telemetry) is one of
//! the most common modern geospatial inputs. Reconstructed tracks feed straight
//! into `vector_to_pmtiles` / H3 binning for web visualization.
//!
//! Points are grouped by `track_field`, sorted by `time_field`, and connected
//! into polylines. A new track segment starts when the gap to the previous point
//! exceeds `time_gap` (seconds) or `distance_gap` (metres for a geographic CRS,
//! CRS units otherwise). Each output track carries `track_id`, point count,
//! start/end time, duration, length, and mean/max speed.
//!
//! With a `dwells` output path the tool also finds **dwell locations**: maximal
//! runs of consecutive points that stay within `dwell_distance` of the run's
//! anchor for at least `dwell_min_duration` seconds — where the mover paused.
//! Each dwell is emitted as a point at the run's centroid with its time span.
//! `time_field` is a numeric field (seconds) or an ISO-8601 timestamp.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct ReconstructTracksTool;

impl Tool for ReconstructTracksTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "reconstruct_tracks",
            display_name: "Reconstruct Tracks",
            summary: "Turn timestamped points into movement track polylines (grouped by id, sorted by time, split on time/distance gaps) with per-track stats, and optionally find dwell locations where the mover paused — like ArcGIS Reconstruct Tracks / Find Dwell Locations.",
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
                    description: "Field identifying each track/mover (e.g. vehicle or animal id).",
                    required: true,
                },
                ToolParamSpec {
                    name: "time_field",
                    description: "Field holding each point's time: numeric seconds or an ISO-8601 timestamp.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output line vector of tracks (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "time_gap",
                    description: "Start a new track segment when the time gap between points exceeds this many seconds. Default: no time split.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distance_gap",
                    description: "Start a new track segment when the distance between points exceeds this (metres for geographic CRS, CRS units otherwise). Default: no distance split.",
                    required: false,
                },
                ToolParamSpec {
                    name: "dwells",
                    description: "Optional output point path for dwell locations (where the mover paused).",
                    required: false,
                },
                ToolParamSpec {
                    name: "dwell_distance",
                    description: "Max radius for a dwell, in the distance units above. Default 50.",
                    required: false,
                },
                ToolParamSpec {
                    name: "dwell_min_duration",
                    description: "Minimum time (seconds) a mover must stay within dwell_distance to count as a dwell. Default 300.",
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
        let dwells_path = parse_optional_str(args, "dwells")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let track_idx = layer.schema.field_index(&prm.track_field).ok_or_else(|| {
            ToolError::Validation(format!("track_field '{}' not found", prm.track_field))
        })?;
        let time_idx = layer.schema.field_index(&prm.time_field).ok_or_else(|| {
            ToolError::Validation(format!("time_field '{}' not found", prm.time_field))
        })?;

        // Geographic distances (metres) when the CRS is lon/lat.
        let geographic = layer.crs_epsg().map(|e| e == 4326).unwrap_or(true);

        // Collect observations grouped by track.
        let mut tracks: BTreeMap<String, Vec<Obs>> = BTreeMap::new();
        let mut skipped = 0usize;
        for feature in layer.iter() {
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
            tracks.entry(id).or_default().push(Obs { x, y, t });
        }

        ctx.progress
            .info(&format!("{} track(s) to reconstruct", tracks.len()));

        // Build the track line layer.
        let mut out = Layer::new("tracks").with_geom_type(GeometryType::LineString);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("track_id", FieldType::Text));
        out.add_field(FieldDef::new("segment", FieldType::Integer));
        out.add_field(FieldDef::new("n_points", FieldType::Integer));
        out.add_field(FieldDef::new("start_t", FieldType::Float));
        out.add_field(FieldDef::new("end_t", FieldType::Float));
        out.add_field(FieldDef::new("duration", FieldType::Float));
        out.add_field(FieldDef::new("length", FieldType::Float));
        out.add_field(FieldDef::new("mean_speed", FieldType::Float));
        out.add_field(FieldDef::new("max_speed", FieldType::Float));

        // Optional dwell layer.
        let mut dwell_layer = Layer::new("dwells").with_geom_type(GeometryType::Point);
        if let Some(epsg) = layer.crs_epsg() {
            dwell_layer = dwell_layer.with_crs_epsg(epsg);
        }
        dwell_layer.add_field(FieldDef::new("track_id", FieldType::Text));
        dwell_layer.add_field(FieldDef::new("start_t", FieldType::Float));
        dwell_layer.add_field(FieldDef::new("end_t", FieldType::Float));
        dwell_layer.add_field(FieldDef::new("duration", FieldType::Float));
        dwell_layer.add_field(FieldDef::new("n_points", FieldType::Integer));

        let mut track_count = 0usize;
        let mut dwell_count = 0usize;
        for (id, mut obs) in tracks {
            obs.sort_by(|a, b| a.t.total_cmp(&b.t));
            // Split into segments on time/distance gaps.
            let segments = split_segments(&obs, &prm, geographic);
            for (seg_no, seg) in segments.iter().enumerate() {
                if seg.len() < 2 {
                    continue;
                }
                let stats = track_stats(seg, geographic);
                let coords: Vec<Coord> = seg.iter().map(|o| Coord::xy(o.x, o.y)).collect();
                out.push(Feature {
                    fid: 0,
                    geometry: Some(Geometry::line_string(coords)),
                    attributes: vec![
                        FieldValue::Text(id.clone()),
                        FieldValue::Integer(seg_no as i64),
                        FieldValue::Integer(seg.len() as i64),
                        FieldValue::Float(stats.start_t),
                        FieldValue::Float(stats.end_t),
                        FieldValue::Float(stats.duration),
                        FieldValue::Float(stats.length),
                        FieldValue::Float(stats.mean_speed),
                        FieldValue::Float(stats.max_speed),
                    ],
                });
                track_count += 1;

                if dwells_path.is_some() {
                    for dw in find_dwells(seg, &prm, geographic) {
                        dwell_layer.push(Feature {
                            fid: 0,
                            geometry: Some(Geometry::point(dw.cx, dw.cy)),
                            attributes: vec![
                                FieldValue::Text(id.clone()),
                                FieldValue::Float(dw.start_t),
                                FieldValue::Float(dw.end_t),
                                FieldValue::Float(dw.end_t - dw.start_t),
                                FieldValue::Integer(dw.n as i64),
                            ],
                        });
                        dwell_count += 1;
                    }
                }
            }
        }

        let out_path = write_or_store_layer(out, output)?;
        let dwells_out = match dwells_path {
            Some(path) => Some(write_or_store_layer(dwell_layer, Some(path))?),
            None => None,
        };

        ctx.progress.info(&format!(
            "{track_count} track segment(s), {dwell_count} dwell(s)"
        ));

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("track_count".to_string(), json!(track_count));
        outputs.insert("dwell_count".to_string(), json!(dwell_count));
        outputs.insert("skipped".to_string(), json!(skipped));
        if let Some(d) = dwells_out {
            outputs.insert("dwells".to_string(), json!(d));
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── Track building ───────────────────────────────────────────────────────────

#[derive(Clone, Copy)]
struct Obs {
    x: f64,
    y: f64,
    t: f64,
}

/// Splits a time-sorted observation list into segments at time/distance gaps.
fn split_segments(obs: &[Obs], prm: &Params, geographic: bool) -> Vec<Vec<Obs>> {
    let mut segments = Vec::new();
    let mut cur: Vec<Obs> = Vec::new();
    for (i, &o) in obs.iter().enumerate() {
        if i > 0 {
            let prev = obs[i - 1];
            let dt = o.t - prev.t;
            let dd = distance(prev.x, prev.y, o.x, o.y, geographic);
            let time_break = prm.time_gap.map(|g| dt > g).unwrap_or(false);
            let dist_break = prm.distance_gap.map(|g| dd > g).unwrap_or(false);
            if (time_break || dist_break) && !cur.is_empty() {
                segments.push(std::mem::take(&mut cur));
            }
        }
        cur.push(o);
    }
    if !cur.is_empty() {
        segments.push(cur);
    }
    segments
}

struct TrackStats {
    start_t: f64,
    end_t: f64,
    duration: f64,
    length: f64,
    mean_speed: f64,
    max_speed: f64,
}

fn track_stats(seg: &[Obs], geographic: bool) -> TrackStats {
    let start_t = seg[0].t;
    let end_t = seg[seg.len() - 1].t;
    let duration = end_t - start_t;
    let mut length = 0.0;
    let mut max_speed = 0.0;
    for w in seg.windows(2) {
        let d = distance(w[0].x, w[0].y, w[1].x, w[1].y, geographic);
        length += d;
        let dt = w[1].t - w[0].t;
        if dt > 0.0 {
            max_speed = f64::max(max_speed, d / dt);
        }
    }
    let mean_speed = if duration > 0.0 {
        length / duration
    } else {
        0.0
    };
    TrackStats {
        start_t,
        end_t,
        duration,
        length,
        mean_speed,
        max_speed,
    }
}

struct Dwell {
    cx: f64,
    cy: f64,
    start_t: f64,
    end_t: f64,
    n: usize,
}

/// Finds maximal runs of points staying within `dwell_distance` of the run's
/// anchor for at least `dwell_min_duration` seconds.
fn find_dwells(seg: &[Obs], prm: &Params, geographic: bool) -> Vec<Dwell> {
    let mut dwells = Vec::new();
    let mut i = 0;
    while i < seg.len() {
        let anchor = seg[i];
        let mut j = i + 1;
        while j < seg.len()
            && distance(anchor.x, anchor.y, seg[j].x, seg[j].y, geographic) <= prm.dwell_distance
        {
            j += 1;
        }
        // Run is seg[i..j].
        let run = &seg[i..j];
        if run.len() >= 2 {
            let duration = run[run.len() - 1].t - run[0].t;
            if duration >= prm.dwell_min_duration {
                let n = run.len();
                let cx = run.iter().map(|o| o.x).sum::<f64>() / n as f64;
                let cy = run.iter().map(|o| o.y).sum::<f64>() / n as f64;
                dwells.push(Dwell {
                    cx,
                    cy,
                    start_t: run[0].t,
                    end_t: run[run.len() - 1].t,
                    n,
                });
                i = j; // continue past the dwell
                continue;
            }
        }
        i += 1;
    }
    dwells
}

// ── Distance ─────────────────────────────────────────────────────────────────

/// Distance between two points: haversine metres for a geographic CRS, planar
/// CRS units otherwise.
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

// ── Time / value parsing ─────────────────────────────────────────────────────

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

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    track_field: String,
    time_field: String,
    time_gap: Option<f64>,
    distance_gap: Option<f64>,
    dwell_distance: f64,
    dwell_min_duration: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let track_field = require_str(args, "track_field")?.to_string();
    let time_field = require_str(args, "time_field")?.to_string();
    let time_gap = opt_pos(args, "time_gap")?;
    let distance_gap = opt_pos(args, "distance_gap")?;
    let dwell_distance = opt_pos(args, "dwell_distance")?.unwrap_or(50.0);
    let dwell_min_duration = opt_f64(args, "dwell_min_duration")?
        .unwrap_or(300.0)
        .max(0.0);
    Ok(Params {
        track_field,
        time_field,
        time_gap,
        distance_gap,
        dwell_distance,
        dwell_min_duration,
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

fn opt_pos(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match opt_f64(args, key)? {
        Some(v) if v > 0.0 && v.is_finite() => Ok(Some(v)),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a positive number"
        ))),
        None => Ok(None),
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

    /// Builds a point layer with (track, time, x, y) rows in a projected CRS.
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

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ReconstructTracksTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn field(layer: &Layer, feat: usize, name: &str) -> f64 {
        let idx = layer.schema.field_index(name).unwrap();
        layer.features[feat].attributes[idx].as_f64().unwrap()
    }

    /// Two movers become two tracks; length and duration are correct.
    #[test]
    fn builds_one_track_per_mover() {
        let rows = [
            ("a", 0.0, 0.0, 0.0),
            ("a", 10.0, 30.0, 40.0), // moved 50 units in 10 s
            ("b", 0.0, 100.0, 0.0),
            ("b", 5.0, 100.0, 10.0), // moved 10 units in 5 s
        ];
        let (out, layer) =
            run(json!({ "input": track_layer(&rows), "track_field": "id", "time_field": "t" }));
        assert_eq!(out.outputs["track_count"], json!(2));
        // Find track 'a' (length 50, duration 10, mean speed 5).
        let tid = layer.schema.field_index("track_id").unwrap();
        let a = (0..2)
            .find(|&i| layer.features[i].attributes[tid].as_str() == Some("a"))
            .unwrap();
        assert!((field(&layer, a, "length") - 50.0).abs() < 1e-6);
        assert!((field(&layer, a, "duration") - 10.0).abs() < 1e-6);
        assert!((field(&layer, a, "mean_speed") - 5.0).abs() < 1e-6);
    }

    /// A time gap splits one mover into two segments.
    #[test]
    fn time_gap_splits_segments() {
        let rows = [
            ("a", 0.0, 0.0, 0.0),
            ("a", 10.0, 10.0, 0.0),
            ("a", 1000.0, 20.0, 0.0), // 990 s gap
            ("a", 1010.0, 30.0, 0.0),
        ];
        let (out, _l) = run(json!({
            "input": track_layer(&rows), "track_field": "id", "time_field": "t", "time_gap": 100.0,
        }));
        assert_eq!(
            out.outputs["track_count"],
            json!(2),
            "a 990s gap should split into 2 segments"
        );
    }

    /// A stationary cluster of points over enough time is a dwell.
    #[test]
    fn detects_a_dwell() {
        // Mover sits near (0,0) for 600s (5 points), then leaves.
        let rows = [
            ("a", 0.0, 0.0, 0.0),
            ("a", 150.0, 1.0, 1.0),
            ("a", 300.0, 0.0, 2.0),
            ("a", 450.0, 2.0, 0.0),
            ("a", 600.0, 1.0, 1.0),
            ("a", 700.0, 500.0, 500.0), // leaves
        ];
        let args: ToolArgs = serde_json::from_value(json!({
            "input": track_layer(&rows), "track_field": "id", "time_field": "t",
            "dwells": null, "dwell_distance": 10.0, "dwell_min_duration": 300.0,
        }))
        .unwrap();
        // Run with an in-memory dwells output.
        let mut v = serde_json::to_value(&args).unwrap();
        v["dwells"] = json!("memory://force"); // any non-empty triggers dwell output
                                               // Simpler: call run() directly and read dwell_count.
        let out = ReconstructTracksTool.run(&args, &ctx());
        // dwells path was null -> dwell_count 0 unless we pass a path; assert the
        // detector via a direct call instead.
        let _ = (v, out);
        let obs: Vec<Obs> = rows
            .iter()
            .map(|(_, t, x, y)| Obs {
                x: *x,
                y: *y,
                t: *t,
            })
            .collect();
        let prm = Params {
            track_field: "id".into(),
            time_field: "t".into(),
            time_gap: None,
            distance_gap: None,
            dwell_distance: 10.0,
            dwell_min_duration: 300.0,
        };
        let dwells = find_dwells(&obs[..5], &prm, false);
        assert_eq!(
            dwells.len(),
            1,
            "the stationary cluster should be one dwell"
        );
        assert!(dwells[0].end_t - dwells[0].start_t >= 300.0);
    }

    #[test]
    fn parses_iso_time() {
        assert_eq!(
            parse_iso8601_seconds("2024-06-01T09:00:00.000Z"),
            Some(1_717_232_400.0)
        );
        let a = parse_iso8601_seconds("2024-06-01T09:00:00Z").unwrap();
        let b = parse_iso8601_seconds("2024-06-01T09:05:00Z").unwrap();
        assert_eq!(b - a, 300.0);
    }

    #[test]
    fn rejects_missing_required() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ReconstructTracksTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "track_field": "id" })).is_err());
        assert!(
            bad(json!({ "input": "a.geojson", "track_field": "id", "time_field": "t" })).is_ok()
        );
    }
}
