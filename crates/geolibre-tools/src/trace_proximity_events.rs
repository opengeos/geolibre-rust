//! GeoLibre tool: contact tracing / proximity events between moving tracks.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Trace Proximity Events* (GeoAnalytics
//! Desktop). `reconstruct_tracks` builds the tracks; nothing analyzed the
//! interactions *between* them. This finds every interval where two movers were
//! within `search_distance` of each other for at least `min_duration`
//! (proximity events — convoy / meeting / contact detection), and optionally
//! traces the transitive downstream spread from a set of seed `entities`
//! (contact tracing) with a generation number.
//!
//! Points are grouped into tracks by `track_field` and time-sorted (as in
//! `reconstruct_tracks`). For each pair of tracks, over the union of their
//! timestamps, both positions are linearly interpolated; on each inter-sample
//! interval the squared separation is a quadratic in time, solved exactly for
//! the sub-interval within `search_distance`. Touching sub-intervals merge into
//! maximal proximity runs. Output is one connector `LineString` per event with
//! the two track ids, start/end time, duration, and minimum distance.
//!
//! With `entities` set, a temporal breadth-first spread runs over the event
//! graph: a seed is "infected" from the start; an event transmits to the other
//! track if it ends at or after the source's infection time, infecting it at
//! `max(event start, source time)` one generation later, up to `depth`. Each
//! event then also carries the generation of each endpoint (`gen_a`/`gen_b`, −1
//! if never reached).

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Times closer than this (map-unit-seconds) are treated as touching.
const T_EPS: f64 = 1e-6;

pub struct TraceProximityEventsTool;

impl Tool for TraceProximityEventsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "trace_proximity_events",
            display_name: "Trace Proximity Events",
            summary: "Find intervals where two moving tracks were within a distance of each other for at least a minimum duration (proximity events), and optionally trace transitive downstream contacts from seed entities, like ArcGIS Trace Proximity Events.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Timestamped point layer (use a projected CRS; distances are in map units).",
                    required: true,
                },
                ToolParamSpec {
                    name: "track_field",
                    description: "Field identifying each moving entity / track.",
                    required: true,
                },
                ToolParamSpec {
                    name: "time_field",
                    description: "Timestamp field (numeric seconds or ISO-8601).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output proximity-event line layer. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "search_distance",
                    description: "Maximum separation for two tracks to be 'in proximity' (map units).",
                    required: true,
                },
                ToolParamSpec {
                    name: "min_duration",
                    description: "Minimum time in proximity to count as an event: seconds or a duration like '5m', '1h' (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "entities",
                    description: "Optional comma-separated seed track ids to trace downstream contacts from (contact tracing).",
                    required: false,
                },
                ToolParamSpec {
                    name: "depth",
                    description: "Maximum tracing generations from the seeds (default unlimited).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "track_field")?;
        require_str(args, "time_field")?;
        if args.get("search_distance").is_none() {
            return Err(ToolError::Validation(
                "missing required parameter 'search_distance'".to_string(),
            ));
        }
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

        // Group into tracks (id -> time-sorted samples).
        let mut tracks: BTreeMap<String, Vec<Sample>> = BTreeMap::new();
        for feature in layer.iter() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let Some((x, y)) = point_xy(geom) else {
                continue;
            };
            let Some(t) = feature.attributes.get(time_idx).and_then(parse_time_value) else {
                continue;
            };
            let id = feature
                .attributes
                .get(track_idx)
                .map(value_string)
                .unwrap_or_default();
            tracks.entry(id).or_default().push(Sample { x, y, t });
        }
        for s in tracks.values_mut() {
            s.sort_by(|a, b| a.t.total_cmp(&b.t));
            s.dedup_by(|a, b| a.t == b.t);
        }
        let ids: Vec<String> = tracks.keys().cloned().collect();
        let series: Vec<&Vec<Sample>> = ids.iter().map(|k| &tracks[k]).collect();
        if ids.len() < 2 {
            return Err(ToolError::Execution(
                "need at least two tracks to find proximity events".to_string(),
            ));
        }

        ctx.progress.info(&format!(
            "scanning {} track pair(s)",
            ids.len() * (ids.len() - 1) / 2
        ));

        // ── Pairwise proximity events ────────────────────────────────────────────
        let d2 = prm.search_distance * prm.search_distance;
        let mut events: Vec<Event> = Vec::new();
        for a in 0..ids.len() {
            for b in (a + 1)..ids.len() {
                for ev in pair_events(series[a], series[b], d2, prm.min_duration) {
                    events.push(Event {
                        a,
                        b,
                        start: ev.start,
                        end: ev.end,
                        min_dist: ev.min_dist,
                        ax: ev.ax,
                        ay: ev.ay,
                        bx: ev.bx,
                        by: ev.by,
                    });
                }
            }
        }

        // ── Optional contact tracing (temporal BFS over the event graph) ─────────
        let generations: Option<Vec<i64>> = prm.entities.as_ref().map(|seeds| {
            let seed_set: Vec<usize> = seeds
                .iter()
                .filter_map(|s| ids.iter().position(|x| x == s))
                .collect();
            trace(&events, ids.len(), &seed_set, prm.depth)
        });

        // ── Build output ─────────────────────────────────────────────────────────
        let mut out = Layer::new("proximity_events").with_geom_type(GeometryType::LineString);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("track_a", FieldType::Text));
        out.add_field(FieldDef::new("track_b", FieldType::Text));
        out.add_field(FieldDef::new("start_t", FieldType::Float));
        out.add_field(FieldDef::new("end_t", FieldType::Float));
        out.add_field(FieldDef::new("duration", FieldType::Float));
        out.add_field(FieldDef::new("min_dist", FieldType::Float));
        if generations.is_some() {
            out.add_field(FieldDef::new("gen_a", FieldType::Integer));
            out.add_field(FieldDef::new("gen_b", FieldType::Integer));
        }

        for ev in &events {
            let geom = Geometry::LineString(vec![Coord::xy(ev.ax, ev.ay), Coord::xy(ev.bx, ev.by)]);
            let mut attrs: Vec<(&str, FieldValue)> = vec![
                ("track_a", FieldValue::Text(ids[ev.a].clone())),
                ("track_b", FieldValue::Text(ids[ev.b].clone())),
                ("start_t", FieldValue::Float(ev.start)),
                ("end_t", FieldValue::Float(ev.end)),
                ("duration", FieldValue::Float(ev.end - ev.start)),
                ("min_dist", FieldValue::Float(ev.min_dist)),
            ];
            if let Some(gen) = &generations {
                attrs.push(("gen_a", FieldValue::Integer(gen[ev.a])));
                attrs.push(("gen_b", FieldValue::Integer(gen[ev.b])));
            }
            out.add_feature(Some(geom), &attrs)
                .map_err(|e| ToolError::Execution(format!("failed writing event: {e}")))?;
        }

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("track_count".to_string(), json!(ids.len()));
        outputs.insert("event_count".to_string(), json!(events.len()));
        if let Some(gen) = &generations {
            let reached = gen.iter().filter(|g| **g >= 0).count();
            outputs.insert("traced_reached".to_string(), json!(reached));
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── Pairwise proximity ─────────────────────────────────────────────────────────

struct Sample {
    x: f64,
    y: f64,
    t: f64,
}

struct Event {
    a: usize,
    b: usize,
    start: f64,
    end: f64,
    min_dist: f64,
    ax: f64,
    ay: f64,
    bx: f64,
    by: f64,
}

struct PairEvent {
    start: f64,
    end: f64,
    min_dist: f64,
    ax: f64,
    ay: f64,
    bx: f64,
    by: f64,
}

/// Interpolate a track's position at time `t` (None outside its span).
fn interp(s: &[Sample], t: f64) -> Option<(f64, f64)> {
    if s.is_empty() || t < s[0].t - T_EPS || t > s[s.len() - 1].t + T_EPS {
        return None;
    }
    // Binary search for the bracketing pair.
    let mut lo = 0usize;
    let mut hi = s.len() - 1;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if s[mid].t <= t {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let (p, q) = (&s[lo], &s[hi]);
    if (q.t - p.t).abs() < T_EPS {
        return Some((p.x, p.y));
    }
    let f = ((t - p.t) / (q.t - p.t)).clamp(0.0, 1.0);
    Some((p.x + (q.x - p.x) * f, p.y + (q.y - p.y) * f))
}

/// All proximity events between two tracks.
fn pair_events(a: &[Sample], b: &[Sample], d2: f64, min_dur: f64) -> Vec<PairEvent> {
    // Common timeline: union of both tracks' timestamps within their overlap.
    let t_lo = a[0].t.max(b[0].t);
    let t_hi = a[a.len() - 1].t.min(b[b.len() - 1].t);
    if t_hi <= t_lo {
        return Vec::new();
    }
    let mut times: Vec<f64> = Vec::new();
    for s in a.iter().chain(b.iter()) {
        if s.t >= t_lo - T_EPS && s.t <= t_hi + T_EPS {
            times.push(s.t);
        }
    }
    times.push(t_lo);
    times.push(t_hi);
    times.sort_by(f64::total_cmp);
    times.dedup_by(|x, y| (*x - *y).abs() < T_EPS);
    if times.len() < 2 {
        return Vec::new();
    }

    // Within-threshold pieces per inter-sample segment, then merge touching ones.
    let mut pieces: Vec<(f64, f64)> = Vec::new(); // (t_start, t_end)
    for w in times.windows(2) {
        let (t0, t1) = (w[0], w[1]);
        let (Some(a0), Some(b0)) = (interp(a, t0), interp(b, t0)) else {
            continue;
        };
        let (Some(a1), Some(b1)) = (interp(a, t1), interp(b, t1)) else {
            continue;
        };
        // Relative position R(s) = R0 + (R1-R0) s, s in [0,1].
        let r0 = (a0.0 - b0.0, a0.1 - b0.1);
        let r1 = (a1.0 - b1.0, a1.1 - b1.1);
        let dr = (r1.0 - r0.0, r1.1 - r0.1);
        // |R|^2 = A s^2 + B s + C <= d2
        let aa = dr.0 * dr.0 + dr.1 * dr.1;
        let bb = 2.0 * (r0.0 * dr.0 + r0.1 * dr.1);
        let cc = r0.0 * r0.0 + r0.1 * r0.1 - d2;
        let (mut s_lo, mut s_hi) = (f64::NAN, f64::NAN);
        if aa.abs() < 1e-15 {
            // Constant separation over the segment.
            if cc <= 0.0 {
                s_lo = 0.0;
                s_hi = 1.0;
            }
        } else {
            let disc = bb * bb - 4.0 * aa * cc;
            if disc >= 0.0 {
                let sq = disc.sqrt();
                let mut lo = (-bb - sq) / (2.0 * aa);
                let mut hi = (-bb + sq) / (2.0 * aa);
                if lo > hi {
                    std::mem::swap(&mut lo, &mut hi);
                }
                s_lo = lo.clamp(0.0, 1.0);
                s_hi = hi.clamp(0.0, 1.0);
                if s_lo >= s_hi {
                    // No within-portion inside [0,1].
                    s_lo = f64::NAN;
                }
            }
        }
        if s_lo.is_finite() && s_hi.is_finite() && s_hi > s_lo {
            let ts = t0 + s_lo * (t1 - t0);
            let te = t0 + s_hi * (t1 - t0);
            pieces.push((ts, te));
        }
    }
    if pieces.is_empty() {
        return Vec::new();
    }

    // Merge touching/overlapping pieces into maximal runs.
    let mut runs: Vec<(f64, f64)> = Vec::new();
    for (ts, te) in pieces {
        match runs.last_mut() {
            Some((_, prev_end)) if ts <= *prev_end + T_EPS => {
                if te > *prev_end {
                    *prev_end = te;
                }
            }
            _ => runs.push((ts, te)),
        }
    }

    // Each run of sufficient duration is an event; sample its minimum distance.
    let mut out = Vec::new();
    for (ts, te) in runs {
        if te - ts < min_dur {
            continue;
        }
        // Minimum distance and its position over the run (sample finely).
        let steps = 32;
        let mut best_d = f64::INFINITY;
        let mut best = (0.0, 0.0, 0.0, 0.0);
        for k in 0..=steps {
            let t = ts + (te - ts) * k as f64 / steps as f64;
            if let (Some(pa), Some(pb)) = (interp(a, t), interp(b, t)) {
                let dd = (pa.0 - pb.0).powi(2) + (pa.1 - pb.1).powi(2);
                if dd < best_d {
                    best_d = dd;
                    best = (pa.0, pa.1, pb.0, pb.1);
                }
            }
        }
        out.push(PairEvent {
            start: ts,
            end: te,
            min_dist: best_d.sqrt(),
            ax: best.0,
            ay: best.1,
            bx: best.2,
            by: best.3,
        });
    }
    out
}

// ── Contact tracing (temporal BFS) ─────────────────────────────────────────────

/// Infection generation per track: seeds are 0; an event transmits from an
/// infected source (infected at time `t_src`) to the other track if it ends at
/// or after `t_src`, infecting it at `max(event.start, t_src)` one generation
/// later, bounded by `depth`. Returns −1 for tracks never reached.
fn trace(events: &[Event], n: usize, seeds: &[usize], depth: Option<usize>) -> Vec<i64> {
    let mut gen = vec![-1i64; n];
    let mut inf_time = vec![f64::INFINITY; n];
    for &s in seeds {
        gen[s] = 0;
        inf_time[s] = f64::NEG_INFINITY;
    }
    let max_gen = depth.map(|d| d as i64).unwrap_or(i64::MAX);
    // Relax until stable (bounded by n passes).
    for _ in 0..=n {
        let mut changed = false;
        for ev in events {
            for (src, dst) in [(ev.a, ev.b), (ev.b, ev.a)] {
                if gen[src] < 0 || gen[src] >= max_gen {
                    continue;
                }
                // Contact must be able to happen after the source was infected.
                if ev.end + T_EPS < inf_time[src] {
                    continue;
                }
                let t_new = ev.start.max(inf_time[src]);
                let g_new = gen[src] + 1;
                // Accept if dst not infected, or reached earlier / at a lower gen.
                if gen[dst] < 0 || t_new < inf_time[dst] - T_EPS {
                    gen[dst] = g_new;
                    inf_time[dst] = t_new;
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
    gen
}

// ── Shared helpers (as in reconstruct_tracks) ──────────────────────────────────

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

fn parse_duration(s: &str) -> Result<f64, ToolError> {
    let s = s.trim();
    if let Ok(v) = s.parse::<f64>() {
        if v >= 0.0 && v.is_finite() {
            return Ok(v);
        }
        return Err(ToolError::Validation(
            "'min_duration' must be >= 0".to_string(),
        ));
    }
    let (num, unit) = s.split_at(s.len().saturating_sub(1));
    let value: f64 = num.trim().parse().map_err(|_| {
        ToolError::Validation(format!("could not parse 'min_duration' value in '{s}'"))
    })?;
    let seconds = match unit {
        "s" => 1.0,
        "m" => 60.0,
        "h" => 3600.0,
        "d" => 86400.0,
        "w" => 604800.0,
        other => {
            return Err(ToolError::Validation(format!(
                "unknown 'min_duration' unit '{other}' (use s/m/h/d/w or a plain number)"
            )))
        }
    };
    Ok(value * seconds)
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    track_field: String,
    time_field: String,
    search_distance: f64,
    min_duration: f64,
    entities: Option<Vec<String>>,
    depth: Option<usize>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let track_field = require_str(args, "track_field")?.to_string();
    let time_field = require_str(args, "time_field")?.to_string();
    let search_distance = match args.get("search_distance") {
        Some(Value::Number(n)) => n.as_f64().unwrap_or(0.0),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'search_distance' must be a number".into()))?,
        _ => {
            return Err(ToolError::Validation(
                "missing required numeric parameter 'search_distance'".to_string(),
            ))
        }
    };
    if search_distance.is_nan() || search_distance <= 0.0 {
        return Err(ToolError::Validation(
            "'search_distance' must be positive".to_string(),
        ));
    }
    let min_duration = match args.get("min_duration") {
        None | Some(Value::Null) => 0.0,
        Some(Value::Number(n)) => n.as_f64().unwrap_or(0.0).max(0.0),
        Some(Value::String(s)) if s.trim().is_empty() => 0.0,
        Some(Value::String(s)) => parse_duration(s)?,
        Some(_) => {
            return Err(ToolError::Validation(
                "'min_duration' must be a number".into(),
            ))
        }
    };
    let entities = parse_optional_str(args, "entities")?.map(|s| {
        s.split(',')
            .map(str::trim)
            .filter(|x| !x.is_empty())
            .map(str::to_string)
            .collect::<Vec<_>>()
    });
    let depth = match args.get("depth") {
        None | Some(Value::Null) => None,
        Some(Value::Number(n)) => n.as_u64().map(|v| v as usize),
        Some(Value::String(s)) if s.trim().is_empty() => None,
        Some(Value::String(s)) => Some(
            s.trim()
                .parse::<usize>()
                .map_err(|_| ToolError::Validation("'depth' must be an integer".into()))?,
        ),
        Some(_) => return Err(ToolError::Validation("'depth' must be a number".into())),
    };
    Ok(Params {
        track_field,
        time_field,
        search_distance,
        min_duration,
        entities,
        depth,
    })
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
    use wbvector::memory_store;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// samples: (track, x, y, t)
    fn layer_of(pts: &[(&str, f64, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("id", FieldType::Text));
        l.add_field(FieldDef::new("t", FieldType::Float));
        for (id, x, y, t) in pts {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(*x, *y))),
                &[("id", (*id).into()), ("t", (*t).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = TraceProximityEventsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// Two tracks that pass close together produce one proximity event.
    #[test]
    fn detects_a_close_pass() {
        // A moves along y=0; B moves along y=2; they are 2 apart the whole time
        // between t=0..10. search_distance 3 -> in proximity throughout.
        let pts = [
            ("A", 0.0, 0.0, 0.0),
            ("A", 10.0, 0.0, 10.0),
            ("B", 0.0, 2.0, 0.0),
            ("B", 10.0, 2.0, 10.0),
        ];
        let input = layer_of(&pts);
        let (out, layer) = run(json!({
            "input": input, "track_field": "id", "time_field": "t", "search_distance": 3.0
        }));
        assert_eq!(out.outputs["event_count"], json!(1));
        let di = layer.schema.field_index("min_dist").unwrap();
        let du = layer.schema.field_index("duration").unwrap();
        let f = layer.iter().next().unwrap();
        assert!((f.attributes[di].as_f64().unwrap() - 2.0).abs() < 1e-6);
        assert!((f.attributes[du].as_f64().unwrap() - 10.0).abs() < 1e-6);
    }

    /// A wide separation yields no events.
    #[test]
    fn far_apart_no_events() {
        let pts = [
            ("A", 0.0, 0.0, 0.0),
            ("A", 10.0, 0.0, 10.0),
            ("B", 0.0, 100.0, 0.0),
            ("B", 10.0, 100.0, 10.0),
        ];
        let input = layer_of(&pts);
        let (out, _l) = run(json!({
            "input": input, "track_field": "id", "time_field": "t", "search_distance": 5.0
        }));
        assert_eq!(out.outputs["event_count"], json!(0));
    }

    /// min_duration filters out brief crossings.
    #[test]
    fn min_duration_filters_brief() {
        // A and B cross paths: within 3 units only briefly around t=5.
        let pts = [
            ("A", 0.0, 0.0, 0.0),
            ("A", 10.0, 0.0, 10.0),
            ("B", 5.0, -10.0, 0.0),
            ("B", 5.0, 10.0, 10.0),
        ];
        let input = layer_of(&pts);
        // Brief pass -> excluded by a long min_duration.
        let (out, _l) = run(json!({
            "input": input, "track_field": "id", "time_field": "t",
            "search_distance": 3.0, "min_duration": 100.0
        }));
        assert_eq!(out.outputs["event_count"], json!(0), "brief pass filtered");
    }

    /// Contact tracing: A meets B (early), B meets C (later) -> C is generation 2.
    #[test]
    fn traces_transitive_contacts() {
        let pts = [
            // A and B together near t=0..2.
            ("A", 0.0, 0.0, 0.0),
            ("A", 0.0, 0.0, 2.0),
            ("A", 0.0, 0.0, 20.0),
            ("B", 1.0, 0.0, 0.0),
            ("B", 1.0, 0.0, 2.0),
            ("B", 100.0, 0.0, 10.0),
            ("B", 100.0, 0.0, 12.0),
            // C meets B later near t=10..12 at (100,0).
            ("C", 101.0, 0.0, 10.0),
            ("C", 101.0, 0.0, 12.0),
        ];
        let input = layer_of(&pts);
        let (out, layer) = run(json!({
            "input": input, "track_field": "id", "time_field": "t",
            "search_distance": 3.0, "entities": "A"
        }));
        assert_eq!(out.outputs["traced_reached"], json!(3), "A, B, C reached");
        // The B-C event should carry gen_b (C) = 2.
        let ga = layer.schema.field_index("gen_a").unwrap();
        let gb = layer.schema.field_index("gen_b").unwrap();
        let ta = layer.schema.field_index("track_a").unwrap();
        let tb = layer.schema.field_index("track_b").unwrap();
        let mut c_gen = -1;
        for f in layer.iter() {
            for (tf, gf) in [(ta, ga), (tb, gb)] {
                if f.attributes[tf].as_str() == Some("C") {
                    c_gen = f.attributes[gf].as_i64().unwrap();
                }
            }
        }
        assert_eq!(c_gen, 2, "C is second generation");
    }

    #[test]
    fn rejects_missing_search_distance() {
        let input = layer_of(&[("A", 0.0, 0.0, 0.0), ("B", 1.0, 0.0, 0.0)]);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input, "track_field": "id", "time_field": "t"
        }))
        .unwrap();
        assert!(TraceProximityEventsTool.validate(&args).is_err());
    }
}
