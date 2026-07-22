//! GeoLibre tool: detect incident episodes along movement tracks.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Detect Incidents* (GeoAnalytics).
//! It completes the movement suite (`reconstruct_tracks`, `calculate_motion_statistics`,
//! `trace_proximity_events`, `find_space_time_matches`): after per-point metrics such
//! as speed or acceleration exist, this tool extracts the *episodes* where a condition
//! holds — speeding intervals, sensor-threshold exceedances, geofence dwell. Nothing
//! comparable is bundled among the ~791 whitebox IDs.
//!
//! Points are grouped by `track_field`, sorted by `time_field`, and a `start_condition`
//! expression over a numeric attribute is evaluated per point. An incident begins at the
//! first point satisfying the start condition. Without an `end_condition` the incident is
//! the maximal run of consecutive points that keep satisfying the start condition. With an
//! `end_condition` the incident stays open (points are "ongoing") until a later point
//! satisfies the end condition (or the track ends). Each incident gets a global id, a
//! per-point status (`start` / `ongoing` / `end`), a within-incident sequence number and
//! the incident's total duration.
//!
//! `mode = points` (default) copies every input point through with the incident fields
//! appended (non-incident points get null ids and an empty status). `mode = segments`
//! emits one polyline per incident that has at least two points, connecting the incident's
//! points in time order, with count / duration / length attributes.
//!
//! The expression grammar is deliberately small: `<field> <op> <number>` where
//! `<op>` is one of `>`, `>=`, `<`, `<=`, `==` (also `=`), `!=`, optionally two clauses
//! joined by `AND` (or `&&`). `time_field` is numeric seconds or an ISO-8601 timestamp.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct DetectIncidentsTool;

impl Tool for DetectIncidentsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "detect_incidents",
            display_name: "Detect Incidents",
            summary: "Flag where a condition starts and ends along each track: group timestamped points by track, sort by time, evaluate a start (and optional end) comparison over a numeric attribute, and emit incident episodes as flagged points or per-episode line segments with ids and durations — like ArcGIS Detect Incidents.",
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
                    name: "start_condition",
                    description: "Condition that opens an incident, e.g. 'speed > 30'. Grammar: '<field> <op> <number>' with op one of > >= < <= == != , optionally two clauses joined by AND.",
                    required: true,
                },
                ToolParamSpec {
                    name: "end_condition",
                    description: "Optional condition that closes an incident (same grammar). If omitted, an incident ends when the start condition stops holding.",
                    required: false,
                },
                ToolParamSpec {
                    name: "mode",
                    description: "Output mode: 'points' (default, every point copied through with incident fields) or 'segments' (one polyline per incident).",
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
        require_str(args, "track_field")?;
        require_str(args, "time_field")?;
        require_str(args, "start_condition")?;
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

        // Resolve condition field indices against the schema.
        let start_cond = prm.start.resolve(&layer)?;
        let end_cond = match &prm.end {
            Some(c) => Some(c.resolve(&layer)?),
            None => None,
        };

        // Collect observations grouped by track (keeping the source feature index).
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
            let start_hit = start_cond.eval(&feature.attributes);
            let end_hit = end_cond.as_ref().map(|c| c.eval(&feature.attributes));
            tracks.entry(id).or_default().push(Obs {
                fid,
                x,
                y,
                t,
                start_hit,
                end_hit,
            });
        }

        // Assign each source feature an incident annotation.
        let mut annotations: BTreeMap<usize, Annotation> = BTreeMap::new();
        let mut incidents: Vec<Incident> = Vec::new();
        for (id, mut obs) in tracks {
            obs.sort_by(|a, b| a.t.total_cmp(&b.t));
            for run in find_incidents(&obs, end_cond.is_some()) {
                let incident_id = incidents.len() as i64 + 1;
                let start_t = obs[run[0]].t;
                let end_t = obs[run[run.len() - 1]].t;
                let duration = end_t - start_t;
                for (seq, &k) in run.iter().enumerate() {
                    let role = if seq == 0 {
                        "start"
                    } else if seq == run.len() - 1 {
                        "end"
                    } else {
                        "ongoing"
                    };
                    annotations.insert(
                        obs[k].fid,
                        Annotation {
                            incident_id,
                            status: role,
                            seq: seq as i64,
                            duration,
                        },
                    );
                }
                incidents.push(Incident {
                    track_id: id.clone(),
                    fids: run.iter().map(|&k| obs[k].fid).collect(),
                    coords: run.iter().map(|&k| (obs[k].x, obs[k].y)).collect(),
                    start_t,
                    end_t,
                    duration,
                });
            }
        }

        let flagged_points = annotations.len();
        let incident_count = incidents.len();
        ctx.progress.info(&format!(
            "{incident_count} incident(s) over {flagged_points} flagged point(s)"
        ));

        let geographic = layer.crs_epsg().map(|e| e == 4326).unwrap_or(true);
        let (out_layer, segment_count) = match prm.mode {
            Mode::Points => (build_points_layer(&layer, &annotations), 0usize),
            Mode::Segments => {
                let (l, n) = build_segments_layer(&layer, &incidents, geographic);
                (l, n)
            }
        };

        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("incident_count".to_string(), json!(incident_count));
        outputs.insert("flagged_points".to_string(), json!(flagged_points));
        outputs.insert("segment_count".to_string(), json!(segment_count));
        outputs.insert("skipped".to_string(), json!(skipped));
        Ok(ToolRunResult { outputs })
    }
}

// ── Incident detection ───────────────────────────────────────────────────────

#[derive(Clone, Copy)]
struct Obs {
    fid: usize,
    x: f64,
    y: f64,
    t: f64,
    start_hit: bool,
    end_hit: Option<bool>,
}

struct Annotation {
    incident_id: i64,
    status: &'static str,
    seq: i64,
    duration: f64,
}

struct Incident {
    track_id: String,
    fids: Vec<usize>,
    coords: Vec<(f64, f64)>,
    start_t: f64,
    end_t: f64,
    duration: f64,
}

/// Returns the incidents in one time-sorted track as lists of indices into `obs`.
///
/// Without an end condition each incident is a maximal run of consecutive points that
/// satisfy the start condition. With an end condition an incident opens at the first
/// start-condition point and stays open (collecting every point) until a later point
/// satisfies the end condition, or until the track ends.
fn find_incidents(obs: &[Obs], has_end: bool) -> Vec<Vec<usize>> {
    let mut incidents = Vec::new();
    if !has_end {
        let mut run: Vec<usize> = Vec::new();
        for (i, o) in obs.iter().enumerate() {
            if o.start_hit {
                run.push(i);
            } else if !run.is_empty() {
                incidents.push(std::mem::take(&mut run));
            }
        }
        if !run.is_empty() {
            incidents.push(run);
        }
    } else {
        let mut run: Vec<usize> = Vec::new();
        let mut open = false;
        for (i, o) in obs.iter().enumerate() {
            if !open {
                if o.start_hit {
                    open = true;
                    run.push(i);
                }
            } else {
                run.push(i);
                if o.end_hit.unwrap_or(false) {
                    incidents.push(std::mem::take(&mut run));
                    open = false;
                }
            }
        }
        if !run.is_empty() {
            incidents.push(run);
        }
    }
    incidents
}

// ── Output layers ────────────────────────────────────────────────────────────

/// Copies every input point through, appending the incident fields.
fn build_points_layer(input: &Layer, annotations: &BTreeMap<usize, Annotation>) -> Layer {
    let mut out = Layer::new("incidents").with_geom_type(GeometryType::Point);
    if let Some(epsg) = input.crs_epsg() {
        out = out.with_crs_epsg(epsg);
    }
    // Preserve the original schema, then append the incident fields.
    for f in input.schema.fields() {
        out.add_field(FieldDef::new(f.name.clone(), f.field_type));
    }
    out.add_field(FieldDef::new("incident_id", FieldType::Integer));
    out.add_field(FieldDef::new("incident_status", FieldType::Text));
    out.add_field(FieldDef::new("incident_seq", FieldType::Integer));
    out.add_field(FieldDef::new("incident_duration", FieldType::Float));

    for (fid, feature) in input.iter().enumerate() {
        let mut attrs = feature.attributes.clone();
        match annotations.get(&fid) {
            Some(a) => {
                attrs.push(FieldValue::Integer(a.incident_id));
                attrs.push(FieldValue::Text(a.status.to_string()));
                attrs.push(FieldValue::Integer(a.seq));
                attrs.push(FieldValue::Float(a.duration));
            }
            None => {
                attrs.push(FieldValue::Null);
                attrs.push(FieldValue::Text(String::new()));
                attrs.push(FieldValue::Null);
                attrs.push(FieldValue::Null);
            }
        }
        out.push(Feature {
            fid: 0,
            geometry: feature.geometry.clone(),
            attributes: attrs,
        });
    }
    out
}

/// Emits one polyline per incident that has at least two points.
fn build_segments_layer(input: &Layer, incidents: &[Incident], geographic: bool) -> (Layer, usize) {
    let mut out = Layer::new("incident_segments").with_geom_type(GeometryType::LineString);
    if let Some(epsg) = input.crs_epsg() {
        out = out.with_crs_epsg(epsg);
    }
    out.add_field(FieldDef::new("incident_id", FieldType::Integer));
    out.add_field(FieldDef::new("track_id", FieldType::Text));
    out.add_field(FieldDef::new("n_points", FieldType::Integer));
    out.add_field(FieldDef::new("start_t", FieldType::Float));
    out.add_field(FieldDef::new("end_t", FieldType::Float));
    out.add_field(FieldDef::new("duration", FieldType::Float));
    out.add_field(FieldDef::new("length", FieldType::Float));

    let mut count = 0usize;
    for (idx, inc) in incidents.iter().enumerate() {
        if inc.coords.len() < 2 {
            continue;
        }
        let coords: Vec<Coord> = inc.coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
        let mut length = 0.0;
        for w in coords.windows(2) {
            length += distance(w[0].x, w[0].y, w[1].x, w[1].y, geographic);
        }
        out.push(Feature {
            fid: 0,
            geometry: Some(Geometry::line_string(coords)),
            attributes: vec![
                FieldValue::Integer(idx as i64 + 1),
                FieldValue::Text(inc.track_id.clone()),
                FieldValue::Integer(inc.fids.len() as i64),
                FieldValue::Float(inc.start_t),
                FieldValue::Float(inc.end_t),
                FieldValue::Float(inc.duration),
                FieldValue::Float(length),
            ],
        });
        count += 1;
    }
    (out, count)
}

// ── Condition expression parsing ─────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Op {
    Gt,
    Ge,
    Lt,
    Le,
    Eq,
    Ne,
}

impl Op {
    fn test(self, lhs: f64, rhs: f64) -> bool {
        match self {
            Op::Gt => lhs > rhs,
            Op::Ge => lhs >= rhs,
            Op::Lt => lhs < rhs,
            Op::Le => lhs <= rhs,
            Op::Eq => lhs == rhs,
            Op::Ne => lhs != rhs,
        }
    }
}

/// A single `<field> <op> <number>` comparison (field name unresolved).
#[derive(Clone)]
struct Clause {
    field: String,
    op: Op,
    value: f64,
}

/// A conjunction (AND) of up to a few clauses, parsed from an expression string.
#[derive(Clone)]
struct Condition {
    clauses: Vec<Clause>,
}

/// A condition with its field names resolved to schema indices.
struct ResolvedCondition {
    clauses: Vec<(usize, Op, f64)>,
}

impl Condition {
    fn resolve(&self, layer: &Layer) -> Result<ResolvedCondition, ToolError> {
        let mut clauses = Vec::with_capacity(self.clauses.len());
        for c in &self.clauses {
            let idx = layer.schema.field_index(&c.field).ok_or_else(|| {
                ToolError::Validation(format!("condition field '{}' not found", c.field))
            })?;
            clauses.push((idx, c.op, c.value));
        }
        Ok(ResolvedCondition { clauses })
    }
}

impl ResolvedCondition {
    /// A point matches only if every clause is satisfiable (numeric) and true.
    fn eval(&self, attrs: &[FieldValue]) -> bool {
        self.clauses.iter().all(|&(idx, op, value)| {
            attrs
                .get(idx)
                .and_then(FieldValue::as_f64)
                .map(|lhs| op.test(lhs, value))
                .unwrap_or(false)
        })
    }
}

fn parse_condition(s: &str) -> Result<Condition, ToolError> {
    let parts = split_and(s);
    let mut clauses = Vec::new();
    for part in parts {
        clauses.push(parse_clause(&part)?);
    }
    if clauses.is_empty() {
        return Err(ToolError::Validation("empty condition expression".into()));
    }
    Ok(Condition { clauses })
}

/// Splits a conjunction on `&&` or a whitespace-delimited `AND` (case-insensitive).
fn split_and(s: &str) -> Vec<String> {
    if s.contains("&&") {
        return s
            .split("&&")
            .map(|p| p.trim().to_string())
            .filter(|p| !p.is_empty())
            .collect();
    }
    let lower = s.to_ascii_lowercase();
    if let Some(pos) = lower.find(" and ") {
        return vec![s[..pos].trim().to_string(), s[pos + 5..].trim().to_string()];
    }
    vec![s.trim().to_string()]
}

fn parse_clause(clause: &str) -> Result<Clause, ToolError> {
    // Two-character operators must be checked before their single-char prefixes.
    const OPS: [(&str, Op); 7] = [
        (">=", Op::Ge),
        ("<=", Op::Le),
        ("==", Op::Eq),
        ("!=", Op::Ne),
        (">", Op::Gt),
        ("<", Op::Lt),
        ("=", Op::Eq),
    ];
    for (tok, op) in OPS {
        if let Some(pos) = clause.find(tok) {
            let field = clause[..pos].trim();
            let rhs = clause[pos + tok.len()..].trim();
            if field.is_empty() {
                return Err(ToolError::Validation(format!(
                    "condition '{clause}' is missing a field name"
                )));
            }
            let value: f64 = rhs.parse().map_err(|_| {
                ToolError::Validation(format!(
                    "condition '{clause}' right-hand side '{rhs}' is not a number"
                ))
            })?;
            return Ok(Clause {
                field: field.to_string(),
                op,
                value,
            });
        }
    }
    Err(ToolError::Validation(format!(
        "condition '{clause}' has no comparison operator (expected one of > >= < <= == != )"
    )))
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

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    Points,
    Segments,
}

struct Params {
    track_field: String,
    time_field: String,
    start: Condition,
    end: Option<Condition>,
    mode: Mode,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let track_field = require_str(args, "track_field")?.to_string();
    let time_field = require_str(args, "time_field")?.to_string();
    let start = parse_condition(require_str(args, "start_condition")?)?;
    let end = match parse_optional_str(args, "end_condition")? {
        Some(s) => Some(parse_condition(s)?),
        None => None,
    };
    let mode = match parse_optional_str(args, "mode")? {
        None => Mode::Points,
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "points" | "point" => Mode::Points,
            "segments" | "segment" | "lines" => Mode::Segments,
            other => {
                return Err(ToolError::Validation(format!(
                    "mode '{other}' must be 'points' or 'segments'"
                )))
            }
        },
    };
    Ok(Params {
        track_field,
        time_field,
        start,
        end,
        mode,
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

    /// Builds a point layer with (track, time, value) rows in a projected CRS.
    fn track_layer(rows: &[(&str, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("id", FieldType::Text));
        l.add_field(FieldDef::new("t", FieldType::Float));
        l.add_field(FieldDef::new("speed", FieldType::Float));
        for (i, (id, t, v)) in rows.iter().enumerate() {
            l.add_feature(
                Some(Geometry::point(i as f64, 0.0)),
                &[
                    ("id", (*id).into()),
                    ("t", (*t).into()),
                    ("speed", (*v).into()),
                ],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = DetectIncidentsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// Two speeding runs in one track become two incidents; every flagged point
    /// satisfies the start condition.
    #[test]
    fn detects_two_episodes_and_flags_only_matching_points() {
        // speed profile: 0 5 6 1 2 7 8 9 1  -> runs [5,6] and [7,8,9] over > 3
        let rows = [
            ("a", 0.0, 0.0),
            ("a", 1.0, 5.0),
            ("a", 2.0, 6.0),
            ("a", 3.0, 1.0),
            ("a", 4.0, 2.0),
            ("a", 5.0, 7.0),
            ("a", 6.0, 8.0),
            ("a", 7.0, 9.0),
            ("a", 8.0, 1.0),
        ];
        let (out, layer) = run(json!({
            "input": track_layer(&rows),
            "track_field": "id",
            "time_field": "t",
            "start_condition": "speed > 3",
        }));
        assert_eq!(out.outputs["incident_count"], json!(2));
        assert_eq!(out.outputs["flagged_points"], json!(5));

        // Every flagged point must satisfy speed > 3; non-flagged must not.
        let iid = layer.schema.field_index("incident_id").unwrap();
        let sid = layer.schema.field_index("speed").unwrap();
        for f in &layer.features {
            let flagged = f.attributes[iid].as_i64().is_some();
            let speed = f.attributes[sid].as_f64().unwrap();
            assert_eq!(flagged, speed > 3.0, "flag disagrees with condition");
        }
    }

    /// With an end condition the incident stays open through sub-threshold points
    /// until the end condition fires (hysteresis).
    #[test]
    fn end_condition_holds_incident_open() {
        // start speed > 8, end speed < 2.  9 5 5 1 -> one incident of 4 points.
        let rows = [
            ("a", 0.0, 9.0),
            ("a", 1.0, 5.0),
            ("a", 2.0, 5.0),
            ("a", 3.0, 1.0),
            ("a", 4.0, 5.0), // after incident closed, no restart (needs >8)
        ];
        let (out, layer) = run(json!({
            "input": track_layer(&rows),
            "track_field": "id",
            "time_field": "t",
            "start_condition": "speed > 8",
            "end_condition": "speed < 2",
        }));
        assert_eq!(out.outputs["incident_count"], json!(1));
        assert_eq!(out.outputs["flagged_points"], json!(4));
        // The closing point (speed 1) carries status "end".
        let iid = layer.schema.field_index("incident_id").unwrap();
        let st = layer.schema.field_index("incident_status").unwrap();
        let closing = layer
            .features
            .iter()
            .find(|f| {
                f.attributes[iid].as_i64() == Some(1) && f.attributes[st].as_str() == Some("end")
            })
            .unwrap();
        // duration 3 (t 0..3)
        let dur = layer.schema.field_index("incident_duration").unwrap();
        assert!((closing.attributes[dur].as_f64().unwrap() - 3.0).abs() < 1e-9);
    }

    /// Two-clause AND condition.
    #[test]
    fn two_clause_and_condition() {
        // flag where speed > 3 AND speed < 8 : values 2 5 9 6 -> flagged [5,6]
        let rows = [
            ("a", 0.0, 2.0),
            ("a", 1.0, 5.0),
            ("a", 2.0, 9.0),
            ("a", 3.0, 6.0),
        ];
        let (out, _l) = run(json!({
            "input": track_layer(&rows),
            "track_field": "id",
            "time_field": "t",
            "start_condition": "speed > 3 AND speed < 8",
        }));
        // 5 and 6 qualify but 9 breaks the run between them -> two 1-point incidents.
        assert_eq!(out.outputs["incident_count"], json!(2));
        assert_eq!(out.outputs["flagged_points"], json!(2));
    }

    /// Segments mode makes one polyline per multi-point incident.
    #[test]
    fn segments_mode_builds_lines() {
        let rows = [
            ("a", 0.0, 0.0),
            ("a", 1.0, 5.0),
            ("a", 2.0, 6.0),
            ("a", 3.0, 7.0),
            ("a", 4.0, 0.0),
        ];
        let (out, layer) = run(json!({
            "input": track_layer(&rows),
            "track_field": "id",
            "time_field": "t",
            "start_condition": "speed > 3",
            "mode": "segments",
        }));
        assert_eq!(out.outputs["incident_count"], json!(1));
        assert_eq!(out.outputs["segment_count"], json!(1));
        assert_eq!(layer.features.len(), 1);
        assert!(matches!(
            layer.features[0].geometry,
            Some(Geometry::LineString(_))
        ));
        let np = layer.schema.field_index("n_points").unwrap();
        assert_eq!(layer.features[0].attributes[np].as_i64(), Some(3));
    }

    #[test]
    fn parses_operators() {
        assert_eq!(parse_clause("speed >= 4").unwrap().op, Op::Ge);
        assert_eq!(parse_clause("x <= 4").unwrap().op, Op::Le);
        assert_eq!(parse_clause("x != 4").unwrap().op, Op::Ne);
        assert_eq!(parse_clause("x == 4").unwrap().op, Op::Eq);
        assert_eq!(parse_clause("x = 4").unwrap().op, Op::Eq);
        assert!(parse_clause("speed 4").is_err());
        assert!(parse_clause("> 4").is_err());
        assert!(parse_clause("speed > fast").is_err());
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            DetectIncidentsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "input": "a.geojson", "track_field": "id", "time_field": "t" })).is_err()
        );
        assert!(bad(json!({
            "input": "a.geojson", "track_field": "id", "time_field": "t",
            "start_condition": "speed !! 3"
        }))
        .is_err());
        assert!(bad(json!({
            "input": "a.geojson", "track_field": "id", "time_field": "t",
            "start_condition": "speed > 3", "mode": "bogus"
        }))
        .is_err());
        assert!(bad(json!({
            "input": "a.geojson", "track_field": "id", "time_field": "t",
            "start_condition": "speed > 3"
        }))
        .is_ok());
    }
}
