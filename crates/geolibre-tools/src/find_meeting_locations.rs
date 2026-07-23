//! GeoLibre tool: find meeting locations from movement tracks.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Find Meeting Locations* (Intelligence
//! / GeoAnalytics). It detects **gatherings**: places and time windows where at
//! least `min_participants` *distinct* tracks are together (within
//! `search_distance`) and stay together for at least `min_meeting_duration`.
//!
//! The authored movement suite covers narrower cases:
//! * `trace_proximity_events` is *pairwise* contact tracing between two tracks,
//!   not N-track gatherings with a participant threshold;
//! * `reconstruct_tracks` finds *single-track* dwell locations, with no
//!   cross-track co-occurrence.
//!
//! Algorithm (deterministic, no RNG):
//! 1. Bin timestamped points into `time_step`-wide slices.
//! 2. In each slice, cluster points by connected components of the
//!    within-`search_distance` graph (union-find over a kd-tree radius search);
//!    keep components holding `>= min_participants` distinct tracks.
//! 3. Link qualifying components in consecutive slices whose centroids are
//!    within `search_distance` — an ongoing meeting that stays put.
//! 4. For each linked meeting compute participants, start/end time, duration,
//!    and emit a convex-hull **area** polygon plus a representative **point**.
//!    Meetings shorter than `min_meeting_duration` (or, if set, longer than
//!    `max_meeting_duration`) are dropped.
//!
//! Distances are in the layer's CRS units (use a projected CRS); times are
//! numeric seconds or ISO-8601 timestamps.

use std::collections::BTreeMap;

use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct FindMeetingLocationsTool;

impl Tool for FindMeetingLocationsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "find_meeting_locations",
            display_name: "Find Meeting Locations",
            summary: "Detect gatherings where multiple distinct tracks are together within a distance for at least a minimum duration (like ArcGIS Find Meeting Locations) — the N-track co-occurrence the authored trace_proximity_events (pairwise) and reconstruct_tracks (single-track dwell) don't cover. Emits area hulls and representative points with participant count, start/end time, and duration.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Timestamped track point layer (projected CRS; distances in map units).",
                    required: true,
                },
                ToolParamSpec {
                    name: "track_field",
                    description: "Field identifying each track/entity.",
                    required: true,
                },
                ToolParamSpec {
                    name: "time_field",
                    description: "Field holding each point's time: numeric seconds or an ISO-8601 timestamp.",
                    required: true,
                },
                ToolParamSpec {
                    name: "search_distance",
                    description: "Maximum distance (map units) for tracks to count as together.",
                    required: true,
                },
                ToolParamSpec {
                    name: "min_meeting_duration",
                    description: "Minimum duration (time units) a gathering must persist to be reported.",
                    required: true,
                },
                ToolParamSpec {
                    name: "max_meeting_duration",
                    description: "Optional maximum duration; longer gatherings are dropped.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_participants",
                    description: "Minimum number of distinct tracks in a gathering (default 2).",
                    required: false,
                },
                ToolParamSpec {
                    name: "time_step",
                    description: "Width of a time slice (time units). Default: min_meeting_duration.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output representative-point layer (one point per meeting). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_area",
                    description: "Optional output area layer (convex-hull polygon per meeting). If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for k in ["input", "track_field", "time_field"] {
            if args
                .get(k)
                .and_then(Value::as_str)
                .map(str::trim)
                .unwrap_or("")
                .is_empty()
            {
                return Err(ToolError::Validation(format!(
                    "missing required string parameter '{k}'"
                )));
            }
        }
        if args.get("search_distance").is_none() {
            return Err(ToolError::Validation(
                "missing required parameter 'search_distance'".to_string(),
            ));
        }
        if args.get("min_meeting_duration").is_none() {
            return Err(ToolError::Validation(
                "missing required parameter 'min_meeting_duration'".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let output_area = parse_optional_str(args, "output_area")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let track_idx = layer.schema.field_index(&prm.track_field).ok_or_else(|| {
            ToolError::Validation(format!("track_field '{}' not found", prm.track_field))
        })?;
        let time_idx = layer.schema.field_index(&prm.time_field).ok_or_else(|| {
            ToolError::Validation(format!("time_field '{}' not found", prm.time_field))
        })?;

        // Parse points: (track, time, x, y).
        let mut recs: Vec<Rec> = Vec::new();
        for feature in layer.features.iter() {
            let Some((x, y)) = feature.geometry.as_ref().and_then(point_xy) else {
                continue;
            };
            let Some(t) = feature.attributes.get(time_idx).and_then(parse_time_value) else {
                continue;
            };
            let track = feature
                .attributes
                .get(track_idx)
                .map(field_key)
                .unwrap_or_default();
            recs.push(Rec { track, t, x, y });
        }
        if recs.is_empty() {
            return Err(ToolError::Execution(
                "no usable timestamped point features".to_string(),
            ));
        }

        let t0 = recs.iter().map(|r| r.t).fold(f64::INFINITY, f64::min);
        ctx.progress
            .info(&format!("{} track point(s), binning by time", recs.len()));

        // ── Bin by time slice ────────────────────────────────────────────────
        let mut by_bin: BTreeMap<i64, Vec<usize>> = BTreeMap::new();
        for (i, r) in recs.iter().enumerate() {
            let bin = ((r.t - t0) / prm.time_step).floor() as i64;
            by_bin.entry(bin).or_default().push(i);
        }

        // ── Per-bin spatial clustering into candidate gathering instances ────
        let mut instances: Vec<Instance> = Vec::new();
        for (&bin, idxs) in &by_bin {
            for comp in components_within(&recs, idxs, prm.search_distance) {
                let tracks: std::collections::BTreeSet<&str> =
                    comp.iter().map(|&i| recs[i].track.as_str()).collect();
                if tracks.len() < prm.min_participants {
                    continue;
                }
                let (cx, cy) = centroid(&recs, &comp);
                instances.push(Instance {
                    bin,
                    members: comp,
                    cx,
                    cy,
                });
            }
        }

        ctx.progress.info(&format!(
            "{} candidate gathering instance(s); linking across time",
            instances.len()
        ));

        // ── Link consecutive-bin instances that stay in place ────────────────
        let mut uf = UnionFind::new(instances.len());
        for a in 0..instances.len() {
            for b in (a + 1)..instances.len() {
                let (ia, ib) = (&instances[a], &instances[b]);
                if (ia.bin - ib.bin).abs() != 1 {
                    continue;
                }
                let d = (ia.cx - ib.cx).hypot(ia.cy - ib.cy);
                if d <= prm.search_distance {
                    uf.union(a, b);
                }
            }
        }

        // Gather instances into meetings.
        let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for i in 0..instances.len() {
            groups.entry(uf.find(i)).or_default().push(i);
        }

        // ── Build meetings, keeping only lingering participants ──────────────
        // A track counts as a participant only if it is present in the gathering
        // for at least `min_meeting_duration` (its own member points span that
        // long) — this excludes passers-by who merely cross the gathering.
        let mut meetings: Vec<Meeting> = Vec::new();
        for insts in groups.values() {
            let mut members: Vec<usize> = Vec::new();
            for &ii in insts {
                members.extend_from_slice(&instances[ii].members);
            }
            members.sort_unstable();
            members.dedup();

            // Per-track time span within this gathering.
            let mut span: BTreeMap<&str, (f64, f64)> = BTreeMap::new();
            for &i in &members {
                let e = span
                    .entry(recs[i].track.as_str())
                    .or_insert((f64::INFINITY, f64::NEG_INFINITY));
                e.0 = e.0.min(recs[i].t);
                e.1 = e.1.max(recs[i].t);
            }
            let participants: std::collections::BTreeSet<&str> = span
                .iter()
                .filter(|(_, &(lo, hi))| hi - lo >= prm.min_meeting_duration)
                .map(|(&t, _)| t)
                .collect();
            if participants.len() < prm.min_participants {
                continue;
            }

            // Keep only the lingering participants' points for the geometry and
            // window.
            let members: Vec<usize> = members
                .into_iter()
                .filter(|&i| participants.contains(recs[i].track.as_str()))
                .collect();
            let start = members
                .iter()
                .map(|&i| recs[i].t)
                .fold(f64::INFINITY, f64::min);
            let end = members
                .iter()
                .map(|&i| recs[i].t)
                .fold(f64::NEG_INFINITY, f64::max);
            let duration = end - start;
            if let Some(maxd) = prm.max_meeting_duration {
                if duration > maxd {
                    continue;
                }
            }
            meetings.push(Meeting {
                members,
                participants: participants.len(),
                start,
                end,
                duration,
            });
        }
        // Stable, deterministic order: earliest start first.
        meetings.sort_by(|a, b| a.start.total_cmp(&b.start));

        ctx.progress
            .info(&format!("{} meeting(s) detected", meetings.len()));

        // ── Emit point and area layers ───────────────────────────────────────
        let mut pts_layer = Layer::new("meeting_points").with_geom_type(GeometryType::Point);
        let mut area_layer = Layer::new("meeting_areas").with_geom_type(GeometryType::Polygon);
        if let Some(epsg) = layer.crs_epsg() {
            pts_layer = pts_layer.with_crs_epsg(epsg);
            area_layer = area_layer.with_crs_epsg(epsg);
        }
        for l in [&mut pts_layer, &mut area_layer] {
            l.add_field(FieldDef::new("meeting_id", FieldType::Integer));
            l.add_field(FieldDef::new("participants", FieldType::Integer));
            l.add_field(FieldDef::new("point_count", FieldType::Integer));
            l.add_field(FieldDef::new("start_time", FieldType::Float));
            l.add_field(FieldDef::new("end_time", FieldType::Float));
            l.add_field(FieldDef::new("duration", FieldType::Float));
        }

        for (mid, m) in meetings.iter().enumerate() {
            let (cx, cy) = centroid(&recs, &m.members);
            let attrs = [
                ("meeting_id", (mid as i64).into()),
                ("participants", (m.participants as i64).into()),
                ("point_count", (m.members.len() as i64).into()),
                ("start_time", m.start.into()),
                ("end_time", m.end.into()),
                ("duration", m.duration.into()),
            ];
            pts_layer
                .add_feature(Some(Geometry::point(cx, cy)), &attrs)
                .map_err(|e| ToolError::Execution(format!("failed writing meeting point: {e}")))?;

            let hull = convex_hull(m.members.iter().map(|&i| (recs[i].x, recs[i].y)).collect());
            if hull.len() >= 3 {
                let ring: Vec<Coord> = hull.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
                area_layer
                    .add_feature(Some(Geometry::polygon(ring, vec![])), &attrs)
                    .map_err(|e| {
                        ToolError::Execution(format!("failed writing meeting area: {e}"))
                    })?;
            }
        }

        let meeting_count = meetings.len();
        let out_points = write_or_store_layer(pts_layer, output)?;
        let out_areas = write_or_store_layer(area_layer, output_area)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_points));
        outputs.insert("output_area".to_string(), json!(out_areas));
        outputs.insert("meeting_count".to_string(), json!(meeting_count));
        Ok(ToolRunResult { outputs })
    }
}

// ── Spatial clustering within a time bin ──────────────────────────────────────

/// Connected components of the "within `radius`" graph over `idxs` (record
/// indices), via a kd-tree radius search + union-find.
fn components_within(recs: &[Rec], idxs: &[usize], radius: f64) -> Vec<Vec<usize>> {
    let n = idxs.len();
    if n == 0 {
        return Vec::new();
    }
    // Local index -> record index; kd-tree over local indices.
    let mut tree: KdTree<f64, usize, [f64; 2]> = KdTree::new(2);
    for (li, &ri) in idxs.iter().enumerate() {
        tree.add([recs[ri].x, recs[ri].y], li).ok();
    }
    let r2 = radius * radius;
    let mut uf = UnionFind::new(n);
    for (li, &ri) in idxs.iter().enumerate() {
        if let Ok(found) = tree.within(&[recs[ri].x, recs[ri].y], r2, &squared_euclidean) {
            for (_d, &lj) in found {
                uf.union(li, lj);
            }
        }
    }
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for li in 0..n {
        groups.entry(uf.find(li)).or_default().push(idxs[li]);
    }
    groups.into_values().collect()
}

fn centroid(recs: &[Rec], members: &[usize]) -> (f64, f64) {
    let n = members.len() as f64;
    let sx: f64 = members.iter().map(|&i| recs[i].x).sum();
    let sy: f64 = members.iter().map(|&i| recs[i].y).sum();
    (sx / n, sy / n)
}

/// Andrew's monotone-chain convex hull. Returns the hull vertices CCW (no
/// closing duplicate). Fewer than 3 unique points → the unique points.
fn convex_hull(mut pts: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    pts.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    pts.dedup();
    let n = pts.len();
    if n < 3 {
        return pts;
    }
    let cross = |o: (f64, f64), a: (f64, f64), b: (f64, f64)| {
        (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
    };
    let mut hull: Vec<(f64, f64)> = Vec::with_capacity(2 * n);
    for &p in &pts {
        while hull.len() >= 2 && cross(hull[hull.len() - 2], hull[hull.len() - 1], p) <= 0.0 {
            hull.pop();
        }
        hull.push(p);
    }
    let lower_len = hull.len() + 1;
    for &p in pts.iter().rev() {
        while hull.len() >= lower_len && cross(hull[hull.len() - 2], hull[hull.len() - 1], p) <= 0.0
        {
            hull.pop();
        }
        hull.push(p);
    }
    hull.pop();
    hull
}

// ── Union-find ────────────────────────────────────────────────────────────────

struct UnionFind {
    parent: Vec<usize>,
}
impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
        }
    }
    fn find(&mut self, x: usize) -> usize {
        let mut r = x;
        while self.parent[r] != r {
            r = self.parent[r];
        }
        let mut c = x;
        while self.parent[c] != r {
            let next = self.parent[c];
            self.parent[c] = r;
            c = next;
        }
        r
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.parent[rb] = ra;
        }
    }
}

// ── Records ───────────────────────────────────────────────────────────────────

struct Rec {
    track: String,
    t: f64,
    x: f64,
    y: f64,
}

struct Instance {
    bin: i64,
    members: Vec<usize>,
    cx: f64,
    cy: f64,
}

struct Meeting {
    members: Vec<usize>,
    participants: usize,
    start: f64,
    end: f64,
    duration: f64,
}

// ── Time / geometry helpers (shared shape with reconstruct_tracks) ────────────

fn field_key(fv: &wbvector::FieldValue) -> String {
    if let Some(i) = fv.as_i64() {
        i.to_string()
    } else if let Some(f) = fv.as_f64() {
        format!("{f}")
    } else {
        fv.as_str().unwrap_or("").to_string()
    }
}

fn parse_time_value(fv: &wbvector::FieldValue) -> Option<f64> {
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

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    track_field: String,
    time_field: String,
    search_distance: f64,
    min_meeting_duration: f64,
    max_meeting_duration: Option<f64>,
    min_participants: usize,
    time_step: f64,
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

fn parse_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(Some(n.as_f64().unwrap_or(f64::NAN))),
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

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let track_field = require_str(args, "track_field")?.to_string();
    let time_field = require_str(args, "time_field")?.to_string();

    let search_distance = parse_f64(args, "search_distance")?
        .ok_or_else(|| ToolError::Validation("missing required 'search_distance'".to_string()))?;
    if search_distance <= 0.0 {
        return Err(ToolError::Validation(
            "'search_distance' must be > 0".to_string(),
        ));
    }
    let min_meeting_duration = parse_f64(args, "min_meeting_duration")?.ok_or_else(|| {
        ToolError::Validation("missing required 'min_meeting_duration'".to_string())
    })?;
    if min_meeting_duration < 0.0 {
        return Err(ToolError::Validation(
            "'min_meeting_duration' must be >= 0".to_string(),
        ));
    }
    let max_meeting_duration = parse_f64(args, "max_meeting_duration")?;
    if let Some(m) = max_meeting_duration {
        if m < min_meeting_duration {
            return Err(ToolError::Validation(
                "'max_meeting_duration' must be >= 'min_meeting_duration'".to_string(),
            ));
        }
    }
    let min_participants = match args.get("min_participants") {
        None | Some(Value::Null) => 2,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(2).max(2) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 2,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'min_participants' must be an integer".into()))?
            .max(2),
        Some(_) => {
            return Err(ToolError::Validation(
                "'min_participants' must be a number".into(),
            ))
        }
    };
    let time_step = match parse_f64(args, "time_step")? {
        None => {
            if min_meeting_duration > 0.0 {
                min_meeting_duration
            } else {
                1.0
            }
        }
        Some(v) if v > 0.0 => v,
        Some(_) => return Err(ToolError::Validation("'time_step' must be > 0".into())),
    };

    Ok(Params {
        track_field,
        time_field,
        search_distance,
        min_meeting_duration,
        max_meeting_duration,
        min_participants,
        time_step,
    })
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

    /// Each record is (track, time, x, y).
    fn layer_of(rows: &[(&str, f64, f64, f64)]) -> String {
        let mut l = Layer::new("track")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("track", FieldType::Text));
        l.add_field(FieldDef::new("t", FieldType::Float));
        for &(tr, t, x, y) in rows {
            l.add_feature(
                Some(Geometry::point(x, y)),
                &[("track", tr.into()), ("t", t.into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = FindMeetingLocationsTool.run(&args, &ctx()).unwrap();
        let pts = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, pts)
    }

    /// Three tracks converge at ~(0,0) over times 0..3, then disperse. One
    /// meeting with 3 participants is detected.
    #[test]
    fn detects_three_way_meeting() {
        let rows = vec![
            // gathering near origin, t = 0,1,2,3
            ("a", 0.0, 0.0, 0.0),
            ("b", 0.0, 1.0, 0.0),
            ("c", 0.0, 0.0, 1.0),
            ("a", 1.0, 0.5, 0.5),
            ("b", 1.0, 1.0, 1.0),
            ("c", 1.0, 0.0, 0.5),
            ("a", 2.0, 0.5, 0.0),
            ("b", 2.0, 1.0, 0.5),
            ("c", 2.0, 0.5, 1.0),
            // then a and b wander off far away, alone
            ("a", 5.0, 500.0, 500.0),
            ("b", 6.0, 900.0, 100.0),
        ];
        let (out, pts) = run(json!({
            "input": layer_of(&rows), "track_field": "track", "time_field": "t",
            "search_distance": 3.0, "min_meeting_duration": 1.0, "time_step": 1.0
        }));
        assert_eq!(out.outputs["meeting_count"], json!(1));
        let pidx = pts.schema.field_index("participants").unwrap();
        assert_eq!(pts.features[0].attributes[pidx].as_i64().unwrap(), 3);
    }

    /// Two tracks that pass each other but never linger produce no meeting when
    /// min_meeting_duration exceeds their brief proximity.
    #[test]
    fn brief_pass_is_not_a_meeting() {
        let rows = vec![
            ("a", 0.0, 0.0, 0.0),
            ("b", 0.0, 100.0, 0.0),
            ("a", 1.0, 1.0, 0.0), // momentarily close at t=1
            ("b", 1.0, 1.5, 0.0),
            ("a", 2.0, 2.0, 0.0),
            ("b", 2.0, 100.0, 0.0), // b leaves
        ];
        let (out, _p) = run(json!({
            "input": layer_of(&rows), "track_field": "track", "time_field": "t",
            "search_distance": 3.0, "min_meeting_duration": 5.0, "time_step": 1.0
        }));
        assert_eq!(
            out.outputs["meeting_count"],
            json!(0),
            "a brief pass is not a meeting"
        );
    }

    /// min_participants gates: two-track gatherings are dropped when 3 required.
    #[test]
    fn min_participants_gate() {
        let rows = vec![
            ("a", 0.0, 0.0, 0.0),
            ("b", 0.0, 1.0, 0.0),
            ("a", 1.0, 0.0, 0.5),
            ("b", 1.0, 1.0, 0.5),
            ("a", 2.0, 0.0, 1.0),
            ("b", 2.0, 1.0, 1.0),
        ];
        let (out, _p) = run(json!({
            "input": layer_of(&rows), "track_field": "track", "time_field": "t",
            "search_distance": 3.0, "min_meeting_duration": 1.0, "min_participants": 3, "time_step": 1.0
        }));
        assert_eq!(out.outputs["meeting_count"], json!(0));
    }

    /// The same points but same-track duplicates never form a meeting (need
    /// distinct tracks).
    #[test]
    fn single_track_is_not_a_meeting() {
        let rows = vec![
            ("a", 0.0, 0.0, 0.0),
            ("a", 1.0, 0.5, 0.5),
            ("a", 2.0, 1.0, 1.0),
        ];
        let (out, _p) = run(json!({
            "input": layer_of(&rows), "track_field": "track", "time_field": "t",
            "search_distance": 3.0, "min_meeting_duration": 1.0, "time_step": 1.0
        }));
        assert_eq!(out.outputs["meeting_count"], json!(0));
    }

    /// A passer-by who crosses the gathering at a single instant is not counted
    /// as a participant (must linger for >= min_meeting_duration).
    #[test]
    fn passerby_is_not_a_participant() {
        let mut rows = vec![
            // three tracks linger near origin over t = 0..3
            ("a", 0.0, 0.0, 0.0),
            ("b", 0.0, 1.0, 0.0),
            ("c", 0.0, 0.0, 1.0),
            ("a", 1.0, 0.5, 0.5),
            ("b", 1.0, 1.0, 1.0),
            ("c", 1.0, 0.0, 0.5),
            ("a", 2.0, 0.5, 0.0),
            ("b", 2.0, 1.0, 0.5),
            ("c", 2.0, 0.5, 1.0),
            ("a", 3.0, 0.0, 0.0),
            ("b", 3.0, 1.0, 0.0),
            ("c", 3.0, 0.0, 1.0),
        ];
        // passer-by 'p' is in the cluster only at t=1, elsewhere far away.
        rows.push(("p", 0.0, 900.0, 900.0));
        rows.push(("p", 1.0, 0.5, 0.5)); // momentarily inside the gathering
        rows.push(("p", 2.0, 900.0, 900.0));
        let (out, pts) = run(json!({
            "input": layer_of(&rows), "track_field": "track", "time_field": "t",
            "search_distance": 3.0, "min_meeting_duration": 2.0, "min_participants": 3, "time_step": 1.0
        }));
        assert_eq!(out.outputs["meeting_count"], json!(1));
        let pidx = pts.schema.field_index("participants").unwrap();
        assert_eq!(
            pts.features[0].attributes[pidx].as_i64().unwrap(),
            3,
            "the passer-by must not be counted"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            FindMeetingLocationsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "input": "a.geojson", "track_field": "t", "time_field": "ts" })).is_err()
        ); // no distance/duration
        assert!(bad(json!({
            "input": "a.geojson", "track_field": "t", "time_field": "ts",
            "search_distance": 10.0, "min_meeting_duration": 60.0
        }))
        .is_ok());
    }
}
