//! GeoLibre tool: detect stationary dwell segments within movement tracks.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Find Dwell Locations* (GeoAnalytics).
//! The repo's movement suite is otherwise complete — `reconstruct_tracks`,
//! `snap_tracks`, `calculate_motion_statistics`, `trace_proximity_events`,
//! `find_meeting_locations`, `detect_incidents` — but the single-track
//! stop/dwell primitive is missing, and dwell segmentation is typically the
//! *first* step of a mobility workflow: it separates travel from stationary
//! periods before anything downstream runs.
//!
//! `find_meeting_locations` solves a different problem: it needs
//! `min_participants` and detects convergence of **multiple distinct tracks**,
//! so a single vehicle idling at a depot for two hours produces no meeting and
//! is invisible to it. `detect_incidents` triggers on attribute conditions, not
//! on spatial persistence.
//!
//! A dwell is a maximal run of consecutive fixes from one track that all stay
//! within `distance_tolerance` of the run's own running centroid, lasting at
//! least `time_tolerance`. Testing against the running centroid rather than the
//! first fix is what keeps slow drift (GPS jitter around a parked vehicle) in
//! one dwell instead of fragmenting it.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// What geometry to emit per dwell.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum OutputType {
    /// Every constituent fix, tagged with its dwell id.
    DwellFeatures,
    /// One point per dwell at its centroid.
    MeanCenters,
    /// One polygon per dwell: the convex hull of its fixes.
    ConvexHulls,
    /// Every input fix, tagged with a dwell id or -1 when moving.
    AllFeatures,
}

pub struct FindDwellLocationsTool;

impl Tool for FindDwellLocationsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "find_dwell_locations",
            display_name: "Find Dwell Locations",
            summary: "Detect where a moving entity stayed put: runs of track fixes within a distance tolerance for at least a time tolerance, emitted as fixes, mean centers or convex hulls, like ArcGIS Find Dwell Locations.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Timestamped point features (track fixes).",
                    required: true,
                },
                ToolParamSpec {
                    name: "track_field",
                    description: "Field identifying the track each fix belongs to.",
                    required: true,
                },
                ToolParamSpec {
                    name: "time_field",
                    description: "Timestamp field (numeric seconds or ISO-8601 text).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distance_tolerance",
                    description: "Maximum distance (map units) a fix may sit from the dwell's running centroid.",
                    required: true,
                },
                ToolParamSpec {
                    name: "time_tolerance",
                    description: "Minimum duration (same units as time_field) for a stationary run to count as a dwell.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output_type",
                    description: "dwell_features (default) | mean_centers | convex_hulls | all_features.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "track_field")?;
        require_str(args, "time_field")?;
        let d = parse_optional_f64(args, "distance_tolerance")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'distance_tolerance'".into())
        })?;
        if !d.is_finite() || d < 0.0 {
            return Err(ToolError::Validation(
                "'distance_tolerance' must be zero or greater".into(),
            ));
        }
        let t = parse_optional_f64(args, "time_tolerance")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'time_tolerance'".into())
        })?;
        if !t.is_finite() || t < 0.0 {
            return Err(ToolError::Validation(
                "'time_tolerance' must be zero or greater".into(),
            ));
        }
        parse_output_type(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let track_field = require_str(args, "track_field")?;
        let time_field = require_str(args, "time_field")?;
        let output = parse_optional_str(args, "output")?;
        let dist_tol = parse_optional_f64(args, "distance_tolerance")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'distance_tolerance'".to_string())
        })?;
        let time_tol = parse_optional_f64(args, "time_tolerance")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'time_tolerance'".to_string())
        })?;
        let output_type = parse_output_type(args)?;

        let layer = load_input_layer(input)?;
        for f in [track_field, time_field] {
            if layer.schema.field_index(f).is_none() {
                return Err(ToolError::Validation(format!(
                    "field '{f}' not found on the input layer"
                )));
            }
        }

        // Group fixes by track, keeping the source row index so all_features can
        // re-emit the original geometry and attributes.
        let mut tracks: BTreeMap<String, Vec<Fix>> = BTreeMap::new();
        for (i, feat) in layer.features.iter().enumerate() {
            let (Some(geom), Ok(tv), Ok(time_v)) = (
                feat.geometry.as_ref(),
                feat.get(&layer.schema, track_field),
                feat.get(&layer.schema, time_field),
            ) else {
                continue;
            };
            let (Some((x, y)), Some(t)) = (point_xy(geom), parse_time_value(time_v)) else {
                continue;
            };
            tracks
                .entry(field_string(tv))
                .or_default()
                .push(Fix { row: i, x, y, t });
        }

        ctx.progress
            .info(&format!("scanning {} track(s) for dwells", tracks.len()));

        // Detect dwells per track. Tracks are independent.
        let mut dwells: Vec<Dwell> = Vec::new();
        for (name, fixes) in tracks.iter_mut() {
            fixes.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
            detect_dwells(name, fixes, dist_tol, time_tol, &mut dwells);
        }
        // Deterministic ordering regardless of map iteration.
        dwells.sort_by(|a, b| {
            a.track.cmp(&b.track).then(
                a.start
                    .partial_cmp(&b.start)
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
        });

        let mut out = Layer::new(layer.name.clone());
        out.crs = layer.crs.clone();
        out.geom_type = Some(match output_type {
            OutputType::ConvexHulls => GeometryType::Polygon,
            // dwell_features / all_features clone the input geometry through
            // verbatim, and point_xy accepts MultiPoint, so the written type is
            // whatever the input was -- inherit it rather than asserting Point.
            OutputType::DwellFeatures | OutputType::AllFeatures => {
                layer.geom_type.unwrap_or(GeometryType::Point)
            }
            OutputType::MeanCenters => GeometryType::Point,
        });
        if matches!(
            output_type,
            OutputType::DwellFeatures | OutputType::AllFeatures
        ) {
            for fd in layer.schema.fields().iter() {
                out.add_field(fd.clone());
            }
        } else {
            out.add_field(FieldDef::new("track", FieldType::Text));
        }
        out.add_field(FieldDef::new("dwell_id", FieldType::Integer));
        out.add_field(FieldDef::new("dwell_start", FieldType::Float));
        out.add_field(FieldDef::new("dwell_end", FieldType::Float));
        out.add_field(FieldDef::new("dwell_duration", FieldType::Float));
        out.add_field(FieldDef::new("dwell_count", FieldType::Integer));

        // Row index -> dwell, for the per-fix output modes.
        let mut row_to_dwell: BTreeMap<usize, usize> = BTreeMap::new();
        for (di, d) in dwells.iter().enumerate() {
            for f in &d.fixes {
                row_to_dwell.insert(f.row, di);
            }
        }

        let mut emitted = 0usize;
        match output_type {
            OutputType::MeanCenters | OutputType::ConvexHulls => {
                for (di, d) in dwells.iter().enumerate() {
                    let geom = if output_type == OutputType::MeanCenters {
                        let (cx, cy) = centroid(&d.fixes);
                        Geometry::Point(Coord::xy(cx, cy))
                    } else {
                        match convex_hull(&d.fixes) {
                            Some(ring) => Geometry::Polygon {
                                exterior: ring,
                                interiors: vec![],
                            },
                            // Fewer than 3 distinct points cannot form a hull;
                            // fall back to the centroid so the dwell is not lost.
                            None => {
                                let (cx, cy) = centroid(&d.fixes);
                                Geometry::Point(Coord::xy(cx, cy))
                            }
                        }
                    };
                    out.add_feature(Some(geom), &dwell_fields(Some(&d.track), di, d))
                        .map_err(|e| ToolError::Execution(format!("failed writing dwell: {e}")))?;
                    emitted += 1;
                }
            }
            OutputType::DwellFeatures | OutputType::AllFeatures => {
                for (i, feat) in layer.features.iter().enumerate() {
                    let member = row_to_dwell.get(&i).copied();
                    if member.is_none() && output_type == OutputType::DwellFeatures {
                        continue;
                    }
                    let mut fields: Vec<(String, FieldValue)> = layer
                        .schema
                        .fields()
                        .iter()
                        .enumerate()
                        .map(|(fi, fd)| (fd.name.clone(), feat.attributes[fi].clone()))
                        .collect();
                    match member {
                        Some(di) => {
                            let d = &dwells[di];
                            fields.push(("dwell_id".into(), FieldValue::Integer(di as i64)));
                            fields.push(("dwell_start".into(), FieldValue::Float(d.start)));
                            fields.push(("dwell_end".into(), FieldValue::Float(d.end)));
                            fields.push((
                                "dwell_duration".into(),
                                FieldValue::Float(d.end - d.start),
                            ));
                            fields.push((
                                "dwell_count".into(),
                                FieldValue::Integer(d.fixes.len() as i64),
                            ));
                        }
                        None => {
                            // Moving fix: -1 marks it explicitly rather than null,
                            // so a filter on dwell_id >= 0 selects stationary fixes.
                            fields.push(("dwell_id".into(), FieldValue::Integer(-1)));
                            for f in ["dwell_start", "dwell_end", "dwell_duration"] {
                                fields.push((f.into(), FieldValue::Null));
                            }
                            fields.push(("dwell_count".into(), FieldValue::Integer(0)));
                        }
                    }
                    let refs: Vec<(&str, FieldValue)> = fields
                        .iter()
                        .map(|(k, v)| (k.as_str(), v.clone()))
                        .collect();
                    out.add_feature(feat.geometry.clone(), &refs)
                        .map_err(|e| ToolError::Execution(format!("failed writing fix: {e}")))?;
                    emitted += 1;
                }
            }
        }

        let total_duration: f64 = dwells.iter().map(|d| d.end - d.start).sum();
        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("dwell_count".to_string(), json!(dwells.len()));
        outputs.insert("track_count".to_string(), json!(tracks.len()));
        outputs.insert("feature_count".to_string(), json!(emitted));
        outputs.insert("total_dwell_duration".to_string(), json!(total_duration));
        Ok(ToolRunResult { outputs })
    }
}

#[derive(Clone)]
struct Fix {
    row: usize,
    x: f64,
    y: f64,
    t: f64,
}

struct Dwell {
    track: String,
    start: f64,
    end: f64,
    fixes: Vec<Fix>,
}

fn dwell_fields<'a>(track: Option<&'a str>, di: usize, d: &Dwell) -> Vec<(&'a str, FieldValue)> {
    let mut v: Vec<(&str, FieldValue)> = Vec::new();
    if let Some(t) = track {
        v.push(("track", FieldValue::Text(t.to_string())));
    }
    v.push(("dwell_id", FieldValue::Integer(di as i64)));
    v.push(("dwell_start", FieldValue::Float(d.start)));
    v.push(("dwell_end", FieldValue::Float(d.end)));
    v.push(("dwell_duration", FieldValue::Float(d.end - d.start)));
    v.push(("dwell_count", FieldValue::Integer(d.fixes.len() as i64)));
    v
}

/// Greedy maximal-run segmentation against the run's running centroid.
fn detect_dwells(track: &str, fixes: &[Fix], dist_tol: f64, time_tol: f64, out: &mut Vec<Dwell>) {
    let n = fixes.len();
    let mut i = 0usize;
    while i < n {
        let (mut sx, mut sy) = (fixes[i].x, fixes[i].y);
        let mut j = i + 1;
        while j < n {
            let k = (j - i) as f64;
            let (cx, cy) = (sx / k, sy / k);
            let dx = fixes[j].x - cx;
            let dy = fixes[j].y - cy;
            if (dx * dx + dy * dy).sqrt() > dist_tol {
                break;
            }
            sx += fixes[j].x;
            sy += fixes[j].y;
            j += 1;
        }
        let start = fixes[i].t;
        let end = fixes[j - 1].t;
        if j - i >= 2 && (end - start) >= time_tol {
            out.push(Dwell {
                track: track.to_string(),
                start,
                end,
                fixes: fixes[i..j].to_vec(),
            });
            i = j; // consume the whole dwell
        } else {
            // Not a dwell: advance one fix so an overlapping run can still start
            // at the next position rather than skipping the whole failed window.
            i += 1;
        }
    }
}

fn centroid(fixes: &[Fix]) -> (f64, f64) {
    let n = fixes.len() as f64;
    (
        fixes.iter().map(|f| f.x).sum::<f64>() / n,
        fixes.iter().map(|f| f.y).sum::<f64>() / n,
    )
}

/// Monotone-chain convex hull. Returns `None` when the points are collinear or
/// fewer than 3 distinct positions exist (no polygon can be formed).
fn convex_hull(fixes: &[Fix]) -> Option<Ring> {
    let mut pts: Vec<(f64, f64)> = fixes.iter().map(|f| (f.x, f.y)).collect();
    pts.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
    });
    pts.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-12 && (a.1 - b.1).abs() < 1e-12);
    if pts.len() < 3 {
        return None;
    }
    let cross = |o: (f64, f64), a: (f64, f64), b: (f64, f64)| {
        (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
    };
    let mut hull: Vec<(f64, f64)> = Vec::with_capacity(pts.len() * 2);
    for &p in pts.iter() {
        while hull.len() >= 2 && cross(hull[hull.len() - 2], hull[hull.len() - 1], p) <= 0.0 {
            hull.pop();
        }
        hull.push(p);
    }
    let lower = hull.len() + 1;
    for &p in pts.iter().rev() {
        while hull.len() >= lower && cross(hull[hull.len() - 2], hull[hull.len() - 1], p) <= 0.0 {
            hull.pop();
        }
        hull.push(p);
    }
    hull.pop(); // last point duplicates the first
    if hull.len() < 3 {
        return None;
    }
    Some(Ring::new(
        hull.into_iter().map(|(x, y)| Coord::xy(x, y)).collect(),
    ))
}

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        _ => None,
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

fn field_string(v: &FieldValue) -> String {
    match v {
        FieldValue::Text(s) | FieldValue::Date(s) | FieldValue::DateTime(s) => s.clone(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(x) => x.to_string(),
        FieldValue::Boolean(b) => b.to_string(),
        FieldValue::Null | FieldValue::Blob(_) => String::new(),
    }
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
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

fn parse_output_type(args: &ToolArgs) -> Result<OutputType, ToolError> {
    match args
        .get("output_type")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some("dwell_features") => Ok(OutputType::DwellFeatures),
        Some("mean_centers") => Ok(OutputType::MeanCenters),
        Some("convex_hulls") => Ok(OutputType::ConvexHulls),
        Some("all_features") => Ok(OutputType::AllFeatures),
        Some(o) => Err(ToolError::Validation(format!(
            "'output_type' must be dwell_features/mean_centers/convex_hulls/all_features, got '{o}'"
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

    /// Fixes as (track, x, y, t).
    fn track_layer(items: Vec<(&str, f64, f64, f64)>) -> String {
        let mut l = Layer::new("fixes")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("trk", FieldType::Text));
        l.add_field(FieldDef::new("t", FieldType::Float));
        for (trk, x, y, t) in items {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(x, y))),
                &[
                    ("trk", FieldValue::Text(trk.to_string())),
                    ("t", FieldValue::Float(t)),
                ],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let mut v = args;
        v["track_field"] = json!("trk");
        v["time_field"] = json!("t");
        let args: ToolArgs = serde_json::from_value(v).unwrap();
        let out = FindDwellLocationsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// A stationary run long enough in time is a dwell; the travel legs are not.
    #[test]
    fn detects_a_stationary_run() {
        let input = track_layer(vec![
            ("a", 0.0, 0.0, 0.0),     // travelling
            ("a", 100.0, 0.0, 10.0),  // travelling
            ("a", 200.0, 0.0, 20.0),  // arrives
            ("a", 200.5, 0.5, 80.0),  // parked
            ("a", 200.2, 0.1, 140.0), // parked
            ("a", 400.0, 0.0, 200.0), // leaves
        ]);
        let (out, layer) = run(json!({
            "input": input, "distance_tolerance": 5.0, "time_tolerance": 60.0
        }));
        assert_eq!(out.outputs["dwell_count"], json!(1));
        assert_eq!(out.outputs["total_dwell_duration"], json!(120.0));
        assert_eq!(out.outputs["feature_count"], json!(3));
        assert_eq!(layer.len(), 3);
    }

    /// A stationary run that is too brief does not qualify.
    #[test]
    fn short_stop_is_not_a_dwell() {
        let input = track_layer(vec![
            ("a", 0.0, 0.0, 0.0),
            ("a", 0.1, 0.1, 5.0),
            ("a", 500.0, 0.0, 100.0),
        ]);
        let (out, _) = run(json!({
            "input": input, "distance_tolerance": 5.0, "time_tolerance": 60.0
        }));
        assert_eq!(out.outputs["dwell_count"], json!(0));
    }

    /// THE distinction from find_meeting_locations: a single track with no
    /// other participants still produces a dwell.
    #[test]
    fn single_track_dwell_is_found() {
        let input = track_layer(vec![
            ("solo", 10.0, 10.0, 0.0),
            ("solo", 10.1, 10.0, 3600.0),
            ("solo", 10.0, 10.2, 7200.0),
        ]);
        let (out, _) = run(json!({
            "input": input, "distance_tolerance": 2.0, "time_tolerance": 600.0
        }));
        assert_eq!(
            out.outputs["dwell_count"],
            json!(1),
            "a lone idling vehicle is a dwell, not a meeting"
        );
    }

    /// Tracks are segmented independently.
    #[test]
    fn tracks_are_independent() {
        let input = track_layer(vec![
            ("a", 0.0, 0.0, 0.0),
            ("a", 0.1, 0.0, 100.0),
            // b sits at the same place but is a different entity.
            ("b", 0.0, 0.0, 0.0),
            ("b", 0.1, 0.0, 100.0),
        ]);
        let (out, _) = run(json!({
            "input": input, "distance_tolerance": 5.0, "time_tolerance": 50.0
        }));
        assert_eq!(out.outputs["dwell_count"], json!(2));
        assert_eq!(out.outputs["track_count"], json!(2));
    }

    /// mean_centers emits one point per dwell at its centroid.
    #[test]
    fn mean_centers_output() {
        let input = track_layer(vec![
            ("a", 0.0, 0.0, 0.0),
            ("a", 2.0, 0.0, 100.0),
            ("a", 1.0, 3.0, 200.0),
        ]);
        let (out, layer) = run(json!({
            "input": input, "distance_tolerance": 5.0,
            "time_tolerance": 50.0, "output_type": "mean_centers"
        }));
        assert_eq!(out.outputs["dwell_count"], json!(1));
        assert_eq!(layer.len(), 1);
        match layer.features[0].geometry.as_ref().unwrap() {
            Geometry::Point(c) => {
                assert!((c.x - 1.0).abs() < 1e-9);
                assert!((c.y - 1.0).abs() < 1e-9);
            }
            other => panic!("expected Point, got {other:?}"),
        }
    }

    /// convex_hulls emits a polygon per dwell.
    #[test]
    fn convex_hulls_output() {
        let input = track_layer(vec![
            ("a", 0.0, 0.0, 0.0),
            ("a", 3.0, 0.0, 100.0),
            ("a", 0.0, 3.0, 200.0),
            ("a", 3.0, 3.0, 300.0),
        ]);
        let (_, layer) = run(json!({
            "input": input, "distance_tolerance": 10.0,
            "time_tolerance": 50.0, "output_type": "convex_hulls"
        }));
        assert_eq!(layer.len(), 1);
        assert!(matches!(
            layer.features[0].geometry.as_ref().unwrap(),
            Geometry::Polygon { .. }
        ));
    }

    /// all_features keeps every fix, tagging moving ones with -1.
    #[test]
    fn all_features_tags_moving_fixes() {
        let input = track_layer(vec![
            ("a", 0.0, 0.0, 0.0),
            ("a", 0.1, 0.0, 100.0),
            ("a", 900.0, 0.0, 200.0),
        ]);
        let (_, layer) = run(json!({
            "input": input, "distance_tolerance": 5.0,
            "time_tolerance": 50.0, "output_type": "all_features"
        }));
        assert_eq!(layer.len(), 3, "every input fix survives");
        let di = layer.schema.field_index("dwell_id").unwrap();
        let ids: Vec<i64> = layer
            .features
            .iter()
            .map(|f| f.attributes[di].as_f64().unwrap() as i64)
            .collect();
        assert_eq!(ids.iter().filter(|v| **v >= 0).count(), 2);
        assert_eq!(ids.iter().filter(|v| **v < 0).count(), 1);
    }

    /// ISO-8601 timestamps are accepted, not just numbers.
    #[test]
    fn accepts_iso8601_timestamps() {
        let mut l = Layer::new("fixes")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("trk", FieldType::Text));
        l.add_field(FieldDef::new("t", FieldType::Text));
        for (t, x) in [("2026-01-01T00:00:00", 0.0), ("2026-01-01T01:00:00", 0.1)] {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(x, 0.0))),
                &[
                    ("trk", FieldValue::Text("a".into())),
                    ("t", FieldValue::Text(t.into())),
                ],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        let (out, _) = run(json!({
            "input": memory_store::make_vector_memory_path(&id),
            "distance_tolerance": 5.0, "time_tolerance": 600.0
        }));
        assert_eq!(out.outputs["dwell_count"], json!(1));
        assert_eq!(out.outputs["total_dwell_duration"], json!(3600.0));
    }

    /// Drift is tracked against the running centroid, so slow GPS wander around
    /// a parked vehicle stays a single dwell rather than fragmenting.
    #[test]
    fn slow_drift_stays_one_dwell() {
        let mut items = Vec::new();
        for i in 0..10 {
            items.push(("a", i as f64 * 0.3, 0.0, i as f64 * 60.0));
        }
        let input = track_layer(items);
        let (out, _) = run(json!({
            "input": input, "distance_tolerance": 5.0, "time_tolerance": 60.0
        }));
        assert_eq!(out.outputs["dwell_count"], json!(1));
    }

    #[test]
    fn rejects_bad_parameters() {
        let p = track_layer(vec![("a", 0.0, 0.0, 0.0)]);
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            FindDwellLocationsTool.validate(&args).is_err()
        };
        assert!(bad(json!({ "track_field": "trk", "time_field": "t",
                           "distance_tolerance": 1, "time_tolerance": 1 })));
        assert!(bad(json!({ "input": p, "time_field": "t",
                           "distance_tolerance": 1, "time_tolerance": 1 })));
        assert!(bad(
            json!({ "input": p, "track_field": "trk", "time_field": "t",
                           "time_tolerance": 1 })
        ));
        assert!(bad(
            json!({ "input": p, "track_field": "trk", "time_field": "t",
                           "distance_tolerance": -1, "time_tolerance": 1 })
        ));
        assert!(bad(
            json!({ "input": p, "track_field": "trk", "time_field": "t",
                           "distance_tolerance": 1, "time_tolerance": 1,
                           "output_type": "bogus" })
        ));
    }
}
