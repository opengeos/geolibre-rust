//! GeoLibre tool: join two point layers by spatial and temporal proximity.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Find Space Time Matches* (Crime
//! Analysis). The bundled `spatial_join` / `near` are space-only, and the
//! shipped `emerging_hot_spot_analysis` / `reconstruct_tracks` handle time
//! within a single layer; nothing joins **two** layers on space × time
//! simultaneously — crimes ↔ calls, sightings ↔ tracks, incidents ↔ sensor
//! events.
//!
//! For every primary feature, the secondary features that fall within
//! `search_distance` AND within `time_window` (before / after / either) are
//! emitted as matched pairs, each carrying the primary id, secondary id,
//! `distance`, and signed `delta_t` (secondary time − primary time, seconds).
//! Output is one straight link `LineString` per match (or the primary points
//! when the geometry is degenerate). Time fields are numeric or ISO-8601.

use std::collections::{BTreeMap, HashMap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

#[derive(Clone, Copy, PartialEq)]
enum Temporal {
    Before,
    After,
    Either,
}

pub struct FindSpaceTimeMatchesTool;

impl Tool for FindSpaceTimeMatchesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "find_space_time_matches",
            display_name: "Find Space Time Matches",
            summary: "Match features between two timestamped point layers that fall within a spatial search distance AND a time window of each other (before/after/either), emitting matched pairs with distance and time offset, like ArcGIS Find Space Time Matches.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "primary",
                    description: "Primary point layer (each feature is matched against the secondary layer).",
                    required: true,
                },
                ToolParamSpec {
                    name: "secondary",
                    description: "Secondary point layer to match against.",
                    required: true,
                },
                ToolParamSpec {
                    name: "time_field",
                    description: "Timestamp field present in BOTH layers (numeric or ISO-8601).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output match layer (one link line per matched pair). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "search_distance",
                    description: "Maximum spatial distance between a primary and secondary feature (map units).",
                    required: true,
                },
                ToolParamSpec {
                    name: "time_window",
                    description: "Maximum time separation: a number of seconds or a duration like '2h', '30m', '7d'.",
                    required: true,
                },
                ToolParamSpec {
                    name: "temporal_relationship",
                    description: "'either' (default), 'before' (secondary before primary), or 'after' (secondary after primary).",
                    required: false,
                },
                ToolParamSpec {
                    name: "primary_id_field",
                    description: "Field naming each primary feature in the output (default: feature index).",
                    required: false,
                },
                ToolParamSpec {
                    name: "secondary_id_field",
                    description: "Field naming each secondary feature in the output (default: feature index).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in [
            "primary",
            "secondary",
            "time_field",
            "search_distance",
            "time_window",
        ] {
            let missing = match args.get(key) {
                None | Some(Value::Null) => true,
                Some(Value::String(s)) => s.trim().is_empty(),
                _ => false,
            };
            if missing {
                return Err(ToolError::Validation(format!(
                    "missing required parameter '{key}'"
                )));
            }
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let primary_path = args.get("primary").and_then(Value::as_str).unwrap();
        let secondary_path = args.get("secondary").and_then(Value::as_str).unwrap();
        let time_field = args.get("time_field").and_then(Value::as_str).unwrap();
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let primary = load_input_layer(primary_path)?;
        let secondary = load_input_layer(secondary_path)?;

        let p_feats = collect(&primary, time_field, prm.primary_id_field.as_deref())?;
        let s_feats = collect(&secondary, time_field, prm.secondary_id_field.as_deref())?;
        if p_feats.is_empty() || s_feats.is_empty() {
            return Err(ToolError::Execution(
                "both layers must contain timestamped point features".to_string(),
            ));
        }

        ctx.progress.info(&format!(
            "matching {} primary × {} secondary features",
            p_feats.len(),
            s_feats.len()
        ));

        // Grid-hash the secondary layer for radius queries.
        let d = prm.search_distance;
        let cell = d.max(1e-9);
        let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
        for (i, s) in s_feats.iter().enumerate() {
            grid.entry(((s.x / cell).floor() as i64, (s.y / cell).floor() as i64))
                .or_default()
                .push(i);
        }

        let mut out = Layer::new("space_time_matches").with_geom_type(GeometryType::LineString);
        if let Some(epsg) = primary.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("primary_id", FieldType::Text));
        out.add_field(FieldDef::new("secondary_id", FieldType::Text));
        out.add_field(FieldDef::new("distance", FieldType::Float));
        out.add_field(FieldDef::new("delta_t", FieldType::Float));

        let d2 = d * d;
        let mut match_count = 0usize;
        let mut matched_primaries = 0usize;
        for p in &p_feats {
            let (gx, gy) = ((p.x / cell).floor() as i64, (p.y / cell).floor() as i64);
            let mut any = false;
            for dx in -1..=1 {
                for dy in -1..=1 {
                    let Some(bucket) = grid.get(&(gx + dx, gy + dy)) else {
                        continue;
                    };
                    for &si in bucket {
                        let s = &s_feats[si];
                        let dist2 = (p.x - s.x).powi(2) + (p.y - s.y).powi(2);
                        if dist2 > d2 {
                            continue;
                        }
                        let delta_t = s.t - p.t; // secondary − primary
                        let in_window = match prm.temporal {
                            Temporal::Either => delta_t.abs() <= prm.time_window,
                            Temporal::Before => delta_t <= 0.0 && -delta_t <= prm.time_window,
                            Temporal::After => delta_t >= 0.0 && delta_t <= prm.time_window,
                        };
                        if !in_window {
                            continue;
                        }
                        let geom =
                            Geometry::LineString(vec![Coord::xy(p.x, p.y), Coord::xy(s.x, s.y)]);
                        out.add_feature(
                            Some(geom),
                            &[
                                ("primary_id", FieldValue::Text(p.id.clone())),
                                ("secondary_id", FieldValue::Text(s.id.clone())),
                                ("distance", FieldValue::Float(dist2.sqrt())),
                                ("delta_t", FieldValue::Float(delta_t)),
                            ],
                        )
                        .map_err(|e| ToolError::Execution(format!("failed writing match: {e}")))?;
                        match_count += 1;
                        any = true;
                    }
                }
            }
            if any {
                matched_primaries += 1;
            }
        }

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("primary_count".to_string(), json!(p_feats.len()));
        outputs.insert("secondary_count".to_string(), json!(s_feats.len()));
        outputs.insert("match_count".to_string(), json!(match_count));
        outputs.insert("matched_primaries".to_string(), json!(matched_primaries));
        Ok(ToolRunResult { outputs })
    }
}

struct Feat {
    x: f64,
    y: f64,
    t: f64,
    id: String,
}

fn collect(
    layer: &Layer,
    time_field: &str,
    id_field: Option<&str>,
) -> Result<Vec<Feat>, ToolError> {
    let t_idx = layer
        .schema
        .field_index(time_field)
        .ok_or_else(|| ToolError::Validation(format!("time_field '{time_field}' not found")))?;
    let id_idx = match id_field {
        Some(f) => Some(
            layer
                .schema
                .field_index(f)
                .ok_or_else(|| ToolError::Validation(format!("id field '{f}' not found")))?,
        ),
        None => None,
    };
    let mut out = Vec::new();
    for (fidx, feature) in layer.features.iter().enumerate() {
        let Some(geom) = feature.geometry.as_ref() else {
            continue;
        };
        let Some((x, y)) = point_xy(geom) else {
            continue;
        };
        let Some(t) = parse_time_value(&feature.attributes[t_idx]) else {
            continue;
        };
        let id = match id_idx {
            Some(i) => value_string(&feature.attributes[i]),
            None => fidx.to_string(),
        };
        out.push(Feat { x, y, t, id });
    }
    Ok(out)
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

// ── Time parsing (shared shape with reconstruct_tracks / emerging_hot_spot) ───

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

/// Parses a duration: a plain number of seconds or `<n><unit>` (s/m/h/d/w/M/y).
fn parse_duration(s: &str) -> Result<f64, ToolError> {
    let s = s.trim();
    if let Ok(v) = s.parse::<f64>() {
        if v > 0.0 && v.is_finite() {
            return Ok(v);
        }
        return Err(ToolError::Validation(
            "'time_window' must be a positive number".to_string(),
        ));
    }
    let (num, unit) = s.split_at(s.len().saturating_sub(1));
    let value: f64 = num.trim().parse().map_err(|_| {
        ToolError::Validation(format!("could not parse 'time_window' value in '{s}'"))
    })?;
    if !(value > 0.0 && value.is_finite()) {
        return Err(ToolError::Validation(
            "'time_window' must be positive".to_string(),
        ));
    }
    let seconds = match unit {
        "s" => 1.0,
        "m" => 60.0,
        "h" => 3600.0,
        "d" => 86400.0,
        "w" => 604800.0,
        "M" => 2_592_000.0,
        "y" => 31_536_000.0,
        other => {
            return Err(ToolError::Validation(format!(
                "unknown 'time_window' unit '{other}' (use s/m/h/d/w/M/y or a plain number)"
            )))
        }
    };
    Ok(value * seconds)
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    search_distance: f64,
    time_window: f64,
    temporal: Temporal,
    primary_id_field: Option<String>,
    secondary_id_field: Option<String>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
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
    let time_window = match args.get("time_window") {
        Some(Value::Number(n)) => {
            let v = n.as_f64().unwrap_or(0.0);
            if v <= 0.0 {
                return Err(ToolError::Validation(
                    "'time_window' must be positive".to_string(),
                ));
            }
            v
        }
        Some(Value::String(s)) => parse_duration(s)?,
        _ => {
            return Err(ToolError::Validation(
                "missing required parameter 'time_window'".to_string(),
            ))
        }
    };
    let temporal =
        match parse_optional_str(args, "temporal_relationship")?.map(|s| s.trim().to_lowercase()) {
            None => Temporal::Either,
            Some(s) if s.is_empty() || s == "either" => Temporal::Either,
            Some(s) if s == "before" => Temporal::Before,
            Some(s) if s == "after" => Temporal::After,
            Some(other) => {
                return Err(ToolError::Validation(format!(
                    "'temporal_relationship' must be either|before|after, got '{other}'"
                )))
            }
        };
    Ok(Params {
        search_distance,
        time_window,
        temporal,
        primary_id_field: parse_optional_str(args, "primary_id_field")?.map(str::to_string),
        secondary_id_field: parse_optional_str(args, "secondary_id_field")?.map(str::to_string),
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

    /// pts: (x, y, time_seconds)
    fn layer_of(pts: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("t", FieldType::Float));
        for (x, y, t) in pts {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(*x, *y))),
                &[("t", (*t).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = FindSpaceTimeMatchesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// A pair close in space and time matches; distant-in-time does not.
    #[test]
    fn matches_within_space_and_time() {
        let p = layer_of(&[(0.0, 0.0, 1000.0)]);
        // s0: 3m away, 60s later (match). s1: 3m away, 10000s later (time-out).
        let s = layer_of(&[(3.0, 0.0, 1060.0), (3.0, 0.0, 11000.0)]);
        let (out, layer) = run(json!({
            "primary": p, "secondary": s, "time_field": "t",
            "search_distance": 10.0, "time_window": 300.0
        }));
        assert_eq!(out.outputs["match_count"], json!(1));
        let di = layer.schema.field_index("distance").unwrap();
        let ti = layer.schema.field_index("delta_t").unwrap();
        let f = layer.iter().next().unwrap();
        assert!((f.attributes[di].as_f64().unwrap() - 3.0).abs() < 1e-6);
        assert!((f.attributes[ti].as_f64().unwrap() - 60.0).abs() < 1e-6);
    }

    /// Spatial radius excludes far-away secondaries even inside the time window.
    #[test]
    fn excludes_spatially_distant() {
        let p = layer_of(&[(0.0, 0.0, 0.0)]);
        let s = layer_of(&[(1000.0, 0.0, 10.0)]);
        let (out, _l) = run(json!({
            "primary": p, "secondary": s, "time_field": "t",
            "search_distance": 50.0, "time_window": 3600.0
        }));
        assert_eq!(out.outputs["match_count"], json!(0));
    }

    /// temporal_relationship=before keeps only secondaries earlier than primary.
    #[test]
    fn before_relationship_filters_direction() {
        let p = layer_of(&[(0.0, 0.0, 1000.0)]);
        // earlier (900) and later (1100), both within 200s and 5m.
        let s = layer_of(&[(1.0, 0.0, 900.0), (1.0, 0.0, 1100.0)]);
        let (out, layer) = run(json!({
            "primary": p, "secondary": s, "time_field": "t",
            "search_distance": 10.0, "time_window": 200.0, "temporal_relationship": "before"
        }));
        assert_eq!(out.outputs["match_count"], json!(1));
        let ti = layer.schema.field_index("delta_t").unwrap();
        // secondary before primary -> delta_t negative.
        assert!(
            layer.iter().next().unwrap().attributes[ti]
                .as_f64()
                .unwrap()
                < 0.0
        );
    }

    /// A duration string like '2h' parses for the time window.
    #[test]
    fn duration_string_window() {
        let p = layer_of(&[(0.0, 0.0, 0.0)]);
        let s = layer_of(&[(1.0, 0.0, 3600.0)]); // 1h later
        let (out, _l) = run(json!({
            "primary": p, "secondary": s, "time_field": "t",
            "search_distance": 10.0, "time_window": "2h"
        }));
        assert_eq!(out.outputs["match_count"], json!(1), "1h within 2h window");
        let (out2, _l) = run(json!({
            "primary": p, "secondary": s, "time_field": "t",
            "search_distance": 10.0, "time_window": "30m"
        }));
        assert_eq!(
            out2.outputs["match_count"],
            json!(0),
            "1h outside 30m window"
        );
    }

    /// ISO-8601 timestamps parse and match.
    #[test]
    fn iso8601_timestamps() {
        let mut lp = Layer::new("p")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        lp.add_field(FieldDef::new("t", FieldType::Text));
        lp.add_feature(
            Some(Geometry::Point(Coord::xy(0.0, 0.0))),
            &[("t", "2024-01-01T12:00:00".into())],
        )
        .unwrap();
        let pid = memory_store::put_vector(lp);
        let mut ls = Layer::new("s")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        ls.add_field(FieldDef::new("t", FieldType::Text));
        ls.add_feature(
            Some(Geometry::Point(Coord::xy(2.0, 0.0))),
            &[("t", "2024-01-01T12:30:00".into())],
        )
        .unwrap();
        let sid = memory_store::put_vector(ls);
        let (out, _l) = run(json!({
            "primary": memory_store::make_vector_memory_path(&pid),
            "secondary": memory_store::make_vector_memory_path(&sid),
            "time_field": "t", "search_distance": 10.0, "time_window": "1h"
        }));
        assert_eq!(
            out.outputs["match_count"],
            json!(1),
            "30 min apart within 1h"
        );
    }

    #[test]
    fn rejects_missing_time_window() {
        let p = layer_of(&[(0.0, 0.0, 0.0)]);
        let s = layer_of(&[(1.0, 0.0, 0.0)]);
        let args: ToolArgs = serde_json::from_value(json!({
            "primary": p, "secondary": s, "time_field": "t", "search_distance": 10.0
        }))
        .unwrap();
        assert!(FindSpaceTimeMatchesTool.validate(&args).is_err());
    }
}
