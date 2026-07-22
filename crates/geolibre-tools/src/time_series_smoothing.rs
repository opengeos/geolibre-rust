//! GeoLibre tool: per-location smoothing of an in-sample time series.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Time Series Smoothing* (Space Time
//! Pattern Mining). Where `time_series_forecast` projects a series *ahead*
//! (Holt + polynomial curve fits), this tool smooths the *observed* series in
//! place and appends the smoothed value to every feature, leaving geometry and
//! all original attributes untouched.
//!
//! Features are grouped by `id_field` (each group is one location's series),
//! ordered by `time_field` (a number or ISO-8601 timestamp), and passed through
//! one of two smoothers:
//!
//! * **moving_average** — the mean of a window of `window` points positioned by
//!   `alignment` (backward / centered / forward). At the series ends the window
//!   is clamped to the available points (a partial window), so every feature
//!   still receives a value.
//! * **local_linear** — an adaptive-bandwidth LOESS (degree 1): for each point
//!   the nearest `bandwidth` points in time are fit by tricube-weighted linear
//!   least squares and the fit is evaluated at that point's time. This
//!   reproduces a genuinely linear trend exactly while damping noise.
//!
//! Features whose value or time cannot be parsed are passed through unchanged
//! with a null smoothed value and are excluded from every series.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct TimeSeriesSmoothingTool;

impl Tool for TimeSeriesSmoothingTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "time_series_smoothing",
            display_name: "Time Series Smoothing",
            summary: "Smooth an in-sample time series per location (like ArcGIS Time Series Smoothing): group features by id_field, order each group by time_field, and append a smoothed value via a backward/centered/forward moving average or adaptive-bandwidth local-linear (LOESS) pass. Geometry and all original attributes are preserved. Complements time_series_forecast, which projects ahead rather than smoothing the observed series.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec { name: "input", description: "Input feature layer of timestamped observations.", required: true },
                ToolParamSpec { name: "value_field", description: "Numeric field whose series is smoothed.", required: true },
                ToolParamSpec { name: "time_field", description: "Field ordering each series: a number or an ISO-8601 timestamp.", required: true },
                ToolParamSpec { name: "id_field", description: "Field identifying each location/series. If omitted, all features form one series.", required: false },
                ToolParamSpec { name: "method", description: "'moving_average' (default) or 'local_linear' (adaptive-bandwidth LOESS).", required: false },
                ToolParamSpec { name: "window", description: "moving_average: number of points in the averaging window (default 3).", required: false },
                ToolParamSpec { name: "bandwidth", description: "local_linear: number of nearest-in-time points used per local fit (default 5).", required: false },
                ToolParamSpec { name: "alignment", description: "moving_average window placement: 'backward', 'centered' (default), or 'forward'.", required: false },
                ToolParamSpec { name: "output_field", description: "Name of the appended smoothed field (default 'smoothed').", required: false },
                ToolParamSpec { name: "output", description: "Output layer path. If omitted, stored in memory.", required: false },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "value_field")?;
        require_str(args, "time_field")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let n = layer.features.len();

        let value_idx = layer.schema.field_index(&prm.value_field).ok_or_else(|| {
            ToolError::Validation(format!("value_field '{}' not found", prm.value_field))
        })?;
        let time_idx = layer.schema.field_index(&prm.time_field).ok_or_else(|| {
            ToolError::Validation(format!("time_field '{}' not found", prm.time_field))
        })?;
        let id_idx = match &prm.id_field {
            Some(f) => Some(
                layer
                    .schema
                    .field_index(f)
                    .ok_or_else(|| ToolError::Validation(format!("id_field '{f}' not found")))?,
            ),
            None => None,
        };
        if layer.schema.field_index(&prm.output_field).is_some() {
            return Err(ToolError::Validation(format!(
                "output_field '{}' already exists on the input; choose another name",
                prm.output_field
            )));
        }

        // Read each feature's group key, time, and value. Features missing a
        // parseable time or value are excluded from all series (smoothed = null).
        let mut group_of: Vec<Option<String>> = vec![None; n];
        for (i, feat) in layer.features.iter().enumerate() {
            // Exclude features whose value or time cannot be parsed.
            if feat
                .attributes
                .get(value_idx)
                .and_then(FieldValue::as_f64)
                .is_none()
            {
                continue;
            }
            if feat
                .attributes
                .get(time_idx)
                .and_then(parse_time_value)
                .is_none()
            {
                continue;
            }
            let key = match id_idx {
                Some(gi) => match feat.attributes.get(gi) {
                    Some(v) if !v.is_null() => group_key(v),
                    _ => continue, // no group key -> excluded
                },
                None => String::new(), // single global series
            };
            group_of[i] = Some(key);
        }

        // Bucket usable feature indices by group key.
        let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for (i, g) in group_of.iter().enumerate() {
            if let Some(key) = g {
                groups.entry(key.clone()).or_default().push(i);
            }
        }

        ctx.progress.info(&format!(
            "{n} feature(s) in {} series; method {}",
            groups.len(),
            prm.method.label()
        ));

        // Compute smoothed values into a per-feature buffer (None = pass-through).
        let mut smoothed: Vec<Option<f64>> = vec![None; n];
        for indices in groups.values() {
            // Order this group's features by time.
            let mut ordered: Vec<usize> = indices.clone();
            ordered.sort_by(|&a, &b| {
                let ta = feature_time(&layer, a, time_idx);
                let tb = feature_time(&layer, b, time_idx);
                ta.partial_cmp(&tb).unwrap_or(std::cmp::Ordering::Equal)
            });
            let times: Vec<f64> = ordered
                .iter()
                .map(|&i| feature_time(&layer, i, time_idx))
                .collect();
            let values: Vec<f64> = ordered
                .iter()
                .map(|&i| {
                    layer.features[i].attributes[value_idx]
                        .as_f64()
                        .unwrap_or(0.0)
                })
                .collect();

            let sm = match prm.method {
                Method::MovingAverage => moving_average(&values, prm.window, prm.alignment),
                Method::LocalLinear => local_linear(&times, &values, prm.bandwidth),
            };
            for (pos, &orig) in ordered.iter().enumerate() {
                smoothed[orig] = Some(sm[pos]);
            }
        }

        let n_smoothed = smoothed.iter().filter(|s| s.is_some()).count();

        layer.add_field(FieldDef::new(&prm.output_field, FieldType::Float));
        for (feat, sm) in layer.features.iter_mut().zip(&smoothed) {
            let fv = match sm {
                Some(v) => FieldValue::Float(*v),
                None => FieldValue::Null,
            };
            feat.attributes.push(fv);
        }

        let out_path = write_or_store_layer(layer, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("series_count".to_string(), json!(groups.len()));
        outputs.insert("smoothed_count".to_string(), json!(n_smoothed));
        outputs.insert("method".to_string(), json!(prm.method.label()));
        Ok(ToolRunResult { outputs })
    }
}

// ── smoothers ────────────────────────────────────────────────────────────────

/// Windowed moving average. `window` is the number of points in the window;
/// `alignment` positions the window relative to the current point. Windows are
/// clamped to `[0, m)` at the ends (partial windows), so every point gets a
/// value even for a series shorter than the window.
fn moving_average(values: &[f64], window: usize, alignment: Alignment) -> Vec<f64> {
    let m = values.len();
    let w = window.max(1);
    // Points taken before / after the current index.
    let (before, after) = match alignment {
        Alignment::Backward => (w - 1, 0),
        Alignment::Forward => (0, w - 1),
        Alignment::Centered => ((w - 1) / 2, w / 2),
    };
    (0..m)
        .map(|j| {
            let lo = j.saturating_sub(before);
            let hi = (j + after).min(m - 1);
            let slice = &values[lo..=hi];
            slice.iter().sum::<f64>() / slice.len() as f64
        })
        .collect()
}

/// Adaptive-bandwidth local-linear smoother (LOESS, degree 1). For each point,
/// the nearest `bandwidth` points in time are weighted by a tricube kernel of
/// their scaled time distance and a weighted line is fit and evaluated at that
/// point's time. Reproduces a perfectly linear series exactly.
fn local_linear(times: &[f64], values: &[f64], bandwidth: usize) -> Vec<f64> {
    let m = times.len();
    if m == 0 {
        return Vec::new();
    }
    let k = bandwidth.clamp(2, m.max(2)).min(m);
    (0..m)
        .map(|j| {
            let t0 = times[j];
            // Distances to every point; take the k nearest.
            let mut dists: Vec<(f64, usize)> = (0..m).map(|i| ((times[i] - t0).abs(), i)).collect();
            dists.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            let neigh = &dists[..k];
            let dmax = neigh.last().map(|d| d.0).unwrap_or(0.0);

            // Tricube weights; degenerate (all-coincident) window -> uniform.
            let mut sw = 0.0;
            let mut swx = 0.0;
            let mut swy = 0.0;
            let mut swxx = 0.0;
            let mut swxy = 0.0;
            for &(d, i) in neigh {
                let w = if dmax > 0.0 {
                    let u = (d / dmax).min(1.0);
                    let t = 1.0 - u * u * u;
                    (t * t * t).max(0.0)
                } else {
                    1.0
                };
                let x = times[i];
                let y = values[i];
                sw += w;
                swx += w * x;
                swy += w * y;
                swxx += w * x * x;
                swxy += w * x * y;
            }
            // Weighted least squares: solve for slope b and intercept a.
            let denom = sw * swxx - swx * swx;
            if denom.abs() < 1e-12 {
                // No spread in x (or zero weight) -> weighted mean.
                if sw > 0.0 {
                    swy / sw
                } else {
                    values[j]
                }
            } else {
                let b = (sw * swxy - swx * swy) / denom;
                let a = (swy - b * swx) / sw;
                a + b * t0
            }
        })
        .collect()
}

// ── param plumbing ───────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    MovingAverage,
    LocalLinear,
}

impl Method {
    fn label(&self) -> &'static str {
        match self {
            Method::MovingAverage => "moving_average",
            Method::LocalLinear => "local_linear",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Alignment {
    Backward,
    Centered,
    Forward,
}

struct Params {
    value_field: String,
    time_field: String,
    id_field: Option<String>,
    method: Method,
    window: usize,
    bandwidth: usize,
    alignment: Alignment,
    output_field: String,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let value_field = require_str(args, "value_field")?.to_string();
    let time_field = require_str(args, "time_field")?.to_string();
    let id_field = parse_optional_str(args, "id_field")?.map(String::from);

    let method = match args
        .get("method")
        .and_then(Value::as_str)
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("") | Some("moving_average") => Method::MovingAverage,
        Some("local_linear") => Method::LocalLinear,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'method' must be moving_average or local_linear, got '{o}'"
            )))
        }
    };

    let alignment = match args
        .get("alignment")
        .and_then(Value::as_str)
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("") | Some("centered") => Alignment::Centered,
        Some("backward") => Alignment::Backward,
        Some("forward") => Alignment::Forward,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'alignment' must be backward, centered, or forward, got '{o}'"
            )))
        }
    };

    let window = opt_u64(args, "window")?.unwrap_or(3).max(1) as usize;
    let bandwidth = opt_u64(args, "bandwidth")?.unwrap_or(5).max(2) as usize;
    let output_field = parse_optional_str(args, "output_field")?
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .unwrap_or("smoothed")
        .to_string();

    Ok(Params {
        value_field,
        time_field,
        id_field,
        method,
        window,
        bandwidth,
        alignment,
        output_field,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn opt_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("'{key}' must be an integer"))),
        _ => Ok(None),
    }
}

/// A stable string key for a group value (numbers formatted deterministically).
fn group_key(v: &FieldValue) -> String {
    match v {
        FieldValue::Text(s) | FieldValue::Date(s) | FieldValue::DateTime(s) => s.clone(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) => format!("{f:.9}"),
        FieldValue::Boolean(b) => b.to_string(),
        other => format!("{other:?}"),
    }
}

fn feature_time(layer: &wbvector::Layer, idx: usize, time_idx: usize) -> f64 {
    layer.features[idx]
        .attributes
        .get(time_idx)
        .and_then(parse_time_value)
        .unwrap_or(0.0)
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

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType};
    use wbvector::{Geometry, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Rows: (lng, lat, id, time, value). Builds an in-memory point layer.
    fn layer_of(rows: &[(f64, f64, &str, f64, Option<f64>)]) -> String {
        let mut l = Layer::new("obs")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("gid", FieldType::Text));
        l.add_field(FieldDef::new("t", FieldType::Float));
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (lng, lat, id, t, v) in rows {
            let vv = match v {
                Some(x) => FieldValue::Float(*x),
                None => FieldValue::Null,
            };
            l.add_feature(
                Some(Geometry::point(*lng, *lat)),
                &[
                    ("gid", FieldValue::Text((*id).to_string())),
                    ("t", (*t).into()),
                    ("v", vv),
                ],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn read_field(path: &str, field: &str) -> Vec<Option<f64>> {
        let l = load_input_layer(path).unwrap();
        let i = l.schema.field_index(field).unwrap();
        l.features
            .iter()
            .map(|f| f.attributes[i].as_f64())
            .collect()
    }

    /// Centered moving average of a constant series is the constant everywhere.
    #[test]
    fn moving_average_preserves_constant() {
        let out = moving_average(&[5.0, 5.0, 5.0, 5.0, 5.0], 3, Alignment::Centered);
        for v in out {
            assert!((v - 5.0).abs() < 1e-9);
        }
    }

    /// Centered MA of a linear ramp reproduces it exactly at interior points.
    #[test]
    fn moving_average_linear_interior() {
        let series: Vec<f64> = (0..9).map(|i| 2.0 + 3.0 * i as f64).collect();
        let out = moving_average(&series, 3, Alignment::Centered);
        for i in 1..8 {
            assert!((out[i] - series[i]).abs() < 1e-9, "interior point {i}");
        }
    }

    /// Backward/forward alignment shift the window as expected.
    #[test]
    fn moving_average_alignment() {
        let s = [0.0, 10.0, 20.0, 30.0];
        // backward window=2 at index 1 = mean(0,10)=5.
        let b = moving_average(&s, 2, Alignment::Backward);
        assert!((b[1] - 5.0).abs() < 1e-9);
        // forward window=2 at index 1 = mean(10,20)=15.
        let f = moving_average(&s, 2, Alignment::Forward);
        assert!((f[1] - 15.0).abs() < 1e-9);
    }

    /// Local-linear reproduces a perfectly linear series exactly, even with
    /// unevenly spaced times.
    #[test]
    fn local_linear_reproduces_line() {
        let times = [0.0, 1.0, 2.5, 4.0, 7.0, 9.0];
        let values: Vec<f64> = times.iter().map(|t| -3.0 + 2.5 * t).collect();
        let out = local_linear(&times, &values, 4);
        for i in 0..times.len() {
            assert!((out[i] - values[i]).abs() < 1e-6, "point {i}: {}", out[i]);
        }
    }

    /// Local-linear damps a spike: variance of the smoothed series drops.
    #[test]
    fn local_linear_reduces_variance() {
        let times: Vec<f64> = (0..11).map(|i| i as f64).collect();
        let mut values = vec![0.0; 11];
        values[5] = 10.0; // lone spike on a flat series
        let out = local_linear(&times, &values, 5);
        let var = |x: &[f64]| {
            let m = x.iter().sum::<f64>() / x.len() as f64;
            x.iter().map(|v| (v - m).powi(2)).sum::<f64>()
        };
        assert!(var(&out) < var(&values), "smoothing should reduce variance");
    }

    /// End-to-end: two independent groups are smoothed separately; every
    /// feature gets a value and original order is preserved.
    #[test]
    fn runs_end_to_end_grouped() {
        // Interleave two groups so ordering by time matters.
        let rows = vec![
            (-100.0, 40.0, "a", 0.0, Some(0.0)),
            (-90.0, 30.0, "b", 0.0, Some(100.0)),
            (-100.0, 40.0, "a", 1.0, Some(10.0)),
            (-90.0, 30.0, "b", 1.0, Some(100.0)),
            (-100.0, 40.0, "a", 2.0, Some(20.0)),
            (-90.0, 30.0, "b", 2.0, Some(100.0)),
        ];
        let path = layer_of(&rows);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": path, "value_field": "v", "time_field": "t", "id_field": "gid",
            "method": "moving_average", "window": 3, "alignment": "centered",
        }))
        .unwrap();
        let out = TimeSeriesSmoothingTool.run(&args, &ctx()).unwrap();
        assert_eq!(out.outputs["series_count"], json!(2));
        assert_eq!(out.outputs["smoothed_count"], json!(6));
        let sm = read_field(out.outputs["output"].as_str().unwrap(), "smoothed");
        // Group b is constant 100 -> stays 100.
        assert!((sm[1].unwrap() - 100.0).abs() < 1e-9);
        // Group a linear 0,10,20 centered MA reproduces interior point (10).
        assert!((sm[2].unwrap() - 10.0).abs() < 1e-9);
        // Group a is not contaminated by group b's high values.
        assert!(sm[0].unwrap() < 20.0);
    }

    /// Features with a null value are passed through with a null smoothed value.
    #[test]
    fn null_values_pass_through() {
        let rows = vec![
            (-100.0, 40.0, "a", 0.0, Some(1.0)),
            (-100.0, 40.0, "a", 1.0, None), // unparseable value -> excluded
            (-100.0, 40.0, "a", 2.0, Some(3.0)),
        ];
        let path = layer_of(&rows);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": path, "value_field": "v", "time_field": "t", "id_field": "gid",
        }))
        .unwrap();
        let out = TimeSeriesSmoothingTool.run(&args, &ctx()).unwrap();
        assert_eq!(out.outputs["smoothed_count"], json!(2));
        let sm = read_field(out.outputs["output"].as_str().unwrap(), "smoothed");
        assert!(sm[1].is_none(), "null value stays null");
        assert!(sm[0].is_some() && sm[2].is_some());
    }

    /// Custom output field name is honored.
    #[test]
    fn custom_output_field() {
        let rows = vec![
            (-100.0, 40.0, "a", 0.0, Some(1.0)),
            (-100.0, 40.0, "a", 1.0, Some(2.0)),
        ];
        let path = layer_of(&rows);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": path, "value_field": "v", "time_field": "t",
            "output_field": "v_smooth",
        }))
        .unwrap();
        let out = TimeSeriesSmoothingTool.run(&args, &ctx()).unwrap();
        let l = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        assert!(l.schema.field_index("v_smooth").is_some());
        // No id_field -> one global series.
        assert_eq!(out.outputs["series_count"], json!(1));
    }

    #[test]
    fn rejects_bad_params() {
        let check = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            TimeSeriesSmoothingTool.validate(&args)
        };
        assert!(check(json!({})).is_err());
        assert!(check(json!({ "input": "a.geojson" })).is_err());
        assert!(check(json!({ "input": "a.geojson", "value_field": "v" })).is_err());
        assert!(check(
            json!({ "input": "a.geojson", "value_field": "v", "time_field": "t", "method": "bogus" })
        )
        .is_err());
        assert!(check(json!({
            "input": "a.geojson", "value_field": "v", "time_field": "t", "alignment": "sideways"
        }))
        .is_err());
        assert!(
            check(json!({ "input": "a.geojson", "value_field": "v", "time_field": "t" })).is_ok()
        );
    }
}
