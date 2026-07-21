//! GeoLibre tool: change-point detection on the time series of a space-time cube.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Change Point Detection* (Space Time
//! Pattern Mining). The repo already builds H3 space-time cubes for
//! `emerging_hot_spot_analysis`, `time_series_clustering`, and
//! `local_outlier_analysis`; this adds temporal segmentation. Nothing bundled
//! performs it — `change_vector_analysis` compares exactly two dates and the
//! Mann-Kendall trend tools detect monotonic trends, not abrupt regime changes.
//!
//! Timestamped points are binned into an H3 cell × time-step cube. Each cell's
//! series is segmented by **binary segmentation**: the split that most reduces
//! the segmentation cost (a shift in `mean` or in linear `slope`) is accepted
//! while its gain exceeds a BIC-style penalty scaled by `sensitivity`, or until
//! `num_change_points` splits are made (`method=defined`). Each H3 cell is output
//! with its change-point count, the year of its largest change, and the segment
//! means on either side.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use h3o::{CellIndex, LatLng, Resolution};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct ChangePointDetectionTool;

impl Tool for ChangePointDetectionTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "change_point_detection",
            display_name: "Change Point Detection",
            summary: "Detect abrupt shifts in the time series at each location of an H3 space-time cube (like ArcGIS Change Point Detection): binary segmentation on a mean- or slope-shift cost with a sensitivity-scaled penalty (or a fixed number of change points), reporting per cell the change count, the largest change's year, and the segment means. The temporal segmentation the two-date change_vector_analysis and monotonic Mann-Kendall tests can't do.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec { name: "input", description: "Input point layer of timestamped events (lon/lat).", required: true },
                ToolParamSpec { name: "time_field", description: "Field holding each point's time: a number or an ISO-8601 timestamp.", required: true },
                ToolParamSpec { name: "output", description: "Output H3 polygon layer with change-point attributes. If omitted, stored in memory.", required: false },
                ToolParamSpec { name: "value_field", description: "Numeric field to aggregate per bin (default: point count).", required: false },
                ToolParamSpec { name: "change_type", description: "'mean' (level shifts; default) or 'slope' (trend shifts).", required: false },
                ToolParamSpec { name: "method", description: "'auto' (penalty-based, default) or 'defined' (exactly num_change_points).", required: false },
                ToolParamSpec { name: "num_change_points", description: "defined method: number of change points to find per location (default 1).", required: false },
                ToolParamSpec { name: "sensitivity", description: "auto method: penalty scale; higher = fewer change points (default 1.0).", required: false },
                ToolParamSpec { name: "time_step", description: "Width of a time step in the time_field's units (default 1).", required: false },
                ToolParamSpec { name: "resolution", description: "H3 resolution 0..15 (default 7).", required: false },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "time_field")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let time_idx = layer.schema.field_index(&prm.time_field).ok_or_else(|| {
            ToolError::Validation(format!("time_field '{}' not found", prm.time_field))
        })?;
        let value_idx =
            match &prm.value_field {
                Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                    ToolError::Validation(format!("value_field '{f}' not found"))
                })?),
                None => None,
            };

        let mut obs: Vec<(CellIndex, f64, f64)> = Vec::new();
        for feature in layer.iter() {
            let Some((lng, lat)) = feature.geometry.as_ref().and_then(point_lnglat) else {
                continue;
            };
            let Some(time) = feature.attributes.get(time_idx).and_then(parse_time_value) else {
                continue;
            };
            let value = match value_idx {
                Some(vi) => match feature.attributes.get(vi).and_then(FieldValue::as_f64) {
                    Some(v) => v,
                    None => continue,
                },
                None => 1.0,
            };
            if let Ok(ll) = LatLng::new(lat, lng) {
                obs.push((ll.to_cell(prm.resolution), time, value));
            }
        }
        if obs.is_empty() {
            return Err(ToolError::Execution("no usable observations".to_string()));
        }

        let t_min = obs.iter().map(|o| o.1).fold(f64::INFINITY, f64::min);
        let t_max = obs.iter().map(|o| o.1).fold(f64::NEG_INFINITY, f64::max);
        let n_times = (((t_max - t_min) / prm.time_step).floor() as usize) + 1;
        if n_times < 4 {
            return Err(ToolError::Execution(
                "need >= 4 time steps for change-point detection (reduce time_step)".to_string(),
            ));
        }
        let time_bin = |t: f64| (((t - t_min) / prm.time_step).floor() as usize).min(n_times - 1);
        let step_year = |ti: usize| t_min + ti as f64 * prm.time_step;

        let cells: Vec<CellIndex> = {
            let set: BTreeSet<u64> = obs.iter().map(|o| u64::from(o.0)).collect();
            set.into_iter()
                .map(|r| CellIndex::try_from(r).unwrap())
                .collect()
        };
        let cell_pos: HashMap<CellIndex, usize> =
            cells.iter().enumerate().map(|(i, &c)| (c, i)).collect();
        let n_cells = cells.len();
        let mut cube = vec![0.0f64; n_cells * n_times];
        for &(cell, time, value) in &obs {
            cube[cell_pos[&cell] * n_times + time_bin(time)] += value;
        }

        ctx.progress.info(&format!(
            "{n_cells} cell(s) x {n_times} step(s); detecting change points"
        ));

        let mut out = Layer::new("change_points")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        out.add_field(FieldDef::new("h3", FieldType::Text));
        out.add_field(FieldDef::new("n_changes", FieldType::Integer));
        out.add_field(FieldDef::new("top_year", FieldType::Float));
        out.add_field(FieldDef::new("mean_before", FieldType::Float));
        out.add_field(FieldDef::new("mean_after", FieldType::Float));

        let mut total_changes = 0i64;
        for ci in 0..n_cells {
            let series = &cube[ci * n_times..ci * n_times + n_times];
            let cps = detect(series, &prm);
            total_changes += cps.len() as i64;

            // Largest change: the split with the biggest mean difference.
            let (top_year, mb, ma) = if let Some(&cp) = cps
                .iter()
                .max_by(|&&a, &&b| shift_size(series, a).total_cmp(&shift_size(series, b)))
            {
                let mb = mean(&series[..cp]);
                let ma = mean(&series[cp..]);
                (step_year(cp), mb, ma)
            } else {
                (f64::NAN, mean(series), mean(series))
            };

            out.add_feature(
                Some(Geometry::polygon(cell_polygon_ring(cells[ci]), Vec::new())),
                &[
                    ("h3", cells[ci].to_string().into()),
                    ("n_changes", (cps.len() as i64).into()),
                    ("top_year", top_year.into()),
                    ("mean_before", mb.into()),
                    ("mean_after", ma.into()),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("cell feature failed: {e}")))?;
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("cell_count".to_string(), json!(n_cells));
        outputs.insert("time_steps".to_string(), json!(n_times));
        outputs.insert("total_change_points".to_string(), json!(total_changes));
        Ok(ToolRunResult { outputs })
    }
}

/// Mean difference across a candidate split at index `cp`.
fn shift_size(series: &[f64], cp: usize) -> f64 {
    if cp == 0 || cp >= series.len() {
        return 0.0;
    }
    (mean(&series[..cp]) - mean(&series[cp..])).abs()
}

/// Binary segmentation returning sorted change-point indices (start of the new
/// segment).
fn detect(series: &[f64], prm: &Params) -> Vec<usize> {
    let n = series.len();
    let total_var = variance(series).max(1e-12);
    let penalty = prm.sensitivity * total_var * (n as f64).ln();
    let mut cps: Vec<usize> = Vec::new();
    // Work list of (start, end) half-open segments.
    let mut stack = vec![(0usize, n)];
    while let Some((a, b)) = stack.pop() {
        if b - a < 4 {
            continue;
        }
        if prm.defined && cps.len() >= prm.num_change_points {
            break;
        }
        if let Some((split, gain)) = best_split(&series[a..b], prm.slope) {
            let accept = if prm.defined {
                gain > 0.0
            } else {
                gain > penalty
            };
            if accept {
                let cp = a + split;
                cps.push(cp);
                stack.push((a, cp));
                stack.push((cp, b));
            }
        }
        if prm.defined && cps.len() >= prm.num_change_points {
            break;
        }
    }
    cps.sort_unstable();
    if prm.defined && cps.len() > prm.num_change_points {
        // Keep the largest-shift change points.
        cps.sort_by(|&x, &y| shift_size(series, y).total_cmp(&shift_size(series, x)));
        cps.truncate(prm.num_change_points);
        cps.sort_unstable();
    }
    cps
}

/// Best single split of a segment: the index maximising cost reduction (SSE for
/// mean, residual SSE for slope). Returns (split index within the segment, gain).
fn best_split(seg: &[f64], slope: bool) -> Option<(usize, f64)> {
    let n = seg.len();
    if n < 4 {
        return None;
    }
    let base = if slope { slope_sse(seg) } else { sse(seg) };
    let mut best = None;
    let mut best_gain = 0.0;
    for k in 2..(n - 1) {
        let cost = if slope {
            slope_sse(&seg[..k]) + slope_sse(&seg[k..])
        } else {
            sse(&seg[..k]) + sse(&seg[k..])
        };
        let gain = base - cost;
        if gain > best_gain {
            best_gain = gain;
            best = Some(k);
        }
    }
    best.map(|k| (k, best_gain))
}

fn mean(v: &[f64]) -> f64 {
    if v.is_empty() {
        0.0
    } else {
        v.iter().sum::<f64>() / v.len() as f64
    }
}

fn sse(v: &[f64]) -> f64 {
    let m = mean(v);
    v.iter().map(|x| (x - m).powi(2)).sum()
}

fn variance(v: &[f64]) -> f64 {
    if v.len() < 2 {
        0.0
    } else {
        sse(v) / v.len() as f64
    }
}

/// Residual SSE of a least-squares line fit over indices 0..n.
fn slope_sse(v: &[f64]) -> f64 {
    let n = v.len();
    if n < 2 {
        return 0.0;
    }
    let nf = n as f64;
    let mx = (nf - 1.0) / 2.0;
    let my = mean(v);
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for (i, &y) in v.iter().enumerate() {
        let dx = i as f64 - mx;
        sxx += dx * dx;
        sxy += dx * (y - my);
    }
    if sxx <= 0.0 {
        return sse(v);
    }
    let slope = sxy / sxx;
    let intercept = my - slope * mx;
    v.iter()
        .enumerate()
        .map(|(i, &y)| {
            let fit = intercept + slope * i as f64;
            (y - fit).powi(2)
        })
        .sum()
}

// ── H3 / time helpers ───────────────────────────────────────────────────────

fn cell_polygon_ring(cell: CellIndex) -> Vec<Coord> {
    cell.boundary()
        .iter()
        .map(|ll| Coord::xy(ll.lng(), ll.lat()))
        .collect()
}

fn point_lnglat(geom: &Geometry) -> Option<(f64, f64)> {
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

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

struct Params {
    time_field: String,
    value_field: Option<String>,
    slope: bool,
    defined: bool,
    num_change_points: usize,
    sensitivity: f64,
    time_step: f64,
    resolution: Resolution,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let time_field = require_str(args, "time_field")?.to_string();
    let value_field = parse_optional_str(args, "value_field")?.map(String::from);
    let slope = match args
        .get("change_type")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some("mean") => false,
        Some("slope") => true,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'change_type' must be 'mean' or 'slope', got '{o}'"
            )))
        }
    };
    let defined = match args.get("method").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("auto") => false,
        Some("defined") => true,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'method' must be 'auto' or 'defined', got '{o}'"
            )))
        }
    };
    let num_change_points = opt_u64(args, "num_change_points")?.unwrap_or(1).max(1) as usize;
    let sensitivity = opt_f64(args, "sensitivity")?.unwrap_or(1.0).max(0.0);
    let time_step = opt_f64(args, "time_step")?.unwrap_or(1.0);
    if time_step <= 0.0 {
        return Err(ToolError::Validation("'time_step' must be positive".into()));
    }
    let res_num = opt_u64(args, "resolution")?.unwrap_or(7);
    let resolution = Resolution::try_from(res_num as u8)
        .map_err(|_| ToolError::Validation("'resolution' must be 0..15".into()))?;
    Ok(Params {
        time_field,
        value_field,
        slope,
        defined,
        num_change_points,
        sensitivity,
        time_step,
        resolution,
    })
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
            .map_err(|_| ToolError::Validation(format!("'{key}' must be a number"))),
        _ => Ok(None),
    }
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

    /// A clear level shift is found at the right index.
    #[test]
    fn finds_a_mean_shift() {
        // low for 6 steps, then high for 6 steps -> change at index 6.
        let series = [1.0, 1.1, 0.9, 1.0, 1.1, 0.9, 5.0, 5.1, 4.9, 5.0, 5.1, 4.9];
        let prm = Params {
            time_field: "t".into(),
            value_field: None,
            slope: false,
            defined: false,
            num_change_points: 1,
            sensitivity: 1.0,
            time_step: 1.0,
            resolution: Resolution::Seven,
        };
        let cps = detect(&series, &prm);
        assert!(!cps.is_empty(), "should find at least one change point");
        assert!(
            cps.iter().any(|&c| (c as i64 - 6).abs() <= 1),
            "change near index 6, got {cps:?}"
        );
    }

    /// A flat series has no change points.
    #[test]
    fn flat_series_has_no_changes() {
        let series = [2.0; 12];
        let prm = Params {
            time_field: "t".into(),
            value_field: None,
            slope: false,
            defined: false,
            num_change_points: 1,
            sensitivity: 1.0,
            time_step: 1.0,
            resolution: Resolution::Seven,
        };
        assert!(detect(&series, &prm).is_empty());
    }

    fn pts(rows: &[(f64, f64, f64, f64)]) -> String {
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("t", FieldType::Float));
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (lng, lat, t, v) in rows {
            l.add_feature(
                Some(Geometry::point(*lng, *lat)),
                &[("t", (*t).into()), ("v", (*v).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    /// End-to-end: cube built, per-cell change attributes written.
    #[test]
    fn runs_end_to_end() {
        let mut rows = Vec::new();
        for t in 0..12 {
            let v = if t < 6 { 1.0 } else { 8.0 };
            for _ in 0..3 {
                rows.push((-100.0, 40.0, t as f64, v));
            }
        }
        let args: ToolArgs = serde_json::from_value(json!({
            "input": pts(&rows), "time_field": "t", "value_field": "v", "resolution": 6,
        }))
        .unwrap();
        let out = ChangePointDetectionTool.run(&args, &ctx()).unwrap();
        assert_eq!(out.outputs["time_steps"], json!(12));
        assert!(out.outputs["total_change_points"].as_i64().unwrap() >= 1);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ChangePointDetectionTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err());
        assert!(
            bad(json!({ "input": "a.geojson", "time_field": "t", "change_type": "var" })).is_err()
        );
        assert!(bad(json!({ "input": "a.geojson", "time_field": "t" })).is_ok());
    }
}
