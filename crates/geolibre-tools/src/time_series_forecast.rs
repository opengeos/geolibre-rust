//! GeoLibre tool: per-location forecasting on an H3 space-time cube.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Exponential Smoothing Forecast* and
//! *Curve Fit Forecast* (Space Time Pattern Mining). Completes the space-time
//! cube suite (`emerging_hot_spot_analysis`, `time_series_clustering`,
//! `change_point_detection`, `local_outlier_analysis`) with per-location
//! forecasting. No bundled equivalent exists — `trend_surface` is a spatial
//! polynomial, not temporal, and nothing produces forward predictions with
//! confidence intervals.
//!
//! Timestamped points are binned into an H3 cell × time-step cube. Each cell's
//! series is forecast `steps` ahead by Holt's exponential smoothing and/or
//! polynomial curve fits (linear, parabolic); `model=auto` picks per cell the
//! one with the lowest hold-out RMSE. The output H3 polygon carries the chosen
//! model, the next-step and final-step forecasts, a 90% confidence half-width,
//! and the validation RMSE.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use h3o::{CellIndex, LatLng, Resolution};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct TimeSeriesForecastTool;

impl Tool for TimeSeriesForecastTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "time_series_forecast",
            display_name: "Time Series Forecast",
            summary: "Per-location forecasting on an H3 space-time cube (like ArcGIS Exponential Smoothing / Curve Fit Forecast): forecast each cell's series ahead by Holt's exponential smoothing or polynomial curve fits, picking per cell the lowest hold-out-RMSE model, with next/final-step forecasts and a 90% confidence half-width. The temporal forecasting the bundled spatial trend_surface can't do.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec { name: "input", description: "Input point layer of timestamped events (lon/lat).", required: true },
                ToolParamSpec { name: "time_field", description: "Field holding each point's time: a number or an ISO-8601 timestamp.", required: true },
                ToolParamSpec { name: "output", description: "Output H3 polygon layer with forecast attributes. If omitted, stored in memory.", required: false },
                ToolParamSpec { name: "value_field", description: "Numeric field to aggregate per bin (default: point count).", required: false },
                ToolParamSpec { name: "steps", description: "Number of time steps to forecast ahead (default 3).", required: false },
                ToolParamSpec { name: "model", description: "'auto' (lowest hold-out RMSE; default), 'exp_smoothing', 'linear', or 'parabolic'.", required: false },
                ToolParamSpec { name: "holdout", description: "Time steps held out to validate/select the model (default 3).", required: false },
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
        if n_times < prm.holdout + 3 {
            return Err(ToolError::Execution(format!(
                "need at least {} time steps (holdout {} + 3); reduce time_step or holdout",
                prm.holdout + 3,
                prm.holdout
            )));
        }
        let time_bin = |t: f64| (((t - t_min) / prm.time_step).floor() as usize).min(n_times - 1);

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
            "{n_cells} cell(s) x {n_times} step(s); forecasting {} step(s)",
            prm.steps
        ));

        let mut out = Layer::new("forecast")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        out.add_field(FieldDef::new("h3", FieldType::Text));
        out.add_field(FieldDef::new("model", FieldType::Text));
        out.add_field(FieldDef::new("forecast_1", FieldType::Float));
        out.add_field(FieldDef::new("forecast_n", FieldType::Float));
        out.add_field(FieldDef::new("ci90", FieldType::Float));
        out.add_field(FieldDef::new("rmse", FieldType::Float));

        for ci in 0..n_cells {
            let series = &cube[ci * n_times..ci * n_times + n_times];
            let fc = forecast_cell(series, &prm);
            out.add_feature(
                Some(Geometry::polygon(cell_polygon_ring(cells[ci]), Vec::new())),
                &[
                    ("h3", cells[ci].to_string().into()),
                    ("model", fc.model.into()),
                    ("forecast_1", fc.first.into()),
                    ("forecast_n", fc.last.into()),
                    ("ci90", fc.ci90.into()),
                    ("rmse", fc.rmse.into()),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("cell feature failed: {e}")))?;
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("cell_count".to_string(), json!(n_cells));
        outputs.insert("time_steps".to_string(), json!(n_times));
        outputs.insert("forecast_steps".to_string(), json!(prm.steps));
        Ok(ToolRunResult { outputs })
    }
}

struct CellForecast {
    model: &'static str,
    first: f64,
    last: f64,
    ci90: f64,
    rmse: f64,
}

/// Selects a model (by hold-out RMSE for `auto`), refits on the full series, and
/// forecasts `steps` ahead.
fn forecast_cell(series: &[f64], prm: &Params) -> CellForecast {
    let n = series.len();
    let train = &series[..n - prm.holdout];
    let test = &series[n - prm.holdout..];

    let candidates: &[Model] = match &prm.model {
        Some(m) => std::slice::from_ref(m),
        None => &[Model::ExpSmoothing, Model::Linear, Model::Parabolic],
    };

    let mut best = (Model::Linear, f64::INFINITY);
    for &m in candidates {
        let pred = m.forecast(train, prm.holdout);
        let rmse = rmse(&pred, test);
        if rmse < best.1 {
            best = (m, rmse);
        }
    }
    let (model, rmse) = best;

    // Refit on the full series and forecast `steps` ahead.
    let full = model.forecast(series, prm.steps);
    // Residual-based 90% CI half-width (1.645 sigma).
    let resid = model.residual_std(series);
    CellForecast {
        model: model.name(),
        first: *full.first().unwrap_or(&0.0),
        last: *full.last().unwrap_or(&0.0),
        ci90: 1.645 * resid,
        rmse,
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Model {
    ExpSmoothing,
    Linear,
    Parabolic,
}

impl Model {
    fn name(&self) -> &'static str {
        match self {
            Model::ExpSmoothing => "exp_smoothing",
            Model::Linear => "linear",
            Model::Parabolic => "parabolic",
        }
    }

    /// Forecast `h` steps beyond the series.
    fn forecast(&self, series: &[f64], h: usize) -> Vec<f64> {
        match self {
            Model::ExpSmoothing => holt(series, h),
            Model::Linear => poly_forecast(series, 1, h),
            Model::Parabolic => poly_forecast(series, 2, h),
        }
    }

    /// Std of in-sample one-step residuals.
    fn residual_std(&self, series: &[f64]) -> f64 {
        let n = series.len();
        if n < 4 {
            return 0.0;
        }
        // One-step-ahead in-sample errors from a refit on a growing window is
        // costly; use the fit residuals of the whole series instead.
        let fit = match self {
            Model::ExpSmoothing => holt_fit(series),
            Model::Linear => poly_fit(series, 1),
            Model::Parabolic => poly_fit(series, 2),
        };
        let sse: f64 = series.iter().zip(&fit).map(|(y, f)| (y - f).powi(2)).sum();
        (sse / (n as f64 - 1.0)).sqrt()
    }
}

/// Holt's linear exponential smoothing; grid-searches alpha/beta, forecasts `h`.
fn holt(series: &[f64], h: usize) -> Vec<f64> {
    let (alpha, beta) = best_holt_params(series);
    let (level, trend) = holt_state(series, alpha, beta);
    (1..=h).map(|k| level + k as f64 * trend).collect()
}

#[allow(clippy::needless_range_loop)]
fn holt_fit(series: &[f64]) -> Vec<f64> {
    let (alpha, beta) = best_holt_params(series);
    let n = series.len();
    let mut level = series[0];
    let mut trend = series[1] - series[0];
    let mut fit = vec![series[0]; n];
    for t in 1..n {
        let prev_level = level;
        fit[t] = level + trend; // one-step-ahead forecast
        level = alpha * series[t] + (1.0 - alpha) * (level + trend);
        trend = beta * (level - prev_level) + (1.0 - beta) * trend;
    }
    fit
}

#[allow(clippy::needless_range_loop)]
fn holt_state(series: &[f64], alpha: f64, beta: f64) -> (f64, f64) {
    let n = series.len();
    let mut level = series[0];
    let mut trend = series[1] - series[0];
    for t in 1..n {
        let prev_level = level;
        level = alpha * series[t] + (1.0 - alpha) * (level + trend);
        trend = beta * (level - prev_level) + (1.0 - beta) * trend;
    }
    (level, trend)
}

fn best_holt_params(series: &[f64]) -> (f64, f64) {
    let mut best = (0.3, 0.1);
    let mut best_sse = f64::INFINITY;
    for ai in 1..=9 {
        for bi in 0..=9 {
            let (a, b) = (ai as f64 / 10.0, bi as f64 / 10.0);
            let sse = holt_insample_sse(series, a, b);
            if sse < best_sse {
                best_sse = sse;
                best = (a, b);
            }
        }
    }
    best
}

#[allow(clippy::needless_range_loop)]
fn holt_insample_sse(series: &[f64], alpha: f64, beta: f64) -> f64 {
    let n = series.len();
    let mut level = series[0];
    let mut trend = series[1] - series[0];
    let mut sse = 0.0;
    for t in 1..n {
        let pred = level + trend;
        sse += (series[t] - pred).powi(2);
        let prev_level = level;
        level = alpha * series[t] + (1.0 - alpha) * (level + trend);
        trend = beta * (level - prev_level) + (1.0 - beta) * trend;
    }
    sse
}

/// Polynomial (degree 1 or 2) least-squares forecast `h` steps ahead.
fn poly_forecast(series: &[f64], degree: usize, h: usize) -> Vec<f64> {
    let coeffs = poly_coeffs(series, degree);
    let n = series.len();
    (0..h)
        .map(|k| {
            let x = (n + k) as f64;
            eval_poly(&coeffs, x)
        })
        .collect()
}

fn poly_fit(series: &[f64], degree: usize) -> Vec<f64> {
    let coeffs = poly_coeffs(series, degree);
    (0..series.len())
        .map(|i| eval_poly(&coeffs, i as f64))
        .collect()
}

/// Least-squares polynomial coefficients (ascending powers) over x = 0..n.
fn poly_coeffs(series: &[f64], degree: usize) -> Vec<f64> {
    let m = degree + 1;
    let n = series.len();
    // Normal equations A^T A c = A^T y.
    let mut ata = vec![vec![0.0f64; m]; m];
    let mut aty = vec![0.0f64; m];
    for (i, &y) in series.iter().enumerate() {
        let x = i as f64;
        let powers: Vec<f64> = (0..m).map(|p| x.powi(p as i32)).collect();
        for r in 0..m {
            aty[r] += powers[r] * y;
            for c in 0..m {
                ata[r][c] += powers[r] * powers[c];
            }
        }
    }
    solve(&mut ata, &mut aty).unwrap_or_else(|| {
        // Fallback: constant = mean.
        let mut v = vec![0.0; m];
        v[0] = series.iter().sum::<f64>() / n as f64;
        v
    })
}

fn eval_poly(coeffs: &[f64], x: f64) -> f64 {
    coeffs
        .iter()
        .enumerate()
        .map(|(p, c)| c * x.powi(p as i32))
        .sum()
}

#[allow(clippy::needless_range_loop)]
fn solve(a: &mut [Vec<f64>], b: &mut [f64]) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        let mut piv = col;
        for r in (col + 1)..n {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-12 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        for r in (col + 1)..n {
            let f = a[r][col] / a[col][col];
            for c in col..n {
                a[r][c] -= f * a[col][c];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let mut s = b[r];
        for c in (r + 1)..n {
            s -= a[r][c] * x[c];
        }
        x[r] = s / a[r][r];
    }
    Some(x)
}

fn rmse(pred: &[f64], actual: &[f64]) -> f64 {
    let n = pred.len().min(actual.len());
    if n == 0 {
        return f64::INFINITY;
    }
    let sse: f64 = (0..n).map(|i| (pred[i] - actual[i]).powi(2)).sum();
    (sse / n as f64).sqrt()
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
    steps: usize,
    model: Option<Model>,
    holdout: usize,
    time_step: f64,
    resolution: Resolution,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let time_field = require_str(args, "time_field")?.to_string();
    let value_field = parse_optional_str(args, "value_field")?.map(String::from);
    let steps = opt_u64(args, "steps")?.unwrap_or(3).max(1) as usize;
    let model = match args.get("model").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("auto") => None,
        Some("exp_smoothing") => Some(Model::ExpSmoothing),
        Some("linear") => Some(Model::Linear),
        Some("parabolic") => Some(Model::Parabolic),
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'model' must be auto/exp_smoothing/linear/parabolic, got '{o}'"
            )))
        }
    };
    let holdout = opt_u64(args, "holdout")?.unwrap_or(3).max(1) as usize;
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
        steps,
        model,
        holdout,
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

    /// A perfectly linear series is forecast (near-)exactly by the linear model.
    #[test]
    fn linear_forecast_is_accurate() {
        // y = 2 + 3t for t=0..9.
        let series: Vec<f64> = (0..10).map(|t| 2.0 + 3.0 * t as f64).collect();
        let f = poly_forecast(&series, 1, 3);
        // Next values at t=10,11,12 -> 32, 35, 38.
        assert!((f[0] - 32.0).abs() < 1e-6);
        assert!((f[2] - 38.0).abs() < 1e-6);
    }

    /// Holt tracks a linear trend forward.
    #[test]
    fn holt_extrapolates_trend() {
        let series: Vec<f64> = (0..12).map(|t| 5.0 + 2.0 * t as f64).collect();
        let f = holt(&series, 2);
        // Should be increasing and well above the last observed value (27).
        assert!(
            f[0] > 27.0 && f[1] > f[0],
            "Holt should extend the upward trend"
        );
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

    /// End-to-end: forecast attributes written per cell; auto picks a model.
    #[test]
    fn runs_end_to_end() {
        let mut rows = Vec::new();
        for t in 0..12 {
            let v = 1.0 + 0.5 * t as f64;
            for _ in 0..2 {
                rows.push((-100.0, 40.0, t as f64, v));
            }
        }
        let args: ToolArgs = serde_json::from_value(json!({
            "input": pts(&rows), "time_field": "t", "value_field": "v", "resolution": 6,
            "steps": 3, "holdout": 3,
        }))
        .unwrap();
        let out = TimeSeriesForecastTool.run(&args, &ctx()).unwrap();
        let l = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        assert!(l.schema.field_index("forecast_1").is_some());
        assert!(l.schema.field_index("model").is_some());
        assert_eq!(out.outputs["forecast_steps"], json!(3));
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            TimeSeriesForecastTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "time_field": "t", "model": "arima" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "time_field": "t" })).is_ok());
    }
}
