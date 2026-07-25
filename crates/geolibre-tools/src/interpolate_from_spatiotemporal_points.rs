//! GeoLibre tool: time-binned interpolation of point observations into a
//! multidimensional raster.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Interpolate From Spatiotemporal
//! Points* (Image Analyst). Every interpolator available today operates on a
//! single snapshot — the bundled `idw_interpolation`,
//! `natural_neighbour_interpolation`, `thin_plate_spline`,
//! `radial_basis_function_interpolation` and the kriging family, plus the
//! shipped `interpolate_with_barriers` / `diffusion_interpolation_with_barriers`
//! / `optimal_interpolation`. Producing a time series of surfaces from a sensor
//! archive means splitting the points by hand, running the interpolator per
//! slice, and reassembling.
//!
//! That reassembled stack is exactly what the shipped time-series suite
//! (`generate_trend_raster`, `multidimensional_anomaly`, `analyze_changes_ccdc`,
//! `time_series_smoothing`) consumes, so this is the bridge from raw
//! observations to all of them.
//!
//! The correctness requirement that makes the output usable downstream is a
//! **single shared grid**: every slice is built on one extent and cell size, so
//! the slices are co-registered and a per-pixel temporal model is meaningful. A
//! per-slice extent derived from that slice's own points would drift and make
//! the stack useless.

use std::collections::BTreeMap;

use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster, RasterConfig};
use wbvector::{FieldValue, Geometry};

use crate::common::{parse_optional_output, write_or_store_output};
use crate::vector_common::{load_input_layer, parse_optional_str};

const OUT_NODATA: f64 = -9999.0;

/// Per-slice interpolation method.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Method {
    Idw,
    Nearest,
    Mean,
    Median,
}

/// Calendar-aware time binning.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum TimeStep {
    Daily,
    Weekly,
    Monthly,
    Quarterly,
    Yearly,
}

pub struct InterpolateFromSpatiotemporalPointsTool;

impl Tool for InterpolateFromSpatiotemporalPointsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "interpolate_from_spatiotemporal_points",
            display_name: "Interpolate From Spatiotemporal Points",
            summary: "Bin timestamped point observations into regular time slices and interpolate each onto one shared grid, producing a co-registered multidimensional raster, like ArcGIS Interpolate From Spatiotemporal Points.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Timestamped point features.",
                    required: true,
                },
                ToolParamSpec {
                    name: "value_field",
                    description: "Numeric field to interpolate.",
                    required: true,
                },
                ToolParamSpec {
                    name: "time_field",
                    description: "Timestamp field: numeric seconds since epoch, or ISO-8601 text.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output raster path for the first slice. Later slices are written alongside with a _t<index> suffix. If omitted, all slices are stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "time_step",
                    description: "daily | weekly | monthly (default) | quarterly | yearly.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Output cell size in map units. Defaults to 1/50th of the larger extent dimension.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "idw (default) | nearest | mean | median.",
                    required: false,
                },
                ToolParamSpec {
                    name: "power",
                    description: "IDW distance exponent (default 2).",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighbors",
                    description: "Number of nearest observations used per cell (default 8).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_points",
                    description: "Minimum observations required for a slice to be interpolated; sparser slices emit all-nodata rather than a misleading surface (default 3).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for k in ["input", "value_field", "time_field"] {
            require_str(args, k)?;
        }
        parse_method(args)?;
        parse_time_step(args)?;
        for (k, lo) in [("cell_size", 0.0), ("power", 0.0)] {
            if let Some(v) = parse_optional_f64(args, k)? {
                if !v.is_finite() || v <= lo {
                    return Err(ToolError::Validation(format!(
                        "'{k}' must be greater than {lo}"
                    )));
                }
            }
        }
        for k in ["neighbors", "min_points"] {
            if let Some(v) = parse_optional_f64(args, k)? {
                if !v.is_finite() || v < 1.0 {
                    return Err(ToolError::Validation(format!("'{k}' must be at least 1")));
                }
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let value_field = require_str(args, "value_field")?;
        let time_field = require_str(args, "time_field")?;
        let output = parse_optional_output(args, "output")?;
        let method = parse_method(args)?;
        let step = parse_time_step(args)?;
        let power = parse_optional_f64(args, "power")?.unwrap_or(2.0);
        let neighbors = parse_optional_f64(args, "neighbors")?.unwrap_or(8.0) as usize;
        let min_points = parse_optional_f64(args, "min_points")?.unwrap_or(3.0) as usize;

        let layer = load_input_layer(input)?;
        for f in [value_field, time_field] {
            if layer.schema.field_index(f).is_none() {
                return Err(ToolError::Validation(format!(
                    "field '{f}' not found on the input layer"
                )));
            }
        }

        // Extract (x, y, value, time).
        let mut obs: Vec<(f64, f64, f64, f64)> = Vec::new();
        for feat in layer.features.iter() {
            let (Some(geom), Ok(vv), Ok(tv)) = (
                feat.geometry.as_ref(),
                feat.get(&layer.schema, value_field),
                feat.get(&layer.schema, time_field),
            ) else {
                continue;
            };
            let (Some((x, y)), Some(v), Some(t)) =
                (point_xy(geom), vv.as_f64(), parse_time_value(&tv))
            else {
                continue;
            };
            if v.is_finite() && t.is_finite() {
                obs.push((x, y, v, t));
            }
        }
        if obs.is_empty() {
            return Err(ToolError::Validation(
                "no usable timestamped point observations found".to_string(),
            ));
        }

        // ONE shared grid over ALL observations, so every slice is co-registered.
        let (min_x, min_y, max_x, max_y) = obs.iter().fold(
            (f64::MAX, f64::MAX, f64::MIN, f64::MIN),
            |(a, b, c, d), (x, y, _, _)| (a.min(*x), b.min(*y), c.max(*x), d.max(*y)),
        );
        let span = (max_x - min_x).max(max_y - min_y).max(f64::EPSILON);
        let cell = parse_optional_f64(args, "cell_size")?.unwrap_or(span / 50.0);
        let cols = (((max_x - min_x) / cell).ceil() as usize + 1).max(1);
        let rows = (((max_y - min_y) / cell).ceil() as usize + 1).max(1);
        if cols * rows > 50_000_000 {
            return Err(ToolError::Validation(format!(
                "cell_size {cell} would produce a {rows}x{cols} grid; supply a larger cell_size"
            )));
        }

        // Bin by calendar-aware slice key.
        let mut slices: BTreeMap<i64, Vec<(f64, f64, f64)>> = BTreeMap::new();
        for (x, y, v, t) in &obs {
            slices
                .entry(bin_key(*t, step))
                .or_default()
                .push((*x, *y, *v));
        }

        ctx.progress.info(&format!(
            "{} observation(s) -> {} slice(s) on a shared {rows}x{cols} grid",
            obs.len(),
            slices.len()
        ));

        let mut written: Vec<String> = Vec::new();
        let mut interpolated = 0usize;
        let mut sparse = 0usize;

        for (si, (key, pts)) in slices.iter().enumerate() {
            let mut data = vec![OUT_NODATA; rows * cols];

            if pts.len() >= min_points {
                let mut tree: KdTree<f64, usize, [f64; 2]> = KdTree::new(2);
                for (i, (x, y, _)) in pts.iter().enumerate() {
                    tree.add([*x, *y], i)
                        .map_err(|e| ToolError::Execution(format!("kd-tree insert: {e:?}")))?;
                }
                let k = neighbors.min(pts.len());
                for r in 0..rows {
                    for c in 0..cols {
                        // Cell centre in map coordinates. Rows run north to
                        // south, matching the raster convention.
                        let px = min_x + (c as f64 + 0.5) * cell;
                        let py = max_y - (r as f64 + 0.5) * cell;
                        let hits = tree
                            .nearest(&[px, py], k, &squared_euclidean)
                            .map_err(|e| ToolError::Execution(format!("kd-tree query: {e:?}")))?;
                        if hits.is_empty() {
                            continue;
                        }
                        let vals: Vec<(f64, f64)> =
                            hits.iter().map(|(d2, &i)| (d2.sqrt(), pts[i].2)).collect();
                        data[r * cols + c] = combine(&vals, method, power);
                    }
                }
                interpolated += 1;
            } else {
                // Too few observations to interpolate honestly: leave the slice
                // empty rather than emit a surface implying information.
                sparse += 1;
            }

            let raster = build_raster(rows, cols, min_x, max_y, cell, data)?;
            let target = match (output, si) {
                (Some(p), 0) => Some(p.to_string()),
                (Some(p), _) => Some(suffixed(p, si)),
                (None, _) => None,
            };
            written.push(write_or_store_output(raster, target.as_deref())?);
            ctx.progress
                .progress((si as f64 + 1.0) / slices.len().max(1) as f64);
            let _ = key;
        }

        let keys: Vec<i64> = slices.keys().copied().collect();
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(written[0]));
        outputs.insert("outputs".to_string(), json!(written));
        outputs.insert("slice_count".to_string(), json!(written.len()));
        outputs.insert("slice_keys".to_string(), json!(keys));
        outputs.insert("interpolated_slices".to_string(), json!(interpolated));
        // Surfaced so an all-nodata slice is visibly "too sparse", not a bug.
        outputs.insert("sparse_slices".to_string(), json!(sparse));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        outputs.insert("cell_size".to_string(), json!(cell));
        Ok(ToolRunResult { outputs })
    }
}

/// Combines the k nearest observations under the chosen method.
fn combine(vals: &[(f64, f64)], method: Method, power: f64) -> f64 {
    match method {
        Method::Nearest => vals
            .iter()
            .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(_, v)| *v)
            .unwrap_or(OUT_NODATA),
        Method::Mean => vals.iter().map(|(_, v)| v).sum::<f64>() / vals.len() as f64,
        Method::Median => {
            let mut vs: Vec<f64> = vals.iter().map(|(_, v)| *v).collect();
            vs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            vs[vs.len() / 2]
        }
        Method::Idw => {
            // An observation sitting exactly on the cell centre would divide by
            // zero; return it directly instead.
            if let Some((_, v)) = vals.iter().find(|(d, _)| *d < 1e-12) {
                return *v;
            }
            let mut num = 0.0;
            let mut den = 0.0;
            for (d, v) in vals {
                let w = 1.0 / d.powf(power);
                num += w * v;
                den += w;
            }
            if den > 0.0 {
                num / den
            } else {
                OUT_NODATA
            }
        }
    }
}

fn build_raster(
    rows: usize,
    cols: usize,
    min_x: f64,
    max_y: f64,
    cell: f64,
    data: Vec<f64>,
) -> Result<Raster, ToolError> {
    let mut r = Raster::new(RasterConfig {
        cols,
        rows,
        bands: 1,
        x_min: min_x,
        y_min: max_y - rows as f64 * cell,
        cell_size: cell,
        cell_size_y: Some(cell),
        nodata: OUT_NODATA,
        data_type: DataType::F32,
        crs: wbraster::CrsInfo::default(),
        metadata: Default::default(),
    });
    for row in 0..rows {
        for col in 0..cols {
            r.set(0, row as isize, col as isize, data[row * cols + col])
                .map_err(|e| ToolError::Execution(format!("failed writing cell: {e}")))?;
        }
    }
    Ok(r)
}

/// Calendar-aware bin key. Monthly/quarterly/yearly use real calendar
/// arithmetic rather than fixed-day approximations, which would drift across
/// month lengths and leap years.
fn bin_key(t: f64, step: TimeStep) -> i64 {
    let secs = t as i64;
    let days = secs.div_euclid(86400);
    match step {
        TimeStep::Daily => days,
        TimeStep::Weekly => days.div_euclid(7),
        _ => {
            let (y, m, _) = civil_from_days(days);
            match step {
                TimeStep::Monthly => y * 12 + (m - 1),
                TimeStep::Quarterly => y * 4 + (m - 1) / 3,
                TimeStep::Yearly => y,
                _ => unreachable!(),
            }
        }
    }
}

/// Inverse of days_from_civil (Howard Hinnant's algorithm).
fn civil_from_days(z: i64) -> (i64, i64, i64) {
    let z = z + 719468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    (if m <= 2 { y + 1 } else { y }, m, d)
}

fn suffixed(path: &str, i: usize) -> String {
    match path.rfind('.') {
        Some(dot) if !path[dot..].contains('/') && !path[dot..].contains('\\') => {
            format!("{}_t{}{}", &path[..dot], i, &path[dot..])
        }
        _ => format!("{path}_t{i}"),
    }
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

// ── parameter parsing ────────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

fn parse_method(args: &ToolArgs) -> Result<Method, ToolError> {
    match args.get("method").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("idw") => Ok(Method::Idw),
        Some("nearest") => Ok(Method::Nearest),
        Some("mean") => Ok(Method::Mean),
        Some("median") => Ok(Method::Median),
        Some(o) => Err(ToolError::Validation(format!(
            "'method' must be idw/nearest/mean/median, got '{o}'"
        ))),
    }
}

fn parse_time_step(args: &ToolArgs) -> Result<TimeStep, ToolError> {
    match args.get("time_step").and_then(Value::as_str).map(str::trim) {
        Some("daily") => Ok(TimeStep::Daily),
        Some("weekly") => Ok(TimeStep::Weekly),
        None | Some("") | Some("monthly") => Ok(TimeStep::Monthly),
        Some("quarterly") => Ok(TimeStep::Quarterly),
        Some("yearly") => Ok(TimeStep::Yearly),
        Some(o) => Err(ToolError::Validation(format!(
            "'time_step' must be daily/weekly/monthly/quarterly/yearly, got '{o}'"
        ))),
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Coord, FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    const DAY: f64 = 86400.0;

    /// Observations as (x, y, value, time-seconds).
    fn obs_layer(items: Vec<(f64, f64, f64, f64)>) -> String {
        let mut l = Layer::new("o")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        l.add_field(FieldDef::new("t", FieldType::Float));
        for (x, y, v, t) in items {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(x, y))),
                &[("v", FieldValue::Float(v)), ("t", FieldValue::Float(t))],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let mut v = args;
        v["value_field"] = json!("v");
        v["time_field"] = json!("t");
        let args: ToolArgs = serde_json::from_value(v).unwrap();
        InterpolateFromSpatiotemporalPointsTool
            .run(&args, &ctx())
            .unwrap()
    }

    fn read(path: &str) -> Raster {
        crate::common::load_input_raster(path).unwrap()
    }

    /// A four-corner grid of observations in two months yields two slices.
    fn two_month_obs() -> Vec<(f64, f64, f64, f64)> {
        let mut v = Vec::new();
        // January: value 10 everywhere.
        for (x, y) in [(0.0, 0.0), (100.0, 0.0), (0.0, 100.0), (100.0, 100.0)] {
            v.push((x, y, 10.0, 0.0));
        }
        // February (about 40 days later): value 50 everywhere.
        for (x, y) in [(0.0, 0.0), (100.0, 0.0), (0.0, 100.0), (100.0, 100.0)] {
            v.push((x, y, 50.0, 40.0 * DAY));
        }
        v
    }

    /// THE property that makes the stack usable downstream: every slice shares
    /// one grid, so a per-pixel temporal model is meaningful.
    #[test]
    fn all_slices_share_one_grid() {
        let out = run(json!({ "input": obs_layer(two_month_obs()), "cell_size": 25.0 }));
        let paths = out.outputs["outputs"].as_array().unwrap();
        assert_eq!(paths.len(), 2);
        let a = read(paths[0].as_str().unwrap());
        let b = read(paths[1].as_str().unwrap());
        assert_eq!((a.rows, a.cols), (b.rows, b.cols));
        assert!((a.x_min - b.x_min).abs() < 1e-9);
        assert!((a.y_min - b.y_min).abs() < 1e-9);
        assert!((a.cell_size_x - b.cell_size_x).abs() < 1e-9);
    }

    /// Each slice carries its own month's values.
    #[test]
    fn slices_carry_their_own_values() {
        let out = run(json!({ "input": obs_layer(two_month_obs()), "cell_size": 50.0 }));
        let paths = out.outputs["outputs"].as_array().unwrap();
        let a = read(paths[0].as_str().unwrap());
        let b = read(paths[1].as_str().unwrap());
        // All four corners are 10 in slice 0, so every cell interpolates to 10.
        assert!(
            (a.get(0, 0, 0) - 10.0).abs() < 1e-6,
            "got {}",
            a.get(0, 0, 0)
        );
        assert!(
            (b.get(0, 0, 0) - 50.0).abs() < 1e-6,
            "got {}",
            b.get(0, 0, 0)
        );
    }

    /// Calendar binning is real, not a fixed-day approximation: Jan 31 and
    /// Feb 1 are one day apart but land in different monthly bins.
    #[test]
    fn monthly_binning_is_calendar_aware() {
        let jan31 = days_from_civil(2026, 1, 31) as f64 * DAY;
        let feb01 = days_from_civil(2026, 2, 1) as f64 * DAY;
        assert_eq!(
            (feb01 - jan31) / DAY,
            1.0,
            "the two timestamps are one day apart"
        );
        assert_ne!(
            bin_key(jan31, TimeStep::Monthly),
            bin_key(feb01, TimeStep::Monthly),
            "consecutive days across a month boundary must bin separately"
        );
        // And the same two days share a bin under yearly.
        assert_eq!(
            bin_key(jan31, TimeStep::Yearly),
            bin_key(feb01, TimeStep::Yearly)
        );
    }

    /// Quarterly binning groups three months.
    #[test]
    fn quarterly_binning_groups_three_months() {
        let jan = days_from_civil(2026, 1, 15) as f64 * DAY;
        let mar = days_from_civil(2026, 3, 15) as f64 * DAY;
        let apr = days_from_civil(2026, 4, 15) as f64 * DAY;
        assert_eq!(
            bin_key(jan, TimeStep::Quarterly),
            bin_key(mar, TimeStep::Quarterly)
        );
        assert_ne!(
            bin_key(mar, TimeStep::Quarterly),
            bin_key(apr, TimeStep::Quarterly)
        );
    }

    /// A slice below min_points emits nodata rather than a misleading surface.
    #[test]
    fn sparse_slice_is_left_empty_and_counted() {
        let mut items = two_month_obs();
        // A third month with a single observation.
        items.push((50.0, 50.0, 999.0, 100.0 * DAY));
        let out = run(json!({
            "input": obs_layer(items), "cell_size": 50.0, "min_points": 3
        }));
        assert_eq!(out.outputs["slice_count"], json!(3));
        assert_eq!(out.outputs["interpolated_slices"], json!(2));
        assert_eq!(out.outputs["sparse_slices"], json!(1));
        let paths = out.outputs["outputs"].as_array().unwrap();
        let sparse = read(paths[2].as_str().unwrap());
        assert_eq!(sparse.get(0, 0, 0), OUT_NODATA);
    }

    /// IDW reproduces an observation sitting on a cell centre exactly, rather
    /// than dividing by a zero distance.
    #[test]
    fn idw_handles_coincident_observation() {
        let vals = [(0.0, 42.0), (10.0, 1.0)];
        assert!((combine(&vals, Method::Idw, 2.0) - 42.0).abs() < 1e-12);
    }

    /// The non-IDW methods behave as named.
    #[test]
    fn alternative_methods() {
        let vals = [(1.0, 10.0), (2.0, 20.0), (3.0, 60.0)];
        assert_eq!(combine(&vals, Method::Nearest, 2.0), 10.0);
        assert_eq!(combine(&vals, Method::Mean, 2.0), 30.0);
        assert_eq!(combine(&vals, Method::Median, 2.0), 20.0);
    }

    /// ISO-8601 timestamps are accepted alongside numeric seconds.
    #[test]
    fn accepts_iso8601_timestamps() {
        let mut l = Layer::new("o")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        l.add_field(FieldDef::new("t", FieldType::Text));
        for (x, y, t) in [
            (0.0, 0.0, "2026-01-05"),
            (10.0, 0.0, "2026-01-06"),
            (0.0, 10.0, "2026-01-07"),
            (0.0, 0.0, "2026-03-05"),
            (10.0, 0.0, "2026-03-06"),
            (0.0, 10.0, "2026-03-07"),
        ] {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(x, y))),
                &[
                    ("v", FieldValue::Float(1.0)),
                    ("t", FieldValue::Text(t.into())),
                ],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        let out = run(json!({
            "input": memory_store::make_vector_memory_path(&id), "cell_size": 5.0
        }));
        assert_eq!(out.outputs["slice_count"], json!(2), "January and March");
    }

    /// The civil-date round trip underpinning calendar binning.
    #[test]
    fn civil_date_round_trip() {
        for (y, m, d) in [(2026, 1, 1), (2026, 2, 28), (2024, 2, 29), (1999, 12, 31)] {
            let days = days_from_civil(y, m, d);
            assert_eq!(civil_from_days(days), (y, m, d), "round trip {y}-{m}-{d}");
        }
    }

    #[test]
    fn rejects_bad_parameters() {
        let p = obs_layer(two_month_obs());
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            InterpolateFromSpatiotemporalPointsTool
                .validate(&args)
                .is_err()
        };
        assert!(bad(json!({ "value_field": "v", "time_field": "t" })));
        assert!(bad(json!({ "input": p, "time_field": "t" })));
        assert!(bad(
            json!({ "input": p, "value_field": "v", "time_field": "t",
                           "method": "kriging" })
        ));
        assert!(bad(
            json!({ "input": p, "value_field": "v", "time_field": "t",
                           "time_step": "hourly" })
        ));
        assert!(bad(
            json!({ "input": p, "value_field": "v", "time_field": "t",
                           "cell_size": 0 })
        ));
        assert!(bad(
            json!({ "input": p, "value_field": "v", "time_field": "t",
                           "neighbors": 0 })
        ));
    }
}
