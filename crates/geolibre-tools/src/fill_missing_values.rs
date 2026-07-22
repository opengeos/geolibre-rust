//! GeoLibre tool: fill missing (null) attribute values from spatial and/or
//! temporal neighbours.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Fill Missing Values*
//! (Space Time Pattern Mining). For every feature whose `fill_field` is null
//! (SQL `NULL` or a non-finite float), the tool gathers neighbouring features
//! that DO have a value and estimates the missing one with a chosen estimator:
//!
//! * `mean` (default), `median`, `min`, `max` — a summary of the neighbours'
//!   values.
//! * `temporal_trend` — a least-squares line fitted to the neighbours'
//!   (time, value) pairs, evaluated at the target feature's own time. Requires
//!   `time_field`; falls back to the neighbour mean when fewer than two distinct
//!   times are available.
//!
//! Neighbours are selected spatially, either as the `k` nearest features
//! (`neighbourhood = knn`, the default) or every feature within
//! `search_radius` (`neighbourhood = distance_band`). When a `time_field` and
//! `time_window` are both supplied, candidate neighbours are additionally
//! restricted to those whose time is within `time_window` of the target — i.e.
//! spatial *and* temporal neighbours.
//!
//! Filled rows are marked with a boolean flag field (default `imputed`) so the
//! imputed values can be told apart downstream. The estimator is deterministic
//! (no RNG); neighbour search is an O(n²) brute-force scan, matching the sizes
//! this crate's other distribution tools target.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Imputes null values in a numeric attribute field from spatial (and optional
/// temporal) neighbours.
pub struct FillMissingValuesTool;

impl Tool for FillMissingValuesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "fill_missing_values",
            display_name: "Fill Missing Values",
            summary: "Estimate and fill null values in a numeric attribute field using spatial neighbours (k-nearest or a distance band), optionally restricted to a temporal window, with a selectable estimator (mean/median/min/max/temporal_trend). Imputed rows are marked with a flag field (ArcGIS Fill Missing Values).",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer with the field to fill.",
                    required: true,
                },
                ToolParamSpec {
                    name: "fill_field",
                    description: "Numeric field whose null values are to be estimated and filled.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "estimator",
                    description: "How to combine neighbour values: 'mean' (default), 'median', 'min', 'max', or 'temporal_trend' (least-squares over time; needs time_field).",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighbourhood",
                    description: "Spatial selection: 'knn' (default, k nearest) or 'distance_band' (all within search_radius).",
                    required: false,
                },
                ToolParamSpec {
                    name: "k",
                    description: "Number of nearest neighbours for 'knn' (default 8).",
                    required: false,
                },
                ToolParamSpec {
                    name: "search_radius",
                    description: "Radius (CRS units) for 'distance_band' selection.",
                    required: false,
                },
                ToolParamSpec {
                    name: "time_field",
                    description: "Optional numeric field giving each feature's time/order (required for 'temporal_trend').",
                    required: false,
                },
                ToolParamSpec {
                    name: "time_window",
                    description: "Optional: restrict neighbours to those within this time distance of the target (needs time_field).",
                    required: false,
                },
                ToolParamSpec {
                    name: "flag_field",
                    description: "Name of the boolean field flagging imputed rows (default 'imputed').",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "fill_field")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let fill_field = require_str(args, "fill_field")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;

        let fill_idx = layer
            .schema
            .field_index(fill_field)
            .ok_or_else(|| ToolError::Validation(format!("fill_field '{fill_field}' not found")))?;
        let fill_is_int = matches!(
            layer.schema.fields()[fill_idx].field_type,
            FieldType::Integer
        );

        let time_idx = match &prm.time_field {
            Some(f) => Some(
                layer
                    .schema
                    .field_index(f)
                    .ok_or_else(|| ToolError::Validation(format!("time_field '{f}' not found")))?,
            ),
            None => None,
        };
        if matches!(prm.estimator, Estimator::TemporalTrend) && time_idx.is_none() {
            return Err(ToolError::Validation(
                "estimator 'temporal_trend' requires 'time_field'".to_string(),
            ));
        }

        // Representative point + value (+ time) for every feature.
        struct Rec {
            xy: Option<(f64, f64)>,
            value: Option<f64>, // finite value, or None if missing
            time: Option<f64>,
        }
        let recs: Vec<Rec> = layer
            .features
            .iter()
            .map(|f| {
                let xy = f.geometry.as_ref().and_then(representative_xy);
                let value = finite_value(&f.attributes, fill_idx);
                let time = time_idx.and_then(|ti| finite_value(&f.attributes, ti));
                Rec { xy, value, time }
            })
            .collect();

        // Source pool: features with a valid value and a location.
        let sources: Vec<usize> = recs
            .iter()
            .enumerate()
            .filter(|(_, r)| r.value.is_some() && r.xy.is_some())
            .map(|(i, _)| i)
            .collect();

        // Build output layer: same schema + a flag field.
        let mut out = layer.clone();
        out.name = "fill_missing_values".to_string();
        let flag_idx = out
            .schema
            .upsert_field(FieldDef::new(prm.flag_field.clone(), FieldType::Integer));
        // Pad existing rows so every feature has a slot for the flag field.
        for feat in out.features.iter_mut() {
            if feat.attributes.len() <= flag_idx {
                feat.attributes.resize(flag_idx + 1, FieldValue::Null);
            }
            feat.attributes[flag_idx] = FieldValue::Integer(0);
        }

        let total = recs.len();
        let mut filled = 0usize;
        let mut still_null = 0usize;

        for i in 0..total {
            // Only impute features that are missing but geolocated.
            if recs[i].value.is_some() {
                continue;
            }
            let Some(target_xy) = recs[i].xy else {
                still_null += 1;
                continue;
            };
            let target_time = recs[i].time;

            // Candidate neighbours: sources satisfying the optional time window.
            let mut cands: Vec<(f64, f64, f64)> = Vec::new(); // (dist, value, time)
            for &s in &sources {
                if s == i {
                    continue;
                }
                let sxy = recs[s].xy.unwrap();
                let sval = recs[s].value.unwrap();
                if let (Some(win), Some(t0), Some(ts)) =
                    (prm.time_window, target_time, recs[s].time)
                {
                    if (ts - t0).abs() > win {
                        continue;
                    }
                }
                let d = (target_xy.0 - sxy.0).hypot(target_xy.1 - sxy.1);
                let t = recs[s].time.unwrap_or(0.0);
                cands.push((d, sval, t));
            }

            // Spatial narrowing.
            let neigh: Vec<(f64, f64)> = match prm.neighbourhood {
                Neighbourhood::Knn => {
                    cands.sort_by(|a, b| a.0.total_cmp(&b.0));
                    cands.iter().take(prm.k).map(|c| (c.1, c.2)).collect()
                }
                Neighbourhood::DistanceBand => cands
                    .iter()
                    .filter(|c| c.0 <= prm.search_radius)
                    .map(|c| (c.1, c.2))
                    .collect(),
            };

            if neigh.is_empty() {
                still_null += 1;
                continue;
            }

            let estimate = match prm.estimator {
                Estimator::Mean => mean(neigh.iter().map(|n| n.0)),
                Estimator::Median => median(neigh.iter().map(|n| n.0).collect()),
                Estimator::Min => neigh.iter().map(|n| n.0).fold(f64::INFINITY, f64::min),
                Estimator::Max => neigh.iter().map(|n| n.0).fold(f64::NEG_INFINITY, f64::max),
                Estimator::TemporalTrend => {
                    let t0 = target_time.unwrap_or(0.0);
                    temporal_trend(&neigh, t0)
                }
            };

            if !estimate.is_finite() {
                still_null += 1;
                continue;
            }

            let fv = if fill_is_int {
                FieldValue::Integer(estimate.round() as i64)
            } else {
                FieldValue::Float(estimate)
            };
            out.features[i].attributes[fill_idx] = fv;
            out.features[i].attributes[flag_idx] = FieldValue::Integer(1);
            filled += 1;
        }

        ctx.progress.info(&format!(
            "{filled} value(s) imputed from {} source feature(s); {still_null} still null",
            sources.len()
        ));

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("total_features".to_string(), json!(total));
        outputs.insert("source_count".to_string(), json!(sources.len()));
        outputs.insert("filled_count".to_string(), json!(filled));
        outputs.insert("remaining_null_count".to_string(), json!(still_null));
        Ok(ToolRunResult { outputs })
    }
}

// ── Estimators ─────────────────────────────────────────────────────────────────

fn mean(vals: impl Iterator<Item = f64>) -> f64 {
    let mut sum = 0.0;
    let mut n = 0.0;
    for v in vals {
        sum += v;
        n += 1.0;
    }
    if n == 0.0 {
        f64::NAN
    } else {
        sum / n
    }
}

fn median(mut vals: Vec<f64>) -> f64 {
    if vals.is_empty() {
        return f64::NAN;
    }
    vals.sort_by(f64::total_cmp);
    let n = vals.len();
    if n % 2 == 1 {
        vals[n / 2]
    } else {
        (vals[n / 2 - 1] + vals[n / 2]) * 0.5
    }
}

/// Least-squares line `value ≈ a + b·time` fitted to the neighbour
/// (value, time) pairs, evaluated at `t0`. Falls back to the neighbour mean
/// when fewer than two distinct times exist (slope undefined).
fn temporal_trend(neigh: &[(f64, f64)], t0: f64) -> f64 {
    let n = neigh.len() as f64;
    let mean_t = neigh.iter().map(|(_, t)| *t).sum::<f64>() / n;
    let mean_v = neigh.iter().map(|(v, _)| *v).sum::<f64>() / n;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for (v, t) in neigh {
        let dt = t - mean_t;
        sxx += dt * dt;
        sxy += dt * (v - mean_v);
    }
    if sxx <= f64::EPSILON {
        return mean_v; // all neighbours share one time -> use their mean
    }
    let slope = sxy / sxx;
    let intercept = mean_v - slope * mean_t;
    intercept + slope * t0
}

// ── Helpers ────────────────────────────────────────────────────────────────────

/// A field value counts as "present" only when it is a finite number.
fn finite_value(attrs: &[FieldValue], idx: usize) -> Option<f64> {
    attrs
        .get(idx)
        .and_then(FieldValue::as_f64)
        .filter(|v| v.is_finite())
}

fn representative_xy(geom: &Geometry) -> Option<(f64, f64)> {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0u64;
    accumulate(geom, &mut sx, &mut sy, &mut n);
    (n > 0).then(|| (sx / n as f64, sy / n as f64))
}

fn accumulate(geom: &Geometry, sx: &mut f64, sy: &mut f64, n: &mut u64) {
    let mut add = |c: &Coord| {
        *sx += c.x;
        *sy += c.y;
        *n += 1;
    };
    match geom {
        Geometry::Point(c) => add(c),
        Geometry::LineString(cs) | Geometry::MultiPoint(cs) => cs.iter().for_each(add),
        Geometry::MultiLineString(lines) => lines.iter().flatten().for_each(add),
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            exterior.coords().iter().for_each(&mut add);
            interiors
                .iter()
                .for_each(|r| r.coords().iter().for_each(&mut add));
        }
        Geometry::MultiPolygon(polys) => {
            for (ext, holes) in polys {
                ext.coords().iter().for_each(&mut add);
                holes
                    .iter()
                    .for_each(|r| r.coords().iter().for_each(&mut add));
            }
        }
        Geometry::GeometryCollection(geoms) => {
            for g in geoms {
                accumulate(g, sx, sy, n);
            }
        }
    }
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

// ── Parameters ──────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Estimator {
    Mean,
    Median,
    Min,
    Max,
    TemporalTrend,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Neighbourhood {
    Knn,
    DistanceBand,
}

struct Params {
    estimator: Estimator,
    neighbourhood: Neighbourhood,
    k: usize,
    search_radius: f64,
    time_field: Option<String>,
    time_window: Option<f64>,
    flag_field: String,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let estimator = match parse_optional_str(args, "estimator")? {
        None => Estimator::Mean,
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "mean" => Estimator::Mean,
            "median" => Estimator::Median,
            "min" => Estimator::Min,
            "max" => Estimator::Max,
            "temporal_trend" | "temporal-trend" | "trend" => Estimator::TemporalTrend,
            other => {
                return Err(ToolError::Validation(format!(
                    "'estimator' must be one of mean/median/min/max/temporal_trend, got '{other}'"
                )))
            }
        },
    };
    let neighbourhood = match parse_optional_str(args, "neighbourhood")? {
        None => Neighbourhood::Knn,
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "knn" | "k_nearest" | "k-nearest" => Neighbourhood::Knn,
            "distance_band" | "distance-band" | "band" => Neighbourhood::DistanceBand,
            other => {
                return Err(ToolError::Validation(format!(
                    "'neighbourhood' must be 'knn' or 'distance_band', got '{other}'"
                )))
            }
        },
    };
    let k = match parse_optional_f64(args, "k")? {
        None => 8,
        Some(v) if v >= 1.0 => v as usize,
        Some(v) => {
            return Err(ToolError::Validation(format!(
                "'k' must be a positive integer, got {v}"
            )))
        }
    };
    let search_radius = parse_optional_f64(args, "search_radius")?.unwrap_or(0.0);
    if neighbourhood == Neighbourhood::DistanceBand && search_radius <= 0.0 {
        return Err(ToolError::Validation(
            "'distance_band' neighbourhood requires a positive 'search_radius'".to_string(),
        ));
    }
    let time_window = match parse_optional_f64(args, "time_window")? {
        None => None,
        Some(v) if v > 0.0 => Some(v),
        Some(v) => {
            return Err(ToolError::Validation(format!(
                "'time_window' must be positive, got {v}"
            )))
        }
    };
    let time_field = parse_optional_str(args, "time_field")?.map(str::to_string);
    if time_window.is_some() && time_field.is_none() {
        return Err(ToolError::Validation(
            "'time_window' requires 'time_field'".to_string(),
        ));
    }
    if estimator == Estimator::TemporalTrend && time_field.is_none() {
        return Err(ToolError::Validation(
            "estimator 'temporal_trend' requires 'time_field'".to_string(),
        ));
    }
    let flag_field = parse_optional_str(args, "flag_field")?
        .unwrap_or("imputed")
        .to_string();

    Ok(Params {
        estimator,
        neighbourhood,
        k,
        search_radius,
        time_field,
        time_window,
        flag_field,
    })
}

/// Parses an optional numeric parameter accepting a JSON number OR a numeric
/// string (host UIs post everything as strings).
fn parse_optional_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s.trim().parse::<f64>().map(Some).map_err(|_| {
            ToolError::Validation(format!("parameter '{key}' must be a number, got '{s}'"))
        }),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a number"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Build a point layer with a numeric `val` field (None -> Null) and an
    /// optional `t` (time) field.
    fn point_layer(rows: &[(f64, f64, Option<f64>, Option<f64>)], with_time: bool) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("val", FieldType::Float));
        if with_time {
            l.add_field(FieldDef::new("t", FieldType::Float));
        }
        let vidx = l.schema.field_index("val").unwrap();
        let tidx = l.schema.field_index("t");
        for (x, y, v, t) in rows {
            let geom = Some(Geometry::point(*x, *y));
            let mut feat = wbvector::Feature {
                fid: 0,
                geometry: geom,
                attributes: vec![FieldValue::Null; l.schema.len()],
            };
            feat.attributes[vidx] = match v {
                Some(v) => FieldValue::Float(*v),
                None => FieldValue::Null,
            };
            if let (Some(ti), Some(t)) = (tidx, t) {
                feat.attributes[ti] = FieldValue::Float(*t);
            }
            l.push(feat);
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = FillMissingValuesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn val(layer: &Layer, i: usize) -> Option<f64> {
        let idx = layer.schema.field_index("val").unwrap();
        layer.features[i].attributes[idx].as_f64()
    }
    fn flag(layer: &Layer, i: usize) -> i64 {
        let idx = layer.schema.field_index("imputed").unwrap();
        layer.features[i].attributes[idx].as_i64().unwrap()
    }

    /// A null value surrounded by identical neighbours is filled with that value
    /// (mean of equal neighbours) and flagged; untouched rows keep flag 0.
    #[test]
    fn fills_null_with_neighbour_mean() {
        // Index 2 (center) is null; neighbours all = 10.
        let input = point_layer(
            &[
                (0.0, 0.0, Some(10.0), None),
                (1.0, 0.0, Some(10.0), None),
                (0.5, 0.5, None, None),
                (0.0, 1.0, Some(10.0), None),
                (1.0, 1.0, Some(10.0), None),
            ],
            false,
        );
        let (out, layer) = run(json!({ "input": input, "fill_field": "val" }));
        assert_eq!(out.outputs["filled_count"], json!(1));
        assert_eq!(out.outputs["remaining_null_count"], json!(0));
        assert!((val(&layer, 2).unwrap() - 10.0).abs() < 1e-9);
        assert_eq!(flag(&layer, 2), 1);
        // A non-null feature is passed through unchanged and unflagged.
        assert_eq!(flag(&layer, 0), 0);
        assert!((val(&layer, 0).unwrap() - 10.0).abs() < 1e-9);
    }

    /// median/min/max estimators pick the expected neighbour statistic.
    #[test]
    fn estimator_variants() {
        // Center null; neighbour values 1,2,3,100 (so mean=26.5, median=2.5).
        let rows = [
            (0.0, 0.0, Some(1.0), None),
            (1.0, 0.0, Some(2.0), None),
            (0.0, 1.0, Some(3.0), None),
            (1.0, 1.0, Some(100.0), None),
            (0.5, 0.5, None, None),
        ];
        let input = point_layer(&rows, false);
        let (_o, med) = run(json!({ "input": input, "fill_field": "val", "estimator": "median" }));
        assert!((val(&med, 4).unwrap() - 2.5).abs() < 1e-9);

        let input = point_layer(&rows, false);
        let (_o, mn) = run(json!({ "input": input, "fill_field": "val", "estimator": "min" }));
        assert!((val(&mn, 4).unwrap() - 1.0).abs() < 1e-9);

        let input = point_layer(&rows, false);
        let (_o, mx) = run(json!({ "input": input, "fill_field": "val", "estimator": "max" }));
        assert!((val(&mx, 4).unwrap() - 100.0).abs() < 1e-9);
    }

    /// temporal_trend extrapolates a linear time series at the target's time.
    #[test]
    fn temporal_trend_predicts_line() {
        // All at the same location; value = 2*t exactly, target at t=5 is null.
        let input = point_layer(
            &[
                (0.0, 0.0, Some(0.0), Some(0.0)),
                (0.0, 0.0, Some(2.0), Some(1.0)),
                (0.0, 0.0, Some(4.0), Some(2.0)),
                (0.0, 0.0, Some(6.0), Some(3.0)),
                (0.0, 0.0, None, Some(5.0)),
            ],
            true,
        );
        let (_o, layer) = run(json!({
            "input": input, "fill_field": "val",
            "estimator": "temporal_trend", "time_field": "t", "k": 4,
        }));
        // Line value=2*t predicts 10 at t=5.
        assert!((val(&layer, 4).unwrap() - 10.0).abs() < 1e-6);
    }

    /// A null feature with no eligible neighbours (empty distance band) stays
    /// null and unflagged, and is counted in remaining_null_count.
    #[test]
    fn no_neighbours_passes_through() {
        let input = point_layer(
            &[
                (0.0, 0.0, Some(5.0), None),
                (1000.0, 1000.0, None, None), // far outside a radius of 10
            ],
            false,
        );
        let (out, layer) = run(json!({
            "input": input, "fill_field": "val",
            "neighbourhood": "distance_band", "search_radius": 10.0,
        }));
        assert_eq!(out.outputs["filled_count"], json!(0));
        assert_eq!(out.outputs["remaining_null_count"], json!(1));
        assert!(val(&layer, 1).is_none());
        assert_eq!(flag(&layer, 1), 0);
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            FillMissingValuesTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err()); // no fill_field
        assert!(
            bad(json!({ "input": "a.geojson", "fill_field": "v", "estimator": "bogus" })).is_err()
        );
        assert!(bad(
            json!({ "input": "a.geojson", "fill_field": "v", "neighbourhood": "distance_band" })
        )
        .is_err()); // no radius
        assert!(bad(
            json!({ "input": "a.geojson", "fill_field": "v", "estimator": "temporal_trend" })
        )
        .is_err()); // no time_field
        assert!(
            bad(json!({ "input": "a.geojson", "fill_field": "v", "time_window": 5.0 })).is_err()
        ); // window w/o time_field
        assert!(bad(json!({ "input": "a.geojson", "fill_field": "v" })).is_ok());
    }
}
