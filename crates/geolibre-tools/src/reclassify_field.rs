//! GeoLibre tool: classify a numeric vector attribute into class bins.
//!
//! The vector-attribute twin of the raster `slice`/`reclass` family: read one
//! numeric `field`, compute class breaks with a chosen classification `method`
//! (equal interval, defined interval, quantile, natural breaks / Jenks,
//! standard deviation, or geometric interval), and write a 1-based integer
//! `class_field` back onto every feature — the standard choropleth-preparation
//! step. Optionally a `break_field` records the upper class limit each feature
//! falls under, so the classification is self-documenting and feeds directly
//! into the repo's `color_polygons` / `render_vector_png` rendering pipeline.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Reclassify Field* (Data Management).
//! Where `transform_fields` (`bin`) offers a subset of these methods as one of
//! many attribute transforms, this tool is the dedicated classifier: it adds
//! natural breaks (Jenks), defined interval, and geometric interval, returns
//! the break vector, and keeps the class column the single, named output.
//!
//! All math is deterministic (no RNG) and null-aware: features whose value is
//! null / non-finite receive a null class rather than poisoning the break fit.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct ReclassifyFieldTool;

impl Tool for ReclassifyFieldTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "reclassify_field",
            display_name: "Reclassify Field",
            summary: "Classify a numeric vector attribute into class bins (like ArcGIS Reclassify Field): equal_interval, defined_interval, quantile, natural_breaks (Jenks), std_dev, or geometric_interval, writing an integer class field plus an optional break-value field.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (or table) with the numeric field to classify.",
                    required: true,
                },
                ToolParamSpec {
                    name: "field",
                    description: "Name of the numeric field to classify.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector layer with the class field added. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "Classification method: equal_interval, defined_interval, quantile, natural_breaks, std_dev, or geometric_interval (default natural_breaks).",
                    required: false,
                },
                ToolParamSpec {
                    name: "classes",
                    description: "Number of classes for equal_interval/quantile/natural_breaks/geometric_interval (default 5). Ignored by defined_interval and std_dev, whose class count follows the data.",
                    required: false,
                },
                ToolParamSpec {
                    name: "interval",
                    description: "Fixed class width for defined_interval (required for that method).",
                    required: false,
                },
                ToolParamSpec {
                    name: "std_dev_interval",
                    description: "Interval size for std_dev, in units of standard deviation (default 1.0; e.g. 0.5 for half-sigma bands).",
                    required: false,
                },
                ToolParamSpec {
                    name: "class_field",
                    description: "Name of the integer class field to add (default 'class').",
                    required: false,
                },
                ToolParamSpec {
                    name: "break_field",
                    description: "If given, also add a float field with the upper class limit each feature falls under.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "field")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let n = layer.len();

        let idx = layer
            .schema
            .field_index(&prm.field)
            .ok_or_else(|| ToolError::Validation(format!("field '{}' not found", prm.field)))?;

        // Read the target column (null-aware).
        let raw: Vec<Option<f64>> = (0..n)
            .map(|fi| {
                layer.features[fi]
                    .attributes
                    .get(idx)
                    .and_then(FieldValue::as_f64)
                    .filter(|v| v.is_finite())
            })
            .collect();
        let vals: Vec<f64> = raw.iter().filter_map(|v| *v).collect();
        let n_valid = vals.len();
        if n_valid == 0 {
            return Err(ToolError::Execution(format!(
                "field '{}' has no finite numeric values to classify",
                prm.field
            )));
        }

        ctx.progress.info(&format!(
            "reclassify_field: '{}' method={} over {n_valid} valid value(s)",
            prm.field,
            prm.method.as_str()
        ));

        // Internal breaks (k-1 for a fixed-k method, variable otherwise).
        let breaks = compute_breaks(&vals, &prm)?;
        let n_classes = breaks.len() + 1;
        // Upper class limits: internal breaks followed by the data maximum.
        let (_lo, hi) = min_max(&vals);
        let mut uppers = breaks.clone();
        uppers.push(hi);

        // Assign classes.
        let class_vals: Vec<FieldValue> = raw
            .iter()
            .map(|v| match v {
                Some(x) => FieldValue::Integer(assign_class(*x, &breaks)),
                None => FieldValue::Null,
            })
            .collect();

        layer.add_field(FieldDef::new(prm.class_field.clone(), FieldType::Integer));
        let break_field_added = prm.break_field.is_some();
        if let Some(bf) = &prm.break_field {
            layer.add_field(FieldDef::new(bf.clone(), FieldType::Float));
        }

        for (fi, feat) in layer.features.iter_mut().enumerate() {
            feat.attributes.push(class_vals[fi].clone());
            if break_field_added {
                let bv = match &class_vals[fi] {
                    FieldValue::Integer(c) => {
                        FieldValue::Float(uppers[(*c as usize - 1).min(uppers.len() - 1)])
                    }
                    _ => FieldValue::Null,
                };
                feat.attributes.push(bv);
            }
        }

        // Per-class feature counts (1-based).
        let mut counts = vec![0i64; n_classes];
        for v in &class_vals {
            if let FieldValue::Integer(c) = v {
                let ci = (*c as usize - 1).min(n_classes - 1);
                counts[ci] += 1;
            }
        }

        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("valid_count".to_string(), json!(n_valid));
        outputs.insert("null_count".to_string(), json!(n - n_valid));
        outputs.insert("method".to_string(), json!(prm.method.as_str()));
        outputs.insert("classes".to_string(), json!(n_classes));
        outputs.insert("breaks".to_string(), json!(breaks));
        outputs.insert("class_field".to_string(), json!(prm.class_field));
        outputs.insert("class_counts".to_string(), json!(counts));
        if let Some(bf) = &prm.break_field {
            outputs.insert("break_field".to_string(), json!(bf));
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── Break computation ────────────────────────────────────────────────────────

/// Returns the internal (upper) class breaks for the chosen method. The number
/// of returned breaks is `classes - 1` for fixed-k methods, or data-dependent
/// for `defined_interval` and `std_dev`.
fn compute_breaks(vals: &[f64], prm: &Params) -> Result<Vec<f64>, ToolError> {
    let (lo, hi) = min_max(vals);
    let k = prm.classes.max(1);

    let breaks = match prm.method {
        Method::EqualInterval => {
            let step = (hi - lo) / k as f64;
            if step > 0.0 {
                (1..k).map(|i| lo + step * i as f64).collect()
            } else {
                Vec::new()
            }
        }
        Method::Quantile => {
            let mut sorted = vals.to_vec();
            sorted.sort_by(|a, b| a.total_cmp(b));
            (1..k)
                .map(|i| percentile(&sorted, 100.0 * i as f64 / k as f64))
                .collect()
        }
        Method::NaturalBreaks => {
            let mut sorted = vals.to_vec();
            sorted.sort_by(|a, b| a.total_cmp(b));
            jenks_breaks(&sorted, k)
        }
        Method::DefinedInterval => {
            let d = prm.interval.ok_or_else(|| {
                ToolError::Validation(
                    "method 'defined_interval' requires the 'interval' parameter".to_string(),
                )
            })?;
            if d <= 0.0 {
                return Err(ToolError::Validation(
                    "'interval' must be a positive number".to_string(),
                ));
            }
            let mut breaks = Vec::new();
            let mut b = lo + d;
            // Cap the class count so a tiny interval on a huge range cannot
            // allocate an unbounded vector.
            while b < hi && breaks.len() < 10_000 {
                breaks.push(b);
                b += d;
            }
            breaks
        }
        Method::StdDev => {
            let (mean, sd) = mean_sd(vals);
            let step = sd * prm.std_dev_interval;
            if step <= 0.0 {
                Vec::new()
            } else {
                let mut breaks = Vec::new();
                // Breaks below the mean (mean - m*step), then the mean, then above.
                let mut m = 1;
                while mean - step * (m as f64) > lo && m < 10_000 {
                    breaks.push(mean - step * m as f64);
                    m += 1;
                }
                if mean > lo && mean < hi {
                    breaks.push(mean);
                }
                let mut m = 1;
                while mean + step * (m as f64) < hi && m < 10_000 {
                    breaks.push(mean + step * m as f64);
                    m += 1;
                }
                breaks.sort_by(|a, b| a.total_cmp(b));
                breaks
            }
        }
        Method::GeometricInterval => {
            // Geometric progression of class widths. Shift the data so its
            // minimum is strictly positive, take a constant multiplier
            // r = (hi'/lo')^(1/k), then shift the breaks back.
            let shift = if lo <= 0.0 { 1.0 - lo } else { 0.0 };
            let lo_p = lo + shift;
            let hi_p = hi + shift;
            if hi_p <= lo_p || lo_p <= 0.0 {
                Vec::new()
            } else {
                let r = (hi_p / lo_p).powf(1.0 / k as f64);
                (1..k).map(|i| lo_p * r.powi(i as i32) - shift).collect()
            }
        }
    };
    Ok(breaks)
}

/// 1-based class index for `x`: class `c` is the smallest `c` such that
/// `x <= breaks[c-1]`; values above the last break fall in the top class. A
/// value exactly equal to a break stays in the lower class.
fn assign_class(x: f64, breaks: &[f64]) -> i64 {
    let mut c = 1i64;
    for &b in breaks {
        if x > b {
            c += 1;
        } else {
            break;
        }
    }
    c
}

/// Fisher–Jenks natural-breaks classification (goodness-of-variance-fit DP).
/// `sorted` must be ascending; returns the `k-1` internal upper breaks.
fn jenks_breaks(sorted: &[f64], n_classes: usize) -> Vec<f64> {
    let n = sorted.len();
    let k = n_classes.min(n);
    if k <= 1 || n < 2 {
        return Vec::new();
    }

    // lower_class_limits[i][j] and variance_combinations[i][j], 1-based over i.
    let mut lcl = vec![vec![0usize; k + 1]; n + 1];
    let mut vc = vec![vec![0f64; k + 1]; n + 1];
    for i in 1..=k {
        lcl[1][i] = 1;
        vc[1][i] = 0.0;
        for row in vc.iter_mut().take(n + 1).skip(2) {
            row[i] = f64::INFINITY;
        }
    }

    for l in 2..=n {
        let mut sum = 0.0;
        let mut sum_sq = 0.0;
        let mut w = 0.0;
        let mut variance = 0.0;
        for m in 1..=l {
            let lower = l - m + 1;
            let val = sorted[lower - 1];
            w += 1.0;
            sum += val;
            sum_sq += val * val;
            variance = sum_sq - (sum * sum) / w;
            let i4 = lower - 1;
            if i4 != 0 {
                for j in 2..=k {
                    let candidate = variance + vc[i4][j - 1];
                    if vc[l][j] >= candidate {
                        lcl[l][j] = lower;
                        vc[l][j] = candidate;
                    }
                }
            }
        }
        lcl[l][1] = 1;
        vc[l][1] = variance;
    }

    // Back-track the class limits.
    let mut kclass = vec![0f64; k + 1];
    kclass[k] = sorted[n - 1];
    kclass[0] = sorted[0];
    let mut kk = n;
    let mut count = k;
    while count > 1 {
        let idx = lcl[kk][count];
        kclass[count - 1] = sorted[idx - 2];
        kk = idx - 1;
        count -= 1;
    }
    kclass[1..k].to_vec()
}

// ── Statistics helpers ───────────────────────────────────────────────────────

fn min_max(vals: &[f64]) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &v in vals {
        if v.is_finite() {
            lo = lo.min(v);
            hi = hi.max(v);
        }
    }
    if !lo.is_finite() {
        (0.0, 0.0)
    } else {
        (lo, hi)
    }
}

fn mean_sd(vals: &[f64]) -> (f64, f64) {
    let n = vals.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    let mean = vals.iter().sum::<f64>() / n as f64;
    let var = if n > 1 {
        vals.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0)
    } else {
        0.0
    };
    (mean, var.sqrt())
}

/// Linear-interpolation percentile (numpy `'linear'` convention) over an
/// already-sorted slice.
fn percentile(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let n = sorted.len();
    if n == 1 {
        return sorted[0];
    }
    let rank = p / 100.0 * (n as f64 - 1.0);
    let lo = rank.floor() as usize;
    let hi = rank.ceil() as usize;
    if lo == hi {
        sorted[lo]
    } else {
        let frac = rank - lo as f64;
        sorted[lo] + (sorted[hi] - sorted[lo]) * frac
    }
}

// ── Parameters ───────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    EqualInterval,
    DefinedInterval,
    Quantile,
    NaturalBreaks,
    StdDev,
    GeometricInterval,
}

impl Method {
    fn as_str(self) -> &'static str {
        match self {
            Method::EqualInterval => "equal_interval",
            Method::DefinedInterval => "defined_interval",
            Method::Quantile => "quantile",
            Method::NaturalBreaks => "natural_breaks",
            Method::StdDev => "std_dev",
            Method::GeometricInterval => "geometric_interval",
        }
    }
}

struct Params {
    field: String,
    method: Method,
    classes: usize,
    interval: Option<f64>,
    std_dev_interval: f64,
    class_field: String,
    break_field: Option<String>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let field = require_str(args, "field")?.to_string();

    let method = match args.get("method").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("natural_breaks") => Method::NaturalBreaks,
        Some("equal_interval") => Method::EqualInterval,
        Some("defined_interval") => Method::DefinedInterval,
        Some("quantile") => Method::Quantile,
        Some("std_dev") => Method::StdDev,
        Some("geometric_interval") => Method::GeometricInterval,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'method' must be equal_interval/defined_interval/quantile/natural_breaks/std_dev/geometric_interval, got '{o}'"
            )))
        }
    };

    let classes = parse_optional_u64(args, "classes")?.unwrap_or(5);
    if classes < 2 {
        return Err(ToolError::Validation(
            "'classes' must be at least 2".to_string(),
        ));
    }

    let interval = parse_optional_f64(args, "interval")?;
    if method == Method::DefinedInterval {
        match interval {
            Some(d) if d > 0.0 => {}
            _ => {
                return Err(ToolError::Validation(
                    "method 'defined_interval' requires a positive 'interval'".to_string(),
                ))
            }
        }
    }

    let std_dev_interval = parse_optional_f64(args, "std_dev_interval")?.unwrap_or(1.0);
    if std_dev_interval <= 0.0 {
        return Err(ToolError::Validation(
            "'std_dev_interval' must be a positive number".to_string(),
        ));
    }

    let class_field = parse_optional_str(args, "class_field")?
        .map(str::to_string)
        .unwrap_or_else(|| "class".to_string());
    let break_field = parse_optional_str(args, "break_field")?.map(str::to_string);

    Ok(Params {
        field,
        method,
        classes: classes as usize,
        interval,
        std_dev_interval,
        class_field,
        break_field,
    })
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

fn parse_optional_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64().or_else(|| n.as_f64().map(|f| f as u64))),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be an integer"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be an integer"
        ))),
    }
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
    use wbvector::{memory_store, Geometry, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn numeric_layer(vals: &[Option<f64>]) -> String {
        let mut l = Layer::new("t")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (i, v) in vals.iter().enumerate() {
            let fv: FieldValue = match v {
                Some(x) => (*x).into(),
                None => FieldValue::Null,
            };
            l.add_feature(Some(Geometry::point(i as f64, 0.0)), &[("v", fv)])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (Layer, ToolRunResult) {
        let a: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ReclassifyFieldTool.run(&a, &ctx()).unwrap();
        let l = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (l, out)
    }

    fn classes(l: &Layer, field: &str) -> Vec<Option<i64>> {
        let i = l.schema.field_index(field).unwrap();
        l.features
            .iter()
            .map(|f| f.attributes[i].as_i64())
            .collect()
    }

    /// Class indices are monotone non-decreasing in the input value: a strictly
    /// larger value never lands in a strictly lower class.
    #[test]
    fn classes_are_monotone_in_value() {
        let vals: Vec<Option<f64>> = (0..40).map(|i| Some((i * i) as f64)).collect();
        for method in [
            "equal_interval",
            "quantile",
            "natural_breaks",
            "geometric_interval",
        ] {
            let (l, _) = run(json!({
                "input": numeric_layer(&vals),
                "field": "v",
                "method": method,
                "classes": 5,
            }));
            let c: Vec<i64> = classes(&l, "class")
                .into_iter()
                .map(|x| x.unwrap())
                .collect();
            for w in c.windows(2) {
                assert!(
                    w[1] >= w[0],
                    "method {method}: class dropped {} -> {}",
                    w[0],
                    w[1]
                );
            }
            let max = *c.iter().max().unwrap();
            assert!(
                (2..=5).contains(&max),
                "method {method}: {max} classes seen"
            );
        }
    }

    /// Quantile classification produces near-equal class counts.
    #[test]
    fn quantile_gives_near_equal_counts() {
        let vals: Vec<Option<f64>> = (0..100).map(|i| Some(i as f64)).collect();
        let (l, out) = run(json!({
            "input": numeric_layer(&vals),
            "field": "v",
            "method": "quantile",
            "classes": 5,
        }));
        assert_eq!(out.outputs["classes"], json!(5));
        let mut counts = [0i64; 6];
        for c in classes(&l, "class").into_iter().flatten() {
            counts[c as usize] += 1;
        }
        for (c, &count) in counts.iter().enumerate().skip(1) {
            assert!(
                (count - 20).abs() <= 1,
                "class {c}: {count} feats, expected ~20"
            );
        }
    }

    /// Equal-interval breaks are evenly spaced between min and max.
    #[test]
    fn equal_interval_breaks_are_even() {
        let vals: Vec<Option<f64>> = (0..=100).map(|i| Some(i as f64)).collect();
        let (_l, out) = run(json!({
            "input": numeric_layer(&vals),
            "field": "v",
            "method": "equal_interval",
            "classes": 5,
        }));
        let breaks: Vec<f64> = serde_json::from_value(out.outputs["breaks"].clone()).unwrap();
        let expected = [20.0, 40.0, 60.0, 80.0];
        assert_eq!(breaks.len(), 4);
        for (b, e) in breaks.iter().zip(expected) {
            assert!((b - e).abs() < 1e-9, "break {b} != {e}");
        }
    }

    /// Natural breaks separates two tight, well-separated clusters cleanly.
    #[test]
    fn natural_breaks_splits_clusters() {
        let mut vals: Vec<Option<f64>> = Vec::new();
        for _ in 0..10 {
            vals.push(Some(1.0));
        }
        for _ in 0..10 {
            vals.push(Some(100.0));
        }
        let (l, out) = run(json!({
            "input": numeric_layer(&vals),
            "field": "v",
            "method": "natural_breaks",
            "classes": 2,
        }));
        let breaks: Vec<f64> = serde_json::from_value(out.outputs["breaks"].clone()).unwrap();
        assert_eq!(breaks.len(), 1);
        assert!(
            breaks[0] >= 1.0 && breaks[0] < 100.0,
            "break at {}",
            breaks[0]
        );
        let c = classes(&l, "class");
        let low = c.iter().take(10).all(|x| *x == Some(1));
        let high = c.iter().skip(10).all(|x| *x == Some(2));
        assert!(low && high, "clusters not split: {c:?}");
    }

    /// The optional break field records the upper class limit; null values stay null.
    #[test]
    fn break_field_and_null_handling() {
        let vals = [Some(1.0), None, Some(50.0), Some(100.0)];
        let (l, out) = run(json!({
            "input": numeric_layer(&vals),
            "field": "v",
            "method": "equal_interval",
            "classes": 4,
            "break_field": "cls_max",
        }));
        assert_eq!(out.outputs["null_count"], json!(1));
        assert_eq!(out.outputs["valid_count"], json!(3));
        let ci = l.schema.field_index("class").unwrap();
        let bi = l.schema.field_index("cls_max").unwrap();
        // the null-valued feature has null class and null break
        assert!(l.features[1].attributes[ci].as_i64().is_none());
        assert!(l.features[1].attributes[bi].as_f64().is_none());
        // finite features have a finite break value >= their own value
        for (fi, v) in vals.iter().enumerate() {
            if let Some(x) = v {
                let b = l.features[fi].attributes[bi].as_f64().unwrap();
                assert!(b + 1e-9 >= *x, "break {b} < value {x}");
            }
        }
    }

    /// defined_interval yields fixed-width classes and a data-driven count.
    #[test]
    fn defined_interval_fixed_width() {
        let vals: Vec<Option<f64>> = (0..=100).map(|i| Some(i as f64)).collect();
        let (_l, out) = run(json!({
            "input": numeric_layer(&vals),
            "field": "v",
            "method": "defined_interval",
            "interval": 25.0,
        }));
        let breaks: Vec<f64> = serde_json::from_value(out.outputs["breaks"].clone()).unwrap();
        assert_eq!(breaks, vec![25.0, 50.0, 75.0]);
        assert_eq!(out.outputs["classes"], json!(4));
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let a: ToolArgs = serde_json::from_value(v).unwrap();
            ReclassifyFieldTool.validate(&a)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err()); // no field
        assert!(bad(json!({ "input": "a.geojson", "field": "v", "method": "kmeans" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "field": "v", "classes": 1 })).is_err());
        assert!(
            bad(json!({ "input": "a.geojson", "field": "v", "method": "defined_interval" }))
                .is_err()
        );
        assert!(bad(json!({
            "input": "a.geojson", "field": "v", "method": "defined_interval", "interval": -5
        }))
        .is_err());
        assert!(bad(json!({ "input": "a.geojson", "field": "v" })).is_ok());
    }

    #[test]
    fn rejects_unknown_field() {
        let vals = [Some(1.0), Some(2.0), Some(3.0)];
        let a: ToolArgs = serde_json::from_value(json!({
            "input": numeric_layer(&vals),
            "field": "nope",
            "method": "quantile",
        }))
        .unwrap();
        assert!(ReclassifyFieldTool.run(&a, &ctx()).is_err());
    }
}
