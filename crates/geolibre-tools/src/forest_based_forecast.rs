//! GeoLibre tool: non-parametric per-location time-series forecasting.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Forest-based Forecast* (Space Time
//! Pattern Mining). The shipped `time_series_forecast` covers `linear`,
//! `parabolic`, `exp_smoothing` and `arima` — all parametric models assuming a
//! specific functional form. Series with threshold effects, regime changes or
//! non-linear seasonality are poorly served by every one of them, and neither
//! registry offers a non-parametric option.
//!
//! The model is a random forest over **lag windows**: predictors are the
//! previous `time_window` values, the response is the next value. Forecasting
//! is recursive — predict one step, append it to the lag vector, repeat.
//!
//! Two properties matter for correctness here:
//!
//! * **Determinism.** The WASM path has no ambient RNG and reproducibility is
//!   required, so bootstrap sampling and feature subsetting are driven by an
//!   explicit `seed` through a small splitmix64 generator rather than any
//!   platform source.
//! * **Honest validation.** With `validation_steps > 0` the forest is refit on
//!   the truncated series before scoring the held-out tail, so the reported
//!   RMSE never sees the data it is judged on. The column is written in the
//!   same shape `time_series_forecast` uses.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Deterministic splitmix64: no platform RNG, identical results everywhere.
struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: usize) -> usize {
        if n == 0 {
            0
        } else {
            (self.next_u64() % n as u64) as usize
        }
    }
}

/// A binary regression tree over lag features.
enum Node {
    Leaf(f64),
    Split {
        feature: usize,
        threshold: f64,
        left: Box<Node>,
        right: Box<Node>,
    },
}

impl Node {
    fn predict(&self, x: &[f64]) -> f64 {
        match self {
            Node::Leaf(v) => *v,
            Node::Split {
                feature,
                threshold,
                left,
                right,
            } => {
                if x[*feature] <= *threshold {
                    left.predict(x)
                } else {
                    right.predict(x)
                }
            }
        }
    }
}

pub struct ForestBasedForecastTool;

impl Tool for ForestBasedForecastTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "forest_based_forecast",
            display_name: "Forest-based Forecast",
            summary: "Forecast a value forward at each location with a random forest trained on sliding lag windows of that location's own history, like ArcGIS Forest-based Forecast.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Timestamped features carrying a location id, a time field and a value field.",
                    required: true,
                },
                ToolParamSpec {
                    name: "location_field",
                    description: "Field identifying the location each observation belongs to.",
                    required: true,
                },
                ToolParamSpec {
                    name: "time_field",
                    description: "Numeric time field used to order each location's series.",
                    required: true,
                },
                ToolParamSpec {
                    name: "value_field",
                    description: "Numeric field to forecast.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output forecast table path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "forecast_steps",
                    description: "Number of steps to forecast forward (default 3).",
                    required: false,
                },
                ToolParamSpec {
                    name: "time_window",
                    description: "Number of prior values used as predictors (default: a quarter of the series length, at least 2).",
                    required: false,
                },
                ToolParamSpec {
                    name: "validation_steps",
                    description: "Trailing steps held out to score RMSE (default 0, no validation).",
                    required: false,
                },
                ToolParamSpec {
                    name: "n_trees",
                    description: "Number of trees in the forest (default 50, max 500).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_leaf_size",
                    description: "Minimum training rows in a leaf (default 2).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_depth",
                    description: "Maximum tree depth (default 8).",
                    required: false,
                },
                ToolParamSpec {
                    name: "seed",
                    description: "Seed for the deterministic RNG driving bootstrapping and feature subsetting (default 42).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for k in ["input", "location_field", "time_field", "value_field"] {
            require_str(args, k)?;
        }
        for (k, lo, hi) in [
            ("forecast_steps", 1.0, 1000.0),
            ("time_window", 1.0, 1000.0),
            ("validation_steps", 0.0, 1000.0),
            ("n_trees", 1.0, 500.0),
            ("min_leaf_size", 1.0, 1000.0),
            ("max_depth", 1.0, 64.0),
        ] {
            if let Some(v) = parse_optional_f64(args, k)? {
                if !v.is_finite() || v < lo || v > hi {
                    return Err(ToolError::Validation(format!(
                        "'{k}' must be between {lo} and {hi}"
                    )));
                }
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let loc_field = require_str(args, "location_field")?;
        let time_field = require_str(args, "time_field")?;
        let value_field = require_str(args, "value_field")?;
        let output = parse_optional_str(args, "output")?;
        let steps = parse_optional_f64(args, "forecast_steps")?.unwrap_or(3.0) as usize;
        let window_opt = parse_optional_f64(args, "time_window")?.map(|v| v as usize);
        let validation = parse_optional_f64(args, "validation_steps")?.unwrap_or(0.0) as usize;
        let n_trees = parse_optional_f64(args, "n_trees")?.unwrap_or(50.0) as usize;
        let min_leaf = parse_optional_f64(args, "min_leaf_size")?.unwrap_or(2.0) as usize;
        let max_depth = parse_optional_f64(args, "max_depth")?.unwrap_or(8.0) as usize;
        let seed = parse_optional_f64(args, "seed")?.unwrap_or(42.0) as u64;

        let layer = load_input_layer(input)?;
        for f in [loc_field, time_field, value_field] {
            if layer.schema.field_index(f).is_none() {
                return Err(ToolError::Validation(format!(
                    "field '{f}' not found on the input layer"
                )));
            }
        }

        // Group observations into per-location series.
        let mut series: BTreeMap<String, Vec<(f64, f64)>> = BTreeMap::new();
        for feat in layer.features.iter() {
            let (Ok(lv), Ok(tv), Ok(vv)) = (
                feat.get(&layer.schema, loc_field),
                feat.get(&layer.schema, time_field),
                feat.get(&layer.schema, value_field),
            ) else {
                continue;
            };
            let (Some(t), Some(v)) = (tv.as_f64(), vv.as_f64()) else {
                continue;
            };
            if !t.is_finite() || !v.is_finite() {
                continue;
            }
            series.entry(field_string(lv)).or_default().push((t, v));
        }

        ctx.progress
            .info(&format!("forecasting {} location series", series.len()));

        let mut out = Layer::new("forest_forecast");
        out.add_field(FieldDef::new("location", FieldType::Text));
        out.add_field(FieldDef::new("step", FieldType::Integer));
        out.add_field(FieldDef::new("time", FieldType::Float));
        out.add_field(FieldDef::new("forecast", FieldType::Float));
        out.add_field(FieldDef::new("rmse", FieldType::Float));
        out.add_field(FieldDef::new("model", FieldType::Text));

        let mut forecast_locations = 0usize;
        let mut skipped = 0usize;

        let location_total = series.len();
        for (li, (name, obs)) in series.iter_mut().enumerate() {
            obs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            let values: Vec<f64> = obs.iter().map(|(_, v)| *v).collect();

            // Default window is a quarter of the series, bounded so there is
            // always at least one training row left over.
            let window = window_opt
                .unwrap_or((values.len() / 4).max(2))
                .clamp(1, values.len().saturating_sub(1).max(1));

            if values.len() < window + 2 {
                skipped += 1;
                continue;
            }

            // Validation: refit on the truncated series so the score never sees
            // the data it is judged on.
            let rmse = if validation > 0 && values.len() > window + validation + 1 {
                let cut = values.len() - validation;
                let train = &values[..cut];
                let truth = &values[cut..];
                let f = Forest::fit(train, window, n_trees, min_leaf, max_depth, seed);
                let pred = f.forecast(train, window, validation);
                let se: f64 = pred
                    .iter()
                    .zip(truth.iter())
                    .map(|(p, t)| (p - t).powi(2))
                    .sum();
                Some((se / truth.len() as f64).sqrt())
            } else {
                None
            };

            let forest = Forest::fit(&values, window, n_trees, min_leaf, max_depth, seed);
            let preds = forest.forecast(&values, window, steps);

            // Extrapolate the time axis using the median observed step, which is
            // robust to an irregular series in a way the mean is not.
            let step_size = median_step(obs);
            let last_t = obs.last().map(|(t, _)| *t).unwrap_or(0.0);

            for (k, p) in preds.iter().enumerate() {
                out.add_feature(
                    None,
                    &[
                        ("location", FieldValue::Text(name.clone())),
                        ("step", FieldValue::Integer(k as i64 + 1)),
                        (
                            "time",
                            FieldValue::Float(last_t + step_size * (k as f64 + 1.0)),
                        ),
                        ("forecast", FieldValue::Float(*p)),
                        (
                            "rmse",
                            rmse.map(FieldValue::Float).unwrap_or(FieldValue::Null),
                        ),
                        ("model", FieldValue::Text("forest".into())),
                    ],
                )
                .map_err(|e| ToolError::Execution(format!("failed writing forecast: {e}")))?;
            }
            forecast_locations += 1;
            ctx.progress
                .progress((li as f64 + 1.0) / location_total.max(1) as f64);
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("location_count".to_string(), json!(series.len()));
        outputs.insert("forecast_locations".to_string(), json!(forecast_locations));
        // Surfaced so a short series silently producing no forecast is visible.
        outputs.insert("skipped_locations".to_string(), json!(skipped));
        Ok(ToolRunResult { outputs })
    }
}

struct Forest {
    trees: Vec<Node>,
}

impl Forest {
    fn fit(
        values: &[f64],
        window: usize,
        n_trees: usize,
        min_leaf: usize,
        max_depth: usize,
        seed: u64,
    ) -> Forest {
        // Build the lag matrix: X[i] = values[i..i+window], y[i] = values[i+window].
        let mut xs: Vec<Vec<f64>> = Vec::new();
        let mut ys: Vec<f64> = Vec::new();
        for i in 0..values.len().saturating_sub(window) {
            xs.push(values[i..i + window].to_vec());
            ys.push(values[i + window]);
        }
        let mut rng = Rng(seed);
        let mut trees = Vec::with_capacity(n_trees);
        for _ in 0..n_trees.max(1) {
            // Bootstrap sample.
            let idx: Vec<usize> = (0..xs.len()).map(|_| rng.below(xs.len().max(1))).collect();
            trees.push(build_tree(&xs, &ys, &idx, min_leaf, max_depth, &mut rng));
        }
        Forest { trees }
    }

    fn predict(&self, x: &[f64]) -> f64 {
        if self.trees.is_empty() {
            return 0.0;
        }
        self.trees.iter().map(|t| t.predict(x)).sum::<f64>() / self.trees.len() as f64
    }

    /// Recursive multi-step forecast: each prediction is appended to the lag
    /// vector and feeds the next step.
    fn forecast(&self, history: &[f64], window: usize, steps: usize) -> Vec<f64> {
        let mut lags: Vec<f64> = history[history.len().saturating_sub(window)..].to_vec();
        while lags.len() < window {
            lags.insert(0, *history.first().unwrap_or(&0.0));
        }
        let mut out = Vec::with_capacity(steps);
        for _ in 0..steps {
            let p = self.predict(&lags);
            out.push(p);
            lags.remove(0);
            lags.push(p);
        }
        out
    }
}

fn build_tree(
    xs: &[Vec<f64>],
    ys: &[f64],
    idx: &[usize],
    min_leaf: usize,
    depth: usize,
    rng: &mut Rng,
) -> Node {
    let mean = if idx.is_empty() {
        0.0
    } else {
        idx.iter().map(|&i| ys[i]).sum::<f64>() / idx.len() as f64
    };
    if depth == 0 || idx.len() <= min_leaf || xs.is_empty() {
        return Node::Leaf(mean);
    }

    let n_features = xs[0].len();
    // Random feature subset, the "random" half of a random forest. sqrt(p) is
    // the standard choice; at least 1.
    let k = ((n_features as f64).sqrt().ceil() as usize)
        .max(1)
        .min(n_features);
    let mut candidates: Vec<usize> = Vec::with_capacity(k);
    for _ in 0..k {
        candidates.push(rng.below(n_features));
    }

    let mut best: Option<(f64, usize, f64)> = None; // (sse, feature, threshold)
    for &f in &candidates {
        let mut vals: Vec<f64> = idx.iter().map(|&i| xs[i][f]).collect();
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        vals.dedup();
        if vals.len() < 2 {
            continue;
        }
        for w in vals.windows(2) {
            let thr = (w[0] + w[1]) / 2.0;
            let (mut ls, mut lc, mut rs, mut rc) = (0.0, 0usize, 0.0, 0usize);
            for &i in idx {
                if xs[i][f] <= thr {
                    ls += ys[i];
                    lc += 1;
                } else {
                    rs += ys[i];
                    rc += 1;
                }
            }
            if lc == 0 || rc == 0 {
                continue;
            }
            let (lm, rm) = (ls / lc as f64, rs / rc as f64);
            let sse: f64 = idx
                .iter()
                .map(|&i| {
                    let m = if xs[i][f] <= thr { lm } else { rm };
                    (ys[i] - m).powi(2)
                })
                .sum();
            if best.is_none_or(|(bs, _, _)| sse < bs) {
                best = Some((sse, f, thr));
            }
        }
    }

    match best {
        Some((_, f, thr)) => {
            let left: Vec<usize> = idx.iter().copied().filter(|&i| xs[i][f] <= thr).collect();
            let right: Vec<usize> = idx.iter().copied().filter(|&i| xs[i][f] > thr).collect();
            if left.is_empty() || right.is_empty() {
                return Node::Leaf(mean);
            }
            Node::Split {
                feature: f,
                threshold: thr,
                left: Box::new(build_tree(xs, ys, &left, min_leaf, depth - 1, rng)),
                right: Box::new(build_tree(xs, ys, &right, min_leaf, depth - 1, rng)),
            }
        }
        None => Node::Leaf(mean),
    }
}

/// Median gap between consecutive observations. Robust to an irregular series
/// where the mean would be dragged by one long gap.
fn median_step(obs: &[(f64, f64)]) -> f64 {
    if obs.len() < 2 {
        return 1.0;
    }
    let mut gaps: Vec<f64> = obs.windows(2).map(|w| w[1].0 - w[0].0).collect();
    gaps.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let m = gaps[gaps.len() / 2];
    if m > 0.0 {
        m
    } else {
        1.0
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Coord, Geometry, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Observations as (location, time, value).
    fn series_layer(items: Vec<(&str, f64, f64)>) -> String {
        let mut l = Layer::new("s")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("loc", FieldType::Text));
        l.add_field(FieldDef::new("t", FieldType::Float));
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (loc, t, v) in items {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(0.0, 0.0))),
                &[
                    ("loc", FieldValue::Text(loc.to_string())),
                    ("t", FieldValue::Float(t)),
                    ("v", FieldValue::Float(v)),
                ],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let mut v = args;
        v["location_field"] = json!("loc");
        v["time_field"] = json!("t");
        v["value_field"] = json!("v");
        let args: ToolArgs = serde_json::from_value(v).unwrap();
        let out = ForestBasedForecastTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn forecasts(layer: &Layer, loc: &str) -> Vec<f64> {
        let li = layer.schema.field_index("location").unwrap();
        let fi = layer.schema.field_index("forecast").unwrap();
        layer
            .features
            .iter()
            .filter(|f| matches!(&f.attributes[li], FieldValue::Text(s) if s == loc))
            .map(|f| f.attributes[fi].as_f64().unwrap())
            .collect()
    }

    /// A repeating square wave is exactly the non-linear pattern a parametric
    /// linear/parabolic fit cannot represent but a lag forest can.
    #[test]
    fn learns_a_repeating_pattern() {
        let mut items = Vec::new();
        for i in 0..40 {
            let v = if (i / 2) % 2 == 0 { 10.0 } else { 90.0 };
            items.push(("a", i as f64, v));
        }
        let (_, layer) = run(json!({
            "input": series_layer(items), "forecast_steps": 4,
            "time_window": 4, "n_trees": 30
        }));
        let preds = forecasts(&layer, "a");
        assert_eq!(preds.len(), 4);
        // Predictions must stay inside the observed range, not run off.
        for p in &preds {
            assert!(
                *p >= 5.0 && *p <= 95.0,
                "forecast {p} escaped the observed range"
            );
        }
    }

    /// A flat series forecasts flat — the basic sanity property.
    #[test]
    fn constant_series_forecasts_the_constant() {
        let items: Vec<(&str, f64, f64)> = (0..20).map(|i| ("a", i as f64, 42.0)).collect();
        let (_, layer) = run(json!({
            "input": series_layer(items), "forecast_steps": 3, "time_window": 3
        }));
        for p in forecasts(&layer, "a") {
            assert!((p - 42.0).abs() < 1e-9, "expected 42, got {p}");
        }
    }

    /// THE reproducibility requirement: the same seed must give identical
    /// output, since the WASM path has no ambient RNG.
    #[test]
    fn seed_makes_output_deterministic() {
        let items: Vec<(&str, f64, f64)> = (0..30)
            .map(|i| ("a", i as f64, ((i * 7) % 13) as f64))
            .collect();
        let input = series_layer(items);
        let (_, a) = run(json!({
            "input": input, "forecast_steps": 3, "time_window": 3, "seed": 7
        }));
        let (_, b) = run(json!({
            "input": input, "forecast_steps": 3, "time_window": 3, "seed": 7
        }));
        assert_eq!(forecasts(&a, "a"), forecasts(&b, "a"));

        let (_, c) = run(json!({
            "input": input, "forecast_steps": 3, "time_window": 3, "seed": 99
        }));
        // A different seed is allowed to differ; the point is that 7 == 7.
        assert_eq!(forecasts(&a, "a").len(), forecasts(&c, "a").len());
    }

    /// Locations are forecast independently.
    #[test]
    fn locations_are_independent() {
        let mut items = Vec::new();
        for i in 0..20 {
            items.push(("low", i as f64, 1.0));
            items.push(("high", i as f64, 100.0));
        }
        let (out, layer) = run(json!({
            "input": series_layer(items), "forecast_steps": 2, "time_window": 3
        }));
        assert_eq!(out.outputs["forecast_locations"], json!(2));
        for p in forecasts(&layer, "low") {
            assert!((p - 1.0).abs() < 1e-9);
        }
        for p in forecasts(&layer, "high") {
            assert!((p - 100.0).abs() < 1e-9);
        }
    }

    /// validation_steps produces an RMSE scored on held-out data.
    #[test]
    fn validation_reports_rmse() {
        let items: Vec<(&str, f64, f64)> = (0..30).map(|i| ("a", i as f64, 5.0)).collect();
        let (_, layer) = run(json!({
            "input": series_layer(items), "forecast_steps": 2,
            "time_window": 3, "validation_steps": 5
        }));
        let ri = layer.schema.field_index("rmse").unwrap();
        let rmse = layer.features[0].attributes[ri].as_f64().unwrap();
        // A constant series is perfectly predictable, so held-out error is ~0.
        assert!(rmse < 1e-6, "expected near-zero RMSE, got {rmse}");
    }

    /// Time is extrapolated by the median observed step.
    #[test]
    fn forecast_times_extend_the_series() {
        let items: Vec<(&str, f64, f64)> = (0..20).map(|i| ("a", i as f64 * 10.0, 3.0)).collect();
        let (_, layer) = run(json!({
            "input": series_layer(items), "forecast_steps": 2, "time_window": 3
        }));
        let ti = layer.schema.field_index("time").unwrap();
        let times: Vec<f64> = layer
            .features
            .iter()
            .map(|f| f.attributes[ti].as_f64().unwrap())
            .collect();
        // Last observation is at t=190, step is 10.
        assert!((times[0] - 200.0).abs() < 1e-9, "got {}", times[0]);
        assert!((times[1] - 210.0).abs() < 1e-9, "got {}", times[1]);
    }

    /// A series too short for the window is skipped, and that is reported
    /// rather than silently producing nothing.
    #[test]
    fn short_series_is_skipped_and_counted() {
        let items = vec![("a", 0.0, 1.0), ("a", 1.0, 2.0)];
        let (out, _) = run(json!({
            "input": series_layer(items), "forecast_steps": 2, "time_window": 5
        }));
        assert_eq!(out.outputs["skipped_locations"], json!(1));
        assert_eq!(out.outputs["forecast_locations"], json!(0));
    }

    #[test]
    fn rejects_bad_parameters() {
        let p = series_layer(vec![("a", 0.0, 1.0)]);
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ForestBasedForecastTool.validate(&args).is_err()
        };
        assert!(bad(
            json!({ "location_field": "loc", "time_field": "t", "value_field": "v" })
        ));
        assert!(bad(
            json!({ "input": p, "time_field": "t", "value_field": "v" })
        ));
        assert!(bad(
            json!({ "input": p, "location_field": "loc", "time_field": "t",
                           "value_field": "v", "n_trees": 0 })
        ));
        assert!(bad(
            json!({ "input": p, "location_field": "loc", "time_field": "t",
                           "value_field": "v", "max_depth": 100 })
        ));
    }

    /// The RNG is a pure function of its seed.
    #[test]
    fn rng_is_reproducible() {
        let a: Vec<u64> = {
            let mut r = Rng(1);
            (0..5).map(|_| r.next_u64()).collect()
        };
        let b: Vec<u64> = {
            let mut r = Rng(1);
            (0..5).map(|_| r.next_u64()).collect()
        };
        assert_eq!(a, b);
        let c: Vec<u64> = {
            let mut r = Rng(2);
            (0..5).map(|_| r.next_u64()).collect()
        };
        assert_ne!(a, c);
    }
}
