//! GeoLibre tool: Kaplan-Meier survival estimate of time until an event.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Estimate Time To Event* (Spatial
//! Statistics). Each feature carries a *duration* (`age_field`) — the length of
//! time it was observed — and an *event indicator* (`event_field`) that records
//! whether the event of interest actually occurred during that time (`1`/true)
//! or whether observation ended before it did (censored, `0`/false). From these
//! the tool builds the non-parametric Kaplan-Meier survival curve
//!
//! ```text
//!   S(t) = ∏_{t_i ≤ t} ( 1 − d_i / n_i )
//! ```
//!
//! where the product runs over distinct event times `t_i`, `d_i` is the number
//! of events at `t_i`, and `n_i` the number still at risk (duration ≥ `t_i`).
//! Censored observations never lower the curve but do shrink the risk set,
//! which is exactly what distinguishes survival analysis from a naive fraction.
//!
//! An optional `stratify_field` fits a separate curve per category (e.g.
//! treatment vs. control), which is the usual comparison in survival studies.
//!
//! Output fields (original attributes preserved):
//! * `km_survival` — the estimated survival probability S(t) at the feature's
//!   own duration, read off its stratum's curve.
//! * `km_median`   — the stratum's estimated median time-to-event (smallest `t`
//!   with S(t) ≤ 0.5); `Null` when the curve never reaches 0.5 ("not reached").
//! * `km_stratum`  — the stratum label the feature was assigned to (`"all"`
//!   when no stratification field is given).
//!
//! The `outputs` map also reports overall event/censor counts and a per-stratum
//! summary (n, events, median) for programmatic use.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct EstimateTimeToEventTool;

impl Tool for EstimateTimeToEventTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "estimate_time_to_event",
            display_name: "Estimate Time To Event",
            summary: "Estimate the time until an event via a Kaplan-Meier survival curve (like ArcGIS Estimate Time To Event): from a duration field and an event/censoring field, compute each feature's survival probability and its stratum's median time-to-event, optionally stratified by a categorical field.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point or polygon layer with a duration field and an event/censoring field.",
                    required: true,
                },
                ToolParamSpec {
                    name: "age_field",
                    description: "Field holding each feature's observed duration (time-to-event or time-to-censoring), a non-negative number.",
                    required: true,
                },
                ToolParamSpec {
                    name: "event_field",
                    description: "Field indicating whether the event occurred (truthy: 1/true) or the observation was censored (falsy: 0/false).",
                    required: true,
                },
                ToolParamSpec {
                    name: "stratify_field",
                    description: "Optional categorical field to fit a separate survival curve per group.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output layer with km_survival, km_median and km_stratum fields added. If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "age_field")?;
        require_str(args, "event_field")?;
        // stratify_field / output are optional but must be strings if present.
        parse_optional_str(args, "stratify_field")?;
        parse_optional_str(args, "output")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let age_field = require_str(args, "age_field")?.to_string();
        let event_field = require_str(args, "event_field")?.to_string();
        let stratify_field = parse_optional_str(args, "stratify_field")?.map(str::to_string);
        let output = parse_optional_str(args, "output")?;

        let mut layer = load_input_layer(input)?;
        let n = layer.features.len();

        let age_idx = layer
            .schema
            .field_index(&age_field)
            .ok_or_else(|| ToolError::Validation(format!("age_field '{age_field}' not found")))?;
        let event_idx = layer.schema.field_index(&event_field).ok_or_else(|| {
            ToolError::Validation(format!("event_field '{event_field}' not found"))
        })?;
        let strat_idx =
            match &stratify_field {
                Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                    ToolError::Validation(format!("stratify_field '{f}' not found"))
                })?),
                None => None,
            };

        // Per-feature observations. `age` is None when the duration is missing or
        // negative (invalid) — such features are excluded from every curve and
        // receive Null outputs, but are still passed through.
        let mut ages: Vec<Option<f64>> = Vec::with_capacity(n);
        let mut events: Vec<bool> = Vec::with_capacity(n);
        let mut strata: Vec<String> = Vec::with_capacity(n);

        for feat in &layer.features {
            let age = feat
                .attributes
                .get(age_idx)
                .and_then(|v| v.as_f64())
                .filter(|a| a.is_finite() && *a >= 0.0);
            let event = feat
                .attributes
                .get(event_idx)
                .map(is_event)
                .unwrap_or(false);
            let stratum = match strat_idx {
                Some(idx) => stratum_key(feat.attributes.get(idx)),
                None => "all".to_string(),
            };
            ages.push(age);
            events.push(event);
            strata.push(stratum);
        }

        // Group valid observations by stratum and fit a Kaplan-Meier curve each.
        let mut grouped: BTreeMap<String, Vec<(f64, bool)>> = BTreeMap::new();
        for i in 0..n {
            if let Some(a) = ages[i] {
                grouped
                    .entry(strata[i].clone())
                    .or_default()
                    .push((a, events[i]));
            }
        }

        let mut curves: BTreeMap<String, KmCurve> = BTreeMap::new();
        for (label, obs) in &grouped {
            curves.insert(label.clone(), KmCurve::fit(obs));
        }

        ctx.progress.info(&format!(
            "{} feature(s), {} stratum/strata",
            n,
            curves.len().max(1)
        ));

        // Add output fields and populate per feature from its stratum's curve.
        layer.add_field(FieldDef::new("km_survival", FieldType::Float));
        layer.add_field(FieldDef::new("km_median", FieldType::Float));
        layer.add_field(FieldDef::new("km_stratum", FieldType::Text));

        for i in 0..n {
            let (surv, median, label) = match ages[i] {
                Some(a) => {
                    let curve = &curves[&strata[i]];
                    (
                        FieldValue::Float(curve.survival_at(a)),
                        curve
                            .median()
                            .map(FieldValue::Float)
                            .unwrap_or(FieldValue::Null),
                        FieldValue::Text(strata[i].clone()),
                    )
                }
                None => (FieldValue::Null, FieldValue::Null, FieldValue::Null),
            };
            let f = &mut layer.features[i];
            f.attributes.push(surv);
            f.attributes.push(median);
            f.attributes.push(label);
        }

        let event_count = (0..n).filter(|&i| ages[i].is_some() && events[i]).count();
        let censored_count = (0..n).filter(|&i| ages[i].is_some() && !events[i]).count();

        // Per-stratum summary for programmatic use / validation.
        let strata_summary: BTreeMap<String, Value> = curves
            .iter()
            .map(|(label, curve)| {
                let obs = &grouped[label];
                let ev = obs.iter().filter(|(_, e)| *e).count();
                (
                    label.clone(),
                    json!({
                        "n": obs.len(),
                        "events": ev,
                        "median": curve.median(),
                    }),
                )
            })
            .collect();

        let out_path = write_or_store_layer(layer, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("event_count".to_string(), json!(event_count));
        outputs.insert("censored_count".to_string(), json!(censored_count));
        outputs.insert("stratum_count".to_string(), json!(curves.len()));
        outputs.insert("strata".to_string(), json!(strata_summary));
        Ok(ToolRunResult { outputs })
    }
}

/// A fitted Kaplan-Meier survival curve: for each distinct *event* time (in
/// ascending order) the survival probability that holds from that time up to
/// (but not including) the next event time.
struct KmCurve {
    /// `(event_time, survival_after)` steps, ascending by time.
    steps: Vec<(f64, f64)>,
}

impl KmCurve {
    /// Fits the curve from `(duration, event_occurred)` observations. Censored
    /// observations (event = false) contribute to the risk set only.
    fn fit(obs: &[(f64, bool)]) -> Self {
        // Distinct times at which at least one event occurred, ascending.
        let mut event_times: Vec<f64> = obs.iter().filter(|(_, e)| *e).map(|(t, _)| *t).collect();
        event_times.sort_by(f64::total_cmp);
        event_times.dedup_by(|a, b| a == b);

        let mut steps = Vec::with_capacity(event_times.len());
        let mut survival = 1.0f64;
        for &t in &event_times {
            // Number at risk: observations with duration ≥ t.
            let n_i = obs.iter().filter(|(dur, _)| *dur >= t).count() as f64;
            // Events exactly at t.
            let d_i = obs.iter().filter(|(dur, e)| *e && *dur == t).count() as f64;
            if n_i > 0.0 {
                survival *= 1.0 - d_i / n_i;
            }
            steps.push((t, survival));
        }
        KmCurve { steps }
    }

    /// Survival probability S(t): the value of the step function at time `t`.
    /// Equals 1.0 before the first event time and is right-continuous at steps.
    fn survival_at(&self, t: f64) -> f64 {
        let mut s = 1.0;
        for &(time, surv) in &self.steps {
            if time <= t {
                s = surv;
            } else {
                break;
            }
        }
        s
    }

    /// Median survival time: the smallest event time at which S(t) ≤ 0.5.
    /// `None` when the curve never falls to 0.5 ("median not reached").
    fn median(&self) -> Option<f64> {
        self.steps
            .iter()
            .find(|(_, surv)| *surv <= 0.5)
            .map(|(t, _)| *t)
    }
}

/// Interprets a field value as an event indicator: numbers ≥ 0.5 or `true`
/// count as an observed event; everything else (0, false, null, text) as
/// censored.
fn is_event(v: &FieldValue) -> bool {
    if let Some(b) = v.as_bool() {
        return b;
    }
    v.as_f64().map(|x| x >= 0.5).unwrap_or(false)
}

/// Stable label for a stratum value (missing/null -> `"<null>"`).
fn stratum_key(v: Option<&FieldValue>) -> String {
    match v {
        None | Some(FieldValue::Null) => "<null>".to_string(),
        Some(fv) => fv.to_string(),
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

    /// Builds a point layer with duration `age`, event `evt`, group `grp`.
    fn layer(rows: &[(f64, f64, i64, &str)]) -> String {
        let mut l = Layer::new("surv")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("age", FieldType::Float));
        l.add_field(FieldDef::new("evt", FieldType::Integer));
        l.add_field(FieldDef::new("grp", FieldType::Text));
        for (i, (age, _y, evt, grp)) in rows.iter().enumerate() {
            l.add_feature(
                Some(Geometry::point(i as f64, 0.0)),
                &[
                    ("age", (*age).into()),
                    ("evt", (*evt).into()),
                    ("grp", (*grp).into()),
                ],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = EstimateTimeToEventTool.run(&args, &ctx()).unwrap();
        let l = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (l, out)
    }

    fn col_f64(l: &Layer, name: &str) -> Vec<Option<f64>> {
        let i = l.schema.field_index(name).unwrap();
        l.features
            .iter()
            .map(|f| match &f.attributes[i] {
                FieldValue::Null => None,
                v => v.as_f64(),
            })
            .collect()
    }

    /// Textbook Kaplan-Meier curve. Durations/events:
    /// (2,1),(3,1),(5,0),(6,1) over n=4.
    ///   t=2: n=4,d=1 -> S=0.75
    ///   t=3: n=3,d=1 -> S=0.50
    ///   t=6: n=1,d=1 -> S=0.00
    /// Median = 3 (first time S ≤ 0.5).
    #[test]
    fn km_curve_matches_hand_computation() {
        let rows = [
            (2.0, 0.0, 1, "a"),
            (3.0, 0.0, 1, "a"),
            (5.0, 0.0, 0, "a"),
            (6.0, 0.0, 1, "a"),
        ];
        let (l, out) = run(json!({
            "input": layer(&rows), "age_field": "age", "event_field": "evt",
        }));
        let surv = col_f64(&l, "km_survival");
        // km_survival is S at each feature's own duration.
        assert!((surv[0].unwrap() - 0.75).abs() < 1e-9); // age 2
        assert!((surv[1].unwrap() - 0.50).abs() < 1e-9); // age 3
        assert!((surv[2].unwrap() - 0.50).abs() < 1e-9); // age 5 (censored) -> S(5)=S(3)
        assert!((surv[3].unwrap() - 0.00).abs() < 1e-9); // age 6
        let median = col_f64(&l, "km_median");
        for m in &median {
            assert!((m.unwrap() - 3.0).abs() < 1e-9);
        }
        assert_eq!(out.outputs["event_count"].as_u64().unwrap(), 3);
        assert_eq!(out.outputs["censored_count"].as_u64().unwrap(), 1);
    }

    /// When the survival curve never drops to 0.5, the median is "not reached"
    /// (Null), and censored observations keep the curve above 0.5.
    #[test]
    fn median_not_reached_is_null() {
        // Only one event early; the rest censored -> S bottoms out at 0.75.
        let rows = [
            (1.0, 0.0, 1, "a"),
            (2.0, 0.0, 0, "a"),
            (3.0, 0.0, 0, "a"),
            (4.0, 0.0, 0, "a"),
        ];
        let (l, _) = run(json!({
            "input": layer(&rows), "age_field": "age", "event_field": "evt",
        }));
        for m in col_f64(&l, "km_median") {
            assert!(m.is_none(), "median should be Null (not reached)");
        }
        // S never below 0.75.
        for s in col_f64(&l, "km_survival") {
            assert!(s.unwrap() >= 0.75 - 1e-9);
        }
    }

    /// Stratification fits an independent curve per group, giving different
    /// medians. Group "fast" events early; group "slow" events late.
    #[test]
    fn stratification_yields_separate_medians() {
        let rows = [
            (1.0, 0.0, 1, "fast"),
            (2.0, 0.0, 1, "fast"),
            (3.0, 0.0, 1, "fast"),
            (10.0, 0.0, 1, "slow"),
            (20.0, 0.0, 1, "slow"),
            (30.0, 0.0, 1, "slow"),
        ];
        let (l, out) = run(json!({
            "input": layer(&rows), "age_field": "age", "event_field": "evt",
            "stratify_field": "grp",
        }));
        let strata = &out.outputs["strata"];
        let fast_med = strata["fast"]["median"].as_f64().unwrap();
        let slow_med = strata["slow"]["median"].as_f64().unwrap();
        assert!(fast_med < slow_med, "fast group must have a smaller median");
        // Median field on a feature matches its own stratum.
        let stratum_idx = l.schema.field_index("km_stratum").unwrap();
        let median = col_f64(&l, "km_median");
        for (i, feat) in l.features.iter().enumerate() {
            let g = feat.attributes[stratum_idx].as_str().unwrap();
            let expected = if g == "fast" { fast_med } else { slow_med };
            assert!((median[i].unwrap() - expected).abs() < 1e-9);
        }
        assert_eq!(out.outputs["stratum_count"].as_u64().unwrap(), 2);
    }

    /// Features with an invalid (negative / missing) duration are preserved but
    /// excluded from the curve and given Null outputs.
    #[test]
    fn invalid_duration_passes_through_as_null() {
        let rows = [
            (2.0, 0.0, 1, "a"),
            (-1.0, 0.0, 1, "a"), // negative duration -> excluded
            (6.0, 0.0, 1, "a"),
        ];
        let (l, out) = run(json!({
            "input": layer(&rows), "age_field": "age", "event_field": "evt",
        }));
        assert_eq!(l.features.len(), 3, "all features preserved");
        let surv = col_f64(&l, "km_survival");
        assert!(
            surv[1].is_none(),
            "invalid-duration feature -> Null survival"
        );
        assert!(col_f64(&l, "km_median")[1].is_none());
        // Only two valid events counted.
        assert_eq!(out.outputs["event_count"].as_u64().unwrap(), 2);
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            EstimateTimeToEventTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "age_field": "age" })).is_err());
        assert!(
            bad(json!({ "input": "a.geojson", "age_field": "age", "event_field": 3 })).is_err()
        );
        assert!(bad(json!({
            "input": "a.geojson", "age_field": "age", "event_field": "evt"
        }))
        .is_ok());
    }
}
