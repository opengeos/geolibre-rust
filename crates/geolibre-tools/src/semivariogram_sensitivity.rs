//! GeoLibre tool: how sensitive kriging predictions are to variogram parameters.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *GA Semivariogram Sensitivity*
//! (Geostatistical Analyst).
//!
//! ## Fit is not the same question as sensitivity
//!
//! The bundled kriging suite is broad — `ordinary_kriging`, `universal_kriging`,
//! `simple_kriging`, `ordinary_cokriging`, `empirical_bayesian_kriging`,
//! `local_kriging`, `spacetime_kriging` — and `kriging_cross_validation`
//! measures **fit**. Nothing measures **sensitivity**: how much the answer moves
//! when the analyst's variogram choices move.
//!
//! Those are independent properties. A model can cross-validate beautifully and
//! still be so parameter-sensitive that a defensible alternative nugget changes
//! the prediction by half its range — which is exactly what a reviewer asks
//! about, and what the analyst currently cannot answer.
//!
//! ## Design
//!
//! A full factorial over perturbed (nugget, partial sill, range) triples, each
//! re-kriged at the evaluation locations, tabulated against the baseline fit.
//! The grid is **enumerated, not sampled**, so there is no RNG and results are
//! reproducible — a hard requirement for the WASM paths.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::args_common::{choice_or, f64_or, req_str, usize_or};
use crate::kriging_common::{fit_variogram, ordinary_kriging, Variogram, VariogramModel};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct SemivariogramSensitivityTool;

impl Tool for SemivariogramSensitivityTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "semivariogram_sensitivity",
            display_name: "Semivariogram Sensitivity",
            summary: "Perturbs the nugget, partial sill and range of a fitted variogram, re-krige at evaluation locations, and reports how far predictions and standard errors move — the standard defensibility check on a kriging model (ArcGIS GA Semivariogram Sensitivity). The bundled kriging suite is broad and kriging_cross_validation measures fit, but nothing measures sensitivity, and the two are independent: a well-cross-validated model can still be highly fragile to the parameters the analyst chose.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Point layer holding the measured values.",
                    required: true,
                },
                ToolParamSpec {
                    name: "field",
                    description: "Numeric field to krige.",
                    required: true,
                },
                ToolParamSpec {
                    name: "locations",
                    description: "Point layer of prediction locations to evaluate the sensitivity at.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Table with one row per perturbed model. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "model",
                    description: "Variogram model: 'exponential' (default), 'spherical', or 'gaussian'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "nugget_span_percent",
                    description: "Perturbation half-width for the nugget, in percent (default 10).",
                    required: false,
                },
                ToolParamSpec {
                    name: "partial_sill_span_percent",
                    description: "Perturbation half-width for the partial sill, in percent (default 10).",
                    required: false,
                },
                ToolParamSpec {
                    name: "range_span_percent",
                    description: "Perturbation half-width for the range, in percent (default 10).",
                    required: false,
                },
                ToolParamSpec {
                    name: "nugget_steps",
                    description: "Samples across the nugget span (default 3; 1 pins it at the fitted value).",
                    required: false,
                },
                ToolParamSpec {
                    name: "partial_sill_steps",
                    description: "Samples across the partial-sill span (default 3).",
                    required: false,
                },
                ToolParamSpec {
                    name: "range_steps",
                    description: "Samples across the range span (default 3).",
                    required: false,
                },
                ToolParamSpec {
                    name: "lag_count",
                    description: "Empirical-variogram lag bins used for the baseline fit (default 12).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        req_str(args, "field")?;
        req_str(args, "locations")?;
        VariogramModel::parse(choice_or(
            args,
            "model",
            &["exponential", "spherical", "gaussian"],
            "exponential",
        )?)?;
        for k in [
            "nugget_span_percent",
            "partial_sill_span_percent",
            "range_span_percent",
        ] {
            let v = f64_or(args, k, 10.0)?;
            if !(0.0..100.0).contains(&v) {
                return Err(ToolError::Validation(format!("'{k}' must be in [0, 100)")));
            }
        }
        for k in ["nugget_steps", "partial_sill_steps", "range_steps"] {
            if usize_or(args, k, 3)? == 0 {
                return Err(ToolError::Validation(format!("'{k}' must be at least 1")));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let field = req_str(args, "field")?.to_string();
        let model = VariogramModel::parse(choice_or(
            args,
            "model",
            &["exponential", "spherical", "gaussian"],
            "exponential",
        )?)?;
        let output = parse_optional_str(args, "output")?;
        let lag_count = usize_or(args, "lag_count", 12)?;

        let (coords, values) = load_points(req_str(args, "input")?, Some(&field))?;
        if coords.len() < 3 {
            return Err(ToolError::Execution(format!(
                "need at least 3 measured points to fit a variogram, got {}",
                coords.len()
            )));
        }
        let (targets, _) = load_points(req_str(args, "locations")?, None)?;
        if targets.is_empty() {
            return Err(ToolError::Execution(
                "'locations' holds no point geometry".to_string(),
            ));
        }

        let baseline = fit_variogram(&coords, &values, model, lag_count);
        ctx.progress.info(&format!(
            "baseline {} variogram: nugget {:.4}, partial sill {:.4}, range {:.4}",
            model.label(),
            baseline.nugget,
            baseline.partial_sill,
            baseline.range
        ));

        // Baseline predictions, the reference every perturbation is measured
        // against.
        let base_pred: Vec<Option<(f64, f64)>> = targets
            .iter()
            .map(|t| ordinary_kriging(&coords, &values, *t, &baseline))
            .collect();

        // Absolute scale for dimensions whose baseline is zero.
        let fallback = (baseline.partial_sill + baseline.nugget).max(1e-12);
        let grid = |value: f64, span_key: &str, step_key: &str| -> Result<Vec<f64>, ToolError> {
            let span = f64_or(args, span_key, 10.0)? / 100.0;
            let steps = usize_or(args, step_key, 3)?;
            if steps == 1 {
                return Ok(vec![value]);
            }
            // A relative span collapses when the baseline is zero (a fitted
            // nugget of 0 is common), which would silently drop that dimension
            // from the factorial. Fall back to an absolute span scaled off the
            // baseline sill so the perturbation stays meaningful.
            let half = if value.abs() > 1e-12 {
                value * span
            } else {
                fallback * span
            };
            Ok((0..steps)
                .map(|i| {
                    let t = i as f64 / (steps - 1) as f64; // 0..1
                    (value - half + 2.0 * half * t).max(0.0)
                })
                .collect())
        };
        let nuggets = grid(baseline.nugget, "nugget_span_percent", "nugget_steps")?;
        let sills = grid(
            baseline.partial_sill,
            "partial_sill_span_percent",
            "partial_sill_steps",
        )?;
        let ranges = grid(baseline.range, "range_span_percent", "range_steps")?;

        let mut out = Layer::new("semivariogram_sensitivity");
        out.add_field(FieldDef::new("NUGGET", FieldType::Float));
        out.add_field(FieldDef::new("PARTIAL_SILL", FieldType::Float));
        out.add_field(FieldDef::new("RANGE", FieldType::Float));
        out.add_field(FieldDef::new("IS_BASELINE", FieldType::Boolean));
        out.add_field(FieldDef::new("MEAN_ABS_DELTA", FieldType::Float));
        out.add_field(FieldDef::new("MAX_ABS_DELTA", FieldType::Float));
        out.add_field(FieldDef::new("MEAN_ABS_SE_DELTA", FieldType::Float));
        out.add_field(FieldDef::new("MEAN_PREDICTION", FieldType::Float));
        out.add_field(FieldDef::new("RESOLVED_LOCATIONS", FieldType::Integer));

        let total_models = nuggets.len() * sills.len() * ranges.len();
        let mut worst_delta = 0.0_f64;
        let mut emitted = 0_u64;
        let mut n_done = 0usize;

        for &nugget in &nuggets {
            for &partial_sill in &sills {
                for &range in &ranges {
                    let vg = Variogram {
                        model,
                        nugget: nugget.max(0.0),
                        partial_sill: partial_sill.max(1e-12),
                        range: range.max(1e-12),
                    };
                    let mut sum_delta = 0.0;
                    let mut max_delta = 0.0_f64;
                    let mut sum_se_delta = 0.0;
                    let mut sum_pred = 0.0;
                    let mut resolved = 0usize;
                    for (k, t) in targets.iter().enumerate() {
                        let (Some((bp, bv)), Some((p, v))) =
                            (base_pred[k], ordinary_kriging(&coords, &values, *t, &vg))
                        else {
                            continue;
                        };
                        let d = (p - bp).abs();
                        sum_delta += d;
                        max_delta = max_delta.max(d);
                        sum_se_delta += (v.sqrt() - bv.sqrt()).abs();
                        sum_pred += p;
                        resolved += 1;
                    }
                    let denom = resolved.max(1) as f64;
                    worst_delta = worst_delta.max(max_delta);

                    // Floating-point identity: the baseline triple appears in
                    // the grid whenever every span is centred, and flagging it
                    // lets callers find the reference row without recomputing.
                    let is_baseline = (nugget - baseline.nugget).abs() <= f64::EPSILON
                        && (partial_sill - baseline.partial_sill).abs() <= f64::EPSILON
                        && (range - baseline.range).abs() <= f64::EPSILON;

                    out.add_feature(
                        None,
                        &[
                            ("NUGGET", FieldValue::Float(vg.nugget)),
                            ("PARTIAL_SILL", FieldValue::Float(vg.partial_sill)),
                            ("RANGE", FieldValue::Float(vg.range)),
                            ("IS_BASELINE", FieldValue::Boolean(is_baseline)),
                            ("MEAN_ABS_DELTA", FieldValue::Float(sum_delta / denom)),
                            ("MAX_ABS_DELTA", FieldValue::Float(max_delta)),
                            ("MEAN_ABS_SE_DELTA", FieldValue::Float(sum_se_delta / denom)),
                            ("MEAN_PREDICTION", FieldValue::Float(sum_pred / denom)),
                            ("RESOLVED_LOCATIONS", FieldValue::Integer(resolved as i64)),
                        ],
                    )
                    .map_err(|e| ToolError::Execution(e.to_string()))?;
                    emitted += 1;
                    n_done += 1;
                    ctx.progress
                        .progress(n_done as f64 / total_models.max(1) as f64);
                }
            }
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("model_count".to_string(), json!(emitted));
        outputs.insert("location_count".to_string(), json!(targets.len()));
        outputs.insert("worst_abs_delta".to_string(), json!(worst_delta));
        outputs.insert("baseline_nugget".to_string(), json!(baseline.nugget));
        outputs.insert(
            "baseline_partial_sill".to_string(),
            json!(baseline.partial_sill),
        );
        outputs.insert("baseline_range".to_string(), json!(baseline.range));
        outputs.insert("model".to_string(), json!(model.label()));
        Ok(ToolRunResult { outputs })
    }
}

/// Reads point coordinates and, optionally, a numeric field.
pub(crate) fn load_points(
    path: &str,
    field: Option<&str>,
) -> Result<(Vec<(f64, f64)>, Vec<f64>), ToolError> {
    let layer = load_input_layer(path)?;
    let idx = match field {
        Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
            ToolError::Validation(format!("field '{f}' not found on layer '{path}'"))
        })?),
        None => None,
    };
    let mut coords = Vec::new();
    let mut values = Vec::new();
    for feature in layer.iter() {
        let Some(Geometry::Point(c)) = feature.geometry.as_ref() else {
            continue;
        };
        if let Some(i) = idx {
            let Some(v) = numeric(&feature.attributes[i]) else {
                // A non-numeric or null measurement cannot enter the solve;
                // skipping is better than substituting a zero.
                continue;
            };
            values.push(v);
        }
        coords.push((c.x, c.y));
    }
    Ok((coords, values))
}

fn numeric(v: &FieldValue) -> Option<f64> {
    match v {
        FieldValue::Float(f) => f.is_finite().then_some(*f),
        FieldValue::Integer(i) => Some(*i as f64),
        FieldValue::Text(s) => s.trim().parse().ok(),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Coord, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn measured(pts: Vec<(f64, f64, f64)>) -> String {
        let mut l = Layer::new("obs");
        l.geom_type = Some(GeometryType::Point);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (x, y, v) in pts {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(x, y))),
                &[("v", FieldValue::Float(v))],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn locations(pts: Vec<(f64, f64)>) -> String {
        let mut l = Layer::new("loc");
        l.geom_type = Some(GeometryType::Point);
        for (x, y) in pts {
            l.add_feature(Some(Geometry::Point(Coord::xy(x, y))), &[])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    /// A smooth deterministic field on a 6x6 grid.
    fn sample_field() -> String {
        let mut pts = Vec::new();
        for i in 0..6 {
            for j in 0..6 {
                let (x, y) = (i as f64 * 10.0, j as f64 * 10.0);
                pts.push((x, y, (x / 25.0).sin() * 3.0 + (y / 25.0).cos() * 2.0));
            }
        }
        measured(pts)
    }

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = SemivariogramSensitivityTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    fn col(layer: &Layer, name: &str) -> Vec<f64> {
        let i = layer.schema.field_index(name).unwrap();
        layer
            .iter()
            .map(|f| match &f.attributes[i] {
                FieldValue::Float(v) => *v,
                FieldValue::Integer(v) => *v as f64,
                other => panic!("expected a number in {name}, got {other:?}"),
            })
            .collect()
    }

    #[test]
    fn the_factorial_has_one_row_per_parameter_triple() {
        let (out, res) = run(json!({
            "input": sample_field(),
            "field": "v",
            "locations": locations(vec![(15.0, 15.0), (35.0, 25.0)]),
            "nugget_steps": 2,
            "partial_sill_steps": 3,
            "range_steps": 4,
        }));
        assert_eq!(res.outputs["model_count"], json!(24));
        assert_eq!(out.iter().count(), 24);
    }

    #[test]
    fn the_baseline_row_has_zero_delta_from_itself() {
        // With odd step counts the centre of each span IS the fitted value.
        let (out, _) = run(json!({
            "input": sample_field(),
            "field": "v",
            "locations": locations(vec![(15.0, 15.0), (35.0, 25.0)]),
        }));
        let base = out.schema.field_index("IS_BASELINE").unwrap();
        let deltas = col(&out, "MEAN_ABS_DELTA");
        let mut found = false;
        for (k, f) in out.iter().enumerate() {
            if f.attributes[base] == FieldValue::Boolean(true) {
                assert!(deltas[k].abs() < 1e-9, "baseline delta {}", deltas[k]);
                found = true;
            }
        }
        assert!(found, "no baseline row was flagged");
    }

    #[test]
    fn a_wider_span_produces_a_larger_worst_case_delta() {
        // The defining behaviour: sensitivity grows with the perturbation.
        let go = |span: f64| {
            let (_, res) = run(json!({
                "input": sample_field(),
                "field": "v",
                "locations": locations(vec![(15.0, 15.0), (35.0, 25.0), (52.0, 8.0)]),
                "range_span_percent": span,
                "nugget_span_percent": span,
                "partial_sill_span_percent": span,
            }));
            res.outputs["worst_abs_delta"].as_f64().unwrap()
        };
        let narrow = go(5.0);
        let wide = go(60.0);
        assert!(wide > narrow, "narrow {narrow} vs wide {wide}");
    }

    #[test]
    fn a_zero_span_collapses_the_grid_to_the_baseline() {
        let (out, res) = run(json!({
            "input": sample_field(),
            "field": "v",
            "locations": locations(vec![(15.0, 15.0)]),
            "nugget_span_percent": 0.0,
            "partial_sill_span_percent": 0.0,
            "range_span_percent": 0.0,
        }));
        // Every row is the same model, so nothing moves.
        assert!(col(&out, "MEAN_ABS_DELTA").iter().all(|d| d.abs() < 1e-9));
        assert!(res.outputs["worst_abs_delta"].as_f64().unwrap() < 1e-9);
    }

    #[test]
    fn single_steps_pin_the_parameters_and_give_one_row() {
        let (_, res) = run(json!({
            "input": sample_field(),
            "field": "v",
            "locations": locations(vec![(15.0, 15.0)]),
            "nugget_steps": 1,
            "partial_sill_steps": 1,
            "range_steps": 1,
        }));
        assert_eq!(res.outputs["model_count"], json!(1));
    }

    #[test]
    fn the_run_is_deterministic() {
        // No RNG anywhere: the factorial is enumerated, not sampled.
        let go = || {
            let (_, res) = run(json!({
                "input": sample_field(),
                "field": "v",
                "locations": locations(vec![(15.0, 15.0), (35.0, 25.0)]),
            }));
            res.outputs["worst_abs_delta"].as_f64().unwrap()
        };
        assert_eq!(go(), go());
    }

    #[test]
    fn every_perturbed_model_stays_physically_valid() {
        let (out, _) = run(json!({
            "input": sample_field(),
            "field": "v",
            "locations": locations(vec![(15.0, 15.0)]),
            "nugget_span_percent": 90.0,
            "partial_sill_span_percent": 90.0,
            "range_span_percent": 90.0,
        }));
        assert!(col(&out, "NUGGET").iter().all(|v| *v >= 0.0));
        assert!(col(&out, "PARTIAL_SILL").iter().all(|v| *v > 0.0));
        assert!(col(&out, "RANGE").iter().all(|v| *v > 0.0));
        assert!(col(&out, "MEAN_ABS_SE_DELTA").iter().all(|v| *v >= 0.0));
    }

    #[test]
    fn rejects_bad_parameters() {
        let obs = sample_field();
        let loc = locations(vec![(1.0, 1.0)]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SemivariogramSensitivityTool.validate(&args).is_err()
        };
        assert!(bad(json!({"field": "v", "locations": loc})));
        assert!(bad(json!({"input": obs, "locations": loc})));
        assert!(bad(
            json!({"input": obs, "field": "v", "locations": loc, "model": "nope"})
        ));
        assert!(bad(
            json!({"input": obs, "field": "v", "locations": loc, "range_span_percent": 150})
        ));
        assert!(bad(
            json!({"input": obs, "field": "v", "locations": loc, "range_steps": 0})
        ));
    }

    #[test]
    fn a_missing_field_is_reported_at_run_time() {
        let args: ToolArgs = serde_json::from_value(json!({
            "input": sample_field(),
            "field": "nope",
            "locations": locations(vec![(1.0, 1.0)]),
        }))
        .unwrap();
        assert!(SemivariogramSensitivityTool.run(&args, &ctx()).is_err());
    }
}
