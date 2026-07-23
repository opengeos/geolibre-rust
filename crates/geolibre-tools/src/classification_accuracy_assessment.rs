//! GeoLibre tool: score a classified map against field-collected reference
//! points and report a confusion matrix with accuracy metrics.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Compute Confusion Matrix* (Spatial
//! Analyst), paired with *Update Accuracy Assessment Points*. The bundled
//! whitebox-wasm suite ships `kappa_index` and
//! `evaluate_object_classification_accuracy`, but those compare **two full
//! rasters** cell-by-cell — they need a completely labelled reference raster.
//! The real-world validation workflow is a **classified raster plus a set of
//! reference sample points**, each carrying a ground-truth class; that
//! point-based accuracy assessment is what this tool provides.
//!
//! For every reference point the tool determines the *predicted* (map) class —
//! either by sampling the classified raster at the point with a **nearest-cell**
//! lookup (classes are categorical, so no interpolation), or by reading a
//! `classified_field` that already carries the prediction — and compares it with
//! the ground-truth class in `class_field`. Class labels are treated as integers
//! (floats are rounded). Cross-tabulating truth × predicted over the sorted set
//! of observed classes yields the confusion matrix `M` (rows = reference/truth,
//! columns = classified/map), from which the tool computes:
//!
//! * `overall_accuracy` = Σ_c M[c][c] / total,
//! * `producers_accuracy[c]` = M[c][c] / (reference-row total for c),
//! * `users_accuracy[c]` = M[c][c] / (map-column total for c),
//! * Cohen's `kappa` = (p₀ − pₑ)/(1 − pₑ), with pₑ the chance-agreement sum.
//!
//! The reference points are echoed to the output layer annotated with
//! `REF_CLASS`, `MAP_CLASS`, and `CORRECT` (1/0). Points that fall outside the
//! raster, on no-data, or that lack a valid prediction get a null `MAP_CLASS`,
//! `CORRECT` = 0, and are excluded from the matrix.

use std::collections::{BTreeMap, BTreeSet};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;
use wbvector::{FieldDef, FieldType, FieldValue, Geometry};

use crate::common::load_input_raster;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct ClassificationAccuracyAssessmentTool;

impl Tool for ClassificationAccuracyAssessmentTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "classification_accuracy_assessment",
            display_name: "Classification Accuracy Assessment",
            summary: "Score a classified raster (or points that already carry a predicted class) against reference points labelled with ground truth, producing a confusion matrix, overall accuracy, per-class producer's/user's accuracy, and Cohen's kappa — like ArcGIS Compute Confusion Matrix.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "points",
                    description: "Reference (accuracy assessment) point layer carrying the ground-truth class field.",
                    required: true,
                },
                ToolParamSpec {
                    name: "class_field",
                    description: "Attribute on the points holding the ground-truth (reference) class label.",
                    required: true,
                },
                ToolParamSpec {
                    name: "input",
                    description: "Classified raster to sample at each reference point (nearest cell). Provide this OR 'classified_field'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "classified_field",
                    description: "Attribute on the points that already holds the predicted (map) class, used instead of sampling a raster. Provide this OR 'input'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band of the classified raster to sample (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path for the reference points annotated with REF_CLASS, MAP_CLASS, CORRECT. If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "points")?;
        require_str(args, "class_field")?;
        let raster = parse_optional_str(args, "input")?;
        let classified_field = parse_optional_str(args, "classified_field")?;
        if raster.is_none() && classified_field.is_none() {
            return Err(ToolError::Validation(
                "provide either 'input' (a classified raster) or 'classified_field'".to_string(),
            ));
        }
        parse_band(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let points_path = require_str(args, "points")?;
        let class_field = require_str(args, "class_field")?;
        let raster_path = parse_optional_str(args, "input")?;
        let classified_field = parse_optional_str(args, "classified_field")?;
        let output = parse_optional_str(args, "output")?;
        let band = parse_band(args)?;

        if raster_path.is_none() && classified_field.is_none() {
            return Err(ToolError::Validation(
                "provide either 'input' (a classified raster) or 'classified_field'".to_string(),
            ));
        }

        ctx.progress.info("reading reference points");
        let mut layer = load_input_layer(points_path)?;

        let truth_idx = layer.schema.field_index(class_field).ok_or_else(|| {
            ToolError::Validation(format!(
                "class_field '{class_field}' not found on the points"
            ))
        })?;

        // When sampling a raster, the raster takes precedence over a supplied
        // classified_field. Otherwise read the prediction straight from a field.
        let dem = match raster_path {
            Some(path) => {
                ctx.progress.info("reading classified raster");
                Some(load_input_raster(path)?)
            }
            None => None,
        };
        let map_idx = match (dem.is_some(), classified_field) {
            (true, _) => None,
            (false, Some(name)) => Some(layer.schema.field_index(name).ok_or_else(|| {
                ToolError::Validation(format!("classified_field '{name}' not found on the points"))
            })?),
            (false, None) => unreachable!("validated above"),
        };

        ctx.progress
            .info(&format!("scoring {} reference point(s)", layer.len()));

        // First pass: resolve (truth, predicted) for every feature.
        let mut resolved: Vec<Resolved> = Vec::with_capacity(layer.len());
        for feature in &layer.features {
            let truth = feature
                .attributes
                .get(truth_idx)
                .and_then(FieldValue::as_f64)
                .filter(|v| v.is_finite())
                .map(round_label);

            let predicted = match &dem {
                Some(r) => feature
                    .geometry
                    .as_ref()
                    .and_then(point_xy)
                    .and_then(|(x, y)| sample_class(r, band, x, y)),
                None => {
                    let idx = map_idx.expect("field index when no raster");
                    feature
                        .attributes
                        .get(idx)
                        .and_then(FieldValue::as_f64)
                        .filter(|v| v.is_finite())
                        .map(round_label)
                }
            };
            resolved.push(Resolved { truth, predicted });
        }

        // Sorted set of every class label observed among scored points.
        let mut class_set: BTreeSet<i64> = BTreeSet::new();
        for r in &resolved {
            if let (Some(t), Some(p)) = (r.truth, r.predicted) {
                class_set.insert(t);
                class_set.insert(p);
            }
        }
        let classes: Vec<i64> = class_set.into_iter().collect();
        let index_of: BTreeMap<i64, usize> =
            classes.iter().enumerate().map(|(i, &c)| (c, i)).collect();
        let n = classes.len();

        // Confusion matrix M[truth][predicted].
        let mut matrix = vec![vec![0u64; n]; n];
        let mut scored = 0usize;
        for r in &resolved {
            if let (Some(t), Some(p)) = (r.truth, r.predicted) {
                let (ti, pi) = (index_of[&t], index_of[&p]);
                matrix[ti][pi] += 1;
                scored += 1;
            }
        }

        let metrics = Metrics::compute(&matrix);

        // Annotate the reference points with REF_CLASS / MAP_CLASS / CORRECT.
        layer.add_field(FieldDef::new("REF_CLASS", FieldType::Integer));
        layer.add_field(FieldDef::new("MAP_CLASS", FieldType::Integer));
        layer.add_field(FieldDef::new("CORRECT", FieldType::Integer));
        for (feature, r) in layer.features.iter_mut().zip(resolved.iter()) {
            feature.attributes.push(match r.truth {
                Some(t) => FieldValue::Integer(t),
                None => FieldValue::Null,
            });
            let both_valid = r.truth.is_some() && r.predicted.is_some();
            feature.attributes.push(match (both_valid, r.predicted) {
                (true, Some(p)) => FieldValue::Integer(p),
                _ => FieldValue::Null,
            });
            let correct = matches!((r.truth, r.predicted), (Some(t), Some(p)) if t == p);
            feature
                .attributes
                .push(FieldValue::Integer(if correct { 1 } else { 0 }));
        }

        ctx.progress.info(&format!(
            "{scored} point(s) scored across {n} class(es); overall accuracy {:.4}",
            metrics.overall_accuracy
        ));

        let out_path = write_or_store_layer(layer, output)?;

        let confusion: Vec<Vec<u64>> = matrix.clone();
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("point_count".to_string(), json!(scored));
        outputs.insert(
            "overall_accuracy".to_string(),
            json!(metrics.overall_accuracy),
        );
        outputs.insert("kappa".to_string(), json!(metrics.kappa));
        outputs.insert("class_count".to_string(), json!(n));
        outputs.insert("classes".to_string(), json!(classes));
        outputs.insert("confusion_matrix".to_string(), json!(confusion));
        outputs.insert(
            "producers_accuracy".to_string(),
            json!(metrics.producers_accuracy),
        );
        outputs.insert("users_accuracy".to_string(), json!(metrics.users_accuracy));
        Ok(ToolRunResult { outputs })
    }
}

/// A per-feature resolution of the ground-truth and predicted class labels.
struct Resolved {
    truth: Option<i64>,
    predicted: Option<i64>,
}

/// Rounds a (possibly floating-point) class label to the nearest integer.
fn round_label(v: f64) -> i64 {
    v.round() as i64
}

/// Extracts the sampling coordinate of a point-like geometry.
fn point_xy(g: &Geometry) -> Option<(f64, f64)> {
    match g {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) => cs.first().map(|c| (c.x, c.y)),
        _ => None,
    }
}

/// Nearest-cell class lookup: returns the rounded class at `(x, y)`, or `None`
/// when the point lies outside the raster or on a no-data / NaN cell.
fn sample_class(raster: &Raster, band: isize, x: f64, y: f64) -> Option<i64> {
    let (col, row) = raster.world_to_pixel(x, y)?;
    if row < 0 || col < 0 || row >= raster.rows as isize || col >= raster.cols as isize {
        return None;
    }
    let v = raster.get(band, row, col);
    if v == raster.nodata || v.is_nan() {
        None
    } else {
        Some(round_label(v))
    }
}

/// Confusion-matrix-derived accuracy metrics.
struct Metrics {
    overall_accuracy: f64,
    kappa: f64,
    producers_accuracy: Vec<f64>,
    users_accuracy: Vec<f64>,
}

impl Metrics {
    /// Computes all metrics from `matrix[truth][predicted]`.
    fn compute(matrix: &[Vec<u64>]) -> Metrics {
        let n = matrix.len();
        let total: u64 = matrix.iter().flat_map(|r| r.iter()).sum();
        let total_f = total as f64;

        let row_totals: Vec<u64> = matrix.iter().map(|r| r.iter().sum()).collect();
        let mut col_totals = vec![0u64; n];
        for row in matrix {
            for (c, &v) in row.iter().enumerate() {
                col_totals[c] += v;
            }
        }

        let diagonal: u64 = (0..n).map(|i| matrix[i][i]).sum();
        let overall_accuracy = if total == 0 {
            0.0
        } else {
            diagonal as f64 / total_f
        };

        // Producer's accuracy: diagonal / reference-row total.
        let producers_accuracy: Vec<f64> = (0..n)
            .map(|i| {
                if row_totals[i] == 0 {
                    0.0
                } else {
                    matrix[i][i] as f64 / row_totals[i] as f64
                }
            })
            .collect();

        // User's accuracy: diagonal / map-column total.
        let users_accuracy: Vec<f64> = (0..n)
            .map(|i| {
                if col_totals[i] == 0 {
                    0.0
                } else {
                    matrix[i][i] as f64 / col_totals[i] as f64
                }
            })
            .collect();

        // Cohen's kappa.
        let p0 = overall_accuracy;
        let pe = if total == 0 {
            0.0
        } else {
            (0..n)
                .map(|i| (row_totals[i] as f64 / total_f) * (col_totals[i] as f64 / total_f))
                .sum::<f64>()
        };
        let kappa = if (1.0 - pe).abs() < f64::EPSILON {
            0.0
        } else {
            (p0 - pe) / (1.0 - pe)
        };

        Metrics {
            overall_accuracy,
            kappa,
            producers_accuracy,
            users_accuracy,
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

/// Parses the optional 1-based `band` (default 1) into a 0-based `isize`.
fn parse_band(args: &ToolArgs) -> Result<isize, ToolError> {
    let band_1based = match args.get("band") {
        None | Some(Value::Null) => 1,
        Some(Value::Number(num)) => num
            .as_i64()
            .or_else(|| num.as_f64().map(|f| f as i64))
            .ok_or_else(|| ToolError::Validation("'band' must be an integer".to_string()))?,
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => s
            .trim()
            .parse::<i64>()
            .map_err(|_| ToolError::Validation("'band' must be an integer".to_string()))?,
        Some(_) => {
            return Err(ToolError::Validation(
                "'band' must be an integer".to_string(),
            ))
        }
    };
    if band_1based < 1 {
        return Err(ToolError::Validation("'band' must be >= 1".to_string()));
    }
    Ok((band_1based - 1) as isize)
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
    use wbvector::{memory_store, Coord, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// A 4x4 classified raster split into two class halves: the left two columns
    /// are class 1, the right two columns are class 2 (cell size 1, origin 0,0).
    fn two_class_raster() -> String {
        let mut r = Raster::new(RasterConfig {
            cols: 4,
            rows: 4,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for row in 0..4 {
            for col in 0..4 {
                let class = if col < 2 { 1.0 } else { 2.0 };
                r.set(0, row as isize, col as isize, class).unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    /// Point layer with a ground-truth field plus an optional predicted field.
    fn points_layer(rows: &[(f64, f64, i64, Option<i64>)], with_predicted: bool) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("truth", FieldType::Integer));
        if with_predicted {
            l.add_field(FieldDef::new("pred", FieldType::Integer));
        }
        for &(x, y, truth, pred) in rows {
            let mut attrs: Vec<(&str, FieldValue)> = vec![("truth", FieldValue::Integer(truth))];
            if with_predicted {
                attrs.push((
                    "pred",
                    pred.map(FieldValue::Integer).unwrap_or(FieldValue::Null),
                ));
            }
            l.add_feature(Some(Geometry::Point(Coord::xy(x, y))), &attrs)
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        ClassificationAccuracyAssessmentTool
            .run(&args, &ctx())
            .unwrap()
    }

    // Cell centres: col c centre x = c + 0.5; row r centre y = y_max - (r+0.5).
    // y_max = 4. Left half (col 0,1) => class 1; right half (col 2,3) => class 2.

    #[test]
    fn all_correct_raster_gives_perfect_scores() {
        let raster = two_class_raster();
        // Truth matches the map everywhere: left points labelled 1, right 2.
        let points = points_layer(
            &[
                (0.5, 3.5, 1, None), // col 0 -> class 1
                (1.5, 2.5, 1, None), // col 1 -> class 1
                (2.5, 1.5, 2, None), // col 2 -> class 2
                (3.5, 0.5, 2, None), // col 3 -> class 2
            ],
            false,
        );
        let out = run(json!({
            "points": points, "class_field": "truth", "input": raster,
        }));
        assert_eq!(out.outputs["point_count"].as_u64().unwrap(), 4);
        assert!((out.outputs["overall_accuracy"].as_f64().unwrap() - 1.0).abs() < 1e-12);
        assert!((out.outputs["kappa"].as_f64().unwrap() - 1.0).abs() < 1e-12);
        assert_eq!(out.outputs["class_count"].as_u64().unwrap(), 2);
        // Confusion matrix is diagonal: [[2,0],[0,2]].
        let m = &out.outputs["confusion_matrix"];
        assert_eq!(m[0][0].as_u64().unwrap(), 2);
        assert_eq!(m[0][1].as_u64().unwrap(), 0);
        assert_eq!(m[1][0].as_u64().unwrap(), 0);
        assert_eq!(m[1][1].as_u64().unwrap(), 2);
    }

    #[test]
    fn two_errors_give_exact_fraction() {
        let raster = two_class_raster();
        // Four points; two mislabelled relative to the map class.
        //   map classes at the sampled cells: 1,1,2,2
        //   truth given:                      1,2,2,1   -> 2 correct of 4
        let points = points_layer(
            &[
                (0.5, 3.5, 1, None), // map 1, truth 1  correct
                (1.5, 2.5, 2, None), // map 1, truth 2  wrong
                (2.5, 1.5, 2, None), // map 2, truth 2  correct
                (3.5, 0.5, 1, None), // map 2, truth 1  wrong
            ],
            false,
        );
        let out = run(json!({
            "points": points, "class_field": "truth", "input": raster,
        }));
        assert_eq!(out.outputs["point_count"].as_u64().unwrap(), 4);
        assert!((out.outputs["overall_accuracy"].as_f64().unwrap() - 0.5).abs() < 1e-12);
        // Matrix rows=truth, cols=map. truth1: {map1:1, map2:1}; truth2:{map1:1,map2:1}.
        // row totals 2,2; col totals 2,2; pe = 0.25+0.25 = 0.5; p0=0.5 => kappa 0.
        assert!((out.outputs["kappa"].as_f64().unwrap()).abs() < 1e-12);
    }

    #[test]
    fn classified_field_path_matches_metrics() {
        // No raster: predictions carried on the points themselves.
        let points = points_layer(
            &[
                (0.0, 0.0, 1, Some(1)),
                (0.0, 0.0, 1, Some(2)),
                (0.0, 0.0, 2, Some(2)),
                (0.0, 0.0, 2, Some(1)),
            ],
            true,
        );
        let out = run(json!({
            "points": points, "class_field": "truth", "classified_field": "pred",
        }));
        assert_eq!(out.outputs["point_count"].as_u64().unwrap(), 4);
        assert!((out.outputs["overall_accuracy"].as_f64().unwrap() - 0.5).abs() < 1e-12);
        assert!((out.outputs["kappa"].as_f64().unwrap()).abs() < 1e-12);
    }

    #[test]
    fn per_class_producer_and_user_accuracy() {
        // Asymmetric error case, three classes, using the classified_field path.
        // Matrix (rows=truth, cols=map):
        //   truth1: map1 x3, map2 x1        row total 4
        //   truth2: map2 x2                  row total 2
        //   truth3: map1 x1, map3 x2         row total 3
        // col totals: map1 = 3+1 = 4, map2 = 1+2 = 3, map3 = 2
        let mut rows: Vec<(f64, f64, i64, Option<i64>)> = Vec::new();
        for _ in 0..3 {
            rows.push((0.0, 0.0, 1, Some(1)));
        }
        rows.push((0.0, 0.0, 1, Some(2)));
        rows.push((0.0, 0.0, 2, Some(2)));
        rows.push((0.0, 0.0, 2, Some(2)));
        rows.push((0.0, 0.0, 3, Some(1)));
        rows.push((0.0, 0.0, 3, Some(3)));
        rows.push((0.0, 0.0, 3, Some(3)));
        let points = points_layer(&rows, true);
        let out = run(json!({
            "points": points, "class_field": "truth", "classified_field": "pred",
        }));
        let classes = out.outputs["classes"].as_array().unwrap();
        assert_eq!(classes.len(), 3);
        // Producer's accuracy = diagonal / row total.
        let pa = out.outputs["producers_accuracy"].as_array().unwrap();
        assert!((pa[0].as_f64().unwrap() - 3.0 / 4.0).abs() < 1e-12); // class 1
        assert!((pa[1].as_f64().unwrap() - 2.0 / 2.0).abs() < 1e-12); // class 2
        assert!((pa[2].as_f64().unwrap() - 2.0 / 3.0).abs() < 1e-12); // class 3
                                                                      // User's accuracy = diagonal / col total.
        let ua = out.outputs["users_accuracy"].as_array().unwrap();
        assert!((ua[0].as_f64().unwrap() - 3.0 / 4.0).abs() < 1e-12); // map 1
        assert!((ua[1].as_f64().unwrap() - 2.0 / 3.0).abs() < 1e-12); // map 2
        assert!((ua[2].as_f64().unwrap() - 2.0 / 2.0).abs() < 1e-12); // map 3
                                                                      // Overall = diagonal sum / total = (3+2+2)/9.
        assert!((out.outputs["overall_accuracy"].as_f64().unwrap() - 7.0 / 9.0).abs() < 1e-12);
    }

    #[test]
    fn annotates_points_and_excludes_out_of_range() {
        let raster = two_class_raster();
        // Third point is outside the 4x4 raster => no prediction, excluded.
        let points = points_layer(
            &[
                (0.5, 3.5, 1, None),
                (2.5, 1.5, 2, None),
                (99.0, 99.0, 1, None),
            ],
            false,
        );
        let out = run(json!({
            "points": points, "class_field": "truth", "input": raster,
        }));
        assert_eq!(out.outputs["point_count"].as_u64().unwrap(), 2);
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        let ref_idx = layer.schema.field_index("REF_CLASS").unwrap();
        let map_idx = layer.schema.field_index("MAP_CLASS").unwrap();
        let cor_idx = layer.schema.field_index("CORRECT").unwrap();
        // First two scored & correct.
        assert_eq!(layer.features[0].attributes[cor_idx].as_i64(), Some(1));
        assert_eq!(layer.features[1].attributes[cor_idx].as_i64(), Some(2 - 1)); // 1
        assert_eq!(layer.features[0].attributes[map_idx].as_i64(), Some(1));
        // Out-of-range point: MAP_CLASS null, CORRECT 0, REF_CLASS still set.
        assert_eq!(layer.features[2].attributes[map_idx], FieldValue::Null);
        assert_eq!(layer.features[2].attributes[cor_idx].as_i64(), Some(0));
        assert_eq!(layer.features[2].attributes[ref_idx].as_i64(), Some(1));
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ClassificationAccuracyAssessmentTool.validate(&args)
        };
        // Missing points.
        assert!(bad(json!({ "class_field": "truth", "input": "r.tif" })).is_err());
        // Missing class_field.
        assert!(bad(json!({ "points": "p.geojson", "input": "r.tif" })).is_err());
        // Neither input nor classified_field.
        assert!(bad(json!({ "points": "p.geojson", "class_field": "truth" })).is_err());
        // Valid: raster path.
        assert!(
            bad(json!({ "points": "p.geojson", "class_field": "truth", "input": "r.tif" })).is_ok()
        );
        // Valid: classified_field path.
        assert!(bad(
            json!({ "points": "p.geojson", "class_field": "truth", "classified_field": "pred" })
        )
        .is_ok());
    }
}
