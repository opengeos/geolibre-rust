//! GeoLibre tool: score object detections against ground truth (precision,
//! recall, average precision, mAP).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Compute Accuracy For Object Detection*
//! (Image Analyst). The repo ships both neighbours of this step —
//! `non_maximum_suppression` cleans up overlapping detections, and
//! `classification_accuracy_assessment` scores a classified raster against
//! reference points — but neither scores *object detections*: instances with a
//! class label and a confidence score, matched to ground truth by overlap.
//! Without it there is no way to tell whether a detection run improved.
//!
//! `classification_accuracy_assessment` cannot be reused: it builds a
//! pixel-level confusion matrix and kappa, which has no notion of instance
//! matching, an IoU threshold, or confidence-ranked average precision. The
//! bundled suite has no object-detection metrics at all.
//!
//! Matching follows the standard detection protocol. Within each class,
//! detections are sorted by descending confidence and greedily matched to the
//! highest-IoU **unmatched** ground-truth instance above `min_iou`. Matched
//! detections are true positives, unmatched are false positives, and unmatched
//! truths are false negatives. Average precision integrates the
//! precision-recall curve in that ranked order using all-point interpolation
//! (the VOC-2010 / COCO convention), which is why confidence order matters and
//! a plain precision/recall pair would be insufficient.

use std::collections::BTreeMap;

use geo::{Area, BooleanOps, BoundingRect, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct ComputeAccuracyForObjectDetectionTool;

impl Tool for ComputeAccuracyForObjectDetectionTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "compute_accuracy_for_object_detection",
            display_name: "Compute Accuracy For Object Detection",
            summary: "Score detected features against ground truth using IoU matching, reporting per-class precision, recall, average precision and overall mAP, like ArcGIS Compute Accuracy For Object Detection.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "detected",
                    description: "Detected features (polygons or bounding boxes).",
                    required: true,
                },
                ToolParamSpec {
                    name: "ground_truth",
                    description: "Reference features to score against.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output accuracy table path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "detected_class_field",
                    description: "Class label field on the detections. Omit to score everything as a single class.",
                    required: false,
                },
                ToolParamSpec {
                    name: "ground_truth_class_field",
                    description: "Class label field on the ground truth. Omit to score everything as a single class.",
                    required: false,
                },
                ToolParamSpec {
                    name: "confidence_field",
                    description: "Confidence score field on the detections, used to rank for average precision. Without it detections are scored in input order.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_iou",
                    description: "Minimum intersection-over-union for a detection to count as a match (default 0.5).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "detected")?;
        require_str(args, "ground_truth")?;
        if let Some(t) = parse_optional_f64(args, "min_iou")? {
            if !t.is_finite() || !(0.0..=1.0).contains(&t) {
                return Err(ToolError::Validation(
                    "'min_iou' must be between 0 and 1".to_string(),
                ));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let det_path = require_str(args, "detected")?;
        let gt_path = require_str(args, "ground_truth")?;
        let output = parse_optional_str(args, "output")?;
        let det_class_field = parse_optional_str(args, "detected_class_field")?.map(String::from);
        let gt_class_field =
            parse_optional_str(args, "ground_truth_class_field")?.map(String::from);
        let conf_field = parse_optional_str(args, "confidence_field")?.map(String::from);
        let min_iou = parse_optional_f64(args, "min_iou")?.unwrap_or(0.5);

        let det_layer = load_input_layer(det_path)?;
        let gt_layer = load_input_layer(gt_path)?;

        for (label, layer, field) in [
            ("detected_class_field", &det_layer, &det_class_field),
            ("ground_truth_class_field", &gt_layer, &gt_class_field),
            ("confidence_field", &det_layer, &conf_field),
        ] {
            if let Some(f) = field {
                if layer.schema.field_index(f).is_none() {
                    return Err(ToolError::Validation(format!(
                        "{label} '{f}' not found on the layer"
                    )));
                }
            }
        }

        // Extract detections with class + confidence, and truths with class.
        let mut dets: Vec<Inst> = det_layer
            .features
            .iter()
            .filter_map(|f| {
                let mp = to_multipolygon(f.geometry.as_ref()?)?;
                let class = det_class_field
                    .as_ref()
                    .map(|c| {
                        f.get(&det_layer.schema, c)
                            .map(field_string)
                            .unwrap_or_default()
                    })
                    .unwrap_or_else(|| "1".to_string());
                let conf = conf_field
                    .as_ref()
                    .and_then(|c| f.get(&det_layer.schema, c).ok())
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0);
                Some(Inst {
                    area: mp.unsigned_area(),
                    bbox: mp.bounding_rect(),
                    geom: mp,
                    class,
                    conf,
                })
            })
            .collect();

        let truths: Vec<Inst> = gt_layer
            .features
            .iter()
            .filter_map(|f| {
                let mp = to_multipolygon(f.geometry.as_ref()?)?;
                let class = gt_class_field
                    .as_ref()
                    .map(|c| {
                        f.get(&gt_layer.schema, c)
                            .map(field_string)
                            .unwrap_or_default()
                    })
                    .unwrap_or_else(|| "1".to_string());
                Some(Inst {
                    area: mp.unsigned_area(),
                    bbox: mp.bounding_rect(),
                    geom: mp,
                    class,
                    conf: 0.0,
                })
            })
            .collect();

        // Highest confidence first; stable so ties keep input order.
        dets.sort_by(|a, b| {
            b.conf
                .partial_cmp(&a.conf)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        ctx.progress.info(&format!(
            "{} detection(s) vs {} truth(s) at IoU >= {min_iou}",
            dets.len(),
            truths.len()
        ));

        let mut classes: Vec<String> = dets
            .iter()
            .map(|d| d.class.clone())
            .chain(truths.iter().map(|t| t.class.clone()))
            .collect();
        classes.sort();
        classes.dedup();

        let mut out = Layer::new("accuracy");
        out.add_field(FieldDef::new("class", FieldType::Text));
        out.add_field(FieldDef::new("truths", FieldType::Integer));
        out.add_field(FieldDef::new("detections", FieldType::Integer));
        out.add_field(FieldDef::new("tp", FieldType::Integer));
        out.add_field(FieldDef::new("fp", FieldType::Integer));
        out.add_field(FieldDef::new("fn", FieldType::Integer));
        out.add_field(FieldDef::new("precision", FieldType::Float));
        out.add_field(FieldDef::new("recall", FieldType::Float));
        out.add_field(FieldDef::new("ap", FieldType::Float));

        let mut ap_sum = 0.0;
        let mut ap_classes = 0usize;
        let mut tp_all = 0usize;
        let mut fp_all = 0usize;
        let mut fn_all = 0usize;

        for class in &classes {
            let cd: Vec<&Inst> = dets.iter().filter(|d| &d.class == class).collect();
            let ct: Vec<&Inst> = truths.iter().filter(|t| &t.class == class).collect();

            let mut used = vec![false; ct.len()];
            // hit[i] = did detection i (in confidence order) match a truth?
            let mut hits: Vec<bool> = Vec::with_capacity(cd.len());

            for d in &cd {
                let mut best: Option<(f64, usize)> = None;
                for (ti, t) in ct.iter().enumerate() {
                    if used[ti] {
                        continue;
                    }
                    // Cheap bbox reject before the exact overlay.
                    if let (Some(a), Some(b)) = (d.bbox, t.bbox) {
                        if a.min().x > b.max().x
                            || a.max().x < b.min().x
                            || a.min().y > b.max().y
                            || a.max().y < b.min().y
                        {
                            continue;
                        }
                    }
                    let iou = iou(d, t);
                    if iou >= min_iou && best.is_none_or(|(bi, _)| iou > bi) {
                        best = Some((iou, ti));
                    }
                }
                match best {
                    Some((_, ti)) => {
                        used[ti] = true;
                        hits.push(true);
                    }
                    None => hits.push(false),
                }
            }

            let tp = hits.iter().filter(|h| **h).count();
            let fp = hits.len() - tp;
            let fneg = ct.len() - tp;
            tp_all += tp;
            fp_all += fp;
            fn_all += fneg;

            let precision = if hits.is_empty() {
                0.0
            } else {
                tp as f64 / hits.len() as f64
            };
            let recall = if ct.is_empty() {
                0.0
            } else {
                tp as f64 / ct.len() as f64
            };
            let ap = average_precision(&hits, ct.len());
            // A class with no ground truth has an undefined AP; excluding it
            // keeps mAP from being dragged down by a spurious 0.
            if !ct.is_empty() {
                ap_sum += ap;
                ap_classes += 1;
            }

            out.add_feature(
                None,
                &[
                    ("class", FieldValue::Text(class.clone())),
                    ("truths", FieldValue::Integer(ct.len() as i64)),
                    ("detections", FieldValue::Integer(hits.len() as i64)),
                    ("tp", FieldValue::Integer(tp as i64)),
                    ("fp", FieldValue::Integer(fp as i64)),
                    ("fn", FieldValue::Integer(fneg as i64)),
                    ("precision", FieldValue::Float(precision)),
                    ("recall", FieldValue::Float(recall)),
                    ("ap", FieldValue::Float(ap)),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed writing class row: {e}")))?;
        }

        let map = if ap_classes > 0 {
            ap_sum / ap_classes as f64
        } else {
            0.0
        };

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("map".to_string(), json!(map));
        outputs.insert("class_count".to_string(), json!(classes.len()));
        outputs.insert("true_positives".to_string(), json!(tp_all));
        outputs.insert("false_positives".to_string(), json!(fp_all));
        outputs.insert("false_negatives".to_string(), json!(fn_all));
        outputs.insert(
            "precision".to_string(),
            json!(if tp_all + fp_all > 0 {
                tp_all as f64 / (tp_all + fp_all) as f64
            } else {
                0.0
            }),
        );
        outputs.insert(
            "recall".to_string(),
            json!(if tp_all + fn_all > 0 {
                tp_all as f64 / (tp_all + fn_all) as f64
            } else {
                0.0
            }),
        );
        Ok(ToolRunResult { outputs })
    }
}

struct Inst {
    geom: MultiPolygon,
    area: f64,
    bbox: Option<geo::Rect<f64>>,
    class: String,
    conf: f64,
}

/// Intersection over union of two instances.
fn iou(a: &Inst, b: &Inst) -> f64 {
    let inter = a.geom.intersection(&b.geom).unsigned_area();
    if inter <= 0.0 {
        return 0.0;
    }
    let union = a.area + b.area - inter;
    if union <= 0.0 {
        0.0
    } else {
        inter / union
    }
}

/// All-point interpolated average precision over a confidence-ranked hit list.
///
/// Walks the ranking accumulating the precision-recall curve, then integrates
/// it after making precision monotonically non-increasing from the right (the
/// standard envelope) so a late precision spike cannot inflate the score.
fn average_precision(hits: &[bool], total_truths: usize) -> f64 {
    if total_truths == 0 || hits.is_empty() {
        return 0.0;
    }
    let mut precisions = Vec::with_capacity(hits.len());
    let mut recalls = Vec::with_capacity(hits.len());
    let mut tp = 0usize;
    for (i, h) in hits.iter().enumerate() {
        if *h {
            tp += 1;
        }
        precisions.push(tp as f64 / (i + 1) as f64);
        recalls.push(tp as f64 / total_truths as f64);
    }
    // Precision envelope, right to left.
    for i in (0..precisions.len().saturating_sub(1)).rev() {
        precisions[i] = precisions[i].max(precisions[i + 1]);
    }
    // Integrate: sum precision at each recall increase.
    let mut ap = 0.0;
    let mut prev_recall = 0.0;
    for i in 0..recalls.len() {
        if recalls[i] > prev_recall {
            ap += (recalls[i] - prev_recall) * precisions[i];
            prev_recall = recalls[i];
        }
    }
    ap
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
    match args.get(key).and_then(Value::as_str) {
        Some(s) if !s.trim().is_empty() => Ok(s),
        _ => Err(ToolError::Validation(format!(
            "missing required string parameter '{key}'"
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

fn to_multipolygon(geom: &Geometry) -> Option<MultiPolygon> {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => Some(MultiPolygon(vec![rings_to_polygon(exterior, interiors)])),
        Geometry::MultiPolygon(parts) => Some(MultiPolygon(
            parts.iter().map(|(e, i)| rings_to_polygon(e, i)).collect(),
        )),
        _ => None,
    }
}

fn rings_to_polygon(exterior: &Ring, interiors: &[Ring]) -> Polygon {
    Polygon::new(
        ring_to_linestring(exterior),
        interiors.iter().map(ring_to_linestring).collect(),
    )
}

fn ring_to_linestring(ring: &Ring) -> LineString {
    LineString::new(
        ring.coords()
            .iter()
            .map(|c| GeoCoord { x: c.x, y: c.y })
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
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

    fn boxg(x0: f64, y0: f64, w: f64, h: f64) -> Geometry {
        Geometry::polygon(
            vec![
                Coord::xy(x0, y0),
                Coord::xy(x0 + w, y0),
                Coord::xy(x0 + w, y0 + h),
                Coord::xy(x0, y0 + h),
            ],
            vec![],
        )
    }

    /// Detections carrying class + confidence.
    fn det_layer(items: Vec<(Geometry, &str, f64)>) -> String {
        let mut l = Layer::new("det")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("cls", FieldType::Text));
        l.add_field(FieldDef::new("conf", FieldType::Float));
        for (g, c, s) in items {
            l.add_feature(
                Some(g),
                &[
                    ("cls", FieldValue::Text(c.to_string())),
                    ("conf", FieldValue::Float(s)),
                ],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn gt_layer(items: Vec<(Geometry, &str)>) -> String {
        let mut l = Layer::new("gt")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("cls", FieldType::Text));
        for (g, c) in items {
            l.add_feature(Some(g), &[("cls", FieldValue::Text(c.to_string()))])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ComputeAccuracyForObjectDetectionTool
            .run(&args, &ctx())
            .unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn row(layer: &Layer, class: &str, field: &str) -> f64 {
        let ci = layer.schema.field_index("class").unwrap();
        let fi = layer.schema.field_index(field).unwrap();
        for f in layer.features.iter() {
            if matches!(&f.attributes[ci], FieldValue::Text(s) if s == class) {
                return f.attributes[fi].as_f64().unwrap();
            }
        }
        panic!("class {class} not in table");
    }

    /// A perfect detection set scores 1.0 across the board.
    #[test]
    fn perfect_detection_scores_one() {
        let g = vec![(boxg(0.0, 0.0, 10.0, 10.0), "car")];
        let d = vec![(boxg(0.0, 0.0, 10.0, 10.0), "car", 0.9)];
        let (out, layer) = run(json!({
            "detected": det_layer(d), "ground_truth": gt_layer(g),
            "detected_class_field": "cls", "ground_truth_class_field": "cls",
            "confidence_field": "conf"
        }));
        assert!((out.outputs["map"].as_f64().unwrap() - 1.0).abs() < 1e-9);
        assert_eq!(row(&layer, "car", "tp"), 1.0);
        assert_eq!(row(&layer, "car", "fp"), 0.0);
        assert_eq!(row(&layer, "car", "fn"), 0.0);
    }

    /// A detection below the IoU threshold is a false positive, and the truth
    /// it failed to cover is a false negative.
    #[test]
    fn low_iou_is_false_positive() {
        // Truth [0,10]^2; detection [9,19]x[0,10] -> IoU = 1/19, well under 0.5.
        let g = vec![(boxg(0.0, 0.0, 10.0, 10.0), "car")];
        let d = vec![(boxg(9.0, 0.0, 10.0, 10.0), "car", 0.9)];
        let (_, layer) = run(json!({
            "detected": det_layer(d), "ground_truth": gt_layer(g),
            "detected_class_field": "cls", "ground_truth_class_field": "cls",
            "confidence_field": "conf"
        }));
        assert_eq!(row(&layer, "car", "tp"), 0.0);
        assert_eq!(row(&layer, "car", "fp"), 1.0);
        assert_eq!(row(&layer, "car", "fn"), 1.0);
    }

    /// One truth can only be matched once; the duplicate is a false positive.
    #[test]
    fn duplicate_detection_is_false_positive() {
        let g = vec![(boxg(0.0, 0.0, 10.0, 10.0), "car")];
        let d = vec![
            (boxg(0.0, 0.0, 10.0, 10.0), "car", 0.9),
            (boxg(0.0, 0.0, 10.0, 10.0), "car", 0.8),
        ];
        let (_, layer) = run(json!({
            "detected": det_layer(d), "ground_truth": gt_layer(g),
            "detected_class_field": "cls", "ground_truth_class_field": "cls",
            "confidence_field": "conf"
        }));
        assert_eq!(row(&layer, "car", "tp"), 1.0);
        assert_eq!(row(&layer, "car", "fp"), 1.0);
    }

    /// Confidence ordering matters: a false positive ranked ABOVE a true
    /// positive depresses AP, which a plain precision/recall pair would miss.
    #[test]
    fn confidence_order_affects_ap() {
        let g = gt_layer(vec![(boxg(0.0, 0.0, 10.0, 10.0), "car")]);
        // Good detection, plus a bogus one far away.
        let good = boxg(0.0, 0.0, 10.0, 10.0);
        let bogus = boxg(500.0, 500.0, 10.0, 10.0);

        // Good ranked first -> AP 1.0
        let (hi, _) = run(json!({
            "detected": det_layer(vec![(good.clone(), "car", 0.9), (bogus.clone(), "car", 0.1)]),
            "ground_truth": g, "detected_class_field": "cls",
            "ground_truth_class_field": "cls", "confidence_field": "conf"
        }));
        // Bogus ranked first -> AP 0.5
        let (lo, _) = run(json!({
            "detected": det_layer(vec![(good, "car", 0.1), (bogus, "car", 0.9)]),
            "ground_truth": g, "detected_class_field": "cls",
            "ground_truth_class_field": "cls", "confidence_field": "conf"
        }));
        let (hi_map, lo_map) = (
            hi.outputs["map"].as_f64().unwrap(),
            lo.outputs["map"].as_f64().unwrap(),
        );
        assert!((hi_map - 1.0).abs() < 1e-9, "got {hi_map}");
        assert!((lo_map - 0.5).abs() < 1e-9, "got {lo_map}");
    }

    /// Cross-class matches must not count: a car detection cannot satisfy a
    /// truck truth even at IoU 1.0.
    #[test]
    fn classes_are_scored_independently() {
        let g = vec![(boxg(0.0, 0.0, 10.0, 10.0), "truck")];
        let d = vec![(boxg(0.0, 0.0, 10.0, 10.0), "car", 0.9)];
        let (out, layer) = run(json!({
            "detected": det_layer(d), "ground_truth": gt_layer(g),
            "detected_class_field": "cls", "ground_truth_class_field": "cls",
            "confidence_field": "conf"
        }));
        assert_eq!(row(&layer, "car", "fp"), 1.0);
        assert_eq!(row(&layer, "truck", "fn"), 1.0);
        // mAP averages only over classes that have ground truth.
        assert!(out.outputs["map"].as_f64().unwrap().abs() < 1e-9);
    }

    /// Without class fields everything is one class.
    #[test]
    fn single_class_mode() {
        let g = vec![(boxg(0.0, 0.0, 10.0, 10.0), "ignored")];
        let d = vec![(boxg(0.0, 0.0, 10.0, 10.0), "ignored", 0.5)];
        let (out, _) = run(json!({
            "detected": det_layer(d), "ground_truth": gt_layer(g)
        }));
        assert_eq!(out.outputs["class_count"], json!(1));
        assert!((out.outputs["map"].as_f64().unwrap() - 1.0).abs() < 1e-9);
    }

    /// min_iou is honoured: the same pair flips as the threshold moves.
    #[test]
    fn min_iou_threshold_is_applied() {
        // Overlap 50 of union 150 -> IoU = 1/3.
        let g = gt_layer(vec![(boxg(0.0, 0.0, 10.0, 10.0), "a")]);
        let d = det_layer(vec![(boxg(5.0, 0.0, 10.0, 10.0), "a", 0.9)]);
        let strict = run(json!({
            "detected": d, "ground_truth": g, "min_iou": 0.5
        }));
        assert_eq!(strict.0.outputs["true_positives"], json!(0));
        let loose = run(json!({
            "detected": d, "ground_truth": g, "min_iou": 0.3
        }));
        assert_eq!(loose.0.outputs["true_positives"], json!(1));
    }

    #[test]
    fn rejects_bad_parameters() {
        let g = gt_layer(vec![(boxg(0.0, 0.0, 1.0, 1.0), "a")]);
        let d = det_layer(vec![(boxg(0.0, 0.0, 1.0, 1.0), "a", 1.0)]);
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ComputeAccuracyForObjectDetectionTool
                .validate(&args)
                .is_err()
        };
        assert!(bad(json!({ "ground_truth": g })));
        assert!(bad(json!({ "detected": d })));
        assert!(bad(
            json!({ "detected": d, "ground_truth": g, "min_iou": 1.5 })
        ));
        assert!(bad(
            json!({ "detected": d, "ground_truth": g, "min_iou": -0.1 })
        ));

        // An unknown field is caught at run time.
        let args: ToolArgs = serde_json::from_value(
            json!({ "detected": d, "ground_truth": g, "confidence_field": "nope" }),
        )
        .unwrap();
        assert!(matches!(
            ComputeAccuracyForObjectDetectionTool
                .run(&args, &ctx())
                .unwrap_err(),
            ToolError::Validation(_)
        ));
    }

    /// The AP integrator itself.
    #[test]
    fn average_precision_math() {
        // All hits -> 1.0
        assert!((average_precision(&[true, true], 2) - 1.0).abs() < 1e-9);
        // Miss then hit, 1 truth: precision at the hit is 1/2 -> AP 0.5
        assert!((average_precision(&[false, true], 1) - 0.5).abs() < 1e-9);
        // No truths -> 0, not NaN
        assert_eq!(average_precision(&[true], 0), 0.0);
    }
}
