//! GeoLibre tool: non-maximum suppression of overlapping detections.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Non Maximum Suppression* (Image
//! Analyst). Object-detection and OBIA outputs routinely emit many overlapping
//! polygons for one real object; NMS keeps, in each overlap cluster, only the
//! highest-`confidence_score_field` detection and discards every lower-scored
//! detection whose intersection-over-union (IoU) with an already-kept one
//! exceeds `max_overlap_ratio`.
//!
//! The bundled/authored suite has no score-based suppressor:
//! * `count_overlapping_features` *counts* overlaps, keeping all of them;
//! * `resolve_building_conflicts` *displaces* features, never removes duplicates;
//! * `find_identical` dedups only *exact* matches, not overlap-by-IoU.
//!
//! Features are sorted by descending score; a greedy pass keeps a detection
//! then rejects later detections that overlap it too much (optionally only
//! within the same `class_value_field`). Bounding-box prefiltering skips pairs
//! that cannot overlap. IoU uses the same `geo` `BooleanOps` machinery as
//! `count_overlapping_features`. Surviving features keep all input attributes.

use std::collections::BTreeMap;

use geo::{Area, BooleanOps, BoundingRect, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldValue, Geometry, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct NonMaximumSuppressionTool;

impl Tool for NonMaximumSuppressionTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "non_maximum_suppression",
            display_name: "Non Maximum Suppression",
            summary: "Keep the highest-confidence detection in each overlap cluster and discard lower-scored polygons whose intersection-over-union exceeds a threshold (like ArcGIS Non Maximum Suppression) — the score-based deduplicator for object-detection/OBIA output that count_overlapping_features (counts, keeps all), resolve_building_conflicts (displaces), and find_identical (exact only) don't provide.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon detections (may overlap).",
                    required: true,
                },
                ToolParamSpec {
                    name: "confidence_score_field",
                    description: "Numeric field holding each detection's confidence/score; higher wins.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon path with the surviving detections. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_overlap_ratio",
                    description: "IoU above which a lower-scored detection is suppressed (0-1, default 0.5).",
                    required: false,
                },
                ToolParamSpec {
                    name: "class_value_field",
                    description: "Optional field; only detections of the same class suppress each other.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args
            .get("input")
            .and_then(Value::as_str)
            .map(str::trim)
            .unwrap_or("")
            .is_empty()
        {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        if args
            .get("confidence_score_field")
            .and_then(Value::as_str)
            .map(str::trim)
            .unwrap_or("")
            .is_empty()
        {
            return Err(ToolError::Validation(
                "missing required string parameter 'confidence_score_field'".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args
            .get("input")
            .and_then(Value::as_str)
            .filter(|s| !s.trim().is_empty())
            .ok_or_else(|| {
                ToolError::Validation("missing required parameter 'input'".to_string())
            })?;
        let score_field = args
            .get("confidence_score_field")
            .and_then(Value::as_str)
            .unwrap();
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let score_idx = layer.schema.field_index(score_field).ok_or_else(|| {
            ToolError::Validation(format!("confidence_score_field '{score_field}' not found"))
        })?;
        let class_idx = match &prm.class_value_field {
            Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                ToolError::Validation(format!("class_value_field '{f}' not found"))
            })?),
            None => None,
        };

        // Collect polygon detections with score, class, geometry, bbox.
        let mut dets: Vec<Detection> = Vec::new();
        for (fi, feature) in layer.features.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let Some(mp) = to_multipolygon(geom) else {
                continue;
            };
            let area = mp.unsigned_area();
            if area <= 0.0 {
                continue;
            }
            let score = feature
                .attributes
                .get(score_idx)
                .and_then(|v| v.as_f64())
                .ok_or_else(|| {
                    ToolError::Execution(format!(
                        "feature {fi} has a missing/non-numeric '{score_field}'"
                    ))
                })?;
            let class = class_idx.map(|ci| field_key(&feature.attributes[ci]));
            let bbox = bbox(&mp);
            dets.push(Detection {
                feat: fi,
                mp,
                area,
                bbox,
                score,
                class,
            });
        }
        let n_in = dets.len();
        if n_in == 0 {
            return Err(ToolError::Execution(
                "no polygon features in input".to_string(),
            ));
        }

        // Greedy NMS: process by descending score, keep a detection, then mark
        // any not-yet-decided detection overlapping it (IoU > threshold) as
        // suppressed. Same-class-only when class_value_field is set.
        ctx.progress
            .info(&format!("suppressing {n_in} detection(s)"));
        let mut order: Vec<usize> = (0..n_in).collect();
        order.sort_by(|&a, &b| dets[b].score.total_cmp(&dets[a].score));

        let mut suppressed = vec![false; n_in];
        let mut kept: Vec<usize> = Vec::new();
        for &i in &order {
            if suppressed[i] {
                continue;
            }
            kept.push(i);
            for &j in &order {
                if j == i || suppressed[j] {
                    continue;
                }
                if dets[i].class != dets[j].class {
                    continue;
                }
                if !bbox_overlap(&dets[i].bbox, &dets[j].bbox) {
                    continue;
                }
                if iou(&dets[i], &dets[j]) > prm.max_overlap_ratio {
                    suppressed[j] = true;
                }
            }
        }

        // Emit survivors in original feature order, preserving all attributes.
        let mut keep_feat = vec![false; layer.len()];
        for &i in &kept {
            keep_feat[dets[i].feat] = true;
        }
        let mut out = Layer::new("nms");
        out.schema = layer.schema.clone();
        for (fi, feature) in layer.features.iter().enumerate() {
            if keep_feat[fi] {
                out.push(feature.clone());
            }
        }

        let n_kept = kept.len();
        let n_suppressed = n_in - n_kept;
        ctx.progress
            .info(&format!("{n_kept} kept, {n_suppressed} suppressed"));

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_count".to_string(), json!(n_in));
        outputs.insert("kept_count".to_string(), json!(n_kept));
        outputs.insert("suppressed_count".to_string(), json!(n_suppressed));
        Ok(ToolRunResult { outputs })
    }
}

// ── Detections / IoU ──────────────────────────────────────────────────────────

struct Detection {
    feat: usize,
    mp: MultiPolygon,
    area: f64,
    bbox: [f64; 4],
    score: f64,
    class: Option<String>,
}

/// Intersection-over-union of two detections.
fn iou(a: &Detection, b: &Detection) -> f64 {
    let inter = a.mp.intersection(&b.mp).unsigned_area();
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

fn bbox(mp: &MultiPolygon) -> [f64; 4] {
    match mp.bounding_rect() {
        Some(r) => [r.min().x, r.min().y, r.max().x, r.max().y],
        None => [
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        ],
    }
}

fn bbox_overlap(a: &[f64; 4], b: &[f64; 4]) -> bool {
    a[0] <= b[2] && b[0] <= a[2] && a[1] <= b[3] && b[1] <= a[3]
}

fn field_key(fv: &FieldValue) -> String {
    if let Some(i) = fv.as_i64() {
        i.to_string()
    } else if let Some(f) = fv.as_f64() {
        format!("{f}")
    } else {
        fv.as_str().unwrap_or("").to_string()
    }
}

// ── geo <-> wbvector conversion (shared shape with count_overlapping_features) ─

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

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    max_overlap_ratio: f64,
    class_value_field: Option<String>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let max_overlap_ratio = match args.get("max_overlap_ratio") {
        None | Some(Value::Null) => 0.5,
        Some(Value::Number(n)) => n.as_f64().unwrap_or(0.5),
        Some(Value::String(s)) if s.trim().is_empty() => 0.5,
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'max_overlap_ratio' must be a number".into()))?,
        Some(_) => {
            return Err(ToolError::Validation(
                "'max_overlap_ratio' must be a number".into(),
            ))
        }
    };
    if !(0.0..=1.0).contains(&max_overlap_ratio) {
        return Err(ToolError::Validation(
            "'max_overlap_ratio' must be in [0, 1]".into(),
        ));
    }
    let class_value_field = parse_optional_str(args, "class_value_field")?.map(str::to_string);
    Ok(Params {
        max_overlap_ratio,
        class_value_field,
    })
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

    fn square(x0: f64, y0: f64, s: f64) -> Geometry {
        Geometry::polygon(
            vec![
                Coord::xy(x0, y0),
                Coord::xy(x0 + s, y0),
                Coord::xy(x0 + s, y0 + s),
                Coord::xy(x0, y0 + s),
            ],
            vec![],
        )
    }

    /// Each detection is (geometry, score, class).
    fn layer_of(dets: Vec<(Geometry, f64, &str)>) -> String {
        let mut l = Layer::new("dets")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(wbvector::FieldDef::new("score", wbvector::FieldType::Float));
        l.add_field(wbvector::FieldDef::new("class", wbvector::FieldType::Text));
        for (g, s, c) in dets {
            l.add_feature(Some(g), &[("score", s.into()), ("class", c.into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = NonMaximumSuppressionTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// Two near-identical boxes (IoU ~1): the higher score survives.
    #[test]
    fn suppresses_duplicate_keeps_highest() {
        let input = layer_of(vec![
            (square(0.0, 0.0, 10.0), 0.9, "car"),
            (square(0.2, 0.2, 10.0), 0.6, "car"), // ~96% IoU with the first
        ]);
        let (out, layer) = run(json!({ "input": input, "confidence_score_field": "score" }));
        assert_eq!(out.outputs["kept_count"], json!(1));
        let sidx = layer.schema.field_index("score").unwrap();
        assert!((layer.features[0].attributes[sidx].as_f64().unwrap() - 0.9).abs() < 1e-9);
    }

    /// Two well-separated boxes (IoU 0): both survive.
    #[test]
    fn keeps_disjoint_detections() {
        let input = layer_of(vec![
            (square(0.0, 0.0, 10.0), 0.9, "car"),
            (square(100.0, 100.0, 10.0), 0.8, "car"),
        ]);
        let (out, _l) = run(json!({ "input": input, "confidence_score_field": "score" }));
        assert_eq!(out.outputs["kept_count"], json!(2));
    }

    /// Partial overlap: suppressed only when IoU exceeds the threshold.
    #[test]
    fn threshold_controls_suppression() {
        // Two 10x10 boxes offset by 5 in x: intersection 50, union 150, IoU=1/3.
        let mk = || {
            layer_of(vec![
                (square(0.0, 0.0, 10.0), 0.9, "car"),
                (square(5.0, 0.0, 10.0), 0.6, "car"),
            ])
        };
        // threshold 0.5 > 0.333 -> both kept.
        let (hi, _l) = run(
            json!({ "input": mk(), "confidence_score_field": "score", "max_overlap_ratio": 0.5 }),
        );
        assert_eq!(hi.outputs["kept_count"], json!(2));
        // threshold 0.2 < 0.333 -> lower one suppressed.
        let (lo, _l) = run(
            json!({ "input": mk(), "confidence_score_field": "score", "max_overlap_ratio": 0.2 }),
        );
        assert_eq!(lo.outputs["kept_count"], json!(1));
    }

    /// class_value_field: overlapping detections of different classes coexist.
    #[test]
    fn class_scoped_suppression() {
        let input = layer_of(vec![
            (square(0.0, 0.0, 10.0), 0.9, "car"),
            (square(0.2, 0.2, 10.0), 0.6, "truck"), // overlaps but different class
        ]);
        let (out, _l) = run(json!({
            "input": input, "confidence_score_field": "score", "class_value_field": "class"
        }));
        assert_eq!(
            out.outputs["kept_count"],
            json!(2),
            "different classes are not suppressed"
        );
    }

    /// A chain of three heavily-overlapping boxes collapses to the top one.
    #[test]
    fn chain_collapses_to_best() {
        let input = layer_of(vec![
            (square(0.0, 0.0, 10.0), 0.5, "car"),
            (square(0.3, 0.3, 10.0), 0.95, "car"),
            (square(0.6, 0.6, 10.0), 0.7, "car"),
        ]);
        let (out, layer) = run(json!({ "input": input, "confidence_score_field": "score" }));
        assert_eq!(out.outputs["kept_count"], json!(1));
        let sidx = layer.schema.field_index("score").unwrap();
        assert!((layer.features[0].attributes[sidx].as_f64().unwrap() - 0.95).abs() < 1e-9);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            NonMaximumSuppressionTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err()); // no score field
        assert!(bad(
            json!({ "input": "a.geojson", "confidence_score_field": "s", "max_overlap_ratio": 1.5 })
        )
        .is_err());
        assert!(bad(json!({ "input": "a.geojson", "confidence_score_field": "s" })).is_ok());
    }
}
