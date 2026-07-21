//! GeoLibre tool: area-weighted attribute transfer between polygon layers.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Apportion Polygon* (Analysis; related
//! to Areal Interpolation). `tabulate_intersection` (#69) reports overlap
//! areas/percentages but stops short of transferring attribute *values* onto the
//! target features. This is the one-step-further tool every census-to-service-
//! area workflow needs — dasymetric reaggregation between incompatible zone
//! systems — and it reuses the same `geo` `BooleanOps` overlay.
//!
//! Each source polygon's numeric field value is split among the target polygons
//! it overlaps, in proportion to the overlap weight: the intersection area
//! (`method = area`) or the intersection area times a target field
//! (`method = weight`). Because each source value is normalised by its total
//! overlap with the targets, **the full value is distributed among the
//! intersecting targets** — so when the targets tile the source layer, the
//! summed apportioned total equals the original total exactly. The output is a
//! copy of the target layer with one apportioned column per field
//! (`<field><suffix>`, default suffix `_app`).

use std::collections::BTreeMap;

use geo::{Area, BooleanOps, BoundingRect, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct ApportionPolygonTool;

impl Tool for ApportionPolygonTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "apportion_polygon",
            display_name: "Apportion Polygon",
            summary: "Transfer numeric attributes from a source polygon layer onto a target polygon layer, apportioned by area of overlap (dasymetric reaggregation between zone systems), like ArcGIS Apportion Polygon.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "target",
                    description: "Target polygon layer (receives the apportioned values; its features and attributes are preserved).",
                    required: true,
                },
                ToolParamSpec {
                    name: "source",
                    description: "Source polygon layer carrying the values to distribute.",
                    required: true,
                },
                ToolParamSpec {
                    name: "fields",
                    description: "Comma-separated numeric field name(s) on the source to apportion.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (a copy of the target with apportioned columns). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'area' (default; split by overlap area) or 'weight' (split by overlap area times the target 'weight_field').",
                    required: false,
                },
                ToolParamSpec {
                    name: "weight_field",
                    description: "Target numeric field used as an extra weight when method=weight.",
                    required: false,
                },
                ToolParamSpec {
                    name: "suffix",
                    description: "Suffix for the apportioned output columns. Default '_app'.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "target")?;
        require_str(args, "source")?;
        require_str(args, "fields")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let target_path = require_str(args, "target")?;
        let source_path = require_str(args, "source")?;
        let fields_arg = require_str(args, "fields")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let mut target = load_input_layer(target_path)?;
        let source = load_input_layer(source_path)?;

        let field_names: Vec<String> = fields_arg
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect();
        if field_names.is_empty() {
            return Err(ToolError::Validation(
                "'fields' must name at least one field".to_string(),
            ));
        }
        let src_field_idx: Vec<usize> = field_names
            .iter()
            .map(|f| {
                source
                    .schema
                    .field_index(f)
                    .ok_or_else(|| ToolError::Validation(format!("source field '{f}' not found")))
            })
            .collect::<Result<_, _>>()?;

        let weight_idx = match (prm.method, &prm.weight_field) {
            (Method::Weight, Some(wf)) => Some(target.schema.field_index(wf).ok_or_else(|| {
                ToolError::Validation(format!("weight_field '{wf}' not found on target"))
            })?),
            (Method::Weight, None) => None, // guarded in parse_params
            (Method::Area, _) => None,
        };

        // Target polygons as geo geometries with bounding boxes.
        let targets: Vec<Option<(MultiPolygon, [f64; 4], f64)>> = target
            .features
            .iter()
            .map(|f| {
                f.geometry.as_ref().and_then(to_multipolygon).map(|mp| {
                    let bb = bbox(&mp);
                    let w = weight_idx
                        .map(|wi| {
                            f.attributes
                                .get(wi)
                                .and_then(FieldValue::as_f64)
                                .unwrap_or(0.0)
                        })
                        .unwrap_or(1.0);
                    (mp, bb, w)
                })
            })
            .collect();

        // Accumulator: one running total per target per field.
        let mut acc: Vec<Vec<f64>> = vec![vec![0.0; field_names.len()]; target.len()];

        ctx.progress.info(&format!(
            "apportioning {} field(s) from {} source polygon(s) onto {} target(s)",
            field_names.len(),
            source.len(),
            target.len()
        ));

        for src in source.features.iter() {
            let Some(src_mp) = src.geometry.as_ref().and_then(to_multipolygon) else {
                continue;
            };
            let src_bb = bbox(&src_mp);
            // Source field values.
            let vals: Vec<f64> = src_field_idx
                .iter()
                .map(|&i| {
                    src.attributes
                        .get(i)
                        .and_then(FieldValue::as_f64)
                        .unwrap_or(0.0)
                })
                .collect();

            // Overlap weight with each intersecting target.
            let mut overlaps: Vec<(usize, f64)> = Vec::new();
            let mut total_w = 0.0;
            for (ti, t) in targets.iter().enumerate() {
                let Some((t_mp, t_bb, t_weight)) = t else {
                    continue;
                };
                if !bbox_overlap(&src_bb, t_bb) {
                    continue;
                }
                let inter_area = src_mp.intersection(t_mp).unsigned_area();
                if inter_area <= 0.0 {
                    continue;
                }
                let w = inter_area * *t_weight;
                if w > 0.0 {
                    overlaps.push((ti, w));
                    total_w += w;
                }
            }
            if total_w <= 0.0 {
                continue;
            }
            // Distribute each field value by the normalised overlap weight.
            for (ti, w) in &overlaps {
                let frac = w / total_w;
                for (fi, v) in vals.iter().enumerate() {
                    acc[*ti][fi] += v * frac;
                }
            }
        }

        // Append apportioned columns to the target.
        for name in &field_names {
            target.add_field(FieldDef::new(
                format!("{name}{}", prm.suffix),
                FieldType::Float,
            ));
        }
        for (ti, feature) in target.features.iter_mut().enumerate() {
            for v in &acc[ti] {
                feature.attributes.push(FieldValue::Float(*v));
            }
        }

        let out_totals: Vec<f64> = (0..field_names.len())
            .map(|fi| acc.iter().map(|row| row[fi]).sum())
            .collect();
        ctx.progress.info(&format!(
            "apportioned totals: {:?}",
            field_names
                .iter()
                .zip(&out_totals)
                .map(|(n, t)| format!("{n}={t:.3}"))
                .collect::<Vec<_>>()
        ));

        let feature_count = target.len();
        let out_path = write_or_store_layer(target, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        for (name, total) in field_names.iter().zip(&out_totals) {
            outputs.insert(format!("total_{name}"), json!(total));
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── geo <-> wbvector conversion ──────────────────────────────────────────────

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

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Area,
    Weight,
}

struct Params {
    method: Method,
    weight_field: Option<String>,
    suffix: String,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let method = match parse_optional_str(args, "method")? {
        None => Method::Area,
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "area" => Method::Area,
            "weight" | "weight_field" => Method::Weight,
            other => {
                return Err(ToolError::Validation(format!(
                    "'method' must be 'area' or 'weight', got '{other}'"
                )))
            }
        },
    };
    let weight_field = parse_optional_str(args, "weight_field")?.map(str::to_string);
    if method == Method::Weight && weight_field.is_none() {
        return Err(ToolError::Validation(
            "method=weight requires 'weight_field'".to_string(),
        ));
    }
    let suffix = parse_optional_str(args, "suffix")?
        .map(str::to_string)
        .unwrap_or_else(|| "_app".to_string());
    Ok(Params {
        method,
        weight_field,
        suffix,
    })
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
    use wbvector::{memory_store, Coord, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn rect(x0: f64, y0: f64, w: f64, h: f64) -> Geometry {
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

    fn poly_layer(feats: &[(Geometry, f64)]) -> String {
        let mut l = Layer::new("polys")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("pop", FieldType::Float));
        for (g, v) in feats {
            l.add_feature(Some(g.clone()), &[("pop", (*v).into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ApportionPolygonTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn apps(layer: &Layer) -> Vec<f64> {
        let idx = layer.schema.field_index("pop_app").unwrap();
        layer
            .iter()
            .map(|f| f.attributes[idx].as_f64().unwrap())
            .collect()
    }

    /// One source (pop 100) fully covered by two equal target halves -> 50/50.
    #[test]
    fn splits_value_by_area() {
        // Source: [0,10]x[0,10] pop 100. Targets: left [0,5], right [5,10].
        let source = poly_layer(&[(rect(0.0, 0.0, 10.0, 10.0), 100.0)]);
        let target = poly_layer(&[
            (rect(0.0, 0.0, 5.0, 10.0), 0.0),
            (rect(5.0, 0.0, 5.0, 10.0), 0.0),
        ]);
        let (out, layer) = run(json!({ "target": target, "source": source, "fields": "pop" }));
        let a = apps(&layer);
        assert!(
            (a[0] - 50.0).abs() < 1e-6 && (a[1] - 50.0).abs() < 1e-6,
            "expected 50/50, got {a:?}"
        );
        // Conservation: total apportioned == source total.
        assert!((out.outputs["total_pop"].as_f64().unwrap() - 100.0).abs() < 1e-6);
    }

    /// Uneven target areas get proportional shares (3:1).
    #[test]
    fn proportional_to_overlap_area() {
        let source = poly_layer(&[(rect(0.0, 0.0, 8.0, 10.0), 80.0)]);
        // Targets [0,6] (area 60) and [6,8] (area 20) -> 3:1 -> 60 and 20.
        let target = poly_layer(&[
            (rect(0.0, 0.0, 6.0, 10.0), 0.0),
            (rect(6.0, 0.0, 2.0, 10.0), 0.0),
        ]);
        let (_o, layer) = run(json!({ "target": target, "source": source, "fields": "pop" }));
        let a = apps(&layer);
        assert!(
            (a[0] - 60.0).abs() < 1e-6 && (a[1] - 20.0).abs() < 1e-6,
            "expected 60/20, got {a:?}"
        );
    }

    /// Full value is distributed among intersecting targets even when the source
    /// extends beyond them (normalisation).
    #[test]
    fn normalizes_over_covered_area() {
        // Source [0,10]x[0,10] pop 100; a single target covering only the left
        // half. The full 100 still lands on that target (it is the only overlap).
        let source = poly_layer(&[(rect(0.0, 0.0, 10.0, 10.0), 100.0)]);
        let target = poly_layer(&[(rect(0.0, 0.0, 5.0, 10.0), 0.0)]);
        let (out, layer) = run(json!({ "target": target, "source": source, "fields": "pop" }));
        assert!((apps(&layer)[0] - 100.0).abs() < 1e-6);
        assert!((out.outputs["total_pop"].as_f64().unwrap() - 100.0).abs() < 1e-6);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ApportionPolygonTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "target": "t.geojson", "source": "s.geojson" })).is_err());
        assert!(bad(json!({ "target": "t.geojson", "source": "s.geojson", "fields": "pop", "method": "weight" })).is_err());
        assert!(
            bad(json!({ "target": "t.geojson", "source": "s.geojson", "fields": "pop" })).is_ok()
        );
    }
}
