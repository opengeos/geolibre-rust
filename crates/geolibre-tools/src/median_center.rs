//! GeoLibre tool: median center (geometric / Weiszfeld median).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Median Center* (Spatial Statistics →
//! Measuring Geographic Distributions). For each group of point features it
//! emits a single point at the **geometric median** — the location that
//! minimizes the sum of (optionally weighted) Euclidean distances to every
//! input point. Unlike the *mean center* (which minimizes summed *squared*
//! distance and is pulled by outliers) the geometric median is robust; unlike
//! the *central feature* the result need not coincide with an input point.
//!
//! The median is found by **Weiszfeld's algorithm**: starting from the weighted
//! mean center, each iteration replaces the estimate with the weighted average
//! of the points, each point inverse-weighted by its distance to the current
//! estimate, until the update falls below a tolerance. A tiny distance floor
//! keeps the update well-defined when the estimate lands on a data point.
//!
//! `weight_field` scales each point's pull; `case_field` partitions the input
//! into groups (one median per distinct value); `attribute_fields` lists numeric
//! fields whose per-group **median value** is appended to the output.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct MedianCenterTool;

impl Tool for MedianCenterTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "median_center",
            display_name: "Median Center",
            summary: "The geometric median (Weiszfeld point that minimizes summed, optionally weighted, Euclidean distance to all input points) — one robust center point per case group, with optional per-field median attributes (ArcGIS Median Center).",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "weight_field",
                    description: "Optional numeric field weighting each point's pull on the median.",
                    required: false,
                },
                ToolParamSpec {
                    name: "case_field",
                    description: "Optional field to group points by; one median per distinct value.",
                    required: false,
                },
                ToolParamSpec {
                    name: "attribute_fields",
                    description: "Optional comma-separated numeric fields; the per-group median of each is appended as '<field>_med'.",
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
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;

        let case_idx = match &prm.case_field {
            Some(f) => Some(
                layer
                    .schema
                    .field_index(f)
                    .ok_or_else(|| ToolError::Validation(format!("case_field '{f}' not found")))?,
            ),
            None => None,
        };
        let weight_idx =
            match &prm.weight_field {
                Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                    ToolError::Validation(format!("weight_field '{f}' not found"))
                })?),
                None => None,
            };
        let attr_idx: Vec<(String, usize)> = prm
            .attribute_fields
            .iter()
            .map(|f| {
                layer
                    .schema
                    .field_index(f)
                    .map(|i| (f.clone(), i))
                    .ok_or_else(|| {
                        ToolError::Validation(format!("attribute field '{f}' not found"))
                    })
            })
            .collect::<Result<_, _>>()?;

        // Group point representatives (x, y, weight) plus row index by case value.
        let mut groups: BTreeMap<String, Vec<(f64, f64, f64, usize)>> = BTreeMap::new();
        for (fidx, feature) in layer.features.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let Some((x, y)) = representative_xy(geom) else {
                continue;
            };
            let w = weight_idx
                .and_then(|wi| feature.attributes.get(wi))
                .and_then(FieldValue::as_f64)
                .filter(|v| v.is_finite() && *v > 0.0)
                .unwrap_or(1.0);
            let key = match case_idx {
                Some(i) => value_string(&feature.attributes[i]),
                None => "ALL".to_string(),
            };
            groups.entry(key).or_default().push((x, y, w, fidx));
        }

        // Build output layer.
        let mut out = Layer::new("median_center").with_geom_type(GeometryType::Point);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("case", FieldType::Text));
        out.add_field(FieldDef::new("n", FieldType::Integer));
        out.add_field(FieldDef::new("sum_dist", FieldType::Float));
        for (name, _) in &attr_idx {
            out.add_field(FieldDef::new(format!("{name}_med"), FieldType::Float));
        }

        let mut count = 0usize;
        for (key, pts) in &groups {
            if pts.is_empty() {
                continue;
            }
            let (mx, my) = weighted_mean(pts);
            let (cx, cy) = weiszfeld(pts, mx, my);
            let sum_dist: f64 = pts
                .iter()
                .map(|&(x, y, w, _)| w * (x - cx).hypot(y - cy))
                .sum();

            let mut attrs = vec![
                FieldValue::Text(key.clone()),
                FieldValue::Integer(pts.len() as i64),
                FieldValue::Float(sum_dist),
            ];
            for (_, idx) in &attr_idx {
                let mut vals: Vec<f64> = pts
                    .iter()
                    .filter_map(|&(_, _, _, fi)| layer.features[fi].attributes.get(*idx))
                    .filter_map(FieldValue::as_f64)
                    .filter(|v| v.is_finite())
                    .collect();
                attrs.push(FieldValue::Float(median(&mut vals)));
            }

            out.push(Feature {
                fid: 0,
                geometry: Some(Geometry::point(cx, cy)),
                attributes: attrs,
            });
            count += 1;
        }

        ctx.progress.info(&format!(
            "{} group(s) -> {count} median center(s)",
            groups.len()
        ));

        let path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(path));
        outputs.insert("group_count".to_string(), json!(groups.len()));
        outputs.insert("feature_count".to_string(), json!(count));
        Ok(ToolRunResult { outputs })
    }
}

// ── Geometric median ──────────────────────────────────────────────────────────

/// Weighted mean center of the points (the Weiszfeld seed).
fn weighted_mean(pts: &[(f64, f64, f64, usize)]) -> (f64, f64) {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sw = 0.0;
    for &(x, y, w, _) in pts {
        sx += w * x;
        sy += w * y;
        sw += w;
    }
    if sw <= 0.0 {
        return (0.0, 0.0);
    }
    (sx / sw, sy / sw)
}

/// Weiszfeld iteration: the weighted geometric median seeded at (sx, sy).
///
/// Each step moves the estimate to the weighted average of the points, each
/// point inverse-weighted by its distance to the current estimate. A tiny
/// distance floor keeps the update finite when the estimate coincides with a
/// point (there the point's inverse weight dominates and pins the median, which
/// is exactly the geometric median when a point carries the majority weight).
fn weiszfeld(pts: &[(f64, f64, f64, usize)], sx: f64, sy: f64) -> (f64, f64) {
    const EPS: f64 = 1e-12;
    const MAX_ITERS: usize = 10_000;
    let mut cx = sx;
    let mut cy = sy;
    for _ in 0..MAX_ITERS {
        let mut num_x = 0.0;
        let mut num_y = 0.0;
        let mut den = 0.0;
        for &(x, y, w, _) in pts {
            let d = ((x - cx).hypot(y - cy)).max(EPS);
            let inv = w / d;
            num_x += inv * x;
            num_y += inv * y;
            den += inv;
        }
        if den <= 0.0 {
            break;
        }
        let nx = num_x / den;
        let ny = num_y / den;
        let step = (nx - cx).hypot(ny - cy);
        cx = nx;
        cy = ny;
        if step < 1e-12 {
            break;
        }
    }
    (cx, cy)
}

fn median(vals: &mut [f64]) -> f64 {
    if vals.is_empty() {
        return f64::NAN;
    }
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = vals.len();
    if n % 2 == 1 {
        vals[n / 2]
    } else {
        0.5 * (vals[n / 2 - 1] + vals[n / 2])
    }
}

// ── Geometry helpers ─────────────────────────────────────────────────────────

/// Representative (x, y) for a geometry — its coordinate mean. For points this
/// is the point itself; for other geometries it degrades to the vertex centroid.
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

fn value_string(fv: &FieldValue) -> String {
    if let Some(i) = fv.as_i64() {
        i.to_string()
    } else if let Some(f) = fv.as_f64() {
        format!("{f}")
    } else {
        fv.as_str().unwrap_or("").to_string()
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    weight_field: Option<String>,
    case_field: Option<String>,
    attribute_fields: Vec<String>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let weight_field = parse_optional_str(args, "weight_field")?.map(str::to_string);
    let case_field = parse_optional_str(args, "case_field")?.map(str::to_string);
    let attribute_fields = parse_optional_str(args, "attribute_fields")?
        .map(|s| {
            s.split(',')
                .map(str::trim)
                .filter(|t| !t.is_empty())
                .map(str::to_string)
                .collect::<Vec<_>>()
        })
        .unwrap_or_default();
    Ok(Params {
        weight_field,
        case_field,
        attribute_fields,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::memory_store;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn point_layer(pts: &[(f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        for (i, (x, y)) in pts.iter().enumerate() {
            l.add_feature(
                Some(Geometry::point(*x, *y)),
                &[("name", format!("p{i}").into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = MedianCenterTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn center(layer: &Layer) -> (f64, f64) {
        match layer.features[0].geometry.as_ref().unwrap() {
            Geometry::Point(c) => (c.x, c.y),
            _ => panic!("expected point"),
        }
    }

    /// Summed distance from the emitted median must be <= that from the mean
    /// center — the defining property of the geometric median.
    #[test]
    fn beats_the_mean_center() {
        // Cluster of four with a lone far outlier that drags the mean away.
        let pts = [
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 1.0),
            (1.0, 1.0),
            (100.0, 100.0),
        ];
        let input = point_layer(&pts);
        let (_o, layer) = run(json!({ "input": input }));
        let (cx, cy) = center(&layer);

        let mean_x = pts.iter().map(|p| p.0).sum::<f64>() / pts.len() as f64;
        let mean_y = pts.iter().map(|p| p.1).sum::<f64>() / pts.len() as f64;
        let sum_from = |x: f64, y: f64| {
            pts.iter()
                .map(|&(px, py)| (px - x).hypot(py - y))
                .sum::<f64>()
        };
        let d_med = sum_from(cx, cy);
        let d_mean = sum_from(mean_x, mean_y);
        assert!(
            d_med <= d_mean + 1e-6,
            "median sum {d_med} should not exceed mean sum {d_mean}"
        );
        // Robust: the median stays inside the cluster, unlike the mean (~20,20).
        assert!(
            cx < 10.0 && cy < 10.0,
            "median dragged by outlier: {cx},{cy}"
        );

        // Reported sum_dist matches.
        let sidx = layer.schema.field_index("sum_dist").unwrap();
        let reported = layer.features[0].attributes[sidx].as_f64().unwrap();
        assert!((reported - d_med).abs() < 1e-6);
    }

    /// For symmetric points the median sits at the centre of symmetry.
    #[test]
    fn symmetric_points_center_exactly() {
        let input = point_layer(&[(-2.0, 0.0), (2.0, 0.0), (0.0, -2.0), (0.0, 2.0)]);
        let (_o, layer) = run(json!({ "input": input }));
        let (cx, cy) = center(&layer);
        assert!(
            cx.abs() < 1e-6 && cy.abs() < 1e-6,
            "expected origin: {cx},{cy}"
        );
    }

    /// A heavy weight pulls the median toward that point.
    #[test]
    fn weight_pulls_the_median() {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("w", FieldType::Float));
        for (x, w) in [(0.0, 1.0), (1.0, 1.0), (10.0, 1.0)] {
            l.add_feature(Some(Geometry::point(x, 0.0)), &[("w", w.into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);

        // Unweighted median of {0,1,10} on a line is the middle point x=1.
        let (_o, unw) = run(json!({ "input": input.clone() }));
        assert!((center(&unw).0 - 1.0).abs() < 1e-6);

        // Heavily weight x=10 -> median jumps to x=10.
        let mut l2 = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l2.add_field(FieldDef::new("w", FieldType::Float));
        for (x, w) in [(0.0, 1.0), (1.0, 1.0), (10.0, 100.0)] {
            l2.add_feature(Some(Geometry::point(x, 0.0)), &[("w", w.into())])
                .unwrap();
        }
        let id2 = memory_store::put_vector(l2);
        let input2 = memory_store::make_vector_memory_path(&id2);
        let (_o2, w) = run(json!({ "input": input2, "weight_field": "w" }));
        assert!(center(&w).0 > 9.0, "weighted median should approach 10");
    }

    /// case_field yields one median per group; attribute_fields appends a median.
    #[test]
    fn grouping_and_attribute_median() {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("grp", FieldType::Text));
        l.add_field(FieldDef::new("v", FieldType::Float));
        let rows = [
            (0.0, 0.0, "a", 10.0),
            (2.0, 0.0, "a", 20.0),
            (4.0, 0.0, "a", 30.0),
            (0.0, 50.0, "b", 5.0),
            (0.0, 52.0, "b", 7.0),
        ];
        for (x, y, g, v) in rows {
            l.add_feature(
                Some(Geometry::point(x, y)),
                &[("grp", g.into()), ("v", v.into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);
        let (out, layer) = run(json!({
            "input": input, "case_field": "grp", "attribute_fields": "v",
        }));
        assert_eq!(out.outputs["feature_count"], json!(2));
        // Group "a" median value is 20 (median of 10,20,30).
        let midx = layer.schema.field_index("v_med").unwrap();
        let cidx = layer.schema.field_index("case").unwrap();
        let a = layer
            .features
            .iter()
            .find(|f| f.attributes[cidx].as_str() == Some("a"))
            .unwrap();
        assert!((a.attributes[midx].as_f64().unwrap() - 20.0).abs() < 1e-9);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            MedianCenterTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "" })).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_ok());
    }

    /// A missing case/weight field is reported at run time.
    #[test]
    fn missing_field_errors() {
        let input = point_layer(&[(0.0, 0.0), (1.0, 1.0)]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "weight_field": "nope" })).unwrap();
        assert!(MedianCenterTool.run(&args, &ctx()).is_err());
    }
}
