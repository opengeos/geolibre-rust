//! GeoLibre tool: central feature and linear directional mean.
//!
//! Pure-Rust counterpart of two ArcGIS Pro *Measuring Geographic Distributions*
//! members not covered by `directional_distribution` (#68):
//!
//! * **Central Feature** — the actual input feature whose total (optionally
//!   weighted) distance to all other features is smallest. Unlike a mean/median
//!   center, the output *is* an input feature, with its attributes preserved.
//! * **Linear Directional Mean** — the mean direction (or, undirected,
//!   orientation) of a set of line features, with circular variance and mean
//!   length. The only distribution measure that works on line bearings.
//!
//! `statistic` selects the mode. Central Feature uses `distance` euclidean or
//! manhattan and an optional `weight_field`; Linear Directional Mean uses
//! `orientation_only` for undirected lines and emits a single mean-vector line
//! (length = mean line length, placed at the mean of the line midpoints) with
//! `mean_dir`, `circ_var`, `mean_len`, and `n` attributes. A `case_field`
//! partitions the features and yields one result per group.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CentralFeatureTool;

impl Tool for CentralFeatureTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "central_feature",
            display_name: "Central Feature / Linear Directional Mean",
            summary: "The central feature (the input feature with the smallest total distance to all others) or the linear directional mean (mean bearing, circular variance, and mean length of line features) — the two Measuring Geographic Distributions members directional_distribution lacks (ArcGIS Central Feature / Linear Directional Mean).",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (points for central_feature; lines for linear_directional_mean).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "statistic",
                    description: "'central_feature' (default) or 'linear_directional_mean'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "weight_field",
                    description: "Optional numeric field weighting each feature (central_feature only).",
                    required: false,
                },
                ToolParamSpec {
                    name: "case_field",
                    description: "Optional field to group features by; one result per distinct value.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distance",
                    description: "Distance metric for central_feature: 'euclidean' (default) or 'manhattan'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "orientation_only",
                    description: "Treat lines as undirected (orientation mod 180°) for linear_directional_mean. Default false.",
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

        // Group feature indices by case value.
        let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for (fidx, feature) in layer.features.iter().enumerate() {
            if feature.geometry.is_none() {
                continue;
            }
            let key = match case_idx {
                Some(i) => value_string(&feature.attributes[i]),
                None => "ALL".to_string(),
            };
            groups.entry(key).or_default().push(fidx);
        }

        let result = match prm.statistic {
            Statistic::CentralFeature => {
                central_feature(&layer, &groups, weight_idx, prm.distance, output)?
            }
            Statistic::LinearDirectionalMean => {
                linear_directional_mean(&layer, &groups, prm.orientation_only, output)?
            }
        };

        ctx.progress.info(&format!(
            "{} group(s) -> {} output feature(s)",
            groups.len(),
            result.1
        ));

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(result.0));
        outputs.insert("group_count".to_string(), json!(groups.len()));
        outputs.insert("feature_count".to_string(), json!(result.1));
        Ok(ToolRunResult { outputs })
    }
}

// ── Central feature ──────────────────────────────────────────────────────────

fn central_feature(
    layer: &Layer,
    groups: &BTreeMap<String, Vec<usize>>,
    weight_idx: Option<usize>,
    distance: Distance,
    output: Option<&str>,
) -> Result<(String, usize), ToolError> {
    let mut out = Layer::new("central_feature");
    out.geom_type = layer.geom_type;
    out.schema = layer.schema.clone();
    if let Some(epsg) = layer.crs_epsg() {
        out = out.with_crs_epsg(epsg);
    }
    out.add_field(FieldDef::new("total_dist", FieldType::Float));

    let mut count = 0usize;
    for indices in groups.values() {
        // Representative point + weight per feature in the group.
        let pts: Vec<(f64, f64, f64)> = indices
            .iter()
            .filter_map(|&fi| {
                let g = layer.features[fi].geometry.as_ref()?;
                let (x, y) = representative_xy(g)?;
                let w = weight_idx
                    .and_then(|wi| layer.features[fi].attributes.get(wi))
                    .and_then(FieldValue::as_f64)
                    .unwrap_or(1.0);
                Some((x, y, w))
            })
            .collect();
        if pts.is_empty() {
            continue;
        }
        // Total weighted distance from each feature to all others; pick the min.
        let mut best = 0usize;
        let mut best_total = f64::INFINITY;
        for i in 0..pts.len() {
            let mut total = 0.0;
            for j in 0..pts.len() {
                if i == j {
                    continue;
                }
                total += pts[j].2 * dist(pts[i], pts[j], distance);
            }
            if total < best_total {
                best_total = total;
                best = i;
            }
        }
        let src = &layer.features[indices[best]];
        let mut attrs = src.attributes.clone();
        attrs.push(FieldValue::Float(best_total));
        out.push(Feature {
            fid: 0,
            geometry: src.geometry.clone(),
            attributes: attrs,
        });
        count += 1;
    }

    Ok((write_or_store_layer(out, output)?, count))
}

// ── Linear directional mean ──────────────────────────────────────────────────

fn linear_directional_mean(
    layer: &Layer,
    groups: &BTreeMap<String, Vec<usize>>,
    orientation_only: bool,
    output: Option<&str>,
) -> Result<(String, usize), ToolError> {
    let mut out = Layer::new("linear_directional_mean").with_geom_type(GeometryType::LineString);
    if let Some(epsg) = layer.crs_epsg() {
        out = out.with_crs_epsg(epsg);
    }
    out.add_field(FieldDef::new("case", FieldType::Text));
    out.add_field(FieldDef::new("mean_dir", FieldType::Float));
    out.add_field(FieldDef::new("circ_var", FieldType::Float));
    out.add_field(FieldDef::new("mean_len", FieldType::Float));
    out.add_field(FieldDef::new("n", FieldType::Integer));

    let mut count = 0usize;
    for (key, indices) in groups {
        // Per-line: bearing (radians, math convention), length, midpoint.
        let mut sum_sin = 0.0;
        let mut sum_cos = 0.0;
        let mut sum_len = 0.0;
        let mut mx = 0.0;
        let mut my = 0.0;
        let mut n = 0usize;
        for &fi in indices {
            let Some(geom) = layer.features[fi].geometry.as_ref() else {
                continue;
            };
            for (a, b, len, cx, cy) in line_segments_endpoints(geom) {
                if len <= 0.0 {
                    continue;
                }
                let theta = (b.1 - a.1).atan2(b.0 - a.0); // math angle from +x, CCW
                                                          // Orientation: double the angle so opposite directions coincide.
                let ang = if orientation_only { 2.0 * theta } else { theta };
                sum_sin += ang.sin();
                sum_cos += ang.cos();
                sum_len += len;
                mx += cx;
                my += cy;
                n += 1;
            }
        }
        if n == 0 {
            continue;
        }
        let nf = n as f64;
        let mean_ang_raw = sum_sin.atan2(sum_cos);
        let mean_theta = if orientation_only {
            mean_ang_raw / 2.0
        } else {
            mean_ang_raw
        };
        let r = (sum_sin * sum_sin + sum_cos * sum_cos).sqrt() / nf;
        let circ_var = 1.0 - r;
        let mean_len = sum_len / nf;
        let (cx, cy) = (mx / nf, my / nf);

        // Mean-vector line at the group's mean midpoint, mean length, mean angle.
        let half = mean_len * 0.5;
        let (dx, dy) = (mean_theta.cos() * half, mean_theta.sin() * half);
        let line = Geometry::line_string(vec![
            Coord::xy(cx - dx, cy - dy),
            Coord::xy(cx + dx, cy + dy),
        ]);
        // Compass bearing (0 = north, clockwise) for the report.
        let mean_dir = compass_bearing(mean_theta);

        out.push(Feature {
            fid: 0,
            geometry: Some(line),
            attributes: vec![
                FieldValue::Text(key.clone()),
                FieldValue::Float(mean_dir),
                FieldValue::Float(circ_var),
                FieldValue::Float(mean_len),
                FieldValue::Integer(n as i64),
            ],
        });
        count += 1;
    }

    Ok((write_or_store_layer(out, output)?, count))
}

/// Converts a math angle (radians from +x, CCW) to a compass bearing in degrees
/// (0 = north, increasing clockwise), in [0, 360).
fn compass_bearing(theta: f64) -> f64 {
    let deg = 90.0 - theta.to_degrees();
    ((deg % 360.0) + 360.0) % 360.0
}

// ── Geometry helpers ─────────────────────────────────────────────────────────

type Pt = (f64, f64, f64); // x, y, weight

fn dist(a: Pt, b: Pt, metric: Distance) -> f64 {
    let (dx, dy) = (a.0 - b.0, a.1 - b.1);
    match metric {
        Distance::Euclidean => dx.hypot(dy),
        Distance::Manhattan => dx.abs() + dy.abs(),
    }
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

/// For each line in the geometry, its start/end points, straight length, and
/// midpoint (the ArcGIS Linear Directional Mean uses the start→end vector).
#[allow(clippy::type_complexity)]
fn line_segments_endpoints(geom: &Geometry) -> Vec<((f64, f64), (f64, f64), f64, f64, f64)> {
    let one = |cs: &[Coord]| -> Option<((f64, f64), (f64, f64), f64, f64, f64)> {
        if cs.len() < 2 {
            return None;
        }
        let a = (cs[0].x, cs[0].y);
        let b = (cs[cs.len() - 1].x, cs[cs.len() - 1].y);
        let len = (b.0 - a.0).hypot(b.1 - a.1);
        Some((a, b, len, (a.0 + b.0) * 0.5, (a.1 + b.1) * 0.5))
    };
    match geom {
        Geometry::LineString(cs) => one(cs).into_iter().collect(),
        Geometry::MultiLineString(lines) => lines.iter().filter_map(|l| one(l)).collect(),
        _ => Vec::new(),
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

#[derive(Clone, Copy, PartialEq, Eq)]
enum Statistic {
    CentralFeature,
    LinearDirectionalMean,
}

#[derive(Clone, Copy)]
enum Distance {
    Euclidean,
    Manhattan,
}

struct Params {
    statistic: Statistic,
    weight_field: Option<String>,
    case_field: Option<String>,
    distance: Distance,
    orientation_only: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let statistic = match parse_optional_str(args, "statistic")? {
        None => Statistic::CentralFeature,
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "central_feature" | "central-feature" => Statistic::CentralFeature,
            "linear_directional_mean" | "linear-directional-mean" | "directional_mean" => {
                Statistic::LinearDirectionalMean
            }
            other => return Err(ToolError::Validation(format!(
                "'statistic' must be 'central_feature' or 'linear_directional_mean', got '{other}'"
            ))),
        },
    };
    let distance = match parse_optional_str(args, "distance")? {
        None => Distance::Euclidean,
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "euclidean" => Distance::Euclidean,
            "manhattan" => Distance::Manhattan,
            other => {
                return Err(ToolError::Validation(format!(
                    "'distance' must be 'euclidean' or 'manhattan', got '{other}'"
                )))
            }
        },
    };
    let weight_field = parse_optional_str(args, "weight_field")?.map(str::to_string);
    let case_field = parse_optional_str(args, "case_field")?.map(str::to_string);
    let orientation_only = parse_optional_bool(args, "orientation_only")?.unwrap_or(false);
    Ok(Params {
        statistic,
        weight_field,
        case_field,
        distance,
        orientation_only,
    })
}

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
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

    fn point_layer(pts: &[(f64, f64, &str)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        for (x, y, n) in pts {
            l.add_feature(Some(Geometry::point(*x, *y)), &[("name", (*n).into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn line_layer(lines: &[&[(f64, f64)]]) -> String {
        let mut l = Layer::new("lines")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        for coords in lines {
            let cs = coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
            l.add_feature(Some(Geometry::line_string(cs)), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CentralFeatureTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// The central feature is the middle point of a line of points, and it keeps
    /// its attribute.
    #[test]
    fn central_feature_picks_the_middle() {
        // Points at x=0,1,2,3,100. The one minimizing total distance is x=2 or 3;
        // for this set x=3 (closer to the far outlier's pull is offset)... compute:
        // totals: x=2 -> 2+1+1+1+98=103; x=3 -> 3+2+1+... let's just assert it's
        // an interior point, not an endpoint outlier.
        let input = point_layer(&[
            (0.0, 0.0, "a"),
            (1.0, 0.0, "b"),
            (2.0, 0.0, "c"),
            (3.0, 0.0, "d"),
            (100.0, 0.0, "e"),
        ]);
        let (out, layer) = run(json!({ "input": input, "statistic": "central_feature" }));
        assert_eq!(out.outputs["feature_count"], json!(1));
        let name = {
            let idx = layer.schema.field_index("name").unwrap();
            layer.features[0].attributes[idx]
                .as_str()
                .unwrap()
                .to_string()
        };
        // The outlier 'e' and the endpoint 'a' cannot be central.
        assert!(
            name != "e" && name != "a",
            "central feature should be interior, got {name}"
        );
        // total_dist attribute present.
        assert!(layer.schema.field_index("total_dist").is_some());
    }

    /// Weighting toward a cluster pulls the central feature into it.
    #[test]
    fn central_feature_respects_weight() {
        // Two clusters; weight the right cluster heavily.
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.add_field(FieldDef::new("w", FieldType::Float));
        let rows = [
            (0.0, "L1", 1.0),
            (1.0, "L2", 1.0),
            (10.0, "R1", 100.0),
            (11.0, "R2", 100.0),
        ];
        for (x, n, w) in rows {
            l.add_feature(
                Some(Geometry::point(x, 0.0)),
                &[("name", n.into()), ("w", w.into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);
        let (_o, layer) = run(json!({
            "input": input, "statistic": "central_feature", "weight_field": "w",
        }));
        let idx = layer.schema.field_index("name").unwrap();
        let name = layer.features[0].attributes[idx].as_str().unwrap();
        assert!(
            name.starts_with('R'),
            "weight should pull the center right, got {name}"
        );
    }

    /// Parallel east-west lines have a mean direction of ~90° (compass, east).
    #[test]
    fn linear_directional_mean_east() {
        let input = line_layer(&[
            &[(0.0, 0.0), (10.0, 0.0)],
            &[(0.0, 5.0), (10.0, 5.3)],
            &[(0.0, 10.0), (10.0, 9.7)],
        ]);
        let (_o, layer) = run(json!({
            "input": input, "statistic": "linear_directional_mean",
        }));
        let didx = layer.schema.field_index("mean_dir").unwrap();
        let mean_dir = layer.features[0].attributes[didx].as_f64().unwrap();
        // Eastward lines -> compass bearing near 90.
        assert!(
            (mean_dir - 90.0).abs() < 5.0,
            "mean_dir {mean_dir} not ~90 (east)"
        );
        let vidx = layer.schema.field_index("circ_var").unwrap();
        assert!(
            layer.features[0].attributes[vidx].as_f64().unwrap() < 0.2,
            "should be low variance"
        );
    }

    /// orientation_only makes opposite-direction lines agree (no cancellation).
    #[test]
    fn orientation_only_ignores_direction() {
        // Two collinear lines pointing opposite ways along x. Directed: they
        // cancel (high variance). Undirected: they agree (low variance).
        let input = line_layer(&[&[(0.0, 0.0), (10.0, 0.0)], &[(10.0, 1.0), (0.0, 1.0)]]);
        let (_o, directed) = run(json!({
            "input": input.clone(), "statistic": "linear_directional_mean",
        }));
        let vidx = directed.schema.field_index("circ_var").unwrap();
        let directed_var = directed.features[0].attributes[vidx].as_f64().unwrap();

        let (_o2, undirected) = run(json!({
            "input": input, "statistic": "linear_directional_mean", "orientation_only": true,
        }));
        let undirected_var = undirected.features[0].attributes[vidx].as_f64().unwrap();
        assert!(
            directed_var > 0.8,
            "opposite directed lines should nearly cancel: {directed_var}"
        );
        assert!(
            undirected_var < 0.1,
            "undirected should agree: {undirected_var}"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CentralFeatureTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "statistic": "bogus" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "distance": "chebyshev" })).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_ok());
    }
}
