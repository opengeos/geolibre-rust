//! GeoLibre tool: place points along 3D polylines using true 3D length.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Generate Points Along 3D Lines*
//! (3D Analyst). The bundled `points_along_lines` measures arc length in **2D**,
//! so on a steep 3D line the emitted points are unevenly spaced in reality: the
//! spacing error is `1 / cos(slope)`, which is already ~15% at a 30-degree
//! grade and worse on pipelines, transmission lines, ski runs and mountain
//! trails. That is exactly the geometry where correct spacing matters, because
//! the points usually drive a profile, an inspection schedule or chainage
//! markers.
//!
//! Z is interpolated linearly within the segment a sample lands on, and
//! `add_chainage` reports cumulative **3D** distance from the line start.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Upper bound on emitted points, so a near-zero spacing fails loudly instead
/// of exhausting memory.
const MAX_POINTS: usize = 20_000_000;

pub struct GeneratePointsAlong3dLinesTool;

impl Tool for GeneratePointsAlong3dLinesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "generate_points_along_3d_lines",
            display_name: "Generate Points Along 3D Lines",
            summary: "Generate evenly spaced points along 3D lines measured in true 3D distance (not the 2D plan length the bundled points_along_lines uses), with interpolated Z, optional end points and a cumulative 3D chainage field. Like ArcGIS Generate Points Along 3D Lines.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input 3D polyline features.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output point path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'distance' (default), 'percentage', or 'distance_field'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distance",
                    description: "Spacing in map units for method 'distance'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "percentage",
                    description: "Spacing as a percentage (0-100) of each line's 3D length for method 'percentage'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distance_field",
                    description: "Field giving a per-feature spacing for method 'distance_field'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "include_end_points",
                    description: "Also emit each line's start and end vertices (default false).",
                    required: false,
                },
                ToolParamSpec {
                    name: "add_chainage",
                    description: "Write a CHAINAGE field holding cumulative 3D distance from the line start (default false).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        let method = parse_method(args)?;
        match method {
            Method::Distance => {
                let d = parse_optional_f64(args, "distance")?.ok_or_else(|| {
                    ToolError::Validation(
                        "'distance' is required when 'method' is 'distance'".to_string(),
                    )
                })?;
                if d <= 0.0 {
                    return Err(ToolError::Validation(
                        "'distance' must be positive".to_string(),
                    ));
                }
            }
            Method::Percentage => {
                let p = parse_optional_f64(args, "percentage")?.ok_or_else(|| {
                    ToolError::Validation(
                        "'percentage' is required when 'method' is 'percentage'".to_string(),
                    )
                })?;
                if p <= 0.0 || p > 100.0 {
                    return Err(ToolError::Validation(
                        "'percentage' must be within (0, 100]".to_string(),
                    ));
                }
            }
            Method::DistanceField => {
                parse_optional_str(args, "distance_field")?.ok_or_else(|| {
                    ToolError::Validation(
                        "'distance_field' is required when 'method' is 'distance_field'"
                            .to_string(),
                    )
                })?;
            }
        }
        parse_optional_bool(args, "include_end_points")?;
        parse_optional_bool(args, "add_chainage")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let method = parse_method(args)?;
        let want_ends = parse_optional_bool(args, "include_end_points")?.unwrap_or(false);
        let want_chainage = parse_optional_bool(args, "add_chainage")?.unwrap_or(false);

        let layer = load_input_layer(input)?;
        if layer.features.is_empty() {
            return Err(ToolError::Execution("input has no features".to_string()));
        }
        let dfield = match method {
            Method::DistanceField => {
                let name = parse_optional_str(args, "distance_field")?.unwrap_or_default();
                Some(layer.schema.field_index(name).ok_or_else(|| {
                    ToolError::Validation(format!("distance_field '{name}' not found"))
                })?)
            }
            _ => None,
        };

        let mut out = Layer::new("points_along_3d_lines").with_geom_type(GeometryType::Point);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        let names: Vec<String> = layer
            .schema
            .fields()
            .iter()
            .map(|f| f.name.clone())
            .collect();
        // Chainage/route layers routinely already carry LINE_ID or CHAINAGE.
        // Appending ours unconditionally would duplicate the name, and the
        // by-name attribute write below would then land in whichever index
        // resolved first. Suffix ours until unique instead.
        let unique = |base: &str| -> String {
            if !names.iter().any(|n| n == base) {
                return base.to_string();
            }
            (1..)
                .map(|k| format!("{base}_{k}"))
                .find(|c| !names.iter().any(|n| n == c))
                .expect("an unused suffix always exists")
        };
        let line_id_name = unique("LINE_ID");
        let seq_name = unique("POINT_SEQ");
        let chainage_name = unique("CHAINAGE");

        for fd in layer.schema.fields() {
            out.add_field(fd.clone());
        }
        out.add_field(FieldDef::new(line_id_name.as_str(), FieldType::Integer));
        out.add_field(FieldDef::new(seq_name.as_str(), FieldType::Integer));
        if want_chainage {
            out.add_field(FieldDef::new(chainage_name.as_str(), FieldType::Float));
        }

        let mut emitted = 0usize;
        let mut skipped = 0usize;
        let mut total_len_3d = 0.0_f64;
        let mut total_len_2d = 0.0_f64;

        for (fid, feat) in layer.iter().enumerate() {
            let Some(g) = &feat.geometry else {
                skipped += 1;
                continue;
            };
            let paths = line_paths(g);
            if paths.is_empty() {
                skipped += 1;
                continue;
            }
            for path in &paths {
                let len3 = path_length_3d(path);
                total_len_3d += len3;
                total_len_2d += path_length_2d(path);
                if len3 <= 0.0 {
                    continue;
                }
                let step = match method {
                    Method::Distance => parse_optional_f64(args, "distance")?.unwrap_or(0.0),
                    Method::Percentage => {
                        len3 * parse_optional_f64(args, "percentage")?.unwrap_or(0.0) / 100.0
                    }
                    Method::DistanceField => feat
                        .attributes
                        .get(dfield.expect("distance_field resolved above"))
                        .and_then(FieldValue::as_f64)
                        .unwrap_or(0.0),
                };
                if step <= 0.0 || !step.is_finite() {
                    skipped += 1;
                    continue;
                }
                if (len3 / step) as usize + emitted > MAX_POINTS {
                    return Err(ToolError::Execution(format!(
                        "spacing {step} would emit more than {MAX_POINTS} points; \
                         increase 'distance' or 'percentage'"
                    )));
                }

                for (seq, (p, chain)) in sample_3d(path, step, want_ends).into_iter().enumerate() {
                    let mut attrs: Vec<(&str, FieldValue)> = names
                        .iter()
                        .enumerate()
                        .map(|(i, nm)| {
                            (
                                nm.as_str(),
                                feat.attributes.get(i).cloned().unwrap_or(FieldValue::Null),
                            )
                        })
                        .collect();
                    attrs.push((line_id_name.as_str(), FieldValue::Integer(fid as i64)));
                    attrs.push((seq_name.as_str(), FieldValue::Integer(seq as i64)));
                    if want_chainage {
                        attrs.push((chainage_name.as_str(), FieldValue::Float(chain)));
                    }
                    out.add_feature(Some(Geometry::point_z(p[0], p[1], p[2])), &attrs)
                        .map_err(|e| ToolError::Execution(format!("failed adding point: {e}")))?;
                    emitted += 1;
                }
            }
        }
        ctx.progress
            .info(&format!("generated {emitted} point(s) along 3D lines"));

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("point_count".to_string(), json!(emitted));
        outputs.insert("skipped_features".to_string(), json!(skipped));
        outputs.insert("total_length_3d".to_string(), json!(total_len_3d));
        outputs.insert("total_length_2d".to_string(), json!(total_len_2d));
        // How much longer the true 3D path is than its plan projection — the
        // number that quantifies what the bundled 2D tool would have got wrong.
        if total_len_2d > 0.0 {
            outputs.insert(
                "length_ratio_3d_2d".to_string(),
                json!(total_len_3d / total_len_2d),
            );
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── Geometry ────────────────────────────────────────────────────────────────

type P3 = [f64; 3];

fn to_p3(c: &Coord) -> P3 {
    [c.x, c.y, c.z.unwrap_or(0.0)]
}

/// Extracts the linear paths from a geometry (points and polygons are ignored).
fn line_paths(g: &Geometry) -> Vec<Vec<P3>> {
    match g {
        Geometry::LineString(cs) if cs.len() >= 2 => vec![cs.iter().map(to_p3).collect()],
        Geometry::MultiLineString(ls) => ls
            .iter()
            .filter(|l| l.len() >= 2)
            .map(|l| l.iter().map(to_p3).collect())
            .collect(),
        Geometry::GeometryCollection(gs) => gs.iter().flat_map(line_paths).collect(),
        _ => Vec::new(),
    }
}

fn seg_len_3d(a: &P3, b: &P3) -> f64 {
    ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2) + (b[2] - a[2]).powi(2)).sqrt()
}

fn path_length_3d(p: &[P3]) -> f64 {
    p.windows(2).map(|w| seg_len_3d(&w[0], &w[1])).sum()
}

fn path_length_2d(p: &[P3]) -> f64 {
    p.windows(2)
        .map(|w| (w[1][0] - w[0][0]).hypot(w[1][1] - w[0][1]))
        .sum()
}

/// Samples a 3D path every `step` units of **3D** distance, returning each
/// point with its cumulative 3D chainage.
fn sample_3d(path: &[P3], step: f64, include_ends: bool) -> Vec<(P3, f64)> {
    let mut out: Vec<(P3, f64)> = Vec::new();
    if path.len() < 2 || step <= 0.0 {
        return out;
    }
    let total = path_length_3d(path);
    if include_ends {
        out.push((path[0], 0.0));
    }

    let mut acc = 0.0_f64;
    let mut next = step;
    for w in path.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        let seg = seg_len_3d(a, b);
        if seg <= 0.0 {
            continue;
        }
        while next <= acc + seg + 1e-9 {
            let t = ((next - acc) / seg).clamp(0.0, 1.0);
            let p = [
                a[0] + t * (b[0] - a[0]),
                a[1] + t * (b[1] - a[1]),
                a[2] + t * (b[2] - a[2]),
            ];
            // With ends requested, a sample landing exactly on the final vertex
            // would duplicate the end point appended below.
            let at_end = (next - total).abs() <= 1e-9;
            if !(include_ends && at_end) {
                out.push((p, next));
            }
            next += step;
        }
        acc += seg;
    }

    if include_ends {
        out.push((path[path.len() - 1], total));
    }
    out
}

// ── Params ──────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Distance,
    Percentage,
    DistanceField,
}

fn parse_method(args: &ToolArgs) -> Result<Method, ToolError> {
    match args
        .get("method")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        None | Some("") | Some("distance") | Some("by_distance") => Ok(Method::Distance),
        Some("percentage") | Some("by_percentage") => Ok(Method::Percentage),
        Some("distance_field") | Some("by_distance_field") => Ok(Method::DistanceField),
        Some(o) => Err(ToolError::Validation(format!(
            "'method' must be 'distance', 'percentage' or 'distance_field', got '{o}'"
        ))),
    }
}

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
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

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn line(pts: &[(f64, f64, f64)], spacing: Option<f64>) -> String {
        let mut l = Layer::new("l")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        if spacing.is_some() {
            l.add_field(FieldDef::new("sp", FieldType::Float));
        }
        let coords: Vec<Coord> = pts.iter().map(|(x, y, z)| Coord::xyz(*x, *y, *z)).collect();
        let a: Vec<(&str, FieldValue)> = match spacing {
            Some(s) => vec![("sp", FieldValue::Float(s))],
            None => vec![],
        };
        l.add_feature(Some(Geometry::line_string(coords)), &a).unwrap();
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = GeneratePointsAlong3dLinesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    #[test]
    fn spacing_is_measured_in_3d_not_plan_distance() {
        // 3-4-5 line: plan length 40, rise 30, true 3D length 50.
        // At spacing 10, a 3D-correct tool emits 5 points (10..50);
        // a 2D tool would emit only 4 (10..40).
        let input = line(&[(0.0, 0.0, 0.0), (40.0, 0.0, 30.0)], None);
        let (out, layer) = run(json!({ "input": input, "distance": 10.0 }));
        assert_eq!(out.outputs["point_count"], json!(5));
        assert!((out.outputs["total_length_3d"].as_f64().unwrap() - 50.0).abs() < 1e-9);
        assert!((out.outputs["total_length_2d"].as_f64().unwrap() - 40.0).abs() < 1e-9);
        assert!((out.outputs["length_ratio_3d_2d"].as_f64().unwrap() - 1.25).abs() < 1e-9);
        // Consecutive points are 10 apart in 3D.
        let pts: Vec<P3> = layer
            .iter()
            .map(|f| match f.geometry.as_ref().unwrap() {
                Geometry::Point(c) => [c.x, c.y, c.z.unwrap_or(0.0)],
                other => panic!("unexpected {other:?}"),
            })
            .collect();
        for w in pts.windows(2) {
            assert!((seg_len_3d(&w[0], &w[1]) - 10.0).abs() < 1e-9);
        }
    }

    #[test]
    fn z_is_interpolated_along_the_segment() {
        let input = line(&[(0.0, 0.0, 0.0), (0.0, 0.0, 100.0)], None);
        let (_o, layer) = run(json!({ "input": input, "distance": 25.0 }));
        let zs: Vec<f64> = layer
            .iter()
            .map(|f| match f.geometry.as_ref().unwrap() {
                Geometry::Point(c) => c.z.unwrap(),
                other => panic!("unexpected {other:?}"),
            })
            .collect();
        assert_eq!(zs.len(), 4);
        for (i, z) in zs.iter().enumerate() {
            assert!((z - (25.0 * (i + 1) as f64)).abs() < 1e-9);
        }
    }

    #[test]
    fn chainage_reports_cumulative_3d_distance() {
        let input = line(&[(0.0, 0.0, 0.0), (40.0, 0.0, 30.0)], None);
        let (_o, layer) = run(json!({
            "input": input, "distance": 10.0, "add_chainage": true
        }));
        let c = layer.schema.field_index("CHAINAGE").unwrap();
        let ch: Vec<f64> = layer
            .iter()
            .map(|f| f.attributes[c].as_f64().unwrap())
            .collect();
        assert_eq!(ch, vec![10.0, 20.0, 30.0, 40.0, 50.0]);
    }

    #[test]
    fn end_points_are_added_without_duplicating_the_last_sample() {
        // Length 50, spacing 10: samples land exactly on the end vertex.
        let input = line(&[(0.0, 0.0, 0.0), (40.0, 0.0, 30.0)], None);
        let (out, _l) = run(json!({
            "input": input, "distance": 10.0, "include_end_points": true
        }));
        // 5 interior/aligned samples, minus the one coincident with the end,
        // plus explicit start and end = 6.
        assert_eq!(out.outputs["point_count"], json!(6));
    }

    #[test]
    fn percentage_method_scales_to_each_line() {
        let input = line(&[(0.0, 0.0, 0.0), (40.0, 0.0, 30.0)], None);
        let (out, _l) = run(json!({
            "input": input, "method": "percentage", "percentage": 25.0
        }));
        // 25% of 50 = 12.5 -> samples at 12.5, 25, 37.5, 50.
        assert_eq!(out.outputs["point_count"], json!(4));
    }

    #[test]
    fn distance_field_gives_per_feature_spacing() {
        let input = line(&[(0.0, 0.0, 0.0), (40.0, 0.0, 30.0)], Some(5.0));
        let (out, _l) = run(json!({
            "input": input, "method": "distance_field", "distance_field": "sp"
        }));
        assert_eq!(out.outputs["point_count"], json!(10));
    }

    #[test]
    fn flat_lines_match_the_plan_length() {
        let input = line(&[(0.0, 0.0, 0.0), (30.0, 0.0, 0.0)], None);
        let (out, _l) = run(json!({ "input": input, "distance": 10.0 }));
        assert_eq!(out.outputs["length_ratio_3d_2d"].as_f64(), Some(1.0));
        assert_eq!(out.outputs["point_count"], json!(3));
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            GeneratePointsAlong3dLinesTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "l.shp" })).is_err()); // distance missing
        assert!(bad(json!({ "input": "l.shp", "distance": 0 })).is_err());
        assert!(bad(json!({ "input": "l.shp", "method": "percentage", "percentage": 0 })).is_err());
        assert!(bad(json!({ "input": "l.shp", "method": "distance_field" })).is_err());
        assert!(bad(json!({ "input": "l.shp", "distance": 5 })).is_ok());
    }
}
