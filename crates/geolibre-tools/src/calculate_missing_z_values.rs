//! GeoLibre tool: fill placeholder / missing Z values along 3D features.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Calculate Missing Z Values* (3D
//! Analyst). Survey data, GPS tracks, and `interpolate_shape` runs over a DEM
//! with holes routinely leave gaps in a feature's Z values, encoded as a
//! sentinel (`-9999`) or as no Z at all. Every downstream 3D tool in the repo —
//! `add_surface_information`, `polygon_volume`, `surface_volume` — then either
//! fails on those vertices or silently returns a wrong answer.
//!
//! Nothing repairs this today. The shipped `fill_missing_values` imputes
//! *attribute* nulls from neighbours and never touches geometry;
//! `adjust_3d_z` rescales Z values that already exist; `repair_geometry` fixes
//! planimetric topology only.
//!
//! Interior gaps are filled by interpolating linearly against **cumulative 2D
//! distance along the feature**, not against vertex index — index interpolation
//! distorts Z badly on unevenly spaced geometry, which is exactly what GPS
//! tracks are. Leading and trailing gaps have only one valid neighbour and are
//! filled by extending it (`extrapolate = true`) or left alone and counted.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Geometry, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// How to fill a run of missing Z values.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Method {
    Linear,
    Nearest,
}

pub struct CalculateMissingZValuesTool;

impl Tool for CalculateMissingZValuesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "calculate_missing_z_values",
            display_name: "Calculate Missing Z Values",
            summary: "Fill placeholder or absent Z values along 3D features by interpolating from the valid Z values on either side, like ArcGIS Calculate Missing Z Values.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input feature layer with 3D vertices.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "placeholder",
                    description: "Sentinel value marking a missing Z (e.g. -9999). Vertices with no Z at all are always treated as missing.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "Fill method: 'linear' (default; interpolate by distance along the feature) or 'nearest'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "extrapolate",
                    description: "Fill leading/trailing gaps by extending the nearest valid vertex (default false; unfilled vertices are counted instead).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_method(args)?;
        parse_optional_f64(args, "placeholder")?;
        parse_optional_bool(args, "extrapolate")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let placeholder = parse_optional_f64(args, "placeholder")?;
        let method = parse_method(args)?;
        let extrapolate = parse_optional_bool(args, "extrapolate")?.unwrap_or(false);

        let layer = load_input_layer(input)?;
        let mut out = Layer::new(layer.name.clone());
        out.crs = layer.crs.clone();
        out.geom_type = layer.geom_type;
        for fd in layer.schema.fields().iter() {
            out.add_field(fd.clone());
        }

        ctx.progress
            .info(&format!("repairing Z on {} feature(s)", layer.len()));

        let mut filled_total = 0usize;
        let mut unfilled_total = 0usize;
        let mut features_touched = 0usize;

        for (fi, feat) in layer.features.iter().enumerate() {
            let geom = match feat.geometry.as_ref() {
                Some(g) => {
                    let mut filled = 0usize;
                    let mut unfilled = 0usize;
                    let repaired = repair_geometry_z(
                        g,
                        placeholder,
                        method,
                        extrapolate,
                        &mut filled,
                        &mut unfilled,
                    );
                    if filled > 0 {
                        features_touched += 1;
                    }
                    filled_total += filled;
                    unfilled_total += unfilled;
                    Some(repaired)
                }
                None => None,
            };
            let fields: Vec<(&str, wbvector::FieldValue)> = layer
                .schema
                .fields()
                .iter()
                .enumerate()
                .map(|(i, fd)| (fd.name.as_str(), feat.attributes[i].clone()))
                .collect();
            out.add_feature(geom, &fields)
                .map_err(|e| ToolError::Execution(format!("failed writing feature: {e}")))?;
            ctx.progress
                .progress((fi as f64 + 1.0) / layer.len().max(1) as f64);
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(layer.len()));
        outputs.insert("features_repaired".to_string(), json!(features_touched));
        outputs.insert("vertices_filled".to_string(), json!(filled_total));
        // Surfaced so a caller can tell "nothing to do" from "could not fill".
        outputs.insert("vertices_unfilled".to_string(), json!(unfilled_total));
        Ok(ToolRunResult { outputs })
    }
}

/// Applies the Z repair to every vertex sequence in a geometry.
fn repair_geometry_z(
    geom: &Geometry,
    placeholder: Option<f64>,
    method: Method,
    extrapolate: bool,
    filled: &mut usize,
    unfilled: &mut usize,
) -> Geometry {
    match geom {
        Geometry::Point(c) => {
            // A lone point has no neighbour to interpolate from.
            if is_missing(c, placeholder) {
                *unfilled += 1;
            }
            Geometry::Point(c.clone())
        }
        Geometry::MultiPoint(cs) => {
            for c in cs {
                if is_missing(c, placeholder) {
                    *unfilled += 1;
                }
            }
            Geometry::MultiPoint(cs.clone())
        }
        Geometry::LineString(cs) => Geometry::LineString(repair_run(
            cs,
            placeholder,
            method,
            extrapolate,
            filled,
            unfilled,
        )),
        Geometry::MultiLineString(parts) => Geometry::MultiLineString(
            parts
                .iter()
                .map(|cs| repair_run(cs, placeholder, method, extrapolate, filled, unfilled))
                .collect(),
        ),
        Geometry::Polygon {
            exterior,
            interiors,
        } => Geometry::Polygon {
            exterior: wbvector::Ring::new(repair_ring(
                exterior.coords(),
                placeholder,
                method,
                extrapolate,
                filled,
                unfilled,
            )),
            interiors: interiors
                .iter()
                .map(|r| {
                    wbvector::Ring::new(repair_ring(
                        r.coords(),
                        placeholder,
                        method,
                        extrapolate,
                        filled,
                        unfilled,
                    ))
                })
                .collect(),
        },
        Geometry::MultiPolygon(parts) => Geometry::MultiPolygon(
            parts
                .iter()
                .map(|(e, hs)| {
                    (
                        wbvector::Ring::new(repair_ring(
                            e.coords(),
                            placeholder,
                            method,
                            extrapolate,
                            filled,
                            unfilled,
                        )),
                        hs.iter()
                            .map(|r| {
                                wbvector::Ring::new(repair_ring(
                                    r.coords(),
                                    placeholder,
                                    method,
                                    extrapolate,
                                    filled,
                                    unfilled,
                                ))
                            })
                            .collect(),
                    )
                })
                .collect(),
        ),
        // Recurse rather than passing through: the other new 3D tools
        // (buffer_3d, near_3d) both descend into collections, and silently
        // reporting "nothing to do" for a collection full of gaps is worse
        // than doing the work.
        Geometry::GeometryCollection(gs) => Geometry::GeometryCollection(
            gs.iter()
                .map(|g| repair_geometry_z(g, placeholder, method, extrapolate, filled, unfilled))
                .collect(),
        ),
    }
}

fn is_missing(c: &Coord, placeholder: Option<f64>) -> bool {
    match c.z {
        None => true,
        Some(z) => !z.is_finite() || placeholder.map(|p| (z - p).abs() < 1e-9).unwrap_or(false),
    }
}

fn repair_ring(
    coords: &[Coord],
    placeholder: Option<f64>,
    method: Method,
    extrapolate: bool,
    filled: &mut usize,
    unfilled: &mut usize,
) -> Vec<Coord> {
    if coords.is_empty() {
        return Vec::new();
    }
    let Some(start) = coords.iter().position(|c| !is_missing(c, placeholder)) else {
        *unfilled += coords.len();
        return coords.to_vec();
    };
    let mut closed: Vec<Coord> = coords[start..]
        .iter()
        .chain(coords[..start].iter())
        .cloned()
        .collect();
    closed.push(closed[0].clone());
    let mut repaired = repair_run(&closed, placeholder, method, extrapolate, filled, unfilled);
    repaired.pop();
    repaired.rotate_right(start);
    repaired
}

/// Fills missing Z across one vertex sequence, interpolating against cumulative
/// 2D distance so unevenly spaced vertices are handled correctly.
fn repair_run(
    coords: &[Coord],
    placeholder: Option<f64>,
    method: Method,
    extrapolate: bool,
    filled: &mut usize,
    unfilled: &mut usize,
) -> Vec<Coord> {
    let n = coords.len();
    let mut out: Vec<Coord> = coords.to_vec();
    if n == 0 {
        return out;
    }

    // Cumulative planimetric distance along the sequence.
    let mut dist = vec![0.0_f64; n];
    for i in 1..n {
        let dx = coords[i].x - coords[i - 1].x;
        let dy = coords[i].y - coords[i - 1].y;
        dist[i] = dist[i - 1] + (dx * dx + dy * dy).sqrt();
    }

    let valid: Vec<usize> = (0..n)
        .filter(|&i| !is_missing(&coords[i], placeholder))
        .collect();
    if valid.is_empty() {
        // Nothing to interpolate from — every missing vertex stays missing.
        *unfilled += n;
        return out;
    }

    for i in 0..n {
        if !is_missing(&coords[i], placeholder) {
            continue;
        }
        // Nearest valid neighbour on each side.
        let split = valid.partition_point(|&v| v < i);
        let before = split.checked_sub(1).map(|j| valid[j]);
        let after = valid.get(split).copied().filter(|&v| v > i);

        let z = match (before, after) {
            (Some(a), Some(b)) => {
                let za = coords[a].z.unwrap();
                let zb = coords[b].z.unwrap();
                match method {
                    Method::Nearest => {
                        if (dist[i] - dist[a]) <= (dist[b] - dist[i]) {
                            Some(za)
                        } else {
                            Some(zb)
                        }
                    }
                    Method::Linear => {
                        let span = dist[b] - dist[a];
                        if span.abs() < f64::EPSILON {
                            Some(za)
                        } else {
                            let t = (dist[i] - dist[a]) / span;
                            Some(za + (zb - za) * t)
                        }
                    }
                }
            }
            // Leading / trailing gap: only one side is available.
            (Some(a), None) if extrapolate => Some(coords[a].z.unwrap()),
            (None, Some(b)) if extrapolate => Some(coords[b].z.unwrap()),
            _ => None,
        };

        match z {
            Some(z) => {
                out[i] = Coord::xyz(coords[i].x, coords[i].y, z);
                out[i].m = coords[i].m;
                *filled += 1;
            }
            None => *unfilled += 1,
        }
    }
    out
}

// ── parameter parsing ────────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    match args.get(key).and_then(Value::as_str) {
        Some(s) if !s.trim().is_empty() => Ok(s),
        _ => Err(ToolError::Validation(format!(
            "missing required string parameter '{key}'"
        ))),
    }
}

fn parse_method(args: &ToolArgs) -> Result<Method, ToolError> {
    match args.get("method").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("linear") => Ok(Method::Linear),
        Some("nearest") => Ok(Method::Nearest),
        Some(o) => Err(ToolError::Validation(format!(
            "'method' must be linear or nearest, got '{o}'"
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

fn parse_optional_bool(args: &ToolArgs, k: &str) -> Result<Option<bool>, ToolError> {
    match args.get(k) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{k}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{k}' must be a boolean"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn line_layer(coords: Vec<Coord>) -> String {
        let mut l = Layer::new("l")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_feature(Some(Geometry::LineString(coords)), &[])
            .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Vec<Coord>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CalculateMissingZValuesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        let cs = match layer.features[0].geometry.as_ref().unwrap() {
            Geometry::LineString(cs) => cs.clone(),
            other => panic!("expected LineString, got {other:?}"),
        };
        (out, cs)
    }

    /// Interior gap fills linearly by distance along the line.
    #[test]
    fn fills_interior_gap_linearly() {
        let cs = vec![
            Coord::xyz(0.0, 0.0, 0.0),
            Coord::xyz(10.0, 0.0, -9999.0),
            Coord::xyz(20.0, 0.0, 100.0),
        ];
        let (out, got) = run(json!({ "input": line_layer(cs), "placeholder": -9999.0 }));
        assert!((got[1].z.unwrap() - 50.0).abs() < 1e-9);
        assert_eq!(out.outputs["vertices_filled"], json!(1));
        assert_eq!(out.outputs["vertices_unfilled"], json!(0));
    }

    /// THE reason to interpolate by distance, not index: unevenly spaced
    /// vertices must weight by geometry.
    #[test]
    fn interpolates_by_distance_not_index() {
        // Gap sits at x=2 of a 0..10 span -> 20% of the way, not 50%.
        let cs = vec![
            Coord::xyz(0.0, 0.0, 0.0),
            Coord::xyz(2.0, 0.0, -9999.0),
            Coord::xyz(10.0, 0.0, 100.0),
        ];
        let (_, got) = run(json!({ "input": line_layer(cs), "placeholder": -9999.0 }));
        assert!(
            (got[1].z.unwrap() - 20.0).abs() < 1e-9,
            "index interpolation would give 50, distance gives 20; got {}",
            got[1].z.unwrap()
        );
    }

    /// Vertices with no Z at all count as missing.
    #[test]
    fn treats_absent_z_as_missing() {
        let cs = vec![
            Coord::xyz(0.0, 0.0, 0.0),
            Coord::xy(10.0, 0.0),
            Coord::xyz(20.0, 0.0, 20.0),
        ];
        let (out, got) = run(json!({ "input": line_layer(cs) }));
        assert!((got[1].z.unwrap() - 10.0).abs() < 1e-9);
        assert_eq!(out.outputs["vertices_filled"], json!(1));
    }

    /// Trailing gaps stay unfilled by default and are reported, not hidden.
    #[test]
    fn trailing_gap_unfilled_by_default() {
        let cs = vec![
            Coord::xyz(0.0, 0.0, 5.0),
            Coord::xyz(10.0, 0.0, 15.0),
            Coord::xyz(20.0, 0.0, -9999.0),
        ];
        let (out, got) = run(json!({ "input": line_layer(cs.clone()), "placeholder": -9999.0 }));
        assert_eq!(out.outputs["vertices_unfilled"], json!(1));
        assert_eq!(got[2].z, Some(-9999.0), "left as the placeholder");

        // With extrapolate it extends the last valid Z.
        let (out2, got2) = run(json!({
            "input": line_layer(cs), "placeholder": -9999.0, "extrapolate": true
        }));
        assert_eq!(out2.outputs["vertices_unfilled"], json!(0));
        assert!((got2[2].z.unwrap() - 15.0).abs() < 1e-9);
    }

    /// nearest snaps to the closer neighbour rather than blending.
    #[test]
    fn nearest_method_snaps() {
        let cs = vec![
            Coord::xyz(0.0, 0.0, 0.0),
            Coord::xyz(1.0, 0.0, -9999.0),
            Coord::xyz(10.0, 0.0, 100.0),
        ];
        let (_, got) = run(json!({
            "input": line_layer(cs), "placeholder": -9999.0, "method": "nearest"
        }));
        assert_eq!(got[1].z, Some(0.0));
    }

    /// A run with no valid Z anywhere cannot be filled and must say so.
    #[test]
    fn all_missing_reports_unfilled() {
        let cs = vec![Coord::xy(0.0, 0.0), Coord::xy(1.0, 0.0)];
        let (out, _) = run(json!({ "input": line_layer(cs) }));
        assert_eq!(out.outputs["vertices_filled"], json!(0));
        assert_eq!(out.outputs["vertices_unfilled"], json!(2));
    }

    /// Valid Z values are never rewritten.
    #[test]
    fn preserves_valid_z() {
        let cs = vec![
            Coord::xyz(0.0, 0.0, 3.5),
            Coord::xyz(10.0, 0.0, -9999.0),
            Coord::xyz(20.0, 0.0, 7.5),
        ];
        let (_, got) = run(json!({ "input": line_layer(cs), "placeholder": -9999.0 }));
        assert_eq!(got[0].z, Some(3.5));
        assert_eq!(got[2].z, Some(7.5));
    }

    #[test]
    fn rejects_bad_parameters() {
        let p = line_layer(vec![Coord::xyz(0.0, 0.0, 1.0)]);
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CalculateMissingZValuesTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({ "input": p, "method": "cubic" })));
        assert!(bad(json!({ "input": p, "placeholder": "abc" })));
        assert!(bad(json!({ "input": p, "extrapolate": "maybe" })));
    }
}
