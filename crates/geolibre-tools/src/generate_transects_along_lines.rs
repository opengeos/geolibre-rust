//! GeoLibre tool: generate perpendicular transects along lines.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Generate Transects Along Lines*
//! (Data Management). The bundled `points_along_lines` places points but nothing
//! generates the perpendicular cross-section lines used for shoreline-change
//! analysis (DSAS-style), riparian surveys, and cross-section extraction. Pairs
//! naturally with `interpolate_shape` for terrain cross-sections.
//!
//! Each input line is walked at a fixed `interval`; at each station the local
//! tangent is taken from the containing segment and a transect of total
//! `length` is emitted perpendicular to it, centred on the line (or offset to
//! one side). Each transect carries its parent feature id, the distance along
//! the line, and the transect bearing in degrees. `include_ends` adds transects
//! exactly at the line's start and end.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct GenerateTransectsAlongLinesTool;

impl Tool for GenerateTransectsAlongLinesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "generate_transects_along_lines",
            display_name: "Generate Transects Along Lines",
            summary: "Walk each line at a fixed interval and emit perpendicular transect lines of a given length (centred or offset), for shoreline-change, riparian, and cross-section sampling — like ArcGIS Generate Transects Along Lines.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input line vector layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output line vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "interval",
                    description: "Spacing between transects along each line, in CRS units. Required.",
                    required: true,
                },
                ToolParamSpec {
                    name: "length",
                    description: "Total transect length, in CRS units. Required.",
                    required: true,
                },
                ToolParamSpec {
                    name: "offset",
                    description: "Signed offset of the transect centre from the line, in CRS units (0 = centred; +/- pushes to the left/right). Default 0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "include_ends",
                    description: "Also place a transect exactly at the start and end of each line. Default false.",
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

        let mut out = Layer::new("transects").with_geom_type(GeometryType::LineString);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("line_id", FieldType::Integer));
        out.add_field(FieldDef::new("dist", FieldType::Float));
        out.add_field(FieldDef::new("bearing", FieldType::Float));

        let mut transect_count = 0usize;
        for (fidx, feature) in layer.features.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            for chain in line_chains(geom) {
                transect_count += emit_transects(&mut out, &chain, fidx as i64, &prm)?;
            }
        }

        ctx.progress
            .info(&format!("generated {transect_count} transect(s)"));

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("transect_count".to_string(), json!(transect_count));
        Ok(ToolRunResult { outputs })
    }
}

/// Emits transects along one polyline. Returns the number produced.
fn emit_transects(
    out: &mut Layer,
    pts: &[P],
    line_id: i64,
    prm: &Params,
) -> Result<usize, ToolError> {
    if pts.len() < 2 {
        return Ok(0);
    }
    // Cumulative distance at each vertex.
    let mut cum = vec![0.0f64; pts.len()];
    for i in 1..pts.len() {
        cum[i] = cum[i - 1] + dist(pts[i - 1], pts[i]);
    }
    let total = *cum.last().unwrap();
    if total <= 0.0 {
        return Ok(0);
    }

    // Station distances along the line.
    let mut stations: Vec<f64> = Vec::new();
    if prm.include_ends {
        stations.push(0.0);
    }
    let mut d = if prm.include_ends { prm.interval } else { 0.0 };
    // When not including ends, ArcGIS still starts at the first interval; place
    // the first interior transect at `interval` and step from there.
    if !prm.include_ends {
        d = prm.interval;
    }
    while d < total - 1e-9 {
        stations.push(d);
        d += prm.interval;
    }
    if prm.include_ends {
        stations.push(total);
    }

    let mut count = 0usize;
    for &s in &stations {
        let Some((center, tangent)) = point_and_tangent(pts, &cum, s) else {
            continue;
        };
        // Perpendicular unit vector (rotate tangent +90°): (-ty, tx).
        let (tx, ty) = tangent;
        let (nx, ny) = (-ty, tx);
        let half = prm.length * 0.5;
        // Centre is offset along the normal by `offset`.
        let cx = center.x + nx * prm.offset;
        let cy = center.y + ny * prm.offset;
        let a = Coord::xy(cx - nx * half, cy - ny * half);
        let b = Coord::xy(cx + nx * half, cy + ny * half);
        let bearing = bearing_deg(nx, ny);
        out.add_feature(
            Some(Geometry::line_string(vec![a, b])),
            &[
                ("line_id", line_id.into()),
                ("dist", wbvector::FieldValue::Float(s)),
                ("bearing", wbvector::FieldValue::Float(bearing)),
            ],
        )
        .map_err(|e| ToolError::Execution(format!("failed writing transect: {e}")))?;
        count += 1;
    }
    Ok(count)
}

/// Point at arc-length `s` along the polyline and the local unit tangent.
fn point_and_tangent(pts: &[P], cum: &[f64], s: f64) -> Option<(P, (f64, f64))> {
    let total = *cum.last()?;
    let s = s.clamp(0.0, total);
    // Find the segment containing s.
    let mut i = 0;
    while i + 1 < pts.len() && cum[i + 1] < s {
        i += 1;
    }
    let (a, b) = (pts[i], pts[(i + 1).min(pts.len() - 1)]);
    let seg_len = dist(a, b);
    if seg_len <= 0.0 {
        return None;
    }
    let t = ((s - cum[i]) / seg_len).clamp(0.0, 1.0);
    let p = P {
        x: a.x + (b.x - a.x) * t,
        y: a.y + (b.y - a.y) * t,
    };
    let tangent = ((b.x - a.x) / seg_len, (b.y - a.y) / seg_len);
    Some((p, tangent))
}

/// Compass bearing (0-360, 0 = north/+y, clockwise) of a direction vector.
fn bearing_deg(dx: f64, dy: f64) -> f64 {
    let deg = dx.atan2(dy).to_degrees();
    (deg + 360.0) % 360.0
}

// ── Geometry helpers ─────────────────────────────────────────────────────────

#[derive(Clone, Copy)]
struct P {
    x: f64,
    y: f64,
}

fn dist(a: P, b: P) -> f64 {
    (a.x - b.x).hypot(a.y - b.y)
}

fn line_chains(geom: &Geometry) -> Vec<Vec<P>> {
    let to_pts = |cs: &[Coord]| -> Vec<P> {
        let mut out: Vec<P> = Vec::with_capacity(cs.len());
        for c in cs {
            let p = P { x: c.x, y: c.y };
            if out.last().is_none_or(|l| dist(*l, p) > 1e-12) {
                out.push(p);
            }
        }
        out
    };
    match geom {
        Geometry::LineString(cs) => vec![to_pts(cs)],
        Geometry::MultiLineString(lines) => lines.iter().map(|l| to_pts(l)).collect(),
        _ => Vec::new(),
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    interval: f64,
    length: f64,
    offset: f64,
    include_ends: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let interval = parse_optional_f64(args, "interval")?
        .ok_or_else(|| ToolError::Validation("required parameter 'interval' is missing".into()))?;
    if !(interval > 0.0 && interval.is_finite()) {
        return Err(ToolError::Validation(
            "'interval' must be a positive number".to_string(),
        ));
    }
    let length = parse_optional_f64(args, "length")?
        .ok_or_else(|| ToolError::Validation("required parameter 'length' is missing".into()))?;
    if !(length > 0.0 && length.is_finite()) {
        return Err(ToolError::Validation(
            "'length' must be a positive number".to_string(),
        ));
    }
    let offset = parse_optional_f64(args, "offset")?.unwrap_or(0.0);
    if !offset.is_finite() {
        return Err(ToolError::Validation(
            "'offset' must be a finite number".to_string(),
        ));
    }
    let include_ends = parse_optional_bool(args, "include_ends")?.unwrap_or(false);
    Ok(Params {
        interval,
        length,
        offset,
        include_ends,
    })
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

    fn line_layer(coords: &[(f64, f64)]) -> String {
        let mut l = Layer::new("lines")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        let cs = coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
        l.add_feature(Some(Geometry::line_string(cs)), &[]).unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = GenerateTransectsAlongLinesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn seg(layer: &Layer, idx: usize) -> (Coord, Coord) {
        match layer.features[idx].geometry.as_ref().unwrap() {
            Geometry::LineString(cs) => (cs[0].clone(), cs[1].clone()),
            other => panic!("expected line, got {other:?}"),
        }
    }

    /// Along a horizontal line, transects are vertical, of the right length,
    /// centred on the line.
    #[test]
    fn perpendicular_transects_on_horizontal_line() {
        // Line from (0,0) to (100,0); interval 25, length 10.
        let input = line_layer(&[(0.0, 0.0), (100.0, 0.0)]);
        let (out, layer) = run(json!({ "input": input, "interval": 25.0, "length": 10.0 }));
        // stations at 25,50,75 (ends excluded) -> 3 transects.
        assert_eq!(out.outputs["transect_count"], json!(3));
        for i in 0..3 {
            let (a, b) = seg(&layer, i);
            // Vertical (x equal), length 10, centred at y=0.
            assert!((a.x - b.x).abs() < 1e-6, "transect not vertical");
            assert!(((a.y - b.y).abs() - 10.0).abs() < 1e-6, "wrong length");
            assert!((a.y + b.y).abs() < 1e-6, "not centred on the line");
        }
    }

    /// include_ends adds transects at distance 0 and total length.
    #[test]
    fn include_ends_adds_endpoint_transects() {
        let input = line_layer(&[(0.0, 0.0), (100.0, 0.0)]);
        let (out, layer) = run(json!({
            "input": input, "interval": 25.0, "length": 10.0, "include_ends": true,
        }));
        // 0,25,50,75,100 -> 5 transects.
        assert_eq!(out.outputs["transect_count"], json!(5));
        let didx = layer.schema.field_index("dist").unwrap();
        let mut dists: Vec<f64> = layer
            .iter()
            .map(|f| f.attributes[didx].as_f64().unwrap())
            .collect();
        dists.sort_by(f64::total_cmp);
        assert!((dists[0]).abs() < 1e-9 && (dists[4] - 100.0).abs() < 1e-9);
    }

    /// offset pushes the transect centre to one side of the line.
    #[test]
    fn offset_shifts_centre() {
        let input = line_layer(&[(0.0, 0.0), (100.0, 0.0)]);
        let (_o, layer) = run(json!({
            "input": input, "interval": 50.0, "length": 10.0, "offset": 3.0,
        }));
        // Tangent +x -> normal (-0,+1)=(0,1); offset +3 shifts centre to y=+3,
        // so the transect spans y in [-2, 8].
        let (a, b) = seg(&layer, 0);
        let lo = a.y.min(b.y);
        let hi = a.y.max(b.y);
        assert!(
            (lo + 2.0).abs() < 1e-6 && (hi - 8.0).abs() < 1e-6,
            "offset wrong: [{lo},{hi}]"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            GenerateTransectsAlongLinesTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "interval": 10 })).is_err()); // no length
        assert!(bad(json!({ "input": "a.geojson", "interval": 0, "length": 10 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "interval": 10, "length": 20 })).is_ok());
    }
}
