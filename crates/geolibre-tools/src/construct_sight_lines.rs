//! GeoLibre tool: build observer-to-target sight lines.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Construct Sight Lines* (3D Analyst).
//! The catalog could already *consume* sight lines but never *produce* them:
//! GeoLibre's `line_of_sight` takes lines and reports visibility along them,
//! and the bundled `viewshed` / `visibility_index` work from observer points
//! against a surface. Nothing generated the observer-to-target geometry, so
//! `line_of_sight` was unusable without hand-building its input layer.
//!
//! This is the missing first step of the visibility pipeline, and it is pure
//! geometry — no surface sampling and no raster work, which keeps it cheap and
//! composable:
//!
//! ```text
//! construct_sight_lines -> line_of_sight -> (visible/obstructed segments)
//! ```
//!
//! Linear and areal targets are sampled along their boundary at
//! `sample_distance`, so a ridgeline or a building footprint yields a fan of
//! sight lines rather than a single line to its centroid.

use std::collections::BTreeMap;
use std::collections::HashMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Guard against an accidental full cross product blowing up memory. ArcGIS
/// silently produces the same explosion; failing loudly is more useful.
const MAX_SIGHT_LINES: usize = 5_000_000;

pub struct ConstructSightLinesTool;

impl Tool for ConstructSightLinesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "construct_sight_lines",
            display_name: "Construct Sight Lines",
            summary: "Build sight-line features between observer points and target features, with observer/target height offsets, sampling along linear and areal targets, optional azimuth and vertical-angle fields, and 2D or 3D distance. Like ArcGIS Construct Sight Lines.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "observers",
                    description: "Observer point features.",
                    required: true,
                },
                ToolParamSpec {
                    name: "targets",
                    description: "Target point, line or polygon features.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output sight-line path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "observer_height_field",
                    description: "Optional field giving each observer's height offset (added to its Z).",
                    required: false,
                },
                ToolParamSpec {
                    name: "target_height_field",
                    description: "Optional field giving each target's height offset (added to its Z).",
                    required: false,
                },
                ToolParamSpec {
                    name: "join_field",
                    description: "Optional field present on both inputs; only observers and targets sharing a value are paired. Without it, every observer pairs with every target.",
                    required: false,
                },
                ToolParamSpec {
                    name: "sample_distance",
                    description: "Spacing in map units for target points generated along line/polygon boundaries (default: one point per vertex).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_direction",
                    description: "Add AZIMUTH and VERT_ANGLE fields (default false).",
                    required: false,
                },
                ToolParamSpec {
                    name: "distance_method",
                    description: "'2d' (default) or '3d' for the DISTANCE field.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "observers")?;
        require_str(args, "targets")?;
        parse_distance_method(args)?;
        parse_optional_bool(args, "output_direction")?;
        if let Some(d) = parse_optional_f64(args, "sample_distance")? {
            if d <= 0.0 {
                return Err(ToolError::Validation(
                    "'sample_distance' must be positive".to_string(),
                ));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let obs_path = require_str(args, "observers")?;
        let tgt_path = require_str(args, "targets")?;
        let output = parse_optional_str(args, "output")?;
        let obs_h_field = parse_optional_str(args, "observer_height_field")?;
        let tgt_h_field = parse_optional_str(args, "target_height_field")?;
        let join_field = parse_optional_str(args, "join_field")?;
        let sample_distance = parse_optional_f64(args, "sample_distance")?;
        let want_dir = parse_optional_bool(args, "output_direction")?.unwrap_or(false);
        let method = parse_distance_method(args)?;

        let observers = load_input_layer(obs_path)?;
        let targets = load_input_layer(tgt_path)?;
        if observers.features.is_empty() {
            return Err(ToolError::Execution(
                "observers layer has no features".to_string(),
            ));
        }
        if targets.features.is_empty() {
            return Err(ToolError::Execution(
                "targets layer has no features".to_string(),
            ));
        }

        let obs_h_idx = optional_field(&observers, obs_h_field, "observer_height_field")?;
        let tgt_h_idx = optional_field(&targets, tgt_h_field, "target_height_field")?;
        let obs_j_idx = optional_field(&observers, join_field, "join_field")?;
        let tgt_j_idx = optional_field(&targets, join_field, "join_field")?;

        // Observers: one point each (non-point observers use their centroid).
        let mut obs_pts: Vec<(usize, [f64; 3], String)> = Vec::new();
        for (i, f) in observers.iter().enumerate() {
            let Some(g) = &f.geometry else { continue };
            let Some(mut p) = first_point(g) else { continue };
            p[2] += obs_h_idx.map_or(0.0, |k| {
                f.attributes.get(k).and_then(FieldValue::as_f64).unwrap_or(0.0)
            });
            let key = obs_j_idx.map_or_else(String::new, |k| key_of(f.attributes.get(k)));
            obs_pts.push((i, p, key));
        }

        // Targets: points pass through; lines/polygons are sampled along the
        // boundary so an areal target yields a fan of sight lines.
        let mut tgt_pts: Vec<(usize, [f64; 3], String)> = Vec::new();
        for (i, f) in targets.iter().enumerate() {
            let Some(g) = &f.geometry else { continue };
            let dz = tgt_h_idx.map_or(0.0, |k| {
                f.attributes.get(k).and_then(FieldValue::as_f64).unwrap_or(0.0)
            });
            let key = tgt_j_idx.map_or_else(String::new, |k| key_of(f.attributes.get(k)));
            for mut p in target_points(g, sample_distance) {
                p[2] += dz;
                tgt_pts.push((i, p, key.clone()));
            }
        }
        if obs_pts.is_empty() || tgt_pts.is_empty() {
            return Err(ToolError::Execution(
                "no usable observer/target geometry found".to_string(),
            ));
        }

        // Pair up: grouped by join value, or the full cross product.
        let mut pairs: Vec<(usize, usize)> = Vec::new();
        if join_field.is_some() {
            let mut by_key: HashMap<&str, Vec<usize>> = HashMap::new();
            for (ti, t) in tgt_pts.iter().enumerate() {
                by_key.entry(t.2.as_str()).or_default().push(ti);
            }
            for (oi, o) in obs_pts.iter().enumerate() {
                if let Some(ts) = by_key.get(o.2.as_str()) {
                    for &ti in ts {
                        pairs.push((oi, ti));
                    }
                }
            }
        } else {
            let total = obs_pts.len().saturating_mul(tgt_pts.len());
            if total > MAX_SIGHT_LINES {
                return Err(ToolError::Execution(format!(
                    "{} observers x {} target points would produce {total} sight lines \
                     (limit {MAX_SIGHT_LINES}); supply 'join_field' or coarsen 'sample_distance'",
                    obs_pts.len(),
                    tgt_pts.len()
                )));
            }
            for oi in 0..obs_pts.len() {
                for ti in 0..tgt_pts.len() {
                    pairs.push((oi, ti));
                }
            }
        }
        ctx.progress
            .info(&format!("constructing {} sight line(s)", pairs.len()));

        let mut out = Layer::new("sight_lines").with_geom_type(GeometryType::LineString);
        if let Some(epsg) = observers.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("OBSERVER_ID", FieldType::Integer));
        out.add_field(FieldDef::new("TARGET_ID", FieldType::Integer));
        out.add_field(FieldDef::new("DISTANCE", FieldType::Float));
        if want_dir {
            out.add_field(FieldDef::new("AZIMUTH", FieldType::Float));
            out.add_field(FieldDef::new("VERT_ANGLE", FieldType::Float));
        }

        for &(oi, ti) in &pairs {
            let (o_fid, o, _) = &obs_pts[oi];
            let (t_fid, t, _) = &tgt_pts[ti];
            let (dx, dy, dz) = (t[0] - o[0], t[1] - o[1], t[2] - o[2]);
            let horiz = dx.hypot(dy);
            let dist = match method {
                DistanceMethod::TwoD => horiz,
                DistanceMethod::ThreeD => (horiz * horiz + dz * dz).sqrt(),
            };
            let mut attrs = vec![
                ("OBSERVER_ID", FieldValue::Integer(*o_fid as i64)),
                ("TARGET_ID", FieldValue::Integer(*t_fid as i64)),
                ("DISTANCE", FieldValue::Float(dist)),
            ];
            if want_dir {
                // Compass azimuth: 0 = north, increasing clockwise.
                let az = dx.atan2(dy).to_degrees().rem_euclid(360.0);
                let vert = if horiz > 0.0 {
                    dz.atan2(horiz).to_degrees()
                } else if dz == 0.0 {
                    0.0
                } else {
                    90.0_f64.copysign(dz)
                };
                attrs.push(("AZIMUTH", FieldValue::Float(az)));
                attrs.push(("VERT_ANGLE", FieldValue::Float(vert)));
            }
            let geom = Geometry::line_string(vec![
                Coord::xyz(o[0], o[1], o[2]),
                Coord::xyz(t[0], t[1], t[2]),
            ]);
            out.add_feature(Some(geom), &attrs)
                .map_err(|e| ToolError::Execution(format!("failed adding sight line: {e}")))?;
        }

        let n_pairs = pairs.len();
        let n_obs = obs_pts.len();
        let n_tgt = tgt_pts.len();
        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("sight_line_count".to_string(), json!(n_pairs));
        outputs.insert("observer_count".to_string(), json!(n_obs));
        outputs.insert("target_point_count".to_string(), json!(n_tgt));
        Ok(ToolRunResult { outputs })
    }
}

/// First representative point of a geometry as `[x, y, z]` (z defaults to 0).
fn first_point(g: &Geometry) -> Option<[f64; 3]> {
    match g {
        Geometry::Point(c) => Some([c.x, c.y, c.z.unwrap_or(0.0)]),
        Geometry::MultiPoint(cs) => cs.first().map(|c| [c.x, c.y, c.z.unwrap_or(0.0)]),
        other => {
            let coords = other.all_coords();
            if coords.is_empty() {
                return None;
            }
            let n = coords.len() as f64;
            Some([
                coords.iter().map(|c| c.x).sum::<f64>() / n,
                coords.iter().map(|c| c.y).sum::<f64>() / n,
                coords.iter().map(|c| c.z.unwrap_or(0.0)).sum::<f64>() / n,
            ])
        }
    }
}

/// Target points for one geometry: points as-is, linear/areal boundaries
/// sampled every `spacing` map units (or at every vertex when `spacing` is
/// `None`).
fn target_points(g: &Geometry, spacing: Option<f64>) -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    match g {
        Geometry::Point(c) => out.push([c.x, c.y, c.z.unwrap_or(0.0)]),
        Geometry::MultiPoint(cs) => out.extend(cs.iter().map(|c| [c.x, c.y, c.z.unwrap_or(0.0)])),
        Geometry::LineString(cs) => sample_path(cs, spacing, &mut out),
        Geometry::MultiLineString(ls) => {
            for l in ls {
                sample_path(l, spacing, &mut out);
            }
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            sample_ring(&exterior.0, spacing, &mut out);
            for r in interiors {
                sample_ring(&r.0, spacing, &mut out);
            }
        }
        Geometry::MultiPolygon(ps) => {
            for (e, hs) in ps {
                sample_ring(&e.0, spacing, &mut out);
                for r in hs {
                    sample_ring(&r.0, spacing, &mut out);
                }
            }
        }
        Geometry::GeometryCollection(gs) => {
            for sub in gs {
                out.extend(target_points(sub, spacing));
            }
        }
    }
    out
}

/// Samples a closed ring, re-appending the first vertex so the closing segment
/// is walked too, then dropping the sample that lands back on the start.
fn sample_ring(coords: &[Coord], spacing: Option<f64>, out: &mut Vec<[f64; 3]>) {
    if coords.is_empty() {
        return;
    }
    let mut closed: Vec<Coord> = coords.to_vec();
    if closed.first().map(|c| (c.x, c.y)) != closed.last().map(|c| (c.x, c.y)) {
        closed.push(closed[0].clone());
    }
    let start = out.len();
    sample_path(&closed, spacing, out);
    // A ring whose perimeter is an exact multiple of the spacing emits a final
    // sample coincident with the first; keep the fan free of that duplicate.
    if out.len() - start >= 2 {
        let (first, last) = (out[start], out[out.len() - 1]);
        if (first[0] - last[0]).hypot(first[1] - last[1]) <= 1e-9 {
            out.pop();
        }
    }
}

/// Walks a polyline emitting a point every `spacing` units (Z interpolated), or
/// one point per vertex when `spacing` is `None`.
///
/// Sampling is driven by cumulative distance along the whole path rather than
/// a per-segment remainder, so a sample landing exactly on a vertex is emitted
/// once instead of being lost between the two segments that share it.
fn sample_path(coords: &[Coord], spacing: Option<f64>, out: &mut Vec<[f64; 3]>) {
    if coords.is_empty() {
        return;
    }
    let Some(step) = spacing else {
        out.extend(coords.iter().map(|c| [c.x, c.y, c.z.unwrap_or(0.0)]));
        return;
    };
    let first = &coords[0];
    out.push([first.x, first.y, first.z.unwrap_or(0.0)]);

    let mut acc = 0.0_f64; // path distance at the current segment's start
    let mut next = step; // path distance of the next sample
    for w in coords.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        let seg = (b.x - a.x).hypot(b.y - a.y);
        if seg <= 0.0 {
            continue;
        }
        let (az, bz) = (a.z.unwrap_or(0.0), b.z.unwrap_or(0.0));
        // Tolerance so a sample sitting on a vertex is not skipped by rounding.
        while next <= acc + seg + 1e-9 {
            let t = ((next - acc) / seg).clamp(0.0, 1.0);
            out.push([
                a.x + t * (b.x - a.x),
                a.y + t * (b.y - a.y),
                az + t * (bz - az),
            ]);
            next += step;
        }
        acc += seg;
    }
}

fn optional_field(
    layer: &Layer,
    name: Option<&str>,
    param: &str,
) -> Result<Option<usize>, ToolError> {
    match name {
        None => Ok(None),
        Some(n) => layer.schema.field_index(n).map(Some).ok_or_else(|| {
            ToolError::Validation(format!("{param} '{n}' not found in the layer"))
        }),
    }
}

fn key_of(v: Option<&FieldValue>) -> String {
    match v {
        None | Some(FieldValue::Null) => "NULL".to_string(),
        Some(FieldValue::Integer(i)) => i.to_string(),
        Some(FieldValue::Float(f)) => format!("{f}"),
        Some(FieldValue::Text(s)) => s.clone(),
        Some(FieldValue::Boolean(b)) => b.to_string(),
        Some(FieldValue::Date(s)) | Some(FieldValue::DateTime(s)) => s.clone(),
        Some(FieldValue::Blob(b)) => format!("blob[{}]", b.len()),
    }
}

// ── Params ──────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum DistanceMethod {
    TwoD,
    ThreeD,
}

fn parse_distance_method(args: &ToolArgs) -> Result<DistanceMethod, ToolError> {
    match args
        .get("distance_method")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        None | Some("") | Some("2d") => Ok(DistanceMethod::TwoD),
        Some("3d") => Ok(DistanceMethod::ThreeD),
        Some(o) => Err(ToolError::Validation(format!(
            "'distance_method' must be '2d' or '3d', got '{o}'"
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

    fn point_layer(pts: &[(f64, f64, f64)], keys: &[&str], hfield: bool) -> String {
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        if !keys.is_empty() {
            l.add_field(FieldDef::new("grp", FieldType::Text));
        }
        if hfield {
            l.add_field(FieldDef::new("ht", FieldType::Float));
        }
        for (i, (x, y, h)) in pts.iter().enumerate() {
            let mut a: Vec<(&str, FieldValue)> = Vec::new();
            if !keys.is_empty() {
                a.push(("grp", FieldValue::Text(keys[i % keys.len()].to_string())));
            }
            if hfield {
                a.push(("ht", FieldValue::Float(*h)));
            }
            l.add_feature(Some(Geometry::point(*x, *y)), &a).unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ConstructSightLinesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    #[test]
    fn cross_product_pairs_every_observer_with_every_target() {
        let obs = point_layer(&[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)], &[], false);
        let tgt = point_layer(&[(0.0, 10.0, 0.0), (5.0, 10.0, 0.0), (9.0, 9.0, 0.0)], &[], false);
        let (out, layer) = run(json!({ "observers": obs, "targets": tgt }));
        assert_eq!(out.outputs["sight_line_count"], json!(6));
        assert_eq!(layer.features.len(), 6);
        // Every emitted feature is a two-vertex line.
        for f in layer.iter() {
            match f.geometry.as_ref().unwrap() {
                Geometry::LineString(cs) => assert_eq!(cs.len(), 2),
                other => panic!("unexpected geometry {other:?}"),
            }
        }
    }

    #[test]
    fn join_field_restricts_pairing() {
        let obs = point_layer(&[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)], &["a", "b"], false);
        let tgt = point_layer(&[(0.0, 10.0, 0.0), (5.0, 10.0, 0.0)], &["a", "b"], false);
        let (out, _l) = run(json!({
            "observers": obs, "targets": tgt, "join_field": "grp"
        }));
        // 2 instead of 4: each observer sees only its own group's target.
        assert_eq!(out.outputs["sight_line_count"], json!(2));
    }

    #[test]
    fn distance_is_2d_by_default_and_3d_on_request() {
        // 3-4-5: horizontal 4, height 3 -> 3D distance 5.
        let obs = point_layer(&[(0.0, 0.0, 0.0)], &[], true);
        let tgt = point_layer(&[(4.0, 0.0, 3.0)], &[], true);
        let (_o, flat) = run(json!({
            "observers": obs.clone(), "targets": tgt.clone(),
            "target_height_field": "ht"
        }));
        let (_o, solid) = run(json!({
            "observers": obs, "targets": tgt,
            "target_height_field": "ht", "distance_method": "3d"
        }));
        let d = |l: &Layer| {
            let i = l.schema.field_index("DISTANCE").unwrap();
            l.features[0].attributes[i].as_f64().unwrap()
        };
        assert!((d(&flat) - 4.0).abs() < 1e-9);
        assert!((d(&solid) - 5.0).abs() < 1e-9);
    }

    #[test]
    fn direction_fields_use_compass_azimuth() {
        // Target due east -> azimuth 90; due north -> 0.
        let obs = point_layer(&[(0.0, 0.0, 0.0)], &[], false);
        let tgt = point_layer(&[(10.0, 0.0, 0.0), (0.0, 10.0, 0.0)], &[], false);
        let (_o, layer) = run(json!({
            "observers": obs, "targets": tgt, "output_direction": true
        }));
        let a = layer.schema.field_index("AZIMUTH").unwrap();
        let az: Vec<f64> = layer
            .iter()
            .map(|f| f.attributes[a].as_f64().unwrap())
            .collect();
        assert!((az[0] - 90.0).abs() < 1e-9, "east azimuth was {}", az[0]);
        assert!(az[1].abs() < 1e-9, "north azimuth was {}", az[1]);
    }

    #[test]
    fn areal_target_is_sampled_into_a_fan() {
        let obs = point_layer(&[(-10.0, 5.0, 0.0)], &[], false);
        let mut l = Layer::new("t")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        // 10x10 square: perimeter 40, sampled every 5 -> 8 boundary points.
        l.add_feature(
            Some(Geometry::polygon(
                vec![
                    Coord::xy(0.0, 0.0),
                    Coord::xy(10.0, 0.0),
                    Coord::xy(10.0, 10.0),
                    Coord::xy(0.0, 10.0),
                    Coord::xy(0.0, 0.0),
                ],
                vec![],
            )),
            &[],
        )
        .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let tgt = wbvector::memory_store::make_vector_memory_path(&id);
        let (out, _l) = run(json!({
            "observers": obs, "targets": tgt, "sample_distance": 5.0
        }));
        assert_eq!(out.outputs["target_point_count"], json!(8));
        assert_eq!(out.outputs["sight_line_count"], json!(8));
    }

    #[test]
    fn height_fields_offset_the_line_endpoints() {
        let obs = point_layer(&[(0.0, 0.0, 1.7)], &[], true);
        let tgt = point_layer(&[(10.0, 0.0, 30.0)], &[], true);
        let (_o, layer) = run(json!({
            "observers": obs, "targets": tgt,
            "observer_height_field": "ht", "target_height_field": "ht"
        }));
        match layer.features[0].geometry.as_ref().unwrap() {
            Geometry::LineString(cs) => {
                assert_eq!(cs[0].z, Some(1.7));
                assert_eq!(cs[1].z, Some(30.0));
            }
            other => panic!("unexpected geometry {other:?}"),
        }
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ConstructSightLinesTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "observers": "o.shp" })).is_err());
        assert!(bad(json!({ "observers": "o.shp", "targets": "t.shp", "sample_distance": 0 })).is_err());
        assert!(bad(json!({ "observers": "o.shp", "targets": "t.shp", "distance_method": "4d" })).is_err());
        assert!(bad(json!({ "observers": "o.shp", "targets": "t.shp" })).is_ok());
    }
}
