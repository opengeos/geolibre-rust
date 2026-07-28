//! GeoLibre tool: per-feature central meridian and standard parallels.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Calculate Central Meridian And
//! Parallels* (Cartography). It completes the repo's map-series projection
//! helpers: `calculate_utm_zone` assigns a UTM zone per feature and
//! `calculate_grid_convergence_angle` gives the grid-vs-true-north rotation,
//! but neither produces the *conic* projection parameters an atlas page needs
//! once its extent is wider than a UTM zone or sits at high latitude, where UTM
//! distortion becomes unacceptable.
//!
//! `grid_index_features` and `strip_map_index_features` (both shipped) generate
//! the page grid; this tool populates each page's projection parameters:
//!
//! * **central meridian** — the midpoint of the feature's longitude range
//! * **standard parallel 1 / 2** — inset from the south and north edges by
//!   `standard_offset` of the latitude range (Esri's 1/6 rule generalized)
//!
//! Features that straddle the antimeridian are handled by unwrapping longitudes
//! before taking the midpoint — averaging raw values there would place the
//! central meridian on the opposite side of the globe. The result is normalized
//! back into [-180, 180].

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CalculateCentralMeridianAndParallelsTool;

impl Tool for CalculateCentralMeridianAndParallelsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "calculate_central_meridian_and_parallels",
            display_name: "Calculate Central Meridian And Parallels",
            summary: "Compute a per-feature central meridian and two standard parallels from each feature's extent, for map-series and atlas pages, like ArcGIS Calculate Central Meridian And Parallels.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input features in geographic coordinates (longitude/latitude), typically a grid_index_features page grid.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "field",
                    description: "Name of the central-meridian field to write (default 'central_meridian').",
                    required: false,
                },
                ToolParamSpec {
                    name: "standard_offset",
                    description: "Fraction of the latitude range used to inset each standard parallel from the feature's edge (default 0.1667, Esri's 1/6 rule). Must be in [0, 0.5).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        if let Some(off) = parse_optional_f64(args, "standard_offset")? {
            if !off.is_finite() || !(0.0..0.5).contains(&off) {
                return Err(ToolError::Validation(
                    "'standard_offset' must be in [0, 0.5)".to_string(),
                ));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let field = parse_optional_str(args, "field")?.unwrap_or("central_meridian");
        let offset = parse_optional_f64(args, "standard_offset")?.unwrap_or(1.0 / 6.0);

        let layer = load_input_layer(input)?;
        if layer
            .crs_epsg()
            .is_some_and(|epsg| epsg != 4326 && epsg != 4269)
        {
            ctx.progress.info(
                "warning: input CRS is not geographic; coordinates will still be treated as longitude/latitude",
            );
        }

        for name in [field, "standard_parallel_1", "standard_parallel_2"] {
            if layer.schema.field_index(name).is_some() {
                return Err(ToolError::Validation(format!(
                    "output field '{name}' conflicts with another field"
                )));
            }
        }
        if field == "standard_parallel_1" || field == "standard_parallel_2" {
            return Err(ToolError::Validation(format!(
                "output field '{field}' conflicts with a standard-parallel field"
            )));
        }

        let mut out = Layer::new(layer.name.clone());
        out.crs = layer.crs.clone();
        out.geom_type = layer.geom_type;
        for fd in layer.schema.fields().iter() {
            out.add_field(fd.clone());
        }
        out.add_field(FieldDef::new(field, FieldType::Float));
        out.add_field(FieldDef::new("standard_parallel_1", FieldType::Float));
        out.add_field(FieldDef::new("standard_parallel_2", FieldType::Float));

        ctx.progress.info(&format!(
            "computing parameters for {} feature(s)",
            layer.len()
        ));

        let mut computed = 0usize;
        let mut antimeridian = 0usize;

        for (fi, feat) in layer.features.iter().enumerate() {
            let params = feat
                .geometry
                .as_ref()
                .and_then(|g| extent_params(g, offset));

            let mut fields: Vec<(String, FieldValue)> = layer
                .schema
                .fields()
                .iter()
                .enumerate()
                .map(|(i, fd)| (fd.name.clone(), feat.attributes[i].clone()))
                .collect();

            match params {
                Some(p) => {
                    if p.crossed_antimeridian {
                        antimeridian += 1;
                    }
                    computed += 1;
                    fields.push((field.to_string(), FieldValue::Float(p.central_meridian)));
                    fields.push((
                        "standard_parallel_1".to_string(),
                        FieldValue::Float(p.standard_parallel_1),
                    ));
                    fields.push((
                        "standard_parallel_2".to_string(),
                        FieldValue::Float(p.standard_parallel_2),
                    ));
                }
                None => {
                    // Empty or non-coordinate geometry: emit nulls rather than a
                    // misleading 0 that would read as the Greenwich meridian.
                    fields.push((field.to_string(), FieldValue::Null));
                    fields.push(("standard_parallel_1".to_string(), FieldValue::Null));
                    fields.push(("standard_parallel_2".to_string(), FieldValue::Null));
                }
            }

            let refs: Vec<(&str, FieldValue)> = fields
                .iter()
                .map(|(k, v)| (k.as_str(), v.clone()))
                .collect();
            out.add_feature(feat.geometry.clone(), &refs)
                .map_err(|e| ToolError::Execution(format!("failed writing feature: {e}")))?;
            ctx.progress
                .progress((fi as f64 + 1.0) / layer.len().max(1) as f64);
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(layer.len()));
        outputs.insert("features_computed".to_string(), json!(computed));
        outputs.insert("antimeridian_features".to_string(), json!(antimeridian));
        Ok(ToolRunResult { outputs })
    }
}

struct Params {
    central_meridian: f64,
    standard_parallel_1: f64,
    standard_parallel_2: f64,
    crossed_antimeridian: bool,
}

/// Derives the projection parameters from a geometry's geographic extent.
fn extent_params(geom: &Geometry, offset: f64) -> Option<Params> {
    let mut lons: Vec<f64> = Vec::new();
    let mut lat_min = f64::INFINITY;
    let mut lat_max = f64::NEG_INFINITY;
    collect(geom, &mut lons, &mut lat_min, &mut lat_max);
    // Drop non-finite longitudes before the min/max folds: a single NaN would
    // poison them and write NaN into the output, and an infinite value would
    // reach normalize_lon.
    lons.retain(|l| l.is_finite());
    if lons.is_empty() || !lat_min.is_finite() || !lat_max.is_finite() {
        return None;
    }

    let (central, crossed) = central_meridian(&lons);
    let span = lat_max - lat_min;
    Some(Params {
        central_meridian: central,
        standard_parallel_1: lat_min + span * offset,
        standard_parallel_2: lat_max - span * offset,
        crossed_antimeridian: crossed,
    })
}

/// Midpoint of a longitude set, unwrapping across the antimeridian.
///
/// A feature spanning 179°E to -179°E has a raw mean of 0° — the wrong side of
/// the planet.
///
/// The crossing test is `raw_span > 180`, which is the actual geometric
/// signature: no non-crossing extent can span more than half the globe, and any
/// crossing one necessarily appears to. Only then is the +360 unwrap applied.
///
/// An earlier version instead compared the raw span against the unwrapped span
/// and took "unwrapped is narrower" as the signal. That is mathematically
/// equivalent but numerically unsafe: for an extent lying wholly in the western
/// hemisphere the two spans are *equal*, and float noise in the two subtraction
/// chains made the shifted one infinitesimally smaller, flagging ordinary
/// extents (e.g. the US South census region, -106.6..-75.1) as crossing.
fn central_meridian(lons: &[f64]) -> (f64, bool) {
    let raw_min = lons.iter().cloned().fold(f64::INFINITY, f64::min);
    let raw_max = lons.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let raw_span = raw_max - raw_min;

    if raw_span <= 180.0 {
        return (normalize_lon((raw_min + raw_max) / 2.0), false);
    }
    if raw_span >= 360.0 - f64::EPSILON {
        return (normalize_lon((raw_min + raw_max) / 2.0), false);
    }

    let shifted: Vec<f64> = lons
        .iter()
        .map(|&l| if l < 0.0 { l + 360.0 } else { l })
        .collect();
    let sh_min = shifted.iter().cloned().fold(f64::INFINITY, f64::min);
    let sh_max = shifted.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    // If unwrapping does not actually narrow the extent, the feature genuinely
    // spans more than half the globe (a whole-world layer); keep the raw answer.
    if sh_max - sh_min < raw_span {
        (normalize_lon((sh_min + sh_max) / 2.0), true)
    } else {
        (normalize_lon((raw_min + raw_max) / 2.0), false)
    }
}

/// Wraps a longitude into [-180, 180].
///
/// Uses `rem_euclid` rather than a subtract-until-in-range loop: the loop form
/// never terminates for a non-finite input (`inf - 360.0 == inf`), which would
/// hang the caller. Non-finite input is returned unchanged so the caller's own
/// finiteness checks stay authoritative.
fn normalize_lon(lon: f64) -> f64 {
    if !lon.is_finite() {
        return lon;
    }
    let wrapped = (lon + 180.0).rem_euclid(360.0) - 180.0;
    // rem_euclid maps exactly +180 to -180; keep the positive representation.
    if wrapped == -180.0 && lon > 0.0 {
        180.0
    } else {
        wrapped
    }
}

fn collect(geom: &Geometry, lons: &mut Vec<f64>, lat_min: &mut f64, lat_max: &mut f64) {
    let mut push = |c: &Coord| {
        if !c.x.is_finite() || !c.y.is_finite() {
            return;
        }
        lons.push(c.x);
        *lat_min = lat_min.min(c.y);
        *lat_max = lat_max.max(c.y);
    };
    match geom {
        Geometry::Point(c) => push(c),
        Geometry::MultiPoint(cs) | Geometry::LineString(cs) => cs.iter().for_each(push),
        Geometry::MultiLineString(parts) => parts.iter().flatten().for_each(push),
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            exterior.coords().iter().for_each(&mut push);
            for r in interiors {
                r.coords().iter().for_each(&mut push);
            }
        }
        Geometry::MultiPolygon(parts) => {
            for (e, hs) in parts {
                e.coords().iter().for_each(&mut push);
                for r in hs {
                    r.coords().iter().for_each(&mut push);
                }
            }
        }
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                collect(g, lons, lat_min, lat_max);
            }
        }
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

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Geometry {
        Geometry::polygon(
            vec![
                Coord::xy(x0, y0),
                Coord::xy(x1, y0),
                Coord::xy(x1, y1),
                Coord::xy(x0, y1),
            ],
            vec![],
        )
    }

    fn layer_of(geoms: Vec<Geometry>) -> String {
        let mut l = Layer::new("pages")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        for g in geoms {
            l.add_feature(Some(g), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CalculateCentralMeridianAndParallelsTool
            .run(&args, &ctx())
            .unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn num(layer: &Layer, row: usize, field: &str) -> f64 {
        let i = layer.schema.field_index(field).unwrap();
        layer.features[row].attributes[i].as_f64().unwrap()
    }

    /// Central meridian is the longitude midpoint; parallels inset by 1/6.
    #[test]
    fn computes_meridian_and_parallels() {
        // Extent 10E..20E, 30N..60N -> CM 15, span 30, inset 5 -> 35N / 55N.
        let (_, layer) = run(json!({ "input": layer_of(vec![rect(10.0, 30.0, 20.0, 60.0)]) }));
        assert!((num(&layer, 0, "central_meridian") - 15.0).abs() < 1e-9);
        assert!((num(&layer, 0, "standard_parallel_1") - 35.0).abs() < 1e-9);
        assert!((num(&layer, 0, "standard_parallel_2") - 55.0).abs() < 1e-9);
    }

    /// THE edge case: a page straddling the antimeridian must not land on 0.
    #[test]
    fn handles_antimeridian_crossing() {
        // 175E .. -175E. Naive mean is 0 (wrong hemisphere); correct is 180.
        let (out, layer) = run(json!({
            "input": layer_of(vec![rect(175.0, 0.0, -175.0, 10.0)])
        }));
        let cm = num(&layer, 0, "central_meridian");
        assert!(
            cm.abs() > 179.0,
            "expected ~180 across the antimeridian, got {cm}"
        );
        assert_eq!(out.outputs["antimeridian_features"], json!(1));
    }

    /// A normal extent must NOT be misdetected as antimeridian-crossing.
    #[test]
    fn does_not_false_positive_on_normal_extent() {
        let (out, layer) = run(json!({
            "input": layer_of(vec![rect(-120.0, 30.0, -100.0, 50.0)])
        }));
        assert!((num(&layer, 0, "central_meridian") - -110.0).abs() < 1e-9);
        assert_eq!(out.outputs["antimeridian_features"], json!(0));
    }

    /// Regression: an extent lying wholly in the western hemisphere has raw and
    /// unwrapped spans that are mathematically EQUAL, so a naive
    /// "unwrapped is narrower" test flips on float noise. Caught on the real US
    /// South census region (-106.6..-75.1), which was wrongly flagged.
    #[test]
    fn western_hemisphere_extent_is_not_flagged_as_crossing() {
        for (lo, hi) in [
            (-106.645646_f64, -75.045448_f64), // US South census region
            (-104.057698, -80.518798),         // US Midwest
            (-80.519891, -66.949895),          // US Northeast
            (-124.733253, -114.039403),        // a single western state
        ] {
            let (_, crossed) = central_meridian(&[lo, hi, (lo + hi) / 2.0]);
            assert!(
                !crossed,
                "extent {lo}..{hi} lies in one hemisphere and must not be flagged"
            );
        }
    }

    /// A genuinely crossing extent (Alaska including the Aleutians) still is.
    #[test]
    fn genuine_crossing_extent_is_flagged() {
        let (cm, crossed) = central_meridian(&[-179.17, 179.77, -130.0]);
        assert!(
            crossed,
            "an extent spanning the antimeridian must be flagged"
        );
        assert!(
            cm.abs() > 90.0,
            "central meridian should sit near 180, got {cm}"
        );
    }

    #[test]
    fn whole_world_extent_is_not_flagged_as_crossing() {
        let (cm, crossed) = central_meridian(&[-180.0, 0.0, 180.0]);
        assert!(!crossed);
        assert!(cm.abs() < 1e-9);
    }

    /// standard_offset = 0 puts the parallels on the edges.
    #[test]
    fn zero_offset_uses_edges() {
        let (_, layer) = run(json!({
            "input": layer_of(vec![rect(0.0, 20.0, 10.0, 40.0)]),
            "standard_offset": 0
        }));
        assert!((num(&layer, 0, "standard_parallel_1") - 20.0).abs() < 1e-9);
        assert!((num(&layer, 0, "standard_parallel_2") - 40.0).abs() < 1e-9);
    }

    /// A custom field name is honoured.
    #[test]
    fn custom_field_name() {
        let (_, layer) = run(json!({
            "input": layer_of(vec![rect(0.0, 0.0, 10.0, 10.0)]),
            "field": "cm"
        }));
        assert!(layer.schema.field_index("cm").is_some());
        assert!((num(&layer, 0, "cm") - 5.0).abs() < 1e-9);
    }

    /// Existing attributes survive.
    #[test]
    fn preserves_input_attributes() {
        let mut l = Layer::new("pages")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("page", FieldType::Integer));
        l.add_feature(
            Some(rect(0.0, 0.0, 10.0, 10.0)),
            &[("page", FieldValue::Integer(7))],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        let (_, layer) = run(json!({ "input": memory_store::make_vector_memory_path(&id) }));
        assert_eq!(num(&layer, 0, "page"), 7.0);
    }

    #[test]
    fn rejects_bad_parameters() {
        let p = layer_of(vec![rect(0.0, 0.0, 1.0, 1.0)]);
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CalculateCentralMeridianAndParallelsTool
                .validate(&args)
                .is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({ "input": p, "standard_offset": 0.5 })));
        assert!(bad(json!({ "input": p, "standard_offset": -0.1 })));
        assert!(bad(json!({ "input": p, "standard_offset": "abc" })));
    }
}
