//! GeoLibre tool: build antenna coverage sectors from cell-tower points.
//!
//! Pure-Rust counterpart of the geometry-generation half of ArcGIS Pro's
//! *Cell Site Records To Feature Class* / *Generate Sector Lines* (Crime
//! Analysis & Safety toolbox). Neither the GeoLibre catalog nor the bundled
//! whitebox suite builds antenna coverage geometry, so this opens a new domain.
//!
//! Each input point is an antenna at a tower. From its **azimuth** (compass
//! bearing of the main lobe, 0 = north, clockwise), **beamwidth** (angular
//! width of the lobe in degrees), and **radius** (coverage distance in CRS
//! units) the tool emits either:
//!   * a *wedge* polygon — the apex at the tower plus an arc sampled at
//!     `segments` steps between `azimuth ± beamwidth/2`, or
//!   * a *line* — the bisector sector line from the tower out to `radius`
//!     along the azimuth.
//!
//! Per-antenna azimuth/beamwidth/radius are read from attribute fields when the
//! `*_field` params are given (falling back to the scalar default when a
//! feature's value is missing/invalid); otherwise the scalar defaults apply to
//! every tower. A beamwidth of ≥ 360 yields an omnidirectional full circle.
//!
//! All original tower attributes are preserved, and each output carries the
//! resolved `azimuth`, `beamwidth`, `radius`, and (for wedges) the polygon
//! `area` so callers can validate coverage — the geometry is pure trigonometry.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CellRecordsToSectorsTool;

impl Tool for CellRecordsToSectorsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "cell_records_to_sectors",
            display_name: "Cell Records To Sectors",
            summary: "Build antenna coverage geometry from cell-tower points: wedge polygons or bisector sector lines from per-tower azimuth, beamwidth, and radius — like ArcGIS Cell Site Records To Feature Class / Generate Sector Lines.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer of antennas/towers.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_type",
                    description: "Geometry to emit: 'wedge' (coverage polygons, default) or 'line' (bisector sector lines).",
                    required: false,
                },
                ToolParamSpec {
                    name: "azimuth_field",
                    description: "Field holding each antenna's azimuth (compass degrees, 0 = north, clockwise). Falls back to 'azimuth' when missing.",
                    required: false,
                },
                ToolParamSpec {
                    name: "beamwidth_field",
                    description: "Field holding each antenna's beamwidth in degrees. Falls back to 'beamwidth' when missing.",
                    required: false,
                },
                ToolParamSpec {
                    name: "radius_field",
                    description: "Field holding each antenna's coverage radius in CRS units. Falls back to 'radius' when missing.",
                    required: false,
                },
                ToolParamSpec {
                    name: "azimuth",
                    description: "Default azimuth in compass degrees when no field is given / a value is missing. Default 0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "beamwidth",
                    description: "Default beamwidth in degrees (0 < bw ≤ 360; ≥ 360 = omnidirectional). Default 65.",
                    required: false,
                },
                ToolParamSpec {
                    name: "radius",
                    description: "Default coverage radius in CRS units (must be positive). Default 1000.",
                    required: false,
                },
                ToolParamSpec {
                    name: "segments",
                    description: "Number of straight segments used to approximate each sector arc (≥ 1). Default 16.",
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

        // Resolve field indices for any *_field params (error if named but absent).
        let az_idx = resolve_field(&layer, prm.azimuth_field.as_deref(), "azimuth_field")?;
        let bw_idx = resolve_field(&layer, prm.beamwidth_field.as_deref(), "beamwidth_field")?;
        let r_idx = resolve_field(&layer, prm.radius_field.as_deref(), "radius_field")?;

        // Build output schema: preserve all input fields, then append computed ones.
        let geom_type = match prm.output_type {
            OutputType::Wedge => GeometryType::Polygon,
            OutputType::Line => GeometryType::LineString,
        };
        let mut out = Layer::new("sectors").with_geom_type(geom_type);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for fd in layer.schema.fields() {
            out.add_field(fd.clone());
        }
        // Computed fields (renamed if they would collide with an input field).
        let az_out = unique_name(&layer, "azimuth");
        let bw_out = unique_name(&layer, "beamwidth");
        let r_out = unique_name(&layer, "radius");
        out.add_field(FieldDef::new(&az_out, FieldType::Float));
        out.add_field(FieldDef::new(&bw_out, FieldType::Float));
        out.add_field(FieldDef::new(&r_out, FieldType::Float));
        let area_out = if matches!(prm.output_type, OutputType::Wedge) {
            let n = unique_name(&layer, "area");
            out.add_field(FieldDef::new(&n, FieldType::Float));
            Some(n)
        } else {
            None
        };

        let mut sector_count = 0usize;
        let mut skipped = 0usize;
        for feature in layer.features.iter() {
            let Some(geom) = feature.geometry.as_ref() else {
                skipped += 1;
                continue;
            };
            let apexes = point_apexes(geom);
            if apexes.is_empty() {
                skipped += 1;
                continue;
            }

            // Resolve per-feature parameters (field value, else scalar default).
            let az = field_or_default(feature, az_idx, prm.azimuth);
            let bw = field_or_default(feature, bw_idx, prm.beamwidth);
            let r = field_or_default(feature, r_idx, prm.radius);
            let valid = bw > 0.0 && bw.is_finite() && r > 0.0 && r.is_finite() && az.is_finite();
            if !valid {
                skipped += 1;
                continue;
            }
            let bw = bw.min(360.0);

            for apex in apexes {
                let geometry = match prm.output_type {
                    OutputType::Wedge => wedge_polygon(apex, az, bw, r, prm.segments),
                    OutputType::Line => bisector_line(apex, az, r),
                };
                // Copy preserved attributes, then set computed ones.
                let mut attrs: Vec<(&str, FieldValue)> = Vec::new();
                for (i, fd) in layer.schema.fields().iter().enumerate() {
                    if let Some(v) = feature.attributes.get(i) {
                        attrs.push((fd.name.as_str(), v.clone()));
                    }
                }
                attrs.push((az_out.as_str(), FieldValue::Float(az)));
                attrs.push((bw_out.as_str(), FieldValue::Float(bw)));
                attrs.push((r_out.as_str(), FieldValue::Float(r)));
                if let Some(name) = &area_out {
                    let area = polygon_area(&geometry);
                    attrs.push((name.as_str(), FieldValue::Float(area)));
                }
                out.add_feature(Some(geometry), &attrs)
                    .map_err(|e| ToolError::Execution(format!("failed writing sector: {e}")))?;
                sector_count += 1;
            }
        }

        ctx.progress.info(&format!(
            "generated {sector_count} sector(s); skipped {skipped} feature(s)"
        ));

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("sector_count".to_string(), json!(sector_count));
        outputs.insert("skipped".to_string(), json!(skipped));
        Ok(ToolRunResult { outputs })
    }
}

// ── Geometry construction ──────────────────────────────────────────────────────

/// Extracts apex coordinates from a point or multipoint geometry.
fn point_apexes(geom: &Geometry) -> Vec<Coord> {
    match geom {
        Geometry::Point(c) => vec![c.clone()],
        Geometry::MultiPoint(cs) => cs.clone(),
        _ => Vec::new(),
    }
}

/// Unit direction vector for a compass bearing (0 = north/+y, clockwise).
fn dir(bearing_deg: f64) -> (f64, f64) {
    let r = bearing_deg.to_radians();
    (r.sin(), r.cos()) // (east/x, north/y)
}

/// Samples the sector arc from `az - bw/2` to `az + bw/2` at `segments` steps.
fn arc_points(apex: &Coord, az: f64, bw: f64, radius: f64, segments: usize) -> Vec<Coord> {
    let n = segments.max(1);
    let start = az - bw * 0.5;
    let mut pts = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let b = start + bw * (i as f64) / (n as f64);
        let (dx, dy) = dir(b);
        pts.push(Coord::xy(apex.x + radius * dx, apex.y + radius * dy));
    }
    pts
}

/// Builds a wedge coverage polygon (apex + arc), or a full circle when bw ≥ 360.
fn wedge_polygon(apex: Coord, az: f64, bw: f64, radius: f64, segments: usize) -> Geometry {
    // Ring stored without the closing duplicate vertex (wbvector convention).
    if bw >= 360.0 {
        // Omnidirectional: closed circle, no apex vertex.
        let n = (segments.max(3)).max(8);
        let mut ring = Vec::with_capacity(n);
        for i in 0..n {
            let b = 360.0 * (i as f64) / (n as f64);
            let (dx, dy) = dir(b);
            ring.push(Coord::xy(apex.x + radius * dx, apex.y + radius * dy));
        }
        return Geometry::polygon(ring, vec![]);
    }
    let mut ring = Vec::with_capacity(segments + 2);
    ring.push(apex.clone());
    ring.extend(arc_points(&apex, az, bw, radius, segments));
    Geometry::polygon(ring, vec![])
}

/// Builds the bisector sector line from the apex out to `radius` along azimuth.
fn bisector_line(apex: Coord, az: f64, radius: f64) -> Geometry {
    let (dx, dy) = dir(az);
    let tip = Coord::xy(apex.x + radius * dx, apex.y + radius * dy);
    Geometry::line_string(vec![apex, tip])
}

/// Shoelace area of a polygon's exterior ring (0 for non-polygons).
fn polygon_area(geom: &Geometry) -> f64 {
    match geom {
        Geometry::Polygon { exterior, .. } => exterior.signed_area().abs(),
        _ => 0.0,
    }
}

// ── Attribute helpers ──────────────────────────────────────────────────────────

/// Reads a numeric field value from a feature, else the scalar default.
fn field_or_default(feature: &wbvector::Feature, idx: Option<usize>, default: f64) -> f64 {
    match idx
        .and_then(|i| feature.attributes.get(i))
        .and_then(FieldValue::as_f64)
    {
        Some(v) if v.is_finite() => v,
        _ => default,
    }
}

/// Resolves a `*_field` name to a schema index, erroring if named but absent.
fn resolve_field(
    layer: &Layer,
    name: Option<&str>,
    param: &str,
) -> Result<Option<usize>, ToolError> {
    match name {
        None => Ok(None),
        Some(n) => layer.schema.field_index(n).map(Some).ok_or_else(|| {
            ToolError::Validation(format!(
                "parameter '{param}' names field '{n}' not in the input"
            ))
        }),
    }
}

/// Returns `base`, or `base_1`, `base_2`, … until it does not collide with an
/// existing input field, so preserved attributes are never overwritten.
fn unique_name(layer: &Layer, base: &str) -> String {
    if layer.schema.field_index(base).is_none() {
        return base.to_string();
    }
    let mut i = 1;
    loop {
        let candidate = format!("{base}_{i}");
        if layer.schema.field_index(&candidate).is_none() {
            return candidate;
        }
        i += 1;
    }
}

// ── Parameters ──────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq)]
enum OutputType {
    Wedge,
    Line,
}

struct Params {
    output_type: OutputType,
    azimuth_field: Option<String>,
    beamwidth_field: Option<String>,
    radius_field: Option<String>,
    azimuth: f64,
    beamwidth: f64,
    radius: f64,
    segments: usize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let output_type = match parse_optional_str(args, "output_type")? {
        None => OutputType::Wedge,
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "wedge" | "polygon" | "sector" => OutputType::Wedge,
            "line" | "bisector" => OutputType::Line,
            other => {
                return Err(ToolError::Validation(format!(
                    "parameter 'output_type' must be 'wedge' or 'line' (got '{other}')"
                )))
            }
        },
    };

    let azimuth = parse_optional_f64(args, "azimuth")?.unwrap_or(0.0);
    if !azimuth.is_finite() {
        return Err(ToolError::Validation(
            "'azimuth' must be a finite number".to_string(),
        ));
    }

    let beamwidth = parse_optional_f64(args, "beamwidth")?.unwrap_or(65.0);
    if !(beamwidth > 0.0 && beamwidth <= 360.0) {
        return Err(ToolError::Validation(
            "'beamwidth' must be in the range (0, 360]".to_string(),
        ));
    }

    let radius = parse_optional_f64(args, "radius")?.unwrap_or(1000.0);
    if !(radius > 0.0 && radius.is_finite()) {
        return Err(ToolError::Validation(
            "'radius' must be a positive number".to_string(),
        ));
    }

    let segments = match parse_optional_f64(args, "segments")? {
        None => 16,
        Some(v) => {
            if !(1.0..=4096.0).contains(&v) {
                return Err(ToolError::Validation(
                    "'segments' must be an integer in [1, 4096]".to_string(),
                ));
            }
            v as usize
        }
    };

    Ok(Params {
        output_type,
        azimuth_field: parse_optional_str(args, "azimuth_field")?.map(str::to_string),
        beamwidth_field: parse_optional_str(args, "beamwidth_field")?.map(str::to_string),
        radius_field: parse_optional_str(args, "radius_field")?.map(str::to_string),
        azimuth,
        beamwidth,
        radius,
        segments,
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;
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

    /// A test tower: (x, y, optional (azimuth, beamwidth, radius) attributes).
    type Tower = (f64, f64, Option<(f64, f64, f64)>);

    /// Builds a point layer with optional (azimuth, beamwidth, radius) fields.
    fn tower_layer(points: &[Tower]) -> String {
        let mut l = Layer::new("towers")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("site", FieldType::Text));
        l.add_field(FieldDef::new("az", FieldType::Float));
        l.add_field(FieldDef::new("bw", FieldType::Float));
        l.add_field(FieldDef::new("rad", FieldType::Float));
        for (i, (x, y, params)) in points.iter().enumerate() {
            let (az, bw, rad) = params.unwrap_or((0.0, 0.0, 0.0));
            l.add_feature(
                Some(Geometry::point(*x, *y)),
                &[
                    ("site", FieldValue::Text(format!("S{i}"))),
                    ("az", FieldValue::Float(az)),
                    ("bw", FieldValue::Float(bw)),
                    ("rad", FieldValue::Float(rad)),
                ],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CellRecordsToSectorsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// Wedge polygon area matches the N-segment circular-sector formula, and one
    /// wedge is produced per tower.
    #[test]
    fn wedge_area_matches_sector_formula() {
        let input = tower_layer(&[(0.0, 0.0, None), (500.0, 0.0, None)]);
        let (out, layer) = run(json!({
            "input": input, "beamwidth": 90.0, "radius": 100.0, "segments": 64,
        }));
        assert_eq!(out.outputs["sector_count"], json!(2));
        assert_eq!(layer.features.len(), 2);
        // N-segment polygon sector area: 0.5 r² N sin(θ/N).
        let theta = 90.0_f64.to_radians();
        let n = 64.0;
        let expected = 0.5 * 100.0 * 100.0 * n * (theta / n).sin();
        let aidx = layer.schema.field_index("area").unwrap();
        for f in layer.iter() {
            let a = f.attributes[aidx].as_f64().unwrap();
            assert!(
                (a - expected).abs() / expected < 1e-6,
                "area {a} != expected {expected}"
            );
        }
    }

    /// The wedge apex sits at the tower and the arc endpoints lie at `radius`,
    /// centred on the azimuth (north here).
    #[test]
    fn wedge_apex_and_arc_geometry() {
        let input = tower_layer(&[(10.0, 20.0, None)]);
        let (_o, layer) = run(json!({
            "input": input, "azimuth": 0.0, "beamwidth": 60.0, "radius": 50.0, "segments": 6,
        }));
        let Geometry::Polygon { exterior, .. } = layer.features[0].geometry.as_ref().unwrap()
        else {
            panic!("expected polygon");
        };
        let ring = exterior.coords();
        // First vertex is the apex (tower location).
        assert!((ring[0].x - 10.0).abs() < 1e-9 && (ring[0].y - 20.0).abs() < 1e-9);
        // Every arc vertex is at `radius` from the apex.
        for c in &ring[1..] {
            let d = (c.x - 10.0).hypot(c.y - 20.0);
            assert!((d - 50.0).abs() < 1e-6, "arc vertex not on radius: {d}");
        }
        // Middle arc vertex points due north (azimuth 0): x == apex x.
        let mid = &ring[1 + 3];
        assert!((mid.x - 10.0).abs() < 1e-6, "sector not centred on azimuth");
        assert!(mid.y > 20.0, "azimuth 0 should point north");
    }

    /// Bisector line output is a 2-vertex line of length `radius` along azimuth.
    #[test]
    fn line_output_is_bisector() {
        let input = tower_layer(&[(0.0, 0.0, None)]);
        let (out, layer) = run(json!({
            "input": input, "output_type": "line", "azimuth": 90.0, "radius": 30.0,
        }));
        assert_eq!(out.outputs["sector_count"], json!(1));
        let Geometry::LineString(cs) = layer.features[0].geometry.as_ref().unwrap() else {
            panic!("expected line");
        };
        assert_eq!(cs.len(), 2);
        // Azimuth 90 = east: tip at (+30, 0).
        assert!((cs[0].x).abs() < 1e-9 && (cs[0].y).abs() < 1e-9);
        assert!((cs[1].x - 30.0).abs() < 1e-6 && cs[1].y.abs() < 1e-6);
    }

    /// Omnidirectional beamwidth (360) yields a full circle whose area ≈ π r².
    #[test]
    fn omnidirectional_is_full_circle() {
        let input = tower_layer(&[(0.0, 0.0, None)]);
        let (_o, layer) = run(json!({
            "input": input, "beamwidth": 360.0, "radius": 100.0, "segments": 256,
        }));
        let aidx = layer.schema.field_index("area").unwrap();
        let a = layer.features[0].attributes[aidx].as_f64().unwrap();
        let circle = PI * 100.0 * 100.0;
        // Inscribed polygon slightly under-estimates; loose upper bound, tight lower.
        assert!(a < circle && a > circle * 0.999, "circle area off: {a}");
    }

    /// Per-tower field values drive each sector; attributes are preserved.
    #[test]
    fn field_driven_and_attributes_preserved() {
        let input = tower_layer(&[
            (0.0, 0.0, Some((0.0, 90.0, 10.0))),
            (100.0, 0.0, Some((0.0, 90.0, 20.0))),
        ]);
        let (_o, layer) = run(json!({
            "input": input,
            "azimuth_field": "az", "beamwidth_field": "bw", "radius_field": "rad",
            "segments": 128,
        }));
        // Original 'site' attribute survives.
        let sidx = layer.schema.field_index("site").unwrap();
        assert_eq!(layer.features[0].attributes[sidx].as_str(), Some("S0"));
        // Areas scale with radius²: second tower r=20 vs first r=10 → 4×.
        let aidx = layer.schema.field_index("area").unwrap();
        let a0 = layer.features[0].attributes[aidx].as_f64().unwrap();
        let a1 = layer.features[1].attributes[aidx].as_f64().unwrap();
        assert!(
            (a1 / a0 - 4.0).abs() < 1e-6,
            "radius scaling wrong: {a0} {a1}"
        );
    }

    /// Features whose resolved beamwidth/radius are invalid are skipped, not emitted.
    #[test]
    fn invalid_features_skipped() {
        // Tower params default to (0,0,0) → bw=0, r=0 which are invalid.
        let input = tower_layer(&[(0.0, 0.0, Some((0.0, 0.0, 0.0)))]);
        let (out, _layer) = run(json!({
            "input": input,
            "azimuth_field": "az", "beamwidth_field": "bw", "radius_field": "rad",
        }));
        assert_eq!(out.outputs["sector_count"], json!(0));
        assert_eq!(out.outputs["skipped"], json!(1));
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CellRecordsToSectorsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err()); // no input
        assert!(bad(json!({ "input": "a.geojson", "beamwidth": 0 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "beamwidth": 400 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "radius": -5 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "output_type": "blob" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "segments": 0 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "beamwidth": 65, "radius": 500 })).is_ok());
    }
}
