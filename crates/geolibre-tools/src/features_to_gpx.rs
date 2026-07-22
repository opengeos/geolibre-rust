//! GeoLibre tool: write point/line GIS features out as a GPX 1.1 file.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Features To GPX* (Conversion): the
//! export direction that pairs with GPX reading. Point features become GPX
//! waypoints (`<wpt>`); line features become GPX tracks (`<trk>` with one
//! `<trkseg>` per line part, exploded into ordered `<trkpt>` vertices).
//!
//! Chosen attribute fields are mapped onto the standard GPX metadata elements:
//! `name_field` -> `<name>`, `description_field` -> `<desc>`,
//! `z_field` -> `<ele>` (elevation, falling back to a geometry Z when the field
//! is absent), and `date_field` -> `<time>` (an ISO-8601 timestamp). Any field
//! left unset is simply omitted, producing schema-valid GPX.
//!
//! GPX coordinates are always lon/lat WGS84, so inputs are assumed EPSG:4326
//! (as the GPX reader emits) and coordinates are written straight through.
//! Polygon / other geometries are skipped and reported as a count. The GPX is
//! serialized with the `quick-xml` writer.

use std::collections::BTreeMap;

use quick_xml::events::{BytesDecl, BytesEnd, BytesStart, BytesText, Event};
use quick_xml::Writer;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldValue, Geometry, Layer, Schema};

use crate::vector_common::{ensure_parent_dir, load_input_layer, parse_optional_str};

pub struct FeaturesToGpxTool;

/// Field names chosen to map onto GPX metadata elements.
struct FieldMap {
    name: Option<String>,
    description: Option<String>,
    z: Option<String>,
    date: Option<String>,
}

/// Counts returned as run outputs.
#[derive(Default)]
struct Counts {
    waypoints: usize,
    tracks: usize,
    track_points: usize,
    skipped: usize,
}

impl Tool for FeaturesToGpxTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "features_to_gpx",
            display_name: "Features To GPX",
            summary: "Write point features as GPX waypoints and line features as GPX tracks, mapping attribute fields to name/description/elevation/time.",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Point or line vector layer to export (assumed WGS84 / EPSG:4326).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output .gpx file path.",
                    required: true,
                },
                ToolParamSpec {
                    name: "name_field",
                    description: "Attribute field mapped to the GPX <name> element (optional).",
                    required: false,
                },
                ToolParamSpec {
                    name: "description_field",
                    description: "Attribute field mapped to the GPX <desc> element (optional).",
                    required: false,
                },
                ToolParamSpec {
                    name: "z_field",
                    description: "Attribute field mapped to the GPX <ele> elevation (optional; falls back to a geometry Z when present).",
                    required: false,
                },
                ToolParamSpec {
                    name: "date_field",
                    description: "Attribute field mapped to the GPX <time> timestamp, an ISO-8601 string (optional).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if str_arg(args, "input").is_empty() {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        if str_arg(args, "output").is_empty() {
            return Err(ToolError::Validation(
                "missing required string parameter 'output'".to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = parse_optional_str(args, "input")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'input'".to_string())
        })?;
        let output = parse_optional_str(args, "output")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'output'".to_string())
        })?;
        if !output.to_ascii_lowercase().ends_with(".gpx") {
            return Err(ToolError::Validation(
                "parameter 'output' must be a .gpx file path".to_string(),
            ));
        }

        let fields = FieldMap {
            name: parse_optional_str(args, "name_field")?.map(str::to_string),
            description: parse_optional_str(args, "description_field")?.map(str::to_string),
            z: parse_optional_str(args, "z_field")?.map(str::to_string),
            date: parse_optional_str(args, "date_field")?.map(str::to_string),
        };

        let layer = load_input_layer(input)?;
        let (gpx, counts) = serialize_gpx(&layer, &fields)
            .map_err(|e| ToolError::Execution(format!("failed serializing GPX: {e}")))?;

        ensure_parent_dir(output)?;
        std::fs::write(output, &gpx)
            .map_err(|e| ToolError::Execution(format!("failed writing GPX file: {e}")))?;

        ctx.progress.info(&format!(
            "wrote {} waypoint(s) and {} track(s) ({} track point(s)) to {output}",
            counts.waypoints, counts.tracks, counts.track_points
        ));

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(output));
        outputs.insert("waypoint_count".to_string(), json!(counts.waypoints));
        outputs.insert("track_count".to_string(), json!(counts.tracks));
        outputs.insert("track_point_count".to_string(), json!(counts.track_points));
        outputs.insert("skipped_count".to_string(), json!(counts.skipped));
        Ok(ToolRunResult { outputs })
    }
}

/// Serializes a layer to a GPX 1.1 document, returning the bytes and feature
/// counts. Points (and MultiPoint parts) become waypoints; lines (and every
/// part of a MultiLineString) become one track with a `<trkseg>` per part.
fn serialize_gpx(layer: &Layer, fields: &FieldMap) -> Result<(Vec<u8>, Counts), quick_xml::Error> {
    let mut w = Writer::new_with_indent(Vec::new(), b' ', 2);
    w.write_event(Event::Decl(BytesDecl::new("1.0", Some("UTF-8"), None)))?;

    let mut gpx = BytesStart::new("gpx");
    gpx.push_attribute(("version", "1.1"));
    gpx.push_attribute(("creator", "geolibre"));
    gpx.push_attribute(("xmlns", "http://www.topografix.com/GPX/1/1"));
    w.write_event(Event::Start(gpx))?;

    let mut counts = Counts::default();
    let schema = &layer.schema;
    for feat in &layer.features {
        let Some(geom) = &feat.geometry else {
            counts.skipped += 1;
            continue;
        };
        let name = field_str(feat, schema, &fields.name);
        let desc = field_str(feat, schema, &fields.description);
        let time = field_str(feat, schema, &fields.date);
        let attr_z = field_f64(feat, schema, &fields.z);

        match geom {
            Geometry::Point(c) => {
                write_waypoint(&mut w, c, &name, &desc, &time, attr_z)?;
                counts.waypoints += 1;
            }
            Geometry::MultiPoint(cs) => {
                for c in cs {
                    write_waypoint(&mut w, c, &name, &desc, &time, attr_z)?;
                    counts.waypoints += 1;
                }
            }
            Geometry::LineString(_) | Geometry::MultiLineString(_) => {
                let parts = line_parts(geom);
                if parts.iter().all(|p| p.len() < 2) {
                    counts.skipped += 1;
                    continue;
                }
                w.write_event(Event::Start(BytesStart::new("trk")))?;
                if let Some(n) = &name {
                    text_el(&mut w, "name", n)?;
                }
                if let Some(d) = &desc {
                    text_el(&mut w, "desc", d)?;
                }
                for part in &parts {
                    if part.len() < 2 {
                        continue;
                    }
                    w.write_event(Event::Start(BytesStart::new("trkseg")))?;
                    for c in part {
                        write_trkpt(&mut w, c, &time, attr_z)?;
                        counts.track_points += 1;
                    }
                    w.write_event(Event::End(BytesEnd::new("trkseg")))?;
                }
                w.write_event(Event::End(BytesEnd::new("trk")))?;
                counts.tracks += 1;
            }
            _ => counts.skipped += 1,
        }
    }

    w.write_event(Event::End(BytesEnd::new("gpx")))?;
    Ok((w.into_inner(), counts))
}

/// Writes one `<wpt>` element. GPX schema order: ele, time, then name, desc.
fn write_waypoint(
    w: &mut Writer<Vec<u8>>,
    c: &Coord,
    name: &Option<String>,
    desc: &Option<String>,
    time: &Option<String>,
    attr_z: Option<f64>,
) -> Result<(), quick_xml::Error> {
    let mut wpt = BytesStart::new("wpt");
    wpt.push_attribute(("lat", fmt_num(c.y).as_str()));
    wpt.push_attribute(("lon", fmt_num(c.x).as_str()));
    w.write_event(Event::Start(wpt))?;
    if let Some(z) = attr_z.or(c.z) {
        text_el(w, "ele", &fmt_num(z))?;
    }
    if let Some(t) = time {
        text_el(w, "time", t)?;
    }
    if let Some(n) = name {
        text_el(w, "name", n)?;
    }
    if let Some(d) = desc {
        text_el(w, "desc", d)?;
    }
    w.write_event(Event::End(BytesEnd::new("wpt")))?;
    Ok(())
}

/// Writes one `<trkpt>` element. GPX schema order within trkpt: ele, time.
fn write_trkpt(
    w: &mut Writer<Vec<u8>>,
    c: &Coord,
    time: &Option<String>,
    attr_z: Option<f64>,
) -> Result<(), quick_xml::Error> {
    let mut pt = BytesStart::new("trkpt");
    pt.push_attribute(("lat", fmt_num(c.y).as_str()));
    pt.push_attribute(("lon", fmt_num(c.x).as_str()));
    // Per-vertex geometry Z is preferred; the z_field is a feature-level fallback.
    let ele = c.z.or(attr_z);
    if ele.is_none() && time.is_none() {
        w.write_event(Event::Empty(pt))?;
        return Ok(());
    }
    w.write_event(Event::Start(pt))?;
    if let Some(z) = ele {
        text_el(w, "ele", &fmt_num(z))?;
    }
    if let Some(t) = time {
        text_el(w, "time", t)?;
    }
    w.write_event(Event::End(BytesEnd::new("trkpt")))?;
    Ok(())
}

/// Writes `<name>text</name>` with automatic XML escaping of `text`.
fn text_el(w: &mut Writer<Vec<u8>>, name: &str, text: &str) -> Result<(), quick_xml::Error> {
    w.write_event(Event::Start(BytesStart::new(name)))?;
    w.write_event(Event::Text(BytesText::new(text)))?;
    w.write_event(Event::End(BytesEnd::new(name)))?;
    Ok(())
}

/// Extracts every line part from a LineString / MultiLineString geometry.
fn line_parts(geom: &Geometry) -> Vec<Vec<Coord>> {
    match geom {
        Geometry::LineString(cs) => vec![cs.clone()],
        Geometry::MultiLineString(ls) => ls.clone(),
        _ => Vec::new(),
    }
}

/// Reads the mapped field as a non-empty string, or `None` when unset/absent.
fn field_str(feat: &wbvector::Feature, schema: &Schema, field: &Option<String>) -> Option<String> {
    let name = field.as_ref()?;
    let s = field_to_string(feat.get(schema, name).ok()?);
    (!s.is_empty()).then_some(s)
}

/// Reads the mapped field as a float (numeric value or a parseable string).
fn field_f64(feat: &wbvector::Feature, schema: &Schema, field: &Option<String>) -> Option<f64> {
    let name = field.as_ref()?;
    let v = feat.get(schema, name).ok()?;
    v.as_f64()
        .or_else(|| v.as_str().and_then(|s| s.trim().parse::<f64>().ok()))
}

/// Renders any attribute value as text (empty for Null / blobs).
fn field_to_string(v: &FieldValue) -> String {
    use wbvector::FieldValue::*;
    match v {
        Null | Blob(_) => String::new(),
        Text(s) | Date(s) | DateTime(s) => s.clone(),
        Integer(i) => i.to_string(),
        Float(f) => fmt_num(*f),
        Boolean(b) => b.to_string(),
    }
}

/// Renders a coordinate/number without a trailing `.0` on whole values.
fn fmt_num(v: f64) -> String {
    if v == v.trunc() && v.is_finite() && v.abs() < 1e15 {
        format!("{}", v as i64)
    } else {
        let mut s = format!("{v:.8}");
        while s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
        s
    }
}

/// Trimmed string value of an argument, or "" when absent/non-string.
fn str_arg(args: &ToolArgs, key: &str) -> String {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .unwrap_or("")
        .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn out_path(tag: &str) -> std::path::PathBuf {
        use std::sync::atomic::{AtomicU64, Ordering};
        static SEQ: AtomicU64 = AtomicU64::new(0);
        let n = SEQ.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir().join(format!("ftgpx_{tag}_{}_{n}.gpx", std::process::id()))
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        FeaturesToGpxTool.run(&args, &ctx()).unwrap()
    }

    fn points_layer() -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.add_field(FieldDef::new("region", FieldType::Text));
        l.add_field(FieldDef::new("elev", FieldType::Float));
        l.add_field(FieldDef::new("when", FieldType::Text));
        l.add_feature(
            Some(Geometry::point(-74.0, 40.0)),
            &[
                ("name", "Alpha & Co".into()),
                ("region", "East".into()),
                ("elev", 12.5.into()),
                ("when", "2026-03-12T10:20:30Z".into()),
            ],
        )
        .unwrap();
        l.add_feature(
            Some(Geometry::point(-118.25, 34.05)),
            &[
                ("name", "Beta".into()),
                ("region", "West".into()),
                ("elev", 89.0.into()),
                ("when", "2026-03-12T11:00:00Z".into()),
            ],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn lines_layer() -> String {
        let mut l = Layer::new("lines")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.add_feature(
            Some(Geometry::LineString(vec![
                Coord::xy(-74.0, 40.0),
                Coord::xy(-74.05, 40.05),
                Coord::xy(-74.1, 40.1),
            ])),
            &[("name", "Trail 1".into())],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    #[test]
    fn points_become_waypoints_with_mapped_fields() {
        let out = out_path("wpt");
        let outs = run(json!({
            "input": points_layer(),
            "output": out.to_str().unwrap(),
            "name_field": "name",
            "description_field": "region",
            "z_field": "elev",
            "date_field": "when",
        }));
        assert_eq!(outs.outputs["waypoint_count"], json!(2));
        assert_eq!(outs.outputs["track_count"], json!(0));

        let text = std::fs::read_to_string(&out).unwrap();
        // Well-formed GPX round-trips through the reader.
        let parsed = wbvector::gpx::parse_str(&text).unwrap();
        assert_eq!(parsed.len(), 2);
        // Waypoint count == input point count.
        let wpts = parsed
            .features
            .iter()
            .filter(|f| matches!(f.geometry, Some(Geometry::Point(_))))
            .count();
        assert_eq!(wpts, 2);
        // Coordinates match (lon/lat preserved).
        let c = match &parsed.features[0].geometry {
            Some(Geometry::Point(c)) => c.clone(),
            _ => panic!("expected point"),
        };
        assert!((c.x - -74.0).abs() < 1e-9 && (c.y - 40.0).abs() < 1e-9);
        // Mapped metadata survived.
        assert_eq!(
            parsed.features[0]
                .get(&parsed.schema, "name")
                .unwrap()
                .as_str(),
            Some("Alpha & Co")
        );
        assert_eq!(
            parsed.features[0]
                .get(&parsed.schema, "desc")
                .unwrap()
                .as_str(),
            Some("East")
        );
        assert!((c.z.unwrap() - 12.5).abs() < 1e-9);
        // XML special char was escaped and decoded back intact.
        assert!(text.contains("Alpha &amp; Co"));
    }

    #[test]
    fn lines_become_tracks() {
        let out = out_path("trk");
        let outs = run(json!({
            "input": lines_layer(),
            "output": out.to_str().unwrap(),
            "name_field": "name",
        }));
        assert_eq!(outs.outputs["track_count"], json!(1));
        assert_eq!(outs.outputs["track_point_count"], json!(3));
        let text = std::fs::read_to_string(&out).unwrap();
        assert!(text.contains("<trk>") && text.contains("<trkseg>"));
        let parsed = wbvector::gpx::parse_str(&text).unwrap();
        assert_eq!(parsed.len(), 1);
        assert!(matches!(
            parsed.features[0].geometry,
            Some(Geometry::LineString(_))
        ));
    }

    #[test]
    fn no_mapping_still_writes_valid_gpx() {
        let out = out_path("plain");
        run(json!({
            "input": points_layer(),
            "output": out.to_str().unwrap(),
        }));
        let text = std::fs::read_to_string(&out).unwrap();
        let parsed = wbvector::gpx::parse_str(&text).unwrap();
        assert_eq!(parsed.len(), 2);
    }

    #[test]
    fn polygons_are_skipped() {
        let mut l = Layer::new("poly")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        l.add_feature(
            Some(Geometry::polygon(
                vec![
                    Coord::xy(0.0, 0.0),
                    Coord::xy(1.0, 0.0),
                    Coord::xy(1.0, 1.0),
                    Coord::xy(0.0, 0.0),
                ],
                vec![],
            )),
            &[],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        let path = memory_store::make_vector_memory_path(&id);
        let out = out_path("poly");
        let outs = run(json!({ "input": path, "output": out.to_str().unwrap() }));
        assert_eq!(outs.outputs["skipped_count"], json!(1));
        assert_eq!(outs.outputs["waypoint_count"], json!(0));
        assert_eq!(outs.outputs["track_count"], json!(0));
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = FeaturesToGpxTool;
        let v = |val: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(val).unwrap();
            tool.validate(&args)
        };
        // Missing input.
        assert!(v(json!({ "output": "x.gpx" })).is_err());
        // Missing output.
        assert!(v(json!({ "input": "x.geojson" })).is_err());
        // Valid.
        assert!(v(json!({ "input": "x.geojson", "output": "x.gpx" })).is_ok());
    }
}
