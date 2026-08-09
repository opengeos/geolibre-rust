//! GeoLibre tool: export a vector layer to KML (or zipped KMZ).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Layer To KML* (Conversion), which
//! also covers *Map To KML*.
//!
//! ## The asymmetry this closes
//!
//! `kml_to_features` shipped; nothing wrote KML back out. Every other exchange
//! format in the catalog is a symmetric pair — `gpx_to_features` /
//! `features_to_gpx`, `gtfs_to_features` / `features_to_gtfs`,
//! `read_geoparquet` / `write_geoparquet`. KML was the one-way street, even
//! though it is still the standard way to hand geometry to Google Earth and to
//! consumer and field viewers, and is frequently the requested deliverable for
//! anything the cartography and 3D suites produce.
//!
//! ## Coordinates are not optional
//!
//! KML is defined on **EPSG:4326 lon/lat**, always. Writing projected
//! coordinates into a KML produces a file that opens without complaint and
//! places everything thousands of kilometres away, which is the single easiest
//! way to get this wrong. A layer declaring a non-geographic EPSG is therefore
//! rejected outright rather than written; a layer with no declared CRS is
//! written with a warning, since its coordinates are all we have to go on.
//!
//! ## Ring orientation
//!
//! KML wants the outer boundary counter-clockwise and closed. Source rings vary,
//! so orientation is normalized from the signed area rather than trusted, and
//! the closing vertex is appended when missing.

use std::collections::BTreeMap;
use std::io::Write as _;

use quick_xml::events::{BytesDecl, BytesEnd, BytesStart, BytesText, Event};
use quick_xml::Writer;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldValue, Geometry, Layer, Ring};

use crate::args_common::choice_or;
use crate::vector_common::{ensure_parent_dir, load_input_layer, parse_optional_str};

const ALTITUDE_MODES: [&str; 3] = ["clamptoground", "relativetoground", "absolute"];

/// KML spells these in camelCase; the parameter is matched case-insensitively.
fn altitude_mode_kml(m: &str) -> &'static str {
    match m {
        "relativetoground" => "relativeToGround",
        "absolute" => "absolute",
        _ => "clampToGround",
    }
}

pub struct FeaturesToKmlTool;

impl Tool for FeaturesToKmlTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "features_to_kml",
            display_name: "Features To KML",
            summary: "Writes a vector layer to a KML 2.2 document (or a zipped KMZ), with attributes as ExtendedData and an optional shared style (ArcGIS Layer To KML / Map To KML). kml_to_features could already read KML but nothing wrote it, unlike the symmetric GPX, GTFS and GeoParquet pairs.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer. Must be in EPSG:4326 lon/lat, which is the only CRS KML defines.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output .kml or .kmz file path.",
                    required: true,
                },
                ToolParamSpec {
                    name: "name_field",
                    description: "Attribute used for each Placemark's <name>.",
                    required: false,
                },
                ToolParamSpec {
                    name: "description_field",
                    description: "Attribute used for <description>. When omitted, all attributes are written as <ExtendedData>.",
                    required: false,
                },
                ToolParamSpec {
                    name: "z_field",
                    description: "Attribute supplying an elevation, written as the third coordinate.",
                    required: false,
                },
                ToolParamSpec {
                    name: "altitude_mode",
                    description: "'clampToGround' (default), 'relativeToGround', or 'absolute'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "color",
                    description: "Line/outline colour as KML aabbggrr hex (e.g. 'ff0000ff' for opaque red).",
                    required: false,
                },
                ToolParamSpec {
                    name: "line_width",
                    description: "Line width in pixels (default 1.0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "fill_color",
                    description: "Polygon fill colour as KML aabbggrr hex. Omit for no fill.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if parse_optional_str(args, "input")?.is_none() {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        let Some(output) = parse_optional_str(args, "output")? else {
            return Err(ToolError::Validation(
                "missing required string parameter 'output'".to_string(),
            ));
        };
        let lower = output.to_ascii_lowercase();
        if !(lower.ends_with(".kml") || lower.ends_with(".kmz")) {
            return Err(ToolError::Validation(
                "parameter 'output' must be a .kml or .kmz file path".to_string(),
            ));
        }
        choice_or(args, "altitude_mode", &ALTITUDE_MODES, "clamptoground")?;
        for key in ["color", "fill_color"] {
            if let Some(c) = parse_optional_str(args, key)? {
                validate_kml_color(key, c)?;
            }
        }
        if let Some(w) = args.get("line_width").and_then(Value::as_f64) {
            if !w.is_finite() || w <= 0.0 {
                return Err(ToolError::Validation(
                    "'line_width' must be a positive number".to_string(),
                ));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = parse_optional_str(args, "input")?
            .ok_or_else(|| ToolError::Validation("missing required parameter 'input'".into()))?;
        let output = parse_optional_str(args, "output")?
            .ok_or_else(|| ToolError::Validation("missing required parameter 'output'".into()))?;
        let style = Style {
            name_field: parse_optional_str(args, "name_field")?.map(str::to_string),
            description_field: parse_optional_str(args, "description_field")?.map(str::to_string),
            z_field: parse_optional_str(args, "z_field")?.map(str::to_string),
            altitude_mode: altitude_mode_kml(choice_or(
                args,
                "altitude_mode",
                &ALTITUDE_MODES,
                "clamptoground",
            )?),
            color: parse_optional_str(args, "color")?.map(str::to_string),
            fill_color: parse_optional_str(args, "fill_color")?.map(str::to_string),
            line_width: args
                .get("line_width")
                .and_then(Value::as_f64)
                .unwrap_or(1.0),
        };

        let layer = load_input_layer(input)?;

        // The single most damaging mistake this tool could make: emit projected
        // coordinates into a format defined on lon/lat.
        let mut crs_warning = None;
        match layer.crs_epsg() {
            Some(4326) => {}
            Some(other) => {
                return Err(ToolError::Validation(format!(
                    "KML is defined on EPSG:4326 lon/lat, but the input layer is EPSG:{other}. \
                     Reproject it to EPSG:4326 first — writing projected coordinates would \
                     produce a KML that opens but places every feature in the wrong place."
                )));
            }
            None => {
                crs_warning = Some(
                    "the input layer declares no CRS; its coordinates were written as lon/lat \
                     without conversion",
                );
                ctx.progress.info(crs_warning.unwrap());
            }
        }
        // A cheap sanity check that catches the common case of an undeclared
        // projected layer, where the numbers are obviously not degrees.
        if let Some(bad) = out_of_degree_range(&layer) {
            return Err(ToolError::Validation(format!(
                "coordinate ({:.3}, {:.3}) is outside the lon/lat range KML requires; the input \
                 looks projected. Reproject it to EPSG:4326 first.",
                bad.0, bad.1
            )));
        }

        let (bytes, counts) = serialize_kml(&layer, &style)
            .map_err(|e| ToolError::Execution(format!("failed serializing KML: {e}")))?;

        ensure_parent_dir(output)?;
        if output.to_ascii_lowercase().ends_with(".kmz") {
            write_kmz(output, &bytes)?;
        } else {
            std::fs::write(output, &bytes)
                .map_err(|e| ToolError::Execution(format!("failed writing KML file: {e}")))?;
        }

        ctx.progress.info(&format!(
            "wrote {} placemark(s) to {output}",
            counts.placemarks
        ));

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(output));
        outputs.insert("placemark_count".to_string(), json!(counts.placemarks));
        outputs.insert("point_count".to_string(), json!(counts.points));
        outputs.insert("line_count".to_string(), json!(counts.lines));
        outputs.insert("polygon_count".to_string(), json!(counts.polygons));
        outputs.insert("skipped_count".to_string(), json!(counts.skipped));
        if let Some(w) = crs_warning {
            outputs.insert("warning".to_string(), json!(w));
        }
        Ok(ToolRunResult { outputs })
    }
}

struct Style {
    name_field: Option<String>,
    description_field: Option<String>,
    z_field: Option<String>,
    altitude_mode: &'static str,
    color: Option<String>,
    fill_color: Option<String>,
    line_width: f64,
}

impl Style {
    fn has_style(&self) -> bool {
        self.color.is_some() || self.fill_color.is_some()
    }
}

#[derive(Default)]
struct Counts {
    placemarks: u64,
    points: u64,
    lines: u64,
    polygons: u64,
    skipped: u64,
}

/// KML colours are `aabbggrr` — 8 hex digits, alpha first and the RGB channels
/// reversed relative to HTML. Accepting a `#rrggbb` here would silently swap
/// red and blue, so only the KML form is allowed.
fn validate_kml_color(key: &str, c: &str) -> Result<(), ToolError> {
    if c.len() == 8 && c.chars().all(|ch| ch.is_ascii_hexdigit()) {
        return Ok(());
    }
    Err(ToolError::Validation(format!(
        "'{key}' must be 8 hex digits in KML aabbggrr order (alpha, blue, green, red), \
         e.g. 'ff0000ff' for opaque red; got '{c}'"
    )))
}

/// Returns the first coordinate outside the lon/lat domain, if any.
fn out_of_degree_range(layer: &Layer) -> Option<(f64, f64)> {
    fn check(coords: &[Coord]) -> Option<(f64, f64)> {
        coords
            .iter()
            .find(|c| !(-180.0..=180.0).contains(&c.x) || !(-90.0..=90.0).contains(&c.y))
            .map(|c| (c.x, c.y))
    }
    for f in layer.iter() {
        let hit = match f.geometry.as_ref() {
            Some(Geometry::Point(p)) => check(std::slice::from_ref(p)),
            Some(Geometry::MultiPoint(ps)) => check(ps),
            Some(Geometry::LineString(cs)) => check(cs),
            Some(Geometry::MultiLineString(ls)) => ls.iter().find_map(|l| check(l)),
            Some(Geometry::Polygon {
                exterior,
                interiors,
            }) => check(&exterior.0).or_else(|| interiors.iter().find_map(|r| check(&r.0))),
            Some(Geometry::MultiPolygon(ps)) => ps
                .iter()
                .find_map(|(e, hs)| check(&e.0).or_else(|| hs.iter().find_map(|r| check(&r.0)))),
            _ => None,
        };
        if hit.is_some() {
            return hit;
        }
    }
    None
}

fn serialize_kml(layer: &Layer, style: &Style) -> Result<(Vec<u8>, Counts), quick_xml::Error> {
    let mut w = Writer::new_with_indent(Vec::new(), b' ', 2);
    w.write_event(Event::Decl(BytesDecl::new("1.0", Some("UTF-8"), None)))?;

    let mut kml = BytesStart::new("kml");
    kml.push_attribute(("xmlns", "http://www.opengis.net/kml/2.2"));
    w.write_event(Event::Start(kml))?;
    w.write_event(Event::Start(BytesStart::new("Document")))?;
    text_el(&mut w, "name", &layer.name)?;

    if style.has_style() {
        w.write_event(Event::Start({
            let mut s = BytesStart::new("Style");
            s.push_attribute(("id", "geolibre"));
            s
        }))?;
        if let Some(c) = &style.color {
            w.write_event(Event::Start(BytesStart::new("LineStyle")))?;
            text_el(&mut w, "color", c)?;
            text_el(&mut w, "width", &fmt_num(style.line_width))?;
            w.write_event(Event::End(BytesEnd::new("LineStyle")))?;
        }
        w.write_event(Event::Start(BytesStart::new("PolyStyle")))?;
        match &style.fill_color {
            Some(c) => {
                text_el(&mut w, "color", c)?;
                text_el(&mut w, "fill", "1")?;
            }
            None => text_el(&mut w, "fill", "0")?,
        }
        text_el(&mut w, "outline", "1")?;
        w.write_event(Event::End(BytesEnd::new("PolyStyle")))?;
        w.write_event(Event::End(BytesEnd::new("Style")))?;
    }

    let mut counts = Counts::default();
    for feature in layer.iter() {
        // Every Geometry variant is writable as KML, so the only feature that
        // cannot be emitted is one carrying no geometry at all.
        let Some(geom) = feature.geometry.as_ref() else {
            counts.skipped += 1;
            continue;
        };
        let z = style
            .z_field
            .as_ref()
            .and_then(|f| layer.schema.field_index(f))
            .and_then(|i| feature.attributes.get(i))
            .and_then(field_as_f64);

        w.write_event(Event::Start(BytesStart::new("Placemark")))?;
        if let Some(name) = field_text(layer, feature, &style.name_field) {
            text_el(&mut w, "name", &name)?;
        }
        if let Some(desc) = field_text(layer, feature, &style.description_field) {
            text_el(&mut w, "description", &desc)?;
        } else {
            write_extended_data(&mut w, layer, feature)?;
        }
        if style.has_style() {
            text_el(&mut w, "styleUrl", "#geolibre")?;
        }
        write_geometry(&mut w, geom, style, z, &mut counts)?;
        w.write_event(Event::End(BytesEnd::new("Placemark")))?;
        counts.placemarks += 1;
    }

    w.write_event(Event::End(BytesEnd::new("Document")))?;
    w.write_event(Event::End(BytesEnd::new("kml")))?;
    Ok((w.into_inner(), counts))
}

fn write_geometry(
    w: &mut Writer<Vec<u8>>,
    geom: &Geometry,
    style: &Style,
    z: Option<f64>,
    counts: &mut Counts,
) -> Result<(), quick_xml::Error> {
    match geom {
        Geometry::Point(p) => {
            w.write_event(Event::Start(BytesStart::new("Point")))?;
            text_el(w, "altitudeMode", style.altitude_mode)?;
            text_el(w, "coordinates", &coord_str(std::slice::from_ref(p), z))?;
            w.write_event(Event::End(BytesEnd::new("Point")))?;
            counts.points += 1;
        }
        Geometry::LineString(cs) => {
            w.write_event(Event::Start(BytesStart::new("LineString")))?;
            text_el(w, "altitudeMode", style.altitude_mode)?;
            text_el(w, "coordinates", &coord_str(cs, z))?;
            w.write_event(Event::End(BytesEnd::new("LineString")))?;
            counts.lines += 1;
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            if ring_is_writable(exterior) {
                write_polygon(w, exterior, interiors, style, z)?;
                counts.polygons += 1;
            } else {
                counts.skipped += 1;
            }
        }
        // Multi-part geometries become a <MultiGeometry> of their parts.
        Geometry::MultiPoint(ps) => {
            w.write_event(Event::Start(BytesStart::new("MultiGeometry")))?;
            for p in ps {
                write_geometry(w, &Geometry::Point(p.clone()), style, z, counts)?;
            }
            w.write_event(Event::End(BytesEnd::new("MultiGeometry")))?;
        }
        Geometry::MultiLineString(ls) => {
            w.write_event(Event::Start(BytesStart::new("MultiGeometry")))?;
            for l in ls {
                write_geometry(w, &Geometry::LineString(l.clone()), style, z, counts)?;
            }
            w.write_event(Event::End(BytesEnd::new("MultiGeometry")))?;
        }
        Geometry::MultiPolygon(polys) => {
            w.write_event(Event::Start(BytesStart::new("MultiGeometry")))?;
            for (ext, holes) in polys {
                if ring_is_writable(ext) {
                    write_polygon(w, ext, holes, style, z)?;
                    counts.polygons += 1;
                } else {
                    counts.skipped += 1;
                }
            }
            w.write_event(Event::End(BytesEnd::new("MultiGeometry")))?;
        }
        Geometry::GeometryCollection(gs) => {
            w.write_event(Event::Start(BytesStart::new("MultiGeometry")))?;
            for g in gs {
                write_geometry(w, g, style, z, counts)?;
            }
            w.write_event(Event::End(BytesEnd::new("MultiGeometry")))?;
        }
    }
    Ok(())
}

/// KML requires a LinearRing to hold at least four coordinate tuples (three
/// distinct vertices plus the closing repeat). Anything shorter produces a
/// document strict validators reject.
fn ring_is_writable(ring: &Ring) -> bool {
    ring.0.len() >= 3
}

fn write_polygon(
    w: &mut Writer<Vec<u8>>,
    exterior: &Ring,
    interiors: &[Ring],
    style: &Style,
    z: Option<f64>,
) -> Result<(), quick_xml::Error> {
    w.write_event(Event::Start(BytesStart::new("Polygon")))?;
    text_el(w, "altitudeMode", style.altitude_mode)?;

    w.write_event(Event::Start(BytesStart::new("outerBoundaryIs")))?;
    w.write_event(Event::Start(BytesStart::new("LinearRing")))?;
    text_el(
        w,
        "coordinates",
        &coord_str(&oriented(&exterior.0, true), z),
    )?;
    w.write_event(Event::End(BytesEnd::new("LinearRing")))?;
    w.write_event(Event::End(BytesEnd::new("outerBoundaryIs")))?;

    for hole in interiors.iter().filter(|r| ring_is_writable(r)) {
        w.write_event(Event::Start(BytesStart::new("innerBoundaryIs")))?;
        w.write_event(Event::Start(BytesStart::new("LinearRing")))?;
        text_el(w, "coordinates", &coord_str(&oriented(&hole.0, false), z))?;
        w.write_event(Event::End(BytesEnd::new("LinearRing")))?;
        w.write_event(Event::End(BytesEnd::new("innerBoundaryIs")))?;
    }

    w.write_event(Event::End(BytesEnd::new("Polygon")))?;
    Ok(())
}

/// Closes the ring and forces the winding KML expects: counter-clockwise for an
/// outer boundary, clockwise for a hole. Source rings vary by format and by
/// producer, so this is normalized rather than trusted.
fn oriented(coords: &[Coord], ccw: bool) -> Vec<Coord> {
    let mut out = coords.to_vec();
    if out.len() >= 2 {
        let (first, last) = (out[0].clone(), out[out.len() - 1].clone());
        if (first.x - last.x).abs() > f64::EPSILON || (first.y - last.y).abs() > f64::EPSILON {
            out.push(first.clone());
        }
    }
    // Shoelace: positive area means counter-clockwise.
    let mut area2 = 0.0;
    for i in 0..out.len() {
        let a = &out[i];
        let b = &out[(i + 1) % out.len()];
        area2 += a.x * b.y - b.x * a.y;
    }
    let is_ccw = area2 > 0.0;
    if is_ccw != ccw {
        out.reverse();
    }
    out
}

/// KML coordinate tuples: `lon,lat[,alt]` comma-separated within a tuple and
/// whitespace-separated between tuples. Getting that inverted is a classic
/// silent corruption, so the two separators are deliberately explicit here.
fn coord_str(coords: &[Coord], z: Option<f64>) -> String {
    let mut s = String::with_capacity(coords.len() * 24);
    for (i, c) in coords.iter().enumerate() {
        if i > 0 {
            s.push(' ');
        }
        s.push_str(&fmt_num(c.x));
        s.push(',');
        s.push_str(&fmt_num(c.y));
        if let Some(z) = z {
            s.push(',');
            s.push_str(&fmt_num(z));
        }
    }
    s
}

fn fmt_num(v: f64) -> String {
    if v == v.trunc() && v.abs() < 1e15 {
        format!("{}", v as i64)
    } else {
        let s = format!("{v:.8}");
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    }
}

fn write_extended_data(
    w: &mut Writer<Vec<u8>>,
    layer: &Layer,
    feature: &Feature,
) -> Result<(), quick_xml::Error> {
    let fields = layer.schema.fields();
    if fields.is_empty() {
        return Ok(());
    }
    w.write_event(Event::Start(BytesStart::new("ExtendedData")))?;
    for (i, f) in fields.iter().enumerate() {
        let mut data = BytesStart::new("Data");
        data.push_attribute(("name", f.name.as_str()));
        w.write_event(Event::Start(data))?;
        let value = feature.attributes.get(i).map(display).unwrap_or_default();
        text_el(w, "value", &value)?;
        w.write_event(Event::End(BytesEnd::new("Data")))?;
    }
    w.write_event(Event::End(BytesEnd::new("ExtendedData")))?;
    Ok(())
}

/// Writes `<name>text</name>`, letting quick-xml escape the text. Attribute
/// values routinely contain `&` and `<`, which would otherwise produce a KML
/// that no parser will open.
fn text_el(w: &mut Writer<Vec<u8>>, name: &str, text: &str) -> Result<(), quick_xml::Error> {
    w.write_event(Event::Start(BytesStart::new(name)))?;
    w.write_event(Event::Text(BytesText::new(text)))?;
    w.write_event(Event::End(BytesEnd::new(name)))?;
    Ok(())
}

fn field_text(layer: &Layer, feature: &Feature, field: &Option<String>) -> Option<String> {
    let name = field.as_ref()?;
    let i = layer.schema.field_index(name)?;
    let s = display(feature.attributes.get(i)?);
    (!s.is_empty()).then_some(s)
}

fn field_as_f64(v: &FieldValue) -> Option<f64> {
    match v {
        FieldValue::Integer(i) => Some(*i as f64),
        FieldValue::Float(f) if f.is_finite() => Some(*f),
        // "NaN" and "inf" both parse successfully; fmt_num would then emit
        // them as the third coordinate and KML readers reject the tuple.
        FieldValue::Text(s) => s.trim().parse::<f64>().ok().filter(|v| v.is_finite()),
        _ => None,
    }
}

fn display(v: &FieldValue) -> String {
    match v {
        FieldValue::Null => String::new(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) => fmt_num(*f),
        FieldValue::Text(s) => s.clone(),
        FieldValue::Boolean(b) => b.to_string(),
        FieldValue::Date(s) | FieldValue::DateTime(s) => s.clone(),
        FieldValue::Blob(b) => format!("blob[{}]", b.len()),
    }
}

/// A KMZ is a zip holding a single `doc.kml`. `kml_to_features` already reads
/// this layout, so the pair round-trips.
fn write_kmz(path: &str, kml: &[u8]) -> Result<(), ToolError> {
    let file = std::fs::File::create(path)
        .map_err(|e| ToolError::Execution(format!("failed creating KMZ file: {e}")))?;
    let mut zip = zip::ZipWriter::new(file);
    let opts: zip::write::FileOptions<'_, ()> =
        zip::write::FileOptions::default().compression_method(zip::CompressionMethod::Deflated);
    zip.start_file("doc.kml", opts)
        .map_err(|e| ToolError::Execution(format!("failed starting KMZ entry: {e}")))?;
    zip.write_all(kml)
        .map_err(|e| ToolError::Execution(format!("failed writing KMZ entry: {e}")))?;
    zip.finish()
        .map_err(|e| ToolError::Execution(format!("failed finalizing KMZ: {e}")))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{FieldDef, FieldType, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn tmp(tag: &str, ext: &str) -> String {
        static N: AtomicU64 = AtomicU64::new(0);
        let n = N.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir()
            .join(format!("ftkml_{tag}_{}_{n}.{ext}", std::process::id()))
            .to_string_lossy()
            .to_string()
    }

    fn store(l: Layer) -> String {
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn points_layer(epsg: Option<u32>) -> String {
        let mut l = Layer::new("sites").with_geom_type(GeometryType::Point);
        if let Some(e) = epsg {
            l = l.with_crs_epsg(e);
        }
        l.add_field(FieldDef::new("label", FieldType::Text));
        l.add_field(FieldDef::new("elev", FieldType::Float));
        l.add_feature(
            Some(Geometry::point(-122.5, 37.7)),
            &[("label", "Ferry & Dock".into()), ("elev", 12.5f64.into())],
        )
        .unwrap();
        store(l)
    }

    fn run(args: Value) -> (String, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = FeaturesToKmlTool.run(&args, &ctx()).unwrap();
        let path = res.outputs["output"].as_str().unwrap().to_string();
        let text = std::fs::read_to_string(&path).unwrap();
        let _ = std::fs::remove_file(&path);
        (text, res)
    }

    #[test]
    fn writes_a_kml_document_with_a_placemark() {
        let out = tmp("basic", "kml");
        let (text, res) = run(json!({"input": points_layer(Some(4326)), "output": out}));
        assert!(text.contains("<kml xmlns=\"http://www.opengis.net/kml/2.2\">"));
        assert!(text.contains("<Placemark>"));
        assert!(text.contains("<Point>"));
        assert_eq!(res.outputs["placemark_count"], json!(1));
    }

    #[test]
    fn coordinates_are_lon_lat_comma_separated() {
        // The classic silent corruption is emitting lat,lon. -122.5 is the
        // longitude and must come first.
        let out = tmp("coords", "kml");
        let (text, _) = run(json!({"input": points_layer(Some(4326)), "output": out}));
        assert!(
            text.contains("<coordinates>-122.5,37.7</coordinates>"),
            "got: {text}"
        );
    }

    #[test]
    fn a_z_field_becomes_the_third_coordinate() {
        let out = tmp("z", "kml");
        let (text, _) = run(json!({
            "input": points_layer(Some(4326)), "output": out,
            "z_field": "elev", "altitude_mode": "absolute",
        }));
        assert!(text.contains("-122.5,37.7,12.5"), "got: {text}");
        assert!(text.contains("<altitudeMode>absolute</altitudeMode>"));
    }

    #[test]
    fn attribute_text_is_xml_escaped() {
        // "Ferry & Dock" must not produce a bare ampersand, which would make
        // the file unparseable.
        let out = tmp("escape", "kml");
        let (text, _) = run(json!({
            "input": points_layer(Some(4326)), "output": out, "name_field": "label",
        }));
        assert!(text.contains("Ferry &amp; Dock"), "got: {text}");
        assert!(!text.contains("Ferry & Dock"));
    }

    #[test]
    fn attributes_land_in_extended_data_when_no_description_field_is_given() {
        let out = tmp("ext", "kml");
        let (text, _) = run(json!({"input": points_layer(Some(4326)), "output": out}));
        assert!(text.contains("<ExtendedData>"));
        assert!(text.contains("name=\"label\""));
    }

    #[test]
    fn a_polygon_outer_ring_is_closed_and_counter_clockwise() {
        // Input is deliberately clockwise and unclosed.
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        l.add_feature(
            Some(Geometry::polygon(
                vec![
                    Coord::xy(0.0, 0.0),
                    Coord::xy(0.0, 1.0),
                    Coord::xy(1.0, 1.0),
                    Coord::xy(1.0, 0.0),
                ],
                Vec::new(),
            )),
            &[],
        )
        .unwrap();
        let out = tmp("poly", "kml");
        let (text, _) = run(json!({"input": store(l), "output": out}));
        let coords = text
            .split("<coordinates>")
            .nth(1)
            .unwrap()
            .split("</coordinates>")
            .next()
            .unwrap();
        let pts: Vec<&str> = coords.split_whitespace().collect();
        assert_eq!(pts.first(), pts.last(), "ring must be closed");
        // Counter-clockwise from (0,0): the second vertex is (1,0), not (0,1).
        assert_eq!(pts[1], "1,0", "ring was not reoriented to CCW: {coords}");
    }

    #[test]
    fn a_hole_is_written_as_an_inner_boundary() {
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        l.add_feature(
            Some(Geometry::polygon(
                vec![
                    Coord::xy(0.0, 0.0),
                    Coord::xy(4.0, 0.0),
                    Coord::xy(4.0, 4.0),
                    Coord::xy(0.0, 4.0),
                    Coord::xy(0.0, 0.0),
                ],
                vec![vec![
                    Coord::xy(1.0, 1.0),
                    Coord::xy(2.0, 1.0),
                    Coord::xy(2.0, 2.0),
                    Coord::xy(1.0, 2.0),
                    Coord::xy(1.0, 1.0),
                ]],
            )),
            &[],
        )
        .unwrap();
        let out = tmp("hole", "kml");
        let (text, _) = run(json!({"input": store(l), "output": out}));
        assert!(text.contains("<outerBoundaryIs>"));
        assert!(text.contains("<innerBoundaryIs>"));
    }

    #[test]
    fn a_degenerate_ring_is_skipped_rather_than_emitting_invalid_kml() {
        // KML requires a LinearRing to carry at least four coordinate tuples;
        // a two-vertex ring produces a document strict validators reject.
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        l.add_feature(
            Some(Geometry::polygon(
                vec![Coord::xy(0.0, 0.0), Coord::xy(1.0, 1.0)],
                Vec::new(),
            )),
            &[],
        )
        .unwrap();
        let out = tmp("degen", "kml");
        let (text, res) = run(json!({"input": store(l), "output": out}));
        assert_eq!(res.outputs["polygon_count"], json!(0));
        assert_eq!(res.outputs["skipped_count"], json!(1));
        assert!(!text.contains("<LinearRing>"), "got: {text}");
    }

    #[test]
    fn a_multi_geometry_wraps_its_parts() {
        let mut l = Layer::new("m")
            .with_geom_type(GeometryType::MultiPoint)
            .with_crs_epsg(4326);
        l.add_feature(
            Some(Geometry::MultiPoint(vec![
                Coord::xy(1.0, 2.0),
                Coord::xy(3.0, 4.0),
            ])),
            &[],
        )
        .unwrap();
        let out = tmp("multi", "kml");
        let (text, res) = run(json!({"input": store(l), "output": out}));
        assert!(text.contains("<MultiGeometry>"));
        assert_eq!(res.outputs["point_count"], json!(2));
    }

    #[test]
    fn a_projected_layer_is_refused_rather_than_silently_misplaced() {
        let mut l = Layer::new("utm")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(32610);
        l.add_feature(Some(Geometry::point(550000.0, 4180000.0)), &[])
            .unwrap();
        let args: ToolArgs = serde_json::from_value(json!({
            "input": store(l), "output": tmp("utm", "kml"),
        }))
        .unwrap();
        let err = FeaturesToKmlTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("EPSG:4326"), "{err}");
    }

    #[test]
    fn undeclared_but_obviously_projected_coordinates_are_caught_too() {
        // No CRS declared, so the EPSG check cannot fire — the range check is
        // the only thing standing between the user and a nonsense KML.
        let mut l = Layer::new("nocrs").with_geom_type(GeometryType::Point);
        l.add_feature(Some(Geometry::point(550000.0, 4180000.0)), &[])
            .unwrap();
        let args: ToolArgs = serde_json::from_value(json!({
            "input": store(l), "output": tmp("nocrs", "kml"),
        }))
        .unwrap();
        let err = FeaturesToKmlTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("looks projected"), "{err}");
    }

    #[test]
    fn a_style_is_emitted_and_referenced() {
        let out = tmp("style", "kml");
        let (text, _) = run(json!({
            "input": points_layer(Some(4326)), "output": out,
            "color": "ff0000ff", "line_width": 3.0, "fill_color": "8000ff00",
        }));
        assert!(text.contains("<Style id=\"geolibre\">"));
        assert!(text.contains("<color>ff0000ff</color>"));
        assert!(text.contains("<width>3</width>"));
        assert!(text.contains("<styleUrl>#geolibre</styleUrl>"));
    }

    #[test]
    fn kmz_output_is_a_zip_that_kml_to_features_can_read_back() {
        let out = tmp("kmz", "kmz");
        let args: ToolArgs = serde_json::from_value(json!({
            "input": points_layer(Some(4326)), "output": out.clone(),
        }))
        .unwrap();
        FeaturesToKmlTool.run(&args, &ctx()).unwrap();
        let bytes = std::fs::read(&out).unwrap();
        assert_eq!(&bytes[..2], b"PK", "KMZ must be a zip");

        // Round-trip through the reader half of the pair.
        let args: ToolArgs = serde_json::from_value(json!({"input": out.clone()})).unwrap();
        let res = crate::kml_to_features::KmlToFeaturesTool
            .run(&args, &ctx())
            .unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        let _ = std::fs::remove_file(&out);
        assert_eq!(layer.features.len(), 1);
        let Some(Geometry::Point(p)) = layer.features[0].geometry.as_ref() else {
            panic!("expected a point back");
        };
        assert!((p.x - -122.5).abs() < 1e-9 && (p.y - 37.7).abs() < 1e-9);
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            FeaturesToKmlTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"input": "a.shp"})));
        assert!(bad(json!({"input": "a.shp", "output": "a.shp"})));
        assert!(bad(
            json!({"input": "a.shp", "output": "a.kml", "altitude_mode": "sky"})
        ));
        // #rrggbb would swap red and blue against KML's aabbggrr order.
        assert!(bad(
            json!({"input": "a.shp", "output": "a.kml", "color": "#ff0000"})
        ));
        assert!(bad(
            json!({"input": "a.shp", "output": "a.kml", "line_width": 0})
        ));
    }
}
