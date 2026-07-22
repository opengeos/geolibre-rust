//! GeoLibre tool: read KML/KMZ Placemarks into point/line/polygon layers.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *KML To Layer* (Conversion): parse a
//! `.kml` file (or a `.kmz`, which is just a zipped KML) and split its
//! `<Placemark>` geometry into up to three `wbvector` layers — points, lines,
//! and polygons — each attributed with the placemark's `name` and
//! `description`.
//!
//! KML is a top-tier public GIS interchange format and is plain XML, so the work
//! (XML pull-parsing + zip unpacking for KMZ) is squarely in the pure-Rust/WASM
//! stack. Neither the repo nor the bundled whitebox-wasm suite reads KML.
//!
//! KML coordinates are always WGS84 and written `lon,lat[,alt]`, so output is
//! EPSG:4326 and the tool honours the KML axis order (longitude first).
//!
//! Geometry inside a `<MultiGeometry>` is flattened: each primitive becomes its
//! own feature carrying the parent placemark's attributes. `<Point>`,
//! `<LineString>`, and `<Polygon>` (with `<innerBoundaryIs>` holes) are read;
//! other overlays (GroundOverlay, Model, …) are ignored.

use std::collections::BTreeMap;
use std::io::Read;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{parse_optional_str, write_or_store_layer};

pub struct KmlToFeaturesTool;

impl Tool for KmlToFeaturesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "kml_to_features",
            display_name: "KML To Features",
            summary: "Read .kml/.kmz Placemarks into point, line, and polygon layers with name/description attributes.",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Path to a .kml file or a .kmz archive (a zipped KML).",
                    required: true,
                },
                ToolParamSpec {
                    name: "points_output",
                    description: "Optional output point vector path (driver from its extension). Written only if the KML has point placemarks.",
                    required: false,
                },
                ToolParamSpec {
                    name: "lines_output",
                    description: "Optional output line vector path (driver from its extension). Written only if the KML has line placemarks.",
                    required: false,
                },
                ToolParamSpec {
                    name: "polygons_output",
                    description: "Optional output polygon vector path (driver from its extension). Written only if the KML has polygon placemarks.",
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
        let points_output = parse_optional_str(args, "points_output")?;
        let lines_output = parse_optional_str(args, "lines_output")?;
        let polygons_output = parse_optional_str(args, "polygons_output")?;

        let kml = read_kml_text(input)?;
        let placemarks = parse_placemarks(&kml)?;
        ctx.progress
            .info(&format!("parsed {} placemark(s)", placemarks.len()));

        // Build one layer per OGC geometry class the KML produced.
        let mut points = new_layer("kml_points", GeometryType::Point);
        let mut lines = new_layer("kml_lines", GeometryType::LineString);
        let mut polygons = new_layer("kml_polygons", GeometryType::Polygon);

        for pm in &placemarks {
            for geom in &pm.geoms {
                let attrs = [
                    ("name", FieldValue::Text(pm.name.clone())),
                    ("description", FieldValue::Text(pm.description.clone())),
                ];
                let layer = match geom.geom_type() {
                    GeometryType::Point => &mut points,
                    GeometryType::LineString => &mut lines,
                    GeometryType::Polygon => &mut polygons,
                    _ => continue,
                };
                layer
                    .add_feature(Some(geom.clone()), &attrs)
                    .map_err(|e| ToolError::Execution(format!("failed adding feature: {e}")))?;
            }
        }

        let point_count = points.len();
        let line_count = lines.len();
        let polygon_count = polygons.len();

        // Write each non-empty layer to its path (always, if a path was given) or
        // stash it in memory when the caller passed no path at all.
        let points_path = emit(points, points_output, point_count)?;
        let lines_path = emit(lines, lines_output, line_count)?;
        let polygons_path = emit(polygons, polygons_output, polygon_count)?;

        let mut outputs = BTreeMap::new();
        // "output" points at the first layer that materialised, for convenience.
        if let Some(p) = points_path
            .as_ref()
            .or(lines_path.as_ref())
            .or(polygons_path.as_ref())
        {
            outputs.insert("output".to_string(), json!(p));
        }
        if let Some(p) = points_path {
            outputs.insert("points_output".to_string(), json!(p));
        }
        if let Some(p) = lines_path {
            outputs.insert("lines_output".to_string(), json!(p));
        }
        if let Some(p) = polygons_path {
            outputs.insert("polygons_output".to_string(), json!(p));
        }
        outputs.insert("placemark_count".to_string(), json!(placemarks.len()));
        outputs.insert("point_count".to_string(), json!(point_count));
        outputs.insert("line_count".to_string(), json!(line_count));
        outputs.insert("polygon_count".to_string(), json!(polygon_count));
        Ok(ToolRunResult { outputs })
    }
}

/// Creates a name/description-attributed WGS84 layer of the given geometry type.
fn new_layer(name: &str, geom: GeometryType) -> Layer {
    let mut layer = Layer::new(name).with_geom_type(geom).with_crs_epsg(4326);
    layer.add_field(FieldDef::new("name", FieldType::Text));
    layer.add_field(FieldDef::new("description", FieldType::Text));
    layer
}

/// Writes `layer` when it has features (to `path`, or to memory when `path` is
/// `None`); returns the resulting path/handle, or `None` when the layer is empty.
fn emit(layer: Layer, path: Option<&str>, count: usize) -> Result<Option<String>, ToolError> {
    if count == 0 && path.is_none() {
        return Ok(None);
    }
    Ok(Some(write_or_store_layer(layer, path)?))
}

/// Reads the KML text from a `.kml` file, or the primary KML entry of a `.kmz`.
fn read_kml_text(path: &str) -> Result<String, ToolError> {
    if path.to_ascii_lowercase().ends_with(".kmz") {
        let file = std::fs::File::open(path)
            .map_err(|e| ToolError::Execution(format!("failed opening {path}: {e}")))?;
        let mut zip = zip::ZipArchive::new(file)
            .map_err(|e| ToolError::Execution(format!("failed reading KMZ archive: {e}")))?;
        // The OGC KMZ convention is a root `doc.kml`; fall back to the first
        // `.kml` entry otherwise.
        let mut chosen: Option<usize> = None;
        for i in 0..zip.len() {
            let name = zip
                .by_index(i)
                .map_err(|e| ToolError::Execution(format!("failed reading KMZ entry: {e}")))?
                .name()
                .to_string();
            let lower = name.to_ascii_lowercase();
            if lower.ends_with(".kml") {
                if lower == "doc.kml" || lower.ends_with("/doc.kml") {
                    chosen = Some(i);
                    break;
                }
                chosen.get_or_insert(i);
            }
        }
        let idx = chosen.ok_or_else(|| {
            ToolError::Execution("KMZ archive contains no .kml entry".to_string())
        })?;
        let mut entry = zip
            .by_index(idx)
            .map_err(|e| ToolError::Execution(format!("failed reading KMZ entry: {e}")))?;
        let mut text = String::new();
        entry
            .read_to_string(&mut text)
            .map_err(|e| ToolError::Execution(format!("failed decoding KML from KMZ: {e}")))?;
        Ok(text)
    } else {
        std::fs::read_to_string(path)
            .map_err(|e| ToolError::Execution(format!("failed reading {path}: {e}")))
    }
}

// ── KML parsing ──────────────────────────────────────────────────────────────

/// One parsed `<Placemark>`: its label, description, and flattened geometry.
struct Placemark {
    name: String,
    description: String,
    geoms: Vec<Geometry>,
}

/// A geometry currently being assembled from nested KML elements.
enum GBuild {
    Point(Option<Coord>),
    Line(Vec<Coord>),
    Poly {
        exterior: Vec<Coord>,
        interiors: Vec<Vec<Coord>>,
        in_inner: bool,
    },
}

impl GBuild {
    /// Finalises the in-progress geometry, dropping degenerate shapes.
    fn finish(self) -> Option<Geometry> {
        match self {
            GBuild::Point(c) => c.map(Geometry::Point),
            GBuild::Line(v) if v.len() >= 2 => Some(Geometry::LineString(v)),
            GBuild::Poly {
                exterior,
                interiors,
                ..
            } if exterior.len() >= 3 => Some(Geometry::polygon(exterior, interiors)),
            _ => None,
        }
    }
}

/// Walks the KML with a streaming pull parser, returning every `<Placemark>`.
fn parse_placemarks(xml: &str) -> Result<Vec<Placemark>, ToolError> {
    use quick_xml::events::Event;
    use quick_xml::Reader;

    let mut reader = Reader::from_str(xml);
    let mut stack: Vec<Vec<u8>> = Vec::new();
    let mut text = String::new();
    let mut pm: Option<Placemark> = None;
    let mut geom: Option<GBuild> = None;
    let mut out: Vec<Placemark> = Vec::new();

    loop {
        let ev = reader
            .read_event()
            .map_err(|e| ToolError::Execution(format!("invalid KML XML: {e}")))?;
        match ev {
            Event::Eof => break,
            Event::Start(e) => {
                let local = e.local_name().as_ref().to_vec();
                match local.as_slice() {
                    b"Placemark" => {
                        pm = Some(Placemark {
                            name: String::new(),
                            description: String::new(),
                            geoms: Vec::new(),
                        });
                        geom = None;
                    }
                    b"Point" if pm.is_some() => geom = Some(GBuild::Point(None)),
                    b"LineString" if pm.is_some() => geom = Some(GBuild::Line(Vec::new())),
                    b"Polygon" if pm.is_some() => {
                        geom = Some(GBuild::Poly {
                            exterior: Vec::new(),
                            interiors: Vec::new(),
                            in_inner: false,
                        })
                    }
                    b"innerBoundaryIs" => {
                        if let Some(GBuild::Poly { in_inner, .. }) = geom.as_mut() {
                            *in_inner = true;
                        }
                    }
                    _ => {}
                }
                stack.push(local);
                text.clear();
            }
            Event::Text(e) => {
                // The reader hands entity references as separate GeneralRef
                // events, so text spans carry no escapes to expand here.
                let raw = e
                    .decode()
                    .map_err(|e| ToolError::Execution(format!("bad KML text encoding: {e}")))?;
                text.push_str(&raw);
            }
            Event::GeneralRef(e) => {
                // Resolve `&amp;`/`&lt;`/`&#48;` etc. into their character(s).
                if let Ok(Some(ch)) = e.resolve_char_ref() {
                    text.push(ch);
                } else if let Ok(name) = e.decode() {
                    if let Some(rep) = quick_xml::escape::resolve_predefined_entity(&name) {
                        text.push_str(rep);
                    }
                }
            }
            Event::CData(e) => {
                let raw = e
                    .decode()
                    .map_err(|e| ToolError::Execution(format!("bad KML CDATA encoding: {e}")))?;
                text.push_str(&raw);
            }
            Event::End(e) => {
                let local = e.local_name().as_ref().to_vec();
                stack.pop();
                let parent = stack.last().map(Vec::as_slice);
                match local.as_slice() {
                    b"coordinates" => {
                        let coords = parse_coords(&text);
                        match geom.as_mut() {
                            Some(GBuild::Point(p)) => *p = coords.into_iter().next(),
                            Some(GBuild::Line(v)) => *v = coords,
                            Some(GBuild::Poly {
                                exterior,
                                interiors,
                                in_inner,
                            }) => {
                                if *in_inner {
                                    interiors.push(coords);
                                } else {
                                    *exterior = coords;
                                }
                            }
                            None => {}
                        }
                    }
                    b"name" if parent == Some(b"Placemark") => {
                        if let Some(pm) = pm.as_mut() {
                            pm.name = text.trim().to_string();
                        }
                    }
                    b"description" if parent == Some(b"Placemark") => {
                        if let Some(pm) = pm.as_mut() {
                            pm.description = text.trim().to_string();
                        }
                    }
                    b"innerBoundaryIs" => {
                        if let Some(GBuild::Poly { in_inner, .. }) = geom.as_mut() {
                            *in_inner = false;
                        }
                    }
                    b"Point" | b"LineString" | b"Polygon" => {
                        if let (Some(pm), Some(g)) = (pm.as_mut(), geom.take()) {
                            if let Some(finished) = g.finish() {
                                pm.geoms.push(finished);
                            }
                        }
                    }
                    b"Placemark" => {
                        if let Some(done) = pm.take() {
                            out.push(done);
                        }
                        geom = None;
                    }
                    _ => {}
                }
                text.clear();
            }
            _ => {}
        }
    }
    Ok(out)
}

/// Parses a KML `<coordinates>` string: whitespace-separated `lon,lat[,alt]`
/// tuples (longitude first, per the KML spec).
fn parse_coords(s: &str) -> Vec<Coord> {
    s.split_whitespace()
        .filter_map(|tok| {
            let mut it = tok.split(',');
            let lon: f64 = it.next()?.trim().parse().ok()?;
            let lat: f64 = it.next()?.trim().parse().ok()?;
            let alt = it.next().and_then(|a| a.trim().parse::<f64>().ok());
            Some(match alt {
                Some(z) => Coord::xyz(lon, lat, z),
                None => Coord::xy(lon, lat),
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::Layer as VLayer;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    const SAMPLE: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <name>Doc name (ignored)</name>
    <Placemark>
      <name>Alpha</name>
      <description>a point</description>
      <Point><coordinates>-74.0,40.0,12</coordinates></Point>
    </Placemark>
    <Placemark>
      <name>Beta &amp; Co</name>
      <LineString><coordinates>-74.0,40.0 -73.9,40.1 -73.8,40.2</coordinates></LineString>
    </Placemark>
    <Placemark>
      <name>Gamma</name>
      <description><![CDATA[<b>rich</b>]]></description>
      <Polygon>
        <outerBoundaryIs><LinearRing><coordinates>
          -74.0,40.0 -73.9,40.0 -73.9,40.1 -74.0,40.1 -74.0,40.0
        </coordinates></LinearRing></outerBoundaryIs>
        <innerBoundaryIs><LinearRing><coordinates>
          -73.97,40.02 -73.93,40.02 -73.93,40.06 -73.97,40.06 -73.97,40.02
        </coordinates></LinearRing></innerBoundaryIs>
      </Polygon>
    </Placemark>
  </Document>
</kml>"#;

    fn write_tmp(name: &str, body: &str) -> std::path::PathBuf {
        use std::sync::atomic::{AtomicU64, Ordering};
        static SEQ: AtomicU64 = AtomicU64::new(0);
        let n = SEQ.fetch_add(1, Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!("kml_test_{}_{n}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        std::fs::write(&path, body).unwrap();
        path
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        KmlToFeaturesTool.run(&args, &ctx()).unwrap()
    }

    fn load(path: &str) -> VLayer {
        crate::vector_common::load_input_layer(path).unwrap()
    }

    #[test]
    fn splits_placemarks_by_geometry_type() {
        let kml = write_tmp("s.kml", SAMPLE);
        let out = run(json!({ "input": kml.to_str().unwrap() }));
        assert_eq!(out.outputs["placemark_count"], json!(3));
        assert_eq!(out.outputs["point_count"], json!(1));
        assert_eq!(out.outputs["line_count"], json!(1));
        assert_eq!(out.outputs["polygon_count"], json!(1));
    }

    #[test]
    fn parses_coordinates_lon_lat_order() {
        let kml = write_tmp("s.kml", SAMPLE);
        let out = run(json!({ "input": kml.to_str().unwrap() }));
        let pts = load(out.outputs["points_output"].as_str().unwrap());
        let f = &pts.features[0];
        // KML is lon,lat[,alt]; x must be the longitude (-74), y the latitude (40).
        if let Some(Geometry::Point(c)) = &f.geometry {
            assert!((c.x - (-74.0)).abs() < 1e-9, "x should be longitude");
            assert!((c.y - 40.0).abs() < 1e-9, "y should be latitude");
            assert_eq!(c.z, Some(12.0), "altitude preserved");
        } else {
            panic!("expected a Point");
        }
        assert_eq!(
            f.get(&pts.schema, "name").unwrap(),
            &FieldValue::Text("Alpha".into())
        );
    }

    #[test]
    fn unescapes_name_and_reads_cdata_description() {
        let kml = write_tmp("s.kml", SAMPLE);
        let out = run(json!({ "input": kml.to_str().unwrap() }));
        let lines = load(out.outputs["lines_output"].as_str().unwrap());
        assert_eq!(
            lines.features[0].get(&lines.schema, "name").unwrap(),
            &FieldValue::Text("Beta & Co".into())
        );
        let polys = load(out.outputs["polygons_output"].as_str().unwrap());
        assert_eq!(
            polys.features[0].get(&polys.schema, "description").unwrap(),
            &FieldValue::Text("<b>rich</b>".into())
        );
    }

    #[test]
    fn polygon_keeps_hole() {
        let kml = write_tmp("s.kml", SAMPLE);
        let out = run(json!({ "input": kml.to_str().unwrap() }));
        let polys = load(out.outputs["polygons_output"].as_str().unwrap());
        if let Some(Geometry::Polygon { interiors, .. }) = &polys.features[0].geometry {
            assert_eq!(interiors.len(), 1, "the inner boundary should be a hole");
        } else {
            panic!("expected a Polygon");
        }
    }

    #[test]
    fn line_vertices_ordered_and_counted() {
        let kml = write_tmp("s.kml", SAMPLE);
        let out = run(json!({ "input": kml.to_str().unwrap() }));
        let lines = load(out.outputs["lines_output"].as_str().unwrap());
        if let Some(Geometry::LineString(cs)) = &lines.features[0].geometry {
            assert_eq!(cs.len(), 3);
            assert!((cs[2].x - (-73.8)).abs() < 1e-9);
        } else {
            panic!("expected a LineString");
        }
    }

    #[test]
    fn rejects_missing_input() {
        let tool = KmlToFeaturesTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "" })).is_err());
        assert!(bad(json!({ "input": "/some/file.kml" })).is_ok());
    }
}
