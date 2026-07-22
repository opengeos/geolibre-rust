//! GeoLibre tool: write point/line GIS features back out as GTFS text files.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Features To GTFS Stops* /
//! *Features To GTFS Shapes* (Public Transit): the reverse of
//! [`gtfs_to_features`](crate::gtfs_to_features). Point features become
//! `stops.txt` (`stop_id`, `stop_name`, `stop_lat`, `stop_lon`); line features
//! become `shapes.txt` (each LineString exploded into ordered
//! `shape_pt_lat`/`shape_pt_lon`/`shape_pt_sequence` rows with a cumulative
//! `shape_dist_traveled` in metres). Together the two tools close the GTFS
//! editor round-trip a transit workflow needs.
//!
//! `stops_input`/`shapes_input` are point/line vector layers (at least one is
//! required); `output_dir` is the folder the `*.txt` files are written to.
//! GTFS coordinates are always WGS84, so inputs are assumed EPSG:4326 (as
//! `gtfs_to_features` emits) and coordinates are written straight through.
//!
//! When an input carries `stop_id`/`stop_name`/`shape_id` attribute fields
//! (exactly what `gtfs_to_features` produces) they are reused, so a feed can be
//! imported, edited, and exported without losing its identifiers; otherwise
//! stable sequential ids are generated.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Geometry, Layer};

use crate::vector_common::{ensure_parent_dir, load_input_layer, parse_optional_str};

pub struct FeaturesToGtfsTool;

impl Tool for FeaturesToGtfsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "features_to_gtfs",
            display_name: "Features To GTFS",
            summary: "Export point features to GTFS stops.txt and line features to GTFS shapes.txt, closing the GTFS round-trip.",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "stops_input",
                    description: "Point vector layer to write as stops.txt (uses stop_id/stop_name attribute fields when present).",
                    required: false,
                },
                ToolParamSpec {
                    name: "shapes_input",
                    description: "Line vector layer to write as shapes.txt (uses a shape_id attribute field when present).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_dir",
                    description: "Directory the GTFS text files (stops.txt / shapes.txt) are written to. Created if missing.",
                    required: true,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if str_arg(args, "output_dir").is_empty() {
            return Err(ToolError::Validation(
                "missing required string parameter 'output_dir'".to_string(),
            ));
        }
        if str_arg(args, "stops_input").is_empty() && str_arg(args, "shapes_input").is_empty() {
            return Err(ToolError::Validation(
                "at least one of 'stops_input' or 'shapes_input' must be provided".to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let output_dir = parse_optional_str(args, "output_dir")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'output_dir'".to_string())
        })?;
        let stops_input = parse_optional_str(args, "stops_input")?;
        let shapes_input = parse_optional_str(args, "shapes_input")?;
        if stops_input.is_none() && shapes_input.is_none() {
            return Err(ToolError::Validation(
                "at least one of 'stops_input' or 'shapes_input' must be provided".to_string(),
            ));
        }

        std::fs::create_dir_all(output_dir)
            .map_err(|e| ToolError::Execution(format!("failed creating output directory: {e}")))?;

        let mut outputs = BTreeMap::new();

        // ── stops.txt ───────────────────────────────────────────────────────
        if let Some(path) = stops_input {
            let layer = load_input_layer(path)?;
            let (csv, count) = build_stops(&layer)?;
            let out = join_dir(output_dir, "stops.txt");
            ensure_parent_dir(&out)?;
            std::fs::write(&out, csv)
                .map_err(|e| ToolError::Execution(format!("failed writing stops.txt: {e}")))?;
            ctx.progress
                .info(&format!("wrote {count} stop(s) to stops.txt"));
            outputs.insert("stops_output".to_string(), json!(out));
            outputs.insert("stop_count".to_string(), json!(count));
        }

        // ── shapes.txt ──────────────────────────────────────────────────────
        if let Some(path) = shapes_input {
            let layer = load_input_layer(path)?;
            let (csv, shape_count, point_count) = build_shapes(&layer)?;
            let out = join_dir(output_dir, "shapes.txt");
            ensure_parent_dir(&out)?;
            std::fs::write(&out, csv)
                .map_err(|e| ToolError::Execution(format!("failed writing shapes.txt: {e}")))?;
            ctx.progress.info(&format!(
                "wrote {shape_count} shape(s), {point_count} vertices to shapes.txt"
            ));
            outputs.insert("shapes_output".to_string(), json!(out));
            outputs.insert("shape_count".to_string(), json!(shape_count));
            outputs.insert("shape_point_count".to_string(), json!(point_count));
        }

        // "output" points at the primary artifact (stops if present, else shapes).
        let primary = outputs
            .get("stops_output")
            .or_else(|| outputs.get("shapes_output"))
            .cloned();
        if let Some(p) = primary {
            outputs.insert("output".to_string(), p);
        }
        outputs.insert("output_dir".to_string(), json!(output_dir));

        Ok(ToolRunResult { outputs })
    }
}

/// Emits `stops.txt` from a point layer. Every point (and every part of a
/// MultiPoint) becomes one stop row.
fn build_stops(layer: &Layer) -> Result<(String, usize), ToolError> {
    let has_id = layer.schema.field_index("stop_id").is_some();
    let has_name = layer.schema.field_index("stop_name").is_some();

    let mut csv = String::from("stop_id,stop_name,stop_lat,stop_lon\n");
    let mut count = 0usize;
    let mut auto = 0usize;
    for feat in &layer.features {
        let Some(geom) = &feat.geometry else { continue };
        let pts = point_coords(geom);
        if pts.is_empty() {
            continue;
        }
        let base_id = if has_id {
            feat.get(&layer.schema, "stop_id")
                .ok()
                .map(field_to_string)
                .filter(|s| !s.is_empty())
        } else {
            None
        };
        let name = if has_name {
            feat.get(&layer.schema, "stop_name")
                .ok()
                .map(field_to_string)
                .unwrap_or_default()
        } else {
            String::new()
        };
        for (i, c) in pts.iter().enumerate() {
            // A single-point feature keeps its id; multipoint parts get a suffix.
            let id = match &base_id {
                Some(id) if pts.len() == 1 => id.clone(),
                Some(id) => format!("{id}_{}", i + 1),
                None => {
                    auto += 1;
                    format!("S{auto}")
                }
            };
            csv.push_str(&csv_field(&id));
            csv.push(',');
            csv.push_str(&csv_field(&name));
            csv.push(',');
            csv.push_str(&fmt_coord(c.y));
            csv.push(',');
            csv.push_str(&fmt_coord(c.x));
            csv.push('\n');
            count += 1;
        }
    }
    Ok((csv, count))
}

/// Emits `shapes.txt` from a line layer. Each LineString (and each part of a
/// MultiLineString) becomes one shape, exploded into ordered vertex rows with a
/// cumulative haversine `shape_dist_traveled` in metres.
fn build_shapes(layer: &Layer) -> Result<(String, usize, usize), ToolError> {
    let has_id = layer.schema.field_index("shape_id").is_some();

    let mut csv =
        String::from("shape_id,shape_pt_lat,shape_pt_lon,shape_pt_sequence,shape_dist_traveled\n");
    let mut shape_count = 0usize;
    let mut point_count = 0usize;
    let mut auto = 0usize;
    for feat in &layer.features {
        let Some(geom) = &feat.geometry else { continue };
        let lines = line_parts(geom);
        if lines.is_empty() {
            continue;
        }
        let base_id = if has_id {
            feat.get(&layer.schema, "shape_id")
                .ok()
                .map(field_to_string)
                .filter(|s| !s.is_empty())
        } else {
            None
        };
        for (li, line) in lines.iter().enumerate() {
            if line.len() < 2 {
                continue;
            }
            let id = match &base_id {
                Some(id) if lines.len() == 1 => id.clone(),
                Some(id) => format!("{id}_{}", li + 1),
                None => {
                    auto += 1;
                    format!("shp{auto}")
                }
            };
            let mut dist = 0.0f64;
            let mut prev: Option<&Coord> = None;
            for (seq, c) in line.iter().enumerate() {
                if let Some(p) = prev {
                    dist += haversine_m(p.y, p.x, c.y, c.x);
                }
                csv.push_str(&csv_field(&id));
                csv.push(',');
                csv.push_str(&fmt_coord(c.y));
                csv.push(',');
                csv.push_str(&fmt_coord(c.x));
                csv.push(',');
                csv.push_str(&(seq + 1).to_string());
                csv.push(',');
                csv.push_str(&fmt_dist(dist));
                csv.push('\n');
                prev = Some(c);
                point_count += 1;
            }
            shape_count += 1;
        }
    }
    Ok((csv, shape_count, point_count))
}

/// Extracts every point coordinate from a Point / MultiPoint geometry.
fn point_coords(geom: &Geometry) -> Vec<Coord> {
    match geom {
        Geometry::Point(c) => vec![c.clone()],
        Geometry::MultiPoint(cs) => cs.clone(),
        _ => Vec::new(),
    }
}

/// Extracts every line part from a LineString / MultiLineString geometry.
fn line_parts(geom: &Geometry) -> Vec<Vec<Coord>> {
    match geom {
        Geometry::LineString(cs) => vec![cs.clone()],
        Geometry::MultiLineString(ls) => ls.clone(),
        _ => Vec::new(),
    }
}

/// Great-circle distance between two WGS84 lat/lon points, in metres.
fn haversine_m(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> f64 {
    const R: f64 = 6_371_000.0;
    let (p1, p2) = (lat1.to_radians(), lat2.to_radians());
    let dphi = (lat2 - lat1).to_radians();
    let dlam = (lon2 - lon1).to_radians();
    let a = (dphi / 2.0).sin().powi(2) + p1.cos() * p2.cos() * (dlam / 2.0).sin().powi(2);
    2.0 * R * a.sqrt().asin()
}

/// Renders a coordinate with enough precision for GTFS (~1 cm) and no trailing
/// zeros, matching how GTFS feeds commonly write lat/lon.
fn fmt_coord(v: f64) -> String {
    trim_num(format!("{v:.7}"))
}

/// Renders a cumulative distance in metres (2 decimals), trimmed.
fn fmt_dist(v: f64) -> String {
    trim_num(format!("{v:.2}"))
}

fn trim_num(mut s: String) -> String {
    if s.contains('.') {
        while s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
    }
    s
}

/// Renders any attribute value as GTFS text.
fn field_to_string(v: &wbvector::FieldValue) -> String {
    use wbvector::FieldValue::*;
    match v {
        Null => String::new(),
        Text(s) | Date(s) | DateTime(s) => s.clone(),
        Integer(i) => i.to_string(),
        Float(f) => trim_num(format!("{f}")),
        Boolean(b) => b.to_string(),
        Blob(_) => String::new(),
    }
}

/// CSV-escapes a field: quote when it contains a comma, quote, or newline.
fn csv_field(s: &str) -> String {
    if s.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// Joins a directory and a file name with a single separator.
fn join_dir(dir: &str, file: &str) -> String {
    std::path::Path::new(dir)
        .join(file)
        .to_string_lossy()
        .into_owned()
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

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        FeaturesToGtfsTool.run(&args, &ctx()).unwrap()
    }

    fn tmp_dir(tag: &str) -> std::path::PathBuf {
        use std::sync::atomic::{AtomicU64, Ordering};
        static SEQ: AtomicU64 = AtomicU64::new(0);
        let n = SEQ.fetch_add(1, Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!("ftg_{tag}_{}_{n}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn stops_layer() -> String {
        let mut l = Layer::new("stops")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("stop_id", FieldType::Text));
        l.add_field(FieldDef::new("stop_name", FieldType::Text));
        l.add_feature(
            Some(Geometry::point(-74.0, 40.0)),
            &[
                ("stop_id", "A".into()),
                ("stop_name", "Alpha, Main St".into()),
            ],
        )
        .unwrap();
        l.add_feature(
            Some(Geometry::point(-74.1, 40.1)),
            &[("stop_id", "B".into()), ("stop_name", "Beta".into())],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn shapes_layer() -> String {
        let mut l = Layer::new("shapes")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("shape_id", FieldType::Text));
        l.add_feature(
            Some(Geometry::LineString(vec![
                Coord::xy(-74.0, 40.0),
                Coord::xy(-74.05, 40.05),
                Coord::xy(-74.1, 40.1),
            ])),
            &[("shape_id", "S1".into())],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    #[test]
    fn writes_stops_txt() {
        let dir = tmp_dir("stops");
        let out = run(json!({
            "stops_input": stops_layer(),
            "output_dir": dir.to_str().unwrap(),
        }));
        assert_eq!(out.outputs["stop_count"], json!(2));
        let txt = std::fs::read_to_string(dir.join("stops.txt")).unwrap();
        let lines: Vec<&str> = txt.lines().collect();
        assert_eq!(lines[0], "stop_id,stop_name,stop_lat,stop_lon");
        // Attribute stop_id reused; lat then lon; comma-bearing name quoted.
        assert_eq!(lines[1], "A,\"Alpha, Main St\",40,-74");
        assert_eq!(lines[2], "B,Beta,40.1,-74.1");
    }

    #[test]
    fn writes_shapes_txt_ordered_with_distance() {
        let dir = tmp_dir("shapes");
        let out = run(json!({
            "shapes_input": shapes_layer(),
            "output_dir": dir.to_str().unwrap(),
        }));
        assert_eq!(out.outputs["shape_count"], json!(1));
        assert_eq!(out.outputs["shape_point_count"], json!(3));
        let txt = std::fs::read_to_string(dir.join("shapes.txt")).unwrap();
        let lines: Vec<&str> = txt.lines().collect();
        assert_eq!(
            lines[0],
            "shape_id,shape_pt_lat,shape_pt_lon,shape_pt_sequence,shape_dist_traveled"
        );
        // First vertex: sequence 1, zero cumulative distance.
        assert!(lines[1].starts_with("S1,40,-74,1,0"));
        // Sequences ascend 1,2,3.
        assert!(lines[2].contains(",2,"));
        assert!(lines[3].contains(",3,"));
        // Cumulative distance strictly increases along the shape.
        let d2: f64 = lines[2].rsplit(',').next().unwrap().parse().unwrap();
        let d3: f64 = lines[3].rsplit(',').next().unwrap().parse().unwrap();
        assert!(d2 > 0.0 && d3 > d2, "distance must accumulate: {d2} {d3}");
    }

    #[test]
    fn both_inputs_write_both_files() {
        let dir = tmp_dir("both");
        run(json!({
            "stops_input": stops_layer(),
            "shapes_input": shapes_layer(),
            "output_dir": dir.to_str().unwrap(),
        }));
        assert!(dir.join("stops.txt").exists());
        assert!(dir.join("shapes.txt").exists());
    }

    #[test]
    fn generates_ids_when_absent() {
        // A layer with no stop_id/stop_name fields still exports valid rows.
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        l.add_feature(Some(Geometry::point(-100.0, 30.0)), &[])
            .unwrap();
        l.add_feature(Some(Geometry::point(-101.0, 31.0)), &[])
            .unwrap();
        let id = memory_store::put_vector(l);
        let path = memory_store::make_vector_memory_path(&id);
        let dir = tmp_dir("auto");
        run(json!({ "stops_input": path, "output_dir": dir.to_str().unwrap() }));
        let txt = std::fs::read_to_string(dir.join("stops.txt")).unwrap();
        let lines: Vec<&str> = txt.lines().collect();
        assert_eq!(lines[1], "S1,,30,-100");
        assert_eq!(lines[2], "S2,,31,-101");
    }

    #[test]
    fn non_matching_geometry_passes_through() {
        // A polygon layer handed to shapes_input yields an empty shapes.txt
        // (header only) rather than erroring.
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
        let dir = tmp_dir("poly");
        let out = run(json!({ "shapes_input": path, "output_dir": dir.to_str().unwrap() }));
        assert_eq!(out.outputs["shape_count"], json!(0));
        let txt = std::fs::read_to_string(dir.join("shapes.txt")).unwrap();
        assert_eq!(txt.lines().count(), 1); // header only
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = FeaturesToGtfsTool;
        let v = |val: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(val).unwrap();
            tool.validate(&args)
        };
        // Missing output_dir.
        assert!(v(json!({ "stops_input": "x.geojson" })).is_err());
        // Missing both inputs.
        assert!(v(json!({ "output_dir": "/tmp/x" })).is_err());
        // Valid.
        assert!(v(json!({ "stops_input": "x.geojson", "output_dir": "/tmp/x" })).is_ok());
    }
}
