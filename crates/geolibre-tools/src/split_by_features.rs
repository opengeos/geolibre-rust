//! GeoLibre tool: split one layer into many outputs using a second layer's
//! field values.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Split* (Analysis). GeoLibre's
//! `split_by_attributes` splits a layer by **its own** attribute, and the
//! bundled `clip` clips against a whole layer at once. Neither performs the
//! standard "one file per county / per tile / per basin" export, which needs
//! the split *geometry* and the split *name* to come from a second layer.
//!
//! Zones that share a `split_field` value are dissolved first, so each distinct
//! value yields exactly one output file — otherwise a county represented by a
//! mainland polygon plus three islands would produce four files that overwrite
//! each other.
//!
//! Zones producing no features are skipped rather than writing empty files;
//! the count is reported so a silent miss is still visible.

use std::collections::BTreeMap;
use std::collections::HashMap;

use geo::{BooleanOps, Coord as GeoCoord, LineString, MultiLineString, MultiPolygon, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldValue, Geometry, Layer, Ring};

use crate::vector_common::{
    ensure_parent_dir, geometry_contains_point, load_input_layer, parse_optional_str,
};

pub struct SplitByFeaturesTool;

impl Tool for SplitByFeaturesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "split_by_features",
            display_name: "Split By Features",
            summary: "Clip an input layer against each zone of a polygon split layer, writing one output file per distinct split-field value. Zones sharing a value are dissolved first, and empty zones are skipped. Like ArcGIS Split (Analysis).",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Features to split (point, line or polygon).",
                    required: true,
                },
                ToolParamSpec {
                    name: "split_features",
                    description: "Polygon layer whose features define the split zones.",
                    required: true,
                },
                ToolParamSpec {
                    name: "split_field",
                    description: "Field on 'split_features' supplying each output's name.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output_dir",
                    description: "Directory receiving one output file per distinct split value.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output_format",
                    description: "Output file extension/driver: 'geojson' (default), 'shp', 'gpkg', 'parquet', 'csv'.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "split_features")?;
        require_str(args, "split_field")?;
        require_str(args, "output_dir")?;
        parse_format(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let split_path = require_str(args, "split_features")?;
        let split_field = require_str(args, "split_field")?;
        let out_dir = require_str(args, "output_dir")?.trim_end_matches('/');
        let ext = parse_format(args)?;

        let layer = load_input_layer(input)?;
        if layer.features.is_empty() {
            return Err(ToolError::Execution("input has no features".to_string()));
        }
        let splits = load_input_layer(split_path)?;
        let sidx = splits.schema.field_index(split_field).ok_or_else(|| {
            ToolError::Validation(format!(
                "split_field '{split_field}' not found in split_features"
            ))
        })?;

        // Group split polygons by value, then dissolve each group.
        let mut order: Vec<String> = Vec::new();
        let mut zones: HashMap<String, MultiPolygon> = HashMap::new();
        for f in splits.iter() {
            let Some(g) = &f.geometry else { continue };
            let Some(mp) = to_multipolygon(g) else { continue };
            let key = key_of(f.attributes.get(sidx));
            match zones.get_mut(&key) {
                Some(acc) => *acc = acc.union(&mp),
                None => {
                    order.push(key.clone());
                    zones.insert(key, mp);
                }
            }
        }
        if order.is_empty() {
            return Err(ToolError::Execution(
                "split_features contains no polygon geometry".to_string(),
            ));
        }
        ctx.progress
            .info(&format!("splitting into {} zone(s)", order.len()));

        let names = safe_file_names(&order);
        let mut written: Vec<Value> = Vec::new();
        let mut empty_zones = 0usize;
        let mut total_out = 0usize;

        for (key, name) in order.iter().zip(names.iter()) {
            let zone = &zones[key];
            // Converted once per zone: `zone_contains` used to rebuild this for
            // every coordinate tested, deep-cloning every ring per point.
            let zone_geom = multipolygon_to_geometry(zone);
            let mut out = Layer::new(name);
            if let Some(gt) = layer.geom_type {
                out = out.with_geom_type(gt);
            }
            if let Some(epsg) = layer.crs_epsg() {
                out = out.with_crs_epsg(epsg);
            }
            for fd in layer.schema.fields() {
                out.add_field(fd.clone());
            }
            let field_names: Vec<String> = layer
                .schema
                .fields()
                .iter()
                .map(|f| f.name.clone())
                .collect();

            for feat in layer.iter() {
                let Some(g) = &feat.geometry else { continue };
                let Some(clipped) = clip_geometry(g, zone, &zone_geom) else {
                    continue;
                };
                let attrs: Vec<(&str, FieldValue)> = field_names
                    .iter()
                    .enumerate()
                    .map(|(fi, nm)| {
                        (
                            nm.as_str(),
                            feat.attributes.get(fi).cloned().unwrap_or(FieldValue::Null),
                        )
                    })
                    .collect();
                out.add_feature(Some(clipped), &attrs)
                    .map_err(|e| ToolError::Execution(format!("failed adding feature: {e}")))?;
            }

            if out.features.is_empty() {
                empty_zones += 1;
                continue;
            }
            let n = out.features.len();
            total_out += n;
            let path = format!("{out_dir}/{name}.{ext}");
            ensure_parent_dir(&path)?;
            let fmt = wbvector::VectorFormat::detect(&path)
                .map_err(|e| ToolError::Validation(format!("unsupported output format: {e}")))?;
            wbvector::write(&out, &path, fmt)
                .map_err(|e| ToolError::Execution(format!("failed writing '{path}': {e}")))?;
            written.push(json!({ "split_value": key, "path": path, "feature_count": n }));
        }

        let mut outputs = BTreeMap::new();
        // 'output' names the directory so generic callers have a single handle.
        outputs.insert("output".to_string(), json!(out_dir));
        outputs.insert("output_dir".to_string(), json!(out_dir));
        outputs.insert("file_count".to_string(), json!(written.len()));
        outputs.insert("zone_count".to_string(), json!(order.len()));
        outputs.insert("empty_zone_count".to_string(), json!(empty_zones));
        outputs.insert("output_feature_count".to_string(), json!(total_out));
        outputs.insert("files".to_string(), json!(written));
        Ok(ToolRunResult { outputs })
    }
}

/// Clips one geometry to `zone`, returning `None` when nothing survives.
///
/// `zone_geom` is the same zone pre-converted to a `wbvector::Geometry` for the
/// point-containment test, so it is built once per zone rather than per point.
fn clip_geometry(g: &Geometry, zone: &MultiPolygon, zone_geom: &Geometry) -> Option<Geometry> {
    match g {
        // Points are kept whole or dropped — clipping cannot subdivide them.
        Geometry::Point(c) => geometry_contains_point(zone_geom, c.x, c.y).then(|| g.clone()),
        Geometry::MultiPoint(cs) => {
            let kept: Vec<Coord> = cs
                .iter()
                .filter(|c| geometry_contains_point(zone_geom, c.x, c.y))
                .cloned()
                .collect();
            match kept.len() {
                0 => None,
                1 => Some(Geometry::Point(kept.into_iter().next()?)),
                _ => Some(Geometry::MultiPoint(kept)),
            }
        }
        Geometry::LineString(cs) => {
            clip_lines(zone, &MultiLineString(vec![coords_to_linestring(cs)]))
        }
        Geometry::MultiLineString(ls) => clip_lines(
            zone,
            &MultiLineString(ls.iter().map(|l| coords_to_linestring(l)).collect()),
        ),
        Geometry::Polygon { .. } | Geometry::MultiPolygon(_) => {
            let mp = to_multipolygon(g)?;
            let cut = zone.intersection(&mp);
            (!cut.0.is_empty()).then(|| multipolygon_to_geometry(&cut))
        }
        Geometry::GeometryCollection(gs) => {
            let kept: Vec<Geometry> = gs
                .iter()
                .filter_map(|s| clip_geometry(s, zone, zone_geom))
                .collect();
            (!kept.is_empty()).then_some(Geometry::GeometryCollection(kept))
        }
    }
}

fn clip_lines(zone: &MultiPolygon, mls: &MultiLineString) -> Option<Geometry> {
    let cut = zone.clip(mls, false);
    let parts: Vec<Vec<Coord>> = cut
        .0
        .iter()
        .filter(|l| l.0.len() >= 2)
        .map(|l| l.0.iter().map(|c| Coord::xy(c.x, c.y)).collect())
        .collect();
    match parts.len() {
        0 => None,
        1 => Some(Geometry::LineString(parts.into_iter().next()?)),
        _ => Some(Geometry::MultiLineString(parts)),
    }
}

/// Turns split values into safe, unique file stems.
///
/// Every emitted stem is recorded, not just the sanitized base: keying on the
/// base alone lets ["west", "west", "west_1"] emit `west_1` twice, silently
/// overwriting the first file — exactly what the dissolve step exists to avoid.
fn safe_file_names(values: &[String]) -> Vec<String> {
    let mut used: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut out = Vec::with_capacity(values.len());
    for v in values {
        let mut base: String = v
            .chars()
            .map(|c| {
                if c.is_ascii_alphanumeric() || c == '-' || c == '_' {
                    c
                } else {
                    '_'
                }
            })
            .collect();
        if base.is_empty() {
            base = "zone".to_string();
        }
        base.truncate(80);
        let name = if used.insert(base.clone()) {
            base
        } else {
            let candidate = (1..)
                .map(|k| format!("{base}_{k}"))
                .find(|c| !used.contains(c))
                .expect("an unused suffix always exists");
            used.insert(candidate.clone());
            candidate
        };
        out.push(name);
    }
    out
}

fn key_of(v: Option<&FieldValue>) -> String {
    match v {
        None | Some(FieldValue::Null) => "NULL".to_string(),
        Some(FieldValue::Integer(i)) => i.to_string(),
        Some(FieldValue::Float(f)) => {
            if f.fract() == 0.0 && f.is_finite() && f.abs() < 1e15 {
                format!("{}", *f as i64)
            } else {
                format!("{f}")
            }
        }
        Some(FieldValue::Text(s)) => s.clone(),
        Some(FieldValue::Boolean(b)) => b.to_string(),
        Some(FieldValue::Date(s)) | Some(FieldValue::DateTime(s)) => s.clone(),
        Some(FieldValue::Blob(b)) => format!("blob[{}]", b.len()),
    }
}

// ── geo <-> wbvector conversion ─────────────────────────────────────────────

fn to_multipolygon(geom: &Geometry) -> Option<MultiPolygon> {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => Some(MultiPolygon(vec![rings_to_polygon(exterior, interiors)])),
        Geometry::MultiPolygon(parts) => Some(MultiPolygon(
            parts
                .iter()
                .map(|(ext, ints)| rings_to_polygon(ext, ints))
                .collect(),
        )),
        _ => None,
    }
}

fn rings_to_polygon(exterior: &Ring, interiors: &[Ring]) -> Polygon {
    Polygon::new(
        ring_to_linestring(exterior),
        interiors.iter().map(ring_to_linestring).collect(),
    )
}

fn ring_to_linestring(ring: &Ring) -> LineString {
    LineString::new(
        ring.coords()
            .iter()
            .map(|c| GeoCoord { x: c.x, y: c.y })
            .collect(),
    )
}

fn coords_to_linestring(cs: &[Coord]) -> LineString {
    LineString::new(cs.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect())
}

fn multipolygon_to_geometry(mp: &MultiPolygon) -> Geometry {
    if mp.0.len() == 1 {
        let (exterior, interiors) = polygon_to_rings(&mp.0[0]);
        Geometry::Polygon {
            exterior,
            interiors,
        }
    } else {
        Geometry::MultiPolygon(mp.0.iter().map(polygon_to_rings).collect())
    }
}

fn polygon_to_rings(poly: &Polygon) -> (Ring, Vec<Ring>) {
    (
        linestring_to_ring(poly.exterior()),
        poly.interiors().iter().map(linestring_to_ring).collect(),
    )
}

fn linestring_to_ring(ls: &LineString) -> Ring {
    let mut coords: Vec<Coord> = ls.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
    if coords.len() >= 2 && coords.first().map(|c| (c.x, c.y)) == coords.last().map(|c| (c.x, c.y))
    {
        coords.pop();
    }
    Ring::new(coords)
}

// ── Params ──────────────────────────────────────────────────────────────────

fn parse_format(args: &ToolArgs) -> Result<String, ToolError> {
    let raw = parse_optional_str(args, "output_format")?
        .unwrap_or("geojson")
        .trim()
        .trim_start_matches('.')
        .to_ascii_lowercase();
    match raw.as_str() {
        "geojson" | "json" | "shp" | "gpkg" | "parquet" | "csv" => Ok(raw),
        o => Err(ToolError::Validation(format!(
            "'output_format' must be one of geojson/shp/gpkg/parquet/csv, got '{o}'"
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
    use wbvector::{FieldDef, FieldType, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    // Test layers use EPSG:4326 on purpose: the GeoJSON writer reprojects to
    // WGS84 per the format spec, so tagging the fixtures 3857 would make the
    // area/length assertions below compare metres against degrees.
    fn tmpdir(tag: &str) -> String {
        let d = std::env::temp_dir().join(format!("geolibre_split_{tag}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d.to_string_lossy().to_string()
    }

    fn square(x0: f64, y0: f64, w: f64) -> Geometry {
        Geometry::polygon(
            vec![
                Coord::xy(x0, y0),
                Coord::xy(x0 + w, y0),
                Coord::xy(x0 + w, y0 + w),
                Coord::xy(x0, y0 + w),
                Coord::xy(x0, y0),
            ],
            vec![],
        )
    }

    /// Two 10x10 zones: "west" at x 0-10, "east" at x 10-20.
    fn zones() -> String {
        let mut l = Layer::new("zones")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("zone", FieldType::Text));
        for (name, x0) in [("west", 0.0), ("east", 10.0)] {
            l.add_feature(
                Some(square(x0, 0.0, 10.0)),
                &[("zone", FieldValue::Text(name.to_string()))],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn points(p: &[(f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("id", FieldType::Integer));
        for (i, (x, y)) in p.iter().enumerate() {
            l.add_feature(
                Some(Geometry::point(*x, *y)),
                &[("id", FieldValue::Integer(i as i64))],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        SplitByFeaturesTool.run(&args, &ctx()).unwrap()
    }

    #[test]
    fn points_are_partitioned_between_zones() {
        let dir = tmpdir("pts");
        let out = run(json!({
            "input": points(&[(1.0, 1.0), (2.0, 2.0), (15.0, 5.0)]),
            "split_features": zones(), "split_field": "zone", "output_dir": dir,
        }));
        assert_eq!(out.outputs["file_count"], json!(2));
        // 2 west + 1 east = 3, i.e. every input landed in exactly one file.
        assert_eq!(out.outputs["output_feature_count"], json!(3));
        let west = load_input_layer(&format!("{dir}/west.geojson")).unwrap();
        let east = load_input_layer(&format!("{dir}/east.geojson")).unwrap();
        assert_eq!(west.features.len(), 2);
        assert_eq!(east.features.len(), 1);
    }

    #[test]
    fn polygons_are_cut_at_the_zone_boundary() {
        // One 20-wide polygon spanning both zones must appear, halved, in each.
        let mut l = Layer::new("big")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        l.add_feature(Some(square(0.0, 0.0, 20.0)), &[]).unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let input = wbvector::memory_store::make_vector_memory_path(&id);

        let dir = tmpdir("poly");
        let out = run(json!({
            "input": input, "split_features": zones(), "split_field": "zone",
            "output_dir": dir,
        }));
        assert_eq!(out.outputs["file_count"], json!(2));
        let area = |p: &str| -> f64 {
            use geo::Area;
            let l = load_input_layer(p).unwrap();
            l.iter()
                .filter_map(|f| f.geometry.as_ref().and_then(to_multipolygon))
                .map(|mp| mp.unsigned_area())
                .sum()
        };
        let w = area(&format!("{dir}/west.geojson"));
        let e = area(&format!("{dir}/east.geojson"));
        // Each zone is 10x10 of the 20x20 input's lower half -> 100 each.
        assert!((w - 100.0).abs() < 1e-6, "west area {w}");
        assert!((e - 100.0).abs() < 1e-6, "east area {e}");
    }

    #[test]
    fn lines_are_clipped_to_each_zone() {
        let mut l = Layer::new("line")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(4326);
        l.add_feature(
            Some(Geometry::line_string(vec![
                Coord::xy(0.0, 5.0),
                Coord::xy(20.0, 5.0),
            ])),
            &[],
        )
        .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let input = wbvector::memory_store::make_vector_memory_path(&id);

        let dir = tmpdir("line");
        run(json!({
            "input": input, "split_features": zones(), "split_field": "zone",
            "output_dir": dir,
        }));
        let length = |p: &str| -> f64 {
            let l = load_input_layer(p).unwrap();
            l.iter()
                .filter_map(|f| f.geometry.clone())
                .map(|g| {
                    g.all_coords()
                        .windows(2)
                        .map(|w| (w[1].x - w[0].x).hypot(w[1].y - w[0].y))
                        .sum::<f64>()
                })
                .sum()
        };
        assert!((length(&format!("{dir}/west.geojson")) - 10.0).abs() < 1e-6);
        assert!((length(&format!("{dir}/east.geojson")) - 10.0).abs() < 1e-6);
    }

    #[test]
    fn zones_sharing_a_value_are_dissolved_into_one_file() {
        // Two disjoint polygons both labelled "west" -> a single west file
        // holding features from both parts, not two files clobbering each other.
        let mut l = Layer::new("zones")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("zone", FieldType::Text));
        for x0 in [0.0, 30.0] {
            l.add_feature(
                Some(square(x0, 0.0, 10.0)),
                &[("zone", FieldValue::Text("west".to_string()))],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        let split = wbvector::memory_store::make_vector_memory_path(&id);

        let dir = tmpdir("dissolve");
        let out = run(json!({
            "input": points(&[(1.0, 1.0), (31.0, 1.0)]),
            "split_features": split, "split_field": "zone", "output_dir": dir,
        }));
        assert_eq!(out.outputs["zone_count"], json!(1));
        assert_eq!(out.outputs["file_count"], json!(1));
        let west = load_input_layer(&format!("{dir}/west.geojson")).unwrap();
        assert_eq!(west.features.len(), 2);
    }

    #[test]
    fn empty_zones_are_skipped_not_written() {
        let dir = tmpdir("empty");
        let out = run(json!({
            "input": points(&[(1.0, 1.0)]),
            "split_features": zones(), "split_field": "zone", "output_dir": dir,
        }));
        assert_eq!(out.outputs["zone_count"], json!(2));
        assert_eq!(out.outputs["file_count"], json!(1));
        assert_eq!(out.outputs["empty_zone_count"], json!(1));
        assert!(!std::path::Path::new(&format!("{dir}/east.geojson")).exists());
    }

    #[test]
    fn attributes_are_preserved_on_the_split_output() {
        let dir = tmpdir("attrs");
        run(json!({
            "input": points(&[(1.0, 1.0)]),
            "split_features": zones(), "split_field": "zone", "output_dir": dir,
        }));
        let west = load_input_layer(&format!("{dir}/west.geojson")).unwrap();
        let i = west.schema.field_index("id").unwrap();
        assert_eq!(west.features[0].attributes[i].as_i64(), Some(0));
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SplitByFeaturesTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.shp", "split_features": "z.shp" })).is_err());
        assert!(bad(json!({
            "input": "a.shp", "split_features": "z.shp", "split_field": "z",
            "output_dir": "/tmp/x", "output_format": "dxf"
        }))
        .is_err());
        assert!(bad(json!({
            "input": "a.shp", "split_features": "z.shp", "split_field": "z",
            "output_dir": "/tmp/x"
        }))
        .is_ok());
    }
}
