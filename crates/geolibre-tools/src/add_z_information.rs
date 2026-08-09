//! GeoLibre tool: append Z-derived geometry statistics to 3D features.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Add Z Information* (3D Analyst).
//!
//! ## The gap
//!
//! `vector_summary_statistics` and `attribute_histogram` summarise
//! *attributes*; `polygon_volume` and `surface_volume` measure a **raster**
//! against a reference plane. Nothing annotates a 3D **vector** feature with
//! measures derived from its own Z — which is the standard first step before
//! any 3D QA, filtering or symbology: how tall is this building, how steep is
//! this face, does this "3D" line actually carry Z at all.
//!
//! ## Noise filtering
//!
//! `noise_filtering` trims the given percentage from each end of the Z
//! distribution before `min_z` / `max_z`, matching ArcGIS. That exists because
//! a single bad vertex — a spike from a bad photogrammetric match — otherwise
//! defines the reported extremes, and the extremes are exactly what people
//! filter on.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::args_common::{f64_or, req_str};
use crate::inside_3d::collect_triangles;
use crate::mesh3d::{mesh_area, mesh_volume, topology, tri_area, tri_normal};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Every property this tool can compute, in output order.
const ALL: &[&str] = &[
    "spot_z",
    "min_z",
    "max_z",
    "mean_z",
    "length_3d",
    "surface_area",
    "volume",
    "min_slope",
    "max_slope",
    "mean_slope",
    "point_count",
    "vertex_count",
];

pub struct AddZInformationTool;

impl Tool for AddZInformationTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "add_z_information",
            display_name: "Add Z Information",
            summary: "Appends Z-derived geometry statistics to 3D features — spot/min/max/mean Z, 3D length, mesh surface area and volume, slope statistics and vertex counts (ArcGIS Add Z Information). vector_summary_statistics and attribute_histogram summarise attributes rather than geometry, and polygon_volume / surface_volume measure a raster against a reference plane, so nothing currently annotates a 3D vector feature from its own Z.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "3D point, line, polygon or multipatch features.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Input features with the requested fields appended. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "properties",
                    description: "Comma-separated subset of spot_z, min_z, max_z, mean_z, length_3d, surface_area, volume, min_slope, max_slope, mean_slope, point_count, vertex_count. Default: all of them (fields not applicable to a geometry are left null).",
                    required: false,
                },
                ToolParamSpec {
                    name: "noise_filtering",
                    description: "Percentage (0-50) trimmed from each end of the Z distribution before min_z / max_z, so a single spike vertex does not define the reported extremes. Default 0.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        parse_properties(args)?;
        let n = f64_or(args, "noise_filtering", 0.0)?;
        if !(0.0..50.0).contains(&n) {
            return Err(ToolError::Validation(
                "'noise_filtering' must be in [0, 50)".to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let props = parse_properties(args)?;
        let noise = f64_or(args, "noise_filtering", 0.0)?;

        let layer = load_input_layer(input)?;
        let mut out = Layer::new("add_z_information");
        out.geom_type = layer.geom_type;
        out.crs = layer.crs.clone();
        for f in layer.schema.fields() {
            out.add_field(f.clone());
        }
        for p in &props {
            let ty = if matches!(p.as_str(), "point_count" | "vertex_count") {
                FieldType::Integer
            } else {
                FieldType::Float
            };
            out.add_field(FieldDef::new(p.to_uppercase(), ty));
        }

        let names: Vec<String> = layer
            .schema
            .fields()
            .iter()
            .map(|f| f.name.clone())
            .collect();

        let mut with_z = 0_u64;
        let total = layer.iter().count().max(1);

        for (i, feature) in layer.iter().enumerate() {
            let stats = feature
                .geometry
                .as_ref()
                .map(|g| ZStats::of(g, noise))
                .unwrap_or_default();
            if stats.has_z {
                with_z += 1;
            }

            let mut attrs: Vec<(String, FieldValue)> = names
                .iter()
                .enumerate()
                .filter_map(|(k, n)| feature.attributes.get(k).map(|v| (n.clone(), v.clone())))
                .collect();
            for p in &props {
                attrs.push((p.to_uppercase(), stats.field(p)));
            }
            let refs: Vec<(&str, FieldValue)> =
                attrs.iter().map(|(n, v)| (n.as_str(), v.clone())).collect();
            out.add_feature(feature.geometry.clone(), &refs)
                .map_err(|e| ToolError::Execution(e.to_string()))?;
            ctx.progress.progress((i as f64 + 1.0) / total as f64);
        }

        ctx.progress.info(&format!(
            "{with_z} of {} feature(s) carry Z",
            layer.iter().count()
        ));

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("features_with_z".to_string(), json!(with_z));
        outputs.insert("properties".to_string(), json!(props));
        Ok(ToolRunResult { outputs })
    }
}

/// Z-derived measures of one geometry.
#[derive(Default)]
struct ZStats {
    has_z: bool,
    spot_z: Option<f64>,
    min_z: Option<f64>,
    max_z: Option<f64>,
    mean_z: Option<f64>,
    length_3d: Option<f64>,
    surface_area: Option<f64>,
    volume: Option<f64>,
    min_slope: Option<f64>,
    max_slope: Option<f64>,
    mean_slope: Option<f64>,
    point_count: Option<i64>,
    vertex_count: Option<i64>,
}

impl ZStats {
    fn of(geom: &Geometry, noise_percent: f64) -> Self {
        let mut s = ZStats::default();
        let verts = vertices(geom);
        s.vertex_count = Some(verts.len() as i64);
        s.point_count = Some(match geom {
            Geometry::Point(_) => 1,
            Geometry::MultiPoint(cs) => cs.len() as i64,
            _ => 0,
        });
        if verts.is_empty() {
            return s;
        }
        s.has_z = verts.iter().any(|v| v[2] != 0.0);
        s.spot_z = Some(verts[0][2]);

        let mut zs: Vec<f64> = verts.iter().map(|v| v[2]).collect();
        zs.sort_by(f64::total_cmp);
        s.mean_z = Some(zs.iter().sum::<f64>() / zs.len() as f64);
        // Trim before reading the extremes: one spike vertex would otherwise
        // define exactly the values callers filter on.
        let trim = ((noise_percent / 100.0) * zs.len() as f64).floor() as usize;
        let lo = trim.min(zs.len().saturating_sub(1));
        let hi = zs.len().saturating_sub(trim).max(lo + 1);
        s.min_z = Some(zs[lo]);
        s.max_z = Some(zs[hi - 1]);

        s.length_3d = Some(length_3d(geom));

        // Mesh measures only make sense for triangle meshes, and volume only
        // for closed ones — reporting a volume for an open surface is the
        // mistake is_closed_3d exists to catch.
        let tris = collect_triangles(geom);
        if !tris.is_empty() {
            s.surface_area = Some(mesh_area(&tris));
            let t = topology(&tris);
            if t.closed && t.consistent_winding {
                s.volume = Some(mesh_volume(&tris));
            }
            let slopes: Vec<f64> = tris
                .iter()
                .filter(|t| tri_area(t) > 0.0)
                .filter_map(|t| {
                    // Slope from the face normal: 0 degrees flat, 90 vertical.
                    tri_normal(t).map(|n| n[2].abs().clamp(-1.0, 1.0).acos().to_degrees())
                })
                .collect();
            if !slopes.is_empty() {
                s.min_slope = Some(slopes.iter().copied().fold(f64::INFINITY, f64::min));
                s.max_slope = Some(slopes.iter().copied().fold(f64::NEG_INFINITY, f64::max));
                s.mean_slope = Some(slopes.iter().sum::<f64>() / slopes.len() as f64);
            }
        }
        s
    }

    fn field(&self, name: &str) -> FieldValue {
        let f = |v: Option<f64>| v.map(FieldValue::Float).unwrap_or(FieldValue::Null);
        let i = |v: Option<i64>| v.map(FieldValue::Integer).unwrap_or(FieldValue::Null);
        match name {
            "spot_z" => f(self.spot_z),
            "min_z" => f(self.min_z),
            "max_z" => f(self.max_z),
            "mean_z" => f(self.mean_z),
            "length_3d" => f(self.length_3d),
            "surface_area" => f(self.surface_area),
            "volume" => f(self.volume),
            "min_slope" => f(self.min_slope),
            "max_slope" => f(self.max_slope),
            "mean_slope" => f(self.mean_slope),
            "point_count" => i(self.point_count),
            "vertex_count" => i(self.vertex_count),
            // `parse_properties` validates every requested name against ALL,
            // so this arm is only reachable if the two lists drift apart.
            _ => FieldValue::Null,
        }
    }
}

fn vertices(geom: &Geometry) -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    collect(geom, &mut out);
    out
}

fn collect(geom: &Geometry, out: &mut Vec<[f64; 3]>) {
    let p = |c: &Coord| [c.x, c.y, c.z.unwrap_or(0.0)];
    match geom {
        Geometry::Point(c) => out.push(p(c)),
        Geometry::MultiPoint(cs) | Geometry::LineString(cs) => out.extend(cs.iter().map(p)),
        Geometry::MultiLineString(parts) => {
            for cs in parts {
                out.extend(cs.iter().map(p));
            }
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            out.extend(exterior.0.iter().map(p));
            for r in interiors {
                out.extend(r.0.iter().map(p));
            }
        }
        Geometry::MultiPolygon(parts) => {
            for (ext, holes) in parts {
                out.extend(ext.0.iter().map(p));
                for r in holes {
                    out.extend(r.0.iter().map(p));
                }
            }
        }
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                collect(g, out);
            }
        }
    }
}

/// True 3D length: only linear geometries have one.
fn length_3d(geom: &Geometry) -> f64 {
    let run = |cs: &[Coord]| {
        cs.windows(2)
            .map(|w| {
                let dx = w[1].x - w[0].x;
                let dy = w[1].y - w[0].y;
                let dz = w[1].z.unwrap_or(0.0) - w[0].z.unwrap_or(0.0);
                (dx * dx + dy * dy + dz * dz).sqrt()
            })
            .sum::<f64>()
    };
    match geom {
        Geometry::LineString(cs) => run(cs),
        Geometry::MultiLineString(parts) => parts.iter().map(|cs| run(cs)).sum(),
        Geometry::GeometryCollection(gs) => gs.iter().map(length_3d).sum(),
        _ => 0.0,
    }
}

fn parse_properties(args: &ToolArgs) -> Result<Vec<String>, ToolError> {
    let Some(spec) = crate::args_common::opt_choice(args, "properties") else {
        return Ok(ALL.iter().map(|s| s.to_string()).collect());
    };
    let mut out = Vec::new();
    for raw in spec.split(',').map(str::trim).filter(|s| !s.is_empty()) {
        let name = ALL.iter().find(|a| **a == raw).ok_or_else(|| {
            ToolError::Validation(format!(
                "unknown property '{raw}'; expected one of {}",
                ALL.join("|")
            ))
        })?;
        if !out.iter().any(|o: &String| o == name) {
            out.push(name.to_string());
        }
    }
    if out.is_empty() {
        return Err(ToolError::Validation("'properties' is empty".to_string()));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, GeometryType};

    use crate::mesh3d::box_mesh;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn layer_of(gt: GeometryType, geoms: Vec<Geometry>) -> String {
        let mut l = Layer::new("in");
        l.geom_type = Some(gt);
        for g in geoms {
            l.add_feature(Some(g), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: Value) -> Layer {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = AddZInformationTool.run(&args, &ctx()).unwrap();
        load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap()
    }

    fn num(layer: &Layer, fid: usize, name: &str) -> f64 {
        let i = layer.schema.field_index(name).unwrap();
        match &layer.iter().nth(fid).unwrap().attributes[i] {
            FieldValue::Float(v) => *v,
            FieldValue::Integer(v) => *v as f64,
            other => panic!("expected a number in {name}, got {other:?}"),
        }
    }

    #[test]
    fn a_3d_line_reports_its_true_slant_length_not_its_plan_length() {
        // A 3-4-5 triangle in the vertical plane: plan length 3, rise 4,
        // true 3D length 5. Reporting 3 would be the classic mistake.
        let line = Geometry::LineString(vec![Coord::xyz(0.0, 0.0, 0.0), Coord::xyz(3.0, 0.0, 4.0)]);
        let out = run(json!({"input": layer_of(GeometryType::LineString, vec![line])}));
        assert!((num(&out, 0, "LENGTH_3D") - 5.0).abs() < 1e-9);
    }

    #[test]
    fn z_extremes_and_mean_come_from_the_vertices() {
        let line = Geometry::LineString(vec![
            Coord::xyz(0.0, 0.0, 10.0),
            Coord::xyz(1.0, 0.0, 20.0),
            Coord::xyz(2.0, 0.0, 30.0),
        ]);
        let out = run(json!({"input": layer_of(GeometryType::LineString, vec![line])}));
        assert!((num(&out, 0, "MIN_Z") - 10.0).abs() < 1e-9);
        assert!((num(&out, 0, "MAX_Z") - 30.0).abs() < 1e-9);
        assert!((num(&out, 0, "MEAN_Z") - 20.0).abs() < 1e-9);
        assert!((num(&out, 0, "SPOT_Z") - 10.0).abs() < 1e-9);
        assert_eq!(num(&out, 0, "VERTEX_COUNT"), 3.0);
    }

    #[test]
    fn noise_filtering_discards_a_spike_vertex() {
        // Nine vertices at z = 1..9 plus one spike at 999. Without trimming
        // MAX_Z is the spike; with it, the real surface height survives.
        let mut cs: Vec<Coord> = (1..=9)
            .map(|i| Coord::xyz(i as f64, 0.0, i as f64))
            .collect();
        cs.push(Coord::xyz(10.0, 0.0, 999.0));
        let path = layer_of(GeometryType::LineString, vec![Geometry::LineString(cs)]);

        let raw = run(json!({"input": path}));
        assert!((num(&raw, 0, "MAX_Z") - 999.0).abs() < 1e-9);

        let trimmed = run(json!({"input": path, "noise_filtering": 10.0}));
        assert!(
            (num(&trimmed, 0, "MAX_Z") - 9.0).abs() < 1e-9,
            "spike survived trimming: {}",
            num(&trimmed, 0, "MAX_Z")
        );
    }

    #[test]
    fn a_closed_mesh_reports_area_and_volume() {
        let out = run(json!({
            "input": layer_of(GeometryType::MultiPolygon, vec![box_mesh([0.0; 3], [2.0, 3.0, 4.0])]),
        }));
        assert!((num(&out, 0, "VOLUME") - 24.0).abs() < 1e-6);
        assert!((num(&out, 0, "SURFACE_AREA") - 52.0).abs() < 1e-6);
    }

    #[test]
    fn an_open_mesh_reports_area_but_no_volume() {
        let tris = collect_triangles(&box_mesh([0.0; 3], [2.0, 2.0, 2.0]));
        let open = crate::mesh3d::triangles_to_geometry(&tris[2..]);
        let out = run(json!({"input": layer_of(GeometryType::MultiPolygon, vec![open])}));
        assert!(num(&out, 0, "SURFACE_AREA") > 0.0);
        let i = out.schema.field_index("VOLUME").unwrap();
        assert_eq!(out.iter().next().unwrap().attributes[i], FieldValue::Null);
    }

    #[test]
    fn box_faces_are_either_flat_or_vertical() {
        // A box has only horizontal (0 degrees) and vertical (90 degrees)
        // faces, which pins the slope convention.
        let out = run(json!({
            "input": layer_of(GeometryType::MultiPolygon, vec![box_mesh([0.0; 3], [1.0, 1.0, 1.0])]),
        }));
        assert!(num(&out, 0, "MIN_SLOPE").abs() < 1e-6);
        assert!((num(&out, 0, "MAX_SLOPE") - 90.0).abs() < 1e-6);
    }

    #[test]
    fn a_property_subset_emits_only_those_fields() {
        let out = run(json!({
            "input": layer_of(GeometryType::MultiPolygon, vec![box_mesh([0.0; 3], [1.0, 1.0, 1.0])]),
            "properties": "min_z,max_z",
        }));
        assert!(out.schema.field_index("MIN_Z").is_some());
        assert!(out.schema.field_index("MAX_Z").is_some());
        assert!(out.schema.field_index("VOLUME").is_none());
    }

    #[test]
    fn rejects_bad_parameters() {
        let path = layer_of(
            GeometryType::MultiPolygon,
            vec![box_mesh([0.0; 3], [1.0, 1.0, 1.0])],
        );
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            AddZInformationTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"input": path, "properties": "nope"})));
        assert!(bad(json!({"input": path, "noise_filtering": 60})));
    }
}
