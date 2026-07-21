//! GeoLibre tool: vector-on-vector zonal summary (tabulate intersection).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Tabulate Intersection* and the core
//! of *Summarize Within* (Analysis): apportion one layer's features across
//! another layer's zones and report how much of each class falls in each zone.
//! Whitebox-wasm's `zonal_statistics` / `cross_tabulation` are raster-only, and
//! a plain spatial join transfers attributes without area-weighting — so "how
//! much of each land-cover class, or how many points, fall inside each district"
//! has had no vector-native answer.
//!
//! For polygon class features the measure is intersected **area** (via `geo`
//! `BooleanOps`); for point class features it is the **count** of points inside
//! the zone. Results are grouped by zone × `class_field` value (or a single
//! "ALL" class when no field is given). Each output row carries:
//!
//! - `zone_id` — the zone identifier (`zone_field` value, or the zone's index);
//! - the class value (under the `class_field` name), when grouping;
//! - `area` or `count` — the intersected measure;
//! - `percentage` — the measure as a percent of the zone's total across all
//!   classes (so a zone's rows sum to 100);
//! - one column per `sum_fields` entry — the class attribute apportioned by the
//!   intersected fraction (`value x intersected_area / class_area` for polygons;
//!   the plain sum over contained points for points).
//!
//! Output geometry is the intersection polygon per zone × class (polygon
//! classes) or the zone polygon (point classes), so the table is also mappable.
//!
//! Scope for v1: line class features (intersected length) are not yet supported
//! — `geo` has no arbitrary line-in-polygon clip; use polygon or point classes.

use std::collections::BTreeMap;

use geo::{
    Area, BooleanOps, Contains, Coord as GeoCoord, LineString, MultiPolygon, Point, Polygon,
};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct TabulateIntersectionTool;

impl Tool for TabulateIntersectionTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "tabulate_intersection",
            display_name: "Tabulate Intersection",
            summary: "Vector-on-vector zonal summary: apportion a class layer (polygons by area, or points by count) across zone polygons, reporting intersected measure, percentage of zone, and area-weighted sum fields.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Zone polygon layer, format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "class_features",
                    description: "Class vector layer to apportion across the zones (polygons summarized by area, or points by count).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "class_field",
                    description: "Optional categorical field in the class layer; results are grouped by zone x this value. Omit to treat all class features as one class.",
                    required: false,
                },
                ToolParamSpec {
                    name: "sum_fields",
                    description: "Optional comma-separated numeric fields in the class layer to apportion into each zone (area-weighted for polygons, summed for points).",
                    required: false,
                },
                ToolParamSpec {
                    name: "zone_field",
                    description: "Optional field identifying each zone in the output. Defaults to the zone feature index.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["input", "class_features"] {
            if args
                .get(key)
                .and_then(Value::as_str)
                .map(str::trim)
                .unwrap_or("")
                .is_empty()
            {
                return Err(ToolError::Validation(format!(
                    "missing required string parameter '{key}'"
                )));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let class_path = require_str(args, "class_features")?;
        let output = parse_optional_str(args, "output")?;
        let class_field = parse_optional_str(args, "class_field")?.map(str::to_string);
        let zone_field = parse_optional_str(args, "zone_field")?.map(str::to_string);
        let sum_fields: Vec<String> = parse_optional_str(args, "sum_fields")?
            .map(|s| {
                s.split(',')
                    .map(str::trim)
                    .filter(|t| !t.is_empty())
                    .map(str::to_string)
                    .collect()
            })
            .unwrap_or_default();

        let zones = load_input_layer(input)?;
        let classes = load_input_layer(class_path)?;
        let class_kind = ClassKind::detect(&classes).ok_or_else(|| {
            ToolError::Validation(
                "class layer has no polygon or point features (line classes are not supported)"
                    .to_string(),
            )
        })?;

        // Pre-extract class features as (geo geometry, class value, sum values).
        let class_schema = classes.schema.clone();
        struct ClassFeat {
            poly: Option<MultiPolygon>,
            point: Option<Point>,
            value: String,
            sums: Vec<f64>,
        }
        let class_feats: Vec<ClassFeat> = classes
            .features
            .iter()
            .filter_map(|feat| {
                let geom = feat.geometry.as_ref()?;
                let value = match &class_field {
                    Some(f) => feat
                        .get(&class_schema, f)
                        .map(field_value_string)
                        .unwrap_or_default(),
                    None => "ALL".to_string(),
                };
                let sums = sum_fields
                    .iter()
                    .map(|f| {
                        feat.get(&class_schema, f)
                            .ok()
                            .and_then(FieldValue::as_f64)
                            .unwrap_or(0.0)
                    })
                    .collect();
                match class_kind {
                    ClassKind::Polygon => to_multipolygon(geom).map(|poly| ClassFeat {
                        poly: Some(poly),
                        point: None,
                        value,
                        sums,
                    }),
                    ClassKind::Point => rep_point(geom).map(|(x, y)| ClassFeat {
                        poly: None,
                        point: Some(Point::new(x, y)),
                        value,
                        sums,
                    }),
                }
            })
            .collect();

        ctx.progress.info(&format!(
            "{} zone(s) x {} class feature(s) ({})",
            zones.len(),
            class_feats.len(),
            class_kind.as_str()
        ));

        // Accumulate per zone -> per class value.
        #[derive(Clone)]
        struct Accum {
            measure: f64,
            geom: MultiPolygon,
            sums: Vec<f64>,
        }
        let new_accum = || Accum {
            measure: 0.0,
            geom: MultiPolygon::new(vec![]),
            sums: vec![0.0; sum_fields.len()],
        };
        let measure_field = class_kind.measure_field();

        let mut out_layer = Layer::new(zones.name.clone());
        out_layer.crs = zones.crs.clone();
        out_layer.add_field(FieldDef::new("zone_id", FieldType::Text));
        if let Some(cf) = &class_field {
            out_layer.add_field(FieldDef::new(cf, FieldType::Text));
        }
        out_layer.add_field(FieldDef::new(measure_field, FieldType::Float));
        out_layer.add_field(FieldDef::new("percentage", FieldType::Float));
        for sf in &sum_fields {
            out_layer.add_field(FieldDef::new(sf, FieldType::Float));
        }
        out_layer.geom_type = Some(match class_kind {
            ClassKind::Polygon => GeometryType::MultiPolygon,
            ClassKind::Point => GeometryType::Polygon,
        });

        let mut rows = 0usize;
        for (zi, zone) in zones.features.iter().enumerate() {
            let Some(zone_mp) = zone.geometry.as_ref().and_then(to_multipolygon) else {
                continue;
            };
            let zone_id = match &zone_field {
                Some(f) => zone
                    .get(&zones.schema, f)
                    .map(field_value_string)
                    .unwrap_or_else(|_| zi.to_string()),
                None => zi.to_string(),
            };

            let mut by_class: BTreeMap<String, Accum> = BTreeMap::new();
            for cf in &class_feats {
                let acc = match class_kind {
                    ClassKind::Polygon => {
                        let cpoly = cf.poly.as_ref().unwrap();
                        let inter = zone_mp.intersection(cpoly);
                        let area = inter.unsigned_area();
                        if area <= 0.0 {
                            continue;
                        }
                        let class_area = cpoly.unsigned_area();
                        let frac = if class_area > 0.0 {
                            area / class_area
                        } else {
                            0.0
                        };
                        let entry = by_class.entry(cf.value.clone()).or_insert_with(new_accum);
                        entry.measure += area;
                        entry.geom = entry.geom.union(&inter);
                        for (s, v) in entry.sums.iter_mut().zip(&cf.sums) {
                            *s += v * frac;
                        }
                        continue;
                    }
                    ClassKind::Point => {
                        let p = cf.point.as_ref().unwrap();
                        if !zone_mp.contains(p) {
                            continue;
                        }
                        by_class.entry(cf.value.clone()).or_insert_with(new_accum)
                    }
                };
                // Point branch accumulation.
                acc.measure += 1.0;
                for (s, v) in acc.sums.iter_mut().zip(&cf.sums) {
                    *s += v;
                }
            }

            let total: f64 = by_class.values().map(|a| a.measure).sum();
            if total <= 0.0 {
                continue;
            }
            for (value, acc) in &by_class {
                let geom = match class_kind {
                    ClassKind::Polygon => multipolygon_to_geometry(&acc.geom),
                    ClassKind::Point => zone.geometry.clone().unwrap(),
                };
                let mut fields: Vec<(&str, FieldValue)> =
                    vec![("zone_id", FieldValue::Text(zone_id.clone()))];
                if let Some(cf) = &class_field {
                    fields.push((cf.as_str(), FieldValue::Text(value.clone())));
                }
                fields.push((measure_field, FieldValue::Float(acc.measure)));
                fields.push(("percentage", FieldValue::Float(acc.measure / total * 100.0)));
                for (sf, v) in sum_fields.iter().zip(&acc.sums) {
                    fields.push((sf.as_str(), FieldValue::Float(*v)));
                }
                out_layer
                    .add_feature(Some(geom), &fields)
                    .map_err(|e| ToolError::Execution(format!("failed writing row: {e}")))?;
                rows += 1;
            }
        }

        ctx.progress
            .info(&format!("wrote {rows} zone x class row(s)"));

        let feature_count = out_layer.len();
        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("zone_count".to_string(), json!(zones.len()));
        outputs.insert("class_kind".to_string(), json!(class_kind.as_str()));
        outputs.insert("row_count".to_string(), json!(rows));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        Ok(ToolRunResult { outputs })
    }
}

// ── Class kind ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum ClassKind {
    Polygon,
    Point,
}

impl ClassKind {
    fn detect(layer: &Layer) -> Option<Self> {
        layer
            .features
            .iter()
            .find_map(|f| match f.geometry.as_ref()? {
                Geometry::Polygon { .. } | Geometry::MultiPolygon(_) => Some(ClassKind::Polygon),
                Geometry::Point(_) | Geometry::MultiPoint(_) => Some(ClassKind::Point),
                _ => None,
            })
    }
    fn as_str(self) -> &'static str {
        match self {
            Self::Polygon => "polygon",
            Self::Point => "point",
        }
    }
    fn measure_field(self) -> &'static str {
        match self {
            Self::Polygon => "area",
            Self::Point => "count",
        }
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

fn rep_point(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) => cs.first().map(|c| (c.x, c.y)),
        _ => None,
    }
}

fn field_value_string(v: &FieldValue) -> String {
    match v {
        FieldValue::Text(s) | FieldValue::Date(s) | FieldValue::DateTime(s) => s.clone(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(x) => x.to_string(),
        FieldValue::Boolean(b) => b.to_string(),
        FieldValue::Null | FieldValue::Blob(_) => String::new(),
    }
}

fn to_multipolygon(geom: &Geometry) -> Option<MultiPolygon> {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => Some(MultiPolygon(vec![rings_to_polygon(exterior, interiors)])),
        Geometry::MultiPolygon(parts) => Some(MultiPolygon(
            parts.iter().map(|(e, i)| rings_to_polygon(e, i)).collect(),
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
    if coords.len() >= 2 && coords.first() == coords.last() {
        coords.pop();
    }
    Ring::new(coords)
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, Layer};

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

    fn store(layer: Layer) -> String {
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = TabulateIntersectionTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn fget(layer: &Layer, idx: usize, name: &str) -> FieldValue {
        layer.features[idx]
            .get(&layer.schema, name)
            .unwrap()
            .clone()
    }
    fn ffloat(layer: &Layer, idx: usize, name: &str) -> f64 {
        FieldValue::as_f64(&fget(layer, idx, name)).unwrap()
    }

    #[test]
    fn polygon_area_composition_sums_to_100() {
        // One zone (0..10 x 0..10, area 100). Two class polygons: class A covers
        // the left 0..4 (area 40), class B the right 4..10 (area 60).
        let mut zones = Layer::new("zones");
        zones
            .add_feature(Some(rect(0.0, 0.0, 10.0, 10.0)), &[])
            .unwrap();
        let zin = store(zones);

        let mut cls = Layer::new("lc");
        cls.add_field(FieldDef::new("class", FieldType::Text));
        cls.add_feature(
            Some(rect(0.0, 0.0, 4.0, 10.0)),
            &[("class", FieldValue::Text("A".into()))],
        )
        .unwrap();
        cls.add_feature(
            Some(rect(4.0, 0.0, 10.0, 10.0)),
            &[("class", FieldValue::Text("B".into()))],
        )
        .unwrap();
        let cin = store(cls);

        let (out, layer) = run_tool(json!({
            "input": zin, "class_features": cin, "class_field": "class"
        }));
        assert_eq!(out.outputs["class_kind"], json!("polygon"));
        assert_eq!(out.outputs["row_count"], json!(2));
        // Rows are (zone A) and (zone B); areas 40 and 60, percentages 40 and 60.
        let mut pairs: Vec<(String, f64, f64)> = (0..2)
            .map(|i| {
                let c = match fget(&layer, i, "class") {
                    FieldValue::Text(s) => s,
                    _ => unreachable!(),
                };
                (
                    c,
                    ffloat(&layer, i, "area"),
                    ffloat(&layer, i, "percentage"),
                )
            })
            .collect();
        pairs.sort_by(|a, b| a.0.cmp(&b.0));
        assert!((pairs[0].1 - 40.0).abs() < 1e-6 && (pairs[0].2 - 40.0).abs() < 1e-6);
        assert!((pairs[1].1 - 60.0).abs() < 1e-6 && (pairs[1].2 - 60.0).abs() < 1e-6);
        let pct_sum: f64 = pairs.iter().map(|p| p.2).sum();
        assert!((pct_sum - 100.0).abs() < 1e-6);
    }

    #[test]
    fn area_weighted_sum_field_apportionment() {
        // A class polygon of area 100 with population 1000 straddles two zones,
        // 30% in zone 0 and 70% in zone 1. Population apportions 300 / 700.
        let mut zones = Layer::new("zones");
        zones
            .add_feature(Some(rect(0.0, 0.0, 3.0, 10.0)), &[])
            .unwrap(); // 30
        zones
            .add_feature(Some(rect(3.0, 0.0, 10.0, 10.0)), &[])
            .unwrap(); // 70
        let zin = store(zones);
        let mut cls = Layer::new("tracts");
        cls.add_field(FieldDef::new("pop", FieldType::Float));
        cls.add_feature(
            Some(rect(0.0, 0.0, 10.0, 10.0)),
            &[("pop", FieldValue::Float(1000.0))],
        )
        .unwrap();
        let cin = store(cls);
        let (_, layer) = run_tool(json!({
            "input": zin, "class_features": cin, "sum_fields": "pop"
        }));
        let mut pops: Vec<f64> = (0..layer.len()).map(|i| ffloat(&layer, i, "pop")).collect();
        pops.sort_by(f64::total_cmp);
        assert!((pops[0] - 300.0).abs() < 1e-6, "zone0 pop {:?}", pops);
        assert!((pops[1] - 700.0).abs() < 1e-6, "zone1 pop {:?}", pops);
    }

    #[test]
    fn point_counts_and_summarize_within() {
        // 3 points, 2 in zone 0 and 1 in zone 1; sum a value field.
        let mut zones = Layer::new("zones");
        zones
            .add_feature(Some(rect(0.0, 0.0, 5.0, 5.0)), &[])
            .unwrap();
        zones
            .add_feature(Some(rect(5.0, 0.0, 10.0, 5.0)), &[])
            .unwrap();
        let zin = store(zones);
        let mut pts = Layer::new("pts");
        pts.add_field(FieldDef::new("val", FieldType::Float));
        for (x, y, v) in [(1.0, 1.0, 10.0), (2.0, 2.0, 20.0), (7.0, 1.0, 100.0)] {
            pts.add_feature(
                Some(Geometry::point(x, y)),
                &[("val", FieldValue::Float(v))],
            )
            .unwrap();
        }
        let cin = store(pts);
        let (out, layer) = run_tool(json!({
            "input": zin, "class_features": cin, "sum_fields": "val"
        }));
        assert_eq!(out.outputs["class_kind"], json!("point"));
        // zone 0: count 2, val 30; zone 1: count 1, val 100.
        let rows: Vec<(String, f64, f64)> = (0..layer.len())
            .map(|i| {
                let z = match fget(&layer, i, "zone_id") {
                    FieldValue::Text(s) => s,
                    _ => unreachable!(),
                };
                (z, ffloat(&layer, i, "count"), ffloat(&layer, i, "val"))
            })
            .collect();
        let z0 = rows.iter().find(|r| r.0 == "0").unwrap();
        let z1 = rows.iter().find(|r| r.0 == "1").unwrap();
        assert_eq!(z0.1, 2.0);
        assert!((z0.2 - 30.0).abs() < 1e-9);
        assert_eq!(z1.1, 1.0);
        assert!((z1.2 - 100.0).abs() < 1e-9);
    }

    #[test]
    fn zones_with_no_class_produce_no_rows() {
        let mut zones = Layer::new("zones");
        zones
            .add_feature(Some(rect(0.0, 0.0, 5.0, 5.0)), &[])
            .unwrap();
        zones
            .add_feature(Some(rect(100.0, 100.0, 105.0, 105.0)), &[])
            .unwrap(); // empty zone
        let zin = store(zones);
        let mut cls = Layer::new("lc");
        cls.add_feature(Some(rect(0.0, 0.0, 5.0, 5.0)), &[])
            .unwrap();
        let cin = store(cls);
        let (out, _) = run_tool(json!({ "input": zin, "class_features": cin }));
        assert_eq!(
            out.outputs["row_count"],
            json!(1),
            "only the overlapping zone yields a row"
        );
    }

    #[test]
    fn rejects_missing_inputs() {
        let tool = TabulateIntersectionTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "input": "z.geojson" })).is_err(),
            "missing class_features"
        );
        assert!(bad(json!({ "input": "z.geojson", "class_features": "c.geojson" })).is_ok());
    }
}
