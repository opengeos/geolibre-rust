//! GeoLibre tool: buffer input features and summarize a second layer within
//! each buffer, in one pass.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Summarize Nearby* (Analysis): draw
//! one or more straight-line distance buffers around every input feature and,
//! for each buffer, summarize the *summary_features* layer that falls inside it
//! — a count of features, the intersected area (polygon summary features), and
//! a sum/mean of one or more numeric fields. It combines what today needs three
//! GeoLibre tools chained by hand: `multiple_ring_buffer` builds the rings but
//! summarizes nothing, `tabulate_intersection` summarizes but wants you to
//! supply the zones, and `neighborhood_summary_statistics` is single-layer.
//!
//! Buffers are **cumulative disks**: distance `d` is the full buffer from the
//! feature out to `d` (`geo`'s `Buffer`, i_overlay backend — the same engine as
//! `BooleanOps`, no GDAL/GEOS). So "summarize nearby within `d`" means every
//! summary feature within `d` of the input feature (including any that overlaps
//! the feature itself). Listing several distances emits one output buffer per
//! input feature × distance, and because the disks nest, the summarized
//! quantities are monotonic in `d`.
//!
//! Summarization mirrors `tabulate_intersection`, dispatched on the summary
//! layer's geometry:
//!
//! - **point** summary features — a feature is counted when the buffer contains
//!   it; `count` is that tally and each `sum_fields` entry is the plain sum of
//!   the contained points' values.
//! - **polygon** summary features — a feature contributes its area intersected
//!   with the buffer; `count` is the number of polygons that overlap the buffer,
//!   `area_within` is the total intersected area, and each `sum_fields` entry is
//!   the field value apportioned by the intersected fraction
//!   (`value × intersected_area / feature_area`), exactly as Summarize Within
//!   apportions.
//!
//! Every output feature carries `input_id` (the value of `id_field`, or the
//! input feature index), the buffer `distance`, `count`, `area_within` (0 for
//! point summaries), and for each summarized field a `sum_<field>` and
//! `mean_<field>` (mean = sum / count, 0 when nothing is nearby). The geometry
//! is the buffer polygon, so the enriched result is directly mappable.
//!
//! Scope for v1: line summary features (intersected length) are not supported —
//! `geo` has no arbitrary line-in-polygon clip, the same limitation
//! `tabulate_intersection` carries; use polygon or point summary layers.
//! Distances are in the layer's CRS units (no on-the-fly unit conversion), input
//! attributes beyond `id_field` are not copied onto the buffers, and the summary
//! statistics are fixed to count/sum/mean.

use std::collections::BTreeMap;

use geo::{
    Area, BooleanOps, Buffer, Contains, Coord as GeoCoord, Geometry as GeoGeometry, LineString,
    MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct SummarizeNearbyTool;

impl Tool for SummarizeNearbyTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "summarize_nearby",
            display_name: "Summarize Nearby",
            summary: "Buffer input features at one or more straight-line distances and summarize a second layer within each buffer: count of features, intersected area (polygons), and area-weighted sum/mean of numeric fields.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer whose features are buffered (points, lines, or polygons), format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "summary_features",
                    description: "Vector layer to summarize within each buffer (polygons summarized by intersected area, or points by count).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distances",
                    description: "Comma-separated list of buffer distances in CRS units, e.g. \"1000,2000,5000\". Sorted ascending; duplicates and non-positive values are dropped.",
                    required: true,
                },
                ToolParamSpec {
                    name: "sum_fields",
                    description: "Optional comma-separated numeric fields in the summary layer to aggregate. Each yields a sum_<field> and mean_<field> column (area-weighted for polygons, plain for points).",
                    required: false,
                },
                ToolParamSpec {
                    name: "id_field",
                    description: "Optional field in the input layer identifying each source feature in the output. Defaults to the input feature index.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["input", "summary_features"] {
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
        parse_distances(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let summary_path = require_str(args, "summary_features")?;
        let output = parse_optional_str(args, "output")?;
        let distances = parse_distances(args)?;
        let id_field = parse_optional_str(args, "id_field")?.map(str::to_string);
        let sum_fields: Vec<String> = parse_optional_str(args, "sum_fields")?
            .map(|s| {
                s.split(',')
                    .map(str::trim)
                    .filter(|t| !t.is_empty())
                    .map(str::to_string)
                    .collect()
            })
            .unwrap_or_default();

        let in_layer = load_input_layer(input)?;
        let summary = load_input_layer(summary_path)?;
        let kind = SummaryKind::detect(&summary).ok_or_else(|| {
            ToolError::Validation(
                "summary_features has no polygon or point features (line summaries are not supported)"
                    .to_string(),
            )
        })?;

        // Pre-extract summary features once: (polygon | point, precomputed area,
        // sum-field values).
        let sum_schema = summary.schema.clone();
        let summary_feats: Vec<SummaryFeat> = summary
            .features
            .iter()
            .filter_map(|feat| {
                let geom = feat.geometry.as_ref()?;
                let sums: Vec<f64> = sum_fields
                    .iter()
                    .map(|f| {
                        feat.get(&sum_schema, f)
                            .ok()
                            .and_then(FieldValue::as_f64)
                            .unwrap_or(0.0)
                    })
                    .collect();
                match kind {
                    SummaryKind::Polygon => to_multipolygon(geom).map(|poly| {
                        let area = poly.unsigned_area();
                        SummaryFeat {
                            poly: Some(poly),
                            point: None,
                            area,
                            sums,
                        }
                    }),
                    SummaryKind::Point => rep_point(geom).map(|(x, y)| SummaryFeat {
                        poly: None,
                        point: Some(Point::new(x, y)),
                        area: 0.0,
                        sums,
                    }),
                }
            })
            .collect();

        ctx.progress.info(&format!(
            "{} input feature(s) x {} distance(s), summarizing {} {} feature(s)",
            in_layer.len(),
            distances.len(),
            summary_feats.len(),
            kind.as_str()
        ));

        // Build output schema.
        let mut out_layer = Layer::new(in_layer.name.clone());
        out_layer.crs = in_layer.crs.clone();
        out_layer.add_field(FieldDef::new("input_id", FieldType::Text));
        out_layer.add_field(FieldDef::new("distance", FieldType::Float));
        out_layer.add_field(FieldDef::new("count", FieldType::Integer));
        out_layer.add_field(FieldDef::new("area_within", FieldType::Float));
        for f in &sum_fields {
            out_layer.add_field(FieldDef::new(format!("sum_{f}"), FieldType::Float));
            out_layer.add_field(FieldDef::new(format!("mean_{f}"), FieldType::Float));
        }
        out_layer.geom_type = Some(GeometryType::MultiPolygon);

        let mut rows = 0usize;
        for (fi, feat) in in_layer.features.iter().enumerate() {
            let Some(src) = feat.geometry.as_ref().and_then(to_geo_geometry) else {
                continue; // geometry we cannot buffer (e.g. nested collection)
            };
            let input_id = match &id_field {
                Some(f) => feat
                    .get(&in_layer.schema, f)
                    .map(field_value_string)
                    .unwrap_or_else(|_| fi.to_string()),
                None => fi.to_string(),
            };
            for &d in &distances {
                let buffer = src.buffer(d);
                if buffer.0.is_empty() {
                    continue;
                }
                let stats = summarize(&buffer, &summary_feats, kind, sum_fields.len());

                let mut fields: Vec<(&str, FieldValue)> = vec![
                    ("input_id", FieldValue::Text(input_id.clone())),
                    ("distance", FieldValue::Float(d)),
                    ("count", FieldValue::Integer(stats.count as i64)),
                    ("area_within", FieldValue::Float(stats.area_within)),
                ];
                let labels: Vec<(String, String)> = sum_fields
                    .iter()
                    .map(|f| (format!("sum_{f}"), format!("mean_{f}")))
                    .collect();
                for (i, (sum_l, mean_l)) in labels.iter().enumerate() {
                    let sum = stats.sums[i];
                    let mean = if stats.count > 0 {
                        sum / stats.count as f64
                    } else {
                        0.0
                    };
                    fields.push((sum_l.as_str(), FieldValue::Float(sum)));
                    fields.push((mean_l.as_str(), FieldValue::Float(mean)));
                }
                out_layer
                    .add_feature(Some(multipolygon_to_geometry(&buffer)), &fields)
                    .map_err(|e| ToolError::Execution(format!("failed writing buffer row: {e}")))?;
                rows += 1;
            }
        }

        ctx.progress
            .info(&format!("wrote {rows} buffer summary row(s)"));

        let feature_count = out_layer.len();
        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_count".to_string(), json!(in_layer.len()));
        outputs.insert("summary_kind".to_string(), json!(kind.as_str()));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("row_count".to_string(), json!(rows));
        outputs.insert("distances".to_string(), json!(distances));
        Ok(ToolRunResult { outputs })
    }
}

// ── Summarization ──────────────────────────────────────────────────────────────

struct SummaryFeat {
    poly: Option<MultiPolygon>,
    point: Option<Point>,
    area: f64,
    sums: Vec<f64>,
}

struct Stats {
    count: usize,
    area_within: f64,
    sums: Vec<f64>,
}

/// Aggregates the summary features that fall within `buffer`.
fn summarize(
    buffer: &MultiPolygon,
    feats: &[SummaryFeat],
    kind: SummaryKind,
    n_sum: usize,
) -> Stats {
    let mut stats = Stats {
        count: 0,
        area_within: 0.0,
        sums: vec![0.0; n_sum],
    };
    for f in feats {
        match kind {
            SummaryKind::Point => {
                let p = f.point.as_ref().unwrap();
                if !buffer.contains(p) {
                    continue;
                }
                stats.count += 1;
                for (s, v) in stats.sums.iter_mut().zip(&f.sums) {
                    *s += v;
                }
            }
            SummaryKind::Polygon => {
                let poly = f.poly.as_ref().unwrap();
                let inter = buffer.intersection(poly);
                let area = inter.unsigned_area();
                if area <= 0.0 {
                    continue;
                }
                stats.count += 1;
                stats.area_within += area;
                let frac = if f.area > 0.0 { area / f.area } else { 0.0 };
                for (s, v) in stats.sums.iter_mut().zip(&f.sums) {
                    *s += v * frac;
                }
            }
        }
    }
    stats
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum SummaryKind {
    Polygon,
    Point,
}

impl SummaryKind {
    fn detect(layer: &Layer) -> Option<Self> {
        layer
            .features
            .iter()
            .find_map(|f| match f.geometry.as_ref()? {
                Geometry::Polygon { .. } | Geometry::MultiPolygon(_) => Some(SummaryKind::Polygon),
                Geometry::Point(_) | Geometry::MultiPoint(_) => Some(SummaryKind::Point),
                _ => None,
            })
    }
    fn as_str(self) -> &'static str {
        match self {
            Self::Polygon => "polygon",
            Self::Point => "point",
        }
    }
}

// ── Parameters ──────────────────────────────────────────────────────────────────

fn parse_distances(args: &ToolArgs) -> Result<Vec<f64>, ToolError> {
    let raw = parse_optional_str(args, "distances")?.ok_or_else(|| {
        ToolError::Validation("missing required parameter 'distances'".to_string())
    })?;
    let mut distances: Vec<f64> = Vec::new();
    for tok in raw.split(',') {
        let t = tok.trim();
        if t.is_empty() {
            continue;
        }
        let v: f64 = t
            .parse()
            .map_err(|_| ToolError::Validation(format!("distance '{t}' is not a number")))?;
        if v > 0.0 && v.is_finite() {
            distances.push(v);
        }
    }
    distances.sort_by(f64::total_cmp);
    distances.dedup();
    if distances.is_empty() {
        return Err(ToolError::Validation(
            "parameter 'distances' must list at least one positive distance".to_string(),
        ));
    }
    Ok(distances)
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
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

fn rep_point(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) => cs.first().map(|c| (c.x, c.y)),
        _ => None,
    }
}

// ── geo <-> wbvector geometry conversion ────────────────────────────────────────

/// Converts a `wbvector` geometry to a `geo` geometry for buffering. Nested
/// `GeometryCollection`s are skipped (returns `None`).
fn to_geo_geometry(g: &Geometry) -> Option<GeoGeometry> {
    let geom = match g {
        Geometry::Point(c) => GeoGeometry::Point(Point::new(c.x, c.y)),
        Geometry::MultiPoint(cs) => GeoGeometry::MultiPoint(MultiPoint(
            cs.iter().map(|c| Point::new(c.x, c.y)).collect(),
        )),
        Geometry::LineString(cs) => GeoGeometry::LineString(coords_to_linestring(cs)),
        Geometry::MultiLineString(ls) => GeoGeometry::MultiLineString(MultiLineString(
            ls.iter().map(|c| coords_to_linestring(c)).collect(),
        )),
        Geometry::Polygon {
            exterior,
            interiors,
        } => GeoGeometry::Polygon(rings_to_polygon(exterior, interiors)),
        Geometry::MultiPolygon(parts) => GeoGeometry::MultiPolygon(MultiPolygon(
            parts.iter().map(|(e, i)| rings_to_polygon(e, i)).collect(),
        )),
        Geometry::GeometryCollection(_) => return None,
    };
    Some(geom)
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

fn coords_to_linestring(coords: &[Coord]) -> LineString {
    LineString::new(coords.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect())
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
    // Drop the closing duplicate vertex `geo` keeps; `Ring` stores it implicitly.
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
        let out = SummarizeNearbyTool.run(&args, &ctx()).unwrap();
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
    fn fint(layer: &Layer, idx: usize, name: &str) -> i64 {
        match fget(layer, idx, name) {
            FieldValue::Integer(i) => i,
            other => panic!("expected integer, got {other:?}"),
        }
    }

    /// One input point buffered at three distances over a field of summary
    /// points; counts must equal the number of points within each radius and be
    /// monotonically non-decreasing in distance.
    #[test]
    fn point_counts_are_monotonic_in_distance() {
        // Input: a single origin point.
        let mut inp = Layer::new("origin");
        inp.add_feature(Some(Geometry::point(0.0, 0.0)), &[])
            .unwrap();
        let iin = store(inp);

        // Summary points at radii 5, 15, 25, 100 from the origin.
        let mut pts = Layer::new("pts");
        for (x, y) in [(5.0, 0.0), (0.0, 15.0), (25.0, 0.0), (100.0, 0.0)] {
            pts.add_feature(Some(Geometry::point(x, y)), &[]).unwrap();
        }
        let sin = store(pts);

        let (out, layer) = run_tool(json!({
            "input": iin, "summary_features": sin, "distances": "10,20,30"
        }));
        assert_eq!(out.outputs["summary_kind"], json!("point"));
        assert_eq!(out.outputs["row_count"], json!(3));
        // Rows in ascending-distance order (single input feature).
        let counts: Vec<i64> = (0..3).map(|i| fint(&layer, i, "count")).collect();
        // d=10 -> {5}; d=20 -> {5,15}; d=30 -> {5,15,25}.
        assert_eq!(counts, vec![1, 2, 3]);
        // Monotonic non-decreasing.
        assert!(counts.windows(2).all(|w| w[1] >= w[0]));
    }

    /// Point summary with a sum field: sum and mean over contained points.
    #[test]
    fn point_sum_and_mean_fields() {
        let mut inp = Layer::new("origin");
        inp.add_feature(Some(Geometry::point(0.0, 0.0)), &[])
            .unwrap();
        let iin = store(inp);

        let mut pts = Layer::new("pts");
        pts.add_field(FieldDef::new("pop", FieldType::Float));
        // Two points within radius 10 (pop 100, 300), one outside (pop 999).
        pts.add_feature(Some(Geometry::point(3.0, 0.0)), &[("pop", 100.0.into())])
            .unwrap();
        pts.add_feature(Some(Geometry::point(0.0, 4.0)), &[("pop", 300.0.into())])
            .unwrap();
        pts.add_feature(Some(Geometry::point(50.0, 0.0)), &[("pop", 999.0.into())])
            .unwrap();
        let sin = store(pts);

        let (_, layer) = run_tool(json!({
            "input": iin, "summary_features": sin, "distances": "10", "sum_fields": "pop"
        }));
        assert_eq!(fint(&layer, 0, "count"), 2);
        assert!((ffloat(&layer, 0, "sum_pop") - 400.0).abs() < 1e-9);
        assert!((ffloat(&layer, 0, "mean_pop") - 200.0).abs() < 1e-9);
    }

    /// Polygon summary features: area_within and area-weighted sum. A summary
    /// polygon of area 100 (pop 1000) sits with its left half inside a buffer
    /// large enough to capture it partially. We use a big square buffer via a
    /// polygon input so the intersection fraction is exact and checkable.
    #[test]
    fn polygon_area_within_and_weighted_sum() {
        // Input polygon: the left half-plane strip 0..5 x -100..100 (buffered by
        // a tiny distance so it stays essentially itself).
        let mut inp = Layer::new("zone");
        inp.add_feature(Some(rect(-100.0, -100.0, 5.0, 100.0)), &[])
            .unwrap();
        let iin = store(inp);

        // Summary polygon 0..10 x 0..10 (area 100), pop 1000. Its left half
        // (x 0..5) lies within the input strip -> 50% inside.
        let mut poly = Layer::new("tracts");
        poly.add_field(FieldDef::new("pop", FieldType::Float));
        poly.add_feature(Some(rect(0.0, 0.0, 10.0, 10.0)), &[("pop", 1000.0.into())])
            .unwrap();
        let sin = store(poly);

        // Buffer distance small (0.001) so the input polygon barely grows.
        let (out, layer) = run_tool(json!({
            "input": iin, "summary_features": sin, "distances": "0.001", "sum_fields": "pop"
        }));
        assert_eq!(out.outputs["summary_kind"], json!("polygon"));
        assert_eq!(fint(&layer, 0, "count"), 1);
        // ~half of the 100-area polygon is inside -> ~50 area, ~500 pop.
        assert!(
            (ffloat(&layer, 0, "area_within") - 50.0).abs() < 0.5,
            "area_within {}",
            ffloat(&layer, 0, "area_within")
        );
        assert!(
            (ffloat(&layer, 0, "sum_pop") - 500.0).abs() < 5.0,
            "sum_pop {}",
            ffloat(&layer, 0, "sum_pop")
        );
    }

    /// A buffer with nothing nearby yields count 0 and zeroed statistics (row
    /// still emitted — the buffer exists even if empty).
    #[test]
    fn empty_buffer_yields_zero_counts() {
        let mut inp = Layer::new("origin");
        inp.add_feature(Some(Geometry::point(0.0, 0.0)), &[])
            .unwrap();
        let iin = store(inp);
        let mut pts = Layer::new("pts");
        pts.add_field(FieldDef::new("v", FieldType::Float));
        pts.add_feature(Some(Geometry::point(500.0, 500.0)), &[("v", 7.0.into())])
            .unwrap();
        let sin = store(pts);
        let (out, layer) = run_tool(json!({
            "input": iin, "summary_features": sin, "distances": "10", "sum_fields": "v"
        }));
        assert_eq!(out.outputs["row_count"], json!(1));
        assert_eq!(fint(&layer, 0, "count"), 0);
        assert_eq!(ffloat(&layer, 0, "sum_v"), 0.0);
        assert_eq!(ffloat(&layer, 0, "mean_v"), 0.0);
    }

    /// Custom id_field is echoed onto every buffer row.
    #[test]
    fn id_field_is_carried_through() {
        let mut inp = Layer::new("stores");
        inp.add_field(FieldDef::new("name", FieldType::Text));
        inp.add_feature(Some(Geometry::point(0.0, 0.0)), &[("name", "A".into())])
            .unwrap();
        let iin = store(inp);
        let mut pts = Layer::new("pts");
        pts.add_feature(Some(Geometry::point(1.0, 0.0)), &[])
            .unwrap();
        let sin = store(pts);
        let (_, layer) = run_tool(json!({
            "input": iin, "summary_features": sin, "distances": "5", "id_field": "name"
        }));
        assert_eq!(fget(&layer, 0, "input_id"), FieldValue::Text("A".into()));
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = SummarizeNearbyTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "input": "a.geojson" })).is_err(),
            "missing summary_features"
        );
        assert!(
            bad(json!({ "input": "a.geojson", "summary_features": "b.geojson" })).is_err(),
            "missing distances"
        );
        assert!(bad(
            json!({ "input": "a.geojson", "summary_features": "b.geojson", "distances": "" })
        )
        .is_err());
        assert!(bad(
            json!({ "input": "a.geojson", "summary_features": "b.geojson", "distances": "-5,0" })
        )
        .is_err());
        assert!(bad(
            json!({ "input": "a.geojson", "summary_features": "b.geojson", "distances": "1,x" })
        )
        .is_err());
        assert!(bad(
            json!({ "input": "a.geojson", "summary_features": "b.geojson", "distances": "1000,2000" })
        )
        .is_ok());
    }
}
