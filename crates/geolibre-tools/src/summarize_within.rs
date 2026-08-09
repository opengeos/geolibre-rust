//! GeoLibre tool: summarize a layer inside an existing polygon layer, splitting
//! partially-overlapping features by area or length.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Summarize Within* (Analysis). The
//! shipped `summarize_nearby` only summarizes inside buffers it generates
//! itself, so it cannot summarize into polygons the user already has;
//! `tabulate_intersection` reports overlap areas but computes no attribute
//! statistics; and the bundled `spatial_join` attaches attributes without any
//! proportional splitting, assigning a straddling feature wholly to one zone.
//!
//! Each summarized feature is weighted by the fraction of itself that falls
//! inside a polygon — intersected area for polygons, intersected length for
//! lines, and 1.0 for a contained point. Statistics are accumulated with those
//! weights, so a parcel half inside a tract contributes half its population.
//! Because the weight is a *fraction of the source feature*, a coverage that
//! tiles the summarized layer reproduces the original totals exactly.
//!
//! `group_field` additionally emits one row per distinct value of a field on the
//! summarized layer, which is the standard "population by category per zone"
//! shape.

use std::collections::BTreeMap;

use geo::{
    Area, BooleanOps, BoundingRect, Coord as GeoCoord, Euclidean, Length, LineString,
    MultiLineString, MultiPolygon, Polygon,
};
use rstar::{RTree, RTreeObject, AABB};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, Layer, Ring};

use crate::vector_common::{
    geometry_contains_point, load_input_layer, parse_optional_str, to_multilinestring,
    write_or_store_layer,
};

/// What kind of geometry the summarized layer holds; decides the weighting rule.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum SummaryKind {
    Polygon,
    Line,
    Point,
}

impl SummaryKind {
    fn as_str(self) -> &'static str {
        match self {
            SummaryKind::Polygon => "polygon",
            SummaryKind::Line => "line",
            SummaryKind::Point => "point",
        }
    }
}

/// One statistic requested for one field.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Stat {
    Sum,
    Mean,
    Min,
    Max,
    StdDev,
    Count,
}

impl Stat {
    fn parse(s: &str) -> Option<Stat> {
        match s.trim().to_ascii_lowercase().as_str() {
            "sum" => Some(Stat::Sum),
            "mean" | "avg" => Some(Stat::Mean),
            "min" => Some(Stat::Min),
            "max" => Some(Stat::Max),
            "stddev" | "std" => Some(Stat::StdDev),
            "count" => Some(Stat::Count),
            _ => None,
        }
    }
    fn label(self) -> &'static str {
        match self {
            Stat::Sum => "sum",
            Stat::Mean => "mean",
            Stat::Min => "min",
            Stat::Max => "max",
            Stat::StdDev => "std",
            Stat::Count => "count",
        }
    }
}

/// Weighted accumulator. Mean/stddev use a weighted Welford update so a long
/// run of small weights does not lose precision the way sum-of-squares does.
#[derive(Default, Clone)]
struct Acc {
    wsum: f64,
    mean: f64,
    m2: f64,
    total: f64,
    min: f64,
    max: f64,
    n: usize,
}

impl Acc {
    fn new() -> Acc {
        Acc {
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
            ..Default::default()
        }
    }
    fn push(&mut self, value: f64, weight: f64) {
        if weight <= 0.0 || !value.is_finite() {
            return;
        }
        self.n += 1;
        self.total += value * weight;
        self.min = self.min.min(value);
        self.max = self.max.max(value);
        let new_w = self.wsum + weight;
        let delta = value - self.mean;
        self.mean += delta * weight / new_w;
        self.m2 += weight * delta * (value - self.mean);
        self.wsum = new_w;
    }
    fn value(&self, stat: Stat) -> f64 {
        match stat {
            Stat::Sum => self.total,
            Stat::Mean => {
                if self.wsum > 0.0 {
                    self.mean
                } else {
                    0.0
                }
            }
            Stat::Min => {
                if self.n > 0 {
                    self.min
                } else {
                    0.0
                }
            }
            Stat::Max => {
                if self.n > 0 {
                    self.max
                } else {
                    0.0
                }
            }
            Stat::StdDev => {
                if self.wsum > 0.0 {
                    (self.m2 / self.wsum).max(0.0).sqrt()
                } else {
                    0.0
                }
            }
            Stat::Count => self.n as f64,
        }
    }
}

pub struct SummarizeWithinTool;

impl Tool for SummarizeWithinTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "summarize_within",
            display_name: "Summarize Within",
            summary: "Summarize a layer's attributes inside an existing polygon layer, weighting partially-overlapping features by intersected area or length, like ArcGIS Summarize Within.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "polygons",
                    description: "Polygon layer to summarize into. Its features and attributes are preserved.",
                    required: true,
                },
                ToolParamSpec {
                    name: "input",
                    description: "Layer being summarized (point, line, or polygon).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output polygon path. If omitted, the result is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "fields",
                    description: "Comma-separated 'field:statistic' pairs, e.g. 'pop:sum,income:mean'. Statistics: sum, mean, min, max, stddev, count.",
                    required: false,
                },
                ToolParamSpec {
                    name: "keep_all",
                    description: "Keep polygons with no intersecting features (default true). When false, empty polygons are dropped.",
                    required: false,
                },
                ToolParamSpec {
                    name: "shape_sum",
                    description: "Append the summed count / length / area of intersected input geometry (default true).",
                    required: false,
                },
                ToolParamSpec {
                    name: "group_field",
                    description: "Optional field on the summarized layer; emits one grouped summary row per distinct value.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "polygons")?;
        require_str(args, "input")?;
        parse_fields(args)?;
        parse_optional_bool(args, "keep_all")?;
        parse_optional_bool(args, "shape_sum")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let poly_path = require_str(args, "polygons")?;
        let input_path = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let field_specs = parse_fields(args)?;
        let keep_all = parse_optional_bool(args, "keep_all")?.unwrap_or(true);
        let shape_sum = parse_optional_bool(args, "shape_sum")?.unwrap_or(true);
        let group_field = parse_optional_str(args, "group_field")?.map(String::from);

        let zones = load_input_layer(poly_path)?;
        let src = load_input_layer(input_path)?;

        for (f, _) in field_specs.iter() {
            if src.schema.field_index(f).is_none() {
                return Err(ToolError::Validation(format!(
                    "field '{f}' not found on the summarized layer"
                )));
            }
        }
        if let Some(g) = &group_field {
            if src.schema.field_index(g).is_none() {
                return Err(ToolError::Validation(format!(
                    "group_field '{g}' not found on the summarized layer"
                )));
            }
        }

        let kind = detect_kind(&src).ok_or_else(|| {
            ToolError::Validation("input layer has no point, line, or polygon features".to_string())
        })?;

        // Pre-extract the summarized features once: geometry in `geo` form, its
        // own total measure (for the weight denominator), field values, group.
        let src_feats: Vec<SrcFeat> = src
            .features
            .iter()
            .filter_map(|f| {
                let geom = f.geometry.as_ref()?;
                let values: Vec<f64> = field_specs
                    .iter()
                    .map(|(name, _)| {
                        f.get(&src.schema, name)
                            .ok()
                            .and_then(FieldValue::as_f64)
                            .unwrap_or(f64::NAN)
                    })
                    .collect();
                let group = group_field.as_ref().map(|g| {
                    f.get(&src.schema, g)
                        .map(field_value_string)
                        .unwrap_or_default()
                });
                build_src_feat(geom, kind, values, group)
            })
            .collect();

        ctx.progress.info(&format!(
            "{} zone(s), {} {} feature(s), {} statistic(s)",
            zones.len(),
            src_feats.len(),
            kind.as_str(),
            field_specs.len()
        ));

        // Output schema: the zone attributes, then the summary columns.
        let mut out = Layer::new(zones.name.clone());
        out.crs = zones.crs.clone();
        out.geom_type = zones.geom_type;
        for fd in zones.schema.fields().iter() {
            out.add_field(fd.clone());
        }
        if group_field.is_some() {
            out.add_field(FieldDef::new("group", FieldType::Text));
        }
        out.add_field(FieldDef::new("count", FieldType::Integer));
        if shape_sum {
            out.add_field(FieldDef::new(shape_sum_label(kind), FieldType::Float));
        }
        let stat_labels: Vec<String> = field_specs
            .iter()
            .map(|(f, s)| format!("{}_{}", f, s.label()))
            .collect();
        for l in &stat_labels {
            out.add_field(FieldDef::new(l.clone(), FieldType::Float));
        }

        struct IndexedFeature {
            index: usize,
            envelope: AABB<[f64; 2]>,
        }
        impl RTreeObject for IndexedFeature {
            type Envelope = AABB<[f64; 2]>;
            fn envelope(&self) -> Self::Envelope {
                self.envelope
            }
        }
        let source_index = RTree::bulk_load(
            src_feats
                .iter()
                .enumerate()
                .filter_map(|(index, sf)| {
                    sf.bbox.map(|bbox| IndexedFeature {
                        index,
                        envelope: AABB::from_corners(
                            [bbox.min().x, bbox.min().y],
                            [bbox.max().x, bbox.max().y],
                        ),
                    })
                })
                .collect(),
        );

        let mut rows = 0usize;
        let mut matched_zones = 0usize;
        for (zi, zf) in zones.features.iter().enumerate() {
            let Some(zgeom) = zf.geometry.as_ref() else {
                continue;
            };
            let Some(zpoly) = to_multipolygon(zgeom) else {
                continue; // non-areal zone feature: nothing to summarize into
            };
            let zbox = zpoly.bounding_rect();

            // group key -> (accumulators, count, shape measure)
            let mut groups: BTreeMap<String, (Vec<Acc>, usize, f64)> = BTreeMap::new();

            let Some(zbox) = zbox else { continue };
            let envelope =
                AABB::from_corners([zbox.min().x, zbox.min().y], [zbox.max().x, zbox.max().y]);
            for indexed in source_index.locate_in_envelope_intersecting(&envelope) {
                let sf = &src_feats[indexed.index];
                let (weight, measure) = match overlap(&zpoly, zgeom, sf, kind) {
                    Some(v) => v,
                    None => continue,
                };
                if weight <= 0.0 {
                    continue;
                }
                let key = sf.group.clone().unwrap_or_default();
                let entry = groups
                    .entry(key)
                    .or_insert_with(|| (vec![Acc::new(); field_specs.len()], 0, 0.0));
                entry.1 += 1;
                entry.2 += measure;
                for (i, v) in sf.values.iter().enumerate() {
                    entry.0[i].push(*v, weight);
                }
            }

            if !groups.is_empty() {
                matched_zones += 1;
            }
            if groups.is_empty() {
                if !keep_all {
                    continue;
                }
                groups.insert(String::new(), (vec![Acc::new(); field_specs.len()], 0, 0.0));
            }

            for (key, (accs, count, measure)) in groups {
                let mut fields: Vec<(String, FieldValue)> = Vec::new();
                for (i, fd) in zones.schema.fields().iter().enumerate() {
                    fields.push((fd.name.clone(), zf.attributes[i].clone()));
                }
                if group_field.is_some() {
                    fields.push(("group".to_string(), FieldValue::Text(key)));
                }
                fields.push(("count".to_string(), FieldValue::Integer(count as i64)));
                if shape_sum {
                    fields.push((
                        shape_sum_label(kind).to_string(),
                        FieldValue::Float(measure),
                    ));
                }
                for (i, (_, stat)) in field_specs.iter().enumerate() {
                    fields.push((
                        stat_labels[i].clone(),
                        FieldValue::Float(accs[i].value(*stat)),
                    ));
                }
                let refs: Vec<(&str, FieldValue)> = fields
                    .iter()
                    .map(|(k, v)| (k.as_str(), v.clone()))
                    .collect();
                out.add_feature(Some(zgeom.clone()), &refs).map_err(|e| {
                    ToolError::Execution(format!("failed writing summary row: {e}"))
                })?;
                rows += 1;
            }
            ctx.progress
                .progress((zi as f64 + 1.0) / zones.len().max(1) as f64);
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(rows));
        outputs.insert("zone_count".to_string(), json!(zones.len()));
        outputs.insert("zones_with_matches".to_string(), json!(matched_zones));
        outputs.insert("summary_kind".to_string(), json!(kind.as_str()));
        Ok(ToolRunResult { outputs })
    }
}

/// A summarized feature, pre-converted and pre-measured.
struct SrcFeat {
    poly: Option<MultiPolygon>,
    lines: Option<MultiLineString>,
    point: Option<(f64, f64)>,
    /// Total measure of the whole feature (area or length); the weight denominator.
    total: f64,
    bbox: Option<geo::Rect<f64>>,
    values: Vec<f64>,
    group: Option<String>,
}

fn build_src_feat(
    geom: &Geometry,
    kind: SummaryKind,
    values: Vec<f64>,
    group: Option<String>,
) -> Option<SrcFeat> {
    match kind {
        SummaryKind::Polygon => {
            let mp = to_multipolygon(geom)?;
            let total = mp.unsigned_area();
            let bbox = mp.bounding_rect();
            Some(SrcFeat {
                poly: Some(mp),
                lines: None,
                point: None,
                total,
                bbox,
                values,
                group,
            })
        }
        SummaryKind::Line => {
            let ml = to_multilinestring(geom)?;
            let total = Euclidean.length(&ml);
            let bbox = ml.bounding_rect();
            Some(SrcFeat {
                poly: None,
                lines: Some(ml),
                point: None,
                total,
                bbox,
                values,
                group,
            })
        }
        SummaryKind::Point => {
            let (x, y) = rep_point(geom)?;
            Some(SrcFeat {
                poly: None,
                lines: None,
                point: Some((x, y)),
                total: 1.0,
                bbox: Some(geo::Rect::new(GeoCoord { x, y }, GeoCoord { x, y })),
                values,
                group,
            })
        }
    }
}

/// Returns `(weight, measure)` for a summarized feature against a zone:
/// weight is the fraction of the feature inside the zone, measure is the
/// absolute intersected area/length (or 1 for a contained point).
fn overlap(
    zone: &MultiPolygon,
    zone_geom: &Geometry,
    sf: &SrcFeat,
    kind: SummaryKind,
) -> Option<(f64, f64)> {
    match kind {
        SummaryKind::Polygon => {
            let mp = sf.poly.as_ref()?;
            if sf.total <= 0.0 {
                return None;
            }
            let inter = zone.intersection(mp);
            let a = inter.unsigned_area();
            if a <= 0.0 {
                return None;
            }
            Some(((a / sf.total).min(1.0), a))
        }
        SummaryKind::Line => {
            let ml = sf.lines.as_ref()?;
            if sf.total <= 0.0 {
                return None;
            }
            // `geo`'s BooleanOps clips a linestring against a polygon; the clipped
            // length divided by the original is the fraction inside.
            let clipped = zone.clip(ml, false);
            let l = Euclidean.length(&clipped);
            if l <= 0.0 {
                return None;
            }
            Some(((l / sf.total).min(1.0), l))
        }
        SummaryKind::Point => {
            let (x, y) = sf.point?;
            if geometry_contains_point(zone_geom, x, y) {
                Some((1.0, 1.0))
            } else {
                None
            }
        }
    }
}

fn shape_sum_label(kind: SummaryKind) -> &'static str {
    match kind {
        SummaryKind::Polygon => "area_within",
        SummaryKind::Line => "length_within",
        SummaryKind::Point => "points_within",
    }
}

fn detect_kind(layer: &Layer) -> Option<SummaryKind> {
    for f in layer.features.iter() {
        if let Some(g) = f.geometry.as_ref() {
            match g {
                Geometry::Polygon { .. } | Geometry::MultiPolygon(_) => {
                    return Some(SummaryKind::Polygon)
                }
                Geometry::LineString(_) | Geometry::MultiLineString(_) => {
                    return Some(SummaryKind::Line)
                }
                Geometry::Point(_) | Geometry::MultiPoint(_) => return Some(SummaryKind::Point),
                _ => continue,
            }
        }
    }
    None
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

// ── parameter parsing ────────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    match args.get(key).and_then(Value::as_str) {
        Some(s) if !s.trim().is_empty() => Ok(s),
        _ => Err(ToolError::Validation(format!(
            "missing required string parameter '{key}'"
        ))),
    }
}

fn parse_optional_bool(args: &ToolArgs, k: &str) -> Result<Option<bool>, ToolError> {
    match args.get(k) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{k}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{k}' must be a boolean"
        ))),
    }
}

/// Parses `fields` as `name:stat` pairs. A bare `name` defaults to `sum`.
fn parse_fields(args: &ToolArgs) -> Result<Vec<(String, Stat)>, ToolError> {
    let raw = match parse_optional_str(args, "fields")? {
        Some(s) => s,
        None => return Ok(Vec::new()),
    };
    let mut out = Vec::new();
    for tok in raw.split(',') {
        let tok = tok.trim();
        if tok.is_empty() {
            continue;
        }
        let (name, stat) = match tok.split_once(':') {
            Some((n, s)) => {
                let parsed = Stat::parse(s).ok_or_else(|| {
                    ToolError::Validation(format!(
                        "unknown statistic '{s}' (use sum, mean, min, max, stddev, count)"
                    ))
                })?;
                (n.trim().to_string(), parsed)
            }
            None => (tok.to_string(), Stat::Sum),
        };
        if name.is_empty() {
            return Err(ToolError::Validation(
                "'fields' entry is missing a field name".to_string(),
            ));
        }
        out.push((name, stat));
    }
    Ok(out)
}

// ── geo <-> wbvector conversion ──────────────────────────────────────────────

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

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Coord, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    #[test]
    fn line_overlap_uses_clipped_length() {
        let zone_geom = Geometry::polygon(
            vec![
                Coord::xy(0.0, 0.0),
                Coord::xy(10.0, 0.0),
                Coord::xy(10.0, 10.0),
                Coord::xy(0.0, 10.0),
            ],
            vec![],
        );
        let zone = to_multipolygon(&zone_geom).unwrap();
        let line = Geometry::LineString(vec![Coord::xy(-5.0, 5.0), Coord::xy(15.0, 5.0)]);
        let sf = build_src_feat(&line, SummaryKind::Line, vec![], None).unwrap();
        let (weight, length) = overlap(&zone, &zone_geom, &sf, SummaryKind::Line).unwrap();
        assert!((weight - 0.5).abs() < 1e-9);
        assert!((length - 10.0).abs() < 1e-9);
    }

    fn rect(x0: f64, y0: f64, w: f64, h: f64) -> Geometry {
        Geometry::polygon(
            vec![
                Coord::xy(x0, y0),
                Coord::xy(x0 + w, y0),
                Coord::xy(x0 + w, y0 + h),
                Coord::xy(x0, y0 + h),
            ],
            vec![],
        )
    }

    fn zones_layer(geoms: Vec<Geometry>) -> String {
        let mut l = Layer::new("zones")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("zid", FieldType::Integer));
        for (i, g) in geoms.into_iter().enumerate() {
            l.add_feature(Some(g), &[("zid", FieldValue::Integer(i as i64))])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    /// Polygons carrying a numeric `pop` field (and optional `cat` group).
    fn src_polys(items: Vec<(Geometry, f64, &str)>) -> String {
        let mut l = Layer::new("src")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("pop", FieldType::Float));
        l.add_field(FieldDef::new("cat", FieldType::Text));
        for (g, pop, cat) in items {
            l.add_feature(
                Some(g),
                &[
                    ("pop", FieldValue::Float(pop)),
                    ("cat", FieldValue::Text(cat.to_string())),
                ],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn src_points(items: Vec<(f64, f64, f64)>) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (x, y, v) in items {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(x, y))),
                &[("v", FieldValue::Float(v))],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SummarizeWithinTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn num(layer: &Layer, row: usize, field: &str) -> f64 {
        let i = layer.schema.field_index(field).unwrap();
        layer.features[row].attributes[i].as_f64().unwrap()
    }

    /// A source polygon exactly half inside a zone contributes half its value.
    #[test]
    fn partial_overlap_is_area_weighted() {
        // Zone [0,10]x[0,10]; source [5,15]x[0,10] with pop=100 -> half inside.
        let zones = zones_layer(vec![rect(0.0, 0.0, 10.0, 10.0)]);
        let src = src_polys(vec![(rect(5.0, 0.0, 10.0, 10.0), 100.0, "a")]);
        let (_, layer) = run(json!({
            "polygons": zones, "input": src, "fields": "pop:sum"
        }));
        assert!((num(&layer, 0, "pop_sum") - 50.0).abs() < 1e-9);
        assert!((num(&layer, 0, "area_within") - 50.0).abs() < 1e-9);
    }

    /// A coverage that tiles the source reproduces the original total exactly —
    /// the invariant that makes apportionment trustworthy.
    #[test]
    fn tiling_zones_conserve_the_total() {
        // Two zones split [0,10]x[0,10] down the middle; one source covers it all.
        let zones = zones_layer(vec![rect(0.0, 0.0, 5.0, 10.0), rect(5.0, 0.0, 5.0, 10.0)]);
        let src = src_polys(vec![(rect(0.0, 0.0, 10.0, 10.0), 80.0, "a")]);
        let (_, layer) = run(json!({
            "polygons": zones, "input": src, "fields": "pop:sum"
        }));
        let total = num(&layer, 0, "pop_sum") + num(&layer, 1, "pop_sum");
        assert!(
            (total - 80.0).abs() < 1e-9,
            "expected 80, got {total} (halves {} + {})",
            num(&layer, 0, "pop_sum"),
            num(&layer, 1, "pop_sum")
        );
    }

    /// Points are counted whole, never split.
    #[test]
    fn points_are_unweighted() {
        let zones = zones_layer(vec![rect(0.0, 0.0, 10.0, 10.0)]);
        let src = src_points(vec![(1.0, 1.0, 5.0), (2.0, 2.0, 7.0), (50.0, 50.0, 99.0)]);
        let (_, layer) = run(json!({
            "polygons": zones, "input": src, "fields": "v:sum"
        }));
        assert_eq!(
            num(&layer, 0, "count"),
            2.0,
            "outside point must be excluded"
        );
        assert!((num(&layer, 0, "v_sum") - 12.0).abs() < 1e-9);
        assert!((num(&layer, 0, "points_within") - 2.0).abs() < 1e-9);
    }

    /// keep_all=false drops zones with no intersecting features.
    #[test]
    fn keep_all_controls_empty_zones() {
        let zones = zones_layer(vec![rect(0.0, 0.0, 5.0, 5.0), rect(100.0, 100.0, 5.0, 5.0)]);
        let src = src_polys(vec![(rect(0.0, 0.0, 5.0, 5.0), 10.0, "a")]);

        let (out_keep, layer_keep) = run(json!({
            "polygons": zones, "input": src, "fields": "pop:sum", "keep_all": true
        }));
        assert_eq!(layer_keep.len(), 2);
        assert_eq!(out_keep.outputs["zones_with_matches"], json!(1));

        let (_, layer_drop) = run(json!({
            "polygons": zones, "input": src, "fields": "pop:sum", "keep_all": false
        }));
        assert_eq!(layer_drop.len(), 1);
    }

    /// group_field emits one row per distinct category within each zone.
    #[test]
    fn group_field_splits_rows() {
        let zones = zones_layer(vec![rect(0.0, 0.0, 10.0, 10.0)]);
        let src = src_polys(vec![
            (rect(0.0, 0.0, 2.0, 2.0), 10.0, "res"),
            (rect(3.0, 3.0, 2.0, 2.0), 20.0, "com"),
            (rect(6.0, 6.0, 2.0, 2.0), 30.0, "res"),
        ]);
        let (_, layer) = run(json!({
            "polygons": zones, "input": src, "fields": "pop:sum", "group_field": "cat"
        }));
        assert_eq!(layer.len(), 2, "expected one row per category");
        let gi = layer.schema.field_index("group").unwrap();
        let mut by_group: BTreeMap<String, f64> = BTreeMap::new();
        for (r, f) in layer.features.iter().enumerate() {
            let g = match &f.attributes[gi] {
                FieldValue::Text(s) => s.clone(),
                _ => String::new(),
            };
            by_group.insert(g, num(&layer, r, "pop_sum"));
        }
        assert!((by_group["res"] - 40.0).abs() < 1e-9);
        assert!((by_group["com"] - 20.0).abs() < 1e-9);
    }

    /// mean is weighted by overlap, not a plain average of the values.
    #[test]
    fn mean_is_weighted() {
        // Zone [0,10]x[0,10]. A: fully inside, value 10. B: 10% inside, value 100.
        let zones = zones_layer(vec![rect(0.0, 0.0, 10.0, 10.0)]);
        let src = src_polys(vec![
            (rect(0.0, 0.0, 10.0, 10.0), 10.0, "a"),
            (rect(9.0, 0.0, 10.0, 10.0), 100.0, "b"),
        ]);
        let (_, layer) = run(json!({
            "polygons": zones, "input": src, "fields": "pop:mean"
        }));
        // weights: A = 1.0, B = 0.1 -> (10*1 + 100*0.1) / 1.1 = 18.18...
        let m = num(&layer, 0, "pop_mean");
        assert!((m - 20.0 / 1.1).abs() < 1e-6, "weighted mean was {m}");
    }

    #[test]
    fn rejects_bad_parameters() {
        let zones = zones_layer(vec![rect(0.0, 0.0, 1.0, 1.0)]);
        let src = src_polys(vec![(rect(0.0, 0.0, 1.0, 1.0), 1.0, "a")]);
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SummarizeWithinTool.validate(&args).is_err()
        };
        assert!(bad(json!({ "input": src })));
        assert!(bad(json!({ "polygons": zones })));
        assert!(bad(
            json!({ "polygons": zones, "input": src, "fields": "pop:bogus" })
        ));
        assert!(bad(
            json!({ "polygons": zones, "input": src, "keep_all": "maybe" })
        ));

        // Unknown field is caught at run time (needs the layer to check).
        let args: ToolArgs = serde_json::from_value(
            json!({ "polygons": zones, "input": src, "fields": "nope:sum" }),
        )
        .unwrap();
        assert!(matches!(
            SummarizeWithinTool.run(&args, &ctx()).unwrap_err(),
            ToolError::Validation(_)
        ));
    }
}
