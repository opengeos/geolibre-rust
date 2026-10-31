//! GeoLibre tool: period-over-period incident change per area.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Summarize Percent Change* (Crime
//! Analysis and Safety). GeoLibre already has the incident family —
//! `detect_incidents`, `eighty_twenty_analysis`, `collect_events`,
//! `trace_proximity_events`, `find_dwell_locations` — and the whole hot-spot
//! suite, but nothing performed the plainest question in that domain:
//! "burglaries this month versus last month, by beat".
//!
//! `summarize_within` aggregates one point layer into polygons, so doing this
//! by hand means running it twice and joining. `emerging_hot_spot_analysis`
//! answers a far more elaborate question (space-time trend significance) and is
//! the wrong instrument for a two-period comparison.
//!
//! One deliberate correctness choice: when the previous count is zero,
//! `PERCENT_CHANGE` is **null**, not 0 and not 100. A ratio against zero has no
//! finite value, and reporting one is the standard way this class of tool
//! misleads. Those areas are labelled `new` instead, which is the information
//! the analyst actually wants.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::vector_common::{
    geometry_contains_point, load_input_layer, parse_optional_str, write_or_store_layer,
};

pub struct SummarizePercentChangeTool;

impl Tool for SummarizePercentChangeTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "summarize_percent_change",
            display_name: "Summarize Percent Change",
            summary: "Count incidents from a current and a previous period within each area feature and report the count difference, percent change and a change class. Percent change is null (class 'new') where the previous count is zero. Like ArcGIS Summarize Percent Change.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Area features (beats, tracts, grid cells) to summarize into.",
                    required: true,
                },
                ToolParamSpec {
                    name: "current_features",
                    description: "Incident point features for the current period.",
                    required: true,
                },
                ToolParamSpec {
                    name: "previous_features",
                    description: "Incident point features for the previous period.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output path for the areas with change fields. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "search_radius",
                    description: "Optional distance in map units; incidents within this distance of an area also count, not only those inside it. Note that a positive radius makes catchments overlap, so an incident can be counted by more than one area (the reported corpus totals stay deduplicated).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "current_features")?;
        require_str(args, "previous_features")?;
        if let Some(r) = parse_optional_f64(args, "search_radius")? {
            if r < 0.0 {
                return Err(ToolError::Validation(
                    "'search_radius' must be non-negative".to_string(),
                ));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let cur_path = require_str(args, "current_features")?;
        let prev_path = require_str(args, "previous_features")?;
        let output = parse_optional_str(args, "output")?;
        let radius = parse_optional_f64(args, "search_radius")?.unwrap_or(0.0);

        let areas = load_input_layer(input)?;
        if areas.features.is_empty() {
            return Err(ToolError::Execution("input has no features".to_string()));
        }
        let current = collect_points(&load_input_layer(cur_path)?);
        let previous = collect_points(&load_input_layer(prev_path)?);
        ctx.progress.info(&format!(
            "summarizing {} current / {} previous incident(s) into {} area(s)",
            current.len(),
            previous.len(),
            areas.features.len()
        ));

        let mut out = Layer::new("percent_change");
        if let Some(gt) = areas.geom_type {
            out = out.with_geom_type(gt);
        }
        if let Some(epsg) = areas.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for fd in areas.schema.fields() {
            out.add_field(fd.clone());
        }
        out.add_field(FieldDef::new("CURRENT_COUNT", FieldType::Integer));
        out.add_field(FieldDef::new("PREVIOUS_COUNT", FieldType::Integer));
        out.add_field(FieldDef::new("COUNT_DIFF", FieldType::Integer));
        out.add_field(FieldDef::new("PERCENT_CHANGE", FieldType::Float));
        out.add_field(FieldDef::new("CHANGE_CLASS", FieldType::Text));

        let names: Vec<String> = areas
            .schema
            .fields()
            .iter()
            .map(|f| f.name.clone())
            .collect();
        // Per-area counts summed would double count any incident that falls in
        // two overlapping catchments — guaranteed once `search_radius` > 0
        // expands each area — and the overall percent change would then be
        // computed from inflated denominators. Corpus totals therefore come from
        // the deduplicated point lists, and the summed-per-area figures are
        // reported separately.
        let (mut sum_cur, mut sum_prev) = (0i64, 0i64);
        let mut class_counts: BTreeMap<&str, usize> = BTreeMap::new();

        for feat in areas.iter() {
            let (cur_n, prev_n) = match &feat.geometry {
                Some(g) => (
                    count_in(g, &current, radius),
                    count_in(g, &previous, radius),
                ),
                None => (0, 0),
            };
            sum_cur += cur_n;
            sum_prev += prev_n;
            let diff = cur_n - prev_n;
            let (pct, class) = classify(cur_n, prev_n);
            *class_counts.entry(class).or_default() += 1;

            let mut attrs: Vec<(&str, FieldValue)> = names
                .iter()
                .enumerate()
                .map(|(fi, nm)| {
                    (
                        nm.as_str(),
                        feat.attributes.get(fi).cloned().unwrap_or(FieldValue::Null),
                    )
                })
                .collect();
            attrs.push(("CURRENT_COUNT", FieldValue::Integer(cur_n)));
            attrs.push(("PREVIOUS_COUNT", FieldValue::Integer(prev_n)));
            attrs.push(("COUNT_DIFF", FieldValue::Integer(diff)));
            attrs.push((
                "PERCENT_CHANGE",
                pct.map_or(FieldValue::Null, FieldValue::Float),
            ));
            attrs.push(("CHANGE_CLASS", FieldValue::Text(class.to_string())));
            out.add_feature(feat.geometry.clone(), &attrs)
                .map_err(|e| ToolError::Execution(format!("failed adding feature: {e}")))?;
        }

        let n_areas = areas.features.len();
        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("area_count".to_string(), json!(n_areas));
        let tot_cur = current.len() as i64;
        let tot_prev = previous.len() as i64;
        outputs.insert("current_total".to_string(), json!(tot_cur));
        outputs.insert("previous_total".to_string(), json!(tot_prev));
        outputs.insert("total_diff".to_string(), json!(tot_cur - tot_prev));
        if tot_prev > 0 {
            outputs.insert(
                "total_percent_change".to_string(),
                json!((tot_cur - tot_prev) as f64 / tot_prev as f64 * 100.0),
            );
        }
        // Summed per-area counts: equal to the corpus totals only when the
        // catchments do not overlap and every incident falls inside one.
        outputs.insert("current_in_areas".to_string(), json!(sum_cur));
        outputs.insert("previous_in_areas".to_string(), json!(sum_prev));
        for (k, v) in class_counts {
            outputs.insert(format!("class_{k}"), json!(v));
        }
        Ok(ToolRunResult { outputs })
    }
}

/// Percent change and change class for one area.
///
/// Returns `None` for the percent when the previous count is zero — see the
/// module docs for why that is deliberately not 0 or 100.
fn classify(cur: i64, prev: i64) -> (Option<f64>, &'static str) {
    match (cur, prev) {
        (0, 0) => (None, "no_change"),
        (c, 0) if c > 0 => (None, "new"),
        (0, p) if p > 0 => (Some(-100.0), "eliminated"),
        (c, p) => {
            let pct = (c - p) as f64 / p as f64 * 100.0;
            let class = if c > p {
                "increase"
            } else if c < p {
                "decrease"
            } else {
                "no_change"
            };
            (Some(pct), class)
        }
    }
}

/// Flattens an incident layer to (x, y) pairs with a bbox for pre-filtering.
fn collect_points(layer: &Layer) -> Vec<(f64, f64)> {
    let mut out = Vec::with_capacity(layer.features.len());
    for f in layer.iter() {
        let Some(g) = &f.geometry else { continue };
        match g {
            Geometry::Point(c) => out.push((c.x, c.y)),
            Geometry::MultiPoint(cs) => out.extend(cs.iter().map(|c| (c.x, c.y))),
            // Non-point incident geometry is summarized at its centroid so a
            // mixed layer still contributes rather than being silently dropped.
            other => {
                if let Some(c) = centroid(other) {
                    out.push(c);
                }
            }
        }
    }
    out
}

/// Counts points inside `geom`, or within `radius` of it when radius > 0.
fn count_in(geom: &Geometry, pts: &[(f64, f64)], radius: f64) -> i64 {
    let Some(bb) = geom.bbox() else { return 0 };
    let (min_x, min_y) = (bb.min_x - radius, bb.min_y - radius);
    let (max_x, max_y) = (bb.max_x + radius, bb.max_y + radius);
    let rings = radius > 0.0;
    let segs = if rings { boundary_segments(geom) } else { Vec::new() };

    let mut n = 0i64;
    for &(x, y) in pts {
        if x < min_x || x > max_x || y < min_y || y > max_y {
            continue;
        }
        // Inside the area, or (with a radius) close enough to its boundary.
        let hit = geometry_contains_point(geom, x, y)
            || (rings && segs.iter().any(|s| point_seg_distance((x, y), s) <= radius));
        if hit {
            n += 1;
        }
    }
    n
}

type Pt = (f64, f64);

fn boundary_segments(geom: &Geometry) -> Vec<[Pt; 2]> {
    fn push_ring(coords: &[Coord], out: &mut Vec<[Pt; 2]>) {
        for w in coords.windows(2) {
            out.push([(w[0].x, w[0].y), (w[1].x, w[1].y)]);
        }
    }
    let mut out = Vec::new();
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            push_ring(&exterior.0, &mut out);
            for r in interiors {
                push_ring(&r.0, &mut out);
            }
        }
        Geometry::MultiPolygon(ps) => {
            for (e, hs) in ps {
                push_ring(&e.0, &mut out);
                for r in hs {
                    push_ring(&r.0, &mut out);
                }
            }
        }
        Geometry::LineString(cs) => push_ring(cs, &mut out),
        Geometry::MultiLineString(ls) => {
            for l in ls {
                push_ring(l, &mut out);
            }
        }
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                out.extend(boundary_segments(g));
            }
        }
        _ => {}
    }
    out
}

fn point_seg_distance(p: Pt, s: &[Pt; 2]) -> f64 {
    let (ax, ay) = s[0];
    let (bx, by) = s[1];
    let (dx, dy) = (bx - ax, by - ay);
    let len2 = dx * dx + dy * dy;
    if len2 <= f64::EPSILON {
        return (p.0 - ax).hypot(p.1 - ay);
    }
    let t = (((p.0 - ax) * dx + (p.1 - ay) * dy) / len2).clamp(0.0, 1.0);
    (p.0 - (ax + t * dx)).hypot(p.1 - (ay + t * dy))
}

/// Vertex-average fallback centroid for non-point incident geometry.
fn centroid(g: &Geometry) -> Option<Pt> {
    let coords = g.all_coords();
    if coords.is_empty() {
        return None;
    }
    let n = coords.len() as f64;
    Some((
        coords.iter().map(|c| c.x).sum::<f64>() / n,
        coords.iter().map(|c| c.y).sum::<f64>() / n,
    ))
}

// ── Params ──────────────────────────────────────────────────────────────────

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
    use wbvector::GeometryType;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Two 10x10 areas side by side: A at x 0-10, B at x 20-30.
    fn areas() -> String {
        let mut l = Layer::new("areas")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        for (name, x0) in [("A", 0.0), ("B", 20.0)] {
            let ring = vec![
                Coord::xy(x0, 0.0),
                Coord::xy(x0 + 10.0, 0.0),
                Coord::xy(x0 + 10.0, 10.0),
                Coord::xy(x0, 10.0),
                Coord::xy(x0, 0.0),
            ];
            l.add_feature(
                Some(Geometry::polygon(ring, vec![])),
                &[("name", FieldValue::Text(name.to_string()))],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn pts(p: &[(f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        for (x, y) in p {
            l.add_feature(Some(Geometry::point(*x, *y)), &[]).unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SummarizePercentChangeTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn col(l: &Layer, name: &str) -> Vec<FieldValue> {
        let i = l.schema.field_index(name).unwrap();
        l.iter().map(|f| f.attributes[i].clone()).collect()
    }

    #[test]
    fn counts_and_percent_change_per_area() {
        // A: 2 -> 3 (+50%). B: 4 -> 2 (-50%).
        let (out, layer) = run(json!({
            "input": areas(),
            "previous_features": pts(&[(1.0, 1.0), (2.0, 2.0), (21.0, 1.0), (22.0, 2.0), (23.0, 3.0), (24.0, 4.0)]),
            "current_features": pts(&[(1.0, 1.0), (2.0, 2.0), (3.0, 3.0), (21.0, 1.0), (22.0, 2.0)]),
        }));
        assert_eq!(out.outputs["current_total"], json!(5));
        assert_eq!(out.outputs["previous_total"], json!(6));
        let pct = col(&layer, "PERCENT_CHANGE");
        assert_eq!(pct[0].as_f64(), Some(50.0));
        assert_eq!(pct[1].as_f64(), Some(-50.0));
        let cls = col(&layer, "CHANGE_CLASS");
        assert_eq!(cls[0].as_str(), Some("increase"));
        assert_eq!(cls[1].as_str(), Some("decrease"));
    }

    #[test]
    fn zero_previous_yields_null_percent_and_new_class() {
        // This is the correctness point: no finite percent exists here.
        let (_o, layer) = run(json!({
            "input": areas(),
            "previous_features": pts(&[(21.0, 1.0)]),
            "current_features": pts(&[(1.0, 1.0), (2.0, 2.0), (21.0, 1.0)]),
        }));
        let pct = col(&layer, "PERCENT_CHANGE");
        let cls = col(&layer, "CHANGE_CLASS");
        assert!(pct[0].is_null());
        assert_eq!(cls[0].as_str(), Some("new"));
    }

    #[test]
    fn dropping_to_zero_is_minus_one_hundred_and_eliminated() {
        let (_o, layer) = run(json!({
            "input": areas(),
            "previous_features": pts(&[(1.0, 1.0), (2.0, 2.0)]),
            "current_features": pts(&[(21.0, 1.0)]),
        }));
        let pct = col(&layer, "PERCENT_CHANGE");
        let cls = col(&layer, "CHANGE_CLASS");
        assert_eq!(pct[0].as_f64(), Some(-100.0));
        assert_eq!(cls[0].as_str(), Some("eliminated"));
    }

    #[test]
    fn empty_both_periods_is_no_change_with_null_percent() {
        let (_o, layer) = run(json!({
            "input": areas(),
            "previous_features": pts(&[(21.0, 1.0)]),
            "current_features": pts(&[(21.0, 1.0)]),
        }));
        let pct = col(&layer, "PERCENT_CHANGE");
        let cls = col(&layer, "CHANGE_CLASS");
        assert!(pct[0].is_null());
        assert_eq!(cls[0].as_str(), Some("no_change"));
    }

    #[test]
    fn search_radius_captures_nearby_incidents() {
        // A point 2 units outside area A's edge. 'current_total' is the corpus
        // size (always 1 here); what the radius changes is how many areas
        // actually count it, reported as 'current_in_areas'.
        let outside = pts(&[(12.0, 5.0)]);
        let (tight, _l) = run(json!({
            "input": areas(), "previous_features": outside.clone(),
            "current_features": outside.clone(),
        }));
        assert_eq!(tight.outputs["current_in_areas"], json!(0));
        assert_eq!(tight.outputs["current_total"], json!(1));
        let (loose, _l) = run(json!({
            "input": areas(), "previous_features": outside.clone(),
            "current_features": outside, "search_radius": 3.0,
        }));
        assert_eq!(loose.outputs["current_in_areas"], json!(1));
    }

    #[test]
    fn corpus_totals_do_not_double_count_overlapping_catchments() {
        // With a radius wide enough that both areas claim the same incident,
        // the per-area sum counts it twice but the corpus total must not.
        let mid = pts(&[(15.0, 5.0)]);
        let (out, _l) = run(json!({
            "input": areas(), "previous_features": mid.clone(),
            "current_features": mid, "search_radius": 6.0,
        }));
        assert_eq!(out.outputs["current_in_areas"], json!(2));
        assert_eq!(out.outputs["current_total"], json!(1));
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SummarizePercentChangeTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.shp", "current_features": "c.shp" })).is_err());
        assert!(bad(json!({
            "input": "a.shp", "current_features": "c.shp", "previous_features": "p.shp",
            "search_radius": -5
        }))
        .is_err());
        assert!(bad(json!({
            "input": "a.shp", "current_features": "c.shp", "previous_features": "p.shp"
        }))
        .is_ok());
    }
}
