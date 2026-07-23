//! GeoLibre tool: assign a shared group id to features within a proximity of
//! one another (single-linkage spatial connected components).
//!
//! Pure-Rust counterpart of ArcGIS GeoAnalytics' *Group By Proximity*. Unlike
//! the density-clustering tools (`hdbscan`, `optics`, bundled `dbscan`) which
//! need density parameters and drop noise, this simply connects features whose
//! geometries lie within a distance (optionally sharing an attribute value) and
//! tags every original feature with a `GROUP_ID`, preserving its geometry.
//!
//! Edges are formed when the minimum distance between two features is within the
//! threshold (`near`) or zero (`intersects`, i.e. touching). Components are
//! resolved with union-find; group ids are assigned deterministically in order
//! of each component's smallest feature index.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct GroupByProximityTool;

impl Tool for GroupByProximityTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "group_by_proximity",
            display_name: "Group By Proximity",
            summary: "Tag features with a GROUP_ID for single-linkage spatial connected components — features within a distance ('near') or touching ('intersects'), optionally sharing an attribute value — like ArcGIS Group By Proximity.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (any geometry type).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "relationship",
                    description: "'near' (within 'spatial_near_distance', default) or 'intersects' (touching, distance 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "spatial_near_distance",
                    description: "Search distance in map units for 'near' (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "attribute_field",
                    description: "Optional field; features only group when they share this attribute value.",
                    required: false,
                },
                ToolParamSpec {
                    name: "group_field",
                    description: "Name of the output group-id field (default 'GROUP_ID').",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_relationship(args)?;
        parse_optional_f64(args, "spatial_near_distance")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let rel = parse_relationship(args)?;
        let distance = match rel {
            Relationship::Intersects => 0.0,
            Relationship::Near => parse_optional_f64(args, "spatial_near_distance")?.unwrap_or(0.0),
        };
        if distance < 0.0 {
            return Err(ToolError::Validation(
                "'spatial_near_distance' must be non-negative".to_string(),
            ));
        }
        let group_field = parse_optional_str(args, "group_field")?.unwrap_or("GROUP_ID");
        let attr_field = parse_optional_str(args, "attribute_field")?;

        let layer = load_input_layer(input)?;
        let n = layer.features.len();
        if n == 0 {
            return Err(ToolError::Execution("input has no features".to_string()));
        }
        let aidx = match attr_field {
            Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                ToolError::Validation(format!("attribute_field '{f}' not found"))
            })?),
            None => None,
        };

        // Segment sets + bboxes per feature.
        let mut segs: Vec<Vec<[Pt; 2]>> = Vec::with_capacity(n);
        let mut bboxes: Vec<Option<[f64; 4]>> = Vec::with_capacity(n);
        for f in &layer.features {
            let s = f
                .geometry
                .as_ref()
                .map(feature_segments)
                .unwrap_or_default();
            let bb = f
                .geometry
                .as_ref()
                .and_then(|g| g.bbox())
                .map(|b| [b.min_x, b.min_y, b.max_x, b.max_y]);
            segs.push(s);
            bboxes.push(bb);
        }

        // Union-find over feature pairs within the threshold.
        let mut uf = UnionFind::new(n);
        ctx.progress.info(&format!("grouping {n} feature(s)"));
        for i in 0..n {
            let Some(bi) = bboxes[i] else { continue };
            if segs[i].is_empty() {
                continue;
            }
            for j in (i + 1)..n {
                let Some(bj) = bboxes[j] else { continue };
                if segs[j].is_empty() {
                    continue;
                }
                // Bbox pre-filter: skip pairs whose bboxes are farther than d.
                if bbox_gap(&bi, &bj) > distance {
                    continue;
                }
                if let (Some(a), Some(b)) = (aidx, aidx) {
                    let va = layer.features[i].attributes.get(a);
                    let vb = layer.features[j].attributes.get(b);
                    if va != vb {
                        continue;
                    }
                }
                if segs_min_distance(&segs[i], &segs[j]) <= distance {
                    uf.union(i, j);
                }
            }
        }

        // Deterministic group ids: root's first-seen order.
        let mut root_to_id: BTreeMap<usize, i64> = BTreeMap::new();
        let mut next_id = 0i64;
        let mut ids = vec![0i64; n];
        for (i, slot) in ids.iter_mut().enumerate() {
            let root = uf.find(i);
            let id = *root_to_id.entry(root).or_insert_with(|| {
                let v = next_id;
                next_id += 1;
                v
            });
            *slot = id;
        }

        // Build output: copy schema + features, append the group field.
        let mut out = Layer::new("grouped");
        if let Some(gt) = layer.geom_type {
            out = out.with_geom_type(gt);
        }
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for fd in layer.schema.fields() {
            out.add_field(fd.clone());
        }
        out.add_field(FieldDef::new(group_field, FieldType::Integer));
        for (i, feat) in layer.features.iter().enumerate() {
            let orig: Vec<(String, FieldValue)> = layer
                .schema
                .fields()
                .iter()
                .enumerate()
                .map(|(fi, fd)| {
                    (
                        fd.name.clone(),
                        feat.attributes.get(fi).cloned().unwrap_or(FieldValue::Null),
                    )
                })
                .collect();
            let mut attrs: Vec<(&str, FieldValue)> =
                orig.iter().map(|(k, v)| (k.as_str(), v.clone())).collect();
            attrs.push((group_field, FieldValue::Integer(ids[i])));
            out.add_feature(feat.geometry.clone(), &attrs)
                .map_err(|e| ToolError::Execution(format!("failed adding feature: {e}")))?;
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("group_count".to_string(), json!(next_id));
        Ok(ToolRunResult { outputs })
    }
}

// ── Geometry helpers (plain (x, y) tuples) ──────────────────────────────────

type Pt = (f64, f64);

/// Flattens a geometry into a set of segments. Points and single-vertex parts
/// become degenerate segments (both endpoints equal) so distance tests still work.
fn feature_segments(geom: &Geometry) -> Vec<[Pt; 2]> {
    fn push_ring(coords: &[Coord], out: &mut Vec<[Pt; 2]>) {
        if coords.len() == 1 {
            let p = (coords[0].x, coords[0].y);
            out.push([p, p]);
        }
        for w in coords.windows(2) {
            out.push([(w[0].x, w[0].y), (w[1].x, w[1].y)]);
        }
    }
    let mut out = Vec::new();
    match geom {
        Geometry::Point(c) => out.push([(c.x, c.y), (c.x, c.y)]),
        Geometry::MultiPoint(cs) => {
            for c in cs {
                out.push([(c.x, c.y), (c.x, c.y)]);
            }
        }
        Geometry::LineString(cs) => push_ring(cs, &mut out),
        Geometry::MultiLineString(ls) => {
            for l in ls {
                push_ring(l, &mut out);
            }
        }
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
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                out.extend(feature_segments(g));
            }
        }
    }
    out
}

/// Minimum distance between two segment sets (early-exit at 0).
fn segs_min_distance(a: &[[Pt; 2]], b: &[[Pt; 2]]) -> f64 {
    let mut best = f64::INFINITY;
    for sa in a {
        for sb in b {
            let d = seg_seg_distance(sa, sb);
            if d < best {
                best = d;
                if best == 0.0 {
                    return 0.0;
                }
            }
        }
    }
    best
}

/// Distance between two segments (0 if they intersect).
fn seg_seg_distance(a: &[Pt; 2], b: &[Pt; 2]) -> f64 {
    if segments_intersect(a, b) {
        return 0.0;
    }
    let d1 = point_seg_distance(a[0], b);
    let d2 = point_seg_distance(a[1], b);
    let d3 = point_seg_distance(b[0], a);
    let d4 = point_seg_distance(b[1], a);
    d1.min(d2).min(d3).min(d4)
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
    let (px, py) = (ax + t * dx, ay + t * dy);
    (p.0 - px).hypot(p.1 - py)
}

fn segments_intersect(a: &[Pt; 2], b: &[Pt; 2]) -> bool {
    let o = |p: Pt, q: Pt, r: Pt| (q.0 - p.0) * (r.1 - p.1) - (q.1 - p.1) * (r.0 - p.0);
    let d1 = o(b[0], b[1], a[0]);
    let d2 = o(b[0], b[1], a[1]);
    let d3 = o(a[0], a[1], b[0]);
    let d4 = o(a[0], a[1], b[1]);
    ((d1 > 0.0) != (d2 > 0.0)) && ((d3 > 0.0) != (d4 > 0.0))
}

/// Gap between two bboxes ([min_x, min_y, max_x, max_y]); 0 if they overlap.
fn bbox_gap(a: &[f64; 4], b: &[f64; 4]) -> f64 {
    let dx = (b[0] - a[2]).max(a[0] - b[2]).max(0.0);
    let dy = (b[1] - a[3]).max(a[1] - b[3]).max(0.0);
    dx.hypot(dy)
}

// ── Union-Find ──────────────────────────────────────────────────────────────

struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }
    fn find(&mut self, x: usize) -> usize {
        let mut root = x;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        let mut cur = x;
        while self.parent[cur] != root {
            let next = self.parent[cur];
            self.parent[cur] = root;
            cur = next;
        }
        root
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        if self.rank[ra] < self.rank[rb] {
            self.parent[ra] = rb;
        } else if self.rank[ra] > self.rank[rb] {
            self.parent[rb] = ra;
        } else {
            self.parent[rb] = ra;
            self.rank[ra] += 1;
        }
    }
}

// ── Params ──────────────────────────────────────────────────────────────────

enum Relationship {
    Near,
    Intersects,
}

fn parse_relationship(args: &ToolArgs) -> Result<Relationship, ToolError> {
    match args
        .get("relationship")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some("near") => Ok(Relationship::Near),
        Some("intersects") => Ok(Relationship::Intersects),
        Some(o) => Err(ToolError::Validation(format!(
            "'relationship' must be 'near' or 'intersects', got '{o}'"
        ))),
    }
}

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

    fn points(pts: &[(f64, f64)], attrs: &[i64]) -> String {
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        let with_attr = !attrs.is_empty();
        if with_attr {
            l.add_field(FieldDef::new("cat", FieldType::Integer));
        }
        for (i, (x, y)) in pts.iter().enumerate() {
            let a: &[(&str, FieldValue)] = if with_attr {
                &[("cat", FieldValue::Integer(attrs[i]))]
            } else {
                &[]
            };
            l.add_feature(Some(Geometry::point(*x, *y)), a).unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = GroupByProximityTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn group_ids(l: &Layer) -> Vec<i64> {
        let idx = l.schema.field_index("GROUP_ID").unwrap();
        l.iter()
            .map(|f| f.attributes[idx].as_i64().unwrap())
            .collect()
    }

    #[test]
    fn two_close_one_far_forms_two_groups() {
        // (0,0)&(1,0) within 1.5; (100,100) alone.
        let input = points(&[(0.0, 0.0), (1.0, 0.0), (100.0, 100.0)], &[]);
        let (out, layer) = run(json!({ "input": input, "spatial_near_distance": 1.5 }));
        assert_eq!(out.outputs["group_count"], json!(2));
        let ids = group_ids(&layer);
        assert_eq!(ids[0], ids[1]);
        assert_ne!(ids[0], ids[2]);
    }

    #[test]
    fn transitive_chaining_links_a_line_of_points() {
        // Single-linkage: 0-1-2-3 each 1 apart, threshold 1.1 -> all one group.
        let input = points(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0)], &[]);
        let (out, layer) = run(json!({ "input": input, "spatial_near_distance": 1.1 }));
        assert_eq!(out.outputs["group_count"], json!(1));
        let ids = group_ids(&layer);
        assert!(ids.iter().all(|&g| g == ids[0]));
    }

    #[test]
    fn attribute_field_splits_groups() {
        // Two coincident-ish points but different category -> stay separate.
        let input = points(&[(0.0, 0.0), (0.5, 0.0)], &[1, 2]);
        let (out, _l) = run(json!({
            "input": input, "spatial_near_distance": 2.0, "attribute_field": "cat"
        }));
        assert_eq!(out.outputs["group_count"], json!(2));
    }

    #[test]
    fn zero_distance_isolates_all() {
        let input = points(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)], &[]);
        let (out, _l) = run(json!({ "input": input, "spatial_near_distance": 0.0 }));
        assert_eq!(out.outputs["group_count"], json!(3));
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            GroupByProximityTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "p.shp", "relationship": "bogus" })).is_err());
        assert!(bad(json!({ "input": "p.shp" })).is_ok());
    }
}
