//! GeoLibre tool: merge matched divided-carriageway pairs into single
//! centerlines.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Merge Divided Roads* (Cartography).
//!
//! This deliberately overlaps a shipped tool, and the issue that requested it
//! said so. GeoLibre already has `collapse_dual_lines_to_centerline`, plus
//! `thin_road_network`, `collapse_road_detail` and `resolve_road_conflicts`.
//! The distinction ArcGIS draws — and the reason this still earns its place —
//! is that Collapse Dual Lines To Centerline is a general dual-line collapse
//! (its natural home is hydrography), whereas Merge Divided Roads is
//! carriageway-aware:
//!
//!   * it pairs lanes using a `merge_field` (road class/name), so only genuine
//!     carriageway pairs merge;
//!   * it **preserves intersection connectivity** instead of collapsing
//!     junctions, which is what keeps the result routable;
//!   * it honours a `character_field` protecting features that must not merge;
//!   * it emits **displacement features**, which feed the shipped
//!     `propagate_displacement`.
//!
//! Junction preservation is the part worth spelling out: nodes where three or
//! more lines meet are pinned, and after a pair collapses to its midline the
//! tool re-attaches each pinned node with a short connector stub. Without that
//! step the merged network would be geometrically prettier and topologically
//! broken.

use std::collections::BTreeMap;
use std::collections::HashMap;

use geo::{Area, BooleanOps, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Coordinate quantum for node identity.
const SNAP: f64 = 1e-6;
/// How far from anti-parallel two carriageways may run and still pair.
const BEARING_TOLERANCE_DEG: f64 = 35.0;
/// Samples used when matching two carriageways along their length.
const SAMPLES: usize = 24;

pub struct MergeDividedRoadsTool;

impl Tool for MergeDividedRoadsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "merge_divided_roads",
            display_name: "Merge Divided Roads",
            summary: "Replace matched pairs of divided road lanes with single centerlines while preserving intersection connectivity, pairing only lanes that share a merge field and run anti-parallel within a separation distance, and emitting displacement polygons for propagate_displacement. Like ArcGIS Merge Divided Roads.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Road centerline features.",
                    required: true,
                },
                ToolParamSpec {
                    name: "merge_field",
                    description: "Field whose matching values make two lines candidate carriageway pairs (e.g. road name or class).",
                    required: true,
                },
                ToolParamSpec {
                    name: "merge_distance",
                    description: "Maximum separation between paired carriageways, in map units.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output merged centerline path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_displacement_features",
                    description: "Optional polygons showing where geometry moved, for propagate_displacement.",
                    required: false,
                },
                ToolParamSpec {
                    name: "character_field",
                    description: "Optional field; features with a non-zero / non-empty value are excluded from merging.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_table",
                    description: "Optional lineage table mapping output features back to their source features.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "merge_field")?;
        let d = parse_optional_f64(args, "merge_distance")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'merge_distance'".to_string())
        })?;
        if d <= 0.0 {
            return Err(ToolError::Validation(
                "'merge_distance' must be positive".to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let merge_field = require_str(args, "merge_field")?;
        let merge_distance = parse_optional_f64(args, "merge_distance")?.unwrap_or(0.0);
        let output = parse_optional_str(args, "output")?;

        let layer = load_input_layer(input)?;
        let n = layer.features.len();
        if n == 0 {
            return Err(ToolError::Execution("input has no features".to_string()));
        }
        let m_idx = layer
            .schema
            .field_index(merge_field)
            .ok_or_else(|| ToolError::Validation(format!("merge_field '{merge_field}' not found")))?;
        let c_idx = match parse_optional_str(args, "character_field")? {
            Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                ToolError::Validation(format!("character_field '{f}' not found"))
            })?),
            None => None,
        };

        // Flatten to single paths, remembering the source feature.
        let mut paths: Vec<(usize, Vec<Coord>)> = Vec::new();
        for (fid, f) in layer.iter().enumerate() {
            let Some(g) = &f.geometry else { continue };
            for p in line_paths(g) {
                paths.push((fid, p));
            }
        }
        if paths.is_empty() {
            return Err(ToolError::Execution(
                "input contains no line geometry".to_string(),
            ));
        }

        // Node degrees over EVERY vertex, not just endpoints: side streets
        // usually T into a carriageway mid-span, and a tool that only pinned
        // endpoint junctions would sever exactly those. A node touched by 3+
        // distinct paths is a junction that must survive the merge.
        let mut degree: HashMap<(i64, i64), usize> = HashMap::new();
        for (_, p) in &paths {
            // A HashSet, not a Vec: densified road geometry has thousands of
            // vertices per feature, and `contains` on a Vec makes this quadratic.
            let mut seen: std::collections::HashSet<(i64, i64)> =
                std::collections::HashSet::new();
            for c in p {
                let k = node_key(c);
                if seen.insert(k) {
                    *degree.entry(k).or_insert(0) += 1;
                }
            }
        }

        // Candidate pairing.
        let protected = |i: usize| -> bool {
            c_idx.is_some_and(|ci| {
                let v = layer.features[paths[i].0].attributes.get(ci);
                match v {
                    Some(FieldValue::Null) | None => false,
                    Some(FieldValue::Integer(x)) => *x != 0,
                    Some(FieldValue::Float(x)) => *x != 0.0,
                    Some(FieldValue::Boolean(b)) => *b,
                    Some(FieldValue::Text(s)) => !s.trim().is_empty(),
                    _ => true,
                }
            })
        };
        // Merge keys and bboxes are loop-invariant: computing the key inside the
        // inner loop re-allocated a String per comparison, and every surviving
        // candidate then ran 25 samples x a full distance-to-path scan.
        let keys: Vec<String> = (0..paths.len())
            .map(|i| key_str(layer.features[paths[i].0].attributes.get(m_idx)))
            .collect();
        let bboxes: Vec<[f64; 4]> = paths.iter().map(|(_, p)| path_bbox(p)).collect();

        let mut partner: Vec<Option<usize>> = vec![None; paths.len()];
        for i in 0..paths.len() {
            if partner[i].is_some() || protected(i) {
                continue;
            }
            let mut best: Option<(usize, f64)> = None;
            for j in (i + 1)..paths.len() {
                if partner[j].is_some() || protected(j) {
                    continue;
                }
                if keys[i] != keys[j] {
                    continue;
                }
                // Cheap reject before the sampling test.
                if bbox_gap(&bboxes[i], &bboxes[j]) > merge_distance {
                    continue;
                }
                let Some(sep) = carriageway_separation(&paths[i].1, &paths[j].1, merge_distance)
                else {
                    continue;
                };
                if best.is_none_or(|(_, bs)| sep < bs) {
                    best = Some((j, sep));
                }
            }
            if let Some((j, _)) = best {
                partner[i] = Some(j);
                partner[j] = Some(i);
            }
        }

        let mut out = Layer::new("merged_roads").with_geom_type(GeometryType::LineString);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for fd in layer.schema.fields() {
            out.add_field(fd.clone());
        }
        out.add_field(FieldDef::new("MERGE_ROLE", FieldType::Text));
        out.add_field(FieldDef::new("SRC_FID_A", FieldType::Integer));
        out.add_field(FieldDef::new("SRC_FID_B", FieldType::Integer));
        let names: Vec<String> = layer
            .schema
            .fields()
            .iter()
            .map(|f| f.name.clone())
            .collect();

        // No declared geometry type: a self-intersecting sweep normalises into
        // more than one part, so `multipolygon_to_geometry` can return a
        // MultiPolygon, which a Polygon-typed layer would reject or coerce.
        let mut displacement = Layer::new("displacement");
        if let Some(epsg) = layer.crs_epsg() {
            displacement = displacement.with_crs_epsg(epsg);
        }
        displacement.add_field(FieldDef::new("SRC_FID", FieldType::Integer));
        displacement.add_field(FieldDef::new("AREA", FieldType::Float));

        let mut lineage = Layer::new("merge_lineage");
        lineage.add_field(FieldDef::new("OUT_INDEX", FieldType::Integer));
        lineage.add_field(FieldDef::new("SRC_FID", FieldType::Integer));
        lineage.add_field(FieldDef::new("ROLE", FieldType::Text));

        let mut emitted = 0usize;
        let mut merged_pairs = 0usize;
        let mut connectors = 0usize;
        let mut displaced_area = 0.0_f64;
        let mut done = vec![false; paths.len()];

        for i in 0..paths.len() {
            if done[i] {
                continue;
            }
            match partner[i] {
                None => {
                    // Unpaired: pass through unchanged.
                    done[i] = true;
                    let fid = paths[i].0;
                    emit(
                        &mut out, &names, &layer, fid,
                        Geometry::line_string(paths[i].1.clone()),
                        "unpaired", fid as i64, -1,
                    )?;
                    add_lineage(&mut lineage, emitted, fid, "unpaired")?;
                    emitted += 1;
                }
                Some(j) => {
                    done[i] = true;
                    done[j] = true;
                    merged_pairs += 1;
                    let mid = midline(&paths[i].1, &paths[j].1);
                    let (fa, fb) = (paths[i].0, paths[j].0);
                    emit(
                        &mut out, &names, &layer, fa,
                        Geometry::line_string(mid.clone()),
                        "merged", fa as i64, fb as i64,
                    )?;
                    add_lineage(&mut lineage, emitted, fa, "merged")?;
                    add_lineage(&mut lineage, emitted, fb, "merged")?;
                    emitted += 1;

                    // Re-attach junctions so cross streets stay connected. Every
                    // vertex is checked, since a T-junction lands mid-span.
                    for k in [i, j] {
                        let verts = paths[k].1.clone();
                        for c in &verts {
                            if degree.get(&node_key(c)).copied().unwrap_or(0) < 3 {
                                continue;
                            }
                            let Some(near) = nearest_vertex(&mid, c) else {
                                continue;
                            };
                            if (near.x - c.x).hypot(near.y - c.y) <= f64::EPSILON {
                                continue;
                            }
                            emit(
                                &mut out, &names, &layer, paths[k].0,
                                Geometry::line_string(vec![c.clone(), near]),
                                "connector", paths[k].0 as i64, -1,
                            )?;
                            add_lineage(&mut lineage, emitted, paths[k].0, "connector")?;
                            emitted += 1;
                            connectors += 1;
                        }
                    }

                    // Displacement: the strip each carriageway swept to reach
                    // the midline.
                    for k in [i, j] {
                        if let Some(mp) = swept_area(&paths[k].1, &mid) {
                            let a = mp.unsigned_area();
                            if a > 0.0 {
                                displaced_area += a;
                                displacement
                                    .add_feature(
                                        Some(multipolygon_to_geometry(&mp)),
                                        &[
                                            ("SRC_FID", FieldValue::Integer(paths[k].0 as i64)),
                                            ("AREA", FieldValue::Float(a)),
                                        ],
                                    )
                                    .map_err(|e| {
                                        ToolError::Execution(format!(
                                            "failed adding displacement: {e}"
                                        ))
                                    })?;
                            }
                        }
                    }
                }
            }
        }

        ctx.progress.info(&format!(
            "merged {merged_pairs} carriageway pair(s); {connectors} junction connector(s)"
        ));
        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_path_count".to_string(), json!(paths.len()));
        outputs.insert("output_feature_count".to_string(), json!(emitted));
        outputs.insert("merged_pair_count".to_string(), json!(merged_pairs));
        outputs.insert("connector_count".to_string(), json!(connectors));
        outputs.insert("displaced_area".to_string(), json!(displaced_area));

        if matches!(args.get("output_displacement_features"), Some(v) if !v.is_null()) {
            let p = parse_optional_str(args, "output_displacement_features")?;
            outputs.insert(
                "output_displacement_features".to_string(),
                json!(write_or_store_layer(displacement, p)?),
            );
        }
        if matches!(args.get("output_table"), Some(v) if !v.is_null()) {
            let p = parse_optional_str(args, "output_table")?;
            outputs.insert("output_table".to_string(), json!(write_or_store_layer(lineage, p)?));
        }
        Ok(ToolRunResult { outputs })
    }
}

#[allow(clippy::too_many_arguments)]
fn emit(
    out: &mut Layer,
    names: &[String],
    layer: &Layer,
    fid: usize,
    geom: Geometry,
    role: &str,
    a: i64,
    b: i64,
) -> Result<(), ToolError> {
    let feat = &layer.features[fid];
    let mut attrs: Vec<(&str, FieldValue)> = names
        .iter()
        .enumerate()
        .map(|(i, nm)| {
            (
                nm.as_str(),
                feat.attributes.get(i).cloned().unwrap_or(FieldValue::Null),
            )
        })
        .collect();
    attrs.push(("MERGE_ROLE", FieldValue::Text(role.to_string())));
    attrs.push(("SRC_FID_A", FieldValue::Integer(a)));
    attrs.push(("SRC_FID_B", FieldValue::Integer(b)));
    out.add_feature(Some(geom), &attrs)
        .map_err(|e| ToolError::Execution(format!("failed adding feature: {e}")))?;
    Ok(())
}

fn add_lineage(
    lineage: &mut Layer,
    out_index: usize,
    fid: usize,
    role: &str,
) -> Result<(), ToolError> {
    lineage
        .add_feature(
            None,
            &[
                ("OUT_INDEX", FieldValue::Integer(out_index as i64)),
                ("SRC_FID", FieldValue::Integer(fid as i64)),
                ("ROLE", FieldValue::Text(role.to_string())),
            ],
        )
        .map_err(|e| ToolError::Execution(format!("failed adding lineage row: {e}")))?;
    Ok(())
}

// ── Pairing and merging ─────────────────────────────────────────────────────

/// Mean separation of two carriageways, or `None` when they are not a plausible
/// divided pair (too far apart, or not running anti-parallel).
fn carriageway_separation(a: &[Coord], b: &[Coord], max_sep: f64) -> Option<f64> {
    let (ba, bb) = (mean_bearing(a)?, mean_bearing(b)?);
    // Carriageways of a divided road run in opposite directions; allow the
    // same direction too, since digitising order is not guaranteed.
    let diff = ((ba - bb).abs() % 360.0).min(360.0 - (ba - bb).abs() % 360.0);
    let anti = (diff - 180.0).abs() <= BEARING_TOLERANCE_DEG;
    let para = diff <= BEARING_TOLERANCE_DEG;
    if !anti && !para {
        return None;
    }
    // Sample A, measure to B, and require the whole length to stay close —
    // two roads that merely cross must not pair.
    let mut total = 0.0;
    let mut worst = 0.0_f64;
    for t in 0..=SAMPLES {
        let p = point_at_fraction(a, t as f64 / SAMPLES as f64)?;
        let d = distance_to_path(b, &p);
        total += d;
        worst = worst.max(d);
    }
    if worst > max_sep {
        return None;
    }
    Some(total / (SAMPLES as f64 + 1.0))
}

/// Midline between two carriageways: sample both by fractional arc length and
/// average. Correspondence comes from the parameterisation, so no skeleton is
/// needed.
fn midline(a: &[Coord], b: &[Coord]) -> Vec<Coord> {
    // Align direction first, or the midline zig-zags end-to-end.
    let reversed: Vec<Coord>;
    let b = if needs_reverse(a, b) {
        reversed = b.iter().rev().cloned().collect();
        &reversed[..]
    } else {
        b
    };
    let mut out = Vec::with_capacity(SAMPLES + 1);
    for t in 0..=SAMPLES {
        let f = t as f64 / SAMPLES as f64;
        let (Some(pa), Some(pb)) = (point_at_fraction(a, f), point_at_fraction(b, f)) else {
            continue;
        };
        out.push(Coord::xy((pa.x + pb.x) / 2.0, (pa.y + pb.y) / 2.0));
    }
    out
}

fn needs_reverse(a: &[Coord], b: &[Coord]) -> bool {
    let (Some(a0), Some(a1)) = (a.first(), a.last()) else {
        return false;
    };
    let (Some(b0), Some(b1)) = (b.first(), b.last()) else {
        return false;
    };
    let same = (a0.x - b0.x).hypot(a0.y - b0.y) + (a1.x - b1.x).hypot(a1.y - b1.y);
    let flip = (a0.x - b1.x).hypot(a0.y - b1.y) + (a1.x - b0.x).hypot(a1.y - b0.y);
    flip < same
}

/// Polygon swept between an original carriageway and the new midline.
fn swept_area(orig: &[Coord], mid: &[Coord]) -> Option<MultiPolygon> {
    if orig.len() < 2 || mid.len() < 2 {
        return None;
    }
    let mut ring: Vec<GeoCoord> = orig.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect();
    ring.extend(mid.iter().rev().map(|c| GeoCoord { x: c.x, y: c.y }));
    if let Some(first) = ring.first().copied() {
        ring.push(first);
    }
    let poly = Polygon::new(LineString::new(ring), vec![]);
    // Self-intersecting sweep rings are common; union with itself lets
    // BooleanOps normalise them into valid polygons.
    let mp = MultiPolygon(vec![poly]);
    Some(mp.union(&MultiPolygon(Vec::new())))
}

// ── Path geometry ───────────────────────────────────────────────────────────

fn line_paths(g: &Geometry) -> Vec<Vec<Coord>> {
    match g {
        Geometry::LineString(cs) if cs.len() >= 2 => vec![cs.clone()],
        Geometry::MultiLineString(ls) => ls.iter().filter(|l| l.len() >= 2).cloned().collect(),
        Geometry::GeometryCollection(gs) => gs.iter().flat_map(line_paths).collect(),
        _ => Vec::new(),
    }
}

/// Axis-aligned bbox of a path as `[min_x, min_y, max_x, max_y]`.
fn path_bbox(p: &[Coord]) -> [f64; 4] {
    let mut bb = [f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY];
    for c in p {
        bb[0] = bb[0].min(c.x);
        bb[1] = bb[1].min(c.y);
        bb[2] = bb[2].max(c.x);
        bb[3] = bb[3].max(c.y);
    }
    bb
}

/// Gap between two bboxes; 0 when they overlap.
fn bbox_gap(a: &[f64; 4], b: &[f64; 4]) -> f64 {
    let dx = (b[0] - a[2]).max(a[0] - b[2]).max(0.0);
    let dy = (b[1] - a[3]).max(a[1] - b[3]).max(0.0);
    dx.hypot(dy)
}

fn path_length(p: &[Coord]) -> f64 {
    p.windows(2)
        .map(|w| (w[1].x - w[0].x).hypot(w[1].y - w[0].y))
        .sum()
}

fn point_at_fraction(p: &[Coord], f: f64) -> Option<Coord> {
    if p.len() < 2 {
        return p.first().cloned();
    }
    let total = path_length(p);
    if total <= 0.0 {
        return p.first().cloned();
    }
    let target = (f.clamp(0.0, 1.0)) * total;
    let mut acc = 0.0;
    for w in p.windows(2) {
        let seg = (w[1].x - w[0].x).hypot(w[1].y - w[0].y);
        if seg <= 0.0 {
            continue;
        }
        if acc + seg >= target {
            let t = (target - acc) / seg;
            return Some(Coord::xy(
                w[0].x + t * (w[1].x - w[0].x),
                w[0].y + t * (w[1].y - w[0].y),
            ));
        }
        acc += seg;
    }
    p.last().cloned()
}

fn mean_bearing(p: &[Coord]) -> Option<f64> {
    let (a, b) = (p.first()?, p.last()?);
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    if dx == 0.0 && dy == 0.0 {
        return None;
    }
    Some(dx.atan2(dy).to_degrees().rem_euclid(360.0))
}

fn distance_to_path(p: &[Coord], q: &Coord) -> f64 {
    let mut best = f64::INFINITY;
    for w in p.windows(2) {
        best = best.min(point_seg_distance(q, &w[0], &w[1]));
    }
    if p.len() == 1 {
        best = (q.x - p[0].x).hypot(q.y - p[0].y);
    }
    best
}

fn point_seg_distance(p: &Coord, a: &Coord, b: &Coord) -> f64 {
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let len2 = dx * dx + dy * dy;
    if len2 <= f64::EPSILON {
        return (p.x - a.x).hypot(p.y - a.y);
    }
    let t = (((p.x - a.x) * dx + (p.y - a.y) * dy) / len2).clamp(0.0, 1.0);
    (p.x - (a.x + t * dx)).hypot(p.y - (a.y + t * dy))
}

fn nearest_vertex(path: &[Coord], q: &Coord) -> Option<Coord> {
    path.iter()
        .min_by(|a, b| {
            (a.x - q.x)
                .hypot(a.y - q.y)
                .total_cmp(&(b.x - q.x).hypot(b.y - q.y))
        })
        .cloned()
}

fn node_key(c: &Coord) -> (i64, i64) {
    ((c.x / SNAP).round() as i64, (c.y / SNAP).round() as i64)
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

fn key_str(v: Option<&FieldValue>) -> String {
    match v {
        None | Some(FieldValue::Null) => "NULL".to_string(),
        Some(FieldValue::Integer(i)) => i.to_string(),
        Some(FieldValue::Float(f)) => format!("{f}"),
        Some(FieldValue::Text(s)) => s.clone(),
        Some(FieldValue::Boolean(b)) => b.to_string(),
        Some(FieldValue::Date(s)) | Some(FieldValue::DateTime(s)) => s.clone(),
        Some(FieldValue::Blob(b)) => format!("blob[{}]", b.len()),
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

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    struct L {
        name: &'static str,
        chr: i64,
        pts: Vec<(f64, f64)>,
    }

    fn roads(items: Vec<L>) -> String {
        let mut l = Layer::new("roads")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.add_field(FieldDef::new("chr", FieldType::Integer));
        for it in items {
            l.add_feature(
                Some(Geometry::line_string(
                    it.pts.iter().map(|(x, y)| Coord::xy(*x, *y)).collect(),
                )),
                &[
                    ("name", FieldValue::Text(it.name.to_string())),
                    ("chr", FieldValue::Integer(it.chr)),
                ],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    /// Two carriageways of "Main St" 10 apart, running opposite ways.
    fn divided_pair() -> Vec<L> {
        vec![
            L {
                name: "Main St",
                chr: 0,
                pts: vec![(0.0, 105.0), (50.0, 105.0), (100.0, 105.0)],
            },
            L {
                name: "Main St",
                chr: 0,
                pts: vec![(100.0, 95.0), (50.0, 95.0), (0.0, 95.0)],
            },
        ]
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = MergeDividedRoadsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn base(input: &str) -> serde_json::Value {
        json!({ "input": input, "merge_field": "name", "merge_distance": 20.0 })
    }

    #[test]
    fn a_divided_pair_collapses_to_one_centerline() {
        let (out, layer) = run(base(&roads(divided_pair())));
        assert_eq!(out.outputs["merged_pair_count"], json!(1));
        let role = layer.schema.field_index("MERGE_ROLE").unwrap();
        let merged: Vec<_> = layer
            .iter()
            .filter(|f| f.attributes[role].as_str() == Some("merged"))
            .collect();
        assert_eq!(merged.len(), 1);
        // The centerline sits midway between the two carriageways (y = 100).
        for c in merged[0].geometry.as_ref().unwrap().all_coords() {
            assert!((c.y - 100.0).abs() < 1e-6, "centerline y = {}", c.y);
        }
    }

    #[test]
    fn roads_with_different_names_do_not_pair() {
        let mut items = divided_pair();
        items[1].name = "Other Ave";
        let (out, _l) = run(base(&roads(items)));
        assert_eq!(out.outputs["merged_pair_count"], json!(0));
    }

    #[test]
    fn roads_beyond_the_merge_distance_do_not_pair() {
        let input = roads(divided_pair());
        let mut args = base(&input);
        args["merge_distance"] = json!(2.0); // carriageways are 10 apart
        let (out, _l) = run(args);
        assert_eq!(out.outputs["merged_pair_count"], json!(0));
    }

    #[test]
    fn crossing_roads_do_not_pair() {
        // Same name, within distance at the crossing point, but perpendicular
        // — a naive nearest-neighbour pairing would wrongly merge these.
        let (out, _l) = run(base(&roads(vec![
            L {
                name: "X",
                chr: 0,
                pts: vec![(0.0, 100.0), (100.0, 100.0)],
            },
            L {
                name: "X",
                chr: 0,
                pts: vec![(50.0, 50.0), (50.0, 150.0)],
            },
        ])));
        assert_eq!(out.outputs["merged_pair_count"], json!(0));
    }

    #[test]
    fn character_field_protects_a_feature_from_merging() {
        let mut items = divided_pair();
        items[1].chr = 1;
        let input = roads(items);
        let mut args = base(&input);
        args["character_field"] = json!("chr");
        let (out, _l) = run(args);
        assert_eq!(out.outputs["merged_pair_count"], json!(0));
    }

    #[test]
    fn junction_connectivity_is_preserved() {
        // A cross street meets the northern carriageway at (50, 105). After the
        // merge the centerline is at y = 100, so a connector must bridge the
        // gap or the network stops being routable.
        let mut items = divided_pair();
        items[0].pts = vec![(0.0, 105.0), (50.0, 105.0), (100.0, 105.0)];
        items.push(L {
            name: "Cross",
            chr: 0,
            pts: vec![(50.0, 105.0), (50.0, 160.0)],
        });
        items.push(L {
            name: "Cross2",
            chr: 0,
            pts: vec![(50.0, 105.0), (10.0, 160.0)],
        });
        let (out, layer) = run(base(&roads(items)));
        assert_eq!(out.outputs["merged_pair_count"], json!(1));
        assert!(
            out.outputs["connector_count"].as_f64().unwrap() > 0.0,
            "no junction connector emitted"
        );
        let role = layer.schema.field_index("MERGE_ROLE").unwrap();
        let conn: Vec<_> = layer
            .iter()
            .filter(|f| f.attributes[role].as_str() == Some("connector"))
            .collect();
        // The connector runs from the old junction to the new centerline.
        let cs = conn[0].geometry.as_ref().unwrap().all_coords();
        assert!((cs[0].y - 105.0).abs() < 1e-6);
        assert!((cs[cs.len() - 1].y - 100.0).abs() < 1e-6);
    }

    #[test]
    fn unpaired_roads_pass_through_unchanged() {
        let (out, layer) = run(base(&roads(vec![L {
            name: "Solo",
            chr: 0,
            pts: vec![(0.0, 0.0), (100.0, 0.0)],
        }])));
        assert_eq!(out.outputs["merged_pair_count"], json!(0));
        assert_eq!(out.outputs["output_feature_count"], json!(1));
        let role = layer.schema.field_index("MERGE_ROLE").unwrap();
        assert_eq!(layer.features[0].attributes[role].as_str(), Some("unpaired"));
        assert_eq!(
            layer.features[0].geometry.as_ref().unwrap().all_coords().len(),
            2
        );
    }

    #[test]
    fn displacement_and_lineage_are_emitted_on_request() {
        let input = roads(divided_pair());
        let mut args = base(&input);
        args["output_displacement_features"] = json!("");
        args["output_table"] = json!("");
        let a: ToolArgs = serde_json::from_value(args).unwrap();
        let out = MergeDividedRoadsTool.run(&a, &ctx()).unwrap();
        let disp =
            load_input_layer(out.outputs["output_displacement_features"].as_str().unwrap())
                .unwrap();
        assert_eq!(disp.features.len(), 2, "one displacement per carriageway");
        assert!(out.outputs["displaced_area"].as_f64().unwrap() > 0.0);
        let lin = load_input_layer(out.outputs["output_table"].as_str().unwrap()).unwrap();
        // Both source features map to the single merged output.
        assert_eq!(lin.features.len(), 2);
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            MergeDividedRoadsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "r.shp", "merge_field": "name" })).is_err());
        assert!(bad(json!({
            "input": "r.shp", "merge_field": "name", "merge_distance": 0
        }))
        .is_err());
        assert!(bad(json!({
            "input": "r.shp", "merge_field": "name", "merge_distance": 20
        }))
        .is_ok());
    }
}
