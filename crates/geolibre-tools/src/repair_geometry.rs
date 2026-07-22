//! GeoLibre tool: detect and fix invalid vector geometry.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Repair Geometry* (and its read-only
//! sibling *Check Geometry*): the standard cleanup step after importing a
//! real-world shapefile or converting a raster to polygons, where the data is
//! peppered with self-intersecting rings, mis-wound rings, duplicate vertices
//! and empty parts that break every downstream overlay/area tool.
//!
//! The ~791 bundled whitebox tools ship only `clean_vector`, which *drops*
//! null/degenerate geometries — it never repairs a self-intersecting polygon,
//! fixes ring orientation, or removes duplicate vertices. This tool does.
//!
//! For each polygonal feature it:
//!
//! - removes consecutive duplicate vertices (and the redundant closing vertex);
//! - drops degenerate rings (fewer than three distinct points) and null/empty
//!   parts, dropping the whole feature when nothing survives;
//! - resolves self-intersections by taking the polygon's **self-union**
//!   (`geo`'s `unary_union`, pure Rust, no GEOS): unioning a bow-tie with
//!   itself splits it into the two clean lobes it implies;
//! - re-winds every ring to the OGC convention (exterior counter-clockwise,
//!   holes clockwise) with `geo`'s `Orient`.
//!
//! Validity of a ring is judged with `geo`'s OGC `Validation` after the cheap
//! vertex cleanup, so a polygon that is merely mis-wound or carries duplicate
//! points is not mislabelled as self-intersecting.
//!
//! Non-polygon geometries (points, lines) pass through untouched.
//!
//! With `check_only=true` nothing is repaired: the tool instead emits the input
//! layer verbatim with two added columns — `problem_code` (e.g.
//! `self_intersection`, `ring_orientation`, `duplicate_vertex`, `null_empty`,
//! or `ok`) and `problem_desc` — mirroring ArcGIS *Check Geometry*'s report
//! table. Per-problem counts are always returned in the tool outputs.

use std::collections::BTreeMap;

use geo::orient::Direction;
use geo::{unary_union, Coord as GeoCoord, LineString, MultiPolygon, Orient, Polygon, Validation};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct RepairGeometryTool;

// Problem codes, ordered by reporting priority.
const CODE_NULL_EMPTY: &str = "null_empty";
const CODE_SELF_INTERSECTION: &str = "self_intersection";
const CODE_RING_ORIENTATION: &str = "ring_orientation";
const CODE_DUPLICATE_VERTEX: &str = "duplicate_vertex";

impl Tool for RepairGeometryTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "repair_geometry",
            display_name: "Repair Geometry",
            summary: "Detect and fix invalid polygon geometry — self-intersections (via self-union), ring winding, duplicate vertices, and null/empty parts — or, with check_only, report a per-feature problem code. Like ArcGIS Repair Geometry / Check Geometry.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon vector file path, format auto-detected (or an in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "check_only",
                    description: "When true, do not modify geometry: emit the input with 'problem_code'/'problem_desc' columns (like ArcGIS Check Geometry). Default false.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_bool(args, "check_only")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let check_only = parse_bool(args, "check_only")?.unwrap_or(false);

        let layer = load_input_layer(input)?;
        let input_count = layer.len();
        let orig_field_count = layer.schema.len();

        // Per-problem tallies (a feature can contribute to several).
        let mut null_empty = 0usize;
        let mut self_intersection = 0usize;
        let mut ring_orientation = 0usize;
        let mut duplicate_vertex = 0usize;
        let mut invalid_features = 0usize;

        if check_only {
            // Report mode: copy the layer verbatim + two diagnostic columns.
            let mut out = Layer::new(layer.name.clone());
            out.geom_type = layer.geom_type;
            out.crs = layer.crs.clone();
            out.schema = layer.schema.clone();
            out.add_field(FieldDef::new("problem_code", FieldType::Text));
            out.add_field(FieldDef::new("problem_desc", FieldType::Text));

            for feature in layer.iter() {
                let analysis = analyze(feature.geometry.as_ref());
                tally(
                    &analysis.codes,
                    &mut null_empty,
                    &mut self_intersection,
                    &mut ring_orientation,
                    &mut duplicate_vertex,
                    &mut invalid_features,
                );
                let mut f = feature.clone();
                f.attributes.resize(orig_field_count, FieldValue::Null);
                f.attributes
                    .push(FieldValue::Text(code_label(&analysis.codes)));
                f.attributes
                    .push(FieldValue::Text(describe(&analysis.codes)));
                out.push(f);
            }

            ctx.progress.info(&format!(
                "checked {input_count} feature(s): {invalid_features} with a geometry problem"
            ));

            let feature_count = out.len();
            let out_path = write_or_store_layer(out, output)?;
            return Ok(counts_result(
                out_path,
                true,
                input_count,
                feature_count,
                0,
                invalid_features,
                null_empty,
                self_intersection,
                ring_orientation,
                duplicate_vertex,
            ));
        }

        // Repair mode: rewrite each geometry, dropping features that become null.
        let mut out = Layer::new(layer.name.clone());
        out.geom_type = layer.geom_type;
        out.crs = layer.crs.clone();
        out.schema = layer.schema.clone();

        let mut dropped = 0usize;
        let mut has_multipolygon = false;
        for feature in layer.features.into_iter() {
            let analysis = analyze(feature.geometry.as_ref());
            tally(
                &analysis.codes,
                &mut null_empty,
                &mut self_intersection,
                &mut ring_orientation,
                &mut duplicate_vertex,
                &mut invalid_features,
            );
            let mut feature = feature;
            match analysis.repaired {
                Some(geom) => {
                    has_multipolygon |= matches!(geom, Geometry::MultiPolygon(_));
                    feature.geometry = Some(geom);
                    feature.fid = out.len() as u64;
                    out.push(feature);
                }
                None => {
                    // Null/empty geometry that could not be repaired: drop it.
                    dropped += 1;
                }
            }
        }
        // A self-union can split a single-part polygon into several parts.
        if has_multipolygon {
            out.geom_type = Some(wbvector::GeometryType::MultiPolygon);
        }

        ctx.progress.info(&format!(
            "repaired {invalid_features} of {input_count} feature(s); {dropped} null/empty dropped"
        ));

        let feature_count = out.len();
        let out_path = write_or_store_layer(out, output)?;
        Ok(counts_result(
            out_path,
            false,
            input_count,
            feature_count,
            dropped,
            invalid_features,
            null_empty,
            self_intersection,
            ring_orientation,
            duplicate_vertex,
        ))
    }
}

// ── Analysis / repair ───────────────────────────────────────────────────────

struct Analysis {
    /// Distinct problem codes found (empty == clean).
    codes: Vec<&'static str>,
    /// Repaired geometry, or `None` when the feature is null/empty.
    repaired: Option<Geometry>,
}

/// Detects problems and (for polygons) repairs the geometry.
fn analyze(geom: Option<&Geometry>) -> Analysis {
    let Some(geom) = geom else {
        return Analysis {
            codes: vec![CODE_NULL_EMPTY],
            repaired: None,
        };
    };
    if geom.is_empty() {
        return Analysis {
            codes: vec![CODE_NULL_EMPTY],
            repaired: None,
        };
    }

    // Split into polygon parts; non-polygons pass through untouched.
    let parts = match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => vec![(exterior.clone(), interiors.clone())],
        Geometry::MultiPolygon(parts) => parts.clone(),
        _ => {
            return Analysis {
                codes: Vec::new(),
                repaired: Some(geom.clone()),
            }
        }
    };

    let mut codes: Vec<&'static str> = Vec::new();
    let mut has_duplicate = false;
    let mut has_orientation = false;

    // Clean every ring (dedupe vertices, check winding, drop degenerate rings).
    let mut clean_parts: Vec<(Vec<Coord>, Vec<Vec<Coord>>)> = Vec::new();
    for (exterior, interiors) in &parts {
        let (ext, ext_dup) = clean_ring(exterior);
        has_duplicate |= ext_dup;
        if distinct_count(&ext) < 3 {
            // Exterior collapsed: the whole part is degenerate/empty.
            continue;
        }
        if ring_signed_area(&ext) <= 0.0 {
            has_orientation = true; // exterior should be CCW (positive area)
        }
        let mut holes: Vec<Vec<Coord>> = Vec::new();
        for hole in interiors {
            let (h, h_dup) = clean_ring(hole);
            has_duplicate |= h_dup;
            if distinct_count(&h) < 3 {
                continue; // drop degenerate hole
            }
            if ring_signed_area(&h) >= 0.0 {
                has_orientation = true; // holes should be CW (negative area)
            }
            holes.push(h);
        }
        clean_parts.push((ext, holes));
    }

    if clean_parts.is_empty() {
        return Analysis {
            codes: vec![CODE_NULL_EMPTY],
            repaired: None,
        };
    }

    // OGC validity on the cleaned geometry: catches self-intersections and
    // ring self-touches (but not mere winding, which OGC does not require).
    let cleaned_mp = clean_parts_to_multipolygon(&clean_parts);
    let self_intersects = !cleaned_mp.is_valid();
    if self_intersects {
        codes.push(CODE_SELF_INTERSECTION);
    }
    // Winding is only meaningful for a topologically valid ring; a
    // self-intersecting ring's winding is ambiguous and is fixed by the
    // self-union + orient below regardless.
    if has_orientation && !self_intersects {
        codes.push(CODE_RING_ORIENTATION);
    }
    if has_duplicate {
        codes.push(CODE_DUPLICATE_VERTEX);
    }

    if codes.is_empty() {
        // Already valid, correctly wound, no duplicates: leave it untouched so
        // valid features are byte-for-byte preserved.
        return Analysis {
            codes,
            repaired: Some(geom.clone()),
        };
    }

    // Repair: self-union resolves self-intersections, then re-orient to OGC.
    let repaired_mp = unary_union(std::iter::once(&cleaned_mp)).orient(Direction::Default);
    let repaired = if repaired_mp.0.is_empty() {
        None
    } else {
        Some(multipolygon_to_geometry(&repaired_mp))
    };
    if repaired.is_none() {
        return Analysis {
            codes: vec![CODE_NULL_EMPTY],
            repaired: None,
        };
    }
    Analysis { codes, repaired }
}

/// Removes consecutive duplicate vertices (and the redundant closing vertex).
/// Returns the cleaned coordinates and whether any interior duplicate was
/// removed (the implicit closing duplicate is stripped silently — it is normal,
/// not a defect).
fn clean_ring(ring: &Ring) -> (Vec<Coord>, bool) {
    let src = ring.coords();
    let mut out: Vec<Coord> = Vec::with_capacity(src.len());
    let mut dropped = false;
    for c in src {
        if let Some(last) = out.last() {
            if coords_eq(last, c) {
                dropped = true;
                continue;
            }
        }
        out.push(c.clone());
    }
    // Strip an explicit closing vertex (first == last) without flagging it.
    while out.len() >= 2 && coords_eq(&out[0], out.last().unwrap()) {
        out.pop();
    }
    (out, dropped)
}

fn coords_eq(a: &Coord, b: &Coord) -> bool {
    a.x == b.x && a.y == b.y
}

/// Count of distinct vertices (a ring needs at least three).
fn distinct_count(coords: &[Coord]) -> usize {
    // Coordinates are already de-duplicated consecutively; count unique points.
    let mut seen: Vec<Coord> = Vec::new();
    for c in coords {
        if !seen.iter().any(|s| coords_eq(s, c)) {
            seen.push(c.clone());
        }
    }
    seen.len()
}

/// Shoelace signed area (positive == counter-clockwise).
fn ring_signed_area(coords: &[Coord]) -> f64 {
    let n = coords.len();
    if n < 3 {
        return 0.0;
    }
    let mut a = 0.0f64;
    for i in 0..n {
        let j = (i + 1) % n;
        a += coords[i].x * coords[j].y - coords[j].x * coords[i].y;
    }
    a * 0.5
}

// ── Outputs ─────────────────────────────────────────────────────────────────

#[allow(clippy::too_many_arguments)]
fn counts_result(
    out_path: String,
    check_only: bool,
    input_count: usize,
    feature_count: usize,
    dropped: usize,
    invalid_features: usize,
    null_empty: usize,
    self_intersection: usize,
    ring_orientation: usize,
    duplicate_vertex: usize,
) -> ToolRunResult {
    let mut outputs = BTreeMap::new();
    outputs.insert("output".to_string(), json!(out_path));
    outputs.insert("check_only".to_string(), json!(check_only));
    outputs.insert("input_count".to_string(), json!(input_count));
    outputs.insert("feature_count".to_string(), json!(feature_count));
    outputs.insert("dropped_count".to_string(), json!(dropped));
    outputs.insert("invalid_count".to_string(), json!(invalid_features));
    outputs.insert("null_empty_count".to_string(), json!(null_empty));
    outputs.insert(
        "self_intersection_count".to_string(),
        json!(self_intersection),
    );
    outputs.insert(
        "ring_orientation_count".to_string(),
        json!(ring_orientation),
    );
    outputs.insert(
        "duplicate_vertex_count".to_string(),
        json!(duplicate_vertex),
    );
    ToolRunResult { outputs }
}

#[allow(clippy::too_many_arguments)]
fn tally(
    codes: &[&str],
    null_empty: &mut usize,
    self_intersection: &mut usize,
    ring_orientation: &mut usize,
    duplicate_vertex: &mut usize,
    invalid_features: &mut usize,
) {
    if !codes.is_empty() {
        *invalid_features += 1;
    }
    for &c in codes {
        match c {
            CODE_NULL_EMPTY => *null_empty += 1,
            CODE_SELF_INTERSECTION => *self_intersection += 1,
            CODE_RING_ORIENTATION => *ring_orientation += 1,
            CODE_DUPLICATE_VERTEX => *duplicate_vertex += 1,
            _ => {}
        }
    }
}

fn code_label(codes: &[&str]) -> String {
    if codes.is_empty() {
        "ok".to_string()
    } else {
        codes.join("+")
    }
}

fn describe(codes: &[&str]) -> String {
    if codes.is_empty() {
        return "ok".to_string();
    }
    codes
        .iter()
        .map(|c| match *c {
            CODE_NULL_EMPTY => "null or empty geometry",
            CODE_SELF_INTERSECTION => "self-intersecting ring(s)",
            CODE_RING_ORIENTATION => "incorrect ring winding (exterior not CCW / hole not CW)",
            CODE_DUPLICATE_VERTEX => "duplicate/consecutive vertices",
            other => other,
        })
        .collect::<Vec<_>>()
        .join("; ")
}

// ── geo <-> wbvector conversion ─────────────────────────────────────────────

fn clean_parts_to_multipolygon(parts: &[(Vec<Coord>, Vec<Vec<Coord>>)]) -> MultiPolygon {
    MultiPolygon(
        parts
            .iter()
            .map(|(ext, holes)| {
                Polygon::new(
                    coords_to_linestring(ext),
                    holes.iter().map(|h| coords_to_linestring(h)).collect(),
                )
            })
            .collect(),
    )
}

fn coords_to_linestring(coords: &[Coord]) -> LineString {
    // `geo` closes rings itself; leaving off the closing vertex is fine.
    LineString::new(coords.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect())
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
    // Drop the closing duplicate `geo` keeps; `Ring` stores it implicitly.
    let mut coords: Vec<Coord> = ls.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
    if coords.len() >= 2 && coords.first() == coords.last() {
        coords.pop();
    }
    Ring::new(coords)
}

// ── Params ──────────────────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

/// Parses an optional boolean, accepting a JSON bool or the strings
/// "true"/"false"/"1"/"0" (host UIs often post booleans as strings).
fn parse_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(Value::Number(n)) => match n.as_i64() {
            Some(0) => Ok(Some(false)),
            Some(1) => Ok(Some(true)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use geo::Area;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = RepairGeometryTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn geo_area(geom: &Geometry) -> f64 {
        match geom {
            Geometry::Polygon {
                exterior,
                interiors,
            } => clean_parts_to_multipolygon(&[(
                exterior.coords().to_vec(),
                interiors.iter().map(|r| r.coords().to_vec()).collect(),
            )])
            .unsigned_area(),
            Geometry::MultiPolygon(parts) => clean_parts_to_multipolygon(
                &parts
                    .iter()
                    .map(|(e, hs)| {
                        (
                            e.coords().to_vec(),
                            hs.iter().map(|r| r.coords().to_vec()).collect(),
                        )
                    })
                    .collect::<Vec<_>>(),
            )
            .unsigned_area(),
            _ => 0.0,
        }
    }

    /// A bow-tie (self-intersecting) polygon is invalid; repair splits it into
    /// the two clean triangular lobes and the result is OGC-valid.
    #[test]
    fn resolves_self_intersection() {
        let mut layer = Layer::new("bowtie").with_geom_type(GeometryType::Polygon);
        // Bow-tie: (0,0)->(2,2)->(2,0)->(0,2)->close. Self-intersects at (1,1).
        layer
            .add_feature(
                Some(Geometry::polygon(
                    vec![
                        Coord::xy(0.0, 0.0),
                        Coord::xy(2.0, 2.0),
                        Coord::xy(2.0, 0.0),
                        Coord::xy(0.0, 2.0),
                    ],
                    vec![],
                )),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input }));
        assert_eq!(out.outputs["self_intersection_count"], json!(1));
        assert_eq!(out.outputs["invalid_count"], json!(1));
        assert_eq!(out.outputs["feature_count"], json!(1));
        // The self-intersection is gone: the repaired geometry is OGC-valid.
        // (geo's unary_union uses a positive/negative fill rule, so a symmetric
        // bow-tie collapses to its dominant lobe — a valid triangle of area 1.)
        let g = layer.features[0].geometry.as_ref().unwrap();
        let (ext, holes): (Vec<Coord>, Vec<Vec<Coord>>) = match g {
            Geometry::Polygon {
                exterior,
                interiors,
            } => (
                exterior.coords().to_vec(),
                interiors.iter().map(|r| r.coords().to_vec()).collect(),
            ),
            _ => panic!("expected a single polygon"),
        };
        let mp = clean_parts_to_multipolygon(&[(ext, holes)]);
        assert!(mp.is_valid(), "repaired geometry must be OGC-valid");
        assert!((geo_area(g) - 1.0).abs() < 1e-9);
    }

    /// Duplicate consecutive vertices are removed; area is preserved exactly.
    #[test]
    fn removes_duplicate_vertices_preserving_area() {
        let mut layer = Layer::new("dups").with_geom_type(GeometryType::Polygon);
        layer
            .add_feature(
                Some(Geometry::polygon(
                    vec![
                        Coord::xy(0.0, 0.0),
                        Coord::xy(0.0, 0.0), // duplicate
                        Coord::xy(4.0, 0.0),
                        Coord::xy(4.0, 3.0),
                        Coord::xy(4.0, 3.0), // duplicate
                        Coord::xy(0.0, 3.0),
                    ],
                    vec![],
                )),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input }));
        assert_eq!(out.outputs["duplicate_vertex_count"], json!(1));
        let g = layer.features[0].geometry.as_ref().unwrap();
        // 4x3 rectangle => area 12, unchanged by dedup.
        assert!((geo_area(g) - 12.0).abs() < 1e-9);
        if let Geometry::Polygon { exterior, .. } = g {
            assert_eq!(exterior.len(), 4, "duplicates removed, 4 distinct vertices");
        } else {
            panic!("expected a single polygon");
        }
    }

    /// A clockwise exterior ring is flagged and re-wound to CCW; area preserved.
    #[test]
    fn fixes_ring_orientation() {
        let mut layer = Layer::new("cw").with_geom_type(GeometryType::Polygon);
        // Clockwise square (negative signed area).
        layer
            .add_feature(
                Some(Geometry::polygon(
                    vec![
                        Coord::xy(0.0, 0.0),
                        Coord::xy(0.0, 2.0),
                        Coord::xy(2.0, 2.0),
                        Coord::xy(2.0, 0.0),
                    ],
                    vec![],
                )),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input }));
        assert_eq!(out.outputs["ring_orientation_count"], json!(1));
        let g = layer.features[0].geometry.as_ref().unwrap();
        if let Geometry::Polygon { exterior, .. } = g {
            assert!(
                ring_signed_area(exterior.coords()) > 0.0,
                "exterior should be CCW after repair"
            );
        } else {
            panic!("expected a single polygon");
        }
        assert!((geo_area(g) - 4.0).abs() < 1e-9);
    }

    /// A valid, correctly-wound polygon passes through unchanged and is not
    /// counted as invalid; points pass through untouched too.
    #[test]
    fn valid_and_non_polygon_pass_through() {
        let mut layer = Layer::new("mixed");
        // Valid CCW square.
        layer
            .add_feature(
                Some(Geometry::polygon(
                    vec![
                        Coord::xy(0.0, 0.0),
                        Coord::xy(2.0, 0.0),
                        Coord::xy(2.0, 2.0),
                        Coord::xy(0.0, 2.0),
                    ],
                    vec![],
                )),
                &[],
            )
            .unwrap();
        layer
            .add_feature(Some(Geometry::point(5.0, 5.0)), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input }));
        assert_eq!(out.outputs["invalid_count"], json!(0));
        assert_eq!(out.outputs["feature_count"], json!(2));
        assert!(layer
            .features
            .iter()
            .any(|f| matches!(f.geometry, Some(Geometry::Point(_)))));
    }

    /// check_only leaves geometry alone and adds problem columns.
    #[test]
    fn check_only_reports_without_repair() {
        let mut layer = Layer::new("bowtie").with_geom_type(GeometryType::Polygon);
        layer
            .add_feature(
                Some(Geometry::polygon(
                    vec![
                        Coord::xy(0.0, 0.0),
                        Coord::xy(2.0, 2.0),
                        Coord::xy(2.0, 0.0),
                        Coord::xy(0.0, 2.0),
                    ],
                    vec![],
                )),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input, "check_only": true }));
        assert_eq!(out.outputs["invalid_count"], json!(1));
        assert_eq!(out.outputs["self_intersection_count"], json!(1));
        // Geometry unchanged: still the original 4-vertex bow-tie polygon.
        let g = layer.features[0].geometry.as_ref().unwrap();
        assert!(matches!(g, Geometry::Polygon { .. }));
        // Diagnostic columns present and populated.
        let code = layer.features[0]
            .get(&layer.schema, "problem_code")
            .unwrap();
        assert_eq!(code, &FieldValue::Text("self_intersection".into()));
    }

    #[test]
    fn drops_null_geometry() {
        let mut layer = Layer::new("nulls").with_geom_type(GeometryType::Polygon);
        layer.add_feature(None, &[]).unwrap();
        layer
            .add_feature(
                Some(Geometry::polygon(
                    vec![
                        Coord::xy(0.0, 0.0),
                        Coord::xy(2.0, 0.0),
                        Coord::xy(2.0, 2.0),
                        Coord::xy(0.0, 2.0),
                    ],
                    vec![],
                )),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input }));
        assert_eq!(out.outputs["null_empty_count"], json!(1));
        assert_eq!(out.outputs["dropped_count"], json!(1));
        assert_eq!(out.outputs["feature_count"], json!(1));
        assert_eq!(layer.len(), 1);
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = RepairGeometryTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input must fail");
        assert!(
            bad(json!({ "input": "x.geojson", "check_only": "maybe" })).is_err(),
            "bad boolean must fail"
        );
        assert!(bad(json!({ "input": "x.geojson" })).is_ok());
        assert!(
            bad(json!({ "input": "x.geojson", "check_only": "true" })).is_ok(),
            "boolean strings ok"
        );
    }
}
