//! GeoLibre tool: delineate built-up areas from building footprints.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Delineate Built-Up Areas*
//! (Cartography): the small-scale-mapping step that generalizes a dense layer
//! of individual building footprints into a handful of settlement / urban-extent
//! polygons. It is the natural sequel to AI building-footprint extraction
//! (GeoAI) — once footprints exist, the common ask is "where are the towns?".
//!
//! The core is a morphological *closing* over the footprints:
//!
//! 1. **Dilate** every footprint outward by `grouping_distance / 2`, unioning
//!    the results. Footprints whose gap is at most `grouping_distance` overlap
//!    once dilated and so fuse into one blob.
//! 2. **Erode** the union back inward by the same radius. Bridges wider than the
//!    erosion survive (keeping fused neighbors joined) while the outer boundary
//!    returns to roughly the footprint extent. Rounded joins (the `geo` buffer
//!    default) give the smooth, generalized boundary cartography wants for free.
//! 3. **Restore** the original footprints via a union, so a closing that nibbles
//!    a tiny isolated building can never shrink a built-up area below the
//!    footprints it must contain (closing is mathematically extensive; this
//!    guards the discrete buffer approximation).
//!
//! Each resulting polygon is one candidate built-up area. Candidates are then
//! filtered: `min_building_count` drops areas that group too few buildings
//! (isolated barns, single houses) and `min_area` drops areas below a size
//! threshold and fills interior holes smaller than it. An optional
//! `simplify_tolerance` thins the vertex-dense rounded boundary with
//! Douglas–Peucker.
//!
//! A building is attributed to the area whose polygon contains its centroid, so
//! `building_count` sums back to the number of grouped footprints. The output is
//! a *new* polygon layer (`building_count`, `area` fields) — unlike in-place
//! cleanup tools, the input footprints are consumed, not carried through;
//! non-polygon features are ignored.
//!
//! Scope for v1: the ArcGIS *edge features* option (snapping built-up
//! boundaries to roads / rivers) and multi-field grouping identifiers are not
//! implemented — grouping is purely geometric.

use std::collections::BTreeMap;

use geo::{
    Area, BooleanOps, Buffer, Centroid, Contains, Coord as GeoCoord, LineString, MultiPolygon,
    Point, Polygon, Simplify, Validation,
};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct DelineateBuiltUpAreasTool;

impl Tool for DelineateBuiltUpAreasTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "delineate_built_up_areas",
            display_name: "Delineate Built-Up Areas",
            summary: "Generalize dense building footprints into smooth settlement / urban-extent polygons by grouping nearby buildings (morphological closing), with minimum building-count and area filters.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input building-footprint polygon layer, format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "grouping_distance",
                    description: "Buildings whose footprints lie within this distance (CRS units) are grouped into one built-up area. Sets the morphological-closing radius to half this value.",
                    required: true,
                },
                ToolParamSpec {
                    name: "min_building_count",
                    description: "Drop built-up areas that contain fewer than this many buildings (default 1 — keep every group).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_area",
                    description: "Drop built-up areas whose area (CRS units squared) is below this value, and fill interior holes smaller than it. Default: keep all.",
                    required: false,
                },
                ToolParamSpec {
                    name: "simplify_tolerance",
                    description: "Douglas–Peucker tolerance (CRS units) applied to the built-up boundary to thin the vertex-dense rounded buffer. Default 0 (no simplification).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args
            .get("input")
            .and_then(Value::as_str)
            .map(str::trim)
            .unwrap_or("")
            .is_empty()
        {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args
            .get("input")
            .and_then(Value::as_str)
            .filter(|s| !s.trim().is_empty())
            .ok_or_else(|| {
                ToolError::Validation("missing required parameter 'input'".to_string())
            })?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let layer_name = layer.name.clone();
        let layer_crs = layer.crs.clone();
        let input_count = layer.len();

        // Collect building footprints as `geo` polygons; remember each one's
        // centroid so a built-up area can be credited the buildings it groups.
        // Non-polygon features are ignored (this tool emits a fresh layer).
        let mut buildings: Vec<Polygon> = Vec::new();
        let mut centroids: Vec<Point> = Vec::new();
        for feature in &layer.features {
            let Some(mp) = feature.geometry.as_ref().and_then(to_multipolygon) else {
                continue;
            };
            if let Some(c) = mp.centroid() {
                centroids.push(c);
            }
            buildings.extend(mp.0);
        }
        let building_count = buildings.len();
        ctx.progress.info(&format!(
            "{input_count} feature(s): {building_count} building footprint(s) to group (radius {:.4})",
            prm.grouping_distance / 2.0
        ));

        // Morphological closing: dilate the whole set, erode back, then union the
        // originals so the result always covers every footprint (see module doc).
        let all_buildings = MultiPolygon(buildings);
        let radius = prm.grouping_distance / 2.0;
        let closed = all_buildings
            .buffer(radius)
            .buffer(-radius)
            .union(&all_buildings);

        // One candidate built-up area per closed part. Attribute buildings to the
        // part containing their centroid, then apply the count / area filters.
        let mut out_layer = Layer::new(layer_name);
        out_layer.crs = layer_crs;
        out_layer.add_field(FieldDef::new("building_count", FieldType::Integer));
        out_layer.add_field(FieldDef::new("area", FieldType::Float));
        out_layer.geom_type = Some(GeometryType::Polygon);

        let mut kept = 0usize;
        let mut dropped = 0usize;
        for part in closed.0 {
            let part = match prm.min_area {
                Some(min) => fill_small_holes(part, min),
                None => part,
            };
            let count = centroids.iter().filter(|c| part.contains(*c)).count();
            let area = part.unsigned_area();

            let too_few = count < prm.min_building_count;
            let too_small = prm.min_area.is_some_and(|min| area < min);
            if too_few || too_small {
                dropped += 1;
                continue;
            }

            let part = match prm.simplify_tolerance {
                Some(tol) if tol > 0.0 => {
                    // Douglas–Peucker can push a boundary across itself; keep the
                    // (always-valid) buffer output for any part it would invalidate.
                    let simplified = part.simplify(tol);
                    if simplified.is_valid() {
                        simplified
                    } else {
                        part
                    }
                }
                _ => part,
            };
            out_layer
                .add_feature(
                    Some(polygon_to_geometry(&part)),
                    &[
                        ("building_count", FieldValue::Integer(count as i64)),
                        ("area", FieldValue::Float(part.unsigned_area())),
                    ],
                )
                .map_err(|e| {
                    ToolError::Execution(format!("failed building output feature: {e}"))
                })?;
            kept += 1;
        }

        ctx.progress.info(&format!(
            "delineated {kept} built-up area(s), dropped {dropped} below the count/area threshold"
        ));

        let feature_count = out_layer.len();
        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_count".to_string(), json!(input_count));
        outputs.insert("building_count".to_string(), json!(building_count));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("area_count".to_string(), json!(kept));
        outputs.insert("dropped_count".to_string(), json!(dropped));
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    grouping_distance: f64,
    min_building_count: usize,
    min_area: Option<f64>,
    simplify_tolerance: Option<f64>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let grouping_distance = parse_optional_f64(args, "grouping_distance")?.ok_or_else(|| {
        ToolError::Validation("missing required parameter 'grouping_distance'".to_string())
    })?;
    if !(grouping_distance > 0.0 && grouping_distance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'grouping_distance' must be a positive number".to_string(),
        ));
    }
    let min_building_count = match parse_optional_f64(args, "min_building_count")? {
        None => 1,
        Some(v) if v.fract() == 0.0 && v >= 1.0 && v.is_finite() => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "parameter 'min_building_count' must be an integer >= 1".to_string(),
            ))
        }
    };
    let min_area = parse_optional_f64(args, "min_area")?;
    if let Some(m) = min_area {
        if !(m > 0.0 && m.is_finite()) {
            return Err(ToolError::Validation(
                "parameter 'min_area' must be a positive number".to_string(),
            ));
        }
    }
    let simplify_tolerance = parse_optional_f64(args, "simplify_tolerance")?;
    if let Some(t) = simplify_tolerance {
        if !(t >= 0.0 && t.is_finite()) {
            return Err(ToolError::Validation(
                "parameter 'simplify_tolerance' must be a non-negative number".to_string(),
            ));
        }
    }
    Ok(Params {
        grouping_distance,
        min_building_count,
        min_area,
        simplify_tolerance,
    })
}

/// Parses an optional numeric parameter, accepting a JSON number or a numeric
/// string (host UIs often post form values as strings).
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

// ── Geometry helpers ──────────────────────────────────────────────────────────

/// Drops interior rings (holes) whose area is below `min_area`, e.g. courtyards
/// or gaps between grouped buildings too small to keep at the output scale.
fn fill_small_holes(poly: Polygon, min_area: f64) -> Polygon {
    let (exterior, interiors) = poly.into_inner();
    let kept: Vec<LineString> = interiors
        .into_iter()
        .filter(|ring| Polygon::new(ring.clone(), vec![]).unsigned_area() >= min_area)
        .collect();
    Polygon::new(exterior, kept)
}

/// Converts a polygonal `wbvector` geometry to a `geo` `MultiPolygon`. Returns
/// `None` for non-polygon geometries (which the tool ignores).
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
    // `geo` closes rings itself; the missing closing vertex in `Ring` is fine.
    LineString::new(
        ring.coords()
            .iter()
            .map(|c| GeoCoord { x: c.x, y: c.y })
            .collect(),
    )
}

/// Converts a single `geo` `Polygon` to a `wbvector` `Geometry::Polygon`.
fn polygon_to_geometry(poly: &Polygon) -> Geometry {
    Geometry::Polygon {
        exterior: linestring_to_ring(poly.exterior()),
        interiors: poly.interiors().iter().map(linestring_to_ring).collect(),
    }
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
    use wbvector::{memory_store, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Axis-aligned square footprint with lower-left corner (x0, y0) and side `s`.
    fn square(x0: f64, y0: f64, s: f64) -> Geometry {
        Geometry::polygon(
            vec![
                Coord::xy(x0, y0),
                Coord::xy(x0 + s, y0),
                Coord::xy(x0 + s, y0 + s),
                Coord::xy(x0, y0 + s),
            ],
            vec![],
        )
    }

    fn layer_of(geoms: Vec<Geometry>) -> String {
        let mut layer = Layer::new("buildings");
        for g in geoms {
            layer.add_feature(Some(g), &[]).unwrap();
        }
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = DelineateBuiltUpAreasTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn count_of(layer: &Layer, idx: usize) -> i64 {
        match layer.features[idx]
            .get(&layer.schema, "building_count")
            .unwrap()
        {
            FieldValue::Integer(i) => *i,
            other => panic!("building_count should be an integer, got {other:?}"),
        }
    }

    /// Two footprints a short gap apart fuse into one area; a distant third
    /// stands alone. Grouping distance (6) exceeds the near gap (2) but not the
    /// far gap (~18).
    #[test]
    fn groups_nearby_buildings_and_isolates_distant_ones() {
        let input = layer_of(vec![
            square(0.0, 0.0, 4.0),  // A: x 0..4
            square(6.0, 0.0, 4.0),  // B: x 6..10  (gap to A = 2)
            square(28.0, 0.0, 4.0), // C: far away
        ]);
        let (out, layer) = run_tool(json!({ "input": input, "grouping_distance": 6.0 }));

        assert_eq!(out.outputs["building_count"], json!(3));
        assert_eq!(out.outputs["feature_count"], json!(2));
        // The two grouped buildings are credited to one area; the loner to the
        // other. Counts sum back to the 3 input footprints.
        let mut counts: Vec<i64> = (0..layer.len()).map(|i| count_of(&layer, i)).collect();
        counts.sort_unstable();
        assert_eq!(counts, vec![1, 2]);
    }

    /// A chain of footprints each within the grouping distance of the next fuses
    /// into a single built-up area crediting all three.
    #[test]
    fn chain_of_close_buildings_fuses() {
        let input = layer_of(vec![
            square(0.0, 0.0, 4.0),  // x 0..4
            square(6.0, 0.0, 4.0),  // x 6..10 (gap 2)
            square(12.0, 0.0, 4.0), // x 12..16 (gap 2)
        ]);
        let (out, layer) = run_tool(json!({ "input": input, "grouping_distance": 6.0 }));
        assert_eq!(out.outputs["feature_count"], json!(1));
        assert_eq!(count_of(&layer, 0), 3);
    }

    /// The output always covers the footprints it groups: a built-up area's
    /// polygon is a superset, so its area exceeds the summed footprint area.
    #[test]
    fn built_up_area_covers_its_buildings() {
        let input = layer_of(vec![square(0.0, 0.0, 4.0), square(6.0, 0.0, 4.0)]);
        let (_, layer) = run_tool(json!({ "input": input, "grouping_distance": 6.0 }));
        assert_eq!(layer.len(), 1);
        let area = match layer.features[0].get(&layer.schema, "area").unwrap() {
            FieldValue::Float(a) => *a,
            other => panic!("area should be a float, got {other:?}"),
        };
        // Two 4x4 footprints = 32; the closed area must strictly exceed that.
        assert!(
            area > 32.0,
            "built-up area {area} does not cover its 32 units of footprint"
        );
    }

    /// `min_building_count` drops areas that group too few buildings.
    #[test]
    fn min_building_count_drops_sparse_areas() {
        let input = layer_of(vec![
            square(0.0, 0.0, 4.0),
            square(6.0, 0.0, 4.0),
            square(28.0, 0.0, 4.0), // isolated -> a 1-building area, must be dropped
        ]);
        let (out, layer) =
            run_tool(json!({ "input": input, "grouping_distance": 6.0, "min_building_count": 2 }));
        assert_eq!(out.outputs["feature_count"], json!(1));
        assert_eq!(out.outputs["dropped_count"], json!(1));
        assert_eq!(count_of(&layer, 0), 2);
    }

    /// `min_area` drops built-up areas below the size threshold.
    #[test]
    fn min_area_drops_small_areas() {
        // One tiny loner far from a fused pair. With a threshold above the
        // loner's area but below the pair's, only the pair survives.
        let input = layer_of(vec![
            square(0.0, 0.0, 4.0),
            square(6.0, 0.0, 4.0),
            square(40.0, 0.0, 1.0), // tiny isolated footprint
        ]);
        let (out, _) =
            run_tool(json!({ "input": input, "grouping_distance": 6.0, "min_area": 20.0 }));
        assert_eq!(out.outputs["feature_count"], json!(1));
        assert_eq!(out.outputs["dropped_count"], json!(1));
    }

    /// Non-polygon features are ignored, not carried into the output.
    #[test]
    fn ignores_non_polygon_features() {
        let mut layer = Layer::new("mixed");
        layer
            .add_feature(Some(Geometry::point(1.0, 1.0)), &[])
            .unwrap();
        layer.add_feature(Some(square(0.0, 0.0, 4.0)), &[]).unwrap();
        layer.add_feature(Some(square(6.0, 0.0, 4.0)), &[]).unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input, "grouping_distance": 6.0 }));
        assert_eq!(out.outputs["building_count"], json!(2));
        assert_eq!(out.outputs["feature_count"], json!(1));
        // Output holds only polygons.
        assert!(layer
            .features
            .iter()
            .all(|f| matches!(f.geometry, Some(Geometry::Polygon { .. }))));
    }

    /// An empty (no-polygon) input yields an empty output, not an error.
    #[test]
    fn empty_input_yields_no_areas() {
        let input = layer_of(vec![Geometry::point(0.0, 0.0)]);
        let (out, _) = run_tool(json!({ "input": input, "grouping_distance": 5.0 }));
        assert_eq!(out.outputs["building_count"], json!(0));
        assert_eq!(out.outputs["feature_count"], json!(0));
    }

    /// Simplifying the dense rounded boundary must never produce a
    /// self-intersecting (invalid) polygon, even with interior holes.
    #[test]
    fn simplify_tolerance_keeps_output_valid() {
        use geo::Validation;
        // A ring of buildings closes into a built-up area with a central hole.
        let mut geoms = Vec::new();
        for i in 0..12 {
            let a = i as f64 * std::f64::consts::TAU / 12.0;
            geoms.push(square(30.0 * a.cos() - 2.0, 30.0 * a.sin() - 2.0, 4.0));
        }
        let input = layer_of(geoms);
        let (out, layer) = run_tool(
            json!({ "input": input, "grouping_distance": 20.0, "simplify_tolerance": 1.0 }),
        );
        assert!(out.outputs["feature_count"].as_i64().unwrap() >= 1);
        for f in &layer.features {
            let mp = to_multipolygon(f.geometry.as_ref().unwrap()).unwrap();
            assert!(mp.is_valid(), "simplified built-up area must be valid");
        }
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = DelineateBuiltUpAreasTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input must fail");
        assert!(
            bad(json!({ "input": "x.geojson" })).is_err(),
            "missing grouping_distance must fail"
        );
        assert!(bad(json!({ "input": "x.geojson", "grouping_distance": 0 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "grouping_distance": -5 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "grouping_distance": "wide" })).is_err());
        assert!(bad(
            json!({ "input": "x.geojson", "grouping_distance": 5, "min_building_count": 0 })
        )
        .is_err());
        assert!(bad(
            json!({ "input": "x.geojson", "grouping_distance": 5, "min_building_count": 2.5 })
        )
        .is_err());
        assert!(
            bad(json!({ "input": "x.geojson", "grouping_distance": 5, "min_area": -1 })).is_err()
        );
        assert!(
            bad(json!({ "input": "x.geojson", "grouping_distance": "5", "min_area": "20" }))
                .is_ok(),
            "numeric strings ok"
        );
    }
}
