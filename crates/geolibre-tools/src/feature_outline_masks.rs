//! GeoLibre tool: cartographic masking polygons around features.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Feature Outline Masks* (Cartography →
//! Masking Tools): generate a polygon at a specified `margin` around every input
//! feature, used to knock out ("halo") other layers around labels, annotation,
//! or symbols so they stay legible. The bundled whitebox-wasm suite has nothing
//! comparable (`unsharp_masking` is an unrelated image filter). Mask polygons
//! pair naturally with the repo's `render_vector_png` and web-map outputs.
//!
//! Three mask shapes are offered, mirroring ArcGIS:
//!
//! - `exact` (default) — the feature outline dilated by the margin: the input
//!   geometry buffered outward by `margin` (rounded joins, `geo`'s `Buffer`).
//!   For points this is a disk, for lines a capsule, for polygons the grown
//!   footprint. This is the tightest mask and follows the feature exactly.
//! - `convex_hull` — the convex hull of the feature's vertices, then dilated by
//!   the margin. Cheaper and smoother than `exact` for busy geometry.
//! - `box` — the feature's bounding envelope expanded by `margin` on every side
//!   (a clean rectangle). The coarsest, fastest mask.
//!
//! Every mask fully contains its source feature (a dilation of a superset of the
//! feature), so `source ⊆ mask` always holds. Each output polygon carries the
//! 0-based index of its source feature in `source_fid`, the `mask_kind`, and —
//! when `id_field` is given — the source id value in `source_id`.
//!
//! **Intersecting Layers Masks.** When a second `masked_layer` is supplied, the
//! output is clipped to only the parts of each mask that overlap the polygonal
//! features of that layer (via `geo`'s `BooleanOps::intersection`), matching the
//! ArcGIS *Intersecting Layers Masks* behaviour: masks are emitted only where
//! the masked layer actually has symbology to knock out. Masks that miss the
//! masked layer entirely are dropped.
//!
//! **Units.** `margin` is read in metres for a geographic layer (EPSG:4326 or an
//! unknown CRS) and converted to degrees at ~111320 m/°; for a projected layer
//! it is taken directly in the layer's CRS units. Buffering happens in the
//! layer's coordinate space, so on geographic data the longitude dilation is
//! slightly compressed away from the equator — containment is unaffected.
//!
//! Scope for v1: source attributes other than the chosen `id_field` are not
//! copied onto the masks, and the geographic metre→degree conversion uses a
//! single spherical scale factor (no per-latitude longitude correction).

use std::collections::BTreeMap;

use geo::{
    BooleanOps, BoundingRect, Buffer, ConvexHull, Coord as GeoCoord, Geometry as GeoGeometry,
    HasDimensions, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Metres per degree of latitude (spherical approximation) used to interpret a
/// metric `margin` on a geographic layer.
const METERS_PER_DEGREE: f64 = 111_320.0;

pub struct FeatureOutlineMasksTool;

impl Tool for FeatureOutlineMasksTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "feature_outline_masks",
            display_name: "Feature Outline Masks",
            summary: "Generate cartographic mask polygons at a margin around features (exact outline, convex hull, or bounding box), optionally clipped to a second layer (Intersecting Layers Masks).",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (any geometry), format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "margin",
                    description: "Mask margin distance. Metres for a geographic CRS (EPSG:4326 or unknown, converted to degrees), otherwise in the layer's CRS units. Must be positive.",
                    required: true,
                },
                ToolParamSpec {
                    name: "mask_kind",
                    description: "Mask shape: 'exact' (default, buffered feature outline), 'convex_hull' (dilated convex hull), or 'box' (bounding envelope + margin).",
                    required: false,
                },
                ToolParamSpec {
                    name: "masked_layer",
                    description: "Optional second polygon layer. When given, masks are clipped to their intersection with it (Intersecting Layers Masks); non-overlapping masks are dropped.",
                    required: false,
                },
                ToolParamSpec {
                    name: "id_field",
                    description: "Optional source field whose value is copied onto each mask as 'source_id'. Defaults to the source feature index.",
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
        let schema = layer.schema.clone();

        // Metric margins are converted to degrees on a geographic layer; projected
        // layers use the value verbatim in their own units.
        let geographic = matches!(layer.crs_epsg(), None | Some(4326));
        let dist = if geographic {
            prm.margin / METERS_PER_DEGREE
        } else {
            prm.margin
        };

        // Optionally load and dissolve the masked layer's polygons once, for the
        // Intersecting Layers Masks clip. Every part is concatenated into one
        // `MultiPolygon` and dissolved in a single boolean pass (union with an
        // empty polygon) — far cheaper than folding a pairwise union over many
        // features, and robust to overlapping masked polygons.
        let clip: Option<MultiPolygon> = match prm.masked_layer.as_deref() {
            Some(path) => {
                let masked = load_input_layer(path)?;
                let mut parts: Vec<Polygon> = Vec::new();
                for f in masked.features.iter() {
                    if let Some(mp) = f.geometry.as_ref().and_then(to_multipolygon) {
                        parts.extend(mp.0);
                    }
                }
                let acc = MultiPolygon(parts).union(&MultiPolygon::<f64>::new(vec![]));
                if acc.is_empty() {
                    ctx.progress.info(
                        "masked_layer has no polygon features; all mask intersections are empty",
                    );
                }
                Some(acc)
            }
            None => None,
        };

        ctx.progress.info(&format!(
            "{input_count} feature(s): building {} masks at margin {}{}",
            prm.mask_kind.as_str(),
            prm.margin,
            if clip.is_some() {
                " (clipped to masked_layer)"
            } else {
                ""
            }
        ));

        let mut out_layer = Layer::new(layer_name);
        out_layer.crs = layer_crs;
        out_layer.geom_type = Some(GeometryType::MultiPolygon);
        out_layer.add_field(FieldDef::new("source_fid", FieldType::Integer));
        out_layer.add_field(FieldDef::new("mask_kind", FieldType::Text));
        let carry_id = prm.id_field.is_some();
        if carry_id {
            out_layer.add_field(FieldDef::new("source_id", FieldType::Text));
        }

        let mut mask_count = 0usize;
        let mut clipped_out = 0usize;
        for (idx, feature) in layer.features.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref().and_then(to_geo_geometry) else {
                continue; // GeometryCollection or empty geometry: skip
            };
            let mut mask = build_mask(&geom, prm.mask_kind, dist);
            if mask.is_empty() {
                continue;
            }
            if let Some(clip) = &clip {
                mask = mask.intersection(clip);
                if mask.is_empty() {
                    clipped_out += 1;
                    continue;
                }
            }

            let source_id = prm.id_field.as_ref().map(|field| {
                feature
                    .get(&schema, field)
                    .map(field_value_string)
                    .unwrap_or_default()
            });
            let mut attrs: Vec<(&str, FieldValue)> = vec![
                ("source_fid", FieldValue::Integer(idx as i64)),
                ("mask_kind", FieldValue::Text(prm.mask_kind.as_str().into())),
            ];
            if let Some(sid) = source_id {
                attrs.push(("source_id", FieldValue::Text(sid)));
            }
            out_layer
                .add_feature(Some(multipolygon_to_geometry(&mask)), &attrs)
                .map_err(|e| ToolError::Execution(format!("failed writing mask: {e}")))?;
            mask_count += 1;
        }

        ctx.progress
            .info(&format!("wrote {mask_count} mask polygon(s)"));

        let feature_count = out_layer.len();
        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_count".to_string(), json!(input_count));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("mask_count".to_string(), json!(mask_count));
        outputs.insert("mask_kind".to_string(), json!(prm.mask_kind.as_str()));
        if clip.is_some() {
            outputs.insert("clipped_out_count".to_string(), json!(clipped_out));
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum MaskKind {
    Exact,
    ConvexHull,
    Box,
}

impl MaskKind {
    fn as_str(self) -> &'static str {
        match self {
            Self::Exact => "exact",
            Self::ConvexHull => "convex_hull",
            Self::Box => "box",
        }
    }
}

struct Params {
    margin: f64,
    mask_kind: MaskKind,
    masked_layer: Option<String>,
    id_field: Option<String>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let margin = opt_f64(args, "margin")?
        .ok_or_else(|| ToolError::Validation("missing required parameter 'margin'".to_string()))?;
    if !(margin > 0.0 && margin.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'margin' must be a positive number".to_string(),
        ));
    }
    let mask_kind = match parse_optional_str(args, "mask_kind")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("exact") => MaskKind::Exact,
        Some("convex_hull") | Some("convexhull") => MaskKind::ConvexHull,
        Some("box") | Some("envelope") => MaskKind::Box,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown mask_kind '{other}' (expected exact, convex_hull, or box)"
            )))
        }
    };
    let masked_layer = parse_optional_str(args, "masked_layer")?.map(str::to_string);
    let id_field = parse_optional_str(args, "id_field")?.map(str::to_string);
    Ok(Params {
        margin,
        mask_kind,
        masked_layer,
        id_field,
    })
}

/// Parses an optional numeric parameter, accepting a JSON number or a numeric
/// string (host UIs often post form values as strings).
fn opt_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
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

// ── Mask construction ─────────────────────────────────────────────────────────

/// Builds the mask polygon for one source geometry under the given kind and
/// dilation distance. Every result contains the source feature.
fn build_mask(geom: &GeoGeometry, kind: MaskKind, dist: f64) -> MultiPolygon {
    match kind {
        MaskKind::Exact => geom.buffer(dist),
        MaskKind::ConvexHull => {
            let pts = MultiPoint(collect_points(geom));
            // convex_hull() of >=3 non-collinear points is an area polygon.
            let buffered = pts.convex_hull().buffer(dist);
            if buffered.is_empty() {
                // Degenerate hull (a single point or collinear vertices has no
                // area): fall back to buffering the source geometry directly,
                // which still fully contains it (a disk / capsule).
                geom.buffer(dist)
            } else {
                buffered
            }
        }
        MaskKind::Box => match geom.bounding_rect() {
            Some(rect) => {
                let min = rect.min();
                let max = rect.max();
                let poly = Polygon::new(
                    LineString::new(vec![
                        GeoCoord {
                            x: min.x - dist,
                            y: min.y - dist,
                        },
                        GeoCoord {
                            x: max.x + dist,
                            y: min.y - dist,
                        },
                        GeoCoord {
                            x: max.x + dist,
                            y: max.y + dist,
                        },
                        GeoCoord {
                            x: min.x - dist,
                            y: max.y + dist,
                        },
                        GeoCoord {
                            x: min.x - dist,
                            y: min.y - dist,
                        },
                    ]),
                    vec![],
                );
                MultiPolygon(vec![poly])
            }
            None => MultiPolygon::new(vec![]),
        },
    }
}

/// Collects every coordinate of a geometry into `geo::Point`s (for the convex
/// hull). Covers all single- and multi-part types.
fn collect_points(geom: &GeoGeometry) -> Vec<Point> {
    let mut pts = Vec::new();
    let mut push = |x: f64, y: f64| pts.push(Point::new(x, y));
    match geom {
        GeoGeometry::Point(p) => push(p.x(), p.y()),
        GeoGeometry::MultiPoint(mp) => mp.iter().for_each(|p| push(p.x(), p.y())),
        GeoGeometry::LineString(ls) => ls.0.iter().for_each(|c| push(c.x, c.y)),
        GeoGeometry::MultiLineString(mls) => mls
            .iter()
            .flat_map(|ls| ls.0.iter())
            .for_each(|c| push(c.x, c.y)),
        GeoGeometry::Polygon(p) => push_polygon_points(p, &mut push),
        GeoGeometry::MultiPolygon(mp) => mp.iter().for_each(|p| push_polygon_points(p, &mut push)),
        _ => {}
    }
    pts
}

fn push_polygon_points(p: &Polygon, push: &mut impl FnMut(f64, f64)) {
    for c in p.exterior().0.iter() {
        push(c.x, c.y);
    }
    for ring in p.interiors() {
        for c in ring.0.iter() {
            push(c.x, c.y);
        }
    }
}

fn field_value_string(v: &FieldValue) -> String {
    match v {
        FieldValue::Text(s) => s.clone(),
        FieldValue::Date(s) | FieldValue::DateTime(s) => s.clone(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) => f.to_string(),
        FieldValue::Boolean(b) => b.to_string(),
        FieldValue::Null | FieldValue::Blob(_) => String::new(),
    }
}

// ── geo <-> wbvector geometry conversion ───────────────────────────────────

/// Converts a `wbvector` geometry to a `geo` geometry. Nested
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

/// Converts a polygonal `wbvector` geometry to a `geo` `MultiPolygon` (for the
/// masked-layer clip). Non-polygon geometries contribute nothing.
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
    use geo::Area;
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

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = FeatureOutlineMasksTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Coord> {
        vec![
            Coord::xy(x0, y0),
            Coord::xy(x1, y0),
            Coord::xy(x1, y1),
            Coord::xy(x0, y1),
        ]
    }

    fn mask_geo(layer: &Layer, idx: usize) -> GeoGeometry {
        to_geo_geometry(layer.features[idx].geometry.as_ref().unwrap()).unwrap()
    }

    fn area_of(layer: &Layer, idx: usize) -> f64 {
        match mask_geo(layer, idx) {
            GeoGeometry::Polygon(p) => p.unsigned_area(),
            GeoGeometry::MultiPolygon(mp) => mp.unsigned_area(),
            _ => 0.0,
        }
    }

    /// A projected point mask is a disk of area ~ pi*margin^2 and the source
    /// point falls inside it.
    #[test]
    fn exact_point_mask_is_a_disk_containing_source() {
        let mut layer = Layer::new("pts").with_crs_epsg(3857);
        layer
            .add_feature(Some(Geometry::point(0.0, 0.0)), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input, "margin": 10.0 }));
        assert_eq!(out.outputs["mask_count"], json!(1));
        let pi = std::f64::consts::PI;
        let a = area_of(&layer, 0);
        assert!((a - 100.0 * pi).abs() < 0.05 * 100.0 * pi, "disk area {a}");
        // Source point (0,0) is inside the mask.
        assert!(geo::Contains::contains(
            &mask_geo(&layer, 0),
            &Point::new(0.0, 0.0)
        ));
    }

    /// A single-point feature has a degenerate convex hull; the mask must still
    /// be a non-empty disk containing the source (regression for the fallback).
    #[test]
    fn convex_hull_of_a_point_is_a_disk() {
        let mut layer = Layer::new("pts").with_crs_epsg(3857);
        layer
            .add_feature(Some(Geometry::point(5.0, 5.0)), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);
        let (out, layer) =
            run_tool(json!({ "input": input, "margin": 4.0, "mask_kind": "convex_hull" }));
        assert_eq!(out.outputs["mask_count"], json!(1));
        let pi = std::f64::consts::PI;
        assert!((area_of(&layer, 0) - 16.0 * pi).abs() < 0.05 * 16.0 * pi);
        assert!(geo::Contains::contains(
            &mask_geo(&layer, 0),
            &Point::new(5.0, 5.0)
        ));
    }

    /// For a polygon input, every mask kind fully contains the source polygon
    /// and the mask area exceeds the source area (source ⊆ mask).
    #[test]
    fn every_kind_contains_the_source_polygon() {
        for kind in ["exact", "convex_hull", "box"] {
            let mut layer = Layer::new("poly").with_crs_epsg(3857);
            layer
                .add_feature(
                    Some(Geometry::polygon(rect(0.0, 0.0, 10.0, 10.0), vec![])),
                    &[],
                )
                .unwrap();
            let id = memory_store::put_vector(layer);
            let input = memory_store::make_vector_memory_path(&id);
            let (_, layer) = run_tool(json!({ "input": input, "margin": 5.0, "mask_kind": kind }));
            let mask = mask_geo(&layer, 0);
            let source = geo::Rect::new(GeoCoord { x: 0.0, y: 0.0 }, GeoCoord { x: 10.0, y: 10.0 })
                .to_polygon();
            assert!(
                geo::Contains::contains(&mask, &source),
                "{kind} mask must contain the source polygon"
            );
            // Mask area exceeds the 100-unit source by roughly the margin ring.
            assert!(area_of(&layer, 0) > 100.0, "{kind} mask larger than source");
        }
    }

    /// A `box` mask around a 10x10 square with margin 5 is a 20x20 rectangle
    /// (area 400).
    #[test]
    fn box_mask_is_the_expanded_envelope() {
        let mut layer = Layer::new("poly").with_crs_epsg(3857);
        layer
            .add_feature(
                Some(Geometry::polygon(rect(0.0, 0.0, 10.0, 10.0), vec![])),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);
        let (_, layer) = run_tool(json!({ "input": input, "margin": 5.0, "mask_kind": "box" }));
        // (0-5 .. 10+5) on each axis -> 20 x 20 = 400.
        assert!((area_of(&layer, 0) - 400.0).abs() < 1e-6);
    }

    /// Intersecting Layers Masks: the mask is clipped to the masked layer; a
    /// feature whose mask misses the masked layer produces no output.
    #[test]
    fn masked_layer_clips_and_drops_non_overlapping() {
        // Two points: one near the masked square, one far away.
        let mut pts = Layer::new("pts").with_crs_epsg(3857);
        pts.add_feature(Some(Geometry::point(0.0, 0.0)), &[])
            .unwrap();
        pts.add_feature(Some(Geometry::point(1000.0, 1000.0)), &[])
            .unwrap();
        let pts_id = memory_store::put_vector(pts);
        let input = memory_store::make_vector_memory_path(&pts_id);

        // Masked layer: a big square covering the origin only.
        let mut masked = Layer::new("masked").with_crs_epsg(3857);
        masked
            .add_feature(
                Some(Geometry::polygon(rect(-50.0, -50.0, 50.0, 50.0), vec![])),
                &[],
            )
            .unwrap();
        let masked_id = memory_store::put_vector(masked);
        let masked_path = memory_store::make_vector_memory_path(&masked_id);

        let (out, layer) = run_tool(json!({
            "input": input,
            "margin": 10.0,
            "masked_layer": masked_path,
        }));
        // Only the origin point's mask overlaps the masked square.
        assert_eq!(out.outputs["mask_count"], json!(1));
        assert_eq!(out.outputs["clipped_out_count"], json!(1));
        // The clipped mask (full disk, inside the square) keeps ~pi*100 area.
        let pi = std::f64::consts::PI;
        assert!((area_of(&layer, 0) - 100.0 * pi).abs() < 0.05 * 100.0 * pi);
    }

    /// `id_field` is copied onto each mask as `source_id`.
    #[test]
    fn id_field_is_carried_onto_masks() {
        let mut layer = Layer::new("pts").with_crs_epsg(3857);
        layer.add_field(FieldDef::new("name", FieldType::Text));
        layer
            .add_feature(Some(Geometry::point(0.0, 0.0)), &[("name", "well".into())])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);
        let (_, layer) = run_tool(json!({ "input": input, "margin": 5.0, "id_field": "name" }));
        let sid = layer.features[0].get(&layer.schema, "source_id").unwrap();
        assert_eq!(sid, &FieldValue::Text("well".into()));
    }

    /// GeometryCollection features are skipped rather than erroring.
    #[test]
    fn skips_geometry_collections() {
        let mut layer = Layer::new("mixed").with_crs_epsg(3857);
        layer
            .add_feature(Some(Geometry::point(0.0, 0.0)), &[])
            .unwrap();
        layer.push(wbvector::Feature {
            fid: 1,
            geometry: Some(Geometry::GeometryCollection(vec![])),
            attributes: vec![],
        });
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);
        let (out, _) = run_tool(json!({ "input": input, "margin": 5.0 }));
        // Only the point produced a mask.
        assert_eq!(out.outputs["mask_count"], json!(1));
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = FeatureOutlineMasksTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input");
        assert!(
            bad(json!({ "input": "x.geojson" })).is_err(),
            "missing margin"
        );
        assert!(bad(json!({ "input": "x.geojson", "margin": 0 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "margin": -5 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "margin": 5, "mask_kind": "blob" })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "margin": 5 })).is_ok());
        assert!(
            bad(json!({ "input": "x.geojson", "margin": "5.0", "mask_kind": "box" })).is_ok(),
            "numeric strings ok"
        );
    }
}
