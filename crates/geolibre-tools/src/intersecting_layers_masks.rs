//! GeoLibre tool: conditional cartographic masks where two layers intersect.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Intersecting Layers Masks*
//! (Cartography → Masking Tools): build a mask polygon at a `margin` around each
//! feature of a `masked_layer`, but **only where that feature actually
//! intersects** a feature of a second `masking_layer`. This is the
//! coverage-aware, pairwise sibling of GeoLibre's `feature_outline_masks` (which
//! masks *every* feature): use it to knock out ("halo") symbology only where two
//! layers conflict — e.g. mask contour labels only where roads cross them, so
//! non-conflicting labels keep their full context.
//!
//! For each feature of the masked layer we test intersection against the masking
//! layer (bounding-box prune, then exact `geo::Intersects`). Only matching
//! features emit a mask; non-intersecting features are dropped. The mask shape
//! mirrors ArcGIS and `feature_outline_masks`:
//!
//! - `exact` (default) — the masked feature's outline dilated by the margin
//!   (`geo`'s `Buffer`, rounded joins). A disk for points, a capsule for lines,
//!   the grown footprint for polygons — the tightest mask.
//! - `convex_hull` — the convex hull of the feature's vertices, dilated by the
//!   margin. Cheaper and smoother for busy geometry.
//! - `box` — the feature's bounding envelope expanded by `margin` on every side
//!   (a clean rectangle). The coarsest, fastest mask.
//!
//! Every mask fully contains its source feature (a dilation of a superset of the
//! feature), so `source ⊆ mask` always holds. Unlike the clip performed by
//! `feature_outline_masks.masked_layer`, the mask here is **not** trimmed to the
//! masking layer — the whole feature is masked once a conflict is detected,
//! which is what ArcGIS does. Each output polygon carries the 0-based index of
//! its source feature in `source_fid`, the `mask_kind`, and — when `id_field` is
//! given — the source id value in `source_id`.
//!
//! **Units.** `margin` is read in metres for a geographic layer (EPSG:4326 or an
//! unknown CRS) and converted to degrees at ~111320 m/°; for a projected layer
//! it is taken directly in the layer's CRS units. Buffering happens in the
//! masked layer's coordinate space.
//!
//! Scope for v1: the two layers are assumed to share a CRS (no reprojection —
//! GeoLibre carries no PROJ); intersection pruning uses a linear bounding-box
//! scan rather than a persistent spatial index; source attributes other than the
//! chosen `id_field` are not copied onto the masks; and the geographic
//! metre→degree conversion uses a single spherical scale factor.

use std::collections::BTreeMap;

use geo::{
    BoundingRect, Buffer, ConvexHull, Coord as GeoCoord, Geometry as GeoGeometry, HasDimensions,
    Intersects, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon, Rect,
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

pub struct IntersectingLayersMasksTool;

impl Tool for IntersectingLayersMasksTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "intersecting_layers_masks",
            display_name: "Intersecting Layers Masks",
            summary: "Generate cartographic mask polygons at a margin around features of a masked layer, but only where they intersect a second masking layer (exact outline, convex hull, or bounding box).",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "masked_layer",
                    description: "Masked vector layer (any geometry): each feature that intersects the masking layer gets a mask. Format auto-detected (or an in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "masking_layer",
                    description: "Masking vector layer (any geometry): a masked feature is masked only where it intersects a feature of this layer.",
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
                    name: "id_field",
                    description: "Optional field on the masked layer whose value is copied onto each mask as 'source_id'. Defaults to the source feature index.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["masked_layer", "masking_layer"] {
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
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let masked_path = req_str(args, "masked_layer")?;
        let masking_path = req_str(args, "masking_layer")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let masked = load_input_layer(masked_path)?;
        let layer_name = masked.name.clone();
        let layer_crs = masked.crs.clone();
        let masked_count = masked.len();
        let schema = masked.schema.clone();

        // Metric margins are converted to degrees on a geographic layer; projected
        // layers use the value verbatim in their own units.
        let geographic = matches!(masked.crs_epsg(), None | Some(4326));
        let dist = if geographic {
            prm.margin / METERS_PER_DEGREE
        } else {
            prm.margin
        };

        // Load the masking layer once and cache each feature's bounding rect for a
        // linear bbox prune before the exact `Intersects` test.
        let masking = load_input_layer(masking_path)?;
        let masking_geoms: Vec<(Rect, GeoGeometry)> = masking
            .features
            .iter()
            .filter_map(|f| f.geometry.as_ref().and_then(to_geo_geometry))
            .filter(|g| !g.is_empty())
            .filter_map(|g| g.bounding_rect().map(|r| (r, g)))
            .collect();
        if masking_geoms.is_empty() {
            ctx.progress
                .info("masking_layer has no usable features; no masks will be produced");
        }

        ctx.progress.info(&format!(
            "{masked_count} masked feature(s) vs {} masking feature(s): building {} masks at margin {}",
            masking_geoms.len(),
            prm.mask_kind.as_str(),
            prm.margin
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
        let mut non_intersecting = 0usize;
        for (idx, feature) in masked.features.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref().and_then(to_geo_geometry) else {
                continue; // GeometryCollection or empty geometry: skip
            };
            let Some(mrect) = geom.bounding_rect() else {
                continue;
            };
            // Conditional test: mask only if the masked feature intersects any
            // masking feature. Bounding-box prune first, then exact geometry test.
            let hit = masking_geoms
                .iter()
                .any(|(rect, mg)| mrect.intersects(rect) && geom.intersects(mg));
            if !hit {
                non_intersecting += 1;
                continue;
            }

            let mask = build_mask(&geom, prm.mask_kind, dist);
            if mask.is_empty() {
                continue;
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

        ctx.progress.info(&format!(
            "wrote {mask_count} mask polygon(s); {non_intersecting} feature(s) had no intersection"
        ));

        let feature_count = out_layer.len();
        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("masked_count".to_string(), json!(masked_count));
        outputs.insert("masking_count".to_string(), json!(masking_geoms.len()));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("mask_count".to_string(), json!(mask_count));
        outputs.insert(
            "non_intersecting_count".to_string(),
            json!(non_intersecting),
        );
        outputs.insert("mask_kind".to_string(), json!(prm.mask_kind.as_str()));
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
    id_field: Option<String>,
}

fn req_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
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
    let id_field = parse_optional_str(args, "id_field")?.map(str::to_string);
    Ok(Params {
        margin,
        mask_kind,
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

    fn store(layer: Layer) -> String {
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = IntersectingLayersMasksTool.run(&args, &ctx()).unwrap();
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

    /// Core property: a mask is emitted only for masked features that intersect
    /// the masking layer. Two points, a masking line crossing only one of them:
    /// exactly one mask, and it wraps the intersecting point.
    #[test]
    fn masks_only_intersecting_features() {
        let mut pts = Layer::new("labels").with_crs_epsg(3857);
        pts.add_feature(Some(Geometry::point(0.0, 0.0)), &[])
            .unwrap();
        pts.add_feature(Some(Geometry::point(1000.0, 1000.0)), &[])
            .unwrap();
        let masked = store(pts);

        // A vertical line through x=0 crosses the origin point but not (1000,1000).
        let mut lines = Layer::new("roads").with_crs_epsg(3857);
        lines
            .add_feature(
                Some(Geometry::LineString(vec![
                    Coord::xy(0.0, -50.0),
                    Coord::xy(0.0, 50.0),
                ])),
                &[],
            )
            .unwrap();
        let masking = store(lines);

        let (out, layer) = run_tool(json!({
            "masked_layer": masked,
            "masking_layer": masking,
            "margin": 10.0,
        }));
        assert_eq!(out.outputs["mask_count"], json!(1));
        assert_eq!(out.outputs["non_intersecting_count"], json!(1));
        assert_eq!(out.outputs["feature_count"], json!(1));
        // The one emitted mask surrounds the origin (its source feature).
        assert_eq!(
            layer.features[0].get(&layer.schema, "source_fid").unwrap(),
            &FieldValue::Integer(0)
        );
        assert!(geo::Contains::contains(
            &mask_geo(&layer, 0),
            &Point::new(0.0, 0.0)
        ));
    }

    /// Non-matching pass-through: if nothing intersects, no masks are produced.
    #[test]
    fn no_intersection_yields_no_masks() {
        let mut pts = Layer::new("labels").with_crs_epsg(3857);
        pts.add_feature(Some(Geometry::point(0.0, 0.0)), &[])
            .unwrap();
        let masked = store(pts);
        let mut far = Layer::new("roads").with_crs_epsg(3857);
        far.add_feature(Some(Geometry::point(9999.0, 9999.0)), &[])
            .unwrap();
        let masking = store(far);
        let (out, _) = run_tool(json!({
            "masked_layer": masked,
            "masking_layer": masking,
            "margin": 5.0,
        }));
        assert_eq!(out.outputs["mask_count"], json!(0));
        assert_eq!(out.outputs["non_intersecting_count"], json!(1));
    }

    /// Every mask kind fully contains its intersecting source polygon
    /// (source ⊆ mask) and exceeds the source area.
    #[test]
    fn every_kind_contains_source() {
        for kind in ["exact", "convex_hull", "box"] {
            let mut polys = Layer::new("poly").with_crs_epsg(3857);
            polys
                .add_feature(
                    Some(Geometry::polygon(rect(0.0, 0.0, 10.0, 10.0), vec![])),
                    &[],
                )
                .unwrap();
            let masked = store(polys);
            // A masking polygon overlapping the source.
            let mut mask_src = Layer::new("mask").with_crs_epsg(3857);
            mask_src
                .add_feature(
                    Some(Geometry::polygon(rect(5.0, 5.0, 20.0, 20.0), vec![])),
                    &[],
                )
                .unwrap();
            let masking = store(mask_src);

            let (out, layer) = run_tool(json!({
                "masked_layer": masked,
                "masking_layer": masking,
                "margin": 5.0,
                "mask_kind": kind,
            }));
            assert_eq!(out.outputs["mask_count"], json!(1), "{kind}");
            let mask = mask_geo(&layer, 0);
            let source = geo::Rect::new(GeoCoord { x: 0.0, y: 0.0 }, GeoCoord { x: 10.0, y: 10.0 })
                .to_polygon();
            assert!(
                geo::Contains::contains(&mask, &source),
                "{kind} contains source"
            );
            assert!(area_of(&layer, 0) > 100.0, "{kind} larger than source");
        }
    }

    /// A `box` mask around a 10x10 square with margin 5 is a 20x20 rectangle
    /// (area 400) when it intersects.
    #[test]
    fn box_mask_is_the_expanded_envelope() {
        let mut polys = Layer::new("poly").with_crs_epsg(3857);
        polys
            .add_feature(
                Some(Geometry::polygon(rect(0.0, 0.0, 10.0, 10.0), vec![])),
                &[],
            )
            .unwrap();
        let masked = store(polys);
        let mut mask_src = Layer::new("mask").with_crs_epsg(3857);
        mask_src
            .add_feature(Some(Geometry::point(5.0, 5.0)), &[])
            .unwrap();
        let masking = store(mask_src);
        let (_, layer) = run_tool(json!({
            "masked_layer": masked,
            "masking_layer": masking,
            "margin": 5.0,
            "mask_kind": "box",
        }));
        assert!((area_of(&layer, 0) - 400.0).abs() < 1e-6);
    }

    /// Touching-only geometries count as intersecting (shared boundary).
    #[test]
    fn touching_counts_as_intersection() {
        let mut a = Layer::new("a").with_crs_epsg(3857);
        a.add_feature(
            Some(Geometry::polygon(rect(0.0, 0.0, 10.0, 10.0), vec![])),
            &[],
        )
        .unwrap();
        let masked = store(a);
        // Shares the edge x=10.
        let mut b = Layer::new("b").with_crs_epsg(3857);
        b.add_feature(
            Some(Geometry::polygon(rect(10.0, 0.0, 20.0, 10.0), vec![])),
            &[],
        )
        .unwrap();
        let masking = store(b);
        let (out, _) = run_tool(json!({
            "masked_layer": masked,
            "masking_layer": masking,
            "margin": 1.0,
        }));
        assert_eq!(out.outputs["mask_count"], json!(1));
    }

    /// `id_field` is copied onto each mask as `source_id`.
    #[test]
    fn id_field_is_carried_onto_masks() {
        let mut pts = Layer::new("labels").with_crs_epsg(3857);
        pts.add_field(FieldDef::new("name", FieldType::Text));
        pts.add_feature(Some(Geometry::point(0.0, 0.0)), &[("name", "peak".into())])
            .unwrap();
        let masked = store(pts);
        let mut m = Layer::new("m").with_crs_epsg(3857);
        m.add_feature(Some(Geometry::point(0.0, 0.0)), &[]).unwrap();
        let masking = store(m);
        let (_, layer) = run_tool(json!({
            "masked_layer": masked,
            "masking_layer": masking,
            "margin": 5.0,
            "id_field": "name",
        }));
        let sid = layer.features[0].get(&layer.schema, "source_id").unwrap();
        assert_eq!(sid, &FieldValue::Text("peak".into()));
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = IntersectingLayersMasksTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing both layers");
        assert!(
            bad(json!({ "masked_layer": "a.geojson", "margin": 5 })).is_err(),
            "missing masking_layer"
        );
        assert!(
            bad(json!({ "masking_layer": "b.geojson", "margin": 5 })).is_err(),
            "missing masked_layer"
        );
        assert!(
            bad(json!({ "masked_layer": "a.geojson", "masking_layer": "b.geojson" })).is_err(),
            "missing margin"
        );
        assert!(bad(
            json!({ "masked_layer": "a.geojson", "masking_layer": "b.geojson", "margin": 0 })
        )
        .is_err());
        assert!(bad(
            json!({ "masked_layer": "a.geojson", "masking_layer": "b.geojson", "margin": -3 })
        )
        .is_err());
        assert!(bad(json!({
            "masked_layer": "a.geojson", "masking_layer": "b.geojson",
            "margin": 5, "mask_kind": "blob"
        }))
        .is_err());
        assert!(bad(json!({
            "masked_layer": "a.geojson", "masking_layer": "b.geojson", "margin": 5
        }))
        .is_ok());
        assert!(
            bad(json!({
                "masked_layer": "a.geojson", "masking_layer": "b.geojson",
                "margin": "5.0", "mask_kind": "box"
            }))
            .is_ok(),
            "numeric strings ok"
        );
    }
}
