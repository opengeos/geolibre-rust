//! GeoLibre tool: flag narrow polygons (slivers and thin necks).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Identify Narrow Polygons*
//! (Topographic Production): a QA step that finds polygons which are thinner
//! than a width tolerance somewhere — full slivers (long, thin, small width but
//! not necessarily small *area*) and pinched necks joining two lobes.
//!
//! The bundled whitebox-wasm suite can drop polygons by area
//! (`filter_vector_features_by_area`), but area alone misses a long, thin
//! polygon: a 500 m × 2 m sliver has a large area yet is clearly narrow. This
//! tool measures *local width* instead, which is what the downstream cleanup
//! tools (`eliminate_polygons`, `collapse_hydro_polygon`) actually key on.
//!
//! ## How narrowness is measured
//!
//! A morphological **opening** at radius `h = tolerance/2`: erode each polygon
//! inward by `h`, then dilate back by `h`. The opening keeps exactly the parts
//! of the polygon wide enough to hold a disk of radius `h` (≥ `tolerance` wide),
//! and drops everything thinner. The area the opening removes — the *narrow
//! area* — is therefore the total area of the too-thin portions:
//!
//! - a fully thin sliver opens to nothing (narrow area ≈ its whole area),
//! - a pinched neck loses the neck, and
//! - a fat polygon loses only minor convex-corner rounding.
//!
//! A polygon is flagged `is_narrow` when its narrow area exceeds
//! `min_narrow_area` (default `tolerance²`, which also absorbs the corner
//! rounding). This is far more stable on real raster-derived coverages than raw
//! part counting: a vertex where the boundary self-touches contributes ~zero
//! area to the opening difference, so it raises no false positive.
//!
//! Each feature is annotated with `is_narrow`, `narrow_area`, and an estimated
//! `min_width = 2·h*` (bisection: the smallest half-distance whose opening
//! removes more than the floor — the width scale of the narrowest sliver or
//! neck). v1 scope: a purely one-sided narrow *protrusion* off an otherwise fat
//! body whose removed area stays under the floor is not flagged; the tool
//! targets slivers and necks, the dominant coverage-QA cases. The opening is
//! `geo`'s `Buffer` (pure Rust, no GEOS). Non-polygon features pass through
//! untouched with `is_narrow = 0`.

use std::collections::BTreeMap;

use geo::{Area, Buffer, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct IdentifyNarrowPolygonsTool;

impl Tool for IdentifyNarrowPolygonsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "identify_narrow_polygons",
            display_name: "Identify Narrow Polygons",
            summary: "Flag polygons that are narrower than a width tolerance somewhere (slivers and thin necks), annotating each with an is_narrow flag and an estimated minimum width.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon vector file path, format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "width_tolerance",
                    description: "Width threshold in CRS units. A polygon is narrow where it is thinner than this. Required.",
                    required: true,
                },
                ToolParamSpec {
                    name: "min_narrow_area",
                    description: "Minimum narrow-portion area (CRS units squared) required to flag a polygon; smaller narrow areas are ignored as noise. Default width_tolerance squared, which also absorbs the opening's minor corner rounding.",
                    required: false,
                },
                ToolParamSpec {
                    name: "narrow_only",
                    description: "When true, drop non-narrow features and output only the flagged narrow polygons. Default false (annotate and keep all features).",
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
        let input_count = layer.len();

        // Extend the schema with the three annotation fields (idempotent: skip
        // if a field of that name already exists so we never clobber user data).
        let mut schema = layer.schema.clone();
        for (name, ty) in [
            ("is_narrow", FieldType::Integer),
            ("min_width", FieldType::Float),
            ("narrow_area", FieldType::Float),
        ] {
            if schema.field_index(name).is_none() {
                schema.add_field(FieldDef::new(name, ty));
            }
        }
        let flag_idx = schema.field_index("is_narrow").unwrap();
        let width_idx = schema.field_index("min_width").unwrap();
        let area_idx = schema.field_index("narrow_area").unwrap();

        ctx.progress.info(&format!(
            "scanning {input_count} feature(s) for narrow polygons"
        ));

        let mut narrow_count = 0usize;
        let mut total_narrow_area = 0.0f64;
        let mut out_features = Vec::with_capacity(input_count);
        for feature in layer.features.into_iter() {
            let mut feature = feature;
            // Grow the attribute vector to the widened schema.
            feature.attributes.resize(schema.len(), FieldValue::Null);

            let info = match feature.geometry.as_ref().and_then(to_multipolygon) {
                Some(mp) => classify(&mp, prm.half_tolerance, prm.min_narrow_area),
                None => NarrowInfo {
                    is_narrow: false,
                    min_width: None,
                    narrow_area: 0.0,
                }, // non-polygon: pass through, not narrow
            };
            if info.is_narrow {
                narrow_count += 1;
                total_narrow_area += info.narrow_area;
            }
            feature.set_by_index(flag_idx, FieldValue::Integer(info.is_narrow as i64));
            feature.set_by_index(
                width_idx,
                match info.min_width {
                    Some(w) => FieldValue::Float(w),
                    None => FieldValue::Null,
                },
            );
            feature.set_by_index(area_idx, FieldValue::Float(info.narrow_area));

            if prm.narrow_only && !info.is_narrow {
                continue;
            }
            feature.fid = out_features.len() as u64;
            out_features.push(feature);
        }

        ctx.progress.info(&format!(
            "flagged {narrow_count} narrow polygon(s) of {input_count} feature(s)"
        ));

        let mut out_layer = wbvector::Layer::new(layer.name);
        out_layer.schema = schema;
        out_layer.crs = layer.crs;
        out_layer.geom_type = layer.geom_type;
        out_layer.features = out_features;

        let feature_count = out_layer.len();
        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_count".to_string(), json!(input_count));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("narrow_count".to_string(), json!(narrow_count));
        outputs.insert("total_narrow_area".to_string(), json!(total_narrow_area));
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    half_tolerance: f64,
    min_narrow_area: f64,
    narrow_only: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let tolerance = parse_optional_f64(args, "width_tolerance")?.ok_or_else(|| {
        ToolError::Validation("missing required parameter 'width_tolerance'".to_string())
    })?;
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'width_tolerance' must be a positive number".to_string(),
        ));
    }
    // Default floor = tolerance²: one tolerance-square of removed area. This also
    // absorbs the minor convex-corner rounding of the morphological opening
    // (~0.86·(tolerance/2)² per fat polygon), so square blocks are not flagged.
    let min_narrow_area =
        parse_optional_f64(args, "min_narrow_area")?.unwrap_or(tolerance * tolerance);
    if min_narrow_area < 0.0 || !min_narrow_area.is_finite() {
        return Err(ToolError::Validation(
            "parameter 'min_narrow_area' must be a non-negative number".to_string(),
        ));
    }
    let narrow_only = parse_optional_bool(args, "narrow_only")?.unwrap_or(false);
    Ok(Params {
        half_tolerance: tolerance / 2.0,
        min_narrow_area,
        narrow_only,
    })
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

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
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
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
}

// ── Narrowness test (morphological opening) ─────────────────────────────────

/// Result of classifying one polygonal geometry.
struct NarrowInfo {
    is_narrow: bool,
    /// Estimated narrowest width in CRS units (only when narrow).
    min_width: Option<f64>,
    /// Area lost to the morphological opening — the total area of thin portions.
    narrow_area: f64,
}

/// Classifies a polygonal geometry using a morphological opening.
///
/// The opening `O = (P ⊖ h) ⊕ h` (erode then dilate by `h = tolerance/2`)
/// removes exactly the portions of `P` too thin to contain a disk of radius `h`
/// — i.e. narrower than the tolerance. The *narrow area* is `area(P) − area(O)`:
/// a fully thin sliver opens to nothing (narrow area ≈ its full area), a pinched
/// neck loses the neck, and a fat polygon loses only minor corner rounding. The
/// polygon is flagged when the narrow area exceeds `min_narrow_area`.
///
/// This is far more stable on real (raster-derived) coverages than raw part
/// counting: a vertex where the boundary self-touches contributes ~zero area to
/// the opening difference, so it does not raise a false positive.
///
/// `min_width` is `2·h*` for the smallest half-distance `h*` (bisection) at
/// which the narrow area first exceeds the floor — the width scale of the
/// narrowest sliver or neck. `narrow_area` is monotone increasing in `h`, so the
/// bisection is well posed.
fn classify(mp: &MultiPolygon, half_tol: f64, min_narrow_area: f64) -> NarrowInfo {
    let base_area = mp.unsigned_area();
    if base_area <= 0.0 {
        return NarrowInfo {
            is_narrow: false,
            min_width: None,
            narrow_area: 0.0,
        };
    }
    let narrow_area = narrow_area_at(mp, half_tol, base_area);
    if narrow_area <= min_narrow_area {
        return NarrowInfo {
            is_narrow: false,
            min_width: None,
            narrow_area,
        };
    }
    // Bisection for the smallest half-distance whose opening removes more than
    // the floor. narrow_area(h) is monotone increasing in h.
    let mut lo = 0.0f64; // opening removes ~nothing
    let mut hi = half_tol; // opening removes > floor
    for _ in 0..20 {
        let mid = 0.5 * (lo + hi);
        if narrow_area_at(mp, mid, base_area) > min_narrow_area {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    NarrowInfo {
        is_narrow: true,
        min_width: Some(2.0 * hi),
        narrow_area,
    }
}

/// Area removed by opening `mp` at radius `h`: `area(mp) − area(opening)`.
fn narrow_area_at(mp: &MultiPolygon, h: f64, base_area: f64) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    let eroded = mp.buffer(-h);
    if eroded.0.is_empty() {
        return base_area; // nothing survives erosion — the whole polygon is thin
    }
    let opened = eroded.buffer(h);
    (base_area - opened.unsigned_area()).max(0.0)
}

// ── geo <-> wbvector geometry conversion ───────────────────────────────────

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

/// Builds a `geo` ring, dropping consecutive-duplicate and collinear vertices.
///
/// Raster-derived coverages carry long runs of collinear rectilinear vertices
/// (median tens, up to ~1300 per ring here). `geo`'s `Buffer` mis-offsets such
/// dense rings — the erode/dilate opening then fails to round-trip and reports
/// spurious "narrow" area even for fat polygons. Reducing each ring to its
/// genuine corners first makes the opening stable without changing the shape.
fn ring_to_linestring(ring: &Ring) -> LineString {
    let src: Vec<GeoCoord> = ring
        .coords()
        .iter()
        .map(|c| GeoCoord { x: c.x, y: c.y })
        .collect();
    LineString::new(drop_collinear(src))
}

/// Removes consecutive duplicates and vertices collinear with their neighbours
/// (relative cross-product test, so it works at any coordinate scale). The ring
/// is treated as closed. Rings that collapse below 3 vertices are returned as-is.
fn drop_collinear(pts: Vec<GeoCoord>) -> Vec<GeoCoord> {
    let n = pts.len();
    if n < 4 {
        return pts;
    }
    let mut out: Vec<GeoCoord> = Vec::with_capacity(n);
    for i in 0..n {
        let prev = last_kept(&out, &pts, i);
        let cur = pts[i];
        let next = pts[(i + 1) % n];
        let (ax, ay) = (cur.x - prev.x, cur.y - prev.y);
        let (bx, by) = (next.x - cur.x, next.y - cur.y);
        let la = ax.hypot(ay);
        let lb = bx.hypot(by);
        if la == 0.0 {
            continue; // coincident with the previous kept vertex
        }
        // Collinear when the normalised cross product is ~0.
        let cross = ax * by - ay * bx;
        if lb > 0.0 && cross.abs() <= 1e-9 * la * lb {
            continue;
        }
        out.push(cur);
    }
    if out.len() < 3 {
        pts
    } else {
        out
    }
}

fn last_kept(out: &[GeoCoord], pts: &[GeoCoord], i: usize) -> GeoCoord {
    match out.last() {
        Some(p) => *p,
        None => pts[(i + pts.len() - 1) % pts.len()],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Coord, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Coord> {
        vec![
            Coord::xy(x0, y0),
            Coord::xy(x1, y0),
            Coord::xy(x1, y1),
            Coord::xy(x0, y1),
        ]
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = IdentifyNarrowPolygonsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn narrow_flag(layer: &Layer, fid: usize) -> i64 {
        layer.features[fid]
            .get(&layer.schema, "is_narrow")
            .unwrap()
            .as_i64()
            .unwrap()
    }

    /// A wide block is not narrow; a long thin sliver of equal (large) area is.
    #[test]
    fn flags_thin_sliver_but_not_wide_block() {
        let mut layer = Layer::new("polys");
        // Wide 20x20 block (area 400).
        layer
            .add_feature(
                Some(Geometry::polygon(rect(0.0, 0.0, 20.0, 20.0), vec![])),
                &[],
            )
            .unwrap();
        // Thin 100x2 sliver (area 200 — large area, tiny width).
        layer
            .add_feature(
                Some(Geometry::polygon(rect(100.0, 0.0, 200.0, 2.0), vec![])),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input, "width_tolerance": 5.0 }));
        assert_eq!(out.outputs["narrow_count"], json!(1));
        assert_eq!(narrow_flag(&layer, 0), 0, "wide block not narrow");
        assert_eq!(narrow_flag(&layer, 1), 1, "thin sliver narrow");
        // The sliver is 2 units wide; the estimate should be close to 2.
        let w = layer.features[1]
            .get(&layer.schema, "min_width")
            .unwrap()
            .as_f64()
            .unwrap();
        assert!((w - 2.0).abs() < 0.6, "estimated width {w} should be ~2");
    }

    /// A dumbbell — two fat lobes joined by a thin neck — splits under erosion
    /// and is flagged even though neither lobe alone is narrow.
    #[test]
    fn flags_thin_neck() {
        let mut layer = Layer::new("polys");
        // Two 10x10 lobes joined by a 1-unit-tall neck along y in [4.5,5.5].
        let coords = vec![
            Coord::xy(0.0, 0.0),
            Coord::xy(10.0, 0.0),
            Coord::xy(10.0, 4.5),
            Coord::xy(20.0, 4.5),
            Coord::xy(20.0, 0.0),
            Coord::xy(30.0, 0.0),
            Coord::xy(30.0, 10.0),
            Coord::xy(20.0, 10.0),
            Coord::xy(20.0, 5.5),
            Coord::xy(10.0, 5.5),
            Coord::xy(10.0, 10.0),
            Coord::xy(0.0, 10.0),
        ];
        layer
            .add_feature(Some(Geometry::polygon(coords, vec![])), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input, "width_tolerance": 3.0 }));
        assert_eq!(out.outputs["narrow_count"], json!(1));
        assert_eq!(narrow_flag(&layer, 0), 1);
    }

    /// A non-polygon passes through, un-flagged.
    #[test]
    fn passes_points_through() {
        let mut layer = Layer::new("mixed");
        layer
            .add_feature(Some(Geometry::point(1.0, 2.0)), &[])
            .unwrap();
        layer
            .add_feature(
                Some(Geometry::polygon(rect(0.0, 0.0, 1.0, 50.0), vec![])),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input, "width_tolerance": 5.0 }));
        assert_eq!(out.outputs["feature_count"], json!(2));
        assert_eq!(narrow_flag(&layer, 0), 0, "point not narrow");
        assert_eq!(narrow_flag(&layer, 1), 1, "1-wide strip narrow");
    }

    /// `narrow_only` drops the non-narrow features.
    #[test]
    fn narrow_only_drops_wide() {
        let mut layer = Layer::new("polys");
        layer
            .add_feature(
                Some(Geometry::polygon(rect(0.0, 0.0, 20.0, 20.0), vec![])),
                &[],
            )
            .unwrap();
        layer
            .add_feature(
                Some(Geometry::polygon(rect(100.0, 0.0, 200.0, 2.0), vec![])),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, _) =
            run_tool(json!({ "input": input, "width_tolerance": 5.0, "narrow_only": true }));
        assert_eq!(out.outputs["feature_count"], json!(1));
        assert_eq!(out.outputs["narrow_count"], json!(1));
    }

    /// `min_narrow_area` suppresses a tiny speck.
    #[test]
    fn min_area_floor_suppresses_speck() {
        let mut layer = Layer::new("polys");
        // 4x0.5 speck, area 2.
        layer
            .add_feature(
                Some(Geometry::polygon(rect(0.0, 0.0, 4.0, 0.5), vec![])),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, _) =
            run_tool(json!({ "input": input, "width_tolerance": 5.0, "min_narrow_area": 10.0 }));
        assert_eq!(
            out.outputs["narrow_count"],
            json!(0),
            "speck below area floor"
        );
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = IdentifyNarrowPolygonsTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input");
        assert!(
            bad(json!({ "input": "x.geojson" })).is_err(),
            "missing width_tolerance"
        );
        assert!(bad(json!({ "input": "x.geojson", "width_tolerance": 0 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "width_tolerance": -1 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "width_tolerance": "5.0" })).is_ok());
        assert!(
            bad(json!({ "input": "x.geojson", "width_tolerance": 5, "min_narrow_area": -1 }))
                .is_err()
        );
    }
}
