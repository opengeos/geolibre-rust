//! GeoLibre tool: buffer features with a choice of end-cap and corner-join
//! styles — the pure-Rust counterpart of ArcGIS Pro's *Graphic Buffer*
//! (Analysis).
//!
//! The bundled `buffer_vector` / `multiple_ring_buffer` tools only ever produce
//! *round* caps and *round* joins (a Minkowski sum with a disk). Graphic
//! Buffer's whole reason to exist is the **cartographic and engineering**
//! buffer geometry that round buffers cannot express: squared-off line ends and
//! sharp, right-angle (mitered) or clipped (beveled) corners.
//!
//! - `cap` — how the ends of line/point features are closed:
//!   - `round` — a semicircular end (a full circle for a point).
//!   - `square` — a squared end that projects half the buffer distance past the
//!     endpoint (a `2·distance` square for a point).
//!   - `butt` / `flat` — a flat end flush with the endpoint (a point vanishes,
//!     having no length to cap).
//! - `join` — how the buffer turns a corner:
//!   - `round` — an arc.
//!   - `miter` — the two offset edges extended until they meet at a sharp
//!     corner, unless the corner is sharper than `miter_limit` allows, in which
//!     case it is beveled.
//!   - `bevel` — the corner cut straight across.
//! - `miter_limit` — the classic ratio of miter length to buffer distance
//!   (ArcGIS/JTS default `10`); a corner whose miter would exceed it is beveled
//!   instead. Only consulted when `join = miter`.
//! - `dissolve` — when true, all buffers are unioned into a single dissolved
//!   (multi)polygon; when false (default) one buffer polygon is emitted per
//!   input feature, carrying that feature's attributes.
//!
//! The offsetting, cap/join construction, and self-union are done by `geo`'s
//! `Buffer` (backed by the pure-Rust `i_overlay` mesh offsetter) — no
//! GDAL/GEOS. `miter_limit` is translated to `i_overlay`'s "minimum sharp
//! angle" threshold via `angle = 2·asin(1 / miter_limit)`, the standard
//! relationship between the miter ratio and the corner half-angle.

use std::collections::BTreeMap;

use geo::algorithm::buffer::{Buffer, BufferStyle, LineCap, LineJoin};
use geo::{
    Area, BooleanOps, Coord as GeoCoord, Geometry as GeoGeometry, LineString as GeoLineString,
    MultiLineString as GeoMultiLineString, MultiPoint as GeoMultiPoint,
    MultiPolygon as GeoMultiPolygon, Point as GeoPoint, Polygon as GeoPolygon,
};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Arc smoothness for `round` caps/joins: the segment-length / radius ratio
/// `i_overlay` uses to tessellate arcs. Matches `geo`'s own default and yields a
/// visually smooth circle (~31 segments over a full turn).
const ROUND_ARC_RATIO: f64 = 0.20;

pub struct GraphicBufferTool;

impl Tool for GraphicBufferTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "graphic_buffer",
            display_name: "Graphic Buffer",
            summary: "Buffer features with a choice of end-cap (round/square/butt) and corner-join (round/miter/bevel) styles, with a miter limit — the squared/mitered buffer geometry that round-only buffers cannot produce.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector file path, format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distance",
                    description: "Buffer distance in the layer's CRS units. Must be positive.",
                    required: true,
                },
                ToolParamSpec {
                    name: "cap",
                    description: "End-cap style for line/point ends: 'round' (default), 'square', or 'butt' (alias 'flat').",
                    required: false,
                },
                ToolParamSpec {
                    name: "join",
                    description: "Corner-join style: 'round' (default), 'miter', or 'bevel'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "miter_limit",
                    description: "Ratio of miter length to buffer distance beyond which a mitered corner is beveled (ArcGIS/JTS default 10). Must be >= 1. Only used when join='miter'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "dissolve",
                    description: "When true, union all buffers into one dissolved (multi)polygon feature. Default false (one buffer per input feature, attributes preserved).",
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
        let schema = layer.schema.clone();
        let layer_name = layer.name.clone();
        let layer_crs = layer.crs.clone();
        let input_count = layer.len();

        ctx.progress.info(&format!(
            "graphic_buffer: {} feature(s) at distance {} (cap={}, join={}{})",
            input_count,
            prm.distance,
            prm.cap.as_str(),
            prm.join.as_str(),
            if prm.join == Join::Miter {
                format!(", miter_limit={}", prm.miter_limit)
            } else {
                String::new()
            }
        ));

        let mut out_features: Vec<Feature> = Vec::with_capacity(input_count);
        let mut empty_count = 0usize;
        // Accumulator for the dissolve case.
        let mut dissolved = GeoMultiPolygon::<f64>::new(Vec::new());

        for feature in &layer.features {
            let Some(geom) = feature.geometry.as_ref() else {
                empty_count += 1;
                continue;
            };
            let Some(geo_geom) = to_geo_geometry(geom) else {
                empty_count += 1;
                continue;
            };
            let buffered = geo_geom.buffer_with_style(prm.style());
            if buffered.0.is_empty() {
                // e.g. a point with a butt cap has no area to enclose.
                empty_count += 1;
                continue;
            }
            if prm.dissolve {
                dissolved = dissolved.union(&buffered);
            } else {
                let mut out = feature.clone();
                out.geometry = Some(multipolygon_to_geometry(&buffered));
                out.fid = out_features.len() as u64;
                out_features.push(out);
            }
        }

        if prm.dissolve && !dissolved.0.is_empty() {
            out_features.push(Feature {
                fid: 0,
                geometry: Some(multipolygon_to_geometry(&dissolved)),
                attributes: Vec::new(),
            });
        }

        let has_multipolygon = out_features
            .iter()
            .any(|f| matches!(f.geometry, Some(Geometry::MultiPolygon(_))));

        let total_area: f64 = out_features
            .iter()
            .filter_map(|f| f.geometry.as_ref())
            .filter_map(to_geo_multipolygon)
            .map(|mp| mp.unsigned_area())
            .sum();

        let mut out_layer = Layer::new(layer_name);
        // When dissolved, attributes are dropped (there is no 1:1 source
        // feature), so start from an empty schema; otherwise keep the input's.
        out_layer.schema = if prm.dissolve {
            wbvector::Schema::default()
        } else {
            schema
        };
        out_layer.crs = layer_crs;
        out_layer.features = out_features;
        out_layer.geom_type = Some(if has_multipolygon {
            GeometryType::MultiPolygon
        } else {
            GeometryType::Polygon
        });

        let feature_count = out_layer.len();
        if empty_count > 0 {
            ctx.progress.info(&format!(
                "{empty_count} feature(s) produced no buffer (empty/unbufferable geometry)"
            ));
        }
        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_count".to_string(), json!(input_count));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("empty_count".to_string(), json!(empty_count));
        outputs.insert("distance".to_string(), json!(prm.distance));
        outputs.insert("cap".to_string(), json!(prm.cap.as_str()));
        outputs.insert("join".to_string(), json!(prm.join.as_str()));
        outputs.insert("dissolve".to_string(), json!(prm.dissolve));
        outputs.insert("total_area".to_string(), json!(total_area));
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Cap {
    Round,
    Square,
    Butt,
}

impl Cap {
    fn as_str(self) -> &'static str {
        match self {
            Self::Round => "round",
            Self::Square => "square",
            Self::Butt => "butt",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Join {
    Round,
    Miter,
    Bevel,
}

impl Join {
    fn as_str(self) -> &'static str {
        match self {
            Self::Round => "round",
            Self::Miter => "miter",
            Self::Bevel => "bevel",
        }
    }
}

struct Params {
    distance: f64,
    cap: Cap,
    join: Join,
    miter_limit: f64,
    dissolve: bool,
}

impl Params {
    /// Builds the `geo` buffer style for these parameters. The miter limit
    /// (miter-length / distance ratio) maps to `i_overlay`'s minimum-sharp-angle
    /// threshold `θ_min = 2·asin(1 / miter_limit)`: a corner is beveled once its
    /// interior angle drops below `θ_min`, which is exactly when the classic
    /// miter ratio `1/sin(θ/2)` exceeds `miter_limit`.
    fn style(&self) -> BufferStyle<f64> {
        let cap = match self.cap {
            Cap::Round => LineCap::Round(ROUND_ARC_RATIO),
            Cap::Square => LineCap::Square,
            Cap::Butt => LineCap::Butt,
        };
        let join = match self.join {
            Join::Round => LineJoin::Round(ROUND_ARC_RATIO),
            Join::Bevel => LineJoin::Bevel,
            Join::Miter => LineJoin::Miter(2.0 * (1.0 / self.miter_limit).asin()),
        };
        BufferStyle::new(self.distance)
            .line_cap(cap)
            .line_join(join)
    }
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let distance = parse_optional_f64(args, "distance")?.ok_or_else(|| {
        ToolError::Validation("missing required numeric parameter 'distance'".to_string())
    })?;
    if !(distance > 0.0 && distance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'distance' must be a positive number".to_string(),
        ));
    }
    let cap = match parse_optional_str(args, "cap")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("round") => Cap::Round,
        Some("square") => Cap::Square,
        Some("butt") | Some("flat") => Cap::Butt,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown cap '{other}' (expected round, square, or butt)"
            )))
        }
    };
    let join = match parse_optional_str(args, "join")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("round") => Join::Round,
        Some("miter") => Join::Miter,
        Some("bevel") => Join::Bevel,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown join '{other}' (expected round, miter, or bevel)"
            )))
        }
    };
    let miter_limit = parse_optional_f64(args, "miter_limit")?.unwrap_or(10.0);
    if !(miter_limit >= 1.0 && miter_limit.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'miter_limit' must be a number >= 1".to_string(),
        ));
    }
    let dissolve = parse_optional_bool(args, "dissolve")?.unwrap_or(false);
    Ok(Params {
        distance,
        cap,
        join,
        miter_limit,
        dissolve,
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

/// Parses an optional boolean parameter, accepting a JSON bool or the strings
/// "true"/"false"/"1"/"0" (host UIs often post checkbox values as strings).
fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
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

// ── geo <-> wbvector geometry conversion ────────────────────────────────────

/// Converts a `wbvector` geometry to a `geo` geometry (dropping any Z/M).
/// Returns `None` for empty geometries that cannot be buffered.
fn to_geo_geometry(geom: &Geometry) -> Option<GeoGeometry<f64>> {
    let g = match geom {
        Geometry::Point(c) => GeoGeometry::Point(GeoPoint::new(c.x, c.y)),
        Geometry::MultiPoint(cs) => {
            if cs.is_empty() {
                return None;
            }
            GeoGeometry::MultiPoint(GeoMultiPoint(
                cs.iter().map(|c| GeoPoint::new(c.x, c.y)).collect(),
            ))
        }
        Geometry::LineString(cs) => {
            if cs.len() < 2 {
                return None;
            }
            GeoGeometry::LineString(coords_to_linestring(cs))
        }
        Geometry::MultiLineString(ls) => {
            let parts: Vec<GeoLineString<f64>> = ls
                .iter()
                .filter(|l| l.len() >= 2)
                .map(|l| coords_to_linestring(l))
                .collect();
            if parts.is_empty() {
                return None;
            }
            GeoGeometry::MultiLineString(GeoMultiLineString(parts))
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            if exterior.len() < 3 {
                return None;
            }
            GeoGeometry::Polygon(rings_to_polygon(exterior, interiors))
        }
        Geometry::MultiPolygon(parts) => {
            let polys: Vec<GeoPolygon<f64>> = parts
                .iter()
                .filter(|(e, _)| e.len() >= 3)
                .map(|(e, hs)| rings_to_polygon(e, hs))
                .collect();
            if polys.is_empty() {
                return None;
            }
            GeoGeometry::MultiPolygon(GeoMultiPolygon(polys))
        }
        Geometry::GeometryCollection(gs) => {
            let members: Vec<GeoGeometry<f64>> = gs.iter().filter_map(to_geo_geometry).collect();
            if members.is_empty() {
                return None;
            }
            GeoGeometry::GeometryCollection(geo::GeometryCollection(members))
        }
    };
    Some(g)
}

/// Converts a `wbvector` polygonal geometry to a `geo` `MultiPolygon` (for the
/// area metric); `None` for non-polygonal geometries.
fn to_geo_multipolygon(geom: &Geometry) -> Option<GeoMultiPolygon<f64>> {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => Some(GeoMultiPolygon(vec![rings_to_polygon(exterior, interiors)])),
        Geometry::MultiPolygon(parts) => Some(GeoMultiPolygon(
            parts
                .iter()
                .map(|(e, hs)| rings_to_polygon(e, hs))
                .collect(),
        )),
        _ => None,
    }
}

fn coords_to_linestring(cs: &[Coord]) -> GeoLineString<f64> {
    GeoLineString::new(cs.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect())
}

fn rings_to_polygon(exterior: &Ring, interiors: &[Ring]) -> GeoPolygon<f64> {
    GeoPolygon::new(
        ring_to_linestring(exterior),
        interiors.iter().map(ring_to_linestring).collect(),
    )
}

fn ring_to_linestring(ring: &Ring) -> GeoLineString<f64> {
    // `geo` closes rings itself; `Ring` stores no closing duplicate.
    GeoLineString::new(
        ring.coords()
            .iter()
            .map(|c| GeoCoord { x: c.x, y: c.y })
            .collect(),
    )
}

/// Converts a `geo` `MultiPolygon` to a `wbvector` geometry: a single part
/// becomes a `Polygon`, multiple parts a `MultiPolygon`.
fn multipolygon_to_geometry(mp: &GeoMultiPolygon<f64>) -> Geometry {
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

fn polygon_to_rings(poly: &GeoPolygon<f64>) -> (Ring, Vec<Ring>) {
    (
        linestring_to_ring(poly.exterior()),
        poly.interiors().iter().map(linestring_to_ring).collect(),
    )
}

fn linestring_to_ring(ls: &GeoLineString<f64>) -> Ring {
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
    use wbvector::{memory_store, FieldDef, FieldType};

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
        let out = GraphicBufferTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn area_of(geom: &Geometry) -> f64 {
        to_geo_multipolygon(geom)
            .map(|mp| mp.unsigned_area())
            .unwrap_or(0.0)
    }

    fn point_layer() -> String {
        let mut layer = Layer::new("pts");
        layer.add_field(FieldDef::new("name", FieldType::Text));
        layer
            .add_feature(Some(Geometry::point(0.0, 0.0)), &[("name", "a".into())])
            .unwrap();
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    /// A square-capped point buffer is a `2d × 2d` square (area `4d²`), which is
    /// strictly larger than the round-capped circle (area `πd²`) — the defining
    /// property of Graphic Buffer's square caps.
    #[test]
    fn square_cap_point_encloses_more_than_round() {
        let d = 5.0;
        let (sq_out, sq) =
            run_tool(json!({ "input": point_layer(), "distance": d, "cap": "square" }));
        let (_rd_out, rd) =
            run_tool(json!({ "input": point_layer(), "distance": d, "cap": "round" }));

        let sq_area = area_of(sq.features[0].geometry.as_ref().unwrap());
        let rd_area = area_of(rd.features[0].geometry.as_ref().unwrap());

        // Square cap area is exactly (2d)^2.
        assert!(
            (sq_area - 4.0 * d * d).abs() < 1e-6,
            "square-cap point area {sq_area} != {}",
            4.0 * d * d
        );
        // Round cap approximates a circle: strictly less than the square.
        assert!(
            rd_area < sq_area,
            "round {rd_area} should be < square {sq_area}"
        );
        assert!(rd_area > 2.5 * d * d, "round {rd_area} should be near πd²");
        assert_eq!(sq_out.outputs["cap"], json!("square"));
        // Attributes are preserved (not dissolved).
        assert_eq!(sq.features.len(), 1);
        assert!(sq
            .features
            .iter()
            .all(|f| matches!(f.geometry, Some(Geometry::Polygon { .. }))));
    }

    /// A butt/flat cap on a point has no length to cap, so it produces no buffer.
    #[test]
    fn butt_cap_point_is_empty() {
        let (out, layer) =
            run_tool(json!({ "input": point_layer(), "distance": 3.0, "cap": "butt" }));
        assert_eq!(out.outputs["feature_count"], json!(0));
        assert_eq!(out.outputs["empty_count"], json!(1));
        assert_eq!(layer.features.len(), 0);
    }

    /// Buffering a square polygon: miter joins keep the outward corners sharp
    /// (right angles), so the mitered buffer encloses strictly more area than the
    /// round-join buffer, which rounds every corner off.
    #[test]
    fn miter_join_encloses_more_than_round_and_bevel() {
        let make = || {
            let mut layer = Layer::new("poly");
            layer
                .add_feature(
                    Some(Geometry::polygon(
                        vec![
                            Coord::xy(0.0, 0.0),
                            Coord::xy(10.0, 0.0),
                            Coord::xy(10.0, 10.0),
                            Coord::xy(0.0, 10.0),
                        ],
                        vec![],
                    )),
                    &[],
                )
                .unwrap();
            let id = memory_store::put_vector(layer);
            memory_store::make_vector_memory_path(&id)
        };
        let d = 2.0f64;
        let (_, miter) = run_tool(json!({ "input": make(), "distance": d, "join": "miter" }));
        let (_, bevel) = run_tool(json!({ "input": make(), "distance": d, "join": "bevel" }));
        let (_, round) = run_tool(json!({ "input": make(), "distance": d, "join": "round" }));

        let ma = area_of(miter.features[0].geometry.as_ref().unwrap());
        let ba = area_of(bevel.features[0].geometry.as_ref().unwrap());
        let ra = area_of(round.features[0].geometry.as_ref().unwrap());

        // A miter-buffered axis-aligned square is itself a bigger square:
        // side 10 + 2d, area (10 + 2d)^2.
        let expected_miter = (10.0 + 2.0 * d).powi(2);
        assert!(
            (ma - expected_miter).abs() < 1e-4,
            "miter area {ma} != {expected_miter}"
        );
        // Miter > bevel (corner triangles clipped) > ... and miter > round.
        assert!(ma > ba, "miter {ma} should exceed bevel {ba}");
        assert!(ma > ra, "miter {ma} should exceed round {ra}");
        // Every output geometry is a valid polygon.
        assert!(miter
            .features
            .iter()
            .all(|f| matches!(f.geometry, Some(Geometry::Polygon { .. }))));
    }

    /// A sharp miter is clipped to a bevel once the corner exceeds `miter_limit`;
    /// a generous limit keeps it sharp, so the tighter limit encloses less area.
    #[test]
    fn miter_limit_clips_sharp_corners() {
        // A thin dart with a very sharp tip.
        let make = || {
            let mut layer = Layer::new("dart");
            layer
                .add_feature(
                    Some(Geometry::line_string(vec![
                        Coord::xy(0.0, 0.0),
                        Coord::xy(20.0, 1.0),
                        Coord::xy(0.0, 2.0),
                    ])),
                    &[],
                )
                .unwrap();
            let id = memory_store::put_vector(layer);
            memory_store::make_vector_memory_path(&id)
        };
        let (_, tight) = run_tool(
            json!({ "input": make(), "distance": 1.0, "join": "miter", "miter_limit": 1.5 }),
        );
        let (_, loose) = run_tool(
            json!({ "input": make(), "distance": 1.0, "join": "miter", "miter_limit": 20.0 }),
        );
        let ta = area_of(tight.features[0].geometry.as_ref().unwrap());
        let la = area_of(loose.features[0].geometry.as_ref().unwrap());
        assert!(
            la > ta,
            "loose miter_limit {la} should enclose more than tight {ta} (clipped spike)"
        );
    }

    /// Dissolve unions overlapping buffers into a single feature.
    #[test]
    fn dissolve_unions_into_one_feature() {
        let mut layer = Layer::new("pts");
        // Two points 4 apart; buffered by 3 they overlap and dissolve to one part.
        layer
            .add_feature(Some(Geometry::point(0.0, 0.0)), &[])
            .unwrap();
        layer
            .add_feature(Some(Geometry::point(4.0, 0.0)), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, dissolved) =
            run_tool(json!({ "input": input, "distance": 3.0, "cap": "round", "dissolve": true }));
        assert_eq!(out.outputs["feature_count"], json!(1));
        assert_eq!(dissolved.features.len(), 1);

        // Without dissolve there are two separate buffer features.
        let (out2, _) = run_tool(json!({ "input": input, "distance": 3.0, "cap": "round" }));
        assert_eq!(out2.outputs["feature_count"], json!(2));
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = GraphicBufferTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input must fail");
        assert!(
            bad(json!({ "input": "x.geojson" })).is_err(),
            "missing distance must fail"
        );
        assert!(bad(json!({ "input": "x.geojson", "distance": 0 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "distance": -5 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "distance": 5, "cap": "bogus" })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "distance": 5, "join": "bogus" })).is_err());
        assert!(
            bad(
                json!({ "input": "x.geojson", "distance": 5, "join": "miter", "miter_limit": 0.5 })
            )
            .is_err(),
            "miter_limit < 1 must fail"
        );
        assert!(
            bad(
                json!({ "input": "x.geojson", "distance": "5.0", "cap": "square", "join": "miter" })
            )
            .is_ok(),
            "numeric strings ok"
        );
        assert!(bad(json!({ "input": "x.geojson", "distance": 5, "cap": "flat" })).is_ok());
    }
}
