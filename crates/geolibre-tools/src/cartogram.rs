//! GeoLibre tool: area cartograms (non-contiguous and Dorling).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Cartogram* toolset (Cartography):
//! distort polygon geometry so that each feature's area is proportional to an
//! attribute value — population, GDP, votes — a web-map-friendly visualization
//! whose output feeds straight into `render_vector_png` or `vector_to_pmtiles`.
//! Nothing comparable exists in the bundled suite.
//!
//! Two methods:
//!
//! - `non_contiguous` (default) — scale each polygon about its own centroid by
//!   `√(target_area / current_area)`, where `target_area` distributes the total
//!   map area in proportion to the value. Shapes and topology are preserved
//!   exactly; polygons shrink or grow in place, opening gaps between them. Total
//!   area is conserved.
//! - `dorling` — replace each polygon with a circle whose area is proportional
//!   to the value, placed at the polygon centroid, then push overlapping circles
//!   apart with a few iterations of force-directed relaxation (a light pull back
//!   toward the original location keeps the layout geographically faithful).
//!
//! Every output feature keeps its attributes. Features with a missing or
//! non-positive value are left undistorted and reported.
//!
//! Scope for v1: the contiguous (Gastner–Newman diffusion) cartogram is not
//! implemented — use `non_contiguous` or `dorling`.

use std::collections::BTreeMap;

use geo::{Area, Centroid, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldValue, Geometry, GeometryType, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Vertices used to approximate each Dorling circle.
const CIRCLE_SEGMENTS: usize = 48;

pub struct CartogramTool;

impl Tool for CartogramTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "cartogram",
            display_name: "Cartogram",
            summary: "Distort polygons so area is proportional to an attribute value: non-contiguous (scale each polygon about its centroid) or Dorling (proportional circles with overlap removal).",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon vector layer, format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "value_field",
                    description: "Numeric field whose value each feature's area should be proportional to.",
                    required: true,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'non_contiguous' (default; scale polygons in place) or 'dorling' (proportional circles).",
                    required: false,
                },
                ToolParamSpec {
                    name: "iterations",
                    description: "Dorling overlap-removal iterations (default 100). Ignored by non_contiguous.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["input", "value_field"] {
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
        let input = require_str(args, "input")?;
        let value_field = require_str(args, "value_field")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let schema = layer.schema.clone();

        // Per feature: geo polygon, centroid, area, value (None if unusable).
        struct Node {
            fi: usize,
            centroid: (f64, f64),
            area: f64,
            value: f64,
        }
        let mut nodes: Vec<Node> = Vec::new();
        for (fi, feature) in layer.features.iter().enumerate() {
            let Some(mp) = feature.geometry.as_ref().and_then(to_multipolygon) else {
                continue;
            };
            let value = feature
                .get(&schema, value_field)
                .ok()
                .and_then(FieldValue::as_f64)
                .filter(|v| v.is_finite() && *v > 0.0);
            let area = mp.unsigned_area();
            let centroid = mp.centroid().map(|c| (c.x(), c.y()));
            match (value, centroid) {
                (Some(value), Some(centroid)) if area > 0.0 => nodes.push(Node {
                    fi,
                    centroid,
                    area,
                    value,
                }),
                _ => {} // left undistorted
            }
        }
        if nodes.is_empty() {
            return Err(ToolError::Execution(
                "no polygon features with a positive value to build a cartogram".to_string(),
            ));
        }
        let total_area: f64 = nodes.iter().map(|nd| nd.area).sum();
        let total_value: f64 = nodes.iter().map(|nd| nd.value).sum();

        ctx.progress.info(&format!(
            "cartogram ({}) over {} feature(s)",
            prm.method.as_str(),
            nodes.len()
        ));

        let mut has_multi = false;
        match prm.method {
            Method::NonContiguous => {
                for nd in &nodes {
                    // target_area / area, area ∝ value, total area conserved.
                    let target = nd.value / total_value * total_area;
                    let s = (target / nd.area).sqrt();
                    let geom = scale_geometry(&layer.features[nd.fi].geometry, nd.centroid, s);
                    has_multi |= matches!(geom, Some(Geometry::MultiPolygon(_)));
                    layer.features[nd.fi].geometry = geom;
                }
            }
            Method::Dorling => {
                // Radius so that π r² ∝ value and Σ π r² = total polygon area.
                let k = total_area / (std::f64::consts::PI * total_value);
                let radii: Vec<f64> = nodes.iter().map(|nd| (k * nd.value).sqrt()).collect();
                let origins: Vec<(f64, f64)> = nodes.iter().map(|nd| nd.centroid).collect();
                let positions = dorling_relax(&origins, &radii, prm.iterations);
                for (i, nd) in nodes.iter().enumerate() {
                    layer.features[nd.fi].geometry =
                        Some(circle(positions[i].0, positions[i].1, radii[i]));
                }
            }
        }
        if matches!(prm.method, Method::Dorling) {
            layer.geom_type = Some(GeometryType::Polygon);
        } else if has_multi {
            layer.geom_type = Some(GeometryType::MultiPolygon);
        }

        let feature_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("method".to_string(), json!(prm.method.as_str()));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("cartogram_count".to_string(), json!(nodes.len()));
        outputs.insert("total_value".to_string(), json!(total_value));
        Ok(ToolRunResult { outputs })
    }
}

// ── Non-contiguous scaling ────────────────────────────────────────────────────

/// Scales a geometry's vertices about `(cx, cy)` by `s`.
fn scale_geometry(geom: &Option<Geometry>, c: (f64, f64), s: f64) -> Option<Geometry> {
    let scale_ring = |ring: &Ring| {
        Ring::new(
            ring.coords()
                .iter()
                .map(|p| Coord::xy(c.0 + s * (p.x - c.0), c.1 + s * (p.y - c.1)))
                .collect(),
        )
    };
    match geom.as_ref()? {
        Geometry::Polygon {
            exterior,
            interiors,
        } => Some(Geometry::Polygon {
            exterior: scale_ring(exterior),
            interiors: interiors.iter().map(scale_ring).collect(),
        }),
        Geometry::MultiPolygon(parts) => Some(Geometry::MultiPolygon(
            parts
                .iter()
                .map(|(e, i)| (scale_ring(e), i.iter().map(scale_ring).collect()))
                .collect(),
        )),
        other => Some(other.clone()),
    }
}

// ── Dorling relaxation ────────────────────────────────────────────────────────

/// Pushes overlapping circles apart with a few force-directed iterations, with a
/// light pull back toward each circle's origin to keep the layout faithful.
fn dorling_relax(origins: &[(f64, f64)], radii: &[f64], iterations: usize) -> Vec<(f64, f64)> {
    let n = origins.len();
    let mut pos = origins.to_vec();
    // Attraction toward the original location each step (gentle).
    let attraction = 0.05;
    for _ in 0..iterations {
        let mut disp = vec![(0.0, 0.0); n];
        for i in 0..n {
            for j in i + 1..n {
                let (dx, dy) = (pos[j].0 - pos[i].0, pos[j].1 - pos[i].1);
                let d = (dx * dx + dy * dy).sqrt();
                let min_d = radii[i] + radii[j];
                if d < min_d && d > 1e-12 {
                    // Split the overlap between the two circles.
                    let push = (min_d - d) * 0.5;
                    let (ux, uy) = (dx / d, dy / d);
                    disp[i].0 -= ux * push;
                    disp[i].1 -= uy * push;
                    disp[j].0 += ux * push;
                    disp[j].1 += uy * push;
                } else if d <= 1e-12 {
                    // Coincident centres: nudge deterministically apart.
                    disp[i].0 -= radii[i] * 0.5;
                    disp[j].0 += radii[j] * 0.5;
                }
            }
        }
        for i in 0..n {
            pos[i].0 += disp[i].0 + attraction * (origins[i].0 - pos[i].0);
            pos[i].1 += disp[i].1 + attraction * (origins[i].1 - pos[i].1);
        }
    }
    pos
}

fn circle(cx: f64, cy: f64, r: f64) -> Geometry {
    let coords: Vec<Coord> = (0..CIRCLE_SEGMENTS)
        .map(|i| {
            let a = 2.0 * std::f64::consts::PI * i as f64 / CIRCLE_SEGMENTS as f64;
            Coord::xy(cx + r * a.cos(), cy + r * a.sin())
        })
        .collect();
    Geometry::Polygon {
        exterior: Ring::new(coords),
        interiors: vec![],
    }
}

// ── Conversion ────────────────────────────────────────────────────────────────

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

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    NonContiguous,
    Dorling,
}

impl Method {
    fn as_str(self) -> &'static str {
        match self {
            Self::NonContiguous => "non_contiguous",
            Self::Dorling => "dorling",
        }
    }
}

struct Params {
    method: Method,
    iterations: usize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let method = match parse_optional_str(args, "method")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("non_contiguous") => Method::NonContiguous,
        Some("dorling") => Method::Dorling,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown method '{other}' (expected non_contiguous or dorling)"
            )))
        }
    };
    let iterations = match parse_optional_f64(args, "iterations")? {
        None => 100,
        Some(v) if v.fract() == 0.0 && (0.0..=100_000.0).contains(&v) => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "parameter 'iterations' must be a non-negative integer".to_string(),
            ))
        }
    };
    Ok(Params { method, iterations })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .map(str::trim)
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
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

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn square(x0: f64, y0: f64, s: f64, val: f64) -> (Geometry, FieldValue) {
        (
            Geometry::polygon(
                vec![
                    Coord::xy(x0, y0),
                    Coord::xy(x0 + s, y0),
                    Coord::xy(x0 + s, y0 + s),
                    Coord::xy(x0, y0 + s),
                ],
                vec![],
            ),
            FieldValue::Float(val),
        )
    }

    fn layer_of(items: Vec<(Geometry, FieldValue)>) -> String {
        let mut layer = Layer::new("polys");
        layer.add_field(FieldDef::new("val", FieldType::Float));
        for (g, v) in items {
            layer.add_feature(Some(g), &[("val", v)]).unwrap();
        }
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CartogramTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn area(layer: &Layer, i: usize) -> f64 {
        to_multipolygon(layer.features[i].geometry.as_ref().unwrap())
            .unwrap()
            .unsigned_area()
    }

    #[test]
    fn non_contiguous_makes_area_proportional_to_value() {
        // Two equal 10x10 squares (area 100) with values 1 and 3. After: areas
        // in ratio 1:3, total area (200) preserved.
        let input = layer_of(vec![
            square(0.0, 0.0, 10.0, 1.0),
            square(20.0, 0.0, 10.0, 3.0),
        ]);
        let (_, layer) = run(json!({ "input": input, "value_field": "val" }));
        let (a0, a1) = (area(&layer, 0), area(&layer, 1));
        assert!(
            (a0 + a1 - 200.0).abs() < 1e-6,
            "total area not conserved: {a0}+{a1}"
        );
        assert!(
            (a1 / a0 - 3.0).abs() < 1e-6,
            "area ratio {} should be 3",
            a1 / a0
        );
        // Attributes preserved.
        assert_eq!(
            layer.features[1].get(&layer.schema, "val").unwrap(),
            &FieldValue::Float(3.0)
        );
    }

    #[test]
    fn non_contiguous_scales_about_centroid() {
        // A single square with any value keeps the same centroid (scaled in place).
        let input = layer_of(vec![square(0.0, 0.0, 10.0, 5.0)]);
        let (_, layer) = run(json!({ "input": input, "value_field": "val" }));
        // Only one feature, so target area == its area -> unchanged.
        assert!((area(&layer, 0) - 100.0).abs() < 1e-6);
    }

    #[test]
    fn dorling_radii_scale_with_sqrt_value_and_do_not_overlap() {
        // Values 1 and 4 -> radius ratio 2. Centroids 100 apart so they start
        // separate; relaxation keeps them non-overlapping.
        let input = layer_of(vec![
            square(0.0, 0.0, 2.0, 1.0),
            square(100.0, 0.0, 2.0, 4.0),
        ]);
        let (out, layer) =
            run(json!({ "input": input, "value_field": "val", "method": "dorling" }));
        assert_eq!(out.outputs["method"], json!("dorling"));
        // Radius from area: r = sqrt(area/pi) per circle polygon.
        let r = |i: usize| (area(&layer, i) / std::f64::consts::PI).sqrt();
        let (r0, r1) = (r(0), r(1));
        assert!(
            (r1 / r0 - 2.0).abs() < 0.05,
            "radius ratio {} should be ~2",
            r1 / r0
        );
        // Circles must be output as polygons.
        assert!(matches!(
            layer.features[0].geometry.as_ref().unwrap(),
            Geometry::Polygon { .. }
        ));
    }

    #[test]
    fn dorling_pushes_overlapping_circles_apart() {
        // Three circles seeded at the SAME centroid must separate after relaxation.
        let input = layer_of(vec![
            square(0.0, 0.0, 4.0, 1.0),
            square(0.0, 0.0, 4.0, 1.0),
            square(0.0, 0.0, 4.0, 1.0),
        ]);
        let (_, layer) = run(
            json!({ "input": input, "value_field": "val", "method": "dorling", "iterations": 300 }),
        );
        // Centres (centroids of the output circles) must be pairwise separated.
        let cen = |i: usize| {
            let mp = to_multipolygon(layer.features[i].geometry.as_ref().unwrap()).unwrap();
            let c = mp.centroid().unwrap();
            (c.x(), c.y())
        };
        let r = (area(&layer, 0) / std::f64::consts::PI).sqrt();
        for a in 0..3 {
            for b in a + 1..3 {
                let (ca, cb) = (cen(a), cen(b));
                let d = ((ca.0 - cb.0).powi(2) + (ca.1 - cb.1).powi(2)).sqrt();
                assert!(
                    d > r,
                    "circles {a},{b} still overlap heavily (d={d}, r={r})"
                );
            }
        }
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = CartogramTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "input": "x.geojson" })).is_err(),
            "missing value_field"
        );
        assert!(
            bad(json!({ "input": "x.geojson", "value_field": "v", "method": "bogus" })).is_err()
        );
        assert!(
            bad(json!({ "input": "x.geojson", "value_field": "v", "iterations": -1 })).is_err()
        );
        assert!(bad(json!({ "input": "x.geojson", "value_field": "v" })).is_ok());
    }
}
