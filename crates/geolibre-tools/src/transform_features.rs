//! GeoLibre tool: apply a plain affine transform to every geometry in a vector
//! layer.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Rotate*, *Mirror*, *Shift* and
//! *Rescale* (Data Management) — bundled here into a single, composable affine
//! step. There is no equivalent whitebox tool: `affine` exists in WhiteboxTools
//! only as a fitted model inside raster GCP warping, and GeoLibre's
//! `warp_raster` / `rubbersheet_features` / `align_features` are raster-warp or
//! control-point conflation rather than a deterministic geometric transform.
//!
//! The transform is composed, about a chosen `anchor` point `(ax, ay)`, as:
//!
//! ```text
//!   p' = R · M · S · (p - anchor) + anchor + (dx, dy)
//! ```
//!
//! where `S` is the per-axis scale `diag(scale_x, scale_y)`, `M` is the optional
//! mirror (reflection across the X or Y axis through the anchor), `R` is a
//! counter-clockwise rotation by `angle` degrees, and `(dx, dy)` is the final
//! shift. Scale, mirror and rotation therefore all pivot about the anchor, so a
//! layer rotated or rescaled about its `CENTROID` stays put; the shift then
//! slides the whole result.
//!
//! `anchor` is one of:
//! - `CENTROID` (default) — the center of the layer's overall bounding box, so
//!   the layer is transformed in place.
//! - `ORIGIN` — the coordinate origin `(0, 0)`.
//! - `XY` — an explicit pivot given by `anchor_x` / `anchor_y`.
//!
//! Every OGC geometry type is handled (points, lines, polygons and their multi-
//! and collection forms), and per-vertex `Z`/`M` values pass through untouched —
//! this is a planar XY transform. Attributes, feature order and the layer schema
//! are preserved exactly; only coordinates move. The transform is affine, so a
//! rotation- or mirror-only transform conserves area and a scale multiplies
//! every area by `|scale_x · scale_y|`, which the tool reports back as
//! `area_scale` for validation.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Geometry, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct TransformFeaturesTool;

impl Tool for TransformFeaturesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "transform_features",
            display_name: "Transform Features",
            summary: "Apply an affine transform (shift, rotate about an anchor, per-axis scale, mirror) to every geometry in a vector layer.",
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
                    name: "dx",
                    description: "Shift applied to X after scale/mirror/rotation, in CRS units. Default 0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "dy",
                    description: "Shift applied to Y after scale/mirror/rotation, in CRS units. Default 0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "angle",
                    description: "Counter-clockwise rotation about the anchor, in degrees. Default 0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "scale_x",
                    description: "Scale factor about the anchor along X. Default 1.",
                    required: false,
                },
                ToolParamSpec {
                    name: "scale_y",
                    description: "Scale factor about the anchor along Y. Default 1.",
                    required: false,
                },
                ToolParamSpec {
                    name: "mirror_axis",
                    description: "Reflect through the anchor: NONE (default), X (flip Y across the horizontal axis), or Y (flip X across the vertical axis).",
                    required: false,
                },
                ToolParamSpec {
                    name: "anchor",
                    description: "Pivot for scale/mirror/rotation: CENTROID (bounding-box center, default), ORIGIN (0,0), or XY (use anchor_x/anchor_y).",
                    required: false,
                },
                ToolParamSpec {
                    name: "anchor_x",
                    description: "Anchor X coordinate; required when anchor=XY.",
                    required: false,
                },
                ToolParamSpec {
                    name: "anchor_y",
                    description: "Anchor Y coordinate; required when anchor=XY.",
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

        let mut layer = load_input_layer(input)?;
        let input_count = layer.len();

        // Resolve the anchor. CENTROID needs the layer's overall extent, so
        // compute it once up front from every vertex.
        let anchor = match prm.anchor {
            Anchor::Origin => (0.0, 0.0),
            Anchor::Xy(x, y) => (x, y),
            Anchor::Centroid => layer_extent_center(&layer).ok_or_else(|| {
                ToolError::Execution(
                    "anchor=CENTROID requires at least one coordinate in the layer".to_string(),
                )
            })?,
        };

        let affine = Affine::compose(&prm, anchor);

        let mut moved = 0usize;
        for feature in &mut layer.features {
            if let Some(geom) = feature.geometry.as_mut() {
                moved += transform_geometry(geom, &affine);
            }
        }

        // An affine transform maps geometry types to themselves, so the declared
        // layer geom_type is still valid and is left as-is.
        ctx.progress.info(&format!(
            "transformed {input_count} feature(s); moved {moved} vertices"
        ));

        let feature_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_count".to_string(), json!(input_count));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("vertices_moved".to_string(), json!(moved));
        outputs.insert("anchor_x".to_string(), json!(anchor.0));
        outputs.insert("anchor_y".to_string(), json!(anchor.1));
        outputs.insert(
            "area_scale".to_string(),
            json!((prm.scale_x * prm.scale_y).abs()),
        );
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Mirror {
    None,
    X,
    Y,
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum Anchor {
    Centroid,
    Origin,
    Xy(f64, f64),
}

struct Params {
    dx: f64,
    dy: f64,
    angle_deg: f64,
    scale_x: f64,
    scale_y: f64,
    mirror: Mirror,
    anchor: Anchor,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let dx = parse_optional_f64(args, "dx")?.unwrap_or(0.0);
    let dy = parse_optional_f64(args, "dy")?.unwrap_or(0.0);
    let angle_deg = parse_optional_f64(args, "angle")?.unwrap_or(0.0);
    let scale_x = parse_optional_f64(args, "scale_x")?.unwrap_or(1.0);
    let scale_y = parse_optional_f64(args, "scale_y")?.unwrap_or(1.0);
    for (name, v) in [("dx", dx), ("dy", dy), ("angle", angle_deg)] {
        if !v.is_finite() {
            return Err(ToolError::Validation(format!(
                "parameter '{name}' must be a finite number"
            )));
        }
    }
    for (name, v) in [("scale_x", scale_x), ("scale_y", scale_y)] {
        if !(v.is_finite() && v != 0.0) {
            return Err(ToolError::Validation(format!(
                "parameter '{name}' must be a non-zero finite number"
            )));
        }
    }

    let mirror = match parse_optional_str(args, "mirror_axis")?
        .map(|s| s.trim().to_ascii_uppercase())
        .as_deref()
    {
        None | Some("NONE") => Mirror::None,
        Some("X") => Mirror::X,
        Some("Y") => Mirror::Y,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown mirror_axis '{other}' (expected NONE, X or Y)"
            )))
        }
    };

    let anchor = match parse_optional_str(args, "anchor")?
        .map(|s| s.trim().to_ascii_uppercase())
        .as_deref()
    {
        None | Some("CENTROID") => Anchor::Centroid,
        Some("ORIGIN") => Anchor::Origin,
        Some("XY") => {
            let x = parse_optional_f64(args, "anchor_x")?.ok_or_else(|| {
                ToolError::Validation("anchor=XY requires 'anchor_x'".to_string())
            })?;
            let y = parse_optional_f64(args, "anchor_y")?.ok_or_else(|| {
                ToolError::Validation("anchor=XY requires 'anchor_y'".to_string())
            })?;
            if !(x.is_finite() && y.is_finite()) {
                return Err(ToolError::Validation(
                    "parameters 'anchor_x'/'anchor_y' must be finite numbers".to_string(),
                ));
            }
            Anchor::Xy(x, y)
        }
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown anchor '{other}' (expected CENTROID, ORIGIN or XY)"
            )))
        }
    };

    Ok(Params {
        dx,
        dy,
        angle_deg,
        scale_x,
        scale_y,
        mirror,
        anchor,
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

// ── Affine transform ────────────────────────────────────────────────────────

/// A 2-D affine map `x' = a·x + b·y + e`, `y' = c·x + d·y + f`, applied to the
/// XY of every coordinate. Z and M pass through unchanged.
struct Affine {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    f: f64,
}

impl Affine {
    /// Builds the composed map `R · M · S · (p - anchor) + anchor + (dx, dy)`.
    fn compose(prm: &Params, anchor: (f64, f64)) -> Self {
        // Linear part before the anchor translation: rotation ∘ mirror ∘ scale.
        // Scale: diag(sx, sy). Mirror flips one axis' sign. Rotation is the 2×2
        // CCW matrix. We multiply them out into a single 2×2 (la..ld).
        let (sx, sy) = (prm.scale_x, prm.scale_y);
        // Mirror X reflects across the X-axis => y -> -y (negate the Y row);
        // Mirror Y reflects across the Y-axis => x -> -x (negate the X row).
        let (mx, my) = match prm.mirror {
            Mirror::None => (1.0, 1.0),
            Mirror::X => (1.0, -1.0),
            Mirror::Y => (-1.0, 1.0),
        };
        // Combined scale+mirror (still diagonal): diag(sx·mx, sy·my).
        let (px, py) = (sx * mx, sy * my);
        let theta = prm.angle_deg.to_radians();
        let (cs, sn) = (theta.cos(), theta.sin());
        // Rotation ∘ diag(px, py):
        //   [cs -sn] [px  0 ]   [cs·px  -sn·py]
        //   [sn  cs] [ 0  py] = [sn·px   cs·py]
        let a = cs * px;
        let b = -sn * py;
        let c = sn * px;
        let d = cs * py;
        // Translation so the linear part pivots about the anchor, then + shift:
        //   p' = L·(p - anchor) + anchor + shift
        //      = L·p + (anchor - L·anchor + shift)
        let (ax, ay) = anchor;
        let e = ax - (a * ax + b * ay) + prm.dx;
        let f = ay - (c * ax + d * ay) + prm.dy;
        Self { a, b, c, d, e, f }
    }

    #[inline]
    fn apply(&self, coord: &mut Coord) {
        let (x, y) = (coord.x, coord.y);
        coord.x = self.a * x + self.b * y + self.e;
        coord.y = self.c * x + self.d * y + self.f;
    }
}

/// Applies `affine` to every coordinate in `geom` in place, returning the number
/// of coordinates moved.
fn transform_geometry(geom: &mut Geometry, affine: &Affine) -> usize {
    match geom {
        Geometry::Point(c) => {
            affine.apply(c);
            1
        }
        Geometry::LineString(cs) | Geometry::MultiPoint(cs) => transform_coords(cs, affine),
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            let mut n = transform_ring(exterior, affine);
            for r in interiors.iter_mut() {
                n += transform_ring(r, affine);
            }
            n
        }
        Geometry::MultiLineString(lines) => {
            lines.iter_mut().map(|l| transform_coords(l, affine)).sum()
        }
        Geometry::MultiPolygon(parts) => {
            let mut n = 0;
            for (ext, ints) in parts.iter_mut() {
                n += transform_ring(ext, affine);
                for r in ints.iter_mut() {
                    n += transform_ring(r, affine);
                }
            }
            n
        }
        Geometry::GeometryCollection(gs) => {
            gs.iter_mut().map(|g| transform_geometry(g, affine)).sum()
        }
    }
}

fn transform_coords(cs: &mut [Coord], affine: &Affine) -> usize {
    for c in cs.iter_mut() {
        affine.apply(c);
    }
    cs.len()
}

fn transform_ring(ring: &mut Ring, affine: &Affine) -> usize {
    transform_coords(&mut ring.0, affine)
}

/// Center of the layer's overall bounding box, or `None` if it has no vertices.
fn layer_extent_center(layer: &wbvector::Layer) -> Option<(f64, f64)> {
    let mut minx = f64::INFINITY;
    let mut miny = f64::INFINITY;
    let mut maxx = f64::NEG_INFINITY;
    let mut maxy = f64::NEG_INFINITY;
    let mut any = false;
    for feature in &layer.features {
        if let Some(geom) = feature.geometry.as_ref() {
            for c in geom.all_coords() {
                any = true;
                minx = minx.min(c.x);
                miny = miny.min(c.y);
                maxx = maxx.max(c.x);
                maxy = maxy.max(c.y);
            }
        }
    }
    any.then_some(((minx + maxx) * 0.5, (miny + maxy) * 0.5))
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

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Coord> {
        vec![
            Coord::xy(x0, y0),
            Coord::xy(x1, y0),
            Coord::xy(x1, y1),
            Coord::xy(x0, y1),
        ]
    }

    /// Unsigned area of a single-ring polygon geometry via the shoelace formula.
    fn poly_area(geom: &Geometry) -> f64 {
        match geom {
            Geometry::Polygon { exterior, .. } => exterior.signed_area().abs(),
            _ => 0.0,
        }
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = TransformFeaturesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn one_poly_layer(coords: Vec<Coord>) -> String {
        let mut layer = Layer::new("shapes");
        layer.add_field(FieldDef::new("name", FieldType::Text));
        layer
            .add_feature(
                Some(Geometry::polygon(coords, vec![])),
                &[("name", "a".into())],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    /// A pure shift moves every vertex by (dx, dy) and preserves area exactly.
    #[test]
    fn shift_translates_and_preserves_area() {
        let input = one_poly_layer(rect(0.0, 0.0, 10.0, 20.0));
        let (out, layer) = run_tool(json!({ "input": input, "dx": 5.0, "dy": -3.0 }));
        assert_eq!(out.outputs["area_scale"], json!(1.0));
        let g = layer.features[0].geometry.as_ref().unwrap();
        assert!((poly_area(g) - 200.0).abs() < 1e-9);
        // Lower-left corner (0,0) -> (5,-3).
        if let Geometry::Polygon { exterior, .. } = g {
            assert!((exterior.0[0].x - 5.0).abs() < 1e-9);
            assert!((exterior.0[0].y + 3.0).abs() < 1e-9);
        } else {
            panic!("expected polygon");
        }
    }

    /// Rotating 90° about the extent centroid conserves area and returns the
    /// layer to (nearly) the same extent (a square maps onto itself).
    #[test]
    fn rotate_about_centroid_conserves_area_and_extent() {
        let input = one_poly_layer(rect(0.0, 0.0, 10.0, 10.0));
        let (_, layer) = run_tool(json!({ "input": input, "angle": 90.0 }));
        let g = layer.features[0].geometry.as_ref().unwrap();
        assert!((poly_area(g) - 100.0).abs() < 1e-6);
        // A square rotated 90° about its center keeps the same bounding box.
        let bb = g.bbox().unwrap();
        assert!((bb.min_x - 0.0).abs() < 1e-6 && (bb.max_x - 10.0).abs() < 1e-6);
        assert!((bb.min_y - 0.0).abs() < 1e-6 && (bb.max_y - 10.0).abs() < 1e-6);
    }

    /// Per-axis scale multiplies area by |scale_x·scale_y|; about ORIGIN the
    /// corner at the origin stays fixed.
    #[test]
    fn scale_multiplies_area_about_origin() {
        let input = one_poly_layer(rect(0.0, 0.0, 4.0, 5.0)); // area 20
        let (out, layer) =
            run_tool(json!({ "input": input, "scale_x": 2.0, "scale_y": 3.0, "anchor": "ORIGIN" }));
        assert_eq!(out.outputs["area_scale"], json!(6.0));
        let g = layer.features[0].geometry.as_ref().unwrap();
        assert!((poly_area(g) - 120.0).abs() < 1e-9);
        let bb = g.bbox().unwrap();
        // Origin corner fixed; opposite corner at (8, 15).
        assert!((bb.min_x).abs() < 1e-9 && (bb.min_y).abs() < 1e-9);
        assert!((bb.max_x - 8.0).abs() < 1e-9 && (bb.max_y - 15.0).abs() < 1e-9);
    }

    /// Mirror across the X axis about ORIGIN negates Y and preserves area.
    #[test]
    fn mirror_x_reflects_y() {
        let input = one_poly_layer(rect(1.0, 2.0, 3.0, 6.0)); // area 8
        let (out, layer) =
            run_tool(json!({ "input": input, "mirror_axis": "X", "anchor": "ORIGIN" }));
        assert_eq!(out.outputs["area_scale"], json!(1.0));
        let g = layer.features[0].geometry.as_ref().unwrap();
        assert!((poly_area(g) - 8.0).abs() < 1e-9);
        let bb = g.bbox().unwrap();
        // Y range 2..6 -> -6..-2, X unchanged.
        assert!((bb.min_y + 6.0).abs() < 1e-9 && (bb.max_y + 2.0).abs() < 1e-9);
        assert!((bb.min_x - 1.0).abs() < 1e-9 && (bb.max_x - 3.0).abs() < 1e-9);
    }

    /// Non-polygon geometries (points, lines) are transformed too; Z/M survive.
    #[test]
    fn transforms_points_and_preserves_z() {
        let mut layer = Layer::new("mixed");
        layer
            .add_feature(Some(Geometry::Point(Coord::xyz(2.0, 3.0, 42.0))), &[])
            .unwrap();
        layer
            .add_feature(
                Some(Geometry::line_string(vec![
                    Coord::xy(0.0, 0.0),
                    Coord::xy(1.0, 1.0),
                ])),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) =
            run_tool(json!({ "input": input, "dx": 10.0, "dy": 20.0, "anchor": "ORIGIN" }));
        assert_eq!(out.outputs["vertices_moved"], json!(3));
        if let Geometry::Point(c) = layer.features[0].geometry.as_ref().unwrap() {
            assert!((c.x - 12.0).abs() < 1e-9 && (c.y - 23.0).abs() < 1e-9);
            assert_eq!(c.z, Some(42.0)); // Z untouched by the planar transform.
        } else {
            panic!("expected point");
        }
    }

    /// Identity transform (all defaults) leaves geometry unchanged.
    #[test]
    fn identity_is_a_no_op() {
        let input = one_poly_layer(rect(0.0, 0.0, 10.0, 10.0));
        let (out, layer) = run_tool(json!({ "input": input }));
        assert_eq!(out.outputs["area_scale"], json!(1.0));
        let g = layer.features[0].geometry.as_ref().unwrap();
        let bb = g.bbox().unwrap();
        assert!((bb.min_x).abs() < 1e-12 && (bb.max_x - 10.0).abs() < 1e-12);
        assert!((bb.min_y).abs() < 1e-12 && (bb.max_y - 10.0).abs() < 1e-12);
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = TransformFeaturesTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input must fail");
        assert!(
            bad(json!({ "input": "x.geojson", "scale_x": 0 })).is_err(),
            "zero scale must fail"
        );
        assert!(
            bad(json!({ "input": "x.geojson", "mirror_axis": "Z" })).is_err(),
            "bad mirror axis"
        );
        assert!(
            bad(json!({ "input": "x.geojson", "anchor": "XY" })).is_err(),
            "anchor XY without coords must fail"
        );
        assert!(bad(json!({ "input": "x.geojson", "anchor": "bogus" })).is_err());
        assert!(
            bad(json!({ "input": "x.geojson", "angle": "45" })).is_ok(),
            "numeric strings ok"
        );
        assert!(
            bad(json!({ "input": "x.geojson", "anchor": "XY", "anchor_x": 1, "anchor_y": 2 }))
                .is_ok()
        );
    }
}
