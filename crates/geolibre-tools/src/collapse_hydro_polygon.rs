//! GeoLibre tool: collapse narrow water polygons to centerlines.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Collapse Hydro Polygon*
//! (Cartography). A river network digitized as polygons becomes unreadable when
//! zoomed out; the cartographic fix is to collapse narrow reaches to a
//! centerline while keeping wide reaches as polygons. The near-misses don't
//! cover it: the GeoLibre `collapse_dual_lines_to_centerline` collapses paired
//! *line* casings, and the bundled `river_centerlines` is raster-based and
//! all-or-nothing (no width threshold). Extends the generalization identity
//! (`thin_road_network`, `aggregate_polygons`, `simplify_shared_edges`).
//!
//! Each water polygon's centerline is traced by walking its principal axis (PCA
//! of the boundary) and, at stations spaced `sample_distance` apart, taking the
//! midpoint of the polygon's perpendicular cross-section — which also measures
//! the local width. A polygon whose median width is at or below `collapse_width`
//! is emitted as a centerline (`output`); wider polygons are passed to an
//! optional `retained` polygon layer. Centerlines shorter than `min_length` are
//! dropped.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{
    Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring,
};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CollapseHydroPolygonTool;

impl Tool for CollapseHydroPolygonTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "collapse_hydro_polygon",
            display_name: "Collapse Hydro Polygon",
            summary: "Collapse narrow water polygons to centerlines while keeping wide reaches as polygons (like ArcGIS Collapse Hydro Polygon): trace each polygon's centerline along its principal axis via perpendicular cross-section midpoints, measure local width, and collapse where the median width is at or below a threshold. The width-thresholded polygon→line generalization the line-casing collapse_dual_lines_to_centerline and raster river_centerlines don't do.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input water polygon layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output centerline layer for collapsed (narrow) polygons. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "collapse_width",
                    description: "Polygons whose median width is at or below this collapse to a centerline.",
                    required: true,
                },
                ToolParamSpec {
                    name: "sample_distance",
                    description: "Spacing of cross-section stations along the axis (default: collapse_width / 2).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_length",
                    description: "Drop collapsed centerlines shorter than this (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "retained",
                    description: "Optional output polygon layer for wide (non-collapsed) reaches.",
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
        let retained_path = parse_optional_str(args, "retained")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;

        let mut lines = Layer::new("centerlines").with_geom_type(GeometryType::LineString);
        let mut retained = Layer::new("retained").with_geom_type(GeometryType::Polygon);
        if let Some(e) = layer.crs_epsg() {
            lines = lines.with_crs_epsg(e);
            retained = retained.with_crs_epsg(e);
        }
        lines.add_field(FieldDef::new("median_w", FieldType::Float));
        lines.add_field(FieldDef::new("length", FieldType::Float));
        retained.add_field(FieldDef::new("median_w", FieldType::Float));

        let mut collapsed = 0usize;
        let mut kept = 0usize;
        for feat in &layer.features {
            let Some(geom) = feat.geometry.as_ref() else {
                continue;
            };
            for ring in exterior_rings(geom) {
                if ring.len() < 3 {
                    continue;
                }
                let Some((center, widths)) = centerline(&ring, prm.sample_distance) else {
                    continue;
                };
                let median_w = median(&mut widths.clone());
                if median_w <= prm.collapse_width && center.len() >= 2 {
                    let length = polyline_length(&center);
                    if length < prm.min_length {
                        continue;
                    }
                    lines.push(Feature {
                        fid: 0,
                        geometry: Some(Geometry::line_string(center)),
                        attributes: vec![FieldValue::Float(median_w), FieldValue::Float(length)],
                    });
                    collapsed += 1;
                } else if retained_path.is_some() {
                    retained.push(Feature {
                        fid: 0,
                        geometry: Some(Geometry::polygon(ring.coords().to_vec(), Vec::new())),
                        attributes: vec![FieldValue::Float(median_w)],
                    });
                    kept += 1;
                } else {
                    kept += 1;
                }
            }
        }

        ctx.progress
            .info(&format!("{collapsed} collapsed, {kept} retained"));

        let out_path = write_or_store_layer(lines, output)?;
        let retained_out = match retained_path {
            Some(p) => Some(write_or_store_layer(retained, Some(p))?),
            None => None,
        };

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("collapsed".to_string(), json!(collapsed));
        outputs.insert("retained".to_string(), json!(kept));
        if let Some(r) = retained_out {
            outputs.insert("retained_output".to_string(), json!(r));
        }
        Ok(ToolRunResult { outputs })
    }
}

/// Traces a polygon ring's centerline along its principal axis, returning the
/// centerline points and the per-station widths.
fn centerline(ring: &Ring, sample_distance: f64) -> Option<(Vec<Coord>, Vec<f64>)> {
    let pts: Vec<(f64, f64)> = ring.coords().iter().map(|c| (c.x, c.y)).collect();
    let n = pts.len();
    if n < 3 {
        return None;
    }
    // Principal axis via PCA of the vertices.
    let (cx, cy) = (
        pts.iter().map(|p| p.0).sum::<f64>() / n as f64,
        pts.iter().map(|p| p.1).sum::<f64>() / n as f64,
    );
    let (mut sxx, mut syy, mut sxy) = (0.0, 0.0, 0.0);
    for &(x, y) in &pts {
        let (dx, dy) = (x - cx, y - cy);
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    let angle = 0.5 * (2.0 * sxy).atan2(sxx - syy);
    let (ux, uy) = (angle.cos(), angle.sin()); // axis direction
    let (vx, vy) = (-uy, ux); // perpendicular

    // Extent along the axis.
    let ts: Vec<f64> = pts
        .iter()
        .map(|&(x, y)| (x - cx) * ux + (y - cy) * uy)
        .collect();
    let tmin = ts.iter().cloned().fold(f64::INFINITY, f64::min);
    let tmax = ts.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let axis_len = tmax - tmin;
    if axis_len <= 0.0 {
        return None;
    }
    let step = if sample_distance > 0.0 {
        sample_distance
    } else {
        (axis_len / 20.0).max(1e-9)
    };
    let n_stations = ((axis_len / step).floor() as usize).max(2);

    // Polygon edges as (a, b).
    let edges: Vec<((f64, f64), (f64, f64))> = (0..n).map(|i| (pts[i], pts[(i + 1) % n])).collect();

    let mut center = Vec::new();
    let mut widths = Vec::new();
    for k in 0..=n_stations {
        let t = tmin + axis_len * (k as f64 / n_stations as f64);
        // Station point on the axis, slightly inset from the ends.
        let px = cx + t * ux;
        let py = cy + t * uy;
        // Intersect the perpendicular line (through P, dir v) with all edges;
        // collect signed offsets along v.
        let mut offs: Vec<f64> = Vec::new();
        for &(a, b) in &edges {
            if let Some(off) = perp_intersection(px, py, vx, vy, a, b) {
                offs.push(off);
            }
        }
        if offs.len() < 2 {
            continue;
        }
        offs.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let lo = offs[0];
        let hi = offs[offs.len() - 1];
        let mid = (lo + hi) / 2.0;
        center.push(Coord::xy(px + mid * vx, py + mid * vy));
        widths.push(hi - lo);
    }
    if center.len() < 2 {
        return None;
    }
    Some((center, widths))
}

/// Offset along `v` from `(px,py)` where the infinite line hits segment `a-b`,
/// or `None`. The line is `P + s·v`; the segment is `a + r·(b-a)`, r ∈ [0,1].
fn perp_intersection(
    px: f64,
    py: f64,
    vx: f64,
    vy: f64,
    a: (f64, f64),
    b: (f64, f64),
) -> Option<f64> {
    let (ex, ey) = (b.0 - a.0, b.1 - a.1);
    let denom = vx * (-ey) - vy * (-ex); // det([v, -e])
    if denom.abs() < 1e-15 {
        return None;
    }
    let (wx, wy) = (a.0 - px, a.1 - py);
    // Solve [v, -e] · [s, r]^T = w.
    let s = (wx * (-ey) - wy * (-ex)) / denom;
    let r = (vx * wy - vy * wx) / denom;
    if (0.0..=1.0).contains(&r) {
        Some(s)
    } else {
        None
    }
}

fn polyline_length(pts: &[Coord]) -> f64 {
    pts.windows(2)
        .map(|w| (w[1].x - w[0].x).hypot(w[1].y - w[0].y))
        .sum()
}

fn median(v: &mut [f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m = v.len() / 2;
    if v.len() % 2 == 1 {
        v[m]
    } else {
        (v[m - 1] + v[m]) / 2.0
    }
}

fn exterior_rings(geom: &Geometry) -> Vec<Ring> {
    match geom {
        Geometry::Polygon { exterior, .. } => vec![exterior.clone()],
        Geometry::MultiPolygon(parts) => parts.iter().map(|(e, _)| e.clone()).collect(),
        _ => Vec::new(),
    }
}

struct Params {
    collapse_width: f64,
    sample_distance: f64,
    min_length: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let collapse_width = match args.get("collapse_width") {
        Some(Value::Number(n)) => n.as_f64().unwrap_or(0.0),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'collapse_width' must be a number".into()))?,
        _ => {
            return Err(ToolError::Validation(
                "required parameter 'collapse_width' is missing".into(),
            ))
        }
    };
    if collapse_width <= 0.0 || collapse_width.is_nan() {
        return Err(ToolError::Validation(
            "'collapse_width' must be positive".into(),
        ));
    }
    let sample_distance = match args.get("sample_distance") {
        None | Some(Value::Null) => collapse_width / 2.0,
        Some(Value::Number(n)) => n
            .as_f64()
            .filter(|v| *v > 0.0)
            .unwrap_or(collapse_width / 2.0),
        Some(Value::String(s)) if s.trim().is_empty() => collapse_width / 2.0,
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'sample_distance' must be a number".into()))?,
        _ => collapse_width / 2.0,
    };
    let min_length = match args.get("min_length") {
        None | Some(Value::Null) => 0.0,
        Some(Value::Number(n)) => n.as_f64().unwrap_or(0.0).max(0.0),
        Some(Value::String(s)) if s.trim().is_empty() => 0.0,
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'min_length' must be a number".into()))?
            .max(0.0),
        _ => 0.0,
    };
    Ok(Params {
        collapse_width,
        sample_distance,
        min_length,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::memory_store;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn poly_layer(rings: &[Vec<(f64, f64)>]) -> String {
        let mut l = Layer::new("water")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("id", FieldType::Integer));
        for (i, r) in rings.iter().enumerate() {
            let coords: Vec<Coord> = r.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
            l.add_feature(
                Some(Geometry::polygon(coords, Vec::new())),
                &[("id", (i as i64).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CollapseHydroPolygonTool.run(&args, &ctx()).unwrap();
        let l = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, l)
    }

    /// A long thin river ribbon (width 2, length 100) collapses to a centerline.
    #[test]
    fn collapses_narrow_ribbon() {
        // Rectangle from (0,0) to (100, 2): width 2.
        let ribbon = vec![(0.0, 0.0), (100.0, 0.0), (100.0, 2.0), (0.0, 2.0)];
        let (out, l) = run(json!({
            "input": poly_layer(&[ribbon]), "collapse_width": 5.0, "sample_distance": 5.0,
        }));
        assert_eq!(
            out.outputs["collapsed"],
            json!(1),
            "narrow ribbon should collapse"
        );
        assert_eq!(l.features.len(), 1);
        // The centerline runs the length of the ribbon near y=1.
        if let Some(Geometry::LineString(cs)) = &l.features[0].geometry {
            let length: f64 = cs
                .windows(2)
                .map(|w| (w[1].x - w[0].x).hypot(w[1].y - w[0].y))
                .sum();
            assert!(
                length > 80.0,
                "centerline should span the ribbon, got {length}"
            );
            assert!(
                cs.iter().all(|c| (c.y - 1.0).abs() < 0.6),
                "centerline near the ribbon centre"
            );
        } else {
            panic!("expected a LineString");
        }
    }

    /// A wide reach (width 40) is retained as a polygon, not collapsed.
    #[test]
    fn retains_wide_reach() {
        let lake = vec![(0.0, 0.0), (100.0, 0.0), (100.0, 40.0), (0.0, 40.0)];
        let (out, _l) = run(json!({
            "input": poly_layer(&[lake]), "collapse_width": 5.0, "sample_distance": 5.0,
        }));
        assert_eq!(
            out.outputs["collapsed"],
            json!(0),
            "wide reach must not collapse"
        );
        assert_eq!(out.outputs["retained"], json!(1));
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CollapseHydroPolygonTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "collapse_width": 0 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "collapse_width": 5 })).is_ok());
    }
}
