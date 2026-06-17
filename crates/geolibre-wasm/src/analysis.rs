//! Geometry + geostatistics helpers backed by the pure-Rust `wbtopology` and
//! `wbspatialstats` engines. Inputs/outputs are flat coordinate/value arrays
//! for easy interop with JavaScript typed arrays.
use wasm_bindgen::prelude::*;
use wbtopology::geom::{Coord, Geometry};
use wbspatialstats::{SpatialWeightsDiagnostics, SpatialWeightsGraph};

/// Convex hull of a 2D point set. Input is `[x0,y0,x1,y1,...]`; output is the
/// hull ring as `[x0,y0,...]` (closed). Needs at least 3 points.
#[wasm_bindgen]
pub fn convex_hull(points_xy: &[f64]) -> Result<Vec<f64>, JsValue> {
    if points_xy.len() % 2 != 0 || points_xy.len() < 6 {
        return Err(JsValue::from_str("need >= 3 points as a flat [x0,y0,x1,y1,...] array"));
    }
    let coords: Vec<Coord> = points_xy
        .chunks_exact(2)
        .map(|c| Coord { x: c[0], y: c[1], z: None })
        .collect();
    let ring = match wbtopology::convex_hull(&coords, 0.0) {
        Geometry::Polygon(p) => p.exterior.coords,
        Geometry::LineString(l) => l.coords,
        Geometry::Point(c) => vec![c],
        Geometry::MultiPoint(v) => v,
        _ => return Err(JsValue::from_str("unexpected hull geometry")),
    };
    let mut out = Vec::with_capacity(ring.len() * 2);
    for c in ring { out.push(c.x); out.push(c.y); }
    Ok(out)
}

/// Global Moran's I spatial autocorrelation for point data, using a binary
/// distance-band spatial weights matrix (neighbors within `distance_threshold`).
///
/// `points_xy` is `[x0,y0,...]`, `values` is one value per point. Returns JSON:
/// `{"ok":true,"morans_i","expected","variance","z_score","p_value","n"}`.
///
/// Builds neighbors in O(n^2); intended for up to a few thousand points.
#[wasm_bindgen]
pub fn morans_i(points_xy: &[f64], values: &[f64], distance_threshold: f64) -> Result<String, JsValue> {
    let n = values.len();
    if points_xy.len() != n * 2 {
        return Err(JsValue::from_str("points_xy length must be 2 * values length"));
    }
    if n < 3 {
        return Err(JsValue::from_str("need >= 3 features"));
    }
    if !(distance_threshold > 0.0) {
        return Err(JsValue::from_str("distance_threshold must be > 0"));
    }
    let pts: Vec<(f64, f64)> = points_xy.chunks_exact(2).map(|c| (c[0], c[1])).collect();
    let t2 = distance_threshold * distance_threshold;
    let mut neighbors: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for i in 0..n {
        for j in 0..n {
            if i == j { continue; }
            let dx = pts[i].0 - pts[j].0;
            let dy = pts[i].1 - pts[j].1;
            if dx * dx + dy * dy <= t2 {
                neighbors[i].push((j, 1.0));
            }
        }
    }
    let counts: Vec<usize> = neighbors.iter().map(|v| v.len()).collect();
    let diagnostics = SpatialWeightsDiagnostics {
        n_features: n,
        n_islands: counts.iter().filter(|&&c| c == 0).count(),
        neighbor_count_min: *counts.iter().min().unwrap_or(&0),
        neighbor_count_mean: counts.iter().sum::<usize>() as f64 / n as f64,
        neighbor_count_max: *counts.iter().max().unwrap_or(&0),
        connected_component_count: wbspatialstats::weights::connected_components(&neighbors),
        row_standardized: false,
        dropped_feature_count: 0,
    };
    let weights = SpatialWeightsGraph { neighbors, diagnostics, warnings: Vec::new() };
    let r = wbspatialstats::autocorrelation::morans_i(values, &weights)
        .map_err(|e| JsValue::from_str(&format!("morans_i: {e}")))?;
    Ok(format!(
        "{{\"ok\":true,\"morans_i\":{},\"expected\":{},\"variance\":{},\"z_score\":{},\"p_value\":{},\"n\":{}}}",
        r.statistic, r.expected_value, r.variance, r.z_score, r.p_value, r.n_features
    ))
}
