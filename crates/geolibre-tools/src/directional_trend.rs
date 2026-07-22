//! GeoLibre tool: 1-D directional trend diagnostic.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Directional Trend* (Geostatistical
//! Analyst) tool. It projects sample points onto a chosen azimuth and fits a
//! 1st-3rd-order polynomial of attribute value vs. distance along that bearing,
//! reporting the fit coefficients and R² (the fraction of variance the trend
//! explains). This is the 1-D anisotropy check a geostatistician runs *before*
//! variography / kriging: a strong directional trend biases the semivariogram.
//!
//! The existing `trend_surface` / `generate_trend_raster` tools fit a 2-D
//! polynomial surface, and `exploratory_regression` is attribute-space
//! regression; none provide this 1-D directional pre-kriging diagnostic.
//!
//! `azimuth` is a compass bearing in degrees (0 = north, clockwise). Passing
//! `determine` (the default) sweeps azimuths 0-179° and returns the bearing of
//! maximum explained variance — the direction along which the attribute trends
//! most strongly. `order` (1-3) is the polynomial degree.
//!
//! The polynomial is fit by least squares (normal equations, Gaussian
//! elimination) in a normalized distance coordinate `t ∈ [-1, 1]` for
//! conditioning; R² is scale-invariant so it is unaffected. The output layer
//! carries each point augmented with its signed projected distance, the fitted
//! value, and the residual — the value-vs-distance table you would scatter-plot.
//!
//! Point and multipoint inputs use their coordinates directly; other geometries
//! use a representative point (the mean of their vertices).

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Azimuth sweep step (degrees) for the `determine` mode.
const SWEEP_STEP_DEG: usize = 1;

pub struct DirectionalTrendTool;

impl Tool for DirectionalTrendTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "directional_trend",
            display_name: "Directional Trend",
            summary: "Project sample points onto a chosen azimuth and fit a 1st-3rd-order polynomial of value-vs-distance to detect a global/anisotropic trend before variography (like ArcGIS's Directional Trend). 'determine' sweeps azimuths for the direction of maximum explained variance.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (points preferred; other geometries use their vertex-mean representative point).",
                    required: true,
                },
                ToolParamSpec {
                    name: "field",
                    description: "Numeric attribute field to analyze the trend of.",
                    required: true,
                },
                ToolParamSpec {
                    name: "azimuth",
                    description: "Compass bearing in degrees (0 = north, clockwise) to project onto, or 'determine' (default) to sweep 0-179° for the bearing of maximum explained variance.",
                    required: false,
                },
                ToolParamSpec {
                    name: "order",
                    description: "Polynomial order 1, 2, or 3 (default 2).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory. Points carry proj_dist, value, fitted, and residual.",
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
        if parse_optional_str(args, "field")?
            .map(str::trim)
            .unwrap_or("")
            .is_empty()
        {
            return Err(ToolError::Validation(
                "missing required string parameter 'field'".to_string(),
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
        let layer_crs = layer.crs.clone();
        let schema = &layer.schema;

        // Collect (x, y, value) observations, skipping missing/non-finite ones.
        let mut obs: Vec<Obs> = Vec::new();
        for feature in &layer.features {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let Some((x, y)) = rep_point(geom) else {
                continue;
            };
            let v = match feature
                .get(schema, &prm.field)
                .ok()
                .and_then(FieldValue::as_f64)
            {
                Some(v) if v.is_finite() => v,
                _ => continue,
            };
            obs.push(Obs { x, y, v });
        }

        let n = obs.len();
        let min_pts = prm.order + 1;
        if n < min_pts {
            return Err(ToolError::Execution(format!(
                "need at least {min_pts} point(s) with a finite '{}' value for an order-{} fit, found {n}",
                prm.field, prm.order
            )));
        }

        // Center coordinates for numerical stability (does not change bearings).
        let cx = obs.iter().map(|o| o.x).sum::<f64>() / n as f64;
        let cy = obs.iter().map(|o| o.y).sum::<f64>() / n as f64;

        // Total variance of the attribute (denominator of R²).
        let vbar = obs.iter().map(|o| o.v).sum::<f64>() / n as f64;
        let ss_tot: f64 = obs.iter().map(|o| (o.v - vbar).powi(2)).sum();
        if ss_tot <= 0.0 {
            return Err(ToolError::Execution(format!(
                "field '{}' is constant across the {n} valid point(s); no trend to fit",
                prm.field
            )));
        }

        // Determine the azimuth: fixed value, or sweep for max explained variance.
        let (azimuth, fit) = match prm.azimuth {
            AzimuthMode::Fixed(az) => {
                let fit = fit_along(&obs, cx, cy, az, prm.order, ss_tot).ok_or_else(|| {
                    ToolError::Execution(
                        "points are coincident along the chosen azimuth; cannot fit".to_string(),
                    )
                })?;
                (az, fit)
            }
            AzimuthMode::Determine => {
                ctx.progress
                    .info("sweeping azimuths for maximum explained variance");
                let mut best: Option<(f64, Fit)> = None;
                let mut az = 0usize;
                while az < 180 {
                    if let Some(fit) = fit_along(&obs, cx, cy, az as f64, prm.order, ss_tot) {
                        if best.as_ref().map(|(_, b)| fit.r2 > b.r2).unwrap_or(true) {
                            best = Some((az as f64, fit));
                        }
                    }
                    az += SWEEP_STEP_DEG;
                }
                best.ok_or_else(|| {
                    ToolError::Execution(
                        "could not fit a trend along any azimuth (points coincident?)".to_string(),
                    )
                })?
            }
        };

        ctx.progress.info(&format!(
            "azimuth {:.0}°, order {}, R² = {:.4}",
            azimuth, prm.order, fit.r2
        ));

        // Build the output point layer: value-vs-distance table.
        let mut out_layer = Layer::new(layer.name.clone());
        out_layer.crs = layer_crs;
        out_layer.geom_type = Some(GeometryType::Point);
        out_layer.add_field(FieldDef::new("proj_dist", FieldType::Float));
        out_layer.add_field(FieldDef::new("value", FieldType::Float));
        out_layer.add_field(FieldDef::new("fitted", FieldType::Float));
        out_layer.add_field(FieldDef::new("residual", FieldType::Float));

        let (rad, dmin, dmax) = azimuth_projection(&obs, cx, cy, azimuth);
        for o in &obs {
            let d = (o.x - cx) * rad.0 + (o.y - cy) * rad.1;
            let fitted = fit.eval(d, dmin, dmax);
            out_layer
                .add_feature(
                    Some(Geometry::point(o.x, o.y)),
                    &[
                        ("proj_dist", FieldValue::Float(d)),
                        ("value", FieldValue::Float(o.v)),
                        ("fitted", FieldValue::Float(fitted)),
                        ("residual", FieldValue::Float(o.v - fitted)),
                    ],
                )
                .map_err(|e| ToolError::Execution(format!("failed writing output feature: {e}")))?;
        }

        let feature_count = out_layer.len();
        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("azimuth".to_string(), json!(azimuth));
        outputs.insert("order".to_string(), json!(prm.order));
        outputs.insert("r_squared".to_string(), json!(fit.r2));
        outputs.insert("rmse".to_string(), json!(fit.rmse));
        outputs.insert(
            "coefficients".to_string(),
            json!(fit.coeffs), // in normalized distance t ∈ [-1, 1]
        );
        outputs.insert("dist_min".to_string(), json!(dmin));
        outputs.insert("dist_max".to_string(), json!(dmax));
        outputs.insert("n_points".to_string(), json!(n));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert(
            "swept".to_string(),
            json!(matches!(prm.azimuth, AzimuthMode::Determine)),
        );
        Ok(ToolRunResult { outputs })
    }
}

// ── Observations & fitting ────────────────────────────────────────────────────

struct Obs {
    x: f64,
    y: f64,
    v: f64,
}

/// A polynomial fit in normalized distance `t = 2*(d - dmin)/(dmax - dmin) - 1`.
struct Fit {
    /// coeffs[k] multiplies t^k; length = order + 1.
    coeffs: Vec<f64>,
    r2: f64,
    rmse: f64,
}

impl Fit {
    /// Evaluates the fitted value at the original signed distance `d`, given the
    /// distance range used to normalize.
    fn eval(&self, d: f64, dmin: f64, dmax: f64) -> f64 {
        let t = normalize(d, dmin, dmax);
        let mut y = 0.0;
        let mut tp = 1.0;
        for &c in &self.coeffs {
            y += c * tp;
            tp *= t;
        }
        y
    }
}

/// Unit vector `(east, north)` of a compass bearing (degrees, clockwise from N).
fn bearing_unit(azimuth_deg: f64) -> (f64, f64) {
    let r = azimuth_deg.to_radians();
    (r.sin(), r.cos())
}

fn normalize(d: f64, dmin: f64, dmax: f64) -> f64 {
    let span = dmax - dmin;
    if span <= 0.0 {
        0.0
    } else {
        2.0 * (d - dmin) / span - 1.0
    }
}

/// Projects observations onto `azimuth` and returns the unit vector plus the
/// min/max signed distance (used to normalize into t-space).
fn azimuth_projection(obs: &[Obs], cx: f64, cy: f64, azimuth: f64) -> ((f64, f64), f64, f64) {
    let unit = bearing_unit(azimuth);
    let mut dmin = f64::INFINITY;
    let mut dmax = f64::NEG_INFINITY;
    for o in obs {
        let d = (o.x - cx) * unit.0 + (o.y - cy) * unit.1;
        dmin = dmin.min(d);
        dmax = dmax.max(d);
    }
    (unit, dmin, dmax)
}

/// Fits an order-`order` polynomial of value vs. distance along `azimuth`.
/// Returns `None` if all points project to the same distance (degenerate).
fn fit_along(
    obs: &[Obs],
    cx: f64,
    cy: f64,
    azimuth: f64,
    order: usize,
    ss_tot: f64,
) -> Option<Fit> {
    let (unit, dmin, dmax) = azimuth_projection(obs, cx, cy, azimuth);
    if dmax - dmin <= 0.0 {
        return None;
    }
    let k = order + 1;
    // Normal equations: (X^T X) c = X^T y, with X_ij = t_i^j.
    let mut ata = vec![vec![0.0f64; k]; k];
    let mut aty = vec![0.0f64; k];
    for o in obs {
        let d = (o.x - cx) * unit.0 + (o.y - cy) * unit.1;
        let t = normalize(d, dmin, dmax);
        // Powers t^0..t^(order).
        let mut powers = vec![1.0; k];
        for p in 1..k {
            powers[p] = powers[p - 1] * t;
        }
        for i in 0..k {
            aty[i] += powers[i] * o.v;
            for j in 0..k {
                ata[i][j] += powers[i] * powers[j];
            }
        }
    }
    let coeffs = solve_linear(ata, aty)?;

    // Residual sum of squares.
    let mut ss_res = 0.0;
    for o in obs {
        let d = (o.x - cx) * unit.0 + (o.y - cy) * unit.1;
        let t = normalize(d, dmin, dmax);
        let mut yhat = 0.0;
        let mut tp = 1.0;
        for &c in &coeffs {
            yhat += c * tp;
            tp *= t;
        }
        ss_res += (o.v - yhat).powi(2);
    }
    let r2 = (1.0 - ss_res / ss_tot).clamp(0.0, 1.0);
    let rmse = (ss_res / obs.len() as f64).sqrt();
    Some(Fit { coeffs, r2, rmse })
}

/// Solves the linear system `a x = b` by Gaussian elimination with partial
/// pivoting. Returns `None` if the matrix is (near-)singular.
#[allow(clippy::needless_range_loop)]
fn solve_linear(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        // Partial pivot.
        let mut pivot = col;
        let mut best = a[col][col].abs();
        for r in (col + 1)..n {
            let v = a[r][col].abs();
            if v > best {
                best = v;
                pivot = r;
            }
        }
        if best < 1e-12 {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        // Eliminate.
        for r in (col + 1)..n {
            let factor = a[r][col] / a[col][col];
            for c in col..n {
                a[r][c] -= factor * a[col][c];
            }
            b[r] -= factor * b[col];
        }
    }
    // Back-substitute.
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in (i + 1)..n {
            s -= a[i][j] * x[j];
        }
        x[i] = s / a[i][i];
    }
    if x.iter().all(|v| v.is_finite()) {
        Some(x)
    } else {
        None
    }
}

/// Representative point of a geometry: the coordinate for a point, otherwise the
/// mean of the geometry's vertices.
fn rep_point(geom: &Geometry) -> Option<(f64, f64)> {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0usize;
    collect_coords(geom, &mut |x, y| {
        sx += x;
        sy += y;
        n += 1;
    });
    (n > 0).then(|| (sx / n as f64, sy / n as f64))
}

fn collect_coords(geom: &Geometry, f: &mut impl FnMut(f64, f64)) {
    match geom {
        Geometry::Point(c) => f(c.x, c.y),
        Geometry::MultiPoint(cs) | Geometry::LineString(cs) => cs.iter().for_each(|c| f(c.x, c.y)),
        Geometry::MultiLineString(ls) => ls.iter().flatten().for_each(|c| f(c.x, c.y)),
        Geometry::Polygon { exterior, .. } => exterior.coords().iter().for_each(|c| f(c.x, c.y)),
        Geometry::MultiPolygon(parts) => parts
            .iter()
            .flat_map(|(e, _)| e.coords())
            .for_each(|c| f(c.x, c.y)),
        Geometry::GeometryCollection(gs) => gs.iter().for_each(|g| collect_coords(g, f)),
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

enum AzimuthMode {
    Fixed(f64),
    Determine,
}

struct Params {
    field: String,
    azimuth: AzimuthMode,
    order: usize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let field = parse_optional_str(args, "field")?
        .map(str::to_string)
        .ok_or_else(|| ToolError::Validation("missing required parameter 'field'".to_string()))?;

    let azimuth = match args.get("azimuth") {
        None | Some(Value::Null) => AzimuthMode::Determine,
        Some(Value::String(s)) if s.trim().is_empty() => AzimuthMode::Determine,
        Some(Value::String(s)) if s.trim().eq_ignore_ascii_case("determine") => {
            AzimuthMode::Determine
        }
        Some(Value::String(s)) => {
            let v = s.trim().parse::<f64>().map_err(|_| {
                ToolError::Validation(
                    "parameter 'azimuth' must be a number in [0, 360) or 'determine'".to_string(),
                )
            })?;
            AzimuthMode::Fixed(normalize_azimuth(v)?)
        }
        Some(Value::Number(nnum)) => {
            let v = nnum.as_f64().ok_or_else(|| {
                ToolError::Validation("parameter 'azimuth' must be a number".to_string())
            })?;
            AzimuthMode::Fixed(normalize_azimuth(v)?)
        }
        Some(_) => {
            return Err(ToolError::Validation(
                "parameter 'azimuth' must be a number or 'determine'".to_string(),
            ))
        }
    };

    let order = match parse_optional_f64(args, "order")? {
        None => 2,
        Some(v) if v.fract() == 0.0 && (1.0..=3.0).contains(&v) => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "parameter 'order' must be 1, 2, or 3".to_string(),
            ))
        }
    };

    Ok(Params {
        field,
        azimuth,
        order,
    })
}

/// Validates an azimuth is finite and wraps it into [0, 360).
fn normalize_azimuth(v: f64) -> Result<f64, ToolError> {
    if !v.is_finite() {
        return Err(ToolError::Validation(
            "parameter 'azimuth' must be finite".to_string(),
        ));
    }
    Ok(v.rem_euclid(360.0))
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
    use wbvector::{memory_store, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a point layer with a single float field `val`.
    fn points_layer(pts: &[(f64, f64, f64)]) -> String {
        let mut layer = Layer::new("pts");
        layer.add_field(FieldDef::new("val", FieldType::Float));
        for &(x, y, v) in pts {
            layer
                .add_feature(
                    Some(Geometry::point(x, y)),
                    &[("val", FieldValue::Float(v))],
                )
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = DirectionalTrendTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn fv(layer: &Layer, idx: usize, name: &str) -> f64 {
        layer.features[idx]
            .get(&layer.schema, name)
            .ok()
            .and_then(FieldValue::as_f64)
            .unwrap()
    }

    #[test]
    fn perfect_linear_trend_along_east_gives_r2_one() {
        // Value = x exactly: a perfect order-1 trend along the east/west line.
        let pts: Vec<(f64, f64, f64)> = (0..20).map(|i| (i as f64, 0.0, i as f64)).collect();
        let (out, layer) = run_tool(json!({
            "input": points_layer(&pts), "field": "val", "azimuth": 90, "order": 1
        }));
        assert!(
            out.outputs["r_squared"].as_f64().unwrap() > 0.9999,
            "R² should be ~1, got {}",
            out.outputs["r_squared"]
        );
        // Residuals are ~0 everywhere.
        for i in 0..layer.len() {
            assert!(fv(&layer, i, "residual").abs() < 1e-6);
        }
        assert_eq!(out.outputs["n_points"], json!(20));
    }

    #[test]
    fn determine_finds_the_trend_bearing() {
        // Value increases toward the north (azimuth 0). The sweep should land
        // near 0° (== 180° as a line) with high R².
        let pts: Vec<(f64, f64, f64)> = (0..20).map(|i| (0.0, i as f64, i as f64)).collect();
        let (out, _) = run_tool(json!({
            "input": points_layer(&pts), "field": "val", "order": 1
        }));
        let az = out.outputs["azimuth"].as_f64().unwrap();
        // North-south line: bearing 0 (swept 0..179 gives 0).
        assert!(
            az < 3.0 || az > 177.0,
            "expected a ~north-south bearing, got {az}"
        );
        assert!(out.outputs["r_squared"].as_f64().unwrap() > 0.999);
        assert_eq!(out.outputs["swept"], json!(true));
    }

    #[test]
    fn no_trend_perpendicular_to_gradient_gives_low_r2() {
        // Value = x, but we project onto the north/south line (azimuth 0),
        // orthogonal to the gradient: R² should be ~0.
        let pts: Vec<(f64, f64, f64)> = (0..10)
            .flat_map(|i| (0..10).map(move |j| (i as f64, j as f64, i as f64)))
            .collect();
        let (out, _) = run_tool(json!({
            "input": points_layer(&pts), "field": "val", "azimuth": 0, "order": 1
        }));
        assert!(
            out.outputs["r_squared"].as_f64().unwrap() < 0.01,
            "R² should be ~0 orthogonal to the gradient, got {}",
            out.outputs["r_squared"]
        );
    }

    #[test]
    fn quadratic_trend_needs_order_two() {
        // Value = x^2 along east: order-1 fit is poor, order-2 is near perfect.
        let pts: Vec<(f64, f64, f64)> = (-10..=10)
            .map(|i| (i as f64, 0.0, (i * i) as f64))
            .collect();
        let (lin, _) = run_tool(json!({
            "input": points_layer(&pts), "field": "val", "azimuth": 90, "order": 1
        }));
        let (quad, _) = run_tool(json!({
            "input": points_layer(&pts), "field": "val", "azimuth": 90, "order": 2
        }));
        let r2_lin = lin.outputs["r_squared"].as_f64().unwrap();
        let r2_quad = quad.outputs["r_squared"].as_f64().unwrap();
        assert!(r2_quad > 0.9999, "order-2 should fit x^2, got {r2_quad}");
        assert!(r2_quad > r2_lin + 0.1, "order-2 must beat order-1");
    }

    #[test]
    fn output_carries_one_point_per_valid_input() {
        let pts: Vec<(f64, f64, f64)> = (0..15).map(|i| (i as f64, 0.0, i as f64)).collect();
        let (out, layer) = run_tool(json!({
            "input": points_layer(&pts), "field": "val", "azimuth": 90, "order": 1
        }));
        assert_eq!(layer.len(), 15);
        assert_eq!(out.outputs["feature_count"], json!(15));
        // Projected distance must be monotonic with x here.
        assert!(fv(&layer, 14, "proj_dist") > fv(&layer, 0, "proj_dist"));
    }

    #[test]
    fn skips_features_with_missing_values() {
        let mut layer = Layer::new("pts");
        layer.add_field(FieldDef::new("val", FieldType::Float));
        for i in 0..10 {
            let fields: Vec<(&str, FieldValue)> = if i == 3 {
                vec![("val", FieldValue::Null)]
            } else {
                vec![("val", FieldValue::Float(i as f64))]
            };
            layer
                .add_feature(Some(Geometry::point(i as f64, 0.0)), &fields)
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);
        let (out, layer) =
            run_tool(json!({ "input": input, "field": "val", "azimuth": 90, "order": 1 }));
        // One point dropped (the null), leaving 9.
        assert_eq!(out.outputs["n_points"], json!(9));
        assert_eq!(layer.len(), 9);
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = DirectionalTrendTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input");
        assert!(
            bad(json!({ "input": "x.geojson" })).is_err(),
            "missing field"
        );
        assert!(
            bad(json!({ "input": "x.geojson", "field": "val", "order": 4 })).is_err(),
            "order out of range"
        );
        assert!(
            bad(json!({ "input": "x.geojson", "field": "val", "order": 1.5 })).is_err(),
            "non-integer order"
        );
        assert!(
            bad(json!({ "input": "x.geojson", "field": "val", "azimuth": "bogus" })).is_err(),
            "non-numeric azimuth"
        );
        assert!(
            bad(json!({ "input": "x.geojson", "field": "val", "azimuth": "determine" })).is_ok()
        );
        assert!(
            bad(json!({ "input": "x.geojson", "field": "val", "azimuth": 45, "order": 3 })).is_ok()
        );
    }

    #[test]
    fn errors_when_too_few_points() {
        // Order 2 needs >= 3 points; give 2.
        let input = points_layer(&[(0.0, 0.0, 1.0), (1.0, 0.0, 2.0)]);
        let args: ToolArgs = serde_json::from_value(
            json!({ "input": input, "field": "val", "azimuth": 90, "order": 2 }),
        )
        .unwrap();
        assert!(DirectionalTrendTool.run(&args, &ctx()).is_err());
    }
}
