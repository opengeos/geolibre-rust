//! GeoLibre tool: geographically weighted regression (GWR).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Geographically Weighted Regression*
//! (Spatial Statistics): fit a separate local linear regression at every feature,
//! weighting nearby observations more heavily with a distance-decay kernel, so
//! the relationship between the dependent and explanatory variables is allowed to
//! vary across space. The bundled suite has only global models (OLS, RF, SVM,
//! KNN); GWR is the flagship spatially-varying model.
//!
//! At each feature `i` the local coefficients solve the weighted least squares
//! `β_i = (Xᵀ W_i X)⁻¹ Xᵀ W_i y`, with `W_i` a diagonal of kernel weights
//! `w_ij = k(d_ij / b_i)`:
//!
//! - `gaussian` — `exp(−½ (d/b)²)` (all neighbours contribute),
//! - `bisquare` — `(1 − (d/b)²)²` for `d < b`, else 0 (compact support).
//!
//! The bandwidth `b_i` is either a **fixed** distance or **adaptive** (the
//! distance to the k-th nearest neighbour, so it shrinks where points are dense).
//! When no bandwidth is given it is chosen by minimising the corrected AIC
//! (AICc) over a search of candidate bandwidths.
//!
//! Each output feature keeps its attributes and gains: `predicted`, `residual`,
//! `local_r2`, and a `b_<name>` coefficient per term (intercept + each
//! explanatory field). The tool report carries the global diagnostics —
//! effective number of parameters (trace of the hat matrix), R², adjusted R²,
//! AICc, and residual sum of squares.
//!
//! The fit is O(n²) in the number of features (a dense per-feature neighbour
//! sweep and a small p×p solve), so it suits moderate feature counts.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct GeographicallyWeightedRegressionTool;

impl Tool for GeographicallyWeightedRegressionTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "geographically_weighted_regression",
            display_name: "Geographically Weighted Regression",
            summary: "Local linear regression with distance-decay kernel weights (GWR): per-feature coefficients, local R², residuals and predictions, plus global AICc / R² diagnostics; fixed or adaptive bandwidth, optionally AICc-optimized.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input feature layer (points, or other geometries via their representative point).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "y_field",
                    description: "Dependent (response) numeric field.",
                    required: true,
                },
                ToolParamSpec {
                    name: "x_fields",
                    description: "Comma-separated explanatory numeric field(s).",
                    required: true,
                },
                ToolParamSpec {
                    name: "kernel",
                    description: "Distance-decay kernel: 'gaussian' (default) or 'bisquare'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "bandwidth_type",
                    description: "'adaptive' (default; bandwidth = k nearest neighbours) or 'fixed' (bandwidth = a distance).",
                    required: false,
                },
                ToolParamSpec {
                    name: "bandwidth",
                    description: "Bandwidth: neighbour count for adaptive, or distance for fixed. Omit to select it by minimizing AICc.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["input", "y_field", "x_fields"] {
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
        let y_field = require_str(args, "y_field")?;
        let x_fields: Vec<String> = require_str(args, "x_fields")?
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(str::to_string)
            .collect();
        if x_fields.is_empty() {
            return Err(ToolError::Validation(
                "'x_fields' must list at least one explanatory field".to_string(),
            ));
        }
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let schema = layer.schema.clone();

        // Collect locations and the design matrix (intercept + explanatory).
        let mut locs: Vec<(f64, f64)> = Vec::new();
        let mut xs: Vec<Vec<f64>> = Vec::new();
        let mut ys: Vec<f64> = Vec::new();
        let mut idx_map: Vec<usize> = Vec::new(); // output feature index per row
        for (fi, feature) in layer.features.iter().enumerate() {
            let Some((px, py)) = feature.geometry.as_ref().and_then(rep_point) else {
                continue;
            };
            let y = feature
                .get(&schema, y_field)
                .ok()
                .and_then(FieldValue::as_f64);
            let mut row = Vec::with_capacity(x_fields.len() + 1);
            row.push(1.0); // intercept
            let mut ok = y.is_some();
            for xf in &x_fields {
                match feature.get(&schema, xf).ok().and_then(FieldValue::as_f64) {
                    Some(v) if v.is_finite() => row.push(v),
                    _ => ok = false,
                }
            }
            if let (true, Some(y)) = (ok && y.unwrap().is_finite(), y) {
                locs.push((px, py));
                xs.push(row);
                ys.push(y);
                idx_map.push(fi);
            }
        }
        let n = ys.len();
        let p = x_fields.len() + 1;
        if n <= p + 1 {
            return Err(ToolError::Execution(format!(
                "need more than {} observations with valid y and x values, found {n}",
                p + 1
            )));
        }

        let data = Data { locs, xs, ys, p };

        // Bandwidth: explicit, or AICc-optimized over candidates.
        let bandwidth = match prm.bandwidth {
            Some(b) => b,
            None => {
                ctx.progress.info("selecting bandwidth by AICc");
                optimize_bandwidth(&data, prm.kernel, prm.adaptive)?
            }
        };
        ctx.progress.info(&format!(
            "GWR: {n} obs, {p} term(s), {} bandwidth {bandwidth:.4}, {} kernel",
            if prm.adaptive { "adaptive" } else { "fixed" },
            prm.kernel.as_str()
        ));

        let fit = gwr_fit(&data, bandwidth, prm.kernel, prm.adaptive).ok_or_else(|| {
            ToolError::Execution(
                "GWR fit failed (singular local system; try a larger bandwidth)".to_string(),
            )
        })?;

        // Append output fields (schema then per-feature values, kept aligned).
        let mut field_names: Vec<String> =
            vec!["predicted".into(), "residual".into(), "local_r2".into()];
        field_names.push("b_intercept".into());
        for xf in &x_fields {
            field_names.push(format!("b_{xf}"));
        }
        for name in &field_names {
            layer.add_field(FieldDef::new(name.clone(), FieldType::Float));
        }
        // Rows -> output features; features without a valid observation get nulls.
        let mut row_for_feature: Vec<Option<usize>> = vec![None; layer.features.len()];
        for (row, &fi) in idx_map.iter().enumerate() {
            row_for_feature[fi] = Some(row);
        }
        let extra = field_names.len();
        for (fi, feature) in layer.features.iter_mut().enumerate() {
            match row_for_feature[fi] {
                Some(r) => {
                    feature.attributes.push(FieldValue::Float(fit.predicted[r]));
                    feature.attributes.push(FieldValue::Float(fit.residual[r]));
                    feature.attributes.push(FieldValue::Float(fit.local_r2[r]));
                    for c in 0..p {
                        feature.attributes.push(FieldValue::Float(fit.coefs[r][c]));
                    }
                }
                None => {
                    for _ in 0..extra {
                        feature.attributes.push(FieldValue::Null);
                    }
                }
            }
        }

        ctx.progress.info(&format!(
            "R2 {:.4}, adjusted R2 {:.4}, AICc {:.2}, effective params {:.2}",
            fit.r2, fit.adj_r2, fit.aicc, fit.trace_s
        ));

        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("observations".to_string(), json!(n));
        outputs.insert("terms".to_string(), json!(p));
        outputs.insert("kernel".to_string(), json!(prm.kernel.as_str()));
        outputs.insert(
            "bandwidth_type".to_string(),
            json!(if prm.adaptive { "adaptive" } else { "fixed" }),
        );
        outputs.insert("bandwidth".to_string(), json!(bandwidth));
        outputs.insert("effective_params".to_string(), json!(fit.trace_s));
        outputs.insert("r2".to_string(), json!(fit.r2));
        outputs.insert("adjusted_r2".to_string(), json!(fit.adj_r2));
        outputs.insert("aicc".to_string(), json!(fit.aicc));
        outputs.insert("residual_ss".to_string(), json!(fit.rss));
        Ok(ToolRunResult { outputs })
    }
}

// ── GWR core ──────────────────────────────────────────────────────────────────

struct Data {
    locs: Vec<(f64, f64)>,
    xs: Vec<Vec<f64>>, // n rows, each length p (intercept first)
    ys: Vec<f64>,
    p: usize,
}

struct Fit {
    coefs: Vec<Vec<f64>>,
    predicted: Vec<f64>,
    residual: Vec<f64>,
    local_r2: Vec<f64>,
    trace_s: f64,
    rss: f64,
    r2: f64,
    adj_r2: f64,
    aicc: f64,
}

/// Fits GWR at the given bandwidth. Returns `None` if every local system is
/// singular (no usable fit).
// Range loops index parallel matrices/vectors (design rows, distances); the
// index form is clearer than zipping several slices here.
#[allow(clippy::needless_range_loop)]
fn gwr_fit(data: &Data, bandwidth: f64, kernel: Kernel, adaptive: bool) -> Option<Fit> {
    let n = data.ys.len();
    let p = data.p;
    let mut coefs = vec![vec![f64::NAN; p]; n];
    let mut predicted = vec![f64::NAN; n];
    let mut residual = vec![0.0; n];
    let mut local_r2 = vec![0.0; n];
    let mut trace_s = 0.0;
    let mut solved = 0usize;

    // Precompute distances lazily per point (O(n) each).
    for i in 0..n {
        let (xi, yi) = data.locs[i];
        // Distances to all points.
        let d: Vec<f64> = (0..n)
            .map(|j| {
                let (xj, yj) = data.locs[j];
                ((xi - xj).powi(2) + (yi - yj).powi(2)).sqrt()
            })
            .collect();
        let b = if adaptive {
            // Bandwidth = distance to the k-th nearest neighbour (k = bandwidth).
            let k = (bandwidth as usize).clamp(p + 1, n);
            let mut sorted = d.clone();
            sorted.sort_by(f64::total_cmp);
            sorted[k - 1].max(f64::MIN_POSITIVE)
        } else {
            bandwidth
        };

        // Weighted normal equations A = XᵀWX (p×p), rhs = XᵀWy.
        let mut a = vec![vec![0.0; p]; p];
        let mut rhs = vec![0.0; p];
        let (mut sw, mut swy, mut swy2) = (0.0, 0.0, 0.0);
        for j in 0..n {
            let w = kernel.weight(d[j], b);
            if w <= 0.0 {
                continue;
            }
            let xj = &data.xs[j];
            let yj = data.ys[j];
            for r in 0..p {
                let wr = w * xj[r];
                for c in 0..p {
                    a[r][c] += wr * xj[c];
                }
                rhs[r] += wr * yj;
            }
            sw += w;
            swy += w * yj;
            swy2 += w * yj * yj;
        }

        let Some(beta) = solve(&a, &rhs) else {
            continue;
        };
        // Hat diagonal S_ii = x_iᵀ A⁻¹ x_i (self-weight is k(0)=1).
        let Some(z) = solve(&a, &data.xs[i]) else {
            continue;
        };
        let s_ii: f64 = (0..p).map(|c| data.xs[i][c] * z[c]).sum();

        let yhat_i: f64 = (0..p).map(|c| data.xs[i][c] * beta[c]).sum();
        predicted[i] = yhat_i;
        residual[i] = data.ys[i] - yhat_i;
        coefs[i] = beta.clone();
        trace_s += s_ii;
        solved += 1;

        // Local weighted R².
        if sw > 0.0 {
            let ybar = swy / sw;
            let tss = swy2 - sw * ybar * ybar; // Σ w (y-ȳ)²
            let mut rss_w = 0.0;
            for j in 0..n {
                let w = kernel.weight(d[j], b);
                if w <= 0.0 {
                    continue;
                }
                let yhat: f64 = (0..p).map(|c| data.xs[j][c] * beta[c]).sum();
                rss_w += w * (data.ys[j] - yhat).powi(2);
            }
            local_r2[i] = if tss > 0.0 {
                (1.0 - rss_w / tss).clamp(0.0, 1.0)
            } else {
                0.0
            };
        }
    }

    if solved == 0 {
        return None;
    }

    // Global diagnostics.
    let ybar = data.ys.iter().sum::<f64>() / n as f64;
    let tss: f64 = data.ys.iter().map(|y| (y - ybar).powi(2)).sum();
    let rss: f64 = residual.iter().map(|r| r * r).sum();
    let r2 = if tss > 0.0 { 1.0 - rss / tss } else { 0.0 };
    let nf = n as f64;
    // GWR AICc (Fotheringham et al.): σ̂ = sqrt(RSS/n).
    let sigma = (rss / nf).sqrt().max(f64::MIN_POSITIVE);
    let aicc = 2.0 * nf * sigma.ln()
        + nf * (2.0 * std::f64::consts::PI).ln()
        + nf * (nf + trace_s) / (nf - 2.0 - trace_s);
    let adj_r2 = if nf - trace_s - 1.0 > 0.0 {
        1.0 - (1.0 - r2) * (nf - 1.0) / (nf - trace_s - 1.0)
    } else {
        r2
    };

    Some(Fit {
        coefs,
        predicted,
        residual,
        local_r2,
        trace_s,
        rss,
        r2,
        adj_r2,
        aicc,
    })
}

/// AICc-minimizing bandwidth over a coarse candidate grid (then a ±1 refine for
/// adaptive), deterministic and dependency-free.
fn optimize_bandwidth(data: &Data, kernel: Kernel, adaptive: bool) -> Result<f64, ToolError> {
    let n = data.ys.len();
    let p = data.p;
    let candidates: Vec<f64> = if adaptive {
        let lo = p + 2;
        let hi = n;
        let steps = 18.min(hi - lo).max(1);
        (0..=steps)
            .map(|s| (lo + s * (hi - lo) / steps) as f64)
            .collect()
    } else {
        // Distances between the median nearest-neighbour distance and the extent.
        let (mut minx, mut miny, mut maxx, mut maxy) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for &(x, y) in &data.locs {
            minx = minx.min(x);
            miny = miny.min(y);
            maxx = maxx.max(x);
            maxy = maxy.max(y);
        }
        let diag = ((maxx - minx).powi(2) + (maxy - miny).powi(2)).sqrt();
        let lo = diag / (n as f64).sqrt() * 0.5;
        let hi = diag;
        (0..=18).map(|s| lo + (hi - lo) * s as f64 / 18.0).collect()
    };

    let mut best = (f64::INFINITY, candidates[0]);
    for &b in &candidates {
        if let Some(fit) = gwr_fit(data, b, kernel, adaptive) {
            if fit.aicc.is_finite() && fit.aicc < best.0 {
                best = (fit.aicc, b);
            }
        }
    }
    if !best.0.is_finite() {
        return Err(ToolError::Execution(
            "bandwidth selection failed (no candidate produced a finite AICc)".to_string(),
        ));
    }
    Ok(best.1)
}

/// Solves the small dense system `a x = b` (Gaussian elimination with partial
/// pivoting). Returns `None` if `a` is singular.
#[allow(clippy::needless_range_loop)]
fn solve(a_in: &[Vec<f64>], b_in: &[f64]) -> Option<Vec<f64>> {
    let n = b_in.len();
    let mut a: Vec<Vec<f64>> = a_in.to_vec();
    let mut b = b_in.to_vec();
    for col in 0..n {
        // Partial pivot.
        let mut piv = col;
        let mut best = a[col][col].abs();
        for r in col + 1..n {
            if a[r][col].abs() > best {
                best = a[r][col].abs();
                piv = r;
            }
        }
        if best < 1e-12 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        let d = a[col][col];
        for r in col + 1..n {
            let f = a[r][col] / d;
            if f == 0.0 {
                continue;
            }
            for c in col..n {
                a[r][c] -= f * a[col][c];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for c in i + 1..n {
            s -= a[i][c] * x[c];
        }
        x[i] = s / a[i][i];
    }
    Some(x)
}

// ── Kernel ────────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kernel {
    Gaussian,
    Bisquare,
}

impl Kernel {
    fn as_str(self) -> &'static str {
        match self {
            Self::Gaussian => "gaussian",
            Self::Bisquare => "bisquare",
        }
    }
    fn weight(self, d: f64, b: f64) -> f64 {
        if b <= 0.0 {
            return 0.0;
        }
        let t = d / b;
        match self {
            Self::Gaussian => (-0.5 * t * t).exp(),
            Self::Bisquare => {
                if t < 1.0 {
                    let u = 1.0 - t * t;
                    u * u
                } else {
                    0.0
                }
            }
        }
    }
}

// ── Helpers / parameters ──────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .map(str::trim)
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

fn rep_point(geom: &Geometry) -> Option<(f64, f64)> {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0usize;
    collect(geom, &mut |x, y| {
        sx += x;
        sy += y;
        n += 1;
    });
    (n > 0).then(|| (sx / n as f64, sy / n as f64))
}

fn collect(geom: &Geometry, f: &mut impl FnMut(f64, f64)) {
    match geom {
        Geometry::Point(c) => f(c.x, c.y),
        Geometry::MultiPoint(cs) | Geometry::LineString(cs) => cs.iter().for_each(|c| f(c.x, c.y)),
        Geometry::MultiLineString(ls) => ls.iter().flatten().for_each(|c| f(c.x, c.y)),
        Geometry::Polygon { exterior, .. } => exterior.coords().iter().for_each(|c| f(c.x, c.y)),
        Geometry::MultiPolygon(parts) => parts
            .iter()
            .flat_map(|(e, _)| e.coords())
            .for_each(|c| f(c.x, c.y)),
        Geometry::GeometryCollection(gs) => gs.iter().for_each(|g| collect(g, f)),
    }
}

struct Params {
    kernel: Kernel,
    adaptive: bool,
    bandwidth: Option<f64>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let kernel = match parse_optional_str(args, "kernel")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("gaussian") => Kernel::Gaussian,
        Some("bisquare") => Kernel::Bisquare,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown kernel '{other}' (expected gaussian or bisquare)"
            )))
        }
    };
    let adaptive = match parse_optional_str(args, "bandwidth_type")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("adaptive") => true,
        Some("fixed") => false,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown bandwidth_type '{other}' (expected adaptive or fixed)"
            )))
        }
    };
    let bandwidth = parse_optional_f64(args, "bandwidth")?;
    if let Some(b) = bandwidth {
        if !(b > 0.0 && b.is_finite()) {
            return Err(ToolError::Validation(
                "parameter 'bandwidth' must be a positive number".to_string(),
            ));
        }
        if adaptive && b.fract() != 0.0 {
            return Err(ToolError::Validation(
                "adaptive 'bandwidth' must be an integer neighbour count".to_string(),
            ));
        }
    }
    Ok(Params {
        kernel,
        adaptive,
        bandwidth,
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

    fn layer_with(rows: &[(f64, f64, f64, f64)]) -> String {
        // (x, y, x1, yval)
        let mut layer = Layer::new("pts");
        layer.add_field(FieldDef::new("x1", FieldType::Float));
        layer.add_field(FieldDef::new("yv", FieldType::Float));
        for &(px, py, x1, yv) in rows {
            layer
                .add_feature(
                    Some(Geometry::point(px, py)),
                    &[("x1", FieldValue::Float(x1)), ("yv", FieldValue::Float(yv))],
                )
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = GeographicallyWeightedRegressionTool
            .run(&args, &ctx())
            .unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn field(layer: &Layer, i: usize, name: &str) -> f64 {
        layer.features[i]
            .get(&layer.schema, name)
            .unwrap()
            .as_f64()
            .unwrap()
    }

    #[test]
    fn recovers_a_global_linear_relationship() {
        // y = 3 + 2*x1 everywhere; GWR should recover it with tiny residuals and
        // R² ≈ 1 regardless of local weighting.
        let mut rows = Vec::new();
        for i in 0..40 {
            let px = (i % 8) as f64;
            let py = (i / 8) as f64;
            let x1 = (i as f64 * 0.37).sin() * 5.0 + px;
            rows.push((px, py, x1, 3.0 + 2.0 * x1));
        }
        let (out, layer) = run(json!({
            "input": layer_with(&rows), "y_field": "yv", "x_fields": "x1",
            "bandwidth_type": "adaptive", "bandwidth": 20
        }));
        assert!(
            out.outputs["r2"].as_f64().unwrap() > 0.999,
            "R2 = {}",
            out.outputs["r2"]
        );
        // Local coefficients ≈ (3, 2).
        assert!((field(&layer, 0, "b_intercept") - 3.0).abs() < 1e-3);
        assert!((field(&layer, 0, "b_x1") - 2.0).abs() < 1e-3);
        assert!(field(&layer, 0, "residual").abs() < 1e-3);
    }

    #[test]
    fn captures_a_spatially_varying_slope() {
        // Slope varies with x-location: left half slope ~1, right half slope ~5.
        // A global OLS can't fit both; GWR's local slopes should differ markedly.
        let mut rows = Vec::new();
        for gx in 0..12 {
            for gy in 0..4 {
                let px = gx as f64;
                let py = gy as f64;
                let slope = if gx < 6 { 1.0 } else { 5.0 };
                let x1 = px + py * 0.1;
                rows.push((px, py, x1, 10.0 + slope * x1));
            }
        }
        let (_, layer) = run(json!({
            "input": layer_with(&rows), "y_field": "yv", "x_fields": "x1",
            "bandwidth_type": "adaptive", "bandwidth": 8
        }));
        // Left-most feature slope near 1, right-most near 5.
        let left = field(&layer, 0, "b_x1");
        let right = field(&layer, layer.len() - 1, "b_x1");
        assert!(left < 2.5, "left slope {left} should be ~1");
        assert!(right > 3.5, "right slope {right} should be ~5");
    }

    #[test]
    fn auto_bandwidth_runs_and_reports_diagnostics() {
        let mut rows = Vec::new();
        for i in 0..50 {
            let px = (i % 10) as f64;
            let py = (i / 10) as f64;
            let x1 = px + 0.5 * py;
            rows.push((px, py, x1, 2.0 + 1.5 * x1 + (i as f64 * 0.9).sin()));
        }
        let (out, _) = run(json!({
            "input": layer_with(&rows), "y_field": "yv", "x_fields": "x1"
        }));
        assert!(out.outputs["bandwidth"].as_f64().unwrap() > 0.0);
        assert!(out.outputs["aicc"].as_f64().unwrap().is_finite());
        assert!(out.outputs["effective_params"].as_f64().unwrap() > 1.0);
        assert!(out.outputs["r2"].as_f64().unwrap() > 0.5);
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = GeographicallyWeightedRegressionTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "input": "p.geojson", "y_field": "yv" })).is_err(),
            "no x_fields"
        );
        assert!(bad(
            json!({ "input": "p.geojson", "y_field": "yv", "x_fields": "x1", "kernel": "tri" })
        )
        .is_err());
        assert!(bad(
            json!({ "input": "p.geojson", "y_field": "yv", "x_fields": "x1", "bandwidth": 0 })
        )
        .is_err());
        assert!(bad(json!({ "input": "p.geojson", "y_field": "yv", "x_fields": "x1", "bandwidth_type": "adaptive", "bandwidth": 5.5 })).is_err());
        assert!(bad(json!({ "input": "p.geojson", "y_field": "yv", "x_fields": "x1" })).is_ok());
    }
}
