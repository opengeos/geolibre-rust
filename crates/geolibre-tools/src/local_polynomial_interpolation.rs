//! GeoLibre tool: local polynomial interpolation (LPI).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Local Polynomial Interpolation*
//! (Geostatistical Analyst). The bundled interpolation suite is broad —
//! `idw_interpolation`, the kriging/variogram family,
//! `natural_neighbour_interpolation`, `radial_basis_function_interpolation`,
//! `thin_plate_spline`, `tin_interpolation` — plus GeoLibre's
//! `interpolate_with_barriers`, `diffusion_interpolation_with_barriers`,
//! `optimal_interpolation` and `empirical_bayesian_kriging`. But the only
//! polynomial fit anywhere in that list is the bundled `trend_surface`, which
//! is a **global** polynomial: one equation for the entire extent.
//!
//! LPI is the local counterpart — an order-1/2/3 polynomial refitted inside a
//! kernel-weighted moving neighbourhood at every cell. It suits surfaces with
//! smooth regional structure plus local variation (contaminant plumes,
//! temperature fields, soil chemistry), and it is the only interpolator in the
//! catalog that also yields a **prediction standard error** and a **condition
//! number** surface, which is what makes a result defensible rather than
//! merely plausible.
//!
//! The local weighted least-squares solve is the same machinery
//! `geographically_weighted_regression` and `mgwr` already use; only the design
//! matrix changes, from covariate columns to polynomial terms in (dx, dy).
//! Fitting in *local* coordinates centred on the target cell both conditions
//! the system far better and makes the prediction simply the intercept.

use std::collections::BTreeMap;

use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};

use crate::common::{parse_optional_output, raster_like_with_data};
use crate::vector_common::{load_input_layer, parse_optional_str};

const NODATA: f64 = -9999.0;

pub struct LocalPolynomialInterpolationTool;

impl Tool for LocalPolynomialInterpolationTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "local_polynomial_interpolation",
            display_name: "Local Polynomial Interpolation",
            summary: "Fit an order-1/2/3 polynomial within a kernel-weighted moving neighbourhood at every cell, producing a prediction surface, a prediction standard error surface, or a condition-number surface. Unlike the bundled global trend_surface, the polynomial is refitted locally. Like ArcGIS Local Polynomial Interpolation.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point features carrying the value to interpolate.",
                    required: true,
                },
                ToolParamSpec {
                    name: "z_field",
                    description: "Numeric field holding the value to interpolate.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output raster path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Output cell size in map units (default: input extent / 200).",
                    required: false,
                },
                ToolParamSpec {
                    name: "order",
                    description: "Polynomial order: 1 (default), 2 or 3.",
                    required: false,
                },
                ToolParamSpec {
                    name: "kernel",
                    description: "'exponential' (default), 'gaussian', 'quartic', 'epanechnikov', 'fifth_order' or 'constant'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "bandwidth",
                    description: "Fixed kernel bandwidth in map units. If omitted, an adaptive bandwidth is used: the distance to the 'neighbors'-th nearest sample.",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighbors",
                    description: "Samples per local fit for the adaptive bandwidth (default: 3x the number of polynomial terms).",
                    required: false,
                },
                ToolParamSpec {
                    name: "weight_field",
                    description: "Optional per-observation measurement weight field.",
                    required: false,
                },
                ToolParamSpec {
                    name: "condition_number",
                    description: "Refuse to predict where the local design's condition number exceeds this threshold (default: no limit).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_type",
                    description: "'prediction' (default), 'standard_error' or 'condition_number'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "epsg",
                    description: "EPSG code stamped on the output raster (default: the input layer's CRS).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "z_field")?;
        let order = parse_order(args)?;
        parse_kernel(args)?;
        parse_output_type(args)?;
        if let Some(b) = parse_optional_f64(args, "bandwidth")? {
            if b <= 0.0 {
                return Err(ToolError::Validation(
                    "'bandwidth' must be positive".to_string(),
                ));
            }
        }
        if let Some(n) = parse_optional_f64(args, "neighbors")? {
            let terms = n_terms(order);
            if (n as usize) < terms {
                return Err(ToolError::Validation(format!(
                    "'neighbors' must be at least {terms} for an order-{order} polynomial"
                )));
            }
        }
        if let Some(c) = parse_optional_f64(args, "cell_size")? {
            if c <= 0.0 {
                return Err(ToolError::Validation(
                    "'cell_size' must be positive".to_string(),
                ));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let z_field = require_str(args, "z_field")?;
        let output = parse_optional_output(args, "output")?;
        let order = parse_order(args)?;
        let kernel = parse_kernel(args)?;
        let out_type = parse_output_type(args)?;
        let fixed_bw = parse_optional_f64(args, "bandwidth")?;
        let terms = n_terms(order);
        let neighbors = parse_optional_f64(args, "neighbors")?
            .map(|v| v as usize)
            .unwrap_or(terms * 3);
        let cond_limit = parse_optional_f64(args, "condition_number")?;
        let weight_field = parse_optional_str(args, "weight_field")?;

        let layer = load_input_layer(input)?;
        let z_idx = layer
            .schema
            .field_index(z_field)
            .ok_or_else(|| ToolError::Validation(format!("z_field '{z_field}' not found")))?;
        let w_idx = match weight_field {
            Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                ToolError::Validation(format!("weight_field '{f}' not found"))
            })?),
            None => None,
        };

        // Collect samples.
        let (mut xs, mut ys, mut zs, mut ws) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for f in layer.iter() {
            let Some(g) = &f.geometry else { continue };
            let coords = g.all_coords();
            let Some(c) = coords.first() else { continue };
            let Some(z) = f.attributes.get(z_idx).and_then(|v| v.as_f64()) else {
                continue;
            };
            if !z.is_finite() {
                continue;
            }
            let w = w_idx
                .and_then(|i| f.attributes.get(i).and_then(|v| v.as_f64()))
                .unwrap_or(1.0)
                .max(0.0);
            xs.push(c.x);
            ys.push(c.y);
            zs.push(z);
            ws.push(w);
        }
        let n = xs.len();
        if n < terms {
            return Err(ToolError::Execution(format!(
                "an order-{order} fit needs at least {terms} samples, found {n}"
            )));
        }

        // Output grid from the sample extent.
        let (min_x, max_x) = bounds(&xs);
        let (min_y, max_y) = bounds(&ys);
        let cell = match parse_optional_f64(args, "cell_size")? {
            Some(c) => c,
            None => ((max_x - min_x).max(max_y - min_y) / 200.0).max(f64::MIN_POSITIVE),
        };
        let cols = (((max_x - min_x) / cell).ceil() as usize).max(1);
        let rows = (((max_y - min_y) / cell).ceil() as usize).max(1);
        if cols.saturating_mul(rows) > 50_000_000 {
            return Err(ToolError::Execution(format!(
                "cell_size {cell} would produce a {rows}x{cols} raster; use a larger cell_size"
            )));
        }

        let mut tree: KdTree<f64, usize, [f64; 2]> = KdTree::new(2);
        for i in 0..n {
            tree.add([xs[i], ys[i]], i).ok();
        }
        ctx.progress.info(&format!(
            "fitting order-{order} local polynomials over {rows}x{cols} cells"
        ));

        let k = neighbors.min(n);
        let mut data = vec![NODATA; rows * cols];
        let mut fitted = 0usize;
        let mut rejected = 0usize;

        for row in 0..rows {
            let y0 = max_y - (row as f64 + 0.5) * cell;
            for col in 0..cols {
                let x0 = min_x + (col as f64 + 0.5) * cell;
                let Ok(found) = tree.nearest(&[x0, y0], k, &squared_euclidean) else {
                    continue;
                };
                if found.len() < terms {
                    continue;
                }
                // Adaptive bandwidth = distance to the farthest of the k
                // neighbours, so the kernel always spans the sample set.
                let bw = fixed_bw.unwrap_or_else(|| {
                    found
                        .iter()
                        .map(|(d2, _)| d2.sqrt())
                        .fold(0.0_f64, f64::max)
                        .max(f64::MIN_POSITIVE)
                });

                let mut design: Vec<Vec<f64>> = Vec::with_capacity(found.len());
                let mut obs: Vec<f64> = Vec::with_capacity(found.len());
                let mut wts: Vec<f64> = Vec::with_capacity(found.len());
                for (d2, idx) in found.iter() {
                    let i = **idx;
                    let w = kernel.weight(d2.sqrt(), bw) * ws[i];
                    if w <= 0.0 {
                        continue;
                    }
                    // Local coordinates: the fit is centred on the target cell,
                    // so the intercept IS the prediction and the system stays
                    // well scaled.
                    design.push(poly_terms(xs[i] - x0, ys[i] - y0, order));
                    obs.push(zs[i]);
                    wts.push(w);
                }
                if design.len() < terms {
                    continue;
                }

                let Some(fit) = weighted_fit(&design, &obs, &wts, terms) else {
                    continue;
                };
                if let Some(limit) = cond_limit {
                    if fit.condition > limit {
                        rejected += 1;
                        continue;
                    }
                }
                let value = match out_type {
                    OutputType::Prediction => fit.beta0,
                    OutputType::StandardError => fit.std_error,
                    OutputType::ConditionNumber => fit.condition,
                };
                if value.is_finite() {
                    data[row * cols + col] = value;
                    fitted += 1;
                }
            }
            ctx.progress.progress((row as f64 + 1.0) / rows as f64);
        }

        let epsg = parse_optional_f64(args, "epsg")?
            .map(|v| v as u32)
            .or_else(|| layer.crs_epsg());
        let crs = epsg.map(CrsInfo::from_epsg).unwrap_or_default();
        let template = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: min_x,
            y_min: max_y - rows as f64 * cell,
            cell_size: cell,
            cell_size_y: Some(cell),
            nodata: NODATA,
            data_type: DataType::F32,
            crs,
            metadata: Default::default(),
        });
        let raster = raster_like_with_data(&template, data, NODATA, DataType::F32)?;
        let out_path = crate::common::write_or_store_output(raster, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("sample_count".to_string(), json!(n));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        outputs.insert("cell_size".to_string(), json!(cell));
        outputs.insert("fitted_cells".to_string(), json!(fitted));
        outputs.insert("rejected_cells".to_string(), json!(rejected));
        outputs.insert("order".to_string(), json!(order));
        outputs.insert("terms".to_string(), json!(terms));
        outputs.insert("output_type".to_string(), json!(out_type.name()));
        Ok(ToolRunResult { outputs })
    }
}

// ── Local weighted least squares ────────────────────────────────────────────

struct Fit {
    /// Intercept = the prediction, because the design is centred on the cell.
    beta0: f64,
    std_error: f64,
    condition: f64,
}

/// Polynomial terms in local coordinates, ordered so index 0 is the intercept.
fn poly_terms(dx: f64, dy: f64, order: usize) -> Vec<f64> {
    let mut t = vec![1.0, dx, dy];
    if order >= 2 {
        t.extend([dx * dx, dx * dy, dy * dy]);
    }
    if order >= 3 {
        t.extend([dx * dx * dx, dx * dx * dy, dx * dy * dy, dy * dy * dy]);
    }
    t
}

fn n_terms(order: usize) -> usize {
    match order {
        1 => 3,
        2 => 6,
        _ => 10,
    }
}

fn weighted_fit(design: &[Vec<f64>], obs: &[f64], wts: &[f64], p: usize) -> Option<Fit> {
    let n = design.len();
    // Normal equations: (X'WX) beta = X'Wy
    let mut ata = vec![0.0; p * p];
    let mut atb = vec![0.0; p];
    for i in 0..n {
        let w = wts[i];
        for a in 0..p {
            for b in 0..p {
                ata[a * p + b] += w * design[i][a] * design[i][b];
            }
            atb[a] += w * design[i][a] * obs[i];
        }
    }
    let condition = condition_number(&ata, p);
    let beta = solve(&ata.clone(), &atb, p)?;

    // Weighted residual sum of squares -> sigma^2 -> Var(beta_0).
    let mut wrss = 0.0;
    for i in 0..n {
        let pred: f64 = (0..p).map(|a| beta[a] * design[i][a]).sum();
        let r = obs[i] - pred;
        wrss += wts[i] * r * r;
    }
    let dof = n as f64 - p as f64;
    let std_error = if dof > 0.0 {
        let sigma2 = wrss / dof;
        // [(X'WX)^-1]_00 via solving (X'WX) v = e0.
        let mut e0 = vec![0.0; p];
        e0[0] = 1.0;
        match solve(&ata, &e0, p) {
            Some(v) if v[0] >= 0.0 => (sigma2 * v[0]).sqrt(),
            _ => f64::NAN,
        }
    } else {
        f64::NAN
    };

    Some(Fit {
        beta0: beta[0],
        std_error,
        condition,
    })
}

/// Gaussian elimination with partial pivoting; `None` when singular.
fn solve(a_in: &[f64], b_in: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut a = a_in.to_vec();
    let mut b = b_in.to_vec();
    for col in 0..n {
        let mut piv = col;
        for r in (col + 1)..n {
            if a[r * n + col].abs() > a[piv * n + col].abs() {
                piv = r;
            }
        }
        if a[piv * n + col].abs() < 1e-14 {
            return None;
        }
        if piv != col {
            for c in 0..n {
                a.swap(col * n + c, piv * n + c);
            }
            b.swap(col, piv);
        }
        let d = a[col * n + col];
        for r in (col + 1)..n {
            let f = a[r * n + col] / d;
            if f == 0.0 {
                continue;
            }
            for c in col..n {
                a[r * n + c] -= f * a[col * n + c];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let mut s = b[row];
        for c in (row + 1)..n {
            s -= a[row * n + c] * x[c];
        }
        x[row] = s / a[row * n + row];
    }
    Some(x)
}

/// 2-norm condition number of the weighted design matrix, i.e.
/// `sqrt(lambda_max / lambda_min)` of the symmetric normal matrix, obtained
/// with a cyclic Jacobi sweep. For 3-10 columns this is cheaper and steadier
/// than a general SVD, and needs no linear-algebra crate.
fn condition_number(ata: &[f64], p: usize) -> f64 {
    let mut a = ata.to_vec();
    for _ in 0..60 {
        // Largest off-diagonal magnitude.
        let mut off = 0.0;
        let (mut pi, mut qi) = (0usize, 1usize);
        for i in 0..p {
            for j in (i + 1)..p {
                let v = a[i * p + j].abs();
                if v > off {
                    off = v;
                    pi = i;
                    qi = j;
                }
            }
        }
        if off < 1e-14 {
            break;
        }
        let (app, aqq, apq) = (a[pi * p + pi], a[qi * p + qi], a[pi * p + qi]);
        let theta = 0.5 * (2.0 * apq).atan2(app - aqq);
        let (c, s) = (theta.cos(), theta.sin());
        for k in 0..p {
            let (akp, akq) = (a[k * p + pi], a[k * p + qi]);
            a[k * p + pi] = c * akp + s * akq;
            a[k * p + qi] = -s * akp + c * akq;
        }
        for k in 0..p {
            let (apk, aqk) = (a[pi * p + k], a[qi * p + k]);
            a[pi * p + k] = c * apk + s * aqk;
            a[qi * p + k] = -s * apk + c * aqk;
        }
    }
    let eig: Vec<f64> = (0..p).map(|i| a[i * p + i].abs()).collect();
    let hi = eig.iter().cloned().fold(0.0_f64, f64::max);
    let lo = eig.iter().cloned().fold(f64::INFINITY, f64::min);
    if lo <= 0.0 || !lo.is_finite() {
        return f64::INFINITY;
    }
    (hi / lo).sqrt()
}

fn bounds(v: &[f64]) -> (f64, f64) {
    v.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &x| {
        (lo.min(x), hi.max(x))
    })
}

// ── Params ──────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kernel {
    Exponential,
    Gaussian,
    Quartic,
    Epanechnikov,
    FifthOrder,
    Constant,
}

impl Kernel {
    fn weight(self, d: f64, h: f64) -> f64 {
        if h <= 0.0 {
            return 0.0;
        }
        let t = d / h;
        match self {
            Kernel::Constant => 1.0,
            Kernel::Exponential => (-t).exp(),
            Kernel::Gaussian => (-0.5 * t * t).exp(),
            Kernel::Quartic => {
                if t < 1.0 {
                    let u = 1.0 - t * t;
                    u * u
                } else {
                    0.0
                }
            }
            Kernel::Epanechnikov => {
                if t < 1.0 {
                    1.0 - t * t
                } else {
                    0.0
                }
            }
            Kernel::FifthOrder => {
                if t < 1.0 {
                    let u = 1.0 - t * t;
                    u * u * u
                } else {
                    0.0
                }
            }
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum OutputType {
    Prediction,
    StandardError,
    ConditionNumber,
}

impl OutputType {
    fn name(self) -> &'static str {
        match self {
            OutputType::Prediction => "prediction",
            OutputType::StandardError => "standard_error",
            OutputType::ConditionNumber => "condition_number",
        }
    }
}

fn parse_order(args: &ToolArgs) -> Result<usize, ToolError> {
    match parse_optional_f64(args, "order")? {
        None => Ok(1),
        Some(v) if (v - 1.0).abs() < 1e-9 => Ok(1),
        Some(v) if (v - 2.0).abs() < 1e-9 => Ok(2),
        Some(v) if (v - 3.0).abs() < 1e-9 => Ok(3),
        Some(v) => Err(ToolError::Validation(format!(
            "'order' must be 1, 2 or 3, got {v}"
        ))),
    }
}

fn parse_kernel(args: &ToolArgs) -> Result<Kernel, ToolError> {
    match args
        .get("kernel")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        None | Some("") | Some("exponential") => Ok(Kernel::Exponential),
        Some("gaussian") => Ok(Kernel::Gaussian),
        Some("quartic") => Ok(Kernel::Quartic),
        Some("epanechnikov") => Ok(Kernel::Epanechnikov),
        Some("fifth_order") => Ok(Kernel::FifthOrder),
        Some("constant") => Ok(Kernel::Constant),
        Some(o) => Err(ToolError::Validation(format!(
            "'kernel' must be exponential/gaussian/quartic/epanechnikov/fifth_order/constant, got '{o}'"
        ))),
    }
}

fn parse_output_type(args: &ToolArgs) -> Result<OutputType, ToolError> {
    match args
        .get("output_type")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        None | Some("") | Some("prediction") => Ok(OutputType::Prediction),
        Some("standard_error") => Ok(OutputType::StandardError),
        Some("condition_number") => Ok(OutputType::ConditionNumber),
        Some(o) => Err(ToolError::Validation(format!(
            "'output_type' must be 'prediction', 'standard_error' or 'condition_number', got '{o}'"
        ))),
    }
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

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Samples on a 9x9 lattice over (0,0)-(80,80), valued by `f(x, y)`.
    fn samples(f: impl Fn(f64, f64) -> f64) -> String {
        let mut l = Layer::new("s")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("z", FieldType::Float));
        for i in 0..9 {
            for j in 0..9 {
                let (x, y) = (i as f64 * 10.0, j as f64 * 10.0);
                l.add_feature(
                    Some(Geometry::point(x, y)),
                    &[("z", FieldValue::Float(f(x, y)))],
                )
                .unwrap();
            }
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = LocalPolynomialInterpolationTool.run(&args, &ctx()).unwrap();
        let r = crate::common::load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, r)
    }

    fn stats(r: &Raster) -> (f64, f64, usize) {
        let (mut lo, mut hi, mut n) = (f64::INFINITY, f64::NEG_INFINITY, 0usize);
        for row in 0..r.rows {
            for col in 0..r.cols {
                let v = r.get(0, row as isize, col as isize);
                if v != r.nodata {
                    lo = lo.min(v);
                    hi = hi.max(v);
                    n += 1;
                }
            }
        }
        (lo, hi, n)
    }

    #[test]
    fn order_one_reproduces_a_planar_field_exactly() {
        // z = 2x + 3y + 5 is in the order-1 model space, so a local linear fit
        // must reproduce it to machine precision.
        let input = samples(|x, y| 2.0 * x + 3.0 * y + 5.0);
        let (out, r) = run(json!({
            "input": input, "z_field": "z", "cell_size": 10.0, "order": 1
        }));
        assert!(out.outputs["fitted_cells"].as_f64().unwrap() > 0.0);
        for row in 0..r.rows {
            for col in 0..r.cols {
                let v = r.get(0, row as isize, col as isize);
                if v == r.nodata {
                    continue;
                }
                let x = r.x_min + (col as f64 + 0.5) * r.cell_size_x;
                let y = r.y_max() - (row as f64 + 0.5) * r.cell_size_y;
                let want = 2.0 * x + 3.0 * y + 5.0;
                assert!((v - want).abs() < 1e-6, "at ({x},{y}) got {v} want {want}");
            }
        }
    }

    #[test]
    fn order_two_beats_order_one_on_a_quadratic_field() {
        // The distinguishing property vs the bundled global trend_surface and
        // vs a linear fit: curvature is captured.
        let f = |x: f64, y: f64| 0.01 * x * x + 0.02 * y * y;
        let input = samples(f);
        let err = |order: i32| -> f64 {
            let (_o, r) = run(json!({
                "input": input.clone(), "z_field": "z", "cell_size": 10.0,
                "order": order, "kernel": "gaussian"
            }));
            let mut worst = 0.0_f64;
            for row in 0..r.rows {
                for col in 0..r.cols {
                    let v = r.get(0, row as isize, col as isize);
                    if v == r.nodata {
                        continue;
                    }
                    let x = r.x_min + (col as f64 + 0.5) * r.cell_size_x;
                    let y = r.y_max() - (row as f64 + 0.5) * r.cell_size_y;
                    worst = worst.max((v - f(x, y)).abs());
                }
            }
            worst
        };
        let (e1, e2) = (err(1), err(2));
        assert!(e2 < e1, "order-2 error {e2} not better than order-1 {e1}");
        assert!(e2 < 1e-6, "order-2 should be near-exact, got {e2}");
    }

    #[test]
    fn constant_kernel_order_one_is_a_local_plane_fit() {
        let input = samples(|x, _y| x);
        let (out, r) = run(json!({
            "input": input, "z_field": "z", "cell_size": 20.0,
            "kernel": "constant", "order": 1
        }));
        assert!(out.outputs["fitted_cells"].as_f64().unwrap() > 0.0);
        let (lo, hi, _n) = stats(&r);
        assert!(lo >= -1e-6 && hi <= 80.0 + 1e-6, "range {lo}..{hi}");
    }

    #[test]
    fn standard_error_surface_is_finite_and_non_negative() {
        // Noise-free input: the local fit is near-perfect, so the standard
        // error should be tiny but well-defined.
        let input = samples(|x, y| x + y);
        let (out, r) = run(json!({
            "input": input, "z_field": "z", "cell_size": 20.0,
            "output_type": "standard_error", "neighbors": 12
        }));
        assert_eq!(out.outputs["output_type"], json!("standard_error"));
        let (lo, hi, n) = stats(&r);
        assert!(n > 0);
        assert!(lo >= 0.0, "negative standard error {lo}");
        assert!(hi.is_finite());
    }

    #[test]
    fn condition_number_surface_is_at_least_one() {
        let input = samples(|x, y| x + y);
        let (_o, r) = run(json!({
            "input": input, "z_field": "z", "cell_size": 20.0,
            "output_type": "condition_number", "neighbors": 12
        }));
        let (lo, _hi, n) = stats(&r);
        assert!(n > 0);
        // A condition number is a ratio of the largest to smallest singular
        // value, so it can never be below 1.
        assert!(lo >= 1.0 - 1e-9, "condition number {lo} below 1");
    }

    #[test]
    fn condition_limit_rejects_cells() {
        let input = samples(|x, y| x + y);
        let (out, _r) = run(json!({
            "input": input, "z_field": "z", "cell_size": 20.0,
            "order": 3, "condition_number": 1.0
        }));
        // An order-3 design essentially never has condition 1, so everything
        // should be refused rather than silently predicted.
        assert!(out.outputs["rejected_cells"].as_f64().unwrap() > 0.0);
        assert_eq!(out.outputs["fitted_cells"], json!(0));
    }

    #[test]
    fn too_few_samples_for_the_order_is_rejected() {
        let mut l = Layer::new("s")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("z", FieldType::Float));
        for i in 0..4 {
            l.add_feature(
                Some(Geometry::point(i as f64, 0.0)),
                &[("z", FieldValue::Float(i as f64))],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        let input = wbvector::memory_store::make_vector_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input, "z_field": "z", "order": 2
        }))
        .unwrap();
        assert!(LocalPolynomialInterpolationTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            LocalPolynomialInterpolationTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "p.shp" })).is_err());
        assert!(bad(json!({ "input": "p.shp", "z_field": "z", "order": 4 })).is_err());
        assert!(bad(json!({ "input": "p.shp", "z_field": "z", "kernel": "tricube" })).is_err());
        assert!(bad(json!({ "input": "p.shp", "z_field": "z", "bandwidth": 0 })).is_err());
        assert!(bad(json!({ "input": "p.shp", "z_field": "z", "neighbors": 2 })).is_err());
        assert!(bad(json!({ "input": "p.shp", "z_field": "z" })).is_ok());
    }
}
