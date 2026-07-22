//! GeoLibre tool: classify the local relationship between two variables.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Local Bivariate Relationships*
//! (Spatial Statistics). It complements the global `bivariate_spatial_association`
//! (Lee's L): rather than one number for the whole map, it labels — for every
//! feature and its local neighbourhood — *how* two variables relate. No bundled
//! tool measures local bivariate form.
//!
//! For each feature the tool gathers its `neighbors` nearest features (plus the
//! feature itself), forming a small local sample of `(x, y)` pairs of the two
//! variable fields. It fits a local **linear** model `y = a + b·x` and a local
//! **quadratic** model `y = a + b·x + c·x²`, and measures dependence with an
//! entropy-reduction statistic: for Gaussian residuals the differential entropy
//! of `y` drops by `ΔH = -½·ln(1 - R²)` nats once `x` is known, so a strong local
//! relationship yields a large `ΔH`. Significance is a **conditional permutation
//! test** — the local `y` values are reshuffled against fixed `x` values with a
//! seeded splitmix64 RNG (reproducible in WASM; no `Date::now`/`rand`) and the
//! statistic recomputed, giving a one-sided pseudo p-value.
//!
//! Significant features are classified from the fitted models into one of ArcGIS's
//! categories: `positive linear`, `negative linear`, `concave`, `convex`, or
//! `undefined` (significant dependence the polynomial forms do not explain).
//! Everything else is `not significant`. Output copies the input and adds
//! `lbr_type`, `lbr_stat` (ΔH), `p_value`, and `lbr_r2` fields. A projected CRS is
//! recommended so neighbour distances are metric.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct LocalBivariateRelationshipsTool;

impl Tool for LocalBivariateRelationshipsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "local_bivariate_relationships",
            display_name: "Local Bivariate Relationships",
            summary: "For each feature, classify how two variables relate in its local neighbourhood (not significant, positive/negative linear, concave, convex, undefined) from local linear and quadratic fits, an entropy-reduction dependence statistic, and a seeded permutation-test p-value — like ArcGIS Local Bivariate Relationships.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer holding the two variable fields.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (a copy of the input with relationship fields). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "field1",
                    description: "First (explanatory) variable field name.",
                    required: true,
                },
                ToolParamSpec {
                    name: "field2",
                    description: "Second (dependent) variable field name.",
                    required: true,
                },
                ToolParamSpec {
                    name: "neighbors",
                    description: "Number of nearest neighbours k defining each local sample (default 12; must be >= 5).",
                    required: false,
                },
                ToolParamSpec {
                    name: "permutations",
                    description: "Number of permutations for the significance test (default 199; 0 disables).",
                    required: false,
                },
                ToolParamSpec {
                    name: "significance",
                    description: "Pseudo p-value threshold for calling a relationship significant (default 0.05).",
                    required: false,
                },
                ToolParamSpec {
                    name: "seed",
                    description: "Seed for the permutation RNG (default 1), for reproducible results.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "field1")?;
        require_str(args, "field2")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let idx1 = layer
            .schema
            .field_index(&prm.field1)
            .ok_or_else(|| ToolError::Validation(format!("field1 '{}' not found", prm.field1)))?;
        let idx2 = layer
            .schema
            .field_index(&prm.field2)
            .ok_or_else(|| ToolError::Validation(format!("field2 '{}' not found", prm.field2)))?;

        let geographic = layer.crs_epsg().map(|e| e == 4326).unwrap_or(true);

        // Collect usable points: geometry present and both variables numeric.
        let mut pts: Vec<(f64, f64)> = Vec::new(); // coordinates
        let mut xv: Vec<f64> = Vec::new();
        let mut yv: Vec<f64> = Vec::new();
        let mut feat_of: Vec<usize> = Vec::new();
        for (fi, feature) in layer.features.iter().enumerate() {
            let Some((px, py)) = feature.geometry.as_ref().and_then(point_xy) else {
                continue;
            };
            let (Some(x), Some(y)) = (
                feature.attributes.get(idx1).and_then(FieldValue::as_f64),
                feature.attributes.get(idx2).and_then(FieldValue::as_f64),
            ) else {
                continue;
            };
            if !x.is_finite() || !y.is_finite() {
                continue;
            }
            pts.push((px, py));
            xv.push(x);
            yv.push(y);
            feat_of.push(fi);
        }

        let n = pts.len();
        let k = prm.neighbors.min(n.saturating_sub(1));
        if n < 6 || k < 5 {
            return Err(ToolError::Execution(format!(
                "need at least 6 usable points and k >= 5 (have {n} points, k={k})"
            )));
        }

        ctx.progress.info(&format!(
            "{n} usable point(s); k={k}, {} permutation(s)",
            prm.permutations
        ));

        // Per feature, gather local sample (itself + k nearest) and classify.
        let mut per_feat: Vec<(f64, f64, f64, &'static str)> =
            vec![(f64::NAN, f64::NAN, f64::NAN, ""); layer.len()];
        let mut counts: BTreeMap<&'static str, usize> = BTreeMap::new();

        for i in 0..n {
            // k nearest neighbours of i (excluding itself), then include i.
            let mut ds: Vec<(f64, usize)> = (0..n)
                .filter(|&j| j != i)
                .map(|j| (dist(pts[i], pts[j], geographic), j))
                .collect();
            ds.sort_by(|a, b| a.0.total_cmp(&b.0));
            ds.truncate(k);

            let mut lx: Vec<f64> = Vec::with_capacity(k + 1);
            let mut ly: Vec<f64> = Vec::with_capacity(k + 1);
            lx.push(xv[i]);
            ly.push(yv[i]);
            for &(_, j) in &ds {
                lx.push(xv[j]);
                ly.push(yv[j]);
            }

            let (stat, p, r2, class) = classify_local(&lx, &mut ly.clone(), &prm);
            *counts.entry(class).or_insert(0) += 1;
            per_feat[feat_of[i]] = (stat, p, r2, class);
        }

        // Write results onto a copy of the input.
        layer.add_field(FieldDef::new("lbr_type", FieldType::Text));
        layer.add_field(FieldDef::new("lbr_stat", FieldType::Float));
        layer.add_field(FieldDef::new("p_value", FieldType::Float));
        layer.add_field(FieldDef::new("lbr_r2", FieldType::Float));
        for (fi, feature) in layer.features.iter_mut().enumerate() {
            let (stat, p, r2, class) = per_feat[fi];
            feature.attributes.push(FieldValue::Text(class.to_string()));
            feature.attributes.push(FieldValue::Float(stat));
            feature.attributes.push(FieldValue::Float(p));
            feature.attributes.push(FieldValue::Float(r2));
        }

        let significant = n
            - counts.get("not significant").copied().unwrap_or(0)
            - counts.get("").copied().unwrap_or(0);
        ctx.progress.info(&format!(
            "{significant} of {n} feature(s) significant; {} positive linear, {} negative linear, {} concave, {} convex",
            counts.get("positive linear").copied().unwrap_or(0),
            counts.get("negative linear").copied().unwrap_or(0),
            counts.get("concave").copied().unwrap_or(0),
            counts.get("convex").copied().unwrap_or(0),
        ));

        let feature_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("analyzed_count".to_string(), json!(n));
        outputs.insert("significant_count".to_string(), json!(significant));
        for cat in [
            "not significant",
            "positive linear",
            "negative linear",
            "concave",
            "convex",
            "undefined",
        ] {
            let key = format!("{}_count", cat.replace(' ', "_"));
            outputs.insert(key, json!(counts.get(cat).copied().unwrap_or(0)));
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── Local classification ──────────────────────────────────────────────────────

/// Fits local linear/quadratic models to `(x, y)`, runs the permutation test and
/// returns `(stat, p_value, r2_quad, class)`.
///
/// `stat` is the entropy-reduction `ΔH = -½·ln(1 - R²_quad)` in nats. `y` is
/// passed as a scratch buffer that the permutation loop shuffles in place.
fn classify_local(x: &[f64], y: &mut [f64], prm: &Params) -> (f64, f64, f64, &'static str) {
    let m = x.len();
    let ss_tot = variance_ss(y);
    let x_var_ss = variance_ss(x);
    // Degenerate local sample: no spread in x or y -> cannot fit.
    if m < 5 || ss_tot <= 1e-12 || x_var_ss <= 1e-12 {
        return (f64::NAN, f64::NAN, f64::NAN, "undefined");
    }

    let lin = fit_linear(x, y);
    let quad = fit_quadratic(x, y);
    let r2_quad = quad.r2.max(lin.r2); // quadratic nests linear; guard numerics
    let obs_stat = entropy_reduction(r2_quad);

    // Conditional permutation test: reshuffle y against fixed x, recompute stat.
    let mut p = f64::NAN;
    if prm.permutations > 0 {
        let mut ge = 1usize; // count sim >= obs (incl. the observed arrangement)
        let mut rng = prm.seed;
        for _ in 0..prm.permutations {
            fisher_yates(y, &mut rng);
            let sim = fit_quadratic(x, y).r2.max(fit_linear(x, y).r2);
            if entropy_reduction(sim) >= obs_stat - 1e-12 {
                ge += 1;
            }
        }
        p = ge as f64 / (prm.permutations as f64 + 1.0);
    }

    // Classify. Without a permutation test, fall back to an R² floor.
    let significant = if prm.permutations > 0 {
        p <= prm.significance
    } else {
        r2_quad >= 0.5
    };
    if !significant {
        return (obs_stat, p, r2_quad, "not significant");
    }

    // Significant but poorly explained by either polynomial form -> undefined.
    if r2_quad < 0.2 {
        return (obs_stat, p, r2_quad, "undefined");
    }

    // Does the quadratic term add materially over the linear fit?
    let improvement = quad.r2 - lin.r2;
    let curved = improvement > 0.10 && quad.curvature.abs() > 1e-12;

    let class = if curved {
        // Vertex of the quadratic (x centred): does the curve turn inside the
        // observed x-range? If so it is genuinely concave/convex; otherwise it is
        // monotonic over the data and better described as linear.
        let vertex = -quad.slope / (2.0 * quad.curvature);
        let (xc_min, xc_max) = centred_range(x);
        if vertex > xc_min && vertex < xc_max {
            if quad.curvature > 0.0 {
                "convex" // opens upward (U-shaped)
            } else {
                "concave" // opens downward (∩-shaped)
            }
        } else if lin.slope >= 0.0 {
            "positive linear"
        } else {
            "negative linear"
        }
    } else if lin.slope >= 0.0 {
        "positive linear"
    } else {
        "negative linear"
    };

    (obs_stat, p, r2_quad, class)
}

/// Entropy reduction (nats) of a Gaussian residual model with coefficient of
/// determination `r2`: `ΔH = -½·ln(1 - r2)`, clamped for numerical safety.
fn entropy_reduction(r2: f64) -> f64 {
    let r2 = r2.clamp(0.0, 0.999_999);
    -0.5 * (1.0 - r2).ln()
}

struct LinFit {
    slope: f64,
    r2: f64,
}

/// Ordinary least squares `y = a + b·x`.
fn fit_linear(x: &[f64], y: &[f64]) -> LinFit {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut syy = 0.0;
    for (&xi, &yi) in x.iter().zip(y) {
        let dx = xi - mx;
        let dy = yi - my;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    if sxx <= 1e-12 || syy <= 1e-12 {
        return LinFit {
            slope: 0.0,
            r2: 0.0,
        };
    }
    let slope = sxy / sxx;
    let r2 = (sxy * sxy) / (sxx * syy);
    LinFit {
        slope,
        r2: r2.clamp(0.0, 1.0),
    }
}

struct QuadFit {
    slope: f64,     // linear coefficient b (x centred)
    curvature: f64, // quadratic coefficient c
    r2: f64,
}

/// Least squares `y = a + b·xc + c·xc²` with `xc = x - mean(x)` (centring keeps
/// the normal equations well-conditioned; it does not change the fit or R²).
fn fit_quadratic(x: &[f64], y: &[f64]) -> QuadFit {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    // Accumulate the normal-equation moments for design [1, xc, xc²].
    let (mut s1, mut s2, mut s3, mut s4) = (0.0, 0.0, 0.0, 0.0); // sums of xc^1..4
    let (mut sy, mut sxy, mut sx2y) = (0.0, 0.0, 0.0);
    for (&xi, &yi) in x.iter().zip(y) {
        let xc = xi - mx;
        let xc2 = xc * xc;
        s1 += xc;
        s2 += xc2;
        s3 += xc2 * xc;
        s4 += xc2 * xc2;
        sy += yi;
        sxy += xc * yi;
        sx2y += xc2 * yi;
    }
    // Symmetric 3x3 system: [[n, s1, s2],[s1, s2, s3],[s2, s3, s4]] · [a,b,c] = [sy, sxy, sx2y].
    let a = [[n, s1, s2], [s1, s2, s3], [s2, s3, s4]];
    let rhs = [sy, sxy, sx2y];
    let Some(sol) = solve3(a, rhs) else {
        // Singular (e.g. degenerate x) -> defer to the linear fit.
        let lin = fit_linear(x, y);
        return QuadFit {
            slope: lin.slope,
            curvature: 0.0,
            r2: lin.r2,
        };
    };
    let (_a, b, c) = (sol[0], sol[1], sol[2]);
    // R² from residuals.
    let my = sy / n;
    let mut ss_res = 0.0;
    let mut ss_tot = 0.0;
    for (&xi, &yi) in x.iter().zip(y) {
        let xc = xi - mx;
        let pred = sol[0] + b * xc + c * xc * xc;
        ss_res += (yi - pred).powi(2);
        ss_tot += (yi - my).powi(2);
    }
    let r2 = if ss_tot <= 1e-12 {
        0.0
    } else {
        (1.0 - ss_res / ss_tot).clamp(0.0, 1.0)
    };
    QuadFit {
        slope: b,
        curvature: c,
        r2,
    }
}

/// Solves a 3x3 linear system by Gaussian elimination with partial pivoting.
fn solve3(mut a: [[f64; 3]; 3], mut b: [f64; 3]) -> Option<[f64; 3]> {
    for col in 0..3 {
        // Partial pivot.
        let mut piv = col;
        for r in (col + 1)..3 {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-12 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        // Eliminate below.
        let pivot_row = a[col];
        let pivot_b = b[col];
        for r in (col + 1)..3 {
            let f = a[r][col] / pivot_row[col];
            a[r].iter_mut()
                .zip(pivot_row.iter())
                .for_each(|(ar, pr)| *ar -= f * pr);
            b[r] -= f * pivot_b;
        }
    }
    // Back-substitution.
    let mut x = [0.0; 3];
    for i in (0..3).rev() {
        let mut s = b[i];
        for j in (i + 1)..3 {
            s -= a[i][j] * x[j];
        }
        x[i] = s / a[i][i];
    }
    Some(x)
}

/// Total sum of squares (variance × n) of a slice.
fn variance_ss(v: &[f64]) -> f64 {
    let n = v.len() as f64;
    if n == 0.0 {
        return 0.0;
    }
    let m = v.iter().sum::<f64>() / n;
    v.iter().map(|&z| (z - m).powi(2)).sum()
}

/// Range of `x` after centring on its mean.
fn centred_range(x: &[f64]) -> (f64, f64) {
    let n = x.len() as f64;
    let m = x.iter().sum::<f64>() / n;
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &xi in x {
        let xc = xi - m;
        lo = lo.min(xc);
        hi = hi.max(xc);
    }
    (lo, hi)
}

// ── Geometry / RNG helpers ────────────────────────────────────────────────────

fn dist(a: (f64, f64), b: (f64, f64), geographic: bool) -> f64 {
    if geographic {
        haversine(a.1, a.0, b.1, b.0)
    } else {
        (a.0 - b.0).hypot(a.1 - b.1)
    }
}

fn haversine(lat0: f64, lon0: f64, lat1: f64, lon1: f64) -> f64 {
    const R: f64 = 6_371_000.0;
    let (p0, p1) = (lat0.to_radians(), lat1.to_radians());
    let dphi = (lat1 - lat0).to_radians();
    let dlmb = (lon1 - lon0).to_radians();
    let a = (dphi / 2.0).sin().powi(2) + p0.cos() * p1.cos() * (dlmb / 2.0).sin().powi(2);
    2.0 * R * a.sqrt().asin()
}

/// In-place Fisher-Yates shuffle with a seeded splitmix64 RNG.
fn fisher_yates<T>(v: &mut [T], rng: &mut u64) {
    let n = v.len();
    for i in (1..n).rev() {
        let j = (next_u64(rng) % (i as u64 + 1)) as usize;
        v.swap(i, j);
    }
}

fn next_u64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => {
            let n = cs.len() as f64;
            Some((
                cs.iter().map(|c| c.x).sum::<f64>() / n,
                cs.iter().map(|c| c.y).sum::<f64>() / n,
            ))
        }
        _ => None,
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    field1: String,
    field2: String,
    neighbors: usize,
    permutations: usize,
    significance: f64,
    seed: u64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let field1 = require_str(args, "field1")?.to_string();
    let field2 = require_str(args, "field2")?.to_string();
    let neighbors = match parse_opt_u64(args, "neighbors")? {
        None => 12,
        Some(v) if v >= 5 => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "'neighbors' must be >= 5".to_string(),
            ))
        }
    };
    let permutations = parse_opt_u64(args, "permutations")?.unwrap_or(199) as usize;
    let significance = match opt_f64(args, "significance")? {
        None => 0.05,
        Some(v) if v > 0.0 && v < 1.0 => v,
        Some(_) => {
            return Err(ToolError::Validation(
                "'significance' must be between 0 and 1".to_string(),
            ))
        }
    };
    let seed = parse_opt_u64(args, "seed")?.unwrap_or(1);
    Ok(Params {
        field1,
        field2,
        neighbors,
        permutations,
        significance,
        seed,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

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

fn parse_opt_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be an integer"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be an integer"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Point layer with (x, y-coord, v1, v2) rows in a projected CRS.
    fn layer_of(rows: &[(f64, f64, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v1", FieldType::Float));
        l.add_field(FieldDef::new("v2", FieldType::Float));
        for &(x, y, v1, v2) in rows {
            l.add_feature(
                Some(Geometry::point(x, y)),
                &[("v1", v1.into()), ("v2", v2.into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = LocalBivariateRelationshipsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// A tight cluster with a strong positive-linear v1->v2 relationship is
    /// classified "positive linear" with a low p-value.
    #[test]
    fn detects_positive_linear() {
        // 20 points on a small grid; v2 = 2*v1 + small deterministic wobble.
        let mut rows = Vec::new();
        for i in 0..20 {
            let x = (i % 5) as f64;
            let y = (i / 5) as f64;
            let v1 = i as f64;
            let wobble = if i % 2 == 0 { 0.3 } else { -0.3 };
            rows.push((x, y, v1, 2.0 * v1 + wobble));
        }
        let (out, layer) = run(json!({
            "input": layer_of(&rows), "field1": "v1", "field2": "v2",
            "neighbors": 8, "permutations": 199, "seed": 7,
        }));
        // The great majority of features should be positive linear.
        let pos = out.outputs["positive_linear_count"].as_u64().unwrap();
        assert!(pos >= 15, "expected mostly positive linear, got {pos}");
        // Every significant feature carries a low p-value.
        let pidx = layer.schema.field_index("p_value").unwrap();
        let tidx = layer.schema.field_index("lbr_type").unwrap();
        for f in &layer.features {
            if f.attributes[tidx].as_str() == Some("positive linear") {
                assert!(f.attributes[pidx].as_f64().unwrap() <= 0.05);
            }
        }
    }

    /// A U-shaped (convex) relationship is detected as convex.
    #[test]
    fn detects_convex() {
        // v2 = (v1 - 5)^2 over a compact cluster.
        let mut rows = Vec::new();
        for i in 0..24 {
            let x = (i % 6) as f64;
            let y = (i / 6) as f64;
            let v1 = (i % 11) as f64; // spans 0..10 so the vertex (5) is inside
            let wob = if i % 2 == 0 { 0.2 } else { -0.2 };
            rows.push((x, y, v1, (v1 - 5.0).powi(2) + wob));
        }
        let (out, _l) = run(json!({
            "input": layer_of(&rows), "field1": "v1", "field2": "v2",
            "neighbors": 12, "permutations": 199, "seed": 3,
        }));
        let convex = out.outputs["convex_count"].as_u64().unwrap();
        assert!(
            convex >= 1,
            "expected some convex classifications, got {convex}"
        );
    }

    /// Pure noise -> predominantly "not significant".
    #[test]
    fn noise_is_not_significant() {
        // Deterministic pseudo-random v1/v2 with no relationship.
        let mut rng = 0x1234_5678u64;
        let mut rows = Vec::new();
        for i in 0..40 {
            let x = (i % 8) as f64;
            let y = (i / 8) as f64;
            let v1 = (next_u64(&mut rng) % 1000) as f64;
            let v2 = (next_u64(&mut rng) % 1000) as f64;
            rows.push((x, y, v1, v2));
        }
        let (out, _l) = run(json!({
            "input": layer_of(&rows), "field1": "v1", "field2": "v2",
            "neighbors": 10, "permutations": 199, "seed": 11,
        }));
        let ns = out.outputs["not_significant_count"].as_u64().unwrap();
        let n = out.outputs["analyzed_count"].as_u64().unwrap();
        assert!(
            ns * 2 >= n,
            "noise should be mostly not significant: {ns}/{n}"
        );
    }

    /// Deterministic: identical seed -> identical per-feature output.
    #[test]
    fn deterministic_with_seed() {
        let mut rows = Vec::new();
        for i in 0..24 {
            let x = (i % 6) as f64;
            let y = (i / 6) as f64;
            let v1 = i as f64;
            rows.push((x, y, v1, 1.5 * v1 + (i % 3) as f64));
        }
        let args = json!({
            "input": layer_of(&rows), "field1": "v1", "field2": "v2",
            "neighbors": 8, "permutations": 149, "seed": 99,
        });
        let (_o1, l1) = run(args.clone());
        let (_o2, l2) = run(args);
        let ti = l1.schema.field_index("lbr_type").unwrap();
        let pi = l1.schema.field_index("p_value").unwrap();
        for (a, b) in l1.features.iter().zip(&l2.features) {
            assert_eq!(a.attributes[ti].as_str(), b.attributes[ti].as_str());
            assert_eq!(
                a.attributes[pi].as_f64().unwrap().to_bits(),
                b.attributes[pi].as_f64().unwrap().to_bits()
            );
        }
    }

    /// Non-point geometry / missing values pass through as empty class.
    #[test]
    fn skips_features_without_values() {
        // One feature has a null v2; it should not be analyzed but still present.
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v1", FieldType::Float));
        l.add_field(FieldDef::new("v2", FieldType::Float));
        for i in 0..10 {
            l.add_feature(
                Some(Geometry::point(i as f64, 0.0)),
                &[("v1", (i as f64).into()), ("v2", (2.0 * i as f64).into())],
            )
            .unwrap();
        }
        // A feature with a null v2.
        l.push(wbvector::Feature {
            fid: 0,
            geometry: Some(Geometry::point(100.0, 100.0)),
            attributes: vec![FieldValue::Float(1.0), FieldValue::Null],
        });
        let id = memory_store::put_vector(l);
        let path = memory_store::make_vector_memory_path(&id);
        let (out, layer) = run(json!({
            "input": path, "field1": "v1", "field2": "v2",
            "neighbors": 5, "permutations": 49, "seed": 1,
        }));
        assert_eq!(out.outputs["feature_count"], json!(11));
        assert_eq!(out.outputs["analyzed_count"], json!(10));
        // The null-value feature keeps an empty class.
        let ti = layer.schema.field_index("lbr_type").unwrap();
        assert_eq!(layer.features[10].attributes[ti].as_str(), Some(""));
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            LocalBivariateRelationshipsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "field1": "v1" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "field1": "v1", "field2": "v2" })).is_ok());
        // neighbors < 5 rejected.
        assert!(bad(
            json!({ "input": "a.geojson", "field1": "v1", "field2": "v2", "neighbors": 3 })
        )
        .is_err());
        // significance out of range rejected.
        assert!(bad(
            json!({ "input": "a.geojson", "field1": "v1", "field2": "v2", "significance": 1.5 })
        )
        .is_err());
    }
}
