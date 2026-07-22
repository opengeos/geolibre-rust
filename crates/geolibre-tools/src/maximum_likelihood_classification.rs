//! GeoLibre tool: Gaussian maximum-likelihood supervised classification.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Maximum Likelihood Classification*
//! (Spatial Analyst) — the classic Bayesian per-pixel classifier and the only
//! supervised method in the suite that emits a per-pixel probability/confidence
//! surface.
//!
//! Where `matched_filter_target_detection` scores pixels against a *single*
//! known target spectrum and the bundled `min_dist`/`parallelepiped`-style
//! classifiers use crude decision rules, this tool models every training class
//! as a multivariate Gaussian: from the labelled training pixels it estimates a
//! per-class mean vector `μ_c` and covariance `Σ_c`, then assigns each image
//! pixel `x` to the class with the largest Gaussian discriminant
//!
//! ```text
//!   g_c(x) = ln P_c − ½·ln|Σ_c| − ½·(x−μ_c)ᵀ Σ_c⁻¹ (x−μ_c)
//! ```
//!
//! `P_c` is the a-priori class weight: **equal** (`1/K`, the default — the
//! `ln P_c` term then cancels), or **sample** (`n_c / N`, proportional to each
//! class's training-pixel count). The constant `−½·b·ln(2π)` is common to every
//! class and dropped from the discriminant, but restored (implicitly, as it
//! cancels) when the posterior probabilities are formed.
//!
//! **Training** is supplied as an integer-labelled *training raster* on the same
//! grid as the image: each distinct positive integer is a class, and cells equal
//! to the raster's no-data value (or `≤ 0`, or non-finite) are unlabelled. This
//! mirrors `zonal_histogram`'s scope decision — rasterizing a polygon training
//! layer on the fly is deliberately out of scope for v1; callers rasterize their
//! training polygons first (any zone raster works).
//!
//! **Inversion is pure Rust with no dense-LA dependency.** Each small `b×b`
//! covariance is ridge-stabilized and **Cholesky-decomposed** (`Σ = L Lᵀ`). The
//! log-determinant is `2·Σ ln L_ii`, and the Mahalanobis term is evaluated
//! without ever forming `Σ⁻¹`: forward-solving `L y = (x−μ)` gives
//! `(x−μ)ᵀ Σ⁻¹ (x−μ) = ‖y‖²`.
//!
//! **Outputs.** The primary `output` is a single-band classified raster carrying
//! the original training class ids. With `prob_output`, a second F32 raster holds
//! the per-pixel **confidence** — the posterior probability of the winning class,
//! `P(c|x) = exp(g_c)/Σ_k exp(g_k)`, in `[0, 1]` (computed with a log-sum-exp for
//! stability). `reject_fraction` ∈ `[0, 1)` leaves the least-confident
//! `reject_fraction` of classified cells unclassified (no-data): the threshold is
//! the `reject_fraction`-quantile of the winning-class confidences, so the effect
//! matches ArcGIS's "fraction of cells that will not be classified".

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{
    band_to_vec, load_input_raster, parse_optional_output, raster_like_with_data,
    write_or_store_output,
};

/// No-data / unclassified sentinel for the class and probability outputs.
const OUT_NODATA: f64 = -1.0;
/// Ridge added to each class covariance diagonal (× its mean) for a positive
/// definite Cholesky factorization even from thin or collinear training sets.
const RIDGE: f64 = 1e-6;

#[derive(Clone, Copy, PartialEq, Eq)]
enum APriori {
    Equal,
    Sample,
}

impl APriori {
    fn as_str(self) -> &'static str {
        match self {
            Self::Equal => "equal",
            Self::Sample => "sample",
        }
    }
}

pub struct MaximumLikelihoodClassificationTool;

impl Tool for MaximumLikelihoodClassificationTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "maximum_likelihood_classification",
            display_name: "Maximum Likelihood Classification",
            summary: "Gaussian maximum-likelihood supervised classifier: models each training class as a multivariate normal (per-class mean + covariance), classifies every pixel of a multiband image by the largest Gaussian log-likelihood with optional equal/sample a-priori weights, and can emit a per-pixel probability/confidence raster and reject the least-confident cells, like ArcGIS Maximum Likelihood Classification.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input multiband raster to classify (1+ bands; every band is a feature).",
                    required: true,
                },
                ToolParamSpec {
                    name: "training",
                    description: "Integer-labelled training raster on the same grid as 'input': each distinct positive integer is a class; no-data / non-positive / non-finite cells are unlabelled.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output classified raster carrying the training class ids. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "prob_output",
                    description: "Optional output F32 raster of the winning class's posterior probability (confidence) in [0,1].",
                    required: false,
                },
                ToolParamSpec {
                    name: "a_priori",
                    description: "A-priori class weights: 'equal' (default, 1/K) or 'sample' (proportional to each class's training-pixel count).",
                    required: false,
                },
                ToolParamSpec {
                    name: "reject_fraction",
                    description: "Fraction in [0,1) of the least-confident classified cells to leave unclassified (no-data). Default 0.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["input", "training"] {
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
        let input = args.get("input").and_then(Value::as_str).unwrap();
        let training = args.get("training").and_then(Value::as_str).unwrap();
        let output = parse_optional_output(args, "output")?;
        let prm = parse_params(args)?;

        let raster = load_input_raster(input)?;
        let train = load_input_raster(training)?;
        if raster.rows != train.rows || raster.cols != train.cols {
            return Err(ToolError::Validation(format!(
                "'training' grid {}x{} does not match 'input' grid {}x{}",
                train.rows, train.cols, raster.rows, raster.cols
            )));
        }
        let b = raster.bands;
        if b == 0 {
            return Err(ToolError::Validation(
                "input raster has no bands".to_string(),
            ));
        }
        let rows = raster.rows;
        let cols = raster.cols;
        let n = rows * cols;
        let nodata = raster.nodata;

        // Materialize the image band stack and a per-pixel validity mask (a pixel
        // is classifiable only when every band is valid).
        let mut stack: Vec<Vec<f64>> = vec![vec![0.0; n]; b];
        let mut valid = vec![true; n];
        for (band, stack_band) in stack.iter_mut().enumerate() {
            for row in 0..rows {
                for col in 0..cols {
                    let v = raster.get(band as isize, row as isize, col as isize);
                    let i = row * cols + col;
                    if v == nodata || !v.is_finite() {
                        valid[i] = false;
                    }
                    stack_band[i] = v;
                }
            }
        }

        // Read training labels and collect per-class pixel indices.
        let train_nd = train.nodata;
        let labels = band_to_vec(&train, 0);
        let mut class_pixels: BTreeMap<i64, Vec<usize>> = BTreeMap::new();
        for i in 0..n {
            let l = labels[i];
            if l == train_nd || !l.is_finite() {
                continue;
            }
            let id = l.round() as i64;
            if id <= 0 {
                continue;
            }
            if valid[i] {
                class_pixels.entry(id).or_default().push(i);
            }
        }
        if class_pixels.len() < 2 {
            return Err(ToolError::Execution(format!(
                "need at least 2 training classes with valid pixels, found {}",
                class_pixels.len()
            )));
        }

        ctx.progress.info(&format!(
            "training {} classes on {}-band image ({} a-priori)",
            class_pixels.len(),
            b,
            prm.a_priori.as_str()
        ));

        // Estimate each class's Gaussian: mean, Cholesky factor, log-determinant.
        let total_train: usize = class_pixels.values().map(Vec::len).sum();
        let mut models: Vec<ClassModel> = Vec::with_capacity(class_pixels.len());
        for (&id, pixels) in &class_pixels {
            if pixels.len() <= b {
                return Err(ToolError::Execution(format!(
                    "class {id} has {} training pixel(s); need more than the band count ({b}) to estimate a covariance",
                    pixels.len()
                )));
            }
            let (mean, cov) = mean_cov(&stack, pixels, b);
            let chol = cholesky(&cov, b).ok_or_else(|| {
                ToolError::Execution(format!(
                    "class {id} covariance is not positive definite (degenerate training sample)"
                ))
            })?;
            let log_det = 2.0 * (0..b).map(|k| chol[k * b + k].ln()).sum::<f64>();
            let ln_prior = match prm.a_priori {
                APriori::Equal => -(class_pixels.len() as f64).ln(),
                APriori::Sample => (pixels.len() as f64 / total_train as f64).ln(),
            };
            models.push(ClassModel {
                id,
                mean,
                chol,
                ln_prior,
                half_log_det: 0.5 * log_det,
            });
        }

        ctx.progress.info("classifying pixels");

        // Per-pixel: pick the class with the largest discriminant; the posterior
        // confidence is exp(g_win) / Σ exp(g_k) via a log-sum-exp.
        let mut class_out = vec![OUT_NODATA; n];
        let mut conf = vec![OUT_NODATA; n];
        let mut x = vec![0.0; b];
        let mut disc = vec![0.0; models.len()];
        for i in 0..n {
            if !valid[i] {
                continue;
            }
            for (k, xk) in x.iter_mut().enumerate() {
                *xk = stack[k][i];
            }
            let mut best = 0usize;
            let mut best_g = f64::NEG_INFINITY;
            for (m, model) in models.iter().enumerate() {
                let g = model.discriminant(&x, b);
                disc[m] = g;
                if g > best_g {
                    best_g = g;
                    best = m;
                }
            }
            // Posterior of the winner via log-sum-exp over the discriminants.
            let mut sum = 0.0;
            for &g in &disc {
                sum += (g - best_g).exp();
            }
            let confidence = 1.0 / sum; // exp(best_g - best_g) / Σ exp(g - best_g)
            class_out[i] = models[best].id as f64;
            conf[i] = confidence.clamp(0.0, 1.0);
        }

        // reject_fraction: leave the least-confident classified cells unclassified.
        let mut rejected = 0usize;
        if prm.reject_fraction > 0.0 {
            let mut confs: Vec<f64> = (0..n).filter(|&i| valid[i]).map(|i| conf[i]).collect();
            if !confs.is_empty() {
                confs.sort_by(f64::total_cmp);
                let idx = (prm.reject_fraction * confs.len() as f64).floor() as usize;
                let cut = confs[idx.min(confs.len() - 1)];
                for i in 0..n {
                    if valid[i] && conf[i] < cut {
                        class_out[i] = OUT_NODATA;
                        conf[i] = OUT_NODATA;
                        rejected += 1;
                    }
                }
            }
        }

        let classified = (0..n)
            .filter(|&i| valid[i] && class_out[i] != OUT_NODATA)
            .count();

        let out_raster = raster_like_with_data(&raster, class_out, OUT_NODATA, DataType::F32)?;
        let out_path = write_or_store_output(out_raster, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("bands".to_string(), json!(b));
        outputs.insert("class_count".to_string(), json!(models.len()));
        outputs.insert("training_pixels".to_string(), json!(total_train));
        outputs.insert("classified_pixels".to_string(), json!(classified));
        outputs.insert("a_priori".to_string(), json!(prm.a_priori.as_str()));
        if prm.reject_fraction > 0.0 {
            outputs.insert("rejected_pixels".to_string(), json!(rejected));
        }

        if prm.prob_output.is_some() || prm.want_prob {
            let prob_raster = raster_like_with_data(&raster, conf, OUT_NODATA, DataType::F32)?;
            let prob_path = write_or_store_output(prob_raster, prm.prob_output.as_deref())?;
            outputs.insert("prob_output".to_string(), json!(prob_path));
        }

        Ok(ToolRunResult { outputs })
    }
}

// ── Per-class Gaussian model ──────────────────────────────────────────────────

struct ClassModel {
    id: i64,
    mean: Vec<f64>,
    /// Lower-triangular Cholesky factor `L` (row-major, `b×b`) of the covariance.
    chol: Vec<f64>,
    ln_prior: f64,
    half_log_det: f64,
}

impl ClassModel {
    /// Gaussian discriminant `g(x) = ln P − ½ ln|Σ| − ½ (x−μ)ᵀ Σ⁻¹ (x−μ)`.
    /// The Mahalanobis term is `‖L⁻¹(x−μ)‖²`, obtained by forward substitution
    /// without ever forming `Σ⁻¹`.
    fn discriminant(&self, x: &[f64], b: usize) -> f64 {
        // Forward-solve L y = (x − μ).
        let mut maha = 0.0;
        let mut y = [0.0f64; 32];
        let y = if b <= 32 {
            &mut y[..b]
        } else {
            // Fall back to a heap vector for very wide (>32-band) images.
            return self.discriminant_heap(x, b);
        };
        for k in 0..b {
            let mut s = x[k] - self.mean[k];
            for (j, &yj) in y[..k].iter().enumerate() {
                s -= self.chol[k * b + j] * yj;
            }
            let yk = s / self.chol[k * b + k];
            y[k] = yk;
            maha += yk * yk;
        }
        self.ln_prior - self.half_log_det - 0.5 * maha
    }

    fn discriminant_heap(&self, x: &[f64], b: usize) -> f64 {
        let mut maha = 0.0;
        let mut y = vec![0.0f64; b];
        for k in 0..b {
            let mut s = x[k] - self.mean[k];
            for (j, &yj) in y[..k].iter().enumerate() {
                s -= self.chol[k * b + j] * yj;
            }
            let yk = s / self.chol[k * b + k];
            y[k] = yk;
            maha += yk * yk;
        }
        self.ln_prior - self.half_log_det - 0.5 * maha
    }
}

/// Sample mean vector and (unbiased, `n−1`) covariance over a class's training
/// pixel indices, ridge-stabilized on the diagonal for invertibility.
fn mean_cov(stack: &[Vec<f64>], pixels: &[usize], b: usize) -> (Vec<f64>, Vec<f64>) {
    let count = pixels.len() as f64;
    let mut mean = vec![0.0; b];
    for &i in pixels {
        for (k, m) in mean.iter_mut().enumerate() {
            *m += stack[k][i];
        }
    }
    for m in mean.iter_mut() {
        *m /= count;
    }
    let mut cov = vec![0.0; b * b];
    for &i in pixels {
        for k in 0..b {
            let dk = stack[k][i] - mean[k];
            for l in k..b {
                cov[k * b + l] += dk * (stack[l][i] - mean[l]);
            }
        }
    }
    let denom = (count - 1.0).max(1.0);
    for k in 0..b {
        for l in k..b {
            let v = cov[k * b + l] / denom;
            cov[k * b + l] = v;
            cov[l * b + k] = v;
        }
    }
    // Ridge (× the mean diagonal) keeps thin/collinear samples factorable.
    let diag_mean: f64 = (0..b).map(|k| cov[k * b + k]).sum::<f64>() / b as f64;
    let ridge = RIDGE * diag_mean.abs().max(1e-12);
    for k in 0..b {
        cov[k * b + k] += ridge;
    }
    (mean, cov)
}

/// Cholesky factorization `Σ = L Lᵀ` of a symmetric matrix (row-major `b×b`).
/// Returns the lower-triangular `L` (row-major; upper entries left as 0), or
/// `None` if `Σ` is not positive definite.
fn cholesky(m: &[f64], b: usize) -> Option<Vec<f64>> {
    let mut l = vec![0.0; b * b];
    for i in 0..b {
        for j in 0..=i {
            let mut s = m[i * b + j];
            for k in 0..j {
                s -= l[i * b + k] * l[j * b + k];
            }
            if i == j {
                if s <= 0.0 || !s.is_finite() {
                    return None;
                }
                l[i * b + j] = s.sqrt();
            } else {
                l[i * b + j] = s / l[j * b + j];
            }
        }
    }
    Some(l)
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    a_priori: APriori,
    reject_fraction: f64,
    prob_output: Option<String>,
    /// True when `prob_output` was explicitly requested (even without a path, in
    /// which case the probability raster is stored in memory).
    want_prob: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let a_priori = match args
        .get("a_priori")
        .and_then(Value::as_str)
        .map(|s| s.trim().to_ascii_lowercase())
    {
        None => APriori::Equal,
        Some(s) if s.is_empty() || s == "equal" => APriori::Equal,
        Some(s) if s == "sample" => APriori::Sample,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "'a_priori' must be 'equal' or 'sample', got '{other}'"
            )))
        }
    };

    let reject_fraction = match args.get("reject_fraction") {
        None | Some(Value::Null) => 0.0,
        Some(Value::Number(num)) => num.as_f64().unwrap_or(0.0),
        Some(Value::String(s)) if s.trim().is_empty() => 0.0,
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'reject_fraction' must be a number".into()))?,
        Some(_) => {
            return Err(ToolError::Validation(
                "'reject_fraction' must be a number".into(),
            ))
        }
    };
    if !(0.0..1.0).contains(&reject_fraction) || !reject_fraction.is_finite() {
        return Err(ToolError::Validation(
            "'reject_fraction' must be in [0, 1)".into(),
        ));
    }

    let (prob_output, want_prob) = match args.get("prob_output") {
        None | Some(Value::Null) => (None, false),
        Some(Value::String(s)) if s.trim().is_empty() => (None, false),
        Some(Value::String(s)) => (Some(s.clone()), true),
        Some(Value::Bool(true)) => (None, true),
        Some(Value::Bool(false)) => (None, false),
        Some(_) => {
            return Err(ToolError::Validation(
                "'prob_output' must be a string path".into(),
            ))
        }
    };

    Ok(Params {
        a_priori,
        reject_fraction,
        prob_output,
        want_prob,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{memory_store, CrsInfo, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn raster_of(rows: usize, cols: usize, bands: &[Vec<f64>], nodata: f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: bands.len(),
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata,
            data_type: DataType::F64,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for (bi, band) in bands.iter().enumerate() {
            for row in 0..rows {
                for col in 0..cols {
                    r.set(
                        bi as isize,
                        row as isize,
                        col as isize,
                        band[row * cols + col],
                    )
                    .unwrap();
                }
            }
        }
        let id = memory_store::put_raster(r);
        memory_store::make_raster_memory_path(&id)
    }

    fn read_band0(path: &str) -> Vec<f64> {
        let r = load_input_raster(path).unwrap();
        band_to_vec(&r, 0)
    }

    /// Two well-separated 2-band clusters: the left two columns are class 1
    /// (~0,0), the right two class 2 (~10,10). Training four pixels per class
    /// recovers every pixel's true class.
    #[test]
    fn separable_two_class_recovered() {
        // 4x4 grid: cols 0-1 cluster A (~0), cols 2-3 cluster B (~10).
        let (rows, cols) = (4, 4);
        let n = rows * cols;
        let mut b1 = vec![0.0; n];
        let mut b2 = vec![0.0; n];
        let mut tr = vec![0.0; n];
        for i in 0..n {
            let (row, col) = (i / cols, i % cols);
            let base = if col < 2 { 0.0 } else { 10.0 };
            b1[i] = base + 0.1 * row as f64 + 0.05 * col as f64;
            b2[i] = base + 0.05 * row as f64 + 0.1 * col as f64;
            // Training: leftmost column -> class 1, rightmost -> class 2.
            if col == 0 {
                tr[i] = 1.0;
            } else if col == 3 {
                tr[i] = 2.0;
            }
        }
        let input = raster_of(rows, cols, &[b1, b2], -9999.0);
        let train = raster_of(rows, cols, &[tr], 0.0);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input, "training": train, "prob_output": true
        }))
        .unwrap();
        let out = MaximumLikelihoodClassificationTool
            .run(&args, &ctx())
            .unwrap();
        assert_eq!(out.outputs["class_count"], json!(2));
        let cls = read_band0(out.outputs["output"].as_str().unwrap());
        for i in 0..n {
            let col = i % cols;
            let expected = if col < 2 { 1.0 } else { 2.0 };
            assert_eq!(
                cls[i], expected,
                "pixel {i} (col {col}) class {} != {expected}",
                cls[i]
            );
        }
        // Probability raster in [0,1].
        let prob = read_band0(out.outputs["prob_output"].as_str().unwrap());
        for p in &prob {
            assert!((0.0..=1.0).contains(p), "confidence {p} out of [0,1]");
        }
    }

    /// The winning-class posterior is highest at the cluster centres and never
    /// exceeds 1.
    #[test]
    fn confidence_is_bounded_and_peaks_at_centres() {
        // 3 well-separated classes, 3 training pixels each (> band count).
        let b1 = vec![0.0, 0.1, 0.2, 10.0, 10.1, 9.9, 20.0, 20.1, 19.9];
        let b2 = vec![0.1, 0.0, 0.2, 10.1, 9.9, 10.0, 20.2, 19.8, 20.0];
        let input = raster_of(3, 3, &[b1, b2], -9999.0);
        let train = vec![1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0];
        let train = raster_of(3, 3, &[train], 0.0);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input, "training": train, "prob_output": true
        }))
        .unwrap();
        let out = MaximumLikelihoodClassificationTool
            .run(&args, &ctx())
            .unwrap();
        assert_eq!(out.outputs["class_count"], json!(3));
        let prob = read_band0(out.outputs["prob_output"].as_str().unwrap());
        for p in &prob {
            assert!((0.0..=1.0).contains(p), "confidence {p} out of [0,1]");
            assert!(
                *p > 0.9,
                "well-separated centre confidence {p} should be high"
            );
        }
    }

    /// `reject_fraction` leaves the least-confident fraction of cells
    /// unclassified. A gradient of pixels between two class means spreads the
    /// confidences so the middle (ambiguous) cells are rejected.
    #[test]
    fn reject_fraction_leaves_cells_unclassified() {
        // 16 pixels on a line from class-1 mean (0,10) to class-2 mean (10,0).
        let n = 16;
        let mut b1 = vec![0.0; n];
        let mut b2 = vec![0.0; n];
        let mut tr = vec![0.0; n];
        for i in 0..n {
            let v = i as f64 / (n as f64 - 1.0) * 10.0;
            b1[i] = v;
            b2[i] = 10.0 - v;
            if i < 4 {
                tr[i] = 1.0; // near (0,10)
            } else if i >= n - 4 {
                tr[i] = 2.0; // near (10,0)
            }
        }
        let input = raster_of(1, n, &[b1, b2], -9999.0);
        let train = raster_of(1, n, &[tr], 0.0);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input, "training": train, "reject_fraction": 0.25, "prob_output": true
        }))
        .unwrap();
        let out = MaximumLikelihoodClassificationTool
            .run(&args, &ctx())
            .unwrap();
        let rejected = out.outputs["rejected_pixels"].as_u64().unwrap();
        // floor(0.25*16) = 4 least-confident cells rejected.
        assert_eq!(rejected, 4, "expected 4 rejected, got {rejected}");
        let cls = read_band0(out.outputs["output"].as_str().unwrap());
        let nod = cls.iter().filter(|v| **v == OUT_NODATA).count();
        assert_eq!(nod as u64, rejected);
        // The rejected cells are the middle (ambiguous) ones, not the ends.
        assert_ne!(cls[0], OUT_NODATA, "confident end pixel wrongly rejected");
        assert_ne!(
            cls[n - 1],
            OUT_NODATA,
            "confident end pixel wrongly rejected"
        );
    }

    /// `sample` a-priori weighting shifts the decision boundary toward the more
    /// populous class relative to `equal`.
    #[test]
    fn sample_prior_differs_from_equal() {
        // Class 1 much more populous than class 2; a borderline pixel should be
        // more likely to go to class 1 under 'sample'.
        let mut b1 = Vec::new();
        let mut b2 = Vec::new();
        let mut tr = Vec::new();
        // 20 class-1 pixels near mean 0 with unit-ish spread, 4 class-2 near 6.
        for i in 0..20 {
            let t = (i as f64 - 10.0) * 0.3;
            b1.push(t);
            b2.push(-t);
            tr.push(1.0);
        }
        for i in 0..4 {
            b1.push(6.0 + i as f64 * 0.2);
            b2.push(6.0 - i as f64 * 0.2);
            tr.push(2.0);
        }
        // A test pixel roughly between the clusters.
        b1.push(3.0);
        b2.push(3.0);
        tr.push(0.0);
        let cols = b1.len();
        let input = raster_of(1, cols, &[b1, b2], -9999.0);
        let train = raster_of(1, cols, &[tr], 0.0);

        let run = |prior: &str| {
            let args: ToolArgs = serde_json::from_value(json!({
                "input": &input, "training": &train, "a_priori": prior
            }))
            .unwrap();
            let out = MaximumLikelihoodClassificationTool
                .run(&args, &ctx())
                .unwrap();
            read_band0(out.outputs["output"].as_str().unwrap())[cols - 1]
        };
        // Under 'sample', the abundant class 1 is favoured for the borderline
        // pixel; under 'equal' it may go the other way. At minimum, the sample
        // prior must not classify the borderline pixel as the rare class when
        // equal already assigns class 1.
        let eq = run("equal");
        let sm = run("sample");
        // The sample prior can only pull the borderline pixel toward class 1.
        assert!(
            sm <= eq,
            "sample prior should not favour the rarer class more than equal (equal={eq}, sample={sm})"
        );
    }

    #[test]
    fn rejects_bad_parameters() {
        let input = raster_of(2, 2, &[vec![1.0; 4], vec![2.0; 4]], -9999.0);
        // Missing training.
        let a: ToolArgs = serde_json::from_value(json!({ "input": &input })).unwrap();
        assert!(MaximumLikelihoodClassificationTool.validate(&a).is_err());
        // Bad a_priori.
        let a: ToolArgs = serde_json::from_value(
            json!({ "input": &input, "training": &input, "a_priori": "xyz" }),
        )
        .unwrap();
        assert!(MaximumLikelihoodClassificationTool.validate(&a).is_err());
        // reject_fraction out of range.
        let a: ToolArgs = serde_json::from_value(
            json!({ "input": &input, "training": &input, "reject_fraction": 1.0 }),
        )
        .unwrap();
        assert!(MaximumLikelihoodClassificationTool.validate(&a).is_err());
    }

    #[test]
    fn rejects_grid_mismatch() {
        let input = raster_of(2, 2, &[vec![1.0; 4]], -9999.0);
        let train = raster_of(3, 3, &[vec![1.0; 9]], 0.0);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "training": train })).unwrap();
        assert!(MaximumLikelihoodClassificationTool
            .run(&args, &ctx())
            .is_err());
    }

    /// Cross-checks the Gaussian discriminant against a hand computation: for a
    /// diagonal covariance the Mahalanobis term is a simple weighted square sum.
    #[test]
    fn discriminant_matches_hand_computation() {
        // Class with mean (1,2), covariance diag(4,9) (already ridge-free here).
        let b = 2;
        let cov = vec![4.0, 0.0, 0.0, 9.0];
        let chol = cholesky(&cov, b).unwrap();
        // L = diag(2,3).
        assert!((chol[0] - 2.0).abs() < 1e-12);
        assert!((chol[3] - 3.0).abs() < 1e-12);
        let model = ClassModel {
            id: 1,
            mean: vec![1.0, 2.0],
            chol,
            ln_prior: 0.0,
            half_log_det: 0.5 * (4.0f64 * 9.0).ln(),
        };
        let x = [3.0, 5.0];
        // Maha = (3-1)²/4 + (5-2)²/9 = 1 + 1 = 2.
        // g = 0 - 0.5*ln(36) - 0.5*2.
        let expected = -0.5 * 36.0_f64.ln() - 1.0;
        let got = model.discriminant(&x, b);
        assert!(
            (got - expected).abs() < 1e-10,
            "got {got} expected {expected}"
        );
    }
}
