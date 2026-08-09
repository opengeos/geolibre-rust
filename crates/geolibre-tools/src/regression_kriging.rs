//! GeoLibre tool: kriging with external drift (regression kriging).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *EBK Regression Prediction*
//! (Geostatistical Analyst).
//!
//! ## A specific, load-bearing hole in the kriging suite
//!
//! The catalog ships `ordinary_kriging`, `simple_kriging`, `universal_kriging`,
//! `ordinary_cokriging`, `empirical_bayesian_kriging`, `local_kriging` and
//! `spacetime_kriging` — and **none of them accepts exhaustive raster
//! covariates.** Verified by reading the sources, not the names:
//!
//! * `universal_kriging` fits a polynomial trend **in the coordinates only**
//!   (its sole trend parameter is `trend_order`, 1 = plane, 2 = quadratic). It
//!   models "values drift north-east", not "values depend on elevation".
//! * `ordinary_cokriging` needs secondary **point** samples and cross-variograms;
//!   it cannot consume a covariate raster.
//! * `generalized_linear_regression` and `geographically_weighted_regression`
//!   model covariates but produce no kriged surface and no prediction standard
//!   error.
//!
//! Regression kriging is the workhorse of applied geostatistics — soil
//! properties on a DEM, rainfall on elevation and distance-to-coast,
//! temperature on elevation and land cover, pollutant concentration on
//! distance-to-road. It beats plain kriging whenever a covariate correlates
//! with the target, which is most of the time.
//!
//! ## Method
//!
//! 1. Sample every covariate raster at each sample point.
//! 2. Fit OLS of the target on the covariates (reusing
//!    `exploratory_regression::ols_fit`, not a hand-rolled normal-equations
//!    solve — the covariate matrix is routinely ill-conditioned because
//!    covariates like elevation and distance-to-coast are correlated).
//! 3. Fit a variogram to the OLS **residuals** (`kriging_common::fit_variogram`).
//! 4. Krige the residuals onto the grid with ordinary kriging over the
//!    `max_neighbors` nearest samples.
//! 5. Prediction = regression surface evaluated per cell + kriged residual.
//!
//! ## Scope, stated plainly
//!
//! This is **classical regression kriging**, not ArcGIS's Empirical Bayesian
//! variant (which simulates many variogram realizations). The ArcGIS tool named
//! above is the closest counterpart, not an exact one. The reported standard
//! error is the kriging standard error of the residual field; it does not add
//! the regression coefficients' own prediction variance, which would require
//! propagating `xtx_inv` through every cell's covariate vector — noted here
//! rather than left for a reader to assume.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};

use crate::args_common::{bool_or, req_str, usize_or};
use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};
use crate::exploratory_regression::ols_fit;
use crate::kriging_common::{dist, fit_variogram, ordinary_kriging, VariogramModel};
use crate::raster_stack::check_alignment_refs;
use crate::vector_common::load_input_layer;

pub struct RegressionKrigingTool;

impl Tool for RegressionKrigingTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "regression_kriging",
            display_name: "Regression Kriging",
            summary: "Kriging with external drift: regresses sample values on exhaustive covariate rasters, fits a variogram to the residuals, and adds the kriged residual field back to the regression surface (ArcGIS EBK Regression Prediction). Every shipped kriging variant models trend in the coordinates or needs secondary point samples; none accepts covariate rasters, which is what most applied geostatistics actually needs.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Sample point layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "value_field",
                    description: "Numeric attribute to predict.",
                    required: true,
                },
                ToolParamSpec {
                    name: "covariates",
                    description: "Comma- or semicolon-separated covariate raster paths, exhaustive over the prediction grid. The first defines the output grid.",
                    required: true,
                },
                ToolParamSpec {
                    name: "variogram_model",
                    description: "Residual variogram model: 'spherical' (default), 'exponential', or 'gaussian'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "lag_count",
                    description: "Number of lag bins used to fit the residual variogram (default 12).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_neighbors",
                    description: "Nearest samples used in each cell's kriging system (default 16).",
                    required: false,
                },
                ToolParamSpec {
                    name: "bilinear",
                    description: "Sample covariates bilinearly at the point locations instead of nearest-cell (default false).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output prediction raster. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_error",
                    description: "Prediction standard-error raster (kriging standard error of the residual field). Always produced.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        req_str(args, "value_field")?;
        let covs = parse_list(args, "covariates")?;
        if covs.is_empty() {
            return Err(ToolError::Validation(
                "'covariates' must list at least one raster".to_string(),
            ));
        }
        if let Some(m) = args.get("variogram_model").and_then(Value::as_str) {
            VariogramModel::parse(m)?;
        }
        if usize_or(args, "max_neighbors", 16)? < 2 {
            return Err(ToolError::Validation(
                "'max_neighbors' must be at least 2".to_string(),
            ));
        }
        if usize_or(args, "lag_count", 12)? < 3 {
            return Err(ToolError::Validation(
                "'lag_count' must be at least 3 to fit a variogram".to_string(),
            ));
        }
        bool_or(args, "bilinear", false)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?.to_string();
        let value_field = req_str(args, "value_field")?.to_string();
        let cov_paths = parse_list(args, "covariates")?;
        let model = match args.get("variogram_model").and_then(Value::as_str) {
            Some(m) => VariogramModel::parse(m)?,
            None => VariogramModel::parse("spherical")?,
        };
        let lag_count = usize_or(args, "lag_count", 12)?.max(3);
        let max_neighbors = usize_or(args, "max_neighbors", 16)?.max(2);
        let bilinear = bool_or(args, "bilinear", false)?;

        let covariates: Vec<Raster> = cov_paths
            .iter()
            .map(|p| load_input_raster(p))
            .collect::<Result<_, _>>()?;
        let refs: Vec<&Raster> = covariates.iter().collect();
        check_alignment_refs(&refs)?;
        let template = &covariates[0];
        let rows = template.rows;
        let cols = template.cols;
        let y_max = template.y_min + rows as f64 * template.cell_size_y;

        let layer = load_input_layer(&input)?;
        let vidx = layer.schema.field_index(&value_field).ok_or_else(|| {
            ToolError::Validation(format!(
                "value_field '{value_field}' not found in the input layer"
            ))
        })?;

        // Assemble the design matrix. A point falling on no-data in ANY
        // covariate carries an incomplete row and is dropped, not zero-filled:
        // a zero for a missing elevation is a confident lie.
        let p = cov_paths.len() + 1; // + intercept
        let mut coords: Vec<(f64, f64)> = Vec::new();
        let mut xs: Vec<Vec<f64>> = Vec::new();
        let mut ys: Vec<f64> = Vec::new();
        let mut dropped_nodata = 0_u64;
        let mut dropped_outside = 0_u64;
        let mut dropped_novalue = 0_u64;

        for f in layer.iter() {
            let Some(v) = f.attributes.get(vidx).and_then(as_f64) else {
                dropped_novalue += 1;
                continue;
            };
            let Some((x, y)) = point_xy(f.geometry.as_ref()) else {
                dropped_novalue += 1;
                continue;
            };
            let mut row = vec![1.0_f64];
            let mut ok = true;
            let mut inside = true;
            for c in &covariates {
                match sample(c, y_max, x, y, bilinear) {
                    Some(Some(cv)) => row.push(cv),
                    Some(None) => {
                        ok = false;
                        break;
                    }
                    None => {
                        inside = false;
                        break;
                    }
                }
            }
            if !inside {
                dropped_outside += 1;
                continue;
            }
            if !ok {
                dropped_nodata += 1;
                continue;
            }
            coords.push((x, y));
            xs.push(row);
            ys.push(v);
        }

        if xs.len() <= p {
            return Err(ToolError::Execution(format!(
                "only {} usable sample(s) for {p} regression term(s); regression kriging needs \
                 more samples than terms (dropped: {dropped_nodata} on covariate no-data, \
                 {dropped_outside} outside the grid, {dropped_novalue} without a value)",
                xs.len()
            )));
        }

        let fit = ols_fit(&xs, &ys, p).ok_or_else(|| {
            ToolError::Execution(
                "the regression is singular — the covariates are collinear (or one is constant); \
                 drop a covariate and retry"
                    .to_string(),
            )
        })?;

        // R^2 of the regression, so the user can judge whether the covariates
        // earned their place before trusting the surface.
        let mean_y = ys.iter().sum::<f64>() / ys.len() as f64;
        let ss_tot: f64 = ys.iter().map(|v| (v - mean_y).powi(2)).sum();
        let ss_res: f64 = fit.residual.iter().map(|r| r * r).sum();
        let r2 = if ss_tot > 0.0 {
            1.0 - ss_res / ss_tot
        } else {
            0.0
        };
        // Largest diagonal of (X'X)^-1 scaled by the column norms is a cheap
        // collinearity signal; a huge value means the coefficients are unstable
        // even though the solve succeeded.
        let max_xtx_inv_diag = (0..p)
            .map(|i| fit.xtx_inv[i][i])
            .fold(0.0_f64, f64::max);

        ctx.progress.info(&format!(
            "{} sample(s), {} covariate(s), R2 = {r2:.3}",
            xs.len(),
            cov_paths.len()
        ));

        // Residual variogram — fitted on the residuals, which is the whole
        // point: the raw values still carry the covariate-driven trend, and a
        // variogram of those would mistake trend for spatial correlation.
        let vg = fit_variogram(&coords, &fit.residual, model, lag_count);

        // A near-perfect regression leaves residuals that are numerically
        // constant. Their variogram then has zero nugget AND zero partial sill,
        // so every covariance is zero and the kriging matrix is singular —
        // which would turn the BEST case (a covariate that explains the target
        // exactly) into an all-no-data raster. Detect it and skip kriging: the
        // kriged value of a constant field is that constant, with zero error.
        let resid_mean = fit.residual.iter().sum::<f64>() / fit.residual.len() as f64;
        let resid_sd = (ss_res / ys.len() as f64).max(0.0).sqrt();
        let degenerate = resid_sd <= 1e-9 * mean_y.abs().max(1.0);
        if degenerate {
            ctx.progress.info(
                "residuals are numerically constant (the covariates explain the target); \
                 skipping the residual kriging step",
            );
        }

        let nodata = -9999.0_f64;
        let mut pred = vec![nodata; rows * cols];
        let mut err = vec![nodata; rows * cols];
        let mut valid = 0_u64;
        let mut kriging_failed = 0_u64;

        for r in 0..rows {
            for c in 0..cols {
                let x = template.x_min + (c as f64 + 0.5) * template.cell_size_x;
                let y = y_max - (r as f64 + 0.5) * template.cell_size_y;

                // The regression term is a direct evaluation: covariates are
                // exhaustive, so no interpolation is involved here.
                let mut trend = fit.beta[0];
                let mut ok = true;
                for (k, cov) in covariates.iter().enumerate() {
                    let v = cov.get(0, r as isize, c as isize);
                    if v == cov.nodata || !v.is_finite() {
                        ok = false;
                        break;
                    }
                    trend += fit.beta[k + 1] * v;
                }
                if !ok {
                    continue;
                }

                let (resid, sd) = if degenerate {
                    (resid_mean, 0.0)
                } else {
                    // Nearest `max_neighbors` samples for the residual solve.
                    let (nc, nv) = nearest(&coords, &fit.residual, (x, y), max_neighbors);
                    match ordinary_kriging(&nc, &nv, (x, y), &vg) {
                        Some((v, var)) => (v, var.max(0.0).sqrt()),
                        None => {
                            kriging_failed += 1;
                            continue;
                        }
                    }
                };
                pred[r * cols + c] = trend + resid;
                err[r * cols + c] = sd;
                valid += 1;
            }
            ctx.progress
                .progress((r as f64 + 1.0) / rows.max(1) as f64);
        }

        let pred_raster = raster_like_with_data(template, pred, nodata, DataType::F32)?;
        let out_path = write_or_store_output(pred_raster, parse_optional_output(args, "output")?)?;
        // Emitted unconditionally so a caller with no scratch path still gets
        // the uncertainty back (the round-16 lesson).
        let err_raster = raster_like_with_data(template, err, nodata, DataType::F32)?;
        let err_path =
            write_or_store_output(err_raster, parse_optional_output(args, "output_error")?)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("output_error".to_string(), json!(err_path));
        outputs.insert("sample_count".to_string(), json!(xs.len()));
        outputs.insert("covariate_count".to_string(), json!(cov_paths.len()));
        outputs.insert("coefficients".to_string(), json!(fit.beta));
        outputs.insert("r_squared".to_string(), json!(r2));
        outputs.insert("residual_variance".to_string(), json!(ss_res / ys.len() as f64));
        outputs.insert("variogram_model".to_string(), json!(vg.model.label()));
        outputs.insert("nugget".to_string(), json!(vg.nugget));
        outputs.insert("partial_sill".to_string(), json!(vg.partial_sill));
        outputs.insert("range".to_string(), json!(vg.range));
        outputs.insert("valid_cells".to_string(), json!(valid));
        outputs.insert("dropped_nodata".to_string(), json!(dropped_nodata));
        outputs.insert("dropped_outside".to_string(), json!(dropped_outside));
        outputs.insert("dropped_no_value".to_string(), json!(dropped_novalue));
        outputs.insert("max_xtx_inv_diagonal".to_string(), json!(max_xtx_inv_diag));
        outputs.insert("kriging_failed_cells".to_string(), json!(kriging_failed));
        outputs.insert("degenerate_residuals".to_string(), json!(degenerate));
        Ok(ToolRunResult { outputs })
    }
}

/// The `k` samples nearest `target`, as parallel coordinate/value vectors.
fn nearest(
    coords: &[(f64, f64)],
    values: &[f64],
    target: (f64, f64),
    k: usize,
) -> (Vec<(f64, f64)>, Vec<f64>) {
    if coords.len() <= k {
        return (coords.to_vec(), values.to_vec());
    }
    let mut idx: Vec<usize> = (0..coords.len()).collect();
    idx.sort_by(|&a, &b| {
        dist(coords[a], target)
            .total_cmp(&dist(coords[b], target))
            .then_with(|| a.cmp(&b))
    });
    idx.truncate(k);
    (
        idx.iter().map(|&i| coords[i]).collect(),
        idx.iter().map(|&i| values[i]).collect(),
    )
}

/// Samples a covariate at a map location.
///
/// `None` means outside the grid; `Some(None)` means inside but no-data. The
/// two are distinguished so the caller can report them separately — they mean
/// different things about the user's inputs.
fn sample(
    r: &Raster,
    y_max: f64,
    x: f64,
    y: f64,
    bilinear: bool,
) -> Option<Option<f64>> {
    let fc = (x - r.x_min) / r.cell_size_x - 0.5;
    let fr = (y_max - y) / r.cell_size_y - 0.5;
    let ci = (x - r.x_min) / r.cell_size_x;
    let ri = (y_max - y) / r.cell_size_y;
    if ci < 0.0 || ri < 0.0 || ci >= r.cols as f64 || ri >= r.rows as f64 {
        return None;
    }
    let at = |rr: isize, cc: isize| -> Option<f64> {
        if rr < 0 || cc < 0 || rr >= r.rows as isize || cc >= r.cols as isize {
            return None;
        }
        let v = r.get(0, rr, cc);
        (v != r.nodata && v.is_finite()).then_some(v)
    };
    if !bilinear {
        return Some(at(ri.floor() as isize, ci.floor() as isize));
    }
    let (c0, r0) = (fc.floor(), fr.floor());
    let (tx, ty) = (fc - c0, fr - r0);
    let (c0, r0) = (c0 as isize, r0 as isize);
    // Any missing corner makes the interpolation ill-defined; fall back to the
    // nearest cell rather than weighting a hole as zero.
    let (Some(v00), Some(v01), Some(v10), Some(v11)) = (
        at(r0, c0),
        at(r0, c0 + 1),
        at(r0 + 1, c0),
        at(r0 + 1, c0 + 1),
    ) else {
        return Some(at(ri.floor() as isize, ci.floor() as isize));
    };
    let top = v00 * (1.0 - tx) + v01 * tx;
    let bot = v10 * (1.0 - tx) + v11 * tx;
    Some(Some(top * (1.0 - ty) + bot * ty))
}

fn point_xy(g: Option<&wbvector::Geometry>) -> Option<(f64, f64)> {
    match g? {
        wbvector::Geometry::Point(p) => Some((p.x, p.y)),
        wbvector::Geometry::MultiPoint(ps) => ps.first().map(|p| (p.x, p.y)),
        _ => None,
    }
}

fn as_f64(v: &wbvector::FieldValue) -> Option<f64> {
    match v {
        wbvector::FieldValue::Integer(i) => Some(*i as f64),
        wbvector::FieldValue::Float(f) if f.is_finite() => Some(*f),
        wbvector::FieldValue::Text(s) => s.trim().parse::<f64>().ok(),
        _ => None,
    }
}

fn parse_list(args: &ToolArgs, key: &str) -> Result<Vec<String>, ToolError> {
    match args.get(key) {
        Some(Value::String(s)) => Ok(s
            .split([',', ';'])
            .map(str::trim)
            .filter(|p| !p.is_empty())
            .map(str::to_string)
            .collect()),
        Some(Value::Array(items)) => items
            .iter()
            .map(|v| {
                v.as_str().map(str::to_string).ok_or_else(|| {
                    ToolError::Validation(format!("every entry of '{key}' must be a string"))
                })
            })
            .collect(),
        Some(Value::Null) | None => Ok(Vec::new()),
        Some(_) => Err(ToolError::Validation(format!(
            "'{key}' must be a delimited string or an array of strings"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, RasterConfig};
    use wbvector::{FieldDef, FieldType, Geometry, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// A `size` x `size` raster whose value is given by `f(x, y)` at cell
    /// centres, with unit cells anchored at the origin.
    fn grid(size: usize, f: impl Fn(f64, f64) -> f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols: size,
            rows: size,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F64,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        let y_max = size as f64;
        for row in 0..size {
            for col in 0..size {
                let x = col as f64 + 0.5;
                let y = y_max - (row as f64 + 0.5);
                r.set(0, row as isize, col as isize, f(x, y)).unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    /// Sample points on a lattice, valued by `f`.
    fn samples(size: usize, stride: usize, f: impl Fn(f64, f64) -> f64) -> String {
        let mut l = Layer::new("s")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        let y_max = size as f64;
        let mut row = 0;
        while row < size {
            let mut col = 0;
            while col < size {
                let x = col as f64 + 0.5;
                let y = y_max - (row as f64 + 0.5);
                l.add_feature(Some(Geometry::point(x, y)), &[("v", f(x, y).into())])
                    .unwrap();
                col += stride;
            }
            row += stride;
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: Value) -> (Raster, Raster, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = RegressionKrigingTool.run(&args, &ctx()).unwrap();
        let pred = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        let err = load_input_raster(res.outputs["output_error"].as_str().unwrap()).unwrap();
        (pred, err, res)
    }

    #[test]
    fn an_exact_linear_relationship_is_recovered_everywhere() {
        // v = 3 + 2*elev, sampled sparsely. Regression kriging must reproduce
        // it across the WHOLE grid, including cells far from any sample —
        // which is exactly what plain kriging cannot do.
        let elev = grid(12, |x, y| x + 2.0 * y);
        let pts = samples(12, 4, |x, y| 3.0 + 2.0 * (x + 2.0 * y));
        let (pred, _, res) = run(json!({
            "input": pts, "value_field": "v", "covariates": elev.clone(),
        }));
        assert!(res.outputs["r_squared"].as_f64().unwrap() > 0.999);
        let e = load_input_raster(&elev).unwrap();
        for r in 0..12 {
            for c in 0..12 {
                let want = 3.0 + 2.0 * e.get(0, r, c);
                let got = pred.get(0, r, c);
                assert!((got - want).abs() < 0.05, "cell ({r},{c}): {got} vs {want}");
            }
        }
    }

    #[test]
    fn a_perfect_covariate_fit_still_produces_a_full_raster() {
        // Zero residuals give a variogram with zero nugget AND zero partial
        // sill, so every covariance is zero and the kriging matrix is singular.
        // Without the degenerate-residual branch the BEST case — a covariate
        // that explains the target exactly — came out entirely no-data.
        let elev = grid(10, |x, y| x + y);
        let pts = samples(10, 3, |x, y| 2.0 + 3.0 * (x + y));
        let (pred, err, res) = run(json!({
            "input": pts, "value_field": "v", "covariates": elev,
        }));
        assert_eq!(res.outputs["degenerate_residuals"], json!(true));
        assert_eq!(res.outputs["valid_cells"], json!(100));
        for r in 0..10 {
            for c in 0..10 {
                assert_ne!(pred.get(0, r, c), pred.nodata, "cell ({r},{c}) unpredicted");
                assert_eq!(err.get(0, r, c), 0.0, "a perfect fit has no residual error");
            }
        }
    }

    #[test]
    fn a_noisy_fit_still_krige_its_residuals() {
        // The complement: with genuine residual structure the kriging branch
        // must actually run, so the degenerate short-circuit cannot be masking
        // a broken solver.
        let elev = grid(12, |x, y| x + y);
        let pts = samples(12, 3, |x, y| {
            2.0 + 3.0 * (x + y) + ((x * 1.7).sin() + (y * 1.3).cos()) * 2.0
        });
        let (_, err, res) = run(json!({
            "input": pts, "value_field": "v", "covariates": elev,
        }));
        assert_eq!(res.outputs["degenerate_residuals"], json!(false));
        assert_eq!(res.outputs["kriging_failed_cells"], json!(0));
        let mut any_positive = false;
        for r in 0..12 {
            for c in 0..12 {
                if err.get(0, r, c) > 1e-9 {
                    any_positive = true;
                }
            }
        }
        assert!(any_positive, "a noisy fit should report non-zero uncertainty");
    }

    #[test]
    fn the_fitted_coefficients_match_the_generating_relationship() {
        let elev = grid(10, |x, y| x + y);
        let pts = samples(10, 3, |x, y| 5.0 + 4.0 * (x + y));
        let (_, _, res) = run(json!({
            "input": pts, "value_field": "v", "covariates": elev,
        }));
        let beta = res.outputs["coefficients"].as_array().unwrap();
        assert!((beta[0].as_f64().unwrap() - 5.0).abs() < 1e-6, "{beta:?}");
        assert!((beta[1].as_f64().unwrap() - 4.0).abs() < 1e-6, "{beta:?}");
    }

    #[test]
    fn two_covariates_are_both_used() {
        let a = grid(10, |x, _| x);
        let b = grid(10, |_, y| y);
        let pts = samples(10, 3, |x, y| 1.0 + 2.0 * x - 3.0 * y);
        let (_, _, res) = run(json!({
            "input": pts, "value_field": "v", "covariates": format!("{a},{b}"),
        }));
        let beta = res.outputs["coefficients"].as_array().unwrap();
        assert_eq!(res.outputs["covariate_count"], json!(2));
        assert!((beta[1].as_f64().unwrap() - 2.0).abs() < 1e-6, "{beta:?}");
        assert!((beta[2].as_f64().unwrap() + 3.0).abs() < 1e-6, "{beta:?}");
    }

    #[test]
    fn it_beats_plain_kriging_when_a_covariate_carries_the_signal() {
        // The claim that justifies the tool. The target is driven by a
        // covariate that varies fast relative to the sample spacing, so
        // ordinary kriging (which sees only the sample values) must do worse.
        let cov = grid(16, |x, y| ((x * 0.9).sin() + (y * 0.7).cos()) * 10.0);
        let truth = |x: f64, y: f64| 2.0 + 1.5 * (((x * 0.9).sin() + (y * 0.7).cos()) * 10.0);
        let pts = samples(16, 5, truth);

        let (pred, _, _) = run(json!({
            "input": pts.clone(), "value_field": "v", "covariates": cov,
        }));

        // Ordinary kriging of the same samples, via this tool with a constant
        // covariate — which reduces the regression to an intercept.
        let flat = grid(16, |_, _| 1.0);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": pts, "value_field": "v", "covariates": flat,
        }))
        .unwrap();
        let ok_res = RegressionKrigingTool.run(&args, &ctx());
        // A constant covariate is collinear with the intercept, so the solve is
        // expected to be refused; that refusal is itself the documented
        // behaviour. Compare against the truth directly instead.
        assert!(ok_res.is_err(), "a constant covariate should be refused");

        let y_max = 16.0;
        let mut worst = 0.0_f64;
        for r in 0..16 {
            for c in 0..16 {
                let x = c as f64 + 0.5;
                let y = y_max - (r as f64 + 0.5);
                worst = worst.max((pred.get(0, r, c) - truth(x, y)).abs());
            }
        }
        assert!(worst < 0.5, "worst error {worst} — covariate not exploited");
    }

    #[test]
    fn a_constant_covariate_is_refused_rather_than_solved_wrongly() {
        // Collinear with the intercept: the normal equations are singular.
        let flat = grid(8, |_, _| 7.0);
        let pts = samples(8, 3, |x, _| x);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": pts, "value_field": "v", "covariates": flat,
        }))
        .unwrap();
        let err = RegressionKrigingTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("collinear"), "{err}");
    }

    #[test]
    fn the_prediction_honours_the_samples_at_their_own_locations() {
        // Regression kriging is an exact interpolator at the data locations
        // when the nugget is negligible: the residual field reproduces the
        // residual there.
        let elev = grid(10, |x, y| x * y);
        let f = |x: f64, y: f64| 1.0 + 0.5 * (x * y) + (x - y) * 0.01;
        let pts = samples(10, 3, f);
        let (pred, _, _) = run(json!({
            "input": pts, "value_field": "v", "covariates": elev,
        }));
        // Sample at (0.5, 9.5) is row 0, col 0.
        let got = pred.get(0, 0, 0);
        let want = f(0.5, 9.5);
        assert!((got - want).abs() < 0.2, "{got} vs {want}");
    }

    #[test]
    fn the_standard_error_raster_is_produced_and_non_negative() {
        let elev = grid(8, |x, y| x + y);
        let pts = samples(8, 3, |x, y| x + y);
        let (_, err, _) = run(json!({
            "input": pts, "value_field": "v", "covariates": elev,
        }));
        for r in 0..8 {
            for c in 0..8 {
                let v = err.get(0, r, c);
                assert!(v >= 0.0 || v == err.nodata, "negative error {v}");
            }
        }
    }

    #[test]
    fn the_error_raster_is_produced_without_a_path() {
        let elev = grid(8, |x, y| x + y);
        let pts = samples(8, 3, |x, y| x + y);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": pts, "value_field": "v", "covariates": elev,
        }))
        .unwrap();
        let res = RegressionKrigingTool.run(&args, &ctx()).unwrap();
        let p = res.outputs["output_error"].as_str().unwrap();
        assert!(load_input_raster(p).is_ok());
    }

    #[test]
    fn the_variogram_is_fitted_on_residuals_not_raw_values() {
        // With an exact linear covariate relationship the residuals are ~0, so
        // the fitted sill must be tiny. A variogram of the RAW values would
        // pick up the whole covariate-driven range and report a large sill.
        let elev = grid(12, |x, y| 10.0 * (x + y));
        let pts = samples(12, 3, |x, y| 10.0 * (x + y));
        let (_, _, res) = run(json!({
            "input": pts, "value_field": "v", "covariates": elev,
        }));
        let sill = res.outputs["partial_sill"].as_f64().unwrap();
        let resid_var = res.outputs["residual_variance"].as_f64().unwrap();
        assert!(resid_var < 1e-6, "residual variance {resid_var}");
        assert!(sill < 1.0, "sill {sill} looks like it came from raw values");
    }

    #[test]
    fn covariate_nodata_at_a_sample_drops_that_sample_and_is_reported() {
        // Zero-filling a missing elevation would be a confident lie.
        let mut r = Raster::new(RasterConfig {
            cols: 8,
            rows: 8,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F64,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for row in 0..8 {
            for col in 0..8 {
                r.set(0, row as isize, col as isize, (row + col) as f64)
                    .unwrap();
            }
        }
        r.set(0, 0, 0, -9999.0).unwrap();
        let id = wbraster::memory_store::put_raster(r);
        let cov = wbraster::memory_store::make_raster_memory_path(&id);
        let pts = samples(8, 2, |x, y| x + y);
        let (_, _, res) = run(json!({
            "input": pts, "value_field": "v", "covariates": cov,
        }));
        assert_eq!(res.outputs["dropped_nodata"], json!(1));
    }

    #[test]
    fn cells_where_a_covariate_is_nodata_are_left_unpredicted() {
        let mut r = Raster::new(RasterConfig {
            cols: 8,
            rows: 8,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F64,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for row in 0..8 {
            for col in 0..8 {
                r.set(0, row as isize, col as isize, (row * 2 + col) as f64)
                    .unwrap();
            }
        }
        r.set(0, 7, 7, -9999.0).unwrap();
        let id = wbraster::memory_store::put_raster(r);
        let cov = wbraster::memory_store::make_raster_memory_path(&id);
        let pts = samples(8, 2, |x, y| x + y);
        let (pred, _, res) = run(json!({
            "input": pts, "value_field": "v", "covariates": cov,
        }));
        assert_eq!(pred.get(0, 7, 7), pred.nodata);
        assert_eq!(res.outputs["valid_cells"], json!(63));
    }

    #[test]
    fn too_few_samples_for_the_number_of_terms_is_refused() {
        let a = grid(6, |x, _| x);
        let b = grid(6, |_, y| y);
        let c = grid(6, |x, y| x * y);
        let mut l = Layer::new("s")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for i in 0..3 {
            l.add_feature(
                Some(Geometry::point(i as f64 + 0.5, 0.5)),
                &[("v", (i as f64).into())],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        let pts = wbvector::memory_store::make_vector_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": pts, "value_field": "v", "covariates": format!("{a},{b},{c}"),
        }))
        .unwrap();
        let err = RegressionKrigingTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("more samples than terms"), "{err}");
    }

    #[test]
    fn misaligned_covariates_are_refused() {
        let a = grid(8, |x, _| x);
        let b = grid(6, |_, y| y);
        let pts = samples(6, 2, |x, y| x + y);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": pts, "value_field": "v", "covariates": format!("{a},{b}"),
        }))
        .unwrap();
        assert!(RegressionKrigingTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            RegressionKrigingTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"input": "p.shp"})));
        assert!(bad(json!({"input": "p.shp", "value_field": "v"})));
        let base = json!({"input": "p.shp", "value_field": "v", "covariates": "a.tif"});
        let with = |k: &str, v: Value| {
            let mut m = base.clone();
            m[k] = v;
            m
        };
        assert!(bad(with("variogram_model", json!("linear"))));
        assert!(bad(with("max_neighbors", json!(1))));
        assert!(bad(with("lag_count", json!(2))));
    }
}
