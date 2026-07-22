//! GeoLibre tool: `optimal_interpolation` — data-assimilation correction of a
//! background raster field using point observations (ArcGIS Pro's *Optimal
//! Interpolation*, Image Analyst).
//!
//! Every interpolator already in the bundled suite — IDW, the kriging family,
//! thin-plate splines, natural-neighbour, EBK — builds a surface *from points
//! alone*. Optimal interpolation (OI) is different in kind: it starts from a
//! prior gridded field (the **background**, e.g. a model forecast, a coarse
//! climatology, or a previous survey) and *nudges it toward observations*,
//! weighting each nudge by how much the background is trusted (its error
//! variance `B`), how much the observations are trusted (their error variance
//! `R`), and how quickly spatial correlations decay (the `correlation_length`
//! `L`). This is the classic sequential data-assimilation analysis step
//!
//! ```text
//! x_a = x_b + B Hᵀ (H B Hᵀ + R)⁻¹ (y − H x_b)
//! ```
//!
//! specialised to a scalar field with a Gaussian spatial correlation model.
//!
//! ## Method (localized OI)
//! For every output cell we solve a small local system rather than one global
//! one (which would be an N×N solve over the whole grid):
//!   1. `H x_b` — bilinearly sample the background at each observation location.
//!   2. Innovation `d_j = y_j − (H x_b)_j`.
//!   3. Gather the observations within `cutoff · L` of the cell (vendored
//!      `kdtree` radius query; `cutoff = 3`, beyond which the Gaussian
//!      correlation is < 0.02), capped at `max_obs` nearest.
//!   4. Background-error covariance among those observations
//!      `P = σ_b² · exp(−½ (r/L)²)` plus observation error `R = diag(σ_o²)`;
//!      solve `(P + R) w = d` by Cholesky (the matrix is SPD).
//!   5. Cross-covariance cell↔obs `c_j = σ_b² · exp(−½ (dist/L)²)`; the analysis
//!      increment is `c · w`, so `x_a = x_b + c · w`.
//!   6. Analysis-error variance `σ_a² = σ_b² − c·(P+R)⁻¹c` (a second solve reuses
//!      the same factorization), written to the optional `analysis_error`
//!      raster. It is ≤ `σ_b²` everywhere (assimilation never increases
//!      uncertainty) and equals `σ_b²` where no observation is in range.
//!
//! Cells with no observation within `cutoff · L` are left exactly equal to the
//! background (pass-through) with analysis error `σ_b²`. Background no-data cells
//! stay no-data.
//!
//! ## Deliberate v1 scope cuts (documented for reviewers)
//! - Isotropic Gaussian correlation only (one scalar `L`); no anisotropy, no
//!   flow-dependent or terrain-following covariances.
//! - Localized (per-cell neighbourhood) solve, not a single global analysis —
//!   standard OI practice for large grids and what keeps this O(cells · k³) with
//!   a small `k` instead of O(N³).
//! - Scalar univariate field: one observation variable, `H` is interpolation of
//!   the background (no multivariate/balance operators).

use std::collections::BTreeMap;

use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};
use wbvector::Geometry;

use crate::common::{load_input_raster, parse_optional_output, raster_like_with_data};
use crate::vector_common::load_input_layer;

/// Correlation cutoff in units of `L`: beyond this the Gaussian weight is < 0.02
/// and is treated as zero (so the radius query and the local solve stay small).
const CUTOFF: f64 = 3.0;

pub struct OptimalInterpolationTool;

impl Tool for OptimalInterpolationTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "optimal_interpolation",
            display_name: "Optimal Interpolation",
            summary: "Correct a background raster field with point observations via the optimal-interpolation / data-assimilation analysis step x_a = x_b + B Hᵀ (H B Hᵀ + R)⁻¹ (y − H x_b), using a Gaussian spatial correlation of the given correlation length and per-observation/background error variances (like ArcGIS's Optimal Interpolation). Unlike IDW/kriging/EBK, which interpolate from points alone, this blends a prior gridded field with observations and can emit an analysis-error-variance raster.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "background",
                    description: "Background (prior) raster field to be corrected.",
                    required: true,
                },
                ToolParamSpec {
                    name: "observations",
                    description: "Observation point layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "field",
                    description: "Numeric field on the observation layer holding the observed value.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output corrected (analysis) raster. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "correlation_length",
                    description: "Spatial correlation length L in CRS units (Gaussian). Default: 5% of the background's longer side.",
                    required: false,
                },
                ToolParamSpec {
                    name: "background_error_variance",
                    description: "Background-error variance σ_b² (default 1.0). Larger => trust observations more.",
                    required: false,
                },
                ToolParamSpec {
                    name: "obs_error_variance",
                    description: "Observation-error variance σ_o² (default 1.0), used when no per-point error field is given. Larger => trust the background more.",
                    required: false,
                },
                ToolParamSpec {
                    name: "error_field",
                    description: "Optional field on the observation layer giving each point's error variance (overrides obs_error_variance where present and > 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based background band to correct (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_obs",
                    description: "Maximum nearest observations used in each cell's local solve (default 50, range 1..=200).",
                    required: false,
                },
                ToolParamSpec {
                    name: "analysis_error",
                    description: "Optional output raster of the analysis-error variance σ_a² (≤ σ_b²).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "background")?;
        require_str(args, "observations")?;
        require_str(args, "field")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let background = require_str(args, "background")?;
        let observations = require_str(args, "observations")?;
        let field = require_str(args, "field")?;
        let output = parse_optional_output(args, "output")?;
        let error_output = parse_optional_output(args, "analysis_error")?;
        let error_field = args
            .get("error_field")
            .and_then(Value::as_str)
            .map(str::trim)
            .filter(|s| !s.is_empty());
        let prm = parse_params(args)?;

        // ── Background raster ────────────────────────────────────────────────
        let bg = load_input_raster(background)?;
        let band_1based = prm.band;
        if (band_1based - 1) as usize >= bg.bands {
            return Err(ToolError::Validation(format!(
                "band {band_1based} out of range (raster has {} band(s))",
                bg.bands
            )));
        }
        let band = (band_1based - 1) as isize;
        let rows = bg.rows;
        let cols = bg.cols;
        let nodata = bg.nodata;
        let x_min = bg.x_min;
        let y_max = bg.y_max();
        let cs_x = bg.cell_size_x;
        let cs_y = bg.cell_size_y;

        let l = prm.correlation_length.unwrap_or_else(|| {
            let w = cols as f64 * cs_x;
            let h = rows as f64 * cs_y;
            (w.max(h) * 0.05).max(cs_x.max(cs_y))
        });

        // ── Observations ─────────────────────────────────────────────────────
        let layer = load_input_layer(observations)?;
        let field_idx = layer
            .schema
            .field_index(field)
            .ok_or_else(|| ToolError::Validation(format!("field '{field}' not found")))?;
        let err_idx =
            match error_field {
                Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                    ToolError::Validation(format!("error_field '{f}' not found"))
                })?),
                None => None,
            };

        let mut obs: Vec<Obs> = Vec::new();
        for feat in layer.iter() {
            let Some((x, y)) = feat.geometry.as_ref().and_then(point_xy) else {
                continue;
            };
            let Some(v) = feat.attributes.get(field_idx).and_then(|f| f.as_f64()) else {
                continue;
            };
            if !x.is_finite() || !y.is_finite() || !v.is_finite() {
                continue;
            }
            // background sampled at the observation location (H x_b)
            let Some(hxb) = bilinear(
                &bg, band, x, y, x_min, y_max, cs_x, cs_y, cols, rows, nodata,
            ) else {
                continue; // observation falls on nodata / outside the grid
            };
            let r = match err_idx {
                Some(ei) => feat
                    .attributes
                    .get(ei)
                    .and_then(|f| f.as_f64())
                    .filter(|r| r.is_finite() && *r > 0.0)
                    .unwrap_or(prm.obs_error_variance),
                None => prm.obs_error_variance,
            };
            obs.push(Obs {
                x,
                y,
                d: v - hxb, // innovation
                r,
            });
        }
        if obs.is_empty() {
            return Err(ToolError::Execution(format!(
                "no observations with a finite '{field}' value fell on valid background cells"
            )));
        }

        // ── kdtree over observations for radius queries ──────────────────────
        let mut tree: KdTree<f64, usize, [f64; 2]> = KdTree::new(2);
        for (i, o) in obs.iter().enumerate() {
            tree.add([o.x, o.y], i)
                .map_err(|e| ToolError::Execution(format!("kdtree insert failed: {e:?}")))?;
        }
        let radius = CUTOFF * l;
        let radius_sq = radius * radius;

        ctx.progress.info(&format!(
            "assimilating {} observation(s) into a {rows}x{cols} background (L={l:.3})",
            obs.len()
        ));

        let sigma_b2 = prm.background_error_variance;
        let n = rows * cols;
        let mut analysis = vec![nodata; n];
        let mut aerr = vec![nodata; n];
        let mut corrected = 0usize;

        for r in 0..rows {
            let cy = y_max - (r as f64 + 0.5) * cs_y;
            for c in 0..cols {
                let xb = bg.get(band, r as isize, c as isize);
                if xb == nodata || !xb.is_finite() {
                    continue; // stays nodata
                }
                let cx = x_min + (c as f64 + 0.5) * cs_x;

                // nearest observations within the correlation cutoff
                let mut found = tree
                    .within(&[cx, cy], radius_sq, &squared_euclidean)
                    .unwrap_or_default();
                if found.is_empty() {
                    analysis[r * cols + c] = xb; // pass-through
                    aerr[r * cols + c] = sigma_b2;
                    continue;
                }
                // `within` returns nearest-first; cap the local solve size.
                if found.len() > prm.max_obs {
                    found.truncate(prm.max_obs);
                }
                let k = found.len();

                // cross-covariance cell<->obs and innovation vector
                let mut cvec = vec![0.0f64; k];
                let mut dvec = vec![0.0f64; k];
                for (a, (dist_sq, &oi)) in found.iter().enumerate() {
                    cvec[a] = sigma_b2 * gauss(*dist_sq, l);
                    dvec[a] = obs[oi].d;
                }

                // (P + R): background covariance among obs + obs error variance
                let mut mat = vec![0.0f64; k * k];
                for a in 0..k {
                    let (_, &oa) = found[a];
                    for b in 0..k {
                        let (_, &ob) = found[b];
                        let dsq = (obs[oa].x - obs[ob].x).powi(2) + (obs[oa].y - obs[ob].y).powi(2);
                        mat[a * k + b] = sigma_b2 * gauss(dsq, l);
                    }
                    mat[a * k + a] += obs[oa].r; // diagonal R
                }

                // Solve (P+R) w = d and (P+R) u = c with one factorization.
                let (incr, var_reduction) = match cholesky(&mat, k) {
                    Some(chol) => {
                        let w = chol.solve(&dvec);
                        let u = chol.solve(&cvec);
                        let incr: f64 = cvec.iter().zip(&w).map(|(cj, wj)| cj * wj).sum();
                        let red: f64 = cvec.iter().zip(&u).map(|(cj, uj)| cj * uj).sum();
                        (incr, red)
                    }
                    None => (0.0, 0.0), // degenerate neighbourhood: fall back to background
                };

                analysis[r * cols + c] = xb + incr;
                aerr[r * cols + c] = (sigma_b2 - var_reduction).max(0.0);
                corrected += 1;
            }
            if rows > 0 {
                ctx.progress.progress((r as f64 + 1.0) / rows as f64);
            }
        }

        // ── Outputs ──────────────────────────────────────────────────────────
        let out_raster = raster_like_with_data(&bg, analysis, nodata, DataType::F32)?;
        let out_path = write_or_store_raster(out_raster, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        if let Some(p) = error_output {
            let e_raster = raster_like_with_data(&bg, aerr, nodata, DataType::F32)?;
            outputs.insert(
                "analysis_error".to_string(),
                json!(write_or_store_raster(e_raster, Some(p))?),
            );
        }
        outputs.insert("observation_count".to_string(), json!(obs.len()));
        outputs.insert("corrected_cells".to_string(), json!(corrected));
        outputs.insert("correlation_length".to_string(), json!(l));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        Ok(ToolRunResult { outputs })
    }
}

/// One observation: location, innovation `d = y − H x_b`, error variance `r`.
#[derive(Clone, Copy)]
struct Obs {
    x: f64,
    y: f64,
    d: f64,
    r: f64,
}

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        _ => None,
    }
}

/// Gaussian spatial correlation `exp(−½ (d/L)²)` from a *squared* distance.
fn gauss(dist_sq: f64, l: f64) -> f64 {
    (-0.5 * dist_sq / (l * l)).exp()
}

/// Bilinear sample of one band at world coordinates `(x, y)`, using cell-centre
/// alignment. Returns `None` if the point is outside the grid or any of the four
/// surrounding cells is no-data (so an observation over no-data is dropped
/// rather than silently corrupting the innovation).
#[allow(clippy::too_many_arguments)]
fn bilinear(
    r: &Raster,
    band: isize,
    x: f64,
    y: f64,
    x_min: f64,
    y_max: f64,
    cs_x: f64,
    cs_y: f64,
    cols: usize,
    rows: usize,
    nodata: f64,
) -> Option<f64> {
    // Fractional cell-centre coordinates.
    let fx = (x - x_min) / cs_x - 0.5;
    let fy = (y_max - y) / cs_y - 0.5;
    if !fx.is_finite() || !fy.is_finite() {
        return None;
    }
    let c0 = fx.floor();
    let r0 = fy.floor();
    let tx = fx - c0;
    let ty = fy - r0;
    let c0 = c0 as isize;
    let r0 = r0 as isize;

    let sample = |rr: isize, cc: isize| -> Option<f64> {
        if rr < 0 || cc < 0 || rr >= rows as isize || cc >= cols as isize {
            return None;
        }
        let v = r.get(band, rr, cc);
        if v == nodata || !v.is_finite() {
            None
        } else {
            Some(v)
        }
    };

    let v00 = sample(r0, c0);
    let v01 = sample(r0, c0 + 1);
    let v10 = sample(r0 + 1, c0);
    let v11 = sample(r0 + 1, c0 + 1);

    match (v00, v01, v10, v11) {
        (Some(a), Some(b), Some(c), Some(d)) => {
            let top = a * (1.0 - tx) + b * tx;
            let bot = c * (1.0 - tx) + d * tx;
            Some(top * (1.0 - ty) + bot * ty)
        }
        // On the grid edge (or one corner is nodata) fall back to nearest valid.
        _ => v00.or(v01).or(v10).or(v11),
    }
}

/// A cached lower-triangular Cholesky factor `L` of an SPD matrix (`A = L Lᵀ`),
/// so several right-hand sides can be solved against one factorization.
struct Cholesky {
    l: Vec<f64>,
    n: usize,
}

impl Cholesky {
    /// Solves `A x = b` via the stored factor.
    fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.n;
        let l = &self.l;
        // Forward solve L y = b.
        let mut y = vec![0.0f64; n];
        for i in 0..n {
            let mut sum = b[i];
            for k in 0..i {
                sum -= l[i * n + k] * y[k];
            }
            y[i] = sum / l[i * n + i];
        }
        // Backward solve Lᵀ x = y.
        let mut x = vec![0.0f64; n];
        for ii in 0..n {
            let i = n - 1 - ii;
            let mut sum = y[i];
            for k in (i + 1)..n {
                sum -= l[k * n + i] * x[k];
            }
            x[i] = sum / l[i * n + i];
        }
        x
    }
}

/// Cholesky-factorizes the row-major SPD matrix `a` (`n×n`). `None` if it is not
/// numerically positive-definite.
fn cholesky(a: &[f64], n: usize) -> Option<Cholesky> {
    let mut l = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i * n + j];
            for k in 0..j {
                sum -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if sum <= 0.0 {
                    return None;
                }
                l[i * n + j] = sum.sqrt();
            } else {
                l[i * n + j] = sum / l[j * n + j];
            }
        }
    }
    Some(Cholesky { l, n })
}

/// Writes a raster to a file path, or stores it in memory and returns a
/// `memory://` handle when no path is given. (Thin wrapper over the shared
/// helper so `run` reads cleanly for both outputs.)
fn write_or_store_raster(raster: Raster, output: Option<&str>) -> Result<String, ToolError> {
    crate::common::write_or_store_output(raster, output)
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    correlation_length: Option<f64>,
    background_error_variance: f64,
    obs_error_variance: f64,
    band: u64,
    max_obs: usize,
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let correlation_length = match parse_f64(args, "correlation_length")? {
        None => None,
        Some(v) if v > 0.0 && v.is_finite() => Some(v),
        Some(_) => {
            return Err(ToolError::Validation(
                "'correlation_length' must be a positive number".to_string(),
            ))
        }
    };

    let background_error_variance = match parse_f64(args, "background_error_variance")? {
        None => 1.0,
        Some(v) if v > 0.0 && v.is_finite() => v,
        Some(_) => {
            return Err(ToolError::Validation(
                "'background_error_variance' must be a positive number".to_string(),
            ))
        }
    };

    let obs_error_variance = match parse_f64(args, "obs_error_variance")? {
        None => 1.0,
        Some(v) if v > 0.0 && v.is_finite() => v,
        Some(_) => {
            return Err(ToolError::Validation(
                "'obs_error_variance' must be a positive number".to_string(),
            ))
        }
    };

    let band = match parse_u64(args, "band")? {
        None => 1,
        Some(v) if v >= 1 => v,
        Some(_) => return Err(ToolError::Validation("'band' must be >= 1".to_string())),
    };

    let max_obs = match parse_u64(args, "max_obs")? {
        None => 50usize,
        Some(v) if (1..=200).contains(&v) => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "'max_obs' must be between 1 and 200".to_string(),
            ))
        }
    };

    Ok(Params {
        correlation_length,
        background_error_variance,
        obs_error_variance,
        band,
        max_obs,
    })
}

fn parse_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
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

fn parse_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64().or_else(|| n.as_f64().map(|f| f as u64))),
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
    use wbraster::{CrsInfo, RasterConfig};
    use wbvector::{memory_store, FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a background raster from a closure `f(x, y)` sampled at cell
    /// centres, over `rows×cols` cells of size 1 starting at the origin.
    fn bg_raster(rows: usize, cols: usize, f: impl Fn(f64, f64) -> f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: f64::NAN,
            data_type: DataType::F32,
            crs: CrsInfo::from_epsg(3857),
            metadata: Vec::new(),
        });
        let y_max = rows as f64;
        for row in 0..rows {
            let y = y_max - (row as f64 + 0.5);
            for col in 0..cols {
                let x = col as f64 + 0.5;
                r.set(0, row as isize, col as isize, f(x, y)).unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn obs_layer(pts: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("obs")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (x, y, v) in pts {
            l.add_feature(Some(Geometry::point(*x, *y)), &[("v", (*v).into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = OptimalInterpolationTool.run(&args, &ctx()).unwrap();
        let r = crate::common::load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, r)
    }

    fn cell(r: &Raster, x: f64, y: f64) -> f64 {
        let col = ((x - r.x_min) / r.cell_size_x).floor() as isize;
        let row = ((r.y_max() - y) / r.cell_size_y).floor() as isize;
        r.get(0, row, col)
    }

    // ── Core: analysis pulls the background toward observations ───────────────

    /// With a constant background of 0 and observations all equal to 10 near
    /// the grid centre, the analysis rises toward 10 at cells near the obs and
    /// stays 0 far away — the defining behaviour of the OI increment.
    #[test]
    fn analysis_pulls_toward_observations() {
        let bg = bg_raster(40, 40, |_, _| 0.0);
        // A cluster of observations near (20, 20).
        let obs = obs_layer(&[
            (19.0, 19.0, 10.0),
            (21.0, 19.0, 10.0),
            (19.0, 21.0, 10.0),
            (21.0, 21.0, 10.0),
            (20.0, 20.0, 10.0),
        ]);
        let (_o, r) = run(json!({
            "background": bg, "observations": obs, "field": "v",
            "correlation_length": 5.0,
            "background_error_variance": 5.0, "obs_error_variance": 0.1,
        }));
        let near = cell(&r, 20.0, 20.0);
        let far = cell(&r, 2.0, 2.0);
        assert!(
            near > 5.0,
            "cell in the observation cluster should be pulled well above the background 0, got {near}"
        );
        assert!(
            far.abs() < 0.5,
            "cell far from any observation should stay ~= background 0, got {far}"
        );
    }

    /// Cells beyond `CUTOFF · L` of every observation are left exactly equal to
    /// the background (pass-through), and analysis error there equals σ_b².
    #[test]
    fn pass_through_far_from_observations() {
        let bg = bg_raster(60, 60, |x, y| 100.0 + x - y);
        let obs = obs_layer(&[(5.0, 5.0, 200.0)]); // one obs in a corner
        let (_o, r) = run(json!({
            "background": bg, "observations": obs, "field": "v",
            "correlation_length": 3.0, "background_error_variance": 4.0,
        }));
        // Far cell (col 55, row 55): centre (55.5, 4.5), background 100+55.5-4.5.
        let far = cell(&r, 55.0, 5.0);
        let bg_far = 100.0 + 55.5 - 4.5;
        assert!(
            (far - bg_far).abs() < 1e-9,
            "cell far from the observation must equal the background {bg_far}, got {far}"
        );
    }

    /// Analysis error variance is ≤ σ_b² everywhere and equals σ_b² where no
    /// observation is in range (assimilation never increases uncertainty).
    #[test]
    fn analysis_error_bounded_by_background_variance() {
        let bg = bg_raster(30, 30, |_, _| 5.0);
        let obs = obs_layer(&[(15.0, 15.0, 9.0), (16.0, 15.0, 9.0), (15.0, 16.0, 9.0)]);
        let sigma_b2 = 3.0;
        let ae_path = std::env::temp_dir()
            .join("geolibre_oi_ae_test.tif")
            .to_string_lossy()
            .into_owned();
        let args: ToolArgs = serde_json::from_value(json!({
            "background": bg, "observations": obs, "field": "v",
            "correlation_length": 4.0, "background_error_variance": sigma_b2,
            "obs_error_variance": 0.5,
            "analysis_error": ae_path,
        }))
        .unwrap();
        let out = OptimalInterpolationTool.run(&args, &ctx()).unwrap();
        let ae = crate::common::load_input_raster(out.outputs["analysis_error"].as_str().unwrap())
            .unwrap();
        let mut saw_reduced = false;
        for row in 0..ae.rows as isize {
            for col in 0..ae.cols as isize {
                let v = ae.get(0, row, col);
                if v.is_nan() {
                    continue;
                }
                assert!(
                    v <= sigma_b2 + 1e-9,
                    "analysis error {v} must not exceed background variance {sigma_b2}"
                );
                assert!(v >= -1e-9, "analysis error {v} must be non-negative");
                if v < sigma_b2 - 1e-6 {
                    saw_reduced = true;
                }
            }
        }
        // The cell at the observation cluster must show a real reduction.
        assert!(
            saw_reduced,
            "at least one near-obs cell should have reduced error"
        );
    }

    /// Trusting observations infinitely more than the background (tiny σ_o²,
    /// large σ_b²) makes the analysis at a lone observation match the observed
    /// value: x_a ≈ y there.
    #[test]
    fn strong_trust_reproduces_observation() {
        let bg = bg_raster(40, 40, |_, _| 0.0);
        let obs = obs_layer(&[(20.5, 20.5, 42.0)]);
        let (_o, r) = run(json!({
            "background": bg, "observations": obs, "field": "v",
            "correlation_length": 6.0,
            "background_error_variance": 1000.0, "obs_error_variance": 1e-6,
        }));
        let at = cell(&r, 20.5, 20.5);
        assert!(
            (at - 42.0).abs() < 1.0,
            "with near-perfect observations the analysis at the obs should be ~42, got {at}"
        );
    }

    /// A known 1-observation OI result, computed by hand, is reproduced exactly.
    /// Single obs: increment = c·(P+R)⁻¹·d where P+R = σ_b² + σ_o², c at the
    /// obs cell = σ_b² (distance 0 → correlation 1), d = y − x_b.
    #[test]
    fn single_obs_matches_closed_form() {
        // Background constant 3, one obs value 8 exactly at a cell centre.
        let bg = bg_raster(21, 21, |_, _| 3.0);
        let obs = obs_layer(&[(10.5, 10.5, 8.0)]);
        let sigma_b2 = 2.0;
        let sigma_o2 = 0.5;
        let (_o, r) = run(json!({
            "background": bg, "observations": obs, "field": "v",
            "correlation_length": 4.0,
            "background_error_variance": sigma_b2, "obs_error_variance": sigma_o2,
        }));
        let d = 8.0 - 3.0;
        let expected = 3.0 + sigma_b2 / (sigma_b2 + sigma_o2) * d; // c=σ_b², P+R=σ_b²+σ_o²
        let got = cell(&r, 10.5, 10.5);
        assert!(
            (got - expected).abs() < 1e-4,
            "closed-form single-obs analysis {expected}, got {got}"
        );
    }

    // ── Cholesky sanity ──────────────────────────────────────────────────────

    #[test]
    fn cholesky_solves_known_system() {
        // A = [[4,1,0],[1,3,1],[0,1,2]] SPD; x = [1,2,3]; b = A x.
        let a = [4.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 2.0];
        let x_true = [1.0, 2.0, 3.0];
        let mut b = [0.0; 3];
        for i in 0..3 {
            for j in 0..3 {
                b[i] += a[i * 3 + j] * x_true[j];
            }
        }
        let chol = cholesky(&a, 3).unwrap();
        let x = chol.solve(&b);
        for i in 0..3 {
            assert!((x[i] - x_true[i]).abs() < 1e-9, "component {i}: {}", x[i]);
        }
    }

    #[test]
    fn cholesky_rejects_non_spd() {
        let a = [1.0, 2.0, 2.0, 1.0]; // indefinite
        assert!(cholesky(&a, 2).is_none());
    }

    // ── Parameter validation ─────────────────────────────────────────────────

    #[test]
    fn rejects_bad_parameters() {
        let bg = bg_raster(10, 10, |_, _| 0.0);
        let obs = obs_layer(&[(5.0, 5.0, 1.0)]);
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            OptimalInterpolationTool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing everything");
        assert!(
            bad(json!({ "background": bg, "observations": obs })).is_err(),
            "missing field"
        );
        assert!(
            bad(json!({ "background": bg, "observations": obs, "field": "v", "correlation_length": -1.0 })).is_err(),
            "negative correlation_length"
        );
        assert!(
            bad(json!({ "background": bg, "observations": obs, "field": "v", "background_error_variance": 0.0 })).is_err(),
            "zero background variance"
        );
        assert!(
            bad(json!({ "background": bg, "observations": obs, "field": "v", "obs_error_variance": -2.0 })).is_err(),
            "negative obs variance"
        );
        assert!(
            bad(json!({ "background": bg, "observations": obs, "field": "v", "max_obs": 0 }))
                .is_err(),
            "max_obs zero"
        );
        assert!(
            bad(json!({ "background": bg, "observations": obs, "field": "v", "band": 0 })).is_err(),
            "band zero"
        );
        assert!(
            bad(json!({ "background": bg, "observations": obs, "field": "v" })).is_ok(),
            "minimal valid args"
        );
    }

    #[test]
    fn rejects_missing_field() {
        let bg = bg_raster(10, 10, |_, _| 0.0);
        let obs = obs_layer(&[(5.0, 5.0, 1.0)]);
        let args: ToolArgs = serde_json::from_value(
            json!({ "background": bg, "observations": obs, "field": "nope" }),
        )
        .unwrap();
        assert!(OptimalInterpolationTool.run(&args, &ctx()).is_err());
    }
}
