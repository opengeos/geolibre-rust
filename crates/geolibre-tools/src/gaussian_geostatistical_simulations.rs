//! GeoLibre tool: Gaussian geostatistical simulation (sequential Gaussian).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Gaussian Geostatistical Simulations*
//! (Geostatistical Analyst). Kriging returns one smoothed surface; this tool
//! produces `num_realizations` equiprobable **conditional** realizations of the
//! field, each honouring the sample data, for uncertainty propagation and
//! Monte-Carlo analysis.
//!
//! The authored suite has `empirical_bayesian_kriging` and the bundled ordinary/
//! universal/simple kriging — all single smoothed predictions — plus the
//! **unconditional** `turning_bands_simulation`. None produces conditional
//! multi-realization output.
//!
//! Method — Sequential Gaussian Simulation (SGS):
//! 1. Fit a semivariogram (`exponential`, `spherical`, or `gaussian`) to the
//!    conditioning points, or use manual `nugget`/`sill`/`range`.
//! 2. For each realization, visit the grid cells in a seeded random order. At
//!    each cell, simple-krige from the nearest `max_neighbors` known values
//!    (conditioning points **and** cells already simulated in this realization)
//!    to get an estimate and variance, then draw N(estimate, variance) with a
//!    seeded RNG and add the cell to the known set.
//!
//! Output is an `num_realizations`-band raster (one band per realization);
//! optional `output_mean` / `output_std` write the ensemble mean and standard
//! deviation. Deterministic: the same `seed` gives identical output (no
//! `Date::now`, no thread RNG — WASM-safe). `cell_size` and coordinates are in
//! the layer's CRS units.

use std::collections::BTreeMap;

use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};

use crate::common::{parse_optional_output, write_or_store_output};
use crate::vector_common::load_input_layer;

const NODATA: f64 = -9999.0;
const MAX_CELLS: usize = 4_000_000;

pub struct GaussianGeostatisticalSimulationsTool;

impl Tool for GaussianGeostatisticalSimulationsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "gaussian_geostatistical_simulations",
            display_name: "Gaussian Geostatistical Simulations",
            summary: "Sequential Gaussian simulation: N equiprobable conditional realizations of a field that honour the sample data (like ArcGIS Gaussian Geostatistical Simulations) — the conditional multi-realization output the single-surface empirical_bayesian_kriging / ordinary kriging and the unconditional turning_bands_simulation don't provide. Seeded and deterministic.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point layer with the conditioning values (projected CRS).",
                    required: true,
                },
                ToolParamSpec {
                    name: "value_field",
                    description: "Numeric field holding the value to simulate.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output multi-band raster (one band per realization). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "num_realizations",
                    description: "Number of realizations / output bands (default 10).",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Output cell size in CRS units. Default: extent's longer side / 100.",
                    required: false,
                },
                ToolParamSpec {
                    name: "variogram_model",
                    description: "'exponential' (default), 'spherical', or 'gaussian'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "nugget",
                    description: "Manual nugget; omit to fit from the data.",
                    required: false,
                },
                ToolParamSpec {
                    name: "sill",
                    description: "Manual partial sill; omit to fit from the data.",
                    required: false,
                },
                ToolParamSpec {
                    name: "range",
                    description: "Manual range; omit to fit from the data.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_neighbors",
                    description: "Max known points used per cell's kriging (default 16).",
                    required: false,
                },
                ToolParamSpec {
                    name: "seed",
                    description: "Seed for the deterministic RNG (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_mean",
                    description: "Optional raster path for the ensemble mean.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_std",
                    description: "Optional raster path for the ensemble standard deviation.",
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
        if args
            .get("value_field")
            .and_then(Value::as_str)
            .map(str::trim)
            .unwrap_or("")
            .is_empty()
        {
            return Err(ToolError::Validation(
                "missing required string parameter 'value_field'".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args.get("input").and_then(Value::as_str).unwrap();
        let value_field = args.get("value_field").and_then(Value::as_str).unwrap();
        let output = parse_optional_output(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let vidx = layer.schema.field_index(value_field).ok_or_else(|| {
            ToolError::Validation(format!("value_field '{value_field}' not found"))
        })?;

        // Conditioning points.
        let mut pts: Vec<Pt> = Vec::new();
        for feature in layer.features.iter() {
            let Some((x, y)) = feature.geometry.as_ref().and_then(point_xy) else {
                continue;
            };
            let Some(v) = feature.attributes.get(vidx).and_then(|f| f.as_f64()) else {
                continue;
            };
            pts.push(Pt { x, y, v });
        }
        if pts.len() < 3 {
            return Err(ToolError::Execution(format!(
                "need at least 3 valued points, found {}",
                pts.len()
            )));
        }

        // Extent & grid.
        let (mut min_x, mut min_y) = (f64::INFINITY, f64::INFINITY);
        let (mut max_x, mut max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
        for p in &pts {
            min_x = min_x.min(p.x);
            min_y = min_y.min(p.y);
            max_x = max_x.max(p.x);
            max_y = max_y.max(p.y);
        }
        let span = (max_x - min_x).max(max_y - min_y).max(1e-9);
        let cell = prm.cell_size.unwrap_or(span / 100.0);
        let cols = (((max_x - min_x) / cell).ceil() as usize + 1).max(1);
        let rows = (((max_y - min_y) / cell).ceil() as usize + 1).max(1);
        if rows.saturating_mul(cols) > MAX_CELLS {
            return Err(ToolError::Validation(format!(
                "grid {rows}x{cols} exceeds {MAX_CELLS} cells; increase cell_size"
            )));
        }

        // Variogram: fit or manual.
        let data_mean = pts.iter().map(|p| p.v).sum::<f64>() / pts.len() as f64;
        let data_var =
            pts.iter().map(|p| (p.v - data_mean).powi(2)).sum::<f64>() / pts.len().max(1) as f64;
        let vg = build_variogram(&pts, &prm, data_var, span);
        ctx.progress.info(&format!(
            "{} points, {rows}x{cols} grid, {} realization(s); variogram {} (nugget {:.3}, sill {:.3}, range {:.3})",
            pts.len(),
            prm.num_realizations,
            vg.model.name(),
            vg.nugget,
            vg.sill,
            vg.range
        ));

        let sill_total = vg.nugget + vg.sill; // C(0)
        let search_r2 = {
            let r = (vg.range * 3.0).max(span);
            r * r
        };

        // Grid cell centres, row 0 = north.
        let cell_xy = |idx: usize| -> (f64, f64) {
            let r = idx / cols;
            let c = idx % cols;
            (
                min_x + (c as f64 + 0.5) * cell,
                max_y - (r as f64 + 0.5) * cell,
            )
        };
        let n_cells = rows * cols;

        // ── Simulate each realization ─────────────────────────────────────────
        let mut bands: Vec<Vec<f64>> = Vec::with_capacity(prm.num_realizations);
        for real in 0..prm.num_realizations {
            let mut rng = Rng::new(prm.seed ^ ((real as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)));

            // Known set starts as the conditioning points; grows as we simulate.
            let mut kx: Vec<f64> = pts.iter().map(|p| p.x).collect();
            let mut ky: Vec<f64> = pts.iter().map(|p| p.y).collect();
            let mut kv: Vec<f64> = pts.iter().map(|p| p.v).collect();
            let mut tree: KdTree<f64, usize, [f64; 2]> = KdTree::new(2);
            for i in 0..kx.len() {
                tree.add([kx[i], ky[i]], i).ok();
            }

            // Seeded random visiting path over grid cells.
            let mut path: Vec<usize> = (0..n_cells).collect();
            fisher_yates(&mut path, &mut rng);

            let mut grid = vec![NODATA; n_cells];
            for &cell_idx in &path {
                let (qx, qy) = cell_xy(cell_idx);
                // Nearest known points within the search radius.
                let found = tree
                    .nearest(&[qx, qy], prm.max_neighbors, &squared_euclidean)
                    .unwrap_or_default();
                let nbrs: Vec<usize> = found
                    .into_iter()
                    .filter(|(d2, _)| *d2 <= search_r2)
                    .map(|(_, &i)| i)
                    .collect();

                let sim = if nbrs.is_empty() {
                    // No neighbours in range: draw from the marginal.
                    data_mean + sill_total.max(0.0).sqrt() * rng.normal()
                } else {
                    let (est, var) =
                        simple_krige(&nbrs, &kx, &ky, &kv, qx, qy, data_mean, sill_total, &vg);
                    est + var.max(0.0).sqrt() * rng.normal()
                };

                grid[cell_idx] = sim;
                // Add the simulated cell to the known set.
                let id = kx.len();
                kx.push(qx);
                ky.push(qy);
                kv.push(sim);
                tree.add([qx, qy], id).ok();
            }
            bands.push(grid);
            ctx.progress
                .progress((real as f64 + 1.0) / prm.num_realizations as f64);
        }

        // ── Assemble outputs ──────────────────────────────────────────────────
        let epsg = layer.crs_epsg();
        let make_raster = |data_bands: &[Vec<f64>]| -> Result<Raster, ToolError> {
            let mut r = Raster::new(RasterConfig {
                cols,
                rows,
                bands: data_bands.len(),
                x_min: min_x - 0.5 * cell,
                y_min: max_y - rows as f64 * cell + 0.5 * cell,
                cell_size: cell,
                cell_size_y: Some(cell),
                nodata: NODATA,
                data_type: DataType::F32,
                crs: CrsInfo {
                    epsg,
                    wkt: None,
                    proj4: None,
                },
                metadata: Vec::new(),
            });
            for (bi, band) in data_bands.iter().enumerate() {
                for row in 0..rows {
                    for col in 0..cols {
                        r.set(
                            bi as isize,
                            row as isize,
                            col as isize,
                            band[row * cols + col],
                        )
                        .map_err(|e| ToolError::Execution(format!("failed writing cell: {e}")))?;
                    }
                }
            }
            Ok(r)
        };

        let out_path = write_or_store_output(make_raster(&bands)?, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("realizations".to_string(), json!(prm.num_realizations));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        outputs.insert("variogram_model".to_string(), json!(vg.model.name()));
        outputs.insert("nugget".to_string(), json!(vg.nugget));
        outputs.insert("sill".to_string(), json!(vg.sill));
        outputs.insert("range".to_string(), json!(vg.range));

        // Optional ensemble mean / std.
        if prm.output_mean.is_some() || prm.output_std.is_some() {
            let mut mean = vec![0.0; n_cells];
            let mut m2 = vec![0.0; n_cells];
            for band in &bands {
                for i in 0..n_cells {
                    mean[i] += band[i];
                }
            }
            for m in mean.iter_mut() {
                *m /= prm.num_realizations as f64;
            }
            for band in &bands {
                for i in 0..n_cells {
                    m2[i] += (band[i] - mean[i]).powi(2);
                }
            }
            let std: Vec<f64> = m2
                .iter()
                .map(|s| (s / prm.num_realizations as f64).sqrt())
                .collect();
            if let Some(path) = &prm.output_mean {
                let p = write_or_store_output(make_raster(&[mean.clone()])?, Some(path.as_str()))?;
                outputs.insert("output_mean".to_string(), json!(p));
            }
            if let Some(path) = &prm.output_std {
                let p = write_or_store_output(make_raster(&[std])?, Some(path.as_str()))?;
                outputs.insert("output_std".to_string(), json!(p));
            }
        }

        Ok(ToolRunResult { outputs })
    }
}

// ── Simple kriging ────────────────────────────────────────────────────────────

/// Simple kriging estimate and variance at (qx, qy) from the neighbour indices.
/// Covariance C(h) = sill_total − γ(h); diagonal uses C(0) = sill_total.
#[allow(clippy::too_many_arguments)]
fn simple_krige(
    nbrs: &[usize],
    kx: &[f64],
    ky: &[f64],
    kv: &[f64],
    qx: f64,
    qy: f64,
    mean: f64,
    sill_total: f64,
    vg: &Vg,
) -> (f64, f64) {
    let k = nbrs.len();
    let cov = |h: f64| -> f64 {
        if h <= 0.0 {
            sill_total
        } else {
            sill_total - gamma(vg, h)
        }
    };
    // Covariance matrix among neighbours (+ tiny ridge) and the RHS to the query.
    let mut a = vec![0.0f64; k * k];
    let mut b = vec![0.0f64; k];
    for i in 0..k {
        let (xi, yi) = (kx[nbrs[i]], ky[nbrs[i]]);
        for j in 0..k {
            let (xj, yj) = (kx[nbrs[j]], ky[nbrs[j]]);
            let h = (xi - xj).hypot(yi - yj);
            let mut c = cov(h);
            if i == j {
                c += 1e-9 * sill_total.max(1.0);
            }
            a[i * k + j] = c;
        }
        let hq = (xi - qx).hypot(yi - qy);
        b[i] = cov(hq);
    }
    let w = match solve(&a, &b, k) {
        Some(w) => w,
        None => {
            // Fall back to the marginal if the system is singular.
            return (mean, sill_total);
        }
    };
    let mut est = mean;
    let mut var = sill_total;
    for i in 0..k {
        est += w[i] * (kv[nbrs[i]] - mean);
        var -= w[i] * b[i];
    }
    (est, var)
}

// ── Variogram ─────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Model {
    Exponential,
    Spherical,
    Gaussian,
}

impl Model {
    fn name(&self) -> &'static str {
        match self {
            Model::Exponential => "exponential",
            Model::Spherical => "spherical",
            Model::Gaussian => "gaussian",
        }
    }
}

#[derive(Clone, Copy)]
struct Vg {
    model: Model,
    nugget: f64,
    sill: f64, // partial sill
    range: f64,
}

/// Structure function g(h): 0 at h=0, →1 far. γ(h) = nugget + sill·g(h).
fn structure(model: Model, h: f64, range: f64) -> f64 {
    let r = range.max(1e-9);
    match model {
        Model::Exponential => 1.0 - (-h / r).exp(),
        Model::Gaussian => 1.0 - (-(h * h) / (r * r)).exp(),
        Model::Spherical => {
            if h >= r {
                1.0
            } else {
                let t = h / r;
                1.5 * t - 0.5 * t * t * t
            }
        }
    }
}

fn gamma(vg: &Vg, h: f64) -> f64 {
    if h <= 0.0 {
        0.0
    } else {
        vg.nugget + vg.sill * structure(vg.model, h, vg.range)
    }
}

/// Builds the variogram from manual params where given, else fits the chosen
/// model to the empirical semivariogram.
fn build_variogram(pts: &[Pt], prm: &Params, data_var: f64, span: f64) -> Vg {
    if let (Some(n), Some(s), Some(r)) = (prm.nugget, prm.sill, prm.range) {
        return Vg {
            model: prm.model,
            nugget: n.max(0.0),
            sill: s.max(1e-12),
            range: r.max(1e-9),
        };
    }
    let mut vg = fit_variogram(pts, prm.model, span);
    // Manual overrides for any individually-specified component.
    if let Some(n) = prm.nugget {
        vg.nugget = n.max(0.0);
    }
    if let Some(s) = prm.sill {
        vg.sill = s.max(1e-12);
    }
    if let Some(r) = prm.range {
        vg.range = r.max(1e-9);
    }
    if !(vg.sill.is_finite() && vg.sill > 0.0) {
        vg.sill = data_var.max(1e-9);
    }
    vg
}

/// Bins the empirical semivariogram, then fits nugget + partial-sill by WLS for
/// each candidate range (golden-section search over range).
fn fit_variogram(pts: &[Pt], model: Model, span: f64) -> Vg {
    let n = pts.len();
    // Empirical pairs (h, semivariance) with pair counts per distance bin.
    let mut pairs: Vec<(f64, f64)> = Vec::with_capacity(n * (n - 1) / 2);
    let mut max_h = 0.0f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let h = (pts[i].x - pts[j].x).hypot(pts[i].y - pts[j].y);
            let g = 0.5 * (pts[i].v - pts[j].v).powi(2);
            pairs.push((h, g));
            max_h = max_h.max(h);
        }
    }
    if pairs.is_empty() || max_h <= 0.0 {
        return Vg {
            model,
            nugget: 0.0,
            sill: 1.0,
            range: span.max(1.0),
        };
    }
    let nbins = 15usize;
    let bw = max_h / nbins as f64;
    let mut sum = vec![0.0f64; nbins];
    let mut cnt = vec![0usize; nbins];
    let mut hsum = vec![0.0f64; nbins];
    for &(h, g) in &pairs {
        let b = ((h / bw) as usize).min(nbins - 1);
        sum[b] += g;
        hsum[b] += h;
        cnt[b] += 1;
    }
    let bins: Vec<(f64, f64, f64)> = (0..nbins)
        .filter(|&b| cnt[b] > 0)
        .map(|b| {
            (
                hsum[b] / cnt[b] as f64,
                sum[b] / cnt[b] as f64,
                cnt[b] as f64,
            )
        })
        .collect();
    if bins.len() < 2 {
        let sill = pairs.iter().map(|p| p.1).sum::<f64>() / pairs.len() as f64;
        return Vg {
            model,
            nugget: 0.0,
            sill: sill.max(1e-9),
            range: (max_h / 3.0).max(1e-9),
        };
    }

    // WLS fit of nugget + sill·g(h; range) for a fixed range.
    let fit_for_range = |range: f64| -> (f64, f64, f64) {
        let (mut s_w, mut s_wx, mut s_wxx, mut s_wy, mut s_wxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for &(h, g, w) in &bins {
            let x = structure(model, h, range);
            s_w += w;
            s_wx += w * x;
            s_wxx += w * x * x;
            s_wy += w * g;
            s_wxy += w * x * g;
        }
        let ridge = 1e-9 * s_w.max(1.0);
        let a = [s_w + ridge, s_wx, s_wx, s_wxx + ridge];
        let b = [s_wy, s_wxy];
        let (nugget, sill) = match solve(&a, &b, 2) {
            Some(v) => (v[0].max(0.0), v[1].max(1e-12)),
            None => (0.0, 1e-6),
        };
        let mut sse = 0.0;
        for &(h, g, w) in &bins {
            let pred = nugget + sill * structure(model, h, range);
            sse += w * (g - pred).powi(2);
        }
        (sse, nugget, sill)
    };

    // Golden-section over range in [max_h/20, max_h].
    let (mut lo, mut hi) = (max_h / 20.0, max_h);
    let gr = (5f64.sqrt() - 1.0) / 2.0;
    let mut c = hi - gr * (hi - lo);
    let mut d = lo + gr * (hi - lo);
    let mut fc = fit_for_range(c).0;
    let mut fd = fit_for_range(d).0;
    for _ in 0..40 {
        if fc < fd {
            hi = d;
            d = c;
            fd = fc;
            c = hi - gr * (hi - lo);
            fc = fit_for_range(c).0;
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = lo + gr * (hi - lo);
            fd = fit_for_range(d).0;
        }
    }
    let range = 0.5 * (lo + hi);
    let (_, nugget, sill) = fit_for_range(range);
    Vg {
        model,
        nugget,
        sill,
        range,
    }
}

// ── Linear solver (Gauss elimination with partial pivoting) ───────────────────

fn solve(a_in: &[f64], b_in: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut a = a_in.to_vec();
    let mut b = b_in.to_vec();
    for col in 0..n {
        // Partial pivot.
        let mut piv = col;
        let mut best = a[col * n + col].abs();
        for r in (col + 1)..n {
            let v = a[r * n + col].abs();
            if v > best {
                best = v;
                piv = r;
            }
        }
        if best < 1e-12 {
            return None;
        }
        if piv != col {
            for k in 0..n {
                a.swap(col * n + k, piv * n + k);
            }
            b.swap(col, piv);
        }
        let d = a[col * n + col];
        for r in 0..n {
            if r == col {
                continue;
            }
            let f = a[r * n + col] / d;
            if f == 0.0 {
                continue;
            }
            for k in col..n {
                a[r * n + k] -= f * a[col * n + k];
            }
            b[r] -= f * b[col];
        }
    }
    Some((0..n).map(|i| b[i] / a[i * n + i]).collect())
}

// ── Deterministic RNG (splitmix64 + Box-Muller) ───────────────────────────────

struct Rng {
    state: u64,
    spare: Option<f64>,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Rng {
            state: seed ^ 0xDEAD_BEEF_CAFE_1234,
            spare: None,
        }
    }
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    /// Standard normal via Box-Muller (caches the second deviate).
    fn normal(&mut self) -> f64 {
        if let Some(s) = self.spare.take() {
            return s;
        }
        let u1 = self.next_f64().max(1e-12);
        let u2 = self.next_f64();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;
        self.spare = Some(r * theta.sin());
        r * theta.cos()
    }
    fn next_usize(&mut self, bound: usize) -> usize {
        (self.next_u64() % bound as u64) as usize
    }
}

fn fisher_yates(v: &mut [usize], rng: &mut Rng) {
    let n = v.len();
    for i in (1..n).rev() {
        let j = rng.next_usize(i + 1);
        v.swap(i, j);
    }
}

// ── Geometry / parameters ─────────────────────────────────────────────────────

struct Pt {
    x: f64,
    y: f64,
    v: f64,
}

fn point_xy(geom: &wbvector::Geometry) -> Option<(f64, f64)> {
    match geom {
        wbvector::Geometry::Point(c) => Some((c.x, c.y)),
        wbvector::Geometry::MultiPoint(cs) if !cs.is_empty() => {
            let n = cs.len() as f64;
            Some((
                cs.iter().map(|c| c.x).sum::<f64>() / n,
                cs.iter().map(|c| c.y).sum::<f64>() / n,
            ))
        }
        _ => None,
    }
}

struct Params {
    num_realizations: usize,
    cell_size: Option<f64>,
    model: Model,
    nugget: Option<f64>,
    sill: Option<f64>,
    range: Option<f64>,
    max_neighbors: usize,
    seed: u64,
    output_mean: Option<String>,
    output_std: Option<String>,
}

fn opt_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(Some(n.as_f64().unwrap_or(f64::NAN))),
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

fn opt_str(args: &ToolArgs, key: &str) -> Option<String> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_string)
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let num_realizations = match args.get("num_realizations") {
        None | Some(Value::Null) => 10,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(10).max(1) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 10,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'num_realizations' must be an integer".into()))?
            .max(1),
        Some(_) => {
            return Err(ToolError::Validation(
                "'num_realizations' must be a number".into(),
            ))
        }
    };
    let cell_size = match opt_f64(args, "cell_size")? {
        None => None,
        Some(v) if v > 0.0 => Some(v),
        Some(_) => return Err(ToolError::Validation("'cell_size' must be > 0".into())),
    };
    let model = match args
        .get("variogram_model")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some("exponential") => Model::Exponential,
        Some("spherical") => Model::Spherical,
        Some("gaussian") => Model::Gaussian,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'variogram_model' must be exponential|spherical|gaussian, got '{o}'"
            )))
        }
    };
    let nugget = opt_f64(args, "nugget")?;
    let sill = opt_f64(args, "sill")?;
    let range = opt_f64(args, "range")?;
    let max_neighbors = match args.get("max_neighbors") {
        None | Some(Value::Null) => 16,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(16).max(1) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 16,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'max_neighbors' must be an integer".into()))?
            .max(1),
        Some(_) => {
            return Err(ToolError::Validation(
                "'max_neighbors' must be a number".into(),
            ))
        }
    };
    let seed = match args.get("seed") {
        None | Some(Value::Null) => 1,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(1),
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
            .map_err(|_| ToolError::Validation("'seed' must be an integer".into()))?,
        Some(_) => return Err(ToolError::Validation("'seed' must be a number".into())),
    };
    Ok(Params {
        num_realizations,
        cell_size,
        model,
        nugget,
        sill,
        range,
        max_neighbors,
        seed,
        output_mean: opt_str(args, "output_mean"),
        output_std: opt_str(args, "output_std"),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, Geometry, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn layer_of(pts: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(32610);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for &(x, y, v) in pts {
            l.add_feature(Some(Geometry::point(x, y)), &[("v", v.into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    /// A planar trend sampled on a grid: realizations honour the data and the
    /// ensemble mean tracks the trend.
    fn trend_points() -> Vec<(f64, f64, f64)> {
        let mut v = Vec::new();
        for i in 0..7 {
            for j in 0..7 {
                let (x, y) = (i as f64 * 10.0, j as f64 * 10.0);
                v.push((x, y, 2.0 + 0.5 * x + 0.3 * y));
            }
        }
        v
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = GaussianGeostatisticalSimulationsTool
            .run(&args, &ctx())
            .unwrap();
        let r = crate::common::load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, r)
    }

    /// Same seed -> identical output; different seed -> different realizations.
    #[test]
    fn deterministic_by_seed() {
        let pts = trend_points();
        let (_o1, r1) = run(json!({
            "input": layer_of(&pts), "value_field": "v", "num_realizations": 2,
            "cell_size": 10.0, "seed": 42
        }));
        let (_o2, r2) = run(json!({
            "input": layer_of(&pts), "value_field": "v", "num_realizations": 2,
            "cell_size": 10.0, "seed": 42
        }));
        let (_o3, r3) = run(json!({
            "input": layer_of(&pts), "value_field": "v", "num_realizations": 2,
            "cell_size": 10.0, "seed": 7
        }));
        let mut same = true;
        let mut diff = false;
        for row in 0..r1.rows as isize {
            for col in 0..r1.cols as isize {
                if (r1.get(0, row, col) - r2.get(0, row, col)).abs() > 1e-9 {
                    same = false;
                }
                if (r1.get(0, row, col) - r3.get(0, row, col)).abs() > 1e-9 {
                    diff = true;
                }
            }
        }
        assert!(same, "same seed must reproduce identical output");
        assert!(diff, "a different seed must change the realization");
    }

    /// Conditional simulation honours the data: at a sample location the mean of
    /// many realizations is close to the sampled value (low nugget).
    #[test]
    fn honours_conditioning_data() {
        let pts = trend_points();
        let (out, r) = run(json!({
            "input": layer_of(&pts), "value_field": "v", "num_realizations": 30,
            "cell_size": 10.0, "nugget": 0.0, "sill": 100.0, "range": 40.0,
            "variogram_model": "exponential", "seed": 3
        }));
        // Sample at (30,30) -> value 2 + 15 + 9 = 26. Find the covering cell.
        let cols = out.outputs["cols"].as_u64().unwrap() as isize;
        let rows = out.outputs["rows"].as_u64().unwrap() as isize;
        // The grid origin: x_min-0.5cell etc; cell centre col c is at 0 + c*10.
        // (30,30) is a data location; the nearest cell centre is col 3, from top.
        let c = 3isize;
        let row = rows - 1 - 3; // y=30 is 3 cells up from the bottom
        let mut sum = 0.0;
        let bands = out.outputs["realizations"].as_u64().unwrap() as isize;
        for b in 0..bands {
            sum += r.get(b, row.clamp(0, rows - 1), c.clamp(0, cols - 1));
        }
        let mean = sum / bands as f64;
        assert!(
            (mean - 26.0).abs() < 6.0,
            "ensemble mean at a data location ({mean}) should track the datum 26"
        );
    }

    /// Ensemble variance is small at/near data locations and larger away from
    /// them — the signature of conditional simulation.
    #[test]
    fn variance_grows_away_from_data() {
        // A single dense cluster near the origin; the far side of the extent has
        // no data, so its simulated variance is larger.
        let mut pts = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                pts.push((i as f64 * 5.0, j as f64 * 5.0, 10.0));
            }
        }
        // one far anchor so the grid spans a wide gap.
        pts.push((300.0, 300.0, 10.0));
        let (out, r) = run(json!({
            "input": layer_of(&pts), "value_field": "v", "num_realizations": 25,
            "cell_size": 15.0, "nugget": 0.0, "sill": 50.0, "range": 20.0, "seed": 5
        }));
        let cols = out.outputs["cols"].as_u64().unwrap() as isize;
        let rows = out.outputs["rows"].as_u64().unwrap() as isize;
        let bands = out.outputs["realizations"].as_u64().unwrap() as isize;
        let var_at = |row: isize, col: isize| -> f64 {
            let vals: Vec<f64> = (0..bands).map(|b| r.get(b, row, col)).collect();
            let m = vals.iter().sum::<f64>() / bands as f64;
            vals.iter().map(|v| (v - m).powi(2)).sum::<f64>() / bands as f64
        };
        // Near the dense cluster (bottom-left) vs the empty middle of the extent.
        let near = var_at(rows - 1, 0);
        let far = var_at(rows / 2, cols / 2);
        assert!(
            far > near,
            "variance should grow away from data (near {near}, far {far})"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            GaussianGeostatisticalSimulationsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err()); // no value_field
        assert!(bad(
            json!({ "input": "a.geojson", "value_field": "v", "variogram_model": "bogus" })
        )
        .is_err());
        assert!(bad(json!({ "input": "a.geojson", "value_field": "v" })).is_ok());
    }
}
