//! GeoLibre tool: slope, aspect, and curvatures from a DEM via a local
//! quadratic/biquadratic surface fit over a chosen neighborhood distance.
//!
//! Pure-Rust counterpart of ArcGIS's *Surface Parameters* (Spatial Analyst /
//! 3D Analyst). The bundled whitebox suite ships fixed 3×3 slope, aspect, and
//! curvature tools; this one fits a least-squares polynomial surface over a
//! neighborhood of a user-chosen radius, so the parameters can be computed at a
//! coarser scale than a single cell — the distinguishing capability of the
//! ArcGIS tool.
//!
//! For each cell the tool gathers the valid cells within `neighborhood_distance`
//! (map units, converted to a cell radius), fits `z = ax² + by² + cxy + dx +
//! ey + f` (`quadratic`) or the 9-term `biquadratic` extension by ordinary
//! least squares, and reads the partial derivatives at the cell center
//! (`fx=d, fy=e, fxx=2a, fyy=2b, fxy=c`). The requested parameter is derived
//! from those with the standard second-order surface formulas. Cells with too
//! few valid neighbors to fit the model, and no-data cells, are left no-data.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};

const NODATA_OUT: f64 = -9999.0;

pub struct SurfaceParametersTool;

impl Tool for SurfaceParametersTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "surface_parameters",
            display_name: "Surface Parameters",
            summary: "Slope, aspect, or curvature (mean/profile/tangential/plan/gaussian/casorati/contour_geodesic_torsion) from a DEM via a local quadratic/biquadratic least-squares surface fit over a chosen neighborhood distance, like ArcGIS Surface Parameters — the scale-adaptive counterpart to the bundled fixed 3×3 terrain tools.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input elevation raster (DEM).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output raster path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "parameter",
                    description: "Which parameter to compute: slope (default), aspect, mean_curvature, profile_curvature, tangential_curvature, plan_curvature, gaussian_curvature, casorati_curvature, contour_geodesic_torsion.",
                    required: false,
                },
                ToolParamSpec {
                    name: "fit",
                    description: "Surface model: 'quadratic' (6-term, default) or 'biquadratic' (9-term).",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighborhood_distance",
                    description: "Neighborhood radius in map units (default = one cell size, i.e. a 3×3 window).",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighborhood_type",
                    description: "'fixed' (default) uses the full radius everywhere; 'adaptive' shrinks the window to 1 cell next to data gaps/edges.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_parameter(args)?;
        parse_fit(args)?;
        parse_neighborhood_type(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_output(args, "output")?;
        let param = parse_parameter(args)?;
        let fit = parse_fit(args)?;
        let adaptive = matches!(parse_neighborhood_type(args)?, Neighborhood::Adaptive);
        let band_1 = parse_band(args)?;

        let raster = load_input_raster(input)?;
        if (band_1 as usize) > raster.bands {
            return Err(ToolError::Validation(format!(
                "band {band_1} out of range (raster has {} band(s))",
                raster.bands
            )));
        }
        let band = (band_1 - 1) as isize;
        let nodata = raster.nodata;
        let rows = raster.rows;
        let cols = raster.cols;
        let dx = raster.cell_size_x.abs();
        let dy = raster.cell_size_y.abs();
        if dx <= 0.0 || dy <= 0.0 {
            return Err(ToolError::Execution(
                "raster has a non-positive cell size".to_string(),
            ));
        }

        let cell = 0.5 * (dx + dy);
        let dist = parse_optional_f64(args, "neighborhood_distance")?.unwrap_or(cell);
        let radius = ((dist / cell).round() as isize).max(1);
        let n_terms = if fit == Fit::Biquadratic { 9 } else { 6 };

        // Dense value buffer for fast neighbor access.
        let mut z = vec![f64::NAN; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let v = raster.get(band, r as isize, c as isize);
                if v != nodata && v.is_finite() {
                    z[r * cols + c] = v;
                }
            }
        }

        ctx.progress
            .info(&format!("{} at radius {radius} cell(s)", param.label()));

        let mut out = vec![NODATA_OUT; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                if z[r * cols + c].is_nan() {
                    continue;
                }
                let rad = if adaptive {
                    adaptive_radius(&z, rows, cols, r, c, radius)
                } else {
                    radius
                };
                // Gather local (x, y, z) samples in map units, origin at the cell.
                let mut xs = Vec::new();
                let mut ys = Vec::new();
                let mut zs = Vec::new();
                for dr in -rad..=rad {
                    for dc in -rad..=rad {
                        let rr = r as isize + dr;
                        let cc = c as isize + dc;
                        if rr < 0 || cc < 0 || rr >= rows as isize || cc >= cols as isize {
                            continue;
                        }
                        let val = z[rr as usize * cols + cc as usize];
                        if val.is_nan() {
                            continue;
                        }
                        xs.push(dc as f64 * dx);
                        ys.push(-dr as f64 * dy); // raster row increases downward
                        zs.push(val);
                    }
                }
                if zs.len() < n_terms {
                    continue;
                }
                let Some(coef) = fit_surface(&xs, &ys, &zs, n_terms) else {
                    continue;
                };
                // Center derivatives: fx=d, fy=e, fxx=2a, fyy=2b, fxy=c.
                let (fxx, fyy, fxy, fx, fy) =
                    (2.0 * coef[0], 2.0 * coef[1], coef[2], coef[3], coef[4]);
                out[r * cols + c] = param.value(fx, fy, fxx, fyy, fxy);
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let out_r = raster_like_with_data(&raster, out, NODATA_OUT, DataType::F32)?;
        let out_path = write_or_store_output(out_r, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("parameter".to_string(), json!(param.label()));
        outputs.insert("radius_cells".to_string(), json!(radius));
        Ok(ToolRunResult { outputs })
    }
}

/// Shrinks the radius to the largest value (≤ `radius`) whose square window
/// contains only valid cells, so parameters stay defined next to data gaps and
/// the raster edge. Returns at least 1.
fn adaptive_radius(
    z: &[f64],
    rows: usize,
    cols: usize,
    r: usize,
    c: usize,
    radius: isize,
) -> isize {
    let mut rad = radius;
    while rad > 1 {
        let mut ok = true;
        'scan: for dr in -rad..=rad {
            for dc in -rad..=rad {
                let rr = r as isize + dr;
                let cc = c as isize + dc;
                if rr < 0
                    || cc < 0
                    || rr >= rows as isize
                    || cc >= cols as isize
                    || z[rr as usize * cols + cc as usize].is_nan()
                {
                    ok = false;
                    break 'scan;
                }
            }
        }
        if ok {
            break;
        }
        rad -= 1;
    }
    rad
}

/// Fits `z ≈ Σ coef_k · term_k(x, y)` by ordinary least squares and returns the
/// coefficient vector ordered `[x², y², xy, x, y, 1, (x²y, xy², x²y²)]`. Returns
/// `None` if the normal equations are singular.
fn fit_surface(xs: &[f64], ys: &[f64], zs: &[f64], n_terms: usize) -> Option<Vec<f64>> {
    let terms = |x: f64, y: f64| -> [f64; 9] {
        [
            x * x,
            y * y,
            x * y,
            x,
            y,
            1.0,
            x * x * y,
            x * y * y,
            x * x * y * y,
        ]
    };
    // Normal equations A^T A c = A^T z.
    let mut ata = vec![vec![0.0_f64; n_terms]; n_terms];
    let mut atz = vec![0.0_f64; n_terms];
    for i in 0..zs.len() {
        let t = terms(xs[i], ys[i]);
        for a in 0..n_terms {
            atz[a] += t[a] * zs[i];
            for b in 0..n_terms {
                ata[a][b] += t[a] * t[b];
            }
        }
    }
    solve(ata, atz)
}

/// Gaussian elimination with partial pivoting for a small dense system.
#[allow(clippy::needless_range_loop)]
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        // Pivot.
        let mut piv = col;
        let mut best = a[col][col].abs();
        for r in (col + 1)..n {
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
        // Eliminate.
        for r in (col + 1)..n {
            let f = a[r][col] / a[col][col];
            if f != 0.0 {
                for k in col..n {
                    a[r][k] -= f * a[col][k];
                }
                b[r] -= f * b[col];
            }
        }
    }
    // Back-substitution.
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let mut s = b[row];
        for k in (row + 1)..n {
            s -= a[row][k] * x[k];
        }
        x[row] = s / a[row][row];
    }
    Some(x)
}

// ── Parameter derivation ────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Parameter {
    Slope,
    Aspect,
    MeanCurvature,
    ProfileCurvature,
    TangentialCurvature,
    PlanCurvature,
    GaussianCurvature,
    CasoratiCurvature,
    ContourGeodesicTorsion,
}

impl Parameter {
    fn label(self) -> &'static str {
        match self {
            Parameter::Slope => "slope",
            Parameter::Aspect => "aspect",
            Parameter::MeanCurvature => "mean_curvature",
            Parameter::ProfileCurvature => "profile_curvature",
            Parameter::TangentialCurvature => "tangential_curvature",
            Parameter::PlanCurvature => "plan_curvature",
            Parameter::GaussianCurvature => "gaussian_curvature",
            Parameter::CasoratiCurvature => "casorati_curvature",
            Parameter::ContourGeodesicTorsion => "contour_geodesic_torsion",
        }
    }

    fn value(self, fx: f64, fy: f64, fxx: f64, fyy: f64, fxy: f64) -> f64 {
        let p = fx;
        let q = fy;
        let pq = p * p + q * q;
        let one_pq = 1.0 + pq;
        match self {
            Parameter::Slope => pq.sqrt().atan().to_degrees(),
            Parameter::Aspect => {
                if pq < 1e-18 {
                    -1.0
                } else {
                    // Downslope azimuth, degrees clockwise from north.
                    let mut a = (-p).atan2(-q).to_degrees();
                    if a < 0.0 {
                        a += 360.0;
                    }
                    a
                }
            }
            Parameter::GaussianCurvature => (fxx * fyy - fxy * fxy) / (one_pq * one_pq),
            Parameter::MeanCurvature => {
                -((1.0 + q * q) * fxx - 2.0 * p * q * fxy + (1.0 + p * p) * fyy)
                    / (2.0 * one_pq.powf(1.5))
            }
            Parameter::ProfileCurvature => {
                if pq < 1e-18 {
                    0.0
                } else {
                    -(fxx * p * p + 2.0 * fxy * p * q + fyy * q * q) / (pq * one_pq.powf(1.5))
                }
            }
            Parameter::PlanCurvature => {
                if pq < 1e-18 {
                    0.0
                } else {
                    -(fxx * q * q - 2.0 * fxy * p * q + fyy * p * p) / pq.powf(1.5)
                }
            }
            Parameter::TangentialCurvature => {
                if pq < 1e-18 {
                    0.0
                } else {
                    -(fxx * q * q - 2.0 * fxy * p * q + fyy * p * p) / (pq * one_pq.sqrt())
                }
            }
            Parameter::CasoratiCurvature => {
                // From principal curvatures via mean (H) and gaussian (K).
                let h = -((1.0 + q * q) * fxx - 2.0 * p * q * fxy + (1.0 + p * p) * fyy)
                    / (2.0 * one_pq.powf(1.5));
                let k = (fxx * fyy - fxy * fxy) / (one_pq * one_pq);
                let disc = (h * h - k).max(0.0).sqrt();
                let k1 = h + disc;
                let k2 = h - disc;
                ((k1 * k1 + k2 * k2) / 2.0).sqrt()
            }
            Parameter::ContourGeodesicTorsion => {
                if pq < 1e-18 {
                    0.0
                } else {
                    ((fyy - fxx) * p * q + fxy * (p * p - q * q)) / (pq * one_pq)
                }
            }
        }
    }
}

#[derive(PartialEq, Eq)]
enum Fit {
    Quadratic,
    Biquadratic,
}

enum Neighborhood {
    Fixed,
    Adaptive,
}

fn parse_parameter(args: &ToolArgs) -> Result<Parameter, ToolError> {
    Ok(
        match args.get("parameter").and_then(Value::as_str).map(str::trim) {
            None | Some("") | Some("slope") => Parameter::Slope,
            Some("aspect") => Parameter::Aspect,
            Some("mean_curvature") => Parameter::MeanCurvature,
            Some("profile_curvature") => Parameter::ProfileCurvature,
            Some("tangential_curvature") => Parameter::TangentialCurvature,
            Some("plan_curvature") => Parameter::PlanCurvature,
            Some("gaussian_curvature") => Parameter::GaussianCurvature,
            Some("casorati_curvature") => Parameter::CasoratiCurvature,
            Some("contour_geodesic_torsion") => Parameter::ContourGeodesicTorsion,
            Some(o) => {
                return Err(ToolError::Validation(format!(
                    "'parameter' must be one of slope|aspect|mean_curvature|profile_curvature|tangential_curvature|plan_curvature|gaussian_curvature|casorati_curvature|contour_geodesic_torsion, got '{o}'"
                )))
            }
        },
    )
}

fn parse_fit(args: &ToolArgs) -> Result<Fit, ToolError> {
    match args.get("fit").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("quadratic") => Ok(Fit::Quadratic),
        Some("biquadratic") => Ok(Fit::Biquadratic),
        Some(o) => Err(ToolError::Validation(format!(
            "'fit' must be 'quadratic' or 'biquadratic', got '{o}'"
        ))),
    }
}

fn parse_neighborhood_type(args: &ToolArgs) -> Result<Neighborhood, ToolError> {
    match args
        .get("neighborhood_type")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some("fixed") => Ok(Neighborhood::Fixed),
        Some("adaptive") => Ok(Neighborhood::Adaptive),
        Some(o) => Err(ToolError::Validation(format!(
            "'neighborhood_type' must be 'fixed' or 'adaptive', got '{o}'"
        ))),
    }
}

fn parse_band(args: &ToolArgs) -> Result<u64, ToolError> {
    match parse_optional_f64(args, "band")? {
        None => Ok(1),
        Some(v) if v.fract() == 0.0 && v >= 1.0 => Ok(v as u64),
        Some(_) => Err(ToolError::Validation(
            "parameter 'band' must be a positive integer".to_string(),
        )),
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
    use wbraster::{Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn raster_path(rows: usize, cols: usize, cell: f64, vals: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: cell,
            cell_size_y: Some(cell),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: Default::default(),
            metadata: Default::default(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, vals[row * cols + col])
                    .unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> Raster {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SurfaceParametersTool.run(&args, &ctx()).unwrap();
        load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    /// A planar surface tilted 1 unit per cell east-west has slope atan(1)=45°.
    #[test]
    fn slope_of_tilted_plane() {
        // z increases by 1 per column, unit cells.
        let vals: Vec<f64> = (0..25).map(|i| (i % 5) as f64).collect();
        let r = run(json!({ "input": raster_path(5, 5, 1.0, &vals), "parameter": "slope" }));
        // Interior cell (2,2): slope should be 45°.
        assert!(
            (r.get(0, 2, 2) - 45.0).abs() < 1e-6,
            "got {}",
            r.get(0, 2, 2)
        );
    }

    /// East-rising plane: steepest descent faces west -> aspect ~270°.
    #[test]
    fn aspect_of_east_rising_plane() {
        let vals: Vec<f64> = (0..25).map(|i| (i % 5) as f64).collect();
        let r = run(json!({ "input": raster_path(5, 5, 1.0, &vals), "parameter": "aspect" }));
        assert!(
            (r.get(0, 2, 2) - 270.0).abs() < 1e-4,
            "got {}",
            r.get(0, 2, 2)
        );
    }

    /// A bowl z = x² + y² is convex: gaussian curvature > 0 at the center.
    #[test]
    fn gaussian_curvature_of_bowl_is_positive() {
        let mut vals = vec![0.0; 25];
        for r in 0..5 {
            for c in 0..5 {
                let x = c as f64 - 2.0;
                let y = 2.0 - r as f64;
                vals[r * 5 + c] = x * x + y * y;
            }
        }
        let out = run(json!({
            "input": raster_path(5, 5, 1.0, &vals),
            "parameter": "gaussian_curvature"
        }));
        assert!(out.get(0, 2, 2) > 0.0, "got {}", out.get(0, 2, 2));
    }

    /// A perfectly flat surface has zero slope on every fully-supported (interior)
    /// cell; edge cells with too few neighbors to fit the quadratic are no-data.
    #[test]
    fn flat_surface_zero_slope() {
        let r = run(json!({ "input": raster_path(5, 5, 1.0, &[7.0; 25]), "parameter": "slope" }));
        for row in 1..4 {
            for col in 1..4 {
                assert!(
                    r.get(0, row, col).abs() < 1e-9,
                    "cell ({row},{col}) = {}",
                    r.get(0, row, col)
                );
            }
        }
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SurfaceParametersTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "x.tif", "parameter": "bogus" })).is_err());
        assert!(bad(json!({ "input": "x.tif", "fit": "cubic" })).is_err());
        assert!(
            bad(json!({ "input": "x.tif", "parameter": "aspect", "fit": "biquadratic" })).is_ok()
        );
    }
}
