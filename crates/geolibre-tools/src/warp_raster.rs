//! GeoLibre tool: georeference a raster from ground control points (GCPs).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Warp* (Data Management). With no GDAL
//! in the stack there is currently no way to georeference scanned maps,
//! historical imagery, or drone frames. The bundled `thin_plate_spline` is a
//! point-interpolation tool, not an image warper, and the GeoLibre
//! `rubbersheet_features` is the vector twin of exactly this operation.
//!
//! GCPs pair a source pixel `(col,row)` with a target world coordinate
//! `(x,y)`. The tool fits an order-1/2/3 polynomial `world = P(col,row)` (to
//! size the output) and its inverse `(col,row) = Q(x,y)` (to resample), both by
//! least squares. Each output cell's world coordinate is mapped back to a source
//! pixel and sampled (`nearest` or `bilinear`). Per-GCP residuals and the total
//! RMS error (world units) are reported.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};

use crate::common::{load_input_raster, parse_optional_output};

pub struct WarpRasterTool;

impl Tool for WarpRasterTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "warp_raster",
            display_name: "Warp Raster",
            summary: "Georeference a raster from ground control points (like ArcGIS Warp): fit an order 1/2/3 polynomial from source-pixel → world-coordinate GCP pairs and resample into a new georeferenced grid (nearest/bilinear), reporting per-GCP residuals and RMS error. The raster half of the conflation story the vector rubbersheet_features already covers.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input raster (its pixel grid is the source coordinate space).",
                    required: true,
                },
                ToolParamSpec {
                    name: "gcps",
                    description: "Ground control points as 'col,row,x,y' quads separated by ';' (>= 3 for order 1, 6 for order 2, 10 for order 3).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output georeferenced raster. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "transform",
                    description: "Polynomial order: 'poly1' (affine; default), 'poly2', or 'poly3'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "resampling",
                    description: "'nearest' (default) or 'bilinear'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Output cell size in world units (default: input width mapped to the target extent).",
                    required: false,
                },
                ToolParamSpec {
                    name: "epsg",
                    description: "EPSG code to tag the output CRS (optional).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to warp (default 1).",
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
        let prm = parse_params(args)?;
        let need = prm.order.terms().len();
        if prm.gcps.len() < need {
            return Err(ToolError::Validation(format!(
                "{} GCPs given; order {} needs at least {}",
                prm.gcps.len(),
                prm.order.degree(),
                need
            )));
        }
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
        let output = parse_optional_output(args, "output")?;
        let prm = parse_params(args)?;

        let raster = load_input_raster(input)?;
        if prm.band < 0 || prm.band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "band {} out of range",
                prm.band + 1
            )));
        }
        let order = prm.order;

        // Fit forward P: (col,row) -> (x,y) for the output extent, and inverse
        // Q: (x,y) -> (col,row) for resampling.
        let fwd_x = fit(&prm.gcps, order, |g| (g.col, g.row, g.x))?;
        let fwd_y = fit(&prm.gcps, order, |g| (g.col, g.row, g.y))?;
        let inv_c = fit(&prm.gcps, order, |g| (g.x, g.y, g.col))?;
        let inv_r = fit(&prm.gcps, order, |g| (g.x, g.y, g.row))?;

        // RMS of forward residuals (world units).
        let mut sse = 0.0;
        let mut residuals = Vec::with_capacity(prm.gcps.len());
        for g in &prm.gcps {
            let px = eval(&fwd_x, g.col, g.row, order);
            let py = eval(&fwd_y, g.col, g.row, order);
            let r = (px - g.x).hypot(py - g.y);
            residuals.push(r);
            sse += r * r;
        }
        let rms = (sse / prm.gcps.len() as f64).sqrt();

        // Output extent from the forward-transformed image corners.
        let w = raster.cols as f64;
        let h = raster.rows as f64;
        let corners = [(0.0, 0.0), (w, 0.0), (0.0, h), (w, h)];
        let (mut xmin, mut ymin, mut xmax, mut ymax) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for (c, r) in corners {
            let x = eval(&fwd_x, c, r, order);
            let y = eval(&fwd_y, c, r, order);
            xmin = xmin.min(x);
            xmax = xmax.max(x);
            ymin = ymin.min(y);
            ymax = ymax.max(y);
        }
        let cell = prm.cell_size.unwrap_or(((xmax - xmin) / w).abs().max(1e-9));
        let out_cols = (((xmax - xmin) / cell).ceil() as usize).max(1);
        let out_rows = (((ymax - ymin) / cell).ceil() as usize).max(1);

        ctx.progress.info(&format!(
            "warp order {} -> {}x{} grid, RMS {:.4}",
            order.degree(),
            out_rows,
            out_cols,
            rms
        ));

        let nodata = raster.nodata;
        let mut out = Raster::new(RasterConfig {
            cols: out_cols,
            rows: out_rows,
            bands: 1,
            x_min: xmin,
            y_min: ymin,
            cell_size: cell,
            cell_size_y: Some(cell),
            nodata,
            data_type: DataType::F32,
            crs: match prm.epsg {
                Some(e) => CrsInfo {
                    epsg: Some(e),
                    wkt: None,
                    proj4: None,
                },
                None => raster.crs.clone(),
            },
            metadata: Vec::new(),
        });

        let src_cols = raster.cols as isize;
        let src_rows = raster.rows as isize;
        for or in 0..out_rows {
            for oc in 0..out_cols {
                // World coordinate of this output cell centre.
                let wx = xmin + (oc as f64 + 0.5) * cell;
                let wy = ymax - (or as f64 + 0.5) * cell;
                // Inverse map to source pixel.
                let sc = eval(&inv_c, wx, wy, order);
                let sr = eval(&inv_r, wx, wy, order);
                let v = sample(
                    &raster,
                    prm.band,
                    sc,
                    sr,
                    src_cols,
                    src_rows,
                    prm.bilinear,
                    nodata,
                );
                out.set(0, or as isize, oc as isize, v)
                    .map_err(|e| ToolError::Execution(format!("write failed: {e}")))?;
            }
            ctx.progress.progress((or as f64 + 1.0) / out_rows as f64);
        }

        let out_path = crate::common::write_or_store_output(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("rms_error".to_string(), json!(rms));
        outputs.insert("residuals".to_string(), json!(residuals));
        outputs.insert("out_cols".to_string(), json!(out_cols));
        outputs.insert("out_rows".to_string(), json!(out_rows));
        Ok(ToolRunResult { outputs })
    }
}

/// Samples the source raster at fractional pixel `(col,row)`.
#[allow(clippy::too_many_arguments)]
fn sample(
    raster: &Raster,
    band: isize,
    c: f64,
    r: f64,
    cols: isize,
    rows: isize,
    bilinear: bool,
    nodata: f64,
) -> f64 {
    if bilinear {
        let c0 = c.floor();
        let r0 = r.floor();
        let fc = c - c0;
        let fr = r - r0;
        let (c0, r0) = (c0 as isize, r0 as isize);
        let g = |rr: isize, cc: isize| -> Option<f64> {
            if rr < 0 || cc < 0 || rr >= rows || cc >= cols {
                return None;
            }
            let v = raster.get(band, rr, cc);
            (v != nodata && v.is_finite()).then_some(v)
        };
        match (g(r0, c0), g(r0, c0 + 1), g(r0 + 1, c0), g(r0 + 1, c0 + 1)) {
            (Some(v00), Some(v01), Some(v10), Some(v11)) => {
                let top = v00 * (1.0 - fc) + v01 * fc;
                let bot = v10 * (1.0 - fc) + v11 * fc;
                top * (1.0 - fr) + bot * fr
            }
            _ => nearest(raster, band, c, r, cols, rows, nodata),
        }
    } else {
        nearest(raster, band, c, r, cols, rows, nodata)
    }
}

fn nearest(
    raster: &Raster,
    band: isize,
    c: f64,
    r: f64,
    cols: isize,
    rows: isize,
    nodata: f64,
) -> f64 {
    let cc = c.round() as isize;
    let rr = r.round() as isize;
    if rr < 0 || cc < 0 || rr >= rows || cc >= cols {
        nodata
    } else {
        raster.get(band, rr, cc)
    }
}

// ── Polynomial least squares ────────────────────────────────────────────────

/// Fits `out = P(inx, iny)` of the given order via normal equations. `select`
/// pulls (inx, iny, out) from each GCP.
fn fit(
    gcps: &[Gcp],
    order: Order,
    select: impl Fn(&Gcp) -> (f64, f64, f64),
) -> Result<Vec<f64>, ToolError> {
    let terms = order.terms();
    let m = terms.len();
    // Normal matrix (m x m) and rhs (m).
    let mut ata = vec![vec![0.0f64; m]; m];
    let mut atb = vec![0.0f64; m];
    for g in gcps {
        let (ix, iy, out) = select(g);
        let d = design(ix, iy, &terms);
        for i in 0..m {
            atb[i] += d[i] * out;
            for j in 0..m {
                ata[i][j] += d[i] * d[j];
            }
        }
    }
    solve(&mut ata, &mut atb).ok_or_else(|| {
        ToolError::Execution("GCP configuration is degenerate (singular fit)".into())
    })
}

fn design(x: f64, y: f64, terms: &[(u32, u32)]) -> Vec<f64> {
    terms
        .iter()
        .map(|&(px, py)| x.powi(px as i32) * y.powi(py as i32))
        .collect()
}

fn eval(coeffs: &[f64], x: f64, y: f64, order: Order) -> f64 {
    let terms = order.terms();
    design(x, y, &terms)
        .iter()
        .zip(coeffs)
        .map(|(d, c)| d * c)
        .sum()
}

/// Gaussian elimination with partial pivoting; solves `a x = b` in place.
fn solve(a: &mut [Vec<f64>], b: &mut [f64]) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        // Pivot.
        let mut piv = col;
        for r in (col + 1)..n {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-12 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        // Eliminate.
        for r in (col + 1)..n {
            let f = a[r][col] / a[col][col];
            for c in col..n {
                a[r][c] -= f * a[col][c];
            }
            b[r] -= f * b[col];
        }
    }
    // Back-substitute.
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let mut s = b[r];
        for c in (r + 1)..n {
            s -= a[r][c] * x[c];
        }
        x[r] = s / a[r][r];
    }
    Some(x)
}

// ── Parameters ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Order {
    Poly1,
    Poly2,
    Poly3,
}

impl Order {
    fn degree(&self) -> u32 {
        match self {
            Order::Poly1 => 1,
            Order::Poly2 => 2,
            Order::Poly3 => 3,
        }
    }
    /// Monomial exponent pairs (px, py) with px + py <= degree.
    fn terms(&self) -> Vec<(u32, u32)> {
        let d = self.degree();
        let mut t = Vec::new();
        for total in 0..=d {
            for px in 0..=total {
                t.push((px, total - px));
            }
        }
        t
    }
}

struct Gcp {
    col: f64,
    row: f64,
    x: f64,
    y: f64,
}

struct Params {
    gcps: Vec<Gcp>,
    order: Order,
    bilinear: bool,
    cell_size: Option<f64>,
    epsg: Option<u32>,
    band: isize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let gcp_str = args
        .get("gcps")
        .and_then(Value::as_str)
        .ok_or_else(|| ToolError::Validation("missing required parameter 'gcps'".to_string()))?;
    let gcps: Vec<Gcp> = gcp_str
        .split(';')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|quad| {
            let nums: Vec<f64> = quad
                .split(',')
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(|n| {
                    n.parse::<f64>()
                        .map_err(|_| format!("bad GCP number '{n}'"))
                })
                .collect::<Result<_, _>>()
                .map_err(ToolError::Validation)?;
            if nums.len() != 4 {
                return Err(ToolError::Validation(format!(
                    "each GCP needs 'col,row,x,y' (4 numbers), got {}",
                    nums.len()
                )));
            }
            Ok(Gcp {
                col: nums[0],
                row: nums[1],
                x: nums[2],
                y: nums[3],
            })
        })
        .collect::<Result<_, _>>()?;
    if gcps.len() < 3 {
        return Err(ToolError::Validation("need at least 3 GCPs".to_string()));
    }
    let order = match args.get("transform").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("poly1") => Order::Poly1,
        Some("poly2") => Order::Poly2,
        Some("poly3") => Order::Poly3,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'transform' must be poly1/poly2/poly3, got '{o}'"
            )))
        }
    };
    let bilinear = match args
        .get("resampling")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some("nearest") => false,
        Some("bilinear") => true,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'resampling' must be nearest/bilinear, got '{o}'"
            )))
        }
    };
    let cell_size = match args.get("cell_size") {
        None | Some(Value::Null) => None,
        Some(Value::Number(n)) => n.as_f64().filter(|v| *v > 0.0),
        Some(Value::String(s)) if s.trim().is_empty() => None,
        Some(Value::String(s)) => Some(
            s.trim()
                .parse::<f64>()
                .map_err(|_| ToolError::Validation("'cell_size' must be a number".into()))?,
        ),
        _ => None,
    };
    let epsg = match args.get("epsg") {
        None | Some(Value::Null) => None,
        Some(Value::Number(n)) => n.as_u64().map(|v| v as u32),
        Some(Value::String(s)) if s.trim().is_empty() => None,
        Some(Value::String(s)) => Some(
            s.trim()
                .parse::<u32>()
                .map_err(|_| ToolError::Validation("'epsg' must be an integer".into()))?,
        ),
        _ => None,
    };
    let band_1based = args.get("band").and_then(Value::as_u64).unwrap_or(1).max(1);
    Ok(Params {
        gcps,
        order,
        bilinear,
        cell_size,
        epsg,
        band: (band_1based - 1) as isize,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn ramp_raster(cols: usize, rows: usize) -> String {
        // value = col*1 + row*10, so we can check resampling picks the right pixel.
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: None,
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(
                    0,
                    row as isize,
                    col as isize,
                    col as f64 + 10.0 * row as f64,
                )
                .unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    /// An affine GCP set (world = pixel scaled by 2, offset by 100/200) fits with
    /// ~0 RMS and produces a georeferenced grid at the right extent.
    #[test]
    fn affine_warp_fits_exactly() {
        let input = ramp_raster(10, 10);
        // world_x = 100 + 2*col ; world_y = 200 + 2*row
        let gcps = "0,0,100,200; 9,0,118,200; 0,9,100,218; 9,9,118,218";
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input, "gcps": gcps, "transform": "poly1", "cell_size": 2.0,
        }))
        .unwrap();
        let out = WarpRasterTool.run(&args, &ctx()).unwrap();
        assert!(
            out.outputs["rms_error"].as_f64().unwrap() < 1e-6,
            "affine GCPs should fit exactly"
        );
        let r = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        assert!(r.cols > 0 && r.rows > 0);
        assert!(
            (r.x_min - 100.0).abs() < 1.0,
            "output x origin near world min"
        );
    }

    /// A rotation-free scaling recovers source values at the mapped locations.
    #[test]
    fn warp_preserves_values() {
        let input = ramp_raster(8, 8);
        // Identity-ish: world == pixel (offset 0, scale 1).
        let gcps = "0,0,0,0; 7,0,7,0; 0,7,0,7; 7,7,7,7";
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input, "gcps": gcps, "cell_size": 1.0, "resampling": "nearest",
        }))
        .unwrap();
        let out = WarpRasterTool.run(&args, &ctx()).unwrap();
        let r = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        // Resampling must pull real source values (ramp spans 0..=7+70=77), so
        // the warped grid should reproduce most of that range.
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for row in 0..r.rows {
            for col in 0..r.cols {
                let v = r.get(0, row as isize, col as isize);
                if v != r.nodata {
                    lo = lo.min(v);
                    hi = hi.max(v);
                }
            }
        }
        assert!(
            lo <= 15.0,
            "min warped value should be near the source min, got {lo}"
        );
        assert!(
            hi >= 70.0,
            "max warped value should be near the source max, got {hi}"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            WarpRasterTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.tif", "gcps": "0,0,1,1" })).is_err()); // <3 gcps
        assert!(bad(
            json!({ "input": "a.tif", "gcps": "0,0,1,1;1,0,2,1;0,1,1,2", "transform": "poly2" })
        )
        .is_err()); // needs 6
        assert!(bad(json!({ "input": "a.tif", "gcps": "0,0,1,1;1,0,2,1;0,1,1,2" })).is_ok());
    }
}
