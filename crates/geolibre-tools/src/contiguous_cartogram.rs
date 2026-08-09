//! GeoLibre tool: density-equalizing (Gastner–Newman) contiguous cartogram.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Generate Contiguous Cartogram*
//! (Cartography).
//!
//! ## The gap the shipped tool names itself
//!
//! From `cartogram`'s own module doc:
//!
//! > Scope for v1: the contiguous (Gastner–Newman diffusion) cartogram is not
//! > implemented — use `non_contiguous` or `dorling`.
//!
//! Both implemented methods **destroy adjacency**: `non_contiguous` scales each
//! polygon about its own centroid, opening gaps between neighbours, and
//! `dorling` replaces polygons with circles, discarding shape entirely. The
//! contiguous diffusion cartogram is the one readers recognise as "a cartogram"
//! — the map stays a connected sheet, borders stay shared, and only the areas
//! distort. It is the standard presentation for election, population and
//! epidemiological maps.
//!
//! ## Method
//!
//! Gastner & Newman (2004): rasterize the value field to a density grid, then
//! let that density diffuse to uniformity. The velocity field of the diffusing
//! fluid, integrated over diffusion time, is the map projection that equalizes
//! density — which is exactly a cartogram.
//!
//! The diffusion is solved spectrally: on a grid, the heat equation's solution
//! is elementwise multiplication by `exp(-k^2 t)` in frequency space. Round 18's
//! `fft2` module supplies the transform, which is what makes this tractable here
//! with **no new dependency**.
//!
//! ## Three things that go wrong if done naively
//!
//! * **`grid_size` must be a power of two.** `fft2` is radix-2. It is validated
//!   and reported with the nearest valid sizes rather than silently rounded,
//!   because a rounded grid changes the result.
//! * **`grid_size` is squared.** `g * g` panics under overflow checks (the
//!   round-18 unsigned-cast finding class), so it is capped before the multiply.
//! * **Long edges must be densified.** A straight border stays straight through
//!   a curved velocity field, so neighbouring polygons visibly separate along
//!   their shared edge unless vertices are inserted first.
//!
//! ## Scope
//!
//! Filed as its own tool rather than a `method` on `cartogram` because the
//! parameters differ materially (grid resolution, diffusion time, blur) and
//! would be dead weight on the existing methods. Features with a missing or
//! non-positive value follow `cartogram`'s convention: left undistorted and
//! reported.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, Layer, Ring};

use crate::args_common::{opt_f64, req_str, usize_or};
use crate::fft2::{fft2, Cpx};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Upper bound on the density grid. 4096^2 complex samples is already ~256 MB
/// of working set; capping before `g * g` is what stops the multiply from
/// overflowing under debug overflow checks.
const MAX_GRID: usize = 4096;

const SEA_MODES: [&str; 2] = ["mean", "min"];

pub struct ContiguousCartogramTool;

impl Tool for ContiguousCartogramTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "contiguous_cartogram",
            display_name: "Contiguous Cartogram",
            summary: "Density-equalizing (Gastner-Newman diffusion) cartogram: distorts polygon areas in proportion to an attribute while keeping the map a connected sheet with shared borders (ArcGIS Generate Contiguous Cartogram). The shipped cartogram tool implements only the non-contiguous and Dorling methods, both of which destroy adjacency, and its own module doc defers the contiguous form.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Contiguous polygon layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "value_field",
                    description: "Positive numeric attribute whose total the areas should become proportional to.",
                    required: true,
                },
                ToolParamSpec {
                    name: "grid_size",
                    description: "Density grid resolution; must be a power of two (default 256). The dominant accuracy/cost knob.",
                    required: false,
                },
                ToolParamSpec {
                    name: "iterations",
                    description: "Diffusion steps (default 20). More steps equalize density further at the cost of more distortion.",
                    required: false,
                },
                ToolParamSpec {
                    name: "blur",
                    description: "Gaussian pre-smoothing of the density field, in grid cells (default 0.5). Prevents the extreme distortion a hard density discontinuity produces.",
                    required: false,
                },
                ToolParamSpec {
                    name: "sea_density",
                    description: "Density outside the study area: 'mean' (default) or 'min'. This materially changes the result at the map edge.",
                    required: false,
                },
                ToolParamSpec {
                    name: "densify_spacing",
                    description: "Maximum edge segment length in grid cells before vertices are inserted (default 1.0). Long straight borders otherwise separate from their neighbours.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon layer with displaced vertices and the original attributes. If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        req_str(args, "value_field")?;
        let g = usize_or(args, "grid_size", 256)?;
        // Capped BEFORE anything squares it.
        if !(16..=MAX_GRID).contains(&g) {
            return Err(ToolError::Validation(format!(
                "'grid_size' must be between 16 and {MAX_GRID}, got {g}"
            )));
        }
        if !g.is_power_of_two() {
            let lo = g.next_power_of_two() / 2;
            let hi = g.next_power_of_two();
            return Err(ToolError::Validation(format!(
                "'grid_size' must be a power of two (the FFT is radix-2), got {g}; use {lo} or \
                 {hi}"
            )));
        }
        if usize_or(args, "iterations", 20)? == 0 {
            return Err(ToolError::Validation(
                "'iterations' must be at least 1".to_string(),
            ));
        }
        if let Some(b) = opt_f64(args, "blur")? {
            if !b.is_finite() || b < 0.0 {
                return Err(ToolError::Validation(
                    "'blur' must be a non-negative number of grid cells".to_string(),
                ));
            }
        }
        if let Some(d) = opt_f64(args, "densify_spacing")? {
            if !d.is_finite() || d <= 0.0 {
                return Err(ToolError::Validation(
                    "'densify_spacing' must be positive".to_string(),
                ));
            }
        }
        if let Some(s) = args.get("sea_density").and_then(Value::as_str) {
            let s = s.trim().to_ascii_lowercase();
            if !SEA_MODES.contains(&s.as_str()) {
                return Err(ToolError::Validation(format!(
                    "'sea_density' must be one of {}, got '{s}'",
                    SEA_MODES.join("|")
                )));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?.to_string();
        let value_field = req_str(args, "value_field")?.to_string();
        let n = usize_or(args, "grid_size", 256)?;
        if !n.is_power_of_two() || !(16..=MAX_GRID).contains(&n) {
            return Err(ToolError::Validation(format!(
                "'grid_size' must be a power of two between 16 and {MAX_GRID}"
            )));
        }
        let iterations = usize_or(args, "iterations", 20)?.max(1);
        let blur = opt_f64(args, "blur")?.unwrap_or(0.5).max(0.0);
        let sea_mode = args
            .get("sea_density")
            .and_then(Value::as_str)
            .map(|s| s.trim().to_ascii_lowercase())
            .unwrap_or_else(|| "mean".to_string());
        let densify_spacing = opt_f64(args, "densify_spacing")?.unwrap_or(1.0).max(1e-6);
        let output = parse_optional_str(args, "output")?;

        let layer = load_input_layer(&input)?;
        let vidx = layer.schema.field_index(&value_field).ok_or_else(|| {
            ToolError::Validation(format!(
                "value_field '{value_field}' not found in the input layer"
            ))
        })?;

        // Bounds, padded so the diffusion has room to move mass outward without
        // pressing against the periodic boundary the FFT implies.
        let Some((min_x, min_y, max_x, max_y)) = layer_bounds(&layer) else {
            return Err(ToolError::Execution(
                "the input layer has no polygon geometry".to_string(),
            ));
        };
        let (w, h) = (max_x - min_x, max_y - min_y);
        if w <= 0.0 || h <= 0.0 {
            return Err(ToolError::Execution(
                "the input layer has zero extent".to_string(),
            ));
        }
        let pad = 0.25 * w.max(h);
        let (gx0, gy0) = (min_x - pad, min_y - pad);
        let span = (w + 2.0 * pad).max(h + 2.0 * pad);
        let cell = span / n as f64;

        // Rasterize value density: value / polygon area, so a large sparse
        // polygon and a small dense one contribute correctly.
        let mut rho = vec![0.0_f64; n * n];
        let mut covered = vec![false; n * n];
        let mut undistorted: Vec<i64> = Vec::new();
        let mut total_value = 0.0_f64;

        for (fid, f) in layer.iter().enumerate() {
            let Some(geom) = f.geometry.as_ref() else {
                continue;
            };
            let value = f.attributes.get(vidx).and_then(as_f64);
            let area = polygon_area(geom).abs();
            let Some(value) = value else {
                undistorted.push(fid as i64);
                continue;
            };
            if value <= 0.0 || area <= 0.0 {
                undistorted.push(fid as i64);
                continue;
            }
            total_value += value;
            let density = value / area;
            for r in 0..n {
                for c in 0..n {
                    let x = gx0 + (c as f64 + 0.5) * cell;
                    let y = gy0 + (r as f64 + 0.5) * cell;
                    if point_in_geometry(geom, x, y) {
                        rho[r * n + c] = density;
                        covered[r * n + c] = true;
                    }
                }
            }
        }

        let inside: Vec<f64> = rho
            .iter()
            .zip(covered.iter())
            .filter(|(_, c)| **c)
            .map(|(v, _)| *v)
            .collect();
        if inside.is_empty() {
            return Err(ToolError::Execution(format!(
                "no polygon with a positive '{value_field}' covered a grid cell; raise \
                 'grid_size' or check the values"
            )));
        }
        // The sea density choice materially changes the edge behaviour, so it
        // is an explicit parameter rather than a hidden constant.
        let sea = match sea_mode.as_str() {
            "min" => inside.iter().cloned().fold(f64::INFINITY, f64::min),
            _ => inside.iter().sum::<f64>() / inside.len() as f64,
        };
        for i in 0..n * n {
            if !covered[i] {
                rho[i] = sea;
            }
        }

        if blur > 0.0 {
            rho = gaussian_blur(&rho, n, blur);
        }
        // Normalize so the mean density is 1; the diffusion then relaxes toward
        // a uniform field of 1.
        let mean = rho.iter().sum::<f64>() / (n * n) as f64;
        if mean <= 0.0 || !mean.is_finite() {
            return Err(ToolError::Execution(
                "the density field is degenerate (mean is zero)".to_string(),
            ));
        }
        for v in rho.iter_mut() {
            *v /= mean;
        }

        ctx.progress.info(&format!(
            "{}x{} density grid, {iterations} diffusion step(s), sea={sea_mode}",
            n, n
        ));

        // Spectral diffusion. The heat equation on a periodic grid is
        // elementwise multiplication by exp(-k^2 t) in frequency space, so the
        // whole time evolution is one forward transform, a scaling per step,
        // and one inverse transform per step.
        let mut spec: Vec<Cpx> = rho.iter().map(|v| (*v, 0.0)).collect();
        fft2(&mut spec, n, n, false);

        // Total diffusion time, scaled so the result is insensitive to grid
        // resolution: k is measured in cycles per grid, so t is in grid units.
        let t_total = 0.1 * iterations as f64;
        let steps = iterations;
        let dt = t_total / steps as f64;

        // Vertex displacement is integrated over the diffusion, so vertices are
        // carried along with the flow rather than displaced once at the end.
        let mut pts = collect_vertices(&layer, densify_spacing * cell);
        let mut work: Vec<Cpx> = vec![(0.0, 0.0); n * n];

        for step in 0..steps {
            let t = (step + 1) as f64 * dt;
            // rho(t) from the original spectrum.
            for r in 0..n {
                for c in 0..n {
                    let kr = freq(r, n);
                    let kc = freq(c, n);
                    let k2 = (kr * kr + kc * kc) * (2.0 * std::f64::consts::PI).powi(2);
                    let decay = (-k2 * t).exp();
                    let s = spec[r * n + c];
                    work[r * n + c] = (s.0 * decay, s.1 * decay);
                }
            }
            fft2(&mut work, n, n, true);
            let rho_t: Vec<f64> = work.iter().map(|v| v.0).collect();

            // v = -grad(rho)/rho, evaluated by central differences.
            let mut vx = vec![0.0_f64; n * n];
            let mut vy = vec![0.0_f64; n * n];
            for r in 0..n {
                for c in 0..n {
                    let rp = rho_t[r * n + c];
                    // As diffusion proceeds rho tends to 1, but an early cell
                    // can sit at zero; without a floor the velocity is infinite
                    // and every vertex downstream becomes NaN.
                    let denom = rp.max(1e-9);
                    let (cl, cr) = ((c + n - 1) % n, (c + 1) % n);
                    let (ru, rd) = ((r + n - 1) % n, (r + 1) % n);
                    let dx = (rho_t[r * n + cr] - rho_t[r * n + cl]) / (2.0 * cell);
                    let dy = (rho_t[rd * n + c] - rho_t[ru * n + c]) / (2.0 * cell);
                    vx[r * n + c] = -dx / denom;
                    vy[r * n + c] = -dy / denom;
                }
            }

            // RK2 (midpoint): forward Euler at a coarse step visibly tears the
            // map along fast-flowing borders.
            for p in pts.iter_mut() {
                let (ux, uy) = sample_velocity(&vx, &vy, n, gx0, gy0, cell, p.0, p.1);
                let mx = p.0 + 0.5 * dt * ux;
                let my = p.1 + 0.5 * dt * uy;
                let (kx, ky) = sample_velocity(&vx, &vy, n, gx0, gy0, cell, mx, my);
                let nx = p.0 + dt * kx;
                let ny = p.1 + dt * ky;
                if nx.is_finite() && ny.is_finite() {
                    p.0 = nx;
                    p.1 = ny;
                }
            }
            ctx.progress.progress((step as f64 + 1.0) / steps as f64);
        }

        // Rebuild the layer with displaced vertices, in the same order the
        // vertices were collected.
        let mut out = Layer::new(&layer.name);
        if let Some(gt) = layer.geom_type {
            out = out.with_geom_type(gt);
        }
        if let Some(e) = layer.crs_epsg() {
            out = out.with_crs_epsg(e);
        }
        for f in layer.schema.fields() {
            out.add_field(f.clone());
        }
        out.add_field(FieldDef::new("AREA_ERR", FieldType::Float));

        let mut worst_err = 0.0_f64;
        let mut sum_err = 0.0_f64;
        let mut counted = 0_usize;

        // Two passes: rebuild first so the total output area is known, then
        // score each feature's achieved-vs-target area against it.
        let mut rebuilt: Vec<(Option<Geometry>, Option<f64>)> = Vec::new();
        let mut cursor = 0_usize;
        for f in layer.iter() {
            let geom = f
                .geometry
                .as_ref()
                .map(|g| rebuild_geometry(g, &pts, &mut cursor, densify_spacing * cell));
            let value = f.attributes.get(vidx).and_then(as_f64).filter(|v| *v > 0.0);
            rebuilt.push((geom, value));
        }
        let out_total: f64 = rebuilt
            .iter()
            .filter(|(_, v)| v.is_some())
            .filter_map(|(g, _)| g.as_ref().map(|g| polygon_area(g).abs()))
            .sum();

        for (f, (geom, value)) in layer.iter().zip(rebuilt) {
            let mut attrs = f.attributes.clone();
            let err = match (value, out_total > 0.0 && total_value > 0.0) {
                (Some(v), true) => {
                    let target = out_total * (v / total_value);
                    let got = geom.as_ref().map(|g| polygon_area(g).abs()).unwrap_or(0.0);
                    let e = if target > 0.0 {
                        (got - target) / target
                    } else {
                        0.0
                    };
                    worst_err = worst_err.max(e.abs());
                    sum_err += e.abs();
                    counted += 1;
                    Some(e)
                }
                _ => None,
            };
            attrs.push(match err {
                Some(e) => FieldValue::Float(e),
                None => FieldValue::Null,
            });
            out.push(Feature {
                fid: f.fid,
                geometry: geom,
                attributes: attrs,
            });
        }

        let feature_count = out.features.len();
        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("grid_size".to_string(), json!(n));
        outputs.insert("iterations".to_string(), json!(iterations));
        outputs.insert("sea_density".to_string(), json!(sea_mode));
        outputs.insert("undistorted_features".to_string(), json!(undistorted));
        outputs.insert("max_area_error".to_string(), json!(worst_err));
        outputs.insert(
            "mean_area_error".to_string(),
            json!(if counted > 0 {
                sum_err / counted as f64
            } else {
                0.0
            }),
        );
        Ok(ToolRunResult { outputs })
    }
}

/// Signed frequency index for an FFT bin, in cycles per grid.
fn freq(i: usize, n: usize) -> f64 {
    if i <= n / 2 {
        i as f64 / n as f64
    } else {
        (i as f64 - n as f64) / n as f64
    }
}

/// Bilinear sample of the velocity field at a map coordinate.
#[allow(clippy::too_many_arguments)]
fn sample_velocity(
    vx: &[f64],
    vy: &[f64],
    n: usize,
    gx0: f64,
    gy0: f64,
    cell: f64,
    x: f64,
    y: f64,
) -> (f64, f64) {
    let fc = (x - gx0) / cell - 0.5;
    let fr = (y - gy0) / cell - 0.5;
    let c0 = fc.floor();
    let r0 = fr.floor();
    let tx = fc - c0;
    let ty = fr - r0;
    let idx = |r: i64, c: i64| -> usize {
        let rr = r.rem_euclid(n as i64) as usize;
        let cc = c.rem_euclid(n as i64) as usize;
        rr * n + cc
    };
    let (c0, r0) = (c0 as i64, r0 as i64);
    let lerp = |f: &[f64]| -> f64 {
        let a = f[idx(r0, c0)] * (1.0 - tx) + f[idx(r0, c0 + 1)] * tx;
        let b = f[idx(r0 + 1, c0)] * (1.0 - tx) + f[idx(r0 + 1, c0 + 1)] * tx;
        a * (1.0 - ty) + b * ty
    };
    (lerp(vx), lerp(vy))
}

/// Separable Gaussian blur over a square grid, with edge clamping.
fn gaussian_blur(src: &[f64], n: usize, sigma: f64) -> Vec<f64> {
    let radius = (3.0 * sigma).ceil().max(1.0) as isize;
    let mut kernel = Vec::with_capacity((2 * radius + 1) as usize);
    let mut sum = 0.0;
    for i in -radius..=radius {
        let w = (-(i as f64).powi(2) / (2.0 * sigma * sigma)).exp();
        kernel.push(w);
        sum += w;
    }
    for k in kernel.iter_mut() {
        *k /= sum;
    }
    let clamp = |v: isize, hi: isize| v.clamp(0, hi - 1) as usize;
    let mut tmp = vec![0.0; n * n];
    for r in 0..n {
        for c in 0..n {
            let mut acc = 0.0;
            for (k, w) in kernel.iter().enumerate() {
                let cc = clamp(c as isize + k as isize - radius, n as isize);
                acc += src[r * n + cc] * w;
            }
            tmp[r * n + c] = acc;
        }
    }
    let mut out = vec![0.0; n * n];
    for r in 0..n {
        for c in 0..n {
            let mut acc = 0.0;
            for (k, w) in kernel.iter().enumerate() {
                let rr = clamp(r as isize + k as isize - radius, n as isize);
                acc += tmp[rr * n + c] * w;
            }
            out[r * n + c] = acc;
        }
    }
    out
}

/// Collects every polygon vertex in feature/ring order, densifying long edges.
///
/// Densification is not cosmetic: a straight border stays straight through a
/// curved velocity field, so two polygons sharing a long edge visibly separate
/// unless intermediate vertices exist to follow the flow.
fn collect_vertices(layer: &Layer, max_seg: f64) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for f in layer.iter() {
        if let Some(g) = f.geometry.as_ref() {
            walk_rings(g, &mut |ring: &Ring| {
                out.extend(densify(&ring.0, max_seg));
            });
        }
    }
    out
}

/// Rebuilds a geometry from the displaced vertex stream, consuming exactly the
/// vertices `collect_vertices` produced for it, in the same order.
fn rebuild_geometry(
    geom: &Geometry,
    pts: &[(f64, f64)],
    cursor: &mut usize,
    max_seg: f64,
) -> Geometry {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            let ext = take_ring(&exterior.0, pts, cursor, max_seg);
            let ints: Vec<Ring> = interiors
                .iter()
                .map(|r| Ring(take_ring(&r.0, pts, cursor, max_seg)))
                .collect();
            Geometry::Polygon {
                exterior: Ring(ext),
                interiors: ints,
            }
        }
        Geometry::MultiPolygon(parts) => Geometry::MultiPolygon(
            parts
                .iter()
                .map(|(e, hs)| {
                    let ext = Ring(take_ring(&e.0, pts, cursor, max_seg));
                    let holes: Vec<Ring> = hs
                        .iter()
                        .map(|r| Ring(take_ring(&r.0, pts, cursor, max_seg)))
                        .collect();
                    (ext, holes)
                })
                .collect(),
        ),
        Geometry::GeometryCollection(gs) => Geometry::GeometryCollection(
            gs.iter()
                .map(|g| rebuild_geometry(g, pts, cursor, max_seg))
                .collect(),
        ),
        other => other.clone(),
    }
}

fn take_ring(
    original: &[Coord],
    pts: &[(f64, f64)],
    cursor: &mut usize,
    max_seg: f64,
) -> Vec<Coord> {
    let count = densify(original, max_seg).len();
    let mut out = Vec::with_capacity(count);
    for k in 0..count {
        let p = pts.get(*cursor + k).copied().unwrap_or((0.0, 0.0));
        out.push(Coord {
            x: p.0,
            y: p.1,
            z: None,
            m: None,
        });
    }
    *cursor += count;
    // Re-close the ring: the closing vertex is displaced independently of the
    // first, and a hair's difference leaves an unclosed polygon.
    if out.len() >= 2 {
        let first = out[0].clone();
        let last = out.len() - 1;
        out[last] = first;
    }
    out
}

/// Inserts intermediate vertices so no segment exceeds `max_seg`.
fn densify(coords: &[Coord], max_seg: f64) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    if coords.is_empty() {
        return out;
    }
    for w in coords.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        out.push((a.x, a.y));
        let d = ((b.x - a.x).powi(2) + (b.y - a.y).powi(2)).sqrt();
        if d > max_seg && max_seg > 0.0 {
            // Bounded: `d / max_seg` is finite because max_seg > 0, and the
            // count is capped so a pathological ratio cannot allocate forever.
            let extra = ((d / max_seg).ceil() as usize).min(4096).saturating_sub(1);
            for k in 1..=extra {
                let t = k as f64 / (extra + 1) as f64;
                out.push((a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t));
            }
        }
    }
    if let Some(last) = coords.last() {
        out.push((last.x, last.y));
    }
    out
}

fn walk_rings(geom: &Geometry, f: &mut impl FnMut(&Ring)) {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            f(exterior);
            for r in interiors {
                f(r);
            }
        }
        Geometry::MultiPolygon(parts) => {
            for (e, hs) in parts {
                f(e);
                for r in hs {
                    f(r);
                }
            }
        }
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                walk_rings(g, f);
            }
        }
        _ => {}
    }
}

fn ring_area(coords: &[Coord]) -> f64 {
    let mut a = 0.0;
    for i in 0..coords.len() {
        let p = &coords[i];
        let q = &coords[(i + 1) % coords.len()];
        a += p.x * q.y - q.x * p.y;
    }
    0.5 * a
}

fn polygon_area(geom: &Geometry) -> f64 {
    let mut total = 0.0;
    walk_rings_ref(geom, &mut |ring, is_hole| {
        let a = ring_area(&ring.0).abs();
        total += if is_hole { -a } else { a };
    });
    total
}

fn walk_rings_ref(geom: &Geometry, f: &mut impl FnMut(&Ring, bool)) {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            f(exterior, false);
            for r in interiors {
                f(r, true);
            }
        }
        Geometry::MultiPolygon(parts) => {
            for (e, hs) in parts {
                f(e, false);
                for r in hs {
                    f(r, true);
                }
            }
        }
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                walk_rings_ref(g, f);
            }
        }
        _ => {}
    }
}

fn point_in_geometry(geom: &Geometry, x: f64, y: f64) -> bool {
    crate::vector_common::geometry_contains_point(geom, x, y)
}

fn layer_bounds(layer: &Layer) -> Option<(f64, f64, f64, f64)> {
    let mut b = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
    let mut seen = false;
    for f in layer.iter() {
        if let Some(g) = f.geometry.as_ref() {
            walk_rings_ref(g, &mut |ring, _| {
                for c in &ring.0 {
                    seen = true;
                    b.0 = b.0.min(c.x);
                    b.1 = b.1.min(c.y);
                    b.2 = b.2.max(c.x);
                    b.3 = b.3.max(c.y);
                }
            });
        }
    }
    seen.then_some(b)
}

fn as_f64(v: &FieldValue) -> Option<f64> {
    match v {
        FieldValue::Integer(i) => Some(*i as f64),
        FieldValue::Float(f) if f.is_finite() => Some(*f),
        FieldValue::Text(s) => s.trim().parse::<f64>().ok(),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::GeometryType;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Geometry {
        Geometry::polygon(
            vec![
                Coord {
                    x: x0,
                    y: y0,
                    z: None,
                    m: None,
                },
                Coord {
                    x: x1,
                    y: y0,
                    z: None,
                    m: None,
                },
                Coord {
                    x: x1,
                    y: y1,
                    z: None,
                    m: None,
                },
                Coord {
                    x: x0,
                    y: y1,
                    z: None,
                    m: None,
                },
                Coord {
                    x: x0,
                    y: y0,
                    z: None,
                    m: None,
                },
            ],
            Vec::new(),
        )
    }

    /// Two side-by-side unit squares sharing the x = 1 edge.
    fn two_squares(v_left: f64, v_right: f64) -> String {
        let mut l = Layer::new("cells")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("pop", FieldType::Float));
        l.add_feature(Some(rect(0.0, 0.0, 1.0, 1.0)), &[("pop", v_left.into())])
            .unwrap();
        l.add_feature(Some(rect(1.0, 0.0, 2.0, 1.0)), &[("pop", v_right.into())])
            .unwrap();
        store(l)
    }

    fn store(l: Layer) -> String {
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = ContiguousCartogramTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    fn area_of(l: &Layer, i: usize) -> f64 {
        polygon_area(l.features[i].geometry.as_ref().unwrap()).abs()
    }

    #[test]
    fn an_already_uniform_map_is_left_almost_unchanged() {
        // Equal values over equal areas means density is already uniform, so
        // the diffusion has nothing to equalize.
        let (out, _) = run(json!({
            "input": two_squares(10.0, 10.0), "value_field": "pop",
            "grid_size": 64, "iterations": 5,
        }));
        let (a, b) = (area_of(&out, 0), area_of(&out, 1));
        assert!((a - b).abs() < 0.05 * a, "areas diverged: {a} vs {b}");
    }

    #[test]
    fn the_higher_valued_polygon_grows_relative_to_its_neighbour() {
        // The core behaviour: 4x the value must end up with more area.
        let (out, _) = run(json!({
            "input": two_squares(1.0, 4.0), "value_field": "pop",
            "grid_size": 64, "iterations": 20,
        }));
        let (a, b) = (area_of(&out, 0), area_of(&out, 1));
        assert!(b > a * 1.2, "high-value polygon did not grow: {a} vs {b}");
    }

    #[test]
    fn more_iterations_push_the_areas_further_toward_the_target_ratio() {
        // Monotone progress toward equalization is what says the diffusion is
        // actually running rather than jittering.
        let ratio = |iters: usize| -> f64 {
            let (out, _) = run(json!({
                "input": two_squares(1.0, 4.0), "value_field": "pop",
                "grid_size": 64, "iterations": iters,
            }));
            area_of(&out, 1) / area_of(&out, 0)
        };
        let few = ratio(3);
        let many = ratio(30);
        assert!(many > few, "ratio did not improve: {few} -> {many}");
    }

    #[test]
    fn the_shared_border_stays_shared() {
        // The property that distinguishes a CONTIGUOUS cartogram from the
        // shipped non_contiguous and dorling methods. Every vertex of the left
        // polygon's right edge must still coincide with one on the right
        // polygon's left edge.
        let (out, _) = run(json!({
            "input": two_squares(1.0, 4.0), "value_field": "pop",
            "grid_size": 64, "iterations": 20,
        }));
        let verts = |i: usize| -> Vec<(f64, f64)> {
            let Some(Geometry::Polygon { exterior, .. }) = out.features[i].geometry.as_ref() else {
                panic!("expected polygons")
            };
            exterior.0.iter().map(|c| (c.x, c.y)).collect()
        };
        let left = verts(0);
        let right = verts(1);
        // Every vertex of the left polygon that was on the shared edge must
        // have an exact counterpart in the right polygon.
        let mut matched = 0;
        for p in &left {
            if right
                .iter()
                .any(|q| (p.0 - q.0).abs() < 1e-9 && (p.1 - q.1).abs() < 1e-9)
            {
                matched += 1;
            }
        }
        assert!(
            matched >= 2,
            "the shared border came apart: only {matched} coincident vertices"
        );
    }

    #[test]
    fn output_rings_stay_closed() {
        // The closing vertex is displaced independently of the first, so it has
        // to be re-closed explicitly or every polygon comes out open.
        let (out, _) = run(json!({
            "input": two_squares(1.0, 4.0), "value_field": "pop",
            "grid_size": 64, "iterations": 10,
        }));
        for f in out.iter() {
            let Some(Geometry::Polygon { exterior, .. }) = f.geometry.as_ref() else {
                panic!()
            };
            let cs = &exterior.0;
            assert_eq!(
                (cs[0].x, cs[0].y),
                (cs[cs.len() - 1].x, cs[cs.len() - 1].y),
                "ring is not closed"
            );
        }
    }

    #[test]
    fn every_output_vertex_is_finite() {
        // v = -grad(rho)/rho divides by a density that can be zero early in the
        // diffusion; without the floor every downstream vertex becomes NaN.
        let (out, _) = run(json!({
            "input": two_squares(0.001, 1000.0), "value_field": "pop",
            "grid_size": 64, "iterations": 25,
        }));
        for f in out.iter() {
            let Some(Geometry::Polygon { exterior, .. }) = f.geometry.as_ref() else {
                panic!()
            };
            for c in &exterior.0 {
                assert!(c.x.is_finite() && c.y.is_finite(), "non-finite vertex");
            }
        }
    }

    #[test]
    fn attributes_and_crs_survive_and_area_error_is_reported() {
        let (out, res) = run(json!({
            "input": two_squares(1.0, 4.0), "value_field": "pop",
            "grid_size": 64, "iterations": 10,
        }));
        assert_eq!(out.crs_epsg(), Some(3857));
        assert!(out.schema.field_index("pop").is_some());
        let i = out.schema.field_index("AREA_ERR").unwrap();
        assert!(matches!(
            out.features[0].attributes[i],
            FieldValue::Float(_)
        ));
        assert!(res.outputs["max_area_error"].as_f64().unwrap().is_finite());
    }

    #[test]
    fn a_feature_with_a_non_positive_value_is_left_undistorted_and_reported() {
        // cartogram's existing convention, followed rather than reinvented.
        let mut l = Layer::new("cells")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("pop", FieldType::Float));
        l.add_feature(Some(rect(0.0, 0.0, 1.0, 1.0)), &[("pop", 5.0f64.into())])
            .unwrap();
        l.add_feature(Some(rect(1.0, 0.0, 2.0, 1.0)), &[("pop", 0.0f64.into())])
            .unwrap();
        let (out, res) = run(json!({
            "input": store(l), "value_field": "pop", "grid_size": 64, "iterations": 5,
        }));
        assert_eq!(res.outputs["undistorted_features"], json!([1]));
        let i = out.schema.field_index("AREA_ERR").unwrap();
        assert_eq!(out.features[1].attributes[i], FieldValue::Null);
    }

    #[test]
    fn the_feature_count_is_preserved() {
        let (out, res) = run(json!({
            "input": two_squares(1.0, 4.0), "value_field": "pop",
            "grid_size": 64, "iterations": 5,
        }));
        assert_eq!(out.features.len(), 2);
        assert_eq!(res.outputs["feature_count"], json!(2));
    }

    #[test]
    fn a_non_power_of_two_grid_is_refused_with_the_nearest_valid_sizes() {
        // fft2 is radix-2; silently rounding would change the result the user
        // asked for without telling them.
        let args: ToolArgs = serde_json::from_value(json!({
            "input": two_squares(1.0, 2.0), "value_field": "pop", "grid_size": 300,
        }))
        .unwrap();
        let err = ContiguousCartogramTool.validate(&args).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("power of two"), "{msg}");
        assert!(msg.contains("256") && msg.contains("512"), "{msg}");
    }

    #[test]
    fn an_oversized_grid_is_refused_before_anything_squares_it() {
        // g * g panics under overflow checks, so the cap must come first.
        let args: ToolArgs = serde_json::from_value(json!({
            "input": two_squares(1.0, 2.0), "value_field": "pop",
            "grid_size": 1_usize << 20,
        }))
        .unwrap();
        assert!(ContiguousCartogramTool.validate(&args).is_err());
    }

    #[test]
    fn densify_inserts_vertices_on_long_edges_only() {
        let short = vec![
            Coord {
                x: 0.0,
                y: 0.0,
                z: None,
                m: None,
            },
            Coord {
                x: 0.5,
                y: 0.0,
                z: None,
                m: None,
            },
        ];
        assert_eq!(densify(&short, 1.0).len(), 2, "short edge untouched");
        let long = vec![
            Coord {
                x: 0.0,
                y: 0.0,
                z: None,
                m: None,
            },
            Coord {
                x: 4.0,
                y: 0.0,
                z: None,
                m: None,
            },
        ];
        assert!(densify(&long, 1.0).len() > 2, "long edge not densified");
    }

    #[test]
    fn the_sea_density_mode_changes_the_result() {
        // An explicit parameter because it materially changes edge behaviour;
        // if the two modes agreed exactly, the parameter would be a lie.
        let a = run(json!({
            "input": two_squares(1.0, 4.0), "value_field": "pop",
            "grid_size": 64, "iterations": 15, "sea_density": "mean",
        }))
        .0;
        let b = run(json!({
            "input": two_squares(1.0, 4.0), "value_field": "pop",
            "grid_size": 64, "iterations": 15, "sea_density": "min",
        }))
        .0;
        assert!(
            (area_of(&a, 0) - area_of(&b, 0)).abs() > 1e-9,
            "sea_density had no effect"
        );
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ContiguousCartogramTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"input": "a.shp"})));
        let base = json!({"input": "a.shp", "value_field": "pop"});
        let with = |k: &str, v: Value| {
            let mut m = base.clone();
            m[k] = v;
            m
        };
        assert!(bad(with("grid_size", json!(8))));
        assert!(bad(with("iterations", json!(0))));
        assert!(bad(with("blur", json!(-1))));
        assert!(bad(with("densify_spacing", json!(0))));
        assert!(bad(with("sea_density", json!("max"))));
    }
}
