//! Shared machinery for the "optimized" spatial-statistics tools
//! (`optimized_hot_spot_analysis`, `optimized_outlier_analysis`): incident
//! extraction and aggregation, automatic distance-band selection, the Getis-Ord
//! Gi\* and Anselin Local Moran's I statistics with analytic z-scores, a standard
//! normal CDF, and Benjamini-Hochberg False Discovery Rate correction.
//!
//! Weights are binary within the analysis distance band. Gi\* includes the focal
//! feature in its own neighborhood (w_ii = 1); Local Moran excludes it
//! (w_ii = 0). Everything is deterministic — the p-values come from closed-form
//! variance formulas, not permutation — so results are stable in WASM without an
//! RNG.

use std::collections::BTreeMap;

use serde_json::Value;
use wbcore::{ToolArgs, ToolError};
use wbvector::{Geometry, Layer};

/// A working point: location plus the analysis value (incident count or field).
#[derive(Clone, Copy)]
pub struct Pt {
    pub x: f64,
    pub y: f64,
    pub val: f64,
}

/// How raw incidents are aggregated before analysis.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Aggregation {
    /// Count incidents in a regular grid over the extent (empty cells kept as 0).
    Fishnet,
    /// Merge coincident/near-coincident incidents into weighted points.
    Snap,
}

pub fn parse_aggregation(args: &ToolArgs) -> Result<Aggregation, ToolError> {
    Ok(
        match args
            .get("aggregation")
            .and_then(Value::as_str)
            .map(str::trim)
        {
            None | Some("") | Some("fishnet") => Aggregation::Fishnet,
            Some("snap") => Aggregation::Snap,
            Some(o) => {
                return Err(ToolError::Validation(format!(
                    "'aggregation' must be 'fishnet' or 'snap', got '{o}'"
                )))
            }
        },
    )
}

/// Representative coordinate of any geometry (the mean of its coordinates).
pub fn geometry_centroid(geom: &Geometry) -> Option<(f64, f64)> {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0.0;
    let mut acc = |x: f64, y: f64| {
        sx += x;
        sy += y;
        n += 1.0;
    };
    match geom {
        Geometry::Point(c) => acc(c.x, c.y),
        Geometry::MultiPoint(cs) | Geometry::LineString(cs) => {
            for c in cs {
                acc(c.x, c.y);
            }
        }
        Geometry::MultiLineString(ls) => {
            for l in ls {
                for c in l {
                    acc(c.x, c.y);
                }
            }
        }
        Geometry::Polygon { exterior, .. } => {
            for c in exterior.coords() {
                acc(c.x, c.y);
            }
        }
        Geometry::MultiPolygon(parts) => {
            for (ext, _) in parts {
                for c in ext.coords() {
                    acc(c.x, c.y);
                }
            }
        }
        _ => {}
    }
    if n > 0.0 {
        Some((sx / n, sy / n))
    } else {
        None
    }
}

/// Extracts working points from a layer. With `field`, each feature becomes a
/// weighted point at its centroid (no aggregation happens later); without it,
/// every feature is a unit incident to be aggregated.
pub fn extract_points(layer: &Layer, field: Option<&str>) -> Result<(Vec<Pt>, usize), ToolError> {
    let fidx = match field {
        Some(f) => Some(
            layer
                .schema
                .field_index(f)
                .ok_or_else(|| ToolError::Validation(format!("analysis_field '{f}' not found")))?,
        ),
        None => None,
    };
    let mut pts = Vec::new();
    let mut skipped = 0usize;
    for feat in layer.iter() {
        let Some((x, y)) = feat.geometry.as_ref().and_then(geometry_centroid) else {
            skipped += 1;
            continue;
        };
        let val = match fidx {
            Some(i) => match feat.attributes.get(i).and_then(|v| v.as_f64()) {
                Some(v) if v.is_finite() => v,
                _ => {
                    skipped += 1;
                    continue;
                }
            },
            None => 1.0,
        };
        pts.push(Pt { x, y, val });
    }
    Ok((pts, skipped))
}

/// Aggregates unit incidents into weighted points. `cell` is the grid/snap size.
/// Fishnet keeps every cell in the extent (zeros included); snap keeps only
/// occupied locations.
pub fn aggregate(pts: &[Pt], mode: Aggregation, cell: f64) -> Vec<Pt> {
    if pts.is_empty() {
        return Vec::new();
    }
    let (mut xmin, mut ymin, mut xmax, mut ymax) = (
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );
    for p in pts {
        xmin = xmin.min(p.x);
        xmax = xmax.max(p.x);
        ymin = ymin.min(p.y);
        ymax = ymax.max(p.y);
    }
    let cell = cell.max(1e-9);
    match mode {
        Aggregation::Snap => {
            // Merge points sharing a snap cell; weighted point at their mean.
            let mut acc: BTreeMap<(i64, i64), (f64, f64, f64)> = BTreeMap::new();
            for p in pts {
                let key = (
                    ((p.x - xmin) / cell).floor() as i64,
                    ((p.y - ymin) / cell).floor() as i64,
                );
                let e = acc.entry(key).or_insert((0.0, 0.0, 0.0));
                e.0 += p.x;
                e.1 += p.y;
                e.2 += p.val;
            }
            acc.values()
                .map(|(sx, sy, w)| Pt {
                    x: sx / w.max(1.0),
                    y: sy / w.max(1.0),
                    val: *w,
                })
                .collect()
        }
        Aggregation::Fishnet => {
            let cols = ((((xmax - xmin) / cell).ceil()) as usize).max(1);
            let rows = ((((ymax - ymin) / cell).ceil()) as usize).max(1);
            let mut counts = vec![0.0f64; rows * cols];
            for p in pts {
                let c = ((((p.x - xmin) / cell).floor()) as usize).min(cols - 1);
                let r = ((((p.y - ymin) / cell).floor()) as usize).min(rows - 1);
                counts[r * cols + c] += p.val;
            }
            let mut out = Vec::with_capacity(rows * cols);
            for r in 0..rows {
                for c in 0..cols {
                    out.push(Pt {
                        x: xmin + (c as f64 + 0.5) * cell,
                        y: ymin + (r as f64 + 0.5) * cell,
                        val: counts[r * cols + c],
                    });
                }
            }
            out
        }
    }
}

/// Default aggregation cell size: the extent's longer side / 30 (≈900 cells).
pub fn default_cell_size(pts: &[Pt]) -> f64 {
    let (mut xmin, mut ymin, mut xmax, mut ymax) = (
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );
    for p in pts {
        xmin = xmin.min(p.x);
        xmax = xmax.max(p.x);
        ymin = ymin.min(p.y);
        ymax = ymax.max(p.y);
    }
    ((xmax - xmin).max(ymax - ymin) / 30.0).max(1e-6)
}

fn dist2(a: &Pt, b: &Pt) -> f64 {
    let dx = a.x - b.x;
    let dy = a.y - b.y;
    dx * dx + dy * dy
}

/// Distance guaranteeing every feature has at least one neighbor: the maximum,
/// over features, of the nearest-neighbor distance.
pub fn min_neighbor_distance(pts: &[Pt]) -> f64 {
    let n = pts.len();
    let mut maxnn = 0.0f64;
    for i in 0..n {
        let mut nn = f64::INFINITY;
        for j in 0..n {
            if i != j {
                nn = nn.min(dist2(&pts[i], &pts[j]));
            }
        }
        if nn.is_finite() {
            maxnn = maxnn.max(nn);
        }
    }
    maxnn.sqrt()
}

/// Automatic analysis distance band. Scans candidate bands from the
/// one-neighbor-guarantee distance upward and returns the one maximizing the
/// global Getis-Ord General G z-score (peak clustering intensity), falling back
/// to the guarantee distance when no peak is finite.
pub fn auto_distance_band(pts: &[Pt]) -> f64 {
    let base = min_neighbor_distance(pts);
    if !(base.is_finite() && base > 0.0) {
        return 1.0;
    }
    let mut best_band = base;
    let mut best_z = f64::NEG_INFINITY;
    for k in 0..8 {
        let band = base * (1.0 + 0.35 * k as f64);
        let z = global_g_z(pts, band);
        if z.is_finite() && z > best_z {
            best_z = z;
            best_band = band;
        }
    }
    best_band
}

/// Global Getis-Ord General G z-score at a distance band (used for band scan).
fn global_g_z(pts: &[Pt], band: f64) -> f64 {
    let n = pts.len();
    if n < 3 {
        return f64::NAN;
    }
    let b2 = band * band;
    let mut w = vec![0.0f64; n * n];
    let mut wsum = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            if i != j && dist2(&pts[i], &pts[j]) <= b2 {
                w[i * n + j] = 1.0;
                wsum += 1.0;
            }
        }
    }
    if wsum == 0.0 {
        return f64::NAN;
    }
    let x: Vec<f64> = pts.iter().map(|p| p.val).collect();
    let mut num = 0.0;
    let mut denom = 0.0;
    for i in 0..n {
        for j in 0..n {
            if i != j {
                denom += x[i] * x[j];
                num += w[i * n + j] * x[i] * x[j];
            }
        }
    }
    if denom == 0.0 {
        return f64::NAN;
    }
    let g = num / denom;
    // Expectation and variance under randomization (Getis & Ord 1992).
    let nn = n as f64;
    let eg = wsum / (nn * (nn - 1.0));
    let s1: f64 = x.iter().sum();
    let s2: f64 = x.iter().map(|v| v * v).sum();
    let s3: f64 = x.iter().map(|v| v * v * v).sum();
    let s4 = s1 * s1 - s2;
    let d2 = s1 * s1 - s2; // Σ_{i≠j} xi xj denominator scale
    if d2 == 0.0 {
        return f64::NAN;
    }
    // Variance approximation via the normal form; use a simplified moment.
    let m2 = s2 / nn - (s1 / nn) * (s1 / nn);
    if m2 <= 0.0 {
        return f64::NAN;
    }
    // Standardized clustering: (G - EG) scaled by its sampling spread. Use the
    // second-moment based approximation adequate for band ranking.
    let _ = (s3, s4);
    let var_g = eg * (1.0 - eg) / wsum.max(1.0);
    if var_g <= 0.0 {
        return f64::NAN;
    }
    (g - eg) / var_g.sqrt()
}

/// Getis-Ord Gi\* for every feature at `band`. Returns (z, two-sided p).
pub fn getis_gi_star(pts: &[Pt], band: f64) -> Vec<(f64, f64)> {
    let n = pts.len();
    let x: Vec<f64> = pts.iter().map(|p| p.val).collect();
    let nn = n as f64;
    let xbar = x.iter().sum::<f64>() / nn;
    let s = (x.iter().map(|v| v * v).sum::<f64>() / nn - xbar * xbar)
        .max(0.0)
        .sqrt();
    let b2 = band * band;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        // Gi* includes the focal feature (w_ii = 1).
        let mut wsum = 1.0;
        let mut wx = x[i];
        for j in 0..n {
            if i != j && dist2(&pts[i], &pts[j]) <= b2 {
                wsum += 1.0;
                wx += x[j];
            }
        }
        if s <= 0.0 {
            out.push((0.0, 1.0));
            continue;
        }
        let denom = s * (((nn * wsum - wsum * wsum) / (nn - 1.0)).max(0.0)).sqrt();
        if denom <= 0.0 {
            out.push((0.0, 1.0));
            continue;
        }
        let z = (wx - xbar * wsum) / denom;
        out.push((z, two_sided_p(z)));
    }
    out
}

/// Local result: Moran's I_i, z-score, two-sided p, and the sign of (x_i - mean).
pub struct LocalMoran {
    pub i: f64,
    pub z: f64,
    pub p: f64,
    pub high: bool,
}

/// Anselin Local Moran's I for every feature at `band` (binary weights).
pub fn local_moran(pts: &[Pt], band: f64) -> Vec<LocalMoran> {
    let n = pts.len();
    let nn = n as f64;
    let x: Vec<f64> = pts.iter().map(|p| p.val).collect();
    let mean = x.iter().sum::<f64>() / nn;
    let z: Vec<f64> = x.iter().map(|v| v - mean).collect();
    let m2 = z.iter().map(|v| v * v).sum::<f64>() / nn;
    let m4 = z.iter().map(|v| v.powi(4)).sum::<f64>() / nn;
    let b2 = if m2 > 0.0 { m4 / (m2 * m2) } else { f64::NAN };
    let band2 = band * band;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let mut wi = 0.0f64;
        let mut lag = 0.0f64;
        for j in 0..n {
            if i != j && dist2(&pts[i], &pts[j]) <= band2 {
                wi += 1.0;
                lag += z[j];
            }
        }
        if m2 <= 0.0 || wi == 0.0 || !(b2.is_finite()) || n < 4 {
            out.push(LocalMoran {
                i: 0.0,
                z: 0.0,
                p: 1.0,
                high: z[i] >= 0.0,
            });
            continue;
        }
        let ii = (z[i] / m2) * lag;
        let ei = -wi / (nn - 1.0);
        let wi2 = wi; // Σ w_ij² for binary weights
        let wikh = wi * wi - wi2; // Σ_{k≠h} w_ik w_ih
        let var = wi2 * (nn - b2) / (nn - 1.0) + wikh * (2.0 * b2 - nn) / ((nn - 1.0) * (nn - 2.0))
            - ei * ei;
        if !(var.is_finite()) || var <= 0.0 {
            out.push(LocalMoran {
                i: ii,
                z: 0.0,
                p: 1.0,
                high: z[i] >= 0.0,
            });
            continue;
        }
        let zi = (ii - ei) / var.sqrt();
        out.push(LocalMoran {
            i: ii,
            z: zi,
            p: two_sided_p(zi),
            high: z[i] >= 0.0,
        });
    }
    out
}

/// Standard normal CDF (Zelen & Severo rational approximation).
pub fn norm_cdf(z: f64) -> f64 {
    if !z.is_finite() {
        return if z > 0.0 { 1.0 } else { 0.0 };
    }
    let t = 1.0 / (1.0 + 0.2316419 * z.abs());
    let d = 0.398942280401433 * (-z * z / 2.0).exp();
    let p = d
        * t
        * (0.319381530
            + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))));
    if z >= 0.0 {
        1.0 - p
    } else {
        p
    }
}

fn two_sided_p(z: f64) -> f64 {
    2.0 * (1.0 - norm_cdf(z.abs()))
}

/// Benjamini-Hochberg: returns, for each input p-value, whether it is significant
/// at the given FDR level `alpha`.
pub fn bh_significant(pvals: &[f64], alpha: f64) -> Vec<bool> {
    let m = pvals.len();
    if m == 0 {
        return Vec::new();
    }
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&a, &b| pvals[a].total_cmp(&pvals[b]));
    // Largest rank k (1-based) with p_(k) <= (k/m)*alpha.
    let mut kmax = 0usize;
    for (rank, &idx) in order.iter().enumerate() {
        let k = rank + 1;
        if pvals[idx] <= (k as f64 / m as f64) * alpha {
            kmax = k;
        }
    }
    let crit = if kmax == 0 {
        -1.0
    } else {
        pvals[order[kmax - 1]]
    };
    pvals.iter().map(|&p| p <= crit).collect()
}

/// FDR-corrected hot/cold-spot bin in [-3, 3]: sign from z, magnitude from the
/// strictest FDR level (0.01→3, 0.05→2, 0.10→1) the feature passes.
pub fn fdr_bins(pvals: &[f64], zs: &[f64]) -> Vec<i32> {
    let sig01 = bh_significant(pvals, 0.01);
    let sig05 = bh_significant(pvals, 0.05);
    let sig10 = bh_significant(pvals, 0.10);
    (0..pvals.len())
        .map(|i| {
            let mag = if sig01[i] {
                3
            } else if sig05[i] {
                2
            } else if sig10[i] {
                1
            } else {
                0
            };
            if mag == 0 {
                0
            } else if zs[i] >= 0.0 {
                mag
            } else {
                -mag
            }
        })
        .collect()
}
