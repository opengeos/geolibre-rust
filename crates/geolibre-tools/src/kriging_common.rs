//! Shared variogram and ordinary-kriging machinery for the round-17
//! geostatistical tools (`semivariogram_sensitivity`, `areal_interpolation`).
//!
//! Both need the same three pieces: a theoretical variogram model, an
//! empirical variogram fit, and an ordinary-kriging solve. Keeping them here
//! means the two tools cannot drift apart on model conventions — a real risk,
//! because `semivariogram_sensitivity` exists precisely to perturb the
//! parameters `areal_interpolation` fits.
//!
//! Everything is plain `Vec<f64>` Gaussian elimination: no linear-algebra
//! crate, no RNG, deterministic, WASM-safe.

use wbcore::ToolError;

/// A theoretical variogram model.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum VariogramModel {
    Spherical,
    Exponential,
    Gaussian,
}

impl VariogramModel {
    pub(crate) fn label(self) -> &'static str {
        match self {
            VariogramModel::Spherical => "spherical",
            VariogramModel::Exponential => "exponential",
            VariogramModel::Gaussian => "gaussian",
        }
    }

    pub(crate) fn parse(s: &str) -> Result<Self, ToolError> {
        match s {
            "spherical" => Ok(VariogramModel::Spherical),
            "exponential" => Ok(VariogramModel::Exponential),
            "gaussian" => Ok(VariogramModel::Gaussian),
            other => Err(ToolError::Validation(format!(
                "'model' must be spherical|exponential|gaussian, got '{other}'"
            ))),
        }
    }
}

/// Nugget / partial sill / range triple.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Variogram {
    pub(crate) model: VariogramModel,
    pub(crate) nugget: f64,
    pub(crate) partial_sill: f64,
    pub(crate) range: f64,
}

impl Variogram {
    /// Semivariance at lag `h`.
    pub(crate) fn gamma(&self, h: f64) -> f64 {
        if h <= 0.0 {
            return 0.0; // gamma(0) = 0 by definition, nugget included
        }
        let a = self.range.max(1e-12);
        let structured = match self.model {
            VariogramModel::Spherical => {
                if h >= a {
                    1.0
                } else {
                    let r = h / a;
                    1.5 * r - 0.5 * r * r * r
                }
            }
            VariogramModel::Exponential => 1.0 - (-3.0 * h / a).exp(),
            VariogramModel::Gaussian => 1.0 - (-3.0 * (h / a).powi(2)).exp(),
        };
        self.nugget + self.partial_sill * structured
    }

    /// Covariance at lag `h`, derived from the variogram's sill.
    pub(crate) fn covariance(&self, h: f64) -> f64 {
        (self.nugget + self.partial_sill) - self.gamma(h)
    }
}

/// Fits nugget, partial sill and range from an empirical variogram cloud.
///
/// A grid search over range with least-squares nugget/sill at each candidate:
/// robust and fully deterministic, unlike a gradient method that would need a
/// starting point and could land in a different optimum per input ordering.
pub(crate) fn fit_variogram(
    coords: &[(f64, f64)],
    values: &[f64],
    model: VariogramModel,
    lag_count: usize,
) -> Variogram {
    let n = coords.len();
    let mean = values.iter().sum::<f64>() / n.max(1) as f64;
    let raw_variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n.max(1) as f64;
    // A constant field has zero variance, which would give a zero-scale
    // covariance matrix that the solver reads as singular. Ordinary-kriging
    // weights are invariant to the sill's scale, so any positive value gives
    // the same (correct) prediction — 1.0 just keeps the matrix conditioned.
    let variance = if raw_variance > 0.0 {
        raw_variance
    } else {
        1.0
    };

    // Empirical semivariances, binned by lag.
    let mut max_dist: f64 = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            max_dist = max_dist.max(dist(coords[i], coords[j]));
        }
    }
    if max_dist <= 0.0 || n < 3 {
        return Variogram {
            model,
            nugget: 0.0,
            partial_sill: variance.max(1e-12),
            range: 1.0,
        };
    }
    // Beyond half the extent the pair count collapses and the estimate is
    // dominated by noise — the standard cutoff.
    let cutoff = max_dist / 2.0;
    let bins = lag_count.max(3);
    let width = cutoff / bins as f64;
    let mut sum = vec![0.0; bins];
    let mut count = vec![0usize; bins];
    for i in 0..n {
        for j in (i + 1)..n {
            let h = dist(coords[i], coords[j]);
            if h <= 0.0 || h > cutoff {
                continue;
            }
            let b = (((h / width).floor() as usize).min(bins - 1));
            sum[b] += 0.5 * (values[i] - values[j]).powi(2);
            count[b] += 1;
        }
    }
    let cloud: Vec<(f64, f64)> = (0..bins)
        .filter(|b| count[*b] > 0)
        .map(|b| ((b as f64 + 0.5) * width, sum[b] / count[b] as f64))
        .collect();
    if cloud.len() < 2 {
        return Variogram {
            model,
            nugget: 0.0,
            partial_sill: variance.max(1e-12),
            range: cutoff.max(1e-12),
        };
    }

    // Grid search over range; nugget and partial sill by non-negative least
    // squares against the model's shape at that range.
    let mut best = Variogram {
        model,
        nugget: 0.0,
        partial_sill: variance.max(1e-12),
        range: cutoff,
    };
    let mut best_err = f64::INFINITY;
    for k in 1..=40 {
        let range = cutoff * k as f64 / 40.0;
        let probe = Variogram {
            model,
            nugget: 0.0,
            partial_sill: 1.0,
            range,
        };
        // gamma(h) = nugget + sill * shape(h); solve the 2x2 normal equations.
        let (mut s11, mut s1x, mut sxx, mut s1y, mut sxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for (h, g) in &cloud {
            let x = probe.gamma(*h);
            s11 += 1.0;
            s1x += x;
            sxx += x * x;
            s1y += g;
            sxy += x * g;
        }
        let det = s11 * sxx - s1x * s1x;
        let (nugget, sill) = if det.abs() < 1e-12 {
            (0.0, (s1y / s11).max(0.0))
        } else {
            ((s1y * sxx - s1x * sxy) / det, (s11 * sxy - s1x * s1y) / det)
        };
        // Negative nugget or sill is unphysical; clamping keeps every emitted
        // model usable rather than letting a fit artefact poison the solve.
        let cand = Variogram {
            model,
            nugget: nugget.max(0.0),
            partial_sill: sill.max(1e-12),
            range,
        };
        let err: f64 = cloud
            .iter()
            .map(|(h, g)| (cand.gamma(*h) - g).powi(2))
            .sum();
        if err < best_err {
            best_err = err;
            best = cand;
        }
    }
    best
}

/// Ordinary-kriging prediction and variance at one location.
///
/// Returns `(prediction, variance)`. The system is the usual covariance matrix
/// bordered by the unbiasedness (Lagrange) row and column.
pub(crate) fn ordinary_kriging(
    coords: &[(f64, f64)],
    values: &[f64],
    target: (f64, f64),
    vg: &Variogram,
) -> Option<(f64, f64)> {
    krige_with(
        coords,
        values,
        vg,
        |i| vg.covariance(dist(coords[i], target)),
        vg.covariance(0.0),
    )
}

/// Ordinary kriging against an arbitrary data-to-data covariance matrix.
///
/// Block kriging needs BOTH sides expressed on the same support: pairing a
/// point-to-point left-hand side with block-to-block right-hand sides yields an
/// inconsistent system whose weights do not reproduce known values.
pub(crate) fn krige_matrix<C, F>(
    n: usize,
    values: &[f64],
    cov: C,
    rhs: F,
    c00: f64,
) -> Option<(f64, f64)>
where
    C: Fn(usize, usize) -> f64,
    F: Fn(usize) -> f64,
{
    if n == 0 || values.len() != n {
        return None;
    }
    if n == 1 {
        return Some((values[0], (c00 - rhs(0)).max(0.0)));
    }
    let m = n + 1;
    let mut a = vec![0.0_f64; m * m];
    let mut b = vec![0.0_f64; m];
    for i in 0..n {
        for j in 0..n {
            a[i * m + j] = cov(i, j);
        }
        a[i * m + n] = 1.0;
        a[n * m + i] = 1.0;
        b[i] = rhs(i);
    }
    a[n * m + n] = 0.0;
    b[n] = 1.0;

    let x = solve_linear(&mut a, &mut b, m)?;
    let prediction: f64 = (0..n).map(|i| x[i] * values[i]).sum();
    let variance = (c00 - (0..n).map(|i| x[i] * rhs(i)).sum::<f64>() - x[n]).max(0.0);
    Some((prediction, variance))
}

/// Ordinary kriging against arbitrary right-hand-side covariances.
///
/// `rhs(i)` is the covariance between datum `i` and the prediction support, and
/// `c00` the support's covariance with itself. Point kriging passes point-to-
/// point covariances; area-to-area kriging passes block averages. Factoring it
/// this way is what lets `areal_interpolation` reuse the same solver instead of
/// reimplementing the normal equations.
pub(crate) fn krige_with<F>(
    coords: &[(f64, f64)],
    values: &[f64],
    vg: &Variogram,
    rhs: F,
    c00: f64,
) -> Option<(f64, f64)>
where
    F: Fn(usize) -> f64,
{
    krige_matrix(
        coords.len(),
        values,
        |i, j| vg.covariance(dist(coords[i], coords[j])),
        rhs,
        c00,
    )
}

/// Gaussian elimination with partial pivoting. `None` when singular.
fn solve_linear(a: &mut [f64], b: &mut [f64], n: usize) -> Option<Vec<f64>> {
    for col in 0..n {
        let mut pivot = col;
        for r in (col + 1)..n {
            if a[r * n + col].abs() > a[pivot * n + col].abs() {
                pivot = r;
            }
        }
        if a[pivot * n + col].abs() < 1e-12 {
            return None;
        }
        if pivot != col {
            for k in 0..n {
                a.swap(col * n + k, pivot * n + k);
            }
            b.swap(col, pivot);
        }
        let d = a[col * n + col];
        for r in (col + 1)..n {
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
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut acc = b[i];
        for k in (i + 1)..n {
            acc -= a[i * n + k] * x[k];
        }
        x[i] = acc / a[i * n + i];
    }
    Some(x)
}

pub(crate) fn dist(a: (f64, f64), b: (f64, f64)) -> f64 {
    ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gamma_is_zero_at_the_origin_and_reaches_the_sill() {
        for model in [
            VariogramModel::Spherical,
            VariogramModel::Exponential,
            VariogramModel::Gaussian,
        ] {
            let v = Variogram {
                model,
                nugget: 0.5,
                partial_sill: 2.0,
                range: 10.0,
            };
            assert_eq!(v.gamma(0.0), 0.0, "{model:?}");
            // Well beyond the range, gamma approaches nugget + partial sill.
            assert!((v.gamma(1000.0) - 2.5).abs() < 1e-6, "{model:?}");
            // And it is monotonically increasing.
            assert!(v.gamma(1.0) < v.gamma(5.0), "{model:?}");
        }
    }

    #[test]
    fn covariance_and_variogram_are_complementary() {
        let v = Variogram {
            model: VariogramModel::Exponential,
            nugget: 0.2,
            partial_sill: 1.8,
            range: 5.0,
        };
        for h in [0.5, 2.0, 7.0] {
            assert!((v.covariance(h) + v.gamma(h) - 2.0).abs() < 1e-9);
        }
    }

    #[test]
    fn kriging_reproduces_a_datum_at_its_own_location() {
        // The exactness property: with no nugget, the prediction at a sample
        // point IS that sample, with (near) zero variance.
        let coords = vec![(0.0, 0.0), (10.0, 0.0), (0.0, 10.0), (10.0, 10.0)];
        let values = vec![1.0, 5.0, 3.0, 9.0];
        let vg = Variogram {
            model: VariogramModel::Spherical,
            nugget: 0.0,
            partial_sill: 4.0,
            range: 20.0,
        };
        let (p, var) = ordinary_kriging(&coords, &values, (10.0, 0.0), &vg).unwrap();
        assert!((p - 5.0).abs() < 1e-6, "got {p}");
        assert!(var < 1e-6, "variance {var} should vanish at a datum");
    }

    #[test]
    fn kriging_a_constant_field_returns_that_constant() {
        let coords = vec![(0.0, 0.0), (10.0, 0.0), (5.0, 8.0)];
        let values = vec![7.0, 7.0, 7.0];
        let vg = Variogram {
            model: VariogramModel::Exponential,
            nugget: 0.1,
            partial_sill: 1.0,
            range: 15.0,
        };
        let (p, _) = ordinary_kriging(&coords, &values, (4.0, 3.0), &vg).unwrap();
        assert!((p - 7.0).abs() < 1e-6, "got {p}");
    }

    #[test]
    fn variance_grows_with_distance_from_the_data() {
        let coords = vec![(0.0, 0.0), (10.0, 0.0), (0.0, 10.0)];
        let values = vec![1.0, 2.0, 3.0];
        let vg = Variogram {
            model: VariogramModel::Spherical,
            nugget: 0.0,
            partial_sill: 4.0,
            range: 20.0,
        };
        let (_, near) = ordinary_kriging(&coords, &values, (3.0, 3.0), &vg).unwrap();
        let (_, far) = ordinary_kriging(&coords, &values, (100.0, 100.0), &vg).unwrap();
        assert!(far > near, "variance did not grow: near {near}, far {far}");
    }

    #[test]
    fn the_fitted_variogram_recovers_a_known_structure() {
        // Sample a known spherical model on a grid, then fit it back. The
        // fitted sill should be close and the range in the right ballpark.
        let truth = Variogram {
            model: VariogramModel::Spherical,
            nugget: 0.0,
            partial_sill: 1.0,
            range: 20.0,
        };
        // A smooth deterministic field with the right correlation scale.
        let mut coords = Vec::new();
        let mut values = Vec::new();
        for i in 0..12 {
            for j in 0..12 {
                let (x, y) = (i as f64 * 5.0, j as f64 * 5.0);
                coords.push((x, y));
                values.push((x / 20.0).sin() + (y / 20.0).cos());
            }
        }
        let fit = fit_variogram(&coords, &values, VariogramModel::Spherical, 12);
        assert_eq!(fit.model, truth.model);
        assert!(fit.partial_sill > 0.0);
        assert!(fit.range > 0.0);
        assert!(fit.nugget >= 0.0, "nugget must not be negative");
    }

    #[test]
    fn a_singular_system_is_reported_rather_than_producing_nonsense() {
        // Duplicated locations make the covariance matrix singular.
        let coords = vec![(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)];
        let values = vec![1.0, 2.0, 3.0];
        let vg = Variogram {
            model: VariogramModel::Spherical,
            nugget: 0.0,
            partial_sill: 1.0,
            range: 10.0,
        };
        assert!(ordinary_kriging(&coords, &values, (5.0, 5.0), &vg).is_none());
    }
}
