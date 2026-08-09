//! Minimal in-place complex FFT (1-D and 2-D), pure Rust, no dependencies.
//!
//! Written for `goldstein_phase_filter`, which needs the power spectrum of a
//! small square window of an interferogram. The crate deliberately avoids
//! pulling a general FFT dependency: the transform sizes here are the
//! user-chosen filter window (typically 32 or 64), always rounded up to a power
//! of two, so a textbook iterative radix-2 Cooley–Tukey kernel is both
//! sufficient and easy to audit.
//!
//! Complex values are stored as `(re, im)` pairs in a flat row-major buffer.

/// A complex sample as `(real, imaginary)`.
pub(crate) type Cpx = (f64, f64);

/// Smallest power of two `>= n` (with `next_pow2(0) == 1`).
pub(crate) fn next_pow2(n: usize) -> usize {
    let mut p = 1usize;
    while p < n {
        p <<= 1;
    }
    p
}

/// In-place radix-2 FFT over `buf`, whose length must be a power of two.
///
/// `inverse` runs the conjugate transform **and** divides by `n`, so
/// `ifft(fft(x)) == x`.
fn fft1(buf: &mut [Cpx], inverse: bool) {
    let n = buf.len();
    debug_assert!(n.is_power_of_two(), "fft1 length must be a power of two");
    if n <= 1 {
        return;
    }

    // Bit-reversal permutation.
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            buf.swap(i, j);
        }
    }

    // Butterfly stages.
    let sign = if inverse { 1.0 } else { -1.0 };
    let mut len = 2usize;
    while len <= n {
        let ang = sign * 2.0 * std::f64::consts::PI / len as f64;
        let (wr, wi) = (ang.cos(), ang.sin());
        let mut i = 0usize;
        while i < n {
            // Recurrence-free twiddle: recomputing per block keeps the error
            // bounded, which matters because the filter squares the magnitudes.
            let (mut cr, mut ci) = (1.0f64, 0.0f64);
            for k in 0..len / 2 {
                let u = buf[i + k];
                let v = buf[i + k + len / 2];
                let vr = v.0 * cr - v.1 * ci;
                let vi = v.0 * ci + v.1 * cr;
                buf[i + k] = (u.0 + vr, u.1 + vi);
                buf[i + k + len / 2] = (u.0 - vr, u.1 - vi);
                let ncr = cr * wr - ci * wi;
                ci = cr * wi + ci * wr;
                cr = ncr;
                let _ = k;
            }
            i += len;
        }
        len <<= 1;
    }

    if inverse {
        let inv = 1.0 / n as f64;
        for c in buf.iter_mut() {
            c.0 *= inv;
            c.1 *= inv;
        }
    }
}

/// In-place 2-D FFT of a `rows x cols` row-major complex buffer.
///
/// Both dimensions must be powers of two. Rows are transformed first, then
/// columns (the transform is separable, so the order does not matter).
pub(crate) fn fft2(buf: &mut [Cpx], rows: usize, cols: usize, inverse: bool) {
    debug_assert_eq!(buf.len(), rows * cols);
    for r in 0..rows {
        fft1(&mut buf[r * cols..(r + 1) * cols], inverse);
    }
    let mut col = vec![(0.0, 0.0); rows];
    for c in 0..cols {
        for (r, slot) in col.iter_mut().enumerate() {
            *slot = buf[r * cols + c];
        }
        fft1(&mut col, inverse);
        for (r, slot) in col.iter().enumerate() {
            buf[r * cols + c] = *slot;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    /// A constant field transforms to a single DC spike.
    #[test]
    fn constant_field_is_dc_only() {
        let (rows, cols) = (8, 8);
        let mut buf = vec![(3.0, 0.0); rows * cols];
        fft2(&mut buf, rows, cols, false);
        assert!(close(buf[0].0, 3.0 * 64.0, 1e-9), "DC = {:?}", buf[0]);
        for (i, c) in buf.iter().enumerate().skip(1) {
            assert!(close(c.0, 0.0, 1e-9) && close(c.1, 0.0, 1e-9), "bin {i} {c:?}");
        }
    }

    /// A single-cycle cosine along x puts all energy in bins (0, ±1).
    #[test]
    fn single_tone_lands_in_one_bin() {
        let (rows, cols) = (4, 16);
        let mut buf = vec![(0.0, 0.0); rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let x = 2.0 * std::f64::consts::PI * c as f64 / cols as f64;
                buf[r * cols + c] = (x.cos(), 0.0);
            }
        }
        fft2(&mut buf, rows, cols, false);
        // Energy concentrated at (row 0, col 1) and its conjugate (0, cols-1).
        let mag = |i: usize| (buf[i].0 * buf[i].0 + buf[i].1 * buf[i].1).sqrt();
        assert!(mag(1) > 1.0, "positive-frequency bin empty");
        assert!(mag(cols - 1) > 1.0, "negative-frequency bin empty");
        assert!(close(mag(0), 0.0, 1e-9), "unexpected DC");
        assert!(close(mag(2), 0.0, 1e-9), "energy leaked to bin 2");
    }

    /// The inverse transform recovers the input.
    #[test]
    fn round_trip() {
        let (rows, cols) = (8, 4);
        let orig: Vec<Cpx> = (0..rows * cols)
            .map(|i| ((i as f64).sin(), (i as f64 * 0.7).cos()))
            .collect();
        let mut buf = orig.clone();
        fft2(&mut buf, rows, cols, false);
        fft2(&mut buf, rows, cols, true);
        for (i, (got, want)) in buf.iter().zip(orig.iter()).enumerate() {
            assert!(
                close(got.0, want.0, 1e-9) && close(got.1, want.1, 1e-9),
                "sample {i}: {got:?} != {want:?}"
            );
        }
    }

    #[test]
    fn pow2_rounding() {
        assert_eq!(next_pow2(0), 1);
        assert_eq!(next_pow2(1), 1);
        assert_eq!(next_pow2(5), 8);
        assert_eq!(next_pow2(32), 32);
        assert_eq!(next_pow2(33), 64);
    }
}
