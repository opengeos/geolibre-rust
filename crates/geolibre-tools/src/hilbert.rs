//! Hilbert-curve spatial ordering for vector features.
//!
//! Sorting features along a Hilbert curve clusters spatially-near records on
//! disk, which is what makes a GeoParquet file's row-group and page statistics
//! (and a bbox covering) effective for spatial pruning. We map each feature's
//! bounding-box center into a `2^ORDER` x `2^ORDER` grid over the dataset
//! extent, then sort by the curve distance of that cell.

/// Grid resolution per axis is `2^ORDER`. 16 gives 65536 cells per axis, the
/// de-facto standard used by GeoParquet writers (e.g. GeoPandas, ogr2ogr).
const ORDER: u32 = 16;

/// Maps grid coordinates `(x, y)` in `[0, 2^ORDER)` to their Hilbert-curve
/// distance `d`.
pub fn xy_to_hilbert(x: u32, y: u32) -> u64 {
    xy_to_hilbert_order(x, y, ORDER)
}

/// Classic iterative xy->d Hilbert transform for an arbitrary grid `order`
/// (grid is `2^order` per axis). Factored out so tests can exhaustively check a
/// small grid.
pub(crate) fn xy_to_hilbert_order(mut x: u32, mut y: u32, order: u32) -> u64 {
    let n: u32 = 1 << order;
    let mut d: u64 = 0;
    let mut s: u32 = n / 2;
    while s > 0 {
        let rx = u32::from((x & s) > 0);
        let ry = u32::from((y & s) > 0);
        d += (s as u64) * (s as u64) * ((3 * rx) ^ ry) as u64;
        // Rotate the quadrant so the curve stays continuous.
        if ry == 0 {
            if rx == 1 {
                x = s.wrapping_sub(1).wrapping_sub(x) & (n - 1);
                y = s.wrapping_sub(1).wrapping_sub(y) & (n - 1);
            }
            std::mem::swap(&mut x, &mut y);
        }
        s /= 2;
    }
    d
}

/// Computes the Hilbert distance for a point at `(px, py)` within the extent
/// `[min_x, max_x] x [min_y, max_y]`. Degenerate (zero-width) axes map to the
/// grid center, so a constant-coordinate dataset still sorts deterministically.
pub fn hilbert_for_point(px: f64, py: f64, min_x: f64, min_y: f64, max_x: f64, max_y: f64) -> u64 {
    let max_cell = ((1u32 << ORDER) - 1) as f64;
    let scale = |v: f64, lo: f64, hi: f64| -> u32 {
        if hi <= lo {
            (max_cell / 2.0) as u32
        } else {
            (((v - lo) / (hi - lo)) * max_cell).round().clamp(0.0, max_cell) as u32
        }
    };
    xy_to_hilbert(scale(px, min_x, max_x), scale(py, min_y, max_y))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_a_bijection_over_a_small_grid() {
        // Over a 4x4 grid the distances must be exactly 0..16 with no repeats.
        let order = 2;
        let n = 1u32 << order;
        let mut seen = vec![false; (n * n) as usize];
        for y in 0..n {
            for x in 0..n {
                let d = xy_to_hilbert_order(x, y, order) as usize;
                assert!(d < seen.len(), "distance out of range");
                assert!(!seen[d], "distance {d} produced twice");
                seen[d] = true;
            }
        }
        assert!(seen.into_iter().all(|s| s), "every distance must be hit");
    }

    #[test]
    fn consecutive_distances_are_grid_adjacent() {
        // The defining Hilbert property: stepping d -> d+1 moves exactly one
        // cell (Manhattan distance 1).
        let order = 2;
        let n = 1u32 << order;
        let mut by_d: Vec<(u32, u32)> = vec![(0, 0); (n * n) as usize];
        for y in 0..n {
            for x in 0..n {
                by_d[xy_to_hilbert_order(x, y, order) as usize] = (x, y);
            }
        }
        for w in by_d.windows(2) {
            let manhattan = w[0].0.abs_diff(w[1].0) + w[0].1.abs_diff(w[1].1);
            assert_eq!(manhattan, 1, "consecutive cells must be adjacent");
        }
    }

    #[test]
    fn degenerate_extent_does_not_panic() {
        let d = hilbert_for_point(5.0, 5.0, 5.0, 5.0, 5.0, 5.0);
        assert!(d > 0);
    }
}
