//! Wang & Liu priority-flood depression filling.
//!
//! Ported from the `fill_depressions_wang_and_liu` core in
//! `whitebox-wasm`'s `wbtools_oss` so that `geolibre-tools` stays free of a
//! `wbtools_oss` dependency (keeping it reusable for native/Python bindings and
//! avoiding a dependency cycle). The algorithm seeds a min-heap from the grid
//! edge and floods inward, raising each cell to at least its lowest already
//! processed neighbor (plus an optional flat increment).

use std::cmp::Ordering;
use std::collections::{BinaryHeap, VecDeque};

/// 8-connectivity offsets used by the flood.
const DX: [isize; 8] = [1, 1, 1, 0, -1, -1, -1, 0];
const DY: [isize; 8] = [-1, 0, 1, 1, 1, 0, -1, -1];

#[inline]
fn in_bounds(r: isize, c: isize, rows: usize, cols: usize) -> bool {
    r >= 0 && c >= 0 && (r as usize) < rows && (c as usize) < cols
}

#[inline]
fn idx(r: usize, c: usize, cols: usize) -> usize {
    r * cols + c
}

/// A min-heap node ordered by ascending elevation (the `BinaryHeap` is a
/// max-heap, so the `Ord` impl is reversed).
struct MinNode {
    elev: f64,
    i: usize,
}

impl PartialEq for MinNode {
    fn eq(&self, other: &Self) -> bool {
        self.elev == other.elev
    }
}
impl Eq for MinNode {}
impl PartialOrd for MinNode {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for MinNode {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reversed so the heap pops the lowest elevation first.
        other.elev.total_cmp(&self.elev)
    }
}

/// Fills depressions in a row-major elevation grid using the Wang & Liu
/// priority-flood. `small` is the flat increment added to enforce a monotone
/// descent across flats (use `0.0` to disable). No-data cells are preserved.
pub fn fill_depressions_wang_and_liu(
    input: &[f64],
    rows: usize,
    cols: usize,
    nodata: f64,
    small: f64,
) -> Vec<f64> {
    let background = (i32::MIN + 1) as f64;
    let mut out = vec![background; rows * cols];
    let mut queue = VecDeque::<(isize, isize)>::new();

    // Seed the flood from one ring outside every edge.
    for r in 0..rows as isize {
        queue.push_back((r, -1));
        queue.push_back((r, cols as isize));
    }
    for c in 0..cols as isize {
        queue.push_back((-1, c));
        queue.push_back((rows as isize, c));
    }

    let mut heap = BinaryHeap::<MinNode>::new();
    while let Some((r, c)) = queue.pop_front() {
        for k in 0..8 {
            let rn = r + DY[k];
            let cn = c + DX[k];
            if !in_bounds(rn, cn, rows, cols) {
                continue;
            }
            let ni = idx(rn as usize, cn as usize, cols);
            if out[ni] != background {
                continue;
            }
            let zin = input[ni];
            if zin == nodata {
                out[ni] = nodata;
                queue.push_back((rn, cn));
            } else {
                out[ni] = zin;
                heap.push(MinNode { elev: zin, i: ni });
            }
        }
    }

    while let Some(cell) = heap.pop() {
        let r = cell.i / cols;
        let c = cell.i % cols;
        let z = out[cell.i];
        for k in 0..8 {
            let rn = r as isize + DY[k];
            let cn = c as isize + DX[k];
            if !in_bounds(rn, cn, rows, cols) {
                continue;
            }
            let ni = idx(rn as usize, cn as usize, cols);
            if out[ni] != background {
                continue;
            }
            let mut zn = input[ni];
            if zn != nodata {
                if zn < z + small {
                    zn = z + small;
                }
                out[ni] = zn;
                heap.push(MinNode { elev: zn, i: ni });
            } else {
                out[ni] = nodata;
            }
        }
    }

    for v in &mut out {
        if *v == background {
            *v = nodata;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fills_a_single_pit() {
        // 3x3 bowl: a low center surrounded by a rim of 10s gets raised to 10.
        let rows = 3;
        let cols = 3;
        let nodata = -9999.0;
        let input = vec![10.0, 10.0, 10.0, 10.0, 1.0, 10.0, 10.0, 10.0, 10.0];
        let out = fill_depressions_wang_and_liu(&input, rows, cols, nodata, 0.0);
        assert_eq!(out[idx(1, 1, cols)], 10.0);
        // Rim is unchanged.
        assert_eq!(out[idx(0, 0, cols)], 10.0);
    }

    #[test]
    fn preserves_nodata() {
        let rows = 2;
        let cols = 2;
        let nodata = -9999.0;
        let input = vec![5.0, nodata, 3.0, 4.0];
        let out = fill_depressions_wang_and_liu(&input, rows, cols, nodata, 0.0);
        assert_eq!(out[1], nodata);
    }

    #[test]
    fn leaves_monotone_surface_unchanged() {
        let rows = 1;
        let cols = 4;
        let nodata = -9999.0;
        let input = vec![4.0, 3.0, 2.0, 1.0];
        let out = fill_depressions_wang_and_liu(&input, rows, cols, nodata, 0.0);
        assert_eq!(out, input);
    }
}
