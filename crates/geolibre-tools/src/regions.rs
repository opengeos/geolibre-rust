//! Connected-component labeling and region properties, ported from the parts of
//! `scipy.ndimage` and `skimage.measure` used by the `lidar` Python package.
//!
//! Two pieces are reproduced faithfully enough to match the Python outputs:
//!
//! * [`region_group`] mirrors `lidar`'s `regionGroup`: label nonzero cells with
//!   4-connectivity (`scipy.ndimage.label`'s default cross structuring element),
//!   drop components with `<= min_size` pixels, then relabel.
//! * [`region_props`] mirrors the `skimage.measure.regionprops` attributes the
//!   package consumes: area, bounding box, first coordinate, intensity
//!   statistics, and the inertia-tensor shape descriptors (major/minor axis,
//!   eccentricity, orientation) plus `skimage`'s 4-connectivity perimeter.

/// 4-connectivity offsets (row, col): up, down, left, right.
const NEIGHBORS_4: [(isize, isize); 4] = [(-1, 0), (1, 0), (0, -1), (0, 1)];

/// Labels connected components of a boolean mask using 4-connectivity.
///
/// Returns `(labels, count)` where `labels[r*cols + c]` is the 1-based component
/// id (0 = background) and `count` is the number of components. This matches
/// `scipy.ndimage.label` with its default (cross) structuring element.
pub fn label_4(mask: &[bool], rows: usize, cols: usize) -> (Vec<u32>, usize) {
    let mut labels = vec![0_u32; rows * cols];
    let mut next: u32 = 0;
    let mut stack: Vec<(usize, usize)> = Vec::new();

    for start_r in 0..rows {
        for start_c in 0..cols {
            let start = start_r * cols + start_c;
            if !mask[start] || labels[start] != 0 {
                continue;
            }
            next += 1;
            labels[start] = next;
            stack.push((start_r, start_c));
            while let Some((r, c)) = stack.pop() {
                for (dr, dc) in NEIGHBORS_4 {
                    let nr = r as isize + dr;
                    let nc = c as isize + dc;
                    if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                        continue;
                    }
                    let ni = nr as usize * cols + nc as usize;
                    if mask[ni] && labels[ni] == 0 {
                        labels[ni] = next;
                        stack.push((nr as usize, nc as usize));
                    }
                }
            }
        }
    }
    (labels, next as usize)
}

/// Ports `lidar`'s `regionGroup`: treat every nonzero cell as foreground, label
/// with 4-connectivity, discard components with `<= min_size` pixels (strictly
/// greater is kept, matching `sizes > min_size`), then relabel the survivors.
///
/// `values` is the working buffer; cells equal to `nodata` are treated as
/// background (0). Returns `(labels, count)` of the relabeled image.
pub fn region_group(values: &[f64], rows: usize, cols: usize, min_size: u64, nodata: f64) -> (Vec<u32>, usize) {
    let mask: Vec<bool> = values
        .iter()
        .map(|&v| v != 0.0 && v != nodata)
        .collect();
    let (labels, count) = label_4(&mask, rows, cols);

    // Count pixels per component.
    let mut sizes = vec![0_u64; count + 1];
    for &l in &labels {
        sizes[l as usize] += 1;
    }

    // Keep components strictly larger than min_size; background (0) is dropped.
    let keep: Vec<bool> = sizes.iter().enumerate().map(|(i, &s)| i != 0 && s > min_size).collect();
    let cleaned: Vec<bool> = labels.iter().map(|&l| keep[l as usize]).collect();
    label_4(&cleaned, rows, cols)
}

/// Region attributes mirroring the `skimage.measure.regionprops` fields used by
/// `lidar`. Geometry is in pixel units; callers scale by spatial resolution.
#[derive(Debug, Clone)]
pub struct RegionProps {
    /// 1-based component id (the `label`).
    pub label: u32,
    /// Pixel count (`area`).
    pub area: u64,
    /// Bounding box `(min_row, min_col, max_row, max_col)`, max exclusive.
    pub bbox: (usize, usize, usize, usize),
    /// First pixel in row-major order `(row, col)` (`coords[0]`).
    pub first_coord: (usize, usize),
    /// Minimum intensity over the region.
    pub min_intensity: f64,
    /// Maximum intensity over the region.
    pub max_intensity: f64,
    /// Sum of intensity over the region (`intensity_image` sums to this).
    pub sum_intensity: f64,
    /// Major axis length in pixels (`4 * sqrt(l1)` from the inertia tensor).
    pub major_axis_length: f64,
    /// Minor axis length in pixels (`4 * sqrt(l2)`).
    pub minor_axis_length: f64,
    /// Eccentricity, `sqrt(1 - l2/l1)`.
    pub eccentricity: f64,
    /// Orientation in radians (`skimage` convention).
    pub orientation: f64,
    /// Perimeter in pixels (`skimage` 4-connectivity weighting).
    pub perimeter: f64,
    /// Filled fraction of the bounding box (`extent` = area / bbox area).
    pub extent: f64,
}

/// Computes [`RegionProps`] for every labeled component (1..=count).
///
/// `labels` and `intensity` are row-major buffers of length `rows*cols`.
/// Background (`label == 0`) is ignored. Results are ordered by ascending label,
/// matching `skimage.measure.regionprops`.
pub fn region_props(
    labels: &[u32],
    count: usize,
    intensity: &[f64],
    rows: usize,
    cols: usize,
) -> Vec<RegionProps> {
    if count == 0 {
        return Vec::new();
    }

    // Per-label accumulators.
    let mut area = vec![0_u64; count + 1];
    let mut min_int = vec![f64::INFINITY; count + 1];
    let mut max_int = vec![f64::NEG_INFINITY; count + 1];
    let mut sum_int = vec![0.0_f64; count + 1];
    let mut min_row = vec![usize::MAX; count + 1];
    let mut min_col = vec![usize::MAX; count + 1];
    let mut max_row = vec![0_usize; count + 1];
    let mut max_col = vec![0_usize; count + 1];
    let mut first = vec![None::<(usize, usize)>; count + 1];
    // Moment sums for centroid and central moments.
    let mut sum_r = vec![0.0_f64; count + 1];
    let mut sum_c = vec![0.0_f64; count + 1];

    for row in 0..rows {
        for col in 0..cols {
            let l = labels[row * cols + col] as usize;
            if l == 0 {
                continue;
            }
            area[l] += 1;
            let v = intensity[row * cols + col];
            if v < min_int[l] {
                min_int[l] = v;
            }
            if v > max_int[l] {
                max_int[l] = v;
            }
            sum_int[l] += v;
            if row < min_row[l] {
                min_row[l] = row;
            }
            if col < min_col[l] {
                min_col[l] = col;
            }
            if row + 1 > max_row[l] {
                max_row[l] = row + 1;
            }
            if col + 1 > max_col[l] {
                max_col[l] = col + 1;
            }
            if first[l].is_none() {
                first[l] = Some((row, col)); // row-major scan => coords[0]
            }
            sum_r[l] += row as f64;
            sum_c[l] += col as f64;
        }
    }

    let mut out = Vec::with_capacity(count);
    for l in 1..=count {
        if area[l] == 0 {
            continue;
        }
        let n = area[l] as f64;
        let cr = sum_r[l] / n;
        let cc = sum_c[l] / n;
        let bbox = (min_row[l], min_col[l], max_row[l], max_col[l]);

        // Second central moments over the binary region (skimage inertia tensor).
        let (mu20, mu02, mu11) = central_moments(labels, l as u32, cr, cc, bbox, cols);
        let a = mu20 / n; // inertia_tensor[0,0]
        let c = mu02 / n; // inertia_tensor[1,1]
        let b = -mu11 / n; // inertia_tensor[0,1] = inertia_tensor[1,0]

        let common = (4.0 * b * b + (a - c) * (a - c)).max(0.0).sqrt();
        let l1 = (a + c) / 2.0 + common / 2.0;
        let l2 = ((a + c) / 2.0 - common / 2.0).max(0.0);
        let major = 4.0 * l1.max(0.0).sqrt();
        let minor = 4.0 * l2.sqrt();
        let eccentricity = if l1 > 0.0 {
            (1.0 - l2 / l1).max(0.0).sqrt()
        } else {
            0.0
        };
        let orientation = if a - c == 0.0 {
            if b < 0.0 {
                -std::f64::consts::FRAC_PI_4
            } else {
                std::f64::consts::FRAC_PI_4
            }
        } else {
            0.5 * (-2.0 * b).atan2(c - a)
        };

        let bbox_area = ((bbox.2 - bbox.0) * (bbox.3 - bbox.1)) as f64;
        let extent = if bbox_area > 0.0 { n / bbox_area } else { 0.0 };
        let perimeter = perimeter_4(labels, l as u32, bbox, cols);

        out.push(RegionProps {
            label: l as u32,
            area: area[l],
            bbox,
            first_coord: first[l].unwrap(),
            min_intensity: min_int[l],
            max_intensity: max_int[l],
            sum_intensity: sum_int[l],
            major_axis_length: major,
            minor_axis_length: minor,
            eccentricity,
            orientation,
            perimeter,
            extent,
        });
    }
    out
}

/// Second-order central moments `(mu20, mu02, mu11)` of a labeled region's
/// binary mask, taken about the centroid `(cr, cc)`.
fn central_moments(
    labels: &[u32],
    label: u32,
    cr: f64,
    cc: f64,
    bbox: (usize, usize, usize, usize),
    cols: usize,
) -> (f64, f64, f64) {
    let (min_row, min_col, max_row, max_col) = bbox;
    let mut mu20 = 0.0;
    let mut mu02 = 0.0;
    let mut mu11 = 0.0;
    for row in min_row..max_row {
        for col in min_col..max_col {
            if labels[row * cols + col] != label {
                continue;
            }
            let dr = row as f64 - cr;
            let dc = col as f64 - cc;
            mu20 += dr * dr;
            mu02 += dc * dc;
            mu11 += dr * dc;
        }
    }
    (mu20, mu02, mu11)
}

/// `skimage.measure.perimeter` with 4-connectivity, evaluated over a region's
/// bounding box.
///
/// The boundary is `mask - erosion(mask)` (border pixels), and each border pixel
/// is weighted by a kernel-coded count of its orthogonal/diagonal border
/// neighbors. Only border pixels (kernel center = 1) carry nonzero weight, so we
/// score them directly instead of materializing the full convolution.
fn perimeter_4(labels: &[u32], label: u32, bbox: (usize, usize, usize, usize), cols: usize) -> f64 {
    let (min_row, min_col, max_row, max_col) = bbox;
    let in_region = |r: isize, c: isize| -> bool {
        if r < 0 || c < 0 {
            return false;
        }
        let (r, c) = (r as usize, c as usize);
        r < (labels.len() / cols) && c < cols && labels[r * cols + c] == label
    };

    // A border pixel is in the region but has at least one 4-neighbor outside
    // it (binary_erosion with a cross structuring element and border_value=0).
    let is_border = |r: isize, c: isize| -> bool {
        if !in_region(r, c) {
            return false;
        }
        for (dr, dc) in NEIGHBORS_4 {
            if !in_region(r + dr, c + dc) {
                return true;
            }
        }
        false
    };

    let sqrt2 = std::f64::consts::SQRT_2;
    let mut total = 0.0;
    for row in min_row..max_row {
        for col in min_col..max_col {
            let r = row as isize;
            let c = col as isize;
            if !is_border(r, c) {
                continue;
            }
            // Count orthogonal and diagonal border neighbors.
            let mut orth = 0;
            for (dr, dc) in NEIGHBORS_4 {
                if is_border(r + dr, c + dc) {
                    orth += 1;
                }
            }
            let mut diag = 0;
            for (dr, dc) in [(-1, -1), (-1, 1), (1, -1), (1, 1)] {
                if is_border(r + dr, c + dc) {
                    diag += 1;
                }
            }
            // Convolution code = 1 (center) + 2*orth + 10*diag, mapped to the
            // skimage perimeter weights.
            let code = 1 + 2 * orth + 10 * diag;
            total += match code {
                5 | 7 | 15 | 17 | 25 | 27 => 1.0,
                21 | 33 => sqrt2,
                13 | 23 => (1.0 + sqrt2) / 2.0,
                _ => 0.0,
            };
        }
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labels_two_separate_blocks() {
        // 4x4 with two 2x2 blocks in opposite corners.
        let rows = 4;
        let cols = 4;
        let mut mask = vec![false; rows * cols];
        for (r, c) in [(0, 0), (0, 1), (1, 0), (1, 1), (2, 2), (2, 3), (3, 2), (3, 3)] {
            mask[r * cols + c] = true;
        }
        let (labels, count) = label_4(&mask, rows, cols);
        assert_eq!(count, 2);
        assert_eq!(labels[0], 1);
        assert_eq!(labels[2 * cols + 2], 2);
    }

    #[test]
    fn region_group_drops_small_components() {
        // One large (>min_size) component and one single-pixel component.
        let rows = 3;
        let cols = 4;
        let mut v = vec![0.0; rows * cols];
        for (r, c) in [(0, 0), (0, 1), (1, 0), (1, 1)] {
            v[r * cols + c] = 5.0;
        }
        v[2 * cols + 3] = 9.0; // lone pixel
        let (labels, count) = region_group(&v, rows, cols, 1, 0.0);
        assert_eq!(count, 1); // 4-pixel block survives min_size=1, lone pixel dropped
        assert_eq!(labels[2 * cols + 3], 0);
        assert_eq!(labels[0], 1);
    }

    #[test]
    fn props_area_bbox_and_intensity() {
        let rows = 3;
        let cols = 3;
        let mut labels = vec![0_u32; rows * cols];
        // A 2x2 block in the top-left.
        for (r, c) in [(0, 0), (0, 1), (1, 0), (1, 1)] {
            labels[r * cols + c] = 1;
        }
        let intensity = vec![1.0, 2.0, 0.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0];
        let props = region_props(&labels, 1, &intensity, rows, cols);
        assert_eq!(props.len(), 1);
        let p = &props[0];
        assert_eq!(p.area, 4);
        assert_eq!(p.bbox, (0, 0, 2, 2));
        assert_eq!(p.first_coord, (0, 0));
        assert_eq!(p.min_intensity, 1.0);
        assert_eq!(p.max_intensity, 4.0);
        assert_eq!(p.sum_intensity, 10.0);
        assert_eq!(p.extent, 1.0);
    }
}
