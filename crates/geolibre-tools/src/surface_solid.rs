//! Shared surface-sampling and shell-building helpers for the round-17
//! surface-to-solid tools (`extrude_between`, `fence_diagram`).
//!
//! Both tools do the same two things: sample a raster surface at arbitrary
//! (x, y) positions, and stitch sampled points into triangle strips. Keeping
//! that here means the watertightness convention (`buffer_3d`'s triangle-per-
//! part `MultiPolygon`) is applied identically by both.

use wbcore::ToolError;
use wbraster::Raster;
use wbvector::Coord;

use crate::inside_3d::Tri;

/// Bilinear sample of a raster at a map position, or `None` off-grid / no-data.
///
/// Bilinear rather than nearest because both callers build *geometry* from the
/// result: nearest-neighbour sampling produces visible stair-stepping along a
/// section trace or a polygon boundary that no downstream smoothing removes.
pub(crate) fn sample_bilinear(raster: &Raster, band: isize, x: f64, y: f64) -> Option<f64> {
    let cols = raster.cols as f64;
    let rows = raster.rows as f64;
    // Continuous cell coordinates with the cell centre at (i + 0.5).
    let fx = (x - raster.x_min) / raster.cell_size_x - 0.5;
    // Rows count down from the top, so y increases as the row index decreases.
    let fy = (rows - 1.0) - ((y - raster.y_min) / raster.cell_size_y - 0.5);
    if !fx.is_finite() || !fy.is_finite() {
        return None;
    }
    let x0 = fx.floor();
    let y0 = fy.floor();
    let tx = fx - x0;
    let ty = fy - y0;

    let mut acc = 0.0;
    let mut wsum = 0.0;
    for (dx, dy, w) in [
        (0.0, 0.0, (1.0 - tx) * (1.0 - ty)),
        (1.0, 0.0, tx * (1.0 - ty)),
        (0.0, 1.0, (1.0 - tx) * ty),
        (1.0, 1.0, tx * ty),
    ] {
        let cx = x0 + dx;
        let cy = y0 + dy;
        if cx < 0.0 || cy < 0.0 || cx >= cols || cy >= rows || w <= 0.0 {
            continue;
        }
        let v = raster.get(band, cy as isize, cx as isize);
        if v == raster.nodata || !v.is_finite() {
            continue;
        }
        acc += w * v;
        wsum += w;
    }
    // Partial coverage still yields a value (renormalised), so a surface edge
    // does not punch holes in an otherwise valid shell.
    (wsum > 0.0).then(|| acc / wsum)
}

/// Stitches two parallel point runs into a triangle strip.
///
/// `lower[i]` and `upper[i]` must correspond. Emits two triangles per quad,
/// wound so the strip's normals face consistently — the precondition for the
/// signed-tetrahedron volume to be meaningful.
pub(crate) fn strip(lower: &[[f64; 3]], upper: &[[f64; 3]], flip: bool) -> Vec<Tri> {
    let n = lower.len().min(upper.len());
    let mut out = Vec::with_capacity(n.saturating_sub(1) * 2);
    for i in 0..n.saturating_sub(1) {
        let (a, b, c, d) = (lower[i], lower[i + 1], upper[i + 1], upper[i]);
        if flip {
            out.push([a, c, b]);
            out.push([a, d, c]);
        } else {
            out.push([a, b, c]);
            out.push([a, c, d]);
        }
    }
    out
}

/// Densifies a coordinate ring/line so no segment exceeds `spacing`.
///
/// Sampling only at the input vertices would let a long straight edge cross a
/// whole valley with no intermediate elevation, which is the difference
/// between a shell that follows the terrain and one that tunnels through it.
pub(crate) fn densify(coords: &[Coord], spacing: f64) -> Vec<(f64, f64)> {
    let mut out: Vec<(f64, f64)> = Vec::new();
    if coords.is_empty() {
        return out;
    }
    out.push((coords[0].x, coords[0].y));
    for w in coords.windows(2) {
        let (x0, y0) = (w[0].x, w[0].y);
        let (x1, y1) = (w[1].x, w[1].y);
        let len = ((x1 - x0).powi(2) + (y1 - y0).powi(2)).sqrt();
        if len <= 0.0 {
            continue;
        }
        let steps = if spacing > 0.0 {
            (len / spacing).ceil().max(1.0) as usize
        } else {
            1
        };
        for k in 1..=steps {
            let t = k as f64 / steps as f64;
            out.push((x0 + t * (x1 - x0), y0 + t * (y1 - y0)));
        }
    }
    out
}

/// Default sampling spacing: the finest cell size among the given rasters.
pub(crate) fn default_spacing(rasters: &[&Raster]) -> Result<f64, ToolError> {
    let s = rasters
        .iter()
        .map(|r| r.cell_size_x.abs().min(r.cell_size_y.abs()))
        .filter(|v| *v > 0.0)
        .fold(f64::INFINITY, f64::min);
    if !s.is_finite() || s <= 0.0 {
        return Err(ToolError::Execution(
            "could not derive a sampling spacing from the input surfaces".to_string(),
        ));
    }
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbraster::{CrsInfo, DataType, RasterConfig};

    fn ramp() -> Raster {
        // 2x2 cells of size 10, values 0, 10 / 20, 30 (row 0 is the top).
        let mut r = Raster::new(RasterConfig {
            cols: 2,
            rows: 2,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 10.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        r.set(0, 0, 0, 0.0).unwrap();
        r.set(0, 0, 1, 10.0).unwrap();
        r.set(0, 1, 0, 20.0).unwrap();
        r.set(0, 1, 1, 30.0).unwrap();
        r
    }

    #[test]
    fn sampling_at_a_cell_centre_returns_that_cell() {
        let r = ramp();
        // Cell (row 1, col 0) is the lower-left, centre (5, 5), value 20.
        assert!((sample_bilinear(&r, 0, 5.0, 5.0).unwrap() - 20.0).abs() < 1e-9);
        // Cell (row 0, col 1) is the upper-right, centre (15, 15), value 10.
        assert!((sample_bilinear(&r, 0, 15.0, 15.0).unwrap() - 10.0).abs() < 1e-9);
    }

    #[test]
    fn sampling_between_centres_interpolates() {
        let r = ramp();
        // Midway between the two lower cells (20 and 30).
        let v = sample_bilinear(&r, 0, 10.0, 5.0).unwrap();
        assert!((v - 25.0).abs() < 1e-9, "got {v}");
    }

    #[test]
    fn sampling_off_the_grid_returns_none() {
        let r = ramp();
        assert!(sample_bilinear(&r, 0, -500.0, -500.0).is_none());
    }

    #[test]
    fn densify_inserts_points_without_exceeding_the_spacing() {
        let line = vec![Coord::xy(0.0, 0.0), Coord::xy(10.0, 0.0)];
        let pts = densify(&line, 2.0);
        assert_eq!(pts.len(), 6); // 0, 2, 4, 6, 8, 10
        for w in pts.windows(2) {
            let d = ((w[1].0 - w[0].0).powi(2) + (w[1].1 - w[0].1).powi(2)).sqrt();
            assert!(d <= 2.0 + 1e-9, "segment of {d} exceeds the spacing");
        }
    }

    #[test]
    fn a_strip_between_two_runs_has_two_triangles_per_quad() {
        let lower = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let upper = vec![[0.0, 0.0, 5.0], [1.0, 0.0, 5.0], [2.0, 0.0, 5.0]];
        let tris = strip(&lower, &upper, false);
        assert_eq!(tris.len(), 4);
    }
}
