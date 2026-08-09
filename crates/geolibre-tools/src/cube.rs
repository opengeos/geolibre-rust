//! A raster cube: an ordered stack of co-registered slices along one dimension.
//!
//! `multidimensional_anomaly` (round 16) established the convention the
//! round-18 multidimensional tools follow — one multiband raster whose bands
//! are slices, or a comma-separated list of co-registered rasters whose bands
//! are concatenated in order. This module is that convention factored out, plus
//! the optional per-slice dimension coordinates (dates, depths, wavelengths)
//! that binning and lagged correlation need.
//!
//! `raster_stack::Stack` is the neighbouring abstraction but answers a
//! different question: it groups bands into *reduction* groups and hands back
//! only the valid values at a cell, discarding which slice each came from. A
//! cube has to preserve slice order and slice identity.

use serde_json::Value;
use wbcore::{ToolArgs, ToolError};
use wbraster::Raster;

use crate::common::load_input_raster;
use crate::raster_stack::{check_alignment_refs, parse_input_paths};

/// An ordered stack of co-registered raster slices.
pub(crate) struct Cube {
    rasters: Vec<Raster>,
    /// One entry per slice, in order: `(raster index, band)`.
    slices: Vec<(usize, isize)>,
    pub(crate) rows: usize,
    pub(crate) cols: usize,
    /// Coordinate of each slice along the dimension, when supplied.
    pub(crate) coords: Option<Vec<f64>>,
    /// Name of the dimension, for reporting.
    pub(crate) dimension: String,
}

impl Cube {
    /// Number of slices.
    pub(crate) fn len(&self) -> usize {
        self.slices.len()
    }

    /// The first raster, used as the output geometry/CRS template.
    pub(crate) fn template(&self) -> &Raster {
        &self.rasters[0]
    }

    /// Value at one slice and cell, or `None` when it is no-data.
    pub(crate) fn get(&self, slice: usize, row: usize, col: usize) -> Option<f64> {
        let (ri, band) = self.slices[slice];
        let r = &self.rasters[ri];
        let v = r.get(band, row as isize, col as isize);
        (v != r.nodata && v.is_finite()).then_some(v)
    }

    /// Fills `out` with the full series at one cell, one entry per slice.
    ///
    /// `out` is resized and overwritten, so callers can reuse a single buffer.
    pub(crate) fn series(&self, row: usize, col: usize, out: &mut Vec<Option<f64>>) {
        out.clear();
        out.reserve(self.slices.len());
        for s in 0..self.slices.len() {
            out.push(self.get(s, row, col));
        }
    }

    /// Coordinate of a slice: the supplied value, or its 1-based index when no
    /// coordinates were given.
    pub(crate) fn coord(&self, slice: usize) -> f64 {
        match &self.coords {
            Some(v) => v[slice],
            None => (slice + 1) as f64,
        }
    }
}

/// Loads a cube from a comma-separated `key` parameter.
///
/// `min_slices` is the smallest usable stack for the calling tool.
///
/// `coords_key` and `dimension_key` name the parameters carrying the optional
/// per-slice coordinates and the dimension's name. They are `Option` because a
/// tool that has no use for one must also not *read* it: reading a key the
/// tool does not declare makes it undiscoverable through the registry and CLI,
/// and lets a supplied value fail validation for a parameter the caller was
/// never told about.
pub(crate) fn load_cube(
    args: &ToolArgs,
    key: &str,
    coords_key: Option<&str>,
    dimension_key: Option<&str>,
    min_slices: usize,
) -> Result<Cube, ToolError> {
    let paths = parse_input_paths(args, key)?;
    let rasters: Vec<Raster> = paths
        .iter()
        .map(|p| load_input_raster(p))
        .collect::<Result<_, _>>()?;
    let refs: Vec<&Raster> = rasters.iter().collect();
    check_alignment_refs(&refs)?;

    let (rows, cols) = (rasters[0].rows, rasters[0].cols);
    let mut slices = Vec::new();
    for (i, r) in rasters.iter().enumerate() {
        for b in 0..r.bands {
            slices.push((i, b as isize));
        }
    }
    if slices.len() < min_slices {
        return Err(ToolError::Validation(format!(
            "'{key}' has {} slice(s) (bands across inputs); at least {min_slices} are needed",
            slices.len()
        )));
    }

    let coords = match coords_key {
        Some(k) => parse_coords(args, k, slices.len())?,
        None => None,
    };
    let dimension = dimension_key
        .and_then(|k| args.get(k))
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .unwrap_or("slice")
        .to_string();

    Ok(Cube {
        rasters,
        slices,
        rows,
        cols,
        coords,
        dimension,
    })
}

/// Parses the per-slice dimension coordinates, checking the count and that they
/// increase.
///
/// Binning and lag both assume the slices are in dimension order; accepting an
/// unsorted list would silently produce bins that interleave.
fn parse_coords(
    args: &ToolArgs,
    key: &str,
    n_slices: usize,
) -> Result<Option<Vec<f64>>, ToolError> {
    let Some(s) = args.get(key).and_then(Value::as_str) else {
        return Ok(None);
    };
    if s.trim().is_empty() {
        return Ok(None);
    }
    let vals: Result<Vec<f64>, ToolError> = s
        .split(',')
        .map(str::trim)
        .filter(|t| !t.is_empty())
        .map(|t| {
            t.parse::<f64>()
                .map_err(|_| ToolError::Validation(format!("'{key}' entry '{t}' is not a number")))
        })
        .collect();
    let vals = vals?;
    if vals.len() != n_slices {
        return Err(ToolError::Validation(format!(
            "'{key}' has {} value(s) but the cube has {n_slices} slice(s)",
            vals.len()
        )));
    }
    if vals.windows(2).any(|w| w[1] <= w[0]) {
        return Err(ToolError::Validation(format!(
            "'{key}' must increase strictly; the slices are assumed to be in dimension order"
        )));
    }
    if vals.iter().any(|v| !v.is_finite()) {
        return Err(ToolError::Validation(format!(
            "'{key}' contains a non-finite value"
        )));
    }
    Ok(Some(vals))
}

#[cfg(test)]
pub(crate) mod test_support {
    use wbraster::{CrsInfo, DataType, Raster, RasterConfig};

    /// Builds an in-memory cube raster from per-slice, row-major buffers.
    pub(crate) fn cube_raster(cols: usize, rows: usize, slices: &[Vec<f64>]) -> String {
        cube_raster_typed(cols, rows, slices, DataType::F32)
    }

    /// As [`cube_raster`], with an explicit data type — needed to test that a
    /// tool preserves precision rather than narrowing to F32.
    pub(crate) fn cube_raster_typed(
        cols: usize,
        rows: usize,
        slices: &[Vec<f64>],
        data_type: DataType,
    ) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: slices.len(),
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for (b, slice) in slices.iter().enumerate() {
            for row in 0..rows {
                for col in 0..cols {
                    r.set(b as isize, row as isize, col as isize, slice[row * cols + col])
                        .unwrap();
                }
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }
}

#[cfg(test)]
mod tests {
    use super::test_support::cube_raster;
    use super::*;
    use serde_json::json;

    fn args(v: Value) -> ToolArgs {
        serde_json::from_value(v).unwrap()
    }

    #[test]
    fn loads_slices_in_band_order() {
        let path = cube_raster(2, 1, &[vec![1.0, 2.0], vec![3.0, 4.0], vec![5.0, 6.0]]);
        let cube = load_cube(&args(json!({"input": path})), "input", Some("dv"), Some("dim"), 1).unwrap();
        assert_eq!(cube.len(), 3);
        assert_eq!(cube.get(0, 0, 0), Some(1.0));
        assert_eq!(cube.get(1, 0, 0), Some(3.0));
        assert_eq!(cube.get(2, 0, 1), Some(6.0));
        // No coordinates supplied -> 1-based slice index.
        assert_eq!(cube.coord(0), 1.0);
        assert_eq!(cube.coord(2), 3.0);
        assert_eq!(cube.dimension, "slice");
    }

    /// Several rasters concatenate into one ordered cube.
    #[test]
    fn concatenates_multiple_inputs() {
        let a = cube_raster(1, 1, &[vec![1.0], vec![2.0]]);
        let b = cube_raster(1, 1, &[vec![3.0]]);
        let cube = load_cube(
            &args(json!({"input": format!("{a},{b}")})),
            "input",
            Some("dv"),
            Some("dim"),
            1,
        )
        .unwrap();
        assert_eq!(cube.len(), 3);
        assert_eq!(cube.get(2, 0, 0), Some(3.0));
    }

    #[test]
    fn no_data_reads_as_none() {
        let path = cube_raster(2, 1, &[vec![1.0, -9999.0]]);
        let cube = load_cube(&args(json!({"input": path})), "input", Some("dv"), Some("dim"), 1).unwrap();
        assert_eq!(cube.get(0, 0, 0), Some(1.0));
        assert_eq!(cube.get(0, 0, 1), None);
        let mut s = Vec::new();
        cube.series(0, 1, &mut s);
        assert_eq!(s, vec![None]);
    }

    #[test]
    fn coordinates_are_validated() {
        let path = cube_raster(1, 1, &[vec![1.0], vec![2.0], vec![3.0]]);
        let load = |v: Value| load_cube(&args(v), "input", Some("dv"), Some("dim"), 1);

        let ok = load(json!({"input": path.clone(), "dv": "2000, 2001, 2002"})).unwrap();
        assert_eq!(ok.coord(1), 2001.0);

        // Wrong count.
        assert!(load(json!({"input": path.clone(), "dv": "2000,2001"})).is_err());
        // Not increasing — binning would interleave.
        assert!(load(json!({"input": path.clone(), "dv": "2000,2002,2001"})).is_err());
        // Not a number.
        assert!(load(json!({"input": path.clone(), "dv": "2000,x,2002"})).is_err());
        // Too few slices for the caller.
        assert!(load_cube(&args(json!({"input": path})), "input", Some("dv"), Some("dim"), 9).is_err());
    }
}
