//! Machinery shared by the SAR tools: unit handling, automatic thresholding,
//! polygon masking, and connected-region extraction with an area filter.
//!
//! The round-17 SAR tools (`multilook`, `sar_coherence`,
//! `apply_radiometric_calibration`, `convert_sar_units`, `compute_sar_indices`)
//! each carried their own copy of the amplitude/intensity/dB conversion. This
//! module is the single definition, extended with the pieces the ocean and
//! water tools need.

use std::collections::HashMap;

use serde_json::{Map, Value};
use wbcore::ToolError;
use wbraster::Raster;
use wbvector::{Coord, Geometry, Layer};

use crate::polygonize::{polygonize_to_geojson, PolygonizeParams};
use crate::vector_common::geometry_contains_point;

/// The unit a SAR raster's samples are in.
///
/// `Dn` and `Amplitude` are field quantities (power is their square);
/// `Intensity`/`Linear` are already power; `Db` is `10*log10(power)`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum SarUnits {
    Dn,
    Amplitude,
    Intensity,
    Db,
}

impl SarUnits {
    pub(crate) fn parse(s: &str) -> Result<Self, ToolError> {
        match s.trim().to_ascii_lowercase().as_str() {
            "" | "intensity" | "linear" | "power" => Ok(SarUnits::Intensity),
            "dn" => Ok(SarUnits::Dn),
            "amplitude" => Ok(SarUnits::Amplitude),
            "db" | "decibel" => Ok(SarUnits::Db),
            other => Err(ToolError::Validation(format!(
                "SAR units must be one of dn|amplitude|intensity|db, got '{other}'"
            ))),
        }
    }

    pub(crate) fn label(self) -> &'static str {
        match self {
            SarUnits::Dn => "dn",
            SarUnits::Amplitude => "amplitude",
            SarUnits::Intensity => "intensity",
            SarUnits::Db => "db",
        }
    }

    /// Converts a sample to linear power, or `None` if the value is not
    /// physically representable (a negative intensity, a non-finite sample).
    pub(crate) fn to_power(self, v: f64) -> Option<f64> {
        if !v.is_finite() {
            return None;
        }
        match self {
            // DN and amplitude are field quantities: power is their square.
            SarUnits::Dn | SarUnits::Amplitude => Some(v * v),
            SarUnits::Intensity => (v >= 0.0).then_some(v),
            SarUnits::Db => Some(10f64.powf(v / 10.0)),
        }
    }
}

/// Linear power to decibels. Non-positive power has no dB value, so callers
/// must map `None` onto no-data rather than emitting `-inf`.
pub(crate) fn power_to_db(p: f64) -> Option<f64> {
    (p > 0.0 && p.is_finite()).then(|| 10.0 * p.log10())
}

/// Otsu's method: the threshold that maximises between-class variance of a
/// histogram of `values`, returned in the units of `values`.
///
/// Used to split a bimodal backscatter histogram (dark water against bright
/// land, or a dark slick against the surrounding sea) without a hand-tuned
/// constant. Returns `None` if there is no spread to split.
#[allow(clippy::needless_range_loop)] // the bin index is the histogram's x axis
pub(crate) fn otsu_threshold(values: &[f64], bins: usize) -> Option<f64> {
    let bins = bins.max(2);
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &v in values {
        if v.is_finite() {
            lo = lo.min(v);
            hi = hi.max(v);
        }
    }
    if !lo.is_finite() || !hi.is_finite() || hi <= lo {
        return None;
    }

    let width = (hi - lo) / bins as f64;
    let mut hist = vec![0usize; bins];
    let mut total = 0usize;
    for &v in values {
        if !v.is_finite() {
            continue;
        }
        let idx = (((v - lo) / width) as usize).min(bins - 1);
        hist[idx] += 1;
        total += 1;
    }
    if total == 0 {
        return None;
    }

    // Bin centres, so the returned threshold is a value not an edge.
    let centre = |i: usize| lo + (i as f64 + 0.5) * width;
    let sum_all: f64 = hist
        .iter()
        .enumerate()
        .map(|(i, &n)| centre(i) * n as f64)
        .sum();

    let (mut w_bg, mut sum_bg) = (0.0f64, 0.0f64);
    let (mut best_var, mut best_idx) = (f64::NEG_INFINITY, 0usize);
    for i in 0..bins {
        w_bg += hist[i] as f64;
        if w_bg == 0.0 {
            continue;
        }
        let w_fg = total as f64 - w_bg;
        if w_fg == 0.0 {
            break;
        }
        sum_bg += centre(i) * hist[i] as f64;
        let mean_bg = sum_bg / w_bg;
        let mean_fg = (sum_all - sum_bg) / w_fg;
        let var = w_bg * w_fg * (mean_bg - mean_fg) * (mean_bg - mean_fg);
        if var > best_var {
            best_var = var;
            best_idx = i;
        }
    }
    (best_var > f64::NEG_INFINITY).then(|| centre(best_idx) + 0.5 * width)
}

/// Which side of a mask polygon set is the analysis area.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum MaskSide {
    /// Polygons delineate land; analyse everything outside them.
    LandPolygon,
    /// Polygons delineate water; analyse everything inside them.
    WaterPolygon,
}

impl MaskSide {
    pub(crate) fn parse(s: &str) -> Result<Self, ToolError> {
        match s.trim().to_ascii_lowercase().as_str() {
            "" | "land_polygon" | "land" => Ok(MaskSide::LandPolygon),
            "water_polygon" | "water" => Ok(MaskSide::WaterPolygon),
            other => Err(ToolError::Validation(format!(
                "mask side must be 'land_polygon' or 'water_polygon', got '{other}'"
            ))),
        }
    }
}

/// Rasterises a polygon mask layer onto `raster`'s grid.
///
/// Returns a per-cell "analyse this cell" flag: for [`MaskSide::LandPolygon`] a
/// cell is analysed when it falls *outside* every polygon, for
/// [`MaskSide::WaterPolygon`] when it falls *inside* at least one. Point-in-
/// polygon runs against cell centres via [`geometry_contains_point`], the same
/// predicate the rest of the crate uses.
pub(crate) fn rasterize_mask(raster: &Raster, layer: &Layer, side: MaskSide) -> Vec<bool> {
    let (rows, cols) = (raster.rows, raster.cols);
    let geoms: Vec<&Geometry> = layer.iter().filter_map(|f| f.geometry.as_ref()).collect();
    let y_max = raster.y_min + rows as f64 * raster.cell_size_y;

    let mut keep = vec![side == MaskSide::LandPolygon; rows * cols];
    if geoms.is_empty() {
        return keep;
    }
    // Envelope per geometry, computed once. Without it every cell walks every
    // ring of every polygon, and a coastline mask of thousands of polygons over
    // a full scene dominates the runtime of all three tools that mask.
    let boxes: Vec<(f64, f64, f64, f64)> = geoms.iter().map(|g| envelope(g)).collect();

    for r in 0..rows {
        let y = y_max - (r as f64 + 0.5) * raster.cell_size_y;
        for c in 0..cols {
            let x = raster.x_min + (c as f64 + 0.5) * raster.cell_size_x;
            let inside = geoms.iter().zip(&boxes).any(|(g, bb)| {
                x >= bb.0 && x <= bb.2 && y >= bb.1 && y <= bb.3
                    && geometry_contains_point(g, x, y)
            });
            keep[r * cols + c] = match side {
                MaskSide::LandPolygon => !inside,
                MaskSide::WaterPolygon => inside,
            };
        }
    }
    keep
}

/// Axis-aligned envelope `(min_x, min_y, max_x, max_y)` of a geometry.
///
/// An empty geometry gets an inverted box, which the prefilter then rejects for
/// every cell — the containment test would say the same, only far more slowly.
/// The variants covered here must therefore track `geometry_contains_point`
/// exactly: a geometry the containment test *can* accept but the envelope skips
/// would be prefiltered away everywhere, masking the whole raster out.
fn envelope(geom: &Geometry) -> (f64, f64, f64, f64) {
    let mut bb = (f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    accumulate_envelope(geom, &mut bb);
    bb
}

fn accumulate_envelope(geom: &Geometry, bb: &mut (f64, f64, f64, f64)) {
    let mut visit = |c: &Coord| {
        bb.0 = bb.0.min(c.x);
        bb.1 = bb.1.min(c.y);
        bb.2 = bb.2.max(c.x);
        bb.3 = bb.3.max(c.y);
    };
    match geom {
        Geometry::Polygon { exterior, .. } => exterior.coords().iter().for_each(&mut visit),
        Geometry::MultiPolygon(parts) => {
            for (ext, _) in parts {
                ext.coords().iter().for_each(&mut visit);
            }
        }
        Geometry::GeometryCollection(members) => {
            for member in members {
                accumulate_envelope(member, bb);
            }
        }
        _ => {}
    }
}

/// One 4-connected run of flagged cells.
pub(crate) struct Region {
    pub(crate) cells: Vec<usize>,
    /// Sum of the source values over the region's cells, for a mean.
    pub(crate) value_sum: f64,
}

impl Region {
    pub(crate) fn mean(&self) -> f64 {
        if self.cells.is_empty() {
            0.0
        } else {
            self.value_sum / self.cells.len() as f64
        }
    }
}

/// Groups flagged cells into 4-connected regions, keeping only those with at
/// least `min_cells` cells.
///
/// `values` supplies the per-cell quantity summed into [`Region::value_sum`]
/// (typically the backscatter that triggered the flag).
pub(crate) fn connected_regions(
    flag: &[bool],
    values: &[f64],
    rows: usize,
    cols: usize,
    min_cells: usize,
) -> Vec<Region> {
    let mut seen = vec![false; rows * cols];
    let mut out = Vec::new();
    let mut stack: Vec<usize> = Vec::new();

    for start in 0..rows * cols {
        if !flag[start] || seen[start] {
            continue;
        }
        seen[start] = true;
        stack.clear();
        stack.push(start);
        let mut cells = Vec::new();
        let mut value_sum = 0.0;
        while let Some(i) = stack.pop() {
            cells.push(i);
            value_sum += values[i];
            let (r, c) = (i / cols, i % cols);
            // 4-connectivity, matching `polygonize`'s component rule so the
            // regions and the traced rings agree.
            let push = |rr: usize, cc: usize, stack: &mut Vec<usize>, seen: &mut Vec<bool>| {
                let j = rr * cols + cc;
                if flag[j] && !seen[j] {
                    seen[j] = true;
                    stack.push(j);
                }
            };
            if r > 0 {
                push(r - 1, c, &mut stack, &mut seen);
            }
            if r + 1 < rows {
                push(r + 1, c, &mut stack, &mut seen);
            }
            if c > 0 {
                push(r, c - 1, &mut stack, &mut seen);
            }
            if c + 1 < cols {
                push(r, c + 1, &mut stack, &mut seen);
            }
        }
        if cells.len() >= min_cells {
            out.push(Region { cells, value_sum });
        }
    }
    out
}

/// Polygonizes a set of regions into one GeoJSON geometry per region, in the
/// raster's CRS. Region *i* becomes label `i + 1`.
///
/// Returns the geometries in region order; a region whose ring trace collapses
/// is dropped, so the result can be shorter than `regions`.
pub(crate) fn regions_to_geometries(
    raster: &Raster,
    regions: &[Region],
    rows: usize,
    cols: usize,
) -> Result<Vec<(usize, Geometry)>, ToolError> {
    let mut labels = vec![0.0f64; rows * cols];
    for (i, reg) in regions.iter().enumerate() {
        for &cell in &reg.cells {
            labels[cell] = (i + 1) as f64;
        }
    }
    let props: HashMap<i64, Map<String, Value>> = HashMap::new();
    let geojson = polygonize_to_geojson(&PolygonizeParams {
        labels: &labels,
        rows,
        cols,
        x_min: raster.x_min,
        y_max: raster.y_min + rows as f64 * raster.cell_size_y,
        cell_size_x: raster.cell_size_x,
        cell_size_y: raster.cell_size_y,
        epsg: raster.crs.epsg,
        props_by_id: &props,
    });
    parse_labelled_polygons(&geojson)
}

/// Parses `polygonize_to_geojson` output into `(region index, geometry)` pairs.
fn parse_labelled_polygons(geojson: &str) -> Result<Vec<(usize, Geometry)>, ToolError> {
    let v: Value = serde_json::from_str(geojson)
        .map_err(|e| ToolError::Execution(format!("polygonize produced invalid GeoJSON: {e}")))?;
    let feats = v
        .get("features")
        .and_then(Value::as_array)
        .ok_or_else(|| ToolError::Execution("polygonize output has no features".to_string()))?;

    let mut out = Vec::new();
    for f in feats {
        let id = f
            .get("properties")
            .and_then(|p| p.get("id"))
            .and_then(Value::as_i64)
            .unwrap_or(0);
        if id <= 0 {
            continue;
        }
        let Some(geom) = f.get("geometry") else {
            continue;
        };
        if let Some(g) = crate::geojson_geom::geometry_from_json(geom) {
            out.push((id as usize - 1, g));
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbraster::{CrsInfo, DataType, RasterConfig};
    use wbvector::{GeometryType, Ring};

    #[test]
    fn unit_conversions_agree_on_one_power() {
        // Amplitude 10, intensity 100 and 20 dB are the same power.
        let a = SarUnits::Amplitude.to_power(10.0).unwrap();
        let i = SarUnits::Intensity.to_power(100.0).unwrap();
        let d = SarUnits::Db.to_power(20.0).unwrap();
        assert!((a - 100.0).abs() < 1e-9);
        assert!((i - 100.0).abs() < 1e-9);
        assert!((d - 100.0).abs() < 1e-9);
    }

    #[test]
    fn negative_intensity_has_no_power() {
        assert!(SarUnits::Intensity.to_power(-1.0).is_none());
        // ...but a negative dB value is perfectly ordinary (power below 1).
        assert!(SarUnits::Db.to_power(-20.0).unwrap() > 0.0);
    }

    #[test]
    fn zero_power_has_no_db() {
        assert!(power_to_db(0.0).is_none());
        assert!(power_to_db(-1.0).is_none());
        assert!((power_to_db(100.0).unwrap() - 20.0).abs() < 1e-9);
    }

    /// Otsu splits a clean bimodal histogram between the two modes.
    #[test]
    fn otsu_splits_two_modes() {
        let mut v: Vec<f64> = std::iter::repeat_n(-20.0, 100).collect();
        v.extend(std::iter::repeat_n(-5.0, 100));
        let t = otsu_threshold(&v, 64).expect("bimodal data must yield a threshold");
        assert!(
            t > -20.0 && t < -5.0,
            "threshold {t} did not land between the modes"
        );
    }

    #[test]
    fn otsu_rejects_constant_data() {
        assert!(otsu_threshold(&[3.0; 50], 64).is_none());
        assert!(otsu_threshold(&[], 64).is_none());
    }

    /// Regions are 4-connected and the area filter drops small ones.
    #[test]
    fn regions_are_four_connected_and_filtered() {
        // 4x4. A 2x2 block top-left, plus a lone diagonal neighbour that must
        // NOT merge into it, plus an isolated single cell.
        let (rows, cols) = (4, 4);
        let mut flag = vec![false; rows * cols];
        for &i in &[0, 1, 4, 5] {
            flag[i] = true;
        }
        flag[10] = true; // diagonal from cell 5 -> separate under 4-connectivity
        let values = vec![2.0; rows * cols];

        let all = connected_regions(&flag, &values, rows, cols, 1);
        assert_eq!(all.len(), 2, "diagonal touch must not merge regions");

        let big = connected_regions(&flag, &values, rows, cols, 4);
        assert_eq!(big.len(), 1, "min_cells must drop the single-cell region");
        assert_eq!(big[0].cells.len(), 4);
        assert!((big[0].mean() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn empty_flag_gives_no_regions() {
        assert!(connected_regions(&[false; 9], &[0.0; 9], 3, 3, 1).is_empty());
    }

    /// A polygon wrapped in a `GeometryCollection` must still mask. The
    /// envelope prefilter used to return an inverted box for a collection, so
    /// every cell was rejected before `geometry_contains_point` — which does
    /// recurse — ever ran.
    #[test]
    fn geometry_collection_mask_is_not_prefiltered_away() {
        // A 4x4 raster over [0,4]x[0,4]; the mask covers the lower-left 2x2.
        let raster = Raster::new(RasterConfig {
            cols: 4,
            rows: 4,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(32610),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        let square = Geometry::Polygon {
            exterior: Ring(vec![
                Coord::xy(0.0, 0.0),
                Coord::xy(2.0, 0.0),
                Coord::xy(2.0, 2.0),
                Coord::xy(0.0, 2.0),
                Coord::xy(0.0, 0.0),
            ]),
            interiors: vec![],
        };
        let collected = Geometry::GeometryCollection(vec![square.clone()]);

        let mask_of = |g: Geometry| {
            let mut layer = Layer::new("mask").with_geom_type(GeometryType::Polygon);
            layer
                .add_feature(Some(g), &[])
                .expect("adding the mask feature must succeed");
            rasterize_mask(&raster, &layer, MaskSide::WaterPolygon)
        };

        let bare = mask_of(square);
        let wrapped = mask_of(collected);
        assert!(bare.iter().any(|&k| k), "the bare polygon must mask cells");
        assert_eq!(
            bare, wrapped,
            "wrapping the polygon in a collection must not change the mask"
        );
    }
}
