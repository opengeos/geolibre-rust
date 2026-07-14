//! ZXY ⇄ PMTiles tile ids and bbox → tile-id range planning.

use crate::{LonLatBounds, PmtilesError};

/// Web Mercator's latitude limit; tiles do not exist beyond it.
const MERCATOR_MAX_LAT: f64 = 85.051_128_78;

/// PMTiles tile id for a `(z, x, y)` tile: the per-zoom base offset plus the
/// Hilbert-curve distance within the zoom. Mirrors the spec's reference
/// `zxyToTileId` exactly (note `s - 1 - x`, signed).
pub fn zxy_to_tile_id(z: u8, x: u32, y: u32) -> u64 {
    let mut acc: u64 = 0;
    for t in 0..z {
        let dim = 1i64 << t;
        acc += (dim * dim) as u64;
    }
    let n = 1i64 << z;
    let (mut xx, mut yy) = (x as i64, y as i64);
    let mut d: i64 = 0;
    let mut s = n / 2;
    while s > 0 {
        let rx = i64::from((xx & s) > 0);
        let ry = i64::from((yy & s) > 0);
        d += s * s * ((3 * rx) ^ ry);
        if ry == 0 {
            if rx == 1 {
                xx = s - 1 - xx;
                yy = s - 1 - yy;
            }
            std::mem::swap(&mut xx, &mut yy);
        }
        s /= 2;
    }
    acc + d as u64
}

/// Slippy-map tile column/row containing `(lon, lat)` at zoom `z`, clamped to
/// the valid grid.
pub fn lonlat_to_tile(lon: f64, lat: f64, z: u8) -> (u32, u32) {
    let n = (1u64 << z) as f64;
    let lat = lat.clamp(-MERCATOR_MAX_LAT, MERCATOR_MAX_LAT);
    let x = ((lon + 180.0) / 360.0 * n).floor();
    let lat_rad = lat.to_radians();
    let y = ((1.0 - (lat_rad.tan() + 1.0 / lat_rad.cos()).ln() / std::f64::consts::PI) / 2.0
        * n)
        .floor();
    let max = (n - 1.0).max(0.0);
    (x.clamp(0.0, max) as u32, y.clamp(0.0, max) as u32)
}

/// The tile column/row rectangle (inclusive) covering `bbox` at zoom `z`.
pub fn bbox_tile_rect(bbox: &LonLatBounds, z: u8) -> (u32, u32, u32, u32) {
    let (x0, y0) = lonlat_to_tile(bbox.min_lon, bbox.max_lat, z); // NW corner
    let (x1, y1) = lonlat_to_tile(bbox.max_lon, bbox.min_lat, z); // SE corner
    (x0, y0, x1, y1)
}

/// Sorted, disjoint, half-open `[start, end)` tile-id ranges covering `bbox`
/// over `min_zoom..=max_zoom`, plus the total addressed-tile count.
///
/// Hilbert ids of a bbox rectangle are not contiguous, so this enumerates each
/// tile in the rectangle per zoom and coalesces consecutive ids. `max_tiles`
/// bounds the enumeration (and the caller's memory) before it starts.
pub fn bbox_tile_ranges(
    bbox: &LonLatBounds,
    min_zoom: u8,
    max_zoom: u8,
    max_tiles: u64,
) -> Result<(Vec<(u64, u64)>, u64), PmtilesError> {
    validate_request(bbox, min_zoom, max_zoom)?;

    // Count first so a too-large request fails before allocating anything.
    let mut total: u64 = 0;
    for z in min_zoom..=max_zoom {
        let (x0, y0, x1, y1) = bbox_tile_rect(bbox, z);
        total += u64::from(x1 - x0 + 1) * u64::from(y1 - y0 + 1);
    }
    if total > max_tiles {
        return Err(PmtilesError::TooManyTiles { requested: total, max: max_tiles });
    }

    let mut ids: Vec<u64> = Vec::with_capacity(total as usize);
    for z in min_zoom..=max_zoom {
        let (x0, y0, x1, y1) = bbox_tile_rect(bbox, z);
        for y in y0..=y1 {
            for x in x0..=x1 {
                ids.push(zxy_to_tile_id(z, x, y));
            }
        }
    }
    ids.sort_unstable();

    let mut ranges: Vec<(u64, u64)> = Vec::new();
    for id in ids {
        match ranges.last_mut() {
            Some((_, end)) if *end == id => *end += 1,
            _ => ranges.push((id, id + 1)),
        }
    }
    Ok((ranges, total))
}

/// Validates a bbox + zoom-range selection, shared by range planning and the
/// extractor so bad requests fail before any bytes are fetched.
pub fn validate_request(
    bbox: &LonLatBounds,
    min_zoom: u8,
    max_zoom: u8,
) -> Result<(), PmtilesError> {
    if !(bbox.min_lon < bbox.max_lon && bbox.min_lat < bbox.max_lat) {
        return Err(PmtilesError::InvalidRequest(format!(
            "bbox must have min < max on both axes (antimeridian-crossing \
             bboxes are not supported): [{}, {}, {}, {}]",
            bbox.min_lon, bbox.min_lat, bbox.max_lon, bbox.max_lat
        )));
    }
    if !(-180.0..=180.0).contains(&bbox.min_lon)
        || !(-180.0..=180.0).contains(&bbox.max_lon)
        || !(-90.0..=90.0).contains(&bbox.min_lat)
        || !(-90.0..=90.0).contains(&bbox.max_lat)
    {
        return Err(PmtilesError::InvalidRequest(
            "bbox out of range (lon in [-180, 180], lat in [-90, 90])".into(),
        ));
    }
    if min_zoom > max_zoom {
        return Err(PmtilesError::InvalidRequest(format!(
            "min_zoom {min_zoom} > max_zoom {max_zoom}"
        )));
    }
    Ok(())
}

/// Intersects `[start, end)` with sorted disjoint `ranges`, invoking `emit`
/// for each overlapping `[seg_start, seg_end)` segment.
pub fn intersect_ranges(
    start: u64,
    end: u64,
    ranges: &[(u64, u64)],
    mut emit: impl FnMut(u64, u64),
) {
    // First range whose end is past `start`.
    let mut i = ranges.partition_point(|&(_, e)| e <= start);
    while let Some(&(s, e)) = ranges.get(i) {
        if s >= end {
            break;
        }
        emit(s.max(start), e.min(end));
        i += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const WORLD: LonLatBounds = LonLatBounds {
        min_lon: -180.0,
        min_lat: -85.0,
        max_lon: 180.0,
        max_lat: 85.0,
    };

    #[test]
    fn tile_ids_match_pmtiles_spec() {
        assert_eq!(zxy_to_tile_id(0, 0, 0), 0);
        assert_eq!(zxy_to_tile_id(1, 0, 0), 1);
        assert_eq!(zxy_to_tile_id(1, 0, 1), 2);
        assert_eq!(zxy_to_tile_id(1, 1, 1), 3);
        assert_eq!(zxy_to_tile_id(1, 1, 0), 4);
        assert_eq!(zxy_to_tile_id(2, 0, 0), 5);
    }

    #[test]
    fn lonlat_to_tile_matches_slippy_reference() {
        assert_eq!(lonlat_to_tile(0.0, 0.0, 0), (0, 0));
        // (0,0) sits on the tile boundary: SE quadrant at z1.
        assert_eq!(lonlat_to_tile(0.0, 0.0, 1), (1, 1));
        assert_eq!(lonlat_to_tile(-180.0, 85.0, 2), (0, 0));
        // Boulder, CO at z10: x = (74.73/360)·1024 = 212.54, y = 387.64.
        assert_eq!(lonlat_to_tile(-105.27, 40.01, 10), (212, 387));
        // Poles clamp into the grid.
        assert_eq!(lonlat_to_tile(179.999, -89.9, 3), (7, 7));
    }

    #[test]
    fn world_bbox_selects_the_full_pyramid() {
        let (ranges, total) = bbox_tile_ranges(&WORLD, 0, 3, 1_000_000).unwrap();
        // 1 + 4 + 16 + 64
        assert_eq!(total, 85);
        // A full pyramid is one contiguous id range.
        assert_eq!(ranges, vec![(0, 85)]);
    }

    #[test]
    fn quadrant_bbox_selects_a_quarter_per_zoom() {
        let ne = LonLatBounds { min_lon: 1.0, min_lat: 1.0, max_lon: 179.0, max_lat: 84.0 };
        let (_, total) = bbox_tile_ranges(&ne, 1, 1, 1_000_000).unwrap();
        assert_eq!(total, 1);
        let (_, total) = bbox_tile_ranges(&ne, 2, 2, 1_000_000).unwrap();
        assert_eq!(total, 4);
    }

    #[test]
    fn tile_cap_is_enforced_before_allocation() {
        let err = bbox_tile_ranges(&WORLD, 0, 15, 1_000).unwrap_err();
        assert!(matches!(err, PmtilesError::TooManyTiles { .. }));
    }

    #[test]
    fn invalid_bboxes_are_rejected() {
        let flipped = LonLatBounds { min_lon: 10.0, min_lat: 0.0, max_lon: -10.0, max_lat: 5.0 };
        assert!(bbox_tile_ranges(&flipped, 0, 1, 100).is_err());
        let out = LonLatBounds { min_lon: -190.0, min_lat: 0.0, max_lon: 0.0, max_lat: 5.0 };
        assert!(bbox_tile_ranges(&out, 0, 1, 100).is_err());
    }

    #[test]
    fn intersect_ranges_emits_clipped_segments() {
        let ranges = [(0u64, 10u64), (20, 30), (40, 50)];
        let mut got = Vec::new();
        intersect_ranges(5, 45, &ranges, |s, e| got.push((s, e)));
        assert_eq!(got, vec![(5, 10), (20, 30), (40, 45)]);
        let mut none = Vec::new();
        intersect_ranges(10, 20, &ranges, |s, e| none.push((s, e)));
        assert!(none.is_empty());
    }
}
