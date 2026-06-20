//! Minimal PMTiles v3 archive writer for PNG raster tiles.
//!
//! Implements just enough of the PMTiles v3 spec to pack a set of `(z, x, y)`
//! PNG tiles into a single clustered archive: a gzip-compressed root directory,
//! gzip-compressed JSON metadata, and the concatenated (deduplicated) tile data.
//! Spec: <https://github.com/protomaps/PMTiles/blob/main/spec/v3/spec.md>.

use std::collections::HashMap;
use std::io::Write;

use flate2::write::GzEncoder;
use flate2::Compression;
use wbcore::ToolError;

const HEADER_LEN: usize = 127;
const TILETYPE_PNG: u8 = 2;
const COMPRESSION_NONE: u8 = 1;
const COMPRESSION_GZIP: u8 = 2;

/// One tile to include in the archive.
pub struct Tile {
    pub z: u8,
    pub x: u32,
    pub y: u32,
    pub data: Vec<u8>,
}

/// Geographic bounds (degrees) for the header.
pub struct LonLatBounds {
    pub min_lon: f64,
    pub min_lat: f64,
    pub max_lon: f64,
    pub max_lat: f64,
}

struct Entry {
    tile_id: u64,
    offset: u64,
    length: u32,
    run_length: u32,
}

/// Builds a complete PMTiles v3 archive from PNG tiles.
///
/// `tiles` need not be sorted. Identical tile blobs are deduplicated. The whole
/// directory is stored in the header's root directory (no leaf directories),
/// which is fine for the bounded pyramids this crate produces.
pub fn build(
    mut tiles: Vec<Tile>,
    bounds: &LonLatBounds,
    min_zoom: u8,
    max_zoom: u8,
) -> Result<Vec<u8>, ToolError> {
    // Address each tile by its PMTiles (Hilbert) tile id and sort ascending.
    tiles.sort_by_key(|t| zxy_to_tile_id(t.z, t.x, t.y));

    // Concatenate tile data, deduplicating identical blobs by content.
    let mut tile_data: Vec<u8> = Vec::new();
    let mut seen: HashMap<Vec<u8>, (u64, u32)> = HashMap::new();
    let mut entries: Vec<Entry> = Vec::with_capacity(tiles.len());
    for t in &tiles {
        let (offset, length) = match seen.get(&t.data) {
            Some(&(off, len)) => (off, len),
            None => {
                let off = tile_data.len() as u64;
                let len = t.data.len() as u32;
                tile_data.extend_from_slice(&t.data);
                seen.insert(t.data.clone(), (off, len));
                (off, len)
            }
        };
        entries.push(Entry {
            tile_id: zxy_to_tile_id(t.z, t.x, t.y),
            offset,
            length,
            run_length: 1,
        });
    }

    let root_dir = gzip(&serialize_directory(&entries))?;
    let metadata = gzip(br#"{"type":"overlay","format":"png"}"#)?;

    let root_dir_offset = HEADER_LEN as u64;
    let root_dir_length = root_dir.len() as u64;
    let metadata_offset = root_dir_offset + root_dir_length;
    let metadata_length = metadata.len() as u64;
    let leaf_dirs_offset = metadata_offset + metadata_length;
    let tile_data_offset = leaf_dirs_offset; // no leaf directories
    let tile_data_length = tile_data.len() as u64;

    let header = build_header(
        root_dir_offset,
        root_dir_length,
        metadata_offset,
        metadata_length,
        leaf_dirs_offset,
        0,
        tile_data_offset,
        tile_data_length,
        tiles.len() as u64,
        entries.len() as u64,
        seen.len() as u64,
        min_zoom,
        max_zoom,
        bounds,
    );

    let mut out = Vec::with_capacity(
        HEADER_LEN + root_dir.len() + metadata.len() + tile_data.len(),
    );
    out.extend_from_slice(&header);
    out.extend_from_slice(&root_dir);
    out.extend_from_slice(&metadata);
    out.extend_from_slice(&tile_data);
    Ok(out)
}

/// PMTiles tile id for a `(z, x, y)` tile: the per-zoom base offset plus the
/// Hilbert-curve distance within the zoom. Mirrors the spec's reference
/// `zxyToTileId` exactly (note `s - 1 - x`, signed).
fn zxy_to_tile_id(z: u8, x: u32, y: u32) -> u64 {
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

/// Serializes directory entries (sorted by tile id) per the v3 spec: counts and
/// columns of varints with run-length / clustered-offset encoding.
fn serialize_directory(entries: &[Entry]) -> Vec<u8> {
    let mut buf = Vec::new();
    write_uvarint(&mut buf, entries.len() as u64);

    let mut last_id = 0u64;
    for e in entries {
        write_uvarint(&mut buf, e.tile_id - last_id);
        last_id = e.tile_id;
    }
    for e in entries {
        write_uvarint(&mut buf, e.run_length as u64);
    }
    for e in entries {
        write_uvarint(&mut buf, e.length as u64);
    }
    for (i, e) in entries.iter().enumerate() {
        if i > 0 && e.offset == entries[i - 1].offset + entries[i - 1].length as u64 {
            write_uvarint(&mut buf, 0);
        } else {
            write_uvarint(&mut buf, e.offset + 1);
        }
    }
    buf
}

#[allow(clippy::too_many_arguments)]
fn build_header(
    root_dir_offset: u64,
    root_dir_length: u64,
    metadata_offset: u64,
    metadata_length: u64,
    leaf_dirs_offset: u64,
    leaf_dirs_length: u64,
    tile_data_offset: u64,
    tile_data_length: u64,
    num_addressed_tiles: u64,
    num_tile_entries: u64,
    num_tile_contents: u64,
    min_zoom: u8,
    max_zoom: u8,
    b: &LonLatBounds,
) -> [u8; HEADER_LEN] {
    let mut h = [0u8; HEADER_LEN];
    h[0..7].copy_from_slice(b"PMTiles");
    h[7] = 3; // version
    let put = |h: &mut [u8; HEADER_LEN], at: usize, v: u64| {
        h[at..at + 8].copy_from_slice(&v.to_le_bytes());
    };
    put(&mut h, 8, root_dir_offset);
    put(&mut h, 16, root_dir_length);
    put(&mut h, 24, metadata_offset);
    put(&mut h, 32, metadata_length);
    put(&mut h, 40, leaf_dirs_offset);
    put(&mut h, 48, leaf_dirs_length);
    put(&mut h, 56, tile_data_offset);
    put(&mut h, 64, tile_data_length);
    put(&mut h, 72, num_addressed_tiles);
    put(&mut h, 80, num_tile_entries);
    put(&mut h, 88, num_tile_contents);
    h[96] = 1; // clustered
    h[97] = COMPRESSION_GZIP; // internal compression (directories + metadata)
    h[98] = COMPRESSION_NONE; // tile compression (PNG is already compressed)
    h[99] = TILETYPE_PNG;
    h[100] = min_zoom;
    h[101] = max_zoom;
    let e7 = |deg: f64| (deg * 1e7).round() as i32;
    let put_i32 = |h: &mut [u8; HEADER_LEN], at: usize, v: i32| {
        h[at..at + 4].copy_from_slice(&v.to_le_bytes());
    };
    put_i32(&mut h, 102, e7(b.min_lon));
    put_i32(&mut h, 106, e7(b.min_lat));
    put_i32(&mut h, 110, e7(b.max_lon));
    put_i32(&mut h, 114, e7(b.max_lat));
    h[118] = min_zoom; // center zoom
    put_i32(&mut h, 119, e7((b.min_lon + b.max_lon) / 2.0));
    put_i32(&mut h, 123, e7((b.min_lat + b.max_lat) / 2.0));
    h
}

fn gzip(data: &[u8]) -> Result<Vec<u8>, ToolError> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder
        .write_all(data)
        .and_then(|_| encoder.finish())
        .map_err(|e| ToolError::Execution(format!("gzip failed: {e}")))
}

/// Unsigned LEB128 varint.
fn write_uvarint(buf: &mut Vec<u8>, mut value: u64) {
    loop {
        let mut byte = (value & 0x7f) as u8;
        value >>= 7;
        if value != 0 {
            byte |= 0x80;
        }
        buf.push(byte);
        if value == 0 {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tile_ids_match_pmtiles_spec() {
        // Known values from the PMTiles v3 reference.
        assert_eq!(zxy_to_tile_id(0, 0, 0), 0);
        assert_eq!(zxy_to_tile_id(1, 0, 0), 1);
        assert_eq!(zxy_to_tile_id(1, 0, 1), 2);
        assert_eq!(zxy_to_tile_id(1, 1, 1), 3);
        assert_eq!(zxy_to_tile_id(1, 1, 0), 4);
        assert_eq!(zxy_to_tile_id(2, 0, 0), 5);
    }

    #[test]
    fn varint_roundtrips_small_and_large() {
        let mut b = Vec::new();
        write_uvarint(&mut b, 0);
        write_uvarint(&mut b, 127);
        write_uvarint(&mut b, 128);
        write_uvarint(&mut b, 300);
        assert_eq!(b, vec![0x00, 0x7f, 0x80, 0x01, 0xac, 0x02]);
    }

    #[test]
    fn build_emits_a_valid_header() {
        let tiles = vec![
            Tile { z: 1, x: 0, y: 0, data: vec![1, 2, 3] },
            Tile { z: 1, x: 1, y: 1, data: vec![1, 2, 3] }, // dup content
        ];
        let bounds = LonLatBounds { min_lon: -10.0, min_lat: -5.0, max_lon: 10.0, max_lat: 5.0 };
        let archive = build(tiles, &bounds, 1, 1).unwrap();
        assert_eq!(&archive[0..7], b"PMTiles");
        assert_eq!(archive[7], 3);
        assert_eq!(archive[99], TILETYPE_PNG);
        // num_tile_contents (offset 88) is 1 (the two tiles share one blob).
        assert_eq!(u64::from_le_bytes(archive[88..96].try_into().unwrap()), 1);
    }
}
