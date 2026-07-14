//! Assembles complete PMTiles v3 archives.
//!
//! [`assemble`] is the general path: caller-provided directory entries over an
//! already-laid-out tile-data section (what [`crate::extract`] produces).
//! [`build_png`] is the historical convenience used by GeoLibre's tiling tools:
//! pack loose PNG tiles into a fresh archive.

use std::collections::HashMap;

use crate::format::{
    compress, serialize_directory, Entry, Header, COMPRESSION_GZIP, COMPRESSION_NONE,
    HEADER_LEN, PRELUDE_LEN, TILETYPE_PNG,
};
use crate::{LonLatBounds, PmtilesError};

/// One tile to include in a [`build_png`] archive.
pub struct Tile {
    pub z: u8,
    pub x: u32,
    pub y: u32,
    pub data: Vec<u8>,
}

/// Header fields the caller decides; everything layout-related is computed.
#[derive(Clone, Debug)]
pub struct ArchiveParams {
    pub tile_type: u8,
    pub tile_compression: u8,
    pub min_zoom: u8,
    pub max_zoom: u8,
    pub bounds: LonLatBounds,
    pub center_zoom: u8,
    pub center_lon: f64,
    pub center_lat: f64,
    pub num_addressed_tiles: u64,
    pub num_tile_contents: u64,
}

/// Builds a complete archive from directory `entries` (sorted by tile id,
/// offsets relative to `tile_data`) plus the tile-data section and
/// *uncompressed* JSON metadata. Directories and metadata are gzipped. When
/// the root directory would overflow the spec's 16 KiB prelude, entries are
/// split into gzipped leaf directories.
pub fn assemble(
    entries: &[Entry],
    tile_data: &[u8],
    metadata: &[u8],
    params: &ArchiveParams,
) -> Result<Vec<u8>, PmtilesError> {
    debug_assert!(entries.windows(2).all(|w| w[0].tile_id < w[1].tile_id));

    let metadata = compress(metadata, COMPRESSION_GZIP)?;
    let (root_dir, leaves) = build_directories(entries)?;

    let root_dir_offset = HEADER_LEN as u64;
    let metadata_offset = root_dir_offset + root_dir.len() as u64;
    let leaf_dirs_offset = metadata_offset + metadata.len() as u64;
    let leaf_dirs_length: u64 = leaves.len() as u64;
    let tile_data_offset = leaf_dirs_offset + leaf_dirs_length;

    let clustered = entries
        .windows(2)
        .all(|w| w[1].offset >= w[0].offset);

    let e7 = |deg: f64| (deg * 1e7).round() as i32;
    let header = Header {
        root_dir_offset,
        root_dir_length: root_dir.len() as u64,
        metadata_offset,
        metadata_length: metadata.len() as u64,
        leaf_dirs_offset,
        leaf_dirs_length,
        tile_data_offset,
        tile_data_length: tile_data.len() as u64,
        num_addressed_tiles: params.num_addressed_tiles,
        num_tile_entries: entries.len() as u64,
        num_tile_contents: params.num_tile_contents,
        clustered,
        internal_compression: COMPRESSION_GZIP,
        tile_compression: params.tile_compression,
        tile_type: params.tile_type,
        min_zoom: params.min_zoom,
        max_zoom: params.max_zoom,
        min_lon_e7: e7(params.bounds.min_lon),
        min_lat_e7: e7(params.bounds.min_lat),
        max_lon_e7: e7(params.bounds.max_lon),
        max_lat_e7: e7(params.bounds.max_lat),
        center_zoom: params.center_zoom,
        center_lon_e7: e7(params.center_lon),
        center_lat_e7: e7(params.center_lat),
    };

    let mut out = Vec::with_capacity(
        HEADER_LEN + root_dir.len() + metadata.len() + leaves.len() + tile_data.len(),
    );
    out.extend_from_slice(&header.to_bytes());
    out.extend_from_slice(&root_dir);
    out.extend_from_slice(&metadata);
    out.extend_from_slice(&leaves);
    out.extend_from_slice(tile_data);
    Ok(out)
}

/// Serializes `entries` as a gzipped root directory, splitting into leaf
/// directories only when the root would not fit the 16 KiB prelude alongside
/// the header. Returns `(root_dir_gzipped, leaf_section)`.
fn build_directories(entries: &[Entry]) -> Result<(Vec<u8>, Vec<u8>), PmtilesError> {
    let root_budget = PRELUDE_LEN - HEADER_LEN;
    let root = compress(&serialize_directory(entries), COMPRESSION_GZIP)?;
    if root.len() <= root_budget {
        return Ok((root, Vec::new()));
    }

    // Split into leaves of `chunk` entries; double until the root of leaf
    // pointers fits. Each doubling roughly halves the root entry count, so
    // this terminates quickly even for millions of entries.
    let mut chunk = 4096usize;
    loop {
        let mut leaves: Vec<u8> = Vec::new();
        let mut root_entries: Vec<Entry> = Vec::with_capacity(entries.len() / chunk + 1);
        for group in entries.chunks(chunk) {
            let leaf = compress(&serialize_directory(group), COMPRESSION_GZIP)?;
            if leaf.len() > u32::MAX as usize {
                return Err(PmtilesError::Corrupt("leaf directory exceeds u32".into()));
            }
            root_entries.push(Entry {
                tile_id: group[0].tile_id,
                offset: leaves.len() as u64,
                length: leaf.len() as u32,
                run_length: 0,
            });
            leaves.extend_from_slice(&leaf);
        }
        let root = compress(&serialize_directory(&root_entries), COMPRESSION_GZIP)?;
        if root.len() <= root_budget {
            return Ok((root, leaves));
        }
        chunk *= 2;
    }
}

/// Builds a complete archive from loose PNG tiles (the original
/// `geolibre-tools` writer). `tiles` need not be sorted; identical blobs are
/// deduplicated by content.
pub fn build_png(
    mut tiles: Vec<Tile>,
    bounds: &LonLatBounds,
    min_zoom: u8,
    max_zoom: u8,
) -> Result<Vec<u8>, PmtilesError> {
    tiles.sort_by_key(|t| crate::tilemath::zxy_to_tile_id(t.z, t.x, t.y));

    let mut tile_data: Vec<u8> = Vec::new();
    let mut seen: HashMap<&[u8], (u64, u32)> = HashMap::new();
    let mut entries: Vec<Entry> = Vec::with_capacity(tiles.len());
    for t in &tiles {
        let (offset, length) = match seen.get(t.data.as_slice()) {
            Some(&(off, len)) => (off, len),
            None => {
                let off = tile_data.len() as u64;
                let len = t.data.len() as u32;
                tile_data.extend_from_slice(&t.data);
                seen.insert(&t.data, (off, len));
                (off, len)
            }
        };
        entries.push(Entry {
            tile_id: crate::tilemath::zxy_to_tile_id(t.z, t.x, t.y),
            offset,
            length,
            run_length: 1,
        });
    }

    let params = ArchiveParams {
        tile_type: TILETYPE_PNG,
        tile_compression: COMPRESSION_NONE, // PNG is already compressed
        min_zoom,
        max_zoom,
        bounds: *bounds,
        center_zoom: min_zoom,
        center_lon: (bounds.min_lon + bounds.max_lon) / 2.0,
        center_lat: (bounds.min_lat + bounds.max_lat) / 2.0,
        num_addressed_tiles: tiles.len() as u64,
        num_tile_contents: seen.len() as u64,
    };
    assemble(&entries, &tile_data, br#"{"type":"overlay","format":"png"}"#, &params)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::format::{decompress, parse_directory};

    fn world() -> LonLatBounds {
        LonLatBounds { min_lon: -10.0, min_lat: -5.0, max_lon: 10.0, max_lat: 5.0 }
    }

    #[test]
    fn build_png_emits_a_valid_deduplicated_archive() {
        let tiles = vec![
            Tile { z: 1, x: 0, y: 0, data: vec![1, 2, 3] },
            Tile { z: 1, x: 1, y: 1, data: vec![1, 2, 3] }, // dup content
        ];
        let archive = build_png(tiles, &world(), 1, 1).unwrap();
        let h = Header::parse(&archive).unwrap();
        assert_eq!(h.tile_type, TILETYPE_PNG);
        assert_eq!(h.num_tile_contents, 1);
        assert_eq!(h.num_addressed_tiles, 2);
        assert!(h.clustered);
        assert_eq!(h.leaf_dirs_length, 0);
        let root = decompress(
            &archive[h.root_dir_offset as usize..(h.root_dir_offset + h.root_dir_length) as usize],
            h.internal_compression,
        )
        .unwrap();
        let entries = parse_directory(&root).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].offset, entries[1].offset);
    }

    #[test]
    fn oversized_roots_split_into_leaves() {
        // Pseudo-random offsets/lengths keep the offset and length varint
        // columns incompressible, so the root genuinely overflows 16 KiB and
        // must split. (A cycling pattern gzips into the root budget.)
        let mut rng = 0x9e3779b97f4a7c15u64;
        let mut next = move || {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            rng
        };
        let entries: Vec<Entry> = (0..40_000u64)
            .map(|i| Entry {
                tile_id: i * 3, // gaps defeat contiguous-id delta encoding
                offset: next() % 1_000_000_007,
                length: 1_000 + (next() as u32 % 100_000),
                run_length: 1,
            })
            .collect();
        let params = ArchiveParams {
            tile_type: TILETYPE_PNG,
            tile_compression: COMPRESSION_NONE,
            min_zoom: 0,
            max_zoom: 8,
            bounds: world(),
            center_zoom: 0,
            center_lon: 0.0,
            center_lat: 0.0,
            num_addressed_tiles: 40_000,
            num_tile_contents: 7,
        };
        let archive = assemble(&entries, &[], b"{}", &params).unwrap();
        let h = Header::parse(&archive).unwrap();
        assert!(h.leaf_dirs_length > 0, "expected leaf directories");
        assert!(h.root_dir_offset + h.root_dir_length <= PRELUDE_LEN as u64);

        // Reassemble the full entry list through root → leaves and compare.
        let root = decompress(
            &archive[h.root_dir_offset as usize..(h.root_dir_offset + h.root_dir_length) as usize],
            h.internal_compression,
        )
        .unwrap();
        let root_entries = parse_directory(&root).unwrap();
        assert!(root_entries.iter().all(|e| e.run_length == 0));
        let mut reassembled: Vec<Entry> = Vec::new();
        for le in &root_entries {
            let start = (h.leaf_dirs_offset + le.offset) as usize;
            let leaf =
                decompress(&archive[start..start + le.length as usize], h.internal_compression)
                    .unwrap();
            reassembled.extend(parse_directory(&leaf).unwrap());
        }
        assert_eq!(reassembled, entries);
    }
}
