//! PMTiles v3 wire format: header, directories, varints, internal compression.

use std::io::{Read, Write};

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

use crate::{LonLatBounds, PmtilesError};

pub const HEADER_LEN: usize = 127;
/// The spec guarantees the header and root directory fit in the first 16 KiB.
pub const PRELUDE_LEN: usize = 16_384;

pub const COMPRESSION_UNKNOWN: u8 = 0;
pub const COMPRESSION_NONE: u8 = 1;
pub const COMPRESSION_GZIP: u8 = 2;

pub const TILETYPE_MVT: u8 = 1;
pub const TILETYPE_PNG: u8 = 2;

/// Parsed 127-byte PMTiles v3 header. Field order mirrors the spec.
#[derive(Clone, Debug, PartialEq)]
pub struct Header {
    pub root_dir_offset: u64,
    pub root_dir_length: u64,
    pub metadata_offset: u64,
    pub metadata_length: u64,
    pub leaf_dirs_offset: u64,
    pub leaf_dirs_length: u64,
    pub tile_data_offset: u64,
    pub tile_data_length: u64,
    pub num_addressed_tiles: u64,
    pub num_tile_entries: u64,
    pub num_tile_contents: u64,
    pub clustered: bool,
    pub internal_compression: u8,
    pub tile_compression: u8,
    pub tile_type: u8,
    pub min_zoom: u8,
    pub max_zoom: u8,
    pub min_lon_e7: i32,
    pub min_lat_e7: i32,
    pub max_lon_e7: i32,
    pub max_lat_e7: i32,
    pub center_zoom: u8,
    pub center_lon_e7: i32,
    pub center_lat_e7: i32,
}

impl Header {
    pub fn parse(bytes: &[u8]) -> Result<Header, PmtilesError> {
        if bytes.len() < HEADER_LEN {
            return Err(PmtilesError::Truncated("header"));
        }
        if &bytes[0..7] != b"PMTiles" {
            return Err(PmtilesError::NotPmtiles);
        }
        if bytes[7] != 3 {
            return Err(PmtilesError::UnsupportedVersion(bytes[7]));
        }
        let u = |at: usize| u64::from_le_bytes(bytes[at..at + 8].try_into().unwrap());
        let i = |at: usize| i32::from_le_bytes(bytes[at..at + 4].try_into().unwrap());
        Ok(Header {
            root_dir_offset: u(8),
            root_dir_length: u(16),
            metadata_offset: u(24),
            metadata_length: u(32),
            leaf_dirs_offset: u(40),
            leaf_dirs_length: u(48),
            tile_data_offset: u(56),
            tile_data_length: u(64),
            num_addressed_tiles: u(72),
            num_tile_entries: u(80),
            num_tile_contents: u(88),
            clustered: bytes[96] == 1,
            internal_compression: bytes[97],
            tile_compression: bytes[98],
            tile_type: bytes[99],
            min_zoom: bytes[100],
            max_zoom: bytes[101],
            min_lon_e7: i(102),
            min_lat_e7: i(106),
            max_lon_e7: i(110),
            max_lat_e7: i(114),
            center_zoom: bytes[118],
            center_lon_e7: i(119),
            center_lat_e7: i(123),
        })
    }

    pub fn to_bytes(&self) -> [u8; HEADER_LEN] {
        let mut h = [0u8; HEADER_LEN];
        h[0..7].copy_from_slice(b"PMTiles");
        h[7] = 3;
        let put = |h: &mut [u8; HEADER_LEN], at: usize, v: u64| {
            h[at..at + 8].copy_from_slice(&v.to_le_bytes());
        };
        let put_i32 = |h: &mut [u8; HEADER_LEN], at: usize, v: i32| {
            h[at..at + 4].copy_from_slice(&v.to_le_bytes());
        };
        put(&mut h, 8, self.root_dir_offset);
        put(&mut h, 16, self.root_dir_length);
        put(&mut h, 24, self.metadata_offset);
        put(&mut h, 32, self.metadata_length);
        put(&mut h, 40, self.leaf_dirs_offset);
        put(&mut h, 48, self.leaf_dirs_length);
        put(&mut h, 56, self.tile_data_offset);
        put(&mut h, 64, self.tile_data_length);
        put(&mut h, 72, self.num_addressed_tiles);
        put(&mut h, 80, self.num_tile_entries);
        put(&mut h, 88, self.num_tile_contents);
        h[96] = u8::from(self.clustered);
        h[97] = self.internal_compression;
        h[98] = self.tile_compression;
        h[99] = self.tile_type;
        h[100] = self.min_zoom;
        h[101] = self.max_zoom;
        put_i32(&mut h, 102, self.min_lon_e7);
        put_i32(&mut h, 106, self.min_lat_e7);
        put_i32(&mut h, 110, self.max_lon_e7);
        put_i32(&mut h, 114, self.max_lat_e7);
        h[118] = self.center_zoom;
        put_i32(&mut h, 119, self.center_lon_e7);
        put_i32(&mut h, 123, self.center_lat_e7);
        h
    }

    pub fn bounds(&self) -> LonLatBounds {
        LonLatBounds {
            min_lon: self.min_lon_e7 as f64 / 1e7,
            min_lat: self.min_lat_e7 as f64 / 1e7,
            max_lon: self.max_lon_e7 as f64 / 1e7,
            max_lat: self.max_lat_e7 as f64 / 1e7,
        }
    }
}

/// One directory entry. `run_length == 0` marks a leaf-directory pointer
/// (offset relative to the leaf section); otherwise a run of `run_length`
/// consecutive tile ids sharing one blob (offset relative to the tile data
/// section).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Entry {
    pub tile_id: u64,
    pub offset: u64,
    pub length: u32,
    pub run_length: u32,
}

/// Serializes directory entries (sorted by tile id) per the v3 spec: an entry
/// count then columns of varints with delta / run-length / clustered-offset
/// encoding.
pub fn serialize_directory(entries: &[Entry]) -> Vec<u8> {
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

/// Parses an (uncompressed) directory. Inverse of [`serialize_directory`].
pub fn parse_directory(bytes: &[u8]) -> Result<Vec<Entry>, PmtilesError> {
    let mut pos = 0usize;
    let n = read_uvarint(bytes, &mut pos)? as usize;
    // Guard against a corrupt count causing a huge allocation: every entry
    // needs at least 4 varint bytes (one per column).
    if n > bytes.len() {
        return Err(PmtilesError::Corrupt(format!(
            "directory claims {n} entries in {} bytes",
            bytes.len()
        )));
    }
    let mut entries = vec![
        Entry { tile_id: 0, offset: 0, length: 0, run_length: 0 };
        n
    ];
    let mut last_id = 0u64;
    for e in entries.iter_mut() {
        last_id = last_id
            .checked_add(read_uvarint(bytes, &mut pos)?)
            .ok_or(PmtilesError::Corrupt("tile id overflow".into()))?;
        e.tile_id = last_id;
    }
    for e in entries.iter_mut() {
        e.run_length = read_uvarint(bytes, &mut pos)? as u32;
    }
    for e in entries.iter_mut() {
        e.length = read_uvarint(bytes, &mut pos)? as u32;
    }
    for i in 0..n {
        let v = read_uvarint(bytes, &mut pos)?;
        entries[i].offset = if v == 0 {
            if i == 0 {
                return Err(PmtilesError::Corrupt(
                    "first directory entry has relative offset".into(),
                ));
            }
            entries[i - 1].offset + entries[i - 1].length as u64
        } else {
            v - 1
        };
    }
    Ok(entries)
}

/// Compresses `data` with a PMTiles internal-compression codec.
pub fn compress(data: &[u8], codec: u8) -> Result<Vec<u8>, PmtilesError> {
    match codec {
        COMPRESSION_NONE => Ok(data.to_vec()),
        COMPRESSION_GZIP => {
            let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
            encoder
                .write_all(data)
                .and_then(|_| encoder.finish())
                .map_err(|e| PmtilesError::Codec(format!("gzip failed: {e}")))
        }
        other => Err(PmtilesError::UnsupportedCompression(other)),
    }
}

/// Decompresses `data` with a PMTiles internal-compression codec.
pub fn decompress(data: &[u8], codec: u8) -> Result<Vec<u8>, PmtilesError> {
    match codec {
        COMPRESSION_NONE => Ok(data.to_vec()),
        COMPRESSION_GZIP => {
            let mut out = Vec::new();
            GzDecoder::new(data)
                .read_to_end(&mut out)
                .map_err(|e| PmtilesError::Codec(format!("gunzip failed: {e}")))?;
            Ok(out)
        }
        other => Err(PmtilesError::UnsupportedCompression(other)),
    }
}

/// Unsigned LEB128 varint.
pub fn write_uvarint(buf: &mut Vec<u8>, mut value: u64) {
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

/// Unsigned LEB128 varint read, advancing `pos`.
pub fn read_uvarint(bytes: &[u8], pos: &mut usize) -> Result<u64, PmtilesError> {
    let mut value = 0u64;
    let mut shift = 0u32;
    loop {
        let byte = *bytes
            .get(*pos)
            .ok_or(PmtilesError::Truncated("varint"))?;
        *pos += 1;
        if shift == 63 && byte > 1 {
            return Err(PmtilesError::Corrupt("varint overflows u64".into()));
        }
        value |= u64::from(byte & 0x7f) << shift;
        if byte & 0x80 == 0 {
            return Ok(value);
        }
        shift += 7;
        if shift > 63 {
            return Err(PmtilesError::Corrupt("varint too long".into()));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn varint_roundtrips() {
        for v in [0u64, 1, 127, 128, 300, 16_383, 16_384, u32::MAX as u64, u64::MAX] {
            let mut buf = Vec::new();
            write_uvarint(&mut buf, v);
            let mut pos = 0;
            assert_eq!(read_uvarint(&buf, &mut pos).unwrap(), v);
            assert_eq!(pos, buf.len());
        }
    }

    #[test]
    fn varint_rejects_truncation_and_overflow() {
        let mut pos = 0;
        assert_eq!(
            read_uvarint(&[0x80], &mut pos),
            Err(PmtilesError::Truncated("varint"))
        );
        let mut pos = 0;
        let too_long = [0x80u8; 11];
        assert!(read_uvarint(&too_long, &mut pos).is_err());
    }

    #[test]
    fn header_roundtrips() {
        let h = Header {
            root_dir_offset: 127,
            root_dir_length: 42,
            metadata_offset: 169,
            metadata_length: 10,
            leaf_dirs_offset: 179,
            leaf_dirs_length: 0,
            tile_data_offset: 179,
            tile_data_length: 1000,
            num_addressed_tiles: 12,
            num_tile_entries: 8,
            num_tile_contents: 5,
            clustered: true,
            internal_compression: COMPRESSION_GZIP,
            tile_compression: COMPRESSION_GZIP,
            tile_type: TILETYPE_MVT,
            min_zoom: 0,
            max_zoom: 15,
            min_lon_e7: -1_800_000_000,
            min_lat_e7: -850_511_287,
            max_lon_e7: 1_800_000_000,
            max_lat_e7: 850_511_287,
            center_zoom: 7,
            center_lon_e7: 123,
            center_lat_e7: -456,
        };
        assert_eq!(Header::parse(&h.to_bytes()).unwrap(), h);
    }

    #[test]
    fn directory_roundtrips_with_runs_and_clustered_offsets() {
        let entries = vec![
            Entry { tile_id: 0, offset: 0, length: 100, run_length: 1 },
            // clustered: offset == previous offset + length → encoded as 0
            Entry { tile_id: 1, offset: 100, length: 50, run_length: 4 },
            // repeated blob (same offset as an earlier entry)
            Entry { tile_id: 10, offset: 0, length: 100, run_length: 1 },
            // leaf pointer
            Entry { tile_id: 20, offset: 9_999, length: 77, run_length: 0 },
        ];
        let bytes = serialize_directory(&entries);
        assert_eq!(parse_directory(&bytes).unwrap(), entries);
    }

    #[test]
    fn parse_rejects_bad_magic_and_version() {
        let mut bytes = [0u8; HEADER_LEN];
        assert_eq!(Header::parse(&bytes), Err(PmtilesError::NotPmtiles));
        bytes[0..7].copy_from_slice(b"PMTiles");
        bytes[7] = 2;
        assert_eq!(Header::parse(&bytes), Err(PmtilesError::UnsupportedVersion(2)));
    }
}
