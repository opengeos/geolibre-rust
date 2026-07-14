//! Sans-IO bbox/zoom extraction from a remote PMTiles archive.
//!
//! [`Extractor`] never performs I/O. It publishes the absolute byte ranges it
//! needs ([`Extractor::wanted`]); the host fetches them however it likes
//! (browser `fetch` with `Range:` headers, reqwest, file reads) and hands the
//! bytes back ([`Extractor::feed`]), in any order and with any concurrency.
//! When [`Extractor::is_done`] turns true, [`Extractor::finish`] assembles a
//! self-contained clustered archive of the selected tiles.
//!
//! Protocol phases (each phase's requests appear in `wanted()` once the
//! previous data arrives):
//!
//! 1. The 16 KiB prelude — spec-guaranteed to contain the header and root
//!    directory, and usually the metadata too.
//! 2. Any leaf directories overlapping the selection (plus the metadata if it
//!    fell outside the prelude).
//! 3. The tile-data ranges, coalesced up to `max_range_gap` so nearby blobs
//!    ride one request.

use std::collections::HashMap;

use crate::format::{decompress, Entry, Header, HEADER_LEN, PRELUDE_LEN};
use crate::tilemath::{bbox_tile_ranges, intersect_ranges, validate_request};
use crate::writer::{assemble, ArchiveParams};
use crate::{LonLatBounds, PmtilesError};

/// An absolute byte range the host must fetch: `offset .. offset + length`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ByteRange {
    pub offset: u64,
    pub length: u64,
}

/// Extraction parameters.
#[derive(Clone, Debug)]
pub struct ExtractOptions {
    pub bbox: LonLatBounds,
    /// Lowest zoom to include. 0 keeps the map usable when zoomed out.
    pub min_zoom: u8,
    /// Highest zoom to include (clamped to what the source archive has).
    pub max_zoom: u8,
    /// Cap on addressed tiles, so an oversized selection fails fast instead of
    /// exhausting memory. Default 2 million (≈ a large metro area at z15).
    pub max_tiles: u64,
    /// Coalesce tile-data requests whose gap is at most this many bytes.
    /// Trades a little overfetch for far fewer HTTP round-trips. Default 64 KiB.
    pub max_range_gap: u64,
}

impl ExtractOptions {
    pub fn new(bbox: LonLatBounds, min_zoom: u8, max_zoom: u8) -> Self {
        ExtractOptions { bbox, min_zoom, max_zoom, max_tiles: 2_000_000, max_range_gap: 65_536 }
    }
}

/// Extraction progress for UIs.
#[derive(Clone, Copy, Debug, Default)]
pub struct Progress {
    /// "header" | "directories" | "data" | "done"
    pub phase: &'static str,
    /// Addressed tiles selected (known once directories are processed).
    pub tiles_selected: u64,
    /// Distinct tile blobs to download.
    pub blobs_total: u64,
    /// Total bytes of planned tile-data requests (includes gap overfetch).
    pub data_bytes_total: u64,
    /// Tile-data bytes received so far.
    pub data_bytes_received: u64,
    /// Rough size of the archive `finish()` will produce.
    pub estimated_output_bytes: u64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Tag {
    Prelude,
    RootDir,
    Leaf,
    Metadata,
    Data { index: usize },
}

pub struct Extractor {
    opts: ExtractOptions,
    pending: Vec<(ByteRange, Tag)>,
    header: Option<Header>,
    /// Effective zoom range after clamping to the source archive.
    eff_zoom: (u8, u8),
    /// Sorted disjoint half-open tile-id ranges wanted by the selection.
    wanted_ids: Vec<(u64, u64)>,
    /// Selected entry segments; `offset` stays source-relative until `finish`.
    selected: Vec<Entry>,
    root_processed: bool,
    leaves_outstanding: usize,
    /// Metadata bytes exactly as stored in the source (still compressed with
    /// the source's internal codec).
    metadata_raw: Option<Vec<u8>>,
    planned: bool,
    /// Merged absolute tile-data ranges, sorted by offset.
    data_ranges: Vec<ByteRange>,
    data: Vec<Option<Vec<u8>>>,
    // Progress accounting.
    blobs_total: u64,
    blob_bytes: u64,
    data_bytes_total: u64,
    data_bytes_received: u64,
}

impl Extractor {
    pub fn new(opts: ExtractOptions) -> Result<Extractor, PmtilesError> {
        validate_request(&opts.bbox, opts.min_zoom, opts.max_zoom)?;
        Ok(Extractor {
            opts,
            pending: vec![(
                ByteRange { offset: 0, length: PRELUDE_LEN as u64 },
                Tag::Prelude,
            )],
            header: None,
            eff_zoom: (0, 0),
            wanted_ids: Vec::new(),
            selected: Vec::new(),
            root_processed: false,
            leaves_outstanding: 0,
            metadata_raw: None,
            planned: false,
            data_ranges: Vec::new(),
            data: Vec::new(),
            blobs_total: 0,
            blob_bytes: 0,
            data_bytes_total: 0,
            data_bytes_received: 0,
        })
    }

    /// The source header, once the prelude has been fed.
    pub fn header(&self) -> Option<&Header> {
        self.header.as_ref()
    }

    /// Adjusts the addressed-tile cap. Effective if set before the prelude
    /// (where the selection is enumerated) is fed.
    pub fn set_max_tiles(&mut self, max_tiles: u64) {
        self.opts.max_tiles = max_tiles;
    }

    /// Adjusts request coalescing. Effective if set before the directories
    /// finish (when tile-data requests are planned).
    pub fn set_max_range_gap(&mut self, max_gap: u64) {
        self.opts.max_range_gap = max_gap;
    }

    /// Byte ranges the host should fetch next. Empty only when done (or when
    /// every outstanding request is already being fetched by the host).
    pub fn wanted(&self) -> Vec<ByteRange> {
        self.pending.iter().map(|(r, _)| *r).collect()
    }

    pub fn is_done(&self) -> bool {
        self.planned && self.pending.is_empty()
    }

    pub fn progress(&self) -> Progress {
        let phase = if self.header.is_none() {
            "header"
        } else if !self.planned {
            "directories"
        } else if !self.pending.is_empty() {
            "data"
        } else {
            "done"
        };
        Progress {
            phase,
            tiles_selected: self.selected.iter().map(|e| e.run_length as u64).sum(),
            blobs_total: self.blobs_total,
            data_bytes_total: self.data_bytes_total,
            data_bytes_received: self.data_bytes_received,
            // Blobs + metadata + a generous directory/header allowance.
            estimated_output_bytes: if self.planned {
                self.blob_bytes
                    + self.metadata_raw.as_ref().map_or(0, |m| m.len() as u64)
                    + PRELUDE_LEN as u64
            } else {
                0
            },
        }
    }

    /// Hands the bytes of one previously-`wanted()` range back to the
    /// extractor. `offset` identifies the request; ranges may arrive in any
    /// order.
    pub fn feed(&mut self, offset: u64, bytes: &[u8]) -> Result<(), PmtilesError> {
        let idx = self
            .pending
            .iter()
            .position(|(r, _)| r.offset == offset)
            .ok_or_else(|| {
                PmtilesError::Protocol(format!("no outstanding request at offset {offset}"))
            })?;
        let (range, tag) = self.pending.swap_remove(idx);

        // The prelude may legitimately come back short (file smaller than
        // 16 KiB); every other range must be complete.
        if tag != Tag::Prelude && (bytes.len() as u64) < range.length {
            self.pending.push((range, tag));
            return Err(PmtilesError::Protocol(format!(
                "range at {offset} returned {} bytes, expected {}",
                bytes.len(),
                range.length
            )));
        }
        let bytes = &bytes[..bytes.len().min(range.length as usize)];

        match tag {
            Tag::Prelude => self.process_prelude(bytes)?,
            Tag::RootDir => {
                self.process_directory(bytes, 0)?;
                self.root_processed = true;
            }
            Tag::Leaf => {
                self.leaves_outstanding -= 1;
                self.process_directory(bytes, 1)?;
            }
            Tag::Metadata => self.metadata_raw = Some(bytes.to_vec()),
            Tag::Data { index } => {
                self.data_bytes_received += bytes.len() as u64;
                self.data[index] = Some(bytes.to_vec());
            }
        }
        self.maybe_plan()?;
        Ok(())
    }

    /// Builds the output archive. Call once `is_done()` is true.
    pub fn finish(mut self) -> Result<Vec<u8>, PmtilesError> {
        if !self.is_done() {
            return Err(PmtilesError::Protocol(
                "finish() called before extraction completed".into(),
            ));
        }
        let header = self.header.clone().expect("done implies header");

        let mut selected = std::mem::take(&mut self.selected);
        selected.sort_unstable_by_key(|e| e.tile_id);

        // Lay out the new tile-data section in tile-id order, deduplicating
        // blobs, and merge segments that are contiguous both in id space and
        // in content back into runs.
        let mut out_data: Vec<u8> = Vec::with_capacity(self.blob_bytes as usize);
        let mut new_offsets: HashMap<(u64, u32), u64> = HashMap::new();
        let mut entries: Vec<Entry> = Vec::with_capacity(selected.len());
        for seg in &selected {
            let new_offset = match new_offsets.get(&(seg.offset, seg.length)) {
                Some(&off) => off,
                None => {
                    let off = out_data.len() as u64;
                    out_data.extend_from_slice(self.blob_bytes_for(&header, seg)?);
                    new_offsets.insert((seg.offset, seg.length), off);
                    off
                }
            };
            match entries.last_mut() {
                Some(last)
                    if last.offset == new_offset
                        && last.length == seg.length
                        && last.tile_id + last.run_length as u64 == seg.tile_id =>
                {
                    last.run_length += seg.run_length;
                }
                _ => entries.push(Entry {
                    tile_id: seg.tile_id,
                    offset: new_offset,
                    length: seg.length,
                    run_length: seg.run_length,
                }),
            }
        }

        let metadata_raw = self.metadata_raw.expect("done implies metadata");
        let metadata = if metadata_raw.is_empty() {
            b"{}".to_vec()
        } else {
            decompress(&metadata_raw, header.internal_compression)?
        };

        // Output bounds: the request clipped to what the source claims to
        // cover (fall back to the request if the source header is degenerate).
        let src = header.bounds();
        let mut bounds = LonLatBounds {
            min_lon: self.opts.bbox.min_lon.max(src.min_lon),
            min_lat: self.opts.bbox.min_lat.max(src.min_lat),
            max_lon: self.opts.bbox.max_lon.min(src.max_lon),
            max_lat: self.opts.bbox.max_lat.min(src.max_lat),
        };
        if bounds.min_lon >= bounds.max_lon || bounds.min_lat >= bounds.max_lat {
            bounds = self.opts.bbox;
        }

        let (min_zoom, max_zoom) = self.eff_zoom;
        let params = ArchiveParams {
            tile_type: header.tile_type,
            tile_compression: header.tile_compression,
            min_zoom,
            max_zoom,
            bounds,
            center_zoom: header.center_zoom.clamp(min_zoom, max_zoom),
            center_lon: (bounds.min_lon + bounds.max_lon) / 2.0,
            center_lat: (bounds.min_lat + bounds.max_lat) / 2.0,
            num_addressed_tiles: entries.iter().map(|e| e.run_length as u64).sum(),
            num_tile_contents: new_offsets.len() as u64,
        };
        assemble(&entries, &out_data, &metadata, &params)
    }

    fn process_prelude(&mut self, bytes: &[u8]) -> Result<(), PmtilesError> {
        if bytes.len() < HEADER_LEN {
            return Err(PmtilesError::Truncated("header"));
        }
        let header = Header::parse(bytes)?;
        // Fail on unsupported internal compression now, with a clear message,
        // rather than deep inside directory parsing.
        decompress(&[], header.internal_compression).and(Ok(())).or_else(|e| {
            if matches!(e, PmtilesError::UnsupportedCompression(_)) { Err(e) } else { Ok(()) }
        })?;

        let eff_min = self.opts.min_zoom.max(header.min_zoom);
        let eff_max = self.opts.max_zoom.min(header.max_zoom);
        if eff_min > eff_max {
            return Err(PmtilesError::InvalidRequest(format!(
                "requested zooms z{}-z{} do not intersect the source archive (z{}-z{})",
                self.opts.min_zoom, self.opts.max_zoom, header.min_zoom, header.max_zoom
            )));
        }
        self.eff_zoom = (eff_min, eff_max);
        let (wanted, _) =
            bbox_tile_ranges(&self.opts.bbox, eff_min, eff_max, self.opts.max_tiles)?;
        self.wanted_ids = wanted;

        // Metadata: empty, already in the prelude, or a follow-up request.
        let (m_off, m_len) = (header.metadata_offset, header.metadata_length);
        if m_len == 0 {
            self.metadata_raw = Some(Vec::new());
        } else if m_off + m_len <= bytes.len() as u64 {
            self.metadata_raw =
                Some(bytes[m_off as usize..(m_off + m_len) as usize].to_vec());
        } else {
            self.pending
                .push((ByteRange { offset: m_off, length: m_len }, Tag::Metadata));
        }

        // Root directory: the spec puts it inside the prelude, but tolerate
        // archives that don't by fetching it explicitly.
        let (r_off, r_len) = (header.root_dir_offset, header.root_dir_length);
        self.header = Some(header);
        if r_off + r_len <= bytes.len() as u64 {
            let root = bytes[r_off as usize..(r_off + r_len) as usize].to_vec();
            self.process_directory(&root, 0)?;
            self.root_processed = true;
        } else {
            self.pending
                .push((ByteRange { offset: r_off, length: r_len }, Tag::RootDir));
        }
        Ok(())
    }

    fn process_directory(&mut self, raw: &[u8], depth: u8) -> Result<(), PmtilesError> {
        let header = self.header.as_ref().expect("directory after header");
        let dir = decompress(raw, header.internal_compression)?;
        let entries = crate::format::parse_directory(&dir)?;
        let leaf_base = header.leaf_dirs_offset;

        let mut new_leaves: Vec<ByteRange> = Vec::new();
        for (i, e) in entries.iter().enumerate() {
            if e.run_length == 0 {
                // Leaf pointer, covering ids up to the next entry.
                if depth > 0 {
                    return Err(PmtilesError::Corrupt("leaf directory inside a leaf".into()));
                }
                if e.length == 0 {
                    return Err(PmtilesError::Corrupt("zero-length leaf directory".into()));
                }
                let coverage_end = entries
                    .get(i + 1)
                    .map(|n| n.tile_id)
                    .unwrap_or(u64::MAX);
                let mut hit = false;
                intersect_ranges(e.tile_id, coverage_end, &self.wanted_ids, |_, _| hit = true);
                if hit {
                    new_leaves.push(ByteRange {
                        offset: leaf_base + e.offset,
                        length: e.length as u64,
                    });
                }
            } else {
                let end = e.tile_id + e.run_length as u64;
                let (offset, length) = (e.offset, e.length);
                intersect_ranges(e.tile_id, end, &self.wanted_ids, |s, seg_end| {
                    self.selected.push(Entry {
                        tile_id: s,
                        offset,
                        length,
                        run_length: (seg_end - s) as u32,
                    });
                });
            }
        }
        self.leaves_outstanding += new_leaves.len();
        self.pending
            .extend(new_leaves.into_iter().map(|r| (r, Tag::Leaf)));
        Ok(())
    }

    /// Once every directory is in, plan the merged tile-data requests.
    fn maybe_plan(&mut self) -> Result<(), PmtilesError> {
        if self.planned || !self.root_processed || self.leaves_outstanding > 0 {
            return Ok(());
        }
        let header = self.header.as_ref().expect("planning after header");
        if self.selected.is_empty() {
            return Err(PmtilesError::InvalidRequest(
                "no tiles matched the requested bbox and zoom range in this archive".into(),
            ));
        }

        // Distinct blobs, sorted by source offset.
        let mut blobs: Vec<(u64, u32)> =
            self.selected.iter().map(|e| (e.offset, e.length)).collect();
        blobs.sort_unstable();
        blobs.dedup();
        self.blobs_total = blobs.len() as u64;
        self.blob_bytes = blobs.iter().map(|&(_, l)| l as u64).sum();

        // Coalesce into absolute ranges, bridging gaps up to max_range_gap.
        let mut ranges: Vec<ByteRange> = Vec::new();
        for &(off, len) in &blobs {
            if len == 0 {
                continue;
            }
            let start = header.tile_data_offset + off;
            let end = start + len as u64;
            match ranges.last_mut() {
                Some(last)
                    if start <= last.offset + last.length + self.opts.max_range_gap =>
                {
                    last.length = last.length.max(end - last.offset);
                }
                _ => ranges.push(ByteRange { offset: start, length: len as u64 }),
            }
        }
        self.data_bytes_total = ranges.iter().map(|r| r.length).sum();
        self.data = vec![None; ranges.len()];
        self.pending.extend(
            ranges
                .iter()
                .enumerate()
                .map(|(index, r)| (*r, Tag::Data { index })),
        );
        self.data_ranges = ranges;
        self.planned = true;
        Ok(())
    }

    /// The bytes of one selected blob, sliced out of the fetched data ranges.
    fn blob_bytes_for<'a>(
        &'a self,
        header: &Header,
        seg: &Entry,
    ) -> Result<&'a [u8], PmtilesError> {
        if seg.length == 0 {
            return Ok(&[]);
        }
        let abs = header.tile_data_offset + seg.offset;
        let idx = self
            .data_ranges
            .partition_point(|r| r.offset <= abs)
            .checked_sub(1)
            .ok_or_else(|| PmtilesError::Corrupt("tile offset before data section".into()))?;
        let range = self.data_ranges[idx];
        if abs + seg.length as u64 > range.offset + range.length {
            return Err(PmtilesError::Corrupt(format!(
                "tile at source offset {} escapes its planned range",
                seg.offset
            )));
        }
        let bytes = self.data[idx]
            .as_ref()
            .expect("done implies all data ranges fetched");
        let rel = (abs - range.offset) as usize;
        Ok(&bytes[rel..rel + seg.length as usize])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::format::{parse_directory, COMPRESSION_NONE, TILETYPE_PNG};
    use crate::tilemath::zxy_to_tile_id;
    use crate::writer::{build_png, Tile};

    /// Distinct, recognizable payload for tile (z, x, y). Padded to a
    /// pseudo-random per-tile length so directory length columns stay
    /// incompressible — otherwise even huge pyramids gzip into a root-only
    /// directory and the leaf paths never exercise.
    fn payload(z: u8, x: u32, y: u32) -> Vec<u8> {
        let mut h = (u64::from(x) << 40) ^ (u64::from(y) << 16) ^ u64::from(z)
            ^ 0x9e37_79b9_7f4a_7c15;
        h ^= h << 13;
        h ^= h >> 7;
        h ^= h << 17;
        let mut data = format!("tile-{z}/{x}/{y}").into_bytes();
        data.resize(data.len() + (h % 700) as usize, 0xAB);
        data
    }

    /// A full z0..=max pyramid archive with unique payloads.
    fn pyramid_archive(max_zoom: u8) -> Vec<u8> {
        let mut tiles = Vec::new();
        for z in 0..=max_zoom {
            let n = 1u32 << z;
            for y in 0..n {
                for x in 0..n {
                    tiles.push(Tile { z, x, y, data: payload(z, x, y) });
                }
            }
        }
        let world =
            LonLatBounds { min_lon: -180.0, min_lat: -85.0, max_lon: 180.0, max_lat: 85.0 };
        build_png(tiles, &world, 0, max_zoom).unwrap()
    }

    /// Drives an extractor to completion against an in-memory archive,
    /// returning the output plus how many fetch round-trips happened.
    fn drive(archive: &[u8], opts: ExtractOptions) -> Result<(Vec<u8>, usize), PmtilesError> {
        let mut ex = Extractor::new(opts)?;
        let mut fetches = 0usize;
        while !ex.is_done() {
            let wants = ex.wanted();
            assert!(!wants.is_empty(), "not done but nothing wanted");
            for r in wants {
                fetches += 1;
                let start = r.offset as usize;
                let end = (r.offset + r.length).min(archive.len() as u64) as usize;
                assert!(start <= archive.len(), "request beyond EOF");
                ex.feed(r.offset, &archive[start..end])?;
            }
        }
        Ok((ex.finish()?, fetches))
    }

    /// Reads one tile out of an archive by walking root → leaf directories.
    fn get_tile(archive: &[u8], z: u8, x: u32, y: u32) -> Option<Vec<u8>> {
        let h = Header::parse(archive).unwrap();
        let want = zxy_to_tile_id(z, x, y);
        let mut dir = decompress(
            &archive[h.root_dir_offset as usize..(h.root_dir_offset + h.root_dir_length) as usize],
            h.internal_compression,
        )
        .unwrap();
        loop {
            let entries = parse_directory(&dir).unwrap();
            // Last entry whose tile_id <= want.
            let idx = entries.partition_point(|e| e.tile_id <= want).checked_sub(1)?;
            let e = entries[idx];
            if e.run_length == 0 {
                let start = (h.leaf_dirs_offset + e.offset) as usize;
                dir = decompress(
                    &archive[start..start + e.length as usize],
                    h.internal_compression,
                )
                .unwrap();
                continue;
            }
            if want >= e.tile_id + e.run_length as u64 {
                return None;
            }
            let start = (h.tile_data_offset + e.offset) as usize;
            return Some(archive[start..start + e.length as usize].to_vec());
        }
    }

    const QUADRANT: LonLatBounds =
        LonLatBounds { min_lon: 5.0, min_lat: 5.0, max_lon: 60.0, max_lat: 55.0 };

    #[test]
    fn extracts_exactly_the_bbox_tiles_with_identical_bytes() {
        let src = pyramid_archive(4);
        let opts = ExtractOptions::new(QUADRANT, 0, 4);
        let (out, _) = drive(&src, opts).unwrap();

        let h = Header::parse(&out).unwrap();
        assert_eq!(h.tile_type, TILETYPE_PNG);
        assert_eq!(h.tile_compression, COMPRESSION_NONE);
        assert_eq!(h.min_zoom, 0);
        assert_eq!(h.max_zoom, 4);
        assert!(h.clustered);

        for z in 0..=4u8 {
            let (x0, y0, x1, y1) = crate::tilemath::bbox_tile_rect(&QUADRANT, z);
            let n = 1u32 << z;
            for y in 0..n {
                for x in 0..n {
                    let inside = (x0..=x1).contains(&x) && (y0..=y1).contains(&y);
                    let got = get_tile(&out, z, x, y);
                    if inside {
                        assert_eq!(
                            got.as_deref(),
                            Some(payload(z, x, y).as_slice()),
                            "tile {z}/{x}/{y} should be present and identical"
                        );
                    } else {
                        assert_eq!(got, None, "tile {z}/{x}/{y} should be absent");
                    }
                }
            }
        }
    }

    #[test]
    fn zoom_subset_and_source_clamping() {
        let src = pyramid_archive(3);
        // Request z2..z9: z9 clamps to the source's max of 3.
        let opts = ExtractOptions::new(QUADRANT, 2, 9);
        let (out, _) = drive(&src, opts).unwrap();
        let h = Header::parse(&out).unwrap();
        assert_eq!((h.min_zoom, h.max_zoom), (2, 3));
        assert_eq!(get_tile(&out, 0, 0, 0), None, "z0 not requested");
        assert!(get_tile(&out, 2, 2, 1).is_some());
    }

    #[test]
    fn counts_and_metadata_survive_extraction() {
        let src = pyramid_archive(3);
        let (out, _) = drive(&src, ExtractOptions::new(QUADRANT, 0, 3)).unwrap();
        let h = Header::parse(&out).unwrap();
        let src_h = Header::parse(&src).unwrap();

        let mut expected: u64 = 0;
        for z in 0..=3u8 {
            let (x0, y0, x1, y1) = crate::tilemath::bbox_tile_rect(&QUADRANT, z);
            expected += u64::from(x1 - x0 + 1) * u64::from(y1 - y0 + 1);
        }
        assert_eq!(h.num_addressed_tiles, expected);

        let meta = decompress(
            &out[h.metadata_offset as usize..(h.metadata_offset + h.metadata_length) as usize],
            h.internal_compression,
        )
        .unwrap();
        let src_meta = decompress(
            &src[src_h.metadata_offset as usize
                ..(src_h.metadata_offset + src_h.metadata_length) as usize],
            src_h.internal_compression,
        )
        .unwrap();
        assert_eq!(meta, src_meta, "metadata must pass through unchanged");
    }

    #[test]
    fn extraction_walks_leaf_directories() {
        // z0..=6 = 5461 entries with non-clustered-friendly payloads still fits
        // a root; z0..=7 (21845) with unique payloads reliably splits.
        let src = pyramid_archive(7);
        assert!(
            Header::parse(&src).unwrap().leaf_dirs_length > 0,
            "test premise: source must have leaf directories"
        );
        let small =
            LonLatBounds { min_lon: 10.0, min_lat: 10.0, max_lon: 12.0, max_lat: 12.0 };
        let (out, _) = drive(&src, ExtractOptions::new(small, 0, 7)).unwrap();
        for z in 0..=7u8 {
            let (x, y) = crate::tilemath::lonlat_to_tile(11.0, 11.0, z);
            assert_eq!(
                get_tile(&out, z, x, y).as_deref(),
                Some(payload(z, x, y).as_slice()),
                "center tile {z}/{x}/{y}"
            );
        }
    }

    #[test]
    fn range_gap_merging_reduces_round_trips() {
        let src = pyramid_archive(5);
        let mut merged = ExtractOptions::new(QUADRANT, 0, 5);
        merged.max_range_gap = 1 << 20;
        let mut unmerged = ExtractOptions::new(QUADRANT, 0, 5);
        unmerged.max_range_gap = 0;
        let (out_a, fetches_a) = drive(&src, merged).unwrap();
        let (out_b, fetches_b) = drive(&src, unmerged).unwrap();
        assert!(fetches_a < fetches_b, "{fetches_a} !< {fetches_b}");
        assert_eq!(out_a, out_b, "gap merging must not change the output");
    }

    #[test]
    fn duplicate_blobs_stay_deduplicated_and_runs_reform() {
        // Every tile shares one payload → source dedups to a single blob and
        // (within each zoom) consecutive ids form runs.
        let mut tiles = Vec::new();
        for z in 0..=3u8 {
            let n = 1u32 << z;
            for y in 0..n {
                for x in 0..n {
                    tiles.push(Tile { z, x, y, data: b"same".to_vec() });
                }
            }
        }
        let world =
            LonLatBounds { min_lon: -180.0, min_lat: -85.0, max_lon: 180.0, max_lat: 85.0 };
        let src = build_png(tiles, &world, 0, 3).unwrap();

        let (out, _) = drive(&src, ExtractOptions::new(world, 0, 3)).unwrap();
        let h = Header::parse(&out).unwrap();
        assert_eq!(h.num_tile_contents, 1, "single shared blob");
        assert_eq!(h.num_addressed_tiles, 85, "full pyramid via runs");
        assert_eq!(h.tile_data_length, 4, "one copy of the blob");
        // The whole pyramid is one contiguous id range sharing one blob, so
        // the directory should collapse to a single run entry.
        let root = decompress(
            &out[h.root_dir_offset as usize..(h.root_dir_offset + h.root_dir_length) as usize],
            h.internal_compression,
        )
        .unwrap();
        let entries = parse_directory(&root).unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].run_length, 85);
    }

    #[test]
    fn empty_selection_is_a_clear_error() {
        let src = pyramid_archive(2);
        // Zooms outside the source range.
        let mut ex =
            Extractor::new(ExtractOptions::new(QUADRANT, 10, 12)).unwrap();
        let r = ex.wanted()[0];
        let err = ex
            .feed(r.offset, &src[..(r.length as usize).min(src.len())])
            .unwrap_err();
        assert!(matches!(err, PmtilesError::InvalidRequest(_)), "{err}");
    }

    #[test]
    fn tile_cap_aborts_early() {
        let src = pyramid_archive(2);
        let mut opts = ExtractOptions::new(QUADRANT, 0, 2);
        opts.max_tiles = 2;
        let mut ex = Extractor::new(opts).unwrap();
        let r = ex.wanted()[0];
        let err = ex
            .feed(r.offset, &src[..(r.length as usize).min(src.len())])
            .unwrap_err();
        assert!(matches!(err, PmtilesError::TooManyTiles { .. }), "{err}");
    }

    #[test]
    fn feed_rejects_unknown_and_short_ranges() {
        let src = pyramid_archive(2);
        let mut ex = Extractor::new(ExtractOptions::new(QUADRANT, 0, 2)).unwrap();
        assert!(matches!(
            ex.feed(999_999, &[0u8; 4]),
            Err(PmtilesError::Protocol(_))
        ));
        // Drive the prelude, then short-change the first follow-up request.
        let r = ex.wanted()[0];
        ex.feed(r.offset, &src[..(r.length as usize).min(src.len())])
            .unwrap();
        if let Some(next) = ex.wanted().first().copied() {
            if next.length > 1 {
                let err = ex
                    .feed(next.offset, &src[next.offset as usize..][..1])
                    .unwrap_err();
                assert!(matches!(err, PmtilesError::Protocol(_)), "{err}");
                // The request must still be outstanding after the bad feed.
                assert!(ex.wanted().iter().any(|w| w.offset == next.offset));
            }
        }
    }

    #[test]
    fn extract_of_extract_is_stable() {
        // Extracting the same bbox from an extract must reproduce every tile.
        let src = pyramid_archive(4);
        let (once, _) = drive(&src, ExtractOptions::new(QUADRANT, 0, 4)).unwrap();
        let (twice, _) = drive(&once, ExtractOptions::new(QUADRANT, 0, 4)).unwrap();
        for z in 0..=4u8 {
            let (x0, y0, x1, y1) = crate::tilemath::bbox_tile_rect(&QUADRANT, z);
            for y in y0..=y1 {
                for x in x0..=x1 {
                    assert_eq!(get_tile(&twice, z, x, y), get_tile(&once, z, x, y));
                }
            }
        }
    }

    #[test]
    fn progress_reports_phases_and_byte_counts() {
        let src = pyramid_archive(3);
        let mut ex = Extractor::new(ExtractOptions::new(QUADRANT, 0, 3)).unwrap();
        assert_eq!(ex.progress().phase, "header");
        let r = ex.wanted()[0];
        ex.feed(r.offset, &src[..(r.length as usize).min(src.len())])
            .unwrap();
        // Small archive: directories resolve straight from the prelude.
        let p = ex.progress();
        assert_eq!(p.phase, "data");
        assert!(p.tiles_selected > 0);
        assert!(p.data_bytes_total > 0);
        assert!(p.estimated_output_bytes > 0);
        for r in ex.wanted() {
            let start = r.offset as usize;
            ex.feed(r.offset, &src[start..start + r.length as usize]).unwrap();
        }
        let p = ex.progress();
        assert_eq!(p.phase, "done");
        assert_eq!(p.data_bytes_received, p.data_bytes_total);
    }
}
