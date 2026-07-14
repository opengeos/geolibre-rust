//! PMTiles bbox extraction over host-driven range requests.
//!
//! Mirrors the [`crate::CogStream`] pattern: the wasm module does **no network
//! I/O**. [`PmtilesExtractor`] publishes the absolute byte ranges it needs;
//! the JS host `fetch`es them (with `Range:` headers, any order, any
//! concurrency) and feeds the bytes back. Typical flow:
//!
//! ```js
//! const ex = new PmtilesExtractor(minLon, minLat, maxLon, maxLat, 0, 15);
//! while (!ex.done) {
//!   const wants = JSON.parse(ex.wanted_json());
//!   await Promise.all(wants.map(async ({ offset, length }) => {
//!     const res = await fetch(url, { headers: {
//!       Range: `bytes=${offset}-${offset + length - 1}` } });
//!     ex.feed(offset, new Uint8Array(await res.arrayBuffer()));
//!   }));
//!   onProgress?.(JSON.parse(ex.progress_json()));
//! }
//! const archive = ex.finish(); // Uint8Array: a complete .pmtiles file
//! ```
//!
//! The output is a self-contained clustered PMTiles v3 archive of every tile
//! intersecting the bbox across the zoom range, with the source's tile type,
//! tile compression, and metadata carried through — suitable for offline use.

use wasm_bindgen::prelude::*;

use geolibre_pmtiles::extract::{ExtractOptions, Extractor};
use geolibre_pmtiles::LonLatBounds;

use crate::jerr;

#[wasm_bindgen]
pub struct PmtilesExtractor {
    inner: Option<Extractor>,
}

impl PmtilesExtractor {
    fn inner(&self) -> Result<&Extractor, JsValue> {
        self.inner
            .as_ref()
            .ok_or_else(|| JsValue::from_str("extractor already finished"))
    }
}

#[wasm_bindgen]
impl PmtilesExtractor {
    /// Plan an extraction of `min_zoom..=max_zoom` tiles intersecting the
    /// WGS84 bbox. Zooms are clamped to what the source archive contains once
    /// its header arrives; `min_zoom` 0 keeps the basemap usable zoomed out.
    #[wasm_bindgen(constructor)]
    pub fn new(
        min_lon: f64,
        min_lat: f64,
        max_lon: f64,
        max_lat: f64,
        min_zoom: u8,
        max_zoom: u8,
    ) -> Result<PmtilesExtractor, JsValue> {
        let bbox = LonLatBounds { min_lon, min_lat, max_lon, max_lat };
        let inner = Extractor::new(ExtractOptions::new(bbox, min_zoom, max_zoom))
            .map_err(jerr("extract"))?;
        Ok(PmtilesExtractor { inner: Some(inner) })
    }

    /// Cap on addressed tiles (default 2,000,000). Raise for huge desktop
    /// extracts; lower to fail fast in memory-constrained embeds.
    pub fn set_max_tiles(&mut self, max_tiles: f64) -> Result<(), JsValue> {
        match self.inner.as_mut() {
            Some(ex) => {
                ex.set_max_tiles(max_tiles as u64);
                Ok(())
            }
            None => Err(JsValue::from_str("extractor already finished")),
        }
    }

    /// Coalesce tile-data requests whose byte gap is at most this (default
    /// 65,536). Larger values trade overfetch for fewer HTTP round-trips.
    pub fn set_max_range_gap(&mut self, max_gap: f64) -> Result<(), JsValue> {
        match self.inner.as_mut() {
            Some(ex) => {
                ex.set_max_range_gap(max_gap as u64);
                Ok(())
            }
            None => Err(JsValue::from_str("extractor already finished")),
        }
    }

    /// Outstanding byte ranges the host should fetch, as a JSON array of
    /// `{"offset":n,"length":n}`. Empty array when nothing is outstanding.
    pub fn wanted_json(&self) -> Result<String, JsValue> {
        let parts: Vec<String> = self
            .inner()?
            .wanted()
            .iter()
            .map(|r| format!("{{\"offset\":{},\"length\":{}}}", r.offset, r.length))
            .collect();
        Ok(format!("[{}]", parts.join(",")))
    }

    /// Hand back the bytes of one `wanted_json()` range, identified by its
    /// offset. Ranges may be fed in any order.
    pub fn feed(&mut self, offset: f64, bytes: &[u8]) -> Result<(), JsValue> {
        self.inner
            .as_mut()
            .ok_or_else(|| JsValue::from_str("extractor already finished"))?
            .feed(offset as u64, bytes)
            .map_err(jerr("extract"))
    }

    /// True once every needed range has been fed; `finish()` is then valid.
    #[wasm_bindgen(getter)]
    pub fn done(&self) -> bool {
        self.inner.as_ref().is_some_and(|ex| ex.is_done())
    }

    /// Source archive header as JSON (`{}` until the first feed): zooms,
    /// bounds, tile type/compression, tile counts. Lets a UI validate the
    /// request and describe the source before committing to the download.
    pub fn header_json(&self) -> Result<String, JsValue> {
        let ex = self.inner()?;
        Ok(match ex.header() {
            None => "{}".to_string(),
            Some(h) => {
                let b = h.bounds();
                format!(
                    "{{\"tile_type\":{},\"tile_compression\":{},\"min_zoom\":{},\"max_zoom\":{},\
                     \"num_addressed_tiles\":{},\"bounds\":[{},{},{},{}]}}",
                    h.tile_type,
                    h.tile_compression,
                    h.min_zoom,
                    h.max_zoom,
                    h.num_addressed_tiles,
                    b.min_lon,
                    b.min_lat,
                    b.max_lon,
                    b.max_lat,
                )
            }
        })
    }

    /// Progress as JSON: `{"phase":"header|directories|data|done",
    /// "tiles_selected":n,"blobs_total":n,"data_bytes_total":n,
    /// "data_bytes_received":n,"estimated_output_bytes":n}`.
    pub fn progress_json(&self) -> Result<String, JsValue> {
        let p = self.inner()?.progress();
        Ok(format!(
            "{{\"phase\":\"{}\",\"tiles_selected\":{},\"blobs_total\":{},\
             \"data_bytes_total\":{},\"data_bytes_received\":{},\
             \"estimated_output_bytes\":{}}}",
            p.phase,
            p.tiles_selected,
            p.blobs_total,
            p.data_bytes_total,
            p.data_bytes_received,
            p.estimated_output_bytes,
        ))
    }

    /// Assemble the extracted archive. Consumes the extractor's buffers; the
    /// returned `Uint8Array` is a complete `.pmtiles` file.
    pub fn finish(&mut self) -> Result<Vec<u8>, JsValue> {
        self.inner
            .take()
            .ok_or_else(|| JsValue::from_str("extractor already finished"))?
            .finish()
            .map_err(jerr("extract"))
    }
}
