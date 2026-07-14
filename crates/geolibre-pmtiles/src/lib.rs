//! Sans-IO PMTiles v3 core for GeoLibre.
//!
//! Implements the subset of the [PMTiles v3 spec] GeoLibre needs, with **no
//! I/O of its own** so the same code drives every binding:
//!
//! - [`format`] — header, directory, and varint encode/decode plus internal
//!   (gzip) compression.
//! - [`tilemath`] — ZXY ⇄ PMTiles (Hilbert) tile ids and bbox → tile-id ranges.
//! - [`writer`] — assembles complete clustered archives, splitting directories
//!   into leaves when the root outgrows the spec's 16 KiB prelude.
//! - [`extract`] — the [`extract::Extractor`] state machine that pulls a
//!   bbox/zoom subset out of a (possibly enormous) remote archive using nothing
//!   but byte ranges the host fetches: browser `fetch` behind wasm-bindgen,
//!   reqwest in a native shell, or plain file reads in tests.
//!
//! [PMTiles v3 spec]: https://github.com/protomaps/PMTiles/blob/main/spec/v3/spec.md

pub mod extract;
pub mod format;
pub mod tilemath;
pub mod writer;

/// Geographic bounds in degrees (WGS84).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LonLatBounds {
    pub min_lon: f64,
    pub min_lat: f64,
    pub max_lon: f64,
    pub max_lat: f64,
}

/// Errors from parsing, planning, or assembling PMTiles data.
#[derive(Debug, Clone, PartialEq)]
pub enum PmtilesError {
    /// The buffer does not start with the PMTiles magic.
    NotPmtiles,
    /// Archive spec version other than 3.
    UnsupportedVersion(u8),
    /// Internal compression this crate cannot decode (brotli, zstd, unknown).
    UnsupportedCompression(u8),
    /// The requested bbox/zoom selection is invalid.
    InvalidRequest(String),
    /// The selection would address more tiles than the configured cap.
    TooManyTiles { requested: u64, max: u64 },
    /// Bytes ended before a complete structure could be read.
    Truncated(&'static str),
    /// The host fed bytes that do not match an outstanding request, or called
    /// a method out of order.
    Protocol(String),
    /// The archive violates the spec.
    Corrupt(String),
    /// (De)compression failure.
    Codec(String),
}

impl std::fmt::Display for PmtilesError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NotPmtiles => write!(f, "not a PMTiles archive (bad magic)"),
            Self::UnsupportedVersion(v) => write!(f, "unsupported PMTiles version {v} (only v3)"),
            Self::UnsupportedCompression(c) => write!(
                f,
                "unsupported internal compression {c} (only none=1 and gzip=2)"
            ),
            Self::InvalidRequest(m) => write!(f, "invalid extract request: {m}"),
            Self::TooManyTiles { requested, max } => write!(
                f,
                "selection addresses {requested} tiles, over the cap of {max}; \
                 shrink the bbox or lower max_zoom"
            ),
            Self::Truncated(what) => write!(f, "truncated {what}"),
            Self::Protocol(m) => write!(f, "protocol error: {m}"),
            Self::Corrupt(m) => write!(f, "corrupt archive: {m}"),
            Self::Codec(m) => write!(f, "compression error: {m}"),
        }
    }
}

impl std::error::Error for PmtilesError {}
