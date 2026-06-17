//! WebAssembly bindings for the pure-Rust Whitebox GeoTIFF engine ([`wbgeotiff`]).
//!
//! The whole stack is pure Rust with no C dependencies (no GDAL/PROJ/HDF5),
//! so it compiles to `wasm32-unknown-unknown` and runs in any browser, Node,
//! Deno, or Wasmtime host. All functions operate on in-memory byte buffers.
//!
//! Two layers are exposed:
//! - **Convenience functions** ([`geotiff_info`], [`geotiff_stats`],
//!   [`geotiff_read_band_f64`]) for one-shot use.
//! - **Stateful classes** ([`GeoTiffReader`] parses once and serves many reads;
//!   [`CogBuilder`] encodes Cloud Optimized GeoTIFFs to bytes).
mod vector;
mod lidar;
mod analysis;
use wasm_bindgen::prelude::*;
use wbgeotiff::{CogLayout, CogWriter, Compression, GeoTiff, GeoTransform, SampleFormat};

/// Install a panic hook so Rust panics surface as readable `console.error`
/// messages instead of an opaque `RuntimeError: unreachable`.
#[wasm_bindgen(start)]
pub fn __start() {
    console_error_panic_hook::set_once();
}

// ───────────────────────── helpers ─────────────────────────

/// Format an `f64` as a JSON value: a finite number, or `null` for NaN /
/// infinity (which are not representable in JSON). Mirrors `JSON.stringify`.
fn json_f64(v: f64) -> String {
    if v.is_finite() { v.to_string() } else { "null".to_string() }
}

fn json_opt_f64(v: Option<f64>) -> String {
    v.map(json_f64).unwrap_or_else(|| "null".to_string())
}

fn err_json(msg: &str) -> String {
    format!("{{\"ok\":false,\"error\":\"{}\"}}", msg.replace('"', "'"))
}

fn jerr<E: std::fmt::Display>(ctx: &str) -> impl Fn(E) -> JsValue + '_ {
    move |e| JsValue::from_str(&format!("{ctx}: {e}"))
}

fn sample_format_str(sf: SampleFormat) -> String {
    format!("{:?}", sf).to_lowercase()
}

/// Axis-aligned bounds `[min_x, min_y, max_x, max_y]` (dataset CRS) from an
/// affine geo-transform and image size.
fn bbox_from_gt(g: &GeoTransform, width: u32, height: u32) -> [f64; 4] {
    let x0 = g.x_origin;
    let x1 = g.x_origin + width as f64 * g.pixel_width;
    let y0 = g.y_origin;
    let y1 = g.y_origin + height as f64 * g.pixel_height;
    [x0.min(x1), y0.min(y1), x0.max(x1), y0.max(y1)]
}

/// Build the source CRS for lon/lat conversion using the pure-Rust
/// `wbprojection` engine (full EPSG support, no PROJ/C). A user-defined
/// projection's PROJ string is preferred (e.g. NLCD's Albers, which has no EPSG
/// code and only reports its geographic base 4326); otherwise the EPSG code.
fn lonlat_crs(epsg: Option<u16>, proj_string: Option<&str>) -> Option<wbprojection::Crs> {
    if let Some(ps) = proj_string {
        if let Ok(c) = wbprojection::from_proj_string(ps) { return Some(c); }
    }
    wbprojection::Crs::from_epsg(epsg? as u32).ok()
}

/// Transform one point from `src` to WGS84 `(lon, lat)`, rejecting implausible
/// results (which catches a wrong/ambiguous CRS tag).
fn project_lonlat(src: &wbprojection::Crs, x: f64, y: f64) -> Option<(f64, f64)> {
    let wgs84 = wbprojection::Crs::wgs84_geographic();
    let (lon, lat) = src.transform_to(x, y, &wgs84).ok()?;
    if lon.is_finite() && lat.is_finite() && lon.abs() <= 180.000_001 && lat.abs() <= 90.000_001 {
        Some((lon, lat))
    } else {
        None
    }
}

/// WGS84 `(lon, lat)` of a single coordinate, or `None` if not convertible.
fn to_lonlat(epsg: Option<u16>, proj_string: Option<&str>, x: f64, y: f64) -> Option<(f64, f64)> {
    project_lonlat(&lonlat_crs(epsg, proj_string)?, x, y)
}

/// `[min_lon, min_lat, max_lon, max_lat]` (WGS84) for a CRS-native bbox.
///
/// The bbox is densified - corners, edge midpoints, and center are reprojected
/// and the min/max taken - so the geographic envelope is correct even for
/// projections where the extremes fall on edges rather than corners. Returns
/// empty if the CRS is not convertible.
fn bounds_lonlat(epsg: Option<u16>, proj_string: Option<&str>, b: [f64; 4]) -> Vec<f64> {
    let src = match lonlat_crs(epsg, proj_string) { Some(c) => c, None => return Vec::new() };
    let (x0, y0, x1, y1) = (b[0], b[1], b[2], b[3]);
    let (mx, my) = ((x0 + x1) / 2.0, (y0 + y1) / 2.0);
    let samples = [
        (x0, y0), (x1, y0), (x0, y1), (x1, y1), // corners
        (mx, y0), (mx, y1), (x0, my), (x1, my), // edge midpoints
        (mx, my),                               // center
    ];
    let (mut min_lon, mut min_lat, mut max_lon, mut max_lat) =
        (f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut any = false;
    for (x, y) in samples {
        if let Some((lon, lat)) = project_lonlat(&src, x, y) {
            any = true;
            min_lon = min_lon.min(lon); max_lon = max_lon.max(lon);
            min_lat = min_lat.min(lat); max_lat = max_lat.max(lat);
        }
    }
    if any { vec![min_lon, min_lat, max_lon, max_lat] } else { Vec::new() }
}

/// Largest full-raster allocation we will attempt, in bytes. WASM is 32-bit
/// (4 GiB linear memory) and the input file already occupies part of it, so we
/// cap whole-raster decodes well below that. Beyond this, callers get a clean
/// error instead of a `capacity overflow` panic (which would poison the module).
const MAX_FULL_READ_BYTES: u64 = 1_200_000_000;

/// Validate that decoding `width * height * bands` elements of `elem_bytes`
/// each will not overflow 32-bit `usize` or blow the memory budget.
fn guard_full_read(width: u32, height: u32, bands: usize, elem_bytes: usize) -> Result<(), String> {
    let cells = (width as u64)
        .checked_mul(height as u64)
        .and_then(|v| v.checked_mul(bands as u64));
    let bytes = cells.and_then(|c| c.checked_mul(elem_bytes as u64));
    match bytes {
        Some(b) if b <= MAX_FULL_READ_BYTES => Ok(()),
        _ => Err(format!(
            "raster too large to fully decode in 32-bit WASM: {}x{} x {} band(s) = {} cells. \
Use geotiff_info for metadata, or stream tiles with CogStream.",
            width, height, bands,
            cells.map(|c| c.to_string()).unwrap_or_else(|| "overflow".into())
        )),
    }
}

/// Compute min/max/mean/valid over band 0, returned as a JSON string.
fn stats_json_from(gt: &GeoTiff) -> String {
    let (w, h, bands) = (gt.width(), gt.height(), gt.band_count());
    let epsg = gt.epsg().map(|e| e.to_string()).unwrap_or_else(|| "null".into());
    if let Err(e) = guard_full_read(w, h, bands, 8) {
        return err_json(&e);
    }
    let values = match gt.read_all_f64() {
        Ok(v) => v,
        Err(e) => return err_json(&format!("read: {e}")),
    };
    let nodata = gt.no_data();
    let (mut min, mut max, mut sum, mut count) = (f64::INFINITY, f64::NEG_INFINITY, 0.0f64, 0u64);
    for &v in values.iter().step_by(bands.max(1)) {
        if v.is_nan() { continue; }
        if let Some(nd) = nodata { if v == nd { continue; } }
        if v < min { min = v; }
        if v > max { max = v; }
        sum += v;
        count += 1;
    }
    if count == 0 {
        return format!(
            "{{\"ok\":true,\"width\":{w},\"height\":{h},\"bands\":{bands},\"epsg\":{epsg},\
\"valid\":0,\"min\":null,\"max\":null,\"mean\":null}}"
        );
    }
    let mean = sum / count as f64;
    format!(
        "{{\"ok\":true,\"width\":{w},\"height\":{h},\"bands\":{bands},\"epsg\":{epsg},\
\"valid\":{count},\"min\":{},\"max\":{},\"mean\":{}}}",
        json_f64(min), json_f64(max), json_f64(mean)
    )
}

// ───────────────────────── convenience functions ─────────────────────────

/// Decode only a GeoTIFF's header and return its metadata as JSON. O(header)
/// memory, so it works on multi-gigabyte rasters that whole-image reads cannot
/// fit in WASM's 4 GiB address space.
///
/// `{"ok":true,"width","height","bands","epsg"|null,"nodata"|null,
///   "bits_per_sample","sample_format","compression","tiled","bigtiff"}`
#[wasm_bindgen]
pub fn geotiff_info(data: &[u8]) -> String {
    let m = match GeoTiff::peek_meta(data) {
        Ok(m) => m,
        Err(e) => return err_json(&format!("decode: {e}")),
    };
    let epsg = m.epsg.map(|e| e.to_string()).unwrap_or_else(|| "null".into());
    let (bbox, center, center_lonlat, bbox_lonlat) = match m.geo_transform.as_ref() {
        Some(g) => {
            let b = bbox_from_gt(g, m.width, m.height);
            let cx = (b[0] + b[2]) / 2.0;
            let cy = (b[1] + b[3]) / 2.0;
            let bbox = format!("[{},{},{},{}]", b[0], b[1], b[2], b[3]);
            let center = format!("[{cx},{cy}]");
            let ps = m.proj_string.as_deref();
            let cll = match to_lonlat(m.epsg, ps, cx, cy) {
                Some((lon, lat)) => format!("[{lon},{lat}]"),
                None => "null".into(),
            };
            let bll = match bounds_lonlat(m.epsg, ps, b).as_slice() {
                [a, c, d, e] => format!("[{a},{c},{d},{e}]"),
                _ => "null".into(),
            };
            (bbox, center, cll, bll)
        }
        None => ("null".into(), "null".into(), "null".into(), "null".into()),
    };
    format!(
        "{{\"ok\":true,\"width\":{},\"height\":{},\"bands\":{},\"epsg\":{},\"nodata\":{},\
\"bits_per_sample\":{},\"sample_format\":\"{}\",\"compression\":\"{:?}\",\"tiled\":{},\"bigtiff\":{},\
\"bbox\":{},\"center\":{},\"center_lonlat\":{},\"bbox_lonlat\":{}}}",
        m.width, m.height, m.bands, epsg, json_opt_f64(m.no_data),
        m.bits_per_sample, sample_format_str(m.sample_format), m.compression, m.tiled, m.is_bigtiff,
        bbox, center, center_lonlat, bbox_lonlat
    )
}

/// Decode a GeoTIFF and return band-0 summary statistics as JSON:
/// `{"ok":true,"width","height","bands","epsg","valid","min","max","mean"}`.
#[wasm_bindgen]
pub fn geotiff_stats(data: &[u8]) -> String {
    match GeoTiff::from_bytes(data) {
        Ok(gt) => stats_json_from(&gt),
        Err(e) => err_json(&format!("decode: {e}")),
    }
}

/// Read a single band of pixel values as an `f64` `Float64Array` (any on-disk
/// sample format is converted), row-major, length `width * height`.
#[wasm_bindgen]
pub fn geotiff_read_band_f64(data: &[u8], band: usize) -> Result<Vec<f64>, JsValue> {
    GeoTiffReader::new(data)?.read_band_f64(band)
}

/// Semantic version of this crate, exposed for runtime feature detection.
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

// ───────────────────────── GeoTiffReader ─────────────────────────

/// A parsed GeoTIFF held in memory. Construct once, then call the accessor and
/// `read_*` methods many times without re-parsing the file.
#[wasm_bindgen]
pub struct GeoTiffReader {
    inner: GeoTiff,
}

#[wasm_bindgen]
impl GeoTiffReader {
    /// Parse a GeoTIFF / BigTIFF / COG from raw bytes.
    #[wasm_bindgen(constructor)]
    pub fn new(data: &[u8]) -> Result<GeoTiffReader, JsValue> {
        let inner = GeoTiff::from_bytes(data).map_err(jerr("decode"))?;
        Ok(GeoTiffReader { inner })
    }

    // ── metadata accessors ──
    #[wasm_bindgen(getter)]
    pub fn width(&self) -> u32 { self.inner.width() }
    #[wasm_bindgen(getter)]
    pub fn height(&self) -> u32 { self.inner.height() }
    #[wasm_bindgen(getter)]
    pub fn bands(&self) -> usize { self.inner.band_count() }
    #[wasm_bindgen(getter)]
    pub fn bits_per_sample(&self) -> u16 { self.inner.bits_per_sample() }
    #[wasm_bindgen(getter)]
    pub fn sample_format(&self) -> String { sample_format_str(self.inner.sample_format()) }
    #[wasm_bindgen(getter)]
    pub fn compression(&self) -> String { format!("{:?}", self.inner.compression()) }
    #[wasm_bindgen(getter)]
    pub fn is_bigtiff(&self) -> bool { self.inner.is_bigtiff }
    /// EPSG code, or `undefined` if the file is not georeferenced by EPSG.
    #[wasm_bindgen(getter)]
    pub fn epsg(&self) -> Option<u32> { self.inner.epsg().map(|e| e as u32) }
    /// No-data sentinel, or `undefined` if none is declared.
    #[wasm_bindgen(getter)]
    pub fn nodata(&self) -> Option<f64> { self.inner.no_data() }

    /// Affine geo-transform as `[x_origin, pixel_width, row_rotation,
    /// y_origin, col_rotation, pixel_height]`, or an empty array if absent.
    pub fn geo_transform(&self) -> Vec<f64> {
        match self.inner.geo_transform() {
            Some(g) => vec![g.x_origin, g.pixel_width, g.row_rotation,
                            g.y_origin, g.col_rotation, g.pixel_height],
            None => Vec::new(),
        }
    }

    /// Bounding box as `[min_x, min_y, max_x, max_y]`, or empty if not georeferenced.
    pub fn bounding_box(&self) -> Vec<f64> {
        match self.inner.bounding_box() {
            Some(b) => vec![b.min_x, b.min_y, b.max_x, b.max_y],
            None => Vec::new(),
        }
    }

    /// Image center `[x, y]` in the dataset CRS, or empty if not georeferenced.
    pub fn center(&self) -> Vec<f64> {
        match self.inner.bounding_box() {
            Some(b) => vec![(b.min_x + b.max_x) / 2.0, (b.min_y + b.max_y) / 2.0],
            None => Vec::new(),
        }
    }

    /// Image center `[lon, lat]` in WGS84 degrees, or empty if not georeferenced
    /// or the CRS is not convertible.
    pub fn center_lonlat(&self) -> Vec<f64> {
        let b = match self.inner.bounding_box() { Some(b) => b, None => return Vec::new() };
        let cx = (b.min_x + b.max_x) / 2.0;
        let cy = (b.min_y + b.max_y) / 2.0;
        match to_lonlat(self.inner.epsg(), self.inner.proj_string().as_deref(), cx, cy) {
            Some((lon, lat)) => vec![lon, lat],
            None => Vec::new(),
        }
    }

    /// Bounds `[min_lon, min_lat, max_lon, max_lat]` in WGS84 degrees, or empty
    /// if not convertible.
    pub fn bounds_lonlat(&self) -> Vec<f64> {
        match self.inner.bounding_box() {
            Some(b) => bounds_lonlat(self.inner.epsg(), self.inner.proj_string().as_deref(),
                                     [b.min_x, b.min_y, b.max_x, b.max_y]),
            None => Vec::new(),
        }
    }

    /// GDAL value transform as `[scale, offset]` (physical = raw*scale+offset),
    /// or empty if none. Apply to `read_*` outputs to get physical values.
    pub fn value_transform(&self) -> Vec<f64> {
        match self.inner.value_transform() {
            Some(t) => vec![t.scale, t.offset],
            None => Vec::new(),
        }
    }

    /// Full metadata as a JSON string (same shape as [`geotiff_info`]).
    pub fn info_json(&self) -> String {
        let epsg = self.inner.epsg().map(|e| e.to_string()).unwrap_or_else(|| "null".into());
        format!(
            "{{\"ok\":true,\"width\":{},\"height\":{},\"bands\":{},\"epsg\":{},\"nodata\":{},\
\"bits_per_sample\":{},\"sample_format\":\"{}\",\"compression\":\"{:?}\"}}",
            self.inner.width(), self.inner.height(), self.inner.band_count(), epsg,
            json_opt_f64(self.inner.no_data()), self.inner.bits_per_sample(),
            sample_format_str(self.inner.sample_format()), self.inner.compression()
        )
    }

    /// Band-0 statistics as a JSON string (same shape as [`geotiff_stats`]).
    pub fn stats_json(&self) -> String { stats_json_from(&self.inner) }

    // ── data reads ──
    fn check_band(&self, band: usize) -> Result<(), JsValue> {
        let n = self.inner.band_count();
        if band >= n {
            return Err(JsValue::from_str(&format!("band {band} out of range (bands={n})")));
        }
        Ok(())
    }

    /// Read a band as `f64`, converting from any on-disk type. `Float64Array`.
    pub fn read_band_f64(&self, band: usize) -> Result<Vec<f64>, JsValue> {
        self.check_band(band)?;
        let bands = self.inner.band_count();
        guard_full_read(self.inner.width(), self.inner.height(), bands, 8)
            .map_err(|e| JsValue::from_str(&e))?;
        let all = self.inner.read_all_f64().map_err(jerr("read"))?;
        if bands == 1 { Ok(all) } else { Ok(all.into_iter().skip(band).step_by(bands).collect()) }
    }

    /// Read every band as `f64`, interleaved per pixel (`band0,band1,...`).
    pub fn read_all_f64(&self) -> Result<Vec<f64>, JsValue> {
        guard_full_read(self.inner.width(), self.inner.height(), self.inner.band_count(), 8)
            .map_err(|e| JsValue::from_str(&e))?;
        self.inner.read_all_f64().map_err(jerr("read"))
    }

    // Guard a native full-image decode (the engine decodes all tiles/strips at
    // `elem_bytes` per sample) against 32-bit overflow / memory blowup.
    fn guard_native(&self, band: usize, elem_bytes: usize) -> Result<(), JsValue> {
        self.check_band(band)?;
        guard_full_read(self.inner.width(), self.inner.height(), self.inner.band_count(), elem_bytes)
            .map_err(|e| JsValue::from_str(&e))
    }

    /// Read a band's raw, undecoded-to-native bytes. `Uint8Array`.
    pub fn read_band_bytes(&self, band: usize) -> Result<Vec<u8>, JsValue> {
        self.guard_native(band, (self.inner.bits_per_sample() as usize + 7) / 8)?;
        self.inner.read_band_bytes(band).map_err(jerr("read"))
    }

    // Native typed reads: these require the band's on-disk type to match and
    // error otherwise. Use `read_band_f64` for format-agnostic access.
    /// Native `u8` band. `Uint8Array`.
    pub fn read_band_u8(&self, band: usize) -> Result<Vec<u8>, JsValue> {
        self.guard_native(band, 1)?; self.inner.read_band_u8(band).map_err(jerr("read"))
    }
    /// Native `i8` band. `Int8Array`.
    pub fn read_band_i8(&self, band: usize) -> Result<Vec<i8>, JsValue> {
        self.guard_native(band, 1)?; self.inner.read_band_i8(band).map_err(jerr("read"))
    }
    /// Native `u16` band. `Uint16Array`.
    pub fn read_band_u16(&self, band: usize) -> Result<Vec<u16>, JsValue> {
        self.guard_native(band, 2)?; self.inner.read_band_u16(band).map_err(jerr("read"))
    }
    /// Native `i16` band. `Int16Array`.
    pub fn read_band_i16(&self, band: usize) -> Result<Vec<i16>, JsValue> {
        self.guard_native(band, 2)?; self.inner.read_band_i16(band).map_err(jerr("read"))
    }
    /// Native `u32` band. `Uint32Array`.
    pub fn read_band_u32(&self, band: usize) -> Result<Vec<u32>, JsValue> {
        self.guard_native(band, 4)?; self.inner.read_band_u32(band).map_err(jerr("read"))
    }
    /// Native `i32` band. `Int32Array`.
    pub fn read_band_i32(&self, band: usize) -> Result<Vec<i32>, JsValue> {
        self.guard_native(band, 4)?; self.inner.read_band_i32(band).map_err(jerr("read"))
    }
    /// Native `f32` band. `Float32Array`.
    pub fn read_band_f32(&self, band: usize) -> Result<Vec<f32>, JsValue> {
        self.guard_native(band, 4)?; self.inner.read_band_f32(band).map_err(jerr("read"))
    }
}

// ───────────────────────── CogBuilder (writer) ─────────────────────────

/// Builder for encoding a Cloud Optimized GeoTIFF (tiled, with overviews and
/// GDAL ghost metadata) to bytes. A COG is also a valid plain GeoTIFF.
///
/// Configure with the `set_*` methods, then call one of `write_*` with the
/// pixel data to get a `Uint8Array` of the encoded file.
#[wasm_bindgen]
pub struct CogBuilder {
    width: u32,
    height: u32,
    bands: u16,
    epsg: Option<u16>,
    nodata: Option<f64>,
    geo_transform: Option<[f64; 6]>,
    compression: Compression,
    tile_size: Option<u32>,
    bigtiff: bool,
    overview_levels: Option<Vec<u32>>,
}

#[wasm_bindgen]
impl CogBuilder {
    /// New builder for a `width` x `height` raster with `bands` bands.
    #[wasm_bindgen(constructor)]
    pub fn new(width: u32, height: u32, bands: u16) -> CogBuilder {
        CogBuilder {
            width, height, bands: bands.max(1),
            epsg: None, nodata: None, geo_transform: None,
            compression: Compression::Deflate, tile_size: None,
            bigtiff: false, overview_levels: None,
        }
    }

    /// Set the EPSG code (1..=65535).
    pub fn set_epsg(&mut self, epsg: u32) { if (1..=65535).contains(&epsg) { self.epsg = Some(epsg as u16); } }
    /// Set the no-data sentinel value.
    pub fn set_nodata(&mut self, v: f64) { self.nodata = Some(v); }
    /// Set the full affine geo-transform:
    /// `[x_origin, pixel_width, row_rotation, y_origin, col_rotation, pixel_height]`.
    pub fn set_geo_transform(&mut self, gt: Vec<f64>) -> Result<(), JsValue> {
        if gt.len() != 6 { return Err(JsValue::from_str("geo_transform needs 6 values")); }
        self.geo_transform = Some([gt[0], gt[1], gt[2], gt[3], gt[4], gt[5]]);
        Ok(())
    }
    /// Convenience: north-up geo-transform from upper-left origin and pixel size.
    pub fn set_origin(&mut self, x_min: f64, y_max: f64, pixel_size: f64) {
        self.geo_transform = Some([x_min, pixel_size, 0.0, y_max, 0.0, -pixel_size]);
    }
    /// Compression: `none`, `lzw`, `deflate`, `packbits`, `webp`, `jpeg`, `jpegxl`.
    pub fn set_compression(&mut self, name: &str) -> Result<(), JsValue> {
        self.compression = match name.to_lowercase().as_str() {
            "none" => Compression::None,
            "lzw" => Compression::Lzw,
            "deflate" | "zip" => Compression::Deflate,
            "packbits" => Compression::PackBits,
            "webp" => Compression::WebP,
            "jpeg" => Compression::Jpeg,
            "jpegxl" | "jxl" => Compression::JpegXl,
            other => return Err(JsValue::from_str(&format!("unknown compression: {other}"))),
        };
        Ok(())
    }
    /// Internal tile size in pixels (default 512).
    pub fn set_tile_size(&mut self, px: u32) { self.tile_size = Some(px); }
    /// Force BigTIFF (64-bit offsets) for very large outputs.
    pub fn set_bigtiff(&mut self, on: bool) { self.bigtiff = on; }
    /// Explicit overview decimation factors (e.g. `[2,4,8]`); empty disables overviews.
    pub fn set_overview_levels(&mut self, levels: Vec<u32>) { self.overview_levels = Some(levels); }

    fn build(&self) -> CogWriter {
        let mut w = CogWriter::new(self.width, self.height, self.bands)
            .compression(self.compression);
        if let Some(e) = self.epsg { w = w.epsg(e); }
        if let Some(nd) = self.nodata { w = w.no_data(nd); }
        if let Some(g) = self.geo_transform {
            w = w.geo_transform(GeoTransform::new(g[0], g[1], g[2], g[3], g[4], g[5]));
        }
        if let Some(ts) = self.tile_size { w = w.tile_size(ts); }
        if self.bigtiff { w = w.bigtiff(true); }
        if let Some(ref lv) = self.overview_levels { w = w.overview_levels(lv.clone()); }
        w
    }

    fn check_len(&self, len: usize) -> Result<(), JsValue> {
        let expected = self.width as usize * self.height as usize * self.bands as usize;
        if len != expected {
            return Err(JsValue::from_str(&format!(
                "data length {len} != width*height*bands {expected}")));
        }
        Ok(())
    }

    /// Encode `u8` pixel data to a COG. `Uint8Array`.
    pub fn write_u8(&self, data: &[u8]) -> Result<Vec<u8>, JsValue> {
        self.check_len(data.len())?;
        self.build().write_u8_to_vec(data).map_err(jerr("write"))
    }
    /// Encode `f32` pixel data to a COG. `Uint8Array`.
    pub fn write_f32(&self, data: &[f32]) -> Result<Vec<u8>, JsValue> {
        self.check_len(data.len())?;
        self.build().write_f32_to_vec(data).map_err(jerr("write"))
    }
    /// Encode `f64` pixel data to a COG. `Uint8Array`.
    pub fn write_f64(&self, data: &[f64]) -> Result<Vec<u8>, JsValue> {
        self.check_len(data.len())?;
        self.build().write_f64_to_vec(data).map_err(jerr("write"))
    }
}

// ───────────────────────── CogStream (range-request reading) ─────────────────────────

/// Range-request reader for a (tiled) Cloud Optimized GeoTIFF.
///
/// The wasm module does no network I/O itself; this class parses the header and
/// tells the JS host exactly which byte ranges to fetch, then decodes the tiles
/// the host fetches. Typical flow:
///
/// 1. Range-fetch the first chunk of the file (e.g. 0..1 MiB) and
///    `new CogStream(headerBytes)`. If it throws "need more header bytes", fetch
///    a larger prefix and retry.
/// 2. Pick a level (0 = full res, higher = overviews) and a pixel window.
/// 3. `tiles_for_window(level, x, y, w, h)` returns the tiles and their byte
///    ranges; range-fetch each, then `decode_tile_f64(level, bytes)`.
#[wasm_bindgen]
pub struct CogStream {
    layout: CogLayout,
}

#[wasm_bindgen]
impl CogStream {
    /// Parse a COG's tile layout from front-of-file header bytes.
    #[wasm_bindgen(constructor)]
    pub fn new(header_bytes: &[u8]) -> Result<CogStream, JsValue> {
        let layout = GeoTiff::parse_cog_layout(header_bytes).map_err(jerr("header"))?;
        Ok(CogStream { layout })
    }

    /// Number of resolution levels (1 + overview count).
    #[wasm_bindgen(getter)]
    pub fn num_levels(&self) -> usize { self.layout.levels.len() }
    /// EPSG code of the full-resolution level, if any.
    #[wasm_bindgen(getter)]
    pub fn epsg(&self) -> Option<u32> { self.layout.epsg.map(|e| e as u32) }
    /// No-data sentinel, if declared.
    #[wasm_bindgen(getter)]
    pub fn nodata(&self) -> Option<f64> { self.layout.no_data }

    /// Level-0 geo-transform `[x_origin, pixel_width, row_rot, y_origin, col_rot,
    /// pixel_height]`, or empty if not georeferenced.
    pub fn geo_transform(&self) -> Vec<f64> {
        match self.layout.geo_transform.as_ref() {
            Some(g) => vec![g.x_origin, g.pixel_width, g.row_rotation,
                            g.y_origin, g.col_rotation, g.pixel_height],
            None => Vec::new(),
        }
    }

    // Full-resolution bbox (dataset CRS) from the level-0 geo-transform.
    fn bbox(&self) -> Option<[f64; 4]> {
        let g = self.layout.geo_transform.as_ref()?;
        let l0 = self.layout.levels.first()?;
        Some(bbox_from_gt(g, l0.width, l0.height))
    }

    /// Bounding box `[min_x, min_y, max_x, max_y]` in the dataset CRS, or empty.
    pub fn bounding_box(&self) -> Vec<f64> {
        self.bbox().map(|b| b.to_vec()).unwrap_or_default()
    }

    /// Image center `[x, y]` in the dataset CRS, or empty.
    pub fn center(&self) -> Vec<f64> {
        match self.bbox() {
            Some(b) => vec![(b[0] + b[2]) / 2.0, (b[1] + b[3]) / 2.0],
            None => Vec::new(),
        }
    }

    /// Image center `[lon, lat]` in WGS84 degrees, or empty if not convertible.
    pub fn center_lonlat(&self) -> Vec<f64> {
        let ps = self.layout.proj_string.as_deref();
        match self.bbox() {
            Some(b) => match to_lonlat(self.layout.epsg, ps, (b[0] + b[2]) / 2.0, (b[1] + b[3]) / 2.0) {
                Some((lon, lat)) => vec![lon, lat],
                None => Vec::new(),
            },
            None => Vec::new(),
        }
    }

    /// Bounds `[min_lon, min_lat, max_lon, max_lat]` in WGS84 degrees, or empty.
    pub fn bounds_lonlat(&self) -> Vec<f64> {
        match self.bbox() {
            Some(b) => bounds_lonlat(self.layout.epsg, self.layout.proj_string.as_deref(), b),
            None => Vec::new(),
        }
    }

    /// JSON array describing every level: `[{level,width,height,tile_width,
    /// tile_height,tiles_x,tiles_y,bands,bits_per_sample,sample_format,compression}]`.
    pub fn levels_json(&self) -> String {
        let mut parts = Vec::with_capacity(self.layout.levels.len());
        for (i, lv) in self.layout.levels.iter().enumerate() {
            parts.push(format!(
                "{{\"level\":{},\"width\":{},\"height\":{},\"tile_width\":{},\"tile_height\":{},\
\"tiles_x\":{},\"tiles_y\":{},\"bands\":{},\"bits_per_sample\":{},\"sample_format\":\"{}\",\"compression\":\"{:?}\"}}",
                i, lv.width, lv.height, lv.tile_width, lv.tile_height, lv.tiles_x, lv.tiles_y,
                lv.samples_per_pixel, lv.bits_per_sample, sample_format_str(lv.sample_format), lv.compression
            ));
        }
        format!("[{}]", parts.join(","))
    }

    fn level(&self, level: usize) -> Result<&wbgeotiff::CogLevel, JsValue> {
        self.layout.levels.get(level)
            .ok_or_else(|| JsValue::from_str(&format!("level {level} out of range (levels={})", self.layout.levels.len())))
    }

    /// `[offset, length]` byte range of the tile at `(col, row)` on `level`.
    pub fn tile_range(&self, level: usize, col: u32, row: u32) -> Result<Vec<f64>, JsValue> {
        let lv = self.level(level)?;
        let (off, len) = lv.tile_range(col, row)
            .ok_or_else(|| JsValue::from_str(&format!("tile ({col},{row}) out of range")))?;
        Ok(vec![off as f64, len as f64])
    }

    /// Tiles covering a pixel window on `level`, as a JSON array of
    /// `{col,row,offset,length}`. Fetch each byte range, then `decode_tile_f64`.
    pub fn tiles_for_window(&self, level: usize, x: u32, y: u32, w: u32, h: u32) -> Result<String, JsValue> {
        let lv = self.level(level)?;
        if w == 0 || h == 0 { return Ok("[]".to_string()); }
        let x1 = (x.saturating_add(w).min(lv.width)).saturating_sub(1);
        let y1 = (y.saturating_add(h).min(lv.height)).saturating_sub(1);
        if x >= lv.width || y >= lv.height {
            return Err(JsValue::from_str("window origin outside image"));
        }
        let c0 = x / lv.tile_width;
        let c1 = x1 / lv.tile_width;
        let r0 = y / lv.tile_height;
        let r1 = y1 / lv.tile_height;
        let mut parts = Vec::new();
        for row in r0..=r1 {
            for col in c0..=c1 {
                if let Some((off, len)) = lv.tile_range(col, row) {
                    parts.push(format!(
                        "{{\"col\":{col},\"row\":{row},\"offset\":{off},\"length\":{len}}}"));
                }
            }
        }
        Ok(format!("[{}]", parts.join(",")))
    }

    /// Decode one tile's fetched (compressed) bytes into an `f64` `Float64Array`,
    /// pixel-interleaved, length `tile_width * tile_height * bands`. Edge tiles
    /// come back full-size; clip to the image/window on the JS side.
    pub fn decode_tile_f64(&self, level: usize, tile_bytes: &[u8]) -> Result<Vec<f64>, JsValue> {
        self.level(level)?.decode_tile_f64(tile_bytes).map_err(jerr("decode tile"))
    }
}
