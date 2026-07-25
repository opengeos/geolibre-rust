//! Vector data: read common formats in memory and emit GeoJSON, with metadata
//! and optional reprojection. Backed by the pure-Rust `wbvector` engine.
//!
//! In-memory formats: GeoJSON, TopoJSON, GML, GPX, KML (text); FlatGeobuf,
//! GeoPackage, GeoParquet, KMZ (binary). Shapefile / MapInfo / OSM-PBF are
//! file-oriented in the engine and not yet exposed here.
//!
//! Every reader here feeds three output paths, all sharing [`read_layer`]:
//! GeoJSON text (this module), deck.gl binary attributes
//! ([`crate::binary`]), and Arrow IPC ([`crate::arrow_ipc`]).
use wasm_bindgen::prelude::*;
use wbvector::Layer;

fn jerr<E: std::fmt::Display>(ctx: &str) -> impl Fn(E) -> JsValue + '_ {
    move |e| JsValue::from_str(&format!("{ctx}: {e}"))
}

/// Parse a vector dataset (bytes + format name) into a `Layer`.
pub(crate) fn read_layer(data: &[u8], format: &str) -> Result<Layer, JsValue> {
    let as_text = || std::str::from_utf8(data)
        .map_err(|_| JsValue::from_str("text vector format requires valid UTF-8"));
    let layer = match format.to_lowercase().as_str() {
        "geojson" | "json" => wbvector::geojson::parse_str(as_text()?),
        "topojson" => wbvector::topojson::parse_str(as_text()?),
        "gml" => wbvector::gml::parse_str(as_text()?),
        "gpx" => wbvector::gpx::parse_str(as_text()?),
        "kml" => wbvector::kml::parse_str(as_text()?),
        "flatgeobuf" | "fgb" => wbvector::flatgeobuf::from_bytes(data),
        "geopackage" | "gpkg" => wbvector::geopackage::from_bytes(data.to_vec()),
        "kmz" => wbvector::kmz::from_bytes(data),
        // Upstream's GeoParquet reader is path-based, which the browser build
        // has no filesystem for; read it from the buffer instead.
        #[cfg(feature = "geoparquet")]
        "geoparquet" | "parquet" => {
            return crate::geoparquet_mem::from_bytes(data)
                .map_err(|e| JsValue::from_str(&format!("read: {e}")))
        }
        #[cfg(not(feature = "geoparquet"))]
        "geoparquet" | "parquet" => return Err(JsValue::from_str(
            "this build was compiled without the 'geoparquet' feature")),
        other => return Err(JsValue::from_str(&format!(
            "unsupported vector format '{other}' (try: {})", vector_formats()))),
    };
    layer.map_err(jerr("read"))
}

/// Apply the reprojection convention shared by every vector entry point:
/// `src_epsg` of `0` means "use the layer's own CRS", falling back to EPSG:4326
/// when it declares none (GeoJSON is WGS84 by RFC 7946).
pub(crate) fn reproject_layer(
    layer: Layer,
    dst_epsg: u32,
    src_epsg: u32,
) -> Result<Layer, JsValue> {
    let src = if src_epsg != 0 {
        src_epsg
    } else {
        layer.crs_epsg().unwrap_or(4326)
    };
    layer
        .reproject_from_to_epsg(src, dst_epsg)
        .map_err(jerr("reproject"))
}

/// Vector formats this build can read from memory (comma-separated).
#[wasm_bindgen]
pub fn vector_formats() -> String {
    let mut formats = "geojson,topojson,gml,gpx,kml,flatgeobuf,geopackage,kmz".to_string();
    if cfg!(feature = "geoparquet") {
        formats.push_str(",geoparquet");
    }
    formats
}

/// Read a vector dataset and return it as a GeoJSON `FeatureCollection` string.
#[wasm_bindgen]
pub fn vector_to_geojson(data: &[u8], format: &str) -> Result<String, JsValue> {
    Ok(wbvector::geojson::to_string(&read_layer(data, format)?))
}

/// Read a vector dataset, reproject it to `dst_epsg`, and return GeoJSON.
/// Uses the bundled pure-Rust projection engine (full EPSG support).
///
/// `src_epsg` overrides the source CRS: pass `0` to use the layer's own CRS, or
/// fall back to EPSG:4326 if it declares none (GeoJSON is WGS84 by RFC 7946).
#[wasm_bindgen]
pub fn vector_to_geojson_reproject(
    data: &[u8],
    format: &str,
    dst_epsg: u32,
    src_epsg: u32,
) -> Result<String, JsValue> {
    let layer = read_layer(data, format)?;
    let rep = reproject_layer(layer, dst_epsg, src_epsg)?;
    Ok(wbvector::geojson::to_string(&rep))
}

/// Read a vector dataset and return metadata as JSON:
/// `{"ok":true,"name","features","geometry","epsg"|null,"fields":[...],
///   "bbox":[min_x,min_y,max_x,max_y]|null}`.
#[wasm_bindgen]
pub fn vector_info(data: &[u8], format: &str) -> Result<String, JsValue> {
    let mut layer = read_layer(data, format)?;
    let count = layer.features.len();
    let geom = layer.geom_type
        .map(|g| format!("\"{g:?}\""))
        .unwrap_or_else(|| "null".into());
    let epsg = layer.crs_epsg().map(|e| e.to_string()).unwrap_or_else(|| "null".into());
    let fields: Vec<String> = layer.schema.fields().iter()
        .map(|f| format!("{{\"name\":\"{}\",\"type\":\"{:?}\"}}", f.name, f.field_type))
        .collect();
    let name = layer.name.replace('"', "'");
    let bbox = match layer.bbox() {
        Some(b) => format!("[{},{},{},{}]", b.min_x, b.min_y, b.max_x, b.max_y),
        None => "null".into(),
    };
    Ok(format!(
        "{{\"ok\":true,\"name\":\"{name}\",\"features\":{count},\"geometry\":{geom},\"epsg\":{epsg},\
\"fields\":[{}],\"bbox\":{bbox}}}",
        fields.join(",")
    ))
}
