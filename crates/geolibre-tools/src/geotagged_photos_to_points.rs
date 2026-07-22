//! GeoLibre tool: build a point layer from a folder of geotagged photos.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *GeoTagged Photos To Points* (Data
//! Management). Field crews return with folders of GPS-stamped phone/camera
//! JPEGs; this tool walks the folder, reads each photo's EXIF header (metadata
//! only — no pixels are decoded), converts the GPS latitude/longitude rationals
//! (degrees/minutes/seconds + N/S/E/W reference) to decimal degrees, and emits
//! one WGS-84 point per photo carrying the file path, capture timestamp,
//! altitude, camera bearing (GPS image direction), and camera make/model.
//!
//! There is no equivalent among the ~791 bundled whitebox-wasm tools. Output
//! points drop straight into `vector_to_pmtiles` / `render_vector_png` for a web
//! map of where each photo was taken.
//!
//! Params:
//! - `input`   — folder path of images (JPEG/TIFF). Required; must exist.
//! - `output`  — output point layer (driver from extension). Optional (memory).
//! - `recursive`      — descend into subfolders. Optional, default false.
//! - `only_geotagged` — skip images with no GPS position. Optional, default true.
//!   When false, non-geotagged images are still emitted with a null geometry so
//!   their path/timestamp/camera metadata is captured.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{parse_optional_str, write_or_store_layer};

pub struct GeotaggedPhotosToPointsTool;

impl Tool for GeotaggedPhotosToPointsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "geotagged_photos_to_points",
            display_name: "GeoTagged Photos To Points",
            summary: "Read a folder of geotagged photos and build a WGS-84 point layer with each photo's path, capture datetime, altitude, camera bearing (GPS image direction), and make/model — parsed from EXIF metadata only (no pixels decoded), like ArcGIS GeoTagged Photos To Points.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Folder path containing photos (JPEG/TIFF) to read EXIF GPS metadata from.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output point vector (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "recursive",
                    description: "Descend into subfolders when scanning for photos. Default false.",
                    required: false,
                },
                ToolParamSpec {
                    name: "only_geotagged",
                    description: "Skip images that carry no GPS position (default true). When false, non-geotagged images are emitted with a null geometry.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        let input = require_str(args, "input")?;
        let path = Path::new(input);
        if !path.exists() {
            return Err(ToolError::Validation(format!(
                "input folder '{input}' does not exist"
            )));
        }
        if !path.is_dir() {
            return Err(ToolError::Validation(format!(
                "input '{input}' is not a folder"
            )));
        }
        opt_bool(args, "recursive")?;
        opt_bool(args, "only_geotagged")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let recursive = opt_bool(args, "recursive")?.unwrap_or(false);
        let only_geotagged = opt_bool(args, "only_geotagged")?.unwrap_or(true);

        // Collect candidate image files (sorted for deterministic output).
        let mut files: Vec<PathBuf> = Vec::new();
        collect_images(Path::new(input), recursive, &mut files)?;
        files.sort();

        ctx.progress
            .info(&format!("scanning {} image file(s)", files.len()));

        let mut out = Layer::new("photos")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        out.add_field(FieldDef::new("path", FieldType::Text));
        out.add_field(FieldDef::new("name", FieldType::Text));
        out.add_field(FieldDef::new("datetime", FieldType::Text));
        out.add_field(FieldDef::new("longitude", FieldType::Float));
        out.add_field(FieldDef::new("latitude", FieldType::Float));
        out.add_field(FieldDef::new("altitude", FieldType::Float));
        out.add_field(FieldDef::new("direction", FieldType::Float));
        out.add_field(FieldDef::new("make", FieldType::Text));
        out.add_field(FieldDef::new("model", FieldType::Text));

        let mut point_count = 0usize;
        let mut skipped = 0usize;
        for file in &files {
            let meta = read_photo_meta(file);
            let name = file
                .file_name()
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_default();
            let path_str = file.to_string_lossy().to_string();

            let has_gps = meta.lon.is_some() && meta.lat.is_some();
            if !has_gps && only_geotagged {
                skipped += 1;
                continue;
            }

            let geometry = match (meta.lon, meta.lat) {
                (Some(lon), Some(lat)) => {
                    point_count += 1;
                    Some(Geometry::point(lon, lat))
                }
                _ => None,
            };

            out.push(Feature {
                fid: 0,
                geometry,
                attributes: vec![
                    FieldValue::Text(path_str),
                    FieldValue::Text(name),
                    opt_text(meta.datetime),
                    opt_float(meta.lon),
                    opt_float(meta.lat),
                    opt_float(meta.altitude),
                    opt_float(meta.direction),
                    opt_text(meta.make),
                    opt_text(meta.model),
                ],
            });
        }

        let out_path = write_or_store_layer(out, output)?;
        ctx.progress.info(&format!(
            "{point_count} geotagged point(s), {skipped} non-geotagged skipped"
        ));

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("photo_count".to_string(), json!(files.len()));
        outputs.insert("point_count".to_string(), json!(point_count));
        outputs.insert("skipped".to_string(), json!(skipped));
        Ok(ToolRunResult { outputs })
    }
}

// ── Folder walk ──────────────────────────────────────────────────────────────

fn is_image(path: &Path) -> bool {
    matches!(
        path.extension()
            .and_then(|e| e.to_str())
            .map(str::to_ascii_lowercase)
            .as_deref(),
        Some("jpg" | "jpeg" | "tif" | "tiff")
    )
}

/// Collects image files under `dir`, optionally recursing into subfolders.
fn collect_images(dir: &Path, recursive: bool, out: &mut Vec<PathBuf>) -> Result<(), ToolError> {
    let entries = std::fs::read_dir(dir).map_err(|e| {
        ToolError::Execution(format!("failed reading folder '{}': {e}", dir.display()))
    })?;
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            if recursive {
                collect_images(&path, recursive, out)?;
            }
        } else if is_image(&path) {
            out.push(path);
        }
    }
    Ok(())
}

// ── EXIF parsing ─────────────────────────────────────────────────────────────

#[derive(Default)]
struct PhotoMeta {
    lat: Option<f64>,
    lon: Option<f64>,
    altitude: Option<f64>,
    direction: Option<f64>,
    datetime: Option<String>,
    make: Option<String>,
    model: Option<String>,
}

/// Reads EXIF metadata from a photo file. Any read/parse failure (or absent
/// tags) yields an empty/partial `PhotoMeta` rather than an error, so one bad
/// file never aborts a whole folder.
fn read_photo_meta(path: &Path) -> PhotoMeta {
    use exif::{In, Tag};

    let mut meta = PhotoMeta::default();
    let Ok(bytes) = std::fs::read(path) else {
        return meta;
    };
    let mut cursor = std::io::Cursor::new(&bytes);
    // `continue_on_error(true)` lets us recover the fields already parsed when a
    // reader hits a non-fatal quirk (e.g. a non-zero sub-IFD "next IFD" pointer,
    // which some EXIF writers emit) instead of discarding the whole header.
    let mut reader = exif::Reader::new();
    reader.continue_on_error(true);
    let exif = match reader.read_from_container(&mut cursor) {
        Ok(e) => e,
        Err(e) => match e.distill_partial_result(|_errors| {}) {
            Ok(e) => e,
            Err(_) => return meta,
        },
    };

    // GPS latitude / longitude: 3 rationals (deg, min, sec) + N/S/E/W reference.
    let lat = exif.get_field(Tag::GPSLatitude, In::PRIMARY).and_then(dms);
    let lat_ref = exif
        .get_field(Tag::GPSLatitudeRef, In::PRIMARY)
        .and_then(ascii_str)
        .unwrap_or_default();
    let lon = exif.get_field(Tag::GPSLongitude, In::PRIMARY).and_then(dms);
    let lon_ref = exif
        .get_field(Tag::GPSLongitudeRef, In::PRIMARY)
        .and_then(ascii_str)
        .unwrap_or_default();
    if let (Some(lat), Some(lon)) = (lat, lon) {
        meta.lat = Some(apply_hemisphere(lat, &lat_ref));
        meta.lon = Some(apply_hemisphere(lon, &lon_ref));
    }

    // Altitude: a single rational; AltitudeRef byte 1 means below sea level.
    if let Some(alt) = exif
        .get_field(Tag::GPSAltitude, In::PRIMARY)
        .and_then(first_rational)
    {
        let below = exif
            .get_field(Tag::GPSAltitudeRef, In::PRIMARY)
            .and_then(first_byte)
            .map(|b| b == 1)
            .unwrap_or(false);
        meta.altitude = Some(if below { -alt } else { alt });
    }

    // GPS image direction (bearing the camera pointed), a single rational.
    meta.direction = exif
        .get_field(Tag::GPSImgDirection, In::PRIMARY)
        .and_then(first_rational);

    // Capture timestamp — prefer the original, fall back to the file datetime.
    meta.datetime = exif
        .get_field(Tag::DateTimeOriginal, In::PRIMARY)
        .and_then(ascii_str)
        .or_else(|| {
            exif.get_field(Tag::DateTime, In::PRIMARY)
                .and_then(ascii_str)
        });

    meta.make = exif.get_field(Tag::Make, In::PRIMARY).and_then(ascii_str);
    meta.model = exif.get_field(Tag::Model, In::PRIMARY).and_then(ascii_str);

    meta
}

/// Converts a GPS coordinate field (3 rationals: deg, min, sec) to unsigned
/// decimal degrees. Returns None if the value is not a 3-part rational.
fn dms(field: &exif::Field) -> Option<f64> {
    if let exif::Value::Rational(ref v) = field.value {
        if v.len() >= 3 {
            return Some(dms_to_decimal(v[0].to_f64(), v[1].to_f64(), v[2].to_f64()));
        }
    }
    None
}

/// Combines degrees/minutes/seconds into unsigned decimal degrees.
fn dms_to_decimal(deg: f64, min: f64, sec: f64) -> f64 {
    deg + min / 60.0 + sec / 3600.0
}

/// Applies the N/S/E/W hemisphere reference: S and W make the value negative.
fn apply_hemisphere(decimal: f64, reference: &str) -> f64 {
    let r = reference.trim().to_ascii_uppercase();
    if r.starts_with('S') || r.starts_with('W') {
        -decimal.abs()
    } else {
        decimal.abs()
    }
}

fn first_rational(field: &exif::Field) -> Option<f64> {
    match field.value {
        exif::Value::Rational(ref v) if !v.is_empty() => Some(v[0].to_f64()),
        exif::Value::SRational(ref v) if !v.is_empty() => Some(v[0].to_f64()),
        _ => None,
    }
}

fn first_byte(field: &exif::Field) -> Option<u8> {
    match field.value {
        exif::Value::Byte(ref v) if !v.is_empty() => Some(v[0]),
        _ => None,
    }
}

fn ascii_str(field: &exif::Field) -> Option<String> {
    if let exif::Value::Ascii(ref parts) = field.value {
        let bytes: Vec<u8> = parts.iter().flatten().copied().collect();
        let s = String::from_utf8_lossy(&bytes).trim().to_string();
        if s.is_empty() {
            None
        } else {
            Some(s)
        }
    } else {
        let s = field.display_value().to_string();
        if s.is_empty() {
            None
        } else {
            Some(s)
        }
    }
}

// ── Attribute helpers ────────────────────────────────────────────────────────

fn opt_text(v: Option<String>) -> FieldValue {
    match v {
        Some(s) => FieldValue::Text(s),
        None => FieldValue::Null,
    }
}

fn opt_float(v: Option<f64>) -> FieldValue {
    match v {
        Some(f) => FieldValue::Float(f),
        None => FieldValue::Null,
    }
}

// ── Parameter helpers ────────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

/// Parses an optional boolean: accepts a JSON bool or the strings
/// "true"/"false"/"1"/"0"/"yes"/"no" (host UIs often post strings).
fn opt_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "true" | "1" | "yes" | "y" => Ok(Some(true)),
            "false" | "0" | "no" | "n" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// The core numeric property: deg/min/sec + hemisphere → signed decimal
    /// degrees, to well within 1e-6 deg.
    #[test]
    fn dms_to_decimal_roundtrips() {
        // 40°26'46" N ≈ 40.446111
        let d = dms_to_decimal(40.0, 26.0, 46.0);
        assert!((d - 40.446111).abs() < 1e-6, "got {d}");
        assert!((apply_hemisphere(d, "N") - 40.446111).abs() < 1e-6);
        // Western/Southern hemispheres flip the sign.
        assert!((apply_hemisphere(d, "S") + 40.446111).abs() < 1e-6);
        assert!((apply_hemisphere(d, "W") + 40.446111).abs() < 1e-6);
        // A whole-degree value with fractional minutes.
        let d2 = dms_to_decimal(122.0, 30.0, 0.0);
        assert!((d2 - 122.5).abs() < 1e-12);
    }

    /// A single fractional decimal degree expressed as d/m/s survives the
    /// conversion at 1e-7 precision (what our synthesized fixtures rely on).
    #[test]
    fn dms_precision_is_tight() {
        // Embed 37.7749295 as deg/min/sec exactly.
        let total = 37.7749295_f64;
        let deg = total.trunc();
        let min_f = (total - deg) * 60.0;
        let min = min_f.trunc();
        let sec = (min_f - min) * 60.0;
        let back = dms_to_decimal(deg, min, sec);
        assert!((back - total).abs() < 1e-9, "got {back}");
    }

    /// End-to-end: an empty folder produces an empty point layer (and does not
    /// error), proving the folder walk + layer construction path.
    #[test]
    fn empty_folder_yields_zero_points() {
        let dir = std::env::temp_dir().join(format!(
            "geotag_test_{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let args: ToolArgs = serde_json::from_value(json!({
            "input": dir.to_string_lossy(),
        }))
        .unwrap();
        let out = GeotaggedPhotosToPointsTool.run(&args, &ctx()).unwrap();
        assert_eq!(out.outputs["point_count"], json!(0));
        assert_eq!(out.outputs["photo_count"], json!(0));
        std::fs::remove_dir_all(&dir).ok();
    }

    /// Non-image files are ignored by the folder scan.
    #[test]
    fn ignores_non_images() {
        assert!(is_image(Path::new("a/b/photo.JPG")));
        assert!(is_image(Path::new("shot.jpeg")));
        assert!(is_image(Path::new("scan.tiff")));
        assert!(!is_image(Path::new("notes.txt")));
        assert!(!is_image(Path::new("data.geojson")));
    }

    #[test]
    fn rejects_missing_or_bad_input() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            GeotaggedPhotosToPointsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "/no/such/folder/xyz123" })).is_err());
        // A file that exists but is not a directory is rejected.
        assert!(bad(json!({ "input": file!() })).is_err());
    }
}
