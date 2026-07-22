//! GeoLibre tool: convert coordinate strings between notations.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Convert Coordinate Notation*. Reads a
//! coordinate from each feature — either the point geometry (interpreted as
//! decimal-degree lon/lat) or a text field in some notation — and writes a new
//! field holding the same location in a target notation:
//!
//! - **DD**  decimal degrees            (`"38.897700, -77.036500"`)
//! - **DMS** degrees-minutes-seconds    (`"38°53'51.720\"N 077°02'11.400\"W"`)
//! - **DDM** degrees-decimal-minutes    (`"38°53.862'N 077°02.190'W"`)
//! - **UTM** Universal Transverse Mercator (`"18N 323480.000 4306479.000"`)
//! - **MGRS** Military Grid Reference System / USNG (`"18SUJ2348006479"`)
//!
//! No PROJ and no new dependency: UTM uses the closed-form Krüger transverse
//! Mercator series on the WGS84 ellipsoid (a = 6378137, f = 1/298.257223563),
//! and MGRS adds the grid-zone/latitude-band and 100km-square lettering plus the
//! standard minimum-northing table to resolve the 2,000,000m northing ambiguity
//! on the way back. Nothing in whitebox-wasm or the catalog converts coordinate
//! *strings* — the projection tools reproject whole rasters/layers.
//!
//! Scope for v1: single combined coordinate field or point geometry as the input
//! source (no separate X/Y-field mode); GARS/GEOREF notations are not supported.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct ConvertCoordinateNotationTool;

impl Tool for ConvertCoordinateNotationTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "convert_coordinate_notation",
            display_name: "Convert Coordinate Notation",
            summary: "Convert each feature's coordinate between notations — decimal degrees (DD), degrees-minutes-seconds (DMS), degrees-decimal-minutes (DDM), UTM, and MGRS/USNG — writing a new field in the target notation, like ArcGIS Convert Coordinate Notation. Pure grid math (Krüger transverse-Mercator + MGRS lettering on WGS84); no PROJ.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point (or table) vector layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "input_notation",
                    description: "Notation of the input coordinate: DD, DMS, DDM, UTM, or MGRS.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output_notation",
                    description: "Target notation to write: DD, DMS, DDM, UTM, or MGRS.",
                    required: true,
                },
                ToolParamSpec {
                    name: "coord_field",
                    description: "Field holding the input coordinate string. If omitted, the point geometry is used (valid only when input_notation=DD; geometry x=lon, y=lat).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_field",
                    description: "Name of the field to write. Default: the lowercased output notation (e.g. 'utm').",
                    required: false,
                },
                ToolParamSpec {
                    name: "precision",
                    description: "MGRS digits per axis (1=10km … 5=1m, up to 9). Default 5. Ignored by other notations.",
                    required: false,
                },
                ToolParamSpec {
                    name: "update_geometry",
                    description: "If true, rebuild each feature's geometry as a lon/lat point in EPSG:4326 from the parsed coordinate. Default false (keep input geometry/CRS).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;

        let coord_idx =
            match &prm.coord_field {
                Some(name) => Some(layer.schema.field_index(name).ok_or_else(|| {
                    ToolError::Validation(format!("coord_field '{name}' not found"))
                })?),
                None => None,
            };

        // Build the output layer: copy schema, add the new field, optionally
        // rebuild point geometry as lon/lat.
        let mut out = Layer::new("converted");
        if prm.update_geometry {
            out = out.with_geom_type(GeometryType::Point).with_crs_epsg(4326);
        } else {
            out.geom_type = layer.geom_type;
            if let Some(epsg) = layer.crs_epsg() {
                out = out.with_crs_epsg(epsg);
            }
        }
        for f in layer.schema.fields() {
            out.add_field(f.clone());
        }
        out.add_field(FieldDef::new(&prm.output_field, FieldType::Text));

        let mut converted = 0usize;
        let mut skipped = 0usize;
        for feature in layer.iter() {
            // Read the source coordinate as (lat, lon).
            let latlon = match coord_idx {
                Some(idx) => feature
                    .attributes
                    .get(idx)
                    .and_then(FieldValue::as_str)
                    .and_then(|s| parse_coord(s, prm.input_notation)),
                None => feature.geometry.as_ref().and_then(point_xy).map(|(x, y)| {
                    // Geometry is decimal-degree lon/lat.
                    (y, x)
                }),
            };
            let Some((lat, lon)) = latlon else {
                skipped += 1;
                // Preserve the row with a null converted value.
                let mut attrs = feature.attributes.clone();
                attrs.push(FieldValue::Null);
                out.push(wbvector::Feature {
                    fid: 0,
                    geometry: feature.geometry.clone(),
                    attributes: attrs,
                });
                continue;
            };

            let text = format_coord(lat, lon, prm.output_notation, prm.precision);
            let mut attrs = feature.attributes.clone();
            attrs.push(FieldValue::Text(text));
            let geometry = if prm.update_geometry {
                Some(Geometry::point(lon, lat))
            } else {
                feature.geometry.clone()
            };
            out.push(wbvector::Feature {
                fid: 0,
                geometry,
                attributes: attrs,
            });
            converted += 1;
        }

        ctx.progress.info(&format!(
            "{converted} converted, {skipped} skipped ({} -> {})",
            prm.input_notation.as_str(),
            prm.output_notation.as_str()
        ));

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("converted".to_string(), json!(converted));
        outputs.insert("skipped".to_string(), json!(skipped));
        outputs.insert("output_field".to_string(), json!(prm.output_field));
        Ok(ToolRunResult { outputs })
    }
}

// ── Notation ─────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Notation {
    Dd,
    Dms,
    Ddm,
    Utm,
    Mgrs,
}

impl Notation {
    fn parse(s: &str) -> Option<Notation> {
        match s.trim().to_ascii_uppercase().as_str() {
            "DD" | "DECIMAL_DEGREES" => Some(Notation::Dd),
            "DMS" => Some(Notation::Dms),
            "DDM" => Some(Notation::Ddm),
            "UTM" => Some(Notation::Utm),
            "MGRS" | "USNG" => Some(Notation::Mgrs),
            _ => None,
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Notation::Dd => "DD",
            Notation::Dms => "DMS",
            Notation::Ddm => "DDM",
            Notation::Utm => "UTM",
            Notation::Mgrs => "MGRS",
        }
    }
}

// ── WGS84 transverse-Mercator (Krüger series, order n^3) ─────────────────────

const A_AXIS: f64 = 6_378_137.0;
const INV_F: f64 = 298.257_223_563;
const K0: f64 = 0.9996;
const FALSE_EASTING: f64 = 500_000.0;
const FALSE_NORTHING_S: f64 = 10_000_000.0;

struct TmConstants {
    n: f64,
    big_a: f64,
    alpha: [f64; 3],
    beta: [f64; 3],
    delta: [f64; 3],
}

fn tm_constants() -> TmConstants {
    let f = 1.0 / INV_F;
    let n = f / (2.0 - f);
    let n2 = n * n;
    let n3 = n2 * n;
    let big_a = A_AXIS / (1.0 + n) * (1.0 + n2 / 4.0 + n2 * n2 / 64.0);
    let alpha = [
        n / 2.0 - 2.0 / 3.0 * n2 + 5.0 / 16.0 * n3,
        13.0 / 48.0 * n2 - 3.0 / 5.0 * n3,
        61.0 / 240.0 * n3,
    ];
    let beta = [
        n / 2.0 - 2.0 / 3.0 * n2 + 37.0 / 96.0 * n3,
        1.0 / 48.0 * n2 + 1.0 / 15.0 * n3,
        17.0 / 480.0 * n3,
    ];
    let delta = [
        2.0 * n - 2.0 / 3.0 * n2 - 2.0 * n3,
        7.0 / 3.0 * n2 - 8.0 / 5.0 * n3,
        56.0 / 15.0 * n3,
    ];
    TmConstants {
        n,
        big_a,
        alpha,
        beta,
        delta,
    }
}

fn utm_zone(lon: f64) -> i32 {
    (((lon + 180.0) / 6.0).floor() as i32 + 1).clamp(1, 60)
}

fn central_meridian(zone: i32) -> f64 {
    (zone * 6 - 183) as f64
}

/// Forward: (lat, lon) degrees -> (zone, is_north, easting, northing).
fn geodetic_to_utm(lat: f64, lon: f64) -> (i32, bool, f64, f64) {
    let c = tm_constants();
    let zone = utm_zone(lon);
    let lon0 = central_meridian(zone);
    let phi = lat.to_radians();
    let dlam = (lon - lon0).to_radians();

    let sphi = phi.sin();
    let two_sqrt_n = 2.0 * c.n.sqrt() / (1.0 + c.n);
    let t = (sphi.atanh() - two_sqrt_n * (two_sqrt_n * sphi).atanh()).sinh();
    let cosl = dlam.cos();
    let xi_p = t.atan2(cosl);
    let eta_p = (dlam.sin() / (t * t + cosl * cosl).sqrt()).asinh();

    let mut xi = xi_p;
    let mut eta = eta_p;
    for j in 1..=3 {
        let jf = j as f64;
        xi += c.alpha[j - 1] * (2.0 * jf * xi_p).sin() * (2.0 * jf * eta_p).cosh();
        eta += c.alpha[j - 1] * (2.0 * jf * xi_p).cos() * (2.0 * jf * eta_p).sinh();
    }
    let easting = FALSE_EASTING + K0 * c.big_a * eta;
    let mut northing = K0 * c.big_a * xi;
    let is_north = lat >= 0.0;
    if !is_north {
        northing += FALSE_NORTHING_S;
    }
    (zone, is_north, easting, northing)
}

/// Inverse: (zone, is_north, easting, northing) -> (lat, lon) degrees.
fn utm_to_geodetic(zone: i32, is_north: bool, easting: f64, northing: f64) -> (f64, f64) {
    let c = tm_constants();
    let lon0 = central_meridian(zone);
    let n_true = if is_north {
        northing
    } else {
        northing - FALSE_NORTHING_S
    };
    let xi = n_true / (K0 * c.big_a);
    let eta = (easting - FALSE_EASTING) / (K0 * c.big_a);

    let mut xi_p = xi;
    let mut eta_p = eta;
    for j in 1..=3 {
        let jf = j as f64;
        xi_p -= c.beta[j - 1] * (2.0 * jf * xi).sin() * (2.0 * jf * eta).cosh();
        eta_p -= c.beta[j - 1] * (2.0 * jf * xi).cos() * (2.0 * jf * eta).sinh();
    }
    let chi = (xi_p.sin() / eta_p.cosh()).asin();
    let mut phi = chi;
    for j in 1..=3 {
        let jf = j as f64;
        phi += c.delta[j - 1] * (2.0 * jf * chi).sin();
    }
    let lam = lon0.to_radians() + eta_p.sinh().atan2(xi_p.cos());
    (phi.to_degrees(), lam.to_degrees())
}

// ── MGRS lettering ───────────────────────────────────────────────────────────

const COL_GROUPS: [&[u8]; 3] = [b"ABCDEFGH", b"JKLMNPQR", b"STUVWXYZ"];
const ROW_LETTERS: &[u8] = b"ABCDEFGHJKLMNPQRSTUV"; // 20, I/O omitted
const BAND8: &[u8] = b"CDEFGHJKLMNPQRSTUVW"; // 19 bands of 8° from -80..72

fn band_letter(lat: f64) -> u8 {
    if lat >= 72.0 {
        return b'X'; // X spans 72..84
    }
    let idx = (((lat + 80.0) / 8.0).floor()).clamp(0.0, 18.0) as usize;
    BAND8[idx]
}

fn band_is_north(band: u8) -> bool {
    matches!(
        band,
        b'N' | b'P' | b'Q' | b'R' | b'S' | b'T' | b'U' | b'V' | b'W' | b'X'
    )
}

/// Minimum true UTM northing (metres) for each latitude band — resolves the
/// 2,000,000m ambiguity when going MGRS -> UTM. Standard MGRS table.
fn min_northing(band: u8) -> f64 {
    match band {
        b'C' => 1_100_000.0,
        b'D' => 2_000_000.0,
        b'E' => 2_800_000.0,
        b'F' => 3_700_000.0,
        b'G' => 4_600_000.0,
        b'H' => 5_500_000.0,
        b'J' => 6_400_000.0,
        b'K' => 7_300_000.0,
        b'L' => 8_200_000.0,
        b'M' => 9_100_000.0,
        b'N' => 0.0,
        b'P' => 800_000.0,
        b'Q' => 1_600_000.0,
        b'R' => 2_400_000.0,
        b'S' => 3_300_000.0,
        b'T' => 4_100_000.0,
        b'U' => 4_900_000.0,
        b'V' => 5_700_000.0,
        b'W' => 6_600_000.0,
        b'X' => 7_500_000.0,
        _ => 0.0,
    }
}

/// Forward: (lat, lon) -> MGRS string with `precision` digits per axis.
fn geodetic_to_mgrs(lat: f64, lon: f64, precision: u32) -> String {
    let (zone, _is_north, easting, northing) = geodetic_to_utm(lat, lon);
    let band = band_letter(lat);

    let col_group = COL_GROUPS[((zone - 1).rem_euclid(3)) as usize];
    let col_num = (easting / 100_000.0).floor() as i64; // 1..8 in practice
    let col_letter = col_group[(col_num - 1).clamp(0, 7) as usize];

    let row_base = (northing / 100_000.0).floor() as i64;
    let row_off = if zone % 2 == 0 { 5 } else { 0 };
    let row_letter = ROW_LETTERS[((row_base + row_off).rem_euclid(20)) as usize];

    let e_num = easting - col_num as f64 * 100_000.0; // 0..100000
    let n_num = northing - row_base as f64 * 100_000.0;
    let scale = 10f64.powi(precision as i32 - 5);
    let width = precision as usize;
    let e_digits = ((e_num * scale).round() as i64).clamp(0, 10i64.pow(precision) - 1);
    let n_digits = ((n_num * scale).round() as i64).clamp(0, 10i64.pow(precision) - 1);

    format!(
        "{}{}{}{}{:0width$}{:0width$}",
        zone,
        band as char,
        col_letter as char,
        row_letter as char,
        e_digits,
        n_digits,
        width = width
    )
}

/// Inverse: MGRS string -> (lat, lon). None on malformed input.
fn mgrs_to_geodetic(s: &str) -> Option<(f64, f64)> {
    let up: String = s.chars().filter(|c| !c.is_whitespace()).collect();
    let up = up.to_ascii_uppercase();
    let bytes = up.as_bytes();
    let mut i = 0;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    if i == 0 || i > 2 || i + 3 > bytes.len() {
        return None;
    }
    let zone: i32 = up[..i].parse().ok()?;
    if !(1..=60).contains(&zone) {
        return None;
    }
    let band = bytes[i];
    let col_letter = bytes[i + 1];
    let row_letter = bytes[i + 2];
    let digits = &up[i + 3..];
    if digits.is_empty()
        || !digits.len().is_multiple_of(2)
        || !digits.bytes().all(|b| b.is_ascii_digit())
    {
        return None;
    }
    let precision = digits.len() / 2;
    let e_str = &digits[..precision];
    let n_str = &digits[precision..];
    let scale = 10f64.powi(5 - precision as i32);
    let e_num = e_str.parse::<f64>().ok()? * scale;
    let n_num = n_str.parse::<f64>().ok()? * scale;

    // Column letter -> easting 100k.
    let col_group = COL_GROUPS[((zone - 1).rem_euclid(3)) as usize];
    let col_pos = col_group.iter().position(|&b| b == col_letter)? as f64;
    let easting = (col_pos + 1.0) * 100_000.0 + e_num;

    // Row letter -> northing 100k (undo even-zone offset), then resolve band.
    let row_pos = ROW_LETTERS.iter().position(|&b| b == row_letter)? as i64;
    let row_off = if zone % 2 == 0 { 5 } else { 0 };
    let row_idx = (row_pos - row_off).rem_euclid(20);
    let mut north100k = row_idx as f64 * 100_000.0;
    let min_n = min_northing(band);
    while north100k < min_n {
        north100k += 2_000_000.0;
    }
    let northing = north100k + n_num;

    let (lat, lon) = utm_to_geodetic(zone, band_is_north(band), easting, northing);
    Some((lat, lon))
}

// ── DD / DMS / DDM formatting & parsing ──────────────────────────────────────

fn format_dd(lat: f64, lon: f64) -> String {
    format!("{lat:.8}, {lon:.8}")
}

/// Degrees / (minutes / seconds) breakdown of an angle magnitude.
fn dms_parts(v: f64) -> (f64, f64, f64) {
    let a = v.abs();
    let deg = a.floor();
    let min_f = (a - deg) * 60.0;
    let min = min_f.floor();
    let sec = (min_f - min) * 60.0;
    (deg, min, sec)
}

fn format_dms_axis(v: f64, is_lat: bool) -> String {
    let hemi = hemisphere(v, is_lat);
    let (mut deg, mut min, mut sec) = dms_parts(v);
    // Guard against rounding 59.9995 -> 60.000.
    if (sec * 1000.0).round() / 1000.0 >= 60.0 {
        sec = 0.0;
        min += 1.0;
    }
    if min >= 60.0 {
        min = 0.0;
        deg += 1.0;
    }
    let deg_width = if is_lat { 2 } else { 3 };
    format!(
        "{:0dw$.0}°{:02.0}'{:06.3}\"{}",
        deg,
        min,
        sec,
        hemi,
        dw = deg_width
    )
}

fn format_ddm_axis(v: f64, is_lat: bool) -> String {
    let hemi = hemisphere(v, is_lat);
    let a = v.abs();
    let mut deg = a.floor();
    let mut min = (a - deg) * 60.0;
    if (min * 1000.0).round() / 1000.0 >= 60.0 {
        min = 0.0;
        deg += 1.0;
    }
    let deg_width = if is_lat { 2 } else { 3 };
    format!("{:0dw$.0}°{:06.3}'{}", deg, min, hemi, dw = deg_width)
}

fn hemisphere(v: f64, is_lat: bool) -> char {
    match (is_lat, v >= 0.0) {
        (true, true) => 'N',
        (true, false) => 'S',
        (false, true) => 'E',
        (false, false) => 'W',
    }
}

fn format_coord(lat: f64, lon: f64, notation: Notation, precision: u32) -> String {
    match notation {
        Notation::Dd => format_dd(lat, lon),
        Notation::Dms => format!(
            "{} {}",
            format_dms_axis(lat, true),
            format_dms_axis(lon, false)
        ),
        Notation::Ddm => format!(
            "{} {}",
            format_ddm_axis(lat, true),
            format_ddm_axis(lon, false)
        ),
        Notation::Utm => {
            let (zone, is_north, e, n) = geodetic_to_utm(lat, lon);
            format!(
                "{}{} {:.3} {:.3}",
                zone,
                if is_north { 'N' } else { 'S' },
                e,
                n
            )
        }
        Notation::Mgrs => geodetic_to_mgrs(lat, lon, precision),
    }
}

/// Parses a coordinate string in `notation` into (lat, lon) degrees.
fn parse_coord(s: &str, notation: Notation) -> Option<(f64, f64)> {
    match notation {
        Notation::Dd => parse_dd(s),
        Notation::Dms | Notation::Ddm => parse_dms_ddm(s),
        Notation::Utm => parse_utm(s),
        Notation::Mgrs => mgrs_to_geodetic(s),
    }
}

/// Extracts every signed float in a string.
fn signed_numbers(s: &str) -> Vec<f64> {
    let mut out = Vec::new();
    let mut cur = String::new();
    for c in s.chars() {
        if c.is_ascii_digit() || c == '.' {
            cur.push(c);
        } else if c == '-' || c == '+' {
            if !cur.is_empty() {
                if let Ok(v) = cur.parse::<f64>() {
                    out.push(v);
                }
                cur.clear();
            }
            cur.push(c);
        } else {
            if !cur.is_empty() {
                if let Ok(v) = cur.parse::<f64>() {
                    out.push(v);
                }
                cur.clear();
            }
        }
    }
    if !cur.is_empty() {
        if let Ok(v) = cur.parse::<f64>() {
            out.push(v);
        }
    }
    out
}

fn parse_dd(s: &str) -> Option<(f64, f64)> {
    // Hemisphere letters override sign if present.
    let up = s.to_ascii_uppercase();
    if up.contains(['N', 'S', 'E', 'W']) {
        return parse_dms_ddm(s);
    }
    let nums = signed_numbers(s);
    if nums.len() < 2 {
        return None;
    }
    Some((nums[0], nums[1])) // lat, lon
}

/// Parses DMS/DDM (and hemisphere-tagged DD) strings: splits into two axes at the
/// first N/S/E/W and reads deg[/min[/sec]] tokens for each.
fn parse_dms_ddm(s: &str) -> Option<(f64, f64)> {
    let up = s.to_ascii_uppercase();
    let idx = up.find(['N', 'S', 'E', 'W'])?;
    let (a, b) = up.split_at(idx + 1);
    let (av, ah) = axis_value(a)?;
    let (bv, bh) = axis_value(b)?;
    let mut lat = None;
    let mut lon = None;
    for (v, h) in [(av, ah), (bv, bh)] {
        match h {
            'N' => lat = Some(v),
            'S' => lat = Some(-v),
            'E' => lon = Some(v),
            'W' => lon = Some(-v),
            _ => {}
        }
    }
    Some((lat?, lon?))
}

/// Reads one axis chunk -> (magnitude in degrees, hemisphere char).
fn axis_value(chunk: &str) -> Option<(f64, char)> {
    let hemi = chunk.chars().rev().find(|c| "NSEW".contains(*c))?;
    let nums = signed_numbers(chunk);
    if nums.is_empty() {
        return None;
    }
    let deg = nums.first().copied().unwrap_or(0.0).abs();
    let min = nums.get(1).copied().unwrap_or(0.0);
    let sec = nums.get(2).copied().unwrap_or(0.0);
    Some((deg + min / 60.0 + sec / 3600.0, hemi))
}

fn parse_utm(s: &str) -> Option<(f64, f64)> {
    let up = s.trim().to_ascii_uppercase();
    let bytes = up.as_bytes();
    let mut i = 0;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    if i == 0 {
        return None;
    }
    let zone: i32 = up[..i].parse().ok()?;
    if !(1..=60).contains(&zone) {
        return None;
    }
    let hemi = up[i..].chars().find(|c| *c == 'N' || *c == 'S')?;
    let nums = signed_numbers(&up[i..]);
    if nums.len() < 2 {
        return None;
    }
    let easting = nums[nums.len() - 2];
    let northing = nums[nums.len() - 1];
    let (lat, lon) = utm_to_geodetic(zone, hemi == 'N', easting, northing);
    Some((lat, lon))
}

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        _ => None,
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    input_notation: Notation,
    output_notation: Notation,
    coord_field: Option<String>,
    output_field: String,
    precision: u32,
    update_geometry: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let input_notation = require_notation(args, "input_notation")?;
    let output_notation = require_notation(args, "output_notation")?;
    let coord_field = parse_optional_str(args, "coord_field")?.map(str::to_string);
    if coord_field.is_none() && input_notation != Notation::Dd {
        return Err(ToolError::Validation(
            "coord_field is required unless input_notation=DD (point geometry as lon/lat)".into(),
        ));
    }
    let output_field = parse_optional_str(args, "output_field")?
        .map(str::to_string)
        .unwrap_or_else(|| output_notation.as_str().to_ascii_lowercase());
    let precision = opt_u32(args, "precision")?.unwrap_or(5).clamp(1, 9);
    let update_geometry = opt_bool(args, "update_geometry")?.unwrap_or(false);
    Ok(Params {
        input_notation,
        output_notation,
        coord_field,
        output_field,
        precision,
        update_geometry,
    })
}

fn require_notation(args: &ToolArgs, key: &str) -> Result<Notation, ToolError> {
    let s = require_str(args, key)?;
    Notation::parse(s).ok_or_else(|| {
        ToolError::Validation(format!(
            "parameter '{key}' must be one of DD, DMS, DDM, UTM, MGRS (got '{s}')"
        ))
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn opt_u32(args: &ToolArgs, key: &str) -> Result<Option<u32>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => n
            .as_u64()
            .map(|v| Some(v as u32))
            .ok_or_else(|| ToolError::Validation(format!("parameter '{key}' must be an integer"))),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<u32>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be an integer"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be an integer"
        ))),
    }
}

fn opt_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
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
    use wbvector::memory_store;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Point layer (lon/lat, EPSG:4326) with a "name" field.
    fn point_layer(rows: &[(&str, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("name", FieldType::Text));
        for (name, lon, lat) in rows {
            l.add_feature(
                Some(Geometry::point(*lon, *lat)),
                &[("name", (*name).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ConvertCoordinateNotationTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn field_str(layer: &Layer, feat: usize, name: &str) -> String {
        let idx = layer.schema.field_index(name).unwrap();
        layer.features[feat].attributes[idx]
            .as_str()
            .unwrap()
            .to_string()
    }

    /// White House ~ (38.8977, -77.0365) -> UTM zone 18N, MGRS square 18SUJ.
    #[test]
    fn white_house_utm_and_mgrs() {
        let (zone, is_north, e, n) = geodetic_to_utm(38.8977, -77.0365);
        assert_eq!(zone, 18);
        assert!(is_north);
        // pyproj EPSG:4326 -> 32618 gives 323394.296 E, 4307395.634 N.
        assert!((e - 323394.296).abs() < 0.01, "easting {e}");
        assert!((n - 4307395.634).abs() < 0.01, "northing {n}");
        let mgrs = geodetic_to_mgrs(38.8977, -77.0365, 5);
        assert!(mgrs.starts_with("18SUJ"), "mgrs {mgrs}");
    }

    /// DD -> UTM -> DD round-trips to well under 1e-6 deg.
    #[test]
    fn roundtrip_dd_utm() {
        for &(lon, lat) in &[
            (-77.0365, 38.8977),
            (2.2945, 48.8584),
            (139.6917, 35.6895),
            (151.2093, -33.8688),
        ] {
            let utm = format_coord(lat, lon, Notation::Utm, 5);
            let (lat2, lon2) = parse_coord(&utm, Notation::Utm).unwrap();
            assert!((lat2 - lat).abs() < 1e-6, "lat {lat} -> {lat2} via {utm}");
            assert!((lon2 - lon).abs() < 1e-6, "lon {lon} -> {lon2} via {utm}");
        }
    }

    /// DD -> MGRS -> DD round-trips within 1e-6 deg at 7-digit precision.
    #[test]
    fn roundtrip_dd_mgrs() {
        for &(lon, lat) in &[(-77.0365, 38.8977), (2.2945, 48.8584), (151.2093, -33.8688)] {
            let mgrs = geodetic_to_mgrs(lat, lon, 7);
            let (lat2, lon2) = mgrs_to_geodetic(&mgrs).unwrap();
            assert!((lat2 - lat).abs() < 1e-6, "lat {lat} -> {lat2} via {mgrs}");
            assert!((lon2 - lon).abs() < 1e-6, "lon {lon} -> {lon2} via {mgrs}");
        }
    }

    /// DMS / DDM formatting round-trips through the parser.
    #[test]
    fn roundtrip_dms_ddm() {
        let (lat, lon) = (38.8977, -77.0365);
        for note in [Notation::Dms, Notation::Ddm] {
            let s = format_coord(lat, lon, note, 5);
            let (lat2, lon2) = parse_coord(&s, note).unwrap();
            assert!((lat2 - lat).abs() < 1e-6, "lat via {s}");
            assert!((lon2 - lon).abs() < 1e-6, "lon via {s}");
        }
    }

    /// End-to-end: DD geometry in, MGRS field out; converted count correct.
    #[test]
    fn adds_mgrs_field_from_geometry() {
        let input = point_layer(&[("wh", -77.0365, 38.8977), ("eiffel", 2.2945, 48.8584)]);
        let (out, layer) = run(json!({
            "input": input,
            "input_notation": "DD",
            "output_notation": "MGRS",
        }));
        assert_eq!(out.outputs["converted"], json!(2));
        assert!(field_str(&layer, 0, "mgrs").starts_with("18SUJ"));
    }

    /// A feature with no geometry passes through with a null converted value.
    #[test]
    fn passes_through_null_geometry() {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.push(wbvector::Feature {
            fid: 0,
            geometry: None,
            attributes: vec![FieldValue::Text("empty".into())],
        });
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);
        let (out, _l) = run(json!({
            "input": input, "input_notation": "DD", "output_notation": "UTM",
        }));
        assert_eq!(out.outputs["converted"], json!(0));
        assert_eq!(out.outputs["skipped"], json!(1));
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ConvertCoordinateNotationTool.validate(&args)
        };
        // Missing notations.
        assert!(bad(json!({ "input": "a.geojson" })).is_err());
        // Unknown notation.
        assert!(bad(json!({
            "input": "a.geojson", "input_notation": "FOO", "output_notation": "DD",
        }))
        .is_err());
        // UTM input needs a coord_field (no geometry parsing for UTM).
        assert!(bad(json!({
            "input": "a.geojson", "input_notation": "UTM", "output_notation": "DD",
        }))
        .is_err());
        // Valid.
        assert!(bad(json!({
            "input": "a.geojson", "input_notation": "DD", "output_notation": "MGRS",
        }))
        .is_ok());
    }
}
