//! GeoLibre tool: scan free-form text for embedded coordinates and emit points.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Extract Locations From Text* (the
//! coordinate/pattern branch — the document/address-locator sibling needs a
//! geocoding service and is out of scope). Where
//! [`convert_coordinate_notation`](crate::convert_coordinate_notation)
//! *reformats* a coordinate already attached to each feature, this tool reads a
//! blob of arbitrary text and finds coordinates *embedded* in it, emitting one
//! WGS84 point feature per hit with the matched substring and its notation.
//!
//! Notations recognised (all detected without any geocoding service):
//! - **DD**  decimal degrees        — `"38.8977, -77.0365"`, `"38.8977°N 77.0365°W"`
//! - **DDM** degrees-decimal-minutes — `"38°53.862'N 077°02.190'W"`
//! - **DMS** degrees-minutes-seconds — `"38°53'51.72\"N 077°02'11.40\"W"`
//! - **MGRS**/USNG grid reference    — `"18SUJ2348006479"`
//!
//! The scanners are hand-written state machines (no `regex` dependency, so no
//! new crate and nothing that pulls C code): a per-axis matcher recognises
//! degree/minute/second tokens ending in an N/S/E/W hemisphere letter and pairs
//! an adjacent latitude token with a longitude token; a separate matcher reads
//! comma-delimited signed-decimal pairs; and an MGRS matcher validates each
//! candidate by running it back through the existing `mgrs_to_geodetic` grid
//! math. Everything is deterministic and WASM-safe.
//!
//! Scope for v1: latitude/longitude embedded coordinates in the four notations
//! above (GARS/GEOREF and UTM-with-hemisphere prose are not scanned); the
//! decimal-degree pair form requires a comma between the two numbers to keep
//! false positives out of ordinary prose.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::convert_coordinate_notation::mgrs_to_geodetic;
use crate::vector_common::{parse_optional_str, write_or_store_layer};

pub struct ExtractLocationsFromTextTool;

impl Tool for ExtractLocationsFromTextTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "extract_locations_from_text",
            display_name: "Extract Locations From Text",
            summary: "Scan free-form text for embedded coordinates — decimal degrees (DD), degrees-decimal-minutes (DDM), degrees-minutes-seconds (DMS), and MGRS/USNG grid references — and emit one WGS84 point per hit, like ArcGIS Extract Locations From Text (coordinate branch). Pure pattern scanning; no geocoding service.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input_text",
                    description: "The free-form text to scan. Either the literal text, or a path to a UTF-8 text file (read when the value resolves to an existing file).",
                    required: true,
                },
                ToolParamSpec {
                    name: "notations",
                    description: "Comma/space separated subset of notations to scan for: DD, DDM, DMS, MGRS (USNG is an alias of MGRS). Default: all four.",
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
        require_str(args, "input_text")?;
        parse_notations(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let raw = require_str(args, "input_text")?;
        let output = parse_optional_str(args, "output")?;
        let enabled = parse_notations(args)?;

        // Accept either literal text or a path to a UTF-8 text file.
        let text = read_text(raw)?;

        let hits = extract_locations(&text, &enabled);
        ctx.progress
            .info(&format!("found {} embedded coordinate(s)", hits.len()));

        let mut out = Layer::new("locations")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        out.add_field(FieldDef::new("match", FieldType::Text));
        out.add_field(FieldDef::new("notation", FieldType::Text));
        out.add_field(FieldDef::new("latitude", FieldType::Float));
        out.add_field(FieldDef::new("longitude", FieldType::Float));
        out.add_field(FieldDef::new("char_offset", FieldType::Integer));

        // Per-notation tallies for the run summary.
        let mut counts: BTreeMap<&'static str, usize> = BTreeMap::new();
        for hit in &hits {
            *counts.entry(hit.notation.as_str()).or_insert(0) += 1;
            out.add_feature(
                Some(Geometry::point(hit.lon, hit.lat)),
                &[
                    ("match", FieldValue::Text(hit.text.clone())),
                    (
                        "notation",
                        FieldValue::Text(hit.notation.as_str().to_string()),
                    ),
                    ("latitude", FieldValue::Float(hit.lat)),
                    ("longitude", FieldValue::Float(hit.lon)),
                    ("char_offset", FieldValue::Integer(hit.offset as i64)),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed building feature: {e}")))?;
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("count".to_string(), json!(hits.len()));
        outputs.insert(
            "counts_by_notation".to_string(),
            json!(counts
                .into_iter()
                .map(|(k, v)| (k.to_string(), v))
                .collect::<BTreeMap<_, _>>()),
        );
        Ok(ToolRunResult { outputs })
    }
}

// ── Notation ─────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Notation {
    // Ordering matters: a paired axis is classified by the "richer" of its two
    // components (Dms > Ddm > Dd), so keep the derived Ord in that order.
    Dd,
    Ddm,
    Dms,
    Mgrs,
}

impl Notation {
    fn parse(s: &str) -> Option<Notation> {
        match s.trim().to_ascii_uppercase().as_str() {
            "DD" | "DECIMAL_DEGREES" => Some(Notation::Dd),
            "DDM" => Some(Notation::Ddm),
            "DMS" => Some(Notation::Dms),
            "MGRS" | "USNG" => Some(Notation::Mgrs),
            _ => None,
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Notation::Dd => "DD",
            Notation::Ddm => "DDM",
            Notation::Dms => "DMS",
            Notation::Mgrs => "MGRS",
        }
    }
}

fn parse_notations(args: &ToolArgs) -> Result<Vec<Notation>, ToolError> {
    let s = match parse_optional_str(args, "notations")? {
        None => {
            return Ok(vec![
                Notation::Dd,
                Notation::Ddm,
                Notation::Dms,
                Notation::Mgrs,
            ])
        }
        Some(s) => s,
    };
    let mut out = Vec::new();
    for tok in s
        .split([',', ' ', ';', '|'])
        .filter(|t| !t.trim().is_empty())
    {
        let n = Notation::parse(tok).ok_or_else(|| {
            ToolError::Validation(format!(
                "unknown notation '{tok}' (expected DD, DDM, DMS, or MGRS)"
            ))
        })?;
        if !out.contains(&n) {
            out.push(n);
        }
    }
    if out.is_empty() {
        return Err(ToolError::Validation(
            "parameter 'notations' selected no valid notations".into(),
        ));
    }
    Ok(out)
}

// ── Scanning ─────────────────────────────────────────────────────────────────

struct Hit {
    offset: usize, // char index of the match start in the text
    text: String,  // the matched substring
    lat: f64,
    lon: f64,
    notation: Notation,
}

const MAX_PAIR_GAP: usize = 25; // chars allowed between a lat token and its lon token

/// Scans `text` for embedded coordinates in the `enabled` notations and returns
/// them ordered by position, with overlapping matches dropped (earliest wins).
fn extract_locations(text: &str, enabled: &[Notation]) -> Vec<Hit> {
    let chars: Vec<char> = text.chars().collect();
    let want = |n: Notation| enabled.contains(&n);
    let want_axis = want(Notation::Dd) || want(Notation::Ddm) || want(Notation::Dms);

    let mut hits: Vec<Hit> = Vec::new();

    // MGRS / USNG grid references.
    if want(Notation::Mgrs) {
        let mut i = 0;
        while i < chars.len() {
            if let Some((end, lat, lon, s)) = match_mgrs(&chars, i) {
                hits.push(Hit {
                    offset: i,
                    text: s,
                    lat,
                    lon,
                    notation: Notation::Mgrs,
                });
                i = end;
            } else {
                i += 1;
            }
        }
    }

    // Hemisphere-tagged DD/DDM/DMS: collect per-axis tokens, then pair an
    // adjacent latitude (N/S) with a longitude (E/W).
    if want_axis {
        let mut tokens: Vec<AxisToken> = Vec::new();
        let mut i = 0;
        while i < chars.len() {
            if let Some(tok) = match_axis(&chars, i) {
                let next = tok.end;
                tokens.push(tok);
                i = next;
            } else {
                i += 1;
            }
        }
        let mut a = 0;
        while a + 1 < tokens.len() {
            let (ta, tb) = (&tokens[a], &tokens[a + 1]);
            let a_lat = matches!(ta.hemi, 'N' | 'S');
            let b_lat = matches!(tb.hemi, 'N' | 'S');
            if a_lat == b_lat || tb.start.saturating_sub(ta.end) > MAX_PAIR_GAP {
                a += 1;
                continue;
            }
            let (lat_tok, lon_tok) = if a_lat { (ta, tb) } else { (tb, ta) };
            let lat = signed(lat_tok.value, lat_tok.hemi);
            let lon = signed(lon_tok.value, lon_tok.hemi);
            let kind = ta.kind.max(tb.kind);
            if want(kind) && lat.abs() <= 90.0 && lon.abs() <= 180.0 {
                let text: String = chars[ta.start..tb.end].iter().collect();
                hits.push(Hit {
                    offset: ta.start,
                    text,
                    lat,
                    lon,
                    notation: kind,
                });
            }
            a += 2; // both tokens consumed by this pair
        }
    }

    // Comma-delimited signed-decimal pairs ("lat, lon").
    if want(Notation::Dd) {
        let mut i = 0;
        while i < chars.len() {
            if let Some((end, lat, lon)) = match_decimal_pair(&chars, i) {
                let text: String = chars[i..end].iter().collect();
                hits.push(Hit {
                    offset: i,
                    text,
                    lat,
                    lon,
                    notation: Notation::Dd,
                });
                i = end;
            } else {
                i += 1;
            }
        }
    }

    // Order by position and drop overlaps (an earlier, longer hemisphere match
    // wins over a bare-decimal match that starts inside it, etc.).
    hits.sort_by_key(|h| (h.offset, usize::MAX - h.text.chars().count()));
    let mut kept: Vec<Hit> = Vec::new();
    let mut covered_to = 0usize; // exclusive char index already consumed
    for hit in hits {
        let start = hit.offset;
        let end = start + hit.text.chars().count();
        if start < covered_to {
            continue;
        }
        covered_to = end;
        kept.push(hit);
    }
    kept
}

/// Applies a hemisphere letter's sign to a positive magnitude.
fn signed(mag: f64, hemi: char) -> f64 {
    match hemi {
        'S' | 'W' => -mag,
        _ => mag,
    }
}

struct AxisToken {
    start: usize,
    end: usize, // exclusive
    value: f64, // positive magnitude in degrees
    hemi: char, // N/S/E/W
    kind: Notation,
}

/// Matches one hemisphere-tagged axis starting at `i`:
/// `DEG [° [MIN ' [SEC "]]] N|S|E|W`, with optional minute/second parts.
fn match_axis(chars: &[char], i: usize) -> Option<AxisToken> {
    // Must start on a boundary so we don't grab digits out of the middle of a
    // longer number/word.
    if i > 0 && (chars[i - 1].is_ascii_digit() || chars[i - 1] == '.') {
        return None;
    }
    let (deg, mut j) = parse_number(chars, i)?;

    // Optional degree symbol (with optional leading whitespace).
    let k = skip_ws(chars, j);
    if k < chars.len() && (chars[k] == '\u{00B0}' || chars[k] == '\u{00BA}') {
        j = k + 1;
    }

    let mut minutes = 0.0;
    let mut min_present = false;
    let mut seconds = 0.0;
    let mut sec_present = false;

    // Minutes: a number immediately followed by a minute marker.
    let m = skip_ws(chars, j);
    if let Some((num, m2)) = parse_number(chars, m) {
        let m3 = skip_ws(chars, m2);
        if m3 < chars.len() && is_minute_marker(chars[m3]) {
            minutes = num;
            min_present = true;
            j = m3 + 1;
            // Seconds: only when minutes were present.
            let s = skip_ws(chars, j);
            if let Some((snum, s2)) = parse_number(chars, s) {
                let s3 = skip_ws(chars, s2);
                if s3 < chars.len() && is_second_marker(chars[s3]) {
                    seconds = snum;
                    sec_present = true;
                    j = s3 + 1;
                } else if s3 + 1 < chars.len()
                    && is_minute_marker(chars[s3])
                    && is_minute_marker(chars[s3 + 1])
                {
                    // Two apostrophes used as a seconds marker.
                    seconds = snum;
                    sec_present = true;
                    j = s3 + 2;
                }
            }
        }
    }

    // Hemisphere letter.
    let h = skip_ws(chars, j);
    if h >= chars.len() {
        return None;
    }
    let hc = chars[h].to_ascii_uppercase();
    if !matches!(hc, 'N' | 'S' | 'E' | 'W') {
        return None;
    }
    // The hemisphere letter must not be part of a longer word.
    if h + 1 < chars.len() && chars[h + 1].is_ascii_alphabetic() {
        return None;
    }
    // Reject nonsensical magnitudes (minutes/seconds out of range).
    if min_present && minutes >= 60.0 {
        return None;
    }
    if sec_present && seconds >= 60.0 {
        return None;
    }

    let value = deg + minutes / 60.0 + seconds / 3600.0;
    let kind = if sec_present {
        Notation::Dms
    } else if min_present {
        Notation::Ddm
    } else {
        Notation::Dd
    };
    Some(AxisToken {
        start: i,
        end: h + 1,
        value,
        hemi: hc,
        kind,
    })
}

/// Matches a comma-delimited signed-decimal pair `LAT , LON` at `i`, where both
/// numbers carry a decimal point and fall in the valid lat/lon ranges. Returns
/// `(end, lat, lon)`.
fn match_decimal_pair(chars: &[char], i: usize) -> Option<(usize, f64, f64)> {
    // Boundary: don't start mid-number/word.
    if i > 0
        && (chars[i - 1].is_ascii_digit()
            || chars[i - 1] == '.'
            || chars[i - 1].is_ascii_alphabetic())
    {
        return None;
    }
    let (lat, j) = parse_signed_decimal(chars, i)?;
    let j = skip_ws(chars, j);
    if j >= chars.len() || chars[j] != ',' {
        return None;
    }
    let j = skip_ws(chars, j + 1);
    let (lon, end) = parse_signed_decimal(chars, j)?;
    // Trailing boundary.
    if end < chars.len() && (chars[end].is_ascii_digit() || chars[end] == '.') {
        return None;
    }
    if lat.abs() > 90.0 || lon.abs() > 180.0 {
        return None;
    }
    Some((end, lat, lon))
}

/// Matches an MGRS/USNG grid reference at `i`. Returns `(end, lat, lon, text)`.
fn match_mgrs(chars: &[char], i: usize) -> Option<(usize, f64, f64, String)> {
    if i > 0 && chars[i - 1].is_ascii_alphanumeric() {
        return None;
    }
    let mut j = i;
    // Zone: 1 or 2 digits.
    let z0 = j;
    while j < chars.len() && chars[j].is_ascii_digit() && j - z0 < 2 {
        j += 1;
    }
    if j == z0 || (j < chars.len() && chars[j].is_ascii_digit()) {
        return None;
    }
    // Band letter, then two grid letters.
    if j >= chars.len() || !is_band_letter(chars[j]) {
        return None;
    }
    j += 1;
    for _ in 0..2 {
        if j >= chars.len() || !is_grid_letter(chars[j]) {
            return None;
        }
        j += 1;
    }
    // An even number of location digits (2..=10 -> precision 1..=5).
    let d0 = j;
    while j < chars.len() && chars[j].is_ascii_digit() {
        j += 1;
    }
    let ndig = j - d0;
    if !(2..=10).contains(&ndig) || !ndig.is_multiple_of(2) {
        return None;
    }
    if j < chars.len() && chars[j].is_ascii_alphanumeric() {
        return None;
    }
    let s: String = chars[i..j].iter().collect();
    let (lat, lon) = mgrs_to_geodetic(&s)?;
    Some((j, lat, lon, s))
}

// ── Small character helpers ──────────────────────────────────────────────────

fn skip_ws(chars: &[char], mut j: usize) -> usize {
    while j < chars.len() && chars[j].is_whitespace() {
        j += 1;
    }
    j
}

/// Parses an unsigned number (digits with at most one decimal point) at `j`.
fn parse_number(chars: &[char], j: usize) -> Option<(f64, usize)> {
    let start = j;
    let mut k = j;
    let mut seen_dot = false;
    let mut has_digit = false;
    while k < chars.len() {
        let c = chars[k];
        if c.is_ascii_digit() {
            has_digit = true;
            k += 1;
        } else if c == '.' && !seen_dot {
            seen_dot = true;
            k += 1;
        } else {
            break;
        }
    }
    if !has_digit {
        return None;
    }
    let s: String = chars[start..k].iter().collect();
    s.parse::<f64>().ok().map(|v| (v, k))
}

/// Parses an optionally-signed decimal number that must contain a decimal point.
fn parse_signed_decimal(chars: &[char], j: usize) -> Option<(f64, usize)> {
    let mut k = j;
    let neg = k < chars.len() && (chars[k] == '-' || chars[k] == '+');
    let minus = k < chars.len() && chars[k] == '-';
    if neg {
        k += 1;
    }
    let start = k;
    let mut seen_dot = false;
    let mut has_digit = false;
    while k < chars.len() {
        let c = chars[k];
        if c.is_ascii_digit() {
            has_digit = true;
            k += 1;
        } else if c == '.' && !seen_dot {
            seen_dot = true;
            k += 1;
        } else {
            break;
        }
    }
    if !has_digit || !seen_dot {
        return None;
    }
    let s: String = chars[start..k].iter().collect();
    let v = s.parse::<f64>().ok()?;
    Some((if minus { -v } else { v }, k))
}

fn is_minute_marker(c: char) -> bool {
    matches!(c, '\'' | '\u{2032}' | '\u{2019}')
}

fn is_second_marker(c: char) -> bool {
    matches!(c, '"' | '\u{2033}' | '\u{201D}')
}

/// MGRS latitude-band letters: C..X excluding I and O.
fn is_band_letter(c: char) -> bool {
    let c = c.to_ascii_uppercase();
    ('C'..='X').contains(&c) && c != 'I' && c != 'O'
}

/// MGRS 100km-square column/row letters: A..Z excluding I and O.
fn is_grid_letter(c: char) -> bool {
    let c = c.to_ascii_uppercase();
    c.is_ascii_uppercase() && c != 'I' && c != 'O'
}

// ── Parameters ────────────────────────────────────────────────────────────────

/// Reads `input_text` either as a literal string or, when it resolves to an
/// existing file, as that file's UTF-8 contents.
fn read_text(raw: &str) -> Result<String, ToolError> {
    let trimmed = raw.trim();
    if !trimmed.is_empty() && !trimmed.contains('\n') && std::path::Path::new(trimmed).is_file() {
        return std::fs::read_to_string(trimmed)
            .map_err(|e| ToolError::Execution(format!("failed reading input_text file: {e}")));
    }
    Ok(raw.to_string())
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
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

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ExtractLocationsFromTextTool.run(&args, &ctx()).unwrap();
        let layer = crate::vector_common::load_input_layer(out.outputs["output"].as_str().unwrap())
            .unwrap();
        (out, layer)
    }

    fn feats(layer: &Layer) -> Vec<(String, String, f64, f64)> {
        let ni = layer.schema.field_index("notation").unwrap();
        let la = layer.schema.field_index("latitude").unwrap();
        let lo = layer.schema.field_index("longitude").unwrap();
        layer
            .features
            .iter()
            .map(|f| {
                (
                    f.attributes[ni].as_str().unwrap().to_string(),
                    match &f.geometry {
                        Some(Geometry::Point(_)) => "Point".to_string(),
                        _ => "?".to_string(),
                    },
                    f.attributes[la].as_f64().unwrap(),
                    f.attributes[lo].as_f64().unwrap(),
                )
            })
            .collect()
    }

    /// The White House in all four notations is found; each lands ~on the same
    /// spot (~38.8977, -77.0365).
    #[test]
    fn finds_all_notations() {
        let text = "Reports place it at 38.8977, -77.0365 (DD), \
                    38°53.862'N 077°02.190'W (DDM), \
                    38°53'51.72\"N 077°02'11.40\"W (DMS), and grid 18SUJ2348006479.";
        let (out, layer) = run(json!({ "input_text": text }));
        assert_eq!(out.outputs["count"], json!(4));
        let mut seen: Vec<String> = feats(&layer).iter().map(|f| f.0.clone()).collect();
        seen.sort();
        assert_eq!(seen, vec!["DD", "DDM", "DMS", "MGRS"]);
        for (_note, geom, lat, lon) in feats(&layer) {
            assert_eq!(geom, "Point");
            assert!((lat - 38.8977).abs() < 0.02, "lat {lat}");
            assert!((lon + 77.0365).abs() < 0.02, "lon {lon}");
        }
    }

    /// The notations filter restricts which patterns are emitted.
    #[test]
    fn notations_filter() {
        let text = "38.8977, -77.0365 and 38°53'51.72\"N 077°02'11.40\"W and 18SUJ2348006479";
        let (out, layer) = run(json!({ "input_text": text, "notations": "MGRS" }));
        assert_eq!(out.outputs["count"], json!(1));
        assert_eq!(feats(&layer)[0].0, "MGRS");
    }

    /// Text with no coordinates yields an empty (but valid) point layer.
    #[test]
    fn no_coordinates_passes_through_empty() {
        let text = "The quick brown fox jumps over the lazy dog. Chapter 3.5, page 12.";
        let (out, layer) = run(json!({ "input_text": text }));
        assert_eq!(out.outputs["count"], json!(0));
        assert_eq!(layer.features.len(), 0);
        // Schema is still present.
        assert!(layer.schema.field_index("notation").is_some());
    }

    /// Hemisphere letters override sign, and a lon-first ordering still pairs.
    #[test]
    fn lon_first_and_hemisphere() {
        // Longitude before latitude, southern/eastern hemisphere (Sydney Opera).
        let text = "Position 151°12'33.5\"E 33°51'24.5\"S was logged.";
        let (out, layer) = run(json!({ "input_text": text, "notations": "DMS" }));
        assert_eq!(out.outputs["count"], json!(1));
        let (_n, _g, lat, lon) = feats(&layer)[0].clone();
        assert!((lat + 33.8568).abs() < 0.01, "lat {lat}");
        assert!((lon - 151.2093).abs() < 0.01, "lon {lon}");
    }

    /// Out-of-range decimal pairs (ordinary prose numbers) are not emitted.
    #[test]
    fn rejects_out_of_range_decimals() {
        let text = "Prices rose 120.5, 250.7 percent over the decade.";
        let (out, _layer) = run(json!({ "input_text": text, "notations": "DD" }));
        assert_eq!(out.outputs["count"], json!(0));
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ExtractLocationsFromTextTool.validate(&args)
        };
        // Missing input_text.
        assert!(bad(json!({ "notations": "DD" })).is_err());
        // Unknown notation.
        assert!(bad(json!({ "input_text": "x", "notations": "FOO" })).is_err());
        // Valid.
        assert!(bad(json!({ "input_text": "38.0, -77.0" })).is_ok());
    }
}
