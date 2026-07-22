//! GeoLibre tool: space-time kernel density on a timestamped point layer.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Space Time Kernel Density* (Spatial
//! Analyst / Space Time Pattern Mining). The bundled `heat_map` tool is a
//! spatial-only 2-D KDE — it flattens all points onto one surface with no notion
//! of *when* each incident happened. This tool adds the temporal axis: it
//! produces a **multiband raster**, one band per time slice, where each band is a
//! density surface estimated from a separable **spatial x temporal** kernel so a
//! slice aggregates contributions from points that are near in both space *and*
//! time.
//!
//! 1. **Bin the time axis.** Points carry a `time_field` (a numeric field in its
//!    own units, or an ISO-8601 date/datetime string parsed to seconds). Slices
//!    are placed every `time_step` across `[t_min, t_max]`; band `b` is centered
//!    at `t_min + b * time_step`. (Same cube-binning shape as
//!    `emerging_hot_spot_analysis` / `reconstruct_tracks`.)
//! 2. **Separable kernel.** Each point deposits mass into a slice weighted by a
//!    spatial kernel (Epanechnikov, default, or quartic; radius
//!    `spatial_bandwidth`) times a temporal kernel (triangular, default, or
//!    Epanechnikov; half-width `temporal_bandwidth`). The kernels are separable,
//!    so the contribution of point *i* to cell *c* in band *b* is
//!    `w_i * S_ic * T_ib`.
//! 3. **Mass-conserving normalization.** For each point the spatial footprint
//!    (the cells within `spatial_bandwidth`) is normalized to sum to 1, and the
//!    temporal weights across slices are normalized to sum to 1. Hence every
//!    point of weight `w_i` deposits exactly `w_i` of mass, split across cells and
//!    slices. The **integrated mass of band `b` equals the weighted count of
//!    points inside that slice's temporal kernel window**, and the total mass
//!    across all bands equals the total point weight — the property the tool is
//!    validated on. Raster cells store a *density* (mass divided by cell area, in
//!    m^-2 for geographic input via a metre-scaled cell area, or CRS-unit^-2 for
//!    projected input), so `sum(band) * cell_area` recovers the band mass.
//!
//! Distances honor the CRS: geographic input (EPSG:4326 or untagged lon/lat) uses
//! haversine metres, projected input uses planar CRS units. Band time ranges are
//! written into the raster metadata (`band_1`, `band_2`, ...) and returned in the
//! run outputs so a UI can label each slice.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
use wbvector::{FieldValue, Geometry};

use crate::common::{parse_optional_output, write_or_store_output};
use crate::vector_common::{load_input_layer, parse_optional_str};

/// Metres per degree of latitude (and of longitude at the equator), spherical
/// approximation used to scale cell areas and the default bandwidth for
/// geographic input.
const M_PER_DEG: f64 = 111_320.0;
/// Mean Earth radius (metres) for haversine distances.
const EARTH_R: f64 = 6_371_000.0;

pub struct SpaceTimeKernelDensityTool;

impl Tool for SpaceTimeKernelDensityTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "space_time_kernel_density",
            display_name: "Space Time Kernel Density",
            summary: "Space-time kernel density of timestamped points: one density band per time slice (a multiband raster) from a separable spatial (Epanechnikov/quartic) x temporal (triangular/Epanechnikov) kernel, so each band aggregates points near in both space and time, like ArcGIS Space Time Kernel Density. The time-aware companion to the bundled 2-D heat_map.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer. Geographic (EPSG:4326/untagged lon-lat) uses haversine metres; projected uses CRS units. Non-point geometries use their first/representative vertex.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output multiband raster path (one band per time slice; driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "time_field",
                    description: "Field holding each point's time: a numeric field (its own units) or an ISO-8601 date/datetime string.",
                    required: true,
                },
                ToolParamSpec {
                    name: "time_step",
                    description: "Slice interval: a plain number (numeric time_field) or a duration like '1w', '7d', '12h', '1M'. Default '1w'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "temporal_bandwidth",
                    description: "Temporal kernel half-width, same units/duration syntax as time_step. Default: equal to time_step.",
                    required: false,
                },
                ToolParamSpec {
                    name: "spatial_bandwidth",
                    description: "Spatial kernel radius, in metres for geographic input or CRS units for projected. Default: max extent / 20.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Output cell size in CRS units (degrees for geographic). Default: max extent / 200.",
                    required: false,
                },
                ToolParamSpec {
                    name: "weight_field",
                    description: "Optional numeric field weighting each point. Default: every point counts as 1.",
                    required: false,
                },
                ToolParamSpec {
                    name: "spatial_kernel",
                    description: "Spatial kernel shape: 'epanechnikov' (default) or 'quartic'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "temporal_kernel",
                    description: "Temporal kernel shape: 'triangular' (default) or 'epanechnikov'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "epsg",
                    description: "EPSG to tag the output raster (default: from the input layer).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "time_field")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_output(args, "output")?;
        let prm = parse_params(args)?;

        ctx.progress.info("reading input points");
        let layer = load_input_layer(input)?;
        let time_idx = layer.schema.field_index(&prm.time_field).ok_or_else(|| {
            ToolError::Validation(format!("time_field '{}' not found", prm.time_field))
        })?;
        let weight_idx =
            match &prm.weight_field {
                Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                    ToolError::Validation(format!("weight_field '{f}' not found"))
                })?),
                None => None,
            };

        // Is the input geographic (haversine) or projected (planar)?
        let geographic = matches!(layer.crs_epsg(), None | Some(4326));

        // ── Collect (x, y, time, weight); track whether time came from strings ──
        let mut pts: Vec<Obs> = Vec::new();
        let mut time_is_iso = false;
        let mut skipped = 0u64;
        for feat in layer.iter() {
            let Some((x, y)) = feat.geometry.as_ref().and_then(point_xy) else {
                skipped += 1;
                continue;
            };
            let Some(fv) = feat.attributes.get(time_idx) else {
                skipped += 1;
                continue;
            };
            let time = match parse_time_value(fv) {
                Some((t, iso)) => {
                    time_is_iso |= iso;
                    t
                }
                None => {
                    skipped += 1;
                    continue;
                }
            };
            let weight = match weight_idx {
                Some(i) => match feat.attributes.get(i).and_then(FieldValue::as_f64) {
                    Some(v) if v.is_finite() => v,
                    _ => {
                        skipped += 1;
                        continue;
                    }
                },
                None => 1.0,
            };
            pts.push(Obs { x, y, time, weight });
        }
        if pts.is_empty() {
            return Err(ToolError::Execution(
                "no usable points (check time_field / coordinates)".to_string(),
            ));
        }

        // ── Time axis: slices every time_step across [t_min, t_max] ─────────────
        let t_min = pts.iter().map(|o| o.time).fold(f64::INFINITY, f64::min);
        let t_max = pts.iter().map(|o| o.time).fold(f64::NEG_INFINITY, f64::max);
        let n_slices = (((t_max - t_min) / prm.time_step).floor() as usize) + 1;
        if n_slices < 2 {
            return Err(ToolError::Execution(format!(
                "only {n_slices} time slice(s) span the data; need >= 2 (reduce time_step)"
            )));
        }
        let slice_time = |b: usize| t_min + b as f64 * prm.time_step;

        // ── Spatial extent (native units) + defaults, padded by the bandwidth ──
        let (mut xmin, mut ymin, mut xmax, mut ymax) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for o in &pts {
            xmin = xmin.min(o.x);
            xmax = xmax.max(o.x);
            ymin = ymin.min(o.y);
            ymax = ymax.max(o.y);
        }
        let lat_mid = 0.5 * (ymin + ymax);
        let coslat = if geographic {
            lat_mid.to_radians().cos().abs().max(1e-6)
        } else {
            1.0
        };
        let (dx, dy) = ((xmax - xmin).max(1e-9), (ymax - ymin).max(1e-9));
        // Longest side expressed in metres, for the default bandwidth.
        let extent_m = if geographic {
            (dx * M_PER_DEG * coslat).max(dy * M_PER_DEG)
        } else {
            dx.max(dy)
        };
        let bandwidth = prm.spatial_bandwidth.unwrap_or((extent_m / 20.0).max(1e-9));
        // Cell size in native units (degrees for geographic).
        let cell = prm.cell_size.unwrap_or((dx.max(dy) / 200.0).max(1e-12));

        // Pad the raster by one bandwidth so kernel footprints are not clipped.
        let (pad_x, pad_y) = if geographic {
            (bandwidth / (M_PER_DEG * coslat), bandwidth / M_PER_DEG)
        } else {
            (bandwidth, bandwidth)
        };
        let x0 = xmin - pad_x;
        let y0 = ymin - pad_y;
        let ext_w = (xmax + pad_x) - x0;
        let ext_h = (ymax + pad_y) - y0;
        let cols = ((ext_w / cell).ceil() as usize).max(1);
        let rows = ((ext_h / cell).ceil() as usize).max(1);
        if rows.saturating_mul(cols).saturating_mul(n_slices) > 300_000_000 {
            return Err(ToolError::Execution(format!(
                "output cube too large ({rows}x{cols}x{n_slices}); increase cell_size or time_step"
            )));
        }

        ctx.progress.info(&format!(
            "{} point(s) -> {n_slices}-band {rows}x{cols} space-time density cube",
            pts.len()
        ));

        // Cell-center coordinate helpers (row 0 = north/top edge).
        let cell_x = |c: usize| x0 + (c as f64 + 0.5) * cell;
        let cell_y = |r: usize| (y0 + ext_h) - (r as f64 + 0.5) * cell;

        // Per-point spatial half-window in cells (bandwidth -> native units).
        let (hw_col, hw_row) = if geographic {
            (
                ((bandwidth / (M_PER_DEG * coslat)) / cell).ceil() as isize,
                ((bandwidth / M_PER_DEG) / cell).ceil() as isize,
            )
        } else {
            let hw = (bandwidth / cell).ceil() as isize;
            (hw, hw)
        };

        // ── Accumulate mass into a band-major cube ─────────────────────────────
        let mut cube = vec![0.0f64; n_slices * rows * cols];
        let mut band_mass = vec![0.0f64; n_slices]; // analytic Σ w_i T_ib
        let mut footprint: Vec<(usize, f64)> = Vec::new();

        for o in &pts {
            // Temporal weights across slices, normalized to sum to 1.
            let mut tw = vec![0.0f64; n_slices];
            let mut tsum = 0.0;
            for (b, twb) in tw.iter_mut().enumerate() {
                let v = prm
                    .temporal_kernel
                    .eval((slice_time(b) - o.time).abs() / prm.temporal_bandwidth);
                *twb = v;
                tsum += v;
            }
            if tsum <= 0.0 {
                // Point beyond the temporal bandwidth of every slice: assign it
                // wholly to the nearest slice so mass is conserved.
                let nearest = (0..n_slices)
                    .min_by(|&a, &b| {
                        (slice_time(a) - o.time)
                            .abs()
                            .total_cmp(&(slice_time(b) - o.time).abs())
                    })
                    .unwrap_or(0);
                tw[nearest] = 1.0;
                tsum = 1.0;
            }
            for twb in tw.iter_mut() {
                *twb /= tsum;
            }

            // Spatial footprint: cells within the bandwidth, normalized to sum 1.
            let cc = (((o.x - x0) / cell).floor() as isize).clamp(0, cols as isize - 1);
            let cr = (((y0 + ext_h - o.y) / cell).floor() as isize).clamp(0, rows as isize - 1);
            footprint.clear();
            let mut ssum = 0.0;
            for r in (cr - hw_row).max(0)..=(cr + hw_row).min(rows as isize - 1) {
                for c in (cc - hw_col).max(0)..=(cc + hw_col).min(cols as isize - 1) {
                    let (r, c) = (r as usize, c as usize);
                    let d = distance(o.x, o.y, cell_x(c), cell_y(r), geographic);
                    let v = prm.spatial_kernel.eval(d / bandwidth);
                    if v > 0.0 {
                        footprint.push((r * cols + c, v));
                        ssum += v;
                    }
                }
            }
            if footprint.is_empty() {
                // Bandwidth narrower than the cell: deposit into the one cell.
                footprint.push((cr as usize * cols + cc as usize, 1.0));
                ssum = 1.0;
            }

            for (b, &twb) in tw.iter().enumerate() {
                if twb <= 0.0 {
                    continue;
                }
                let band_off = b * rows * cols;
                let m = o.weight * twb;
                band_mass[b] += m;
                for &(idx, sv) in &footprint {
                    cube[band_off + idx] += m * (sv / ssum);
                }
            }
        }

        // ── Convert per-cell mass to density (mass / cell area) ────────────────
        // Geographic: metre-scaled area varies with latitude; projected: constant.
        ctx.progress
            .info("normalizing to density and writing bands");
        for r in 0..rows {
            let area = if geographic {
                let lat = cell_y(r);
                (cell * M_PER_DEG * lat.to_radians().cos().abs().max(1e-6)) * (cell * M_PER_DEG)
            } else {
                cell * cell
            };
            let inv = 1.0 / area;
            for b in 0..n_slices {
                let band_off = b * rows * cols + r * cols;
                for c in 0..cols {
                    cube[band_off + c] *= inv;
                }
            }
        }

        // ── Build the multiband raster + band metadata ─────────────────────────
        let nodata = -9999.0f64;
        let crs = match prm.epsg.or_else(|| layer.crs_epsg()) {
            Some(e) => CrsInfo {
                epsg: Some(e),
                wkt: None,
                proj4: None,
            },
            None => CrsInfo {
                epsg: None,
                wkt: None,
                proj4: None,
            },
        };
        let mut metadata: Vec<(String, String)> = Vec::new();
        let mut band_labels: Vec<Value> = Vec::new();
        for (b, &mass) in band_mass.iter().enumerate() {
            let (s, e) = (slice_time(b), slice_time(b) + prm.time_step);
            let label = if time_is_iso {
                format!("{}..{}", format_epoch(s), format_epoch(e))
            } else {
                format!("{s}..{e}")
            };
            metadata.push((format!("band_{}", b + 1), label.clone()));
            band_labels.push(json!({
                "band": b + 1,
                "time_range": label,
                "point_mass": mass,
            }));
        }

        let mut out = Raster::new(RasterConfig {
            cols,
            rows,
            bands: n_slices,
            x_min: x0,
            y_min: y0,
            cell_size: cell,
            cell_size_y: Some(cell),
            nodata,
            data_type: DataType::F32,
            crs,
            metadata,
        });
        for b in 0..n_slices {
            let band_off = b * rows * cols;
            for r in 0..rows {
                for c in 0..cols {
                    out.set(
                        b as isize,
                        r as isize,
                        c as isize,
                        cube[band_off + r * cols + c],
                    )
                    .map_err(|e| ToolError::Execution(format!("write failed: {e}")))?;
                }
            }
        }

        let out_path = write_or_store_output(out, output)?;

        let total_mass: f64 = band_mass.iter().sum();
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("point_count".to_string(), json!(pts.len()));
        outputs.insert("skipped".to_string(), json!(skipped));
        outputs.insert("bands".to_string(), json!(n_slices));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        outputs.insert("cell_size".to_string(), json!(cell));
        outputs.insert("spatial_bandwidth".to_string(), json!(bandwidth));
        outputs.insert("geographic".to_string(), json!(geographic));
        outputs.insert("total_point_mass".to_string(), json!(total_mass));
        outputs.insert("bands_info".to_string(), json!(band_labels));
        Ok(ToolRunResult { outputs })
    }
}

// ── Types & kernels ────────────────────────────────────────────────────────────

struct Obs {
    x: f64,
    y: f64,
    time: f64,
    weight: f64,
}

#[derive(Clone, Copy)]
enum SpatialKernel {
    Epanechnikov,
    Quartic,
}

impl SpatialKernel {
    /// Kernel profile at normalized radius `u = d / bandwidth`; zero for `u >= 1`.
    /// The leading normalization constant is dropped because the spatial
    /// footprint is renormalized to unit mass per point.
    fn eval(&self, u: f64) -> f64 {
        if u.is_nan() || u >= 1.0 {
            return 0.0;
        }
        let a = 1.0 - u * u;
        match self {
            SpatialKernel::Epanechnikov => a,
            SpatialKernel::Quartic => a * a,
        }
    }
}

#[derive(Clone, Copy)]
enum TemporalKernel {
    Triangular,
    Epanechnikov,
}

impl TemporalKernel {
    /// Kernel profile at normalized lag `u = |Δt| / temporal_bandwidth`.
    fn eval(&self, u: f64) -> f64 {
        if u.is_nan() || u >= 1.0 {
            return 0.0;
        }
        match self {
            TemporalKernel::Triangular => 1.0 - u,
            TemporalKernel::Epanechnikov => 1.0 - u * u,
        }
    }
}

/// Distance between two points: haversine metres for geographic input, planar
/// Euclidean in CRS units otherwise. Arguments are (x=lon, y=lat) for geographic.
fn distance(x1: f64, y1: f64, x2: f64, y2: f64, geographic: bool) -> f64 {
    if geographic {
        let (lat1, lat2) = (y1.to_radians(), y2.to_radians());
        let dlat = lat2 - lat1;
        let dlon = (x2 - x1).to_radians();
        let a = (dlat * 0.5).sin().powi(2) + lat1.cos() * lat2.cos() * (dlon * 0.5).sin().powi(2);
        2.0 * EARTH_R * a.sqrt().clamp(0.0, 1.0).asin()
    } else {
        (x2 - x1).hypot(y2 - y1)
    }
}

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        Geometry::LineString(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        Geometry::Polygon { exterior, .. } if !exterior.coords().is_empty() => {
            let c = &exterior.coords()[0];
            Some((c.x, c.y))
        }
        _ => None,
    }
}

// ── Time / value parsing (shared shape with emerging_hot_spot_analysis) ─────────

/// Parses a field value as a time coordinate, returning `(seconds, is_iso)`.
/// A numeric field is used directly; a string is parsed as ISO-8601 (seconds
/// since the epoch).
fn parse_time_value(fv: &FieldValue) -> Option<(f64, bool)> {
    if let Some(n) = fv.as_f64() {
        return Some((n, false));
    }
    fv.as_str()
        .and_then(parse_iso8601_seconds)
        .map(|s| (s, true))
}

/// Minimal ISO-8601 parser: `YYYY-MM-DD` with optional `THH:MM:SS` and a trailing
/// `Z`/offset (offset ignored). Returns seconds since the Unix epoch.
fn parse_iso8601_seconds(s: &str) -> Option<f64> {
    let s = s.trim();
    let bytes = s.as_bytes();
    if bytes.len() < 10 {
        return None;
    }
    let year: i64 = s.get(0..4)?.parse().ok()?;
    if bytes[4] != b'-' {
        return None;
    }
    let month: i64 = s.get(5..7)?.parse().ok()?;
    if bytes[7] != b'-' {
        return None;
    }
    let day: i64 = s.get(8..10)?.parse().ok()?;
    let (mut hh, mut mm, mut ss) = (0i64, 0i64, 0i64);
    if bytes.len() >= 19 && (bytes[10] == b'T' || bytes[10] == b' ') {
        hh = s.get(11..13)?.parse().ok()?;
        mm = s.get(14..16)?.parse().ok()?;
        ss = s.get(17..19)?.parse().ok()?;
    }
    let days = days_from_civil(year, month, day);
    Some((days * 86400 + hh * 3600 + mm * 60 + ss) as f64)
}

/// Days from 1970-01-01 to the given civil date (Howard Hinnant's algorithm).
fn days_from_civil(y: i64, m: i64, d: i64) -> i64 {
    let y = if m <= 2 { y - 1 } else { y };
    let era = if y >= 0 { y } else { y - 399 } / 400;
    let yoe = y - era * 400;
    let doy = (153 * (if m > 2 { m - 3 } else { m + 9 }) + 2) / 5 + d - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146097 + doe - 719468
}

/// Civil date `(year, month, day)` from days since 1970-01-01 (inverse of the
/// above).
fn civil_from_days(z: i64) -> (i64, i64, i64) {
    let z = z + 719468;
    let era = (if z >= 0 { z } else { z - 146096 }) / 146097;
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    (if m <= 2 { y + 1 } else { y }, m, d)
}

/// Formats epoch seconds as `YYYY-MM-DD` (with `THH:MM:SS` when not midnight).
fn format_epoch(secs: f64) -> String {
    let total = secs.floor() as i64;
    let days = total.div_euclid(86400);
    let rem = total.rem_euclid(86400);
    let (y, m, d) = civil_from_days(days);
    if rem == 0 {
        format!("{y:04}-{m:02}-{d:02}")
    } else {
        let (hh, mm, ss) = (rem / 3600, (rem % 3600) / 60, rem % 60);
        format!("{y:04}-{m:02}-{d:02}T{hh:02}:{mm:02}:{ss:02}")
    }
}

/// Reads an optional duration parameter accepting a JSON number (plain units /
/// seconds) or a string (plain number or a duration like `1w`, `12h`).
fn opt_duration(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => parse_duration(s, key).map(Some),
        Some(Value::Number(n)) => match n.as_f64() {
            Some(v) if v > 0.0 && v.is_finite() => Ok(Some(v)),
            _ => Err(ToolError::Validation(format!("'{key}' must be positive"))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "'{key}' must be a number or duration string"
        ))),
    }
}

/// Parses `time_step` / `temporal_bandwidth`: a plain number or a duration like
/// `1w`, `7d`, `12h`, `30m`, `1M`, `1y` (returned in seconds; M≈30d, y≈365d).
fn parse_duration(s: &str, key: &str) -> Result<f64, ToolError> {
    let s = s.trim();
    if let Ok(v) = s.parse::<f64>() {
        if v > 0.0 && v.is_finite() {
            return Ok(v);
        }
        return Err(ToolError::Validation(format!("'{key}' must be positive")));
    }
    let (num, unit) = s.split_at(s.len().saturating_sub(1));
    let value: f64 = num
        .trim()
        .parse()
        .map_err(|_| ToolError::Validation(format!("could not parse '{key}' value in '{s}'")))?;
    if !(value > 0.0 && value.is_finite()) {
        return Err(ToolError::Validation(format!("'{key}' must be positive")));
    }
    let seconds = match unit {
        "s" => 1.0,
        "m" => 60.0,
        "h" => 3600.0,
        "d" => 86400.0,
        "w" => 604800.0,
        "M" => 2_592_000.0,
        "y" => 31_536_000.0,
        other => {
            return Err(ToolError::Validation(format!(
                "unknown '{key}' unit '{other}' (use s/m/h/d/w/M/y or a plain number)"
            )))
        }
    };
    Ok(value * seconds)
}

// ── Parameters ──────────────────────────────────────────────────────────────────

struct Params {
    time_field: String,
    time_step: f64,
    temporal_bandwidth: f64,
    spatial_bandwidth: Option<f64>,
    cell_size: Option<f64>,
    weight_field: Option<String>,
    spatial_kernel: SpatialKernel,
    temporal_kernel: TemporalKernel,
    epsg: Option<u32>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let time_field = require_str(args, "time_field")?.to_string();
    let time_step = opt_duration(args, "time_step")?.unwrap_or(604800.0); // '1w'
    let temporal_bandwidth = opt_duration(args, "temporal_bandwidth")?.unwrap_or(time_step);
    let spatial_bandwidth = opt_pos_f64(args, "spatial_bandwidth")?;
    let cell_size = opt_pos_f64(args, "cell_size")?;
    let weight_field = parse_optional_str(args, "weight_field")?.map(str::to_string);
    let spatial_kernel = match args
        .get("spatial_kernel")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some("epanechnikov") => SpatialKernel::Epanechnikov,
        Some("quartic") | Some("biweight") => SpatialKernel::Quartic,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'spatial_kernel' must be epanechnikov/quartic, got '{o}'"
            )))
        }
    };
    let temporal_kernel = match args
        .get("temporal_kernel")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some("triangular") => TemporalKernel::Triangular,
        Some("epanechnikov") => TemporalKernel::Epanechnikov,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'temporal_kernel' must be triangular/epanechnikov, got '{o}'"
            )))
        }
    };
    let epsg = opt_u32(args, "epsg")?;
    Ok(Params {
        time_field,
        time_step,
        temporal_bandwidth,
        spatial_bandwidth,
        cell_size,
        weight_field,
        spatial_kernel,
        temporal_kernel,
        epsg,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

/// Parses an optional positive f64 accepting a JSON number or numeric string.
fn opt_pos_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    let v = match args.get(key) {
        None | Some(Value::Null) => return Ok(None),
        Some(Value::Number(n)) => n.as_f64(),
        Some(Value::String(s)) if s.trim().is_empty() => return Ok(None),
        Some(Value::String(s)) => s.trim().parse::<f64>().ok(),
        Some(_) => None,
    };
    match v {
        Some(v) if v > 0.0 && v.is_finite() => Ok(Some(v)),
        _ => Err(ToolError::Validation(format!(
            "'{key}' must be a positive number"
        ))),
    }
}

fn opt_u32(args: &ToolArgs, key: &str) -> Result<Option<u32>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<u32>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("'{key}' must be an integer"))),
        Some(Value::Number(n)) => n
            .as_u64()
            .filter(|v| *v <= u32::MAX as u64)
            .map(|v| Some(v as u32))
            .ok_or_else(|| ToolError::Validation(format!("'{key}' must be an integer"))),
        Some(_) => Err(ToolError::Validation(format!(
            "'{key}' must be a number when provided"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::load_input_raster;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Sum of all cells of a band times the (constant, projected) cell area.
    fn band_mass(r: &Raster, band: usize, cell_area: f64) -> f64 {
        let mut s = 0.0;
        for row in 0..r.rows {
            for col in 0..r.cols {
                let v = r.get(band as isize, row as isize, col as isize);
                if v != r.nodata {
                    s += v;
                }
            }
        }
        s * cell_area
    }

    fn point_layer(pts: &[(f64, f64, f64, f64)], epsg: u32, weighted: bool) -> String {
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(epsg);
        l.add_field(FieldDef::new("t", FieldType::Float));
        l.add_field(FieldDef::new("w", FieldType::Float));
        for (x, y, t, w) in pts {
            let _ = weighted;
            l.add_feature(
                Some(Geometry::point(*x, *y)),
                &[("t", (*t).into()), ("w", (*w).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(v: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(v).unwrap();
        let out = SpaceTimeKernelDensityTool.run(&args, &ctx()).unwrap();
        let r = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, r)
    }

    /// Core property: with a projected grid (constant cell area), the integrated
    /// mass of each band equals the analytic Σ weights in its temporal window, and
    /// the total mass across bands equals the total point weight.
    #[test]
    fn mass_is_conserved_per_band() {
        // Three points, one per time step (t = 0, 10, 20), unit weight.
        let input = point_layer(
            &[
                (0.0, 0.0, 0.0, 1.0),
                (100.0, 0.0, 10.0, 1.0),
                (0.0, 100.0, 20.0, 1.0),
            ],
            32633,
            false,
        );
        let (out, r) = run(json!({
            "input": input,
            "time_field": "t",
            "time_step": "10",
            "temporal_bandwidth": "10",
            "spatial_bandwidth": "30",
            "cell_size": "5",
        }));
        assert_eq!(out.outputs["bands"], json!(3));
        assert!(!out.outputs["geographic"].as_bool().unwrap());
        let cell_area = 5.0 * 5.0;
        // Each of the 3 slices should carry ~1 unit of mass (one point each).
        let mut total = 0.0;
        for b in 0..3 {
            let m = band_mass(&r, b, cell_area);
            total += m;
            assert!(
                (m - 1.0).abs() < 1e-3,
                "band {b} mass {m} should be ~1 (one point in the window)"
            );
        }
        assert!(
            (total - 3.0).abs() < 1e-3,
            "total mass {total} should equal total weight 3"
        );
        // The tool's own analytic accounting must agree.
        assert!((out.outputs["total_point_mass"].as_f64().unwrap() - 3.0).abs() < 1e-9);
    }

    /// A weight field scales each point's deposited mass.
    #[test]
    fn weight_field_scales_mass() {
        let input = point_layer(
            &[(0.0, 0.0, 0.0, 2.0), (100.0, 0.0, 10.0, 5.0)],
            32633,
            true,
        );
        let (out, r) = run(json!({
            "input": input,
            "time_field": "t",
            "weight_field": "w",
            "time_step": "10",
            "temporal_bandwidth": "10",
            "spatial_bandwidth": "30",
            "cell_size": "5",
        }));
        let cell_area = 25.0;
        assert!((band_mass(&r, 0, cell_area) - 2.0).abs() < 1e-3);
        assert!((band_mass(&r, 1, cell_area) - 5.0).abs() < 1e-3);
        assert!((out.outputs["total_point_mass"].as_f64().unwrap() - 7.0).abs() < 1e-9);
    }

    /// Temporal smoothing: a point midway between two slices splits its mass
    /// roughly evenly between them (triangular kernel).
    #[test]
    fn point_between_slices_splits_temporally() {
        // Points at t=0 and t=20 anchor the axis to two steps of width 10, so
        // slices sit at t=0,10,20. A point at t=5 is midway between slice 0 and 1.
        let input = point_layer(
            &[
                (0.0, 0.0, 0.0, 1.0),
                (0.0, 0.0, 20.0, 1.0),
                (0.0, 0.0, 5.0, 1.0),
            ],
            32633,
            false,
        );
        let (_out, r) = run(json!({
            "input": input,
            "time_field": "t",
            "time_step": "10",
            "temporal_bandwidth": "10",
            "spatial_bandwidth": "30",
            "cell_size": "5",
        }));
        let ca = 25.0;
        let m0 = band_mass(&r, 0, ca);
        let m1 = band_mass(&r, 1, ca);
        // slice0: point@0 (1.0) + half of point@5 (~0.5) = ~1.5
        // slice1: half of point@5 (~0.5) = ~0.5
        assert!((m0 - 1.5).abs() < 1e-2, "band0 mass {m0} expected ~1.5");
        assert!((m1 - 0.5).abs() < 1e-2, "band1 mass {m1} expected ~0.5");
    }

    /// Geographic input goes through the haversine path and still conserves total
    /// mass (checked via the tool's analytic accounting).
    #[test]
    fn geographic_input_conserves_total_mass() {
        let input = point_layer(
            &[
                (-77.0, 39.0, 0.0, 1.0),
                (-77.01, 39.01, 10.0, 1.0),
                (-76.99, 38.99, 20.0, 1.0),
            ],
            4326,
            false,
        );
        let (out, _r) = run(json!({
            "input": input,
            "time_field": "t",
            "time_step": "10",
            "temporal_bandwidth": "10",
        }));
        assert!(out.outputs["geographic"].as_bool().unwrap());
        assert!((out.outputs["total_point_mass"].as_f64().unwrap() - 3.0).abs() < 1e-9);
    }

    #[test]
    fn iso_timestamps_and_band_labels() {
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("ts", FieldType::Text));
        for wk in 0..4 {
            l.add_feature(
                Some(Geometry::point(-77.0 + wk as f64 * 0.001, 39.0)),
                &[("ts", format!("2024-01-{:02}T00:00:00Z", 1 + wk * 7).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input,
            "time_field": "ts",
            "time_step": "1w",
            "temporal_bandwidth": "1w",
        }))
        .unwrap();
        let out = SpaceTimeKernelDensityTool.run(&args, &ctx()).unwrap();
        assert_eq!(out.outputs["bands"], json!(4));
        let info = out.outputs["bands_info"].as_array().unwrap();
        // ISO labels should be date-formatted, not raw epoch seconds.
        let lbl = info[0]["time_range"].as_str().unwrap();
        assert!(lbl.starts_with("2024-01-01"), "unexpected band label {lbl}");
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SpaceTimeKernelDensityTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err()); // no time_field
        assert!(bad(json!({ "input": "a.geojson", "time_field": "t", "time_step": "0" })).is_err());
        assert!(bad(
            json!({ "input": "a.geojson", "time_field": "t", "spatial_kernel": "gaussian" })
        )
        .is_err());
        assert!(
            bad(json!({ "input": "a.geojson", "time_field": "t", "spatial_bandwidth": "-5" }))
                .is_err()
        );
        assert!(bad(json!({ "input": "a.geojson", "time_field": "t" })).is_ok());
    }
}
