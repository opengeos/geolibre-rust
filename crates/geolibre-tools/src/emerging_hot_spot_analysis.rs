//! GeoLibre tool: emerging hot spot analysis on an H3 space-time cube.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Emerging Hot Spot Analysis*
//! (Space Time Pattern Mining), combined with *Create Space Time Cube*. The
//! bundled `getis_ord_gi_star` is spatial-only — it has no notion of time. This
//! tool adds the temporal dimension:
//!
//! 1. **Bin** timestamped points into a space-time cube: an H3 cell (spatial
//!    axis, via [`h3o`], as in `vector_to_h3`) × a time step (temporal axis).
//!    Each bin holds a count of points, or the sum of a value field.
//! 2. **Gi\*** per bin: the Getis-Ord Gi\* z-score at every (cell, time) using a
//!    space-time neighborhood — the spatial k-ring of the cell crossed with a
//!    window of ± time steps — against the mean and standard deviation of the
//!    whole cube. A high positive z means the bin sits in a space-time cluster
//!    of high values (a hot spot); a low negative z, a cold spot.
//! 3. **Trend**: a Mann-Kendall test on each cell's time series of Gi\* z-scores
//!    tells whether its clustering is intensifying or diminishing over time.
//! 4. **Classify** each cell into the Esri categories — new / consecutive /
//!    intensifying / persistent / diminishing / sporadic / oscillating /
//!    historical hot (and cold) spot, or no pattern — from the pattern of
//!    significant bins plus the trend.
//!
//! The whole computation is deterministic: no Monte-Carlo, no RNG. Output is one
//! H3 polygon per cell (renderable straight through `h3_to_vector`/PMTiles) with
//! its category, latest Gi\* z-score, and Mann-Kendall z/p.
//!
//! The `time_field` may be a numeric field (its own units) or an ISO-8601
//! date/datetime string (parsed to seconds); `time_step` is then a plain number
//! or a duration like `1w`, `7d`, `12h`. The study area is the set of occupied
//! cells; empty bins inside it count as zeros (as a space-time cube requires).

use std::collections::{BTreeMap, BTreeSet, HashMap};

use h3o::{CellIndex, LatLng, Resolution};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Default H3 resolution: 7 (~5 km² cells, ~1.4 km edge) suits city-scale
/// point sets where multiple points share a cell each time step.
const DEFAULT_RESOLUTION: u8 = 7;
/// Gi* significance threshold (two-sided 95%).
const Z_CRIT: f64 = 1.96;
/// A cell is "persistently" hot/cold when at least this fraction of its time
/// steps are significant (Esri's 90%).
const PERSIST_FRACTION: f64 = 0.9;

pub struct EmergingHotSpotAnalysisTool;

impl Tool for EmergingHotSpotAnalysisTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "emerging_hot_spot_analysis",
            display_name: "Emerging Hot Spot Analysis",
            summary: "Space-time hot spot trends on an H3 space-time cube: bin timestamped points into H3 cells x time steps, compute Getis-Ord Gi* per bin, and classify each cell (new/intensifying/persistent/diminishing/... hot or cold spot) via a Mann-Kendall trend test, like ArcGIS Emerging Hot Spot Analysis.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer (geographic lon/lat, EPSG:4326). Other geometries use their vertex-mean representative point.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (driver from extension). One H3 polygon per cell. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "time_field",
                    description: "Field holding each point's time: a numeric field (its own units) or an ISO-8601 date/datetime string.",
                    required: true,
                },
                ToolParamSpec {
                    name: "time_step",
                    description: "Width of a time step: a plain number (numeric time_field) or a duration like '1w', '7d', '12h', '1M' for a date field. Default '1w'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "value_field",
                    description: "Optional numeric field summed into each bin. Default: count of points per bin.",
                    required: false,
                },
                ToolParamSpec {
                    name: "resolution",
                    description: "H3 resolution 0 (coarsest) to 15 (finest); default 7 (~5 km² cells).",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighborhood",
                    description: "Spatial neighborhood radius in H3 k-rings (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "time_window",
                    description: "Temporal neighborhood radius in time steps (default 1).",
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
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        ctx.progress.info("reading input points");
        let layer = load_input_layer(input)?;
        if let Some(epsg) = layer.crs_epsg() {
            if epsg != 4326 {
                return Err(ToolError::Validation(format!(
                    "input CRS is EPSG:{epsg}; reproject to EPSG:4326 (lon/lat) before H3 binning"
                )));
            }
        }
        let time_idx = layer.schema.field_index(&prm.time_field).ok_or_else(|| {
            ToolError::Validation(format!("time_field '{}' not found", prm.time_field))
        })?;
        let value_idx =
            match &prm.value_field {
                Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                    ToolError::Validation(format!("value_field '{f}' not found"))
                })?),
                None => None,
            };

        // ── Collect (cell, time, value) observations ──────────────────────────
        ctx.progress.info("binning points into the space-time cube");
        let mut obs: Vec<(CellIndex, f64, f64)> = Vec::new();
        let mut skipped = 0u64;
        for feature in layer.iter() {
            let Some(geom) = feature.geometry.as_ref() else {
                skipped += 1;
                continue;
            };
            let Some((lng, lat)) = representative_lnglat(geom) else {
                skipped += 1;
                continue;
            };
            let Some(time) = feature.attributes.get(time_idx).and_then(parse_time_value) else {
                skipped += 1;
                continue;
            };
            let value = match value_idx {
                Some(vi) => match feature.attributes.get(vi).and_then(FieldValue::as_f64) {
                    Some(v) => v,
                    None => {
                        skipped += 1;
                        continue;
                    }
                },
                None => 1.0,
            };
            match LatLng::new(lat, lng) {
                Ok(ll) => obs.push((ll.to_cell(prm.resolution), time, value)),
                Err(_) => skipped += 1,
            }
        }
        if obs.is_empty() {
            return Err(ToolError::Execution(
                "no usable observations (check time_field / coordinates)".to_string(),
            ));
        }

        // ── Time axis: bin each observation into a step index 0..n_times ──────
        let t_min = obs.iter().map(|o| o.1).fold(f64::INFINITY, f64::min);
        let t_max = obs.iter().map(|o| o.1).fold(f64::NEG_INFINITY, f64::max);
        let n_times = (((t_max - t_min) / prm.time_step).floor() as usize) + 1;
        if n_times < 2 {
            return Err(ToolError::Execution(format!(
                "only {n_times} time step(s) span the data; need >= 2 (reduce time_step)"
            )));
        }
        let time_bin =
            |t: f64| -> usize { (((t - t_min) / prm.time_step).floor() as usize).min(n_times - 1) };

        // ── Cube: study cells × times, summing values (missing bins = 0) ──────
        let cells: Vec<CellIndex> = {
            let set: BTreeSet<u64> = obs.iter().map(|o| u64::from(o.0)).collect();
            set.into_iter()
                .map(|r| CellIndex::try_from(r).unwrap())
                .collect()
        };
        let cell_pos: HashMap<CellIndex, usize> =
            cells.iter().enumerate().map(|(i, &c)| (c, i)).collect();
        let n_cells = cells.len();
        let mut cube = vec![0.0f64; n_cells * n_times];
        for &(cell, time, value) in &obs {
            let ci = cell_pos[&cell];
            let ti = time_bin(time);
            cube[ci * n_times + ti] += value;
        }

        ctx.progress.info(&format!(
            "cube {n_cells} cell(s) x {n_times} time step(s); computing space-time Gi*"
        ));

        // ── Global mean / std over all cube bins (the Gi* reference) ──────────
        let n = (n_cells * n_times) as f64;
        let sum: f64 = cube.iter().sum();
        let mean = sum / n;
        let var = cube.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
        let std = var.sqrt();

        // Precompute each study cell's in-study k-ring neighbors.
        let neighbors: Vec<Vec<usize>> = cells
            .iter()
            .map(|&c| {
                c.grid_disk::<Vec<CellIndex>>(prm.neighborhood)
                    .into_iter()
                    .filter_map(|nc| cell_pos.get(&nc).copied())
                    .collect()
            })
            .collect();

        // ── Gi* z-score per (cell, time) ──────────────────────────────────────
        let mut giz = vec![0.0f64; n_cells * n_times];
        for ci in 0..n_cells {
            for ti in 0..n_times {
                let t0 = ti.saturating_sub(prm.time_window);
                let t1 = (ti + prm.time_window).min(n_times - 1);
                let mut local_sum = 0.0;
                let mut w = 0.0;
                for &nj in &neighbors[ci] {
                    for ts in t0..=t1 {
                        local_sum += cube[nj * n_times + ts];
                        w += 1.0;
                    }
                }
                giz[ci * n_times + ti] = gi_star(local_sum, w, mean, std, n);
            }
        }

        // ── Per-cell Mann-Kendall trend + Esri classification ─────────────────
        ctx.progress.info("classifying cells");
        let mut categories: BTreeMap<&'static str, usize> = BTreeMap::new();
        let mut out = Layer::new("emerging_hot_spots")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        out.add_field(FieldDef::new("h3", FieldType::Text));
        out.add_field(FieldDef::new("category", FieldType::Text));
        out.add_field(FieldDef::new("gi_z_last", FieldType::Float));
        out.add_field(FieldDef::new("mk_z", FieldType::Float));
        out.add_field(FieldDef::new("mk_p", FieldType::Float));
        out.add_field(FieldDef::new("n_hot", FieldType::Integer));
        out.add_field(FieldDef::new("n_cold", FieldType::Integer));

        for ci in 0..n_cells {
            let series = &giz[ci * n_times..ci * n_times + n_times];
            let mk = mann_kendall(series);
            let trend_up = mk.p < 0.05 && mk.z > 0.0;
            let trend_down = mk.p < 0.05 && mk.z < 0.0;
            let category = classify(series, trend_up, trend_down);
            *categories.entry(category).or_default() += 1;

            let n_hot = series.iter().filter(|z| **z > Z_CRIT).count() as i64;
            let n_cold = series.iter().filter(|z| **z < -Z_CRIT).count() as i64;
            let ring = cell_polygon_ring(cells[ci]);
            out.add_feature(
                Some(Geometry::polygon(ring, Vec::new())),
                &[
                    ("h3", cells[ci].to_string().into()),
                    ("category", category.into()),
                    ("gi_z_last", series[n_times - 1].into()),
                    ("mk_z", mk.z.into()),
                    ("mk_p", mk.p.into()),
                    ("n_hot", n_hot.into()),
                    ("n_cold", n_cold.into()),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed building cell feature: {e}")))?;
        }

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("cell_count".to_string(), json!(n_cells));
        outputs.insert("time_steps".to_string(), json!(n_times));
        outputs.insert("skipped".to_string(), json!(skipped));
        for (cat, count) in &categories {
            outputs.insert(format!("category_{cat}"), json!(count));
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── Getis-Ord Gi* ──────────────────────────────────────────────────────────

/// Gi* z-score for a neighborhood: `local_sum` is the sum of values over the
/// `w` space-time neighbor bins; `mean`/`std` are over all `n` cube bins. Binary
/// weights, so the S1 term equals `w`.
fn gi_star(local_sum: f64, w: f64, mean: f64, std: f64, n: f64) -> f64 {
    if std <= 0.0 || w <= 0.0 || n <= 1.0 {
        return 0.0;
    }
    let denom = std * ((n * w - w * w) / (n - 1.0)).max(0.0).sqrt();
    if denom <= 0.0 {
        return 0.0;
    }
    (local_sum - mean * w) / denom
}

// ── Mann-Kendall trend test ─────────────────────────────────────────────────

struct MannKendall {
    z: f64,
    p: f64,
}

/// Mann-Kendall trend test with tie correction on the continuity-corrected z.
fn mann_kendall(series: &[f64]) -> MannKendall {
    let n = series.len();
    if n < 3 {
        return MannKendall { z: 0.0, p: 1.0 };
    }
    let mut s = 0i64;
    for i in 0..n {
        for j in i + 1..n {
            s += (series[j] - series[i]).signum() as i64;
        }
    }
    // Variance with a tie correction over groups of equal values.
    let mut tie_term = 0.0;
    let mut counts: HashMap<u64, u64> = HashMap::new();
    for v in series {
        *counts.entry(v.to_bits()).or_default() += 1;
    }
    for &t in counts.values() {
        if t > 1 {
            let t = t as f64;
            tie_term += t * (t - 1.0) * (2.0 * t + 5.0);
        }
    }
    let nf = n as f64;
    let var = (nf * (nf - 1.0) * (2.0 * nf + 5.0) - tie_term) / 18.0;
    if var <= 0.0 {
        return MannKendall { z: 0.0, p: 1.0 };
    }
    let z = if s > 0 {
        (s as f64 - 1.0) / var.sqrt()
    } else if s < 0 {
        (s as f64 + 1.0) / var.sqrt()
    } else {
        0.0
    };
    let p = 2.0 * (1.0 - normal_cdf(z.abs()));
    MannKendall { z, p }
}

/// Standard normal CDF via the erf rational approximation (A&S 7.1.26).
fn normal_cdf(x: f64) -> f64 {
    0.5 * erfc(-x / std::f64::consts::SQRT_2)
}

fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let ans = t
        * (-z * z - 1.26551223
            + t * (1.00002368
                + t * (0.37409196
                    + t * (0.09678418
                        + t * (-0.18628806
                            + t * (0.27886807
                                + t * (-1.13520398
                                    + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
            .exp();
    if x >= 0.0 {
        ans
    } else {
        2.0 - ans
    }
}

// ── Esri emerging hot-spot classification ────────────────────────────────────

/// Classify a cell's Gi* z-score time series into an Esri emerging hot/cold-spot
/// category, given whether the Mann-Kendall trend rises or falls significantly.
fn classify(z: &[f64], trend_up: bool, trend_down: bool) -> &'static str {
    let t = z.len();
    if t == 0 {
        return "no_pattern";
    }
    let hot: Vec<bool> = z.iter().map(|v| *v > Z_CRIT).collect();
    let cold: Vec<bool> = z.iter().map(|v| *v < -Z_CRIT).collect();
    let nhot = hot.iter().filter(|b| **b).count();
    let ncold = cold.iter().filter(|b| **b).count();
    if nhot == 0 && ncold == 0 {
        return "no_pattern";
    }
    let final_hot = hot[t - 1];
    let final_cold = cold[t - 1];
    // Polarity: the final step decides; else whichever type is more present.
    let hot_side = if final_hot {
        true
    } else if final_cold {
        false
    } else {
        nhot >= ncold
    };
    if hot_side {
        // For hot spots, intensifying means the trend rises.
        classify_side(&hot, &cold, nhot, trend_up, trend_down, true)
    } else {
        // For cold spots, "intensifying" means the trend falls (more negative).
        classify_side(&cold, &hot, ncold, trend_down, trend_up, false)
    }
}

#[allow(clippy::too_many_arguments)]
fn classify_side(
    same: &[bool],
    opp: &[bool],
    nsame: usize,
    trend_intensify: bool,
    trend_diminish: bool,
    hot: bool,
) -> &'static str {
    let t = same.len();
    let pct = nsame as f64 / t as f64;
    let final_same = same[t - 1];
    let ever_opp = opp.iter().any(|b| *b);
    let ever_same_before = same[..t - 1].iter().any(|b| *b);

    if final_same {
        if ever_opp {
            return pick(hot, "oscillating_hot_spot", "oscillating_cold_spot");
        }
        if pct >= PERSIST_FRACTION {
            if trend_intensify {
                return pick(hot, "intensifying_hot_spot", "intensifying_cold_spot");
            }
            if trend_diminish {
                return pick(hot, "diminishing_hot_spot", "diminishing_cold_spot");
            }
            return pick(hot, "persistent_hot_spot", "persistent_cold_spot");
        }
        if !ever_same_before {
            return pick(hot, "new_hot_spot", "new_cold_spot");
        }
        if consecutive_run_to_end(same) {
            return pick(hot, "consecutive_hot_spot", "consecutive_cold_spot");
        }
        pick(hot, "sporadic_hot_spot", "sporadic_cold_spot")
    } else {
        // Final step is not this type.
        if pct >= PERSIST_FRACTION {
            return pick(hot, "historical_hot_spot", "historical_cold_spot");
        }
        if !ever_opp && nsame > 0 {
            return pick(hot, "sporadic_hot_spot", "sporadic_cold_spot");
        }
        "no_pattern"
    }
}

fn pick(hot: bool, hot_cat: &'static str, cold_cat: &'static str) -> &'static str {
    if hot {
        hot_cat
    } else {
        cold_cat
    }
}

/// True when the trailing run of `true`s reaches the end and there is no earlier
/// `true` before that run (a single uninterrupted run ending at the final step).
fn consecutive_run_to_end(same: &[bool]) -> bool {
    let t = same.len();
    if !same[t - 1] {
        return false;
    }
    let mut i = t;
    while i > 0 && same[i - 1] {
        i -= 1;
    }
    !same[..i].iter().any(|b| *b)
}

// ── Time / value parsing ─────────────────────────────────────────────────────

/// Parses a field value as a time coordinate: a numeric field is used directly;
/// a string is parsed as an ISO-8601 date/datetime and returned as seconds since
/// the Unix epoch.
fn parse_time_value(fv: &FieldValue) -> Option<f64> {
    if let Some(n) = fv.as_f64() {
        return Some(n);
    }
    fv.as_str().and_then(parse_iso8601_seconds)
}

/// Minimal ISO-8601 parser: `YYYY-MM-DD` with an optional `THH:MM:SS` and a
/// trailing `Z` or offset (offset ignored). Returns seconds since the epoch.
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

/// Parses `time_step`: a plain number (numeric field units) or a duration like
/// `1w`, `7d`, `12h`, `30m`, `1M`, `1y` (returned in seconds; M≈30d, y≈365d).
fn parse_time_step(s: &str) -> Result<f64, ToolError> {
    let s = s.trim();
    if let Ok(v) = s.parse::<f64>() {
        if v > 0.0 && v.is_finite() {
            return Ok(v);
        }
        return Err(ToolError::Validation(
            "'time_step' must be a positive number".to_string(),
        ));
    }
    let (num, unit) = s.split_at(s.len().saturating_sub(1));
    let value: f64 = num.trim().parse().map_err(|_| {
        ToolError::Validation(format!("could not parse 'time_step' value in '{s}'"))
    })?;
    if !(value > 0.0 && value.is_finite()) {
        return Err(ToolError::Validation(
            "'time_step' must be positive".to_string(),
        ));
    }
    let seconds = match unit {
        "s" => 1.0,
        "m" => 60.0,
        "h" => 3600.0,
        "d" => 86400.0,
        "w" => 604800.0,
        "M" => 2_592_000.0,  // 30 days
        "y" => 31_536_000.0, // 365 days
        other => {
            return Err(ToolError::Validation(format!(
                "unknown 'time_step' unit '{other}' (use s/m/h/d/w/M/y or a plain number)"
            )))
        }
    };
    Ok(value * seconds)
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    time_field: String,
    time_step: f64,
    value_field: Option<String>,
    resolution: Resolution,
    neighborhood: u32,
    time_window: usize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let time_field = require_str(args, "time_field")?.to_string();
    let time_step = match parse_optional_str(args, "time_step")? {
        Some(s) => parse_time_step(s)?,
        None => 604800.0, // 1 week, matching the default '1w'
    };
    let value_field = parse_optional_str(args, "value_field")?.map(str::to_string);
    let res_u8 = parse_optional_u32(args, "resolution")?
        .map(|v| v as u8)
        .unwrap_or(DEFAULT_RESOLUTION);
    let resolution = Resolution::try_from(res_u8)
        .map_err(|_| ToolError::Validation(format!("'resolution' must be 0-15, got {res_u8}")))?;
    let neighborhood = parse_optional_u32(args, "neighborhood")?.unwrap_or(1);
    let time_window = parse_optional_u32(args, "time_window")?.unwrap_or(1) as usize;
    Ok(Params {
        time_field,
        time_step,
        value_field,
        resolution,
        neighborhood,
        time_window,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn parse_optional_u32(args: &ToolArgs, key: &str) -> Result<Option<u32>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => {
            s.trim().parse::<u32>().map(Some).map_err(|_| {
                ToolError::Validation(format!("'{key}' must be a non-negative integer"))
            })
        }
        Some(Value::Number(n)) => n
            .as_u64()
            .filter(|v| *v <= u32::MAX as u64)
            .map(|v| Some(v as u32))
            .ok_or_else(|| {
                ToolError::Validation(format!("'{key}' must be a non-negative integer"))
            }),
        Some(_) => Err(ToolError::Validation(format!(
            "'{key}' must be a number when provided"
        ))),
    }
}

// ── Geometry helpers (shared shape with vector_to_h3) ─────────────────────────

fn representative_lnglat(geom: &Geometry) -> Option<(f64, f64)> {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0u64;
    accumulate_coords(geom, &mut sx, &mut sy, &mut n);
    (n > 0).then(|| (sx / n as f64, sy / n as f64))
}

fn accumulate_coords(geom: &Geometry, sx: &mut f64, sy: &mut f64, n: &mut u64) {
    let mut add = |c: &Coord| {
        *sx += c.x;
        *sy += c.y;
        *n += 1;
    };
    match geom {
        Geometry::Point(c) => add(c),
        Geometry::LineString(cs) | Geometry::MultiPoint(cs) => cs.iter().for_each(add),
        Geometry::MultiLineString(lines) => lines.iter().flatten().for_each(add),
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            exterior.coords().iter().for_each(&mut add);
            interiors
                .iter()
                .for_each(|r| r.coords().iter().for_each(&mut add));
        }
        Geometry::MultiPolygon(polys) => {
            for (ext, holes) in polys {
                ext.coords().iter().for_each(&mut add);
                holes
                    .iter()
                    .for_each(|r| r.coords().iter().for_each(&mut add));
            }
        }
        Geometry::GeometryCollection(geoms) => {
            for g in geoms {
                accumulate_coords(g, sx, sy, n);
            }
        }
    }
}

fn cell_polygon_ring(cell: CellIndex) -> Vec<Coord> {
    cell.boundary()
        .iter()
        .map(|ll| Coord::xy(ll.lng(), ll.lat()))
        .collect()
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

    #[test]
    fn iso8601_parses_dates() {
        assert_eq!(parse_iso8601_seconds("1970-01-01"), Some(0.0));
        assert_eq!(parse_iso8601_seconds("1970-01-02"), Some(86400.0));
        assert_eq!(
            parse_iso8601_seconds("2024-01-01T00:00:00Z"),
            Some(1_704_067_200.0)
        );
        // one week apart
        let a = parse_iso8601_seconds("2024-01-01T00:00:00Z").unwrap();
        let b = parse_iso8601_seconds("2024-01-08T00:00:00Z").unwrap();
        assert_eq!(b - a, 604800.0);
        assert_eq!(parse_iso8601_seconds("not-a-date"), None);
    }

    #[test]
    fn time_step_parses_durations() {
        assert_eq!(parse_time_step("1w").unwrap(), 604800.0);
        assert_eq!(parse_time_step("7d").unwrap(), 604800.0);
        assert_eq!(parse_time_step("12h").unwrap(), 43200.0);
        assert_eq!(parse_time_step("3").unwrap(), 3.0);
        assert!(parse_time_step("0").is_err());
        assert!(parse_time_step("1x").is_err());
    }

    #[test]
    fn mann_kendall_detects_monotone_trend() {
        let up = mann_kendall(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert!(
            up.z > 0.0 && up.p < 0.05,
            "expected significant up trend: {:?}",
            (up.z, up.p)
        );
        let down = mann_kendall(&[6.0, 5.0, 4.0, 3.0, 2.0, 1.0]);
        assert!(down.z < 0.0 && down.p < 0.05);
        let flat = mann_kendall(&[1.0, 1.0, 1.0, 1.0]);
        assert!(flat.p > 0.05);
    }

    #[test]
    fn classify_categories() {
        let h = 3.0; // significant hot
        let c = -3.0; // significant cold
        let o = 0.0; // not significant
                     // New: hot only in the final step.
        assert_eq!(
            classify(&[o, o, o, o, o, o, o, o, o, o, o, h], false, false),
            "new_hot_spot"
        );
        // Intensifying: hot throughout with a rising trend.
        assert_eq!(
            classify(&[h, h, h, h, h, h, h, h, h, h, h, h], true, false),
            "intensifying_hot_spot"
        );
        // Persistent: hot throughout, no trend.
        assert_eq!(
            classify(&[h, h, h, h, h, h, h, h, h, h, h, h], false, false),
            "persistent_hot_spot"
        );
        // Diminishing: hot throughout, falling trend.
        assert_eq!(
            classify(&[h, h, h, h, h, h, h, h, h, h, h, h], false, true),
            "diminishing_hot_spot"
        );
        // Historical: hot most of the way but not at the end.
        assert_eq!(
            classify(&[h, h, h, h, h, h, h, h, h, h, h, o], false, false),
            "historical_hot_spot"
        );
        // Oscillating: hot at the end, cold earlier.
        assert_eq!(
            classify(&[c, c, o, o, o, o, o, o, o, o, o, h], false, false),
            "oscillating_hot_spot"
        );
        // No pattern: never significant.
        assert_eq!(classify(&[o, o, o, o, o, o], false, false), "no_pattern");
        // Cold side mirrors.
        assert_eq!(
            classify(&[o, o, o, o, o, o, o, o, o, o, o, c], false, false),
            "new_cold_spot"
        );
    }

    #[test]
    fn end_to_end_on_a_synthetic_cube() {
        // A tight cluster of points whose count grows every week (should read as
        // an intensifying hot spot) plus a scatter of background points.
        let mut layer = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        layer.add_field(FieldDef::new("t", FieldType::Integer));
        // Hot cluster near (0,0): more points each week.
        for week in 0..12 {
            for k in 0..(2 + week * 3) {
                let jitter = k as f64 * 1e-4;
                layer
                    .add_feature(
                        Some(Geometry::point(0.001 + jitter, 0.001)),
                        &[("t", (week as i64).into())],
                    )
                    .unwrap();
            }
        }
        // Background scatter far away, steady low counts.
        for week in 0..12 {
            layer
                .add_feature(
                    Some(Geometry::point(0.5 + week as f64 * 1e-3, 0.5)),
                    &[("t", (week as i64).into())],
                )
                .unwrap();
        }
        let id = wbvector::memory_store::put_vector(layer);
        let input = wbvector::memory_store::make_vector_memory_path(&id);

        let args: ToolArgs = serde_json::from_value(json!({
            "input": input,
            "time_field": "t",
            "time_step": "1",
            "resolution": "9",
        }))
        .unwrap();
        let out = EmergingHotSpotAnalysisTool.run(&args, &ctx()).unwrap();
        assert_eq!(out.outputs["time_steps"], json!(12));
        // The growing cluster must read as some kind of hot spot, while the
        // steady far-away background must not (no pattern).
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        let idx = layer.schema.field_index("category").unwrap();
        let cats: Vec<String> = layer
            .iter()
            .map(|f| f.attributes[idx].as_str().unwrap_or("").to_string())
            .collect();
        assert!(
            cats.iter().any(|c| c.ends_with("hot_spot")),
            "expected the growing cluster to read as a hot spot, got {cats:?}"
        );
        assert!(
            cats.iter().any(|c| c == "no_pattern"),
            "expected the steady background to read as no pattern, got {cats:?}"
        );
    }

    #[test]
    fn rejects_missing_required() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            EmergingHotSpotAnalysisTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "x.geojson" })).is_err()); // no time_field
        assert!(bad(json!({ "input": "x.geojson", "time_field": "t", "time_step": "0" })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "time_field": "t" })).is_ok());
    }
}
