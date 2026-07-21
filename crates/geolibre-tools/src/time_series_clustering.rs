//! GeoLibre tool: cluster space-time-cube locations by temporal profile.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Time Series Clustering* (Space Time
//! Pattern Mining). `emerging_hot_spot_analysis` builds the same H3 space-time
//! cube but classifies each cell independently; nothing groups cells by their
//! *whole temporal trajectory*. This bins timestamped points into an H3 × time
//! cube (identical to `emerging_hot_spot_analysis`), then clusters the per-cell
//! time series with deterministic **k-medoids** (PAM) under one of three
//! distances:
//!
//! * `value` — Euclidean distance on the raw value series;
//! * `profile` — Euclidean distance on each series after z-normalization (mean
//!   0, std 1), so cells with the same *shape* at different levels group together;
//! * `correlation` — `1 − Pearson r`, grouping cells whose series move together.
//!
//! Output is one H3 polygon per cell (renderable straight through PMTiles) with
//! `cluster_id` and an `is_medoid` flag. Deterministic: a seeded splitmix64 RNG
//! chooses the initial medoids, so runs are reproducible in WASM.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use h3o::{CellIndex, LatLng, Resolution};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

const DEFAULT_RESOLUTION: u8 = 7;
const MAX_PAM_ITERS: usize = 50;
/// k-medoids restarts from distinct seeded inits; the lowest-cost result wins,
/// so the clustering is robust to a single unlucky initialization.
const PAM_RESTARTS: usize = 12;

#[derive(Clone, Copy, PartialEq)]
enum Characteristic {
    Value,
    Profile,
    Correlation,
}

pub struct TimeSeriesClusteringTool;

impl Tool for TimeSeriesClusteringTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "time_series_clustering",
            display_name: "Time Series Clustering",
            summary: "Cluster the cells of an H3 space-time cube by the similarity of their time series (value, z-normalized profile, or correlation) with deterministic k-medoids, so places that evolve alike group together, like ArcGIS Time Series Clustering.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point layer in EPSG:4326 (lon/lat) with a timestamp field.",
                    required: true,
                },
                ToolParamSpec {
                    name: "time_field",
                    description: "Timestamp field (numeric or ISO-8601 date/datetime).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output H3 polygon layer with cluster_id/is_medoid. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "num_clusters",
                    description: "Number of clusters k (default 4).",
                    required: false,
                },
                ToolParamSpec {
                    name: "characteristic",
                    description: "Time-series distance: 'value' (raw), 'profile' (z-normalized, default), or 'correlation' (1 − r).",
                    required: false,
                },
                ToolParamSpec {
                    name: "time_step",
                    description: "Width of a time step: a plain number (numeric time_field) or a duration like '1w', '7d', '12h'. Default '1w'.",
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
                    name: "seed",
                    description: "Random seed for reproducible medoid initialization (default 1).",
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

        // ── Bin points into an H3 × time cube (as in emerging_hot_spot_analysis) ──
        ctx.progress.info("binning points into the space-time cube");
        let mut obs: Vec<(CellIndex, f64, f64)> = Vec::new();
        for feature in layer.iter() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let Some((lng, lat)) = representative_lnglat(geom) else {
                continue;
            };
            let Some(time) = feature.attributes.get(time_idx).and_then(parse_time_value) else {
                continue;
            };
            let value = match value_idx {
                Some(vi) => match feature.attributes.get(vi).and_then(FieldValue::as_f64) {
                    Some(v) => v,
                    None => continue,
                },
                None => 1.0,
            };
            if let Ok(ll) = LatLng::new(lat, lng) {
                obs.push((ll.to_cell(prm.resolution), time, value));
            }
        }
        if obs.is_empty() {
            return Err(ToolError::Execution(
                "no usable observations (check time_field / coordinates)".to_string(),
            ));
        }

        let t_min = obs.iter().map(|o| o.1).fold(f64::INFINITY, f64::min);
        let t_max = obs.iter().map(|o| o.1).fold(f64::NEG_INFINITY, f64::max);
        let n_times = (((t_max - t_min) / prm.time_step).floor() as usize) + 1;
        if n_times < 2 {
            return Err(ToolError::Execution(format!(
                "only {n_times} time step(s) span the data; need >= 2 (reduce time_step)"
            )));
        }
        let time_bin = |t: f64| ((t - t_min) / prm.time_step).floor() as usize % n_times.max(1);

        let cells: Vec<CellIndex> = {
            let set: BTreeSet<u64> = obs.iter().map(|o| u64::from(o.0)).collect();
            set.into_iter()
                .map(|r| CellIndex::try_from(r).unwrap())
                .collect()
        };
        let cell_pos: HashMap<CellIndex, usize> =
            cells.iter().enumerate().map(|(i, &c)| (c, i)).collect();
        let n_cells = cells.len();
        if n_cells < prm.num_clusters {
            return Err(ToolError::Execution(format!(
                "only {n_cells} cell(s) but num_clusters={}; reduce k or the resolution",
                prm.num_clusters
            )));
        }
        let mut cube = vec![0.0f64; n_cells * n_times];
        for &(cell, time, value) in &obs {
            let ci = cell_pos[&cell];
            let ti = time_bin(time).min(n_times - 1);
            cube[ci * n_times + ti] += value;
        }

        // ── Per-cell feature vectors for the chosen characteristic ──────────────
        let features: Vec<Vec<f64>> = (0..n_cells)
            .map(|ci| {
                let series = &cube[ci * n_times..ci * n_times + n_times];
                match prm.characteristic {
                    Characteristic::Value => series.to_vec(),
                    Characteristic::Profile | Characteristic::Correlation => znorm(series),
                }
            })
            .collect();

        ctx.progress.info(&format!(
            "clustering {n_cells} cell series into {} cluster(s)",
            prm.num_clusters
        ));

        // Pairwise distance closure.
        let dist = |a: &[f64], b: &[f64]| -> f64 {
            match prm.characteristic {
                Characteristic::Value | Characteristic::Profile => euclidean(a, b),
                Characteristic::Correlation => 1.0 - pearson(a, b),
            }
        };

        let (assign, medoids) = k_medoids(&features, prm.num_clusters, prm.seed, &dist);

        // ── Build output H3 polygons ────────────────────────────────────────────
        let mut out = Layer::new("time_series_clusters")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        out.add_field(FieldDef::new("h3", FieldType::Text));
        out.add_field(FieldDef::new("cluster_id", FieldType::Integer));
        out.add_field(FieldDef::new("is_medoid", FieldType::Integer));

        let medoid_set: BTreeSet<usize> = medoids.iter().copied().collect();
        let mut cluster_sizes = vec![0usize; prm.num_clusters];
        for ci in 0..n_cells {
            cluster_sizes[assign[ci]] += 1;
            let ring = cell_polygon_ring(cells[ci]);
            out.add_feature(
                Some(Geometry::polygon(ring, Vec::new())),
                &[
                    ("h3", cells[ci].to_string().into()),
                    ("cluster_id", (assign[ci] as i64).into()),
                    ("is_medoid", (medoid_set.contains(&ci) as i64).into()),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed writing cell: {e}")))?;
        }

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("cell_count".to_string(), json!(n_cells));
        outputs.insert("time_steps".to_string(), json!(n_times));
        outputs.insert("num_clusters".to_string(), json!(prm.num_clusters));
        outputs.insert("cluster_sizes".to_string(), json!(cluster_sizes));
        Ok(ToolRunResult { outputs })
    }
}

// ── k-medoids (PAM) ────────────────────────────────────────────────────────────

/// Deterministic k-medoids with multiple seeded restarts; the partition with
/// the lowest total point-to-medoid cost wins.
fn k_medoids(
    features: &[Vec<f64>],
    k: usize,
    seed: u64,
    dist: &dyn Fn(&[f64], &[f64]) -> f64,
) -> (Vec<usize>, Vec<usize>) {
    let n = features.len();
    // Precompute the full distance matrix (n is the H3 cell count, modest).
    let mut d = vec![0.0f64; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            let v = dist(&features[i], &features[j]);
            d[i * n + j] = v;
            d[j * n + i] = v;
        }
    }

    let mut state = seed ^ 0x9E37_79B9_7F4A_7C15;
    let mut best: Option<(f64, Vec<usize>, Vec<usize>)> = None;
    for _ in 0..PAM_RESTARTS {
        let (assign, medoids) = pam_once(&d, n, k, &mut state);
        let cost: f64 = (0..n).map(|i| d[i * n + medoids[assign[i]]]).sum();
        if best.as_ref().is_none_or(|(bc, _, _)| cost < *bc) {
            best = Some((cost, assign, medoids));
        }
    }
    let (_, assign, medoids) = best.expect("at least one restart");
    (assign, medoids)
}

/// One PAM run from a fresh seeded init: alternate assign / medoid-update until
/// stable.
fn pam_once(d: &[f64], n: usize, k: usize, state: &mut u64) -> (Vec<usize>, Vec<usize>) {
    let mut medoids: Vec<usize> = Vec::with_capacity(k);
    while medoids.len() < k {
        let cand = (splitmix(state) % n as u64) as usize;
        if !medoids.contains(&cand) {
            medoids.push(cand);
        }
    }

    let mut assign = vec![0usize; n];
    for _ in 0..MAX_PAM_ITERS {
        let mut changed = false;
        for i in 0..n {
            let mut best = 0usize;
            let mut best_d = f64::INFINITY;
            for (m, &med) in medoids.iter().enumerate() {
                let dd = d[i * n + med];
                if dd < best_d {
                    best_d = dd;
                    best = m;
                }
            }
            if assign[i] != best {
                assign[i] = best;
                changed = true;
            }
        }
        let mut new_medoids = medoids.clone();
        for (m, med) in medoids.iter().enumerate() {
            let members: Vec<usize> = (0..n).filter(|&i| assign[i] == m).collect();
            if members.is_empty() {
                continue;
            }
            let mut best = *med;
            let mut best_cost = f64::INFINITY;
            for &cand in &members {
                let cost: f64 = members.iter().map(|&i| d[cand * n + i]).sum();
                if cost < best_cost {
                    best_cost = cost;
                    best = cand;
                }
            }
            new_medoids[m] = best;
        }
        if new_medoids == medoids && !changed {
            break;
        }
        medoids = new_medoids;
    }
    (assign, medoids)
}

fn splitmix(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

// ── Time-series distances ──────────────────────────────────────────────────────

fn euclidean(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
}

fn znorm(s: &[f64]) -> Vec<f64> {
    let n = s.len() as f64;
    let mean = s.iter().sum::<f64>() / n;
    let var = s.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
    let sd = var.sqrt();
    if sd <= 0.0 {
        vec![0.0; s.len()]
    } else {
        s.iter().map(|v| (v - mean) / sd).collect()
    }
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let ma = a.iter().sum::<f64>() / n;
    let mb = b.iter().sum::<f64>() / n;
    let mut cov = 0.0;
    let mut va = 0.0;
    let mut vb = 0.0;
    for (x, y) in a.iter().zip(b) {
        cov += (x - ma) * (y - mb);
        va += (x - ma).powi(2);
        vb += (y - mb).powi(2);
    }
    if va <= 0.0 || vb <= 0.0 {
        0.0
    } else {
        cov / (va.sqrt() * vb.sqrt())
    }
}

// ── H3 / geometry / time helpers (shared shape with emerging_hot_spot) ────────

fn cell_polygon_ring(cell: CellIndex) -> Vec<Coord> {
    cell.boundary()
        .iter()
        .map(|ll| Coord::xy(ll.lng(), ll.lat()))
        .collect()
}

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

fn parse_time_value(fv: &FieldValue) -> Option<f64> {
    if let Some(n) = fv.as_f64() {
        return Some(n);
    }
    fv.as_str().and_then(parse_iso8601_seconds)
}

fn parse_iso8601_seconds(s: &str) -> Option<f64> {
    let s = s.trim();
    if s.len() < 10 {
        return None;
    }
    let b = s.as_bytes();
    let year: i64 = s.get(0..4)?.parse().ok()?;
    if b[4] != b'-' {
        return None;
    }
    let month: i64 = s.get(5..7)?.parse().ok()?;
    let day: i64 = s.get(8..10)?.parse().ok()?;
    let (mut hh, mut mm, mut ss) = (0i64, 0i64, 0i64);
    if s.len() >= 19 && (b[10] == b'T' || b[10] == b' ') {
        hh = s.get(11..13)?.parse().ok()?;
        mm = s.get(14..16)?.parse().ok()?;
        ss = s.get(17..19)?.parse().ok()?;
    }
    Some((days_from_civil(year, month, day) * 86400 + hh * 3600 + mm * 60 + ss) as f64)
}

fn days_from_civil(y: i64, m: i64, d: i64) -> i64 {
    let y = if m <= 2 { y - 1 } else { y };
    let era = if y >= 0 { y } else { y - 399 } / 400;
    let yoe = y - era * 400;
    let doy = (153 * (if m > 2 { m - 3 } else { m + 9 }) + 2) / 5 + d - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146097 + doe - 719468
}

fn parse_time_step(s: &str) -> Result<f64, ToolError> {
    let s = s.trim();
    if let Ok(v) = s.parse::<f64>() {
        if v > 0.0 && v.is_finite() {
            return Ok(v);
        }
        return Err(ToolError::Validation(
            "'time_step' must be positive".to_string(),
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
        "M" => 2_592_000.0,
        "y" => 31_536_000.0,
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
    num_clusters: usize,
    characteristic: Characteristic,
    seed: u64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let time_field = require_str(args, "time_field")?.to_string();
    let time_step = match parse_optional_str(args, "time_step")? {
        Some(s) => parse_time_step(s)?,
        None => 604800.0,
    };
    let value_field = parse_optional_str(args, "value_field")?.map(str::to_string);
    let res_u8 = parse_optional_u32(args, "resolution")?
        .map(|v| v as u8)
        .unwrap_or(DEFAULT_RESOLUTION);
    let resolution = Resolution::try_from(res_u8)
        .map_err(|_| ToolError::Validation(format!("'resolution' must be 0-15, got {res_u8}")))?;
    let num_clusters = parse_optional_u32(args, "num_clusters")?
        .unwrap_or(4)
        .max(1) as usize;
    let characteristic =
        match parse_optional_str(args, "characteristic")?.map(|s| s.trim().to_lowercase()) {
            None => Characteristic::Profile,
            Some(s) if s.is_empty() || s == "profile" => Characteristic::Profile,
            Some(s) if s == "value" => Characteristic::Value,
            Some(s) if s == "correlation" => Characteristic::Correlation,
            Some(other) => {
                return Err(ToolError::Validation(format!(
                    "'characteristic' must be value|profile|correlation, got '{other}'"
                )))
            }
        };
    let seed = parse_optional_u32(args, "seed")?.unwrap_or(1) as u64;
    Ok(Params {
        time_field,
        time_step,
        value_field,
        resolution,
        num_clusters,
        characteristic,
        seed,
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
        Some(_) => Err(ToolError::Validation(format!("'{key}' must be a number"))),
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

    /// Points at (lng, lat, day-of-year seconds, value).
    fn layer_of(pts: &[(f64, f64, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("t", FieldType::Float));
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (lng, lat, t, v) in pts {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(*lng, *lat))),
                &[("t", (*t).into()), ("v", (*v).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = TimeSeriesClusteringTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// Two spatially separated groups with opposite temporal trends split into
    /// two clusters.
    #[test]
    fn separates_opposite_trends() {
        let mut pts = Vec::new();
        // Group A near (0,0): rising over 4 weeks.
        for (wk, v) in [(0.0, 1.0), (1.0, 2.0), (2.0, 3.0), (3.0, 10.0)] {
            for _ in 0..(v as usize) {
                pts.push((0.0, 0.0, wk * 604800.0, 1.0));
            }
        }
        // Group B near (40,40): falling over 4 weeks.
        for (wk, v) in [(0.0, 10.0), (1.0, 3.0), (2.0, 2.0), (3.0, 1.0)] {
            for _ in 0..(v as usize) {
                pts.push((40.0, 40.0, wk * 604800.0, 1.0));
            }
        }
        let input = layer_of(&pts);
        let (out, layer) = run(json!({
            "input": input, "time_field": "t", "num_clusters": 2,
            "characteristic": "profile", "resolution": 3
        }));
        // Two H3 cells, two clusters.
        assert_eq!(out.outputs["num_clusters"], json!(2));
        let gi = layer.schema.field_index("cluster_id").unwrap();
        let clusters: Vec<i64> = layer
            .iter()
            .map(|f| f.attributes[gi].as_i64().unwrap())
            .collect();
        assert_eq!(clusters.len(), 2);
        assert_ne!(
            clusters[0], clusters[1],
            "opposite trends -> different clusters"
        );
    }

    /// Deterministic: same seed reproduces the cluster assignment.
    #[test]
    fn deterministic_by_seed() {
        let mut pts = Vec::new();
        for (lng, lat) in [(0.0, 0.0), (10.0, 10.0), (20.0, 20.0), (30.0, 30.0)] {
            for wk in 0..4 {
                pts.push((lng, lat, wk as f64 * 604800.0, (wk + 1) as f64));
            }
        }
        let input = layer_of(&pts);
        let clusters = |seed: u64| -> Vec<i64> {
            let (_o, l) = run(json!({
                "input": input, "time_field": "t", "value_field": "v",
                "num_clusters": 2, "resolution": 2, "seed": seed
            }));
            let gi = l.schema.field_index("cluster_id").unwrap();
            l.iter()
                .map(|f| f.attributes[gi].as_i64().unwrap())
                .collect()
        };
        assert_eq!(clusters(5), clusters(5), "same seed reproducible");
    }

    /// Exactly one medoid per cluster is flagged.
    #[test]
    fn one_medoid_per_cluster() {
        let mut pts = Vec::new();
        for (lng, lat) in [(0.0, 0.0), (15.0, 15.0), (30.0, 30.0), (45.0, 5.0)] {
            for wk in 0..4 {
                pts.push((lng, lat, wk as f64 * 604800.0, 1.0));
            }
        }
        let input = layer_of(&pts);
        let (out, layer) = run(json!({
            "input": input, "time_field": "t", "num_clusters": 2, "resolution": 2
        }));
        let mi = layer.schema.field_index("is_medoid").unwrap();
        let medoids: i64 = layer
            .iter()
            .map(|f| f.attributes[mi].as_i64().unwrap())
            .sum();
        assert_eq!(medoids, out.outputs["num_clusters"].as_i64().unwrap());
    }

    #[test]
    fn rejects_missing_time_field() {
        let input = layer_of(&[(0.0, 0.0, 0.0, 1.0)]);
        let args: ToolArgs = serde_json::from_value(json!({ "input": input })).unwrap();
        assert!(TimeSeriesClusteringTool.validate(&args).is_err());
    }
}
