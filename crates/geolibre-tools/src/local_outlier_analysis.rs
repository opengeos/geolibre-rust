//! GeoLibre tool: space-time Local Outlier Analysis (Anselin Local Moran's I on
//! an H3 space-time cube).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Local Outlier Analysis* (Space Time
//! Pattern Mining). The bundled `local_morans_i_lisa` is spatial-only; with the
//! H3 space-time cube machinery already built for `emerging_hot_spot_analysis`
//! and `time_series_clustering`, extending LISA into the time dimension answers
//! "when and where did this location behave unlike its neighbours" — which the
//! hot-spot trend categories don't capture.
//!
//! Timestamped points are binned into an H3 cell × time-step cube (as in
//! `emerging_hot_spot_analysis`). Each bin's value is standardised over the whole
//! cube, and the local Moran's I of every bin is computed against its space-time
//! neighbourhood (spatial `kring` × ± `time_window` steps). A seeded permutation
//! test gives a pseudo p-value, and each bin is classified High-High / Low-Low
//! cluster or High-Low / Low-High outlier. The output is one polygon per H3 cell
//! summarising, over time, how many bins were clusters vs. outliers.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use h3o::{CellIndex, LatLng, Resolution};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct LocalOutlierAnalysisTool;

impl Tool for LocalOutlierAnalysisTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "local_outlier_analysis",
            display_name: "Local Outlier Analysis",
            summary: "Space-time Local Outlier Analysis (Anselin Local Moran's I on an H3 space-time cube): classify each cell x time bin as a High-High/Low-Low cluster or High-Low/Low-High outlier with a seeded permutation p-value, then summarise per cell over time — like ArcGIS Local Outlier Analysis. The space-time extension of the bundled spatial-only local_morans_i_lisa.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec { name: "input", description: "Input point layer of timestamped events (lon/lat).", required: true },
                ToolParamSpec { name: "time_field", description: "Field holding each point's time: a number or an ISO-8601 timestamp.", required: true },
                ToolParamSpec { name: "output", description: "Output H3 polygon summary layer. If omitted, stored in memory.", required: false },
                ToolParamSpec { name: "value_field", description: "Numeric field to aggregate per bin (default: point count).", required: false },
                ToolParamSpec { name: "time_step", description: "Width of a time step in the time_field's units (default 1).", required: false },
                ToolParamSpec { name: "resolution", description: "H3 resolution 0..15 (default 7).", required: false },
                ToolParamSpec { name: "kring", description: "Spatial neighbourhood k-ring radius (default 1).", required: false },
                ToolParamSpec { name: "time_window", description: "Temporal neighbourhood half-width in steps (default 1).", required: false },
                ToolParamSpec { name: "permutations", description: "Permutations for the significance test (default 99; 0 disables).", required: false },
                ToolParamSpec { name: "seed", description: "Seed for the permutation RNG (default 1).", required: false },
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

        // Bin observations into (cell, time, value).
        let mut obs: Vec<(CellIndex, f64, f64)> = Vec::new();
        for feature in layer.iter() {
            let Some((lng, lat)) = feature.geometry.as_ref().and_then(point_lnglat) else {
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
            return Err(ToolError::Execution("no usable observations".to_string()));
        }

        let t_min = obs.iter().map(|o| o.1).fold(f64::INFINITY, f64::min);
        let t_max = obs.iter().map(|o| o.1).fold(f64::NEG_INFINITY, f64::max);
        let n_times = (((t_max - t_min) / prm.time_step).floor() as usize) + 1;
        if n_times < 2 {
            return Err(ToolError::Execution(
                "need >= 2 time steps (reduce time_step)".to_string(),
            ));
        }
        let time_bin = |t: f64| (((t - t_min) / prm.time_step).floor() as usize).min(n_times - 1);

        let cells: Vec<CellIndex> = {
            let set: BTreeSet<u64> = obs.iter().map(|o| u64::from(o.0)).collect();
            set.into_iter()
                .map(|r| CellIndex::try_from(r).unwrap())
                .collect()
        };
        let cell_pos: HashMap<CellIndex, usize> =
            cells.iter().enumerate().map(|(i, &c)| (c, i)).collect();
        let n_cells = cells.len();
        let n = n_cells * n_times;
        let mut cube = vec![0.0f64; n];
        for &(cell, time, value) in &obs {
            cube[cell_pos[&cell] * n_times + time_bin(time)] += value;
        }

        // Standardise the cube (z-scores over all bins).
        let mean = cube.iter().sum::<f64>() / n as f64;
        let var = cube.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n as f64;
        let sd = var.sqrt().max(1e-12);
        let z: Vec<f64> = cube.iter().map(|v| (v - mean) / sd).collect();

        // Space-time neighbours per cell (spatial k-ring, indices into `cells`).
        let spatial: Vec<Vec<usize>> = cells
            .iter()
            .map(|c| {
                c.grid_disk::<Vec<CellIndex>>(prm.kring)
                    .into_iter()
                    .filter(|nc| nc != c)
                    .filter_map(|nc| cell_pos.get(&nc).copied())
                    .collect()
            })
            .collect();

        ctx.progress.info(&format!(
            "cube {n_cells} cell(s) x {n_times} step(s); local Moran with {} permutation(s)",
            prm.permutations
        ));

        // Spatio-temporal lag (row-standardised mean of neighbour z) per bin.
        let lag_of = |zz: &[f64], ci: usize, ti: usize| -> f64 {
            let mut sum = 0.0;
            let mut cnt = 0usize;
            for &nc in &spatial[ci] {
                for dt in -(prm.time_window as isize)..=(prm.time_window as isize) {
                    let nt = ti as isize + dt;
                    if nt < 0 || nt >= n_times as isize {
                        continue;
                    }
                    sum += zz[nc * n_times + nt as usize];
                    cnt += 1;
                }
            }
            // Include same-cell other-time neighbours.
            for dt in -(prm.time_window as isize)..=(prm.time_window as isize) {
                if dt == 0 {
                    continue;
                }
                let nt = ti as isize + dt;
                if nt < 0 || nt >= n_times as isize {
                    continue;
                }
                sum += zz[ci * n_times + nt as usize];
                cnt += 1;
            }
            if cnt > 0 {
                sum / cnt as f64
            } else {
                0.0
            }
        };

        // Observed local I per bin.
        let mut local_i = vec![0.0f64; n];
        for ci in 0..n_cells {
            for ti in 0..n_times {
                let idx = ci * n_times + ti;
                local_i[idx] = z[idx] * lag_of(&z, ci, ti);
            }
        }

        // Permutation p-values (conditional-ish: shuffle all z, recompute).
        let mut ge = vec![1usize; n];
        if prm.permutations > 0 {
            let mut perm = z.clone();
            let mut rng = prm.seed;
            for _ in 0..prm.permutations {
                fisher_yates(&mut perm, &mut rng);
                for ci in 0..n_cells {
                    for ti in 0..n_times {
                        let idx = ci * n_times + ti;
                        let pi = perm[idx] * lag_of(&perm, ci, ti);
                        if pi.abs() >= local_i[idx].abs() {
                            ge[idx] += 1;
                        }
                    }
                }
            }
        }
        let denom = (prm.permutations + 1) as f64;

        // Per-cell summary over time.
        let mut out = Layer::new("local_outliers").with_geom_type(GeometryType::Polygon);
        out = out.with_crs_epsg(4326);
        use wbvector::{FieldDef, FieldType};
        out.add_field(FieldDef::new("h3", FieldType::Text));
        out.add_field(FieldDef::new("n_hh", FieldType::Integer));
        out.add_field(FieldDef::new("n_ll", FieldType::Integer));
        out.add_field(FieldDef::new("n_hl", FieldType::Integer));
        out.add_field(FieldDef::new("n_lh", FieldType::Integer));
        out.add_field(FieldDef::new("n_outlier", FieldType::Integer));
        out.add_field(FieldDef::new("dominant", FieldType::Text));
        out.add_field(FieldDef::new("mean_i", FieldType::Float));

        let mut total_outliers = 0i64;
        #[allow(clippy::needless_range_loop)]
        for ci in 0..n_cells {
            let mut counts = [0i64; 4]; // HH, LL, HL, LH
            let mut mean_i = 0.0;
            for ti in 0..n_times {
                let idx = ci * n_times + ti;
                mean_i += local_i[idx];
                let p = ge[idx] as f64 / denom;
                if prm.permutations > 0 && p > 0.05 {
                    continue;
                }
                let zi = z[idx];
                let lag = lag_of(&z, ci, ti);
                let k = match (zi >= 0.0, lag >= 0.0) {
                    (true, true) => 0,   // HH
                    (false, false) => 1, // LL
                    (true, false) => 2,  // HL outlier
                    (false, true) => 3,  // LH outlier
                };
                counts[k] += 1;
            }
            mean_i /= n_times as f64;
            let n_out = counts[2] + counts[3];
            total_outliers += n_out;
            let dominant = dominant_label(&counts);
            let ring = cell_polygon_ring(cells[ci]);
            out.add_feature(
                Some(Geometry::polygon(ring, Vec::new())),
                &[
                    ("h3", cells[ci].to_string().into()),
                    ("n_hh", counts[0].into()),
                    ("n_ll", counts[1].into()),
                    ("n_hl", counts[2].into()),
                    ("n_lh", counts[3].into()),
                    ("n_outlier", n_out.into()),
                    ("dominant", dominant.into()),
                    ("mean_i", mean_i.into()),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("cell feature failed: {e}")))?;
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("cell_count".to_string(), json!(n_cells));
        outputs.insert("time_steps".to_string(), json!(n_times));
        outputs.insert("total_outlier_bins".to_string(), json!(total_outliers));
        Ok(ToolRunResult { outputs })
    }
}

fn dominant_label(counts: &[i64; 4]) -> &'static str {
    let labels = ["HH", "LL", "HL", "LH"];
    let (mut best, mut bi) = (0, usize::MAX);
    for (i, &c) in counts.iter().enumerate() {
        if c > best {
            best = c;
            bi = i;
        }
    }
    if bi == usize::MAX {
        "none"
    } else {
        labels[bi]
    }
}

fn cell_polygon_ring(cell: CellIndex) -> Vec<Coord> {
    cell.boundary()
        .iter()
        .map(|ll| Coord::xy(ll.lng(), ll.lat()))
        .collect()
}

fn point_lnglat(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        _ => None,
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

fn fisher_yates(v: &mut [f64], rng: &mut u64) {
    for i in (1..v.len()).rev() {
        let j = (next_u64(rng) % (i as u64 + 1)) as usize;
        v.swap(i, j);
    }
}

fn next_u64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

struct Params {
    time_field: String,
    value_field: Option<String>,
    time_step: f64,
    resolution: Resolution,
    kring: u32,
    time_window: u32,
    permutations: usize,
    seed: u64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let time_field = require_str(args, "time_field")?.to_string();
    let value_field = parse_optional_str(args, "value_field")?.map(String::from);
    let time_step = opt_f64(args, "time_step")?.unwrap_or(1.0);
    if time_step <= 0.0 {
        return Err(ToolError::Validation("'time_step' must be positive".into()));
    }
    let res_num = opt_u64(args, "resolution")?.unwrap_or(7);
    let resolution = Resolution::try_from(res_num as u8)
        .map_err(|_| ToolError::Validation("'resolution' must be 0..15".into()))?;
    let kring = opt_u64(args, "kring")?.unwrap_or(1) as u32;
    let time_window = opt_u64(args, "time_window")?.unwrap_or(1) as u32;
    let permutations = opt_u64(args, "permutations")?.unwrap_or(99) as usize;
    let seed = opt_u64(args, "seed")?.unwrap_or(1);
    Ok(Params {
        time_field,
        value_field,
        time_step,
        resolution,
        kring,
        time_window,
        permutations,
        seed,
    })
}

fn opt_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("'{key}' must be a number"))),
        _ => Ok(None),
    }
}

fn opt_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("'{key}' must be an integer"))),
        _ => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn pts(rows: &[(f64, f64, f64, f64)]) -> String {
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        l.add_field(FieldDef::new("t", FieldType::Float));
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (lng, lat, t, v) in rows {
            l.add_feature(
                Some(Geometry::point(*lng, *lat)),
                &[("t", (*t).into()), ("v", (*v).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = LocalOutlierAnalysisTool.run(&args, &ctx()).unwrap();
        let l = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, l)
    }

    /// A dense grid of events with one anomalously-high cell yields the expected
    /// cube shape and an outlier/cluster summary per cell.
    #[test]
    fn builds_cube_and_summary() {
        // 3x3 block of H3-adjacent-ish points over 4 time steps, plus a spike.
        let mut rows = Vec::new();
        for gx in 0..4 {
            for gy in 0..4 {
                for t in 0..4 {
                    let lng = -100.0 + gx as f64 * 0.05;
                    let lat = 40.0 + gy as f64 * 0.05;
                    rows.push((lng, lat, t as f64, 1.0));
                }
            }
        }
        // Spike at one cell/time.
        rows.push((-100.0, 40.0, 1.0, 50.0));
        let (out, l) = run(json!({
            "input": pts(&rows), "time_field": "t", "value_field": "v",
            "resolution": 6, "kring": 1, "time_window": 1, "permutations": 49, "seed": 3,
        }));
        assert!(out.outputs["cell_count"].as_i64().unwrap() > 1);
        assert_eq!(out.outputs["time_steps"], json!(4));
        assert!(l.schema.field_index("n_outlier").is_some());
        assert!(l.schema.field_index("dominant").is_some());
        assert_eq!(
            l.features.len() as i64,
            out.outputs["cell_count"].as_i64().unwrap()
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            LocalOutlierAnalysisTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "time_field": "t", "resolution": 99 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "time_field": "t" })).is_ok());
    }
}
