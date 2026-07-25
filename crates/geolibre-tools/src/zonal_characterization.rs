//! GeoLibre tool: multi-raster zonal statistics in a single pass.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Zonal Characterization* (Spatial
//! Analyst). The bundled `zonal_statistics` handles exactly one value raster
//! with a fixed statistic, so characterizing zones against a stack of
//! covariates (elevation, slope, NDVI, precipitation) means running it once per
//! raster and joining the results by hand — re-reading and re-indexing the zone
//! grid every time.
//!
//! This does it in one pass with a per-raster choice of statistic, which is the
//! shape needed for building a training or summary table. The shipped
//! `zonal_histogram` is a binned *distribution* per zone, not a statistics
//! matrix, and `zonal_geometry` characterizes zone shape rather than values.
//!
//! Running statistics (mean/sum/min/max/stddev via Welford) need no storage;
//! order statistics (median, percentile, majority, minority, variety) require
//! retaining each zone's values, so those buffers are allocated **only** for the
//! rasters whose statistic demands them.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Layer};

use crate::common::load_input_raster;
use crate::vector_common::{parse_optional_str, write_or_store_layer};

/// A statistic requested for one value raster.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Stat {
    Mean,
    Sum,
    Min,
    Max,
    Range,
    StdDev,
    Median,
    Majority,
    Minority,
    Variety,
    Count,
    Percentile,
}

impl Stat {
    fn parse(s: &str) -> Option<Stat> {
        match s.trim().to_ascii_lowercase().as_str() {
            "mean" => Some(Stat::Mean),
            "sum" => Some(Stat::Sum),
            "min" | "minimum" => Some(Stat::Min),
            "max" | "maximum" => Some(Stat::Max),
            "range" => Some(Stat::Range),
            "stddev" | "std" => Some(Stat::StdDev),
            "median" => Some(Stat::Median),
            "majority" => Some(Stat::Majority),
            "minority" => Some(Stat::Minority),
            "variety" => Some(Stat::Variety),
            "count" => Some(Stat::Count),
            "percentile" => Some(Stat::Percentile),
            _ => None,
        }
    }
    fn label(self) -> &'static str {
        match self {
            Stat::Mean => "mean",
            Stat::Sum => "sum",
            Stat::Min => "min",
            Stat::Max => "max",
            Stat::Range => "range",
            Stat::StdDev => "std",
            Stat::Median => "median",
            Stat::Majority => "majority",
            Stat::Minority => "minority",
            Stat::Variety => "variety",
            Stat::Count => "count",
            Stat::Percentile => "pct",
        }
    }
    /// Does this statistic need every value retained?
    fn needs_values(self) -> bool {
        matches!(
            self,
            Stat::Median | Stat::Majority | Stat::Minority | Stat::Variety | Stat::Percentile
        )
    }
}

/// Running accumulator for the streaming statistics.
#[derive(Clone)]
struct Acc {
    n: usize,
    sum: f64,
    mean: f64,
    m2: f64,
    min: f64,
    max: f64,
    values: Option<Vec<f64>>,
}

impl Acc {
    fn new(keep: bool) -> Acc {
        Acc {
            n: 0,
            sum: 0.0,
            mean: 0.0,
            m2: 0.0,
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
            values: if keep { Some(Vec::new()) } else { None },
        }
    }
    fn push(&mut self, v: f64) {
        self.n += 1;
        self.sum += v;
        let d = v - self.mean;
        self.mean += d / self.n as f64;
        self.m2 += d * (v - self.mean);
        self.min = self.min.min(v);
        self.max = self.max.max(v);
        if let Some(vs) = self.values.as_mut() {
            vs.push(v);
        }
    }
}

pub struct ZonalCharacterizationTool;

impl Tool for ZonalCharacterizationTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "zonal_characterization",
            display_name: "Zonal Characterization",
            summary: "Compute statistics for several value rasters over one zone raster in a single pass, with a per-raster statistic and arbitrary percentiles, like ArcGIS Zonal Characterization.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "zones",
                    description: "Zone raster; each distinct value defines a zone.",
                    required: true,
                },
                ToolParamSpec {
                    name: "rasters",
                    description: "Comma-separated 'path:statistic[:name]' entries, e.g. 'dem.tif:mean,ndvi.tif:median'. Statistics: mean, sum, min, max, range, stddev, median, majority, minority, variety, count, percentile.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output statistics table path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "percentile",
                    description: "Percentile (0-100) applied to any raster whose statistic is 'percentile' (default 50).",
                    required: false,
                },
                ToolParamSpec {
                    name: "zone_band",
                    description: "1-based band to read from the zone raster (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "ignore_nodata",
                    description: "Skip no-data cells in the value rasters (default true). When false, any no-data in a zone makes that zone's statistic null.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "zones")?;
        let specs = parse_rasters(args)?;
        if specs.is_empty() {
            return Err(ToolError::Validation(
                "'rasters' must name at least one raster".to_string(),
            ));
        }
        if let Some(p) = parse_optional_f64(args, "percentile")? {
            if !p.is_finite() || !(0.0..=100.0).contains(&p) {
                return Err(ToolError::Validation(
                    "'percentile' must be between 0 and 100".to_string(),
                ));
            }
        }
        parse_optional_bool(args, "ignore_nodata")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let zones_path = require_str(args, "zones")?;
        let specs = parse_rasters(args)?;
        let output = parse_optional_str(args, "output")?;
        let percentile = parse_optional_f64(args, "percentile")?.unwrap_or(50.0);
        let zone_band = parse_band(args, "zone_band")?;
        let ignore_nodata = parse_optional_bool(args, "ignore_nodata")?.unwrap_or(true);

        let zones = load_input_raster(zones_path)?;
        if zone_band < 0 || zone_band as usize >= zones.bands {
            return Err(ToolError::Validation(format!(
                "zone_band {} out of range ({} band(s))",
                zone_band + 1,
                zones.bands
            )));
        }
        let (rows, cols) = (zones.rows, zones.cols);

        // Load every value raster once and check alignment up front, so a
        // mismatch is an error rather than a silently cropped result.
        let mut value_rasters = Vec::with_capacity(specs.len());
        for spec in &specs {
            let r = load_input_raster(&spec.path)?;
            if r.rows != rows || r.cols != cols {
                return Err(ToolError::Validation(format!(
                    "raster '{}' is {}x{}, expected {rows}x{cols} to match the zone raster",
                    spec.path, r.rows, r.cols
                )));
            }
            value_rasters.push(r);
        }

        ctx.progress.info(&format!(
            "characterizing {rows}x{cols} zones against {} raster(s)",
            specs.len()
        ));

        // zone id -> per-raster accumulators (+ cell count for the zone itself)
        let mut zone_acc: BTreeMap<i64, (Vec<Acc>, usize)> = BTreeMap::new();
        // Zones invalidated by a no-data value cell when ignore_nodata = false.
        let mut poisoned: BTreeMap<i64, Vec<bool>> = BTreeMap::new();

        for r in 0..rows {
            for c in 0..cols {
                let z = zones.get(zone_band, r as isize, c as isize);
                if z == zones.nodata || !z.is_finite() {
                    continue;
                }
                let zid = z as i64;
                let entry = zone_acc.entry(zid).or_insert_with(|| {
                    (
                        specs
                            .iter()
                            .map(|s| Acc::new(s.stat.needs_values()))
                            .collect(),
                        0,
                    )
                });
                entry.1 += 1;
                let pois = poisoned
                    .entry(zid)
                    .or_insert_with(|| vec![false; specs.len()]);
                for (i, vr) in value_rasters.iter().enumerate() {
                    let v = vr.get(0, r as isize, c as isize);
                    if v == vr.nodata || !v.is_finite() {
                        if !ignore_nodata {
                            pois[i] = true;
                        }
                        continue;
                    }
                    entry.0[i].push(v);
                }
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let mut out = Layer::new("zonal_characterization");
        out.add_field(FieldDef::new("zone", FieldType::Integer));
        out.add_field(FieldDef::new("cell_count", FieldType::Integer));
        let labels: Vec<String> = specs
            .iter()
            .map(|s| {
                s.name
                    .clone()
                    .unwrap_or_else(|| format!("{}_{}", s.stem(), s.stat.label()))
            })
            .collect();
        for l in &labels {
            out.add_field(FieldDef::new(l.clone(), FieldType::Float));
        }

        for (zid, (accs, cells)) in &zone_acc {
            let pois = poisoned.get(zid);
            let mut fields: Vec<(String, FieldValue)> = vec![
                ("zone".into(), FieldValue::Integer(*zid)),
                ("cell_count".into(), FieldValue::Integer(*cells as i64)),
            ];
            for (i, spec) in specs.iter().enumerate() {
                let invalid = pois.map(|p| p[i]).unwrap_or(false);
                let v = if invalid {
                    None
                } else {
                    compute(&accs[i], spec.stat, percentile)
                };
                fields.push((
                    labels[i].clone(),
                    v.map(FieldValue::Float).unwrap_or(FieldValue::Null),
                ));
            }
            let refs: Vec<(&str, FieldValue)> = fields
                .iter()
                .map(|(k, v)| (k.as_str(), v.clone()))
                .collect();
            out.add_feature(None, &refs)
                .map_err(|e| ToolError::Execution(format!("failed writing zone row: {e}")))?;
        }

        let zone_count = zone_acc.len();
        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("zone_count".to_string(), json!(zone_count));
        outputs.insert("raster_count".to_string(), json!(specs.len()));
        Ok(ToolRunResult { outputs })
    }
}

struct Spec {
    path: String,
    stat: Stat,
    name: Option<String>,
}

impl Spec {
    /// File stem, used to build a default column name.
    fn stem(&self) -> String {
        std::path::Path::new(&self.path)
            .file_stem()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| "raster".to_string())
            .replace(['-', ' ', '.'], "_")
    }
}

fn compute(acc: &Acc, stat: Stat, percentile: f64) -> Option<f64> {
    if acc.n == 0 {
        return None;
    }
    Some(match stat {
        Stat::Mean => acc.mean,
        Stat::Sum => acc.sum,
        Stat::Min => acc.min,
        Stat::Max => acc.max,
        Stat::Range => acc.max - acc.min,
        Stat::StdDev => (acc.m2 / acc.n as f64).max(0.0).sqrt(),
        Stat::Count => acc.n as f64,
        Stat::Median => quantile(acc.values.as_ref()?, 50.0),
        Stat::Percentile => quantile(acc.values.as_ref()?, percentile),
        Stat::Variety => distinct_count(acc.values.as_ref()?) as f64,
        Stat::Majority => modal(acc.values.as_ref()?, true)?,
        Stat::Minority => modal(acc.values.as_ref()?, false)?,
    })
}

/// Linear-interpolated quantile over a copy of the values.
fn quantile(values: &[f64], pct: f64) -> f64 {
    let mut v = values.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if v.len() == 1 {
        return v[0];
    }
    let pos = (pct / 100.0) * (v.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    if lo == hi {
        v[lo]
    } else {
        v[lo] + (v[hi] - v[lo]) * (pos - lo as f64)
    }
}

/// Most (or least) frequent value. Ties break on the smaller value so the
/// result is deterministic across runs.
fn modal(values: &[f64], most: bool) -> Option<f64> {
    let mut counts: BTreeMap<u64, (f64, usize)> = BTreeMap::new();
    for v in values {
        let key = v.to_bits();
        counts.entry(key).or_insert((*v, 0)).1 += 1;
    }
    let mut best: Option<(f64, usize)> = None;
    for (val, n) in counts.values() {
        let take = match best {
            None => true,
            Some((bv, bn)) => {
                if most {
                    *n > bn || (*n == bn && *val < bv)
                } else {
                    *n < bn || (*n == bn && *val < bv)
                }
            }
        };
        if take {
            best = Some((*val, *n));
        }
    }
    best.map(|(v, _)| v)
}

fn distinct_count(values: &[f64]) -> usize {
    let mut seen: Vec<u64> = values.iter().map(|v| v.to_bits()).collect();
    seen.sort_unstable();
    seen.dedup();
    seen.len()
}

// ── parameter parsing ────────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

/// Parses `path:statistic[:name]` entries. Paths may contain `:` only in the
/// Windows-drive position, so split from the right on the known statistic set.
fn parse_rasters(args: &ToolArgs) -> Result<Vec<Spec>, ToolError> {
    let raw = require_str(args, "rasters")?;
    let mut out = Vec::new();
    for tok in raw.split(',') {
        let tok = tok.trim();
        if tok.is_empty() {
            continue;
        }
        let parts: Vec<&str> = tok.split(':').collect();
        let (path, stat, name) = match parts.len() {
            1 => (parts[0], Stat::Mean, None),
            2 => (
                parts[0],
                Stat::parse(parts[1]).ok_or_else(|| {
                    ToolError::Validation(format!("unknown statistic '{}'", parts[1]))
                })?,
                None,
            ),
            _ => {
                // Last field is the alias, second-to-last the statistic; anything
                // before rejoins as the path (tolerates a drive letter).
                let alias = parts[parts.len() - 1];
                let stat_s = parts[parts.len() - 2];
                let stat = Stat::parse(stat_s).ok_or_else(|| {
                    ToolError::Validation(format!("unknown statistic '{stat_s}'"))
                })?;
                (
                    tok.get(..tok.len() - alias.len() - stat_s.len() - 2)
                        .unwrap_or(parts[0]),
                    stat,
                    Some(alias.to_string()),
                )
            }
        };
        if path.trim().is_empty() {
            return Err(ToolError::Validation(
                "'rasters' entry is missing a path".to_string(),
            ));
        }
        out.push(Spec {
            path: path.trim().to_string(),
            stat,
            name,
        });
    }
    Ok(out)
}

fn parse_band(args: &ToolArgs, key: &str) -> Result<isize, ToolError> {
    let n = match args.get(key) {
        None | Some(Value::Null) => 1u64,
        Some(Value::Number(n)) => n
            .as_u64()
            .ok_or_else(|| ToolError::Validation(format!("'{key}' must be a positive integer")))?,
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
            .map_err(|_| ToolError::Validation(format!("'{key}' must be a positive integer")))?,
        Some(_) => return Err(ToolError::Validation(format!("'{key}' must be an integer"))),
    };
    Ok((n.max(1) - 1) as isize)
}

fn parse_optional_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be a number"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a number"
        ))),
    }
}

fn parse_optional_bool(args: &ToolArgs, k: &str) -> Result<Option<bool>, ToolError> {
    match args.get(k) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{k}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{k}' must be a boolean"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{memory_store, DataType, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    const ND: f64 = -9999.0;

    fn raster_of(rows: usize, cols: usize, data: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: ND,
            data_type: DataType::F32,
            crs: wbraster::CrsInfo::default(),
            metadata: Default::default(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, data[row * cols + col])
                    .unwrap();
            }
        }
        let id = memory_store::put_raster(r);
        memory_store::make_raster_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ZonalCharacterizationTool.run(&args, &ctx()).unwrap();
        let layer = crate::vector_common::load_input_layer(out.outputs["output"].as_str().unwrap())
            .unwrap();
        (out, layer)
    }

    fn cell(layer: &Layer, zone: i64, field: &str) -> Option<f64> {
        let zi = layer.schema.field_index("zone").unwrap();
        let fi = layer.schema.field_index(field).unwrap();
        for f in layer.features.iter() {
            if f.attributes[zi].as_f64() == Some(zone as f64) {
                return f.attributes[fi].as_f64();
            }
        }
        panic!("zone {zone} not in table");
    }

    /// THE point of the tool: several rasters, each with its own statistic, in
    /// a single pass over one zone grid.
    #[test]
    fn multiple_rasters_each_with_own_statistic() {
        // Zones: left column 1, right column 2.
        let zones = raster_of(2, 2, &[1.0, 2.0, 1.0, 2.0]);
        let a = raster_of(2, 2, &[10.0, 100.0, 20.0, 200.0]);
        let b = raster_of(2, 2, &[1.0, 5.0, 3.0, 9.0]);
        let (out, layer) = run(json!({
            "zones": zones,
            "rasters": format!("{a}:mean:a_mean,{b}:max:b_max")
        }));
        assert_eq!(out.outputs["zone_count"], json!(2));
        assert_eq!(out.outputs["raster_count"], json!(2));
        assert_eq!(cell(&layer, 1, "a_mean"), Some(15.0));
        assert_eq!(cell(&layer, 2, "a_mean"), Some(150.0));
        assert_eq!(cell(&layer, 1, "b_max"), Some(3.0));
        assert_eq!(cell(&layer, 2, "b_max"), Some(9.0));
    }

    /// Streaming statistics are exact on a known set.
    #[test]
    fn streaming_statistics_are_correct() {
        let zones = raster_of(1, 4, &[1.0, 1.0, 1.0, 1.0]);
        let v = raster_of(1, 4, &[2.0, 4.0, 4.0, 6.0]);
        let (_, layer) = run(json!({
            "zones": zones,
            "rasters": format!("{v}:mean:m,{v}:sum:s,{v}:stddev:sd,{v}:range:rg,{v}:count:n")
        }));
        assert_eq!(cell(&layer, 1, "m"), Some(4.0));
        assert_eq!(cell(&layer, 1, "s"), Some(16.0));
        // population sd of [2,4,4,6] = sqrt(2) = 1.4142...
        assert!((cell(&layer, 1, "sd").unwrap() - 2.0_f64.sqrt()).abs() < 1e-9);
        assert_eq!(cell(&layer, 1, "rg"), Some(4.0));
        assert_eq!(cell(&layer, 1, "n"), Some(4.0));
    }

    /// Order statistics need the retained-value path.
    #[test]
    fn order_statistics() {
        let zones = raster_of(1, 5, &[1.0; 5]);
        let v = raster_of(1, 5, &[1.0, 2.0, 2.0, 3.0, 100.0]);
        let (_, layer) = run(json!({
            "zones": zones,
            "rasters": format!("{v}:median:med,{v}:majority:maj,{v}:variety:var,{v}:percentile:p")
            , "percentile": 100
        }));
        assert_eq!(cell(&layer, 1, "med"), Some(2.0));
        assert_eq!(cell(&layer, 1, "maj"), Some(2.0));
        assert_eq!(cell(&layer, 1, "var"), Some(4.0));
        assert_eq!(cell(&layer, 1, "p"), Some(100.0));
    }

    /// ignore_nodata=true skips no-data cells; false nulls the whole zone.
    #[test]
    fn nodata_handling() {
        let zones = raster_of(1, 3, &[1.0, 1.0, 1.0]);
        let v = raster_of(1, 3, &[10.0, ND, 20.0]);

        let (_, skip) = run(json!({
            "zones": zones, "rasters": format!("{v}:mean:m"), "ignore_nodata": true
        }));
        assert_eq!(skip.schema.field_index("m").map(|_| ()), Some(()));
        assert_eq!(cell(&skip, 1, "m"), Some(15.0));

        let (_, strict) = run(json!({
            "zones": zones, "rasters": format!("{v}:mean:m"), "ignore_nodata": false
        }));
        assert_eq!(
            cell(&strict, 1, "m"),
            None,
            "a no-data cell must invalidate the zone when ignore_nodata is false"
        );
    }

    /// No-data in the ZONE raster excludes those cells entirely.
    #[test]
    fn zone_nodata_is_excluded() {
        let zones = raster_of(1, 3, &[1.0, ND, 2.0]);
        let v = raster_of(1, 3, &[5.0, 999.0, 7.0]);
        let (out, layer) = run(json!({
            "zones": zones, "rasters": format!("{v}:mean:m")
        }));
        assert_eq!(out.outputs["zone_count"], json!(2));
        assert_eq!(cell(&layer, 1, "m"), Some(5.0));
        assert_eq!(cell(&layer, 2, "m"), Some(7.0));
    }

    /// Mismatched grids are a clear validation error, not a silent crop.
    #[test]
    fn rejects_mismatched_grids() {
        let zones = raster_of(2, 2, &[1.0; 4]);
        let v = raster_of(1, 2, &[1.0, 2.0]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "zones": zones, "rasters": format!("{v}:mean") }))
                .unwrap();
        assert!(matches!(
            ZonalCharacterizationTool.run(&args, &ctx()).unwrap_err(),
            ToolError::Validation(_)
        ));
    }

    #[test]
    fn rejects_bad_parameters() {
        let z = raster_of(1, 1, &[1.0]);
        let v = raster_of(1, 1, &[1.0]);
        let bad = |val: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(val).unwrap();
            ZonalCharacterizationTool.validate(&args).is_err()
        };
        assert!(bad(json!({ "rasters": format!("{v}:mean") })));
        assert!(bad(json!({ "zones": z })));
        assert!(bad(json!({ "zones": z, "rasters": format!("{v}:bogus") })));
        assert!(bad(
            json!({ "zones": z, "rasters": format!("{v}:mean"), "percentile": 150 })
        ));
    }

    /// The quantile and modal helpers themselves.
    #[test]
    fn helper_math() {
        assert_eq!(quantile(&[1.0, 2.0, 3.0], 50.0), 2.0);
        assert_eq!(quantile(&[1.0, 3.0], 50.0), 2.0); // interpolated
        assert_eq!(quantile(&[5.0], 90.0), 5.0);
        assert_eq!(modal(&[1.0, 2.0, 2.0], true), Some(2.0));
        assert_eq!(modal(&[1.0, 2.0, 2.0], false), Some(1.0));
        assert_eq!(distinct_count(&[1.0, 1.0, 2.0]), 2);
    }
}
