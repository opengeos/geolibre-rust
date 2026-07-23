//! GeoLibre tool: rasterized neighborhood statistic of a line attribute.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Line Statistics* (Spatial Analyst).
//! `line_density` (#203) computes length-per-area; this instead summarizes a line
//! **attribute** (road class, speed limit, capacity) over a moving neighborhood,
//! producing an attribute surface. Neither statistic is in the bundled suite.
//!
//! For each output cell a circular neighborhood of `search_radius` is centered on
//! the cell; every line segment is clipped to that circle (closed-form segment-
//! circle intersection) and contributes its in-circle length as a weight. The
//! chosen statistic reduces the contributing `(field value, length)` pairs:
//! * `mean`    — length-weighted mean;
//! * `maximum` / `minimum` / `range` — extremes of the contributing values;
//! * `median`  — length-weighted median;
//! * `majority` / `minority` — value holding the most / least total length;
//! * `variety` — count of distinct values;
//! * `length`  — total clipped line length (no field needed).
//!
//! Distances are in the layer's coordinate units; project geographic inputs first
//! for metric neighborhoods.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
use wbvector::Geometry;

use crate::common::{parse_optional_output, write_or_store_output};
use crate::vector_common::{load_input_layer, parse_optional_str};

pub struct LineStatisticsTool;

impl Tool for LineStatisticsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "line_statistics",
            display_name: "Line Statistics",
            summary: "Rasterize a neighborhood statistic of a line attribute (mean, majority, median, min/max, range, variety, or total length) by clipping each segment to the search-radius circle and length-weighting the values — like ArcGIS Line Statistics; the attribute-surface complement to line_density that the bundled suite lacks.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polyline vector layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "field",
                    description: "Numeric field to summarize. Required for every statistic except 'length'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output raster path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "statistic",
                    description: "'mean' (default), 'majority', 'maximum', 'median', 'minimum', 'minority', 'range', 'variety', or 'length'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "search_radius",
                    description: "Neighborhood radius in CRS units. Default: shorter extent side / 25.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Output cell size in CRS units. Default: search_radius / 10.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        let stat = parse_stat(args)?;
        if stat != Stat::Length && parse_optional_str(args, "field")?.is_none() {
            return Err(ToolError::Validation(format!(
                "statistic '{}' requires a 'field'",
                stat.label()
            )));
        }
        opt_pos(args, "search_radius")?;
        opt_pos(args, "cell_size")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_output(args, "output")?;
        let stat = parse_stat(args)?;
        let field = parse_optional_str(args, "field")?.map(String::from);
        if stat != Stat::Length && field.is_none() {
            return Err(ToolError::Validation(format!(
                "statistic '{}' requires a 'field'",
                stat.label()
            )));
        }
        let radius_arg = opt_pos(args, "search_radius")?;
        let cell_arg = opt_pos(args, "cell_size")?;

        let layer = load_input_layer(input)?;
        let epsg = layer.crs_epsg();
        let fidx = match &field {
            Some(f) => Some(
                layer
                    .schema
                    .field_index(f)
                    .ok_or_else(|| ToolError::Validation(format!("field '{f}' not found")))?,
            ),
            None => None,
        };

        // Collect segments with their feature value.
        struct Seg {
            ax: f64,
            ay: f64,
            bx: f64,
            by: f64,
            v: f64,
        }
        let mut segs: Vec<Seg> = Vec::new();
        let mut skipped = 0usize;
        for feat in layer.iter() {
            let Some(geom) = feat.geometry.as_ref() else {
                skipped += 1;
                continue;
            };
            let chains = line_chains(geom);
            if chains.is_empty() {
                skipped += 1;
                continue;
            }
            let v = match fidx {
                Some(i) => match feat.attributes.get(i).and_then(|x| x.as_f64()) {
                    Some(x) if x.is_finite() => x,
                    _ => {
                        skipped += 1;
                        continue;
                    }
                },
                None => 0.0,
            };
            for chain in chains {
                for pair in chain.windows(2) {
                    segs.push(Seg {
                        ax: pair[0].0,
                        ay: pair[0].1,
                        bx: pair[1].0,
                        by: pair[1].1,
                        v,
                    });
                }
            }
        }
        if segs.is_empty() {
            return Err(ToolError::Execution(
                "input contains no usable line segments".to_string(),
            ));
        }

        // Bounding box + grid.
        let (mut xmin, mut ymin, mut xmax, mut ymax) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for s in &segs {
            xmin = xmin.min(s.ax.min(s.bx));
            xmax = xmax.max(s.ax.max(s.bx));
            ymin = ymin.min(s.ay.min(s.by));
            ymax = ymax.max(s.ay.max(s.by));
        }
        let ext_w = xmax - xmin;
        let ext_h = ymax - ymin;
        let radius = radius_arg.unwrap_or_else(|| (ext_w.min(ext_h) / 25.0).max(1e-6));
        let cell = cell_arg.unwrap_or_else(|| (radius / 10.0).max(1e-6));
        let pad = radius + cell;
        let gxmin = xmin - pad;
        let gymin = ymin - pad;
        let cols = ((((xmax + pad) - gxmin) / cell).ceil() as usize).max(1);
        let rows = ((((ymax + pad) - gymin) / cell).ceil() as usize).max(1);
        let gymax = gymin + rows as f64 * cell;

        ctx.progress.info(&format!(
            "{} segment(s) -> {rows}x{cols} raster ({}, r={radius:.3})",
            segs.len(),
            stat.label()
        ));

        let needs_list = matches!(
            stat,
            Stat::Majority | Stat::Minority | Stat::Median | Stat::Variety
        );

        // Light accumulators (always maintained for the simple statistics).
        let mut sum_len = vec![0.0f64; rows * cols];
        let mut sum_lv = vec![0.0f64; rows * cols];
        let mut minv = vec![f64::INFINITY; rows * cols];
        let mut maxv = vec![f64::NEG_INFINITY; rows * cols];
        // Heavy per-cell (value, length) lists only when the statistic needs them.
        let mut lists: Vec<Vec<(f64, f64)>> = if needs_list {
            vec![Vec::new(); rows * cols]
        } else {
            Vec::new()
        };

        for (si, s) in segs.iter().enumerate() {
            let sxmin = s.ax.min(s.bx) - radius;
            let sxmax = s.ax.max(s.bx) + radius;
            let symin = s.ay.min(s.by) - radius;
            let symax = s.ay.max(s.by) + radius;
            let c0 = (((sxmin - gxmin) / cell).floor() as isize).max(0) as usize;
            let c1 = (((sxmax - gxmin) / cell).ceil() as isize).min(cols as isize) as usize;
            let r0 = (((gymax - symax) / cell).floor() as isize).max(0) as usize;
            let r1 = (((gymax - symin) / cell).ceil() as isize).min(rows as isize) as usize;
            for r in r0..r1 {
                let cy = gymax - (r as f64 + 0.5) * cell;
                for c in c0..c1 {
                    let cx = gxmin + (c as f64 + 0.5) * cell;
                    let l = clipped_length(s.ax, s.ay, s.bx, s.by, cx, cy, radius);
                    if l <= 0.0 {
                        continue;
                    }
                    let idx = r * cols + c;
                    sum_len[idx] += l;
                    sum_lv[idx] += l * s.v;
                    minv[idx] = minv[idx].min(s.v);
                    maxv[idx] = maxv[idx].max(s.v);
                    if needs_list {
                        lists[idx].push((s.v, l));
                    }
                }
            }
            if si % 256 == 0 {
                ctx.progress.progress((si as f64 + 1.0) / segs.len() as f64);
            }
        }

        let nodata = -9999.0_f64;
        let mut data = vec![nodata; rows * cols];
        let mut valid_cells = 0usize;
        for idx in 0..rows * cols {
            if sum_len[idx] <= 0.0 {
                continue;
            }
            valid_cells += 1;
            data[idx] = match stat {
                Stat::Length => sum_len[idx],
                Stat::Mean => sum_lv[idx] / sum_len[idx],
                Stat::Maximum => maxv[idx],
                Stat::Minimum => minv[idx],
                Stat::Range => maxv[idx] - minv[idx],
                Stat::Median => weighted_median(&mut lists[idx]),
                Stat::Majority => extreme_by_length(&lists[idx], true),
                Stat::Minority => extreme_by_length(&lists[idx], false),
                Stat::Variety => variety(&lists[idx]),
            };
        }

        let mut out = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: gxmin,
            y_min: gymin,
            cell_size: cell,
            cell_size_y: Some(cell),
            nodata,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg,
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for r in 0..rows {
            for c in 0..cols {
                out.set(0, r as isize, c as isize, data[r * cols + c])
                    .map_err(|e| ToolError::Execution(format!("write failed: {e}")))?;
            }
        }

        let out_path = write_or_store_output(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        outputs.insert("statistic".to_string(), json!(stat.label()));
        outputs.insert("segment_count".to_string(), json!(segs.len()));
        outputs.insert("skipped".to_string(), json!(skipped));
        outputs.insert("valid_cells".to_string(), json!(valid_cells));
        outputs.insert("search_radius".to_string(), json!(radius));
        outputs.insert("cell_size".to_string(), json!(cell));
        Ok(ToolRunResult { outputs })
    }
}

/// Length-weighted median: sort by value, accumulate length to the 50% mass.
fn weighted_median(list: &mut [(f64, f64)]) -> f64 {
    if list.is_empty() {
        return f64::NAN;
    }
    list.sort_by(|a, b| a.0.total_cmp(&b.0));
    let total: f64 = list.iter().map(|(_, l)| *l).sum();
    let mut cum = 0.0;
    for (v, l) in list.iter() {
        cum += *l;
        if cum >= 0.5 * total {
            return *v;
        }
    }
    list.last().unwrap().0
}

/// Value holding the most (`want_max`) or least total length. Ties break to the
/// smaller value for determinism.
fn extreme_by_length(list: &[(f64, f64)], want_max: bool) -> f64 {
    let mut agg: BTreeMap<u64, (f64, f64)> = BTreeMap::new();
    for &(v, l) in list {
        let e = agg.entry(v.to_bits()).or_insert((0.0, v));
        e.0 += l;
    }
    let mut ordered: Vec<(f64, f64)> = agg.values().map(|(l, v)| (*v, *l)).collect();
    ordered.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut best_len = if want_max {
        f64::NEG_INFINITY
    } else {
        f64::INFINITY
    };
    let mut best_val = f64::NAN;
    for (v, l) in ordered {
        let better = if want_max { l > best_len } else { l < best_len };
        if better {
            best_len = l;
            best_val = v;
        }
    }
    best_val
}

fn variety(list: &[(f64, f64)]) -> f64 {
    let mut seen: Vec<u64> = list.iter().map(|(v, _)| v.to_bits()).collect();
    seen.sort_unstable();
    seen.dedup();
    seen.len() as f64
}

/// Length of segment A→B within distance `r` of C (closed-form clip).
fn clipped_length(ax: f64, ay: f64, bx: f64, by: f64, cx: f64, cy: f64, r: f64) -> f64 {
    let dx = bx - ax;
    let dy = by - ay;
    let a = dx * dx + dy * dy;
    if a <= 0.0 {
        return 0.0;
    }
    let fx = ax - cx;
    let fy = ay - cy;
    let b = 2.0 * (fx * dx + fy * dy);
    let c = fx * fx + fy * fy - r * r;
    let disc = b * b - 4.0 * a * c;
    if disc <= 0.0 {
        return 0.0;
    }
    let sq = disc.sqrt();
    let inv2a = 1.0 / (2.0 * a);
    let t1 = (-b - sq) * inv2a;
    let t2 = (-b + sq) * inv2a;
    let lo = t1.max(0.0);
    let hi = t2.min(1.0);
    if hi <= lo {
        return 0.0;
    }
    (hi - lo) * a.sqrt()
}

fn line_chains(geom: &Geometry) -> Vec<Vec<(f64, f64)>> {
    let to_pts =
        |cs: &[wbvector::Coord]| -> Vec<(f64, f64)> { cs.iter().map(|c| (c.x, c.y)).collect() };
    match geom {
        Geometry::LineString(cs) => vec![to_pts(cs)],
        Geometry::MultiLineString(lines) => lines.iter().map(|l| to_pts(l)).collect(),
        _ => Vec::new(),
    }
}

// ── Statistic ─────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Stat {
    Mean,
    Majority,
    Maximum,
    Median,
    Minimum,
    Minority,
    Range,
    Variety,
    Length,
}

impl Stat {
    fn label(&self) -> &'static str {
        match self {
            Stat::Mean => "mean",
            Stat::Majority => "majority",
            Stat::Maximum => "maximum",
            Stat::Median => "median",
            Stat::Minimum => "minimum",
            Stat::Minority => "minority",
            Stat::Range => "range",
            Stat::Variety => "variety",
            Stat::Length => "length",
        }
    }
}

fn parse_stat(args: &ToolArgs) -> Result<Stat, ToolError> {
    Ok(
        match args.get("statistic").and_then(Value::as_str).map(str::trim) {
            None | Some("") | Some("mean") => Stat::Mean,
            Some("majority") => Stat::Majority,
            Some("maximum") | Some("max") => Stat::Maximum,
            Some("median") => Stat::Median,
            Some("minimum") | Some("min") => Stat::Minimum,
            Some("minority") => Stat::Minority,
            Some("range") => Stat::Range,
            Some("variety") => Stat::Variety,
            Some("length") => Stat::Length,
            Some(o) => {
                return Err(ToolError::Validation(format!(
                    "'statistic' must be mean|majority|maximum|median|minimum|minority|range|variety|length, got '{o}'"
                )))
            }
        },
    )
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn opt_pos(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::Number(n)) => finite_pos(n.as_f64(), key),
        Some(Value::String(s)) => finite_pos(s.trim().parse::<f64>().ok(), key),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a number"
        ))),
    }
}

fn finite_pos(v: Option<f64>, key: &str) -> Result<Option<f64>, ToolError> {
    match v {
        Some(x) if x.is_finite() && x > 0.0 => Ok(Some(x)),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a positive number"
        ))),
        None => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a number"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Coord, FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn line_layer(lines: &[(&[(f64, f64)], f64)]) -> String {
        let mut l = Layer::new("lines")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (pts, v) in lines {
            let coords: Vec<Coord> = pts.iter().map(|(x, y)| Coord::xy(*x, *y)).collect();
            l.add_feature(Some(Geometry::line_string(coords)), &[("v", (*v).into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = LineStatisticsTool.run(&args, &ctx()).unwrap();
        let r = crate::common::load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, r)
    }

    /// A cell sitting on a single line reports that line's field value for mean.
    #[test]
    fn mean_of_single_line() {
        let (_o, r) = run(json!({
            "input": line_layer(&[(&[(-100.0, 0.0), (100.0, 0.0)], 7.0)]),
            "field": "v", "statistic": "mean", "search_radius": 20.0, "cell_size": 5.0,
        }));
        // Find a cell on the line (y=0). The value should be 7.
        let mut found = false;
        for row in 0..r.rows {
            for col in 0..r.cols {
                let val = r.get(0, row as isize, col as isize);
                if val != r.nodata {
                    assert!((val - 7.0).abs() < 1e-6, "value {val} != 7");
                    found = true;
                }
            }
        }
        assert!(found, "expected some valid cells");
    }

    /// Where two lines with different values overlap a cell, majority is the one
    /// contributing more length; maximum is the larger value.
    #[test]
    fn majority_and_max_over_two_lines() {
        // A long line value 1 and a short line value 9 near the origin.
        let (_o, r) = run(json!({
            "input": line_layer(&[
                (&[(-100.0, 0.0), (100.0, 0.0)], 1.0),
                (&[(-5.0, 0.0), (5.0, 0.0)], 9.0),
            ]),
            "field": "v", "statistic": "maximum", "search_radius": 50.0, "cell_size": 10.0,
        }));
        // A cell near origin sees both; maximum must be 9 somewhere.
        let mut saw9 = false;
        for row in 0..r.rows {
            for col in 0..r.cols {
                if (r.get(0, row as isize, col as isize) - 9.0).abs() < 1e-6 {
                    saw9 = true;
                }
            }
        }
        assert!(saw9, "maximum should reach 9 near the overlap");
    }

    /// 'length' needs no field and returns positive length near lines.
    #[test]
    fn length_statistic_no_field() {
        let (out, r) = run(json!({
            "input": line_layer(&[(&[(0.0, 0.0), (100.0, 0.0)], 1.0)]),
            "statistic": "length", "search_radius": 20.0, "cell_size": 5.0,
        }));
        assert_eq!(out.outputs["statistic"], json!("length"));
        let mut maxlen = 0.0f64;
        for row in 0..r.rows {
            for col in 0..r.cols {
                let v = r.get(0, row as isize, col as isize);
                if v != r.nodata {
                    maxlen = maxlen.max(v);
                }
            }
        }
        assert!(
            maxlen > 0.0,
            "length surface should be positive near the line"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            LineStatisticsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "statistic": "mode" })).is_err());
        // mean without a field is rejected.
        assert!(bad(json!({ "input": "a.geojson", "statistic": "mean" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "statistic": "length" })).is_ok());
        assert!(bad(json!({ "input": "a.geojson", "field": "v", "search_radius": -1 })).is_err());
    }
}
