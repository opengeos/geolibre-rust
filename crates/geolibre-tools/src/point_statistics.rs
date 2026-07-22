//! GeoLibre tool: rasterize a moving-window statistic of a point attribute.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Point Statistics* (Spatial Analyst).
//! For every output cell a neighborhood (circle or rectangle) is centered on the
//! cell and a statistic is computed over the values of the input point attribute
//! that fall inside that neighborhood. This differs from the bundled density
//! tools (`line_density`, `heat_map`), which measure *density* rather than a
//! statistic of a point *attribute*, and from `neighborhood_summary_statistics`,
//! which is vector→vector rather than a raster.
//!
//! Rather than sweeping every cell against every point, each point is *scattered*
//! onto the cells whose neighborhood contains it (a disc/box of cells around the
//! point), appending its value to that cell's value list. A final pass reduces
//! each cell's value list to the requested statistic. Cells whose neighborhood
//! contains no points are written as no-data.
//!
//! For a geographic (EPSG:4326) input the geometry is projected to a local
//! equirectangular metre frame centered on the extent, so the radius / rectangle
//! size and cell size are all true metres; the output raster's georeferencing is
//! converted back to degrees (with distinct x/y cell sizes) so it still overlays
//! the input. For a projected input everything is in the CRS's native units.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};

use crate::common::{parse_optional_output, write_or_store_output};
use crate::vector_common::load_input_layer;

/// Mean Earth radius (metres) for the local equirectangular projection.
const EARTH_R: f64 = 6_371_000.0;

pub struct PointStatisticsTool;

impl Tool for PointStatisticsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "point_statistics",
            display_name: "Point Statistics",
            summary: "Rasterize a moving-window statistic of a point-feature attribute: for each output cell, a neighborhood (circle or rectangle) collects the values of the chosen numeric field from the points inside it and reduces them to mean/majority/maximum/median/minimum/minority/range/std/sum/variety. Like ArcGIS Point Statistics; complements the density-only line_density/heat_map and the vector-only neighborhood_summary_statistics.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer (Point / MultiPoint).",
                    required: true,
                },
                ToolParamSpec {
                    name: "field",
                    description: "Numeric attribute field to summarize.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output raster path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "statistic",
                    description: "Statistic per neighborhood: mean (default), majority, maximum, median, minimum, minority, range, std, sum, or variety.",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighborhood",
                    description: "Neighborhood shape: 'circle' (default, uses 'radius') or 'rectangle' (a square of side 2*radius).",
                    required: false,
                },
                ToolParamSpec {
                    name: "radius",
                    description: "Neighborhood radius / rectangle half-side (metres for a geographic CRS, CRS units otherwise). Default: shorter extent side / 20.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Output cell size (same units as radius). Default: radius / 5.",
                    required: false,
                },
                ToolParamSpec {
                    name: "epsg",
                    description: "Override the input CRS EPSG code (e.g. when the layer is unlabeled). 4326 triggers the metre projection.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "field")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let field = require_str(args, "field")?;
        let output = parse_optional_output(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let fidx = layer
            .schema
            .field_index(field)
            .ok_or_else(|| ToolError::Validation(format!("field '{field}' not found")))?;

        // Collect (x, y, value) for every point with a finite numeric value.
        let mut raw: Vec<(f64, f64, f64)> = Vec::new();
        let mut skipped = 0usize;
        for feat in &layer.features {
            let Some(geom) = feat.geometry.as_ref() else {
                skipped += 1;
                continue;
            };
            let v = feat.attributes.get(fidx).and_then(|a| a.as_f64());
            let Some(v) = v.filter(|v| v.is_finite()) else {
                skipped += 1;
                continue;
            };
            match geom {
                wbvector::Geometry::Point(c) => raw.push((c.x, c.y, v)),
                wbvector::Geometry::MultiPoint(cs) => {
                    for c in cs {
                        raw.push((c.x, c.y, v));
                    }
                }
                _ => skipped += 1,
            }
        }
        if raw.is_empty() {
            return Err(ToolError::Execution(
                "input contains no point features with a finite value in 'field'".to_string(),
            ));
        }

        // Native bounding box of the points.
        let (mut nxmin, mut nymin, mut nxmax, mut nymax) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for (x, y, _) in &raw {
            nxmin = nxmin.min(*x);
            nxmax = nxmax.max(*x);
            nymin = nymin.min(*y);
            nymax = nymax.max(*y);
        }

        // Local metre projection for a geographic CRS; identity otherwise.
        let epsg = prm.epsg.or_else(|| layer.crs_epsg());
        let geographic = epsg == Some(4326);
        let (kx, ky) = if geographic {
            let lat0 = 0.5 * (nymin + nymax);
            let ky = EARTH_R * std::f64::consts::PI / 180.0;
            let kx = ky * lat0.to_radians().cos().max(1e-9);
            (kx, ky)
        } else {
            (1.0, 1.0)
        };
        let (lon0, lat0) = (nxmin, nymin); // projection origin (native)
        let fwd_x = |x: f64| (x - lon0) * kx;
        let fwd_y = |y: f64| (y - lat0) * ky;

        // Working-frame (metre / native) points.
        let pts: Vec<(f64, f64, f64)> = raw
            .iter()
            .map(|(x, y, v)| (fwd_x(*x), fwd_y(*y), *v))
            .collect();

        // Working-frame extent + defaults for radius / cell size.
        let ext_w = (nxmax - nxmin) * kx;
        let ext_h = (nymax - nymin) * ky;
        let radius = prm
            .radius
            .unwrap_or_else(|| (ext_w.min(ext_h) / 20.0).max(1e-6));
        let cell = prm.cell_size.unwrap_or_else(|| (radius / 5.0).max(1e-6));

        // Grid extent = point bbox padded by the neighborhood reach so every cell
        // whose neighborhood can contain a point exists in the grid.
        let pad = radius + cell;
        let gxmin = -pad;
        let gymin = -pad;
        let gxmax = ext_w + pad;
        let gymax_raw = ext_h + pad;
        let cols = (((gxmax - gxmin) / cell).ceil() as usize).max(1);
        let rows = (((gymax_raw - gymin) / cell).ceil() as usize).max(1);
        let gymax = gymin + rows as f64 * cell; // snap top edge to whole cells

        if rows.saturating_mul(cols) > 40_000_000 {
            return Err(ToolError::Execution(format!(
                "output grid {rows}x{cols} too large; increase cell_size or reduce radius"
            )));
        }

        ctx.progress.info(&format!(
            "{} point(s) -> {rows}x{cols} {} raster (r={radius:.3}, cell={cell:.3}, stat={})",
            pts.len(),
            prm.neighborhood.name(),
            prm.statistic.name()
        ));

        // Scatter each point onto every cell whose neighborhood contains it.
        let mut lists: Vec<Vec<f64>> = vec![Vec::new(); rows * cols];
        let r2 = radius * radius;
        for (pi, &(px, py, v)) in pts.iter().enumerate() {
            let c0 = (((px - radius - gxmin) / cell).floor() as isize).max(0) as usize;
            let c1 = (((px + radius - gxmin) / cell).ceil() as isize).min(cols as isize) as usize;
            let r0 = (((gymax - (py + radius)) / cell).floor() as isize).max(0) as usize;
            let r1 = (((gymax - (py - radius)) / cell).ceil() as isize).min(rows as isize) as usize;
            for r in r0..r1 {
                let cy = gymax - (r as f64 + 0.5) * cell;
                for c in c0..c1 {
                    let cx = gxmin + (c as f64 + 0.5) * cell;
                    let inside = match prm.neighborhood {
                        Neighborhood::Circle => {
                            let dx = cx - px;
                            let dy = cy - py;
                            dx * dx + dy * dy <= r2
                        }
                        Neighborhood::Rectangle => {
                            (cx - px).abs() <= radius && (cy - py).abs() <= radius
                        }
                    };
                    if inside {
                        lists[r * cols + c].push(v);
                    }
                }
            }
            if pi % 256 == 0 {
                ctx.progress.progress((pi as f64 + 1.0) / pts.len() as f64);
            }
        }

        // Reduce each cell's value list to the requested statistic.
        let nodata = -9999.0f64;
        let mut nonempty = 0usize;
        let (mut vmin, mut vmax) = (f64::INFINITY, f64::NEG_INFINITY);
        let data: Vec<f64> = lists
            .iter_mut()
            .map(|vals| {
                if vals.is_empty() {
                    return nodata;
                }
                nonempty += 1;
                let s = prm.statistic.reduce(vals);
                if s.is_finite() {
                    vmin = vmin.min(s);
                    vmax = vmax.max(s);
                }
                s
            })
            .collect();

        // Georeference the output back to native units.
        let (out_cell_x, out_cell_y) = if geographic {
            (cell / kx, cell / ky)
        } else {
            (cell, cell)
        };
        let out_xmin = lon0 + gxmin / kx;
        let out_ymin = lat0 + gymin / ky;

        let crs = CrsInfo {
            epsg,
            wkt: None,
            proj4: None,
        };
        let mut out = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: out_xmin,
            y_min: out_ymin,
            cell_size: out_cell_x,
            cell_size_y: Some(out_cell_y),
            nodata,
            data_type: DataType::F32,
            crs,
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
        outputs.insert("point_count".to_string(), json!(pts.len()));
        outputs.insert("skipped".to_string(), json!(skipped));
        outputs.insert("statistic".to_string(), json!(prm.statistic.name()));
        outputs.insert("neighborhood".to_string(), json!(prm.neighborhood.name()));
        outputs.insert("radius".to_string(), json!(radius));
        outputs.insert("cell_size".to_string(), json!(cell));
        outputs.insert("valued_cells".to_string(), json!(nonempty));
        outputs.insert(
            "value_min".to_string(),
            json!(if vmin.is_finite() { vmin } else { 0.0 }),
        );
        outputs.insert(
            "value_max".to_string(),
            json!(if vmax.is_finite() { vmax } else { 0.0 }),
        );
        Ok(ToolRunResult { outputs })
    }
}

// ── Neighborhood ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Neighborhood {
    Circle,
    Rectangle,
}

impl Neighborhood {
    fn name(self) -> &'static str {
        match self {
            Neighborhood::Circle => "circle",
            Neighborhood::Rectangle => "rectangle",
        }
    }
}

// ── Statistic ─────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Statistic {
    Mean,
    Majority,
    Maximum,
    Median,
    Minimum,
    Minority,
    Range,
    Std,
    Sum,
    Variety,
}

impl Statistic {
    fn name(self) -> &'static str {
        match self {
            Statistic::Mean => "mean",
            Statistic::Majority => "majority",
            Statistic::Maximum => "maximum",
            Statistic::Median => "median",
            Statistic::Minimum => "minimum",
            Statistic::Minority => "minority",
            Statistic::Range => "range",
            Statistic::Std => "std",
            Statistic::Sum => "sum",
            Statistic::Variety => "variety",
        }
    }

    /// Reduces a non-empty value list to this statistic.
    fn reduce(self, vals: &mut [f64]) -> f64 {
        let n = vals.len() as f64;
        match self {
            Statistic::Mean => vals.iter().sum::<f64>() / n,
            Statistic::Sum => vals.iter().sum::<f64>(),
            Statistic::Maximum => vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            Statistic::Minimum => vals.iter().cloned().fold(f64::INFINITY, f64::min),
            Statistic::Range => {
                let mx = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                let mn = vals.iter().cloned().fold(f64::INFINITY, f64::min);
                mx - mn
            }
            Statistic::Std => {
                let mean = vals.iter().sum::<f64>() / n;
                let var = vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
                var.sqrt()
            }
            Statistic::Median => {
                vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                let m = vals.len() / 2;
                if vals.len() % 2 == 1 {
                    vals[m]
                } else {
                    0.5 * (vals[m - 1] + vals[m])
                }
            }
            Statistic::Majority | Statistic::Minority => {
                let (majority, minority) = mode_extremes(vals);
                if self == Statistic::Majority {
                    majority
                } else {
                    minority
                }
            }
            Statistic::Variety => distinct_count(vals) as f64,
        }
    }
}

/// Groups values (by bit pattern) and returns (most-frequent, least-frequent).
/// Ties are broken by the smaller value, matching ArcGIS's deterministic pick.
fn mode_extremes(vals: &mut [f64]) -> (f64, f64) {
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut best_major = (0usize, f64::INFINITY); // (count, value)
    let mut best_minor = (usize::MAX, f64::INFINITY);
    let mut i = 0;
    while i < vals.len() {
        let v = vals[i];
        let mut j = i + 1;
        while j < vals.len() && vals[j] == v {
            j += 1;
        }
        let count = j - i;
        if count > best_major.0 {
            best_major = (count, v);
        }
        if count < best_minor.0 {
            best_minor = (count, v);
        }
        i = j;
    }
    (best_major.1, best_minor.1)
}

/// Counts distinct values (by exact equality on sorted data).
fn distinct_count(vals: &mut [f64]) -> usize {
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut n = 0;
    let mut i = 0;
    while i < vals.len() {
        let v = vals[i];
        let mut j = i + 1;
        while j < vals.len() && vals[j] == v {
            j += 1;
        }
        n += 1;
        i = j;
    }
    n
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    statistic: Statistic,
    neighborhood: Neighborhood,
    radius: Option<f64>,
    cell_size: Option<f64>,
    epsg: Option<u32>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let statistic = match args.get("statistic").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("mean") | Some("MEAN") => Statistic::Mean,
        Some("majority") | Some("MAJORITY") => Statistic::Majority,
        Some("maximum") | Some("MAXIMUM") | Some("max") => Statistic::Maximum,
        Some("median") | Some("MEDIAN") => Statistic::Median,
        Some("minimum") | Some("MINIMUM") | Some("min") => Statistic::Minimum,
        Some("minority") | Some("MINORITY") => Statistic::Minority,
        Some("range") | Some("RANGE") => Statistic::Range,
        Some("std") | Some("STD") | Some("stddev") => Statistic::Std,
        Some("sum") | Some("SUM") => Statistic::Sum,
        Some("variety") | Some("VARIETY") => Statistic::Variety,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'statistic' must be one of mean/majority/maximum/median/minimum/minority/range/std/sum/variety, got '{o}'"
            )))
        }
    };
    let neighborhood = match args
        .get("neighborhood")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some("circle") | Some("CIRCLE") => Neighborhood::Circle,
        Some("rectangle") | Some("RECTANGLE") | Some("square") => Neighborhood::Rectangle,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'neighborhood' must be 'circle' or 'rectangle', got '{o}'"
            )))
        }
    };
    let radius = opt_pos(args, "radius")?;
    let cell_size = opt_pos(args, "cell_size")?;
    let epsg = match opt_f64(args, "epsg")? {
        Some(v) if v > 0.0 => Some(v as u32),
        _ => None,
    };
    Ok(Params {
        statistic,
        neighborhood,
        radius,
        cell_size,
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

fn opt_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
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

fn opt_pos(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match opt_f64(args, key)? {
        Some(v) if v > 0.0 && v.is_finite() => Ok(Some(v)),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a positive number"
        ))),
        None => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Coord, FieldDef, FieldType, Geometry, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a projected point layer from (x, y, value) triples.
    fn point_layer(pts: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (x, y, v) in pts {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(*x, *y))),
                &[("v", (*v).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = PointStatisticsTool.run(&args, &ctx()).unwrap();
        let r = crate::common::load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, r)
    }

    /// Reads the value of the output cell nearest working-frame (px, py), where
    /// the point bbox min is the projection origin. Here the layer is projected
    /// (identity), so native == working coords.
    fn cell_at(r: &Raster, out: &ToolRunResult, px: f64, py: f64) -> f64 {
        let cell = r.cell_size_x; // projected: x == y
        let ymax = r.y_min + r.rows as f64 * cell;
        let col = ((px - r.x_min) / cell).floor() as isize;
        let row = ((ymax - py) / cell).floor() as isize;
        let _ = out;
        r.get(0, row, col)
    }

    /// A cell centered on a lone point returns that point's value for mean/sum,
    /// and empty far-away cells are no-data.
    #[test]
    fn single_point_mean_and_nodata() {
        let (out, r) = run(json!({
            "input": point_layer(&[(0.0, 0.0, 42.0)]),
            "field": "v", "radius": 100.0, "cell_size": 10.0, "statistic": "mean",
        }));
        let v = cell_at(&r, &out, 0.0, 0.0);
        assert!(
            (v - 42.0).abs() < 1e-6,
            "cell on point should be 42, got {v}"
        );
        // A corner cell far outside the neighborhood is no-data.
        assert_eq!(r.get(0, 0, 0), r.nodata);
    }

    /// Two points inside the same neighborhood: mean averages, sum adds,
    /// max/min/range pick extremes, variety counts distinct values.
    #[test]
    fn multi_point_statistics() {
        let layer = point_layer(&[(0.0, 0.0, 10.0), (20.0, 0.0, 30.0)]);
        let common = |stat: &str| {
            run(json!({
                "input": layer.clone(), "field": "v",
                "radius": 100.0, "cell_size": 10.0, "statistic": stat,
            }))
        };
        // Cell centered at the midpoint (10,0) sees both points.
        let mid = |r: &Raster, o: &ToolRunResult| cell_at(r, o, 10.0, 0.0);

        let (o, r) = common("mean");
        assert!((mid(&r, &o) - 20.0).abs() < 1e-6);
        let (o, r) = common("sum");
        assert!((mid(&r, &o) - 40.0).abs() < 1e-6);
        let (o, r) = common("maximum");
        assert!((mid(&r, &o) - 30.0).abs() < 1e-6);
        let (o, r) = common("minimum");
        assert!((mid(&r, &o) - 10.0).abs() < 1e-6);
        let (o, r) = common("range");
        assert!((mid(&r, &o) - 20.0).abs() < 1e-6);
        let (o, r) = common("variety");
        assert!((mid(&r, &o) - 2.0).abs() < 1e-6);
        let (o, r) = common("median");
        assert!((mid(&r, &o) - 20.0).abs() < 1e-6);
    }

    /// Majority / minority pick the most and least frequent values.
    #[test]
    fn majority_minority() {
        // Three coincident-ish points: two 5s and one 9 near the origin.
        let layer = point_layer(&[(0.0, 0.0, 5.0), (1.0, 0.0, 5.0), (-1.0, 0.0, 9.0)]);
        let (o, r) = run(json!({
            "input": layer.clone(), "field": "v",
            "radius": 100.0, "cell_size": 10.0, "statistic": "majority",
        }));
        assert!((cell_at(&r, &o, 0.0, 0.0) - 5.0).abs() < 1e-6);
        let (o, r) = run(json!({
            "input": layer, "field": "v",
            "radius": 100.0, "cell_size": 10.0, "statistic": "minority",
        }));
        assert!((cell_at(&r, &o, 0.0, 0.0) - 9.0).abs() < 1e-6);
    }

    /// std of {10,30} over a full-window cell is 10 (population std).
    #[test]
    fn std_statistic() {
        let (o, r) = run(json!({
            "input": point_layer(&[(0.0, 0.0, 10.0), (20.0, 0.0, 30.0)]),
            "field": "v", "radius": 100.0, "cell_size": 10.0, "statistic": "std",
        }));
        assert!((cell_at(&r, &o, 10.0, 0.0) - 10.0).abs() < 1e-6);
    }

    /// A point just outside the circular radius is excluded, but a rectangle of
    /// the same half-side (its corner reaches farther) includes it.
    #[test]
    fn circle_vs_rectangle() {
        // Point at (14,14): distance from origin ~19.8 > 15 (excluded by circle)
        // but within |dx|,|dy| <= 15 (included by rectangle).
        let layer = point_layer(&[(14.0, 14.0, 7.0)]);
        let (o, r) = run(json!({
            "input": layer.clone(), "field": "v",
            "radius": 15.0, "cell_size": 2.0, "neighborhood": "circle", "statistic": "mean",
        }));
        assert_eq!(cell_at(&r, &o, 0.0, 0.0), r.nodata);
        let (o, r) = run(json!({
            "input": layer, "field": "v",
            "radius": 15.0, "cell_size": 2.0, "neighborhood": "rectangle", "statistic": "mean",
        }));
        assert!((cell_at(&r, &o, 0.0, 0.0) - 7.0).abs() < 1e-6);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            PointStatisticsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err()); // missing field
        assert!(bad(json!({ "input": "a.geojson", "field": "v", "radius": -1.0 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "field": "v", "statistic": "mode" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "field": "v", "neighborhood": "hex" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "field": "v", "statistic": "median" })).is_ok());
    }
}
