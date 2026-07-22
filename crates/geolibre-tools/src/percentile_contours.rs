//! GeoLibre tool: percentile contour polygons of a surface.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Value Percentile Contours* and
//! *Volume Percentile Contours* (Spatial Analyst), unified under one `mode`
//! parameter. Given a continuous raster (a density/KDE surface, suitability,
//! elevation), it emits nested polygons that enclose the most extreme cells —
//! the standard way to turn a kernel-density surface into 50/90/95 % home-range
//! or hotspot footprints.
//!
//! The bundled `contours_from_raster` traces fixed z-levels; nothing computes a
//! percentile-of-value or percentile-of-cumulative-mass threshold. The two
//! modes:
//!
//! - `value` (default): rank cells by their own value; a `percentile` of `p`
//!   selects cells at or above the p-th percentile value — the top `(100 − p)%`
//!   of cells by count.
//! - `volume`: rank cells by value but threshold on **cumulative magnitude**;
//!   `p` selects the highest cells that together carry the top `(100 − p)%` of
//!   the total. Use for "the region accounting for the top X % of the mass".
//!   Assumes non-negative values — pair with `ignore_negative`.
//!
//! Each requested percentile is thresholded into a mask and vectorized with the
//! repo's own `polygonize` (cell-edge rings, holes preserved), so a percentile
//! footprint may be several disjoint polygons. Percentiles are nested and each
//! output feature carries its `percentile` and the `threshold` value; higher
//! percentiles are the innermost, most-extreme rings. Optional `smooth` applies
//! Douglas–Peucker (`geo`'s pure-Rust `Simplify`) to soften the staircase.

use std::collections::{BTreeMap, HashMap};

use geo::{Coord as GeoCoord, LineString, Polygon, Simplify};
use serde_json::{json, Map, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{
    Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring,
};

use crate::common::{band_to_vec, load_input_raster};
use crate::polygonize::{polygonize_to_geojson, PolygonizeParams};
use crate::vector_common::{parse_optional_str, write_or_store_layer};

pub struct PercentileContoursTool;

impl Tool for PercentileContoursTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "percentile_contours",
            display_name: "Percentile Contours",
            summary: "Emit nested polygons enclosing the most extreme cells of a surface, thresholded by value percentile or by cumulative-mass (volume) percentile.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input continuous surface raster (density, suitability, elevation).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output polygon vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "percentiles",
                    description: "Comma-separated percentiles in [0,100] (e.g. \"50,90,95\"). Default \"90\".",
                    required: false,
                },
                ToolParamSpec {
                    name: "mode",
                    description: "'value' (default): cells at/above the p-th percentile value. 'volume': highest cells carrying the top (100-p)% of cumulative value.",
                    required: false,
                },
                ToolParamSpec {
                    name: "ignore_negative",
                    description: "Exclude negative cell values from ranking and output. Default false.",
                    required: false,
                },
                ToolParamSpec {
                    name: "smooth",
                    description: "Simplify the cell-edge polygons with Douglas-Peucker (tolerance ~1.5 cells) to soften the staircase. Default false.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args
            .get("input")
            .and_then(Value::as_str)
            .map(str::trim)
            .unwrap_or("")
            .is_empty()
        {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args
            .get("input")
            .and_then(Value::as_str)
            .filter(|s| !s.trim().is_empty())
            .ok_or_else(|| {
                ToolError::Validation("missing required parameter 'input'".to_string())
            })?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let raster = load_input_raster(input)?;
        if prm.band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "band out of range (raster has {} band(s))",
                raster.bands
            )));
        }
        let rows = raster.rows;
        let cols = raster.cols;
        let nodata = raster.nodata;
        let data = band_to_vec(&raster, prm.band);

        // Collect valid cell values for ranking.
        let mut valid: Vec<f64> = Vec::new();
        for &v in &data {
            if v == nodata || !v.is_finite() {
                continue;
            }
            if prm.ignore_negative && v < 0.0 {
                continue;
            }
            valid.push(v);
        }
        if valid.is_empty() {
            return Err(ToolError::Execution(
                "raster band has no valid cells to contour".to_string(),
            ));
        }
        let mut sorted = valid.clone();
        sorted.sort_by(|a, b| a.total_cmp(b));
        let total: f64 = sorted.iter().sum();

        let cs_x = raster.cell_size_x;
        let cs_y = raster.cell_size_y;
        let cell_area = cs_x * cs_y;
        let y_max = raster.y_min + rows as f64 * cs_y;
        let epsg = raster.crs.epsg;

        let mut out_layer = Layer::new("percentile_contours");
        out_layer.geom_type = Some(GeometryType::Polygon);
        if let Some(e) = epsg {
            out_layer = out_layer.with_crs_epsg(e);
        }
        out_layer.add_field(FieldDef::new("id", FieldType::Integer));
        out_layer.add_field(FieldDef::new("percentile", FieldType::Float));
        out_layer.add_field(FieldDef::new("threshold", FieldType::Float));
        out_layer.add_field(FieldDef::new("mode", FieldType::Text));

        let mut per_pct: Vec<Value> = Vec::new();
        let mut fid = 0u64;
        for &p in &prm.percentiles {
            let threshold = match prm.mode {
                Mode::Value => value_at_percentile(&sorted, p),
                Mode::Volume => volume_threshold(&sorted, total, p),
            };
            // Build the selection mask as a 0/1 label raster.
            let mut labels = vec![0.0f64; rows * cols];
            let mut selected = 0usize;
            for (i, &v) in data.iter().enumerate() {
                if v == nodata || !v.is_finite() {
                    continue;
                }
                if prm.ignore_negative && v < 0.0 {
                    continue;
                }
                if v >= threshold {
                    labels[i] = 1.0;
                    selected += 1;
                }
            }
            ctx.progress.info(&format!(
                "percentile {p}: threshold {threshold:.4}, {selected} cell(s)"
            ));

            let props: HashMap<i64, Map<String, Value>> = HashMap::new();
            let geojson = polygonize_to_geojson(&PolygonizeParams {
                labels: &labels,
                rows,
                cols,
                x_min: raster.x_min,
                y_max,
                cell_size_x: cs_x,
                cell_size_y: cs_y,
                epsg,
                props_by_id: &props,
            });
            let mut polys = parse_polygons(&geojson)?;
            if prm.smooth {
                let tol = 1.5 * cs_x.max(cs_y);
                polys = polys
                    .into_iter()
                    .map(|g| simplify_geometry(g, tol))
                    .collect();
            }
            for g in polys {
                let mut f = Feature::with_geometry(fid, g, out_layer.schema.len());
                f.set_by_index(0, FieldValue::Integer(fid as i64));
                f.set_by_index(1, FieldValue::Float(p));
                f.set_by_index(2, FieldValue::Float(threshold));
                f.set_by_index(3, FieldValue::Text(prm.mode.as_str().to_string()));
                out_layer.push(f);
                fid += 1;
            }

            per_pct.push(json!({
                "percentile": p,
                "threshold": threshold,
                "cell_count": selected,
                "area": selected as f64 * cell_area,
            }));
        }

        let feature_count = out_layer.len();
        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("valid_cells".to_string(), json!(valid.len()));
        outputs.insert("mode".to_string(), json!(prm.mode.as_str()));
        outputs.insert("percentiles".to_string(), json!(per_pct));
        Ok(ToolRunResult { outputs })
    }
}

/// Value at percentile `p` (0..100) of an ascending-sorted slice, by
/// nearest-rank. `p = 0` returns the minimum (selecting all cells), `p = 100`
/// the maximum.
fn value_at_percentile(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    if n == 1 {
        return sorted[0];
    }
    let idx = ((p / 100.0) * (n as f64 - 1.0)).round() as usize;
    sorted[idx.min(n - 1)]
}

/// Cumulative-mass threshold: the smallest value below which the ascending
/// cumulative sum first reaches `p%` of `total`. Cells at or above it carry the
/// top `(100 − p)%` of the mass.
fn volume_threshold(sorted: &[f64], total: f64, p: f64) -> f64 {
    if sorted.len() == 1 || total <= 0.0 {
        return sorted[0];
    }
    let target = (p / 100.0) * total;
    let mut acc = 0.0;
    for &v in sorted {
        acc += v;
        if acc >= target {
            return v;
        }
    }
    *sorted.last().unwrap()
}

// ── GeoJSON polygon parsing + smoothing ─────────────────────────────────────

/// Parses the `polygonize` `FeatureCollection` into `wbvector` polygon
/// geometries (dropping the closing duplicate vertex of each ring).
fn parse_polygons(geojson: &str) -> Result<Vec<Geometry>, ToolError> {
    let v: Value = serde_json::from_str(geojson)
        .map_err(|e| ToolError::Execution(format!("failed parsing polygonize output: {e}")))?;
    let mut out = Vec::new();
    if let Some(features) = v.get("features").and_then(Value::as_array) {
        for f in features {
            let coords = match f.pointer("/geometry/coordinates").and_then(Value::as_array) {
                Some(c) => c,
                None => continue,
            };
            let mut rings = coords.iter().filter_map(ring_from_json);
            let Some(exterior) = rings.next() else {
                continue;
            };
            let interiors: Vec<Ring> = rings.collect();
            out.push(Geometry::Polygon {
                exterior,
                interiors,
            });
        }
    }
    Ok(out)
}

fn ring_from_json(ring: &Value) -> Option<Ring> {
    let pts = ring.as_array()?;
    let mut coords: Vec<Coord> = pts
        .iter()
        .filter_map(|p| {
            let a = p.as_array()?;
            Some(Coord::xy(a.first()?.as_f64()?, a.get(1)?.as_f64()?))
        })
        .collect();
    if coords.len() >= 2 && coords.first() == coords.last() {
        coords.pop();
    }
    if coords.len() < 3 {
        return None;
    }
    Some(Ring::new(coords))
}

/// Douglas–Peucker simplification of a polygon geometry via `geo`.
fn simplify_geometry(geom: Geometry, tol: f64) -> Geometry {
    let Geometry::Polygon {
        exterior,
        interiors,
    } = &geom
    else {
        return geom;
    };
    let poly = Polygon::new(
        ring_to_ls(exterior),
        interiors.iter().map(ring_to_ls).collect(),
    );
    let s = poly.simplify(tol);
    let ext = ls_to_ring(s.exterior());
    let ints: Vec<Ring> = s.interiors().iter().map(ls_to_ring).collect();
    Geometry::Polygon {
        exterior: ext,
        interiors: ints,
    }
}

fn ring_to_ls(r: &Ring) -> LineString {
    LineString::new(
        r.coords()
            .iter()
            .map(|c| GeoCoord { x: c.x, y: c.y })
            .collect(),
    )
}

fn ls_to_ring(ls: &LineString) -> Ring {
    let mut coords: Vec<Coord> = ls.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
    if coords.len() >= 2 && coords.first() == coords.last() {
        coords.pop();
    }
    Ring::new(coords)
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    Value,
    Volume,
}

impl Mode {
    fn as_str(self) -> &'static str {
        match self {
            Mode::Value => "value",
            Mode::Volume => "volume",
        }
    }
}

struct Params {
    percentiles: Vec<f64>,
    mode: Mode,
    ignore_negative: bool,
    smooth: bool,
    band: isize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    // Accept a comma-separated string ("50,90,95"), or a bare number (host UIs
    // and the CLI coerce a single numeric-looking value to a JSON number).
    let raw = match args.get("percentiles") {
        None | Some(Value::Null) => "90".to_string(),
        Some(Value::String(s)) if s.trim().is_empty() => "90".to_string(),
        Some(Value::String(s)) => s.clone(),
        Some(Value::Number(n)) => n.to_string(),
        Some(_) => {
            return Err(ToolError::Validation(
                "parameter 'percentiles' must be a number or comma-separated string".to_string(),
            ))
        }
    };
    let mut percentiles = Vec::new();
    for tok in raw.split(',').map(str::trim).filter(|t| !t.is_empty()) {
        let p: f64 = tok
            .parse()
            .map_err(|_| ToolError::Validation(format!("percentile '{tok}' is not a number")))?;
        if !(0.0..=100.0).contains(&p) {
            return Err(ToolError::Validation(format!(
                "percentile {p} out of range [0, 100]"
            )));
        }
        percentiles.push(p);
    }
    if percentiles.is_empty() {
        return Err(ToolError::Validation(
            "'percentiles' has no valid values".to_string(),
        ));
    }
    percentiles.sort_by(|a, b| a.total_cmp(b));
    percentiles.dedup();
    let mode = match parse_optional_str(args, "mode")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("value") => Mode::Value,
        Some("volume") => Mode::Volume,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown mode '{other}' (expected value or volume)"
            )))
        }
    };
    let ignore_negative = parse_optional_bool(args, "ignore_negative")?.unwrap_or(false);
    let smooth = parse_optional_bool(args, "smooth")?.unwrap_or(false);
    let band = parse_optional_u64(args, "band")?.unwrap_or(1).max(1) as isize - 1;
    Ok(Params {
        percentiles,
        mode,
        ignore_negative,
        smooth,
        band,
    })
}

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
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

fn parse_optional_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be an integer"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be an integer"
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

    fn raster_from(rows: usize, cols: usize, data: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: Default::default(),
            metadata: vec![],
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
        let out = PercentileContoursTool.run(&args, &ctx()).unwrap();
        let layer = crate::vector_common::load_input_layer(out.outputs["output"].as_str().unwrap())
            .unwrap();
        (out, layer)
    }

    fn total_area(layer: &Layer) -> f64 {
        use geo::Area;
        layer
            .features
            .iter()
            .filter_map(|f| f.geometry.as_ref())
            .filter_map(|g| match g {
                Geometry::Polygon {
                    exterior,
                    interiors,
                } => Some(Polygon::new(
                    ring_to_ls(exterior),
                    interiors.iter().map(ring_to_ls).collect(),
                )),
                _ => None,
            })
            .map(|p| p.unsigned_area())
            .sum()
    }

    /// A radial bump: value mode at p=90 selects ~10% of cells; the polygon area
    /// equals that cell count exactly (rectilinear cell boundaries).
    #[test]
    fn value_percentile_area_matches_cell_count() {
        // 10x10 grid with a single tall peak at the centre and a gradient.
        let rows = 10;
        let cols = 10;
        let mut data = vec![0.0; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let dr = r as f64 - 4.5;
                let dc = c as f64 - 4.5;
                data[r * cols + c] = 100.0 - (dr * dr + dc * dc); // peak at centre
            }
        }
        let input = raster_from(rows, cols, &data);
        let (out, layer) = run(json!({ "input": input, "percentiles": "90", "mode": "value" }));
        let pj = &out.outputs["percentiles"][0];
        let cells = pj["cell_count"].as_u64().unwrap();
        // ~10% of 100 cells.
        assert!((5..=15).contains(&cells), "selected {cells} cells");
        // Polygon area equals the selected cell count (cell size 1).
        assert!((total_area(&layer) - cells as f64).abs() < 1e-6);
    }

    /// Nested: a higher percentile encloses a smaller area than a lower one.
    #[test]
    fn higher_percentile_is_smaller() {
        let rows = 20;
        let cols = 20;
        let mut data = vec![0.0; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let dr = r as f64 - 9.5;
                let dc = c as f64 - 9.5;
                data[r * cols + c] = 200.0 - (dr * dr + dc * dc);
            }
        }
        let input = raster_from(rows, cols, &data);
        let (out, _) = run(json!({ "input": input, "percentiles": "50,90", "mode": "value" }));
        let a50 = out.outputs["percentiles"][0]["cell_count"]
            .as_u64()
            .unwrap();
        let a90 = out.outputs["percentiles"][1]["cell_count"]
            .as_u64()
            .unwrap();
        assert!(
            a90 < a50,
            "p90 ({a90}) should select fewer cells than p50 ({a50})"
        );
    }

    /// Volume mode selects fewer high-value cells than value mode when the mass
    /// is concentrated in a few tall cells.
    #[test]
    fn volume_mode_concentrates_on_mass() {
        // Mostly 1s with a few very tall spikes carrying most of the mass.
        let rows = 10;
        let cols = 10;
        let mut data = vec![1.0; rows * cols];
        data[0] = 1000.0;
        data[99] = 1000.0;
        let input = raster_from(rows, cols, &data);
        let (out, _) = run(json!({ "input": input, "percentiles": "90", "mode": "volume" }));
        // The 2 spikes carry ~2000 of ~2098 total -> the top-10%-of-mass line
        // sits just below the spikes, so only a couple of cells are selected.
        let cells = out.outputs["percentiles"][0]["cell_count"]
            .as_u64()
            .unwrap();
        assert!(
            cells <= 5,
            "volume p90 selected {cells} cells (expected the few spikes)"
        );
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = PercentileContoursTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "x.tif", "percentiles": "150" })).is_err());
        assert!(bad(json!({ "input": "x.tif", "mode": "bogus" })).is_err());
        assert!(bad(json!({ "input": "x.tif", "percentiles": "50,90,95" })).is_ok());
    }
}
