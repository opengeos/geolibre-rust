//! GeoLibre tool: point-to-point line of sight over a DEM.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Line Of Sight* / *Construct Sight
//! Lines* (3D Analyst). The bundled suite answers "what area can this observer
//! see" (`viewshed`, `visibility_index`) and "what is the horizon profile"
//! (`skyline_analysis`), but not the vector question **"can A see B, and where
//! exactly does the line break"** — a different output shape: classified line
//! segments, not a raster.
//!
//! For each observer→target pair the sight line is sampled across the DEM (one
//! sample per cell along the segment). Walking outward from the observer we keep
//! the running maximum *vertical angle* `(ground − observer_z) / distance`; a
//! ground sample is **visible** when its own angle is at least the running
//! maximum (nothing closer rises above the line of sight to it) and
//! **obstructed** otherwise. The target — lifted by its offset — is visible when
//! its angle clears the running maximum of the intervening terrain. The first
//! obstructed sample is the obstruction point.
//!
//! Output is a line layer: each pair contributes one LineString per contiguous
//! visible/obstructed run (so it renders green/red directly), carrying the
//! observer and target ids, the run's `visible` flag, and the pair's overall
//! `tgt_vis` (target visible) flag. Observer and target heights above ground are
//! `observer_offset` / `target_offset`. Pairs are formed all-to-all, or matched
//! by a shared `pair_field`. DEM and points must share a projected CRS (distances
//! and heights in the same units); a geographic DEM is rejected.

use std::collections::{BTreeMap, HashMap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::common::load_input_raster;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Safety cap on all-to-all pair expansion (observers × targets).
const DEFAULT_MAX_PAIRS: usize = 100_000;
/// Angular slack so a sample exactly on the sight line reads as visible.
const ANGLE_EPS: f64 = 1e-9;

pub struct LineOfSightTool;

impl Tool for LineOfSightTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "line_of_sight",
            display_name: "Line Of Sight",
            summary: "Point-to-point visibility over a DEM: for each observer-target pair emit the sight line split into visible and obstructed segments plus a target-visible flag, like ArcGIS Line Of Sight.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "dem",
                    description: "Input elevation raster (projected CRS; distances/heights in its units).",
                    required: true,
                },
                ToolParamSpec {
                    name: "observers",
                    description: "Observer point vector layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "targets",
                    description: "Target point vector layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output line vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "observer_offset",
                    description: "Observer height above ground, in DEM units (default 1.75).",
                    required: false,
                },
                ToolParamSpec {
                    name: "target_offset",
                    description: "Target height above ground, in DEM units (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "pair_field",
                    description: "Optional field present on both layers; observers and targets are paired where this field matches (instead of all-to-all).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "DEM band to sample (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "dem")?;
        require_str(args, "observers")?;
        require_str(args, "targets")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let dem_path = require_str(args, "dem")?;
        let obs_path = require_str(args, "observers")?;
        let tgt_path = require_str(args, "targets")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        ctx.progress.info("reading DEM and points");
        let dem = load_input_raster(dem_path)?;
        // A geographic DEM would mix degrees (distance) and metres (height).
        if dem.crs.epsg == Some(4326) {
            return Err(ToolError::Validation(
                "DEM is geographic (EPSG:4326); reproject to a projected CRS so distances and heights share units".to_string(),
            ));
        }
        let obs_layer = load_input_layer(obs_path)?;
        let tgt_layer = load_input_layer(tgt_path)?;

        let observers = collect_points(&obs_layer, prm.pair_field.as_deref());
        let targets = collect_points(&tgt_layer, prm.pair_field.as_deref());
        if observers.is_empty() || targets.is_empty() {
            return Err(ToolError::Execution(
                "no usable observer or target points".to_string(),
            ));
        }

        // ── Build the observer→target pair list ───────────────────────────────
        let pairs: Vec<(usize, usize)> = match &prm.pair_field {
            Some(_) => {
                let mut by_key: HashMap<String, Vec<usize>> = HashMap::new();
                for (ti, t) in targets.iter().enumerate() {
                    by_key.entry(t.key.clone()).or_default().push(ti);
                }
                let mut v = Vec::new();
                for (oi, o) in observers.iter().enumerate() {
                    if let Some(tis) = by_key.get(&o.key) {
                        for &ti in tis {
                            v.push((oi, ti));
                        }
                    }
                }
                v
            }
            None => {
                let mut v = Vec::with_capacity(observers.len() * targets.len());
                for oi in 0..observers.len() {
                    for ti in 0..targets.len() {
                        v.push((oi, ti));
                    }
                }
                v
            }
        };
        if pairs.len() > DEFAULT_MAX_PAIRS {
            return Err(ToolError::Execution(format!(
                "{} observer-target pairs exceed the {DEFAULT_MAX_PAIRS} cap; use a pair_field or fewer points",
                pairs.len()
            )));
        }

        ctx.progress
            .info(&format!("computing {} sight line(s)", pairs.len()));

        let mut out = Layer::new("line_of_sight")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(dem.crs.epsg.unwrap_or(0));
        out.add_field(FieldDef::new("obs_id", FieldType::Text));
        out.add_field(FieldDef::new("tgt_id", FieldType::Text));
        out.add_field(FieldDef::new("visible", FieldType::Integer));
        out.add_field(FieldDef::new("tgt_vis", FieldType::Integer));

        let mut visible_pairs = 0usize;
        let mut skipped = 0usize;
        for &(oi, ti) in &pairs {
            let o = &observers[oi];
            let t = &targets[ti];
            let Some(los) =
                compute_los(&dem, prm.band, o, t, prm.observer_offset, prm.target_offset)
            else {
                skipped += 1;
                continue;
            };
            if los.target_visible {
                visible_pairs += 1;
            }
            for run in &los.runs {
                if run.pts.len() < 2 {
                    continue;
                }
                let coords: Vec<Coord> = run.pts.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
                out.add_feature(
                    Some(Geometry::line_string(coords)),
                    &[
                        ("obs_id", o.key.clone().into()),
                        ("tgt_id", t.key.clone().into()),
                        ("visible", (run.visible as i64).into()),
                        ("tgt_vis", (los.target_visible as i64).into()),
                    ],
                )
                .map_err(|e| ToolError::Execution(format!("failed building sight line: {e}")))?;
            }
        }

        let feature_count = out.len();
        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("pair_count".to_string(), json!(pairs.len()));
        outputs.insert("visible_pairs".to_string(), json!(visible_pairs));
        outputs.insert("skipped_pairs".to_string(), json!(skipped));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        Ok(ToolRunResult { outputs })
    }
}

// ── Line-of-sight computation ────────────────────────────────────────────────

struct PointObs {
    x: f64,
    y: f64,
    key: String,
}

/// A contiguous run of samples with the same visibility, as a polyline.
struct Run {
    visible: bool,
    pts: Vec<(f64, f64)>,
}

struct Los {
    target_visible: bool,
    runs: Vec<Run>,
}

/// Computes the line of sight from `o` to `t` over the DEM. Returns `None` when
/// either endpoint has no DEM value or the two points coincide.
fn compute_los(
    dem: &Raster,
    band: isize,
    o: &PointObs,
    t: &PointObs,
    obs_off: f64,
    tgt_off: f64,
) -> Option<Los> {
    let obs_ground = sample_dem(dem, band, o.x, o.y)?;
    let tgt_ground = sample_dem(dem, band, t.x, t.y)?;
    let obs_z = obs_ground + obs_off;
    let tgt_z = tgt_ground + tgt_off;
    let dx = t.x - o.x;
    let dy = t.y - o.y;
    let total = dx.hypot(dy);
    if total <= 0.0 {
        return None;
    }
    // One sample per cell along the segment.
    let cell = dem.cell_size_x.min(dem.cell_size_y).max(f64::MIN_POSITIVE);
    let steps = ((total / cell).ceil() as usize).max(2);

    let mut max_angle = f64::NEG_INFINITY;
    let mut runs: Vec<Run> = Vec::new();
    // Start the first run at the observer, visible by definition.
    let mut cur = Run {
        visible: true,
        pts: vec![(o.x, o.y)],
    };
    // Intermediate ground samples (exclude the observer at i=0 and target at
    // i=steps; the target is handled with its offset below).
    for i in 1..steps {
        let f = i as f64 / steps as f64;
        let x = o.x + dx * f;
        let y = o.y + dy * f;
        let d = total * f;
        let (px, py) = (x, y);
        let vis = match sample_dem(dem, band, x, y) {
            Some(g) => {
                let angle = (g - obs_z) / d;
                let v = angle >= max_angle - ANGLE_EPS;
                if angle > max_angle {
                    max_angle = angle;
                }
                v
            }
            // No terrain data here: treat as non-blocking (visible), don't raise
            // the horizon.
            None => true,
        };
        if vis == cur.visible {
            cur.pts.push((px, py));
        } else {
            // Close the current run at this breakpoint and start a new one that
            // shares the vertex, so the segments join without a gap.
            cur.pts.push((px, py));
            let finished = std::mem::replace(
                &mut cur,
                Run {
                    visible: vis,
                    pts: vec![(px, py)],
                },
            );
            runs.push(finished);
        }
    }
    // The target sample, lifted by its offset. The final segment (last sample →
    // target) takes the target's visibility.
    let target_angle = (tgt_z - obs_z) / total;
    let target_visible = target_angle >= max_angle - ANGLE_EPS;
    let last_pt = *cur.pts.last().unwrap();
    if target_visible == cur.visible {
        cur.pts.push((t.x, t.y));
        runs.push(cur);
    } else {
        runs.push(cur);
        runs.push(Run {
            visible: target_visible,
            pts: vec![last_pt, (t.x, t.y)],
        });
    }

    Some(Los {
        target_visible,
        runs,
    })
}

/// Samples the DEM at world `(x, y)`, returning `None` for outside-extent or
/// no-data cells.
fn sample_dem(dem: &Raster, band: isize, x: f64, y: f64) -> Option<f64> {
    let (col, row) = dem.world_to_pixel(x, y)?;
    let v = dem.get(band, row, col);
    if v == dem.nodata || v.is_nan() {
        None
    } else {
        Some(v)
    }
}

// ── Point collection ─────────────────────────────────────────────────────────

/// Collects representative points and a pairing key (the `pair_field` value, or
/// the feature index) from a layer.
fn collect_points(layer: &Layer, pair_field: Option<&str>) -> Vec<PointObs> {
    let field_idx = pair_field.and_then(|f| layer.schema.field_index(f));
    let mut out = Vec::new();
    for (i, feature) in layer.features.iter().enumerate() {
        let Some(geom) = feature.geometry.as_ref() else {
            continue;
        };
        let Some((x, y)) = representative_xy(geom) else {
            continue;
        };
        let key = match field_idx {
            Some(fi) => field_key(&feature.attributes[fi]),
            None => i.to_string(),
        };
        out.push(PointObs { x, y, key });
    }
    out
}

/// A stable string key for a field value (so numeric and text ids both pair).
fn field_key(fv: &FieldValue) -> String {
    if let Some(i) = fv.as_i64() {
        i.to_string()
    } else if let Some(f) = fv.as_f64() {
        format!("{f}")
    } else {
        fv.as_str().unwrap_or("").to_string()
    }
}

fn representative_xy(geom: &Geometry) -> Option<(f64, f64)> {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0u64;
    accumulate(geom, &mut sx, &mut sy, &mut n);
    (n > 0).then(|| (sx / n as f64, sy / n as f64))
}

fn accumulate(geom: &Geometry, sx: &mut f64, sy: &mut f64, n: &mut u64) {
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
                accumulate(g, sx, sy, n);
            }
        }
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    observer_offset: f64,
    target_offset: f64,
    pair_field: Option<String>,
    band: isize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let observer_offset = parse_optional_f64(args, "observer_offset")?.unwrap_or(1.75);
    let target_offset = parse_optional_f64(args, "target_offset")?.unwrap_or(0.0);
    if !observer_offset.is_finite() || !target_offset.is_finite() {
        return Err(ToolError::Validation(
            "offsets must be finite numbers".to_string(),
        ));
    }
    let pair_field = parse_optional_str(args, "pair_field")?.map(str::to_string);
    // Users pass a 1-based band; `Raster::get` is 0-based.
    let band_1based = parse_optional_f64(args, "band")?
        .map(|v| v as i64)
        .unwrap_or(1);
    if band_1based < 1 {
        return Err(ToolError::Validation("'band' must be >= 1".to_string()));
    }
    Ok(Params {
        observer_offset,
        target_offset,
        pair_field,
        band: (band_1based - 1) as isize,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
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

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{DataType, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a small projected DEM from a row-major elevation buffer.
    fn dem_from(cols: usize, rows: usize, cell: f64, data: Vec<f64>) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: cell,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: wbraster::CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, data[row * cols + col])
                    .unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn point_layer(name: &str, pts: &[(f64, f64, &str)]) -> String {
        let mut l = Layer::new(name)
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("id", FieldType::Text));
        for (x, y, id) in pts {
            l.add_feature(Some(Geometry::point(*x, *y)), &[("id", (*id).into())])
                .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = LineOfSightTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// A flat plane: the target is always visible.
    #[test]
    fn flat_terrain_is_visible() {
        // 10x1 flat DEM at elevation 100, 1 unit cells.
        let dem = dem_from(10, 1, 1.0, vec![100.0; 10]);
        let obs = point_layer("obs", &[(0.5, 0.5, "A")]);
        let tgt = point_layer("tgt", &[(9.5, 0.5, "A")]);
        let (out, _layer) = run(json!({
            "dem": dem, "observers": obs, "targets": tgt,
            "observer_offset": 2.0, "target_offset": 2.0,
        }));
        assert_eq!(out.outputs["pair_count"], json!(1));
        assert_eq!(out.outputs["visible_pairs"], json!(1));
    }

    /// A ridge between observer and target blocks the view.
    #[test]
    fn ridge_blocks_the_view() {
        // Elevations: flat 100 except a tall 300 spike in the middle (col 5).
        let mut data = vec![100.0; 11];
        data[5] = 300.0;
        let dem = dem_from(11, 1, 1.0, data);
        let obs = point_layer("obs", &[(0.5, 0.5, "A")]);
        let tgt = point_layer("tgt", &[(10.5, 0.5, "A")]);
        let (out, layer) = run(json!({
            "dem": dem, "observers": obs, "targets": tgt,
            "observer_offset": 2.0, "target_offset": 2.0,
        }));
        assert_eq!(out.outputs["visible_pairs"], json!(0), "ridge should block");
        // There must be both a visible run (before the ridge) and an obstructed
        // run (after it).
        let vidx = layer.schema.field_index("visible").unwrap();
        let vis_vals: Vec<i64> = layer
            .iter()
            .map(|f| f.attributes[vidx].as_i64().unwrap())
            .collect();
        assert!(vis_vals.contains(&1), "expected a visible near segment");
        assert!(vis_vals.contains(&0), "expected an obstructed far segment");
    }

    /// Raising the observer high above the ridge restores the view.
    #[test]
    fn tall_observer_sees_over_ridge() {
        let mut data = vec![100.0; 11];
        data[5] = 150.0;
        let dem = dem_from(11, 1, 1.0, data);
        let obs = point_layer("obs", &[(0.5, 0.5, "A")]);
        let tgt = point_layer("tgt", &[(10.5, 0.5, "A")]);
        let (out, _l) = run(json!({
            "dem": dem, "observers": obs, "targets": tgt,
            "observer_offset": 500.0, "target_offset": 2.0,
        }));
        assert_eq!(out.outputs["visible_pairs"], json!(1));
    }

    /// pair_field matches observers and targets one-to-one.
    #[test]
    fn pair_field_matches_by_key() {
        let dem = dem_from(10, 10, 1.0, vec![100.0; 100]);
        let obs = point_layer("obs", &[(0.5, 0.5, "A"), (0.5, 9.5, "B")]);
        let tgt = point_layer("tgt", &[(9.5, 0.5, "A"), (9.5, 9.5, "B")]);
        let (out, _l) = run(json!({
            "dem": dem, "observers": obs, "targets": tgt, "pair_field": "id",
        }));
        // Two matched pairs (A-A, B-B), not the four of an all-to-all join.
        assert_eq!(out.outputs["pair_count"], json!(2));
    }

    #[test]
    fn rejects_missing_required() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            LineOfSightTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "dem": "d.tif", "observers": "o.geojson" })).is_err());
        assert!(
            bad(json!({ "dem": "d.tif", "observers": "o.geojson", "targets": "t.geojson" }))
                .is_ok()
        );
    }
}
