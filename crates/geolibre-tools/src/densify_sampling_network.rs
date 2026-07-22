//! GeoLibre tool: `densify_sampling_network` — propose new sample sites at the
//! cells of highest kriging prediction error, honoring a minimum-spacing
//! (inhibition) distance.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Densify Sampling Network*
//! (Geostatistical Analyst). GeoLibre already generates sample locations with
//! `create_spatially_balanced_points` and `create_spatial_sampling_locations`,
//! but none of them are *error-guided*: they spread points to fill space, not
//! to target the least-known areas of a surface. This tool consumes a
//! prediction-standard-error raster — the kind produced by ordinary kriging or
//! GeoLibre's `empirical_bayesian_kriging` — and greedily proposes where to add
//! the next field samples so that the largest remaining uncertainty is measured
//! first.
//!
//! ## Method
//! 1. Read the (single-band) prediction-error raster and collect every valid
//!    (non-nodata, finite) cell as a candidate keyed by its standard-error
//!    value and cell-center map coordinate.
//! 2. Optionally restrict candidates to those whose center falls inside a
//!    `mask` polygon layer (e.g. the study area or accessible terrain).
//! 3. Sort candidates by error descending, then greedily accept the largest
//!    remaining error whose cell center is at least `inhibition_distance` map
//!    units from every already-accepted point. Repeat until `count` points are
//!    placed or no admissible candidate remains.
//!
//! Because acceptance is max-error-first, any prefix of the output is itself the
//! set of highest-priority sites, and the inhibition distance prevents new
//! samples from clustering in one hot spot. When the inhibition distance is 0
//! (the default), the tool simply returns the `count` highest-error cells.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::common::load_input_raster;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Proposes new sample locations at the cells of highest kriging prediction
/// error, spaced at least `inhibition_distance` apart.
pub struct DensifySamplingNetworkTool;

impl Tool for DensifySamplingNetworkTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "densify_sampling_network",
            display_name: "Densify Sampling Network",
            summary: "Propose new sample locations at the cells of highest kriging prediction (standard-error) surface, honoring a minimum-spacing inhibition distance, like ArcGIS Densify Sampling Network.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "prediction_error",
                    description: "Input prediction standard-error raster (e.g. a kriging variance/error surface). New samples target its highest-error cells.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output point vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "count",
                    description: "Number of new sample locations to propose.",
                    required: true,
                },
                ToolParamSpec {
                    name: "inhibition_distance",
                    description: "Minimum spacing between proposed samples in map units (default 0 = no inhibition; just the highest-error cells).",
                    required: false,
                },
                ToolParamSpec {
                    name: "mask",
                    description: "Optional polygon layer; only cells whose center falls inside are eligible candidates.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band of the prediction-error raster to read (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args
            .get("prediction_error")
            .and_then(Value::as_str)
            .map(str::trim)
            .unwrap_or("")
            .is_empty()
        {
            return Err(ToolError::Validation(
                "missing required string parameter 'prediction_error'".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let error_path = args
            .get("prediction_error")
            .and_then(Value::as_str)
            .unwrap();
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let raster = load_input_raster(error_path)?;
        let band = (prm.band_1based - 1) as isize;
        if band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "band {} out of range (raster has {} band(s))",
                prm.band_1based, raster.bands
            )));
        }

        // Optional candidate mask polygons.
        let mask_polys: Option<Vec<Region>> = match &prm.mask_path {
            Some(path) => {
                let layer = load_input_layer(path)?;
                let mut polys = Vec::new();
                for feature in &layer.features {
                    if let Some(geom) = feature.geometry.as_ref() {
                        polys.extend(polygon_regions(geom));
                    }
                }
                if polys.is_empty() {
                    return Err(ToolError::Execution(
                        "mask layer has no usable polygon area".to_string(),
                    ));
                }
                Some(polys)
            }
            None => None,
        };

        // Collect valid candidate cells (value, x, y).
        ctx.progress.info("collecting candidate cells");
        let nodata = raster.nodata;
        let rows = raster.rows as isize;
        let cols = raster.cols as isize;
        let mut candidates: Vec<Candidate> = Vec::new();
        for row in 0..rows {
            for col in 0..cols {
                let v = raster.get(band, row, col);
                if v == nodata || !v.is_finite() {
                    continue;
                }
                let x = raster.col_center_x(col);
                let y = raster.row_center_y(row);
                if let Some(polys) = &mask_polys {
                    if !point_in_any(x, y, polys) {
                        continue;
                    }
                }
                candidates.push(Candidate { value: v, x, y });
            }
        }

        if candidates.is_empty() {
            return Err(ToolError::Execution(
                "prediction-error raster has no valid candidate cells (check nodata / mask)"
                    .to_string(),
            ));
        }

        // Sort by error descending (NaNs already filtered out).
        candidates.sort_by(|a, b| b.value.partial_cmp(&a.value).unwrap());

        // Greedy max-error-first selection with an inhibition radius.
        ctx.progress
            .info(&format!("selecting up to {} sample site(s)", prm.count));
        let inhib_sq = prm.inhibition_distance * prm.inhibition_distance;
        let mut picked: Vec<Candidate> = Vec::with_capacity(prm.count);
        for cand in &candidates {
            if picked.len() >= prm.count {
                break;
            }
            if prm.inhibition_distance > 0.0 {
                let too_close = picked.iter().any(|p| {
                    let dx = p.x - cand.x;
                    let dy = p.y - cand.y;
                    dx * dx + dy * dy < inhib_sq
                });
                if too_close {
                    continue;
                }
            }
            picked.push(*cand);
        }

        if picked.is_empty() {
            return Err(ToolError::Execution(
                "no sample sites could be placed — check the inhibition distance".to_string(),
            ));
        }

        // Build the output point layer.
        let mut out = Layer::new("densified_samples").with_geom_type(GeometryType::Point);
        if let Some(epsg) = raster.crs.epsg {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("rank", FieldType::Integer));
        out.add_field(FieldDef::new("pred_error", FieldType::Float));

        for (i, p) in picked.iter().enumerate() {
            out.add_feature(
                Some(Geometry::Point(Coord::xy(p.x, p.y))),
                &[
                    ("rank", FieldValue::Integer((i + 1) as i64)),
                    ("pred_error", FieldValue::Float(p.value)),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed writing point: {e}")))?;
        }

        // Minimum realized spacing between the chosen sites (validation-friendly).
        let mut min_spacing = f64::INFINITY;
        for i in 0..picked.len() {
            for j in (i + 1)..picked.len() {
                let dx = picked[i].x - picked[j].x;
                let dy = picked[i].y - picked[j].y;
                min_spacing = min_spacing.min((dx * dx + dy * dy).sqrt());
            }
        }
        let min_spacing = if min_spacing.is_finite() {
            min_spacing
        } else {
            0.0
        };

        let max_error = picked[0].value;
        let min_error = picked.iter().map(|p| p.value).fold(f64::INFINITY, f64::min);

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("requested".to_string(), json!(prm.count));
        outputs.insert("generated".to_string(), json!(picked.len()));
        outputs.insert("candidates".to_string(), json!(candidates.len()));
        outputs.insert("min_spacing".to_string(), json!(min_spacing));
        outputs.insert("max_pred_error".to_string(), json!(max_error));
        outputs.insert("min_pred_error".to_string(), json!(min_error));
        Ok(ToolRunResult { outputs })
    }
}

/// A candidate cell: its prediction error and cell-center map coordinate.
#[derive(Clone, Copy)]
struct Candidate {
    value: f64,
    x: f64,
    y: f64,
}

/// One polygon part: a list of rings (exterior first, then holes), each ring a
/// chain of `(x, y)` vertices.
type Region = Vec<Vec<(f64, f64)>>;

// ── Geometry helpers (mask polygons) ─────────────────────────────────────────

fn point_in_any(x: f64, y: f64, polys: &[Region]) -> bool {
    polys.iter().any(|rings| point_in_rings(x, y, rings))
}

fn point_in_rings(x: f64, y: f64, rings: &[Vec<(f64, f64)>]) -> bool {
    if rings.is_empty() || !point_in_ring(x, y, &rings[0]) {
        return false;
    }
    !rings[1..].iter().any(|h| point_in_ring(x, y, h))
}

fn point_in_ring(x: f64, y: f64, ring: &[(f64, f64)]) -> bool {
    let n = ring.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = ring[i];
        let (xj, yj) = ring[j];
        if (yi > y) != (yj > y) {
            let xcross = (xj - xi) * (y - yi) / (yj - yi) + xi;
            if x < xcross {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

/// Each polygon part becomes one region: a Vec of rings (exterior first).
fn polygon_regions(geom: &Geometry) -> Vec<Region> {
    let ring_pts =
        |ring: &Ring| -> Vec<(f64, f64)> { ring.coords().iter().map(|c| (c.x, c.y)).collect() };
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            let mut rings = vec![ring_pts(exterior)];
            rings.extend(interiors.iter().map(&ring_pts));
            vec![rings]
        }
        Geometry::MultiPolygon(parts) => parts
            .iter()
            .map(|(ext, holes)| {
                let mut rings = vec![ring_pts(ext)];
                rings.extend(holes.iter().map(&ring_pts));
                rings
            })
            .collect(),
        _ => Vec::new(),
    }
}

// ── Parameters ───────────────────────────────────────────────────────────────

struct Params {
    count: usize,
    inhibition_distance: f64,
    mask_path: Option<String>,
    band_1based: u64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let count = match args.get("count") {
        Some(Value::Number(n)) => n.as_u64().unwrap_or(0) as usize,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'count' must be an integer".into()))?,
        _ => {
            return Err(ToolError::Validation(
                "missing required integer parameter 'count'".to_string(),
            ))
        }
    };
    if count == 0 {
        return Err(ToolError::Validation(
            "'count' must be at least 1".to_string(),
        ));
    }

    let inhibition_distance = match args.get("inhibition_distance") {
        None | Some(Value::Null) => 0.0,
        Some(Value::Number(n)) => n.as_f64().unwrap_or(0.0),
        Some(Value::String(s)) if s.trim().is_empty() => 0.0,
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'inhibition_distance' must be a number".into()))?,
        Some(_) => {
            return Err(ToolError::Validation(
                "'inhibition_distance' must be a number".into(),
            ))
        }
    };
    if inhibition_distance < 0.0 || !inhibition_distance.is_finite() {
        return Err(ToolError::Validation(
            "'inhibition_distance' must be a non-negative, finite number".into(),
        ));
    }

    let band_1based = match args.get("band") {
        None | Some(Value::Null) => 1,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(1).max(1),
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
            .map_err(|_| ToolError::Validation("'band' must be an integer".into()))?
            .max(1),
        Some(_) => return Err(ToolError::Validation("'band' must be a number".into())),
    };

    Ok(Params {
        count,
        inhibition_distance,
        mask_path: parse_optional_str(args, "mask")?.map(str::to_string),
        band_1based,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{memory_store as rmem, DataType, Raster, RasterConfig};
    use wbvector::memory_store as vmem;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Build a `rows × cols` single-band raster from a row-major value slice,
    /// origin (0,0), unit cell size, store it in memory, return its path.
    fn raster_from(rows: usize, cols: usize, cell: f64, values: &[f64], nodata: f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            data_type: DataType::F64,
            cell_size: cell,
            nodata,
            ..Default::default()
        });
        r.crs.epsg = Some(3857);
        for row in 0..rows as isize {
            for col in 0..cols as isize {
                let v = values[row as usize * cols + col as usize];
                r.set(0, row, col, v).unwrap();
            }
        }
        let id = rmem::put_raster(r);
        rmem::make_raster_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = DensifySamplingNetworkTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    // Core property: with no inhibition, the tool returns exactly the N
    // highest-error cells, ranked descending.
    #[test]
    fn picks_highest_error_cells_ranked() {
        // 3x3 with a clear ordering; max at (row2,col2)=9, then 8, 7...
        #[rustfmt::skip]
        let vals = vec![
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0,
        ];
        let path = raster_from(3, 3, 1.0, &vals, -9999.0);
        let (res, layer) = run(json!({ "prediction_error": path, "count": 3 }));
        assert_eq!(res.outputs["generated"].as_u64().unwrap(), 3);
        let errs: Vec<f64> = layer
            .features
            .iter()
            .map(|f| f.get_by_index(1).unwrap().as_f64().unwrap())
            .collect();
        assert_eq!(errs, vec![9.0, 8.0, 7.0]);
        // rank 1 is the global max, at cell center (2.5, 0.5).
        let g = layer.features[0].geometry.as_ref().unwrap();
        if let Geometry::Point(c) = g {
            assert!((c.x - 2.5).abs() < 1e-9);
            assert!((c.y - 0.5).abs() < 1e-9);
        } else {
            panic!("expected point");
        }
    }

    // Inhibition distance is honored: no two chosen points closer than it.
    #[test]
    fn respects_inhibition_distance() {
        // 5x5 grid, unit cells; make top-right corner the hottest region.
        let mut vals = vec![0.0; 25];
        for (i, v) in vals.iter_mut().enumerate() {
            let row = (i / 5) as f64;
            let col = (i % 5) as f64;
            // higher toward top-right
            *v = col - row;
        }
        let path = raster_from(5, 5, 1.0, &vals, -9999.0);
        let (res, layer) = run(json!({
            "prediction_error": path,
            "count": 4,
            "inhibition_distance": 2.5
        }));
        let min_spacing = res.outputs["min_spacing"].as_f64().unwrap();
        assert!(
            min_spacing >= 2.5,
            "min spacing {min_spacing} below inhibition distance"
        );
        // pairwise verify from geometry too
        let pts: Vec<(f64, f64)> = layer
            .features
            .iter()
            .map(|f| match f.geometry.as_ref().unwrap() {
                Geometry::Point(c) => (c.x, c.y),
                _ => panic!(),
            })
            .collect();
        for i in 0..pts.len() {
            for j in (i + 1)..pts.len() {
                let d = ((pts[i].0 - pts[j].0).powi(2) + (pts[i].1 - pts[j].1).powi(2)).sqrt();
                assert!(d >= 2.5 - 1e-9);
            }
        }
    }

    // nodata cells are never chosen (pass-through / exclusion behavior).
    #[test]
    fn skips_nodata_cells() {
        // The single highest raw value is nodata; the real max is 5.0.
        #[rustfmt::skip]
        let vals = vec![
            5.0, 2.0,
            3.0, 999.0, // 999 marked nodata below
        ];
        let path = raster_from(2, 2, 1.0, &vals, 999.0);
        let (res, layer) = run(json!({ "prediction_error": path, "count": 1 }));
        assert_eq!(res.outputs["candidates"].as_u64().unwrap(), 3);
        let e = layer.features[0].get_by_index(1).unwrap().as_f64().unwrap();
        assert_eq!(e, 5.0);
    }

    // Large inhibition distance means fewer points than requested (as many as fit).
    #[test]
    fn caps_when_inhibition_too_large() {
        #[rustfmt::skip]
        let vals = vec![
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0,
        ];
        let path = raster_from(3, 3, 1.0, &vals, -9999.0);
        // Inhibition bigger than the whole extent -> only one point can be placed.
        let (res, _layer) = run(json!({
            "prediction_error": path,
            "count": 5,
            "inhibition_distance": 100.0
        }));
        assert_eq!(res.outputs["requested"].as_u64().unwrap(), 5);
        assert_eq!(res.outputs["generated"].as_u64().unwrap(), 1);
    }

    // Mask restricts candidates to inside the polygon.
    #[test]
    fn mask_restricts_candidates() {
        #[rustfmt::skip]
        let vals = vec![
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0,
        ];
        let path = raster_from(3, 3, 1.0, &vals, -9999.0);
        // Mask covers only the top-left cell center (~0.5, 2.5) -> value 1.0.
        let mut m = Layer::new("m")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        m.add_feature(
            Some(Geometry::polygon(
                vec![
                    Coord::xy(0.0, 2.0),
                    Coord::xy(1.0, 2.0),
                    Coord::xy(1.0, 3.0),
                    Coord::xy(0.0, 3.0),
                ],
                vec![],
            )),
            &[],
        )
        .unwrap();
        let mid = vmem::put_vector(m);
        let mpath = vmem::make_vector_memory_path(&mid);
        let (res, layer) = run(json!({
            "prediction_error": path,
            "count": 5,
            "mask": mpath
        }));
        assert_eq!(res.outputs["candidates"].as_u64().unwrap(), 1);
        assert_eq!(layer.features.len(), 1);
        let e = layer.features[0].get_by_index(1).unwrap().as_f64().unwrap();
        assert_eq!(e, 1.0);
    }

    #[test]
    fn rejects_bad_parameters() {
        let path = raster_from(2, 2, 1.0, &[1.0, 2.0, 3.0, 4.0], -9999.0);
        // missing count
        let a: ToolArgs = serde_json::from_value(json!({ "prediction_error": path })).unwrap();
        assert!(DensifySamplingNetworkTool.validate(&a).is_err());
        // zero count
        let a: ToolArgs =
            serde_json::from_value(json!({ "prediction_error": path, "count": 0 })).unwrap();
        assert!(DensifySamplingNetworkTool.validate(&a).is_err());
        // negative inhibition
        let a: ToolArgs = serde_json::from_value(
            json!({ "prediction_error": path, "count": 1, "inhibition_distance": -1.0 }),
        )
        .unwrap();
        assert!(DensifySamplingNetworkTool.validate(&a).is_err());
        // missing prediction_error
        let a: ToolArgs = serde_json::from_value(json!({ "count": 1 })).unwrap();
        assert!(DensifySamplingNetworkTool.validate(&a).is_err());
    }
}
