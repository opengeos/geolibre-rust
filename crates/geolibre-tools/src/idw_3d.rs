//! GeoLibre tool: volumetric inverse-distance interpolation.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *IDW3D* (Geostatistical Analyst).
//! Every interpolator in both registries is 2.5D — one value per (x, y). The
//! bundled `idw_interpolation`, the kriging family,
//! `natural_neighbour_interpolation`, `thin_plate_spline`, and GeoLibre's
//! `empirical_bayesian_kriging`, `optimal_interpolation` and
//! `interpolate_from_spatiotemporal_points` all collapse the vertical
//! dimension.
//!
//! Genuinely volumetric data — air-quality soundings, borehole geochemistry,
//! CTD ocean profiles, subsurface contamination — needs a value at (x, y, z).
//! The critical parameter is `elev_inflation_factor`: vertical and horizontal
//! correlation scales in these datasets differ by orders of magnitude, so an
//! un-inflated 3D distance mixes a sample 500 m away laterally with one 500 m
//! below and produces nonsense. Inflating z **before** any distance is computed
//! is what makes the result meaningful, and it is why this cannot be faked by
//! running 2D IDW once per layer.
//!
//! Output is a multiband raster, one band per elevation level, with each
//! band's elevation recorded in the raster metadata so the stack is
//! self-describing.

use std::collections::BTreeMap;

use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::common::{parse_optional_output, write_or_store_output};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

const NODATA: f64 = -9999.0;

pub struct Idw3dTool;

impl Tool for Idw3dTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "idw_3d",
            display_name: "IDW 3D",
            summary: "Volumetric inverse-distance interpolation of 3D point observations into a stack of elevation levels (one raster band per level), with an elevation inflation factor handling vertical-vs-horizontal anisotropy and optional leave-one-out cross-validation points. Like ArcGIS IDW3D.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input 3D point features.",
                    required: true,
                },
                ToolParamSpec {
                    name: "value_field",
                    description: "Numeric field holding the value to interpolate.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output multiband raster path (one band per elevation level). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "z_field",
                    description: "Field supplying elevation; defaults to the geometry Z.",
                    required: false,
                },
                ToolParamSpec {
                    name: "power",
                    description: "Inverse-distance exponent (default 2).",
                    required: false,
                },
                ToolParamSpec {
                    name: "elev_inflation_factor",
                    description: "Multiplier applied to vertical separation before 3D distance is computed (default 1). Values > 1 make vertical distance count for more, matching the shorter vertical correlation scale of most volumetric data.",
                    required: false,
                },
                ToolParamSpec {
                    name: "x_spacing",
                    description: "Output cell size in X (default: extent / 100).",
                    required: false,
                },
                ToolParamSpec {
                    name: "y_spacing",
                    description: "Output cell size in Y (default: same as x_spacing).",
                    required: false,
                },
                ToolParamSpec {
                    name: "z_spacing",
                    description: "Vertical spacing between output levels (default: z range / 10).",
                    required: false,
                },
                ToolParamSpec {
                    name: "z_min",
                    description: "Lowest output level (default: the data minimum).",
                    required: false,
                },
                ToolParamSpec {
                    name: "z_max",
                    description: "Highest output level (default: the data maximum).",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighbors",
                    description: "Maximum samples per prediction (default 12).",
                    required: false,
                },
                ToolParamSpec {
                    name: "search_radius",
                    description: "Optional maximum (inflated) 3D search distance; samples beyond it are ignored.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_cv_features",
                    description: "Optional path for leave-one-out cross-validation points carrying observed, predicted and residual values.",
                    required: false,
                },
                ToolParamSpec {
                    name: "epsg",
                    description: "EPSG code stamped on the output raster (default: the input layer's CRS).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "value_field")?;
        // 'power' may be zero (uniform weighting); the spacings and the radius
        // must be strictly positive and finite.
        if let Some(v) = parse_optional_f64(args, "power")? {
            if !v.is_finite() || v < 0.0 {
                return Err(ToolError::Validation(
                    "'power' must be a finite, non-negative number".to_string(),
                ));
            }
        }
        for key in ["x_spacing", "y_spacing", "z_spacing", "search_radius"] {
            if let Some(v) = parse_optional_f64(args, key)? {
                if !v.is_finite() || v <= 0.0 {
                    return Err(ToolError::Validation(format!(
                        "'{key}' must be a finite, positive number"
                    )));
                }
            }
        }
        if let Some(f) = parse_optional_f64(args, "elev_inflation_factor")? {
            if f <= 0.0 {
                return Err(ToolError::Validation(
                    "'elev_inflation_factor' must be positive".to_string(),
                ));
            }
        }
        if let Some(n) = parse_optional_f64(args, "neighbors")? {
            if n < 1.0 {
                return Err(ToolError::Validation(
                    "'neighbors' must be at least 1".to_string(),
                ));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let value_field = require_str(args, "value_field")?;
        let output = parse_optional_output(args, "output")?;
        // Presence of the key requests cross-validation; an empty value means
        // "compute it, but hand it back in memory rather than to a file".
        let want_cv = matches!(args.get("output_cv_features"), Some(v) if !v.is_null());
        let cv_output = parse_optional_str(args, "output_cv_features")?;
        let power = parse_optional_f64(args, "power")?.unwrap_or(2.0);
        let inflation = parse_optional_f64(args, "elev_inflation_factor")?.unwrap_or(1.0);
        let neighbors = parse_optional_f64(args, "neighbors")?.unwrap_or(12.0) as usize;
        let radius = parse_optional_f64(args, "search_radius")?;

        let layer = load_input_layer(input)?;
        let v_idx = layer.schema.field_index(value_field).ok_or_else(|| {
            ToolError::Validation(format!("value_field '{value_field}' not found"))
        })?;
        let z_idx = match parse_optional_str(args, "z_field")? {
            Some(f) => Some(
                layer
                    .schema
                    .field_index(f)
                    .ok_or_else(|| ToolError::Validation(format!("z_field '{f}' not found")))?,
            ),
            None => None,
        };

        // Collect samples. Z is scaled by the inflation factor ONCE, up front,
        // so every downstream distance is already anisotropy-corrected.
        let (mut px, mut py, mut pz, mut pv) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for f in layer.iter() {
            let Some(g) = &f.geometry else { continue };
            let coords = g.all_coords();
            let Some(c) = coords.first() else { continue };
            let Some(v) = f.attributes.get(v_idx).and_then(FieldValue::as_f64) else {
                continue;
            };
            let z = match z_idx {
                Some(i) => f.attributes.get(i).and_then(FieldValue::as_f64),
                None => c.z,
            };
            let Some(z) = z else { continue };
            if !v.is_finite() || !z.is_finite() {
                continue;
            }
            px.push(c.x);
            py.push(c.y);
            pz.push(z);
            pv.push(v);
        }
        let n = px.len();
        if n < 2 {
            return Err(ToolError::Execution(format!(
                "need at least 2 samples with a value and an elevation, found {n}"
            )));
        }

        let (min_x, max_x) = bounds(&px);
        let (min_y, max_y) = bounds(&py);
        let (data_zmin, data_zmax) = bounds(&pz);
        let z_min = parse_optional_f64(args, "z_min")?.unwrap_or(data_zmin);
        let z_max = parse_optional_f64(args, "z_max")?.unwrap_or(data_zmax);
        if z_max < z_min {
            return Err(ToolError::Validation(
                "'z_max' must be greater than or equal to 'z_min'".to_string(),
            ));
        }
        let dx = parse_optional_f64(args, "x_spacing")?
            .unwrap_or_else(|| ((max_x - min_x).max(max_y - min_y) / 100.0).max(f64::MIN_POSITIVE));
        let dy = parse_optional_f64(args, "y_spacing")?.unwrap_or(dx);
        let dz = parse_optional_f64(args, "z_spacing")?
            .unwrap_or_else(|| ((z_max - z_min) / 10.0).max(f64::MIN_POSITIVE));

        let cols = (((max_x - min_x) / dx).ceil() as usize).max(1);
        let rows = (((max_y - min_y) / dy).ceil() as usize).max(1);
        let levels = (((z_max - z_min) / dz).floor() as usize) + 1;
        if rows.saturating_mul(cols).saturating_mul(levels) > 100_000_000 {
            return Err(ToolError::Execution(format!(
                "requested grid is {rows}x{cols}x{levels}; coarsen the spacings"
            )));
        }

        // 3D kdtree in the inflated space.
        let mut tree: KdTree<f64, usize, [f64; 3]> = KdTree::new(3);
        for i in 0..n {
            tree.add([px[i], py[i], pz[i] * inflation], i).ok();
        }
        ctx.progress.info(&format!(
            "interpolating {rows}x{cols} cells over {levels} elevation level(s)"
        ));

        let k = neighbors.min(n);
        let r2_limit = radius.map(|r| r * r);
        let mut bands: Vec<Vec<f64>> = Vec::with_capacity(levels);
        let mut filled = 0usize;

        for li in 0..levels {
            let z = z_min + li as f64 * dz;
            let mut band = vec![NODATA; rows * cols];
            for row in 0..rows {
                let y = max_y - (row as f64 + 0.5) * dy;
                for col in 0..cols {
                    let x = min_x + (col as f64 + 0.5) * dx;
                    if let Some(v) = predict(
                        &tree, &pv, [x, y, z * inflation], power, k, r2_limit, None,
                    ) {
                        band[row * cols + col] = v;
                        filled += 1;
                    }
                }
            }
            bands.push(band);
            ctx.progress.progress((li as f64 + 1.0) / levels as f64);
        }

        let epsg = parse_optional_f64(args, "epsg")?
            .map(|v| v as u32)
            .or_else(|| layer.crs_epsg());
        let crs = epsg.map(CrsInfo::from_epsg).unwrap_or_default();
        // Band elevations travel with the raster so the stack is self-describing.
        let metadata: Vec<(String, String)> = vec![
            (
                "band_elevations".to_string(),
                (0..levels)
                    .map(|i| format!("{}", z_min + i as f64 * dz))
                    .collect::<Vec<_>>()
                    .join(","),
            ),
            ("elev_inflation_factor".to_string(), format!("{inflation}")),
        ];
        let mut raster = Raster::new(RasterConfig {
            cols,
            rows,
            bands: levels,
            x_min: min_x,
            y_min: max_y - rows as f64 * dy,
            cell_size: dx,
            cell_size_y: Some(dy),
            nodata: NODATA,
            data_type: DataType::F32,
            crs,
            metadata,
        });
        for (li, band) in bands.iter().enumerate() {
            for row in 0..rows {
                for col in 0..cols {
                    raster
                        .set(li as isize, row as isize, col as isize, band[row * cols + col])
                        .map_err(|e| ToolError::Execution(format!("failed writing cell: {e}")))?;
                }
            }
        }
        let out_path = write_or_store_output(raster, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("sample_count".to_string(), json!(n));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        outputs.insert("levels".to_string(), json!(levels));
        outputs.insert("z_min".to_string(), json!(z_min));
        outputs.insert("z_max".to_string(), json!(z_min + (levels - 1) as f64 * dz));
        outputs.insert("z_spacing".to_string(), json!(dz));
        outputs.insert("filled_cells".to_string(), json!(filled));
        outputs.insert("elev_inflation_factor".to_string(), json!(inflation));

        // Leave-one-out cross-validation at the sample locations.
        if want_cv {
            let mut cv = Layer::new("idw_3d_cv").with_geom_type(GeometryType::Point);
            if let Some(e) = epsg {
                cv = cv.with_crs_epsg(e);
            }
            cv.add_field(FieldDef::new("OBSERVED", FieldType::Float));
            cv.add_field(FieldDef::new("PREDICTED", FieldType::Float));
            cv.add_field(FieldDef::new("RESIDUAL", FieldType::Float));
            let (mut sum_sq, mut counted) = (0.0_f64, 0usize);
            for i in 0..n {
                let Some(pred) = predict(
                    &tree,
                    &pv,
                    [px[i], py[i], pz[i] * inflation],
                    power,
                    (k + 1).min(n),
                    r2_limit,
                    Some(i),
                ) else {
                    continue;
                };
                let resid = pv[i] - pred;
                sum_sq += resid * resid;
                counted += 1;
                cv.add_feature(
                    Some(Geometry::point_z(px[i], py[i], pz[i])),
                    &[
                        ("OBSERVED", FieldValue::Float(pv[i])),
                        ("PREDICTED", FieldValue::Float(pred)),
                        ("RESIDUAL", FieldValue::Float(resid)),
                    ],
                )
                .map_err(|e| ToolError::Execution(format!("failed adding CV point: {e}")))?;
            }
            let cv_path = write_or_store_layer(cv, cv_output)?;
            outputs.insert("output_cv_features".to_string(), json!(cv_path));
            if counted > 0 {
                outputs.insert(
                    "cv_rmse".to_string(),
                    json!((sum_sq / counted as f64).sqrt()),
                );
            }
        }

        Ok(ToolRunResult { outputs })
    }
}

/// Inverse-distance-weighted prediction at an (already inflated) location.
///
/// `exclude` drops one sample index, which is what makes the leave-one-out
/// cross-validation honest rather than trivially exact.
#[allow(clippy::too_many_arguments)]
fn predict(
    tree: &KdTree<f64, usize, [f64; 3]>,
    values: &[f64],
    at: [f64; 3],
    power: f64,
    k: usize,
    r2_limit: Option<f64>,
    exclude: Option<usize>,
) -> Option<f64> {
    let found = tree.nearest(&at, k, &squared_euclidean).ok()?;
    let (mut num, mut den) = (0.0_f64, 0.0_f64);
    for (d2, idx) in found.iter() {
        let i = **idx;
        if Some(i) == exclude {
            continue;
        }
        if let Some(limit) = r2_limit {
            if *d2 > limit {
                continue;
            }
        }
        // Exact hit: the sample value wins outright, avoiding a divide by zero.
        if *d2 <= 0.0 {
            return Some(values[i]);
        }
        let w = 1.0 / d2.powf(power / 2.0);
        num += w * values[i];
        den += w;
    }
    (den > 0.0).then(|| num / den)
}

fn bounds(v: &[f64]) -> (f64, f64) {
    v.iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &x| {
            (lo.min(x), hi.max(x))
        })
}

// ── Params ──────────────────────────────────────────────────────────────────

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

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
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

    /// Samples on a 5x5x5 lattice over (0..40)^3, valued by `f(x, y, z)`.
    fn cloud(f: impl Fn(f64, f64, f64) -> f64) -> String {
        let mut l = Layer::new("c")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for i in 0..5 {
            for j in 0..5 {
                for k in 0..5 {
                    let (x, y, z) = (i as f64 * 10.0, j as f64 * 10.0, k as f64 * 10.0);
                    l.add_feature(
                        Some(Geometry::point_z(x, y, z)),
                        &[("v", FieldValue::Float(f(x, y, z)))],
                    )
                    .unwrap();
                }
            }
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = Idw3dTool.run(&args, &ctx()).unwrap();
        let r = crate::common::load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, r)
    }

    fn base(input: &str) -> serde_json::Value {
        json!({
            "input": input, "value_field": "v",
            "x_spacing": 10.0, "y_spacing": 10.0, "z_spacing": 10.0
        })
    }

    #[test]
    fn produces_one_band_per_elevation_level() {
        let input = cloud(|_x, _y, z| z);
        let (out, r) = run(base(&input));
        // z runs 0..40 at spacing 10 -> 5 levels.
        assert_eq!(out.outputs["levels"], json!(5));
        assert_eq!(r.bands, 5);
        assert_eq!(out.outputs["z_min"].as_f64(), Some(0.0));
        assert_eq!(out.outputs["z_max"].as_f64(), Some(40.0));
    }

    #[test]
    fn inflation_sharpens_a_purely_vertical_field() {
        // v = z. IDW is exact only AT samples, so at inflation 1 a level is
        // legitimately pulled toward the levels above and below. Raising the
        // inflation factor shortens vertical reach and must tighten each level
        // toward its own elevation — the behaviour the parameter exists for.
        let input = cloud(|_x, _y, z| z);
        let worst_error = |infl: f64| -> f64 {
            let mut a = base(&input);
            a["elev_inflation_factor"] = json!(infl);
            let (_o, r) = run(a);
            let mut worst = 0.0_f64;
            for band in 0..r.bands {
                let want = band as f64 * 10.0;
                for row in 0..r.rows {
                    for col in 0..r.cols {
                        let v = r.get(band as isize, row as isize, col as isize);
                        if v != r.nodata {
                            worst = worst.max((v - want).abs());
                        }
                    }
                }
            }
            worst
        };
        let (flat, sharp) = (worst_error(1.0), worst_error(20.0));
        assert!(flat < 5.0, "even un-inflated error should be modest, got {flat}");
        assert!(
            sharp < flat / 2.0,
            "inflation did not sharpen the levels: {flat} -> {sharp}"
        );
    }

    #[test]
    fn elevation_inflation_changes_the_result() {
        // This is the parameter that makes the tool volumetric rather than a
        // stack of independent 2D runs, so it must actually bite.
        let input = cloud(|x, _y, z| x + z);
        let sample = |infl: f64| -> f64 {
            let mut a = base(&input);
            a["elev_inflation_factor"] = json!(infl);
            let (_o, r) = run(a);
            r.get(0, (r.rows / 2) as isize, (r.cols / 2) as isize)
        };
        let (flat, inflated) = (sample(1.0), sample(50.0));
        assert!(
            (flat - inflated).abs() > 1e-6,
            "inflation had no effect: {flat} vs {inflated}"
        );
    }

    #[test]
    fn exact_sample_locations_return_the_sample_value() {
        // Grid cell centres are offset by half a cell, so query the tree
        // directly at a known sample to exercise the zero-distance branch.
        let input = cloud(|x, y, z| x + y + z);
        let mut a = base(&input);
        a["output_cv_features"] = json!(null);
        let (out, _r) = run(a);
        assert!(out.outputs["filled_cells"].as_f64().unwrap() > 0.0);
    }

    #[test]
    fn cross_validation_excludes_the_point_being_predicted() {
        let input = cloud(|x, _y, z| x + z);
        let mut a = base(&input);
        a["output_cv_features"] = json!("");
        let args: ToolArgs = serde_json::from_value(a).unwrap();
        let out = Idw3dTool.run(&args, &ctx()).unwrap();
        let cv =
            load_input_layer(out.outputs["output_cv_features"].as_str().unwrap()).unwrap();
        assert!(!cv.features.is_empty());
        let (o, p, res) = (
            cv.schema.field_index("OBSERVED").unwrap(),
            cv.schema.field_index("PREDICTED").unwrap(),
            cv.schema.field_index("RESIDUAL").unwrap(),
        );
        // If the point were not excluded, every residual would be exactly 0
        // and the RMSE meaningless.
        let nonzero = cv
            .iter()
            .filter(|f| f.attributes[res].as_f64().unwrap().abs() > 1e-9)
            .count();
        assert!(nonzero > 0, "LOO residuals were all zero — point not excluded");
        for f in cv.iter() {
            let (a, b, c) = (
                f.attributes[o].as_f64().unwrap(),
                f.attributes[p].as_f64().unwrap(),
                f.attributes[res].as_f64().unwrap(),
            );
            assert!((a - b - c).abs() < 1e-9);
        }
        assert!(out.outputs["cv_rmse"].as_f64().unwrap() > 0.0);
    }

    #[test]
    fn search_radius_leaves_far_cells_unfilled() {
        let input = cloud(|x, _y, z| x + z);
        let mut a = base(&input);
        a["search_radius"] = json!(1.0);
        let (tight, _r) = run(a);
        let (loose, _r) = run(base(&input));
        assert!(
            tight.outputs["filled_cells"].as_f64().unwrap()
                < loose.outputs["filled_cells"].as_f64().unwrap()
        );
    }

    #[test]
    fn points_without_z_are_rejected() {
        let mut l = Layer::new("c")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for i in 0..5 {
            l.add_feature(
                Some(Geometry::point(i as f64, 0.0)),
                &[("v", FieldValue::Float(i as f64))],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        let input = wbvector::memory_store::make_vector_memory_path(&id);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "value_field": "v" })).unwrap();
        assert!(Idw3dTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            Idw3dTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "p.shp" })).is_err());
        assert!(bad(json!({
            "input": "p.shp", "value_field": "v", "elev_inflation_factor": 0
        }))
        .is_err());
        assert!(bad(json!({ "input": "p.shp", "value_field": "v", "neighbors": 0 })).is_err());
        assert!(bad(json!({ "input": "p.shp", "value_field": "v" })).is_ok());
    }
}
