//! GeoLibre tool: model-based reaggregation between incompatible polygon
//! geographies.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Areal Interpolation Layer To
//! Polygons* (Geostatistical Analyst).
//!
//! ## Why `apportion_polygon` is not this
//!
//! `apportion_polygon` (round 3) does **area-weighted** apportionment: it
//! assumes the underlying quantity is spread uniformly inside each source
//! polygon. That is the classic modifiable-areal-unit-problem failure — it
//! puts as much population in a census tract's lake as in its housing.
//!
//! Areal interpolation instead treats the observations as *areal averages of a
//! smooth underlying surface*, fits a variogram to that surface, and
//! re-integrates it over the targets. Two things follow that area weighting
//! cannot give: neighbouring source polygons inform each target, and every
//! prediction carries a **standard error**.
//!
//! The bundled kriging suite (`ordinary_kriging`, `empirical_bayesian_kriging`,
//! `local_kriging`, …) is all point-support; none accepts polygon support.
//!
//! ## Method
//!
//! Each polygon is discretised into a deterministic point set. Area-to-area
//! covariances are then averages of point-to-point covariances from the fitted
//! variogram, so the whole problem reduces to the ordinary-kriging normal
//! equations `kriging_common` already solves — that is the whole trick.
//!
//! The discretisation is a **regular lattice clipped to the polygon**, seeded
//! by the polygon's own bounding box. No RNG at all, so results are
//! reproducible and the WASM constraint is satisfied by construction rather
//! than by seeding a generator.
//!
//! ## Mass preservation
//!
//! For `field_type = count`, predictions are rescaled so the target total
//! matches the source total wherever the two geographies cover the same area.
//! Kriging is not mass-preserving on its own, and a population reaggregation
//! that quietly invents or loses people is worse than useless — so the
//! adjustment is applied and its factor reported rather than left implicit.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::args_common::{bool_or, choice_or, req_str, usize_or};
use crate::kriging_common::{fit_variogram, krige_matrix, Variogram, VariogramModel};
use crate::vector_common::{
    geometry_contains_point, load_input_layer, parse_optional_str, write_or_store_layer,
};

pub struct ArealInterpolationTool;

impl Tool for ArealInterpolationTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "areal_interpolation",
            display_name: "Areal Interpolation",
            summary: "Predicts aggregated values for target polygons from source polygons with a different geography, by fitting a smooth underlying surface and re-integrating it rather than weighting by area (ArcGIS Areal Interpolation Layer To Polygons). apportion_polygon assumes the quantity is uniform inside each source polygon — the classic modifiable-areal-unit-problem failure — and gives no standard errors; the bundled kriging suite is entirely point-support and cannot accept polygon observations.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Source polygons carrying the measured field.",
                    required: true,
                },
                ToolParamSpec {
                    name: "target",
                    description: "Target polygons to predict values for.",
                    required: true,
                },
                ToolParamSpec {
                    name: "field",
                    description: "Numeric value field on the source polygons.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Target polygons with PREDICTED, PRED_STD_ERR and the diagnostics appended. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "field_type",
                    description: "'count' (default): an extensive total, rescaled to preserve the source sum. 'average': an intensive mean, predicted directly with no rescaling.",
                    required: false,
                },
                ToolParamSpec {
                    name: "discretization",
                    description: "Approximate number of interior points used to represent each polygon's area support (default 25). Higher is more accurate and much slower.",
                    required: false,
                },
                ToolParamSpec {
                    name: "model",
                    description: "Variogram model: 'exponential' (default), 'spherical', or 'gaussian'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "preserve_total",
                    description: "Rescale count predictions so the target total matches the source total (default true for 'count', always off for 'average').",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        req_str(args, "target")?;
        req_str(args, "field")?;
        choice_or(args, "field_type", &["count", "average"], "count")?;
        VariogramModel::parse(choice_or(
            args,
            "model",
            &["exponential", "spherical", "gaussian"],
            "exponential",
        )?)?;
        let d = usize_or(args, "discretization", 25)?;
        if !(1..=2000).contains(&d) {
            return Err(ToolError::Validation(
                "'discretization' must be between 1 and 2000".to_string(),
            ));
        }
        bool_or(args, "preserve_total", true)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let field = req_str(args, "field")?.to_string();
        let is_count = choice_or(args, "field_type", &["count", "average"], "count")? == "count";
        let model = VariogramModel::parse(choice_or(
            args,
            "model",
            &["exponential", "spherical", "gaussian"],
            "exponential",
        )?)?;
        let discretization = usize_or(args, "discretization", 25)?;
        let preserve_total = is_count && bool_or(args, "preserve_total", true)?;
        let output = parse_optional_str(args, "output")?;

        let source = load_input_layer(req_str(args, "input")?)?;
        let field_idx = source.schema.field_index(&field).ok_or_else(|| {
            ToolError::Validation(format!("field '{field}' not found on the source layer"))
        })?;

        // Source supports: a point cloud per polygon, plus the observed value.
        // For counts the spatially smooth quantity is the *density*, not the
        // total — a big tract and a small one with the same population are not
        // the same observation — so counts are divided by AREA and multiplied
        // back by area at the end.
        //
        // Per unit area, not per support point: the point count is an artefact
        // of the discretisation lattice, so normalising by it would make the
        // answer depend on `discretization` instead of on geography.
        let mut supports: Vec<Vec<(f64, f64)>> = Vec::new();
        let mut observed: Vec<f64> = Vec::new();
        let mut totals: Vec<f64> = Vec::new();
        for feature in source.iter() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let Some(value) = numeric(&feature.attributes[field_idx]) else {
                continue;
            };
            let pts = discretize(geom, discretization);
            let area = polygon_area(geom);
            if pts.is_empty() || (is_count && area <= 0.0) {
                continue;
            }
            let intensity = if is_count { value / area } else { value };
            supports.push(pts);
            observed.push(intensity);
            totals.push(value);
        }
        if supports.len() < 2 {
            return Err(ToolError::Execution(format!(
                "need at least 2 usable source polygons, got {}",
                supports.len()
            )));
        }

        // Fit the variogram on the support centroids: the observations are
        // areal, so this is an approximation, and saying so is better than
        // implying a point-support fit the data cannot supply.
        let centroids: Vec<(f64, f64)> = supports.iter().map(|p| centroid(p)).collect();
        let vg = fit_variogram(&centroids, &observed, model, 12);
        let c00 = vg.covariance(0.0);
        // Data-to-data covariances on BLOCK support, computed once. Pairing a
        // point-support left-hand side with block-support right-hand sides
        // would make the system inconsistent, and the weights would then fail
        // to reproduce even a source polygon's own value.
        let n_src = supports.len();
        let mut data_cov = vec![0.0_f64; n_src * n_src];
        for i in 0..n_src {
            for j in i..n_src {
                let c = block_covariance(&supports[i], &supports[j], &vg);
                data_cov[i * n_src + j] = c;
                data_cov[j * n_src + i] = c;
            }
        }
        ctx.progress.info(&format!(
            "{} source polygon(s), {} variogram: nugget {:.4}, sill {:.4}, range {:.4}",
            supports.len(),
            model.label(),
            vg.nugget,
            vg.partial_sill,
            vg.range
        ));

        let target = load_input_layer(req_str(args, "target")?)?;
        let mut out = Layer::new("areal_interpolation");
        out.geom_type = target.geom_type;
        out.crs = target.crs.clone();
        for f in target.schema.fields() {
            out.add_field(f.clone());
        }
        out.add_field(FieldDef::new("PREDICTED", FieldType::Float));
        out.add_field(FieldDef::new("PRED_STD_ERR", FieldType::Float));
        out.add_field(FieldDef::new("AREA", FieldType::Float));

        let names: Vec<String> = target
            .schema
            .fields()
            .iter()
            .map(|f| f.name.clone())
            .collect();

        // First pass: predict.
        let mut rows: Vec<(Option<Geometry>, Vec<FieldValue>, Option<(f64, f64, f64)>)> =
            Vec::new();
        let total_n = target.iter().count().max(1);
        for (k, feature) in target.iter().enumerate() {
            let attrs: Vec<FieldValue> = names
                .iter()
                .enumerate()
                .map(|(i, _)| {
                    feature
                        .attributes
                        .get(i)
                        .cloned()
                        .unwrap_or(FieldValue::Null)
                })
                .collect();
            let prediction = feature.geometry.as_ref().and_then(|g| {
                let pts = discretize(g, discretization);
                let area = polygon_area(g);
                if pts.is_empty() || (is_count && area <= 0.0) {
                    return None;
                }
                // Block-to-block covariance: the average point-to-point
                // covariance between the two supports. This is the step that
                // makes polygon observations usable at all.
                let rhs = |i: usize| block_covariance(&supports[i], &pts, &vg);
                let c_block = block_covariance(&pts, &pts, &vg).min(c00);
                krige_matrix(
                    n_src,
                    &observed,
                    |i, j| data_cov[i * n_src + j],
                    rhs,
                    c_block,
                )
                .map(|(p, var)| (p, var.sqrt(), area))
            });
            rows.push((feature.geometry.clone(), attrs, prediction));
            ctx.progress.progress((k as f64 + 1.0) / total_n as f64);
        }

        // Convert densities back to totals for counts.
        let mut raw: Vec<Option<f64>> = rows
            .iter()
            .map(|(_, _, p)| p.map(|(v, _, area)| if is_count { v * area } else { v }))
            .collect();

        // Mass preservation. Kriging is not mass-preserving, and a population
        // reaggregation that invents or loses people is worse than useless.
        let source_total: f64 = totals.iter().sum();
        let predicted_total: f64 = raw.iter().flatten().sum();
        let scale = if preserve_total && predicted_total.abs() > 1e-12 {
            source_total / predicted_total
        } else {
            1.0
        };
        if (scale - 1.0).abs() > f64::EPSILON {
            for v in raw.iter_mut().flatten() {
                *v *= scale;
            }
        }

        let mut resolved = 0_u64;
        for (k, (geom, attrs, pred)) in rows.into_iter().enumerate() {
            let mut a: Vec<(&str, FieldValue)> = names
                .iter()
                .zip(attrs.iter())
                .map(|(n, v)| (n.as_str(), v.clone()))
                .collect();
            match (raw[k], pred) {
                (Some(value), Some((_, se, area))) => {
                    resolved += 1;
                    a.push(("PREDICTED", FieldValue::Float(value)));
                    // The standard error must carry the same conversion the
                    // prediction did, or it would be in the wrong units.
                    let se_scaled = if is_count { se * area * scale } else { se };
                    a.push(("PRED_STD_ERR", FieldValue::Float(se_scaled)));
                    a.push(("AREA", FieldValue::Float(area)));
                }
                _ => {
                    a.push(("PREDICTED", FieldValue::Null));
                    a.push(("PRED_STD_ERR", FieldValue::Null));
                    a.push(("AREA", FieldValue::Float(0.0)));
                }
            }
            out.add_feature(geom, &a)
                .map_err(|e| ToolError::Execution(e.to_string()))?;
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("source_count".to_string(), json!(supports.len()));
        outputs.insert("resolved_targets".to_string(), json!(resolved));
        outputs.insert("source_total".to_string(), json!(source_total));
        outputs.insert(
            "predicted_total".to_string(),
            json!(raw.iter().flatten().sum::<f64>()),
        );
        outputs.insert("mass_scale".to_string(), json!(scale));
        outputs.insert("model".to_string(), json!(model.label()));
        outputs.insert("range".to_string(), json!(vg.range));
        Ok(ToolRunResult { outputs })
    }
}

/// Mean covariance between two point supports.
fn block_covariance(a: &[(f64, f64)], b: &[(f64, f64)], vg: &Variogram) -> f64 {
    if a.is_empty() || b.is_empty() {
        return 0.0;
    }
    let mut acc = 0.0;
    for p in a {
        for q in b {
            acc += vg.covariance(crate::kriging_common::dist(*p, *q));
        }
    }
    acc / (a.len() * b.len()) as f64
}

/// Discretises a polygon into roughly `target` interior points.
///
/// A regular lattice clipped to the polygon, sized from the polygon's own
/// bounding box: fully deterministic, so no RNG is needed and the WASM
/// constraint holds by construction. Falls back to the centroid for polygons
/// too small or thin for the lattice to catch.
fn discretize(geom: &Geometry, target: usize) -> Vec<(f64, f64)> {
    let Some((min_x, min_y, max_x, max_y)) = bbox(geom) else {
        return Vec::new();
    };
    let (w, h) = (max_x - min_x, max_y - min_y);
    if w <= 0.0 || h <= 0.0 {
        return Vec::new();
    }
    // Aspect-aware lattice so a long thin polygon is not sampled by a single
    // column of points.
    let aspect = (w / h).sqrt();
    let nx = ((target as f64).sqrt() * aspect).ceil().max(1.0) as usize;
    let ny = ((target as f64).sqrt() / aspect).ceil().max(1.0) as usize;

    let mut pts = Vec::with_capacity(target);
    for j in 0..ny {
        let y = min_y + h * (j as f64 + 0.5) / ny as f64;
        for i in 0..nx {
            let x = min_x + w * (i as f64 + 0.5) / nx as f64;
            if geometry_contains_point(geom, x, y) {
                pts.push((x, y));
            }
        }
    }
    if pts.is_empty() {
        // A sliver the lattice missed entirely still deserves a support.
        pts.push((min_x + w / 2.0, min_y + h / 2.0));
    }
    pts
}

fn bbox(geom: &Geometry) -> Option<(f64, f64, f64, f64)> {
    let mut min = (f64::INFINITY, f64::INFINITY);
    let mut max = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut seen = false;
    let mut take = |cs: &Vec<Coord>| {
        for c in cs {
            seen = true;
            min.0 = min.0.min(c.x);
            min.1 = min.1.min(c.y);
            max.0 = max.0.max(c.x);
            max.1 = max.1.max(c.y);
        }
    };
    match geom {
        Geometry::Polygon { exterior, .. } => take(&exterior.0),
        Geometry::MultiPolygon(parts) => {
            for (ext, _) in parts {
                take(&ext.0);
            }
        }
        _ => return None,
    }
    seen.then_some((min.0, min.1, max.0, max.1))
}

/// Planar polygon area (exterior minus holes) by the shoelace formula.
fn polygon_area(geom: &Geometry) -> f64 {
    let ring_area = |cs: &Vec<Coord>| -> f64 {
        let n = cs.len();
        if n < 3 {
            return 0.0;
        }
        let mut acc = 0.0;
        for i in 0..n {
            let a = &cs[i];
            let b = &cs[(i + 1) % n];
            acc += a.x * b.y - b.x * a.y;
        }
        (acc / 2.0).abs()
    };
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => (ring_area(&exterior.0) - interiors.iter().map(|r| ring_area(&r.0)).sum::<f64>())
            .max(0.0),
        Geometry::MultiPolygon(parts) => parts
            .iter()
            .map(|(ext, holes)| {
                (ring_area(&ext.0) - holes.iter().map(|r| ring_area(&r.0)).sum::<f64>()).max(0.0)
            })
            .sum(),
        _ => 0.0,
    }
}

fn centroid(pts: &[(f64, f64)]) -> (f64, f64) {
    let n = pts.len().max(1) as f64;
    (
        pts.iter().map(|p| p.0).sum::<f64>() / n,
        pts.iter().map(|p| p.1).sum::<f64>() / n,
    )
}

fn numeric(v: &FieldValue) -> Option<f64> {
    match v {
        FieldValue::Float(f) => f.is_finite().then_some(*f),
        FieldValue::Integer(i) => Some(*i as f64),
        FieldValue::Text(s) => s.trim().parse().ok(),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, GeometryType, Ring};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Geometry {
        Geometry::Polygon {
            exterior: Ring::new(vec![
                Coord::xy(x0, y0),
                Coord::xy(x1, y0),
                Coord::xy(x1, y1),
                Coord::xy(x0, y1),
            ]),
            interiors: vec![],
        }
    }

    /// Polygons with a value field.
    fn valued(cells: Vec<(Geometry, f64)>) -> String {
        let mut l = Layer::new("src");
        l.geom_type = Some(GeometryType::Polygon);
        l.add_field(FieldDef::new("pop", FieldType::Float));
        for (g, v) in cells {
            l.add_feature(Some(g), &[("pop", FieldValue::Float(v))])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn plain(geoms: Vec<Geometry>) -> String {
        let mut l = Layer::new("tgt");
        l.geom_type = Some(GeometryType::Polygon);
        for g in geoms {
            l.add_feature(Some(g), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = ArealInterpolationTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    fn col(layer: &Layer, name: &str) -> Vec<Option<f64>> {
        let i = layer.schema.field_index(name).unwrap();
        layer
            .iter()
            .map(|f| match &f.attributes[i] {
                FieldValue::Float(v) => Some(*v),
                FieldValue::Integer(v) => Some(*v as f64),
                _ => None,
            })
            .collect()
    }

    /// A 4x1 strip of source cells with a linear gradient.
    fn gradient_source() -> String {
        valued(
            (0..4)
                .map(|i| {
                    let x = i as f64 * 10.0;
                    (rect(x, 0.0, x + 10.0, 10.0), 100.0 + i as f64 * 50.0)
                })
                .collect(),
        )
    }

    #[test]
    fn counts_are_mass_preserving_over_a_matching_extent() {
        // The property that makes this usable for population: the target total
        // must equal the source total when both cover the same ground.
        let (_, res) = run(json!({
            "input": gradient_source(),
            "target": plain(vec![rect(0.0, 0.0, 20.0, 10.0), rect(20.0, 0.0, 40.0, 10.0)]),
            "field": "pop",
        }));
        let src = res.outputs["source_total"].as_f64().unwrap();
        let got = res.outputs["predicted_total"].as_f64().unwrap();
        assert!((got - src).abs() < 1e-6, "source {src}, predicted {got}");
    }

    #[test]
    fn the_prediction_follows_the_underlying_gradient() {
        // The western half must come out lower than the eastern half.
        let (out, _) = run(json!({
            "input": gradient_source(),
            "target": plain(vec![rect(0.0, 0.0, 20.0, 10.0), rect(20.0, 0.0, 40.0, 10.0)]),
            "field": "pop",
        }));
        let p = col(&out, "PREDICTED");
        let (west, east) = (p[0].unwrap(), p[1].unwrap());
        assert!(west < east, "west {west} should be below east {east}");
    }

    #[test]
    fn every_prediction_carries_a_standard_error() {
        // The thing apportion_polygon cannot give at all.
        let (out, _) = run(json!({
            "input": gradient_source(),
            "target": plain(vec![rect(0.0, 0.0, 20.0, 10.0), rect(20.0, 0.0, 40.0, 10.0)]),
            "field": "pop",
        }));
        let se = col(&out, "PRED_STD_ERR");
        assert!(se
            .iter()
            .all(|v| v.is_some_and(|x| x >= 0.0 && x.is_finite())));
    }

    #[test]
    fn a_target_matching_a_source_polygon_recovers_its_value() {
        // With identical geographies the answer must be the input.
        let (out, _) = run(json!({
            "input": gradient_source(),
            "target": plain((0..4).map(|i| {
                let x = i as f64 * 10.0;
                rect(x, 0.0, x + 10.0, 10.0)
            }).collect::<Vec<_>>()),
            "field": "pop",
        }));
        let p = col(&out, "PREDICTED");
        for (i, want) in [100.0, 150.0, 200.0, 250.0].iter().enumerate() {
            let got = p[i].unwrap();
            assert!(
                (got - want).abs() < want * 0.25,
                "cell {i}: got {got}, expected near {want}"
            );
        }
    }

    #[test]
    fn average_field_type_is_intensive_and_not_rescaled() {
        // An average must not be inflated to match a source "total".
        let (out, res) = run(json!({
            "input": valued(vec![
                (rect(0.0, 0.0, 10.0, 10.0), 20.0),
                (rect(10.0, 0.0, 20.0, 10.0), 20.0),
                (rect(20.0, 0.0, 30.0, 10.0), 20.0),
            ]),
            "target": plain(vec![rect(0.0, 0.0, 30.0, 10.0)]),
            "field": "pop",
            "field_type": "average",
        }));
        assert_eq!(res.outputs["mass_scale"], json!(1.0));
        let p = col(&out, "PREDICTED")[0].unwrap();
        assert!((p - 20.0).abs() < 1e-3, "constant average became {p}");
    }

    #[test]
    fn counts_split_a_uniform_field_in_proportion_to_area() {
        // A constant density over three equal cells, reaggregated to one
        // double-width and one single-width target, must split 2:1.
        let (out, _) = run(json!({
            "input": valued(vec![
                (rect(0.0, 0.0, 10.0, 10.0), 90.0),
                (rect(10.0, 0.0, 20.0, 10.0), 90.0),
                (rect(20.0, 0.0, 30.0, 10.0), 90.0),
            ]),
            "target": plain(vec![rect(0.0, 0.0, 20.0, 10.0), rect(20.0, 0.0, 30.0, 10.0)]),
            "field": "pop",
        }));
        let p = col(&out, "PREDICTED");
        let (big, small) = (p[0].unwrap(), p[1].unwrap());
        assert!((big + small - 270.0).abs() < 1e-6, "total drifted");
        assert!(
            (big / small - 2.0).abs() < 0.2,
            "expected a 2:1 split, got {big} vs {small}"
        );
    }

    #[test]
    fn preserve_total_can_be_turned_off() {
        let (_, res) = run(json!({
            "input": gradient_source(),
            "target": plain(vec![rect(0.0, 0.0, 10.0, 10.0)]),
            "field": "pop",
            "preserve_total": false,
        }));
        assert_eq!(res.outputs["mass_scale"], json!(1.0));
    }

    #[test]
    fn the_run_is_deterministic() {
        // The discretisation is a lattice, not a random sample.
        let go = || {
            let (out, _) = run(json!({
                "input": gradient_source(),
                "target": plain(vec![rect(5.0, 2.0, 27.0, 9.0)]),
                "field": "pop",
            }));
            col(&out, "PREDICTED")[0].unwrap()
        };
        assert_eq!(go(), go());
    }

    #[test]
    fn a_non_polygon_target_is_left_unresolved_rather_than_guessed() {
        let mut l = Layer::new("tgt");
        l.geom_type = Some(GeometryType::Point);
        l.add_feature(Some(Geometry::Point(Coord::xy(5.0, 5.0))), &[])
            .unwrap();
        let id = memory_store::put_vector(l);
        let point_target = memory_store::make_vector_memory_path(&id);

        let (out, res) = run(json!({
            "input": gradient_source(),
            "target": point_target,
            "field": "pop",
        }));
        assert_eq!(res.outputs["resolved_targets"], json!(0));
        assert!(col(&out, "PREDICTED")[0].is_none());
    }

    #[test]
    fn rejects_bad_parameters() {
        let src = gradient_source();
        let tgt = plain(vec![rect(0.0, 0.0, 10.0, 10.0)]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ArealInterpolationTool.validate(&args).is_err()
        };
        assert!(bad(json!({"target": tgt, "field": "pop"})));
        assert!(bad(json!({"input": src, "field": "pop"})));
        assert!(bad(
            json!({"input": src, "target": tgt, "field": "pop", "field_type": "nope"})
        ));
        assert!(bad(
            json!({"input": src, "target": tgt, "field": "pop", "discretization": 0})
        ));
    }

    #[test]
    fn a_missing_field_is_reported() {
        let args: ToolArgs = serde_json::from_value(json!({
            "input": gradient_source(),
            "target": plain(vec![rect(0.0, 0.0, 10.0, 10.0)]),
            "field": "nope",
        }))
        .unwrap();
        assert!(ArealInterpolationTool.run(&args, &ctx()).is_err());
    }
}
