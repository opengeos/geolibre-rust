//! GeoLibre tool: residuals and RMSE for a transformation fitted from
//! displacement links.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Calculate Transformation Errors*
//! (Editing). GeoLibre ships an unusually complete conflation family —
//! `transform_features`, `rubbersheet_features`, `edgematch_features`,
//! `align_features`, `propagate_displacement`, `integrate`, `snap_tracks` —
//! but none of them tell you whether the control you fed them was any good.
//!
//! This closes that loop: fit the transformation, then report how far each
//! link's source point lands from its own target under the fitted model. A
//! single mistyped control point dominates the RMSE and shows up immediately
//! in `RESIDUAL`, so it can be dropped *before* the transformation is applied
//! rather than discovered by eye afterwards.
//!
//! Each displacement link is a two-or-more-vertex line: its first vertex is
//! the source (from) location, its last vertex the target (to) location.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CalculateTransformationErrorsTool;

impl Tool for CalculateTransformationErrorsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "calculate_transformation_errors",
            display_name: "Calculate Transformation Errors",
            summary: "Fit a similarity/affine/projective transformation from displacement links and report each link's residual error plus the overall RMSE, so bad control points can be found and removed before transforming. Like ArcGIS Calculate Transformation Errors.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Displacement link lines: first vertex is the source (from) point, last vertex the target (to) point.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output path for the link table with residual fields appended. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'affine' (default, 6 parameters), 'similarity' (4) or 'projective' (8).",
                    required: false,
                },
                ToolParamSpec {
                    name: "keep_geometry",
                    description: "Keep each link's geometry on the output (default true). Set false for a pure table.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_method(args)?;
        parse_optional_bool(args, "keep_geometry")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let method = parse_method(args)?;
        let keep_geom = parse_optional_bool(args, "keep_geometry")?.unwrap_or(true);

        let layer = load_input_layer(input)?;
        if layer.features.is_empty() {
            return Err(ToolError::Execution("input has no features".to_string()));
        }

        // Extract (source, target) pairs from the link geometries.
        let mut links: Vec<(usize, [f64; 2], [f64; 2])> = Vec::new();
        for (i, f) in layer.features.iter().enumerate() {
            if let Some(g) = &f.geometry {
                if let Some((a, b)) = link_endpoints(g) {
                    links.push((i, a, b));
                }
            }
        }
        let n = links.len();
        let min_links = method.min_links();
        if n < min_links {
            return Err(ToolError::Execution(format!(
                "'{}' transformation needs at least {min_links} displacement link(s), found {n}",
                method.name()
            )));
        }
        ctx.progress
            .info(&format!("fitting {} from {n} link(s)", method.name()));

        let coeffs = fit(method, &links)?;

        // Residuals: where each source lands under the fitted model vs its target.
        let mut residuals = vec![(0.0_f64, 0.0_f64, 0.0_f64); layer.features.len()];
        let mut sum_sq = 0.0;
        let mut max_res = 0.0_f64;
        for &(i, src, dst) in &links {
            let p = apply(method, &coeffs, src);
            let (dx, dy) = (p[0] - dst[0], p[1] - dst[1]);
            let r = dx.hypot(dy);
            residuals[i] = (dx, dy, r);
            sum_sq += r * r;
            max_res = max_res.max(r);
        }
        let rmse = (sum_sq / n as f64).sqrt();

        // Build output: original attributes + residual fields.
        let mut out = Layer::new("transformation_errors");
        if keep_geom {
            if let Some(gt) = layer.geom_type {
                out = out.with_geom_type(gt);
            }
        }
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for fd in layer.schema.fields() {
            out.add_field(fd.clone());
        }
        out.add_field(FieldDef::new("RESIDUAL_X", FieldType::Float));
        out.add_field(FieldDef::new("RESIDUAL_Y", FieldType::Float));
        out.add_field(FieldDef::new("RESIDUAL", FieldType::Float));
        out.add_field(FieldDef::new("ERR_SHARE", FieldType::Float));

        let names: Vec<String> = layer
            .schema
            .fields()
            .iter()
            .map(|f| f.name.clone())
            .collect();
        let mut outliers = 0usize;
        // A link is flagged when its residual exceeds 3x the RMSE — the usual
        // "this control point is wrong" threshold.
        let outlier_cut = 3.0 * rmse;
        for (i, feat) in layer.features.iter().enumerate() {
            let (dx, dy, r) = residuals[i];
            if rmse > 0.0 && r > outlier_cut {
                outliers += 1;
            }
            let share = if sum_sq > 0.0 { r * r / sum_sq } else { 0.0 };
            let mut attrs: Vec<(&str, FieldValue)> = names
                .iter()
                .enumerate()
                .map(|(fi, nm)| {
                    (
                        nm.as_str(),
                        feat.attributes.get(fi).cloned().unwrap_or(FieldValue::Null),
                    )
                })
                .collect();
            attrs.push(("RESIDUAL_X", FieldValue::Float(dx)));
            attrs.push(("RESIDUAL_Y", FieldValue::Float(dy)));
            attrs.push(("RESIDUAL", FieldValue::Float(r)));
            attrs.push(("ERR_SHARE", FieldValue::Float(share)));
            let geom = if keep_geom { feat.geometry.clone() } else { None };
            out.add_feature(geom, &attrs)
                .map_err(|e| ToolError::Execution(format!("failed adding feature: {e}")))?;
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("method".to_string(), json!(method.name()));
        outputs.insert("link_count".to_string(), json!(n));
        outputs.insert("rmse".to_string(), json!(rmse));
        outputs.insert("max_residual".to_string(), json!(max_res));
        outputs.insert("outlier_count".to_string(), json!(outliers));
        outputs.insert("coefficients".to_string(), json!(coeffs));
        Ok(ToolRunResult { outputs })
    }
}

// ── Transformations ─────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Similarity,
    Affine,
    Projective,
}

impl Method {
    fn name(self) -> &'static str {
        match self {
            Method::Similarity => "similarity",
            Method::Affine => "affine",
            Method::Projective => "projective",
        }
    }
    /// Degrees of freedom / 2, i.e. the minimum control-point count.
    fn min_links(self) -> usize {
        match self {
            Method::Similarity => 2,
            Method::Affine => 3,
            Method::Projective => 4,
        }
    }
}

/// Least-squares fit of the chosen model to the (source -> target) pairs.
fn fit(method: Method, links: &[(usize, [f64; 2], [f64; 2])]) -> Result<Vec<f64>, ToolError> {
    match method {
        // x' =  a*x - b*y + tx
        // y' =  b*x + a*y + ty        (4 unknowns: a, b, tx, ty)
        Method::Similarity => {
            let mut ata = vec![0.0; 4 * 4];
            let mut atb = vec![0.0; 4];
            for &(_, s, d) in links {
                let (x, y) = (s[0], s[1]);
                // Row for x': [x, -y, 1, 0]
                accumulate(&mut ata, &mut atb, &[x, -y, 1.0, 0.0], d[0]);
                // Row for y': [y,  x, 0, 1]
                accumulate(&mut ata, &mut atb, &[y, x, 0.0, 1.0], d[1]);
            }
            solve(&mut ata, &mut atb, 4)
        }
        // Two independent 3-parameter solves (x and y are decoupled in affine).
        Method::Affine => {
            let mut ata = vec![0.0; 3 * 3];
            let mut atbx = vec![0.0; 3];
            let mut atby = vec![0.0; 3];
            for &(_, s, d) in links {
                let row = [s[0], s[1], 1.0];
                // Shared design matrix, two right-hand sides.
                for (i, ri) in row.iter().enumerate() {
                    for (j, rj) in row.iter().enumerate() {
                        ata[i * 3 + j] += ri * rj;
                    }
                    atbx[i] += ri * d[0];
                    atby[i] += ri * d[1];
                }
            }
            let mut ata2 = ata.clone();
            let cx = solve(&mut ata, &mut atbx, 3)?;
            let cy = solve(&mut ata2, &mut atby, 3)?;
            Ok(vec![cx[0], cx[1], cx[2], cy[0], cy[1], cy[2]])
        }
        // Direct linear transform, normalising h_33 = 1:
        //   x' = (a*x + b*y + c) / (g*x + h*y + 1)
        //   y' = (d*x + e*y + f) / (g*x + h*y + 1)
        // rearranged to linear form in [a,b,c,d,e,f,g,h].
        Method::Projective => {
            let mut ata = vec![0.0; 8 * 8];
            let mut atb = vec![0.0; 8];
            for &(_, s, d) in links {
                let (x, y) = (s[0], s[1]);
                let (xp, yp) = (d[0], d[1]);
                accumulate(
                    &mut ata,
                    &mut atb,
                    &[x, y, 1.0, 0.0, 0.0, 0.0, -x * xp, -y * xp],
                    xp,
                );
                accumulate(
                    &mut ata,
                    &mut atb,
                    &[0.0, 0.0, 0.0, x, y, 1.0, -x * yp, -y * yp],
                    yp,
                );
            }
            solve(&mut ata, &mut atb, 8)
        }
    }
}

/// Adds one observation row to the normal equations `AtA p = Atb`.
fn accumulate(ata: &mut [f64], atb: &mut [f64], row: &[f64], rhs: f64) {
    let n = row.len();
    for i in 0..n {
        for j in 0..n {
            ata[i * n + j] += row[i] * row[j];
        }
        atb[i] += row[i] * rhs;
    }
}

/// Gaussian elimination with partial pivoting on an `n x n` system.
fn solve(a: &mut [f64], b: &mut [f64], n: usize) -> Result<Vec<f64>, ToolError> {
    for col in 0..n {
        let mut piv = col;
        for r in (col + 1)..n {
            if a[r * n + col].abs() > a[piv * n + col].abs() {
                piv = r;
            }
        }
        if a[piv * n + col].abs() < 1e-12 {
            return Err(ToolError::Execution(
                "displacement links are degenerate (collinear or coincident); \
                 the transformation is not determined"
                    .to_string(),
            ));
        }
        if piv != col {
            for c in 0..n {
                a.swap(col * n + c, piv * n + c);
            }
            b.swap(col, piv);
        }
        let d = a[col * n + col];
        for r in (col + 1)..n {
            let factor = a[r * n + col] / d;
            if factor == 0.0 {
                continue;
            }
            for c in col..n {
                a[r * n + c] -= factor * a[col * n + c];
            }
            b[r] -= factor * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let mut sum = b[row];
        for c in (row + 1)..n {
            sum -= a[row * n + c] * x[c];
        }
        x[row] = sum / a[row * n + row];
    }
    Ok(x)
}

/// Applies a fitted model to a source point.
fn apply(method: Method, c: &[f64], p: [f64; 2]) -> [f64; 2] {
    let (x, y) = (p[0], p[1]);
    match method {
        Method::Similarity => [
            c[0] * x - c[1] * y + c[2],
            c[1] * x + c[0] * y + c[3],
        ],
        Method::Affine => [
            c[0] * x + c[1] * y + c[2],
            c[3] * x + c[4] * y + c[5],
        ],
        Method::Projective => {
            let w = c[6] * x + c[7] * y + 1.0;
            // A degenerate denominator means the point maps to infinity; return
            // it unmoved so the residual reports the failure rather than NaN.
            if w.abs() < 1e-12 {
                return p;
            }
            [
                (c[0] * x + c[1] * y + c[2]) / w,
                (c[3] * x + c[4] * y + c[5]) / w,
            ]
        }
    }
}

/// First and last vertex of a link geometry.
fn link_endpoints(g: &Geometry) -> Option<([f64; 2], [f64; 2])> {
    fn ends(cs: &[Coord]) -> Option<([f64; 2], [f64; 2])> {
        if cs.len() < 2 {
            return None;
        }
        let a = cs.first()?;
        let b = cs.last()?;
        Some(([a.x, a.y], [b.x, b.y]))
    }
    match g {
        Geometry::LineString(cs) => ends(cs),
        Geometry::MultiLineString(ls) => ls.iter().find_map(|l| ends(l)),
        Geometry::MultiPoint(cs) => ends(cs),
        _ => None,
    }
}

// ── Params ──────────────────────────────────────────────────────────────────

fn parse_method(args: &ToolArgs) -> Result<Method, ToolError> {
    match args
        .get("method")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        None | Some("") | Some("affine") => Ok(Method::Affine),
        Some("similarity") => Ok(Method::Similarity),
        Some("projective") => Ok(Method::Projective),
        Some(o) => Err(ToolError::Validation(format!(
            "'method' must be 'affine', 'similarity' or 'projective', got '{o}'"
        ))),
    }
}

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
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
    use wbvector::GeometryType;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds link lines from (from_xy, to_xy) pairs.
    fn links(pairs: &[([f64; 2], [f64; 2])]) -> String {
        let mut l = Layer::new("links")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        for (a, b) in pairs {
            l.add_feature(
                Some(Geometry::line_string(vec![
                    Coord::xy(a[0], a[1]),
                    Coord::xy(b[0], b[1]),
                ])),
                &[],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CalculateTransformationErrorsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    #[test]
    fn exact_affine_fit_has_zero_residual() {
        // Targets generated by x' = 2x + 10, y' = 2y - 5 exactly.
        let f = |p: [f64; 2]| [2.0 * p[0] + 10.0, 2.0 * p[1] - 5.0];
        let src = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [3.0, 4.0]];
        let pairs: Vec<_> = src.iter().map(|&s| (s, f(s))).collect();
        let (out, _l) = run(json!({ "input": links(&pairs), "method": "affine" }));
        assert!(out.outputs["rmse"].as_f64().unwrap() < 1e-9);
        assert_eq!(out.outputs["outlier_count"], json!(0));
    }

    #[test]
    fn one_bad_link_dominates_the_error_share() {
        // Three clean links on x' = x + 100, plus one control point typo.
        let mut pairs: Vec<([f64; 2], [f64; 2])> = vec![
            ([0.0, 0.0], [100.0, 0.0]),
            ([10.0, 0.0], [110.0, 0.0]),
            ([0.0, 10.0], [100.0, 10.0]),
            ([10.0, 10.0], [110.0, 10.0]),
        ];
        pairs.push(([5.0, 5.0], [900.0, 900.0])); // typo
        let (out, layer) = run(json!({ "input": links(&pairs), "method": "affine" }));
        let share = layer.schema.field_index("ERR_SHARE").unwrap();
        let shares: Vec<f64> = layer
            .iter()
            .map(|f| f.attributes[share].as_f64().unwrap())
            .collect();
        // The bad link is the last one and must carry most of the squared error.
        let worst = shares
            .iter()
            .cloned()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .unwrap();
        assert_eq!(worst.0, 4);
        assert!(worst.1 > 0.5, "bad link share was {}", worst.1);
        assert!(out.outputs["rmse"].as_f64().unwrap() > 1.0);
    }

    #[test]
    fn similarity_recovers_rotation_and_scale() {
        // 90-degree rotation about the origin, scale 2.
        let f = |p: [f64; 2]| [-2.0 * p[1], 2.0 * p[0]];
        let src = [[1.0, 0.0], [0.0, 1.0], [2.0, 3.0]];
        let pairs: Vec<_> = src.iter().map(|&s| (s, f(s))).collect();
        let (out, _l) = run(json!({ "input": links(&pairs), "method": "similarity" }));
        assert!(out.outputs["rmse"].as_f64().unwrap() < 1e-9);
        let c = out.outputs["coefficients"].as_array().unwrap();
        // a = scale*cos(theta) = 0, b = scale*sin(theta) = 2
        assert!((c[0].as_f64().unwrap() - 0.0).abs() < 1e-9);
        assert!((c[1].as_f64().unwrap() - 2.0).abs() < 1e-9);
    }

    #[test]
    fn projective_fits_a_perspective_warp() {
        // Unit square -> irregular quad: affine cannot fit this, projective can.
        let pairs = [
            ([0.0, 0.0], [0.0, 0.0]),
            ([1.0, 0.0], [1.0, 0.0]),
            ([1.0, 1.0], [0.8, 1.2]),
            ([0.0, 1.0], [0.1, 1.0]),
        ];
        let (proj, _l) = run(json!({ "input": links(&pairs), "method": "projective" }));
        let (aff, _l) = run(json!({ "input": links(&pairs), "method": "affine" }));
        assert!(proj.outputs["rmse"].as_f64().unwrap() < 1e-9);
        assert!(aff.outputs["rmse"].as_f64().unwrap() > 1e-3);
    }

    #[test]
    fn too_few_links_is_rejected() {
        let pairs = [([0.0, 0.0], [1.0, 1.0]), ([1.0, 0.0], [2.0, 1.0])];
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": links(&pairs), "method": "projective" }))
                .unwrap();
        assert!(CalculateTransformationErrorsTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn collinear_links_are_reported_as_degenerate() {
        // All sources on one line: an affine fit is underdetermined.
        let pairs = [
            ([0.0, 0.0], [0.0, 0.0]),
            ([1.0, 0.0], [1.0, 0.0]),
            ([2.0, 0.0], [2.0, 0.0]),
        ];
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": links(&pairs), "method": "affine" })).unwrap();
        assert!(CalculateTransformationErrorsTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CalculateTransformationErrorsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "l.shp", "method": "helmert" })).is_err());
        assert!(bad(json!({ "input": "l.shp" })).is_ok());
    }
}
