//! GeoLibre tool: descriptive geographic distribution statistics.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Measuring Geographic Distributions*
//! toolset, exposed as one tool with a `statistic` selector:
//!
//! - `mean_center` — the (optionally weighted) average location.
//! - `median_center` — the point minimizing total distance to all features,
//!   found by Weiszfeld iteration (robust to outliers).
//! - `central_feature` — the input feature with the smallest total distance to
//!   all others (an actual location, not a synthetic one).
//! - `standard_distance` — the mean center plus a circle whose radius is the
//!   standard distance (root-mean-square spread), times `n_std`.
//! - `standard_deviational_ellipse` — the directional-distribution ellipse from
//!   the coordinate covariance matrix (its eigenvectors give the axis
//!   directions, the square-root eigenvalues the axis lengths), times `n_std`.
//!
//! Whitebox-wasm covers the inferential spatial statistics (Moran's I, Getis-Ord
//! Gi*, nearest-neighbour index) but none of these descriptive measures.
//!
//! A `weight_field` weights each feature (e.g. population); a `case_field`
//! partitions the features and emits one result per group. Point and multipoint
//! inputs use their coordinates directly; other geometries use a representative
//! point (the mean of their vertices). Output attributes carry the numbers a
//! caller would report (centre coordinates, standard distance, ellipse rotation
//! and axis lengths).
//!
//! Scope for v1: ArcGIS's *Linear Directional Mean* (circular statistics over
//! line orientation) is not included — this tool covers the point-distribution
//! measures.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Vertices used to approximate the standard-distance circle and the ellipse.
const ELLIPSE_SEGMENTS: usize = 120;

pub struct DirectionalDistributionTool;

impl Tool for DirectionalDistributionTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "directional_distribution",
            display_name: "Directional Distribution",
            summary: "Descriptive geographic distribution statistics: mean center, median center, central feature, standard distance (circle), and the standard deviational ellipse — with optional weighting and grouping.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (points preferred; other geometries use their vertex-mean representative point).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "statistic",
                    description: "Which measure to compute: 'mean_center', 'median_center', 'central_feature', 'standard_distance', or 'standard_deviational_ellipse'.",
                    required: true,
                },
                ToolParamSpec {
                    name: "weight_field",
                    description: "Optional numeric field to weight each feature (e.g. population). Default: every feature weighted 1.",
                    required: false,
                },
                ToolParamSpec {
                    name: "case_field",
                    description: "Optional field to group features by; one output feature is produced per distinct value.",
                    required: false,
                },
                ToolParamSpec {
                    name: "n_std",
                    description: "Standard-deviation multiplier (1, 2, or 3) for the standard-distance circle and the deviational ellipse. Default 1.",
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

        let layer = load_input_layer(input)?;
        let layer_crs = layer.crs.clone();
        let schema = &layer.schema;

        // Collect (group, observation) pairs.
        let mut groups: BTreeMap<String, Vec<Obs>> = BTreeMap::new();
        for feature in &layer.features {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let Some((x, y)) = rep_point(geom) else {
                continue;
            };
            let w = match &prm.weight_field {
                Some(f) => match feature.get(schema, f).ok().and_then(FieldValue::as_f64) {
                    Some(v) if v.is_finite() && v > 0.0 => v,
                    _ => continue, // skip features with a missing / non-positive weight
                },
                None => 1.0,
            };
            let key = match &prm.case_field {
                Some(f) => feature
                    .get(schema, f)
                    .map(field_value_string)
                    .unwrap_or_default(),
                None => String::new(),
            };
            groups.entry(key).or_default().push(Obs { x, y, w });
        }

        ctx.progress.info(&format!(
            "{} feature(s): computing {} over {} group(s)",
            layer.len(),
            prm.statistic.as_str(),
            groups.len().max(1)
        ));

        let mut out_layer = Layer::new(layer.name.clone());
        out_layer.crs = layer_crs;
        let grouped = prm.case_field.is_some();
        if grouped {
            out_layer.add_field(FieldDef::new("group", FieldType::Text));
        }
        for &name in prm.statistic.output_fields() {
            out_layer.add_field(FieldDef::new(name, FieldType::Float));
        }
        out_layer.geom_type = Some(prm.statistic.geom_type());

        let mut produced = 0usize;
        for (key, obs) in &groups {
            let Some((geom, mut fields)) = prm.statistic.compute(obs, prm.n_std) else {
                continue;
            };
            if grouped {
                fields.insert(0, ("group", FieldValue::Text(key.clone())));
            }
            out_layer
                .add_feature(Some(geom), &fields)
                .map_err(|e| ToolError::Execution(format!("failed writing output feature: {e}")))?;
            produced += 1;
        }

        let feature_count = out_layer.len();
        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("statistic".to_string(), json!(prm.statistic.as_str()));
        outputs.insert("group_count".to_string(), json!(produced));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        Ok(ToolRunResult { outputs })
    }
}

// ── Observations & statistics ─────────────────────────────────────────────────

struct Obs {
    x: f64,
    y: f64,
    w: f64,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Statistic {
    MeanCenter,
    MedianCenter,
    CentralFeature,
    StandardDistance,
    StandardDeviationalEllipse,
}

impl Statistic {
    fn as_str(self) -> &'static str {
        match self {
            Self::MeanCenter => "mean_center",
            Self::MedianCenter => "median_center",
            Self::CentralFeature => "central_feature",
            Self::StandardDistance => "standard_distance",
            Self::StandardDeviationalEllipse => "standard_deviational_ellipse",
        }
    }

    fn geom_type(self) -> GeometryType {
        match self {
            Self::MeanCenter | Self::MedianCenter | Self::CentralFeature => GeometryType::Point,
            Self::StandardDistance | Self::StandardDeviationalEllipse => GeometryType::Polygon,
        }
    }

    /// Output attribute names (all Float) for this statistic, in order.
    fn output_fields(self) -> &'static [&'static str] {
        match self {
            Self::MeanCenter | Self::MedianCenter | Self::CentralFeature => &["x", "y"],
            Self::StandardDistance => &["center_x", "center_y", "std_dist", "n_std", "radius"],
            Self::StandardDeviationalEllipse => &[
                "center_x",
                "center_y",
                "rotation",
                "semi_major",
                "semi_minor",
                "n_std",
            ],
        }
    }

    /// Computes the statistic for one group, returning its geometry and the
    /// float attribute values (matching `output_fields`). `None` if the group
    /// is degenerate (e.g. empty).
    fn compute(
        self,
        obs: &[Obs],
        n_std: f64,
    ) -> Option<(Geometry, Vec<(&'static str, FieldValue)>)> {
        if obs.is_empty() {
            return None;
        }
        let sw: f64 = obs.iter().map(|o| o.w).sum();
        if sw <= 0.0 {
            return None;
        }
        let (mx, my) = (
            obs.iter().map(|o| o.w * o.x).sum::<f64>() / sw,
            obs.iter().map(|o| o.w * o.y).sum::<f64>() / sw,
        );
        match self {
            Self::MeanCenter => Some((Geometry::point(mx, my), vec![("x", f(mx)), ("y", f(my))])),
            Self::MedianCenter => {
                let (cx, cy) = weiszfeld(obs, mx, my);
                Some((Geometry::point(cx, cy), vec![("x", f(cx)), ("y", f(cy))]))
            }
            Self::CentralFeature => {
                let (cx, cy) = central_feature(obs);
                Some((Geometry::point(cx, cy), vec![("x", f(cx)), ("y", f(cy))]))
            }
            Self::StandardDistance => {
                let var = obs
                    .iter()
                    .map(|o| o.w * ((o.x - mx).powi(2) + (o.y - my).powi(2)))
                    .sum::<f64>()
                    / sw;
                let sd = var.sqrt();
                let radius = sd * n_std;
                let geom = circle(mx, my, radius);
                Some((
                    geom,
                    vec![
                        ("center_x", f(mx)),
                        ("center_y", f(my)),
                        ("std_dist", f(sd)),
                        ("n_std", f(n_std)),
                        ("radius", f(radius)),
                    ],
                ))
            }
            Self::StandardDeviationalEllipse => {
                let (a, b, theta) = ellipse_axes(obs, mx, my)?;
                let (semi_major, semi_minor) = (a * n_std, b * n_std);
                let geom = ellipse(mx, my, semi_major, semi_minor, theta);
                Some((
                    geom,
                    vec![
                        ("center_x", f(mx)),
                        ("center_y", f(my)),
                        ("rotation", f(theta.to_degrees())),
                        ("semi_major", f(semi_major)),
                        ("semi_minor", f(semi_minor)),
                        ("n_std", f(n_std)),
                    ],
                ))
            }
        }
    }
}

fn f(v: f64) -> FieldValue {
    FieldValue::Float(v)
}

/// Weiszfeld iteration for the weighted geometric median, seeded at the mean.
fn weiszfeld(obs: &[Obs], mut x: f64, mut y: f64) -> (f64, f64) {
    for _ in 0..256 {
        let (mut num_x, mut num_y, mut den) = (0.0, 0.0, 0.0);
        let mut coincident = false;
        for o in obs {
            let d = ((o.x - x).powi(2) + (o.y - y).powi(2)).sqrt();
            if d < 1e-12 {
                coincident = true;
                break;
            }
            let wd = o.w / d;
            num_x += wd * o.x;
            num_y += wd * o.y;
            den += wd;
        }
        if coincident || den <= 0.0 {
            break;
        }
        let (nx, ny) = (num_x / den, num_y / den);
        let step = ((nx - x).powi(2) + (ny - y).powi(2)).sqrt();
        x = nx;
        y = ny;
        if step < 1e-9 {
            break;
        }
    }
    (x, y)
}

/// The observation minimizing the weighted sum of distances to all others.
fn central_feature(obs: &[Obs]) -> (f64, f64) {
    let mut best = (obs[0].x, obs[0].y);
    let mut best_cost = f64::INFINITY;
    for a in obs {
        let cost: f64 = obs
            .iter()
            .map(|b| b.w * ((a.x - b.x).powi(2) + (a.y - b.y).powi(2)).sqrt())
            .sum();
        if cost < best_cost {
            best_cost = cost;
            best = (a.x, a.y);
        }
    }
    best
}

/// Semi-major, semi-minor, and rotation (radians, CCW from +x) of the standard
/// deviational ellipse: the eigen-decomposition of the weighted 2x2 coordinate
/// covariance matrix. `None` if degenerate (all points coincident).
fn ellipse_axes(obs: &[Obs], mx: f64, my: f64) -> Option<(f64, f64, f64)> {
    let sw: f64 = obs.iter().map(|o| o.w).sum();
    let mut sxx = 0.0;
    let mut syy = 0.0;
    let mut sxy = 0.0;
    for o in obs {
        let dx = o.x - mx;
        let dy = o.y - my;
        sxx += o.w * dx * dx;
        syy += o.w * dy * dy;
        sxy += o.w * dx * dy;
    }
    sxx /= sw;
    syy /= sw;
    sxy /= sw;

    let trace = sxx + syy;
    let det = sxx * syy - sxy * sxy;
    let disc = ((trace * 0.5).powi(2) - det).max(0.0).sqrt();
    let l1 = trace * 0.5 + disc; // larger eigenvalue -> major axis
    let l2 = (trace * 0.5 - disc).max(0.0);
    if l1 <= 0.0 {
        return None;
    }
    // Eigenvector for l1.
    let theta = if sxy.abs() > 1e-12 {
        (l1 - sxx).atan2(sxy)
    } else if sxx >= syy {
        0.0
    } else {
        std::f64::consts::FRAC_PI_2
    };
    Some((l1.sqrt(), l2.sqrt(), theta))
}

// ── Geometry builders ─────────────────────────────────────────────────────────

fn circle(cx: f64, cy: f64, r: f64) -> Geometry {
    ellipse(cx, cy, r, r, 0.0)
}

/// A closed ellipse polygon centered at `(cx, cy)` with the given semi-axes,
/// rotated `theta` radians CCW.
fn ellipse(cx: f64, cy: f64, a: f64, b: f64, theta: f64) -> Geometry {
    let (ct, st) = (theta.cos(), theta.sin());
    let coords: Vec<Coord> = (0..ELLIPSE_SEGMENTS)
        .map(|i| {
            let phi = 2.0 * std::f64::consts::PI * i as f64 / ELLIPSE_SEGMENTS as f64;
            let (ex, ey) = (a * phi.cos(), b * phi.sin());
            Coord::xy(cx + ex * ct - ey * st, cy + ex * st + ey * ct)
        })
        .collect();
    Geometry::Polygon {
        exterior: Ring::new(coords),
        interiors: vec![],
    }
}

/// Representative point of a geometry: the coordinate for a point, otherwise the
/// mean of the geometry's vertices.
fn rep_point(geom: &Geometry) -> Option<(f64, f64)> {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0usize;
    collect_coords(geom, &mut |x, y| {
        sx += x;
        sy += y;
        n += 1;
    });
    (n > 0).then(|| (sx / n as f64, sy / n as f64))
}

fn collect_coords(geom: &Geometry, f: &mut impl FnMut(f64, f64)) {
    match geom {
        Geometry::Point(c) => f(c.x, c.y),
        Geometry::MultiPoint(cs) | Geometry::LineString(cs) => cs.iter().for_each(|c| f(c.x, c.y)),
        Geometry::MultiLineString(ls) => ls.iter().flatten().for_each(|c| f(c.x, c.y)),
        Geometry::Polygon { exterior, .. } => exterior.coords().iter().for_each(|c| f(c.x, c.y)),
        Geometry::MultiPolygon(parts) => parts
            .iter()
            .flat_map(|(e, _)| e.coords())
            .for_each(|c| f(c.x, c.y)),
        Geometry::GeometryCollection(gs) => gs.iter().for_each(|g| collect_coords(g, f)),
    }
}

fn field_value_string(v: &FieldValue) -> String {
    match v {
        FieldValue::Text(s) | FieldValue::Date(s) | FieldValue::DateTime(s) => s.clone(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(x) => x.to_string(),
        FieldValue::Boolean(b) => b.to_string(),
        FieldValue::Null | FieldValue::Blob(_) => String::new(),
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    statistic: Statistic,
    weight_field: Option<String>,
    case_field: Option<String>,
    n_std: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let statistic = match parse_optional_str(args, "statistic")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        Some("mean_center") => Statistic::MeanCenter,
        Some("median_center") => Statistic::MedianCenter,
        Some("central_feature") => Statistic::CentralFeature,
        Some("standard_distance") => Statistic::StandardDistance,
        Some("standard_deviational_ellipse") | Some("ellipse") => {
            Statistic::StandardDeviationalEllipse
        }
        None => {
            return Err(ToolError::Validation(
                "missing required parameter 'statistic'".to_string(),
            ))
        }
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown statistic '{other}' (expected mean_center, median_center, central_feature, standard_distance, or standard_deviational_ellipse)"
            )))
        }
    };
    let weight_field = parse_optional_str(args, "weight_field")?.map(str::to_string);
    let case_field = parse_optional_str(args, "case_field")?.map(str::to_string);
    let n_std = match parse_optional_f64(args, "n_std")? {
        None => 1.0,
        Some(v) if v.fract() == 0.0 && (1.0..=3.0).contains(&v) => v,
        Some(_) => {
            return Err(ToolError::Validation(
                "parameter 'n_std' must be 1, 2, or 3".to_string(),
            ))
        }
    };
    Ok(Params {
        statistic,
        weight_field,
        case_field,
        n_std,
    })
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
    use wbvector::{memory_store, FieldDef, FieldType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn points_layer(pts: &[(f64, f64)]) -> String {
        let mut layer = Layer::new("pts");
        for &(x, y) in pts {
            layer.add_feature(Some(Geometry::point(x, y)), &[]).unwrap();
        }
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = DirectionalDistributionTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn fv(layer: &Layer, idx: usize, name: &str) -> f64 {
        layer.features[idx]
            .get(&layer.schema, name)
            .ok()
            .and_then(FieldValue::as_f64)
            .unwrap()
    }

    #[test]
    fn mean_center_of_a_symmetric_cross() {
        let input = points_layer(&[(-2.0, 0.0), (2.0, 0.0), (0.0, -2.0), (0.0, 2.0)]);
        let (_, layer) = run_tool(json!({ "input": input, "statistic": "mean_center" }));
        assert_eq!(layer.len(), 1);
        assert!(fv(&layer, 0, "x").abs() < 1e-9 && fv(&layer, 0, "y").abs() < 1e-9);
        match layer.features[0].geometry.as_ref().unwrap() {
            Geometry::Point(c) => assert!(c.x.abs() < 1e-9 && c.y.abs() < 1e-9),
            other => panic!("expected point, got {other:?}"),
        }
    }

    #[test]
    fn weighted_mean_center_shifts_toward_weight() {
        let mut layer = Layer::new("pts");
        layer.add_field(FieldDef::new("pop", FieldType::Float));
        layer
            .add_feature(
                Some(Geometry::point(0.0, 0.0)),
                &[("pop", FieldValue::Float(1.0))],
            )
            .unwrap();
        layer
            .add_feature(
                Some(Geometry::point(10.0, 0.0)),
                &[("pop", FieldValue::Float(3.0))],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);
        let (_, layer) =
            run_tool(json!({ "input": input, "statistic": "mean_center", "weight_field": "pop" }));
        // Weighted mean x = (0*1 + 10*3)/4 = 7.5.
        assert!((fv(&layer, 0, "x") - 7.5).abs() < 1e-9);
    }

    #[test]
    fn standard_distance_circle_radius() {
        // Four points at distance 2 from the origin: RMS distance = 2.
        let input = points_layer(&[(-2.0, 0.0), (2.0, 0.0), (0.0, -2.0), (0.0, 2.0)]);
        let (_, layer) = run_tool(json!({ "input": input, "statistic": "standard_distance" }));
        assert!((fv(&layer, 0, "std_dist") - 2.0).abs() < 1e-9);
        assert!((fv(&layer, 0, "radius") - 2.0).abs() < 1e-9);
        // n_std=2 doubles the radius.
        let (_, layer2) = run_tool(
            json!({ "input": points_layer(&[(-2.0, 0.0), (2.0, 0.0), (0.0, -2.0), (0.0, 2.0)]), "statistic": "standard_distance", "n_std": 2 }),
        );
        assert!((fv(&layer2, 0, "radius") - 4.0).abs() < 1e-9);
    }

    #[test]
    fn ellipse_aligns_with_the_spread_axis() {
        // Points spread along x far more than y -> major axis ~horizontal
        // (rotation near 0 degrees) and semi_major > semi_minor.
        let pts: Vec<(f64, f64)> = (-5..=5)
            .map(|i| (i as f64 * 4.0, if i % 2 == 0 { 1.0 } else { -1.0 }))
            .collect();
        let input = points_layer(&pts);
        let (_, layer) =
            run_tool(json!({ "input": input, "statistic": "standard_deviational_ellipse" }));
        let (major, minor, rot) = (
            fv(&layer, 0, "semi_major"),
            fv(&layer, 0, "semi_minor"),
            fv(&layer, 0, "rotation"),
        );
        assert!(
            major > minor * 3.0,
            "major {major} should dominate minor {minor}"
        );
        let rot_norm = rot.rem_euclid(180.0);
        assert!(
            !(5.0..=175.0).contains(&rot_norm),
            "major axis should be ~horizontal, rotation was {rot}"
        );
        // Output is a closed polygon.
        assert!(matches!(
            layer.features[0].geometry.as_ref().unwrap(),
            Geometry::Polygon { .. }
        ));
    }

    #[test]
    fn central_feature_is_an_input_point() {
        // A tight cluster plus one outlier; the central feature is in the cluster.
        let input = points_layer(&[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (100.0, 100.0)]);
        let (_, layer) = run_tool(json!({ "input": input, "statistic": "central_feature" }));
        let (x, y) = (fv(&layer, 0, "x"), fv(&layer, 0, "y"));
        assert!(
            x < 50.0 && y < 50.0,
            "central feature ({x},{y}) should be in the cluster"
        );
        // It coincides with one of the inputs.
        let inputs = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (100.0, 100.0)];
        assert!(inputs
            .iter()
            .any(|&(ix, iy)| (ix - x).abs() < 1e-9 && (iy - y).abs() < 1e-9));
    }

    #[test]
    fn case_field_produces_one_feature_per_group() {
        let mut layer = Layer::new("pts");
        layer.add_field(FieldDef::new("region", FieldType::Text));
        for (x, y, r) in [
            (0.0, 0.0, "A"),
            (2.0, 0.0, "A"),
            (100.0, 0.0, "B"),
            (102.0, 0.0, "B"),
        ] {
            layer
                .add_feature(
                    Some(Geometry::point(x, y)),
                    &[("region", FieldValue::Text(r.into()))],
                )
                .unwrap();
        }
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);
        let (out, layer) =
            run_tool(json!({ "input": input, "statistic": "mean_center", "case_field": "region" }));
        assert_eq!(out.outputs["group_count"], json!(2));
        assert_eq!(layer.len(), 2);
        // Group A mean x = 1, group B mean x = 101.
        let xs: Vec<f64> = (0..2).map(|i| fv(&layer, i, "x")).collect();
        assert!(xs.iter().any(|&x| (x - 1.0).abs() < 1e-9));
        assert!(xs.iter().any(|&x| (x - 101.0).abs() < 1e-9));
    }

    #[test]
    fn median_center_resists_an_outlier() {
        // Cluster near origin plus a far outlier: the median center stays near
        // the cluster, unlike the mean.
        let input = points_layer(&[
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 1.0),
            (1.0, 1.0),
            (100.0, 100.0),
        ]);
        let (_, med) = run_tool(json!({ "input": input.clone(), "statistic": "median_center" }));
        let (_, mean) = run_tool(json!({ "input": input, "statistic": "mean_center" }));
        assert!(
            fv(&med, 0, "x") < fv(&mean, 0, "x"),
            "median should resist the outlier"
        );
        assert!(fv(&med, 0, "x") < 10.0 && fv(&med, 0, "y") < 10.0);
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = DirectionalDistributionTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "input": "x.geojson" })).is_err(),
            "missing statistic"
        );
        assert!(bad(json!({ "input": "x.geojson", "statistic": "bogus" })).is_err());
        assert!(
            bad(json!({ "input": "x.geojson", "statistic": "mean_center", "n_std": 4 })).is_err()
        );
        assert!(
            bad(json!({ "input": "x.geojson", "statistic": "mean_center", "n_std": 1.5 })).is_err()
        );
        assert!(bad(json!({ "input": "x.geojson", "statistic": "mean_center" })).is_ok());
    }
}
