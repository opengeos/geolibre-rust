//! GeoLibre tool: incremental spatial autocorrelation (Moran's I vs distance).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Incremental Spatial Autocorrelation*
//! (Spatial Statistics). The bundled `global_morans_i` runs at one distance
//! conceptualization; users picking a neighbourhood band for
//! `getis_ord_gi_star` otherwise guess. This computes global Moran's I at a
//! series of increasing fixed-distance bands and reports the z-score curve plus
//! its **first** and **maximum** peaks — the standard way to choose the scale of
//! clustering.
//!
//! At each band distance `d` the spatial weight `w_ij` is 1 when two features
//! are within `d` (binary, symmetric, no self-weight). Moran's I, its expected
//! value `-1/(n-1)`, and the **randomization** variance (Esri's default, using
//! S0/S1/S2 and the data kurtosis) give a z-score and p-value per band. Peaks in
//! the z-curve mark the distances at which the process is most clustered.
//!
//! `begin_distance`, `increment`, and `num_bands` default to a sweep anchored on
//! the average nearest-neighbour distance (so every feature has at least one
//! neighbour in the first band). Output is a table `distance, morans_i,
//! expected_i, variance, z_score, p_value`, written to CSV when `output` is
//! given and always returned in the result. Pairwise, so this suits moderate
//! feature counts. Use a projected CRS (distances in its units).

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldValue, Geometry};

use crate::common::write_text_output;
use crate::vector_common::{load_input_layer, parse_optional_str};

pub struct IncrementalSpatialAutocorrelationTool;

impl Tool for IncrementalSpatialAutocorrelationTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "incremental_spatial_autocorrelation",
            display_name: "Incremental Spatial Autocorrelation",
            summary: "Compute global Moran's I across a series of increasing distance bands and report the z-score curve with its first and maximum peaks — the defensible way to pick a clustering distance, like ArcGIS Incremental Spatial Autocorrelation.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer (other geometries use their vertex-mean representative point). Use a projected CRS.",
                    required: true,
                },
                ToolParamSpec {
                    name: "field",
                    description: "Numeric field to test for spatial autocorrelation.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional CSV path for the per-band table. Always returned in the result.",
                    required: false,
                },
                ToolParamSpec {
                    name: "begin_distance",
                    description: "First band distance, in CRS units. Default: the maximum nearest-neighbour distance (every feature has >=1 neighbour).",
                    required: false,
                },
                ToolParamSpec {
                    name: "increment",
                    description: "Distance added at each band, in CRS units. Default: the begin distance.",
                    required: false,
                },
                ToolParamSpec {
                    name: "num_bands",
                    description: "Number of distance bands to evaluate (default 10).",
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
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let fidx = layer
            .schema
            .field_index(field)
            .ok_or_else(|| ToolError::Validation(format!("field '{field}' not found")))?;

        // Collect points and values.
        let mut pts: Vec<(f64, f64)> = Vec::new();
        let mut vals: Vec<f64> = Vec::new();
        for feature in layer.iter() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let Some((x, y)) = representative_xy(geom) else {
                continue;
            };
            let Some(v) = feature.attributes.get(fidx).and_then(FieldValue::as_f64) else {
                continue;
            };
            pts.push((x, y));
            vals.push(v);
        }
        let n = pts.len();
        if n < 4 {
            return Err(ToolError::Execution(format!(
                "need at least 4 valued features, found {n}"
            )));
        }

        // Precompute the pairwise distance matrix (upper triangle).
        ctx.progress.info("computing pairwise distances");
        let mut nearest = vec![f64::INFINITY; n];
        let dmat: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        if i == j {
                            0.0
                        } else {
                            let d = ((pts[i].0 - pts[j].0).powi(2) + (pts[i].1 - pts[j].1).powi(2))
                                .sqrt();
                            if d < nearest[i] {
                                nearest[i] = d;
                            }
                            d
                        }
                    })
                    .collect()
            })
            .collect();

        // Default band sweep from the nearest-neighbour distances.
        let max_nn = nearest
            .iter()
            .copied()
            .filter(|d| d.is_finite())
            .fold(0.0, f64::max);
        let begin = prm.begin_distance.unwrap_or(max_nn.max(f64::MIN_POSITIVE));
        let increment = prm.increment.unwrap_or(begin.max(f64::MIN_POSITIVE));
        let num_bands = prm.num_bands;

        // Value statistics that are constant across bands.
        let mean = vals.iter().sum::<f64>() / n as f64;
        let dev: Vec<f64> = vals.iter().map(|v| v - mean).collect();
        let m2 = dev.iter().map(|d| d * d).sum::<f64>(); // Σ(x-x̄)²
        let m4 = dev.iter().map(|d| d.powi(4)).sum::<f64>();
        if m2 <= 0.0 {
            return Err(ToolError::Execution(
                "field is constant; Moran's I is undefined".to_string(),
            ));
        }
        let nf = n as f64;
        let kurt = nf * m4 / (m2 * m2); // K in Esri's variance formula
        let e_i = -1.0 / (nf - 1.0);

        ctx.progress
            .info(&format!("evaluating {num_bands} distance band(s)"));

        let mut csv = String::from("distance,morans_i,expected_i,variance,z_score,p_value\n");
        let mut rows: Vec<BandResult> = Vec::new();
        for b in 0..num_bands {
            let d = begin + increment * b as f64;
            let r = morans_i_band(&dmat, &dev, m2, nf, kurt, e_i, d);
            csv.push_str(&format!(
                "{:.4},{:.6},{:.6},{:.6e},{:.6},{:.6}\n",
                d, r.i, e_i, r.variance, r.z, r.p
            ));
            rows.push(BandResult { distance: d, ..r });
        }

        // First peak (first local maximum of z) and maximum peak.
        let max_peak = rows
            .iter()
            .filter(|r| r.z.is_finite())
            .max_by(|a, b| a.z.total_cmp(&b.z))
            .map(|r| r.distance);
        let first_peak = first_local_max(&rows);

        if let Some(path) = output {
            write_text_output(&csv, path)?;
        }

        ctx.progress.info(&format!(
            "max-peak distance {:?}, first-peak distance {:?}",
            max_peak, first_peak
        ));

        let table: Vec<Value> = rows
            .iter()
            .map(|r| {
                json!({
                    "distance": r.distance,
                    "morans_i": r.i,
                    "expected_i": e_i,
                    "variance": r.variance,
                    "z_score": r.z,
                    "p_value": r.p,
                })
            })
            .collect();

        let mut outputs = BTreeMap::new();
        if let Some(path) = output {
            outputs.insert("output".to_string(), json!(path));
        }
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("num_bands".to_string(), json!(num_bands));
        outputs.insert("begin_distance".to_string(), json!(begin));
        outputs.insert("increment".to_string(), json!(increment));
        outputs.insert("max_peak_distance".to_string(), json!(max_peak));
        outputs.insert("first_peak_distance".to_string(), json!(first_peak));
        outputs.insert("table".to_string(), json!(table));
        Ok(ToolRunResult { outputs })
    }
}

// ── Moran's I with randomization variance ────────────────────────────────────

#[derive(Clone)]
struct BandResult {
    distance: f64,
    i: f64,
    variance: f64,
    z: f64,
    p: f64,
}

/// Moran's I at a fixed distance band, with the randomization-null z-score.
fn morans_i_band(
    dmat: &[Vec<f64>],
    dev: &[f64],
    m2: f64,
    nf: f64,
    kurt: f64,
    e_i: f64,
    d: f64,
) -> BandResult {
    let n = dev.len();
    // Weighted cross-product and the weight sums S0, S2 (S1 = 2 S0 for symmetric
    // binary weights).
    let mut cross = 0.0;
    let mut s0 = 0.0;
    let mut degree = vec![0.0f64; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            if dmat[i][j] <= d {
                cross += dev[i] * dev[j];
                s0 += 1.0;
                degree[i] += 1.0;
            }
        }
    }
    if s0 <= 0.0 {
        return BandResult {
            distance: d,
            i: f64::NAN,
            variance: f64::NAN,
            z: f64::NAN,
            p: f64::NAN,
        };
    }
    let i_val = (nf / s0) * (cross / m2);

    // S1 = 2 S0 (symmetric binary), S2 = Σ (2 k_i)² = 4 Σ k_i².
    let s1 = 2.0 * s0;
    let s2 = 4.0 * degree.iter().map(|k| k * k).sum::<f64>();

    // Esri/randomization variance.
    let n2 = nf * nf;
    let a = nf * ((n2 - 3.0 * nf + 3.0) * s1 - nf * s2 + 3.0 * s0 * s0);
    let b = kurt * ((n2 - nf) * s1 - 2.0 * nf * s2 + 6.0 * s0 * s0);
    let denom = (nf - 1.0) * (nf - 2.0) * (nf - 3.0) * s0 * s0;
    let var = if denom != 0.0 {
        (a - b) / denom - e_i * e_i
    } else {
        f64::NAN
    };
    let (z, p) = if var > 0.0 {
        let z = (i_val - e_i) / var.sqrt();
        (z, 2.0 * (1.0 - normal_cdf(z.abs())))
    } else {
        (f64::NAN, f64::NAN)
    };
    BandResult {
        distance: d,
        i: i_val,
        variance: var,
        z,
        p,
    }
}

/// Distance of the first local maximum of the z-curve (a rise then a fall).
fn first_local_max(rows: &[BandResult]) -> Option<f64> {
    for w in rows.windows(3) {
        if w[1].z.is_finite() && w[1].z > w[0].z && w[1].z >= w[2].z {
            return Some(w[1].distance);
        }
    }
    // Monotone increasing: the last band; monotone decreasing: the first.
    rows.iter().find(|r| r.z.is_finite()).map(|r| r.distance)
}

fn normal_cdf(x: f64) -> f64 {
    0.5 * erfc(-x / std::f64::consts::SQRT_2)
}

fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let ans = t
        * (-z * z - 1.26551223
            + t * (1.00002368
                + t * (0.37409196
                    + t * (0.09678418
                        + t * (-0.18628806
                            + t * (0.27886807
                                + t * (-1.13520398
                                    + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
            .exp();
    if x >= 0.0 {
        ans
    } else {
        2.0 - ans
    }
}

// ── Geometry / parameters ────────────────────────────────────────────────────

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

struct Params {
    begin_distance: Option<f64>,
    increment: Option<f64>,
    num_bands: usize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let begin_distance = parse_optional_f64(args, "begin_distance")?;
    if let Some(v) = begin_distance {
        if !(v > 0.0 && v.is_finite()) {
            return Err(ToolError::Validation(
                "'begin_distance' must be a positive number".to_string(),
            ));
        }
    }
    let increment = parse_optional_f64(args, "increment")?;
    if let Some(v) = increment {
        if !(v > 0.0 && v.is_finite()) {
            return Err(ToolError::Validation(
                "'increment' must be a positive number".to_string(),
            ));
        }
    }
    let num_bands = match parse_optional_f64(args, "num_bands")? {
        None => 10,
        Some(v) if v.fract() == 0.0 && (2.0..=200.0).contains(&v) => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "'num_bands' must be an integer between 2 and 200".to_string(),
            ))
        }
    };
    Ok(Params {
        begin_distance,
        increment,
        num_bands,
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
    use wbvector::{memory_store, FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a grid of points with the given values.
    fn grid_layer(side: usize, value: impl Fn(usize, usize) -> f64) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for r in 0..side {
            for c in 0..side {
                l.add_feature(
                    Some(Geometry::point(c as f64, r as f64)),
                    &[("v", value(r, c).into())],
                )
                .unwrap();
            }
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        IncrementalSpatialAutocorrelationTool
            .run(&args, &ctx())
            .unwrap()
    }

    /// A smooth spatial gradient is strongly positively autocorrelated: Moran's I
    /// is high and z is very positive at short bands.
    #[test]
    fn gradient_is_positively_autocorrelated() {
        // value = row + col: neighbours have similar values.
        let input = grid_layer(10, |r, c| (r + c) as f64);
        let out = run(json!({
            "input": input, "field": "v", "begin_distance": 1.5, "increment": 1.0, "num_bands": 5,
        }));
        let table = out.outputs["table"].as_array().unwrap();
        let first = &table[0];
        assert!(
            first["morans_i"].as_f64().unwrap() > 0.3,
            "gradient Moran's I should be strongly positive, got {}",
            first["morans_i"]
        );
        assert!(
            first["z_score"].as_f64().unwrap() > 2.0,
            "z should be significant, got {}",
            first["z_score"]
        );
    }

    /// A checkerboard is negatively autocorrelated at the 1-cell band (Moran's I
    /// < 0).
    #[test]
    fn checkerboard_is_negatively_autocorrelated() {
        let input = grid_layer(10, |r, c| if (r + c) % 2 == 0 { 1.0 } else { 0.0 });
        let out = run(json!({
            "input": input, "field": "v", "begin_distance": 1.1, "increment": 1.0, "num_bands": 3,
        }));
        let table = out.outputs["table"].as_array().unwrap();
        // At the rook-neighbour band, alternating values -> negative I.
        assert!(
            table[0]["morans_i"].as_f64().unwrap() < 0.0,
            "checkerboard should be negatively autocorrelated, got {}",
            table[0]["morans_i"]
        );
    }

    /// The tool reports a peak distance and a well-formed table.
    #[test]
    fn reports_peaks_and_table() {
        let input = grid_layer(8, |r, c| ((r as f64) * 0.7 + (c as f64) * 0.3).sin());
        let out = run(json!({ "input": input, "field": "v", "num_bands": 6 }));
        assert_eq!(out.outputs["num_bands"], json!(6));
        assert!(out.outputs["max_peak_distance"].is_number());
        assert_eq!(out.outputs["table"].as_array().unwrap().len(), 6);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            IncrementalSpatialAutocorrelationTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "field": "v", "num_bands": 1 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "field": "v" })).is_ok());
    }
}
