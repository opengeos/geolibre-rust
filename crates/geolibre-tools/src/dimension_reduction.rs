//! GeoLibre tool: principal-component analysis over feature attributes.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Dimension Reduction* (Spatial
//! Statistics). The bundled `principal_component_analysis` operates on raster
//! imagery bands only — there is no PCA over *vector feature attributes*, yet
//! collinear attribute inputs quietly distort `calculate_composite_index`,
//! `similarity_search`, and `build_balanced_zones`. This tool reduces a set of
//! numeric fields to a smaller set of orthogonal principal components.
//!
//! Selected fields are optionally standardized (z-score, sample SD), the
//! covariance (raw) or correlation (standardized) matrix is formed, and its
//! symmetric eigen-problem is solved with a hand-rolled cyclic **Jacobi
//! rotation** (no linear-algebra crate). Features are projected onto the leading
//! components, written back as `PC1`, `PC2`, … score fields. A companion report
//! table lists each component's eigenvalue, variance explained, cumulative
//! variance, and per-variable loadings (the eigenvector coefficients).
//!
//! The number of retained components comes from `num_components`, or from the
//! smallest count reaching `min_variance` cumulative variance, else all of them.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue};

use crate::common::write_text_output;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct DimensionReductionTool;

impl Tool for DimensionReductionTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "dimension_reduction",
            display_name: "Dimension Reduction",
            summary: "Principal-component analysis over numeric feature attributes (like ArcGIS Dimension Reduction): standardize selected fields, eigen-decompose the correlation/covariance matrix with a hand-rolled Jacobi rotation, and write PC1..PCk component scores back to the features plus an eigenvalue/variance/loadings report — the vector-attribute counterpart of the bundled band-only principal_component_analysis.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (or table) with numeric attribute fields.",
                    required: true,
                },
                ToolParamSpec {
                    name: "fields",
                    description: "Comma-separated list of at least two numeric fields to reduce.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector layer with PC score fields appended. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "table",
                    description: "Optional CSV path for the components report (eigenvalue, variance explained, cumulative, per-variable loadings). Always returned in the result.",
                    required: false,
                },
                ToolParamSpec {
                    name: "num_components",
                    description: "Number of leading components to keep. Default: enough to reach 'min_variance', else all fields.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_variance",
                    description: "Keep the fewest components whose cumulative variance explained reaches this fraction (0-1). Ignored if 'num_components' is set.",
                    required: false,
                },
                ToolParamSpec {
                    name: "standardize",
                    description: "Standardize each field (z-score) before PCA so scales are comparable — analyse the correlation matrix. Default true; false analyses the covariance matrix.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "fields")?;
        parse_params(args)?;
        Ok(())
    }

    #[allow(clippy::needless_range_loop)] // dense matrix math reads clearest with explicit indices
    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let table_path = parse_optional_str(args, "table")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let n = layer.features.len();
        let p = prm.fields.len();
        if n < 2 {
            return Err(ToolError::Validation(
                "at least 2 features are required for PCA".to_string(),
            ));
        }

        // Resolve field indices.
        let mut idx = Vec::with_capacity(p);
        for f in &prm.fields {
            let i = layer
                .schema
                .field_index(f)
                .ok_or_else(|| ToolError::Validation(format!("field '{f}' not found")))?;
            idx.push(i);
        }

        // Read raw values per column (None = missing).
        let mut raw: Vec<Vec<Option<f64>>> = vec![vec![None; n]; p];
        for (fi, feat) in layer.features.iter().enumerate() {
            for (ci, &i) in idx.iter().enumerate() {
                raw[ci][fi] = feat
                    .attributes
                    .get(i)
                    .and_then(|v| v.as_f64())
                    .filter(|x| x.is_finite());
            }
        }

        // Column mean and sample SD over present values; impute missing with mean.
        let mut means = vec![0.0f64; p];
        let mut sds = vec![0.0f64; p];
        let mut imputed = 0usize;
        for ci in 0..p {
            let present: Vec<f64> = raw[ci].iter().filter_map(|v| *v).collect();
            let m = present.len();
            let mean = if m == 0 {
                0.0
            } else {
                present.iter().sum::<f64>() / m as f64
            };
            let var = if m > 1 {
                present.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (m as f64 - 1.0)
            } else {
                0.0
            };
            means[ci] = mean;
            sds[ci] = var.sqrt();
            imputed += n - m;
        }

        // Build the centered (and optionally standardized) data matrix X (n x p).
        let mut x = vec![vec![0.0f64; p]; n];
        for ci in 0..p {
            let scale = if prm.standardize && sds[ci] > 0.0 {
                1.0 / sds[ci]
            } else {
                1.0
            };
            for fi in 0..n {
                let raw_val = raw[ci][fi].unwrap_or(means[ci]);
                x[fi][ci] = (raw_val - means[ci]) * scale;
            }
        }

        // Symmetric p x p covariance/correlation matrix S = X^T X / (n-1).
        let mut s = vec![vec![0.0f64; p]; p];
        for a in 0..p {
            for b in a..p {
                let mut acc = 0.0;
                for fi in 0..n {
                    acc += x[fi][a] * x[fi][b];
                }
                let v = acc / (n as f64 - 1.0);
                s[a][b] = v;
                s[b][a] = v;
            }
        }

        ctx.progress
            .info(&format!("PCA over {p} field(s), {n} feature(s)"));

        // Eigen-decompose, then sort components by descending eigenvalue.
        let (eigvals, eigvecs) = jacobi_eigen(&s);
        let mut order: Vec<usize> = (0..p).collect();
        order.sort_by(|&a, &b| eigvals[b].total_cmp(&eigvals[a]));

        // Sorted eigenvalues and sign-fixed eigenvectors (largest |loading| positive).
        let mut vals = Vec::with_capacity(p);
        let mut vecs: Vec<Vec<f64>> = Vec::with_capacity(p); // vecs[comp][var]
        for &o in &order {
            let ev = eigvals[o].max(0.0);
            let mut vec_comp: Vec<f64> = (0..p).map(|var| eigvecs[var][o]).collect();
            let lead = vec_comp
                .iter()
                .cloned()
                .enumerate()
                .max_by(|a, b| a.1.abs().total_cmp(&b.1.abs()))
                .map(|(i, _)| i)
                .unwrap_or(0);
            if vec_comp[lead] < 0.0 {
                for c in vec_comp.iter_mut() {
                    *c = -*c;
                }
            }
            vals.push(ev);
            vecs.push(vec_comp);
        }

        let total: f64 = vals.iter().sum();
        let mut ratio = vec![0.0f64; p];
        let mut cumulative = vec![0.0f64; p];
        let mut run = 0.0;
        for i in 0..p {
            ratio[i] = if total > 0.0 { vals[i] / total } else { 0.0 };
            run += ratio[i];
            cumulative[i] = run;
        }

        // How many components to keep.
        let keep = match (prm.num_components, prm.min_variance) {
            (Some(k), _) => k.clamp(1, p),
            (None, Some(mv)) => {
                let mut k = p;
                for i in 0..p {
                    if cumulative[i] >= mv - 1e-12 {
                        k = i + 1;
                        break;
                    }
                }
                k
            }
            (None, None) => p,
        };

        // Project features onto the kept components: score = X . eigenvector.
        let mut scores = vec![vec![0.0f64; keep]; n];
        for fi in 0..n {
            for c in 0..keep {
                let mut acc = 0.0;
                for var in 0..p {
                    acc += x[fi][var] * vecs[c][var];
                }
                scores[fi][c] = acc;
            }
        }

        // Append PC score fields.
        for c in 0..keep {
            layer.add_field(FieldDef::new(format!("PC{}", c + 1), FieldType::Float));
        }
        for fi in 0..n {
            let feat = &mut layer.features[fi];
            for c in 0..keep {
                feat.attributes.push(FieldValue::Float(scores[fi][c]));
            }
        }

        // Build the components report (all p components; `kept` flags retention).
        let mut header =
            String::from("component,eigenvalue,variance_explained,cumulative_variance,kept");
        for f in &prm.fields {
            header.push_str(&format!(",loading_{f}"));
        }
        let mut csv = header.clone();
        csv.push('\n');
        let mut report: Vec<Value> = Vec::with_capacity(p);
        for i in 0..p {
            csv.push_str(&format!(
                "PC{},{:.10},{:.10},{:.10},{}",
                i + 1,
                vals[i],
                ratio[i],
                cumulative[i],
                if i < keep { 1 } else { 0 },
            ));
            let mut loadings = serde_json::Map::new();
            for (vi, f) in prm.fields.iter().enumerate() {
                csv.push_str(&format!(",{:.10}", vecs[i][vi]));
                loadings.insert(f.clone(), json!(vecs[i][vi]));
            }
            csv.push('\n');
            report.push(json!({
                "component": format!("PC{}", i + 1),
                "eigenvalue": vals[i],
                "variance_explained": ratio[i],
                "cumulative_variance": cumulative[i],
                "kept": i < keep,
                "loadings": Value::Object(loadings),
            }));
        }

        if let Some(path) = table_path {
            write_text_output(&csv, path)?;
        }

        let out_path = write_or_store_layer(layer, output)?;

        ctx.progress.info(&format!(
            "kept {keep} of {p} component(s); {:.1}% variance",
            cumulative[keep - 1] * 100.0
        ));

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("num_fields".to_string(), json!(p));
        outputs.insert("num_components".to_string(), json!(keep));
        outputs.insert("eigenvalues".to_string(), json!(vals));
        outputs.insert("variance_explained".to_string(), json!(ratio));
        outputs.insert("cumulative_variance".to_string(), json!(cumulative));
        outputs.insert("cumulative_kept".to_string(), json!(cumulative[keep - 1]));
        outputs.insert("imputed_values".to_string(), json!(imputed));
        outputs.insert("report".to_string(), json!(report));
        if let Some(path) = table_path {
            outputs.insert("table".to_string(), json!(path));
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── Eigen-decomposition (cyclic Jacobi rotation for a symmetric matrix) ────────

/// Solves the symmetric eigen-problem `S v = λ v` via cyclic Jacobi rotations.
/// Returns `(eigenvalues, eigenvectors)` where `eigenvectors[row][col]` holds the
/// `col`-th eigenvector (a column). No sorting is applied here.
#[allow(clippy::needless_range_loop)] // Jacobi rotations index rows/cols directly
fn jacobi_eigen(s: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = s.len();
    let mut a: Vec<Vec<f64>> = s.to_vec();
    let mut v = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        v[i][i] = 1.0;
    }
    if n <= 1 {
        let vals = (0..n).map(|i| a[i][i]).collect();
        return (vals, v);
    }

    for _sweep in 0..100 {
        // Sum of squared off-diagonal magnitudes; stop once negligible.
        let mut off = 0.0;
        for p in 0..n {
            for q in (p + 1)..n {
                off += a[p][q] * a[p][q];
            }
        }
        if off < 1e-28 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p][q];
                if apq.abs() < 1e-300 {
                    continue;
                }
                let app = a[p][p];
                let aqq = a[q][q];
                let theta = (aqq - app) / (2.0 * apq);
                let t = if theta == 0.0 {
                    1.0
                } else {
                    let sign = if theta > 0.0 { 1.0 } else { -1.0 };
                    sign / (theta.abs() + (theta * theta + 1.0).sqrt())
                };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let sn = t * c;
                // Rotate columns p,q of A.
                for k in 0..n {
                    let akp = a[k][p];
                    let akq = a[k][q];
                    a[k][p] = c * akp - sn * akq;
                    a[k][q] = sn * akp + c * akq;
                }
                // Rotate rows p,q of A.
                for k in 0..n {
                    let apk = a[p][k];
                    let aqk = a[q][k];
                    a[p][k] = c * apk - sn * aqk;
                    a[q][k] = sn * apk + c * aqk;
                }
                // Accumulate rotation into V.
                for k in 0..n {
                    let vkp = v[k][p];
                    let vkq = v[k][q];
                    v[k][p] = c * vkp - sn * vkq;
                    v[k][q] = sn * vkp + c * vkq;
                }
            }
        }
    }

    let vals = (0..n).map(|i| a[i][i]).collect();
    (vals, v)
}

// ── Parameters ─────────────────────────────────────────────────────────────────

struct Params {
    fields: Vec<String>,
    num_components: Option<usize>,
    min_variance: Option<f64>,
    standardize: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let fields: Vec<String> = require_str(args, "fields")?
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect();
    if fields.len() < 2 {
        return Err(ToolError::Validation(
            "'fields' must name at least two numeric fields".to_string(),
        ));
    }

    let num_components = match opt_f64(args, "num_components")? {
        Some(v) if v >= 1.0 && v.is_finite() => Some(v.round() as usize),
        Some(_) => {
            return Err(ToolError::Validation(
                "'num_components' must be a positive integer".to_string(),
            ))
        }
        None => None,
    };
    let min_variance = match opt_f64(args, "min_variance")? {
        Some(v) if (0.0..=1.0).contains(&v) => Some(v),
        Some(_) => {
            return Err(ToolError::Validation(
                "'min_variance' must be between 0 and 1".to_string(),
            ))
        }
        None => None,
    };
    let standardize = opt_bool(args, "standardize")?.unwrap_or(true);

    Ok(Params {
        fields,
        num_components,
        min_variance,
        standardize,
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

fn opt_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" | "y" => Ok(Some(true)),
            "false" | "0" | "no" | "n" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Geometry, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a point layer with three numeric fields a, b, c.
    fn layer(rows: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("z")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("a", FieldType::Float));
        l.add_field(FieldDef::new("b", FieldType::Float));
        l.add_field(FieldDef::new("c", FieldType::Float));
        for (i, (a, b, c)) in rows.iter().enumerate() {
            l.add_feature(
                Some(Geometry::point(i as f64, 0.0)),
                &[("a", (*a).into()), ("b", (*b).into()), ("c", (*c).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = DimensionReductionTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn col(l: &Layer, name: &str) -> Vec<f64> {
        let i = l.schema.field_index(name).unwrap();
        l.features
            .iter()
            .map(|f| f.attributes[i].as_f64().unwrap())
            .collect()
    }

    fn corr(a: &[f64], b: &[f64]) -> f64 {
        let n = a.len() as f64;
        let ma = a.iter().sum::<f64>() / n;
        let mb = b.iter().sum::<f64>() / n;
        let mut num = 0.0;
        let mut da = 0.0;
        let mut db = 0.0;
        for i in 0..a.len() {
            num += (a[i] - ma) * (b[i] - mb);
            da += (a[i] - ma).powi(2);
            db += (b[i] - mb).powi(2);
        }
        if da == 0.0 || db == 0.0 {
            0.0
        } else {
            num / (da.sqrt() * db.sqrt())
        }
    }

    /// Eigenvalues of a correlation matrix sum to the field count, variance
    /// ratios sum to 1, and successive component scores are orthogonal
    /// (uncorrelated) — the defining PCA properties.
    #[test]
    fn pca_properties_hold() {
        // a and b are strongly collinear; c is nearly independent.
        let rows = [
            (1.0, 2.1, 9.0),
            (2.0, 4.0, 3.0),
            (3.0, 5.9, 7.0),
            (4.0, 8.2, 1.0),
            (5.0, 9.8, 6.0),
            (6.0, 12.1, 2.0),
        ];
        let (out, layer) = run(json!({ "input": layer(&rows), "fields": "a,b,c" }));
        let eig = out.outputs["eigenvalues"].as_array().unwrap();
        let sum: f64 = eig.iter().map(|v| v.as_f64().unwrap()).sum();
        assert!((sum - 3.0).abs() < 1e-6, "corr-matrix eigenvalues sum to p");
        let ratios = out.outputs["variance_explained"].as_array().unwrap();
        let rsum: f64 = ratios.iter().map(|v| v.as_f64().unwrap()).sum();
        assert!((rsum - 1.0).abs() < 1e-9, "variance ratios sum to 1");
        // Descending order.
        let e0 = eig[0].as_f64().unwrap();
        let e1 = eig[1].as_f64().unwrap();
        assert!(e0 >= e1);
        // PC1 and PC2 are orthogonal.
        let pc1 = col(&layer, "PC1");
        let pc2 = col(&layer, "PC2");
        assert!(corr(&pc1, &pc2).abs() < 1e-9, "components are uncorrelated");
        // First component captures the a/b collinearity: dominant variance.
        assert!(e0 > 1.5, "leading component captures the collinear pair");
    }

    /// A perfectly collinear pair collapses to essentially one component.
    #[test]
    fn collinear_pair_one_component() {
        // b = 2*a exactly; c = -a. All three are perfectly correlated.
        let rows = [
            (1.0, 2.0, -1.0),
            (2.0, 4.0, -2.0),
            (3.0, 6.0, -3.0),
            (4.0, 8.0, -4.0),
            (5.0, 10.0, -5.0),
        ];
        let (out, _l) = run(json!({ "input": layer(&rows), "fields": "a,b,c" }));
        let ratios = out.outputs["variance_explained"].as_array().unwrap();
        assert!(
            ratios[0].as_f64().unwrap() > 0.999,
            "one component explains ~all variance for collinear inputs"
        );
    }

    /// min_variance selects the fewest components reaching the threshold.
    #[test]
    fn min_variance_selects_components() {
        let rows = [
            (1.0, 2.0, -1.0),
            (2.0, 4.0, -2.0),
            (3.0, 6.0, -3.0),
            (4.0, 8.0, -4.0),
            (5.0, 10.0, -5.0),
        ];
        let (out, layer) = run(json!({
            "input": layer(&rows), "fields": "a,b,c", "min_variance": 0.9,
        }));
        assert_eq!(
            out.outputs["num_components"],
            json!(1),
            "collinear inputs need only 1 component to reach 90%"
        );
        assert!(layer.schema.field_index("PC1").is_some());
        assert!(layer.schema.field_index("PC2").is_none());
    }

    /// num_components caps the kept components and PC fields.
    #[test]
    fn num_components_caps_output() {
        let rows = [
            (1.0, 2.1, 9.0),
            (2.0, 4.0, 3.0),
            (3.0, 5.9, 7.0),
            (4.0, 8.2, 1.0),
        ];
        let (out, layer) = run(json!({
            "input": layer(&rows), "fields": "a,b,c", "num_components": 2,
        }));
        assert_eq!(out.outputs["num_components"], json!(2));
        assert!(layer.schema.field_index("PC2").is_some());
        assert!(layer.schema.field_index("PC3").is_none());
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            DimensionReductionTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err()); // no fields
        assert!(bad(json!({ "input": "a.geojson", "fields": "a" })).is_err()); // need >=2
        assert!(
            bad(json!({ "input": "a.geojson", "fields": "a,b", "min_variance": 2.0 })).is_err()
        );
        assert!(
            bad(json!({ "input": "a.geojson", "fields": "a,b", "num_components": 0 })).is_err()
        );
        assert!(bad(json!({ "input": "a.geojson", "fields": "a,b" })).is_ok());
    }
}
