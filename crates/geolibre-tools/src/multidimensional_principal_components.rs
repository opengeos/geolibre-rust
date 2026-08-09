//! GeoLibre tool: principal components of a raster cube.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Multidimensional Principal
//! Components* (Image Analyst).
//!
//! ## Why the catalog needs it
//!
//! A cube of 200 monthly NDVI slices is not 200 independent pictures — nearly
//! all of it is one seasonal cycle plus a trend. PCA along the slice dimension
//! finds those few real modes: the first components map *where* each mode acts,
//! and the eigenvalues say how much of the record each explains. It is the
//! standard way to compress a long series, to separate a seasonal signal from
//! an interannual one, and to spot the handful of pixels that do not follow the
//! regional pattern.
//!
//! The bundled `principal_component_analysis` cannot do this: it decomposes
//! *band space within a single scene* (the classic multispectral PCA), so its
//! components mix spectral bands, not time. `dimension_reduction` operates on
//! vector attribute tables. `multidimensional_anomaly` removes a baseline but
//! finds no modes, and `time_series_clustering` groups whole series rather than
//! decomposing them.
//!
//! ## The two reductions, and why one eigendecomposition serves both
//!
//! * **`dimension_reduction`** treats each slice as a variable and each pixel
//!   as an observation. Components are new *rasters* (spatial patterns); the
//!   loadings say how each input slice contributes to each.
//! * **`spatial_reduction`** treats each pixel as a variable and each slice as
//!   an observation. Components are *series* over the dimension; the loadings
//!   are rasters.
//!
//! These are transposes of one another, and both follow from the eigenvectors
//! of the same `n_slices x n_slices` covariance matrix. That matters
//! practically: the pixel-space covariance of even a small cube would be
//! millions of rows square, while the slice-space one is a few hundred at
//! worst and is solved exactly by the shared cyclic-Jacobi routine.
//!
//! Only cells valid in **every** slice take part, because a covariance built
//! from different pixel subsets per slice pair is not a consistent matrix; the
//! count of participating cells is reported.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;
use wbvector::{FieldDef, FieldType, FieldValue, Layer};

use crate::args_common::{bool_or, choice_or, opt_usize};
use crate::common::{parse_optional_output, write_or_store_output};
use crate::cube::{load_cube, Cube};
use crate::dimension_reduction::jacobi_eigen;
use crate::raster_stack::raster_like_multiband;
use crate::vector_common::write_or_store_layer;

pub struct MultidimensionalPrincipalComponentsTool;

impl Tool for MultidimensionalPrincipalComponentsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "multidimensional_principal_components",
            display_name: "Multidimensional Principal Components",
            summary: "Decomposes a raster cube into the few real modes underlying a long time series, emitting component rasters, loadings and eigenvalues, in either dimension-reduction or spatial-reduction form (ArcGIS Multidimensional Principal Components). The bundled principal_component_analysis decomposes band space within one scene, so its components mix spectral bands rather than time; dimension_reduction works on vector attribute tables. Both reductions come from one slice-space covariance, since the pixel-space equivalent would be millions of rows square.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "One multiband raster (each band is a slice) or a comma-separated list of co-registered rasters.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output multiband raster: component scores for 'dimension_reduction', component loadings for 'spatial_reduction'. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_loadings",
                    description: "Output table of the per-slice loadings (or per-slice component series under 'spatial_reduction'). Always produced; stored in memory when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_eigenvalues",
                    description: "Output table of eigenvalues and the variance each explains. Always produced; stored in memory when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "mode",
                    description: "'dimension_reduction' (default; components are spatial patterns) or 'spatial_reduction' (components are series over the dimension).",
                    required: false,
                },
                ToolParamSpec {
                    name: "dimension",
                    description: "Name of the dimension, used only in the report (default 'slice').",
                    required: false,
                },
                ToolParamSpec {
                    name: "dimension_values",
                    description: "Comma-separated coordinate of each slice, strictly increasing. Defaults to the 1-based slice index.",
                    required: false,
                },
                ToolParamSpec {
                    name: "number_of_pc",
                    description: "How many components to keep (default: all slices).",
                    required: false,
                },
                ToolParamSpec {
                    name: "correlation",
                    description: "Decompose the correlation matrix instead of the covariance (default false). Use when slices are on different scales.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        crate::raster_stack::parse_input_paths(args, "input")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let prm = parse_params(args)?;
        let output = parse_optional_output(args, "output")?;
        let out_loadings = parse_optional_output(args, "output_loadings")?;
        let out_eigen = parse_optional_output(args, "output_eigenvalues")?;

        let cube = load_cube(args, "input", "dimension_values", "dimension", 2)?;
        let (rows, cols, n) = (cube.rows, cube.cols, cube.len());

        // Complete cases only: a covariance whose entries came from different
        // pixel subsets is not a consistent matrix and can fail to be positive
        // semi-definite.
        let complete: Vec<usize> = (0..rows * cols)
            .filter(|&i| (0..n).all(|s| cube.get(s, i / cols, i % cols).is_some()))
            .collect();
        if complete.len() < 2 {
            return Err(ToolError::Execution(format!(
                "only {} cell(s) are valid in every slice; PCA needs at least 2",
                complete.len()
            )));
        }

        ctx.progress.info(&format!(
            "{n} slice(s) over '{}', {} complete cell(s), mode {}",
            cube.dimension,
            complete.len(),
            prm.mode.label()
        ));

        // Slice means and the slice-space covariance (or correlation).
        let m = complete.len() as f64;
        let mut means = vec![0.0f64; n];
        for s in 0..n {
            let sum: f64 = complete
                .iter()
                .map(|&i| cube.get(s, i / cols, i % cols).unwrap())
                .sum();
            means[s] = sum / m;
        }
        let mut cov = vec![vec![0.0f64; n]; n];
        for &i in &complete {
            let (r, c) = (i / cols, i % cols);
            let dev: Vec<f64> = (0..n)
                .map(|s| cube.get(s, r, c).unwrap() - means[s])
                .collect();
            for a in 0..n {
                for b in a..n {
                    cov[a][b] += dev[a] * dev[b];
                }
            }
        }
        for a in 0..n {
            for b in a..n {
                cov[a][b] /= m;
                cov[b][a] = cov[a][b];
            }
        }

        // Standard deviations, kept for the correlation form and for reporting.
        let sds: Vec<f64> = (0..n).map(|s| cov[s][s].max(0.0).sqrt()).collect();
        if prm.correlation {
            for a in 0..n {
                for b in 0..n {
                    // A slice with no variance cannot be standardised; leave it
                    // at zero so it contributes to no component rather than
                    // producing a NaN that would poison every eigenvector.
                    let d = sds[a] * sds[b];
                    cov[a][b] = if d > 0.0 { cov[a][b] / d } else { 0.0 };
                }
            }
        }

        let (mut values, vectors) = jacobi_eigen(&cov);
        // Order components by explained variance, descending.
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| values[b].total_cmp(&values[a]));
        let keep = prm.number_of_pc.unwrap_or(n).min(n);

        // Sign convention: the largest-magnitude loading of each component is
        // positive. An eigenvector's sign is arbitrary, so without this the
        // same input can produce components that flip between runs of a
        // different length and comparisons become meaningless.
        let mut loadings: Vec<Vec<f64>> = Vec::with_capacity(keep);
        for k in 0..keep {
            let col = order[k];
            let mut v: Vec<f64> = (0..n).map(|r| vectors[r][col]).collect();
            let dominant = v
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.abs().total_cmp(&b.1.abs()))
                .map(|(i, _)| i)
                .unwrap_or(0);
            if v[dominant] < 0.0 {
                for x in v.iter_mut() {
                    *x = -*x;
                }
            }
            loadings.push(v);
        }
        let eigenvalues: Vec<f64> = (0..keep).map(|k| values[order[k]].max(0.0)).collect();
        let total: f64 = values.iter().map(|v| v.max(0.0)).sum();
        values.truncate(n); // keep clippy from flagging the unused tail

        // Project every complete cell onto the retained components.
        let nodata = -9999.0_f64;
        let mut bands: Vec<Vec<f64>> = vec![vec![nodata; rows * cols]; keep];
        let mut series: Vec<Vec<f64>> = vec![vec![0.0; n]; keep];
        for &i in &complete {
            let (r, c) = (i / cols, i % cols);
            let dev: Vec<f64> = (0..n)
                .map(|s| {
                    let d = cube.get(s, r, c).unwrap() - means[s];
                    // The correlation form decomposes standardised variables,
                    // so the projection must standardise too.
                    if prm.correlation && sds[s] > 0.0 {
                        d / sds[s]
                    } else {
                        d
                    }
                })
                .collect();
            for k in 0..keep {
                let score: f64 = dev.iter().zip(&loadings[k]).map(|(d, l)| d * l).sum();
                bands[k][i] = score;
                // Spatial reduction wants the component's series over the
                // dimension; accumulate the score-weighted deviations.
                if prm.mode == Mode::Spatial {
                    for (s, d) in dev.iter().enumerate() {
                        series[k][s] += score * d;
                    }
                }
            }
        }

        // Under spatial reduction the roles swap: the rasters carry the
        // loadings and the table carries the component series.
        let (raster_bands, table_rows): (Vec<Vec<f64>>, Vec<Vec<f64>>) = match prm.mode {
            Mode::Dimension => (bands, loadings.clone()),
            Mode::Spatial => {
                // Normalise each component's series so it is comparable across
                // components rather than scaled by the cell count.
                let normalised: Vec<Vec<f64>> = series
                    .iter()
                    .map(|v| {
                        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
                        if norm > 0.0 {
                            v.iter().map(|x| x / norm).collect()
                        } else {
                            v.clone()
                        }
                    })
                    .collect();
                (bands, normalised)
            }
        };

        let out_raster =
            raster_like_multiband(cube.template(), &raster_bands, nodata, DataType::F32)?;
        let out_path = write_or_store_output(out_raster, output)?;

        let loadings_layer = loadings_table(&cube, &table_rows, prm.mode)?;
        let loadings_path = write_or_store_layer(loadings_layer, out_loadings)?;

        let eigen_layer = eigen_table(&eigenvalues, total)?;
        let eigen_path = write_or_store_layer(eigen_layer, out_eigen)?;

        let explained: Vec<f64> = eigenvalues
            .iter()
            .map(|v| if total > 0.0 { 100.0 * v / total } else { 0.0 })
            .collect();

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("output_loadings".to_string(), json!(loadings_path));
        outputs.insert("output_eigenvalues".to_string(), json!(eigen_path));
        outputs.insert("mode".to_string(), json!(prm.mode.label()));
        outputs.insert("dimension".to_string(), json!(cube.dimension));
        outputs.insert("slices".to_string(), json!(n));
        outputs.insert("components".to_string(), json!(keep));
        outputs.insert("complete_cells".to_string(), json!(complete.len()));
        outputs.insert("eigenvalues".to_string(), json!(eigenvalues));
        outputs.insert("variance_explained".to_string(), json!(explained));
        Ok(ToolRunResult { outputs })
    }
}

/// One row per slice, one column per retained component.
fn loadings_table(cube: &Cube, rows_data: &[Vec<f64>], mode: Mode) -> Result<Layer, ToolError> {
    let mut layer = Layer::new(match mode {
        Mode::Dimension => "pc_loadings",
        Mode::Spatial => "pc_series",
    });
    // The dimension's own name is used for its column, so a caller who asked
    // for "year" gets a `year` column rather than a generic one.
    let dim_col = cube.dimension.clone();
    let pc_cols: Vec<String> = (0..rows_data.len())
        .map(|k| format!("pc{}", k + 1))
        .collect();

    layer.add_field(FieldDef::new("slice", FieldType::Integer));
    layer.add_field(FieldDef::new(&dim_col, FieldType::Float));
    for name in &pc_cols {
        layer.add_field(FieldDef::new(name, FieldType::Float));
    }
    for s in 0..cube.len() {
        let mut attrs: Vec<(&str, FieldValue)> = vec![
            ("slice", FieldValue::Integer(s as i64 + 1)),
            (dim_col.as_str(), FieldValue::Float(cube.coord(s))),
        ];
        for (name, row) in pc_cols.iter().zip(rows_data) {
            attrs.push((
                name.as_str(),
                FieldValue::Float(row.get(s).copied().unwrap_or(0.0)),
            ));
        }
        layer
            .add_feature(None, &attrs)
            .map_err(|e| ToolError::Execution(format!("writing the loadings table: {e}")))?;
    }
    Ok(layer)
}

/// One row per retained component.
fn eigen_table(eigenvalues: &[f64], total: f64) -> Result<Layer, ToolError> {
    let mut layer = Layer::new("pc_eigenvalues");
    layer.add_field(FieldDef::new("component", FieldType::Integer));
    layer.add_field(FieldDef::new("eigenvalue", FieldType::Float));
    layer.add_field(FieldDef::new("variance_explained", FieldType::Float));
    layer.add_field(FieldDef::new("cumulative_explained", FieldType::Float));

    let mut cumulative = 0.0;
    for (k, &v) in eigenvalues.iter().enumerate() {
        let pct = if total > 0.0 { 100.0 * v / total } else { 0.0 };
        cumulative += pct;
        layer
            .add_feature(
                None,
                &[
                    ("component", FieldValue::Integer(k as i64 + 1)),
                    ("eigenvalue", FieldValue::Float(v)),
                    ("variance_explained", FieldValue::Float(pct)),
                    ("cumulative_explained", FieldValue::Float(cumulative)),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("writing the eigenvalue table: {e}")))?;
    }
    Ok(layer)
}

// ── Parameters ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    Dimension,
    Spatial,
}

impl Mode {
    fn label(self) -> &'static str {
        match self {
            Mode::Dimension => "dimension_reduction",
            Mode::Spatial => "spatial_reduction",
        }
    }
}

struct Params {
    mode: Mode,
    number_of_pc: Option<usize>,
    correlation: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let mode = match choice_or(
        args,
        "mode",
        &["dimension_reduction", "spatial_reduction"],
        "dimension_reduction",
    )? {
        "spatial_reduction" => Mode::Spatial,
        _ => Mode::Dimension,
    };
    let number_of_pc = match opt_usize(args, "number_of_pc")? {
        None => None,
        Some(k) if k >= 1 => Some(k),
        Some(_) => {
            return Err(ToolError::Validation(
                "'number_of_pc' must be at least 1".to_string(),
            ))
        }
    };
    let correlation = bool_or(args, "correlation", false)?;
    Ok(Params {
        mode,
        number_of_pc,
        correlation,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::load_input_raster;
    use crate::cube::test_support::cube_raster;
    use crate::vector_common::load_input_layer;
    use wbvector::Feature;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::Raster;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn run(args: Value) -> (Raster, BTreeMap<String, Value>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = MultidimensionalPrincipalComponentsTool
            .run(&args, &ctx())
            .unwrap();
        let r = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (r, out.outputs)
    }

    /// A cube whose slices are all multiples of one spatial pattern has exactly
    /// one real mode: PC1 must capture essentially all the variance, and the
    /// remaining components must be empty.
    #[test]
    fn one_mode_cube_concentrates_in_pc1() {
        let (rows, cols) = (6, 6);
        // Spatial pattern: a diagonal gradient.
        let pattern: Vec<f64> = (0..rows * cols)
            .map(|i| (i / cols) as f64 + (i % cols) as f64)
            .collect();
        // Four slices, each the pattern scaled by a different amount.
        let slices: Vec<Vec<f64>> = [1.0, 2.0, -1.0, 3.0]
            .iter()
            .map(|k| pattern.iter().map(|v| k * v).collect())
            .collect();
        let (_, outputs) = run(json!({ "input": cube_raster(cols, rows, &slices) }));

        let explained = outputs["variance_explained"].as_array().unwrap();
        let pc1 = explained[0].as_f64().unwrap();
        assert!(
            pc1 > 99.9,
            "a single-mode cube should put everything in PC1, got {pc1}%"
        );
        assert_eq!(outputs["components"].as_u64().unwrap(), 4);
        assert_eq!(outputs["complete_cells"].as_u64().unwrap(), 36);
    }

    /// Two independent spatial patterns produce two real modes, and the
    /// loadings identify which slices carry which.
    #[test]
    fn two_modes_are_separated() {
        let (rows, cols) = (8, 8);
        // Pattern A varies along rows, pattern B along columns — orthogonal.
        let a: Vec<f64> = (0..rows * cols)
            .map(|i| (i / cols) as f64 - 3.5)
            .collect();
        let b: Vec<f64> = (0..rows * cols)
            .map(|i| (i % cols) as f64 - 3.5)
            .collect();
        // Slices 0,1 are pattern A; slices 2,3 are pattern B.
        let slices: Vec<Vec<f64>> = vec![
            a.clone(),
            a.iter().map(|v| 2.0 * v).collect(),
            b.clone(),
            b.iter().map(|v| 2.0 * v).collect(),
        ];
        let (_, outputs) = run(json!({ "input": cube_raster(cols, rows, &slices) }));
        let explained = outputs["variance_explained"].as_array().unwrap();
        let two = explained[0].as_f64().unwrap() + explained[1].as_f64().unwrap();
        assert!(
            two > 99.9,
            "two orthogonal patterns should need exactly two components, got {two}%"
        );

        // The loadings table separates them: PC1 loads on one pair of slices,
        // PC2 on the other.
        let path = outputs["output_loadings"].as_str().unwrap();
        let layer = load_input_layer(path).unwrap();
        assert_eq!(layer.len(), 4, "one row per slice");
        let pc1 = layer.schema.field_index("pc1").unwrap();
        let pc2 = layer.schema.field_index("pc2").unwrap();
        let get = |f: &Feature, i: usize| match f.attributes[i] {
            FieldValue::Float(v) => v,
            _ => panic!("loading must be a float"),
        };
        let feats: Vec<&Feature> = layer.iter().collect();
        // Slices 0/1 load on one component and 2/3 on the other; which is which
        // depends on their variance, so test the block structure instead.
        let a_on_1 = get(feats[0], pc1).abs() + get(feats[1], pc1).abs();
        let b_on_1 = get(feats[2], pc1).abs() + get(feats[3], pc1).abs();
        let a_on_2 = get(feats[0], pc2).abs() + get(feats[1], pc2).abs();
        let b_on_2 = get(feats[2], pc2).abs() + get(feats[3], pc2).abs();
        assert!(
            (a_on_1 > 0.9 && b_on_1 < 0.1 && b_on_2 > 0.9 && a_on_2 < 0.1)
                || (b_on_1 > 0.9 && a_on_1 < 0.1 && a_on_2 > 0.9 && b_on_2 < 0.1),
            "loadings did not separate the two patterns: {a_on_1} {b_on_1} {a_on_2} {b_on_2}"
        );
    }

    /// Eigenvalues are ordered and the explained variance accumulates to 100%.
    #[test]
    fn eigenvalue_table_is_ordered_and_complete() {
        let (rows, cols) = (5, 5);
        let slices: Vec<Vec<f64>> = (0..4)
            .map(|k| {
                (0..rows * cols)
                    .map(|i| ((i * 7 + k * 13) % 11) as f64)
                    .collect()
            })
            .collect();
        let (_, outputs) = run(json!({ "input": cube_raster(cols, rows, &slices) }));
        let path = outputs["output_eigenvalues"].as_str().unwrap();
        let layer = load_input_layer(path).unwrap();
        assert_eq!(layer.len(), 4);

        let ev = layer.schema.field_index("eigenvalue").unwrap();
        let cum = layer.schema.field_index("cumulative_explained").unwrap();
        let vals: Vec<f64> = layer
            .iter()
            .map(|f| match f.attributes[ev] {
                FieldValue::Float(v) => v,
                _ => panic!(),
            })
            .collect();
        assert!(
            vals.windows(2).all(|w| w[0] >= w[1] - 1e-9),
            "eigenvalues must be descending, got {vals:?}"
        );
        let last = layer.iter().last().unwrap();
        let FieldValue::Float(total) = last.attributes[cum] else {
            panic!()
        };
        assert!(
            (total - 100.0).abs() < 1e-6,
            "cumulative variance should reach 100%, got {total}"
        );
    }

    /// `number_of_pc` truncates the output.
    #[test]
    fn number_of_pc_truncates() {
        let (rows, cols) = (5, 5);
        let slices: Vec<Vec<f64>> = (0..5)
            .map(|k| {
                (0..rows * cols)
                    .map(|i| ((i * 3 + k * 5) % 7) as f64)
                    .collect()
            })
            .collect();
        let (out, outputs) = run(json!({
            "input": cube_raster(cols, rows, &slices), "number_of_pc": 2
        }));
        assert_eq!(out.bands, 2);
        assert_eq!(outputs["components"].as_u64().unwrap(), 2);
    }

    /// The correlation form standardises slices, so a slice with a much larger
    /// scale no longer dominates PC1 purely by virtue of its units.
    #[test]
    fn correlation_mode_removes_scale_dominance() {
        let (rows, cols) = (6, 6);
        let p: Vec<f64> = (0..rows * cols).map(|i| (i % 5) as f64).collect();
        let q: Vec<f64> = (0..rows * cols).map(|i| (i / cols) as f64).collect();
        // Slice 0 is 1000x the scale of slice 1, but carries no more structure.
        let slices = vec![
            p.iter().map(|v| 1000.0 * v).collect::<Vec<f64>>(),
            q.clone(),
        ];
        let src = cube_raster(cols, rows, &slices);

        let cov = run(json!({ "input": src.clone() })).1;
        let cov_pc1 = cov["variance_explained"][0].as_f64().unwrap();
        let cor = run(json!({ "input": src, "correlation": true })).1;
        let cor_pc1 = cor["variance_explained"][0].as_f64().unwrap();
        assert!(
            cov_pc1 > 99.9,
            "on covariance the large-scale slice should dominate, got {cov_pc1}%"
        );
        assert!(
            cor_pc1 < cov_pc1 - 20.0,
            "correlation mode should break that dominance: {cov_pc1}% -> {cor_pc1}%"
        );
    }

    /// Spatial reduction emits the component series over the dimension.
    #[test]
    fn spatial_reduction_emits_series() {
        let (rows, cols) = (5, 5);
        let pattern: Vec<f64> = (0..rows * cols).map(|i| (i % 4) as f64).collect();
        let slices: Vec<Vec<f64>> = [1.0, 3.0, 2.0]
            .iter()
            .map(|k| pattern.iter().map(|v| k * v).collect())
            .collect();
        let (_, outputs) = run(json!({
            "input": cube_raster(cols, rows, &slices),
            "mode": "spatial_reduction", "dimension": "year",
            "dimension_values": "2001,2002,2003"
        }));
        assert_eq!(outputs["mode"].as_str().unwrap(), "spatial_reduction");
        let layer = load_input_layer(outputs["output_loadings"].as_str().unwrap()).unwrap();
        assert_eq!(layer.len(), 3, "one row per slice");
        // The dimension coordinate is carried through under its own name.
        let yi = layer.schema.field_index("year").unwrap();
        let FieldValue::Float(y0) = layer.iter().next().unwrap().attributes[yi] else {
            panic!("dimension coordinate must be a float")
        };
        assert!((y0 - 2001.0).abs() < 1e-9);
    }

    /// Cells missing in any slice are excluded, and the count says so.
    #[test]
    fn incomplete_cells_are_excluded() {
        let (rows, cols) = (3, 3);
        let mut s0: Vec<f64> = (0..9).map(|i| i as f64).collect();
        let s1: Vec<f64> = (0..9).map(|i| 2.0 * i as f64).collect();
        s0[4] = -9999.0; // one cell missing in slice 0
        let (out, outputs) = run(json!({ "input": cube_raster(cols, rows, &[s0, s1]) }));
        assert_eq!(outputs["complete_cells"].as_u64().unwrap(), 8);
        assert_eq!(
            out.get(0, 1, 1),
            -9999.0,
            "the incomplete cell must stay no-data in the scores"
        );
    }

    /// Too few complete cells is an error, not a degenerate answer.
    #[test]
    fn errors_when_nothing_is_complete() {
        let s0 = vec![-9999.0; 9];
        let s1: Vec<f64> = (0..9).map(|i| i as f64).collect();
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": cube_raster(3, 3, &[s0, s1]) })).unwrap();
        let err = MultidimensionalPrincipalComponentsTool
            .run(&args, &ctx())
            .unwrap_err();
        assert!(
            format!("{err:?}").contains("valid in every slice"),
            "expected a complete-case error, got {err:?}"
        );
    }

    /// The sign convention is stable: the dominant loading of each component is
    /// positive, so runs are comparable.
    #[test]
    fn component_signs_are_deterministic() {
        let (rows, cols) = (5, 5);
        let pattern: Vec<f64> = (0..rows * cols).map(|i| (i % 6) as f64).collect();
        // Negative scalings would naturally yield a negative-leaning PC1.
        let slices: Vec<Vec<f64>> = [-1.0, -2.0, -3.0]
            .iter()
            .map(|k| pattern.iter().map(|v| k * v).collect())
            .collect();
        let (_, outputs) = run(json!({ "input": cube_raster(cols, rows, &slices) }));
        let layer = load_input_layer(outputs["output_loadings"].as_str().unwrap()).unwrap();
        let pc1 = layer.schema.field_index("pc1").unwrap();
        let loads: Vec<f64> = layer
            .iter()
            .map(|f| match f.attributes[pc1] {
                FieldValue::Float(v) => v,
                _ => panic!(),
            })
            .collect();
        let dominant = loads
            .iter()
            .cloned()
            .fold(0.0f64, |acc, v| if v.abs() > acc.abs() { v } else { acc });
        assert!(
            dominant > 0.0,
            "the largest loading should be positive by convention, got {loads:?}"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            MultidimensionalPrincipalComponentsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({"input": "a.tif", "mode": "temporal"})).is_err());
        assert!(bad(json!({"input": "a.tif", "number_of_pc": 0})).is_err());
        assert!(bad(json!({"input": "a.tif", "correlation": true})).is_ok());
    }
}
