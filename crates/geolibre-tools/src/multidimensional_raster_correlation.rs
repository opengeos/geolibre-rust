//! GeoLibre tool: per-pixel correlation between two raster cubes.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Multidimensional Raster Correlation*
//! (Image Analyst).
//!
//! ## Why the catalog needs it
//!
//! "Where does rainfall drive greenness, and with what delay?" is a per-pixel
//! question: the answer is a *map* of correlation coefficients, not one number
//! for the scene. Together with the lag map it is how a teleconnection, a
//! drought response time, or a lake's reaction to snowmelt is actually located.
//!
//! Nothing in either registry produces that map. `image_correlation` and
//! `image_autocorrelation` collapse a whole image pair to a single scalar;
//! `time_series_cross_correlation` (round 17) works on *vector* features;
//! `attribute_correlation` is a table operation; `bivariate_spatial_association`
//! compares two single-band rasters, so it has no time axis to correlate along.
//!
//! ## Method
//!
//! For every cell the two cubes' series are paired slice by slice and reduced
//! to one coefficient:
//!
//! * **Pearson** — linear association.
//! * **Spearman** — rank association, so a monotone but curved relationship
//!   still scores 1 and outliers do not dominate.
//! * **Kendall** — concordant minus discordant pairs (tau-b, tie-corrected);
//!   the most robust and the most expensive, at O(n^2) in the slice count.
//!
//! A fixed `lag` shifts the second cube forward before pairing. With
//! `cross_correlation` the tool instead searches every lag in
//! `-max_lag ..= max_lag` and writes both the strongest coefficient and the lag
//! at which it occurred — which is the delay itself, the thing the analysis is
//! usually after.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::args_common::{bool_or, choice_or, opt_usize, usize_or};
use crate::common::{parse_optional_output, raster_like_with_data, write_or_store_output};
use crate::cube::load_cube;
use crate::raster_stack::check_alignment_refs;

pub struct MultidimensionalRasterCorrelationTool;

impl Tool for MultidimensionalRasterCorrelationTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "multidimensional_raster_correlation",
            display_name: "Multidimensional Raster Correlation",
            summary: "Correlates two raster cubes cell by cell along their shared dimension, producing a map of Pearson, Spearman or Kendall coefficients and — in cross-correlation mode — the lag at which each cell's association peaks (ArcGIS Multidimensional Raster Correlation). Nothing in either registry maps this: image_correlation collapses an image pair to one scalar, round 17's time_series_cross_correlation is vector-only, attribute_correlation is a table operation, and bivariate_spatial_association compares two single-band rasters with no time axis.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input1",
                    description: "First cube: one multiband raster (each band a slice) or a comma-separated list of co-registered rasters.",
                    required: true,
                },
                ToolParamSpec {
                    name: "input2",
                    description: "Second cube, co-registered with the first and with the same number of slices.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output correlation raster in [-1, 1]. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_lag",
                    description: "Output raster of the best lag per cell (cross-correlation mode); the fixed lag everywhere otherwise. Always produced; stored in memory when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_count",
                    description: "Output raster of the paired-sample count behind each coefficient. Always produced; stored in memory when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'pearson' (default), 'spearman', or 'kendall'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "lag",
                    description: "Fixed number of slices to shift the second cube forward before pairing (default 0). Ignored in cross-correlation mode.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cross_correlation",
                    description: "Search every lag in -max_lag..=max_lag and keep the strongest coefficient (default false).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_lag",
                    description: "Half-width of the lag search in slices (default 3). Cross-correlation mode only.",
                    required: false,
                },
                ToolParamSpec {
                    name: "use_absolute",
                    description: "In cross-correlation mode, rank lags by absolute correlation so a strong negative association wins (default true).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_valid",
                    description: "Minimum paired samples before a coefficient is written (default 3).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        crate::raster_stack::parse_input_paths(args, "input1")?;
        crate::raster_stack::parse_input_paths(args, "input2")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let prm = parse_params(args)?;
        let output = parse_optional_output(args, "output")?;
        let out_lag = parse_optional_output(args, "output_lag")?;
        let out_count = parse_optional_output(args, "output_count")?;

        let a = load_cube(args, "input1", "dimension_values", "dimension", 2)?;
        let b = load_cube(args, "input2", "dimension_values", "dimension", 2)?;
        check_alignment_refs(&[a.template(), b.template()])?;
        if a.len() != b.len() {
            return Err(ToolError::Validation(format!(
                "the cubes have different slice counts ({} and {}); they must be paired slice by \
                 slice",
                a.len(),
                b.len()
            )));
        }

        let (rows, cols, n) = (a.rows, a.cols, a.len());
        let lags: Vec<isize> = if prm.cross {
            let m = prm.max_lag as isize;
            (-m..=m).collect()
        } else {
            vec![prm.lag]
        };
        if lags.iter().all(|&l| l.unsigned_abs() >= n) {
            return Err(ToolError::Validation(format!(
                "every requested lag is at least the slice count ({n}); no pairs would remain"
            )));
        }

        ctx.progress.info(&format!(
            "{rows}x{cols}, {n} slice(s), {} correlation over {} lag(s)",
            prm.method.label(),
            lags.len()
        ));

        let nodata = -9999.0_f64;
        let mut corr = vec![nodata; rows * cols];
        let mut lag_out = vec![nodata; rows * cols];
        let mut count_out = vec![nodata; rows * cols];

        let mut sa: Vec<Option<f64>> = Vec::with_capacity(n);
        let mut sb: Vec<Option<f64>> = Vec::with_capacity(n);
        let mut xs: Vec<f64> = Vec::with_capacity(n);
        let mut ys: Vec<f64> = Vec::with_capacity(n);

        for r in 0..rows {
            for c in 0..cols {
                a.series(r, c, &mut sa);
                b.series(r, c, &mut sb);

                let mut best: Option<(f64, isize, usize)> = None;
                for &lag in &lags {
                    xs.clear();
                    ys.clear();
                    // Pair slice i of cube A with slice i+lag of cube B.
                    for i in 0..n {
                        let j = i as isize + lag;
                        if j < 0 || j as usize >= n {
                            continue;
                        }
                        if let (Some(x), Some(y)) = (sa[i], sb[j as usize]) {
                            xs.push(x);
                            ys.push(y);
                        }
                    }
                    if xs.len() < prm.min_valid {
                        continue;
                    }
                    let Some(v) = correlate(&mut xs, &mut ys, prm.method) else {
                        continue;
                    };
                    let score = if prm.cross && prm.use_absolute {
                        v.abs()
                    } else {
                        v
                    };
                    let better = match best {
                        None => true,
                        Some((bv, blag, _)) => {
                            let bscore = if prm.cross && prm.use_absolute {
                                bv.abs()
                            } else {
                                bv
                            };
                            // Ties resolve to the smaller |lag|. Lags are
                            // searched from -max_lag upwards, so taking the
                            // first winner would systematically report the most
                            // negative of several equally good lags.
                            score > bscore
                                || (score == bscore && lag.abs() < blag.abs())
                        }
                    };
                    if better {
                        best = Some((v, lag, xs.len()));
                    }
                }

                if let Some((v, lag, count)) = best {
                    let i = r * cols + c;
                    corr[i] = v;
                    lag_out[i] = lag as f64;
                    count_out[i] = count as f64;
                }
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let template = a.template();
        let out_path = write_or_store_output(
            raster_like_with_data(template, corr, nodata, DataType::F32)?,
            output,
        )?;
        let lag_path = write_or_store_output(
            raster_like_with_data(template, lag_out, nodata, DataType::F32)?,
            out_lag,
        )?;
        let count_path = write_or_store_output(
            raster_like_with_data(template, count_out, nodata, DataType::F32)?,
            out_count,
        )?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("output_lag".to_string(), json!(lag_path));
        outputs.insert("output_count".to_string(), json!(count_path));
        outputs.insert("method".to_string(), json!(prm.method.label()));
        outputs.insert("slices".to_string(), json!(n));
        outputs.insert("cross_correlation".to_string(), json!(prm.cross));
        outputs.insert("lags_searched".to_string(), json!(lags.len()));
        Ok(ToolRunResult { outputs })
    }
}

/// One correlation coefficient. `xs`/`ys` may be reordered.
fn correlate(xs: &mut [f64], ys: &mut [f64], method: Method) -> Option<f64> {
    match method {
        Method::Pearson => pearson(xs, ys),
        Method::Spearman => {
            // Spearman is Pearson on the ranks; average ranks handle ties.
            let rx = ranks(xs);
            let ry = ranks(ys);
            pearson(&rx, &ry)
        }
        Method::Kendall => kendall_tau_b(xs, ys),
    }
}

fn pearson(xs: &[f64], ys: &[f64]) -> Option<f64> {
    let n = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
    for (x, y) in xs.iter().zip(ys) {
        let (dx, dy) = (x - mx, y - my);
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    // A constant series has no variance, so no correlation is defined; report
    // no-data rather than 0, which would read as "measured, and unrelated".
    if sxx <= 0.0 || syy <= 0.0 {
        return None;
    }
    Some((sxy / (sxx * syy).sqrt()).clamp(-1.0, 1.0))
}

/// Fractional ranks, averaging over ties.
fn ranks(v: &[f64]) -> Vec<f64> {
    let n = v.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| v[a].total_cmp(&v[b]));
    let mut out = vec![0.0f64; n];
    let mut i = 0usize;
    while i < n {
        let mut j = i + 1;
        while j < n && v[idx[j]] == v[idx[i]] {
            j += 1;
        }
        // Mean of the 1-based ranks spanned by this tie group.
        let mean_rank = ((i + 1 + j) as f64) / 2.0;
        for &k in &idx[i..j] {
            out[k] = mean_rank;
        }
        i = j;
    }
    out
}

/// Kendall's tau-b, which corrects for ties in either series.
fn kendall_tau_b(xs: &[f64], ys: &[f64]) -> Option<f64> {
    let n = xs.len();
    let (mut concordant, mut discordant) = (0i64, 0i64);
    let (mut tied_x, mut tied_y) = (0i64, 0i64);
    for i in 0..n {
        for j in i + 1..n {
            let dx = xs[i] - xs[j];
            let dy = ys[i] - ys[j];
            if dx == 0.0 && dy == 0.0 {
                // Tied in both: counted in neither denominator term.
                continue;
            }
            if dx == 0.0 {
                tied_x += 1;
            } else if dy == 0.0 {
                tied_y += 1;
            } else if (dx > 0.0) == (dy > 0.0) {
                concordant += 1;
            } else {
                discordant += 1;
            }
        }
    }
    let n0 = concordant + discordant;
    let denom = (((n0 + tied_x) as f64) * ((n0 + tied_y) as f64)).sqrt();
    if denom <= 0.0 {
        return None;
    }
    Some((((concordant - discordant) as f64) / denom).clamp(-1.0, 1.0))
}

// ── Parameters ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Pearson,
    Spearman,
    Kendall,
}

impl Method {
    fn label(self) -> &'static str {
        match self {
            Method::Pearson => "pearson",
            Method::Spearman => "spearman",
            Method::Kendall => "kendall",
        }
    }
}

struct Params {
    method: Method,
    lag: isize,
    cross: bool,
    max_lag: usize,
    use_absolute: bool,
    min_valid: usize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let method = match choice_or(
        args,
        "method",
        &["pearson", "spearman", "kendall"],
        "pearson",
    )? {
        "spearman" => Method::Spearman,
        "kendall" => Method::Kendall,
        _ => Method::Pearson,
    };

    // A lag is signed, so it cannot go through the usize parsers.
    let lag = match args.get("lag") {
        None | Some(Value::Null) => 0isize,
        Some(Value::Number(n)) => n
            .as_i64()
            .ok_or_else(|| ToolError::Validation("'lag' must be a whole number".to_string()))?
            as isize,
        Some(Value::String(s)) if s.trim().is_empty() => 0,
        Some(Value::String(s)) => s
            .trim()
            .parse::<isize>()
            .map_err(|_| ToolError::Validation("'lag' must be a whole number".to_string()))?,
        Some(_) => {
            return Err(ToolError::Validation(
                "'lag' must be a whole number".to_string(),
            ))
        }
    };

    let cross = bool_or(args, "cross_correlation", false)?;
    let max_lag = usize_or(args, "max_lag", 3)?;
    if cross && max_lag == 0 {
        return Err(ToolError::Validation(
            "'max_lag' must be at least 1 in cross-correlation mode".to_string(),
        ));
    }
    let use_absolute = bool_or(args, "use_absolute", true)?;

    let min_valid = match opt_usize(args, "min_valid")? {
        None => 3,
        Some(v) if v >= 2 => v,
        Some(v) => {
            return Err(ToolError::Validation(format!(
                "'min_valid' must be at least 2; a correlation needs two paired samples, got {v}"
            )))
        }
    };

    Ok(Params {
        method,
        lag,
        cross,
        max_lag,
        use_absolute,
        min_valid,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::load_input_raster;
    use crate::cube::test_support::cube_raster;
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

    fn series_cube(vals: &[f64]) -> String {
        let slices: Vec<Vec<f64>> = vals.iter().map(|&v| vec![v]).collect();
        cube_raster(1, 1, &slices)
    }

    fn run(args: Value) -> (Raster, BTreeMap<String, Value>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = MultidimensionalRasterCorrelationTool.run(&args, &ctx()).unwrap();
        let r = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (r, out.outputs)
    }

    /// Perfect positive and negative linear relationships hit exactly +/-1.
    #[test]
    fn pearson_extremes() {
        let up = series_cube(&[1.0, 2.0, 3.0, 4.0, 5.0]);
        let same = series_cube(&[2.0, 4.0, 6.0, 8.0, 10.0]);
        let down = series_cube(&[5.0, 4.0, 3.0, 2.0, 1.0]);
        assert!((run(json!({"input1": up.clone(), "input2": same})).0.get(0, 0, 0) - 1.0).abs() < 1e-6);
        assert!((run(json!({"input1": up, "input2": down})).0.get(0, 0, 0) + 1.0).abs() < 1e-6);
    }

    /// Spearman scores a monotone but curved relationship 1 where Pearson does
    /// not — the reason both are offered.
    #[test]
    fn spearman_handles_a_monotone_curve() {
        let x = series_cube(&[1.0, 2.0, 3.0, 4.0, 5.0]);
        let y = series_cube(&[1.0, 4.0, 9.0, 16.0, 250.0]);
        let p = run(json!({"input1": x.clone(), "input2": y.clone()}))
            .0
            .get(0, 0, 0);
        let s = run(json!({"input1": x.clone(), "input2": y.clone(), "method": "spearman"}))
            .0
            .get(0, 0, 0);
        let k = run(json!({"input1": x, "input2": y, "method": "kendall"}))
            .0
            .get(0, 0, 0);
        assert!((s - 1.0).abs() < 1e-6, "spearman should be exactly 1, got {s}");
        assert!((k - 1.0).abs() < 1e-6, "kendall should be exactly 1, got {k}");
        assert!(
            p < 0.95,
            "pearson should be dragged down by the curvature, got {p}"
        );
    }

    /// The headline capability: recovering the delay between two series.
    #[test]
    fn cross_correlation_recovers_the_lag() {
        // y is x delayed by two slices: y[i] = x[i-2]. The tool pairs A[i]
        // with B[i+lag], so the match is at lag = +2 — B's copy of the event
        // A saw at slice i turns up two slices later.
        let x = series_cube(&[1.0, 5.0, 2.0, 8.0, 3.0, 9.0, 4.0, 7.0]);
        let y = series_cube(&[0.0, 0.0, 1.0, 5.0, 2.0, 8.0, 3.0, 9.0]);
        let args: ToolArgs = serde_json::from_value(json!({
            "input1": x, "input2": y, "cross_correlation": true, "max_lag": 3, "min_valid": 4
        }))
        .unwrap();
        let out = MultidimensionalRasterCorrelationTool.run(&args, &ctx()).unwrap();
        let corr = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        let lag = load_input_raster(out.outputs["output_lag"].as_str().unwrap()).unwrap();
        assert_eq!(
            lag.get(0, 0, 0),
            2.0,
            "a two-slice delay should show up as lag +2"
        );
        assert!(
            corr.get(0, 0, 0) > 0.99,
            "at the right lag the series match exactly, got {}",
            corr.get(0, 0, 0)
        );
    }

    /// A fixed lag is applied verbatim rather than searched.
    #[test]
    fn fixed_lag_is_used_as_given() {
        let x = series_cube(&[1.0, 5.0, 2.0, 8.0, 3.0, 9.0]);
        let y = series_cube(&[0.0, 1.0, 5.0, 2.0, 8.0, 3.0]);
        // y is x delayed by one slice, so the alignment is at lag +1.
        let (out, outputs) = run(json!({
            "input1": x, "input2": y, "lag": 1, "min_valid": 4
        }));
        assert_eq!(outputs["lags_searched"].as_u64().unwrap(), 1);
        assert!(
            out.get(0, 0, 0) > 0.99,
            "shifting by one slice should align the series exactly"
        );
    }

    /// The output is a map: different cells get different answers.
    #[test]
    fn correlation_varies_per_cell() {
        // Cell 0 correlates positively, cell 1 negatively.
        let a = cube_raster(2, 1, &[vec![1.0, 1.0], vec![2.0, 2.0], vec![3.0, 3.0]]);
        let b = cube_raster(2, 1, &[vec![10.0, 30.0], vec![20.0, 20.0], vec![30.0, 10.0]]);
        let (out, _) = run(json!({"input1": a, "input2": b}));
        assert!((out.get(0, 0, 0) - 1.0).abs() < 1e-6);
        assert!((out.get(0, 0, 1) + 1.0).abs() < 1e-6);
    }

    /// A constant series has no variance, so no coefficient exists — report
    /// no-data rather than 0, which would read as "measured, and unrelated".
    #[test]
    fn constant_series_is_nodata_not_zero() {
        let x = series_cube(&[1.0, 2.0, 3.0, 4.0]);
        let flat = series_cube(&[7.0, 7.0, 7.0, 7.0]);
        let (out, _) = run(json!({"input1": x, "input2": flat}));
        assert_eq!(out.get(0, 0, 0), -9999.0);
    }

    /// No-data slices are dropped from the pairing and the surviving count is
    /// reported.
    #[test]
    fn pairs_skip_nodata_and_report_the_count() {
        let x = series_cube(&[1.0, 2.0, 3.0, 4.0, 5.0]);
        let y = series_cube(&[2.0, -9999.0, 6.0, 8.0, 10.0]);
        let args: ToolArgs =
            serde_json::from_value(json!({"input1": x, "input2": y})).unwrap();
        let out = MultidimensionalRasterCorrelationTool.run(&args, &ctx()).unwrap();
        let corr = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        let count = load_input_raster(out.outputs["output_count"].as_str().unwrap()).unwrap();
        assert_eq!(count.get(0, 0, 0), 4.0, "one pair should be dropped");
        assert!((corr.get(0, 0, 0) - 1.0).abs() < 1e-6);
    }

    /// `min_valid` suppresses coefficients built from too few pairs.
    #[test]
    fn min_valid_suppresses_thin_series() {
        let x = series_cube(&[1.0, 2.0, 3.0]);
        let y = series_cube(&[2.0, 4.0, 6.0]);
        let (out, _) = run(json!({"input1": x, "input2": y, "min_valid": 5}));
        assert_eq!(out.get(0, 0, 0), -9999.0);
    }

    /// Mismatched cube lengths are an error, not a silent truncation.
    #[test]
    fn mismatched_slice_counts_are_rejected() {
        let x = series_cube(&[1.0, 2.0, 3.0]);
        let y = series_cube(&[1.0, 2.0]);
        let args: ToolArgs =
            serde_json::from_value(json!({"input1": x, "input2": y})).unwrap();
        let err = MultidimensionalRasterCorrelationTool
            .run(&args, &ctx())
            .unwrap_err();
        assert!(
            format!("{err:?}").contains("slice counts"),
            "expected a slice-count error, got {err:?}"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            MultidimensionalRasterCorrelationTool.validate(&args)
        };
        assert!(bad(json!({"input2": "b.tif"})).is_err());
        assert!(bad(json!({"input1": "a.tif"})).is_err());
        let base = json!({"input1": "a.tif", "input2": "b.tif"});
        assert!(bad(base.clone()).is_ok());
        let with = |k: &str, v: Value| {
            let mut m = base.as_object().unwrap().clone();
            m.insert(k.into(), v);
            Value::Object(m)
        };
        assert!(bad(with("method", json!("kruskal"))).is_err());
        assert!(bad(with("lag", json!("x"))).is_err());
        assert!(bad(with("min_valid", json!(1))).is_err());
        assert!(bad(with("method", json!("kendall"))).is_ok());
        assert!(bad(with("lag", json!(-2))).is_ok());
    }
}
