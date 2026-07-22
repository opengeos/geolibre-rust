//! GeoLibre tool: derive criterion weights from an AHP pairwise-comparison matrix.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Assign Weights By Pairwise Comparison*
//! (Spatial Analyst / Suitability Modeler). The catalog has a full suitability
//! stack — `fuzzy_overlay`, `calculate_composite_index`, bundled
//! `weighted_overlay`/`weighted_sum` — but no defensible way to *derive* the
//! weights those tools consume. The Analytic Hierarchy Process (AHP) is the
//! standard method and is absent from the ~791 bundled tool IDs.
//!
//! The input is an `n x n` reciprocal comparison matrix on Saaty's 1-9 scale
//! (`a[i][j]` = how many times more important criterion `i` is than `j`), given
//! as a JSON 2-D array (`matrix`) and/or a CSV table (`input`, labeled or plain).
//! The criterion weights are the normalized principal eigenvector of that matrix,
//! recovered by power iteration (no linear-algebra crate). The tool also reports
//! the consistency of the judgements: the principal eigenvalue `lambda_max`, the
//! consistency index `CI = (lambda_max - n)/(n - 1)`, and the consistency ratio
//! `CR = CI / RI[n]` against Saaty's random-index table — emitting a warning when
//! `CR > 0.1` (the accepted inconsistency threshold).
//!
//! Output is a geometry-less table `criterion, weight, rank` (or a CSV when the
//! path ends in `.csv`) that feeds straight into `calculate_composite_index` /
//! weighted overlay.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Feature, FieldDef, FieldType, FieldValue, Layer};

use crate::common::write_text_output;
use crate::vector_common::{parse_optional_str, write_or_store_layer};

pub struct PairwiseComparisonWeightsTool;

/// Saaty's random-index table (average CI of random reciprocal matrices),
/// indexed by matrix order `n`. Index 0 is a placeholder.
const RANDOM_INDEX: [f64; 11] = [
    0.0, 0.0, 0.0, 0.58, 0.90, 1.12, 1.24, 1.32, 1.41, 1.45, 1.49,
];

impl Tool for PairwiseComparisonWeightsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "pairwise_comparison_weights",
            display_name: "Pairwise Comparison Weights",
            summary: "Derive criterion weights from an AHP pairwise-comparison matrix (Saaty 1-9 scale) as the normalized principal eigenvector, with a consistency check (lambda_max, consistency index, consistency ratio) — like ArcGIS Assign Weights By Pairwise Comparison.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "matrix",
                    description: "The n x n reciprocal comparison matrix as a JSON 2-D array, e.g. [[1,3,5],[0.333,1,2],[0.2,0.5,1]]. Provide this or `input`.",
                    required: false,
                },
                ToolParamSpec {
                    name: "input",
                    description: "CSV file holding the comparison matrix (plain numeric rows, or labeled with a header row and a leading name column). Provide this or `matrix`.",
                    required: false,
                },
                ToolParamSpec {
                    name: "criteria",
                    description: "Comma-separated criterion names (must match the matrix order). Overrides CSV labels; defaults to C1..Cn.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output table path — a CSV (extension .csv) or a geometry-less vector table (criterion, weight, rank). If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        let inputs = ParsedInputs::from_args(args)?;
        // Building the model validates the matrix (square, positive, finite) and names.
        AhpModel::build(inputs.matrix, inputs.names)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let output = parse_optional_str(args, "output")?;
        let inputs = ParsedInputs::from_args(args)?;
        let model = AhpModel::build(inputs.matrix, inputs.names)?;

        let result = derive_weights(&model.matrix);
        let n = model.matrix.len();

        ctx.progress.info(&format!(
            "AHP over {n} criteria: lambda_max={:.4}, CR={:.4}",
            result.lambda_max, result.consistency_ratio
        ));

        // Rank criteria by descending weight (1 = most important).
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| result.weights[b].total_cmp(&result.weights[a]));
        let mut rank = vec![0i64; n];
        for (r, &idx) in order.iter().enumerate() {
            rank[idx] = (r + 1) as i64;
        }

        // ── Emit the weight table ─────────────────────────────────────────────
        let mut table = Layer::new("pairwise_comparison_weights");
        table.add_field(FieldDef::new("criterion", FieldType::Text));
        table.add_field(FieldDef::new("weight", FieldType::Float));
        table.add_field(FieldDef::new("rank", FieldType::Integer));

        let mut csv = String::from("criterion,weight,rank\n");
        for (i, name) in model.names.iter().enumerate() {
            table.push(Feature {
                fid: 0,
                geometry: None,
                attributes: vec![
                    FieldValue::Text(name.clone()),
                    FieldValue::Float(result.weights[i]),
                    FieldValue::Integer(rank[i]),
                ],
            });
            csv.push_str(&format!(
                "{},{},{}\n",
                escape_csv(name),
                result.weights[i],
                rank[i]
            ));
        }

        let out_path = match output {
            Some(path) if path.to_ascii_lowercase().ends_with(".csv") => {
                write_text_output(&csv, path)?;
                path.to_string()
            }
            other => write_or_store_layer(table, other)?,
        };

        let consistent = result.consistency_ratio <= 0.1;
        let warning = if consistent {
            None
        } else {
            Some(format!(
                "consistency ratio {:.4} exceeds 0.1 — the pairwise judgements are inconsistent; revise the comparisons",
                result.consistency_ratio
            ))
        };
        if let Some(w) = &warning {
            ctx.progress.info(w);
        }

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("n".to_string(), json!(n));
        outputs.insert("criteria".to_string(), json!(model.names));
        outputs.insert("weights".to_string(), json!(result.weights));
        outputs.insert("lambda_max".to_string(), json!(result.lambda_max));
        outputs.insert(
            "consistency_index".to_string(),
            json!(result.consistency_index),
        );
        outputs.insert(
            "consistency_ratio".to_string(),
            json!(result.consistency_ratio),
        );
        outputs.insert("random_index".to_string(), json!(result.random_index));
        outputs.insert("consistent".to_string(), json!(consistent));
        if let Some(w) = warning {
            outputs.insert("warning".to_string(), json!(w));
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── AHP computation ────────────────────────────────────────────────────────────

struct WeightResult {
    weights: Vec<f64>,
    lambda_max: f64,
    consistency_index: f64,
    consistency_ratio: f64,
    random_index: f64,
}

/// Derives the normalized principal eigenvector (criterion weights) of a
/// positive reciprocal matrix by power iteration, plus the consistency metrics.
fn derive_weights(a: &[Vec<f64>]) -> WeightResult {
    let n = a.len();

    // Power iteration for the principal eigenvector.
    let mut v = vec![1.0 / n as f64; n];
    for _ in 0..1000 {
        let mut w = mat_vec(a, &v);
        let sum: f64 = w.iter().sum();
        if sum > 0.0 {
            for x in &mut w {
                *x /= sum;
            }
        }
        let delta = w
            .iter()
            .zip(&v)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        v = w;
        if delta < 1e-15 {
            break;
        }
    }
    // Ensure the weights sum to exactly 1.
    let sum: f64 = v.iter().sum();
    if sum > 0.0 {
        for x in &mut v {
            *x /= sum;
        }
    }

    // lambda_max = mean over i of (A w)_i / w_i.
    let aw = mat_vec(a, &v);
    let lambda_max = if n == 0 {
        0.0
    } else {
        let s: f64 = (0..n)
            .map(|i| if v[i] > 0.0 { aw[i] / v[i] } else { 0.0 })
            .sum();
        s / n as f64
    };

    let consistency_index = if n > 1 {
        (lambda_max - n as f64) / (n as f64 - 1.0)
    } else {
        0.0
    };
    let random_index = if n < RANDOM_INDEX.len() {
        RANDOM_INDEX[n]
    } else {
        // Beyond the tabulated orders, RI plateaus; use the largest value.
        RANDOM_INDEX[RANDOM_INDEX.len() - 1]
    };
    let consistency_ratio = if random_index > 0.0 {
        (consistency_index / random_index).max(0.0)
    } else {
        0.0
    };

    WeightResult {
        weights: v,
        lambda_max,
        consistency_index,
        consistency_ratio,
        random_index,
    }
}

fn mat_vec(a: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    a.iter()
        .map(|row| row.iter().zip(v).map(|(x, y)| x * y).sum())
        .collect()
}

// ── Input parsing ──────────────────────────────────────────────────────────────

struct ParsedInputs {
    matrix: Vec<Vec<f64>>,
    names: Option<Vec<String>>,
}

impl ParsedInputs {
    fn from_args(args: &ToolArgs) -> Result<Self, ToolError> {
        let mut names_from_csv: Option<Vec<String>> = None;
        let matrix = match args.get("matrix") {
            Some(Value::Null) | None => match parse_optional_str(args, "input")? {
                Some(path) => {
                    let (m, names) = parse_csv_matrix(path)?;
                    names_from_csv = names;
                    m
                }
                None => {
                    return Err(ToolError::Validation(
                        "provide either 'matrix' (JSON 2-D array) or 'input' (CSV path)".into(),
                    ))
                }
            },
            Some(v) => parse_matrix_value(v)?,
        };

        // Explicit `criteria` overrides CSV-derived names.
        let names = match parse_optional_str(args, "criteria")? {
            Some(s) => Some(
                s.split(',')
                    .map(|t| t.trim().to_string())
                    .filter(|t| !t.is_empty())
                    .collect::<Vec<_>>(),
            ),
            None => names_from_csv,
        };

        Ok(ParsedInputs { matrix, names })
    }
}

/// Parses the `matrix` parameter, which may be a JSON array value or a JSON
/// string (host UIs post scalar params as strings).
fn parse_matrix_value(v: &Value) -> Result<Vec<Vec<f64>>, ToolError> {
    let value: Value = match v {
        Value::String(s) => serde_json::from_str(s.trim())
            .map_err(|e| ToolError::Validation(format!("'matrix' is not valid JSON: {e}")))?,
        other => other.clone(),
    };
    let rows = value
        .as_array()
        .ok_or_else(|| ToolError::Validation("'matrix' must be a JSON 2-D array".into()))?;
    let mut m = Vec::with_capacity(rows.len());
    for row in rows {
        let cells = row.as_array().ok_or_else(|| {
            ToolError::Validation("every row of 'matrix' must be a JSON array".into())
        })?;
        let mut r = Vec::with_capacity(cells.len());
        for c in cells {
            let x = c
                .as_f64()
                .ok_or_else(|| ToolError::Validation("'matrix' entries must be numbers".into()))?;
            r.push(x);
        }
        m.push(r);
    }
    Ok(m)
}

/// A parsed comparison matrix and, if the source was labeled, its criterion names.
type ParsedMatrix = (Vec<Vec<f64>>, Option<Vec<String>>);

/// Reads a CSV comparison matrix. Supports a plain numeric grid, or a labeled
/// grid with a header row (criterion names) and a leading name column.
fn parse_csv_matrix(path: &str) -> Result<ParsedMatrix, ToolError> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| ToolError::Execution(format!("cannot read '{path}': {e}")))?;
    let rows: Vec<Vec<String>> = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .map(|l| l.split(',').map(|c| c.trim().to_string()).collect())
        .collect();
    if rows.is_empty() {
        return Err(ToolError::Validation(format!("'{path}' is empty")));
    }

    // Plain numeric path: every cell parses and the grid is square.
    let all_numeric = rows
        .iter()
        .all(|r| r.iter().all(|c| c.parse::<f64>().is_ok()));
    let n = rows.len();
    if all_numeric && rows.iter().all(|r| r.len() == n) {
        let m = rows
            .iter()
            .map(|r| r.iter().map(|c| c.parse::<f64>().unwrap()).collect())
            .collect();
        return Ok((m, None));
    }

    // Labeled path: header row + leading name column.
    if rows.len() < 2 {
        return Err(ToolError::Validation(format!(
            "'{path}' is not a valid comparison matrix (need a square numeric grid or a labeled grid)"
        )));
    }
    let data = &rows[1..];
    let dn = data.len();
    let mut names = Vec::with_capacity(dn);
    let mut m = Vec::with_capacity(dn);
    for row in data {
        if row.len() != dn + 1 {
            return Err(ToolError::Validation(format!(
                "labeled row '{}' should have 1 label + {dn} values",
                row.first().cloned().unwrap_or_default()
            )));
        }
        names.push(row[0].clone());
        let mut r = Vec::with_capacity(dn);
        for c in &row[1..] {
            let x = c.parse::<f64>().map_err(|_| {
                ToolError::Validation(format!("non-numeric matrix cell '{c}' in '{path}'"))
            })?;
            r.push(x);
        }
        m.push(r);
    }
    Ok((m, Some(names)))
}

// ── Model / validation ─────────────────────────────────────────────────────────

struct AhpModel {
    matrix: Vec<Vec<f64>>,
    names: Vec<String>,
}

impl AhpModel {
    fn build(matrix: Vec<Vec<f64>>, names: Option<Vec<String>>) -> Result<Self, ToolError> {
        let n = matrix.len();
        if n == 0 {
            return Err(ToolError::Validation(
                "the comparison matrix is empty".into(),
            ));
        }
        for (i, row) in matrix.iter().enumerate() {
            if row.len() != n {
                return Err(ToolError::Validation(format!(
                    "the comparison matrix must be square: row {i} has {} entries, expected {n}",
                    row.len()
                )));
            }
            for &x in row {
                if !x.is_finite() || x <= 0.0 {
                    return Err(ToolError::Validation(
                        "matrix entries must be positive, finite numbers (Saaty scale)".into(),
                    ));
                }
            }
        }
        let names = match names {
            Some(ns) => {
                if ns.len() != n {
                    return Err(ToolError::Validation(format!(
                        "{} criterion name(s) given for a {n}x{n} matrix",
                        ns.len()
                    )));
                }
                ns
            }
            None => (1..=n).map(|i| format!("C{i}")).collect(),
        };
        Ok(AhpModel { matrix, names })
    }
}

fn escape_csv(s: &str) -> String {
    if s.contains([',', '"', '\n']) {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
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

    fn run(args: serde_json::Value) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        PairwiseComparisonWeightsTool.run(&args, &ctx()).unwrap()
    }

    /// A perfectly consistent matrix A[i][j] = w[i]/w[j] recovers w exactly,
    /// the weights sum to 1, and CR is ~0.
    #[test]
    fn recovers_consistent_weights() {
        let w = [0.5, 0.3, 0.2];
        let matrix: Vec<Vec<f64>> = (0..3)
            .map(|i| (0..3).map(|j| w[i] / w[j]).collect())
            .collect();
        let out = run(json!({ "matrix": matrix, "criteria": "a,b,c" }));
        let weights = out.outputs["weights"].as_array().unwrap();
        let got: Vec<f64> = weights.iter().map(|v| v.as_f64().unwrap()).collect();
        let sum: f64 = got.iter().sum();
        assert!((sum - 1.0).abs() < 1e-9, "weights must sum to 1, got {sum}");
        for (g, e) in got.iter().zip(w.iter()) {
            assert!((g - e).abs() < 1e-6, "weight {g} != expected {e}");
        }
        assert!(
            out.outputs["consistency_ratio"].as_f64().unwrap() < 1e-6,
            "a consistent matrix must have CR ~ 0"
        );
        assert_eq!(out.outputs["consistent"], json!(true));
        // lambda_max should equal n for a consistent matrix.
        assert!((out.outputs["lambda_max"].as_f64().unwrap() - 3.0).abs() < 1e-6);
    }

    /// The most-important criterion (largest weight) gets rank 1.
    #[test]
    fn ranks_by_weight() {
        // C2 clearly dominates.
        let matrix = json!([[1.0, 0.2, 0.5], [5.0, 1.0, 3.0], [2.0, 0.3333, 1.0]]);
        let out = run(json!({ "matrix": matrix }));
        let names = out.outputs["criteria"].as_array().unwrap();
        assert_eq!(names[1], json!("C2"));
        let weights: Vec<f64> = out.outputs["weights"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap())
            .collect();
        let top = (0..3)
            .max_by(|&a, &b| weights[a].total_cmp(&weights[b]))
            .unwrap();
        assert_eq!(top, 1, "C2 must carry the largest weight");
    }

    /// A cyclic, strongly inconsistent matrix triggers the CR > 0.1 warning.
    #[test]
    fn flags_inconsistent_matrix() {
        let matrix = json!([[1.0, 0.3333, 3.0], [3.0, 1.0, 0.3333], [0.3333, 3.0, 1.0]]);
        let out = run(json!({ "matrix": matrix }));
        assert!(
            out.outputs["consistency_ratio"].as_f64().unwrap() > 0.1,
            "cyclic judgements should be inconsistent"
        );
        assert_eq!(out.outputs["consistent"], json!(false));
        assert!(out.outputs.contains_key("warning"));
    }

    /// A labeled CSV matrix is parsed, and its row labels become criterion names.
    #[test]
    fn parses_labeled_csv() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("ahp_test_{}.csv", std::process::id()));
        let csv = ",cost,slope,access\ncost,1,3,5\nslope,0.3333,1,2\naccess,0.2,0.5,1\n";
        std::fs::write(&path, csv).unwrap();
        let out = run(json!({ "input": path.to_str().unwrap() }));
        assert_eq!(out.outputs["n"], json!(3));
        assert_eq!(out.outputs["criteria"], json!(["cost", "slope", "access"]));
        // cost should be the most important.
        let weights: Vec<f64> = out.outputs["weights"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap())
            .collect();
        assert!(weights[0] > weights[1] && weights[1] > weights[2]);
        std::fs::remove_file(&path).ok();
    }

    /// Single-criterion matrix: weight 1, CR 0.
    #[test]
    fn single_criterion() {
        let out = run(json!({ "matrix": [[1.0]], "criteria": "only" }));
        assert_eq!(out.outputs["weights"], json!([1.0]));
        assert_eq!(out.outputs["consistency_ratio"], json!(0.0));
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            PairwiseComparisonWeightsTool.validate(&args)
        };
        // No matrix and no input.
        assert!(bad(json!({})).is_err());
        // Non-square matrix.
        assert!(bad(json!({ "matrix": [[1.0, 2.0], [0.5, 1.0, 3.0]] })).is_err());
        // Non-positive entry.
        assert!(bad(json!({ "matrix": [[1.0, -2.0], [0.5, 1.0]] })).is_err());
        // Wrong number of criterion names.
        assert!(bad(json!({ "matrix": [[1.0, 2.0], [0.5, 1.0]], "criteria": "a,b,c" })).is_err());
        // Valid.
        assert!(bad(json!({ "matrix": [[1.0, 2.0], [0.5, 1.0]] })).is_ok());
    }
}
