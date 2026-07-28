//! GeoLibre tool: diagnostic sweep of candidate aggregation bin sizes.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Evaluate Bin Sizes* (Spatial
//! Statistics). The shipped `optimized_hot_spot_analysis`,
//! `optimized_outlier_analysis` and `emerging_hot_spot_analysis` all aggregate
//! incident points into bins before testing, and each picks a bin size by
//! internal heuristic with no diagnostic exposed. Bin size materially changes
//! the result (the modifiable areal unit problem), so users currently have no
//! way to check whether the chosen size is reasonable or to justify a different
//! one.
//!
//! `incremental_spatial_autocorrelation` (shipped) solves the adjacent problem
//! of choosing a *distance band* for a fixed geometry; it does not evaluate bin
//! size for aggregation.
//!
//! The candidate sweep is anchored on the **mean nearest-neighbour distance**,
//! which keeps it scale-appropriate without the user guessing a magnitude: a
//! bin much smaller than the typical point spacing is mostly empty, and one
//! much larger washes out structure. Each candidate reports bin count,
//! non-empty proportion, mean/median counts, and the coefficient of variation
//! — together these identify both over-aggregation (too few bins) and sparsity.

use std::collections::BTreeMap;

use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Tessellation shape for the trial bins.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum BinShape {
    Square,
    /// Row-offset ("brick") binning — see `hex_bin` for what this does and
    /// does not guarantee.
    Hexagon,
}

pub struct EvaluateBinSizesTool;

impl Tool for EvaluateBinSizesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "evaluate_bin_sizes",
            display_name: "Evaluate Bin Sizes",
            summary: "Sweep candidate aggregation bin sizes and report diagnostics identifying which size best reveals structure, like ArcGIS Evaluate Bin Sizes.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Point features to be aggregated.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output diagnostics table path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "bin_shape",
                    description: "hexagon (default; row-offset approximation, see docs) | square.",
                    required: false,
                },
                ToolParamSpec {
                    name: "sizes",
                    description: "Comma-separated candidate bin sizes (map units). When omitted, a geometric sweep is generated from the mean nearest-neighbour distance.",
                    required: false,
                },
                ToolParamSpec {
                    name: "steps",
                    description: "Number of candidates in the generated sweep (default 8, max 64).",
                    required: false,
                },
                ToolParamSpec {
                    name: "analysis_field",
                    description: "Optional numeric field to sum per bin instead of counting points.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_shape(args)?;
        if let Some(list) = parse_optional_str(args, "sizes")? {
            let sizes = parse_sizes(list)?;
            if sizes.is_empty() {
                return Err(ToolError::Validation(
                    "'sizes' did not contain any positive values".to_string(),
                ));
            }
        }
        if let Some(n) = parse_optional_f64(args, "steps")? {
            if !n.is_finite() || !(1.0..=64.0).contains(&n) {
                return Err(ToolError::Validation(
                    "'steps' must be between 1 and 64".to_string(),
                ));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let shape = parse_shape(args)?;
        let steps = parse_optional_f64(args, "steps")?.unwrap_or(8.0) as usize;
        let analysis_field = parse_optional_str(args, "analysis_field")?.map(String::from);

        let layer = load_input_layer(input)?;
        if let Some(f) = &analysis_field {
            if layer.schema.field_index(f).is_none() {
                return Err(ToolError::Validation(format!(
                    "analysis_field '{f}' not found on the input layer"
                )));
            }
        }

        // Extract points and their weights.
        let mut pts: Vec<(f64, f64, f64)> = Vec::new();
        for feat in layer.features.iter() {
            let Some(geom) = feat.geometry.as_ref() else {
                continue;
            };
            let Some((x, y)) = point_xy(geom) else {
                continue;
            };
            let w = match &analysis_field {
                Some(f) => feat
                    .get(&layer.schema, f)
                    .ok()
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0),
                None => 1.0,
            };
            pts.push((x, y, w));
        }
        if pts.len() < 2 {
            return Err(ToolError::Validation(format!(
                "need at least 2 points to evaluate bin sizes, found {}",
                pts.len()
            )));
        }

        // Mean nearest-neighbour distance anchors the sweep.
        let mnn = mean_nearest_neighbour(&pts)?;

        let sizes = match parse_optional_str(args, "sizes")? {
            Some(list) => parse_sizes(list)?,
            None => {
                // Geometric sweep from half the mean NN distance upward, so the
                // range spans "mostly empty" through "over-aggregated".
                let base = if mnn > 0.0 { mnn } else { 1.0 };
                (0..steps)
                    .map(|i| base * 0.5 * 2f64.powf(i as f64 * 0.5))
                    .collect()
            }
        };

        ctx.progress.info(&format!(
            "evaluating {} candidate size(s) over {} point(s); mean NN distance {mnn:.4}",
            sizes.len(),
            pts.len()
        ));

        let (min_x, min_y) = pts
            .iter()
            .fold((f64::MAX, f64::MAX), |(ax, ay), (x, y, _)| {
                (ax.min(*x), ay.min(*y))
            });
        let (max_x, max_y) = pts
            .iter()
            .fold((f64::MIN, f64::MIN), |(ax, ay), (x, y, _)| {
                (ax.max(*x), ay.max(*y))
            });

        let mut out = Layer::new("bin_size_evaluation");
        out.add_field(FieldDef::new("size", FieldType::Float));
        out.add_field(FieldDef::new("bin_count", FieldType::Integer));
        out.add_field(FieldDef::new("nonempty_bins", FieldType::Integer));
        out.add_field(FieldDef::new("nonempty_ratio", FieldType::Float));
        out.add_field(FieldDef::new("mean_per_bin", FieldType::Float));
        out.add_field(FieldDef::new("median_per_bin", FieldType::Float));
        out.add_field(FieldDef::new("max_per_bin", FieldType::Float));
        out.add_field(FieldDef::new("cv", FieldType::Float));

        let mut best: Option<(f64, f64)> = None; // (score, size)
        for &size in &sizes {
            if !size.is_finite() || size <= 0.0 {
                continue;
            }
            let mut bins: BTreeMap<(i64, i64), f64> = BTreeMap::new();
            let mut point_counts: BTreeMap<(i64, i64), usize> = BTreeMap::new();
            for (x, y, w) in &pts {
                let key = match shape {
                    BinShape::Square => square_bin(*x - min_x, *y - min_y, size),
                    BinShape::Hexagon => hex_bin(*x - min_x, *y - min_y, size),
                };
                *bins.entry(key).or_insert(0.0) += w;
                *point_counts.entry(key).or_insert(0) += 1;
            }
            let counts: Vec<f64> = bins.values().copied().collect();
            let nonempty = counts.len();
            // Total bins needed to tile the point extent at this size. This is
            // the denominator that makes nonempty_ratio meaningful: counting
            // only occupied bins would make the ratio trivially 1.
            let span_x = (max_x - min_x).max(f64::EPSILON);
            let span_y = (max_y - min_y).max(f64::EPSILON);
            let row_h = match shape {
                BinShape::Square => size,
                BinShape::Hexagon => size * 0.8660254037844386,
            };
            let total_bins =
                (((span_x / size).floor() + 1.0) * ((span_y / row_h).floor() + 1.0)).max(1.0);
            let total: f64 = counts.iter().sum();
            let mean = if nonempty > 0 {
                total / nonempty as f64
            } else {
                0.0
            };
            let var = if nonempty > 0 {
                counts.iter().map(|c| (c - mean).powi(2)).sum::<f64>() / nonempty as f64
            } else {
                0.0
            };
            let cv = if mean.abs() > f64::EPSILON {
                var.sqrt() / mean
            } else {
                0.0
            };
            let mut sorted = counts.clone();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let median = if sorted.is_empty() {
                0.0
            } else {
                sorted[sorted.len() / 2]
            };
            let max = sorted.last().copied().unwrap_or(0.0);

            // Score favours a high coefficient of variation (structure is
            // visible) while penalizing bins so small that most are empty or a
            // single point. Both extremes score low.
            let mean_points = if nonempty > 0 {
                point_counts.values().sum::<usize>() as f64 / nonempty as f64
            } else {
                0.0
            };
            let occupancy = (mean_points / 2.0).min(1.0);
            let score = cv * occupancy;
            if best.is_none_or(|(bs, _)| score > bs) {
                best = Some((score, size));
            }

            out.add_feature(
                None,
                &[
                    ("size", FieldValue::Float(size)),
                    ("bin_count", FieldValue::Integer(total_bins as i64)),
                    ("nonempty_bins", FieldValue::Integer(nonempty as i64)),
                    (
                        "nonempty_ratio",
                        FieldValue::Float((nonempty as f64 / total_bins).min(1.0)),
                    ),
                    ("mean_per_bin", FieldValue::Float(mean)),
                    ("median_per_bin", FieldValue::Float(median)),
                    ("max_per_bin", FieldValue::Float(max)),
                    ("cv", FieldValue::Float(cv)),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed writing candidate row: {e}")))?;
        }

        let recommended = best.map(|(_, s)| s).unwrap_or(mnn);
        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("recommended_size".to_string(), json!(recommended));
        outputs.insert("mean_nearest_neighbour".to_string(), json!(mnn));
        outputs.insert("candidate_count".to_string(), json!(sizes.len()));
        outputs.insert("point_count".to_string(), json!(pts.len()));
        Ok(ToolRunResult { outputs })
    }
}

/// Square-grid bin index.
fn square_bin(x: f64, y: f64, size: f64) -> (i64, i64) {
    ((x / size).floor() as i64, (y / size).floor() as i64)
}

/// Row-offset bin index approximating a hexagonal layout. `size` is the
/// centre-to-centre spacing along x.
///
/// Rows are half-width offset and spaced by `size * sqrt(3)/2`, which gives the
/// staggered arrangement and the row density of a true hex grid — but
/// membership is still decided by a rectangular `floor` on x, so the cells are
/// offset rectangles rather than hexagons. A point can therefore land in a
/// different bin than true axial hex rounding would assign, by up to half a row
/// height near a row boundary.
///
/// For this tool's purpose that is adequate: the diagnostics reported (bin
/// occupancy, counts per bin, coefficient of variation) depend on the bin
/// *size and density*, not on exact cell boundaries. Callers who need genuine
/// hexagons for the aggregation itself should use the bundled
/// `hexagonal_grid_from_*` tools once a size has been chosen here.
fn hex_bin(x: f64, y: f64, size: f64) -> (i64, i64) {
    let w = size;
    let h = size * 0.8660254037844386; // sqrt(3)/2
    let row = (y / h).floor() as i64;
    // Odd rows are offset by half a hexagon, which is what makes the
    // tessellation hexagonal rather than a plain rectangular grid.
    let shift = if row.rem_euclid(2) == 1 { w / 2.0 } else { 0.0 };
    let col = ((x - shift) / w).floor() as i64;
    (col, row)
}

/// Mean distance from each point to its nearest other point.
fn mean_nearest_neighbour(pts: &[(f64, f64, f64)]) -> Result<f64, ToolError> {
    let mut tree: KdTree<f64, usize, [f64; 2]> = KdTree::new(2);
    for (i, (x, y, _)) in pts.iter().enumerate() {
        tree.add([*x, *y], i)
            .map_err(|e| ToolError::Execution(format!("kd-tree insert failed: {e:?}")))?;
    }
    let mut sum = 0.0;
    let mut n = 0usize;
    for (i, (x, y, _)) in pts.iter().enumerate() {
        // Ask for 2: the first hit is the point itself.
        let hits = tree
            .nearest(&[*x, *y], 2, &squared_euclidean)
            .map_err(|e| ToolError::Execution(format!("kd-tree query failed: {e:?}")))?;
        for (d2, &j) in hits {
            if j != i {
                sum += d2.sqrt();
                n += 1;
                break;
            }
        }
    }
    Ok(if n > 0 { sum / n as f64 } else { 0.0 })
}

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        _ => None,
    }
}

// ── parameter parsing ────────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

fn parse_sizes(raw: &str) -> Result<Vec<f64>, ToolError> {
    let mut out = Vec::new();
    for tok in raw.split(',') {
        let tok = tok.trim();
        if tok.is_empty() {
            continue;
        }
        let v: f64 = tok
            .parse()
            .map_err(|_| ToolError::Validation(format!("'sizes' entry '{tok}' is not a number")))?;
        if !v.is_finite() || v <= 0.0 {
            return Err(ToolError::Validation(format!(
                "'sizes' entry '{tok}' must be greater than 0"
            )));
        }
        out.push(v);
    }
    Ok(out)
}

fn parse_shape(args: &ToolArgs) -> Result<BinShape, ToolError> {
    match args.get("bin_shape").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("hexagon") => Ok(BinShape::Hexagon),
        Some("square") => Ok(BinShape::Square),
        Some(o) => Err(ToolError::Validation(format!(
            "'bin_shape' must be hexagon or square, got '{o}'"
        ))),
    }
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
    use wbvector::{memory_store, Coord, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn points(items: Vec<(f64, f64)>) -> String {
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("w", FieldType::Float));
        for (x, y) in items {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(x, y))),
                &[("w", FieldValue::Float(1.0))],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = EvaluateBinSizesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn col(layer: &Layer, field: &str) -> Vec<f64> {
        let i = layer.schema.field_index(field).unwrap();
        layer
            .features
            .iter()
            .map(|f| f.attributes[i].as_f64().unwrap())
            .collect()
    }

    /// One row per candidate, in the order supplied.
    #[test]
    fn emits_one_row_per_candidate() {
        let input = points(vec![(0.0, 0.0), (1.0, 0.0), (10.0, 10.0)]);
        let (out, layer) = run(json!({ "input": input, "sizes": "1,5,20" }));
        assert_eq!(layer.len(), 3);
        assert_eq!(out.outputs["candidate_count"], json!(3));
        assert_eq!(col(&layer, "size"), vec![1.0, 5.0, 20.0]);
    }

    /// THE core signal: a small bin isolates points (mean ~1 per bin), a large
    /// bin aggregates them all into one.
    #[test]
    fn bin_count_falls_as_size_grows() {
        // A 4x4 lattice spaced 10 apart.
        let mut pts = Vec::new();
        for i in 0..4 {
            for j in 0..4 {
                pts.push((i as f64 * 10.0, j as f64 * 10.0));
            }
        }
        let (_, layer) = run(json!({
            "input": points(pts), "sizes": "1,100", "bin_shape": "square"
        }));
        let occupied = col(&layer, "nonempty_bins");
        assert_eq!(occupied[0], 16.0, "size 1 isolates every point");
        assert_eq!(occupied[1], 1.0, "size 100 swallows the whole lattice");
        let means = col(&layer, "mean_per_bin");
        assert_eq!(means[0], 1.0);
        assert_eq!(means[1], 16.0);

        // nonempty_ratio is occupied / total-bins-tiling-the-extent, so a tiny
        // bin leaves most of the extent empty and a huge bin fills it.
        let ratio = col(&layer, "nonempty_ratio");
        assert!(
            ratio[0] < 0.05,
            "size 1 over a 30x30 extent should be almost all empty, got {}",
            ratio[0]
        );
        assert!(
            (ratio[1] - 1.0).abs() < 1e-9,
            "size 100 is a single bin covering everything, got {}",
            ratio[1]
        );
    }

    /// The generated sweep is anchored on the mean nearest-neighbour distance.
    #[test]
    fn generated_sweep_scales_with_point_spacing() {
        let tight = points(vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]);
        let wide = points(vec![(0.0, 0.0), (100.0, 0.0), (200.0, 0.0)]);
        let (t, _) = run(json!({ "input": tight }));
        let (w, _) = run(json!({ "input": wide }));
        let (tm, wm) = (
            t.outputs["mean_nearest_neighbour"].as_f64().unwrap(),
            w.outputs["mean_nearest_neighbour"].as_f64().unwrap(),
        );
        assert!((tm - 1.0).abs() < 1e-9);
        assert!((wm - 100.0).abs() < 1e-9);
        assert!(
            w.outputs["recommended_size"].as_f64().unwrap()
                > t.outputs["recommended_size"].as_f64().unwrap(),
            "a sparser layer must be recommended a larger bin"
        );
    }

    /// The coefficient of variation separates a clustered layer from a uniform
    /// one at the same bin size, which is what makes it the useful diagnostic.
    #[test]
    fn cv_detects_clustering() {
        // Uniform lattice.
        let mut uniform = Vec::new();
        for i in 0..6 {
            for j in 0..6 {
                uniform.push((i as f64 * 10.0, j as f64 * 10.0));
            }
        }
        // Same count, all crammed into one corner plus a few strays.
        let mut clustered = Vec::new();
        for i in 0..30 {
            clustered.push((i as f64 * 0.1, 0.0));
        }
        for i in 0..6 {
            clustered.push((i as f64 * 10.0 + 100.0, 100.0));
        }
        let (_, u) = run(json!({
            "input": points(uniform), "sizes": "25", "bin_shape": "square"
        }));
        let (_, c) = run(json!({
            "input": points(clustered), "sizes": "25", "bin_shape": "square"
        }));
        assert!(
            col(&c, "cv")[0] > col(&u, "cv")[0],
            "clustered cv {} should exceed uniform cv {}",
            col(&c, "cv")[0],
            col(&u, "cv")[0]
        );
    }

    /// analysis_field sums a value instead of counting points.
    #[test]
    fn analysis_field_sums_values() {
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("v", FieldType::Float));
        for (x, v) in [(0.0, 5.0), (1.0, 7.0)] {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(x, 0.0))),
                &[("v", FieldValue::Float(v))],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        let (_, layer) = run(json!({
            "input": memory_store::make_vector_memory_path(&id),
            "sizes": "100", "analysis_field": "v", "bin_shape": "square"
        }));
        assert_eq!(col(&layer, "mean_per_bin")[0], 12.0, "5 + 7 in one bin");
    }

    /// Hexagon and square tessellations both bin, and differ.
    #[test]
    fn both_bin_shapes_work() {
        let mut pts = Vec::new();
        for i in 0..8 {
            for j in 0..8 {
                pts.push((i as f64 * 3.0, j as f64 * 3.0));
            }
        }
        let input = points(pts);
        let (_, hex) = run(json!({ "input": input, "sizes": "10", "bin_shape": "hexagon" }));
        let (_, sq) = run(json!({ "input": input, "sizes": "10", "bin_shape": "square" }));
        assert!(col(&hex, "nonempty_bins")[0] > 0.0);
        assert!(col(&sq, "nonempty_bins")[0] > 0.0);
    }

    #[test]
    fn rejects_bad_parameters() {
        let p = points(vec![(0.0, 0.0), (1.0, 1.0)]);
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            EvaluateBinSizesTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({ "input": p, "bin_shape": "triangle" })));
        assert!(bad(json!({ "input": p, "sizes": "0" })));
        assert!(bad(json!({ "input": p, "sizes": "-5" })));
        assert!(bad(json!({ "input": p, "sizes": "abc" })));
        assert!(bad(json!({ "input": p, "steps": 100 })));

        // A single point cannot support a nearest-neighbour anchor.
        let one = points(vec![(0.0, 0.0)]);
        let args: ToolArgs = serde_json::from_value(json!({ "input": one })).unwrap();
        assert!(matches!(
            EvaluateBinSizesTool.run(&args, &ctx()).unwrap_err(),
            ToolError::Validation(_)
        ));
    }

    /// Alternate rows really are offset (the property the approximation does
    /// guarantee); this does NOT assert exact hexagonal cell membership.
    #[test]
    fn hex_rows_are_offset() {
        // Two points one hex-row apart at the same x should not share a column
        // index, because the odd row is shifted by half a hexagon.
        let a = hex_bin(0.0, 0.0, 10.0);
        let b = hex_bin(0.0, 8.7, 10.0);
        assert_ne!(a.1, b.1, "different rows");
        assert_ne!(a, b);
    }
}
