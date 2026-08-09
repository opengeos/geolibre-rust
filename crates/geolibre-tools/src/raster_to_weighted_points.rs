//! GeoLibre tool: points distributed in proportion to a raster's values.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Raster To Weighted Points* (Spatial
//! Analyst).
//!
//! ## Why the catalog needs it
//!
//! Turning a density surface back into discrete points is how a continuous
//! field becomes something you can route to, cluster, or sample: population
//! density into candidate service locations, a habitat suitability surface into
//! survey sites, a cost surface into agent start positions. The requirement is
//! that point *density* track cell value, so that a cell twice as populous gets
//! twice as many points.
//!
//! Nothing in either registry does that. `random_sample` and
//! `random_points_in_polygon` are uniform — they ignore the surface entirely.
//! `create_spatially_balanced_points` (round 2) spreads points *evenly*, which
//! is the opposite requirement. `raster_to_vector_points` emits one point per
//! cell, so it neither thins nor weights.
//!
//! ## Determinism
//!
//! There is no RNG here. The per-cell allocation uses the **largest-remainder**
//! method — exact integer quotas that sum to the requested total — and the
//! points inside a cell are placed on a deterministic pattern. That matters
//! because the WASM builds have no random source the native build shares, so a
//! sampled result has to be reproducible from the inputs alone.
//!
//! Placement patterns, following ArcGIS's list:
//!
//! * `fibonacci_lattice` (default) — a golden-ratio lattice, the lowest-
//!   discrepancy way to fill a square, so points neither clump nor grid-align;
//! * `fibonacci_spiral` — the same golden angle in polar form, which reads
//!   naturally for radial phenomena;
//! * `circular` — concentric rings;
//! * `centroid` — every point at the cell centre, for a purely aggregate use.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::args_common::{band_index, bool_or, choice_or, opt_f64, req_str, usize_or};
use crate::common::{load_input_raster, parse_optional_output};
use crate::vector_common::write_or_store_layer;

/// The golden ratio conjugate, the irrational step that makes a Fibonacci
/// lattice low-discrepancy.
const GOLDEN: f64 = 0.618_033_988_749_894_9;

pub struct RasterToWeightedPointsTool;

impl Tool for RasterToWeightedPointsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "raster_to_weighted_points",
            display_name: "Raster To Weighted Points",
            summary: "Distributes a fixed number of points across a raster so their density is proportional to cell value, placing them on a deterministic Fibonacci lattice, spiral, circular or centroid pattern (ArcGIS Raster To Weighted Points). Nothing in either registry weights by a surface: random_sample and random_points_in_polygon are uniform, create_spatially_balanced_points deliberately spreads points evenly, and raster_to_vector_points emits one per cell without thinning or weighting. Allocation uses the largest-remainder method and needs no RNG, so results are reproducible in WASM.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input surface raster; cell values are the weights.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output point layer. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_number_of_points",
                    description: "Total points to distribute (default 1000).",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "In-cell placement: 'fibonacci_lattice' (default), 'fibonacci_spiral', 'circular', or 'centroid'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_value",
                    description: "Cells at or below this weight get no points (default 0). Negative weights are always excluded.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_points_per_cell",
                    description: "Cap the points any single cell may receive, so one hot cell cannot absorb the whole budget.",
                    required: false,
                },
                ToolParamSpec {
                    name: "include_weight",
                    description: "Write the source cell value onto each point (default true).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band supplying the weights (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input_path = req_str(args, "input")?.to_string();
        let prm = parse_params(args)?;
        let band = band_index(args, "band")?;
        let output = parse_optional_output(args, "output")?;

        let raster = load_input_raster(&input_path)?;
        let (rows, cols) = (raster.rows, raster.cols);

        // Eligible cells and their weights.
        let mut cells: Vec<(usize, f64)> = Vec::new();
        let mut total_weight = 0.0;
        for r in 0..rows {
            for c in 0..cols {
                let v = raster.get(band, r as isize, c as isize);
                if v == raster.nodata || !v.is_finite() || v <= prm.min_value || v <= 0.0 {
                    continue;
                }
                cells.push((r * cols + c, v));
                total_weight += v;
            }
        }
        if cells.is_empty() || total_weight <= 0.0 {
            return Err(ToolError::Execution(
                "no cell carries a positive weight, so no points can be distributed".to_string(),
            ));
        }

        let counts = allocate(&cells, total_weight, prm.total, prm.max_per_cell);
        let placed: usize = counts.iter().sum();
        ctx.progress.info(&format!(
            "{} eligible cell(s), total weight {total_weight:.4}, placing {placed} point(s) by {}",
            cells.len(),
            prm.method.label()
        ));

        let mut layer = Layer::new("weighted_points");
        layer.geom_type = Some(GeometryType::Point);
        if let Some(e) = raster.crs.epsg {
            layer = layer.with_crs_epsg(e);
        }
        layer.add_field(FieldDef::new("id", FieldType::Integer));
        layer.add_field(FieldDef::new("row", FieldType::Integer));
        layer.add_field(FieldDef::new("col", FieldType::Integer));
        if prm.include_weight {
            layer.add_field(FieldDef::new("weight", FieldType::Float));
        }

        let (csx, csy) = (raster.cell_size_x, raster.cell_size_y);
        let y_max = raster.y_min + rows as f64 * csy;
        let mut fid = 0u64;

        for (ci, &(index, weight)) in cells.iter().enumerate() {
            let n = counts[ci];
            if n == 0 {
                continue;
            }
            let (r, c) = (index / cols, index % cols);
            let x0 = raster.x_min + c as f64 * csx;
            let y0 = y_max - (r as f64 + 1.0) * csy;

            for k in 0..n {
                let (u, v) = prm.method.offset(k, n);
                let x = x0 + u * csx;
                let y = y0 + v * csy;
                let mut f = Feature::with_geometry(
                    fid,
                    Geometry::Point(Coord::xy(x, y)),
                    layer.schema.len(),
                );
                f.set_by_index(0, FieldValue::Integer(fid as i64));
                f.set_by_index(1, FieldValue::Integer(r as i64));
                f.set_by_index(2, FieldValue::Integer(c as i64));
                if prm.include_weight {
                    f.set_by_index(3, FieldValue::Float(weight));
                }
                layer.push(f);
                fid += 1;
            }
            ctx.progress
                .progress((ci as f64 + 1.0) / cells.len() as f64);
        }

        let point_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("point_count".to_string(), json!(point_count));
        outputs.insert("requested_points".to_string(), json!(prm.total));
        outputs.insert("eligible_cells".to_string(), json!(cells.len()));
        outputs.insert("total_weight".to_string(), json!(total_weight));
        outputs.insert("method".to_string(), json!(prm.method.label()));
        Ok(ToolRunResult { outputs })
    }
}

/// Distributes `total` points across cells in proportion to their weight, using
/// the largest-remainder method.
///
/// Proportional shares are rarely whole numbers. Rounding each independently
/// would miss or overshoot the requested total, sometimes badly; the
/// largest-remainder method floors every share and then hands the leftover
/// points to the cells with the biggest fractional parts, so the result is both
/// exact and deterministic. Ties break on cell order, which is fixed.
fn allocate(
    cells: &[(usize, f64)],
    total_weight: f64,
    total: usize,
    max_per_cell: Option<usize>,
) -> Vec<usize> {
    let mut counts = vec![0usize; cells.len()];
    let mut remainders: Vec<(f64, usize)> = Vec::with_capacity(cells.len());
    let mut assigned = 0usize;

    for (i, &(_, w)) in cells.iter().enumerate() {
        let share = total as f64 * w / total_weight;
        let mut whole = share.floor() as usize;
        if let Some(cap) = max_per_cell {
            whole = whole.min(cap);
        }
        counts[i] = whole;
        assigned += whole;
        remainders.push((share - share.floor(), i));
    }

    // Hand out the leftovers, largest fractional part first.
    remainders.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
    let mut leftover = total.saturating_sub(assigned);
    for &(_, i) in &remainders {
        if leftover == 0 {
            break;
        }
        if let Some(cap) = max_per_cell {
            if counts[i] >= cap {
                continue;
            }
        }
        counts[i] += 1;
        leftover -= 1;
    }
    counts
}

/// Where within a cell the `k`-th of `n` points sits, in unit coordinates.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    FibonacciLattice,
    FibonacciSpiral,
    Circular,
    Centroid,
}

impl Method {
    fn label(self) -> &'static str {
        match self {
            Method::FibonacciLattice => "fibonacci_lattice",
            Method::FibonacciSpiral => "fibonacci_spiral",
            Method::Circular => "circular",
            Method::Centroid => "centroid",
        }
    }

    /// Unit-square offset of point `k` of `n`, kept strictly inside the cell so
    /// no point lands on a shared boundary.
    fn offset(self, k: usize, n: usize) -> (f64, f64) {
        if n == 1 || self == Method::Centroid {
            return (0.5, 0.5);
        }
        match self {
            Method::Centroid => (0.5, 0.5),
            Method::FibonacciLattice => {
                // Stratified in one axis, golden-ratio stepped in the other:
                // the standard low-discrepancy fill of a square.
                let y = (k as f64 + 0.5) / n as f64;
                let x = ((k as f64 * GOLDEN) % 1.0).clamp(0.0, 1.0);
                (clamp_inside(x), clamp_inside(y))
            }
            Method::FibonacciSpiral => {
                // Golden-angle spiral with sqrt radius, so area density is even.
                let theta = 2.0 * PI * GOLDEN * k as f64;
                let radius = 0.5 * ((k as f64 + 0.5) / n as f64).sqrt();
                (
                    clamp_inside(0.5 + radius * theta.cos()),
                    clamp_inside(0.5 + radius * theta.sin()),
                )
            }
            Method::Circular => {
                // Concentric rings: ring r holds roughly 6r points, which keeps
                // the spacing between rings and along them comparable.
                let mut remaining = k;
                let mut ring = 0usize;
                loop {
                    let capacity = if ring == 0 { 1 } else { 6 * ring };
                    if remaining < capacity {
                        break;
                    }
                    remaining -= capacity;
                    ring += 1;
                }
                if ring == 0 {
                    return (0.5, 0.5);
                }
                let capacity = 6 * ring;
                let theta = 2.0 * PI * remaining as f64 / capacity as f64;
                // Rings are spaced to fit inside the cell however many there are.
                let max_ring = ((n as f64 / 6.0).sqrt().ceil() as usize).max(1);
                let radius = 0.5 * ring as f64 / (max_ring as f64 + 0.5);
                (
                    clamp_inside(0.5 + radius * theta.cos()),
                    clamp_inside(0.5 + radius * theta.sin()),
                )
            }
        }
    }
}

/// Keeps an offset strictly inside the unit cell.
fn clamp_inside(v: f64) -> f64 {
    v.clamp(1e-6, 1.0 - 1e-6)
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    total: usize,
    method: Method,
    min_value: f64,
    max_per_cell: Option<usize>,
    include_weight: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let total = usize_or(args, "max_number_of_points", 1000)?;
    if total == 0 {
        return Err(ToolError::Validation(
            "'max_number_of_points' must be at least 1".to_string(),
        ));
    }
    let method = match choice_or(
        args,
        "method",
        &[
            "fibonacci_lattice",
            "fibonacci_spiral",
            "circular",
            "centroid",
        ],
        "fibonacci_lattice",
    )? {
        "fibonacci_spiral" => Method::FibonacciSpiral,
        "circular" => Method::Circular,
        "centroid" => Method::Centroid,
        _ => Method::FibonacciLattice,
    };
    let min_value = opt_f64(args, "min_value")?.unwrap_or(0.0);
    if !min_value.is_finite() {
        return Err(ToolError::Validation(
            "'min_value' must be finite".to_string(),
        ));
    }
    let max_per_cell = match crate::args_common::opt_usize(args, "max_points_per_cell")? {
        None => None,
        Some(v) if v >= 1 => Some(v),
        Some(_) => {
            return Err(ToolError::Validation(
                "'max_points_per_cell' must be at least 1".to_string(),
            ))
        }
    };
    let include_weight = bool_or(args, "include_weight", true)?;
    Ok(Params {
        total,
        method,
        min_value,
        max_per_cell,
        include_weight,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vector_common::load_input_layer;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
    use serde_json::Value;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn raster_of(cols: usize, rows: usize, vals: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 10.0,
            cell_size_y: Some(10.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(32610),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, vals[row * cols + col])
                    .unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, BTreeMap<String, Value>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = RasterToWeightedPointsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (layer, out.outputs)
    }

    /// Points per cell counted from the output.
    fn per_cell(layer: &Layer) -> BTreeMap<(i64, i64), usize> {
        let ri = layer.schema.field_index("row").unwrap();
        let ci = layer.schema.field_index("col").unwrap();
        let mut out: BTreeMap<(i64, i64), usize> = BTreeMap::new();
        for f in layer.iter() {
            let (FieldValue::Integer(r), FieldValue::Integer(c)) =
                (&f.attributes[ri], &f.attributes[ci])
            else {
                panic!("row/col must be integers")
            };
            *out.entry((*r, *c)).or_default() += 1;
        }
        out
    }

    /// The defining property: a cell worth three times another gets three times
    /// the points. This is exactly what `random_sample` cannot do.
    #[test]
    fn density_is_proportional_to_weight() {
        // Two cells, weights 1 and 3, 400 points -> 100 and 300.
        let (layer, outputs) = run(json!({
            "input": raster_of(2, 1, &[1.0, 3.0]), "max_number_of_points": 400
        }));
        assert_eq!(outputs["point_count"].as_u64().unwrap(), 400);
        let counts = per_cell(&layer);
        assert_eq!(counts[&(0, 0)], 100);
        assert_eq!(counts[&(0, 1)], 300);
    }

    /// The largest-remainder allocation hits the requested total exactly, even
    /// when the shares are all fractional — independent rounding would not.
    #[test]
    fn total_is_exact_with_fractional_shares() {
        // Three equal cells and 10 points: shares are 3.33 each.
        let (layer, outputs) = run(json!({
            "input": raster_of(3, 1, &[1.0, 1.0, 1.0]), "max_number_of_points": 10
        }));
        assert_eq!(layer.len(), 10, "must place exactly the requested total");
        assert_eq!(outputs["point_count"].as_u64().unwrap(), 10);
        let counts = per_cell(&layer);
        let mut sizes: Vec<usize> = counts.values().copied().collect();
        sizes.sort_unstable();
        assert_eq!(sizes, vec![3, 3, 4], "the leftover goes to one cell");
    }

    /// Results are reproducible: the same input gives the same points, with no
    /// RNG anywhere.
    #[test]
    fn output_is_deterministic() {
        let v = vec![1.0, 5.0, 2.0, 9.0];
        let first = run(json!({
            "input": raster_of(2, 2, &v), "max_number_of_points": 37
        }))
        .0;
        let second = run(json!({
            "input": raster_of(2, 2, &v), "max_number_of_points": 37
        }))
        .0;
        assert_eq!(first.len(), second.len());
        let coords = |l: &Layer| -> Vec<(u64, u64)> {
            l.iter()
                .map(|f| match f.geometry.as_ref() {
                    Some(Geometry::Point(p)) => (p.x.to_bits(), p.y.to_bits()),
                    _ => panic!("expected points"),
                })
                .collect()
        };
        assert_eq!(coords(&first), coords(&second), "runs must be identical");
    }

    /// Points land inside the cell they were allocated to.
    #[test]
    fn points_fall_inside_their_own_cell() {
        let (layer, _) = run(json!({
            "input": raster_of(3, 3, &[1.0; 9]), "max_number_of_points": 60
        }));
        let ri = layer.schema.field_index("row").unwrap();
        let ci = layer.schema.field_index("col").unwrap();
        let y_max = 3.0 * 10.0;
        for f in layer.iter() {
            let (FieldValue::Integer(r), FieldValue::Integer(c)) =
                (&f.attributes[ri], &f.attributes[ci])
            else {
                panic!()
            };
            let Some(Geometry::Point(p)) = f.geometry.as_ref() else {
                panic!()
            };
            let x0 = *c as f64 * 10.0;
            let y0 = y_max - (*r as f64 + 1.0) * 10.0;
            assert!(
                p.x > x0 && p.x < x0 + 10.0,
                "x {} outside cell column {c}",
                p.x
            );
            assert!(
                p.y > y0 && p.y < y0 + 10.0,
                "y {} outside cell row {r}",
                p.y
            );
        }
    }

    /// Every placement pattern is usable and spreads points within a cell —
    /// except `centroid`, which deliberately stacks them.
    #[test]
    fn placement_methods_behave_as_described() {
        for method in ["fibonacci_lattice", "fibonacci_spiral", "circular"] {
            let (layer, _) = run(json!({
                "input": raster_of(1, 1, &[1.0]), "max_number_of_points": 20, "method": method
            }));
            assert_eq!(layer.len(), 20, "{method} should place all 20");
            let pts: Vec<(f64, f64)> = layer
                .iter()
                .map(|f| match f.geometry.as_ref() {
                    Some(Geometry::Point(p)) => (p.x, p.y),
                    _ => panic!(),
                })
                .collect();
            let distinct: std::collections::HashSet<(u64, u64)> =
                pts.iter().map(|p| (p.0.to_bits(), p.1.to_bits())).collect();
            assert_eq!(
                distinct.len(),
                20,
                "{method} should spread points, not stack them"
            );
        }

        let (centroid, _) = run(json!({
            "input": raster_of(1, 1, &[1.0]), "max_number_of_points": 20, "method": "centroid"
        }));
        let distinct: std::collections::HashSet<(u64, u64)> = centroid
            .iter()
            .map(|f| match f.geometry.as_ref() {
                Some(Geometry::Point(p)) => (p.x.to_bits(), p.y.to_bits()),
                _ => panic!(),
            })
            .collect();
        assert_eq!(distinct.len(), 1, "centroid stacks every point at the centre");
    }

    /// Zero, negative and no-data cells get nothing.
    #[test]
    fn nonpositive_and_nodata_cells_are_skipped() {
        let (layer, outputs) = run(json!({
            "input": raster_of(4, 1, &[0.0, -5.0, -9999.0, 8.0]),
            "max_number_of_points": 50
        }));
        assert_eq!(outputs["eligible_cells"].as_u64().unwrap(), 1);
        let counts = per_cell(&layer);
        assert_eq!(counts.len(), 1);
        assert_eq!(counts[&(0, 3)], 50);
    }

    /// `min_value` raises the eligibility bar.
    #[test]
    fn min_value_excludes_low_cells() {
        let (_, outputs) = run(json!({
            "input": raster_of(4, 1, &[1.0, 2.0, 3.0, 4.0]),
            "max_number_of_points": 40, "min_value": 2.0
        }));
        // Only the 3 and 4 cells clear a min_value of 2.
        assert_eq!(outputs["eligible_cells"].as_u64().unwrap(), 2);
    }

    /// The per-cell cap stops one hot cell absorbing the whole budget.
    #[test]
    fn max_points_per_cell_is_respected() {
        let (layer, _) = run(json!({
            "input": raster_of(2, 1, &[1.0, 99.0]),
            "max_number_of_points": 100, "max_points_per_cell": 10
        }));
        let counts = per_cell(&layer);
        for (cell, n) in &counts {
            assert!(*n <= 10, "cell {cell:?} got {n} points, over the cap");
        }
    }

    /// The source weight travels with each point, and can be suppressed.
    #[test]
    fn weight_field_is_optional() {
        let (with, _) = run(json!({
            "input": raster_of(1, 1, &[7.5]), "max_number_of_points": 3
        }));
        let wi = with.schema.field_index("weight").unwrap();
        assert!(matches!(
            with.iter().next().unwrap().attributes[wi],
            FieldValue::Float(v) if (v - 7.5).abs() < 1e-6
        ));

        let (without, _) = run(json!({
            "input": raster_of(1, 1, &[7.5]), "max_number_of_points": 3,
            "include_weight": false
        }));
        assert!(without.schema.field_index("weight").is_none());
    }

    /// An all-zero surface has nothing to weight by; say so.
    #[test]
    fn empty_surface_is_an_error() {
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": raster_of(3, 3, &[0.0; 9]) })).unwrap();
        let err = RasterToWeightedPointsTool.run(&args, &ctx()).unwrap_err();
        assert!(
            format!("{err:?}").contains("positive weight"),
            "expected a no-weight error, got {err:?}"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            RasterToWeightedPointsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({"input": "a.tif", "max_number_of_points": 0})).is_err());
        assert!(bad(json!({"input": "a.tif", "method": "poisson"})).is_err());
        assert!(bad(json!({"input": "a.tif", "max_points_per_cell": 0})).is_err());
        assert!(bad(json!({"input": "a.tif", "method": "circular"})).is_ok());
    }
}
