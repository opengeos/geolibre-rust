//! GeoLibre tool: collapse exactly-coincident incident points into weighted
//! points carrying an `ICOUNT` count field.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Collect Events* (Spatial Statistics).
//! Raw incident data typically stacks many events at the same address/junction;
//! hot-spot tools want one point per location weighted by how many events
//! landed there. `eliminate_coincident_points` *removes* the duplicates and
//! `aggregate_points` merges by a distance threshold into polygons — neither
//! yields the count-bearing point collapse (`ICOUNT`) those downstream tools
//! expect.
//!
//! A single pass hashes each point's coordinate (snapped to `tolerance` map
//! units so genuine duplicates coincide despite float noise) into a map that
//! counts occurrences and remembers the first-seen location. Each unique
//! location becomes one output point whose `ICOUNT` is the number of coincident
//! events there. Deterministic: keys are emitted in a stable sorted order, and
//! `sum(ICOUNT)` always equals the input point count.

use std::collections::{BTreeMap, HashMap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Default snapping tolerance (map units) when the user does not supply one.
/// Small enough to be "exact" for real data yet coarse enough to absorb the
/// last-bit float noise of re-parsed decimal coordinates.
const DEFAULT_TOLERANCE: f64 = 1e-8;

/// Collapses exactly-coincident points into one point per location with an
/// `ICOUNT` field counting the events there.
pub struct CollectEventsTool;

impl Tool for CollectEventsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "collect_events",
            display_name: "Collect Events",
            summary: "Collapse coincident incident points into unique locations, each carrying an ICOUNT count of the events there, like ArcGIS Collect Events.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer (Point or MultiPoint).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output point vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Snapping tolerance in map units; points whose coordinates round to the same cell are treated as coincident (default 1e-8, i.e. effectively exact).",
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
        parse_tolerance(args)?;
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
        let tolerance = parse_tolerance(args)?;

        let layer = load_input_layer(input)?;

        // ── Hash each point's snapped coordinate, counting occurrences and
        //    remembering the first-seen (unsnapped) location as the emitted point.
        let mut order: Vec<(i64, i64)> = Vec::new();
        let mut counts: HashMap<(i64, i64), (usize, f64, f64)> = HashMap::new();
        let mut input_points = 0usize;
        for feature in layer.features.iter() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            match geom {
                Geometry::Point(c) => accumulate(
                    &mut order,
                    &mut counts,
                    &mut input_points,
                    tolerance,
                    c.x,
                    c.y,
                ),
                Geometry::MultiPoint(cs) => {
                    for c in cs {
                        accumulate(
                            &mut order,
                            &mut counts,
                            &mut input_points,
                            tolerance,
                            c.x,
                            c.y,
                        );
                    }
                }
                _ => {}
            }
        }
        if input_points == 0 {
            return Err(ToolError::Execution(
                "no point features in input".to_string(),
            ));
        }

        ctx.progress.info(&format!(
            "collapsing {} point(s) into {} location(s)",
            input_points,
            counts.len()
        ));

        // ── Build the output point layer.
        let mut out = Layer::new("collected_events").with_geom_type(GeometryType::Point);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("ICOUNT", FieldType::Integer));

        // Deterministic order: keys are pushed in first-seen order; sort them so
        // output ordering is independent of input traversal quirks.
        order.sort_unstable();
        let mut max_icount = 0i64;
        let mut sum_icount = 0i64;
        for key in &order {
            let (count, x, y) = counts[key];
            let icount = count as i64;
            max_icount = max_icount.max(icount);
            sum_icount += icount;
            out.add_feature(
                Some(Geometry::Point(Coord::xy(x, y))),
                &[("ICOUNT", icount.into())],
            )
            .map_err(|e| ToolError::Execution(format!("failed writing point: {e}")))?;
        }

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_points".to_string(), json!(input_points));
        outputs.insert("output_points".to_string(), json!(order.len()));
        outputs.insert("sum_icount".to_string(), json!(sum_icount));
        outputs.insert("max_icount".to_string(), json!(max_icount));
        Ok(ToolRunResult { outputs })
    }
}

/// Snaps `(x, y)` to the tolerance grid and folds it into the running tally,
/// keeping the first-seen coordinate as the representative location.
fn accumulate(
    order: &mut Vec<(i64, i64)>,
    counts: &mut HashMap<(i64, i64), (usize, f64, f64)>,
    input_points: &mut usize,
    tolerance: f64,
    x: f64,
    y: f64,
) {
    let key = (
        (x / tolerance).round() as i64,
        (y / tolerance).round() as i64,
    );
    counts
        .entry(key)
        .and_modify(|e| e.0 += 1)
        .or_insert_with(|| {
            order.push(key);
            (1, x, y)
        });
    *input_points += 1;
}

/// Parses the optional `tolerance` parameter (JSON number or numeric string),
/// defaulting to [`DEFAULT_TOLERANCE`] and rejecting non-positive values.
fn parse_tolerance(args: &ToolArgs) -> Result<f64, ToolError> {
    let tol = match args.get("tolerance") {
        None | Some(Value::Null) => DEFAULT_TOLERANCE,
        Some(Value::Number(n)) => n.as_f64().unwrap_or(DEFAULT_TOLERANCE),
        Some(Value::String(s)) if s.trim().is_empty() => DEFAULT_TOLERANCE,
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'tolerance' must be a number".into()))?,
        Some(_) => {
            return Err(ToolError::Validation(
                "'tolerance' must be a number".to_string(),
            ))
        }
    };
    if !tol.is_finite() || tol <= 0.0 {
        return Err(ToolError::Validation(
            "'tolerance' must be a positive number".to_string(),
        ));
    }
    Ok(tol)
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::memory_store;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn layer_of(pts: &[(f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        for (x, y) in pts {
            l.add_feature(Some(Geometry::Point(Coord::xy(*x, *y))), &[])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CollectEventsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn icounts(layer: &Layer) -> Vec<i64> {
        let idx = layer.schema.field_index("ICOUNT").unwrap();
        layer
            .iter()
            .map(|f| f.attributes[idx].as_i64().unwrap())
            .collect()
    }

    /// Three events at one spot + one lone event -> two locations, ICOUNT 3 and 1.
    #[test]
    fn collapses_coincident_points() {
        let pts = [(0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (10.0, 10.0)];
        let input = layer_of(&pts);
        let (out, layer) = run(json!({ "input": input }));
        assert_eq!(out.outputs["input_points"], json!(4));
        assert_eq!(out.outputs["output_points"], json!(2));
        let mut ic = icounts(&layer);
        ic.sort_unstable();
        assert_eq!(ic, vec![1, 3]);
    }

    /// The core invariant: sum(ICOUNT) == number of input points.
    #[test]
    fn sum_icount_equals_input_count() {
        let pts = [
            (1.0, 1.0),
            (1.0, 1.0),
            (2.0, 2.0),
            (2.0, 2.0),
            (2.0, 2.0),
            (5.0, 9.0),
        ];
        let input = layer_of(&pts);
        let (out, layer) = run(json!({ "input": input }));
        assert_eq!(out.outputs["sum_icount"], json!(6));
        assert_eq!(out.outputs["max_icount"], json!(3));
        let total: i64 = icounts(&layer).iter().sum();
        assert_eq!(total, 6);
    }

    /// All-distinct points pass through unchanged, every ICOUNT == 1.
    #[test]
    fn distinct_points_pass_through() {
        let pts = [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)];
        let input = layer_of(&pts);
        let (out, layer) = run(json!({ "input": input }));
        assert_eq!(out.outputs["output_points"], json!(3));
        assert!(icounts(&layer).iter().all(|&c| c == 1));
    }

    /// A coarse tolerance snaps near-but-not-equal points together.
    #[test]
    fn tolerance_snaps_near_points() {
        let pts = [(0.0, 0.0), (0.4, 0.3), (100.0, 100.0)];
        let input = layer_of(&pts);
        let (out, _l) = run(json!({ "input": input, "tolerance": 1.0 }));
        assert_eq!(out.outputs["output_points"], json!(2));
        // Numeric string accepted too.
        let (out2, _l2) = run(json!({ "input": input, "tolerance": "1.0" }));
        assert_eq!(out2.outputs["output_points"], json!(2));
    }

    #[test]
    fn rejects_bad_parameters() {
        let input = layer_of(&[(0.0, 0.0)]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "tolerance": -1.0 })).unwrap();
        assert!(CollectEventsTool.validate(&args).is_err());
        let args: ToolArgs = serde_json::from_value(json!({ "tolerance": 1.0 })).unwrap();
        assert!(CollectEventsTool.validate(&args).is_err());
    }
}
