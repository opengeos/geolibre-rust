//! GeoLibre tool: thin a vector point layer to a minimum spacing.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Reduce Point Density* (Bathymetry).
//!
//! ## Why the catalog could only add points
//!
//! `lidar_thin` and `lidar_thin_high_density` are **LAS-only** — they cannot
//! touch a shapefile or GeoJSON point layer. `densify_sampling_network` adds
//! points; `create_spatially_balanced_points` and
//! `create_spatial_sampling_locations` generate new ones. Nothing decimated an
//! existing layer, and `find_identical` only catches exact coincidences, not a
//! dense cluster of distinct-but-nearby points.
//!
//! That gap has teeth because dense point layers break the interpolators.
//! Bathymetric soundings, GPS traces and sensor logs arrive with tens of
//! thousands of near-coincident points; feeding those to `idw_interpolation`,
//! the kriging family, `thin_plate_spline` or `natural_neighbour_interpolation`
//! is slow at best and singular at worst, because near-duplicate points make the
//! kriging system ill-conditioned. Thinning is the standard prep step.
//!
//! ## Safety: `keep_field`
//!
//! In bathymetry, dropping the shallowest sounding in a cluster is a
//! navigational hazard, not a cosmetic loss. Points flagged by `keep_field` are
//! inserted first and are never displaced, so a shoal survives thinning even
//! when it sits inside a dense cluster.

use std::collections::BTreeMap;

use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldValue, Geometry, Layer};

use crate::args_common::{choice_or, opt_positive_f64};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

const METHODS: [&str; 2] = ["spacing", "bin"];
const ORDERS: [&str; 2] = ["descending", "ascending"];

pub struct ReducePointDensityTool;

impl Tool for ReducePointDensityTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "reduce_point_density",
            display_name: "Reduce Point Density",
            summary: "Thins a vector point layer so no two retained points are closer than a minimum spacing, or keeps one point per grid cell (ArcGIS Reduce Point Density). lidar_thin is LAS-only and find_identical catches only exact duplicates, so dense sounding, GPS and sensor layers could not be decimated before interpolation.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point (or multipoint) layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'spacing' (default): greedy minimum-separation thinning. 'bin': keep one point per grid cell.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_distance",
                    description: "Minimum separation between retained points, in layer units. Required by method 'spacing'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "bin_size",
                    description: "Grid cell size in layer units. Required by method 'bin'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "sort_field",
                    description: "Numeric attribute deciding which point wins a conflict. Without it, input order decides.",
                    required: false,
                },
                ToolParamSpec {
                    name: "sort_order",
                    description: "'descending' (default, highest sort_field wins) or 'ascending'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "keep_field",
                    description: "Attribute marking points that must always be retained (non-zero / true / non-empty). Use for hazards such as shoal soundings.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output layer of retained points. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_removed",
                    description: "Optional layer receiving the dropped points, for review.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if parse_optional_str(args, "input")?.is_none() {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        let method = choice_or(args, "method", &METHODS, "spacing")?;
        choice_or(args, "sort_order", &ORDERS, "descending")?;
        // opt_positive_f64 rejects zero and negatives, which would otherwise
        // silently retain everything (spacing) or divide by zero (bin).
        let min_distance = opt_positive_f64(args, "min_distance")?;
        let bin_size = opt_positive_f64(args, "bin_size")?;
        match method {
            "bin" if bin_size.is_none() => Err(ToolError::Validation(
                "method 'bin' requires a positive 'bin_size'".to_string(),
            )),
            "spacing" if min_distance.is_none() => Err(ToolError::Validation(
                "method 'spacing' requires a positive 'min_distance'".to_string(),
            )),
            _ => Ok(()),
        }
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = parse_optional_str(args, "input")?
            .ok_or_else(|| ToolError::Validation("missing required parameter 'input'".into()))?;
        let method = choice_or(args, "method", &METHODS, "spacing")?;
        let descending = choice_or(args, "sort_order", &ORDERS, "descending")? == "descending";
        let min_distance = opt_positive_f64(args, "min_distance")?;
        let bin_size = opt_positive_f64(args, "bin_size")?;
        let sort_field = parse_optional_str(args, "sort_field")?.map(str::to_string);
        let keep_field = parse_optional_str(args, "keep_field")?.map(str::to_string);
        let output = parse_optional_str(args, "output")?;
        let output_removed = parse_optional_str(args, "output_removed")?;

        let layer = load_input_layer(input)?;

        let sort_idx = match &sort_field {
            Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                ToolError::Validation(format!("sort_field '{f}' not found in the input layer"))
            })?),
            None => None,
        };
        let keep_idx = match &keep_field {
            Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                ToolError::Validation(format!("keep_field '{f}' not found in the input layer"))
            })?),
            None => None,
        };

        // Collect candidate points. A non-point geometry cannot be thinned by
        // spacing in any meaningful way; passing it through would make the
        // output type ambiguous, so it is an error rather than a silent copy.
        let mut pts: Vec<Pt> = Vec::with_capacity(layer.features.len());
        for (fid, f) in layer.iter().enumerate() {
            let (x, y) = match f.geometry.as_ref() {
                Some(Geometry::Point(p)) => (p.x, p.y),
                Some(Geometry::MultiPoint(ps)) if ps.len() == 1 => (ps[0].x, ps[0].y),
                Some(_) => {
                    return Err(ToolError::Validation(
                        "reduce_point_density expects a point layer; feature geometry is not a \
                         single point"
                            .to_string(),
                    ))
                }
                None => continue,
            };
            if !x.is_finite() || !y.is_finite() {
                continue;
            }
            pts.push(Pt {
                fid,
                x,
                y,
                sort: sort_idx
                    .and_then(|i| f.attributes.get(i))
                    .and_then(as_f64)
                    .unwrap_or(0.0),
                pinned: keep_idx
                    .and_then(|i| f.attributes.get(i))
                    .is_some_and(truthy),
            });
        }
        if pts.is_empty() {
            return Err(ToolError::Execution(
                "the input layer has no usable point geometry".to_string(),
            ));
        }
        ctx.progress
            .info(&format!("{} point(s), method {method}", pts.len()));

        // Pinned points first, then by sort priority, then by input order so
        // the result is deterministic for ties.
        pts.sort_by(|a, b| {
            b.pinned
                .cmp(&a.pinned)
                .then_with(|| {
                    if sort_idx.is_none() {
                        std::cmp::Ordering::Equal
                    } else if descending {
                        b.sort.total_cmp(&a.sort)
                    } else {
                        a.sort.total_cmp(&b.sort)
                    }
                })
                .then_with(|| a.fid.cmp(&b.fid))
        });

        let kept: Vec<usize> = if method == "bin" {
            let size = bin_size.ok_or_else(|| {
                ToolError::Validation("method 'bin' requires 'bin_size'".to_string())
            })?;
            thin_by_bin(&pts, size)
        } else {
            let d = min_distance.ok_or_else(|| {
                ToolError::Validation("method 'spacing' requires 'min_distance'".to_string())
            })?;
            thin_by_spacing(&pts, d)?
        };

        let mut keep_flag = vec![false; layer.features.len()];
        for &fid in &kept {
            keep_flag[fid] = true;
        }

        let mut out = clone_schema(&layer);
        let mut removed = clone_schema(&layer);
        for (fid, f) in layer.iter().enumerate() {
            if keep_flag[fid] {
                out.push(f.clone());
            } else {
                removed.push(f.clone());
            }
        }

        let kept_n = out.features.len();
        let removed_n = removed.features.len();
        let out_path = write_or_store_layer(out, output)?;
        // Emitted unconditionally: a secondary output gated on a supplied path
        // silently vanishes for callers with no scratch directory.
        let removed_path = write_or_store_layer(removed, output_removed)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("output_removed".to_string(), json!(removed_path));
        outputs.insert("input_count".to_string(), json!(layer.features.len()));
        outputs.insert("kept_count".to_string(), json!(kept_n));
        outputs.insert("removed_count".to_string(), json!(removed_n));
        outputs.insert("method".to_string(), json!(method));
        Ok(ToolRunResult { outputs })
    }
}

struct Pt {
    fid: usize,
    x: f64,
    y: f64,
    sort: f64,
    pinned: bool,
}

/// Greedy minimum-separation thinning: walk the points in priority order and
/// accept one only when nothing already accepted lies within `min_distance`.
///
/// The accepted set is held in a k-d tree so each test is a radius query rather
/// than a scan of everything kept so far.
fn thin_by_spacing(pts: &[Pt], min_distance: f64) -> Result<Vec<usize>, ToolError> {
    let mut tree: KdTree<f64, usize, [f64; 2]> = KdTree::new(2);
    let mut kept = Vec::new();
    // The kdtree crate's `within` works in the metric its distance function
    // returns, so the radius has to be squared to match squared_euclidean.
    let r2 = min_distance * min_distance;
    // `within` is inclusive, but a point sitting *exactly* min_distance away is
    // not "closer than min_distance" and must be kept — otherwise a regular
    // grid at precisely the requested spacing would be decimated by half. The
    // returned squared distances are re-checked strictly, with a relative slack
    // so a boundary pair is not split by one ulp of rounding.
    let cutoff = r2 * (1.0 - 1e-9);
    for p in pts {
        let neighbours = tree
            .within(&[p.x, p.y], r2, &squared_euclidean)
            .map_err(|e| ToolError::Execution(format!("kd-tree query failed: {e:?}")))?;
        let too_close = neighbours.iter().any(|(d2, _)| *d2 < cutoff);
        if !too_close {
            tree.add([p.x, p.y], p.fid)
                .map_err(|e| ToolError::Execution(format!("kd-tree insert failed: {e:?}")))?;
            kept.push(p.fid);
        }
    }
    Ok(kept)
}

/// Grid thinning: keep the highest-priority point in each `bin_size` cell.
///
/// Points arrive already sorted by priority, so the first point to claim a cell
/// is the winner and later ones in that cell are dropped.
fn thin_by_bin(pts: &[Pt], bin_size: f64) -> Vec<usize> {
    let mut seen: BTreeMap<(i64, i64), usize> = BTreeMap::new();
    for p in pts {
        let key = (
            (p.x / bin_size).floor() as i64,
            (p.y / bin_size).floor() as i64,
        );
        seen.entry(key).or_insert(p.fid);
    }
    let mut kept: Vec<usize> = seen.into_values().collect();
    kept.sort_unstable();
    kept
}

fn clone_schema(layer: &Layer) -> Layer {
    let mut out = Layer::new(&layer.name);
    if let Some(gt) = layer.geom_type {
        out = out.with_geom_type(gt);
    }
    if let Some(e) = layer.crs_epsg() {
        out = out.with_crs_epsg(e);
    }
    for f in layer.schema.fields() {
        out.add_field(f.clone());
    }
    out
}

fn as_f64(v: &FieldValue) -> Option<f64> {
    match v {
        FieldValue::Integer(i) => Some(*i as f64),
        FieldValue::Float(f) if f.is_finite() => Some(*f),
        FieldValue::Boolean(b) => Some(if *b { 1.0 } else { 0.0 }),
        FieldValue::Text(s) => s.trim().parse::<f64>().ok(),
        _ => None,
    }
}

fn truthy(v: &FieldValue) -> bool {
    match v {
        FieldValue::Boolean(b) => *b,
        FieldValue::Integer(i) => *i != 0,
        FieldValue::Float(f) => *f != 0.0,
        FieldValue::Text(s) => {
            let t = s.trim();
            !t.is_empty() && !t.eq_ignore_ascii_case("false") && t != "0"
        }
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{FieldDef, FieldType, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Points on a line at x = 0, 1, 2, ... with a `depth` and `shoal` field.
    fn line_points(n: usize) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("depth", FieldType::Float));
        l.add_field(FieldDef::new("shoal", FieldType::Integer));
        for i in 0..n {
            l.add_feature(
                Some(Geometry::point(i as f64, 0.0)),
                &[("depth", (i as f64).into()), ("shoal", 0i64.into())],
            )
            .unwrap();
        }
        store(l)
    }

    fn store(l: Layer) -> String {
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = ReducePointDensityTool.run(&args, &ctx()).unwrap();
        let kept = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        let removed = load_input_layer(res.outputs["output_removed"].as_str().unwrap()).unwrap();
        (kept, removed, res)
    }

    fn xs(l: &Layer) -> Vec<f64> {
        l.iter()
            .map(|f| match f.geometry.as_ref() {
                Some(Geometry::Point(p)) => p.x,
                _ => panic!("expected points"),
            })
            .collect()
    }

    #[test]
    fn spacing_keeps_only_points_at_least_min_distance_apart() {
        // Points at 0..5 with min_distance 2 leaves 0, 2, 4.
        let (kept, _, res) = run(json!({
            "input": line_points(6), "min_distance": 2.0,
        }));
        assert_eq!(xs(&kept), vec![0.0, 2.0, 4.0]);
        assert_eq!(res.outputs["kept_count"], json!(3));
        assert_eq!(res.outputs["removed_count"], json!(3));
    }

    #[test]
    fn every_retained_pair_really_is_far_enough_apart() {
        // The invariant the tool exists to guarantee.
        let mut l = Layer::new("cluster").with_geom_type(GeometryType::Point);
        for i in 0..40 {
            let a = i as f64 * 0.37;
            l.add_feature(Some(Geometry::point(a.cos() * 5.0, a.sin() * 5.0)), &[])
                .unwrap();
        }
        let (kept, _, _) = run(json!({"input": store(l), "min_distance": 2.0}));
        let pts: Vec<(f64, f64)> = kept
            .iter()
            .map(|f| match f.geometry.as_ref() {
                Some(Geometry::Point(p)) => (p.x, p.y),
                _ => panic!(),
            })
            .collect();
        for i in 0..pts.len() {
            for j in (i + 1)..pts.len() {
                let d = ((pts[i].0 - pts[j].0).powi(2) + (pts[i].1 - pts[j].1).powi(2)).sqrt();
                assert!(d >= 2.0 - 1e-9, "kept two points only {d} apart");
            }
        }
        assert!(!pts.is_empty());
    }

    #[test]
    fn points_exactly_min_distance_apart_are_both_kept() {
        // A regular grid at precisely the requested spacing must survive
        // intact; an inclusive radius test would decimate it by half.
        let (kept, _, _) = run(json!({
            "input": line_points(4), "min_distance": 1.0,
        }));
        assert_eq!(xs(&kept), vec![0.0, 1.0, 2.0, 3.0]);
    }

    #[test]
    fn nothing_is_dropped_when_the_layer_is_already_sparse() {
        let (kept, removed, _) = run(json!({
            "input": line_points(4), "min_distance": 0.5,
        }));
        assert_eq!(kept.features.len(), 4);
        assert_eq!(removed.features.len(), 0);
    }

    #[test]
    fn sort_field_decides_which_point_survives_a_conflict() {
        // Two coincident-ish points; the deeper one must win under descending.
        let mut l = Layer::new("p").with_geom_type(GeometryType::Point);
        l.add_field(FieldDef::new("depth", FieldType::Float));
        l.add_feature(Some(Geometry::point(0.0, 0.0)), &[("depth", 1.0f64.into())])
            .unwrap();
        l.add_feature(Some(Geometry::point(0.1, 0.0)), &[("depth", 9.0f64.into())])
            .unwrap();
        let path = store(l);
        let (kept, _, _) = run(json!({
            "input": path.clone(), "min_distance": 1.0, "sort_field": "depth",
        }));
        assert_eq!(xs(&kept), vec![0.1], "highest depth should win");

        let (kept, _, _) = run(json!({
            "input": path, "min_distance": 1.0, "sort_field": "depth",
            "sort_order": "ascending",
        }));
        assert_eq!(xs(&kept), vec![0.0], "lowest depth should win");
    }

    #[test]
    fn a_flagged_hazard_survives_a_dense_cluster() {
        // The safety case: the shoal is last in input order and has the worst
        // sort value, so only the pin can save it.
        let mut l = Layer::new("p").with_geom_type(GeometryType::Point);
        l.add_field(FieldDef::new("depth", FieldType::Float));
        l.add_field(FieldDef::new("shoal", FieldType::Integer));
        for i in 0..5 {
            l.add_feature(
                Some(Geometry::point(i as f64 * 0.1, 0.0)),
                &[("depth", (9.0 - i as f64).into()), ("shoal", 0i64.into())],
            )
            .unwrap();
        }
        l.add_feature(
            Some(Geometry::point(0.45, 0.0)),
            &[("depth", 0.0f64.into()), ("shoal", 1i64.into())],
        )
        .unwrap();
        let (kept, _, _) = run(json!({
            "input": store(l), "min_distance": 5.0,
            "sort_field": "depth", "keep_field": "shoal",
        }));
        assert_eq!(xs(&kept), vec![0.45], "the pinned hazard must be retained");
    }

    #[test]
    fn bin_keeps_one_point_per_cell() {
        // 0..5 with bin_size 2: cells [0,2) [2,4) [4,6) -> 3 points.
        let (kept, _, res) = run(json!({
            "input": line_points(6), "method": "bin", "bin_size": 2.0,
        }));
        assert_eq!(kept.features.len(), 3);
        assert_eq!(xs(&kept), vec![0.0, 2.0, 4.0]);
        assert_eq!(res.outputs["method"], json!("bin"));
    }

    #[test]
    fn bin_respects_the_sort_priority_within_a_cell() {
        let (kept, _, _) = run(json!({
            "input": line_points(6), "method": "bin", "bin_size": 2.0,
            "sort_field": "depth",
        }));
        // Descending depth: the higher x in each cell wins.
        assert_eq!(xs(&kept), vec![1.0, 3.0, 5.0]);
    }

    #[test]
    fn attributes_and_crs_survive_the_thinning() {
        let (kept, _, _) = run(json!({
            "input": line_points(6), "min_distance": 3.0,
        }));
        assert_eq!(kept.crs_epsg(), Some(3857));
        assert!(kept.schema.field_index("depth").is_some());
        let i = kept.schema.field_index("depth").unwrap();
        assert_eq!(kept.features[0].attributes[i], FieldValue::Float(0.0));
    }

    #[test]
    fn kept_and_removed_partition_the_input_exactly() {
        let (kept, removed, res) = run(json!({
            "input": line_points(20), "min_distance": 4.0,
        }));
        assert_eq!(kept.features.len() + removed.features.len(), 20);
        assert_eq!(res.outputs["input_count"], json!(20));
    }

    #[test]
    fn the_removed_layer_is_produced_without_a_path() {
        let args: ToolArgs = serde_json::from_value(json!({
            "input": line_points(6), "min_distance": 2.0,
        }))
        .unwrap();
        let res = ReducePointDensityTool.run(&args, &ctx()).unwrap();
        let p = res.outputs["output_removed"].as_str().unwrap();
        assert!(load_input_layer(p).is_ok());
    }

    #[test]
    fn a_non_point_layer_is_refused() {
        let mut l = Layer::new("lines").with_geom_type(GeometryType::LineString);
        l.add_feature(
            Some(Geometry::line_string(vec![
                wbvector::Coord::xy(0.0, 0.0),
                wbvector::Coord::xy(1.0, 1.0),
            ])),
            &[],
        )
        .unwrap();
        let args: ToolArgs = serde_json::from_value(json!({
            "input": store(l), "min_distance": 1.0,
        }))
        .unwrap();
        let err = ReducePointDensityTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("point layer"), "{err}");
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ReducePointDensityTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        // Zero or negative would silently retain everything.
        assert!(bad(json!({"input": "p.shp", "min_distance": 0})));
        assert!(bad(json!({"input": "p.shp", "min_distance": -1})));
        assert!(bad(json!({"input": "p.shp"})));
        assert!(bad(json!({"input": "p.shp", "method": "bin"})));
        assert!(bad(json!({"input": "p.shp", "method": "grid", "bin_size": 1})));
    }
}
