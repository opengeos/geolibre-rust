//! GeoLibre tool: Pareto (80/20) concentration of incidents across locations.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Eighty Twenty Analysis* (Crime
//! Analysis and Safety). It answers "what small fraction of locations account
//! for a large fraction of incidents?" — a statistic that the clustering /
//! hot-spot tools (`aggregate_points`, `colocation_analysis`,
//! `emerging_hot_spot_analysis`) do not produce.
//!
//! Steps (all deterministic):
//!   1. Snap coincident / near-coincident incident points into weighted
//!      *locations* using a grid-hashed single-link union-find over the
//!      `cluster_tolerance` distance (mirrors `aggregate_points`' clustering).
//!      Each incident contributes 1, or the value of `weight_field` when given.
//!   2. Sort locations by incident count descending (ties broken by centroid
//!      coordinates for a stable order).
//!   3. Walk the sorted list accumulating a running incident total, tagging
//!      each location with its 1-based `rank` and `cumulative_pct` (share of all
//!      incidents at this location and every higher-ranked one).
//!   4. The shortest prefix whose `cumulative_pct` reaches `threshold`
//!      (default 80%) is the Pareto "head"; those locations are tagged
//!      `band = "head"`, the rest `band = "tail"`.
//!
//! Output is one point per weighted location (placed at the cluster centroid)
//! carrying `location_id`, `incident_count`, `rank`, `cumulative_pct`, and
//! `band`. The run result reports the crossover — e.g. how few locations
//! account for the threshold share of incidents.

use std::collections::{BTreeMap, HashMap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct EightyTwentyAnalysisTool;

impl Tool for EightyTwentyAnalysisTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "eighty_twenty_analysis",
            display_name: "Eighty Twenty Analysis",
            summary: "Snap incident points into weighted locations and compute Pareto concentration ('X% of incidents occur at Y% of locations'), tagging each location with its rank, cumulative-incident percent, and head/tail band, like ArcGIS Eighty Twenty Analysis.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input incident point vector layer (Point or MultiPoint).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output point vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cluster_tolerance",
                    description: "Distance within which incidents snap into one location (map units, default 0 = only exactly coincident points).",
                    required: false,
                },
                ToolParamSpec {
                    name: "weight_field",
                    description: "Optional numeric field giving each incident's weight (default: each incident counts as 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "threshold",
                    description: "Concentration threshold percent that defines the Pareto 'head' (default 80).",
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

        // Resolve the weight field index, if requested.
        let weight_idx =
            match &prm.weight_field {
                Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                    ToolError::Validation(format!("weight_field '{f}' not found"))
                })?),
                None => None,
            };

        // Collect incident points: (x, y, weight).
        let mut pts: Vec<(f64, f64, f64)> = Vec::new();
        for feature in layer.features.iter() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let w = match weight_idx {
                Some(i) => feature.attributes[i].as_f64().unwrap_or(0.0),
                None => 1.0,
            };
            match geom {
                Geometry::Point(c) => pts.push((c.x, c.y, w)),
                Geometry::MultiPoint(cs) => {
                    for c in cs {
                        pts.push((c.x, c.y, w));
                    }
                }
                _ => {}
            }
        }
        if pts.is_empty() {
            return Err(ToolError::Execution(
                "no point features in input".to_string(),
            ));
        }

        ctx.progress
            .info(&format!("snapping {} incident point(s)", pts.len()));

        // ── Grid-hashed single-link clustering (union-find) ──────────────────────
        // A cell size == tolerance means only points within one cell (and its 8
        // neighbours) are compared; tolerance 0 collapses only exact duplicates.
        let tol = prm.cluster_tolerance;
        let cell = tol.max(1e-9);
        let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
        for (i, (x, y, _)) in pts.iter().enumerate() {
            grid.entry(((x / cell).floor() as i64, (y / cell).floor() as i64))
                .or_default()
                .push(i);
        }
        let mut uf = UnionFind::new(pts.len());
        let tol2 = tol * tol;
        for (i, (x, y, _)) in pts.iter().enumerate() {
            let (gx, gy) = ((x / cell).floor() as i64, (y / cell).floor() as i64);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    if let Some(bucket) = grid.get(&(gx + dx, gy + dy)) {
                        for &j in bucket {
                            if j <= i {
                                continue;
                            }
                            let (xj, yj, _) = pts[j];
                            if (x - xj).powi(2) + (y - yj).powi(2) <= tol2 {
                                uf.union(i, j);
                            }
                        }
                    }
                }
            }
        }

        // Gather clusters by root -> member point indices.
        let mut clusters: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..pts.len() {
            clusters.entry(uf.find(i)).or_default().push(i);
        }

        // Reduce each cluster to a weighted location at its centroid.
        let mut locations: Vec<Location> = clusters
            .values()
            .map(|members| {
                let n = members.len() as f64;
                let sx: f64 = members.iter().map(|&i| pts[i].0).sum();
                let sy: f64 = members.iter().map(|&i| pts[i].1).sum();
                let count: f64 = members.iter().map(|&i| pts[i].2).sum();
                Location {
                    x: sx / n,
                    y: sy / n,
                    count,
                }
            })
            .collect();

        let total: f64 = locations.iter().map(|l| l.count).sum();
        if total <= 0.0 {
            return Err(ToolError::Execution(
                "total incident weight is zero; check 'weight_field' values".to_string(),
            ));
        }

        // Sort by incident count descending; break ties by coordinates so the
        // ranking (and therefore every cumulative value) is deterministic.
        locations.sort_by(|a, b| {
            b.count
                .partial_cmp(&a.count)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.x.partial_cmp(&b.x).unwrap_or(std::cmp::Ordering::Equal))
                .then(a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal))
        });

        // ── Build the output point layer ─────────────────────────────────────────
        let mut out = Layer::new("eighty_twenty_locations").with_geom_type(GeometryType::Point);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("location_id", FieldType::Integer));
        out.add_field(FieldDef::new("incident_count", FieldType::Float));
        out.add_field(FieldDef::new("rank", FieldType::Integer));
        out.add_field(FieldDef::new("cumulative_pct", FieldType::Float));
        out.add_field(FieldDef::new("band", FieldType::Text));

        let mut cumulative = 0.0f64;
        let mut head_count = 0usize;
        let mut reached = false; // has cumulative_pct crossed the threshold yet?
        for (i, loc) in locations.iter().enumerate() {
            cumulative += loc.count;
            let cumulative_pct = cumulative / total * 100.0;
            // The head is the shortest prefix reaching the threshold: every
            // location up to and including the first that crosses it.
            let band = if !reached { "head" } else { "tail" };
            if band == "head" {
                head_count += 1;
            }
            if cumulative_pct >= prm.threshold {
                reached = true;
            }
            out.add_feature(
                Some(Geometry::Point(Coord::xy(loc.x, loc.y))),
                &[
                    ("location_id", FieldValue::Integer(i as i64)),
                    ("incident_count", FieldValue::Float(loc.count)),
                    ("rank", FieldValue::Integer(i as i64 + 1)),
                    ("cumulative_pct", FieldValue::Float(cumulative_pct)),
                    ("band", FieldValue::Text(band.to_string())),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed writing location: {e}")))?;
        }

        let location_count = locations.len();
        let head_location_pct = head_count as f64 / location_count as f64 * 100.0;
        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_points".to_string(), json!(pts.len()));
        outputs.insert("location_count".to_string(), json!(location_count));
        outputs.insert("total_incidents".to_string(), json!(total));
        outputs.insert("threshold_pct".to_string(), json!(prm.threshold));
        outputs.insert("head_location_count".to_string(), json!(head_count));
        outputs.insert("head_location_pct".to_string(), json!(head_location_pct));
        outputs.insert(
            "pareto_summary".to_string(),
            json!(format!(
                "{:.0}% of incidents occur at {:.1}% of locations ({} of {})",
                prm.threshold, head_location_pct, head_count, location_count
            )),
        );
        Ok(ToolRunResult { outputs })
    }
}

/// A weighted location: cluster centroid and its total incident count.
struct Location {
    x: f64,
    y: f64,
    count: f64,
}

// ── Union-Find ─────────────────────────────────────────────────────────────────

struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn new(n: usize) -> UnionFind {
        UnionFind {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }
    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        match self.rank[ra].cmp(&self.rank[rb]) {
            std::cmp::Ordering::Less => self.parent[ra] = rb,
            std::cmp::Ordering::Greater => self.parent[rb] = ra,
            std::cmp::Ordering::Equal => {
                self.parent[rb] = ra;
                self.rank[ra] += 1;
            }
        }
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    cluster_tolerance: f64,
    weight_field: Option<String>,
    threshold: f64,
}

/// Parses a numeric parameter accepting a JSON number or a numeric string.
fn parse_num(args: &ToolArgs, key: &str, default: f64) -> Result<f64, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(default),
        Some(Value::Number(n)) => Ok(n.as_f64().unwrap_or(default)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(default),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation(format!("'{key}' must be a number"))),
        Some(_) => Err(ToolError::Validation(format!("'{key}' must be a number"))),
    }
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let cluster_tolerance = parse_num(args, "cluster_tolerance", 0.0)?;
    if cluster_tolerance.is_nan() || cluster_tolerance < 0.0 {
        return Err(ToolError::Validation(
            "'cluster_tolerance' must be zero or positive".to_string(),
        ));
    }
    let threshold = parse_num(args, "threshold", 80.0)?;
    if threshold.is_nan() || threshold <= 0.0 || threshold > 100.0 {
        return Err(ToolError::Validation(
            "'threshold' must be in the range (0, 100]".to_string(),
        ));
    }
    let weight_field = parse_optional_str(args, "weight_field")?.map(str::to_string);
    Ok(Params {
        cluster_tolerance,
        weight_field,
        threshold,
    })
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

    /// Builds an in-memory point layer of (x, y, weight) rows.
    fn layer_of(pts: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("incidents")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("w", FieldType::Float));
        for (x, y, w) in pts {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(*x, *y))),
                &[("w", (*w).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = EightyTwentyAnalysisTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// Weighted counts 70/10/10/10 (total 100): at threshold 80 the head is the
    /// two highest-ranked locations (70% then 80%), i.e. 50% of locations.
    #[test]
    fn pareto_head_crossover() {
        let input = layer_of(&[
            (0.0, 0.0, 70.0),
            (100.0, 0.0, 10.0),
            (200.0, 0.0, 10.0),
            (300.0, 0.0, 10.0),
        ]);
        let (out, layer) = run(json!({
            "input": input, "weight_field": "w", "threshold": 80.0
        }));
        assert_eq!(out.outputs["location_count"], json!(4));
        assert_eq!(out.outputs["total_incidents"], json!(100.0));
        assert_eq!(out.outputs["head_location_count"], json!(2));
        assert!(
            (out.outputs["head_location_pct"].as_f64().unwrap() - 50.0).abs() < 1e-9,
            "expected 50% of locations in the head"
        );
        // Rank 1 is the dominant location; its cumulative share is 70%.
        let rank_idx = layer.schema.field_index("rank").unwrap();
        let cum_idx = layer.schema.field_index("cumulative_pct").unwrap();
        let band_idx = layer.schema.field_index("band").unwrap();
        let mut by_rank: Vec<_> = layer.iter().collect();
        by_rank.sort_by_key(|f| f.attributes[rank_idx].as_i64().unwrap());
        assert!((by_rank[0].attributes[cum_idx].as_f64().unwrap() - 70.0).abs() < 1e-9);
        assert!((by_rank[1].attributes[cum_idx].as_f64().unwrap() - 80.0).abs() < 1e-9);
        assert_eq!(by_rank[0].attributes[band_idx].as_str().unwrap(), "head");
        assert_eq!(by_rank[1].attributes[band_idx].as_str().unwrap(), "head");
        assert_eq!(by_rank[2].attributes[band_idx].as_str().unwrap(), "tail");
        // Final cumulative percent must reach 100.
        assert!((by_rank[3].attributes[cum_idx].as_f64().unwrap() - 100.0).abs() < 1e-9);
    }

    /// Coincident incidents (no weight field) snap into one weighted location
    /// whose count is the number of incidents there.
    #[test]
    fn coincident_points_snap_and_count() {
        // Three incidents at the origin, one elsewhere.
        let input = layer_of(&[
            (0.0, 0.0, 1.0),
            (0.0, 0.0, 1.0),
            (0.0, 0.0, 1.0),
            (500.0, 500.0, 1.0),
        ]);
        let (out, layer) = run(json!({ "input": input })); // tolerance defaults to 0
        assert_eq!(out.outputs["input_points"], json!(4));
        assert_eq!(out.outputs["location_count"], json!(2));
        assert_eq!(out.outputs["total_incidents"], json!(4.0));
        let cnt_idx = layer.schema.field_index("incident_count").unwrap();
        let rank_idx = layer.schema.field_index("rank").unwrap();
        let top = layer
            .iter()
            .find(|f| f.attributes[rank_idx].as_i64().unwrap() == 1)
            .unwrap();
        assert!((top.attributes[cnt_idx].as_f64().unwrap() - 3.0).abs() < 1e-9);
    }

    /// Distinct, non-coincident points at tolerance 0 stay separate: one
    /// location per incident, and the ranking is a proper 1..=n permutation.
    #[test]
    fn distinct_points_pass_through() {
        let input = layer_of(&[(0.0, 0.0, 1.0), (100.0, 0.0, 1.0), (200.0, 0.0, 1.0)]);
        let (out, layer) = run(json!({ "input": input }));
        assert_eq!(out.outputs["location_count"], json!(3));
        let rank_idx = layer.schema.field_index("rank").unwrap();
        let mut ranks: Vec<i64> = layer
            .iter()
            .map(|f| f.attributes[rank_idx].as_i64().unwrap())
            .collect();
        ranks.sort_unstable();
        assert_eq!(ranks, vec![1, 2, 3]);
    }

    /// A tolerance that spans the gap merges near points into one location.
    #[test]
    fn tolerance_snaps_near_points() {
        let input = layer_of(&[(0.0, 0.0, 1.0), (5.0, 0.0, 1.0), (1000.0, 0.0, 1.0)]);
        let (out, _l) = run(json!({ "input": input, "cluster_tolerance": 20.0 }));
        assert_eq!(out.outputs["location_count"], json!(2));
    }

    #[test]
    fn rejects_bad_parameters() {
        let input = layer_of(&[(0.0, 0.0, 1.0)]);
        // Threshold out of range.
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "threshold": 150.0 })).unwrap();
        assert!(EightyTwentyAnalysisTool.validate(&args).is_err());
        // Negative tolerance.
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "cluster_tolerance": -1.0 })).unwrap();
        assert!(EightyTwentyAnalysisTool.validate(&args).is_err());
        // Missing input.
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(EightyTwentyAnalysisTool.validate(&args).is_err());
    }

    #[test]
    fn rejects_unknown_weight_field() {
        let input = layer_of(&[(0.0, 0.0, 1.0)]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "weight_field": "nope" })).unwrap();
        assert!(EightyTwentyAnalysisTool.run(&args, &ctx()).is_err());
    }
}
