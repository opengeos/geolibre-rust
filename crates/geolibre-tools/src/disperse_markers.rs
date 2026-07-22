//! GeoLibre tool: disperse coincident/overlapping point markers.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Disperse Markers* (Cartography).
//! Point datasets frequently carry exact or near-duplicate coordinates
//! (multiple incidents at one address, co-located wells); at map scale the
//! symbols stack and only the top one is visible. This tool nudges each
//! coincident/overlapping group apart to a minimum center-to-center spacing
//! while leaving well-separated points and non-point geometries untouched.
//!
//! ## Clustering
//! Points closer than `min_spacing` are grouped with a grid-hashed
//! union-find — the same single-link idiom `aggregate_points` uses. A
//! bucket grid sized to `min_spacing` keeps the neighbour search close to
//! linear for the modestly sized clusters this tool targets; a full kdtree
//! would only pay off for very large single clusters, which are out of
//! scope for v1.
//!
//! ## Placement
//! Each cluster's centroid is held fixed and members are laid out on one of
//! four patterns, then re-centred so the new centroid exactly matches the
//! original:
//! * `square` — a centred grid with `min_spacing` row/column pitch (every
//!   pair is at least `min_spacing` apart: axis-aligned neighbours are
//!   exactly `min_spacing`, diagonal neighbours are `min_spacing * sqrt(2)`).
//! * `cross` — four cardinal arms (N/E/S/W); members fill the arms in
//!   round-robin, radius growing by `min_spacing` every 4th member.
//! * `ring` — one circle sized so adjacent members are exactly
//!   `min_spacing` apart; since chord length is maximal between angular
//!   neighbours on an evenly spaced circle, that also bounds every other
//!   pair.
//! * `expanded` (default) — concentric rings filled outward using the same
//!   chord sizing as `ring`, alternating each ring's start angle, then a
//!   bounded (60 iterations, O(n^2) each) local-repulsion pass nudges any
//!   pair still under `min_spacing` apart until the cluster clears it or the
//!   iteration budget runs out. Pathologically large single clusters (e.g.
//!   thousands of exactly-coincident points) may not fully converge within
//!   the budget; `square`/`cross`/`ring` give an exact guarantee for any
//!   cluster size and are the safer choice there.
//!
//! Singleton points (no neighbour within `min_spacing`) and non-point
//! geometries pass through untouched. `orig_x`/`orig_y` preserve the
//! pre-dispersal position and `displaced` records the move distance (0 for
//! untouched points, null for non-point features). Feature order in the
//! output mirrors the input; clustering itself processes points sorted by
//! feature id for determinism. `seed` only seeds a tiny inline hash used to
//! pick each cluster's starting angle and to break ties when members start
//! out exactly coincident (direction otherwise undefined) — no wall-clock or
//! unseeded randomness anywhere.

use std::collections::{BTreeMap, HashMap};
use std::f64::consts::{PI, TAU};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

#[derive(Clone, Copy, PartialEq)]
enum Pattern {
    Expanded,
    Ring,
    Cross,
    Square,
}

pub struct DisperseMarkersTool;

impl Tool for DisperseMarkersTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "disperse_markers",
            display_name: "Disperse Markers",
            summary: "Spread coincident or overlapping point symbols apart to a minimum center-to-center spacing (expanded/ring/cross/square patterns) so stacked points stay individually visible, like ArcGIS Disperse Markers.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output point vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_spacing",
                    description: "Minimum center-to-center spacing to enforce between clustered points, in map units.",
                    required: true,
                },
                ToolParamSpec {
                    name: "pattern",
                    description: "'expanded' (default, concentric rings + local repulsion), 'ring' (single circle), 'cross' (four arms), or 'square' (grid).",
                    required: false,
                },
                ToolParamSpec {
                    name: "seed",
                    description: "Integer seed for deterministic tie-breaking when members start out exactly coincident (default 0).",
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

        // Collect point features (fid, feature index, x, y), sorted by fid
        // for deterministic clustering.
        let mut pts: Vec<(u64, usize, f64, f64)> = layer
            .features
            .iter()
            .enumerate()
            .filter_map(|(fidx, f)| match f.geometry.as_ref() {
                Some(Geometry::Point(c)) => Some((f.fid, fidx, c.x, c.y)),
                _ => None,
            })
            .collect();
        pts.sort_by_key(|&(fid, fidx, _, _)| (fid, fidx));

        ctx.progress
            .info(&format!("dispersing {} point(s)", pts.len()));

        // Grid-hashed single-link clustering (union-find), same idiom as
        // `aggregate_points`. Points strictly closer than `min_spacing` join
        // the same cluster.
        let cell = prm.min_spacing.max(1e-9);
        let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
        for (i, &(_, _, x, y)) in pts.iter().enumerate() {
            grid.entry(((x / cell).floor() as i64, (y / cell).floor() as i64))
                .or_default()
                .push(i);
        }
        let mut uf = UnionFind::new(pts.len());
        let d2 = prm.min_spacing * prm.min_spacing;
        for (i, &(_, _, x, y)) in pts.iter().enumerate() {
            let (gx, gy) = ((x / cell).floor() as i64, (y / cell).floor() as i64);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    if let Some(bucket) = grid.get(&(gx + dx, gy + dy)) {
                        for &j in bucket {
                            if j <= i {
                                continue;
                            }
                            let (_, _, xj, yj) = pts[j];
                            if (x - xj).powi(2) + (y - yj).powi(2) < d2 {
                                uf.union(i, j);
                            }
                        }
                    }
                }
            }
        }

        let mut clusters: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for i in 0..pts.len() {
            clusters.entry(uf.find(i)).or_default().push(i);
        }

        // New (x, y) per `pts` index; defaults to the original position
        // (singletons stay untouched).
        let mut new_pos: Vec<(f64, f64)> = pts.iter().map(|&(_, _, x, y)| (x, y)).collect();
        let mut cluster_count = 0usize;
        let mut moved_points = 0usize;
        for members in clusters.values() {
            if members.len() < 2 {
                continue;
            }
            cluster_count += 1;
            moved_points += members.len();
            let cluster_key = pts[members[0]].0; // smallest fid in the cluster (members is fid-sorted)
            let cx = members.iter().map(|&i| pts[i].2).sum::<f64>() / members.len() as f64;
            let cy = members.iter().map(|&i| pts[i].3).sum::<f64>() / members.len() as f64;

            let mut offsets = place_cluster(
                prm.pattern,
                members.len(),
                prm.min_spacing,
                prm.seed,
                cluster_key,
            );
            // Re-centre so the new arrangement's centroid exactly matches the
            // original centroid, regardless of any placement asymmetry.
            let mean_x = offsets.iter().map(|p| p.0).sum::<f64>() / offsets.len() as f64;
            let mean_y = offsets.iter().map(|p| p.1).sum::<f64>() / offsets.len() as f64;
            for (ox, oy) in offsets.iter_mut() {
                *ox -= mean_x;
                *oy -= mean_y;
            }
            for (k, &i) in members.iter().enumerate() {
                new_pos[i] = (cx + offsets[k].0, cy + offsets[k].1);
            }
        }

        // Map fid-sorted `pts` results back onto original feature index.
        let mut per_feature: Vec<Option<(f64, f64, f64, f64)>> = vec![None; layer.features.len()];
        for (i, &(_, fidx, ox, oy)) in pts.iter().enumerate() {
            let (nx, ny) = new_pos[i];
            per_feature[fidx] = Some((ox, oy, nx, ny));
        }

        // Build the output layer: input schema + orig_x/orig_y/displaced.
        let mut out = wbvector::Layer::new(layer.name.clone());
        out.geom_type = layer.geom_type;
        out.crs = layer.crs.clone();
        out.schema = layer.schema.clone();
        out.add_field(FieldDef::new("orig_x", FieldType::Float));
        out.add_field(FieldDef::new("orig_y", FieldType::Float));
        out.add_field(FieldDef::new("displaced", FieldType::Float));

        let mut non_point_features = 0usize;
        for (fidx, feature) in layer.features.iter().enumerate() {
            let mut f = feature.clone();
            match per_feature[fidx] {
                Some((ox, oy, nx, ny)) => {
                    let displaced = ((nx - ox).powi(2) + (ny - oy).powi(2)).sqrt();
                    f.geometry = Some(Geometry::point(nx, ny));
                    f.attributes.push(FieldValue::Float(ox));
                    f.attributes.push(FieldValue::Float(oy));
                    f.attributes.push(FieldValue::Float(displaced));
                }
                None => {
                    non_point_features += 1;
                    f.attributes.push(FieldValue::Null);
                    f.attributes.push(FieldValue::Null);
                    f.attributes.push(FieldValue::Null);
                }
            }
            out.push(f);
        }

        let feature_count = out.len();
        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("input_points".to_string(), json!(pts.len()));
        outputs.insert("non_point_features".to_string(), json!(non_point_features));
        outputs.insert("cluster_count".to_string(), json!(cluster_count));
        outputs.insert("moved_points".to_string(), json!(moved_points));
        Ok(ToolRunResult { outputs })
    }
}

/// Offsets (relative to an arbitrary common origin) for `n` cluster members
/// on the chosen pattern, sized against `spacing`.
fn place_cluster(
    pattern: Pattern,
    n: usize,
    spacing: f64,
    seed: u64,
    cluster_key: u64,
) -> Vec<(f64, f64)> {
    let angle0 = cluster_angle_offset(seed, cluster_key);
    match pattern {
        Pattern::Square => square_offsets(n, spacing),
        Pattern::Cross => cross_offsets(n, spacing),
        Pattern::Ring => ring_offsets(n, spacing, angle0),
        Pattern::Expanded => {
            let mut offsets = expanded_offsets(n, spacing, angle0);
            relax(&mut offsets, spacing, seed, cluster_key);
            offsets
        }
    }
}

/// Centred grid, `spacing` pitch. Axis-aligned neighbours are exactly
/// `spacing` apart; diagonal neighbours are `spacing * sqrt(2)`.
fn square_offsets(n: usize, spacing: f64) -> Vec<(f64, f64)> {
    let cols = (n as f64).sqrt().ceil().max(1.0) as usize;
    let rows = n.div_ceil(cols);
    (0..n)
        .map(|i| {
            let row = i / cols;
            let col = i % cols;
            let x = (col as f64 - (cols as f64 - 1.0) / 2.0) * spacing;
            let y = (row as f64 - (rows as f64 - 1.0) / 2.0) * spacing;
            (x, y)
        })
        .collect()
}

const ARMS: [(f64, f64); 4] = [(0.0, 1.0), (1.0, 0.0), (0.0, -1.0), (-1.0, 0.0)];

/// Four cardinal arms filled round-robin; radius grows by `spacing` every
/// 4th member, so same-arm neighbours are exactly `spacing` apart and
/// perpendicular same-radius neighbours are `spacing * sqrt(2)`.
fn cross_offsets(n: usize, spacing: f64) -> Vec<(f64, f64)> {
    (0..n)
        .map(|i| {
            let ring = (i / 4 + 1) as f64;
            let (dx, dy) = ARMS[i % 4];
            (dx * ring * spacing, dy * ring * spacing)
        })
        .collect()
}

/// A single circle sized so adjacent members are exactly `spacing` apart.
fn ring_offsets(n: usize, spacing: f64, angle0: f64) -> Vec<(f64, f64)> {
    if n <= 1 {
        return vec![(0.0, 0.0); n];
    }
    let radius = spacing / (2.0 * (PI / n as f64).sin());
    (0..n)
        .map(|i| {
            let theta = angle0 + TAU * i as f64 / n as f64;
            (radius * theta.cos(), radius * theta.sin())
        })
        .collect()
}

/// Concentric rings filled outward: ring `k` (radius `k * spacing`) holds as
/// many evenly spaced members as fit while keeping adjacent chord length
/// `>= spacing`. Alternates each ring's start angle to reduce radial
/// alignment between consecutive rings.
fn expanded_offsets(n: usize, spacing: f64, angle0: f64) -> Vec<(f64, f64)> {
    let mut offsets = Vec::with_capacity(n);
    let mut remaining = n;
    let mut ring_idx = 1u32;
    while remaining > 0 {
        let radius = ring_idx as f64 * spacing;
        // radius >= spacing (ring_idx >= 1) so this ratio is always <= 0.5.
        let ratio = spacing / (2.0 * radius);
        let cap = if ratio >= 1.0 {
            1
        } else {
            ((PI / ratio.asin()).floor() as usize).max(1)
        };
        let take = cap.min(remaining);
        let stagger = if ring_idx.is_multiple_of(2) {
            PI / take as f64
        } else {
            0.0
        };
        for i in 0..take {
            let theta = angle0 + stagger + TAU * i as f64 / take as f64;
            offsets.push((radius * theta.cos(), radius * theta.sin()));
        }
        remaining -= take;
        ring_idx += 1;
    }
    offsets
}

/// Bounded local-repulsion refinement: while any pair is still closer than
/// `spacing`, push each such pair apart by half the deficit. Deterministic;
/// exactly coincident points get a hash-derived push direction instead of
/// dividing by zero.
fn relax(offsets: &mut [(f64, f64)], spacing: f64, seed: u64, cluster_key: u64) {
    let n = offsets.len();
    if n < 2 {
        return;
    }
    const MAX_ITERS: u32 = 60;
    for _ in 0..MAX_ITERS {
        let mut moved = false;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = offsets[j].0 - offsets[i].0;
                let dy = offsets[j].1 - offsets[i].1;
                let dist = (dx * dx + dy * dy).sqrt();
                if dist < spacing - 1e-9 {
                    moved = true;
                    let (ux, uy) = if dist > 1e-9 {
                        (dx / dist, dy / dist)
                    } else {
                        let theta = tie_break_angle(seed, cluster_key, i as u64, j as u64);
                        (theta.cos(), theta.sin())
                    };
                    let push = (spacing - dist) / 2.0 + 1e-9;
                    offsets[i].0 -= ux * push;
                    offsets[i].1 -= uy * push;
                    offsets[j].0 += ux * push;
                    offsets[j].1 += uy * push;
                }
            }
        }
        if !moved {
            break;
        }
    }
}

/// A tiny deterministic 64-bit mix (splitmix64-style) — no RNG crate needed,
/// no wall-clock, reproducible for a given `seed`.
fn hash64(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d049bb133111eb);
    x ^= x >> 31;
    x
}

fn cluster_angle_offset(seed: u64, cluster_key: u64) -> f64 {
    let h = hash64(seed ^ hash64(cluster_key));
    (h as f64 / u64::MAX as f64) * TAU
}

fn tie_break_angle(seed: u64, cluster_key: u64, i: u64, j: u64) -> f64 {
    let h = hash64(seed ^ hash64(cluster_key ^ (i << 32) ^ j));
    (h as f64 / u64::MAX as f64) * TAU
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
    min_spacing: f64,
    pattern: Pattern,
    seed: u64,
}

fn parse_optional_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64().or_else(|| n.as_f64().map(|f| f as u64))),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be an integer"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be an integer"
        ))),
    }
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let min_spacing = match args.get("min_spacing") {
        Some(Value::Number(n)) => n.as_f64().unwrap_or(0.0),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'min_spacing' must be a number".into()))?,
        _ => {
            return Err(ToolError::Validation(
                "missing required numeric parameter 'min_spacing'".to_string(),
            ))
        }
    };
    if min_spacing.is_nan() || min_spacing <= 0.0 {
        return Err(ToolError::Validation(
            "'min_spacing' must be positive".to_string(),
        ));
    }
    let pattern = match parse_optional_str(args, "pattern")?.map(|s| s.trim().to_lowercase()) {
        None => Pattern::Expanded,
        Some(s) if s.is_empty() || s == "expanded" => Pattern::Expanded,
        Some(s) if s == "ring" => Pattern::Ring,
        Some(s) if s == "cross" => Pattern::Cross,
        Some(s) if s == "square" => Pattern::Square,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "'pattern' must be one of 'expanded', 'ring', 'cross', 'square', got '{other}'"
            )))
        }
    };
    let seed = parse_optional_u64(args, "seed")?.unwrap_or(0);
    Ok(Params {
        min_spacing,
        pattern,
        seed,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Coord, FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn layer_of_points(pts: &[(f64, f64)]) -> String {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("val", FieldType::Integer));
        for (i, (x, y)) in pts.iter().enumerate() {
            l.add_feature(Some(Geometry::point(*x, *y)), &[("val", (i as i64).into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, wbvector::Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = DisperseMarkersTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn min_pairwise_dist(pts: &[(f64, f64)]) -> f64 {
        let mut m = f64::INFINITY;
        for i in 0..pts.len() {
            for j in (i + 1)..pts.len() {
                let d = ((pts[i].0 - pts[j].0).powi(2) + (pts[i].1 - pts[j].1).powi(2)).sqrt();
                m = m.min(d);
            }
        }
        m
    }

    fn point_xy(layer: &wbvector::Layer) -> Vec<(f64, f64)> {
        layer
            .iter()
            .map(|f| match f.geometry.as_ref().unwrap() {
                Geometry::Point(c) => (c.x, c.y),
                other => panic!("expected point, got {other:?}"),
            })
            .collect()
    }

    /// A cluster of exactly-coincident points separates so every pair clears
    /// `min_spacing`, for every pattern.
    #[test]
    fn overlapping_points_separate_to_min_spacing() {
        for pattern in ["expanded", "ring", "cross", "square"] {
            let pts = [(10.0, 10.0); 7];
            let input = layer_of_points(&pts);
            let (out, layer) = run(json!({
                "input": input, "min_spacing": 5.0, "pattern": pattern
            }));
            assert_eq!(out.outputs["cluster_count"], json!(1));
            assert_eq!(out.outputs["moved_points"], json!(7));
            let after = point_xy(&layer);
            let min_after = min_pairwise_dist(&after);
            assert!(
                min_after >= 5.0 - 1e-6,
                "pattern {pattern}: min pairwise distance {min_after} < 5.0"
            );
            // Point count is preserved.
            assert_eq!(after.len(), 7);
        }
    }

    /// A near-coincident (but not exactly identical) cluster also separates.
    #[test]
    fn near_coincident_points_separate() {
        let pts = [(0.0, 0.0), (0.01, 0.0), (0.0, 0.01), (0.01, 0.01)];
        let before_min = min_pairwise_dist(&pts);
        assert!(before_min < 10.0);
        let input = layer_of_points(&pts);
        let (_out, layer) = run(json!({ "input": input, "min_spacing": 10.0 }));
        let after = point_xy(&layer);
        assert!(min_pairwise_dist(&after) >= 10.0 - 1e-6);
    }

    /// A lone point far from everything else is left exactly where it was.
    #[test]
    fn singleton_untouched() {
        let pts = [(0.0, 0.0), (5000.0, 5000.0)];
        let input = layer_of_points(&pts);
        let (out, layer) = run(json!({ "input": input, "min_spacing": 5.0 }));
        assert_eq!(out.outputs["cluster_count"], json!(0));
        let after = point_xy(&layer);
        assert_eq!(after, pts);
        let didx = layer.schema.field_index("displaced").unwrap();
        for f in layer.iter() {
            assert_eq!(f.attributes[didx].as_f64().unwrap(), 0.0);
        }
    }

    /// Non-point geometries pass through unchanged (geometry and attributes).
    #[test]
    fn non_point_geometry_passes_through() {
        let mut l = Layer::new("mixed")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("val", FieldType::Integer));
        l.add_feature(
            Some(Geometry::line_string(vec![
                Coord::xy(0.0, 0.0),
                Coord::xy(1.0, 1.0),
            ])),
            &[("val", 1i64.into())],
        )
        .unwrap();
        l.add_feature(Some(Geometry::point(0.0, 0.0)), &[("val", 2i64.into())])
            .unwrap();
        l.add_feature(Some(Geometry::point(0.0, 0.0)), &[("val", 3i64.into())])
            .unwrap();
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run(json!({ "input": input, "min_spacing": 5.0 }));
        assert_eq!(out.outputs["non_point_features"], json!(1));
        assert_eq!(out.outputs["cluster_count"], json!(1));
        let line_feature = layer.iter().next().unwrap();
        match line_feature.geometry.as_ref().unwrap() {
            Geometry::LineString(cs) => {
                assert_eq!(cs, &vec![Coord::xy(0.0, 0.0), Coord::xy(1.0, 1.0)])
            }
            other => panic!("expected line string, got {other:?}"),
        }
        let oidx = layer.schema.field_index("orig_x").unwrap();
        assert!(line_feature.attributes[oidx].is_null());
    }

    #[test]
    fn rejects_bad_parameters() {
        let input = layer_of_points(&[(0.0, 0.0)]);
        // Missing min_spacing.
        let args: ToolArgs = serde_json::from_value(json!({ "input": input })).unwrap();
        assert!(DisperseMarkersTool.validate(&args).is_err());
        // Non-positive min_spacing.
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "min_spacing": 0.0 })).unwrap();
        assert!(DisperseMarkersTool.validate(&args).is_err());
        // Unknown pattern.
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input, "min_spacing": 5.0, "pattern": "spiral"
        }))
        .unwrap();
        assert!(DisperseMarkersTool.validate(&args).is_err());
    }
}
