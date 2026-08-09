//! GeoLibre tool: spatial weights from network travel distance.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Generate Network Spatial Weights*
//! (Spatial Statistics).
//!
//! ## The gap
//!
//! `generate_spatial_weights_matrix` offers `knn`, `fixed_distance_band`,
//! `inverse_distance`, `contiguity_edges`, `contiguity_edges_corners` and
//! `delaunay` — every one of them Euclidean or topological. For anything
//! constrained to a network (crime along streets, retail catchments, stream
//! ecology, accessibility) straight-line weights connect features that cannot
//! actually reach each other: two addresses either side of a motorway are
//! 30 m apart and 3 km apart, and the Euclidean matrix says 30 m.
//!
//! That bias then propagates into every statistic built on the matrix —
//! `global_morans_i`, `local_morans_i_lisa`, `getis_ord_gi_star`.
//!
//! The bundled network suite (`network_od_cost_matrix`, `shortest_path_network`,
//! `network_service_area`) can produce the distances but nothing turns them
//! into a weights matrix.
//!
//! ## Interchange
//!
//! The output uses **exactly** the schema `generate_spatial_weights_matrix`
//! emits (`origin`, `neighbor`, `weight`), so every existing consumer —
//! including `compare_spatial_weights` — accepts it with no change. That is the
//! point: this is a new way to build the same object, not a new object.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap};

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::args_common::{bool_or, choice_or, f64_or, opt_positive_f64, opt_usize, req_str};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct GenerateNetworkSwmTool;

impl Tool for GenerateNetworkSwmTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "generate_network_swm",
            display_name: "Generate Network SWM",
            summary: "Builds a spatial weights matrix in which neighbours and weights come from travel distance along a line network rather than straight-line geometry (ArcGIS Generate Network Spatial Weights). generate_spatial_weights_matrix offers only Euclidean and topological conceptualizations (knn, fixed distance band, inverse distance, contiguity, Delaunay), so for network-constrained phenomena it links features that cannot reach each other and biases every downstream statistic. Output uses the same table schema, so global_morans_i, local_morans_i_lisa, getis_ord_gi_star and compare_spatial_weights consume it unchanged.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Point features to weight.",
                    required: true,
                },
                ToolParamSpec {
                    name: "network",
                    description: "Line network the travel distances are measured along.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Weights table (CSV by extension, or a geometry-less vector table), in the same origin/neighbor/weight schema generate_spatial_weights_matrix emits. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "id_field",
                    description: "Field holding each feature's unique id. Default: the 0-based feature index, matching generate_spatial_weights_matrix.",
                    required: false,
                },
                ToolParamSpec {
                    name: "impedance_field",
                    description: "Network cost field. Default: each segment's geometric length.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distance_cutoff",
                    description: "Maximum network impedance a neighbour may lie at. Default: unlimited.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_neighbors",
                    description: "Cap on neighbours per origin, keeping the nearest. Default: unlimited.",
                    required: false,
                },
                ToolParamSpec {
                    name: "conceptualization",
                    description: "'fixed' (default, binary weight 1 within the cutoff) or 'inverse' (weight = 1 / impedance^exponent).",
                    required: false,
                },
                ToolParamSpec {
                    name: "exponent",
                    description: "Decay exponent for 'inverse' (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "row_standardization",
                    description: "Divide each origin's weights by their row sum so every origin totals 1. Default false, matching generate_spatial_weights_matrix.",
                    required: false,
                },
                ToolParamSpec {
                    name: "snap_tolerance",
                    description: "Distance within which two network vertices are welded into one node, in CRS units (default 1e-6). A network whose segments almost meet would otherwise disconnect silently, and every affected origin would emit zero rows.",
                    required: false,
                },
                ToolParamSpec {
                    name: "search_tolerance",
                    description: "Maximum snapping distance from a point to the network. Points further away are reported as unsnapped rather than silently attached to a distant node.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        req_str(args, "network")?;
        choice_or(args, "conceptualization", &["fixed", "inverse"], "fixed")?;
        opt_positive_f64(args, "distance_cutoff")?;
        opt_positive_f64(args, "search_tolerance")?;
        opt_positive_f64(args, "snap_tolerance")?;
        let e = f64_or(args, "exponent", 1.0)?;
        if e <= 0.0 {
            return Err(ToolError::Validation("'exponent' must be > 0".to_string()));
        }
        bool_or(args, "row_standardization", false)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?;
        let network_path = req_str(args, "network")?;
        let output = parse_optional_str(args, "output")?;
        let inverse =
            choice_or(args, "conceptualization", &["fixed", "inverse"], "fixed")? == "inverse";
        let exponent = f64_or(args, "exponent", 1.0)?;
        let cutoff = opt_positive_f64(args, "distance_cutoff")?.unwrap_or(f64::INFINITY);
        let max_neighbors = opt_usize(args, "max_neighbors")?.unwrap_or(usize::MAX);
        let row_standardize = bool_or(args, "row_standardization", false)?;
        let tolerance = opt_positive_f64(args, "search_tolerance")?.unwrap_or(f64::INFINITY);
        let snap = opt_positive_f64(args, "snap_tolerance")?.unwrap_or(1e-6);

        let network = load_input_layer(network_path)?;
        let impedance_field = parse_optional_str(args, "impedance_field")?;
        let mut graph = Graph::build(&network, impedance_field, snap)?;
        if graph.nodes.is_empty() {
            return Err(ToolError::Execution(
                "'network' holds no line geometry".to_string(),
            ));
        }

        let points = load_input_layer(input)?;
        let id_field = parse_optional_str(args, "id_field")?;
        let id_idx = match id_field {
            Some(f) => Some(points.schema.field_index(f).ok_or_else(|| {
                ToolError::Validation(format!("id_field '{f}' not found on the input layer"))
            })?),
            None => None,
        };

        // Snap each point to its nearest network node.
        let mut ids: Vec<String> = Vec::new();
        let mut snapped: Vec<Option<usize>> = Vec::new();
        let mut unsnapped = 0_u64;
        for (fid, feature) in points.iter().enumerate() {
            let id = match id_idx {
                Some(i) => field_to_string(&feature.attributes[i]),
                None => fid.to_string(),
            };
            let node = feature
                .geometry
                .as_ref()
                .and_then(point_xy)
                .and_then(|(x, y)| graph.snap_point(x, y, tolerance));
            if node.is_none() {
                unsnapped += 1;
            }
            ids.push(id);
            snapped.push(node);
        }
        if snapped.iter().all(Option::is_none) {
            return Err(ToolError::Execution(
                "no input point snapped to the network within 'search_tolerance'".to_string(),
            ));
        }

        ctx.progress.info(&format!(
            "{} node(s), {} point(s), {unsnapped} unsnapped",
            graph.nodes.len(),
            ids.len()
        ));

        let mut out = Layer::new("generate_network_swm");
        out.add_field(FieldDef::new("origin", FieldType::Text));
        out.add_field(FieldDef::new("neighbor", FieldType::Text));
        out.add_field(FieldDef::new("weight", FieldType::Float));

        // Which node each destination point sits on, so one Dijkstra per
        // origin serves every destination at once.
        let mut by_node: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (i, n) in snapped.iter().enumerate() {
            if let Some(n) = n {
                by_node.entry(*n).or_default().push(i);
            }
        }

        let mut pairs = 0_u64;
        let total = ids.len().max(1);
        for (i, origin_node) in snapped.iter().enumerate() {
            let Some(origin_node) = origin_node else {
                continue;
            };
            let dist = graph.dijkstra(*origin_node, cutoff);

            // Collect reachable destinations, nearest first.
            // Scan the snapped destinations, not the reachability map: `dist`
            // holds every node within the cutoff (up to the whole network),
            // while `by_node` holds at most one entry per input point.
            let mut found: Vec<(usize, f64)> = Vec::new();
            for (node, members) in &by_node {
                let Some(d) = dist.get(node) else { continue };
                for &j in members {
                    if j != i {
                        found.push((j, *d));
                    }
                }
            }
            found.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
            if found.len() > max_neighbors {
                found.truncate(max_neighbors);
            }

            let mut weights: Vec<(usize, f64)> = found
                .into_iter()
                .map(|(j, d)| {
                    let w = if inverse {
                        // A coincident neighbour would divide by zero; treat
                        // zero impedance as maximal proximity instead.
                        if d <= 0.0 {
                            1.0
                        } else {
                            1.0 / d.powf(exponent)
                        }
                    } else {
                        1.0
                    };
                    (j, w)
                })
                .collect();

            if row_standardize {
                let sum: f64 = weights.iter().map(|(_, w)| *w).sum();
                if sum > 0.0 {
                    for (_, w) in weights.iter_mut() {
                        *w /= sum;
                    }
                }
            }

            for (j, w) in weights {
                out.add_feature(
                    None,
                    &[
                        ("origin", FieldValue::Text(ids[i].clone())),
                        ("neighbor", FieldValue::Text(ids[j].clone())),
                        ("weight", FieldValue::Float(w)),
                    ],
                )
                .map_err(|e| ToolError::Execution(e.to_string()))?;
                pairs += 1;
            }
            ctx.progress.progress((i as f64 + 1.0) / total as f64);
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("pair_count".to_string(), json!(pairs));
        outputs.insert("point_count".to_string(), json!(ids.len()));
        outputs.insert("unsnapped_count".to_string(), json!(unsnapped));
        outputs.insert("node_count".to_string(), json!(graph.nodes.len()));
        Ok(ToolRunResult { outputs })
    }
}

/// Undirected weighted graph built from a line layer.
///
/// Edges are kept alongside the adjacency so an input point can be snapped
/// onto the interior of a segment (splitting it) rather than jumping to the
/// nearest vertex — on a single 0..10 edge, points at x=2 and x=8 must be 6
/// apart, not 10.
struct Graph {
    nodes: Vec<(f64, f64)>,
    /// `(a, b, cost)` per undirected edge, indices into `nodes`.
    edges: Vec<(usize, usize, f64)>,
    adjacency: Vec<Vec<(usize, f64)>>,
    /// Weld grid, reused when a split introduces a node.
    snap: f64,
}

impl Graph {
    fn build(layer: &Layer, impedance_field: Option<&str>, snap: f64) -> Result<Graph, ToolError> {
        let imp_idx = match impedance_field {
            Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                ToolError::Validation(format!(
                    "impedance_field '{f}' not found on the network layer"
                ))
            })?),
            None => None,
        };

        let mut edges: Vec<(usize, usize, f64)> = Vec::new();
        let mut nodes: Vec<(f64, f64)> = Vec::new();
        let mut index: BTreeMap<(i64, i64), Vec<usize>> = BTreeMap::new();
        let mut adjacency: Vec<Vec<(usize, f64)>> = Vec::new();
        // Weld vertices onto a `snap`-sized grid so two segments meeting at a
        // shared endpoint connect. Neighbouring buckets are probed too:
        // quantisation alone still splits any pair straddling a bucket edge,
        // which would disconnect the network silently. A 1e-9 grid — the
        // previous behaviour — is effectively exact coordinate equality and
        // tolerates no float noise at all.
        let cell = snap.max(f64::MIN_POSITIVE);
        let mut node_of = |x: f64,
                           y: f64,
                           nodes: &mut Vec<(f64, f64)>,
                           adjacency: &mut Vec<Vec<(usize, f64)>>| {
            let kx = (x / cell).round() as i64;
            let ky = (y / cell).round() as i64;
            for dx in -1..=1 {
                for dy in -1..=1 {
                    // Every candidate in the bucket, not just the latest: two
                    // vertices farther apart than `snap` can still quantise
                    // together, and keeping only one would strand the other.
                    if let Some(bucket) = index.get(&(kx + dx, ky + dy)) {
                        for &existing in bucket {
                            let (nx, ny) = nodes[existing];
                            if ((nx - x).powi(2) + (ny - y).powi(2)).sqrt() <= snap {
                                return existing;
                            }
                        }
                    }
                }
            }
            nodes.push((x, y));
            adjacency.push(Vec::new());
            let id = nodes.len() - 1;
            index.entry((kx, ky)).or_default().push(id);
            id
        };

        for feature in layer.iter() {
            let parts: Vec<&Vec<Coord>> = match feature.geometry.as_ref() {
                Some(Geometry::LineString(cs)) => vec![cs],
                Some(Geometry::MultiLineString(ps)) => ps.iter().collect(),
                _ => continue,
            };
            // A feature-level impedance is apportioned across its segments by
            // length, so a long polyline is not charged the same as a short one.
            let declared = imp_idx
                .and_then(|i| feature.attributes.get(i))
                .and_then(field_to_f64);
            // Dijkstra is undefined for non-positive edges. On a cyclic network
            // — the normal case for streets and streams — relaxation keeps
            // succeeding around the cycle until the cost underflows, so the
            // tool would hang on user data rather than return an error.
            if let Some(d) = declared {
                if d <= 0.0 || !d.is_finite() {
                    return Err(ToolError::Validation(format!(
                        "impedance value {d} is not positive; network costs must be > 0"
                    )));
                }
            }
            // Apportion a FEATURE-level impedance across the feature's WHOLE
            // length, not each part's: otherwise a two-part feature declared as
            // 100 contributes 200 of traversal cost.
            let total: f64 = parts
                .iter()
                .map(|cs| cs.windows(2).map(seg_len).sum::<f64>())
                .sum();
            for cs in parts {
                for w in cs.windows(2) {
                    let len = seg_len(w);
                    if len <= 0.0 {
                        continue;
                    }
                    let cost = match declared {
                        Some(d) if total > 0.0 => d * (len / total),
                        _ => len,
                    };
                    let a = node_of(w[0].x, w[0].y, &mut nodes, &mut adjacency);
                    let b = node_of(w[1].x, w[1].y, &mut nodes, &mut adjacency);
                    if a != b {
                        edges.push((a, b, cost));
                    }
                }
            }
        }
        let mut g = Graph {
            nodes,
            edges,
            adjacency,
            snap,
        };
        g.rebuild_adjacency();
        Ok(g)
    }

    /// Recomputes adjacency from the current edge list.
    fn rebuild_adjacency(&mut self) {
        self.adjacency = vec![Vec::new(); self.nodes.len()];
        for &(a, b, cost) in &self.edges {
            self.adjacency[a].push((b, cost));
            self.adjacency[b].push((a, cost));
        }
    }

    /// Snaps a point onto the nearest segment within `tolerance`, splitting
    /// that segment at the projection when the point falls in its interior.
    ///
    /// Snapping to the nearest *vertex* instead would make two points on the
    /// same long edge measure the full edge length apart rather than their
    /// true along-edge separation, which distorts both cutoff inclusion and
    /// inverse-distance weights.
    fn snap_point(&mut self, x: f64, y: f64, tolerance: f64) -> Option<usize> {
        let mut best: Option<(usize, f64, f64)> = None; // (edge, t, distance)
        for (ei, &(a, b, _)) in self.edges.iter().enumerate() {
            let (ax, ay) = self.nodes[a];
            let (bx, by) = self.nodes[b];
            let (dx, dy) = (bx - ax, by - ay);
            let len2 = dx * dx + dy * dy;
            if len2 <= 0.0 {
                continue;
            }
            let t = (((x - ax) * dx + (y - ay) * dy) / len2).clamp(0.0, 1.0);
            let (px, py) = (ax + t * dx, ay + t * dy);
            let d = ((px - x).powi(2) + (py - y).powi(2)).sqrt();
            if d <= tolerance && best.is_none_or(|(_, _, bd)| d < bd) {
                best = Some((ei, t, d));
            }
        }
        let (ei, t, _) = best?;
        let (a, b, cost) = self.edges[ei];
        let (ax, ay) = self.nodes[a];
        let (bx, by) = self.nodes[b];
        let (px, py) = (ax + t * (bx - ax), ay + t * (by - ay));

        // Land on an existing endpoint when the projection is within the weld
        // tolerance of one, so a point at a junction does not split anything.
        if ((ax - px).powi(2) + (ay - py).powi(2)).sqrt() <= self.snap {
            return Some(a);
        }
        if ((bx - px).powi(2) + (by - py).powi(2)).sqrt() <= self.snap {
            return Some(b);
        }

        // Split: the edge's cost is apportioned along its length.
        self.nodes.push((px, py));
        let mid = self.nodes.len() - 1;
        self.edges[ei] = (a, mid, cost * t);
        self.edges.push((mid, b, cost * (1.0 - t)));
        self.rebuild_adjacency();
        Some(mid)
    }

    #[cfg(test)]
    fn nearest_node(&self, x: f64, y: f64, tolerance: f64) -> Option<usize> {
        let mut best: Option<(usize, f64)> = None;
        for (i, (nx, ny)) in self.nodes.iter().enumerate() {
            let d = ((nx - x).powi(2) + (ny - y).powi(2)).sqrt();
            if d <= tolerance && best.is_none_or(|(_, bd)| d < bd) {
                best = Some((i, d));
            }
        }
        best.map(|(i, _)| i)
    }

    /// Shortest network impedance from `start` to every node within `cutoff`.
    fn dijkstra(&self, start: usize, cutoff: f64) -> BTreeMap<usize, f64> {
        let mut dist: BTreeMap<usize, f64> = BTreeMap::new();
        let mut heap: BinaryHeap<Step> = BinaryHeap::new();
        dist.insert(start, 0.0);
        heap.push(Step {
            cost: 0.0,
            node: start,
        });
        while let Some(Step { cost, node }) = heap.pop() {
            if dist.get(&node).is_some_and(|d| cost > *d) {
                continue; // stale entry
            }
            for &(next, w) in &self.adjacency[node] {
                let nd = cost + w;
                if nd > cutoff {
                    continue;
                }
                if dist.get(&next).is_none_or(|d| nd < *d) {
                    dist.insert(next, nd);
                    heap.push(Step {
                        cost: nd,
                        node: next,
                    });
                }
            }
        }
        dist
    }
}

struct Step {
    cost: f64,
    node: usize,
}
impl PartialEq for Step {
    fn eq(&self, other: &Self) -> bool {
        self.cost == other.cost && self.node == other.node
    }
}
impl Eq for Step {}
impl Ord for Step {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .cost
            .total_cmp(&self.cost)
            .then_with(|| other.node.cmp(&self.node))
    }
}
impl PartialOrd for Step {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn seg_len(w: &[Coord]) -> f64 {
    ((w[1].x - w[0].x).powi(2) + (w[1].y - w[0].y).powi(2)).sqrt()
}

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) => cs.first().map(|c| (c.x, c.y)),
        _ => None,
    }
}

fn field_to_string(v: &FieldValue) -> String {
    match v {
        FieldValue::Text(s) => s.clone(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) => f.to_string(),
        FieldValue::Boolean(b) => b.to_string(),
        other => format!("{other:?}"),
    }
}

fn field_to_f64(v: &FieldValue) -> Option<f64> {
    match v {
        FieldValue::Float(f) => Some(*f),
        FieldValue::Integer(i) => Some(*i as f64),
        FieldValue::Text(s) => s.trim().parse().ok(),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn points(ps: Vec<(f64, f64)>) -> String {
        let mut l = Layer::new("pts");
        l.geom_type = Some(GeometryType::Point);
        for (x, y) in ps {
            l.add_feature(Some(Geometry::Point(Coord::xy(x, y))), &[])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn network(lines: Vec<Vec<(f64, f64)>>) -> String {
        let mut l = Layer::new("net");
        l.geom_type = Some(GeometryType::LineString);
        for cs in lines {
            l.add_feature(
                Some(Geometry::LineString(
                    cs.into_iter().map(|(x, y)| Coord::xy(x, y)).collect(),
                )),
                &[],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = GenerateNetworkSwmTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    /// All (origin, neighbor, weight) rows.
    fn rows(layer: &Layer) -> Vec<(String, String, f64)> {
        let o = layer.schema.field_index("origin").unwrap();
        let n = layer.schema.field_index("neighbor").unwrap();
        let w = layer.schema.field_index("weight").unwrap();
        layer
            .iter()
            .map(|f| {
                (
                    field_to_string(&f.attributes[o]),
                    field_to_string(&f.attributes[n]),
                    field_to_f64(&f.attributes[w]).unwrap(),
                )
            })
            .collect()
    }

    #[test]
    fn the_output_uses_the_interchangeable_weights_schema() {
        // The whole point: existing consumers must accept it unchanged.
        let (out, _) = run(json!({
            "input": points(vec![(0.0, 0.0), (10.0, 0.0)]),
            "network": network(vec![vec![(0.0, 0.0), (10.0, 0.0)]]),
        }));
        for f in ["origin", "neighbor", "weight"] {
            assert!(out.schema.field_index(f).is_some(), "missing field {f}");
        }
    }

    #[test]
    fn network_distance_beats_euclidean_when_they_disagree() {
        // Two points 2 units apart in a straight line but 22 units apart along
        // a U-shaped detour. A cutoff of 5 must find NO neighbours — which a
        // Euclidean matrix would never conclude.
        let pts = points(vec![(0.0, 0.0), (2.0, 0.0)]);
        let detour = network(vec![vec![(0.0, 0.0), (0.0, 10.0), (2.0, 10.0), (2.0, 0.0)]]);
        let (out, _) = run(json!({
            "input": pts, "network": detour, "distance_cutoff": 5.0,
        }));
        assert!(
            rows(&out).is_empty(),
            "euclidean-style neighbours leaked through: {:?}",
            rows(&out)
        );

        // Widening the cutoff past the true detour length finds them.
        let (out, _) = run(json!({
            "input": pts, "network": detour, "distance_cutoff": 25.0,
        }));
        assert_eq!(rows(&out).len(), 2);
    }

    #[test]
    fn a_disconnected_component_yields_no_pairs_across_it() {
        // Two separate line segments: points on one cannot reach the other.
        let (out, _) = run(json!({
            "input": points(vec![(0.0, 0.0), (100.0, 0.0)]),
            "network": network(vec![
                vec![(0.0, 0.0), (5.0, 0.0)],
                vec![(100.0, 0.0), (105.0, 0.0)],
            ]),
        }));
        assert!(rows(&out).is_empty());
    }

    #[test]
    fn fixed_conceptualization_gives_unit_weights() {
        let (out, _) = run(json!({
            "input": points(vec![(0.0, 0.0), (10.0, 0.0)]),
            "network": network(vec![vec![(0.0, 0.0), (10.0, 0.0)]]),
        }));
        assert!(rows(&out).iter().all(|(_, _, w)| (*w - 1.0).abs() < 1e-12));
    }

    #[test]
    fn inverse_conceptualization_decays_with_network_distance() {
        // Points at 0, 10 and 30 on a line. From origin 0, the nearer
        // neighbour must weigh more, and by the documented 1/d law.
        let (out, _) = run(json!({
            "input": points(vec![(0.0, 0.0), (10.0, 0.0), (30.0, 0.0)]),
            "network": network(vec![vec![(0.0, 0.0), (10.0, 0.0), (30.0, 0.0)]]),
            "conceptualization": "inverse",
        }));
        let from0: Vec<(String, f64)> = rows(&out)
            .into_iter()
            .filter(|(o, _, _)| o == "0")
            .map(|(_, n, w)| (n, w))
            .collect();
        let w1 = from0.iter().find(|(n, _)| n == "1").unwrap().1;
        let w2 = from0.iter().find(|(n, _)| n == "2").unwrap().1;
        assert!((w1 - 1.0 / 10.0).abs() < 1e-9, "got {w1}");
        assert!((w2 - 1.0 / 30.0).abs() < 1e-9, "got {w2}");
        assert!(w1 > w2);
    }

    #[test]
    fn row_standardization_makes_each_origin_total_one() {
        let (out, _) = run(json!({
            "input": points(vec![(0.0, 0.0), (10.0, 0.0), (30.0, 0.0)]),
            "network": network(vec![vec![(0.0, 0.0), (10.0, 0.0), (30.0, 0.0)]]),
            "conceptualization": "inverse",
            "row_standardization": true,
        }));
        let mut sums: BTreeMap<String, f64> = BTreeMap::new();
        for (o, _, w) in rows(&out) {
            *sums.entry(o).or_insert(0.0) += w;
        }
        for (o, s) in sums {
            assert!((s - 1.0).abs() < 1e-9, "origin {o} totals {s}");
        }
    }

    #[test]
    fn max_neighbors_keeps_the_nearest() {
        let (out, _) = run(json!({
            "input": points(vec![(0.0, 0.0), (10.0, 0.0), (30.0, 0.0)]),
            "network": network(vec![vec![(0.0, 0.0), (10.0, 0.0), (30.0, 0.0)]]),
            "max_neighbors": 1,
        }));
        let from0: Vec<String> = rows(&out)
            .into_iter()
            .filter(|(o, _, _)| o == "0")
            .map(|(_, n, _)| n)
            .collect();
        assert_eq!(from0, vec!["1".to_string()]);
    }

    #[test]
    fn an_impedance_field_overrides_geometric_length() {
        // A 10-unit segment declared as costing 100 must fall outside a
        // 50-unit cutoff, proving the field is honoured.
        let mut l = Layer::new("net");
        l.geom_type = Some(GeometryType::LineString);
        l.add_field(FieldDef::new("minutes", FieldType::Float));
        l.add_feature(
            Some(Geometry::LineString(vec![
                Coord::xy(0.0, 0.0),
                Coord::xy(10.0, 0.0),
            ])),
            &[("minutes", FieldValue::Float(100.0))],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        let net = memory_store::make_vector_memory_path(&id);

        let (out, _) = run(json!({
            "input": points(vec![(0.0, 0.0), (10.0, 0.0)]),
            "network": net,
            "impedance_field": "minutes",
            "distance_cutoff": 50.0,
        }));
        assert!(rows(&out).is_empty(), "impedance field ignored");
    }

    #[test]
    fn points_beyond_the_search_tolerance_are_reported_unsnapped() {
        let (_, res) = run(json!({
            "input": points(vec![(0.0, 0.0), (0.0, 900.0)]),
            "network": network(vec![vec![(0.0, 0.0), (10.0, 0.0)]]),
            "search_tolerance": 1.0,
        }));
        assert_eq!(res.outputs["unsnapped_count"], json!(1));
    }

    #[test]
    fn a_custom_id_field_labels_the_rows() {
        let mut l = Layer::new("pts");
        l.geom_type = Some(GeometryType::Point);
        l.add_field(FieldDef::new("name", FieldType::Text));
        for (x, n) in [(0.0, "west"), (10.0, "east")] {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(x, 0.0))),
                &[("name", FieldValue::Text(n.into()))],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        let pts = memory_store::make_vector_memory_path(&id);

        let (out, _) = run(json!({
            "input": pts,
            "network": network(vec![vec![(0.0, 0.0), (10.0, 0.0)]]),
            "id_field": "name",
        }));
        assert!(rows(&out)
            .iter()
            .any(|(o, n, _)| o == "west" && n == "east"));
    }

    #[test]
    fn a_multilinestring_impedance_is_shared_across_its_parts() {
        // Regression: `total` was computed per part, so each part received the
        // whole declared impedance and a two-part feature cost double.
        let mut l = Layer::new("net");
        l.geom_type = Some(GeometryType::MultiLineString);
        l.add_field(FieldDef::new("minutes", FieldType::Float));
        l.add_feature(
            Some(Geometry::MultiLineString(vec![
                vec![Coord::xy(0.0, 0.0), Coord::xy(5.0, 0.0)],
                vec![Coord::xy(5.0, 0.0), Coord::xy(10.0, 0.0)],
            ])),
            &[("minutes", FieldValue::Float(10.0))],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        let net = memory_store::make_vector_memory_path(&id);

        // The whole feature costs 10, so end-to-end is reachable at a cutoff
        // of 12. Charging each part 10 would make it 20 and find nothing.
        let (out, _) = run(json!({
            "input": points(vec![(0.0, 0.0), (10.0, 0.0)]),
            "network": net,
            "impedance_field": "minutes",
            "distance_cutoff": 12.0,
        }));
        assert_eq!(rows(&out).len(), 2, "feature impedance was not apportioned");
    }

    #[test]
    fn a_non_positive_impedance_is_rejected_rather_than_hanging() {
        // Dijkstra is undefined for negative edges; on a cyclic network the
        // relaxation would loop until the cost underflows.
        let mut l = Layer::new("net");
        l.geom_type = Some(GeometryType::LineString);
        l.add_field(FieldDef::new("minutes", FieldType::Float));
        l.add_feature(
            Some(Geometry::LineString(vec![
                Coord::xy(0.0, 0.0),
                Coord::xy(10.0, 0.0),
            ])),
            &[("minutes", FieldValue::Float(-5.0))],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        let net = memory_store::make_vector_memory_path(&id);

        let args: ToolArgs = serde_json::from_value(json!({
            "input": points(vec![(0.0, 0.0), (10.0, 0.0)]),
            "network": net,
            "impedance_field": "minutes",
        }))
        .unwrap();
        assert!(GenerateNetworkSwmTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn points_snap_onto_a_segment_not_to_its_nearest_vertex() {
        // On a single 0..10 edge, points at x=2 and x=8 are 6 apart along the
        // edge. Snapping each to the nearest VERTEX would report 10 and change
        // both cutoff inclusion and inverse-distance weights.
        let (out, _) = run(json!({
            "input": points(vec![(2.0, 0.0), (8.0, 0.0)]),
            "network": network(vec![vec![(0.0, 0.0), (10.0, 0.0)]]),
            "conceptualization": "inverse",
        }));
        let w = rows(&out)[0].2;
        assert!(
            (w - 1.0 / 6.0).abs() < 1e-9,
            "expected a 6-unit separation, got weight {w} (=1/{})",
            1.0 / w
        );
    }

    #[test]
    fn a_point_at_a_junction_does_not_split_the_edge() {
        let (out, _) = run(json!({
            "input": points(vec![(0.0, 0.0), (10.0, 0.0)]),
            "network": network(vec![vec![(0.0, 0.0), (10.0, 0.0)]]),
            "conceptualization": "inverse",
        }));
        let w = rows(&out)[0].2;
        assert!((w - 1.0 / 10.0).abs() < 1e-9, "got {w}");
    }

    #[test]
    fn rejects_bad_parameters() {
        let pts = points(vec![(0.0, 0.0)]);
        let net = network(vec![vec![(0.0, 0.0), (1.0, 0.0)]]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            GenerateNetworkSwmTool.validate(&args).is_err()
        };
        assert!(bad(json!({"network": net})));
        assert!(bad(
            json!({"input": pts, "network": net, "conceptualization": "nope"})
        ));
        assert!(bad(json!({"input": pts, "network": net, "exponent": 0})));
    }
}
