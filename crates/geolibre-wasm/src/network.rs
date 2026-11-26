//! Shortest-path searches over a routable line network, for GeoLibre's local
//! network tools (shortest path, closest facility, service area).
//!
//! GeoLibre builds the graph in JavaScript (topology, snapping, point
//! placement) and hands it over as compact typed arrays: node coordinates,
//! edge endpoints and per-direction costs, the source feature of each edge, and
//! the node adjacency in CSR form. This module then runs the searches:
//!
//! - **Node search**: Dijkstra, or A* toward a target with a haversine
//!   heuristic scaled by the cheapest cost per metre in the graph, which never
//!   overestimates, so A* returns the same optimum as Dijkstra.
//! - **Directed-edge search** once turn rules are set: a search state is "has
//!   just traversed edge `e` in direction `d`" (state `2e` runs `a -> b`, `2e+1`
//!   runs `b -> a`), so turn table entries, left/right turn costs and U-turn
//!   policies can depend on the edge a route arrived by.
//!
//! Both searches support a cost cap, reverse searches (following edges against
//! travel, giving costs *to* the sources), a stop condition over weighted
//! nodes (closest facility: stop once enough facilities are settled), and a
//! shared adjacency work budget. A [`RoutingSearch`] runs in chunks of bounded
//! work, so a caller can yield to its event loop and cancel between chunks.
//!
//! The search order mirrors GeoLibre's JavaScript engine exactly (the same
//! binary heap, tie-breaking and relaxation order), so the two engines return
//! the same costs and the same paths and can be checked against each other.
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::rc::Rc;

use wasm_bindgen::prelude::*;

const DEG: f64 = std::f64::consts::PI / 180.0;

const BUDGET_EXCEEDED: &str =
    "Local network adjacency work budget exceeded; use a smaller network or selection";

/// Where U-turns (reversing onto the edge just traversed) are allowed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum UTurnPolicy {
    /// Anywhere.
    Allow,
    /// Anywhere but the middle of a line (nodes of degree 2).
    Junctions,
    /// Only at dead ends (nodes of degree 0 or 1).
    DeadEnds,
    /// Nowhere.
    Never,
}

impl UTurnPolicy {
    /// Decode the policy code used across the JS boundary (0 allow,
    /// 1 junctions, 2 dead ends, 3 never).
    pub fn from_code(code: u8) -> Result<Self, String> {
        match code {
            0 => Ok(Self::Allow),
            1 => Ok(Self::Junctions),
            2 => Ok(Self::DeadEnds),
            3 => Ok(Self::Never),
            _ => Err(format!("invalid U-turn policy code {code}")),
        }
    }
}

/// Signed shortest longitude offset, including the dateline seam.
fn longitude_delta(from: f64, to: f64) -> f64 {
    if (to - from).abs() <= 180.0 {
        return to - from;
    }
    ((((to - from + 180.0) % 360.0) + 360.0) % 360.0) - 180.0
}

/// Great-circle distance between two lon/lat positions on a sphere of
/// `radius` metres.
fn haversine_meters(lon1: f64, lat1: f64, lon2: f64, lat2: f64, radius: f64) -> f64 {
    let s_lat = ((lat2 - lat1) * DEG / 2.0).sin();
    let s_lon = ((lon2 - lon1) * DEG / 2.0).sin();
    let h = s_lat * s_lat + (lat1 * DEG).cos() * (lat2 * DEG).cos() * s_lon * s_lon;
    2.0 * radius * h.sqrt().min(1.0).asin()
}

/// A routable network: nodes, edges with per-direction costs, and adjacency.
#[derive(Debug)]
pub struct Graph {
    lon: Vec<f64>,
    lat: Vec<f64>,
    edge_a: Vec<u32>,
    edge_b: Vec<u32>,
    cost_ab: Vec<f64>,
    cost_ba: Vec<f64>,
    feature: Vec<u32>,
    adj_offsets: Vec<u32>,
    adj_edges: Vec<u32>,
}

impl Graph {
    /// Build and validate a graph from flat arrays.
    ///
    /// `adj_offsets` has one entry per node plus one: the edges incident to
    /// node `n` are `adj_edges[adj_offsets[n]..adj_offsets[n + 1]]`, in the
    /// order the searches expand them. A cost of infinity closes that direction.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        lon: Vec<f64>,
        lat: Vec<f64>,
        edge_a: Vec<u32>,
        edge_b: Vec<u32>,
        cost_ab: Vec<f64>,
        cost_ba: Vec<f64>,
        feature: Vec<u32>,
        adj_offsets: Vec<u32>,
        adj_edges: Vec<u32>,
    ) -> Result<Self, String> {
        let n = lon.len();
        let m = edge_a.len();
        if lat.len() != n {
            return Err("lon and lat must have the same length".into());
        }
        if edge_b.len() != m || cost_ab.len() != m || cost_ba.len() != m || feature.len() != m {
            return Err("edge arrays must all have the same length".into());
        }
        // Node and state ids are stored as i32 in the predecessor arrays.
        if n > (i32::MAX as usize) || 2 * m > (i32::MAX as usize) {
            return Err("network is too large".into());
        }
        if edge_a.iter().chain(edge_b.iter()).any(|&v| v as usize >= n) {
            return Err("edge endpoint out of range".into());
        }
        // A self-loop has no direction for the directed-edge states to tell
        // apart, and never shortens a route; GeoLibre never builds one.
        if edge_a.iter().zip(&edge_b).any(|(a, b)| a == b) {
            return Err("self-loop edges are not supported".into());
        }
        if cost_ab.iter().chain(cost_ba.iter()).any(|c| c.is_nan() || *c < 0.0) {
            return Err("edge costs must be non-negative or infinity".into());
        }
        if adj_offsets.len() != n + 1 || adj_offsets[0] != 0 {
            return Err("adj_offsets must have node count + 1 entries starting at 0".into());
        }
        if adj_offsets.windows(2).any(|w| w[1] < w[0]) || adj_offsets[n] as usize != adj_edges.len() {
            return Err("adj_offsets must be non-decreasing and end at adj_edges length".into());
        }
        if adj_edges.iter().any(|&e| e as usize >= m) {
            return Err("adjacency edge out of range".into());
        }
        // Each adjacency entry must name an edge incident to its node, or the
        // searches would step from a node along an edge that does not touch it.
        for node in 0..n {
            let range = adj_offsets[node] as usize..adj_offsets[node + 1] as usize;
            for &e in &adj_edges[range] {
                let (a, b) = (edge_a[e as usize] as usize, edge_b[e as usize] as usize);
                if a != node && b != node {
                    return Err(format!("adjacency of node {node} lists edge {e}, which does not touch it"));
                }
            }
        }
        Ok(Self { lon, lat, edge_a, edge_b, cost_ab, cost_ba, feature, adj_offsets, adj_edges })
    }

    /// Number of nodes.
    pub fn node_count(&self) -> usize {
        self.lon.len()
    }

    /// Number of edges.
    pub fn edge_count(&self) -> usize {
        self.edge_a.len()
    }

    fn adj(&self, node: usize) -> &[u32] {
        &self.adj_edges[self.adj_offsets[node] as usize..self.adj_offsets[node + 1] as usize]
    }

    /// The node at the other end of `edge` from `node`, with the travel cost.
    fn step(&self, edge: usize, node: usize, reverse: bool) -> (usize, f64) {
        if self.edge_a[edge] as usize == node {
            (self.edge_b[edge] as usize, if reverse { self.cost_ba[edge] } else { self.cost_ab[edge] })
        } else {
            (self.edge_a[edge] as usize, if reverse { self.cost_ab[edge] } else { self.cost_ba[edge] })
        }
    }

    /// Start node of a directed-edge state.
    fn state_from(&self, state: usize) -> usize {
        let e = state >> 1;
        (if state & 1 == 1 { self.edge_b[e] } else { self.edge_a[e] }) as usize
    }

    /// End node of a directed-edge state.
    fn state_to(&self, state: usize) -> usize {
        let e = state >> 1;
        (if state & 1 == 1 { self.edge_a[e] } else { self.edge_b[e] }) as usize
    }

    /// Travel cost of a directed-edge state (infinity when not allowed).
    fn state_cost(&self, state: usize) -> f64 {
        let e = state >> 1;
        if state & 1 == 1 { self.cost_ba[e] } else { self.cost_ab[e] }
    }

    /// Signed deflection at `node` from arriving along `in_state` to leaving
    /// along `out_state`, in degrees: positive turns left, negative right.
    pub fn turn_angle(&self, in_state: usize, node: usize, out_state: usize) -> f64 {
        let from = self.state_from(in_state);
        let to = self.state_to(out_state);
        let k = (self.lat[node] * DEG).cos();
        let x1 = longitude_delta(self.lon[from], self.lon[node]) * k;
        let y1 = self.lat[node] - self.lat[from];
        let x2 = longitude_delta(self.lon[node], self.lon[to]) * k;
        let y2 = self.lat[to] - self.lat[node];
        (x1 * y2 - y1 * x2).atan2(x1 * x2 + y1 * y2) / DEG
    }
}

/// Hasher for the turn table's packed `(from, to)` feature keys: one
/// multiply instead of SipHash, since the keys are not attacker-chosen and the
/// lookup runs once per expanded turn.
#[derive(Default)]
struct PairHasher(u64);

impl Hasher for PairHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0.rotate_left(8) ^ u64::from(b)).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        }
    }

    fn write_u64(&mut self, value: u64) {
        self.0 = (value ^ (value >> 29)).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    }
}

type PairMap = HashMap<u64, f64, BuildHasherDefault<PairHasher>>;

/// Packs a `(from, to)` feature pair into one key.
fn pair_key(from: u32, to: u32) -> u64 {
    (u64::from(from) << 32) | u64::from(to)
}

/// Non-negative finite cost, or 0.
fn turn_cost_value(value: f64) -> f64 {
    if value.is_finite() && value > 0.0 { value } else { 0.0 }
}

/// Turn restrictions and costs, applied by the directed-edge search.
#[derive(Debug, Clone)]
pub struct TurnRules {
    table: PairMap,
    left_cost: f64,
    right_cost: f64,
    uturn_cost: f64,
    uturns: UTurnPolicy,
    straight_angle: f64,
}

impl TurnRules {
    /// Compile turn rules. Table entries from a feature to itself are ignored;
    /// a NaN or negative cost prohibits the turn (infinity), and when several
    /// entries name the same pair the largest cost wins, so a prohibition
    /// overrides a cost. Deflections up to `straight_angle` degrees either way
    /// are straight on (a non-finite angle means the default of 30).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        left_cost: f64,
        right_cost: f64,
        uturn_cost: f64,
        uturns: UTurnPolicy,
        straight_angle: f64,
        table_from: &[u32],
        table_to: &[u32],
        table_cost: &[f64],
    ) -> Result<Self, String> {
        if table_to.len() != table_from.len() || table_cost.len() != table_from.len() {
            return Err("turn table arrays must all have the same length".into());
        }
        let mut table = PairMap::with_capacity_and_hasher(table_from.len(), Default::default());
        for i in 0..table_from.len() {
            if table_from[i] == table_to[i] {
                continue;
            }
            let raw = table_cost[i];
            let cost = if raw.is_nan() || raw < 0.0 { f64::INFINITY } else { raw };
            let slot = table.entry(pair_key(table_from[i], table_to[i])).or_insert(0.0);
            if cost > *slot {
                *slot = cost;
            }
        }
        Ok(Self {
            table,
            left_cost: turn_cost_value(left_cost),
            right_cost: turn_cost_value(right_cost),
            uturn_cost: turn_cost_value(uturn_cost),
            uturns,
            straight_angle: if straight_angle.is_finite() { straight_angle.clamp(0.0, 180.0) } else { 30.0 },
        })
    }

    /// Cost of moving at `node` from traversal `in_state` onto `out_state`:
    /// infinity for a prohibited turn, else the turn table cost plus any left
    /// or right turn cost, or the U-turn cost when `out_state` reverses along
    /// the same edge.
    pub fn turn_cost(&self, graph: &Graph, in_state: usize, node: usize, out_state: usize) -> f64 {
        let in_edge = in_state >> 1;
        let out_edge = out_state >> 1;
        let degree = graph.adj(node).len();
        if in_edge == out_edge {
            let allowed = match self.uturns {
                UTurnPolicy::Allow => true,
                UTurnPolicy::Junctions => degree != 2,
                UTurnPolicy::DeadEnds => degree <= 1,
                UTurnPolicy::Never => false,
            };
            return if allowed { self.uturn_cost } else { f64::INFINITY };
        }
        let mut cost = 0.0;
        let from_feature = graph.feature[in_edge];
        let to_feature = graph.feature[out_edge];
        if from_feature != to_feature && !self.table.is_empty() {
            if let Some(&entry) = self.table.get(&pair_key(from_feature, to_feature)) {
                if !entry.is_finite() {
                    return f64::INFINITY;
                }
                cost += entry;
            }
        }
        // Only a junction offers a choice, and only leaving the line is a
        // turn: a bend in one road is free even where a side road joins there.
        if from_feature != to_feature && degree >= 3 && (self.left_cost > 0.0 || self.right_cost > 0.0) {
            let angle = graph.turn_angle(in_state, node, out_state);
            if angle > self.straight_angle {
                cost += self.left_cost;
            } else if angle < -self.straight_angle {
                cost += self.right_cost;
            }
        }
        cost
    }
}

/// Min binary heap of (priority, item) pairs with lazy deletion. The sift
/// order matches GeoLibre's JavaScript heap, so equal priorities pop in the
/// same order in both engines.
#[derive(Debug, Default)]
struct MinHeap {
    keys: Vec<f64>,
    vals: Vec<u32>,
}

impl MinHeap {
    fn push(&mut self, key: f64, val: u32) {
        let mut i = self.keys.len();
        self.keys.push(key);
        self.vals.push(val);
        while i > 0 {
            let p = (i - 1) >> 1;
            if self.keys[p] <= key {
                break;
            }
            self.keys[i] = self.keys[p];
            self.vals[i] = self.vals[p];
            i = p;
        }
        self.keys[i] = key;
        self.vals[i] = val;
    }

    fn pop(&mut self) -> Option<u32> {
        let top = *self.vals.first()?;
        let last_key = self.keys.pop()?;
        let last_val = self.vals.pop()?;
        let n = self.keys.len();
        if n > 0 {
            let mut i = 0;
            loop {
                let l = 2 * i + 1;
                if l >= n {
                    break;
                }
                let r = l + 1;
                let c = if r < n && self.keys[r] < self.keys[l] { r } else { l };
                if self.keys[c] >= last_key {
                    break;
                }
                self.keys[i] = self.keys[c];
                self.vals[i] = self.vals[c];
                i = c;
            }
            self.keys[i] = last_key;
            self.vals[i] = last_val;
        }
        Some(top)
    }
}

/// Inputs for one search.
#[derive(Debug, Clone)]
pub struct SearchOptions {
    /// Start nodes (cost 0).
    pub sources: Vec<u32>,
    /// Stop once this node is settled.
    pub target: Option<u32>,
    /// Use the A* heuristic toward `target`.
    pub astar: bool,
    /// Do not expand past this cost.
    pub max_cost: f64,
    /// Follow edges against their allowed direction.
    pub reverse: bool,
    /// Sphere radius in metres for the heuristic's haversine distance.
    pub radius_m: f64,
    /// Smallest finite cost per metre over all edges (the heuristic scale).
    pub min_cost_per_meter: f64,
    /// Per-node weights: settling a weighted node records it, and the search
    /// stops once the settled weight reaches `stop_count`.
    pub stop_weights: Option<Vec<u32>>,
    /// Settled weight that stops the search (with `stop_weights`).
    pub stop_count: u64,
    /// Adjacency entries this search may inspect before failing.
    pub budget: u64,
}

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
enum Phase {
    Start,
    Main,
    Done,
}

/// Directed-edge state arrays, present when the search applies turn rules.
#[derive(Debug)]
struct TurnState {
    rules: Rc<TurnRules>,
    entry: Vec<f64>,
    prev_state: Vec<i32>,
    best_state: Vec<i32>,
}

/// A resumable search: [`Search::run`] advances it by a bounded amount of
/// work, so a caller can interleave other work and cancel between calls.
#[derive(Debug)]
pub struct Search {
    graph: Rc<Graph>,
    opts: SearchOptions,
    phase: Phase,
    dist: Vec<f64>,
    prev_edge: Vec<i32>,
    prev_node: Vec<i32>,
    settled: Vec<u8>,
    heap: MinHeap,
    turns: Option<TurnState>,
    use_heuristic: bool,
    target_lon: f64,
    target_lat: f64,
    remaining: u64,
    visits: u64,
    found: u64,
    settled_targets: Vec<u32>,
}

impl Search {
    /// Prepare a search; no work is done until [`Search::run`].
    pub fn new(graph: Rc<Graph>, rules: Option<Rc<TurnRules>>, opts: SearchOptions) -> Result<Self, String> {
        let n = graph.node_count();
        if opts.sources.iter().any(|&s| s as usize >= n) {
            return Err("source node out of range".into());
        }
        if opts.target.is_some_and(|t| t as usize >= n) {
            return Err("target node out of range".into());
        }
        if opts.stop_weights.as_ref().is_some_and(|w| w.len() != n) {
            return Err("stop_weights must have one entry per node".into());
        }
        if opts.max_cost.is_nan() {
            return Err("max_cost must be a number".into());
        }
        if opts.stop_weights.is_some() && opts.stop_count == 0 {
            return Err("stop_count must be at least 1 with stop_weights".into());
        }
        // A weighted stop needs settle order by true cost, so it runs Dijkstra
        // even when a target is given.
        let use_heuristic = opts.astar
            && opts.target.is_some()
            && opts.stop_weights.is_none()
            && opts.min_cost_per_meter.is_finite()
            && opts.min_cost_per_meter > 0.0;
        let (target_lon, target_lat) = match opts.target {
            Some(t) => (graph.lon[t as usize], graph.lat[t as usize]),
            None => (0.0, 0.0),
        };
        let states = rules.as_ref().map_or(n, |_| 2 * graph.edge_count());
        let turns = rules.map(|rules| TurnState {
            rules,
            entry: vec![f64::INFINITY; states],
            prev_state: vec![-1; states],
            best_state: vec![-1; n],
        });
        Ok(Self {
            dist: vec![f64::INFINITY; n],
            prev_edge: vec![-1; n],
            prev_node: vec![-1; n],
            settled: vec![0; states],
            heap: MinHeap::default(),
            turns,
            use_heuristic,
            target_lon,
            target_lat,
            remaining: opts.budget,
            visits: 0,
            found: 0,
            settled_targets: Vec::new(),
            phase: Phase::Start,
            graph,
            opts,
        })
    }

    /// Heuristic cost from `node` to the target (0 without A*).
    fn h(&self, node: usize) -> f64 {
        if !self.use_heuristic {
            return 0.0;
        }
        let g = &self.graph;
        haversine_meters(g.lon[node], g.lat[node], self.target_lon, self.target_lat, self.opts.radius_m)
            * self.opts.min_cost_per_meter
    }

    /// Charge one adjacency entry against the budget, even one skipped later.
    fn charge(&mut self) -> Result<(), String> {
        if self.remaining == 0 {
            self.phase = Phase::Done;
            return Err(BUDGET_EXCEEDED.into());
        }
        self.remaining -= 1;
        self.visits += 1;
        Ok(())
    }

    /// Records a settled node against the stop condition; true to stop.
    fn settle(&mut self, node: usize) -> bool {
        let Some(weights) = &self.opts.stop_weights else { return false };
        let w = weights[node];
        if w > 0 {
            self.settled_targets.push(node as u32);
            self.found += u64::from(w);
        }
        self.found >= self.opts.stop_count
    }

    /// Advance the search by about `max_visits` adjacency entries (it always
    /// finishes the node or state it is expanding).
    ///
    /// Returns true once the search is complete. A search that exceeds its
    /// budget returns an error and is complete.
    pub fn run(&mut self, max_visits: u64) -> Result<bool, String> {
        let stop_at = self.visits.saturating_add(max_visits.max(1));
        if self.turns.is_some() {
            self.run_turns(stop_at)
        } else {
            self.run_nodes(stop_at)
        }
    }

    fn run_nodes(&mut self, stop_at: u64) -> Result<bool, String> {
        if self.phase == Phase::Start {
            for i in 0..self.opts.sources.len() {
                let s = self.opts.sources[i] as usize;
                if self.dist[s] == 0.0 {
                    continue;
                }
                self.dist[s] = 0.0;
                let key = self.h(s);
                self.heap.push(key, s as u32);
            }
            self.phase = Phase::Main;
        }
        let reverse = self.opts.reverse;
        while self.phase == Phase::Main {
            if self.visits >= stop_at {
                return Ok(false);
            }
            let Some(u) = self.heap.pop() else { break };
            let u = u as usize;
            if self.settled[u] != 0 {
                continue;
            }
            self.settled[u] = 1;
            if self.settle(u) || self.opts.target == Some(u as u32) {
                break;
            }
            let graph = Rc::clone(&self.graph);
            for &e in graph.adj(u) {
                self.charge()?;
                let (v, cost) = graph.step(e as usize, u, reverse);
                if !cost.is_finite() || self.settled[v] != 0 {
                    continue;
                }
                let nd = self.dist[u] + cost;
                if nd > self.opts.max_cost {
                    continue;
                }
                if nd < self.dist[v] {
                    self.dist[v] = nd;
                    self.prev_edge[v] = e as i32;
                    self.prev_node[v] = u as i32;
                    let key = nd + self.h(v);
                    self.heap.push(key, v as u32);
                }
            }
        }
        self.phase = Phase::Done;
        Ok(true)
    }

    /// The node a state reaches in search order.
    fn reached(&self, state: usize) -> usize {
        if self.opts.reverse { self.graph.state_from(state) } else { self.graph.state_to(state) }
    }

    /// The state leaving (or, reverse, entering) `node` along edge `e`.
    fn state_at(&self, e: usize, node: usize) -> usize {
        let g = &self.graph;
        let odd = if self.opts.reverse { g.edge_b[e] as usize != node } else { g.edge_a[e] as usize != node };
        2 * e + usize::from(odd)
    }

    fn run_turns(&mut self, stop_at: u64) -> Result<bool, String> {
        let graph = Rc::clone(&self.graph);
        let rules = Rc::clone(&self.turns.as_ref().expect("turn state").rules);
        let reverse = self.opts.reverse;
        let max_cost = self.opts.max_cost;
        if self.phase == Phase::Start {
            // Distinct sources in first-seen order (tie-breaking depends on it).
            let mut seen = vec![false; graph.node_count()];
            let mut starts: Vec<usize> = Vec::with_capacity(self.opts.sources.len());
            for &s in &self.opts.sources {
                if !std::mem::replace(&mut seen[s as usize], true) {
                    starts.push(s as usize);
                }
            }
            for &s in &starts {
                self.dist[s] = 0.0;
            }
            for &s in &starts {
                if self.settle(s) || self.opts.target == Some(s as u32) {
                    self.phase = Phase::Done;
                    return Ok(true);
                }
            }
            for &s in &starts {
                for &e in graph.adj(s) {
                    self.charge()?;
                    let state = self.state_at(e as usize, s);
                    let cost = graph.state_cost(state);
                    let turns = self.turns.as_mut().expect("turn state");
                    if !cost.is_finite() || turns.entry[state] <= 0.0 {
                        continue;
                    }
                    turns.entry[state] = 0.0;
                    if cost <= max_cost {
                        let key = cost + self.h(self.reached(state));
                        self.heap.push(key, state as u32);
                    }
                }
            }
            self.phase = Phase::Main;
        }
        while self.phase == Phase::Main {
            if self.visits >= stop_at {
                return Ok(false);
            }
            let Some(state) = self.heap.pop() else { break };
            let state = state as usize;
            if self.settled[state] != 0 {
                continue;
            }
            self.settled[state] = 1;
            let cost = self.turns.as_ref().expect("turn state").entry[state] + graph.state_cost(state);
            let node = self.reached(state);
            // The first state to reach a node is its cheapest: states pop in
            // cost order, and the heuristic is the same at one node.
            if self.dist[node] == f64::INFINITY {
                self.dist[node] = cost;
                self.turns.as_mut().expect("turn state").best_state[node] = state as i32;
                self.prev_edge[node] = (state >> 1) as i32;
                self.prev_node[node] =
                    (if reverse { graph.state_to(state) } else { graph.state_from(state) }) as i32;
                if self.settle(node) || self.opts.target == Some(node as u32) {
                    break;
                }
            }
            for &e in graph.adj(node) {
                self.charge()?;
                let next = self.state_at(e as usize, node);
                if self.settled[next] != 0 {
                    continue;
                }
                let next_cost = graph.state_cost(next);
                if !next_cost.is_finite() {
                    continue;
                }
                let turn = if reverse {
                    rules.turn_cost(&graph, next, node, state)
                } else {
                    rules.turn_cost(&graph, state, node, next)
                };
                if !turn.is_finite() {
                    continue;
                }
                let start = cost + turn;
                let turns = self.turns.as_mut().expect("turn state");
                if start < turns.entry[next] {
                    turns.entry[next] = start;
                    turns.prev_state[next] = state as i32;
                    if start + next_cost <= max_cost {
                        let key = start + next_cost + self.h(self.reached(next));
                        self.heap.push(key, next as u32);
                    }
                }
            }
        }
        self.phase = Phase::Done;
        Ok(true)
    }

    /// Whether the search has finished.
    pub fn done(&self) -> bool {
        self.phase == Phase::Done
    }

    /// Per node, the cheapest cost from (or, reverse, to) the sources.
    pub fn dist(&self) -> &[f64] {
        &self.dist
    }

    /// Per node, the edge on the best path (-1 for sources and unreached nodes).
    pub fn prev_edge(&self) -> &[i32] {
        &self.prev_edge
    }

    /// Per node, the previous node on the best path (-1 when none).
    pub fn prev_node(&self) -> &[i32] {
        &self.prev_node
    }

    /// Per directed-edge state, the cost of starting it (forward) or the cost
    /// still to go once it ends (reverse); empty without turn rules.
    pub fn entry(&self) -> &[f64] {
        self.turns.as_ref().map_or(&[], |t| &t.entry)
    }

    /// Per state, the neighboring state toward the sources (-1 when none).
    pub fn prev_state(&self) -> &[i32] {
        self.turns.as_ref().map_or(&[], |t| &t.prev_state)
    }

    /// Per node, the state that first reached it (-1 for sources).
    pub fn best_state(&self) -> &[i32] {
        self.turns.as_ref().map_or(&[], |t| &t.best_state)
    }

    /// Whether this search ran over directed-edge states.
    pub fn has_turns(&self) -> bool {
        self.turns.is_some()
    }

    /// Weighted nodes settled, in settle order (with `stop_weights`).
    pub fn settled_targets(&self) -> &[u32] {
        &self.settled_targets
    }

    /// Budget left after the work done so far.
    pub fn remaining_budget(&self) -> u64 {
        self.remaining
    }

    /// Adjacency entries inspected so far.
    pub fn visits(&self) -> u64 {
        self.visits
    }
}

// ───────────────────────── wasm-bindgen API ─────────────────────────

fn js_err(msg: String) -> JsValue {
    JsValue::from_str(&msg)
}

/// Clamp a JS number to a non-negative integer count (NaN counts as 0;
/// infinity and anything past `u64::MAX` saturate).
fn js_count(value: f64) -> u64 {
    if value.is_nan() || value <= 0.0 { 0 } else { value as u64 }
}

/// A routable network held in WebAssembly memory. Build it once per network,
/// optionally set turn rules, then start any number of [`RoutingSearch`]es.
#[wasm_bindgen]
pub struct RoutingGraph {
    graph: Rc<Graph>,
    rules: Option<Rc<TurnRules>>,
}

#[wasm_bindgen]
impl RoutingGraph {
    /// Build a graph from typed arrays: node `lon`/`lat`; per edge its end
    /// nodes, the cost from `a` to `b` and from `b` to `a` (`Infinity` closes a
    /// direction), and its source feature; and the node adjacency in CSR form
    /// (`adj_offsets` has node count + 1 entries).
    #[wasm_bindgen(constructor)]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        lon: Vec<f64>,
        lat: Vec<f64>,
        edge_a: Vec<u32>,
        edge_b: Vec<u32>,
        cost_ab: Vec<f64>,
        cost_ba: Vec<f64>,
        edge_feature: Vec<u32>,
        adj_offsets: Vec<u32>,
        adj_edges: Vec<u32>,
    ) -> Result<RoutingGraph, JsValue> {
        let graph = Graph::new(lon, lat, edge_a, edge_b, cost_ab, cost_ba, edge_feature, adj_offsets, adj_edges)
            .map_err(js_err)?;
        Ok(RoutingGraph { graph: Rc::new(graph), rules: None })
    }

    /// Number of nodes.
    #[wasm_bindgen(getter)]
    pub fn node_count(&self) -> u32 {
        self.graph.node_count() as u32
    }

    /// Number of edges.
    #[wasm_bindgen(getter)]
    pub fn edge_count(&self) -> u32 {
        self.graph.edge_count() as u32
    }

    /// Whether later searches apply turn rules.
    #[wasm_bindgen(getter)]
    pub fn has_turn_rules(&self) -> bool {
        self.rules.is_some()
    }

    /// Apply turn rules to later searches. `uturn_policy`: 0 allow, 1 at
    /// junctions, 2 at dead ends, 3 never. Table rows name features by index;
    /// a cost of `Infinity` (or NaN / negative) prohibits the turn.
    #[allow(clippy::too_many_arguments)]
    pub fn set_turn_rules(
        &mut self,
        left_cost: f64,
        right_cost: f64,
        uturn_cost: f64,
        uturn_policy: u8,
        straight_angle: f64,
        table_from: &[u32],
        table_to: &[u32],
        table_cost: &[f64],
    ) -> Result<(), JsValue> {
        let policy = UTurnPolicy::from_code(uturn_policy).map_err(js_err)?;
        let rules = TurnRules::new(
            left_cost, right_cost, uturn_cost, policy, straight_angle, table_from, table_to, table_cost,
        )
        .map_err(js_err)?;
        self.rules = Some(Rc::new(rules));
        Ok(())
    }

    /// Remove turn rules; later searches run over nodes.
    pub fn clear_turn_rules(&mut self) {
        self.rules = None;
    }

    /// Start a search. `target` is a node or -1 for none. `stop_weights`, when
    /// given, has one entry per node: the search stops once the settled weight
    /// reaches `stop_count`. `budget` caps the adjacency entries inspected.
    #[allow(clippy::too_many_arguments)]
    pub fn start_search(
        &self,
        sources: Vec<u32>,
        target: i32,
        astar: bool,
        max_cost: f64,
        reverse: bool,
        radius_m: f64,
        min_cost_per_meter: f64,
        stop_weights: Option<Vec<u32>>,
        stop_count: f64,
        budget: f64,
    ) -> Result<RoutingSearch, JsValue> {
        let opts = SearchOptions {
            sources,
            target: u32::try_from(target).ok(),
            astar,
            max_cost,
            reverse,
            radius_m,
            min_cost_per_meter,
            stop_weights,
            stop_count: js_count(stop_count),
            budget: js_count(budget),
        };
        let inner = Search::new(Rc::clone(&self.graph), self.rules.clone(), opts).map_err(js_err)?;
        Ok(RoutingSearch { inner })
    }
}

/// One search started by [`RoutingGraph::start_search`].
#[wasm_bindgen]
pub struct RoutingSearch {
    inner: Search,
}

#[wasm_bindgen]
impl RoutingSearch {
    /// Advance by about `max_visits` adjacency entries; true once complete.
    /// Throws when the budget runs out.
    pub fn run(&mut self, max_visits: f64) -> Result<bool, JsValue> {
        self.inner.run(js_count(max_visits)).map_err(js_err)
    }

    /// Whether the search has finished.
    #[wasm_bindgen(getter)]
    pub fn done(&self) -> bool {
        self.inner.done()
    }

    /// Whether the search ran over directed-edge states (turn rules).
    #[wasm_bindgen(getter)]
    pub fn has_turns(&self) -> bool {
        self.inner.has_turns()
    }

    /// Budget left after the work done so far.
    #[wasm_bindgen(getter)]
    pub fn remaining_budget(&self) -> f64 {
        self.inner.remaining_budget() as f64
    }

    /// Adjacency entries inspected so far.
    #[wasm_bindgen(getter)]
    pub fn visits(&self) -> f64 {
        self.inner.visits() as f64
    }

    /// Per node, the cheapest cost (`Infinity` when unreached).
    pub fn dist(&self) -> Vec<f64> {
        self.inner.dist().to_vec()
    }

    /// Per node, the edge used to reach it (-1 when none).
    pub fn prev_edge(&self) -> Vec<i32> {
        self.inner.prev_edge().to_vec()
    }

    /// Per node, the previous node on the best path (-1 when none).
    pub fn prev_node(&self) -> Vec<i32> {
        self.inner.prev_node().to_vec()
    }

    /// Per directed-edge state, its entry cost (empty without turn rules).
    pub fn entry(&self) -> Vec<f64> {
        self.inner.entry().to_vec()
    }

    /// Per state, the previous state on the best path (empty without turn rules).
    pub fn prev_state(&self) -> Vec<i32> {
        self.inner.prev_state().to_vec()
    }

    /// Per node, the state that first reached it (empty without turn rules).
    pub fn best_state(&self) -> Vec<i32> {
        self.inner.best_state().to_vec()
    }

    /// Weighted nodes settled, in settle order.
    pub fn settled_targets(&self) -> Vec<u32> {
        self.inner.settled_targets().to_vec()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const R: f64 = 6_371_008.8;

    /// Degrees spanning `m` metres along a meridian (or the equator).
    fn deg(m: f64) -> f64 {
        m / (R * DEG)
    }

    /// Builds a graph from metre coordinates and `(a, b, cost_ab, cost_ba,
    /// feature)` edges, with adjacency in edge order.
    fn graph(nodes: &[(f64, f64)], edges: &[(u32, u32, f64, f64, u32)]) -> Rc<Graph> {
        let n = nodes.len();
        let mut adj: Vec<Vec<u32>> = vec![Vec::new(); n];
        for (i, e) in edges.iter().enumerate() {
            adj[e.0 as usize].push(i as u32);
            adj[e.1 as usize].push(i as u32);
        }
        let mut offsets = vec![0u32];
        let mut flat = Vec::new();
        for list in &adj {
            flat.extend_from_slice(list);
            offsets.push(flat.len() as u32);
        }
        Rc::new(
            Graph::new(
                nodes.iter().map(|p| deg(p.0)).collect(),
                nodes.iter().map(|p| deg(p.1)).collect(),
                edges.iter().map(|e| e.0).collect(),
                edges.iter().map(|e| e.1).collect(),
                edges.iter().map(|e| e.2).collect(),
                edges.iter().map(|e| e.3).collect(),
                edges.iter().map(|e| e.4).collect(),
                offsets,
                flat,
            )
            .unwrap(),
        )
    }

    fn opts(sources: &[u32]) -> SearchOptions {
        SearchOptions {
            sources: sources.to_vec(),
            target: None,
            astar: false,
            max_cost: f64::INFINITY,
            reverse: false,
            radius_m: R,
            min_cost_per_meter: 1.0,
            stop_weights: None,
            stop_count: 0,
            budget: 50_000_000,
        }
    }

    fn finish(g: &Rc<Graph>, rules: Option<TurnRules>, o: SearchOptions) -> Search {
        let mut s = Search::new(Rc::clone(g), rules.map(Rc::new), o).unwrap();
        assert!(s.run(u64::MAX).unwrap());
        s
    }

    /// Nodes from the source to `target`, rebuilt like GeoLibre's pathTo.
    fn node_path(g: &Graph, s: &Search, target: usize) -> Vec<usize> {
        if s.has_turns() {
            let mut states = Vec::new();
            let mut st = s.best_state()[target];
            while st >= 0 {
                states.push(st as usize);
                st = s.prev_state()[st as usize];
            }
            states.reverse();
            let mut nodes = vec![g.state_from(states[0])];
            nodes.extend(states.iter().map(|&st| g.state_to(st)));
            return nodes;
        }
        let mut nodes = vec![target];
        let mut cur = target;
        while s.prev_node()[cur] >= 0 {
            cur = s.prev_node()[cur] as usize;
            nodes.push(cur);
        }
        nodes.reverse();
        nodes
    }

    /// A 3x3 grid, 100 m spacing, unit cost per metre, plus a costly diagonal.
    fn grid() -> Rc<Graph> {
        let mut nodes = Vec::new();
        for y in 0..3 {
            for x in 0..3 {
                nodes.push((x as f64 * 100.0, y as f64 * 100.0));
            }
        }
        let mut edges = Vec::new();
        let mut f = 0;
        for y in 0..3u32 {
            for x in 0..3u32 {
                let i = y * 3 + x;
                if x < 2 {
                    edges.push((i, i + 1, 100.0, 100.0, f));
                    f += 1;
                }
                if y < 2 {
                    edges.push((i, i + 3, 100.0, 100.0, f));
                    f += 1;
                }
            }
        }
        edges.push((0, 8, 1000.0, 1000.0, f));
        graph(&nodes, &edges)
    }

    #[test]
    fn dijkstra_and_astar_agree_on_cost() {
        let g = grid();
        let d = finish(&g, None, opts(&[0]));
        assert_eq!(d.dist()[8], 400.0);
        assert_eq!(d.dist()[4], 200.0);
        let mut o = opts(&[0]);
        o.target = Some(8);
        o.astar = true;
        o.min_cost_per_meter = 0.99;
        let a = finish(&g, None, o);
        assert_eq!(a.dist()[8], 400.0);
        assert_eq!(node_path(&g, &a, 8).len(), 5);
        // A* settles no more than Dijkstra did on the way to the target.
        assert!(a.visits() <= d.visits());
    }

    #[test]
    fn one_way_edges_and_reverse_search() {
        // 0 -> 1 only, 1 <-> 2.
        let g = graph(&[(0.0, 0.0), (100.0, 0.0), (200.0, 0.0)], &[
            (0, 1, 100.0, f64::INFINITY, 0),
            (1, 2, 100.0, 100.0, 1),
        ]);
        let fwd = finish(&g, None, opts(&[0]));
        assert_eq!(fwd.dist()[2], 200.0);
        let back = finish(&g, None, opts(&[2]));
        assert!(back.dist()[0].is_infinite());
        // Reverse from 2 gives costs *to* 2, so 0 reaches it.
        let mut o = opts(&[2]);
        o.reverse = true;
        let rev = finish(&g, None, o);
        assert_eq!(rev.dist()[0], 200.0);
    }

    #[test]
    fn max_cost_caps_the_search() {
        let g = grid();
        let mut o = opts(&[0]);
        o.max_cost = 150.0;
        let s = finish(&g, None, o);
        assert_eq!(s.dist()[1], 100.0);
        assert!(s.dist()[2].is_infinite());
    }

    #[test]
    fn stop_weights_stop_after_enough_targets() {
        let g = grid();
        let mut w = vec![0u32; 9];
        w[8] = 1;
        w[2] = 1;
        w[1] = 2;
        let mut o = opts(&[0]);
        o.stop_weights = Some(w);
        o.stop_count = 3;
        let s = finish(&g, None, o);
        // Node 1 (weight 2, cost 100) then node 2 (cost 200) reach 3.
        assert_eq!(s.settled_targets(), &[1, 2]);
        assert!(s.dist()[8].is_infinite() || s.dist()[8] >= 300.0);
    }

    #[test]
    fn chunked_run_matches_one_shot() {
        let g = grid();
        let whole = finish(&g, None, opts(&[0]));
        let mut s = Search::new(Rc::clone(&g), None, opts(&[0])).unwrap();
        let mut chunks = 0;
        while !s.run(1).unwrap() {
            chunks += 1;
        }
        assert!(chunks > 3);
        assert_eq!(s.dist(), whole.dist());
        assert_eq!(s.prev_edge(), whole.prev_edge());
    }

    #[test]
    fn budget_is_enforced_and_reported() {
        let g = grid();
        let mut o = opts(&[0]);
        o.budget = 3;
        let mut s = Search::new(Rc::clone(&g), None, o).unwrap();
        let err = s.run(u64::MAX).unwrap_err();
        assert_eq!(err, BUDGET_EXCEEDED);
        assert_eq!(s.remaining_budget(), 0);
        assert!(s.done());
    }

    /// A plus-shaped junction: 0 west, 1 centre, 2 east, 3 north, 4 south,
    /// each arm its own feature, plus a loop 3 - 5 - 2 around the north-east.
    fn plus() -> Rc<Graph> {
        graph(
            &[(0.0, 0.0), (100.0, 0.0), (200.0, 0.0), (100.0, 100.0), (100.0, -100.0), (200.0, 100.0)],
            &[
                (0, 1, 100.0, 100.0, 0),
                (1, 2, 100.0, 100.0, 1),
                (1, 3, 100.0, 100.0, 2),
                (1, 4, 100.0, 100.0, 3),
                (3, 5, 100.0, 100.0, 4),
                (5, 2, 100.0, 100.0, 5),
            ],
        )
    }

    #[test]
    fn turn_table_prohibition_forces_a_detour() {
        let g = plus();
        // No left turn from the west arm (feature 0) onto the north arm (2).
        // U-turns are barred so the detour around the block is the only route.
        let rules = TurnRules::new(0.0, 0.0, 0.0, UTurnPolicy::Never, f64::NAN, &[0], &[2], &[f64::INFINITY]).unwrap();
        let mut o = opts(&[0]);
        o.target = Some(3);
        let s = finish(&g, Some(rules), o);
        // West -> centre -> east -> NE -> north = 400 instead of 200.
        assert_eq!(s.dist()[3], 400.0);
        assert_eq!(node_path(&g, &s, 3), vec![0, 1, 2, 5, 3]);
    }

    #[test]
    fn left_and_right_turn_costs_use_the_deflection() {
        let g = plus();
        let in_state = 0; // edge 0, west -> centre
        // Edge 2 runs centre -> north: a left turn.
        assert!(g.turn_angle(in_state, 1, 4) > 80.0);
        // Edge 3 runs centre -> south: a right turn.
        assert!(g.turn_angle(in_state, 1, 6) < -80.0);
        let rules = TurnRules::new(7.0, 3.0, 0.0, UTurnPolicy::Allow, 30.0, &[], &[], &[]).unwrap();
        assert_eq!(rules.turn_cost(&g, in_state, 1, 4), 7.0);
        assert_eq!(rules.turn_cost(&g, in_state, 1, 6), 3.0);
        assert_eq!(rules.turn_cost(&g, in_state, 1, 2), 0.0);
        let mut o = opts(&[0]);
        o.target = Some(4);
        let s = finish(&g, Some(rules), o);
        assert_eq!(s.dist()[4], 203.0);
    }

    #[test]
    fn uturn_policy_at_dead_ends() {
        // A dead-end spur 1 - 2 off a line 0 - 1 - 3; reaching 0 from the spur
        // tip needs no U-turn, but a route that must reverse can only do so at
        // a dead end.
        let g = graph(&[(0.0, 0.0), (100.0, 0.0), (100.0, 100.0), (200.0, 0.0)], &[
            (0, 1, 100.0, 100.0, 0),
            (1, 2, 100.0, 100.0, 1),
            (1, 3, 100.0, 100.0, 2),
        ]);
        let never = TurnRules::new(0.0, 0.0, 0.0, UTurnPolicy::Never, 30.0, &[], &[], &[]).unwrap();
        assert!(never.turn_cost(&g, 2, 2, 3).is_infinite());
        let dead = TurnRules::new(0.0, 0.0, 5.0, UTurnPolicy::DeadEnds, 30.0, &[], &[], &[]).unwrap();
        // At node 2 (degree 1) a U-turn costs 5; at node 1 (degree 3) it is barred.
        assert_eq!(dead.turn_cost(&g, 2, 2, 3), 5.0);
        assert!(dead.turn_cost(&g, 0, 1, 1).is_infinite());
        let junction = TurnRules::new(0.0, 0.0, 5.0, UTurnPolicy::Junctions, 30.0, &[], &[], &[]).unwrap();
        assert_eq!(junction.turn_cost(&g, 0, 1, 1), 5.0);
    }

    #[test]
    fn turn_search_reverse_matches_forward_costs() {
        let g = plus();
        let rules = || {
            Some(TurnRules::new(7.0, 3.0, 50.0, UTurnPolicy::Allow, 30.0, &[0], &[2], &[f64::INFINITY]).unwrap())
        };
        let fwd = finish(&g, rules(), opts(&[0]));
        for t in 1..6u32 {
            let mut o = opts(&[t]);
            o.reverse = true;
            let rev = finish(&g, rules(), o);
            assert_eq!(rev.dist()[0], fwd.dist()[t as usize], "node {t}");
        }
    }

    #[test]
    fn turn_search_stops_at_a_source_target() {
        let g = plus();
        let rules = TurnRules::new(1.0, 1.0, 0.0, UTurnPolicy::Allow, 30.0, &[], &[], &[]).unwrap();
        let mut o = opts(&[2, 2]);
        o.target = Some(2);
        let s = finish(&g, Some(rules), o);
        assert_eq!(s.dist()[2], 0.0);
        assert_eq!(s.visits(), 0);
    }

    #[test]
    fn graph_validation_rejects_bad_input() {
        let ok = Graph::new(vec![0.0, 1.0], vec![0.0, 0.0], vec![0], vec![1], vec![1.0], vec![1.0], vec![0], vec![0, 1, 2], vec![0, 0]);
        assert!(ok.is_ok());
        let bad_end = Graph::new(vec![0.0, 1.0], vec![0.0, 0.0], vec![0], vec![2], vec![1.0], vec![1.0], vec![0], vec![0, 1, 2], vec![0, 0]);
        assert!(bad_end.is_err());
        let bad_cost = Graph::new(vec![0.0, 1.0], vec![0.0, 0.0], vec![0], vec![1], vec![-1.0], vec![1.0], vec![0], vec![0, 1, 2], vec![0, 0]);
        assert!(bad_cost.is_err());
        let bad_adj = Graph::new(vec![0.0, 1.0, 2.0], vec![0.0; 3], vec![0], vec![1], vec![1.0], vec![1.0], vec![0], vec![0, 1, 2, 3], vec![0, 0, 0]);
        assert!(bad_adj.is_err());
        let self_loop = Graph::new(vec![0.0, 1.0], vec![0.0, 0.0], vec![0, 0], vec![1, 0], vec![1.0, 1.0], vec![1.0, 1.0], vec![0, 1], vec![0, 3, 4], vec![0, 1, 1, 0]);
        assert!(self_loop.is_err());
        let g = grid();
        assert!(Search::new(Rc::clone(&g), None, opts(&[99])).is_err());
        let mut zero_stop = opts(&[0]);
        zero_stop.stop_weights = Some(vec![1; 9]);
        assert!(Search::new(Rc::clone(&g), None, zero_stop).is_err());
    }

    #[test]
    fn weighted_stop_ignores_the_heuristic() {
        // A line A(-100) - 0 - B(150) - T(300): facility A costs 100 and B 150
        // from node 0. A* toward T would pop B first (150 + 150 left) before
        // A (100 + 400 left); the stop must still see A, the nearer one.
        let g = graph(&[(0.0, 0.0), (-100.0, 0.0), (150.0, 0.0), (300.0, 0.0)], &[
            (0, 1, 100.0, 100.0, 0),
            (0, 2, 150.0, 150.0, 1),
            (2, 3, 150.0, 150.0, 2),
        ]);
        let mut o = opts(&[0]);
        o.stop_weights = Some(vec![0, 1, 1, 0]);
        o.stop_count = 1;
        o.target = Some(3);
        o.astar = true;
        o.min_cost_per_meter = 1.0;
        let s = finish(&g, None, o);
        assert_eq!(s.settled_targets(), &[1]);
    }

    #[test]
    fn longitude_delta_wraps_the_seam() {
        assert_eq!(longitude_delta(179.0, -179.0), 2.0);
        assert_eq!(longitude_delta(-179.0, 179.0), -2.0);
        assert_eq!(longitude_delta(10.0, 20.0), 10.0);
    }
}
