//! Fill-Spill-Merge (FSM) core numerics — clean-room Rust implementation.
//!
//! Distributes a finite quantity of surface water across a DEM to produce
//! realistic standing-water (lake / inundation) extents, rather than
//! unconditionally filling every depression to its spill point. Water flows
//! downhill into pits; each depression fills only as far as its available water
//! volume allows (partial fill); when a depression overflows its sill the excess
//! **spills** into the neighbouring depression, and adjacent flooded depressions
//! **merge** into larger lakes.
//!
//! This is a from-scratch reimplementation of the algorithm described in:
//!
//! > Barnes, R., Callaghan, K.L., Wickert, A.D. (2020). *Computing water flow
//! > through complex landscapes – Part 2: Finding hydrologic connectivity and
//! > dividing surface water*. Earth Surface Dynamics 8, 431–445.
//! > <https://doi.org/10.5194/esurf-8-431-2020>
//!
//! Written from the published paper and RichDEM's documentation — **no RichDEM
//! source code is copied** (RichDEM is GPL-3; this crate is MIT). The three
//! phases mirror the paper: build a depression hierarchy, move surface water
//! into pits, then move water through the hierarchy (spill/merge) and spread it
//! back out as standing water.

use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet, VecDeque};

/// Column offsets for the eight D8 neighbours, indexed 1..=8 (index 0 unused).
const DCOL: [isize; 9] = [0, -1, -1, 0, 1, 1, 1, 0, -1];
/// Row offsets for the eight D8 neighbours, indexed 1..=8.
const DROW: [isize; 9] = [0, 0, -1, -1, -1, 0, 1, 1, 1];
/// Direction leading from neighbour `n` back to the central cell.
const D8_INVERSE: [u8; 9] = [0, 5, 6, 7, 8, 1, 2, 3, 4];
/// Flow-direction sentinel for cells with no downhill neighbour (pits/ocean).
const NO_FLOW: u8 = 0;

/// Label of the ocean "depression" (the ultimate sink).
pub const OCEAN: u32 = 0;
/// Label for a cell not yet assigned to any depression.
const NO_DEP: u32 = u32::MAX;
/// Sentinel for an absent depression id (parent/child/overflow links).
const NO_VALUE: u32 = u32::MAX;
/// Sentinel for an absent cell index.
const NO_CELL: usize = usize::MAX;

/// Floating-point comparison tolerance (matches the reference algorithm).
const FP_ERR: f64 = 1e-6;

#[inline]
fn fp_eq(a: f64, b: f64) -> bool {
    (a - b).abs() < FP_ERR
}
#[inline]
fn fp_le(a: f64, b: f64) -> bool {
    a < b || (a - b).abs() < FP_ERR
}

/// One node of the depression hierarchy — a leaf depression or a metadepression
/// (the merge of two children). Volumes are computed for the whole subtree.
#[derive(Clone)]
struct Depression {
    /// Lowest cell in the depression (flat index), or [`NO_CELL`].
    pit_cell: usize,
    /// Outlet / sill cell through which the depression first overflows.
    out_cell: usize,
    /// Parent metadepression, or [`NO_VALUE`].
    parent: u32,
    /// Logical overflow depression (the sibling this one spills into).
    odep: u32,
    /// Geographic leaf into which overflow is initially routed.
    geolink: u32,
    /// Elevation of the pit cell.
    pit_elev: f64,
    /// Elevation of the outlet cell.
    out_elev: f64,
    /// Left / right children (both set, or both [`NO_VALUE`]).
    lchild: u32,
    rchild: u32,
    /// True when this depression spills directly to the ocean.
    ocean_parent: bool,
    /// Depressions that reach the ocean *through* this one (spill in, not up).
    ocean_linked: Vec<u32>,
    /// Cells contained in this depression and its children.
    cell_count: u32,
    /// Water-holding volume up to the outlet (Water-Level Equation).
    dep_vol: f64,
    /// Water currently held by this depression and its children.
    water_vol: f64,
    /// Sum of the elevations of the contained cells.
    total_elevation: f64,
}

impl Depression {
    fn new() -> Self {
        Depression {
            pit_cell: NO_CELL,
            out_cell: NO_CELL,
            parent: NO_VALUE,
            odep: NO_VALUE,
            geolink: NO_VALUE,
            pit_elev: f64::INFINITY,
            out_elev: f64::INFINITY,
            lchild: NO_VALUE,
            rchild: NO_VALUE,
            ocean_parent: false,
            ocean_linked: Vec::new(),
            cell_count: 0,
            dep_vol: 0.0,
            water_vol: 0.0,
            total_elevation: 0.0,
        }
    }
}

/// A priority-queue item ordered so the queue pops the lowest elevation first,
/// breaking ties in favour of the most recently inserted cell (needed to make a
/// single wavefront claim an entire flat).
struct PqItem {
    elev: f64,
    seq: u64,
    cell: usize,
}
impl PartialEq for PqItem {
    fn eq(&self, other: &Self) -> bool {
        self.elev == other.elev && self.seq == other.seq
    }
}
impl Eq for PqItem {}
impl PartialOrd for PqItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for PqItem {
    fn cmp(&self, other: &Self) -> Ordering {
        // `BinaryHeap` is a max-heap; the item that should pop first must be the
        // greatest. Lower elevation => greater; on a tie, higher seq => greater.
        match other.elev.total_cmp(&self.elev) {
            Ordering::Equal => self.seq.cmp(&other.seq),
            ord => ord,
        }
    }
}

/// Union-find over depression labels. `merge_a_into_b` makes `b`'s root the
/// parent of `a`'s root, so `find` returns the top-most metadepression.
struct Dsu {
    parent: Vec<u32>,
}
impl Dsu {
    fn with_len(n: usize) -> Self {
        Dsu {
            parent: (0..n as u32).collect(),
        }
    }
    fn ensure(&mut self, n: u32) {
        while (self.parent.len() as u32) <= n {
            let l = self.parent.len() as u32;
            self.parent.push(l);
        }
    }
    fn find(&mut self, x: u32) -> u32 {
        let mut x = x;
        while self.parent[x as usize] != x {
            let p = self.parent[x as usize];
            self.parent[x as usize] = self.parent[p as usize];
            x = p;
        }
        x
    }
    fn merge_a_into_b(&mut self, a: u32, b: u32) {
        let ra = self.find(a);
        let rb = self.find(b);
        self.parent[ra as usize] = rb;
    }
}

#[inline]
fn in_bounds(r: isize, c: isize, rows: usize, cols: usize) -> bool {
    r >= 0 && c >= 0 && (r as usize) < rows && (c as usize) < cols
}

/// Result of a Fill-Spill-Merge run.
pub struct FsmResult {
    /// Standing water depth per cell (`>= 0`; `0` where dry). Row-major.
    pub wtd: Vec<f64>,
    /// Total water volume that drained off the grid into the ocean sink.
    pub ocean_volume: f64,
    /// Number of leaf + meta depressions found (including the ocean).
    pub depression_count: usize,
}

/// Builds the ocean label array. A cell is ocean when it is NoData, on the grid
/// border (`edge_outlet`), or connected to the border through cells at or below
/// `ocean_level`. Everything else is [`NO_DEP`].
fn label_ocean(
    dem: &[f64],
    rows: usize,
    cols: usize,
    nodata: f64,
    ocean_level: Option<f64>,
    edge_outlet: bool,
) -> Vec<u32> {
    let n = rows * cols;
    let mut label = vec![NO_DEP; n];
    let mut queue: VecDeque<usize> = VecDeque::new();

    for r in 0..rows {
        for c in 0..cols {
            let i = r * cols + c;
            let is_border = r == 0 || c == 0 || r == rows - 1 || c == cols - 1;
            let is_nodata = dem[i] == nodata;
            if is_nodata || (edge_outlet && is_border) {
                // First (and only) labeling pass, so no cell is OCEAN yet.
                label[i] = OCEAN;
                queue.push_back(i);
            }
        }
    }

    // Bucket-fill inland through cells at or below the ocean level (coastal
    // inundation baseline). Without an ocean level this loop is a no-op.
    if let Some(level) = ocean_level {
        while let Some(i) = queue.pop_front() {
            let r = (i / cols) as isize;
            let c = (i % cols) as isize;
            for n in 1..=8 {
                let rn = r + DROW[n];
                let cn = c + DCOL[n];
                if !in_bounds(rn, cn, rows, cols) {
                    continue;
                }
                let ni = rn as usize * cols + cn as usize;
                if label[ni] == OCEAN {
                    continue;
                }
                if dem[ni] == nodata || dem[ni] <= level {
                    label[ni] = OCEAN;
                    queue.push_back(ni);
                }
            }
        }
    }

    label
}

/// Phase 1: build the depression hierarchy. Fills `label` (leaf depression of
/// each cell, `OCEAN` for freely-draining cells) and returns the depressions and
/// the D8 flow-direction grid. `label` must already mark ocean cells.
fn get_depression_hierarchy(
    dem: &[f64],
    rows: usize,
    cols: usize,
    label: &mut [u32],
) -> Result<(Vec<Depression>, Vec<u8>), String> {
    let n = rows * cols;
    let mut flowdirs = vec![NO_FLOW; n];
    let mut deps: Vec<Depression> = Vec::new();

    // Depression 0 is always the ocean.
    let mut ocean = Depression::new();
    ocean.pit_elev = f64::NEG_INFINITY;
    ocean.pit_cell = NO_CELL;
    deps.push(ocean);

    let mut pq: BinaryHeap<PqItem> = BinaryHeap::new();
    let mut seq: u64 = 0;
    let mut ocean_cells = 0usize;

    // Seed the queue with ocean cells that border a non-ocean cell.
    for r in 0..rows {
        for c in 0..cols {
            let i = r * cols + c;
            if label[i] != OCEAN {
                continue;
            }
            let mut borders_land = false;
            for k in 1..=8 {
                let rn = r as isize + DROW[k];
                let cn = c as isize + DCOL[k];
                if in_bounds(rn, cn, rows, cols) && label[rn as usize * cols + cn as usize] != OCEAN
                {
                    borders_land = true;
                    break;
                }
            }
            if borders_land {
                pq.push(PqItem {
                    elev: dem[i],
                    seq,
                    cell: i,
                });
                seq += 1;
                ocean_cells += 1;
            }
        }
    }
    if ocean_cells == 0 {
        return Err("no ocean/outlet cells found; set ocean_level or enable edge_outlet".into());
    }

    // Seed the queue with pit cells (no strictly-lower neighbour).
    for r in 0..rows {
        for c in 0..cols {
            let i = r * cols + c;
            if label[i] == OCEAN {
                continue;
            }
            let my_elev = dem[i];
            let mut has_lower = false;
            for k in 1..=8 {
                let rn = r as isize + DROW[k];
                let cn = c as isize + DCOL[k];
                if !in_bounds(rn, cn, rows, cols) {
                    continue;
                }
                if dem[rn as usize * cols + cn as usize] < my_elev {
                    has_lower = true;
                    break;
                }
            }
            if !has_lower {
                pq.push(PqItem {
                    elev: my_elev,
                    seq,
                    cell: i,
                });
                seq += 1;
            }
        }
    }

    // Outlets keyed by the (min,max) pair of the two depressions they connect;
    // we retain the lowest outlet between any pair.
    let mut outlet_db: HashMap<(u32, u32), (usize, f64)> = HashMap::new();

    while let Some(item) = pq.pop() {
        let ci = item.cell;
        let celev = item.elev;
        let cr = (ci / cols) as isize;
        let cc = (ci % cols) as isize;
        let mut clabel = label[ci];

        if clabel == NO_DEP {
            // The only NO_DEP cells popped are pit seeds not yet claimed by a
            // flat wavefront: each starts a new leaf depression.
            clabel = deps.len() as u32;
            let mut nd = Depression::new();
            nd.pit_cell = ci;
            nd.pit_elev = celev;
            deps.push(nd);
            label[ci] = clabel;
        }

        for k in 1..=8 {
            let rn = cr + DROW[k];
            let cn = cc + DCOL[k];
            if !in_bounds(rn, cn, rows, cols) {
                continue;
            }
            let ni = rn as usize * cols + cn as usize;
            let nlabel = label[ni];
            if nlabel == NO_DEP {
                label[ni] = clabel;
                flowdirs[ni] = D8_INVERSE[k];
                pq.push(PqItem {
                    elev: dem[ni],
                    seq,
                    cell: ni,
                });
                seq += 1;
            } else if nlabel == clabel {
                // Same depression — nothing to do.
            } else {
                // Found a link between two depressions. The outlet is the higher
                // of the two cells.
                let mut out_cell = ci;
                let mut out_elev = celev;
                if dem[ni] > out_elev {
                    out_cell = ni;
                    out_elev = dem[ni];
                }
                let key = if clabel < nlabel {
                    (clabel, nlabel)
                } else {
                    (nlabel, clabel)
                };
                outlet_db
                    .entry(key)
                    .and_modify(|e| {
                        if e.1 > out_elev {
                            *e = (out_cell, out_elev);
                        }
                    })
                    .or_insert((out_cell, out_elev));
            }
        }
    }

    // Build the hierarchy: visit outlets from lowest to highest, joining
    // depressions into metadepressions with a union-find.
    let mut outlets: Vec<(u32, u32, usize, f64)> = outlet_db
        .into_iter()
        .map(|((a, b), (cell, elev))| (a, b, cell, elev))
        .collect();
    outlets.sort_by(|x, y| x.3.total_cmp(&y.3));

    let mut djset = Dsu::with_len(deps.len());
    for outlet in &outlets {
        let (mut depa, mut depb, out_cell, out_elev) = *outlet;
        let mut depa_set = djset.find(depa);
        let mut depb_set = djset.find(depb);
        if depa_set == depb_set {
            continue;
        }

        if depa_set == OCEAN || depb_set == OCEAN {
            // Exactly one links to the ocean; make depb the ocean side.
            if depa_set == OCEAN {
                std::mem::swap(&mut depa, &mut depb);
                std::mem::swap(&mut depa_set, &mut depb_set);
            }
            let dep = &mut deps[depa_set as usize];
            dep.parent = depb;
            dep.out_elev = out_elev;
            dep.out_cell = out_cell;
            dep.odep = NO_VALUE;
            dep.ocean_parent = true;
            dep.geolink = depb;
            deps[depb as usize].ocean_linked.push(depa_set);
            djset.merge_a_into_b(depa_set, OCEAN);
        } else {
            // Neither has found the ocean: merge into a new metadepression.
            let newlabel = deps.len() as u32;
            let depa_pit = deps[depa_set as usize].pit_cell;
            {
                let da = &mut deps[depa_set as usize];
                da.parent = newlabel;
                da.out_cell = out_cell;
                da.out_elev = out_elev;
                da.odep = depb_set;
                da.geolink = depb;
            }
            {
                let db = &mut deps[depb_set as usize];
                db.parent = newlabel;
                db.out_cell = out_cell;
                db.out_elev = out_elev;
                db.odep = depa_set;
                db.geolink = depa;
            }
            let mut nd = Depression::new();
            nd.lchild = depa_set;
            nd.rchild = depb_set;
            nd.pit_cell = depa_pit;
            deps.push(nd);
            djset.ensure(newlabel);
            djset.merge_a_into_b(depa_set, newlabel);
            djset.merge_a_into_b(depb_set, newlabel);
        }
    }

    calculate_marginal_volumes(&mut deps, dem, label);
    calculate_total_volumes(&mut deps);

    Ok((deps, flowdirs))
}

/// Assigns each cell's elevation to the deepest metadepression whose outlet lies
/// above the cell (its "marginal" contribution). Cells that drain to the ocean
/// contribute nothing.
fn calculate_marginal_volumes(deps: &mut [Depression], dem: &[f64], label: &[u32]) {
    let mut counts = vec![0u32; deps.len()];
    let mut elevs = vec![0.0f64; deps.len()];
    for i in 0..dem.len() {
        let my_elev = dem[i];
        let mut clabel = label[i];
        while clabel != OCEAN {
            let out_elev = deps[clabel as usize].out_elev;
            if my_elev <= out_elev {
                break;
            }
            if deps[clabel as usize].ocean_parent {
                clabel = OCEAN;
                break;
            }
            clabel = deps[clabel as usize].parent;
        }
        if clabel == OCEAN {
            continue;
        }
        counts[clabel as usize] += 1;
        elevs[clabel as usize] += my_elev;
    }
    for (d, dep) in deps.iter_mut().enumerate() {
        dep.cell_count += counts[d];
        dep.total_elevation += elevs[d];
    }
}

/// Rolls child cell counts and elevations up into their parents, then computes
/// each depression's water-holding volume via the Water-Level Equation.
fn calculate_total_volumes(deps: &mut [Depression]) {
    for d in 0..deps.len() {
        if deps[d].lchild != NO_VALUE {
            let (lc, rc) = (deps[d].lchild as usize, deps[d].rchild as usize);
            let (lcount, lelev) = (deps[lc].cell_count, deps[lc].total_elevation);
            let (rcount, relev) = (deps[rc].cell_count, deps[rc].total_elevation);
            deps[d].cell_count += lcount + rcount;
            deps[d].total_elevation += lelev + relev;
        }
        if d as u32 == OCEAN {
            continue;
        }
        deps[d].dep_vol = deps[d].cell_count as f64 * deps[d].out_elev - deps[d].total_elevation;
    }
}

/// Phase 2: route each cell's surface water down the flow field into its leaf
/// depression's pit, accumulating into `deps[*].water_vol`. Afterwards every
/// `wtd` value is `0` (all surface water has moved into the hierarchy).
fn move_water_into_pits(
    rows: usize,
    cols: usize,
    label: &[u32],
    flowdirs: &[u8],
    deps: &mut [Depression],
    wtd: &mut [f64],
) {
    let n = rows * cols;
    let mut dependencies = vec![0u32; n];
    for r in 0..rows {
        for c in 0..cols {
            let i = r * cols + c;
            for k in 1..=8 {
                let rn = r as isize + DROW[k];
                let cn = c as isize + DCOL[k];
                if !in_bounds(rn, cn, rows, cols) {
                    continue;
                }
                let ni = rn as usize * cols + cn as usize;
                if flowdirs[ni] == D8_INVERSE[k] {
                    dependencies[i] += 1;
                }
            }
        }
    }

    let mut q: VecDeque<usize> = VecDeque::new();
    for (i, &d) in dependencies.iter().enumerate() {
        if d == 0 {
            q.push_back(i);
        }
    }

    while let Some(c) = q.pop_front() {
        let dir = flowdirs[c];
        if dir == NO_FLOW {
            // Pit (or ocean) cell: deposit its water into the depression.
            if wtd[c] > 0.0 {
                deps[label[c] as usize].water_vol += wtd[c];
                wtd[c] = 0.0;
            }
        } else {
            let cr = (c / cols) as isize;
            let cc = (c % cols) as isize;
            let nr = cr + DROW[dir as usize];
            let nc = cc + DCOL[dir as usize];
            let ncell = nr as usize * cols + nc as usize;
            if wtd[c] > 0.0 {
                wtd[ncell] += wtd[c];
                wtd[c] = 0.0;
            }
            dependencies[ncell] -= 1;
            if dependencies[ncell] == 0 {
                q.push_back(ncell);
            }
        }
    }
}

/// Returns the traversal children of a depression: its ocean-linked depressions
/// followed by its two structural children (skipping absent ones).
fn children_of(dep: &Depression) -> Vec<u32> {
    let mut kids = dep.ocean_linked.clone();
    if dep.lchild != NO_VALUE {
        kids.push(dep.lchild);
    }
    if dep.rchild != NO_VALUE {
        kids.push(dep.rchild);
    }
    kids
}

/// Phase 3a: walk the hierarchy bottom-up; when a depression holds more water
/// than it can contain, overflow the excess into its neighbour and, if needed,
/// its parent (spill/merge). Implemented as an explicit post-order traversal to
/// avoid deep recursion on large grids.
fn move_water_in_dep_hier(deps: &mut [Depression]) {
    let mut jump_table: HashMap<u32, u32> = HashMap::new();
    let mut stack: Vec<(u32, usize)> = vec![(OCEAN, 0)];
    while let Some(&(node, ci)) = stack.last() {
        let kids = children_of(&deps[node as usize]);
        if ci < kids.len() {
            stack.last_mut().unwrap().1 += 1;
            stack.push((kids[ci], 0));
            continue;
        }
        stack.pop();

        if node == OCEAN {
            continue;
        }

        let lchild = deps[node as usize].lchild;
        let rchild = deps[node as usize].rchild;
        if lchild != NO_VALUE {
            let lfull = deps[lchild as usize].water_vol == deps[lchild as usize].dep_vol;
            let rfull = deps[rchild as usize].water_vol == deps[rchild as usize].dep_vol;
            if lfull && rfull && deps[node as usize].water_vol == 0.0 {
                deps[node as usize].water_vol +=
                    deps[lchild as usize].water_vol + deps[rchild as usize].water_vol;
            }
        }

        if deps[node as usize].water_vol > deps[node as usize].dep_vol {
            let parent = deps[node as usize].parent;
            overflow_into(node, parent, deps, &mut jump_table, 0.0);
        }
    }
}

/// Moves `extra_water` starting at depression `root`, stashing it in `root`, its
/// overflow neighbour, or its parent, chaining until it reaches `stop_node` or
/// the ocean. Iterative rewrite of the reference's recursive `OverflowInto`; the
/// `jump_table` maps every visited depression straight to the final destination
/// so the whole traversal stays near-linear.
fn overflow_into(
    root: u32,
    stop_node: u32,
    deps: &mut [Depression],
    jump_table: &mut HashMap<u32, u32>,
    extra_water: f64,
) -> u32 {
    let mut root = root;
    let mut extra_water = extra_water;
    let mut chain: Vec<u32> = Vec::new();

    let dest = loop {
        // Absorb any of this depression's own excess.
        if deps[root as usize].water_vol > deps[root as usize].dep_vol {
            extra_water += deps[root as usize].water_vol - deps[root as usize].dep_vol;
            deps[root as usize].water_vol = deps[root as usize].dep_vol;
        }

        if root == stop_node || root == OCEAN {
            deps[root as usize].water_vol += extra_water;
            break root;
        }

        // (1) Stash water in this depression.
        if deps[root as usize].water_vol < deps[root as usize].dep_vol {
            let capacity = deps[root as usize].dep_vol - deps[root as usize].water_vol;
            if extra_water < capacity {
                deps[root as usize].water_vol =
                    (deps[root as usize].water_vol + extra_water).min(deps[root as usize].dep_vol);
                extra_water = 0.0;
            } else {
                deps[root as usize].water_vol = deps[root as usize].dep_vol;
                extra_water -= capacity;
            }
        }
        if fp_eq(extra_water, 0.0) {
            break root;
        }

        // Jump-table shortcut past already-filled depressions.
        if let Some(&j) = jump_table.get(&root) {
            chain.push(root);
            root = j;
            continue;
        }

        // (2) Stash water in the overflow neighbour.
        let odep = deps[root as usize].odep;
        let geolink = deps[root as usize].geolink;
        if odep != NO_VALUE {
            let ow = deps[odep as usize].water_vol;
            let ov = deps[odep as usize].dep_vol;
            if ow < ov {
                chain.push(root);
                root = geolink;
                continue;
            } else if ow > ov {
                extra_water += ow - ov;
                deps[odep as usize].water_vol = ov;
            }
        }

        // (3) Pass water up to the parent.
        let parent = deps[root as usize].parent;
        let ocean_parent = deps[root as usize].ocean_parent;
        if deps[parent as usize].water_vol == 0.0 && !ocean_parent {
            let add = deps[root as usize].water_vol
                + if odep != NO_VALUE {
                    deps[odep as usize].water_vol
                } else {
                    0.0
                };
            deps[parent as usize].water_vol += add;
        }
        chain.push(root);
        root = parent;
    };

    for r in chain {
        jump_table.insert(r, dest);
    }
    dest
}

/// Info bubbled up while deciding which depressions to spread water across.
#[derive(Clone)]
struct SubtreeInfo {
    leaf_label: u32,
    my_labels: HashSet<u32>,
}
impl SubtreeInfo {
    fn empty() -> Self {
        SubtreeInfo {
            leaf_label: NO_VALUE,
            my_labels: HashSet::new(),
        }
    }
}

/// Phase 3b: find each partially-filled (meta)depression and spread its water as
/// standing surface water. Explicit post-order traversal; each node reads its
/// children's subtree info from `info`.
fn find_depressions_to_fill(
    dem: &[f64],
    cols: usize,
    label: &[u32],
    deps: &[Depression],
    wtd: &mut [f64],
) {
    let mut info: HashMap<u32, SubtreeInfo> = HashMap::new();
    let mut stack: Vec<(u32, usize)> = vec![(OCEAN, 0)];
    while let Some(&(node, ci)) = stack.last() {
        let kids = children_of(&deps[node as usize]);
        if ci < kids.len() {
            stack.last_mut().unwrap().1 += 1;
            stack.push((kids[ci], 0));
            continue;
        }
        stack.pop();

        if node == OCEAN {
            info.insert(node, SubtreeInfo::empty());
            continue;
        }

        let this = &deps[node as usize];
        let left = if this.lchild != NO_VALUE {
            info.get(&this.lchild)
                .cloned()
                .unwrap_or_else(SubtreeInfo::empty)
        } else {
            SubtreeInfo::empty()
        };
        let right = if this.rchild != NO_VALUE {
            info.get(&this.rchild)
                .cloned()
                .unwrap_or_else(SubtreeInfo::empty)
        } else {
            SubtreeInfo::empty()
        };

        let mut combined = SubtreeInfo::empty();
        combined.my_labels.insert(node);
        for l in left.my_labels {
            combined.my_labels.insert(l);
        }
        for l in right.my_labels {
            combined.my_labels.insert(l);
        }
        combined.leaf_label = if left.leaf_label != NO_VALUE {
            left.leaf_label
        } else {
            node
        };

        let parent_water = if this.parent != NO_VALUE {
            deps[this.parent as usize].water_vol
        } else {
            0.0
        };
        let can_fill_here = this.water_vol < this.dep_vol
            || this.ocean_parent
            || (this.water_vol == this.dep_vol && parent_water == 0.0);

        if can_fill_here {
            let pit_cell = deps[combined.leaf_label as usize].pit_cell;
            fill_depressions(
                pit_cell,
                this.out_cell,
                &combined.my_labels,
                this.water_vol,
                dem,
                cols,
                label,
                wtd,
            );
            info.insert(node, SubtreeInfo::empty());
        } else {
            info.insert(node, combined);
        }
    }
}

/// Volume dammed behind a sill of the given elevation: `count*sill - Σelev`.
#[inline]
fn depression_volume(sill_elev: f64, cells: usize, total_elev: f64) -> f64 {
    cells as f64 * sill_elev - total_elev
}

/// The Lake-Level Equation: elevation the water surface settles at given the
/// water volume, the sill elevation, and the cells spread across.
fn determine_water_level(water_vol: f64, sill_elev: f64, cells: usize, total_elev: f64) -> f64 {
    let current = depression_volume(sill_elev, cells, total_elev);
    if water_vol > current {
        // Excess beyond the above-ground volume; the surface is at the sill.
        sill_elev
    } else if fp_eq(water_vol, current) {
        sill_elev
    } else {
        let nominal = (water_vol + total_elev) / cells as f64;
        if fp_eq(nominal, sill_elev) {
            sill_elev
        } else {
            nominal
        }
    }
}

/// Raises the water table of every affected cell to `water_level`.
fn backfill_depression(water_level: f64, dem: &[f64], wtd: &mut [f64], cells: &[usize]) {
    for &c in cells {
        wtd[c] = (water_level - dem[c]).max(0.0);
    }
}

/// Spreads `water_vol` across a (meta)depression by priority-flooding from its
/// pit up to the level the volume supports, writing standing water into `wtd`.
#[allow(clippy::too_many_arguments)]
fn fill_depressions(
    pit_cell: usize,
    out_cell: usize,
    dep_labels: &HashSet<u32>,
    water_vol: f64,
    dem: &[f64],
    cols: usize,
    label: &[u32],
    wtd: &mut [f64],
) {
    if water_vol == 0.0 || pit_cell == NO_CELL {
        return;
    }
    let out_elev = dem[out_cell];

    let mut visited: HashSet<usize> = HashSet::new();
    let mut flood: BinaryHeap<PqItem> = BinaryHeap::new();
    let mut seq: u64 = 0;

    flood.push(PqItem {
        elev: dem[pit_cell],
        seq,
        cell: pit_cell,
    });
    seq += 1;
    visited.insert(pit_cell);

    let mut water_vol = water_vol;
    let mut cells_affected: Vec<usize> = Vec::new();
    let mut total_elevation = 0.0;

    while let Some(item) = flood.pop() {
        let c = item.cell;
        let cz = dem[c];
        let current_volume = depression_volume(cz, cells_affected.len(), total_elevation);

        if fp_le(water_vol, current_volume - wtd[c]) {
            let mut water_level =
                determine_water_level(water_vol, cz, cells_affected.len(), total_elevation);
            if fp_eq(water_level, out_elev) {
                water_level = out_elev;
            }
            backfill_depression(water_level, dem, wtd, &cells_affected);
            return;
        }

        if c != out_cell {
            cells_affected.push(c);
            water_vol += wtd[c]; // wtd <= 0 here (surface water already moved out)
            wtd[c] = 0.0;
            total_elevation += cz;
        }

        let cr = (c / cols) as isize;
        let cc = (c % cols) as isize;
        let rows = dem.len() / cols;
        for k in 1..=8 {
            let rn = cr + DROW[k];
            let cn = cc + DCOL[k];
            if !in_bounds(rn, cn, rows, cols) {
                continue;
            }
            let ni = rn as usize * cols + cn as usize;
            if !dep_labels.contains(&label[ni]) && ni != out_cell {
                continue;
            }
            if dem[ni] > out_elev {
                continue;
            }
            if !visited.contains(&ni) {
                flood.push(PqItem {
                    elev: dem[ni],
                    seq,
                    cell: ni,
                });
                seq += 1;
                visited.insert(ni);
            }
        }

        if flood.is_empty() {
            flood.push(PqItem {
                elev: dem[out_cell],
                seq,
                cell: out_cell,
            });
            seq += 1;
            visited.insert(out_cell);
        }
    }
}

/// Runs Fill-Spill-Merge on a row-major DEM.
///
/// `initial_wtd` is the surface water depth applied to each cell (uniform or
/// per-cell, all `>= 0`). Returns the standing water depth per cell plus summary
/// figures. NoData cells and grid-edge/ocean cells act as free outlets.
#[allow(clippy::too_many_arguments)]
pub fn fill_spill_merge(
    dem: &[f64],
    rows: usize,
    cols: usize,
    nodata: f64,
    initial_wtd: &[f64],
    ocean_level: Option<f64>,
    edge_outlet: bool,
) -> Result<FsmResult, String> {
    let n = rows * cols;
    if dem.len() != n || initial_wtd.len() != n {
        return Err("DEM and water buffers must match rows*cols".into());
    }

    let mut label = label_ocean(dem, rows, cols, nodata, ocean_level, edge_outlet);

    // Start water table: the requested surface water on land, zero on ocean and
    // NoData cells.
    let mut wtd = vec![0.0f64; n];
    for i in 0..n {
        if label[i] != OCEAN && dem[i] != nodata {
            wtd[i] = initial_wtd[i].max(0.0);
        }
    }

    let (mut deps, flowdirs) = get_depression_hierarchy(dem, rows, cols, &mut label)?;
    move_water_into_pits(rows, cols, &label, &flowdirs, &mut deps, &mut wtd);
    move_water_in_dep_hier(&mut deps);
    find_depressions_to_fill(dem, cols, &label, &deps, &mut wtd);

    Ok(FsmResult {
        ocean_volume: deps[OCEAN as usize].water_vol,
        depression_count: deps.len(),
        wtd,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    // Sum of standing water across the grid.
    fn total_water(wtd: &[f64]) -> f64 {
        wtd.iter().sum()
    }

    #[test]
    fn single_pit_partial_fill() {
        // 5x5 bowl draining nowhere internally; a single deep pit at the centre.
        // Rim = 10, interior ring = 5, pit = 0. Border cells are outlets.
        let rows = 5;
        let cols = 5;
        #[rustfmt::skip]
        let dem = vec![
            10.0,10.0,10.0,10.0,10.0,
            10.0, 5.0, 5.0, 5.0,10.0,
            10.0, 5.0, 0.0, 5.0,10.0,
            10.0, 5.0, 5.0, 5.0,10.0,
            10.0,10.0,10.0,10.0,10.0,
        ];
        // Add a little water everywhere; not enough to overflow the rim.
        let water = vec![0.2; rows * cols];
        let res = fill_spill_merge(&dem, rows, cols, -9999.0, &water, None, true).unwrap();

        // The interior 3x3 depression collects water; it should pool, not vanish.
        assert!(total_water(&res.wtd) > 0.0);
        // Water only accumulates on interior cells (the border is the outlet).
        for r in 0..rows {
            for c in 0..cols {
                if r == 0 || c == 0 || r == rows - 1 || c == cols - 1 {
                    assert_eq!(res.wtd[r * cols + c], 0.0, "border must stay dry");
                }
            }
        }
        // The deepest cell must hold the most water.
        let center = 2 * cols + 2;
        assert!(res.wtd[center] > 0.0);
    }

    #[test]
    fn excess_water_drains_to_ocean() {
        // The same bowl, but flood it with far more water than it can hold. The
        // depression fills to its rim and the rest is lost to the ocean sink.
        let rows = 5;
        let cols = 5;
        #[rustfmt::skip]
        let dem = vec![
            10.0,10.0,10.0,10.0,10.0,
            10.0, 5.0, 5.0, 5.0,10.0,
            10.0, 5.0, 0.0, 5.0,10.0,
            10.0, 5.0, 5.0, 5.0,10.0,
            10.0,10.0,10.0,10.0,10.0,
        ];
        let water = vec![100.0; rows * cols];
        let res = fill_spill_merge(&dem, rows, cols, -9999.0, &water, None, true).unwrap();

        // A lot of water should have drained off the grid.
        assert!(res.ocean_volume > 0.0);
        // No interior cell should hold water above the rim elevation (10).
        for r in 1..rows - 1 {
            for c in 1..cols - 1 {
                let i = r * cols + c;
                let surface = dem[i] + res.wtd[i];
                assert!(
                    surface <= 10.0 + 1e-3,
                    "cell {i} surface {surface} exceeds rim"
                );
            }
        }
    }

    #[test]
    fn flat_plane_holds_no_water() {
        // A perfectly flat, freely draining plane: all water runs off the edges.
        let rows = 6;
        let cols = 6;
        let dem = vec![5.0; rows * cols];
        let water = vec![1.0; rows * cols];
        let res = fill_spill_merge(&dem, rows, cols, -9999.0, &water, None, true).unwrap();
        assert!(
            total_water(&res.wtd) < 1e-6,
            "flat plane should retain no standing water"
        );
    }

    #[test]
    fn mass_is_conserved() {
        // Total input water == standing water + water lost to the ocean.
        let rows = 7;
        let cols = 7;
        let mut dem = vec![20.0; rows * cols];
        // Two separate pits of differing depth.
        dem[2 * cols + 2] = 0.0;
        dem[2 * cols + 3] = 8.0;
        dem[3 * cols + 2] = 8.0;
        dem[3 * cols + 3] = 8.0;
        dem[4 * cols + 4] = 2.0;
        dem[4 * cols + 5] = 9.0;
        dem[5 * cols + 4] = 9.0;
        let depth = 0.5;
        let water = vec![depth; rows * cols];
        let res = fill_spill_merge(&dem, rows, cols, -9999.0, &water, None, true).unwrap();

        // Input water lands only on non-border cells (border is instant outlet).
        let interior = (rows - 2) * (cols - 2);
        let input = interior as f64 * depth;
        let standing = total_water(&res.wtd);
        let accounted = standing + res.ocean_volume;
        assert!(
            (accounted - input).abs() < 1e-3,
            "mass not conserved: input={input} standing={standing} ocean={} accounted={accounted}",
            res.ocean_volume
        );
    }

    #[test]
    fn rejects_mismatched_buffers() {
        let dem = vec![1.0; 9];
        let water = vec![1.0; 4];
        assert!(fill_spill_merge(&dem, 3, 3, -9999.0, &water, None, true).is_err());
    }
}
