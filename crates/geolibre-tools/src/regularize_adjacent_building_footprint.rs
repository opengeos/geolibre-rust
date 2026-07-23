//! GeoLibre tool: regularize a *group* of adjacent building footprints jointly.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Regularize Adjacent Building
//! Footprint* (3D Analyst). The authored `regularize_building_footprints`
//! regularizes each polygon **independently**, so it picks a per-building
//! dominant direction and adjacent footprints that share a wall can drift out
//! of alignment. This tool regularizes each **group** to one shared
//! orientation so shared walls stay parallel and (via offset snapping)
//! collinear.
//!
//! Per group:
//! 1. A shared dominant orientation θ — the wall-length-weighted circular mean
//!    of every edge direction (taken mod 90°, since a rectangular building's
//!    walls lie at θ and θ+90).
//! 2. Each footprint is regularized against the **group** axis set: `right_angles`
//!    uses {θ, θ+90}; `right_angles_and_diagonals` adds {θ+45, θ+135}. Every
//!    edge snaps to the nearest allowed direction; the wall is rebuilt as the
//!    line through its (grid-snapped) offset, and new vertices are the
//!    intersections of consecutive wall lines.
//! 3. Wall offsets are snapped to a shared `precision` grid in the group frame,
//!    so near-coincident parallel walls of neighbouring buildings land on the
//!    same line (collinear shared edges).
//!
//! Groups come from a `group` field, or (absent one) from boundary adjacency:
//! footprints whose vertices come within `adjacency_distance` are one group.
//! Short edges (< `tolerance`) are pruned first as segmentation noise.
//! Distances are in the layer's CRS units.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldValue, Geometry, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct RegularizeAdjacentBuildingFootprintTool;

impl Tool for RegularizeAdjacentBuildingFootprintTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "regularize_adjacent_building_footprint",
            display_name: "Regularize Adjacent Building Footprint",
            summary: "Regularize a group of adjacent building footprints to one shared orientation so shared walls stay parallel and collinear (like ArcGIS Regularize Adjacent Building Footprint) — the group-consistent counterpart of the per-building regularize_building_footprints, which lets neighbouring footprints drift out of alignment.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input building-footprint polygon layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "group",
                    description: "Optional field grouping footprints to regularize together. Without it, groups are formed by boundary adjacency.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'right_angles' (default) or 'right_angles_and_diagonals'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Minimum wall length; shorter edges are pruned as noise (CRS units). Default 1.0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "precision",
                    description: "Grid step for snapping wall offsets so nearby parallel walls become collinear (CRS units). Default: tolerance.",
                    required: false,
                },
                ToolParamSpec {
                    name: "adjacency_distance",
                    description: "When no 'group' field is given, footprints whose vertices come within this distance are one group (CRS units). Default: 2 x tolerance.",
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

        let mut layer = load_input_layer(input)?;
        let group_idx =
            match &prm.group {
                Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                    ToolError::Validation(format!("group field '{f}' not found"))
                })?),
                None => None,
            };

        // Collect polygon footprints (exterior ring only for orientation; holes
        // are preserved but not group-aligned).
        let mut foot: Vec<Footprint> = Vec::new();
        for (fi, feature) in layer.features.iter().enumerate() {
            if let Some((ext, holes)) = polygon_rings(feature.geometry.as_ref()) {
                let gk = group_idx.map(|gi| field_key(&feature.attributes[gi]));
                foot.push(Footprint {
                    feat: fi,
                    ext,
                    holes,
                    group: gk,
                });
            }
        }
        if foot.is_empty() {
            return Err(ToolError::Execution(
                "no polygon footprints in input".to_string(),
            ));
        }

        // ── Group footprints ──────────────────────────────────────────────────
        let groups = build_groups(&foot, &prm);
        ctx.progress.info(&format!(
            "{} footprint(s) in {} group(s)",
            foot.len(),
            groups.len()
        ));

        // ── Regularize each group to a shared orientation ─────────────────────
        let mut new_geom: BTreeMap<usize, Geometry> = BTreeMap::new();
        let mut changed = 0usize;
        for members in groups.values() {
            let theta = group_orientation(&foot, members);
            for &m in members {
                let f = &foot[m];
                let ext = regularize_ring(&f.ext, theta, &prm);
                if let Some(ext) = ext {
                    let holes: Vec<Ring> = f
                        .holes
                        .iter()
                        .filter_map(|h| regularize_ring(h, theta, &prm))
                        .collect();
                    changed += 1;
                    new_geom.insert(
                        f.feat,
                        Geometry::Polygon {
                            exterior: ext,
                            interiors: holes,
                        },
                    );
                } else {
                    // Keep original geometry when regularization degenerates.
                }
            }
        }

        // Write regularized geometries back onto a copy of the input.
        for (fi, feature) in layer.features.iter_mut().enumerate() {
            if let Some(g) = new_geom.remove(&fi) {
                feature.geometry = Some(g);
            }
        }

        let feature_count = layer.len();
        let group_count = groups.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("group_count".to_string(), json!(group_count));
        outputs.insert("changed_count".to_string(), json!(changed));
        outputs.insert("method".to_string(), json!(prm.method.as_str()));
        Ok(ToolRunResult { outputs })
    }
}

// ── Footprints & grouping ─────────────────────────────────────────────────────

struct Footprint {
    feat: usize,
    ext: Vec<(f64, f64)>,
    holes: Vec<Vec<(f64, f64)>>,
    group: Option<String>,
}

/// Groups footprints by the `group` field, or by boundary adjacency when no
/// field is given. Returns group-key -> member indices.
fn build_groups(foot: &[Footprint], prm: &Params) -> BTreeMap<String, Vec<usize>> {
    if foot.iter().any(|f| f.group.is_some()) {
        let mut g: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for (i, f) in foot.iter().enumerate() {
            let key = f.group.clone().unwrap_or_else(|| format!("__nogroup_{i}"));
            g.entry(key).or_default().push(i);
        }
        return g;
    }
    // Adjacency grouping via union-find on vertex proximity.
    let n = foot.len();
    let mut uf = UnionFind::new(n);
    let bboxes: Vec<[f64; 4]> = foot.iter().map(|f| bbox(&f.ext)).collect();
    let d = prm.adjacency_distance;
    for a in 0..n {
        for b in (a + 1)..n {
            if !bbox_within(&bboxes[a], &bboxes[b], d) {
                continue;
            }
            if rings_touch(&foot[a].ext, &foot[b].ext, d) {
                uf.union(a, b);
            }
        }
    }
    let mut g: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        g.entry(format!("adj_{}", uf.find(i))).or_default().push(i);
    }
    g
}

fn rings_touch(a: &[(f64, f64)], b: &[(f64, f64)], d: f64) -> bool {
    let d2 = d * d;
    for &(ax, ay) in a {
        for &(bx, by) in b {
            let (dx, dy) = (ax - bx, ay - by);
            if dx * dx + dy * dy <= d2 {
                return true;
            }
        }
    }
    false
}

// ── Orientation & regularization ──────────────────────────────────────────────

/// Wall-length-weighted dominant orientation θ (radians, in [0, π/2)) shared by
/// a group. Uses the 4θ circular mean because right-angle walls at θ and θ+90°
/// are equivalent under a π/2 period.
fn group_orientation(foot: &[Footprint], members: &[usize]) -> f64 {
    let (mut sc, mut ss) = (0.0f64, 0.0f64);
    for &m in members {
        let ring = &foot[m].ext;
        let n = ring.len();
        for i in 0..n {
            let (x0, y0) = ring[i];
            let (x1, y1) = ring[(i + 1) % n];
            let (dx, dy) = (x1 - x0, y1 - y0);
            let len = dx.hypot(dy);
            if len <= 0.0 {
                continue;
            }
            let a = dy.atan2(dx); // edge angle
            sc += len * (4.0 * a).cos();
            ss += len * (4.0 * a).sin();
        }
    }
    if sc == 0.0 && ss == 0.0 {
        return 0.0;
    }
    let theta = ss.atan2(sc) / 4.0;
    // Normalize into [0, π/2).
    let half = std::f64::consts::FRAC_PI_2;
    ((theta % half) + half) % half
}

/// Regularize one ring against the group orientation. Returns None if the
/// rebuilt ring degenerates (< 3 vertices).
fn regularize_ring(ring: &[(f64, f64)], theta: f64, prm: &Params) -> Option<Ring> {
    // 1. Prune short edges (segmentation noise): drop a vertex whose incoming
    //    edge is shorter than tolerance.
    let pruned = prune_short_edges(ring, prm.tolerance);
    if pruned.len() < 3 {
        return None;
    }
    let n = pruned.len();

    // 2. Allowed wall directions (angles mod π) for the chosen method.
    let allowed = allowed_angles(theta, prm.method);

    // 3. Each edge -> a wall line: snapped direction + grid-snapped offset.
    let mut lines: Vec<Line> = Vec::with_capacity(n);
    for i in 0..n {
        let (x0, y0) = pruned[i];
        let (x1, y1) = pruned[(i + 1) % n];
        let a_edge = (y1 - y0).atan2(x1 - x0);
        let a = snap_angle(a_edge, &allowed);
        // Unit normal to direction a.
        let (nx, ny) = (-a.sin(), a.cos());
        let mid = ((x0 + x1) * 0.5, (y0 + y1) * 0.5);
        let mut offset = mid.0 * nx + mid.1 * ny;
        // Snap offset to the shared precision grid for collinearity.
        offset = (offset / prm.precision).round() * prm.precision;
        lines.push(Line { nx, ny, offset });
    }

    // 4. New vertices = intersections of consecutive wall lines.
    let mut out: Vec<Coord> = Vec::with_capacity(n);
    for i in 0..n {
        let prev = &lines[(i + n - 1) % n];
        let cur = &lines[i];
        match intersect(prev, cur) {
            Some((x, y)) => out.push(Coord::xy(x, y)),
            None => return None, // parallel consecutive walls -> degenerate
        }
    }
    // Drop near-duplicate consecutive vertices.
    dedup_ring(&mut out);
    if out.len() < 3 {
        return None;
    }
    Some(Ring::new(out))
}

struct Line {
    nx: f64,
    ny: f64,
    offset: f64,
}

/// Intersection of two lines p·n = offset. None if (near) parallel.
fn intersect(a: &Line, b: &Line) -> Option<(f64, f64)> {
    let det = a.nx * b.ny - a.ny * b.nx;
    if det.abs() < 1e-12 {
        return None;
    }
    let x = (a.offset * b.ny - a.ny * b.offset) / det;
    let y = (a.nx * b.offset - a.offset * b.nx) / det;
    Some((x, y))
}

/// Allowed wall angles for the method: {θ, θ+90} or, with diagonals, also
/// {θ+45, θ+135}. Angles are returned mod π (line orientations).
fn allowed_angles(theta: f64, method: Method) -> Vec<f64> {
    let mut v = vec![theta, theta + std::f64::consts::FRAC_PI_2];
    if method == Method::RightAnglesAndDiagonals {
        v.push(theta + std::f64::consts::FRAC_PI_4);
        v.push(theta + 3.0 * std::f64::consts::FRAC_PI_4);
    }
    v.iter()
        .map(|a| a.rem_euclid(std::f64::consts::PI))
        .collect()
}

/// Snaps an edge angle to the nearest allowed line orientation (mod π).
fn snap_angle(a_edge: f64, allowed: &[f64]) -> f64 {
    let pi = std::f64::consts::PI;
    let a = a_edge.rem_euclid(pi);
    let mut best = allowed[0];
    let mut best_d = f64::INFINITY;
    for &c in allowed {
        // Circular distance mod π.
        let mut diff = (a - c).abs() % pi;
        if diff > pi / 2.0 {
            diff = pi - diff;
        }
        if diff < best_d {
            best_d = diff;
            best = c;
        }
    }
    best
}

/// Removes vertices whose incoming edge is shorter than `tolerance`, keeping the
/// ring closed. Never reduces below 3 vertices.
fn prune_short_edges(ring: &[(f64, f64)], tolerance: f64) -> Vec<(f64, f64)> {
    // Strip an explicit closing duplicate if present.
    let mut pts: Vec<(f64, f64)> = ring.to_vec();
    if pts.len() >= 2 && pts.first() == pts.last() {
        pts.pop();
    }
    if pts.len() <= 3 {
        return pts;
    }
    let mut kept: Vec<(f64, f64)> = Vec::with_capacity(pts.len());
    for &p in &pts {
        if let Some(&last) = kept.last() {
            let (dx, dy) = (p.0 - last.0, p.1 - last.1);
            if dx.hypot(dy) < tolerance {
                continue;
            }
        }
        kept.push(p);
    }
    // Also check the closing edge (last -> first).
    if kept.len() > 3 {
        if let (Some(&first), Some(&last)) = (kept.first(), kept.last()) {
            let (dx, dy) = (first.0 - last.0, first.1 - last.1);
            if dx.hypot(dy) < tolerance {
                kept.pop();
            }
        }
    }
    if kept.len() < 3 {
        pts
    } else {
        kept
    }
}

fn dedup_ring(pts: &mut Vec<Coord>) {
    let eps = 1e-6;
    let mut out: Vec<Coord> = Vec::with_capacity(pts.len());
    for p in pts.iter() {
        if let Some(last) = out.last() {
            if (last.x - p.x).abs() < eps && (last.y - p.y).abs() < eps {
                continue;
            }
        }
        out.push(Coord::xy(p.x, p.y));
    }
    if out.len() >= 2 {
        let (fx, fy) = (out[0].x, out[0].y);
        let l = &out[out.len() - 1];
        if (fx - l.x).abs() < eps && (fy - l.y).abs() < eps {
            out.pop();
        }
    }
    *pts = out;
}

// ── Geometry helpers ──────────────────────────────────────────────────────────

fn polygon_rings(geom: Option<&Geometry>) -> Option<(Vec<(f64, f64)>, Vec<Vec<(f64, f64)>>)> {
    match geom? {
        Geometry::Polygon {
            exterior,
            interiors,
        } => Some((
            ring_coords(exterior),
            interiors.iter().map(ring_coords).collect(),
        )),
        Geometry::MultiPolygon(parts) if !parts.is_empty() => {
            // Regularize the largest part; keep others as holes-free parts is
            // out of scope, so take the first part's rings.
            let (e, i) = &parts[0];
            Some((ring_coords(e), i.iter().map(ring_coords).collect()))
        }
        _ => None,
    }
}

fn ring_coords(r: &Ring) -> Vec<(f64, f64)> {
    r.coords().iter().map(|c| (c.x, c.y)).collect()
}

fn bbox(ring: &[(f64, f64)]) -> [f64; 4] {
    let mut b = [
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    ];
    for &(x, y) in ring {
        b[0] = b[0].min(x);
        b[1] = b[1].min(y);
        b[2] = b[2].max(x);
        b[3] = b[3].max(y);
    }
    b
}

fn bbox_within(a: &[f64; 4], b: &[f64; 4], d: f64) -> bool {
    a[0] - d <= b[2] && b[0] - d <= a[2] && a[1] - d <= b[3] && b[1] - d <= a[3]
}

fn field_key(fv: &FieldValue) -> String {
    if let Some(i) = fv.as_i64() {
        i.to_string()
    } else if let Some(f) = fv.as_f64() {
        format!("{f}")
    } else {
        fv.as_str().unwrap_or("").to_string()
    }
}

struct UnionFind {
    parent: Vec<usize>,
}
impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
        }
    }
    fn find(&mut self, x: usize) -> usize {
        let mut r = x;
        while self.parent[r] != r {
            r = self.parent[r];
        }
        let mut c = x;
        while self.parent[c] != r {
            let next = self.parent[c];
            self.parent[c] = r;
            c = next;
        }
        r
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.parent[rb] = ra;
        }
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    RightAngles,
    RightAnglesAndDiagonals,
}

impl Method {
    fn as_str(&self) -> &'static str {
        match self {
            Method::RightAngles => "right_angles",
            Method::RightAnglesAndDiagonals => "right_angles_and_diagonals",
        }
    }
}

struct Params {
    group: Option<String>,
    method: Method,
    tolerance: f64,
    precision: f64,
    adjacency_distance: f64,
}

fn parse_optional_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(Some(n.as_f64().unwrap_or(f64::NAN))),
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

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let group = parse_optional_str(args, "group")?.map(str::to_string);
    let method = match parse_optional_str(args, "method")?.map(str::trim) {
        None | Some("") | Some("right_angles") => Method::RightAngles,
        Some("right_angles_and_diagonals") => Method::RightAnglesAndDiagonals,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'method' must be 'right_angles' or 'right_angles_and_diagonals', got '{o}'"
            )))
        }
    };
    let tolerance = parse_optional_f64(args, "tolerance")?.unwrap_or(1.0);
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "'tolerance' must be a positive number".to_string(),
        ));
    }
    let precision = match parse_optional_f64(args, "precision")? {
        None => tolerance,
        Some(v) if v > 0.0 => v,
        Some(_) => {
            return Err(ToolError::Validation(
                "'precision' must be a positive number".to_string(),
            ))
        }
    };
    let adjacency_distance = match parse_optional_f64(args, "adjacency_distance")? {
        None => 2.0 * tolerance,
        Some(v) if v > 0.0 => v,
        Some(_) => {
            return Err(ToolError::Validation(
                "'adjacency_distance' must be a positive number".to_string(),
            ))
        }
    };
    Ok(Params {
        group,
        method,
        tolerance,
        precision,
        adjacency_distance,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Rotate a point by angle (radians) about the origin.
    fn rot(x: f64, y: f64, a: f64) -> (f64, f64) {
        (x * a.cos() - y * a.sin(), x * a.sin() + y * a.cos())
    }

    /// Builds a layer from rings; each ring is (vertices, group).
    fn layer_of(rings: Vec<(Vec<(f64, f64)>, &str)>) -> String {
        let mut l = Layer::new("b")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(32618);
        l.add_field(FieldDef::new("grp", FieldType::Text));
        for (verts, g) in rings {
            let ring: Vec<Coord> = verts.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
            l.add_feature(Some(Geometry::polygon(ring, vec![])), &[("grp", g.into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = RegularizeAdjacentBuildingFootprintTool
            .run(&args, &ctx())
            .unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// The set of distinct wall line-orientations (mod 180°, degrees, rounded)
    /// over a ring — a rectilinear building has exactly two, θ and θ+90.
    fn wall_orientations(g: &Geometry) -> Vec<i64> {
        let (ext, _) = polygon_rings(Some(g)).unwrap();
        let n = ext.len();
        let mut set = std::collections::BTreeSet::new();
        for i in 0..n {
            let (x0, y0) = ext[i];
            let (x1, y1) = ext[(i + 1) % n];
            let a = (y1 - y0).atan2(x1 - x0).to_degrees().rem_euclid(180.0);
            set.insert(a.round() as i64 % 180);
        }
        set.into_iter().collect()
    }

    /// Two near-rectangular footprints sharing a wall, both rotated ~20°, become
    /// exactly rectilinear against a single shared orientation.
    #[test]
    fn shared_orientation_makes_walls_parallel() {
        let a = 20f64.to_radians();
        // Two 10x10 boxes side by side (share the x=10 wall), slightly noisy,
        // then rotated by 20°.
        let mk = |x0: f64| {
            vec![
                rot(x0, 0.0, a),
                rot(x0 + 10.2, 0.3, a),
                rot(x0 + 9.8, 10.1, a),
                rot(x0, 9.9, a),
            ]
        };
        let input = layer_of(vec![(mk(0.0), "blk"), (mk(10.0), "blk")]);
        let (out, layer) = run(json!({
            "input": input, "group": "grp", "method": "right_angles", "tolerance": 0.5
        }));
        assert_eq!(out.outputs["group_count"], json!(1));
        // Both footprints share the same two wall orientations.
        let w0 = wall_orientations(layer.features[0].geometry.as_ref().unwrap());
        let w1 = wall_orientations(layer.features[1].geometry.as_ref().unwrap());
        assert_eq!(w0, w1, "the two footprints must share wall orientations");
        assert_eq!(
            w0.len(),
            2,
            "a rectilinear building has 2 wall orientations"
        );
    }

    /// Auto-adjacency grouping: two touching footprints (no group field) are
    /// found to be one group and aligned together.
    #[test]
    fn auto_adjacency_grouping() {
        let a = 15f64.to_radians();
        let mk = |x0: f64| {
            vec![
                rot(x0, 0.0, a),
                rot(x0 + 10.0, 0.2, a),
                rot(x0 + 10.1, 10.0, a),
                rot(x0, 9.8, a),
            ]
        };
        // share the x=10 edge (touching within adjacency distance).
        let input = layer_of(vec![(mk(0.0), ""), (mk(10.0), "")]);
        let (out, _l) = run(json!({
            "input": input, "method": "right_angles", "tolerance": 0.5, "adjacency_distance": 1.0
        }));
        assert_eq!(
            out.outputs["group_count"],
            json!(1),
            "touching footprints form one group"
        );
        assert_eq!(out.outputs["changed_count"], json!(2));
    }

    /// Area is roughly preserved (regularization shifts walls by < tolerance).
    #[test]
    fn area_roughly_preserved() {
        let verts = vec![(0.0, 0.0), (20.3, 0.2), (19.8, 10.1), (0.1, 9.9)];
        let orig_area = shoelace(&verts);
        let input = layer_of(vec![(verts, "g")]);
        let (_o, layer) = run(json!({ "input": input, "group": "grp", "tolerance": 0.5 }));
        let (ext, _) = polygon_rings(layer.features[0].geometry.as_ref()).unwrap();
        let new_area = shoelace(&ext);
        assert!(
            (new_area - orig_area).abs() / orig_area < 0.1,
            "area changed too much: {orig_area} -> {new_area}"
        );
    }

    fn shoelace(pts: &[(f64, f64)]) -> f64 {
        let n = pts.len();
        let mut s = 0.0;
        for i in 0..n {
            let (x0, y0) = pts[i];
            let (x1, y1) = pts[(i + 1) % n];
            s += x0 * y1 - x1 * y0;
        }
        s.abs() / 2.0
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            RegularizeAdjacentBuildingFootprintTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "method": "curvy" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "tolerance": -1 })).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_ok());
    }
}
