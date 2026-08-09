//! Shared triangle-mesh machinery for GeoLibre's 3D tools.
//!
//! `inside_3d` introduced the mesh representation (`Tri`, `Solid`,
//! `collect_triangles`) and `union_3d` added exact volumes plus voxel
//! occupancy. Round 17 adds eight more 3D tools that all need the same few
//! primitives — an edge map, boundary-loop extraction, triangle emission, and
//! occupancy sampling — so they live here instead of being copied eight times.
//!
//! ## The edge map is the centre of this module
//!
//! Almost every question a triangle mesh raises reduces to "how many triangles
//! share this edge, and in which direction":
//!
//! * exactly two, everywhere → the mesh is **closed** (`is_closed_3d`);
//! * exactly one → that edge is on a **boundary loop**, which is what
//!   `enclose_multipatch` caps;
//! * two, traversed in opposite directions → **consistent winding**, without
//!   which the signed-tetrahedron volume is meaningless;
//! * two whose normals disagree against a light direction → a **silhouette**
//!   edge, which is what `sun_shadow_volume` extrudes.
//!
//! Vertices are quantised before keying so float noise from mesh construction
//! does not split an edge that is geometrically shared — the same 1e-9
//! quantisation `inside_3d::is_closed` uses.

use std::collections::BTreeMap;

use wbvector::{Coord, Geometry, Ring};

use crate::inside_3d::{Solid, Tri};

/// Quantised vertex key. Coordinates are snapped to a 1e-9 grid: fine enough
/// not to merge distinct vertices, coarse enough to survive construction noise.
pub(crate) type VKey = (i64, i64, i64);

/// Undirected edge key, always stored with the smaller vertex first.
pub(crate) type EdgeKey = (VKey, VKey);

pub(crate) fn vkey(p: [f64; 3]) -> VKey {
    (
        (p[0] / 1e-9).round() as i64,
        (p[1] / 1e-9).round() as i64,
        (p[2] / 1e-9).round() as i64,
    )
}

pub(crate) fn edge_key(a: [f64; 3], b: [f64; 3]) -> EdgeKey {
    let (ka, kb) = (vkey(a), vkey(b));
    if ka <= kb {
        (ka, kb)
    } else {
        (kb, ka)
    }
}

/// One triangle's use of an edge: which triangle, and whether it traversed the
/// edge in the key's canonical direction.
#[derive(Clone, Copy)]
pub(crate) struct EdgeUse {
    pub(crate) tri: usize,
    /// True when this triangle walks the edge low-key → high-key.
    pub(crate) forward: bool,
}

/// Maps every undirected edge to the triangles using it.
pub(crate) fn edge_map(tris: &[Tri]) -> BTreeMap<EdgeKey, Vec<EdgeUse>> {
    let mut map: BTreeMap<EdgeKey, Vec<EdgeUse>> = BTreeMap::new();
    for (ti, t) in tris.iter().enumerate() {
        for k in 0..3 {
            let p = t[k];
            let q = t[(k + 1) % 3];
            let key = edge_key(p, q);
            let forward = vkey(p) <= vkey(q);
            map.entry(key)
                .or_default()
                .push(EdgeUse { tri: ti, forward });
        }
    }
    map
}

/// Watertightness and orientation diagnostics for a mesh.
pub(crate) struct MeshTopology {
    /// Edges used by exactly one triangle — the mesh's open boundary.
    pub(crate) open_edges: usize,
    /// Edges used by three or more triangles (non-manifold).
    pub(crate) nonmanifold_edges: usize,
    /// True when every edge is used by exactly two triangles.
    pub(crate) closed: bool,
    /// True when every shared edge is traversed in opposite directions by its
    /// two triangles — the condition that makes the signed volume meaningful.
    pub(crate) consistent_winding: bool,
}

pub(crate) fn topology(tris: &[Tri]) -> MeshTopology {
    let map = edge_map(tris);
    let mut open = 0usize;
    let mut nonmanifold = 0usize;
    let mut consistent = true;
    for uses in map.values() {
        match uses.len() {
            0 => {}
            1 => open += 1,
            2 => {
                // Opposite traversal directions mean the two triangles agree on
                // which side is "outside".
                if uses[0].forward == uses[1].forward {
                    consistent = false;
                }
            }
            _ => nonmanifold += 1,
        }
    }
    MeshTopology {
        closed: !map.is_empty() && open == 0 && nonmanifold == 0,
        open_edges: open,
        nonmanifold_edges: nonmanifold,
        consistent_winding: consistent && !map.is_empty(),
    }
}

/// Chains the mesh's open (single-use) edges into closed boundary loops.
///
/// Returns loops as vertex-position rings. An edge whose partner is missing —
/// a genuinely dangling boundary — simply ends its loop; the loop is only
/// emitted when it closes back on its start, so `enclose_multipatch` never
/// tries to cap a chain that is not actually a hole.
pub(crate) fn boundary_loops(tris: &[Tri]) -> Vec<Vec<[f64; 3]>> {
    // Collect directed boundary edges, keeping real coordinates alongside keys
    // so the emitted loop carries unquantised positions.
    let mut next: BTreeMap<VKey, Vec<([f64; 3], [f64; 3])>> = BTreeMap::new();
    for (key, uses) in edge_map(tris) {
        if uses.len() != 1 {
            continue;
        }
        let t = &tris[uses[0].tri];
        // Recover the directed edge as this triangle walked it.
        let mut directed = None;
        for k in 0..3 {
            let (p, q) = (t[k], t[(k + 1) % 3]);
            if edge_key(p, q) == key {
                directed = Some((p, q));
                break;
            }
        }
        if let Some((p, q)) = directed {
            next.entry(vkey(p)).or_default().push((p, q));
        }
    }

    let mut loops = Vec::new();
    while let Some((&start_key, _)) = next.iter().find(|(_, v)| !v.is_empty()) {
        let mut ring: Vec<[f64; 3]> = Vec::new();
        let mut cursor = start_key;
        loop {
            let Some(bucket) = next.get_mut(&cursor) else {
                break;
            };
            let Some((p, q)) = bucket.pop() else { break };
            ring.push(p);
            cursor = vkey(q);
            if cursor == start_key {
                // Closed the loop.
                if ring.len() >= 3 {
                    loops.push(std::mem::take(&mut ring));
                }
                break;
            }
            // Guard against a pathological walk revisiting forever.
            if ring.len() > tris.len() * 3 + 3 {
                break;
            }
        }
        next.retain(|_, v| !v.is_empty());
        if next.is_empty() {
            break;
        }
    }
    loops
}

/// Triangulates a (possibly non-planar) boundary loop with a centroid fan.
///
/// A fan from the *centroid* rather than from vertex 0 is what makes this work
/// on non-planar loops: vertex-0 fans on a saddle-shaped loop produce
/// self-intersecting triangles, while the centroid fan stays well-formed and,
/// crucially, contributes the correct signed volume.
///
/// `reverse` flips the emitted winding so the cap faces outward relative to the
/// mesh it closes.
pub(crate) fn fan_triangulate(ring: &[[f64; 3]], reverse: bool) -> Vec<Tri> {
    if ring.len() < 3 {
        return Vec::new();
    }
    let n = ring.len() as f64;
    let mut c = [0.0; 3];
    for p in ring {
        for k in 0..3 {
            c[k] += p[k] / n;
        }
    }
    let mut out = Vec::with_capacity(ring.len());
    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        out.push(if reverse { [c, b, a] } else { [c, a, b] });
    }
    out
}

/// Emits triangles as a `MultiPolygon` with one triangle per part — the
/// convention `buffer_3d`, `minimum_bounding_volume` and `voxel_isosurface`
/// already share, so the result is valid input to every other 3D tool.
///
/// The ring carries **three** coordinates and is deliberately not closed with a
/// repeated first vertex: `inside_3d::collect_triangles` fans each ring from
/// its first coordinate, so a closing vertex would turn every triangle into two
/// — one real and one degenerate — and silently double the mesh on round trip.
pub(crate) fn triangles_to_geometry(tris: &[Tri]) -> Geometry {
    Geometry::MultiPolygon(
        tris.iter()
            .map(|t| {
                (
                    Ring::new(
                        t.iter()
                            .map(|p| Coord::xyz(p[0], p[1], p[2]))
                            .collect::<Vec<_>>(),
                    ),
                    Vec::new(),
                )
            })
            .collect(),
    )
}

/// Exact volume of a closed triangle mesh by signed-tetrahedron summation.
pub(crate) fn mesh_volume(tris: &[Tri]) -> f64 {
    let mut v = 0.0;
    for t in tris {
        let (a, b, c) = (t[0], t[1], t[2]);
        v += a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]);
    }
    (v / 6.0).abs()
}

/// Total surface area of a triangle mesh.
pub(crate) fn mesh_area(tris: &[Tri]) -> f64 {
    tris.iter().map(tri_area).sum()
}

/// Area of one triangle (half the cross-product magnitude).
pub(crate) fn tri_area(t: &Tri) -> f64 {
    let u = sub3(t[1], t[0]);
    let v = sub3(t[2], t[0]);
    let n = cross3(u, v);
    0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
}

/// Unit normal of a triangle, or `None` for a degenerate one.
pub(crate) fn tri_normal(t: &Tri) -> Option<[f64; 3]> {
    let n = cross3(sub3(t[1], t[0]), sub3(t[2], t[0]));
    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    (len > 1e-15).then(|| [n[0] / len, n[1] / len, n[2] / len])
}

pub(crate) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub(crate) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub(crate) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Bounding box shared by two solids, or `None` when they do not overlap.
pub(crate) fn intersect_bbox(a: &Solid, b: &Solid) -> Option<([f64; 3], [f64; 3])> {
    let mut min = [0.0_f64; 3];
    let mut max = [0.0_f64; 3];
    for k in 0..3 {
        min[k] = a.min[k].max(b.min[k]);
        max[k] = a.max[k].min(b.max[k]);
        if max[k] <= min[k] {
            return None;
        }
    }
    Some((min, max))
}

pub(crate) fn bbox_overlap(a: &Solid, b: &Solid) -> bool {
    (0..3).all(|k| a.min[k] <= b.max[k] && a.max[k] >= b.min[k])
}

/// Builds a sampling grid over a box: the longest axis gets `resolution` cells
/// and the others are scaled to keep voxels cubic, so accuracy does not depend
/// on the box's aspect ratio.
pub(crate) fn grid_for(
    min: [f64; 3],
    max: [f64; 3],
    resolution: usize,
) -> (usize, usize, usize, f64, [f64; 3]) {
    let span = [max[0] - min[0], max[1] - min[1], max[2] - min[2]];
    let longest = span[0].max(span[1]).max(span[2]);
    // A zero resolution would divide by zero; every tool validates it, but the
    // helper is crate-visible so it defends itself.
    if longest <= 0.0 || resolution == 0 {
        return (0, 0, 0, 0.0, [0.0; 3]);
    }
    let h = longest / resolution as f64;
    let n: Vec<usize> = span
        .iter()
        .map(|s| ((s / h).ceil() as usize).max(1))
        .collect();
    let step = [
        span[0] / n[0] as f64,
        span[1] / n[1] as f64,
        span[2] / n[2] as f64,
    ];
    (n[0], n[1], n[2], step[0] * step[1] * step[2], step)
}

/// Volume of the region satisfying `predicate`, sampled over `min..max`.
///
/// This is the generalisation of `union_3d`'s occupancy counter: the caller
/// supplies the set membership test, so union (any), intersection (all) and
/// difference (first and not the rest) all share one sampler and therefore one
/// accuracy story.
pub(crate) fn occupancy_volume<F>(
    min: [f64; 3],
    max: [f64; 3],
    resolution: usize,
    predicate: F,
) -> f64
where
    F: Fn(f64, f64, f64) -> bool,
{
    let (nx, ny, nz, cell_vol, step) = grid_for(min, max, resolution);
    if cell_vol <= 0.0 {
        return 0.0;
    }
    let mut occupied = 0_u64;
    for k in 0..nz {
        let z = min[2] + (k as f64 + 0.5) * step[2];
        for j in 0..ny {
            let y = min[1] + (j as f64 + 0.5) * step[1];
            for i in 0..nx {
                let x = min[0] + (i as f64 + 0.5) * step[0];
                if predicate(x, y, z) {
                    occupied += 1;
                }
            }
        }
    }
    occupied as f64 * cell_vol
}

/// Möller-Trumbore **segment**/triangle intersection.
///
/// Returns the parameter in `[0, 1]` along `a -> b` where the segment meets the
/// triangle. Unlike `inside_3d`'s infinite-ray version there is no parity
/// counting here, so the shared-edge double-hit hazard does not apply: callers
/// want "is anything in the way", not a crossing count.
pub(crate) fn segment_triangle(a: [f64; 3], b: [f64; 3], tri: &Tri) -> Option<f64> {
    const EPS: f64 = 1e-12;
    let dir = sub3(b, a);
    let e1 = sub3(tri[1], tri[0]);
    let e2 = sub3(tri[2], tri[0]);
    let p = cross3(dir, e2);
    let det = dot3(e1, p);
    if det.abs() < EPS {
        return None; // parallel to the triangle's plane
    }
    let inv = 1.0 / det;
    let t_vec = sub3(a, tri[0]);
    let u = dot3(t_vec, p) * inv;
    if !(-EPS..=1.0 + EPS).contains(&u) {
        return None;
    }
    let q = cross3(t_vec, e1);
    let v = dot3(dir, q) * inv;
    if v < -EPS || u + v > 1.0 + EPS {
        return None;
    }
    let t = dot3(e2, q) * inv;
    (-EPS..=1.0 + EPS).contains(&t).then_some(t)
}

/// Axis-aligned box as a closed, outward-wound triangle mesh.
///
/// Originally `inside_3d`'s test helper; promoted here because `intersect_3d`
/// emits one as an intersection's bounding solid at run time, and every 3D
/// test in the crate builds its fixtures from it.
pub(crate) fn box_mesh(min: [f64; 3], max: [f64; 3]) -> Geometry {
    let (x0, y0, z0) = (min[0], min[1], min[2]);
    let (x1, y1, z1) = (max[0], max[1], max[2]);
    let v = [
        [x0, y0, z0],
        [x1, y0, z0],
        [x1, y1, z0],
        [x0, y1, z0],
        [x0, y0, z1],
        [x1, y0, z1],
        [x1, y1, z1],
        [x0, y1, z1],
    ];
    // Outward-facing triangles for all six faces.
    let faces: [[usize; 3]; 12] = [
        [0, 2, 1],
        [0, 3, 2], // bottom
        [4, 5, 6],
        [4, 6, 7], // top
        [0, 1, 5],
        [0, 5, 4], // front
        [1, 2, 6],
        [1, 6, 5], // right
        [2, 3, 7],
        [2, 7, 6], // back
        [3, 0, 4],
        [3, 4, 7], // left
    ];
    Geometry::MultiPolygon(
        faces
            .iter()
            .map(|f| {
                (
                    Ring::new(
                        f.iter()
                            .map(|i| Coord::xyz(v[*i][0], v[*i][1], v[*i][2]))
                            .collect::<Vec<_>>(),
                    ),
                    Vec::new(),
                )
            })
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    fn box_tris(min: [f64; 3], max: [f64; 3]) -> Vec<Tri> {
        crate::inside_3d::collect_triangles(&box_mesh(min, max))
    }

    #[test]
    fn a_closed_box_reports_closed_and_consistently_wound() {
        let t = topology(&box_tris([0.0; 3], [1.0, 1.0, 1.0]));
        assert!(t.closed);
        assert_eq!(t.open_edges, 0);
        assert_eq!(t.nonmanifold_edges, 0);
        assert!(t.consistent_winding);
    }

    #[test]
    fn removing_a_face_opens_the_mesh_and_exposes_a_loop() {
        let mut tris = box_tris([0.0; 3], [1.0, 1.0, 1.0]);
        // Drop the two triangles forming the bottom face.
        tris.drain(0..2);
        let t = topology(&tris);
        assert!(!t.closed);
        // A square hole has four boundary edges.
        assert_eq!(t.open_edges, 4);

        let loops = boundary_loops(&tris);
        assert_eq!(loops.len(), 1, "expected one boundary loop");
        assert_eq!(loops[0].len(), 4);
    }

    #[test]
    fn capping_the_boundary_loop_restores_watertightness_and_volume() {
        // The property enclose_multipatch depends on.
        let full = box_tris([0.0; 3], [2.0, 3.0, 4.0]);
        let expected = mesh_volume(&full);

        let mut open = full.clone();
        open.drain(0..2);
        assert!(!topology(&open).closed);

        // `reverse` is not cosmetic: boundary_loops reports each edge in its
        // owning triangle's direction, so the cap must traverse it the other
        // way. Capping the same way round yields a watertight but
        // inconsistently oriented mesh whose signed volume cancels to zero.
        for ring in boundary_loops(&open) {
            open.extend(fan_triangulate(&ring, true));
        }
        let capped = topology(&open);
        assert!(capped.closed, "cap left {} open edges", capped.open_edges);
        assert!(
            capped.consistent_winding,
            "cap closed the mesh but flipped its orientation"
        );
        assert!((mesh_volume(&open) - expected).abs() < 1e-9);

        // And the wrong winding really does break it, so the check above is
        // not vacuous.
        let mut wrong = full.clone();
        wrong.drain(0..2);
        for ring in boundary_loops(&wrong) {
            wrong.extend(fan_triangulate(&ring, false));
        }
        assert!(!topology(&wrong).consistent_winding);
    }

    #[test]
    fn mesh_volume_and_area_match_a_known_box() {
        let tris = box_tris([0.0; 3], [2.0, 3.0, 4.0]);
        assert!((mesh_volume(&tris) - 24.0).abs() < 1e-9);
        // Surface area 2*(2*3 + 3*4 + 2*4) = 52.
        assert!((mesh_area(&tris) - 52.0).abs() < 1e-9);
    }

    #[test]
    fn segment_triangle_hits_only_within_the_segment() {
        let tri: Tri = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        // Straight through the middle of the triangle.
        assert!(segment_triangle([0.25, 0.25, -1.0], [0.25, 0.25, 1.0], &tri).is_some());
        // Same line, but the segment stops short of the plane.
        assert!(segment_triangle([0.25, 0.25, -1.0], [0.25, 0.25, -0.5], &tri).is_none());
        // Passes outside the triangle's extent.
        assert!(segment_triangle([5.0, 5.0, -1.0], [5.0, 5.0, 1.0], &tri).is_none());
    }

    #[test]
    fn triangles_survive_a_geometry_round_trip() {
        let tris = box_tris([0.0; 3], [1.0, 1.0, 1.0]);
        let geom = triangles_to_geometry(&tris);
        let back = crate::inside_3d::collect_triangles(&geom);
        assert_eq!(back.len(), tris.len());
        assert!((mesh_volume(&back) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn occupancy_recovers_a_box_volume() {
        let tris = box_tris([0.0; 3], [4.0, 4.0, 4.0]);
        let solid = Solid::new(0, tris);
        let v = occupancy_volume([0.0; 3], [4.0, 4.0, 4.0], 32, |x, y, z| {
            solid.contains(x, y, z)
        });
        assert!((v - 64.0).abs() < 1.0, "occupancy gave {v}, expected ~64");
    }
}
