//! Raster-to-vector polygonization, a pure-Rust port of the `gdal.Polygonize`
//! step in `lidar` (`filling.py` / `slicing.py`).
//!
//! Each connected group of equal-valued, nonzero cells (4-connectivity, as GDAL
//! uses by default) becomes one polygon feature. Cell-edge boundaries are traced
//! into rings, inner rings are emitted as holes, and the per-region attributes
//! are joined onto each feature. Output is a GeoJSON `FeatureCollection` in the
//! raster's own CRS, which is what GeoLibre / MapLibre consume directly.

use std::collections::HashMap;

use serde_json::{json, Map, Value};

/// Geometry/attribute inputs for one polygonization pass.
pub struct PolygonizeParams<'a> {
    pub labels: &'a [f64],
    pub rows: usize,
    pub cols: usize,
    pub x_min: f64,
    pub y_max: f64,
    pub cell_size_x: f64,
    pub cell_size_y: f64,
    pub epsg: Option<u32>,
    /// Attribute table keyed by integer id; merged into each feature's
    /// properties (alongside the always-present `id`).
    pub props_by_id: &'a HashMap<i64, Map<String, Value>>,
}

/// A vertex on the cell-corner grid: `(col, row)`, with `0..=cols` / `0..=rows`.
type Vertex = (i64, i64);

/// Polygonizes a label raster into a GeoJSON `FeatureCollection` string.
pub fn polygonize_to_geojson(p: &PolygonizeParams) -> String {
    let groups = connected_groups(p.labels, p.rows, p.cols);

    let mut features: Vec<Value> = Vec::new();
    for group in groups {
        let rings = trace_rings(&group.cells, p.cols);
        if rings.is_empty() {
            continue;
        }
        let geometry = rings_to_geometry(rings, p);

        let mut props = p
            .props_by_id
            .get(&group.id)
            .cloned()
            .unwrap_or_default();
        props.entry("id".to_string()).or_insert(json!(group.id));

        features.push(json!({
            "type": "Feature",
            "properties": Value::Object(props),
            "geometry": geometry,
        }));
    }

    let mut fc = Map::new();
    fc.insert("type".to_string(), json!("FeatureCollection"));
    if let Some(epsg) = p.epsg {
        fc.insert(
            "crs".to_string(),
            json!({"type": "name", "properties": {"name": format!("EPSG:{epsg}")}}),
        );
    }
    fc.insert("features".to_string(), Value::Array(features));
    Value::Object(fc).to_string()
}

/// A connected component of equal-valued, nonzero cells.
struct Group {
    id: i64,
    cells: Vec<usize>,
}

/// Flood-fills the raster into 4-connected groups of equal nonzero value
/// (matching GDAL's per-value connected polygons).
fn connected_groups(labels: &[f64], rows: usize, cols: usize) -> Vec<Group> {
    let mut visited = vec![false; rows * cols];
    let mut groups = Vec::new();
    let mut stack: Vec<usize> = Vec::new();

    for start in 0..rows * cols {
        if visited[start] || labels[start] == 0.0 {
            continue;
        }
        let id = labels[start] as i64;
        let mut cells = Vec::new();
        stack.push(start);
        visited[start] = true;
        while let Some(i) = stack.pop() {
            cells.push(i);
            let r = i / cols;
            let c = i % cols;
            for (dr, dc) in [(-1isize, 0isize), (1, 0), (0, -1), (0, 1)] {
                let nr = r as isize + dr;
                let nc = c as isize + dc;
                if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                    continue;
                }
                let ni = nr as usize * cols + nc as usize;
                if !visited[ni] && labels[ni] as i64 == id && labels[ni] != 0.0 {
                    visited[ni] = true;
                    stack.push(ni);
                }
            }
        }
        groups.push(Group { id, cells });
    }
    groups
}

/// Traces the boundary of a cell group into one or more closed rings of grid
/// vertices. Interior cell edges cancel; only the outer/inner boundaries remain.
fn trace_rings(cells: &[usize], cols: usize) -> Vec<Vec<Vertex>> {
    let cellset: std::collections::HashSet<usize> = cells.iter().copied().collect();
    let in_group = |r: i64, c: i64| -> bool {
        if r < 0 || c < 0 {
            return false;
        }
        cellset.contains(&(r as usize * cols + c as usize))
    };

    // Directed boundary edges (interior on a consistent side); shared edges
    // between two in-group cells appear in both directions and are omitted.
    let mut outgoing: HashMap<Vertex, Vec<Vertex>> = HashMap::new();
    let mut push_edge = |a: Vertex, b: Vertex| outgoing.entry(a).or_default().push(b);

    for &i in cells {
        let r = (i / cols) as i64;
        let c = (i % cols) as i64;
        let a = (c, r);
        let b = (c + 1, r);
        let d = (c, r + 1);
        let e = (c + 1, r + 1);
        if !in_group(r - 1, c) {
            push_edge(a, b); // top
        }
        if !in_group(r, c + 1) {
            push_edge(b, e); // right
        }
        if !in_group(r + 1, c) {
            push_edge(e, d); // bottom
        }
        if !in_group(r, c - 1) {
            push_edge(d, a); // left
        }
    }

    // Walk directed edges into closed rings.
    let mut rings = Vec::new();
    let mut vertices: Vec<Vertex> = outgoing.keys().copied().collect();
    vertices.sort();
    for start in vertices {
        while outgoing.get(&start).map(|v| !v.is_empty()).unwrap_or(false) {
            let mut ring = vec![start];
            let mut cur = start;
            let mut prev_dir: Option<(i64, i64)> = None;
            loop {
                let next = {
                    let outs = outgoing.get_mut(&cur).expect("vertex has edges");
                    let idx = choose_next(outs, cur, prev_dir);
                    outs.swap_remove(idx)
                };
                prev_dir = Some((next.0 - cur.0, next.1 - cur.1));
                cur = next;
                if cur == start {
                    break;
                }
                ring.push(cur);
            }
            if ring.len() >= 3 {
                rings.push(simplify_collinear(ring));
            }
        }
    }
    rings
}

/// Picks the next edge when walking a ring. At a simple (degree-2) vertex there
/// is one choice; at a checkerboard pinch we prefer the left-most turn so the
/// two loops separate cleanly.
fn choose_next(outs: &[Vertex], cur: Vertex, prev_dir: Option<(i64, i64)>) -> usize {
    if outs.len() == 1 || prev_dir.is_none() {
        return 0;
    }
    let (dx, dy) = prev_dir.unwrap();
    let mut best = 0;
    let mut best_score = i64::MAX;
    for (k, &target) in outs.iter().enumerate() {
        // Edges are emitted so the region's interior is on the right of each
        // directed edge. At a checkerboard pinch (a degree-4 vertex), taking the
        // most-clockwise turn (minimum cross product) hugs the interior and
        // splits the boundary into two non-crossing loops.
        let (ex, ey) = (target.0 - cur.0, target.1 - cur.1);
        let cross = dx * ey - dy * ex;
        if cross < best_score {
            best_score = cross;
            best = k;
        }
    }
    best
}

/// Removes vertices that are collinear with their neighbors (rings are
/// rectilinear, so this merges runs along the same axis).
fn simplify_collinear(ring: Vec<Vertex>) -> Vec<Vertex> {
    let n = ring.len();
    if n < 3 {
        return ring;
    }
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let prev = ring[(i + n - 1) % n];
        let cur = ring[i];
        let next = ring[(i + 1) % n];
        let d1 = (cur.0 - prev.0, cur.1 - prev.1);
        let d2 = (next.0 - cur.0, next.1 - cur.1);
        // Keep the vertex only if the direction changes (cross product != 0).
        if d1.0 * d2.1 - d1.1 * d2.0 != 0 {
            out.push(cur);
        }
    }
    out
}

/// Converts traced rings to a GeoJSON `Polygon`, classifying the largest ring as
/// the exterior and the rest as holes, with RFC 7946 winding (exterior CCW,
/// holes CW).
fn rings_to_geometry(rings: Vec<Vec<Vertex>>, p: &PolygonizeParams) -> Value {
    let to_world = |v: Vertex| -> [f64; 2] {
        [
            p.x_min + v.0 as f64 * p.cell_size_x,
            p.y_max - v.1 as f64 * p.cell_size_y,
        ]
    };

    // World-space rings with their signed area.
    let mut world: Vec<(Vec<[f64; 2]>, f64)> = rings
        .into_iter()
        .map(|ring| {
            let coords: Vec<[f64; 2]> = ring.iter().map(|&v| to_world(v)).collect();
            let area = signed_area(&coords);
            (coords, area)
        })
        .collect();

    // Exterior = largest absolute area.
    let exterior_idx = world
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.1.abs().total_cmp(&b.1.abs()))
        .map(|(i, _)| i)
        .unwrap_or(0);

    let mut coordinates: Vec<Value> = Vec::with_capacity(world.len());
    // Emit exterior first (CCW), then holes (CW).
    let exterior = std::mem::take(&mut world[exterior_idx]);
    coordinates.push(close_ring(orient(exterior.0, exterior.1, true)));
    for (i, (coords, area)) in world.into_iter().enumerate() {
        if i == exterior_idx || coords.is_empty() {
            continue;
        }
        coordinates.push(close_ring(orient(coords, area, false)));
    }

    json!({"type": "Polygon", "coordinates": coordinates})
}

/// Shoelace signed area (positive = counterclockwise).
fn signed_area(ring: &[[f64; 2]]) -> f64 {
    let n = ring.len();
    let mut a = 0.0;
    for i in 0..n {
        let p1 = ring[i];
        let p2 = ring[(i + 1) % n];
        a += p1[0] * p2[1] - p2[0] * p1[1];
    }
    a / 2.0
}

/// Orients a ring CCW (`want_ccw = true`) or CW, given its current signed area.
fn orient(mut ring: Vec<[f64; 2]>, area: f64, want_ccw: bool) -> Vec<[f64; 2]> {
    let is_ccw = area > 0.0;
    if is_ccw != want_ccw {
        ring.reverse();
    }
    ring
}

/// Closes a ring (repeats the first point) and converts to JSON.
fn close_ring(ring: Vec<[f64; 2]>) -> Value {
    let mut pts: Vec<Value> = ring.iter().map(|p| json!([p[0], p[1]])).collect();
    if let Some(first) = ring.first() {
        pts.push(json!([first[0], first[1]]));
    }
    Value::Array(pts)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(s: &str) -> Value {
        serde_json::from_str(s).unwrap()
    }

    #[test]
    fn square_block_becomes_one_polygon() {
        // 4x4 grid, a 2x2 block of id 1 in the top-left.
        let rows = 4;
        let cols = 4;
        let mut labels = vec![0.0; rows * cols];
        for (r, c) in [(0, 0), (0, 1), (1, 0), (1, 1)] {
            labels[r * cols + c] = 1.0;
        }
        let props: HashMap<i64, Map<String, Value>> = HashMap::new();
        let geojson = polygonize_to_geojson(&PolygonizeParams {
            labels: &labels,
            rows,
            cols,
            x_min: 0.0,
            y_max: 4.0,
            cell_size_x: 1.0,
            cell_size_y: 1.0,
            epsg: Some(4326),
            props_by_id: &props,
        });
        let v = parse(&geojson);
        let features = v["features"].as_array().unwrap();
        assert_eq!(features.len(), 1);
        let coords = features[0]["geometry"]["coordinates"].as_array().unwrap();
        // One exterior ring, no holes.
        assert_eq!(coords.len(), 1);
        let ring = coords[0].as_array().unwrap();
        // 4 corners after collinear merge, plus the closing point.
        assert_eq!(ring.len(), 5);
        assert_eq!(features[0]["properties"]["id"], json!(1));
    }

    #[test]
    fn donut_has_a_hole() {
        // 5x5 ring of id 1 around a hole at the center.
        let rows = 5;
        let cols = 5;
        let mut labels = vec![0.0; rows * cols];
        for r in 1..4 {
            for c in 1..4 {
                labels[r * cols + c] = 1.0;
            }
        }
        labels[2 * cols + 2] = 0.0; // punch the hole
        let props: HashMap<i64, Map<String, Value>> = HashMap::new();
        let geojson = polygonize_to_geojson(&PolygonizeParams {
            labels: &labels,
            rows,
            cols,
            x_min: 0.0,
            y_max: 5.0,
            cell_size_x: 1.0,
            cell_size_y: 1.0,
            epsg: None,
            props_by_id: &props,
        });
        let v = parse(&geojson);
        let coords = v["features"][0]["geometry"]["coordinates"]
            .as_array()
            .unwrap();
        // Exterior + one hole.
        assert_eq!(coords.len(), 2);
    }

    #[test]
    fn winding_is_rfc7946() {
        let rows = 3;
        let cols = 3;
        let mut labels = vec![0.0; rows * cols];
        labels[0] = 1.0;
        let props: HashMap<i64, Map<String, Value>> = HashMap::new();
        let geojson = polygonize_to_geojson(&PolygonizeParams {
            labels: &labels,
            rows,
            cols,
            x_min: 0.0,
            y_max: 3.0,
            cell_size_x: 1.0,
            cell_size_y: 1.0,
            epsg: None,
            props_by_id: &props,
        });
        let v = parse(&geojson);
        let ring = v["features"][0]["geometry"]["coordinates"][0]
            .as_array()
            .unwrap();
        let pts: Vec<[f64; 2]> = ring
            .iter()
            .map(|p| [p[0].as_f64().unwrap(), p[1].as_f64().unwrap()])
            .collect();
        // Exterior ring must be counterclockwise (positive signed area).
        assert!(signed_area(&pts[..pts.len() - 1]) > 0.0);
    }
}
