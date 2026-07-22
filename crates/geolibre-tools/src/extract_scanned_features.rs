//! GeoLibre tool: vectorize scanned / binarized map imagery into features.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Extract Scanned Lines* and
//! *Extract Scanned Polygons* (Conversion, the ArcScan successor): turn a binary
//! raster (a scanned or thresholded map) into clean vector **lines** or
//! **polygons**, the natural front end to the repo's raster-to-vector cleanup
//! pipeline (`regularize_building_footprints`, `smooth_natural_features`, …).
//!
//! The bundled `raster_to_vector_lines` traces raw cell edges only. This tool
//! does the full scanned-map workflow, all classical image processing (no deep
//! learning):
//!
//! 1. **Binarize** — cells equal to `foreground_value` (or ≥ `threshold`, or
//!    simply non-zero/valid) are foreground.
//! 2. **Denoise** — drop foreground components smaller than `noise_size` cells
//!    and fill background holes smaller than `hole_size` cells.
//! 3. **Lines**: **Zhang–Suen thinning** to a 1-px skeleton, then trace the
//!    skeleton into polylines that break only at junctions (degree ≠ 2), close
//!    small gaps between nearly-collinear endpoints within `gap_distance`, and
//!    simplify with Douglas–Peucker (`simplify_tolerance`).
//! 4. **Polygons**: vectorize the cleaned foreground with the repo's own
//!    `polygonize` (cell-edge rings, holes preserved) and simplify.
//!
//! Output carries the raster's CRS, ready for the downstream regularize/smooth
//! tools. `max_line_width` filtering is implicit in thinning (everything reduces
//! to a 1-px centreline); an explicit width cut is future work.

use std::collections::{BTreeMap, HashMap, HashSet};

use geo::{Coord as GeoCoord, LineString, Simplify};
use serde_json::{json, Map, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{
    Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring,
};

use crate::common::{band_to_vec, load_input_raster};
use crate::polygonize::{polygonize_to_geojson, PolygonizeParams};
use crate::vector_common::{parse_optional_str, write_or_store_layer};

pub struct ExtractScannedFeaturesTool;

impl Tool for ExtractScannedFeaturesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "extract_scanned_features",
            display_name: "Extract Scanned Features",
            summary: "Vectorize a binary (scanned/thresholded) raster into clean lines (skeleton tracing) or polygons, with noise removal, hole filling, gap closure, and simplification.",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input binary/classified raster to vectorize.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "feature_type",
                    description: "'lines' (default): skeletonize and trace centrelines. 'polygons': vectorize filled regions.",
                    required: false,
                },
                ToolParamSpec {
                    name: "foreground_value",
                    description: "Cell value that marks foreground. If omitted, cells >= 'threshold' (or, absent that, any non-zero valid cell) are foreground.",
                    required: false,
                },
                ToolParamSpec {
                    name: "threshold",
                    description: "Foreground is cells with value >= this (used when 'foreground_value' is not given).",
                    required: false,
                },
                ToolParamSpec {
                    name: "noise_size",
                    description: "Remove foreground components smaller than this many cells (default 4).",
                    required: false,
                },
                ToolParamSpec {
                    name: "hole_size",
                    description: "Fill background holes smaller than this many cells (default 4).",
                    required: false,
                },
                ToolParamSpec {
                    name: "gap_distance",
                    description: "For lines: bridge gaps up to this distance (CRS units) between nearly-collinear line endpoints. Default 0 (off).",
                    required: false,
                },
                ToolParamSpec {
                    name: "simplify_tolerance",
                    description: "Douglas-Peucker tolerance (CRS units) for the output geometry. Default ~1 cell.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read (default 1).",
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

        let raster = load_input_raster(input)?;
        if prm.band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "band out of range (raster has {} band(s))",
                raster.bands
            )));
        }
        let rows = raster.rows;
        let cols = raster.cols;
        let nodata = raster.nodata;
        let data = band_to_vec(&raster, prm.band);
        let csx = raster.cell_size_x;
        let csy = raster.cell_size_y;
        let x_min = raster.x_min;
        let y_max = raster.y_min + rows as f64 * csy;
        let simplify_tol = prm.simplify_tolerance.unwrap_or(csx.max(csy));

        // 1. Binarize.
        let mut fg = vec![false; rows * cols];
        for (i, &v) in data.iter().enumerate() {
            if v == nodata || !v.is_finite() {
                continue;
            }
            fg[i] = match (prm.foreground_value, prm.threshold) {
                (Some(fv), _) => (v - fv).abs() < 0.5,
                (None, Some(t)) => v >= t,
                (None, None) => v != 0.0,
            };
        }
        let fg_before = fg.iter().filter(|&&b| b).count();

        // 2. Denoise: drop tiny foreground blobs, fill tiny background holes.
        remove_small_components(&mut fg, rows, cols, prm.noise_size, true);
        fill_small_holes(&mut fg, rows, cols, prm.hole_size);
        let fg_after = fg.iter().filter(|&&b| b).count();
        ctx.progress.info(&format!(
            "foreground {fg_before} -> {fg_after} cell(s) after denoise"
        ));

        let epsg = raster.crs.epsg;
        let px = |c: i64| x_min + (c as f64 + 0.5) * csx;
        let py = |r: i64| y_max - (r as f64 + 0.5) * csy;

        let (geoms, kind) = match prm.feature_type {
            FeatureType::Lines => {
                let mut skel = fg.clone();
                zhang_suen_thin(&mut skel, rows, cols);
                let skel_px = skel.iter().filter(|&&b| b).count();
                ctx.progress
                    .info(&format!("skeleton: {skel_px} centreline cell(s)"));
                let mut lines = trace_skeleton(&skel, rows, cols, &px, &py);
                if prm.gap_distance > 0.0 {
                    close_gaps(&mut lines, prm.gap_distance);
                }
                let geoms = lines
                    .into_iter()
                    .filter(|l| l.len() >= 2)
                    .map(|l| simplify_line(l, simplify_tol))
                    .collect::<Vec<_>>();
                (geoms, GeometryType::LineString)
            }
            FeatureType::Polygons => {
                let labels: Vec<f64> = fg.iter().map(|&b| if b { 1.0 } else { 0.0 }).collect();
                let props: HashMap<i64, Map<String, Value>> = HashMap::new();
                let geojson = polygonize_to_geojson(&PolygonizeParams {
                    labels: &labels,
                    rows,
                    cols,
                    x_min,
                    y_max,
                    cell_size_x: csx,
                    cell_size_y: csy,
                    epsg,
                    props_by_id: &props,
                });
                let geoms = parse_polygons(&geojson)?
                    .into_iter()
                    .map(|g| simplify_polygon(g, simplify_tol))
                    .collect::<Vec<_>>();
                (geoms, GeometryType::Polygon)
            }
        };

        // Build the output layer.
        let mut layer = Layer::new("scanned_features");
        layer.geom_type = Some(kind);
        if let Some(e) = epsg {
            layer = layer.with_crs_epsg(e);
        }
        layer.add_field(FieldDef::new("id", FieldType::Integer));
        layer.add_field(FieldDef::new("measure", FieldType::Float)); // length or area
        let mut total_measure = 0.0;
        for (i, g) in geoms.into_iter().enumerate() {
            let m = measure_of(&g);
            total_measure += m;
            let mut f = Feature::with_geometry(i as u64, g, 2);
            f.set_by_index(0, FieldValue::Integer(i as i64));
            f.set_by_index(1, FieldValue::Float(m));
            layer.push(f);
        }

        let feature_count = layer.len();
        ctx.progress
            .info(&format!("extracted {feature_count} feature(s)"));
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("feature_type".to_string(), json!(prm.feature_type.as_str()));
        outputs.insert("foreground_cells".to_string(), json!(fg_after));
        outputs.insert("total_measure".to_string(), json!(total_measure));
        Ok(ToolRunResult { outputs })
    }
}

// ── Connected-component denoise ──────────────────────────────────────────────

/// Removes 8-connected `target` components smaller than `min_size` cells.
fn remove_small_components(
    grid: &mut [bool],
    rows: usize,
    cols: usize,
    min_size: usize,
    target: bool,
) {
    if min_size <= 1 {
        return;
    }
    let mut visited = vec![false; rows * cols];
    let mut stack = Vec::new();
    for start in 0..rows * cols {
        if visited[start] || grid[start] != target {
            continue;
        }
        let mut comp = Vec::new();
        stack.push(start);
        visited[start] = true;
        while let Some(i) = stack.pop() {
            comp.push(i);
            let (r, c) = ((i / cols) as i64, (i % cols) as i64);
            for (dr, dc) in NEIGH8 {
                let (nr, nc) = (r + dr, c + dc);
                if nr < 0 || nc < 0 || nr >= rows as i64 || nc >= cols as i64 {
                    continue;
                }
                let ni = nr as usize * cols + nc as usize;
                if !visited[ni] && grid[ni] == target {
                    visited[ni] = true;
                    stack.push(ni);
                }
            }
        }
        if comp.len() < min_size {
            for i in comp {
                grid[i] = !target;
            }
        }
    }
}

/// Fills background holes (4-connected `false` regions not touching the border)
/// smaller than `max_size` cells.
fn fill_small_holes(grid: &mut [bool], rows: usize, cols: usize, max_size: usize) {
    if max_size == 0 {
        return;
    }
    let mut visited = vec![false; rows * cols];
    let mut stack = Vec::new();
    for start in 0..rows * cols {
        if visited[start] || grid[start] {
            continue;
        }
        let mut comp = Vec::new();
        let mut touches_border = false;
        stack.push(start);
        visited[start] = true;
        while let Some(i) = stack.pop() {
            comp.push(i);
            let (r, c) = ((i / cols) as i64, (i % cols) as i64);
            if r == 0 || c == 0 || r == rows as i64 - 1 || c == cols as i64 - 1 {
                touches_border = true;
            }
            for (dr, dc) in NEIGH4 {
                let (nr, nc) = (r + dr, c + dc);
                if nr < 0 || nc < 0 || nr >= rows as i64 || nc >= cols as i64 {
                    continue;
                }
                let ni = nr as usize * cols + nc as usize;
                if !visited[ni] && !grid[ni] {
                    visited[ni] = true;
                    stack.push(ni);
                }
            }
        }
        if !touches_border && comp.len() <= max_size {
            for i in comp {
                grid[i] = true;
            }
        }
    }
}

// ── Zhang–Suen thinning ──────────────────────────────────────────────────────

/// Reduces the foreground to a 1-pixel-wide skeleton in place (Zhang & Suen
/// 1984). Two sub-iterations per pass until no pixel changes.
fn zhang_suen_thin(grid: &mut [bool], rows: usize, cols: usize) {
    let at = |g: &[bool], r: i64, c: i64| -> bool {
        r >= 0 && c >= 0 && r < rows as i64 && c < cols as i64 && g[r as usize * cols + c as usize]
    };
    loop {
        let mut changed = false;
        for sub in 0..2 {
            let mut to_clear = Vec::new();
            for r in 0..rows as i64 {
                for c in 0..cols as i64 {
                    if !at(grid, r, c) {
                        continue;
                    }
                    // 8 neighbours P2..P9 clockwise from north.
                    let p = [
                        at(grid, r - 1, c),     // P2 N
                        at(grid, r - 1, c + 1), // P3 NE
                        at(grid, r, c + 1),     // P4 E
                        at(grid, r + 1, c + 1), // P5 SE
                        at(grid, r + 1, c),     // P6 S
                        at(grid, r + 1, c - 1), // P7 SW
                        at(grid, r, c - 1),     // P8 W
                        at(grid, r - 1, c - 1), // P9 NW
                    ];
                    let bp: usize = p.iter().filter(|&&x| x).count();
                    if !(2..=6).contains(&bp) {
                        continue;
                    }
                    // A(P1): number of 0->1 transitions in P2..P9,P2.
                    let mut a = 0;
                    for k in 0..8 {
                        if !p[k] && p[(k + 1) % 8] {
                            a += 1;
                        }
                    }
                    if a != 1 {
                        continue;
                    }
                    let (c1, c2) = if sub == 0 {
                        // P2*P4*P6 and P4*P6*P8
                        (p[0] && p[2] && p[4], p[2] && p[4] && p[6])
                    } else {
                        // P2*P4*P8 and P2*P6*P8
                        (p[0] && p[2] && p[6], p[0] && p[4] && p[6])
                    };
                    if c1 || c2 {
                        continue;
                    }
                    to_clear.push(r as usize * cols + c as usize);
                }
            }
            if !to_clear.is_empty() {
                changed = true;
                for i in to_clear {
                    grid[i] = false;
                }
            }
        }
        if !changed {
            break;
        }
    }
}

// ── Skeleton tracing ─────────────────────────────────────────────────────────

/// Traces a 1-px skeleton into polylines that break at junctions (8-neighbour
/// degree ≠ 2). Isolated loops are broken at an arbitrary pixel.
fn trace_skeleton(
    skel: &[bool],
    rows: usize,
    cols: usize,
    px: &impl Fn(i64) -> f64,
    py: &impl Fn(i64) -> f64,
) -> Vec<Vec<Coord>> {
    let idx = |r: i64, c: i64| r as usize * cols + c as usize;
    let inb = |r: i64, c: i64| r >= 0 && c >= 0 && r < rows as i64 && c < cols as i64;
    let neighbours = |r: i64, c: i64| -> Vec<(i64, i64)> {
        let mut v = Vec::new();
        for (dr, dc) in NEIGH8 {
            let (nr, nc) = (r + dr, c + dc);
            if inb(nr, nc) && skel[idx(nr, nc)] {
                v.push((nr, nc));
            }
        }
        v
    };
    let degree = |r: i64, c: i64| neighbours(r, c).len();

    let mut used_edges: HashSet<(usize, usize)> = HashSet::new();
    let edge = |a: usize, b: usize| -> (usize, usize) {
        if a < b {
            (a, b)
        } else {
            (b, a)
        }
    };
    let mut lines = Vec::new();

    // Trace from every node (degree != 2) along each unused incident edge.
    let world = |r: i64, c: i64| Coord::xy(px(c), py(r));
    for r in 0..rows as i64 {
        for c in 0..cols as i64 {
            if !skel[idx(r, c)] || degree(r, c) == 2 {
                continue;
            }
            for (nr, nc) in neighbours(r, c) {
                let e = edge(idx(r, c), idx(nr, nc));
                if used_edges.contains(&e) {
                    continue;
                }
                // Walk the path until the next node.
                let mut path = vec![world(r, c)];
                used_edges.insert(e);
                let (mut pr, mut pc) = (r, c);
                let (mut cr, mut cc) = (nr, nc);
                loop {
                    path.push(world(cr, cc));
                    if degree(cr, cc) != 2 {
                        break;
                    }
                    // Continue to the neighbour that is not where we came from.
                    let next = neighbours(cr, cc)
                        .into_iter()
                        .find(|&(a, b)| !(a == pr && b == pc));
                    match next {
                        Some((a, b)) => {
                            let e2 = edge(idx(cr, cc), idx(a, b));
                            if used_edges.contains(&e2) {
                                break;
                            }
                            used_edges.insert(e2);
                            (pr, pc) = (cr, cc);
                            (cr, cc) = (a, b);
                        }
                        None => break,
                    }
                }
                if path.len() >= 2 {
                    lines.push(path);
                }
            }
        }
    }

    // Remaining pixels form pure loops (all degree 2): break each at a pixel.
    let mut seen: HashSet<usize> = HashSet::new();
    for r in 0..rows as i64 {
        for c in 0..cols as i64 {
            let i = idx(r, c);
            if !skel[i] || degree(r, c) != 2 || seen.contains(&i) {
                continue;
            }
            // Only start a loop if none of its edges were traced.
            if neighbours(r, c)
                .iter()
                .all(|&(a, b)| !used_edges.contains(&edge(i, idx(a, b))))
            {
                let mut path = vec![world(r, c)];
                seen.insert(i);
                let (mut pr, mut pc) = (r, c);
                let start = (r, c);
                let (mut cr, mut cc) = neighbours(r, c)[0];
                loop {
                    let e = edge(idx(pr, pc), idx(cr, cc));
                    if used_edges.contains(&e) {
                        break;
                    }
                    used_edges.insert(e);
                    path.push(world(cr, cc));
                    seen.insert(idx(cr, cc));
                    if (cr, cc) == start {
                        break;
                    }
                    let next = neighbours(cr, cc)
                        .into_iter()
                        .find(|&(a, b)| !(a == pr && b == pc));
                    match next {
                        Some(n) => {
                            (pr, pc) = (cr, cc);
                            (cr, cc) = n;
                        }
                        None => break,
                    }
                }
                if path.len() >= 2 {
                    lines.push(path);
                }
            }
        }
    }
    lines
}

/// Bridges gaps between nearly-collinear endpoints of distinct lines within
/// `max_gap`, greedily connecting each free endpoint to its best partner.
fn close_gaps(lines: &mut Vec<Vec<Coord>>, max_gap: f64) {
    // Endpoints: (line index, is_tail).
    let n = lines.len();
    let mut joined = vec![false; n];
    // Greedy: for each line's endpoints, find the nearest compatible endpoint.
    for i in 0..n {
        if joined[i] || lines[i].len() < 2 {
            continue;
        }
        for tail_i in [true, false] {
            let ep = endpoint(&lines[i], tail_i);
            let dir_i = direction(&lines[i], tail_i);
            let mut best: Option<(usize, bool, f64)> = None;
            for (j, lj) in lines.iter().enumerate() {
                if j == i || joined[j] || lj.len() < 2 {
                    continue;
                }
                for tail_j in [true, false] {
                    let q = endpoint(lj, tail_j);
                    let d = (ep.x - q.x).hypot(ep.y - q.y);
                    if d > max_gap || d == 0.0 {
                        continue;
                    }
                    // Prefer roughly collinear continuations (angle < ~60deg).
                    let dir_j = direction(lj, tail_j);
                    let dot = dir_i.0 * (-dir_j.0) + dir_i.1 * (-dir_j.1);
                    if dot < 0.5 {
                        continue;
                    }
                    if best.map(|b| d < b.2).unwrap_or(true) {
                        best = Some((j, tail_j, d));
                    }
                }
            }
            if let Some((j, tail_j, _)) = best {
                // Append line j onto line i, oriented to continue.
                let mut other = lines[j].clone();
                if tail_j {
                    other.reverse();
                }
                if tail_i {
                    lines[i].extend(other);
                } else {
                    let mut new = other;
                    new.reverse();
                    new.extend(lines[i].clone());
                    lines[i] = new;
                }
                joined[j] = true;
            }
        }
    }
    let mut out = Vec::new();
    for (i, l) in lines.drain(..).enumerate() {
        if !joined[i] {
            out.push(l);
        }
    }
    *lines = out;
}

fn endpoint(line: &[Coord], tail: bool) -> Coord {
    if tail {
        line.last().unwrap().clone()
    } else {
        line.first().unwrap().clone()
    }
}

/// Unit direction pointing outward at the given endpoint.
fn direction(line: &[Coord], tail: bool) -> (f64, f64) {
    let (a, b) = if tail {
        (&line[line.len() - 2], &line[line.len() - 1])
    } else {
        (&line[1], &line[0])
    };
    let (dx, dy) = (b.x - a.x, b.y - a.y);
    let l = dx.hypot(dy).max(f64::MIN_POSITIVE);
    (dx / l, dy / l)
}

// ── Simplify + measures ──────────────────────────────────────────────────────

fn simplify_line(coords: Vec<Coord>, tol: f64) -> Geometry {
    let ls = LineString::new(coords.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect());
    let s = ls.simplify(tol);
    Geometry::LineString(s.0.iter().map(|c| Coord::xy(c.x, c.y)).collect())
}

fn simplify_polygon(geom: Geometry, tol: f64) -> Geometry {
    let Geometry::Polygon {
        exterior,
        interiors,
    } = &geom
    else {
        return geom;
    };
    let simp = |r: &Ring| -> Ring {
        let ls = LineString::new(
            r.coords()
                .iter()
                .map(|c| GeoCoord { x: c.x, y: c.y })
                .collect(),
        );
        let s = ls.simplify(tol);
        let mut cs: Vec<Coord> = s.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
        if cs.len() >= 2 && cs.first() == cs.last() {
            cs.pop();
        }
        Ring::new(cs)
    };
    Geometry::Polygon {
        exterior: simp(exterior),
        interiors: interiors.iter().map(simp).collect(),
    }
}

fn measure_of(g: &Geometry) -> f64 {
    match g {
        Geometry::LineString(cs) => cs
            .windows(2)
            .map(|w| (w[1].x - w[0].x).hypot(w[1].y - w[0].y))
            .sum(),
        Geometry::Polygon { exterior, .. } => {
            let n = exterior.coords().len();
            let c = exterior.coords();
            let mut a = 0.0;
            for i in 0..n {
                let j = (i + 1) % n;
                a += c[i].x * c[j].y - c[j].x * c[i].y;
            }
            (a * 0.5).abs()
        }
        _ => 0.0,
    }
}

// ── GeoJSON polygon parsing (from polygonize) ────────────────────────────────

fn parse_polygons(geojson: &str) -> Result<Vec<Geometry>, ToolError> {
    let v: Value = serde_json::from_str(geojson)
        .map_err(|e| ToolError::Execution(format!("failed parsing polygonize output: {e}")))?;
    let mut out = Vec::new();
    if let Some(features) = v.get("features").and_then(Value::as_array) {
        for f in features {
            let Some(coords) = f.pointer("/geometry/coordinates").and_then(Value::as_array) else {
                continue;
            };
            let mut rings = coords.iter().filter_map(ring_from_json);
            let Some(exterior) = rings.next() else {
                continue;
            };
            out.push(Geometry::Polygon {
                exterior,
                interiors: rings.collect(),
            });
        }
    }
    Ok(out)
}

fn ring_from_json(ring: &Value) -> Option<Ring> {
    let pts = ring.as_array()?;
    let mut coords: Vec<Coord> = pts
        .iter()
        .filter_map(|p| {
            let a = p.as_array()?;
            Some(Coord::xy(a.first()?.as_f64()?, a.get(1)?.as_f64()?))
        })
        .collect();
    if coords.len() >= 2 && coords.first() == coords.last() {
        coords.pop();
    }
    if coords.len() < 3 {
        return None;
    }
    Some(Ring::new(coords))
}

const NEIGH8: [(i64, i64); 8] = [
    (-1, -1),
    (-1, 0),
    (-1, 1),
    (0, -1),
    (0, 1),
    (1, -1),
    (1, 0),
    (1, 1),
];
const NEIGH4: [(i64, i64); 4] = [(-1, 0), (1, 0), (0, -1), (0, 1)];

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum FeatureType {
    Lines,
    Polygons,
}

impl FeatureType {
    fn as_str(self) -> &'static str {
        match self {
            FeatureType::Lines => "lines",
            FeatureType::Polygons => "polygons",
        }
    }
}

struct Params {
    feature_type: FeatureType,
    foreground_value: Option<f64>,
    threshold: Option<f64>,
    noise_size: usize,
    hole_size: usize,
    gap_distance: f64,
    simplify_tolerance: Option<f64>,
    band: isize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let feature_type = match parse_optional_str(args, "feature_type")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("lines") | Some("line") => FeatureType::Lines,
        Some("polygons") | Some("polygon") => FeatureType::Polygons,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "unknown feature_type '{o}' (lines or polygons)"
            )))
        }
    };
    let foreground_value = parse_optional_f64(args, "foreground_value")?;
    let threshold = parse_optional_f64(args, "threshold")?;
    let noise_size = parse_optional_u64(args, "noise_size")?.unwrap_or(4) as usize;
    let hole_size = parse_optional_u64(args, "hole_size")?.unwrap_or(4) as usize;
    let gap_distance = parse_optional_f64(args, "gap_distance")?.unwrap_or(0.0);
    if gap_distance < 0.0 {
        return Err(ToolError::Validation(
            "'gap_distance' must be non-negative".to_string(),
        ));
    }
    let simplify_tolerance = parse_optional_f64(args, "simplify_tolerance")?;
    if let Some(t) = simplify_tolerance {
        if t < 0.0 {
            return Err(ToolError::Validation(
                "'simplify_tolerance' must be non-negative".to_string(),
            ));
        }
    }
    let band = parse_optional_u64(args, "band")?.unwrap_or(1).max(1) as isize - 1;
    Ok(Params {
        feature_type,
        foreground_value,
        threshold,
        noise_size,
        hole_size,
        gap_distance,
        simplify_tolerance,
        band,
    })
}

fn parse_optional_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
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

fn parse_optional_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64().or_else(|| n.as_f64().map(|f| f.max(0.0) as u64))),
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

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{memory_store, DataType, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn raster_from(rows: usize, cols: usize, data: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: Default::default(),
            metadata: vec![],
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, data[row * cols + col])
                    .unwrap();
            }
        }
        let id = memory_store::put_raster(r);
        memory_store::make_raster_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ExtractScannedFeaturesTool.run(&args, &ctx()).unwrap();
        let layer = crate::vector_common::load_input_layer(out.outputs["output"].as_str().unwrap())
            .unwrap();
        (out, layer)
    }

    /// A thick horizontal bar thins to a single centreline polyline.
    #[test]
    fn thick_bar_becomes_one_centreline() {
        // 7 rows x 20 cols; a 3-row-thick horizontal bar in the middle.
        let rows = 7;
        let cols = 20;
        let mut d = vec![0.0; rows * cols];
        for r in 2..5 {
            for c in 2..18 {
                d[r * cols + c] = 1.0;
            }
        }
        let input = raster_from(rows, cols, &d);
        let (out, layer) =
            run(json!({ "input": input, "feature_type": "lines", "simplify_tolerance": 0.0 }));
        assert_eq!(out.outputs["feature_type"], json!("lines"));
        assert_eq!(layer.len(), 1, "one centreline");
        // Length ~ bar length (16 cells).
        let m = out.outputs["total_measure"].as_f64().unwrap();
        assert!((10.0..=17.0).contains(&m), "centreline length {m}");
    }

    /// A tiny noise speck is removed; only the real bar survives.
    #[test]
    fn noise_speck_is_removed() {
        let rows = 9;
        let cols = 20;
        let mut d = vec![0.0; rows * cols];
        for c in 2..18 {
            d[4 * cols + c] = 1.0; // 1-px bar
        }
        d[0] = 1.0; // isolated speck
        let input = raster_from(rows, cols, &d);
        let (out, layer) = run(json!({ "input": input, "feature_type": "lines", "noise_size": 3 }));
        // Speck removed -> one line only.
        assert_eq!(layer.len(), 1);
        let _ = out;
    }

    /// A filled block extracts as one polygon in polygon mode.
    #[test]
    fn polygon_mode_extracts_regions() {
        let rows = 10;
        let cols = 10;
        let mut d = vec![0.0; rows * cols];
        for r in 2..8 {
            for c in 2..8 {
                d[r * cols + c] = 1.0;
            }
        }
        let input = raster_from(rows, cols, &d);
        let (out, layer) =
            run(json!({ "input": input, "feature_type": "polygons", "simplify_tolerance": 0.0 }));
        assert_eq!(out.outputs["feature_type"], json!("polygons"));
        assert_eq!(layer.len(), 1);
        // Area ~ 36 (6x6 block).
        let m = out.outputs["total_measure"].as_f64().unwrap();
        assert!((m - 36.0).abs() < 1e-6, "polygon area {m}");
    }

    /// A T-junction traces into three arms meeting at the node.
    #[test]
    fn junction_splits_into_arms() {
        // vertical bar + horizontal bar forming a T (1-px lines).
        let rows = 11;
        let cols = 11;
        let mut d = vec![0.0; rows * cols];
        for r in 1..10 {
            d[r * cols + 5] = 1.0; // vertical
        }
        for c in 1..10 {
            d[5 * cols + c] = 1.0; // horizontal
        }
        let input = raster_from(rows, cols, &d);
        let (_, layer) =
            run(json!({ "input": input, "feature_type": "lines", "simplify_tolerance": 0.0 }));
        // A plus-sign has 4 arms from the centre node.
        assert!(layer.len() >= 3, "expected >=3 arms, got {}", layer.len());
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = ExtractScannedFeaturesTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "x.tif", "feature_type": "bogus" })).is_err());
        assert!(bad(json!({ "input": "x.tif", "feature_type": "polygons" })).is_ok());
    }
}
