//! GeoLibre tool: seamline cut polygons for overlapping raster tiles.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Build Seamlines* (Data Management).
//! The bundled `mosaic` and `mosaic_with_feathering` decide overlap by simple
//! ordering and blending — they never compute **where** the cut ought to go.
//! That works for radiometrically consistent tiles and fails visibly on
//! anything else: date-mismatched orthophotos, differing sun angles, seasonal
//! change and cloud edges all produce an obvious straight join running through
//! buildings and fields. Seamline generation is the standard fix and had no
//! equivalent anywhere in either registry.
//!
//! The repo has no mosaic-dataset concept, so the interface is a raster list
//! plus optional footprints rather than a managed dataset. Every cell of the
//! mosaic extent is assigned to exactly one input image, and each image's
//! assigned cells are polygonised into its seamline polygon — so the outputs
//! tile the extent once, with no overlaps and no gaps, ready to hand to
//! `mosaic_with_feathering` as per-image clip masks.
//!
//! Assignment methods:
//!   * `voronoi` (default) — nearest footprint centroid; the classic
//!     equal-distance partition.
//!   * `order` — priority order wins, reproducing "copy footprint" behaviour.
//!   * `radiometry` — the image whose value at that cell is closest to the
//!     consensus (mean) of all covering images, which puts the cut where the
//!     images already agree.
//!   * `edge_detection` — radiometry weighted by local gradient, so the seam
//!     prefers to run *along* real edges (field boundaries, roads) rather than
//!     across them.
//!
//! **Scope for v1:** footprints, when not supplied, are derived as each
//! raster's valid-data bounding rectangle rather than a traced NoData outline,
//! and `min_thinness_ratio` / `max_sliver_size` are not implemented.
//! `min_region_size` is.

use std::collections::BTreeMap;

use geo::{Area, BooleanOps, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::common::load_input_raster;
use crate::vector_common::{
    geometry_contains_point, load_input_layer, parse_optional_str, write_or_store_layer,
};

pub struct BuildSeamlinesTool;

impl Tool for BuildSeamlinesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "build_seamlines",
            display_name: "Build Seamlines",
            summary: "Given a set of overlapping rasters (and optional footprints), compute where each image should be cut so the mosaic joins along a minimally visible seam, emitting one non-overlapping seamline polygon per image via a Voronoi, priority-order, radiometric or edge-aware assignment. Like ArcGIS Build Seamlines.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "inputs",
                    description: "Comma/semicolon-separated raster paths to mosaic.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output seamline polygon path (one feature per input raster). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "footprints",
                    description: "Optional polygon layer giving each raster's valid extent, in input order. Derived from valid data when omitted.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'voronoi' (default), 'order', 'radiometry' or 'edge_detection'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "sort_ascending",
                    description: "For method 'order': true (default) gives earlier rasters priority, false gives later ones priority.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Analysis cell size (default: the finest input cell size).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band used for the radiometric methods (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_region_size",
                    description: "Reassign assigned regions smaller than this many cells to a neighbouring image (default 0, i.e. keep all).",
                    required: false,
                },
                ToolParamSpec {
                    name: "blend_width",
                    description: "Feather width in map units, written as an attribute for downstream mosaicking (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "blend_type",
                    description: "'both' (default), 'inside' or 'outside'; written as an attribute.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        let inputs = split_list(require_str(args, "inputs")?);
        if inputs.len() < 2 {
            return Err(ToolError::Validation(
                "'inputs' must list at least 2 rasters".to_string(),
            ));
        }
        parse_method(args)?;
        parse_blend_type(args)?;
        for key in ["cell_size", "blend_width", "min_region_size"] {
            if let Some(v) = parse_optional_f64(args, key)? {
                if v < 0.0 {
                    return Err(ToolError::Validation(format!(
                        "'{key}' must be non-negative"
                    )));
                }
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let paths = split_list(require_str(args, "inputs")?);
        let output = parse_optional_str(args, "output")?;
        let method = parse_method(args)?;
        let blend_type = parse_blend_type(args)?;
        let blend_width = parse_optional_f64(args, "blend_width")?.unwrap_or(0.0);
        let ascending = parse_optional_bool(args, "sort_ascending")?.unwrap_or(true);
        let band = (parse_optional_f64(args, "band")?.unwrap_or(1.0).max(1.0) as isize) - 1;
        let min_region = parse_optional_f64(args, "min_region_size")?.unwrap_or(0.0) as usize;

        ctx.progress
            .info(&format!("loading {} raster(s)", paths.len()));
        let mut rasters: Vec<Raster> = Vec::with_capacity(paths.len());
        for p in &paths {
            rasters.push(load_input_raster(p)?);
        }
        let m = rasters.len();

        // Footprints: supplied, or each raster's valid-data rectangle.
        let supplied: Option<Vec<Geometry>> = match parse_optional_str(args, "footprints")? {
            Some(p) => {
                let l = load_input_layer(p)?;
                let g: Vec<Geometry> = l.iter().filter_map(|f| f.geometry.clone()).collect();
                if g.len() < m {
                    return Err(ToolError::Validation(format!(
                        "'footprints' has {} feature(s) but {m} raster(s) were supplied",
                        g.len()
                    )));
                }
                Some(g)
            }
            None => None,
        };
        let footprints: Vec<Geometry> = (0..m)
            .map(|i| match &supplied {
                Some(g) => g[i].clone(),
                None => valid_data_rect(&rasters[i], band),
            })
            .collect();

        // Analysis grid over the union extent.
        let cell = parse_optional_f64(args, "cell_size")?.unwrap_or_else(|| {
            rasters
                .iter()
                .map(|r| r.cell_size_x.min(r.cell_size_y))
                .fold(f64::INFINITY, f64::min)
        });
        if !(cell > 0.0 && cell.is_finite()) {
            return Err(ToolError::Execution(
                "could not determine a positive analysis cell size".to_string(),
            ));
        }
        let (mut min_x, mut min_y) = (f64::INFINITY, f64::INFINITY);
        let (mut max_x, mut max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
        for g in &footprints {
            if let Some(b) = g.bbox() {
                min_x = min_x.min(b.min_x);
                min_y = min_y.min(b.min_y);
                max_x = max_x.max(b.max_x);
                max_y = max_y.max(b.max_y);
            }
        }
        if !min_x.is_finite() {
            return Err(ToolError::Execution(
                "input rasters have no valid extent".to_string(),
            ));
        }
        let cols = (((max_x - min_x) / cell).ceil() as usize).max(1);
        let rows = (((max_y - min_y) / cell).ceil() as usize).max(1);
        if rows.saturating_mul(cols) > 40_000_000 {
            return Err(ToolError::Execution(format!(
                "analysis grid would be {rows}x{cols}; supply a coarser 'cell_size'"
            )));
        }

        let centroids: Vec<(f64, f64)> = footprints
            .iter()
            .map(|g| centroid(g).unwrap_or((0.0, 0.0)))
            .collect();

        // Assign every cell to exactly one image.
        ctx.progress.info("assigning cells to images");
        let mut owner = vec![-1_i32; rows * cols];
        // Cells covered by exactly two images, and which pair — the radiometric
        // seam is solved per pair.
        let mut pair_of: Vec<Option<(usize, usize)>> = vec![None; rows * cols];
        let mut covered_cells = 0usize;
        let mut overlap_cells = 0usize;
        for row in 0..rows {
            let y = max_y - (row as f64 + 0.5) * cell;
            for col in 0..cols {
                let x = min_x + (col as f64 + 0.5) * cell;
                let mut covering: Vec<usize> = Vec::new();
                for i in 0..m {
                    if geometry_contains_point(&footprints[i], x, y)
                        && sample(&rasters[i], band, x, y).is_some()
                    {
                        covering.push(i);
                    }
                }
                if covering.is_empty() {
                    continue;
                }
                covered_cells += 1;
                if covering.len() > 1 {
                    overlap_cells += 1;
                }
                if covering.len() == 2 {
                    pair_of[row * cols + col] = Some((covering[0], covering[1]));
                }
                let pick = match method {
                    Method::Order => {
                        if ascending {
                            covering[0]
                        } else {
                            covering[covering.len() - 1]
                        }
                    }
                    Method::Voronoi => *covering
                        .iter()
                        .min_by(|&&a, &&b| {
                            let da = (centroids[a].0 - x).hypot(centroids[a].1 - y);
                            let db = (centroids[b].0 - x).hypot(centroids[b].1 - y);
                            da.total_cmp(&db)
                        })
                        .expect("covering is non-empty"),
                    // Radiometric methods need the whole overlap at once, not a
                    // per-cell decision, so they are resolved in a second pass
                    // below. Seed with the Voronoi answer, which is also the
                    // final answer wherever three or more images overlap.
                    Method::Radiometry | Method::EdgeDetection => *covering
                        .iter()
                        .min_by(|&&a, &&b| {
                            let da = (centroids[a].0 - x).hypot(centroids[a].1 - y);
                            let db = (centroids[b].0 - x).hypot(centroids[b].1 - y);
                            da.total_cmp(&db)
                        })
                        .expect("covering is non-empty"),
                };
                owner[row * cols + col] = pick as i32;
            }
            ctx.progress.progress((row as f64 + 1.0) / rows as f64);
        }
        if covered_cells == 0 {
            return Err(ToolError::Execution(
                "no analysis cell was covered by any input raster; check CRS and extents"
                    .to_string(),
            ));
        }

        if matches!(method, Method::Radiometry | Method::EdgeDetection) {
            ctx.progress.info("solving least-cost radiometric seams");
            radiometric_seams(
                &mut owner, &pair_of, &rasters, &centroids, band, rows, cols, min_x, max_y, cell,
                method == Method::EdgeDetection,
            );
        }

        let reassigned = if min_region > 0 {
            drop_small_regions(&mut owner, rows, cols, min_region)
        } else {
            0
        };

        // Polygonise each image's assigned cells.
        ctx.progress.info("polygonising seamline regions");
        let mut out = Layer::new("seamlines").with_geom_type(GeometryType::Polygon);
        if let Some(e) = rasters[0].crs.epsg {
            out = out.with_crs_epsg(e);
        }
        out.add_field(FieldDef::new("IMAGE_ID", FieldType::Integer));
        out.add_field(FieldDef::new("SOURCE", FieldType::Text));
        out.add_field(FieldDef::new("CELL_COUNT", FieldType::Integer));
        out.add_field(FieldDef::new("AREA", FieldType::Float));
        out.add_field(FieldDef::new("BLEND_WIDTH", FieldType::Float));
        out.add_field(FieldDef::new("BLEND_TYPE", FieldType::Text));

        let mut total_area = 0.0;
        let mut emitted = 0usize;
        for (i, source) in paths.iter().enumerate() {
            let count = owner.iter().filter(|&&o| o == i as i32).count();
            if count == 0 {
                continue;
            }
            let mp = polygonise(&owner, rows, cols, i as i32, min_x, max_y, cell);
            if mp.0.is_empty() {
                continue;
            }
            let area = mp.unsigned_area();
            total_area += area;
            emitted += 1;
            out.add_feature(
                Some(multipolygon_to_geometry(&mp)),
                &[
                    ("IMAGE_ID", FieldValue::Integer(i as i64)),
                    ("SOURCE", FieldValue::Text(source.clone())),
                    ("CELL_COUNT", FieldValue::Integer(count as i64)),
                    ("AREA", FieldValue::Float(area)),
                    ("BLEND_WIDTH", FieldValue::Float(blend_width)),
                    ("BLEND_TYPE", FieldValue::Text(blend_type.to_string())),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed adding seamline: {e}")))?;
        }
        if emitted == 0 {
            return Err(ToolError::Execution(
                "no seamline region could be built".to_string(),
            ));
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("raster_count".to_string(), json!(m));
        outputs.insert("seamline_count".to_string(), json!(emitted));
        outputs.insert("covered_cells".to_string(), json!(covered_cells));
        outputs.insert("overlap_cells".to_string(), json!(overlap_cells));
        outputs.insert("reassigned_cells".to_string(), json!(reassigned));
        outputs.insert("total_area".to_string(), json!(total_area));
        outputs.insert("method".to_string(), json!(method.name()));
        outputs.insert("cell_size".to_string(), json!(cell));
        Ok(ToolRunResult { outputs })
    }
}

// ── Assignment helpers ──────────────────────────────────────────────────────

/// Each raster's valid-data bounding rectangle, as a polygon.
fn valid_data_rect(r: &Raster, band: isize) -> Geometry {
    let (mut lo_r, mut lo_c) = (usize::MAX, usize::MAX);
    let (mut hi_r, mut hi_c) = (0usize, 0usize);
    let mut any = false;
    for row in 0..r.rows {
        for col in 0..r.cols {
            let v = r.get(band.max(0), row as isize, col as isize);
            if v == r.nodata || !v.is_finite() {
                continue;
            }
            any = true;
            lo_r = lo_r.min(row);
            lo_c = lo_c.min(col);
            hi_r = hi_r.max(row);
            hi_c = hi_c.max(col);
        }
    }
    if !any {
        lo_r = 0;
        lo_c = 0;
        hi_r = r.rows.saturating_sub(1);
        hi_c = r.cols.saturating_sub(1);
    }
    let x0 = r.x_min + lo_c as f64 * r.cell_size_x;
    let x1 = r.x_min + (hi_c + 1) as f64 * r.cell_size_x;
    let y1 = r.y_max() - lo_r as f64 * r.cell_size_y;
    let y0 = r.y_max() - (hi_r + 1) as f64 * r.cell_size_y;
    Geometry::polygon(
        vec![
            Coord::xy(x0, y0),
            Coord::xy(x1, y0),
            Coord::xy(x1, y1),
            Coord::xy(x0, y1),
            Coord::xy(x0, y0),
        ],
        vec![],
    )
}

fn sample(r: &Raster, band: isize, x: f64, y: f64) -> Option<f64> {
    let (col, row) = r.world_to_pixel(x, y)?;
    if row < 0 || col < 0 || row >= r.rows as isize || col >= r.cols as isize {
        return None;
    }
    let v = r.get(band.max(0), row, col);
    (v != r.nodata && v.is_finite()).then_some(v)
}

/// Central-difference gradient magnitude at a map location.
fn gradient(r: &Raster, band: isize, x: f64, y: f64, step: f64) -> f64 {
    let e = sample(r, band, x + step, y);
    let w = sample(r, band, x - step, y);
    let n = sample(r, band, x, y + step);
    let s = sample(r, band, x, y - step);
    match (e, w, n, s) {
        (Some(e), Some(w), Some(n), Some(s)) => ((e - w) / 2.0).hypot((n - s) / 2.0),
        _ => 0.0,
    }
}

/// Reassigns connected regions smaller than `min_cells` to the most common
/// neighbouring owner, so tiny specks do not become their own seamline parts.
fn drop_small_regions(owner: &mut [i32], rows: usize, cols: usize, min_cells: usize) -> usize {
    let mut visited = vec![false; rows * cols];
    let mut reassigned = 0usize;
    for start in 0..(rows * cols) {
        if visited[start] || owner[start] < 0 {
            continue;
        }
        let id = owner[start];
        let mut comp = Vec::new();
        let mut stack = vec![start];
        visited[start] = true;
        while let Some(u) = stack.pop() {
            comp.push(u);
            let (r, c) = (u / cols, u % cols);
            for (dr, dc) in [(-1_isize, 0_isize), (1, 0), (0, -1), (0, 1)] {
                let (nr, nc) = (r as isize + dr, c as isize + dc);
                if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                    continue;
                }
                let v = nr as usize * cols + nc as usize;
                if !visited[v] && owner[v] == id {
                    visited[v] = true;
                    stack.push(v);
                }
            }
        }
        if comp.len() >= min_cells {
            continue;
        }
        // Majority vote among the region's differing neighbours.
        let mut votes: BTreeMap<i32, usize> = BTreeMap::new();
        for &u in &comp {
            let (r, c) = (u / cols, u % cols);
            for (dr, dc) in [(-1_isize, 0_isize), (1, 0), (0, -1), (0, 1)] {
                let (nr, nc) = (r as isize + dr, c as isize + dc);
                if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                    continue;
                }
                let o = owner[nr as usize * cols + nc as usize];
                if o >= 0 && o != id {
                    *votes.entry(o).or_default() += 1;
                }
            }
        }
        if let Some((&winner, _)) = votes.iter().max_by_key(|&(_, v)| *v) {
            for &u in &comp {
                owner[u] = winner;
            }
            reassigned += comp.len();
        }
    }
    reassigned
}

/// Polygonises the cells owned by `id`.
///
/// Cells are grouped into horizontal runs before unioning — a per-cell union
/// would be correct but pathologically slow, and the run rectangles produce the
/// identical outline.
fn polygonise(
    owner: &[i32],
    rows: usize,
    cols: usize,
    id: i32,
    x_min: f64,
    y_max: f64,
    cell: f64,
) -> MultiPolygon {
    let mut parts: Vec<Polygon> = Vec::new();
    for row in 0..rows {
        let mut col = 0usize;
        while col < cols {
            if owner[row * cols + col] != id {
                col += 1;
                continue;
            }
            let start = col;
            while col < cols && owner[row * cols + col] == id {
                col += 1;
            }
            let x0 = x_min + start as f64 * cell;
            let x1 = x_min + col as f64 * cell;
            let y1 = y_max - row as f64 * cell;
            let y0 = y_max - (row + 1) as f64 * cell;
            parts.push(Polygon::new(
                LineString::new(vec![
                    GeoCoord { x: x0, y: y0 },
                    GeoCoord { x: x1, y: y0 },
                    GeoCoord { x: x1, y: y1 },
                    GeoCoord { x: x0, y: y1 },
                    GeoCoord { x: x0, y: y0 },
                ]),
                vec![],
            ));
        }
    }
    if parts.is_empty() {
        return MultiPolygon(Vec::new());
    }
    // Balanced pairwise union: far cheaper than folding one run at a time.
    let mut level: Vec<MultiPolygon> = parts.into_iter().map(|p| MultiPolygon(vec![p])).collect();
    while level.len() > 1 {
        let mut next: Vec<MultiPolygon> = Vec::with_capacity(level.len().div_ceil(2));
        let mut it = level.into_iter();
        while let Some(a) = it.next() {
            match it.next() {
                Some(b) => next.push(a.union(&b)),
                None => next.push(a),
            }
        }
        level = next;
    }
    level.into_iter().next().unwrap_or(MultiPolygon(Vec::new()))
}

/// Resolves every two-image overlap by cutting along the least-cost seam.
///
/// A per-cell "which image is closer to the consensus" rule cannot work for a
/// pair: both values are equidistant from their own mean, so the test is always
/// a tie and the seam degenerates. The meaningful question is *where to put the
/// cut*, which is a path problem: find the route across the overlap that
/// minimises the radiometric mismatch |v_a - v_b| it has to cross, then give
/// each side of that route to the nearer image.
#[allow(clippy::too_many_arguments)]
fn radiometric_seams(
    owner: &mut [i32],
    pair_of: &[Option<(usize, usize)>],
    rasters: &[Raster],
    centroids: &[(f64, f64)],
    band: isize,
    rows: usize,
    cols: usize,
    x_min: f64,
    y_max: f64,
    cell: f64,
    edge_aware: bool,
) {
    // Group the overlap cells by which image pair covers them.
    let mut groups: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    for (idx, p) in pair_of.iter().enumerate() {
        if let Some(pair) = p {
            groups.entry(*pair).or_default().push(idx);
        }
    }

    for ((a, b), cells) in groups {
        if cells.len() < 2 {
            continue;
        }
        let in_group: std::collections::HashSet<usize> = cells.iter().copied().collect();
        let at = |idx: usize| -> (f64, f64) {
            let (r, c) = (idx / cols, idx % cols);
            (
                x_min + (c as f64 + 0.5) * cell,
                y_max - (r as f64 + 0.5) * cell,
            )
        };
        // Per-cell mismatch cost; edge-aware mode discounts it where the image
        // already has a strong gradient, so the seam prefers to follow real
        // features (field edges, roads) instead of cutting across them.
        let cost_at = |idx: usize| -> f64 {
            let (x, y) = at(idx);
            let va = sample(&rasters[a], band, x, y).unwrap_or(0.0);
            let vb = sample(&rasters[b], band, x, y).unwrap_or(0.0);
            let mut c = (va - vb).abs();
            if edge_aware {
                let g = gradient(&rasters[a], band, x, y, cell)
                    .max(gradient(&rasters[b], band, x, y, cell));
                c /= 1.0 + g;
            }
            c + 1e-9
        };

        // The seam runs across the overlap, perpendicular to the line joining
        // the two image centroids.
        let (ax, ay) = centroids[a];
        let (bx, by) = centroids[b];
        let (dx, dy) = (bx - ax, by - ay);
        let len = dx.hypot(dy);
        if len <= 0.0 {
            continue;
        }
        let perp = (-dy / len, dx / len);
        let proj = |idx: usize| -> f64 {
            let (x, y) = at(idx);
            x * perp.0 + y * perp.1
        };
        let lo = cells.iter().map(|&i| proj(i)).fold(f64::INFINITY, f64::min);
        let hi = cells
            .iter()
            .map(|&i| proj(i))
            .fold(f64::NEG_INFINITY, f64::max);
        let sources: Vec<usize> = cells
            .iter()
            .copied()
            .filter(|&i| proj(i) <= lo + cell)
            .collect();
        let targets: std::collections::HashSet<usize> = cells
            .iter()
            .copied()
            .filter(|&i| proj(i) >= hi - cell)
            .collect();
        if sources.is_empty() || targets.is_empty() {
            continue;
        }

        // Dijkstra confined to the overlap.
        let mut dist: BTreeMap<usize, f64> = BTreeMap::new();
        let mut back: BTreeMap<usize, usize> = BTreeMap::new();
        let mut heap: std::collections::BinaryHeap<(std::cmp::Reverse<OrderedF64>, usize)> =
            std::collections::BinaryHeap::new();
        for &s in &sources {
            dist.insert(s, cost_at(s));
            heap.push((std::cmp::Reverse(OrderedF64(cost_at(s))), s));
        }
        let mut reached = None;
        while let Some((std::cmp::Reverse(OrderedF64(d)), u)) = heap.pop() {
            if d > *dist.get(&u).unwrap_or(&f64::INFINITY) {
                continue;
            }
            if targets.contains(&u) {
                reached = Some(u);
                break;
            }
            let (r, c) = (u / cols, u % cols);
            for (dr, dc) in [
                (-1_isize, -1_isize), (-1, 0), (-1, 1),
                (0, -1),                       (0, 1),
                (1, -1),             (1, 0),   (1, 1),
            ] {
                let (nr, nc) = (r as isize + dr, c as isize + dc);
                if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                    continue;
                }
                let v = nr as usize * cols + nc as usize;
                if !in_group.contains(&v) {
                    continue;
                }
                let nd = d + cost_at(v);
                if nd < *dist.get(&v).unwrap_or(&f64::INFINITY) {
                    dist.insert(v, nd);
                    back.insert(v, u);
                    heap.push((std::cmp::Reverse(OrderedF64(nd)), v));
                }
            }
        }
        let Some(end) = reached else { continue };

        // Walk the seam back to its source.
        let mut seam: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut cur = end;
        seam.insert(cur);
        while let Some(&prev) = back.get(&cur) {
            cur = prev;
            seam.insert(cur);
        }

        // Flood the overlap from the side nearest image `a`, stopping at the
        // seam; everything the flood reaches belongs to `a`, the rest to `b`.
        let start = cells
            .iter()
            .copied()
            .filter(|i| !seam.contains(i))
            .min_by(|&i, &j| {
                let (xi, yi) = at(i);
                let (xj, yj) = at(j);
                (xi - ax)
                    .hypot(yi - ay)
                    .total_cmp(&(xj - ax).hypot(yj - ay))
            });
        let Some(start) = start else { continue };
        let mut side_a: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut stack = vec![start];
        side_a.insert(start);
        while let Some(u) = stack.pop() {
            let (r, c) = (u / cols, u % cols);
            for (dr, dc) in [(-1_isize, 0_isize), (1, 0), (0, -1), (0, 1)] {
                let (nr, nc) = (r as isize + dr, c as isize + dc);
                if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                    continue;
                }
                let v = nr as usize * cols + nc as usize;
                if in_group.contains(&v) && !seam.contains(&v) && side_a.insert(v) {
                    stack.push(v);
                }
            }
        }
        for &idx in &cells {
            owner[idx] = if side_a.contains(&idx) || seam.contains(&idx) {
                a as i32
            } else {
                b as i32
            };
        }
    }
}

/// Total-order wrapper so `f64` costs can live in a `BinaryHeap`.
#[derive(PartialEq)]
struct OrderedF64(f64);
impl Eq for OrderedF64 {}
impl PartialOrd for OrderedF64 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for OrderedF64 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0)
    }
}

fn centroid(g: &Geometry) -> Option<(f64, f64)> {
    let mut coords = g.all_coords();
    if coords.is_empty() {
        return None;
    }
    // A closed ring repeats its first vertex; averaging it twice drags the
    // centroid toward that corner (on a rectangle, by a full 10% of the span).
    if coords.len() > 2 {
        let (first, last) = (coords[0], coords[coords.len() - 1]);
        if (first.x - last.x).abs() < 1e-12 && (first.y - last.y).abs() < 1e-12 {
            coords.pop();
        }
    }
    let n = coords.len() as f64;
    Some((
        coords.iter().map(|c| c.x).sum::<f64>() / n,
        coords.iter().map(|c| c.y).sum::<f64>() / n,
    ))
}

fn multipolygon_to_geometry(mp: &MultiPolygon) -> Geometry {
    if mp.0.len() == 1 {
        let (exterior, interiors) = polygon_to_rings(&mp.0[0]);
        Geometry::Polygon {
            exterior,
            interiors,
        }
    } else {
        Geometry::MultiPolygon(mp.0.iter().map(polygon_to_rings).collect())
    }
}

fn polygon_to_rings(poly: &Polygon) -> (Ring, Vec<Ring>) {
    (
        linestring_to_ring(poly.exterior()),
        poly.interiors().iter().map(linestring_to_ring).collect(),
    )
}

fn linestring_to_ring(ls: &LineString) -> Ring {
    let mut coords: Vec<Coord> = ls.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
    if coords.len() >= 2 && coords.first().map(|c| (c.x, c.y)) == coords.last().map(|c| (c.x, c.y))
    {
        coords.pop();
    }
    Ring::new(coords)
}

// ── Params ──────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Voronoi,
    Order,
    Radiometry,
    EdgeDetection,
}

impl Method {
    fn name(self) -> &'static str {
        match self {
            Method::Voronoi => "voronoi",
            Method::Order => "order",
            Method::Radiometry => "radiometry",
            Method::EdgeDetection => "edge_detection",
        }
    }
}

fn parse_method(args: &ToolArgs) -> Result<Method, ToolError> {
    match args
        .get("method")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        None | Some("") | Some("voronoi") => Ok(Method::Voronoi),
        Some("order") | Some("copy_footprint") => Ok(Method::Order),
        Some("radiometry") => Ok(Method::Radiometry),
        Some("edge_detection") => Ok(Method::EdgeDetection),
        Some(o) => Err(ToolError::Validation(format!(
            "'method' must be 'voronoi', 'order', 'radiometry' or 'edge_detection', got '{o}'"
        ))),
    }
}

fn parse_blend_type(args: &ToolArgs) -> Result<&'static str, ToolError> {
    match args
        .get("blend_type")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        None | Some("") | Some("both") => Ok("both"),
        Some("inside") => Ok("inside"),
        Some("outside") => Ok("outside"),
        Some(o) => Err(ToolError::Validation(format!(
            "'blend_type' must be 'both', 'inside' or 'outside', got '{o}'"
        ))),
    }
}

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
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

fn split_list(s: &str) -> Vec<String> {
    s.split([',', ';'])
        .map(str::trim)
        .filter(|p| !p.is_empty())
        .map(str::to_string)
        .collect()
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// A raster whose lower-left corner is (x0, y0), `n` cells square, cell 1.
    fn tile(x0: f64, y0: f64, n: usize, f: impl Fn(f64, f64) -> f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols: n,
            rows: n,
            bands: 1,
            x_min: x0,
            y_min: y0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo::from_epsg(3857),
            metadata: Default::default(),
        });
        for row in 0..n {
            for col in 0..n {
                let x = x0 + 0.5 + col as f64;
                let y = y0 + n as f64 - 0.5 - row as f64;
                r.set(0, row as isize, col as isize, f(x, y)).unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    /// Two 20x20 tiles overlapping in x 10..20.
    fn pair() -> (String, String) {
        (
            tile(0.0, 0.0, 20, |_x, _y| 10.0),
            tile(10.0, 0.0, 20, |_x, _y| 20.0),
        )
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = BuildSeamlinesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    #[test]
    fn seamlines_tile_the_extent_without_overlapping() {
        // The whole point: outputs must partition the covered extent exactly
        // once — no double-counted pixels, no holes.
        let (a, b) = pair();
        let (out, layer) = run(json!({ "inputs": format!("{a},{b}"), "cell_size": 1.0 }));
        assert_eq!(out.outputs["seamline_count"], json!(2));
        assert!(out.outputs["overlap_cells"].as_f64().unwrap() > 0.0);

        let mps: Vec<MultiPolygon> = layer
            .iter()
            .map(|f| to_mp(f.geometry.as_ref().unwrap()))
            .collect();
        // Pairwise intersections must be empty.
        let inter = mps[0].intersection(&mps[1]).unsigned_area();
        assert!(inter < 1e-6, "seamlines overlap by {inter}");
        // Union area equals the covered extent (30 x 20 = 600).
        let union = mps[0].union(&mps[1]).unsigned_area();
        assert!((union - 600.0).abs() < 1e-6, "union area {union}");
    }

    /// AREA attribute of the seamline belonging to `img`, or 0 when absent.
    fn area_for(layer: &Layer, img: i64) -> f64 {
        let id = layer.schema.field_index("IMAGE_ID").unwrap();
        let ar = layer.schema.field_index("AREA").unwrap();
        layer
            .iter()
            .find(|f| f.attributes[id].as_i64() == Some(img))
            .map(|f| f.attributes[ar].as_f64().unwrap())
            .unwrap_or(0.0)
    }

    fn to_mp(g: &Geometry) -> MultiPolygon {
        match g {
            Geometry::Polygon {
                exterior,
                interiors,
            } => MultiPolygon(vec![Polygon::new(
                ring_ls(exterior),
                interiors.iter().map(ring_ls).collect(),
            )]),
            Geometry::MultiPolygon(ps) => MultiPolygon(
                ps.iter()
                    .map(|(e, hs)| Polygon::new(ring_ls(e), hs.iter().map(ring_ls).collect()))
                    .collect(),
            ),
            other => panic!("unexpected geometry {other:?}"),
        }
    }
    fn ring_ls(r: &Ring) -> LineString {
        LineString::new(r.coords().iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect())
    }

    #[test]
    fn voronoi_splits_the_overlap_near_the_middle() {
        let (a, b) = pair();
        let (_o, layer) = run(json!({ "inputs": format!("{a},{b}"), "cell_size": 1.0 }));
        // Tile A spans x 0..20, tile B x 10..30, so the bisector sits at x=15.
        let id = layer.schema.field_index("IMAGE_ID").unwrap();
        let first = layer
            .iter()
            .find(|f| f.attributes[id].as_i64() == Some(0))
            .unwrap();
        let max_x = first
            .geometry
            .as_ref()
            .unwrap()
            .all_coords()
            .iter()
            .map(|c| c.x)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!((max_x - 15.0).abs() < 1.5, "A extends to x = {max_x}");
    }

    #[test]
    fn order_method_gives_the_whole_overlap_to_one_image() {
        let (a, b) = pair();
        let area_of = |img: i64, asc: bool| -> f64 {
            let (_o, layer) = run(json!({
                "inputs": format!("{a},{b}"), "cell_size": 1.0,
                "method": "order", "sort_ascending": asc
            }));
            area_for(&layer, img)
        };
        // Ascending: A (first) keeps its full 20x20; B gets only 10x20.
        assert!((area_of(0, true) - 400.0).abs() < 1e-6);
        assert!((area_of(1, true) - 200.0).abs() < 1e-6);
        // Descending flips it.
        assert!((area_of(0, false) - 200.0).abs() < 1e-6);
        assert!((area_of(1, false) - 400.0).abs() < 1e-6);
    }

    #[test]
    fn radiometry_differs_from_voronoi_when_values_disagree() {
        // A ramps, B is flat: the radiometric cut should not land on the
        // geometric bisector.
        let a = tile(0.0, 0.0, 20, |x, _y| x);
        let b = tile(10.0, 0.0, 20, |_x, _y| 25.0);
        let bounds = |method: &str| -> f64 {
            let (_o, layer) = run(json!({
                "inputs": format!("{a},{b}"), "cell_size": 1.0, "method": method
            }));
            area_for(&layer, 0)
        };
        assert!(
            (bounds("voronoi") - bounds("radiometry")).abs() > 1.0,
            "radiometry produced the same cut as voronoi"
        );
    }

    #[test]
    fn blend_settings_are_carried_as_attributes() {
        let (a, b) = pair();
        let (_o, layer) = run(json!({
            "inputs": format!("{a},{b}"), "cell_size": 1.0,
            "blend_width": 7.5, "blend_type": "inside"
        }));
        let (w, t) = (
            layer.schema.field_index("BLEND_WIDTH").unwrap(),
            layer.schema.field_index("BLEND_TYPE").unwrap(),
        );
        assert_eq!(layer.features[0].attributes[w].as_f64(), Some(7.5));
        assert_eq!(layer.features[0].attributes[t].as_str(), Some("inside"));
    }

    #[test]
    fn min_region_size_absorbs_isolated_specks() {
        // 5x5 owned by image 0, with a single stray cell of image 1 in the
        // middle. A one-cell region must be absorbed by its neighbours.
        let (rows, cols) = (5usize, 5usize);
        let mut owner = vec![0_i32; rows * cols];
        owner[2 * cols + 2] = 1;
        let moved = drop_small_regions(&mut owner, rows, cols, 4);
        assert_eq!(moved, 1);
        assert!(owner.iter().all(|&o| o == 0), "speck was not absorbed");

        // A region at or above the threshold is left alone.
        let mut owner2 = vec![0_i32; rows * cols];
        for r in 0..2 {
            for c in 0..2 {
                owner2[r * cols + c] = 1;
            }
        }
        let moved2 = drop_small_regions(&mut owner2, rows, cols, 4);
        assert_eq!(moved2, 0);
        assert_eq!(owner2.iter().filter(|&&o| o == 1).count(), 4);
    }

    #[test]
    fn min_region_size_defaults_to_no_reassignment() {
        let (a, b) = pair();
        let (out, _l) = run(json!({
            "inputs": format!("{a},{b}"), "cell_size": 1.0, "method": "radiometry"
        }));
        assert_eq!(out.outputs["reassigned_cells"], json!(0));
    }

    #[test]
    fn non_overlapping_tiles_each_keep_their_own_extent() {
        let a = tile(0.0, 0.0, 10, |_x, _y| 1.0);
        let b = tile(50.0, 0.0, 10, |_x, _y| 2.0);
        let (out, _l) = run(json!({ "inputs": format!("{a},{b}"), "cell_size": 1.0 }));
        assert_eq!(out.outputs["overlap_cells"], json!(0));
        assert!((out.outputs["total_area"].as_f64().unwrap() - 200.0).abs() < 1e-6);
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            BuildSeamlinesTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "inputs": "a.tif" })).is_err());
        assert!(bad(json!({ "inputs": "a.tif,b.tif", "method": "disparity" })).is_err());
        assert!(bad(json!({ "inputs": "a.tif,b.tif", "blend_type": "sideways" })).is_err());
        assert!(bad(json!({ "inputs": "a.tif,b.tif" })).is_ok());
    }
}
