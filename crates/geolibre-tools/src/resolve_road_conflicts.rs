//! GeoLibre tool: displace roads whose symbolized widths graphically conflict.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Resolve Road Conflicts* (Cartography),
//! which also subsumes *Propagate Displacement*. It is the road-vs-road sibling of
//! the shipped `resolve_building_conflicts` (buildings vs. roads) and reuses the
//! localized vertex-displacement idea behind `rubbersheet_features`. No bundled
//! whitebox tool performs displacement cartography.
//!
//! At a target map scale, each road is drawn as a strip of some *symbolized width*.
//! Two roads conflict when their strips overlap or crowd (their centreline
//! separation is smaller than half the sum of their symbol widths plus an optional
//! `gap`). The tool iteratively pushes the *lower-hierarchy* road's vertices away
//! from the higher-hierarchy road, along the local conflict normal and by the
//! overlap amount, until conflicts clear or a max-iteration cap is hit. Endpoints
//! are pinned by default and the per-vertex displacement field is lightly smoothed,
//! so junctions/connectivity are preserved and shifts taper along each line.
//!
//! Symbol widths are read from `symbol_width_field` (per feature) or a fixed
//! `symbol_width`. With `scale` supplied the widths are treated as page
//! millimetres and converted to ground units (`mm * scale / 1000`); without it
//! they are taken directly as map units. Geographic layers (CRS None or 4326) are
//! projected to a local equirectangular metre space for the computation and
//! projected back on output, so distances stay metric.
//!
//! Hierarchy is read from `hierarchy_field` (numeric; **lower value = more
//! important**). In a conflict the more-important road stays fixed and the
//! less-important road moves; equal ranks share the displacement. An optional
//! `links` output writes displacement links (from → to lines) that
//! `rubbersheet_features` can consume to propagate the shift onto companion layers.
//!
//! v1 scope: pairwise, Jacobi-style iterative displacement with endpoint pinning
//! and displacement-field smoothing. Conflict counting uses exact segment-to-
//! segment distance (truthful), but displacement is applied only at existing
//! vertices — very sparsely digitized lines may need densifying for full
//! resolution. It is not a full least-squares network adjustment, roads are not
//! split, and a road boxed between two more-important roads (or two roads meeting
//! at a pinned junction) may retain a residual conflict — reported honestly in
//! `conflicts_after`, never silently hidden.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Metres per degree of latitude (local equirectangular scaling for geographic CRS).
const DEG_M: f64 = 111_320.0;
/// Fraction of the computed overlap applied each iteration (Jacobi damping).
const DAMPING: f64 = 0.5;
/// Overlap below which a pair is considered clear (working units).
const TOL: f64 = 1e-4;

pub struct ResolveRoadConflictsTool;

impl Tool for ResolveRoadConflictsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "resolve_road_conflicts",
            display_name: "Resolve Road Conflicts",
            summary: "Displace roads whose symbolized (drawn-width) representations overlap or crowd at a target scale: lower-hierarchy roads are pushed clear of higher-hierarchy ones with endpoint-pinned, smoothed vertex shifts, optionally emitting displacement links, like ArcGIS Resolve Road Conflicts / Propagate Displacement.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input road line layer (LineString / MultiLineString).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output displaced road layer with 'status' and 'shift' fields. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "symbol_width",
                    description: "Fixed symbolized width for all roads, in map units (or page mm when 'scale' is set). Used when 'symbol_width_field' is absent.",
                    required: false,
                },
                ToolParamSpec {
                    name: "symbol_width_field",
                    description: "Numeric field giving each road's symbolized width (map units, or page mm when 'scale' is set). Falls back to 'symbol_width' when a value is missing.",
                    required: false,
                },
                ToolParamSpec {
                    name: "hierarchy_field",
                    description: "Numeric field ranking road importance (lower value = more important, held fixed). Absent: all roads share displacement equally.",
                    required: false,
                },
                ToolParamSpec {
                    name: "scale",
                    description: "Target map scale denominator. When set, symbol widths are page millimetres converted to ground units (mm * scale / 1000). Default: widths taken directly as map units.",
                    required: false,
                },
                ToolParamSpec {
                    name: "gap",
                    description: "Extra clear spacing to keep between symbol strips, in map units. Default 0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_iter",
                    description: "Maximum displacement iterations. Default 50.",
                    required: false,
                },
                ToolParamSpec {
                    name: "pin_endpoints",
                    description: "Keep the first/last vertex of every road fixed to preserve junctions/connectivity (default true).",
                    required: false,
                },
                ToolParamSpec {
                    name: "links",
                    description: "Optional output path for displacement links (from→to lines) consumable by rubbersheet_features.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;

        // Geographic layers are projected to local metres so symbol widths (which
        // are metric once scaled) and centreline distances share the same units.
        let geographic = matches!(layer.crs_epsg(), None | Some(4326));
        let proj = Projector::build(&layer, geographic);

        // ── Load roads into working coordinates ───────────────────────────────
        let width_idx = prm
            .symbol_width_field
            .as_deref()
            .and_then(|f| layer.schema.field_index(f));
        let rank_idx = prm
            .hierarchy_field
            .as_deref()
            .and_then(|f| layer.schema.field_index(f));

        let mut roads: Vec<Road> = Vec::new();
        for (fidx, f) in layer.features.iter().enumerate() {
            let Some(geom) = f.geometry.as_ref() else {
                continue;
            };
            let (parts, is_multi) = match geom {
                Geometry::LineString(cs) if cs.len() >= 2 => {
                    (vec![cs.iter().map(|c| proj.fwd(c.x, c.y)).collect()], false)
                }
                Geometry::MultiLineString(ls) => {
                    let parts: Vec<Vec<(f64, f64)>> = ls
                        .iter()
                        .filter(|cs| cs.len() >= 2)
                        .map(|cs| cs.iter().map(|c| proj.fwd(c.x, c.y)).collect())
                        .collect();
                    if parts.is_empty() {
                        continue;
                    }
                    (parts, true)
                }
                _ => continue,
            };

            let width = width_idx
                .and_then(|i| f.attributes.get(i))
                .and_then(FieldValue::as_f64)
                .filter(|w| w.is_finite() && *w > 0.0)
                .or(prm.symbol_width)
                .unwrap_or(0.0);
            let half = 0.5 * width * prm.width_scale;
            let rank = rank_idx
                .and_then(|i| f.attributes.get(i))
                .and_then(FieldValue::as_f64)
                .filter(|r| r.is_finite())
                .unwrap_or(0.0);

            let orig = parts.clone();
            let pin_first = parts.iter().map(|_| false).collect();
            let pin_last = parts.iter().map(|_| false).collect();
            roads.push(Road {
                src: fidx,
                parts,
                orig,
                is_multi,
                half,
                rank,
                pin_first,
                pin_last,
                max_shift: 0.0,
            });
        }
        if roads.is_empty() {
            return Err(ToolError::Execution(
                "no road line features in input".to_string(),
            ));
        }

        // Junctions: endpoint coordinates shared by two or more road-part ends.
        // Only these are pinned (when pin_endpoints), so connectivity survives
        // while free dangling ends stay free to move.
        if prm.pin_endpoints {
            let mut counts: BTreeMap<(i64, i64), u32> = BTreeMap::new();
            let key = |x: f64, y: f64| ((x * 1000.0).round() as i64, (y * 1000.0).round() as i64);
            for r in &roads {
                for p in &r.parts {
                    *counts.entry(key(p[0].0, p[0].1)).or_default() += 1;
                    let last = *p.last().unwrap();
                    *counts.entry(key(last.0, last.1)).or_default() += 1;
                }
            }
            for r in roads.iter_mut() {
                for (pi, p) in r.parts.iter().enumerate() {
                    let last = *p.last().unwrap();
                    r.pin_first[pi] = counts[&key(p[0].0, p[0].1)] >= 2;
                    r.pin_last[pi] = counts[&key(last.0, last.1)] >= 2;
                }
            }
        }

        let conflicts_before = conflict_count(&roads, prm.gap);
        ctx.progress.info(&format!(
            "{} road(s); {} symbol-overlap conflict(s) before",
            roads.len(),
            conflicts_before
        ));

        // ── Iterative displacement ────────────────────────────────────────────
        let mut iters = 0usize;
        for _ in 0..prm.max_iter {
            iters += 1;
            let mut disp: Vec<Vec<Vec<(f64, f64)>>> = roads
                .iter()
                .map(|r| r.parts.iter().map(|p| vec![(0.0, 0.0); p.len()]).collect())
                .collect();
            let mut max_overlap = 0.0f64;

            for a in 0..roads.len() {
                for b in 0..roads.len() {
                    if a == b || roads[a].half <= 0.0 || roads[b].half <= 0.0 {
                        continue;
                    }
                    let w = move_weight(roads[a].rank, roads[b].rank);
                    if w <= 0.0 {
                        continue;
                    }
                    let req = roads[a].half + roads[b].half + prm.gap;
                    for (pi, part) in roads[a].parts.iter().enumerate() {
                        for (ki, &(x, y)) in part.iter().enumerate() {
                            let (d, bp, bdir) = nearest_on_road(x, y, &roads[b].parts);
                            if d < req {
                                let overlap = req - d;
                                max_overlap = max_overlap.max(overlap);
                                let (nx, ny) = away_normal(x, y, bp, bdir);
                                disp[a][pi][ki].0 += nx * overlap * w;
                                disp[a][pi][ki].1 += ny * overlap * w;
                            }
                        }
                    }
                }
            }

            if max_overlap <= TOL {
                break;
            }

            // Smooth the displacement field along each line, then apply (damped),
            // keeping pinned junction endpoints fixed.
            for (road, disp_r) in roads.iter_mut().zip(disp.iter()) {
                for (pi, part) in road.parts.iter_mut().enumerate() {
                    let n = part.len();
                    let raw = &disp_r[pi];
                    let mut sm = raw.clone();
                    for k in 1..n.saturating_sub(1) {
                        sm[k].0 = 0.5 * raw[k].0 + 0.25 * raw[k - 1].0 + 0.25 * raw[k + 1].0;
                        sm[k].1 = 0.5 * raw[k].1 + 0.25 * raw[k - 1].1 + 0.25 * raw[k + 1].1;
                    }
                    let pin_first = road.pin_first[pi];
                    let pin_last = road.pin_last[pi];
                    for (k, (v, s)) in part.iter_mut().zip(sm.iter()).enumerate() {
                        if (k == 0 && pin_first) || (k == n - 1 && pin_last) {
                            continue;
                        }
                        v.0 += s.0 * DAMPING;
                        v.1 += s.1 * DAMPING;
                    }
                }
            }
        }

        // Per-road maximum vertex shift (working units == metres).
        for r in roads.iter_mut() {
            let mut m = 0.0f64;
            for (part, orig) in r.parts.iter().zip(r.orig.iter()) {
                for (&(x, y), &(ox, oy)) in part.iter().zip(orig.iter()) {
                    m = m.max((x - ox).hypot(y - oy));
                }
            }
            r.max_shift = m;
        }

        let conflicts_after = conflict_count(&roads, prm.gap);
        let displaced = roads.iter().filter(|r| r.max_shift > TOL).count();
        let total_shift: f64 = roads.iter().map(|r| r.max_shift).sum();
        let mean_shift = total_shift / roads.len() as f64;
        let max_shift = roads.iter().map(|r| r.max_shift).fold(0.0, f64::max);

        ctx.progress.info(&format!(
            "{iters} iteration(s); {conflicts_after} conflict(s) after; {displaced} road(s) displaced"
        ));

        // ── Build displaced road output ───────────────────────────────────────
        let mut out = Layer::new("resolved_roads");
        if let Some(gt) = layer.geom_type {
            out = out.with_geom_type(gt);
        }
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for field in layer.schema.fields() {
            out.add_field(field.clone());
        }
        out.add_field(FieldDef::new("status", FieldType::Text));
        out.add_field(FieldDef::new("shift", FieldType::Float));

        // Map source feature index -> road index for in-order emission.
        let mut road_of: BTreeMap<usize, usize> = BTreeMap::new();
        for (ri, r) in roads.iter().enumerate() {
            road_of.insert(r.src, ri);
        }

        for (fidx, f) in layer.features.iter().enumerate() {
            let mut attrs = f.attributes.clone();
            match road_of.get(&fidx) {
                Some(&ri) => {
                    let r = &roads[ri];
                    let status = if r.max_shift > TOL {
                        "displaced"
                    } else {
                        "unchanged"
                    };
                    attrs.push(FieldValue::Text(status.to_string()));
                    attrs.push(FieldValue::Float(r.max_shift));
                    out.push(Feature {
                        fid: 0,
                        geometry: Some(road_to_geometry(r, &proj)),
                        attributes: attrs,
                    });
                }
                None => {
                    // Non-line (or degenerate) feature: pass through unchanged.
                    attrs.push(FieldValue::Text("unchanged".to_string()));
                    attrs.push(FieldValue::Float(0.0));
                    out.push(Feature {
                        fid: 0,
                        geometry: f.geometry.clone(),
                        attributes: attrs,
                    });
                }
            }
        }

        let out_path = write_or_store_layer(out, output)?;

        // ── Optional displacement links ───────────────────────────────────────
        let mut link_count = 0usize;
        let mut link_path: Option<String> = None;
        if let Some(links_out) = prm.links.as_deref() {
            let mut ll = Layer::new("displacement_links").with_geom_type(GeometryType::LineString);
            if let Some(epsg) = layer.crs_epsg() {
                ll = ll.with_crs_epsg(epsg);
            }
            ll.add_field(FieldDef::new("road_fid", FieldType::Integer));
            ll.add_field(FieldDef::new("shift", FieldType::Float));
            for r in &roads {
                for (part, orig) in r.parts.iter().zip(r.orig.iter()) {
                    for (&(x, y), &(ox, oy)) in part.iter().zip(orig.iter()) {
                        let s = (x - ox).hypot(y - oy);
                        if s > TOL {
                            let (fx, fy) = proj.inv(ox, oy);
                            let (tx, ty) = proj.inv(x, y);
                            ll.push(Feature {
                                fid: 0,
                                geometry: Some(Geometry::LineString(vec![
                                    Coord::xy(fx, fy),
                                    Coord::xy(tx, ty),
                                ])),
                                attributes: vec![
                                    FieldValue::Integer(r.src as i64),
                                    FieldValue::Float(s),
                                ],
                            });
                            link_count += 1;
                        }
                    }
                }
            }
            link_path = Some(write_or_store_layer(ll, Some(links_out))?);
        }

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("road_count".to_string(), json!(roads.len()));
        outputs.insert("conflicts_before".to_string(), json!(conflicts_before));
        outputs.insert("conflicts_after".to_string(), json!(conflicts_after));
        outputs.insert("displaced".to_string(), json!(displaced));
        outputs.insert("iterations".to_string(), json!(iters));
        outputs.insert("mean_shift".to_string(), json!(mean_shift));
        outputs.insert("max_shift".to_string(), json!(max_shift));
        if let Some(lp) = link_path {
            outputs.insert("links".to_string(), json!(lp));
            outputs.insert("link_count".to_string(), json!(link_count));
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── Road model ───────────────────────────────────────────────────────────────

struct Road {
    src: usize,
    parts: Vec<Vec<(f64, f64)>>,
    orig: Vec<Vec<(f64, f64)>>,
    is_multi: bool,
    half: f64,
    rank: f64,
    /// Whether the first/last vertex of each part is a pinned junction.
    pin_first: Vec<bool>,
    pin_last: Vec<bool>,
    max_shift: f64,
}

/// Fraction of a conflict's overlap that road `a` (rank `ra`) absorbs relative to
/// road `b` (rank `rb`). Lower rank == more important == stays fixed.
fn move_weight(ra: f64, rb: f64) -> f64 {
    if ra > rb {
        1.0 // a is less important -> a moves fully
    } else if ra < rb {
        0.0 // a is more important -> a stays
    } else {
        0.5 // equal -> share
    }
}

/// Unit vector pushing point `(x, y)` away from its nearest point `bp` on the other
/// road; if it lies on the road, push perpendicular to that road's local segment.
fn away_normal(x: f64, y: f64, bp: (f64, f64), bdir: (f64, f64)) -> (f64, f64) {
    let (mut nx, mut ny) = (x - bp.0, y - bp.1);
    let len = nx.hypot(ny);
    if len < 1e-9 {
        let (dx, dy) = bdir;
        let l = dx.hypot(dy).max(1e-9);
        return (-dy / l, dx / l);
    }
    nx /= len;
    ny /= len;
    (nx, ny)
}

/// Minimum distance from `(x, y)` to any segment of `parts`, plus the nearest point
/// and the direction of the nearest segment.
fn nearest_on_road(x: f64, y: f64, parts: &[Vec<(f64, f64)>]) -> (f64, (f64, f64), (f64, f64)) {
    let mut best = f64::INFINITY;
    let mut bp = (x, y);
    let mut bdir = (1.0, 0.0);
    for part in parts {
        for w in part.windows(2) {
            let (ax, ay) = w[0];
            let (bx, by) = w[1];
            let (px, py) = nearest_on_seg(x, y, ax, ay, bx, by);
            let d = (x - px).hypot(y - py);
            if d < best {
                best = d;
                bp = (px, py);
                bdir = (bx - ax, by - ay);
            }
        }
    }
    (best, bp, bdir)
}

fn nearest_on_seg(px: f64, py: f64, ax: f64, ay: f64, bx: f64, by: f64) -> (f64, f64) {
    let dx = bx - ax;
    let dy = by - ay;
    let len2 = dx * dx + dy * dy;
    let t = if len2 <= 0.0 {
        0.0
    } else {
        (((px - ax) * dx + (py - ay) * dy) / len2).clamp(0.0, 1.0)
    };
    (ax + t * dx, ay + t * dy)
}

/// Exact minimum distance between two segments (0 if they cross), so conflict
/// counting is truthful even where the closest approach is interior-to-interior
/// on sparse near-parallel lines (a vertex-only test would under-report).
fn seg_seg_dist(a0: (f64, f64), a1: (f64, f64), b0: (f64, f64), b1: (f64, f64)) -> f64 {
    if segments_cross(a0, a1, b0, b1) {
        return 0.0;
    }
    let pd = |p: (f64, f64), s0: (f64, f64), s1: (f64, f64)| {
        let (px, py) = nearest_on_seg(p.0, p.1, s0.0, s0.1, s1.0, s1.1);
        (p.0 - px).hypot(p.1 - py)
    };
    pd(a0, b0, b1)
        .min(pd(a1, b0, b1))
        .min(pd(b0, a0, a1))
        .min(pd(b1, a0, a1))
}

fn orient(p: (f64, f64), q: (f64, f64), r: (f64, f64)) -> f64 {
    (q.0 - p.0) * (r.1 - p.1) - (q.1 - p.1) * (r.0 - p.0)
}

fn segments_cross(a0: (f64, f64), a1: (f64, f64), b0: (f64, f64), b1: (f64, f64)) -> bool {
    let d1 = orient(a0, a1, b0);
    let d2 = orient(a0, a1, b1);
    let d3 = orient(b0, b1, a0);
    let d4 = orient(b0, b1, a1);
    ((d1 > 0.0) != (d2 > 0.0)) && ((d3 > 0.0) != (d4 > 0.0))
}

/// Exact minimum centreline distance between two roads.
fn min_road_dist(a: &Road, b: &Road) -> f64 {
    let mut m = f64::INFINITY;
    for pa in &a.parts {
        for wa in pa.windows(2) {
            for pb in &b.parts {
                for wb in pb.windows(2) {
                    m = m.min(seg_seg_dist(wa[0], wa[1], wb[0], wb[1]));
                    if m == 0.0 {
                        return 0.0;
                    }
                }
            }
        }
    }
    m
}

/// Number of unordered road pairs whose symbol strips overlap or crowd.
fn conflict_count(roads: &[Road], gap: f64) -> usize {
    let mut c = 0;
    for i in 0..roads.len() {
        for j in (i + 1)..roads.len() {
            if roads[i].half <= 0.0 || roads[j].half <= 0.0 {
                continue;
            }
            let req = roads[i].half + roads[j].half + gap;
            if min_road_dist(&roads[i], &roads[j]) < req - TOL {
                c += 1;
            }
        }
    }
    c
}

fn road_to_geometry(r: &Road, proj: &Projector) -> Geometry {
    if r.is_multi {
        Geometry::MultiLineString(
            r.parts
                .iter()
                .map(|p| {
                    p.iter()
                        .map(|&(x, y)| {
                            let (lx, ly) = proj.inv(x, y);
                            Coord::xy(lx, ly)
                        })
                        .collect()
                })
                .collect(),
        )
    } else {
        Geometry::LineString(
            r.parts[0]
                .iter()
                .map(|&(x, y)| {
                    let (lx, ly) = proj.inv(x, y);
                    Coord::xy(lx, ly)
                })
                .collect(),
        )
    }
}

// ── Projection (local equirectangular for geographic layers) ──────────────────

struct Projector {
    geographic: bool,
    lon0: f64,
    lat0: f64,
    cos0: f64,
}

impl Projector {
    fn build(layer: &Layer, geographic: bool) -> Self {
        if !geographic {
            return Projector {
                geographic: false,
                lon0: 0.0,
                lat0: 0.0,
                cos0: 1.0,
            };
        }
        let (mut sx, mut sy, mut n) = (0.0, 0.0, 0u64);
        for f in layer.features.iter() {
            if let Some(g) = f.geometry.as_ref() {
                for c in g.all_coords() {
                    sx += c.x;
                    sy += c.y;
                    n += 1;
                }
            }
        }
        let lon0 = if n > 0 { sx / n as f64 } else { 0.0 };
        let lat0 = if n > 0 { sy / n as f64 } else { 0.0 };
        Projector {
            geographic: true,
            lon0,
            lat0,
            cos0: (lat0.to_radians()).cos().max(1e-6),
        }
    }

    fn fwd(&self, x: f64, y: f64) -> (f64, f64) {
        if self.geographic {
            ((x - self.lon0) * self.cos0 * DEG_M, (y - self.lat0) * DEG_M)
        } else {
            (x, y)
        }
    }

    fn inv(&self, x: f64, y: f64) -> (f64, f64) {
        if self.geographic {
            (self.lon0 + x / (self.cos0 * DEG_M), self.lat0 + y / DEG_M)
        } else {
            (x, y)
        }
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    symbol_width: Option<f64>,
    symbol_width_field: Option<String>,
    hierarchy_field: Option<String>,
    width_scale: f64,
    gap: f64,
    max_iter: usize,
    pin_endpoints: bool,
    links: Option<String>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let symbol_width = opt_f64(args, "symbol_width")?;
    if let Some(w) = symbol_width {
        if !(w.is_finite() && w > 0.0) {
            return Err(ToolError::Validation(
                "'symbol_width' must be a positive number".to_string(),
            ));
        }
    }
    let symbol_width_field = parse_optional_str(args, "symbol_width_field")?.map(str::to_string);
    let hierarchy_field = parse_optional_str(args, "hierarchy_field")?.map(str::to_string);
    if symbol_width.is_none() && symbol_width_field.is_none() {
        return Err(ToolError::Validation(
            "provide 'symbol_width' or 'symbol_width_field' to size the road symbols".to_string(),
        ));
    }
    let scale = opt_f64(args, "scale")?;
    let width_scale = match scale {
        None => 1.0,
        Some(s) if s.is_finite() && s > 0.0 => s / 1000.0,
        Some(_) => {
            return Err(ToolError::Validation(
                "'scale' must be a positive number".to_string(),
            ))
        }
    };
    let gap = opt_f64(args, "gap")?.unwrap_or(0.0);
    if !(gap.is_finite() && gap >= 0.0) {
        return Err(ToolError::Validation(
            "'gap' must be a non-negative number".to_string(),
        ));
    }
    let max_iter = match opt_f64(args, "max_iter")? {
        None => 50,
        Some(v) if v >= 1.0 && v.is_finite() => v as usize,
        Some(_) => {
            return Err(ToolError::Validation(
                "'max_iter' must be an integer >= 1".to_string(),
            ))
        }
    };
    let pin_endpoints = opt_bool(args, "pin_endpoints")?.unwrap_or(true);
    let links = parse_optional_str(args, "links")?.map(str::to_string);

    Ok(Params {
        symbol_width,
        symbol_width_field,
        hierarchy_field,
        width_scale,
        gap,
        max_iter,
        pin_endpoints,
        links,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn opt_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
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

fn opt_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
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

    /// (polyline points, symbol width, hierarchy rank)
    type RoadSpec = (Vec<(f64, f64)>, f64, f64);

    /// Build a road layer from a list of road specs.
    fn road_layer(roads: &[RoadSpec], epsg: Option<u32>) -> String {
        let mut l = Layer::new("roads").with_geom_type(GeometryType::LineString);
        if let Some(e) = epsg {
            l = l.with_crs_epsg(e);
        }
        l.add_field(FieldDef::new("w", FieldType::Float));
        l.add_field(FieldDef::new("rank", FieldType::Integer));
        for (pts, w, rank) in roads {
            l.add_feature(
                Some(Geometry::LineString(
                    pts.iter().map(|(x, y)| Coord::xy(*x, *y)).collect(),
                )),
                &[("w", (*w).into()), ("rank", (*rank as i64).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ResolveRoadConflictsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn geom_pts(g: &Geometry) -> Vec<(f64, f64)> {
        match g {
            Geometry::LineString(cs) => cs.iter().map(|c| (c.x, c.y)).collect(),
            _ => panic!("expected linestring"),
        }
    }

    /// Two near-parallel equal-rank roads closer than their combined symbol width
    /// are pushed apart until the conflict clears; each moves (share).
    #[test]
    fn separates_parallel_equal_rank_roads() {
        // Road A at y=0, road B at y=6, both width 10 (half 5 each -> need >=10 apart).
        let a: Vec<(f64, f64)> = (0..=10).map(|i| (i as f64 * 5.0, 0.0)).collect();
        let b: Vec<(f64, f64)> = (0..=10).map(|i| (i as f64 * 5.0, 6.0)).collect();
        let path = road_layer(&[(a, 10.0, 0.0), (b, 10.0, 0.0)], Some(3857));
        let (out, layer) = run(json!({
            "input": path, "symbol_width_field": "w", "hierarchy_field": "rank", "gap": 0.0
        }));
        assert_eq!(out.outputs["conflicts_before"], json!(1));
        assert_eq!(out.outputs["conflicts_after"], json!(0));
        // Both roads moved (shared displacement); midpoints separated by >= 10.
        let ay = geom_pts(layer.features[0].geometry.as_ref().unwrap())[5].1;
        let by = geom_pts(layer.features[1].geometry.as_ref().unwrap())[5].1;
        assert!((by - ay) >= 10.0 - 1e-3, "gap {} too small", by - ay);
        assert!(ay < 0.0 && by > 6.0, "roads should move apart");
    }

    /// The higher-hierarchy road (lower rank) stays fixed; the lower one moves.
    #[test]
    fn fixes_higher_hierarchy_road() {
        let a: Vec<(f64, f64)> = (0..=10).map(|i| (i as f64 * 5.0, 0.0)).collect();
        let b: Vec<(f64, f64)> = (0..=10).map(|i| (i as f64 * 5.0, 6.0)).collect();
        // A rank 0 (important, fixed), B rank 5 (moves).
        let path = road_layer(&[(a, 10.0, 0.0), (b, 10.0, 5.0)], Some(3857));
        let (out, layer) = run(json!({
            "input": path, "symbol_width_field": "w", "hierarchy_field": "rank"
        }));
        assert_eq!(out.outputs["conflicts_after"], json!(0));
        let si = layer.schema.field_index("status").unwrap();
        assert_eq!(layer.features[0].attributes[si].as_str(), Some("unchanged"));
        assert_eq!(layer.features[1].attributes[si].as_str(), Some("displaced"));
        // Road A interior vertex unmoved.
        let ay = geom_pts(layer.features[0].geometry.as_ref().unwrap())[5].1;
        assert!(ay.abs() < 1e-9, "important road moved: {ay}");
        // Road B pushed up clear.
        let by = geom_pts(layer.features[1].geometry.as_ref().unwrap())[5].1;
        assert!(by >= 10.0 - 1e-3, "B not clear: {by}");
    }

    /// Roads that already clear their symbol widths are left unchanged.
    #[test]
    fn leaves_clear_roads_unchanged() {
        let a: Vec<(f64, f64)> = vec![(0.0, 0.0), (50.0, 0.0)];
        let b: Vec<(f64, f64)> = vec![(0.0, 100.0), (50.0, 100.0)];
        let path = road_layer(&[(a, 4.0, 0.0), (b, 4.0, 0.0)], Some(3857));
        let (out, _l) = run(json!({ "input": path, "symbol_width_field": "w" }));
        assert_eq!(out.outputs["conflicts_before"], json!(0));
        assert_eq!(out.outputs["displaced"], json!(0));
        assert_eq!(out.outputs["max_shift"], json!(0.0));
    }

    /// Displacement links are emitted for moved vertices.
    #[test]
    fn emits_displacement_links() {
        let a: Vec<(f64, f64)> = (0..=10).map(|i| (i as f64 * 5.0, 0.0)).collect();
        let b: Vec<(f64, f64)> = (0..=10).map(|i| (i as f64 * 5.0, 6.0)).collect();
        let path = road_layer(&[(a, 10.0, 0.0), (b, 10.0, 5.0)], Some(3857));
        let links_path =
            std::env::temp_dir().join(format!("rrc_links_{}.geojson", std::process::id()));
        let links_str = links_path.to_string_lossy().to_string();
        let (out, _l) = run(json!({
            "input": path, "symbol_width_field": "w", "hierarchy_field": "rank",
            "links": links_str
        }));
        // link_count present and positive.
        assert!(
            out.outputs
                .get("link_count")
                .and_then(|v| v.as_u64())
                .unwrap_or(0)
                > 0
        );
        let ll = load_input_layer(out.outputs["links"].as_str().unwrap()).unwrap();
        assert!(!ll.features.is_empty());
        let _ = std::fs::remove_file(&links_path);
    }

    /// Geographic (EPSG:4326) input is handled metrically and conflicts clear.
    #[test]
    fn handles_geographic_crs() {
        // Two E-W roads ~6 m apart near the equator; width 10 m needs >=10 m sep.
        let dy = 6.0 / DEG_M; // ~6 metres in degrees latitude
        let a: Vec<(f64, f64)> = (0..=10).map(|i| (i as f64 * 0.0005, 0.0)).collect();
        let b: Vec<(f64, f64)> = (0..=10).map(|i| (i as f64 * 0.0005, dy)).collect();
        let path = road_layer(&[(a, 10.0, 0.0), (b, 10.0, 0.0)], Some(4326));
        let (out, _l) = run(json!({ "input": path, "symbol_width_field": "w" }));
        assert_eq!(out.outputs["conflicts_before"], json!(1));
        assert_eq!(out.outputs["conflicts_after"], json!(0));
    }

    #[test]
    fn rejects_missing_input() {
        let args: ToolArgs = serde_json::from_value(json!({ "symbol_width": 5.0 })).unwrap();
        assert!(ResolveRoadConflictsTool.validate(&args).is_err());
    }

    #[test]
    fn rejects_missing_symbol_width() {
        let path = road_layer(&[(vec![(0.0, 0.0), (10.0, 0.0)], 5.0, 0.0)], Some(3857));
        let args: ToolArgs = serde_json::from_value(json!({ "input": path })).unwrap();
        assert!(ResolveRoadConflictsTool.validate(&args).is_err());
    }
}
