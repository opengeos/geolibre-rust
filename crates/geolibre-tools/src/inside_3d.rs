//! GeoLibre tool: true 3D containment of features within closed volumes.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Inside 3D* (3D Analyst).
//!
//! Rounds 14-15 built a real 3D suite — `buffer_3d`, `near_3d`, `idw_3d`,
//! `simplify_3d_line`, `minimum_bounding_volume`, `generate_points_along_3d_lines`,
//! `calculate_missing_z_values`, `adjust_3d_z`. What that suite still could not
//! answer is the most basic volumetric question: **is this thing inside that
//! thing**.
//!
//! Every containment predicate in either registry is 2D. `select_by_location`,
//! `spatial_join`, `clip` and the point-in-polygon helpers all project to the XY
//! plane, so a sensor 200 m above a building footprint tests as "inside the
//! building" and a utility line beneath a parcel tests as "inside the parcel".
//! For subsurface utilities, airspace volumes, building interiors and plume
//! envelopes that answer is simply wrong.
//!
//! Solids follow the convention `buffer_3d` and `minimum_bounding_volume`
//! established: a `MultiPolygon` whose parts are Z-bearing triangles.
//!
//! Containment uses ray casting with Möller-Trumbore triangle intersection. The
//! classic failure is a ray striking an edge shared by two triangles and being
//! counted twice, which flips the parity — and it is not exotic: any
//! axis-aligned ray through the middle of a quad face hits that face's diagonal.
//! Hits are therefore collected with an inclusive edge test and then
//! **deduplicated by ray parameter**, collapsing each shared-edge hit back to
//! the single surface crossing it physically is. The ray direction is derived
//! deterministically from the query point's own coordinates (**not** an RNG, so
//! WASM and native agree exactly).
//!
//! Scope note: for solid *targets*, `complex` mode reports the vertex-based
//! verdict and the fraction of vertices inside — not an exact intersection
//! volume. Exact solid-solid intersection is `union_3d`'s problem.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Tests which 3D features fall inside closed volumetric containers.
pub struct Inside3dTool;

pub(crate) type Tri = [[f64; 3]; 3];

impl Tool for Inside3dTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "inside_3d",
            display_name: "Inside 3D",
            summary: "Determines which 3D features (points, polylines or solids) fall inside closed volumetric containers, optionally reporting the contained length or vertex fraction (ArcGIS Inside 3D). Every containment predicate in either registry projects to the XY plane, so a sensor above a footprint or a pipe beneath a parcel currently tests as inside it. Containers follow the triangle-mesh convention of the shipped buffer_3d and minimum_bounding_volume.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "target",
                    description: "3D features to test: points, polylines, or triangle-mesh solids.",
                    required: true,
                },
                ToolParamSpec {
                    name: "container",
                    description: "Closed 3D features (triangle-mesh solids) acting as containers.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output table path of target/container pairs. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "mode",
                    description: "'simple' (default) reports whether the target is inside; 'complex' adds the contained length and inside fraction.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_features",
                    description: "Optional path for the clipped inside portions of line targets (complex mode). If omitted, stored in memory (still returned).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["target", "container"] {
            if args
                .get(key)
                .and_then(Value::as_str)
                .map(str::trim)
                .unwrap_or("")
                .is_empty()
            {
                return Err(ToolError::Validation(format!(
                    "missing required string parameter '{key}'"
                )));
            }
        }
        parse_mode(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let target_path = required_str(args, "target")?;
        let container_path = required_str(args, "container")?;
        let output = parse_optional_str(args, "output")?;
        let output_features = parse_optional_str(args, "output_features")?;
        let complex = parse_mode(args)?;

        let targets = load_input_layer(target_path)?;
        let containers = load_input_layer(container_path)?;

        // ── Prepare containers ───────────────────────────────────────────────
        let mut solids: Vec<Solid> = Vec::new();
        let mut open_meshes = 0_u64;
        for (cid, feature) in containers.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let tris = collect_triangles(geom);
            if tris.is_empty() {
                continue;
            }
            let solid = Solid::new(cid, tris);
            // A mesh with boundary edges bounds no volume, so ray-cast parity
            // against it is arbitrary. Skip it and report the count, matching
            // union_3d — keeping it would emit exactly the garbage rows the
            // check exists to prevent.
            if !solid.closed {
                open_meshes += 1;
                continue;
            }
            solids.push(solid);
        }
        if solids.is_empty() {
            return Err(ToolError::Execution(format!(
                "container layer holds no CLOSED triangle-mesh solids (expected MultiPolygon parts with Z); {open_meshes} open mesh(es) were skipped"
            )));
        }
        ctx.progress.info(&format!(
            "{} container solid(s), {} not closed",
            solids.len(),
            open_meshes
        ));

        // ── Test targets ─────────────────────────────────────────────────────
        let mut out = Layer::new("inside_3d");
        out.add_field(FieldDef::new("target_fid", FieldType::Integer));
        out.add_field(FieldDef::new("container_fid", FieldType::Integer));
        out.add_field(FieldDef::new("inside", FieldType::Integer));
        if complex {
            out.add_field(FieldDef::new("inside_length", FieldType::Float));
            out.add_field(FieldDef::new("inside_fraction", FieldType::Float));
        }

        let mut clipped = Layer::new("inside_3d_parts");
        clipped.add_field(FieldDef::new("target_fid", FieldType::Integer));
        clipped.add_field(FieldDef::new("container_fid", FieldType::Integer));
        clipped.crs = targets.crs.clone();

        let mut pair_count = 0_u64;
        let mut inside_count = 0_u64;

        for (tid, feature) in targets.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            for solid in &solids {
                // Broad phase: bounding boxes must overlap at all.
                if !solid.bbox_may_contain(geom) {
                    continue;
                }

                let verdict = match geom {
                    Geometry::Point(c) => {
                        let inside = solid.contains(c.x, c.y, c.z.unwrap_or(0.0));
                        Some((inside, 0.0, if inside { 1.0 } else { 0.0 }, Vec::new()))
                    }
                    Geometry::MultiPoint(cs) => {
                        let n = cs.len().max(1) as f64;
                        let hits = cs
                            .iter()
                            .filter(|c| solid.contains(c.x, c.y, c.z.unwrap_or(0.0)))
                            .count() as f64;
                        Some((hits > 0.0, 0.0, hits / n, Vec::new()))
                    }
                    Geometry::LineString(cs) => Some(line_verdict(cs, solid)),
                    Geometry::MultiLineString(parts) => {
                        let mut total_in = 0.0;
                        let mut total_len = 0.0;
                        let mut spans = Vec::new();
                        for cs in parts {
                            let (_, l_in, _, mut sp) = line_verdict(cs, solid);
                            total_in += l_in;
                            total_len += polyline_length_3d(cs);
                            spans.append(&mut sp);
                        }
                        let frac = if total_len > 0.0 {
                            total_in / total_len
                        } else {
                            0.0
                        };
                        Some((total_in > 0.0, total_in, frac, spans))
                    }
                    // Solid targets: vertex-based verdict only (see module docs).
                    Geometry::MultiPolygon(_) | Geometry::Polygon { .. } => {
                        let verts = collect_vertices(geom);
                        if verts.is_empty() {
                            None
                        } else {
                            let hits = verts
                                .iter()
                                .filter(|v| solid.contains(v[0], v[1], v[2]))
                                .count() as f64;
                            let frac = hits / verts.len() as f64;
                            Some((hits > 0.0, 0.0, frac, Vec::new()))
                        }
                    }
                    _ => None,
                };

                let Some((inside, inside_len, frac, spans)) = verdict else {
                    continue;
                };
                // Every pair that reached a verdict counts as evaluated; only
                // the containing ones count as inside.
                pair_count += 1;
                if !inside {
                    continue;
                }
                inside_count += 1;

                let mut attrs: Vec<(&str, FieldValue)> = vec![
                    ("target_fid", FieldValue::Integer(tid as i64)),
                    ("container_fid", FieldValue::Integer(solid.fid as i64)),
                    ("inside", FieldValue::Integer(1)),
                ];
                if complex {
                    attrs.push(("inside_length", FieldValue::Float(inside_len)));
                    attrs.push(("inside_fraction", FieldValue::Float(frac)));
                }
                out.add_feature(None, &attrs)
                    .map_err(|e| ToolError::Execution(format!("failed writing row: {e}")))?;

                if complex && !spans.is_empty() {
                    clipped
                        .add_feature(
                            Some(Geometry::MultiLineString(spans)),
                            &[
                                ("target_fid", FieldValue::Integer(tid as i64)),
                                ("container_fid", FieldValue::Integer(solid.fid as i64)),
                            ],
                        )
                        .map_err(|e| {
                            ToolError::Execution(format!("failed writing clipped part: {e}"))
                        })?;
                }
            }
            ctx.progress
                .progress((tid as f64 + 1.0) / targets.len().max(1) as f64);
        }

        let out_path = write_or_store_layer(out, output)?;
        let parts_path = write_or_store_layer(clipped, output_features)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("output_features".to_string(), json!(parts_path));
        outputs.insert("pair_count".to_string(), json!(pair_count));
        outputs.insert("inside_count".to_string(), json!(inside_count));
        outputs.insert("container_count".to_string(), json!(solids.len()));
        outputs.insert("open_container_count".to_string(), json!(open_meshes));
        Ok(ToolRunResult { outputs })
    }
}

/// Classifies a polyline against a solid, returning
/// `(inside, inside_length, fraction, inside_spans)`.
fn line_verdict(cs: &[Coord], solid: &Solid) -> (bool, f64, f64, Vec<Vec<Coord>>) {
    let total = polyline_length_3d(cs);
    let mut inside_len = 0.0;
    let mut spans: Vec<Vec<Coord>> = Vec::new();

    for w in cs.windows(2) {
        let a = [w[0].x, w[0].y, w[0].z.unwrap_or(0.0)];
        let b = [w[1].x, w[1].y, w[1].z.unwrap_or(0.0)];
        // Parameters where the segment crosses the mesh, plus the endpoints.
        let mut ts = solid.segment_crossings(a, b);
        ts.push(0.0);
        ts.push(1.0);
        ts.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        ts.dedup_by(|x, y| (*x - *y).abs() < 1e-12);

        for pair in ts.windows(2) {
            let (t0, t1) = (pair[0], pair[1]);
            if t1 - t0 < 1e-12 {
                continue;
            }
            // Classify the span by a single interior sample.
            let tm = (t0 + t1) / 2.0;
            let p = lerp3(a, b, tm);
            if !solid.contains(p[0], p[1], p[2]) {
                continue;
            }
            let p0 = lerp3(a, b, t0);
            let p1 = lerp3(a, b, t1);
            inside_len += dist3(p0, p1);
            spans.push(vec![
                Coord::xyz(p0[0], p0[1], p0[2]),
                Coord::xyz(p1[0], p1[1], p1[2]),
            ]);
        }
    }

    let frac = if total > 0.0 { inside_len / total } else { 0.0 };
    (inside_len > 0.0, inside_len, frac, spans)
}

// ── Solid ─────────────────────────────────────────────────────────────────────

pub(crate) struct Solid {
    pub(crate) fid: usize,
    pub(crate) tris: Vec<Tri>,
    pub(crate) min: [f64; 3],
    pub(crate) max: [f64; 3],
    pub(crate) closed: bool,
}

impl Solid {
    pub(crate) fn new(fid: usize, tris: Vec<Tri>) -> Self {
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for t in &tris {
            for v in t {
                for k in 0..3 {
                    min[k] = min[k].min(v[k]);
                    max[k] = max[k].max(v[k]);
                }
            }
        }
        let closed = is_closed(&tris);
        Self {
            fid,
            tris,
            min,
            max,
            closed,
        }
    }

    /// Cheap rejection: does the target's bounding box overlap the solid's?
    fn bbox_may_contain(&self, geom: &Geometry) -> bool {
        let verts = collect_vertices(geom);
        if verts.is_empty() {
            return false;
        }
        let mut tmin = [f64::INFINITY; 3];
        let mut tmax = [f64::NEG_INFINITY; 3];
        for v in &verts {
            for k in 0..3 {
                tmin[k] = tmin[k].min(v[k]);
                tmax[k] = tmax[k].max(v[k]);
            }
        }
        (0..3).all(|k| tmin[k] <= self.max[k] && tmax[k] >= self.min[k])
    }

    /// Ray-casting containment test. Odd crossing count means inside.
    pub(crate) fn contains(&self, x: f64, y: f64, z: f64) -> bool {
        let origin = [x, y, z];
        if (0..3).any(|k| origin[k] < self.min[k] || origin[k] > self.max[k]) {
            return false;
        }
        let dir = ray_direction(origin);
        self.crossing_params(origin, dir, f64::INFINITY).len() % 2 == 1
    }

    /// Forward ray/mesh crossing parameters, deduplicated.
    ///
    /// A ray striking an edge shared by two triangles is reported by *both* of
    /// them at the same `t`. Counting those separately breaks parity — which is
    /// exactly what happens for an axis-aligned ray through the middle of a
    /// quad face, since the face's two triangles meet on its diagonal. Merging
    /// hits at (near-)equal `t` collapses each such shared-edge crossing back to
    /// the single surface crossing it physically is.
    fn crossing_params(&self, origin: [f64; 3], dir: [f64; 3], t_max: f64) -> Vec<f64> {
        let mut ts: Vec<f64> = Vec::new();
        for tri in &self.tris {
            if let Some(t) = ray_triangle(origin, dir, tri) {
                if t > 1e-12 && t < t_max {
                    ts.push(t);
                }
            }
        }
        ts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        ts.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
        ts
    }

    /// Parameters in `(0, 1)` where the segment crosses the mesh.
    fn segment_crossings(&self, a: [f64; 3], b: [f64; 3]) -> Vec<f64> {
        let dir = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let len = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
        if len < 1e-15 {
            return Vec::new();
        }
        self.crossing_params(a, dir, 1.0 - 1e-12)
    }
}

/// A mesh bounds a volume only if every edge is shared by exactly two triangles.
fn is_closed(tris: &[Tri]) -> bool {
    let mut edges: BTreeMap<(u64, u64, u64, u64, u64, u64), usize> = BTreeMap::new();
    for t in tris {
        for k in 0..3 {
            let p = t[k];
            let q = t[(k + 1) % 3];
            // Undirected edge key, quantised so float noise does not split it.
            let a = (quant(p[0]), quant(p[1]), quant(p[2]));
            let b = (quant(q[0]), quant(q[1]), quant(q[2]));
            let key = if a <= b {
                (a.0, a.1, a.2, b.0, b.1, b.2)
            } else {
                (b.0, b.1, b.2, a.0, a.1, a.2)
            };
            *edges.entry(key).or_insert(0) += 1;
        }
    }
    !edges.is_empty() && edges.values().all(|c| *c == 2)
}

fn quant(v: f64) -> u64 {
    // 1e-9 quantisation: fine enough not to merge distinct vertices, coarse
    // enough to survive the float noise of mesh construction.
    ((v / 1e-9).round() as i64) as u64
}

/// Möller-Trumbore ray/triangle intersection, returning the ray parameter of a
/// forward hit.
///
/// The barycentric test is deliberately *inclusive* at the edges: a hit exactly
/// on a shared edge is reported by both adjacent triangles, and the caller
/// dedupes by `t`. Excluding edge hits instead would drop the crossing
/// altogether whenever a ray passes through a face diagonal.
fn ray_triangle(origin: [f64; 3], dir: [f64; 3], tri: &Tri) -> Option<f64> {
    const EPS: f64 = 1e-12;
    const EDGE: f64 = 1e-9;
    let e1 = sub3(tri[1], tri[0]);
    let e2 = sub3(tri[2], tri[0]);
    let p = cross3(dir, e2);
    let det = dot3(e1, p);
    if det.abs() < EPS {
        // Ray parallel to the triangle plane.
        return None;
    }
    let inv = 1.0 / det;
    let tvec = sub3(origin, tri[0]);
    let u = dot3(tvec, p) * inv;
    if !(-EDGE..=1.0 + EDGE).contains(&u) {
        return None;
    }
    let q = cross3(tvec, e1);
    let v = dot3(dir, q) * inv;
    if v < -EDGE || u + v > 1.0 + EDGE {
        return None;
    }
    let t = dot3(e2, q) * inv;
    if t <= EPS {
        return None;
    }
    Some(t)
}

/// Deterministic pseudo-random ray direction derived from the query point, so
/// results are reproducible on every platform (no RNG, WASM-safe).
///
/// A single cast suffices because shared-edge hits are resolved by deduplicating
/// on the ray parameter rather than by re-casting.
fn ray_direction(p: [f64; 3]) -> [f64; 3] {
    let seed = p[0] * 12.9898 + p[1] * 78.233 + p[2] * 37.719;
    let a = fract(seed.sin() * 43758.5453) * std::f64::consts::TAU;
    let b = fract((seed * 1.618).cos() * 24634.6345) * std::f64::consts::PI;
    [b.sin() * a.cos(), b.sin() * a.sin(), b.cos()]
}

fn fract(v: f64) -> f64 {
    v - v.floor()
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + t * (b[0] - a[0]),
        a[1] + t * (b[1] - a[1]),
        a[2] + t * (b[2] - a[2]),
    ]
}

fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub3(a, b);
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

fn polyline_length_3d(cs: &[Coord]) -> f64 {
    cs.windows(2)
        .map(|w| {
            dist3(
                [w[0].x, w[0].y, w[0].z.unwrap_or(0.0)],
                [w[1].x, w[1].y, w[1].z.unwrap_or(0.0)],
            )
        })
        .sum()
}

/// Extracts triangles from a triangle-mesh geometry (each MultiPolygon part is
/// one triangle, as `buffer_3d` and `minimum_bounding_volume` emit). Rings with
/// more than three vertices are fan-triangulated.
pub(crate) fn collect_triangles(geom: &Geometry) -> Vec<Tri> {
    let mut out = Vec::new();
    let mut push_ring = |cs: &[Coord]| {
        if cs.len() < 3 {
            return;
        }
        let p: Vec<[f64; 3]> = cs.iter().map(|c| [c.x, c.y, c.z.unwrap_or(0.0)]).collect();
        for k in 1..(p.len() - 1) {
            out.push([p[0], p[k], p[k + 1]]);
        }
    };
    match geom {
        Geometry::MultiPolygon(parts) => {
            for (ext, _) in parts {
                push_ring(&ext.0);
            }
        }
        Geometry::Polygon { exterior, .. } => push_ring(&exterior.0),
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                out.extend(collect_triangles(g));
            }
        }
        _ => {}
    }
    out
}

fn collect_vertices(geom: &Geometry) -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    let push = |c: &Coord, out: &mut Vec<[f64; 3]>| out.push([c.x, c.y, c.z.unwrap_or(0.0)]);
    match geom {
        Geometry::Point(c) => push(c, &mut out),
        Geometry::MultiPoint(cs) | Geometry::LineString(cs) => {
            for c in cs {
                push(c, &mut out);
            }
        }
        Geometry::MultiLineString(parts) => {
            for cs in parts {
                for c in cs {
                    push(c, &mut out);
                }
            }
        }
        Geometry::Polygon { exterior, .. } => {
            for c in &exterior.0 {
                push(c, &mut out);
            }
        }
        Geometry::MultiPolygon(parts) => {
            for (ext, _) in parts {
                for c in &ext.0 {
                    push(c, &mut out);
                }
            }
        }
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                out.extend(collect_vertices(g));
            }
        }
    }
    out
}

fn parse_mode(args: &ToolArgs) -> Result<bool, ToolError> {
    match parse_optional_str(args, "mode")? {
        None => Ok(false),
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "simple" => Ok(false),
            "complex" => Ok(true),
            other => Err(ToolError::Validation(format!(
                "unknown mode '{other}' (expected 'simple' or 'complex')"
            ))),
        },
    }
}

fn required_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

// `box_mesh` — the axis-aligned closed test box — now lives in `mesh3d`, where
// the round-17 3D tools can use it outside `cfg(test)` too (`intersect_3d`
// emits one as an intersection's bounding solid).
#[cfg(test)]
pub(crate) use crate::mesh3d::box_mesh;

#[cfg(test)]
mod tests {
    use super::*;
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

    fn layer_of(name: &str, geoms: Vec<Geometry>) -> String {
        let mut l = Layer::new(name);
        l.geom_type = Some(GeometryType::MultiPolygon);
        for g in geoms {
            l.add_feature(Some(g), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn unit_box() -> String {
        layer_of("box", vec![box_mesh([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])])
    }

    fn run(target: String, container: String, extra: Value) -> (Layer, ToolRunResult) {
        let mut obj = serde_json::Map::new();
        obj.insert("target".to_string(), json!(target));
        obj.insert("container".to_string(), json!(container));
        if let Value::Object(m) = extra {
            for (k, v) in m {
                obj.insert(k, v);
            }
        }
        let args: ToolArgs = serde_json::from_value(Value::Object(obj)).unwrap();
        let res = Inside3dTool.run(&args, &ctx()).unwrap();
        let table = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (table, res)
    }

    /// The core case 2D containment gets wrong: a point directly above the
    /// footprint is NOT inside the volume.
    #[test]
    fn point_above_the_footprint_is_not_inside() {
        let inside_pt = Geometry::Point(Coord::xyz(5.0, 5.0, 5.0));
        let above_pt = Geometry::Point(Coord::xyz(5.0, 5.0, 200.0));
        let target = layer_of("pts", vec![inside_pt, above_pt]);
        let (table, res) = run(target, unit_box(), json!({}));

        assert_eq!(res.outputs["inside_count"], json!(1));
        assert_eq!(table.len(), 1);
        let i = table.schema.field_index("target_fid").unwrap();
        assert_eq!(
            table.iter().next().unwrap().attributes[i],
            FieldValue::Integer(0),
            "only the point at z=5 is inside; the one at z=200 is above the box"
        );
    }

    /// A point below the volume is likewise outside.
    #[test]
    fn point_below_the_volume_is_not_inside() {
        let target = layer_of("pts", vec![Geometry::Point(Coord::xyz(5.0, 5.0, -50.0))]);
        let (_, res) = run(target, unit_box(), json!({}));
        assert_eq!(res.outputs["inside_count"], json!(0));
    }

    /// Points outside in XY are outside regardless of Z.
    #[test]
    fn point_outside_in_plan_is_not_inside() {
        let target = layer_of("pts", vec![Geometry::Point(Coord::xyz(50.0, 5.0, 5.0))]);
        let (_, res) = run(target, unit_box(), json!({}));
        assert_eq!(res.outputs["inside_count"], json!(0));
    }

    /// A line passing through the box reports the contained length. The box is
    /// 10 wide, so a line crossing it along X has exactly 10 units inside.
    #[test]
    fn complex_mode_measures_contained_length() {
        let line = Geometry::LineString(vec![
            Coord::xyz(-10.0, 5.0, 5.0),
            Coord::xyz(20.0, 5.0, 5.0),
        ]);
        let target = layer_of("lines", vec![line]);
        let (table, _) = run(target, unit_box(), json!({ "mode": "complex" }));

        assert_eq!(table.len(), 1);
        let li = table.schema.field_index("inside_length").unwrap();
        let fi = table.schema.field_index("inside_fraction").unwrap();
        let f = table.iter().next().unwrap();
        let FieldValue::Float(len) = f.attributes[li] else {
            panic!("expected float length");
        };
        let FieldValue::Float(frac) = f.attributes[fi] else {
            panic!("expected float fraction");
        };
        assert!(
            (len - 10.0).abs() < 1e-6,
            "expected 10 units inside the box, got {len}"
        );
        // The line is 30 long overall, so a third of it is inside.
        assert!((frac - 1.0 / 3.0).abs() < 1e-6, "got fraction {frac}");
    }

    /// A line entirely below the box has nothing inside, even though it passes
    /// directly under it in plan.
    #[test]
    fn line_under_the_box_has_no_contained_length() {
        let line = Geometry::LineString(vec![
            Coord::xyz(-10.0, 5.0, -5.0),
            Coord::xyz(20.0, 5.0, -5.0),
        ]);
        let target = layer_of("lines", vec![line]);
        let (_, res) = run(target, unit_box(), json!({ "mode": "complex" }));
        assert_eq!(res.outputs["inside_count"], json!(0));
    }

    /// Multiple containers are each reported separately.
    #[test]
    fn reports_one_row_per_container_hit() {
        let containers = layer_of(
            "boxes",
            vec![
                box_mesh([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]),
                // A second box overlapping the first around the query point.
                box_mesh([4.0, 4.0, 4.0], [20.0, 20.0, 20.0]),
            ],
        );
        let target = layer_of("pts", vec![Geometry::Point(Coord::xyz(5.0, 5.0, 5.0))]);
        let (table, res) = run(target, containers, json!({}));
        assert_eq!(table.len(), 2, "the point is inside both boxes");
        assert_eq!(res.outputs["container_count"], json!(2));
    }

    /// An open (non-closed) mesh bounds no volume, so it is skipped and counted
    /// rather than ray-cast — parity against an open mesh is arbitrary, which is
    /// exactly the garbage the check exists to prevent.
    #[test]
    fn open_meshes_are_skipped_not_used() {
        // A single triangle bounds no volume.
        let open = || {
            Geometry::MultiPolygon(vec![(
                wbvector::Ring::new(vec![
                    Coord::xyz(0.0, 0.0, 0.0),
                    Coord::xyz(1.0, 0.0, 0.0),
                    Coord::xyz(0.0, 1.0, 0.0),
                ]),
                Vec::new(),
            )])
        };
        let target = layer_of("pts", vec![Geometry::Point(Coord::xyz(5.0, 5.0, 5.0))]);

        // Alongside a real solid: the open mesh is counted but contributes no row.
        let mixed = layer_of(
            "mixed",
            vec![box_mesh([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]), open()],
        );
        let (table, res) = run(target.clone(), mixed, json!({}));
        assert_eq!(res.outputs["open_container_count"], json!(1));
        assert_eq!(
            res.outputs["container_count"],
            json!(1),
            "only the closed box is a container"
        );
        assert_eq!(table.len(), 1, "the open mesh must not produce a row");

        // On its own it leaves nothing to test against, which is an error.
        let only_open = layer_of("open", vec![open()]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "target": target, "container": only_open })).unwrap();
        assert!(Inside3dTool.run(&args, &ctx()).is_err());
    }

    /// The containment test is deterministic: repeated runs agree exactly,
    /// which is what the seeded (non-RNG) ray direction guarantees.
    #[test]
    fn results_are_deterministic() {
        let target = layer_of(
            "pts",
            (0..25)
                .map(|i| {
                    Geometry::Point(Coord::xyz((i % 5) as f64 * 3.0, (i / 5) as f64 * 3.0, 5.0))
                })
                .collect(),
        );
        let (_, a) = run(target.clone(), unit_box(), json!({}));
        let (_, b) = run(target, unit_box(), json!({}));
        assert_eq!(a.outputs["inside_count"], b.outputs["inside_count"]);
    }

    #[test]
    fn rejects_bad_parameters() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(Inside3dTool.validate(&args).is_err());

        let t = layer_of("pts", vec![Geometry::Point(Coord::xyz(0.0, 0.0, 0.0))]);
        let args: ToolArgs = serde_json::from_value(
            json!({ "target": t.clone(), "container": unit_box(), "mode": "fuzzy" }),
        )
        .unwrap();
        assert!(Inside3dTool.validate(&args).is_err());

        // A container layer with no meshes is a run-time error.
        let empty = layer_of("empty", vec![Geometry::Point(Coord::xyz(0.0, 0.0, 0.0))]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "target": t, "container": empty })).unwrap();
        assert!(Inside3dTool.run(&args, &ctx()).is_err());
    }
}
