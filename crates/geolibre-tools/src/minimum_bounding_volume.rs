//! GeoLibre tool: minimum enclosing volume around 3D features.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Minimum Bounding Volume* (3D
//! Analyst). GeoLibre's 3D family is now substantial — `buffer_3d`, `near_3d`,
//! `adjust_3d_z`, `simplify_3d_line`, `calculate_missing_z_values`,
//! `interpolate_shape` — but every bounding-geometry tool in either registry is
//! strictly 2D: the bundled `minimum_bounding_box`, `minimum_bounding_circle`
//! and `minimum_bounding_envelope` all discard Z. There was no way to ask for
//! the *volume* occupied by a point-cloud cluster, a plume, a building's lidar
//! returns or a set of 3D tracks.
//!
//! The vector model has no multipatch type, so the enclosing solid is written
//! as its triangulated faces — one Z-valued polygon per face, all sharing an
//! `MBV_ID`. That is renderable, measurable, and consistent with how the rest
//! of the crate emits 3D results.
//!
//! **Scope for v1:** `convex_hull`, `sphere` and `envelope` are implemented.
//! `concave_hull` (an alpha shape over a 3D Delaunay tetrahedralisation) is
//! deliberately deferred rather than approximated — a wrong concave hull is
//! worse than an absent one, because its volume looks plausible.

use std::collections::BTreeMap;
use std::collections::HashMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct MinimumBoundingVolumeTool;

impl Tool for MinimumBoundingVolumeTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "minimum_bounding_volume",
            display_name: "Minimum Bounding Volume",
            summary: "Compute the minimum enclosing volume (convex hull, sphere or envelope) around 3D features, per feature, per group, or over all input, emitting the solid as Z-valued triangular faces with volume and surface-area attributes. Like ArcGIS Minimum Bounding Volume.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input 3D point, multipoint or line features.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output path for the volume faces. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "z_value",
                    description: "'geometry' (default) to take Z from the geometry, or the name of a numeric field supplying Z.",
                    required: false,
                },
                ToolParamSpec {
                    name: "geometry_type",
                    description: "'convex_hull' (default), 'sphere' or 'envelope'. ('concave_hull' is not implemented in this version.)",
                    required: false,
                },
                ToolParamSpec {
                    name: "group_option",
                    description: "'none' (one volume per feature, default), 'all' (a single volume over all input) or 'list' (one per 'group_field' value).",
                    required: false,
                },
                ToolParamSpec {
                    name: "group_field",
                    description: "Field defining groups when 'group_option' is 'list'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "mbv_fields",
                    description: "Append MBV_VOLUME and MBV_AREA fields (default false).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_shape(args)?;
        let group = parse_group(args)?;
        if group == GroupOption::List && parse_optional_str(args, "group_field")?.is_none() {
            return Err(ToolError::Validation(
                "'group_field' is required when 'group_option' is 'list'".to_string(),
            ));
        }
        parse_optional_bool(args, "mbv_fields")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let shape = parse_shape(args)?;
        let group = parse_group(args)?;
        let group_field = parse_optional_str(args, "group_field")?;
        let want_fields = parse_optional_bool(args, "mbv_fields")?.unwrap_or(false);
        let z_value = parse_optional_str(args, "z_value")?.unwrap_or("geometry");

        let layer = load_input_layer(input)?;
        if layer.features.is_empty() {
            return Err(ToolError::Execution("input has no features".to_string()));
        }
        let z_idx = if z_value.eq_ignore_ascii_case("geometry") {
            None
        } else {
            Some(layer.schema.field_index(z_value).ok_or_else(|| {
                ToolError::Validation(format!("z_value field '{z_value}' not found"))
            })?)
        };
        let g_idx = match (group, group_field) {
            (GroupOption::List, Some(f)) => Some(layer.schema.field_index(f).ok_or_else(|| {
                ToolError::Validation(format!("group_field '{f}' not found"))
            })?),
            _ => None,
        };

        // Bucket the point cloud according to the grouping option.
        let mut labels: Vec<String> = Vec::new();
        let mut buckets: Vec<Vec<P3>> = Vec::new();
        let mut pos: HashMap<String, usize> = HashMap::new();
        for (fid, feat) in layer.iter().enumerate() {
            let Some(g) = &feat.geometry else { continue };
            let zoverride = z_idx.and_then(|i| feat.attributes.get(i).and_then(FieldValue::as_f64));
            let pts: Vec<P3> = g
                .all_coords()
                .iter()
                .map(|c| [c.x, c.y, zoverride.unwrap_or_else(|| c.z.unwrap_or(0.0))])
                .collect();
            if pts.is_empty() {
                continue;
            }
            let key = match group {
                GroupOption::None => fid.to_string(),
                GroupOption::All => "ALL".to_string(),
                GroupOption::List => key_of(feat.attributes.get(g_idx.expect("resolved above"))),
            };
            let b = *pos.entry(key.clone()).or_insert_with(|| {
                labels.push(key.clone());
                buckets.push(Vec::new());
                buckets.len() - 1
            });
            buckets[b].extend(pts);
        }
        if buckets.is_empty() {
            return Err(ToolError::Execution(
                "no usable geometry found in the input".to_string(),
            ));
        }
        ctx.progress
            .info(&format!("building {} bounding volume(s)", buckets.len()));

        let mut out = Layer::new("minimum_bounding_volume").with_geom_type(GeometryType::Polygon);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("MBV_ID", FieldType::Text));
        out.add_field(FieldDef::new("FACE_ID", FieldType::Integer));
        out.add_field(FieldDef::new("POINT_COUNT", FieldType::Integer));
        if want_fields {
            out.add_field(FieldDef::new("MBV_VOLUME", FieldType::Float));
            out.add_field(FieldDef::new("MBV_AREA", FieldType::Float));
        }

        let mut total_volume = 0.0;
        let mut solids = 0usize;
        let mut degenerate = 0usize;
        let mut per_volume: Vec<Value> = Vec::new();

        for (label, pts) in labels.iter().zip(buckets.iter()) {
            let Some(solid) = build_solid(shape, pts) else {
                degenerate += 1;
                continue;
            };
            solids += 1;
            total_volume += solid.volume;
            per_volume.push(json!({
                "mbv_id": label,
                "point_count": pts.len(),
                "volume": solid.volume,
                "area": solid.area,
                "face_count": solid.faces.len(),
            }));

            for (fi, tri) in solid.faces.iter().enumerate() {
                let ring: Vec<Coord> = tri
                    .iter()
                    .map(|p| Coord::xyz(p[0], p[1], p[2]))
                    .collect();
                let mut attrs = vec![
                    ("MBV_ID", FieldValue::Text(label.clone())),
                    ("FACE_ID", FieldValue::Integer(fi as i64)),
                    ("POINT_COUNT", FieldValue::Integer(pts.len() as i64)),
                ];
                if want_fields {
                    attrs.push(("MBV_VOLUME", FieldValue::Float(solid.volume)));
                    attrs.push(("MBV_AREA", FieldValue::Float(solid.area)));
                }
                out.add_feature(Some(Geometry::polygon(ring, vec![])), &attrs)
                    .map_err(|e| ToolError::Execution(format!("failed adding face: {e}")))?;
            }
        }
        if solids == 0 {
            return Err(ToolError::Execution(
                "every group was degenerate (fewer than 4 non-coplanar points); \
                 a bounding volume is not defined"
                    .to_string(),
            ));
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("volume_count".to_string(), json!(solids));
        outputs.insert("degenerate_count".to_string(), json!(degenerate));
        outputs.insert("total_volume".to_string(), json!(total_volume));
        outputs.insert("geometry_type".to_string(), json!(shape.name()));
        outputs.insert("volumes".to_string(), json!(per_volume));
        Ok(ToolRunResult { outputs })
    }
}

// ── Solids ──────────────────────────────────────────────────────────────────

type P3 = [f64; 3];

struct Solid {
    faces: Vec<[P3; 3]>,
    volume: f64,
    area: f64,
}

fn build_solid(shape: Shape, pts: &[P3]) -> Option<Solid> {
    match shape {
        Shape::Envelope => envelope(pts),
        Shape::Sphere => sphere(pts),
        Shape::ConvexHull => convex_hull(pts),
    }
}

/// Axis-aligned bounding box as 12 triangles.
fn envelope(pts: &[P3]) -> Option<Solid> {
    if pts.is_empty() {
        return None;
    }
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in pts {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let (dx, dy, dz) = (hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
    let corner = |i: usize| -> P3 {
        [
            if i & 1 == 0 { lo[0] } else { hi[0] },
            if i & 2 == 0 { lo[1] } else { hi[1] },
            if i & 4 == 0 { lo[2] } else { hi[2] },
        ]
    };
    // Six quads, each split into two triangles, wound outward.
    const QUADS: [[usize; 4]; 6] = [
        [0, 2, 3, 1], // z = lo
        [4, 5, 7, 6], // z = hi
        [0, 1, 5, 4], // y = lo
        [2, 6, 7, 3], // y = hi
        [0, 4, 6, 2], // x = lo
        [1, 3, 7, 5], // x = hi
    ];
    let mut faces = Vec::with_capacity(12);
    for q in QUADS {
        let c: Vec<P3> = q.iter().map(|&i| corner(i)).collect();
        faces.push([c[0], c[1], c[2]]);
        faces.push([c[0], c[2], c[3]]);
    }
    Some(Solid {
        faces,
        volume: dx * dy * dz,
        area: 2.0 * (dx * dy + dy * dz + dz * dx),
    })
}

/// Minimum enclosing sphere: Ritter's approximation followed by expansion
/// passes until every point is inside, then tessellated to triangles.
fn sphere(pts: &[P3]) -> Option<Solid> {
    if pts.is_empty() {
        return None;
    }
    // Ritter seed: the pair of points farthest apart along each axis extreme.
    let mut min_i = [0usize; 3];
    let mut max_i = [0usize; 3];
    for (i, p) in pts.iter().enumerate() {
        for k in 0..3 {
            if p[k] < pts[min_i[k]][k] {
                min_i[k] = i;
            }
            if p[k] > pts[max_i[k]][k] {
                max_i[k] = i;
            }
        }
    }
    // Compare all pairs among the six axis extremes, not just same-axis pairs:
    // on a cube the same-axis pairing picks an edge (length s) instead of a
    // face diagonal, which seeds a needlessly large sphere.
    let extremes: Vec<usize> = {
        let mut v: Vec<usize> = min_i.iter().chain(max_i.iter()).copied().collect();
        v.sort_unstable();
        v.dedup();
        v
    };
    let mut best = (extremes[0], extremes[0], -1.0_f64);
    for (i, &a) in extremes.iter().enumerate() {
        for &b in &extremes[i + 1..] {
            let d = dist(&pts[a], &pts[b]);
            if d > best.2 {
                best = (a, b, d);
            }
        }
    }
    let (a, b) = (pts[best.0], pts[best.1]);
    let mut center = [
        (a[0] + b[0]) / 2.0,
        (a[1] + b[1]) / 2.0,
        (a[2] + b[2]) / 2.0,
    ];

    // Badoiu-Clarkson refinement: repeatedly step the centre a shrinking
    // fraction of the way toward the currently farthest point. This converges
    // on the 1-centre, where Ritter's expansion alone stalls well wide of it
    // (~35% excess radius on a cube).
    for i in 0..128 {
        let Some(far) = pts
            .iter()
            .max_by(|p, q| dist(&center, p).total_cmp(&dist(&center, q)))
        else {
            break;
        };
        let step = 1.0 / (i as f64 + 2.0);
        for k in 0..3 {
            center[k] += (far[k] - center[k]) * step;
        }
    }

    // Enclosure is guaranteed by construction: take the radius from the final
    // centre rather than carrying one forward from the refinement.
    let radius = pts
        .iter()
        .map(|p| dist(&center, p))
        .fold(0.0_f64, f64::max);
    if radius <= 0.0 || !radius.is_finite() {
        return None;
    }

    // UV tessellation; 12 stacks x 24 sectors is smooth enough to render and
    // cheap enough not to bloat the output.
    const STACKS: usize = 12;
    const SECTORS: usize = 24;
    let at = |i: usize, j: usize| -> P3 {
        let phi = std::f64::consts::PI * i as f64 / STACKS as f64;
        let theta = std::f64::consts::TAU * j as f64 / SECTORS as f64;
        [
            center[0] + radius * phi.sin() * theta.cos(),
            center[1] + radius * phi.sin() * theta.sin(),
            center[2] + radius * phi.cos(),
        ]
    };
    let mut faces = Vec::new();
    for i in 0..STACKS {
        for j in 0..SECTORS {
            let (p00, p01, p10, p11) = (at(i, j), at(i, j + 1), at(i + 1, j), at(i + 1, j + 1));
            if i > 0 {
                faces.push([p00, p10, p01]);
            }
            if i + 1 < STACKS {
                faces.push([p01, p10, p11]);
            }
        }
    }
    Some(Solid {
        faces,
        volume: 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3),
        area: 4.0 * std::f64::consts::PI * radius * radius,
    })
}

/// Incremental 3D convex hull (quickhull-style): seed a tetrahedron, then for
/// each remaining point delete the faces it can see and re-triangulate the
/// horizon.
fn convex_hull(pts: &[P3]) -> Option<Solid> {
    let n = pts.len();
    if n < 4 {
        return None;
    }
    // Seed: extreme point, farthest from it, farthest from that line, farthest
    // from that plane. Any failure means the cloud is degenerate (coplanar).
    let i0 = (0..n).min_by(|&a, &b| pts[a][0].total_cmp(&pts[b][0]))?;
    let i1 = (0..n).max_by(|&a, &b| dist(&pts[i0], &pts[a]).total_cmp(&dist(&pts[i0], &pts[b])))?;
    if dist(&pts[i0], &pts[i1]) < 1e-12 {
        return None;
    }
    let i2 = (0..n).max_by(|&a, &b| {
        line_dist(&pts[i0], &pts[i1], &pts[a]).total_cmp(&line_dist(&pts[i0], &pts[i1], &pts[b]))
    })?;
    if line_dist(&pts[i0], &pts[i1], &pts[i2]) < 1e-12 {
        return None;
    }
    let nrm = cross(&sub(&pts[i1], &pts[i0]), &sub(&pts[i2], &pts[i0]));
    let i3 = (0..n).max_by(|&a, &b| {
        dot(&nrm, &sub(&pts[a], &pts[i0]))
            .abs()
            .total_cmp(&dot(&nrm, &sub(&pts[b], &pts[i0])).abs())
    })?;
    if dot(&nrm, &sub(&pts[i3], &pts[i0])).abs() < 1e-12 {
        return None;
    }

    let interior = [
        (pts[i0][0] + pts[i1][0] + pts[i2][0] + pts[i3][0]) / 4.0,
        (pts[i0][1] + pts[i1][1] + pts[i2][1] + pts[i3][1]) / 4.0,
        (pts[i0][2] + pts[i1][2] + pts[i2][2] + pts[i3][2]) / 4.0,
    ];
    let mut faces: Vec<[usize; 3]> = vec![
        [i0, i1, i2],
        [i0, i1, i3],
        [i0, i2, i3],
        [i1, i2, i3],
    ];
    for f in &mut faces {
        orient(f, pts, &interior);
    }

    for (pi, p) in pts.iter().enumerate() {
        if pi == i0 || pi == i1 || pi == i2 || pi == i3 {
            continue;
        }
        let visible: Vec<usize> = (0..faces.len())
            .filter(|&fi| face_sees(&faces[fi], pts, p))
            .collect();
        if visible.is_empty() {
            continue;
        }
        // Horizon = edges belonging to exactly one visible face.
        let mut edge_count: HashMap<(usize, usize), i32> = HashMap::new();
        for &fi in &visible {
            let f = faces[fi];
            for e in [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])] {
                let key = if e.0 < e.1 { e } else { (e.1, e.0) };
                *edge_count.entry(key).or_insert(0) += 1;
            }
        }
        let horizon: Vec<(usize, usize)> = edge_count
            .into_iter()
            .filter(|&(_, c)| c == 1)
            .map(|(e, _)| e)
            .collect();

        let vis: std::collections::HashSet<usize> = visible.into_iter().collect();
        let kept: Vec<[usize; 3]> = faces
            .iter()
            .enumerate()
            .filter(|(i, _)| !vis.contains(i))
            .map(|(_, f)| *f)
            .collect();
        faces = kept;
        for (a, b) in horizon {
            let mut f = [a, b, pi];
            orient(&mut f, pts, &interior);
            faces.push(f);
        }
    }
    if faces.len() < 4 {
        return None;
    }

    // Volume via the divergence theorem, area as the triangle-area sum.
    let mut volume = 0.0;
    let mut area = 0.0;
    let mut tris = Vec::with_capacity(faces.len());
    for f in &faces {
        let (a, b, c) = (pts[f[0]], pts[f[1]], pts[f[2]]);
        volume += dot(&a, &cross(&b, &c)) / 6.0;
        area += norm(&cross(&sub(&b, &a), &sub(&c, &a))) / 2.0;
        tris.push([a, b, c]);
    }
    Some(Solid {
        faces: tris,
        volume: volume.abs(),
        area,
    })
}

/// Flips a face so its normal points away from the interior reference point.
fn orient(f: &mut [usize; 3], pts: &[P3], interior: &P3) {
    let (a, b, c) = (pts[f[0]], pts[f[1]], pts[f[2]]);
    let n = cross(&sub(&b, &a), &sub(&c, &a));
    if dot(&n, &sub(interior, &a)) > 0.0 {
        f.swap(1, 2);
    }
}

fn face_sees(f: &[usize; 3], pts: &[P3], p: &P3) -> bool {
    let (a, b, c) = (pts[f[0]], pts[f[1]], pts[f[2]]);
    let n = cross(&sub(&b, &a), &sub(&c, &a));
    let len = norm(&n);
    if len <= 0.0 {
        return false;
    }
    dot(&n, &sub(p, &a)) / len > 1e-10
}

fn sub(a: &P3, b: &P3) -> P3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross(a: &P3, b: &P3) -> P3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: &P3, b: &P3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(a: &P3) -> f64 {
    dot(a, a).sqrt()
}
fn dist(a: &P3, b: &P3) -> f64 {
    norm(&sub(a, b))
}
fn line_dist(a: &P3, b: &P3, p: &P3) -> f64 {
    let ab = sub(b, a);
    let len = norm(&ab);
    if len <= 0.0 {
        return dist(a, p);
    }
    norm(&cross(&ab, &sub(p, a))) / len
}

// ── Params ──────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Shape {
    ConvexHull,
    Sphere,
    Envelope,
}

impl Shape {
    fn name(self) -> &'static str {
        match self {
            Shape::ConvexHull => "convex_hull",
            Shape::Sphere => "sphere",
            Shape::Envelope => "envelope",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum GroupOption {
    None,
    All,
    List,
}

fn parse_shape(args: &ToolArgs) -> Result<Shape, ToolError> {
    match args
        .get("geometry_type")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        None | Some("") | Some("convex_hull") => Ok(Shape::ConvexHull),
        Some("sphere") => Ok(Shape::Sphere),
        Some("envelope") => Ok(Shape::Envelope),
        Some("concave_hull") => Err(ToolError::Validation(
            "'concave_hull' is not implemented in this version; use 'convex_hull', \
             'sphere' or 'envelope'"
                .to_string(),
        )),
        Some(o) => Err(ToolError::Validation(format!(
            "'geometry_type' must be 'convex_hull', 'sphere' or 'envelope', got '{o}'"
        ))),
    }
}

fn parse_group(args: &ToolArgs) -> Result<GroupOption, ToolError> {
    match args
        .get("group_option")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        None | Some("") | Some("none") => Ok(GroupOption::None),
        Some("all") => Ok(GroupOption::All),
        Some("list") => Ok(GroupOption::List),
        Some(o) => Err(ToolError::Validation(format!(
            "'group_option' must be 'none', 'all' or 'list', got '{o}'"
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

fn key_of(v: Option<&FieldValue>) -> String {
    match v {
        None | Some(FieldValue::Null) => "NULL".to_string(),
        Some(FieldValue::Integer(i)) => i.to_string(),
        Some(FieldValue::Float(f)) => format!("{f}"),
        Some(FieldValue::Text(s)) => s.clone(),
        Some(FieldValue::Boolean(b)) => b.to_string(),
        Some(FieldValue::Date(s)) | Some(FieldValue::DateTime(s)) => s.clone(),
        Some(FieldValue::Blob(b)) => format!("blob[{}]", b.len()),
    }
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

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn cloud(pts: &[(f64, f64, f64)], groups: &[&str]) -> String {
        let mut l = Layer::new("c")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        if !groups.is_empty() {
            l.add_field(FieldDef::new("grp", FieldType::Text));
        }
        for (i, (x, y, z)) in pts.iter().enumerate() {
            let a: Vec<(&str, FieldValue)> = if groups.is_empty() {
                vec![]
            } else {
                vec![("grp", FieldValue::Text(groups[i % groups.len()].to_string()))]
            };
            l.add_feature(Some(Geometry::point_z(*x, *y, *z)), &a).unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    /// The 8 corners of an axis-aligned cube of side `s`, plus an interior point.
    fn cube(s: f64) -> Vec<(f64, f64, f64)> {
        let mut v = Vec::new();
        for i in 0..8 {
            v.push((
                if i & 1 == 0 { 0.0 } else { s },
                if i & 2 == 0 { 0.0 } else { s },
                if i & 4 == 0 { 0.0 } else { s },
            ));
        }
        v.push((s / 2.0, s / 2.0, s / 2.0));
        v
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = MinimumBoundingVolumeTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    #[test]
    fn convex_hull_of_a_cube_has_the_cube_volume() {
        // 9 points (8 corners + centre) grouped into one hull of side 10.
        let input = cloud(&cube(10.0), &[]);
        let (out, _l) = run(json!({
            "input": input, "group_option": "all", "geometry_type": "convex_hull"
        }));
        let v = out.outputs["total_volume"].as_f64().unwrap();
        assert!((v - 1000.0).abs() < 1e-6, "hull volume was {v}");
        let area = out.outputs["volumes"].as_array().unwrap()[0]["area"]
            .as_f64()
            .unwrap();
        assert!((area - 600.0).abs() < 1e-6, "hull area was {area}");
    }

    #[test]
    fn envelope_matches_the_axis_aligned_extent() {
        // Points spanning 2 x 3 x 4 -> volume 24, area 2*(6+12+8) = 52.
        let input = cloud(
            &[(0.0, 0.0, 0.0), (2.0, 3.0, 4.0), (1.0, 1.0, 1.0)],
            &[],
        );
        let (out, _l) = run(json!({
            "input": input, "group_option": "all", "geometry_type": "envelope"
        }));
        assert!((out.outputs["total_volume"].as_f64().unwrap() - 24.0).abs() < 1e-9);
        let area = out.outputs["volumes"].as_array().unwrap()[0]["area"]
            .as_f64()
            .unwrap();
        assert!((area - 52.0).abs() < 1e-9);
    }

    #[test]
    fn sphere_encloses_every_input_point() {
        let pts = cube(10.0);
        let input = cloud(&pts, &[]);
        let (out, _l) = run(json!({
            "input": input, "group_option": "all", "geometry_type": "sphere"
        }));
        // A sphere around the cube must be at least the cube's circumsphere:
        // r = sqrt(3)*10/2 -> V = 4/3 pi r^3 ~= 2721.
        let v = out.outputs["total_volume"].as_f64().unwrap();
        let r_min = 3.0_f64.sqrt() * 10.0 / 2.0;
        let v_min = 4.0 / 3.0 * std::f64::consts::PI * r_min.powi(3);
        assert!(v >= v_min - 1e-6, "sphere volume {v} < circumsphere {v_min}");
        // ...and not wildly larger than it.
        assert!(v < v_min * 1.3, "sphere volume {v} too loose");
    }

    #[test]
    fn hull_volume_never_exceeds_the_envelope() {
        // Fundamental ordering: convex hull <= envelope for the same cloud.
        let pts = [
            (0.0, 0.0, 0.0),
            (10.0, 0.0, 0.0),
            (0.0, 10.0, 0.0),
            (0.0, 0.0, 10.0),
        ];
        let input = cloud(&pts, &[]);
        let (hull, _l) = run(json!({
            "input": input.clone(), "group_option": "all", "geometry_type": "convex_hull"
        }));
        let (env, _l) = run(json!({
            "input": input, "group_option": "all", "geometry_type": "envelope"
        }));
        let (h, e) = (
            hull.outputs["total_volume"].as_f64().unwrap(),
            env.outputs["total_volume"].as_f64().unwrap(),
        );
        // The tetrahedron is exactly 1/6 of its bounding cube.
        assert!((h - 1000.0 / 6.0).abs() < 1e-6, "tetra volume {h}");
        assert!(h < e);
    }

    #[test]
    fn group_field_builds_one_volume_per_group() {
        let mut pts = cube(10.0);
        pts.extend(cube(10.0).iter().map(|(x, y, z)| (x + 100.0, *y, *z)));
        let groups: Vec<&str> = std::iter::repeat_n("a", 9)
            .chain(std::iter::repeat_n("b", 9))
            .collect();
        let mut l = Layer::new("c")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("grp", FieldType::Text));
        for (p, g) in pts.iter().zip(groups.iter()) {
            l.add_feature(
                Some(Geometry::point_z(p.0, p.1, p.2)),
                &[("grp", FieldValue::Text((*g).to_string()))],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        let input = wbvector::memory_store::make_vector_memory_path(&id);

        let (out, _l) = run(json!({
            "input": input, "group_option": "list", "group_field": "grp"
        }));
        assert_eq!(out.outputs["volume_count"], json!(2));
        // Two separate 10-cubes, not one hull spanning the 100-unit gap.
        assert!((out.outputs["total_volume"].as_f64().unwrap() - 2000.0).abs() < 1e-6);
    }

    #[test]
    fn faces_carry_z_and_share_an_mbv_id() {
        let input = cloud(&cube(4.0), &[]);
        let (_o, layer) = run(json!({
            "input": input, "group_option": "all", "mbv_fields": true
        }));
        assert!(!layer.features.is_empty());
        let id = layer.schema.field_index("MBV_ID").unwrap();
        let vol = layer.schema.field_index("MBV_VOLUME").unwrap();
        let first = layer.features[0].attributes[id].as_str().unwrap().to_string();
        assert!(layer
            .iter()
            .all(|f| f.attributes[id].as_str() == Some(first.as_str())));
        assert!(layer.iter().all(|f| f.attributes[vol].as_f64().is_some()));
        for f in layer.iter() {
            match f.geometry.as_ref().unwrap() {
                Geometry::Polygon { exterior, .. } => {
                    assert_eq!(exterior.0.len(), 3);
                    assert!(exterior.0.iter().all(|c| c.z.is_some()));
                }
                other => panic!("unexpected geometry {other:?}"),
            }
        }
    }

    #[test]
    fn coplanar_cloud_is_reported_as_degenerate() {
        // All z = 0: no enclosing volume exists.
        let input = cloud(
            &[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.0)],
            &[],
        );
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input, "group_option": "all", "geometry_type": "convex_hull"
        }))
        .unwrap();
        assert!(MinimumBoundingVolumeTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            MinimumBoundingVolumeTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "p.shp", "geometry_type": "concave_hull" })).is_err());
        assert!(bad(json!({ "input": "p.shp", "geometry_type": "blob" })).is_err());
        assert!(bad(json!({ "input": "p.shp", "group_option": "list" })).is_err());
        assert!(bad(json!({ "input": "p.shp" })).is_ok());
    }
}
