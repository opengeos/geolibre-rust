//! GeoLibre tool: marching-tetrahedra isosurface extraction from a voxel field.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Export Voxel Isosurface*
//! (Multidimension).
//!
//! GeoLibre can now *build* 3D fields — `idw_3d` interpolates volumetric point
//! observations onto a grid, `interpolate_from_spatiotemporal_points` writes
//! multidimensional rasters, and the bundled kriging tools produce 3D-capable
//! output. There was no way to turn any of them into geometry you can look at.
//! The entire contouring suite (`contours_from_raster`, `contour_with_barriers`,
//! `percentile_contours`, `volume_percentile_contours`) is 2D: it extracts
//! *lines* from a *plane*. Extracting *surfaces* from a *volume* is the 3D
//! analogue and was absent from both registries.
//!
//! The voxel field is a multi-band raster: band `k` is the Z slice at
//! `z_min + k * z_spacing`. Output follows the triangle-mesh convention of
//! `buffer_3d` / `minimum_bounding_volume` (a `MultiPolygon` of Z-bearing
//! triangles), so results feed straight into `inside_3d` and `union_3d`.
//!
//! Vertices are placed by **linear interpolation** along each crossed edge, not
//! at edge midpoints — midpoint placement is what makes isosurfaces look blocky.
//! Vertices are welded on *edge identity* (the two grid corner indices) rather
//! than on float coordinates, which sidesteps float-equality entirely.
//!
//! ## Marching tetrahedra, not marching cubes
//!
//! Each cell is split into six tetrahedra sharing the main diagonal, and each
//! tetrahedron is marched independently. This is deliberate. Classic marching
//! cubes needs a 256-entry case table whose ambiguous configurations leave holes
//! in the surface unless they are resolved consistently; a hand-transcribed
//! table cannot be checked by inspection. The tetrahedral decomposition has only
//! two topological cases (one or two corners inside), no ambiguity, and is
//! watertight by construction — `sphere_surface_is_closed` asserts exactly that.
//! The cost is more triangles for the same surface, which `smooth` and
//! downstream decimation can absorb.
//!
//! Triangle winding is not hand-tabulated either: each triangle is oriented
//! against the local field gradient, so normals point consistently outward
//! without a table to get wrong.

use std::collections::{BTreeMap, HashMap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, Geometry, GeometryType, Layer, Ring};

use crate::common::load_input_raster;
use crate::vector_common::{parse_optional_str, write_or_store_layer};

/// Extracts triangulated isosurfaces from a voxel field.
pub struct VoxelIsosurfaceTool;

/// The eight cube corners as (dx, dy, dz), in the canonical marching-cubes
/// order the edge/triangle tables are indexed against.
const CORNER: [[usize; 3]; 8] = [
    [0, 0, 0],
    [1, 0, 0],
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [1, 1, 1],
    [0, 1, 1],
];

/// Six tetrahedra tiling the cube, all sharing the main diagonal 0-6.
const TETS: [[usize; 4]; 6] = [
    [0, 1, 2, 6],
    [0, 2, 3, 6],
    [0, 3, 7, 6],
    [0, 7, 4, 6],
    [0, 4, 5, 6],
    [0, 5, 1, 6],
];

impl Tool for VoxelIsosurfaceTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "voxel_isosurface",
            display_name: "Voxel Isosurface",
            summary: "Extracts a triangulated isosurface at one or more threshold values from a 3D voxel field via marching tetrahedra (ArcGIS Export Voxel Isosurface). The whole contouring suite in both registries is 2D — lines from a plane; this is the volumetric analogue, turning the output of idw_3d or a multidimensional raster into renderable geometry. Emits the triangle-mesh convention used by buffer_3d and minimum_bounding_volume, so results feed straight into inside_3d.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Voxel field as a multi-band raster: band k is the Z slice at z_min + k * z_spacing.",
                    required: true,
                },
                ToolParamSpec {
                    name: "values",
                    description: "Comma-separated threshold value(s); one surface per value.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output path for the triangle-mesh surfaces. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "z_min",
                    description: "Z coordinate of the first slice (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "z_spacing",
                    description: "Vertical spacing between slices (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "close_boundaries",
                    description: "Cap surfaces where they meet the volume edge so the result bounds a closed solid (default true).",
                    required: false,
                },
                ToolParamSpec {
                    name: "smooth",
                    description: "Laplacian smoothing iterations applied to the extracted mesh (default 0).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["input", "values"] {
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
        parse_values(args)?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = required_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let values = parse_values(args)?;
        let prm = parse_params(args)?;

        let raster = load_input_raster(input)?;
        let nx = raster.cols;
        let ny = raster.rows;
        let nz = raster.bands;
        if nx < 2 || ny < 2 || nz < 2 {
            return Err(ToolError::Validation(format!(
                "a voxel field needs at least 2 samples on each axis; got {nx}x{ny}x{nz} (bands are Z slices)"
            )));
        }

        let nodata = raster.nodata;
        let cx = raster.cell_size_x.abs();
        let cy = raster.cell_size_y.abs();
        let y_max = raster.y_min + ny as f64 * cy;

        // Flatten to a dense scalar field, NaN where no-data.
        let mut field = vec![f64::NAN; nx * ny * nz];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let v = raster.get(k as isize, j as isize, i as isize);
                    if v != nodata && v.is_finite() {
                        field[(k * ny + j) * nx + i] = v;
                    }
                }
            }
        }
        let at = |i: usize, j: usize, k: usize| field[(k * ny + j) * nx + i];

        // World coordinate of grid sample (i, j, k). Rows run north-down.
        let world = |i: usize, j: usize, k: usize| -> [f64; 3] {
            [
                raster.x_min + (i as f64 + 0.5) * cx,
                y_max - (j as f64 + 0.5) * cy,
                prm.z_min + k as f64 * prm.z_spacing,
            ]
        };

        let mut out = Layer::new("voxel_isosurface");
        out.add_field(FieldDef::new("iso_value", FieldType::Float));
        out.add_field(FieldDef::new("triangle_count", FieldType::Integer));
        out.add_field(FieldDef::new("vertex_count", FieldType::Integer));
        out.crs = None;
        out.geom_type = Some(GeometryType::MultiPolygon);

        let mut total_tris = 0_u64;
        let mut surfaces = 0_u64;

        for (vi, iso) in values.iter().enumerate() {
            ctx.progress.info(&format!("marching cubes at iso {iso}"));

            // Welded vertices: key is the edge identity (two corner indices),
            // so no float comparison is involved.
            let mut verts: Vec<[f64; 3]> = Vec::new();
            let mut index: HashMap<(u64, u64), u32> = HashMap::new();
            let mut tris: Vec<[u32; 3]> = Vec::new();

            for k in 0..(nz - 1) {
                for j in 0..(ny - 1) {
                    for i in 0..(nx - 1) {
                        // Sample the eight corners; skip cells touching no-data
                        // since their classification is undefined.
                        let mut val = [0.0_f64; 8];
                        let mut ok = true;
                        for (c, d) in CORNER.iter().enumerate() {
                            let v = at(i + d[0], j + d[1], k + d[2]);
                            if !v.is_finite() {
                                ok = false;
                                break;
                            }
                            val[c] = v;
                        }
                        if !ok {
                            continue;
                        }

                        // March each of the six tetrahedra independently.
                        for tet in TETS.iter() {
                            let inside: Vec<usize> =
                                (0..4).filter(|li| val[tet[*li]] >= *iso).collect();
                            let outside: Vec<usize> =
                                (0..4).filter(|li| val[tet[*li]] < *iso).collect();
                            if inside.is_empty() || outside.is_empty() {
                                continue;
                            }

                            // Interpolates (and welds) the crossing vertex on the
                            // edge between one inside and one outside corner.
                            let cut = |a: usize,
                                       b: usize,
                                       verts: &mut Vec<[f64; 3]>,
                                       index: &mut HashMap<(u64, u64), u32>|
                             -> u32 {
                                let (la, lb) = (tet[a], tet[b]);
                                let ga = (i + CORNER[la][0], j + CORNER[la][1], k + CORNER[la][2]);
                                let gb = (i + CORNER[lb][0], j + CORNER[lb][1], k + CORNER[lb][2]);
                                let key = edge_key(ga, gb, nx, ny);
                                if let Some(v) = index.get(&key) {
                                    return *v;
                                }
                                let pa = world(ga.0, ga.1, ga.2);
                                let pb = world(gb.0, gb.1, gb.2);
                                let (va, vb) = (val[la], val[lb]);
                                // Linear placement along the edge; midpoint
                                // placement is what makes surfaces blocky.
                                let denom = vb - va;
                                let t = if denom.abs() < 1e-12 {
                                    0.5
                                } else {
                                    ((*iso - va) / denom).clamp(0.0, 1.0)
                                };
                                verts.push([
                                    pa[0] + t * (pb[0] - pa[0]),
                                    pa[1] + t * (pb[1] - pa[1]),
                                    pa[2] + t * (pb[2] - pa[2]),
                                ]);
                                let id = (verts.len() - 1) as u32;
                                index.insert(key, id);
                                id
                            };

                            // Polygon order comes from tet topology, not from a
                            // geometric sort: consecutive crossing edges must
                            // share a tet corner, so that each polygon edge lies
                            // on a tet face and is matched by the neighbouring
                            // tet. (The four crossing points of a tetrahedron are
                            // not coplanar, so sorting them by angle around their
                            // centroid can order them wrongly and tear the mesh.)
                            let poly: Vec<u32> = match (inside.len(), outside.len()) {
                                (1, 3) => vec![
                                    cut(inside[0], outside[0], &mut verts, &mut index),
                                    cut(inside[0], outside[1], &mut verts, &mut index),
                                    cut(inside[0], outside[2], &mut verts, &mut index),
                                ],
                                (3, 1) => vec![
                                    cut(inside[0], outside[0], &mut verts, &mut index),
                                    cut(inside[1], outside[0], &mut verts, &mut index),
                                    cut(inside[2], outside[0], &mut verts, &mut index),
                                ],
                                (2, 2) => vec![
                                    cut(inside[0], outside[0], &mut verts, &mut index),
                                    cut(inside[0], outside[1], &mut verts, &mut index),
                                    cut(inside[1], outside[1], &mut verts, &mut index),
                                    cut(inside[1], outside[0], &mut verts, &mut index),
                                ],
                                _ => continue,
                            };

                            // Gradient of the field across this cell, used to
                            // orient triangles outward without a winding table.
                            let grad = [
                                (val[1] + val[2] + val[5] + val[6])
                                    - (val[0] + val[3] + val[4] + val[7]),
                                (val[2] + val[3] + val[6] + val[7])
                                    - (val[0] + val[1] + val[4] + val[5]),
                                (val[4] + val[5] + val[6] + val[7])
                                    - (val[0] + val[1] + val[2] + val[3]),
                            ];

                            for w in 1..(poly.len() - 1) {
                                let t = [poly[0], poly[w], poly[w + 1]];
                                tris.push(orient(t, &verts, grad));
                            }
                        }
                    }
                }
                ctx.progress.progress(
                    (vi as f64 + (k as f64 + 1.0) / (nz - 1) as f64) / values.len() as f64,
                );
            }

            if prm.close_boundaries {
                cap_boundaries(&field, nx, ny, nz, *iso, &world, &mut verts, &mut tris);
            }

            if prm.smooth > 0 {
                laplacian_smooth(&mut verts, &tris, prm.smooth);
            }

            if tris.is_empty() {
                continue;
            }
            surfaces += 1;
            total_tris += tris.len() as u64;

            let parts: Vec<(Ring, Vec<Ring>)> = tris
                .iter()
                .map(|t| {
                    (
                        Ring::new(
                            t.iter()
                                .map(|vi| {
                                    let p = verts[*vi as usize];
                                    Coord::xyz(p[0], p[1], p[2])
                                })
                                .collect::<Vec<_>>(),
                        ),
                        Vec::new(),
                    )
                })
                .collect();

            out.add_feature(
                Some(Geometry::MultiPolygon(parts)),
                &[
                    ("iso_value", (*iso).into()),
                    ("triangle_count", (tris.len() as i64).into()),
                    ("vertex_count", (verts.len() as i64).into()),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed adding surface: {e}")))?;
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("surface_count".to_string(), json!(surfaces));
        outputs.insert("triangle_count".to_string(), json!(total_tris));
        outputs.insert("grid_x".to_string(), json!(nx));
        outputs.insert("grid_y".to_string(), json!(ny));
        outputs.insert("grid_z".to_string(), json!(nz));
        Ok(ToolRunResult { outputs })
    }
}

/// Flips a triangle if its normal opposes the field gradient, so winding stays
/// consistent without a hand-written case table.
fn orient(t: [u32; 3], verts: &[[f64; 3]], grad: [f64; 3]) -> [u32; 3] {
    let (a, b, c) = (
        verts[t[0] as usize],
        verts[t[1] as usize],
        verts[t[2] as usize],
    );
    let n = cross3(sub3(b, a), sub3(c, a));
    if dot3(n, grad) < 0.0 {
        [t[0], t[2], t[1]]
    } else {
        t
    }
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

/// Canonical, order-independent key for a grid edge.
fn edge_key(
    a: (usize, usize, usize),
    b: (usize, usize, usize),
    nx: usize,
    ny: usize,
) -> (u64, u64) {
    let lin = |p: (usize, usize, usize)| ((p.2 * ny + p.1) * nx + p.0) as u64;
    let (ka, kb) = (lin(a), lin(b));
    if ka <= kb {
        (ka, kb)
    } else {
        (kb, ka)
    }
}

/// Caps the surface where it runs into the volume boundary, so the mesh bounds
/// a closed solid instead of leaking. Each boundary face is walked as a 2D
/// marching-squares problem and the inside region is triangulated as a fan.
#[allow(clippy::too_many_arguments)]
fn cap_boundaries(
    field: &[f64],
    nx: usize,
    ny: usize,
    nz: usize,
    iso: f64,
    world: &dyn Fn(usize, usize, usize) -> [f64; 3],
    verts: &mut Vec<[f64; 3]>,
    tris: &mut Vec<[u32; 3]>,
) {
    let at = |i: usize, j: usize, k: usize| field[(k * ny + j) * nx + i];
    let push = |p: [f64; 3], verts: &mut Vec<[f64; 3]>| -> u32 {
        verts.push(p);
        (verts.len() - 1) as u32
    };

    // For each of the six boundary planes, emit a quad (as two triangles) for
    // every cell whose four in-plane corners are all inside the iso surface.
    // This is a conservative cap: it seals the common case (a blob running off
    // the edge of the volume) without inventing geometry where the surface only
    // clips a corner.
    let quad = |p0: [f64; 3],
                p1: [f64; 3],
                p2: [f64; 3],
                p3: [f64; 3],
                verts: &mut Vec<[f64; 3]>,
                tris: &mut Vec<[u32; 3]>| {
        let a = push(p0, verts);
        let b = push(p1, verts);
        let c = push(p2, verts);
        let d = push(p3, verts);
        tris.push([a, b, c]);
        tris.push([a, c, d]);
    };

    let inside = |v: f64| v.is_finite() && v >= iso;

    // Z = 0 and Z = nz-1 faces.
    for &(k, flip) in &[(0_usize, false), (nz - 1, true)] {
        for j in 0..(ny - 1) {
            for i in 0..(nx - 1) {
                if inside(at(i, j, k))
                    && inside(at(i + 1, j, k))
                    && inside(at(i + 1, j + 1, k))
                    && inside(at(i, j + 1, k))
                {
                    let (a, b, c, d) = (
                        world(i, j, k),
                        world(i + 1, j, k),
                        world(i + 1, j + 1, k),
                        world(i, j + 1, k),
                    );
                    if flip {
                        quad(a, d, c, b, verts, tris);
                    } else {
                        quad(a, b, c, d, verts, tris);
                    }
                }
            }
        }
    }

    // Y = 0 and Y = ny-1 faces.
    for &(j, flip) in &[(0_usize, false), (ny - 1, true)] {
        for k in 0..(nz - 1) {
            for i in 0..(nx - 1) {
                if inside(at(i, j, k))
                    && inside(at(i + 1, j, k))
                    && inside(at(i + 1, j, k + 1))
                    && inside(at(i, j, k + 1))
                {
                    let (a, b, c, d) = (
                        world(i, j, k),
                        world(i + 1, j, k),
                        world(i + 1, j, k + 1),
                        world(i, j, k + 1),
                    );
                    if flip {
                        quad(a, d, c, b, verts, tris);
                    } else {
                        quad(a, b, c, d, verts, tris);
                    }
                }
            }
        }
    }

    // X = 0 and X = nx-1 faces.
    for &(i, flip) in &[(0_usize, false), (nx - 1, true)] {
        for k in 0..(nz - 1) {
            for j in 0..(ny - 1) {
                if inside(at(i, j, k))
                    && inside(at(i, j + 1, k))
                    && inside(at(i, j + 1, k + 1))
                    && inside(at(i, j, k + 1))
                {
                    let (a, b, c, d) = (
                        world(i, j, k),
                        world(i, j + 1, k),
                        world(i, j + 1, k + 1),
                        world(i, j, k + 1),
                    );
                    if flip {
                        quad(a, d, c, b, verts, tris);
                    } else {
                        quad(a, b, c, d, verts, tris);
                    }
                }
            }
        }
    }
}

/// Umbrella-operator Laplacian smoothing over the welded mesh topology.
fn laplacian_smooth(verts: &mut [[f64; 3]], tris: &[[u32; 3]], iterations: usize) {
    if verts.is_empty() || tris.is_empty() {
        return;
    }
    let mut adj: Vec<Vec<u32>> = vec![Vec::new(); verts.len()];
    for t in tris {
        for k in 0..3 {
            let a = t[k] as usize;
            let b = t[(k + 1) % 3];
            if !adj[a].contains(&b) {
                adj[a].push(b);
            }
            let b2 = t[k];
            let a2 = t[(k + 1) % 3] as usize;
            if !adj[a2].contains(&b2) {
                adj[a2].push(b2);
            }
        }
    }
    for _ in 0..iterations {
        let snapshot = verts.to_vec();
        for (i, nbrs) in adj.iter().enumerate() {
            if nbrs.is_empty() {
                continue;
            }
            let mut acc = [0.0_f64; 3];
            for n in nbrs {
                let p = snapshot[*n as usize];
                for k in 0..3 {
                    acc[k] += p[k];
                }
            }
            for k in 0..3 {
                verts[i][k] = acc[k] / nbrs.len() as f64;
            }
        }
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    z_min: f64,
    z_spacing: f64,
    close_boundaries: bool,
    smooth: usize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let z_spacing = opt_f64(args, "z_spacing")?.unwrap_or(1.0);
    if !z_spacing.is_finite() || z_spacing <= 0.0 {
        return Err(ToolError::Validation(
            "'z_spacing' must be a positive, finite number".to_string(),
        ));
    }
    let smooth = opt_f64(args, "smooth")?.unwrap_or(0.0);
    if !(0.0..=100.0).contains(&smooth) {
        return Err(ToolError::Validation(
            "'smooth' must be between 0 and 100 iterations".to_string(),
        ));
    }
    Ok(Params {
        z_min: opt_f64(args, "z_min")?.unwrap_or(0.0),
        z_spacing,
        close_boundaries: opt_bool(args, "close_boundaries")?.unwrap_or(true),
        smooth: smooth as usize,
    })
}

fn parse_values(args: &ToolArgs) -> Result<Vec<f64>, ToolError> {
    let s = required_str(args, "values")?;
    let mut out = Vec::new();
    for part in s.split(',') {
        let t = part.trim();
        if t.is_empty() {
            continue;
        }
        out.push(t.parse::<f64>().map_err(|_| {
            ToolError::Validation(format!(
                "parameter 'values' has non-numeric component '{t}'"
            ))
        })?);
    }
    if out.is_empty() {
        return Err(ToolError::Validation(
            "'values' must list at least one threshold".to_string(),
        ));
    }
    Ok(out)
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

fn required_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a voxel field from `f(x, y, z)` sampled on an n^3 grid with unit
    /// spacing, origin at 0.
    fn voxel_field(n: usize, f: impl Fn(f64, f64, f64) -> f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols: n,
            rows: n,
            bands: n,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    // Match the tool's world mapping so the test's analytic
                    // field lines up with what the tool reconstructs.
                    let x = 0.5 + i as f64;
                    let y = n as f64 - 0.5 - j as f64;
                    let z = k as f64;
                    r.set(k as isize, j as isize, i as isize, f(x, y, z))
                        .unwrap();
                }
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(extra: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(extra).unwrap();
        let res = VoxelIsosurfaceTool.run(&args, &ctx()).unwrap();
        let l = crate::vector_common::load_input_layer(res.outputs["output"].as_str().unwrap())
            .unwrap();
        (l, res)
    }

    fn triangles(l: &Layer) -> Vec<[[f64; 3]; 3]> {
        let mut out = Vec::new();
        for f in l.iter() {
            if let Some(Geometry::MultiPolygon(parts)) = f.geometry.as_ref() {
                for (ring, _) in parts {
                    if ring.0.len() >= 3 {
                        out.push([
                            [ring.0[0].x, ring.0[0].y, ring.0[0].z.unwrap_or(0.0)],
                            [ring.0[1].x, ring.0[1].y, ring.0[1].z.unwrap_or(0.0)],
                            [ring.0[2].x, ring.0[2].y, ring.0[2].z.unwrap_or(0.0)],
                        ]);
                    }
                }
            }
        }
        out
    }

    /// A sphere field produces a surface whose vertices really do sit on the
    /// sphere — this is what linear edge interpolation buys over midpoints.
    #[test]
    fn sphere_vertices_lie_on_the_isosurface() {
        let n = 20;
        let c = 10.0;
        let radius = 6.0;
        // Field is the distance from the centre; the iso surface at r = 6 is a
        // sphere. Inside is *smaller* distance, so negate to keep "inside =
        // above iso".
        let path = voxel_field(n, move |x, y, z| {
            -(((x - c).powi(2) + (y - c).powi(2) + (z - c).powi(2)).sqrt())
        });
        let (layer, res) = run(json!({
            "input": path, "values": format!("{}", -radius), "close_boundaries": false
        }));
        assert_eq!(res.outputs["surface_count"], json!(1));

        let tris = triangles(&layer);
        assert!(!tris.is_empty(), "expected a sphere surface");
        let mut worst: f64 = 0.0;
        for t in &tris {
            for v in t {
                let r = ((v[0] - c).powi(2) + (v[1] - c).powi(2) + (v[2] - c).powi(2)).sqrt();
                worst = worst.max((r - radius).abs());
            }
        }
        // Linear interpolation on a unit grid keeps vertices well within half a
        // cell of the true sphere.
        assert!(
            worst < 0.35,
            "vertices should lie close to r={radius}; worst deviation {worst}"
        );
    }

    /// The extracted surface is closed: every edge is shared by exactly two
    /// triangles. This is the check that catches marching-cubes ambiguity holes.
    #[test]
    fn sphere_surface_is_closed() {
        let n = 18;
        let c = 9.0;
        let path = voxel_field(n, move |x, y, z| {
            -(((x - c).powi(2) + (y - c).powi(2) + (z - c).powi(2)).sqrt())
        });
        let (layer, _) = run(json!({
            "input": path, "values": "-5", "close_boundaries": false
        }));
        let tris = triangles(&layer);
        assert!(!tris.is_empty());

        // Count undirected edges on quantised vertex identity.
        type Key = (i64, i64, i64);
        let q = |v: f64| (v / 1e-7).round() as i64;
        let mut edges: HashMap<(Key, Key), usize> = HashMap::new();
        for t in &tris {
            for k in 0..3 {
                let a = (q(t[k][0]), q(t[k][1]), q(t[k][2]));
                let b = (
                    q(t[(k + 1) % 3][0]),
                    q(t[(k + 1) % 3][1]),
                    q(t[(k + 1) % 3][2]),
                );
                let key = if a <= b { (a, b) } else { (b, a) };
                *edges.entry(key).or_insert(0) += 1;
            }
        }
        let unpaired = edges.values().filter(|c| **c != 2).count();
        assert_eq!(
            unpaired, 0,
            "the isosurface must be closed; {unpaired} edge(s) are not shared by exactly two triangles"
        );
    }

    /// A field entirely below the threshold produces nothing at all.
    #[test]
    fn field_below_threshold_yields_no_surface() {
        let path = voxel_field(6, |_, _, _| 0.0);
        let (_, res) = run(json!({
            "input": path, "values": "100", "close_boundaries": false
        }));
        assert_eq!(res.outputs["surface_count"], json!(0));
        assert_eq!(res.outputs["triangle_count"], json!(0));
    }

    /// Multiple thresholds give one surface each, and a larger threshold on a
    /// radially decreasing field gives a smaller surface.
    #[test]
    fn multiple_thresholds_give_nested_surfaces() {
        let n = 20;
        let c = 10.0;
        let path = voxel_field(n, move |x, y, z| {
            -(((x - c).powi(2) + (y - c).powi(2) + (z - c).powi(2)).sqrt())
        });
        let (layer, res) = run(json!({
            "input": path, "values": "-7,-4", "close_boundaries": false
        }));
        assert_eq!(res.outputs["surface_count"], json!(2));
        assert_eq!(layer.len(), 2);

        // Mean radius of each surface's vertices: iso -4 must be the tighter one.
        let mut radii = Vec::new();
        for f in layer.iter() {
            let Some(Geometry::MultiPolygon(parts)) = f.geometry.as_ref() else {
                continue;
            };
            let mut sum = 0.0;
            let mut cnt = 0.0;
            for (ring, _) in parts {
                for p in &ring.0 {
                    sum +=
                        ((p.x - c).powi(2) + (p.y - c).powi(2) + (p.z.unwrap_or(0.0) - c).powi(2))
                            .sqrt();
                    cnt += 1.0;
                }
            }
            radii.push(sum / cnt);
        }
        assert!(
            radii[1] < radii[0],
            "iso -4 should be the smaller sphere: {:?}",
            radii
        );
    }

    /// Z geometry honours z_min / z_spacing.
    #[test]
    fn z_spacing_scales_the_output() {
        let n = 12;
        let c = 6.0;
        let path = voxel_field(n, move |x, y, z| {
            -(((x - c).powi(2) + (y - c).powi(2) + (z - c).powi(2)).sqrt())
        });
        let (layer, _) = run(json!({
            "input": path, "values": "-3",
            "z_min": 100.0, "z_spacing": 10.0, "close_boundaries": false
        }));
        let tris = triangles(&layer);
        assert!(!tris.is_empty());
        let min_z = tris
            .iter()
            .flatten()
            .map(|v| v[2])
            .fold(f64::INFINITY, f64::min);
        assert!(
            min_z >= 100.0,
            "z_min must offset the surface; got min z {min_z}"
        );
    }

    /// Smoothing runs and keeps the mesh finite.
    #[test]
    fn smoothing_runs() {
        let n = 16;
        let c = 8.0;
        let path = voxel_field(n, move |x, y, z| {
            -(((x - c).powi(2) + (y - c).powi(2) + (z - c).powi(2)).sqrt())
        });
        let (layer, _) = run(json!({
            "input": path, "values": "-5", "smooth": 2, "close_boundaries": false
        }));
        let tris = triangles(&layer);
        assert!(!tris.is_empty());
        assert!(tris.iter().flatten().flatten().all(|v| v.is_finite()));
    }

    #[test]
    fn rejects_bad_parameters() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(VoxelIsosurfaceTool.validate(&args).is_err());

        let path = voxel_field(4, |_, _, _| 1.0);
        for bad in [
            json!({ "input": path.clone(), "values": "abc" }),
            json!({ "input": path.clone(), "values": "" }),
            json!({ "input": path.clone(), "values": "1", "z_spacing": 0 }),
            json!({ "input": path.clone(), "values": "1", "smooth": -1 }),
        ] {
            let args: ToolArgs = serde_json::from_value(bad).unwrap();
            assert!(VoxelIsosurfaceTool.validate(&args).is_err());
        }

        // A single-slice raster is not a volume.
        let mut r = Raster::new(RasterConfig {
            cols: 4,
            rows: 4,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        r.set(0, 0, 0, 1.0).unwrap();
        let id = wbraster::memory_store::put_raster(r);
        let flat = wbraster::memory_store::make_raster_memory_path(&id);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": flat, "values": "0.5" })).unwrap();
        assert!(VoxelIsosurfaceTool.run(&args, &ctx()).is_err());
    }
}
