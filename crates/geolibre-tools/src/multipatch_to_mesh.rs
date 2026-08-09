//! GeoLibre tool: export 3D solids to OBJ or glTF 2.0 / GLB.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Multipatch To Collada* (Conversion).
//!
//! ## The 3D suite had no way out
//!
//! Rounds 16-18 built a real 3D solid-modeling suite — `union_3d`,
//! `difference_3d`, `intersect_3d`, `intersect_3d_line_with_surface`,
//! `enclose_multipatch`, `is_closed_3d`, `inside_3d`, `minimum_bounding_volume`,
//! `multipatch_footprint`, `voxel_isosurface`, `sun_shadow_volume`,
//! `extrude_between`, `fence_diagram`, `buffer_3d`. **None of it could be
//! exported to a 3D viewer**: no `gltf`, `obj`, `collada` or `mesh` tool id
//! existed anywhere in the ~1,100-tool registry, so every result was trapped in
//! a vector format only this catalog can read.
//!
//! That contradicts the repo's stated identity. GeoLibre ships `write_pmtiles`,
//! `vector_to_pmtiles`, `raster_to_h3`, `render_raster_png` and
//! `render_vector_png` for 2D web output; glTF is the equivalent for 3D — the
//! format three.js, deck.gl, Cesium and every browser viewer load natively.
//!
//! ## No new dependency
//!
//! glTF 2.0 is JSON plus a binary buffer, writable with `serde_json` and raw
//! little-endian byte packing. GLB is that JSON plus the buffer in a 12-byte
//! container. OBJ is plain text. Triangulation reuses
//! `inside_3d::collect_triangles`, which already understands this codebase's
//! multipatch layout (`Geometry::MultiPolygon` parts are `(Ring, Vec<Ring>)`
//! tuples; there is no `wbvector::Polygon` struct).
//!
//! ## The precision trap
//!
//! glTF positions are **`f32`**, which carries about 7 significant digits. A
//! UTM easting of 500000.0 therefore quantizes to roughly 0.03 m, and an ECEF
//! coordinate far worse — the mesh arrives visibly jittered. Translating to a
//! local origin is a correctness requirement here, not a convenience, so it is
//! **on by default** and the applied offset is reported so the model can be
//! placed back in the world.
//!
//! glTF carries no CRS, so the source EPSG is reported too; it cannot be
//! embedded in the file.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::FieldValue;

use crate::args_common::{choice_or, req_str};
use crate::common::{write_bytes, write_text_output};
use crate::inside_3d::collect_triangles;
use crate::vector_common::{load_input_layer, parse_optional_str};

const FORMATS: [&str; 3] = ["glb", "gltf", "obj"];
const ORIGINS: [&str; 3] = ["centroid", "min_corner", "none"];

/// glTF component type for `f32`.
const GL_FLOAT: u64 = 5126;
/// glTF component type for `u32`.
const GL_UNSIGNED_INT: u64 = 5125;
/// glTF primitive mode for triangles.
const GL_TRIANGLES: u64 = 4;
/// glTF bufferView target for vertex attributes.
const ARRAY_BUFFER: u64 = 34962;
/// glTF bufferView target for indices.
const ELEMENT_ARRAY_BUFFER: u64 = 34963;

pub struct MultipatchToMeshTool;

impl Tool for MultipatchToMeshTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "multipatch_to_mesh",
            display_name: "Multipatch To Mesh",
            summary: "Exports multipatch and 3D solid geometry to OBJ or glTF 2.0 / GLB for any 3D or web viewer (ArcGIS Multipatch To Collada). The union_3d / difference_3d / enclose_multipatch / voxel_isosurface suite had no export path at all — no mesh format existed anywhere in the registry — so 3D results could only be read back by this catalog.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Multipatch / 3D polygon layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output .glb, .gltf or .obj file path. The extension selects the format unless 'format' is given.",
                    required: true,
                },
                ToolParamSpec {
                    name: "format",
                    description: "'glb' (default for .glb), 'gltf' (JSON with an embedded base64 buffer), or 'obj'. Defaults from the output extension.",
                    required: false,
                },
                ToolParamSpec {
                    name: "origin",
                    description: "Local-origin translation: 'centroid' (default), 'min_corner', or 'none'. Required for precision — glTF positions are 32-bit floats.",
                    required: false,
                },
                ToolParamSpec {
                    name: "y_up",
                    description: "Swap to glTF's Y-up convention from geospatial Z-up (default true for glTF/GLB, false for OBJ).",
                    required: false,
                },
                ToolParamSpec {
                    name: "name_field",
                    description: "Attribute naming each feature's mesh in the output.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        let output = req_str(args, "output")?;
        let explicit = args.get("format").and_then(Value::as_str);
        if explicit.is_none() && format_from_path(output).is_none() {
            return Err(ToolError::Validation(
                "cannot infer the mesh format from 'output'; use a .glb, .gltf or .obj extension, \
                 or set 'format'"
                    .to_string(),
            ));
        }
        if explicit.is_some() {
            choice_or(args, "format", &FORMATS, "glb")?;
        }
        choice_or(args, "origin", &ORIGINS, "centroid")?;
        if let Some(v) = args.get("y_up") {
            if !v.is_boolean() && !v.is_null() {
                return Err(ToolError::Validation(
                    "'y_up' must be a boolean".to_string(),
                ));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?.to_string();
        let output = req_str(args, "output")?.to_string();
        let format = match args.get("format").and_then(Value::as_str) {
            Some(_) => choice_or(args, "format", &FORMATS, "glb")?.to_string(),
            None => format_from_path(&output)
                .ok_or_else(|| {
                    ToolError::Validation("cannot infer the mesh format from 'output'".into())
                })?
                .to_string(),
        };
        let origin_mode = choice_or(args, "origin", &ORIGINS, "centroid")?;
        // OBJ is plain text with f64 precision and no fixed axis convention, so
        // the Y-up swap only earns its keep for glTF.
        let y_up = args
            .get("y_up")
            .and_then(Value::as_bool)
            .unwrap_or(format != "obj");
        let name_field = parse_optional_str(args, "name_field")?.map(str::to_string);

        let layer = load_input_layer(&input)?;
        let name_idx = match &name_field {
            Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                ToolError::Validation(format!("name_field '{f}' not found in the input layer"))
            })?),
            None => None,
        };

        // One mesh per feature, triangulated via the existing 3D machinery.
        let mut meshes: Vec<Mesh> = Vec::new();
        let mut degenerate = 0_u64;
        let mut skipped = 0_u64;
        for (fid, feature) in layer.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                skipped += 1;
                continue;
            };
            let tris = collect_triangles(geom);
            if tris.is_empty() {
                skipped += 1;
                continue;
            }
            let mut kept = Vec::with_capacity(tris.len());
            for t in tris {
                // union_3d / difference_3d output can contain slivers; a
                // zero-area triangle has no normal and renders as a black seam.
                if triangle_area(&t) > 1e-12 {
                    kept.push(t);
                } else {
                    degenerate += 1;
                }
            }
            if kept.is_empty() {
                skipped += 1;
                continue;
            }
            let name = name_idx
                .and_then(|i| feature.attributes.get(i))
                .map(display)
                .filter(|s| !s.is_empty())
                .unwrap_or_else(|| format!("feature_{fid}"));
            meshes.push(Mesh { name, tris: kept });
        }

        if meshes.is_empty() {
            return Err(ToolError::Execution(
                "no feature yielded a triangle mesh; multipatch_to_mesh expects 3D polygon or \
                 multipatch geometry (see enclose_multipatch, union_3d, voxel_isosurface)"
                    .to_string(),
            ));
        }

        // Local origin. Computed over every vertex before any conversion, so
        // both formats place the model identically.
        let offset = compute_origin(&meshes, origin_mode);
        let tri_count: usize = meshes.iter().map(|m| m.tris.len()).sum();
        ctx.progress.info(&format!(
            "{} mesh(es), {tri_count} triangle(s), format {format}",
            meshes.len()
        ));

        let bytes_written = match format.as_str() {
            "obj" => {
                let (obj, mtl) = write_obj(&meshes, offset, y_up);
                write_text_output(&obj, &output)?;
                let mtl_path = replace_ext(&output, "mtl");
                write_text_output(&mtl, &mtl_path)?;
                obj.len()
            }
            "gltf" => {
                let json = build_gltf(&meshes, offset, y_up, true)?;
                let text = serde_json::to_string_pretty(&json).map_err(|e| {
                    ToolError::Execution(format!("failed serializing glTF JSON: {e}"))
                })?;
                write_text_output(&text, &output)?;
                text.len()
            }
            _ => {
                let glb = build_glb(&meshes, offset, y_up)?;
                write_bytes(&output, &glb)?;
                glb.len()
            }
        };

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(output));
        outputs.insert("format".to_string(), json!(format));
        outputs.insert("mesh_count".to_string(), json!(meshes.len()));
        outputs.insert("triangle_count".to_string(), json!(tri_count));
        outputs.insert(
            "vertex_count".to_string(),
            json!(meshes.iter().map(|m| m.tris.len() * 3).sum::<usize>()),
        );
        outputs.insert("origin_offset".to_string(), json!(offset));
        outputs.insert("y_up".to_string(), json!(y_up));
        outputs.insert("degenerate_triangles".to_string(), json!(degenerate));
        outputs.insert("skipped_features".to_string(), json!(skipped));
        outputs.insert("bytes_written".to_string(), json!(bytes_written));
        // glTF carries no CRS, so it can only be reported, never embedded.
        outputs.insert("source_epsg".to_string(), json!(layer.crs_epsg()));
        Ok(ToolRunResult { outputs })
    }
}

struct Mesh {
    name: String,
    tris: Vec<[[f64; 3]; 3]>,
}

fn format_from_path(p: &str) -> Option<&'static str> {
    let lower = p.to_ascii_lowercase();
    if lower.ends_with(".glb") {
        Some("glb")
    } else if lower.ends_with(".gltf") {
        Some("gltf")
    } else if lower.ends_with(".obj") {
        Some("obj")
    } else {
        None
    }
}

fn replace_ext(path: &str, ext: &str) -> String {
    match path.rfind('.') {
        Some(i) => format!("{}.{ext}", &path[..i]),
        None => format!("{path}.{ext}"),
    }
}

fn triangle_area(t: &[[f64; 3]; 3]) -> f64 {
    let u = sub(t[1], t[0]);
    let v = sub(t[2], t[0]);
    let c = cross(u, v);
    0.5 * (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Unit face normal, or `[0, 0, 1]` for a degenerate triangle (which the caller
/// has already filtered, so this is only a safety net).
fn face_normal(t: &[[f64; 3]; 3]) -> [f64; 3] {
    let n = cross(sub(t[1], t[0]), sub(t[2], t[0]));
    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if len <= 0.0 {
        [0.0, 0.0, 1.0]
    } else {
        [n[0] / len, n[1] / len, n[2] / len]
    }
}

fn compute_origin(meshes: &[Mesh], mode: &str) -> [f64; 3] {
    if mode == "none" {
        return [0.0; 3];
    }
    let mut min = [f64::MAX; 3];
    let mut max = [f64::MIN; 3];
    let mut sum = [0.0_f64; 3];
    let mut n = 0_usize;
    for m in meshes {
        for t in &m.tris {
            for v in t {
                for k in 0..3 {
                    min[k] = min[k].min(v[k]);
                    max[k] = max[k].max(v[k]);
                    sum[k] += v[k];
                }
                n += 1;
            }
        }
    }
    if n == 0 {
        return [0.0; 3];
    }
    match mode {
        "min_corner" => min,
        // Bounding-box centre rather than the vertex mean: a mesh with many
        // vertices clustered at one end would otherwise be offset lopsidedly
        // and lose the precision the translation exists to protect.
        _ => [
            0.5 * (min[0] + max[0]),
            0.5 * (min[1] + max[1]),
            0.5 * (min[2] + max[2]),
        ],
    }
}

/// Applies the origin translation and, when requested, the geospatial Z-up to
/// glTF Y-up conversion.
///
/// The swap is `(x, y, z) -> (x, z, -y)`, a right-handed rotation. Simply
/// swapping Y and Z would mirror the model, which is the classic way an
/// exported building ends up inside-out.
fn transform(v: [f64; 3], offset: [f64; 3], y_up: bool) -> [f64; 3] {
    let p = [v[0] - offset[0], v[1] - offset[1], v[2] - offset[2]];
    if y_up {
        [p[0], p[2], -p[1]]
    } else {
        p
    }
}

fn write_obj(meshes: &[Mesh], offset: [f64; 3], y_up: bool) -> (String, String) {
    let mut obj = String::new();
    obj.push_str("# generated by GeoLibre multipatch_to_mesh\n");
    obj.push_str(&format!("mtllib {}\n", "mesh.mtl"));
    let mut normals = String::new();
    let mut faces = String::new();
    let mut vi = 1_usize; // OBJ indices are 1-based
    let mut ni = 1_usize;
    let mut verts = String::new();

    for m in meshes {
        faces.push_str(&format!("o {}\nusemtl geolibre\n", sanitize(&m.name)));
        for t in &m.tris {
            let n = face_normal(t);
            let n = transform(
                [n[0] + offset[0], n[1] + offset[1], n[2] + offset[2]],
                offset,
                y_up,
            );
            normals.push_str(&format!("vn {} {} {}\n", n[0], n[1], n[2]));
            for v in t {
                let p = transform(*v, offset, y_up);
                verts.push_str(&format!("v {} {} {}\n", p[0], p[1], p[2]));
            }
            faces.push_str(&format!(
                "f {}//{} {}//{} {}//{}\n",
                vi,
                ni,
                vi + 1,
                ni,
                vi + 2,
                ni
            ));
            vi += 3;
            ni += 1;
        }
    }
    obj.push_str(&verts);
    obj.push_str(&normals);
    obj.push_str(&faces);

    let mtl = "# generated by GeoLibre multipatch_to_mesh\nnewmtl geolibre\nKa 0.2 0.2 0.2\n\
               Kd 0.8 0.8 0.8\nKs 0.0 0.0 0.0\nd 1.0\nillum 2\n"
        .to_string();
    (obj, mtl)
}

fn sanitize(s: &str) -> String {
    s.chars()
        .map(|c| if c.is_whitespace() { '_' } else { c })
        .collect()
}

/// Packs positions, normals and indices for one mesh into little-endian bytes.
///
/// Vertices are not shared between triangles: each face carries its own normal,
/// and welding vertices would average normals across faces and round off every
/// hard edge a building has.
struct Packed {
    positions: Vec<u8>,
    normals: Vec<u8>,
    indices: Vec<u8>,
    count: usize,
    min: [f32; 3],
    max: [f32; 3],
}

fn pack(mesh: &Mesh, offset: [f64; 3], y_up: bool) -> Packed {
    let mut positions = Vec::with_capacity(mesh.tris.len() * 36);
    let mut normals = Vec::with_capacity(mesh.tris.len() * 36);
    let mut indices = Vec::with_capacity(mesh.tris.len() * 12);
    let mut min = [f32::MAX; 3];
    let mut max = [f32::MIN; 3];
    let mut idx = 0_u32;

    for t in &mesh.tris {
        let n = face_normal(t);
        let n = if y_up { [n[0], n[2], -n[1]] } else { n };
        for v in t {
            let p = transform(*v, offset, y_up);
            for (k, lo) in min.iter_mut().enumerate() {
                let f = p[k] as f32;
                *lo = lo.min(f);
                max[k] = max[k].max(f);
                positions.extend_from_slice(&f.to_le_bytes());
            }
            for nk in n {
                normals.extend_from_slice(&(nk as f32).to_le_bytes());
            }
            indices.extend_from_slice(&idx.to_le_bytes());
            idx += 1;
        }
    }
    Packed {
        positions,
        normals,
        indices,
        count: idx as usize,
        min,
        max,
    }
}

/// Assembles the glTF JSON plus the binary buffer.
///
/// `embed` base64-encodes the buffer into a data URI (the `.gltf` form);
/// otherwise the buffer is left external for GLB to append.
fn build_gltf_parts(
    meshes: &[Mesh],
    offset: [f64; 3],
    y_up: bool,
) -> Result<(Value, Vec<u8>), ToolError> {
    let mut buffer: Vec<u8> = Vec::new();
    let mut views = Vec::new();
    let mut accessors = Vec::new();
    let mut gltf_meshes = Vec::new();
    let mut nodes = Vec::new();

    for (i, mesh) in meshes.iter().enumerate() {
        let p = pack(mesh, offset, y_up);

        // Every bufferView offset must be 4-byte aligned. Positions, normals
        // and u32 indices are all 4-byte scalars, so alignment holds as long as
        // each section starts on a multiple of 4 — asserted rather than assumed,
        // because a misaligned view parses but renders nothing.
        let mut push_view = |data: &[u8], target: u64, buffer: &mut Vec<u8>| -> usize {
            while !buffer.len().is_multiple_of(4) {
                buffer.push(0);
            }
            let offset = buffer.len();
            buffer.extend_from_slice(data);
            views.push(json!({
                "buffer": 0,
                "byteOffset": offset,
                "byteLength": data.len(),
                "target": target,
            }));
            views.len() - 1
        };

        let pos_view = push_view(&p.positions, ARRAY_BUFFER, &mut buffer);
        let nrm_view = push_view(&p.normals, ARRAY_BUFFER, &mut buffer);
        let idx_view = push_view(&p.indices, ELEMENT_ARRAY_BUFFER, &mut buffer);

        // POSITION requires min/max; viewers use it for culling and framing,
        // and omitting it is a validation error.
        accessors.push(json!({
            "bufferView": pos_view,
            "componentType": GL_FLOAT,
            "count": p.count,
            "type": "VEC3",
            "min": [p.min[0], p.min[1], p.min[2]],
            "max": [p.max[0], p.max[1], p.max[2]],
        }));
        let pos_acc = accessors.len() - 1;
        accessors.push(json!({
            "bufferView": nrm_view,
            "componentType": GL_FLOAT,
            "count": p.count,
            "type": "VEC3",
        }));
        let nrm_acc = accessors.len() - 1;
        accessors.push(json!({
            "bufferView": idx_view,
            "componentType": GL_UNSIGNED_INT,
            "count": p.count,
            "type": "SCALAR",
        }));
        let idx_acc = accessors.len() - 1;

        gltf_meshes.push(json!({
            "name": mesh.name,
            "primitives": [{
                "attributes": {"POSITION": pos_acc, "NORMAL": nrm_acc},
                "indices": idx_acc,
                "mode": GL_TRIANGLES,
                "material": 0,
            }],
        }));
        nodes.push(json!({"mesh": i, "name": mesh.name}));
    }

    // The buffer's own length must also be 4-byte aligned for GLB.
    while !buffer.len().is_multiple_of(4) {
        buffer.push(0);
    }

    let json_doc = json!({
        "asset": {"version": "2.0", "generator": "GeoLibre multipatch_to_mesh"},
        "scene": 0,
        "scenes": [{"nodes": (0..nodes.len()).collect::<Vec<_>>()}],
        "nodes": nodes,
        "meshes": gltf_meshes,
        "materials": [{
            "name": "geolibre",
            "pbrMetallicRoughness": {
                "baseColorFactor": [0.8, 0.8, 0.8, 1.0],
                "metallicFactor": 0.0,
                "roughnessFactor": 0.9,
            },
            "doubleSided": true,
        }],
        "accessors": accessors,
        "bufferViews": views,
        "buffers": [{"byteLength": buffer.len()}],
    });
    Ok((json_doc, buffer))
}

fn build_gltf(
    meshes: &[Mesh],
    offset: [f64; 3],
    y_up: bool,
    embed: bool,
) -> Result<Value, ToolError> {
    let (mut doc, buffer) = build_gltf_parts(meshes, offset, y_up)?;
    if embed {
        let uri = format!("data:application/octet-stream;base64,{}", base64(&buffer));
        doc["buffers"] = json!([{"byteLength": buffer.len(), "uri": uri}]);
    }
    Ok(doc)
}

/// GLB container: a 12-byte header followed by a JSON chunk and a BIN chunk,
/// each chunk itself 4-byte aligned and padded (JSON with spaces, BIN with
/// zeros, as the specification requires).
fn build_glb(meshes: &[Mesh], offset: [f64; 3], y_up: bool) -> Result<Vec<u8>, ToolError> {
    let (doc, buffer) = build_gltf_parts(meshes, offset, y_up)?;
    let mut json_bytes = serde_json::to_vec(&doc)
        .map_err(|e| ToolError::Execution(format!("failed serializing glTF JSON: {e}")))?;
    while !json_bytes.len().is_multiple_of(4) {
        json_bytes.push(b' ');
    }
    let mut bin = buffer;
    while !bin.len().is_multiple_of(4) {
        bin.push(0);
    }

    let total = 12 + 8 + json_bytes.len() + 8 + bin.len();
    let mut out = Vec::with_capacity(total);
    out.extend_from_slice(b"glTF");
    out.extend_from_slice(&2_u32.to_le_bytes());
    out.extend_from_slice(&(total as u32).to_le_bytes());

    out.extend_from_slice(&(json_bytes.len() as u32).to_le_bytes());
    out.extend_from_slice(b"JSON");
    out.extend_from_slice(&json_bytes);

    out.extend_from_slice(&(bin.len() as u32).to_le_bytes());
    out.extend_from_slice(b"BIN\0");
    out.extend_from_slice(&bin);
    Ok(out)
}

/// Minimal base64 encoder — avoids adding a dependency for one data URI.
fn base64(data: &[u8]) -> String {
    const T: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(data.len().div_ceil(3) * 4);
    for chunk in data.chunks(3) {
        let b = [
            chunk[0],
            *chunk.get(1).unwrap_or(&0),
            *chunk.get(2).unwrap_or(&0),
        ];
        let n = ((b[0] as u32) << 16) | ((b[1] as u32) << 8) | b[2] as u32;
        out.push(T[(n >> 18) as usize & 63] as char);
        out.push(T[(n >> 12) as usize & 63] as char);
        out.push(if chunk.len() > 1 {
            T[(n >> 6) as usize & 63] as char
        } else {
            '='
        });
        out.push(if chunk.len() > 2 {
            T[n as usize & 63] as char
        } else {
            '='
        });
    }
    out
}

fn display(v: &FieldValue) -> String {
    match v {
        FieldValue::Null => String::new(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) => format!("{f}"),
        FieldValue::Text(s) | FieldValue::Date(s) | FieldValue::DateTime(s) => s.clone(),
        FieldValue::Boolean(b) => b.to_string(),
        FieldValue::Blob(b) => format!("blob[{}]", b.len()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{Coord, FieldDef, FieldType, Geometry, GeometryType, Layer, Ring};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn tmp(tag: &str, ext: &str) -> String {
        static N: AtomicU64 = AtomicU64::new(0);
        let n = N.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir()
            .join(format!("m2m_{tag}_{}_{n}.{ext}", std::process::id()))
            .to_string_lossy()
            .to_string()
    }

    fn tri(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> (Ring, Vec<Ring>) {
        let ring = Ring(vec![
            Coord {
                x: a[0],
                y: a[1],
                z: Some(a[2]),
                m: None,
            },
            Coord {
                x: b[0],
                y: b[1],
                z: Some(b[2]),
                m: None,
            },
            Coord {
                x: c[0],
                y: c[1],
                z: Some(c[2]),
                m: None,
            },
            Coord {
                x: a[0],
                y: a[1],
                z: Some(a[2]),
                m: None,
            },
        ]);
        (ring, Vec::new())
    }

    /// A two-triangle square in the z = 0 plane, offset by `dx`.
    fn square_layer(dx: f64, epsg: Option<u32>) -> String {
        let mut l = Layer::new("solids").with_geom_type(GeometryType::MultiPolygon);
        if let Some(e) = epsg {
            l = l.with_crs_epsg(e);
        }
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.add_feature(
            Some(Geometry::MultiPolygon(vec![
                tri([dx, 0.0, 0.0], [dx + 1.0, 0.0, 0.0], [dx + 1.0, 1.0, 0.0]),
                tri([dx, 0.0, 0.0], [dx + 1.0, 1.0, 0.0], [dx, 1.0, 0.0]),
            ])),
            &[("name", "slab".into())],
        )
        .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: Value) -> (Vec<u8>, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = MultipatchToMeshTool.run(&args, &ctx()).unwrap();
        let p = res.outputs["output"].as_str().unwrap();
        let bytes = std::fs::read(p).unwrap();
        let _ = std::fs::remove_file(p);
        let _ = std::fs::remove_file(replace_ext(p, "mtl"));
        (bytes, res)
    }

    #[test]
    fn writes_a_glb_with_the_right_magic_and_length() {
        let (bytes, res) = run(json!({
            "input": square_layer(0.0, Some(3857)), "output": tmp("glb", "glb"),
        }));
        assert_eq!(&bytes[0..4], b"glTF");
        assert_eq!(u32::from_le_bytes(bytes[4..8].try_into().unwrap()), 2);
        // The header's declared length must match the file exactly, or every
        // loader rejects it.
        assert_eq!(
            u32::from_le_bytes(bytes[8..12].try_into().unwrap()) as usize,
            bytes.len()
        );
        assert_eq!(res.outputs["triangle_count"], json!(2));
        assert_eq!(res.outputs["format"], json!("glb"));
    }

    #[test]
    fn the_glb_chunks_are_four_byte_aligned_and_correctly_tagged() {
        // A misaligned bufferView parses but renders nothing, so this is the
        // silent-failure mode worth pinning down.
        let (bytes, _) = run(json!({
            "input": square_layer(0.0, Some(3857)), "output": tmp("align", "glb"),
        }));
        let json_len = u32::from_le_bytes(bytes[12..16].try_into().unwrap()) as usize;
        assert_eq!(&bytes[16..20], b"JSON");
        assert_eq!(json_len % 4, 0, "JSON chunk must be 4-byte aligned");
        let bin_start = 20 + json_len;
        let bin_len =
            u32::from_le_bytes(bytes[bin_start..bin_start + 4].try_into().unwrap()) as usize;
        assert_eq!(&bytes[bin_start + 4..bin_start + 8], b"BIN\0");
        assert_eq!(bin_len % 4, 0, "BIN chunk must be 4-byte aligned");
        assert_eq!(bin_start + 8 + bin_len, bytes.len());
    }

    #[test]
    fn the_glb_json_parses_and_declares_a_complete_mesh() {
        let (bytes, _) = run(json!({
            "input": square_layer(0.0, Some(3857)), "output": tmp("json", "glb"),
        }));
        let json_len = u32::from_le_bytes(bytes[12..16].try_into().unwrap()) as usize;
        let doc: Value = serde_json::from_slice(&bytes[20..20 + json_len]).unwrap();
        assert_eq!(doc["asset"]["version"], json!("2.0"));
        assert_eq!(doc["meshes"].as_array().unwrap().len(), 1);
        let prim = &doc["meshes"][0]["primitives"][0];
        assert_eq!(prim["mode"], json!(4));
        assert!(prim["attributes"]["POSITION"].is_number());
        // Without normals most viewers render the mesh flat black.
        assert!(prim["attributes"]["NORMAL"].is_number());
        assert!(prim["indices"].is_number());
        // POSITION accessors require min/max; omitting them is a validation
        // error and breaks camera framing.
        let pos_acc = prim["attributes"]["POSITION"].as_u64().unwrap() as usize;
        let acc = &doc["accessors"][pos_acc];
        assert_eq!(acc["count"], json!(6));
        assert_eq!(acc["min"].as_array().unwrap().len(), 3);
        assert_eq!(acc["max"].as_array().unwrap().len(), 3);
    }

    #[test]
    fn buffer_views_start_on_four_byte_boundaries() {
        let (bytes, _) = run(json!({
            "input": square_layer(0.0, Some(3857)), "output": tmp("views", "glb"),
        }));
        let json_len = u32::from_le_bytes(bytes[12..16].try_into().unwrap()) as usize;
        let doc: Value = serde_json::from_slice(&bytes[20..20 + json_len]).unwrap();
        for v in doc["bufferViews"].as_array().unwrap() {
            assert_eq!(
                v["byteOffset"].as_u64().unwrap() % 4,
                0,
                "bufferView not aligned: {v}"
            );
        }
    }

    #[test]
    fn the_declared_buffer_length_matches_the_bin_chunk() {
        let (bytes, _) = run(json!({
            "input": square_layer(0.0, Some(3857)), "output": tmp("buflen", "glb"),
        }));
        let json_len = u32::from_le_bytes(bytes[12..16].try_into().unwrap()) as usize;
        let doc: Value = serde_json::from_slice(&bytes[20..20 + json_len]).unwrap();
        let declared = doc["buffers"][0]["byteLength"].as_u64().unwrap() as usize;
        let bin_start = 20 + json_len;
        let bin_len =
            u32::from_le_bytes(bytes[bin_start..bin_start + 4].try_into().unwrap()) as usize;
        assert_eq!(declared, bin_len);
    }

    #[test]
    fn a_far_from_origin_mesh_keeps_its_shape_in_f32() {
        // The precision trap: raw UTM eastings quantize to centimetres in f32
        // and the mesh arrives visibly jittered. With the default centroid
        // origin the 1-unit square must survive intact.
        let (bytes, res) = run(json!({
            "input": square_layer(500_000.0, Some(32610)), "output": tmp("prec", "glb"),
        }));
        let offset = res.outputs["origin_offset"].as_array().unwrap();
        assert!(
            offset[0].as_f64().unwrap() > 400_000.0,
            "origin not applied"
        );

        let json_len = u32::from_le_bytes(bytes[12..16].try_into().unwrap()) as usize;
        let doc: Value = serde_json::from_slice(&bytes[20..20 + json_len]).unwrap();
        let pos_acc = doc["meshes"][0]["primitives"][0]["attributes"]["POSITION"]
            .as_u64()
            .unwrap() as usize;
        let acc = &doc["accessors"][pos_acc];
        let span = acc["max"][0].as_f64().unwrap() - acc["min"][0].as_f64().unwrap();
        assert!(
            (span - 1.0).abs() < 1e-5,
            "the unit square lost its width in f32: span {span}"
        );
    }

    #[test]
    fn origin_none_leaves_world_coordinates_in_place() {
        let (_, res) = run(json!({
            "input": square_layer(500_000.0, Some(32610)),
            "output": tmp("noorigin", "glb"), "origin": "none",
        }));
        assert_eq!(res.outputs["origin_offset"], json!([0.0, 0.0, 0.0]));
    }

    #[test]
    fn min_corner_origin_places_the_mesh_at_the_positive_octant() {
        let (bytes, _) = run(json!({
            "input": square_layer(10.0, Some(3857)),
            "output": tmp("mincorner", "glb"), "origin": "min_corner", "y_up": false,
        }));
        let json_len = u32::from_le_bytes(bytes[12..16].try_into().unwrap()) as usize;
        let doc: Value = serde_json::from_slice(&bytes[20..20 + json_len]).unwrap();
        let pos_acc = doc["meshes"][0]["primitives"][0]["attributes"]["POSITION"]
            .as_u64()
            .unwrap() as usize;
        let min = &doc["accessors"][pos_acc]["min"];
        for k in 0..3 {
            assert!(
                min[k].as_f64().unwrap().abs() < 1e-6,
                "min not at origin: {min}"
            );
        }
    }

    #[test]
    fn the_y_up_swap_is_a_rotation_not_a_mirror() {
        // Swapping Y and Z outright would flip handedness and turn an exported
        // building inside-out. (x, y, z) -> (x, z, -y) preserves it.
        let v = transform([1.0, 2.0, 3.0], [0.0; 3], true);
        assert_eq!(v, [1.0, 3.0, -2.0]);
        // A right-handed basis must stay right-handed.
        let ex = transform([1.0, 0.0, 0.0], [0.0; 3], true);
        let ey = transform([0.0, 1.0, 0.0], [0.0; 3], true);
        let ez = transform([0.0, 0.0, 1.0], [0.0; 3], true);
        let c = cross(ex, ey);
        let dot = c[0] * ez[0] + c[1] * ez[1] + c[2] * ez[2];
        assert!(dot > 0.0, "handedness flipped: {dot}");
    }

    #[test]
    fn obj_output_has_vertices_normals_and_faces() {
        let (bytes, res) = run(json!({
            "input": square_layer(0.0, Some(3857)), "output": tmp("obj", "obj"),
        }));
        let text = String::from_utf8(bytes).unwrap();
        assert_eq!(res.outputs["format"], json!("obj"));
        assert_eq!(text.lines().filter(|l| l.starts_with("v ")).count(), 6);
        assert_eq!(text.lines().filter(|l| l.starts_with("vn ")).count(), 2);
        assert_eq!(text.lines().filter(|l| l.starts_with("f ")).count(), 2);
        // OBJ indices are 1-based; a 0 index makes the file unloadable.
        assert!(text.contains("f 1//1 2//1 3//1"), "got: {text}");
    }

    #[test]
    fn obj_writes_a_material_sidecar() {
        let out = tmp("mtl", "obj");
        let args: ToolArgs = serde_json::from_value(json!({
            "input": square_layer(0.0, Some(3857)), "output": out.clone(),
        }))
        .unwrap();
        MultipatchToMeshTool.run(&args, &ctx()).unwrap();
        let mtl = replace_ext(&out, "mtl");
        let text = std::fs::read_to_string(&mtl).unwrap();
        let _ = std::fs::remove_file(&out);
        let _ = std::fs::remove_file(&mtl);
        assert!(text.contains("newmtl geolibre"));
    }

    #[test]
    fn gltf_output_embeds_the_buffer_as_a_data_uri() {
        let (bytes, res) = run(json!({
            "input": square_layer(0.0, Some(3857)), "output": tmp("gltf", "gltf"),
        }));
        assert_eq!(res.outputs["format"], json!("gltf"));
        let doc: Value = serde_json::from_slice(&bytes).unwrap();
        let uri = doc["buffers"][0]["uri"].as_str().unwrap();
        assert!(uri.starts_with("data:application/octet-stream;base64,"));
        // A .gltf with an unresolvable external buffer loads to nothing.
        assert!(uri.len() > 40);
    }

    #[test]
    fn the_base64_encoder_matches_known_vectors() {
        assert_eq!(base64(b""), "");
        assert_eq!(base64(b"f"), "Zg==");
        assert_eq!(base64(b"fo"), "Zm8=");
        assert_eq!(base64(b"foo"), "Zm9v");
        assert_eq!(base64(b"foob"), "Zm9vYg==");
        assert_eq!(base64(b"foobar"), "Zm9vYmFy");
    }

    #[test]
    fn degenerate_triangles_are_culled_and_counted() {
        let mut l = Layer::new("s").with_geom_type(GeometryType::MultiPolygon);
        l.add_feature(
            Some(Geometry::MultiPolygon(vec![
                tri([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]),
                // Collinear: zero area, no normal, renders as a black seam.
                tri([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]),
            ])),
            &[],
        )
        .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let path = wbvector::memory_store::make_vector_memory_path(&id);
        let (_, res) = run(json!({"input": path, "output": tmp("degen", "glb")}));
        assert_eq!(res.outputs["degenerate_triangles"], json!(1));
        assert_eq!(res.outputs["triangle_count"], json!(1));
    }

    #[test]
    fn each_feature_becomes_its_own_named_mesh() {
        let mut l = Layer::new("s").with_geom_type(GeometryType::MultiPolygon);
        l.add_field(FieldDef::new("name", FieldType::Text));
        for (i, n) in ["tower", "annex"].iter().enumerate() {
            let dx = i as f64 * 5.0;
            l.add_feature(
                Some(Geometry::MultiPolygon(vec![tri(
                    [dx, 0.0, 0.0],
                    [dx + 1.0, 0.0, 0.0],
                    [dx + 1.0, 1.0, 0.0],
                )])),
                &[("name", (*n).into())],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        let path = wbvector::memory_store::make_vector_memory_path(&id);
        let (bytes, res) = run(json!({
            "input": path, "output": tmp("named", "glb"), "name_field": "name",
        }));
        assert_eq!(res.outputs["mesh_count"], json!(2));
        let json_len = u32::from_le_bytes(bytes[12..16].try_into().unwrap()) as usize;
        let doc: Value = serde_json::from_slice(&bytes[20..20 + json_len]).unwrap();
        assert_eq!(doc["meshes"][0]["name"], json!("tower"));
        assert_eq!(doc["meshes"][1]["name"], json!("annex"));
        assert_eq!(doc["scenes"][0]["nodes"], json!([0, 1]));
    }

    #[test]
    fn the_source_crs_is_reported_since_gltf_cannot_carry_it() {
        let (_, res) = run(json!({
            "input": square_layer(0.0, Some(32610)), "output": tmp("crs", "glb"),
        }));
        assert_eq!(res.outputs["source_epsg"], json!(32610));
    }

    #[test]
    fn a_layer_with_no_3d_geometry_is_refused() {
        let mut l = Layer::new("pts").with_geom_type(GeometryType::Point);
        l.add_feature(Some(Geometry::point(0.0, 0.0)), &[]).unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let path = wbvector::memory_store::make_vector_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": path, "output": tmp("nope", "glb"),
        }))
        .unwrap();
        let err = MultipatchToMeshTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("triangle mesh"), "{err}");
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            MultipatchToMeshTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"input": "a.shp"})));
        // An unknown extension with no explicit format is ambiguous.
        assert!(bad(json!({"input": "a.shp", "output": "mesh.dae"})));
        assert!(bad(
            json!({"input": "a.shp", "output": "m.glb", "format": "collada"})
        ));
        assert!(bad(
            json!({"input": "a.shp", "output": "m.glb", "origin": "corner"})
        ));
        assert!(bad(
            json!({"input": "a.shp", "output": "m.glb", "y_up": "yes"})
        ));
    }
}
