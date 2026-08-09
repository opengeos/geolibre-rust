//! GeoLibre tool: split 3D lines where they cross a surface or a multipatch.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Intersect 3D Line With Surface* and
//! *Intersect 3D Line With Multi Patch* (3D Analyst).
//!
//! ## Why the catalog needs it
//!
//! Almost every 3D line means something different above and below the ground.
//! A sightline is blocked where it enters terrain; a proposed pipeline is
//! trenched where it runs under the surface and on piers where it does not; a
//! drill trace changes formation at each contact. The primitive all of those
//! need is the same: split the line at every crossing and label each piece.
//!
//! Round 17 added the 3D overlay suite, but none of it answers this:
//!
//! * `intersect_3d` computes the shared **volume of two closed solids** — it
//!   takes no lines at all;
//! * `inside_3d` is a containment *predicate*: it reports which parts are
//!   inside, but does not cut the geometry or find the crossing points;
//! * `line_of_sight` returns visibility between an observer and a target, not
//!   the geometry of a general 3D line against a surface;
//! * `interpolate_shape` drapes a line **onto** a surface, discarding the
//!   line's own z — the very quantity the comparison needs.
//!
//! ## Method
//!
//! Each line is densified in plan view to about one cell (or `spacing`), and at
//! every vertex the surface height is sampled bilinearly and compared with the
//! line's own interpolated z. A sign change in that difference is a crossing;
//! its position is found by linear interpolation between the bracketing
//! vertices, which is exact for the piecewise-linear geometry involved. The
//! line is cut there and each piece is tagged `above` or `below`.
//!
//! A `multipatch` may be supplied instead of a raster, in which case each
//! densified segment is intersected against the mesh triangles with the shared
//! `mesh3d::segment_triangle` routine — the *Intersect 3D Line With Multi
//! Patch* case.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::args_common::{band_index, bool_or, opt_positive_f64, req_str};
use crate::common::{load_input_raster, parse_optional_output};
use crate::inside_3d::{collect_triangles, Tri};
use crate::mesh3d::segment_triangle;
use crate::surface_solid::sample_bilinear;
use crate::vector_common::{load_input_layer, write_or_store_layer};

/// Upper bound on the samples one segment may be densified into. Only a
/// pathologically small `spacing` reaches it; the point is that it is reached
/// rather than exhausting memory.
const MAX_DENSIFY_STEPS: usize = 1_000_000;

pub struct Intersect3dLineWithSurfaceTool;

impl Tool for Intersect3dLineWithSurfaceTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "intersect_3d_line_with_surface",
            display_name: "Intersect 3D Line With Surface",
            summary: "Splits 3D lines at every point where they cross an elevation surface or a multipatch, tagging each piece above or below and emitting the crossing points (ArcGIS Intersect 3D Line With Surface / Intersect 3D Line With Multi Patch). Round 17's intersect_3d takes closed solids rather than lines, inside_3d is a containment predicate that neither cuts the geometry nor locates the crossings, and interpolate_shape drapes a line onto a surface, discarding the line z the comparison depends on.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "3D line layer, with z-bearing vertices.",
                    required: true,
                },
                ToolParamSpec {
                    name: "surface",
                    description: "Elevation raster to cut against. Required unless 'multipatch' is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "multipatch",
                    description: "Multipatch layer (triangle multipolygons) to cut against, instead of a raster surface.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output line layer split at every crossing. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_points",
                    description: "Output point layer of the crossings themselves. Always produced; stored in memory when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "spacing",
                    description: "Plan-view densification step in map units (default: one surface cell). Smaller resolves crossings more finely at more cost.",
                    required: false,
                },
                ToolParamSpec {
                    name: "keep_above",
                    description: "Keep the pieces above the surface (default true).",
                    required: false,
                },
                ToolParamSpec {
                    name: "keep_below",
                    description: "Keep the pieces below the surface (default true).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band of the surface raster (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        let prm = parse_params(args)?;
        let has_surface = args
            .get("surface")
            .and_then(Value::as_str)
            .is_some_and(|s| !s.trim().is_empty());
        let has_mesh = args
            .get("multipatch")
            .and_then(Value::as_str)
            .is_some_and(|s| !s.trim().is_empty());
        if !has_surface && !has_mesh {
            return Err(ToolError::Validation(
                "one of 'surface' or 'multipatch' is required".to_string(),
            ));
        }
        if has_surface && has_mesh {
            return Err(ToolError::Validation(
                "'surface' and 'multipatch' are alternatives; supply exactly one".to_string(),
            ));
        }
        if !prm.keep_above && !prm.keep_below {
            return Err(ToolError::Validation(
                "'keep_above' and 'keep_below' cannot both be false; nothing would be written"
                    .to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        self.validate(args)?;
        let input_path = req_str(args, "input")?.to_string();
        let prm = parse_params(args)?;
        let band = band_index(args, "band")?;
        let output = parse_optional_output(args, "output")?;
        let out_points = parse_optional_output(args, "output_points")?;

        let lines = load_input_layer(&input_path)?;

        // Either a raster surface or a triangle mesh.
        let surface = match args.get("surface").and_then(Value::as_str) {
            Some(p) if !p.trim().is_empty() => Some(load_input_raster(p.trim())?),
            _ => None,
        };
        let mesh = match args.get("multipatch").and_then(Value::as_str) {
            Some(p) if !p.trim().is_empty() => Some(layer_triangles(&load_input_layer(p.trim())?)),
            _ => None,
        };

        let spacing = match prm.spacing {
            Some(s) => s,
            None => match &surface {
                Some(r) => r.cell_size_x.min(r.cell_size_y),
                // With no raster to take a scale from, fall back to the mesh's
                // own extent rather than an arbitrary constant.
                None => mesh_spacing(mesh.as_deref()),
            },
        };
        if spacing <= 0.0 || !spacing.is_finite() {
            return Err(ToolError::Execution(
                "could not determine a densification spacing; supply 'spacing'".to_string(),
            ));
        }

        let mesh_z = mesh_z_range(mesh.as_deref().unwrap_or(&[]));

        ctx.progress.info(&format!(
            "{} line(s) against a {}, spacing {spacing}",
            lines.len(),
            if surface.is_some() {
                "raster surface"
            } else {
                "multipatch"
            }
        ));

        let mut out_lines = Layer::new("split_lines");
        out_lines.geom_type = Some(GeometryType::LineString);
        if let Some(e) = lines.crs_epsg() {
            out_lines = out_lines.with_crs_epsg(e);
        }
        out_lines.add_field(FieldDef::new("id", FieldType::Integer));
        out_lines.add_field(FieldDef::new("source_fid", FieldType::Integer));
        out_lines.add_field(FieldDef::new("part", FieldType::Integer));
        out_lines.add_field(FieldDef::new("position", FieldType::Text));
        out_lines.add_field(FieldDef::new("length_3d", FieldType::Float));
        out_lines.add_field(FieldDef::new("min_clearance", FieldType::Float));
        out_lines.add_field(FieldDef::new("max_clearance", FieldType::Float));

        let mut out_pts = Layer::new("crossings");
        out_pts.geom_type = Some(GeometryType::Point);
        if let Some(e) = lines.crs_epsg() {
            out_pts = out_pts.with_crs_epsg(e);
        }
        out_pts.add_field(FieldDef::new("id", FieldType::Integer));
        out_pts.add_field(FieldDef::new("source_fid", FieldType::Integer));
        out_pts.add_field(FieldDef::new("z", FieldType::Float));
        out_pts.add_field(FieldDef::new("direction", FieldType::Text));

        let mut line_fid = 0u64;
        let mut point_fid = 0u64;
        let mut crossings_total = 0usize;

        for (n, feat) in lines.iter().enumerate() {
            let Some(coords) = line_coords(feat.geometry.as_ref()) else {
                continue;
            };
            if coords.len() < 2 {
                continue;
            }

            // Densify in plan view and evaluate the line-minus-surface
            // clearance at each sample.
            let samples = densify_3d(&coords, spacing);
            // Split into runs at every sample the surface has no value for.
            // Dropping those samples and treating the rest as continuous makes
            // the sample before a gap adjacent to the one after it, so a sign
            // change across the gap reads as a crossing and its position is
            // interpolated between two points that may be far apart.
            let mut runs: Vec<Vec<(Coord, f64)>> = Vec::new();
            let mut run: Vec<(Coord, f64)> = Vec::new();
            for p in samples {
                match surface_height(&surface, band, mesh.as_deref(), p.x, p.y, mesh_z) {
                    Some(sz) => run.push((p.clone(), p.z.unwrap_or(0.0) - sz)),
                    None => {
                        if run.len() >= 2 {
                            runs.push(std::mem::take(&mut run));
                        } else {
                            run.clear();
                        }
                    }
                }
            }
            if run.len() >= 2 {
                runs.push(run);
            }
            if runs.is_empty() {
                continue;
            }
            for evaluated in runs {
                // Cut at every sign change. Each part records the inclusive range
                // of *real* samples it covers, because `coords` also carries the
                // interpolated cut vertices and those have no clearance value.
                let mut parts: Vec<(Vec<Coord>, bool, usize, usize)> = Vec::new();
                let mut current: Vec<Coord> = vec![evaluated[0].0.clone()];
                let mut above = evaluated[0].1 >= 0.0;
                let mut sample_start = 0usize;
                for (si, w) in evaluated.windows(2).enumerate() {
                    let (a, da) = (&w[0].0, w[0].1);
                    let (b, db) = (&w[1].0, w[1].1);
                    if (da >= 0.0) != (db >= 0.0) && (da - db).abs() > 0.0 {
                        // Linear interpolation is exact here: between two densified
                        // samples both the line and the sampled surface are linear,
                        // so their difference is too.
                        let t = da / (da - db);
                        let x = a.x + t * (b.x - a.x);
                        let y = a.y + t * (b.y - a.y);
                        let z = a.z.unwrap_or(0.0) + t * (b.z.unwrap_or(0.0) - a.z.unwrap_or(0.0));
                        let cut = Coord::xyz(x, y, z);

                        current.push(cut.clone());
                        // This part owns samples `sample_start ..= si`; the sample
                        // at `si + 1` is already on the other side of the surface.
                        parts.push((std::mem::take(&mut current), above, sample_start, si));
                        current.push(cut.clone());
                        above = db >= 0.0;
                        sample_start = si + 1;
                        crossings_total += 1;

                        let mut pf = Feature::with_geometry(
                            point_fid,
                            Geometry::Point(cut),
                            out_pts.schema.len(),
                        );
                        pf.set_by_index(0, FieldValue::Integer(point_fid as i64));
                        pf.set_by_index(1, FieldValue::Integer(feat.fid as i64));
                        pf.set_by_index(2, FieldValue::Float(z));
                        pf.set_by_index(
                            3,
                            FieldValue::Text(
                                if db >= 0.0 { "emerging" } else { "entering" }.to_string(),
                            ),
                        );
                        out_pts.push(pf);
                        point_fid += 1;
                    }
                    current.push(b.clone());
                }
                if current.len() >= 2 {
                    parts.push((current, above, sample_start, evaluated.len() - 1));
                }

                // Clearance statistics per part, from the real samples it spans.
                let mut part_index = 0i64;
                for (coords, is_above, lo, hi) in parts {
                    let stats: Vec<f64> =
                        evaluated[lo..=hi.max(lo)].iter().map(|(_, d)| *d).collect();

                    let keep = if is_above {
                        prm.keep_above
                    } else {
                        prm.keep_below
                    };
                    if !keep {
                        part_index += 1;
                        continue;
                    }
                    let length = length_3d(&coords);
                    let (mn, mx) = if stats.is_empty() {
                        (0.0, 0.0)
                    } else {
                        (
                            stats.iter().copied().fold(f64::INFINITY, f64::min),
                            stats.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                        )
                    };

                    let mut lf = Feature::with_geometry(
                        line_fid,
                        Geometry::LineString(coords),
                        out_lines.schema.len(),
                    );
                    lf.set_by_index(0, FieldValue::Integer(line_fid as i64));
                    lf.set_by_index(1, FieldValue::Integer(feat.fid as i64));
                    lf.set_by_index(2, FieldValue::Integer(part_index));
                    lf.set_by_index(
                        3,
                        FieldValue::Text(if is_above { "above" } else { "below" }.to_string()),
                    );
                    lf.set_by_index(4, FieldValue::Float(length));
                    lf.set_by_index(5, FieldValue::Float(mn));
                    lf.set_by_index(6, FieldValue::Float(mx));
                    out_lines.push(lf);
                    line_fid += 1;
                    part_index += 1;
                }
            }
            ctx.progress
                .progress((n as f64 + 1.0) / lines.len().max(1) as f64);
        }

        let line_count = out_lines.len();
        let point_count = out_pts.len();
        let out_path = write_or_store_layer(out_lines, output)?;
        let pts_path = write_or_store_layer(out_pts, out_points)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("output_points".to_string(), json!(pts_path));
        outputs.insert("part_count".to_string(), json!(line_count));
        outputs.insert("crossing_count".to_string(), json!(point_count));
        outputs.insert("crossings_found".to_string(), json!(crossings_total));
        outputs.insert("spacing".to_string(), json!(spacing));
        Ok(ToolRunResult { outputs })
    }
}

/// The coordinate list of a line feature, if it is one.
fn line_coords(geom: Option<&Geometry>) -> Option<Vec<Coord>> {
    match geom? {
        Geometry::LineString(c) => Some(c.clone()),
        // A multi-part line is flattened: the parts are treated as one
        // sequence, which is what splitting by z-crossing means for a route.
        Geometry::MultiLineString(parts) => {
            Some(parts.iter().flat_map(|p| p.iter().cloned()).collect())
        }
        _ => None,
    }
}

/// Densifies a 3D polyline in plan view, interpolating z along the way.
fn densify_3d(coords: &[Coord], spacing: f64) -> Vec<Coord> {
    let mut out: Vec<Coord> = Vec::new();
    for w in coords.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        let (za, zb) = (a.z.unwrap_or(0.0), b.z.unwrap_or(0.0));
        let dx = b.x - a.x;
        let dy = b.y - a.y;
        let plan = (dx * dx + dy * dy).sqrt();
        // Clamped as well as guarded: `as usize` saturates rather than
        // overflowing, so an accidentally tiny spacing would otherwise turn
        // this into an unbounded push loop.
        let steps = ((plan / spacing).ceil() as usize).clamp(1, MAX_DENSIFY_STEPS);
        for s in 0..steps {
            let t = s as f64 / steps as f64;
            out.push(Coord::xyz(a.x + t * dx, a.y + t * dy, za + t * (zb - za)));
        }
    }
    if let Some(last) = coords.last() {
        out.push(Coord::xyz(last.x, last.y, last.z.unwrap_or(0.0)));
    }
    out
}

/// 3D length of a coordinate run.
fn length_3d(coords: &[Coord]) -> f64 {
    coords
        .windows(2)
        .map(|w| {
            let dx = w[1].x - w[0].x;
            let dy = w[1].y - w[0].y;
            let dz = w[1].z.unwrap_or(0.0) - w[0].z.unwrap_or(0.0);
            (dx * dx + dy * dy + dz * dz).sqrt()
        })
        .sum()
}

/// Surface height under a plan position, from a raster or a mesh.
///
/// For a mesh the height is the highest triangle hit by a vertical ray, which
/// is the surface an overflying line would meet first.
fn surface_height(
    surface: &Option<Raster>,
    band: isize,
    mesh: Option<&[Tri]>,
    x: f64,
    y: f64,
    // Precomputed once per run: the range does not change, and recomputing it
    // here scanned every vertex of the mesh for every densified sample.
    z_range: (f64, f64),
) -> Option<f64> {
    if let Some(r) = surface {
        return sample_bilinear(r, band, x, y);
    }
    let tris = mesh?;
    // Shoot a tall vertical segment and keep the highest hit.
    let (lo, hi) = z_range;
    let pad = (hi - lo).abs().max(1.0);
    let a = [x, y, lo - pad];
    let b = [x, y, hi + pad];
    let mut best: Option<f64> = None;
    for t in tris {
        if let Some(param) = segment_triangle(a, b, t) {
            let z = a[2] + param * (b[2] - a[2]);
            best = Some(match best {
                None => z,
                Some(prev) => prev.max(z),
            });
        }
    }
    best
}

/// Triangles of every feature in a multipatch layer, via the shared
/// `inside_3d::collect_triangles` (the crate's convention is that each part of
/// a `MultiPolygon` is one triangle).
fn layer_triangles(layer: &Layer) -> Vec<Tri> {
    layer
        .iter()
        .filter_map(|f| f.geometry.as_ref())
        .flat_map(collect_triangles)
        .collect()
}

/// Lowest and highest z over a triangle list; `(0, 0)` when it is empty.
fn mesh_z_range(tris: &[Tri]) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for t in tris {
        for v in t {
            lo = lo.min(v[2]);
            hi = hi.max(v[2]);
        }
    }
    if lo.is_finite() {
        (lo, hi)
    } else {
        (0.0, 0.0)
    }
}

/// A default densification step for a mesh: a fiftieth of its plan diagonal.
fn mesh_spacing(mesh: Option<&[Tri]>) -> f64 {
    let Some(tris) = mesh else { return 0.0 };
    let (mut x0, mut y0) = (f64::INFINITY, f64::INFINITY);
    let (mut x1, mut y1) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for t in tris {
        for v in t {
            x0 = x0.min(v[0]);
            x1 = x1.max(v[0]);
            y0 = y0.min(v[1]);
            y1 = y1.max(v[1]);
        }
    }
    if !x0.is_finite() {
        return 0.0;
    }
    let diag = ((x1 - x0).powi(2) + (y1 - y0).powi(2)).sqrt();
    // A mesh with no plan extent has no scale to derive a step from. Returning
    // `f64::MIN_POSITIVE` instead would clear the caller's `> 0` guard, make
    // `plan / spacing` infinite in `densify_3d`, saturate the `as usize` cast
    // to `usize::MAX`, and push coordinates until memory ran out. Return 0 and
    // let the caller report that it needs an explicit `spacing`.
    diag / 50.0
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    spacing: Option<f64>,
    keep_above: bool,
    keep_below: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    Ok(Params {
        spacing: opt_positive_f64(args, "spacing")?,
        keep_above: bool_or(args, "keep_above", true)?,
        keep_below: bool_or(args, "keep_below", true)?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, RasterConfig};
    use wbvector::Ring;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// A flat surface at `height`, 10 m cells over a 10x10 grid.
    fn flat_surface(height: f64) -> String {
        let (rows, cols) = (10usize, 10usize);
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 10.0,
            cell_size_y: Some(10.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(32610),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, height).unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn line_layer(pts: &[(f64, f64, f64)]) -> String {
        let mut layer = Layer::new("lines");
        layer.geom_type = Some(GeometryType::LineString);
        layer = layer.with_crs_epsg(32610);
        layer.add_field(FieldDef::new("id", FieldType::Integer));
        let coords: Vec<Coord> = pts.iter().map(|&(x, y, z)| Coord::xyz(x, y, z)).collect();
        let mut f = Feature::with_geometry(0, Geometry::LineString(coords), layer.schema.len());
        f.set_by_index(0, FieldValue::Integer(0));
        layer.push(f);
        write_or_store_layer(layer, None).unwrap()
    }

    fn run(args: Value) -> (Layer, BTreeMap<String, Value>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = Intersect3dLineWithSurfaceTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (layer, out.outputs)
    }

    fn text(layer: &Layer, f: &Feature, name: &str) -> String {
        match &f.attributes[layer.schema.field_index(name).unwrap()] {
            FieldValue::Text(t) => t.clone(),
            other => panic!("{name} should be text, got {other:?}"),
        }
    }

    /// A line rising through a flat surface splits into exactly two parts at the
    /// analytically known crossing point.
    #[test]
    fn splits_at_a_known_crossing() {
        // Surface at z = 50. The line runs from (5, 50, 0) to (95, 50, 100),
        // so z = 50 at exactly the midpoint, x = 50.
        let (layer, outputs) = run(json!({
            "input": line_layer(&[(5.0, 50.0, 0.0), (95.0, 50.0, 100.0)]),
            "surface": flat_surface(50.0), "spacing": 1.0
        }));
        assert_eq!(outputs["crossing_count"].as_u64().unwrap(), 1);
        assert_eq!(layer.len(), 2, "one crossing gives two parts");

        let parts: Vec<&Feature> = layer.iter().collect();
        assert_eq!(text(&layer, parts[0], "position"), "below");
        assert_eq!(text(&layer, parts[1], "position"), "above");

        // The crossing point is at x = 50, z = 50.
        let args: ToolArgs = serde_json::from_value(json!({
            "input": line_layer(&[(5.0, 50.0, 0.0), (95.0, 50.0, 100.0)]),
            "surface": flat_surface(50.0), "spacing": 1.0
        }))
        .unwrap();
        let out = Intersect3dLineWithSurfaceTool.run(&args, &ctx()).unwrap();
        let pts = load_input_layer(out.outputs["output_points"].as_str().unwrap()).unwrap();
        let Some(Geometry::Point(p)) = pts.iter().next().unwrap().geometry.as_ref() else {
            panic!("expected a point")
        };
        assert!((p.x - 50.0).abs() < 0.5, "crossing x {} != 50", p.x);
        assert!(
            (p.z.unwrap_or(0.0) - 50.0).abs() < 0.5,
            "crossing z {:?} != 50",
            p.z
        );
    }

    /// A line that dips under and comes back out yields three parts and two
    /// crossings, correctly labelled.
    #[test]
    fn handles_multiple_crossings() {
        let (layer, outputs) = run(json!({
            "input": line_layer(&[
                (5.0, 50.0, 60.0), (35.0, 50.0, 20.0), (65.0, 50.0, 20.0), (95.0, 50.0, 60.0)
            ]),
            "surface": flat_surface(40.0), "spacing": 1.0
        }));
        assert_eq!(outputs["crossing_count"].as_u64().unwrap(), 2);
        assert_eq!(layer.len(), 3);
        let positions: Vec<String> = layer.iter().map(|f| text(&layer, f, "position")).collect();
        assert_eq!(positions, vec!["above", "below", "above"]);

        // The crossing points are tagged with their direction of travel.
        let args: ToolArgs = serde_json::from_value(json!({
            "input": line_layer(&[
                (5.0, 50.0, 60.0), (35.0, 50.0, 20.0), (65.0, 50.0, 20.0), (95.0, 50.0, 60.0)
            ]),
            "surface": flat_surface(40.0), "spacing": 1.0
        }))
        .unwrap();
        let out = Intersect3dLineWithSurfaceTool.run(&args, &ctx()).unwrap();
        let pts = load_input_layer(out.outputs["output_points"].as_str().unwrap()).unwrap();
        let di = pts.schema.field_index("direction").unwrap();
        let dirs: Vec<String> = pts
            .iter()
            .map(|f| match &f.attributes[di] {
                FieldValue::Text(t) => t.clone(),
                _ => panic!(),
            })
            .collect();
        assert_eq!(dirs, vec!["entering", "emerging"]);
    }

    /// A line entirely above the surface is not cut at all.
    #[test]
    fn line_clear_of_the_surface_is_untouched() {
        let (layer, outputs) = run(json!({
            "input": line_layer(&[(5.0, 50.0, 200.0), (95.0, 50.0, 210.0)]),
            "surface": flat_surface(50.0), "spacing": 5.0
        }));
        assert_eq!(outputs["crossing_count"].as_u64().unwrap(), 0);
        assert_eq!(layer.len(), 1);
        assert_eq!(
            text(&layer, layer.iter().next().unwrap(), "position"),
            "above"
        );
    }

    /// Clearance statistics describe how far the piece runs from the surface —
    /// the quantity a pipeline or sightline study is actually after.
    #[test]
    fn reports_clearance_per_part() {
        let (layer, _) = run(json!({
            "input": line_layer(&[(5.0, 50.0, 60.0), (95.0, 50.0, 80.0)]),
            "surface": flat_surface(50.0), "spacing": 5.0
        }));
        let f = layer.iter().next().unwrap();
        let mn = match f.attributes[layer.schema.field_index("min_clearance").unwrap()] {
            FieldValue::Float(v) => v,
            _ => panic!(),
        };
        let mx = match f.attributes[layer.schema.field_index("max_clearance").unwrap()] {
            FieldValue::Float(v) => v,
            _ => panic!(),
        };
        assert!((mn - 10.0).abs() < 1.0, "min clearance {mn} should be ~10");
        assert!((mx - 30.0).abs() < 1.0, "max clearance {mx} should be ~30");
    }

    /// The keep filters select which side survives.
    #[test]
    fn keep_filters_select_a_side() {
        let src = json!({
            "input": line_layer(&[(5.0, 50.0, 0.0), (95.0, 50.0, 100.0)]),
            "surface": flat_surface(50.0), "spacing": 2.0
        });
        let mut above_only = src.as_object().unwrap().clone();
        above_only.insert("keep_below".into(), json!(false));
        let (layer, _) = run(Value::Object(above_only));
        assert_eq!(layer.len(), 1);
        assert_eq!(
            text(&layer, layer.iter().next().unwrap(), "position"),
            "above"
        );
    }

    /// The multipatch form cuts against a mesh instead of a raster — the
    /// *Intersect 3D Line With Multi Patch* case.
    #[test]
    fn cuts_against_a_multipatch() {
        // A flat square slab at z = 30, spanning x,y in [0, 100].
        let mut mesh = Layer::new("slab");
        mesh.geom_type = Some(GeometryType::MultiPolygon);
        mesh = mesh.with_crs_epsg(32610);
        mesh.add_field(FieldDef::new("id", FieldType::Integer));
        let tri = |a: (f64, f64), b: (f64, f64), c: (f64, f64)| {
            Ring::new(vec![
                Coord::xyz(a.0, a.1, 30.0),
                Coord::xyz(b.0, b.1, 30.0),
                Coord::xyz(c.0, c.1, 30.0),
            ])
        };
        // A MultiPolygon part is an (exterior, interiors) tuple.
        let parts = vec![
            (tri((0.0, 0.0), (100.0, 0.0), (100.0, 100.0)), Vec::new()),
            (tri((0.0, 0.0), (100.0, 100.0), (0.0, 100.0)), Vec::new()),
        ];
        let mut f = Feature::with_geometry(0, Geometry::MultiPolygon(parts), mesh.schema.len());
        f.set_by_index(0, FieldValue::Integer(0));
        mesh.push(f);
        let mesh_path = write_or_store_layer(mesh, None).unwrap();

        let (layer, outputs) = run(json!({
            "input": line_layer(&[(10.0, 50.0, 0.0), (90.0, 50.0, 60.0)]),
            "multipatch": mesh_path, "spacing": 1.0
        }));
        assert_eq!(outputs["crossing_count"].as_u64().unwrap(), 1);
        assert_eq!(layer.len(), 2);
        let positions: Vec<String> = layer.iter().map(|f| text(&layer, f, "position")).collect();
        assert_eq!(positions, vec!["below", "above"]);
    }

    /// Crossing points are emitted even without a path.
    #[test]
    fn crossing_points_emitted_without_a_path() {
        let args: ToolArgs = serde_json::from_value(json!({
            "input": line_layer(&[(5.0, 50.0, 0.0), (95.0, 50.0, 100.0)]),
            "surface": flat_surface(50.0), "spacing": 2.0
        }))
        .unwrap();
        let out = Intersect3dLineWithSurfaceTool.run(&args, &ctx()).unwrap();
        let path = out.outputs["output_points"].as_str().unwrap();
        assert!(!path.is_empty());
        assert_eq!(load_input_layer(path).unwrap().len(), 1);
    }

    /// Regression: `coords` carries the interpolated cut vertices, which have
    /// no clearance value, so taking `coords.len()` samples read one value from
    /// the other side of the surface. An `above` part could then report a
    /// negative minimum clearance and a `below` part a positive maximum.
    #[test]
    fn clearance_signs_match_each_part() {
        let (layer, _) = run(json!({
            "input": line_layer(&[
                (5.0, 50.0, 60.0), (35.0, 50.0, 20.0), (65.0, 50.0, 20.0), (95.0, 50.0, 60.0)
            ]),
            "surface": flat_surface(40.0), "spacing": 1.0
        }));
        assert_eq!(layer.len(), 3);
        let mi = layer.schema.field_index("min_clearance").unwrap();
        let xi = layer.schema.field_index("max_clearance").unwrap();
        for f in layer.iter() {
            let (FieldValue::Float(mn), FieldValue::Float(mx)) =
                (&f.attributes[mi], &f.attributes[xi])
            else {
                panic!("clearances must be floats")
            };
            match text(&layer, f, "position").as_str() {
                "above" => assert!(
                    *mn >= 0.0 && *mx >= 0.0,
                    "an above part must not report a negative clearance ({mn}, {mx})"
                ),
                "below" => assert!(
                    *mn <= 0.0 && *mx <= 0.0,
                    "a below part must not report a positive clearance ({mn}, {mx})"
                ),
                other => panic!("unexpected position {other}"),
            }
        }
    }

    /// A mesh with no plan extent yields no spacing, and must be reported
    /// rather than densified into an unbounded loop.
    #[test]
    fn degenerate_multipatch_is_reported() {
        let mut mesh = Layer::new("degenerate");
        mesh.geom_type = Some(GeometryType::MultiPolygon);
        mesh = mesh.with_crs_epsg(32610);
        mesh.add_field(FieldDef::new("id", FieldType::Integer));
        // Three vertices sharing one plan position: zero plan diagonal.
        let ring = Ring::new(vec![
            Coord::xyz(50.0, 50.0, 0.0),
            Coord::xyz(50.0, 50.0, 10.0),
            Coord::xyz(50.0, 50.0, 20.0),
        ]);
        let mut f = Feature::with_geometry(
            0,
            Geometry::MultiPolygon(vec![(ring, Vec::new())]),
            mesh.schema.len(),
        );
        f.set_by_index(0, FieldValue::Integer(0));
        mesh.push(f);
        let mesh_path = write_or_store_layer(mesh, None).unwrap();

        let args: ToolArgs = serde_json::from_value(json!({
            "input": line_layer(&[(10.0, 50.0, 0.0), (90.0, 50.0, 60.0)]),
            "multipatch": mesh_path
        }))
        .unwrap();
        let err = Intersect3dLineWithSurfaceTool
            .run(&args, &ctx())
            .unwrap_err();
        assert!(
            format!("{err:?}").contains("spacing"),
            "expected a spacing error, got {err:?}"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            Intersect3dLineWithSurfaceTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        // Neither cutting geometry supplied.
        assert!(bad(json!({"input": "l.shp"})).is_err());
        // Both supplied is ambiguous.
        assert!(bad(json!({"input": "l.shp", "surface": "s.tif", "multipatch": "m.shp"})).is_err());
        assert!(bad(json!({
            "input": "l.shp", "surface": "s.tif", "keep_above": false, "keep_below": false
        }))
        .is_err());
        assert!(bad(json!({"input": "l.shp", "surface": "s.tif", "spacing": -1})).is_err());
        assert!(bad(json!({"input": "l.shp", "surface": "s.tif"})).is_ok());
    }
}
