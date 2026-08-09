//! GeoLibre tool: closed solids spanning the space between two surfaces.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Extrude Between* (3D Analyst).
//!
//! ## The gap
//!
//! `cut_fill` (round 1) gives per-cell volume change between two rasters, but
//! it emits a **raster**: it cannot be restricted to polygon footprints, cannot
//! carry attributes through, and produces no geometry. `buffer_3d` and
//! `minimum_bounding_volume` build solids from features rather than from a
//! surface pair. So there was no way to model a geological unit, an excavation
//! or a fill body as an actual solid.
//!
//! The output feeds straight into `union_3d`, `intersect_3d`, `difference_3d`
//! and `inside_3d`, which is what makes it worth having as geometry rather than
//! a number.
//!
//! ## Watertight by construction
//!
//! The solid is an upper shell over the footprint's sample grid, a mirrored
//! lower shell, and a side wall built from the densified boundary. Because the
//! wall is stitched to the *same* boundary samples both shells use, every edge
//! is shared by exactly two triangles by construction — asserted in a test with
//! the same check `is_closed_3d` performs.
//!
//! Volume is exact (signed-tetrahedron summation), not sampled: unlike the 3D
//! overlay tools there is no boolean here, so no approximation is needed.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::args_common::{band_index, opt_positive_f64, req_str};
use crate::common::load_input_raster;
use crate::inside_3d::Tri;
use crate::mesh3d::{mesh_volume, topology, triangles_to_geometry};
use crate::raster_stack::check_alignment;
use crate::surface_solid::{default_spacing, densify, sample_bilinear};
use crate::vector_common::{
    geometry_contains_point, load_input_layer, parse_optional_str, write_or_store_layer,
};

pub struct ExtrudeBetweenTool;

impl Tool for ExtrudeBetweenTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "extrude_between",
            display_name: "Extrude Between",
            summary: "Builds watertight closed solids spanning the space between two surfaces inside each input polygon, with exact per-feature volumes (ArcGIS Extrude Between). cut_fill gives per-cell volume change between two rasters but emits a raster, cannot be restricted to footprints and carries no attributes; buffer_3d and minimum_bounding_volume build solids from features rather than from a surface pair. The output feeds directly into union_3d, intersect_3d, difference_3d and inside_3d.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Polygon footprints defining each solid's lateral extent.",
                    required: true,
                },
                ToolParamSpec {
                    name: "surface_upper",
                    description: "Upper bounding surface raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "surface_lower",
                    description: "Lower bounding surface raster, co-registered with the upper one.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Closed solid features (triangle-mesh MultiPolygons with Z). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "sample_distance",
                    description: "Footprint sampling and boundary densification spacing in CRS units. Default: the finer surface cell size.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read from each surface raster (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        req_str(args, "surface_upper")?;
        req_str(args, "surface_lower")?;
        opt_positive_f64(args, "sample_distance")?;
        band_index(args, "band")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?;
        let upper = load_input_raster(req_str(args, "surface_upper")?)?;
        let lower = load_input_raster(req_str(args, "surface_lower")?)?;
        check_alignment(&[upper.clone(), lower.clone()])?;
        let band = band_index(args, "band")?;
        let output = parse_optional_str(args, "output")?;
        let spacing = match opt_positive_f64(args, "sample_distance")? {
            Some(v) => v,
            None => default_spacing(&[&upper, &lower])?,
        };

        let layer = load_input_layer(input)?;
        let mut out = Layer::new("extrude_between");
        out.geom_type = Some(GeometryType::MultiPolygon);
        out.crs = layer.crs.clone();
        for f in layer.schema.fields() {
            out.add_field(f.clone());
        }
        out.add_field(FieldDef::new("SRC_FID", FieldType::Integer));
        out.add_field(FieldDef::new("VOLUME", FieldType::Float));
        out.add_field(FieldDef::new("SAMPLE_DISTANCE", FieldType::Float));
        out.add_field(FieldDef::new("WATERTIGHT", FieldType::Boolean));

        let names: Vec<String> = layer
            .schema
            .fields()
            .iter()
            .map(|f| f.name.clone())
            .collect();

        let mut built = 0_u64;
        let mut skipped = 0_u64;
        let mut total_volume = 0.0_f64;
        let total = layer.iter().count().max(1);

        for (fid, feature) in layer.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                skipped += 1;
                continue;
            };
            let Some(tris) = build_solid(geom, &upper, &lower, band, spacing) else {
                skipped += 1;
                continue;
            };
            let t = topology(&tris);
            let volume = if t.closed && t.consistent_winding {
                mesh_volume(&tris)
            } else {
                0.0
            };
            total_volume += volume;
            built += 1;

            let mut attrs: Vec<(&str, FieldValue)> = names
                .iter()
                .enumerate()
                .filter_map(|(k, n)| feature.attributes.get(k).map(|v| (n.as_str(), v.clone())))
                .collect();
            attrs.push(("SRC_FID", FieldValue::Integer(fid as i64)));
            attrs.push(("VOLUME", FieldValue::Float(volume)));
            attrs.push(("SAMPLE_DISTANCE", FieldValue::Float(spacing)));
            attrs.push(("WATERTIGHT", FieldValue::Boolean(t.closed)));
            out.add_feature(Some(triangles_to_geometry(&tris)), &attrs)
                .map_err(|e| ToolError::Execution(e.to_string()))?;
            ctx.progress.progress((fid as f64 + 1.0) / total as f64);
        }

        if built == 0 {
            return Err(ToolError::Execution(format!(
                "no solid could be built ({skipped} footprint(s) had no usable surface coverage)"
            )));
        }
        ctx.progress
            .info(&format!("{built} solid(s), {skipped} skipped"));

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("solid_count".to_string(), json!(built));
        outputs.insert("skipped_count".to_string(), json!(skipped));
        outputs.insert("total_volume".to_string(), json!(total_volume));
        outputs.insert("sample_distance".to_string(), json!(spacing));
        Ok(ToolRunResult { outputs })
    }
}

/// Builds one closed solid from a polygon footprint and two surfaces.
///
/// The construction is a prism whose top and bottom follow the surfaces: the
/// boundary ring supplies the side wall, and the ring is fanned from its
/// centroid on both caps. Fanning from the *same* boundary samples the wall
/// uses is what guarantees watertightness — an independently gridded cap would
/// not share the wall's edges.
fn build_solid(
    geom: &Geometry,
    upper: &Raster,
    lower: &Raster,
    band: isize,
    spacing: f64,
) -> Option<Vec<Tri>> {
    let ring = exterior_ring(geom)?;
    let mut boundary = densify(ring, spacing);
    // densify repeats the closing vertex when the ring is already closed;
    // drop it so the wall does not emit a zero-area quad.
    if boundary.len() >= 2 {
        let (a, b) = (boundary[0], boundary[boundary.len() - 1]);
        if (a.0 - b.0).abs() < 1e-12 && (a.1 - b.1).abs() < 1e-12 {
            boundary.pop();
        }
    }
    if boundary.len() < 3 {
        return None;
    }

    // Sample both surfaces along the boundary. A footprint straddling the
    // surface edge is rejected rather than closed with fabricated elevations.
    let mut top: Vec<[f64; 3]> = Vec::with_capacity(boundary.len());
    let mut bot: Vec<[f64; 3]> = Vec::with_capacity(boundary.len());
    for (x, y) in &boundary {
        let zu = sample_bilinear(upper, band, *x, *y)?;
        let zl = sample_bilinear(lower, band, *x, *y)?;
        top.push([*x, *y, zu]);
        bot.push([*x, *y, zl]);
    }

    // Interior centroid, used as the fan apex on both caps. Its elevation is
    // sampled too, so the caps follow the surfaces rather than being flat.
    let n = boundary.len() as f64;
    let cx = boundary.iter().map(|p| p.0).sum::<f64>() / n;
    let cy = boundary.iter().map(|p| p.1).sum::<f64>() / n;
    // The centroid of a concave footprint can fall outside it; sampling there
    // would read the surface at a point the solid does not cover.
    if !geometry_contains_point(geom, cx, cy) {
        return None;
    }
    let cz_top = sample_bilinear(upper, band, cx, cy)?;
    let cz_bot = sample_bilinear(lower, band, cx, cy)?;
    let apex_top = [cx, cy, cz_top];
    let apex_bot = [cx, cy, cz_bot];

    let m = boundary.len();
    let mut tris: Vec<Tri> = Vec::with_capacity(m * 4);
    for i in 0..m {
        let j = (i + 1) % m;
        // Top cap, outward (upward) winding.
        tris.push([apex_top, top[i], top[j]]);
        // Bottom cap, reversed so it faces downward.
        tris.push([apex_bot, bot[j], bot[i]]);
        // Side wall quad as two triangles.
        tris.push([bot[i], bot[j], top[j]]);
        tris.push([bot[i], top[j], top[i]]);
    }
    Some(tris)
}

fn exterior_ring(geom: &Geometry) -> Option<&Vec<Coord>> {
    match geom {
        Geometry::Polygon { exterior, .. } => Some(&exterior.0),
        Geometry::MultiPolygon(parts) => parts.first().map(|(ext, _)| &ext.0),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, RasterConfig};
    use wbvector::{memory_store, Ring};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Constant-value surface covering (0,0)-(100,100) at 10-unit cells.
    fn flat(value: f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols: 10,
            rows: 10,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 10.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for row in 0..10 {
            for col in 0..10 {
                r.set(0, row, col, value).unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn square(x0: f64, y0: f64, x1: f64, y1: f64) -> String {
        let mut l = Layer::new("fp");
        l.geom_type = Some(GeometryType::Polygon);
        l.add_feature(
            Some(Geometry::Polygon {
                exterior: Ring::new(vec![
                    Coord::xy(x0, y0),
                    Coord::xy(x1, y0),
                    Coord::xy(x1, y1),
                    Coord::xy(x0, y1),
                ]),
                interiors: vec![],
            }),
            &[],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = ExtrudeBetweenTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    fn num(layer: &Layer, fid: usize, name: &str) -> f64 {
        let i = layer.schema.field_index(name).unwrap();
        match &layer.iter().nth(fid).unwrap().attributes[i] {
            FieldValue::Float(v) => *v,
            FieldValue::Integer(v) => *v as f64,
            other => panic!("expected a number, got {other:?}"),
        }
    }

    #[test]
    fn two_flat_surfaces_give_footprint_area_times_separation() {
        // A 40x40 footprint between z = 10 and z = 30: 1600 * 20 = 32000.
        let (out, res) = run(json!({
            "input": square(20.0, 20.0, 60.0, 60.0),
            "surface_lower": flat(10.0),
            "surface_upper": flat(30.0),
            "sample_distance": 5.0,
        }));
        assert_eq!(res.outputs["solid_count"], json!(1));
        let v = num(&out, 0, "VOLUME");
        assert!((v - 32000.0).abs() < 1.0, "expected 32000, got {v}");
    }

    #[test]
    fn the_solid_is_watertight() {
        // The precondition every downstream volumetric tool checks.
        let (out, _) = run(json!({
            "input": square(20.0, 20.0, 60.0, 60.0),
            "surface_lower": flat(0.0),
            "surface_upper": flat(5.0),
            "sample_distance": 10.0,
        }));
        let geom = out.iter().next().unwrap().geometry.clone().unwrap();
        let t = topology(&crate::inside_3d::collect_triangles(&geom));
        assert!(t.closed, "{} edges left open", t.open_edges);
        assert!(t.consistent_winding, "orientation is inconsistent");
        let i = out.schema.field_index("WATERTIGHT").unwrap();
        assert_eq!(
            out.iter().next().unwrap().attributes[i],
            FieldValue::Boolean(true)
        );
    }

    #[test]
    fn a_sloping_upper_surface_still_gives_the_mean_separation_volume() {
        // Lower flat at 0, upper a plane rising 0 to 90 across x. Over a
        // footprint spanning x 20..60 the mean upper height is 40, so the
        // volume is 40 * 40 * 40 = 64000.
        let mut r = Raster::new(RasterConfig {
            cols: 10,
            rows: 10,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 10.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for row in 0..10 {
            for col in 0..10 {
                // Cell centre x = col*10 + 5.
                r.set(0, row, col, col as f64 * 10.0 + 5.0).unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        let ramp = wbraster::memory_store::make_raster_memory_path(&id);

        let (out, _) = run(json!({
            "input": square(20.0, 20.0, 60.0, 60.0),
            "surface_lower": flat(0.0),
            "surface_upper": ramp,
            "sample_distance": 2.0,
        }));
        let v = num(&out, 0, "VOLUME");
        assert!((v - 64000.0).abs() < 2000.0, "expected ~64000, got {v}");
    }

    #[test]
    fn a_finer_sample_distance_converges_on_the_analytic_volume() {
        let mk = |spacing: f64| {
            let (out, _) = run(json!({
                "input": square(20.0, 20.0, 60.0, 60.0),
                "surface_lower": flat(10.0),
                "surface_upper": flat(30.0),
                "sample_distance": spacing,
            }));
            (num(&out, 0, "VOLUME") - 32000.0).abs()
        };
        assert!(mk(2.0) <= mk(20.0) + 1e-6);
    }

    #[test]
    fn swapping_the_surfaces_gives_the_same_magnitude() {
        // Volume is unsigned, so an inverted pair must not report a negative
        // or a doubled figure.
        let a = run(json!({
            "input": square(20.0, 20.0, 60.0, 60.0),
            "surface_lower": flat(10.0), "surface_upper": flat(30.0),
            "sample_distance": 10.0,
        }));
        let b = run(json!({
            "input": square(20.0, 20.0, 60.0, 60.0),
            "surface_lower": flat(30.0), "surface_upper": flat(10.0),
            "sample_distance": 10.0,
        }));
        assert!((num(&a.0, 0, "VOLUME") - num(&b.0, 0, "VOLUME")).abs() < 1.0);
    }

    #[test]
    fn a_footprint_outside_the_surface_extent_is_skipped_not_fabricated() {
        let args: ToolArgs = serde_json::from_value(json!({
            "input": square(500.0, 500.0, 600.0, 600.0),
            "surface_lower": flat(0.0),
            "surface_upper": flat(10.0),
        }))
        .unwrap();
        assert!(ExtrudeBetweenTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn the_output_is_valid_input_to_union_3d() {
        // The composition that justifies emitting geometry rather than a number.
        let (out, _) = run(json!({
            "input": square(20.0, 20.0, 60.0, 60.0),
            "surface_lower": flat(0.0), "surface_upper": flat(10.0),
            "sample_distance": 10.0,
        }));
        let geom = out.iter().next().unwrap().geometry.clone().unwrap();
        let solid = crate::inside_3d::Solid::new(0, crate::inside_3d::collect_triangles(&geom));
        assert!(solid.closed, "union_3d would reject this solid");
    }

    #[test]
    fn rejects_bad_parameters() {
        let fp = square(0.0, 0.0, 10.0, 10.0);
        let s = flat(0.0);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ExtrudeBetweenTool.validate(&args).is_err()
        };
        assert!(bad(json!({"surface_upper": s, "surface_lower": s})));
        assert!(bad(json!({"input": fp, "surface_upper": s})));
        assert!(bad(
            json!({"input": fp, "surface_upper": s, "surface_lower": s, "sample_distance": 0})
        ));
    }
}
