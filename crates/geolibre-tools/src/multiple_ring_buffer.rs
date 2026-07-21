//! GeoLibre tool: concentric buffers at multiple distances.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Multiple Ring Buffer* (Analysis):
//! buffer the input features at a list of distances and return either nested
//! disks or non-overlapping rings (donuts), one polygon per distance band. The
//! everyday proximity-zone tool — drive-up catchments, impact bands, setback
//! zones — where whitebox-wasm ships only a single-distance buffer.
//!
//! For each distance `d` the features are buffered with `geo`'s `Buffer`
//! (i_overlay backend, same as `BooleanOps`; rounded joins, no GDAL/GEOS). With
//! `ring_type = disks` each band is the full buffer out to `d` (nested,
//! overlapping disks). With `ring_type = rings` (default) each band is the
//! buffer at `d` minus the buffer at the previous distance, giving
//! non-overlapping rings; for polygon inputs the innermost ring also excludes
//! the source polygon area (outside-only), so the rings tile the space around
//! the features.
//!
//! `dissolve = per_ring` unions every feature's band at a given distance into a
//! single polygon (the classic dissolved concentric-ring output); `dissolve =
//! none` (default) keeps one band polygon per source feature. Every output
//! feature carries the band distance in `distance_field` (default `distance`).
//! Input geometry of any type (points, lines, polygons, and their multi
//! variants) is accepted; nested `GeometryCollection`s are skipped.
//!
//! Scope for v1: source attributes are not copied onto the output bands (each
//! band carries only its distance) and distances are interpreted directly in
//! the layer's CRS units (no on-the-fly unit conversion).

use std::collections::BTreeMap;

use geo::{
    BooleanOps, Buffer, Coord as GeoCoord, Geometry as GeoGeometry, HasDimensions, LineString,
    MultiLineString, MultiPoint, MultiPolygon, Point, Polygon,
};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct MultipleRingBufferTool;

impl Tool for MultipleRingBufferTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "multiple_ring_buffer",
            display_name: "Multiple Ring Buffer",
            summary: "Buffer features at a list of distances into concentric bands — nested disks or non-overlapping rings — optionally dissolved per distance, with the band distance stored as an attribute.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (points, lines, or polygons), format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distances",
                    description: "Comma-separated list of buffer distances in CRS units, e.g. \"100,250,500\". Sorted ascending; duplicates and non-positive values are dropped.",
                    required: true,
                },
                ToolParamSpec {
                    name: "ring_type",
                    description: "'rings' (default) for non-overlapping donut bands, or 'disks' for nested cumulative buffers.",
                    required: false,
                },
                ToolParamSpec {
                    name: "dissolve",
                    description: "'none' (default) keeps one band per source feature; 'per_ring' unions all features' bands at each distance into one polygon.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distance_field",
                    description: "Name of the output attribute holding each band's distance. Default 'distance'.",
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

        let layer = load_input_layer(input)?;
        let layer_name = layer.name.clone();
        let layer_crs = layer.crs.clone();
        let input_count = layer.len();

        // Convert each usable source feature to a `geo` geometry once.
        let sources: Vec<GeoGeometry> = layer
            .features
            .iter()
            .filter_map(|f| f.geometry.as_ref().and_then(to_geo_geometry))
            .collect();

        ctx.progress.info(&format!(
            "{input_count} feature(s): buffering {} source geometr(ies) at {} distance(s)",
            sources.len(),
            prm.distances.len()
        ));

        let mut out_layer = Layer::new(layer_name);
        out_layer.crs = layer_crs;
        out_layer.add_field(FieldDef::new(&prm.distance_field, FieldType::Float));
        out_layer.geom_type = Some(GeometryType::MultiPolygon);

        let mut band_count = 0usize;
        if prm.dissolve == Dissolve::PerRing {
            // Union all sources at each distance, then difference successive
            // dissolved disks for rings.
            let mut prev = MultiPolygon::<f64>::new(vec![]);
            for &d in &prm.distances {
                let disk = union_all(sources.iter().map(|g| g.buffer(d)));
                let band = match prm.ring_type {
                    RingType::Disks => disk.clone(),
                    RingType::Rings => disk.difference(&prev),
                };
                if prm.ring_type == RingType::Rings {
                    prev = disk;
                }
                if emit_band(&mut out_layer, &band, &prm.distance_field, d)? {
                    band_count += 1;
                }
            }
        } else {
            // One independent set of bands per source feature.
            for g in &sources {
                let base = polygonal_base(g);
                let mut prev = base;
                for &d in &prm.distances {
                    let disk = g.buffer(d);
                    let band = match prm.ring_type {
                        RingType::Disks => disk.clone(),
                        RingType::Rings => disk.difference(&prev),
                    };
                    if prm.ring_type == RingType::Rings {
                        prev = disk;
                    }
                    if emit_band(&mut out_layer, &band, &prm.distance_field, d)? {
                        band_count += 1;
                    }
                }
            }
        }

        ctx.progress
            .info(&format!("wrote {band_count} buffer band(s)"));

        let feature_count = out_layer.len();
        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_count".to_string(), json!(input_count));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("band_count".to_string(), json!(band_count));
        outputs.insert("distances".to_string(), json!(prm.distances));
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum RingType {
    Rings,
    Disks,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Dissolve {
    None,
    PerRing,
}

struct Params {
    distances: Vec<f64>,
    ring_type: RingType,
    dissolve: Dissolve,
    distance_field: String,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let raw = parse_optional_str(args, "distances")?.ok_or_else(|| {
        ToolError::Validation("missing required parameter 'distances'".to_string())
    })?;
    let mut distances: Vec<f64> = Vec::new();
    for tok in raw.split(',') {
        let t = tok.trim();
        if t.is_empty() {
            continue;
        }
        let v: f64 = t
            .parse()
            .map_err(|_| ToolError::Validation(format!("distance '{t}' is not a number")))?;
        if v > 0.0 && v.is_finite() {
            distances.push(v);
        }
    }
    distances.sort_by(f64::total_cmp);
    distances.dedup();
    if distances.is_empty() {
        return Err(ToolError::Validation(
            "parameter 'distances' must list at least one positive distance".to_string(),
        ));
    }
    let ring_type = match parse_optional_str(args, "ring_type")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("rings") => RingType::Rings,
        Some("disks") => RingType::Disks,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown ring_type '{other}' (expected rings or disks)"
            )))
        }
    };
    let dissolve = match parse_optional_str(args, "dissolve")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("none") => Dissolve::None,
        Some("per_ring") => Dissolve::PerRing,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown dissolve '{other}' (expected none or per_ring)"
            )))
        }
    };
    let distance_field = parse_optional_str(args, "distance_field")?
        .map(str::to_string)
        .unwrap_or_else(|| "distance".to_string());
    Ok(Params {
        distances,
        ring_type,
        dissolve,
        distance_field,
    })
}

// ── Geometry helpers ──────────────────────────────────────────────────────────

/// Unions an iterator of `MultiPolygon`s into one.
fn union_all(parts: impl Iterator<Item = MultiPolygon>) -> MultiPolygon {
    let mut acc = MultiPolygon::<f64>::new(vec![]);
    for p in parts {
        acc = acc.union(&p);
    }
    acc
}

/// The polygonal footprint of a source geometry (for outside-only rings). Empty
/// for points and lines, which have no area to exclude.
fn polygonal_base(g: &GeoGeometry) -> MultiPolygon {
    match g {
        GeoGeometry::Polygon(p) => MultiPolygon(vec![p.clone()]),
        GeoGeometry::MultiPolygon(mp) => mp.clone(),
        _ => MultiPolygon::new(vec![]),
    }
}

/// Appends a band polygon (skipping empty ones) to the output layer.
fn emit_band(
    layer: &mut Layer,
    band: &MultiPolygon,
    field: &str,
    distance: f64,
) -> Result<bool, ToolError> {
    if band.is_empty() {
        return Ok(false);
    }
    layer
        .add_feature(
            Some(multipolygon_to_geometry(band)),
            &[(field, FieldValue::Float(distance))],
        )
        .map_err(|e| ToolError::Execution(format!("failed writing buffer band: {e}")))?;
    Ok(true)
}

/// Converts a `wbvector` geometry to a `geo` geometry for buffering. Nested
/// `GeometryCollection`s are skipped (returns `None`).
fn to_geo_geometry(g: &Geometry) -> Option<GeoGeometry> {
    let geom = match g {
        Geometry::Point(c) => GeoGeometry::Point(Point::new(c.x, c.y)),
        Geometry::MultiPoint(cs) => GeoGeometry::MultiPoint(MultiPoint(
            cs.iter().map(|c| Point::new(c.x, c.y)).collect(),
        )),
        Geometry::LineString(cs) => GeoGeometry::LineString(coords_to_linestring(cs)),
        Geometry::MultiLineString(ls) => GeoGeometry::MultiLineString(MultiLineString(
            ls.iter().map(|c| coords_to_linestring(c)).collect(),
        )),
        Geometry::Polygon {
            exterior,
            interiors,
        } => GeoGeometry::Polygon(rings_to_polygon(exterior, interiors)),
        Geometry::MultiPolygon(parts) => GeoGeometry::MultiPolygon(MultiPolygon(
            parts.iter().map(|(e, i)| rings_to_polygon(e, i)).collect(),
        )),
        Geometry::GeometryCollection(_) => return None,
    };
    Some(geom)
}

fn coords_to_linestring(coords: &[Coord]) -> LineString {
    LineString::new(coords.iter().map(|c| GeoCoord { x: c.x, y: c.y }).collect())
}

fn rings_to_polygon(exterior: &Ring, interiors: &[Ring]) -> Polygon {
    Polygon::new(
        ring_to_linestring(exterior),
        interiors.iter().map(ring_to_linestring).collect(),
    )
}

fn ring_to_linestring(ring: &Ring) -> LineString {
    LineString::new(
        ring.coords()
            .iter()
            .map(|c| GeoCoord { x: c.x, y: c.y })
            .collect(),
    )
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
    if coords.len() >= 2 && coords.first() == coords.last() {
        coords.pop();
    }
    Ring::new(coords)
}

#[cfg(test)]
mod tests {
    use super::*;
    use geo::Area;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn point_layer(pts: &[(f64, f64)]) -> String {
        let mut layer = Layer::new("pts");
        for &(x, y) in pts {
            layer.add_feature(Some(Geometry::point(x, y)), &[]).unwrap();
        }
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = MultipleRingBufferTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn dist_of(layer: &Layer, idx: usize) -> f64 {
        match layer.features[idx].get(&layer.schema, "distance").unwrap() {
            FieldValue::Float(f) => *f,
            other => panic!("distance should be float, got {other:?}"),
        }
    }

    fn area_of(layer: &Layer, idx: usize) -> f64 {
        to_geo_geometry(layer.features[idx].geometry.as_ref().unwrap())
            .and_then(|g| match g {
                GeoGeometry::Polygon(p) => Some(p.unsigned_area()),
                GeoGeometry::MultiPolygon(mp) => Some(mp.unsigned_area()),
                _ => None,
            })
            .unwrap()
    }

    #[test]
    fn rings_are_nonoverlapping_annuli() {
        let input = point_layer(&[(0.0, 0.0)]);
        let (out, layer) = run_tool(json!({ "input": input, "distances": "10,20,30" }));
        assert_eq!(out.outputs["band_count"], json!(3));
        assert_eq!(layer.len(), 3);
        // Distances ascending, one band each.
        for (i, d) in [10.0, 20.0, 30.0].iter().enumerate() {
            assert_eq!(dist_of(&layer, i), *d);
        }
        // Ring k area ~ pi*(dk^2 - d(k-1)^2): 100pi, 300pi, 500pi.
        let pi = std::f64::consts::PI;
        for (i, expect) in [100.0 * pi, 300.0 * pi, 500.0 * pi].iter().enumerate() {
            let a = area_of(&layer, i);
            assert!(
                (a - expect).abs() < 0.05 * expect,
                "ring {i} area {a} vs {expect}"
            );
        }
    }

    #[test]
    fn disks_are_cumulative() {
        let input = point_layer(&[(0.0, 0.0)]);
        let (_, layer) =
            run_tool(json!({ "input": input, "distances": "10,20", "ring_type": "disks" }));
        let pi = std::f64::consts::PI;
        // Disk areas ~ pi*d^2: 100pi then 400pi.
        assert!((area_of(&layer, 0) - 100.0 * pi).abs() < 0.05 * 100.0 * pi);
        assert!((area_of(&layer, 1) - 400.0 * pi).abs() < 0.05 * 400.0 * pi);
    }

    #[test]
    fn dissolve_per_ring_merges_overlapping_features() {
        // Two points 10 apart; at distance 20 their disks overlap. per_ring
        // dissolves them into a single band; none keeps two.
        let input = point_layer(&[(0.0, 0.0), (10.0, 0.0)]);
        let (none, none_layer) =
            run_tool(json!({ "input": input.clone(), "distances": "20", "dissolve": "none" }));
        let (diss, diss_layer) =
            run_tool(json!({ "input": input, "distances": "20", "dissolve": "per_ring" }));
        assert_eq!(
            none.outputs["band_count"],
            json!(2),
            "none -> one band per point"
        );
        assert_eq!(
            diss.outputs["band_count"],
            json!(1),
            "per_ring -> one merged band"
        );
        // Dissolved area is less than the sum of two disks (they overlap).
        let pi = std::f64::consts::PI;
        assert!(area_of(&diss_layer, 0) < 2.0 * 400.0 * pi);
        assert_eq!(none_layer.len(), 2);
    }

    #[test]
    fn polygon_rings_exclude_the_source_area() {
        // A 10x10 square buffered outward by 5: the first ring is the outside
        // band only, so its area excludes the 100-unit source.
        let mut layer = Layer::new("poly");
        layer
            .add_feature(
                Some(Geometry::polygon(
                    vec![
                        Coord::xy(0.0, 0.0),
                        Coord::xy(10.0, 0.0),
                        Coord::xy(10.0, 10.0),
                        Coord::xy(0.0, 10.0),
                    ],
                    vec![],
                )),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);
        let (_, layer) = run_tool(json!({ "input": input, "distances": "5" }));
        // Outside band area ~ perimeter*d + pi*d^2 = 40*5 + 25pi ~= 278.5; must
        // NOT include the 100-unit interior.
        let a = area_of(&layer, 0);
        let expect = 40.0 * 5.0 + std::f64::consts::PI * 25.0;
        assert!(
            (a - expect).abs() < 0.1 * expect,
            "outside ring area {a} vs {expect}"
        );
    }

    #[test]
    fn custom_distance_field_name() {
        let input = point_layer(&[(0.0, 0.0)]);
        let (_, layer) =
            run_tool(json!({ "input": input, "distances": "10", "distance_field": "band_m" }));
        assert!(layer.schema.field("band_m").is_some());
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = MultipleRingBufferTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "input": "x.geojson" })).is_err(),
            "missing distances"
        );
        assert!(bad(json!({ "input": "x.geojson", "distances": "" })).is_err());
        assert!(
            bad(json!({ "input": "x.geojson", "distances": "-5,0" })).is_err(),
            "no positive distance"
        );
        assert!(bad(json!({ "input": "x.geojson", "distances": "10,abc" })).is_err());
        assert!(
            bad(json!({ "input": "x.geojson", "distances": "10", "ring_type": "square" })).is_err()
        );
        assert!(
            bad(json!({ "input": "x.geojson", "distances": "10", "dissolve": "all" })).is_err()
        );
        assert!(bad(json!({ "input": "x.geojson", "distances": "10,20" })).is_ok());
    }
}
