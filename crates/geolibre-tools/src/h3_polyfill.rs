//! GeoLibre tool: cover polygons with H3 cells (polyfill).
//!
//! For each polygon in the input layer, finds the H3 cells at a chosen
//! resolution that cover it and emits one polygon per distinct cell with its
//! `h3` index. Coverage uses [`ContainmentMode::Covers`], so the returned cells
//! fully blanket each input polygon (the complete-coverage choice, as opposed to
//! centroid-only binning).
//!
//! This complements [`crate::vector_to_h3`] (which bins features *into* cells by
//! a representative point): polyfill instead tiles an *area* with cells, the H3
//! analogue of a polygon-to-grid overlay.
//!
//! H3 work is done by [`h3o`](https://docs.rs/h3o) (pure Rust, no C deps); the
//! polygon tiler lives behind its `geo` feature. See issue #22.

use std::collections::BTreeSet;

use h3o::geom::{ContainmentMode, TilerBuilder};
use h3o::{CellIndex, Resolution};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Default H3 resolution (see [`crate::vector_to_h3`]); 8 is ~0.7 km² cells.
const DEFAULT_RESOLUTION: u8 = 8;

/// Covers input polygons with H3 cells at a chosen resolution.
pub struct H3PolyfillTool;

impl Tool for H3PolyfillTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "h3_polyfill",
            display_name: "H3 Polyfill",
            summary: "Cover input polygons with H3 discrete-global-grid cells at a chosen resolution.",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon vector file path (must be geographic lon/lat, EPSG:4326) or in-memory handle.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path; the driver is taken from its extension (.geojson, .fgb, .parquet, ...). If omitted, the layer is stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "resolution",
                    description: "H3 resolution 0 (coarsest) to 15 (finest); default 8 (~0.7 km² cells).",
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
        if let Some(res) = parse_optional_resolution(args)? {
            Resolution::try_from(res).map_err(|_| {
                ToolError::Validation(format!("'resolution' must be 0-15, got {res}"))
            })?;
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args.get("input").and_then(Value::as_str).ok_or_else(|| {
            ToolError::Validation("missing required parameter 'input'".to_string())
        })?;
        let output = parse_optional_str(args, "output")?;
        let res_u8 = parse_optional_resolution(args)?.unwrap_or(DEFAULT_RESOLUTION);
        let resolution = Resolution::try_from(res_u8).map_err(|_| {
            ToolError::Validation(format!("'resolution' must be 0-15, got {res_u8}"))
        })?;

        ctx.progress.info("reading input vector");
        let layer = load_input_layer(input)?;
        // H3 indexes geographic coordinates; refuse a non-4326 CRS rather than
        // silently producing garbage cells. A missing CRS is assumed to be 4326.
        if let Some(epsg) = layer.crs_epsg() {
            if epsg != 4326 {
                return Err(ToolError::Validation(format!(
                    "input CRS is EPSG:{epsg}; reproject to EPSG:4326 (lon/lat) before H3 polyfill"
                )));
            }
        }
        let feature_count = layer.len();

        ctx.progress.info("covering polygons with H3 cells");
        let mut cells: BTreeSet<u64> = BTreeSet::new();
        let mut skipped = 0u64;
        for feature in layer.iter() {
            let Some(geom) = feature.geometry.as_ref() else {
                skipped += 1;
                continue;
            };
            let polys = geo_polygons(geom);
            if polys.is_empty() {
                skipped += 1; // not an areal geometry
                continue;
            }
            for poly in polys {
                // A degenerate/invalid ring fails to tile; skip it rather than
                // aborting the whole run.
                let mut tiler = TilerBuilder::new(resolution)
                    .containment_mode(ContainmentMode::Covers)
                    .build();
                if tiler.add(poly).is_err() {
                    skipped += 1;
                    continue;
                }
                cells.extend(tiler.into_coverage().map(u64::from));
            }
        }

        ctx.progress.info("building H3 cell polygons");
        let mut out = Layer::new("h3_polyfill")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        out.add_field(FieldDef::new("h3", FieldType::Text));
        for &cell_raw in &cells {
            let cell = CellIndex::try_from(cell_raw)
                .map_err(|e| ToolError::Execution(format!("invalid H3 cell: {e}")))?;
            out.add_feature(
                Some(Geometry::polygon(cell_polygon_ring(cell), Vec::new())),
                &[("h3", FieldValue::Text(cell.to_string()))],
            )
            .map_err(|e| ToolError::Execution(format!("failed building H3 feature: {e}")))?;
        }

        ctx.progress.info("writing output vector");
        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("resolution".to_string(), json!(res_u8));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("cell_count".to_string(), json!(cells.len()));
        outputs.insert("skipped".to_string(), json!(skipped));
        Ok(ToolRunResult { outputs })
    }
}

/// Parses the optional `resolution` parameter (a JSON number or numeric string).
fn parse_optional_resolution(args: &ToolArgs) -> Result<Option<u8>, ToolError> {
    match args.get("resolution") {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s.trim().parse::<u8>().map(Some).map_err(|_| {
            ToolError::Validation(format!("'resolution' must be an integer, got '{s}'"))
        }),
        Some(Value::Number(n)) => n
            .as_u64()
            .filter(|v| *v <= u8::MAX as u64)
            .map(|v| Some(v as u8))
            .ok_or_else(|| {
                ToolError::Validation("'resolution' must be an integer 0-15".to_string())
            }),
        Some(_) => Err(ToolError::Validation(
            "'resolution' must be a number when provided".to_string(),
        )),
    }
}

/// Converts an input geometry into the `geo::Polygon`s the tiler accepts. Only
/// (multi)polygons contribute; point/line geometries yield an empty list and
/// are counted as skipped by the caller.
fn geo_polygons(geom: &Geometry) -> Vec<geo::Polygon> {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => vec![to_geo_polygon(exterior, interiors)],
        Geometry::MultiPolygon(polys) => polys
            .iter()
            .map(|(ext, holes)| to_geo_polygon(ext, holes))
            .collect(),
        Geometry::GeometryCollection(geoms) => geoms.iter().flat_map(geo_polygons).collect(),
        _ => Vec::new(),
    }
}

fn to_geo_polygon(exterior: &Ring, interiors: &[Ring]) -> geo::Polygon {
    geo::Polygon::new(
        ring_to_linestring(exterior),
        interiors.iter().map(ring_to_linestring).collect(),
    )
}

fn ring_to_linestring(ring: &Ring) -> geo::LineString {
    // `geo::Polygon::new` closes rings itself, so the missing closing vertex in
    // `Ring` is fine.
    geo::LineString::new(
        ring.coords()
            .iter()
            .map(|c| geo::Coord { x: c.x, y: c.y })
            .collect(),
    )
}

/// Builds the exterior ring (lon/lat coordinates) of an H3 cell's boundary,
/// without the closing duplicate vertex (matching `Ring`).
fn cell_polygon_ring(cell: CellIndex) -> Vec<Coord> {
    cell.boundary()
        .iter()
        .map(|ll| Coord::xy(ll.lng(), ll.lat()))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;
    use wbcore::{AllowAllCapabilities, ProgressSink};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// A ~0.2° square around San Francisco (lon/lat), CCW, unclosed ring.
    fn square_layer() -> Layer {
        let mut l = Layer::new("poly")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        let ring = vec![
            Coord::xy(-122.5, 37.7),
            Coord::xy(-122.3, 37.7),
            Coord::xy(-122.3, 37.9),
            Coord::xy(-122.5, 37.9),
        ];
        l.add_feature(Some(Geometry::polygon(ring, Vec::new())), &[])
            .unwrap();
        l
    }

    #[test]
    fn covers_polygon_with_cells() {
        let id = wbvector::memory_store::put_vector(square_layer());
        let input = wbvector::memory_store::make_vector_memory_path(&id);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "resolution": "7" })).unwrap();

        let result = H3PolyfillTool.run(&args, &ctx()).unwrap();
        assert_eq!(result.outputs["feature_count"], json!(1));
        assert_eq!(result.outputs["skipped"], json!(0));
        let n = result.outputs["cell_count"].as_u64().unwrap();
        assert!(n > 1, "a ~0.2-degree square should need many res-7 cells, got {n}");

        let out = load_input_layer(result.outputs["output"].as_str().unwrap()).unwrap();
        assert_eq!(out.len() as u64, n);
        assert_eq!(out.geom_type, Some(GeometryType::Polygon));
        // Every output cell id parses back to a res-7 cell.
        let hidx = out.schema.field_index("h3").unwrap();
        for f in out.iter() {
            let id = f.attributes[hidx].as_str().unwrap();
            let cell = h3o::CellIndex::from_str(id).unwrap();
            assert_eq!(cell.resolution(), Resolution::Seven);
        }
    }

    #[test]
    fn finer_resolution_yields_more_cells() {
        let run = |res: &str| {
            let id = wbvector::memory_store::put_vector(square_layer());
            let input = wbvector::memory_store::make_vector_memory_path(&id);
            let args: ToolArgs =
                serde_json::from_value(json!({ "input": input, "resolution": res })).unwrap();
            H3PolyfillTool.run(&args, &ctx()).unwrap().outputs["cell_count"]
                .as_u64()
                .unwrap()
        };
        assert!(run("8") > run("6"), "finer resolution must yield more cells");
    }

    #[test]
    fn non_polygon_features_are_skipped() {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        l.add_feature(Some(Geometry::point(-122.4, 37.8)), &[]).unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let input = wbvector::memory_store::make_vector_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({ "input": input })).unwrap();
        let r = H3PolyfillTool.run(&args, &ctx()).unwrap();
        assert_eq!(r.outputs["skipped"], json!(1));
        assert_eq!(r.outputs["cell_count"], json!(0));
    }
}
