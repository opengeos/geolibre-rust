//! GeoLibre tool: bin a vector layer into H3 cells (DGGS binning).
//!
//! Aggregates the features of an input layer into the cells of Uber's H3
//! discrete global grid at a chosen resolution, emitting one polygon per
//! occupied cell with a `h3` index column and a `count` of features that fell
//! in it. This is the discrete-global-grid analogue of the suite's existing
//! `vector_hex_binning`, but with stable, hierarchical H3 cell ids instead of an
//! arbitrary local hex grid.
//!
//! H3 work is done by [`h3o`](https://docs.rs/h3o), a pure-Rust reimplementation
//! of H3 with no C dependencies, so it compiles to the same wasm targets as the
//! rest of the suite. See issue #22.
//!
//! Each feature is reduced to a single representative lon/lat point (the mean of
//! its coordinates) before lookup, so points bin exactly and lines/polygons bin
//! by their centroid — adequate for density/heatmap-style binning. Polygon
//! coverage (polyfill) is a separate tool tracked in #22.

use std::collections::BTreeMap;

use h3o::{LatLng, Resolution};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Default H3 resolution (0 = coarsest ~1000 km, 15 = finest ~0.5 m). 8 is a
/// neighbourhood-scale cell (~0.7 km²), a sensible middle ground for binning.
const DEFAULT_RESOLUTION: u8 = 8;

/// Bins a vector layer into H3 cells and counts features per cell.
pub struct VectorToH3Tool;

impl Tool for VectorToH3Tool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "vector_to_h3",
            display_name: "Vector to H3 Bins",
            summary: "Bin a vector layer into H3 discrete-global-grid cells, counting features per cell.",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector file path (must be geographic lon/lat, EPSG:4326) or in-memory handle.",
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
                    "input CRS is EPSG:{epsg}; reproject to EPSG:4326 (lon/lat) before H3 binning"
                )));
            }
        }
        let feature_count = layer.len();

        ctx.progress.info("binning features into H3 cells");
        let mut counts: BTreeMap<u64, u64> = BTreeMap::new();
        let mut skipped = 0u64;
        for feature in layer.iter() {
            let Some(geom) = feature.geometry.as_ref() else {
                skipped += 1;
                continue;
            };
            let Some((lng, lat)) = representative_lnglat(geom) else {
                skipped += 1;
                continue;
            };
            match LatLng::new(lat, lng) {
                Ok(ll) => {
                    let cell = ll.to_cell(resolution);
                    *counts.entry(u64::from(cell)).or_insert(0) += 1;
                }
                Err(_) => skipped += 1, // out-of-range coordinate
            }
        }

        ctx.progress.info("building H3 cell polygons");
        let mut out = Layer::new("h3_bins")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        out.add_field(FieldDef::new("h3", FieldType::Text));
        out.add_field(FieldDef::new("count", FieldType::Integer));
        for (&cell_raw, &count) in &counts {
            let cell = h3o::CellIndex::try_from(cell_raw)
                .map_err(|e| ToolError::Execution(format!("invalid H3 cell: {e}")))?;
            let exterior = cell_polygon_ring(cell);
            out.add_feature(
                Some(Geometry::polygon(exterior, Vec::new())),
                &[
                    ("h3", cell.to_string().into()),
                    ("count", (count as i64).into()),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed building H3 feature: {e}")))?;
        }

        ctx.progress.info("writing output vector");
        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("resolution".to_string(), json!(res_u8));
        outputs.insert("cell_count".to_string(), json!(counts.len()));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("skipped".to_string(), json!(skipped));
        Ok(ToolRunResult { outputs })
    }
}

/// Parses the optional `resolution` parameter (accepts a JSON number or a
/// numeric string, as the CLI passes everything as strings).
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

/// Reduces a geometry to a single representative (lon, lat) point: the mean of
/// all its coordinates. Exact for points; the centroid-ish mean for everything
/// else, which is what binning needs. Returns `None` for empty geometries.
fn representative_lnglat(geom: &Geometry) -> Option<(f64, f64)> {
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut n = 0u64;
    accumulate_coords(geom, &mut sum_x, &mut sum_y, &mut n);
    (n > 0).then(|| (sum_x / n as f64, sum_y / n as f64))
}

fn accumulate_coords(geom: &Geometry, sx: &mut f64, sy: &mut f64, n: &mut u64) {
    let mut add = |c: &Coord| {
        *sx += c.x;
        *sy += c.y;
        *n += 1;
    };
    match geom {
        Geometry::Point(c) => add(c),
        Geometry::LineString(cs) | Geometry::MultiPoint(cs) => cs.iter().for_each(add),
        Geometry::MultiLineString(lines) => lines.iter().flatten().for_each(add),
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            exterior.coords().iter().for_each(&mut add);
            interiors
                .iter()
                .for_each(|r| r.coords().iter().for_each(&mut add));
        }
        Geometry::MultiPolygon(polys) => {
            for (ext, holes) in polys {
                ext.coords().iter().for_each(&mut add);
                holes
                    .iter()
                    .for_each(|r| r.coords().iter().for_each(&mut add));
            }
        }
        Geometry::GeometryCollection(geoms) => {
            for g in geoms {
                accumulate_coords(g, sx, sy, n);
            }
        }
    }
}

/// Builds the exterior ring (lon/lat coordinates) of an H3 cell's boundary.
/// The ring is stored without the closing duplicate vertex, matching `Ring`.
fn cell_polygon_ring(cell: h3o::CellIndex) -> Vec<Coord> {
    cell.boundary()
        .iter()
        .map(|ll| Coord::xy(ll.lng(), ll.lat()))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink, ToolContext};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn points_layer() -> Layer {
        let mut l = Layer::new("pts")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        // Three points: two essentially co-located (same H3 cell at res 8) and
        // one far away (a different cell).
        l.add_feature(Some(Geometry::point(-122.4194, 37.7749)), &[])
            .unwrap();
        l.add_feature(Some(Geometry::point(-122.4195, 37.7750)), &[])
            .unwrap();
        l.add_feature(Some(Geometry::point(2.3522, 48.8566)), &[])
            .unwrap();
        l
    }

    #[test]
    fn bins_points_into_cells_with_counts() {
        let id = wbvector::memory_store::put_vector(points_layer());
        let input = wbvector::memory_store::make_vector_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input,
            "resolution": "8",
        }))
        .unwrap();

        let ctx = ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        };
        let result = VectorToH3Tool.run(&args, &ctx).unwrap();
        assert_eq!(result.outputs["feature_count"], json!(3));
        // SF's two near-identical points collapse to one cell; Paris is another.
        assert_eq!(result.outputs["cell_count"], json!(2));
        assert_eq!(result.outputs["skipped"], json!(0));

        // The stored output layer should have 2 polygon features, counts summing to 3.
        let out_path = result.outputs["output"].as_str().unwrap();
        let out = load_input_layer(out_path).unwrap();
        assert_eq!(out.len(), 2);
        assert_eq!(out.geom_type, Some(GeometryType::Polygon));
        let total: i64 = out
            .iter()
            .map(|f| {
                let idx = out.schema.field_index("count").unwrap();
                f.attributes[idx].as_i64().unwrap()
            })
            .sum();
        assert_eq!(total, 3);
        // Each cell boundary is a hexagon (or pentagon): >= 5 vertices.
        for f in out.iter() {
            if let Some(Geometry::Polygon { exterior, .. }) = &f.geometry {
                assert!(
                    exterior.len() >= 5,
                    "H3 cell ring too small: {}",
                    exterior.len()
                );
            } else {
                panic!("expected polygon geometry");
            }
        }
    }

    #[test]
    fn rejects_out_of_range_resolution() {
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": "x", "resolution": "42" })).unwrap();
        assert!(VectorToH3Tool.validate(&args).is_err());
    }
}
