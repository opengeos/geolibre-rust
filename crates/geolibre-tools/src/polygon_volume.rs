//! GeoLibre tool: per-polygon volume and 3D surface area of a surface raster
//! above and/or below a per-feature reference plane.
//!
//! Pure-Rust counterpart of ArcGIS 3D Analyst's *Polygon Volume* — the zonal
//! counterpart of the repo's `surface_volume` (which integrates a whole raster
//! against one plane). For every polygon, using its own `height_field` as the
//! reference elevation, it accumulates over the cells whose centers fall inside:
//!   * `volume` — Σ |z − plane| · cellArea for the qualifying side(s);
//!   * `surface_area` — draped 3D area, each cell footprint scaled by the local
//!     finite-difference tilt `sqrt(1 + (dz/dx)² + (dz/dy)²)`;
//!   * `area_2d` — planimetric footprint of the qualifying cells.
//!
//! Neither `cut_fill` (which differences two surfaces), `storage_capacity`
//! (basin stage curves), nor the bundled `zonal_statistics` performs this
//! above/below-plane volumetric integration per polygon. No-data cells are
//! ignored.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::common::load_input_raster;
use crate::vector_common::{
    geometry_contains_point, load_input_layer, parse_optional_str, write_or_store_layer,
};

#[derive(Clone, Copy, PartialEq, Eq)]
enum Direction {
    Above,
    Below,
    Both,
}

pub struct PolygonVolumeTool;

impl Tool for PolygonVolumeTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "polygon_volume",
            display_name: "Polygon Volume",
            summary: "Per-polygon volume and 3D surface area of a surface raster above/below each feature's reference-plane height, like ArcGIS Polygon Volume — the zonal counterpart of surface_volume.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "surface",
                    description: "Input elevation surface raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "input",
                    description: "Polygon vector layer (zones).",
                    required: true,
                },
                ToolParamSpec {
                    name: "height_field",
                    description: "Attribute field giving each polygon's reference-plane elevation.",
                    required: true,
                },
                ToolParamSpec {
                    name: "direction",
                    description: "Which side of each plane to integrate: 'above' (default), 'below', or 'both'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read from the surface (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "surface")?;
        require_str(args, "input")?;
        require_str(args, "height_field")?;
        parse_direction(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let surface_path = require_str(args, "surface")?;
        let input_path = require_str(args, "input")?;
        let height_field = require_str(args, "height_field")?;
        let direction = parse_direction(args)?;
        let output = parse_optional_str(args, "output")?;
        let band_1 = parse_band(args)?;

        let surface = load_input_raster(surface_path)?;
        if (band_1 as usize) > surface.bands {
            return Err(ToolError::Validation(format!(
                "band {band_1} out of range (surface has {} band(s))",
                surface.bands
            )));
        }
        let band = (band_1 - 1) as isize;
        let layer = load_input_layer(input_path)?;
        let hidx = layer.schema.field_index(height_field).ok_or_else(|| {
            ToolError::Validation(format!("height_field '{height_field}' not found in input"))
        })?;

        let cell_x = surface.cell_size_x.abs();
        let cell_y = surface.cell_size_y.abs();
        let cell_area = cell_x * cell_y;
        if cell_area <= 0.0 {
            return Err(ToolError::Execution(
                "surface has a non-positive cell size".to_string(),
            ));
        }
        let nodata = surface.nodata;

        // Output layer: copy input schema, append the volume/area fields.
        let mut out = Layer::new("polygon_volume");
        if let Some(gt) = layer.geom_type {
            out = out.with_geom_type(gt);
        }
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for fd in layer.schema.fields() {
            out.add_field(fd.clone());
        }
        let new_fields: Vec<&str> = match direction {
            Direction::Above | Direction::Below => vec!["volume", "area_2d", "surface_area"],
            Direction::Both => vec![
                "volume_above",
                "volume_below",
                "area_above",
                "area_below",
                "surface_area_above",
                "surface_area_below",
            ],
        };
        for f in &new_fields {
            out.add_field(FieldDef::new(*f, FieldType::Float));
        }

        let mut total_volume = 0.0_f64;
        let n = layer.features.len();
        for (i, feat) in layer.features.iter().enumerate() {
            let Some(geom) = feat.geometry.as_ref() else {
                continue;
            };
            let plane = feat
                .attributes
                .get(hidx)
                .and_then(FieldValue::as_f64)
                .unwrap_or(f64::NAN);

            let (above, below) = if plane.is_finite() {
                integrate(
                    &surface, band, geom, plane, cell_x, cell_y, cell_area, nodata,
                )
            } else {
                (Totals::zero(), Totals::zero())
            };

            // Original attributes + new fields.
            let orig: Vec<(String, FieldValue)> = layer
                .schema
                .fields()
                .iter()
                .enumerate()
                .map(|(fi, fd)| {
                    (
                        fd.name.clone(),
                        feat.attributes.get(fi).cloned().unwrap_or(FieldValue::Null),
                    )
                })
                .collect();
            let mut attrs: Vec<(&str, FieldValue)> =
                orig.iter().map(|(k, v)| (k.as_str(), v.clone())).collect();
            match direction {
                Direction::Above => {
                    attrs.push(("volume", FieldValue::Float(above.volume)));
                    attrs.push(("area_2d", FieldValue::Float(above.area_2d)));
                    attrs.push(("surface_area", FieldValue::Float(above.area_3d)));
                    total_volume += above.volume;
                }
                Direction::Below => {
                    attrs.push(("volume", FieldValue::Float(below.volume)));
                    attrs.push(("area_2d", FieldValue::Float(below.area_2d)));
                    attrs.push(("surface_area", FieldValue::Float(below.area_3d)));
                    total_volume += below.volume;
                }
                Direction::Both => {
                    attrs.push(("volume_above", FieldValue::Float(above.volume)));
                    attrs.push(("volume_below", FieldValue::Float(below.volume)));
                    attrs.push(("area_above", FieldValue::Float(above.area_2d)));
                    attrs.push(("area_below", FieldValue::Float(below.area_2d)));
                    attrs.push(("surface_area_above", FieldValue::Float(above.area_3d)));
                    attrs.push(("surface_area_below", FieldValue::Float(below.area_3d)));
                    total_volume += above.volume + below.volume;
                }
            }
            out.add_feature(Some(geom.clone()), &attrs)
                .map_err(|e| ToolError::Execution(format!("failed adding feature: {e}")))?;
            ctx.progress.progress((i as f64 + 1.0) / n.max(1) as f64);
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(n));
        outputs.insert("total_volume".to_string(), json!(total_volume));
        Ok(ToolRunResult { outputs })
    }
}

struct Totals {
    area_2d: f64,
    area_3d: f64,
    volume: f64,
}

impl Totals {
    fn zero() -> Totals {
        Totals {
            area_2d: 0.0,
            area_3d: 0.0,
            volume: 0.0,
        }
    }
}

/// Integrates the surface within `geom` above and below `plane`.
#[allow(clippy::too_many_arguments)]
fn integrate(
    surface: &Raster,
    band: isize,
    geom: &Geometry,
    plane: f64,
    cell_x: f64,
    cell_y: f64,
    cell_area: f64,
    nodata: f64,
) -> (Totals, Totals) {
    let mut above = Totals::zero();
    let mut below = Totals::zero();
    let Some(bb) = geom.bbox() else {
        return (above, below);
    };
    let cols = surface.cols as isize;
    let rows = surface.rows as isize;
    // Cell window covering the polygon bbox.
    let c0 = (((bb.min_x - surface.x_min) / cell_x).floor() as isize).max(0);
    let c1 = (((bb.max_x - surface.x_min) / cell_x).ceil() as isize).min(cols);
    let top = ((rows as f64 - (bb.max_y - surface.y_min) / cell_y).floor() as isize).max(0);
    let bot = ((rows as f64 - (bb.min_y - surface.y_min) / cell_y).ceil() as isize).min(rows);

    for row in top..bot {
        for col in c0..c1 {
            let z = surface.get(band, row, col);
            if z == nodata || !z.is_finite() {
                continue;
            }
            let x = surface.x_min + (col as f64 + 0.5) * cell_x;
            let y = surface.y_min + (rows as f64 - 1.0 - row as f64 + 0.5) * cell_y;
            if !geometry_contains_point(geom, x, y) {
                continue;
            }
            let factor = surface_factor(surface, band, row, col, cell_x, cell_y, nodata);
            let cell_3d = cell_area * factor;
            if z >= plane {
                above.area_2d += cell_area;
                above.area_3d += cell_3d;
                above.volume += (z - plane) * cell_area;
            }
            if z <= plane {
                below.area_2d += cell_area;
                below.area_3d += cell_3d;
                below.volume += (plane - z) * cell_area;
            }
        }
    }
    (above, below)
}

/// Per-cell 3D-area scale factor from a central finite-difference gradient
/// (one-sided at edges; no-data neighbors treated as flat).
fn surface_factor(
    r: &Raster,
    band: isize,
    row: isize,
    col: isize,
    dx: f64,
    dy: f64,
    nodata: f64,
) -> f64 {
    let center = r.get(band, row, col);
    let sample = |rr: isize, cc: isize| -> Option<f64> {
        if rr < 0 || cc < 0 || rr >= r.rows as isize || cc >= r.cols as isize {
            return None;
        }
        let v = r.get(band, rr, cc);
        (v != nodata && v.is_finite()).then_some(v)
    };
    let grad = |a: Option<f64>, b: Option<f64>, step: f64| match (a, b) {
        (Some(a), Some(b)) => (b - a) / (2.0 * step),
        (Some(a), None) => (center - a) / step,
        (None, Some(b)) => (b - center) / step,
        (None, None) => 0.0,
    };
    let gx = grad(sample(row, col - 1), sample(row, col + 1), dx);
    let gy = grad(sample(row - 1, col), sample(row + 1, col), dy);
    (1.0 + gx * gx + gy * gy).sqrt()
}

fn parse_direction(args: &ToolArgs) -> Result<Direction, ToolError> {
    match args.get("direction").and_then(Value::as_str).map(str::trim) {
        None | Some("") | Some("above") => Ok(Direction::Above),
        Some("below") => Ok(Direction::Below),
        Some("both") => Ok(Direction::Both),
        Some(o) => Err(ToolError::Validation(format!(
            "'direction' must be 'above', 'below', or 'both', got '{o}'"
        ))),
    }
}

fn parse_band(args: &ToolArgs) -> Result<u64, ToolError> {
    match args.get("band") {
        None | Some(Value::Null) => Ok(1),
        Some(Value::Number(n)) => match n.as_f64() {
            Some(v) if v.fract() == 0.0 && v >= 1.0 => Ok(v as u64),
            _ => Err(ToolError::Validation(
                "'band' must be a positive integer".to_string(),
            )),
        },
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
            .map_err(|_| ToolError::Validation("'band' must be a positive integer".to_string()))
            .and_then(|v| {
                if v >= 1 {
                    Ok(v)
                } else {
                    Err(ToolError::Validation(
                        "'band' must be a positive integer".to_string(),
                    ))
                }
            }),
        Some(_) => Err(ToolError::Validation(
            "'band' must be a positive integer".to_string(),
        )),
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
    use wbraster::{DataType, RasterConfig};
    use wbvector::{Coord, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Flat surface at z, unit cells, cols x rows.
    fn surface(rows: usize, cols: usize, z: f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: Default::default(),
            metadata: Default::default(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, z).unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    /// One rectangular polygon covering [x0,x1]x[y0,y1] with a height attribute.
    fn poly_layer(x0: f64, y0: f64, x1: f64, y1: f64, height: f64) -> String {
        let mut l = Layer::new("p").with_geom_type(GeometryType::Polygon);
        l.add_field(FieldDef::new("h", FieldType::Float));
        let ring = vec![
            Coord::xy(x0, y0),
            Coord::xy(x1, y0),
            Coord::xy(x1, y1),
            Coord::xy(x0, y1),
            Coord::xy(x0, y0),
        ];
        l.add_feature(
            Some(Geometry::polygon(ring, vec![])),
            &[("h", FieldValue::Float(height))],
        )
        .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = PolygonVolumeTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    #[test]
    fn flat_surface_volume_matches_height_times_area() {
        // 4x4 flat surface at z=5; polygon over the central 2x2 (x in [1,3], y in [1,3]).
        let s = surface(4, 4, 5.0);
        let p = poly_layer(1.0, 1.0, 3.0, 3.0, 0.0);
        let (out, layer) =
            run(json!({ "surface": s, "input": p, "height_field": "h", "direction": "above" }));
        // 4 cells inside, each contributes (5-0)*1 = 5 -> volume 20, area 4.
        let f = &layer.features[0];
        let vol = f.get(&layer.schema, "volume").unwrap().as_f64().unwrap();
        let area = f.get(&layer.schema, "area_2d").unwrap().as_f64().unwrap();
        assert!((vol - 20.0).abs() < 1e-9, "vol {vol}");
        assert!((area - 4.0).abs() < 1e-9, "area {area}");
        assert!((out.outputs["total_volume"].as_f64().unwrap() - 20.0).abs() < 1e-9);
    }

    #[test]
    fn plane_above_surface_gives_below_volume() {
        let s = surface(4, 4, 5.0);
        let p = poly_layer(1.0, 1.0, 3.0, 3.0, 8.0);
        let (_o, layer) =
            run(json!({ "surface": s, "input": p, "height_field": "h", "direction": "below" }));
        let f = &layer.features[0];
        // plane 8 above surface 5 -> below volume = (8-5)*4 = 12.
        let vol = f.get(&layer.schema, "volume").unwrap().as_f64().unwrap();
        assert!((vol - 12.0).abs() < 1e-9, "vol {vol}");
    }

    #[test]
    fn both_directions_emit_split_fields() {
        let s = surface(4, 4, 5.0);
        let p = poly_layer(1.0, 1.0, 3.0, 3.0, 5.0);
        let (_o, layer) =
            run(json!({ "surface": s, "input": p, "height_field": "h", "direction": "both" }));
        let f = &layer.features[0];
        assert!(f.get(&layer.schema, "volume_above").is_ok());
        assert!(f.get(&layer.schema, "volume_below").is_ok());
        // plane == surface -> both volumes zero.
        assert!(
            f.get(&layer.schema, "volume_above")
                .unwrap()
                .as_f64()
                .unwrap()
                .abs()
                < 1e-9
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            PolygonVolumeTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "surface": "s.tif", "input": "p.shp" })).is_err(),
            "needs height_field"
        );
        assert!(bad(
            json!({ "surface": "s.tif", "input": "p.shp", "height_field": "h", "direction": "x" })
        )
        .is_err());
        assert!(bad(json!({ "surface": "s.tif", "input": "p.shp", "height_field": "h" })).is_ok());
    }
}
