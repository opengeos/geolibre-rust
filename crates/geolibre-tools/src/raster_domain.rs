//! GeoLibre tool: the footprint of a raster's valid data.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Raster Domain* (3D Analyst).
//!
//! ## Why the catalog needs it
//!
//! A satellite scene is a rotated parallelogram inside a north-up rectangle; a
//! clipped DEM is a watershed inside a bounding box; a mosaic tile has ragged
//! no-data margins where the swaths did not reach. In every case the useful
//! extent is the *valid-data* region, and it is what you need in order to build
//! a seamline, index a tile set, clip a neighbour, or draw an honest coverage
//! map. The bounding rectangle is not that region and can be several times its
//! area.
//!
//! The bundled `layer_footprint_raster` returns exactly that rectangle — its
//! own summary says "a polygon footprint representing the full **extent**" —
//! so it cannot answer the question. `polygonize` could, but only after the
//! caller has separately built a 0/1 validity mask, which is the whole job.
//!
//! ## Method
//!
//! Cells are classified valid or not — no-data, non-finite, and optionally
//! outside a `value_range` or equal to a nominated `ignore_value` — and the
//! valid set is traced with the shared `polygonize` machinery. Interior holes
//! (a no-data lake inside a scene) are preserved as polygon holes, so the
//! footprint really is the domain and not just its outline. `geometry_type:
//! line` emits the same boundaries as closed polylines instead.

use std::collections::BTreeMap;
use std::collections::HashMap;

use serde_json::{json, Map, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::args_common::{band_index, bool_or, choice_or, opt_f64, opt_positive_f64, req_str};
use crate::common::{load_input_raster, parse_optional_output};
use crate::geojson_geom::geometry_from_json;
use crate::polygonize::{polygonize_to_geojson, PolygonizeParams};
use crate::vector_common::write_or_store_layer;

pub struct RasterDomainTool;

impl Tool for RasterDomainTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "raster_domain",
            display_name: "Raster Domain",
            summary: "Traces the boundary of a raster's valid-data region as polygons (or closed polylines), preserving interior no-data holes (ArcGIS Raster Domain). The bundled layer_footprint_raster returns the full rectangular extent, which for a rotated scene, a clipped watershed or a ragged mosaic tile can be several times the actual coverage; polygonize could do it but only from a 0/1 validity mask the caller has to build first, which is the entire task.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input raster whose valid-data footprint is wanted.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon or line layer. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "geometry_type",
                    description: "'polygon' (default) or 'line' for the boundary as closed polylines.",
                    required: false,
                },
                ToolParamSpec {
                    name: "ignore_value",
                    description: "Treat this value as invalid too, alongside the raster's own no-data (commonly 0 for imagery with a black border).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_value",
                    description: "Treat values below this as invalid.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_value",
                    description: "Treat values above this as invalid.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_area",
                    description: "Discard footprint parts smaller than this area in map units squared, which removes speckle around a ragged margin.",
                    required: false,
                },
                ToolParamSpec {
                    name: "fill_holes",
                    description: "Drop interior holes and return the solid footprint (default false).",
                    required: false,
                },
                ToolParamSpec {
                    name: "all_bands",
                    description: "Require a cell to be valid in every band rather than just the selected one (default false).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to test (default 1). Ignored when 'all_bands' is set.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input_path = req_str(args, "input")?.to_string();
        let prm = parse_params(args)?;
        let band = band_index(args, "band")?;
        let output = parse_optional_output(args, "output")?;

        let raster = load_input_raster(&input_path)?;
        let (rows, cols) = (raster.rows, raster.cols);

        let valid = |r: usize, c: usize, b: isize| -> bool {
            let v = raster.get(b, r as isize, c as isize);
            if v == raster.nodata || !v.is_finite() {
                return false;
            }
            if let Some(ig) = prm.ignore_value {
                if v == ig {
                    return false;
                }
            }
            if let Some(mn) = prm.min_value {
                if v < mn {
                    return false;
                }
            }
            if let Some(mx) = prm.max_value {
                if v > mx {
                    return false;
                }
            }
            true
        };

        let mut labels = vec![0.0f64; rows * cols];
        let mut valid_cells = 0usize;
        for r in 0..rows {
            for c in 0..cols {
                let ok = if prm.all_bands {
                    (0..raster.bands).all(|b| valid(r, c, b as isize))
                } else {
                    valid(r, c, band)
                };
                if ok {
                    labels[r * cols + c] = 1.0;
                    valid_cells += 1;
                }
            }
        }
        if valid_cells == 0 {
            return Err(ToolError::Execution(
                "the raster has no valid cells, so it has no domain".to_string(),
            ));
        }

        let cell_area = raster.cell_size_x * raster.cell_size_y;
        ctx.progress.info(&format!(
            "{rows}x{cols}, {valid_cells} valid cell(s) = {:.0} of {:.0} map units squared",
            valid_cells as f64 * cell_area,
            (rows * cols) as f64 * cell_area
        ));

        let props: HashMap<i64, Map<String, Value>> = HashMap::new();
        let geojson = polygonize_to_geojson(&PolygonizeParams {
            labels: &labels,
            rows,
            cols,
            x_min: raster.x_min,
            y_max: raster.y_min + rows as f64 * raster.cell_size_y,
            cell_size_x: raster.cell_size_x,
            cell_size_y: raster.cell_size_y,
            epsg: raster.crs.epsg,
            props_by_id: &props,
        });
        let parsed: Value = serde_json::from_str(&geojson).map_err(|e| {
            ToolError::Execution(format!("polygonize produced invalid GeoJSON: {e}"))
        })?;
        let feats = parsed
            .get("features")
            .and_then(Value::as_array)
            .map(Vec::as_slice)
            .unwrap_or_default();

        let mut layer = Layer::new("raster_domain");
        layer.geom_type = Some(if prm.as_lines {
            GeometryType::LineString
        } else {
            GeometryType::Polygon
        });
        if let Some(e) = raster.crs.epsg {
            layer = layer.with_crs_epsg(e);
        }
        layer.add_field(FieldDef::new("id", FieldType::Integer));
        layer.add_field(FieldDef::new("area", FieldType::Float));
        layer.add_field(FieldDef::new("hole_count", FieldType::Integer));
        layer.add_field(FieldDef::new("ring", FieldType::Text));

        let mut fid = 0u64;
        let mut total_area = 0.0;
        let mut parts = 0usize;
        for f in feats {
            let Some(geom) = f.get("geometry").and_then(geometry_from_json) else {
                continue;
            };
            let Geometry::Polygon {
                exterior,
                interiors,
            } = geom
            else {
                continue;
            };
            // Net area: the exterior less its holes, so a scene with a lake
            // does not report the lake as coverage.
            let outer = ring_area(exterior.coords());
            let holes: f64 = interiors.iter().map(|r| ring_area(r.coords())).sum();
            let net = if prm.fill_holes { outer } else { outer - holes };
            if let Some(min) = prm.min_area {
                if net < min {
                    continue;
                }
            }
            parts += 1;
            total_area += net;

            if prm.as_lines {
                // One closed polyline per ring, exterior first.
                let mut emit = |coords: &[Coord], kind: &str, fid: &mut u64| {
                    let mut cs = coords.to_vec();
                    if let Some(first) = cs.first().cloned() {
                        cs.push(first);
                    }
                    let mut feat = Feature::with_geometry(
                        *fid,
                        Geometry::LineString(cs),
                        layer.schema.len(),
                    );
                    feat.set_by_index(0, FieldValue::Integer(*fid as i64));
                    feat.set_by_index(1, FieldValue::Float(net));
                    // Match the polygon branch: with `fill_holes` the holes
                    // are not emitted, so the count must read 0 in both forms.
                    feat.set_by_index(
                        2,
                        FieldValue::Integer(if prm.fill_holes {
                            0
                        } else {
                            interiors.len() as i64
                        }),
                    );
                    feat.set_by_index(3, FieldValue::Text(kind.to_string()));
                    layer.push(feat);
                    *fid += 1;
                };
                emit(exterior.coords(), "exterior", &mut fid);
                if !prm.fill_holes {
                    for r in &interiors {
                        emit(r.coords(), "hole", &mut fid);
                    }
                }
            } else {
                let geom = Geometry::Polygon {
                    exterior,
                    interiors: if prm.fill_holes {
                        Vec::new()
                    } else {
                        interiors.clone()
                    },
                };
                let mut feat = Feature::with_geometry(fid, geom, layer.schema.len());
                feat.set_by_index(0, FieldValue::Integer(fid as i64));
                feat.set_by_index(1, FieldValue::Float(net));
                feat.set_by_index(
                    2,
                    FieldValue::Integer(if prm.fill_holes {
                        0
                    } else {
                        interiors.len() as i64
                    }),
                );
                feat.set_by_index(3, FieldValue::Text("exterior".to_string()));
                layer.push(feat);
                fid += 1;
            }
        }

        let feature_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("part_count".to_string(), json!(parts));
        outputs.insert("valid_cells".to_string(), json!(valid_cells));
        outputs.insert("domain_area".to_string(), json!(total_area));
        outputs.insert(
            "extent_area".to_string(),
            json!((rows * cols) as f64 * cell_area),
        );
        Ok(ToolRunResult { outputs })
    }
}

/// Absolute shoelace area of an unclosed ring.
fn ring_area(coords: &[Coord]) -> f64 {
    let n = coords.len();
    if n < 3 {
        return 0.0;
    }
    let mut a = 0.0;
    for i in 0..n {
        let p = &coords[i];
        let q = &coords[(i + 1) % n];
        a += p.x * q.y - q.x * p.y;
    }
    (a / 2.0).abs()
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    as_lines: bool,
    ignore_value: Option<f64>,
    min_value: Option<f64>,
    max_value: Option<f64>,
    min_area: Option<f64>,
    fill_holes: bool,
    all_bands: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let as_lines = choice_or(args, "geometry_type", &["polygon", "line"], "polygon")? == "line";
    let ignore_value = opt_f64(args, "ignore_value")?;
    let min_value = opt_f64(args, "min_value")?;
    let max_value = opt_f64(args, "max_value")?;
    if let (Some(lo), Some(hi)) = (min_value, max_value) {
        if lo > hi {
            return Err(ToolError::Validation(format!(
                "'min_value' ({lo}) exceeds 'max_value' ({hi}); no cell could be valid"
            )));
        }
    }
    let min_area = opt_positive_f64(args, "min_area")?;
    let fill_holes = bool_or(args, "fill_holes", false)?;
    let all_bands = bool_or(args, "all_bands", false)?;
    Ok(Params {
        as_lines,
        ignore_value,
        min_value,
        max_value,
        min_area,
        fill_holes,
        all_bands,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vector_common::load_input_layer;
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

    fn raster_of(cols: usize, rows: usize, bands: usize, vals: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands,
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
        for b in 0..bands {
            for row in 0..rows {
                for col in 0..cols {
                    r.set(
                        b as isize,
                        row as isize,
                        col as isize,
                        vals[b * rows * cols + row * cols + col],
                    )
                    .unwrap();
                }
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, BTreeMap<String, Value>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = RasterDomainTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (layer, out.outputs)
    }

    /// The point of the tool: for a diagonal swath inside a north-up raster the
    /// domain is far smaller than the rectangular extent that
    /// `layer_footprint_raster` would return.
    #[test]
    fn footprint_is_smaller_than_the_extent() {
        let (rows, cols) = (10, 10);
        // A diagonal band of valid cells; everything else is no-data.
        let mut v = vec![-9999.0; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                if (c as i32 - r as i32).abs() <= 1 {
                    v[r * cols + c] = 1.0;
                }
            }
        }
        let (_, outputs) = run(json!({ "input": raster_of(cols, rows, 1, &v) }));
        let domain = outputs["domain_area"].as_f64().unwrap();
        let extent = outputs["extent_area"].as_f64().unwrap();
        assert!((extent - 10_000.0).abs() < 1e-6, "extent {extent}");
        assert!(
            domain < 0.4 * extent,
            "the diagonal swath domain {domain} should be far under the extent {extent}"
        );
        // 28 valid cells of 100 m^2.
        assert!((domain - 2800.0).abs() < 1e-6, "domain area {domain}");
    }

    /// An interior no-data lake becomes a polygon hole, so the reported area is
    /// the real coverage rather than the outline.
    #[test]
    fn interior_nodata_becomes_a_hole() {
        let (rows, cols) = (8, 8);
        let mut v = vec![5.0; rows * cols];
        for r in 3..5 {
            for c in 3..5 {
                v[r * cols + c] = -9999.0;
            }
        }
        let (layer, outputs) = run(json!({ "input": raster_of(cols, rows, 1, &v) }));
        assert_eq!(layer.len(), 1);
        let f = layer.iter().next().unwrap();
        let hi = layer.schema.field_index("hole_count").unwrap();
        assert!(
            matches!(f.attributes[hi], FieldValue::Integer(1)),
            "the lake should be one hole"
        );
        // 64 cells less the 4-cell lake, at 100 m^2 each.
        assert!((outputs["domain_area"].as_f64().unwrap() - 6000.0).abs() < 1e-6);

        // fill_holes returns the solid outline instead.
        let (filled, filled_out) = run(json!({
            "input": raster_of(cols, rows, 1, &v), "fill_holes": true
        }));
        let ff = filled.iter().next().unwrap();
        assert!(matches!(ff.attributes[hi], FieldValue::Integer(0)));
        assert!((filled_out["domain_area"].as_f64().unwrap() - 6400.0).abs() < 1e-6);
    }

    /// A fully valid raster gives the whole rectangle — the case that used to
    /// panic inside `polygonize` before the row-wrap fix.
    #[test]
    fn fully_valid_raster_traces_the_whole_rectangle() {
        let (rows, cols) = (5, 5);
        let (layer, outputs) = run(json!({ "input": raster_of(cols, rows, 1, &[3.0; 25]) }));
        assert_eq!(layer.len(), 1);
        assert!((outputs["domain_area"].as_f64().unwrap() - 2500.0).abs() < 1e-6);
        assert!(
            (outputs["domain_area"].as_f64().unwrap()
                - outputs["extent_area"].as_f64().unwrap())
            .abs()
                < 1e-6,
            "a fully valid raster's domain is its extent"
        );
    }

    /// `ignore_value` handles the black border that imagery usually carries
    /// instead of a proper no-data flag.
    #[test]
    fn ignore_value_strips_a_zero_border() {
        let (rows, cols) = (6, 6);
        let mut v = vec![0.0; rows * cols];
        for r in 1..5 {
            for c in 1..5 {
                v[r * cols + c] = 42.0;
            }
        }
        let (_, plain) = run(json!({ "input": raster_of(cols, rows, 1, &v) }));
        assert!((plain["domain_area"].as_f64().unwrap() - 3600.0).abs() < 1e-6);

        let (_, stripped) = run(json!({
            "input": raster_of(cols, rows, 1, &v), "ignore_value": 0.0
        }));
        // 16 interior cells of 100 m^2.
        assert!((stripped["domain_area"].as_f64().unwrap() - 1600.0).abs() < 1e-6);
    }

    /// A value range narrows what counts as valid.
    #[test]
    fn value_range_restricts_the_domain() {
        let (rows, cols) = (4, 4);
        let v: Vec<f64> = (0..16).map(|i| i as f64).collect();
        let (_, outputs) = run(json!({
            "input": raster_of(cols, rows, 1, &v), "min_value": 8.0
        }));
        // Values 8..15 = 8 cells.
        assert_eq!(outputs["valid_cells"].as_u64().unwrap(), 8);
    }

    /// Disconnected parts are separate features, and `min_area` drops specks.
    #[test]
    fn min_area_drops_speckle() {
        let (rows, cols) = (10, 10);
        let mut v = vec![-9999.0; rows * cols];
        for r in 1..5 {
            for c in 1..5 {
                v[r * cols + c] = 1.0; // 16 cells = 1600 m^2
            }
        }
        v[9 * cols + 9] = 1.0; // a single speck = 100 m^2
        let (all, _) = run(json!({ "input": raster_of(cols, rows, 1, &v) }));
        assert_eq!(all.len(), 2);
        let (big, outputs) = run(json!({
            "input": raster_of(cols, rows, 1, &v), "min_area": 500.0
        }));
        assert_eq!(big.len(), 1);
        assert_eq!(outputs["part_count"].as_u64().unwrap(), 1);
    }

    /// `all_bands` requires validity everywhere, not just in one band.
    #[test]
    fn all_bands_intersects_the_valid_masks() {
        let (rows, cols) = (4, 4);
        // Band 0 valid on the left half, band 1 valid on the top half.
        let mut b0 = vec![-9999.0; 16];
        let mut b1 = vec![-9999.0; 16];
        for r in 0..4 {
            for c in 0..4 {
                if c < 2 {
                    b0[r * cols + c] = 1.0;
                }
                if r < 2 {
                    b1[r * cols + c] = 1.0;
                }
            }
        }
        let mut both = b0.clone();
        both.extend(b1);
        let src = || raster_of(cols, rows, 2, &both);

        let (_, single) = run(json!({ "input": src(), "band": 1 }));
        assert_eq!(single["valid_cells"].as_u64().unwrap(), 8);
        let (_, all) = run(json!({ "input": src(), "all_bands": true }));
        assert_eq!(
            all["valid_cells"].as_u64().unwrap(),
            4,
            "only the top-left quadrant is valid in both bands"
        );
    }

    /// The line form emits closed boundaries, one per ring.
    #[test]
    fn line_form_emits_closed_rings() {
        let (rows, cols) = (8, 8);
        let mut v = vec![5.0; rows * cols];
        for r in 3..5 {
            for c in 3..5 {
                v[r * cols + c] = -9999.0;
            }
        }
        let (layer, _) = run(json!({
            "input": raster_of(cols, rows, 1, &v), "geometry_type": "line"
        }));
        assert_eq!(layer.len(), 2, "one exterior plus one hole boundary");
        for f in layer.iter() {
            let Some(Geometry::LineString(cs)) = f.geometry.as_ref() else {
                panic!("expected a line")
            };
            assert_eq!(
                (cs.first().map(|c| (c.x, c.y)), cs.last().map(|c| (c.x, c.y))).0,
                (cs.first().map(|c| (c.x, c.y)), cs.last().map(|c| (c.x, c.y))).1,
                "boundary polylines must close"
            );
        }
    }

    /// A raster with nothing valid has no domain; say so.
    #[test]
    fn empty_raster_is_an_error() {
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": raster_of(3, 3, 1, &[-9999.0; 9]) }))
                .unwrap();
        let err = RasterDomainTool.run(&args, &ctx()).unwrap_err();
        assert!(
            format!("{err:?}").contains("no valid cells"),
            "expected an empty-raster error, got {err:?}"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            RasterDomainTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({"input": "a.tif", "geometry_type": "point"})).is_err());
        assert!(bad(json!({"input": "a.tif", "min_value": 10, "max_value": 5})).is_err());
        assert!(bad(json!({"input": "a.tif", "min_area": -1})).is_err());
        assert!(bad(json!({"input": "a.tif", "geometry_type": "line"})).is_ok());
    }
}
