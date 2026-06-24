//! GeoLibre tool: aggregate raster values into H3 cells (DGGS binning).
//!
//! The raster analogue of [`crate::vector_to_h3`]: each pixel is reduced to its
//! centre (lon/lat), assigned to the H3 cell containing it at a chosen
//! resolution, and its band value folded into that cell's aggregate (mean, sum,
//! min, max, count, or median). One polygon per occupied cell is emitted with
//! the aggregate `value` and the pixel `count`.
//!
//! H3 work is done by [`h3o`](https://docs.rs/h3o) (pure Rust, no C deps), so it
//! compiles to the same wasm targets as the rest of the suite. See issue #26.

use std::collections::BTreeMap;

use h3o::{LatLng, Resolution};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::common::load_input_raster;
use crate::vector_common::{parse_optional_str, write_or_store_layer};

/// Default H3 resolution (see [`crate::vector_to_h3`]); 8 is ~0.7 km² cells.
const DEFAULT_RESOLUTION: u8 = 8;

/// How to combine the pixel values that fall in a cell.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Aggregate {
    Mean,
    Sum,
    Min,
    Max,
    Count,
    Median,
}

impl Aggregate {
    fn parse(name: &str) -> Result<Self, ToolError> {
        match name.trim().to_ascii_lowercase().as_str() {
            "mean" => Ok(Self::Mean),
            "sum" => Ok(Self::Sum),
            "min" => Ok(Self::Min),
            "max" => Ok(Self::Max),
            "count" => Ok(Self::Count),
            "median" => Ok(Self::Median),
            other => Err(ToolError::Validation(format!(
                "unknown aggregate '{other}' (expected mean, sum, min, max, count, or median)"
            ))),
        }
    }
}

/// Running per-cell accumulator. `values` is only populated for `median`, which
/// needs the full sample; the other aggregates stay O(1) per cell.
struct Acc {
    count: u64,
    sum: f64,
    min: f64,
    max: f64,
    values: Vec<f64>,
}

impl Acc {
    fn new() -> Self {
        Self {
            count: 0,
            sum: 0.0,
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
            values: Vec::new(),
        }
    }

    fn push(&mut self, v: f64, keep_values: bool) {
        self.count += 1;
        self.sum += v;
        if v < self.min {
            self.min = v;
        }
        if v > self.max {
            self.max = v;
        }
        if keep_values {
            self.values.push(v);
        }
    }

    fn finalize(&mut self, agg: Aggregate) -> f64 {
        match agg {
            Aggregate::Mean => self.sum / self.count as f64,
            Aggregate::Sum => self.sum,
            Aggregate::Min => self.min,
            Aggregate::Max => self.max,
            Aggregate::Count => self.count as f64,
            Aggregate::Median => {
                self.values.sort_by(|a, b| a.total_cmp(b));
                let n = self.values.len();
                if n % 2 == 1 {
                    self.values[n / 2]
                } else {
                    (self.values[n / 2 - 1] + self.values[n / 2]) / 2.0
                }
            }
        }
    }
}

/// Bins raster pixel values into H3 cells.
pub struct RasterToH3Tool;

impl Tool for RasterToH3Tool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "raster_to_h3",
            display_name: "Raster to H3 Bins",
            summary: "Aggregate raster band values into H3 discrete-global-grid cells (mean, sum, min, max, count, median).",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input raster file path (must be geographic lon/lat, EPSG:4326).",
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
                ToolParamSpec {
                    name: "band",
                    description: "1-based raster band to aggregate (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "aggregate",
                    description: "How to combine pixel values per cell: mean, sum, min, max, count, or median (default mean).",
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
        if let Some(agg) = parse_optional_str(args, "aggregate")? {
            Aggregate::parse(agg)?;
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
        let aggregate = match parse_optional_str(args, "aggregate")? {
            Some(a) => Aggregate::parse(a)?,
            None => Aggregate::Mean,
        };
        let band_1based = args.get("band").and_then(Value::as_u64).unwrap_or(1).max(1);
        let band = (band_1based - 1) as isize;

        ctx.progress.info("reading input raster");
        let raster = load_input_raster(input)?;
        // H3 indexes geographic coordinates; refuse a non-4326 CRS rather than
        // silently producing garbage cells. A missing CRS is assumed to be 4326.
        if let Some(epsg) = raster.crs.epsg {
            if epsg != 4326 {
                return Err(ToolError::Validation(format!(
                    "input CRS is EPSG:{epsg}; reproject to EPSG:4326 (lon/lat) before H3 binning"
                )));
            }
        }
        if band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "band {band_1based} out of range (raster has {} band(s))",
                raster.bands
            )));
        }

        ctx.progress.info("aggregating pixels into H3 cells");
        let rows = raster.rows;
        let cols = raster.cols;
        let keep_values = aggregate == Aggregate::Median;
        let mut cells: BTreeMap<u64, Acc> = BTreeMap::new();
        let mut pixel_count = 0u64;
        let mut skipped = 0u64;
        for row in 0..rows {
            // Pixel-centre latitude (rows run top-down from the north edge).
            let lat = raster.y_min + (rows - row) as f64 * raster.cell_size_y
                - 0.5 * raster.cell_size_y;
            for col in 0..cols {
                let Some(v) = raster.get_opt(band, row as isize, col as isize) else {
                    skipped += 1; // nodata
                    continue;
                };
                let lng = raster.x_min + col as f64 * raster.cell_size_x + 0.5 * raster.cell_size_x;
                match LatLng::new(lat, lng) {
                    Ok(ll) => {
                        let cell = ll.to_cell(resolution);
                        cells
                            .entry(u64::from(cell))
                            .or_insert_with(Acc::new)
                            .push(v, keep_values);
                        pixel_count += 1;
                    }
                    Err(_) => skipped += 1, // out-of-range coordinate
                }
            }
        }

        ctx.progress.info("building H3 cell polygons");
        let mut out = Layer::new("h3_bins")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(4326);
        out.add_field(FieldDef::new("h3", FieldType::Text));
        out.add_field(FieldDef::new("value", FieldType::Float));
        out.add_field(FieldDef::new("count", FieldType::Integer));
        for (&cell_raw, acc) in cells.iter_mut() {
            let cell = h3o::CellIndex::try_from(cell_raw)
                .map_err(|e| ToolError::Execution(format!("invalid H3 cell: {e}")))?;
            let value = acc.finalize(aggregate);
            out.add_feature(
                Some(Geometry::polygon(cell_polygon_ring(cell), Vec::new())),
                &[
                    ("h3", FieldValue::Text(cell.to_string())),
                    ("value", FieldValue::Float(value)),
                    ("count", FieldValue::Integer(acc.count as i64)),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed building H3 feature: {e}")))?;
        }

        ctx.progress.info("writing output vector");
        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("resolution".to_string(), json!(res_u8));
        outputs.insert("cell_count".to_string(), json!(cells.len()));
        outputs.insert("pixel_count".to_string(), json!(pixel_count));
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

/// Builds the exterior ring (lon/lat coordinates) of an H3 cell's boundary,
/// without the closing duplicate vertex (matching `Ring`).
fn cell_polygon_ring(cell: h3o::CellIndex) -> Vec<Coord> {
    cell.boundary()
        .iter()
        .map(|ll| Coord::xy(ll.lng(), ll.lat()))
        .collect()
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

    /// A 2x2 lon/lat raster of ~11 m cells over San Francisco. Such a tiny
    /// extent sits well inside a single res-8 H3 cell (~461 m edge), so the
    /// three valid pixels deterministically share one cell. Top-down rows.
    /// Values: [[1, 2], [3, 100]] with 100 marked nodata.
    fn small_raster() -> Raster {
        let cfg = RasterConfig {
            cols: 2,
            rows: 2,
            bands: 1,
            x_min: -122.4194,
            y_min: 37.7749,
            cell_size: 0.0001,
            cell_size_y: Some(0.0001),
            nodata: 100.0,
            data_type: DataType::F64,
            crs: CrsInfo::from_epsg(4326),
            metadata: Vec::new(),
        };
        // data is band-major, row-major, top-down.
        Raster::from_data(cfg, vec![1.0, 2.0, 3.0, 100.0]).unwrap()
    }

    fn run_agg(agg: &str, res: &str) -> ToolRunResult {
        let id = wbraster::memory_store::put_raster(small_raster());
        let input = wbraster::memory_store::make_raster_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input, "resolution": res, "aggregate": agg,
        }))
        .unwrap();
        RasterToH3Tool.run(&args, &ctx()).unwrap()
    }

    #[test]
    fn skips_nodata_and_counts_pixels() {
        // At a coarse resolution all 3 valid pixels share one cell; nodata skipped.
        let r = run_agg("mean", "8");
        assert_eq!(r.outputs["pixel_count"], json!(3));
        assert_eq!(r.outputs["skipped"], json!(1));
    }

    #[test]
    fn aggregates_are_correct_when_pixels_share_a_cell() {
        use crate::vector_common::load_input_layer;
        // Coarse res: the 3 valid pixels (1, 2, 3) fall in one H3 cell.
        for (agg, expected) in [("mean", 2.0), ("sum", 6.0), ("min", 1.0), ("max", 3.0), ("count", 3.0), ("median", 2.0)] {
            let r = run_agg(agg, "8");
            assert_eq!(r.outputs["cell_count"], json!(1), "{agg}: expected a single shared cell");
            let out = load_input_layer(r.outputs["output"].as_str().unwrap()).unwrap();
            let vidx = out.schema.field_index("value").unwrap();
            let cidx = out.schema.field_index("count").unwrap();
            let f = out.iter().next().unwrap();
            assert_eq!(f.attributes[vidx].as_f64().unwrap(), expected, "aggregate {agg}");
            assert_eq!(f.attributes[cidx].as_i64().unwrap(), 3);
        }
    }

    #[test]
    fn rejects_non_4326_and_bad_args() {
        // Non-4326 CRS rejected.
        let mut r = small_raster();
        r.crs = CrsInfo::from_epsg(3857);
        let id = wbraster::memory_store::put_raster(r);
        let input = wbraster::memory_store::make_raster_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({ "input": input })).unwrap();
        assert!(RasterToH3Tool.run(&args, &ctx()).is_err());

        // Bad resolution / aggregate rejected at validate().
        let a: ToolArgs = serde_json::from_value(json!({ "input": "x", "resolution": "42" })).unwrap();
        assert!(RasterToH3Tool.validate(&a).is_err());
        let b: ToolArgs = serde_json::from_value(json!({ "input": "x", "aggregate": "nope" })).unwrap();
        assert!(RasterToH3Tool.validate(&b).is_err());
    }
}
