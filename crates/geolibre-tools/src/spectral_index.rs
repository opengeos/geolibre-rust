//! GeoLibre tool: compute a named spectral index from a multi-band raster.
//!
//! The whitebox suite ships a few hardcoded derived products but no general
//! index calculator. This computes NDVI / NDWI / NDBI / NBR / EVI / SAVI from
//! band numbers the caller maps to red / nir / green / blue / swir.

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata,
    ToolParamSpec, ToolRunResult,
};
use wbraster::DataType;

use crate::common::{load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output};

/// Nodata value written for cells where any required band is nodata.
const INDEX_NODATA: f64 = -9999.0;

/// Which sensor band each index needs (1-based band numbers from the input).
#[derive(Clone, Copy)]
enum Index {
    Ndvi,
    Ndwi,
    Ndbi,
    Nbr,
    Evi,
    Savi,
}

impl Index {
    fn parse(name: &str) -> Result<Self, ToolError> {
        match name.trim().to_ascii_lowercase().as_str() {
            "ndvi" => Ok(Self::Ndvi),
            "ndwi" => Ok(Self::Ndwi),
            "ndbi" => Ok(Self::Ndbi),
            "nbr" => Ok(Self::Nbr),
            "evi" => Ok(Self::Evi),
            "savi" => Ok(Self::Savi),
            other => Err(ToolError::Validation(format!(
                "unknown index '{other}' (expected ndvi, ndwi, ndbi, nbr, evi, or savi)"
            ))),
        }
    }
}

/// Computes a normalized or enhanced vegetation/water/built-up index.
pub struct SpectralIndexTool;

impl Tool for SpectralIndexTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "spectral_index",
            display_name: "Spectral Index",
            summary: "Compute a spectral index (NDVI, NDWI, NDBI, NBR, EVI, SAVI) from a multi-band raster.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input multi-band raster file path.",
                    required: true,
                },
                ToolParamSpec {
                    name: "index",
                    description: "Index to compute: ndvi, ndwi, ndbi, nbr, evi, or savi.",
                    required: true,
                },
                ToolParamSpec { name: "red", description: "1-based band number for the red band.", required: false },
                ToolParamSpec { name: "nir", description: "1-based band number for the near-infrared band.", required: false },
                ToolParamSpec { name: "green", description: "1-based band number for the green band.", required: false },
                ToolParamSpec { name: "blue", description: "1-based band number for the blue band.", required: false },
                ToolParamSpec { name: "swir", description: "1-based band number for the shortwave-infrared band.", required: false },
                ToolParamSpec {
                    name: "soil_factor",
                    description: "Soil-adjustment factor L for SAVI (default 0.5).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output raster path. If omitted, the result is stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args.get("input").and_then(Value::as_str).map(str::trim).unwrap_or("").is_empty() {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        match args.get("index").and_then(Value::as_str) {
            Some(name) => Index::parse(name).map(|_| ()),
            None => Err(ToolError::Validation(
                "missing required parameter 'index'".to_string(),
            )),
        }
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args
            .get("input")
            .and_then(Value::as_str)
            .ok_or_else(|| ToolError::Validation("missing required parameter 'input'".to_string()))?;
        let index = Index::parse(
            args.get("index")
                .and_then(Value::as_str)
                .ok_or_else(|| ToolError::Validation("missing required parameter 'index'".to_string()))?,
        )?;
        let output = parse_optional_output(args, "output")?;
        let soil_factor = args.get("soil_factor").and_then(Value::as_f64).unwrap_or(0.5);

        let raster = load_input_raster(input)?;
        // Resolve the 0-based band indices each formula needs.
        let band = |key: &str| -> Result<isize, ToolError> {
            let n = args.get(key).and_then(Value::as_u64).ok_or_else(|| {
                ToolError::Validation(format!(
                    "index '{}' requires the '{key}' band number (1-based)",
                    args.get("index").and_then(Value::as_str).unwrap_or("?")
                ))
            })?;
            if n == 0 {
                return Err(ToolError::Validation(format!(
                    "'{key}' band number is 1-based and must be >= 1"
                )));
            }
            let idx = (n - 1) as isize;
            if idx as usize >= raster.bands {
                return Err(ToolError::Validation(format!(
                    "band {n} out of range (raster has {} band(s))",
                    raster.bands
                )));
            }
            Ok(idx)
        };

        let (red, nir, green, blue, swir) = match index {
            Index::Ndvi => (Some(band("red")?), Some(band("nir")?), None, None, None),
            Index::Savi => (Some(band("red")?), Some(band("nir")?), None, None, None),
            Index::Evi => (Some(band("red")?), Some(band("nir")?), None, Some(band("blue")?), None),
            Index::Ndwi => (None, Some(band("nir")?), Some(band("green")?), None, None),
            Index::Ndbi => (None, Some(band("nir")?), None, None, Some(band("swir")?)),
            Index::Nbr => (None, Some(band("nir")?), None, None, Some(band("swir")?)),
        };

        let nodata = raster.nodata;
        let rows = raster.rows as isize;
        let cols = raster.cols as isize;
        let valid = |b: Option<isize>, row: isize, col: isize| -> Option<f64> {
            let b = b?;
            let v = raster.get(b, row, col);
            if v == nodata || !v.is_finite() {
                None
            } else {
                Some(v)
            }
        };

        ctx.progress.info("computing spectral index");
        let mut data = vec![INDEX_NODATA; (rows * cols) as usize];
        for row in 0..rows {
            for col in 0..cols {
                let r = valid(red, row, col);
                let n = valid(nir, row, col);
                let g = valid(green, row, col);
                let bl = valid(blue, row, col);
                let sw = valid(swir, row, col);
                let value = match index {
                    Index::Ndvi => norm_diff(n, r),
                    Index::Ndwi => norm_diff(g, n),
                    Index::Ndbi => norm_diff(sw, n),
                    Index::Nbr => norm_diff(n, sw),
                    Index::Evi => match (n, r, bl) {
                        (Some(n), Some(r), Some(b)) => {
                            let denom = n + 6.0 * r - 7.5 * b + 1.0;
                            if denom == 0.0 { None } else { Some(2.5 * (n - r) / denom) }
                        }
                        _ => None,
                    },
                    Index::Savi => match (n, r) {
                        (Some(n), Some(r)) => {
                            let denom = n + r + soil_factor;
                            if denom == 0.0 { None } else { Some((1.0 + soil_factor) * (n - r) / denom) }
                        }
                        _ => None,
                    },
                };
                if let Some(v) = value.filter(|v| v.is_finite()) {
                    data[(row * cols + col) as usize] = v;
                }
            }
            ctx.progress.progress((row as f64 + 1.0) / rows as f64);
        }

        let out = raster_like_with_data(&raster, data, INDEX_NODATA, DataType::F32)?;
        let out_path = write_or_store_output(out, output)?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        Ok(ToolRunResult { outputs })
    }
}

/// Normalized difference `(a - b) / (a + b)`; `None` if either is missing or the
/// sum is zero.
fn norm_diff(a: Option<f64>, b: Option<f64>) -> Option<f64> {
    match (a, b) {
        (Some(a), Some(b)) => {
            let sum = a + b;
            if sum == 0.0 {
                None
            } else {
                Some((a - b) / sum)
            }
        }
        _ => None,
    }
}
