//! GeoLibre tool: weighted Voronoi (dominance / market-area) allocation raster.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Generate Weighted Voronoi* (Spatial
//! Analyst). The bundled `voronoi_diagram` is unweighted — every site has equal
//! pull. Weighted Voronoi models unequal facilities (store market areas by size,
//! service dominance by capacity) and produces curved boundaries no exact-geometry
//! library in the stack computes; a raster allocation sidesteps that cleanly.
//!
//! Every cell is assigned to the site with the smallest *weighted* distance:
//! * **multiplicative** — `d / w` (Apollonius; larger weight = larger territory).
//! * **additive** — `d - w` (weight is a distance handicap/head-start).
//! * **power** — `d² - w²` (power/Laguerre diagram).
//!
//! The output is a categorical raster of 1-based site indices over the points'
//! extent (padded by `margin`). Polygonise it with `raster_to_vector_polygons`
//! for vector market areas.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
use wbvector::Geometry;

use crate::common::parse_optional_output;
use crate::vector_common::load_input_layer;

pub struct WeightedVoronoiTool;

impl Tool for WeightedVoronoiTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "weighted_voronoi",
            display_name: "Weighted Voronoi",
            summary: "Weighted Voronoi allocation raster (market/dominance areas): assign each cell to the site with the smallest weighted distance — multiplicative (d/w, Apollonius), additive (d-w), or power (d²-w²) — like ArcGIS Generate Weighted Voronoi. The unequal-site version of the bundled unweighted voronoi_diagram.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point layer of sites.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output allocation raster (1-based site index per cell). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "weight_field",
                    description: "Numeric field giving each site's weight/attractiveness (default: all sites weight 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "weight_type",
                    description: "'multiplicative' (d/w; default), 'additive' (d-w), or 'power' (d²-w²).",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Output cell size in CRS units (default: extent / 400).",
                    required: false,
                },
                ToolParamSpec {
                    name: "margin",
                    description: "Fraction of the site extent to pad the raster by (default 0.1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "epsg",
                    description: "EPSG to tag the output (default: from the input layer).",
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
        let output = parse_optional_output(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let widx =
            match &prm.weight_field {
                Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                    ToolError::Validation(format!("weight_field '{f}' not found"))
                })?),
                None => None,
            };

        // Collect sites.
        let mut sites: Vec<Site> = Vec::new();
        for feat in &layer.features {
            let Some((x, y)) = feat.geometry.as_ref().and_then(point_xy) else {
                continue;
            };
            let w = match widx {
                Some(i) => feat
                    .attributes
                    .get(i)
                    .and_then(|v| v.as_f64())
                    .unwrap_or(1.0),
                None => 1.0,
            };
            sites.push(Site {
                x,
                y,
                w: w.max(1e-9),
            });
        }
        if sites.is_empty() {
            return Err(ToolError::Execution("no point sites in input".to_string()));
        }

        // Extent from sites + margin.
        let (mut xmin, mut ymin, mut xmax, mut ymax) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for s in &sites {
            xmin = xmin.min(s.x);
            xmax = xmax.max(s.x);
            ymin = ymin.min(s.y);
            ymax = ymax.max(s.y);
        }
        let (mut dx, mut dy) = (xmax - xmin, ymax - ymin);
        if dx <= 0.0 {
            dx = 1.0;
        }
        if dy <= 0.0 {
            dy = 1.0;
        }
        xmin -= dx * prm.margin;
        xmax += dx * prm.margin;
        ymin -= dy * prm.margin;
        ymax += dy * prm.margin;
        let ext_w = xmax - xmin;
        let ext_h = ymax - ymin;
        let cell = prm
            .cell_size
            .unwrap_or((ext_w.max(ext_h) / 400.0).max(1e-9));
        let cols = ((ext_w / cell).ceil() as usize).max(1);
        let rows = ((ext_h / cell).ceil() as usize).max(1);

        ctx.progress.info(&format!(
            "{} site(s) -> {rows}x{cols} allocation raster ({})",
            sites.len(),
            prm.weight_type.label()
        ));

        let nodata = -1.0f64;
        let mut data = vec![nodata; rows * cols];
        for r in 0..rows {
            let cy = ymax - (r as f64 + 0.5) * cell;
            for c in 0..cols {
                let cx = xmin + (c as f64 + 0.5) * cell;
                let mut best = f64::INFINITY;
                let mut best_i = 0usize;
                for (i, s) in sites.iter().enumerate() {
                    let d = (cx - s.x).hypot(cy - s.y);
                    let wd = match prm.weight_type {
                        WeightType::Multiplicative => d / s.w,
                        WeightType::Additive => d - s.w,
                        WeightType::Power => d * d - s.w * s.w,
                    };
                    if wd < best {
                        best = wd;
                        best_i = i;
                    }
                }
                data[r * cols + c] = (best_i + 1) as f64;
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let out = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: xmin,
            y_min: ymin,
            cell_size: cell,
            cell_size_y: Some(cell),
            nodata,
            data_type: DataType::F32,
            crs: match prm.epsg.or_else(|| layer.crs_epsg()) {
                Some(e) => CrsInfo {
                    epsg: Some(e),
                    wkt: None,
                    proj4: None,
                },
                None => CrsInfo {
                    epsg: None,
                    wkt: None,
                    proj4: None,
                },
            },
            metadata: Vec::new(),
        });
        let mut out = out;
        for r in 0..rows {
            for c in 0..cols {
                out.set(0, r as isize, c as isize, data[r * cols + c])
                    .map_err(|e| ToolError::Execution(format!("write failed: {e}")))?;
            }
        }

        let out_path = crate::common::write_or_store_output(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("site_count".to_string(), json!(sites.len()));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        Ok(ToolRunResult { outputs })
    }
}

struct Site {
    x: f64,
    y: f64,
    w: f64,
}

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        _ => None,
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum WeightType {
    Multiplicative,
    Additive,
    Power,
}

impl WeightType {
    fn label(&self) -> &'static str {
        match self {
            WeightType::Multiplicative => "multiplicative",
            WeightType::Additive => "additive",
            WeightType::Power => "power",
        }
    }
}

struct Params {
    weight_field: Option<String>,
    weight_type: WeightType,
    cell_size: Option<f64>,
    margin: f64,
    epsg: Option<u32>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let weight_field =
        crate::vector_common::parse_optional_str(args, "weight_field")?.map(String::from);
    let weight_type = match args
        .get("weight_type")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some("multiplicative") => WeightType::Multiplicative,
        Some("additive") => WeightType::Additive,
        Some("power") => WeightType::Power,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'weight_type' must be multiplicative/additive/power, got '{o}'"
            )))
        }
    };
    let cell_size = match args.get("cell_size") {
        None | Some(Value::Null) => None,
        Some(Value::Number(n)) => n.as_f64().filter(|v| *v > 0.0),
        Some(Value::String(s)) if s.trim().is_empty() => None,
        Some(Value::String(s)) => Some(
            s.trim()
                .parse::<f64>()
                .map_err(|_| ToolError::Validation("'cell_size' must be a number".into()))?,
        ),
        _ => None,
    };
    let margin = match args.get("margin") {
        None | Some(Value::Null) => 0.1,
        Some(Value::Number(n)) => n.as_f64().unwrap_or(0.1).max(0.0),
        Some(Value::String(s)) if s.trim().is_empty() => 0.1,
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'margin' must be a number".into()))?
            .max(0.0),
        _ => 0.1,
    };
    let epsg = match args.get("epsg") {
        None | Some(Value::Null) => None,
        Some(Value::Number(n)) => n.as_u64().map(|v| v as u32),
        Some(Value::String(s)) if s.trim().is_empty() => None,
        Some(Value::String(s)) => Some(
            s.trim()
                .parse::<u32>()
                .map_err(|_| ToolError::Validation("'epsg' must be an integer".into()))?,
        ),
        _ => None,
    };
    Ok(Params {
        weight_field,
        weight_type,
        cell_size,
        margin,
        epsg,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn site_layer(pts: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("s")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("w", FieldType::Float));
        for (x, y, w) in pts {
            l.add_feature(Some(Geometry::point(*x, *y)), &[("w", (*w).into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = WeightedVoronoiTool.run(&args, &ctx()).unwrap();
        let r = crate::common::load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, r)
    }

    fn frac(r: &Raster, site: f64) -> f64 {
        let mut n = 0usize;
        let mut tot = 0usize;
        for row in 0..r.rows {
            for col in 0..r.cols {
                let v = r.get(0, row as isize, col as isize);
                if v != r.nodata {
                    tot += 1;
                    if v == site {
                        n += 1;
                    }
                }
            }
        }
        n as f64 / tot as f64
    }

    /// Equal weights -> the standard (unweighted) Voronoi: two symmetric sites
    /// split the field roughly 50/50.
    #[test]
    fn equal_weights_split_evenly() {
        let (_o, r) = run(json!({
            "input": site_layer(&[(0.0, 0.0, 1.0), (10.0, 0.0, 1.0)]),
            "weight_type": "multiplicative", "cell_size": 0.5, "margin": 0.2,
        }));
        let f1 = frac(&r, 1.0);
        assert!(
            (f1 - 0.5).abs() < 0.08,
            "equal weights should split ~50/50, got {f1}"
        );
    }

    /// A heavier multiplicative weight claims a larger territory.
    #[test]
    fn heavier_site_wins_more_area() {
        let (_o, r) = run(json!({
            "input": site_layer(&[(0.0, 0.0, 1.0), (10.0, 0.0, 3.0)]),
            "weight_field": "w", "weight_type": "multiplicative", "cell_size": 0.5, "margin": 0.2,
        }));
        let f2 = frac(&r, 2.0); // the weight-3 site
        assert!(f2 > 0.6, "the 3x-weight site should dominate, got {f2}");
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            WeightedVoronoiTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "weight_type": "log" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "weight_type": "power" })).is_ok());
    }
}
