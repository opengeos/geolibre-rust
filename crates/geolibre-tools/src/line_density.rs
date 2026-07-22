//! GeoLibre tool: line density raster (length of lines per unit area).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Line Density* (Spatial Analyst). The
//! bundled `heat_map` is a point KDE and `edge_density` is a landscape metric
//! over a categorical raster — neither computes vector-line density, one of the
//! most common Spatial Analyst requests (road density, stream density, fault
//! density).
//!
//! For every output cell a circular neighborhood of `search_radius` is centered
//! on the cell. The tool sums the length of every line segment (times an optional
//! per-feature weight) that falls inside that circle, then divides by the
//! neighborhood area (`pi * r^2`). Segments are clipped to the circle with a
//! closed-form segment-circle intersection, so partial overlaps contribute their
//! exact in-circle length. The result is a density raster in `length / area`
//! units over the input's bounding box padded by the search radius.
//!
//! Work is scattered per segment: each segment only touches the cells inside its
//! bounding box expanded by the radius (a capsule region), which keeps the cost
//! proportional to the data rather than `cells * segments`.
//!
//! For a geographic (EPSG:4326) input the geometry is projected to a local
//! equirectangular metre frame centered on the extent so radii, lengths, and the
//! neighborhood area are true metres; the output raster's georeferencing is
//! converted back to degrees (with distinct x/y cell sizes) so it still overlays
//! the input. For a projected input everything is in the CRS's native units.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
use wbvector::Geometry;

use crate::common::{parse_optional_output, write_or_store_output};
use crate::vector_common::{load_input_layer, parse_optional_str};

/// Mean Earth radius (metres) for the local equirectangular projection.
const EARTH_R: f64 = 6_371_000.0;

pub struct LineDensityTool;

impl Tool for LineDensityTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "line_density",
            display_name: "Line Density",
            summary: "Density of linear features (length per unit area) in a circular neighborhood around each raster cell — road, stream, or fault density — by clipping each line segment to the search-radius circle (closed-form) and dividing the summed weighted length by pi*r^2. Like ArcGIS Line Density; the vector-line counterpart the bundled point-only heat_map and categorical edge_density lack.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polyline vector layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output density raster path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "weight_field",
                    description: "Numeric field weighting each line (e.g. lanes, traffic, population). Default: every line weight 1.",
                    required: false,
                },
                ToolParamSpec {
                    name: "search_radius",
                    description: "Neighborhood radius (metres for a geographic CRS, CRS units otherwise). Default: shorter extent side / 25.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Output cell size (same units as search_radius). Default: search_radius / 10.",
                    required: false,
                },
                ToolParamSpec {
                    name: "area_units",
                    description: "Area unit of the density values: 'square_map_units' (default), 'square_meters', 'square_kilometers', 'square_miles', or 'square_feet'.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
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

        // Local metre projection for a geographic CRS; identity otherwise.
        let geographic = layer.crs_epsg() == Some(4326);
        let epsg = layer.crs_epsg();

        // Collect line segments (in native coords) with their weight.
        struct RawSeg {
            ax: f64,
            ay: f64,
            bx: f64,
            by: f64,
            w: f64,
        }
        let mut raw: Vec<RawSeg> = Vec::new();
        let mut skipped = 0usize;
        for feat in layer.iter() {
            let Some(geom) = feat.geometry.as_ref() else {
                skipped += 1;
                continue;
            };
            let chains = line_chains(geom);
            if chains.is_empty() {
                skipped += 1;
                continue;
            }
            let w = match widx {
                Some(i) => feat
                    .attributes
                    .get(i)
                    .and_then(|v| v.as_f64())
                    .unwrap_or(1.0),
                None => 1.0,
            };
            if !(w.is_finite() && w > 0.0) {
                continue;
            }
            for chain in chains {
                for pair in chain.windows(2) {
                    raw.push(RawSeg {
                        ax: pair[0].0,
                        ay: pair[0].1,
                        bx: pair[1].0,
                        by: pair[1].1,
                        w,
                    });
                }
            }
        }
        if raw.is_empty() {
            return Err(ToolError::Execution(
                "input contains no line segments".to_string(),
            ));
        }

        // Native bounding box of the input lines.
        let (mut nxmin, mut nymin, mut nxmax, mut nymax) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for s in &raw {
            nxmin = nxmin.min(s.ax.min(s.bx));
            nxmax = nxmax.max(s.ax.max(s.bx));
            nymin = nymin.min(s.ay.min(s.by));
            nymax = nymax.max(s.ay.max(s.by));
        }

        // Projection factors: metres per native x/y unit at the extent center.
        let (kx, ky) = if geographic {
            let lat0 = 0.5 * (nymin + nymax);
            let ky = EARTH_R * std::f64::consts::PI / 180.0;
            let kx = ky * lat0.to_radians().cos().max(1e-9);
            (kx, ky)
        } else {
            (1.0, 1.0)
        };
        let (lon0, lat0) = (nxmin, nymin); // projection origin (native)
        let fwd_x = |x: f64| (x - lon0) * kx;
        let fwd_y = |y: f64| (y - lat0) * ky;

        // Working-frame (metre) segments.
        let segs: Vec<Seg> = raw
            .iter()
            .map(|s| Seg {
                ax: fwd_x(s.ax),
                ay: fwd_y(s.ay),
                bx: fwd_x(s.bx),
                by: fwd_y(s.by),
                w: s.w,
            })
            .collect();

        // Working-frame extent + defaults for radius / cell size.
        let ext_w = (nxmax - nxmin) * kx;
        let ext_h = (nymax - nymin) * ky;
        let radius = prm
            .search_radius
            .unwrap_or_else(|| (ext_w.min(ext_h) / 25.0).max(1e-6));
        let cell = prm.cell_size.unwrap_or_else(|| (radius / 10.0).max(1e-6));

        // Grid extent = line bbox padded by radius + one cell (so every circle
        // around a line point is fully covered → mass is conserved).
        let pad = radius + cell;
        let gxmin = -pad;
        let gymin = -pad;
        let gxmax = ext_w + pad;
        let gymax = ext_h + pad;
        let cols = (((gxmax - gxmin) / cell).ceil() as usize).max(1);
        let rows = (((gymax - gymin) / cell).ceil() as usize).max(1);
        let gymax = gymin + rows as f64 * cell; // snap top edge to whole cells

        let nbhd_area = std::f64::consts::PI * radius * radius;
        ctx.progress.info(&format!(
            "{} segment(s) -> {rows}x{cols} density raster (r={radius:.3}, cell={cell:.3})",
            segs.len()
        ));

        // Scatter each segment onto the cells within its expanded bounding box.
        let mut length = vec![0.0f64; rows * cols];
        for (si, s) in segs.iter().enumerate() {
            let sxmin = s.ax.min(s.bx) - radius;
            let sxmax = s.ax.max(s.bx) + radius;
            let symin = s.ay.min(s.by) - radius;
            let symax = s.ay.max(s.by) + radius;
            let c0 = (((sxmin - gxmin) / cell).floor() as isize).max(0) as usize;
            let c1 = (((sxmax - gxmin) / cell).ceil() as isize).min(cols as isize) as usize;
            let r0 = (((gymax - symax) / cell).floor() as isize).max(0) as usize;
            let r1 = (((gymax - symin) / cell).ceil() as isize).min(rows as isize) as usize;
            for r in r0..r1 {
                let cy = gymax - (r as f64 + 0.5) * cell;
                for c in c0..c1 {
                    let cx = gxmin + (c as f64 + 0.5) * cell;
                    let l = clipped_length(s.ax, s.ay, s.bx, s.by, cx, cy, radius);
                    if l > 0.0 {
                        length[r * cols + c] += l * s.w;
                    }
                }
            }
            if si % 256 == 0 {
                ctx.progress.progress((si as f64 + 1.0) / segs.len() as f64);
            }
        }

        // Density = weighted length / neighborhood area, scaled to area units.
        let unit_scale = prm.area_units.scale();
        let mut max_density = 0.0f64;
        let mut sum_len = 0.0f64;
        let data: Vec<f64> = length
            .iter()
            .map(|&l| {
                sum_len += l;
                let d = (l / nbhd_area) * unit_scale;
                if d > max_density {
                    max_density = d;
                }
                d
            })
            .collect();

        // Georeference the output back to native units.
        let (out_cell_x, out_cell_y) = if geographic {
            (cell / kx, cell / ky)
        } else {
            (cell, cell)
        };
        let out_xmin = lon0 + gxmin / kx;
        let out_ymin = lat0 + gymin / ky;

        let crs = CrsInfo {
            epsg,
            wkt: None,
            proj4: None,
        };
        let mut out = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: out_xmin,
            y_min: out_ymin,
            cell_size: out_cell_x,
            cell_size_y: Some(out_cell_y),
            nodata: -9999.0,
            data_type: DataType::F32,
            crs,
            metadata: Vec::new(),
        });
        for r in 0..rows {
            for c in 0..cols {
                out.set(0, r as isize, c as isize, data[r * cols + c])
                    .map_err(|e| ToolError::Execution(format!("write failed: {e}")))?;
            }
        }

        // Mass check: integral of density over the raster (in working units)
        // recovers the total weighted line length inside the extent.
        let integral_len = sum_len; // = Σ (clipped_len * w) already
        let out_path = write_or_store_output(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        outputs.insert("segment_count".to_string(), json!(segs.len()));
        outputs.insert("skipped".to_string(), json!(skipped));
        outputs.insert("search_radius".to_string(), json!(radius));
        outputs.insert("cell_size".to_string(), json!(cell));
        outputs.insert("neighborhood_area".to_string(), json!(nbhd_area));
        outputs.insert("max_density".to_string(), json!(max_density));
        // Σ clipped weighted length scattered across the grid (working units);
        // divided by (πr²)·cell² this equals Σ density·cell_area (mass check).
        outputs.insert("scattered_length".to_string(), json!(integral_len));
        Ok(ToolRunResult { outputs })
    }
}

/// A line segment in the working (metre / native) frame with its weight.
struct Seg {
    ax: f64,
    ay: f64,
    bx: f64,
    by: f64,
    w: f64,
}

/// Length of the portion of segment A→B that lies within distance `r` of C.
///
/// Closed-form segment-circle intersection: parametrize P(t)=A+t(B−A), solve the
/// quadratic |P(t)−C|²=r² for the entry/exit parameters, and intersect [t1,t2]
/// with [0,1]. The clipped length is (t_hi−t_lo)·|B−A|.
fn clipped_length(ax: f64, ay: f64, bx: f64, by: f64, cx: f64, cy: f64, r: f64) -> f64 {
    let dx = bx - ax;
    let dy = by - ay;
    let a = dx * dx + dy * dy;
    if a <= 0.0 {
        return 0.0; // zero-length segment contributes no length
    }
    let fx = ax - cx;
    let fy = ay - cy;
    let b = 2.0 * (fx * dx + fy * dy);
    let c = fx * fx + fy * fy - r * r;
    let disc = b * b - 4.0 * a * c;
    if disc <= 0.0 {
        return 0.0; // misses the circle (or grazes it: zero length)
    }
    let sq = disc.sqrt();
    let inv2a = 1.0 / (2.0 * a);
    let t1 = (-b - sq) * inv2a;
    let t2 = (-b + sq) * inv2a;
    let lo = t1.max(0.0);
    let hi = t2.min(1.0);
    if hi <= lo {
        return 0.0;
    }
    (hi - lo) * a.sqrt()
}

/// Extracts the point chains of a (multi)line geometry; non-line geometry yields
/// no chains (and is skipped by the caller).
fn line_chains(geom: &Geometry) -> Vec<Vec<(f64, f64)>> {
    let to_pts =
        |cs: &[wbvector::Coord]| -> Vec<(f64, f64)> { cs.iter().map(|c| (c.x, c.y)).collect() };
    match geom {
        Geometry::LineString(cs) => vec![to_pts(cs)],
        Geometry::MultiLineString(lines) => lines.iter().map(|l| to_pts(l)).collect(),
        _ => Vec::new(),
    }
}

// ── Area units ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum AreaUnits {
    MapUnits,
    Meters,
    Kilometers,
    Miles,
    Feet,
}

impl AreaUnits {
    /// Multiplier converting density in per-working-unit to per-target-unit.
    /// (density is 1/length, so per-km = per-metre × metres-per-km, etc.)
    fn scale(self) -> f64 {
        match self {
            AreaUnits::MapUnits | AreaUnits::Meters => 1.0,
            AreaUnits::Kilometers => 1000.0,
            AreaUnits::Miles => 1609.344,
            AreaUnits::Feet => 0.3048,
        }
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    weight_field: Option<String>,
    search_radius: Option<f64>,
    cell_size: Option<f64>,
    area_units: AreaUnits,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let weight_field = parse_optional_str(args, "weight_field")?.map(String::from);
    let search_radius = opt_pos(args, "search_radius")?;
    let cell_size = opt_pos(args, "cell_size")?;
    let area_units = match args
        .get("area_units")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some("square_map_units") | Some("map_units") => AreaUnits::MapUnits,
        Some("square_meters") | Some("meters") => AreaUnits::Meters,
        Some("square_kilometers") | Some("kilometers") => AreaUnits::Kilometers,
        Some("square_miles") | Some("miles") => AreaUnits::Miles,
        Some("square_feet") | Some("feet") => AreaUnits::Feet,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'area_units' must be one of square_map_units/square_meters/square_kilometers/square_miles/square_feet, got '{o}'"
            )))
        }
    };
    Ok(Params {
        weight_field,
        search_radius,
        cell_size,
        area_units,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn opt_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be a number"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a number"
        ))),
    }
}

fn opt_pos(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match opt_f64(args, key)? {
        Some(v) if v > 0.0 && v.is_finite() => Ok(Some(v)),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a positive number"
        ))),
        None => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Coord, FieldDef, FieldType, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Builds a projected line layer from polylines (each a list of (x,y)).
    fn line_layer(lines: &[(&[(f64, f64)], f64)]) -> String {
        let mut l = Layer::new("lines")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("w", FieldType::Float));
        for (pts, w) in lines {
            let coords: Vec<Coord> = pts.iter().map(|(x, y)| Coord::xy(*x, *y)).collect();
            l.add_feature(Some(Geometry::line_string(coords)), &[("w", (*w).into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Raster) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = LineDensityTool.run(&args, &ctx()).unwrap();
        let r = crate::common::load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, r)
    }

    /// Integrating density over the raster (Σ value · cell²) recovers the total
    /// line length — mass is conserved because the extent is padded by the radius.
    #[test]
    fn conserves_total_length() {
        // A single 1000-unit horizontal segment.
        let (out, r) = run(json!({
            "input": line_layer(&[(&[(0.0, 0.0), (1000.0, 0.0)], 1.0)]),
            "search_radius": 100.0, "cell_size": 10.0,
        }));
        let cell = out.outputs["cell_size"].as_f64().unwrap();
        let mut integral = 0.0;
        for row in 0..r.rows {
            for col in 0..r.cols {
                let v = r.get(0, row as isize, col as isize);
                if v != r.nodata {
                    integral += v * cell * cell;
                }
            }
        }
        assert!(
            (integral - 1000.0).abs() / 1000.0 < 0.02,
            "integrated density {integral} should recover length 1000"
        );
    }

    /// Peak density near a line is ~ (2·clip)/(π r²); a cell right on the line
    /// sees a full chord of length 2r, so density ≈ 2r/(π r²) = 2/(π r).
    #[test]
    fn peak_density_matches_formula() {
        let r_search = 50.0;
        let (out, r) = run(json!({
            "input": line_layer(&[(&[(-500.0, 0.0), (500.0, 0.0)], 1.0)]),
            "search_radius": r_search, "cell_size": 5.0,
        }));
        let expected = 2.0 / (std::f64::consts::PI * r_search);
        let got = out.outputs["max_density"].as_f64().unwrap();
        assert!(
            (got - expected).abs() / expected < 0.15,
            "peak density {got} should be near 2/(pi*r) = {expected}"
        );
        // A far-away cell must be exactly zero.
        assert_eq!(r.get(0, 0, 0), 0.0);
    }

    /// Doubling a line's weight doubles the density it produces.
    #[test]
    fn weight_scales_density() {
        let base = run(json!({
            "input": line_layer(&[(&[(0.0, 0.0), (400.0, 0.0)], 1.0)]),
            "search_radius": 60.0, "cell_size": 6.0,
        }))
        .0;
        let heavy = run(json!({
            "input": line_layer(&[(&[(0.0, 0.0), (400.0, 0.0)], 3.0)]),
            "weight_field": "w", "search_radius": 60.0, "cell_size": 6.0,
        }))
        .0;
        let b = base.outputs["max_density"].as_f64().unwrap();
        let h = heavy.outputs["max_density"].as_f64().unwrap();
        assert!((h / b - 3.0).abs() < 0.05, "weight 3 should triple density");
    }

    /// Non-line geometry is skipped, not counted.
    #[test]
    fn skips_non_line_geometry() {
        let mut l = Layer::new("mix")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("w", FieldType::Float));
        l.add_feature(
            Some(Geometry::line_string(vec![
                Coord::xy(0.0, 0.0),
                Coord::xy(200.0, 0.0),
            ])),
            &[("w", (1.0f64).into())],
        )
        .unwrap();
        l.add_feature(Some(Geometry::point(0.0, 0.0)), &[("w", (1.0f64).into())])
            .unwrap();
        let id = memory_store::put_vector(l);
        let path = memory_store::make_vector_memory_path(&id);
        let out = run(json!({ "input": path, "search_radius": 40.0, "cell_size": 4.0 })).0;
        assert_eq!(out.outputs["segment_count"], json!(1));
        assert_eq!(out.outputs["skipped"], json!(1));
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            LineDensityTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "search_radius": -1.0 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "area_units": "square_furlongs" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "area_units": "square_kilometers" })).is_ok());
    }
}
