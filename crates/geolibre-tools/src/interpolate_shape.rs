//! GeoLibre tool: drape vector features on a surface and add surface metrics.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Interpolate Shape* + *Add Surface
//! Information* (3D Analyst). Nothing bundled drapes arbitrary vectors on a
//! raster — point-sampling tools handle points only, and none compute surface
//! length or slope along a line. This is the bridge between the repo's raster
//! (terrain) and vector halves: trail profiles, pipeline lengths, slope-aware
//! routing prep.
//!
//! Each feature's geometry is densified at `sample_distance` (so long segments
//! follow the terrain), every vertex is sampled from the surface (bilinear or
//! nearest) and given a Z value, and per-feature surface metrics are written as
//! attributes:
//!
//! * `z_min` / `z_max` / `z_mean` — elevation statistics over the sampled
//!   vertices;
//! * `surf_len` — the true 3D (over-the-surface) length, summed over all line
//!   and ring segments;
//! * `avg_slope` — the length-weighted mean slope in degrees.
//!
//! `attributes` selects which of these columns to add (default: all). The output
//! geometry carries the draped Z. No-data samples are filled from the nearest
//! valid neighbour when possible; a vertex with no usable surface value keeps no
//! Z and is skipped in the statistics. The vector and raster must share a CRS.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Ring};

use crate::common::load_input_raster;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct InterpolateShapeTool;

impl Tool for InterpolateShapeTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "interpolate_shape",
            display_name: "Interpolate Shape",
            summary: "Drape points/lines/polygons on a surface raster: densify, sample Z per vertex (bilinear or nearest), write 3D geometry, and add surface metrics (z_min/z_max/z_mean, 3D surface length, average slope), like ArcGIS Interpolate Shape / Add Surface Information.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer (points, lines, or polygons) sharing the surface's CRS.",
                    required: true,
                },
                ToolParamSpec {
                    name: "surface",
                    description: "Elevation/surface raster to sample.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "sample_distance",
                    description: "Densification interval in CRS units; long segments are split so they follow the surface. Default: the raster cell size.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "Surface sampling: 'bilinear' (default) or 'nearest'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "attributes",
                    description: "Comma list of metric columns to add: z_min,z_max,z_mean,surf_len,avg_slope. Default: all.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based surface band to sample (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "surface")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let surface = require_str(args, "surface")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        ctx.progress.info("reading surface and vector");
        let dem = load_input_raster(surface)?;
        let mut layer = load_input_layer(input)?;

        let sample_distance = prm
            .sample_distance
            .unwrap_or_else(|| dem.cell_size_x.min(dem.cell_size_y))
            .max(f64::MIN_POSITIVE);

        // Add the requested metric fields.
        for m in &prm.attributes {
            layer.add_field(FieldDef::new(m.field_name(), FieldType::Float));
        }

        ctx.progress
            .info(&format!("draping {} feature(s)", layer.len()));

        let mut draped = 0usize;
        let mut skipped = 0usize;
        for feature in layer.features.iter_mut() {
            let Some(geom) = feature.geometry.as_ref() else {
                skipped += 1;
                continue;
            };
            let (new_geom, metrics) =
                drape_geometry(geom, &dem, prm.band, sample_distance, prm.method);
            feature.geometry = Some(new_geom);
            // Append metric values in the same order they were added as fields.
            for m in &prm.attributes {
                feature.attributes.push(metric_value(m, &metrics));
            }
            if metrics.count > 0 {
                draped += 1;
            } else {
                skipped += 1;
            }
        }
        layer.extent = None;

        ctx.progress
            .info(&format!("{draped} draped, {skipped} without surface data"));

        let feature_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("draped_count".to_string(), json!(draped));
        outputs.insert("skipped_count".to_string(), json!(skipped));
        outputs.insert("sample_distance".to_string(), json!(sample_distance));
        Ok(ToolRunResult { outputs })
    }
}

// ── Draping ──────────────────────────────────────────────────────────────────

#[derive(Default)]
struct Metrics {
    count: usize,
    z_sum: f64,
    z_min: f64,
    z_max: f64,
    surf_len: f64,
    slope_len_weighted: f64,
    horiz_len: f64,
}

impl Metrics {
    fn new() -> Self {
        Metrics {
            z_min: f64::INFINITY,
            z_max: f64::NEG_INFINITY,
            ..Default::default()
        }
    }
    fn add_vertex(&mut self, z: f64) {
        self.count += 1;
        self.z_sum += z;
        self.z_min = self.z_min.min(z);
        self.z_max = self.z_max.max(z);
    }
    /// Accumulates one segment's 3D length and slope (both endpoints must have Z).
    fn add_segment(&mut self, horiz: f64, dz: f64) {
        if horiz <= 0.0 {
            return;
        }
        let seg3d = (horiz * horiz + dz * dz).sqrt();
        self.surf_len += seg3d;
        let slope_deg = (dz.abs() / horiz).atan().to_degrees();
        self.slope_len_weighted += slope_deg * horiz;
        self.horiz_len += horiz;
    }
    fn z_mean(&self) -> f64 {
        if self.count > 0 {
            self.z_sum / self.count as f64
        } else {
            f64::NAN
        }
    }
    fn avg_slope(&self) -> f64 {
        if self.horiz_len > 0.0 {
            self.slope_len_weighted / self.horiz_len
        } else {
            f64::NAN
        }
    }
}

/// Drapes a geometry: returns the Z-valued geometry and its surface metrics.
fn drape_geometry(
    geom: &Geometry,
    dem: &Raster,
    band: isize,
    step: f64,
    method: Method,
) -> (Geometry, Metrics) {
    let mut m = Metrics::new();
    let g = match geom {
        Geometry::Point(c) => {
            let z = sample(dem, band, c.x, c.y, method);
            if let Some(z) = z {
                m.add_vertex(z);
                Geometry::Point(Coord::xyz(c.x, c.y, z))
            } else {
                Geometry::Point(Coord::xy(c.x, c.y))
            }
        }
        Geometry::MultiPoint(cs) => {
            let out = cs
                .iter()
                .map(|c| match sample(dem, band, c.x, c.y, method) {
                    Some(z) => {
                        m.add_vertex(z);
                        Coord::xyz(c.x, c.y, z)
                    }
                    None => Coord::xy(c.x, c.y),
                })
                .collect();
            Geometry::MultiPoint(out)
        }
        Geometry::LineString(cs) => {
            Geometry::LineString(drape_line(cs, dem, band, step, method, &mut m))
        }
        Geometry::MultiLineString(lines) => Geometry::MultiLineString(
            lines
                .iter()
                .map(|l| drape_line(l, dem, band, step, method, &mut m))
                .collect(),
        ),
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            let ext = drape_ring(exterior, dem, band, step, method, &mut m);
            let holes = interiors
                .iter()
                .map(|r| drape_ring(r, dem, band, step, method, &mut m))
                .collect();
            Geometry::Polygon {
                exterior: ext,
                interiors: holes,
            }
        }
        Geometry::MultiPolygon(parts) => Geometry::MultiPolygon(
            parts
                .iter()
                .map(|(ext, holes)| {
                    let e = drape_ring(ext, dem, band, step, method, &mut m);
                    let h = holes
                        .iter()
                        .map(|r| drape_ring(r, dem, band, step, method, &mut m))
                        .collect();
                    (e, h)
                })
                .collect(),
        ),
        other => other.clone(),
    };
    (g, m)
}

fn drape_line(
    cs: &[Coord],
    dem: &Raster,
    band: isize,
    step: f64,
    method: Method,
    m: &mut Metrics,
) -> Vec<Coord> {
    let dense = densify(cs, step, false);
    drape_chain(&dense, dem, band, method, m)
}

fn drape_ring(
    ring: &Ring,
    dem: &Raster,
    band: isize,
    step: f64,
    method: Method,
    m: &mut Metrics,
) -> Ring {
    let dense = densify(ring.coords(), step, true);
    Ring::new(drape_chain(&dense, dem, band, method, m))
}

/// Samples Z for each vertex of a chain, accumulating vertex and segment metrics.
fn drape_chain(
    coords: &[Coord],
    dem: &Raster,
    band: isize,
    method: Method,
    m: &mut Metrics,
) -> Vec<Coord> {
    let mut out: Vec<Coord> = Vec::with_capacity(coords.len());
    let mut prev: Option<(f64, f64, f64)> = None;
    for c in coords {
        match sample(dem, band, c.x, c.y, method) {
            Some(z) => {
                m.add_vertex(z);
                if let Some((px, py, pz)) = prev {
                    let horiz = (c.x - px).hypot(c.y - py);
                    m.add_segment(horiz, z - pz);
                }
                prev = Some((c.x, c.y, z));
                out.push(Coord::xyz(c.x, c.y, z));
            }
            None => {
                // Keep the 2D vertex; break the segment chain across the gap.
                prev = None;
                out.push(Coord::xy(c.x, c.y));
            }
        }
    }
    out
}

/// Densifies a coordinate chain so no segment is longer than `max_len`.
fn densify(coords: &[Coord], max_len: f64, closed: bool) -> Vec<Coord> {
    let n = coords.len();
    if n < 2 || max_len <= 0.0 {
        return coords.to_vec();
    }
    let edges = if closed { n } else { n - 1 };
    let mut out = Vec::with_capacity(n * 2);
    for i in 0..edges {
        let (ax, ay) = (coords[i].x, coords[i].y);
        let (bx, by) = (coords[(i + 1) % n].x, coords[(i + 1) % n].y);
        out.push(Coord::xy(ax, ay));
        let d = (bx - ax).hypot(by - ay);
        let pieces = (d / max_len).ceil().max(1.0) as usize;
        for j in 1..pieces {
            let t = j as f64 / pieces as f64;
            out.push(Coord::xy(ax + (bx - ax) * t, ay + (by - ay) * t));
        }
    }
    if !closed {
        out.push(Coord::xy(coords[n - 1].x, coords[n - 1].y));
    }
    out
}

// ── Surface sampling ─────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Bilinear,
    Nearest,
}

/// Samples the surface at world `(x, y)`. Bilinear blends the four surrounding
/// cell centres (falling back to the nearest valid corner when some are
/// no-data); nearest returns the containing cell. Returns `None` when no usable
/// value is available.
fn sample(dem: &Raster, band: isize, x: f64, y: f64, method: Method) -> Option<f64> {
    match method {
        Method::Nearest => sample_nearest(dem, band, x, y),
        Method::Bilinear => sample_bilinear(dem, band, x, y),
    }
}

fn cell_value(dem: &Raster, band: isize, row: isize, col: isize) -> Option<f64> {
    if row < 0 || col < 0 || row >= dem.rows as isize || col >= dem.cols as isize {
        return None;
    }
    let v = dem.get(band, row, col);
    if v == dem.nodata || v.is_nan() {
        None
    } else {
        Some(v)
    }
}

fn sample_nearest(dem: &Raster, band: isize, x: f64, y: f64) -> Option<f64> {
    let (col, row) = dem.world_to_pixel(x, y)?;
    cell_value(dem, band, row, col)
}

fn sample_bilinear(dem: &Raster, band: isize, x: f64, y: f64) -> Option<f64> {
    // Fractional position in cell-center space (col/row of the cell whose centre
    // is just north-west of the point).
    let fx = (x - dem.x_min) / dem.cell_size_x - 0.5;
    let fy = (dem.y_max() - y) / dem.cell_size_y - 0.5;
    let col0 = fx.floor() as isize;
    let row0 = fy.floor() as isize;
    let tx = fx - col0 as f64;
    let ty = fy - row0 as f64;

    let v00 = cell_value(dem, band, row0, col0);
    let v01 = cell_value(dem, band, row0, col0 + 1);
    let v10 = cell_value(dem, band, row0 + 1, col0);
    let v11 = cell_value(dem, band, row0 + 1, col0 + 1);

    // All four valid: full bilinear blend.
    if let (Some(a), Some(b), Some(c), Some(d)) = (v00, v01, v10, v11) {
        let top = a * (1.0 - tx) + b * tx;
        let bot = c * (1.0 - tx) + d * tx;
        return Some(top * (1.0 - ty) + bot * ty);
    }
    // Otherwise fall back to the nearest valid corner (by bilinear weight).
    let candidates = [
        (v00, (1.0 - tx) * (1.0 - ty)),
        (v01, tx * (1.0 - ty)),
        (v10, (1.0 - tx) * ty),
        (v11, tx * ty),
    ];
    candidates
        .iter()
        .filter_map(|(v, w)| v.map(|v| (v, *w)))
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(v, _)| v)
}

// ── Metric selection ─────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Metric {
    ZMin,
    ZMax,
    ZMean,
    SurfLen,
    AvgSlope,
}

impl Metric {
    fn field_name(&self) -> &'static str {
        match self {
            Metric::ZMin => "z_min",
            Metric::ZMax => "z_max",
            Metric::ZMean => "z_mean",
            Metric::SurfLen => "surf_len",
            Metric::AvgSlope => "avg_slope",
        }
    }
    fn parse(s: &str) -> Option<Metric> {
        match s.trim().to_ascii_lowercase().as_str() {
            "z_min" => Some(Metric::ZMin),
            "z_max" => Some(Metric::ZMax),
            "z_mean" => Some(Metric::ZMean),
            "surf_len" | "surface_length" => Some(Metric::SurfLen),
            "avg_slope" => Some(Metric::AvgSlope),
            _ => None,
        }
    }
}

fn metric_value(m: &Metric, metrics: &Metrics) -> FieldValue {
    let v = match m {
        Metric::ZMin if metrics.count > 0 => metrics.z_min,
        Metric::ZMax if metrics.count > 0 => metrics.z_max,
        Metric::ZMean => metrics.z_mean(),
        Metric::SurfLen => metrics.surf_len,
        Metric::AvgSlope => metrics.avg_slope(),
        _ => f64::NAN,
    };
    FieldValue::Float(v)
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    sample_distance: Option<f64>,
    method: Method,
    attributes: Vec<Metric>,
    band: isize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let sample_distance = parse_optional_f64(args, "sample_distance")?;
    if let Some(v) = sample_distance {
        if !(v > 0.0 && v.is_finite()) {
            return Err(ToolError::Validation(
                "'sample_distance' must be a positive number".to_string(),
            ));
        }
    }
    let method = match parse_optional_str(args, "method")? {
        None => Method::Bilinear,
        Some(s) => match s.trim().to_ascii_lowercase().as_str() {
            "bilinear" => Method::Bilinear,
            "nearest" => Method::Nearest,
            other => {
                return Err(ToolError::Validation(format!(
                    "'method' must be 'bilinear' or 'nearest', got '{other}'"
                )))
            }
        },
    };
    let attributes = match parse_optional_str(args, "attributes")? {
        None => vec![
            Metric::ZMin,
            Metric::ZMax,
            Metric::ZMean,
            Metric::SurfLen,
            Metric::AvgSlope,
        ],
        Some(s) => {
            let mut v = Vec::new();
            for part in s.split(',').filter(|p| !p.trim().is_empty()) {
                let m = Metric::parse(part)
                    .ok_or_else(|| ToolError::Validation(format!("unknown attribute '{part}'")))?;
                if !v.contains(&m) {
                    v.push(m);
                }
            }
            v
        }
    };
    let band_1based = parse_optional_f64(args, "band")?
        .map(|v| v as i64)
        .unwrap_or(1);
    if band_1based < 1 {
        return Err(ToolError::Validation("'band' must be >= 1".to_string()));
    }
    Ok(Params {
        sample_distance,
        method,
        attributes,
        band: (band_1based - 1) as isize,
    })
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn parse_optional_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, Raster, RasterConfig};
    use wbvector::{memory_store, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// A tilted plane surface: elevation = 10 * col (so slope is known).
    fn ramp_dem(cols: usize, rows: usize, cell: f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: cell,
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
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, 10.0 * col as f64)
                    .unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn line_layer(coords: &[(f64, f64)]) -> String {
        let mut l = Layer::new("lines")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("name", FieldType::Text));
        let cs = coords.iter().map(|&(x, y)| Coord::xy(x, y)).collect();
        l.add_feature(Some(Geometry::line_string(cs)), &[("name", "l".into())])
            .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = InterpolateShapeTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    #[test]
    fn drapes_a_point_with_z() {
        // 10x10 ramp, cell 1: elevation at col c is 10*c. Point at x=5.0 (col 4
        // or 5) sampled bilinearly ~ 10*(5.0-0.5)=45 at cell-center space.
        let dem = ramp_dem(10, 10, 1.0);
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_feature(Some(Geometry::point(5.0, 5.0)), &[]).unwrap();
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);
        let (_o, layer) = run(json!({ "input": input, "surface": dem }));
        match layer.features[0].geometry.as_ref().unwrap() {
            Geometry::Point(c) => {
                assert!(c.z.is_some(), "point should be draped with Z");
                let z = c.z.unwrap();
                assert!(
                    (z - 45.0).abs() < 1e-6,
                    "expected z≈45 on the ramp, got {z}"
                );
            }
            other => panic!("expected point, got {other:?}"),
        }
    }

    #[test]
    fn line_gets_surface_length_and_slope() {
        // Ramp: elevation rises 10 per column. A horizontal line along y at
        // x 0.5..8.5 climbs 10 per unit x -> slope atan(10)=84.3 deg, and the
        // 3D length is sqrt(1+100) per unit horizontal.
        let dem = ramp_dem(10, 10, 1.0);
        let input = line_layer(&[(0.5, 5.0), (8.5, 5.0)]);
        let (_o, layer) = run(json!({
            "input": input, "surface": dem, "method": "bilinear",
        }));
        let g = |name: &str| {
            let idx = layer.schema.field_index(name).unwrap();
            layer.features[0].attributes[idx].as_f64().unwrap()
        };
        let horiz = 8.0;
        let expected_surf = horiz * (1.0 + 100.0f64).sqrt();
        assert!(
            (g("surf_len") - expected_surf).abs() < 0.5,
            "surf_len {} vs expected {}",
            g("surf_len"),
            expected_surf
        );
        let slope = g("avg_slope");
        assert!((slope - 84.3).abs() < 1.0, "avg_slope {slope} vs ~84.3");
        // z rises from ~0 to ~80 across the ramp.
        assert!(
            g("z_min") < 10.0 && g("z_max") >= 79.0,
            "z range off: {}..{}",
            g("z_min"),
            g("z_max")
        );
    }

    #[test]
    fn selects_attribute_subset() {
        let dem = ramp_dem(10, 10, 1.0);
        let input = line_layer(&[(0.5, 5.0), (8.5, 5.0)]);
        let (_o, layer) = run(json!({
            "input": input, "surface": dem, "attributes": "z_mean,surf_len",
        }));
        assert!(layer.schema.field_index("z_mean").is_some());
        assert!(layer.schema.field_index("surf_len").is_some());
        assert!(
            layer.schema.field_index("avg_slope").is_none(),
            "avg_slope not requested"
        );
        assert!(layer.schema.field_index("z_min").is_none());
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            InterpolateShapeTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err());
        assert!(
            bad(json!({ "input": "a.geojson", "surface": "d.tif", "method": "spline" })).is_err()
        );
        assert!(
            bad(json!({ "input": "a.geojson", "surface": "d.tif", "attributes": "z_bogus" }))
                .is_err()
        );
        assert!(bad(json!({ "input": "a.geojson", "surface": "d.tif" })).is_ok());
    }
}
