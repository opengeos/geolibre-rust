//! GeoLibre tool: append surface-derived attributes to vector features.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Add Surface Information* (3D Analyst).
//! Where `interpolate_shape` *drapes* geometry into 3D and rewrites the shapes,
//! this tool leaves the 2D geometry untouched and only **appends attribute
//! columns** describing the surface underneath each feature — the way a GIS
//! analyst annotates an existing trail, parcel, or footprint layer with terrain
//! statistics without changing the shapes.
//!
//! The requested `properties` (ArcGIS keyword groups) map onto fields:
//!
//! * `Z` → `Z` — the surface value at a point feature.
//! * `SURFACE_LENGTH` → `SLength` — the true 3D (over-the-surface) length of a
//!   line, or the 3D perimeter of a polygon.
//! * `SURFACE_AREA` → `SArea` — the 3D surface area draped over a polygon's
//!   footprint (∑ cell planimetric area · √(1+fx²+fy²), the surface-area-ratio
//!   integral).
//! * `MIN_MAX_MEAN_Z` → `Min_Z` / `Max_Z` / `Mean_Z` — elevation statistics
//!   (over line vertices, or over the raster cells inside a polygon).
//! * `MIN_MAX_AVG_SLOPE` → `Min_Slope` / `Max_Slope` / `Avg_Slope` — slope in
//!   degrees (per line segment, or per cell inside a polygon).
//!
//! Line elevation/length come from densifying each edge to `sample_distance`
//! and bilinearly (or nearest) sampling the surface at every vertex. Polygon
//! area/z/slope statistics come from the raster cells whose centre falls inside
//! the polygon (holes respected); `Avg_Slope` for lines is length-weighted, for
//! polygons it is the per-cell mean. Which fields apply to which geometry kind
//! follows ArcGIS: points get `Z`, lines get length + z + slope, polygons get
//! area + length + z + slope. The vector and raster must share a CRS.

use std::collections::BTreeMap;

use geo::{Area, Contains};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Ring};

use crate::common::load_input_raster;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct AddSurfaceInformationTool;

impl Tool for AddSurfaceInformationTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "add_surface_information",
            display_name: "Add Surface Information",
            summary: "Append per-feature surface statistics from an elevation raster — Z, 3D surface length, 3D surface area, min/max/mean Z, and min/max/avg slope — as attribute fields without altering the input geometry, like ArcGIS Add Surface Information.",
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
                    name: "properties",
                    description: "Comma list of surface property groups to append: Z, SURFACE_LENGTH, SURFACE_AREA, MIN_MAX_MEAN_Z, MIN_MAX_AVG_SLOPE. Default: all applicable to the layer's geometry.",
                    required: false,
                },
                ToolParamSpec {
                    name: "sample_distance",
                    description: "Densification interval in CRS units for line/perimeter sampling; long segments are split so they follow the surface. Default: the raster cell size.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "Vertex surface sampling: 'bilinear' (default) or 'nearest'.",
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

        let step = prm
            .sample_distance
            .unwrap_or_else(|| dem.cell_size_x.min(dem.cell_size_y))
            .max(f64::MIN_POSITIVE);

        // Determine which geometry kinds are present so we only add applicable
        // fields (ArcGIS adds Z only for points, SArea only for polygons, etc.).
        let mut kinds = Kinds::default();
        for f in &layer.features {
            if let Some(g) = f.geometry.as_ref() {
                kinds.observe(g);
            }
        }

        // Fields to add, in a stable order, filtered by request and applicability.
        let fields = resolve_fields(&prm.properties, &kinds);
        if fields.is_empty() {
            return Err(ToolError::Validation(
                "no requested 'properties' apply to the input geometry".to_string(),
            ));
        }
        for f in &fields {
            layer.add_field(FieldDef::new(f.name(), FieldType::Float));
        }

        ctx.progress
            .info(&format!("annotating {} feature(s)", layer.len()));

        let mut annotated = 0usize;
        let mut skipped = 0usize;
        for feature in layer.features.iter_mut() {
            let info = match feature.geometry.as_ref() {
                Some(g) => compute_surface_info(g, &dem, prm.band, step, prm.method),
                None => SurfaceInfo::default(),
            };
            for f in &fields {
                feature.attributes.push(field_value(*f, &info));
            }
            if info.any_computed() {
                annotated += 1;
            } else {
                skipped += 1;
            }
        }

        ctx.progress.info(&format!(
            "{annotated} annotated, {skipped} without surface data"
        ));

        let feature_count = layer.len();
        let field_names: Vec<&str> = fields.iter().map(|f| f.name()).collect();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("annotated_count".to_string(), json!(annotated));
        outputs.insert("skipped_count".to_string(), json!(skipped));
        outputs.insert("fields_added".to_string(), json!(field_names));
        outputs.insert("sample_distance".to_string(), json!(step));
        Ok(ToolRunResult { outputs })
    }
}

// ── Geometry-kind bookkeeping ─────────────────────────────────────────────────

#[derive(Default)]
struct Kinds {
    point: bool,
    line: bool,
    polygon: bool,
}

impl Kinds {
    fn observe(&mut self, g: &Geometry) {
        match g {
            Geometry::Point(_) | Geometry::MultiPoint(_) => self.point = true,
            Geometry::LineString(_) | Geometry::MultiLineString(_) => self.line = true,
            Geometry::Polygon { .. } | Geometry::MultiPolygon(_) => self.polygon = true,
            _ => {}
        }
    }
}

// ── Per-feature surface statistics ────────────────────────────────────────────

#[derive(Default)]
struct SurfaceInfo {
    z: Option<f64>,
    surface_length: Option<f64>,
    surface_area: Option<f64>,
    min_z: Option<f64>,
    max_z: Option<f64>,
    mean_z: Option<f64>,
    min_slope: Option<f64>,
    max_slope: Option<f64>,
    avg_slope: Option<f64>,
}

impl SurfaceInfo {
    fn any_computed(&self) -> bool {
        self.z.is_some()
            || self.surface_length.is_some()
            || self.surface_area.is_some()
            || self.mean_z.is_some()
            || self.avg_slope.is_some()
    }
}

fn compute_surface_info(
    geom: &Geometry,
    dem: &Raster,
    band: isize,
    step: f64,
    method: Method,
) -> SurfaceInfo {
    match geom {
        Geometry::Point(c) => point_info(std::slice::from_ref(c), dem, band, method),
        Geometry::MultiPoint(cs) => point_info(cs, dem, band, method),
        Geometry::LineString(cs) => line_info(
            std::slice::from_ref(&cs.as_slice()),
            dem,
            band,
            step,
            method,
        ),
        Geometry::MultiLineString(lines) => {
            let refs: Vec<&[Coord]> = lines.iter().map(|l| l.as_slice()).collect();
            line_info(&refs, dem, band, step, method)
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => polygon_info(
            std::slice::from_ref(&(exterior, interiors)),
            dem,
            band,
            step,
            method,
        ),
        Geometry::MultiPolygon(parts) => {
            let refs: Vec<(&Ring, &Vec<Ring>)> = parts.iter().map(|(e, h)| (e, h)).collect();
            polygon_info(&refs, dem, band, step, method)
        }
        _ => SurfaceInfo::default(),
    }
}

/// Point/MultiPoint: sample the surface at each point.
fn point_info(coords: &[Coord], dem: &Raster, band: isize, method: Method) -> SurfaceInfo {
    let mut acc = ZAcc::new();
    for c in coords {
        if let Some(z) = sample(dem, band, c.x, c.y, method) {
            acc.add(z);
        }
    }
    let mut info = SurfaceInfo::default();
    if acc.count > 0 {
        info.z = Some(acc.mean());
        info.min_z = Some(acc.min);
        info.max_z = Some(acc.max);
        info.mean_z = Some(acc.mean());
    }
    info
}

/// Line/MultiLine: densify, sample vertices, accumulate 3D length + z + slope.
fn line_info(
    lines: &[&[Coord]],
    dem: &Raster,
    band: isize,
    step: f64,
    method: Method,
) -> SurfaceInfo {
    let mut z = ZAcc::new();
    let mut slope = SlopeAcc::new();
    let mut surf_len = 0.0;
    for cs in lines {
        let dense = densify(cs, step, false);
        let mut prev: Option<(f64, f64, f64)> = None;
        for c in &dense {
            match sample(dem, band, c.x, c.y, method) {
                Some(zc) => {
                    z.add(zc);
                    if let Some((px, py, pz)) = prev {
                        let horiz = (c.x - px).hypot(c.y - py);
                        if horiz > 0.0 {
                            let seg3d = (horiz * horiz + (zc - pz).powi(2)).sqrt();
                            surf_len += seg3d;
                            let s = ((zc - pz).abs() / horiz).atan().to_degrees();
                            slope.add_weighted(s, horiz);
                        }
                    }
                    prev = Some((c.x, c.y, zc));
                }
                None => prev = None,
            }
        }
    }
    let mut info = SurfaceInfo::default();
    if z.count > 0 {
        info.surface_length = Some(surf_len);
        info.min_z = Some(z.min);
        info.max_z = Some(z.max);
        info.mean_z = Some(z.mean());
    }
    if slope.horiz_len > 0.0 {
        info.min_slope = Some(slope.min);
        info.max_slope = Some(slope.max);
        info.avg_slope = Some(slope.weighted_mean());
    }
    info
}

/// Polygon/MultiPolygon: 3D perimeter from the boundary rings; z/slope/area
/// statistics from the raster cells whose centre lies inside the footprint.
fn polygon_info(
    parts: &[(&Ring, &Vec<Ring>)],
    dem: &Raster,
    band: isize,
    step: f64,
    method: Method,
) -> SurfaceInfo {
    // 3D perimeter (surface length) over all rings.
    let mut surf_len = 0.0;
    let mut perim_z = ZAcc::new();
    for (ext, holes) in parts {
        surf_len += ring_surface_length(ext.coords(), dem, band, step, method, &mut perim_z);
        for h in holes.iter() {
            surf_len += ring_surface_length(h.coords(), dem, band, step, method, &mut perim_z);
        }
    }

    // Cell-based statistics inside the footprint.
    let geo_polys = to_geo_polygons(parts);
    let (mut cell_z, mut cell_slope, mut surface_area) = (ZAcc::new(), SlopeAcc::new(), 0.0);
    let cell_area = dem.cell_size_x * dem.cell_size_y;

    // Bounding box of all parts in world coordinates.
    if let Some((min_x, min_y, max_x, max_y)) = parts_bbox(parts) {
        // Map bbox to inclusive row/col ranges.
        let col_lo = ((min_x - dem.x_min) / dem.cell_size_x).floor() as isize - 1;
        let col_hi = ((max_x - dem.x_min) / dem.cell_size_x).ceil() as isize + 1;
        let row_lo = ((dem.y_max() - max_y) / dem.cell_size_y).floor() as isize - 1;
        let row_hi = ((dem.y_max() - min_y) / dem.cell_size_y).ceil() as isize + 1;
        let col_lo = col_lo.max(0);
        let row_lo = row_lo.max(0);
        let col_hi = col_hi.min(dem.cols as isize - 1);
        let row_hi = row_hi.min(dem.rows as isize - 1);
        for row in row_lo..=row_hi {
            let cy = dem.row_center_y(row);
            for col in col_lo..=col_hi {
                let cx = dem.col_center_x(col);
                let pt = geo::Point::new(cx, cy);
                if !geo_polys.iter().any(|p| p.contains(&pt)) {
                    continue;
                }
                let Some(zc) = cell_value(dem, band, row, col) else {
                    continue;
                };
                cell_z.add(zc);
                let (fx, fy) = cell_gradient(dem, band, row, col, zc);
                let g2 = fx * fx + fy * fy;
                let factor = (1.0 + g2).sqrt();
                surface_area += cell_area * factor;
                let slope_deg = g2.sqrt().atan().to_degrees();
                cell_slope.add_weighted(slope_deg, 1.0);
            }
        }
    }

    let mut info = SurfaceInfo::default();
    if perim_z.count > 0 {
        info.surface_length = Some(surf_len);
    }
    if cell_z.count > 0 {
        info.surface_area = Some(surface_area);
        info.min_z = Some(cell_z.min);
        info.max_z = Some(cell_z.max);
        info.mean_z = Some(cell_z.mean());
        info.min_slope = Some(cell_slope.min);
        info.max_slope = Some(cell_slope.max);
        info.avg_slope = Some(cell_slope.weighted_mean());
    } else {
        // Footprint smaller than a cell: fall back to the centroid so tiny
        // polygons still get z/slope/area rather than nulls.
        if let Some((cx, cy)) = centroid(parts) {
            if let Some((col, row)) = dem.world_to_pixel(cx, cy) {
                if let Some(zc) = cell_value(dem, band, row, col) {
                    let (fx, fy) = cell_gradient(dem, band, row, col, zc);
                    let g2 = fx * fx + fy * fy;
                    let factor = (1.0 + g2).sqrt();
                    let slope_deg = g2.sqrt().atan().to_degrees();
                    let planimetric: f64 = geo_polys.iter().map(|p| p.unsigned_area()).sum();
                    info.surface_area = Some(planimetric * factor);
                    info.min_z = Some(zc);
                    info.max_z = Some(zc);
                    info.mean_z = Some(zc);
                    info.min_slope = Some(slope_deg);
                    info.max_slope = Some(slope_deg);
                    info.avg_slope = Some(slope_deg);
                }
            }
        }
    }
    info
}

/// 3D length of one ring (densified), also feeding perimeter z samples.
fn ring_surface_length(
    coords: &[Coord],
    dem: &Raster,
    band: isize,
    step: f64,
    method: Method,
    z: &mut ZAcc,
) -> f64 {
    let dense = densify(coords, step, true);
    let mut len = 0.0;
    let mut prev: Option<(f64, f64, f64)> = None;
    for c in &dense {
        match sample(dem, band, c.x, c.y, method) {
            Some(zc) => {
                z.add(zc);
                if let Some((px, py, pz)) = prev {
                    let horiz = (c.x - px).hypot(c.y - py);
                    if horiz > 0.0 {
                        len += (horiz * horiz + (zc - pz).powi(2)).sqrt();
                    }
                }
                prev = Some((c.x, c.y, zc));
            }
            None => prev = None,
        }
    }
    len
}

/// Central-difference surface gradient (dz/dx, dz/dy) at a cell, in z-units per
/// CRS unit. No-data neighbours fall back to the centre value (one-sided).
fn cell_gradient(dem: &Raster, band: isize, row: isize, col: isize, center: f64) -> (f64, f64) {
    let e = cell_value(dem, band, row, col + 1).unwrap_or(center);
    let w = cell_value(dem, band, row, col - 1).unwrap_or(center);
    let n = cell_value(dem, band, row - 1, col).unwrap_or(center);
    let s = cell_value(dem, band, row + 1, col).unwrap_or(center);
    let fx = (e - w) / (2.0 * dem.cell_size_x);
    // row increases southward, so dz/dy uses (north - south).
    let fy = (n - s) / (2.0 * dem.cell_size_y);
    (fx, fy)
}

// ── Small accumulators ────────────────────────────────────────────────────────

struct ZAcc {
    count: usize,
    sum: f64,
    min: f64,
    max: f64,
}
impl ZAcc {
    fn new() -> Self {
        ZAcc {
            count: 0,
            sum: 0.0,
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
        }
    }
    fn add(&mut self, v: f64) {
        self.count += 1;
        self.sum += v;
        self.min = self.min.min(v);
        self.max = self.max.max(v);
    }
    fn mean(&self) -> f64 {
        self.sum / self.count as f64
    }
}

struct SlopeAcc {
    horiz_len: f64,
    weighted: f64,
    min: f64,
    max: f64,
}
impl SlopeAcc {
    fn new() -> Self {
        SlopeAcc {
            horiz_len: 0.0,
            weighted: 0.0,
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
        }
    }
    fn add_weighted(&mut self, slope_deg: f64, weight: f64) {
        self.horiz_len += weight;
        self.weighted += slope_deg * weight;
        self.min = self.min.min(slope_deg);
        self.max = self.max.max(slope_deg);
    }
    fn weighted_mean(&self) -> f64 {
        self.weighted / self.horiz_len
    }
}

// ── geo conversion / geometry helpers ─────────────────────────────────────────

fn to_geo_polygons(parts: &[(&Ring, &Vec<Ring>)]) -> Vec<geo::Polygon<f64>> {
    parts
        .iter()
        .map(|(ext, holes)| {
            let shell = ring_to_geo(ext.coords());
            let interiors = holes.iter().map(|h| ring_to_geo(h.coords())).collect();
            geo::Polygon::new(shell, interiors)
        })
        .collect()
}

fn ring_to_geo(coords: &[Coord]) -> geo::LineString<f64> {
    geo::LineString(
        coords
            .iter()
            .map(|c| geo::Coord { x: c.x, y: c.y })
            .collect(),
    )
}

fn parts_bbox(parts: &[(&Ring, &Vec<Ring>)]) -> Option<(f64, f64, f64, f64)> {
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    let mut any = false;
    for (ext, _) in parts {
        for c in ext.coords() {
            any = true;
            min_x = min_x.min(c.x);
            min_y = min_y.min(c.y);
            max_x = max_x.max(c.x);
            max_y = max_y.max(c.y);
        }
    }
    any.then_some((min_x, min_y, max_x, max_y))
}

fn centroid(parts: &[(&Ring, &Vec<Ring>)]) -> Option<(f64, f64)> {
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut n = 0.0;
    for (ext, _) in parts {
        for c in ext.coords() {
            sx += c.x;
            sy += c.y;
            n += 1.0;
        }
    }
    (n > 0.0).then(|| (sx / n, sy / n))
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

// ── Surface sampling (bilinear / nearest) ─────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    Bilinear,
    Nearest,
}

fn sample(dem: &Raster, band: isize, x: f64, y: f64, method: Method) -> Option<f64> {
    match method {
        Method::Nearest => {
            let (col, row) = dem.world_to_pixel(x, y)?;
            cell_value(dem, band, row, col)
        }
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

fn sample_bilinear(dem: &Raster, band: isize, x: f64, y: f64) -> Option<f64> {
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

    if let (Some(a), Some(b), Some(c), Some(d)) = (v00, v01, v10, v11) {
        let top = a * (1.0 - tx) + b * tx;
        let bot = c * (1.0 - tx) + d * tx;
        return Some(top * (1.0 - ty) + bot * ty);
    }
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

// ── Property groups → output fields ───────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Property {
    Z,
    SurfaceLength,
    SurfaceArea,
    MinMaxMeanZ,
    MinMaxAvgSlope,
}

impl Property {
    fn parse(s: &str) -> Option<Property> {
        match s
            .trim()
            .to_ascii_uppercase()
            .replace([' ', '-'], "_")
            .as_str()
        {
            "Z" => Some(Property::Z),
            "SURFACE_LENGTH" | "SLENGTH" => Some(Property::SurfaceLength),
            "SURFACE_AREA" | "SAREA" => Some(Property::SurfaceArea),
            "MIN_MAX_MEAN_Z" | "Z_STATS" => Some(Property::MinMaxMeanZ),
            "MIN_MAX_AVG_SLOPE" | "SLOPE" | "SLOPE_STATS" => Some(Property::MinMaxAvgSlope),
            _ => None,
        }
    }

    fn all() -> [Property; 5] {
        [
            Property::Z,
            Property::SurfaceLength,
            Property::SurfaceArea,
            Property::MinMaxMeanZ,
            Property::MinMaxAvgSlope,
        ]
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Field {
    Z,
    SLength,
    SArea,
    MinZ,
    MaxZ,
    MeanZ,
    MinSlope,
    MaxSlope,
    AvgSlope,
}

impl Field {
    fn name(&self) -> &'static str {
        match self {
            Field::Z => "Z",
            Field::SLength => "SLength",
            Field::SArea => "SArea",
            Field::MinZ => "Min_Z",
            Field::MaxZ => "Max_Z",
            Field::MeanZ => "Mean_Z",
            Field::MinSlope => "Min_Slope",
            Field::MaxSlope => "Max_Slope",
            Field::AvgSlope => "Avg_Slope",
        }
    }
}

/// Expands the requested property groups into concrete fields, keeping only
/// those applicable to a geometry kind present in the layer.
fn resolve_fields(props: &[Property], kinds: &Kinds) -> Vec<Field> {
    let mut out = Vec::new();
    let mut push = |f: Field, ok: bool| {
        if ok && !out.contains(&f) {
            out.push(f);
        }
    };
    for p in props {
        match p {
            Property::Z => push(Field::Z, kinds.point),
            Property::SurfaceLength => push(Field::SLength, kinds.line || kinds.polygon),
            Property::SurfaceArea => push(Field::SArea, kinds.polygon),
            Property::MinMaxMeanZ => {
                let ok = kinds.point || kinds.line || kinds.polygon;
                push(Field::MinZ, ok);
                push(Field::MaxZ, ok);
                push(Field::MeanZ, ok);
            }
            Property::MinMaxAvgSlope => {
                let ok = kinds.line || kinds.polygon;
                push(Field::MinSlope, ok);
                push(Field::MaxSlope, ok);
                push(Field::AvgSlope, ok);
            }
        }
    }
    out
}

fn field_value(f: Field, info: &SurfaceInfo) -> FieldValue {
    let v = match f {
        Field::Z => info.z,
        Field::SLength => info.surface_length,
        Field::SArea => info.surface_area,
        Field::MinZ => info.min_z,
        Field::MaxZ => info.max_z,
        Field::MeanZ => info.mean_z,
        Field::MinSlope => info.min_slope,
        Field::MaxSlope => info.max_slope,
        Field::AvgSlope => info.avg_slope,
    };
    match v {
        Some(x) if x.is_finite() => FieldValue::Float(x),
        _ => FieldValue::Null,
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    properties: Vec<Property>,
    sample_distance: Option<f64>,
    method: Method,
    band: isize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let properties = match parse_optional_str(args, "properties")? {
        None => Property::all().to_vec(),
        Some(s) => {
            let mut v = Vec::new();
            for part in s.split(',').filter(|p| !p.trim().is_empty()) {
                let p = Property::parse(part).ok_or_else(|| {
                    ToolError::Validation(format!("unknown surface property '{part}'"))
                })?;
                if !v.contains(&p) {
                    v.push(p);
                }
            }
            if v.is_empty() {
                Property::all().to_vec()
            } else {
                v
            }
        }
    };

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

    let band_1based = parse_optional_f64(args, "band")?
        .map(|v| v as i64)
        .unwrap_or(1);
    if band_1based < 1 {
        return Err(ToolError::Validation("'band' must be >= 1".to_string()));
    }

    Ok(Params {
        properties,
        sample_distance,
        method,
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

    /// A tilted plane: elevation = slope_per_col * col (units of z per cell).
    fn ramp_dem(cols: usize, rows: usize, cell: f64, per_col: f64) -> String {
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
                r.set(0, row as isize, col as isize, per_col * col as f64)
                    .unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    /// A perfectly flat surface at constant elevation.
    fn flat_dem(cols: usize, rows: usize, cell: f64, z: f64) -> String {
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
                r.set(0, row as isize, col as isize, z).unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> Layer {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = AddSurfaceInformationTool.run(&args, &ctx()).unwrap();
        load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    fn get(layer: &Layer, feat: usize, name: &str) -> Option<f64> {
        let idx = layer.schema.field_index(name)?;
        layer.features[feat].attributes[idx].as_f64()
    }

    #[test]
    fn point_gets_z() {
        // per_col=10, cell=1: elevation at col c = 10*c. Point x=5.0 sampled
        // bilinearly in cell-center space -> 10*(5.0-0.5)=45.
        let dem = ramp_dem(10, 10, 1.0, 10.0);
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_feature(Some(Geometry::point(5.0, 5.0)), &[]).unwrap();
        let id = memory_store::put_vector(l);
        let input = memory_store::make_vector_memory_path(&id);
        let layer = run(json!({ "input": input, "surface": dem }));
        let z = get(&layer, 0, "Z").expect("Z field present");
        assert!((z - 45.0).abs() < 1e-6, "expected Z≈45, got {z}");
        // Geometry must be unchanged (still 2D, no draped Z).
        match layer.features[0].geometry.as_ref().unwrap() {
            Geometry::Point(c) => assert!(c.z.is_none(), "geometry should remain 2D"),
            other => panic!("expected point, got {other:?}"),
        }
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

    #[test]
    fn line_gets_length_z_and_slope() {
        // Ramp climbs 10 per column (per unit x). Horizontal line at y=5 from
        // x=0.5 to x=8.5: slope atan(10)=84.3deg, 3D length 8*sqrt(1+100).
        let dem = ramp_dem(10, 10, 1.0, 10.0);
        let input = line_layer(&[(0.5, 5.0), (8.5, 5.0)]);
        let layer = run(json!({ "input": input, "surface": dem, "method": "bilinear" }));
        let expected = 8.0 * (1.0 + 100.0f64).sqrt();
        let sl = get(&layer, 0, "SLength").unwrap();
        assert!((sl - expected).abs() < 0.5, "SLength {sl} vs {expected}");
        let avg = get(&layer, 0, "Avg_Slope").unwrap();
        assert!((avg - 84.3).abs() < 1.0, "Avg_Slope {avg} vs ~84.3");
        assert!(get(&layer, 0, "Min_Z").unwrap() < 10.0);
        assert!(get(&layer, 0, "Max_Z").unwrap() >= 79.0);
        // A line layer must NOT receive an area or Z field.
        assert!(layer.schema.field_index("SArea").is_none());
        assert!(layer.schema.field_index("Z").is_none());
    }

    fn square_polygon(x0: f64, y0: f64, x1: f64, y1: f64) -> String {
        let mut l = Layer::new("poly")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        let ring = Ring::new(vec![
            Coord::xy(x0, y0),
            Coord::xy(x1, y0),
            Coord::xy(x1, y1),
            Coord::xy(x0, y1),
            Coord::xy(x0, y0),
        ]);
        l.add_feature(
            Some(Geometry::Polygon {
                exterior: ring,
                interiors: vec![],
            }),
            &[],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    #[test]
    fn polygon_flat_surface_area_equals_footprint() {
        // Flat surface: 3D surface area == planimetric area. 10x10 square over a
        // 20x20 grid, cell 1 -> area 100, slope 0, SLength perimeter 40.
        let dem = flat_dem(20, 20, 1.0, 100.0);
        let input = square_polygon(3.0, 3.0, 13.0, 13.0);
        let layer = run(json!({ "input": input, "surface": dem }));
        let area = get(&layer, 0, "SArea").unwrap();
        assert!((area - 100.0).abs() < 6.0, "flat SArea {area} vs ~100");
        assert!(
            (get(&layer, 0, "Avg_Slope").unwrap()).abs() < 1e-6,
            "flat slope ~0"
        );
        assert!((get(&layer, 0, "Mean_Z").unwrap() - 100.0).abs() < 1e-6);
        let sl = get(&layer, 0, "SLength").unwrap();
        assert!((sl - 40.0).abs() < 0.5, "perimeter SLength {sl} vs 40");
    }

    #[test]
    fn polygon_tilted_area_matches_secant_of_slope() {
        // Ramp of 1 z-unit per cell -> gradient 1, slope 45deg, surface-area
        // factor sqrt(1+1)=sqrt(2). A 10x10 footprint -> area ~100*sqrt(2).
        let dem = ramp_dem(20, 20, 1.0, 1.0);
        let input = square_polygon(4.0, 4.0, 14.0, 14.0);
        let layer = run(json!({ "input": input, "surface": dem }));
        let area = get(&layer, 0, "SArea").unwrap();
        let expected = 100.0 * 2.0f64.sqrt();
        assert!(
            (area - expected).abs() < 12.0,
            "tilted SArea {area} vs ~{expected}"
        );
        let avg = get(&layer, 0, "Avg_Slope").unwrap();
        assert!((avg - 45.0).abs() < 1.0, "Avg_Slope {avg} vs ~45");
    }

    #[test]
    fn selects_property_subset() {
        let dem = ramp_dem(20, 20, 1.0, 1.0);
        let input = square_polygon(4.0, 4.0, 14.0, 14.0);
        let layer = run(json!({
            "input": input, "surface": dem, "properties": "SURFACE_AREA,MIN_MAX_MEAN_Z",
        }));
        assert!(layer.schema.field_index("SArea").is_some());
        assert!(layer.schema.field_index("Mean_Z").is_some());
        assert!(
            layer.schema.field_index("Avg_Slope").is_none(),
            "slope not requested"
        );
        assert!(layer.schema.field_index("SLength").is_none());
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            AddSurfaceInformationTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson" })).is_err());
        assert!(
            bad(json!({ "input": "a.geojson", "surface": "d.tif", "method": "spline" })).is_err()
        );
        assert!(
            bad(json!({ "input": "a.geojson", "surface": "d.tif", "properties": "BOGUS" }))
                .is_err()
        );
        assert!(
            bad(json!({ "input": "a.geojson", "surface": "d.tif", "sample_distance": -1 }))
                .is_err()
        );
        assert!(bad(json!({ "input": "a.geojson", "surface": "d.tif" })).is_ok());
    }
}
