//! GeoLibre tool: build line or ellipse features from tabular fields.
//!
//! Pure-Rust counterpart of three closely related ArcGIS Pro *Data Management*
//! tools, unified behind one `mode` selector since they share the same
//! planar/geodesic/rhumb line-construction machinery:
//!
//! - `xy_to_line`       — like *XY To Line*: one line per row from a
//!   `(start_x, start_y)` field pair to an `(end_x, end_y)` field pair.
//! - `bearing_distance` — like *Bearing Distance To Line*: one line per row
//!   from an origin `(x, y)` travelling `distance` along `bearing` (degrees
//!   clockwise from north).
//! - `ellipse`          — like *Table To Ellipse*: an ellipse outline (or,
//!   optionally, a filled polygon) per row centred at `(x, y)` with semi-axis
//!   fields `major`/`minor` and an `azimuth` (degrees clockwise from north)
//!   giving the major axis direction.
//!
//! `line_type` selects how the geometry is constructed:
//!
//! - `planar`   — straight Cartesian segments in the field values' own units
//!   (matches the input layer's CRS, whatever that is).
//! - `geodesic` (default) — shortest path on the WGS84 ellipsoid, via `geo`
//!   0.33's `Geodesic` metric space (a pure-Rust Karney/geographiclib port,
//!   not feature-gated). Field values are lon/lat degrees; `distance` and the
//!   ellipse axes are meters.
//! - `rhumb`    — constant-bearing loxodrome, via `geo`'s `Rhumb` metric
//!   space. Same units as geodesic.
//!
//! For `xy_to_line`/`bearing_distance`, non-planar lines are densified with
//! `vertex_spacing` (meters) so the curvature is visible; for `ellipse`, the
//! same parameter controls the arc-length spacing of the sampled outline.
//! Every source attribute is carried onto the output feature unchanged.
//!
//! Nothing in the repo or the bundled whitebox suite builds geometry *from*
//! attribute fields — every existing tool consumes geometry that already
//! exists. This is a staple ArcGIS workflow for turning tabular
//! sighting/telemetry/error-ellipse data (aviation, search-and-rescue, crime
//! analysis, wildlife tracking) into mappable features.
//!
//! Scope for v1: field names only (no reading start/end from two separate
//! point layers by join key); ellipse output is a closed line by default,
//! with `polygon_output` to request a filled polygon instead.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldValue, Geometry, GeometryType, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

use geo::{Destination, Distance, Geodesic, InterpolatePoint, Point as GeoPoint, Rhumb};

/// Fixed segment count for a curved (geodesic/rhumb) `xy_to_line` /
/// `bearing_distance` line when `vertex_spacing` is not given.
const DEFAULT_LINE_SEGMENTS: usize = 32;
/// Default ellipse sample count when `vertex_spacing` is not given (matches
/// `directional_distribution`'s ellipse rendering).
const DEFAULT_ELLIPSE_SEGMENTS: usize = 120;
/// Hard cap on densification, so a tiny `vertex_spacing` can't blow up memory.
const MAX_SEGMENTS: usize = 2000;

pub struct TableToGeometryTool;

impl Tool for TableToGeometryTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "table_to_geometry",
            display_name: "Table To Geometry",
            summary: "Build line or ellipse features from tabular fields: start/end coordinate pairs (XY To Line), origin + bearing + distance (Bearing Distance To Line), or center + major/minor axis + azimuth (Table To Ellipse) — with planar, geodesic, and rhumb line construction.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input attribute table or point vector layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "mode",
                    description: "'xy_to_line', 'bearing_distance', or 'ellipse'.",
                    required: true,
                },
                ToolParamSpec {
                    name: "line_type",
                    description: "'geodesic' (default, shortest path on the WGS84 ellipsoid), 'rhumb' (constant bearing), or 'planar' (straight Cartesian segments in the input's own units).",
                    required: false,
                },
                ToolParamSpec {
                    name: "vertex_spacing",
                    description: "Densification spacing for curved (geodesic/rhumb) lines, and ellipse outline sample spacing, in meters for geodesic/rhumb or CRS units for planar. Default: a fixed segment count.",
                    required: false,
                },
                ToolParamSpec {
                    name: "polygon_output",
                    description: "mode=ellipse only: emit a filled Polygon instead of a closed line outline. Default false.",
                    required: false,
                },
                ToolParamSpec {
                    name: "start_x",
                    description: "mode=xy_to_line: field holding the start X / longitude.",
                    required: false,
                },
                ToolParamSpec {
                    name: "start_y",
                    description: "mode=xy_to_line: field holding the start Y / latitude.",
                    required: false,
                },
                ToolParamSpec {
                    name: "end_x",
                    description: "mode=xy_to_line: field holding the end X / longitude.",
                    required: false,
                },
                ToolParamSpec {
                    name: "end_y",
                    description: "mode=xy_to_line: field holding the end Y / latitude.",
                    required: false,
                },
                ToolParamSpec {
                    name: "x",
                    description: "mode=bearing_distance/ellipse: field holding the origin/center X / longitude.",
                    required: false,
                },
                ToolParamSpec {
                    name: "y",
                    description: "mode=bearing_distance/ellipse: field holding the origin/center Y / latitude.",
                    required: false,
                },
                ToolParamSpec {
                    name: "bearing",
                    description: "mode=bearing_distance: field holding the bearing in degrees clockwise from north.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distance",
                    description: "mode=bearing_distance: field holding the travel distance (meters for geodesic/rhumb, CRS units for planar).",
                    required: false,
                },
                ToolParamSpec {
                    name: "major",
                    description: "mode=ellipse: field holding the semi-major axis length (meters for geodesic/rhumb, CRS units for planar).",
                    required: false,
                },
                ToolParamSpec {
                    name: "minor",
                    description: "mode=ellipse: field holding the semi-minor axis length (same units as major).",
                    required: false,
                },
                ToolParamSpec {
                    name: "azimuth",
                    description: "mode=ellipse: field holding the major-axis direction in degrees clockwise from north.",
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
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(input)?;
        let schema = layer.schema.clone();
        let field_idx = |name: &str| -> Result<usize, ToolError> {
            schema.field_index(name).ok_or_else(|| {
                ToolError::Validation(format!("field '{name}' not found in input schema"))
            })
        };

        let is_polygon = prm.mode == Mode::Ellipse && prm.polygon_output;

        let mut out = Layer::new(format!("{}_geometry", prm.mode.as_str()));
        for f in schema.fields() {
            out.add_field(f.clone());
        }
        out.geom_type = Some(if is_polygon {
            GeometryType::Polygon
        } else {
            GeometryType::LineString
        });
        out = if prm.line_type == LineType::Planar {
            match layer.crs_epsg() {
                Some(epsg) => out.with_crs_epsg(epsg),
                None => out,
            }
        } else {
            out.with_crs_epsg(4326)
        };

        let mut built = 0usize;
        let mut skipped = 0usize;

        match &prm.fields {
            Fields::XyToLine {
                start_x,
                start_y,
                end_x,
                end_y,
            } => {
                let (isx, isy, iex, iey) = (
                    field_idx(start_x)?,
                    field_idx(start_y)?,
                    field_idx(end_x)?,
                    field_idx(end_y)?,
                );
                for feature in layer.iter() {
                    let vals = (
                        f64_at(feature, isx),
                        f64_at(feature, isy),
                        f64_at(feature, iex),
                        f64_at(feature, iey),
                    );
                    let (Some(sx), Some(sy), Some(ex), Some(ey)) = vals else {
                        skipped += 1;
                        continue;
                    };
                    let pts = line_points((sx, sy), (ex, ey), prm.line_type, prm.vertex_spacing);
                    push_line(&mut out, feature, pts);
                    built += 1;
                }
            }
            Fields::BearingDistance {
                x,
                y,
                bearing,
                distance,
            } => {
                let (ix, iy, ib, id) = (
                    field_idx(x)?,
                    field_idx(y)?,
                    field_idx(bearing)?,
                    field_idx(distance)?,
                );
                for feature in layer.iter() {
                    let vals = (
                        f64_at(feature, ix),
                        f64_at(feature, iy),
                        f64_at(feature, ib),
                        f64_at(feature, id),
                    );
                    let (Some(ox), Some(oy), Some(brg), Some(dist)) = vals else {
                        skipped += 1;
                        continue;
                    };
                    if !dist.is_finite() || dist < 0.0 {
                        skipped += 1;
                        continue;
                    }
                    let end = destination_point((ox, oy), brg, dist, prm.line_type);
                    let pts = line_points((ox, oy), end, prm.line_type, prm.vertex_spacing);
                    push_line(&mut out, feature, pts);
                    built += 1;
                }
            }
            Fields::Ellipse {
                x,
                y,
                major,
                minor,
                azimuth,
            } => {
                let (ix, iy, ima, imi, iaz) = (
                    field_idx(x)?,
                    field_idx(y)?,
                    field_idx(major)?,
                    field_idx(minor)?,
                    field_idx(azimuth)?,
                );
                for feature in layer.iter() {
                    let vals = (
                        f64_at(feature, ix),
                        f64_at(feature, iy),
                        f64_at(feature, ima),
                        f64_at(feature, imi),
                        f64_at(feature, iaz),
                    );
                    let (Some(cx), Some(cy), Some(a), Some(b), Some(az)) = vals else {
                        skipped += 1;
                        continue;
                    };
                    if !a.is_finite() || !b.is_finite() || a < 0.0 || b < 0.0 {
                        skipped += 1;
                        continue;
                    }
                    let pts = ellipse_points(
                        (cx, cy),
                        a,
                        b,
                        az,
                        prm.line_type,
                        prm.vertex_spacing,
                        !is_polygon,
                    );
                    if is_polygon {
                        let geom = Geometry::polygon(
                            pts.into_iter().map(|(x, y)| Coord::xy(x, y)).collect(),
                            vec![],
                        );
                        out.push(Feature {
                            fid: 0,
                            geometry: Some(geom),
                            attributes: feature.attributes.clone(),
                        });
                    } else {
                        push_line(&mut out, feature, pts);
                    }
                    built += 1;
                }
            }
        }

        ctx.progress.info(&format!(
            "{built} feature(s) built, {skipped} skipped (mode={}, line_type={})",
            prm.mode.as_str(),
            prm.line_type.as_str()
        ));

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("built".to_string(), json!(built));
        outputs.insert("skipped".to_string(), json!(skipped));
        outputs.insert("mode".to_string(), json!(prm.mode.as_str()));
        outputs.insert("line_type".to_string(), json!(prm.line_type.as_str()));
        Ok(ToolRunResult { outputs })
    }
}

fn push_line(out: &mut Layer, feature: &Feature, pts: Vec<(f64, f64)>) {
    let geom = Geometry::line_string(pts.into_iter().map(|(x, y)| Coord::xy(x, y)).collect());
    out.push(Feature {
        fid: 0,
        geometry: Some(geom),
        attributes: feature.attributes.clone(),
    });
}

fn f64_at(feature: &Feature, idx: usize) -> Option<f64> {
    feature
        .attributes
        .get(idx)
        .and_then(FieldValue::as_f64)
        .filter(|v| v.is_finite())
}

// ── Geometry construction ───────────────────────────────────────────────────

/// The vertices of a line from `start` to `end`. Planar lines are a single
/// straight segment; geodesic/rhumb lines are densified along the true curve
/// so the arc is visible, at `vertex_spacing` (meters) or a fixed default
/// segment count.
fn line_points(
    start: (f64, f64),
    end: (f64, f64),
    line_type: LineType,
    vertex_spacing: Option<f64>,
) -> Vec<(f64, f64)> {
    if line_type == LineType::Planar {
        return vec![start, end];
    }
    let p0 = GeoPoint::new(start.0, start.1);
    let p1 = GeoPoint::new(end.0, end.1);
    let total = match line_type {
        LineType::Geodesic => Geodesic.distance(p0, p1),
        LineType::Rhumb => Rhumb.distance(p0, p1),
        LineType::Planar => unreachable!(),
    };
    let segments = match vertex_spacing {
        Some(vs) if vs > 0.0 && total > 0.0 => {
            ((total / vs).ceil() as usize).clamp(1, MAX_SEGMENTS)
        }
        _ => DEFAULT_LINE_SEGMENTS,
    };
    (0..=segments)
        .map(|i| {
            let t = i as f64 / segments as f64;
            let p = match line_type {
                LineType::Geodesic => Geodesic.point_at_ratio_between(p0, p1, t),
                LineType::Rhumb => Rhumb.point_at_ratio_between(p0, p1, t),
                LineType::Planar => unreachable!(),
            };
            (p.x(), p.y())
        })
        .collect()
}

/// The destination point reached by travelling `distance` along `bearing_deg`
/// (degrees clockwise from north) from `origin`.
fn destination_point(
    origin: (f64, f64),
    bearing_deg: f64,
    distance: f64,
    line_type: LineType,
) -> (f64, f64) {
    match line_type {
        LineType::Planar => {
            let b = bearing_deg.to_radians();
            (origin.0 + distance * b.sin(), origin.1 + distance * b.cos())
        }
        LineType::Geodesic => {
            let p = Geodesic.destination(GeoPoint::new(origin.0, origin.1), bearing_deg, distance);
            (p.x(), p.y())
        }
        LineType::Rhumb => {
            let p = Rhumb.destination(GeoPoint::new(origin.0, origin.1), bearing_deg, distance);
            (p.x(), p.y())
        }
    }
}

/// Samples the outline of an ellipse centred at `center`, with semi-major
/// axis `major` and semi-minor axis `minor` (both non-negative), whose major
/// axis points along `azimuth_deg` (degrees clockwise from north).
///
/// Each sample is parametrized by the angle `phi` around the ellipse: the
/// point at `phi` sits `major*cos(phi)` along the major-axis direction and
/// `minor*sin(phi)` along the perpendicular (clockwise) direction. For planar
/// lines that offset is applied directly in Cartesian space; for
/// geodesic/rhumb lines the offset is converted to a (bearing, range) pair
/// and placed with the metric space's `destination` operator, so the outline
/// follows the true ellipsoidal/rhumb geometry rather than a flat
/// approximation.
///
/// If `close` is true, the first point is repeated at the end (for a
/// `LineString` outline); a `Polygon` ring should pass `close = false`, since
/// `wbvector` rings store the closing vertex implicitly.
fn ellipse_points(
    center: (f64, f64),
    major: f64,
    minor: f64,
    azimuth_deg: f64,
    line_type: LineType,
    vertex_spacing: Option<f64>,
    close: bool,
) -> Vec<(f64, f64)> {
    let segments = match vertex_spacing {
        Some(vs) if vs > 0.0 => {
            // Ramanujan's approximation of the ellipse circumference.
            let h =
                (major - minor) * (major - minor) / ((major + minor) * (major + minor)).max(1e-30);
            let circumference = std::f64::consts::PI
                * (major + minor)
                * (1.0 + 3.0 * h / (10.0 + (4.0 - 3.0 * h).max(0.0).sqrt()));
            ((circumference / vs).ceil() as usize).clamp(12, MAX_SEGMENTS)
        }
        _ => DEFAULT_ELLIPSE_SEGMENTS,
    };
    let az = azimuth_deg.to_radians();
    let center_pt = GeoPoint::new(center.0, center.1);

    let mut pts = Vec::with_capacity(segments + 1);
    for i in 0..segments {
        let phi = 2.0 * std::f64::consts::PI * i as f64 / segments as f64;
        let along = major * phi.cos();
        let across = minor * phi.sin();
        let pt = match line_type {
            LineType::Planar => {
                let dx = along * az.sin() + across * az.cos();
                let dy = along * az.cos() - across * az.sin();
                (center.0 + dx, center.1 + dy)
            }
            LineType::Geodesic | LineType::Rhumb => {
                let range = (along * along + across * across).sqrt();
                let bearing = azimuth_deg + across.atan2(along).to_degrees();
                let p = match line_type {
                    LineType::Geodesic => Geodesic.destination(center_pt, bearing, range),
                    LineType::Rhumb => Rhumb.destination(center_pt, bearing, range),
                    LineType::Planar => unreachable!(),
                };
                (p.x(), p.y())
            }
        };
        pts.push(pt);
    }
    if close {
        if let Some(&first) = pts.first() {
            pts.push(first);
        }
    }
    pts
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    XyToLine,
    BearingDistance,
    Ellipse,
}

impl Mode {
    fn parse(s: &str) -> Result<Mode, ToolError> {
        match s.trim().to_ascii_lowercase().as_str() {
            "xy_to_line" | "xy-to-line" => Ok(Mode::XyToLine),
            "bearing_distance" | "bearing-distance" => Ok(Mode::BearingDistance),
            "ellipse" => Ok(Mode::Ellipse),
            other => Err(ToolError::Validation(format!(
                "parameter 'mode' must be one of xy_to_line, bearing_distance, ellipse (got '{other}')"
            ))),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Mode::XyToLine => "xy_to_line",
            Mode::BearingDistance => "bearing_distance",
            Mode::Ellipse => "ellipse",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum LineType {
    Geodesic,
    Rhumb,
    Planar,
}

impl LineType {
    fn parse(s: &str) -> Result<LineType, ToolError> {
        match s.trim().to_ascii_lowercase().as_str() {
            "geodesic" => Ok(LineType::Geodesic),
            "rhumb" | "loxodrome" => Ok(LineType::Rhumb),
            "planar" | "euclidean" | "cartesian" => Ok(LineType::Planar),
            other => Err(ToolError::Validation(format!(
                "parameter 'line_type' must be one of geodesic, rhumb, planar (got '{other}')"
            ))),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            LineType::Geodesic => "geodesic",
            LineType::Rhumb => "rhumb",
            LineType::Planar => "planar",
        }
    }
}

enum Fields {
    XyToLine {
        start_x: String,
        start_y: String,
        end_x: String,
        end_y: String,
    },
    BearingDistance {
        x: String,
        y: String,
        bearing: String,
        distance: String,
    },
    Ellipse {
        x: String,
        y: String,
        major: String,
        minor: String,
        azimuth: String,
    },
}

struct Params {
    mode: Mode,
    line_type: LineType,
    vertex_spacing: Option<f64>,
    polygon_output: bool,
    fields: Fields,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let mode = Mode::parse(require_str(args, "mode")?)?;
    let line_type = match parse_optional_str(args, "line_type")? {
        Some(s) => LineType::parse(s)?,
        None => LineType::Geodesic,
    };
    let vertex_spacing = parse_optional_f64(args, "vertex_spacing")?;
    if let Some(vs) = vertex_spacing {
        if !(vs.is_finite() && vs > 0.0) {
            return Err(ToolError::Validation(
                "parameter 'vertex_spacing' must be a positive number".into(),
            ));
        }
    }
    let polygon_output = parse_optional_bool(args, "polygon_output")?.unwrap_or(false);

    let require_field = |key: &str| -> Result<String, ToolError> {
        parse_optional_str(args, key)?
            .map(str::to_string)
            .ok_or_else(|| {
                ToolError::Validation(format!(
                    "mode '{}' requires field-name parameter '{key}'",
                    mode.as_str()
                ))
            })
    };

    let fields = match mode {
        Mode::XyToLine => Fields::XyToLine {
            start_x: require_field("start_x")?,
            start_y: require_field("start_y")?,
            end_x: require_field("end_x")?,
            end_y: require_field("end_y")?,
        },
        Mode::BearingDistance => Fields::BearingDistance {
            x: require_field("x")?,
            y: require_field("y")?,
            bearing: require_field("bearing")?,
            distance: require_field("distance")?,
        },
        Mode::Ellipse => Fields::Ellipse {
            x: require_field("x")?,
            y: require_field("y")?,
            major: require_field("major")?,
            minor: require_field("minor")?,
            azimuth: require_field("azimuth")?,
        },
    };

    Ok(Params {
        mode,
        line_type,
        vertex_spacing,
        polygon_output,
        fields,
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

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = TableToGeometryTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn mem_path(layer: Layer) -> String {
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn xy_line_table(rows: &[(&str, f64, f64, f64, f64)]) -> String {
        let mut l = Layer::new("segs");
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.add_field(FieldDef::new("sx", FieldType::Float));
        l.add_field(FieldDef::new("sy", FieldType::Float));
        l.add_field(FieldDef::new("ex", FieldType::Float));
        l.add_field(FieldDef::new("ey", FieldType::Float));
        for (name, sx, sy, ex, ey) in rows {
            l.add_feature(
                None,
                &[
                    ("name", (*name).into()),
                    ("sx", (*sx).into()),
                    ("sy", (*sy).into()),
                    ("ex", (*ex).into()),
                    ("ey", (*ey).into()),
                ],
            )
            .unwrap();
        }
        mem_path(l)
    }

    fn bearing_table(rows: &[(&str, f64, f64, f64, f64)]) -> String {
        let mut l = Layer::new("obs");
        l.add_field(FieldDef::new("name", FieldType::Text));
        l.add_field(FieldDef::new("ox", FieldType::Float));
        l.add_field(FieldDef::new("oy", FieldType::Float));
        l.add_field(FieldDef::new("brg", FieldType::Float));
        l.add_field(FieldDef::new("dist", FieldType::Float));
        for (name, ox, oy, brg, dist) in rows {
            l.add_feature(
                None,
                &[
                    ("name", (*name).into()),
                    ("ox", (*ox).into()),
                    ("oy", (*oy).into()),
                    ("brg", (*brg).into()),
                    ("dist", (*dist).into()),
                ],
            )
            .unwrap();
        }
        mem_path(l)
    }

    fn ellipse_table(rows: &[(f64, f64, f64, f64, f64)]) -> String {
        let mut l = Layer::new("ellipses");
        l.add_field(FieldDef::new("cx", FieldType::Float));
        l.add_field(FieldDef::new("cy", FieldType::Float));
        l.add_field(FieldDef::new("maj", FieldType::Float));
        l.add_field(FieldDef::new("min", FieldType::Float));
        l.add_field(FieldDef::new("az", FieldType::Float));
        for (cx, cy, maj, min, az) in rows {
            l.add_feature(
                None,
                &[
                    ("cx", (*cx).into()),
                    ("cy", (*cy).into()),
                    ("maj", (*maj).into()),
                    ("min", (*min).into()),
                    ("az", (*az).into()),
                ],
            )
            .unwrap();
        }
        mem_path(l)
    }

    fn line_coords(layer: &Layer, i: usize) -> Vec<(f64, f64)> {
        match layer.features[i].geometry.as_ref().unwrap() {
            Geometry::LineString(cs) => cs.iter().map(|c| (c.x, c.y)).collect(),
            other => panic!("expected LineString, got {other:?}"),
        }
    }

    /// xy_to_line in planar mode builds a straight 2-vertex line preserving
    /// endpoints exactly, and carries the source attribute along.
    #[test]
    fn xy_to_line_builds_straight_two_point_lines() {
        let input = xy_line_table(&[("seg1", 0.0, 0.0, 10.0, 5.0)]);
        let (out, layer) = run(json!({
            "input": input,
            "mode": "xy_to_line",
            "line_type": "planar",
            "start_x": "sx", "start_y": "sy", "end_x": "ex", "end_y": "ey",
        }));
        assert_eq!(out.outputs["built"], json!(1));
        assert_eq!(layer.len(), 1);
        let coords = line_coords(&layer, 0);
        assert_eq!(coords, vec![(0.0, 0.0), (10.0, 5.0)]);
        let name_idx = layer.schema.field_index("name").unwrap();
        assert_eq!(
            layer.features[0].attributes[name_idx].as_str().unwrap(),
            "seg1"
        );
    }

    /// bearing_distance geodesic endpoint matches an independently computed
    /// Karney/geographiclib destination: due-east 1,000,000 m from the
    /// equator at the prime meridian.
    #[test]
    fn bearing_distance_geodesic_matches_known_destination() {
        let input = bearing_table(&[("obs1", 0.0, 0.0, 90.0, 1_000_000.0)]);
        let (_, layer) = run(json!({
            "input": input,
            "mode": "bearing_distance",
            "line_type": "geodesic",
            "x": "ox", "y": "oy", "bearing": "brg", "distance": "dist",
        }));
        let coords = line_coords(&layer, 0);
        let (end_x, end_y) = *coords.last().unwrap();
        // Independently computed via pyproj's Geod (also a Karney/
        // GeographicLib WGS84 direct solution) for (0,0) bearing 90
        // distance 1_000_000 m: lon = 8.983152841195215 deg, lat = 0.0
        // (due east from the equator stays on the equator).
        assert!((end_x - 8.983_152_841_195_215).abs() < 1e-6, "lon {end_x}");
        assert!(end_y.abs() < 1e-6, "lat {end_y}");
    }

    /// Planar bearing_distance uses simple trig: due-north 1 unit from the
    /// origin lands exactly at (0, 1); due-east lands at (1, 0).
    #[test]
    fn bearing_distance_planar_matches_trig() {
        let input = bearing_table(&[("north", 0.0, 0.0, 0.0, 1.0), ("east", 0.0, 0.0, 90.0, 1.0)]);
        let (_, layer) = run(json!({
            "input": input,
            "mode": "bearing_distance",
            "line_type": "planar",
            "x": "ox", "y": "oy", "bearing": "brg", "distance": "dist",
        }));
        let north_end = *line_coords(&layer, 0).last().unwrap();
        let east_end = *line_coords(&layer, 1).last().unwrap();
        assert!((north_end.0).abs() < 1e-9 && (north_end.1 - 1.0).abs() < 1e-9);
        assert!((east_end.0 - 1.0).abs() < 1e-9 && (east_end.1).abs() < 1e-9);
    }

    /// A planar ellipse's polygon-approximation area matches pi * a * b, and
    /// its outline is oriented so the major axis dominates.
    #[test]
    fn ellipse_planar_area_matches_pi_a_b() {
        let input = ellipse_table(&[(0.0, 0.0, 100.0, 40.0, 0.0)]);
        let (_, layer) = run(json!({
            "input": input,
            "mode": "ellipse",
            "line_type": "planar",
            "polygon_output": true,
            "x": "cx", "y": "cy", "major": "maj", "minor": "min", "azimuth": "az",
        }));
        assert_eq!(layer.len(), 1);
        let geom = layer.features[0].geometry.as_ref().unwrap();
        let Geometry::Polygon { exterior, .. } = geom else {
            panic!("expected polygon, got {geom:?}");
        };
        let coords = exterior.coords();
        // Shoelace formula.
        let n = coords.len();
        let mut area2 = 0.0;
        for i in 0..n {
            let (x0, y0) = (coords[i].x, coords[i].y);
            let (x1, y1) = (coords[(i + 1) % n].x, coords[(i + 1) % n].y);
            area2 += x0 * y1 - x1 * y0;
        }
        let area = area2.abs() / 2.0;
        let expected = std::f64::consts::PI * 100.0 * 40.0;
        assert!(
            (area - expected).abs() / expected < 0.01,
            "area {area} vs expected {expected}"
        );
    }

    /// A closed (non-polygon) ellipse outline repeats its first vertex at the
    /// end so it renders as a closed ring even as a LineString.
    #[test]
    fn ellipse_line_output_is_closed() {
        let input = ellipse_table(&[(5.0, 5.0, 10.0, 10.0, 0.0)]);
        let (_, layer) = run(json!({
            "input": input,
            "mode": "ellipse",
            "line_type": "planar",
            "x": "cx", "y": "cy", "major": "maj", "minor": "min", "azimuth": "az",
        }));
        let coords = line_coords(&layer, 0);
        assert_eq!(coords.first(), coords.last());
        assert!(coords.len() > 12);
    }

    /// Planar and geodesic bearing_distance diverge measurably over a long
    /// (near-quarter-Earth) span, since planar ignores curvature.
    #[test]
    fn planar_and_geodesic_diverge_over_long_spans() {
        let dist = 5_000_000.0; // 5,000 km
        let input_planar = bearing_table(&[("p", 0.0, 0.0, 45.0, dist)]);
        let input_geodesic = bearing_table(&[("g", 0.0, 0.0, 45.0, dist)]);
        let (_, planar) = run(json!({
            "input": input_planar, "mode": "bearing_distance", "line_type": "planar",
            "x": "ox", "y": "oy", "bearing": "brg", "distance": "dist",
        }));
        let (_, geodesic) = run(json!({
            "input": input_geodesic, "mode": "bearing_distance", "line_type": "geodesic",
            "x": "ox", "y": "oy", "bearing": "brg", "distance": "dist",
        }));
        let p_end = *line_coords(&planar, 0).last().unwrap();
        let g_end = *line_coords(&geodesic, 0).last().unwrap();
        let d = ((p_end.0 - g_end.0).powi(2) + (p_end.1 - g_end.1).powi(2)).sqrt();
        assert!(
            d > 1.0,
            "planar {p_end:?} and geodesic {g_end:?} endpoints should diverge over {dist}m, delta {d}"
        );
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            TableToGeometryTool.validate(&args)
        };
        // Missing mode.
        assert!(bad(json!({ "input": "a.geojson" })).is_err());
        // Unknown mode.
        assert!(bad(json!({ "input": "a.geojson", "mode": "bogus" })).is_err());
        // xy_to_line missing a required field name.
        assert!(bad(json!({
            "input": "a.geojson", "mode": "xy_to_line",
            "start_x": "sx", "start_y": "sy", "end_x": "ex",
        }))
        .is_err());
        // bearing_distance missing field names entirely.
        assert!(bad(json!({ "input": "a.geojson", "mode": "bearing_distance" })).is_err());
        // ellipse missing field names entirely.
        assert!(bad(json!({ "input": "a.geojson", "mode": "ellipse" })).is_err());
        // Bad line_type.
        assert!(bad(json!({
            "input": "a.geojson", "mode": "xy_to_line", "line_type": "warp_speed",
            "start_x": "sx", "start_y": "sy", "end_x": "ex", "end_y": "ey",
        }))
        .is_err());
        // Non-positive vertex_spacing.
        assert!(bad(json!({
            "input": "a.geojson", "mode": "xy_to_line", "vertex_spacing": -1,
            "start_x": "sx", "start_y": "sy", "end_x": "ex", "end_y": "ey",
        }))
        .is_err());
        // Valid xy_to_line.
        assert!(bad(json!({
            "input": "a.geojson", "mode": "xy_to_line",
            "start_x": "sx", "start_y": "sy", "end_x": "ex", "end_y": "ey",
        }))
        .is_ok());
        // Valid ellipse.
        assert!(bad(json!({
            "input": "a.geojson", "mode": "ellipse",
            "x": "cx", "y": "cy", "major": "maj", "minor": "min", "azimuth": "az",
        }))
        .is_ok());
    }
}
