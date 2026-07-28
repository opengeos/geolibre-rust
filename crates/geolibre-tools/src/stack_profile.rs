//! GeoLibre tool: profile a line across several surfaces at once.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Stack Profile* (3D Analyst). Two
//! nearby tools exist and neither does this: the bundled `profile` samples
//! **one** surface along a line, and `image_stack_profile` samples a raster
//! **stack at points**, not along a line.
//!
//! The common request — bare-earth vs first-return DSM vs a proposed design
//! surface, plotted as one cross-section — needs a line sampled against N
//! surfaces on a *shared* distance axis. Running `profile` N times gives N
//! tables on N independently-sampled axes, which cannot then be joined without
//! resampling. Doing it in one pass is the whole point: the line is densified
//! once, and every surface is read at exactly the same stations.
//!
//! Output is long-form (one row per station per surface), which is the shape a
//! plotting library wants and which keeps the schema fixed no matter how many
//! surfaces are supplied. Stations landing on NoData in a given surface are
//! emitted with a null Z rather than dropped, so the distance axis stays
//! aligned across surfaces.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::common::load_input_raster;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Cap on emitted rows so a tiny `sample_distance` fails loudly.
const MAX_ROWS: usize = 20_000_000;

pub struct StackProfileTool;

impl Tool for StackProfileTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "stack_profile",
            display_name: "Stack Profile",
            summary: "Sample several surface rasters along the same polyline and emit one long-form table (LINE_ID, FIRST_DIST, FIRST_Z, SRC_ID, SRC_NAME) so the surfaces can be compared in cross-section on a shared distance axis. Like ArcGIS Stack Profile.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input line features to profile.",
                    required: true,
                },
                ToolParamSpec {
                    name: "surfaces",
                    description: "Comma/semicolon-separated list of surface raster paths, sampled in the order given.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output profile table path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "sample_distance",
                    description: "Station spacing in map units (default: the finest input cell size).",
                    required: false,
                },
                ToolParamSpec {
                    name: "line_id_field",
                    description: "Optional field supplying LINE_ID; defaults to the feature index.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "Surface sampling: 'bilinear' (default) or 'nearest'.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        if split_list(require_str(args, "surfaces")?).is_empty() {
            return Err(ToolError::Validation(
                "'surfaces' must list at least one raster".to_string(),
            ));
        }
        parse_method(args)?;
        if let Some(d) = parse_optional_f64(args, "sample_distance")? {
            if d <= 0.0 {
                return Err(ToolError::Validation(
                    "'sample_distance' must be positive".to_string(),
                ));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let surface_paths = split_list(require_str(args, "surfaces")?);
        let output = parse_optional_str(args, "output")?;
        let method = parse_method(args)?;
        let id_field = parse_optional_str(args, "line_id_field")?;

        let layer = load_input_layer(input)?;
        if layer.features.is_empty() {
            return Err(ToolError::Execution("input has no features".to_string()));
        }
        let id_idx = match id_field {
            Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                ToolError::Validation(format!("line_id_field '{f}' not found"))
            })?),
            None => None,
        };

        ctx.progress
            .info(&format!("loading {} surface(s)", surface_paths.len()));
        let mut surfaces: Vec<(String, Raster)> = Vec::with_capacity(surface_paths.len());
        for p in &surface_paths {
            let r = load_input_raster(p)?;
            surfaces.push((surface_name(p), r));
        }

        // Warn (do not fail) on mixed CRS: the station axis is line-driven, so
        // a mismatch produces NoData rather than silently wrong elevations, but
        // the user should know why their column came back empty.
        let epsgs: Vec<Option<u32>> = surfaces.iter().map(|(_, r)| raster_epsg(r)).collect();
        if epsgs.iter().flatten().collect::<std::collections::BTreeSet<_>>().len() > 1 {
            ctx.progress
                .info("surfaces do not share a CRS; stations outside a surface will be null");
        }

        let step = match parse_optional_f64(args, "sample_distance")? {
            Some(d) => d,
            None => surfaces
                .iter()
                .map(|(_, r)| r.cell_size_x.min(r.cell_size_y))
                .fold(f64::INFINITY, f64::min),
        };
        if !(step > 0.0 && step.is_finite()) {
            return Err(ToolError::Execution(
                "could not determine a positive sample distance from the surfaces".to_string(),
            ));
        }

        let mut out = Layer::new("stack_profile");
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("LINE_ID", FieldType::Text));
        out.add_field(FieldDef::new("FIRST_DIST", FieldType::Float));
        out.add_field(FieldDef::new("FIRST_Z", FieldType::Float));
        out.add_field(FieldDef::new("SRC_ID", FieldType::Integer));
        out.add_field(FieldDef::new("SRC_NAME", FieldType::Text));

        let mut rows = 0usize;
        let mut stations_total = 0usize;
        let mut null_samples = 0usize;
        let mut per_surface_valid = vec![0usize; surfaces.len()];

        for (fid, feat) in layer.iter().enumerate() {
            let Some(g) = &feat.geometry else { continue };
            let line_id = match id_idx {
                Some(i) => key_of(feat.attributes.get(i)),
                None => fid.to_string(),
            };
            for path in line_paths(g) {
                let stations = densify(&path, step);
                if stations.is_empty() {
                    continue;
                }
                stations_total += stations.len();
                if rows + stations.len() * surfaces.len() > MAX_ROWS {
                    return Err(ToolError::Execution(format!(
                        "sample_distance {step} would emit more than {MAX_ROWS} rows; \
                         increase 'sample_distance'"
                    )));
                }
                // Sample every surface at the SAME station list — the point of
                // the tool.
                for (si, (name, raster)) in surfaces.iter().enumerate() {
                    for &(x, y, dist) in &stations {
                        let z = sample(raster, 0, x, y, method);
                        if z.is_some() {
                            per_surface_valid[si] += 1;
                        } else {
                            null_samples += 1;
                        }
                        out.add_feature(
                            None,
                            &[
                                ("LINE_ID", FieldValue::Text(line_id.clone())),
                                ("FIRST_DIST", FieldValue::Float(dist)),
                                ("FIRST_Z", z.map_or(FieldValue::Null, FieldValue::Float)),
                                ("SRC_ID", FieldValue::Integer(si as i64)),
                                ("SRC_NAME", FieldValue::Text(name.clone())),
                            ],
                        )
                        .map_err(|e| ToolError::Execution(format!("failed adding row: {e}")))?;
                        rows += 1;
                    }
                }
            }
        }
        if rows == 0 {
            return Err(ToolError::Execution(
                "no line geometry produced any profile stations".to_string(),
            ));
        }
        ctx.progress
            .info(&format!("profiled {stations_total} station(s) x {} surface(s)", surfaces.len()));

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("row_count".to_string(), json!(rows));
        outputs.insert("station_count".to_string(), json!(stations_total));
        outputs.insert("surface_count".to_string(), json!(surfaces.len()));
        outputs.insert("sample_distance".to_string(), json!(step));
        outputs.insert("null_sample_count".to_string(), json!(null_samples));
        outputs.insert(
            "surface_names".to_string(),
            json!(surfaces.iter().map(|(n, _)| n.clone()).collect::<Vec<_>>()),
        );
        outputs.insert("valid_per_surface".to_string(), json!(per_surface_valid));
        Ok(ToolRunResult { outputs })
    }
}

/// Densifies a path into `(x, y, cumulative_distance)` stations every `step`
/// units, always including the start and end vertices so the profile spans the
/// whole line.
fn densify(path: &[Coord], step: f64) -> Vec<(f64, f64, f64)> {
    let mut out = Vec::new();
    if path.len() < 2 {
        return out;
    }
    let total: f64 = path
        .windows(2)
        .map(|w| (w[1].x - w[0].x).hypot(w[1].y - w[0].y))
        .sum();
    if total <= 0.0 {
        return out;
    }
    out.push((path[0].x, path[0].y, 0.0));

    let mut acc = 0.0_f64;
    let mut next = step;
    for w in path.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        let seg = (b.x - a.x).hypot(b.y - a.y);
        if seg <= 0.0 {
            continue;
        }
        while next <= acc + seg + 1e-9 {
            let t = ((next - acc) / seg).clamp(0.0, 1.0);
            // Skip a station that coincides with the end vertex appended below.
            if (next - total).abs() > 1e-9 {
                out.push((a.x + t * (b.x - a.x), a.y + t * (b.y - a.y), next));
            }
            next += step;
        }
        acc += seg;
    }
    let last = &path[path.len() - 1];
    out.push((last.x, last.y, total));
    out
}

fn line_paths(g: &Geometry) -> Vec<Vec<Coord>> {
    match g {
        Geometry::LineString(cs) if cs.len() >= 2 => vec![cs.clone()],
        Geometry::MultiLineString(ls) => ls.iter().filter(|l| l.len() >= 2).cloned().collect(),
        // Polygon boundaries profile fine and are a legitimate cross-section.
        Geometry::Polygon { exterior, .. } if exterior.0.len() >= 2 => vec![exterior.0.clone()],
        Geometry::GeometryCollection(gs) => gs.iter().flat_map(line_paths).collect(),
        _ => Vec::new(),
    }
}

fn surface_name(path: &str) -> String {
    std::path::Path::new(path)
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| path.to_string())
}

fn raster_epsg(r: &Raster) -> Option<u32> {
    r.crs.epsg
}

// ── Surface sampling ────────────────────────────────────────────────────────

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
    // Partial coverage at the raster edge: fall back to the nearest valid
    // contributor rather than emitting null for the whole station.
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

// ── Params ──────────────────────────────────────────────────────────────────

fn parse_method(args: &ToolArgs) -> Result<Method, ToolError> {
    match args
        .get("method")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        None | Some("") | Some("bilinear") => Ok(Method::Bilinear),
        Some("nearest") => Ok(Method::Nearest),
        Some(o) => Err(ToolError::Validation(format!(
            "'method' must be 'bilinear' or 'nearest', got '{o}'"
        ))),
    }
}

fn key_of(v: Option<&FieldValue>) -> String {
    match v {
        None | Some(FieldValue::Null) => "NULL".to_string(),
        Some(FieldValue::Integer(i)) => i.to_string(),
        Some(FieldValue::Float(f)) => format!("{f}"),
        Some(FieldValue::Text(s)) => s.clone(),
        Some(FieldValue::Boolean(b)) => b.to_string(),
        Some(FieldValue::Date(s)) | Some(FieldValue::DateTime(s)) => s.clone(),
        Some(FieldValue::Blob(b)) => format!("blob[{}]", b.len()),
    }
}

fn split_list(s: &str) -> Vec<String> {
    s.split([',', ';'])
        .map(str::trim)
        .filter(|p| !p.is_empty())
        .map(str::to_string)
        .collect()
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
    use wbvector::GeometryType;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// 10x10 raster over (0,0)-(10,10), cell 1, filled by `f(x, y)`.
    fn raster(f: impl Fn(f64, f64) -> f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols: 10,
            rows: 10,
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
        for row in 0..10 {
            for col in 0..10 {
                let x = 0.5 + col as f64;
                let y = 9.5 - row as f64;
                r.set(0, row as isize, col as isize, f(x, y)).unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn line() -> String {
        let mut l = Layer::new("l")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_feature(
            Some(Geometry::line_string(vec![
                Coord::xy(0.5, 5.0),
                Coord::xy(9.5, 5.0),
            ])),
            &[],
        )
        .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = StackProfileTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    #[test]
    fn all_surfaces_share_one_distance_axis() {
        // Two surfaces: ground = x, and a DSM 10 above it.
        let ground = raster(|x, _y| x);
        let dsm = raster(|x, _y| x + 10.0);
        let (out, layer) = run(json!({
            "input": line(), "surfaces": format!("{ground},{dsm}"), "sample_distance": 1.0
        }));
        assert_eq!(out.outputs["surface_count"], json!(2));
        let stations = out.outputs["station_count"].as_f64().unwrap() as usize;
        assert_eq!(out.outputs["row_count"].as_f64().unwrap() as usize, stations * 2);

        // The two surfaces must report identical FIRST_DIST sequences — the
        // property that makes the output joinable/plottable.
        let (d, s) = (
            layer.schema.field_index("FIRST_DIST").unwrap(),
            layer.schema.field_index("SRC_ID").unwrap(),
        );
        let axis = |src: i64| -> Vec<f64> {
            layer
                .iter()
                .filter(|f| f.attributes[s].as_i64() == Some(src))
                .map(|f| f.attributes[d].as_f64().unwrap())
                .collect()
        };
        assert_eq!(axis(0), axis(1));
    }

    #[test]
    fn z_values_track_each_surface() {
        let ground = raster(|x, _y| x);
        let dsm = raster(|x, _y| x + 10.0);
        let (_o, layer) = run(json!({
            "input": line(), "surfaces": format!("{ground},{dsm}"), "sample_distance": 3.0
        }));
        let (z, s) = (
            layer.schema.field_index("FIRST_Z").unwrap(),
            layer.schema.field_index("SRC_ID").unwrap(),
        );
        let zs = |src: i64| -> Vec<f64> {
            layer
                .iter()
                .filter(|f| f.attributes[s].as_i64() == Some(src))
                .map(|f| f.attributes[z].as_f64().unwrap())
                .collect()
        };
        let (g, d) = (zs(0), zs(1));
        assert_eq!(g.len(), d.len());
        for (a, b) in g.iter().zip(d.iter()) {
            assert!((b - a - 10.0).abs() < 1e-6, "offset was {}", b - a);
        }
    }

    #[test]
    fn surface_names_come_from_the_paths() {
        let a = raster(|x, _| x);
        let (out, _l) = run(json!({
            "input": line(), "surfaces": a, "sample_distance": 5.0
        }));
        let names = out.outputs["surface_names"].as_array().unwrap();
        assert_eq!(names.len(), 1);
    }

    #[test]
    fn stations_include_both_endpoints() {
        let a = raster(|x, _| x);
        let (_o, layer) = run(json!({
            "input": line(), "surfaces": a, "sample_distance": 2.0
        }));
        let d = layer.schema.field_index("FIRST_DIST").unwrap();
        let dists: Vec<f64> = layer
            .iter()
            .map(|f| f.attributes[d].as_f64().unwrap())
            .collect();
        assert_eq!(dists.first(), Some(&0.0));
        // The line runs 0.5 -> 9.5, i.e. 9.0 long.
        assert!((dists.last().unwrap() - 9.0).abs() < 1e-9);
    }

    #[test]
    fn out_of_range_stations_are_null_not_dropped() {
        // A line running well outside the raster still yields rows, with nulls.
        let mut l = Layer::new("l")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_feature(
            Some(Geometry::line_string(vec![
                Coord::xy(100.0, 100.0),
                Coord::xy(110.0, 100.0),
            ])),
            &[],
        )
        .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let far = wbvector::memory_store::make_vector_memory_path(&id);
        let a = raster(|x, _| x);
        let (out, layer) = run(json!({
            "input": far, "surfaces": a, "sample_distance": 5.0
        }));
        assert!(out.outputs["row_count"].as_f64().unwrap() > 0.0);
        let z = layer.schema.field_index("FIRST_Z").unwrap();
        assert!(layer.iter().all(|f| f.attributes[z].is_null()));
        assert_eq!(out.outputs["valid_per_surface"], json!([0]));
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            StackProfileTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "l.shp" })).is_err());
        assert!(bad(json!({ "input": "l.shp", "surfaces": "a.tif", "sample_distance": 0 })).is_err());
        assert!(bad(json!({ "input": "l.shp", "surfaces": "a.tif", "method": "cubic" })).is_err());
        assert!(bad(json!({ "input": "l.shp", "surfaces": "a.tif,b.tif" })).is_ok());
    }
}
