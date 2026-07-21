//! GeoLibre tool: reservoir stage-area-volume curves from a DEM.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Storage Capacity* (Spatial Analyst).
//! The bundled `impoundment_size_index` evaluates dam *siting*; nothing produces
//! the elevation → surface-area → storage-volume table a given basin needs for
//! reservoir and detention design. Strong hydrology-identity fit alongside
//! `cut_fill` and the depression/sink tools.
//!
//! For each analysis zone (a polygon from `zones`, or the whole DEM when none is
//! given), a series of water-surface elevations is swept from the zone minimum
//! to its maximum (`num_levels` steps, or a fixed `increment`). At each level
//! `L` the tool accumulates, over cells whose terrain `z <= L`:
//!
//! * `area  = Σ cell_area`
//! * `volume = Σ (L − z) · cell_area`
//!
//! reusing `cut_fill`'s volume accounting. The output is a CSV of
//! `zone_id, level, elevation, area, volume`; a per-zone minimum/maximum
//! elevation is reported for reference.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;
use wbvector::{Geometry, Ring};

use crate::common::{load_input_raster, write_text_output};
use crate::vector_common::{load_input_layer, parse_optional_str};

const DEFAULT_LEVELS: usize = 20;

pub struct StorageCapacityTool;

impl Tool for StorageCapacityTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "storage_capacity",
            display_name: "Storage Capacity",
            summary: "Sweep water-surface elevations over a DEM and report the flooded surface area and storage volume at each level (the stage-area-volume curve) per analysis zone or over the whole DEM, like ArcGIS Storage Capacity.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "dem",
                    description: "Input DEM raster (elevations).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output CSV path (zone_id, level, elevation, area, volume).",
                    required: true,
                },
                ToolParamSpec {
                    name: "zones",
                    description: "Optional polygon layer; one stage-area-volume curve per polygon. If omitted, the whole DEM is one zone.",
                    required: false,
                },
                ToolParamSpec {
                    name: "zone_id_field",
                    description: "Field naming each zone in the output (default: the feature index).",
                    required: false,
                },
                ToolParamSpec {
                    name: "num_levels",
                    description: "Number of elevation steps between each zone's min and max (default 20). Ignored if 'increment' is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "increment",
                    description: "Fixed elevation step between levels (map units). Overrides 'num_levels'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_elevation",
                    description: "Optional lower elevation bound for the sweep (default: each zone's minimum terrain).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_elevation",
                    description: "Optional upper elevation bound for the sweep (default: each zone's maximum terrain).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based DEM band (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["dem", "output"] {
            if args
                .get(key)
                .and_then(Value::as_str)
                .map(str::trim)
                .unwrap_or("")
                .is_empty()
            {
                return Err(ToolError::Validation(format!(
                    "missing required string parameter '{key}'"
                )));
            }
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let dem_path = args.get("dem").and_then(Value::as_str).unwrap();
        let out_path = args.get("output").and_then(Value::as_str).unwrap();
        let prm = parse_params(args)?;

        let dem = load_input_raster(dem_path)?;
        if prm.band as usize >= dem.bands {
            return Err(ToolError::Validation(format!(
                "band {} out of range (DEM has {} band(s))",
                prm.band + 1,
                dem.bands
            )));
        }
        let cell_area = dem.cell_size_x * dem.cell_size_y;
        let nodata = dem.nodata;

        // Collect the analysis zones as lists of (row, col) cell indices.
        let zones: Vec<Zone> = match &prm.zones_path {
            Some(path) => {
                let layer = load_input_layer(path)?;
                let id_idx = match &prm.zone_id_field {
                    Some(f) => Some(layer.schema.field_index(f).ok_or_else(|| {
                        ToolError::Validation(format!("zone_id_field '{f}' not found"))
                    })?),
                    None => None,
                };
                let mut zs = Vec::new();
                for (fidx, feature) in layer.features.iter().enumerate() {
                    let Some(geom) = feature.geometry.as_ref() else {
                        continue;
                    };
                    let rings = polygon_rings(geom);
                    if rings.is_empty() {
                        continue;
                    }
                    let cells = cells_in_polygon(&dem, &rings, prm.band, nodata);
                    let id = match id_idx {
                        Some(i) => field_key(&feature.attributes[i]),
                        None => fidx.to_string(),
                    };
                    zs.push(Zone { id, cells });
                }
                zs
            }
            None => {
                // Whole DEM is one zone.
                let mut cells = Vec::new();
                for row in 0..dem.rows as isize {
                    for col in 0..dem.cols as isize {
                        let z = dem.get(prm.band, row, col);
                        if z != nodata && z.is_finite() {
                            cells.push((row, col, z));
                        }
                    }
                }
                vec![Zone {
                    id: "0".to_string(),
                    cells,
                }]
            }
        };

        if zones.iter().all(|z| z.cells.is_empty()) {
            return Err(ToolError::Execution(
                "no valid DEM cells fall within the analysis zone(s)".to_string(),
            ));
        }

        ctx.progress
            .info(&format!("sweeping {} zone(s)", zones.len()));

        let mut csv = String::from("zone_id,level,elevation,area,volume\n");
        let mut total_rows = 0usize;
        let mut max_volume = 0.0_f64;
        for zone in &zones {
            if zone.cells.is_empty() {
                continue;
            }
            let zmin = zone.cells.iter().map(|c| c.2).fold(f64::INFINITY, f64::min);
            let zmax = zone
                .cells
                .iter()
                .map(|c| c.2)
                .fold(f64::NEG_INFINITY, f64::max);
            let lo = prm.min_elevation.unwrap_or(zmin);
            let hi = prm.max_elevation.unwrap_or(zmax);
            if hi <= lo {
                // Flat or inverted range: single level at hi.
                let (area, volume) = accumulate(&zone.cells, hi, cell_area);
                csv.push_str(&format!("{},{},{hi},{area},{volume}\n", zone.id, 0));
                total_rows += 1;
                max_volume = max_volume.max(volume);
                continue;
            }

            let levels = level_series(lo, hi, prm.num_levels, prm.increment);
            for (i, level) in levels.iter().enumerate() {
                let (area, volume) = accumulate(&zone.cells, *level, cell_area);
                csv.push_str(&format!("{},{i},{level},{area},{volume}\n", zone.id));
                total_rows += 1;
                max_volume = max_volume.max(volume);
            }
        }

        write_text_output(&csv, out_path)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("zone_count".to_string(), json!(zones.len()));
        outputs.insert("row_count".to_string(), json!(total_rows));
        outputs.insert("max_volume".to_string(), json!(max_volume));
        Ok(ToolRunResult { outputs })
    }
}

struct Zone {
    id: String,
    /// (row, col, elevation) of the zone's valid DEM cells.
    cells: Vec<(isize, isize, f64)>,
}

/// Surface area and storage volume at water-surface elevation `level`.
fn accumulate(cells: &[(isize, isize, f64)], level: f64, cell_area: f64) -> (f64, f64) {
    let mut area = 0.0;
    let mut volume = 0.0;
    for &(_, _, z) in cells {
        if z <= level {
            area += cell_area;
            volume += (level - z) * cell_area;
        }
    }
    (area, volume)
}

/// The elevation levels to sweep: fixed `increment` (inclusive of `hi`) or
/// `num_levels` evenly spaced from `lo` to `hi`.
fn level_series(lo: f64, hi: f64, num_levels: usize, increment: Option<f64>) -> Vec<f64> {
    if let Some(step) = increment {
        let step = step.abs().max(f64::MIN_POSITIVE);
        let mut levels = Vec::new();
        let mut l = lo;
        while l < hi {
            levels.push(l);
            l += step;
        }
        levels.push(hi);
        levels
    } else {
        let n = num_levels.max(1);
        (0..=n)
            .map(|i| lo + (hi - lo) * i as f64 / n as f64)
            .collect()
    }
}

/// Every valid DEM cell whose center falls inside the polygon (bbox-prefiltered).
fn cells_in_polygon(
    dem: &Raster,
    rings: &[Vec<(f64, f64)>],
    band: isize,
    nodata: f64,
) -> Vec<(isize, isize, f64)> {
    // Zone bounding box in world coords.
    let (mut minx, mut miny, mut maxx, mut maxy) = (
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );
    for ring in rings {
        for &(x, y) in ring {
            minx = minx.min(x);
            miny = miny.min(y);
            maxx = maxx.max(x);
            maxy = maxy.max(y);
        }
    }
    // Restrict to the pixel window overlapping the bbox.
    let mut cells = Vec::new();
    for row in 0..dem.rows as isize {
        let cy = dem.row_center_y(row);
        if cy < miny || cy > maxy {
            continue;
        }
        for col in 0..dem.cols as isize {
            let cx = dem.col_center_x(col);
            if cx < minx || cx > maxx {
                continue;
            }
            if !point_in_rings(cx, cy, rings) {
                continue;
            }
            let z = dem.get(band, row, col);
            if z != nodata && z.is_finite() {
                cells.push((row, col, z));
            }
        }
    }
    cells
}

/// Even-odd test across an exterior ring minus holes (rings[0] = exterior).
fn point_in_rings(x: f64, y: f64, rings: &[Vec<(f64, f64)>]) -> bool {
    if rings.is_empty() || !point_in_ring(x, y, &rings[0]) {
        return false;
    }
    !rings[1..].iter().any(|h| point_in_ring(x, y, h))
}

fn point_in_ring(x: f64, y: f64, ring: &[(f64, f64)]) -> bool {
    let n = ring.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = ring[i];
        let (xj, yj) = ring[j];
        if (yi > y) != (yj > y) {
            let xcross = (xj - xi) * (y - yi) / (yj - yi) + xi;
            if x < xcross {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

/// Polygon rings as (x,y) chains (exterior first) without the closing duplicate.
fn polygon_rings(geom: &Geometry) -> Vec<Vec<(f64, f64)>> {
    let ring_pts =
        |ring: &Ring| -> Vec<(f64, f64)> { ring.coords().iter().map(|c| (c.x, c.y)).collect() };
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            let mut out = vec![ring_pts(exterior)];
            out.extend(interiors.iter().map(&ring_pts));
            out
        }
        Geometry::MultiPolygon(parts) => {
            // Treat every part's exterior as its own analysis region within one zone.
            let mut out = Vec::new();
            for (ext, holes) in parts {
                out.push(ring_pts(ext));
                out.extend(holes.iter().map(&ring_pts));
            }
            out
        }
        _ => Vec::new(),
    }
}

fn field_key(fv: &wbvector::FieldValue) -> String {
    if let Some(i) = fv.as_i64() {
        i.to_string()
    } else if let Some(f) = fv.as_f64() {
        format!("{f}")
    } else {
        fv.as_str().unwrap_or("").to_string()
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    zones_path: Option<String>,
    zone_id_field: Option<String>,
    num_levels: usize,
    increment: Option<f64>,
    min_elevation: Option<f64>,
    max_elevation: Option<f64>,
    band: isize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let num_levels = match args.get("num_levels") {
        None | Some(Value::Null) => DEFAULT_LEVELS,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(DEFAULT_LEVELS as u64).max(1) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => DEFAULT_LEVELS,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'num_levels' must be an integer".into()))?
            .max(1),
        Some(_) => {
            return Err(ToolError::Validation(
                "'num_levels' must be a number".into(),
            ))
        }
    };
    let increment = parse_f64(args, "increment")?;
    if let Some(i) = increment {
        if i.is_nan() || i <= 0.0 {
            return Err(ToolError::Validation("'increment' must be positive".into()));
        }
    }
    let band = args.get("band").and_then(Value::as_u64).unwrap_or(1).max(1) as isize - 1;
    Ok(Params {
        zones_path: parse_optional_str(args, "zones")?.map(str::to_string),
        zone_id_field: parse_optional_str(args, "zone_id_field")?.map(str::to_string),
        num_levels,
        increment,
        min_elevation: parse_f64(args, "min_elevation")?,
        max_elevation: parse_f64(args, "max_elevation")?,
        band,
    })
}

fn parse_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_f64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("'{key}' must be a number"))),
        Some(_) => Err(ToolError::Validation(format!("'{key}' must be a number"))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{memory_store, CrsInfo, DataType, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Build a rows×cols DEM (1 m cells, top-down) from a row-major buffer.
    fn dem_of(rows: usize, cols: usize, vals: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
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
                r.set(0, row as isize, col as isize, vals[row * cols + col])
                    .unwrap();
            }
        }
        let id = memory_store::put_raster(r);
        memory_store::make_raster_memory_path(&id)
    }

    fn parse_csv(text: &str) -> Vec<Vec<String>> {
        text.lines()
            .skip(1)
            .filter(|l| !l.is_empty())
            .map(|l| l.split(',').map(str::to_string).collect())
            .collect()
    }

    fn run_csv(args: serde_json::Value) -> Vec<Vec<String>> {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = StorageCapacityTool.run(&args, &ctx()).unwrap();
        let path = out.outputs["output"].as_str().unwrap();
        parse_csv(&std::fs::read_to_string(path).unwrap())
    }

    /// A flat basin at z=0 with 1 m² cells: at level L, area = ncells, volume = L·ncells.
    #[test]
    fn flat_basin_area_and_volume() {
        // 3x3 flat DEM at elevation 0 -> 9 cells, each 1 m2.
        let dem = dem_of(3, 3, &[0.0; 9]);
        let tmp = std::env::temp_dir().join("storage_flat.csv");
        let rows = run_csv(json!({
            "dem": dem, "output": tmp.to_str().unwrap(),
            "min_elevation": 0.0, "max_elevation": 4.0, "num_levels": 4
        }));
        // Levels 0,1,2,3,4. At level L: area=9, volume=9*L.
        for r in &rows {
            let level: f64 = r[2].parse().unwrap();
            let area: f64 = r[3].parse().unwrap();
            let volume: f64 = r[4].parse().unwrap();
            assert!((area - 9.0).abs() < 1e-9, "area at level {level}");
            assert!(
                (volume - 9.0 * level).abs() < 1e-6,
                "volume at level {level}"
            );
        }
    }

    /// A bowl: volume is monotonic non-decreasing with level; area too.
    #[test]
    fn bowl_curve_is_monotonic() {
        // center low, rim high.
        let dem = dem_of(
            3,
            3,
            &[
                2.0, 2.0, 2.0, //
                2.0, 0.0, 2.0, //
                2.0, 2.0, 2.0,
            ],
        );
        let tmp = std::env::temp_dir().join("storage_bowl.csv");
        let rows = run_csv(json!({
            "dem": dem, "output": tmp.to_str().unwrap(), "num_levels": 4
        }));
        let mut last_area = -1.0;
        let mut last_vol = -1.0;
        for r in &rows {
            let area: f64 = r[3].parse().unwrap();
            let vol: f64 = r[4].parse().unwrap();
            assert!(area >= last_area - 1e-9, "area monotonic");
            assert!(vol >= last_vol - 1e-9, "volume monotonic");
            last_area = area;
            last_vol = vol;
        }
        // At the top level (2.0), all 9 cells flooded; volume = (2-0)*1 for the
        // center + 0 for the 8 rim cells = 2.
        let top = rows.last().unwrap();
        assert!((top[3].parse::<f64>().unwrap() - 9.0).abs() < 1e-9);
        assert!((top[4].parse::<f64>().unwrap() - 2.0).abs() < 1e-6);
    }

    /// A zone polygon restricts the sweep to cells inside it.
    #[test]
    fn zone_polygon_masks_cells() {
        use wbvector::{Coord, FieldDef, FieldType, GeometryType, Layer};
        let dem = dem_of(4, 4, &[0.0; 16]);
        // Polygon covering the top-left 2x2 block: x in [0,2], y in [2,4]
        // (row 0-1, col 0-1). DEM y_max = 4.
        let mut zl = Layer::new("z")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        zl.add_field(FieldDef::new("name", FieldType::Text));
        zl.add_feature(
            Some(Geometry::polygon(
                vec![
                    Coord::xy(0.0, 2.0),
                    Coord::xy(2.0, 2.0),
                    Coord::xy(2.0, 4.0),
                    Coord::xy(0.0, 4.0),
                ],
                vec![],
            )),
            &[("name", "nw".into())],
        )
        .unwrap();
        let zid = wbvector::memory_store::put_vector(zl);
        let zpath = wbvector::memory_store::make_vector_memory_path(&zid);

        let tmp = std::env::temp_dir().join("storage_zone.csv");
        let rows = run_csv(json!({
            "dem": dem, "output": tmp.to_str().unwrap(), "zones": zpath,
            "zone_id_field": "name", "min_elevation": 0.0, "max_elevation": 1.0, "num_levels": 1
        }));
        // 4 cells inside the polygon; at level 1 area=4.
        let top = rows.last().unwrap();
        assert_eq!(top[0], "nw", "zone id preserved");
        assert!(
            (top[3].parse::<f64>().unwrap() - 4.0).abs() < 1e-9,
            "4 masked cells"
        );
    }

    #[test]
    fn rejects_missing_output() {
        let dem = dem_of(2, 2, &[0.0; 4]);
        let args: ToolArgs = serde_json::from_value(json!({ "dem": dem })).unwrap();
        assert!(StorageCapacityTool.validate(&args).is_err());
    }
}
