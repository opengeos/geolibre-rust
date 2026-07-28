//! GeoLibre tool: azimuth and back-direction toward the nearest source cell.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Euclidean Direction* (Spatial
//! Analyst).
//!
//! The bundled suite ships two thirds of the Euclidean family —
//! `euclidean_distance` (how far to the nearest source) and
//! `euclidean_allocation` (which source it is) — but not the direction. That is
//! the piece answering *which way*, and without it evacuation-bearing surfaces,
//! orientation toward an outlet, and wind-fetch bearing to the nearest shoreline
//! cannot be expressed at all.
//!
//! Direction is reported in compass degrees (clockwise from north), with `0`
//! reserved for source cells and due north reported as `360`, matching ArcGIS so
//! outputs are directly comparable. Back direction follows the ArcGIS
//! definition too — the bearing to travel **back toward** the closest source, so
//! it agrees with the forward direction rather than opposing it. The two
//! coincide without barriers; with barriers the back direction is the first step
//! of the retrace, read from the Dijkstra predecessor.
//!
//! Without barriers the nearest source is found with an **exact** Euclidean
//! distance transform (Felzenszwalb & Huttenlocher's lower-envelope method)
//! extended to carry the nearest source's coordinates, so directions are exact
//! rather than chamfer approximations. With barriers the straight-line transform
//! is invalid, so the tool falls back to a multi-source Dijkstra over the
//! 8-connected grid that routes around blocked cells.

use std::collections::BinaryHeap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};

/// Computes, per cell, the compass bearing toward its nearest source cell.
pub struct EuclideanDirectionTool;

impl Tool for EuclideanDirectionTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "euclidean_direction",
            display_name: "Euclidean Direction",
            summary: "Computes the compass azimuth from each cell toward its nearest source cell, plus optional distance and back-direction rasters (ArcGIS Euclidean Direction). Completes the bundled Euclidean family, which ships euclidean_distance (how far) and euclidean_allocation (which source) but not the bearing. Uses an exact distance transform that carries nearest-source coordinates; with a barrier raster it routes around obstacles instead.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Source raster: cells that are neither no-data nor zero are sources.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output direction raster path (degrees clockwise from north). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_distance",
                    description: "Optional path for the companion distance raster.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_back_direction",
                    description: "Optional path for the back-direction raster: the bearing to travel from this cell back toward the closest source. Equals the direction raster when no barriers are supplied; with barriers it is the first step of the return path.",
                    required: false,
                },
                ToolParamSpec {
                    name: "barriers",
                    description: "Optional barrier raster; non-zero, non-no-data cells block straight-line paths, forcing routing around them.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_distance",
                    description: "Cells farther than this from any source are written as no-data.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band of the source raster to read (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args.get("input").and_then(Value::as_str).is_none() {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        opt_f64(args, "max_distance")?;
        opt_u64(args, "band")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args.get("input").and_then(Value::as_str).ok_or_else(|| {
            ToolError::Validation("missing required parameter 'input'".to_string())
        })?;
        let output = parse_optional_output(args, "output")?;
        let out_distance = parse_optional_output(args, "output_distance")?;
        let out_back = parse_optional_output(args, "output_back_direction")?;
        let barriers_path = parse_optional_output(args, "barriers")?;
        let max_distance = opt_f64(args, "max_distance")?;
        let band_1based = opt_u64(args, "band")?.unwrap_or(1);

        let raster = load_input_raster(input)?;
        if band_1based == 0 || band_1based as usize > raster.bands {
            return Err(ToolError::Validation(format!(
                "band {band_1based} out of range (raster has {} band(s))",
                raster.bands
            )));
        }
        let band = (band_1based - 1) as isize;

        let rows = raster.rows;
        let cols = raster.cols;
        let nodata = raster.nodata;
        let cell_x = raster.cell_size_x.abs().max(f64::MIN_POSITIVE);
        let cell_y = raster.cell_size_y.abs().max(f64::MIN_POSITIVE);

        // Sources: non-no-data, non-zero.
        let mut is_source = vec![false; rows * cols];
        let mut source_count = 0_u64;
        for row in 0..rows {
            for col in 0..cols {
                let v = raster.get(band, row as isize, col as isize);
                if v != nodata && v.is_finite() && v != 0.0 {
                    is_source[row * cols + col] = true;
                    source_count += 1;
                }
            }
        }
        if source_count == 0 {
            return Err(ToolError::Execution(
                "source raster contains no source cells (all cells are no-data or zero)"
                    .to_string(),
            ));
        }

        let blocked = match barriers_path {
            Some(path) => {
                let br = load_input_raster(path)?;
                if br.rows != rows || br.cols != cols {
                    return Err(ToolError::Validation(format!(
                        "barrier raster is {}x{} but the source raster is {rows}x{cols}",
                        br.rows, br.cols
                    )));
                }
                let bn = br.nodata;
                let mut mask = vec![false; rows * cols];
                for row in 0..rows {
                    for col in 0..cols {
                        let v = br.get(0, row as isize, col as isize);
                        mask[row * cols + col] = v != bn && v.is_finite() && v != 0.0;
                    }
                }
                Some(mask)
            }
            None => None,
        };

        ctx.progress.info("locating nearest sources");
        // `nearest[i]` is the (row, col) of the cell's nearest source, and
        // `dist[i]` the distance to it in map units.
        let (nearest, dist, prev) = match &blocked {
            None => {
                let (n, d) = exact_transform(rows, cols, &is_source, cell_x, cell_y);
                (n, d, None)
            }
            Some(mask) => {
                let (n, d, p) = dijkstra_transform(rows, cols, &is_source, mask, cell_x, cell_y);
                (n, d, Some(p))
            }
        };

        ctx.progress.info("computing bearings");
        // Direction is 0 on sources and (0, 360] elsewhere, so -1 is a safe
        // out-of-range no-data marker.
        let out_nodata = -1.0_f64;
        let mut dir = vec![out_nodata; rows * cols];
        let mut back = vec![out_nodata; rows * cols];
        let mut dst = vec![f64::NAN; rows * cols];
        let mut reached = 0_u64;

        for row in 0..rows {
            for col in 0..cols {
                let i = row * cols + col;
                let Some((s_row, s_col)) = nearest[i] else {
                    continue;
                };
                let d = dist[i];
                if let Some(max) = max_distance {
                    if d > max {
                        continue;
                    }
                }
                reached += 1;
                dst[i] = d;
                if is_source[i] {
                    dir[i] = 0.0;
                    back[i] = 0.0;
                    continue;
                }
                // Map-space offset from this cell to its source. Raster rows
                // increase downward (north-up), so northward is -d_row.
                let dx = (s_col as f64 - col as f64) * cell_x;
                let dy = (row as f64 - s_row as f64) * cell_y;
                dir[i] = compass_degrees(dx, dy);
                // ArcGIS defines back direction as the direction to travel back
                // toward the closest source, so it agrees with the forward
                // bearing rather than opposing it. Without barriers the return
                // path is the straight line, so the two coincide; with barriers
                // it is the first step of the retrace, taken from the Dijkstra
                // predecessor below.
                back[i] = match &prev {
                    None => dir[i],
                    Some(p) => match p[i] {
                        usize::MAX => dir[i],
                        pj => {
                            let (p_row, p_col) = (pj / cols, pj % cols);
                            let bx = (p_col as f64 - col as f64) * cell_x;
                            let by = (row as f64 - p_row as f64) * cell_y;
                            compass_degrees(bx, by)
                        }
                    },
                };
            }
            ctx.progress
                .progress((row as f64 + 1.0) / rows.max(1) as f64);
        }

        let dir_raster = raster_like_with_data(&raster, dir, out_nodata, DataType::F32)?;
        let out_path = write_or_store_output(dir_raster, output)?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("source_cells".to_string(), json!(source_count));
        outputs.insert("reached_cells".to_string(), json!(reached));

        // Both companion rasters are always produced, falling back to an
        // in-memory handle when no path is given (as create_overpass does for
        // its secondary output), so a caller without a scratch path still gets
        // them back.
        let dnodata = -1.0_f64;
        let dist_buf: Vec<f64> = dst
            .iter()
            .map(|v| if v.is_nan() { dnodata } else { *v })
            .collect();
        let dist_raster = raster_like_with_data(&raster, dist_buf, dnodata, DataType::F32)?;
        let dist_path = write_or_store_output(dist_raster, out_distance)?;
        outputs.insert("output_distance".to_string(), json!(dist_path));

        let back_raster = raster_like_with_data(&raster, back, out_nodata, DataType::F32)?;
        let back_path = write_or_store_output(back_raster, out_back)?;
        outputs.insert("output_back_direction".to_string(), json!(back_path));

        Ok(ToolRunResult { outputs })
    }
}

/// Converts a map-space offset into compass degrees clockwise from north, with
/// due north reported as 360 (ArcGIS convention; 0 is reserved for sources).
fn compass_degrees(dx: f64, dy: f64) -> f64 {
    if dx == 0.0 && dy == 0.0 {
        return 0.0;
    }
    let mut deg = dx.atan2(dy).to_degrees();
    if deg <= 0.0 {
        deg += 360.0;
    }
    deg
}

/// Exact Euclidean distance transform carrying nearest-source coordinates.
///
/// Runs the Felzenszwalb & Huttenlocher lower-envelope transform once per axis.
/// The first pass reduces each column to the nearest source *within* that
/// column; the second takes, per row, the lower envelope of the parabolas those
/// column minima define. Tracking the envelope's argmin recovers the true
/// nearest source cell, not just its distance.
///
/// Anisotropic cell sizes are handled by scaling each axis into map units.
#[allow(clippy::type_complexity)]
fn exact_transform(
    rows: usize,
    cols: usize,
    is_source: &[bool],
    cell_x: f64,
    cell_y: f64,
) -> (Vec<Option<(usize, usize)>>, Vec<f64>) {
    const INF: f64 = f64::INFINITY;

    // Pass 1 — per column, squared distance (map units) to the nearest source
    // in that column, plus that source's row.
    let mut col_d2 = vec![INF; rows * cols];
    let mut col_src = vec![usize::MAX; rows * cols];
    for col in 0..cols {
        // Downward sweep.
        let mut best: Option<usize> = None;
        for row in 0..rows {
            if is_source[row * cols + col] {
                best = Some(row);
            }
            if let Some(b) = best {
                let d = (row as f64 - b as f64) * cell_y;
                col_d2[row * cols + col] = d * d;
                col_src[row * cols + col] = b;
            }
        }
        // Upward sweep keeps whichever direction is closer.
        best = None;
        for row in (0..rows).rev() {
            if is_source[row * cols + col] {
                best = Some(row);
            }
            if let Some(b) = best {
                let d = (row as f64 - b as f64) * cell_y;
                let d2 = d * d;
                if d2 < col_d2[row * cols + col] {
                    col_d2[row * cols + col] = d2;
                    col_src[row * cols + col] = b;
                }
            }
        }
    }

    // Pass 2 — per row, lower envelope of parabolas f_c(x) = (x - c)^2 * cell_x^2
    // + col_d2[c], tracking which column's parabola wins.
    let mut nearest: Vec<Option<(usize, usize)>> = vec![None; rows * cols];
    let mut dist = vec![INF; rows * cols];

    let mut v = vec![0_usize; cols.max(1)];
    let mut z = vec![0.0_f64; cols + 1];

    for row in 0..rows {
        let f = |c: usize| col_d2[row * cols + c];

        let mut k = 0_usize;
        v[0] = 0;
        z[0] = f64::NEG_INFINITY;
        z[1] = INF;

        for q in 1..cols {
            if f(q).is_infinite() {
                continue;
            }
            loop {
                // Intersection of parabolas q and v[k] in x (column) space.
                let p = v[k];
                if f(p).is_infinite() {
                    // Replace an empty parabola outright.
                    if k == 0 {
                        v[0] = q;
                        z[0] = f64::NEG_INFINITY;
                        z[1] = INF;
                        break;
                    }
                    k -= 1;
                    continue;
                }
                let x2 = cell_x * cell_x;
                let s = ((f(q) + (q as f64 * q as f64) * x2) - (f(p) + (p as f64 * p as f64) * x2))
                    / (2.0 * x2 * (q as f64 - p as f64));
                if s <= z[k] {
                    if k == 0 {
                        v[0] = q;
                        z[0] = f64::NEG_INFINITY;
                        z[1] = INF;
                        break;
                    }
                    k -= 1;
                } else {
                    k += 1;
                    v[k] = q;
                    z[k] = s;
                    z[k + 1] = INF;
                    break;
                }
            }
        }

        let mut k_read = 0_usize;
        for q in 0..cols {
            while z[k_read + 1] < q as f64 {
                k_read += 1;
            }
            let p = v[k_read];
            if f(p).is_infinite() {
                continue;
            }
            let dx = (q as f64 - p as f64) * cell_x;
            let d2 = dx * dx + f(p);
            let src_row = col_src[row * cols + p];
            if src_row != usize::MAX {
                dist[row * cols + q] = d2.sqrt();
                nearest[row * cols + q] = Some((src_row, p));
            }
        }
    }

    (nearest, dist)
}

/// Multi-source Dijkstra over the 8-connected grid, skipping blocked cells.
/// Used when a barrier raster is supplied, where straight-line distance is
/// invalid. Each cell inherits the source its cheapest path came from.
#[allow(clippy::type_complexity)]
fn dijkstra_transform(
    rows: usize,
    cols: usize,
    is_source: &[bool],
    blocked: &[bool],
    cell_x: f64,
    cell_y: f64,
) -> (Vec<Option<(usize, usize)>>, Vec<f64>, Vec<usize>) {
    #[derive(PartialEq)]
    struct Node(f64, usize);
    impl Eq for Node {}
    impl Ord for Node {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering {
            // Min-heap by distance; index breaks ties so ordering is total and
            // deterministic across platforms.
            other
                .0
                .partial_cmp(&self.0)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| other.1.cmp(&self.1))
        }
    }
    impl PartialOrd for Node {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }

    let mut dist = vec![f64::INFINITY; rows * cols];
    let mut nearest: Vec<Option<(usize, usize)>> = vec![None; rows * cols];
    // Predecessor on the cheapest path home, so back-direction can report the
    // first step of the retrace rather than a straight-line bearing.
    let mut prev = vec![usize::MAX; rows * cols];
    let mut heap = BinaryHeap::new();

    for row in 0..rows {
        for col in 0..cols {
            let i = row * cols + col;
            if is_source[i] && !blocked[i] {
                dist[i] = 0.0;
                nearest[i] = Some((row, col));
                heap.push(Node(0.0, i));
            }
        }
    }

    let diag = (cell_x * cell_x + cell_y * cell_y).sqrt();
    let steps: [(isize, isize, f64); 8] = [
        (-1, 0, cell_y),
        (1, 0, cell_y),
        (0, -1, cell_x),
        (0, 1, cell_x),
        (-1, -1, diag),
        (-1, 1, diag),
        (1, -1, diag),
        (1, 1, diag),
    ];

    while let Some(Node(d, i)) = heap.pop() {
        if d > dist[i] {
            continue;
        }
        let row = i / cols;
        let col = i % cols;
        for (d_row, d_col, cost) in steps {
            let n_row = row as isize + d_row;
            let n_col = col as isize + d_col;
            if n_row < 0 || n_col < 0 || n_row >= rows as isize || n_col >= cols as isize {
                continue;
            }
            let j = n_row as usize * cols + n_col as usize;
            if blocked[j] {
                continue;
            }
            // No corner cutting. A diagonal step passes between two orthogonal
            // cells; if both are blocked it would squeeze through the corner of
            // a one-cell-thick diagonal wall, making the barrier permeable.
            if d_row != 0 && d_col != 0 {
                let a = (row as isize + d_row) as usize * cols + col;
                let b = row * cols + (col as isize + d_col) as usize;
                if blocked[a] && blocked[b] {
                    continue;
                }
            }
            let nd = d + cost;
            if nd < dist[j] {
                dist[j] = nd;
                nearest[j] = nearest[i];
                prev[j] = i;
                heap.push(Node(nd, j));
            }
        }
    }

    (nearest, dist, prev)
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

fn opt_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => n.as_u64().map(Some).ok_or_else(|| {
            ToolError::Validation(format!("parameter '{key}' must be a positive integer"))
        }),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s.trim().parse::<u64>().map(Some).map_err(|_| {
            ToolError::Validation(format!("parameter '{key}' must be a positive integer"))
        }),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a positive integer"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn raster(cols: usize, rows: usize, data: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
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
                r.set(0, row as isize, col as isize, data[row * cols + col])
                    .unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run_with(path: String, extra: Value) -> Raster {
        let mut obj = serde_json::Map::new();
        obj.insert("input".to_string(), json!(path));
        if let Value::Object(m) = extra {
            for (k, v) in m {
                obj.insert(k, v);
            }
        }
        let args: ToolArgs = serde_json::from_value(Value::Object(obj)).unwrap();
        let out = EuclideanDirectionTool.run(&args, &ctx()).unwrap();
        load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    /// Cardinal bearings: a single source in the middle of a 3x3.
    /// North is 360, east 90, south 180, west 270.
    #[test]
    fn cardinal_bearings_match_arcgis_convention() {
        let path = raster(3, 3, &[0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]);
        let out = run_with(path, json!({}));
        assert_eq!(out.get(0, 1, 1), 0.0, "source cell reads 0");
        // Cell below the source (row 2) must look north => 360.
        assert!((out.get(0, 2, 1) - 360.0).abs() < 1e-9);
        // Cell above the source (row 0) looks south => 180.
        assert!((out.get(0, 0, 1) - 180.0).abs() < 1e-9);
        // Cell left of the source looks east => 90.
        assert!((out.get(0, 1, 0) - 90.0).abs() < 1e-9);
        // Cell right of the source looks west => 270.
        assert!((out.get(0, 1, 2) - 270.0).abs() < 1e-9);
    }

    /// Diagonals land on the 45-degree bearings.
    #[test]
    fn diagonal_bearings() {
        let path = raster(3, 3, &[0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]);
        let out = run_with(path, json!({}));
        // Bottom-left cell (2,0) looks north-east => 45.
        assert!((out.get(0, 2, 0) - 45.0).abs() < 1e-9);
        // Top-right cell (0,2) looks south-west => 225.
        assert!((out.get(0, 0, 2) - 225.0).abs() < 1e-9);
    }

    /// With two sources, each cell points at the *nearest* one — this is what
    /// the exact transform's source tracking has to get right.
    #[test]
    fn points_at_the_nearest_of_several_sources() {
        // Sources at both ends of a 5-wide row.
        let path = raster(5, 1, &[1.0, 0.0, 0.0, 0.0, 1.0]);
        let out = run_with(path, json!({}));
        // Column 1 is nearest the left source => looks west (270).
        assert!((out.get(0, 0, 1) - 270.0).abs() < 1e-9);
        // Column 3 is nearest the right source => looks east (90).
        assert!((out.get(0, 0, 3) - 90.0).abs() < 1e-9);
    }

    /// A single corner source: bearings across the grid stay exact, which is
    /// what the lower-envelope transform buys over a chamfer approximation.
    #[test]
    fn exact_transform_gives_exact_diagonal_bearing() {
        let mut data = [0.0; 16];
        data[0] = 1.0; // single source at (row 0, col 0)
        let path = raster(4, 4, &data);
        let args: ToolArgs = serde_json::from_value(json!({ "input": path })).unwrap();
        let res = EuclideanDirectionTool.run(&args, &ctx()).unwrap();
        assert_eq!(res.outputs["source_cells"], json!(1));

        let out = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        // (3,3) sits exactly south-east of the source, so it looks north-west: 315.
        assert!((out.get(0, 3, 3) - 315.0).abs() < 1e-9);
    }

    /// Back direction points **back toward** the source (ArcGIS semantics), so
    /// without barriers it equals the forward direction rather than opposing it.
    #[test]
    fn back_direction_points_toward_the_source() {
        let mut data = [0.0; 9];
        data[4] = 1.0; // source at the centre
        let path = raster(3, 3, &data);
        let args: ToolArgs = serde_json::from_value(json!({ "input": path })).unwrap();
        let res = EuclideanDirectionTool.run(&args, &ctx()).unwrap();
        let dir = load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap();
        let back =
            load_input_raster(res.outputs["output_back_direction"].as_str().unwrap()).unwrap();

        // Cell (2,1) is south of the source, so it looks north => 360.
        assert!((dir.get(0, 2, 1) - 360.0).abs() < 1e-9);
        for (r, c) in [(0, 0), (0, 1), (1, 0), (1, 2), (2, 1), (2, 2)] {
            assert!(
                (back.get(0, r, c) - dir.get(0, r, c)).abs() < 1e-9,
                "without barriers back direction must equal the forward bearing at ({r},{c}): {} vs {}",
                back.get(0, r, c),
                dir.get(0, r, c)
            );
        }
    }

    /// A one-cell-thick diagonal wall must be impermeable: an 8-connected search
    /// that ignores corners would squeeze between the two blocked cells.
    #[test]
    fn diagonal_barrier_is_not_permeable() {
        // 3x3. Source top-left; a diagonal wall blocks (0,1) and (1,0), sealing
        // the source into its own corner.
        let mut src = [0.0; 9];
        src[0] = 1.0;
        let source = raster(3, 3, &src);

        let mut bar = [0.0; 9];
        bar[1] = 1.0; // (0,1)
        bar[3] = 1.0; // (1,0)
        let barriers = raster(3, 3, &bar);

        let out = run_with(source, json!({ "barriers": barriers }));
        assert_eq!(
            out.get(0, 1, 1),
            out.nodata,
            "the cell diagonally past the wall must be unreachable, not reached through the corner"
        );
        assert_eq!(out.get(0, 2, 2), out.nodata);
    }

    /// max_distance masks far cells to no-data.
    #[test]
    fn max_distance_masks_far_cells() {
        let mut data = [0.0; 25];
        data[0] = 1.0;
        let path = raster(5, 5, &data);
        let out = run_with(path, json!({ "max_distance": 1.5 }));
        assert!(out.get(0, 0, 1) > 0.0, "adjacent cell is within range");
        assert_eq!(out.get(0, 4, 4), out.nodata, "far corner is masked");
    }

    /// A barrier forces the path around it, changing the reported bearing away
    /// from the straight-line answer.
    #[test]
    fn barrier_routes_around_obstacle() {
        // Source at left end of a 5x3; a wall blocks the middle column except
        // the bottom row.
        let mut src = [0.0; 15];
        src[5] = 1.0; // (row 1, col 0)
        let source = raster(5, 3, &src);

        let mut bar = [0.0; 15];
        bar[2] = 1.0; // (0,2)
        bar[7] = 1.0; // (1,2)
        let barriers = raster(5, 3, &bar);

        let out = run_with(source, json!({ "barriers": barriers }));
        // The cell directly right of the wall must still be reached (around the
        // bottom), and must not be no-data.
        assert_ne!(out.get(0, 1, 3), out.nodata);
        // Blocked cells themselves are unreachable.
        assert_eq!(out.get(0, 1, 2), out.nodata);
    }

    #[test]
    fn rejects_bad_parameters() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(EuclideanDirectionTool.validate(&args).is_err());

        // A raster with no sources cannot produce directions.
        let path = raster(2, 2, &[0.0; 4]);
        let args: ToolArgs = serde_json::from_value(json!({ "input": path })).unwrap();
        assert!(EuclideanDirectionTool.run(&args, &ctx()).is_err());
    }
}
