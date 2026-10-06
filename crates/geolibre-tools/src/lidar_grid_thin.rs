//! GeoLibre tool: thin a LiDAR point cloud to one point per grid cell, **streaming**.
//!
//! The bundled `lidar_thin` / `lidar_thin_high_density` read the entire cloud
//! into a `Vec<PointRecord>` first. A decoded point is 336 bytes, so a
//! 17-million-point USGS 3DEP tile needs ~5.7 GB — more than a 32-bit
//! WebAssembly tab can ever allocate (2 GB per allocation, 4 GB total). Those
//! tools therefore abort in the browser exactly on the data one most wants to
//! thin.
//!
//! This tool never holds the input. It reads points one at a time (LAS, LAZ) or
//! one COPC node at a time, keeps the best candidate per cell of a global grid,
//! and only the survivors are materialized, so peak memory scales with the
//! **output** point count. The grid is anchored at the origin
//! (`floor(x / cell_size)`), which needs no bounding-box pass over the data.
//!
//! `method` picks the survivor in each cell: `nearest_center` (neutral, the
//! default), `lowest` (ground-biased, for bare-earth DEMs), or `highest`
//! (surface-biased, for DSMs and canopy). Ties keep the first point read, so the
//! result is deterministic. Points flagged withheld are skipped.

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wblidar::copc::CopcReader;
use wblidar::frontend::LidarFormat;
use wblidar::las::LasReader;
use wblidar::laz::LazReader;
use wblidar::{memory_store, PointCloud, PointReader, PointRecord};

use crate::args_common::{choice_or, opt_positive_f64, usize_or};
use crate::lidar_common::write_or_store_cloud;

/// Cap on surviving points on 32-bit targets (wasm32). Measured natively, each
/// survivor costs ~0.7 KB at peak (its 344 B cell-table entry plus the 336 B copy
/// in the output cloud; 1.4 GB for 2M points), and a wasm32 tab holds 4 GB with at
/// most 2 GB per allocation, so fail with an actionable error well before that.
/// Unlimited elsewhere.
#[cfg(target_pointer_width = "32")]
const MAX_POINTS: usize = 2_000_000;
#[cfg(not(target_pointer_width = "32"))]
const MAX_POINTS: usize = usize::MAX;

/// Largest accepted `points_per_cell`. Each offer is O(points_per_cell) (a sorted
/// insert of 344-byte entries), so an unbounded value would make dense cells
/// quadratic; no thinning workflow needs more than this per cell.
const MAX_POINTS_PER_CELL: usize = 1000;

/// Largest cell index magnitude accepted (well inside `i64`, exactly representable in `f64`).
const MAX_CELL_INDEX: f64 = 4.0e15;

const METHODS: [&str; 3] = ["nearest_center", "lowest", "highest"];

/// LAS "withheld" bit in `PointRecord::flags`. wblidar's LAS reader packs the top
/// three bits of the classification byte into `flags` (synthetic=1, key-point=2,
/// withheld=4), the same bit `wbtools_oss` tests in its own `is_withheld`.
const FLAG_WITHHELD: u8 = 0x04;

pub struct LidarGridThinTool;

impl Tool for LidarGridThinTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "lidar_grid_thin",
            display_name: "LiDAR Grid Thin",
            summary: "Thin a LiDAR point cloud to at most `points_per_cell` points (default one) per grid cell by streaming it, so peak memory scales with the output, not the input. Handles LAS/LAZ/COPC tiles with tens of millions of points that the bundled lidar_thin cannot load in a 32-bit WebAssembly tab.",
            category: ToolCategory::Lidar,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input LAS/LAZ/COPC point cloud (read as a stream, never fully loaded).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output point cloud path (.las, .laz or .copc.laz by extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Grid cell size in CRS units; at most points_per_cell points are kept per cell. Default 1.0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "Which points survive in each cell: nearest_center (default), lowest, or highest.",
                    required: false,
                },
                ToolParamSpec {
                    name: "points_per_cell",
                    description: "Most points to keep in each cell, best first by method (default 1, at most 1000). Cells with fewer points keep them all.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        input_path(args)?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = input_path(args)?;
        let prm = parse_params(args)?;
        let output = crate::common::parse_optional_output(args, "output")?;

        let mut grid = GridThinner::new(prm.cell_size, prm.method, prm.points_per_cell, MAX_POINTS);
        let crs = if memory_store::lidar_is_memory_path(input) {
            // Borrow the stored cloud; `load_input_cloud` would deep-copy it.
            let id = memory_store::lidar_path_to_id(input).ok_or_else(|| {
                ToolError::Validation("malformed in-memory lidar path".to_string())
            })?;
            let cloud = memory_store::get_lidar_arc_by_id(id).ok_or_else(|| {
                ToolError::Validation(format!("unknown in-memory lidar id '{id}'"))
            })?;
            for p in &cloud.points {
                if !grid.offer(p) {
                    break;
                }
            }
            cloud.crs.clone()
        } else {
            stream_points(input, |p| grid.offer(p))?
        };

        if grid.bad_cell {
            return Err(ToolError::Execution(format!(
                "cell_size {} is too small for these coordinates (cell index out of range); increase cell_size",
                prm.cell_size
            )));
        }
        if grid.overflowed {
            return Err(ToolError::Execution(format!(
                "more than {} points would be kept (cell_size {}, points_per_cell {}), which exceeds the memory of this build; increase cell_size or lower points_per_cell",
                grid.max_points, prm.cell_size, prm.points_per_cell
            )));
        }
        if crs.is_none() {
            ctx.progress
                .info("input has no readable CRS metadata; the output will have none");
        }
        let points_in = grid.seen;
        if grid.kept == 0 {
            return Err(ToolError::Execution(
                "no points to thin (input is empty or all points are withheld)".to_string(),
            ));
        }
        let points_out = grid.kept;
        ctx.progress.info(&format!(
            "kept {points_out} of {points_in} points ({:.1}%)",
            100.0 * points_out as f64 / points_in.max(1) as f64
        ));

        let cloud = PointCloud {
            points: grid.into_points(),
            crs,
        };
        let out_path = write_or_store_cloud(cloud, output)?;

        let mut outputs = std::collections::BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("points_in".to_string(), json!(points_in));
        outputs.insert("points_out".to_string(), json!(points_out));
        Ok(ToolRunResult { outputs })
    }
}

/// Keeps the best `per_cell` points seen so far for every occupied cell of an
/// origin-anchored grid.
struct GridThinner {
    cell: f64,
    method: Method,
    /// Most points kept in one cell.
    per_cell: usize,
    /// Occupied cell -> its kept `(rank, point)` pairs, best (lowest rank) first.
    /// Lower rank wins (see [`Self::rank`]).
    cells: HashMap<(i64, i64), Vec<(f64, PointRecord)>>,
    /// Points offered, excluding withheld ones.
    seen: usize,
    /// Points currently kept across all cells.
    kept: usize,
    /// Most points that may be kept before new ones are refused.
    max_points: usize,
    /// Set when `max_points` was exceeded; the result is then discarded.
    overflowed: bool,
    /// Set when a coordinate divided by `cell` does not fit a cell index.
    bad_cell: bool,
}

impl GridThinner {
    fn new(cell: f64, method: Method, per_cell: usize, max_points: usize) -> Self {
        Self {
            cell,
            method,
            per_cell: per_cell.max(1),
            cells: HashMap::new(),
            seen: 0,
            kept: 0,
            max_points,
            overflowed: false,
            bad_cell: false,
        }
    }

    /// Lower is better: squared distance to the cell centre, `z`, or `-z`.
    fn rank(&self, p: &PointRecord, cx: i64, cy: i64) -> f64 {
        match self.method {
            Method::NearestCenter => {
                let mx = (cx as f64 + 0.5) * self.cell;
                let my = (cy as f64 + 0.5) * self.cell;
                (p.x - mx).powi(2) + (p.y - my).powi(2)
            }
            Method::Lowest => p.z,
            Method::Highest => -p.z,
        }
    }

    /// Considers one point. Returns `false` once the memory cap is hit, so the
    /// caller can stop reading instead of decoding the rest of the file.
    fn offer(&mut self, p: &PointRecord) -> bool {
        if p.flags & FLAG_WITHHELD != 0 || !(p.x.is_finite() && p.y.is_finite()) {
            return true;
        }
        // A NaN z gives a NaN rank, which never compares less than anything, so
        // it could pin a cell forever; z only matters to the z-based methods.
        if self.method != Method::NearestCenter && !p.z.is_finite() {
            return true;
        }
        self.seen += 1;
        let (fx, fy) = ((p.x / self.cell).floor(), (p.y / self.cell).floor());
        // `as i64` saturates, which would silently merge distant cells into one.
        if fx.abs() >= MAX_CELL_INDEX || fy.abs() >= MAX_CELL_INDEX {
            self.bad_cell = true;
            return false;
        }
        let key = (fx as i64, fy as i64);
        let rank = self.rank(p, key.0, key.1);
        let slot = self.cells.entry(key).or_default();
        // Ties go after existing entries, so the first point read wins.
        let at = slot.partition_point(|e| e.0 <= rank);
        if slot.len() < self.per_cell {
            if self.kept >= self.max_points {
                self.overflowed = true;
                if slot.is_empty() {
                    self.cells.remove(&key);
                }
                return false;
            }
            // A fresh `Vec` of 336-byte records would otherwise jump to capacity 4
            // (~1.4 KB per cell even for one point), so grow one slot at a time.
            slot.reserve_exact(1);
            slot.insert(at, (rank, *p));
            self.kept += 1;
        } else if at < self.per_cell {
            // Drop the worst first so the insert reuses its slot: inserting into a
            // full `Vec` would grow it (to capacity 4, ~1.4 KB per cell).
            slot.pop();
            slot.insert(at, (rank, *p));
        }
        true
    }

    /// Survivors in a stable order (cell row-major, best point first within a
    /// cell) so output is reproducible.
    fn into_points(self) -> Vec<PointRecord> {
        let mut cells: Vec<((i64, i64), Vec<(f64, PointRecord)>)> =
            self.cells.into_iter().collect();
        cells.sort_unstable_by_key(|(k, _)| (k.1, k.0));
        // Exact capacity avoids a doubling reallocation (a transient 1.5-2x peak),
        // and each cell's Vec is freed as it is drained, so total memory stays
        // near one copy of the survivors.
        let mut out = Vec::with_capacity(self.kept);
        for (_, pts) in cells {
            out.extend(pts.into_iter().map(|(_, p)| p));
        }
        out
    }
}

/// Streams every point of a LAS/LAZ/COPC file to `visit` without materializing
/// the cloud, and returns the file's CRS. COPC is decoded one node at a time.
/// `visit` returns `false` to stop early (e.g. once the memory cap is hit).
fn stream_points<F: FnMut(&PointRecord) -> bool>(
    path: &str,
    mut visit: F,
) -> Result<Option<wblidar::Crs>, ToolError> {
    let open = || {
        File::open(path)
            .map_err(|e| ToolError::Execution(format!("failed opening input lidar: {e}")))
    };
    let read_err =
        |e: wblidar::Error| ToolError::Execution(format!("failed reading input lidar: {e}"));
    let format = LidarFormat::detect(std::path::Path::new(path)).map_err(read_err)?;

    // The LAS header + VLRs carry the CRS for LAS, LAZ and COPC alike.
    let crs = LasReader::new(BufReader::new(open()?))
        .ok()
        .and_then(|r| r.crs().cloned());

    match format {
        LidarFormat::Las => {
            let mut reader = LasReader::new(BufReader::new(open()?)).map_err(read_err)?;
            let mut p = PointRecord::default();
            while reader.read_point(&mut p).map_err(read_err)? {
                if !visit(&p) {
                    break;
                }
            }
        }
        LidarFormat::Laz => {
            let mut reader = LazReader::new(BufReader::new(open()?)).map_err(read_err)?;
            let mut p = PointRecord::default();
            while reader.read_point(&mut p).map_err(read_err)? {
                if !visit(&p) {
                    break;
                }
            }
        }
        LidarFormat::Copc => {
            let mut reader = CopcReader::new(BufReader::new(open()?)).map_err(read_err)?;
            let mut node: Vec<PointRecord> = Vec::new();
            'nodes: for key in reader.data_node_keys() {
                node.clear();
                reader.read_node(key, &mut node).map_err(read_err)?;
                for p in &node {
                    if !visit(p) {
                        break 'nodes;
                    }
                }
            }
        }
        _ => {
            return Err(ToolError::Validation(
                "input must be a LAS, LAZ or COPC file".to_string(),
            ))
        }
    }
    Ok(crs)
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    NearestCenter,
    Lowest,
    Highest,
}

struct Params {
    cell_size: f64,
    method: Method,
    points_per_cell: usize,
}

fn input_path(args: &ToolArgs) -> Result<&str, ToolError> {
    args.get("input")
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation("missing required parameter 'input'".to_string()))
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let cell_size = opt_positive_f64(args, "cell_size")?.unwrap_or(1.0);
    if !cell_size.is_finite() {
        return Err(ToolError::Validation(
            "'cell_size' must be a finite number".to_string(),
        ));
    }
    let method = match choice_or(args, "method", &METHODS, "nearest_center")? {
        "lowest" => Method::Lowest,
        "highest" => Method::Highest,
        _ => Method::NearestCenter,
    };
    let points_per_cell = usize_or(args, "points_per_cell", 1)?;
    if points_per_cell == 0 || points_per_cell > MAX_POINTS_PER_CELL {
        return Err(ToolError::Validation(format!(
            "'points_per_cell' must be between 1 and {MAX_POINTS_PER_CELL}"
        )));
    }
    Ok(Params {
        cell_size,
        method,
        points_per_cell,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn pt(x: f64, y: f64, z: f64) -> PointRecord {
        let mut p = PointRecord::default();
        p.x = x;
        p.y = y;
        p.z = z;
        p
    }

    fn cloud_path(pts: &[PointRecord]) -> String {
        let mut cloud = PointCloud::default();
        cloud.points = pts.to_vec();
        memory_store::make_lidar_memory_path(&memory_store::put_lidar(cloud))
    }

    fn run(v: serde_json::Value) -> (ToolRunResult, Vec<PointRecord>) {
        let args: ToolArgs = serde_json::from_value(v).unwrap();
        let out = LidarGridThinTool.run(&args, &ctx()).unwrap();
        let path = out.outputs["output"].as_str().unwrap().to_string();
        let cloud = crate::lidar_common::load_input_cloud(&path).unwrap();
        (out, cloud.points)
    }

    /// Four points in one 10 m cell collapse to one; a point in another cell survives.
    #[test]
    fn keeps_one_point_per_cell() {
        let input = cloud_path(&[
            pt(1.0, 1.0, 5.0),
            pt(2.0, 2.0, 6.0),
            pt(9.0, 9.0, 7.0),
            pt(4.0, 6.0, 8.0),
            pt(15.0, 5.0, 9.0),
        ]);
        let (out, pts) = run(json!({ "input": input, "cell_size": 10.0 }));
        assert_eq!(out.outputs["points_in"], json!(5));
        assert_eq!(out.outputs["points_out"], json!(2));
        assert_eq!(pts.len(), 2);
    }

    #[test]
    fn nearest_center_picks_closest_to_cell_centre() {
        // Cell (0,0) with size 10 has centre (5,5); (4.9,5.1) is closest.
        let input = cloud_path(&[pt(1.0, 1.0, 1.0), pt(4.9, 5.1, 2.0), pt(9.0, 9.0, 3.0)]);
        let (_o, pts) = run(json!({ "input": input, "cell_size": 10.0 }));
        assert_eq!(pts.len(), 1);
        assert_eq!(pts[0].z, 2.0);
    }

    #[test]
    fn lowest_and_highest_pick_extreme_z() {
        let pts_in = [pt(1.0, 1.0, 3.0), pt(2.0, 2.0, 9.0), pt(3.0, 3.0, 5.0)];
        let (_o, low) =
            run(json!({ "input": cloud_path(&pts_in), "cell_size": 10.0, "method": "lowest" }));
        let (_o, high) =
            run(json!({ "input": cloud_path(&pts_in), "cell_size": 10.0, "method": "highest" }));
        assert_eq!(low[0].z, 3.0);
        assert_eq!(high[0].z, 9.0);
    }

    /// Negative coordinates fall in their own cells (floor, not truncation toward zero).
    #[test]
    fn negative_coordinates_use_floor() {
        let input = cloud_path(&[pt(-0.5, 0.5, 1.0), pt(0.5, 0.5, 2.0)]);
        let (out, _p) = run(json!({ "input": input, "cell_size": 1.0 }));
        assert_eq!(out.outputs["points_out"], json!(2));
    }

    #[test]
    fn withheld_points_are_skipped() {
        let mut w = pt(1.0, 1.0, 100.0);
        w.flags = FLAG_WITHHELD;
        let input = cloud_path(&[w, pt(2.0, 2.0, 4.0)]);
        let (out, pts) = run(json!({ "input": input, "cell_size": 10.0, "method": "highest" }));
        assert_eq!(out.outputs["points_in"], json!(1));
        assert_eq!(pts[0].z, 4.0);
    }

    /// File round trip through the streaming LAS reader, preserving the CRS.
    #[test]
    fn streams_a_las_file() {
        let dir = std::env::temp_dir().join(format!("lgt_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let src = dir.join("in.las");
        let dst = dir.join("out.las");
        let mut cloud = PointCloud::default();
        cloud.crs = Some(wblidar::Crs {
            epsg: Some(32610),
            wkt: None,
        });
        for i in 0..100 {
            cloud.points.push(pt(i as f64 * 0.1, 0.5, i as f64));
        }
        cloud.write(&src).unwrap();
        let args: ToolArgs = serde_json::from_value(json!({
            "input": src.to_str().unwrap(), "output": dst.to_str().unwrap(), "cell_size": 1.0
        }))
        .unwrap();
        let out = LidarGridThinTool.run(&args, &ctx()).unwrap();
        assert_eq!(out.outputs["points_in"], json!(100));
        assert_eq!(out.outputs["points_out"], json!(10));
        let back = PointCloud::read(&dst).unwrap();
        assert_eq!(back.points.len(), 10);
        assert_eq!(back.crs.and_then(|c| c.epsg), Some(32610));
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Exceeding the cell cap is flagged instead of growing without bound.
    #[test]
    fn cell_cap_sets_overflow() {
        let mut g = GridThinner::new(1.0, Method::NearestCenter, 1, 2);
        for i in 0..5 {
            g.offer(&pt(i as f64 + 0.5, 0.5, 0.0));
        }
        assert!(g.overflowed);
        assert_eq!(g.kept, 2);
        assert_eq!(g.cells.len(), 2);
    }

    /// Keeps the N best per cell, best first, and never more than the cell holds.
    #[test]
    fn keeps_n_best_points_per_cell() {
        let pts_in = [
            pt(1.0, 1.0, 5.0),
            pt(2.0, 2.0, 3.0),
            pt(3.0, 3.0, 9.0),
            pt(4.0, 4.0, 1.0),
            pt(15.0, 5.0, 7.0), // a second cell with only one point
        ];
        let (out, low) = run(
            json!({ "input": cloud_path(&pts_in), "cell_size": 10.0, "method": "lowest", "points_per_cell": 2 }),
        );
        assert_eq!(out.outputs["points_out"], json!(3));
        // Cell (0,0): the two lowest, lowest first; then the lone point of cell (1,0).
        let z: Vec<f64> = low.iter().map(|p| p.z).collect();
        assert_eq!(z, vec![1.0, 3.0, 7.0]);
        let (_o, high) = run(
            json!({ "input": cloud_path(&pts_in), "cell_size": 10.0, "method": "highest", "points_per_cell": 2 }),
        );
        let z: Vec<f64> = high.iter().map(|p| p.z).collect();
        assert_eq!(z, vec![9.0, 5.0, 7.0]);
    }

    /// Equal ranks keep the first points read.
    #[test]
    fn ties_keep_first_read() {
        let mut a = pt(1.0, 1.0, 4.0);
        a.intensity = 1;
        let mut b = pt(2.0, 2.0, 4.0);
        b.intensity = 2;
        let mut c = pt(3.0, 3.0, 4.0);
        c.intensity = 3;
        let (_o, pts) = run(
            json!({ "input": cloud_path(&[a, b, c]), "cell_size": 10.0, "method": "lowest", "points_per_cell": 2 }),
        );
        let ids: Vec<u16> = pts.iter().map(|p| p.intensity).collect();
        assert_eq!(ids, vec![1, 2]);
    }

    /// The cap counts every kept point, not just occupied cells.
    #[test]
    fn point_cap_counts_all_kept_points() {
        let mut g = GridThinner::new(10.0, Method::Lowest, 3, 2);
        for i in 0..6 {
            g.offer(&pt(1.0, 1.0, i as f64));
        }
        assert!(g.overflowed);
        assert_eq!(g.kept, 2);
    }

    /// `offer` reports `false` as soon as the cap is hit so streaming can stop.
    #[test]
    fn offer_signals_stop_at_cap() {
        let mut g = GridThinner::new(1.0, Method::NearestCenter, 1, 2);
        assert!(g.offer(&pt(0.5, 0.5, 0.0)));
        assert!(g.offer(&pt(1.5, 0.5, 0.0)));
        assert!(!g.offer(&pt(2.5, 0.5, 0.0)));
    }

    /// A NaN z must neither pin a cell nor displace a valid point.
    #[test]
    fn nan_z_is_skipped_for_z_methods() {
        let input = cloud_path(&[pt(1.0, 1.0, f64::NAN), pt(2.0, 2.0, 4.0), pt(3.0, 3.0, 2.0)]);
        let (out, pts) = run(json!({ "input": input, "cell_size": 10.0, "method": "lowest" }));
        assert_eq!(out.outputs["points_in"], json!(2));
        assert_eq!(pts[0].z, 2.0);
    }

    /// The CRS survives a LAZ and a COPC round trip, not just uncompressed LAS.
    #[test]
    fn crs_survives_compressed_inputs() {
        let dir = std::env::temp_dir().join(format!("lgt_crs_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        for name in ["in.laz", "in.copc.laz"] {
            let src = dir.join(name);
            let dst = dir.join(format!("out_{name}.las"));
            let mut cloud = PointCloud::default();
            cloud.crs = Some(wblidar::Crs {
                epsg: Some(32610),
                wkt: None,
            });
            for i in 0..50 {
                cloud
                    .points
                    .push(pt(500000.0 + i as f64, 4000000.0, i as f64));
            }
            cloud.write(&src).unwrap();
            let args: ToolArgs = serde_json::from_value(json!({
                "input": src.to_str().unwrap(), "output": dst.to_str().unwrap(), "cell_size": 10.0
            }))
            .unwrap();
            let out = LidarGridThinTool.run(&args, &ctx()).unwrap();
            assert_eq!(out.outputs["points_in"], json!(50), "{name}");
            let back = PointCloud::read(&dst).unwrap();
            assert_eq!(back.crs.and_then(|c| c.epsg), Some(32610), "{name}");
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Withheld points are recognized through the real LAS decoder, not just the
    /// constant: write flags=4 to a file and let the streaming reader decode it.
    #[test]
    fn withheld_flag_matches_las_decoding() {
        let dir = std::env::temp_dir().join(format!("lgt_wh_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let src = dir.join("in.las");
        let mut w = pt(1.0, 1.0, 100.0);
        w.flags = 0x04;
        let mut cloud = PointCloud::default();
        cloud.points = vec![w, pt(2.0, 2.0, 4.0), pt(3.0, 3.0, 5.0)];
        cloud.write(&src).unwrap();
        let args: ToolArgs = serde_json::from_value(json!({
            "input": src.to_str().unwrap(), "cell_size": 10.0, "method": "highest"
        }))
        .unwrap();
        let out = LidarGridThinTool.run(&args, &ctx()).unwrap();
        assert_eq!(out.outputs["points_in"], json!(2));
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A cell size far too small for the coordinates is an error, not silent merging.
    #[test]
    fn tiny_cell_size_is_rejected_not_saturated() {
        let input = cloud_path(&[pt(500000.0, 4000000.0, 1.0), pt(500001.0, 4000001.0, 2.0)]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": input, "cell_size": 1e-12 })).unwrap();
        let err = LidarGridThinTool.run(&args, &ctx()).unwrap_err();
        assert!(err.to_string().contains("too small"), "{err}");
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            LidarGridThinTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.laz", "cell_size": 0 })).is_err());
        assert!(bad(json!({ "input": "a.laz", "cell_size": -2 })).is_err());
        assert!(bad(json!({ "input": "a.laz", "method": "bogus" })).is_err());
        assert!(bad(json!({ "input": "a.laz", "points_per_cell": 0 })).is_err());
        assert!(bad(json!({ "input": "a.laz", "points_per_cell": 1001 })).is_err());
        assert!(bad(json!({ "input": "a.laz", "points_per_cell": 1000 })).is_ok());
        assert!(bad(json!({ "input": "a.laz", "points_per_cell": 1.5 })).is_err());
        assert!(bad(json!({ "input": "a.laz", "points_per_cell": "3" })).is_ok());
        assert!(bad(json!({ "input": "a.laz" })).is_ok());
        assert!(bad(json!({ "input": "a.laz", "method": "lowest", "cell_size": "2.5" })).is_ok());
    }
}
