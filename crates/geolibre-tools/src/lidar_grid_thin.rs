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

use crate::args_common::{choice_or, opt_positive_f64};
use crate::lidar_common::write_or_store_cloud;

/// Cap on surviving points on 32-bit targets (wasm32). Each survivor costs ~360 B
/// in the cell table, and a single allocation above 2 GB aborts the tab, so fail
/// with an actionable error well before that. Unlimited elsewhere.
#[cfg(target_pointer_width = "32")]
const MAX_CELLS: usize = 2_000_000;
#[cfg(not(target_pointer_width = "32"))]
const MAX_CELLS: usize = usize::MAX;

const METHODS: [&str; 3] = ["nearest_center", "lowest", "highest"];

/// LAS "withheld" bit in the packed `flags` byte (synthetic=1, key-point=2, withheld=4).
const FLAG_WITHHELD: u8 = 0x04;

pub struct LidarGridThinTool;

impl Tool for LidarGridThinTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "lidar_grid_thin",
            display_name: "LiDAR Grid Thin",
            summary: "Thin a LiDAR point cloud to one point per grid cell by streaming it, so peak memory scales with the output, not the input. Handles LAS/LAZ/COPC tiles with tens of millions of points that the bundled lidar_thin cannot load in a 32-bit WebAssembly tab.",
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
                    description: "Grid cell size in CRS units; one point is kept per cell. Default 1.0.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "Which point survives in each cell: nearest_center (default), lowest, or highest.",
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

        let mut grid = GridThinner::new(prm.cell_size, prm.method, MAX_CELLS);
        let crs = if memory_store::lidar_is_memory_path(input) {
            let cloud = crate::lidar_common::load_input_cloud(input)?;
            for p in &cloud.points {
                grid.offer(p);
            }
            cloud.crs
        } else {
            stream_points(input, |p| grid.offer(p))?
        };

        if grid.overflowed {
            return Err(ToolError::Execution(format!(
                "more than {} cells would be kept at cell_size {}, which exceeds the memory of this build; increase cell_size",
                grid.max_cells, prm.cell_size
            )));
        }
        let points_in = grid.seen;
        if grid.cells.is_empty() {
            return Err(ToolError::Execution(
                "no points to thin (input is empty or all points are withheld)".to_string(),
            ));
        }
        let points_out = grid.cells.len();
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

/// Keeps the best point seen so far for every occupied cell of an origin-anchored grid.
struct GridThinner {
    cell: f64,
    method: Method,
    /// Occupied cell -> (rank, point). Lower rank wins (see [`Self::rank`]).
    cells: HashMap<(i64, i64), (f64, PointRecord)>,
    /// Points offered, excluding withheld ones.
    seen: usize,
    /// Most cells that may be occupied before new cells are refused.
    max_cells: usize,
    /// Set when `max_cells` was exceeded; the result is then discarded.
    overflowed: bool,
}

impl GridThinner {
    fn new(cell: f64, method: Method, max_cells: usize) -> Self {
        Self {
            cell,
            method,
            cells: HashMap::new(),
            seen: 0,
            max_cells,
            overflowed: false,
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

    fn offer(&mut self, p: &PointRecord) {
        if p.flags & FLAG_WITHHELD != 0 || !(p.x.is_finite() && p.y.is_finite()) {
            return;
        }
        self.seen += 1;
        let key = (
            (p.x / self.cell).floor() as i64,
            (p.y / self.cell).floor() as i64,
        );
        let rank = self.rank(p, key.0, key.1);
        if let Some(slot) = self.cells.get_mut(&key) {
            if rank < slot.0 {
                *slot = (rank, *p);
            }
        } else if self.cells.len() >= self.max_cells {
            self.overflowed = true;
        } else {
            self.cells.insert(key, (rank, *p));
        }
    }

    /// Survivors in a stable (cell-row-major) order so output is reproducible.
    fn into_points(self) -> Vec<PointRecord> {
        let mut entries: Vec<((i64, i64), PointRecord)> =
            self.cells.into_iter().map(|(k, (_, p))| (k, p)).collect();
        entries.sort_unstable_by_key(|(k, _)| (k.1, k.0));
        entries.into_iter().map(|(_, p)| p).collect()
    }
}

/// Streams every point of a LAS/LAZ/COPC file to `visit` without materializing
/// the cloud, and returns the file's CRS. COPC is decoded one node at a time.
fn stream_points<F: FnMut(&PointRecord)>(
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
                visit(&p);
            }
        }
        LidarFormat::Laz => {
            let mut reader = LazReader::new(BufReader::new(open()?)).map_err(read_err)?;
            let mut p = PointRecord::default();
            while reader.read_point(&mut p).map_err(read_err)? {
                visit(&p);
            }
        }
        LidarFormat::Copc => {
            let mut reader = CopcReader::new(BufReader::new(open()?)).map_err(read_err)?;
            let mut node: Vec<PointRecord> = Vec::new();
            for key in reader.data_node_keys() {
                node.clear();
                reader.read_node(key, &mut node).map_err(read_err)?;
                for p in &node {
                    visit(p);
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
    Ok(Params { cell_size, method })
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
        let mut g = GridThinner::new(1.0, Method::NearestCenter, 2);
        for i in 0..5 {
            g.offer(&pt(i as f64 + 0.5, 0.5, 0.0));
        }
        assert!(g.overflowed);
        assert_eq!(g.cells.len(), 2);
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
        assert!(bad(json!({ "input": "a.laz" })).is_ok());
        assert!(bad(json!({ "input": "a.laz", "method": "lowest", "cell_size": "2.5" })).is_ok());
    }
}
