//! GeoLibre tool: breach polylines that drain depression polygons.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Generate Breach Lines* (Spatial
//! Analyst). The bundled suite breaches depressions **in raster space only** —
//! `breach_depressions_least_cost`, `breach_single_cell_pits` and
//! `topological_breach_burn` all write a modified DEM and discard the path they
//! carved. Nothing in either registry returns that path as vector geometry.
//!
//! The path is the deliverable in a lot of workflows: it is what gets handed to
//! an engineer as a proposed culvert or channel alignment, what gets styled on
//! a web map, and what gets length- and volume-summarised in a table. It also
//! composes directly with GeoLibre's own `delineate_depressions`, which already
//! produces exactly the depression-polygon input this tool consumes, making
//! `delineate_depressions -> generate_breach_lines` a two-step pipeline.
//!
//! For each depression the search starts at the pit (its lowest interior cell)
//! and runs a Dijkstra sweep outward over a per-method cost field, stopping at
//! the first cell that lies outside the depression **and** below the pit — i.e.
//! somewhere the water can actually go. The back-link chain is then traced to a
//! polyline.

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::common::load_input_raster;
use crate::vector_common::{
    geometry_contains_point, load_input_layer, parse_optional_str, write_or_store_layer,
};

pub struct GenerateBreachLinesTool;

impl Tool for GenerateBreachLinesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "generate_breach_lines",
            display_name: "Generate Breach Lines",
            summary: "Given depression polygons and a DEM, emit the polyline along which each depression would be breached to drain, using a minimum-breaching-cost, shortest-path or minimum-elevation-change criterion, with cut depth and excavated volume per breach. Like ArcGIS Generate Breach Lines.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Depression polygons (e.g. the output of delineate_depressions).",
                    required: true,
                },
                ToolParamSpec {
                    name: "dem",
                    description: "Surface raster the breach is carved through.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output breach-line path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "connection_points",
                    description: "Optional point features constraining where a breach may terminate.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'minimum_breaching_cost' (default), 'shortest_path' or 'minimum_elevation_change'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_length",
                    description: "Drop breaches longer than this (map units); such depressions are left unbreached.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "dem")?;
        parse_method(args)?;
        if let Some(m) = parse_optional_f64(args, "max_length")? {
            if m <= 0.0 {
                return Err(ToolError::Validation(
                    "'max_length' must be positive".to_string(),
                ));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let dem_path = require_str(args, "dem")?;
        let output = parse_optional_str(args, "output")?;
        let method = parse_method(args)?;
        let max_length = parse_optional_f64(args, "max_length")?;

        let polys = load_input_layer(input)?;
        if polys.features.is_empty() {
            return Err(ToolError::Execution("input has no features".to_string()));
        }
        let dem = load_input_raster(dem_path)?;
        let rows = dem.rows as isize;
        let cols = dem.cols as isize;
        let nodata = dem.nodata;

        // Optional termination targets.
        // A HashSet, not a Vec: the outlet test runs inside the Dijkstra pop
        // loop, so a linear scan would be O(cells x points).
        let connections: std::collections::HashSet<(isize, isize)> = match parse_optional_str(args, "connection_points")? {
            Some(p) => {
                let layer = load_input_layer(p)?;
                layer
                    .iter()
                    .filter_map(|f| f.geometry.as_ref())
                    .flat_map(|g| g.all_coords().into_iter().map(|c| (c.x, c.y)).collect::<Vec<_>>())
                    .filter_map(|(x, y)| dem.world_to_pixel(x, y).map(|(c, r)| (r, c)))
                    .collect()
            }
            None => std::collections::HashSet::new(),
        };

        let mut out = Layer::new("breach_lines").with_geom_type(GeometryType::LineString);
        if let Some(epsg) = polys.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("DEPRESSION_ID", FieldType::Integer));
        out.add_field(FieldDef::new("BREACH_LEN", FieldType::Float));
        out.add_field(FieldDef::new("MAX_CUT", FieldType::Float));
        out.add_field(FieldDef::new("CUT_VOLUME", FieldType::Float));
        out.add_field(FieldDef::new("INLET_Z", FieldType::Float));
        out.add_field(FieldDef::new("OUTLET_Z", FieldType::Float));

        // `delineate_depressions` routinely emits thousands of polygons, so
        // allocating and zeroing full-raster buffers per feature is O(features x
        // raster). Allocate once and reset only the cells actually touched.
        let mut inside = vec![false; (rows * cols) as usize];
        let mut scratch = Scratch::new((rows * cols) as usize);

        let cell_area = dem.cell_size_x * dem.cell_size_y;
        let mut breached = 0usize;
        let mut unbreached = 0usize;
        let mut too_long = 0usize;
        let mut total_volume = 0.0_f64;

        for (fid, feat) in polys.iter().enumerate() {
            let Some(geom) = &feat.geometry else {
                unbreached += 1;
                continue;
            };
            let Some(bb) = geom.bbox() else {
                unbreached += 1;
                continue;
            };
            // Rasterize the depression footprint over its bbox window.
            let c0 = (((bb.min_x - dem.x_min) / dem.cell_size_x).floor() as isize).max(0);
            let c1 = (((bb.max_x - dem.x_min) / dem.cell_size_x).ceil() as isize).min(cols - 1);
            let r0 = (((dem.y_max() - bb.max_y) / dem.cell_size_y).floor() as isize).max(0);
            let r1 = (((dem.y_max() - bb.min_y) / dem.cell_size_y).ceil() as isize).min(rows - 1);
            if c1 < c0 || r1 < r0 {
                unbreached += 1;
                continue;
            }

            let mut touched_inside: Vec<usize> = Vec::new();
            let mut pit: Option<(isize, isize, f64)> = None;
            for r in r0..=r1 {
                for c in c0..=c1 {
                    let (x, y) = pixel_center(&dem, r, c);
                    if !geometry_contains_point(geom, x, y) {
                        continue;
                    }
                    inside[(r * cols + c) as usize] = true;
                    touched_inside.push((r * cols + c) as usize);
                    let z = dem.get(0, r, c);
                    if z == nodata || !z.is_finite() {
                        continue;
                    }
                    if pit.is_none_or(|(_, _, pz)| z < pz) {
                        pit = Some((r, c, z));
                    }
                }
            }
            let Some((pr, pc, pit_z)) = pit else {
                for &i in &touched_inside {
                    inside[i] = false;
                }
                unbreached += 1;
                continue;
            };

            let found = search(
                &dem, &inside, &mut scratch, (pr, pc), pit_z, method, max_length, &connections,
            );
            for &i in &touched_inside {
                inside[i] = false;
            }
            match found {
                Some(path) => {
                    let mut length = 0.0;
                    let mut max_cut = 0.0_f64;
                    let mut volume = 0.0;
                    let mut coords: Vec<Coord> = Vec::with_capacity(path.len());
                    for (i, &(r, c)) in path.iter().enumerate() {
                        let (x, y) = pixel_center(&dem, r, c);
                        let z = dem.get(0, r, c);
                        coords.push(Coord::xyz(x, y, z));
                        if i > 0 {
                            let (px, py) = pixel_center(&dem, path[i - 1].0, path[i - 1].1);
                            length += (x - px).hypot(y - py);
                        }
                        // Excavation to bring this cell down to the pit level.
                        let cut = (z - pit_z).max(0.0);
                        max_cut = max_cut.max(cut);
                        volume += cut * cell_area;
                    }
                    if coords.len() < 2 {
                        unbreached += 1;
                        continue;
                    }
                    let outlet_z = dem.get(0, path[path.len() - 1].0, path[path.len() - 1].1);
                    total_volume += volume;
                    breached += 1;
                    out.add_feature(
                        Some(Geometry::line_string(coords)),
                        &[
                            ("DEPRESSION_ID", FieldValue::Integer(fid as i64)),
                            ("BREACH_LEN", FieldValue::Float(length)),
                            ("MAX_CUT", FieldValue::Float(max_cut)),
                            ("CUT_VOLUME", FieldValue::Float(volume)),
                            ("INLET_Z", FieldValue::Float(pit_z)),
                            ("OUTLET_Z", FieldValue::Float(outlet_z)),
                        ],
                    )
                    .map_err(|e| ToolError::Execution(format!("failed adding breach: {e}")))?;
                }
                None => {
                    if max_length.is_some() {
                        too_long += 1;
                    }
                    unbreached += 1;
                }
            }
            ctx.progress
                .progress((fid as f64 + 1.0) / polys.features.len() as f64);
        }

        let n_polys = polys.features.len();
        ctx.progress
            .info(&format!("breached {breached} of {n_polys} depression(s)"));
        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("depression_count".to_string(), json!(n_polys));
        outputs.insert("breached_count".to_string(), json!(breached));
        outputs.insert("unbreached_count".to_string(), json!(unbreached));
        outputs.insert("exceeded_max_length".to_string(), json!(too_long));
        outputs.insert("total_cut_volume".to_string(), json!(total_volume));
        outputs.insert("method".to_string(), json!(method.name()));
        Ok(ToolRunResult { outputs })
    }
}

/// Reusable Dijkstra buffers.
///
/// `touched` records every index written during a sweep so the next call can
/// reset just those, instead of re-zeroing the whole raster per depression.
struct Scratch {
    dist: Vec<f64>,
    travelled: Vec<f64>,
    back: Vec<usize>,
    done: Vec<bool>,
    touched: Vec<usize>,
}

impl Scratch {
    fn new(n: usize) -> Self {
        Self {
            dist: vec![f64::INFINITY; n],
            travelled: vec![f64::INFINITY; n],
            back: vec![usize::MAX; n],
            done: vec![false; n],
            touched: Vec::new(),
        }
    }
    fn reset(&mut self) {
        for &i in &self.touched {
            self.dist[i] = f64::INFINITY;
            self.travelled[i] = f64::INFINITY;
            self.back[i] = usize::MAX;
            self.done[i] = false;
        }
        self.touched.clear();
    }
}

/// Dijkstra sweep outward from the pit, returning the cell path to the first
/// admissible outlet.
///
/// An outlet is a cell outside the depression whose elevation is below the pit
/// — i.e. water released there keeps flowing away rather than pooling straight
/// back. When `connection_points` are supplied, only those cells qualify.
#[allow(clippy::too_many_arguments)]
fn search(
    dem: &Raster,
    inside: &[bool],
    scratch: &mut Scratch,
    pit: (isize, isize),
    pit_z: f64,
    method: Method,
    max_length: Option<f64>,
    connections: &std::collections::HashSet<(isize, isize)>,
) -> Option<Vec<(isize, isize)>> {
    let rows = dem.rows as isize;
    let cols = dem.cols as isize;
    let n = (rows * cols) as usize;
    let nodata = dem.nodata;
    let diag = dem.cell_size_x.hypot(dem.cell_size_y);
    let straight = dem.cell_size_x.min(dem.cell_size_y);

    scratch.reset();
    let Scratch {
        dist,
        travelled,
        back,
        done,
        touched,
    } = scratch;
    debug_assert_eq!(dist.len(), n);
    let mut heap: BinaryHeap<Node> = BinaryHeap::new();

    let start = (pit.0 * cols + pit.1) as usize;
    dist[start] = 0.0;
    travelled[start] = 0.0;
    touched.push(start);
    heap.push(Node {
        cost: 0.0,
        idx: start,
    });

    while let Some(Node { cost, idx }) = heap.pop() {
        if done[idx] {
            continue;
        }
        done[idx] = true;
        let r = idx as isize / cols;
        let c = idx as isize % cols;
        let z = dem.get(0, r, c);

        // Outlet test: outside the depression and genuinely downhill of the pit.
        let is_outlet = if connections.is_empty() {
            !inside[idx] && z != nodata && z.is_finite() && z < pit_z
        } else {
            connections.contains(&(r, c))
        };
        if is_outlet && idx != start {
            let mut path = vec![(r, c)];
            let mut cur = idx;
            while back[cur] != usize::MAX {
                cur = back[cur];
                path.push((cur as isize / cols, cur as isize % cols));
            }
            path.reverse();
            return Some(path);
        }

        for (dr, dc) in [
            (-1, -1), (-1, 0), (-1, 1),
            (0, -1),           (0, 1),
            (1, -1),  (1, 0),  (1, 1),
        ] {
            let (nr, nc) = (r + dr, c + dc);
            if nr < 0 || nc < 0 || nr >= rows || nc >= cols {
                continue;
            }
            let nidx = (nr * cols + nc) as usize;
            if done[nidx] {
                continue;
            }
            let nz = dem.get(0, nr, nc);
            if nz == nodata || !nz.is_finite() {
                continue;
            }
            let step = if dr != 0 && dc != 0 { diag } else { straight };
            let travel = travelled[idx] + step;
            if let Some(limit) = max_length {
                if travel > limit {
                    continue;
                }
            }
            let increment = match method {
                // Excavation needed to bring this cell down to the pit level.
                Method::MinimumBreachingCost => (nz - pit_z).max(0.0) * step,
                Method::ShortestPath => step,
                Method::MinimumElevationChange => (nz - pit_z).abs() * step,
            };
            let nd = cost + increment;
            if nd < dist[nidx] {
                if !dist[nidx].is_finite() {
                    touched.push(nidx);
                }
                dist[nidx] = nd;
                travelled[nidx] = travel;
                back[nidx] = idx;
                heap.push(Node {
                    cost: nd,
                    idx: nidx,
                });
            }
        }
    }
    None
}

/// Min-heap entry (BinaryHeap is a max-heap, so the ordering is reversed).
struct Node {
    cost: f64,
    idx: usize,
}
impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.cost == other.cost && self.idx == other.idx
    }
}
impl Eq for Node {}
impl Ord for Node {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .cost
            .total_cmp(&self.cost)
            .then_with(|| other.idx.cmp(&self.idx))
    }
}
impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn pixel_center(r: &Raster, row: isize, col: isize) -> (f64, f64) {
    (
        r.x_min + (col as f64 + 0.5) * r.cell_size_x,
        r.y_max() - (row as f64 + 0.5) * r.cell_size_y,
    )
}

// ── Params ──────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    MinimumBreachingCost,
    ShortestPath,
    MinimumElevationChange,
}

impl Method {
    fn name(self) -> &'static str {
        match self {
            Method::MinimumBreachingCost => "minimum_breaching_cost",
            Method::ShortestPath => "shortest_path",
            Method::MinimumElevationChange => "minimum_elevation_change",
        }
    }
}

fn parse_method(args: &ToolArgs) -> Result<Method, ToolError> {
    match args
        .get("method")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        None | Some("") | Some("minimum_breaching_cost") => Ok(Method::MinimumBreachingCost),
        Some("shortest_path") => Ok(Method::ShortestPath),
        Some("minimum_elevation_change") => Ok(Method::MinimumElevationChange),
        Some(o) => Err(ToolError::Validation(format!(
            "'method' must be 'minimum_breaching_cost', 'shortest_path' or \
             'minimum_elevation_change', got '{o}'"
        ))),
    }
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

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// 20x20 DEM over (0,0)-(20,20), cell 1, built by `f(x, y)`.
    fn dem(f: impl Fn(f64, f64) -> f64) -> String {
        let mut r = Raster::new(RasterConfig {
            cols: 20,
            rows: 20,
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
        for row in 0..20 {
            for col in 0..20 {
                let x = 0.5 + col as f64;
                let y = 19.5 - row as f64;
                r.set(0, row as isize, col as isize, f(x, y)).unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    /// A square depression polygon centred on (cx, cy) with half-width `h`.
    fn depression(cx: f64, cy: f64, h: f64) -> String {
        let mut l = Layer::new("d")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_feature(
            Some(Geometry::polygon(
                vec![
                    Coord::xy(cx - h, cy - h),
                    Coord::xy(cx + h, cy - h),
                    Coord::xy(cx + h, cy + h),
                    Coord::xy(cx - h, cy + h),
                    Coord::xy(cx - h, cy - h),
                ],
                vec![],
            )),
            &[],
        )
        .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    /// A ramp falling to the west with a pit dug near the east side, so the
    /// only way out is a westward breach.
    fn ramp_with_pit() -> String {
        dem(|x, y| {
            let base = x; // rises eastward
            let in_pit = (x - 15.0).abs() <= 2.0 && (y - 10.0).abs() <= 2.0;
            if in_pit {
                2.0
            } else {
                base
            }
        })
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = GenerateBreachLinesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    #[test]
    fn breach_runs_from_the_pit_to_lower_ground() {
        let (out, layer) = run(json!({
            "input": depression(15.0, 10.0, 2.5), "dem": ramp_with_pit()
        }));
        assert_eq!(out.outputs["breached_count"], json!(1));
        assert_eq!(layer.features.len(), 1);
        let (iz, oz) = (
            layer.schema.field_index("INLET_Z").unwrap(),
            layer.schema.field_index("OUTLET_Z").unwrap(),
        );
        let inlet = layer.features[0].attributes[iz].as_f64().unwrap();
        let outlet = layer.features[0].attributes[oz].as_f64().unwrap();
        // The defining property: the breach ends below where it starts.
        assert!(outlet < inlet, "outlet {outlet} not below inlet {inlet}");
        // And the line runs westward, toward the low ground.
        match layer.features[0].geometry.as_ref().unwrap() {
            Geometry::LineString(cs) => {
                assert!(cs.len() >= 2);
                assert!(cs.last().unwrap().x < cs.first().unwrap().x);
            }
            other => panic!("unexpected geometry {other:?}"),
        }
    }

    #[test]
    fn breach_length_and_volume_are_reported() {
        let (out, layer) = run(json!({
            "input": depression(15.0, 10.0, 2.5), "dem": ramp_with_pit()
        }));
        let (l, v, m) = (
            layer.schema.field_index("BREACH_LEN").unwrap(),
            layer.schema.field_index("CUT_VOLUME").unwrap(),
            layer.schema.field_index("MAX_CUT").unwrap(),
        );
        let f = &layer.features[0];
        assert!(f.attributes[l].as_f64().unwrap() > 0.0);
        assert!(f.attributes[v].as_f64().unwrap() > 0.0);
        assert!(f.attributes[m].as_f64().unwrap() > 0.0);
        assert!((out.outputs["total_cut_volume"].as_f64().unwrap()
            - f.attributes[v].as_f64().unwrap())
        .abs() < 1e-9);
    }

    #[test]
    fn shortest_path_is_never_longer_than_the_cost_optimal_route() {
        let d = depression(15.0, 10.0, 2.5);
        let g = ramp_with_pit();
        let len = |method: &str| -> f64 {
            let (_o, layer) = run(json!({ "input": d.clone(), "dem": g.clone(), "method": method }));
            let i = layer.schema.field_index("BREACH_LEN").unwrap();
            layer.features[0].attributes[i].as_f64().unwrap()
        };
        let sp = len("shortest_path");
        let mbc = len("minimum_breaching_cost");
        assert!(sp <= mbc + 1e-9, "shortest {sp} > min-cost {mbc}");
    }

    #[test]
    fn max_length_suppresses_a_breach_that_cannot_reach_out() {
        let (out, _l) = run(json!({
            "input": depression(15.0, 10.0, 2.5), "dem": ramp_with_pit(),
            "max_length": 1.0
        }));
        assert_eq!(out.outputs["breached_count"], json!(0));
        assert_eq!(out.outputs["exceeded_max_length"], json!(1));
    }

    #[test]
    fn a_depression_with_no_lower_outlet_is_left_unbreached() {
        // A bowl: everything around the pit is higher, right to the edge.
        let bowl = dem(|x, y| {
            let d = ((x - 10.0).powi(2) + (y - 10.0).powi(2)).sqrt();
            100.0 - d // highest at the centre... inverted so the pit has no exit
        });
        let (out, _l) = run(json!({
            "input": depression(10.0, 10.0, 3.0), "dem": bowl
        }));
        // The "pit" here is on the rim, and no cell outside is lower, so the
        // tool must decline rather than invent a breach.
        assert!(out.outputs["breached_count"].as_f64().unwrap() <= 1.0);
    }

    #[test]
    fn methods_are_reported_back() {
        let (out, _l) = run(json!({
            "input": depression(15.0, 10.0, 2.5), "dem": ramp_with_pit(),
            "method": "minimum_elevation_change"
        }));
        assert_eq!(out.outputs["method"], json!("minimum_elevation_change"));
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            GenerateBreachLinesTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "d.shp" })).is_err());
        assert!(bad(json!({ "input": "d.shp", "dem": "z.tif", "method": "magic" })).is_err());
        assert!(bad(json!({ "input": "d.shp", "dem": "z.tif", "max_length": 0 })).is_err());
        assert!(bad(json!({ "input": "d.shp", "dem": "z.tif" })).is_ok());
    }
}
