//! GeoLibre tool: least-cost network connecting multiple sites over a cost
//! surface.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Cost Connectivity* / *Optimal Region
//! Connections* (Spatial Analyst). The bundled suite computes cost distance from
//! one source set (`cost_distance`, `cost_pathway`) and GeoLibre ships
//! `path_distance` and `corridor`, but nothing produces the least-cost *network*
//! among N sites — the standard output for wildlife-corridor design and
//! infrastructure planning (connect all reserves/sites at minimum total cost).
//!
//! Every source point seeds a multi-source Dijkstra over the cost grid, so each
//! cell learns its nearest source (allocation), accumulated cost, and back-link.
//! Where two sources' allocation regions touch, the crossing cell with the
//! smallest summed cost defines the least-cost path between that pair. The pair
//! costs form a graph; `connections=mst` returns its minimum spanning tree
//! (connect all sites at minimum total cost), `all_neighbors` returns every
//! adjacent-region path. Paths are emitted as polylines with from/to ids and
//! accumulated cost.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::common::load_input_raster;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CostConnectivityTool;

impl Tool for CostConnectivityTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "cost_connectivity",
            display_name: "Cost Connectivity",
            summary: "Least-cost network connecting multiple sites over a cost surface (like ArcGIS Cost Connectivity / Optimal Region Connections): multi-source Dijkstra allocation, pairwise least-cost paths between adjacent regions, and their minimum spanning tree — the least-cost network the bundled single-source cost_distance and path_distance can't build.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "sources",
                    description: "Point layer of sites/regions to connect.",
                    required: true,
                },
                ToolParamSpec {
                    name: "cost",
                    description: "Cost/friction raster (per-cell traversal cost; no-data cells are barriers).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polyline layer of least-cost paths. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "connections",
                    description: "'mst' (minimum spanning tree connecting all sites; default) or 'all_neighbors' (every adjacent-region path).",
                    required: false,
                },
                ToolParamSpec {
                    name: "id_field",
                    description: "Optional field to label sources (default: 1-based index).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band of the cost raster (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "sources")?;
        require_str(args, "cost")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let sources_path = require_str(args, "sources")?;
        let cost_path = require_str(args, "cost")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let cost = load_input_raster(cost_path)?;
        if prm.band < 0 || prm.band as usize >= cost.bands {
            return Err(ToolError::Validation(format!(
                "band {} out of range",
                prm.band + 1
            )));
        }
        let layer = load_input_layer(sources_path)?;
        let id_idx = match &prm.id_field {
            Some(f) => Some(
                layer
                    .schema
                    .field_index(f)
                    .ok_or_else(|| ToolError::Validation(format!("id_field '{f}' not found")))?,
            ),
            None => None,
        };

        let rows = cost.rows;
        let cols = cost.cols;
        let nodata = cost.nodata;
        let cx = cost.cell_size_x;
        let cy = cost.cell_size_y;
        let x0 = cost.x_min;
        let ymax = cost.y_min + rows as f64 * cy;

        // Read cost into a flat grid.
        let mut cgrid = vec![f64::NAN; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let v = cost.get(prm.band, r as isize, c as isize);
                cgrid[r * cols + c] = if v != nodata && v.is_finite() && v >= 0.0 {
                    v
                } else {
                    f64::NAN
                };
            }
        }

        // Map source points to cells.
        let mut sources: Vec<Source> = Vec::new();
        for (fi, feat) in layer.features.iter().enumerate() {
            let Some((wx, wy)) = feat.geometry.as_ref().and_then(point_xy) else {
                continue;
            };
            let col = (((wx - x0) / cx).floor() as isize).clamp(0, cols as isize - 1) as usize;
            let row = (((ymax - wy) / cy).floor() as isize).clamp(0, rows as isize - 1) as usize;
            let idx = row * cols + col;
            if cgrid[idx].is_nan() {
                continue; // source on a barrier
            }
            let label = match id_idx {
                Some(i) => feat.attributes.get(i).map(value_string).unwrap_or_default(),
                None => (fi + 1).to_string(),
            };
            sources.push(Source { cell: idx, label });
        }
        if sources.len() < 2 {
            return Err(ToolError::Execution(
                "need at least 2 valid sources on non-barrier cells".to_string(),
            ));
        }

        ctx.progress.info(&format!(
            "{} source(s); multi-source Dijkstra",
            sources.len()
        ));

        // Multi-source Dijkstra: allocation, accumulated cost, back-link.
        let n = rows * cols;
        let mut acc = vec![f64::INFINITY; n];
        let mut region = vec![usize::MAX; n];
        let mut back = vec![usize::MAX; n];
        let mut heap: BinaryHeap<State> = BinaryHeap::new();
        for (si, s) in sources.iter().enumerate() {
            acc[s.cell] = 0.0;
            region[s.cell] = si;
            heap.push(State {
                cost: 0.0,
                cell: s.cell,
            });
        }
        let neigh: [(isize, isize); 8] = [
            (-1, -1),
            (-1, 0),
            (-1, 1),
            (0, -1),
            (0, 1),
            (1, -1),
            (1, 0),
            (1, 1),
        ];
        while let Some(State { cost: cc, cell }) = heap.pop() {
            if cc > acc[cell] {
                continue;
            }
            let r = (cell / cols) as isize;
            let c = (cell % cols) as isize;
            for (dr, dc) in neigh {
                let nr = r + dr;
                let nc = c + dc;
                if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                    continue;
                }
                let nidx = nr as usize * cols + nc as usize;
                if cgrid[nidx].is_nan() {
                    continue;
                }
                let dist = if dr != 0 && dc != 0 {
                    (cx * cx + cy * cy).sqrt()
                } else if dr != 0 {
                    cy
                } else {
                    cx
                };
                let step = 0.5 * (cgrid[cell] + cgrid[nidx]) * dist;
                let nc_cost = cc + step;
                if nc_cost < acc[nidx] {
                    acc[nidx] = nc_cost;
                    region[nidx] = region[cell];
                    back[nidx] = cell;
                    heap.push(State {
                        cost: nc_cost,
                        cell: nidx,
                    });
                }
            }
        }

        // Find the min-cost crossing per region pair.
        let mut best_pair: BTreeMap<(usize, usize), Crossing> = BTreeMap::new();
        for r in 0..rows {
            for c in 0..cols {
                let a = r * cols + c;
                if region[a] == usize::MAX {
                    continue;
                }
                for (dr, dc) in [(0isize, 1isize), (1, 0), (1, 1), (1, -1)] {
                    let nr = r as isize + dr;
                    let nc = c as isize + dc;
                    if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                        continue;
                    }
                    let b = nr as usize * cols + nc as usize;
                    if region[b] == usize::MAX || region[a] == region[b] {
                        continue;
                    }
                    let dist = if dr != 0 && dc != 0 {
                        (cx * cx + cy * cy).sqrt()
                    } else if dr != 0 {
                        cy
                    } else {
                        cx
                    };
                    let link = 0.5 * (cgrid[a] + cgrid[b]) * dist;
                    let total = acc[a] + acc[b] + link;
                    let key = order2(region[a], region[b]);
                    let e = best_pair.entry(key).or_insert(Crossing {
                        total: f64::INFINITY,
                        a,
                        b,
                    });
                    if total < e.total {
                        *e = Crossing { total, a, b };
                    }
                }
            }
        }

        // Choose the edge set: MST (Kruskal) or all neighbours.
        let edges: Vec<((usize, usize), Crossing)> = {
            let mut v: Vec<_> = best_pair.into_iter().collect();
            v.sort_by(|x, y| x.1.total.total_cmp(&y.1.total));
            if prm.mst {
                let mut uf = UnionFind::new(sources.len());
                v.into_iter()
                    .filter(|((ra, rb), _)| uf.union(*ra, *rb))
                    .collect()
            } else {
                v
            }
        };

        // Build the output polylines.
        let mut out = Layer::new("cost_paths").with_geom_type(GeometryType::LineString);
        if let Some(e) = cost.crs.epsg {
            out = out.with_crs_epsg(e);
        }
        out.add_field(FieldDef::new("from_id", FieldType::Text));
        out.add_field(FieldDef::new("to_id", FieldType::Text));
        out.add_field(FieldDef::new("cost", FieldType::Float));

        let cell_xy = |cell: usize| -> Coord {
            let r = cell / cols;
            let c = cell % cols;
            Coord::xy(x0 + (c as f64 + 0.5) * cx, ymax - (r as f64 + 0.5) * cy)
        };

        let mut total_cost = 0.0;
        for ((ra, rb), cr) in &edges {
            // Path: source(ra) <- ... <- a  ++  b -> ... -> source(rb).
            let mut left = trace(cr.a, &back);
            left.reverse(); // now source(a) .. a
            let right = trace(cr.b, &back); // b .. source(b)
            let mut coords: Vec<Coord> = left.iter().map(|&cl| cell_xy(cl)).collect();
            coords.extend(right.iter().map(|&cl| cell_xy(cl)));
            if coords.len() < 2 {
                continue;
            }
            total_cost += cr.total;
            out.push(Feature {
                fid: 0,
                geometry: Some(Geometry::line_string(coords)),
                attributes: vec![
                    FieldValue::Text(sources[*ra].label.clone()),
                    FieldValue::Text(sources[*rb].label.clone()),
                    FieldValue::Float(cr.total),
                ],
            });
        }
        let path_count = out.features.len();

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("source_count".to_string(), json!(sources.len()));
        outputs.insert("path_count".to_string(), json!(path_count));
        outputs.insert("total_cost".to_string(), json!(total_cost));
        Ok(ToolRunResult { outputs })
    }
}

/// Traces back-links from `cell` to its source, returning [cell .. source].
fn trace(mut cell: usize, back: &[usize]) -> Vec<usize> {
    let mut path = vec![cell];
    while back[cell] != usize::MAX {
        cell = back[cell];
        path.push(cell);
    }
    path
}

struct Source {
    cell: usize,
    label: String,
}

#[derive(Clone, Copy)]
struct Crossing {
    total: f64,
    a: usize,
    b: usize,
}

/// Dijkstra heap state (min-heap via reversed ordering).
struct State {
    cost: f64,
    cell: usize,
}
impl PartialEq for State {
    fn eq(&self, o: &Self) -> bool {
        self.cost == o.cost
    }
}
impl Eq for State {}
impl PartialOrd for State {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for State {
    fn cmp(&self, o: &Self) -> Ordering {
        o.cost.total_cmp(&self.cost).then(self.cell.cmp(&o.cell))
    }
}

struct UnionFind {
    parent: Vec<usize>,
}
impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
        }
    }
    fn find(&mut self, x: usize) -> usize {
        let mut r = x;
        while self.parent[r] != r {
            r = self.parent[r];
        }
        let mut c = x;
        while self.parent[c] != r {
            let next = self.parent[c];
            self.parent[c] = r;
            c = next;
        }
        r
    }
    fn union(&mut self, a: usize, b: usize) -> bool {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            false
        } else {
            self.parent[ra] = rb;
            true
        }
    }
}

fn order2(a: usize, b: usize) -> (usize, usize) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        _ => None,
    }
}

fn value_string(fv: &FieldValue) -> String {
    if let Some(i) = fv.as_i64() {
        i.to_string()
    } else if let Some(f) = fv.as_f64() {
        format!("{f}")
    } else {
        fv.as_str().unwrap_or("").to_string()
    }
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

struct Params {
    mst: bool,
    id_field: Option<String>,
    band: isize,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let mst = match args
        .get("connections")
        .and_then(Value::as_str)
        .map(str::trim)
    {
        None | Some("") | Some("mst") => true,
        Some("all_neighbors") => false,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "'connections' must be 'mst' or 'all_neighbors', got '{o}'"
            )))
        }
    };
    let id_field = parse_optional_str(args, "id_field")?.map(String::from);
    let band_1based = args.get("band").and_then(Value::as_u64).unwrap_or(1).max(1);
    Ok(Params {
        mst,
        id_field,
        band: (band_1based - 1) as isize,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn cost_raster(cols: usize, rows: usize, data: Vec<f64>) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -1.0,
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

    fn point_layer(pts: &[(f64, f64)]) -> String {
        let mut l = Layer::new("s")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("id", FieldType::Integer));
        for (i, (x, y)) in pts.iter().enumerate() {
            l.add_feature(Some(Geometry::point(*x, *y)), &[("id", (i as i64).into())])
                .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    /// Three sites -> an MST has exactly 2 edges connecting all of them.
    #[test]
    fn mst_connects_all_sites() {
        // 10x10 uniform cost 1.
        let cost = cost_raster(10, 10, vec![1.0; 100]);
        // three points near three corners (world y is up; row0 is top).
        let src = point_layer(&[(0.5, 9.5), (9.5, 9.5), (5.0, 0.5)]);
        let args: ToolArgs = serde_json::from_value(json!({
            "sources": src, "cost": cost, "connections": "mst",
        }))
        .unwrap();
        let out = CostConnectivityTool.run(&args, &ctx()).unwrap();
        assert_eq!(out.outputs["source_count"], json!(3));
        assert_eq!(
            out.outputs["path_count"],
            json!(2),
            "MST of 3 nodes has 2 edges"
        );
    }

    /// A barrier of high cost makes the connecting path more expensive than the
    /// straight-line distance would suggest.
    #[test]
    fn barrier_raises_path_cost() {
        // 5x5: a vertical wall of high cost in the middle column.
        let mut data = vec![1.0; 25];
        for r in 0..5 {
            data[r * 5 + 2] = 100.0;
        }
        let cost = cost_raster(5, 5, data);
        let src = point_layer(&[(0.5, 2.5), (4.5, 2.5)]); // opposite sides of the wall
        let args: ToolArgs =
            serde_json::from_value(json!({ "sources": src, "cost": cost })).unwrap();
        let out = CostConnectivityTool.run(&args, &ctx()).unwrap();
        assert_eq!(out.outputs["path_count"], json!(1));
        assert!(
            out.outputs["total_cost"].as_f64().unwrap() > 50.0,
            "crossing the wall must be costly"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CostConnectivityTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "sources": "a.geojson" })).is_err());
        assert!(
            bad(json!({ "sources": "a.geojson", "cost": "c.tif", "connections": "star" })).is_err()
        );
        assert!(bad(json!({ "sources": "a.geojson", "cost": "c.tif" })).is_ok());
    }
}
