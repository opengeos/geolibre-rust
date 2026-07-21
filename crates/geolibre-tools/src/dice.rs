//! GeoLibre tool: split pathologically large geometries into pieces below a
//! vertex limit, using an adaptive grid.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Dice* (Data Management). Million-vertex
//! geometries (dissolved coastlines, land cover) choke boolean overlay and tile
//! clipping. The nearest tools don't solve this: `subdivide_polygon` divides by
//! target *area* (not vertex count) and handles only polygons; `split_with_lines`
//! needs a user-supplied cutter. Dice is the standard safety valve before
//! `vector_to_pmtiles` or heavy overlay.
//!
//! Each feature with more than `vertex_limit` vertices has its bounding box
//! recursively quartered; the geometry is intersected with each quadrant
//! (polygon parts via `geo`'s `BooleanOps`, line parts via rectangle clipping)
//! and the quadrant recurses until every piece is under the limit. Features
//! already under the limit pass through untouched. Attributes are copied to
//! every output piece.

use std::collections::BTreeMap;

use geo::{BooleanOps, Coord as GeoCoord, LineString, MultiPolygon, Polygon, Rect};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Hard recursion cap so a degenerate geometry can't loop forever.
const MAX_DEPTH: u32 = 24;

pub struct DiceTool;

impl Tool for DiceTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "dice",
            display_name: "Dice",
            summary: "Split polygons or polylines with more than a vertex limit into a grid of smaller pieces (like ArcGIS Dice), so downstream overlay and tiling stay fast — an adaptive quadtree intersects the geometry with each quadrant until every piece is under the limit.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon or polyline layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output layer of diced pieces. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "vertex_limit",
                    description: "Maximum vertices per output feature (default 10000). Features at or under this pass through unchanged.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args
            .get("input")
            .and_then(Value::as_str)
            .map(str::trim)
            .unwrap_or("")
            .is_empty()
        {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        parse_vertex_limit(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = args
            .get("input")
            .and_then(Value::as_str)
            .filter(|s| !s.trim().is_empty())
            .ok_or_else(|| {
                ToolError::Validation("missing required parameter 'input'".to_string())
            })?;
        let output = parse_optional_str(args, "output")?;
        let limit = parse_vertex_limit(args)?;

        let layer = load_input_layer(input)?;

        let mut out = Layer::new(layer.name.clone());
        out.geom_type = layer.geom_type;
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for f in layer.schema.fields() {
            out.add_field(f.clone());
        }

        let mut diced = 0usize;
        let mut pieces_added = 0usize;
        for feat in &layer.features {
            let Some(geom) = feat.geometry.as_ref() else {
                out.push(feat.clone());
                pieces_added += 1;
                continue;
            };
            if vertex_count(geom) <= limit {
                out.push(feat.clone());
                pieces_added += 1;
                continue;
            }
            diced += 1;
            let parts = dice_geometry(geom, limit);
            for g in parts {
                out.push(Feature {
                    fid: 0,
                    geometry: Some(g),
                    attributes: feat.attributes.clone(),
                });
                pieces_added += 1;
            }
        }

        ctx.progress.info(&format!(
            "{diced} feature(s) diced into {pieces_added} piece(s)"
        ));

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("features_diced".to_string(), json!(diced));
        outputs.insert("output_features".to_string(), json!(pieces_added));
        Ok(ToolRunResult { outputs })
    }
}

/// Total vertex count across all rings/parts of a geometry.
fn vertex_count(geom: &Geometry) -> usize {
    match geom {
        Geometry::Point(_) => 1,
        Geometry::MultiPoint(cs) => cs.len(),
        Geometry::LineString(cs) => cs.len(),
        Geometry::MultiLineString(ls) => ls.iter().map(|l| l.len()).sum(),
        Geometry::Polygon {
            exterior,
            interiors,
        } => exterior.len() + interiors.iter().map(|r| r.len()).sum::<usize>(),
        Geometry::MultiPolygon(parts) => parts
            .iter()
            .map(|(e, ints)| e.len() + ints.iter().map(|r| r.len()).sum::<usize>())
            .sum(),
        Geometry::GeometryCollection(gs) => gs.iter().map(vertex_count).sum(),
    }
}

fn dice_geometry(geom: &Geometry, limit: usize) -> Vec<Geometry> {
    let Some(bb) = geom.bbox() else {
        return vec![geom.clone()];
    };
    let rect = (bb.min_x, bb.min_y, bb.max_x, bb.max_y);
    if let Some(mp) = to_multipolygon(geom) {
        return dice_polygon(&mp, rect, limit, 0)
            .into_iter()
            .map(|p| multipolygon_to_geometry(&p))
            .collect();
    }
    if let Some(lines) = to_lines(geom) {
        // Each returned piece is one contiguous polyline.
        return dice_lines(&lines, rect, limit, 0)
            .into_iter()
            .filter(|l| l.len() >= 2)
            .map(|l| lines_to_geometry(vec![l]))
            .collect();
    }
    vec![geom.clone()]
}

// ── Polygon dicing ──────────────────────────────────────────────────────────

fn dice_polygon(
    mp: &MultiPolygon,
    (xmin, ymin, xmax, ymax): (f64, f64, f64, f64),
    limit: usize,
    depth: u32,
) -> Vec<MultiPolygon> {
    let verts = mp_vertices(mp);
    if verts <= limit || depth >= MAX_DEPTH || (xmax - xmin) <= 0.0 || (ymax - ymin) <= 0.0 {
        return if mp.0.is_empty() {
            vec![]
        } else {
            vec![mp.clone()]
        };
    }
    let (mx, my) = ((xmin + xmax) / 2.0, (ymin + ymax) / 2.0);
    let quads = [
        (xmin, ymin, mx, my),
        (mx, ymin, xmax, my),
        (xmin, my, mx, ymax),
        (mx, my, xmax, ymax),
    ];
    let mut out = Vec::new();
    for q in quads {
        let clip = rect_polygon(q);
        let piece = mp.intersection(&clip);
        if piece.0.is_empty() {
            continue;
        }
        out.extend(dice_polygon(&piece, q, limit, depth + 1));
    }
    out
}

fn mp_vertices(mp: &MultiPolygon) -> usize {
    mp.0.iter()
        .map(|p| p.exterior().0.len() + p.interiors().iter().map(|r| r.0.len()).sum::<usize>())
        .sum()
}

fn rect_polygon((xmin, ymin, xmax, ymax): (f64, f64, f64, f64)) -> MultiPolygon {
    let rect = Rect::new(GeoCoord { x: xmin, y: ymin }, GeoCoord { x: xmax, y: ymax });
    MultiPolygon(vec![rect.to_polygon()])
}

// ── Line dicing ─────────────────────────────────────────────────────────────

type Line = Vec<(f64, f64)>;

fn dice_lines(
    lines: &[Line],
    (xmin, ymin, xmax, ymax): (f64, f64, f64, f64),
    limit: usize,
    depth: u32,
) -> Vec<Line> {
    let verts: usize = lines.iter().map(|l| l.len()).sum();
    if verts <= limit || depth >= MAX_DEPTH || (xmax - xmin) <= 0.0 || (ymax - ymin) <= 0.0 {
        return lines.to_vec();
    }
    let (mx, my) = ((xmin + xmax) / 2.0, (ymin + ymax) / 2.0);
    let quads = [
        (xmin, ymin, mx, my),
        (mx, ymin, xmax, my),
        (xmin, my, mx, ymax),
        (mx, my, xmax, ymax),
    ];
    let mut out = Vec::new();
    for q in quads {
        let clipped: Vec<Line> = lines
            .iter()
            .flat_map(|l| clip_line_to_rect(l, q))
            .filter(|l| l.len() >= 2)
            .collect();
        if clipped.is_empty() {
            continue;
        }
        out.extend(dice_lines(&clipped, q, limit, depth + 1));
    }
    out
}

/// Clips a polyline to a rectangle, returning the contiguous inside pieces.
fn clip_line_to_rect(line: &Line, (xmin, ymin, xmax, ymax): (f64, f64, f64, f64)) -> Vec<Line> {
    let mut pieces: Vec<Line> = Vec::new();
    let mut cur: Line = Vec::new();
    for w in line.windows(2) {
        if let Some(((ax, ay), (bx, by))) = liang_barsky(w[0], w[1], xmin, ymin, xmax, ymax) {
            if cur.is_empty() {
                cur.push((ax, ay));
            } else if cur.last() != Some(&(ax, ay)) {
                // Discontinuity: the clipped segment starts away from the last
                // point -> begin a new piece.
                pieces.push(std::mem::take(&mut cur));
                cur.push((ax, ay));
            }
            cur.push((bx, by));
        } else if !cur.is_empty() {
            pieces.push(std::mem::take(&mut cur));
        }
    }
    if cur.len() >= 2 {
        pieces.push(cur);
    }
    pieces
}

/// Liang–Barsky segment clip; returns the clipped endpoints inside the rect.
#[allow(clippy::too_many_arguments)]
fn liang_barsky(
    (x0, y0): (f64, f64),
    (x1, y1): (f64, f64),
    xmin: f64,
    ymin: f64,
    xmax: f64,
    ymax: f64,
) -> Option<((f64, f64), (f64, f64))> {
    let dx = x1 - x0;
    let dy = y1 - y0;
    let mut t0 = 0.0f64;
    let mut t1 = 1.0f64;
    let checks = [
        (-dx, x0 - xmin),
        (dx, xmax - x0),
        (-dy, y0 - ymin),
        (dy, ymax - y0),
    ];
    for (p, q) in checks {
        if p == 0.0 {
            if q < 0.0 {
                return None; // parallel and outside
            }
        } else {
            let r = q / p;
            if p < 0.0 {
                if r > t1 {
                    return None;
                }
                if r > t0 {
                    t0 = r;
                }
            } else {
                if r < t0 {
                    return None;
                }
                if r < t1 {
                    t1 = r;
                }
            }
        }
    }
    Some(((x0 + t0 * dx, y0 + t0 * dy), (x0 + t1 * dx, y0 + t1 * dy)))
}

// ── Geometry <-> geo conversions ────────────────────────────────────────────

fn to_multipolygon(geom: &Geometry) -> Option<MultiPolygon> {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => Some(MultiPolygon(vec![rings_to_polygon(exterior, interiors)])),
        Geometry::MultiPolygon(parts) => Some(MultiPolygon(
            parts.iter().map(|(e, i)| rings_to_polygon(e, i)).collect(),
        )),
        _ => None,
    }
}

fn rings_to_polygon(exterior: &Ring, interiors: &[Ring]) -> Polygon {
    Polygon::new(
        ring_to_linestring(exterior),
        interiors.iter().map(ring_to_linestring).collect(),
    )
}

fn ring_to_linestring(ring: &Ring) -> LineString {
    LineString::new(
        ring.coords()
            .iter()
            .map(|c| GeoCoord { x: c.x, y: c.y })
            .collect(),
    )
}

fn multipolygon_to_geometry(mp: &MultiPolygon) -> Geometry {
    if mp.0.len() == 1 {
        let (exterior, interiors) = polygon_to_rings(&mp.0[0]);
        Geometry::Polygon {
            exterior,
            interiors,
        }
    } else {
        Geometry::MultiPolygon(mp.0.iter().map(polygon_to_rings).collect())
    }
}

fn polygon_to_rings(poly: &Polygon) -> (Ring, Vec<Ring>) {
    (
        linestring_to_ring(poly.exterior()),
        poly.interiors().iter().map(linestring_to_ring).collect(),
    )
}

fn linestring_to_ring(ls: &LineString) -> Ring {
    let mut coords: Vec<Coord> = ls.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
    if coords.len() >= 2 && coords.first() == coords.last() {
        coords.pop();
    }
    Ring::new(coords)
}

fn to_lines(geom: &Geometry) -> Option<Vec<Line>> {
    match geom {
        Geometry::LineString(cs) => Some(vec![cs.iter().map(|c| (c.x, c.y)).collect()]),
        Geometry::MultiLineString(ls) => Some(
            ls.iter()
                .map(|l| l.iter().map(|c| (c.x, c.y)).collect())
                .collect(),
        ),
        _ => None,
    }
}

fn lines_to_geometry(lines: Vec<Line>) -> Geometry {
    if lines.len() == 1 {
        Geometry::LineString(
            lines
                .into_iter()
                .next()
                .unwrap()
                .iter()
                .map(|&(x, y)| Coord::xy(x, y))
                .collect(),
        )
    } else {
        Geometry::MultiLineString(
            lines
                .into_iter()
                .map(|l| l.iter().map(|&(x, y)| Coord::xy(x, y)).collect())
                .collect(),
        )
    }
}

// ── Parameters ──────────────────────────────────────────────────────────────

fn parse_vertex_limit(args: &ToolArgs) -> Result<usize, ToolError> {
    let limit = match args.get("vertex_limit") {
        None | Some(Value::Null) => 10_000,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(10_000) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 10_000,
        Some(Value::String(s)) => s.trim().parse::<usize>().map_err(|_| {
            ToolError::Validation("'vertex_limit' must be a positive integer".into())
        })?,
        _ => {
            return Err(ToolError::Validation(
                "'vertex_limit' must be an integer".into(),
            ))
        }
    };
    if limit < 4 {
        return Err(ToolError::Validation(
            "'vertex_limit' must be at least 4".to_string(),
        ));
    }
    Ok(limit)
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

    /// A big rectangular polygon with many densified boundary vertices gets
    /// split into pieces, each under the vertex limit, covering the same area.
    #[test]
    fn dices_a_dense_polygon() {
        // Build a square with many collinear points on each edge (high vertex
        // count, simple shape).
        let n = 200usize;
        let mut ring: Vec<Coord> = Vec::new();
        for i in 0..n {
            ring.push(Coord::xy(100.0 * i as f64 / n as f64, 0.0));
        }
        for i in 0..n {
            ring.push(Coord::xy(100.0, 100.0 * i as f64 / n as f64));
        }
        for i in 0..n {
            ring.push(Coord::xy(100.0 - 100.0 * i as f64 / n as f64, 100.0));
        }
        for i in 0..n {
            ring.push(Coord::xy(0.0, 100.0 - 100.0 * i as f64 / n as f64));
        }
        let geom = Geometry::Polygon {
            exterior: Ring::new(ring),
            interiors: vec![],
        };
        let mut l = Layer::new("big")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("id", FieldType::Integer));
        l.add_feature(Some(geom), &[("id", 1i64.into())]).unwrap();
        let id = memory_store::put_vector(l);
        let path = memory_store::make_vector_memory_path(&id);

        let args: ToolArgs =
            serde_json::from_value(json!({ "input": path, "vertex_limit": 100 })).unwrap();
        let out = DiceTool.run(&args, &ctx()).unwrap();
        let res = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();

        assert!(
            res.features.len() > 1,
            "dense polygon should split into pieces"
        );
        for f in &res.features {
            let vc = vertex_count(f.geometry.as_ref().unwrap());
            assert!(
                vc <= 100,
                "each piece must be under the vertex limit, got {vc}"
            );
            // attributes preserved
            assert_eq!(f.attributes[0].as_i64(), Some(1));
        }
        // Total area is conserved by the intersection tiling.
        let total: f64 = res
            .features
            .iter()
            .filter_map(|f| to_multipolygon(f.geometry.as_ref().unwrap()))
            .map(|mp| geo::Area::unsigned_area(&mp))
            .sum();
        assert!(
            (total - 10_000.0).abs() < 1.0,
            "diced area must equal the original 100x100"
        );
    }

    /// A small feature under the limit passes through unchanged.
    #[test]
    fn passes_small_features_through() {
        let geom = Geometry::Polygon {
            exterior: Ring::new(vec![
                Coord::xy(0.0, 0.0),
                Coord::xy(1.0, 0.0),
                Coord::xy(1.0, 1.0),
                Coord::xy(0.0, 1.0),
            ]),
            interiors: vec![],
        };
        let mut l = Layer::new("s")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("id", FieldType::Integer));
        l.add_feature(Some(geom), &[("id", 7i64.into())]).unwrap();
        let id = memory_store::put_vector(l);
        let path = memory_store::make_vector_memory_path(&id);
        let args: ToolArgs = serde_json::from_value(json!({ "input": path })).unwrap();
        let out = DiceTool.run(&args, &ctx()).unwrap();
        assert_eq!(out.outputs["features_diced"], json!(0));
        assert_eq!(out.outputs["output_features"], json!(1));
    }

    /// A long dense polyline gets diced and stays within the limit.
    #[test]
    fn dices_a_polyline() {
        let coords: Vec<Coord> = (0..500)
            .map(|i| Coord::xy(i as f64, (i as f64 * 0.1).sin() * 10.0))
            .collect();
        let geom = Geometry::LineString(coords);
        let mut l = Layer::new("ln")
            .with_geom_type(GeometryType::LineString)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("id", FieldType::Integer));
        l.add_feature(Some(geom), &[("id", 3i64.into())]).unwrap();
        let id = memory_store::put_vector(l);
        let path = memory_store::make_vector_memory_path(&id);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": path, "vertex_limit": 50 })).unwrap();
        let out = DiceTool.run(&args, &ctx()).unwrap();
        let res = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        assert!(res.features.len() > 1);
        for f in &res.features {
            assert!(vertex_count(f.geometry.as_ref().unwrap()) <= 50);
        }
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            DiceTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "vertex_limit": 2 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "vertex_limit": 5000 })).is_ok());
    }
}
