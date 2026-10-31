//! GeoLibre tool: contours at an explicit list of values.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Contour List* (Spatial Analyst /
//! 3D Analyst).
//!
//! ## Why the catalog needs it
//!
//! Most contour sets that matter are not evenly spaced. Nautical charts use the
//! IHO depth series (2, 5, 10, 20, 30, 50, 100 m); flood mapping wants the 1-,
//! 2- and 5-percent annual-chance stages; an aviation chart wants the minimum
//! safe altitudes; a soils map wants published class breaks. None of those is
//! an arithmetic progression.
//!
//! The bundled `contours_from_raster` takes only `interval` and `base`, so it
//! can generate an evenly-spaced series and nothing else. Reaching an irregular
//! set through it means one run per value followed by a merge, and even that
//! cannot produce a single layer whose features are labelled by level.
//! `percentile_contours` chooses its own levels from the distribution rather
//! than accepting given ones.
//!
//! ## Method
//!
//! Marching squares over the grid of cell centres. For each requested level,
//! every 2x2 block of centres is classified by which corners sit above it, and
//! the crossing points on the block's edges are interpolated linearly. The
//! resulting segments are chained end to end into polylines, so a contour comes
//! out as a few long lines rather than thousands of two-point stubs. The
//! saddle case (two opposite corners above the level) is resolved with the
//! block's mean, which is the standard disambiguation and keeps contours from
//! crossing each other.
//!
//! `geometry_type: polygon` instead returns the filled region at or above each
//! level, traced with the shared `polygonize` machinery.

use std::collections::BTreeMap;
use std::collections::HashMap;

use serde_json::{json, Map, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::args_common::{band_index, bool_or, choice_or, opt_positive_f64, req_str};
use crate::common::{load_input_raster, parse_optional_output};
use crate::geojson_geom::geometry_from_json;
use crate::polygonize::{polygonize_to_geojson, PolygonizeParams};
use crate::vector_common::write_or_store_layer;

pub struct ContourListTool;

impl Tool for ContourListTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "contour_list",
            display_name: "Contour List",
            summary: "Generates contours at an explicit list of values rather than a fixed interval, as chained polylines or filled polygons, each labelled with its level (ArcGIS Contour List). The bundled contours_from_raster accepts only 'interval' and 'base', so the irregular series that real charts use — IHO depth steps, flood recurrence stages, published class breaks — need one run per value and a merge, and still cannot come out as one layer labelled by level.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input surface raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "contour_values",
                    description: "Comma-separated list of levels to contour, in the raster's own z units.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output contour layer. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "geometry_type",
                    description: "'polyline' (default) for contour lines, or 'polygon' for the filled region at or above each level.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_length",
                    description: "Discard contour lines shorter than this, in map units, which removes the specks a noisy surface produces.",
                    required: false,
                },
                ToolParamSpec {
                    name: "close_rings",
                    description: "Repeat the first vertex on contours that close on themselves (default true).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to contour (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        parse_levels(args)?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input_path = req_str(args, "input")?.to_string();
        let levels = parse_levels(args)?;
        let prm = parse_params(args)?;
        let band = band_index(args, "band")?;
        let output = parse_optional_output(args, "output")?;

        let raster = load_input_raster(&input_path)?;
        let (rows, cols) = (raster.rows, raster.cols);
        if rows < 2 || cols < 2 {
            return Err(ToolError::Validation(
                "contouring needs a raster at least 2x2".to_string(),
            ));
        }

        let mut z = vec![f64::NAN; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let v = raster.get(band, r as isize, c as isize);
                if v != raster.nodata && v.is_finite() {
                    z[r * cols + c] = v;
                }
            }
        }

        ctx.progress.info(&format!(
            "{rows}x{cols}, {} level(s), {} output",
            levels.len(),
            if prm.polygons { "polygon" } else { "polyline" }
        ));

        let mut layer = Layer::new("contour_list");
        layer.geom_type = Some(if prm.polygons {
            GeometryType::Polygon
        } else {
            GeometryType::LineString
        });
        if let Some(e) = raster.crs.epsg {
            layer = layer.with_crs_epsg(e);
        }
        layer.add_field(FieldDef::new("id", FieldType::Integer));
        layer.add_field(FieldDef::new("contour", FieldType::Float));
        layer.add_field(FieldDef::new("length", FieldType::Float));
        layer.add_field(FieldDef::new("closed", FieldType::Boolean));

        // Cell centres, which is where the sampled values live.
        let y_max = raster.y_min + rows as f64 * raster.cell_size_y;
        let px = |c: f64| raster.x_min + (c + 0.5) * raster.cell_size_x;
        let py = |r: f64| y_max - (r + 0.5) * raster.cell_size_y;

        let mut fid = 0u64;
        let mut per_level: Vec<Value> = Vec::new();
        for (li, &level) in levels.iter().enumerate() {
            let mut emitted = 0usize;

            if prm.polygons {
                let mut labels = vec![0.0f64; rows * cols];
                for i in 0..rows * cols {
                    if z[i].is_finite() && z[i] >= level {
                        labels[i] = 1.0;
                    }
                }
                let props: HashMap<i64, Map<String, Value>> = HashMap::new();
                let geojson = polygonize_to_geojson(&PolygonizeParams {
                    labels: &labels,
                    rows,
                    cols,
                    x_min: raster.x_min,
                    y_max,
                    cell_size_x: raster.cell_size_x,
                    cell_size_y: raster.cell_size_y,
                    epsg: raster.crs.epsg,
                    props_by_id: &props,
                });
                let parsed: Value = serde_json::from_str(&geojson).map_err(|e| {
                    ToolError::Execution(format!("polygonize produced invalid GeoJSON: {e}"))
                })?;
                if let Some(feats) = parsed.get("features").and_then(Value::as_array) {
                    for f in feats {
                        let Some(geom) = f.get("geometry").and_then(geometry_from_json) else {
                            continue;
                        };
                        let mut feat =
                            Feature::with_geometry(fid, geom, layer.schema.len());
                        feat.set_by_index(0, FieldValue::Integer(fid as i64));
                        feat.set_by_index(1, FieldValue::Float(level));
                        feat.set_by_index(2, FieldValue::Float(0.0));
                        feat.set_by_index(3, FieldValue::Boolean(true));
                        layer.push(feat);
                        fid += 1;
                        emitted += 1;
                    }
                }
            } else {
                let segments = marching_squares(&z, rows, cols, level);
                for chain in chain_segments(segments) {
                    let coords: Vec<Coord> = chain
                        .iter()
                        .map(|&(c, r)| Coord::xy(px(c), py(r)))
                        .collect();
                    if coords.len() < 2 {
                        continue;
                    }
                    let closed = coords.len() > 2
                        && near(coords[0].x, coords[coords.len() - 1].x)
                        && near(coords[0].y, coords[coords.len() - 1].y);
                    let length = polyline_length(&coords);
                    if let Some(min) = prm.min_length {
                        if length < min {
                            continue;
                        }
                    }
                    let mut coords = coords;
                    if closed && prm.close_rings {
                        let first = coords[0].clone();
                        if !near(coords[coords.len() - 1].x, first.x)
                            || !near(coords[coords.len() - 1].y, first.y)
                        {
                            coords.push(first);
                        }
                    }

                    let mut feat = Feature::with_geometry(
                        fid,
                        Geometry::LineString(coords),
                        layer.schema.len(),
                    );
                    feat.set_by_index(0, FieldValue::Integer(fid as i64));
                    feat.set_by_index(1, FieldValue::Float(level));
                    feat.set_by_index(2, FieldValue::Float(length));
                    feat.set_by_index(3, FieldValue::Boolean(closed));
                    layer.push(feat);
                    fid += 1;
                    emitted += 1;
                }
            }

            per_level.push(json!({ "contour": level, "feature_count": emitted }));
            ctx.progress
                .progress((li as f64 + 1.0) / levels.len() as f64);
        }

        let feature_count = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("level_count".to_string(), json!(levels.len()));
        outputs.insert("levels".to_string(), Value::Array(per_level));
        Ok(ToolRunResult { outputs })
    }
}

/// A point on the cell-centre grid, in fractional `(col, row)`.
type Pt = (f64, f64);

/// How close two crossing points must be, in cell widths, to be treated as the
/// same vertex when chaining. Small enough never to merge distinct crossings,
/// large enough to absorb both float noise and the [`nudge_level`] split.
const WELD_TOLERANCE: f64 = 1e-5;

fn near(a: f64, b: f64) -> bool {
    (a - b).abs() <= WELD_TOLERANCE
}

/// Nudges a level off any value the surface actually attains.
///
/// A corner sitting *exactly* on the level is the classic marching-squares
/// degeneracy: the interpolated crossing collapses onto that corner, adjacent
/// blocks emit zero-length segments there, and the ring fragments. Shifting the
/// level by a hair — far below any plotting precision, and scaled to the
/// surface's own magnitude — removes the case instead of special-casing it
/// downstream. Integer surfaces contoured at integer levels hit it constantly.
///
/// The shift is **downward**, so a cell exactly on the level counts as inside.
/// Nudging upward instead carves those cells out of the enclosed region, and
/// every one that happens to lie on the contour's own path sprouts a spurious
/// micro-loop around itself — a cone contoured at a level its grid attains
/// exactly produced four of them. Downward also keeps this agreeing with the
/// polygon form, which selects cells with `value >= level`.
fn nudge_level(z: &[f64], level: f64) -> f64 {
    let scale = z
        .iter()
        .filter(|v| v.is_finite())
        .fold(0.0f64, |acc, v| acc.max(v.abs()))
        .max(level.abs())
        .max(1.0);
    let eps = scale * 1e-9;
    if z.iter().any(|&v| v.is_finite() && v == level) {
        level - eps
    } else {
        level
    }
}

/// Marching squares over the grid of cell centres, returning unordered
/// segments in fractional `(col, row)` space.
fn marching_squares(z: &[f64], rows: usize, cols: usize, level: f64) -> Vec<(Pt, Pt)> {
    let level = nudge_level(z, level);
    let mut out = Vec::new();
    for r in 0..rows - 1 {
        for c in 0..cols - 1 {
            // Corners of the block, clockwise from top-left.
            let v = [
                z[r * cols + c],
                z[r * cols + c + 1],
                z[(r + 1) * cols + c + 1],
                z[(r + 1) * cols + c],
            ];
            // A block touching no-data cannot be interpolated across.
            if v.iter().any(|x| !x.is_finite()) {
                continue;
            }
            let corner = [
                (c as f64, r as f64),
                (c as f64 + 1.0, r as f64),
                (c as f64 + 1.0, r as f64 + 1.0),
                (c as f64, r as f64 + 1.0),
            ];
            // Bit per corner at or above the level.
            let mut case = 0usize;
            for (k, &val) in v.iter().enumerate() {
                if val >= level {
                    case |= 1 << k;
                }
            }
            if case == 0 || case == 15 {
                continue;
            }

            // Crossing point on edge k (corner k -> corner k+1).
            let cross = |k: usize| -> Pt {
                let a = v[k];
                let b = v[(k + 1) % 4];
                let (pa, pb) = (corner[k], corner[(k + 1) % 4]);
                // The edge is only queried when the two ends straddle the
                // level, so the denominator cannot vanish.
                let t = ((level - a) / (b - a)).clamp(0.0, 1.0);
                (pa.0 + t * (pb.0 - pa.0), pa.1 + t * (pb.1 - pa.1))
            };

            // Edges 0=top, 1=right, 2=bottom, 3=left.
            //
            // Segments shorter than a millionth of a cell are discarded. Where
            // the surface touches the level exactly at a grid point, the
            // epsilon nudge above leaves a hair-thin sliver cutting that corner
            // — a few times 1e-8 of a cell across. It is a rounding artefact,
            // not geometry: it carries no cartographic information, and left in
            // it becomes a stray two-point feature that breaks the "one ring
            // per closed contour" the caller expects. No real contour segment
            // is anywhere near that short.
            let mut push = |e1: usize, e2: usize| {
                let (a, b) = (cross(e1), cross(e2));
                if (a.0 - b.0).hypot(a.1 - b.1) > 1e-6 {
                    out.push((a, b));
                }
            };
            match case {
                1 | 14 => push(3, 0),
                2 | 13 => push(0, 1),
                3 | 12 => push(3, 1),
                4 | 11 => push(1, 2),
                6 | 9 => push(0, 2),
                7 | 8 => push(3, 2),
                // Saddles: two opposite corners above the level. Which pair of
                // crossings joins up is genuinely ambiguous from the corners
                // alone; the block mean is the standard tie-breaker and is what
                // keeps contours of different levels from crossing.
                5 => {
                    let mean = v.iter().sum::<f64>() / 4.0;
                    if mean >= level {
                        push(3, 2);
                        push(0, 1);
                    } else {
                        push(3, 0);
                        push(1, 2);
                    }
                }
                10 => {
                    let mean = v.iter().sum::<f64>() / 4.0;
                    if mean >= level {
                        push(3, 0);
                        push(1, 2);
                    } else {
                        push(3, 2);
                        push(0, 1);
                    }
                }
                _ => {}
            }
        }
    }
    out
}

/// Chains unordered segments end to end into polylines.
///
/// Endpoints are matched within [`WELD_TOLERANCE`] rather than exactly. Two
/// blocks sharing an edge compute the same crossing by different arithmetic, so
/// their results differ by an ulp; and where the surface touches the level
/// exactly at a grid point, the epsilon nudge in [`nudge_level`] leaves the
/// contour's single passage through that point as two crossings on adjacent
/// edges a few times 1e-7 of a cell apart. Exact matching leaves the ring in
/// pieces — a cone contoured at a level its grid attains exactly came out as
/// four disconnected quadrant arcs.
///
/// Candidates are gathered from the endpoint's own bucket and all eight
/// neighbours, then distance-checked, so a pair straddling a bucket boundary is
/// still found.
fn chain_segments(segments: Vec<(Pt, Pt)>) -> Vec<Vec<Pt>> {
    // Bucket size equals the weld tolerance, so a true match is always in the
    // 3x3 neighbourhood of buckets.
    let bucket = |p: Pt| -> (i64, i64) {
        (
            (p.0 / WELD_TOLERANCE).floor() as i64,
            (p.1 / WELD_TOLERANCE).floor() as i64,
        )
    };
    let close = |a: Pt, b: Pt| (a.0 - b.0).hypot(a.1 - b.1) <= WELD_TOLERANCE;

    let mut buckets: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    for (i, s) in segments.iter().enumerate() {
        buckets.entry(bucket(s.0)).or_default().push(i);
        buckets.entry(bucket(s.1)).or_default().push(i);
    }

    let mut used = vec![false; segments.len()];
    let mut chains: Vec<Vec<Pt>> = Vec::new();

    for start in 0..segments.len() {
        if used[start] {
            continue;
        }
        used[start] = true;
        let mut chain = vec![segments[start].0, segments[start].1];

        // Extend forward, reverse, extend again: an open contour is followed to
        // both of its ends, and a closed one comes back to where it began.
        for _pass in 0..2 {
            loop {
                let tail = *chain.last().unwrap();
                let (bx, by) = bucket(tail);
                let mut advanced = false;
                'search: for dx in -1..=1 {
                    for dy in -1..=1 {
                        let Some(cands) = buckets.get(&(bx + dx, by + dy)) else {
                            continue;
                        };
                        for &si in cands {
                            if used[si] {
                                continue;
                            }
                            let s = segments[si];
                            let next = if close(s.0, tail) {
                                s.1
                            } else if close(s.1, tail) {
                                s.0
                            } else {
                                continue;
                            };
                            used[si] = true;
                            chain.push(next);
                            advanced = true;
                            break 'search;
                        }
                    }
                }
                if !advanced {
                    break;
                }
            }
            chain.reverse();
        }
        chains.push(chain);
    }
    chains
}

fn polyline_length(coords: &[Coord]) -> f64 {
    coords
        .windows(2)
        .map(|w| ((w[1].x - w[0].x).powi(2) + (w[1].y - w[0].y).powi(2)).sqrt())
        .sum()
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    polygons: bool,
    min_length: Option<f64>,
    close_rings: bool,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let polygons =
        choice_or(args, "geometry_type", &["polyline", "polygon"], "polyline")? == "polygon";
    Ok(Params {
        polygons,
        min_length: opt_positive_f64(args, "min_length")?,
        close_rings: bool_or(args, "close_rings", true)?,
    })
}

/// The requested contour levels, de-duplicated and sorted.
fn parse_levels(args: &ToolArgs) -> Result<Vec<f64>, ToolError> {
    let raw = args
        .get("contour_values")
        .and_then(Value::as_str)
        .ok_or_else(|| {
            ToolError::Validation("missing required parameter 'contour_values'".to_string())
        })?;
    let mut out: Vec<f64> = Vec::new();
    for part in raw.split(',').map(str::trim).filter(|p| !p.is_empty()) {
        let v: f64 = part.parse().map_err(|_| {
            ToolError::Validation(format!("'contour_values' entry '{part}' is not a number"))
        })?;
        if !v.is_finite() {
            return Err(ToolError::Validation(format!(
                "'contour_values' entry '{part}' is not finite"
            )));
        }
        out.push(v);
    }
    if out.is_empty() {
        return Err(ToolError::Validation(
            "'contour_values' must list at least one level".to_string(),
        ));
    }
    out.sort_by(f64::total_cmp);
    out.dedup_by(|a, b| a == b);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vector_common::load_input_layer;
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

    fn raster_of(cols: usize, rows: usize, vals: &[f64]) -> String {
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
                epsg: Some(32610),
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
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, BTreeMap<String, Value>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ContourListTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (layer, out.outputs)
    }

    /// Absolute shoelace area of an unclosed ring.
    fn ring_area(coords: &[Coord]) -> f64 {
        let n = coords.len();
        if n < 3 {
            return 0.0;
        }
        let mut a = 0.0;
        for i in 0..n {
            let p = &coords[i];
            let q = &coords[(i + 1) % n];
            a += p.x * q.y - q.x * p.y;
        }
        (a / 2.0).abs()
    }

    /// A planar ramp: the 25 contour must be a straight line at exactly the
    /// right place, which is the analytic anchor for the interpolation.
    #[test]
    fn contour_of_a_ramp_is_at_the_right_place() {
        // z = 10 * x, sampled at cell centres x = 0.5, 1.5, ... so z = 5, 15, 25...
        let (rows, cols) = (5, 6);
        let mut v = vec![0.0; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                v[r * cols + c] = 10.0 * (c as f64 + 0.5);
            }
        }
        let (layer, outputs) = run(json!({
            "input": raster_of(cols, rows, &v), "contour_values": "25"
        }));
        assert_eq!(outputs["level_count"].as_u64().unwrap(), 1);
        assert_eq!(layer.len(), 1, "a ramp has exactly one contour line");

        let f = layer.iter().next().unwrap();
        let Some(Geometry::LineString(cs)) = f.geometry.as_ref() else {
            panic!("expected a line")
        };
        // z = 25 at x = 2.5 in map units.
        for c in cs {
            assert!(
                (c.x - 2.5).abs() < 1e-6,
                "contour vertex at x = {} should be 2.5",
                c.x
            );
        }
        // It should span the full height of the raster.
        let ys: Vec<f64> = cs.iter().map(|c| c.y).collect();
        let span = ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
            - ys.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(span >= 3.9, "contour should cross the raster, span {span}");
    }

    /// The whole point: an irregular list, each feature labelled with its own
    /// level. `contours_from_raster` cannot express this.
    #[test]
    fn irregular_levels_each_get_their_own_features() {
        let (rows, cols) = (6, 20);
        let mut v = vec![0.0; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                v[r * cols + c] = 5.0 * (c as f64 + 0.5);
            }
        }
        // The IHO-style depth series is not an arithmetic progression.
        let (layer, outputs) = run(json!({
            "input": raster_of(cols, rows, &v), "contour_values": "10, 20, 50, 90"
        }));
        assert_eq!(outputs["level_count"].as_u64().unwrap(), 4);
        assert_eq!(layer.len(), 4, "one line per requested level");

        let ci = layer.schema.field_index("contour").unwrap();
        let mut levels: Vec<f64> = layer
            .iter()
            .map(|f| match f.attributes[ci] {
                FieldValue::Float(v) => v,
                _ => panic!("contour must be a float"),
            })
            .collect();
        levels.sort_by(f64::total_cmp);
        assert_eq!(levels, vec![10.0, 20.0, 50.0, 90.0]);
    }

    /// A cone produces a closed contour ring, flagged as closed.
    #[test]
    fn closed_contours_are_detected() {
        let (rows, cols) = (21, 21);
        let mut v = vec![0.0; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let dx = c as f64 - 10.0;
                let dy = r as f64 - 10.0;
                v[r * cols + c] = 100.0 - (dx * dx + dy * dy).sqrt() * 5.0;
            }
        }
        let (layer, _) = run(json!({
            "input": raster_of(cols, rows, &v), "contour_values": "75"
        }));
        assert_eq!(layer.len(), 1, "a cone gives one ring at each level");
        let f = layer.iter().next().unwrap();
        let closed_idx = layer.schema.field_index("closed").unwrap();
        assert!(
            matches!(f.attributes[closed_idx], FieldValue::Boolean(true)),
            "the ring should be flagged closed"
        );
        // z = 75 is 5 units of radius out from the peak: circumference ~ 31.4.
        let li = layer.schema.field_index("length").unwrap();
        let FieldValue::Float(len) = f.attributes[li] else {
            panic!()
        };
        assert!(
            (len - 31.4).abs() < 3.0,
            "ring length {len} should be near the 31.4 circumference"
        );
    }

    /// Segments are chained, not emitted as thousands of stubs.
    #[test]
    fn segments_are_chained_into_polylines() {
        let (rows, cols) = (10, 10);
        let mut v = vec![0.0; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                v[r * cols + c] = c as f64;
            }
        }
        let (layer, _) = run(json!({
            "input": raster_of(cols, rows, &v), "contour_values": "4.5"
        }));
        assert_eq!(layer.len(), 1, "one chained polyline, not nine stubs");
        let Some(Geometry::LineString(cs)) = layer.iter().next().unwrap().geometry.as_ref() else {
            panic!()
        };
        assert!(cs.len() >= 9, "expected a chain of vertices, got {}", cs.len());
    }

    /// The polygon form returns filled regions at or above each level, nested
    /// as the levels rise.
    #[test]
    fn polygon_form_returns_filled_regions() {
        let (rows, cols) = (11, 11);
        let mut v = vec![0.0; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let dx = c as f64 - 5.0;
                let dy = r as f64 - 5.0;
                v[r * cols + c] = 50.0 - (dx * dx + dy * dy).sqrt() * 5.0;
            }
        }
        let (layer, _) = run(json!({
            "input": raster_of(cols, rows, &v), "contour_values": "20, 40",
            "geometry_type": "polygon"
        }));
        assert_eq!(layer.len(), 2);
        // The higher level encloses less area.
        let ci = layer.schema.field_index("contour").unwrap();
        let mut by_level: Vec<(f64, f64)> = layer
            .iter()
            .map(|f| {
                let FieldValue::Float(level) = f.attributes[ci] else {
                    panic!()
                };
                let Some(Geometry::Polygon { exterior, .. }) = f.geometry.as_ref() else {
                    panic!("expected a polygon")
                };
                (level, ring_area(exterior.coords()))
            })
            .collect();
        by_level.sort_by(|a, b| a.0.total_cmp(&b.0));
        assert!(
            by_level[1].1 < by_level[0].1,
            "the 40 region ({}) should be inside the 20 region ({})",
            by_level[1].1,
            by_level[0].1
        );
    }

    /// A level outside the surface's range yields nothing rather than failing.
    #[test]
    fn level_outside_the_range_is_empty() {
        let (rows, cols) = (4, 4);
        let (layer, outputs) = run(json!({
            "input": raster_of(cols, rows, &[5.0; 16]), "contour_values": "100"
        }));
        assert_eq!(layer.len(), 0);
        assert_eq!(outputs["levels"][0]["feature_count"].as_u64().unwrap(), 0);
    }

    /// No-data blocks are skipped rather than interpolated across.
    #[test]
    fn nodata_blocks_are_skipped() {
        let (rows, cols) = (5, 5);
        let mut v: Vec<f64> = (0..25).map(|i| (i % 5) as f64 * 10.0).collect();
        v[12] = -9999.0; // a hole in the middle
        let (layer, _) = run(json!({
            "input": raster_of(cols, rows, &v), "contour_values": "15"
        }));
        // The contour is interrupted by the hole, so it comes out in pieces
        // rather than crossing invented ground.
        assert!(!layer.is_empty());
        for f in layer.iter() {
            let Some(Geometry::LineString(cs)) = f.geometry.as_ref() else {
                panic!()
            };
            assert!(cs.len() >= 2);
        }
    }

    /// Duplicate levels collapse to one.
    #[test]
    fn duplicate_levels_are_deduplicated() {
        let (rows, cols) = (5, 5);
        let v: Vec<f64> = (0..25).map(|i| (i % 5) as f64 * 10.0).collect();
        let (_, outputs) = run(json!({
            "input": raster_of(cols, rows, &v), "contour_values": "15, 15, 15"
        }));
        assert_eq!(outputs["level_count"].as_u64().unwrap(), 1);
    }

    /// `min_length` drops short specks.
    #[test]
    fn min_length_drops_specks() {
        let (rows, cols) = (12, 12);
        let mut v = vec![0.0; rows * cols];
        // A long ridge plus a single-cell bump far away.
        for r in 0..rows {
            v[r * cols + 5] = 10.0;
        }
        v[10 * cols + 10] = 10.0;
        let src = raster_of(cols, rows, &v);
        let (all, _) = run(json!({ "input": src, "contour_values": "5" }));
        let (long, _) = run(json!({
            "input": raster_of(cols, rows, &v), "contour_values": "5", "min_length": 5.0
        }));
        assert!(
            long.len() < all.len(),
            "min_length should drop the bump: {} vs {}",
            long.len(),
            all.len()
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ContourListTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({"input": "a.tif"})).is_err());
        assert!(bad(json!({"input": "a.tif", "contour_values": ""})).is_err());
        assert!(bad(json!({"input": "a.tif", "contour_values": "10,x"})).is_err());
        assert!(bad(json!({"input": "a.tif", "contour_values": "10", "geometry_type": "point"}))
            .is_err());
        assert!(bad(json!({"input": "a.tif", "contour_values": "10", "min_length": -1})).is_err());
        assert!(bad(json!({"input": "a.tif", "contour_values": "2,5,10,20"})).is_ok());
    }
}


