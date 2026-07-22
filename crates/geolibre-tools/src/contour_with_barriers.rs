//! GeoLibre tool: contour a surface with barriers the lines must not cross.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Contour With Barriers* (Spatial
//! Analyst): trace iso-value lines of a raster while treating a barrier vector
//! layer (faults, cliffs, coastlines) as hard discontinuities that contour
//! lines terminate at instead of crossing. The bundled `contours_from_raster`
//! ignores barriers, so contours bleed across features that should break them;
//! this complements the already-shipped `interpolate_with_barriers`.
//!
//! ## Method
//!
//! Standard marching squares over the cell-centre grid: for each contour level
//! and each 2×2 block of neighbouring cells, the level's crossings on the four
//! edges are linearly interpolated and connected per the 16-case table. Barrier
//! geometry is rasterized onto the grid — every cell a barrier passes through is
//! flagged — and any marching-squares block touching a flagged cell (or a
//! no-data cell) is skipped, so contours stop cleanly at the barrier rather than
//! jumping across it. The resulting segments are chained into polylines by
//! snapping shared endpoints.
//!
//! Levels come from an `interval` (+ optional `base`) or an explicit `levels`
//! list. Output is a line layer with a `level` attribute; barriers as
//! LineString, MultiLineString, Polygon, or MultiPolygon boundaries are all
//! honoured. v1 emits contour lines (not filled polygons).

use std::collections::{BTreeMap, HashMap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::common::{band_to_vec, load_input_raster};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct ContourWithBarriersTool;

impl Tool for ContourWithBarriersTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "contour_with_barriers",
            display_name: "Contour With Barriers",
            summary: "Trace iso-value contour lines of a surface that terminate at barrier features instead of crossing them.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input continuous surface raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output line vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "barriers",
                    description: "Optional barrier vector layer (lines or polygons). Contours will not cross these features.",
                    required: false,
                },
                ToolParamSpec {
                    name: "interval",
                    description: "Contour interval in surface units. Levels are base + k*interval within the data range. Ignored if 'levels' is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "base",
                    description: "Base value the interval is measured from (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "levels",
                    description: "Explicit comma-separated contour values (e.g. \"10,20,30\"). Overrides interval/base.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to read (default 1).",
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
        let band = parse_optional_u64(args, "band")?.unwrap_or(1).max(1) as isize - 1;

        let raster = load_input_raster(input)?;
        if band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "band out of range (raster has {} band(s))",
                raster.bands
            )));
        }
        let rows = raster.rows;
        let cols = raster.cols;
        let nodata = raster.nodata;
        let data = band_to_vec(&raster, band);
        let csx = raster.cell_size_x;
        let csy = raster.cell_size_y;
        let x_min = raster.x_min;
        let y_max = raster.y_min + rows as f64 * csy;

        // Data range for interval-based levels.
        let (mut vmin, mut vmax) = (f64::INFINITY, f64::NEG_INFINITY);
        for &v in &data {
            if v != nodata && v.is_finite() {
                vmin = vmin.min(v);
                vmax = vmax.max(v);
            }
        }
        if !vmin.is_finite() {
            return Err(ToolError::Execution(
                "raster band has no valid cells".to_string(),
            ));
        }
        let levels = resolve_levels(args, vmin, vmax)?;

        // Rasterize barriers to a per-cell blocked mask.
        let mut blocked = vec![false; rows * cols];
        let mut barrier_cells = 0usize;
        if let Some(bpath) = parse_optional_str(args, "barriers")? {
            let barriers = load_input_layer(bpath)?;
            for f in &barriers.features {
                if let Some(g) = &f.geometry {
                    rasterize_barrier(g, x_min, y_max, csx, csy, rows, cols, &mut blocked);
                }
            }
            barrier_cells = blocked.iter().filter(|&&b| b).count();
        }

        ctx.progress.info(&format!(
            "contouring {} level(s) over {rows}x{cols} cells, {barrier_cells} barrier cell(s)",
            levels.len()
        ));

        let mut out_layer = Layer::new("contours");
        out_layer.geom_type = Some(GeometryType::LineString);
        if let Some(e) = raster.crs.epsg {
            out_layer = out_layer.with_crs_epsg(e);
        }
        out_layer.add_field(FieldDef::new("level", FieldType::Float));

        // Cell-centre world coordinate.
        let px = |c: usize| x_min + (c as f64 + 0.5) * csx;
        let py = |r: usize| y_max - (r as f64 + 0.5) * csy;
        let val = |r: usize, c: usize| data[r * cols + c];
        let ok = |r: usize, c: usize| {
            let v = val(r, c);
            v != nodata && v.is_finite() && !blocked[r * cols + c]
        };

        let mut fid = 0u64;
        let mut total_segments = 0usize;
        for &z in &levels {
            let mut segs: Vec<(Coord, Coord)> = Vec::new();
            for r in 0..rows.saturating_sub(1) {
                for c in 0..cols.saturating_sub(1) {
                    // Skip a block touching nodata or a barrier cell.
                    if !(ok(r, c) && ok(r, c + 1) && ok(r + 1, c) && ok(r + 1, c + 1)) {
                        continue;
                    }
                    let tl = val(r, c);
                    let tr = val(r, c + 1);
                    let bl = val(r + 1, c);
                    let br = val(r + 1, c + 1);
                    let (xl, xr) = (px(c), px(c + 1));
                    let (yt, yb) = (py(r), py(r + 1));
                    march_cell(z, tl, tr, br, bl, xl, xr, yt, yb, &mut segs);
                }
            }
            total_segments += segs.len();
            // Chain segments into polylines and emit.
            for line in chain_segments(segs, csx.min(csy)) {
                let mut f = Feature::with_geometry(fid, Geometry::LineString(line), 1);
                f.set_by_index(0, FieldValue::Float(z));
                out_layer.push(f);
                fid += 1;
            }
        }

        ctx.progress.info(&format!(
            "traced {} contour line(s) from {total_segments} segment(s)",
            out_layer.len()
        ));

        let feature_count = out_layer.len();
        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("num_levels".to_string(), json!(levels.len()));
        outputs.insert("barrier_cells".to_string(), json!(barrier_cells));
        outputs.insert("segments".to_string(), json!(total_segments));
        Ok(ToolRunResult { outputs })
    }
}

// ── Marching squares ────────────────────────────────────────────────────────

/// Emits the contour segment(s) for one 2×2 block at level `z`. Corners are
/// top-left, top-right, bottom-right, bottom-left with the given world extents.
#[allow(clippy::too_many_arguments)]
fn march_cell(
    z: f64,
    tl: f64,
    tr: f64,
    br: f64,
    bl: f64,
    xl: f64,
    xr: f64,
    yt: f64,
    yb: f64,
    out: &mut Vec<(Coord, Coord)>,
) {
    // Bit per corner that is at or above the level.
    let case = (tl >= z) as u8
        | (((tr >= z) as u8) << 1)
        | (((br >= z) as u8) << 2)
        | (((bl >= z) as u8) << 3);
    if case == 0 || case == 15 {
        return;
    }
    // Edge crossing points (interpolated). T: top (tl-tr), R: right (tr-br),
    // B: bottom (br-bl), L: left (bl-tl).
    let t = || Coord::xy(lerp(xl, xr, frac(tl, tr, z)), yt);
    let r = || Coord::xy(xr, lerp(yt, yb, frac(tr, br, z)));
    let b = || Coord::xy(lerp(xl, xr, frac(bl, br, z)), yb);
    let l = || Coord::xy(xl, lerp(yt, yb, frac(tl, bl, z)));
    match case {
        1 | 14 => out.push((l(), t())),
        2 | 13 => out.push((t(), r())),
        3 | 12 => out.push((l(), r())),
        4 | 11 => out.push((r(), b())),
        6 | 9 => out.push((t(), b())),
        7 | 8 => out.push((l(), b())),
        5 => {
            // saddle: connect L-T and R-B
            out.push((l(), t()));
            out.push((r(), b()));
        }
        10 => {
            out.push((t(), r()));
            out.push((b(), l()));
        }
        _ => {}
    }
}

/// Fraction along an edge from `a` to `b` at which the value crosses `z`.
fn frac(a: f64, b: f64, z: f64) -> f64 {
    if (b - a).abs() < f64::EPSILON {
        0.5
    } else {
        ((z - a) / (b - a)).clamp(0.0, 1.0)
    }
}

fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

// ── Segment chaining ─────────────────────────────────────────────────────────

/// Chains 2-point segments into polylines by matching endpoints snapped to a
/// grid of cell `snap`. Open chains and closed loops are both produced.
fn chain_segments(segs: Vec<(Coord, Coord)>, snap: f64) -> Vec<Vec<Coord>> {
    let tol = (snap * 1e-3).max(f64::MIN_POSITIVE);
    let key = |c: &Coord| ((c.x / tol).round() as i64, (c.y / tol).round() as i64);

    // Adjacency: node key -> list of (segment index, which end).
    let mut adj: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    for (i, (a, b)) in segs.iter().enumerate() {
        adj.entry(key(a)).or_default().push(i);
        adj.entry(key(b)).or_default().push(i);
    }
    let mut used = vec![false; segs.len()];
    let mut lines = Vec::new();

    for start in 0..segs.len() {
        if used[start] {
            continue;
        }
        used[start] = true;
        let mut chain = vec![segs[start].0.clone(), segs[start].1.clone()];
        // Extend forward from the tail.
        loop {
            let tail = chain.last().unwrap().clone();
            match next_seg(&adj, &segs, &used, key(&tail), tol) {
                Some((ni, other)) => {
                    used[ni] = true;
                    chain.push(other);
                }
                None => break,
            }
        }
        // Extend backward from the head.
        loop {
            let head = chain.first().unwrap().clone();
            match next_seg(&adj, &segs, &used, key(&head), tol) {
                Some((ni, other)) => {
                    used[ni] = true;
                    chain.insert(0, other);
                }
                None => break,
            }
        }
        if chain.len() >= 2 {
            lines.push(chain);
        }
    }
    lines
}

/// Finds an unused segment incident to node `k`, returning its index and the
/// far endpoint.
fn next_seg(
    adj: &HashMap<(i64, i64), Vec<usize>>,
    segs: &[(Coord, Coord)],
    used: &[bool],
    k: (i64, i64),
    tol: f64,
) -> Option<(usize, Coord)> {
    let cands = adj.get(&k)?;
    for &i in cands {
        if used[i] {
            continue;
        }
        let (a, b) = &segs[i];
        let ka = ((a.x / tol).round() as i64, (a.y / tol).round() as i64);
        if ka == k {
            return Some((i, b.clone()));
        } else {
            return Some((i, a.clone()));
        }
    }
    None
}

// ── Barrier rasterization ────────────────────────────────────────────────────

/// Marks every grid cell a barrier geometry passes through as blocked, by
/// walking each segment and sampling it at ~half-cell spacing.
#[allow(clippy::too_many_arguments)]
fn rasterize_barrier(
    geom: &Geometry,
    x_min: f64,
    y_max: f64,
    csx: f64,
    csy: f64,
    rows: usize,
    cols: usize,
    blocked: &mut [bool],
) {
    let mut mark = |x: f64, y: f64| {
        let c = ((x - x_min) / csx - 0.5).round();
        let r = ((y_max - y) / csy - 0.5).round();
        if r >= 0.0 && c >= 0.0 && (r as usize) < rows && (c as usize) < cols {
            blocked[r as usize * cols + c as usize] = true;
        }
    };
    let step = 0.5 * csx.min(csy);
    let mut walk = |pts: &[Coord]| {
        for w in pts.windows(2) {
            let (a, b) = (&w[0], &w[1]);
            let d = ((b.x - a.x).powi(2) + (b.y - a.y).powi(2)).sqrt();
            let n = (d / step).ceil().max(1.0) as usize;
            for k in 0..=n {
                let t = k as f64 / n as f64;
                mark(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t);
            }
        }
    };
    match geom {
        Geometry::LineString(cs) => walk(cs),
        Geometry::MultiLineString(parts) => parts.iter().for_each(|p| walk(p)),
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            walk(exterior.coords());
            interiors.iter().for_each(|r| walk(r.coords()));
        }
        Geometry::MultiPolygon(parts) => {
            for (ext, ints) in parts {
                walk(ext.coords());
                ints.iter().for_each(|r| walk(r.coords()));
            }
        }
        _ => {}
    }
}

// ── Levels ───────────────────────────────────────────────────────────────────

fn resolve_levels(args: &ToolArgs, vmin: f64, vmax: f64) -> Result<Vec<f64>, ToolError> {
    if let Some(s) = parse_optional_str(args, "levels")? {
        let mut v = Vec::new();
        for tok in s.split(',').map(str::trim).filter(|t| !t.is_empty()) {
            v.push(
                tok.parse::<f64>()
                    .map_err(|_| ToolError::Validation(format!("level '{tok}' is not a number")))?,
            );
        }
        if v.is_empty() {
            return Err(ToolError::Validation("'levels' has no values".to_string()));
        }
        v.sort_by(|a, b| a.total_cmp(b));
        v.dedup();
        return Ok(v);
    }
    let interval = parse_optional_f64(args, "interval")?.unwrap_or_else(|| {
        // Default: ~10 contours across the range.
        let span = vmax - vmin;
        if span > 0.0 {
            span / 10.0
        } else {
            1.0
        }
    });
    if !(interval > 0.0 && interval.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'interval' must be a positive number".to_string(),
        ));
    }
    let base = parse_optional_f64(args, "base")?.unwrap_or(0.0);
    // First level >= vmin.
    let mut k = ((vmin - base) / interval).ceil();
    let mut levels = Vec::new();
    loop {
        let z = base + k * interval;
        if z > vmax {
            break;
        }
        if z >= vmin {
            levels.push(z);
        }
        k += 1.0;
        if levels.len() > 100_000 {
            break; // runaway guard
        }
    }
    if levels.is_empty() {
        return Err(ToolError::Execution(
            "no contour levels fall within the data range".to_string(),
        ));
    }
    Ok(levels)
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

fn parse_optional_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64()),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
            .map(Some)
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be an integer"))),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be an integer"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{memory_store as rmem, DataType, Raster, RasterConfig};
    use wbvector::{memory_store as vmem, Layer as VLayer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// A pure west-east gradient: contours are vertical lines.
    fn gradient_raster(rows: usize, cols: usize) -> String {
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
            crs: Default::default(),
            metadata: vec![],
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, col as f64).unwrap();
            }
        }
        let id = rmem::put_raster(r);
        rmem::make_raster_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, VLayer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ContourWithBarriersTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn total_length(layer: &VLayer) -> f64 {
        layer
            .features
            .iter()
            .filter_map(|f| f.geometry.as_ref())
            .map(|g| match g {
                Geometry::LineString(cs) => cs
                    .windows(2)
                    .map(|w| ((w[1].x - w[0].x).powi(2) + (w[1].y - w[0].y).powi(2)).sqrt())
                    .sum(),
                _ => 0.0,
            })
            .sum()
    }

    /// A west-east gradient produces the requested number of contour levels.
    #[test]
    fn contours_a_gradient() {
        let input = gradient_raster(20, 20);
        let (out, layer) = run(json!({ "input": input, "levels": "5,10,15" }));
        assert_eq!(out.outputs["num_levels"], json!(3));
        // Each level should appear (3 vertical contour lines, one per level).
        assert!(layer.len() >= 3);
        // level attribute present
        assert!(layer.features.iter().all(|f| f
            .get(&layer.schema, "level")
            .map(|v| v.as_f64().is_some())
            .unwrap_or(false)));
    }

    /// A vertical barrier at x=10 shortens (breaks) the horizontal contour that
    /// would otherwise cross it.
    #[test]
    fn barrier_breaks_contours() {
        let input = gradient_raster(20, 20);
        let (no_barrier, nb_layer) = run(json!({ "input": input, "levels": "10" }));

        // Build a horizontal barrier spanning the whole width at mid-height,
        // which every vertical contour must cross.
        let mut barrier = VLayer::new("wall");
        barrier
            .add_feature(
                Some(Geometry::LineString(vec![
                    Coord::xy(-1.0, 10.0),
                    Coord::xy(25.0, 10.0),
                ])),
                &[],
            )
            .unwrap();
        let bid = vmem::put_vector(barrier);
        let bpath = vmem::make_vector_memory_path(&bid);

        let (with, wb_layer) = run(json!({ "input": input, "levels": "10", "barriers": bpath }));
        assert!(with.outputs["barrier_cells"].as_u64().unwrap() > 0);
        // The barrier removes contour length near y=10.
        assert!(
            total_length(&wb_layer) < total_length(&nb_layer),
            "barrier should shorten total contour length ({} vs {})",
            total_length(&wb_layer),
            total_length(&nb_layer)
        );
        let _ = no_barrier;
    }

    #[test]
    fn rejects_missing_input() {
        let tool = ContourWithBarriersTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "x.tif" })).is_ok());
    }
}
