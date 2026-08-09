//! GeoLibre tool: CFAR ship/target detection in SAR imagery over water.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Detect Bright Ocean Objects* (Image
//! Analyst).
//!
//! ## Why the catalog needs it
//!
//! Metal vessels are corner reflectors: against the near-specular surface of
//! open water they stand 15–25 dB above the background, which is why SAR is the
//! workhorse for maritime domain awareness — it sees at night and through
//! cloud, and it is the only routinely-available way to find vessels that have
//! switched off their AIS transponders.
//!
//! Nothing in either registry can do it. A global threshold fails immediately,
//! because sea clutter brightens with wind and with incidence angle, so any
//! constant that finds small vessels in calm water floods a rough sea with
//! false alarms. `matched_filter_target_detection` is a *spectral* detector for
//! optical imagery and has no notion of a clutter distribution;
//! `detect_image_anomalies` works on a global statistical model rather than a
//! local sliding one; `image_segmentation` groups homogeneous regions, which is
//! the opposite of what a point target is.
//!
//! ## CFAR
//!
//! The standard answer is a **constant false-alarm-rate** detector. For every
//! cell, the sea-clutter mean and standard deviation are estimated from a
//! *ring* of background cells, separated from the cell under test by a guard
//! band wide enough that a vessel's own return cannot contaminate its own
//! background estimate. The cell is declared a target when
//!
//! ```text
//! (x - mu_background) / sigma_background >= threshold
//! ```
//!
//! Because the statistics are local, the detector automatically adapts to wind
//! streaks, swell and the across-swath brightness gradient — hence "constant
//! false-alarm rate".
//!
//! Detections are grouped into 4-connected regions, sized, and emitted as
//! oriented bounding boxes (or their traced outline) with length, width,
//! heading and peak backscatter, so a length filter can discard whitecaps and
//! keep vessels.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;
use wbvector::{Coord, Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::args_common::{band_index, choice_or, f64_or, opt_f64, req_str, usize_or};
use crate::common::{load_input_raster, parse_optional_output};
use crate::sar_common::{
    connected_regions, power_to_db, rasterize_mask, regions_to_geometries, MaskSide, Region,
    SarUnits,
};
use crate::vector_common::{load_input_layer, write_or_store_layer};

pub struct DetectBrightOceanObjectsTool;

impl Tool for DetectBrightOceanObjectsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "detect_bright_ocean_objects",
            display_name: "Detect Bright Ocean Objects",
            summary: "Finds vessels and other bright point targets on water in SAR imagery with a two-parameter CFAR detector, and emits them as oriented bounding boxes with length, width, heading and peak backscatter (ArcGIS Detect Bright Ocean Objects). Neither registry has a local-clutter detector: a global threshold cannot work because sea clutter brightens with wind and incidence angle, matched_filter_target_detection is a spectral optical detector, and detect_image_anomalies uses a global rather than a sliding-window model.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "SAR backscatter raster.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon layer of detected objects. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "geometry_type",
                    description: "'bounding_box' (default) for the oriented minimum rectangle, or 'perimeter' for the traced detection outline.",
                    required: false,
                },
                ToolParamSpec {
                    name: "threshold",
                    description: "CFAR detection threshold in background standard deviations (default 5.0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "guard_size",
                    description: "Half-width in cells of the guard band excluded from the background estimate (default 3). Must be smaller than background_size.",
                    required: false,
                },
                ToolParamSpec {
                    name: "background_size",
                    description: "Half-width in cells of the background ring (default 12).",
                    required: false,
                },
                ToolParamSpec {
                    name: "mask_features",
                    description: "Polygon layer masking land or water. Cells excluded by the mask are neither tested nor used as background.",
                    required: false,
                },
                ToolParamSpec {
                    name: "mask_type",
                    description: "'land_polygon' (default; analyse outside the polygons) or 'water_polygon' (analyse inside them).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_object_length",
                    description: "Discard objects shorter than this, in map units.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_object_length",
                    description: "Discard objects longer than this, in map units.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_object_width",
                    description: "Discard objects narrower than this, in map units.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_object_width",
                    description: "Discard objects wider than this, in map units.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_cells",
                    description: "Discard detections smaller than this many cells (default 2), which suppresses isolated speckle spikes.",
                    required: false,
                },
                ToolParamSpec {
                    name: "input_units",
                    description: "'intensity' (default), 'dn', 'amplitude', or 'db'. Detection runs on linear power, where the clutter statistics are meaningful.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band to search (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input_path = req_str(args, "input")?.to_string();
        let prm = parse_params(args)?;
        let band = band_index(args, "band")?;
        let output = parse_optional_output(args, "output")?;

        let raster = load_input_raster(&input_path)?;
        let (rows, cols) = (raster.rows, raster.cols);

        // Linear power, NaN where invalid.
        let mut power = vec![f64::NAN; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let v = raster.get(band, r as isize, c as isize);
                if v != raster.nodata && v.is_finite() {
                    if let Some(p) = prm.input_units.to_power(v) {
                        power[r * cols + c] = p;
                    }
                }
            }
        }

        // Mask: cells outside the analysis area are neither tested nor allowed
        // to pollute a background estimate — land is far brighter than sea and
        // would raise the local mean enough to hide a vessel next to a coast.
        let analyse = match args.get("mask_features").and_then(Value::as_str) {
            Some(p) if !p.trim().is_empty() => {
                let layer = load_input_layer(p.trim())?;
                rasterize_mask(&raster, &layer, prm.mask_side)
            }
            _ => vec![true; rows * cols],
        };

        ctx.progress.info(&format!(
            "{rows}x{cols}, CFAR threshold {} sigma, guard {} background {}",
            prm.threshold, prm.guard, prm.background
        ));

        let flag = cfar(&power, &analyse, rows, cols, &prm, ctx);
        let regions = connected_regions(&flag, &power, rows, cols, prm.min_cells);
        ctx.progress
            .info(&format!("{} candidate detection(s)", regions.len()));

        let layer = build_layer(&raster, &power, &regions, rows, cols, &prm)?;
        let detected = layer.len();
        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("object_count".to_string(), json!(detected));
        outputs.insert("candidate_count".to_string(), json!(regions.len()));
        outputs.insert(
            "detected_cells".to_string(),
            json!(flag.iter().filter(|f| **f).count()),
        );
        outputs.insert("threshold".to_string(), json!(prm.threshold));
        Ok(ToolRunResult { outputs })
    }
}

/// Two-parameter CFAR over a guard/background ring.
fn cfar(
    power: &[f64],
    analyse: &[bool],
    rows: usize,
    cols: usize,
    prm: &Params,
    ctx: &ToolContext,
) -> Vec<bool> {
    let g = prm.guard as isize;
    let b = prm.background as isize;
    let mut flag = vec![false; rows * cols];

    for r in 0..rows as isize {
        for c in 0..cols as isize {
            let i = r as usize * cols + c as usize;
            if !analyse[i] || !power[i].is_finite() {
                continue;
            }
            // Background ring: inside the background box, outside the guard box.
            let mut sum = 0.0;
            let mut sum_sq = 0.0;
            let mut n = 0usize;
            for dr in -b..=b {
                for dc in -b..=b {
                    if dr.abs() <= g && dc.abs() <= g {
                        continue;
                    }
                    let (rr, cc) = (r + dr, c + dc);
                    if rr < 0 || cc < 0 || rr >= rows as isize || cc >= cols as isize {
                        continue;
                    }
                    let j = rr as usize * cols + cc as usize;
                    if !analyse[j] || !power[j].is_finite() {
                        continue;
                    }
                    sum += power[j];
                    sum_sq += power[j] * power[j];
                    n += 1;
                }
            }
            // Too little background to estimate a distribution from.
            if n < 8 {
                continue;
            }
            let mean = sum / n as f64;
            let var = (sum_sq / n as f64 - mean * mean).max(0.0);
            let sd = var.sqrt();
            if sd <= 0.0 {
                continue;
            }
            if (power[i] - mean) / sd >= prm.threshold {
                flag[i] = true;
            }
        }
        ctx.progress.progress((r as f64 + 1.0) / rows as f64);
    }
    flag
}

/// Oriented extent of a set of cells: length, width, heading, and the four
/// corners of the minimum-area-ish rectangle aligned to the principal axis.
struct Oriented {
    length: f64,
    width: f64,
    /// Compass heading of the long axis, degrees from north, in [0, 180).
    heading: f64,
    corners: [(f64, f64); 4],
}

/// Fits an oriented box to a region by principal-axis analysis of its cell
/// centres, in world coordinates.
///
/// A single-cell region has no principal axis, so it falls back to the cell
/// footprint — returning a zero-size box would make every length filter reject
/// exactly the smallest targets the detector is most likely to find.
fn oriented_extent(raster: &Raster, cells: &[usize], cols: usize) -> Oriented {
    let (csx, csy) = (raster.cell_size_x, raster.cell_size_y);
    let y_max = raster.y_min + raster.rows as f64 * csy;
    let pts: Vec<(f64, f64)> = cells
        .iter()
        .map(|&i| {
            let (r, c) = (i / cols, i % cols);
            (
                raster.x_min + (c as f64 + 0.5) * csx,
                y_max - (r as f64 + 0.5) * csy,
            )
        })
        .collect();

    let n = pts.len() as f64;
    let cx = pts.iter().map(|p| p.0).sum::<f64>() / n;
    let cy = pts.iter().map(|p| p.1).sum::<f64>() / n;

    // Covariance of the centred points.
    let (mut sxx, mut sxy, mut syy) = (0.0, 0.0, 0.0);
    for p in &pts {
        let (dx, dy) = (p.0 - cx, p.1 - cy);
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    sxx /= n;
    sxy /= n;
    syy /= n;

    // Principal axis = dominant eigenvector of the 2x2 covariance.
    let theta = if sxy.abs() < 1e-15 && (sxx - syy).abs() < 1e-15 {
        0.0
    } else {
        0.5 * (2.0 * sxy).atan2(sxx - syy)
    };
    let (ct, st) = (theta.cos(), theta.sin());

    // Extents along and across the principal axis, padded by the half-cell
    // footprint so a one-cell detection still has the size of a cell.
    let (mut amin, mut amax, mut bmin, mut bmax) =
        (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
    for p in &pts {
        let (dx, dy) = (p.0 - cx, p.1 - cy);
        let a = dx * ct + dy * st;
        let b = -dx * st + dy * ct;
        amin = amin.min(a);
        amax = amax.max(a);
        bmin = bmin.min(b);
        bmax = bmax.max(b);
    }
    let pad_a = 0.5 * (csx * ct.abs() + csy * st.abs());
    let pad_b = 0.5 * (csx * st.abs() + csy * ct.abs());
    amin -= pad_a;
    amax += pad_a;
    bmin -= pad_b;
    bmax += pad_b;

    let mut length = amax - amin;
    let mut width = bmax - bmin;
    let mut axis = theta;
    if width > length {
        std::mem::swap(&mut length, &mut width);
        axis += PI / 2.0;
    }

    let corner = |a: f64, b: f64| (cx + a * ct - b * st, cy + a * st + b * ct);
    let corners = [
        corner(amin, bmin),
        corner(amax, bmin),
        corner(amax, bmax),
        corner(amin, bmax),
    ];

    // Compass heading: measured from north, clockwise. `axis` is a mathematical
    // angle from east, counter-clockwise.
    let mut heading = (90.0 - axis.to_degrees()).rem_euclid(180.0);
    if !heading.is_finite() {
        heading = 0.0;
    }

    Oriented {
        length,
        width,
        heading,
        corners,
    }
}

/// Builds the output layer, applying the size filters.
fn build_layer(
    raster: &Raster,
    power: &[f64],
    regions: &[Region],
    rows: usize,
    cols: usize,
    prm: &Params,
) -> Result<Layer, ToolError> {
    let mut layer = Layer::new("bright_ocean_objects");
    layer.geom_type = Some(GeometryType::Polygon);
    if let Some(e) = raster.crs.epsg {
        layer = layer.with_crs_epsg(e);
    }
    layer.add_field(FieldDef::new("id", FieldType::Integer));
    layer.add_field(FieldDef::new("cell_count", FieldType::Integer));
    layer.add_field(FieldDef::new("area", FieldType::Float));
    layer.add_field(FieldDef::new("length", FieldType::Float));
    layer.add_field(FieldDef::new("width", FieldType::Float));
    layer.add_field(FieldDef::new("heading", FieldType::Float));
    layer.add_field(FieldDef::new("mean_power", FieldType::Float));
    layer.add_field(FieldDef::new("max_power", FieldType::Float));
    layer.add_field(FieldDef::new("max_db", FieldType::Float));

    // Traced outlines, only needed for the perimeter form.
    let traced = if prm.perimeter {
        regions_to_geometries(raster, regions, rows, cols)?
    } else {
        Vec::new()
    };

    let cell_area = raster.cell_size_x * raster.cell_size_y;
    let mut fid = 0u64;
    for (i, reg) in regions.iter().enumerate() {
        let ext = oriented_extent(raster, &reg.cells, cols);
        if let Some(v) = prm.min_length {
            if ext.length < v {
                continue;
            }
        }
        if let Some(v) = prm.max_length {
            if ext.length > v {
                continue;
            }
        }
        if let Some(v) = prm.min_width {
            if ext.width < v {
                continue;
            }
        }
        if let Some(v) = prm.max_width {
            if ext.width > v {
                continue;
            }
        }

        let geom = if prm.perimeter {
            match traced.iter().find(|(idx, _)| *idx == i) {
                Some((_, g)) => g.clone(),
                // A region whose ring trace collapsed has no outline to emit;
                // dropping it is better than substituting a different shape.
                None => continue,
            }
        } else {
            let coords: Vec<Coord> = ext
                .corners
                .iter()
                .map(|&(x, y)| Coord::xy(x, y))
                .collect();
            Geometry::Polygon {
                exterior: Ring::new(coords),
                interiors: Vec::new(),
            }
        };

        let max_power = reg
            .cells
            .iter()
            .map(|&c| power[c])
            .fold(f64::NEG_INFINITY, f64::max);

        let mut f = Feature::with_geometry(fid, geom, layer.schema.len());
        f.set_by_index(0, FieldValue::Integer(fid as i64));
        f.set_by_index(1, FieldValue::Integer(reg.cells.len() as i64));
        f.set_by_index(2, FieldValue::Float(reg.cells.len() as f64 * cell_area));
        f.set_by_index(3, FieldValue::Float(ext.length));
        f.set_by_index(4, FieldValue::Float(ext.width));
        f.set_by_index(5, FieldValue::Float(ext.heading));
        f.set_by_index(6, FieldValue::Float(reg.mean()));
        f.set_by_index(7, FieldValue::Float(max_power));
        f.set_by_index(
            8,
            match power_to_db(max_power) {
                Some(v) => FieldValue::Float(v),
                None => FieldValue::Null,
            },
        );
        layer.push(f);
        fid += 1;
    }
    Ok(layer)
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    threshold: f64,
    guard: usize,
    background: usize,
    min_cells: usize,
    perimeter: bool,
    mask_side: MaskSide,
    input_units: SarUnits,
    min_length: Option<f64>,
    max_length: Option<f64>,
    min_width: Option<f64>,
    max_width: Option<f64>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let threshold = f64_or(args, "threshold", 5.0)?;
    if !threshold.is_finite() || threshold <= 0.0 {
        return Err(ToolError::Validation(format!(
            "'threshold' must be positive, got {threshold}"
        )));
    }
    let guard = usize_or(args, "guard_size", 3)?;
    let background = usize_or(args, "background_size", 12)?;
    if background <= guard {
        return Err(ToolError::Validation(format!(
            "'background_size' ({background}) must exceed 'guard_size' ({guard}), \
             otherwise the background ring is empty"
        )));
    }
    let min_cells = usize_or(args, "min_cells", 2)?.max(1);
    let perimeter = choice_or(
        args,
        "geometry_type",
        &["bounding_box", "perimeter"],
        "bounding_box",
    )? == "perimeter";
    let mask_side = MaskSide::parse(
        args.get("mask_type")
            .and_then(Value::as_str)
            .unwrap_or(""),
    )?;
    let input_units = SarUnits::parse(args.get("input_units").and_then(Value::as_str).unwrap_or(""))?;

    let bound = |key: &str| -> Result<Option<f64>, ToolError> {
        match opt_f64(args, key)? {
            None => Ok(None),
            Some(v) if v >= 0.0 && v.is_finite() => Ok(Some(v)),
            Some(v) => Err(ToolError::Validation(format!(
                "'{key}' must be a non-negative distance, got {v}"
            ))),
        }
    };
    let min_length = bound("min_object_length")?;
    let max_length = bound("max_object_length")?;
    let min_width = bound("min_object_width")?;
    let max_width = bound("max_object_width")?;
    if let (Some(lo), Some(hi)) = (min_length, max_length) {
        if lo > hi {
            return Err(ToolError::Validation(
                "'min_object_length' exceeds 'max_object_length'".to_string(),
            ));
        }
    }
    if let (Some(lo), Some(hi)) = (min_width, max_width) {
        if lo > hi {
            return Err(ToolError::Validation(
                "'min_object_width' exceeds 'max_object_width'".to_string(),
            ));
        }
    }

    Ok(Params {
        threshold,
        guard,
        background,
        min_cells,
        perimeter,
        mask_side,
        input_units,
        min_length,
        max_length,
        min_width,
        max_width,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Deterministic clutter in [0, 1) — no RNG, so WASM matches native.
    fn clutter(i: usize) -> f64 {
        let mut x = (i as u64).wrapping_mul(6364136223846793005).wrapping_add(1);
        x ^= x >> 33;
        x = x.wrapping_mul(0xff51afd7ed558ccd);
        x ^= x >> 33;
        (x % 10007) as f64 / 10007.0
    }

    fn raster_of(cols: usize, rows: usize, vals: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 10.0,
            cell_size_y: Some(10.0),
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
        let out = DetectBrightOceanObjectsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (layer, out.outputs)
    }

    /// A calm sea with two bright targets: both are found and nothing else is.
    #[test]
    fn finds_bright_targets_in_clutter() {
        let (rows, cols) = (60, 60);
        let mut v: Vec<f64> = (0..rows * cols).map(|i| 0.05 + 0.02 * clutter(i)).collect();
        // Two 2x2 "vessels", far apart and away from the border.
        for &(r0, c0) in &[(20usize, 20usize), (40, 42)] {
            for dr in 0..2 {
                for dc in 0..2 {
                    v[(r0 + dr) * cols + c0 + dc] = 3.0;
                }
            }
        }
        let (layer, outputs) = run(json!({
            "input": raster_of(cols, rows, &v),
            "threshold": 6.0, "guard_size": 3, "background_size": 10
        }));
        assert_eq!(
            layer.len(),
            2,
            "expected exactly the two vessels, got {} ({} candidates)",
            layer.len(),
            outputs["candidate_count"]
        );
        // Each is 2x2 cells of 10 m -> about 20 m across.
        let li = layer.schema.field_index("length").unwrap();
        for f in layer.iter() {
            let FieldValue::Float(len) = f.attributes[li] else {
                panic!("length must be a float")
            };
            assert!(
                (15.0..35.0).contains(&len),
                "implausible vessel length {len}"
            );
        }
    }

    /// This is the property a global threshold cannot have: the same target
    /// must be found on both a calm and a rough sea, because the detector
    /// rescales to the local clutter.
    #[test]
    fn adapts_to_the_local_clutter_level() {
        let (rows, cols) = (48, 48);
        let build = |sea: f64| {
            let mut v: Vec<f64> = (0..rows * cols)
                .map(|i| sea * (0.8 + 0.4 * clutter(i)))
                .collect();
            for dr in 0..2 {
                for dc in 0..2 {
                    // Target sits 30x above its own local sea state.
                    v[(24 + dr) * cols + 24 + dc] = sea * 30.0;
                }
            }
            v
        };
        for sea in [0.02, 0.5] {
            let (layer, _) = run(json!({
                "input": raster_of(cols, rows, &build(sea)),
                "threshold": 6.0, "guard_size": 3, "background_size": 10
            }));
            assert_eq!(
                layer.len(),
                1,
                "sea state {sea}: expected 1 detection, got {}",
                layer.len()
            );
        }
    }

    /// The size filters actually discard objects.
    #[test]
    fn length_filter_rejects_small_objects() {
        let (rows, cols) = (48, 48);
        let mut v: Vec<f64> = (0..rows * cols).map(|i| 0.05 + 0.02 * clutter(i)).collect();
        // A 1x4 "long" object and a 2x2 "small" one.
        for dc in 0..4 {
            v[16 * cols + 10 + dc] = 3.0;
        }
        for dr in 0..2 {
            for dc in 0..2 {
                v[(32 + dr) * cols + 30 + dc] = 3.0;
            }
        }
        let src = raster_of(cols, rows, &v);
        let (all, _) = run(json!({
            "input": src.clone(), "threshold": 6.0, "guard_size": 3, "background_size": 10
        }));
        assert_eq!(all.len(), 2, "both objects should be found unfiltered");

        // The 4-cell object spans about 50 m; the 2x2 about 20 m.
        let (long_only, _) = run(json!({
            "input": src, "threshold": 6.0, "guard_size": 3, "background_size": 10,
            "min_object_length": 35.0
        }));
        assert_eq!(
            long_only.len(),
            1,
            "length filter should keep only the long object"
        );
    }

    /// Land must be excluded, both from testing and from the background — a
    /// bright coastline next to open water would otherwise both trigger
    /// detections and mask real ones nearby.
    #[test]
    fn mask_excludes_land() {
        let (rows, cols) = (40, 40);
        // Left half is bright "land", right half calm sea with one vessel.
        let mut v = vec![0.0; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                v[r * cols + c] = if c < 20 {
                    2.0 + clutter(r * cols + c)
                } else {
                    0.05 + 0.02 * clutter(r * cols + c)
                };
            }
        }
        for dr in 0..2 {
            for dc in 0..2 {
                v[(20 + dr) * cols + 30 + dc] = 3.0;
            }
        }
        let src = raster_of(cols, rows, &v);

        // Land polygon covers the left half (x < 200 in map units).
        let mut land = Layer::new("land");
        land.geom_type = Some(GeometryType::Polygon);
        land = land.with_crs_epsg(32610);
        land.add_field(FieldDef::new("id", FieldType::Integer));
        let ring = Ring::new(vec![
            Coord::xy(0.0, 0.0),
            Coord::xy(200.0, 0.0),
            Coord::xy(200.0, 400.0),
            Coord::xy(0.0, 400.0),
        ]);
        let mut f = Feature::with_geometry(
            0,
            Geometry::Polygon {
                exterior: ring,
                interiors: Vec::new(),
            },
            land.schema.len(),
        );
        f.set_by_index(0, FieldValue::Integer(1));
        land.push(f);
        let land_path = write_or_store_layer(land, None).unwrap();

        let (masked, _) = run(json!({
            "input": src, "threshold": 6.0, "guard_size": 3, "background_size": 8,
            "mask_features": land_path, "mask_type": "land_polygon"
        }));
        assert_eq!(
            masked.len(),
            1,
            "only the vessel on open water should survive masking, got {}",
            masked.len()
        );
    }

    /// The perimeter form emits the traced detection outline instead of a box.
    #[test]
    fn perimeter_geometry_is_available() {
        let (rows, cols) = (40, 40);
        let mut v: Vec<f64> = (0..rows * cols).map(|i| 0.05 + 0.02 * clutter(i)).collect();
        for dr in 0..3 {
            for dc in 0..3 {
                v[(20 + dr) * cols + 20 + dc] = 4.0;
            }
        }
        let (layer, _) = run(json!({
            "input": raster_of(cols, rows, &v), "threshold": 6.0,
            "guard_size": 3, "background_size": 9, "geometry_type": "perimeter"
        }));
        assert_eq!(layer.len(), 1);
        let Some(Geometry::Polygon { exterior, .. }) = layer.iter().next().unwrap().geometry.as_ref()
        else {
            panic!("expected a polygon");
        };
        // A traced 3x3 block is a square ring of 4 corners, not the 4-corner
        // oriented box — both have 4 points, so check the extent instead.
        let xs: Vec<f64> = exterior.coords().iter().map(|c| c.x).collect();
        let span = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
            - xs.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(
            (span - 30.0).abs() < 1e-6,
            "traced outline should span exactly 3 cells (30 m), got {span}"
        );
    }

    /// An empty sea produces no detections rather than an error.
    #[test]
    fn calm_sea_yields_nothing() {
        let (rows, cols) = (40, 40);
        let v: Vec<f64> = (0..rows * cols).map(|i| 0.05 + 0.02 * clutter(i)).collect();
        let (layer, outputs) = run(json!({
            "input": raster_of(cols, rows, &v),
            "threshold": 8.0, "guard_size": 3, "background_size": 10
        }));
        assert_eq!(layer.len(), 0);
        assert_eq!(outputs["object_count"].as_u64().unwrap(), 0);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            DetectBrightOceanObjectsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({"input": "a.tif", "threshold": 0.0})).is_err());
        // A background ring no larger than the guard band is empty.
        assert!(bad(json!({"input": "a.tif", "guard_size": 5, "background_size": 5})).is_err());
        assert!(bad(json!({"input": "a.tif", "mask_type": "sky"})).is_err());
        assert!(bad(json!({"input": "a.tif", "geometry_type": "circle"})).is_err());
        assert!(
            bad(json!({"input": "a.tif", "min_object_length": 50, "max_object_length": 10}))
                .is_err()
        );
        assert!(bad(json!({"input": "a.tif"})).is_ok());
    }
}
