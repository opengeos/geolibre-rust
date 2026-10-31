//! GeoLibre tool: which observers can see each cell, and how high a target
//! would have to be raised to be seen.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Observer Points* and the `OBSERVERS`
//! analysis type of *Viewshed 2* (3D Analyst / Spatial Analyst).
//!
//! ## Why the catalog needs it
//!
//! The bundled `viewshed` "counts how many stations can see each cell" — it
//! returns a scalar frequency and nothing else. That answers "how exposed is
//! this spot?" but not the two questions siting analysis actually asks:
//!
//! * **Which** observer covers this parcel? A count of 1 does not say whether
//!   it is the tower you are about to decommission.
//! * **How high** would an antenna have to be here to be seen at all? A count
//!   of 0 gives no way to tell "two metres short" from "behind a mountain".
//!
//! Neither is recoverable from a frequency raster, and re-running `viewshed`
//! once per observer to recover the first still cannot produce the second.
//!
//! ## Outputs
//!
//! * `output` — a per-observer **bit mask** (observer *k* sets bit *k*), or a
//!   plain frequency count under `analysis_type: frequency`. The mask is
//!   written as 64-bit so all 32 observers stay exactly representable; a 32-bit
//!   float would start rounding at observer 25.
//! * `output_agl` — the **above-ground level** raster: the smallest height a
//!   target at that cell would have to be raised to become visible to at least
//!   one observer. Zero where the ground is already visible.
//! * `output_table` — one row per observer with the cells and area it sees.
//!
//! Visibility uses the standard line-of-sight test: walk the ray from observer
//! to target, tracking the largest elevation angle met so far; the target is
//! visible when its own angle exceeds that horizon. Earth curvature and
//! atmospheric refraction are applied to the target's apparent elevation when
//! `curvature` is on.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::{DataType, Raster};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, Layer};

use crate::args_common::{band_index, bool_or, choice_or, f64_or, opt_f64, req_str};
use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};
use crate::vector_common::{load_input_layer, write_or_store_layer};

/// ArcGIS caps the observer bit mask at 32; so does this, and for the same
/// reason — beyond that there are no bits left in the cell value.
const MAX_OBSERVERS: usize = 32;

/// Mean Earth radius in metres, for the curvature correction.
const EARTH_RADIUS: f64 = 6_371_000.0;

pub struct ObserverPointsTool;

impl Tool for ObserverPointsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "observer_points",
            display_name: "Observer Points",
            summary: "Identifies which of up to 32 observers can see each cell, as a per-observer bit mask, and writes the above-ground-level raster giving the height a target would need to be raised to become visible (ArcGIS Observer Points / Viewshed 2 in OBSERVERS mode). The bundled viewshed returns only a count of visible stations, from which neither answer can be recovered: a count of 1 does not say which tower covers a parcel, and a count of 0 does not distinguish two metres short from behind a mountain.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input elevation raster (DEM or DSM).",
                    required: true,
                },
                ToolParamSpec {
                    name: "observers",
                    description: "Point layer of observer positions, at most 32.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output visibility raster: an observer bit mask, or a count under analysis_type 'frequency'. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_agl",
                    description: "Output above-ground-level raster: the height a target must be raised to become visible. Always produced; stored in memory when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output_table",
                    description: "Output table with one row per observer and the cells and area it sees. Always produced; stored in memory when no path is given.",
                    required: false,
                },
                ToolParamSpec {
                    name: "analysis_type",
                    description: "'observers' (default; a bit mask) or 'frequency' (a count, matching the bundled viewshed).",
                    required: false,
                },
                ToolParamSpec {
                    name: "observer_offset",
                    description: "Height of each observer above the surface, in z units (default 2.0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "observer_offset_field",
                    description: "Attribute holding a per-observer height, overriding 'observer_offset'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "target_offset",
                    description: "Height of the target above the surface, in z units (default 0.0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_distance",
                    description: "Ignore cells farther than this from an observer, in map units. Unlimited by default.",
                    required: false,
                },
                ToolParamSpec {
                    name: "curvature",
                    description: "Apply the Earth-curvature and refraction correction (default false, i.e. a flat earth).",
                    required: false,
                },
                ToolParamSpec {
                    name: "refractivity_coefficient",
                    description: "Atmospheric refraction coefficient used with 'curvature' (default 0.13).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band of the elevation raster (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        req_str(args, "observers")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let dem_path = req_str(args, "input")?.to_string();
        let obs_path = req_str(args, "observers")?.to_string();
        let prm = parse_params(args)?;
        let band = band_index(args, "band")?;
        let output = parse_optional_output(args, "output")?;
        let out_agl = parse_optional_output(args, "output_agl")?;
        let out_table = parse_optional_output(args, "output_table")?;

        let dem = load_input_raster(&dem_path)?;
        let (rows, cols) = (dem.rows, dem.cols);
        let mut z = vec![f64::NAN; rows * cols];
        for r in 0..rows {
            for c in 0..cols {
                let v = dem.get(band, r as isize, c as isize);
                if v != dem.nodata && v.is_finite() {
                    z[r * cols + c] = v;
                }
            }
        }

        let layer = load_input_layer(&obs_path)?;
        let observers = collect_observers(&dem, &layer, &prm, band)?;
        if observers.is_empty() {
            return Err(ToolError::Validation(
                "no observer points fell inside the elevation raster".to_string(),
            ));
        }
        if observers.len() > MAX_OBSERVERS {
            return Err(ToolError::Validation(format!(
                "{} observers were supplied; at most {MAX_OBSERVERS} are allowed, because each \
                 gets one bit of the output cell value",
                observers.len()
            )));
        }

        ctx.progress.info(&format!(
            "{rows}x{cols}, {} observer(s), {} mode{}",
            observers.len(),
            prm.analysis.label(),
            if prm.curvature { ", curvature on" } else { "" }
        ));

        // Bit mask (or count) plus the running minimum required height.
        let nodata = -9999.0_f64;
        let mut mask = vec![0.0f64; rows * cols];
        let mut freq = vec![0.0f64; rows * cols];
        let mut agl = vec![f64::INFINITY; rows * cols];
        let mut per_observer = vec![0usize; observers.len()];

        for (k, obs) in observers.iter().enumerate() {
            // With a range limit, only the cells inside that radius can be
            // visible, so scanning the whole grid is wasted work — the full
            // sweep is O(observers x rows x cols x max(rows, cols)) and a
            // routine DEM makes the tool look hung.
            let (r0, r1, c0, c1) = match prm.max_distance {
                None => (0, rows, 0, cols),
                Some(d) => {
                    let dr = (d / dem.cell_size_y).ceil() as usize + 1;
                    let dc = (d / dem.cell_size_x).ceil() as usize + 1;
                    (
                        obs.row.saturating_sub(dr),
                        (obs.row + dr + 1).min(rows),
                        obs.col.saturating_sub(dc),
                        (obs.col + dc + 1).min(cols),
                    )
                }
            };
            let scanned_rows = r1.saturating_sub(r0).max(1);
            for r in r0..r1 {
                // Report inside the pass as well as between passes, so a long
                // single-observer sweep stays observable.
                if r % 64 == 0 {
                    ctx.progress.progress(
                        (k as f64 + (r - r0) as f64 / scanned_rows as f64)
                            / observers.len() as f64,
                    );
                }
                for c in c0..c1 {
                    let i = r * cols + c;
                    if !z[i].is_finite() {
                        continue;
                    }
                    let Some(req) = required_height(&z, rows, cols, &dem, obs, r, c, &prm) else {
                        continue;
                    };
                    if req <= 0.0 {
                        mask[i] += (1u64 << k) as f64;
                        freq[i] += 1.0;
                        per_observer[k] += 1;
                        agl[i] = 0.0;
                    } else if req < agl[i] {
                        agl[i] = req;
                    }
                }
            }
            ctx.progress
                .progress((k as f64 + 1.0) / observers.len() as f64);
        }

        // Cells with no elevation carry no answer at all.
        let mut visible_cells = 0usize;
        for i in 0..rows * cols {
            if !z[i].is_finite() {
                mask[i] = nodata;
                freq[i] = nodata;
                agl[i] = nodata;
                continue;
            }
            if freq[i] > 0.0 {
                visible_cells += 1;
            }
            if !agl[i].is_finite() {
                // Out of range of every observer, so no height would help.
                agl[i] = nodata;
            }
        }

        let (values, dtype) = match prm.analysis {
            Analysis::Observers => (mask, DataType::F64),
            Analysis::Frequency => (freq, DataType::F32),
        };
        let out_path = write_or_store_output(
            raster_like_with_data(&dem, values, nodata, dtype)?,
            output,
        )?;
        let agl_path = write_or_store_output(
            raster_like_with_data(&dem, agl, nodata, DataType::F32)?,
            out_agl,
        )?;

        let cell_area = dem.cell_size_x * dem.cell_size_y;
        let table = observer_table(&observers, &per_observer, cell_area)?;
        let table_path = write_or_store_layer(table, out_table)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("output_agl".to_string(), json!(agl_path));
        outputs.insert("output_table".to_string(), json!(table_path));
        outputs.insert("observer_count".to_string(), json!(observers.len()));
        outputs.insert("visible_cells".to_string(), json!(visible_cells));
        outputs.insert("analysis_type".to_string(), json!(prm.analysis.label()));
        outputs.insert(
            "cells_per_observer".to_string(),
            json!(per_observer),
        );
        Ok(ToolRunResult { outputs })
    }
}

/// One observer, resolved onto the raster grid.
struct Observer {
    row: usize,
    col: usize,
    /// Eye elevation: terrain plus the observer's own offset.
    eye: f64,
    x: f64,
    y: f64,
}

/// Resolves the observer layer onto the grid, dropping points that fall outside
/// or land on no-data.
fn collect_observers(
    dem: &Raster,
    layer: &Layer,
    prm: &Params,
    band: isize,
) -> Result<Vec<Observer>, ToolError> {
    let (rows, cols) = (dem.rows, dem.cols);
    let y_max = dem.y_min + rows as f64 * dem.cell_size_y;

    let offset_idx = match &prm.offset_field {
        None => None,
        Some(name) => Some(layer.schema.field_index(name).ok_or_else(|| {
            ToolError::Validation(format!(
                "'observer_offset_field' references field '{name}', which the observer layer does \
                 not have"
            ))
        })?),
    };

    let mut out = Vec::new();
    for f in layer.iter() {
        let Some(Geometry::Point(p)) = f.geometry.as_ref() else {
            continue;
        };
        let col_f = (p.x - dem.x_min) / dem.cell_size_x;
        let row_f = (y_max - p.y) / dem.cell_size_y;
        if col_f < 0.0 || row_f < 0.0 {
            continue;
        }
        let (col, row) = (col_f as usize, row_f as usize);
        if row >= rows || col >= cols {
            continue;
        }
        // The same band `run` built the elevation grid from; reading band 0
        // here would put the eye height on a different surface than the
        // line-of-sight walk.
        let ground = dem.get(band, row as isize, col as isize);
        if ground == dem.nodata || !ground.is_finite() {
            continue;
        }
        let offset = match offset_idx {
            None => prm.observer_offset,
            Some(i) => match &f.attributes[i] {
                FieldValue::Float(v) => *v,
                FieldValue::Integer(v) => *v as f64,
                // A row with no height falls back to the global default rather
                // than being dropped silently.
                _ => prm.observer_offset,
            },
        };
        out.push(Observer {
            row,
            col,
            eye: ground + offset,
            x: p.x,
            y: p.y,
        });
    }
    Ok(out)
}

/// How much a target at `(tr, tc)` must be raised to be visible from `obs`.
///
/// Returns `0.0` when the ground itself is already visible, a positive height
/// when it is not, or `None` when the cell is out of range or unusable.
///
/// The ray is walked in equal steps of about one cell, sampling the terrain
/// underneath and tracking the largest elevation angle seen so far. A target is
/// visible exactly when its own angle clears that horizon.
#[allow(clippy::too_many_arguments)]
fn required_height(
    z: &[f64],
    rows: usize,
    cols: usize,
    dem: &Raster,
    obs: &Observer,
    tr: usize,
    tc: usize,
    prm: &Params,
) -> Option<f64> {
    let (csx, csy) = (dem.cell_size_x, dem.cell_size_y);
    let dx = (tc as f64 - obs.col as f64) * csx;
    let dy = (tr as f64 - obs.row as f64) * csy;
    let dist = (dx * dx + dy * dy).sqrt();

    if let Some(max) = prm.max_distance {
        if dist > max {
            return None;
        }
    }
    let target_ground = z[tr * cols + tc];
    if !target_ground.is_finite() {
        return None;
    }
    // The observer's own cell is trivially visible.
    if dist == 0.0 {
        return Some(0.0);
    }

    // Apparent drop from curvature, offset by refraction bending the ray back
    // down toward the surface.
    let drop = |d: f64| -> f64 {
        if prm.curvature {
            (1.0 - prm.refractivity) * d * d / (2.0 * EARTH_RADIUS)
        } else {
            0.0
        }
    };

    // Step along the ray at roughly one cell at a time.
    let steps = ((tr as f64 - obs.row as f64)
        .abs()
        .max((tc as f64 - obs.col as f64).abs())
        .ceil() as usize)
        .max(1);

    let mut horizon = f64::NEG_INFINITY;
    for s in 1..steps {
        let t = s as f64 / steps as f64;
        let rr = obs.row as f64 + t * (tr as f64 - obs.row as f64);
        let cc = obs.col as f64 + t * (tc as f64 - obs.col as f64);
        let ri = rr.round() as isize;
        let ci = cc.round() as isize;
        if ri < 0 || ci < 0 || ri as usize >= rows || ci as usize >= cols {
            continue;
        }
        let zi = z[ri as usize * cols + ci as usize];
        if !zi.is_finite() {
            // A no-data gap blocks nothing; treating it as ground height zero
            // would invent a valley that lets the ray through.
            continue;
        }
        let d = t * dist;
        if d <= 0.0 {
            continue;
        }
        let angle = (zi - drop(d) - obs.eye) / d;
        if angle > horizon {
            horizon = angle;
        }
    }

    // Elevation the target must reach for its angle to clear the horizon.
    let needed_z = obs.eye + horizon.max(f64::NEG_INFINITY) * dist + drop(dist);
    let have = target_ground + prm.target_offset;
    if horizon == f64::NEG_INFINITY || have >= needed_z {
        Some(0.0)
    } else {
        Some(needed_z - have)
    }
}

/// One row per observer: its position and how much it sees.
fn observer_table(
    observers: &[Observer],
    counts: &[usize],
    cell_area: f64,
) -> Result<Layer, ToolError> {
    let mut layer = Layer::new("observer_visibility");
    layer.add_field(FieldDef::new("observer", FieldType::Integer));
    layer.add_field(FieldDef::new("bit", FieldType::Integer));
    layer.add_field(FieldDef::new("x", FieldType::Float));
    layer.add_field(FieldDef::new("y", FieldType::Float));
    layer.add_field(FieldDef::new("eye_z", FieldType::Float));
    layer.add_field(FieldDef::new("visible_cells", FieldType::Integer));
    layer.add_field(FieldDef::new("visible_area", FieldType::Float));

    for (k, obs) in observers.iter().enumerate() {
        layer
            .add_feature(
                None,
                &[
                    ("observer", FieldValue::Integer(k as i64)),
                    ("bit", FieldValue::Integer(1i64 << k)),
                    ("x", FieldValue::Float(obs.x)),
                    ("y", FieldValue::Float(obs.y)),
                    ("eye_z", FieldValue::Float(obs.eye)),
                    ("visible_cells", FieldValue::Integer(counts[k] as i64)),
                    (
                        "visible_area",
                        FieldValue::Float(counts[k] as f64 * cell_area),
                    ),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("writing the observer table: {e}")))?;
    }
    Ok(layer)
}

// ── Parameters ──────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Analysis {
    Observers,
    Frequency,
}

impl Analysis {
    fn label(self) -> &'static str {
        match self {
            Analysis::Observers => "observers",
            Analysis::Frequency => "frequency",
        }
    }
}

struct Params {
    analysis: Analysis,
    observer_offset: f64,
    offset_field: Option<String>,
    target_offset: f64,
    max_distance: Option<f64>,
    curvature: bool,
    refractivity: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let analysis = match choice_or(
        args,
        "analysis_type",
        &["observers", "frequency"],
        "observers",
    )? {
        "frequency" => Analysis::Frequency,
        _ => Analysis::Observers,
    };
    let observer_offset = f64_or(args, "observer_offset", 2.0)?;
    if !observer_offset.is_finite() {
        return Err(ToolError::Validation(
            "'observer_offset' must be finite".to_string(),
        ));
    }
    let offset_field = args
        .get("observer_offset_field")
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(String::from);
    let target_offset = f64_or(args, "target_offset", 0.0)?;
    if !target_offset.is_finite() {
        return Err(ToolError::Validation(
            "'target_offset' must be finite".to_string(),
        ));
    }
    let max_distance = match opt_f64(args, "max_distance")? {
        None => None,
        Some(v) if v > 0.0 && v.is_finite() => Some(v),
        Some(v) => {
            return Err(ToolError::Validation(format!(
                "'max_distance' must be positive, got {v}"
            )))
        }
    };
    let curvature = bool_or(args, "curvature", false)?;
    let refractivity = f64_or(args, "refractivity_coefficient", 0.13)?;
    if !(0.0..1.0).contains(&refractivity) {
        return Err(ToolError::Validation(format!(
            "'refractivity_coefficient' must be in [0, 1), got {refractivity}"
        )));
    }

    Ok(Params {
        analysis,
        observer_offset,
        offset_field,
        target_offset,
        max_distance,
        curvature,
        refractivity,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, RasterConfig};
    use wbvector::{Coord, Feature, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn dem_of(cols: usize, rows: usize, vals: &[f64]) -> String {
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

    /// Observer points at the centre of the given (row, col) cells of a
    /// `rows`-high grid with 10 m cells.
    fn observers_at(rows: usize, cells: &[(usize, usize)]) -> String {
        let mut layer = Layer::new("observers");
        layer.geom_type = Some(GeometryType::Point);
        layer = layer.with_crs_epsg(32610);
        layer.add_field(FieldDef::new("id", FieldType::Integer));
        let y_max = rows as f64 * 10.0;
        for (i, &(r, c)) in cells.iter().enumerate() {
            let x = (c as f64 + 0.5) * 10.0;
            let y = y_max - (r as f64 + 0.5) * 10.0;
            let mut f = Feature::with_geometry(
                i as u64,
                Geometry::Point(Coord::xy(x, y)),
                layer.schema.len(),
            );
            f.set_by_index(0, FieldValue::Integer(i as i64));
            layer.push(f);
        }
        write_or_store_layer(layer, None).unwrap()
    }

    fn run(args: Value) -> (Raster, BTreeMap<String, Value>) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ObserverPointsTool.run(&args, &ctx()).unwrap();
        let r = load_input_raster(out.outputs["output"].as_str().unwrap()).unwrap();
        (r, out.outputs)
    }

    /// The capability the bundled `viewshed` lacks: on a flat plain split by a
    /// ridge, the mask says *which* observer sees each cell, so the two sides
    /// are distinguishable. A frequency raster would read 1 on both.
    #[test]
    fn mask_identifies_which_observer_sees_what() {
        let (rows, cols) = (1, 11);
        // Flat ground with a tall wall in the middle cell.
        let mut z = vec![0.0; cols];
        z[5] = 100.0;
        let dem = dem_of(cols, rows, &z);
        // One observer on each side of the wall.
        let obs = observers_at(rows, &[(0, 0), (0, 10)]);

        let (out, outputs) = run(json!({
            "input": dem, "observers": obs, "observer_offset": 2.0
        }));
        assert_eq!(outputs["observer_count"].as_u64().unwrap(), 2);

        // Bit 0 = the west observer, bit 1 = the east one.
        let west_side = out.get(0, 0, 2) as u64;
        let east_side = out.get(0, 0, 8) as u64;
        assert_eq!(west_side & 1, 1, "west cell should be seen by observer 0");
        assert_eq!(
            west_side & 2,
            0,
            "the wall should hide the west cell from observer 1"
        );
        assert_eq!(east_side & 2, 2, "east cell should be seen by observer 1");
        assert_eq!(
            east_side & 1,
            0,
            "the wall should hide the east cell from observer 0"
        );

        // The frequency form collapses exactly this distinction — which is why
        // the bundled viewshed cannot answer the question.
        let (freq, _) = run(json!({
            "input": dem_of(cols, rows, &z), "observers": observers_at(rows, &[(0, 0), (0, 10)]),
            "analysis_type": "frequency"
        }));
        assert_eq!(freq.get(0, 0, 2), 1.0);
        assert_eq!(freq.get(0, 0, 8), 1.0);
    }

    /// The second capability: AGL says how much higher a hidden target would
    /// have to be, which a count of zero cannot express.
    #[test]
    fn agl_gives_the_height_needed_to_be_seen() {
        let (rows, cols) = (1, 11);
        let mut z = vec![0.0; cols];
        z[5] = 100.0; // a 100 m wall at x = 55 m
        let dem = dem_of(cols, rows, &z);
        let obs = observers_at(rows, &[(0, 0)]);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": dem, "observers": obs, "observer_offset": 0.0
        }))
        .unwrap();
        let out = ObserverPointsTool.run(&args, &ctx()).unwrap();
        let agl = load_input_raster(out.outputs["output_agl"].as_str().unwrap()).unwrap();

        // Visible ground has zero requirement.
        assert_eq!(agl.get(0, 0, 2), 0.0, "cell before the wall is visible");

        // Behind the wall the sightline from the observer's eye (0 m at x = 5)
        // over the wall top (100 m at x = 55) rises at 2 m per metre, so at
        // x = 105 m (cell 10) the ray is at 200 m.
        let need = agl.get(0, 0, 10) as f64;
        assert!(
            (need - 200.0).abs() < 1.0,
            "expected about 200 m of required height behind the wall, got {need}"
        );
    }

    /// Flat ground with nothing to block: everyone sees everything.
    #[test]
    fn flat_terrain_is_fully_visible() {
        let (rows, cols) = (5, 5);
        let dem = dem_of(cols, rows, &[10.0; 25]);
        let obs = observers_at(rows, &[(2, 2)]);
        let (out, outputs) = run(json!({ "input": dem, "observers": obs }));
        assert_eq!(outputs["visible_cells"].as_u64().unwrap(), 25);
        for r in 0..5 {
            for c in 0..5 {
                assert_eq!(out.get(0, r, c), 1.0, "cell ({r},{c}) should be visible");
            }
        }
    }

    /// A raised observer sees over an obstacle that blocks a ground-level one.
    #[test]
    fn observer_offset_changes_what_is_visible() {
        let (rows, cols) = (1, 9);
        let mut z = vec![0.0; cols];
        z[4] = 20.0;
        let build = |offset: f64| {
            let args: ToolArgs = serde_json::from_value(json!({
                "input": dem_of(cols, rows, &z),
                "observers": observers_at(rows, &[(0, 0)]),
                "observer_offset": offset
            }))
            .unwrap();
            let out = ObserverPointsTool.run(&args, &ctx()).unwrap();
            out.outputs["visible_cells"].as_u64().unwrap()
        };
        let low = build(0.0);
        let high = build(1000.0);
        assert!(
            high > low,
            "a much higher observer must see more, got {low} then {high}"
        );
    }

    /// `max_distance` limits the analysis and leaves far cells unanswered.
    #[test]
    fn max_distance_limits_the_search() {
        let (rows, cols) = (1, 11);
        let dem = dem_of(cols, rows, &[0.0; 11]);
        let obs = observers_at(rows, &[(0, 0)]);
        let (out, outputs) = run(json!({
            "input": dem, "observers": obs, "max_distance": 35.0
        }));
        // Cells 0..3 are within 35 m of the observer at x = 5 m.
        assert_eq!(out.get(0, 0, 0), 1.0);
        assert_eq!(out.get(0, 0, 3), 1.0);
        assert_eq!(out.get(0, 0, 9), 0.0, "cell beyond the range is not visible");
        assert!(outputs["visible_cells"].as_u64().unwrap() < 11);
    }

    /// Per-observer heights can come from an attribute.
    #[test]
    fn offset_field_overrides_the_default() {
        let (rows, cols) = (1, 9);
        let mut z = vec![0.0; cols];
        z[4] = 20.0;
        let mut layer = Layer::new("observers");
        layer.geom_type = Some(GeometryType::Point);
        layer = layer.with_crs_epsg(32610);
        layer.add_field(FieldDef::new("id", FieldType::Integer));
        layer.add_field(FieldDef::new("mast", FieldType::Float));
        let mut f = Feature::with_geometry(
            0,
            Geometry::Point(Coord::xy(5.0, 5.0)),
            layer.schema.len(),
        );
        f.set_by_index(0, FieldValue::Integer(0));
        f.set_by_index(1, FieldValue::Float(1000.0));
        layer.push(f);
        let obs = write_or_store_layer(layer, None).unwrap();

        let (_, tall) = run(json!({
            "input": dem_of(cols, rows, &z), "observers": obs.clone(),
            "observer_offset": 0.0, "observer_offset_field": "mast"
        }));
        let (_, short) = run(json!({
            "input": dem_of(cols, rows, &z), "observers": obs, "observer_offset": 0.0
        }));
        assert!(
            tall["visible_cells"].as_u64().unwrap() > short["visible_cells"].as_u64().unwrap(),
            "the mast height from the attribute was not applied"
        );
    }

    /// The observer table records per-observer coverage.
    #[test]
    fn table_reports_per_observer_coverage() {
        let (rows, cols) = (3, 3);
        let dem = dem_of(cols, rows, &[0.0; 9]);
        let obs = observers_at(rows, &[(0, 0), (2, 2)]);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": dem, "observers": obs })).unwrap();
        let out = ObserverPointsTool.run(&args, &ctx()).unwrap();
        let table = load_input_layer(out.outputs["output_table"].as_str().unwrap()).unwrap();
        assert_eq!(table.len(), 2);
        let bi = table.schema.field_index("bit").unwrap();
        let ci = table.schema.field_index("visible_cells").unwrap();
        let ai = table.schema.field_index("visible_area").unwrap();
        let rows_out: Vec<&Feature> = table.iter().collect();
        assert!(matches!(rows_out[0].attributes[bi], FieldValue::Integer(1)));
        assert!(matches!(rows_out[1].attributes[bi], FieldValue::Integer(2)));
        // Flat ground: each observer sees all nine cells of 100 m^2.
        assert!(matches!(rows_out[0].attributes[ci], FieldValue::Integer(9)));
        let FieldValue::Float(area) = rows_out[0].attributes[ai] else {
            panic!("area must be a float")
        };
        assert!((area - 900.0).abs() < 1e-6);
    }

    /// Curvature lowers a distant target's apparent elevation, so switching it
    /// on can only reduce what is visible over a long flat sightline.
    #[test]
    fn curvature_reduces_long_range_visibility() {
        // 200 cells of 10 m = a 2 km line; the curvature drop at 2 km is about
        // 0.27 m, so a flat surface at exactly eye level goes out of sight.
        let (rows, cols) = (1, 200);
        let z = vec![0.0; cols];
        let build = |curv: bool| {
            let args: ToolArgs = serde_json::from_value(json!({
                "input": dem_of(cols, rows, &z),
                "observers": observers_at(rows, &[(0, 0)]),
                "observer_offset": 0.0, "curvature": curv
            }))
            .unwrap();
            ObserverPointsTool
                .run(&args, &ctx())
                .unwrap()
                .outputs["visible_cells"]
                .as_u64()
                .unwrap()
        };
        assert!(
            build(true) < build(false),
            "curvature must cut off the far end of a flat 2 km sightline"
        );
    }

    /// Beyond 32 observers there are no bits left; fail rather than silently
    /// wrapping two observers onto the same bit.
    #[test]
    fn rejects_more_than_32_observers() {
        let (rows, cols) = (1, 40);
        let dem = dem_of(cols, rows, &[0.0; 40]);
        let cells: Vec<(usize, usize)> = (0..33).map(|c| (0, c)).collect();
        let obs = observers_at(rows, &cells);
        let args: ToolArgs =
            serde_json::from_value(json!({ "input": dem, "observers": obs })).unwrap();
        let err = ObserverPointsTool.run(&args, &ctx()).unwrap_err();
        assert!(
            format!("{err:?}").contains("at most 32"),
            "expected an observer-count error, got {err:?}"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ObserverPointsTool.validate(&args)
        };
        assert!(bad(json!({"observers": "o.shp"})).is_err());
        assert!(bad(json!({"input": "d.tif"})).is_err());
        let base = json!({"input": "d.tif", "observers": "o.shp"});
        assert!(bad(base.clone()).is_ok());
        let with = |k: &str, v: Value| {
            let mut m = base.as_object().unwrap().clone();
            m.insert(k.into(), v);
            Value::Object(m)
        };
        assert!(bad(with("analysis_type", json!("agl"))).is_err());
        assert!(bad(with("max_distance", json!(-1))).is_err());
        assert!(bad(with("refractivity_coefficient", json!(1.5))).is_err());
        assert!(bad(with("analysis_type", json!("frequency"))).is_ok());
    }
}
