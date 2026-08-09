//! GeoLibre tool: 3D shadow volumes cast by solids for a given sun position.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Sun Shadow Volume* (3D Analyst).
//!
//! ## The gap
//!
//! `shadow_image`, `shadow_animation` and `time_in_daylight` all operate on a
//! **raster DSM** and emit rasters; `solar_radiation` (round 2) computes
//! irradiance. None of them produces shadow *geometry*.
//!
//! Geometry is what urban planning and solar-access analysis actually need,
//! because the question is rarely "is this ground cell shaded" but "is this
//! window shaded", "does this new tower shade the park", "how much of this
//! plot's airspace is in shadow at 3pm in December". Those are 3D containment
//! and overlay questions, and the answer composes with `inside_3d`,
//! `intersect_3d` and `union_3d` — none of which a shadow raster can feed.
//!
//! ## Construction
//!
//! For each solid, the **silhouette** edges against the sun direction (edges
//! whose two adjoining faces disagree on whether they face the sun) are swept
//! along the sun vector down to the ground plane. That sweep, plus the original
//! silhouette and its ground-plane projection, closes the shadow volume.
//!
//! Using the silhouette rather than every edge is what keeps the result a
//! single clean prism instead of a self-intersecting tangle of per-face sweeps.
//!
//! ## Determinism
//!
//! The sun position comes from the `datetime` **parameter** — never from a
//! clock — so the tool is deterministic and WASM-safe. No RNG is involved.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, GeometryType, Layer};

use crate::args_common::{bool_or, f64_or, opt_f64, req_str};
use crate::inside_3d::{collect_triangles, Tri};
use crate::mesh3d::{edge_map, mesh_volume, topology, tri_normal, triangles_to_geometry};
use crate::solar_radiation::{declination, sun_position};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct SunShadowVolumeTool;

impl Tool for SunShadowVolumeTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "sun_shadow_volume",
            display_name: "Sun Shadow Volume",
            summary: "Extrudes 3D solids along the solar vector to the ground to produce closed shadow-volume solids for a given date and time (ArcGIS Sun Shadow Volume). shadow_image, shadow_animation and time_in_daylight all work on a raster DSM and emit rasters, and solar_radiation computes irradiance — none produces shadow geometry, so questions like 'is this window shaded' or 'how much of this plot's airspace is in shadow' cannot be answered. The output composes with inside_3d, intersect_3d and union_3d.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "3D building solids or extruded footprints (triangle-mesh MultiPolygons with Z).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Closed shadow-volume solids. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "datetime",
                    description: "Local date-time of the sun position as 'YYYY-MM-DDTHH:MM' or 'YYYY-MM-DD HH:MM'. Taken as a parameter, never from a clock, so results are reproducible.",
                    required: true,
                },
                ToolParamSpec {
                    name: "latitude",
                    description: "Latitude in degrees, used for the solar geometry. Required: it cannot be derived without a projection engine.",
                    required: true,
                },
                ToolParamSpec {
                    name: "utc_offset",
                    description: "Hours from UTC of the supplied local time (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "adjusted_for_dst",
                    description: "Treat the supplied time as daylight-saving-adjusted, subtracting an hour before the solar calculation (default false), matching ArcGIS.",
                    required: false,
                },
                ToolParamSpec {
                    name: "ground_elevation",
                    description: "Constant ground Z the shadow is projected onto. Default: each feature's own minimum Z.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_length",
                    description: "Cap on the sweep distance in CRS units, so a low sun does not produce a kilometres-long prism. Default: unlimited.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        parse_datetime(req_str(args, "datetime")?)?;
        let lat = opt_f64(args, "latitude")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'latitude' (degrees)".to_string())
        })?;
        if !(-90.0..=90.0).contains(&lat) {
            return Err(ToolError::Validation(
                "'latitude' must be in [-90, 90]".to_string(),
            ));
        }
        crate::args_common::opt_positive_f64(args, "max_length")?;
        bool_or(args, "adjusted_for_dst", false)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let (day_of_year, hour) = parse_datetime(req_str(args, "datetime")?)?;
        let latitude = opt_f64(args, "latitude")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'latitude'".to_string())
        })?;
        let utc_offset = f64_or(args, "utc_offset", 0.0)?;
        let dst = bool_or(args, "adjusted_for_dst", false)?;
        let ground_param = opt_f64(args, "ground_elevation")?;
        let max_length =
            crate::args_common::opt_positive_f64(args, "max_length")?.unwrap_or(f64::INFINITY);

        // Convert the supplied local clock time to solar hour. utc_offset moves
        // it to UTC; the DST adjustment undoes the civil hour shift.
        let solar_hour = hour - utc_offset - if dst { 1.0 } else { 0.0 } + longitude_hour(latitude);
        let (altitude, azimuth) =
            sun_position(latitude.to_radians(), declination(day_of_year), solar_hour);

        if altitude <= 0.0 {
            return Err(ToolError::Execution(format!(
                "the sun is below the horizon at that date and time (altitude {:.2} degrees); \
                 there is no shadow volume to build",
                altitude.to_degrees()
            )));
        }

        // Direction light travels: horizontally away from the sun's bearing,
        // and downward at the sun's altitude.
        let dir = [
            -azimuth.sin() * altitude.cos(),
            -azimuth.cos() * altitude.cos(),
            -altitude.sin(),
        ];
        ctx.progress.info(&format!(
            "sun altitude {:.2} deg, azimuth {:.2} deg",
            altitude.to_degrees(),
            azimuth.to_degrees()
        ));

        let layer = load_input_layer(input)?;
        let mut out = Layer::new("sun_shadow_volume");
        out.geom_type = Some(GeometryType::MultiPolygon);
        out.crs = layer.crs.clone();
        out.add_field(FieldDef::new("SRC_FID", FieldType::Integer));
        out.add_field(FieldDef::new("VOLUME", FieldType::Float));
        out.add_field(FieldDef::new("SUN_AZIMUTH", FieldType::Float));
        out.add_field(FieldDef::new("SUN_ALTITUDE", FieldType::Float));
        out.add_field(FieldDef::new("GROUND_Z", FieldType::Float));
        out.add_field(FieldDef::new("WATERTIGHT", FieldType::Boolean));

        let mut built = 0_u64;
        let mut skipped = 0_u64;
        let mut total_volume = 0.0_f64;
        let total = layer.iter().count().max(1);

        for (fid, feature) in layer.iter().enumerate() {
            let tris = feature
                .geometry
                .as_ref()
                .map(collect_triangles)
                .unwrap_or_default();
            if tris.is_empty() {
                skipped += 1;
                continue;
            }
            let ground = ground_param.unwrap_or_else(|| {
                tris.iter()
                    .flatten()
                    .map(|v| v[2])
                    .fold(f64::INFINITY, f64::min)
            });

            let Some(shadow) = sweep(&tris, dir, ground, max_length) else {
                skipped += 1;
                continue;
            };
            // The sweep is closed by construction (lit cap, projected floor,
            // swept silhouette wall), so the signed-tetrahedron sum is the
            // shadow volume by the divergence theorem. Watertightness is
            // reported as a diagnostic rather than used as a gate: a sliver
            // dropped from a near-degenerate face should not silently turn a
            // real shadow into a reported volume of zero.
            let t = topology(&shadow);
            let volume = mesh_volume(&shadow);
            total_volume += volume;
            built += 1;

            out.add_feature(
                Some(triangles_to_geometry(&shadow)),
                &[
                    ("SRC_FID", FieldValue::Integer(fid as i64)),
                    ("VOLUME", FieldValue::Float(volume)),
                    ("SUN_AZIMUTH", FieldValue::Float(azimuth.to_degrees())),
                    ("SUN_ALTITUDE", FieldValue::Float(altitude.to_degrees())),
                    ("GROUND_Z", FieldValue::Float(ground)),
                    ("WATERTIGHT", FieldValue::Boolean(t.closed)),
                ],
            )
            .map_err(|e| ToolError::Execution(e.to_string()))?;
            ctx.progress.progress((fid as f64 + 1.0) / total as f64);
        }

        if built == 0 {
            return Err(ToolError::Execution(format!(
                "no shadow volume could be built ({skipped} feature(s) skipped)"
            )));
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("shadow_count".to_string(), json!(built));
        outputs.insert("skipped_count".to_string(), json!(skipped));
        outputs.insert("total_volume".to_string(), json!(total_volume));
        outputs.insert("sun_azimuth".to_string(), json!(azimuth.to_degrees()));
        outputs.insert("sun_altitude".to_string(), json!(altitude.to_degrees()));
        Ok(ToolRunResult { outputs })
    }
}

/// Sweeps a solid's silhouette along `dir` to the `ground` plane.
///
/// Returns a closed mesh: the sun-facing faces (the cap), the swept silhouette
/// wall, and the ground-plane projection (the floor). Faces pointing away from
/// the sun are dropped — they are inside the shadow volume, not on its
/// boundary — which is precisely what makes the result a single prism.
fn sweep(tris: &[Tri], dir: [f64; 3], ground: f64, max_length: f64) -> Option<Vec<Tri>> {
    // Project a point along `dir` until it reaches the ground plane.
    let project = |p: [f64; 3]| -> [f64; 3] {
        if dir[2] >= -1e-12 {
            return p; // sun at or below the horizon; caller already guarded
        }
        let t = ((p[2] - ground) / -dir[2]).max(0.0).min(max_length);
        [p[0] + dir[0] * t, p[1] + dir[1] * t, p[2] + dir[2] * t]
    };

    // Sun-facing test per triangle: the face is lit when its outward normal
    // opposes the light direction.
    let lit: Vec<bool> = tris
        .iter()
        .map(|t| tri_normal(t).is_some_and(|n| dot(n, dir) < 0.0))
        .collect();
    if !lit.iter().any(|b| *b) {
        return None;
    }

    let coincident = |a: [f64; 3], b: [f64; 3]| {
        (a[0] - b[0]).abs() < 1e-9 && (a[1] - b[1]).abs() < 1e-9 && (a[2] - b[2]).abs() < 1e-9
    };

    let mut out: Vec<Tri> = Vec::new();
    // Counts geometry that was actually displaced. Without any, the "shadow"
    // is an open cap, and the signed-tetrahedron sum of an open mesh is an
    // artefact rather than a volume — so that case reports no shadow at all.
    let mut swept = 0usize;
    // Cap: the lit faces themselves.
    for (i, t) in tris.iter().enumerate() {
        if lit[i] {
            out.push(*t);
        }
    }
    // Floor: the lit faces projected to the ground, wound the other way so the
    // prism closes. A lit face already lying ON the ground projects onto
    // itself, which would duplicate the cap and break edge pairing, so it
    // contributes no floor triangle — the cap already bounds that side.
    for (i, t) in tris.iter().enumerate() {
        if !lit[i] {
            continue;
        }
        let p = [project(t[0]), project(t[1]), project(t[2])];
        if (0..3).all(|k| coincident(p[k], t[k])) {
            continue;
        }
        out.push([p[0], p[2], p[1]]);
        swept += 1;
    }
    // Wall: sweep every silhouette edge — an edge with exactly one lit
    // neighbour. Interior lit-lit and unlit-unlit edges are not on the
    // shadow's boundary and must not be swept, or the prism self-intersects.
    for (_key, uses) in edge_map(tris) {
        let lit_count = uses.iter().filter(|u| lit[u.tri]).count();
        let is_silhouette = match uses.len() {
            1 => lit_count == 1,
            2 => lit_count == 1,
            _ => false,
        };
        if !is_silhouette {
            continue;
        }
        // Take the edge as the LIT triangle walks it, so the wall's winding
        // agrees with the cap it grows from.
        let Some(u) = uses.iter().find(|u| lit[u.tri]) else {
            continue;
        };
        let t = &tris[u.tri];
        let mut edge = None;
        for k in 0..3 {
            let (p, q) = (t[k], t[(k + 1) % 3]);
            if crate::mesh3d::edge_key(p, q) == _key {
                edge = Some((p, q));
                break;
            }
        }
        let Some((p, q)) = edge else { continue };
        let (pp, qp) = (project(p), project(q));
        // An edge already sitting on the ground sweeps to itself: the wall
        // quad would be two zero-area triangles carrying a degenerate p-p
        // edge, which breaks the edge map. Skipping it is not a hole — the cap
        // edge and the floor edge are then the same edge and pair directly.
        if coincident(p, pp) && coincident(q, qp) {
            continue;
        }
        out.push([p, pp, qp]);
        out.push([p, qp, q]);
        swept += 1;
    }
    // Drop any remaining slivers so they cannot leave unpaired edges behind.
    out.retain(|t| crate::mesh3d::tri_area(t) > 1e-12);
    (swept > 0 && !out.is_empty()).then_some(out)
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Longitude correction is not available without a projection engine, so the
/// supplied time is treated as *local solar* time once the UTC offset and DST
/// are removed. Returning zero here documents that explicitly rather than
/// silently pretending to a precision the inputs do not support.
fn longitude_hour(_latitude: f64) -> f64 {
    0.0
}

/// Parses `YYYY-MM-DD[T ]HH:MM[:SS]` into (day of year, fractional hour).
fn parse_datetime(s: &str) -> Result<(i64, f64), ToolError> {
    let bad = || {
        ToolError::Validation(format!(
            "'datetime' must look like 'YYYY-MM-DDTHH:MM' (got '{s}')"
        ))
    };
    let (date, time) = s
        .split_once('T')
        .or_else(|| s.split_once(' '))
        .ok_or_else(bad)?;
    let d: Vec<&str> = date.split('-').collect();
    if d.len() != 3 {
        return Err(bad());
    }
    let year: i64 = d[0].parse().map_err(|_| bad())?;
    let month: i64 = d[1].parse().map_err(|_| bad())?;
    let day: i64 = d[2].parse().map_err(|_| bad())?;
    if !(1..=12).contains(&month) {
        return Err(bad());
    }
    let leap = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
    let month_len = [
        31,
        if leap { 29 } else { 28 },
        31,
        30,
        31,
        30,
        31,
        31,
        30,
        31,
        30,
        31,
    ][(month - 1) as usize];
    if !(1..=month_len).contains(&day) {
        return Err(bad());
    }

    let t: Vec<&str> = time.split(':').collect();
    if t.len() < 2 {
        return Err(bad());
    }
    let hh: f64 = t[0].parse().map_err(|_| bad())?;
    let mm: f64 = t[1].parse().map_err(|_| bad())?;
    let ss: f64 = if t.len() > 2 {
        t[2].parse().unwrap_or(0.0)
    } else {
        0.0
    };
    if !(0.0..24.0).contains(&hh) || !(0.0..60.0).contains(&mm) {
        return Err(bad());
    }

    let cumulative = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    let mut doy = cumulative[(month - 1) as usize] + day;
    if leap && month > 2 {
        doy += 1;
    }
    Ok((doy, hh + mm / 60.0 + ss / 3600.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Geometry};

    use crate::mesh3d::box_mesh;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn layer_of(geoms: Vec<Geometry>) -> String {
        let mut l = Layer::new("in");
        l.geom_type = Some(GeometryType::MultiPolygon);
        for g in geoms {
            l.add_feature(Some(g), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = SunShadowVolumeTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    fn num(layer: &Layer, fid: usize, name: &str) -> f64 {
        let i = layer.schema.field_index(name).unwrap();
        match &layer.iter().nth(fid).unwrap().attributes[i] {
            FieldValue::Float(v) => *v,
            FieldValue::Integer(v) => *v as f64,
            other => panic!("expected a number, got {other:?}"),
        }
    }

    /// Midsummer noon at mid-latitude: a high sun.
    fn noon() -> Value {
        json!({"datetime": "2026-06-21T12:00", "latitude": 45.0})
    }

    #[test]
    fn a_building_casts_a_shadow_volume_with_positive_extent() {
        let mut a = noon();
        a["input"] = json!(layer_of(vec![box_mesh(
            [0.0, 0.0, 0.0],
            [10.0, 10.0, 20.0]
        )]));
        let (out, res) = run(a);
        assert_eq!(res.outputs["shadow_count"], json!(1));
        assert!(num(&out, 0, "VOLUME") > 0.0);
        assert!(num(&out, 0, "SUN_ALTITUDE") > 0.0);
    }

    #[test]
    fn the_shadow_geometry_extends_below_the_building_top() {
        // A shadow must reach the ground plane, not float at roof level.
        let mut a = noon();
        a["input"] = json!(layer_of(vec![box_mesh(
            [0.0, 0.0, 0.0],
            [10.0, 10.0, 20.0]
        )]));
        let (out, _) = run(a);
        let geom = out.iter().next().unwrap().geometry.clone().unwrap();
        let zs: Vec<f64> = collect_triangles(&geom)
            .iter()
            .flatten()
            .map(|v| v[2])
            .collect();
        let min_z = zs.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(
            min_z <= 1e-6,
            "shadow bottoms out at {min_z}, not the ground"
        );
    }

    #[test]
    fn an_overhead_sun_gives_exactly_the_vertical_prism_under_the_roof() {
        // The one case with a closed-form answer. On the June solstice at
        // latitude 23.45 (the Tropic of Cancer) solar noon puts the sun at the
        // zenith, so the shadow is the vertical prism under the roof:
        // 10 x 10 footprint x 20 height = 2000.
        let mut a = json!({"datetime": "2026-06-21T12:00", "latitude": 23.45});
        a["input"] = json!(layer_of(vec![box_mesh(
            [0.0, 0.0, 0.0],
            [10.0, 10.0, 20.0]
        )]));
        let (out, res) = run(a);
        assert!(
            res.outputs["sun_altitude"].as_f64().unwrap() > 89.0,
            "sun was not near the zenith: {:?}",
            res.outputs["sun_altitude"]
        );
        let v = num(&out, 0, "VOLUME");
        assert!((v - 2000.0).abs() < 1.0, "expected 2000, got {v}");
    }

    #[test]
    fn a_ground_plane_at_roof_level_yields_no_shadow_at_all() {
        // Projecting onto the roof itself is a zero-length sweep. The right
        // answer is "no shadow", not "an open cap with an artefact volume":
        // the signed-tetrahedron sum only means something on a closed mesh.
        // Uses the zenith sun so the roof is the only lit face.
        let mut a = json!({"datetime": "2026-06-21T12:00", "latitude": 23.45});
        a["input"] = json!(layer_of(vec![box_mesh(
            [0.0, 0.0, 0.0],
            [10.0, 10.0, 20.0]
        )]));
        a["ground_elevation"] = json!(20.0);
        let args: ToolArgs = serde_json::from_value(a).unwrap();
        assert!(SunShadowVolumeTool.run(&args, &ctx()).is_err());
    }

    #[test]
    fn a_lower_sun_casts_a_longer_shadow() {
        // The core physical behaviour. Late afternoon vs noon.
        let vol = |time: &str| {
            let mut a = json!({"datetime": time, "latitude": 45.0});
            a["input"] = json!(layer_of(vec![box_mesh(
                [0.0, 0.0, 0.0],
                [10.0, 10.0, 20.0]
            )]));
            let (out, _) = run(a);
            (num(&out, 0, "VOLUME"), num(&out, 0, "SUN_ALTITUDE"))
        };
        let (v_noon, alt_noon) = vol("2026-06-21T12:00");
        let (v_late, alt_late) = vol("2026-06-21T17:00");
        assert!(alt_late < alt_noon, "sun did not descend");
        assert!(
            v_late > v_noon,
            "lower sun ({alt_late:.1} deg) gave a smaller shadow: {v_late} vs {v_noon}"
        );
    }

    #[test]
    fn the_sun_below_the_horizon_is_an_explicit_error() {
        let mut a = json!({"datetime": "2026-12-21T23:00", "latitude": 60.0});
        a["input"] = json!(layer_of(vec![box_mesh([0.0; 3], [1.0, 1.0, 1.0])]));
        let args: ToolArgs = serde_json::from_value(a).unwrap();
        let err = SunShadowVolumeTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err:?}").contains("horizon"), "got {err:?}");
    }

    #[test]
    fn max_length_caps_the_sweep() {
        let vol = |cap: Option<f64>| {
            let mut a = json!({"datetime": "2026-06-21T17:30", "latitude": 55.0});
            a["input"] = json!(layer_of(vec![box_mesh(
                [0.0, 0.0, 0.0],
                [10.0, 10.0, 30.0]
            )]));
            if let Some(c) = cap {
                a["max_length"] = json!(c);
            }
            let (out, _) = run(a);
            num(&out, 0, "VOLUME")
        };
        assert!(vol(Some(5.0)) < vol(None));
    }

    #[test]
    fn a_ground_elevation_raises_the_shadow_floor() {
        // Projecting onto a higher ground plane means a shorter sweep.
        let vol = |g: f64| {
            let mut a = noon();
            a["input"] = json!(layer_of(vec![box_mesh(
                [0.0, 0.0, 0.0],
                [10.0, 10.0, 20.0]
            )]));
            a["ground_elevation"] = json!(g);
            let (out, _) = run(a);
            num(&out, 0, "VOLUME")
        };
        assert!(vol(15.0) < vol(0.0));
    }

    #[test]
    fn results_are_reproducible_for_the_same_parameters() {
        // No clock, no RNG: the same inputs must give the same answer.
        let go = || {
            let mut a = noon();
            a["input"] = json!(layer_of(vec![box_mesh([0.0, 0.0, 0.0], [7.0, 3.0, 12.0])]));
            let (out, _) = run(a);
            num(&out, 0, "VOLUME")
        };
        assert_eq!(go(), go());
    }

    #[test]
    fn datetime_parsing_accepts_both_separators_and_rejects_junk() {
        assert!(parse_datetime("2026-06-21T12:00").is_ok());
        assert!(parse_datetime("2026-06-21 12:00:30").is_ok());
        // 2024 is a leap year, so 1 March is day 61 rather than 60.
        assert_eq!(parse_datetime("2024-03-01T00:00").unwrap().0, 61);
        assert_eq!(parse_datetime("2026-03-01T00:00").unwrap().0, 60);
        assert!(parse_datetime("21/06/2026 12:00").is_err());
        assert!(parse_datetime("2026-06-21T99:00").is_err());
        assert!(parse_datetime("2026-13-01T12:00").is_err());
    }

    #[test]
    fn rejects_bad_parameters() {
        let path = layer_of(vec![box_mesh([0.0; 3], [1.0, 1.0, 1.0])]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SunShadowVolumeTool.validate(&args).is_err()
        };
        assert!(bad(
            json!({"datetime": "2026-06-21T12:00", "latitude": 45.0})
        ));
        assert!(bad(json!({"input": path, "latitude": 45.0})));
        assert!(bad(json!({"input": path, "datetime": "2026-06-21T12:00"})));
        assert!(bad(
            json!({"input": path, "datetime": "2026-06-21T12:00", "latitude": 200.0})
        ));
    }
}
