//! GeoLibre tool: polarimetric SAR vegetation and structure indices.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Compute SAR Indices* (Image Analyst).
//!
//! ## Why `spectral_index` does not cover this
//!
//! `spectral_index` is optical band math: its formulas, band conventions and
//! naming are all reflectance-domain. Polarimetric indices are ratios between
//! *polarization channels* (VV, VH, HH, HV) of a radar acquisition, which
//! `spectral_index` has no concept of.
//!
//! These indices are the cloud-independent counterpart to NDVI, and they are
//! the analysis endpoint that makes the rest of the SAR chain worth having:
//! `multilook` → `apply_radiometric_calibration` → a speckle filter → an index.
//! Without them the SAR tools produce inputs with nothing to consume them.
//!
//! ## The dB trap
//!
//! Every index here is a **ratio**, and ratios of logarithms are not logarithms
//! of ratios — computing RVI on dB input yields a number with no physical
//! meaning that still looks plausible on a map. Because SAR products commonly
//! ship in dB, `input_units` is explicit and dB input is converted to linear
//! before any ratio arithmetic rather than being silently accepted.

use std::collections::BTreeMap;

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::DataType;

use crate::args_common::{choice_or, opt_choice, req_str};
use crate::common::{
    load_input_raster, parse_optional_output, raster_like_with_data, write_or_store_output,
};

pub struct ComputeSarIndicesTool;

impl Tool for ComputeSarIndicesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "compute_sar_indices",
            display_name: "Compute SAR Indices",
            summary: "Derives polarimetric SAR indices — RVI (Radar Vegetation Index), RFDI (Radar Forest Degradation Index), CSI (Canopy Structure Index) and DPSVI (Dual-Pol SAR Vegetation Index) — from multi-polarization backscatter (ArcGIS Compute SAR Indices). spectral_index covers optical band math only: its formulas and band conventions are reflectance-domain and it has no notion of polarization channels, so the catalog's SAR chain (multilook, sar_coherence, speckle filters) currently has no analysis endpoint.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Multi-polarization calibrated SAR raster (one band per polarization).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output index raster. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "index",
                    description: "'rvi' (default), 'rfdi', 'csi', or 'dpsvi'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "polarization_bands",
                    description: "1-based band mapping, e.g. 'vv=1,vh=2' or 'hh=1,hv=2,vv=3'. Default by band count: 2 bands = 'vv=1,vh=2', 3+ bands = 'hh=1,hv=2,vv=3'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "input_units",
                    description: "'linear' (default) or 'db'. dB input is converted to linear power before the ratio arithmetic, because ratios of decibels are not meaningful.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        parse_index(args)?;
        choice_or(args, "input_units", &["linear", "db"], "linear")?;
        if let Some(spec) = opt_choice(args, "polarization_bands") {
            parse_band_map(&spec)?;
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?.to_string();
        let index = parse_index(args)?;
        let db_input = choice_or(args, "input_units", &["linear", "db"], "linear")? == "db";
        let output = parse_optional_output(args, "output")?;

        let raster = load_input_raster(&input)?;
        let (rows, cols) = (raster.rows, raster.cols);

        let map = match opt_choice(args, "polarization_bands") {
            Some(spec) => parse_band_map(&spec)?,
            None => default_band_map(raster.bands)?,
        };
        // Fail up front on a missing channel rather than producing an
        // all-no-data raster that looks like a data problem.
        let satisfied = index
            .channel_forms()
            .iter()
            .find(|form| form.iter().all(|pol| map.contains_key(*pol)))
            .ok_or_else(|| {
                ToolError::Validation(format!(
                    "index '{}' needs {}; supply the channels via 'polarization_bands' (e.g. '{}')",
                    index.label(),
                    index
                        .channel_forms()
                        .iter()
                        .map(|f| f.join("+"))
                        .collect::<Vec<_>>()
                        .join(" or "),
                    index.example_mapping()
                ))
            })?;
        for pol in satisfied.iter() {
            let band = map[*pol];
            if band >= raster.bands {
                return Err(ToolError::Validation(format!(
                    "polarization {pol} maps to band {} but '{input}' has {} band(s)",
                    band + 1,
                    raster.bands
                )));
            }
        }

        ctx.progress.info(&format!(
            "{rows}x{cols}, index {}, {} input",
            index.label(),
            if db_input { "dB" } else { "linear" }
        ));

        let nodata = -9999.0_f64;
        let mut out = vec![nodata; rows * cols];
        let fetch = |pol: &str, r: usize, c: usize| -> Option<f64> {
            let band = *map.get(pol)? as isize;
            let v = raster.get(band, r as isize, c as isize);
            if v == raster.nodata || !v.is_finite() {
                return None;
            }
            // Convert first: every formula below is a ratio of linear powers.
            Some(if db_input { 10.0_f64.powf(v / 10.0) } else { v })
        };

        for r in 0..rows {
            for c in 0..cols {
                if let Some(v) = index.compute(&|pol| fetch(pol, r, c)) {
                    out[r * cols + c] = v;
                }
            }
            ctx.progress.progress((r as f64 + 1.0) / rows as f64);
        }

        let out_r = raster_like_with_data(&raster, out, nodata, DataType::F32)?;
        let out_path = write_or_store_output(out_r, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("index".to_string(), json!(index.label()));
        outputs.insert(
            "polarization_bands".to_string(),
            json!(describe_map(&map)),
        );
        outputs.insert("rows".to_string(), json!(rows));
        outputs.insert("cols".to_string(), json!(cols));
        Ok(ToolRunResult { outputs })
    }
}

// ── Indices ─────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Index {
    Rvi,
    Rfdi,
    Csi,
    Dpsvi,
}

impl Index {
    fn label(self) -> &'static str {
        match self {
            Index::Rvi => "rvi",
            Index::Rfdi => "rfdi",
            Index::Csi => "csi",
            Index::Dpsvi => "dpsvi",
        }
    }

    /// The acceptable channel sets, most specific first.
    ///
    /// RVI has two published forms — quad-pol (HH/VV/HV) and dual-pol (VV/VH) —
    /// and either is valid input, so this returns alternatives rather than a
    /// single required set. Checking only the dual-pol minimum would reject
    /// perfectly good quad-pol data, which has no VH channel at all.
    fn channel_forms(self) -> &'static [&'static [&'static str]] {
        match self {
            Index::Rvi => &[&["hh", "vv", "hv"], &["vv", "vh"]],
            Index::Rfdi => &[&["hh", "hv"]],
            Index::Csi => &[&["vv", "hh"]],
            Index::Dpsvi => &[&["vv", "vh"]],
        }
    }

    fn example_mapping(self) -> &'static str {
        match self {
            Index::Rvi | Index::Dpsvi => "vv=1,vh=2",
            Index::Rfdi => "hh=1,hv=2",
            Index::Csi => "hh=1,vv=2",
        }
    }

    /// Evaluates the index from a channel accessor. `None` when any required
    /// channel is no-data or the denominator vanishes.
    fn compute(self, get: &dyn Fn(&str) -> Option<f64>) -> Option<f64> {
        let ratio = |num: f64, den: f64| (den.abs() > f64::EPSILON).then_some(num / den);
        match self {
            Index::Rvi => {
                // Quad-pol form when the full basis is available, dual-pol
                // otherwise. Both are the standard published definitions.
                match (get("hh"), get("vv"), get("hv")) {
                    (Some(hh), Some(vv), Some(hv)) => ratio(8.0 * hv, hh + vv + 2.0 * hv),
                    _ => {
                        let (vv, vh) = (get("vv")?, get("vh")?);
                        ratio(4.0 * vh, vv + vh)
                    }
                }
            }
            Index::Rfdi => {
                let (hh, hv) = (get("hh")?, get("hv")?);
                ratio(hh - hv, hh + hv)
            }
            Index::Csi => {
                let (vv, hh) = (get("vv")?, get("hh")?);
                ratio(vv, vv + hh)
            }
            Index::Dpsvi => {
                let (vv, vh) = (get("vv")?, get("vh")?);
                ratio(vh * (vv + vh), vv)
            }
        }
    }
}

fn parse_index(args: &ToolArgs) -> Result<Index, ToolError> {
    Ok(
        match choice_or(args, "index", &["rvi", "rfdi", "csi", "dpsvi"], "rvi")? {
            "rvi" => Index::Rvi,
            "rfdi" => Index::Rfdi,
            "csi" => Index::Csi,
            _ => Index::Dpsvi,
        },
    )
}

// ── Polarization band mapping ───────────────────────────────────────────────

type BandMap = BTreeMap<String, usize>;

/// Parses `'vv=1,vh=2'` into 0-based band indices.
fn parse_band_map(spec: &str) -> Result<BandMap, ToolError> {
    let mut map = BandMap::new();
    for part in spec.split(',').map(str::trim).filter(|s| !s.is_empty()) {
        let (pol, band) = part.split_once('=').ok_or_else(|| {
            ToolError::Validation(format!(
                "'polarization_bands' entry '{part}' must look like 'vv=1'"
            ))
        })?;
        let pol = pol.trim().to_ascii_lowercase();
        if !["hh", "hv", "vh", "vv"].contains(&pol.as_str()) {
            return Err(ToolError::Validation(format!(
                "unknown polarization '{pol}'; expected hh, hv, vh or vv"
            )));
        }
        let band: usize = band.trim().parse().map_err(|_| {
            ToolError::Validation(format!("band for '{pol}' must be a whole number"))
        })?;
        if band == 0 {
            return Err(ToolError::Validation(format!(
                "'polarization_bands' is 1-based; use 1 for the first band (got {pol}=0)"
            )));
        }
        map.insert(pol, band - 1);
    }
    if map.is_empty() {
        return Err(ToolError::Validation(
            "'polarization_bands' is empty".to_string(),
        ));
    }
    Ok(map)
}

fn default_band_map(bands: usize) -> Result<BandMap, ToolError> {
    let mut map = BandMap::new();
    match bands {
        0 | 1 => {
            return Err(ToolError::Validation(format!(
                "SAR indices need at least two polarization bands, input has {bands}"
            )))
        }
        2 => {
            // The near-universal Sentinel-1 dual-pol ordering.
            map.insert("vv".into(), 0);
            map.insert("vh".into(), 1);
        }
        _ => {
            map.insert("hh".into(), 0);
            map.insert("hv".into(), 1);
            map.insert("vv".into(), 2);
        }
    }
    Ok(map)
}

fn describe_map(map: &BandMap) -> String {
    map.iter()
        .map(|(pol, band)| format!("{pol}={}", band + 1))
        .collect::<Vec<_>>()
        .join(",")
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Output rasters are F32, so value comparisons use a **relative**
    /// tolerance — 0.8 and 0.4 are not exactly representable in f32.
    fn close(actual: f64, expect: f64) -> bool {
        (actual - expect).abs() <= 1e-6 * expect.abs().max(1.0)
    }

    /// A 1x1 raster with one band per supplied value.
    fn cell(bands: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols: 1,
            rows: 1,
            bands: bands.len(),
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F32,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for (b, v) in bands.iter().enumerate() {
            r.set(b as isize, 0, 0, *v).unwrap();
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: Value) -> Raster {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = ComputeSarIndicesTool.run(&args, &ctx()).unwrap();
        load_input_raster(res.outputs["output"].as_str().unwrap()).unwrap()
    }

    #[test]
    fn dual_pol_rvi_matches_the_published_formula() {
        // VV=0.2, VH=0.05 gives 4*0.05 / (0.2+0.05) = 0.8
        let out = run(json!({"input": cell(&[0.2, 0.05]), "index": "rvi"}));
        assert!(close(out.get(0, 0, 0), 0.8));
    }

    #[test]
    fn quad_pol_rvi_is_used_when_the_full_basis_is_present() {
        // HH=0.3, HV=0.1, VV=0.2 gives 8*0.1 / (0.3+0.2+0.2) = 0.8/0.7.
        // Quad-pol data has no VH channel at all, so requiring VH would
        // wrongly reject this input.
        let out = run(json!({
            "input": cell(&[0.3, 0.1, 0.2]),
            "index": "rvi",
            "polarization_bands": "hh=1,hv=2,vv=3",
        }));
        assert!(close(out.get(0, 0, 0), 0.8 / 0.7));
    }

    #[test]
    fn rfdi_and_csi_match_their_formulas() {
        // RFDI = (HH - HV)/(HH + HV) = (0.3-0.1)/0.4 = 0.5
        let out = run(json!({
            "input": cell(&[0.3, 0.1]), "index": "rfdi", "polarization_bands": "hh=1,hv=2",
        }));
        assert!(close(out.get(0, 0, 0), 0.5));
        // CSI = VV/(VV + HH) = 0.2/(0.2+0.3) = 0.4
        let out = run(json!({
            "input": cell(&[0.3, 0.2]), "index": "csi", "polarization_bands": "hh=1,vv=2",
        }));
        assert!(close(out.get(0, 0, 0), 0.4));
    }

    #[test]
    fn dpsvi_matches_its_formula() {
        // DPSVI = VH*(VV + VH)/VV = 0.05*0.25/0.2 = 0.0625
        let out = run(json!({"input": cell(&[0.2, 0.05]), "index": "dpsvi"}));
        assert!(close(out.get(0, 0, 0), 0.0625));
    }

    #[test]
    fn db_input_is_converted_before_the_ratio_not_after() {
        // The trap the tool exists to prevent. The same physical values
        // expressed in dB must give the SAME index.
        let linear = run(json!({"input": cell(&[0.2, 0.05]), "index": "rvi"}));
        let db_vv = 10.0 * 0.2_f64.log10();
        let db_vh = 10.0 * 0.05_f64.log10();
        let from_db = run(json!({
            "input": cell(&[db_vv, db_vh]), "index": "rvi", "input_units": "db",
        }));
        assert!(
            close(from_db.get(0, 0, 0), linear.get(0, 0, 0)),
            "dB path gave {} but linear gave {}",
            from_db.get(0, 0, 0),
            linear.get(0, 0, 0)
        );
        // Treating dB as linear gives something quite different, so the check
        // above is not vacuous.
        let mistaken = run(json!({"input": cell(&[db_vv, db_vh]), "index": "rvi"}));
        assert!((mistaken.get(0, 0, 0) - linear.get(0, 0, 0)).abs() > 0.1);
    }

    #[test]
    fn nodata_in_any_required_channel_yields_nodata() {
        let out = run(json!({"input": cell(&[-9999.0, 0.05]), "index": "rvi"}));
        assert_eq!(out.get(0, 0, 0), out.nodata);
    }

    #[test]
    fn a_vanishing_denominator_yields_nodata_not_infinity() {
        // VV = -VH makes the RVI denominator zero.
        let out = run(json!({"input": cell(&[0.05, -0.05]), "index": "rvi"}));
        assert_eq!(out.get(0, 0, 0), out.nodata);
    }

    #[test]
    fn missing_required_channel_is_rejected_up_front() {
        // The default 2-band mapping is vv/vh, so RFDI (hh/hv) cannot be met.
        let args: ToolArgs =
            serde_json::from_value(json!({"input": cell(&[0.3, 0.1]), "index": "rfdi"})).unwrap();
        let err = ComputeSarIndicesTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err:?}").contains("hh"), "unhelpful error: {err:?}");
    }

    #[test]
    fn rejects_bad_parameters() {
        let src = cell(&[0.2, 0.05]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ComputeSarIndicesTool.validate(&args).is_err()
        };
        assert!(bad(json!({"input": src, "index": "nope"})));
        assert!(bad(json!({"input": src, "input_units": "watts"})));
        assert!(bad(json!({"input": src, "polarization_bands": "xx=1"})));
        assert!(bad(json!({"input": src, "polarization_bands": "vv=0"})));
        assert!(bad(json!({"input": src, "polarization_bands": "vv"})));
    }
}
