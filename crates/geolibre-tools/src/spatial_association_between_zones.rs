//! GeoLibre tool: spatial association between two categorical zone rasters.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Spatial Association Between Zones*
//! (Spatial Statistics): measures how strongly two categorical maps of the same
//! area correspond — e.g. land cover vs. ecoregions, or a predicted vs. an
//! observed classification. The bundled `cross_tabulation` produces only the
//! raw contingency table between two rasters; it computes no association
//! strength and no per-zone diagnostics. This tool turns that table into the
//! information-theoretic **V-measure** (Rosenberg & Hirschberg 2007) plus
//! per-zone association scores.
//!
//! For every cell where both rasters are valid, the joint class counts form a
//! contingency table, from which:
//!
//! - **homogeneity** `h = 1 − H(Z1 | Z2) / H(Z1)` — how completely each zone-2
//!   class falls inside a single zone-1 class (0 = independent, 1 = every
//!   zone-2 patch is pure in zone 1),
//! - **completeness** `c = 1 − H(Z2 | Z1) / H(Z2)` — the symmetric direction, and
//! - **V-measure** `V = 2hc / (h + c)` — their harmonic mean, the overall
//!   association (0 = independent partitions, 1 = identical up to relabelling).
//!
//! Per zone-1 class `k`, an association score `1 − H(Z2 | Z1=k) / H(Z2)` says
//! how tightly that zone concentrates within the other map, alongside its
//! dominant zone-2 class and that class's share. Entropies use natural logs, so
//! the ratios (and hence h, c, V) are base-independent and match scikit-learn's
//! `homogeneity_completeness_v_measure` exactly.
//!
//! v1 takes two categorical **rasters** on the same grid (like the bundled
//! `cross_tabulation`); overlaying two polygon zonings is future work.

use std::collections::{BTreeMap, HashMap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};

use crate::common::{band_to_vec, load_input_raster, parse_optional_output, write_text_output};

pub struct SpatialAssociationBetweenZonesTool;

impl Tool for SpatialAssociationBetweenZonesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "spatial_association_between_zones",
            display_name: "Spatial Association Between Zones",
            summary: "Measure how strongly two categorical zone rasters correspond via the V-measure (homogeneity/completeness) plus per-zone association scores.",
            category: ToolCategory::Raster,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "zones1",
                    description: "First categorical zone raster (e.g. land cover). Same grid as zones2.",
                    required: true,
                },
                ToolParamSpec {
                    name: "zones2",
                    description: "Second categorical zone raster (e.g. ecoregions) on the same grid as zones1.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional per-zone-1-class CSV table (class, count, association, dominant zone-2 class and share). Overall metrics are always returned in the result.",
                    required: false,
                },
                ToolParamSpec {
                    name: "band1",
                    description: "1-based band of zones1 to read (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band2",
                    description: "1-based band of zones2 to read (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        for key in ["zones1", "zones2"] {
            if args
                .get(key)
                .and_then(Value::as_str)
                .map(str::trim)
                .unwrap_or("")
                .is_empty()
            {
                return Err(ToolError::Validation(format!(
                    "missing required string parameter '{key}'"
                )));
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let z1_path = req_str(args, "zones1")?;
        let z2_path = req_str(args, "zones2")?;
        let output = parse_optional_output(args, "output")?;
        let band1 = band_index(args, "band1")?;
        let band2 = band_index(args, "band2")?;

        let r1 = load_input_raster(z1_path)?;
        let r2 = load_input_raster(z2_path)?;
        if r1.rows != r2.rows || r1.cols != r2.cols {
            return Err(ToolError::Validation(format!(
                "zones1 ({}x{}) and zones2 ({}x{}) must share the same grid",
                r1.rows, r1.cols, r2.rows, r2.cols
            )));
        }
        if band1 as usize >= r1.bands {
            return Err(ToolError::Validation(format!(
                "band1 out of range (zones1 has {} band(s))",
                r1.bands
            )));
        }
        if band2 as usize >= r2.bands {
            return Err(ToolError::Validation(format!(
                "band2 out of range (zones2 has {} band(s))",
                r2.bands
            )));
        }

        let a = band_to_vec(&r1, band1);
        let b = band_to_vec(&r2, band2);
        let (nd1, nd2) = (r1.nodata, r2.nodata);

        ctx.progress.info("building contingency table");
        // Joint and marginal counts over cells valid in both rasters.
        let mut joint: HashMap<(i64, i64), u64> = HashMap::new();
        let mut m1: BTreeMap<i64, u64> = BTreeMap::new();
        let mut m2: BTreeMap<i64, u64> = BTreeMap::new();
        let mut n: u64 = 0;
        for (&v1, &v2) in a.iter().zip(b.iter()) {
            if v1 == nd1 || v2 == nd2 || !v1.is_finite() || !v2.is_finite() {
                continue;
            }
            let c1 = v1.round() as i64;
            let c2 = v2.round() as i64;
            *joint.entry((c1, c2)).or_insert(0) += 1;
            *m1.entry(c1).or_insert(0) += 1;
            *m2.entry(c2).or_insert(0) += 1;
            n += 1;
        }
        if n == 0 {
            return Err(ToolError::Execution(
                "no cells are valid in both rasters".to_string(),
            ));
        }

        let nf = n as f64;
        let h1 = entropy(m1.values().map(|&c| c as f64 / nf));
        let h2 = entropy(m2.values().map(|&c| c as f64 / nf));
        // Conditional entropies from the joint table.
        // H(Z1|Z2) = -sum p(c1,c2) log( p(c1,c2)/p(c2) ).
        let mut h1_given2 = 0.0;
        let mut h2_given1 = 0.0;
        let mut mutual = 0.0;
        for (&(c1, c2), &cnt) in &joint {
            let pij = cnt as f64 / nf;
            let p1 = m1[&c1] as f64 / nf;
            let p2 = m2[&c2] as f64 / nf;
            h1_given2 += -pij * (pij / p2).ln();
            h2_given1 += -pij * (pij / p1).ln();
            mutual += pij * (pij / (p1 * p2)).ln();
        }

        // Rosenberg–Hirschberg conventions: perfectly homogeneous/complete when
        // the conditioning entropy is zero (including the degenerate single-class
        // cases), matching scikit-learn.
        let homogeneity = if h1 == 0.0 { 1.0 } else { 1.0 - h1_given2 / h1 };
        let completeness = if h2 == 0.0 { 1.0 } else { 1.0 - h2_given1 / h2 };
        let v_measure = if homogeneity + completeness == 0.0 {
            0.0
        } else {
            2.0 * homogeneity * completeness / (homogeneity + completeness)
        };

        ctx.progress.info(&format!(
            "V-measure {v_measure:.4} (homogeneity {homogeneity:.4}, completeness {completeness:.4}) over {n} cells"
        ));

        // Per-zone-1 association: 1 - H(Z2 | Z1=k)/H(Z2), plus dominant zone-2.
        let mut zone_rows: Vec<ZoneRow> = Vec::new();
        for (&c1, &m1c) in &m1 {
            let total = m1c as f64;
            let mut cond = 0.0; // H(Z2 | Z1=c1)
            let mut dominant = (0i64, 0u64);
            for (&(a1, a2), &cnt) in &joint {
                if a1 != c1 {
                    continue;
                }
                let p = cnt as f64 / total;
                if p > 0.0 {
                    cond += -p * p.ln();
                }
                if cnt > dominant.1 {
                    dominant = (a2, cnt);
                }
            }
            let assoc = if h2 == 0.0 { 1.0 } else { 1.0 - cond / h2 };
            zone_rows.push(ZoneRow {
                class: c1,
                count: m1c,
                association: assoc,
                dominant_class: dominant.0,
                dominant_share: dominant.1 as f64 / total,
            });
        }

        let mut outputs = BTreeMap::new();
        if let Some(path) = output {
            let mut csv =
                String::from("zone1_class,count,association,dominant_zone2_class,dominant_share\n");
            for z in &zone_rows {
                csv.push_str(&format!(
                    "{},{},{:.6},{},{:.6}\n",
                    z.class, z.count, z.association, z.dominant_class, z.dominant_share
                ));
            }
            write_text_output(&csv, path)?;
            outputs.insert("output".to_string(), json!(path));
        }

        outputs.insert("n_cells".to_string(), json!(n));
        outputs.insert("num_zones1".to_string(), json!(m1.len()));
        outputs.insert("num_zones2".to_string(), json!(m2.len()));
        outputs.insert("homogeneity".to_string(), json!(homogeneity));
        outputs.insert("completeness".to_string(), json!(completeness));
        outputs.insert("v_measure".to_string(), json!(v_measure));
        outputs.insert("mutual_information".to_string(), json!(mutual));
        Ok(ToolRunResult { outputs })
    }
}

struct ZoneRow {
    class: i64,
    count: u64,
    association: f64,
    dominant_class: i64,
    dominant_share: f64,
}

/// Shannon entropy (natural log) of a probability distribution.
fn entropy(probs: impl Iterator<Item = f64>) -> f64 {
    let mut h = 0.0;
    for p in probs {
        if p > 0.0 {
            h += -p * p.ln();
        }
    }
    h
}

fn req_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required parameter '{key}'")))
}

/// Parses a 1-based band parameter into a 0-based index (default band 1 -> 0).
fn band_index(args: &ToolArgs, key: &str) -> Result<isize, ToolError> {
    let one_based = match args.get(key) {
        None | Some(Value::Null) => 1,
        Some(Value::Number(nu)) => nu.as_u64().unwrap_or(1),
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
            .map_err(|_| ToolError::Validation(format!("parameter '{key}' must be an integer")))?,
        Some(_) => {
            return Err(ToolError::Validation(format!(
                "parameter '{key}' must be an integer"
            )))
        }
    };
    Ok(one_based.max(1) as isize - 1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{memory_store, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn raster_from(rows: usize, cols: usize, data: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: Some(1.0),
            nodata: -9999.0,
            data_type: wbraster::DataType::F32,
            crs: Default::default(),
            metadata: vec![],
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, data[row * cols + col])
                    .unwrap();
            }
        }
        let id = memory_store::put_raster(r);
        memory_store::make_raster_memory_path(&id)
    }

    fn run(z1: &str, z2: &str) -> ToolRunResult {
        let args: ToolArgs = serde_json::from_value(json!({ "zones1": z1, "zones2": z2 })).unwrap();
        SpatialAssociationBetweenZonesTool
            .run(&args, &ctx())
            .unwrap()
    }

    /// Identical zonings -> perfect association (h = c = V = 1).
    #[test]
    fn identical_zonings_are_perfectly_associated() {
        let z = raster_from(2, 2, &[1.0, 1.0, 2.0, 2.0]);
        let out = run(&z, &z);
        assert!((out.outputs["v_measure"].as_f64().unwrap() - 1.0).abs() < 1e-9);
        assert!((out.outputs["homogeneity"].as_f64().unwrap() - 1.0).abs() < 1e-9);
        assert!((out.outputs["completeness"].as_f64().unwrap() - 1.0).abs() < 1e-9);
    }

    /// A finer zoning nested inside a coarser one: completeness 1 (fine
    /// determines coarse) but homogeneity < 1 (coarse does not determine fine).
    #[test]
    fn nested_zoning_is_complete_but_not_homogeneous() {
        // coarse: two halves; fine: four quarters (each nested in a half).
        let coarse = raster_from(2, 2, &[1.0, 1.0, 2.0, 2.0]);
        let fine = raster_from(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        // zones1 = fine (classes), zones2 = coarse (partition).
        let out = run(&fine, &coarse);
        let h = out.outputs["homogeneity"].as_f64().unwrap();
        let c = out.outputs["completeness"].as_f64().unwrap();
        // Each coarse class is split across two fine classes -> homogeneity < 1.
        assert!(h < 0.99, "homogeneity {h} should be < 1");
        // Each fine class lies entirely in one coarse class -> completeness = 1.
        assert!((c - 1.0).abs() < 1e-9, "completeness {c} should be 1");
    }

    /// Independent (checkerboard vs. stripes) -> low association.
    #[test]
    fn independent_zonings_have_low_v() {
        // z1 stripes by row, z2 stripes by column: knowing one says nothing
        // about the other -> mutual information 0, V-measure 0.
        let z1 = raster_from(2, 2, &[1.0, 1.0, 2.0, 2.0]);
        let z2 = raster_from(2, 2, &[1.0, 2.0, 1.0, 2.0]);
        let out = run(&z1, &z2);
        assert!(out.outputs["v_measure"].as_f64().unwrap().abs() < 1e-9);
        assert!(out.outputs["mutual_information"].as_f64().unwrap().abs() < 1e-9);
    }

    #[test]
    fn rejects_mismatched_grids() {
        let z1 = raster_from(2, 2, &[1.0, 1.0, 2.0, 2.0]);
        let z2 = raster_from(2, 3, &[1.0, 1.0, 2.0, 2.0, 1.0, 2.0]);
        let args: ToolArgs = serde_json::from_value(json!({ "zones1": z1, "zones2": z2 })).unwrap();
        assert!(SpatialAssociationBetweenZonesTool
            .run(&args, &ctx())
            .is_err());
    }

    #[test]
    fn rejects_missing_input() {
        let tool = SpatialAssociationBetweenZonesTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "zones1": "a.tif" })).is_err());
        assert!(bad(json!({ "zones1": "a.tif", "zones2": "b.tif" })).is_ok());
    }
}
