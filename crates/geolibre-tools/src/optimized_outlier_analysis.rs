//! GeoLibre tool: optimized outlier analysis (Anselin Local Moran's I).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Optimized Outlier Analysis* (Spatial
//! Statistics). It is the Local Moran companion to `optimized_hot_spot_analysis`:
//! aggregate raw incidents, auto-select the analysis distance band, run Anselin
//! Local Moran's I, and report False Discovery Rate–corrected cluster/outlier
//! types (HH, LL, HL, LH).
//!
//! The registry exposes the raw `local_morans_i_lisa` (whitebox) but no optimized
//! wrapper, and its two existing outlier tools solve different problems:
//! `local_outlier_analysis` (#164) is the *space-time* Local Moran on a
//! space-time cube, and `spatial_outlier_detection` (#159) is a Local Outlier
//! Factor (density) score — neither is the non-temporal Anselin Local Moran with
//! automatic aggregation, band selection, and FDR correction.
//!
//! With an `analysis_field` the input features are analyzed directly; without one,
//! incidents are aggregated (fishnet/snap). Output is a point layer carrying
//! `VALUE`, `LMiIndex`, `LMiZScore`, `LMiPValue`, and `COType` (HH/LL/HL/LH or
//! empty when not FDR-significant). Shares [`crate::hotspot_common`].

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, Geometry, GeometryType, Layer};

use crate::hotspot_common::{
    aggregate, auto_distance_band, bh_significant, default_cell_size, extract_points, local_moran,
    parse_aggregation,
};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct OptimizedOutlierAnalysisTool;

impl Tool for OptimizedOutlierAnalysisTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "optimized_outlier_analysis",
            display_name: "Optimized Outlier Analysis",
            summary: "Auto-aggregate incidents, auto-select the analysis distance band, run Anselin Local Moran's I, and report FDR-corrected cluster/outlier types (HH/LL/HL/LH) — like ArcGIS Optimized Outlier Analysis; the non-temporal Local Moran wrapper the raw bundled local_morans_i_lisa and the space-time local_outlier_analysis don't provide.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input feature layer (incident points, or features with an analysis field).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output point layer with Local Moran results. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "analysis_field",
                    description: "Numeric field to analyze (weighted, no aggregation). Omit to count incidents.",
                    required: false,
                },
                ToolParamSpec {
                    name: "aggregation",
                    description: "Incident aggregation when no analysis_field: 'fishnet' (grid count; default) or 'snap' (merge coincident).",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_size",
                    description: "Aggregation cell/snap size in CRS units. Default: extent's longer side / 30.",
                    required: false,
                },
                ToolParamSpec {
                    name: "distance_band",
                    description: "Analysis distance band in CRS units. Default: auto-selected by peak clustering intensity.",
                    required: false,
                },
                ToolParamSpec {
                    name: "fdr_alpha",
                    description: "False Discovery Rate level for COType significance (default 0.05).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_aggregation(args)?;
        opt_pos(args, "cell_size")?;
        opt_pos(args, "distance_band")?;
        parse_alpha(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let field = parse_optional_str(args, "analysis_field")?.map(String::from);
        let mode = parse_aggregation(args)?;
        let cell_arg = opt_pos(args, "cell_size")?;
        let band_arg = opt_pos(args, "distance_band")?;
        let alpha = parse_alpha(args)?;

        let layer = load_input_layer(input)?;
        let epsg = layer.crs_epsg();
        let (raw, skipped) = extract_points(&layer, field.as_deref())?;
        if raw.len() < 4 {
            return Err(ToolError::Execution(format!(
                "need at least 4 usable features, got {}",
                raw.len()
            )));
        }

        let cell = cell_arg.unwrap_or_else(|| default_cell_size(&raw));
        let pts = if field.is_some() {
            raw
        } else {
            aggregate(&raw, mode, cell)
        };
        if pts.len() < 4 {
            return Err(ToolError::Execution(
                "aggregation produced fewer than 4 features; try a smaller cell_size".to_string(),
            ));
        }

        let band = band_arg.unwrap_or_else(|| auto_distance_band(&pts));
        ctx.progress.info(&format!(
            "{} feature(s) after aggregation, distance band {band:.4}",
            pts.len()
        ));

        let lm = local_moran(&pts, band);
        let ps: Vec<f64> = lm.iter().map(|r| r.p).collect();
        let sig = bh_significant(&ps, alpha);

        // COType from the sign of z (cluster vs outlier) and the feature's value.
        let mut counts: BTreeMap<&'static str, usize> = BTreeMap::new();
        let cotypes: Vec<&'static str> = lm
            .iter()
            .enumerate()
            .map(|(k, r)| {
                if !sig[k] || r.z == 0.0 {
                    return "";
                }
                let t = if r.z > 0.0 {
                    // Positive local autocorrelation -> cluster (HH or LL).
                    if r.high {
                        "HH"
                    } else {
                        "LL"
                    }
                } else {
                    // Negative -> spatial outlier (HL or LH).
                    if r.high {
                        "HL"
                    } else {
                        "LH"
                    }
                };
                *counts.entry(t).or_insert(0) += 1;
                t
            })
            .collect();

        let mut out = Layer::new("optimized_outliers").with_geom_type(GeometryType::Point);
        if let Some(e) = epsg {
            out = out.with_crs_epsg(e);
        }
        out.add_field(FieldDef::new("VALUE", FieldType::Float));
        out.add_field(FieldDef::new("LMiIndex", FieldType::Float));
        out.add_field(FieldDef::new("LMiZScore", FieldType::Float));
        out.add_field(FieldDef::new("LMiPValue", FieldType::Float));
        out.add_field(FieldDef::new("COType", FieldType::Text));
        for (k, p) in pts.iter().enumerate() {
            out.add_feature(
                Some(Geometry::point(p.x, p.y)),
                &[
                    ("VALUE", p.val.into()),
                    ("LMiIndex", lm[k].i.into()),
                    ("LMiZScore", lm[k].z.into()),
                    ("LMiPValue", lm[k].p.into()),
                    ("COType", cotypes[k].to_string().into()),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed adding feature: {e}")))?;
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(pts.len()));
        outputs.insert("skipped".to_string(), json!(skipped));
        outputs.insert("distance_band".to_string(), json!(band));
        outputs.insert("cell_size".to_string(), json!(cell));
        outputs.insert("aggregated".to_string(), json!(field.is_none()));
        for t in ["HH", "LL", "HL", "LH"] {
            outputs.insert(t.to_string(), json!(counts.get(t).copied().unwrap_or(0)));
        }
        Ok(ToolRunResult { outputs })
    }
}

fn parse_alpha(args: &ToolArgs) -> Result<f64, ToolError> {
    match args.get("fdr_alpha") {
        None | Some(Value::Null) => Ok(0.05),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(0.05),
        v => {
            let x = match v {
                Some(Value::Number(n)) => n.as_f64(),
                Some(Value::String(s)) => s.trim().parse::<f64>().ok(),
                _ => None,
            };
            match x {
                Some(a) if a > 0.0 && a < 1.0 => Ok(a),
                _ => Err(ToolError::Validation(
                    "'fdr_alpha' must be in (0, 1)".to_string(),
                )),
            }
        }
    }
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

fn opt_pos(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::Number(n)) => check_pos(n.as_f64(), key),
        Some(Value::String(s)) => check_pos(s.trim().parse::<f64>().ok(), key),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a number"
        ))),
    }
}

fn check_pos(v: Option<f64>, key: &str) -> Result<Option<f64>, ToolError> {
    match v {
        Some(x) if x.is_finite() && x > 0.0 => Ok(Some(x)),
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a positive number"
        ))),
        None => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a number"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Coord, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// A point layer with a numeric analysis field.
    fn valued_layer(pts: &[(f64, f64, f64)]) -> String {
        let mut l = Layer::new("v")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("val", FieldType::Float));
        for (x, y, v) in pts {
            l.add_feature(
                Some(Geometry::Point(Coord::xy(*x, *y))),
                &[("val", (*v).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = OptimizedOutlierAnalysisTool.run(&args, &ctx()).unwrap();
        let l = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, l)
    }

    /// A single high value surrounded by low values is a High-Low spatial
    /// outlier; the output carries the five result fields.
    #[test]
    fn detects_high_low_outlier() {
        // 5x5 grid, all value 1 except a spike of 100 at the center.
        let mut pts = Vec::new();
        for r in 0..5 {
            for c in 0..5 {
                let v = if r == 2 && c == 2 { 100.0 } else { 1.0 };
                pts.push((c as f64, r as f64, v));
            }
        }
        let (out, l) = run(json!({
            "input": valued_layer(&pts), "analysis_field": "val", "distance_band": 1.5,
        }));
        for f in ["VALUE", "LMiIndex", "LMiZScore", "LMiPValue", "COType"] {
            assert!(l.schema.field_index(f).is_some(), "missing field {f}");
        }
        // With a strong spike the tool should flag at least one HL outlier.
        assert!(
            out.outputs["HL"].as_i64().unwrap() >= 1,
            "expected a High-Low outlier, got {:?}",
            out.outputs
        );
    }

    /// Two adjacent high-value clusters produce High-High clustering.
    #[test]
    fn detects_high_high_cluster() {
        let mut pts = Vec::new();
        for r in 0..6 {
            for c in 0..6 {
                // Left half high, right half low.
                let v = if c < 3 { 50.0 } else { 1.0 };
                pts.push((c as f64, r as f64, v));
            }
        }
        let (out, _l) = run(json!({
            "input": valued_layer(&pts), "analysis_field": "val", "distance_band": 1.5,
        }));
        assert!(
            out.outputs["HH"].as_i64().unwrap() >= 1 || out.outputs["LL"].as_i64().unwrap() >= 1,
            "expected clustering, got {:?}",
            out.outputs
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            OptimizedOutlierAnalysisTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "fdr_alpha": 1.5 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "aggregation": "grid" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "fdr_alpha": 0.1 })).is_ok());
    }
}
