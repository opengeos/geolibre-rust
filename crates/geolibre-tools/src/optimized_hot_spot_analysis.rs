//! GeoLibre tool: optimized hot spot analysis (Getis-Ord Gi\*).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Optimized Hot Spot Analysis* (Spatial
//! Statistics). The registry exposes only the raw statistic (whitebox's
//! `getis_ord_gi_star`), which makes the user pre-aggregate incidents, hand-pick a
//! distance band, and read raw z-scores. This tool is the workflow wrapper: it
//! aggregates raw incident points (fishnet grid count or snap-to-weighted-points),
//! auto-selects the analysis distance band, runs Gi\*, and reports a False
//! Discovery Rate–corrected hot/cold-spot bin — the automation that makes hot spot
//! analysis reliable.
//!
//! With an `analysis_field` the input features are analyzed directly (weighted, no
//! aggregation); without one, each feature is a unit incident to aggregate. The
//! output is a point layer carrying `VALUE` (the analyzed value), `GiZScore`,
//! `GiPValue`, and `Gi_Bin` (−3…3; ±1/±2/±3 = 90/95/99% FDR-significant
//! cold/hot spots, 0 = not significant). Reuses the shared incident/statistics
//! machinery in [`crate::hotspot_common`], as does `optimized_outlier_analysis`.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, Geometry, GeometryType, Layer};

use crate::hotspot_common::{
    aggregate, auto_distance_band, default_cell_size, extract_points, fdr_bins, getis_gi_star,
    parse_aggregation,
};
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct OptimizedHotSpotAnalysisTool;

impl Tool for OptimizedHotSpotAnalysisTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "optimized_hot_spot_analysis",
            display_name: "Optimized Hot Spot Analysis",
            summary: "Auto-aggregate incidents, auto-select the analysis distance band, run Getis-Ord Gi*, and report an FDR-corrected hot/cold-spot bin (like ArcGIS Optimized Hot Spot Analysis) — the workflow wrapper the raw bundled getis_ord_gi_star lacks.",
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
                    description: "Output point layer with Gi* results. If omitted, stored in memory.",
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
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_aggregation(args)?;
        opt_pos(args, "cell_size")?;
        opt_pos(args, "distance_band")?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let field = parse_optional_str(args, "analysis_field")?.map(String::from);
        let mode = parse_aggregation(args)?;
        let cell_arg = opt_pos(args, "cell_size")?;
        let band_arg = opt_pos(args, "distance_band")?;

        let layer = load_input_layer(input)?;
        let epsg = layer.crs_epsg();
        let (raw, skipped) = extract_points(&layer, field.as_deref())?;
        if raw.len() < 3 {
            return Err(ToolError::Execution(format!(
                "need at least 3 usable features, got {}",
                raw.len()
            )));
        }

        // Aggregate incidents only when counting (no analysis field).
        let cell = cell_arg.unwrap_or_else(|| default_cell_size(&raw));
        let pts = if field.is_some() {
            raw
        } else {
            aggregate(&raw, mode, cell)
        };
        if pts.len() < 3 {
            return Err(ToolError::Execution(
                "aggregation produced fewer than 3 features; try a smaller cell_size".to_string(),
            ));
        }

        let band = band_arg.unwrap_or_else(|| auto_distance_band(&pts));
        ctx.progress.info(&format!(
            "{} feature(s) after aggregation, distance band {band:.4}",
            pts.len()
        ));

        let gi = getis_gi_star(&pts, band);
        let zs: Vec<f64> = gi.iter().map(|(z, _)| *z).collect();
        let ps: Vec<f64> = gi.iter().map(|(_, p)| *p).collect();
        let bins = fdr_bins(&ps, &zs);

        let hot = bins.iter().filter(|&&b| b > 0).count();
        let cold = bins.iter().filter(|&&b| b < 0).count();

        // Build the output point layer.
        let mut out = Layer::new("optimized_hot_spots").with_geom_type(GeometryType::Point);
        if let Some(e) = epsg {
            out = out.with_crs_epsg(e);
        }
        out.add_field(FieldDef::new("VALUE", FieldType::Float));
        out.add_field(FieldDef::new("GiZScore", FieldType::Float));
        out.add_field(FieldDef::new("GiPValue", FieldType::Float));
        out.add_field(FieldDef::new("Gi_Bin", FieldType::Integer));
        for (k, p) in pts.iter().enumerate() {
            out.add_feature(
                Some(Geometry::point(p.x, p.y)),
                &[
                    ("VALUE", p.val.into()),
                    ("GiZScore", zs[k].into()),
                    ("GiPValue", ps[k].into()),
                    ("Gi_Bin", (bins[k] as i64).into()),
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
        outputs.insert("hot_spots".to_string(), json!(hot));
        outputs.insert("cold_spots".to_string(), json!(cold));
        Ok(ToolRunResult { outputs })
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
    use wbvector::{memory_store, Coord, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// A point layer of incidents (no attributes).
    fn point_layer(pts: &[(f64, f64)]) -> String {
        let mut l = Layer::new("inc")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        for (x, y) in pts {
            l.add_feature(Some(Geometry::Point(Coord::xy(*x, *y))), &[])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = OptimizedHotSpotAnalysisTool.run(&args, &ctx()).unwrap();
        let l = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, l)
    }

    /// A dense cluster of incidents against a sparse background yields a positive
    /// Gi* hot-spot bin somewhere.
    #[test]
    fn detects_hot_cluster() {
        let mut pts = Vec::new();
        // Dense 6x6 cluster near the origin.
        for i in 0..6 {
            for j in 0..6 {
                pts.push((i as f64 * 0.5, j as f64 * 0.5));
            }
        }
        // Sparse scatter far away.
        for k in 0..6 {
            pts.push((100.0 + k as f64 * 10.0, 100.0 + k as f64 * 10.0));
        }
        let (out, l) = run(json!({ "input": point_layer(&pts), "aggregation": "fishnet" }));
        assert!(
            out.outputs["hot_spots"].as_i64().unwrap() > 0,
            "expected at least one hot-spot cell"
        );
        // Output carries the four result fields.
        for f in ["VALUE", "GiZScore", "GiPValue", "Gi_Bin"] {
            assert!(l.schema.field_index(f).is_some(), "missing field {f}");
        }
    }

    /// distance_band and cell_size overrides are honoured and reported.
    #[test]
    fn honours_overrides() {
        let pts: Vec<(f64, f64)> = (0..20).map(|i| (i as f64, (i % 3) as f64)).collect();
        let (out, _l) = run(json!({
            "input": point_layer(&pts),
            "cell_size": 2.0, "distance_band": 5.0,
        }));
        assert_eq!(out.outputs["cell_size"].as_f64().unwrap(), 2.0);
        assert_eq!(out.outputs["distance_band"].as_f64().unwrap(), 5.0);
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            OptimizedHotSpotAnalysisTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "aggregation": "kmeans" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "distance_band": -1 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "aggregation": "snap" })).is_ok());
    }
}
