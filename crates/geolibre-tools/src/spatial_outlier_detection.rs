//! GeoLibre tool: score points by their Local Outlier Factor (LOF).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Spatial Outlier Detection* (Spatial
//! Statistics). The bundled suite offers density-based *clustering* (`dbscan`,
//! and the GeoLibre `hdbscan`) but no continuous per-point outlier score:
//! DBSCAN's binary "noise" label has no degree, and `lidar_remove_outliers` is
//! elevation-specific. LOF gives each point a score for how isolated it is
//! relative to the local density of its neighbours — the standard way to spot
//! anomalous locations.
//!
//! For each point the tool finds its `neighbors` nearest neighbours, then
//! computes the k-distance, reachability distance, local reachability density
//! (LRD), and finally `LOF = mean(LRD(neighbour) / LRD(point))`. LOF ≈ 1 is an
//! inlier; LOF ≫ 1 is an outlier. Points are flagged as outliers when their LOF
//! ranks in the top `percent_outlier`% (or exceeds an explicit `threshold`).
//! Distances are haversine metres for a geographic CRS, CRS units otherwise.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct SpatialOutlierDetectionTool;

impl Tool for SpatialOutlierDetectionTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "spatial_outlier_detection",
            display_name: "Spatial Outlier Detection",
            summary: "Score points by their Local Outlier Factor (LOF) — how isolated each point is relative to the local density of its k nearest neighbours — and flag the top outliers, like ArcGIS Spatial Outlier Detection. Fills the gap DBSCAN's binary noise label leaves.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input point vector layer.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output points with LOF score and outlier flag. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighbors",
                    description: "Number of nearest neighbours k (default 20).",
                    required: false,
                },
                ToolParamSpec {
                    name: "percent_outlier",
                    description: "Flag the top this-percent of points by LOF as outliers (default 5). Ignored if 'threshold' is set.",
                    required: false,
                },
                ToolParamSpec {
                    name: "threshold",
                    description: "Explicit LOF cutoff: points with LOF above this are outliers (overrides percent_outlier).",
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
        parse_params(args)?;
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
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let geographic = layer.crs_epsg().map(|e| e == 4326).unwrap_or(true);

        // Collect point coordinates; features without a point get LOF 0 / not flagged.
        let pts: Vec<Option<(f64, f64)>> = layer
            .features
            .iter()
            .map(|f| f.geometry.as_ref().and_then(point_xy))
            .collect();
        let idx: Vec<usize> = (0..pts.len()).filter(|&i| pts[i].is_some()).collect();
        let m = idx.len();
        let k = prm.neighbors.min(m.saturating_sub(1)).max(1);

        ctx.progress.info(&format!("LOF over {m} point(s), k={k}"));

        // For every point: sorted neighbour list (index, distance).
        let coords: Vec<(f64, f64)> = idx.iter().map(|&i| pts[i].unwrap()).collect();
        let mut neigh: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
        let mut kdist = vec![0.0f64; m];
        for a in 0..m {
            let mut d: Vec<(usize, f64)> = (0..m)
                .filter(|&b| b != a)
                .map(|b| (b, distance(coords[a], coords[b], geographic)))
                .collect();
            d.sort_by(|x, y| x.1.total_cmp(&y.1));
            d.truncate(k);
            kdist[a] = d.last().map(|(_, dd)| *dd).unwrap_or(0.0);
            neigh[a] = d;
        }

        // Local reachability density.
        let mut lrd = vec![0.0f64; m];
        for a in 0..m {
            let sum_reach: f64 = neigh[a].iter().map(|&(b, dab)| kdist[b].max(dab)).sum();
            lrd[a] = if sum_reach > 0.0 {
                neigh[a].len() as f64 / sum_reach
            } else {
                f64::INFINITY // identical coincident points -> infinite density
            };
        }

        // LOF = mean(lrd(neighbour) / lrd(point)).
        let mut lof = vec![1.0f64; m];
        for a in 0..m {
            if neigh[a].is_empty() {
                continue;
            }
            let s: f64 = neigh[a]
                .iter()
                .map(|&(b, _)| safe_ratio(lrd[b], lrd[a]))
                .sum();
            lof[a] = s / neigh[a].len() as f64;
        }

        // Decide the outlier flag.
        let cutoff = match prm.threshold {
            Some(t) => t,
            None => {
                let mut sorted = lof.clone();
                sorted.sort_by(|x, y| y.total_cmp(x)); // descending
                let n_flag = ((prm.percent_outlier / 100.0) * m as f64).ceil() as usize;
                if n_flag == 0 || sorted.is_empty() {
                    f64::INFINITY
                } else {
                    sorted[n_flag.min(sorted.len()) - 1]
                }
            }
        };

        // Attach lof + is_outlier back to the source features.
        let mut lof_all = vec![0.0f64; layer.features.len()];
        let mut flag_all = vec![0i64; layer.features.len()];
        let mut outliers = 0usize;
        for (a, &orig) in idx.iter().enumerate() {
            lof_all[orig] = lof[a];
            let is_out = lof[a] >= cutoff && cutoff.is_finite();
            flag_all[orig] = is_out as i64;
            if is_out {
                outliers += 1;
            }
        }

        layer.add_field(FieldDef::new("lof", FieldType::Float));
        layer.add_field(FieldDef::new("is_outlier", FieldType::Integer));
        for (i, f) in layer.features.iter_mut().enumerate() {
            f.attributes.push(FieldValue::Float(lof_all[i]));
            f.attributes.push(FieldValue::Integer(flag_all[i]));
        }

        let out_path = write_or_store_layer(layer, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("point_count".to_string(), json!(m));
        outputs.insert("outliers".to_string(), json!(outliers));
        Ok(ToolRunResult { outputs })
    }
}

fn safe_ratio(num: f64, den: f64) -> f64 {
    if den == 0.0 {
        if num == 0.0 {
            1.0
        } else {
            f64::INFINITY
        }
    } else if num.is_infinite() && den.is_infinite() {
        1.0
    } else {
        num / den
    }
}

fn distance(a: (f64, f64), b: (f64, f64), geographic: bool) -> f64 {
    if geographic {
        haversine(a.1, a.0, b.1, b.0)
    } else {
        (a.0 - b.0).hypot(a.1 - b.1)
    }
}

fn haversine(lat0: f64, lon0: f64, lat1: f64, lon1: f64) -> f64 {
    const R: f64 = 6_371_000.0;
    let (p0, p1) = (lat0.to_radians(), lat1.to_radians());
    let dphi = (lat1 - lat0).to_radians();
    let dlmb = (lon1 - lon0).to_radians();
    let a = (dphi / 2.0).sin().powi(2) + p0.cos() * p1.cos() * (dlmb / 2.0).sin().powi(2);
    2.0 * R * a.sqrt().asin()
}

fn point_xy(geom: &Geometry) -> Option<(f64, f64)> {
    match geom {
        Geometry::Point(c) => Some((c.x, c.y)),
        Geometry::MultiPoint(cs) if !cs.is_empty() => Some((cs[0].x, cs[0].y)),
        _ => None,
    }
}

struct Params {
    neighbors: usize,
    percent_outlier: f64,
    threshold: Option<f64>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let neighbors = opt_usize(args, "neighbors")?.unwrap_or(20).max(1);
    let percent_outlier = opt_f64(args, "percent_outlier")?.unwrap_or(5.0);
    if !(0.0..=100.0).contains(&percent_outlier) {
        return Err(ToolError::Validation(
            "'percent_outlier' must be between 0 and 100".to_string(),
        ));
    }
    let threshold = opt_f64(args, "threshold")?;
    if let Some(t) = threshold {
        if t <= 0.0 {
            return Err(ToolError::Validation(
                "'threshold' must be positive".to_string(),
            ));
        }
    }
    Ok(Params {
        neighbors,
        percent_outlier,
        threshold,
    })
}

fn opt_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
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

fn opt_usize(args: &ToolArgs, key: &str) -> Result<Option<usize>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64().map(|v| v as usize)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
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
    use wbvector::{memory_store, GeometryType, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn point_layer(pts: &[(f64, f64)]) -> String {
        let mut l = Layer::new("p")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("id", FieldType::Integer));
        for (i, (x, y)) in pts.iter().enumerate() {
            l.add_feature(Some(Geometry::point(*x, *y)), &[("id", (i as i64).into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SpatialOutlierDetectionTool.run(&args, &ctx()).unwrap();
        let l = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, l)
    }

    /// A point far from a tight cluster gets the highest LOF and is flagged.
    #[test]
    fn detects_the_obvious_outlier() {
        let mut pts: Vec<(f64, f64)> = Vec::new();
        // Tight 6x6 grid cluster near the origin.
        for r in 0..6 {
            for c in 0..6 {
                pts.push((c as f64, r as f64));
            }
        }
        // One far-away outlier.
        pts.push((100.0, 100.0));
        let outlier_i = pts.len() - 1;
        let (out, l) = run(json!({
            "input": point_layer(&pts), "neighbors": 5, "percent_outlier": 5,
        }));
        assert!(out.outputs["outliers"].as_i64().unwrap() >= 1);
        let lof = l.schema.field_index("lof").unwrap();
        let flag = l.schema.field_index("is_outlier").unwrap();
        let lofs: Vec<f64> = l
            .features
            .iter()
            .map(|f| f.attributes[lof].as_f64().unwrap())
            .collect();
        // The far point has the maximum LOF...
        let argmax = (0..lofs.len())
            .max_by(|&a, &b| lofs[a].total_cmp(&lofs[b]))
            .unwrap();
        assert_eq!(
            argmax, outlier_i,
            "the far point should have the highest LOF"
        );
        // ...and is flagged.
        assert_eq!(l.features[outlier_i].attributes[flag].as_i64(), Some(1));
    }

    /// A uniform grid has LOF near 1 for all points (no outliers by structure).
    #[test]
    fn uniform_grid_has_low_lof() {
        let mut pts = Vec::new();
        for r in 0..8 {
            for c in 0..8 {
                pts.push((c as f64, r as f64));
            }
        }
        let (_out, l) =
            run(json!({ "input": point_layer(&pts), "neighbors": 8, "threshold": 2.0 }));
        let lof = l.schema.field_index("lof").unwrap();
        let interior_ok = l
            .features
            .iter()
            .all(|f| f.attributes[lof].as_f64().unwrap() < 3.0);
        assert!(
            interior_ok,
            "a uniform grid should not have extreme LOF values"
        );
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            SpatialOutlierDetectionTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "percent_outlier": 150 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "threshold": -1 })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "neighbors": 10 })).is_ok());
    }
}
