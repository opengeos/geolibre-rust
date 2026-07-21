//! GeoLibre tool: Lee's L bivariate spatial association between two variables.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Bivariate Spatial Association*
//! (Spatial Statistics). The bundled suite covers univariate spatial
//! autocorrelation thoroughly (`global_morans_i`, `local_morans_i_lisa`,
//! `getis_ord_gi_star`) but nothing measures the spatial association *between
//! two* continuous variables. The GeoLibre `colocation_analysis` handles
//! categorical point colocations; this is its continuous-field counterpart
//! (e.g. income vs pollution across tracts).
//!
//! Lee's L combines Pearson correlation with spatial smoothing. Using
//! row-standardised weights over each feature's `neighbors` nearest neighbours,
//! let `zx`, `zy` be the mean-centred variables and `lx`, `ly` their spatial
//! lags. Global `L = Σ lx·ly / (√Σzx² · √Σzy²)`; the local
//! `L_i = n·lx_i·ly_i / (√Σzx² · √Σzy²)` shows where the two variables co-vary.
//! Each feature is classified High-High / High-Low / Low-High / Low-Low from the
//! sign of its own `zx` and its lagged `zy`. Significance is a seeded
//! permutation test (reproducible in WASM).

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct BivariateSpatialAssociationTool;

impl Tool for BivariateSpatialAssociationTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "bivariate_spatial_association",
            display_name: "Bivariate Spatial Association",
            summary: "Lee's L global and local statistic measuring where two continuous variables co-vary spatially, with a seeded permutation-test p-value and a High-High/High-Low/Low-High/Low-Low class per feature — like ArcGIS Bivariate Spatial Association. The continuous-field counterpart of colocation_analysis.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer with two numeric fields.",
                    required: true,
                },
                ToolParamSpec {
                    name: "x_field",
                    description: "First numeric variable field.",
                    required: true,
                },
                ToolParamSpec {
                    name: "y_field",
                    description: "Second numeric variable field.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output layer with local L, p-value, and class fields. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "neighbors",
                    description: "Number of nearest neighbours for the spatial weights (default 8).",
                    required: false,
                },
                ToolParamSpec {
                    name: "permutations",
                    description: "Number of permutations for the significance test (default 199; 0 disables).",
                    required: false,
                },
                ToolParamSpec {
                    name: "seed",
                    description: "Seed for the permutation RNG (default 1), for reproducible results.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "x_field")?;
        require_str(args, "y_field")?;
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let xi = layer
            .schema
            .field_index(&prm.x_field)
            .ok_or_else(|| ToolError::Validation(format!("x_field '{}' not found", prm.x_field)))?;
        let yi = layer
            .schema
            .field_index(&prm.y_field)
            .ok_or_else(|| ToolError::Validation(format!("y_field '{}' not found", prm.y_field)))?;

        // Gather feature centroid + (x, y) values; skip incomplete rows.
        let mut used: Vec<usize> = Vec::new();
        let mut cx = Vec::new();
        let mut cy = Vec::new();
        let mut xv = Vec::new();
        let mut yv = Vec::new();
        for (fi, f) in layer.features.iter().enumerate() {
            let (Some(c), Some(x), Some(y)) = (
                f.geometry.as_ref().and_then(centroid),
                f.attributes.get(xi).and_then(|v| v.as_f64()),
                f.attributes.get(yi).and_then(|v| v.as_f64()),
            ) else {
                continue;
            };
            used.push(fi);
            cx.push(c.0);
            cy.push(c.1);
            xv.push(x);
            yv.push(y);
        }
        let n = used.len();
        if n < 3 {
            return Err(ToolError::Execution(
                "need at least 3 features with geometry and both values".to_string(),
            ));
        }
        let k = prm.neighbors.min(n - 1).max(1);

        // Row-standardised k-NN weights as neighbour index lists (weight = 1/k).
        let neighbors: Vec<Vec<usize>> = (0..n)
            .map(|a| {
                let mut d: Vec<(f64, usize)> = (0..n)
                    .filter(|&b| b != a)
                    .map(|b| ((cx[a] - cx[b]).hypot(cy[a] - cy[b]), b))
                    .collect();
                d.sort_by(|p, q| p.0.total_cmp(&q.0));
                d.truncate(k);
                d.into_iter().map(|(_, b)| b).collect()
            })
            .collect();

        ctx.progress.info(&format!(
            "Lee's L over {n} feature(s), k={k}, {} permutation(s)",
            prm.permutations
        ));

        // Mean-centre.
        let mx = xv.iter().sum::<f64>() / n as f64;
        let my = yv.iter().sum::<f64>() / n as f64;
        let zx: Vec<f64> = xv.iter().map(|v| v - mx).collect();
        let zy: Vec<f64> = yv.iter().map(|v| v - my).collect();
        let ssx: f64 = zx.iter().map(|v| v * v).sum();
        let ssy: f64 = zy.iter().map(|v| v * v).sum();
        let denom = (ssx.sqrt() * ssy.sqrt()).max(1e-300);

        let (local, global) = lee_l(&zx, &zy, &neighbors, denom, n);

        // Classification from own zx and lagged zy.
        let classes: Vec<&str> = (0..n)
            .map(|a| {
                let lag_y = mean_lag(&zy, &neighbors[a]);
                match (zx[a] >= 0.0, lag_y >= 0.0) {
                    (true, true) => "HH",
                    (true, false) => "HL",
                    (false, true) => "LH",
                    (false, false) => "LL",
                }
            })
            .collect();

        // Permutation test: reshuffle y, recompute local L magnitudes.
        let mut p_local = vec![f64::NAN; n];
        let mut global_p = f64::NAN;
        if prm.permutations > 0 {
            let mut ge_local = vec![1usize; n]; // count |perm| >= |obs| (obs itself)
            let mut ge_global = 1usize;
            let mut perm_zy = zy.clone();
            let mut rng = prm.seed;
            for _ in 0..prm.permutations {
                fisher_yates(&mut perm_zy, &mut rng);
                let (ploc, pglob) = lee_l(&zx, &perm_zy, &neighbors, denom, n);
                if pglob.abs() >= global.abs() {
                    ge_global += 1;
                }
                for a in 0..n {
                    if ploc[a].abs() >= local[a].abs() {
                        ge_local[a] += 1;
                    }
                }
            }
            let d = (prm.permutations + 1) as f64;
            for a in 0..n {
                p_local[a] = ge_local[a] as f64 / d;
            }
            global_p = ge_global as f64 / d;
        }

        // Write fields.
        layer.add_field(FieldDef::new("local_l", FieldType::Float));
        layer.add_field(FieldDef::new("lee_class", FieldType::Text));
        layer.add_field(FieldDef::new("lee_p", FieldType::Float));
        let mut per_used = 0usize;
        // Build lookup from original feature index -> position in `used`.
        let mut pos = vec![usize::MAX; layer.features.len()];
        for (a, &fi) in used.iter().enumerate() {
            pos[fi] = a;
        }
        for (fi, f) in layer.features.iter_mut().enumerate() {
            let a = pos[fi];
            if a == usize::MAX {
                f.attributes.push(FieldValue::Float(0.0));
                f.attributes.push(FieldValue::Text("n/a".to_string()));
                f.attributes.push(FieldValue::Float(f64::NAN));
            } else {
                f.attributes.push(FieldValue::Float(local[a]));
                f.attributes.push(FieldValue::Text(classes[a].to_string()));
                f.attributes.push(FieldValue::Float(p_local[a]));
                per_used += 1;
            }
        }

        let out_path = write_or_store_layer(layer, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(per_used));
        outputs.insert("global_l".to_string(), json!(global));
        outputs.insert("global_p".to_string(), json!(global_p));
        Ok(ToolRunResult { outputs })
    }
}

/// Local and global Lee's L for centred variables `zx`, `zy`.
fn lee_l(
    zx: &[f64],
    zy: &[f64],
    neighbors: &[Vec<usize>],
    denom: f64,
    n: usize,
) -> (Vec<f64>, f64) {
    let mut local = vec![0.0f64; n];
    let mut sum = 0.0f64;
    for a in 0..n {
        let lx = mean_lag(zx, &neighbors[a]);
        let ly = mean_lag(zy, &neighbors[a]);
        let prod = lx * ly;
        local[a] = n as f64 * prod / denom;
        sum += prod;
    }
    let global = sum / denom;
    (local, global)
}

/// Row-standardised spatial lag = mean of the values over the neighbour set.
fn mean_lag(z: &[f64], nbrs: &[usize]) -> f64 {
    if nbrs.is_empty() {
        return 0.0;
    }
    nbrs.iter().map(|&b| z[b]).sum::<f64>() / nbrs.len() as f64
}

fn centroid(geom: &Geometry) -> Option<(f64, f64)> {
    let bb = geom.bbox()?;
    Some(((bb.min_x + bb.max_x) / 2.0, (bb.min_y + bb.max_y) / 2.0))
}

fn fisher_yates<T>(v: &mut [T], rng: &mut u64) {
    let n = v.len();
    for i in (1..n).rev() {
        let j = (next_u64(rng) % (i as u64 + 1)) as usize;
        v.swap(i, j);
    }
}

fn next_u64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

struct Params {
    x_field: String,
    y_field: String,
    neighbors: usize,
    permutations: usize,
    seed: u64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let x_field = require_str(args, "x_field")?.to_string();
    let y_field = require_str(args, "y_field")?.to_string();
    let neighbors = opt_usize(args, "neighbors")?.unwrap_or(8).max(1);
    let permutations = opt_usize(args, "permutations")?.unwrap_or(199);
    let seed = opt_usize(args, "seed")?.unwrap_or(1) as u64;
    Ok(Params {
        x_field,
        y_field,
        neighbors,
        permutations,
        seed,
    })
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

    fn layer(rows: &[(f64, f64, f64, f64)]) -> String {
        let mut l = Layer::new("z")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("x", FieldType::Float));
        l.add_field(FieldDef::new("y", FieldType::Float));
        for (px, py, x, y) in rows {
            l.add_feature(
                Some(Geometry::point(*px, *py)),
                &[("x", (*x).into()), ("y", (*y).into())],
            )
            .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = BivariateSpatialAssociationTool.run(&args, &ctx()).unwrap();
        let l = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, l)
    }

    /// Two variables with the same spatial gradient are strongly positively
    /// associated (global L > 0).
    #[test]
    fn positive_association_when_aligned() {
        // 6x6 grid; both x and y increase to the right -> co-vary spatially.
        let mut rows = Vec::new();
        for r in 0..6 {
            for c in 0..6 {
                let v = c as f64;
                rows.push((c as f64, r as f64, v, v + 0.1 * r as f64));
            }
        }
        let (out, _l) = run(json!({
            "input": layer(&rows), "x_field": "x", "y_field": "y",
            "neighbors": 4, "permutations": 99, "seed": 7,
        }));
        let l_global = out.outputs["global_l"].as_f64().unwrap();
        assert!(
            l_global > 0.3,
            "aligned gradients should give strong positive L, got {l_global}"
        );
    }

    /// Opposed gradients give a negative global L.
    #[test]
    fn negative_association_when_opposed() {
        let mut rows = Vec::new();
        for r in 0..6 {
            for c in 0..6 {
                rows.push((c as f64, r as f64, c as f64, (5 - c) as f64));
            }
        }
        let (out, _l) = run(json!({
            "input": layer(&rows), "x_field": "x", "y_field": "y",
            "neighbors": 4, "permutations": 0,
        }));
        assert!(out.outputs["global_l"].as_f64().unwrap() < 0.0);
    }

    /// Output carries a class label and local L per feature.
    #[test]
    fn writes_local_fields() {
        let mut rows = Vec::new();
        for r in 0..5 {
            for c in 0..5 {
                rows.push((c as f64, r as f64, c as f64, c as f64));
            }
        }
        let (_out, l) = run(json!({
            "input": layer(&rows), "x_field": "x", "y_field": "y", "permutations": 49,
        }));
        assert!(l.schema.field_index("local_l").is_some());
        assert!(l.schema.field_index("lee_class").is_some());
        assert!(l.schema.field_index("lee_p").is_some());
    }

    #[test]
    fn rejects_bad_params() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            BivariateSpatialAssociationTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "a.geojson", "x_field": "x" })).is_err());
        assert!(bad(json!({ "input": "a.geojson", "x_field": "x", "y_field": "y" })).is_ok());
    }
}
