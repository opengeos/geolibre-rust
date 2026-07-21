//! GeoLibre tool: spatially balanced (quasi-random) point sampling.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Create Spatially Balanced Points*
//! (Data Management). The bundled `random_points_in_polygon` is plain uniform
//! sampling — it clumps and leaves voids — and nothing does balanced /
//! quasi-random or inclusion-probability-weighted sampling, the standard design
//! for field surveys and monitoring networks. Pairs with
//! `generate_transects_along_lines` for survey design.
//!
//! Candidate locations come from a 2-D **Halton sequence** (bases 2 and 3) over
//! the constraint bounding box, seeded by skipping a deterministic number of
//! leading terms (WASM-safe, no system RNG). A candidate is accepted when it
//! falls inside the constraint polygons and — if a `probability` raster is
//! given — when a seeded uniform draw is below its inclusion probability
//! (rejection sampling, `probability` normalized by its maximum). Because the
//! Halton order is space-filling, the `sample_order` field means any prefix of
//! the output is itself spatially balanced.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbraster::Raster;
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::common::load_input_raster;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Give up after this many Halton candidates per requested point.
const MAX_ATTEMPTS_PER_POINT: usize = 2000;

pub struct CreateSpatiallyBalancedPointsTool;

impl Tool for CreateSpatiallyBalancedPointsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "create_spatially_balanced_points",
            display_name: "Create Spatially Balanced Points",
            summary: "Generate spatially well-spread (quasi-random Halton) sample points within a constraint polygon, optionally weighted by an inclusion-probability raster, each tagged with a balanced sample order, like ArcGIS Create Spatially Balanced Points.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "constraint",
                    description: "Constraint polygon layer; points are generated inside these polygons.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output point vector path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "num_points",
                    description: "Number of sample points to generate.",
                    required: true,
                },
                ToolParamSpec {
                    name: "probability",
                    description: "Optional inclusion-probability raster; higher cell values are more likely to be sampled (normalized by its max).",
                    required: false,
                },
                ToolParamSpec {
                    name: "seed",
                    description: "Random seed for reproducible sampling (default 1).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if args
            .get("constraint")
            .and_then(Value::as_str)
            .map(str::trim)
            .unwrap_or("")
            .is_empty()
        {
            return Err(ToolError::Validation(
                "missing required string parameter 'constraint'".to_string(),
            ));
        }
        parse_params(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let constraint = args.get("constraint").and_then(Value::as_str).unwrap();
        let output = parse_optional_str(args, "output")?;
        let prm = parse_params(args)?;

        let layer = load_input_layer(constraint)?;

        // Collect constraint polygons as (x,y) ring chains and the overall bbox.
        let mut polys: Vec<Vec<Vec<(f64, f64)>>> = Vec::new();
        let (mut minx, mut miny, mut maxx, mut maxy) = (
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        );
        for feature in &layer.features {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            for rings in polygon_regions(geom) {
                for ring in &rings {
                    for &(x, y) in ring {
                        minx = minx.min(x);
                        miny = miny.min(y);
                        maxx = maxx.max(x);
                        maxy = maxy.max(y);
                    }
                }
                polys.push(rings);
            }
        }
        if polys.is_empty() || !(maxx > minx && maxy > miny) {
            return Err(ToolError::Execution(
                "constraint layer has no usable polygon area".to_string(),
            ));
        }

        // Optional inclusion-probability raster.
        let prob = match &prm.probability_path {
            Some(path) => {
                let r = load_input_raster(path)?;
                let mut pmax = f64::NEG_INFINITY;
                for row in 0..r.rows as isize {
                    for col in 0..r.cols as isize {
                        let v = r.get(0, row, col);
                        if v != r.nodata && v.is_finite() {
                            pmax = pmax.max(v);
                        }
                    }
                }
                if pmax.is_nan() || pmax <= 0.0 {
                    return Err(ToolError::Execution(
                        "probability raster has no positive values".to_string(),
                    ));
                }
                Some((r, pmax))
            }
            None => None,
        };

        ctx.progress
            .info(&format!("sampling {} balanced point(s)", prm.num_points));

        // Seed: skip a deterministic number of leading Halton terms, and seed the
        // acceptance RNG.
        let mut rng = prm.seed ^ 0x2545_F491_4F6C_DD1D;
        let skip = 100 + (splitmix(&mut rng) % 4096) as usize;
        let mut accept_state = splitmix(&mut rng);

        let mut out = Layer::new("balanced_points").with_geom_type(GeometryType::Point);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("sample_order", FieldType::Integer));

        let max_attempts = prm
            .num_points
            .saturating_mul(MAX_ATTEMPTS_PER_POINT)
            .max(10_000);
        let mut accepted = 0usize;
        let mut i = 0usize;
        while accepted < prm.num_points && i < max_attempts {
            let idx = skip + i;
            i += 1;
            let hx = halton(idx as u64, 2);
            let hy = halton(idx as u64, 3);
            let x = minx + hx * (maxx - minx);
            let y = miny + hy * (maxy - miny);
            if !point_in_any(x, y, &polys) {
                continue;
            }
            if let Some((r, pmax)) = &prob {
                let p = sample_probability(r, x, y);
                let u = next_uniform(&mut accept_state);
                if u > p / pmax {
                    continue;
                }
            }
            out.add_feature(
                Some(Geometry::Point(Coord::xy(x, y))),
                &[("sample_order", FieldValue::Integer(accepted as i64))],
            )
            .map_err(|e| ToolError::Execution(format!("failed writing point: {e}")))?;
            accepted += 1;
        }

        if accepted == 0 {
            return Err(ToolError::Execution(
                "generated no points — check the constraint area or probability raster".to_string(),
            ));
        }

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("requested".to_string(), json!(prm.num_points));
        outputs.insert("generated".to_string(), json!(accepted));
        outputs.insert("attempts".to_string(), json!(i));
        Ok(ToolRunResult { outputs })
    }
}

// ── Halton sequence ───────────────────────────────────────────────────────────

/// The `index`-th term of the van der Corput / Halton sequence in `base`.
fn halton(mut index: u64, base: u64) -> f64 {
    let mut f = 1.0;
    let mut r = 0.0;
    while index > 0 {
        f /= base as f64;
        r += f * (index % base) as f64;
        index /= base;
    }
    r
}

fn splitmix(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn next_uniform(state: &mut u64) -> f64 {
    (splitmix(state) >> 11) as f64 / (1u64 << 53) as f64
}

// ── Geometry ──────────────────────────────────────────────────────────────────

fn sample_probability(r: &Raster, x: f64, y: f64) -> f64 {
    match r.world_to_pixel(x, y) {
        Some((col, row)) => {
            let v = r.get(0, row, col);
            if v == r.nodata || !v.is_finite() {
                0.0
            } else {
                v.max(0.0)
            }
        }
        None => 0.0,
    }
}

fn point_in_any(x: f64, y: f64, polys: &[Vec<Vec<(f64, f64)>>]) -> bool {
    polys.iter().any(|rings| point_in_rings(x, y, rings))
}

fn point_in_rings(x: f64, y: f64, rings: &[Vec<(f64, f64)>]) -> bool {
    if rings.is_empty() || !point_in_ring(x, y, &rings[0]) {
        return false;
    }
    !rings[1..].iter().any(|h| point_in_ring(x, y, h))
}

fn point_in_ring(x: f64, y: f64, ring: &[(f64, f64)]) -> bool {
    let n = ring.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = ring[i];
        let (xj, yj) = ring[j];
        if (yi > y) != (yj > y) {
            let xcross = (xj - xi) * (y - yi) / (yj - yi) + xi;
            if x < xcross {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

/// Each polygon part becomes one region: a Vec of rings (exterior first).
fn polygon_regions(geom: &Geometry) -> Vec<Vec<Vec<(f64, f64)>>> {
    let ring_pts =
        |ring: &Ring| -> Vec<(f64, f64)> { ring.coords().iter().map(|c| (c.x, c.y)).collect() };
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            let mut rings = vec![ring_pts(exterior)];
            rings.extend(interiors.iter().map(&ring_pts));
            vec![rings]
        }
        Geometry::MultiPolygon(parts) => parts
            .iter()
            .map(|(ext, holes)| {
                let mut rings = vec![ring_pts(ext)];
                rings.extend(holes.iter().map(&ring_pts));
                rings
            })
            .collect(),
        _ => Vec::new(),
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    num_points: usize,
    probability_path: Option<String>,
    seed: u64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let num_points = match args.get("num_points") {
        Some(Value::Number(n)) => n.as_u64().unwrap_or(0) as usize,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'num_points' must be an integer".into()))?,
        _ => {
            return Err(ToolError::Validation(
                "missing required integer parameter 'num_points'".to_string(),
            ))
        }
    };
    if num_points == 0 {
        return Err(ToolError::Validation(
            "'num_points' must be at least 1".to_string(),
        ));
    }
    let seed = match args.get("seed") {
        None | Some(Value::Null) => 1,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(1),
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
            .map_err(|_| ToolError::Validation("'seed' must be an integer".into()))?,
        Some(_) => return Err(ToolError::Validation("'seed' must be a number".into())),
    };
    Ok(Params {
        num_points,
        probability_path: parse_optional_str(args, "probability")?.map(str::to_string),
        seed,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::memory_store;

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn square_layer(x0: f64, y0: f64, s: f64) -> String {
        let mut l = Layer::new("c")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_feature(
            Some(Geometry::polygon(
                vec![
                    Coord::xy(x0, y0),
                    Coord::xy(x0 + s, y0),
                    Coord::xy(x0 + s, y0 + s),
                    Coord::xy(x0, y0 + s),
                ],
                vec![],
            )),
            &[],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CreateSpatiallyBalancedPointsTool
            .run(&args, &ctx())
            .unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// Generates the requested count, all inside the constraint, with 0..n order.
    #[test]
    fn generates_requested_points_inside() {
        let c = square_layer(0.0, 0.0, 100.0);
        let (out, layer) = run(json!({ "constraint": c, "num_points": 50 }));
        assert_eq!(out.outputs["generated"], json!(50));
        let oi = layer.schema.field_index("sample_order").unwrap();
        let mut orders = Vec::new();
        for f in layer.iter() {
            if let Some(Geometry::Point(cc)) = &f.geometry {
                assert!(cc.x >= 0.0 && cc.x <= 100.0 && cc.y >= 0.0 && cc.y <= 100.0);
            }
            orders.push(f.attributes[oi].as_i64().unwrap());
        }
        orders.sort_unstable();
        assert_eq!(orders, (0..50).collect::<Vec<_>>());
    }

    /// Same seed -> identical points; different seed -> different set.
    #[test]
    fn deterministic_by_seed() {
        let c = square_layer(0.0, 0.0, 100.0);
        let coords = |args: serde_json::Value| -> Vec<(u64, u64)> {
            let (_o, l) = run(args);
            l.iter()
                .filter_map(|f| match &f.geometry {
                    Some(Geometry::Point(c)) => Some((c.x.to_bits(), c.y.to_bits())),
                    _ => None,
                })
                .collect()
        };
        let a = coords(json!({ "constraint": c, "num_points": 30, "seed": 7 }));
        let b = coords(json!({ "constraint": c, "num_points": 30, "seed": 7 }));
        let d = coords(json!({ "constraint": c, "num_points": 30, "seed": 8 }));
        assert_eq!(a, b, "same seed reproducible");
        assert_ne!(a, d, "different seed differs");
    }

    /// The Halton sample is more spread than clumpy: every 10x10 subcell of a
    /// 100x100 square gets at least one of 100 points (uniform random rarely does).
    #[test]
    fn sample_is_well_spread() {
        let c = square_layer(0.0, 0.0, 100.0);
        let (_o, layer) = run(json!({ "constraint": c, "num_points": 100, "seed": 3 }));
        let mut occupied = std::collections::HashSet::new();
        for f in layer.iter() {
            if let Some(Geometry::Point(cc)) = &f.geometry {
                let cx = (cc.x / 10.0).floor() as i32;
                let cy = (cc.y / 10.0).floor() as i32;
                occupied.insert((cx.clamp(0, 9), cy.clamp(0, 9)));
            }
        }
        // Halton(2,3) fills the grid far more evenly than uniform random, whose
        // expected coverage for 100 points over 100 cells is only ~63.
        assert!(
            occupied.len() >= 75,
            "coverage {} of 100 cells",
            occupied.len()
        );
    }

    #[test]
    fn rejects_missing_num_points() {
        let c = square_layer(0.0, 0.0, 10.0);
        let args: ToolArgs = serde_json::from_value(json!({ "constraint": c })).unwrap();
        assert!(CreateSpatiallyBalancedPointsTool.validate(&args).is_err());
    }
}
