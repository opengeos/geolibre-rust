//! GeoLibre tool: generate stratified reference points for accuracy assessment.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Create Accuracy Assessment Points*
//! (Spatial Analyst / Image Analyst).
//!
//! ## The missing first half
//!
//! `classification_accuracy_assessment` **consumes** a set of reference points,
//! each carrying a ground-truth class, and returns the confusion matrix, kappa,
//! and per-class producer's/user's accuracy. Nothing produced those points, so
//! the validation workflow started with hand-authoring a point layer.
//!
//! `create_spatial_sampling_locations` looks like it should cover this and does
//! not: its strata are **polygon features or an attribute field**, so it cannot
//! stratify by the class values of a classified *raster*, and it has no
//! equalized allocation. Equalized allocation is the whole point here — a class
//! covering 1% of the map needs as many validation points as one covering 60%,
//! because per-class accuracy is estimated from the per-class sample.
//! The bundled `random_points_in_polygon` is unstratified.
//!
//! ## Output contract
//!
//! The output schema is exactly what `classification_accuracy_assessment`
//! expects, so the two chain directly:
//!
//! ```text
//! create_accuracy_assessment_points(input=classified.tif) -> points
//! # ... user fills GrndTruth in the field ...
//! classification_accuracy_assessment(points=points, class_field="GrndTruth",
//!                                    input=classified.tif)
//! ```
//!
//! `CLASSVALUE` is prefilled from the map. `GrndTruth` is written as the
//! sentinel `-1`, meaning "not yet collected" — a null would be indistinguishable
//! from a class the user genuinely could not determine.
//!
//! ## Determinism
//!
//! All randomness is a seeded splitmix64 stream, matching
//! `create_spatial_sampling_locations`. No wall-clock and no thread RNG, so the
//! WASM builds behave identically to native.

use std::collections::{BTreeMap, HashSet};

use serde_json::json;
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer};

use crate::args_common::{band_index, choice_or, req_str, usize_or};
use crate::common::load_input_raster;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

const STRATEGIES: [&str; 3] = ["stratified_random", "equalized_stratified", "random"];

/// Ground-truth sentinel for "not yet collected".
const UNCOLLECTED: i64 = -1;

pub struct CreateAccuracyAssessmentPointsTool;

impl Tool for CreateAccuracyAssessmentPointsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "create_accuracy_assessment_points",
            display_name: "Create Accuracy Assessment Points",
            summary: "Generates randomly distributed reference points stratified by the classes of a classified raster or polygon layer, prefilled with the map class and ready for ground-truth labelling (ArcGIS Create Accuracy Assessment Points). classification_accuracy_assessment consumes such points but nothing produced them, and create_spatial_sampling_locations stratifies by polygon feature or field rather than by raster class and has no equalized allocation.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Classified raster, or a classified polygon layer with 'class_field'.",
                    required: true,
                },
                ToolParamSpec {
                    name: "class_field",
                    description: "Class attribute, required when 'input' is a vector layer. Ignored for raster input, which uses the cell value.",
                    required: false,
                },
                ToolParamSpec {
                    name: "num_points",
                    description: "Total number of reference points (default 500).",
                    required: false,
                },
                ToolParamSpec {
                    name: "sampling_strategy",
                    description: "One of 'stratified_random' (default, allocate proportional to class area), 'equalized_stratified' (equal count per class), 'random' (ignore class).",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_points_per_class",
                    description: "Floor applied under 'stratified_random' so rare classes are not dropped (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "band",
                    description: "1-based band of a raster input (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "seed",
                    description: "Seed for the deterministic RNG (default 1).",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output point layer with CLASSVALUE (map class), GrndTruth (-1 until collected) and Stratum. If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        req_str(args, "input")?;
        choice_or(args, "sampling_strategy", &STRATEGIES, "stratified_random")?;
        band_index(args, "band")?;
        if usize_or(args, "num_points", 500)? == 0 {
            return Err(ToolError::Validation(
                "'num_points' must be at least 1".to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = req_str(args, "input")?.to_string();
        let class_field = parse_optional_str(args, "class_field")?.map(str::to_string);
        let num_points = usize_or(args, "num_points", 500)?;
        let strategy = choice_or(args, "sampling_strategy", &STRATEGIES, "stratified_random")?;
        let min_per_class = usize_or(args, "min_points_per_class", 1)?;
        let band = band_index(args, "band")?;
        let seed = args.get("seed").and_then(|v| v.as_u64()).unwrap_or(1);
        let output = parse_optional_str(args, "output")?;

        // Candidate locations, bucketed by class. Raster and vector inputs
        // differ only in how the buckets are filled.
        let (buckets, epsg, source_kind) =
            collect_candidates(&input, class_field.as_deref(), band)?;
        if buckets.is_empty() {
            return Err(ToolError::Execution(
                "the input contained no valid classified cells or features".to_string(),
            ));
        }
        let total_candidates: usize = buckets.values().map(|v| v.len()).sum();
        ctx.progress.info(&format!(
            "{source_kind}: {} class(es), {total_candidates} candidate location(s)",
            buckets.len()
        ));

        let allocation = allocate(&buckets, num_points, strategy, min_per_class);

        let mut rng = Rng::new(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ 0xDEAD_BEEF);
        let mut layer =
            Layer::new("accuracy_assessment_points").with_geom_type(GeometryType::Point);
        if let Some(e) = epsg {
            layer = layer.with_crs_epsg(e);
        }
        layer.add_field(FieldDef::new("CLASSVALUE", FieldType::Integer));
        layer.add_field(FieldDef::new("GrndTruth", FieldType::Integer));
        layer.add_field(FieldDef::new("Stratum", FieldType::Integer));

        let mut exhausted: Vec<i64> = Vec::new();
        let mut per_class: BTreeMap<String, usize> = BTreeMap::new();
        for (class, want) in &allocation {
            let pool = &buckets[class];
            // A class with fewer candidates than its allocation is sampled
            // exhaustively and reported. Looping until `want` is reached would
            // never terminate.
            let take = (*want).min(pool.len());
            if take < *want {
                exhausted.push(*class);
            }
            for &(x, y) in sample_without_replacement(pool, take, &mut rng).iter() {
                layer.push(Feature {
                    fid: 0,
                    geometry: Some(Geometry::point(x, y)),
                    attributes: vec![
                        FieldValue::Integer(*class),
                        FieldValue::Integer(UNCOLLECTED),
                        FieldValue::Integer(*class),
                    ],
                });
            }
            per_class.insert(class.to_string(), take);
        }

        let point_count = layer.features.len();
        let out_path = write_or_store_layer(layer, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("point_count".to_string(), json!(point_count));
        outputs.insert("class_count".to_string(), json!(buckets.len()));
        outputs.insert("sampling_strategy".to_string(), json!(strategy));
        outputs.insert("points_per_class".to_string(), json!(per_class));
        outputs.insert("classes_exhausted".to_string(), json!(exhausted));
        Ok(ToolRunResult { outputs })
    }
}

/// Distributes `num_points` across the classes.
///
/// `stratified_random` allocates in proportion to class extent with a
/// `min_per_class` floor; `equalized_stratified` gives every class the same
/// count; `random` pools everything into a single pseudo-class so the draw is
/// unstratified. Every branch corrects the rounding residual so the totals sum
/// to `num_points` exactly (when the candidates allow it).
fn allocate(
    buckets: &BTreeMap<i64, Vec<(f64, f64)>>,
    num_points: usize,
    strategy: &str,
    min_per_class: usize,
) -> Vec<(i64, usize)> {
    let classes: Vec<i64> = buckets.keys().copied().collect();
    let n = classes.len();

    let mut alloc: Vec<(i64, usize)> = match strategy {
        "equalized_stratified" => {
            let base = num_points / n;
            let mut rem = num_points % n;
            classes
                .iter()
                .map(|c| {
                    let extra = usize::from(rem > 0);
                    rem = rem.saturating_sub(1);
                    (*c, base + extra)
                })
                .collect()
        }
        "random" => {
            // Unstratified: allocate strictly in proportion to candidate count
            // with no floor, which is the same distribution a single pooled
            // uniform draw would produce in expectation.
            let total: usize = buckets.values().map(|v| v.len()).sum();
            classes
                .iter()
                .map(|c| {
                    let share = buckets[c].len() as f64 / total as f64;
                    (*c, (share * num_points as f64).floor() as usize)
                })
                .collect()
        }
        _ => {
            let total: usize = buckets.values().map(|v| v.len()).sum();
            classes
                .iter()
                .map(|c| {
                    let share = buckets[c].len() as f64 / total as f64;
                    let want = (share * num_points as f64).round() as usize;
                    (*c, want.max(min_per_class))
                })
                .collect()
        }
    };

    // Correct the residual so the total lands on num_points. Give (or take)
    // one point at a time, largest class first, and never push a class below
    // its floor or above its candidate count.
    let floor = if strategy == "stratified_random" {
        min_per_class
    } else {
        0
    };
    let mut order: Vec<usize> = (0..alloc.len()).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(buckets[&alloc[i].0].len()));

    let mut guard = 0_usize;
    loop {
        let sum: usize = alloc.iter().map(|(_, v)| *v).sum();
        if sum == num_points {
            break;
        }
        // Bounded by construction: each pass either changes the sum or exits.
        guard += 1;
        if guard > alloc.len() * num_points + alloc.len() + 1 {
            break;
        }
        let mut changed = false;
        if sum < num_points {
            for &i in &order {
                if alloc[i].1 < buckets[&alloc[i].0].len() {
                    alloc[i].1 += 1;
                    changed = true;
                    break;
                }
            }
        } else {
            for &i in order.iter().rev() {
                if alloc[i].1 > floor {
                    alloc[i].1 -= 1;
                    changed = true;
                    break;
                }
            }
        }
        // No class can absorb the difference (every class is at its candidate
        // ceiling, or every class is at its floor). Stop rather than spin.
        if !changed {
            break;
        }
    }
    alloc
}

/// Draws `take` distinct entries from `pool` using partial Fisher-Yates over an
/// index set, so no candidate is selected twice and the cost is O(take) rather
/// than O(pool) shuffling.
fn sample_without_replacement(pool: &[(f64, f64)], take: usize, rng: &mut Rng) -> Vec<(f64, f64)> {
    if take >= pool.len() {
        return pool.to_vec();
    }
    let mut chosen: HashSet<usize> = HashSet::with_capacity(take);
    let mut out = Vec::with_capacity(take);
    // Rejection sampling is fine while take is well below pool.len(); the
    // swap-based fallback covers the dense case without a full shuffle.
    if take * 2 <= pool.len() {
        while out.len() < take {
            let i = (rng.f64() * pool.len() as f64) as usize;
            let i = i.min(pool.len() - 1);
            if chosen.insert(i) {
                out.push(pool[i]);
            }
        }
    } else {
        let mut idx: Vec<usize> = (0..pool.len()).collect();
        for k in 0..take {
            let j = k + ((rng.f64() * (pool.len() - k) as f64) as usize).min(pool.len() - k - 1);
            idx.swap(k, j);
            out.push(pool[idx[k]]);
        }
    }
    out
}

/// Builds the per-class candidate buckets from a raster or vector input.
///
/// Raster: every valid cell is a candidate at its cell centre, keyed by the
/// rounded cell value (classes are categorical). Vector: every polygon feature
/// contributes its representative interior point, keyed by `class_field`.
#[allow(clippy::type_complexity)]
fn collect_candidates(
    input: &str,
    class_field: Option<&str>,
    band: isize,
) -> Result<(BTreeMap<i64, Vec<(f64, f64)>>, Option<u32>, &'static str), ToolError> {
    let mut buckets: BTreeMap<i64, Vec<(f64, f64)>> = BTreeMap::new();

    if let Ok(raster) = load_input_raster(input) {
        if band as usize >= raster.bands {
            return Err(ToolError::Validation(format!(
                "band {} out of range (raster has {} band(s))",
                band + 1,
                raster.bands
            )));
        }
        let rows = raster.rows;
        let cols = raster.cols;
        let y_max = raster.y_min + rows as f64 * raster.cell_size_y;
        for r in 0..rows {
            for c in 0..cols {
                let v = raster.get(band, r as isize, c as isize);
                if v == raster.nodata || !v.is_finite() {
                    continue;
                }
                let x = raster.x_min + (c as f64 + 0.5) * raster.cell_size_x;
                let y = y_max - (r as f64 + 0.5) * raster.cell_size_y;
                buckets.entry(v.round() as i64).or_default().push((x, y));
            }
        }
        return Ok((buckets, raster.crs.epsg, "raster"));
    }

    let layer = load_input_layer(input)?;
    let field = class_field.ok_or_else(|| {
        ToolError::Validation(
            "'class_field' is required when 'input' is a vector layer".to_string(),
        )
    })?;
    let fidx = layer.schema.field_index(field).ok_or_else(|| {
        ToolError::Validation(format!(
            "class_field '{field}' not found in the input layer"
        ))
    })?;
    for feature in layer.iter() {
        let Some(class) = feature.attributes.get(fidx).and_then(field_as_i64) else {
            continue;
        };
        if let Some((x, y)) = representative_point(feature.geometry.as_ref()) {
            buckets.entry(class).or_default().push((x, y));
        }
    }
    Ok((buckets, layer.crs_epsg(), "vector"))
}

/// A point guaranteed to lie on (or in) the geometry. Polygons use the centroid
/// when it falls inside and otherwise the first vertex, which keeps the sample
/// on the feature for concave shapes where the centroid escapes.
fn representative_point(geom: Option<&Geometry>) -> Option<(f64, f64)> {
    match geom? {
        Geometry::Point(p) => Some((p.x, p.y)),
        Geometry::MultiPoint(ps) => ps.first().map(|p| (p.x, p.y)),
        Geometry::Polygon { exterior, .. } => {
            let coords = &exterior.0;
            let (mut sx, mut sy) = (0.0, 0.0);
            for c in coords {
                sx += c.x;
                sy += c.y;
            }
            let n = coords.len().max(1) as f64;
            let centroid = (sx / n, sy / n);
            if crate::vector_common::ring_contains(coords, centroid.0, centroid.1) {
                Some(centroid)
            } else {
                coords.first().map(|c| (c.x, c.y))
            }
        }
        Geometry::MultiPolygon(polys) => polys.first().and_then(|(ext, _)| {
            representative_point(Some(&Geometry::Polygon {
                exterior: ext.clone(),
                interiors: Vec::new(),
            }))
        }),
        _ => None,
    }
}

fn field_as_i64(v: &FieldValue) -> Option<i64> {
    match v {
        FieldValue::Integer(i) => Some(*i),
        FieldValue::Float(f) if f.is_finite() => Some(f.round() as i64),
        FieldValue::Text(s) => s.trim().parse::<i64>().ok(),
        _ => None,
    }
}

// ── Deterministic RNG (splitmix64), matching create_spatial_sampling_locations ──

struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Self(seed)
    }
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbraster::{CrsInfo, DataType, Raster, RasterConfig};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn raster(rows: usize, cols: usize, data: &[f64]) -> String {
        let mut r = Raster::new(RasterConfig {
            cols,
            rows,
            bands: 1,
            x_min: 0.0,
            y_min: 0.0,
            cell_size: 1.0,
            cell_size_y: None,
            nodata: -9999.0,
            data_type: DataType::F64,
            crs: CrsInfo {
                epsg: Some(3857),
                wkt: None,
                proj4: None,
            },
            metadata: Vec::new(),
        });
        for row in 0..rows {
            for col in 0..cols {
                r.set(0, row as isize, col as isize, data[row * cols + col])
                    .unwrap();
            }
        }
        let id = wbraster::memory_store::put_raster(r);
        wbraster::memory_store::make_raster_memory_path(&id)
    }

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = CreateAccuracyAssessmentPointsTool
            .run(&args, &ctx())
            .unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    /// 10x10, class 1 covers 90 cells and class 2 covers 10.
    fn skewed() -> String {
        let mut data = vec![1.0; 100];
        for d in data.iter_mut().take(10) {
            *d = 2.0;
        }
        raster(10, 10, &data)
    }

    fn class_counts(layer: &Layer) -> BTreeMap<i64, usize> {
        let idx = layer.schema.field_index("CLASSVALUE").unwrap();
        let mut m = BTreeMap::new();
        for f in layer.iter() {
            if let FieldValue::Integer(c) = f.attributes[idx] {
                *m.entry(c).or_insert(0) += 1;
            }
        }
        m
    }

    #[test]
    fn produces_the_requested_number_of_points() {
        let (layer, res) = run(json!({"input": skewed(), "num_points": 40}));
        assert_eq!(layer.features.len(), 40);
        assert_eq!(res.outputs["point_count"], json!(40));
    }

    #[test]
    fn equalized_gives_every_class_the_same_count() {
        // The whole reason this tool exists: the rare class (10 cells of 100)
        // gets the same sample size as the dominant one, so its per-class
        // accuracy is estimable. 16 points is within both classes' ceilings.
        let (layer, _) = run(json!({
            "input": skewed(), "num_points": 16,
            "sampling_strategy": "equalized_stratified",
        }));
        let counts = class_counts(&layer);
        assert_eq!(counts[&1], 8);
        assert_eq!(counts[&2], 8);
    }

    #[test]
    fn equalized_still_caps_at_a_class_candidate_ceiling() {
        // Class 2 holds only 10 cells, so an equal share of 20 is impossible.
        // Capping (and reporting it) is right; inventing points is not.
        let (layer, res) = run(json!({
            "input": skewed(), "num_points": 40,
            "sampling_strategy": "equalized_stratified",
        }));
        let counts = class_counts(&layer);
        assert_eq!(counts[&1], 20);
        assert_eq!(counts[&2], 10, "capped at the number of class-2 cells");
        assert_eq!(res.outputs["classes_exhausted"], json!([2]));
    }

    #[test]
    fn stratified_random_follows_class_area() {
        let (layer, _) = run(json!({
            "input": skewed(), "num_points": 100,
            "sampling_strategy": "stratified_random",
        }));
        let counts = class_counts(&layer);
        // 90/10 split of the map, so roughly 90/10 of the sample.
        assert_eq!(counts[&1], 90);
        assert_eq!(counts[&2], 10);
    }

    #[test]
    fn min_points_per_class_rescues_a_rare_class() {
        // Class 2 holds 1 cell in 100. A pure proportional split of 10 points
        // rounds it to zero, which would make its accuracy unestimable.
        let mut data = vec![1.0; 100];
        data[0] = 2.0;
        let (layer, _) = run(json!({
            "input": raster(10, 10, &data), "num_points": 10,
            "min_points_per_class": 1,
        }));
        let counts = class_counts(&layer);
        assert_eq!(*counts.get(&2).unwrap_or(&0), 1, "rare class must survive");
    }

    #[test]
    fn ground_truth_starts_uncollected_and_stratum_matches_the_map_class() {
        let (layer, _) = run(json!({"input": skewed(), "num_points": 10}));
        let gt = layer.schema.field_index("GrndTruth").unwrap();
        let cv = layer.schema.field_index("CLASSVALUE").unwrap();
        let st = layer.schema.field_index("Stratum").unwrap();
        for f in layer.iter() {
            assert_eq!(f.attributes[gt], FieldValue::Integer(-1));
            assert_eq!(f.attributes[cv], f.attributes[st]);
        }
    }

    #[test]
    fn points_land_on_cells_of_the_class_they_claim() {
        // A point whose CLASSVALUE disagrees with the raster underneath it
        // would silently corrupt every downstream accuracy figure.
        let path = skewed();
        let r = load_input_raster(&path).unwrap();
        let (layer, _) = run(json!({"input": path, "num_points": 50}));
        let cv = layer.schema.field_index("CLASSVALUE").unwrap();
        let y_max = r.y_min + r.rows as f64 * r.cell_size_y;
        for f in layer.iter() {
            let Some(Geometry::Point(p)) = f.geometry.as_ref() else {
                panic!("expected points")
            };
            let col = ((p.x - r.x_min) / r.cell_size_x).floor() as isize;
            let row = ((y_max - p.y) / r.cell_size_y).floor() as isize;
            let under = r.get(0, row, col).round() as i64;
            assert_eq!(f.attributes[cv], FieldValue::Integer(under));
        }
    }

    #[test]
    fn no_candidate_cell_is_sampled_twice() {
        // Ask for every cell of a tiny raster; duplicates would show up as a
        // short distinct-coordinate set.
        let (layer, _) = run(json!({
            "input": raster(2, 2, &[1.0, 1.0, 1.0, 1.0]), "num_points": 4,
        }));
        let mut seen = HashSet::new();
        for f in layer.iter() {
            let Some(Geometry::Point(p)) = f.geometry.as_ref() else {
                panic!()
            };
            assert!(
                seen.insert((p.x.to_bits(), p.y.to_bits())),
                "duplicate sample location"
            );
        }
        assert_eq!(seen.len(), 4);
    }

    #[test]
    fn asking_for_more_points_than_cells_is_capped_and_reported() {
        let (layer, res) = run(json!({
            "input": raster(2, 2, &[1.0, 1.0, 1.0, 1.0]), "num_points": 999,
        }));
        assert_eq!(layer.features.len(), 4);
        assert_eq!(res.outputs["classes_exhausted"], json!([1]));
    }

    #[test]
    fn nodata_cells_are_never_sampled() {
        let (layer, res) = run(json!({
            "input": raster(1, 4, &[1.0, -9999.0, -9999.0, 2.0]), "num_points": 2,
        }));
        assert_eq!(res.outputs["class_count"], json!(2));
        for f in layer.iter() {
            let Some(Geometry::Point(p)) = f.geometry.as_ref() else {
                panic!()
            };
            // Only cells 0 and 3 are valid, at x = 0.5 and 3.5.
            assert!(p.x < 1.0 || p.x > 3.0, "sampled a no-data cell at {}", p.x);
        }
    }

    #[test]
    fn the_same_seed_reproduces_the_same_points() {
        let path = skewed();
        let a = run(json!({"input": path.clone(), "num_points": 20, "seed": 7})).0;
        let b = run(json!({"input": path.clone(), "num_points": 20, "seed": 7})).0;
        let c = run(json!({"input": path, "num_points": 20, "seed": 8})).0;
        let coords = |l: &Layer| -> Vec<(u64, u64)> {
            l.iter()
                .map(|f| match f.geometry.as_ref() {
                    Some(Geometry::Point(p)) => (p.x.to_bits(), p.y.to_bits()),
                    _ => panic!(),
                })
                .collect()
        };
        assert_eq!(coords(&a), coords(&b), "same seed must be reproducible");
        assert_ne!(coords(&a), coords(&c), "a different seed must differ");
    }

    #[test]
    fn a_vector_input_stratifies_by_class_field() {
        let mut l = Layer::new("polys")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("cls", FieldType::Integer));
        for i in 0..6_i64 {
            let x = i as f64;
            l.add_feature(
                Some(Geometry::polygon(
                    vec![
                        wbvector::Coord::xy(x, 0.0),
                        wbvector::Coord::xy(x + 1.0, 0.0),
                        wbvector::Coord::xy(x + 1.0, 1.0),
                        wbvector::Coord::xy(x, 1.0),
                        wbvector::Coord::xy(x, 0.0),
                    ],
                    Vec::new(),
                )),
                &[("cls", (i % 2).into())],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        let path = wbvector::memory_store::make_vector_memory_path(&id);
        let (layer, res) = run(json!({
            "input": path, "class_field": "cls", "num_points": 4,
            "sampling_strategy": "equalized_stratified",
        }));
        assert_eq!(res.outputs["class_count"], json!(2));
        let counts = class_counts(&layer);
        assert_eq!(counts[&0], 2);
        assert_eq!(counts[&1], 2);
    }

    #[test]
    fn rejects_bad_parameters() {
        let r = raster(1, 2, &[1.0, 2.0]);
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CreateAccuracyAssessmentPointsTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"input": r.clone(), "num_points": 0})));
        assert!(bad(
            json!({"input": r.clone(), "sampling_strategy": "nope"})
        ));
        // band is 1-based; 0 is a common off-by-one and must not silently pass.
        assert!(bad(json!({"input": r, "band": 0})));
    }
}
