//! GeoLibre tool: create spatial sampling locations over a study area.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Create Spatial Sampling Locations*
//! (Data Management): generate sample points over a polygon study area with the
//! classic survey-design schemes. The repo already has GRTS spatially-balanced
//! sampling (`create_spatially_balanced_points`) and the bundled suite has
//! `random_points_in_polygon`, but there is no stratified, systematic, or
//! cluster sampling with allocation rules — the everyday design toolkit.
//!
//! Methods:
//!
//! - `simple_random` — uniform points anywhere in the study area (rejection
//!   sampling in the bounding box), with an optional `min_distance` spacing.
//! - `stratified` — partition the area into strata (each input feature, or
//!   features grouped by `strata_field`), allocate `num_samples` across strata
//!   by `allocation` (`equal`, `proportional` to area, or by a `population_field`
//!   total), then sample uniformly within each.
//! - `systematic` — a regular lattice of points at `bin_size` spacing clipped to
//!   the study area, in a `square`, `hexagon`, or `triangle` arrangement.
//! - `cluster` — tessellate the extent into `bin_size` square bins, randomly
//!   pick `num_clusters` bins that overlap the area, and sample points within
//!   each (tagged with a `cluster` id).
//!
//! All randomness is a seeded splitmix64 stream, so output is deterministic
//! (WASM-safe — no wall-clock or thread RNG). Point-in-polygon tests use `geo`'s
//! `Contains` (pure Rust). Output is a point layer tagged with the sampling
//! `method`, plus `stratum`/`cluster` where relevant.

use std::collections::BTreeMap;

use geo::{Area, Contains, Coord as GeoCoord, LineString, MultiPolygon, Point, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Feature, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CreateSpatialSamplingLocationsTool;

impl Tool for CreateSpatialSamplingLocationsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "create_spatial_sampling_locations",
            display_name: "Create Spatial Sampling Locations",
            summary: "Generate sample points over a polygon study area by simple-random, stratified, systematic (square/hexagon/triangle), or cluster design, with seeded reproducible randomness.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Study-area polygon vector file path (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output point vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "method",
                    description: "'simple_random' (default), 'stratified', 'systematic', or 'cluster'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "num_samples",
                    description: "Total number of sample points for simple_random, stratified, and cluster. Default 100.",
                    required: false,
                },
                ToolParamSpec {
                    name: "strata_field",
                    description: "For stratified: group input features into strata by this attribute. If omitted, each input feature is its own stratum.",
                    required: false,
                },
                ToolParamSpec {
                    name: "allocation",
                    description: "For stratified: 'proportional' to stratum area (default), 'equal' per stratum, or 'population_field' proportional to a field total.",
                    required: false,
                },
                ToolParamSpec {
                    name: "population_field",
                    description: "For allocation=population_field: the numeric field whose per-stratum total drives allocation.",
                    required: false,
                },
                ToolParamSpec {
                    name: "bin_shape",
                    description: "For systematic: 'square' (default), 'hexagon', or 'triangle' lattice.",
                    required: false,
                },
                ToolParamSpec {
                    name: "bin_size",
                    description: "Spacing (CRS units) for systematic lattices and cluster bins. Required for those methods.",
                    required: false,
                },
                ToolParamSpec {
                    name: "num_clusters",
                    description: "For cluster: how many bins to randomly select. Default 10.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_distance",
                    description: "Minimum spacing between simple_random / stratified / cluster points (CRS units). Default 0 (no constraint).",
                    required: false,
                },
                ToolParamSpec {
                    name: "seed",
                    description: "Seed for the deterministic RNG (default 1).",
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

        let layer = load_input_layer(input)?;
        let schema = layer.schema.clone();

        // Collect polygon features (each remembered so strata can be built).
        let mut parts: Vec<(usize, MultiPolygon)> = Vec::new();
        for (idx, f) in layer.features.iter().enumerate() {
            if let Some(mp) = f.geometry.as_ref().and_then(to_multipolygon) {
                if mp.unsigned_area() > 0.0 {
                    parts.push((idx, mp));
                }
            }
        }
        if parts.is_empty() {
            return Err(ToolError::Validation(
                "study area has no polygon features".to_string(),
            ));
        }

        let mut rng = Rng::new(prm.seed);
        let mut points: Vec<SamplePoint> = Vec::new();

        match prm.method {
            Method::SimpleRandom => {
                let whole = MultiPolygon(parts.iter().flat_map(|(_, mp)| mp.0.clone()).collect());
                sample_in(
                    &whole,
                    prm.num_samples,
                    prm.min_distance,
                    &mut rng,
                    None,
                    None,
                    &mut points,
                );
            }
            Method::Stratified => {
                let strata = build_strata(&parts, &layer, &schema, &prm)?;
                let allocations = allocate(&strata, prm.num_samples, &prm)?;
                for (stratum, n) in strata.iter().zip(allocations) {
                    sample_in(
                        &stratum.geom,
                        n,
                        prm.min_distance,
                        &mut rng,
                        Some(stratum.name.clone()),
                        None,
                        &mut points,
                    );
                }
            }
            Method::Systematic => {
                let whole = MultiPolygon(parts.iter().flat_map(|(_, mp)| mp.0.clone()).collect());
                systematic(&whole, prm.bin_size, prm.bin_shape, &mut points);
            }
            Method::Cluster => {
                let whole = MultiPolygon(parts.iter().flat_map(|(_, mp)| mp.0.clone()).collect());
                cluster(
                    &whole,
                    prm.bin_size,
                    prm.num_clusters,
                    prm.num_samples,
                    prm.min_distance,
                    &mut rng,
                    &mut points,
                );
            }
        }

        ctx.progress
            .info(&format!("placed {} sample point(s)", points.len()));

        // Build the output layer.
        let mut out = Layer::new("sampling_locations");
        out.geom_type = Some(GeometryType::Point);
        out.crs = layer.crs.clone();
        out.add_field(FieldDef::new("id", FieldType::Integer));
        out.add_field(FieldDef::new("method", FieldType::Text));
        out.add_field(FieldDef::new("stratum", FieldType::Text));
        out.add_field(FieldDef::new("cluster", FieldType::Integer));
        for (i, p) in points.iter().enumerate() {
            let mut f = Feature::with_geometry(i as u64, Geometry::point(p.x, p.y), 4);
            f.set_by_index(0, FieldValue::Integer(i as i64));
            f.set_by_index(1, FieldValue::Text(prm.method.as_str().to_string()));
            f.set_by_index(
                2,
                p.stratum
                    .clone()
                    .map(FieldValue::Text)
                    .unwrap_or(FieldValue::Null),
            );
            f.set_by_index(
                3,
                p.cluster
                    .map(FieldValue::Integer)
                    .unwrap_or(FieldValue::Null),
            );
            out.push(f);
        }

        let feature_count = out.len();
        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("method".to_string(), json!(prm.method.as_str()));
        Ok(ToolRunResult { outputs })
    }
}

/// A placed sample point with optional stratum / cluster tags.
struct SamplePoint {
    x: f64,
    y: f64,
    stratum: Option<String>,
    cluster: Option<i64>,
}

// ── Sampling primitives ──────────────────────────────────────────────────────

/// Rejection-samples up to `n` uniform points inside `mp`, honouring
/// `min_distance` (a grid-hash spacing check). Tags each with `stratum`/`cluster`.
#[allow(clippy::too_many_arguments)]
fn sample_in(
    mp: &MultiPolygon,
    n: usize,
    min_distance: f64,
    rng: &mut Rng,
    stratum: Option<String>,
    cluster: Option<i64>,
    out: &mut Vec<SamplePoint>,
) {
    if n == 0 {
        return;
    }
    let Some((minx, miny, maxx, maxy)) = bbox(mp) else {
        return;
    };
    let mut spacing = SpacingGuard::new(min_distance);
    let max_attempts = (n as u64 + 8) * 400;
    let mut placed = 0usize;
    let mut attempts = 0u64;
    while placed < n && attempts < max_attempts {
        attempts += 1;
        let x = minx + rng.f64() * (maxx - minx);
        let y = miny + rng.f64() * (maxy - miny);
        if !mp.contains(&Point::new(x, y)) {
            continue;
        }
        if !spacing.accept(x, y) {
            continue;
        }
        out.push(SamplePoint {
            x,
            y,
            stratum: stratum.clone(),
            cluster,
        });
        placed += 1;
    }
}

/// A regular lattice of points at `bin_size` spacing clipped to `mp`.
fn systematic(mp: &MultiPolygon, bin_size: f64, shape: BinShape, out: &mut Vec<SamplePoint>) {
    let Some((minx, miny, maxx, maxy)) = bbox(mp) else {
        return;
    };
    let s = bin_size;
    let push_if_in = |x: f64, y: f64, out: &mut Vec<SamplePoint>| {
        if mp.contains(&Point::new(x, y)) {
            out.push(SamplePoint {
                x,
                y,
                stratum: None,
                cluster: None,
            });
        }
    };
    match shape {
        BinShape::Square => {
            let mut y = miny + s * 0.5;
            while y <= maxy {
                let mut x = minx + s * 0.5;
                while x <= maxx {
                    push_if_in(x, y, out);
                    x += s;
                }
                y += s;
            }
        }
        BinShape::Hexagon | BinShape::Triangle => {
            // Rows offset by half a step; row spacing = s*sqrt(3)/2 gives a
            // hexagonal/triangular lattice of equal nearest-neighbour distance.
            let row_h = s * (3.0_f64).sqrt() / 2.0;
            let mut y = miny + s * 0.5;
            let mut row = 0usize;
            while y <= maxy {
                let offset = if row % 2 == 1 { s * 0.5 } else { 0.0 };
                let mut x = minx + s * 0.5 + offset;
                while x <= maxx {
                    push_if_in(x, y, out);
                    x += s;
                }
                y += row_h;
                row += 1;
            }
        }
    }
}

/// Cluster sampling: tessellate the extent into `bin_size` square bins that
/// overlap `mp`, randomly pick `num_clusters` of them, and sample points inside.
fn cluster(
    mp: &MultiPolygon,
    bin_size: f64,
    num_clusters: usize,
    num_samples: usize,
    min_distance: f64,
    rng: &mut Rng,
    out: &mut Vec<SamplePoint>,
) {
    let Some((minx, miny, maxx, maxy)) = bbox(mp) else {
        return;
    };
    let s = bin_size;
    // Candidate bins whose centre lies in the study area (a cheap overlap proxy).
    let mut bins: Vec<(f64, f64)> = Vec::new();
    let mut y = miny + s * 0.5;
    while y <= maxy + s * 0.5 {
        let mut x = minx + s * 0.5;
        while x <= maxx + s * 0.5 {
            if mp.contains(&Point::new(x, y)) {
                bins.push((x, y));
            }
            x += s;
        }
        y += s;
    }
    if bins.is_empty() {
        return;
    }
    // Random selection without replacement (partial Fisher–Yates).
    let k = num_clusters.min(bins.len());
    for i in 0..k {
        let j = i + (rng.next_u64() as usize) % (bins.len() - i);
        bins.swap(i, j);
    }
    let per = (num_samples as f64 / k as f64).ceil() as usize;
    for (ci, &(bx, by)) in bins.iter().take(k).enumerate() {
        // The bin as a small MultiPolygon intersected implicitly by contains().
        let bin_poly = square_polygon(bx, by, s);
        // Sample within the bin, keeping only points inside the study area.
        let mut tmp = Vec::new();
        let bin_mp = MultiPolygon(vec![bin_poly]);
        sample_in(
            &bin_mp,
            per * 3,
            min_distance,
            rng,
            None,
            Some(ci as i64),
            &mut tmp,
        );
        let mut kept = 0usize;
        for p in tmp {
            if kept >= per {
                break;
            }
            if mp.contains(&Point::new(p.x, p.y)) {
                out.push(p);
                kept += 1;
            }
        }
    }
}

fn square_polygon(cx: f64, cy: f64, s: f64) -> Polygon {
    let h = s * 0.5;
    Polygon::new(
        LineString::new(vec![
            GeoCoord {
                x: cx - h,
                y: cy - h,
            },
            GeoCoord {
                x: cx + h,
                y: cy - h,
            },
            GeoCoord {
                x: cx + h,
                y: cy + h,
            },
            GeoCoord {
                x: cx - h,
                y: cy + h,
            },
            GeoCoord {
                x: cx - h,
                y: cy - h,
            },
        ]),
        vec![],
    )
}

/// Enforces a minimum spacing via a grid hash of accepted points.
struct SpacingGuard {
    min_d: f64,
    cell: f64,
    grid: std::collections::HashMap<(i64, i64), Vec<(f64, f64)>>,
}

impl SpacingGuard {
    fn new(min_d: f64) -> Self {
        Self {
            min_d,
            cell: if min_d > 0.0 { min_d } else { 1.0 },
            grid: std::collections::HashMap::new(),
        }
    }
    fn accept(&mut self, x: f64, y: f64) -> bool {
        if self.min_d <= 0.0 {
            return true;
        }
        let (cx, cy) = (
            (x / self.cell).floor() as i64,
            (y / self.cell).floor() as i64,
        );
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(pts) = self.grid.get(&(cx + dx, cy + dy)) {
                    for &(px, py) in pts {
                        if (px - x).hypot(py - y) < self.min_d {
                            return false;
                        }
                    }
                }
            }
        }
        self.grid.entry((cx, cy)).or_default().push((x, y));
        true
    }
}

// ── Strata ───────────────────────────────────────────────────────────────────

struct Stratum {
    name: String,
    geom: MultiPolygon,
    area: f64,
    population: f64,
}

fn build_strata(
    parts: &[(usize, MultiPolygon)],
    layer: &Layer,
    schema: &wbvector::Schema,
    prm: &Params,
) -> Result<Vec<Stratum>, ToolError> {
    let pop_idx = match &prm.population_field {
        Some(name) => Some(schema.field_index(name).ok_or_else(|| {
            ToolError::Validation(format!("population_field '{name}' not found"))
        })?),
        None => None,
    };
    let strata_idx =
        match &prm.strata_field {
            Some(name) => Some(schema.field_index(name).ok_or_else(|| {
                ToolError::Validation(format!("strata_field '{name}' not found"))
            })?),
            None => None,
        };

    let mut map: BTreeMap<String, Stratum> = BTreeMap::new();
    for (idx, mp) in parts {
        let f = &layer.features[*idx];
        let name = match strata_idx {
            Some(si) => f
                .attributes
                .get(si)
                .map(field_value_string)
                .unwrap_or_default(),
            None => format!("feature_{idx}"),
        };
        let pop = pop_idx
            .and_then(|pi| f.attributes.get(pi).and_then(FieldValue::as_f64))
            .unwrap_or(0.0);
        let e = map.entry(name.clone()).or_insert_with(|| Stratum {
            name,
            geom: MultiPolygon(vec![]),
            area: 0.0,
            population: 0.0,
        });
        e.geom.0.extend(mp.0.clone());
        e.area += mp.unsigned_area();
        e.population += pop;
    }
    Ok(map.into_values().collect())
}

/// Allocates `total` samples across strata by the chosen rule.
fn allocate(strata: &[Stratum], total: usize, prm: &Params) -> Result<Vec<usize>, ToolError> {
    let n = strata.len();
    if n == 0 {
        return Ok(vec![]);
    }
    let weights: Vec<f64> = match prm.allocation {
        Allocation::Equal => vec![1.0; n],
        Allocation::Proportional => strata.iter().map(|s| s.area).collect(),
        Allocation::Population => {
            if prm.population_field.is_none() {
                return Err(ToolError::Validation(
                    "allocation=population_field requires 'population_field'".to_string(),
                ));
            }
            strata.iter().map(|s| s.population).collect()
        }
    };
    let sum: f64 = weights.iter().sum();
    if sum <= 0.0 {
        // Fall back to equal split.
        return Ok(largest_remainder(&vec![1.0; n], total));
    }
    Ok(largest_remainder(&weights, total))
}

/// Largest-remainder apportionment of `total` across positive `weights`.
fn largest_remainder(weights: &[f64], total: usize) -> Vec<usize> {
    let sum: f64 = weights.iter().sum();
    if sum <= 0.0 {
        return vec![0; weights.len()];
    }
    let raw: Vec<f64> = weights.iter().map(|w| w / sum * total as f64).collect();
    let mut counts: Vec<usize> = raw.iter().map(|r| r.floor() as usize).collect();
    let assigned: usize = counts.iter().sum();
    let mut remainders: Vec<(usize, f64)> = raw
        .iter()
        .enumerate()
        .map(|(i, r)| (i, r - r.floor()))
        .collect();
    remainders.sort_by(|a, b| b.1.total_cmp(&a.1));
    let mut left = total.saturating_sub(assigned);
    for (i, _) in remainders {
        if left == 0 {
            break;
        }
        counts[i] += 1;
        left -= 1;
    }
    counts
}

// ── Geometry helpers ─────────────────────────────────────────────────────────

fn bbox(mp: &MultiPolygon) -> Option<(f64, f64, f64, f64)> {
    let mut b = (
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );
    for poly in mp {
        for c in poly.exterior() {
            b.0 = b.0.min(c.x);
            b.1 = b.1.min(c.y);
            b.2 = b.2.max(c.x);
            b.3 = b.3.max(c.y);
        }
    }
    if b.0.is_finite() {
        Some(b)
    } else {
        None
    }
}

fn to_multipolygon(geom: &Geometry) -> Option<MultiPolygon> {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => Some(MultiPolygon(vec![rings_to_polygon(exterior, interiors)])),
        Geometry::MultiPolygon(parts) => Some(MultiPolygon(
            parts
                .iter()
                .map(|(ext, ints)| rings_to_polygon(ext, ints))
                .collect(),
        )),
        _ => None,
    }
}

fn rings_to_polygon(exterior: &Ring, interiors: &[Ring]) -> Polygon {
    Polygon::new(
        ring_to_linestring(exterior),
        interiors.iter().map(ring_to_linestring).collect(),
    )
}

fn ring_to_linestring(ring: &Ring) -> LineString {
    LineString::new(
        ring.coords()
            .iter()
            .map(|c| GeoCoord { x: c.x, y: c.y })
            .collect(),
    )
}

fn field_value_string(v: &FieldValue) -> String {
    match v {
        FieldValue::Text(s) | FieldValue::Date(s) | FieldValue::DateTime(s) => s.clone(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) => f.to_string(),
        FieldValue::Boolean(b) => b.to_string(),
        FieldValue::Null | FieldValue::Blob(_) => String::new(),
    }
}

// ── Deterministic RNG (splitmix64) ──────────────────────────────────────────

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

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Method {
    SimpleRandom,
    Stratified,
    Systematic,
    Cluster,
}

impl Method {
    fn as_str(self) -> &'static str {
        match self {
            Method::SimpleRandom => "simple_random",
            Method::Stratified => "stratified",
            Method::Systematic => "systematic",
            Method::Cluster => "cluster",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum BinShape {
    Square,
    Hexagon,
    Triangle,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Allocation {
    Equal,
    Proportional,
    Population,
}

struct Params {
    method: Method,
    num_samples: usize,
    strata_field: Option<String>,
    allocation: Allocation,
    population_field: Option<String>,
    bin_shape: BinShape,
    bin_size: f64,
    num_clusters: usize,
    min_distance: f64,
    seed: u64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let method = match parse_optional_str(args, "method")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("simple_random") | Some("random") => Method::SimpleRandom,
        Some("stratified") => Method::Stratified,
        Some("systematic") => Method::Systematic,
        Some("cluster") => Method::Cluster,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "unknown method '{o}' (simple_random, stratified, systematic, cluster)"
            )))
        }
    };
    let num_samples = parse_optional_u64(args, "num_samples")?.unwrap_or(100) as usize;
    let strata_field = parse_optional_str(args, "strata_field")?.map(str::to_string);
    let allocation = match parse_optional_str(args, "allocation")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("proportional") => Allocation::Proportional,
        Some("equal") => Allocation::Equal,
        Some("population_field") | Some("population") => Allocation::Population,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "unknown allocation '{o}' (proportional, equal, population_field)"
            )))
        }
    };
    let population_field = parse_optional_str(args, "population_field")?.map(str::to_string);
    let bin_shape = match parse_optional_str(args, "bin_shape")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("square") => BinShape::Square,
        Some("hexagon") | Some("hex") => BinShape::Hexagon,
        Some("triangle") => BinShape::Triangle,
        Some(o) => {
            return Err(ToolError::Validation(format!(
                "unknown bin_shape '{o}' (square, hexagon, triangle)"
            )))
        }
    };
    let bin_size = parse_optional_f64(args, "bin_size")?.unwrap_or(0.0);
    let num_clusters = parse_optional_u64(args, "num_clusters")?.unwrap_or(10) as usize;
    let min_distance = parse_optional_f64(args, "min_distance")?.unwrap_or(0.0);
    let seed = parse_optional_u64(args, "seed")?.unwrap_or(1);

    if matches!(method, Method::Systematic | Method::Cluster)
        && !(bin_size.is_finite() && bin_size > 0.0)
    {
        return Err(ToolError::Validation(format!(
            "method '{}' requires a positive 'bin_size'",
            method.as_str()
        )));
    }
    if min_distance < 0.0 {
        return Err(ToolError::Validation(
            "'min_distance' must be non-negative".to_string(),
        ));
    }
    Ok(Params {
        method,
        num_samples,
        strata_field,
        allocation,
        population_field,
        bin_shape,
        bin_size,
        num_clusters,
        min_distance,
        seed,
    })
}

fn parse_optional_f64(args: &ToolArgs, key: &str) -> Result<Option<f64>, ToolError> {
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

fn parse_optional_u64(args: &ToolArgs, key: &str) -> Result<Option<u64>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64().or_else(|| n.as_f64().map(|f| f.max(0.0) as u64))),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => s
            .trim()
            .parse::<u64>()
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
    use wbvector::{memory_store, Coord, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Vec<Coord> {
        vec![
            Coord::xy(x0, y0),
            Coord::xy(x1, y0),
            Coord::xy(x1, y1),
            Coord::xy(x0, y1),
        ]
    }

    fn run(input: &str, args: serde_json::Value) -> (ToolRunResult, Layer) {
        let mut m = args.as_object().unwrap().clone();
        m.insert("input".to_string(), json!(input));
        let args: ToolArgs = serde_json::from_value(Value::Object(m)).unwrap();
        let out = CreateSpatialSamplingLocationsTool
            .run(&args, &ctx())
            .unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn square_layer() -> String {
        let mut layer = Layer::new("area");
        layer
            .add_feature(
                Some(Geometry::polygon(rect(0.0, 0.0, 100.0, 100.0), vec![])),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn inside(layer: &Layer) -> bool {
        let poly = MultiPolygon(vec![Polygon::new(
            LineString::new(
                rect(0.0, 0.0, 100.0, 100.0)
                    .iter()
                    .map(|c| GeoCoord { x: c.x, y: c.y })
                    .collect(),
            ),
            vec![],
        )]);
        layer.features.iter().all(|f| {
            if let Some(Geometry::Point(c)) = &f.geometry {
                poly.contains(&Point::new(c.x, c.y))
            } else {
                false
            }
        })
    }

    /// Simple random: correct count, all inside, and deterministic by seed.
    #[test]
    fn simple_random_inside_and_deterministic() {
        let a = square_layer();
        let (out, layer) = run(
            &a,
            json!({ "method": "simple_random", "num_samples": 50, "seed": 7 }),
        );
        assert_eq!(out.outputs["feature_count"], json!(50));
        assert!(inside(&layer), "all points inside the study area");
        // Same seed -> identical first point.
        let (_, layer2) = run(
            &a,
            json!({ "method": "simple_random", "num_samples": 50, "seed": 7 }),
        );
        let p0 = |l: &Layer| match &l.features[0].geometry {
            Some(Geometry::Point(c)) => (c.x, c.y),
            _ => (0.0, 0.0),
        };
        assert_eq!(p0(&layer), p0(&layer2), "deterministic for a fixed seed");
    }

    /// min_distance is respected (all pairwise distances >= the constraint).
    #[test]
    fn min_distance_is_respected() {
        let a = square_layer();
        let (_, layer) = run(
            &a,
            json!({ "method": "simple_random", "num_samples": 30, "min_distance": 10.0, "seed": 3 }),
        );
        let pts: Vec<(f64, f64)> = layer
            .features
            .iter()
            .filter_map(|f| match &f.geometry {
                Some(Geometry::Point(c)) => Some((c.x, c.y)),
                _ => None,
            })
            .collect();
        for i in 0..pts.len() {
            for j in i + 1..pts.len() {
                let d = (pts[i].0 - pts[j].0).hypot(pts[i].1 - pts[j].1);
                assert!(d >= 10.0 - 1e-9, "pair {i},{j} too close: {d}");
            }
        }
    }

    /// Systematic square lattice: regular spacing, count ~ area / bin_size^2.
    #[test]
    fn systematic_square_grid() {
        let a = square_layer();
        let (out, _) = run(
            &a,
            json!({ "method": "systematic", "bin_shape": "square", "bin_size": 10.0 }),
        );
        // 100x100 area, 10-unit spacing -> ~10x10 = 100 points.
        let n = out.outputs["feature_count"].as_u64().unwrap();
        assert!(
            (90..=110).contains(&n),
            "systematic grid produced {n} points"
        );
    }

    /// Stratified proportional: two strata of area 3:1 get ~3:1 samples.
    #[test]
    fn stratified_proportional_allocation() {
        let mut layer = Layer::new("area");
        layer.add_field(FieldDef::new("zone", FieldType::Text));
        // Big stratum (area 300) and small (area 100).
        layer
            .add_feature(
                Some(Geometry::polygon(rect(0.0, 0.0, 30.0, 10.0), vec![])),
                &[("zone", "big".into())],
            )
            .unwrap();
        layer
            .add_feature(
                Some(Geometry::polygon(rect(0.0, 20.0, 10.0, 30.0), vec![])),
                &[("zone", "small".into())],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run(
            &input,
            json!({ "method": "stratified", "strata_field": "zone", "allocation": "proportional", "num_samples": 40, "seed": 5 }),
        );
        assert_eq!(out.outputs["feature_count"], json!(40));
        let count = |name: &str| {
            layer
                .features
                .iter()
                .filter(|f| {
                    f.get(&layer.schema, "stratum").unwrap() == &FieldValue::Text(name.into())
                })
                .count()
        };
        // 3:1 area -> ~30 vs ~10.
        assert_eq!(count("big"), 30);
        assert_eq!(count("small"), 10);
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = CreateSpatialSamplingLocationsTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(
            bad(json!({ "input": "x.geojson", "method": "systematic" })).is_err(),
            "systematic needs bin_size"
        );
        assert!(bad(json!({ "input": "x.geojson", "method": "bogus" })).is_err());
        assert!(
            bad(json!({ "input": "x.geojson", "method": "systematic", "bin_size": 5 })).is_ok()
        );
    }
}
