//! GeoLibre tool: eliminate sliver polygons by merging them into a neighbor.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Eliminate* (Data Management): the
//! classic cleanup step after an overlay or a raster-to-vector conversion,
//! where the output is peppered with tiny sliver polygons that should be
//! absorbed into an adjacent, "real" polygon rather than deleted outright.
//!
//! In ArcGIS you first *select* the polygons to eliminate; here the selection
//! is expressed declaratively so the tool is a single, reproducible step:
//!
//! - `max_area` — a polygon is a sliver candidate when its area (in CRS units²)
//!   is at or below this threshold.
//! - `where` — a polygon is a candidate when it matches a simple attribute
//!   condition (`FIELD OP VALUE`, e.g. `class = 0` or `gridcode <= 2`).
//! - When both are given a polygon must satisfy *both* (small **and** matching).
//! - `exclude` — a second condition that protects features from elimination
//!   even if they would otherwise qualify (they stay, and can host merges).
//!
//! Each selected sliver is merged (unioned) into one neighboring non-sliver
//! polygon, chosen by `strategy`:
//!
//! - `longest_border` (default) — the neighbor sharing the longest boundary.
//! - `largest_area` — the neighbor with the largest area.
//!
//! The receiving polygon keeps its own attributes (the sliver's are dropped),
//! matching the ArcGIS tool. Slivers are absorbed iteratively, so a chain of
//! adjacent slivers collapses into the real polygon at the end of the chain: a
//! sliver touching only other slivers is merged once one of them has itself
//! merged into a host. A sliver with no reachable non-sliver neighbor (e.g. an
//! isolated island, or the case where every polygon was selected) is left
//! unchanged and reported as unmerged. Non-polygon features pass through
//! untouched.
//!
//! Shared-boundary length is measured by accumulating the collinear overlap of
//! boundary segments within `tolerance` (CRS units), which is exact for a
//! clean coverage (adjacent polygons share identical edges) and still robust
//! when overlay output splits a shared edge into several collinear pieces.
//!
//! The union itself is `geo`'s `BooleanOps` (pure Rust, no GEOS), so a merged
//! host that ends up with disjoint parts becomes a `MultiPolygon`.

use std::collections::{BTreeMap, HashSet};

use geo::{Area, BooleanOps, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, FieldValue, Geometry, GeometryType, Ring, Schema};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct EliminatePolygonsTool;

impl Tool for EliminatePolygonsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "eliminate_polygons",
            display_name: "Eliminate Polygons",
            summary: "Merge sliver polygons (by max area and/or an attribute query) into a neighboring polygon, chosen by longest shared border or largest area.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon vector file path, format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "max_area",
                    description: "Select polygons whose area (in CRS units squared) is at or below this value. May be combined with 'where' (both must match).",
                    required: false,
                },
                ToolParamSpec {
                    name: "where",
                    description: "Select polygons matching a simple attribute condition 'FIELD OP VALUE', where OP is one of = != < <= > >= (e.g. \"class = 0\"). May be combined with 'max_area'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "exclude",
                    description: "Protect features matching this condition (same 'FIELD OP VALUE' syntax) from elimination; they are retained and can receive merges.",
                    required: false,
                },
                ToolParamSpec {
                    name: "strategy",
                    description: "Which neighbor a sliver merges into: 'longest_border' (default) or 'largest_area'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Maximum perpendicular distance (in CRS units) for two boundary segments to count as a shared edge. Default 1e-6.",
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
        let layer_name = layer.name.clone();
        let layer_crs = layer.crs.clone();
        let input_count = layer.len();

        // Classify every feature: polygon candidate (sliver), polygon host, or
        // pass-through (non-polygon). Areas come from the geo geometry so the
        // area threshold and the 'largest_area' strategy agree.
        let mut sliver_geom: BTreeMap<usize, MultiPolygon> = BTreeMap::new();
        let mut host_acc: BTreeMap<usize, MultiPolygon> = BTreeMap::new();
        let mut candidate_order: Vec<usize> = Vec::new();
        for (idx, feature) in layer.features.iter().enumerate() {
            let Some(mp) = feature.geometry.as_ref().and_then(to_multipolygon) else {
                continue; // non-polygon: pass through untouched
            };
            let excluded = prm
                .exclude
                .as_ref()
                .is_some_and(|c| c.matches(feature, &schema));
            let selected = !excluded && prm.selects(&mp, feature, &schema);
            if selected {
                candidate_order.push(idx);
                sliver_geom.insert(idx, mp);
            } else {
                host_acc.insert(idx, mp);
            }
        }

        ctx.progress.info(&format!(
            "{} feature(s): {} sliver candidate(s), {} host polygon(s)",
            input_count,
            candidate_order.len(),
            host_acc.len()
        ));
        if host_acc.is_empty() && !candidate_order.is_empty() {
            ctx.progress.info(
                "warning: no host polygons available; slivers cannot be merged and are kept as-is",
            );
        }

        // Absorb slivers iteratively. Each pass, a sliver that shares a border
        // with some host is unioned into the best host (which then grows, so a
        // neighboring sliver can attach on a later pass). Stop when a pass
        // makes no progress; whatever remains has no reachable host.
        let mut remaining = candidate_order.clone();
        let mut merged_hosts: HashSet<usize> = HashSet::new();
        loop {
            let mut still: Vec<usize> = Vec::new();
            for &s in &remaining {
                let sliver = &sliver_geom[&s];
                let sliver_segs = Segments::of(sliver);
                let mut best: Option<Choice> = None;
                for (&h, hgeom) in &host_acc {
                    if !sliver_segs
                        .bbox
                        .intersects(&Segments::bbox_of(hgeom), prm.tolerance)
                    {
                        continue;
                    }
                    let shared = sliver_segs.shared_border_len(hgeom, prm.tolerance);
                    if shared <= prm.tolerance {
                        continue;
                    }
                    let cand = Choice {
                        host: h,
                        shared,
                        area: hgeom.unsigned_area(),
                    };
                    if best.as_ref().is_none_or(|b| cand.beats(b, prm.strategy)) {
                        best = Some(cand);
                    }
                }
                match best {
                    Some(choice) => {
                        let grown = host_acc[&choice.host].union(sliver);
                        host_acc.insert(choice.host, grown);
                        merged_hosts.insert(choice.host);
                    }
                    None => still.push(s),
                }
            }
            if still.len() == remaining.len() {
                remaining = still;
                break;
            }
            remaining = still;
        }

        let eliminated: HashSet<usize> = candidate_order
            .iter()
            .copied()
            .filter(|idx| !remaining.contains(idx))
            .collect();
        let unmerged = remaining.len();

        // Rebuild the layer in original feature order, dropping eliminated
        // slivers and swapping in the grown geometry for hosts that received a
        // merge. Everything else (retained slivers, untouched hosts,
        // non-polygons) is copied verbatim.
        let mut has_multipolygon = false;
        let mut out_features: Vec<Feature> = Vec::with_capacity(input_count - eliminated.len());
        for (idx, feature) in layer.features.into_iter().enumerate() {
            if eliminated.contains(&idx) {
                continue;
            }
            let mut feature = feature;
            if merged_hosts.contains(&idx) {
                let geom = multipolygon_to_geometry(&host_acc[&idx]);
                has_multipolygon |= matches!(geom, Geometry::MultiPolygon(_));
                feature.geometry = Some(geom);
            } else if matches!(feature.geometry, Some(Geometry::MultiPolygon(_))) {
                has_multipolygon = true;
            }
            feature.fid = out_features.len() as u64;
            out_features.push(feature);
        }

        let mut out_layer = wbvector::Layer::new(layer_name);
        out_layer.schema = schema;
        out_layer.crs = layer_crs;
        out_layer.features = out_features;
        // A union of adjacent polygons can produce a MultiPolygon even when the
        // input was single-part; widen the declared type so strict writers
        // (e.g. Shapefile) accept the output.
        out_layer.geom_type = if has_multipolygon {
            Some(GeometryType::MultiPolygon)
        } else {
            Some(GeometryType::Polygon)
        };

        ctx.progress.info(&format!(
            "eliminated {} sliver(s), {} kept (no host neighbor)",
            eliminated.len(),
            unmerged
        ));

        let feature_count = out_layer.len();
        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_count".to_string(), json!(input_count));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("eliminated_count".to_string(), json!(eliminated.len()));
        outputs.insert("unmerged_count".to_string(), json!(unmerged));
        outputs.insert("strategy".to_string(), json!(prm.strategy.as_str()));
        Ok(ToolRunResult { outputs })
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Strategy {
    LongestBorder,
    LargestArea,
}

impl Strategy {
    fn as_str(self) -> &'static str {
        match self {
            Self::LongestBorder => "longest_border",
            Self::LargestArea => "largest_area",
        }
    }
}

struct Params {
    max_area: Option<f64>,
    where_cond: Option<Condition>,
    exclude: Option<Condition>,
    strategy: Strategy,
    tolerance: f64,
}

impl Params {
    /// A polygon is selected when it satisfies every provided criterion. With
    /// both `max_area` and `where`, it must be small *and* match the query.
    fn selects(&self, mp: &MultiPolygon, feature: &Feature, schema: &Schema) -> bool {
        let area_ok = self.max_area.map(|m| mp.unsigned_area() <= m);
        let where_ok = self.where_cond.as_ref().map(|c| c.matches(feature, schema));
        match (area_ok, where_ok) {
            (Some(a), Some(w)) => a && w,
            (Some(a), None) => a,
            (None, Some(w)) => w,
            // parse_params guarantees at least one selector is present.
            (None, None) => false,
        }
    }
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let max_area = parse_optional_f64(args, "max_area")?;
    if let Some(m) = max_area {
        if !(m > 0.0 && m.is_finite()) {
            return Err(ToolError::Validation(
                "parameter 'max_area' must be a positive number".to_string(),
            ));
        }
    }
    let where_cond = parse_optional_str(args, "where")?
        .map(Condition::parse)
        .transpose()?;
    let exclude = parse_optional_str(args, "exclude")?
        .map(Condition::parse)
        .transpose()?;
    if max_area.is_none() && where_cond.is_none() {
        return Err(ToolError::Validation(
            "provide 'max_area' and/or 'where' to select the polygons to eliminate".to_string(),
        ));
    }
    let strategy = match parse_optional_str(args, "strategy")?
        .map(|s| s.trim().to_ascii_lowercase())
        .as_deref()
    {
        None | Some("longest_border") => Strategy::LongestBorder,
        Some("largest_area") => Strategy::LargestArea,
        Some(other) => {
            return Err(ToolError::Validation(format!(
                "unknown strategy '{other}' (expected longest_border or largest_area)"
            )))
        }
    };
    let tolerance = parse_optional_f64(args, "tolerance")?.unwrap_or(1e-6);
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'tolerance' must be a positive number".to_string(),
        ));
    }
    Ok(Params {
        max_area,
        where_cond,
        exclude,
        strategy,
        tolerance,
    })
}

/// Parses an optional numeric parameter, accepting a JSON number or a numeric
/// string (host UIs often post form values as strings).
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

// ── Attribute condition (FIELD OP VALUE) ───────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Op {
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
}

/// A single `FIELD OP VALUE` attribute filter. Deliberately minimal — one
/// condition, no boolean combinators — so it is unambiguous to parse and test;
/// richer selection can be done upstream with `extract_by_attribute`.
struct Condition {
    field: String,
    op: Op,
    value: String,
}

impl Condition {
    fn parse(raw: &str) -> Result<Self, ToolError> {
        // Two-char operators first so `<=` is not read as `<`.
        let (field, op, value) = if let Some((l, r)) = raw.split_once("!=") {
            (l, Op::Ne, r)
        } else if let Some((l, r)) = raw.split_once("<=") {
            (l, Op::Le, r)
        } else if let Some((l, r)) = raw.split_once(">=") {
            (l, Op::Ge, r)
        } else if let Some((l, r)) = raw.split_once("==") {
            (l, Op::Eq, r)
        } else if let Some((l, r)) = raw.split_once('=') {
            (l, Op::Eq, r)
        } else if let Some((l, r)) = raw.split_once('<') {
            (l, Op::Lt, r)
        } else if let Some((l, r)) = raw.split_once('>') {
            (l, Op::Gt, r)
        } else {
            return Err(ToolError::Validation(format!(
                "condition '{raw}' must be 'FIELD OP VALUE' with OP one of = != < <= > >="
            )));
        };
        let field = field.trim();
        let value = value.trim().trim_matches(|c| c == '\'' || c == '"');
        if field.is_empty() || value.is_empty() {
            return Err(ToolError::Validation(format!(
                "condition '{raw}' must have a non-empty field and value"
            )));
        }
        Ok(Self {
            field: field.to_string(),
            op,
            value: value.to_string(),
        })
    }

    /// Evaluates the condition against a feature. A missing field, or an
    /// ordering comparison on non-numeric operands, yields `false`.
    fn matches(&self, feature: &Feature, schema: &Schema) -> bool {
        let Ok(fv) = feature.get(schema, &self.field) else {
            return false;
        };
        // Numeric comparison when both sides look numeric.
        if let (Some(lhs), Ok(rhs)) = (fv.as_f64(), self.value.parse::<f64>()) {
            return match self.op {
                Op::Eq => lhs == rhs,
                Op::Ne => lhs != rhs,
                Op::Lt => lhs < rhs,
                Op::Le => lhs <= rhs,
                Op::Gt => lhs > rhs,
                Op::Ge => lhs >= rhs,
            };
        }
        // Text comparison otherwise; ordering falls back to lexicographic.
        let lhs = field_value_string(fv);
        match self.op {
            Op::Eq => lhs == self.value,
            Op::Ne => lhs != self.value,
            Op::Lt => lhs.as_str() < self.value.as_str(),
            Op::Le => lhs.as_str() <= self.value.as_str(),
            Op::Gt => lhs.as_str() > self.value.as_str(),
            Op::Ge => lhs.as_str() >= self.value.as_str(),
        }
    }
}

fn field_value_string(v: &FieldValue) -> String {
    match v {
        FieldValue::Text(s) => s.clone(),
        FieldValue::Date(s) | FieldValue::DateTime(s) => s.clone(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) => f.to_string(),
        FieldValue::Boolean(b) => b.to_string(),
        FieldValue::Null | FieldValue::Blob(_) => String::new(),
    }
}

// ── Neighbor choice ─────────────────────────────────────────────────────────

struct Choice {
    host: usize,
    shared: f64,
    area: f64,
}

impl Choice {
    /// True when `self` is a better host than `other` under the strategy.
    /// Ties break on the secondary measure, then on the lower feature index so
    /// the result is deterministic regardless of map iteration order.
    fn beats(&self, other: &Choice, strategy: Strategy) -> bool {
        let (a, b) = match strategy {
            Strategy::LongestBorder => (
                (self.shared, self.area, other.host),
                (other.shared, other.area, self.host),
            ),
            Strategy::LargestArea => (
                (self.area, self.shared, other.host),
                (other.area, other.shared, self.host),
            ),
        };
        a.0 > b.0
            || (a.0 == b.0 && a.1 > b.1)
            || (a.0 == b.0 && a.1 == b.1 && self.host < other.host)
    }
}

// ── Shared-border measurement ───────────────────────────────────────────────

#[derive(Clone, Copy)]
struct BBox {
    minx: f64,
    miny: f64,
    maxx: f64,
    maxy: f64,
}

impl BBox {
    fn empty() -> Self {
        Self {
            minx: f64::INFINITY,
            miny: f64::INFINITY,
            maxx: f64::NEG_INFINITY,
            maxy: f64::NEG_INFINITY,
        }
    }
    fn expand(&mut self, x: f64, y: f64) {
        self.minx = self.minx.min(x);
        self.miny = self.miny.min(y);
        self.maxx = self.maxx.max(x);
        self.maxy = self.maxy.max(y);
    }
    /// True when the boxes come within `pad` of each other on both axes.
    fn intersects(&self, other: &BBox, pad: f64) -> bool {
        self.minx <= other.maxx + pad
            && self.maxx >= other.minx - pad
            && self.miny <= other.maxy + pad
            && self.maxy >= other.miny - pad
    }
}

/// Boundary segments of a multipolygon (every ring of every part), with a
/// cached bounding box for cheap pruning.
struct Segments {
    segs: Vec<(GeoCoord, GeoCoord)>,
    bbox: BBox,
}

impl Segments {
    fn of(mp: &MultiPolygon) -> Self {
        let mut segs = Vec::new();
        let mut bbox = BBox::empty();
        for poly in mp {
            for ring in std::iter::once(poly.exterior()).chain(poly.interiors()) {
                push_ring_segments(ring, &mut segs, &mut bbox);
            }
        }
        Self { segs, bbox }
    }

    fn bbox_of(mp: &MultiPolygon) -> BBox {
        let mut bbox = BBox::empty();
        for poly in mp {
            for ring in std::iter::once(poly.exterior()).chain(poly.interiors()) {
                for c in &ring.0 {
                    bbox.expand(c.x, c.y);
                }
            }
        }
        bbox
    }

    /// Total length over which these segments are collinear (within `tol`) and
    /// overlapping with the boundary of `other`.
    fn shared_border_len(&self, other: &MultiPolygon, tol: f64) -> f64 {
        let mut other_segs = Vec::new();
        let mut other_bbox = BBox::empty();
        for poly in other {
            for ring in std::iter::once(poly.exterior()).chain(poly.interiors()) {
                push_ring_segments(ring, &mut other_segs, &mut other_bbox);
            }
        }
        let mut total = 0.0;
        for &(p1, p2) in &self.segs {
            for &(q1, q2) in &other_segs {
                total += collinear_overlap(p1, p2, q1, q2, tol);
            }
        }
        total
    }
}

fn push_ring_segments(ring: &LineString, segs: &mut Vec<(GeoCoord, GeoCoord)>, bbox: &mut BBox) {
    let pts = &ring.0;
    for w in pts.windows(2) {
        segs.push((w[0], w[1]));
        bbox.expand(w[0].x, w[0].y);
    }
    if let Some(last) = pts.last() {
        bbox.expand(last.x, last.y);
    }
}

/// Length over which segment `p1p2` and segment `q1q2` are collinear (both `q`
/// endpoints within `tol` of the infinite line through `p`) and overlapping.
/// Direction-agnostic, so a shared edge traversed in opposite windings (as two
/// adjacent polygons do) is still measured.
fn collinear_overlap(p1: GeoCoord, p2: GeoCoord, q1: GeoCoord, q2: GeoCoord, tol: f64) -> f64 {
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let len = dx.hypot(dy);
    if len <= tol {
        return 0.0;
    }
    let (ux, uy) = (dx / len, dy / len);
    // Perpendicular distance of each q endpoint from the p-line.
    let perp = |q: GeoCoord| ((q.x - p1.x) * uy - (q.y - p1.y) * ux).abs();
    if perp(q1) > tol || perp(q2) > tol {
        return 0.0;
    }
    // Project onto the unit direction; p spans [0, len].
    let proj = |q: GeoCoord| (q.x - p1.x) * ux + (q.y - p1.y) * uy;
    let (tq1, tq2) = (proj(q1), proj(q2));
    let lo = 0.0f64.max(tq1.min(tq2));
    let hi = len.min(tq1.max(tq2));
    (hi - lo).max(0.0)
}

// ── geo <-> wbvector geometry conversion ───────────────────────────────────

/// Converts a polygonal `wbvector` geometry to a `geo` `MultiPolygon`. Returns
/// `None` for non-polygon geometries (which the tool passes through untouched).
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
    // `geo` closes rings itself; the missing closing vertex in `Ring` is fine.
    LineString::new(
        ring.coords()
            .iter()
            .map(|c| GeoCoord { x: c.x, y: c.y })
            .collect(),
    )
}

/// Converts a `geo` `MultiPolygon` back to a `wbvector` geometry: a single part
/// becomes a `Polygon`, multiple parts a `MultiPolygon`.
fn multipolygon_to_geometry(mp: &MultiPolygon) -> Geometry {
    if mp.0.len() == 1 {
        let (exterior, interiors) = polygon_to_rings(&mp.0[0]);
        Geometry::Polygon {
            exterior,
            interiors,
        }
    } else {
        Geometry::MultiPolygon(mp.0.iter().map(polygon_to_rings).collect())
    }
}

fn polygon_to_rings(poly: &Polygon) -> (Ring, Vec<Ring>) {
    (
        linestring_to_ring(poly.exterior()),
        poly.interiors().iter().map(linestring_to_ring).collect(),
    )
}

fn linestring_to_ring(ls: &LineString) -> Ring {
    // Drop the closing duplicate vertex `geo` keeps; `Ring` stores it implicitly.
    let mut coords: Vec<Coord> = ls.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
    if coords.len() >= 2 && coords.first() == coords.last() {
        coords.pop();
    }
    Ring::new(coords)
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, FieldDef, FieldType, Layer};

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

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = EliminatePolygonsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn area_of(geom: &Geometry) -> f64 {
        to_multipolygon(geom)
            .map(|mp| mp.unsigned_area())
            .unwrap_or(0.0)
    }

    /// Two big rectangles side by side (share the x=10 edge) plus a thin sliver
    /// riding on the right block's right edge. Eliminating by area should merge
    /// the sliver into the right block.
    #[test]
    fn merges_sliver_into_only_neighbor() {
        let mut layer = Layer::new("parcels");
        layer.add_field(FieldDef::new("name", FieldType::Text));
        // left block: 0..10 x 0..10 (area 100)
        layer
            .add_feature(
                Some(Geometry::polygon(rect(0.0, 0.0, 10.0, 10.0), vec![])),
                &[("name", "left".into())],
            )
            .unwrap();
        // right block: 10..20 x 0..10 (area 100)
        layer
            .add_feature(
                Some(Geometry::polygon(rect(10.0, 0.0, 20.0, 10.0), vec![])),
                &[("name", "right".into())],
            )
            .unwrap();
        // sliver: 20..20.5 x 0..10 (area 5), touches right block on x=20
        layer
            .add_feature(
                Some(Geometry::polygon(rect(20.0, 0.0, 20.5, 10.0), vec![])),
                &[("name", "sliver".into())],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input, "max_area": 10.0 }));
        assert_eq!(out.outputs["eliminated_count"], json!(1));
        assert_eq!(out.outputs["feature_count"], json!(2));
        // The right block absorbed the sliver: area 100 -> 105.
        let right = layer
            .features
            .iter()
            .find(|f| f.get(&layer.schema, "name").unwrap() == &FieldValue::Text("right".into()))
            .unwrap();
        assert!((area_of(right.geometry.as_ref().unwrap()) - 105.0).abs() < 1e-6);
    }

    /// A sliver touching two blocks: a long shared edge with A, a short one with
    /// B. `longest_border` picks A; `largest_area` picks the bigger block B.
    #[test]
    fn strategy_picks_longest_border_vs_largest_area() {
        // A shares a length-4 edge with the sliver but is small; B shares only
        // a length-2 edge but is by far the larger polygon. The two strategies
        // therefore pick different hosts.
        let build2 = || {
            let mut layer = Layer::new("parcels");
            layer.add_field(FieldDef::new("name", FieldType::Text));
            // A: 0..10 x 0..10 (area 100), shares x=10 over y 0..4 (len 4).
            layer
                .add_feature(
                    Some(Geometry::polygon(rect(0.0, 0.0, 10.0, 10.0), vec![])),
                    &[("name", "A".into())],
                )
                .unwrap();
            // B: 10..60 x -20..0 (area 1000), shares y=0 over x 10..12 (len 2).
            layer
                .add_feature(
                    Some(Geometry::polygon(rect(10.0, -20.0, 60.0, 0.0), vec![])),
                    &[("name", "B".into())],
                )
                .unwrap();
            // Sliver: 10..12 x 0..4 (area 8). Edge with A: x=10,y0..4 -> len 4.
            //                                  Edge with B: y=0,x10..12 -> len 2.
            layer
                .add_feature(
                    Some(Geometry::polygon(rect(10.0, 0.0, 12.0, 4.0), vec![])),
                    &[("name", "S".into())],
                )
                .unwrap();
            let id = memory_store::put_vector(layer);
            memory_store::make_vector_memory_path(&id)
        };

        let host_name = |layer: &Layer, other: &str| -> f64 {
            let f = layer
                .features
                .iter()
                .find(|f| f.get(&layer.schema, "name").unwrap() == &FieldValue::Text(other.into()))
                .unwrap();
            area_of(f.geometry.as_ref().unwrap())
        };

        let (_, la) =
            run_tool(json!({ "input": build2(), "max_area": 10.0, "strategy": "longest_border" }));
        // A (longer border) grew by 8 -> 108; B unchanged at 1000.
        assert!(
            (host_name(&la, "A") - 108.0).abs() < 1e-6,
            "longest_border should merge into A"
        );
        assert!((host_name(&la, "B") - 1000.0).abs() < 1e-6);

        let (_, lb) =
            run_tool(json!({ "input": build2(), "max_area": 10.0, "strategy": "largest_area" }));
        // B (larger area) grew by 8 -> 1008; A unchanged at 100.
        assert!(
            (host_name(&lb, "B") - 1008.0).abs() < 1e-6,
            "largest_area should merge into B"
        );
        assert!((host_name(&lb, "A") - 100.0).abs() < 1e-6);
    }

    /// A chain of two adjacent slivers between a host and open space: the sliver
    /// touching the host merges first, then the second sliver merges into the
    /// grown host on the next pass.
    #[test]
    fn absorbs_a_chain_of_slivers_iteratively() {
        let mut layer = Layer::new("parcels");
        // Host: 0..10 x 0..10 (area 100).
        layer
            .add_feature(
                Some(Geometry::polygon(rect(0.0, 0.0, 10.0, 10.0), vec![])),
                &[],
            )
            .unwrap();
        // Sliver 1: 10..12 x 0..10 (area 20) touches host at x=10.
        layer
            .add_feature(
                Some(Geometry::polygon(rect(10.0, 0.0, 12.0, 10.0), vec![])),
                &[],
            )
            .unwrap();
        // Sliver 2: 12..13 x 0..10 (area 10) touches sliver 1 at x=12 only.
        layer
            .add_feature(
                Some(Geometry::polygon(rect(12.0, 0.0, 13.0, 10.0), vec![])),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input, "max_area": 25.0 }));
        assert_eq!(out.outputs["eliminated_count"], json!(2));
        assert_eq!(out.outputs["feature_count"], json!(1));
        assert_eq!(out.outputs["unmerged_count"], json!(0));
        assert!((area_of(layer.features[0].geometry.as_ref().unwrap()) - 130.0).abs() < 1e-6);
    }

    /// `where` selects by attribute; `exclude` protects a matching feature.
    #[test]
    fn where_selects_and_exclude_protects() {
        let make = || {
            let mut layer = Layer::new("landcover");
            layer.add_field(FieldDef::new("class", FieldType::Integer));
            // Host class 1 (big).
            layer
                .add_feature(
                    Some(Geometry::polygon(rect(0.0, 0.0, 10.0, 10.0), vec![])),
                    &[("class", 1i64.into())],
                )
                .unwrap();
            // class 0 sliver touching host — should be eliminated by `where`.
            layer
                .add_feature(
                    Some(Geometry::polygon(rect(10.0, 0.0, 11.0, 10.0), vec![])),
                    &[("class", 0i64.into())],
                )
                .unwrap();
            let id = memory_store::put_vector(layer);
            memory_store::make_vector_memory_path(&id)
        };

        let (out, _) = run_tool(json!({ "input": make(), "where": "class = 0" }));
        assert_eq!(out.outputs["eliminated_count"], json!(1));

        // Excluding class 0 keeps it.
        let (out2, _) =
            run_tool(json!({ "input": make(), "where": "class = 0", "exclude": "class = 0" }));
        assert_eq!(out2.outputs["eliminated_count"], json!(0));
        assert_eq!(out2.outputs["feature_count"], json!(2));
    }

    /// A sliver with no polygon neighbor cannot be merged and is retained.
    #[test]
    fn isolated_sliver_is_retained() {
        let mut layer = Layer::new("parcels");
        layer
            .add_feature(
                Some(Geometry::polygon(rect(0.0, 0.0, 10.0, 10.0), vec![])),
                &[],
            )
            .unwrap();
        // Detached sliver far away.
        layer
            .add_feature(
                Some(Geometry::polygon(rect(100.0, 100.0, 101.0, 101.0), vec![])),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, _) = run_tool(json!({ "input": input, "max_area": 5.0 }));
        assert_eq!(out.outputs["eliminated_count"], json!(0));
        assert_eq!(out.outputs["unmerged_count"], json!(1));
        assert_eq!(out.outputs["feature_count"], json!(2));
    }

    #[test]
    fn passes_non_polygons_through() {
        let mut layer = Layer::new("mixed");
        layer
            .add_feature(Some(Geometry::point(1.0, 2.0)), &[])
            .unwrap();
        layer
            .add_feature(
                Some(Geometry::polygon(rect(0.0, 0.0, 10.0, 10.0), vec![])),
                &[],
            )
            .unwrap();
        layer
            .add_feature(
                Some(Geometry::polygon(rect(10.0, 0.0, 11.0, 10.0), vec![])),
                &[],
            )
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input, "max_area": 15.0 }));
        // Only the sliver polygon is eliminated; the point survives.
        assert_eq!(out.outputs["eliminated_count"], json!(1));
        assert!(layer
            .features
            .iter()
            .any(|f| matches!(f.geometry, Some(Geometry::Point(_)))));
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = EliminatePolygonsTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input must fail");
        assert!(
            bad(json!({ "input": "x.geojson" })).is_err(),
            "no selector must fail"
        );
        assert!(bad(json!({ "input": "x.geojson", "max_area": 0 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "max_area": 5, "strategy": "bogus" })).is_err());
        assert!(
            bad(json!({ "input": "x.geojson", "where": "class" })).is_err(),
            "malformed condition"
        );
        assert!(bad(json!({ "input": "x.geojson", "where": "class = 0" })).is_ok());
        assert!(
            bad(json!({ "input": "x.geojson", "max_area": "5.0" })).is_ok(),
            "numeric strings ok"
        );
    }
}
