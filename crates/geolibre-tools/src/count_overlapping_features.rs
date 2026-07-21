//! GeoLibre tool: flatten overlapping polygons into counted disjoint regions.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Count Overlapping Features*
//! (Analysis). The bundled `overlaps` is a boolean predicate test and `count` is
//! a raster/statistics op — nothing flattens a set of overlapping polygons into
//! the disjoint regions, each attributed with how many inputs cover it, that
//! buffer-coverage analysis, service-area redundancy, and imagery-footprint
//! dedup need.
//!
//! The input polygons are overlaid incrementally into a set of mutually
//! disjoint regions, each carrying the set of input feature ids that cover it
//! (so `count` is the coverage depth). For each input polygon `P` and each
//! existing region `R`:
//!
//! * `R ∩ P` becomes a region with `R`'s ids plus `P`;
//! * `R \ P` keeps `R`'s ids;
//! * whatever of `P` was not already covered by any region becomes a new
//!   depth-1 region.
//!
//! Bounding-box prefiltering skips region/polygon pairs that cannot overlap.
//! `min_count` drops regions below a coverage depth (e.g. 2 to keep only
//! overlaps); `report_ids` writes a region→feature-id CSV. Uses the same `geo`
//! `BooleanOps` machinery as `tabulate_intersection`.

use std::collections::BTreeMap;

use geo::{Area, BooleanOps, BoundingRect, Coord as GeoCoord, LineString, MultiPolygon, Polygon};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, FieldDef, FieldType, FieldValue, Geometry, GeometryType, Layer, Ring};

use crate::common::write_text_output;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Regions smaller than this fraction of the mean input area are dropped as
/// numerical slivers.
const SLIVER_AREA_EPS: f64 = 1e-9;

pub struct CountOverlappingFeaturesTool;

impl Tool for CountOverlappingFeaturesTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "count_overlapping_features",
            display_name: "Count Overlapping Features",
            summary: "Flatten overlapping polygons into disjoint regions, each attributed with the number (and optionally the ids) of input features covering it, like ArcGIS Count Overlapping Features.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon vector layer (features may overlap).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon vector path (driver from extension) with a 'count' field. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "min_count",
                    description: "Keep only regions covered by at least this many features (default 1; use 2 for overlaps only).",
                    required: false,
                },
                ToolParamSpec {
                    name: "id_field",
                    description: "Optional field identifying each input feature (used in the ids column / report). Default: the feature index.",
                    required: false,
                },
                ToolParamSpec {
                    name: "report_ids",
                    description: "Optional CSV path for a region_id,feature_id join table (one row per covering feature).",
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
        let id_idx = match &prm.id_field {
            Some(f) => Some(
                layer
                    .schema
                    .field_index(f)
                    .ok_or_else(|| ToolError::Validation(format!("id_field '{f}' not found")))?,
            ),
            None => None,
        };

        // Collect input polygons as geo MultiPolygons with their ids.
        let mut inputs: Vec<(String, MultiPolygon)> = Vec::new();
        for (fidx, feature) in layer.features.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let Some(mp) = to_multipolygon(geom) else {
                continue;
            };
            if mp.unsigned_area() <= 0.0 {
                continue;
            }
            let id = match id_idx {
                Some(i) => field_key(&feature.attributes[i]),
                None => fidx.to_string(),
            };
            inputs.push((id, mp));
        }
        if inputs.is_empty() {
            return Err(ToolError::Execution(
                "no polygon features in input".to_string(),
            ));
        }

        ctx.progress
            .info(&format!("overlaying {} polygon(s)", inputs.len()));

        // ── Incremental disjoint-region overlay ───────────────────────────────
        let mut regions: Vec<Region> = Vec::new();
        for (id, poly) in &inputs {
            let p_bbox = bbox(poly);
            let mut remaining = poly.clone();
            let mut next: Vec<Region> = Vec::with_capacity(regions.len() + 1);
            for r in regions.drain(..) {
                if !bbox_overlap(&r.bbox, &p_bbox) {
                    next.push(r);
                    continue;
                }
                let inter = r.geom.intersection(poly);
                if inter.unsigned_area() > SLIVER_AREA_EPS {
                    let mut ids = r.ids.clone();
                    ids.push(id.clone());
                    next.push(Region::new(inter, ids));
                    let diff = r.geom.difference(poly);
                    if diff.unsigned_area() > SLIVER_AREA_EPS {
                        next.push(Region::new(diff, r.ids));
                    }
                    remaining = remaining.difference(&r.geom);
                } else {
                    next.push(r);
                }
            }
            if remaining.unsigned_area() > SLIVER_AREA_EPS {
                next.push(Region::new(remaining, vec![id.clone()]));
            }
            regions = next;
        }

        ctx.progress
            .info(&format!("{} disjoint region(s)", regions.len()));

        // ── Build the output layer ────────────────────────────────────────────
        let mut out = Layer::new("overlap_counts").with_geom_type(GeometryType::MultiPolygon);
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new("region_id", FieldType::Integer));
        out.add_field(FieldDef::new("count", FieldType::Integer));
        out.add_field(FieldDef::new("ids", FieldType::Text));

        let mut csv = String::from("region_id,feature_id\n");
        let mut kept = 0usize;
        let mut max_count = 0usize;
        for region in &regions {
            let count = region.ids.len();
            if count < prm.min_count {
                continue;
            }
            max_count = max_count.max(count);
            let region_id = kept as i64;
            let geom = multipolygon_to_geometry(&region.geom);
            let ids_joined = region.ids.join(";");
            out.add_feature(
                Some(geom),
                &[
                    ("region_id", region_id.into()),
                    ("count", (count as i64).into()),
                    ("ids", ids_joined.into()),
                ],
            )
            .map_err(|e| ToolError::Execution(format!("failed writing region: {e}")))?;
            if prm.report_ids.is_some() {
                for fid in &region.ids {
                    csv.push_str(&format!("{region_id},{fid}\n"));
                }
            }
            kept += 1;
        }

        if let Some(path) = &prm.report_ids {
            write_text_output(&csv, path)?;
        }

        let out_path = write_or_store_layer(out, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_count".to_string(), json!(inputs.len()));
        outputs.insert("region_count".to_string(), json!(kept));
        outputs.insert("max_overlap".to_string(), json!(max_count));
        if let Some(path) = &prm.report_ids {
            outputs.insert("report_ids".to_string(), json!(path));
        }
        Ok(ToolRunResult { outputs })
    }
}

// ── Regions ──────────────────────────────────────────────────────────────────

struct Region {
    geom: MultiPolygon,
    ids: Vec<String>,
    bbox: [f64; 4],
}

impl Region {
    fn new(geom: MultiPolygon, ids: Vec<String>) -> Region {
        let bbox = bbox(&geom);
        Region { geom, ids, bbox }
    }
}

fn bbox(mp: &MultiPolygon) -> [f64; 4] {
    match mp.bounding_rect() {
        Some(r) => [r.min().x, r.min().y, r.max().x, r.max().y],
        None => [
            f64::INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
        ],
    }
}

fn bbox_overlap(a: &[f64; 4], b: &[f64; 4]) -> bool {
    a[0] <= b[2] && b[0] <= a[2] && a[1] <= b[3] && b[1] <= a[3]
}

fn field_key(fv: &FieldValue) -> String {
    if let Some(i) = fv.as_i64() {
        i.to_string()
    } else if let Some(f) = fv.as_f64() {
        format!("{f}")
    } else {
        fv.as_str().unwrap_or("").to_string()
    }
}

// ── geo <-> wbvector conversion (shared shape with tabulate_intersection) ─────

fn to_multipolygon(geom: &Geometry) -> Option<MultiPolygon> {
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => Some(MultiPolygon(vec![rings_to_polygon(exterior, interiors)])),
        Geometry::MultiPolygon(parts) => Some(MultiPolygon(
            parts.iter().map(|(e, i)| rings_to_polygon(e, i)).collect(),
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

fn multipolygon_to_geometry(mp: &MultiPolygon) -> Geometry {
    Geometry::MultiPolygon(mp.0.iter().map(polygon_to_rings).collect())
}

fn polygon_to_rings(poly: &Polygon) -> (Ring, Vec<Ring>) {
    (
        linestring_to_ring(poly.exterior()),
        poly.interiors().iter().map(linestring_to_ring).collect(),
    )
}

fn linestring_to_ring(ls: &LineString) -> Ring {
    let mut coords: Vec<Coord> = ls.0.iter().map(|c| Coord::xy(c.x, c.y)).collect();
    if coords.len() >= 2 && coords.first() == coords.last() {
        coords.pop();
    }
    Ring::new(coords)
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    min_count: usize,
    id_field: Option<String>,
    report_ids: Option<String>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let min_count = match args.get("min_count") {
        None | Some(Value::Null) => 1,
        Some(Value::Number(n)) => n.as_u64().unwrap_or(1).max(1) as usize,
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => s
            .trim()
            .parse::<usize>()
            .map_err(|_| ToolError::Validation("'min_count' must be an integer".to_string()))?
            .max(1),
        Some(_) => {
            return Err(ToolError::Validation(
                "'min_count' must be a number".to_string(),
            ))
        }
    };
    let id_field = parse_optional_str(args, "id_field")?.map(str::to_string);
    let report_ids = parse_optional_str(args, "report_ids")?.map(str::to_string);
    Ok(Params {
        min_count,
        id_field,
        report_ids,
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

    fn square(x0: f64, y0: f64, s: f64) -> Geometry {
        Geometry::polygon(
            vec![
                Coord::xy(x0, y0),
                Coord::xy(x0 + s, y0),
                Coord::xy(x0 + s, y0 + s),
                Coord::xy(x0, y0 + s),
            ],
            vec![],
        )
    }

    fn layer_of(geoms: Vec<Geometry>) -> String {
        let mut l = Layer::new("polys")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        for g in geoms {
            l.add_feature(Some(g), &[]).unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CountOverlappingFeaturesTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn area(geom: &Geometry) -> f64 {
        to_multipolygon(geom)
            .map(|m| m.unsigned_area())
            .unwrap_or(0.0)
    }

    /// Two overlapping unit squares -> three regions: two depth-1, one depth-2.
    #[test]
    fn two_overlapping_squares() {
        // A: [0,10]^2, B: [5,15]x[0,10]; overlap is [5,10]x[0,10] (area 50).
        let a = square(0.0, 0.0, 10.0);
        let b = square(5.0, 0.0, 10.0);
        let input = layer_of(vec![a, b]);
        let (out, layer) = run(json!({ "input": input }));
        assert_eq!(out.outputs["max_overlap"], json!(2));

        let cidx = layer.schema.field_index("count").unwrap();
        let mut total_overlap = 0.0;
        let mut count2 = 0;
        for f in layer.iter() {
            let c = f.attributes[cidx].as_i64().unwrap();
            if c == 2 {
                count2 += 1;
                total_overlap += area(f.geometry.as_ref().unwrap());
            }
        }
        assert_eq!(count2, 1, "exactly one depth-2 region");
        assert!(
            (total_overlap - 50.0).abs() < 1e-6,
            "overlap area {total_overlap} != 50"
        );
    }

    /// min_count=2 keeps only the overlap regions.
    #[test]
    fn min_count_filters_to_overlaps() {
        let a = square(0.0, 0.0, 10.0);
        let b = square(5.0, 0.0, 10.0);
        let input = layer_of(vec![a, b]);
        let (out, layer) = run(json!({ "input": input, "min_count": 2 }));
        assert_eq!(out.outputs["region_count"], json!(1));
        let cidx = layer.schema.field_index("count").unwrap();
        for f in layer.iter() {
            assert!(f.attributes[cidx].as_i64().unwrap() >= 2);
        }
    }

    /// Three mutually overlapping squares reach depth 3 in the common core, and
    /// the total flattened area equals the area of their union.
    #[test]
    fn triple_overlap_depth_and_union_area() {
        let a = square(0.0, 0.0, 10.0);
        let b = square(4.0, 0.0, 10.0);
        let c = square(2.0, 3.0, 10.0);
        let union_area = {
            let ma = to_multipolygon(&a).unwrap();
            let mb = to_multipolygon(&b).unwrap();
            let mc = to_multipolygon(&c).unwrap();
            ma.union(&mb).union(&mc).unsigned_area()
        };
        let input = layer_of(vec![a, b, c]);
        let (out, layer) = run(json!({ "input": input }));
        assert_eq!(
            out.outputs["max_overlap"],
            json!(3),
            "core should be covered 3x"
        );
        let flattened: f64 = layer
            .iter()
            .map(|f| area(f.geometry.as_ref().unwrap()))
            .sum();
        assert!(
            (flattened - union_area).abs() < 1e-4,
            "flattened area {flattened} != union {union_area}"
        );
        // Regions are disjoint: sum of area == union (no double counting).
    }

    /// Non-overlapping squares stay depth 1.
    #[test]
    fn disjoint_squares_stay_depth_one() {
        let a = square(0.0, 0.0, 5.0);
        let b = square(100.0, 100.0, 5.0);
        let input = layer_of(vec![a, b]);
        let (out, _l) = run(json!({ "input": input }));
        assert_eq!(out.outputs["max_overlap"], json!(1));
        assert_eq!(out.outputs["region_count"], json!(2));
    }

    #[test]
    fn rejects_missing_input() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(CountOverlappingFeaturesTool.validate(&args).is_err());
    }
}
