//! GeoLibre tool: trim short dangling line segments (overshoots and spurs).
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Trim Line* (Editing toolset): remove
//! *dangling* line segments — those that hang free at one or both ends — when
//! they are shorter than a `dangle_length` tolerance. These short dangles are
//! the overshoots and spurs left behind by raster-to-vector line extraction,
//! topology import, or careless digitizing; deleting them is a standard cleanup
//! before a network is usable.
//!
//! The bundled Whitebox tools do not cover this case: `fix_dangling_arcs`
//! *snaps* undershoot/overshoot dangles to the nearest arc rather than deleting
//! them; `prune_vector_streams` / `remove_short_streams` are stream-network
//! specific; `remove_spurs` operates on rasters. This tool works on any
//! vectorized line layer.
//!
//! ## How it decides
//!
//! Endpoints are snapped onto a grid (`snap_tolerance`) to build a node graph;
//! the *degree* of a node is how many line-ends touch it. For each segment we
//! count its **free ends** — endpoints whose node has degree 1 (nothing else
//! meets there):
//!
//! - **0 free ends** — an interior segment wired into the network at both ends;
//!   never removed, however short (it is not a dangle).
//! - **1 free end** — a true dangle (overshoot / spur) attached to the network
//!   at one end; removed when its length < `dangle_length`.
//! - **2 free ends** — an isolated segment touching nothing. Removed when short
//!   *only* if `keep_short` is false; by default (`keep_short` = true) these are
//!   preserved, matching ArcGIS's KEEP_SHORT default for features unconnected at
//!   either end.
//!
//! Decisions use the original topology in a single pass (they are not
//! re-evaluated after removals cascade), matching ArcGIS *Trim Line*. Non-line
//! features pass through untouched; `MultiLineString` inputs are exploded into
//! their parts, each judged independently.

use std::collections::{BTreeMap, HashMap};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Coord, Feature, Geometry, GeometryType};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct TrimLineTool;

impl Tool for TrimLineTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "trim_line",
            display_name: "Trim Line",
            summary: "Delete short dangling line segments (overshoots and spurs) shorter than a length tolerance, keeping true intersections and long lines.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input line vector file path, format auto-detected (or in-memory handle).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output vector path (driver from its extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "dangle_length",
                    description: "Length tolerance in CRS units. A dangling segment shorter than this is trimmed; longer dangles are kept. Must be positive.",
                    required: true,
                },
                ToolParamSpec {
                    name: "keep_short",
                    description: "When true (default), short segments that are unconnected at both ends (isolated) are kept; when false they are also deleted. Segments dangling at exactly one end are always trimmed when short.",
                    required: false,
                },
                ToolParamSpec {
                    name: "snap_tolerance",
                    description: "Endpoints within this distance (CRS units) are treated as the same node. Default 1e-6.",
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

        // Explode features into individual polylines, remembering their source
        // feature (for attributes). Non-line features are collected for
        // pass-through.
        let mut lines: Vec<LineRef> = Vec::new();
        let mut passthrough: Vec<Feature> = Vec::new();
        for (fidx, feature) in layer.features.iter().enumerate() {
            match &feature.geometry {
                Some(Geometry::LineString(cs)) => {
                    if cs.len() >= 2 {
                        lines.push(LineRef {
                            coords: cs.clone(),
                            feat: fidx,
                        });
                    }
                }
                Some(Geometry::MultiLineString(parts)) => {
                    for cs in parts {
                        if cs.len() >= 2 {
                            lines.push(LineRef {
                                coords: cs.clone(),
                                feat: fidx,
                            });
                        }
                    }
                }
                _ => passthrough.push(feature.clone()),
            }
        }
        let input_line_count = lines.len();

        // Assign a node id to every distinct endpoint (snapped) and tally the
        // degree of each node (number of incident line-ends).
        let mut nodes = NodeIndex::new(prm.snap_tolerance);
        let mut ends: Vec<(usize, usize)> = Vec::with_capacity(lines.len()); // (start_node, end_node)
        for l in &lines {
            let a = nodes.id_of(&l.coords[0]);
            let b = nodes.id_of(l.coords.last().unwrap());
            ends.push((a, b));
        }
        let mut degree: HashMap<usize, usize> = HashMap::new();
        for &(a, b) in &ends {
            *degree.entry(a).or_insert(0) += 1;
            *degree.entry(b).or_insert(0) += 1;
        }

        ctx.progress.info(&format!(
            "{input_line_count} line segment(s), {} node(s)",
            degree.len()
        ));

        // Decide, per segment, whether to keep it.
        let mut removed_count = 0usize;
        let mut isolated_removed = 0usize;
        let mut keep = vec![true; lines.len()];
        for (li, l) in lines.iter().enumerate() {
            let (a, b) = ends[li];
            let free_ends = (degree.get(&a).copied().unwrap_or(0) == 1) as u8
                + (degree.get(&b).copied().unwrap_or(0) == 1) as u8;
            if free_ends == 0 {
                continue; // interior segment, never a dangle
            }
            let len = polyline_length(&l.coords);
            if len >= prm.dangle_length {
                continue; // dangle, but long enough to keep
            }
            if free_ends == 2 && prm.keep_short {
                continue; // isolated short line preserved by default
            }
            keep[li] = false;
            removed_count += 1;
            if free_ends == 2 {
                isolated_removed += 1;
            }
        }

        // Build the output layer: surviving line segments (each as a LineString
        // carrying its source attributes), then pass-through non-line features.
        let mut out_features: Vec<Feature> = Vec::with_capacity(lines.len() + passthrough.len());
        for (li, l) in lines.iter().enumerate() {
            if !keep[li] {
                continue;
            }
            let mut f = layer.features[l.feat].clone();
            f.geometry = Some(Geometry::LineString(l.coords.clone()));
            f.fid = out_features.len() as u64;
            out_features.push(f);
        }
        for mut f in passthrough {
            f.fid = out_features.len() as u64;
            out_features.push(f);
        }

        ctx.progress.info(&format!(
            "trimmed {removed_count} short dangle(s) ({isolated_removed} isolated); {} line(s) remain",
            input_line_count - removed_count
        ));

        let mut out_layer = wbvector::Layer::new(layer.name);
        out_layer.schema = layer.schema.clone();
        out_layer.crs = layer.crs;
        out_layer.geom_type = Some(GeometryType::LineString);
        out_layer.features = out_features;

        let feature_count = out_layer.len();
        let out_path = write_or_store_layer(out_layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_line_count".to_string(), json!(input_line_count));
        outputs.insert("feature_count".to_string(), json!(feature_count));
        outputs.insert("removed_count".to_string(), json!(removed_count));
        outputs.insert("isolated_removed".to_string(), json!(isolated_removed));
        Ok(ToolRunResult { outputs })
    }
}

/// A single polyline together with its source feature index.
struct LineRef {
    coords: Vec<Coord>,
    feat: usize,
}

/// Total planar (Euclidean) length of a polyline.
fn polyline_length(coords: &[Coord]) -> f64 {
    coords
        .windows(2)
        .map(|w| {
            let dx = w[1].x - w[0].x;
            let dy = w[1].y - w[0].y;
            (dx * dx + dy * dy).sqrt()
        })
        .sum()
}

// ── Node index (endpoint snapping) ──────────────────────────────────────────

/// Maps snapped endpoint coordinates to integer node ids. Coordinates are
/// quantized onto a grid of cell `tol`, and each new grid cell (and its already
/// seen 8-neighbours, so points straddling a cell edge still merge) yields one
/// node id.
struct NodeIndex {
    tol: f64,
    map: HashMap<(i64, i64), usize>,
    next: usize,
}

impl NodeIndex {
    fn new(tol: f64) -> Self {
        Self {
            tol: tol.max(f64::MIN_POSITIVE),
            map: HashMap::new(),
            next: 0,
        }
    }

    fn cell(&self, c: &Coord) -> (i64, i64) {
        (
            (c.x / self.tol).round() as i64,
            (c.y / self.tol).round() as i64,
        )
    }

    fn id_of(&mut self, c: &Coord) -> usize {
        let (cx, cy) = self.cell(c);
        // Reuse a node from this cell or any of the 8 neighbours (handles points
        // that snap to adjacent grid cells).
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(&id) = self.map.get(&(cx + dx, cy + dy)) {
                    self.map.insert((cx, cy), id);
                    return id;
                }
            }
        }
        let id = self.next;
        self.next += 1;
        self.map.insert((cx, cy), id);
        id
    }
}

// ── Parameters ──────────────────────────────────────────────────────────────

struct Params {
    dangle_length: f64,
    keep_short: bool,
    snap_tolerance: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let dangle_length = parse_optional_f64(args, "dangle_length")?.ok_or_else(|| {
        ToolError::Validation("missing required numeric parameter 'dangle_length'".to_string())
    })?;
    if !(dangle_length > 0.0 && dangle_length.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'dangle_length' must be a positive number".to_string(),
        ));
    }
    let keep_short = parse_optional_bool(args, "keep_short")?.unwrap_or(true);
    let snap_tolerance = parse_optional_f64(args, "snap_tolerance")?.unwrap_or(1e-6);
    if !(snap_tolerance > 0.0 && snap_tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "parameter 'snap_tolerance' must be a positive number".to_string(),
        ));
    }
    Ok(Params {
        dangle_length,
        keep_short,
        snap_tolerance,
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

fn parse_optional_bool(args: &ToolArgs, key: &str) -> Result<Option<bool>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) if s.trim().is_empty() => Ok(None),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "true" | "1" | "yes" | "y" => Ok(Some(true)),
            "false" | "0" | "no" | "n" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{key}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{key}' must be a boolean"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{memory_store, Layer};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn line(pts: &[(f64, f64)]) -> Geometry {
        Geometry::LineString(pts.iter().map(|&(x, y)| Coord::xy(x, y)).collect())
    }

    fn run_tool(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = TrimLineTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// A long backbone with one short spur hanging off a mid-node: the spur is
    /// trimmed, the backbone survives.
    #[test]
    fn trims_a_short_spur() {
        let mut layer = Layer::new("net");
        // Backbone: (0,0)->(10,0)->(20,0), two segments meeting at (10,0).
        layer
            .add_feature(Some(line(&[(0.0, 0.0), (10.0, 0.0)])), &[])
            .unwrap();
        layer
            .add_feature(Some(line(&[(10.0, 0.0), (20.0, 0.0)])), &[])
            .unwrap();
        // Short spur hanging off the (10,0) node, length 1.
        layer
            .add_feature(Some(line(&[(10.0, 0.0), (10.0, 1.0)])), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input, "dangle_length": 2.0 }));
        assert_eq!(out.outputs["removed_count"], json!(1));
        assert_eq!(out.outputs["feature_count"], json!(2));
        // No surviving segment should be the vertical spur.
        assert!(layer.features.iter().all(|f| {
            !matches!(&f.geometry, Some(Geometry::LineString(cs)) if cs.last().unwrap().y == 1.0)
        }));
    }

    /// A long dangle (over the tolerance) is kept.
    #[test]
    fn keeps_long_dangle() {
        let mut layer = Layer::new("net");
        layer
            .add_feature(Some(line(&[(0.0, 0.0), (10.0, 0.0)])), &[])
            .unwrap();
        layer
            .add_feature(Some(line(&[(10.0, 0.0), (20.0, 0.0)])), &[])
            .unwrap();
        // Dangle length 5 > tolerance 2 -> kept.
        layer
            .add_feature(Some(line(&[(10.0, 0.0), (10.0, 5.0)])), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, _) = run_tool(json!({ "input": input, "dangle_length": 2.0 }));
        assert_eq!(out.outputs["removed_count"], json!(0));
        assert_eq!(out.outputs["feature_count"], json!(3));
    }

    /// An interior segment shorter than the tolerance is NOT removed (both ends
    /// connected -> not a dangle).
    #[test]
    fn keeps_short_interior_segment() {
        let mut layer = Layer::new("net");
        // A short middle segment (length 1) wired between two long arms.
        layer
            .add_feature(Some(line(&[(0.0, 0.0), (10.0, 0.0)])), &[])
            .unwrap();
        layer
            .add_feature(Some(line(&[(10.0, 0.0), (11.0, 0.0)])), &[]) // short interior
            .unwrap();
        layer
            .add_feature(Some(line(&[(11.0, 0.0), (21.0, 0.0)])), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, _) = run_tool(json!({ "input": input, "dangle_length": 5.0 }));
        assert_eq!(out.outputs["removed_count"], json!(0));
        assert_eq!(out.outputs["feature_count"], json!(3));
    }

    /// An isolated short line: kept by default, deleted when keep_short=false.
    #[test]
    fn isolated_short_line_honours_keep_short() {
        let build = || {
            let mut layer = Layer::new("net");
            // Connected backbone.
            layer
                .add_feature(Some(line(&[(0.0, 0.0), (10.0, 0.0)])), &[])
                .unwrap();
            layer
                .add_feature(Some(line(&[(10.0, 0.0), (20.0, 0.0)])), &[])
                .unwrap();
            // Isolated short line, both ends free, length 1.
            layer
                .add_feature(Some(line(&[(0.0, 50.0), (1.0, 50.0)])), &[])
                .unwrap();
            memory_store::make_vector_memory_path(&memory_store::put_vector(layer))
        };

        // Default keep_short=true -> isolated line kept.
        let (kept, _) = run_tool(json!({ "input": build(), "dangle_length": 2.0 }));
        assert_eq!(kept.outputs["removed_count"], json!(0));
        assert_eq!(kept.outputs["feature_count"], json!(3));

        // keep_short=false -> isolated short line deleted.
        let (del, _) =
            run_tool(json!({ "input": build(), "dangle_length": 2.0, "keep_short": false }));
        assert_eq!(del.outputs["removed_count"], json!(1));
        assert_eq!(del.outputs["isolated_removed"], json!(1));
        assert_eq!(del.outputs["feature_count"], json!(2));
    }

    /// Non-line features pass through untouched.
    #[test]
    fn passes_points_through() {
        let mut layer = Layer::new("mixed");
        layer
            .add_feature(Some(Geometry::point(5.0, 5.0)), &[])
            .unwrap();
        layer
            .add_feature(Some(line(&[(0.0, 0.0), (10.0, 0.0)])), &[])
            .unwrap();
        layer
            .add_feature(Some(line(&[(10.0, 0.0), (20.0, 0.0)])), &[])
            .unwrap();
        // Short spur to be trimmed.
        layer
            .add_feature(Some(line(&[(10.0, 0.0), (10.0, 1.0)])), &[])
            .unwrap();
        let id = memory_store::put_vector(layer);
        let input = memory_store::make_vector_memory_path(&id);

        let (out, layer) = run_tool(json!({ "input": input, "dangle_length": 2.0 }));
        // spur removed; point + 2 backbone segments remain.
        assert_eq!(out.outputs["removed_count"], json!(1));
        assert_eq!(out.outputs["feature_count"], json!(3));
        assert!(layer
            .features
            .iter()
            .any(|f| matches!(f.geometry, Some(Geometry::Point(_)))));
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = TrimLineTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err()); // no input
        assert!(bad(json!({ "input": "x.geojson" })).is_err()); // no dangle_length
        assert!(bad(json!({ "input": "x.geojson", "dangle_length": 0 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "dangle_length": -1 })).is_err());
        assert!(bad(json!({ "input": "x.geojson", "dangle_length": 5.0 })).is_ok());
        assert!(bad(json!({ "input": "x.geojson", "dangle_length": "5.0" })).is_ok());
        assert!(
            bad(json!({ "input": "x.geojson", "dangle_length": 5.0, "keep_short": "false" }))
                .is_ok()
        );
        assert!(
            bad(json!({ "input": "x.geojson", "dangle_length": 5.0, "snap_tolerance": 0 }))
                .is_err()
        );
    }
}
