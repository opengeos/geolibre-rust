//! GeoLibre tool: populate directional neighbour-sheet fields on a map-series grid.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Calculate Adjacent Fields*
//! (Cartography). GeoLibre's `grid_index_features` / `strip_map_index_features`
//! build map-book index grids and `polygon_neighbors` reports adjacency counts,
//! but nothing writes the eight "continued on sheet …" margin labels that a map
//! book needs. This tool fills that gap: for every grid tile it locates the
//! neighbouring sheet in each of the eight compass directions and writes that
//! neighbour's page name into a directional field (`N`, `NE`, `E`, `SE`, `S`,
//! `SW`, `W`, `NW`).
//!
//! **How neighbours are found.** Every polygon boundary is decomposed into
//! undirected edges (endpoint-keyed, optionally snapped to a grid so
//! near-coincident borders match) exactly as `polygon_neighbors` does. Two tiles
//! that share an edge are *rook* neighbours and get a cardinal slot (N/E/S/W);
//! two tiles that touch only at a corner (a shared vertex, no shared edge) are
//! *diagonal* neighbours and get a corner slot (NE/SE/SW/NW). The specific slot
//! comes from the sign / dominant axis of the vector between the two tiles'
//! bounding-box centres, which is exact for the axis-aligned rectangles a map
//! series is made of and robust to any tile aspect ratio. If more than one
//! neighbour lands in the same slot (irregular grids) the nearest one wins.
//!
//! The output layer preserves the input geometry and attributes and appends the
//! directional text fields (four cardinal ones only when `include_diagonal` is
//! false).

use std::collections::{BTreeMap, HashMap, HashSet};

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Geometry, Ring};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct CalculateAdjacentFieldsTool;

/// Directional slots, indexed 0..8 (N, NE, E, SE, S, SW, W, NW).
const DIR_NAMES: [&str; 8] = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"];
const CARDINAL: [usize; 4] = [0, 2, 4, 6];
const DIAGONAL: [usize; 4] = [1, 3, 5, 7];

impl Tool for CalculateAdjacentFieldsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "calculate_adjacent_fields",
            display_name: "Calculate Adjacent Fields",
            summary: "Populate eight directional neighbour-sheet fields (N, NE, E, SE, S, SW, W, NW) on each tile of a map-series index grid with the adjoining page's name, for 'continued on sheet …' margin labels — like ArcGIS Calculate Adjacent Fields.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input polygon map-series index grid (e.g. output of grid_index_features).",
                    required: true,
                },
                ToolParamSpec {
                    name: "page_name_field",
                    description: "Field whose value labels each page (written into the neighbours' directional fields).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output polygon vector (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "include_diagonal",
                    description: "Also populate the four diagonal fields (NE, SE, SW, NW) from corner-touching tiles. Default true.",
                    required: false,
                },
                ToolParamSpec {
                    name: "snap_tolerance",
                    description: "Quantize vertices onto a grid of this size (CRS units) before matching borders. Default 0 (exact coordinates).",
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
        if parse_optional_str(args, "page_name_field")?.is_none() {
            return Err(ToolError::Validation(
                "missing required string parameter 'page_name_field'".to_string(),
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
        let page_field = parse_optional_str(args, "page_name_field")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'page_name_field'".to_string())
        })?;
        let prm = parse_params(args)?;

        let mut layer = load_input_layer(input)?;
        let name_idx = layer.schema.field_index(page_field).ok_or_else(|| {
            ToolError::Validation(format!("page_name_field '{page_field}' not found"))
        })?;

        let n = layer.features.len();

        // ── Per-feature page name + bounding-box centre ───────────────────────
        let mut names: Vec<String> = Vec::with_capacity(n);
        let mut centers: Vec<Option<(f64, f64)>> = Vec::with_capacity(n);
        for feature in &layer.features {
            names.push(field_string(&feature.attributes[name_idx]));
            centers.push(feature.geometry.as_ref().and_then(bbox_center));
        }

        // ── Edge → features and vertex → features maps (polygon_neighbors) ────
        let mut edge_feats: HashMap<(Key, Key), HashSet<usize>> = HashMap::new();
        let mut vert_feats: HashMap<Key, HashSet<usize>> = HashMap::new();
        let mut poly_count = 0usize;

        for (fidx, feature) in layer.features.iter().enumerate() {
            let Some(geom) = feature.geometry.as_ref() else {
                continue;
            };
            let rings = polygon_rings(geom, prm.snap_tolerance);
            if rings.is_empty() {
                continue;
            }
            poly_count += 1;
            for ring in &rings {
                let m = ring.len();
                for i in 0..m {
                    let a = ring[i];
                    let b = ring[(i + 1) % m];
                    vert_feats.entry(key(a)).or_default().insert(fidx);
                    if key(a) == key(b) {
                        continue;
                    }
                    edge_feats.entry(edge_key(a, b)).or_default().insert(fidx);
                }
            }
        }

        // ── Per-pair adjacency: edge-sharing (rook) vs node-only (diagonal) ───
        // shared_edge[(a,b)] = true when they share at least one full edge.
        let mut shared_edge: HashMap<(usize, usize), bool> = HashMap::new();
        for feats in edge_feats.values() {
            if feats.len() < 2 {
                continue;
            }
            let list: Vec<usize> = feats.iter().copied().collect();
            for i in 0..list.len() {
                for j in (i + 1)..list.len() {
                    shared_edge.insert(order(list[i], list[j]), true);
                }
            }
        }
        // node-only pairs: touch at a vertex but share no edge.
        let mut pairs: HashSet<(usize, usize)> = shared_edge.keys().copied().collect();
        for feats in vert_feats.values() {
            if feats.len() < 2 {
                continue;
            }
            let list: Vec<usize> = feats.iter().copied().collect();
            for i in 0..list.len() {
                for j in (i + 1)..list.len() {
                    pairs.insert(order(list[i], list[j]));
                }
            }
        }

        // ── Assign each neighbour to a directional slot (nearest wins) ────────
        // slot[fidx][dir] = (distance², neighbour idx)
        let mut slots: Vec<[Option<(f64, usize)>; 8]> = vec![[None; 8]; n];
        let assign = |src: usize, dst: usize, slots: &mut Vec<[Option<(f64, usize)>; 8]>| {
            let (Some((sx, sy)), Some((dx_, dy_))) = (centers[src], centers[dst]) else {
                return;
            };
            let dx = dx_ - sx;
            let dy = dy_ - sy;
            let d2 = dx * dx + dy * dy;
            let is_rook = shared_edge.contains_key(&order(src, dst));
            let Some(dir) = direction_slot(dx, dy, is_rook) else {
                return;
            };
            if !prm.include_diagonal && DIAGONAL.contains(&dir) {
                return;
            }
            match slots[src][dir] {
                Some((best, _)) if best <= d2 => {}
                _ => slots[src][dir] = Some((d2, dst)),
            }
        };
        for &(a, b) in &pairs {
            assign(a, b, &mut slots);
            assign(b, a, &mut slots);
        }

        // ── Append directional fields and fill them ───────────────────────────
        let active: Vec<usize> = if prm.include_diagonal {
            (0..8).collect()
        } else {
            CARDINAL.to_vec()
        };
        // Directional field names, disambiguated against existing schema names.
        let field_of: Vec<String> = active
            .iter()
            .map(|&d| unique_field_name(&layer, DIR_NAMES[d]))
            .collect();
        for name in &field_of {
            layer.add_field(FieldDef::new(name.clone(), FieldType::Text));
        }

        let mut populated = 0usize;
        for (fidx, feature) in layer.features.iter_mut().enumerate() {
            for &dir in &active {
                let value = match slots[fidx][dir] {
                    Some((_, nbr)) => {
                        populated += 1;
                        FieldValue::Text(names[nbr].clone())
                    }
                    None => FieldValue::Text(String::new()),
                };
                feature.attributes.push(value);
            }
        }

        ctx.progress.info(&format!(
            "{poly_count} page(s); {} neighbour pair(s); {populated} directional label(s) written",
            pairs.len()
        ));

        let out_path = write_or_store_layer(layer, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("page_count".to_string(), json!(poly_count));
        outputs.insert("neighbor_pairs".to_string(), json!(pairs.len()));
        outputs.insert("labels_written".to_string(), json!(populated));
        outputs.insert("direction_fields".to_string(), json!(field_of));
        Ok(ToolRunResult { outputs })
    }
}

// ── Direction classification ─────────────────────────────────────────────────

/// Map a centre-to-centre delta to a directional slot (0..8). Edge-sharing
/// (`rook`) pairs land in a cardinal slot by dominant axis; corner-only pairs
/// land in a diagonal slot by the sign quadrant. Returns None when the delta is
/// degenerate (zero or on an axis where a diagonal is expected).
fn direction_slot(dx: f64, dy: f64, rook: bool) -> Option<usize> {
    if rook {
        if dx == 0.0 && dy == 0.0 {
            return None;
        }
        if dy.abs() >= dx.abs() {
            Some(if dy > 0.0 { 0 } else { 4 }) // N / S
        } else {
            Some(if dx > 0.0 { 2 } else { 6 }) // E / W
        }
    } else {
        match (dx > 0.0, dy > 0.0, dx == 0.0, dy == 0.0) {
            (_, _, true, _) | (_, _, _, true) => None,
            (true, true, ..) => Some(1),   // NE
            (true, false, ..) => Some(3),  // SE
            (false, false, ..) => Some(5), // SW
            (false, true, ..) => Some(7),  // NW
        }
    }
}

// ── Field-name helpers ───────────────────────────────────────────────────────

/// Return `base`, or `base_2`, `base_3`, … if `base` already exists in the
/// schema, so appended direction fields never collide with input attributes.
fn unique_field_name(layer: &wbvector::Layer, base: &str) -> String {
    if layer.schema.field_index(base).is_none() {
        return base.to_string();
    }
    let mut i = 2;
    loop {
        let candidate = format!("{base}_{i}");
        if layer.schema.field_index(&candidate).is_none() {
            return candidate;
        }
        i += 1;
    }
}

// ── Geometry keys (shared with polygon_neighbors) ────────────────────────────

#[derive(Clone, Copy)]
struct P {
    x: f64,
    y: f64,
}

type Key = (u64, u64);

fn key(p: P) -> Key {
    (p.x.to_bits(), p.y.to_bits())
}

fn edge_key(a: P, b: P) -> (Key, Key) {
    let (ka, kb) = (key(a), key(b));
    if ka <= kb {
        (ka, kb)
    } else {
        (kb, ka)
    }
}

fn order(a: usize, b: usize) -> (usize, usize) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

fn canonical(x: f64, y: f64, snap: f64) -> P {
    if snap > 0.0 {
        P {
            x: (x / snap).round() * snap,
            y: (y / snap).round() * snap,
        }
    } else {
        P { x, y }
    }
}

fn bbox_center(geom: &Geometry) -> Option<(f64, f64)> {
    geom.bbox().map(|b| b.center())
}

/// All rings (exterior + interiors) of a polygon geometry as canonical vertex
/// chains without the closing duplicate.
fn polygon_rings(geom: &Geometry, snap: f64) -> Vec<Vec<P>> {
    let ring_pts = |ring: &Ring| -> Vec<P> {
        let mut pts: Vec<P> = Vec::with_capacity(ring.len());
        for c in ring.coords() {
            let p = canonical(c.x, c.y, snap);
            if pts.last().is_none_or(|l| key(*l) != key(p)) {
                pts.push(p);
            }
        }
        while pts.len() >= 2 && key(pts[0]) == key(*pts.last().unwrap()) {
            pts.pop();
        }
        pts
    };
    let mut out = Vec::new();
    match geom {
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            out.push(ring_pts(exterior));
            out.extend(interiors.iter().map(&ring_pts));
        }
        Geometry::MultiPolygon(parts) => {
            for (ext, holes) in parts {
                out.push(ring_pts(ext));
                out.extend(holes.iter().map(&ring_pts));
            }
        }
        _ => {}
    }
    out.retain(|r| r.len() >= 3);
    out
}

fn field_string(fv: &FieldValue) -> String {
    if let Some(s) = fv.as_str() {
        s.to_string()
    } else if let Some(i) = fv.as_i64() {
        i.to_string()
    } else if let Some(f) = fv.as_f64() {
        format!("{f}")
    } else {
        String::new()
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

struct Params {
    include_diagonal: bool,
    snap_tolerance: f64,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let include_diagonal = parse_optional_bool(args, "include_diagonal")?.unwrap_or(true);
    let snap_tolerance = match args.get("snap_tolerance") {
        None | Some(Value::Null) => 0.0,
        Some(Value::Number(num)) => num.as_f64().unwrap_or(0.0),
        Some(Value::String(s)) if s.trim().is_empty() => 0.0,
        Some(Value::String(s)) => s
            .trim()
            .parse::<f64>()
            .map_err(|_| ToolError::Validation("'snap_tolerance' must be a number".into()))?,
        Some(_) => {
            return Err(ToolError::Validation(
                "'snap_tolerance' must be a number".into(),
            ))
        }
    };
    if !(snap_tolerance >= 0.0 && snap_tolerance.is_finite()) {
        return Err(ToolError::Validation(
            "'snap_tolerance' must be a non-negative number".to_string(),
        ));
    }
    Ok(Params {
        include_diagonal,
        snap_tolerance,
    })
}

fn parse_optional_bool(args: &ToolArgs, k: &str) -> Result<Option<bool>, ToolError> {
    match args.get(k) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Bool(b)) => Ok(Some(*b)),
        Some(Value::String(s)) => match s.trim().to_ascii_lowercase().as_str() {
            "" => Ok(None),
            "true" | "1" | "yes" => Ok(Some(true)),
            "false" | "0" | "no" => Ok(Some(false)),
            _ => Err(ToolError::Validation(format!(
                "parameter '{k}' must be a boolean"
            ))),
        },
        Some(_) => Err(ToolError::Validation(format!(
            "parameter '{k}' must be a boolean"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;
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

    /// Build a layer of named square tiles and return its memory:// path.
    fn layer_of(named: &[(&str, Geometry)]) -> String {
        let mut l = Layer::new("grid")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("page_name", FieldType::Text));
        for (name, g) in named {
            l.add_feature(Some(g.clone()), &[("page_name", (*name).into())])
                .unwrap();
        }
        let id = memory_store::put_vector(l);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = CalculateAdjacentFieldsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    /// Read a directional field's value for the tile whose page_name == `page`.
    fn dir(layer: &Layer, page: &str, field: &str) -> String {
        let pi = layer.schema.field_index("page_name").unwrap();
        let di = layer.schema.field_index(field).unwrap();
        for f in layer.iter() {
            if f.attributes[pi].as_str() == Some(page) {
                return f.attributes[di].as_str().unwrap_or("").to_string();
            }
        }
        panic!("page {page} not found");
    }

    /// A 3x3 grid of unit tiles. Names laid out (row 1 = top / north):
    ///   A1 B1 C1
    ///   A2 B2 C2
    ///   A3 B3 C3
    fn grid3x3() -> String {
        let s = 10.0;
        let mut v = Vec::new();
        // (name, col 0..3, row 0..3 top->bottom); y grows north so row0 is top.
        let cells = [
            ("A1", 0, 0),
            ("B1", 1, 0),
            ("C1", 2, 0),
            ("A2", 0, 1),
            ("B2", 1, 1),
            ("C2", 2, 1),
            ("A3", 0, 2),
            ("B3", 1, 2),
            ("C3", 2, 2),
        ];
        // Keep geometries alive: build owned squares first.
        let squares: Vec<(&str, Geometry)> = cells
            .iter()
            .map(|&(n, c, r)| {
                let x0 = c as f64 * s;
                let y0 = (2 - r) as f64 * s; // top row highest y
                (n, square(x0, y0, s))
            })
            .collect();
        for (n, g) in &squares {
            v.push((*n, g.clone()));
        }
        layer_of(&v)
    }

    /// The centre tile B2 has all eight neighbours in the right directions.
    #[test]
    fn center_tile_has_all_eight_neighbors() {
        let input = grid3x3();
        let (out, layer) = run(json!({ "input": input, "page_name_field": "page_name" }));
        assert_eq!(out.outputs["page_count"], json!(9));
        assert_eq!(dir(&layer, "B2", "N"), "B1");
        assert_eq!(dir(&layer, "B2", "S"), "B3");
        assert_eq!(dir(&layer, "B2", "E"), "C2");
        assert_eq!(dir(&layer, "B2", "W"), "A2");
        assert_eq!(dir(&layer, "B2", "NE"), "C1");
        assert_eq!(dir(&layer, "B2", "NW"), "A1");
        assert_eq!(dir(&layer, "B2", "SE"), "C3");
        assert_eq!(dir(&layer, "B2", "SW"), "A3");
    }

    /// A corner tile (A1, top-left) has neighbours only to its E, S, SE and
    /// empty strings elsewhere.
    fn field_val(layer: &Layer, page: &str, field: &str) -> String {
        dir(layer, page, field)
    }
    #[test]
    fn corner_tile_edges_are_empty() {
        let input = grid3x3();
        let (_o, layer) = run(json!({ "input": input, "page_name_field": "page_name" }));
        assert_eq!(field_val(&layer, "A1", "E"), "B1");
        assert_eq!(field_val(&layer, "A1", "S"), "A2");
        assert_eq!(field_val(&layer, "A1", "SE"), "B2");
        assert_eq!(field_val(&layer, "A1", "N"), "");
        assert_eq!(field_val(&layer, "A1", "W"), "");
        assert_eq!(field_val(&layer, "A1", "NW"), "");
    }

    /// include_diagonal=false adds only the four cardinal fields.
    #[test]
    fn diagonal_disabled_drops_corner_fields() {
        let input = grid3x3();
        let (_o, layer) = run(json!({
            "input": input,
            "page_name_field": "page_name",
            "include_diagonal": false
        }));
        assert!(layer.schema.field_index("N").is_some());
        assert!(layer.schema.field_index("E").is_some());
        assert!(layer.schema.field_index("NE").is_none());
        assert!(layer.schema.field_index("SW").is_none());
        assert_eq!(dir(&layer, "B2", "N"), "B1");
    }

    /// Input geometry and existing attributes are preserved on output.
    #[test]
    fn preserves_input_features_and_attributes() {
        let input = grid3x3();
        let (out, layer) = run(json!({ "input": input, "page_name_field": "page_name" }));
        assert_eq!(layer.len(), 9);
        assert_eq!(out.outputs["labels_written"], json!(40)); // 9 tiles, 40 adjacencies
                                                              // page_name column untouched.
        let pi = layer.schema.field_index("page_name").unwrap();
        let mut names: Vec<String> = layer
            .iter()
            .map(|f| f.attributes[pi].as_str().unwrap().to_string())
            .collect();
        names.sort();
        assert_eq!(
            names,
            vec!["A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"]
        );
        // every feature carries geometry.
        assert!(layer.iter().all(|f| f.geometry.is_some()));
    }

    /// A field name that collides with an existing attribute is disambiguated.
    #[test]
    fn colliding_field_name_is_suffixed() {
        let s = 10.0;
        let mut l = Layer::new("grid")
            .with_geom_type(GeometryType::Polygon)
            .with_crs_epsg(3857);
        l.add_field(FieldDef::new("page_name", FieldType::Text));
        l.add_field(FieldDef::new("N", FieldType::Text)); // pre-existing "N"
        l.add_feature(
            Some(square(0.0, 0.0, s)),
            &[("page_name", "A".into()), ("N", "x".into())],
        )
        .unwrap();
        l.add_feature(
            Some(square(0.0, s, s)),
            &[("page_name", "B".into()), ("N", "y".into())],
        )
        .unwrap();
        let id = memory_store::put_vector(l);
        let path = memory_store::make_vector_memory_path(&id);
        let (_o, layer) = run(json!({ "input": path, "page_name_field": "page_name" }));
        assert!(layer.schema.field_index("N_2").is_some());
        // original "N" preserved.
        let pi = layer.schema.field_index("page_name").unwrap();
        let n2 = layer.schema.field_index("N_2").unwrap();
        for f in layer.iter() {
            if f.attributes[pi].as_str() == Some("A") {
                assert_eq!(f.attributes[n2].as_str(), Some("B")); // B is north of A
            }
        }
    }

    #[test]
    fn rejects_missing_input() {
        let args: ToolArgs = serde_json::from_value(json!({})).unwrap();
        assert!(CalculateAdjacentFieldsTool.validate(&args).is_err());
    }

    #[test]
    fn rejects_missing_page_name_field() {
        let args: ToolArgs = serde_json::from_value(json!({ "input": "memory://x" })).unwrap();
        assert!(CalculateAdjacentFieldsTool.validate(&args).is_err());
    }

    #[test]
    fn rejects_bad_parameters() {
        let input = grid3x3();
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input,
            "page_name_field": "page_name",
            "snap_tolerance": "nope"
        }))
        .unwrap();
        assert!(CalculateAdjacentFieldsTool.validate(&args).is_err());
    }
}
