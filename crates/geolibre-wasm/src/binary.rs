//! Vector data as flat typed arrays, in deck.gl's binary-attribute layout.
//!
//! [`crate::vector::vector_to_geojson`] hands JavaScript a GeoJSON *string*.
//! The consumer then pays for `JSON.parse` and for re-encoding the parsed
//! objects into the typed arrays a GPU actually wants — twice the work, on a
//! representation that is roughly an order of magnitude larger than the
//! coordinates it carries. For a few hundred thousand features that dominates
//! load time.
//!
//! [`vector_to_binary`] skips both steps: geometry comes back as
//! `Float64Array` positions plus `Uint32Array` index arrays laid out exactly as
//! `@deck.gl/layers`' binary `GeoJsonLayer` expects, and attributes come back
//! as columnar arrays (one per field, one entry per feature) instead of one JS
//! object per feature.
//!
//! Layout, per geometry class (`points` / `lines` / `polygons`):
//! - `positions` — `[x, y, (z), ...]`, `position_size()` values per vertex
//! - `feature_ids` — per vertex, index into that class's feature list
//! - `global_feature_ids` — per vertex, index into the whole layer
//! - `feature_index` — per class feature, its index in the whole layer
//! - lines: `path_indices` — vertex offsets delimiting each path, `n_paths + 1`
//!   entries, ending at the total vertex count
//! - polygons: `primitive_polygon_indices` (ring offsets) and `polygon_indices`
//!   (polygon offsets), same convention
//!
//! Rings are emitted closed (first vertex repeated last) so outlines draw
//! correctly even when the source format stores them open.

use wasm_bindgen::prelude::*;
use wbvector::feature::{FieldType, FieldValue, Layer};
use wbvector::geometry::{Coord, Geometry, Ring};

/// One geometry class's buffers.
#[derive(Default)]
struct ClassBuf {
    positions: Vec<f64>,
    feature_ids: Vec<u32>,
    global_feature_ids: Vec<u32>,
    /// Class-local feature index -> index in the whole layer.
    feature_index: Vec<u32>,
    /// Vertex offsets delimiting paths (lines) or rings (polygons).
    part_indices: Vec<u32>,
    /// Vertex offsets delimiting whole polygons (polygons only).
    polygon_indices: Vec<u32>,
}

impl ClassBuf {
    /// Vertices written so far.
    fn vertex_count(&self, size: usize) -> u32 {
        (self.positions.len() / size) as u32
    }

    /// Claim a class-local id for the layer feature `global`.
    fn claim(&mut self, global: u32) -> u32 {
        self.feature_index.push(global);
        (self.feature_index.len() - 1) as u32
    }

    fn push_coord(&mut self, c: &Coord, size: usize, local: u32, global: u32) {
        self.positions.push(c.x);
        self.positions.push(c.y);
        if size == 3 {
            self.positions.push(c.z.unwrap_or(0.0));
        }
        self.feature_ids.push(local);
        self.global_feature_ids.push(global);
    }

    /// Append a ring, repeating the first vertex when the source leaves it open.
    fn push_ring(&mut self, ring: &Ring, size: usize, local: u32, global: u32) {
        let coords = ring.coords();
        if coords.is_empty() {
            return;
        }
        for c in coords {
            self.push_coord(c, size, local, global);
        }
        let first = &coords[0];
        let last = &coords[coords.len() - 1];
        if first.x != last.x || first.y != last.y {
            self.push_coord(first, size, local, global);
        }
    }
}

/// A layer's geometry and attributes as flat typed arrays.
///
/// Produced by [`vector_to_binary`]. Every accessor copies into a fresh
/// JavaScript typed array, so read each one once and keep the reference.
#[wasm_bindgen]
pub struct VectorBinary {
    feature_count: usize,
    position_size: usize,
    epsg: Option<u32>,
    bounds: Option<[f64; 4]>,
    points: ClassBuf,
    lines: ClassBuf,
    polygons: ClassBuf,
    fields: Vec<(String, FieldType)>,
    /// Per field, per feature. Parallel to `fields`.
    columns: Vec<Vec<FieldValue>>,
}

#[wasm_bindgen]
impl VectorBinary {
    /// Number of features in the source layer.
    #[wasm_bindgen(getter)]
    pub fn feature_count(&self) -> usize {
        self.feature_count
    }

    /// Values per vertex in every `positions` array: 2, or 3 when any geometry
    /// carries a Z coordinate.
    #[wasm_bindgen(getter)]
    pub fn position_size(&self) -> usize {
        self.position_size
    }

    /// EPSG code of the layer's CRS, if it declares one.
    #[wasm_bindgen(getter)]
    pub fn epsg(&self) -> Option<u32> {
        self.epsg
    }

    /// Bounds as `[min_x, min_y, max_x, max_y]`, or an empty array when the
    /// layer has no geometry.
    #[wasm_bindgen(getter)]
    pub fn bbox(&self) -> Vec<f64> {
        self.bounds.map(|b| b.to_vec()).unwrap_or_default()
    }

    /// Attribute schema as JSON: `[{"name":...,"type":...}, ...]`. Field values
    /// themselves come from [`Self::numeric_column`] / [`Self::text_column`].
    #[wasm_bindgen(getter)]
    pub fn schema_json(&self) -> String {
        let fields: Vec<String> = self
            .fields
            .iter()
            .map(|(n, t)| {
                format!(
                    "{{\"name\":\"{}\",\"type\":\"{}\"}}",
                    n.replace('"', "'"),
                    t.as_str()
                )
            })
            .collect();
        format!("[{}]", fields.join(","))
    }

    /// Number of attribute fields.
    #[wasm_bindgen(getter)]
    pub fn field_count(&self) -> usize {
        self.fields.len()
    }

    // ── points ──────────────────────────────────────────────────────────────

    /// Point positions, `position_size` values per vertex.
    pub fn point_positions(&self) -> Vec<f64> {
        self.points.positions.clone()
    }
    /// Per point vertex: index into the point feature list.
    pub fn point_feature_ids(&self) -> Vec<u32> {
        self.points.feature_ids.clone()
    }
    /// Per point vertex: index into the whole layer.
    pub fn point_global_feature_ids(&self) -> Vec<u32> {
        self.points.global_feature_ids.clone()
    }
    /// Per point feature: its index in the whole layer.
    pub fn point_feature_index(&self) -> Vec<u32> {
        self.points.feature_index.clone()
    }

    // ── lines ───────────────────────────────────────────────────────────────

    /// Line positions, `position_size` values per vertex.
    pub fn line_positions(&self) -> Vec<f64> {
        self.lines.positions.clone()
    }
    /// Vertex offsets delimiting each path; ends at the total vertex count.
    pub fn line_path_indices(&self) -> Vec<u32> {
        terminated(&self.lines, self.position_size)
    }
    /// Per line vertex: index into the line feature list.
    pub fn line_feature_ids(&self) -> Vec<u32> {
        self.lines.feature_ids.clone()
    }
    /// Per line vertex: index into the whole layer.
    pub fn line_global_feature_ids(&self) -> Vec<u32> {
        self.lines.global_feature_ids.clone()
    }
    /// Per line feature: its index in the whole layer.
    pub fn line_feature_index(&self) -> Vec<u32> {
        self.lines.feature_index.clone()
    }

    // ── polygons ────────────────────────────────────────────────────────────

    /// Polygon positions, `position_size` values per vertex.
    pub fn polygon_positions(&self) -> Vec<f64> {
        self.polygons.positions.clone()
    }
    /// Vertex offsets delimiting each ring (exterior and holes alike).
    pub fn polygon_primitive_polygon_indices(&self) -> Vec<u32> {
        terminated(&self.polygons, self.position_size)
    }
    /// Vertex offsets delimiting each whole polygon.
    pub fn polygon_indices(&self) -> Vec<u32> {
        let mut v = self.polygons.polygon_indices.clone();
        if !v.is_empty() {
            v.push(self.polygons.vertex_count(self.position_size));
        }
        v
    }
    /// Per polygon vertex: index into the polygon feature list.
    pub fn polygon_feature_ids(&self) -> Vec<u32> {
        self.polygons.feature_ids.clone()
    }
    /// Per polygon vertex: index into the whole layer.
    pub fn polygon_global_feature_ids(&self) -> Vec<u32> {
        self.polygons.global_feature_ids.clone()
    }
    /// Per polygon feature: its index in the whole layer.
    pub fn polygon_feature_index(&self) -> Vec<u32> {
        self.polygons.feature_index.clone()
    }

    // ── attributes ──────────────────────────────────────────────────────────

    /// Field `index` as one `f64` per feature: integers and floats as-is,
    /// booleans as 1/0, everything else (including nulls) as `NaN`.
    pub fn numeric_column(&self, index: usize) -> Result<Vec<f64>, JsValue> {
        self.numeric_column_impl(index)
            .map_err(|e| JsValue::from_str(&e))
    }

    /// Field `index` as UTF-8 bytes for every feature, concatenated. Split it
    /// with [`Self::text_column_offsets`]. Non-text values are stringified;
    /// nulls are empty.
    pub fn text_column(&self, index: usize) -> Result<Vec<u8>, JsValue> {
        self.text_column_impl(index)
            .map_err(|e| JsValue::from_str(&e))
    }

    /// Byte offsets into [`Self::text_column`], `feature_count + 1` entries.
    pub fn text_column_offsets(&self, index: usize) -> Result<Vec<u32>, JsValue> {
        self.text_column_offsets_impl(index)
            .map_err(|e| JsValue::from_str(&e))
    }
}

// Plain-Rust bodies, kept off the `#[wasm_bindgen]` impl block so they can be
// unit-tested natively: constructing a `JsValue` panics off wasm32.
impl VectorBinary {
    fn numeric_column_impl(&self, index: usize) -> Result<Vec<f64>, String> {
        let col = self.column(index)?;
        Ok(col
            .iter()
            .map(|v| match v {
                FieldValue::Integer(i) => *i as f64,
                FieldValue::Float(f) => *f,
                FieldValue::Boolean(b) => {
                    if *b {
                        1.0
                    } else {
                        0.0
                    }
                }
                _ => f64::NAN,
            })
            .collect())
    }

    fn text_column_impl(&self, index: usize) -> Result<Vec<u8>, String> {
        let col = self.column(index)?;
        let mut out = Vec::new();
        for v in col {
            out.extend_from_slice(field_text(v).as_bytes());
        }
        Ok(out)
    }

    fn text_column_offsets_impl(&self, index: usize) -> Result<Vec<u32>, String> {
        let col = self.column(index)?;
        let mut out = Vec::with_capacity(col.len() + 1);
        let mut at = 0u32;
        out.push(at);
        for v in col {
            at += field_text(v).len() as u32;
            out.push(at);
        }
        Ok(out)
    }

    fn column(&self, index: usize) -> Result<&Vec<FieldValue>, String> {
        self.columns.get(index).ok_or_else(|| {
            format!(
                "field index {index} out of range (layer has {} fields)",
                self.fields.len()
            )
        })
    }
}

/// A part-offset array with its terminating total-vertex-count entry. Empty
/// stays empty so consumers can test for "this class has no geometry".
fn terminated(buf: &ClassBuf, size: usize) -> Vec<u32> {
    let mut v = buf.part_indices.clone();
    if !v.is_empty() {
        v.push(buf.vertex_count(size));
    }
    v
}

/// Display form of an attribute value, shared with the Arrow encoder.
pub(crate) fn field_text(v: &FieldValue) -> String {
    match v {
        FieldValue::Null => String::new(),
        FieldValue::Text(s) | FieldValue::Date(s) | FieldValue::DateTime(s) => s.clone(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) => f.to_string(),
        FieldValue::Boolean(b) => b.to_string(),
        FieldValue::Blob(b) => format!("<{} bytes>", b.len()),
    }
}

/// True if any coordinate in the layer carries a Z value.
fn needs_z(layer: &Layer) -> bool {
    fn geom_has_z(g: &Geometry) -> bool {
        match g {
            Geometry::Point(c) => c.z.is_some(),
            Geometry::LineString(cs) | Geometry::MultiPoint(cs) => cs.iter().any(|c| c.z.is_some()),
            Geometry::MultiLineString(ls) => ls.iter().flatten().any(|c| c.z.is_some()),
            Geometry::Polygon {
                exterior,
                interiors,
            } => ring_has_z(exterior) || interiors.iter().any(ring_has_z),
            Geometry::MultiPolygon(ps) => ps
                .iter()
                .any(|(ext, ints)| ring_has_z(ext) || ints.iter().any(ring_has_z)),
            Geometry::GeometryCollection(gs) => gs.iter().any(geom_has_z),
        }
    }
    fn ring_has_z(r: &Ring) -> bool {
        r.coords().iter().any(|c| c.z.is_some())
    }
    layer
        .features
        .iter()
        .filter_map(|f| f.geometry.as_ref())
        .any(geom_has_z)
}

/// Flatten one layer into the binary layout.
pub fn build(layer: &Layer) -> VectorBinary {
    let size = if needs_z(layer) { 3 } else { 2 };
    let mut points = ClassBuf::default();
    let mut lines = ClassBuf::default();
    let mut polygons = ClassBuf::default();
    // Class-local ids are claimed lazily, so a feature that contributes to two
    // classes (a GeometryCollection) gets one id in each.
    let mut bounds: Option<[f64; 4]> = None;

    for (gi, feature) in layer.features.iter().enumerate() {
        let Some(geom) = feature.geometry.as_ref() else {
            continue;
        };
        let gi = gi as u32;
        let mut ids = (None, None, None);
        write_geometry(
            geom,
            size,
            gi,
            &mut ids,
            &mut points,
            &mut lines,
            &mut polygons,
        );
    }

    for chunk in points
        .positions
        .chunks_exact(size)
        .chain(lines.positions.chunks_exact(size))
        .chain(polygons.positions.chunks_exact(size))
    {
        let (x, y) = (chunk[0], chunk[1]);
        bounds = Some(match bounds {
            None => [x, y, x, y],
            Some([a, b, c, d]) => [a.min(x), b.min(y), c.max(x), d.max(y)],
        });
    }

    let fields: Vec<(String, FieldType)> = layer
        .schema
        .fields()
        .iter()
        .map(|f| (f.name.clone(), f.field_type))
        .collect();
    let columns: Vec<Vec<FieldValue>> = (0..fields.len())
        .map(|i| {
            layer
                .features
                .iter()
                .map(|f| f.get_by_index(i).cloned().unwrap_or(FieldValue::Null))
                .collect()
        })
        .collect();

    VectorBinary {
        feature_count: layer.features.len(),
        position_size: size,
        epsg: layer.crs_epsg(),
        bounds,
        points,
        lines,
        polygons,
        fields,
        columns,
    }
}

/// Per-feature class-local ids, claimed on first use: `(point, line, polygon)`.
type LocalIds = (Option<u32>, Option<u32>, Option<u32>);

fn write_geometry(
    geom: &Geometry,
    size: usize,
    global: u32,
    ids: &mut LocalIds,
    points: &mut ClassBuf,
    lines: &mut ClassBuf,
    polygons: &mut ClassBuf,
) {
    match geom {
        Geometry::Point(c) => {
            let local = *ids.0.get_or_insert_with(|| points.claim(global));
            points.push_coord(c, size, local, global);
        }
        Geometry::MultiPoint(cs) => {
            let local = *ids.0.get_or_insert_with(|| points.claim(global));
            for c in cs {
                points.push_coord(c, size, local, global);
            }
        }
        Geometry::LineString(cs) => {
            let local = *ids.1.get_or_insert_with(|| lines.claim(global));
            write_path(cs, size, local, global, lines);
        }
        Geometry::MultiLineString(paths) => {
            let local = *ids.1.get_or_insert_with(|| lines.claim(global));
            for cs in paths {
                write_path(cs, size, local, global, lines);
            }
        }
        Geometry::Polygon {
            exterior,
            interiors,
        } => {
            let local = *ids.2.get_or_insert_with(|| polygons.claim(global));
            write_polygon(exterior, interiors, size, local, global, polygons);
        }
        Geometry::MultiPolygon(parts) => {
            let local = *ids.2.get_or_insert_with(|| polygons.claim(global));
            for (exterior, interiors) in parts {
                write_polygon(exterior, interiors, size, local, global, polygons);
            }
        }
        Geometry::GeometryCollection(gs) => {
            for g in gs {
                write_geometry(g, size, global, ids, points, lines, polygons);
            }
        }
    }
}

fn write_path(coords: &[Coord], size: usize, local: u32, global: u32, buf: &mut ClassBuf) {
    if coords.is_empty() {
        return;
    }
    buf.part_indices.push(buf.vertex_count(size));
    for c in coords {
        buf.push_coord(c, size, local, global);
    }
}

fn write_polygon(
    exterior: &Ring,
    interiors: &[Ring],
    size: usize,
    local: u32,
    global: u32,
    buf: &mut ClassBuf,
) {
    if exterior.coords().is_empty() {
        return;
    }
    buf.polygon_indices.push(buf.vertex_count(size));
    buf.part_indices.push(buf.vertex_count(size));
    buf.push_ring(exterior, size, local, global);
    for hole in interiors {
        if hole.coords().is_empty() {
            continue;
        }
        buf.part_indices.push(buf.vertex_count(size));
        buf.push_ring(hole, size, local, global);
    }
}

/// Read a vector dataset and return its geometry and attributes as flat typed
/// arrays — no GeoJSON string, no `JSON.parse`.
///
/// `format` accepts the same names as [`crate::vector::vector_formats`].
#[wasm_bindgen]
pub fn vector_to_binary(data: &[u8], format: &str) -> Result<VectorBinary, JsValue> {
    Ok(build(&crate::vector::read_layer(data, format)?))
}

/// Read a vector dataset, reproject it to `dst_epsg`, and return flat typed
/// arrays. `src_epsg` of `0` uses the layer's own CRS (falling back to
/// EPSG:4326), matching [`crate::vector::vector_to_geojson_reproject`].
#[wasm_bindgen]
pub fn vector_to_binary_reproject(
    data: &[u8],
    format: &str,
    dst_epsg: u32,
    src_epsg: u32,
) -> Result<VectorBinary, JsValue> {
    let layer = crate::vector::read_layer(data, format)?;
    Ok(build(&crate::vector::reproject_layer(
        layer, dst_epsg, src_epsg,
    )?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbvector::feature::{FieldDef, FieldType};
    use wbvector::geometry::GeometryType;

    fn ring(coords: &[(f64, f64)]) -> Ring {
        Ring::new(coords.iter().map(|(x, y)| Coord::xy(*x, *y)).collect())
    }

    #[test]
    fn points_flatten_to_positions_and_ids() {
        let mut layer = Layer::new("pts").with_geom_type(GeometryType::Point);
        layer.add_field(FieldDef::new("pop", FieldType::Integer));
        layer
            .add_feature(Some(Geometry::point(1.0, 2.0)), &[("pop", 10i64.into())])
            .unwrap();
        layer
            .add_feature(Some(Geometry::point(3.0, 4.0)), &[("pop", 20i64.into())])
            .unwrap();

        let b = build(&layer);
        assert_eq!(b.position_size(), 2);
        assert_eq!(b.point_positions(), vec![1.0, 2.0, 3.0, 4.0]);
        assert_eq!(b.point_feature_ids(), vec![0, 1]);
        assert_eq!(b.point_global_feature_ids(), vec![0, 1]);
        assert_eq!(b.point_feature_index(), vec![0, 1]);
        assert_eq!(b.numeric_column(0).unwrap(), vec![10.0, 20.0]);
        assert_eq!(b.bbox(), vec![1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn multi_line_paths_are_delimited_and_terminated() {
        let mut layer = Layer::new("lines");
        let geom = Geometry::MultiLineString(vec![
            vec![Coord::xy(0.0, 0.0), Coord::xy(1.0, 1.0)],
            vec![
                Coord::xy(2.0, 2.0),
                Coord::xy(3.0, 3.0),
                Coord::xy(4.0, 4.0),
            ],
        ]);
        layer.add_feature(Some(geom), &[]).unwrap();

        let b = build(&layer);
        // Two paths of 2 and 3 vertices: offsets 0, 2, terminated at 5.
        assert_eq!(b.line_path_indices(), vec![0, 2, 5]);
        assert_eq!(b.line_positions().len(), 10);
        assert_eq!(b.line_feature_ids(), vec![0; 5]);
    }

    #[test]
    fn polygon_rings_close_and_holes_get_their_own_offsets() {
        let mut layer = Layer::new("polys");
        // Deliberately open ring: the builder must close it.
        let exterior = ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]);
        let hole = ring(&[(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0), (1.0, 1.0)]);
        layer
            .add_feature(
                Some(Geometry::Polygon {
                    exterior,
                    interiors: vec![hole],
                }),
                &[],
            )
            .unwrap();

        let b = build(&layer);
        // Exterior closes to 5 vertices, hole is already closed at 5.
        assert_eq!(b.polygon_primitive_polygon_indices(), vec![0, 5, 10]);
        assert_eq!(b.polygon_indices(), vec![0, 10]);
        let pos = b.polygon_positions();
        assert_eq!(pos.len(), 20);
        assert_eq!(&pos[0..2], &[0.0, 0.0]);
        assert_eq!(&pos[8..10], &[0.0, 0.0]); // ring closed back to the start
    }

    #[test]
    fn geometry_collection_feeds_every_class_with_its_own_local_id() {
        let mut layer = Layer::new("mixed");
        layer
            .add_feature(Some(Geometry::point(9.0, 9.0)), &[])
            .unwrap();
        layer
            .add_feature(
                Some(Geometry::GeometryCollection(vec![
                    Geometry::point(0.0, 0.0),
                    Geometry::LineString(vec![Coord::xy(0.0, 0.0), Coord::xy(1.0, 0.0)]),
                ])),
                &[],
            )
            .unwrap();

        let b = build(&layer);
        assert_eq!(b.point_feature_index(), vec![0, 1]);
        assert_eq!(b.point_feature_ids(), vec![0, 1]);
        assert_eq!(b.line_feature_index(), vec![1]);
        assert_eq!(b.line_feature_ids(), vec![0, 0]);
        assert_eq!(b.line_path_indices(), vec![0, 2]);
    }

    #[test]
    fn z_coordinates_widen_every_class_to_size_three() {
        let mut layer = Layer::new("3d");
        layer
            .add_feature(Some(Geometry::point_z(1.0, 2.0, 3.0)), &[])
            .unwrap();
        layer
            .add_feature(Some(Geometry::point(4.0, 5.0)), &[])
            .unwrap();

        let b = build(&layer);
        assert_eq!(b.position_size(), 3);
        assert_eq!(b.point_positions(), vec![1.0, 2.0, 3.0, 4.0, 5.0, 0.0]);
    }

    #[test]
    fn text_columns_split_on_byte_offsets() {
        let mut layer = Layer::new("named");
        layer.add_field(FieldDef::new("name", FieldType::Text));
        layer
            .add_feature(
                Some(Geometry::point(0.0, 0.0)),
                &[("name", "Zürich".into())],
            )
            .unwrap();
        layer
            .add_feature(Some(Geometry::point(1.0, 1.0)), &[("name", "Oslo".into())])
            .unwrap();

        let b = build(&layer);
        let bytes = b.text_column(0).unwrap();
        let offsets = b.text_column_offsets(0).unwrap();
        assert_eq!(offsets.len(), 3);
        let first = std::str::from_utf8(&bytes[offsets[0] as usize..offsets[1] as usize]).unwrap();
        let second = std::str::from_utf8(&bytes[offsets[1] as usize..offsets[2] as usize]).unwrap();
        assert_eq!(first, "Zürich"); // multi-byte char: offsets are bytes, not chars
        assert_eq!(second, "Oslo");
        assert!(b.numeric_column_impl(1).is_err());
    }

    #[test]
    fn empty_classes_stay_empty() {
        let mut layer = Layer::new("pts");
        layer
            .add_feature(Some(Geometry::point(0.0, 0.0)), &[])
            .unwrap();
        let b = build(&layer);
        assert!(b.line_path_indices().is_empty());
        assert!(b.polygon_indices().is_empty());
        assert!(b.polygon_primitive_polygon_indices().is_empty());
    }
}
