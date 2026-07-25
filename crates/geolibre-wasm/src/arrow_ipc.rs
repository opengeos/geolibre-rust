//! Vector data as a GeoArrow-encoded Arrow IPC stream.
//!
//! Where [`crate::binary`] targets deck.gl's own attribute layout, this targets
//! the wider Arrow ecosystem: the bytes returned here go straight into
//! `apache-arrow`'s `tableFromIPC` (zero-copy), into DuckDB-WASM via
//! `insertArrowFromIPCStream`, or into `@geoarrow/deck.gl-layers`. No GeoJSON
//! string, no `JSON.parse`, and one representation shared by the query engine
//! and the renderer.
//!
//! Geometry uses GeoArrow's native (interleaved-coordinate) encoding when every
//! feature in the layer shares one geometry type, which is the case for
//! essentially all real layers, and falls back to `geoarrow.wkb` for mixed
//! layers. Attributes become one Arrow column per field.

use std::sync::Arc;

use arrow_array::{
    Array, ArrayRef, BinaryArray, BooleanArray, FixedSizeListArray, Float64Array, Int64Array,
    ListArray, RecordBatch, StringArray,
};
use arrow_buffer::{NullBuffer, OffsetBuffer, ScalarBuffer};
use arrow_ipc::writer::StreamWriter;
use arrow_schema::{DataType, Field, Schema};
use wasm_bindgen::prelude::*;
use wbvector::feature::{FieldType, FieldValue, Layer};
use wbvector::geometry::{Coord, Geometry, Ring};

const EXT_NAME: &str = "ARROW:extension:name";
const EXT_META: &str = "ARROW:extension:metadata";

/// The GeoArrow encoding chosen for a layer.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Encoding {
    Point,
    LineString,
    Polygon,
    MultiPoint,
    MultiLineString,
    MultiPolygon,
    /// Mixed or collection geometries: serialized as WKB.
    Wkb,
}

impl Encoding {
    fn extension_name(self) -> &'static str {
        match self {
            Self::Point => "geoarrow.point",
            Self::LineString => "geoarrow.linestring",
            Self::Polygon => "geoarrow.polygon",
            Self::MultiPoint => "geoarrow.multipoint",
            Self::MultiLineString => "geoarrow.multilinestring",
            Self::MultiPolygon => "geoarrow.multipolygon",
            Self::Wkb => "geoarrow.wkb",
        }
    }
}

/// Pick the encoding: the shared native type when every geometry agrees, WKB
/// otherwise. A layer with no geometry at all is treated as WKB (all null).
fn choose_encoding(layer: &Layer) -> Encoding {
    let mut chosen: Option<Encoding> = None;
    for geom in layer.features.iter().filter_map(|f| f.geometry.as_ref()) {
        let e = match geom {
            Geometry::Point(_) => Encoding::Point,
            Geometry::LineString(_) => Encoding::LineString,
            Geometry::Polygon { .. } => Encoding::Polygon,
            Geometry::MultiPoint(_) => Encoding::MultiPoint,
            Geometry::MultiLineString(_) => Encoding::MultiLineString,
            Geometry::MultiPolygon(_) => Encoding::MultiPolygon,
            Geometry::GeometryCollection(_) => return Encoding::Wkb,
        };
        match chosen {
            None => chosen = Some(e),
            Some(prev) if prev == e => {}
            Some(_) => return Encoding::Wkb,
        }
    }
    chosen.unwrap_or(Encoding::Wkb)
}

/// True if any coordinate carries Z, in which case every coordinate is written
/// as XYZ (missing Z filled with 0) so the fixed-size list stays uniform.
fn coord_dims(layer: &Layer) -> i32 {
    fn ring_z(r: &Ring) -> bool {
        r.coords().iter().any(|c| c.z.is_some())
    }
    fn geom_z(g: &Geometry) -> bool {
        match g {
            Geometry::Point(c) => c.z.is_some(),
            Geometry::LineString(cs) | Geometry::MultiPoint(cs) => cs.iter().any(|c| c.z.is_some()),
            Geometry::MultiLineString(ls) => ls.iter().flatten().any(|c| c.z.is_some()),
            Geometry::Polygon {
                exterior,
                interiors,
            } => ring_z(exterior) || interiors.iter().any(ring_z),
            Geometry::MultiPolygon(ps) => ps.iter().any(|(e, i)| ring_z(e) || i.iter().any(ring_z)),
            Geometry::GeometryCollection(gs) => gs.iter().any(geom_z),
        }
    }
    let has_z = layer
        .features
        .iter()
        .filter_map(|f| f.geometry.as_ref())
        .any(geom_z);
    if has_z {
        3
    } else {
        2
    }
}

/// Growing interleaved coordinate buffer.
struct Coords {
    values: Vec<f64>,
    dims: i32,
}

impl Coords {
    fn new(dims: i32) -> Self {
        Self {
            values: Vec::new(),
            dims,
        }
    }
    fn len(&self) -> i32 {
        (self.values.len() / self.dims as usize) as i32
    }
    fn push(&mut self, c: &Coord) {
        self.values.push(c.x);
        self.values.push(c.y);
        if self.dims == 3 {
            self.values.push(c.z.unwrap_or(0.0));
        }
    }
    /// Push a ring, closing it when the source stores it open.
    fn push_ring(&mut self, r: &Ring) {
        let coords = r.coords();
        if coords.is_empty() {
            return;
        }
        for c in coords {
            self.push(c);
        }
        let (first, last) = (&coords[0], &coords[coords.len() - 1]);
        if first.x != last.x || first.y != last.y {
            self.push(first);
        }
    }
    /// Finish as a `FixedSizeList<Float64>[dims]` with GeoArrow's field naming.
    fn finish(self) -> FixedSizeListArray {
        self.finish_with_nulls(None)
    }

    fn finish_with_nulls(self, nulls: Option<NullBuffer>) -> FixedSizeListArray {
        let name = if self.dims == 3 { "xyz" } else { "xy" };
        let field = Arc::new(Field::new(name, DataType::Float64, false));
        FixedSizeListArray::new(
            field,
            self.dims,
            Arc::new(Float64Array::from(self.values)),
            nulls,
        )
    }
}

fn list_of(
    name: &str,
    values: ArrayRef,
    offsets: Vec<i32>,
    nulls: Option<NullBuffer>,
) -> ListArray {
    let field = Arc::new(Field::new(name, values.data_type().clone(), false));
    ListArray::new(
        field,
        OffsetBuffer::new(ScalarBuffer::from(offsets)),
        values,
        nulls,
    )
}

/// Validity for the geometry column: null wherever a feature has no geometry.
fn geometry_nulls(layer: &Layer) -> Option<NullBuffer> {
    let has_null = layer.features.iter().any(|f| f.geometry.is_none());
    has_null.then(|| NullBuffer::from_iter(layer.features.iter().map(|f| f.geometry.is_some())))
}

/// Build the geometry column for `encoding`.
///
/// Every native branch walks the features once, appending coordinates and
/// recording the offsets that delimit each nesting level.
fn build_geometry(layer: &Layer, encoding: Encoding, dims: i32) -> Result<ArrayRef, String> {
    let nulls = geometry_nulls(layer);
    let geoms = || layer.features.iter().map(|f| f.geometry.as_ref());

    Ok(match encoding {
        Encoding::Point => {
            let mut coords = Coords::new(dims);
            for g in geoms() {
                match g {
                    Some(Geometry::Point(c)) => coords.push(c),
                    // A null geometry still needs a slot in the fixed-size list.
                    _ => coords.push(&Coord::xy(f64::NAN, f64::NAN)),
                }
            }
            Arc::new(coords.finish_with_nulls(nulls))
        }
        Encoding::LineString | Encoding::MultiPoint => {
            let vertex_name = if encoding == Encoding::MultiPoint {
                "points"
            } else {
                "vertices"
            };
            let mut coords = Coords::new(dims);
            let mut offsets = vec![0i32];
            for g in geoms() {
                match g {
                    Some(Geometry::LineString(cs)) | Some(Geometry::MultiPoint(cs)) => {
                        for c in cs {
                            coords.push(c);
                        }
                    }
                    _ => {}
                }
                offsets.push(coords.len());
            }
            Arc::new(list_of(
                vertex_name,
                Arc::new(coords.finish()),
                offsets,
                nulls,
            ))
        }
        Encoding::Polygon | Encoding::MultiLineString => {
            let (outer_name, inner_name) = if encoding == Encoding::Polygon {
                ("rings", "vertices")
            } else {
                ("linestrings", "vertices")
            };
            let mut coords = Coords::new(dims);
            let mut part_offsets = vec![0i32];
            let mut feat_offsets = vec![0i32];
            for g in geoms() {
                match g {
                    Some(Geometry::Polygon {
                        exterior,
                        interiors,
                    }) => {
                        coords.push_ring(exterior);
                        part_offsets.push(coords.len());
                        for hole in interiors {
                            coords.push_ring(hole);
                            part_offsets.push(coords.len());
                        }
                    }
                    Some(Geometry::MultiLineString(paths)) => {
                        for cs in paths {
                            for c in cs {
                                coords.push(c);
                            }
                            part_offsets.push(coords.len());
                        }
                    }
                    _ => {}
                }
                feat_offsets.push((part_offsets.len() - 1) as i32);
            }
            let inner = list_of(inner_name, Arc::new(coords.finish()), part_offsets, None);
            Arc::new(list_of(outer_name, Arc::new(inner), feat_offsets, nulls))
        }
        Encoding::MultiPolygon => {
            let mut coords = Coords::new(dims);
            let mut ring_offsets = vec![0i32];
            let mut poly_offsets = vec![0i32];
            let mut feat_offsets = vec![0i32];
            for g in geoms() {
                if let Some(Geometry::MultiPolygon(parts)) = g {
                    for (exterior, interiors) in parts {
                        coords.push_ring(exterior);
                        ring_offsets.push(coords.len());
                        for hole in interiors {
                            coords.push_ring(hole);
                            ring_offsets.push(coords.len());
                        }
                        poly_offsets.push((ring_offsets.len() - 1) as i32);
                    }
                }
                feat_offsets.push((poly_offsets.len() - 1) as i32);
            }
            let rings = list_of("vertices", Arc::new(coords.finish()), ring_offsets, None);
            let polys = list_of("rings", Arc::new(rings), poly_offsets, None);
            Arc::new(list_of("polygons", Arc::new(polys), feat_offsets, nulls))
        }
        Encoding::Wkb => {
            let wkb: Vec<Option<Vec<u8>>> = geoms().map(|g| g.map(|g| g.to_wkb())).collect();
            let refs: Vec<Option<&[u8]>> = wkb.iter().map(|b| b.as_deref()).collect();
            Arc::new(BinaryArray::from_opt_vec(refs))
        }
    })
}

/// GeoArrow extension metadata for the geometry field: the layer's CRS, in the
/// two shapes the spec allows us to state confidently.
fn extension_metadata(layer: &Layer) -> String {
    if let Some(epsg) = layer.crs_epsg() {
        return format!("{{\"crs\":\"EPSG:{epsg}\",\"crs_type\":\"authority_code\"}}");
    }
    match layer.crs_wkt() {
        // The WKT flavour is unknown, so pass it through without a crs_type
        // claim rather than mislabelling it.
        Some(wkt) => format!("{{\"crs\":{}}}", json_string(wkt)),
        None => "{}".to_string(),
    }
}

/// Minimal JSON string escaping — enough for CRS WKT, which is ASCII with
/// quotes and the occasional backslash.
fn json_string(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// Build one Arrow column from attribute field `index`.
fn build_attribute(layer: &Layer, index: usize, field_type: FieldType) -> ArrayRef {
    let values = || layer.features.iter().map(move |f| f.get_by_index(index));
    match field_type {
        FieldType::Integer => Arc::new(Int64Array::from_iter(
            values().map(|v| v.and_then(FieldValue::as_i64)),
        )) as ArrayRef,
        FieldType::Float => Arc::new(Float64Array::from_iter(
            values().map(|v| v.and_then(FieldValue::as_f64)),
        )),
        FieldType::Boolean => Arc::new(BooleanArray::from_iter(
            values().map(|v| v.and_then(FieldValue::as_bool)),
        )),
        FieldType::Blob => {
            let owned: Vec<Option<Vec<u8>>> = values()
                .map(|v| v.and_then(FieldValue::as_blob).map(|b| b.to_vec()))
                .collect();
            let refs: Vec<Option<&[u8]>> = owned.iter().map(|b| b.as_deref()).collect();
            Arc::new(BinaryArray::from_opt_vec(refs))
        }
        // Text, Date, DateTime and Json all land in Arrow as UTF-8. Dates keep
        // their ISO-8601 text rather than being reinterpreted as a temporal
        // type, which would need a parse that can fail on real-world data.
        _ => Arc::new(StringArray::from_iter(values().map(|v| match v {
            Some(FieldValue::Text(s))
            | Some(FieldValue::Date(s))
            | Some(FieldValue::DateTime(s)) => Some(s.clone()),
            Some(FieldValue::Null) | None => None,
            Some(other) => Some(crate::binary::field_text(other)),
        }))),
    }
}

/// Serialize a layer as a GeoArrow-encoded Arrow IPC stream.
pub fn to_ipc(layer: &Layer) -> Result<Vec<u8>, String> {
    let encoding = choose_encoding(layer);
    let dims = coord_dims(layer);
    let geometry = build_geometry(layer, encoding, dims)?;

    let mut fields: Vec<Field> = Vec::with_capacity(layer.schema.len() + 1);
    let mut columns: Vec<ArrayRef> = Vec::with_capacity(layer.schema.len() + 1);

    for (i, def) in layer.schema.fields().iter().enumerate() {
        let col = build_attribute(layer, i, def.field_type);
        fields.push(Field::new(def.name.clone(), col.data_type().clone(), true));
        columns.push(col);
    }

    let geom_field = Field::new("geometry", geometry.data_type().clone(), true).with_metadata(
        [
            (EXT_NAME.to_string(), encoding.extension_name().to_string()),
            (EXT_META.to_string(), extension_metadata(layer)),
        ]
        .into_iter()
        .collect(),
    );
    fields.push(geom_field);
    columns.push(geometry);

    let schema = Arc::new(
        Schema::new(fields).with_metadata(
            [("geolibre:layer".to_string(), layer.name.clone())]
                .into_iter()
                .collect(),
        ),
    );
    let batch = RecordBatch::try_new(schema.clone(), columns)
        .map_err(|e| format!("failed building record batch: {e}"))?;

    let mut buf: Vec<u8> = Vec::new();
    {
        let mut writer = StreamWriter::try_new(&mut buf, &schema)
            .map_err(|e| format!("failed opening IPC stream: {e}"))?;
        writer
            .write(&batch)
            .map_err(|e| format!("failed writing record batch: {e}"))?;
        writer
            .finish()
            .map_err(|e| format!("failed finishing IPC stream: {e}"))?;
    }
    Ok(buf)
}

/// Read a vector dataset and return it as a GeoArrow-encoded Arrow IPC stream.
///
/// `format` accepts the same names as [`crate::vector::vector_formats`]. Feed
/// the bytes to `tableFromIPC` (apache-arrow) or
/// `insertArrowFromIPCStream` (DuckDB-WASM).
#[wasm_bindgen]
pub fn vector_to_arrow_ipc(data: &[u8], format: &str) -> Result<Vec<u8>, JsValue> {
    let layer = crate::vector::read_layer(data, format)?;
    to_ipc(&layer).map_err(|e| JsValue::from_str(&e))
}

/// Read a vector dataset, reproject it to `dst_epsg`, and return a GeoArrow
/// Arrow IPC stream. `src_epsg` of `0` uses the layer's own CRS (falling back
/// to EPSG:4326).
#[wasm_bindgen]
pub fn vector_to_arrow_ipc_reproject(
    data: &[u8],
    format: &str,
    dst_epsg: u32,
    src_epsg: u32,
) -> Result<Vec<u8>, JsValue> {
    let layer = crate::vector::read_layer(data, format)?;
    let layer = crate::vector::reproject_layer(layer, dst_epsg, src_epsg)?;
    to_ipc(&layer).map_err(|e| JsValue::from_str(&e))
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow_ipc::reader::StreamReader;
    use wbvector::feature::{FieldDef, FieldType};
    use wbvector::geometry::GeometryType;

    fn read_back(bytes: &[u8]) -> RecordBatch {
        let mut reader = StreamReader::try_new(std::io::Cursor::new(bytes.to_vec()), None).unwrap();
        reader.next().unwrap().unwrap()
    }

    fn ring(coords: &[(f64, f64)]) -> Ring {
        Ring::new(coords.iter().map(|(x, y)| Coord::xy(*x, *y)).collect())
    }

    #[test]
    fn points_round_trip_as_geoarrow_point() {
        let mut layer = Layer::new("cities")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        layer.add_field(FieldDef::new("name", FieldType::Text));
        layer.add_field(FieldDef::new("pop", FieldType::Integer));
        layer
            .add_feature(
                Some(Geometry::point(-0.1, 51.5)),
                &[("name", "London".into()), ("pop", 9_000_000i64.into())],
            )
            .unwrap();
        layer
            .add_feature(
                Some(Geometry::point(10.7, 59.9)),
                &[("name", "Oslo".into()), ("pop", 700_000i64.into())],
            )
            .unwrap();

        let batch = read_back(&to_ipc(&layer).unwrap());
        assert_eq!(batch.num_rows(), 2);
        assert_eq!(batch.num_columns(), 3);

        let geom_field = batch.schema().field(2).clone();
        assert_eq!(geom_field.name(), "geometry");
        assert_eq!(
            geom_field.metadata().get(EXT_NAME).unwrap(),
            "geoarrow.point"
        );
        assert!(geom_field
            .metadata()
            .get(EXT_META)
            .unwrap()
            .contains("EPSG:4326"));

        let coords = batch
            .column(2)
            .as_any()
            .downcast_ref::<FixedSizeListArray>()
            .unwrap();
        assert_eq!(coords.value_length(), 2);
        let flat = coords
            .values()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();
        assert_eq!(flat.values(), &[-0.1, 51.5, 10.7, 59.9]);

        let names = batch
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(names.value(0), "London");
        let pops = batch
            .column(1)
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap();
        assert_eq!(pops.value(1), 700_000);
    }

    #[test]
    fn polygons_nest_rings_and_close_them() {
        let mut layer = Layer::new("polys");
        layer
            .add_feature(
                Some(Geometry::Polygon {
                    // Open ring: the writer must close it.
                    exterior: ring(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]),
                    interiors: vec![ring(&[(0.5, 0.5), (1.0, 0.5), (1.0, 1.0), (0.5, 0.5)])],
                }),
                &[],
            )
            .unwrap();

        let batch = read_back(&to_ipc(&layer).unwrap());
        assert_eq!(
            batch.schema().field(0).metadata().get(EXT_NAME).unwrap(),
            "geoarrow.polygon"
        );
        let polys = batch
            .column(0)
            .as_any()
            .downcast_ref::<ListArray>()
            .unwrap();
        let rings = polys
            .value(0)
            .as_any()
            .downcast_ref::<ListArray>()
            .unwrap()
            .clone();
        assert_eq!(rings.len(), 2); // exterior + one hole
        assert_eq!(rings.value(0).len(), 5); // closed
        assert_eq!(rings.value(1).len(), 4); // already closed
    }

    #[test]
    fn mixed_geometry_types_fall_back_to_wkb() {
        let mut layer = Layer::new("mixed");
        layer
            .add_feature(Some(Geometry::point(0.0, 0.0)), &[])
            .unwrap();
        layer
            .add_feature(
                Some(Geometry::LineString(vec![
                    Coord::xy(0.0, 0.0),
                    Coord::xy(1.0, 1.0),
                ])),
                &[],
            )
            .unwrap();

        assert_eq!(choose_encoding(&layer), Encoding::Wkb);
        let batch = read_back(&to_ipc(&layer).unwrap());
        assert_eq!(
            batch.schema().field(0).metadata().get(EXT_NAME).unwrap(),
            "geoarrow.wkb"
        );
        let wkb = batch
            .column(0)
            .as_any()
            .downcast_ref::<BinaryArray>()
            .unwrap();
        // Both geometries survive the round trip through WKB.
        assert!(matches!(
            Geometry::from_wkb(wkb.value(0)).unwrap(),
            Geometry::Point(_)
        ));
        assert!(matches!(
            Geometry::from_wkb(wkb.value(1)).unwrap(),
            Geometry::LineString(_)
        ));
    }

    #[test]
    fn missing_geometry_becomes_a_null_slot() {
        let mut layer = Layer::new("sparse");
        layer
            .add_feature(Some(Geometry::point(1.0, 1.0)), &[])
            .unwrap();
        layer.add_feature(None, &[]).unwrap();

        let batch = read_back(&to_ipc(&layer).unwrap());
        let geom = batch.column(0);
        assert!(!geom.is_null(0));
        assert!(geom.is_null(1));
    }

    #[test]
    fn z_coordinates_widen_the_fixed_size_list() {
        let mut layer = Layer::new("3d");
        layer
            .add_feature(Some(Geometry::point_z(1.0, 2.0, 3.0)), &[])
            .unwrap();
        layer
            .add_feature(Some(Geometry::point(4.0, 5.0)), &[])
            .unwrap();

        let batch = read_back(&to_ipc(&layer).unwrap());
        let coords = batch
            .column(0)
            .as_any()
            .downcast_ref::<FixedSizeListArray>()
            .unwrap();
        assert_eq!(coords.value_length(), 3);
        let flat = coords
            .values()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();
        assert_eq!(flat.values(), &[1.0, 2.0, 3.0, 4.0, 5.0, 0.0]);
    }

    #[test]
    fn wkt_crs_is_passed_through_without_a_type_claim() {
        let layer = Layer::new("wkt").with_crs_wkt(r#"GEOGCS["WGS 84","quoted"]"#);
        let meta = extension_metadata(&layer);
        assert!(meta.contains("\\\"quoted\\\""), "unescaped WKT: {meta}");
        assert!(!meta.contains("crs_type"));
    }
}
