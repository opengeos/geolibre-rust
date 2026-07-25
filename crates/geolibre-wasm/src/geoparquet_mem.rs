//! In-memory GeoParquet reader.
//!
//! `wbvector::geoparquet::read` takes a `Path`, which is unusable in the
//! browser: the wasm32-unknown-unknown build has no filesystem at all. This
//! module reads the same format straight from a byte buffer, so a fetched or
//! drag-and-dropped `.parquet` can be decoded client-side and fed to every
//! other vector entry point (GeoJSON, binary attributes, Arrow IPC).
//!
//! Scope matches the upstream reader: WKB geometry in a binary column named by
//! the `geo` file metadata (default `geometry`), scalar attribute columns, and
//! the GeoParquet 1.1 `bbox` covering column skipped on read.

use std::collections::HashMap;

use parquet::basic::{ConvertedType, Type as PhysicalType};
use parquet::file::metadata::KeyValue;
use parquet::file::reader::{FileReader, SerializedFileReader};
use parquet::record::{Field, Row};
use serde_json::Value;
use wbvector::feature::{FieldDef, FieldType, FieldValue, Layer};
use wbvector::geometry::Geometry;

/// Name of the GeoParquet 1.1 bbox covering column (a `struct<xmin, ymin,
/// xmax, ymax>`), written next to the geometry column and skipped on read.
const BBOX_COL: &str = "bbox";
/// Key under which the upstream writer records exact wbvector field types,
/// which Parquet's own type system cannot represent (Date vs DateTime vs Text).
const WBVECTOR_FIELD_TYPES_KEY: &str = "wbvector_field_types";

/// Read a GeoParquet dataset from bytes into a [`Layer`].
pub fn from_bytes(data: &[u8]) -> Result<Layer, String> {
    let reader = SerializedFileReader::new(bytes::Bytes::from(data.to_vec()))
        .map_err(|e| format!("failed opening parquet bytes: {e}"))?;

    let file_meta = reader.metadata().file_metadata();
    let kv_meta = file_meta.key_value_metadata().map(|v| v.as_slice());
    let geo_meta = parse_geo_metadata(kv_meta)?;
    let declared_types = parse_wbvector_field_types(kv_meta);
    let geom_col = geo_meta
        .primary_column
        .clone()
        .unwrap_or_else(|| "geometry".to_owned());
    let schema_types = infer_types_from_schema(file_meta, &geom_col);

    // Column order comes from the file schema; any column that only shows up in
    // the row data (nested/edge cases) is appended as it is encountered.
    let mut ordered_attr_names: Vec<String> = file_meta
        .schema_descr()
        .columns()
        .iter()
        .filter_map(|c| {
            // The bbox covering is a struct, so its leaves are nested (path
            // ["bbox", "xmin"]); a user attribute named "bbox" is a top-level
            // scalar, so requiring nesting keeps it from being dropped.
            let parts = c.path().parts();
            let root = parts.first().map(String::as_str);
            let is_covering_leaf = root == Some(BBOX_COL) && parts.len() >= 2;
            if root == Some(geom_col.as_str()) || is_covering_leaf {
                None
            } else {
                Some(c.name().to_owned())
            }
        })
        .collect();

    let rows: Vec<Row> = reader
        .get_row_iter(None)
        .map_err(|e| format!("failed to iterate rows: {e}"))?
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| format!("failed reading parquet row: {e}"))?;

    let mut inferred_types: HashMap<String, FieldType> = HashMap::new();
    for name in &ordered_attr_names {
        if let Some(t) = declared_types.get(name).or_else(|| schema_types.get(name)) {
            inferred_types.insert(name.clone(), *t);
        }
    }
    for row in &rows {
        for (name, field) in row.get_column_iter() {
            if name.as_str() == geom_col.as_str() || is_bbox_covering(name, field) {
                continue;
            }
            if !ordered_attr_names.iter().any(|n| n == name) {
                ordered_attr_names.push(name.clone());
            }
            let inferred = declared_types
                .get(name)
                .copied()
                .or_else(|| schema_types.get(name).copied())
                .unwrap_or_else(|| infer_field_type(field));
            let entry = inferred_types.entry(name.clone()).or_insert(inferred);
            *entry = FieldValue::widen_type(*entry, inferred);
        }
    }

    let mut layer = Layer::new("layer");
    if let Some(epsg) = geo_meta.epsg {
        layer = layer.with_crs_epsg(epsg);
    }
    if let Some(wkt) = geo_meta.wkt {
        layer = layer.with_crs_wkt(wkt);
    }
    for name in &ordered_attr_names {
        let ty = inferred_types.get(name).copied().unwrap_or(FieldType::Text);
        layer.add_field(FieldDef::new(name, ty));
    }

    for (idx, row) in rows.into_iter().enumerate() {
        let mut geom = None;
        let mut attrs = vec![FieldValue::Null; layer.schema.len()];

        for (name, field) in row.get_column_iter() {
            if is_bbox_covering(name, field) {
                continue;
            }
            if name.as_str() == geom_col.as_str() {
                geom = geometry_from_field(field)?;
            } else if let Some(i) = layer.schema.field_index(name) {
                let hinted = layer.schema.fields()[i].field_type;
                attrs[i] = field_to_value_with_hint(field, hinted);
            }
        }

        if layer.geom_type.is_none() {
            if let Some(g) = &geom {
                layer.geom_type = Some(g.geom_type());
            }
        }

        layer
            .add_feature(geom, &[])
            .map_err(|e| format!("failed adding feature: {e}"))?;
        if let Some(f) = layer.features.get_mut(idx) {
            f.fid = idx as u64;
            f.attributes = attrs;
        }
    }

    Ok(layer)
}

fn is_bbox_covering(name: &str, field: &Field) -> bool {
    name == BBOX_COL && matches!(field, Field::Group(_))
}

fn geometry_from_field(field: &Field) -> Result<Option<Geometry>, String> {
    match field {
        Field::Null => Ok(None),
        Field::Bytes(bytes) => Geometry::from_wkb(bytes.data())
            .map(Some)
            .map_err(|e| format!("invalid WKB geometry: {e}")),
        other => Err(format!(
            "geometry column must be binary WKB, found {other:?}"
        )),
    }
}

fn infer_field_type(field: &Field) -> FieldType {
    match field {
        Field::Bool(_) => FieldType::Boolean,
        Field::Byte(_) | Field::Short(_) | Field::Int(_) | Field::Long(_) => FieldType::Integer,
        Field::UByte(_) | Field::UShort(_) | Field::UInt(_) | Field::ULong(_) => FieldType::Integer,
        Field::Float(_) | Field::Double(_) => FieldType::Float,
        Field::Bytes(_) => FieldType::Blob,
        _ => FieldType::Text,
    }
}

fn field_to_value(field: &Field) -> FieldValue {
    match field {
        Field::Null => FieldValue::Null,
        Field::Bool(v) => FieldValue::Boolean(*v),
        Field::Byte(v) => FieldValue::Integer(*v as i64),
        Field::Short(v) => FieldValue::Integer(*v as i64),
        Field::Int(v) => FieldValue::Integer(*v as i64),
        Field::Long(v) => FieldValue::Integer(*v),
        Field::UByte(v) => FieldValue::Integer(*v as i64),
        Field::UShort(v) => FieldValue::Integer(*v as i64),
        Field::UInt(v) => FieldValue::Integer(*v as i64),
        Field::ULong(v) => FieldValue::Integer(*v as i64),
        Field::Float(v) => FieldValue::Float(*v as f64),
        Field::Double(v) => FieldValue::Float(*v),
        Field::Str(v) => FieldValue::Text(v.clone()),
        Field::Bytes(v) => FieldValue::Blob(v.data().to_vec()),
        other => FieldValue::Text(format!("{other:?}")),
    }
}

fn field_to_value_with_hint(field: &Field, hint: FieldType) -> FieldValue {
    match (hint, field_to_value(field)) {
        (FieldType::Date, FieldValue::Text(s)) => FieldValue::Date(s),
        (FieldType::DateTime, FieldValue::Text(s)) => FieldValue::DateTime(s),
        (_, v) => v,
    }
}

fn parse_wbvector_field_types(kv: Option<&[KeyValue]>) -> HashMap<String, FieldType> {
    let mut out = HashMap::new();
    let Some(raw) = kv.and_then(|pairs| {
        pairs
            .iter()
            .find(|p| p.key == WBVECTOR_FIELD_TYPES_KEY)
            .and_then(|p| p.value.clone())
    }) else {
        return out;
    };
    // A malformed hint block is not fatal: fall back to schema inference.
    let Ok(Value::Object(obj)) = serde_json::from_str::<Value>(&raw) else {
        return out;
    };
    for (name, value) in obj {
        if let Some(t) = value.as_str().and_then(parse_field_type_name) {
            out.insert(name, t);
        }
    }
    out
}

fn parse_field_type_name(s: &str) -> Option<FieldType> {
    Some(match s {
        "Integer" => FieldType::Integer,
        "Float" => FieldType::Float,
        "Text" => FieldType::Text,
        "Boolean" => FieldType::Boolean,
        "Blob" => FieldType::Blob,
        "Date" => FieldType::Date,
        "DateTime" => FieldType::DateTime,
        "Json" => FieldType::Json,
        _ => return None,
    })
}

fn infer_types_from_schema(
    file_meta: &parquet::file::metadata::FileMetaData,
    geom_col: &str,
) -> HashMap<String, FieldType> {
    let mut out = HashMap::new();
    for col in file_meta.schema_descr().columns() {
        let name = col.name();
        if name == geom_col {
            continue;
        }
        let ty = match col.converted_type() {
            ConvertedType::JSON => FieldType::Json,
            ConvertedType::UTF8 => FieldType::Text,
            ConvertedType::DATE => FieldType::Date,
            _ => match col.physical_type() {
                PhysicalType::BOOLEAN => FieldType::Boolean,
                PhysicalType::INT32 | PhysicalType::INT64 => FieldType::Integer,
                PhysicalType::FLOAT | PhysicalType::DOUBLE => FieldType::Float,
                PhysicalType::BYTE_ARRAY | PhysicalType::FIXED_LEN_BYTE_ARRAY => FieldType::Blob,
                _ => FieldType::Text,
            },
        };
        out.insert(name.to_owned(), ty);
    }
    out
}

#[derive(Debug, Default)]
struct GeoMeta {
    primary_column: Option<String>,
    epsg: Option<u32>,
    wkt: Option<String>,
}

fn parse_geo_metadata(kv: Option<&[KeyValue]>) -> Result<GeoMeta, String> {
    let mut meta = GeoMeta::default();
    let Some(raw) = kv.and_then(|pairs| {
        pairs
            .iter()
            .find(|p| p.key == "geo")
            .and_then(|p| p.value.clone())
    }) else {
        return Ok(meta);
    };

    let v: Value =
        serde_json::from_str(&raw).map_err(|e| format!("invalid 'geo' metadata JSON: {e}"))?;
    meta.primary_column = v
        .get("primary_column")
        .and_then(|x| x.as_str())
        .map(ToOwned::to_owned);

    if let Some(pc) = meta.primary_column.clone() {
        if let Some(col) = v.get("columns").and_then(|c| c.get(&pc)) {
            parse_crs_hint(col, &mut meta);
        }
    }
    Ok(meta)
}

/// Resolve the column's `crs` entry to an EPSG code where possible. Handles the
/// three shapes GeoParquet writers emit: an SRS reference string
/// (`"EPSG:4326"`, `"OGC:CRS84"`), a PROJJSON object with an `id`, and raw WKT.
fn parse_crs_hint(col_meta: &Value, out: &mut GeoMeta) {
    let Some(crs_v) = col_meta.get("crs") else {
        return;
    };

    if let Some(s) = crs_v.as_str() {
        match epsg_from_srs_reference(s) {
            Some(epsg) => out.epsg = Some(epsg),
            None => match epsg_from_wkt_lenient(s) {
                Some(epsg) => out.epsg = Some(epsg),
                None => out.wkt = Some(s.to_owned()),
            },
        }
        return;
    }

    if let Some(obj) = crs_v.as_object() {
        if let Some(wkt) = obj.get("wkt").and_then(|x| x.as_str()) {
            out.wkt = Some(wkt.to_owned());
            out.epsg = epsg_from_wkt_lenient(wkt);
            return;
        }
        if let Some(id) = obj.get("id") {
            let authority = id
                .get("authority")
                .and_then(|a| a.as_str())
                .unwrap_or_default();
            if authority.eq_ignore_ascii_case("EPSG") {
                out.epsg = id.get("code").and_then(|c| c.as_u64()).map(|c| c as u32);
            }
        }
    }
}

/// `"EPSG:4326"` / `"urn:ogc:def:crs:EPSG::4326"` / `"OGC:CRS84"` -> EPSG code.
fn epsg_from_srs_reference(s: &str) -> Option<u32> {
    let t = s.trim();
    if t.eq_ignore_ascii_case("OGC:CRS84")
        || t.eq_ignore_ascii_case("urn:ogc:def:crs:OGC:1.3:CRS84")
    {
        // CRS84 is WGS84 with lon/lat axis order, which is the order every
        // reader here already uses.
        return Some(4326);
    }
    let upper = t.to_ascii_uppercase();
    let code = upper.rsplit_once("EPSG:")?.1.trim_start_matches(':');
    code.parse().ok()
}

/// Pull the EPSG code out of the trailing `AUTHORITY["EPSG","4326"]` (WKT1) or
/// `ID["EPSG",4326]` (WKT2) node without a full WKT parse.
fn epsg_from_wkt_lenient(wkt: &str) -> Option<u32> {
    let upper = wkt.to_ascii_uppercase();
    let start = upper.rfind("AUTHORITY[").or_else(|| upper.rfind("ID["))?;
    let tail = &wkt[start..];
    let inner = tail.split_once('[')?.1;
    let inner = inner.split(']').next()?;
    let mut parts = inner.split(',');
    let authority = parts.next()?.trim().trim_matches('"');
    if !authority.eq_ignore_ascii_case("EPSG") {
        return None;
    }
    parts.next()?.trim().trim_matches('"').parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbvector::feature::FieldDef;
    use wbvector::geometry::GeometryType;

    /// End-to-end check against the format the upstream writer actually emits
    /// (WKB geometry, `geo` metadata, bbox covering column, field-type hints).
    #[test]
    fn reads_geoparquet_written_by_the_upstream_writer() {
        let mut layer = Layer::new("cities")
            .with_geom_type(GeometryType::Point)
            .with_crs_epsg(4326);
        layer.add_field(FieldDef::new("name", FieldType::Text));
        layer.add_field(FieldDef::new("pop", FieldType::Integer));
        layer.add_field(FieldDef::new("area", FieldType::Float));
        layer
            .add_feature(
                Some(Geometry::point(-0.1278, 51.5074)),
                &[
                    ("name", "London".into()),
                    ("pop", 9_000_000i64.into()),
                    ("area", 1572.0.into()),
                ],
            )
            .unwrap();
        layer
            .add_feature(
                Some(Geometry::point(10.7522, 59.9139)),
                &[
                    ("name", "Oslo".into()),
                    ("pop", 700_000i64.into()),
                    ("area", 454.0.into()),
                ],
            )
            .unwrap();

        let path = std::env::temp_dir().join("geolibre_geoparquet_mem_roundtrip.parquet");
        wbvector::geoparquet::write(&layer, &path).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        let _ = std::fs::remove_file(&path);

        let back = from_bytes(&bytes).unwrap();
        assert_eq!(back.len(), 2);
        assert_eq!(back.crs_epsg(), Some(4326));
        assert_eq!(back.geom_type, Some(GeometryType::Point));

        let names: Vec<&str> = back
            .schema
            .fields()
            .iter()
            .map(|f| f.name.as_str())
            .collect();
        assert_eq!(names, vec!["name", "pop", "area"]);
        // The bbox covering column must not leak in as an attribute.
        assert!(!names.contains(&"bbox"));

        let first = &back.features[0];
        assert_eq!(
            first.get(&back.schema, "name").unwrap().as_str(),
            Some("London")
        );
        assert_eq!(
            first.get(&back.schema, "pop").unwrap().as_i64(),
            Some(9_000_000)
        );
        assert_eq!(
            first.get(&back.schema, "area").unwrap().as_f64(),
            Some(1572.0)
        );
        match first.geometry.as_ref().unwrap() {
            Geometry::Point(c) => {
                assert!((c.x - -0.1278).abs() < 1e-9);
                assert!((c.y - 51.5074).abs() < 1e-9);
            }
            other => panic!("expected a point, got {other:?}"),
        }
    }

    #[test]
    fn rejects_bytes_that_are_not_parquet() {
        assert!(from_bytes(b"definitely not parquet").is_err());
    }

    #[test]
    fn srs_reference_forms_resolve_to_epsg() {
        assert_eq!(epsg_from_srs_reference("EPSG:3857"), Some(3857));
        assert_eq!(epsg_from_srs_reference("epsg:4326"), Some(4326));
        assert_eq!(
            epsg_from_srs_reference("urn:ogc:def:crs:EPSG::32617"),
            Some(32617)
        );
        assert_eq!(epsg_from_srs_reference("OGC:CRS84"), Some(4326));
        assert_eq!(epsg_from_srs_reference("not a crs"), None);
    }

    #[test]
    fn wkt_authority_node_resolves_to_epsg() {
        let wkt1 = r#"GEOGCS["WGS 84",DATUM["WGS_1984"],AUTHORITY["EPSG","4326"]]"#;
        assert_eq!(epsg_from_wkt_lenient(wkt1), Some(4326));
        let wkt2 = r#"PROJCRS["WGS 84 / UTM zone 17N",ID["EPSG",32617]]"#;
        assert_eq!(epsg_from_wkt_lenient(wkt2), Some(32617));
        let other = r#"PROJCRS["local",ID["ESRI",102008]]"#;
        assert_eq!(epsg_from_wkt_lenient(other), None);
    }

    #[test]
    fn crs_hint_reads_projjson_id() {
        let col: Value =
            serde_json::from_str(r#"{"crs":{"id":{"authority":"EPSG","code":3857}}}"#).unwrap();
        let mut meta = GeoMeta::default();
        parse_crs_hint(&col, &mut meta);
        assert_eq!(meta.epsg, Some(3857));
    }

    #[test]
    fn field_type_hints_round_trip_through_names() {
        for (name, ty) in [
            ("Integer", FieldType::Integer),
            ("Float", FieldType::Float),
            ("DateTime", FieldType::DateTime),
        ] {
            assert_eq!(parse_field_type_name(name), Some(ty));
        }
        assert_eq!(parse_field_type_name("Nonsense"), None);
    }
}
