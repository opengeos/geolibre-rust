//! GeoLibre tool: read an Excel/OpenDocument worksheet into an attribute table.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Excel To Table* (Conversion): import
//! one worksheet from an `.xlsx`/`.xls`/`.xlsb`/`.ods` workbook into a
//! non-spatial attribute table (a `wbvector::Layer` with no geometry) for
//! downstream joins and field transforms.
//!
//! Neither the repo nor the bundled whitebox-wasm suite has ever read a
//! spreadsheet — analysts who receive attribute/lookup tables as Excel had to
//! pre-convert to CSV. The `calamine` crate is pure Rust with no native
//! dependencies (it reads through `std::fs`, so it runs on the WASI `/work`
//! filesystem exactly like `gtfs_to_features`), keeping the tool inside the
//! GDAL/GEOS/PROJ-free stack.
//!
//! Behaviour:
//!   - `sheet` selects a worksheet by name or 0-based index (default: first).
//!   - `cell_range` (A1 notation, e.g. `"B2:E40"`) restricts the imported block;
//!     omitted, the sheet's whole used range is read.
//!   - `field_names_row` is the 1-based row *within the selected range* that
//!     holds the column names (default 1). Rows above it are ignored and rows
//!     below it are data. Set it to `0` to treat every row as data and
//!     auto-generate `Field1`, `Field2`, … names.
//!
//! Per column, the cell types are inspected across the data rows and one
//! `FieldType` is inferred: all-boolean → `Boolean`, all-integer → `Integer`,
//! all-numeric (any float present) → `Float`, otherwise `Text` (dates and mixed
//! columns fall back to their textual form). Empty cells become `Null`.
//!
//! Output is a table layer; with no `output` path it is stored in memory. For a
//! file path, use a container that carries non-spatial rows — GeoParquet
//! (`.parquet`) is the natural table format; GeoJSON writes null-geometry
//! features.

use std::collections::BTreeMap;

use calamine::{open_workbook_auto, Data, Reader};
use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Feature, FieldDef, FieldType, FieldValue, Layer};

use crate::vector_common::{parse_optional_str, write_or_store_layer};

pub struct ExcelToTableTool;

impl Tool for ExcelToTableTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "excel_to_table",
            display_name: "Excel To Table",
            summary: "Read a worksheet from an .xlsx/.xls/.xlsb/.ods workbook into an attribute table, inferring numeric vs text columns.",
            category: ToolCategory::Conversion,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input Excel/OpenDocument workbook path (.xlsx, .xls, .xlsb, or .ods).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output table path (driver from its extension; GeoParquet .parquet recommended). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "sheet",
                    description: "Worksheet to read, by name or 0-based index. Default: the first sheet.",
                    required: false,
                },
                ToolParamSpec {
                    name: "cell_range",
                    description: "Optional A1-notation range to import, e.g. 'B2:E40'. Default: the sheet's whole used range.",
                    required: false,
                },
                ToolParamSpec {
                    name: "field_names_row",
                    description: "1-based row within the selected range holding the column names (default 1). Rows above it are skipped; 0 means no header row (auto-generate Field1, Field2, ...).",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        parse_field_names_row(args)?;
        if let Some(cr) = parse_optional_str(args, "cell_range")? {
            parse_cell_range(cr)?;
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let sheet_arg = parse_optional_str(args, "sheet")?;
        let cell_range = parse_optional_str(args, "cell_range")?
            .map(parse_cell_range)
            .transpose()?;
        let field_names_row = parse_field_names_row(args)?;

        let mut workbook = open_workbook_auto(input)
            .map_err(|e| ToolError::Execution(format!("failed opening workbook '{input}': {e}")))?;

        // Resolve the worksheet: prefer an exact name match, then a 0-based index.
        let names = workbook.sheet_names().to_vec();
        if names.is_empty() {
            return Err(ToolError::Execution(
                "workbook contains no worksheets".to_string(),
            ));
        }
        let sheet_name = match sheet_arg {
            None => names[0].clone(),
            Some(s) => {
                if names.iter().any(|n| n == s) {
                    s.to_string()
                } else if let Ok(idx) = s.parse::<usize>() {
                    names.get(idx).cloned().ok_or_else(|| {
                        ToolError::Validation(format!(
                            "sheet index {idx} out of range (workbook has {} sheet(s))",
                            names.len()
                        ))
                    })?
                } else {
                    return Err(ToolError::Validation(format!(
                        "sheet '{s}' not found; available: {}",
                        names.join(", ")
                    )));
                }
            }
        };

        let full = workbook.worksheet_range(&sheet_name).map_err(|e| {
            ToolError::Execution(format!("failed reading sheet '{sheet_name}': {e}"))
        })?;
        let range = match cell_range {
            Some((start, end)) => full.range(start, end),
            None => full,
        };

        let rows: Vec<&[Data]> = range.rows().collect();
        let ncols = range.width();
        if rows.is_empty() || ncols == 0 {
            return Err(ToolError::Execution(format!(
                "sheet '{sheet_name}' (selected range) is empty"
            )));
        }

        // Split header vs data rows. field_names_row is 1-based within the range;
        // 0 means "no header".
        let (field_names, data_rows): (Vec<String>, &[&[Data]]) = if field_names_row == 0 {
            let names = (1..=ncols).map(|c| format!("Field{c}")).collect();
            (names, &rows[..])
        } else {
            let hidx = field_names_row - 1;
            if hidx >= rows.len() {
                return Err(ToolError::Validation(format!(
                    "field_names_row {field_names_row} is beyond the selected range ({} row(s))",
                    rows.len()
                )));
            }
            let header = header_names(rows[hidx], ncols);
            (header, &rows[hidx + 1..])
        };

        // Infer one FieldType per column across the data rows.
        let field_types: Vec<FieldType> =
            (0..ncols).map(|c| infer_field_type(data_rows, c)).collect();

        let mut layer = Layer::new(sheet_name.clone());
        for (name, ftype) in field_names.iter().zip(&field_types) {
            layer.add_field(FieldDef::new(name.clone(), *ftype));
        }

        let mut row_count = 0usize;
        for row in data_rows {
            // Skip rows that are entirely empty (trailing blank lines Excel keeps).
            if (0..ncols).all(|c| matches!(row.get(c), Some(Data::Empty) | None)) {
                continue;
            }
            let attributes: Vec<FieldValue> = (0..ncols)
                .map(|c| cell_to_value(row.get(c).unwrap_or(&Data::Empty), field_types[c]))
                .collect();
            layer.push(Feature {
                fid: 0,
                geometry: None,
                attributes,
            });
            row_count += 1;
        }

        ctx.progress.info(&format!(
            "imported {row_count} row(s) × {ncols} column(s) from sheet '{sheet_name}'"
        ));

        let out_path = write_or_store_layer(layer, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("sheet".to_string(), json!(sheet_name));
        outputs.insert("row_count".to_string(), json!(row_count));
        outputs.insert("field_count".to_string(), json!(ncols));
        outputs.insert("field_names".to_string(), json!(field_names));
        Ok(ToolRunResult { outputs })
    }
}

/// Builds unique, non-empty field names from a header row, padding to `ncols`.
fn header_names(header: &[Data], ncols: usize) -> Vec<String> {
    let mut names: Vec<String> = Vec::with_capacity(ncols);
    let mut seen: BTreeMap<String, usize> = BTreeMap::new();
    for c in 0..ncols {
        let raw = header.get(c).map(cell_to_string).unwrap_or_default();
        let mut name = raw.trim().to_string();
        if name.is_empty() {
            name = format!("Field{}", c + 1);
        }
        // Disambiguate duplicate column names deterministically.
        let count = seen.entry(name.clone()).or_insert(0);
        if *count > 0 {
            name = format!("{name}_{count}");
        }
        *seen.entry(name.clone()).or_insert(0) += 1;
        names.push(name);
    }
    names
}

/// Infers the storage type for column `c` by scanning every data-row cell.
fn infer_field_type(rows: &[&[Data]], c: usize) -> FieldType {
    let mut any = false;
    let mut all_int = true;
    let mut all_num = true;
    let mut all_bool = true;
    for row in rows {
        match row.get(c) {
            None | Some(Data::Empty) => continue,
            Some(Data::Int(_)) => {
                any = true;
                all_bool = false;
            }
            Some(Data::Float(f)) => {
                any = true;
                all_bool = false;
                // A whole-valued float can still stay integer-typed.
                if f.fract() != 0.0 {
                    all_int = false;
                }
            }
            Some(Data::Bool(_)) => {
                any = true;
                all_int = false;
                all_num = false;
            }
            Some(_) => {
                any = true;
                all_int = false;
                all_num = false;
                all_bool = false;
            }
        }
    }
    if !any {
        FieldType::Text
    } else if all_bool {
        FieldType::Boolean
    } else if all_int {
        FieldType::Integer
    } else if all_num {
        FieldType::Float
    } else {
        FieldType::Text
    }
}

/// Converts one spreadsheet cell to a `FieldValue` in the column's inferred type.
fn cell_to_value(cell: &Data, ftype: FieldType) -> FieldValue {
    if matches!(cell, Data::Empty) {
        return FieldValue::Null;
    }
    match ftype {
        FieldType::Integer => match cell {
            Data::Int(i) => FieldValue::Integer(*i),
            Data::Float(f) => FieldValue::Integer(*f as i64),
            _ => FieldValue::Null,
        },
        FieldType::Float => match cell {
            Data::Int(i) => FieldValue::Float(*i as f64),
            Data::Float(f) => FieldValue::Float(*f),
            _ => FieldValue::Null,
        },
        FieldType::Boolean => match cell {
            Data::Bool(b) => FieldValue::Boolean(*b),
            _ => FieldValue::Null,
        },
        _ => FieldValue::Text(cell_to_string(cell)),
    }
}

/// Renders a cell as text (used for header names and Text columns).
fn cell_to_string(cell: &Data) -> String {
    match cell {
        Data::Empty => String::new(),
        Data::String(s) => s.clone(),
        Data::Bool(b) => b.to_string(),
        Data::Int(i) => i.to_string(),
        Data::Float(f) => f.to_string(),
        other => other.to_string(),
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

/// Parses `field_names_row` (JSON number or numeric string), default 1.
fn parse_field_names_row(args: &ToolArgs) -> Result<usize, ToolError> {
    let n: i64 = match args.get("field_names_row") {
        None | Some(Value::Null) => 1,
        Some(Value::Number(num)) => num.as_i64().ok_or_else(bad_row_num)?,
        Some(Value::String(s)) if s.trim().is_empty() => 1,
        Some(Value::String(s)) => s.trim().parse::<i64>().map_err(|_| bad_row_num())?,
        Some(_) => return Err(bad_row_num()),
    };
    if n < 0 {
        return Err(bad_row_num());
    }
    Ok(n as usize)
}

fn bad_row_num() -> ToolError {
    ToolError::Validation("parameter 'field_names_row' must be a non-negative integer".to_string())
}

/// Inclusive 0-based `((start_row, start_col), (end_row, end_col))` cell window.
type CellRange = ((u32, u32), (u32, u32));

/// Parses an A1-notation range like `"B2:E40"` to inclusive 0-based
/// `(start_row, start_col)` / `(end_row, end_col)` cell coordinates.
fn parse_cell_range(s: &str) -> Result<CellRange, ToolError> {
    let s = s.trim();
    let (a, b) = s
        .split_once(':')
        .ok_or_else(|| ToolError::Validation(format!("cell_range '{s}' must be like 'A1:D20'")))?;
    let start = parse_a1(a.trim())?;
    let end = parse_a1(b.trim())?;
    let lo = (start.0.min(end.0), start.1.min(end.1));
    let hi = (start.0.max(end.0), start.1.max(end.1));
    Ok((lo, hi))
}

/// Parses a single A1 cell reference (e.g. `"B2"`) to 0-based `(row, col)`.
fn parse_a1(s: &str) -> Result<(u32, u32), ToolError> {
    let bad = || ToolError::Validation(format!("invalid A1 cell reference '{s}'"));
    let split = s.find(|c: char| c.is_ascii_digit()).ok_or_else(bad)?;
    let (letters, digits) = s.split_at(split);
    if letters.is_empty() || digits.is_empty() {
        return Err(bad());
    }
    let mut col: u32 = 0;
    for ch in letters.chars() {
        if !ch.is_ascii_alphabetic() {
            return Err(bad());
        }
        col = col
            .checked_mul(26)
            .and_then(|v| v.checked_add((ch.to_ascii_uppercase() as u32) - ('A' as u32) + 1))
            .ok_or_else(bad)?;
    }
    let row: u32 = digits.parse().map_err(|_| bad())?;
    if col == 0 || row == 0 {
        return Err(bad());
    }
    Ok((row - 1, col - 1))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_xlsxwriter::Workbook;
    use std::sync::atomic::{AtomicU64, Ordering};
    use wbcore::{AllowAllCapabilities, ProgressSink};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Cell value for the fixture writer.
    enum V<'a> {
        S(&'a str),
        I(i32),
        F(f64),
        B(bool),
        Blank,
    }

    /// Writes a two-sheet .xlsx fixture to a unique temp path and returns it.
    ///
    /// Sheet "people" carries the primary test table; sheet "extra" exists so
    /// sheet selection (by name and by index) can be exercised.
    fn write_fixture(rows: &[Vec<V>]) -> std::path::PathBuf {
        static SEQ: AtomicU64 = AtomicU64::new(0);
        let n = SEQ.fetch_add(1, Ordering::Relaxed);
        let path =
            std::env::temp_dir().join(format!("excel_to_table_{}_{n}.xlsx", std::process::id()));
        let mut wb = Workbook::new();
        let ws = wb.add_worksheet();
        ws.set_name("people").unwrap();
        for (r, row) in rows.iter().enumerate() {
            for (c, v) in row.iter().enumerate() {
                let (r, c) = (r as u32, c as u16);
                match v {
                    V::S(s) => {
                        ws.write_string(r, c, *s).unwrap();
                    }
                    V::I(i) => {
                        ws.write_number(r, c, *i as f64).unwrap();
                    }
                    V::F(f) => {
                        ws.write_number(r, c, *f).unwrap();
                    }
                    V::B(b) => {
                        ws.write_boolean(r, c, *b).unwrap();
                    }
                    V::Blank => {}
                }
            }
        }
        let ws2 = wb.add_worksheet();
        ws2.set_name("extra").unwrap();
        ws2.write_string(0, 0, "only").unwrap();
        ws2.write_number(1, 0, 42.0).unwrap();
        wb.save(&path).unwrap();
        path
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = ExcelToTableTool.run(&args, &ctx()).unwrap();
        let layer = crate::vector_common::load_input_layer(out.outputs["output"].as_str().unwrap())
            .unwrap();
        (out, layer)
    }

    fn people_fixture() -> std::path::PathBuf {
        write_fixture(&[
            vec![V::S("name"), V::S("pop"), V::S("area"), V::S("coastal")],
            vec![V::S("Alpha"), V::I(100), V::F(12.5), V::B(true)],
            vec![V::S("Beta"), V::I(200), V::F(3.0), V::B(false)],
            vec![V::S("Gamma"), V::I(300), V::F(9.25), V::B(true)],
        ])
    }

    /// Core property: every data row is imported with the header names, and each
    /// column's type is inferred (Text / Integer / Float / Boolean) with values
    /// preserved.
    #[test]
    fn imports_rows_and_infers_column_types() {
        let path = people_fixture();
        let (out, layer) = run(json!({ "input": path.to_str().unwrap() }));
        assert_eq!(out.outputs["row_count"], json!(3));
        assert_eq!(out.outputs["field_count"], json!(4));
        assert_eq!(out.outputs["sheet"], json!("people"));
        assert_eq!(layer.len(), 3);
        assert!(layer.geom_type.is_none(), "table has no geometry");

        let ftype = |name: &str| layer.schema.field(name).unwrap().field_type;
        assert_eq!(ftype("name"), FieldType::Text);
        assert_eq!(ftype("pop"), FieldType::Integer);
        assert_eq!(ftype("area"), FieldType::Float);
        assert_eq!(ftype("coastal"), FieldType::Boolean);

        // Values round-trip.
        let f0 = &layer.features[0];
        assert_eq!(
            f0.get(&layer.schema, "name").unwrap(),
            &FieldValue::Text("Alpha".into())
        );
        assert_eq!(
            f0.get(&layer.schema, "pop").unwrap(),
            &FieldValue::Integer(100)
        );
        assert_eq!(
            f0.get(&layer.schema, "area").unwrap(),
            &FieldValue::Float(12.5)
        );
        assert_eq!(
            f0.get(&layer.schema, "coastal").unwrap(),
            &FieldValue::Boolean(true)
        );
        let _ = std::fs::remove_file(&path);
    }

    /// A column mixing text and numbers falls back to Text; an empty cell in a
    /// numeric column becomes Null (edge cases).
    #[test]
    fn mixed_column_is_text_and_blanks_are_null() {
        let path = write_fixture(&[
            vec![V::S("id"), V::S("val")],
            vec![V::S("a"), V::I(1)],
            vec![V::S("b"), V::Blank],
            vec![V::S("c"), V::S("n/a")], // makes column "val" non-numeric -> Text
        ]);
        let (_, layer) = run(json!({ "input": path.to_str().unwrap() }));
        let vi = layer.schema.field_index("val").unwrap();
        assert_eq!(
            layer.schema.field("val").unwrap().field_type,
            FieldType::Text
        );
        // Row "a": integer rendered as text.
        assert_eq!(
            layer.features[0].attributes[vi],
            FieldValue::Text("1".into())
        );
        // Row "b": blank -> Null.
        assert_eq!(layer.features[1].attributes[vi], FieldValue::Null);
        let _ = std::fs::remove_file(&path);
    }

    /// field_names_row=0 treats the first row as data and auto-names columns.
    #[test]
    fn no_header_autogenerates_field_names() {
        let path = people_fixture();
        let (out, layer) = run(json!({
            "input": path.to_str().unwrap(),
            "field_names_row": 0,
        }));
        // Header row is now data, so 4 rows total.
        assert_eq!(out.outputs["row_count"], json!(4));
        assert!(layer.schema.field_index("Field1").is_some());
        assert!(layer.schema.field_index("name").is_none());
        let _ = std::fs::remove_file(&path);
    }

    /// Sheet can be selected by name or by 0-based index; the second sheet holds
    /// a single value.
    #[test]
    fn selects_sheet_by_name_and_index() {
        let path = people_fixture();
        let by_name = run(json!({ "input": path.to_str().unwrap(), "sheet": "extra" }));
        assert_eq!(by_name.0.outputs["row_count"], json!(1));
        assert_eq!(by_name.0.outputs["sheet"], json!("extra"));

        let by_index = run(json!({ "input": path.to_str().unwrap(), "sheet": "1" }));
        assert_eq!(by_index.0.outputs["sheet"], json!("extra"));
        let _ = std::fs::remove_file(&path);
    }

    /// cell_range restricts the imported block to the requested A1 window.
    #[test]
    fn cell_range_restricts_import() {
        let path = people_fixture();
        // A1:B3 -> header row (name,pop) + 2 data rows.
        let (out, layer) = run(json!({
            "input": path.to_str().unwrap(),
            "cell_range": "A1:B3",
        }));
        assert_eq!(out.outputs["field_count"], json!(2));
        assert_eq!(out.outputs["row_count"], json!(2));
        assert!(layer.schema.field_index("area").is_none());
        assert!(layer.schema.field_index("pop").is_some());
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn parse_a1_maps_letters_to_columns() {
        assert_eq!(parse_a1("A1").unwrap(), (0, 0));
        assert_eq!(parse_a1("B2").unwrap(), (1, 1));
        assert_eq!(parse_a1("AA1").unwrap(), (0, 26));
        assert_eq!(parse_a1("Z10").unwrap(), (9, 25));
        assert!(parse_a1("1A").is_err());
        assert!(parse_a1("A0").is_err());
        assert!(parse_a1("").is_err());
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            ExcelToTableTool.validate(&args)
        };
        // Missing input.
        assert!(bad(json!({})).is_err());
        // Malformed cell_range (no colon).
        assert!(bad(json!({ "input": "a.xlsx", "cell_range": "A1" })).is_err());
        // Invalid A1 reference.
        assert!(bad(json!({ "input": "a.xlsx", "cell_range": "foo:bar" })).is_err());
        // Negative field_names_row.
        assert!(bad(json!({ "input": "a.xlsx", "field_names_row": -1 })).is_err());
        // Valid.
        assert!(
            bad(json!({ "input": "a.xlsx", "cell_range": "B2:E40", "field_names_row": 2 })).is_ok()
        );
    }
}
