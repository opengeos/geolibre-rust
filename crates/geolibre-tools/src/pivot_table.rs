//! GeoLibre tool: reshape a long-form attribute table to wide.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Pivot Table* (Data Management). Rows
//! that repeat per entity — one row per (tract, year), per (zone, land-cover
//! class), per (station, parameter) — are collapsed into a single row per
//! entity with one column per distinct value of the pivot field.
//!
//! The catalog can already *aggregate* long-form data (`summary_statistics`
//! does GROUP BY; the bundled `cross_tabulation` does raster class tables) but
//! it cannot *reshape* it. That is a hard blocker in practice, because
//! `multivariate_clustering`, `dimension_reduction`,
//! `spatially_constrained_multivariate_clustering`, `similarity_search` and
//! `calculate_composite_index` all require one row per feature with one column
//! per variable, while real inputs arrive one row per entity-variable pair.
//!
//! Output is a pure attribute table (no geometry), emitted in first-seen key
//! order so runs are reproducible and diffs stay small.

use std::collections::BTreeMap;
use std::collections::HashMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

/// Hard cap on generated columns; a high-cardinality pivot field is a user
/// error, and silently emitting 50k columns would be worse than failing.
const MAX_PIVOT_COLUMNS: usize = 512;

pub struct PivotTableTool;

impl Tool for PivotTableTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "pivot_table",
            display_name: "Pivot Table",
            summary: "Reshape a long-form attribute table to wide: one output row per unique combination of the identity fields, and one column per distinct value of the pivot field, filled from the value field. Like ArcGIS Pivot Table.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer or attribute table in long form.",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output table path (driver from extension). If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "fields",
                    description: "Comma/semicolon-separated field name(s) identifying an output row; carried through unchanged.",
                    required: true,
                },
                ToolParamSpec {
                    name: "pivot_field",
                    description: "Field whose distinct values become new columns.",
                    required: true,
                },
                ToolParamSpec {
                    name: "value_field",
                    description: "Field supplying the cell values.",
                    required: true,
                },
                ToolParamSpec {
                    name: "aggregate",
                    description: "Applied when one row/column cell has multiple source records: 'first' (default), 'sum', 'mean', 'min', 'max', 'count'.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        let fields = split_list(require_str(args, "fields")?);
        if fields.is_empty() {
            return Err(ToolError::Validation(
                "'fields' must name at least one field".to_string(),
            ));
        }
        require_str(args, "pivot_field")?;
        require_str(args, "value_field")?;
        parse_aggregate(args)?;
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let id_names = split_list(require_str(args, "fields")?);
        let pivot_name = require_str(args, "pivot_field")?;
        let value_name = require_str(args, "value_field")?;
        let agg = parse_aggregate(args)?;

        let layer = load_input_layer(input)?;
        if layer.features.is_empty() {
            return Err(ToolError::Execution("input has no features".to_string()));
        }

        let field_idx = |name: &str| -> Result<usize, ToolError> {
            layer
                .schema
                .field_index(name)
                .ok_or_else(|| ToolError::Validation(format!("field '{name}' not found in input")))
        };
        let id_idx: Vec<usize> = id_names
            .iter()
            .map(|n| field_idx(n))
            .collect::<Result<_, _>>()?;
        let pivot_idx = field_idx(pivot_name)?;
        let value_idx = field_idx(value_name)?;

        // Pass 1: distinct pivot values, in first-seen order.
        let mut pivot_values: Vec<String> = Vec::new();
        let mut pivot_pos: HashMap<String, usize> = HashMap::new();
        for feat in &layer.features {
            let key = value_key(feat.attributes.get(pivot_idx));
            if !pivot_pos.contains_key(&key) {
                if pivot_values.len() >= MAX_PIVOT_COLUMNS {
                    return Err(ToolError::Execution(format!(
                        "pivot_field '{pivot_name}' has more than {MAX_PIVOT_COLUMNS} distinct \
                         values; it is probably not a categorical field"
                    )));
                }
                pivot_pos.insert(key.clone(), pivot_values.len());
                pivot_values.push(key);
            }
        }
        let n_cols = pivot_values.len();
        ctx.progress
            .info(&format!("pivoting into {n_cols} column(s)"));

        // Pass 2: accumulate cells, grouped by the identity-field tuple.
        let mut row_order: Vec<String> = Vec::new();
        let mut row_pos: HashMap<String, usize> = HashMap::new();
        let mut row_ids: Vec<Vec<FieldValue>> = Vec::new();
        let mut cells: Vec<Vec<Accum>> = Vec::new();

        for feat in &layer.features {
            let ids: Vec<FieldValue> = id_idx
                .iter()
                .map(|&i| feat.attributes.get(i).cloned().unwrap_or(FieldValue::Null))
                .collect();
            let rkey = ids
                .iter()
                .map(|v| value_key(Some(v)))
                .collect::<Vec<_>>()
                .join("\u{1f}");
            let r = *row_pos.entry(rkey.clone()).or_insert_with(|| {
                row_order.push(rkey.clone());
                row_ids.push(ids.clone());
                cells.push(vec![Accum::default(); n_cols]);
                row_order.len() - 1
            });
            let c = pivot_pos[&value_key(feat.attributes.get(pivot_idx))];
            cells[r][c].push(feat.attributes.get(value_idx), agg);
        }

        // Column naming: sanitize, then disambiguate collisions.
        let col_names = column_names(&pivot_values);

        // Numeric aggregates always emit Float; 'count' emits Integer; 'first'
        // preserves the source field's type where it is representable.
        let out_type = match agg {
            Aggregate::Count => FieldType::Integer,
            Aggregate::First => layer.schema.fields()[value_idx].field_type,
            _ => FieldType::Float,
        };

        let mut out = Layer::new("pivot");
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        for (name, &idx) in id_names.iter().zip(id_idx.iter()) {
            out.add_field(FieldDef::new(
                name.as_str(),
                layer.schema.fields()[idx].field_type,
            ));
        }
        for name in &col_names {
            out.add_field(FieldDef::new(name.as_str(), out_type));
        }

        let mut filled = 0usize;
        for (r, ids) in row_ids.iter().enumerate() {
            let mut attrs: Vec<(&str, FieldValue)> = Vec::with_capacity(id_names.len() + n_cols);
            for (name, v) in id_names.iter().zip(ids.iter()) {
                attrs.push((name.as_str(), v.clone()));
            }
            for (c, name) in col_names.iter().enumerate() {
                let v = cells[r][c].finish(agg, out_type);
                if !v.is_null() {
                    filled += 1;
                }
                attrs.push((name.as_str(), v));
            }
            out.add_feature(None, &attrs)
                .map_err(|e| ToolError::Execution(format!("failed adding row: {e}")))?;
        }

        let n_rows = row_ids.len();
        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("row_count".to_string(), json!(n_rows));
        outputs.insert("column_count".to_string(), json!(n_cols));
        outputs.insert("filled_cells".to_string(), json!(filled));
        outputs.insert(
            "empty_cells".to_string(),
            json!(n_rows.saturating_mul(n_cols) - filled),
        );
        Ok(ToolRunResult { outputs })
    }
}

// ── Aggregation ─────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Aggregate {
    First,
    Sum,
    Mean,
    Min,
    Max,
    Count,
}

/// Running accumulator for one output cell.
#[derive(Clone, Default)]
struct Accum {
    n: u64,
    sum: f64,
    min: Option<f64>,
    max: Option<f64>,
    first: Option<FieldValue>,
}

impl Accum {
    fn push(&mut self, v: Option<&FieldValue>, agg: Aggregate) {
        let Some(v) = v else { return };
        if v.is_null() {
            return;
        }
        if agg == Aggregate::Count {
            self.n += 1;
            return;
        }
        if self.first.is_none() {
            self.first = Some(v.clone());
        }
        if let Some(x) = v.as_f64() {
            self.n += 1;
            self.sum += x;
            self.min = Some(self.min.map_or(x, |m: f64| m.min(x)));
            self.max = Some(self.max.map_or(x, |m: f64| m.max(x)));
        }
    }

    fn finish(&self, agg: Aggregate, out_type: FieldType) -> FieldValue {
        match agg {
            Aggregate::Count => FieldValue::Integer(self.n as i64),
            Aggregate::First => match &self.first {
                Some(v) => coerce(v, out_type),
                None => FieldValue::Null,
            },
            Aggregate::Sum if self.n > 0 => FieldValue::Float(self.sum),
            Aggregate::Mean if self.n > 0 => FieldValue::Float(self.sum / self.n as f64),
            Aggregate::Min => self.min.map_or(FieldValue::Null, FieldValue::Float),
            Aggregate::Max => self.max.map_or(FieldValue::Null, FieldValue::Float),
            _ => FieldValue::Null,
        }
    }
}

/// Coerces a preserved `first` value into the declared output field type.
fn coerce(v: &FieldValue, ty: FieldType) -> FieldValue {
    match ty {
        FieldType::Integer => v.as_i64().map_or(FieldValue::Null, FieldValue::Integer),
        FieldType::Float => v.as_f64().map_or(FieldValue::Null, FieldValue::Float),
        FieldType::Text => FieldValue::Text(value_key(Some(v))),
        _ => v.clone(),
    }
}

fn parse_aggregate(args: &ToolArgs) -> Result<Aggregate, ToolError> {
    match args
        .get("aggregate")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        None | Some("") | Some("first") => Ok(Aggregate::First),
        Some("sum") => Ok(Aggregate::Sum),
        Some("mean") | Some("avg") => Ok(Aggregate::Mean),
        Some("min") => Ok(Aggregate::Min),
        Some("max") => Ok(Aggregate::Max),
        Some("count") => Ok(Aggregate::Count),
        Some(o) => Err(ToolError::Validation(format!(
            "'aggregate' must be one of first/sum/mean/min/max/count, got '{o}'"
        ))),
    }
}

// ── Naming / keys ───────────────────────────────────────────────────────────

/// Canonical text key for a field value. Floats are formatted without trailing
/// noise so `1.0` and `1` group together, which is what a pivot field wants.
fn value_key(v: Option<&FieldValue>) -> String {
    match v {
        None | Some(FieldValue::Null) => "NULL".to_string(),
        Some(FieldValue::Integer(i)) => i.to_string(),
        Some(FieldValue::Float(f)) => {
            if f.fract() == 0.0 && f.is_finite() && f.abs() < 1e15 {
                format!("{}", *f as i64)
            } else {
                format!("{f}")
            }
        }
        Some(FieldValue::Text(s)) => s.clone(),
        Some(FieldValue::Boolean(b)) => b.to_string(),
        Some(FieldValue::Date(s)) | Some(FieldValue::DateTime(s)) => s.clone(),
        Some(FieldValue::Blob(b)) => format!("blob[{}]", b.len()),
    }
}

/// Turns pivot values into safe, unique column names.
fn column_names(values: &[String]) -> Vec<String> {
    let mut used: HashMap<String, usize> = HashMap::new();
    let mut out = Vec::with_capacity(values.len());
    for v in values {
        let mut base: String = v
            .chars()
            .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
            .collect();
        if base.is_empty() {
            base = "COL".to_string();
        }
        if base.chars().next().is_some_and(|c| c.is_ascii_digit()) {
            base.insert(0, '_');
        }
        let name = match used.get_mut(&base) {
            None => {
                used.insert(base.clone(), 1);
                base
            }
            Some(n) => {
                let candidate = format!("{base}_{n}");
                *n += 1;
                candidate
            }
        };
        out.push(name);
    }
    out
}

fn split_list(s: &str) -> Vec<String> {
    s.split([',', ';'])
        .map(str::trim)
        .filter(|p| !p.is_empty())
        .map(str::to_string)
        .collect()
}

fn require_str<'a>(args: &'a ToolArgs, key: &str) -> Result<&'a str, ToolError> {
    args.get(key)
        .and_then(Value::as_str)
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .ok_or_else(|| ToolError::Validation(format!("missing required string parameter '{key}'")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    /// Long-form table: (site, year, value).
    fn long_table(rows: &[(&str, i64, f64)]) -> String {
        let mut l = Layer::new("long");
        l.add_field(FieldDef::new("site", FieldType::Text));
        l.add_field(FieldDef::new("year", FieldType::Integer));
        l.add_field(FieldDef::new("val", FieldType::Float));
        for (site, year, val) in rows {
            l.add_feature(
                None,
                &[
                    ("site", FieldValue::Text((*site).to_string())),
                    ("year", FieldValue::Integer(*year)),
                    ("val", FieldValue::Float(*val)),
                ],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = PivotTableTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    #[test]
    fn reshapes_long_to_wide() {
        let input = long_table(&[
            ("a", 2020, 1.0),
            ("a", 2021, 2.0),
            ("b", 2020, 3.0),
            ("b", 2021, 4.0),
        ]);
        let (out, layer) = run(json!({
            "input": input, "fields": "site", "pivot_field": "year", "value_field": "val"
        }));
        assert_eq!(out.outputs["row_count"], json!(2));
        assert_eq!(out.outputs["column_count"], json!(2));
        assert_eq!(out.outputs["empty_cells"], json!(0));
        // Columns are named after the pivot values, in first-seen order.
        let c2020 = layer.schema.field_index("_2020").unwrap();
        let c2021 = layer.schema.field_index("_2021").unwrap();
        let site = layer.schema.field_index("site").unwrap();
        assert_eq!(layer.features[0].attributes[site].as_str(), Some("a"));
        assert_eq!(layer.features[0].attributes[c2020].as_f64(), Some(1.0));
        assert_eq!(layer.features[1].attributes[c2021].as_f64(), Some(4.0));
    }

    #[test]
    fn missing_combinations_are_null_not_zero() {
        // 'b' has no 2021 record — that cell must stay null, since 0 would be
        // a fabricated observation.
        let input = long_table(&[("a", 2020, 1.0), ("a", 2021, 2.0), ("b", 2020, 3.0)]);
        let (out, layer) = run(json!({
            "input": input, "fields": "site", "pivot_field": "year", "value_field": "val"
        }));
        assert_eq!(out.outputs["empty_cells"], json!(1));
        let c2021 = layer.schema.field_index("_2021").unwrap();
        assert!(layer.features[1].attributes[c2021].is_null());
    }

    #[test]
    fn aggregate_combines_duplicate_cells() {
        let input = long_table(&[("a", 2020, 1.0), ("a", 2020, 3.0)]);
        let (_o, first) = run(json!({
            "input": input.clone(), "fields": "site", "pivot_field": "year", "value_field": "val"
        }));
        let (_o, summed) = run(json!({
            "input": input.clone(), "fields": "site", "pivot_field": "year",
            "value_field": "val", "aggregate": "sum"
        }));
        let (_o, meaned) = run(json!({
            "input": input, "fields": "site", "pivot_field": "year",
            "value_field": "val", "aggregate": "mean"
        }));
        let c = first.schema.field_index("_2020").unwrap();
        assert_eq!(first.features[0].attributes[c].as_f64(), Some(1.0));
        assert_eq!(summed.features[0].attributes[c].as_f64(), Some(4.0));
        assert_eq!(meaned.features[0].attributes[c].as_f64(), Some(2.0));
    }

    #[test]
    fn multiple_identity_fields_key_the_row() {
        // Keying on (site, year) with year also carried leaves one row each.
        let input = long_table(&[("a", 2020, 1.0), ("a", 2021, 2.0)]);
        let (out, _l) = run(json!({
            "input": input, "fields": "site,year", "pivot_field": "year", "value_field": "val"
        }));
        assert_eq!(out.outputs["row_count"], json!(2));
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            PivotTableTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        assert!(bad(json!({ "input": "t.csv", "pivot_field": "y", "value_field": "v" })).is_err());
        assert!(bad(json!({
            "input": "t.csv", "fields": "s", "pivot_field": "y", "value_field": "v",
            "aggregate": "median"
        }))
        .is_err());
        assert!(bad(json!({
            "input": "t.csv", "fields": "s", "pivot_field": "y", "value_field": "v"
        }))
        .is_ok());
    }

    #[test]
    fn unknown_field_is_rejected_at_run() {
        let input = long_table(&[("a", 2020, 1.0)]);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input, "fields": "nope", "pivot_field": "year", "value_field": "val"
        }))
        .unwrap();
        assert!(PivotTableTool.run(&args, &ctx()).is_err());
    }
}
