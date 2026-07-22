//! GeoLibre tool: SQL-style `GROUP BY` attribute aggregation.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Summary Statistics* (Analysis)
//! tool. It groups the input attribute table by one or more *case fields* and,
//! for each unique case-field tuple, computes a configurable list of
//! `(field, statistic)` aggregates plus a group record `count`.
//!
//! Supported statistics: `count` (non-null values), `sum`, `mean`, `min`,
//! `max`, `std` (sample standard deviation), `first`, and `last`. `mean` and
//! `std` use Welford's online algorithm so they stay numerically stable over a
//! single streaming pass — no second pass, no catastrophic cancellation.
//!
//! The output is a **pure attribute table** (no geometry): one row per unique
//! case-field combination, ordered deterministically by the case-field key.
//! Passing no case field collapses to a single summary row over the whole
//! table (which is how you get a one-line grand total). The degenerate
//! count-only configuration reproduces ArcGIS **Frequency**.
//!
//! The bundled whitebox `vector_summary_statistics` summarizes a *single*
//! field over the *whole* layer with no grouping; `cross_tabulation` only does
//! a two-field contingency count. This tool is the daily-driver
//! `GROUP BY(fields) -> statistics` workflow neither covers.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{Feature, FieldDef, FieldType, FieldValue, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct SummaryStatisticsTool;

impl Tool for SummaryStatisticsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "summary_statistics",
            display_name: "Summary Statistics",
            summary: "SQL-style GROUP BY over one or more case fields: for each unique case-field combination, compute a list of (field, statistic) aggregates (count/sum/mean/min/max/std/first/last) plus a group record count. Output is a pure attribute table. No case field yields a single grand-total row (like ArcGIS's Summary Statistics; count-only reproduces Frequency).",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer or attribute table to aggregate.",
                    required: true,
                },
                ToolParamSpec {
                    name: "statistics",
                    description: "Semicolon/comma-separated list of 'field statistic' pairs (also accepts 'field:statistic'), e.g. 'population sum; population mean; elev max'. Statistics: count, sum, mean, min, max, std, first, last.",
                    required: true,
                },
                ToolParamSpec {
                    name: "case_fields",
                    description: "Optional comma/semicolon-separated case (GROUP BY) field name(s). If omitted, a single summary row is produced over the whole table.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Optional output table path (driver from its extension; GeoParquet or CSV). If omitted, the table is stored in memory. Columns: the case fields, a record 'count', then one '<field>_<statistic>' column per requested pair.",
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
        let schema = &layer.schema;

        // Resolve requested field names to positional indices up front.
        for field in prm
            .case_fields
            .iter()
            .chain(prm.stats.iter().map(|s| &s.field))
        {
            if schema.field_index(field).is_none() {
                return Err(ToolError::Execution(format!(
                    "field '{field}' not found in the input table"
                )));
            }
        }

        // Accumulate one group per unique case-field tuple.
        //
        // `groups` maps the joined string key -> its accumulator; `order` keeps
        // the insertion order so the final rows can be sorted deterministically
        // by the case-field key.
        let mut groups: BTreeMap<String, Group> = BTreeMap::new();

        for feature in &layer.features {
            let key = group_key(feature, schema, &prm.case_fields);
            let group = groups.entry(key).or_insert_with(|| Group {
                case_values: prm
                    .case_fields
                    .iter()
                    .map(|f| feature.get(schema, f).cloned().unwrap_or(FieldValue::Null))
                    .collect(),
                count: 0,
                accs: vec![Acc::default(); prm.stats.len()],
            });
            group.count += 1;
            for (acc, stat) in group.accs.iter_mut().zip(prm.stats.iter()) {
                let value = feature.get(schema, &stat.field).ok();
                acc.push(value);
            }
        }

        // Build the output table schema: case fields, record count, then one
        // column per (field, statistic) pair.
        let mut table = Layer::new("summary_statistics");
        for field in &prm.case_fields {
            let ty = schema
                .field(field)
                .map(|d| d.field_type)
                .unwrap_or(FieldType::Text);
            table.add_field(FieldDef::new(field.clone(), ty));
        }
        table.add_field(FieldDef::new("count", FieldType::Integer));
        for stat in &prm.stats {
            table.add_field(FieldDef::new(stat.column_name(), stat.kind.output_type()));
        }

        // Emit one row per group. `BTreeMap` already iterates in sorted key
        // order, giving deterministic output.
        let mut fid = 0u64;
        for group in groups.values() {
            let mut attributes: Vec<FieldValue> = group.case_values.clone();
            attributes.push(FieldValue::Integer(group.count as i64));
            for (acc, stat) in group.accs.iter().zip(prm.stats.iter()) {
                attributes.push(acc.finish(stat.kind));
            }
            table.push(Feature {
                fid,
                geometry: None,
                attributes,
            });
            fid += 1;
        }

        let group_count = table.features.len();
        ctx.progress.info(&format!(
            "aggregated {} record(s) into {} group(s)",
            layer.features.len(),
            group_count
        ));

        let out_path = write_or_store_layer(table, output)?;

        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("group_count".to_string(), json!(group_count));
        outputs.insert("input_records".to_string(), json!(layer.features.len()));
        outputs.insert("case_fields".to_string(), json!(prm.case_fields));
        outputs.insert(
            "statistics".to_string(),
            json!(prm.stats.iter().map(Stat::column_name).collect::<Vec<_>>()),
        );
        Ok(ToolRunResult { outputs })
    }
}

// ── Grouping ──────────────────────────────────────────────────────────────────

/// Field-value separator for the composite group key. Unit Separator (0x1F) is
/// a control char that will not appear in real attribute text.
const KEY_SEP: char = '\u{1f}';

/// Builds the composite group key for a feature from its case-field values.
/// With no case fields this is the empty string, yielding one grand-total group.
fn group_key(feature: &Feature, schema: &wbvector::Schema, case_fields: &[String]) -> String {
    let mut key = String::new();
    for (i, field) in case_fields.iter().enumerate() {
        if i > 0 {
            key.push(KEY_SEP);
        }
        let val = feature.get(schema, field).ok();
        key.push_str(&value_key(val));
    }
    key
}

/// Canonical string form of a value for keying. `Null` gets a distinct sentinel
/// so nulls group together and never collide with the text "null".
fn value_key(val: Option<&FieldValue>) -> String {
    match val {
        None | Some(FieldValue::Null) => "\u{0}<null>".to_string(),
        Some(FieldValue::Integer(v)) => format!("i{v}"),
        Some(FieldValue::Float(v)) => format!("f{v}"),
        Some(FieldValue::Boolean(b)) => format!("b{b}"),
        Some(FieldValue::Text(s)) | Some(FieldValue::Date(s)) | Some(FieldValue::DateTime(s)) => {
            format!("s{s}")
        }
        Some(FieldValue::Blob(b)) => format!("x{}", b.len()),
    }
}

struct Group {
    /// The case-field values (first feature in the group), used for output.
    case_values: Vec<FieldValue>,
    /// Number of records in the group (the `count` column / ArcGIS FREQUENCY).
    count: usize,
    /// One accumulator per requested `(field, statistic)` pair.
    accs: Vec<Acc>,
}

// ── Statistic accumulators ────────────────────────────────────────────────────

/// A single-pass accumulator carrying everything any statistic needs. Numeric
/// mean/variance use Welford's algorithm (stable one-pass mean and M2).
#[derive(Clone, Default)]
struct Acc {
    /// Count of non-null values seen.
    n_nonnull: u64,
    /// Count of finite numeric values (denominator for mean/std).
    n_num: u64,
    /// Running mean (Welford).
    mean: f64,
    /// Running sum of squared deviations from the mean (Welford M2).
    m2: f64,
    /// Plain running sum.
    sum: f64,
    /// Min / max of numeric values.
    min: f64,
    max: f64,
    /// First / last non-null value seen (original type preserved).
    first: Option<FieldValue>,
    last: Option<FieldValue>,
}

impl Acc {
    fn push(&mut self, value: Option<&FieldValue>) {
        let Some(value) = value else { return };
        if value.is_null() {
            return;
        }
        self.n_nonnull += 1;
        if self.first.is_none() {
            self.first = Some(value.clone());
        }
        self.last = Some(value.clone());

        if let Some(x) = value.as_f64() {
            if x.is_finite() {
                self.n_num += 1;
                // Welford update.
                let delta = x - self.mean;
                self.mean += delta / self.n_num as f64;
                self.m2 += delta * (x - self.mean);
                self.sum += x;
                if self.n_num == 1 {
                    self.min = x;
                    self.max = x;
                } else {
                    self.min = self.min.min(x);
                    self.max = self.max.max(x);
                }
            }
        }
    }

    /// Materializes the requested statistic as a `FieldValue`. Statistics with
    /// no data (e.g. `mean` over zero numeric values, `std` over fewer than two)
    /// yield `Null`.
    fn finish(&self, kind: StatKind) -> FieldValue {
        match kind {
            StatKind::Count => FieldValue::Integer(self.n_nonnull as i64),
            StatKind::Sum => {
                if self.n_num == 0 {
                    FieldValue::Null
                } else {
                    FieldValue::Float(self.sum)
                }
            }
            StatKind::Mean => {
                if self.n_num == 0 {
                    FieldValue::Null
                } else {
                    FieldValue::Float(self.mean)
                }
            }
            StatKind::Min => {
                if self.n_num == 0 {
                    FieldValue::Null
                } else {
                    FieldValue::Float(self.min)
                }
            }
            StatKind::Max => {
                if self.n_num == 0 {
                    FieldValue::Null
                } else {
                    FieldValue::Float(self.max)
                }
            }
            StatKind::Std => {
                if self.n_num < 2 {
                    FieldValue::Null
                } else {
                    // Sample standard deviation (divide by n-1).
                    FieldValue::Float((self.m2 / (self.n_num - 1) as f64).sqrt())
                }
            }
            StatKind::First => self.first.clone().unwrap_or(FieldValue::Null),
            StatKind::Last => self.last.clone().unwrap_or(FieldValue::Null),
        }
    }
}

// ── Parameters ────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum StatKind {
    Count,
    Sum,
    Mean,
    Min,
    Max,
    Std,
    First,
    Last,
}

impl StatKind {
    fn parse(s: &str) -> Option<StatKind> {
        match s.trim().to_ascii_lowercase().as_str() {
            "count" => Some(StatKind::Count),
            "sum" => Some(StatKind::Sum),
            "mean" | "avg" | "average" => Some(StatKind::Mean),
            "min" | "minimum" => Some(StatKind::Min),
            "max" | "maximum" => Some(StatKind::Max),
            "std" | "stddev" | "std_dev" => Some(StatKind::Std),
            "first" => Some(StatKind::First),
            "last" => Some(StatKind::Last),
            _ => None,
        }
    }

    fn suffix(self) -> &'static str {
        match self {
            StatKind::Count => "count",
            StatKind::Sum => "sum",
            StatKind::Mean => "mean",
            StatKind::Min => "min",
            StatKind::Max => "max",
            StatKind::Std => "std",
            StatKind::First => "first",
            StatKind::Last => "last",
        }
    }

    /// The output column type. Numeric summaries are Float; count is Integer;
    /// first/last echo the source and default to Text.
    fn output_type(self) -> FieldType {
        match self {
            StatKind::Count => FieldType::Integer,
            StatKind::First | StatKind::Last => FieldType::Text,
            _ => FieldType::Float,
        }
    }
}

struct Stat {
    field: String,
    kind: StatKind,
}

impl Stat {
    fn column_name(&self) -> String {
        format!("{}_{}", self.field, self.kind.suffix())
    }
}

struct Params {
    case_fields: Vec<String>,
    stats: Vec<Stat>,
}

fn parse_params(args: &ToolArgs) -> Result<Params, ToolError> {
    let case_fields = parse_optional_str(args, "case_fields")?
        .map(|s| {
            s.split([',', ';'])
                .map(str::trim)
                .filter(|t| !t.is_empty())
                .map(str::to_string)
                .collect::<Vec<_>>()
        })
        .unwrap_or_default();

    let raw = parse_optional_str(args, "statistics")?.ok_or_else(|| {
        ToolError::Validation("missing required parameter 'statistics'".to_string())
    })?;

    let mut stats = Vec::new();
    for pair in raw.split([';', ',', '\n']) {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }
        // Within a pair, "field:stat" takes precedence, else split on the last
        // run of whitespace so "field name stat" keeps a spaced field name.
        let (field, stat) = if let Some((f, s)) = pair.split_once(':') {
            (f.trim(), s.trim())
        } else {
            pair.rsplit_once(char::is_whitespace)
                .map(|(f, s)| (f.trim(), s.trim()))
                .ok_or_else(|| {
                    ToolError::Validation(format!(
                        "statistics entry '{pair}' must be 'field statistic' or 'field:statistic'"
                    ))
                })?
        };
        if field.is_empty() {
            return Err(ToolError::Validation(format!(
                "statistics entry '{pair}' is missing a field name"
            )));
        }
        let kind = StatKind::parse(stat).ok_or_else(|| {
            ToolError::Validation(format!(
                "unknown statistic '{stat}' (expected count/sum/mean/min/max/std/first/last)"
            ))
        })?;
        stats.push(Stat {
            field: field.to_string(),
            kind,
        });
    }

    if stats.is_empty() {
        return Err(ToolError::Validation(
            "parameter 'statistics' listed no valid (field, statistic) pairs".to_string(),
        ));
    }

    Ok(Params { case_fields, stats })
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

    /// Builds a geometry-less table with a Text `grp` field and a Float `val`.
    fn table(rows: &[(&str, Option<f64>)]) -> String {
        let mut layer = Layer::new("t");
        layer.add_field(FieldDef::new("grp", FieldType::Text));
        layer.add_field(FieldDef::new("val", FieldType::Float));
        for &(g, v) in rows {
            let val = v.map(FieldValue::Float).unwrap_or(FieldValue::Null);
            layer.push(Feature {
                fid: 0,
                geometry: None,
                attributes: vec![FieldValue::Text(g.to_string()), val],
            });
        }
        let id = memory_store::put_vector(layer);
        memory_store::make_vector_memory_path(&id)
    }

    fn run(args: serde_json::Value) -> Layer {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = SummaryStatisticsTool.run(&args, &ctx()).unwrap();
        load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap()
    }

    fn cell(layer: &Layer, row: usize, name: &str) -> FieldValue {
        layer.features[row]
            .get(&layer.schema, name)
            .unwrap()
            .clone()
    }

    /// Locates the output row whose `grp` equals `g`.
    fn row_of(layer: &Layer, g: &str) -> usize {
        (0..layer.features.len())
            .find(|&i| cell(layer, i, "grp").as_str() == Some(g))
            .unwrap()
    }

    #[test]
    fn group_by_sum_mean_count() {
        let input = table(&[
            ("a", Some(1.0)),
            ("a", Some(3.0)),
            ("b", Some(10.0)),
            ("b", Some(20.0)),
            ("b", Some(30.0)),
        ]);
        let layer = run(json!({
            "input": input,
            "case_fields": "grp",
            "statistics": "val sum; val mean; val count",
        }));
        assert_eq!(layer.features.len(), 2);
        let a = row_of(&layer, "a");
        let b = row_of(&layer, "b");
        assert_eq!(cell(&layer, a, "count"), FieldValue::Integer(2));
        assert_eq!(cell(&layer, a, "val_sum").as_f64().unwrap(), 4.0);
        assert_eq!(cell(&layer, a, "val_mean").as_f64().unwrap(), 2.0);
        assert_eq!(cell(&layer, b, "val_sum").as_f64().unwrap(), 60.0);
        assert_eq!(cell(&layer, b, "val_mean").as_f64().unwrap(), 20.0);
        assert_eq!(cell(&layer, b, "val_count"), FieldValue::Integer(3));
    }

    #[test]
    fn std_is_sample_stddev() {
        // Values 2,4,4,4,5,5,7,9 -> sample std = 2.13809...
        let rows: Vec<(&str, Option<f64>)> = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]
            .iter()
            .map(|&v| ("g", Some(v)))
            .collect();
        let input = table(&rows);
        let layer = run(json!({
            "input": input, "case_fields": "grp", "statistics": "val std; val min; val max",
        }));
        let std = cell(&layer, 0, "val_std").as_f64().unwrap();
        assert!((std - 2.138089935).abs() < 1e-6, "got {std}");
        assert_eq!(cell(&layer, 0, "val_min").as_f64().unwrap(), 2.0);
        assert_eq!(cell(&layer, 0, "val_max").as_f64().unwrap(), 9.0);
    }

    #[test]
    fn no_case_field_is_single_summary_row() {
        let input = table(&[("a", Some(1.0)), ("b", Some(2.0)), ("c", Some(3.0))]);
        let layer = run(json!({ "input": input, "statistics": "val sum; val count" }));
        assert_eq!(layer.features.len(), 1);
        assert_eq!(cell(&layer, 0, "count"), FieldValue::Integer(3));
        assert_eq!(cell(&layer, 0, "val_sum").as_f64().unwrap(), 6.0);
        // No case column exists in the degenerate output.
        assert!(layer.schema.field_index("grp").is_none());
    }

    #[test]
    fn count_ignores_nulls_but_record_count_does_not() {
        let input = table(&[("a", Some(1.0)), ("a", None), ("a", Some(3.0))]);
        let layer = run(json!({
            "input": input, "case_fields": "grp", "statistics": "val count; val mean",
        }));
        // Record count = 3, but non-null value count = 2, mean = 2.0.
        assert_eq!(cell(&layer, 0, "count"), FieldValue::Integer(3));
        assert_eq!(cell(&layer, 0, "val_count"), FieldValue::Integer(2));
        assert_eq!(cell(&layer, 0, "val_mean").as_f64().unwrap(), 2.0);
    }

    #[test]
    fn first_and_last_preserve_order() {
        let input = table(&[("a", Some(10.0)), ("a", Some(20.0)), ("a", Some(30.0))]);
        let layer = run(json!({
            "input": input, "case_fields": "grp", "statistics": "val first; val last",
        }));
        assert_eq!(cell(&layer, 0, "val_first").as_f64().unwrap(), 10.0);
        assert_eq!(cell(&layer, 0, "val_last").as_f64().unwrap(), 30.0);
    }

    #[test]
    fn empty_group_std_and_mean_are_null() {
        // Group "a" has only nulls -> mean/std null; count column = 0.
        let input = table(&[("a", None), ("a", None)]);
        let layer = run(json!({
            "input": input, "case_fields": "grp", "statistics": "val mean; val std; val count",
        }));
        assert!(cell(&layer, 0, "val_mean").is_null());
        assert!(cell(&layer, 0, "val_std").is_null());
        assert_eq!(cell(&layer, 0, "val_count"), FieldValue::Integer(0));
    }

    #[test]
    fn colon_and_whitespace_pair_syntax_both_work() {
        let input = table(&[("a", Some(1.0)), ("a", Some(3.0))]);
        let layer = run(json!({
            "input": input, "case_fields": "grp", "statistics": "val:sum, val:mean",
        }));
        assert_eq!(cell(&layer, 0, "val_sum").as_f64().unwrap(), 4.0);
        assert_eq!(cell(&layer, 0, "val_mean").as_f64().unwrap(), 2.0);
    }

    #[test]
    fn rejects_bad_parameters() {
        let tool = SummaryStatisticsTool;
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            tool.validate(&args)
        };
        assert!(bad(json!({})).is_err(), "missing input");
        assert!(
            bad(json!({ "input": "x.geojson" })).is_err(),
            "missing statistics"
        );
        assert!(
            bad(json!({ "input": "x.geojson", "statistics": "val bogus" })).is_err(),
            "unknown statistic"
        );
        assert!(
            bad(json!({ "input": "x.geojson", "statistics": "val" })).is_err(),
            "pair missing statistic"
        );
        assert!(
            bad(json!({ "input": "x.geojson", "statistics": "val sum" })).is_ok(),
            "valid single pair"
        );
    }

    #[test]
    fn unknown_field_errors_at_run() {
        let input = table(&[("a", Some(1.0))]);
        let args: ToolArgs = serde_json::from_value(json!({
            "input": input, "case_fields": "grp", "statistics": "does_not_exist sum",
        }))
        .unwrap();
        assert!(SummaryStatisticsTool.run(&args, &ctx()).is_err());
    }
}
