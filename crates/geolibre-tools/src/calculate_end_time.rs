//! GeoLibre tool: fill an end-time field from the next record's start time.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Calculate End Time* (Data Management).
//!
//! ## Small, and genuinely missing
//!
//! Several of this crate's temporal tools reason about an *interval* — how long
//! a state held, how long an entity stayed — but observational data arrives as
//! instants: one row per reading, with a timestamp and no duration. Deriving
//! each interval's end from the next observation's start, within an entity
//! grouping, is the standard bridge and nothing in either registry did it.
//!
//! Downstream consumers: `find_dwell_locations`, `space_time_kernel_density`,
//! `emerging_hot_spot_analysis`, `estimate_time_to_event`, `reconstruct_tracks`,
//! `trace_proximity_events`, `find_space_time_matches`.
//!
//! `transform_fields` does per-value numeric transforms (zscore, log, binning),
//! not row-ordered temporal fill; `summary_statistics` aggregates groups without
//! writing back per row; `calculate_adjacent_fields` works on polygon adjacency,
//! not record order.
//!
//! ## Ordering, and why output order is not sort order
//!
//! Rows are ordered by `(group, start)` to compute the fill, but **indices** are
//! sorted rather than the features themselves, so the output preserves input row
//! order and geometry is never shuffled. A tool that silently reordered a layer
//! would break every join the caller had already made against it.
//!
//! ## Unparseable starts
//!
//! A row whose start time is null or unparseable is excluded from the ordering
//! entirely and gets a null end, and the count is reported. Sorting such rows to
//! the front (which is what a `unwrap_or(0.0)` would do) would corrupt the first
//! real row's interval by making it start at the epoch.

use std::collections::BTreeMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Layer};

use crate::args_common::{choice_or, opt_f64};
use crate::find_dwell_locations::parse_time_value;
use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

const LAST_RECORD: [&str; 3] = ["null", "same_as_start", "duration"];

/// How the end value is written, derived from the target column's declared
/// type so the written values always agree with the schema.
#[derive(Clone, Copy, PartialEq, Eq)]
enum EndKind {
    Text,
    Integer,
    Float,
}

pub struct CalculateEndTimeTool;

impl Tool for CalculateEndTimeTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "calculate_end_time",
            display_name: "Calculate End Time",
            summary: "Fills an end-time field from the start time of the next record within each entity group, turning instant-in-time observations into time-extent records (ArcGIS Calculate End Time). transform_fields does per-value transforms and summary_statistics aggregates groups, but nothing performed a row-ordered temporal fill.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input vector layer or table.",
                    required: true,
                },
                ToolParamSpec {
                    name: "start_field",
                    description: "Timestamp field: a numeric epoch, or an ISO-8601 date/date-time string.",
                    required: true,
                },
                ToolParamSpec {
                    name: "end_field",
                    description: "Field receiving the end time (default 'END_TIME'). Created if absent, overwritten if present.",
                    required: false,
                },
                ToolParamSpec {
                    name: "id_fields",
                    description: "Comma- or semicolon-separated grouping fields. End times are derived within a group only, never across group boundaries.",
                    required: false,
                },
                ToolParamSpec {
                    name: "last_record",
                    description: "What to write for each group's final row: 'null' (default), 'same_as_start', or 'duration' (start + default_duration).",
                    required: false,
                },
                ToolParamSpec {
                    name: "default_duration",
                    description: "Seconds added to the start for last_record 'duration'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output layer. If omitted, stored in memory.",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        if parse_optional_str(args, "input")?.is_none() {
            return Err(ToolError::Validation(
                "missing required string parameter 'input'".to_string(),
            ));
        }
        if parse_optional_str(args, "start_field")?.is_none() {
            return Err(ToolError::Validation(
                "missing required string parameter 'start_field'".to_string(),
            ));
        }
        parse_list(args, "id_fields")?;
        let last = choice_or(args, "last_record", &LAST_RECORD, "null")?;
        let dur = opt_f64(args, "default_duration")?;
        if last == "duration" {
            match dur {
                None => {
                    return Err(ToolError::Validation(
                        "last_record 'duration' requires 'default_duration' in seconds".to_string(),
                    ))
                }
                Some(d) if !d.is_finite() || d < 0.0 => {
                    return Err(ToolError::Validation(
                        "'default_duration' must be a non-negative number of seconds".to_string(),
                    ))
                }
                _ => {}
            }
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = parse_optional_str(args, "input")?
            .ok_or_else(|| ToolError::Validation("missing required parameter 'input'".into()))?;
        let start_field = parse_optional_str(args, "start_field")?.ok_or_else(|| {
            ToolError::Validation("missing required parameter 'start_field'".into())
        })?;
        let end_name = parse_optional_str(args, "end_field")?
            .unwrap_or("END_TIME")
            .to_string();
        let id_fields = parse_list(args, "id_fields")?;
        let last_record = choice_or(args, "last_record", &LAST_RECORD, "null")?;
        let default_duration = opt_f64(args, "default_duration")?.unwrap_or(0.0);
        let output = parse_optional_str(args, "output")?;

        let layer = load_input_layer(input)?;
        let start_idx = layer.schema.field_index(start_field).ok_or_else(|| {
            ToolError::Validation(format!(
                "start_field '{start_field}' not found in the input layer"
            ))
        })?;
        let mut id_idx = Vec::with_capacity(id_fields.len());
        for f in &id_fields {
            id_idx.push(layer.schema.field_index(f).ok_or_else(|| {
                ToolError::Validation(format!("id_field '{f}' not found in the input layer"))
            })?);
        }

        // The start field's declared type decides how the end is written back:
        // a text timestamp column gets an ISO-8601 string, a numeric one gets a
        // number. Writing an epoch integer into a column of ISO strings would
        // make the layer self-inconsistent.
        // When the end column already exists, ITS declared type decides how the
        // value is written; only a freshly created column follows the start
        // field. Writing a Float into a Text column (or the reverse) produces a
        // layer whose values disagree with its own schema, which the shapefile
        // and GeoParquet writers both expect to be consistent.
        let type_source = layer
            .schema
            .field_index(&end_name)
            .map(|i| layer.schema.fields()[i].field_type)
            .unwrap_or_else(|| layer.schema.fields()[start_idx].field_type);
        // Three cases, not two: collapsing to text-vs-not writes a Float into
        // an existing Integer column holding whole epoch seconds, which is the
        // same schema/value mismatch the Text case fixed.
        let out_kind = match type_source {
            FieldType::Text | FieldType::Date | FieldType::DateTime => EndKind::Text,
            FieldType::Integer => EndKind::Integer,
            _ => EndKind::Float,
        };

        let n = layer.features.len();
        let mut starts: Vec<Option<f64>> = Vec::with_capacity(n);
        for f in layer.iter() {
            starts.push(
                f.attributes
                    .get(start_idx)
                    .filter(|v| !matches!(v, FieldValue::Null))
                    .and_then(parse_time_value),
            );
        }

        // Group key per row, materialized once. Building it inside the
        // comparator would allocate a String on every comparison, so sorting
        // alone would cost O(n log n) allocations.
        let keys: Vec<String> = (0..n)
            .map(|i| {
                id_idx
                    .iter()
                    .map(|&k| {
                        layer.features[i]
                            .attributes
                            .get(k)
                            .map(key_string)
                            .unwrap_or_default()
                    })
                    .collect::<Vec<_>>()
                    .join("\u{1}")
            })
            .collect();

        // Sort INDICES, not features: the output must keep input row order.
        let mut order: Vec<usize> = (0..n).filter(|&i| starts[i].is_some()).collect();
        order.sort_by(|&a, &b| {
            keys[a]
                .cmp(&keys[b])
                .then_with(|| starts[a].unwrap().total_cmp(&starts[b].unwrap()))
                .then_with(|| a.cmp(&b))
        });

        let mut ends: Vec<Option<f64>> = vec![None; n];
        let mut groups = 0_u64;
        let mut i = 0;
        while i < order.len() {
            let key = &keys[order[i]];
            let mut j = i;
            while j < order.len() && &keys[order[j]] == key {
                j += 1;
            }
            groups += 1;
            // Within the group, each row's end is the next row's start; the
            // final row follows `last_record`.
            for k in i..j - 1 {
                ends[order[k]] = starts[order[k + 1]];
            }
            let last = order[j - 1];
            ends[last] = match last_record {
                "same_as_start" => starts[last],
                "duration" => starts[last].map(|s| s + default_duration),
                _ => None,
            };
            i = j;
        }

        let skipped = n - order.len();
        let filled = ends.iter().filter(|e| e.is_some()).count();
        ctx.progress.info(&format!(
            "{n} row(s), {groups} group(s), {filled} end time(s) filled, {skipped} skipped"
        ));

        // Build the output, appending the end field only when it is new.
        let mut out = Layer::new(&layer.name);
        if let Some(gt) = layer.geom_type {
            out = out.with_geom_type(gt);
        }
        if let Some(e) = layer.crs_epsg() {
            out = out.with_crs_epsg(e);
        }
        for f in layer.schema.fields() {
            out.add_field(f.clone());
        }
        let end_idx = match layer.schema.field_index(&end_name) {
            Some(i) => i,
            None => {
                out.add_field(FieldDef::new(
                    end_name.as_str(),
                    match out_kind {
                        EndKind::Text => FieldType::Text,
                        EndKind::Integer => FieldType::Integer,
                        EndKind::Float => FieldType::Float,
                    },
                ));
                out.schema.fields().len() - 1
            }
        };

        for (i, f) in layer.iter().enumerate() {
            let mut attrs = f.attributes.clone();
            // A pre-existing end field is at end_idx; a fresh one needs the
            // vector grown to match the widened schema.
            while attrs.len() <= end_idx {
                attrs.push(FieldValue::Null);
            }
            attrs[end_idx] = match (ends[i], out_kind) {
                (None, _) => FieldValue::Null,
                (Some(v), EndKind::Text) => FieldValue::Text(format_iso8601(v)),
                (Some(v), EndKind::Integer) => {
                    // Rounding would move the endpoint: a next start of 10.4
                    // stored as 10 leaves an interval that no longer reaches
                    // the following record. Refuse rather than silently
                    // shortening it.
                    if !v.is_finite() || v.fract() != 0.0 || v.abs() > i64::MAX as f64 {
                        return Err(ToolError::Validation(format!(
                            "end time {v} is not a whole number of seconds, but '{end_name}' is \
                             an Integer column; storing it would move the endpoint. Use a Float \
                             end_field, or supply whole-second starts and default_duration"
                        )));
                    }
                    FieldValue::Integer(v as i64)
                }
                (Some(v), EndKind::Float) => FieldValue::Float(v),
            };
            out.push(wbvector::Feature {
                fid: f.fid,
                geometry: f.geometry.clone(),
                attributes: attrs,
            });
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("row_count".to_string(), json!(n));
        outputs.insert("group_count".to_string(), json!(groups));
        outputs.insert("filled_count".to_string(), json!(filled));
        outputs.insert("skipped_count".to_string(), json!(skipped));
        outputs.insert("end_field".to_string(), json!(end_name));
        Ok(ToolRunResult { outputs })
    }
}

/// Formats epoch seconds back as `YYYY-MM-DDTHH:MM:SS`, the inverse of the
/// parser in `find_dwell_locations`.
fn format_iso8601(secs: f64) -> String {
    let total = secs.floor() as i64;
    let days = total.div_euclid(86_400);
    let rem = total.rem_euclid(86_400);
    let (y, m, d) = civil_from_days(days);
    format!(
        "{y:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}",
        rem / 3600,
        (rem % 3600) / 60,
        rem % 60
    )
}

/// Inverse of `find_dwell_locations::days_from_civil` (Howard Hinnant's
/// civil-calendar algorithm).
fn civil_from_days(z: i64) -> (i64, i64, i64) {
    let z = z + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    (if m <= 2 { y + 1 } else { y }, m, d)
}

fn key_string(v: &FieldValue) -> String {
    match v {
        FieldValue::Null => "\u{0}NULL".to_string(),
        FieldValue::Integer(i) => i.to_string(),
        FieldValue::Float(f) => format!("{f:.12}"),
        FieldValue::Text(s) | FieldValue::Date(s) | FieldValue::DateTime(s) => s.clone(),
        FieldValue::Boolean(b) => b.to_string(),
        FieldValue::Blob(b) => format!("blob[{}]", b.len()),
    }
}

/// Accepts a delimited string or a JSON array, matching `transpose_fields`.
///
/// Reading only strings meant `["id"]` yielded an empty list, every row landed
/// in one group, and end times were filled straight across entity boundaries —
/// the exact failure the module doc warns about, with no error to the caller.
fn parse_list(args: &ToolArgs, key: &str) -> Result<Vec<String>, ToolError> {
    match args.get(key) {
        None | Some(Value::Null) => Ok(Vec::new()),
        Some(Value::String(s)) => Ok(s
            .split([',', ';'])
            .map(str::trim)
            .filter(|p| !p.is_empty())
            .map(str::to_string)
            .collect()),
        Some(Value::Array(items)) => items
            .iter()
            .map(|v| {
                v.as_str().map(str::to_string).ok_or_else(|| {
                    ToolError::Validation(format!("every entry of '{key}' must be a string"))
                })
            })
            .collect(),
        Some(other) => Err(ToolError::Validation(format!(
            "'{key}' must be a delimited string or an array of strings; got {other}"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use wbcore::{AllowAllCapabilities, ProgressSink};
    use wbvector::{Feature, Geometry, GeometryType};

    struct NullProgress;
    impl ProgressSink for NullProgress {}

    fn ctx() -> ToolContext<'static> {
        ToolContext {
            progress: &NullProgress,
            capabilities: &AllowAllCapabilities,
        }
    }

    fn store(l: Layer) -> String {
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    /// Numeric-epoch readings for two entities, deliberately interleaved and
    /// out of time order.
    fn readings() -> String {
        let mut l = Layer::new("obs").with_geom_type(GeometryType::Point);
        l.add_field(FieldDef::new("id", FieldType::Text));
        l.add_field(FieldDef::new("t", FieldType::Float));
        for (id, t) in [
            ("a", 30.0),
            ("b", 100.0),
            ("a", 10.0),
            ("b", 300.0),
            ("a", 20.0),
        ] {
            l.add_feature(
                Some(Geometry::point(0.0, 0.0)),
                &[("id", id.into()), ("t", t.into())],
            )
            .unwrap();
        }
        store(l)
    }

    fn run(args: Value) -> (Layer, ToolRunResult) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let res = CalculateEndTimeTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(res.outputs["output"].as_str().unwrap()).unwrap();
        (layer, res)
    }

    fn col(l: &Layer, name: &str) -> Vec<FieldValue> {
        let i = l.schema.field_index(name).unwrap();
        l.iter().map(|f| f.attributes[i].clone()).collect()
    }

    #[test]
    fn each_row_gets_the_next_start_within_its_group() {
        let (layer, res) = run(json!({
            "input": readings(), "start_field": "t", "id_fields": "id",
        }));
        // Input order preserved: a30, b100, a10, b300, a20.
        assert_eq!(
            col(&layer, "END_TIME"),
            vec![
                FieldValue::Null,         // a30 is a's last
                FieldValue::Float(300.0), // b100 -> b300
                FieldValue::Float(20.0),  // a10 -> a20
                FieldValue::Null,         // b300 is b's last
                FieldValue::Float(30.0),  // a20 -> a30
            ]
        );
        assert_eq!(res.outputs["group_count"], json!(2));
        assert_eq!(res.outputs["filled_count"], json!(3));
    }

    #[test]
    fn output_row_order_matches_input_row_order() {
        // Sorting the features rather than the indices would shuffle the layer
        // and break every join already made against it.
        let (layer, _) = run(json!({
            "input": readings(), "start_field": "t", "id_fields": "id",
        }));
        assert_eq!(
            col(&layer, "t"),
            vec![
                FieldValue::Float(30.0),
                FieldValue::Float(100.0),
                FieldValue::Float(10.0),
                FieldValue::Float(300.0),
                FieldValue::Float(20.0),
            ]
        );
    }

    #[test]
    fn end_times_never_cross_a_group_boundary() {
        // Without grouping, a10 would chain to b100. With it, a's chain stops.
        let (grouped, _) = run(json!({
            "input": readings(), "start_field": "t", "id_fields": "id",
        }));
        let (ungrouped, res) = run(json!({"input": readings(), "start_field": "t"}));
        assert_eq!(res.outputs["group_count"], json!(1));
        // a30 (row 0) is the global maximum only when ungrouped... it is not:
        // ungrouped the order is 10,20,30,100,300, so a30 -> 100.
        assert_eq!(col(&ungrouped, "END_TIME")[0], FieldValue::Float(100.0));
        assert_eq!(col(&grouped, "END_TIME")[0], FieldValue::Null);
    }

    #[test]
    fn last_record_same_as_start_gives_a_zero_length_interval() {
        let (layer, _) = run(json!({
            "input": readings(), "start_field": "t", "id_fields": "id",
            "last_record": "same_as_start",
        }));
        assert_eq!(col(&layer, "END_TIME")[0], FieldValue::Float(30.0));
        assert_eq!(col(&layer, "END_TIME")[3], FieldValue::Float(300.0));
    }

    #[test]
    fn last_record_duration_extends_the_final_row() {
        let (layer, _) = run(json!({
            "input": readings(), "start_field": "t", "id_fields": "id",
            "last_record": "duration", "default_duration": 60.0,
        }));
        assert_eq!(col(&layer, "END_TIME")[0], FieldValue::Float(90.0));
        assert_eq!(col(&layer, "END_TIME")[3], FieldValue::Float(360.0));
    }

    #[test]
    fn iso8601_timestamps_round_trip_as_iso8601() {
        // A text timestamp column must not come back as an epoch number.
        let mut l = Layer::new("obs");
        l.add_field(FieldDef::new("t", FieldType::Text));
        for t in ["2024-03-01T00:00:00", "2024-03-01T01:30:00"] {
            l.add_feature(None, &[("t", t.into())]).unwrap();
        }
        let (layer, _) = run(json!({"input": store(l), "start_field": "t"}));
        assert_eq!(
            col(&layer, "END_TIME")[0],
            FieldValue::Text("2024-03-01T01:30:00".into())
        );
    }

    #[test]
    fn the_iso_formatter_inverts_the_parser() {
        // Round-trip across a leap day and a year boundary, which is where a
        // hand-rolled civil-calendar conversion goes wrong.
        for s in [
            "2024-02-29T23:59:59",
            "2000-01-01T00:00:00",
            "1970-01-01T00:00:00",
            "2038-12-31T12:00:00",
        ] {
            let secs = crate::find_dwell_locations::parse_iso8601_seconds(s).unwrap();
            assert_eq!(format_iso8601(secs), s, "round trip failed for {s}");
        }
    }

    #[test]
    fn a_null_start_is_skipped_and_reported_not_sorted_to_the_front() {
        // Treating a null start as epoch 0 would make it the group's first row
        // and hand the next real row a bogus interval.
        let mut l = Layer::new("obs");
        l.add_field(FieldDef::new("t", FieldType::Float));
        l.push(Feature {
            fid: 0,
            geometry: None,
            attributes: vec![FieldValue::Null],
        });
        l.add_feature(None, &[("t", 100.0f64.into())]).unwrap();
        l.add_feature(None, &[("t", 200.0f64.into())]).unwrap();
        let (layer, res) = run(json!({"input": store(l), "start_field": "t"}));
        assert_eq!(res.outputs["skipped_count"], json!(1));
        let ends = col(&layer, "END_TIME");
        assert_eq!(ends[0], FieldValue::Null, "the null row gets a null end");
        assert_eq!(ends[1], FieldValue::Float(200.0), "100 chains to 200");
    }

    #[test]
    fn an_unparseable_timestamp_is_skipped_rather_than_read_as_zero() {
        let mut l = Layer::new("obs");
        l.add_field(FieldDef::new("t", FieldType::Text));
        for t in ["not a date", "2024-01-01T00:00:00", "2024-01-02T00:00:00"] {
            l.add_feature(None, &[("t", t.into())]).unwrap();
        }
        let (_, res) = run(json!({"input": store(l), "start_field": "t"}));
        assert_eq!(res.outputs["skipped_count"], json!(1));
        assert_eq!(res.outputs["filled_count"], json!(1));
    }

    #[test]
    fn identical_starts_produce_a_zero_length_interval_rather_than_an_error() {
        let mut l = Layer::new("obs");
        l.add_field(FieldDef::new("t", FieldType::Float));
        for t in [5.0, 5.0, 9.0] {
            l.add_feature(None, &[("t", t.into())]).unwrap();
        }
        let (layer, _) = run(json!({"input": store(l), "start_field": "t"}));
        let ends = col(&layer, "END_TIME");
        assert_eq!(ends[0], FieldValue::Float(5.0), "tie -> zero-length");
        assert_eq!(ends[1], FieldValue::Float(9.0));
    }

    #[test]
    fn an_existing_end_field_is_overwritten_not_duplicated() {
        let mut l = Layer::new("obs");
        l.add_field(FieldDef::new("t", FieldType::Float));
        l.add_field(FieldDef::new("END_TIME", FieldType::Float));
        l.add_feature(None, &[("t", 1.0f64.into()), ("END_TIME", 999.0f64.into())])
            .unwrap();
        l.add_feature(None, &[("t", 2.0f64.into()), ("END_TIME", 999.0f64.into())])
            .unwrap();
        let (layer, _) = run(json!({"input": store(l), "start_field": "t"}));
        assert_eq!(layer.schema.fields().len(), 2, "no duplicate column");
        assert_eq!(col(&layer, "END_TIME")[0], FieldValue::Float(2.0));
    }

    #[test]
    fn a_fractional_end_time_is_refused_for_an_integer_column() {
        // Rounding 10.4 to 10 leaves an interval that stops short of the next
        // record, which is a silent data change rather than a formatting one.
        let mut l = Layer::new("obs");
        l.add_field(FieldDef::new("t", FieldType::Float));
        l.add_field(FieldDef::new("END_TIME", FieldType::Integer));
        for t in [1.0, 10.4] {
            l.add_feature(None, &[("t", t.into()), ("END_TIME", 0i64.into())])
                .unwrap();
        }
        let args: ToolArgs =
            serde_json::from_value(json!({"input": store(l), "start_field": "t"})).unwrap();
        let err = CalculateEndTimeTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("whole number"), "{err}");
    }

    #[test]
    fn whole_second_endpoints_are_written_to_an_integer_column() {
        let mut l = Layer::new("obs");
        l.add_field(FieldDef::new("t", FieldType::Float));
        l.add_field(FieldDef::new("END_TIME", FieldType::Integer));
        for t in [1.0, 10.0] {
            l.add_feature(None, &[("t", t.into()), ("END_TIME", 0i64.into())])
                .unwrap();
        }
        let (layer, _) = run(json!({"input": store(l), "start_field": "t"}));
        assert_eq!(col(&layer, "END_TIME")[0], FieldValue::Integer(10));
    }

    #[test]
    fn a_custom_end_field_name_is_honoured() {
        let (layer, res) = run(json!({
            "input": readings(), "start_field": "t", "id_fields": "id",
            "end_field": "stop",
        }));
        assert!(layer.schema.field_index("stop").is_some());
        assert_eq!(res.outputs["end_field"], json!("stop"));
    }

    #[test]
    fn geometry_and_attributes_are_preserved() {
        let (layer, _) = run(json!({
            "input": readings(), "start_field": "t", "id_fields": "id",
        }));
        assert!(layer.iter().all(|f| f.geometry.is_some()));
        assert_eq!(col(&layer, "id")[0], FieldValue::Text("a".into()));
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            CalculateEndTimeTool.validate(&args).is_err()
        };
        assert!(bad(json!({})));
        assert!(bad(json!({"input": "a.shp"})));
        assert!(bad(
            json!({"input": "a.shp", "start_field": "t", "last_record": "next"})
        ));
        // 'duration' without a duration would silently mean 'same_as_start'.
        assert!(bad(
            json!({"input": "a.shp", "start_field": "t", "last_record": "duration"})
        ));
        assert!(bad(json!({
            "input": "a.shp", "start_field": "t", "last_record": "duration",
            "default_duration": -5,
        })));
    }

    #[test]
    fn an_unknown_field_name_fails_loudly() {
        let args: ToolArgs =
            serde_json::from_value(json!({"input": readings(), "start_field": "nope"})).unwrap();
        let err = CalculateEndTimeTool.run(&args, &ctx()).unwrap_err();
        assert!(format!("{err}").contains("nope"), "{err}");
    }
}
