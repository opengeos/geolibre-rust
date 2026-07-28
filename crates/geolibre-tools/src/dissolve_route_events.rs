//! GeoLibre tool: merge adjacent or overlapping linear-referencing events.
//!
//! Pure-Rust counterpart of ArcGIS Pro's *Dissolve Route Events* (Linear
//! Referencing). GeoLibre has built out most of that family — `create_routes`,
//! `locate_lines_along_routes`, `transform_route_events` — on top of the
//! bundled `route_event_overlay` and `locate_points_along_routes`. The dissolve
//! was the missing cleanup step.
//!
//! It matters because it is what you run immediately after an overlay:
//! overlaying two event tables fragments every span at each boundary, so the
//! result is a shattered table full of adjacent rows with identical
//! attributes. Dissolving reassembles them into the smallest set of spans that
//! carries the same information.
//!
//! Two modes:
//!   * `dissolve` — only merge touching/overlapping spans whose dissolve-field
//!     values all match.
//!   * `concatenate` — merge any touching/overlapping spans on a route and join
//!     the distinct dissolve-field values with a separator.
//!
//! Point events (`from == to`) fall out of the same sweep correctly.

use std::collections::BTreeMap;
use std::collections::HashMap;

use serde_json::{json, Value};
use wbcore::{
    LicenseTier, Tool, ToolArgs, ToolCategory, ToolContext, ToolError, ToolMetadata, ToolParamSpec,
    ToolRunResult,
};
use wbvector::{FieldDef, FieldType, FieldValue, Layer};

use crate::vector_common::{load_input_layer, parse_optional_str, write_or_store_layer};

pub struct DissolveRouteEventsTool;

impl Tool for DissolveRouteEventsTool {
    fn metadata(&self) -> ToolMetadata {
        ToolMetadata {
            id: "dissolve_route_events",
            display_name: "Dissolve Route Events",
            summary: "Merge route events that touch or overlap along the same route into single measure spans, either only when the dissolve fields match ('dissolve') or by concatenating differing values ('concatenate'). Like ArcGIS Dissolve Route Events.",
            category: ToolCategory::Vector,
            license_tier: LicenseTier::Open,
            params: vec![
                ToolParamSpec {
                    name: "input",
                    description: "Input event table (line events with from/to measures).",
                    required: true,
                },
                ToolParamSpec {
                    name: "output",
                    description: "Output event table path. If omitted, stored in memory.",
                    required: false,
                },
                ToolParamSpec {
                    name: "route_id_field",
                    description: "Field holding the route identifier.",
                    required: true,
                },
                ToolParamSpec {
                    name: "from_measure_field",
                    description: "Field holding each event's from-measure.",
                    required: true,
                },
                ToolParamSpec {
                    name: "to_measure_field",
                    description: "Field holding each event's to-measure. For point events, use the same field as 'from_measure_field'.",
                    required: true,
                },
                ToolParamSpec {
                    name: "dissolve_fields",
                    description: "Comma/semicolon-separated field name(s) that must agree for two spans to merge (or whose values are concatenated).",
                    required: true,
                },
                ToolParamSpec {
                    name: "mode",
                    description: "'dissolve' (default) or 'concatenate'.",
                    required: false,
                },
                ToolParamSpec {
                    name: "tolerance",
                    description: "Measure gap below which two spans count as adjacent (default 0).",
                    required: false,
                },
                ToolParamSpec {
                    name: "separator",
                    description: "Separator used to join values in 'concatenate' mode (default ', ').",
                    required: false,
                },
            ],
        }
    }

    fn validate(&self, args: &ToolArgs) -> Result<(), ToolError> {
        require_str(args, "input")?;
        require_str(args, "route_id_field")?;
        require_str(args, "from_measure_field")?;
        require_str(args, "to_measure_field")?;
        if split_list(require_str(args, "dissolve_fields")?).is_empty() {
            return Err(ToolError::Validation(
                "'dissolve_fields' must name at least one field".to_string(),
            ));
        }
        parse_mode(args)?;
        let tol = parse_optional_f64(args, "tolerance")?.unwrap_or(0.0);
        if tol < 0.0 {
            return Err(ToolError::Validation(
                "'tolerance' must be non-negative".to_string(),
            ));
        }
        Ok(())
    }

    fn run(&self, args: &ToolArgs, ctx: &ToolContext) -> Result<ToolRunResult, ToolError> {
        let input = require_str(args, "input")?;
        let output = parse_optional_str(args, "output")?;
        let route_name = require_str(args, "route_id_field")?;
        let from_name = require_str(args, "from_measure_field")?;
        let to_name = require_str(args, "to_measure_field")?;
        let diss_names = split_list(require_str(args, "dissolve_fields")?);
        let mode = parse_mode(args)?;
        let tol = parse_optional_f64(args, "tolerance")?.unwrap_or(0.0);
        let sep = parse_optional_str(args, "separator")?.unwrap_or(", ");

        let layer = load_input_layer(input)?;
        if layer.features.is_empty() {
            return Err(ToolError::Execution("input has no features".to_string()));
        }
        let idx = |name: &str| -> Result<usize, ToolError> {
            layer
                .schema
                .field_index(name)
                .ok_or_else(|| ToolError::Validation(format!("field '{name}' not found in input")))
        };
        let r_i = idx(route_name)?;
        let f_i = idx(from_name)?;
        let t_i = idx(to_name)?;
        let d_i: Vec<usize> = diss_names
            .iter()
            .map(|n| idx(n))
            .collect::<Result<_, _>>()?;

        // Collect events, keyed by route (plus the dissolve tuple in 'dissolve'
        // mode so non-matching spans never merge).
        let mut groups: Vec<Group> = Vec::new();
        let mut pos: HashMap<String, usize> = HashMap::new();
        let mut route_keys: std::collections::HashSet<String> = std::collections::HashSet::new();
        let mut skipped = 0usize;
        for feat in &layer.features {
            let (Some(from), Some(to)) = (
                feat.attributes.get(f_i).and_then(FieldValue::as_f64),
                feat.attributes.get(t_i).and_then(FieldValue::as_f64),
            ) else {
                skipped += 1;
                continue;
            };
            if !from.is_finite() || !to.is_finite() {
                skipped += 1;
                continue;
            }
            // Tolerate reversed measures rather than dropping the row.
            let (lo, hi) = if from <= to { (from, to) } else { (to, from) };
            let route = key_of(feat.attributes.get(r_i));
            route_keys.insert(route.clone());
            let vals: Vec<String> = d_i
                .iter()
                .map(|&i| key_of(feat.attributes.get(i)))
                .collect();
            let gkey = match mode {
                Mode::Dissolve => format!("{route}\u{1f}{}", vals.join("\u{1f}")),
                Mode::Concatenate => route.clone(),
            };
            let g = *pos.entry(gkey).or_insert_with(|| {
                groups.push(Group {
                    route: feat.attributes.get(r_i).cloned().unwrap_or(FieldValue::Null),
                    spans: Vec::new(),
                });
                groups.len() - 1
            });
            groups[g].spans.push(Span {
                lo,
                hi,
                vals,
                count: 1,
            });
        }

        let n_in = layer.features.len() - skipped;
        if n_in == 0 {
            return Err(ToolError::Execution(
                "no input rows had usable from/to measures".to_string(),
            ));
        }
        ctx.progress
            .info(&format!("dissolving {n_in} event(s) on {} group(s)", groups.len()));

        // Sweep-merge each group.
        let mut out = Layer::new("dissolved_events");
        if let Some(epsg) = layer.crs_epsg() {
            out = out.with_crs_epsg(epsg);
        }
        out.add_field(FieldDef::new(
            route_name,
            layer.schema.fields()[r_i].field_type,
        ));
        out.add_field(FieldDef::new(from_name, FieldType::Float));
        out.add_field(FieldDef::new(to_name, FieldType::Float));
        // Values travel through the sweep as strings (they are also the group
        // key), so both modes emit Text. Declaring the source type here while
        // writing Text would leave the schema and the values disagreeing, which
        // downstream writers either reject or silently coerce.
        for name in &diss_names {
            out.add_field(FieldDef::new(name.as_str(), FieldType::Text));
        }
        out.add_field(FieldDef::new("EVENT_COUNT", FieldType::Integer));

        let mut n_out = 0usize;
        for group in &mut groups {
            group.spans.sort_by(|a, b| a.lo.total_cmp(&b.lo));
            let mut merged: Vec<Span> = Vec::new();
            for span in group.spans.drain(..) {
                match merged.last_mut() {
                    // Adjacent (within tolerance) or overlapping -> extend.
                    Some(cur) if span.lo <= cur.hi + tol => {
                        cur.hi = cur.hi.max(span.hi);
                        cur.count += span.count;
                        for (i, v) in span.vals.into_iter().enumerate() {
                            if !cur.vals[i].split(sep).any(|p| p == v) {
                                cur.vals[i].push_str(sep);
                                cur.vals[i].push_str(&v);
                            }
                        }
                    }
                    _ => merged.push(span),
                }
            }
            for span in merged {
                let mut attrs: Vec<(&str, FieldValue)> = vec![
                    (route_name, group.route.clone()),
                    (from_name, FieldValue::Float(span.lo)),
                    (to_name, FieldValue::Float(span.hi)),
                ];
                for (name, v) in diss_names.iter().zip(span.vals.iter()) {
                    attrs.push((name.as_str(), FieldValue::Text(v.clone())));
                }
                attrs.push(("EVENT_COUNT", FieldValue::Integer(span.count)));
                out.add_feature(None, &attrs)
                    .map_err(|e| ToolError::Execution(format!("failed adding row: {e}")))?;
                n_out += 1;
            }
        }

        let out_path = write_or_store_layer(out, output)?;
        let mut outputs = BTreeMap::new();
        outputs.insert("output".to_string(), json!(out_path));
        outputs.insert("input_event_count".to_string(), json!(n_in));
        outputs.insert("output_event_count".to_string(), json!(n_out));
        outputs.insert("group_count".to_string(), json!(groups.len()));
        // In 'dissolve' mode a group is (route x dissolve tuple), so the group
        // count is not the route count; report the distinct routes separately.
        outputs.insert("route_count".to_string(), json!(route_keys.len()));
        outputs.insert("skipped_rows".to_string(), json!(skipped));
        Ok(ToolRunResult { outputs })
    }
}

struct Group {
    route: FieldValue,
    spans: Vec<Span>,
}

struct Span {
    lo: f64,
    hi: f64,
    vals: Vec<String>,
    count: i64,
}

// ── Params ──────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    Dissolve,
    Concatenate,
}

fn parse_mode(args: &ToolArgs) -> Result<Mode, ToolError> {
    match args
        .get("mode")
        .and_then(Value::as_str)
        .map(str::trim)
        .map(str::to_ascii_lowercase)
        .as_deref()
    {
        None | Some("") | Some("dissolve") => Ok(Mode::Dissolve),
        Some("concatenate") => Ok(Mode::Concatenate),
        Some(o) => Err(ToolError::Validation(format!(
            "'mode' must be 'dissolve' or 'concatenate', got '{o}'"
        ))),
    }
}

fn key_of(v: Option<&FieldValue>) -> String {
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

fn split_list(s: &str) -> Vec<String> {
    s.split([',', ';'])
        .map(str::trim)
        .filter(|p| !p.is_empty())
        .map(str::to_string)
        .collect()
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

    /// Event table rows: (route, from, to, surface).
    fn events(rows: &[(&str, f64, f64, &str)]) -> String {
        let mut l = Layer::new("events");
        l.add_field(FieldDef::new("RID", FieldType::Text));
        l.add_field(FieldDef::new("FMEAS", FieldType::Float));
        l.add_field(FieldDef::new("TMEAS", FieldType::Float));
        l.add_field(FieldDef::new("SURF", FieldType::Text));
        for (rid, f, t, s) in rows {
            l.add_feature(
                None,
                &[
                    ("RID", FieldValue::Text((*rid).to_string())),
                    ("FMEAS", FieldValue::Float(*f)),
                    ("TMEAS", FieldValue::Float(*t)),
                    ("SURF", FieldValue::Text((*s).to_string())),
                ],
            )
            .unwrap();
        }
        let id = wbvector::memory_store::put_vector(l);
        wbvector::memory_store::make_vector_memory_path(&id)
    }

    fn base(input: &str) -> serde_json::Value {
        json!({
            "input": input, "route_id_field": "RID", "from_measure_field": "FMEAS",
            "to_measure_field": "TMEAS", "dissolve_fields": "SURF"
        })
    }

    fn run(args: serde_json::Value) -> (ToolRunResult, Layer) {
        let args: ToolArgs = serde_json::from_value(args).unwrap();
        let out = DissolveRouteEventsTool.run(&args, &ctx()).unwrap();
        let layer = load_input_layer(out.outputs["output"].as_str().unwrap()).unwrap();
        (out, layer)
    }

    fn spans(l: &Layer) -> Vec<(String, f64, f64)> {
        let (r, f, t) = (
            l.schema.field_index("RID").unwrap(),
            l.schema.field_index("FMEAS").unwrap(),
            l.schema.field_index("TMEAS").unwrap(),
        );
        l.iter()
            .map(|x| {
                (
                    x.attributes[r].as_str().unwrap().to_string(),
                    x.attributes[f].as_f64().unwrap(),
                    x.attributes[t].as_f64().unwrap(),
                )
            })
            .collect()
    }

    #[test]
    fn adjacent_matching_spans_merge() {
        // The classic post-overlay fragmentation: 0-1, 1-2, 2-3 all asphalt.
        let input = events(&[
            ("A", 0.0, 1.0, "asphalt"),
            ("A", 1.0, 2.0, "asphalt"),
            ("A", 2.0, 3.0, "asphalt"),
        ]);
        let (out, layer) = run(base(&input));
        assert_eq!(out.outputs["output_event_count"], json!(1));
        assert_eq!(spans(&layer), vec![("A".to_string(), 0.0, 3.0)]);
        let c = layer.schema.field_index("EVENT_COUNT").unwrap();
        assert_eq!(layer.features[0].attributes[c].as_i64(), Some(3));
    }

    #[test]
    fn differing_attribute_blocks_the_merge() {
        let input = events(&[
            ("A", 0.0, 1.0, "asphalt"),
            ("A", 1.0, 2.0, "gravel"),
            ("A", 2.0, 3.0, "asphalt"),
        ]);
        let (out, layer) = run(base(&input));
        assert_eq!(out.outputs["output_event_count"], json!(3));
        // The two asphalt spans are non-adjacent, so they stay separate.
        let mut got = spans(&layer);
        got.sort_by(|a, b| a.1.total_cmp(&b.1));
        assert_eq!(got[0], ("A".to_string(), 0.0, 1.0));
    }

    #[test]
    fn different_routes_never_merge() {
        let input = events(&[("A", 0.0, 1.0, "asphalt"), ("B", 1.0, 2.0, "asphalt")]);
        let (out, _l) = run(base(&input));
        assert_eq!(out.outputs["output_event_count"], json!(2));
        assert_eq!(out.outputs["route_count"], json!(2));
    }

    #[test]
    fn overlapping_spans_merge_to_their_union() {
        let input = events(&[("A", 0.0, 5.0, "x"), ("A", 3.0, 9.0, "x")]);
        let (_o, layer) = run(base(&input));
        assert_eq!(spans(&layer), vec![("A".to_string(), 0.0, 9.0)]);
    }

    #[test]
    fn tolerance_bridges_a_small_gap() {
        let input = events(&[("A", 0.0, 1.0, "x"), ("A", 1.05, 2.0, "x")]);
        let (tight, _l) = run(base(&input));
        assert_eq!(tight.outputs["output_event_count"], json!(2));
        let mut with_tol = base(&input);
        with_tol["tolerance"] = json!(0.1);
        let (loose, _l) = run(with_tol);
        assert_eq!(loose.outputs["output_event_count"], json!(1));
    }

    #[test]
    fn concatenate_mode_joins_differing_values() {
        let input = events(&[("A", 0.0, 1.0, "asphalt"), ("A", 1.0, 2.0, "gravel")]);
        let mut args = base(&input);
        args["mode"] = json!("concatenate");
        let (out, layer) = run(args);
        assert_eq!(out.outputs["output_event_count"], json!(1));
        let s = layer.schema.field_index("SURF").unwrap();
        assert_eq!(
            layer.features[0].attributes[s].as_str(),
            Some("asphalt, gravel")
        );
    }

    #[test]
    fn point_events_collapse_when_coincident() {
        // from == to; two at the same measure merge, a distant one does not.
        let input = events(&[("A", 5.0, 5.0, "x"), ("A", 5.0, 5.0, "x"), ("A", 9.0, 9.0, "x")]);
        let (out, _l) = run(base(&input));
        assert_eq!(out.outputs["output_event_count"], json!(2));
    }

    #[test]
    fn rows_without_measures_are_skipped_not_fatal() {
        let mut l = Layer::new("events");
        l.add_field(FieldDef::new("RID", FieldType::Text));
        l.add_field(FieldDef::new("FMEAS", FieldType::Float));
        l.add_field(FieldDef::new("TMEAS", FieldType::Float));
        l.add_field(FieldDef::new("SURF", FieldType::Text));
        l.add_feature(
            None,
            &[
                ("RID", FieldValue::Text("A".into())),
                ("FMEAS", FieldValue::Null),
                ("TMEAS", FieldValue::Float(1.0)),
                ("SURF", FieldValue::Text("x".into())),
            ],
        )
        .unwrap();
        l.add_feature(
            None,
            &[
                ("RID", FieldValue::Text("A".into())),
                ("FMEAS", FieldValue::Float(0.0)),
                ("TMEAS", FieldValue::Float(1.0)),
                ("SURF", FieldValue::Text("x".into())),
            ],
        )
        .unwrap();
        let id = wbvector::memory_store::put_vector(l);
        let input = wbvector::memory_store::make_vector_memory_path(&id);
        let (out, _l) = run(base(&input));
        assert_eq!(out.outputs["skipped_rows"], json!(1));
        assert_eq!(out.outputs["output_event_count"], json!(1));
    }

    #[test]
    fn rejects_bad_parameters() {
        let bad = |v: serde_json::Value| {
            let args: ToolArgs = serde_json::from_value(v).unwrap();
            DissolveRouteEventsTool.validate(&args)
        };
        assert!(bad(json!({})).is_err());
        let mut ok = base("e.csv");
        assert!(bad(ok.clone()).is_ok());
        ok["mode"] = json!("bogus");
        assert!(bad(ok.clone()).is_err());
        ok["mode"] = json!("dissolve");
        ok["tolerance"] = json!(-1);
        assert!(bad(ok).is_err());
    }
}
