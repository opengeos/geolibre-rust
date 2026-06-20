//! GeoLibre WASI tool runner.
//!
//! A thin command-line adapter over the `whitebox_next_gen` tool registry,
//! compiled to `wasm32-wasip1` and executed in the browser through a WASI shim
//! ([`@bjorn3/browser_wasi_shim`]) with an in-memory filesystem. The JS glue
//! (`npm/tools.mjs`) preopens `/work`, writes input files there, invokes this
//! binary, and reads back any files the tool wrote.
//!
//! Invocation contract (argv after the program name):
//!
//! ```text
//! list                 -> print every tool id, one per line
//! manifests            -> print a JSON array of every tool manifest
//! manifest <id>        -> print one tool manifest as JSON
//! version              -> print the runner version
//! <tool> [--k=v ...]   -> run a tool; reads/writes files under /work via std::fs
//! ```
//!
//! Examples:
//!
//! ```text
//! geolibre list
//! geolibre slope --input=/work/dem.tif --output=/work/slope.tif --units=degrees
//! ```
//!
//! Exit codes: `0` success, `1` tool/execution error, `2` usage error.

use std::collections::{BTreeMap, HashSet};
use std::process::ExitCode;

use serde_json::Value;
use wbcore::{AllowAllCapabilities, ProgressSink, ToolArgs, ToolContext};
use wbtools_oss::{register_default_tools, ToolRegistry};

/// Builds the registry with whitebox's default (OSS) tools plus GeoLibre's own
/// tools that extend the suite.
fn build_registry() -> ToolRegistry {
    let mut registry = ToolRegistry::new();
    register_default_tools(&mut registry);
    for tool in geolibre_tools::geolibre_tools() {
        registry.register(tool);
    }
    registry
}

/// The ids of GeoLibre-authored tools, so manifests can be tagged with their
/// provenance (`"geolibre"` vs `"whitebox"`).
fn geolibre_tool_ids() -> HashSet<String> {
    geolibre_tools::geolibre_tools()
        .iter()
        .map(|t| t.metadata().id.to_string())
        .collect()
}

/// Adds a `"source"` field (`"geolibre"` or `"whitebox"`) to a manifest JSON
/// object based on its `id`. The upstream `ToolManifest` has no such field, so
/// it is injected here at serialization time.
fn tag_source(manifest: &mut Value, geolibre_ids: &HashSet<String>) {
    if let Some(obj) = manifest.as_object_mut() {
        let is_geolibre = obj
            .get("id")
            .and_then(Value::as_str)
            .map(|id| geolibre_ids.contains(id))
            .unwrap_or(false);
        obj.insert(
            "source".to_string(),
            Value::String(if is_geolibre { "geolibre" } else { "whitebox" }.to_string()),
        );
    }
}

fn main() -> ExitCode {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let Some(command) = argv.first() else {
        eprintln!("usage: geolibre <list|manifests|manifest <id>|version|<tool> [--k=v ...]>");
        return ExitCode::from(2);
    };

    match command.as_str() {
        "list" => {
            let registry = build_registry();
            for meta in registry.list() {
                println!("{}", meta.id);
            }
            ExitCode::SUCCESS
        }
        "manifests" => {
            let registry = build_registry();
            let geolibre_ids = geolibre_tool_ids();
            match serde_json::to_value(registry.manifests()) {
                Ok(mut value) => {
                    if let Some(arr) = value.as_array_mut() {
                        for m in arr.iter_mut() {
                            tag_source(m, &geolibre_ids);
                        }
                    }
                    println!("{value}");
                    ExitCode::SUCCESS
                }
                Err(e) => {
                    eprintln!("failed to serialize manifests: {e}");
                    ExitCode::from(1)
                }
            }
        }
        "manifest" => {
            let Some(id) = argv.get(1) else {
                eprintln!("usage: geolibre manifest <tool-id>");
                return ExitCode::from(2);
            };
            let registry = build_registry();
            match registry.manifests().into_iter().find(|m| &m.id == id) {
                Some(manifest) => match serde_json::to_value(&manifest) {
                    Ok(mut value) => {
                        tag_source(&mut value, &geolibre_tool_ids());
                        println!("{value}");
                        ExitCode::SUCCESS
                    }
                    Err(e) => {
                        eprintln!("failed to serialize manifest: {e}");
                        ExitCode::from(1)
                    }
                },
                None => {
                    eprintln!("tool not found: {id}");
                    ExitCode::from(1)
                }
            }
        }
        "version" | "--version" | "-V" => {
            println!("geolibre-cli {}", env!("CARGO_PKG_VERSION"));
            ExitCode::SUCCESS
        }
        tool_id => run_tool(tool_id, &argv[1..]),
    }
}

/// Runs a single tool, mapping `--key=value` CLI flags onto the JSON [`ToolArgs`]
/// the registry expects.
fn run_tool(tool_id: &str, flags: &[String]) -> ExitCode {
    let args = parse_args(flags);
    let registry = build_registry();

    let progress = StdoutProgress::default();
    let capabilities = AllowAllCapabilities;
    let ctx = ToolContext {
        progress: &progress,
        capabilities: &capabilities,
    };

    match registry.run(tool_id, &args, &ctx) {
        Ok(result) => {
            // The JS side reads results from new files under /work; emit the
            // structured outputs to stdout too so callers can inspect them.
            if let Ok(json) = serde_json::to_string(&result.outputs) {
                println!("{json}");
            }
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("{e}");
            ExitCode::from(1)
        }
    }
}

/// Parses `--key=value`, `--key value`, and bare `--flag` tokens into
/// [`ToolArgs`].
///
/// Values are type-inferred to match the JSON shapes tools expect: `true`/`false`
/// become booleans, integers and floats become numbers, and everything else
/// (including file paths) stays a string. A bare `--flag` with no value is
/// treated as `true`.
fn parse_args(flags: &[String]) -> ToolArgs {
    let mut args: BTreeMap<String, Value> = BTreeMap::new();
    let mut i = 0;
    while i < flags.len() {
        let token = &flags[i];
        let Some(stripped) = token.strip_prefix("--") else {
            i += 1;
            continue;
        };

        if let Some((key, value)) = stripped.split_once('=') {
            args.insert(key.to_string(), infer_value(value));
            i += 1;
        } else if i + 1 < flags.len() && !flags[i + 1].starts_with("--") {
            args.insert(stripped.to_string(), infer_value(&flags[i + 1]));
            i += 2;
        } else {
            // Bare flag, e.g. `--zero-background`.
            args.insert(stripped.to_string(), Value::Bool(true));
            i += 1;
        }
    }
    args
}

/// Infers a JSON value from a raw CLI string.
fn infer_value(raw: &str) -> Value {
    match raw {
        "true" => return Value::Bool(true),
        "false" => return Value::Bool(false),
        _ => {}
    }
    if let Ok(n) = raw.parse::<i64>() {
        return Value::from(n);
    }
    if let Ok(f) = raw.parse::<f64>() {
        if f.is_finite() {
            return Value::from(f);
        }
    }
    Value::String(raw.to_string())
}

/// Progress sink that forwards tool messages to stdout (captured by the WASI
/// shim and surfaced to JS). Per-percent progress is dropped to keep stdout
/// quiet; flip `EMIT_PROGRESS` if you want it.
#[derive(Default)]
struct StdoutProgress;

const EMIT_PROGRESS: bool = false;

impl ProgressSink for StdoutProgress {
    fn info(&self, msg: &str) {
        println!("{msg}");
    }

    fn progress(&self, pct: f64) {
        if EMIT_PROGRESS {
            println!("progress: {:.0}%", (pct.clamp(0.0, 1.0) * 100.0));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn infers_value_types() {
        assert_eq!(infer_value("true"), Value::Bool(true));
        assert_eq!(infer_value("42"), Value::from(42i64));
        assert_eq!(infer_value("3.5"), Value::from(3.5f64));
        assert_eq!(
            infer_value("/work/dem.tif"),
            Value::String("/work/dem.tif".to_string())
        );
    }

    #[test]
    fn parses_equals_space_and_bare_flags() {
        let flags = vec![
            "--input=/work/dem.tif".to_string(),
            "--units".to_string(),
            "degrees".to_string(),
            "--zero-background".to_string(),
        ];
        let args = parse_args(&flags);
        assert_eq!(args["input"], Value::String("/work/dem.tif".to_string()));
        assert_eq!(args["units"], Value::String("degrees".to_string()));
        assert_eq!(args["zero-background"], Value::Bool(true));
    }

    #[test]
    fn registry_lists_tools() {
        let registry = build_registry();
        assert!(!registry.list().is_empty());
    }
}
