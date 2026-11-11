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

use geolibre_cli::wbcore::{AllowAllCapabilities, ProgressSink, ToolArgs, ToolContext};
use geolibre_cli::{build_registry, manifest_json, manifests_json, ToolRegistry};
use serde_json::Value;

/// Stack for the command on native targets. Windows gives the main thread only
/// 1 MB, which a debug build can overflow while building the manifests, and some
/// tools recurse deeply on large rasters; Linux and macOS default to 8 MB.
#[cfg(not(target_family = "wasm"))]
const STACK_BYTES: usize = 64 * 1024 * 1024;

fn main() -> ExitCode {
    // WASI has no threads: run directly there.
    #[cfg(target_family = "wasm")]
    {
        run_command()
    }
    #[cfg(not(target_family = "wasm"))]
    {
        match std::thread::Builder::new()
            .stack_size(STACK_BYTES)
            .spawn(run_command)
        {
            // A panic has already printed its message (and release builds abort
            // instead of unwinding); exit like any other failure.
            Ok(handle) => handle.join().unwrap_or_else(|_| ExitCode::from(1)),
            // No thread available: fall back to the main thread's stack.
            Err(_) => run_command(),
        }
    }
}

fn run_command() -> ExitCode {
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
            // The I/O-enriched manifests (each param's `schema`, `io_role` and
            // `data_kind`), so host UIs can route raster/vector/lidar inputs and
            // render widgets without a separate catalog.
            let registry = build_registry();
            println!("{}", Value::Array(manifests_json(&registry)));
            ExitCode::SUCCESS
        }
        "manifest" => {
            let Some(id) = argv.get(1) else {
                eprintln!("usage: geolibre manifest <tool-id>");
                return ExitCode::from(2);
            };
            match manifest_json(&build_registry(), id) {
                Some(value) => {
                    println!("{value}");
                    ExitCode::SUCCESS
                }
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
    let mut args = parse_args(flags);
    let registry = build_registry();
    decode_json_args(&mut args, || string_default_keys(&registry, tool_id));

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
/// (including file paths) stays a string. JSON lists and objects are decoded
/// afterwards by [`decode_json_args`], which needs the tool's manifest. A bare `--flag` with no value is
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

/// Turns `[...]` and `{...}` values that parse as JSON into that JSON, so list
/// and object parameters (`--scales=[8,16]`, `--background_clr=[0,0,0,255]`)
/// reach tools as arrays and objects rather than strings.
///
/// A parameter whose manifest default is a string keeps its text, because such
/// tools parse the JSON themselves (e.g. `fit_variogram --lags_json=[...]`).
/// `string_keys` builds that set and runs only when a value looks like JSON, so
/// ordinary runs don't pay for building the manifests.
///
/// Args:
///     args: Parsed CLI arguments, updated in place.
///     string_keys: Returns the names of parameters whose default is a string.
fn decode_json_args<F>(args: &mut ToolArgs, string_keys: F)
where
    F: FnOnce() -> HashSet<String>,
{
    let mut string_keys = Some(string_keys);
    let mut keep_text: Option<HashSet<String>> = None;
    for (key, value) in args.iter_mut() {
        let Value::String(text) = value else { continue };
        let Some(parsed) = parse_json_text(text) else {
            continue;
        };
        let keep =
            keep_text.get_or_insert_with(|| string_keys.take().map(|f| f()).unwrap_or_default());
        if !keep.contains(key) {
            *value = parsed;
        }
    }
}

/// Parses `raw` as a JSON array or object, or returns `None` when it is not
/// one (paths, expressions and other free text stay strings).
fn parse_json_text(raw: &str) -> Option<Value> {
    let trimmed = raw.trim_start();
    if !(trimmed.starts_with('[') || trimmed.starts_with('{')) {
        return None;
    }
    match serde_json::from_str::<Value>(raw) {
        Ok(v @ (Value::Array(_) | Value::Object(_))) => Some(v),
        _ => None,
    }
}

/// Names of `tool_id`'s parameters whose manifest default is a string.
fn string_default_keys(registry: &ToolRegistry, tool_id: &str) -> HashSet<String> {
    registry
        .manifests()
        .into_iter()
        .find(|m| m.id == tool_id)
        .map(|m| {
            m.defaults
                .into_iter()
                .filter(|(_, v)| v.is_string())
                .map(|(k, _)| k)
                .collect()
        })
        .unwrap_or_default()
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
    fn decodes_json_lists_and_objects() {
        let flags = vec![
            "--scales=[8,16]".to_string(),
            "--background_clr=[0,0,0,255]".to_string(),
            "--library".to_string(),
            r#"{"name":"water"}"#.to_string(),
            "--input=/work/[dem].tif".to_string(),
            "--expression={dem} * 2".to_string(),
        ];
        let mut args = parse_args(&flags);
        decode_json_args(&mut args, HashSet::new);
        assert_eq!(args["scales"], serde_json::json!([8, 16]));
        assert_eq!(args["background_clr"], serde_json::json!([0, 0, 0, 255]));
        assert_eq!(args["library"], serde_json::json!({"name": "water"}));
        // Not JSON: left as text.
        assert_eq!(args["input"], Value::from("/work/[dem].tif"));
        assert_eq!(args["expression"], Value::from("{dem} * 2"));
    }

    #[test]
    fn keeps_json_text_for_string_parameters() {
        let mut args = parse_args(&["--lags_json=[{\"lag\":1}]".to_string()]);
        decode_json_args(&mut args, || HashSet::from(["lags_json".to_string()]));
        assert_eq!(args["lags_json"], Value::from(r#"[{"lag":1}]"#));
    }

    #[test]
    fn skips_manifest_lookup_without_json_values() {
        let mut args = parse_args(&["--input=dem.tif".to_string(), "--z=2".to_string()]);
        decode_json_args(&mut args, || panic!("manifests built for plain args"));
        assert_eq!(args["input"], Value::from("dem.tif"));
    }

    #[test]
    fn string_default_keys_come_from_the_manifest() {
        let registry = build_registry();
        let keys = string_default_keys(&registry, "fit_variogram");
        assert!(keys.contains("lags_json"), "{keys:?}");
        assert!(!string_default_keys(&registry, "frangi_filter").contains("scales"));
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
}
