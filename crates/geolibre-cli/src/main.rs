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
use wbcore::{
    AllowAllCapabilities, ProgressSink, ToolArgs, ToolContext, ToolManifest, ToolParamDescriptor,
    ToolParamSchema,
};
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

/// Parameters for the two tools that declare none anywhere -- not in their
/// manifest, not in their metadata, and not in a schema table -- so there is
/// nothing to backfill them from.
///
/// Both are the hillshade pair in `wbtools_oss`'
/// `geomorphometry::basic_terrain_tools`: `hillshade_metadata`/`hillshade_manifest`
/// (and the `multidirectional_` twins) are written with `params: vec![]`, while
/// `run_shade_core` reads `input`, `output`, `z_factor`, `altitude` and `azimuth`
/// (plus `full_360_mode` for the multidirectional variant). Names and optionality
/// below mirror that function: `input` is the only required one
/// (`parse_raster_path_arg(args, "input")`), the rest have defaults.
///
/// Keep this list as short as the upstream gap requires -- a tool with a schema
/// table or metadata params needs no entry.
const HILLSHADE_PARAMS: &[(&str, &str, bool)] = &[
    ("input", "Input DEM raster", true),
    (
        "z_factor",
        "Vertical exaggeration applied to elevations",
        false,
    ),
    (
        "altitude",
        "Illumination altitude above the horizon in degrees",
        false,
    ),
];

fn fallback_params(tool_id: &str) -> Option<Vec<ToolParamDescriptor>> {
    let extra: &[(&str, &str, bool)] = match tool_id {
        "hillshade" => &[(
            "azimuth",
            "Illumination azimuth in degrees clockwise from north",
            false,
        )],
        "multidirectional_hillshade" => &[(
            "full_360_mode",
            "Combine 8 light sources over the full 360 degrees instead of the default 4",
            false,
        )],
        _ => return None,
    };
    Some(
        HILLSHADE_PARAMS
            .iter()
            .chain(extra)
            .chain(std::iter::once(&("output", "Output raster path", false)))
            .map(|(name, description, required)| ToolParamDescriptor {
                name: (*name).to_string(),
                description: (*description).to_string(),
                required: *required,
            })
            .collect(),
    )
}

/// Tools whose schema table disagrees with what their implementation reads, so
/// the backfill takes their `metadata()` params instead.
///
/// The table is right for most tools (see [`schema_table_params`]), but for
/// these seven it is stale while the metadata matches the code:
///
/// - `downslope_flowpath_length`, `unnest_basins` -- read a pointer raster via
///   `parse_pointer_input` (`d8_pntr`), not the `dem` the table lists, and
///   `unnest_basins` also needs `pour_points`.
/// - `num_inflowing_neighbours` -- reads a DEM via `parse_dem_and_output`, not
///   the table's `d8_pntr`.
/// - `longest_flowpath`, `max_upslope_value` -- need a second raster the table
///   omits (`basins` and `values` respectively).
/// - `hypsometric_analysis`, `slope_vs_elev_plot` -- take a raster *list*
///   (`inputs`, plural) plus `watershed`; the table lists neither.
///
/// Each was found by running the tool with the params its manifest advertised
/// and catching `missing required parameter`; `manifest_params_satisfy_the_runner`
/// re-runs them so the list cannot rot silently.
const PREFER_METADATA_PARAMS: &[&str] = &[
    "downslope_flowpath_length",
    "hypsometric_analysis",
    "longest_flowpath",
    "max_upslope_value",
    "num_inflowing_neighbours",
    "slope_vs_elev_plot",
    "unnest_basins",
];

/// Sort key that renders a form in the order a user fills it: inputs, then
/// options, then outputs. A schema table is a `BTreeMap`, so without this the
/// params come out alphabetically and `output` lands in the middle of the
/// options.
fn param_order(schema: &ToolParamSchema) -> u8 {
    match schema {
        ToolParamSchema::Input(_) => 0,
        ToolParamSchema::Output(_) => 2,
        _ => 1,
    }
}

/// The params a tool declares in its schema table, if it has one.
///
/// `wbtools_oss::tools::tool_param_schemas` is the same per-tool table the
/// runner's own arg parsing follows, and it is more reliable than `metadata()`
/// for the tools this backfill targets: 30 of them are generated by a shared
/// macro (`create_stream_tool_impl!`) whose metadata hands every tool the same
/// three placeholder params. `extract_valleys`, for one, is declared
/// `d8_pntr`/`streams_raster`/`output` by that macro while its implementation
/// reads `dem`, `line_thin`, `variant`, `filter_size` and `output` -- which is
/// exactly what the schema table says.
///
/// Descriptions and the required flag come from the companion tables
/// (`tool_param_descriptions`/`tool_param_required`, both sourced from
/// whitebox's generated param docs), with the tool's own metadata as a fallback
/// for a param the docs do not cover.
fn schema_table_params(
    tool_id: &str,
    metadata_params: &[ToolParamDescriptor],
) -> Option<Vec<ToolParamDescriptor>> {
    let schemas = wbtools_oss::tools::tool_param_schemas(tool_id)?;
    if schemas.is_empty() {
        return None;
    }
    let descriptions = wbtools_oss::tools::tool_param_descriptions(tool_id).unwrap_or_default();
    let required = wbtools_oss::tools::tool_param_required(tool_id).unwrap_or_default();
    let by_name: BTreeMap<&str, &ToolParamDescriptor> = metadata_params
        .iter()
        .map(|param| (param.name.as_str(), param))
        .collect();

    let mut params: Vec<(u8, ToolParamDescriptor)> = schemas
        .iter()
        .map(|(name, schema)| {
            let from_metadata = by_name.get(name.as_str());
            let descriptor = ToolParamDescriptor {
                name: name.clone(),
                description: descriptions
                    .get(name)
                    .cloned()
                    .or_else(|| from_metadata.map(|param| param.description.clone()))
                    .unwrap_or_default(),
                required: required
                    .get(name)
                    .copied()
                    .or_else(|| from_metadata.map(|param| param.required))
                    .unwrap_or(false),
            };
            (param_order(schema), descriptor)
        })
        .collect();
    params.sort_by_key(|(order, _)| *order);
    Some(params.into_iter().map(|(_, param)| param).collect())
}

/// Every tool manifest, with an empty `params` list backfilled from whatever the
/// tool does declare elsewhere.
///
/// A tool that does not override [`wbcore::Tool::manifest`] gets the trait
/// default, which derives its params from `metadata()`. A tool that *does*
/// override it hand-writes the whole `ToolManifest`, and 138 of them (every
/// Hydrology tool -- `d8_pointer`, `fill_depressions`, `basins`, ... -- plus
/// `aspect`, `sky_view_factor` and more) write `params: vec![]` while their
/// schema table and `metadata()` still describe the real parameters and the
/// runner still requires them. The manifest then says "this tool takes no
/// parameters" and the engine answers `validation error: missing required
/// parameter '...'`, so a host that builds its form or its CLI args from the
/// manifest cannot run the tool at all (opengeos/geolibre-rust#327; GeoLibre
/// rendered "This tool has no parameters." for all 138).
///
/// Sources, in order: the tool's schema table ([`schema_table_params`]), its
/// metadata, then [`fallback_params`]. An empty list is treated as missing
/// metadata rather than as a parameterless tool -- none of the 138 is genuinely
/// parameterless. A manifest that declares params is left untouched, so a
/// hand-written list still wins over all three.
fn manifests_with_metadata_params(registry: &ToolRegistry) -> Vec<ToolManifest> {
    let metadata_params: BTreeMap<String, Vec<ToolParamDescriptor>> = registry
        .list()
        .into_iter()
        .map(|meta| {
            (
                meta.id.to_string(),
                meta.params.iter().map(ToolParamDescriptor::from).collect(),
            )
        })
        .collect();
    registry
        .manifests()
        .into_iter()
        .map(|mut manifest| {
            if manifest.params.is_empty() {
                let from_metadata = metadata_params
                    .get(&manifest.id)
                    .map(Vec::as_slice)
                    .unwrap_or_default();
                let from_schema_table = if PREFER_METADATA_PARAMS.contains(&manifest.id.as_str()) {
                    None
                } else {
                    schema_table_params(&manifest.id, from_metadata)
                };
                if let Some(params) = from_schema_table
                    .or_else(|| Some(from_metadata.to_vec()).filter(|p| !p.is_empty()))
                    .or_else(|| fallback_params(&manifest.id))
                {
                    manifest.params = params;
                }
            }
            manifest
        })
        .collect()
}

/// Ids whose manifest declares no params, so [`manifests_with_metadata_params`]
/// has to reconstruct them.
fn backfilled_tool_ids(registry: &ToolRegistry) -> HashSet<String> {
    registry
        .manifests()
        .into_iter()
        .filter(|manifest| manifest.params.is_empty())
        .map(|manifest| manifest.id)
        .collect()
}

/// Serializes a manifest with per-param I/O schema. GeoLibre-authored tools use
/// their explicit schemas ([`geolibre_tools::geolibre_param_schemas`]); a
/// whitebox tool whose params were backfilled uses its own schema table
/// ([`wbtools_oss::tools::tool_param_schemas`]), since that table is where those
/// params came from and it types them exactly; everything else keeps wbcore's
/// name/description-based inference.
///
/// The table is deliberately *not* applied to a tool that shipped its own
/// params. It types 291 of them differently from the inference, and while most
/// of those differences are corrections (`bilateral_filter`'s `sigma_int` is a
/// number, not the raster input the inference guessed), a few are losses:
/// `global_morans_i`'s `output_html`/`output_csv`/`output_json` are real output
/// paths that the table calls plain strings, and a host would stop naming and
/// collecting those files. Re-typing tools that already work is a separate
/// change from making 138 broken ones work, and it needs its own verification.
fn enriched_manifest(manifest: &wbcore::ToolManifest, backfilled: &HashSet<String>) -> Value {
    match geolibre_tools::geolibre_param_schemas(&manifest.id).or_else(|| {
        backfilled
            .contains(&manifest.id)
            .then(|| wbtools_oss::tools::tool_param_schemas(&manifest.id))
            .flatten()
    }) {
        Some(schemas) => wbcore::manifest_with_param_schema_json(manifest, &schemas),
        None => wbcore::manifest_with_io_schema_json(manifest),
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
            let backfilled = backfilled_tool_ids(&registry);
            // Emit the I/O-enriched manifest (each param's `schema`, `io_role`,
            // and `data_kind`) rather than the bare manifest, so host UIs can
            // route raster/vector/lidar inputs and render widgets without a
            // separate catalog. Mirrors what the ArcGIS catalog snapshot carries.
            let arr: Vec<Value> = manifests_with_metadata_params(&registry)
                .iter()
                .map(|m| {
                    let mut value = enriched_manifest(m, &backfilled);
                    tag_source(&mut value, &geolibre_ids);
                    value
                })
                .collect();
            println!("{}", Value::Array(arr));
            ExitCode::SUCCESS
        }
        "manifest" => {
            let Some(id) = argv.get(1) else {
                eprintln!("usage: geolibre manifest <tool-id>");
                return ExitCode::from(2);
            };
            let registry = build_registry();
            let backfilled = backfilled_tool_ids(&registry);
            match manifests_with_metadata_params(&registry)
                .into_iter()
                .find(|m| &m.id == id)
            {
                Some(manifest) => {
                    let mut value = enriched_manifest(&manifest, &backfilled);
                    tag_source(&mut value, &geolibre_tool_ids());
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

    #[test]
    fn no_manifest_is_empty_while_its_metadata_has_params() {
        // opengeos/geolibre-rust#327: 138 tools hand-write a `ToolManifest` with
        // `params: vec![]` while `metadata()` declares the real ones and
        // `validate()` requires them, so `manifests` advertised parameterless
        // tools that no host could run. A tool that declares params anywhere must
        // declare them in its manifest.
        //
        // Only the *empty* case is checked. A manifest that lists params is left
        // as written, deliberate omissions included -- e.g.
        // terrain_corrected_optical_analytics keeps a legacy `output` in its
        // metadata ("Deprecated -- use output_prefix. Ignored.") that its
        // manifest drops on purpose.
        let registry = build_registry();
        let manifests = manifests_with_metadata_params(&registry);
        let by_id: BTreeMap<_, _> = manifests.iter().map(|m| (m.id.clone(), m)).collect();

        let mut empty: Vec<&str> = Vec::new();
        for meta in registry.list() {
            if meta.params.is_empty() {
                continue;
            }
            let manifest = by_id
                .get(meta.id)
                .unwrap_or_else(|| panic!("manifest for {} present", meta.id));
            if manifest.params.is_empty() {
                empty.push(meta.id);
            }
        }
        assert!(
            empty.is_empty(),
            "manifests report no params for tools whose metadata declares them: {empty:?}"
        );
    }

    #[test]
    fn no_tool_advertises_itself_as_parameterless() {
        // Stronger than the check above: after the metadata backfill, every tool
        // in the registry declares at least one parameter -- none of them takes
        // zero. A new tool landing here with an empty list is almost certainly
        // the #327 bug again; give it params in its manifest (or, if both its
        // manifest and metadata are empty upstream, an entry in fallback_params).
        let registry = build_registry();
        let empty: Vec<String> = manifests_with_metadata_params(&registry)
            .into_iter()
            .filter(|m| m.params.is_empty())
            .map(|m| m.id)
            .collect();
        assert!(empty.is_empty(), "tools declaring no params: {empty:?}");
    }

    #[test]
    fn the_hillshade_pair_falls_back_to_curated_params() {
        // Neither hillshade nor multidirectional_hillshade declares params in its
        // manifest *or* its metadata, so the backfill has nothing to copy; these
        // come from fallback_params, which mirrors upstream's run_shade_core.
        let registry = build_registry();
        let manifests = manifests_with_metadata_params(&registry);
        for (id, extra) in [
            ("hillshade", "azimuth"),
            ("multidirectional_hillshade", "full_360_mode"),
        ] {
            let manifest = manifests
                .iter()
                .find(|m| m.id == id)
                .unwrap_or_else(|| panic!("{id} present"));
            let names: Vec<&str> = manifest.params.iter().map(|p| p.name.as_str()).collect();
            assert!(names.contains(&"input"), "{id} needs its raster input");
            assert!(names.contains(&"output"), "{id} needs its output path");
            assert!(names.contains(&extra), "{id} needs its {extra} control");
            assert!(
                manifest
                    .params
                    .iter()
                    .any(|p| p.name == "input" && p.required),
                "{id}'s input is the only required param",
            );
        }
    }

    const SAMPLE_DEM: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/sample.tif");

    /// A scratch directory holding the fixtures the run-based test feeds in: a
    /// point GeoJSON and a real D8 pointer raster derived from the sample DEM.
    /// The pointer matters -- several tools decode `d8_pntr` cell by cell and
    /// walk the flow network, and handing them an arbitrary DEM instead sends
    /// `basins` into a non-terminating traversal rather than an error.
    fn fixtures(registry: &ToolRegistry) -> (std::path::PathBuf, std::path::PathBuf) {
        let scratch = std::env::temp_dir().join("geolibre-cli-manifest-params");
        std::fs::create_dir_all(&scratch).expect("scratch dir");
        let vector = scratch.join("points.geojson");
        std::fs::write(
            &vector,
            r#"{"type":"FeatureCollection","features":[{"type":"Feature","properties":{"id":1},"geometry":{"type":"Point","coordinates":[0.5,0.5]}}]}"#,
        )
        .expect("scratch vector");

        let pointer = scratch.join("d8_pointer.tif");
        let mut args = ToolArgs::new();
        args.insert("dem".to_string(), Value::String(SAMPLE_DEM.to_string()));
        args.insert(
            "output".to_string(),
            Value::String(pointer.to_string_lossy().into_owned()),
        );
        let progress = StdoutProgress;
        let capabilities = AllowAllCapabilities;
        let ctx = ToolContext {
            progress: &progress,
            capabilities: &capabilities,
        };
        registry
            .run("d8_pointer", &args, &ctx)
            .expect("d8_pointer builds the pointer fixture");
        (vector, pointer)
    }

    /// Runs `manifest`'s tool with a value for every param it declares, and
    /// returns the error text (empty when the tool succeeded). Inputs point at
    /// the fixtures and outputs at a scratch directory, so the tools that can
    /// complete on a 64x48 DEM do, and the ones that cannot fail on the *data*
    /// rather than on a parameter they were never offered.
    fn run_with_manifest_params(
        registry: &ToolRegistry,
        manifest: &ToolManifest,
        vector: &std::path::Path,
        pointer: &std::path::Path,
    ) -> String {
        let scratch = std::env::temp_dir().join("geolibre-cli-manifest-params");
        let json = enriched_manifest(manifest, &backfilled_tool_ids(registry));
        let mut args = ToolArgs::new();
        for param in json["params"].as_array().into_iter().flatten() {
            let name = param["name"].as_str().unwrap_or_default().to_string();
            let data_kind = param["data_kind"].as_str().unwrap_or_default();
            let value = match param["io_role"].as_str() {
                Some("input") if data_kind == "vector" => {
                    Value::String(vector.to_string_lossy().into_owned())
                }
                Some("input") if name.contains("pntr") || name.contains("pointer") => {
                    Value::String(pointer.to_string_lossy().into_owned())
                }
                Some("input") => Value::String(SAMPLE_DEM.to_string()),
                Some("output") => Value::String(
                    scratch
                        .join(format!("{}_{name}.tif", manifest.id))
                        .to_string_lossy()
                        .into_owned(),
                ),
                _ => match param["schema"]["kind"].as_str() {
                    Some("enum") => param["schema"]["options"][0]["value"].clone(),
                    _ if data_kind == "bool" => Value::Bool(false),
                    _ if data_kind == "number" => Value::from(1),
                    // A required free-form string on these tools is a path (or a
                    // ';'-separated list of them), e.g. hypsometric_analysis'
                    // `inputs`; a user types one, so the test does too.
                    _ if param["required"].as_bool().unwrap_or(false) => {
                        Value::String(SAMPLE_DEM.to_string())
                    }
                    _ => continue,
                },
            };
            args.insert(name, value);
        }

        let progress = StdoutProgress;
        let capabilities = AllowAllCapabilities;
        let ctx = ToolContext {
            progress: &progress,
            capabilities: &capabilities,
        };
        match registry.run(&manifest.id, &args, &ctx) {
            Ok(_) => String::new(),
            Err(err) => err.to_string(),
        }
    }

    #[test]
    fn manifest_params_satisfy_the_runner() {
        // The manifests are only useful if a host can drive the tool from them,
        // so actually run the tools whose params this file has to reconstruct:
        // the two hillshades (fallback_params), the seven whose schema table is
        // stale (PREFER_METADATA_PARAMS), and a few plain schema-table backfills.
        // None may come back asking for a parameter it never advertised -- that
        // is the #327 failure, and it is what would silently return if any of
        // these lists drifted from upstream.
        let registry = build_registry();
        let manifests = manifests_with_metadata_params(&registry);
        let checked: Vec<&str> = ["hillshade", "multidirectional_hillshade"]
            .into_iter()
            .chain(PREFER_METADATA_PARAMS.iter().copied())
            .chain([
                "d8_pointer",
                "fill_depressions",
                "aspect",
                "basins",
                "extract_valleys",
                "sky_view_factor",
            ])
            .collect();

        let (vector, pointer) = fixtures(&registry);
        let mut unsatisfied: Vec<String> = Vec::new();
        for id in checked {
            let manifest = manifests
                .iter()
                .find(|m| m.id == id)
                .unwrap_or_else(|| panic!("{id} present"));
            let error = run_with_manifest_params(&registry, manifest, &vector, &pointer);
            if error.contains("missing required") {
                unsatisfied.push(format!("{id}: {error}"));
            }
        }
        assert!(
            unsatisfied.is_empty(),
            "tools asked for a parameter their manifest does not declare: {unsatisfied:#?}"
        );
    }

    #[test]
    fn a_backfilled_manifest_names_every_parameter_the_tool_takes() {
        // The names are the contract a host builds its args from, so they must be
        // the full set the tool reads -- ordered inputs, options, output.
        // Descriptions and `required` come from whitebox's generated param docs
        // (richer than the terse metadata prose), so they are not compared here.
        let registry = build_registry();
        let manifests = manifests_with_metadata_params(&registry);
        let by_id: BTreeMap<_, _> = manifests.iter().map(|m| (m.id.clone(), m)).collect();

        let expected = [
            ("d8_pointer", vec!["dem", "esri_pntr", "output"]),
            (
                "fill_depressions",
                vec![
                    "dem",
                    "fix_flats",
                    "flat_increment",
                    "flat_resolution",
                    "max_depth",
                    "output",
                ],
            ),
            ("aspect", vec!["input", "z_factor", "output"]),
            ("basins", vec!["d8_pntr", "esri_pntr", "output"]),
            // Generated by create_stream_tool_impl!, whose metadata gives every
            // tool the same d8_pntr/streams_raster/output placeholders; the real
            // params are the ones its run arm reads.
            (
                "extract_valleys",
                vec!["dem", "filter_size", "line_thin", "variant", "output"],
            ),
        ];
        for (id, names) in expected {
            let manifest = by_id.get(id).unwrap_or_else(|| panic!("{id} manifest"));
            assert_eq!(
                manifest
                    .params
                    .iter()
                    .map(|p| p.name.as_str())
                    .collect::<Vec<_>>(),
                names,
                "{id} params",
            );
        }
    }

    #[test]
    fn d8_pointer_manifest_exposes_its_parameters() {
        // The tool from the issue report: its manifest declared no params, so the
        // host rendered no form and the run failed validation on the first
        // required parameter. Named params, kinds and roles must all survive.
        let registry = build_registry();
        let manifest = manifests_with_metadata_params(&registry)
            .into_iter()
            .find(|m| m.id == "d8_pointer")
            .expect("d8_pointer present");
        assert_eq!(
            manifest
                .params
                .iter()
                .map(|p| p.name.as_str())
                .collect::<Vec<_>>(),
            ["dem", "esri_pntr", "output"],
        );
        assert!(
            manifest
                .params
                .iter()
                .any(|p| p.name == "dem" && p.required),
            "the DEM input stays required"
        );

        let json = enriched_manifest(&manifest, &backfilled_tool_ids(&registry));
        assert_eq!(param(&json, "dem")["schema"]["kind"], "input");
        assert_eq!(param(&json, "dem")["data_kind"], "raster");
        assert_eq!(param(&json, "output")["schema"]["kind"], "output");
    }

    // Looks up an enriched-manifest param by name; panics if it's absent.
    fn param<'a>(manifest: &'a Value, name: &str) -> &'a Value {
        manifest["params"]
            .as_array()
            .expect("params array")
            .iter()
            .find(|p| p["name"] == name)
            .unwrap_or_else(|| panic!("param '{name}' not found"))
    }

    #[test]
    fn spatial_join_manifest_renders_correct_controls() {
        // Regression for the broken demo inputs (opengeos/geolibre-rust#17). These
        // controls now come from wbcore's improved manifest inference rather than a
        // curated override: layer paths are file inputs, strategy/predicate are
        // enums (not a free-form box or a file picker), and distance is numeric.
        let registry = build_registry();
        let manifest = registry
            .manifests()
            .into_iter()
            .find(|m| m.id == "spatial_join")
            .expect("spatial_join present");
        let json = enriched_manifest(&manifest, &backfilled_tool_ids(&registry));

        for layer in ["target", "join"] {
            assert_eq!(param(&json, layer)["schema"]["kind"], "input");
            assert_eq!(param(&json, layer)["io_role"], "input");
        }
        assert_eq!(param(&json, "output")["schema"]["kind"], "output");
        assert_eq!(param(&json, "distance")["schema"]["kind"], "scalar");

        for enum_param in ["predicate", "strategy"] {
            let schema = &param(&json, enum_param)["schema"];
            assert_eq!(schema["kind"], "enum", "{enum_param} should be an enum");
            assert!(
                schema["options"].as_array().is_some_and(|o| !o.is_empty()),
                "{enum_param} should carry its option list"
            );
        }
    }
}
