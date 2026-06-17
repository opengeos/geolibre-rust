# geolibre-rust

The [`whitebox_next_gen`](https://github.com/jblindsay/whitebox_next_gen) pure-Rust
geospatial tool suite compiled to WebAssembly (WASI) for in-browser execution in
[GeoLibre](https://github.com/opengeos/GeoLibre).

It builds a single WASI binary (`geolibre-cli.wasm`) that exposes whitebox's
built-in tool registry **plus GeoLibre's own new tools** over a small
command-line contract, and ships a thin JS wrapper (`tools.mjs`) that runs it
through [`@bjorn3/browser_wasi_shim`](https://github.com/bjorn3/browser_wasi_shim)
with an in-memory filesystem. No server, no Python, no native install.

New tools live in the `geolibre-tools` crate and are registered alongside
whitebox's, so GeoLibre sees them through the same interface as the built-ins.

## Architecture

```
JS (browser/Node)                WASI binary (geolibre-cli.wasm)
-----------------                --------------------------------
tools.mjs                        crates/geolibre-cli (main.rs)
  write inputs -> /work    -->     argv -> ToolArgs (JSON)
  argv ["slope", "--..."]  -->     ToolRegistry::run
  read new files from /work <--      register_default_tools (whitebox)
                                     + geolibre_tools (new tools)
                                   tool writes via std::fs to /work
```

## Adding a new tool

1. Add a module with a `wbcore::Tool` impl under `crates/geolibre-tools/src/`
   (see `raster_normalize.rs` for the template: `metadata` / `validate` / `run`,
   reading and writing rasters by path).
2. Push it into the list returned by `geolibre_tools()` in
   `crates/geolibre-tools/src/lib.rs`.
3. Rebuild (`./build.sh`); it appears in `listTools()` automatically.

The crate depends only on `wbcore` plus the data crates a tool needs (e.g.
`wbraster`), so the same tools can later back a native CLI or the Python sidecar,
not just WASM.

The data boundary is the WASI virtual filesystem: inputs are placed under
`/work`, tools read/write there via ordinary `std::fs`, and any new file is
returned to JS as a `Uint8Array`. Raster outputs are Cloud Optimized GeoTIFFs.

## CLI contract

```text
geolibre list                 # print every tool id
geolibre manifests            # print all tool manifests as JSON (param schemas)
geolibre manifest <id>        # print one manifest as JSON
geolibre version              # print the runner version
geolibre <tool> [--k=v ...]   # run a tool over /work
```

`--key=value`, `--key value`, and bare `--flag` are all accepted. Values are
type-inferred: `true`/`false` -> bool, numbers -> number, everything else
(including `/work/...` paths) -> string.

## Build

```bash
rustup target add wasm32-wasip1
sudo apt-get install -y binaryen   # provides wasm-opt
./build.sh                         # -> npm/geolibre-cli.wasm
```

The `whitebox_next_gen` crates are referenced as path dependencies in
`crates/geolibre-cli/Cargo.toml`. Switch them to git or published versions before
releasing (note `wbtools_oss` is `publish = false`, so git or vendoring is
required).

### TODO: remove the vendored `kdtree` patch once `kdtree 0.8.1` ships

`vendor/kdtree/` and the `[patch.crates-io]` block in the root `Cargo.toml` work
around a bug in the published `kdtree 0.8.0` (it declares `criterion` as a normal
dependency, which pulls `rayon` and breaks the WASI build). The fix is already on
`kdtree-rs` `master` (PRs #70 and #89) but unreleased. Tracking issue:
https://github.com/mrhooray/kdtree-rs/issues/91

When `kdtree 0.8.1` (or later) is published, delete `vendor/kdtree/` and the
`[patch.crates-io]` block, then rebuild to confirm the WASI build stays green.

## Use from JavaScript

```js
import { runTool, listTools } from "geolibre-rust/tools";

const tools = await listTools();

const { files } = await runTool("slope", {
  args: ["--input=/work/dem.tif", "--output=/work/slope.tif", "--units=degrees"],
  input: { "dem.tif": demBytes }, // Uint8Array
});
const slopeCog = files["slope.tif"]; // Uint8Array (COG GeoTIFF)
```

## GeoLibre integration

The interface is byte-compatible with the existing `whitebox-wasm/tools` client:

1. Add `geolibre-rust` to `packages/processing/package.json`.
2. Add it to `optimizeDeps.exclude` in `apps/geolibre-desktop/vite.config.ts`
   (required for the `new URL("./*.wasm", import.meta.url)` glue).
3. Point `packages/processing/src/wasm-client.ts`'s lazy
   `import("whitebox-wasm/tools")` at `geolibre-rust/tools`, or add a sibling
   client and a source toggle in `ProcessingDialog.tsx`.

`listManifests()` is a value-add over the legacy package: it lets GeoLibre build
tool dialogs (parameter schemas, `raster_in`/`vector_in` roles) fully offline,
without the Python sidecar.

## License

MIT OR Apache-2.0
