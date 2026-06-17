# geolibre-wasm

The [`whitebox_next_gen`](https://github.com/jblindsay/whitebox_next_gen) pure-Rust
geospatial tool suite, **plus new [GeoLibre](https://github.com/opengeos/GeoLibre)
tools**, compiled to WebAssembly (WASI) and runnable entirely in the browser,
Node, Deno, or any bundler. No server, no Python, no native install.

Source repo: [opengeos/geolibre-rust](https://github.com/opengeos/geolibre-rust).

## Install

```bash
npm install geolibre-wasm
```

It has a single runtime dependency, [`@bjorn3/browser_wasi_shim`](https://github.com/bjorn3/browser_wasi_shim),
which runs the WASI binary over an in-memory filesystem.

## Usage

```js
import { runTool, listTools, listManifests } from "geolibre-wasm/tools";

// every available tool id
const tools = await listTools();

// run a tool: inputs go into an in-memory /work dir; new files come back out
const { exitCode, stdout, files } = await runTool("slope", {
  args: ["--input=/work/dem.tif", "--output=/work/slope.tif", "--units=degrees"],
  input: { "dem.tif": demBytes }, // Uint8Array
});
const slopeCog = files["slope.tif"]; // Uint8Array (Cloud Optimized GeoTIFF)
```

Inputs are placed under `/work` (keyed by filename); any file a tool writes is
returned in `files`. Raster outputs are Cloud Optimized GeoTIFFs; vector outputs
are GeoJSON.

## API

| Export | Description |
|---|---|
| `initTools(source?)` | Compile the WASI runner once. Omit `source` in browsers/bundlers; pass wasm bytes or a URL/Response in Node. |
| `listTools(): Promise<string[]>` | Every available tool id. |
| `listManifests(): Promise<ToolManifest[]>` | All tool manifests (parameter schemas), for building UIs offline. |
| `runTool(tool, { args?, input? }): Promise<ToolResult>` | Run one tool over the in-memory filesystem. |

`ToolResult` is `{ exitCode: number, stdout: string[], files: Record<string, Uint8Array> }`.

## Notes

- Bundler users (Vite, etc.): exclude this package from dependency pre-bundling
  so the `new URL("./geolibre-cli.wasm", import.meta.url)` reference is preserved
  (e.g. Vite's `optimizeDeps.exclude`).
- Bounded by WebAssembly's ~4 GiB memory and single-threaded execution; use a
  server-side path for very large data.

## License

MIT OR Apache-2.0
