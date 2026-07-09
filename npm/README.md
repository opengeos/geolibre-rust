# geolibre-wasm

[![npm version](https://img.shields.io/npm/v/geolibre-wasm.svg)](https://www.npmjs.com/package/geolibre-wasm)
[![npm downloads](https://img.shields.io/npm/dm/geolibre-wasm.svg)](https://www.npmjs.com/package/geolibre-wasm)
[![license](https://img.shields.io/npm/l/geolibre-wasm.svg)](https://github.com/opengeos/geolibre-rust#license)

A pure-Rust geospatial toolkit compiled to WebAssembly, runnable entirely in the
browser, Node, Deno, or any bundler. No server, no Python, no native install.

It is a **superset of [`whitebox-wasm`](https://www.npmjs.com/package/whitebox-wasm)**:
everything that package offers, plus new
[GeoLibre](https://github.com/opengeos/GeoLibre) tools. Built on
[`opengeos/whitebox-wasm`](https://github.com/opengeos/whitebox-wasm), the
WASM-ready fork of [`whitebox_next_gen`](https://github.com/jblindsay/whitebox_next_gen).

Source repo: [opengeos/geolibre-rust](https://github.com/opengeos/geolibre-rust).

Two layers, two entry points:

- **`geolibre-wasm`** (the `.` export) -- a `wasm-bindgen` browser library with
  typed in-memory APIs for GeoTIFF/COG read+write, projections, vector, and LiDAR.
- **`geolibre-wasm/tools`** (the `./tools` export) -- a WASI tool runner exposing
  the whitebox tool suite plus GeoLibre's own tools.

## Install

```bash
npm install geolibre-wasm
```

`geolibre-wasm/tools` uses [`@bjorn3/browser_wasi_shim`](https://github.com/bjorn3/browser_wasi_shim)
(a runtime dependency) to run the WASI binary over an in-memory filesystem.

## Library usage (`.`)

```js
import init, { GeoTiffReader, CogBuilder, geotiff_info, version } from "geolibre-wasm";

await init(); // load the wasm-bindgen module (browsers/bundlers)

const r = new GeoTiffReader(tiffBytes); // Uint8Array
console.log(r.width, r.height, r.bands, r.epsg);
const band0 = r.read_band_f64(0);        // Float64Array

// header-only metadata (works on multi-GB files)
const meta = JSON.parse(geotiff_info(tiffBytes));
```

Classes include `GeoTiffReader` (parse once, read many), `CogBuilder` (encode
Cloud Optimized GeoTIFFs), and `CogStream` (range-request tiled COG reading).

## Tools usage (`./tools`)

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

COG subsets can be extracted from a local file or an HTTP source. HTTP sources
use byte-range requests and do not download the full object:

```js
const { files } = await runTool("extract_cog_subset", {
  args: [
    "--url=https://example.com/scene.cog.tif",
    "--bbox=-122.55,37.70,-122.35,37.84",
    "--bbox_crs=4326",
    "--output_crs=3857",
    "--resolution=10",
    "--nodata=-9999",
    "--output=/work/subset.tif",
  ],
});
const subsetCog = files["subset.tif"];
```

WMS subsets can be requested from a GetMap endpoint that supports GeoTIFF
responses:

```js
const { files } = await runTool("extract_wms_subset", {
  args: [
    "--url=https://example.com/geoserver/wms",
    "--layers=workspace:layer",
    "--bbox=-122.55,37.70,-122.35,37.84",
    "--bbox_crs=4326",
    "--resolution=0.0001",
    "--format=image/geotiff",
    "--output=/work/wms_subset.tif",
  ],
});
```

For a local COG, upload/provide bytes through `input` and point `--input` at
the `/work` path:

```js
const { files } = await runTool("extract_cog_subset", {
  args: [
    "--input=/work/local.cog.tif",
    "--bbox=-122.55,37.70,-122.35,37.84",
    "--bbox_crs=4326",
    "--output=/work/subset.tif",
  ],
  input: { "local.cog.tif": localCogBytes },
});
```

Inputs are placed under `/work` (keyed by filename); any file a tool writes is
returned in `files`. Raster outputs are Cloud Optimized GeoTIFFs; vector outputs
are GeoJSON.

When `output_crs` is set, `resolution` is interpreted in output CRS units and
the subset is reprojected with nearest-neighbor resampling. Without
`output_crs`, the output usually stays on the source COG grid. Sources with
user-defined projection strings that cannot be written as EPSG metadata, such as
the NLCD sample, default to `bbox_crs` output so GIS viewers can place the
result. Outputs are Deflate-compressed by default. `nodata` overrides the output
nodata metadata and is used to fill pixels introduced during reprojection.
The subset writer preserves the source sample type for supported COG types
(`uint8`, `float32`, and `float64`); palette-indexed `uint8` rasters retain
their source ColorMap.

## API

| Export | Description |
|---|---|
| `initTools(source?)` | Compile the WASI runner once. Omit `source` in browsers/bundlers; pass wasm bytes or a URL/Response in Node. |
| `listTools(): Promise<string[]>` | Every available tool id. |
| `listManifests(): Promise<ToolManifest[]>` | All tool manifests (parameter schemas), for building UIs offline. |
| `runTool(tool, { args?, input? }): Promise<ToolResult>` | Run one tool over the in-memory filesystem. |
| `extractCogSubset(source, opts): Promise<Uint8Array>` | Extract a local or HTTP COG subset directly. |
| `extractWmsSubset(url, opts): Promise<Uint8Array>` | Request a WMS GeoTIFF subset and encode it as a COG. |

`ToolResult` is `{ exitCode: number, stdout: string[], files: Record<string, Uint8Array> }`.

## Notes

- Bundler users (Vite, etc.): exclude this package from dependency pre-bundling
  so the `new URL("./geolibre-cli.wasm", import.meta.url)` reference is preserved
  (e.g. Vite's `optimizeDeps.exclude`).
- Bounded by WebAssembly's ~4 GiB memory and single-threaded execution; use a
  server-side path for very large data.

## License

MIT
