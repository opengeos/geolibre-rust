// Node smoke test for both package exports:
//   - the browser library (.)      : GeoTIFF/projection/vector/lidar/topology API
//   - the WASI tool runner (./tools): the whitebox + GeoLibre tool suite
//
//   node examples/node-demo.mjs path/to/dem.tif
import { readFile } from "node:fs/promises";
import initLib, { version, geotiff_info } from "../npm/geolibre_wasm.js";
import { initTools, listTools, runTool } from "../npm/tools.mjs";

const demPath = process.argv[2] ?? new URL("./sample.tif", import.meta.url);
const dem = new Uint8Array(await readFile(demPath));

// ── library export (.) ──
await initLib(await readFile(new URL("../npm/geolibre_wasm_bg.wasm", import.meta.url)));
console.log(`library version: ${version()}`);
const info = JSON.parse(geotiff_info(dem));
if (!info.ok) throw new Error("geotiff_info failed");
console.log(`library geotiff_info: ${info.width}x${info.height}, epsg ${info.epsg}`);

// ── tools export (./tools) ──
await initTools(await readFile(new URL("../npm/geolibre-cli.wasm", import.meta.url)));

const tools = await listTools();
console.log(`tools available: ${tools.length}`);
if (tools.length === 0) throw new Error("expected a non-empty tool list");

const { exitCode, stdout, files } = await runTool("slope", {
  args: ["--input=/work/dem.tif", "--output=/work/slope.tif", "--units=degrees"],
  input: { "dem.tif": dem },
});

console.log("exitCode:", exitCode);
console.log("stdout:", stdout.join("\n"));
console.log("output files:", Object.keys(files));
if (exitCode !== 0) process.exit(1);
if (!files["slope.tif"]) throw new Error("expected slope.tif output");
console.log(`slope.tif: ${files["slope.tif"].length} bytes`);
