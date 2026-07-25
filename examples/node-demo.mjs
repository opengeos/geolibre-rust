// Node smoke test for both package exports:
//   - the browser library (.)      : GeoTIFF/projection/vector/lidar/topology API
//   - the WASI tool runner (./tools): the whitebox + GeoLibre tool suite
//
//   node examples/node-demo.mjs path/to/dem.tif
import { readFile } from "node:fs/promises";
import initLib, {
  version,
  geotiff_info,
  vector_formats,
  vector_to_binary,
  vector_to_arrow_ipc,
} from "../npm/geolibre_wasm.js";
import { initTools, listTools, runTool } from "../npm/tools.mjs";

const demPath = process.argv[2] ?? new URL("./sample.tif", import.meta.url);
const dem = new Uint8Array(await readFile(demPath));

// ── library export (.) ──
await initLib(await readFile(new URL("../npm/geolibre_wasm_bg.wasm", import.meta.url)));
console.log(`library version: ${version()}`);
const info = JSON.parse(geotiff_info(dem));
if (!info.ok) throw new Error("geotiff_info failed");
console.log(`library geotiff_info: ${info.width}x${info.height}, epsg ${info.epsg}`);

// ── binary vector interop: typed arrays and Arrow IPC, no GeoJSON string ──
const sampleGeoJson = new TextEncoder().encode(
  JSON.stringify({
    type: "FeatureCollection",
    features: [
      { type: "Feature", geometry: { type: "Point", coordinates: [-0.1278, 51.5074] },
        properties: { name: "London", pop: 9000000 } },
      { type: "Feature", geometry: { type: "Point", coordinates: [10.7522, 59.9139] },
        properties: { name: "Oslo", pop: 700000 } },
    ],
  }),
);
console.log(`library vector_formats: ${vector_formats()}`);

const bin = vector_to_binary(sampleGeoJson, "geojson");
const positions = bin.point_positions();
const schema = JSON.parse(bin.schema_json);
const popIndex = schema.findIndex((f) => f.name === "pop");
const pops = bin.numeric_column(popIndex);
if (!(positions instanceof Float64Array)) throw new Error("expected Float64Array positions");
if (positions.length !== 4) throw new Error(`expected 4 position values, got ${positions.length}`);
if (pops[0] !== 9000000) throw new Error(`expected London's pop, got ${pops[0]}`);
console.log(
  `library vector_to_binary: ${bin.feature_count} features, ` +
    `${positions.length / bin.position_size} vertices, fields ${schema.map((f) => f.name).join()}`,
);
bin.free();

const ipc = vector_to_arrow_ipc(sampleGeoJson, "geojson");
// Arrow IPC streams start with the 0xFFFFFFFF continuation marker.
const marker = new DataView(ipc.buffer, ipc.byteOffset, 4).getUint32(0, true);
if (marker !== 0xffffffff) throw new Error(`not an Arrow IPC stream (marker ${marker})`);
console.log(`library vector_to_arrow_ipc: ${ipc.length} bytes of GeoArrow IPC`);

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
