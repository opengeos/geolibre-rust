// Node smoke test: load the WASI runner, list tools, and run one tool over a
// sample DEM placed in the in-memory /work dir.
//
//   node examples/node-demo.mjs path/to/dem.tif
import { readFile } from "node:fs/promises";
import { initTools, listTools, runTool } from "../npm/tools.mjs";

const wasmBytes = await readFile(new URL("../npm/geolibre-cli.wasm", import.meta.url));
await initTools(wasmBytes);

const tools = await listTools();
console.log(`tools available: ${tools.length}`);
if (tools.length === 0) throw new Error("expected a non-empty tool list");

const demPath = process.argv[2] ?? new URL("./sample.tif", import.meta.url);

const dem = new Uint8Array(await readFile(demPath));
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
