// geolibre-rust - run the whitebox_next_gen geospatial tool suite from
// JavaScript. The tools are the WASI binary `geolibre-cli.wasm`; this module
// executes them through a WASI shim with an in-memory filesystem, so they run in
// browsers, Node, Deno, and bundlers without a real disk. Raster outputs are
// Cloud Optimized GeoTIFFs.
//
//   import { runTool, listTools } from "geolibre-rust/tools";
//   const { files } = await runTool("slope", {
//     args: ["--input=/work/dem.tif", "--output=/work/slope.tif", "--units=degrees"],
//     input: { "dem.tif": demBytes },   // Uint8Array, placed under /work
//   });
//   const slopeCog = files["slope.tif"];  // Uint8Array
import { WASI, File, OpenFile, ConsoleStdout, PreopenDirectory } from "@bjorn3/browser_wasi_shim";

let _module = null;

/**
 * Compile the WASI tool runner once. In browsers/bundlers it loads the bundled
 * `geolibre-cli.wasm` relative to this module. In Node (no fetch of file URLs),
 * pass the wasm bytes or a URL/Response explicitly.
 * @param {URL|Response|BufferSource|string} [source]
 * @returns {Promise<WebAssembly.Module>}
 */
export async function initTools(source) {
  if (_module) return _module;
  if (!source) source = new URL("./geolibre-cli.wasm", import.meta.url);
  if (source instanceof Uint8Array || source instanceof ArrayBuffer) {
    _module = await WebAssembly.compile(source);
  } else if (source instanceof Response) {
    _module = await WebAssembly.compileStreaming(source);
  } else {
    _module = await WebAssembly.compileStreaming(fetch(source));
  }
  return _module;
}

async function exec(argv, inputFiles) {
  const mod = await initTools();
  const inNames = new Set(Object.keys(inputFiles));
  const contents = new Map(
    Object.entries(inputFiles).map(([k, v]) => [k, new File(new Uint8Array(v))]));
  const work = new PreopenDirectory("/work", contents);
  const stdout = [];
  const fds = [
    new OpenFile(new File(new Uint8Array())),
    ConsoleStdout.lineBuffered((s) => stdout.push(s)),
    ConsoleStdout.lineBuffered((s) => stdout.push(s)),
    work,
  ];
  const wasi = new WASI(["geolibre", ...argv], [], fds, { debug: false });
  const inst = await WebAssembly.instantiate(mod, { wasi_snapshot_preview1: wasi.wasiImport });
  let exitCode = 0;
  try { exitCode = wasi.start(inst); }
  catch (e) { if (e && e.constructor && e.constructor.name === "WASIProcExit") exitCode = e.code; else throw e; }
  const files = {};
  for (const [name, entry] of work.dir.contents) {
    if (entry.data && !inNames.has(name)) files[name] = entry.data;
  }
  return { exitCode, stdout, files };
}

/**
 * List every available tool id.
 * @returns {Promise<string[]>}
 */
export async function listTools() {
  const { stdout } = await exec(["list"], {});
  return stdout.map((s) => s.trim()).filter(Boolean);
}

/**
 * Fetch every tool manifest (id, display name, parameter schema, license tier).
 * Lets a host build tool dialogs fully offline, without a server.
 * @returns {Promise<object[]>}
 */
export async function listManifests() {
  const { stdout } = await exec(["manifests"], {});
  return JSON.parse(stdout.join(""));
}

/**
 * Run one tool over an in-memory filesystem.
 * @param {string} tool  tool id, e.g. "slope" (see {@link listTools})
 * @param {object} [opts]
 * @param {string[]} [opts.args]  CLI args, e.g. ["--input=/work/dem.tif","--output=/work/out.tif","--units=degrees"]
 * @param {Object<string, Uint8Array>} [opts.input]  files placed under /work (key = filename)
 * @returns {Promise<{exitCode:number, stdout:string[], files:Object<string,Uint8Array>}>}
 *   `files` contains any new files the tool wrote (e.g. the --output path).
 */
export async function runTool(tool, opts = {}) {
  const { args = [], input = {} } = opts;
  return exec([tool, ...args], input);
}
