// geolibre-wasm - run the whitebox_next_gen geospatial tool suite from
// JavaScript. The tools are the WASI binary `geolibre-cli.wasm`; this module
// executes them through a WASI shim with an in-memory filesystem, so they run in
// browsers, Node, Deno, and bundlers without a real disk. Raster outputs are
// Cloud Optimized GeoTIFFs.
//
//   import { runTool, listTools } from "geolibre-wasm/tools";
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

// Resolve one input value to bytes: a Uint8Array/ArrayBuffer as-is, or an
// http(s) URL string to fetch. Format-agnostic (raster or vector).
async function materializeInput(value) {
  if (typeof value === "string") {
    if (!/^https?:\/\//i.test(value))
      throw new Error(`input string must be an http(s) URL, got: ${value}`);
    // A User-Agent helps with CDNs that reject non-browser agents (browsers
    // ignore this header and send their own; Node/undici honors it).
    const resp = await fetch(value, { headers: { "User-Agent": "Mozilla/5.0 (geolibre-wasm)" } });
    return new Uint8Array(await resp.arrayBuffer());
  }
  return new Uint8Array(value);
}

async function exec(argv, inputFiles) {
  const mod = await initTools();
  const inNames = new Set(Object.keys(inputFiles));
  const entries = await Promise.all(
    Object.entries(inputFiles).map(async ([k, v]) => [k, new File(await materializeInput(v))]));
  const contents = new Map(entries);
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
  // Collect every new file under /work, recursing into subdirectories so tools
  // that write a tree (e.g. raster_to_tiles' {z}/{x}/{y}.png) are surfaced too.
  // Keys are paths relative to /work (nested files use "/" separators).
  const files = {};
  const walk = (dir, prefix) => {
    for (const [name, entry] of dir.contents) {
      const rel = prefix ? `${prefix}/${name}` : name;
      if (entry && entry.contents) {
        walk(entry, rel); // subdirectory
      } else if (entry && entry.data && !(prefix === "" && inNames.has(name))) {
        files[rel] = entry.data; // new file (top-level inputs excluded)
      }
    }
  };
  walk(work.dir, "");
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
