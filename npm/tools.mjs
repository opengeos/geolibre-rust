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
import initGeoLibre, { CogBuilder, CogStream, GeoTiffReader, PmtilesExtractor, transform_bbox_epsg, transform_points_epsg } from "./geolibre_wasm.js";

let _module = null;
let _libraryReady = null;
const COG_SUBSET_TOOL_ID = "extract_cog_subset";
const WMS_SUBSET_TOOL_ID = "extract_wms_subset";
const XYZ_SUBSET_TOOL_ID = "extract_xyz_tile_subset";
const PMTILES_EXTRACT_TOOL_ID = "pmtiles_extract";
// The `url` + `bbox_crs` params the JS interception adds to the WASI
// `pmtiles_extract` tool's own manifest, so a host can extract from a remote
// archive by byte-range (no whole-file download) and offer the same
// "Use map extent" / "Draw on map" bbox affordance the COG/WMS/XYZ extractors
// get. The native `geolibre` CLI keeps its local file→file behaviour.
const PMTILES_EXTRACT_EXTRA_PARAMS = [
  { name: "url", description: "HTTP(S) PMTiles URL. Provide either url or input; HTTP sources use byte-range requests (the host must allow cross-origin Range reads).", required: false, schema: { type: "string" } },
  { name: "bbox_crs", description: "EPSG code of the bbox coordinates. Defaults to 4326 (WGS84 lon/lat); other CRSs are reprojected to 4326.", required: false, schema: { kind: "scalar", scalar: "integer" } },
];
const COG_SUBSET_MANIFEST = {
  id: COG_SUBSET_TOOL_ID,
  display_name: "Extract COG Subset",
  summary: "Extract a bbox subset from a local or HTTP Cloud Optimized GeoTIFF.",
  category: "Raster",
  license_tier: "Open",
  source: "geolibre",
  params: [
    { name: "input", description: "Local COG file. Provide either input or url.", required: false, schema: { kind: "input", data_kind: "raster", sample_default: false } },
    { name: "url", description: "HTTP(S) COG URL. Provide either url or input; HTTP sources use byte-range requests.", required: false, schema: { type: "string" } },
    { name: "bbox", description: "Bounding box as minX,minY,maxX,maxY in bbox_crs.", required: true, schema: { type: "string" } },
    { name: "bbox_crs", description: "EPSG code of bbox coordinates.", required: true, schema: { kind: "scalar", scalar: "integer" } },
    { name: "output", description: "Output COG path.", required: false, schema: { type: "string" } },
    { name: "level", description: "COG overview level to read; 0 is full resolution.", required: false, schema: { kind: "scalar", scalar: "integer" } },
    { name: "resolution", description: "Target output pixel size. Uses output_crs units when output_crs is set; otherwise bbox_crs units. Selects the closest COG overview when level is omitted.", required: false, schema: { kind: "scalar", scalar: "float" } },
    { name: "output_crs", description: "Optional output EPSG code. When set, the subset is reprojected to this CRS. For source CRSs that cannot be written as EPSG metadata, defaults to bbox_crs.", required: false, schema: { kind: "scalar", scalar: "integer" } },
    { name: "nodata", description: "Optional output nodata value. Used as reprojection fill and written as output nodata metadata.", required: false, schema: { kind: "scalar", scalar: "float" } },
  ],
};
const WMS_SUBSET_MANIFEST = {
  id: WMS_SUBSET_TOOL_ID,
  display_name: "Extract WMS Subset",
  summary: "Request a bbox subset from a WMS GetMap endpoint and write it as a Deflate COG.",
  category: "Raster",
  license_tier: "Open",
  source: "geolibre",
  params: [
    { name: "url", description: "WMS endpoint URL.", required: true, schema: { type: "string" } },
    { name: "layers", description: "WMS layer name(s), comma-separated.", required: true, schema: { type: "string" } },
    { name: "styles", description: "WMS style name(s), comma-separated. Defaults to empty.", required: false, schema: { type: "string" } },
    { name: "bbox", description: "Bounding box as minX,minY,maxX,maxY in bbox_crs.", required: true, schema: { type: "string" } },
    { name: "bbox_crs", description: "EPSG code of bbox coordinates.", required: true, schema: { kind: "scalar", scalar: "integer" } },
    { name: "output", description: "Output COG path.", required: false, schema: { type: "string" } },
    { name: "resolution", description: "Target output pixel size in output_crs units when output_crs is set; otherwise bbox_crs units.", required: false, schema: { kind: "scalar", scalar: "float" } },
    { name: "width", description: "Output request width in pixels. Used when resolution is omitted.", required: false, schema: { kind: "scalar", scalar: "integer" } },
    { name: "height", description: "Output request height in pixels. Used when resolution is omitted.", required: false, schema: { kind: "scalar", scalar: "integer" } },
    { name: "output_crs", description: "Optional output/request EPSG code. Defaults to bbox_crs.", required: false, schema: { kind: "scalar", scalar: "integer" } },
    { name: "format", description: "WMS image format. Defaults to image/geotiff.", required: false, schema: { type: "string" } },
    { name: "version", description: "WMS version. Defaults to 1.1.1.", required: false, schema: { type: "string" } },
    { name: "nodata", description: "Optional output nodata value.", required: false, schema: { kind: "scalar", scalar: "float" } },
  ],
};
const XYZ_SUBSET_MANIFEST = {
  id: XYZ_SUBSET_TOOL_ID,
  display_name: "Extract XYZ Tile Subset",
  summary: "Fetch XYZ raster tiles for a bbox, mosaic them, and write a Deflate RGB COG.",
  category: "Raster",
  license_tier: "Open",
  source: "geolibre",
  params: [
    { name: "url", description: "XYZ tile URL template, e.g. https://.../{z}/{x}/{y}.png.", required: true, schema: { type: "string" } },
    { name: "zoom", description: "XYZ zoom level.", required: true, schema: { kind: "scalar", scalar: "integer" } },
    { name: "bbox", description: "Bounding box as minX,minY,maxX,maxY in bbox_crs.", required: true, schema: { type: "string" } },
    { name: "bbox_crs", description: "EPSG code of bbox coordinates.", required: true, schema: { kind: "scalar", scalar: "integer" } },
    { name: "output", description: "Output COG path.", required: false, schema: { type: "string" } },
    { name: "resolution", description: "Target output pixel size in output_crs units. Defaults to native tile resolution.", required: false, schema: { kind: "scalar", scalar: "float" } },
    { name: "width", description: "Output width in pixels. Used when resolution is omitted.", required: false, schema: { kind: "scalar", scalar: "integer" } },
    { name: "height", description: "Output height in pixels. Used when resolution is omitted.", required: false, schema: { kind: "scalar", scalar: "integer" } },
    { name: "output_crs", description: "Optional output EPSG code. Defaults to bbox_crs.", required: false, schema: { kind: "scalar", scalar: "integer" } },
    { name: "tile_size", description: "Tile size in pixels. Defaults to 256.", required: false, schema: { kind: "scalar", scalar: "integer" } },
    { name: "subdomains", description: "Optional subdomain letters for {s}, e.g. abc.", required: false, schema: { type: "string" } },
    { name: "nodata", description: "Optional output nodata value used for missing/transparent pixels.", required: false, schema: { kind: "scalar", scalar: "float" } },
  ],
};

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

async function initLibrary() {
  if (!_libraryReady) _libraryReady = initGeoLibre();
  return _libraryReady;
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
  const tools = stdout.map((s) => s.trim()).filter(Boolean);
  if (!tools.includes(COG_SUBSET_TOOL_ID)) tools.push(COG_SUBSET_TOOL_ID);
  if (!tools.includes(WMS_SUBSET_TOOL_ID)) tools.push(WMS_SUBSET_TOOL_ID);
  if (!tools.includes(XYZ_SUBSET_TOOL_ID)) tools.push(XYZ_SUBSET_TOOL_ID);
  return tools;
}

/**
 * Fetch every tool manifest (id, display name, parameter schema, license tier).
 * Lets a host build tool dialogs fully offline, without a server.
 * @returns {Promise<object[]>}
 */
export async function listManifests() {
  const { stdout } = await exec(["manifests"], {});
  const manifests = JSON.parse(stdout.join(""));
  if (!manifests.some((m) => m.id === COG_SUBSET_TOOL_ID)) {
    manifests.push(COG_SUBSET_MANIFEST);
  }
  if (!manifests.some((m) => m.id === WMS_SUBSET_TOOL_ID)) {
    manifests.push(WMS_SUBSET_MANIFEST);
  }
  if (!manifests.some((m) => m.id === XYZ_SUBSET_TOOL_ID)) {
    manifests.push(XYZ_SUBSET_MANIFEST);
  }
  // pmtiles_extract is a real WASI tool, so its manifest already lists
  // input/output/bbox/min_zoom/max_zoom/max_tiles. Augment it with the `url`
  // and `bbox_crs` params the JS interception handles, so hosts can extract
  // from a remote archive and reproject the bbox.
  const pmtiles = manifests.find((m) => m.id === PMTILES_EXTRACT_TOOL_ID);
  if (pmtiles && Array.isArray(pmtiles.params)) {
    pmtiles.source = "geolibre";
    for (const extra of PMTILES_EXTRACT_EXTRA_PARAMS) {
      if (!pmtiles.params.some((p) => p.name === extra.name)) {
        pmtiles.params.push({ ...extra });
      }
    }
  }
  return manifests;
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
  if (tool === COG_SUBSET_TOOL_ID) return runCogSubsetTool(args, input);
  if (tool === WMS_SUBSET_TOOL_ID) return runWmsSubsetTool(args);
  if (tool === XYZ_SUBSET_TOOL_ID) return runXyzTileSubsetTool(args);
  if (tool === PMTILES_EXTRACT_TOOL_ID) return runPmtilesExtractTool(args, input);
  return exec([tool, ...args], input);
}

function parseFlagArgs(args) {
  const out = {};
  for (let i = 0; i < args.length; i++) {
    const token = args[i];
    if (!token.startsWith("--")) continue;
    const stripped = token.slice(2);
    if (stripped.includes("=")) {
      const [key, ...rest] = stripped.split("=");
      out[key] = rest.join("=");
    } else if (i + 1 < args.length && !args[i + 1].startsWith("--")) {
      out[stripped] = args[++i];
    } else {
      out[stripped] = true;
    }
  }
  return out;
}

function parseBbox(raw) {
  const bbox = String(raw || "").split(",").map((v) => Number(v.trim()));
  if (bbox.length !== 4 || bbox.some((v) => !Number.isFinite(v)) || bbox[0] >= bbox[2] || bbox[1] >= bbox[3]) {
    throw new Error("--bbox must be ordered as minX,minY,maxX,maxY");
  }
  return bbox;
}

function parseOptionalNumber(raw, name) {
  if (raw == null || raw === true || String(raw).trim() === "") return undefined;
  const value = Number(raw);
  if (!Number.isFinite(value)) throw new Error(`--${name} must be a number`);
  return value;
}

function parseOptionalInteger(raw, name) {
  if (raw == null || raw === true || String(raw).trim() === "") return undefined;
  const value = Number(raw);
  if (!Number.isInteger(value) || value <= 0) throw new Error(`--${name} must be a positive integer`);
  return value;
}

function parseOptionalEpsg(raw, name) {
  if (raw == null || raw === true || String(raw).trim() === "") return undefined;
  const value = Number(raw);
  if (!Number.isInteger(value) || value <= 0) throw new Error(`--${name} must be a positive EPSG code`);
  return value;
}

function outputKey(path) {
  if (!path || path === true) return "subset.tif";
  const normalized = String(path).replace(/^\/work\/?/, "");
  return normalized || "subset.tif";
}

async function runCogSubsetTool(args, inputFiles) {
  const flags = parseFlagArgs(args);
  const url = flags.url;
  const inputPath = flags.input;
  const bbox = parseBbox(flags.bbox);
  const bboxCrs = Number(flags.bbox_crs ?? flags.bboxCrs ?? flags.crs);
  const level = parseOptionalNumber(flags.level, "level");
  const resolution = parseOptionalNumber(flags.resolution, "resolution");
  const outputCrs = parseOptionalEpsg(flags.output_crs ?? flags.outputCrs, "output_crs");
  const nodata = parseOptionalNumber(flags.nodata, "nodata");
  const key = outputKey(flags.output);
  const stdout = [];

  try {
    const source = await resolveCogSubsetSource({ url, inputPath, inputFiles });
    const bytes = await extractCogSubset(source, { bbox, bboxCrs, level, resolution, outputCrs, nodata });
    stdout.push(JSON.stringify({ output: `/work/${key}`, bytes: bytes.byteLength }));
    return { exitCode: 0, stdout, files: { [key]: bytes } };
  } catch (error) {
    stdout.push(String(error?.message || error));
    return { exitCode: 1, stdout, files: {} };
  }
}

async function runWmsSubsetTool(args) {
  const flags = parseFlagArgs(args);
  const url = flags.url;
  const layers = flags.layers;
  const styles = flags.styles;
  const bbox = parseBbox(flags.bbox);
  const bboxCrs = Number(flags.bbox_crs ?? flags.bboxCrs ?? flags.crs);
  const resolution = parseOptionalNumber(flags.resolution, "resolution");
  const width = parseOptionalInteger(flags.width, "width");
  const height = parseOptionalInteger(flags.height, "height");
  const outputCrs = parseOptionalEpsg(flags.output_crs ?? flags.outputCrs, "output_crs");
  const nodata = parseOptionalNumber(flags.nodata, "nodata");
  const format = flags.format == null || flags.format === true ? undefined : String(flags.format);
  const version = flags.version == null || flags.version === true ? undefined : String(flags.version);
  const key = outputKey(flags.output);
  const stdout = [];

  try {
    const bytes = await extractWmsSubset(url, {
      layers, styles, bbox, bboxCrs, resolution, width, height, outputCrs, nodata, format, version,
    });
    stdout.push(JSON.stringify({ output: `/work/${key}`, bytes: bytes.byteLength }));
    return { exitCode: 0, stdout, files: { [key]: bytes } };
  } catch (error) {
    stdout.push(String(error?.message || error));
    return { exitCode: 1, stdout, files: {} };
  }
}

async function runXyzTileSubsetTool(args) {
  const flags = parseFlagArgs(args);
  const url = flags.url;
  const zoom = parseOptionalInteger(flags.zoom, "zoom");
  const bbox = parseBbox(flags.bbox);
  const bboxCrs = Number(flags.bbox_crs ?? flags.bboxCrs ?? flags.crs);
  const resolution = parseOptionalNumber(flags.resolution, "resolution");
  const width = parseOptionalInteger(flags.width, "width");
  const height = parseOptionalInteger(flags.height, "height");
  const outputCrs = parseOptionalEpsg(flags.output_crs ?? flags.outputCrs, "output_crs");
  const tileSize = parseOptionalInteger(flags.tile_size ?? flags.tileSize, "tile_size");
  const nodata = parseOptionalNumber(flags.nodata, "nodata");
  const subdomains = flags.subdomains == null || flags.subdomains === true ? undefined : String(flags.subdomains);
  const key = outputKey(flags.output);
  const stdout = [];

  try {
    const bytes = await extractXyzTileSubset(url, {
      zoom, bbox, bboxCrs, resolution, width, height, outputCrs, tileSize, subdomains, nodata,
    });
    stdout.push(JSON.stringify({ output: `/work/${key}`, bytes: bytes.byteLength }));
    return { exitCode: 0, stdout, files: { [key]: bytes } };
  } catch (error) {
    stdout.push(String(error?.message || error));
    return { exitCode: 1, stdout, files: {} };
  }
}

/**
 * Extract a bbox/zoom subset of a PMTiles archive (local or remote) into a new
 * self-contained archive, driving the wasm `PmtilesExtractor` with byte-range
 * reads. Mirrors the COG/WMS/XYZ interception so a remote planet build is
 * subset by range instead of downloaded whole. The bbox is reprojected from
 * `bbox_crs` to WGS84 (what PMTiles address); `bbox_crs` defaults to 4326.
 */
async function runPmtilesExtractTool(args, inputFiles) {
  const flags = parseFlagArgs(args);
  const url = flags.url;
  const inputPath = flags.input;
  const rawBbox = parseBbox(flags.bbox);
  const bboxCrs = parseOptionalEpsg(flags.bbox_crs ?? flags.bboxCrs ?? flags.crs, "bbox_crs");
  const minZoom = clampZoom(flags.min_zoom ?? flags.minZoom, 0);
  const maxZoom = clampZoom(flags.max_zoom ?? flags.maxZoom, 30);
  const maxTiles = parseOptionalInteger(flags.max_tiles ?? flags.maxTiles, "max_tiles");
  const key = pmtilesOutputKey(flags.output);
  const stdout = [];

  try {
    await initLibrary();
    // PMTiles are addressed in WGS84 lon/lat; reproject if the caller gave a
    // bbox in another CRS (the map extent affordance always sends 4326).
    const bbox =
      bboxCrs && bboxCrs !== 4326
        ? Array.from(transform_bbox_epsg(bboxCrs, 4326, Float64Array.from(rawBbox)))
        : rawBbox;

    const reader = makeSourceReader(
      await resolvePmtilesSource({ url, inputPath, inputFiles }),
    );
    const archive = await drivePmtilesExtractor(reader, {
      bbox,
      minZoom,
      maxZoom,
      maxTiles,
    });
    stdout.push(
      JSON.stringify({ output: `/work/${key}`, bytes: archive.byteLength }),
    );
    return { exitCode: 0, stdout, files: { [key]: archive } };
  } catch (error) {
    stdout.push(String(error?.message || error));
    return { exitCode: 1, stdout, files: {} };
  }
}

/**
 * Drive a `PmtilesExtractor` to completion over a range reader, fetching each
 * round of wanted ranges concurrently.
 * @param {{range: (offset:number, length:number) => Promise<Uint8Array>}} reader
 * @param {{bbox:number[], minZoom:number, maxZoom:number, maxTiles?:number}} opts
 * @returns {Promise<Uint8Array>}
 */
async function drivePmtilesExtractor(reader, { bbox, minZoom, maxZoom, maxTiles }) {
  const extractor = new PmtilesExtractor(
    bbox[0], bbox[1], bbox[2], bbox[3], minZoom, maxZoom,
  );
  try {
    if (maxTiles != null) extractor.set_max_tiles(maxTiles);
    while (!extractor.done) {
      const wanted = JSON.parse(extractor.wanted_json());
      if (wanted.length === 0) {
        throw new Error("PMTiles extractor stalled: not done but nothing wanted");
      }
      const chunks = await Promise.all(
        wanted.map(({ offset, length }) => reader.range(offset, length)),
      );
      wanted.forEach(({ offset }, i) => extractor.feed(offset, chunks[i]));
    }
    return extractor.finish();
  } finally {
    extractor.free();
  }
}

/** Resolve a PMTiles source to a URL string or local bytes for makeSourceReader. */
async function resolvePmtilesSource({ url, inputPath, inputFiles }) {
  const hasUrl = url && url !== true && String(url).trim() !== "";
  if (!hasUrl && !inputPath) {
    throw new Error("provide either --url=<http pmtiles> or --input=/work/local.pmtiles");
  }
  if (hasUrl && inputPath) {
    throw new Error("provide only one source: --url or --input");
  }
  if (hasUrl) return String(url).trim();

  const key = String(inputPath).replace(/^\/work\/?/, "");
  if (!inputFiles || !(key in inputFiles)) {
    throw new Error(`input file not found in /work: ${inputPath}`);
  }
  return materializeInput(inputFiles[key]);
}

/** Clamp/parse a zoom flag to an integer in [0, 30], falling back to a default. */
function clampZoom(raw, fallback) {
  if (raw == null || raw === true || String(raw).trim() === "") return fallback;
  const value = Number(raw);
  if (!Number.isInteger(value) || value < 0 || value > 30) {
    throw new Error("zoom levels must be whole numbers in 0..=30");
  }
  return value;
}

/** Default a PMTiles output path to extract.pmtiles, stripping the /work prefix. */
function pmtilesOutputKey(path) {
  if (!path || path === true) return "extract.pmtiles";
  const normalized = String(path).replace(/^\/work\/?/, "");
  return normalized || "extract.pmtiles";
}

async function resolveCogSubsetSource({ url, inputPath, inputFiles }) {
  if ((url == null || url === true || String(url).trim() === "") && !inputPath) {
    throw new Error("provide either --url=<http COG> or --input=/work/local.tif");
  }
  if (url && url !== true && inputPath) {
    throw new Error("provide only one source: --url or --input");
  }
  if (url && url !== true) return String(url).trim();

  const key = outputKey(inputPath);
  if (!inputFiles || !(key in inputFiles)) {
    throw new Error(`input file not found in /work: ${inputPath}`);
  }
  return materializeInput(inputFiles[key]);
}

async function fetchRange(url, offset, length, fetchOptions) {
  const end = offset + length - 1;
  const headers = new Headers(fetchOptions?.headers || {});
  headers.set("Range", `bytes=${offset}-${end}`);
  try {
    if (!headers.has("User-Agent")) headers.set("User-Agent", "Mozilla/5.0 (geolibre-wasm)");
  } catch {
    // Browsers treat User-Agent as a forbidden header; Node accepts it.
  }

  const resp = await fetch(url, { ...fetchOptions, headers });
  if (resp.status !== 206) {
    throw new Error(`server must support HTTP range requests (expected 206, got ${resp.status})`);
  }
  return new Uint8Array(await resp.arrayBuffer());
}

function makeSourceReader(source, fetchOptions) {
  if (typeof source === "string") {
    if (!/^https?:\/\//i.test(source)) throw new Error(`url must be HTTP(S), got: ${source}`);
    return {
      type: "http",
      async range(offset, length) {
        return fetchRange(source, offset, length, fetchOptions);
      },
    };
  }

  const bytes = new Uint8Array(source);
  return {
    type: "local",
    async range(offset, length) {
      if (offset < 0 || length < 0 || offset >= bytes.byteLength) {
        throw new Error(`requested byte range ${offset}-${offset + length - 1} exceeds local COG size`);
      }
      return bytes.slice(offset, Math.min(bytes.byteLength, offset + length));
    },
  };
}

async function openCogStream(reader, options) {
  const maxHeaderBytes = options.maxHeaderBytes ?? 8 * 1024 * 1024;
  let headerBytes = options.initialHeaderBytes ?? 256 * 1024;
  let lastError = null;

  while (headerBytes <= maxHeaderBytes) {
    const prefix = await reader.range(0, headerBytes);
    try {
      return { stream: new CogStream(prefix), headerBytes, header: prefix };
    } catch (error) {
      lastError = error;
      const message = String(error?.message || error);
      if (!/(need more header bytes|failed to fill whole buffer)/i.test(message)) throw error;
      headerBytes *= 2;
    }
  }

  throw new Error(`could not parse COG header within ${maxHeaderBytes} bytes: ${lastError}`);
}

function parseLevels(stream) {
  return JSON.parse(stream.levels_json());
}

function tiffAccess(bytes) {
  if (bytes[0] !== 0x49 || bytes[1] !== 0x49) throw new Error("only little-endian TIFF metadata is supported");
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const magic = dv.getUint16(2, true);
  const big = magic === 43;
  if (!big && magic !== 42) throw new Error("not a TIFF header");
  return {
    dv,
    big,
    inlineBytes: big ? 8 : 4,
    firstIfdOffset: Number(big ? dv.getBigUint64(8, true) : dv.getUint32(4, true)),
    readOffset(pos) { return Number(big ? dv.getBigUint64(pos, true) : dv.getUint32(pos, true)); },
    readCount(pos) { return Number(big ? dv.getBigUint64(pos, true) : dv.getUint16(pos, true)); },
    writeFirstIfd(out, offset) {
      if (big) new DataView(out.buffer).setBigUint64(8, BigInt(offset), true);
      else new DataView(out.buffer).setUint32(4, offset, true);
    },
  };
}

const TIFF_TYPE_BYTES = { 1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 12: 8, 16: 8 };

function readTiffIfd(bytes, offset, big) {
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const count = Number(big ? dv.getBigUint64(offset, true) : dv.getUint16(offset, true));
  const countBytes = big ? 8 : 2;
  const entryBytes = big ? 20 : 12;
  const inlineBytes = big ? 8 : 4;
  const entries = [];
  for (let i = 0; i < count; i++) {
    const pos = offset + countBytes + i * entryBytes;
    const tag = dv.getUint16(pos, true);
    const type = dv.getUint16(pos + 2, true);
    const valueCount = Number(big ? dv.getBigUint64(pos + 4, true) : dv.getUint32(pos + 4, true));
    const bytesLen = valueCount * (TIFF_TYPE_BYTES[type] || 1);
    const valuePos = pos + (big ? 12 : 8);
    const valueOffset = bytesLen <= inlineBytes
      ? valuePos
      : Number(big ? dv.getBigUint64(valuePos, true) : dv.getUint32(valuePos, true));
    entries.push({ tag, type, count: valueCount, bytesLen, valuePos, valueOffset });
  }
  const nextOffsetPos = offset + countBytes + count * entryBytes;
  const nextOffset = Number(big ? dv.getBigUint64(nextOffsetPos, true) : dv.getUint32(nextOffsetPos, true));
  return { count, entries, nextOffset };
}

function readShortTag(bytes, ifd, tag) {
  const entry = ifd.entries.find((e) => e.tag === tag);
  if (!entry || entry.type !== 3 || entry.count < 1) return undefined;
  return new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength).getUint16(entry.valueOffset, true);
}

function readTagBytes(bytes, ifd, tag) {
  const entry = ifd.entries.find((e) => e.tag === tag);
  if (!entry || entry.valueOffset + entry.bytesLen > bytes.byteLength) return undefined;
  return { type: entry.type, count: entry.count, bytes: bytes.slice(entry.valueOffset, entry.valueOffset + entry.bytesLen) };
}

function parseTiffPalette(headerBytes) {
  try {
    const tiff = tiffAccess(headerBytes);
    const ifd = readTiffIfd(headerBytes, tiff.firstIfdOffset, tiff.big);
    if (readShortTag(headerBytes, ifd, 262) !== 3) return null;
    const colorMap = readTagBytes(headerBytes, ifd, 320);
    if (!colorMap || colorMap.type !== 3 || colorMap.count < 3 || colorMap.count % 3 !== 0) return null;
    return colorMap;
  } catch {
    return null;
  }
}

function writeTiffEntry(out, pos, big, tag, type, count, valueBytes, valueDataOffset) {
  const dv = new DataView(out.buffer);
  dv.setUint16(pos, tag, true);
  dv.setUint16(pos + 2, type, true);
  if (big) {
    dv.setBigUint64(pos + 4, BigInt(count), true);
    if (valueBytes.byteLength <= 8) out.set(valueBytes, pos + 12);
    else dv.setBigUint64(pos + 12, BigInt(valueDataOffset), true);
  } else {
    dv.setUint32(pos + 4, count, true);
    if (valueBytes.byteLength <= 4) out.set(valueBytes, pos + 8);
    else dv.setUint32(pos + 8, valueDataOffset, true);
  }
}

function patchTiffPalette(bytes, palette) {
  if (!palette) return bytes;
  const tiff = tiffAccess(bytes);
  const oldIfd = readTiffIfd(bytes, tiff.firstIfdOffset, tiff.big);
  const entries = oldIfd.entries.filter((e) => e.tag !== 320).sort((a, b) => a.tag - b.tag);
  const photo = entries.find((e) => e.tag === 262);
  if (!photo || photo.type !== 3 || photo.count !== 1) return bytes;

  const insertAt = entries.findIndex((e) => e.tag > 320);
  const paletteEntry = { tag: 320, type: palette.type, count: palette.count, bytes: palette.bytes };
  const ordered = entries.map((e) => ({ ...e, bytes: bytes.slice(e.valueOffset, e.valueOffset + e.bytesLen) }));
  ordered.splice(insertAt === -1 ? ordered.length : insertAt, 0, paletteEntry);

  const countBytes = tiff.big ? 8 : 2;
  const entryBytes = tiff.big ? 20 : 12;
  const nextBytes = tiff.big ? 8 : 4;
  const newIfdOffset = bytes.byteLength;
  const ifdBytes = countBytes + ordered.length * entryBytes + nextBytes;
  const extraStart = newIfdOffset + ifdBytes;
  let extraLen = 0;
  for (const e of ordered) {
    const inline = e.bytes.byteLength <= tiff.inlineBytes;
    if (!inline) extraLen += e.bytes.byteLength + (extraLen % 2);
  }

  const out = new Uint8Array(bytes.byteLength + ifdBytes + extraLen);
  out.set(bytes, 0);
  tiff.writeFirstIfd(out, newIfdOffset);
  const dv = new DataView(out.buffer);
  if (tiff.big) dv.setBigUint64(newIfdOffset, BigInt(ordered.length), true);
  else dv.setUint16(newIfdOffset, ordered.length, true);

  let extraOffset = extraStart;
  for (let i = 0; i < ordered.length; i++) {
    const e = ordered[i];
    const entryPos = newIfdOffset + countBytes + i * entryBytes;
    let valueBytes = e.bytes;
    if (e.tag === 262) valueBytes = new Uint8Array([3, 0]);
    let valueDataOffset = 0;
    if (valueBytes.byteLength > tiff.inlineBytes) {
      if ((extraOffset - extraStart) % 2) extraOffset++;
      valueDataOffset = extraOffset;
      out.set(valueBytes, valueDataOffset);
      extraOffset += valueBytes.byteLength;
    }
    writeTiffEntry(out, entryPos, tiff.big, e.tag, e.type, e.count, valueBytes, valueDataOffset);
  }

  const nextPos = newIfdOffset + countBytes + ordered.length * entryBytes;
  if (tiff.big) dv.setBigUint64(nextPos, BigInt(oldIfd.nextOffset), true);
  else dv.setUint32(nextPos, oldIfd.nextOffset, true);
  return out;
}

function selectLevelForResolution(levels, gt, datasetBbox, bbox, resolution) {
  if (resolution == null) return 0;
  if (!Number.isFinite(resolution) || resolution <= 0) {
    throw new Error("resolution must be a positive number");
  }

  const datasetWidth = Math.abs(datasetBbox[2] - datasetBbox[0]);
  const datasetHeight = Math.abs(datasetBbox[3] - datasetBbox[1]);
  const bboxWidth = Math.abs(bbox[2] - bbox[0]);
  const bboxHeight = Math.abs(bbox[3] - bbox[1]);
  const scaleX = bboxWidth > 0 ? datasetWidth / bboxWidth : 1;
  const scaleY = bboxHeight > 0 ? datasetHeight / bboxHeight : 1;
  const targetX = resolution * scaleX;
  const targetY = resolution * scaleY;

  let best = 0;
  let bestScore = Infinity;
  for (let i = 0; i < levels.length; i++) {
    const level = levels[i];
    const levelScaleX = level.width / levels[0].width;
    const levelScaleY = level.height / levels[0].height;
    const px = Math.abs(gt[1] / levelScaleX);
    const py = Math.abs(gt[5] / levelScaleY);
    const score = Math.abs(Math.log(px / targetX)) + Math.abs(Math.log(py / targetY));
    if (score < bestScore) {
      best = i;
      bestScore = score;
    }
  }
  return best;
}

function outputBboxForCrs(bbox, bboxCrs, outputCrs) {
  if (outputCrs == null || outputCrs === bboxCrs) return bbox.slice();
  return Array.from(transform_bbox_epsg(bboxCrs, outputCrs, bbox));
}

function reprojectSubsetNearest(stream, source, src, dst, outputCrs, nodata) {
  const fill = nodata ?? Number.NaN;
  const out = new Float64Array(dst.width * dst.height * src.bands);
  out.fill(fill);

  const batchRows = 32;
  for (let row0 = 0; row0 < dst.height; row0 += batchRows) {
    const row1 = Math.min(dst.height, row0 + batchRows);
    const coords = new Array((row1 - row0) * dst.width * 2);
    let k = 0;
    for (let row = row0; row < row1; row++) {
      const y = dst.y0 + (row + 0.5) * dst.pixelHeight;
      for (let col = 0; col < dst.width; col++) {
        coords[k++] = dst.x0 + (col + 0.5) * dst.pixelWidth;
        coords[k++] = y;
      }
    }

    const srcCoords = stream.points_to_dataset_crs(outputCrs, coords);
    k = 0;
    for (let row = row0; row < row1; row++) {
      for (let col = 0; col < dst.width; col++) {
        const x = srcCoords[k++];
        const y = srcCoords[k++];
        if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
        const srcCol = Math.floor((x - src.x0) / src.pixelWidth);
        const srcRow = Math.floor((y - src.y0) / src.pixelHeight);
        if (srcCol < 0 || srcRow < 0 || srcCol >= src.width || srcRow >= src.height) continue;
        const srcPixel = (srcRow * src.width + srcCol) * src.bands;
        const dstPixel = ((row * dst.width) + col) * src.bands;
        for (let band = 0; band < src.bands; band++) out[dstPixel + band] = source[srcPixel + band];
      }
    }
  }

  return out;
}

function readGeotiffInterleaved(reader) {
  const bands = reader.bands;
  const pixels = reader.width * reader.height;
  const out = new Float64Array(pixels * bands);
  for (let band = 0; band < bands; band++) {
    const values = reader.read_band_f64(band);
    for (let i = 0; i < pixels; i++) out[i * bands + band] = values[i];
  }
  return out;
}

function writeTypedCog({ data, width, height, bands, sampleFormat, bitsPerSample, geoTransform, epsg, nodata, palette }) {
  const builder = new CogBuilder(width, height, bands);
  builder.set_geo_transform(geoTransform);
  builder.set_compression("deflate");
  if (epsg != null) builder.set_epsg(epsg);
  if (nodata != null) builder.set_nodata(nodata);
  if (sampleFormat === "uint" && bitsPerSample === 8) {
    const u8 = new Uint8Array(data.length);
    const fill = nodata == null ? 0 : Math.max(0, Math.min(255, Math.round(nodata)));
    for (let i = 0; i < data.length; i++) {
      const v = data[i];
      u8[i] = Number.isFinite(v) ? Math.max(0, Math.min(255, Math.round(v))) : fill;
    }
    const bytes = builder.write_u8(u8);
    return bands === 1 && palette ? patchTiffPalette(bytes, palette) : bytes;
  }
  if (sampleFormat === "ieeefloat" && bitsPerSample === 32) {
    return builder.write_f32(Float32Array.from(data));
  }
  if (sampleFormat === "ieeefloat" && bitsPerSample === 64) {
    return builder.write_f64(data);
  }
  throw new Error(`preserving source sample type is not yet supported for ${sampleFormat}/${bitsPerSample}-bit rasters`);
}

function windowFromBbox(gt, baseLevel, level, bbox) {
  const [x0, pixelWidth, rowRotation, y0, colRotation, pixelHeight] = gt;
  if (Math.abs(rowRotation) > 1e-12 || Math.abs(colRotation) > 1e-12) {
    throw new Error("rotated/skewed COG geo-transforms are not supported");
  }
  if (!(pixelWidth > 0) || !(pixelHeight < 0)) {
    throw new Error("only north-up COGs with positive pixel width and negative pixel height are supported");
  }

  const scaleX = level.width / baseLevel.width;
  const scaleY = level.height / baseLevel.height;
  const px = pixelWidth / scaleX;
  const py = pixelHeight / scaleY;

  const minCol = Math.floor((bbox[0] - x0) / px);
  const maxCol = Math.ceil((bbox[2] - x0) / px);
  const minRow = Math.floor((bbox[3] - y0) / py);
  const maxRow = Math.ceil((bbox[1] - y0) / py);

  const x = Math.max(0, Math.min(level.width, minCol));
  const y = Math.max(0, Math.min(level.height, minRow));
  const x2 = Math.max(0, Math.min(level.width, maxCol));
  const y2 = Math.max(0, Math.min(level.height, maxRow));
  if (x2 <= x || y2 <= y) throw new Error("bbox does not intersect the COG extent");
  return { x, y, width: x2 - x, height: y2 - y, pixelWidth: px, pixelHeight: py };
}

function dimensionsForBbox(bbox, resolution, width, height) {
  const bboxWidth = bbox[2] - bbox[0];
  const bboxHeight = bbox[3] - bbox[1];
  if (!(bboxWidth > 0) || !(bboxHeight > 0)) throw new Error("bbox dimensions must be positive");
  if (resolution != null) {
    if (!(resolution > 0)) throw new Error("resolution must be positive");
    return {
      width: Math.max(1, Math.ceil(bboxWidth / resolution)),
      height: Math.max(1, Math.ceil(bboxHeight / resolution)),
    };
  }
  if (width != null && height != null) return { width, height };
  if (width != null) return { width, height: Math.max(1, Math.round(width * bboxHeight / bboxWidth)) };
  if (height != null) return { width: Math.max(1, Math.round(height * bboxWidth / bboxHeight)), height };
  const maxDim = 1024;
  return bboxWidth >= bboxHeight
    ? { width: maxDim, height: Math.max(1, Math.round(maxDim * bboxHeight / bboxWidth)) }
    : { width: Math.max(1, Math.round(maxDim * bboxWidth / bboxHeight)), height: maxDim };
}

const WEB_MERCATOR_EXTENT = 20037508.342789244;

function lonLatToTile(lon, lat, zoom) {
  const n = 2 ** zoom;
  const clampedLat = Math.max(-85.05112878, Math.min(85.05112878, lat));
  const latRad = clampedLat * Math.PI / 180;
  return {
    x: Math.floor(((lon + 180) / 360) * n),
    y: Math.floor((1 - Math.log(Math.tan(latRad) + 1 / Math.cos(latRad)) / Math.PI) / 2 * n),
  };
}

function tileMercatorBounds(x, y, zoom) {
  const n = 2 ** zoom;
  const span = (WEB_MERCATOR_EXTENT * 2) / n;
  const minX = -WEB_MERCATOR_EXTENT + x * span;
  const maxY = WEB_MERCATOR_EXTENT - y * span;
  return [minX, maxY - span, minX + span, maxY];
}

function xyzTileUrl(template, x, y, z, subdomains) {
  const s = subdomains ? subdomains[Math.abs(x + y + z) % subdomains.length] : "";
  return template
    .replaceAll("{x}", String(x))
    .replaceAll("{y}", String(y))
    .replaceAll("{z}", String(z))
    .replaceAll("{s}", s);
}

function canvasFor(width, height) {
  if (typeof OffscreenCanvas !== "undefined") return new OffscreenCanvas(width, height);
  if (typeof document !== "undefined" && document.createElement) {
    const canvas = document.createElement("canvas");
    canvas.width = width;
    canvas.height = height;
    return canvas;
  }
  throw new Error("extract_xyz_tile_subset requires browser image decoding support");
}

async function decodeImageRgba(bytes) {
  if (typeof createImageBitmap === "undefined" || typeof Blob === "undefined") {
    throw new Error("extract_xyz_tile_subset requires createImageBitmap support to decode PNG/JPEG tiles");
  }
  const bitmap = await createImageBitmap(new Blob([bytes]));
  const canvas = canvasFor(bitmap.width, bitmap.height);
  const ctx = canvas.getContext("2d", { willReadFrequently: true });
  ctx.drawImage(bitmap, 0, 0);
  const image = ctx.getImageData(0, 0, bitmap.width, bitmap.height);
  if (bitmap.close) bitmap.close();
  return { width: image.width, height: image.height, data: image.data };
}

async function fetchTileRgba(url, fetchOptions) {
  const resp = await fetch(url, fetchOptions);
  if (!resp.ok) throw new Error(`tile request failed (${resp.status}): ${url}`);
  return decodeImageRgba(new Uint8Array(await resp.arrayBuffer()));
}

function buildWmsGetMapUrl(endpoint, opts) {
  const u = new URL(endpoint);
  const version = opts.version || "1.1.1";
  const crsParam = version.startsWith("1.3") ? "CRS" : "SRS";
  const params = {
    SERVICE: "WMS",
    VERSION: version,
    REQUEST: "GetMap",
    LAYERS: opts.layers,
    STYLES: opts.styles ?? "",
    FORMAT: opts.format,
    TRANSPARENT: "FALSE",
    WIDTH: String(opts.width),
    HEIGHT: String(opts.height),
    BBOX: opts.bbox.join(","),
    [crsParam]: `EPSG:${opts.crs}`,
  };
  for (const [key, value] of Object.entries(params)) u.searchParams.set(key, value);
  return u.toString();
}

/**
 * Request a bbox subset from a WMS GetMap endpoint and write it as a Deflate COG.
 *
 * The WMS response must be a GeoTIFF-compatible format, typically
 * `image/geotiff` or `image/tiff`. The request CRS defaults to `bboxCrs` unless
 * `outputCrs` is supplied.
 *
 * @param {string} url WMS endpoint URL.
 * @param {object} opts
 * @param {string} opts.layers WMS layer name(s), comma-separated.
 * @param {string} [opts.styles] WMS style name(s), comma-separated.
 * @param {[number, number, number, number]} opts.bbox [minX,minY,maxX,maxY].
 * @param {number} opts.bboxCrs EPSG code of `bbox`.
 * @param {number} [opts.resolution] Target output pixel size in output CRS units.
 * @param {number} [opts.width] Request width in pixels; used when resolution is omitted.
 * @param {number} [opts.height] Request height in pixels; used when resolution is omitted.
 * @param {number} [opts.outputCrs] Optional output/request EPSG code.
 * @param {string} [opts.format="image/geotiff"] WMS response format.
 * @param {string} [opts.version="1.1.1"] WMS version.
 * @param {number} [opts.nodata] Optional output nodata value.
 * @param {RequestInit} [opts.fetchOptions] Extra fetch options for the GetMap request.
 * @returns {Promise<Uint8Array>}
 */
export async function extractWmsSubset(url, opts) {
  opts = opts || {};
  await initLibrary();
  const { layers, styles, bbox, bboxCrs, resolution, width, height, nodata } = opts;
  const outputCrs = opts.outputCrs ?? bboxCrs;
  const format = opts.format || "image/geotiff";
  const version = opts.version || "1.1.1";
  if (!/^https?:\/\//i.test(url)) throw new Error(`url must be HTTP(S), got: ${url}`);
  if (!layers || String(layers).trim() === "") throw new Error("opts.layers is required");
  if (!/tiff|geotiff|geotif/i.test(format)) {
    throw new Error("extract_wms_subset requires a GeoTIFF WMS format such as image/geotiff or image/tiff");
  }
  if (!Array.isArray(bbox) || bbox.length !== 4 || bbox.some((v) => !Number.isFinite(v)) || bbox[0] >= bbox[2] || bbox[1] >= bbox[3]) {
    throw new Error("opts.bbox must be finite and ordered [minX,minY,maxX,maxY]");
  }
  if (!Number.isInteger(bboxCrs) || bboxCrs <= 0) throw new Error("opts.bboxCrs must be a positive EPSG code");
  if (!Number.isInteger(outputCrs) || outputCrs <= 0) throw new Error("opts.outputCrs must be a positive EPSG code");
  if (nodata != null && !Number.isFinite(nodata)) throw new Error("opts.nodata must be a finite number");

  const requestBbox = outputBboxForCrs(bbox, bboxCrs, outputCrs);
  const dims = dimensionsForBbox(requestBbox, resolution, width, height);
  const requestUrl = buildWmsGetMapUrl(url, {
    layers: String(layers),
    styles: styles == null ? "" : String(styles),
    bbox: requestBbox,
    crs: outputCrs,
    width: dims.width,
    height: dims.height,
    format,
    version,
  });
  const resp = await fetch(requestUrl, opts.fetchOptions);
  const bytes = new Uint8Array(await resp.arrayBuffer());
  const contentType = resp.headers.get("content-type") || "";
  if (!resp.ok) throw new Error(`WMS GetMap failed (${resp.status}): ${new TextDecoder().decode(bytes.slice(0, 512))}`);
  if (/xml|text/i.test(contentType)) {
    throw new Error(`WMS returned an exception instead of GeoTIFF: ${new TextDecoder().decode(bytes.slice(0, 1024))}`);
  }

  const reader = new GeoTiffReader(bytes);
  const data = readGeotiffInterleaved(reader);
  const palette = parseTiffPalette(bytes);
  return writeTypedCog({
    data,
    width: reader.width,
    height: reader.height,
    bands: reader.bands,
    sampleFormat: reader.sample_format,
    bitsPerSample: reader.bits_per_sample,
    geoTransform: [requestBbox[0], (requestBbox[2] - requestBbox[0]) / reader.width, 0, requestBbox[3], 0, -(requestBbox[3] - requestBbox[1]) / reader.height],
    epsg: outputCrs,
    nodata: nodata ?? reader.nodata,
    palette,
  });
}

/**
 * Fetch XYZ raster tiles for a bbox, mosaic them, and write a Deflate RGB COG.
 *
 * Tiles are assumed to be in the standard XYZ Web Mercator grid (EPSG:3857).
 * The output CRS defaults to `bboxCrs`.
 *
 * @param {string} url XYZ tile URL template with `{z}`, `{x}`, `{y}` and optional `{s}`.
 * @param {object} opts
 * @param {number} opts.zoom XYZ zoom level.
 * @param {[number, number, number, number]} opts.bbox [minX,minY,maxX,maxY].
 * @param {number} opts.bboxCrs EPSG code of `bbox`.
 * @param {number} [opts.resolution] Target output pixel size in output CRS units.
 * @param {number} [opts.width] Output width in pixels; used when resolution is omitted.
 * @param {number} [opts.height] Output height in pixels; used when resolution is omitted.
 * @param {number} [opts.outputCrs] Optional output EPSG code; defaults to bboxCrs.
 * @param {number} [opts.tileSize=256] Tile size in pixels.
 * @param {string} [opts.subdomains] Optional subdomain letters for `{s}`.
 * @param {number} [opts.nodata] Optional output nodata/fill value for missing or transparent pixels.
 * @param {RequestInit} [opts.fetchOptions] Extra fetch options for tile requests.
 * @returns {Promise<Uint8Array>}
 */
export async function extractXyzTileSubset(url, opts) {
  opts = opts || {};
  await initLibrary();
  const { zoom, bbox, bboxCrs, resolution, width, height, nodata, subdomains } = opts;
  const outputCrs = opts.outputCrs ?? bboxCrs;
  const tileSize = opts.tileSize ?? 256;
  if (!/^https?:\/\//i.test(url)) throw new Error(`url must be HTTP(S), got: ${url}`);
  if (!Number.isInteger(zoom) || zoom < 0 || zoom > 30) throw new Error("opts.zoom must be an integer 0-30");
  if (!Number.isInteger(tileSize) || tileSize <= 0) throw new Error("opts.tileSize must be a positive integer");
  if (!Array.isArray(bbox) || bbox.length !== 4 || bbox.some((v) => !Number.isFinite(v)) || bbox[0] >= bbox[2] || bbox[1] >= bbox[3]) {
    throw new Error("opts.bbox must be finite and ordered [minX,minY,maxX,maxY]");
  }
  if (!Number.isInteger(bboxCrs) || bboxCrs <= 0) throw new Error("opts.bboxCrs must be a positive EPSG code");
  if (!Number.isInteger(outputCrs) || outputCrs <= 0) throw new Error("opts.outputCrs must be a positive EPSG code");
  if (nodata != null && !Number.isFinite(nodata)) throw new Error("opts.nodata must be a finite number");

  const bbox4326 = bboxCrs === 4326 ? bbox : Array.from(transform_bbox_epsg(bboxCrs, 4326, bbox));
  const n = 2 ** zoom;
  const nw = lonLatToTile(bbox4326[0], bbox4326[3], zoom);
  const se = lonLatToTile(bbox4326[2], bbox4326[1], zoom);
  const minTileX = Math.max(0, Math.min(n - 1, nw.x));
  const maxTileX = Math.max(0, Math.min(n - 1, se.x));
  const minTileY = Math.max(0, Math.min(n - 1, nw.y));
  const maxTileY = Math.max(0, Math.min(n - 1, se.y));
  if (maxTileX < minTileX || maxTileY < minTileY) throw new Error("bbox does not intersect the XYZ tile grid");

  const tileCols = maxTileX - minTileX + 1;
  const tileRows = maxTileY - minTileY + 1;
  if (tileCols * tileRows > 512) throw new Error(`bbox intersects too many tiles at zoom ${zoom}: ${tileCols * tileRows}`);
  const mosaicWidth = tileCols * tileSize;
  const mosaicHeight = tileRows * tileSize;
  const fill = nodata == null ? 0 : Math.max(0, Math.min(255, Math.round(nodata)));
  const mosaic = new Uint8Array(mosaicWidth * mosaicHeight * 4);
  for (let i = 0; i < mosaic.length; i += 4) {
    mosaic[i] = fill;
    mosaic[i + 1] = fill;
    mosaic[i + 2] = fill;
    mosaic[i + 3] = 0;
  }

  for (let ty = minTileY; ty <= maxTileY; ty++) {
    for (let tx = minTileX; tx <= maxTileX; tx++) {
      const tileUrl = xyzTileUrl(url, tx, ty, zoom, subdomains);
      const tile = await fetchTileRgba(tileUrl, opts.fetchOptions);
      const copyW = Math.min(tileSize, tile.width);
      const copyH = Math.min(tileSize, tile.height);
      const dstX = (tx - minTileX) * tileSize;
      const dstY = (ty - minTileY) * tileSize;
      for (let row = 0; row < copyH; row++) {
        for (let col = 0; col < copyW; col++) {
          const src = (row * tile.width + col) * 4;
          const dst = ((dstY + row) * mosaicWidth + dstX + col) * 4;
          mosaic[dst] = tile.data[src];
          mosaic[dst + 1] = tile.data[src + 1];
          mosaic[dst + 2] = tile.data[src + 2];
          mosaic[dst + 3] = tile.data[src + 3];
        }
      }
    }
  }

  const mb0 = tileMercatorBounds(minTileX, maxTileY, zoom);
  const mb1 = tileMercatorBounds(maxTileX, minTileY, zoom);
  const mosaicBbox3857 = [mb0[0], mb0[1], mb1[2], mb1[3]];
  const mosaicPx = (mosaicBbox3857[2] - mosaicBbox3857[0]) / mosaicWidth;
  const mosaicPy = -(mosaicBbox3857[3] - mosaicBbox3857[1]) / mosaicHeight;

  const outBbox = outputBboxForCrs(bbox, bboxCrs, outputCrs);
  const nativeResolution = outputCrs === 3857
    ? Math.abs(mosaicPx)
    : Math.max((outBbox[2] - outBbox[0]) / Math.max(1, Math.round((transform_bbox_epsg(outputCrs, 3857, outBbox)[2] - transform_bbox_epsg(outputCrs, 3857, outBbox)[0]) / Math.abs(mosaicPx))), 1e-12);
  const dims = dimensionsForBbox(outBbox, resolution ?? nativeResolution, width, height);
  const out = new Float64Array(dims.width * dims.height * 3);
  const outPx = (outBbox[2] - outBbox[0]) / dims.width;
  const outPy = -(outBbox[3] - outBbox[1]) / dims.height;

  const batchRows = 32;
  for (let row0 = 0; row0 < dims.height; row0 += batchRows) {
    const row1 = Math.min(dims.height, row0 + batchRows);
    const coords = new Array((row1 - row0) * dims.width * 2);
    let k = 0;
    for (let row = row0; row < row1; row++) {
      const y = outBbox[3] + (row + 0.5) * outPy;
      for (let col = 0; col < dims.width; col++) {
        coords[k++] = outBbox[0] + (col + 0.5) * outPx;
        coords[k++] = y;
      }
    }
    const merc = outputCrs === 3857 ? coords : Array.from(transform_points_epsg(outputCrs, 3857, coords));
    k = 0;
    for (let row = row0; row < row1; row++) {
      for (let col = 0; col < dims.width; col++) {
        const x = merc[k++];
        const y = merc[k++];
        const srcCol = Math.floor((x - mosaicBbox3857[0]) / mosaicPx);
        const srcRow = Math.floor((y - mosaicBbox3857[3]) / mosaicPy);
        const dst = (row * dims.width + col) * 3;
        if (srcCol < 0 || srcRow < 0 || srcCol >= mosaicWidth || srcRow >= mosaicHeight) {
          out[dst] = fill; out[dst + 1] = fill; out[dst + 2] = fill;
          continue;
        }
        const src = (srcRow * mosaicWidth + srcCol) * 4;
        if (mosaic[src + 3] === 0) {
          out[dst] = fill; out[dst + 1] = fill; out[dst + 2] = fill;
        } else {
          out[dst] = mosaic[src];
          out[dst + 1] = mosaic[src + 1];
          out[dst + 2] = mosaic[src + 2];
        }
      }
    }
  }

  return writeTypedCog({
    data: out,
    width: dims.width,
    height: dims.height,
    bands: 3,
    sampleFormat: "uint",
    bitsPerSample: 8,
    geoTransform: [outBbox[0], outPx, 0, outBbox[3], 0, outPy],
    epsg: outputCrs,
    nodata,
  });
}

/**
 * Extract a bbox subset from a local or HTTP Cloud Optimized GeoTIFF. HTTP
 * sources are read with byte-range requests, without downloading the full COG.
 *
 * The returned bytes are a new COG containing all bands from the selected
 * source level, preserving supported source sample types. `bboxCrs` is an EPSG
 * code for `bbox`; it is reprojected to the COG CRS before selecting tiles. If
 * `resolution` is set and `level` is omitted, the closest available COG overview
 * level is selected. If `outputCrs` is set, the extracted source window is
 * reprojected to that EPSG CRS with nearest-neighbor resampling. Sources with
 * user-defined projection strings default to `bboxCrs` output so the result can
 * be written with standard EPSG metadata.
 *
 * @param {string|Uint8Array|ArrayBuffer} source HTTP(S) COG URL or local COG bytes.
 * @param {object} opts
 * @param {[number, number, number, number]} opts.bbox [minX,minY,maxX,maxY].
 * @param {number} opts.bboxCrs EPSG code of `bbox`.
 * @param {number} [opts.level] COG overview level to read; 0 is full res.
 * @param {number} [opts.resolution] Target output pixel size in outputCrs units when outputCrs is set; otherwise bboxCrs units.
 * @param {number} [opts.outputCrs] Optional output EPSG code.
 * @param {number} [opts.nodata] Optional output nodata value.
 * @param {RequestInit} [opts.fetchOptions] Extra fetch options for all requests.
 * @param {number} [opts.initialHeaderBytes=262144] Initial COG header prefix size.
 * @param {number} [opts.maxHeaderBytes=8388608] Maximum COG header prefix size.
 * @returns {Promise<Uint8Array>}
 */
export async function extractCogSubset(source, opts) {
  opts = opts || {};
  await initLibrary();
  const { bbox, bboxCrs, resolution, nodata } = opts || {};
  let { level, outputCrs } = opts || {};
  if (!Array.isArray(bbox) || bbox.length !== 4) {
    throw new Error("opts.bbox must be [minX,minY,maxX,maxY]");
  }
  if (bbox.some((v) => !Number.isFinite(v)) || bbox[0] >= bbox[2] || bbox[1] >= bbox[3]) {
    throw new Error("opts.bbox must be finite and ordered min < max");
  }
  if (!Number.isInteger(bboxCrs) || bboxCrs <= 0) {
    throw new Error("opts.bboxCrs must be a positive EPSG code");
  }
  if (level != null && (!Number.isInteger(level) || level < 0)) {
    throw new Error("opts.level must be a non-negative integer");
  }
  if (outputCrs != null && (!Number.isInteger(outputCrs) || outputCrs <= 0)) {
    throw new Error("opts.outputCrs must be a positive EPSG code");
  }
  if (nodata != null && !Number.isFinite(nodata)) {
    throw new Error("opts.nodata must be a finite number");
  }

  const reader = makeSourceReader(source, opts.fetchOptions);
  const { stream, header } = await openCogStream(reader, opts || {});
  const levels = parseLevels(stream);
  const sourcePalette = parseTiffPalette(header);

  if (outputCrs == null && stream.has_projection_string) outputCrs = bboxCrs;
  const datasetBbox = Array.from(stream.bbox_to_dataset_crs(bboxCrs, bbox));
  const requestedOutputBbox = outputBboxForCrs(bbox, bboxCrs, outputCrs);
  const gt = Array.from(stream.geo_transform());
  if (gt.length !== 6) throw new Error("COG has no affine geo-transform");
  if (level == null) level = selectLevelForResolution(levels, gt, datasetBbox, requestedOutputBbox, resolution);
  const selected = levels[level];
  if (!selected) throw new Error(`level ${level} out of range`);
  const win = windowFromBbox(gt, levels[0], selected, datasetBbox);

  const tileSpecs = JSON.parse(stream.tiles_for_window(level, win.x, win.y, win.width, win.height));
  const out = new Float64Array(win.width * win.height * selected.bands);
  const tileStride = selected.tile_width * selected.tile_height * selected.bands;

  for (const tile of tileSpecs) {
    const bytes = await reader.range(tile.offset, tile.length);
    const decoded = stream.decode_tile_f64(level, bytes);
    if (decoded.length !== tileStride) {
      throw new Error(`decoded tile size mismatch for tile ${tile.col},${tile.row}`);
    }

    const tileX0 = tile.col * selected.tile_width;
    const tileY0 = tile.row * selected.tile_height;
    const copyX0 = Math.max(win.x, tileX0);
    const copyY0 = Math.max(win.y, tileY0);
    const copyX1 = Math.min(win.x + win.width, tileX0 + selected.tile_width, selected.width);
    const copyY1 = Math.min(win.y + win.height, tileY0 + selected.tile_height, selected.height);

    for (let row = copyY0; row < copyY1; row++) {
      for (let col = copyX0; col < copyX1; col++) {
        const srcPixel = ((row - tileY0) * selected.tile_width + (col - tileX0)) * selected.bands;
        const dstPixel = ((row - win.y) * win.width + (col - win.x)) * selected.bands;
        for (let band = 0; band < selected.bands; band++) {
          out[dstPixel + band] = decoded[srcPixel + band];
        }
      }
    }
  }

  const subsetX0 = gt[0] + win.x * win.pixelWidth;
  const subsetY0 = gt[3] + win.y * win.pixelHeight;

  let finalData = out;
  let finalWidth = win.width;
  let finalHeight = win.height;
  let finalX0 = subsetX0;
  let finalY0 = subsetY0;
  let finalPixelWidth = win.pixelWidth;
  let finalPixelHeight = win.pixelHeight;
  let finalEpsg = !stream.has_projection_string ? stream.epsg : undefined;
  const outputNodata = nodata ?? stream.nodata;

  if (outputCrs != null) {
    const outWidth = resolution == null
      ? win.width
      : Math.max(1, Math.ceil((requestedOutputBbox[2] - requestedOutputBbox[0]) / resolution));
    const outHeight = resolution == null
      ? win.height
      : Math.max(1, Math.ceil((requestedOutputBbox[3] - requestedOutputBbox[1]) / resolution));
    const outPixelWidth = (requestedOutputBbox[2] - requestedOutputBbox[0]) / outWidth;
    const outPixelHeight = -(requestedOutputBbox[3] - requestedOutputBbox[1]) / outHeight;

    finalData = reprojectSubsetNearest(
      stream,
      out,
      { x0: subsetX0, y0: subsetY0, pixelWidth: win.pixelWidth, pixelHeight: win.pixelHeight, width: win.width, height: win.height, bands: selected.bands },
      { x0: requestedOutputBbox[0], y0: requestedOutputBbox[3], pixelWidth: outPixelWidth, pixelHeight: outPixelHeight, width: outWidth, height: outHeight },
      outputCrs,
      outputNodata,
    );
    finalWidth = outWidth;
    finalHeight = outHeight;
    finalX0 = requestedOutputBbox[0];
    finalY0 = requestedOutputBbox[3];
    finalPixelWidth = outPixelWidth;
    finalPixelHeight = outPixelHeight;
    finalEpsg = outputCrs;
  }

  return writeTypedCog({
    data: finalData,
    width: finalWidth,
    height: finalHeight,
    bands: selected.bands,
    sampleFormat: selected.sample_format,
    bitsPerSample: selected.bits_per_sample,
    geoTransform: [finalX0, finalPixelWidth, 0, finalY0, 0, finalPixelHeight],
    epsg: finalEpsg,
    nodata: outputNodata,
    palette: sourcePalette,
  });
}
