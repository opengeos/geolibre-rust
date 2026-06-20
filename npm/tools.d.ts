/** Result of running a tool. */
export interface ToolResult {
  /** Process exit code (0 = success). */
  exitCode: number;
  /** Captured stdout/stderr lines. */
  stdout: string[];
  /** New files the tool wrote, keyed by path relative to /work. Top-level
   *  outputs use their basename (e.g. "slope.tif"); tools that write a tree
   *  (e.g. raster_to_tiles) use nested keys like "tiles/15/4779/16383.png". */
  files: Record<string, Uint8Array>;
}

export interface RunToolOptions {
  /** CLI args, e.g. ["--input=/work/dem.tif", "--output=/work/out.tif", "--units=degrees"]. */
  args?: string[];
  /** Input files placed under /work, keyed by filename. */
  input?: Record<string, Uint8Array>;
}

/** A single parameter in a tool manifest. */
export interface ToolParam {
  name: string;
  [key: string]: unknown;
}

/** A tool's metadata and parameter schema. */
export interface ToolManifest {
  id: string;
  display_name: string;
  summary: string;
  params: ToolParam[];
  [key: string]: unknown;
}

/** Compile the WASI tool runner once. Omit `source` in browsers/bundlers; pass
 *  the wasm bytes or a URL/Response in Node. */
export function initTools(source?: URL | Response | BufferSource | string): Promise<WebAssembly.Module>;

/** List every available tool id. */
export function listTools(): Promise<string[]>;

/** Fetch every tool manifest (parameter schemas), for building UIs offline. */
export function listManifests(): Promise<ToolManifest[]>;

/** Run one tool over an in-memory filesystem. */
export function runTool(tool: string, opts?: RunToolOptions): Promise<ToolResult>;
