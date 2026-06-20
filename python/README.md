# geolibre-wasm (Python)

Run the [`geolibre-rust`](https://github.com/opengeos/geolibre-rust) geospatial
tool suite (the `whitebox_next_gen` tools plus GeoLibre's own) from Python. The
tools are a single WebAssembly (WASI) module executed in-process via
[`wasmtime`](https://github.com/bytecodealliance/wasmtime-py), so there is **no
native install, no GDAL, and no server** — just `pip install`.

This mirrors the JavaScript `geolibre-wasm/tools` API (`list_tools`,
`list_manifests`, `run_tool`), so the two stay in sync.

> The import package is `geolibre_wasm` (the distribution is `geolibre-wasm`),
> matching the npm package name and avoiding a clash with the separate
> `geolibre` application package.

## Install

```bash
pip install geolibre-wasm
```

On first use the runtime (`geolibre-cli.wasm`, ~20 MB) is downloaded from the
matching GitHub release and cached under `~/.cache/geolibre/`. To use a local
copy instead, set `GEOLIBRE_WASM=/path/to/geolibre-cli.wasm` or pass
`wasm_path=` to any call.

## Usage

Inputs are passed as `bytes` under `/work`; the tool reads/writes there and any
new files come back as `bytes`.

```python
import geolibre_wasm as gl

# Discover tools (each manifest carries a "source": "geolibre" | "whitebox")
tools = gl.list_tools()
manifests = gl.list_manifests()

# Raster: compute slope from a DEM
dem = open("dem.tif", "rb").read()
res = gl.run_tool(
    "slope",
    args=["--input=/work/dem.tif", "--output=/work/slope.tif", "--units=degrees"],
    input={"dem.tif": dem},
)
assert res.exit_code == 0
open("slope.tif", "wb").write(res.files["slope.tif"])

# Reproject (warp) to a target EPSG
res = gl.run_tool(
    "reproject_raster",
    args=["--input=/work/dem.tif", "--output=/work/wgs84.tif", "--epsg=4326"],
    input={"dem.tif": dem},
)

# Vector: GeoJSON -> GeoParquet (Hilbert-sorted, bbox covering, ZSTD by default)
gj = open("cities.geojson", "rb").read()
res = gl.run_tool(
    "write_geoparquet",
    args=["--input=/work/in.geojson", "--output=/work/out.parquet"],
    input={"in.geojson": gj},
)
open("cities.parquet", "wb").write(res.files["out.parquet"])
```

Tools that write a directory tree (e.g. `raster_to_tiles`) return nested keys:

```python
res = gl.run_tool(
    "raster_to_tiles",
    args=["--input=/work/dem.tif", "--output_dir=/work/tiles", "--min_zoom=16", "--max_zoom=18"],
    input={"dem.tif": dem},
)
for path, data in res.files.items():
    # e.g. "tiles/16/9559/32767.png"
    ...
```

## API

- `list_tools(wasm_path=None) -> list[str]`
- `list_manifests(wasm_path=None) -> list[dict]`
- `run_tool(tool, args=None, input=None, wasm_path=None) -> ToolResult`
- `ToolResult(exit_code: int, stdout: list[str], files: dict[str, bytes])`
- `runtime_path(wasm_path=None) -> str` — resolve the runtime (explicit > `GEOLIBRE_WASM` > cached download)
- `download_runtime(dest=None) -> str` — fetch the runtime ahead of time

## License

MIT
