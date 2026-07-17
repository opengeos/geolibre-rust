"""GeoLibre: run the whitebox_next_gen + GeoLibre geospatial tool suite from Python.

The tools are compiled to a single WebAssembly (WASI) module and executed
in-process via wasmtime, so there is no native install, GDAL, or server. The API
mirrors the JavaScript ``geolibre-wasm/tools`` package.

Example:
    >>> import geolibre_wasm as gl
    >>> dem = open("dem.tif", "rb").read()
    >>> result = gl.run_tool(
    ...     "slope",
    ...     args=["--input=/work/dem.tif", "--output=/work/slope.tif", "--units=degrees"],
    ...     input={"dem.tif": dem},
    ... )
    >>> result.exit_code
    0
    >>> open("slope.tif", "wb").write(result.files["slope.tif"])
"""

from ._core import (
    RUNTIME_VERSION,
    ToolResult,
    download_runtime,
    list_manifests,
    list_tools,
    run_tool,
    runtime_path,
)

__all__ = [
    "RUNTIME_VERSION",
    "ToolResult",
    "download_runtime",
    "list_manifests",
    "list_tools",
    "run_tool",
    "runtime_path",
]

__version__ = "0.9.0"
