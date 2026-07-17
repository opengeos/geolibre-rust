"""Run the GeoLibre WASI tool runner from Python via an in-process wasmtime.

The tools are the single ``geolibre-cli.wasm`` (a ``wasm32-wasip1`` module). Each
call runs it over a private temporary directory preopened as ``/work``: input
files are written there, the tool reads/writes via ordinary ``std::fs``, and any
new files are returned as ``bytes``. This mirrors the JavaScript ``tools.mjs``
API (``list_tools`` / ``list_manifests`` / ``run_tool``) so the two stay in sync.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import urllib.request
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

import wasmtime

#: Magic header of every WebAssembly module ("\0asm").
_WASM_MAGIC = b"\x00asm"
#: Network timeout (seconds) for the one-time runtime download.
_DOWNLOAD_TIMEOUT = 120
#: Sent on HTTP fetches; some CDNs reject the default ``Python-urllib`` agent.
_USER_AGENT = "Mozilla/5.0 (geolibre-wasm)"


def _http_get(url: str, timeout: int) -> bytes:
    """GET a URL to bytes with a browser-like User-Agent and a timeout."""
    request = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310
        return response.read()


#: Release whose ``geolibre-cli.wasm`` asset this wrapper downloads by default.
#: Kept in sync with the package version's ``vMAJOR.MINOR.PATCH`` tag.
RUNTIME_VERSION = "v0.9.0"
_ASSET = f"geolibre-cli-{RUNTIME_VERSION}.wasm"
_RELEASE_URL = (
    "https://github.com/opengeos/geolibre-rust/releases/download/"
    f"{RUNTIME_VERSION}/{_ASSET}"
)
_STDOUT_CAPTURE = ".geolibre-stdout"
#: Network timeout (seconds) for fetching http(s) input files.
_INPUT_TIMEOUT = 120

PathLike = Union[str, os.PathLike]
#: An input file's contents: raw bytes, an http(s) URL to fetch, or a path to a
#: local file to read. Works the same for raster and vector inputs.
InputSource = Union[bytes, bytearray, memoryview, str, os.PathLike]


def _materialize_input(value: InputSource) -> bytes:
    """Resolves an input value to the bytes written under ``/work``.

    Args:
        value: Raw bytes, an ``http(s)`` URL to download, or a path to a local
            file to read.

    Returns:
        The file contents as ``bytes``.
    """
    if isinstance(value, (bytes, bytearray, memoryview)):
        return bytes(value)
    if isinstance(value, str) and value.startswith(("http://", "https://")):
        return _http_get(value, _INPUT_TIMEOUT)
    if isinstance(value, (str, os.PathLike)):
        path = Path(value)
        if path.is_file():
            return path.read_bytes()
    raise TypeError(
        "input values must be bytes, an http(s) URL, or a path to an existing "
        f"file; got {value!r}"
    )


@dataclass
class ToolResult:
    """Result of running a tool.

    Attributes:
        exit_code: Process exit code (0 = success).
        stdout: Captured stdout/stderr lines.
        files: New files the tool wrote, keyed by path relative to ``/work``.
            Tools that write a tree (e.g. ``raster_to_tiles``) use nested keys
            such as ``"tiles/15/4779/16383.png"``.
    """

    exit_code: int
    stdout: list[str]
    files: dict[str, bytes]


_engine: Optional[wasmtime.Engine] = None
_module_cache: dict[str, "wasmtime.Module"] = {}


def _get_engine() -> "wasmtime.Engine":
    global _engine
    if _engine is None:
        _engine = wasmtime.Engine()
    return _engine


def _cache_path() -> Path:
    base = os.environ.get("XDG_CACHE_HOME") or (Path.home() / ".cache")
    return Path(base) / "geolibre" / _ASSET


def download_runtime(dest: Optional[PathLike] = None) -> str:
    """Download the ``geolibre-cli.wasm`` runtime from the GitHub release.

    Args:
        dest: Where to write the file. Defaults to the per-user cache
            (``$XDG_CACHE_HOME/geolibre/`` or ``~/.cache/geolibre/``).

    Returns:
        The path to the downloaded runtime.
    """
    target = Path(dest) if dest is not None else _cache_path()
    target.parent.mkdir(parents=True, exist_ok=True)
    data = _http_get(_RELEASE_URL, _DOWNLOAD_TIMEOUT)
    # Guard against truncated or error-page downloads: every wasm module starts
    # with the "\0asm" magic.
    if not data.startswith(_WASM_MAGIC):
        raise RuntimeError(
            f"downloaded runtime from {_RELEASE_URL} is not a valid WASM module"
        )
    tmp = target.with_name(target.name + ".download")
    tmp.write_bytes(data)
    tmp.replace(target)
    return str(target)


def runtime_path(wasm_path: Optional[PathLike] = None) -> str:
    """Resolve the runtime ``.wasm`` to use.

    Resolution order: the explicit ``wasm_path`` argument, the ``GEOLIBRE_WASM``
    environment variable, then the cached download (fetched on first use).

    Args:
        wasm_path: Explicit path to a ``geolibre-cli.wasm``; takes precedence.

    Returns:
        Filesystem path to the runtime module.
    """
    if wasm_path is not None:
        return str(wasm_path)
    env = os.environ.get("GEOLIBRE_WASM")
    if env:
        return env
    cache = _cache_path()
    if not cache.exists():
        download_runtime(cache)
    return str(cache)


def _load_module(path: str) -> "wasmtime.Module":
    module = _module_cache.get(path)
    if module is None:
        module = wasmtime.Module.from_file(_get_engine(), path)
        _module_cache[path] = module
    return module


def _exec(
    argv: Sequence[str],
    inputs: Optional[Mapping[str, InputSource]] = None,
    wasm_path: Optional[PathLike] = None,
) -> ToolResult:
    inputs = inputs or {}
    module = _load_module(runtime_path(wasm_path))
    engine = _get_engine()
    work = Path(tempfile.mkdtemp(prefix="geolibre-"))
    try:
        for name, value in inputs.items():
            dest = work / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(_materialize_input(value))

        stdout_file = work / _STDOUT_CAPTURE
        store = wasmtime.Store(engine)
        wasi = wasmtime.WasiConfig()
        wasi.argv = ["geolibre", *argv]
        wasi.preopen_dir(str(work), "/work")
        # Merge stdout + stderr into one capture, matching tools.mjs.
        wasi.stdout_file = str(stdout_file)
        wasi.stderr_file = str(stdout_file)
        store.set_wasi(wasi)

        linker = wasmtime.Linker(engine)
        linker.define_wasi()
        instance = linker.instantiate(store, module)
        start = instance.exports(store)["_start"]

        exit_code = 0
        try:
            start(store)  # WASI returns by calling proc_exit -> ExitTrap
        except wasmtime.ExitTrap as exit_trap:
            exit_code = exit_trap.code

        stdout: list[str] = []
        if stdout_file.exists():
            stdout = stdout_file.read_text(errors="replace").splitlines()

        input_names = set(inputs)
        files: dict[str, bytes] = {}
        for path in sorted(work.rglob("*")):
            if not path.is_file():
                continue
            rel = path.relative_to(work).as_posix()
            if rel == _STDOUT_CAPTURE or rel in input_names:
                continue
            files[rel] = path.read_bytes()

        return ToolResult(exit_code=exit_code, stdout=stdout, files=files)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def list_tools(wasm_path: Optional[PathLike] = None) -> list[str]:
    """List every available tool id.

    Args:
        wasm_path: Optional explicit runtime path (see :func:`runtime_path`).

    Returns:
        The tool ids, one per registered tool.
    """
    result = _exec(["list"], wasm_path=wasm_path)
    return [line.strip() for line in result.stdout if line.strip()]


def list_manifests(wasm_path: Optional[PathLike] = None) -> list[dict]:
    """Fetch every tool manifest (id, parameters, category, provenance, ...).

    Args:
        wasm_path: Optional explicit runtime path (see :func:`runtime_path`).

    Returns:
        A list of manifest dicts, each including a ``"source"`` field of
        ``"geolibre"`` or ``"whitebox"``.
    """
    result = _exec(["manifests"], wasm_path=wasm_path)
    return json.loads("".join(result.stdout))


def run_tool(
    tool: str,
    args: Optional[Sequence[str]] = None,
    input: Optional[Mapping[str, InputSource]] = None,
    wasm_path: Optional[PathLike] = None,
) -> ToolResult:
    """Run one tool over an in-memory ``/work`` filesystem.

    Args:
        tool: Tool id, e.g. ``"slope"`` (see :func:`list_tools`).
        args: CLI args, e.g.
            ``["--input=/work/dem.tif", "--output=/work/slope.tif", "--units=degrees"]``.
        input: Files placed under ``/work`` before the run, keyed by filename.
            Each value may be ``bytes``, an ``http(s)`` URL to download, or a path
            to a local file to read. This is format-agnostic: it works for raster
            inputs (``cog.tif``) and vector inputs (``data.parquet``,
            ``data.geojson``, ``data.fgb``, ...) alike. The whole file is fetched
            (no HTTP range reads).
        wasm_path: Optional explicit runtime path (see :func:`runtime_path`).

    Returns:
        A :class:`ToolResult` with the exit code, captured output, and any new
        files the tool wrote (keyed by path relative to ``/work``).
    """
    return _exec([tool, *(args or [])], inputs=input, wasm_path=wasm_path)
