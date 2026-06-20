"""Smoke tests for the geolibre Python wrapper.

These run the real WASI tool runner. Point ``GEOLIBRE_WASM`` at a local
``geolibre-cli.wasm`` to avoid a network download (the build stages one into
``npm/``); otherwise the runtime is fetched from the GitHub release on first use.
"""

from __future__ import annotations

import json

import pytest

import geolibre_wasm as geolibre
from geolibre_wasm._core import _materialize_input


def test_list_tools():
    tools = geolibre.list_tools()
    assert len(tools) > 100
    # The GeoLibre-authored tools are registered alongside whitebox's.
    for tool in ("slope", "reproject_raster", "write_geoparquet", "read_geoparquet"):
        assert tool in tools


def test_manifests_carry_provenance():
    manifests = geolibre.list_manifests()
    by_id = {m["id"]: m for m in manifests}
    assert by_id["reproject_raster"]["source"] == "geolibre"
    assert by_id["slope"]["source"] == "whitebox"


def test_geoparquet_roundtrip():
    geojson = json.dumps(
        {
            "type": "FeatureCollection",
            "features": [
                {
                    "type": "Feature",
                    "properties": {"name": "A", "val": 1},
                    "geometry": {"type": "Point", "coordinates": [-122.3, 47.6]},
                },
                {
                    "type": "Feature",
                    "properties": {"name": "B", "val": 2},
                    "geometry": {"type": "Point", "coordinates": [-122.4, 47.7]},
                },
            ],
        }
    ).encode()

    written = geolibre.run_tool(
        "write_geoparquet",
        args=["--input=/work/in.geojson", "--output=/work/out.parquet"],
        input={"in.geojson": geojson},
    )
    assert written.exit_code == 0
    parquet = written.files["out.parquet"]
    assert parquet[:4] == b"PAR1"

    back = geolibre.run_tool(
        "read_geoparquet",
        args=["--input=/work/x.parquet", "--output=/work/back.geojson"],
        input={"x.parquet": parquet},
    )
    assert back.exit_code == 0
    fc = json.loads(back.files["back.geojson"])
    assert len(fc["features"]) == 2


def test_materialize_input_sources(tmp_path):
    # bytes pass through
    assert _materialize_input(b"abc") == b"abc"
    assert _materialize_input(bytearray(b"xy")) == b"xy"
    # a local file path is read
    p = tmp_path / "f.bin"
    p.write_bytes(b"hello")
    assert _materialize_input(p) == b"hello"
    assert _materialize_input(str(p)) == b"hello"
    # a non-URL, non-file string is rejected
    with pytest.raises(TypeError):
        _materialize_input("not-a-url-or-file")


def test_run_tool_accepts_a_local_path_input(tmp_path):
    # run_tool's `input` values may be a path, not just bytes. (URL inputs use
    # the same code path; they are exercised in test_materialize_input_sources
    # without hitting the network.)
    import os

    src = os.path.join("examples", "sample.tif")  # repo-relative; present in CI
    if not os.path.isfile(src):
        pytest.skip("sample raster not available")
    sample = tmp_path / "dem.tif"
    sample.write_bytes(open(src, "rb").read())
    res = geolibre.run_tool(
        "slope",
        args=["--input=/work/dem.tif", "--output=/work/slope.tif", "--units=degrees"],
        input={"dem.tif": sample},  # a pathlib.Path, not bytes
    )
    assert res.exit_code == 0
    assert "slope.tif" in res.files
