"""Smoke tests for the geolibre Python wrapper.

These run the real WASI tool runner. Point ``GEOLIBRE_WASM`` at a local
``geolibre-cli.wasm`` to avoid a network download (the build stages one into
``npm/``); otherwise the runtime is fetched from the GitHub release on first use.
"""

from __future__ import annotations

import json

import geolibre_wasm as geolibre


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
