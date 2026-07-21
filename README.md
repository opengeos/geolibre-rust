# geolibre-rust

[![npm version](https://img.shields.io/npm/v/geolibre-wasm.svg)](https://www.npmjs.com/package/geolibre-wasm)
[![PyPI version](https://img.shields.io/pypi/v/geolibre-wasm.svg)](https://pypi.org/project/geolibre-wasm/)
[![npm downloads](https://img.shields.io/npm/dm/geolibre-wasm.svg)](https://www.npmjs.com/package/geolibre-wasm)
[![CI](https://github.com/opengeos/geolibre-rust/actions/workflows/ci.yml/badge.svg)](https://github.com/opengeos/geolibre-rust/actions/workflows/ci.yml)
[![license](https://img.shields.io/npm/l/geolibre-wasm.svg)](https://github.com/opengeos/geolibre-rust#license)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/opengeos/geolibre-rust/blob/main/examples/geolibre_wasm.ipynb)

A pure-Rust geospatial toolkit for [GeoLibre](https://github.com/opengeos/GeoLibre),
built on [`opengeos/whitebox-wasm`](https://github.com/opengeos/whitebox-wasm)
(the WASM-ready fork of
[`whitebox_next_gen`](https://github.com/jblindsay/whitebox_next_gen)) and
compiled to WebAssembly. It is a **superset of `whitebox-wasm`**: everything that
package offers, plus new GeoLibre-authored tools.

The published npm package (`geolibre-wasm`) ships two layers:

- **Browser library** (`.` export, `wasm-bindgen`) -- typed in-memory APIs for
  GeoTIFF/COG read+write, projections, vector, LiDAR, and topology
  (`GeoTiffReader`, `CogBuilder`, `CogStream`, ...). Same surface as
  `whitebox-wasm`.
- **Tool runner** (`./tools` export, WASI) -- the whitebox tool registry **plus
  GeoLibre's own tools**, run over an in-memory `/work` filesystem via
  [`@bjorn3/browser_wasi_shim`](https://github.com/bjorn3/browser_wasi_shim).

No server, no GDAL, no native install. Use it from JavaScript (npm
`geolibre-wasm`) or Python (PyPI `geolibre-wasm`). New tools live in the
`geolibre-tools` crate and are registered alongside whitebox's, so GeoLibre sees
them through the same interface as the built-ins.

## Try it in the browser

`demo/index.html` is a self-contained page that loads every tool manifest,
renders a parameter form for whichever tool you pick, and runs it on a sample DEM
(or your own GeoTIFF) entirely in the browser via the WASI runner.

```bash
./build.sh          # once, to produce npm/geolibre-cli.wasm and npm/tools.mjs
./demo/serve.sh     # serve on http://localhost:8000 (pass a port to override)
```

Open the printed URL, filter the tool list, fill in the auto-generated form, and
click **Run** to see the exit code, stdout, output files, and a download link.
`serve.sh` stages the runtime (`npm/tools.mjs`, `npm/geolibre-cli.wasm`) and the
sample raster (`examples/sample.tif`) next to the page in a temp directory, so the
repo's `demo/` stays clean; Ctrl-C stops the server and cleans up.

### Self-host with Docker

The same demo ships as a container image so you can host it yourself. The image
is a static site (nginx) — every tool still runs in the visitor's browser, so
there's no server-side compute, GDAL, or database.

Pull the published image with Docker Compose:

```bash
docker compose up -d        # serves on http://localhost:8080
```

Or build and run it straight from source (needs only Docker — the Rust/WASM
toolchain lives inside the build):

```bash
docker build -t geolibre-rust .
docker run --rm -p 8080:80 geolibre-rust
```

Images are published to `ghcr.io/opengeos/geolibre-rust` on each release (and on
`v*` tags); `:latest` tracks the most recent release. To build from source via
Compose instead of pulling, uncomment `build: .` in `docker-compose.yml`.

## Architecture

```
crates/geolibre-wasm   wasm-bindgen browser library  -> geolibre_wasm{.js,_bg.wasm,.d.ts}  (npm ".")
crates/geolibre-cli    WASI tool runner              -> geolibre-cli.wasm + tools.mjs       (npm "./tools")
crates/geolibre-tools  new Tool impls (raster_normalize, ...), registered by geolibre-cli

JS (browser/Node)                WASI binary (geolibre-cli.wasm)
-----------------                --------------------------------
tools.mjs                        crates/geolibre-cli (main.rs)
  write inputs -> /work    -->     argv -> ToolArgs (JSON)
  argv ["slope", "--..."]  -->     ToolRegistry::run
  read new files from /work <--      register_default_tools (whitebox)
                                     + geolibre_tools (new tools)
                                   tool writes via std::fs to /work
```

## GeoLibre-authored tools

In addition to the whitebox suite, `geolibre-tools` ships cloud-native I/O and
rendering tools that the whitebox suite lacks (all pure-Rust, running in WASM):

| Tool id | What it does |
|---|---|
| `reproject_raster` | Reproject (warp) a raster into a target EPSG CRS, with selectable resampling. |
| `assign_projection_raster` | Assign an EPSG CRS to a raster's metadata without warping its cells (for data whose coordinates are already in that CRS but carry a missing/wrong projection tag). |
| `assign_projection_vector` | Assign an EPSG CRS to a vector layer without reprojecting its geometries. |
| `assign_projection_lidar` | Assign an EPSG CRS to a LiDAR point cloud (LAS/LAZ/COPC) without reprojecting its points. |
| `render_raster_png` | Render a raster band to a PNG through a colormap (viridis/magma/turbo/terrain/grayscale); no-data becomes transparent. |
| `raster_to_tiles` | Slice a raster into a Web Mercator (EPSG:3857) XYZ PNG tile pyramid for web maps. |
| `write_pmtiles` | Render a raster into a single PMTiles archive (the Web Mercator PNG pyramid as one file). |
| `vector_to_pmtiles` | Pack a vector layer (GeoJSON, Shapefile, GeoPackage, FlatGeobuf, GeoParquet, ...) into a single PMTiles archive of Mapbox Vector Tiles, ready to style in MapLibre. Clipping, per-zoom simplification and MVT encoding come from [freestiler](https://walker-data.com/freestiler/). |
| `pmtiles_extract` | Extract a bbox/zoom subset of a PMTiles archive into a new self-contained archive (e.g. an offline basemap from a Protomaps planet build). The browser library exposes the same engine as `PmtilesExtractor`, driven by host `fetch` range requests. |
| `spectral_index` | Compute a spectral index (NDVI, NDWI, NDBI, NBR, EVI, SAVI) from a multi-band raster. |
| `regularize_building_footprints` | Normalize noisy building footprint polygons into regular shapes (like ArcGIS's *Regularize Building Footprint*): snap walls to right angles, right angles + 45° diagonals, straighten at any angle, or fit a best-fit circle; features that can't be regularized within the tolerance keep their original shape and are flagged in a `status` field. |
| `smooth_natural_features` | Smooth pixelated polygons and jagged lines (raster-to-vector outputs such as land cover, water bodies, or vegetation masks) into natural-looking curves — like [Smoothify](https://github.com/DPIRD-DMA/Smoothify): Douglas–Peucker de-noising plus Chaikin corner cutting, with each polygon's original area restored afterwards. |
| `eliminate_polygons` | Merge sliver polygons into a neighbor (like ArcGIS's *Eliminate*): select slivers by maximum area and/or a simple attribute query (with an optional `exclude` filter to protect features), then dissolve each into the neighbor sharing the longest border or having the largest area. The classic cleanup after an overlay or raster-to-vector conversion; slivers with no polygon neighbor are kept and reported. |
| `simplify_shared_edges` | Simplify a polygon coverage while keeping boundaries shared between adjacent polygons coincident (like ArcGIS's *Simplify Shared Edges*): build an arc–node topology, Douglas–Peucker each shared arc exactly once, and reassemble every polygon — so no gaps or slivers open up the way per-feature `simplify_features` would. Optionally leave the outer boundary untouched, and snap nearly-coincident vertices first. |
| `smooth_shared_edges` | Smooth a polygon coverage while keeping boundaries shared between adjacent polygons coincident (like ArcGIS's *Smooth Shared Edges*): the smoothing twin of `simplify_shared_edges` — build the same arc–node topology, smooth each shared arc exactly once (`paek` Gaussian-kernel smoothing or `bezier` Chaikin corner-cutting) with junction nodes pinned, and reassemble — so curving a coverage never opens gaps or slivers the way per-feature `smooth_natural_features` would. Optionally leave the outer boundary untouched, and snap nearly-coincident vertices first. |
| `cartogram` | Distort polygons so area is proportional to an attribute value (like ArcGIS's *Cartogram* toolset): `non_contiguous` (scale each polygon about its centroid, preserving shape and conserving total area) or `dorling` (proportional circles placed at centroids with force-directed overlap removal). Output feeds straight into `render_vector_png` / `vector_to_pmtiles`. |
| `build_balanced_zones` | Group contiguous polygons into balanced zones (like ArcGIS's *Build Balanced Zones* / *Spatially Constrained Multivariate Clustering*): a SKATER spanning-tree partition (contiguity graph → minimum spanning forest → greedy edge cutting) that keeps every zone connected while balancing feature count, an attribute sum, or attribute homogeneity. Rook/queen contiguity; deterministic, no RNG. For districting, sales territories, and ecological regionalization. |
| `geographically_weighted_regression` | Local linear regression with distance-decay kernel weights (like ArcGIS's *Geographically Weighted Regression*): fit a separate weighted least-squares model at every feature (gaussian or bisquare kernel, fixed or adaptive bandwidth, optionally AICc-optimized), producing per-feature coefficients, local R², residuals and predictions, plus global AICc / R² diagnostics. |
| `ripleys_k` | Multi-distance point-pattern analysis (like ArcGIS's *Multi-Distance Spatial Cluster Analysis / Ripley's K*): the K/L function across distance bands with Monte-Carlo complete-spatial-randomness envelopes (deterministic seeded RNG), to detect clustering (L above the envelope) or dispersion (below) across scales. Outputs a distance/observed/expected/envelope table. |
| `emerging_hot_spot_analysis` | Space-time hot spot trends on an H3 space-time cube (like ArcGIS's *Emerging Hot Spot Analysis* + *Create Space Time Cube*): bin timestamped points into H3 cells × time steps, compute the Getis-Ord Gi\* z-score per bin over a space-time neighborhood (spatial k-ring × ± time window), run a Mann-Kendall trend test on each cell's Gi\* series, and classify every cell into the Esri categories (new / consecutive / intensifying / persistent / diminishing / sporadic / oscillating / historical hot or cold spot, or no pattern). Adds the temporal dimension the bundled spatial-only `getis_ord_gi_star` lacks; deterministic, no RNG; output H3 polygons render straight through `h3_to_vector` / PMTiles. |
| `cut_fill` | Volumetric change between two DEM surfaces (like ArcGIS's *Cut Fill* / *Surface Volume*): a signed elevation-change raster (Δz = after − before, or surface − a reference plane) plus cut, fill, and net volumes (Σ \|Δz\| × cell area), with optional contiguous cut/fill region labelling and a per-region volume CSV. For earthworks, erosion/deposition, stockpiles, and lidar change detection. |
| `line_of_sight` | Point-to-point visibility over a DEM (like ArcGIS's *Line Of Sight* / *Construct Sight Lines*): for each observer→target pair, walk the DEM along the segment tracking the running maximum vertical angle, and emit the sight line split into contiguous **visible** and **obstructed** LineString segments plus a target-visible flag. Answers "can A see B, and where does the view break" — the per-pair vector question the raster `viewshed` / `skyline_analysis` can't. Observer/target heights, all-to-all or `pair_field` matching. |
| `corridor` | Least-cost corridor between two accumulated-cost surfaces (like ArcGIS's *Corridor* / Least Cost Corridor): sum two `cost_distance` rasters into a corridor surface whose every cell is the cost of the cheapest A→cell→B path through it, then optionally threshold (absolute or percent-above-minimum) to the swath of near-optimal routes the bundled cost suite (`cost_distance`/`cost_pathway`, which give a single path) can't. Direct (two cost rasters) or convenience mode (one friction raster + two source rasters, accumulated internally). For wildlife corridors and route planning. |
| `interpolate_shape` | Drape points/lines/polygons on a surface raster (like ArcGIS's *Interpolate Shape* + *Add Surface Information*): densify each geometry at a sample interval, sample Z per vertex (bilinear or nearest), write true 3D geometry, and add surface metrics — `z_min`/`z_max`/`z_mean`, 3D surface length, and length-weighted average slope. The bridge between the repo's raster (terrain) and vector halves — trail profiles, pipeline lengths, slope-aware routing — that no bundled point-sampling tool provides for lines/polygons. |
| `thin_road_network` | Generalize a road network for small-scale display (like ArcGIS's *Thin Road Network*): hide short, low-hierarchy roads while preserving connectivity — a candidate road is thinned only if it is not a bridge in the current network, so the visible network keeps its connected components. Non-destructive (flags each road visible/thinned) or filtered. Pairs with `download_osm_vector`. |
| `tabulate_intersection` | Vector-on-vector zonal summary (like ArcGIS's *Tabulate Intersection* and the core of *Summarize Within*): apportion a class layer across zone polygons — polygons summarized by intersected area, points by count — reporting the measure, its percentage of each zone, and area-weighted `sum_fields`. Answers "how much of each class, or how many points, fall inside each zone". |
| `directional_distribution` | Descriptive geographic-distribution statistics (like ArcGIS's *Measuring Geographic Distributions*): mean center, median center (Weiszfeld), central feature, standard distance (circle), and the standard deviational ellipse (from the coordinate covariance eigenvectors) — with optional weighting and per-group `case_field` output. |
| `multiple_ring_buffer` | Buffer features at a list of distances into concentric bands (like ArcGIS's *Multiple Ring Buffer*): non-overlapping rings (donuts) or nested disks, optionally dissolved across features per distance, with the band distance stored as an attribute. The everyday proximity-zone tool (catchments, impact bands, setbacks). |
| `aggregate_polygons` | Combine polygons within an aggregation distance into larger polygons (like ArcGIS's *Aggregate Polygons*): a morphological closing (dilate–erode by half the `aggregation_distance`) fuses nearby polygons, with minimum-area and minimum-hole-size filters and an optional `barrier` layer (lines or polygons) that aggregation may not cross. Pairs with `regularize_building_footprints`; the general-purpose sibling of `delineate_built_up_areas`. Each output carries a `part_count` of the source polygons it merges. |
| `delineate_built_up_areas` | Generalize dense building footprints into smooth settlement / urban-extent polygons (like ArcGIS's *Delineate Built-Up Areas*): a morphological closing (dilate–erode by half the `grouping_distance`) fuses footprints within the grouping distance into one built-up area, with rounded, cartographic boundaries. Filter the result by a minimum building count and minimum area, fill small interior holes, and optionally simplify the boundary. The natural sequel to AI building-footprint extraction. |
| `write_geoparquet` | Convert any supported vector format to GeoParquet, Hilbert-sorted with a bbox covering column and ZSTD compression by default. |
| `read_geoparquet` | Read GeoParquet and convert it to another vector format (or store it in memory). |
| `vector_convert` | Convert a vector dataset between formats (the output extension picks the driver). |
| `render_vector_png` | Draw a vector layer (points/lines/polygons) to a PNG map image. |

It also ships pure-Rust ports of the DEM depression/mount algorithms from
[`opengeos/lidar`](https://github.com/opengeos/lidar) (no GDAL, RichDEM, SciPy,
or scikit-image dependency; they run in WASM):

| Tool id | Source | What it does |
|---|---|---|
| `dem_filter` | `filtering.py` | Mean / median / Gaussian smoothing of a DEM. |
| `extract_sinks` | `filling.py` | Wang & Liu fill, then group filled cells into sinks larger than `min_size`; emits sink/region/depth/filled rasters, an attribute CSV, and region polygons (`vector_output`, GeoJSON). |
| `delineate_depressions` | `slicing.py` | Level-set slicing of a sink raster into a nested-depression hierarchy; emits id/level rasters, a CSV, and depression polygons (`vector_output`, GeoJSON). |
| `delineate_mounts` | `mounts.py` | Flip the DEM, then run the sink + depression pipeline to delineate nested elevated features (rasters, CSV, and GeoJSON). |

Typical chain (over the WASI `/work` filesystem or via paths):

```text
extract_sinks --input=dem.tif --output=sink.tif --min_size=100
delineate_depressions --input=sink.tif --output=dep_id.tif --level_output=dep_level.tif
```

Depression filling reuses a port of whitebox's Wang & Liu priority-flood (kept
inside `geolibre-tools` so the crate stays free of a `wbtools_oss` dependency).
The morphological attributes (perimeter, axes, eccentricity, orientation) mirror
`scikit-image`'s `regionprops`. The `vector_output` parameter polygonizes the
label raster into GeoJSON (one feature per connected component, holes preserved,
RFC 7946 winding) in the source CRS, with the attribute table joined onto each
feature -- a pure-Rust replacement for the `gdal.Polygonize` + GeoPackage join
in the Python original.

## Adding a new tool

1. Add a module with a `wbcore::Tool` impl under `crates/geolibre-tools/src/`
   (see `raster_normalize.rs` for the template: `metadata` / `validate` / `run`,
   reading and writing rasters by path).
2. Push it into the list returned by `geolibre_tools()` in
   `crates/geolibre-tools/src/lib.rs`.
3. **Declare each parameter's schema** in `geolibre_param_schemas()` (same file):
   add a match arm mapping the tool id to its params, e.g.
   `input` -> `ToolParamSchema::input_raster()`, `output` ->
   `ToolParamSchema::output_raster()` (or `output(ToolDatasetSchema::File)` for a
   non-dataset file), a count -> `scalar_integer()`, a factor -> `scalar_float()`,
   a flag -> `bool()`, a fixed choice -> `enum_values(&[...])`. This is what makes
   `geolibre manifests` emit an accurate `io_role`/`data_kind`/`schema` per param,
   so host UIs route raster/vector inputs and render the right widget. Without it,
   the manifest falls back to keyword inference, which mis-types scalars whose
   description mentions a dataset (a flag that "sorts features" would read as a
   vector input). A unit test (`every_tool_has_explicit_param_schemas`) fails if a
   tool's param is left without a schema.
4. Rebuild (`./build.sh`); it appears in `listTools()` automatically.

The crate depends only on `wbcore` plus the data crates a tool needs (e.g.
`wbraster`), so the same tools can later back a native CLI or the Python sidecar,
not just WASM.

The data boundary is the WASI virtual filesystem: inputs are placed under
`/work`, tools read/write there via ordinary `std::fs`, and any new file is
returned to JS as a `Uint8Array`. Raster outputs are Cloud Optimized GeoTIFFs.

## CLI contract

```text
geolibre list                 # print every tool id
geolibre manifests            # print all tool manifests as JSON (param schemas)
geolibre manifest <id>        # print one manifest as JSON
geolibre version              # print the runner version
geolibre <tool> [--k=v ...]   # run a tool over /work
```

`--key=value`, `--key value`, and bare `--flag` are all accepted. Values are
type-inferred: `true`/`false` -> bool, numbers -> number, everything else
(including `/work/...` paths) -> string.

## Build

```bash
rustup target add wasm32-wasip1
sudo apt-get install -y binaryen   # provides wasm-opt
./build.sh                         # -> npm/geolibre-cli.wasm
```

The `whitebox_next_gen` crates are referenced as path dependencies in
`crates/geolibre-cli/Cargo.toml`. Switch them to git or published versions before
releasing (note `wbtools_oss` is `publish = false`, so git or vendoring is
required).

### TODO: remove the vendored `kdtree` patch once `kdtree 0.8.1` ships

`vendor/kdtree/` and the `[patch.crates-io]` block in the root `Cargo.toml` work
around a bug in the published `kdtree 0.8.0` (it declares `criterion` as a normal
dependency, which pulls `rayon` and breaks the WASI build). The fix is already on
`kdtree-rs` `master` (PRs #70 and #89) but unreleased. Tracking issue:
https://github.com/mrhooray/kdtree-rs/issues/91

When `kdtree 0.8.1` (or later) is published, delete `vendor/kdtree/` and the
`[patch.crates-io]` block, then rebuild to confirm the WASI build stays green.

## Use from JavaScript

> Note: the repository is `geolibre-rust` (the Rust source), but the published
> npm package is **`geolibre-wasm`** (the WASM artifact), mirroring `whitebox-wasm`.

```bash
npm install geolibre-wasm
```

Browser library (the `.` export) -- typed GeoTIFF/projection/vector/LiDAR APIs:

```js
import init, { GeoTiffReader, CogBuilder, version } from "geolibre-wasm";

await init(); // load the wasm-bindgen module
const r = new GeoTiffReader(tiffBytes);   // Uint8Array
console.log(r.width, r.height, r.bands, r.epsg);
const band0 = r.read_band_f64(0);          // Float64Array
```

Tool runner (the `./tools` export) -- the whitebox + GeoLibre tool suite:

```js
import { runTool, listTools } from "geolibre-wasm/tools";

const tools = await listTools();

const { files } = await runTool("slope", {
  args: ["--input=/work/dem.tif", "--output=/work/slope.tif", "--units=degrees"],
  input: { "dem.tif": demBytes }, // Uint8Array
});
const slopeCog = files["slope.tif"]; // Uint8Array (COG GeoTIFF)
```

An `input` value may also be an `http(s)` URL string, fetched for you (whole
file, no range reads) -- the same for raster and vector inputs:

```js
await runTool("write_geoparquet", {
  args: ["--input=/work/in.geojson", "--output=/work/out.parquet"],
  input: { "in.geojson": "https://example.com/data/cities.geojson" },
});
```

## Use from Python

The `python/` package (`geolibre-wasm` on PyPI, `import geolibre_wasm`) runs the
same WASI tool runner in-process via `wasmtime`, mirroring the JS `./tools` API.
No native install, GDAL, or server.

Try it in your browser, no setup:
[**Open the quickstart notebook in Google Colab**](https://colab.research.google.com/github/opengeos/geolibre-rust/blob/main/examples/geolibre_wasm.ipynb)
([`examples/geolibre_wasm.ipynb`](examples/geolibre_wasm.ipynb)) -- it reads and
processes a real DEM and building footprints end to end.

```bash
pip install geolibre-wasm
```

```python
import geolibre_wasm as gl

tools = gl.list_tools()                 # every tool id
manifests = gl.list_manifests()         # schemas + "source": geolibre|whitebox

res = gl.run_tool(
    "slope",
    # Paths in `args` refer to the tool's sandbox (/work), NOT your host disk.
    # `input` files are placed at /work/<name>; `res.files` keys are relative
    # to /work. Mixing in host paths (e.g. /content on Colab) will not work.
    args=["--input=/work/dem.tif", "--output=/work/slope.tif", "--units=degrees"],
    input={"dem.tif": open("dem.tif", "rb").read()},   # -> /work/dem.tif
)
assert res.exit_code == 0, res.stdout                  # surfaces tool errors
open("slope.tif", "wb").write(res.files["slope.tif"])  # key is relative to /work
```

Each `input` value may be `bytes`, an `http(s)` URL (downloaded for you), or a
local file path -- the same for raster and vector inputs:

```python
gl.run_tool(
    "write_geoparquet",
    args=["--input=/work/in.geojson", "--output=/work/out.parquet"],
    input={"in.geojson": "https://example.com/data/cities.geojson"},
)
```

The runtime `.wasm` is downloaded from the matching release on first use (or set
`GEOLIBRE_WASM`). See [`python/README.md`](python/README.md) for details.

## Recipes: reading and processing various formats

The examples below use the Python API; the JavaScript `runTool` takes the same
`args` and `input` (just `camelCase`). Each `input` value can be `bytes`, an
`http(s)` URL, or a local path. Output files come back in `result.files` keyed by
their `/work`-relative path. These all run against the real tool suite.

### Vector (GeoParquet, GeoJSON, FlatGeobuf, Shapefile, GeoPackage, ...)

```python
import geolibre_wasm as gl

# Convert GeoJSON -> GeoParquet (Hilbert-sorted, bbox covering, ZSTD by default)
gj = open("cities.geojson", "rb").read()
gl.run_tool("write_geoparquet",
            args=["--input=/work/in.geojson", "--output=/work/out.parquet"],
            input={"in.geojson": gj})

# Read GeoParquet -> any vector format (driver picked from the output extension:
# .geojson, .fgb, .shp, .gpkg, ...). Omit --output to keep it in memory.
res = gl.run_tool("read_geoparquet",
                  args=["--input=/work/in.parquet", "--output=/work/out.fgb"],
                  input={"in.parquet": "https://example.com/data.parquet"})
open("out.fgb", "wb").write(res.files["out.fgb"])

# Buffer features, then simplify (Douglas-Peucker)
res = gl.run_tool("buffer_vector",
                  args=["--input=/work/in.geojson", "--distance=25", "--output=/work/buf.geojson"],
                  input={"in.geojson": gj})
res = gl.run_tool("simplify_features",
                  args=["--input=/work/buf.geojson", "--tolerance=5", "--output=/work/simple.geojson"],
                  input={"buf.geojson": res.files["buf.geojson"]})

# Add geometry attributes (area / perimeter / centroid, ...)
gl.run_tool("add_geometry_attributes",
            args=["--input=/work/in.geojson", "--area=true", "--centroid=true",
                  "--output=/work/attrs.geojson"],
            input={"in.geojson": gj})
```

`reproject_vector` works the same, but the input must carry a source CRS (a
Shapefile `.prj`, or a GeoParquet/GeoPackage with CRS metadata), e.g.
`args=["--input=/work/in.fgb", "--epsg=3857", "--output=/work/out.fgb"]`.

### LiDAR point clouds (LAS / LAZ)

```python
import geolibre_wasm as gl

cloud = open("cloud.las", "rb").read()        # or a .laz, or an http(s) URL

# Summary report (point count, bounds, density, ...). Output must be .txt/.html.
res = gl.run_tool("lidar_info",
                  args=["--input=/work/cloud.las", "--output=/work/info.txt"],
                  input={"cloud.las": cloud})
print(res.files["info.txt"].decode())

# Rasterize to a DEM (IDW) -> Cloud Optimized GeoTIFF
res = gl.run_tool("lidar_idw_interpolation",
                  args=["--input=/work/cloud.las", "--resolution=1.0", "--output=/work/dtm.tif"],
                  input={"cloud.las": cloud})
open("dtm.tif", "wb").write(res.files["dtm.tif"])

# Drop unwanted classes (comma-delimited list; e.g. exclude 1=unclassified, 7=noise)
gl.run_tool("filter_lidar_classes",
            args=["--input=/work/cloud.las", "--excluded_classes=1,7", "--output=/work/clean.las"],
            input={"cloud.las": cloud})

# Export points to a Shapefile
gl.run_tool("las_to_shapefile",
            args=["--input=/work/cloud.las", "--output=/work/points.shp"],
            input={"cloud.las": cloud})
```

### Raster (GeoTIFF / COG)

```python
dem = open("dem.tif", "rb").read()            # or an http(s) URL to a COG

# Warp to Web Mercator, then render a PNG preview through a colormap
gl.run_tool("reproject_raster",
            args=["--input=/work/dem.tif", "--epsg=3857", "--output=/work/merc.tif"],
            input={"dem.tif": dem})
gl.run_tool("render_raster_png",
            args=["--input=/work/dem.tif", "--colormap=terrain", "--output=/work/preview.png"],
            input={"dem.tif": dem})
```

Run `gl.list_tools()` for all 740+ tool ids and `gl.list_manifests()` for each
tool's parameters and provenance (`"source": "geolibre" | "whitebox"`).

## GeoLibre integration

The interface is byte-compatible with the existing `whitebox-wasm/tools` client:

1. Add `geolibre-wasm` to `packages/processing/package.json`.
2. Add it to `optimizeDeps.exclude` in `apps/geolibre-desktop/vite.config.ts`
   (required for the `new URL("./*.wasm", import.meta.url)` glue).
3. Point `packages/processing/src/wasm-client.ts`'s lazy
   `import("whitebox-wasm/tools")` at `geolibre-wasm/tools`, or add a sibling
   client and a source toggle in `ProcessingDialog.tsx`.

`listManifests()` is a value-add over the legacy package: it lets GeoLibre build
tool dialogs (parameter schemas, `raster_in`/`vector_in` roles) fully offline,
without the Python sidecar.

## License

MIT. See [LICENSE](LICENSE).
