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
| `integrate` | Snap all vertices across all features within a cluster tolerance to a shared location so nearly-coincident shared boundaries become exactly coincident (like ArcGIS's *Integrate* + the editing *Snap* tool): grid-hashed union-find vertex clustering → move each vertex to its cluster centroid, then optionally insert T-junction vertices where one feature's vertex lands on another's segment. The topology-cleaning precondition the coverage-safe family (`*_shared_edges`, `polygon_neighbors`, `eliminate_polygons`) all assume — the bundled `snap_endnodes` only does polyline endpoints. |
| `eliminate_polygons` | Merge sliver polygons into a neighbor (like ArcGIS's *Eliminate*): select slivers by maximum area and/or a simple attribute query (with an optional `exclude` filter to protect features), then dissolve each into the neighbor sharing the longest border or having the largest area. The classic cleanup after an overlay or raster-to-vector conversion; slivers with no polygon neighbor are kept and reported. |
| `simplify_shared_edges` | Simplify a polygon coverage while keeping boundaries shared between adjacent polygons coincident (like ArcGIS's *Simplify Shared Edges*): build an arc–node topology, Douglas–Peucker each shared arc exactly once, and reassemble every polygon — so no gaps or slivers open up the way per-feature `simplify_features` would. Optionally leave the outer boundary untouched, and snap nearly-coincident vertices first. |
| `smooth_shared_edges` | Smooth a polygon coverage while keeping boundaries shared between adjacent polygons coincident (like ArcGIS's *Smooth Shared Edges*): the smoothing twin of `simplify_shared_edges` — build the same arc–node topology, smooth each shared arc exactly once (`paek` Gaussian-kernel smoothing or `bezier` Chaikin corner-cutting) with junction nodes pinned, and reassemble — so curving a coverage never opens gaps or slivers the way per-feature `smooth_natural_features` would. Optionally leave the outer boundary untouched, and snap nearly-coincident vertices first. |
| `cartogram` | Distort polygons so area is proportional to an attribute value (like ArcGIS's *Cartogram* toolset): `non_contiguous` (scale each polygon about its centroid, preserving shape and conserving total area) or `dorling` (proportional circles placed at centroids with force-directed overlap removal). Output feeds straight into `render_vector_png` / `vector_to_pmtiles`. |
| `build_balanced_zones` | Group contiguous polygons into balanced zones (like ArcGIS's *Build Balanced Zones* / *Spatially Constrained Multivariate Clustering*): a SKATER spanning-tree partition (contiguity graph → minimum spanning forest → greedy edge cutting) that keeps every zone connected while balancing feature count, an attribute sum, or attribute homogeneity. Rook/queen contiguity; deterministic, no RNG. For districting, sales territories, and ecological regionalization. |
| `similarity_search` | Rank candidate features by attribute similarity to one or more reference features (like ArcGIS's *Similarity Search*): z-standardize the chosen numeric fields over the combined distribution, build the reference profile, and score every candidate by Euclidean distance or cosine similarity — most / least / both ends, with per-field standardized differences. "Find the tracts most like this one" for site selection and market analysis; the attribute-space ranking the bundled classifiers (which predict classes) don't do. |
| `geographically_weighted_regression` | Local linear regression with distance-decay kernel weights (like ArcGIS's *Geographically Weighted Regression*): fit a separate weighted least-squares model at every feature (gaussian or bisquare kernel, fixed or adaptive bandwidth, optionally AICc-optimized), producing per-feature coefficients, local R², residuals and predictions, plus global AICc / R² diagnostics. |
| `hdbscan` | Hierarchical density-based clustering — HDBSCAN* (like the HDBSCAN option of ArcGIS's *Density-based Clustering*): core-distance density estimate → mutual-reachability minimum spanning tree → single-linkage dendrogram → condensed cluster tree → excess-of-mass selection, with per-point `cluster_id` (−1 noise), membership `probability`, and `outlier_score`. Handles variable density with no epsilon, unlike the bundled `dbscan`. Deterministic. |
| `colocation_analysis` | Local colocation quotient between two point categories (like ArcGIS's *Colocation Analysis*): for each category-A point, the kernel-weighted fraction of its k nearest neighbours that are category B over B's global share (Leslie & Kronenfeld) — `CLQ>1` = A drawn toward B, `<1` = avoids — with a seeded conditional-permutation p-value and a colocated/isolated class. The asymmetric two-population association the single-population `ripleys_k`/`nearest_neighbour_index` can't measure. Deterministic. |
| `ripleys_k` | Multi-distance point-pattern analysis (like ArcGIS's *Multi-Distance Spatial Cluster Analysis / Ripley's K*): the K/L function across distance bands with Monte-Carlo complete-spatial-randomness envelopes (deterministic seeded RNG), to detect clustering (L above the envelope) or dispersion (below) across scales. Outputs a distance/observed/expected/envelope table. |
| `incremental_spatial_autocorrelation` | Global Moran's I across a series of increasing distance bands, with the z-score curve and its first/maximum peaks (like ArcGIS's *Incremental Spatial Autocorrelation*): binary fixed-distance weights, Esri's randomization variance (S0/S1/S2 + kurtosis), per-band `distance, morans_i, expected_i, variance, z, p` table. The defensible way to pick a clustering distance band for `getis_ord_gi_star`, which the single-scale bundled `global_morans_i` can't give. |
| `central_feature` | The two *Measuring Geographic Distributions* members `directional_distribution` lacks (ArcGIS's *Central Feature* / *Linear Directional Mean*): **central_feature** returns the actual input feature (with its attributes) minimizing total (optionally weighted, euclidean/manhattan) distance to all others; **linear_directional_mean** gives the circular mean bearing (or undirected orientation), circular variance, and mean length of a set of lines as a single mean-vector line. `case_field` grouping. |
| `calculate_motion_statistics` | Annotate timestamped track points with motion statistics (like ArcGIS's *Calculate Motion Statistics*): group points by `track_field`, sort by `time_field`, and add per-point `seq`, segment and cumulative distance, `dt`, elapsed time, instantaneous `speed`, trailing-window `avg_speed`, `accel`, `bearing` (degrees from north), and an `idle` flag (mover barely moved over the look-back window). Distances are haversine metres for a geographic CRS, CRS units otherwise; original attributes are preserved. The per-point movement derivative the reconstruct/snap/trace track tools don't emit. |
| `sort_features` | Reorder a vector layer by attribute fields or along a Hilbert space-filling curve (like ArcGIS's *Sort*, including its spatial-sort methods): `method=hilbert` sorts by the Hilbert-curve distance of each feature's bounding-box centre over the dataset extent, clustering spatially-near records so `write_geoparquet` row-group/bbox statistics prune better and `vector_to_pmtiles` packs more local features per tile; `method=attribute` sorts by a `field:asc/desc` list with a spatial tiebreak. Optionally writes the curve index as a field. Reuses the same Hilbert mapping the GeoParquet writer already uses; the bundled suite only sorts LiDAR. |
| `calculate_composite_index` | Combine several numeric attributes into a single composite index (like ArcGIS's *Calculate Composite Index*): per-variable scaling (min-max / z-score / percentile) with an optional `:reverse` for variables where high = worse, a weighted combination (mean / sum / geometric mean), and output rescaling (0–1, 0–100, or z-score). Emits the `index`, its `index_rank` and `index_pctl`, and a `<field>_scaled` column per variable. The vector-side sibling of the raster `fuzzy_overlay` — the standard build for vulnerability / deprivation / SDG indices, which the bundled `weighted_overlay`/`weighted_sum` can only do on rasters. |
| `calculate_rates` | Compute smoothed rates from count and population fields (like ArcGIS's *Calculate Rates*): `crude` (`count/pop×per`), `eb_global` (Marshall global empirical Bayes — shrink each crude rate toward the global mean, more for smaller populations), or `eb_spatial` (shrink toward a local reference rate over each area's k-nearest neighbours). Emits `crude_rate`, `smooth_rate`, and a Poisson `rate_se`. The rate-stabilization step the hot-spot statistics (`getis_ord_gi_star`, `local_morans_i_lisa`) assume but the bundled suite can't produce — small-population noise no longer dominates the map. |
| `color_polygons` | Assign each polygon a small integer colour index so no two adjacent polygons share a value (like ArcGIS's *Calculate Color Theorem Field*): build shared-edge (rook) or shared-corner (queen) contiguity, then colour with the DSATUR heuristic, staying within 4–6 colours on planar maps. Writes a 1-based `color_id` field for instant choropleth-safe styling of parcels, admin units, or `build_balanced_zones` output through `render_vector_png` / PMTiles. Reuses the shared-edge adjacency machinery from `polygon_neighbors`; `snap_tolerance` matches near-coincident borders in unclean coverages. |
| `dice` | Split polygons or polylines with more than a `vertex_limit` (default 10000) into a grid of smaller pieces (like ArcGIS's *Dice*), so downstream overlay and tiling don't choke on million-vertex geometries: an adaptive quadtree recursively quarters each oversized feature's bounding box and intersects the geometry with every quadrant (polygon parts via `geo` `BooleanOps`, line parts via Liang–Barsky clipping) until each piece is under the limit. Features under the limit pass through untouched; attributes copy to every piece. The vertex-count safety valve the area-based `subdivide_polygon` and cutter-based `split_with_lines` don't provide, e.g. before `vector_to_pmtiles`. |
| `spatial_outlier_detection` | Score points by their Local Outlier Factor (LOF) — how isolated each point is relative to the local density of its k nearest neighbours (like ArcGIS's *Spatial Outlier Detection*): computes k-distance, reachability distance, local reachability density, then `LOF = mean(LRD(neighbour)/LRD(point))`, and flags the top `percent_outlier`% (or points above an explicit `threshold`). LOF ≈ 1 is an inlier, ≫ 1 an outlier. The continuous per-point outlier score DBSCAN's binary noise label and the elevation-only `lidar_remove_outliers` can't give. |
| `bivariate_spatial_association` | Lee's L global and local statistic for where two continuous variables co-vary spatially (like ArcGIS's *Bivariate Spatial Association*): row-standardised k-nearest-neighbour weights, mean-centred variables and their spatial lags give global `L = Σ lx·ly / (√Σzx²·√Σzy²)` and local `L_i = n·lx_i·ly_i / (…)`, with a High-High/High-Low/Low-High/Low-Low class per feature and a seeded permutation-test p-value. The continuous-field counterpart of the categorical `colocation_analysis`, and the bivariate complement to the bundled univariate `global_morans_i`/`local_morans_i_lisa`. |
| `generate_trend_raster` | Fit a per-pixel temporal trend across a time series of co-registered rasters (like ArcGIS's *Generate Trend Raster*): `linear` gives OLS slope, intercept, and R²; `mann_kendall` gives the non-parametric Mann-Kendall trend test with Sen's slope and a two-sided p-value (tie- and continuity-corrected normal approximation). Pixels with fewer than `min_valid` valid observations become no-data. The per-pixel *temporal* trend the bundled `trend_surface` (spatial polynomial) and `change_vector_analysis` (two-date) can't produce — the workhorse of NDVI/temperature change monitoring, pairing with `spectral_index` and `detect_image_anomalies`. |
| `warp_raster` | Georeference a raster from ground control points (like ArcGIS's *Warp*): fit an order 1/2/3 polynomial from `col,row,x,y` GCP pairs (source pixel → world coordinate) by least squares, then resample the image into a new georeferenced grid (`nearest`/`bilinear`), reporting per-GCP residuals and RMS error. The raster half of the conflation/registration story the vector `rubbersheet_features` already covers — and the only way to georeference scanned maps or drone frames in a stack with no GDAL (the bundled `thin_plate_spline` interpolates points, it doesn't warp images). |
| `weighted_voronoi` | Weighted Voronoi (dominance / market-area) allocation raster (like ArcGIS's *Generate Weighted Voronoi*): assign each cell to the site with the smallest *weighted* distance — `multiplicative` (d/w, Apollonius: larger weight → larger territory), `additive` (d-w), or `power` (d²-w²). Output is a categorical raster of 1-based site indices over the sites' padded extent; polygonise with `raster_to_vector_polygons` for vector market areas. The unequal-site version of the bundled `voronoi_diagram`, producing the curved boundaries no exact-geometry library in the stack can. |
| `pycnophylactic_interpolation` | Tobler's mass-preserving areal interpolation (like ArcGIS's *Areal Interpolation*): turn zone-aggregated counts (e.g. census population) into a smooth density raster whose per-zone cell sums still equal the input totals. Zones are rasterised and seeded uniformly, then iterated — 3×3 mean smoothing followed by a per-zone additive mass correction with non-negativity redistribution — until convergence. The smooth alternative to the uniform-density assumption of `apportion_polygon`; absent from the bundled suite (kriging interpolates point samples, not areal aggregates). |
| `cost_connectivity` | Least-cost network connecting multiple sites over a cost surface (like ArcGIS's *Cost Connectivity* / *Optimal Region Connections*): a multi-source Dijkstra allocates every cell to its nearest source (accumulated cost + back-link); the min-cost crossing between each pair of allocation regions defines their least-cost path; `connections=mst` returns the minimum spanning tree that connects all sites at minimum total cost, `all_neighbors` every adjacent-region path. Paths are polylines with from/to ids and cost. The least-cost *network* the bundled single-source `cost_distance` and GeoLibre `path_distance`/`corridor` can't build — for wildlife corridors and infrastructure planning. |
| `locate_regions` | Find the best contiguous region(s) of a target area from a suitability raster (like ArcGIS's *Locate Regions*): best-first region growing from the highest-suitability seeds, where each candidate cell's score blends suitability with a compactness penalty (`shape` 0..1 → rounder regions), grown to a target cell count; the next region seeds outside a `min_distance` buffer of the ones already chosen. Output is a raster of 1-based region ids with per-region area and mean suitability. The siting step that turns a `fuzzy_overlay`/`weighted_overlay` surface into actual regions — which `clump` (labels existing regions) and thresholding (fragmented blobs) can't. |
| `edgematch_features` | Connect line datasets across a tile/sheet boundary (like ArcGIS's *Generate Edgematch Links* + *Edgematch Features*): match dangling endpoints (line ends not shared with another feature) one-to-one within a `tolerance` by distance — optionally disambiguated by attribute agreement on `match_fields` — then reconcile each pair (`midpoint` moves both ends to their midpoint, `move_endpoint` snaps the second onto the first), with an optional `links` layer for QA. The cross-feature one-to-one matching the blind bundled `snap_endnodes` lacks; completes the conflation suite (`integrate`, `rubbersheet_features`, `detect_feature_changes`). |
| `landtrendr` | LandTrendr temporal segmentation of a yearly image series (like ArcGIS's *Analyze Changes Using LandTrendr*): per pixel, despike then fit a piecewise-linear trajectory by greedy vertex insertion (up to `max_segments`), and report the greatest disturbance — a drop (`direction=loss`) or rise (`gain`) in a vegetation index — as its **year** (primary output), **magnitude**, and **duration**. The per-pixel change-history segmentation the two-date `change_vector_analysis` and single-trend `generate_trend_raster` can't do; feed it `spectral_index` NBR/NDVI stacks for forest-disturbance mapping. |
| `local_outlier_analysis` | Space-time Local Outlier Analysis — Anselin Local Moran's I on an H3 space-time cube (like ArcGIS's *Local Outlier Analysis*): bin timestamped points into H3 cells × time steps, standardise, and compute each bin's local Moran's I against its space-time neighbourhood (spatial `kring` × ± `time_window`), with a seeded permutation p-value; classify each bin High-High/Low-Low cluster or High-Low/Low-High outlier, and emit one H3 polygon per cell summarising cluster vs. outlier bins over time. The space-time extension of the bundled spatial-only `local_morans_i_lisa`; reuses the cube machinery from `emerging_hot_spot_analysis`. |
| `collapse_hydro_polygon` | Collapse narrow water polygons to centerlines while keeping wide reaches as polygons (like ArcGIS's *Collapse Hydro Polygon*): trace each polygon's centerline along its principal axis (PCA of the boundary), taking the midpoint of the polygon's perpendicular cross-section at stations spaced `sample_distance` apart — which also measures the local width — and collapse to a line where the median width is at or below `collapse_width`, routing wider reaches to an optional `retained` polygon layer. The width-thresholded polygon→line generalization the line-casing `collapse_dual_lines_to_centerline` and raster `river_centerlines` don't do. |
| `change_point_detection` | Detect abrupt shifts in the time series at each location of an H3 space-time cube (like ArcGIS's *Change Point Detection*): bin timestamped points into H3 cells × time steps, then segment each cell's series by binary segmentation on a mean-shift or linear-slope-shift cost — accepting splits whose gain beats a BIC-style penalty scaled by `sensitivity` (`method=auto`) or taking a fixed `num_change_points` (`method=defined`). Each H3 cell reports its change count, the year of its largest change, and the segment means. The temporal segmentation the two-date `change_vector_analysis` and monotonic Mann-Kendall tests can't do; reuses the cube machinery from `emerging_hot_spot_analysis`. |
| `time_series_forecast` | Per-location forecasting on an H3 space-time cube (like ArcGIS's *Exponential Smoothing Forecast* / *Curve Fit Forecast*): bin timestamped points into H3 cells × time steps, then forecast each cell's series `steps` ahead by Holt's exponential smoothing or polynomial (linear/parabolic) curve fits — `model=auto` picks per cell the lowest hold-out-RMSE model. Each H3 cell reports the chosen model, next-step and final-step forecasts, a 90% confidence half-width, and the validation RMSE. Completes the space-time cube suite with the temporal forecasting the bundled spatial `trend_surface` can't do. |
| `reconstruct_tracks` | Turn timestamped points into movement track polylines (like ArcGIS's *Reconstruct Tracks* / *Find Dwell Locations*): group by `track_field`, sort by `time_field`, split on time/distance gaps, and emit per-track stats (duration, length, mean/max speed, haversine distances for geographic CRS). With a `dwells` output it also finds dwell locations — where the mover stayed within a radius for a minimum time. The first trajectory/movement tool in the suite; tracks feed straight into `vector_to_pmtiles` / H3 binning. |
| `emerging_hot_spot_analysis` | Space-time hot spot trends on an H3 space-time cube (like ArcGIS's *Emerging Hot Spot Analysis* + *Create Space Time Cube*): bin timestamped points into H3 cells × time steps, compute the Getis-Ord Gi\* z-score per bin over a space-time neighborhood (spatial k-ring × ± time window), run a Mann-Kendall trend test on each cell's Gi\* series, and classify every cell into the Esri categories (new / consecutive / intensifying / persistent / diminishing / sporadic / oscillating / historical hot or cold spot, or no pattern). Adds the temporal dimension the bundled spatial-only `getis_ord_gi_star` lacks; deterministic, no RNG; output H3 polygons render straight through `h3_to_vector` / PMTiles. |
| `expand_shrink` | Grow (expand) or shrink selected classes of a categorical raster by N cells, leaving other classes intact (like ArcGIS's *Expand* / *Shrink*): iterative 8-connected dilation/erosion where boundary cells adopt the most frequent selected (expand) or non-selected (shrink) neighbour class, with no-data and the raster edge as barriers. The class-aware morphology the binary `buffer_raster`/`nibble` can't do — strengthens the raster→vector cleanup pipeline before `polygonize`. |
| `boundary_clean` | Smooth categorical-raster zone boundaries and remove classification speckle (like ArcGIS's *Boundary Clean* / *Majority Filter*): a **majority** mode replaces each cell with the most frequent value among its 4- or 8-connected neighbours when that value reaches a `majority`/`half` threshold (removing isolated speckle), and an **expand_shrink** mode smooths boundaries with a priority-ordered expansion pass followed by a shrink pass (larger, smaller, or unsorted zones win). No-data cells are barriers; ties break toward the smaller class. The generalization step the bundled `clump`/`nibble` and GeoLibre `expand_shrink` don't cover — cleans classified rasters before `raster_to_vector_polygons` → `regularize_building_footprints` / `smooth_natural_features`. |
| `solar_radiation` | Incoming solar radiation (Wh/m²) over a DEM for a date range (like ArcGIS's *Raster Solar Radiation*): per-cell slope/aspect (Horn), 16-sector horizon shading + sky-view factor, and integration of direct-normal irradiance `I0·E0·τ^m` projected onto the slope plus an isotropic diffuse term, over sampled days × time steps. Optional direct/diffuse component rasters. The flagship insolation tool the bundled `horizon_angle`/`shadow_image` building blocks never integrate into energy units. Deterministic (dates are parameters). |
| `cut_fill` | Volumetric change between two DEM surfaces (like ArcGIS's *Cut Fill* / *Surface Volume*): a signed elevation-change raster (Δz = after − before, or surface − a reference plane) plus cut, fill, and net volumes (Σ \|Δz\| × cell area), with optional contiguous cut/fill region labelling and a per-region volume CSV. For earthworks, erosion/deposition, stockpiles, and lidar change detection. |
| `line_of_sight` | Point-to-point visibility over a DEM (like ArcGIS's *Line Of Sight* / *Construct Sight Lines*): for each observer→target pair, walk the DEM along the segment tracking the running maximum vertical angle, and emit the sight line split into contiguous **visible** and **obstructed** LineString segments plus a target-visible flag. Answers "can A see B, and where does the view break" — the per-pair vector question the raster `viewshed` / `skyline_analysis` can't. Observer/target heights, all-to-all or `pair_field` matching. |
| `corridor` | Least-cost corridor between two accumulated-cost surfaces (like ArcGIS's *Corridor* / Least Cost Corridor): sum two `cost_distance` rasters into a corridor surface whose every cell is the cost of the cheapest A→cell→B path through it, then optionally threshold (absolute or percent-above-minimum) to the swath of near-optimal routes the bundled cost suite (`cost_distance`/`cost_pathway`, which give a single path) can't. Direct (two cost rasters) or convenience mode (one friction raster + two source rasters, accumulated internally). For wildlife corridors and route planning. |
| `interpolate_shape` | Drape points/lines/polygons on a surface raster (like ArcGIS's *Interpolate Shape* + *Add Surface Information*): densify each geometry at a sample interval, sample Z per vertex (bilinear or nearest), write true 3D geometry, and add surface metrics — `z_min`/`z_max`/`z_mean`, 3D surface length, and length-weighted average slope. The bridge between the repo's raster (terrain) and vector halves — trail profiles, pipeline lengths, slope-aware routing — that no bundled point-sampling tool provides for lines/polygons. |
| `generate_transects_along_lines` | Walk each line at a fixed `interval` and emit perpendicular transect lines of a given `length`, centred or offset (like ArcGIS's *Generate Transects Along Lines*): the standard sampling structure for shoreline-change (DSAS-style), riparian surveys, and cross-section extraction — the bundled `points_along_lines` only places points. Each transect carries its parent line id, distance-along-line, and bearing; pairs with `interpolate_shape` for terrain cross-sections. |
| `collapse_dual_lines_to_centerline` | Collapse paired dual carriageways into single centerlines (like ArcGIS's *Collapse Dual Lines To Centerline* / *Merge Divided Roads*): detect roughly parallel line pairs within a width band (directed-overlap parallelism test on densified vertices, optional attribute match), replace each pair with the ordered midpoint centerline, and snap the new endpoints back to surviving lines so the network stays routable. Completes the road-generalization arc with `thin_road_network`; the bundled `river_centerlines` only works from polygons/rasters, not paired line features. |
| `rubbersheet_features` | Warp a vector layer to align with a target using displacement links (like ArcGIS's *Rubbersheet Features* + *Generate Rubbersheet Links*): a Delaunay-TIN piecewise-affine transform (self-contained Bowyer–Watson) that lands control points exactly on their targets and deforms everything between smoothly, with IDW falloff outside the hull (or IDW everywhere). Links come from a `links` line layer or are auto-generated by matching `input` to a `target`. The conflation *transform* step `detect_feature_changes` was built to feed — absent from the bundled suite. |
| `detect_feature_changes` | Match two line datasets and classify each feature's change (like ArcGIS's *Detect Feature Changes*): each update line is matched to the nearest base line by symmetric discrete Hausdorff distance within `search_distance`, then labelled `unchanged` / `spatial` / `attribute` / `spatial_attribute` / `new`, with unmatched base lines emitted as `deleted`. Carries `change_type`, `match_id`, and `match_dist`. The vector change-detection / conflation entry point absent from the bundled suite (its `change_vector_analysis` is raster). |
| `snap_tracks` | Sequence-aware map matching of GPS tracks to a road network (like ArcGIS's *Snap Tracks*): a per-track Viterbi dynamic program trading off snap distance (emission) against route continuity (transition = |along-network distance − GPS movement|, with a big penalty for hopping to a disconnected edge), so matched tracks follow one plausible path instead of zig-zagging between parallel streets like the bundled per-point `snap_points_to_network`. Continues the movement arc of `reconstruct_tracks`. |
| `remove_overlap_multiple` | Reallocate every overlap among a set of polygons to exactly one feature, producing a gap-free, overlap-free partition of their union (like ArcGIS's *Remove Overlap (Multiple)*): the incremental `geo` `BooleanOps` overlay from `count_overlapping_features` finds each shared region, then a grid divides it either by centre line (each point goes to the feature it lies deepest inside) or by nearest generator centroid (Thiessen), conserving total area exactly. `count_overlapping_features` only *counts* overlaps — nothing in the bundled suite resolves them, as trade areas and service territories require. |
| `fuzzy_overlay` | Rescale a raster to a 0..1 fuzzy-membership surface (linear/gaussian/small/large/ms_small/ms_large) or combine several membership rasters with fuzzy AND/OR/PRODUCT/SUM/GAMMA (like ArcGIS's *Fuzzy Membership* + *Fuzzy Overlay*): the standard multi-criteria suitability workflow — site selection, habitat modeling — as pure cell-wise math with no-data propagation. The bundled `weighted_overlay`/`weighted_sum` are crisp reclass-and-add and `fuzzy_knn_classification` is a classifier; nothing does fuzzy suitability. |
| `aggregate_points` | Group points that fall within an aggregation distance of each other into one polygon per cluster (like ArcGIS's *Aggregate Points*): a grid-hashed union-find single-link clustering (equivalent to DBSCAN with `min_samples=1`) whose clusters become convex hulls or buffered-union footprints, each carrying `point_count` and optional per-cluster field sums. The point analogue of `aggregate_polygons`/`delineate_built_up_areas`; the bundled `vector_hex_binning` imposes an arbitrary grid and `concave_hull` produces a single hull for the whole layer. |
| `generate_od_links` | Draw straight origin-destination desire lines from origin points to destination points (like ArcGIS's *Generate Origin-Destination Links*): pair them by a shared `id_field`, within a `search_distance`, and/or to the `num_nearest` (rules combine), each link carrying origin/destination ids and length. The bundled OD tools (`network_od_cost_matrix`, `multimodal_od_cost_matrix`) output cost *tables*, not the geometry used for flow maps and catchment spider diagrams; pairs with `render_vector_png`/`vector_to_pmtiles`. |
| `neighborhood_summary_statistics` | For each feature, summarize numeric fields over its neighbours — k nearest, within a distance band, or shared-edge (rook) contiguity — adding `<field>_nbr_mean/median/std/min/max/sum` and a neighbour count, optionally inverse-distance weighted (like ArcGIS's *Neighborhood Summary Statistics*). The "spatial lag" columns every weights-based workflow (spatial-regression prep, smoothing, anomaly screening) starts from; the bundled suite computes global/local indices (`global_morans_i`, `getis_ord_gi_star`) but never exports the neighbour statistics themselves. |
| `storage_capacity` | Sweep water-surface elevations over a DEM and report the flooded surface area and storage volume at each level — the stage-area-volume curve used in reservoir and detention design — per polygon zone or over the whole DEM (like ArcGIS's *Storage Capacity*). One masked pass per level accumulates `area = Σ cell_area` and `volume = Σ(level − z)·cell_area`, output as a CSV. The bundled `impoundment_size_index` only evaluates dam *siting*; nothing produced the elevation→area→volume table for a given basin. Hydrology-identity fit alongside `cut_fill` and the depression/sink tools. |
| `find_space_time_matches` | Match features between two timestamped point layers that fall within a spatial `search_distance` **and** a `time_window` of each other (before/after/either), emitting matched pairs with `distance` and signed `delta_t` (like ArcGIS's *Find Space Time Matches*): crimes↔calls, sightings↔tracks, incidents↔sensor events. The bundled `spatial_join`/`near` are space-only and `emerging_hot_spot_analysis`/`reconstruct_tracks` handle time within one layer; nothing joined two layers on space × time. Time fields are numeric or ISO-8601. |
| `create_spatially_balanced_points` | Generate spatially well-spread (quasi-random Halton) sample points inside a constraint polygon, optionally weighted by an inclusion-probability raster, each tagged with a balanced `sample_order` so any prefix is itself balanced (like ArcGIS's *Create Spatially Balanced Points*). The bundled `random_points_in_polygon` is plain uniform (clumps and voids); nothing did balanced or probability-weighted sampling — the standard design for field surveys and monitoring networks. Deterministic (seeded, WASM-safe); pairs with `generate_transects_along_lines`. |
| `find_identical` | Find features that are identical on geometry and/or chosen fields, grouping duplicates (like ArcGIS's *Find Identical* / *Delete Identical*): a canonical geometry key (all vertices, optionally snapped to `xy_tolerance`) plus field-value keys, single-pass hashed into groups. `report` adds `dup_group`/`dup_seq` columns; `delete` keeps the first of each group. The bundled `remove_duplicates` is LiDAR-only and `eliminate_coincident_points` is points-only — nothing deduplicated arbitrary vector features after a merge/append. |
| `path_distance` | Accumulated least-cost distance where each 8-connected step pays the true 3-D **surface distance** from a DEM (`√(planar² + Δz²)`) times a slope-dependent **vertical factor** times friction (like ArcGIS's *Path Distance*): Tobler's hiking function (default), linear/sym/inverse-linear, or a binary max-slope cutoff. The bundled `cost_distance` is planar and slope-blind; this extends `corridor`'s Dijkstra with elevation-aware step costs for realistic hiking/wildlife/access travel surfaces. |
| `time_series_clustering` | Cluster the cells of an H3 space-time cube by the similarity of their time series — raw `value`, z-normalized `profile`, or `correlation` (1 − r) — with deterministic multi-restart k-medoids, so "places that evolve alike" group together (like ArcGIS's *Time Series Clustering*). Reuses `emerging_hot_spot_analysis`' H3×time binning, which classifies each cell independently; nothing grouped cells by their whole temporal trajectory. Output is H3 polygons with `cluster_id` and an `is_medoid` flag. |
| `trace_proximity_events` | Find intervals where two moving tracks were within `search_distance` of each other for at least `min_duration` (proximity events — convoy/meeting/contact detection), and optionally trace transitive downstream contacts from seed `entities` with a generation number (like ArcGIS's *Trace Proximity Events*). Positions are linearly interpolated on the union timeline and the squared separation solved exactly (a quadratic per interval) for the within-distance runs; tracing is a temporal BFS over the event graph. `reconstruct_tracks` builds the tracks — this analyzes the interactions between them. |
| `detect_image_anomalies` | Score every pixel of a multiband image by its squared Mahalanobis distance to the scene's (or a moving window's) band statistics — the unsupervised Reed–Xiaoli (RX) anomaly detector — with an optional percentile-threshold mask (like ArcGIS's *Detect Image Anomalies*). The bundled remote-sensing suite *transforms* imagery (`principal_component_analysis`, `minimum_noise_fraction`, `linear_spectral_unmixing`) but never scored anomalies; RX needs no training data and complements `spectral_index`. Covariance inverted with a hand-rolled Gauss–Jordan solve. |
| `resolve_building_conflicts` | Displace, shrink, or hide building footprints that graphically conflict with symbolized road barriers (and each other) for small-scale mapping, tagging each with a `status` (like ArcGIS's *Resolve Building Conflicts*): barriers are buffered by `barrier_width`/2 + `gap`, conflicting buildings are pushed away from the nearest barrier and relaxed apart over several passes, and any that still cannot be placed are shrunk toward `min_size` or hidden. The displacement-cartography piece completing the generalization family (`regularize_building_footprints`, `delineate_built_up_areas`, `thin_road_network`, `collapse_dual_lines_to_centerline`). |
| `thin_road_network` | Generalize a road network for small-scale display (like ArcGIS's *Thin Road Network*): hide short, low-hierarchy roads while preserving connectivity — a candidate road is thinned only if it is not a bridge in the current network, so the visible network keeps its connected components. Non-destructive (flags each road visible/thinned) or filtered. Pairs with `download_osm_vector`. |
| `subdivide_polygon` | Divide each polygon into equal-area parts (a set number, or a target area each) using straight parallel cuts at a given angle (like ArcGIS's *Subdivide Polygon*): rotate into the cut frame, binary-search each cut position where the area to its left hits the target (monotonic → fast bisection), and clip strips with `geo` `BooleanOps`. No bundled equivalent — for parcel pre-division, sampling frames, and splitting oversized polygons for tiling/parallel processing. |
| `split_by_attributes` | Split a vector layer into one output file per unique value (or combination) of the given field(s) (like ArcGIS's *Split By Attributes*): group features by key, write `output_dir/<sanitized_value>.<format>` (geojson/fgb/parquet/shp) preserving schema, geometry type, and CRS, plus a `split_summary.csv`. The per-class / per-year / per-county partitioning chore no bundled tool does. |
| `polygon_neighbors` | Build a polygon adjacency (contiguity) table (like ArcGIS's *Polygon Neighbors*): decompose every boundary into endpoint-keyed edges, and for each pair of polygons emit their shared-border `length` and `node_count` (point-touches) — `length>0` = edge/rook neighbours, `length==0, node_count>0` = corner-only neighbours. `both_sides` for one or two rows per pair, optional snapping. The contiguity data (for spatial weights, redistricting QA, region-merge) the bundled topology-rule checks never export. |
| `count_overlapping_features` | Flatten overlapping polygons into disjoint regions, each attributed with the number of features covering it (like ArcGIS's *Count Overlapping Features*): an incremental `geo` `BooleanOps` overlay (intersection/difference) that tracks coverage depth and the covering feature ids per region, with an optional `min_count` filter and a region→id CSV. For buffer-coverage analysis, service-area redundancy, and imagery-footprint dedup — the bundled `overlaps` is only a boolean predicate. |
| `apportion_polygon` | Transfer numeric attributes from a source polygon layer onto a target polygon layer, apportioned by area of overlap (like ArcGIS's *Apportion Polygon* / Areal Interpolation): each source value is split among the targets it overlaps in proportion to intersection area (or area × a target weight field), normalised so the full value is distributed — dasymetric reaggregation between incompatible zone systems. The value-transfer step `tabulate_intersection` stops short of; reuses the same `geo` `BooleanOps` overlay. |
| `tabulate_intersection` | Vector-on-vector zonal summary (like ArcGIS's *Tabulate Intersection* and the core of *Summarize Within*): apportion a class layer across zone polygons — polygons summarized by intersected area, points by count — reporting the measure, its percentage of each zone, and area-weighted `sum_fields`. Answers "how much of each class, or how many points, fall inside each zone". |
| `directional_distribution` | Descriptive geographic-distribution statistics (like ArcGIS's *Measuring Geographic Distributions*): mean center, median center (Weiszfeld), central feature, standard distance (circle), and the standard deviational ellipse (from the coordinate covariance eigenvectors) — with optional weighting and per-group `case_field` output. |
| `multiple_ring_buffer` | Buffer features at a list of distances into concentric bands (like ArcGIS's *Multiple Ring Buffer*): non-overlapping rings (donuts) or nested disks, optionally dissolved across features per distance, with the band distance stored as an attribute. The everyday proximity-zone tool (catchments, impact bands, setbacks). |
| `aggregate_polygons` | Combine polygons within an aggregation distance into larger polygons (like ArcGIS's *Aggregate Polygons*): a morphological closing (dilate–erode by half the `aggregation_distance`) fuses nearby polygons, with minimum-area and minimum-hole-size filters and an optional `barrier` layer (lines or polygons) that aggregation may not cross. Pairs with `regularize_building_footprints`; the general-purpose sibling of `delineate_built_up_areas`. Each output carries a `part_count` of the source polygons it merges. |
| `delineate_built_up_areas` | Generalize dense building footprints into smooth settlement / urban-extent polygons (like ArcGIS's *Delineate Built-Up Areas*): a morphological closing (dilate–erode by half the `grouping_distance`) fuses footprints within the grouping distance into one built-up area, with rounded, cartographic boundaries. Filter the result by a minimum building count and minimum area, fill small interior holes, and optionally simplify the boundary. The natural sequel to AI building-footprint extraction. |
| `write_geoparquet` | Convert any supported vector format to GeoParquet, Hilbert-sorted with a bbox covering column and ZSTD compression by default. |
| `read_geoparquet` | Read GeoParquet and convert it to another vector format (or store it in memory). |
| `vector_convert` | Convert a vector dataset between formats (the output extension picks the driver). |
| `render_vector_png` | Draw a vector layer (points/lines/polygons) to a PNG map image. |
| `find_argument_statistics` | Per-pixel *argument* statistics over a raster stack — one multiband raster or a list of co-registered rasters (like ArcGIS's *Find Argument Statistics*): for each pixel find the slice index (or `dates` value, e.g. day-of-year) of the maximum (`argmax`), minimum (`argmin`), or ordered-median (`median_position`), or the count of slices past a `threshold`/`comparison` (`duration`) or the longest consecutive run that pass (`longest_run`). Answers "when does NDVI peak", "how many weeks below X", "longest dry spell" — the argument-position questions the value-aggregating bundled cell-statistics and `image_stack_profile` can't, complementing the per-pixel `generate_trend_raster` and `landtrendr`. |
| `detect_incidents` | Flag where a condition starts and ends along each track (like ArcGIS's *Detect Incidents*): group timestamped points by `track_field`, sort by `time_field`, and evaluate a small `start_condition` comparison (`<field> <op> <number>`, optional two-clause `AND`) over a numeric attribute — an incident opens at the first matching point and, without an `end_condition`, spans the maximal run that keeps matching, or with one stays open (points marked `ongoing`) until the end condition fires. Emits every point copied through with incident id / status / sequence / duration (`mode=points`) or one polyline per episode (`mode=segments`). Turns the per-point speed/acceleration from `calculate_motion_statistics` into discrete episodes; nothing bundled extracts condition intervals along a track. |
| `kernel_density_ratio` | Relative-risk ratio of two kernel density surfaces (like ArcGIS's *Calculate Kernel Density Ratio*): evaluate a numerator point layer (cases) and a denominator point layer (population at risk) with the same quartic/biweight kernel and a shared `bandwidth`, then divide cell-by-cell over the padded union extent. A `denominator_floor` writes no-data where the population is too thin so no `inf`/`NaN` leaks, and `log_ratio` natural-log-transforms the surface so over/under-representation is symmetric around 0; haversine metres for a geographic CRS. The density-*ratio* the single-layer bundled `heat_map` cannot express — the standard epidemiology/crime relative-risk map. |
| `pairwise_comparison_weights` | Derive criterion weights from an Analytic Hierarchy Process (AHP) pairwise-comparison matrix on Saaty's 1-9 scale (like ArcGIS's *Assign Weights By Pairwise Comparison*): the weights are the normalized principal eigenvector recovered by power iteration (no linear-algebra crate), emitted as a `criterion, weight, rank` table with a consistency check — principal eigenvalue `lambda_max`, consistency index `CI = (lambda_max−n)/(n−1)`, and consistency ratio `CR = CI/RI[n]` against Saaty's random-index table, warning when `CR > 0.1`. The `matrix` is a JSON 2-D array or a labeled/plain CSV (`input`). The missing front-end of the suitability stack — it *derives* the weights `calculate_composite_index` / `fuzzy_overlay` / bundled `weighted_overlay` consume. |
| `line_density` | Density of linear features — length per unit area — in a circular neighborhood around each raster cell (like ArcGIS's *Line Density*): for every output cell, clip each line segment to the `search_radius` circle with a closed-form segment-circle intersection, sum the clipped lengths (× an optional `weight_field`), and divide by the neighborhood area (πr²). Output is a density raster over the input's radius-padded extent, in per-map-unit or km²/mi²/m² units (`area_units`), with geographic input handled in a local metre frame. The vector-line density the bundled point-only `heat_map` (KDE) and categorical `edge_density` can't produce — for road, stream, and fault density. |
| `feature_outline_masks` | Generate cartographic mask polygons at a `margin` around features (like ArcGIS's *Feature Outline Masks* / *Intersecting Layers Masks*): `exact` buffers the feature outline, `convex_hull` dilates the convex hull, and `box` expands the bounding envelope — each mask fully containing its source and carrying `source_fid`/`mask_kind`. An optional `masked_layer` clips masks to only where a second polygon layer overlaps. The masking toolset the bundled suite lacks (`unsharp_masking` is an image filter); masks feed straight into `render_vector_png` label haloes. |
| `dimension_reduction` | Principal-component analysis over numeric feature attributes (like ArcGIS's *Dimension Reduction*): optionally z-score the selected `fields`, form the correlation (standardized) or covariance matrix, and eigen-decompose it with a hand-rolled cyclic Jacobi rotation — no linear-algebra crate — then project features onto the leading components, appending `PC1..PCk` scores and emitting an eigenvalue / variance-explained / cumulative / per-variable-loadings report table. Keep count comes from `num_components`, a `min_variance` cumulative target, or all components. The bundled `principal_component_analysis` is imagery-band-only; this decorrelates the attribute inputs that feed `calculate_composite_index`, `similarity_search`, and `build_balanced_zones`. |
| `local_bivariate_relationships` | For each feature, classify *how* two variables relate across its local neighbourhood — not significant, positive/negative linear, concave, convex, or undefined (like ArcGIS's *Local Bivariate Relationships*): gather each feature's `neighbors` nearest points, fit local linear `y=a+b·x` and quadratic `y=a+b·x+c·x²` models, measure dependence by the Gaussian entropy reduction `ΔH=-½·ln(1−R²)`, and test significance with a seeded conditional permutation test before classifying the form from the fits. The local, form-aware complement to the global `bivariate_spatial_association` (Lee's L) and the continuous-variable counterpart of the categorical `colocation_analysis`; no bundled tool maps local bivariate form. |
| `grid_index_features` | Build a map-book / atlas index grid of page-sized rectangles with page names (like ArcGIS's *Grid Index Features*, plus a `mode=strip` variant of *Strip Map Index Features*): tile an extent (explicit or from the input layer) into pages sized by `tile_width`/`tile_height` or derived from a paper `page_size` preset × `map_scale`, aligned to an `origin`, named `alphanumeric` (column letter + row number, `A1`/`B3`) or `sequential`, with `intersect_only` dropping pages that miss the data (bbox prefilter + `geo` intersection). Strip mode walks a route line placing overlapping rectangles rotated to the local bearing. Unlike the bundled bare-fishnet `rectangular_grid_from_*`, it carries real page semantics (`page_name`, `row`, `col`, `page_number`) for cartographic map-series output. |
| `repair_geometry` | Detect and fix invalid polygon geometry (like ArcGIS's *Repair Geometry* / *Check Geometry*): remove duplicate/consecutive vertices, drop degenerate rings and null/empty parts, re-wind rings to the OGC convention (exterior CCW, holes CW) with `geo`'s `Orient`, and resolve self-intersections by taking each polygon's `geo` `unary_union` self-union — a bow-tie splits into its two clean lobes. Validity is judged by `geo`'s OGC `Validation`. With `check_only` it reports a per-feature `problem_code`/`problem_desc` instead of editing. The bundled `clean_vector` only *drops* bad geometries; this actually repairs them. |
| `convert_coordinate_notation` | Convert each feature's coordinate between decimal degrees (DD), degrees-minutes-seconds (DMS), degrees-decimal-minutes (DDM), UTM, and MGRS/USNG, writing a new field in the target notation (like ArcGIS's *Convert Coordinate Notation*): read the point geometry as lon/lat or parse a coordinate string field, and optionally rebuild geometry as an EPSG:4326 point. Pure grid math — a closed-form Krüger transverse-Mercator series (WGS84) for UTM plus MGRS grid-zone/latitude-band and 100km-square lettering — so no PROJ and no new dependency; nothing in the bundled suite converts coordinate *strings*, its projection tools reproject whole rasters/layers. |
| `interpolate_with_barriers` | Interpolate point measurements to a raster while respecting absolute barriers — shorelines, ridges, faults, walls (like ArcGIS's *Kernel Interpolation With Barriers* / *Spline With Barriers*): rasterize the barrier polylines/polygons into an impassable mask, replace the straight-line sample→cell distance with a **cost/geodesic distance** from a per-source 8-connected Dijkstra over the free space (the least-cost engine of `path_distance`/`cost_connectivity`/`corridor`), then apply IDW (`w=1/d^power`) or a first-order local-polynomial Gaussian-kernel fit on those geodesic distances; cells unreachable within `radius` become no-data. The barrier-aware interpolator the barrier-blind bundled `idw_interpolation`, kriging, `natural_neighbour_interpolation`, and `thin_plate_spline` can't be. |
| `generalized_linear_regression` | Fit one global regression of a dependent field on explanatory fields by iteratively reweighted least squares (like ArcGIS's *Generalized Linear Regression*): `gaussian` (OLS), `poisson` (log link, counts), or `logistic` (logit link, 0/1). Writes `glr_estimated` / `glr_residual` / `glr_std_resid` per feature and a full diagnostics report — per-term coefficient, standard error, t/z, probability and VIF for multicollinearity, plus AICc, R²/deviance & McFadden pseudo-R², and the studentized Koenker (Breusch–Pagan) heteroskedasticity test (optional `report` CSV). The global companion to the local `geographically_weighted_regression`; the bundled regressions are raster image classifiers, not attribute regression with inferential statistics. |

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
