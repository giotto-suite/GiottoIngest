# giottoingest

An optional Rust backend for writing GiottoDisk tile stores. It replaces the slow, R-heavy parts of the generic tile writer and reproduces the R writer's output tile for tile. It's a prototype: validated against GiottoDisk, but not yet wired in. It's called from bench scripts, and GiottoDisk is unchanged.

Also bundled: `prototypes/ateratxrs`, a separate R package with a Rust decoder for Atera v1 `transcripts.zarr.zip`. It is the zarr input, still to be folded into the shared core.

## Scope and decisions
- **Always optional.** GiottoDisk keeps an R-native method for every ingest step. That method is the reference, the fallback, and the default when giottoingest isn't installed. Nothing in GiottoDisk may require Rust. The planned switch is `giottodisk.ingest_backend = c("auto", "rust", "r")`: `"auto"` uses Rust only when it's installed, and `"r"` forces the R path.
- **Goal:** much faster ingestion, and far less memory. The R path's parallelism is limited by memory, so the gain is largest on smaller machines.
- **One package,** with pluggable inputs and a shared core. Coverage follows the measured bottlenecks; it doesn't try to cover everything. `ateratxrs` folds in as the zarr input. `cosmxscan` comes in only if its wide-CSV path becomes a bottleneck worth sharing, and only with its author's agreement.
- **Transcripts and polygons** across the major platforms. Expression matrices only where they're a measured bottleneck.
- **Tiles are identical to R.** The tile plan is part of GiottoDisk's format, so it's reproduced exactly; only the counting is faster.
- **WKB is optional, in Rust only.** `gi_write_point_tiles(geom = FALSE)` writes coordinates without the WKB `geom` column. GiottoDisk's point stores still expect `geom` (store validity, the SedonaDB/DuckDB engines in `spatRelate` and `aggregate`), so the default stays `TRUE` until GiottoDisk adopts it.

| Platform | Transcripts | Polygons |
|---|---|---|
| Xenium (parquet) | validated | validated |
| Atera preview (Xenium layout) | validated | validated |
| Atera v1 zarr | validated (`ateratxrs`) | not yet tested |
| VisiumHD, Stereo-seq | bin points validated | flat stores, not tiled |
| MERSCOPE | same path as Xenium, untested (needs the data) | in-memory path, untested |
| CosMx | no disk transcript path in GiottoDisk | GiottoDisk reader fails on the polygon CSV (see Findings) |

## Status

All runs: Apple M4 Max (16 cores, 128 GB), GiottoDisk `dev` at `a41f8ca`, 8 threads. "Identical" means the same tile plan (bounds and counts), the same tile directories, the same rows with identical values in every column, and the same geoparquet metadata.

| Data | Size | Result | Rust | R (same step) |
|---|---|---|---|---|
| Xenium lung transcripts (parquet) | 10.6 M rows kept, 15 tiles | identical | 1.6 s | about 13–18 s |
| Atera breast transcripts (parquet) | 624 M rows kept, 853 tiles | identical plan; all 853 tile counts equal; 40 sampled tiles identical | 2.8 min, 2.0 GB | 10.4 min, 34 GB (8 workers) |
| VisiumHD 2 µm bin points | 9.0 M bins, 11 tiles | identical | 0.5 s | about 4.8 s |
| Stereo-seq testis bin points | 36.7 M bins, 43 tiles | identical | 2.5 s | about 26.5 s |
| Polygons, Xenium lung and Atera breast, cells and nuclei | 157–170 k polygons, 11 tiles each | identical as written | plan and routing 0.2–0.4 s, plus R build 5.2–9.0 s | 9.5–13.0 s |
| Atera kidney transcripts (`ateratxrs`, zarr) | 1.08 B rows kept, 1,169 tiles | identical, by the same checks | 2.2 min on 4 threads, 3.1 GB | about 18 min, single decoder (16 GB simulated) |

- **Points** (transcripts, bin points): Rust does everything: extent, the counting passes for the plan, row routing, and writing each tile.
- **Polygons:** Rust makes the plan and assigns polygons to tiles. R still builds the polygons with GiottoDisk's own tile code, because reproducing terra/GEOS validity repair isn't realistic. What's left of the polygon time is that R build, run one tile at a time.

**Without WKB** (`geom = FALSE`), the Rust write gets only slightly faster; the gain is in size. On Atera breast the store shrinks from 19.0 to 14.9 GB and the write from 106 to 100 s. Xenium lung shrinks from 497 to 379 MB, with identical rows in every remaining column.

## What the R writer imposes (equivalence rules)
- **The plan is part of the format.** `tilework::quadtreePlan`'s control flow is kept: the R loop, `.subdivide_extent` and `.merge_rect_tiles`. Only the counting moves to Rust, one parallel pass per depth instead of one scan per tile.
- **Edges:**
  - Left and top edges are inclusive.
  - Bottom and right edges are inclusive only on the outer boundary of the tile set used *in that pass*. That's the set being counted, or the non-empty tiles being written, never the global extent. Inclusion is decided on the double bounds.
- **float32 columns:** Arrow rounds each bound to float32 before comparing, so Rust compares against the rounded bound. Filter values are rounded the same way.
- **Types after the R round trip:** float → double, uint8 and int16 → int32, uint32 → double; the reader's string casts (`transcript_id` etc.) are kept.
- **Polygons:**
  - The plan counts vertices.
  - Each polygon goes to the tile holding its box center. With float32 columns, that center is `(xmin + xmax) / 2` with the sum done in float32.
  - Vertices are written in `row_index` order.
- **Not reproducible, so not a target:** the row order inside a tile, and the `row_index` derived from it. Two R runs on the same file give different orders, because Arrow's threaded scan returns batches as they finish. Rows are compared by key (`transcript_id`, `feat_ID`). The Rust writer keeps source order, so it's deterministic.

## Reader profile: where R spends its time

Each reader call ran into a fresh project. "seq" is a sequential plan, profiled; "w8" uses 8 mirai workers. Memory is the peak macOS footprint.

| Job | Input | seq | w8 | Peak (w8) | Where the time goes (seq) |
|---|---|---|---|---|---|
| Atera breast transcripts | parquet, 740 M rows | 34.0 min | 10.4 min | 34.4 GB | R tile work 41%, quadtree count scans 22%, per-tile scans 22%, write 11% |
| Xenium lung transcripts | parquet, 12.2 M rows | 18.0 s | 13.5 s | 18.6 GB | R tile work 57%, write 21%, scans 16% |
| Xenium lung transcripts | zarr | 22.0 s | 16.3 s | 20.7 GB | as above, plus zarr→parquet 12% |
| Xenium lung cell polygons | parquet, 4.1 M vertices | 9.0 s | 13.8 s | 15.4 GB | polygon build 45%, plan scans 18%, future overhead 14%, tile scans 12% |
| Atera breast cell polygons | parquet, 4.2 M vertices | 10.8 s | 14.7 s | 15.2 GB | similar |
| Atera kidney expression | zarr | 43.4 s | 45.6 s | 5.2 GB | parquet write 31%, R zarr decode about 28% |
| CosMx expression, `cosmxscan` | csv.gz | 21.6 s | 23.0 s | 2.2 GB | scanner 58%, write 20% |
| CosMx expression, R fallback | csv | 80.8 s | 84.0 s | 3.2 GB | `readLines` 31%, data.table work about 30% |
| VisiumHD 2 µm expression | h5 | 25.1 s | 22.3 s | 13.5 GB | h5 read 37%, write 20% |
| VisiumHD 2 µm bin points | h5 + parquet | 31.5 s | 46.0 s | 15.7 GB | expression ingest 54%, tile path about 30% |
| Stereo-seq testis expression | gef | 10.7 s | 11.4 s | 3.4 GB | write loop 97% |
| Stereo-seq testis bin points | gef | 136.8 s | 142.2 s | 27.2 GB | GEF ingest in R (`paste0`, joins) about 70%, tile path 19% |

## Findings
1. **Workers hurt small inputs.** Polygons and VisiumHD bin points are slower with 8 workers than without, and the 8 workers alone take about 8 GB. Threads in one process avoid both costs.
2. **Bug: CosMx polygon CSV.** Arrow reads a column as int64, but the file contains `1e+05` (row 513,945, column 6). `.cosmx_poly_disk` (`R/convenience-cosmx.R`, around line 458) fails on this vendor file.
3. **Possible lost polygons.** The plan counts vertices but polygons are assigned by box center, so a polygon whose center falls in a tile with no vertices would be dropped. This didn't occur in the four layers tested. `gi_route_polygons()` reports a count of such polygons (`unassigned`).
4. **Stereo-seq's remaining cost is the GEF matrix ingest** (about 98 s), not tiling. It may be fixable in R before reaching for Rust and HDF5.

## Design and integration
- **One package with pluggable inputs feeding a shared core.** The core is the plan counts, the row router, the point writer, polygon routing and, later, an expression writer.
- **Inputs, in order:**
  1. parquet (done)
  2. Atera and Xenium zarr (`ateratxrs`)
  3. CSV and csv.gz, including wide matrices (absorbing or wrapping `cosmxscan`)
  4. HDF5 only if profiling justifies the dependency
- **GiottoDisk integration:**
  - The GiottoKernels pattern: an optional Suggest, an option `giottodisk.ingest_backend = c("auto", "rust", "r")`, and the R path kept as the reference and fallback.
  - The hook points are `storeWrite(parquetGeomTileStore, queryableStore)` for points and `storeWrite(parquetGeomTileStore, parquetStore)` for polygons, in `R/methods-storeWrite.R`.
  - The R side keeps tilework's control flow and calls `gi_count_tiles()` per depth, then `gi_write_point_tiles()` or `gi_route_polygons()`.
- **Distribution:** r-universe binaries. Shared crate versions: extendr 0.9, arrow/parquet 55.

## Next (ranked by measured cost)
1. **MERSCOPE transcripts and polygons:** the tutorial reports 31 min, nearly all of it tiling 151 M transcripts. Needs the data.
2. **Fold `ateratxrs` into the core** as the zarr input (about 18 → 2 min on Atera kidney).
3. **Parallel R polygon build** across tiles: 75–95% of what's left of polygon time.
4. **Cache the filtered coordinates** across the counting passes, about 1 min on Atera breast.
5. **Stereo-seq GEF matrix ingest** (about 98 s). Look at an R-side fix first.
6. **Wide-CSV matrices** (MERSCOPE `cell_by_gene.csv`; CosMx without cosmxscan, 81 s against 22 s). This is where cosmxscan would fold in.
7. **CosMx:** the polygon CSV bug is a GiottoDisk fix; transcripts need a GiottoDisk disk reader first.
8. **Wire into GiottoDisk** behind the option, keeping the R path, validated with the bench scripts.

## Requirements
- **Rust toolchain:** Cargo and rustc ≥ 1.70. The first build fetches crates and takes a few minutes.
- **R:**
  - GiottoDisk `dev`, at `a41f8ca` or later. The bench scripts call some of its internal functions, so other versions may differ.
  - tilework, terra, data.table, future, future.mirai, mirai.
  - arrow built with zstd support.
  - Optional: cosmxscan, for the CosMx matrix job.
- **Platform:** developed on macOS. `fp_sampler.sh` (memory sampling) uses macOS `footprint`; the rest is portable.

## Install
```sh
R CMD INSTALL giottoingest
R CMD INSTALL giottoingest/prototypes/ateratxrs   # optional, Atera zarr decoder
```

## Datasets
The bench jobs read public vendor datasets from one root folder, `GI_DATA` (default `~/Datasets`), in this layout:

| Job prefix | Folder under `GI_DATA` |
|---|---|
| `xenium_lung_*` | `Xenium/LungCancer_10X/Xenium_V1_humanLung_Cancer_FFPE_outs` |
| `xenium_mini_*` | `Cropped_datasets/Xenium/xenium_mini_lung_cancer` |
| `atera_bc_*` | `Atera/WTA_Preview_FFPE_Breast_Cancer_outs` ([10x dataset page](https://www.10xgenomics.com/datasets/atera-wta-ffpe-human-breast-cancer)) |
| `atera_kid_*` | `Atera/Atera_v1_WTA_Kidney_non_diseased_FFPE_outs` ([10x dataset page](https://www.10xgenomics.com/datasets/atera-wta-ffpe-human-nd-kidney)) |
| `cosmx_*` | `CosMx/Breast_cancer/Flatfiles_RNA/flatFiles/BreastCancer` |
| `vhd_*` | `visiumHD/Human_Colon_Cancer_v4/outs/binned_outputs` |
| `vhd_mini_*` | `Cropped_datasets/VisiumHD/visiumhd_mini_crc/binned_outputs` |
| `st_*` | `Stereoseq/Mouse_testis` |
| `st_mini_*` | `Cropped_datasets/Stereoseq/stereoseq_mini_mouse_eyeball` |

Edit the paths at the top of `inst/bench/readers/jobs.R` if your folders differ. Start with the mini datasets.

## Workflow
Every comparison needs an R reference store, written by GiottoDisk's own reader. The R writer is compared against Rust on the same input.

```sh
B=giottoingest/inst/bench
export GI_DATA=~/Datasets

# 1. transcripts: run the R reader once (sequential), then compare
Rscript $B/readers/profile_job.R xenium_mini_tx out/xen_mini seq
R_ROOT=$(ls -d out/xen_mini/project/artifacts/*/data/source_id=*)
TX=$GI_DATA/Cropped_datasets/Xenium/xenium_mini_lung_cancer/transcripts.parquet
Rscript $B/points_vs_r.R $TX "$R_ROOT" out/xen_mini_rs 8

# 2. bin points: capture the bins parquet the reader writes, then compare (kind "bins")
Rscript $B/readers/capture_bins.R vhd_mini_tx out/cap_vhd
#   points_vs_r.R out/cap_vhd/bins.parquet <the cap project's points store> out/vhd_rs 8 bins

# 3. polygons: capture the polygon writer's inputs, then compare
Rscript $B/readers/capture_polys.R xenium_mini_poly out/cap_poly
Rscript $B/polygons_vs_r.R out/cap_poly/poly_args.rds out/poly_rs 8
```

- **`points_vs_r.R TX_PARQUET R_TILE_ROOT OUT [threads] [xenium|bins]`:** compares the plan and tiles with the R store. `SKIP_R=1` runs Rust only, for timing; `SAMPLE=n` checks n random tiles row by row on large data; `GEOM=0` writes without WKB and compares every other column.
- **`polygons_vs_r.R`:** the Rust plan and routing, plus R's own polygon build, compared with the captured R store.
- **Reader profiling:** `readers/run_jobs.sh OUT seq|w8 JOB...` times each reader with memory sampling. Job names are in `jobs.R`.
- **Atera zarr decoder:** `prototypes/ateratxrs/inst/bench/` (`run_cell.R` is the R decoder on one grid cell, `run_cell_rs.R` is Rust, `compare*.R` checks the outputs).
