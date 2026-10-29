# giottoingest

An optional Rust backend that reads spatial omics inputs fast and writes GiottoDisk point tile stores, tile for tile the same as GiottoDisk's R writer. Everything after the write stays in the Giotto Suite. It's a prototype: validated against GiottoDisk, not yet wired in (GiottoDisk is unchanged; the bench scripts call it).

## Scope and decisions
- **Always optional.** GiottoDisk keeps an R-native method for every ingest step. That method is the reference, the fallback, and the default when giottoingest isn't installed. Nothing in GiottoDisk may require Rust. The planned switch is `giottodisk.ingest_backend = c("auto", "rust", "r")`.
- **Read fast, write tiles, stop.** The backend reads, filters, plans (through tilework) and writes; GiottoDisk builds the store object and does everything after.
- **Transcripts first.** Point data (transcripts, bin points) is where the gains are. Polygon routing works and is validated, but it is frozen: R still builds the polygons (terra/GEOS repair), so the end-to-end gain is small.
- **One package, one crate:** a shared core with one module per input format.
- **Tiles are identical to R.** The tile plan is part of GiottoDisk's format: it is tilework's own `quadtreePlan()`, only the counting is in Rust.
- **WKB is optional, in Rust only.** `geom = FALSE` writes coordinates without the WKB `geom` column. GiottoDisk's point stores still expect `geom`, so the default stays `TRUE`.

## Use

```r
library(giottoingest)
src <- gi_source("transcripts.parquet", "x_location", "y_location", flip_y = TRUE,
    filters = list(list(col = "qv", op = ">=", value = 20)))
ex <- gi_extent(src, threads = 8)                      # n and extent
plan <- tilework::quadtreePlan(src, tiles = seed, threshold = threshold,
    threads = 8)                                       # GiottoDisk's plan, Rust counts
gi_write_points(src, plan, root, id_col = "feature_name",
    to_string = c("transcript_id", "cell_id", "feature_name"), threads = 8)
```

`seed` and `threshold` are GiottoDisk's (the starting `tilePlan` over `ex$ext`, and `.auto_threshold(ex$n, "points")`). `root` is the store's `source_id=<uid>` directory. Tiles go to `tile_index=<k>/part-0.parquet`, numbered by plan index (empty tiles leave gaps, as in GiottoDisk).

| Function | What it does |
|---|---|
| `gi_source()` | Files, coordinate columns, y flip, numeric filters (`>=`, `>`, `<=`, `<`, `==`, `!=`, `in`) |
| `gi_extent()` | Rows and extent of the filtered source |
| `tileApply()` method | Counts every tile of a pass in one Rust scan, so `tilework::quadtreePlan(src, ...)` works |
| `gi_write_points()` | Writes the store's tiles from a plan |
| `gi_atera_tx_cells()` | The zarr input: Atera v1 grid cells, on GiottoDisk's Atera plan |
| `gi_route_polygons()` | Frozen: polygon vertices to tiles by envelope centroid; R builds the polygons |
| `gi_count_tiles()`, `gi_write_point_tiles()` | The low-level passes the above use |

## Layout

```
src/rust/src/
  lib.rs              R entry points (extendr), mimalloc allocator
  core/               shared by every input
    morton.rs         GiottoDisk's per-tile row order
    tile_write.rs     point tile part: row_index, x/y_index, WKB geom, parquet writer
    geo.rs            geoparquet metadata, jsonlite's number format exactly
  input/              one module per format
    parquet.rs        parquet sources (filters, flip, several files)
    zarr/             Atera v1 transcripts.zarr.zip (zip + blosc, grid cell decoder)
                      next: csv (delimited, wide matrices)
  tiles.rs            tile membership under getTile's contiguous edge rules
  points.rs           extent, counts, routed write with spill
  polygons.rs         frozen polygon routing
```

## Status

All runs: Apple M2 Ultra (24 cores, 192 GB), macOS 14, GiottoDisk `dev` `b29361a`, tilework `dev` `447503b`, 8 threads (Atera zarr: 4). "Identical" means the same plan (bounds and counts), the same tile directories, the same rows with identical values in every column but `row_index` (R's is not reproducible), and the same geo metadata. Rust time is the whole Rust path (extent, plan, write) in a fresh process; memory is the process's peak footprint. R is GiottoDisk's own reader for the same step, 8 mirai workers unless marked.

| Data | Size | Result | Rust | R reader |
|---|---|---|---|---|
| Xenium lung transcripts (parquet) | 10.6 M rows, 15 tiles | identical | 2.3 s, 3.5 GB | 25.7 s, 2.6 GB (sequential) |
| MERSCOPE colon transcripts (parquet) | 151 M rows, 209 tiles (one empty plan tile) | identical, every tile | 23.6 s, 9.9 GB | 5.7 min, 32.4 GB |
| Atera breast transcripts (parquet) | 624 M rows, 853 tiles | identical plan, all 853 counts equal, 40 sampled tiles identical | 3.7 min, 4.3 GB | 9.7 min, 34.3 GB |
| Atera kidney transcripts (zarr) | 1.08 B rows, 1,169 tiles | all 1,169 counts equal, 40 sampled tiles identical | 2.6 min, 4.0 GB (4 threads) | 4.2 min, 30.6 GB |
| VisiumHD 2 µm bin points | 9.0 M bins, 11 tiles | identical | 0.9 s, 2.2 GB | tiling step 1.0 s (of 45 s, mostly h5) |
| Stereo-seq bin points (tissue GEF) | 5.0 M bins, 7 tiles | identical | 0.4 s, 1.8 GB | tiling step 0.8 s |
| Polygons: Xenium lung cells and nuclei, Atera breast cells | 157–170 k polygons, 11 tiles each | identical as written | plan + routing 0.4–0.5 s, + R build 8.3–12.4 s | 11–15 s (plan 2.8–3.2 s + build) |

- **Points** (transcripts, bin points): Rust does everything: extent, the counting passes (through `quadtreePlan()`), routing, and writing each tile.
- **Polygons:** Rust plans and routes; R builds the polygons, which is most of the time, so the step barely changes.
- **Without WKB** (`geom = FALSE`) the write is only slightly faster; the gain is store size (Atera breast 19.0 to 14.9 GB).

## What the R writer imposes (equivalence rules)
- **The plan is part of the format.** It is `tilework::quadtreePlan()` itself (subdivide, merge), driven through a `tileApply()` method; only the counting is in Rust, one parallel pass per depth instead of one scan per tile.
- **Edges:** left and top inclusive; bottom and right inclusive only on the outer boundary of the tile set *of that pass* (the tiles being counted, or the non-empty tiles being written). So only non-empty tiles are passed to the writer.
- **Tile dirs** are numbered by plan index; empty plan tiles leave gaps.
- **float32 columns:** Arrow rounds each bound (and filter value) to float32 before comparing; so does Rust.
- **Nulls:** a null coordinate is in no tile and not in the extent, but counts in `n` (R's `nrow()`), which sets the threshold.
- **Several files:** the first file's schema, as `open_dataset()`: later files cast to it, missing columns null.
- **Types after the R round trip, per tile:** float to double; small integers to int32; uint32/uint64 to int32 when the tile's values fit, else double; int64 to int32 when they fit, else int64; strings (also large and view) to string; factors stay dictionaries with R's index type (int8 below 127 levels, int16 below 32,767, else int32); the reader's string casts are kept.
- **Factor level order** is the one known difference: R restores `levels()` from its own `"r"` schema metadata (only present on R-written files); a parquet dictionary holds levels in order of first appearance. Values and types are identical.
- **Geo metadata:** `.geoparquet_metadata()` via jsonlite, `digits = 4` (scientific below 1e-5 and above 1e16, round half to even).
- **Polygons:** the plan counts vertices; each polygon goes to the tile holding its box center (float32 sum for float32 columns); vertices in `row_index` order; a null id is one polygon.
- **Not reproducible, so not a target:** R's row order inside a tile (and `row_index`), which changes between runs because Arrow's threaded scan returns batches as they finish. Rows are compared by key. Rust keeps source order, so it's deterministic.

## Findings
1. **macOS allocator.** With the system allocator, the writers' large reallocations become kernel `vm_copy` calls that contend across threads: the same Xenium lung write took 3–363 s (363 s was 2,590 s of kernel time). mimalloc as the Rust global allocator makes it a steady 2.1 s.
2. **Workers hurt small inputs** in the R path: polygons and VisiumHD bin points are slower with 8 workers, and the workers alone take about 8 GB.
3. **Possible lost polygons** (R and Rust alike): the plan counts vertices but polygons go by box center, so a polygon whose center falls in a tile with no vertices is dropped. Not seen in the layers tested; `gi_route_polygons()` reports it as `unassigned`.
4. **GiottoDisk/tilework:** the starting grid fails for a zero-width or zero-height extent (a single point, or points on one line). Fixed on branches `fix/seed-grid-zero-extent` in both repos (not yet merged).
5. **CosMx polygon CSV** (GiottoDisk): Arrow reads a column as int64 but the file has `1e+05` (row 513,945, column 6); `.cosmx_poly_disk` fails on this vendor file.
6. **Stereo-seq's remaining cost is the GEF matrix ingest**, not tiling.

## Next
1. **Wire into GiottoDisk** behind the option: `storeWrite(parquetGeomTileStore, queryableStore)` builds a `gi_source` from the reader's query, plans with `quadtreePlan(src, ...)`, writes with `gi_write_points()`; the Atera zarr reader calls `gi_atera_tx_cells()`.
2. **csv input** (delimited transcripts, wide matrices: MERSCOPE `cell_by_gene.csv`, CosMx without cosmxscan).
3. **Cache the filtered coordinates** across the counting passes (the plan is about a third of Atera breast's time).
4. **Parallel R polygon build** across tiles, in GiottoDisk, if polygons need to be faster.

## Requirements
- **Rust:** Cargo and rustc ≥ 1.70, and a C compiler (blosc). The first build fetches crates and takes a few minutes.
- **R:** tilework (Imports). For tests and benches: arrow (with zstd), GiottoDisk `dev`, terra, dplyr, future, rlang, testthat; data.table, future.mirai, mirai for the reader profiles.
- **Platform:** developed on macOS. `fp_sampler.sh` (memory sampling) uses macOS `footprint`; the rest is portable.

## Install
```sh
R CMD INSTALL giottoingest
```

## Tests
`tests/testthat` compares every pass with GiottoDisk's own code on synthetic data (edge rules, float32, filters, nulls, types, several files, spill, threads, Morton order, the plan, tile naming, geo metadata), plus the tiny Atera zarr dataset against GiottoDisk's R decoder (skipped without it; set `GI_DATA`). Rust unit tests: `cargo test` in `src/rust` (jsonlite number format).

## Datasets
The bench jobs read public vendor datasets from `GI_DATA` (default `~/Documents/Datasets`):

| Job prefix | Folder under `GI_DATA` |
|---|---|
| `xenium_lung_*` | `Xenium/LungCancer_10X/Xenium_V1_humanLung_Cancer_FFPE_outs` |
| `xenium_mini_*` | `Xenium/LungCancer_10X_small/Xenium_V1_humanLung_Cancer_FFPE_outs` |
| `atera_bc_*` | `Atera/WTA_Preview_FFPE_Breast_Cancer_outs` ([10x](https://www.10xgenomics.com/datasets/atera-wta-ffpe-human-breast-cancer)) |
| `atera_kid_*` | `Atera/Atera_v1_WTA_Kidney_non_diseased_FFPE_outs` ([10x](https://www.10xgenomics.com/datasets/atera-wta-ffpe-human-nd-kidney)) |
| `atera_tiny_*` | `Atera/tiny_atera_dataset_outs` |
| `merscope_*` | `Merscope/Colon_cancer/Vizgen_Colon_Cancer_V2` |
| `vhd_*` | `Visium_HD/Human_Colorectal_Cancer` |
| `st_*` | `StereoSeq/C04687E314_backup/outs` |
| `cosmx_*` | `CosMx/Breast_cancer/Flatfiles_RNA/flatFiles/BreastCancer` |

Edit the paths at the top of `inst/bench/readers/jobs.R` if your folders differ.

## Workflow
Every comparison needs an R reference store, written by GiottoDisk's own reader.

```sh
B=giottoingest/inst/bench
rroot() { for d in $1/project/artifacts/*/data/source_id=*; do ls -d $d/tile_index=* >/dev/null 2>&1 && { echo $d; return; }; done; }

# 1. transcripts: the R reader once, then compare (kinds: xenium, merscope, bins)
$B/readers/run_jobs.sh out seq xenium_lung_tx
Rscript $B/points_vs_r.R $GI_DATA/Xenium/LungCancer_10X/Xenium_V1_humanLung_Cancer_FFPE_outs/transcripts.parquet \
  "$(rroot out/xenium_lung_tx.seq)" out/xen_rs 8

# 2. bin points: capture the bins parquet the reader writes, then compare
Rscript $B/readers/capture_bins.R vhd_tx out/cap_vhd
Rscript $B/points_vs_r.R out/cap_vhd/bins.parquet "$(rroot out/cap_vhd)" out/vhd_rs 8 bins

# 3. polygons: capture the polygon writer's inputs, then compare
Rscript $B/readers/capture_polys.R xenium_lung_poly out/cap_poly
Rscript $B/polygons_vs_r.R out/cap_poly/poly_args.rds out/poly_rs 8

# 4. Atera zarr: the R reader, then the Rust decoder on every grid cell, then compare
$B/readers/run_jobs.sh out w8 atera_kid_tx
Rscript $B/zarr/run_cell_rs.R $GI_DATA/Atera/Atera_v1_WTA_Kidney_non_diseased_FFPE_outs/transcripts.zarr.zip all out/kid_rs 16 4
Rscript $B/zarr/compare_full.R "$(rroot out/atera_kid_tx.w8)" out/kid_rs/store 40
```

- **`points_vs_r.R TX R_TILE_ROOT OUT [threads] [xenium|merscope|bins]`:** the extent, plan and tiles against the R store. `SKIP_R=1` runs Rust only (timing, memory); `SAMPLE=n` checks n random tiles row by row; `GEOM=0` writes without WKB.
- **`polygons_vs_r.R`:** the Rust plan and routing, plus R's own polygon build, against the captured R store.
- **`readers/run_jobs.sh OUT seq|w8 JOB...`:** times each reader with memory sampling. Job names are in `jobs.R`.
- **`zarr/`:** `run_cell.R` (R decoder, one grid cell, profiled), `run_cell_rs.R` (Rust), `compare.R` / `compare_full.R`.
