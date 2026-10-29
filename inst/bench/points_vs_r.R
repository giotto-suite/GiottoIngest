#!/usr/bin/env Rscript
# Rust point tile write vs GiottoDisk's R writer.
#   points_vs_r.R TX_PARQUET R_TILE_ROOT OUT_DIR [threads] [kind]
# kind "xenium" (default): 10x transcripts query; "merscope": Vizgen
# detected_transcripts query (reader defaults); "bins": the bins parquet of
# .binpoints_store() (plain open_dataset, feat_ID/x/y).
# R_TILE_ROOT: the `source_id=...` directory of a store written by the R
# reader (`importXeniumDisk(...)$load_transcripts()`), for the tile check.
suppressPackageStartupMessages({
  library(arrow)
  library(GiottoDisk)
  library(GiottoIngest)
})
a <- commandArgs(TRUE)
tx <- a[[1]]; r_root <- a[[2]]; out <- a[[3]]
threads <- if (length(a) >= 4) as.integer(a[[4]]) else 8L
kind <- if (length(a) >= 5) a[[5]] else "xenium"
skip_r <- nzchar(Sys.getenv("SKIP_R"))
drop <- character()
with_geom <- !identical(Sys.getenv("GEOM"), "0")   # GEOM=0: no WKB column
n_sample <- if (nzchar(Sys.getenv("SAMPLE"))) as.integer(Sys.getenv("SAMPLE")) else Inf
gd <- asNamespace("GiottoDisk"); tw <- asNamespace("tilework")
unlink(out, recursive = TRUE); dir.create(out, recursive = TRUE)
tm <- function(expr) { t0 <- Sys.time(); v <- expr; list(v = v, s = as.numeric(Sys.time() - t0, units = "secs")) }

if (kind == "xenium") {
  # the reader's lazy query (.xenium_transcript_disk read_fun)
  read_fun <- function(x, schema = NULL) {
    sc <- arrow::schema(arrow::open_dataset(x))
    cast_cols <- intersect(c("transcript_id", "cell_id", "feature_name"), names(sc))
    schema <- do.call(arrow::schema, lapply(names(sc), function(nm) {
      if (nm %in% cast_cols) arrow::field(nm, arrow::string()) else sc$GetFieldByName(nm)
    }))
    arrow::open_dataset(sources = x, schema = schema) |>
      dplyr::filter(qv >= 20) |>
      dplyr::mutate(y_location = -y_location)
  }
  sx <- "x_location"; sy <- "y_location"; id_col <- "feature_name"; key <- "transcript_id"
  src <- gi_source(tx, sx, sy, flip_y = TRUE,
                   filters = list(list(col = "qv", op = ">=", value = 20)))
  to_string <- intersect(c("transcript_id", "cell_id", "feature_name"),
                         names(arrow::schema(arrow::open_dataset(tx))))
} else if (kind == "merscope") {
  # the reader's lazy query (.merscope_transcript_disk read_fun), defaults
  read_fun <- function(x, schema = NULL) gd$.merscope_tx_query(x,
    dropcols = c("x", "y", "barcode_id", "transcript_id"), flip_vertical = TRUE)
  sx <- "global_x"; sy <- "global_y"; id_col <- "gene"; key <- NA
  src <- gi_source(tx, sx, sy, flip_y = TRUE)
  to_string <- "cell_id"
  drop <- c("", "x", "y", "barcode_id", "transcript_id")
} else {
  read_fun <- function(path, ...) arrow::open_dataset(path, ...)
  sx <- "x"; sy <- "y"; id_col <- "feat_ID"; key <- "feat_ID"
  src <- gi_source(tx, sx, sy)
  to_string <- character()
}
qs <- as(fileStore(path = tx, read_fun = read_fun), "queryableStore")

# ---- extent ----------------------------------------------------------------
r_ex <- if (skip_r) NULL else tm({ aq <- storeRead(qs); list(n = gd$.dplyr_nrow(aq),
            e = tw$.ext_to_num_vec(gd$.dplyr_ext(aq, sx, sy))) })
rs_ex <- tm(gi_extent(src, threads))
if (skip_r) cat(sprintf("extent: Rust %.1f s\n", rs_ex$s)) else
cat(sprintf("extent: R %.1f s, Rust %.1f s; n equal %s, extent equal %s\n", r_ex$s, rs_ex$s,
    identical(r_ex$v$n, rs_ex$v$n), identical(r_ex$v$e, rs_ex$v$ext)))

# ---- plan: GiottoDisk's own, and tilework's on the gi_source (Rust counts) --
n <- rs_ex$v$n
threshold <- gd$.auto_threshold(n, type = "points")
seed <- function(e) {
  tiles <- tilework::tilePlan("spatial")
  terra::ext(tiles) <- terra::ext(e[1:2], e[3:4])
  erange <- range(terra::ext(tiles))
  length(tiles) <- round(max(erange) / min(erange)) * 4L
  tiles
}
r_plan <- if (skip_r) NULL else tm({
  st <- storeCreate(path = file.path(out, "r_plan"), type = "parquetGeomTile")
  gd$.tile_plan_quadtree(st, qs, threshold = threshold, tiles = seed(rs_ex$v$ext), i = NULL,
    dry_run = FALSE, verbose = FALSE,
    qtree_args = list(sdimx = sx, sdimy = sy, contiguous = TRUE))$store@tiles
})

# the same plan through tilework itself, counts in Rust (tileApply method)
rs_plan <- tm({
  fp <- tilework::quadtreePlan(src, tiles = seed(rs_ex$v$ext), threshold = threshold,
    threads = threads, verbose = FALSE)
  list(plan = fp, bounds = fp$bounds, n = as.numeric(fp$n_records))
})
if (skip_r) {
  cat(sprintf("plan: Rust %.1f s; %d tiles\n", rs_plan$s, nrow(rs_plan$v$bounds)))
} else {
  same_b <- identical(unname(r_plan$v@bounds), unname(rs_plan$v$bounds))
  same_n <- identical(as.numeric(r_plan$v@metadata$n_records), as.numeric(rs_plan$v$n))
  cat(sprintf("plan: R %.1f s, Rust %.1f s; %d tiles; bounds equal %s, counts equal %s\n",
      r_plan$s, rs_plan$s, nrow(rs_plan$v$bounds), same_b, same_n))
}

# ---- write ---------------------------------------------------------------
rs_root <- file.path(out, "rs_tiles")
rs_w <- tm(gi_write_points(src, rs_plan$v$plan, rs_root, id_col = id_col, to_string = to_string,
  drop = drop, geom = with_geom, threads = threads))
cat(sprintf("write: Rust %.1f s, %s rows in %d tiles\n", rs_w$s,
    format(sum(rs_w$v$counts), big.mark = ","), length(rs_w$v$tile_index)))
cat(sprintf("Rust total %.1f s (extent + plan + write)\n", rs_ex$s + rs_plan$s + rs_w$s))
cat(sprintf("store size %.0f MB (geom %s)\n", sum(file.size(list.files(rs_root, recursive = TRUE,
    full.names = TRUE))) / 1e6, if (with_geom) "written" else "omitted"))
fp <- system2("/usr/bin/footprint", c("-p", Sys.getpid()), stdout = TRUE)
cat(grep("phys_footprint_peak", fp, value = TRUE), "\n")
if (skip_r) quit(save = "no")

# ---- tile-for-tile check against the R reader's store ----------------------
rd <- sort(list.dirs(r_root, recursive = FALSE)); rd <- rd[grepl("tile_index=", rd)]
sd <- sort(list.dirs(rs_root, recursive = FALSE)); sd <- sd[grepl("tile_index=", sd)]
cat("tile dirs equal:", identical(basename(rd), basename(sd)), "\n")
nr <- function(d) vapply(d, function(p) nrow(open_dataset(p)), numeric(1))
cat("tiles with equal row counts:", sum(nr(rd) == nr(sd)), "of", length(rd), "\n")
set.seed(1)
pick <- if (is.finite(n_sample) && n_sample < length(rd)) sort(sample(seq_along(rd), n_sample)) else seq_along(rd)
# R's per-tile row order (and so row_index) changes from run to run: arrow
# returns scan batches as threads finish. Compare rows by their key, every
# column but row_index.
bad <- 0L
for (k in pick) {
  x <- as.data.frame(open_dataset(rd[k])); y <- as.data.frame(open_dataset(sd[k]))
  skip <- c("row_index", if (!with_geom) "geom")
  # no unique key (MERSCOPE): order by every column but the WKB, which is
  # derived from x_index / y_index
  if (is.na(key)) skip <- c(skip, "geom")
  x <- x[, !names(x) %in% skip]; y <- y[, !names(y) %in% skip]
  ox <- if (is.na(key)) do.call(order, unname(as.list(x))) else order(x[[key]])
  oy <- if (is.na(key)) do.call(order, unname(as.list(y))) else order(y[[key]])
  x <- x[ox, ]; y <- y[oy, ]
  rownames(x) <- rownames(y) <- NULL
  if (!identical(x, y)) {
    bad <- bad + 1L
    cols <- if (identical(dim(x), dim(y))) names(x)[!mapply(identical, x, y)] else
      sprintf("rows R %d, Rust %d", nrow(x), nrow(y))
    cat("DIFF", basename(rd[k]), ":", paste(cols, collapse = ", "), "\n")
  }
}
cat(sprintf("%d of %d checked tiles hold identical rows (all columns but row_index)\n",
    length(pick) - bad, length(pick)))
f1 <- function(d) list.files(d, pattern = "parquet$", full.names = TRUE)[1]
if (with_geom) cat("geo metadata equal:", identical(
  ParquetFileReader$create(f1(rd[1]))$GetSchema()$metadata$geo,
  ParquetFileReader$create(f1(sd[1]))$GetSchema()$metadata$geo), "\n")
