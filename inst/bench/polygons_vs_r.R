#!/usr/bin/env Rscript
# Rust polygon plan + routing, R polygon building, vs GiottoDisk's R writer.
#   polygons_vs_r.R POLY_ARGS_RDS OUT_DIR [threads]
# POLY_ARGS_RDS: written by capture_polys.R (the writer's parquetStore and
# arguments); the R-written tile store sits in the same capture project.
suppressPackageStartupMessages({
  library(arrow)
  library(GiottoDisk)
  library(giottoingest)
})
options(tilework.warn_sequential = FALSE)
a <- commandArgs(TRUE)
cap <- readRDS(a[[1]]); out <- a[[2]]
threads <- if (length(a) >= 3) as.integer(a[[3]]) else 8L
gd <- asNamespace("GiottoDisk"); tw <- asNamespace("tilework")
unlink(out, recursive = TRUE); dir.create(out, recursive = TRUE)
tm <- function(expr) { t0 <- Sys.time(); v <- expr; list(v = v, s = as.numeric(Sys.time() - t0, units = "secs")) }

data <- cap$data   # vertex rows; its read_fun already applies the flip
files <- sort(list.files(data@path, pattern = "parquet$", recursive = TRUE, full.names = TRUE))
src <- gi_source(files, cap$sdimx, cap$sdimy, flip_y = isTRUE(cap$flip_vertical))
threshold <- cap$threshold
seed <- function(e) {
  tiles <- tilework::tilePlan("spatial")
  terra::ext(tiles) <- terra::ext(e[1:2], e[3:4])
  erange <- range(terra::ext(tiles))
  length(tiles) <- round(max(erange) / min(erange)) * 4L
  tiles
}

# ---- plan: R's, and tilework's on the gi_source (Rust counts) ----------------
rs_ex <- tm(gi_extent(src, threads))
r_plan <- tm({
  st <- storeCreate(path = file.path(out, "r_plan"), type = "parquetGeomTile")
  gd$.tile_plan_quadtree(st, data, threshold = threshold, tiles = seed(rs_ex$v$ext), i = NULL,
    dry_run = FALSE, verbose = FALSE,
    qtree_args = list(sdimx = cap$sdimx, sdimy = cap$sdimy, contiguous = cap$contiguous))
})
rs_plan <- tm({
  fp <- tilework::quadtreePlan(src, tiles = seed(rs_ex$v$ext), threshold = threshold,
    threads = threads, verbose = FALSE)
  list(bounds = fp$bounds, n = as.numeric(fp$n_records))
})
rfp <- r_plan$v$store@tiles
cat(sprintf("plan: R %.1f s, Rust %.1f s; %d tiles; bounds equal %s, counts equal %s\n",
    r_plan$s, rs_ex$s + rs_plan$s, nrow(rs_plan$v$bounds),
    identical(unname(rfp@bounds), unname(rs_plan$v$bounds)),
    identical(as.numeric(rfp@metadata$n_records), as.numeric(rs_plan$v$n))))

# ---- route polygons to tiles (Rust) ------------------------------------------
keep <- which(rs_plan$v$n > 0)
b <- rs_plan$v$bounds
tfiles <- sprintf("t%03d.parquet", seq_along(keep))
route <- tm(gi_route_polygons(src, b[keep, , drop = FALSE], group = cap$group_col,
  out_dir = file.path(out, "routed"), tile_files = tfiles, threads = threads))
cat(sprintf("route: Rust %.1f s; %s polygons, %s in no tile; %s vertex rows routed\n", route$s,
    format(route$v$groups, big.mark = ","), format(route$v$unassigned, big.mark = ","),
    format(sum(route$v$counts), big.mark = ",")))

# ---- build each tile in R, as the R writer's tile FUN does --------------------
st <- r_plan$v$store
st@geomtype <- "polygons"
st@path <- file.path(out, "rs_store")
build <- tm(for (I in seq_along(keep)) {
  f <- file.path(out, "routed", tfiles[I])
  if (!file.exists(f)) next
  mem_data <- arrow::open_dataset(f) |>
    dplyr::arrange(dplyr::across(dplyr::any_of(c("source_id", "row_index")))) |>
    dplyr::collect() |>
    dplyr::select(-dplyr::any_of(c("row_index", "source_id")))
  # tile dirs are numbered by plan index, as in R (empty tiles leave gaps)
  tile_store <- storeCreate(path = gd$.tilepath(gd$.idpath(st@path, st@uid), idx = keep[I]),
                            type = "parquetGeom", uid = st@uid)
  fmt <- if (isTRUE(cap$split_geom)) sub("%d", paste0(keep[I], "_%d"), cap$split_geom_fmt) else cap$split_geom_fmt
  do.call(storeWrite, c(list(tile_store, mem_data, type = "polygons", uid_partition = FALSE,
    tile_idx = NULL, split_geom = cap$split_geom, split_geom_fmt = fmt,
    split_geom_sourcename = cap$split_geom_sourcename,
    .arrow_meta = gd$.arrow_meta_add_geoparquet(st)),
    list(id_col = cap$id_col, sdimx = cap$sdimx, sdimy = cap$sdimy, part_col = cap$part_col),
    cap$write_param))
})
cat(sprintf("build: R %.1f s over %d tiles\n", build$s, length(keep)))
cat(sprintf("Rust plan + route + R build: %.1f s\n", rs_ex$s + rs_plan$s + route$s + build$s))

# ---- compare with the R writer's store --------------------------------------
r_root <- gd$.idpath(cap$store@path, cap$store@uid)
s_root <- gd$.idpath(st@path, st@uid)
rd <- sort(list.dirs(r_root, recursive = FALSE)); rd <- rd[grepl("tile_index=", rd)]
sd <- sort(list.dirs(s_root, recursive = FALSE)); sd <- sd[grepl("tile_index=", sd)]
cat("tile dirs equal:", identical(basename(rd), basename(sd)), "\n")
same <- 0L; same_sorted <- 0L
for (k in seq_along(rd)) {
  x <- as.data.frame(open_dataset(rd[k])); y <- as.data.frame(open_dataset(sd[k]))
  if (identical(x, y)) { same <- same + 1L; same_sorted <- same_sorted + 1L; next }
  id <- intersect(c("poly_ID", "feat_ID"), names(x))[1]
  xs <- x[order(x[[id]]), ]; ys <- y[order(y[[id]]), ]
  rownames(xs) <- rownames(ys) <- NULL
  if (identical(xs, ys)) same_sorted <- same_sorted + 1L else
    cat("DIFF", basename(rd[k]), ": rows R", nrow(x), "Rust", nrow(y), ";",
        if (identical(dim(xs), dim(ys))) paste(names(xs)[!mapply(identical, xs, ys)], collapse = ", "), "\n")
}
cat(sprintf("%d of %d tiles identical as written; %d identical by polygon id\n", same, length(rd), same_sorted))
