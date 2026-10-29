#!/usr/bin/env Rscript
# The Rust decoder (ateratxrs) on grid cells, sized like run_cell.R.
#   run_cell_rs.R ZARR TILES OUT_DIR [free_gb] [threads]
# TILES: one grid cell ("5,3"), a comma-free list "5,3;1,2", or "all".
suppressPackageStartupMessages({
  library(arrow)
  library(GiottoDisk)
  library(ateratxrs)
})
a <- commandArgs(TRUE)
zarr <- a[[1]]; tiles_arg <- a[[2]]; out <- a[[3]]
free_gb <- if (length(a) >= 4) as.numeric(a[[4]]) else 16
threads <- if (length(a) >= 5) as.integer(a[[5]]) else 1L
ns <- asNamespace("GiottoDisk")
assignInNamespace(".sc_free_ram", local({b <- free_gb * 1024^3; function() b}),
                  "GiottoDisk")

src <- ns$.zarr_open(zarr)
n_total <- as.numeric(ns$.zarr_attrs(src, "")$number_rnas)
threshold <- ns$.auto_threshold(n_total, type = "points")
plan <- ns$.atera_tx_tile_plan(src, threshold)
lk <- ns$.atera_tx_lookups(src)
ns$.zarr_close(src)
sz <- ns$.atera_tx_sizes(1L, 1L)
if (nzchar(Sys.getenv("WIN"))) sz$window_rows <- as.numeric(Sys.getenv("WIN"))
if (nzchar(Sys.getenv("BUF"))) sz$buffer_rows <- as.numeric(Sys.getenv("BUF"))

rows <- attr(plan, "rows")
cells <- if (identical(tiles_arg, "all")) {
  names(sort(rows, decreasing = TRUE))
} else {
  strsplit(tiles_arg, ";", fixed = TRUE)[[1]]
}
ords <- match(cells, unique(plan$tile)) - 1L
s <- vapply(cells, function(t) plan$s[plan$tile == t][1L], numeric(1))
ti <- lapply(cells, function(t) plan$tile_index[plan$tile == t])
cat(sprintf("%d cell(s), %s rows, threads %d, window %s, buffer %s\n", length(cells),
            format(sum(rows[cells]), big.mark = ","), threads,
            format(sz$window_rows, big.mark = ","), format(sz$buffer_rows, big.mark = ",")))

unlink(out, recursive = TRUE); dir.create(out, recursive = TRUE)
t0 <- Sys.time()
res <- atera_tx_cells(zarr, cells, ords, s, ti, lk$gene_names, lk$cw_gene,
  root = file.path(out, "store"), qv_threshold = 20, flip = TRUE,
  window_rows = sz$window_rows, buffer_rows = sz$buffer_rows, threads = threads)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
kept <- sum(vapply(res, function(r) sum(r$counts), numeric(1)))
cat(sprintf("wall %.1f s, %.2f M rows/s read, %s rows kept\n", secs,
            sum(rows[cells]) / secs / 1e6, format(kept, big.mark = ",")))
fp <- system2("/usr/bin/footprint", c("-p", Sys.getpid()), stdout = TRUE)
cat(grep("phys_footprint", fp, value = TRUE), sep = "\n")
# same shape as the R decoder's result, for compare.R
r1 <- res[[1L]]
r1$bounds <- unname(r1$bounds)
saveRDS(list(res = r1, all = res, secs = secs, n = sum(rows[cells])),
        file.path(out, "result.rds"))
