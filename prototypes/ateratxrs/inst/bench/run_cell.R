#!/usr/bin/env Rscript
# Run the Atera transcript decoder on one grid cell, in-session, profiled.
#   run_cell.R ZARR TILE OUT_DIR [free_gb] [patch.R]
# Sizes window/buffer as the reader would for a single decoder with
# `free_gb` free (default 16). patch.R, if given, is sourced into the
# GiottoDisk namespace first (prototype replacements of internals).
suppressPackageStartupMessages({
  library(arrow)
  library(GiottoDisk)
  library(data.table)
})
a <- commandArgs(TRUE)
zarr <- a[[1]]; tile <- a[[2]]; out <- a[[3]]
free_gb <- if (length(a) >= 4) as.numeric(a[[4]]) else 16
patch <- if (length(a) >= 5) a[[5]] else NA
ns <- asNamespace("GiottoDisk")
if (!is.na(patch)) {
  env <- new.env(parent = ns)
  sys.source(patch, envir = env)
  for (f in ls(env, all.names = TRUE)) assignInNamespace(f, get(f, env), "GiottoDisk")
  cat("patched:", ls(env, all.names = TRUE), "\n")
}
assignInNamespace(".sc_free_ram", local({b <- free_gb * 1024^3; function() b}),
                  "GiottoDisk")

src <- ns$.zarr_open(zarr)
n_total <- as.numeric(ns$.zarr_attrs(src, "")$number_rnas)
threshold <- ns$.auto_threshold(n_total, type = "points")
plan <- ns$.atera_tx_tile_plan(src, threshold)
ns$.zarr_close(src)
sz <- ns$.atera_tx_sizes(1L, 1L)
if (nzchar(Sys.getenv("WIN"))) sz$window_rows <- as.numeric(Sys.getenv("WIN"))
if (nzchar(Sys.getenv("BUF"))) sz$buffer_rows <- as.numeric(Sys.getenv("BUF"))
n <- attr(plan, "rows")[[tile]]
cat(sprintf("tile %s: %s rows, %d store tiles, window %s, buffer %s\n", tile,
            format(n, big.mark = ","), sum(plan$tile == tile),
            format(sz$window_rows, big.mark = ","), format(sz$buffer_rows, big.mark = ",")))

unlink(out, recursive = TRUE); dir.create(out, recursive = TRUE)
root <- file.path(out, "store")
prof <- file.path(out, "Rprof.out")
gc(reset = TRUE)
pool <- arrow::default_memory_pool()
t0 <- Sys.time()
Rprof(prof, interval = 0.005, event = "elapsed", line.profiling = FALSE)
res <- ns$.atera_tx_cell_worker(zarr, tile, plan, root, "uid0",
  qv_threshold = 20, flip = TRUE, window_rows = sz$window_rows,
  buffer_rows = sz$buffer_rows, compression = ns$.parquet_compression())
Rprof(NULL)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
g <- gc()
cat(sprintf("wall %.1f s, %.2f M rows/s read, %s rows kept, R heap max %.0f MB, arrow max %.0f MB\n",
            secs, n / secs / 1e6, format(sum(res$counts), big.mark = ","),
            sum(g[, ncol(g)]), pool$max_memory / 1e6))
fp <- system2("/usr/bin/footprint", c("-p", Sys.getpid()), stdout = TRUE)
cat(grep("phys_footprint", fp, value = TRUE), sep = "\n")
saveRDS(list(res = res, secs = secs, n = n), file.path(out, "result.rds"))

# inclusive time per stage, from the R call stacks
s <- summaryRprof(prof)$by.total
stage <- c(
  "1 zip read"          = ".zarr_read_raw",
  "2 blosc decompress"  = ".zarr_blosc_decompress",
  "3 decode to vectors" = ".zarr_decode_values",
  "3 column assembly"   = ".zarr_columns",
  "4 unit table"        = ".atera_tx_unit_table",
  "5 tile assignment"   = ".quantile_slot",
  "6 flush/write"       = "flush"
)
tt <- vapply(stage, function(f) {
  k <- paste0('"', f, '"'); if (k %in% rownames(s)) s[k, "total.time"] else 0
}, numeric(1))
total <- summaryRprof(prof)$sampling.time
cat("\ninclusive seconds (nested: decode within columns within unit table):\n")
print(round(data.frame(sec = tt, share = tt / total), 3))
cat("\ntop self time:\n")
print(head(summaryRprof(prof)$by.self, 25))
