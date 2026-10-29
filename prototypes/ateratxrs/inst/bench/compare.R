#!/usr/bin/env Rscript
# Are two decoder outputs identical? Compares, per store tile, every column
# of the rows ordered by row_index, plus the returned counts and bounds.
#   compare.R OUT_A OUT_B
suppressPackageStartupMessages(library(arrow))
a <- commandArgs(TRUE)
ra <- readRDS(file.path(a[[1]], "result.rds"))$res
rb <- readRDS(file.path(a[[2]], "result.rds"))$res
stopifnot(identical(ra$tile_index, rb$tile_index))
cat("counts identical:", identical(ra$counts, rb$counts),
    " bounds identical:", identical(ra$bounds, rb$bounds), "\n")
tiles <- function(o) {
  d <- list.dirs(file.path(o, "store"), recursive = TRUE)
  d[grepl("tile_index=", basename(d)) | grepl("^[0-9]+$", basename(d))]
}
ta <- tiles(a[[1]]); tb <- tiles(a[[2]])
stopifnot(identical(basename(ta), basename(tb)))
bad <- 0L
for (i in seq_along(ta)) {
  read <- function(p) {
    x <- as.data.frame(open_dataset(p, format = "parquet"))
    x[order(x$row_index), , drop = FALSE]
  }
  x <- read(ta[i]); y <- read(tb[i])
  rownames(x) <- rownames(y) <- NULL
  if (!identical(x, y)) {
    bad <- bad + 1L
    cat("DIFF in", basename(ta[i]), ":",
        paste(names(x)[!mapply(identical, x, y)], collapse = ", "), "\n")
  }
}
cat(sprintf("%d of %d tiles identical\n", length(ta) - bad, length(ta)))
