#!/usr/bin/env Rscript
# Whole-store comparison: row counts of every tile, full equality on a sample.
#   compare_full.R R_STORE RS_STORE [n_sample]
suppressPackageStartupMessages(library(arrow))
a <- commandArgs(TRUE)
ns <- if (length(a) >= 3) as.integer(a[[3]]) else 40L
ta <- sort(list.dirs(a[[1]], recursive = FALSE)); ta <- ta[grepl("tile_index=", ta)]
tb <- sort(list.dirs(a[[2]], recursive = FALSE)); tb <- tb[grepl("tile_index=", tb)]
cat("tiles:", length(ta), length(tb), " same names:", identical(basename(ta), basename(tb)), "\n")
na <- vapply(ta, function(d) nrow(open_dataset(d)), numeric(1))
nb <- vapply(tb, function(d) nrow(open_dataset(d)), numeric(1))
cat(sprintf("rows: %s vs %s, tiles with equal counts: %d of %d\n",
            format(sum(na), big.mark = ","), format(sum(nb), big.mark = ","),
            sum(na == nb), length(na)))
set.seed(1)
pick <- sort(sample(seq_along(ta), min(ns, length(ta))))
bad <- 0L
for (i in pick) {
  read <- function(p) {
    x <- as.data.frame(open_dataset(p))
    x <- x[order(x$row_index), , drop = FALSE]; rownames(x) <- NULL; x
  }
  x <- read(ta[i]); y <- read(tb[i])
  if (!identical(x, y)) {
    bad <- bad + 1L
    cat("DIFF", basename(ta[i]), ":", paste(names(x)[!mapply(identical, x, y)], collapse = ", "), "\n")
  }
}
cat(sprintf("sampled %d tiles (%s rows): %d identical\n", length(pick),
            format(sum(na[pick]), big.mark = ","), length(pick) - bad))
