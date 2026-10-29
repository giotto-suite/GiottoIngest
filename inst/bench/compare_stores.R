#!/usr/bin/env Rscript
# Two point tile stores (e.g. GiottoDisk's R path and its GiottoIngest path):
# tile dirs, rows per tile, rows in every column but row_index (as sets:
# R's order within a tile is not reproducible) and the file-level geo
# metadata (DuckDB's parquet_kv_metadata).
#   compare_stores.R ROOT_A ROOT_B [n_sample]
suppressPackageStartupMessages(library(arrow))
a <- commandArgs(TRUE)
n_sample <- if (length(a) >= 3) as.integer(a[[3]]) else Inf
dirs <- function(r) { d <- list.dirs(r, recursive = FALSE); sort(d[grepl("^tile_index=", basename(d))]) }
da <- dirs(a[[1]]); db <- dirs(a[[2]])
cat("tile dirs equal:", identical(basename(da), basename(db)), "(", length(da), ")\n")
nr <- function(d) vapply(d, function(p) nrow(open_dataset(p)), numeric(1))
na <- nr(da); nb <- nr(db)
cat(sprintf("rows: %s vs %s; tiles with equal counts: %d of %d\n", format(sum(na), big.mark = ","),
    format(sum(nb), big.mark = ","), sum(na == nb), length(na)))
set.seed(1)
pick <- if (is.finite(n_sample) && n_sample < length(da)) sort(sample(seq_along(da), n_sample)) else seq_along(da)
rows <- function(d) {
  x <- as.data.frame(open_dataset(d))
  x <- x[, setdiff(names(x), c("row_index", "geom")), drop = FALSE]
  x <- x[do.call(order, unname(as.list(x))), , drop = FALSE]
  rownames(x) <- NULL
  x
}
bad <- 0L
for (k in pick) {
  x <- rows(da[k]); y <- rows(db[k])
  if (!identical(x, y)) { bad <- bad + 1L; cat("DIFF", basename(da[k]), "\n") }
}
cat(sprintf("%d of %d checked tiles hold identical rows (all columns but row_index)\n",
    length(pick) - bad, length(pick)))
con <- DBI::dbConnect(duckdb::duckdb())
geo <- function(d) {
  f <- list.files(d, pattern = "parquet$", full.names = TRUE)[1]
  kv <- DBI::dbGetQuery(con, sprintf(
    "SELECT decode(key) AS key, decode(value) AS value FROM parquet_kv_metadata(%s)",
    DBI::dbQuoteString(con, f)))
  kv$value[kv$key == "geo"]
}
cat("file geo metadata equal:", identical(geo(da[1]), geo(db[1])), "\n")
DBI::dbDisconnect(con, shutdown = TRUE)
