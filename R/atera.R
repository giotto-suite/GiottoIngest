#' Stream Atera grid cells into tile partitions
#'
#' The zarr input: Rust counterpart of GiottoDisk's per-grid-cell Atera v1
#' transcript writer (`.atera_tx_cell_worker()`). The caller supplies
#' GiottoDisk's Atera tile plan and lookups (`.atera_tx_tile_plan()`,
#' `.atera_tx_lookups()`). Each
#' grid cell is decoded, routed to its `s` x `s` store tiles and written as
#' `tile_index=<k>/part-NNNNN-0.parquet` under `root`; `threads` cells are
#' decoded at once in this process.
#'
#' @param zarr_path path to `transcripts.zarr.zip`
#' @param tiles grid cell names (`"i,j"`), largest first
#' @param ords zero-based grid cell ordinals, for `transcript_id`
#' @param s tiles per side of each grid cell's split
#' @param tile_index list of store tile indices per grid cell (`s^2` each)
#' @param gene_names,cw_gene the archive's gene names and codeword to gene
#'   mapping
#' @param root directory the `tile_index=` partitions are written under
#' @param qv_threshold minimum quality score kept (`NULL` keeps all)
#' @param flip negate y
#' @param window_rows rows decoded per read
#' @param buffer_rows rows buffered per grid cell before a write
#' @param threads grid cells decoded at once
#' @returns one list per grid cell: `tile_index`, `counts`, and `bounds`
#'   (`s^2` x 4 matrix, xmin xmax ymin ymax)
#' @export
gi_atera_tx_cells <- function(zarr_path, tiles, ords, s, tile_index, gene_names,
    cw_gene, root, qv_threshold = 20, flip = TRUE, window_rows = 4194304,
    buffer_rows = 4194304, threads = 1L) {
    n <- length(tiles)
    if (length(ords) != n || length(s) != n || length(tile_index) != n) {
        stop("[GiottoIngest] `tiles`, `ords`, `s` and `tile_index` need one entry per grid cell",
            call. = FALSE)
    }
    if (any(lengths(tile_index) != as.integer(s)^2L)) {
        stop("[GiottoIngest] each grid cell needs s^2 tile indices", call. = FALSE)
    }
    res <- atera_tx_cells_rs(normalizePath(zarr_path), as.character(tiles),
        as.integer(ords), as.integer(s), lapply(tile_index, as.integer),
        as.character(gene_names), as.integer(cw_gene), path.expand(root),
        if (is.null(qv_threshold)) NA_real_ else as.numeric(qv_threshold),
        isTRUE(flip), as.numeric(window_rows), as.numeric(buffer_rows),
        .threads(threads))
    lapply(res, function(r) {
        r$bounds <- matrix(r$bounds, ncol = 4L, byrow = TRUE,
            dimnames = list(NULL, c("xmin", "xmax", "ymin", "ymax")))
        r
    })
}
