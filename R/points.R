#' Describe a parquet point source
#'
#' The lazy query a GiottoDisk reader builds, reduced to what the Rust passes
#' need: files, coordinate columns, a y flip and numeric row filters.
#'
#' @param files parquet file paths, in the order the reader scans them
#' @param x,y coordinate columns
#' @param flip_y negate y
#' @param filters list of `list(col, op, value)`; `op` one of `>=`, `>`,
#'   `<=`, `<`, `==`, `!=`
#' @returns a list used by the `gi_*` functions
#' @export
gi_source <- function(files, x, y, flip_y = FALSE, filters = list()) {
    list(files = normalizePath(files), x = x, y = y, flip_y = isTRUE(flip_y),
        f_cols = vapply(filters, `[[`, "", "col"),
        f_ops = vapply(filters, `[[`, "", "op"),
        f_vals = as.numeric(vapply(filters, `[[`, 0, "value")))
}

.src_args <- function(src) {
    src[c("files", "x", "y", "flip_y", "f_cols", "f_ops", "f_vals")]
}

.bounds_vec <- function(bounds) as.numeric(t(as.matrix(bounds)))

#' Rows and extent of a point source
#' @param src from [gi_source()]
#' @param threads decoding threads
#' @returns `list(n, ext)`, `ext` as xmin, xmax, ymin, ymax
#' @export
gi_extent <- function(src, threads = 1L) {
    do.call(gi_extent_rs, c(.src_args(src), threads = as.integer(threads)))
}

#' Rows of a point source per tile
#'
#' Tile membership follows GiottoDisk's contiguous `getTile()` over the given
#' set of tiles.
#' @param src from [gi_source()]
#' @param bounds matrix, one row per tile: xmin, xmax, ymin, ymax
#' @param threads decoding threads
#' @returns numeric counts, one per tile
#' @export
gi_count_tiles <- function(src, bounds, threads = 1L) {
    do.call(gi_count_tiles_rs, c(.src_args(src),
        list(bounds_v = .bounds_vec(bounds), threads = as.integer(threads))))
}

#' Write a point tile store's tiles
#'
#' @param src from [gi_source()]
#' @param bounds matrix of the tiles to write (non-empty plan tiles)
#' @param tile_dirs directory name of each tile under `root`
#' @param id_col feature id column, written as `feat_ID`
#' @param to_string columns cast to string, as the reader's schema does
#' @param drop columns left out
#' @param root directory the tiles are written under
#' @param geo_json geoparquet metadata written in every tile
#' @param spill_rows rows buffered in memory before spilling to disk
#' @param geom write the WKB `geom` column and geoparquet metadata, as
#'   GiottoDisk's point stores have them. `FALSE` writes coordinates only
#'   (`x_index`, `y_index`, and z among the metadata columns): faster and
#'   smaller, but not yet a format GiottoDisk reads.
#' @param threads threads
#' @returns numeric rows per tile
#' @export
gi_write_point_tiles <- function(src, bounds, tile_dirs, id_col, to_string = character(),
    drop = character(), root, geo_json, spill_rows = 4e6, geom = TRUE, threads = 1L) {
    spill <- file.path(root, paste0(".spill-", Sys.getpid()))
    do.call(gi_write_point_tiles_rs, c(.src_args(src), list(
        bounds_v = .bounds_vec(bounds), tile_dirs = as.character(tile_dirs),
        id_col = id_col, to_string = as.character(to_string),
        drop_cols = as.character(drop), root = path.expand(root),
        geo_json = geo_json, spill_dir = spill,
        spill_rows = as.numeric(spill_rows), geom = isTRUE(geom),
        threads = as.integer(threads))))
}

#' Route polygon vertices to tiles
#'
#' Each polygon goes to the tile holding its envelope centroid, as
#' GiottoDisk's polygon tile write selects it; its vertex rows are written to
#' one parquet file per tile in `row_index` order.
#' @param src from [gi_source()], over the vertex rows
#' @param bounds matrix of the tiles written (non-empty plan tiles)
#' @param group polygon id column
#' @param out_dir directory for the per-tile files
#' @param tile_files file name per tile
#' @param threads threads
#' @returns `list(counts, groups, unassigned)`: vertex rows per tile, polygons,
#'   and polygons whose centroid fell in no tile
#' @export
gi_route_polygons <- function(src, bounds, group, out_dir, tile_files, threads = 1L) {
    do.call(gi_route_polygons_rs, c(.src_args(src), list(
        bounds_v = .bounds_vec(bounds), group = group, out_dir = path.expand(out_dir),
        tile_files = as.character(tile_files), threads = as.integer(threads))))
}
