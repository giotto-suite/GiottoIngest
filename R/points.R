#' Describe a parquet point source
#'
#' The lazy query a GiottoDisk reader builds, reduced to what the Rust passes
#' need: files, coordinate columns, a y flip and numeric row filters.
#'
#' @param files parquet file paths, in the order the reader scans them
#' @param x,y coordinate columns
#' @param flip_y negate y
#' @param filters list of `list(col, op, value)` on numeric columns; `op`
#'   one of `>=`, `>`, `<=`, `<`, `==`, `!=`, or `in` with a vector `value`
#'   (R's `%in%`). Rows where the column is NA are dropped, as in arrow.
#' @returns a list used by the `gi_*` functions
#' @export
gi_source <- function(files, x, y, flip_y = FALSE, filters = list()) {
    ops <- c(">=", ">", "<=", "<", "==", "!=", "in")
    for (f in filters) {
        if (!is.list(f) || !all(c("col", "op", "value") %in% names(f))) {
            stop("[GiottoIngest] each filter is list(col, op, value)", call. = FALSE)
        }
        if (!f$op %in% ops) {
            stop("[GiottoIngest] unsupported filter op `", f$op, "`", call. = FALSE)
        }
        if (f$op != "in" && length(f$value) != 1L) {
            stop("[GiottoIngest] filter `", f$col, " ", f$op, "` needs one value", call. = FALSE)
        }
    }
    s <- list(files = normalizePath(files), x = x, y = y, flip_y = isTRUE(flip_y),
        f_cols = vapply(filters, `[[`, "", "col"),
        f_ops = vapply(filters, `[[`, "", "op"),
        f_vals = lapply(filters, function(f) as.numeric(f$value)))
    class(s) <- "gi_source"
    s
}

.src_args <- function(src) {
    src[c("files", "x", "y", "flip_y", "f_cols", "f_ops", "f_vals")]
}

.bounds_vec <- function(bounds) {
    b <- as.matrix(bounds)
    if (length(b) == 0L) return(numeric(0))
    if (ncol(b) != 4L || !is.numeric(b) || anyNA(b)) {
        stop("[GiottoIngest] `bounds` must be a numeric matrix with columns ",
            "xmin, xmax, ymin, ymax and no NA", call. = FALSE)
    }
    as.numeric(t(b))
}

.n_bounds <- function(bounds) if (length(bounds) == 0L) 0L else nrow(as.matrix(bounds))

.threads <- function(threads) {
    threads <- suppressWarnings(as.integer(threads))
    if (length(threads) != 1L || is.na(threads)) {
        stop("[GiottoIngest] `threads` must be a single integer", call. = FALSE)
    }
    max(threads, 1L)
}

#' Rows and extent of a point source
#' @param src from [gi_source()]
#' @param threads decoding threads
#' @returns `list(n, ext)`, `ext` as xmin, xmax, ymin, ymax
#' @export
gi_extent <- function(src, threads = 1L) {
    do.call(gi_extent_rs, c(.src_args(src), threads = .threads(threads)))
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
        list(bounds_v = .bounds_vec(bounds), threads = .threads(threads))))
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
    bv <- .bounds_vec(bounds)
    if (length(tile_dirs) != .n_bounds(bounds)) {
        stop("[GiottoIngest] `tile_dirs` needs one name per tile (", .n_bounds(bounds), ")",
            call. = FALSE)
    }
    threads <- .threads(threads)
    root <- path.expand(root)
    # a spill dir per call (concurrent calls in one process do not collide);
    # removed afterwards, and on error so are the tiles this call created
    spill <- tempfile(paste0(".spill-", Sys.getpid(), "-"), tmpdir = root)
    dirs <- file.path(root, tile_dirs)
    pre <- dir.exists(dirs)
    root_pre <- dir.exists(root)
    ok <- FALSE
    on.exit({
        unlink(spill, recursive = TRUE)
        if (!ok) {
            unlink(dirs[!pre], recursive = TRUE)
            if (!root_pre && length(list.files(root, all.files = TRUE, no.. = TRUE)) == 0L) {
                unlink(root, recursive = TRUE)
            }
        }
    }, add = TRUE)
    res <- do.call(gi_write_point_tiles_rs, c(.src_args(src), list(
        bounds_v = bv, tile_dirs = as.character(tile_dirs),
        id_col = id_col, to_string = as.character(to_string),
        drop_cols = as.character(drop), root = root,
        geo_json = geo_json, spill_dir = spill,
        spill_rows = as.numeric(spill_rows), geom = isTRUE(geom),
        threads = threads)))
    ok <- TRUE
    res
}

#' Write a point tile store from a plan
#'
#' The one call a reader needs after planning: writes every non-empty tile
#' of `plan` as GiottoDisk's point store writer does, into
#' `root/tile_index=<k>/part-0.parquet` with `k` the tile's index in the plan
#' (empty tiles leave gaps, as in GiottoDisk), and the store's geoparquet
#' metadata in every tile.
#'
#' Only the non-empty tiles are passed on, as R writes them: the contiguous
#' edge rules (which tiles' right and bottom edges are inclusive) are taken
#' over that set.
#'
#' @param src from [gi_source()]
#' @param plan a `freeTilePlan` from [tilework::quadtreePlan()] (with its
#'   `n_records`), or a list with `bounds` (one row per tile: xmin, xmax,
#'   ymin, ymax) and `n` (rows per tile)
#' @param root the store's `source_id=` directory the tiles go under
#' @param id_col feature id column, written as `feat_ID`
#' @param to_string columns cast to string, as the reader's schema does
#' @param drop columns left out
#' @param geom write the WKB `geom` column and geoparquet metadata
#'   (GiottoDisk's current point format)
#' @param threads threads
#' @param spill_rows rows buffered in memory before spilling to disk
#' @returns `list(tile_index, counts)`: the plan indices written and their
#'   row counts
#' @export
gi_write_points <- function(src, plan, root, id_col, to_string = character(),
    drop = character(), geom = TRUE, threads = 1L, spill_rows = 4e6) {
    if (methods::is(plan, "freeTilePlan")) {
        bounds <- plan$bounds
        n <- plan$n_records
    } else {
        bounds <- plan$bounds
        n <- plan$n
    }
    bounds <- as.matrix(bounds)
    if (is.null(n) || length(n) != nrow(bounds)) {
        stop("[GiottoIngest] `plan` needs one row count per tile", call. = FALSE)
    }
    keep <- which(n > 0)
    # the store's extent: the whole plan, empty tiles included
    ext <- c(min(bounds[, 1L]), max(bounds[, 2L]), min(bounds[, 3L]), max(bounds[, 4L]))
    counts <- gi_write_point_tiles(src, bounds[keep, , drop = FALSE],
        tile_dirs = sprintf("tile_index=%03d", keep), id_col = id_col,
        to_string = to_string, drop = drop, root = root,
        geo_json = gi_geo_json_rs(as.numeric(ext)), spill_rows = spill_rows,
        geom = geom, threads = threads)
    list(tile_index = keep, counts = counts)
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
    bv <- .bounds_vec(bounds)
    if (length(tile_files) != .n_bounds(bounds)) {
        stop("[GiottoIngest] `tile_files` needs one name per tile (", .n_bounds(bounds), ")",
            call. = FALSE)
    }
    do.call(gi_route_polygons_rs, c(.src_args(src), list(
        bounds_v = bv, group = group, out_dir = path.expand(out_dir),
        tile_files = as.character(tile_files), threads = .threads(threads))))
}
