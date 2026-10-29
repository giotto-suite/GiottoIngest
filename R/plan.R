#' Tile planning with tilework, counted in Rust
#'
#' A `tileApply()` method for [gi_source()] objects: every tile of a pass is
#' counted in one Rust scan of the source, under GiottoDisk's contiguous
#' `getTile()` edge rules. With it, tilework's own
#' [tilework::quadtreePlan()] plans a point store from a source, so the plan
#' is GiottoDisk's by construction:
#'
#' ```r
#' src <- gi_source(files, "x_location", "y_location", flip_y = TRUE)
#' plan <- tilework::quadtreePlan(src, tiles = seed, threshold = threshold,
#'     threads = 8L)
#' ```
#'
#' `seed` is the starting `tilePlan` GiottoDisk builds from the source
#' extent ([gi_extent()]). Only counting is supported: `FUN` is ignored, and
#' `sdimx`, `sdimy` and `contiguous` (the source carries its coordinates;
#' the edge rules are always the contiguous ones).
#'
#' @param x a [gi_source()]
#' @param tiles a `freeTilePlan`
#' @param threads decoding threads
#' @param ... ignored
#' @returns a list with one count per tile
#' @name gi_tileApply
#' @aliases tileApply,gi_source,missing,ANY-method
#' @exportMethod tileApply
#' @importFrom tilework tileApply
NULL

methods::setOldClass("gi_source")

methods::setMethod("tileApply", methods::signature("gi_source", "missing", "ANY"),
    function(x, tiles, threads = 1L, ...) {
        as.list(gi_count_tiles(x, tiles$bounds, threads))
    })
