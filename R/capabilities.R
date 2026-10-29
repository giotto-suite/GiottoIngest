#' Capabilities of this build
#'
#' What this GiottoIngest build can do. GiottoDisk gates each fast path on
#' its capability, so an older build without it falls back to R:
#'
#' - `"points_parquet"`: point tile stores from parquet sources
#'   ([gi_source()], [gi_write_points()]).
#' - `"points_atera_zarr"`: Atera v1 `transcripts.zarr.zip` grid cells
#'   ([gi_atera_tx_cells()]).
#'
#' @param name capability name
#' @returns `gi_capabilities()`: character vector; `has_capability()`: `TRUE` or
#'   `FALSE`
#' @export
gi_capabilities <- function() capabilities_rs()

#' @rdname gi_capabilities
#' @export
has_capability <- function(name) {
    isTRUE(name %in% capabilities_rs())
}
