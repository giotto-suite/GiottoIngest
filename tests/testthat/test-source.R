# Argument handling and error paths. Errors must surface as R errors,
# never as a crash of the R process.

pts_file <- function() {
    skip_if_not_installed("arrow")
    write_pq(synth_points(500))
}

test_that("gi_source builds the spec it documents", {
    f <- pts_file()
    s <- gi_source(f, "x", "y", flip_y = TRUE,
        filters = list(list(col = "qv", op = ">=", value = 20)))
    expect_identical(s$f_cols, "qv")
    expect_identical(s$f_ops, ">=")
    expect_identical(s$f_vals, list(20))
    expect_true(s$flip_y)
    expect_identical(s$files, normalizePath(f))
})

test_that("gi_source accepts no filters", {
    s <- gi_source(pts_file(), "x", "y")
    expect_length(s$f_cols, 0L)
    expect_length(s$f_vals, 0L)
    expect_type(s$f_vals, "list")
})

test_that("gi_source accepts an integer filter value", {
    # `value = 20L` is a natural way to write a qv threshold
    f <- pts_file()
    expect_no_error(s <- gi_source(f, "x", "y",
        filters = list(list(col = "qv", op = ">=", value = 20L))))
    expect_identical(s$f_vals, list(20))
})

test_that("a missing input file is an R error", {
    s <- gi_source(pts_file(), "x", "y")
    s$files <- file.path(tempdir(), "does-not-exist.parquet")
    expect_error(gi_extent(s), "does-not-exist")
})

test_that("a missing coordinate column is an R error", {
    expect_error(gi_extent(gi_source(pts_file(), "x", "nope")), "nope")
    expect_error(gi_count_tiles(gi_source(pts_file(), "nope", "y"),
        matrix(c(0, 1, 0, 1), 1)), "nope")
})

test_that("a missing filter column is an R error", {
    s <- gi_source(pts_file(), "x", "y", filters = list(list(col = "zz", op = ">", value = 1)))
    expect_error(gi_extent(s), "zz")
})

test_that("an unsupported filter op is an R error", {
    expect_error(gi_source(pts_file(), "x", "y",
        filters = list(list(col = "qv", op = "%in%", value = 1))), "unsupported filter op")
    expect_error(gi_source(pts_file(), "x", "y",
        filters = list(list(col = "qv", op = ">=", value = c(1, 2)))), "needs one value")
})

test_that("bounds whose length is not a multiple of 4 are an R error", {
    s <- gi_source(pts_file(), "x", "y")
    expect_error(gi_count_tiles(s, c(0, 1, 0)))
})

test_that("fewer tile_dirs than tiles is an R error, and writes nothing partial", {
    s <- gi_source(pts_file(), "x", "y")
    root <- tempfile()
    b <- rbind(c(0, 50, 0, 50), c(50, 100, 0, 50))
    expect_error(gi_write_point_tiles(s, b, tile_dirs = "tile_index=1", id_col = "gene",
        root = root, geo_json = "{}"))
    # no orphaned spill parts (.spill-<pid>) or partial tiles left in the store root
    expect_identical(list.files(root, all.files = TRUE, no.. = TRUE), character(0))
})

test_that("fewer tile_files than tiles is an R error", {
    s <- gi_source(pts_file(), "x", "y")
    b <- rbind(c(0, 50, 0, 50), c(50, 100, 0, 50))
    expect_error(gi_route_polygons(s, b, group = "gene", out_dir = tempfile(),
        tile_files = "t1.parquet"))
})

test_that("threads = 0 or negative falls back to one thread", {
    s <- gi_source(pts_file(), "x", "y")
    ref <- gi_extent(s, 1L)
    expect_identical(gi_extent(s, 0L), ref)
    expect_identical(gi_extent(s, -1L), ref)
    b <- matrix(c(0, 100, 0, 50), 1)
    expect_identical(gi_count_tiles(s, b, -1L), gi_count_tiles(s, b, 1L))
})

test_that("threads = NA is an R error or one thread, not a crash", {
    s <- gi_source(pts_file(), "x", "y")
    res <- tryCatch(gi_extent(s, NA_integer_), error = function(e) "error")
    expect_true(identical(res, "error") || identical(res, gi_extent(s, 1L)))
})
