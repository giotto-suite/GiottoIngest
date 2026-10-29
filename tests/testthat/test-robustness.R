# What GiottoDisk relies on to use the Rust path safely: unsupported input is
# a classed error (GiottoDisk then writes in R), forked workers work after
# the parent used Rust, and empty or degenerate sources do not fail.

odd_file <- function() {
    skip_if_not_installed("arrow")
    df <- data.frame(x = 1:10 + 0.5, y = 1:10 + 0.5, gene = "g", lab = letters[1:10],
        ts = as.POSIXct("2020-01-01", tz = "UTC") + 1:10)
    tbl <- arrow::arrow_table(df)
    tbl$lst <- arrow::Array$create(lapply(1:10, function(i) 1:2))
    tbl$dec <- arrow::Array$create(1:10 + 0.25)$cast(arrow::decimal128(10, 2))
    write_pq(tbl)
}

write_one <- function(src, drop, root = tempfile()) {
    gi_write_point_tiles(src, matrix(c(0, 20, 0, 20), 1), "tile_index=001", "gene",
        drop = drop, root = root, geo_json = "{}")
}

test_that("columns of types R's writer would not round trip are giottoingest_unsupported", {
    f <- odd_file()
    src <- gi_source(f, "x", "y")
    others <- c("ts", "lst", "dec")
    for (col in others) {
        expect_error(write_one(src, drop = setdiff(others, col)), class = "giottoingest_unsupported",
            label = col)
    }
    # dropped, they are fine
    expect_identical(write_one(src, drop = others), 10)
})

test_that("a filter on a non-numeric column is giottoingest_unsupported, not 0 rows", {
    f <- odd_file()
    s <- gi_source(f, "x", "y", filters = list(list(col = "lab", op = ">=", value = 1)))
    expect_error(gi_extent(s), class = "giottoingest_unsupported")
    expect_error(gi_count_tiles(s, matrix(c(0, 20, 0, 20), 1)), class = "giottoingest_unsupported")
})

test_that("non-numeric filter values are giottoingest_unsupported", {
    f <- odd_file()
    expect_error(gi_source(f, "x", "y",
        filters = list(list(col = "lab", op = "in", value = c("a", "b")))),
        class = "giottoingest_unsupported")
})

test_that("a classed error leaves no tiles or spill behind", {
    f <- odd_file()
    root <- tempfile()
    expect_error(write_one(gi_source(f, "x", "y"), drop = c("lst", "dec"), root = root),
        class = "giottoingest_unsupported")
    expect_identical(list.files(root, all.files = TRUE, no.. = TRUE), character(0))
})

test_that("forked workers write after the parent used Rust", {
    skip_on_os("windows")
    f <- write_pq(synth_points(4000, digits = 4), chunk_size = 300L)
    src <- gi_source(f, "x", "y")
    ref <- gi_extent(src, 4L)   # parent: thread pool and allocator in use
    b <- rbind(c(0, 50, 0, 50), c(50, 100, 0, 50))
    jobs <- lapply(1:3, function(i) parallel::mcparallel({
        w <- gi_write_points(src, list(bounds = b, n = c(1, 1)), tempfile(), "gene",
            threads = 2L)
        list(ext = gi_extent(src, 2L), counts = w$counts)
    }))
    res <- parallel::mccollect(jobs, wait = TRUE, timeout = 120)
    expect_length(res, 3L)
    for (r in res) {
        expect_identical(r$ext, ref)
        expect_identical(sum(r$counts), ref$n)
    }
})

test_that("a source with every row filtered out has n 0 and writes nothing", {
    f <- write_pq(synth_points(500))
    src <- gi_source(f, "x", "y", filters = list(list(col = "qv", op = ">", value = 1e9)))
    ex <- gi_extent(src)
    expect_identical(ex$n, 0)
    root <- tempfile()
    w <- gi_write_points(src, list(bounds = rbind(c(0, 100, 0, 50)), n = 0), root, "gene")
    expect_identical(w$tile_index, integer(0))
    expect_false(dir.exists(file.path(root, "tile_index=001")))
})

test_that("a single point and points on one line count and write", {
    one <- write_pq(data.frame(x = 5, y = 3, gene = "g", tid = 1L))
    src <- gi_source(one, "x", "y")
    expect_identical(gi_extent(src)$ext, c(5, 5, 3, 3))
    b <- rbind(c(5, 5, 3, 3))
    expect_identical(gi_count_tiles(src, b), 1)
    out <- tempfile()
    expect_identical(gi_write_points(src, list(bounds = b, n = 1), out, "gene")$counts, 1)
    expect_identical(nrow(read_tile(file.path(out, "tile_index=001"))), 1L)

    line <- write_pq(data.frame(x = 5, y = 0:5 + 0.5, gene = "g", tid = 1:6))
    src <- gi_source(line, "x", "y")
    b <- rbind(c(5, 5, 0.5, 5.5))
    expect_identical(gi_count_tiles(src, b), 6)
    out <- tempfile()
    expect_identical(gi_write_points(src, list(bounds = b, n = 6), out, "gene")$counts, 6)
})
