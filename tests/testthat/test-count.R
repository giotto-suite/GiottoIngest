# gi_count_tiles: per-tile counts under GiottoDisk's contiguous getTile edge
# rules over the given set of tiles, against R's own per-tile filtered scan.

expect_counts_match <- function(x, bounds, sx = "x", sy = "y", flip = FALSE, filters = list(),
                                chunk_size = 250L) {
    f <- if (is.character(x)) x else write_pq(x, chunk_size = chunk_size)
    ref <- r_count_tiles(make_qs(f, sx, sy, flip, filters), bounds, sx, sy)
    src <- gi_source(f, sx, sy, flip_y = flip, filters = filters)
    for (th in c(1L, 4L)) {
        expect_identical(gi_count_tiles(src, bounds, th), ref)
    }
    invisible(ref)
}

grid_bounds <- function(xs, ys) {
    g <- expand.grid(i = seq_len(length(xs) - 1L), j = seq_len(length(ys) - 1L))
    unname(cbind(xs[g$i], xs[g$i + 1L], ys[g$j], ys[g$j + 1L]))
}

# points on every grid line and corner, plus interior points
lattice_points <- function(step = 5, w = 100, h = 50) {
    g <- expand.grid(x = seq(0, w, by = step / 2), y = seq(0, h, by = step / 2))
    data.frame(x = g$x, y = g$y, gene = "g", qv = 30, tid = seq_len(nrow(g)))
}

test_that("points on shared and outer edges count as in R", {
    skip_if_no_reference()
    b <- grid_bounds(seq(0, 100, 25), seq(0, 50, 12.5))
    ref <- expect_counts_match(lattice_points(), b)
    # every point lands in exactly one tile of a contiguous partition
    expect_equal(sum(ref), nrow(lattice_points()))
})

test_that("counts match R with a y flip", {
    skip_if_no_reference()
    b <- grid_bounds(seq(0, 100, 25), seq(-50, 0, 12.5))
    expect_counts_match(lattice_points(), b, flip = TRUE)
})

test_that("a subset of tiles uses the subset's own outer edges", {
    skip_if_no_reference()
    # only the left half of the grid: its right edge (x = 50) becomes inclusive
    b <- grid_bounds(seq(0, 50, 25), seq(0, 50, 12.5))
    expect_counts_match(lattice_points(), b)
})

test_that("irregular quadtree-like tiles count as in R", {
    skip_if_no_reference()
    b <- rbind(c(0, 50, 25, 50), c(50, 100, 25, 50), c(0, 25, 0, 25), c(25, 50, 0, 12.5),
        c(25, 50, 12.5, 25), c(50, 100, 0, 25))
    expect_counts_match(lattice_points(2.5), b)
})

test_that("overlapping tiles count shared points in both, as R does", {
    skip_if_no_reference()
    b <- rbind(c(0, 60, 0, 50), c(40, 100, 0, 50))
    expect_counts_match(lattice_points(), b)
})

test_that("points outside every tile are not counted", {
    skip_if_no_reference()
    b <- grid_bounds(c(20, 40, 60), c(10, 30))
    expect_counts_match(synth_points(3000), b)
})

test_that("float32 coordinates compare against float32-rounded bounds, as arrow does", {
    skip_if_no_reference()
    # coordinates that are exactly float32 values right around bounds that
    # are not representable in float32 (1/3 steps)
    set.seed(3)
    xs <- seq(0, 100, length.out = 7)  # 16.666.. bounds
    ys <- seq(0, 50, length.out = 4)
    near <- c(xs, xs + 1e-5, xs - 1e-5)
    g <- expand.grid(x = near, y = c(ys, ys + 1e-5, ys - 1e-5, runif(5, 0, 50)))
    df <- data.frame(x = g$x, y = g$y, gene = "g", qv = 30, tid = seq_len(nrow(g)))
    tbl <- arrow::arrow_table(df, schema = arrow::schema(x = arrow::float32(),
        y = arrow::float32(), gene = arrow::utf8(), qv = arrow::float64(), tid = arrow::int32()))
    expect_counts_match(write_pq(tbl), grid_bounds(xs, ys))
    expect_counts_match(write_pq(tbl), grid_bounds(xs, -rev(ys)), flip = TRUE)
})

test_that("counts match R with filters applied", {
    skip_if_no_reference()
    b <- grid_bounds(seq(0, 100, 20), seq(0, 50, 10))
    expect_counts_match(synth_points(4000), b,
        filters = list(list(col = "qv", op = ">=", value = 20)))
})

test_that("rows with null coordinates are not counted, as in R", {
    skip_if_no_reference()
    df <- synth_points(1000)
    df$x[c(3, 400, 401)] <- NA
    df$y[c(10, 700)] <- NA
    expect_counts_match(df, grid_bounds(seq(0, 100, 25), seq(0, 50, 25)))
})

test_that("an empty set of bounds returns no counts", {
    skip_if_no_reference()
    f <- write_pq(synth_points(100))
    expect_identical(gi_count_tiles(gi_source(f, "x", "y"), matrix(numeric(0), 0, 4)), numeric(0))
})

test_that("many tiles (grid index path) count as in R", {
    skip_if_no_reference()
    b <- grid_bounds(seq(0, 100, 5), seq(0, 50, 5))   # 200 tiles
    expect_counts_match(synth_points(20000, digits = 0), b, chunk_size = 2000L)
})
