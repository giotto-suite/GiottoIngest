# gi_route_polygons: each polygon to the tile holding its envelope centroid,
# all its vertex rows in row_index order, against GiottoDisk's own per-tile
# selection (getTile on the written tiles, envelope = TRUE, contiguous).

r_route <- function(qs, bounds, group, sx = "x", sy = "y") {
    fp <- tilework::freeTilePlan()
    fp$bounds <- bounds
    lapply(seq_len(nrow(bounds)), function(k) {
        GiottoDisk::getTile(qs, fp, i = k, contiguous = TRUE, sdimx = sx, sdimy = sy,
            group_col = group, envelope = TRUE)[[1]] |>
            dplyr::arrange(dplyr::across(dplyr::any_of("row_index"))) |>
            dplyr::collect() |>
            as.data.frame()
    })
}

rs_route <- function(src, bounds, group, threads = 2L) {
    out <- tempfile()
    files <- sprintf("t%03d.parquet", seq_len(nrow(bounds)))
    res <- gi_route_polygons(src, bounds, group, out, files, threads)
    tiles <- lapply(files, function(f) {
        p <- file.path(out, f)
        if (file.exists(p)) as.data.frame(arrow::read_parquet(p)) else NULL
    })
    list(res = res, tiles = tiles)
}

expect_routes_match <- function(x, bounds, group = "poly_ID", sx = "x", sy = "y", flip = FALSE) {
    f <- if (is.character(x)) x else write_pq(x, chunk_size = 200L)
    ref <- r_route(make_qs(f, sx, sy, flip), bounds, group, sx, sy)
    for (th in c(1L, 4L)) {
        got <- rs_route(gi_source(f, sx, sy, flip_y = flip), bounds, group, th)
        for (k in seq_along(ref)) {
            r <- ref[[k]]; s <- got$tiles[[k]]
            if (nrow(r) == 0L) {
                expect_true(is.null(s) || nrow(s) == 0L, label = sprintf("tile %d empty", k))
            } else {
                rownames(r) <- NULL
                expect_identical(s, r, label = sprintf("tile %d (threads %d)", k, th))
            }
        }
        expect_identical(got$res$counts, vapply(ref, nrow, numeric(1)))
    }
    invisible(got)
}

# random star-shaped polygons, vertex rows shuffled (source order != row_index)
synth_polys <- function(np = 300, seed = 2, w = 100, h = 50) {
    set.seed(seed)
    ctr <- data.frame(id = sprintf("p%03d", seq_len(np)), cx = stats::runif(np, 0, w),
        cy = stats::runif(np, 0, h))
    v <- do.call(rbind, lapply(seq_len(np), function(k) {
        a <- sort(stats::runif(6, 0, 2 * pi)); r <- stats::runif(6, 0.5, 2)
        data.frame(poly_ID = ctr$id[k], x = ctr$cx[k] + r * cos(a), y = ctr$cy[k] + r * sin(a))
    }))
    v$row_index <- seq_len(nrow(v))
    v[sample(nrow(v)), ]
}

quad_bounds <- rbind(c(-5, 50, -5, 25), c(50, 105, -5, 25), c(-5, 50, 25, 55), c(50, 105, 25, 55))

test_that("random polygons route as R selects them", {
    skip_if_no_reference()
    got <- expect_routes_match(synth_polys(), quad_bounds)
    expect_identical(got$res$groups, 300)
    expect_identical(got$res$unassigned, 0)
})

test_that("routing matches R with a y flip", {
    skip_if_no_reference()
    b <- quad_bounds; b[, 3:4] <- -b[, 4:3]
    expect_routes_match(synth_polys(), b, flip = TRUE)
})

test_that("envelope centroids exactly on shared and outer edges route as in R", {
    skip_if_no_reference()
    sq <- function(id, x0, x1, y0, y1) data.frame(poly_ID = id, x = c(x0, x1, x1, x0),
        y = c(y0, y0, y1, y1))
    v <- rbind(
        sq("on_x50", 48, 52, 10, 12),     # centroid x = 50: shared vertical edge
        sq("on_y25", 10, 12, 23, 27),     # centroid y = 25: shared horizontal edge
        sq("on_both", 49, 51, 24, 26),    # centroid on the inner corner
        sq("on_right", 103, 107, 5, 7),   # centroid on the outer right edge (x = 105)
        sq("on_bottom", 30, 32, -7, -3),  # centroid on the outer bottom edge (y = -5)
        sq("on_top", 30, 32, 53, 57),     # centroid on the outer top edge (y = 55)
        sq("on_left", -7, -3, 30, 32))    # centroid on the outer left edge (x = -5)
    v$row_index <- seq_len(nrow(v))
    expect_routes_match(v, quad_bounds)
})

test_that("float32 coordinates: centroid sum in float32 routes as in R", {
    skip_if_no_reference()
    # (xmin + xmax) / 2 lands just either side of x = 50 depending on whether
    # the sum is done in float32 or double
    set.seed(5)
    eps <- c(-3e-6, -1e-6, 0, 1e-6, 3e-6)
    v <- do.call(rbind, lapply(seq_along(eps), function(k) {
        lo <- 49.3 + eps[k]; hi <- 50.7 + eps[k] * 2
        data.frame(poly_ID = sprintf("e%d", k), x = c(lo, hi, hi, lo), y = c(10, 10, 11, 11) + k)
    }))
    v <- rbind(v, synth_polys(100, seed = 6)[, c("poly_ID", "x", "y")])
    v$row_index <- seq_len(nrow(v))
    tbl <- arrow::arrow_table(v, schema = arrow::schema(poly_ID = arrow::utf8(),
        x = arrow::float32(), y = arrow::float32(), row_index = arrow::int32()))
    expect_routes_match(write_pq(tbl), quad_bounds)
})

test_that("integer polygon ids route as in R", {
    skip_if_no_reference()
    v <- synth_polys(200, seed = 7)
    v$poly_ID <- as.integer(sub("p", "", v$poly_ID))
    expect_routes_match(v, quad_bounds)
})

test_that("rows with a null polygon id are handled as in R", {
    skip_if_no_reference()
    v <- synth_polys(100, seed = 8)
    v$poly_ID[v$poly_ID == "p010"] <- NA
    expect_routes_match(v, quad_bounds)
})

test_that("without a row_index column, tiles hold R's rows, in source order", {
    skip_if_no_reference()
    # R's row order is not reproducible here (its envelope selection is a
    # join), so compare the rows as sets; Rust keeps source order.
    v <- synth_polys(100, seed = 9)
    v <- v[order(v$row_index), names(v) != "row_index"]
    f <- write_pq(v, chunk_size = 200L)
    ref <- r_route(make_qs(f, "x", "y"), quad_bounds, "poly_ID")
    got <- rs_route(gi_source(f, "x", "y"), quad_bounds, "poly_ID")
    key <- function(d) paste(d$poly_ID, d$x, d$y)
    for (k in seq_along(ref)) {
        s <- key(got$tiles[[k]])
        expect_setequal(s, key(ref[[k]]))
        expect_identical(s, key(v)[key(v) %in% s])
    }
})

test_that("a polygon whose centroid is in no written tile is reported, as R loses it too", {
    skip_if_no_reference()
    # a C shape: vertices in the left and right tiles only; its envelope
    # centroid is in the middle tile, which holds no vertices and so is not
    # among the tiles written (the plan counts vertices)
    # C spans x 1..29 (vertices at x 1, 3 and 29); centroid x = 15
    v <- data.frame(poly_ID = rep(c("C", "R"), c(8, 4)),
        x = c(1, 29, 29, 3, 3, 29, 29, 1, 21, 28, 28, 21),
        y = c(1, 1, 3, 3, 7, 7, 9, 9, 1, 1, 9, 9))
    v$row_index <- seq_len(nrow(v))
    b <- rbind(c(0, 10, 0, 10), c(20, 30, 0, 10))  # middle tile 10..20 omitted
    got <- expect_routes_match(v, b)
    expect_identical(got$res$unassigned, 1)
})

test_that("no polygon is lost when written tiles are the plan's non-empty tiles", {
    skip_if_no_reference()
    # README finding 3, through the real plan: every polygon must reach a tile
    v <- data.frame(poly_ID = "C", x = c(0, 100, 100, 2, 2, 100, 100, 0),
        y = c(0, 0, 2, 2, 48, 48, 50, 50))
    fill <- synth_polys(60, seed = 10)
    fill <- fill[fill$x < 30 | fill$x > 70, ]
    v <- rbind(v, fill[, c("poly_ID", "x", "y")])
    v$row_index <- seq_len(nrow(v))
    src <- gi_source(write_pq(v), "x", "y")
    plan <- rs_plan(src, threshold = 40)
    keep <- which(plan$n > 0)
    got <- rs_route(src, plan$bounds[keep, , drop = FALSE], "poly_ID")
    expect_identical(got$res$unassigned, 0)
})
