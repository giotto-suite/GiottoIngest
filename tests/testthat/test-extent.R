# gi_extent: row count and extent of the filtered (and flipped) source,
# against GiottoDisk's own arrow query (.dplyr_nrow / .dplyr_ext).

r_extent <- function(qs, sx, sy) {
    a <- GiottoDisk::storeRead(qs)
    list(n = gd()$.dplyr_nrow(a),
        ext = tw()$.ext_to_num_vec(gd()$.dplyr_ext(a, sx, sy)))
}

expect_extent_matches <- function(df, sx = "x", sy = "y", flip = FALSE, filters = list(),
                                  threads = 2L, chunk_size = 250L) {
    f <- write_pq(df, chunk_size = chunk_size)
    ref <- r_extent(make_qs(f, sx, sy, flip, filters), sx, sy)
    got <- gi_extent(gi_source(f, sx, sy, flip_y = flip, filters = filters), threads)
    expect_identical(got$n, as.numeric(ref$n))
    expect_identical(got$ext, ref$ext)
}

test_that("extent and n match R for float64 coordinates", {
    skip_if_no_reference()
    expect_extent_matches(synth_points(3000, digits = 6))
})

test_that("extent matches R with a y flip", {
    skip_if_no_reference()
    expect_extent_matches(synth_points(3000, digits = 6), flip = TRUE)
})

test_that("extent matches R for float32 coordinates", {
    skip_if_no_reference()
    df <- synth_points(3000, digits = 6)
    tbl <- arrow::arrow_table(df, schema = arrow::schema(x = arrow::float32(),
        y = arrow::float32(), gene = arrow::utf8(), qv = arrow::float32(), tid = arrow::int32()))
    expect_extent_matches(tbl, flip = TRUE)
})

for (op in c(">=", ">", "<=", "<", "==", "!=")) {
    test_that(sprintf("filter op `%s` keeps the same rows as R", op), {
        skip_if_no_reference()
        df <- synth_points(3000)
        df$qv <- round(df$qv)   # ties on the threshold
        expect_extent_matches(df, filters = list(list(col = "qv", op = op, value = 20)))
    })
}

test_that("a float32 filter column compares like arrow (threshold not exact in f32)", {
    skip_if_no_reference()
    df <- synth_points(3000)
    df$qv <- 0.1 * (0:2999 %% 3)  # 0, 0.1, 0.2 ...: 0.1 is not exact in float32
    tbl <- arrow::arrow_table(df, schema = arrow::schema(x = arrow::float64(),
        y = arrow::float64(), gene = arrow::utf8(), qv = arrow::float32(), tid = arrow::int32()))
    for (op in c(">=", ">", "==")) {
        expect_extent_matches(tbl, filters = list(list(col = "qv", op = op, value = 0.1)))
    }
})

test_that("filter op `in` keeps the same rows as R's %in%", {
    skip_if_no_reference()
    # MERSCOPE's reader keeps a set of FOVs; NA rows are dropped
    df <- synth_points(3000)
    df$fov <- df$tid %% 40L
    df$fov[c(3, 77)] <- NA
    expect_extent_matches(df, filters = list(list(col = "fov", op = "in", value = c(1, 5, 7, 38))))
    tbl <- arrow::arrow_table(df)
    tbl$fov <- tbl$fov$cast(arrow::float32())
    expect_extent_matches(tbl, filters = list(list(col = "fov", op = "in", value = c(0.1, 5))))
})

test_that("several filters combine with AND", {
    skip_if_no_reference()
    expect_extent_matches(synth_points(3000), filters = list(
        list(col = "qv", op = ">=", value = 10), list(col = "tid", op = "<", value = 2000)))
})

test_that("rows with a null filter value are dropped, as in R", {
    skip_if_no_reference()
    df <- synth_points(1000)
    df$qv[c(5, 500, 999)] <- NA
    expect_extent_matches(df, filters = list(list(col = "qv", op = ">=", value = 0)))
})

test_that("null coordinates do not enter the extent", {
    skip_if_no_reference()
    # all data far from 0: a null read as 0.0 would pull xmin / ymin to 0
    df <- synth_points(1000)
    df$x <- df$x + 1000; df$y <- df$y + 1000
    df$x[c(3, 400)] <- NA
    df$y[c(10, 700)] <- NA
    expect_extent_matches(df)
})

test_that("NaN coordinates do not enter the extent", {
    skip_if_no_reference()
    df <- synth_points(1000)
    df$x[7] <- NaN
    f <- write_pq(df)
    got <- gi_extent(gi_source(f, "x", "y"))
    expect_identical(got$ext[1:2], range(df$x, na.rm = TRUE))
})

test_that("several files and many row groups give the same extent", {
    skip_if_no_reference()
    df <- synth_points(4000, digits = 6)
    fs <- c(write_pq(df[1:1500, ], chunk_size = 100L), write_pq(df[1501:4000, ], chunk_size = 333L))
    ref <- r_extent(make_qs(fs, "x", "y"), "x", "y")
    for (th in c(1L, 3L, 8L)) {
        got <- gi_extent(gi_source(fs, "x", "y"), th)
        expect_identical(got$n, as.numeric(ref$n))
        expect_identical(got$ext, ref$ext)
    }
})
