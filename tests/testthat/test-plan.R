# The high-level path: tilework::quadtreePlan() on a gi_source (counts in
# Rust through the tileApply method), then gi_write_points(), against
# GiottoDisk's own storeWrite() of the same source.

test_that("quadtreePlan() on a gi_source is GiottoDisk's plan", {
    skip_if_no_reference()
    df <- synth_points(6000, digits = 4)
    f <- write_pq(df)
    flt <- list(list(col = "qv", op = ">=", value = 20))
    r <- r_point_store(make_qs(f, "x", "y", TRUE, flt), tempfile("r_"), 400, "gene", "x", "y")
    src <- gi_source(f, "x", "y", flip_y = TRUE, filters = flt)
    for (th in c(1L, 4L)) {
        fp <- tilework::quadtreePlan(src, tiles = seed_tiles(gi_extent(src)$ext),
            threshold = 400, threads = th, verbose = FALSE)
        expect_identical(unname(fp$bounds), r$bounds)
        expect_identical(as.numeric(fp$n_records), r$n)
    }
})

test_that("gi_write_points() writes GiottoDisk's store from a quadtree plan", {
    skip_if_no_reference()
    df <- synth_points(6000, digits = 4)
    df$cell_id <- sprintf("c%04d", df$tid %% 97)
    f <- write_pq(df)
    flt <- list(list(col = "qv", op = ">=", value = 20))
    ts <- c("tid", "cell_id", "gene")
    r <- r_point_store(make_qs(f, "x", "y", TRUE, flt, ts), tempfile("r_"), 400, "gene", "x", "y")
    src <- gi_source(f, "x", "y", flip_y = TRUE, filters = flt)
    fp <- tilework::quadtreePlan(src, tiles = seed_tiles(gi_extent(src)$ext), threshold = 400,
        threads = 2L, verbose = FALSE)
    out <- tempfile("rs_")
    w <- gi_write_points(src, fp, out, "gene", to_string = ts, threads = 2L)
    expect_identical(w$counts, r$n[r$n > 0])
    expect_identical(compare_tile_roots(r$root, out, "tid", "row_index"), character(0))
    geo <- function(root) {
        f <- list.files(tile_dirs(root)[1], pattern = "parquet$", full.names = TRUE)[1]
        arrow::ParquetFileReader$create(f)$GetSchema()$metadata$geo
    }
    expect_identical(geo(out), geo(r$root))
})

test_that("empty plan tiles leave gaps in the tile numbering, as in GiottoDisk", {
    skip_if_no_reference()
    df <- synth_points(3000, digits = 4)
    df <- df[df$x < 30 | df$x > 70, ]   # nothing in the middle tile
    f <- write_pq(df)
    qs <- make_qs(f, "x", "y")
    b <- rbind(c(0, 30, 0, 50), c(30, 70, 0, 50), c(70, 100, 0, 50))
    fp <- tilework::freeTilePlan()
    fp$bounds <- b
    fp$n_records <- r_count_tiles(qs, b, "x", "y")
    expect_identical(fp$n_records, c(fp$n_records[1], 0, fp$n_records[3]))
    st <- GiottoDisk::storeCreate(path = tempfile("r_"), type = "parquetGeomTile")
    st <- GiottoDisk::storeWrite(st, qs, tiles = fp, id_col = "gene", sdimx = "x", sdimy = "y",
        type = "points")
    r_root <- gd()$.idpath(st@path, st@uid)
    out <- tempfile("rs_")
    w <- gi_write_points(gi_source(f, "x", "y"), fp, out, "gene")
    expect_identical(w$tile_index, c(1L, 3L))
    expect_identical(basename(tile_dirs(out)), c("tile_index=001", "tile_index=003"))
    expect_identical(compare_tile_roots(r_root, out, "tid", "row_index"), character(0))
})

test_that("gi_write_points() checks the plan", {
    f <- write_pq(synth_points(100))
    expect_error(gi_write_points(gi_source(f, "x", "y"),
        list(bounds = rbind(c(0, 1, 0, 1)), n = c(1, 2)), tempfile(), "gene"), "row count per tile")
})
