# gi_write_point_tiles, end to end against GiottoDisk's own point tile
# store: same plan, same tile dirs, same rows (by key) in every column, same
# geoparquet metadata; plus the per-tile frame rules (row_index, Morton).

# Run R's writer and the Rust plan + write on the same source; return both.
both_stores <- function(x, sx = "x", sy = "y", flip = FALSE, filters = list(),
                        to_string = character(), id_col = "gene", threshold = 400,
                        threads = 2L, spill_rows = 4e6, geom = TRUE, chunk_size = 250L) {
    f <- if (is.character(x)) x else write_pq(x, chunk_size = chunk_size)
    qs <- make_qs(f, sx, sy, flip, filters, to_string)
    r <- r_point_store(qs, tempfile("r_"), threshold, id_col, sx, sy)
    src <- gi_source(f, sx, sy, flip_y = flip, filters = filters)
    plan <- rs_plan(src, threshold, threads)
    s <- rs_point_store(src, plan, tempfile("rs_"), id_col, to_string = to_string,
        geom = geom, threads = threads, spill_rows = spill_rows)
    list(r = r, plan = plan, s = s, src = src, files = f)
}

expect_same_store <- function(b, key = "tid", skip = "row_index") {
    expect_identical(b$plan$bounds, b$r$bounds)
    expect_identical(b$plan$n, b$r$n)
    expect_identical(compare_tile_roots(b$r$root, b$s$root, key, skip), character(0))
    expect_identical(b$s$counts, b$r$n[b$r$n > 0])
}

geo_meta <- function(root) {
    f <- list.files(tile_dirs(root)[1], pattern = "parquet$", full.names = TRUE)[1]
    arrow::ParquetFileReader$create(f)$GetSchema()$metadata$geo
}

test_that("Xenium-like source (filter, flip, string casts) matches R tile for tile", {
    skip_if_no_reference()
    df <- synth_points(6000, digits = 4)
    df$cell_id <- sprintf("c%04d", df$tid %% 97)
    b <- both_stores(df, flip = TRUE, filters = list(list(col = "qv", op = ">=", value = 20)),
        to_string = c("tid", "cell_id", "gene"))
    expect_gt(nrow(b$plan$bounds), 3L)
    expect_same_store(b)
    expect_identical(geo_meta(b$s$root), geo_meta(b$r$root))
})

test_that("float32 coordinates and narrow integer types match R", {
    skip_if_no_reference()
    df <- synth_points(5000, digits = 6)
    df$u8 <- df$tid %% 200L; df$i16 <- df$tid %% 3000L - 1500L
    df$flag <- df$tid %% 2L == 0L
    tbl <- arrow::arrow_table(df, schema = arrow::schema(x = arrow::float32(), y = arrow::float32(),
        gene = arrow::utf8(), qv = arrow::float32(), tid = arrow::int32(), u8 = arrow::uint8(),
        i16 = arrow::int16(), flag = arrow::boolean()))
    b <- both_stores(write_pq(tbl), flip = TRUE,
        filters = list(list(col = "qv", op = ">", value = 12.3)))
    expect_same_store(b)
})

test_that("uint32 columns get the type R gives them", {
    skip_if_no_reference()
    # R converts uint32 to integer when a tile's values fit in int32 and to
    # double otherwise, so the written type depends on the values
    df <- synth_points(5000, digits = 6)
    df$u32 <- df$tid * 3L
    tbl <- arrow::arrow_table(df)
    tbl$u32 <- tbl$u32$cast(arrow::uint32())
    expect_same_store(both_stores(write_pq(tbl)))
})

test_that("each tile is in Morton order of its own extent, row_index in source order", {
    skip_if_no_reference()
    b <- both_stores(synth_points(6000, digits = 3), flip = TRUE)
    for (d in tile_dirs(b$s$root)) {
        t <- read_tile(d)
        expect_true(in_morton_order(t), label = basename(d))
        expect_setequal(t$row_index, seq_len(nrow(t)))
        # row_index numbers the tile's rows in source order
        expect_identical(t$tid[order(t$row_index)], sort(t$tid))
    }
})

test_that("Morton keys with the high bit set (upper half of y) order as in R", {
    skip_if_no_reference()
    # one tile: R's bitwShiftL wraps to negative int32 for y in the upper half
    df <- synth_points(3000, digits = 6)
    b <- both_stores(df, threshold = Inf)
    t <- read_tile(tile_dirs(b$s$root)[1])
    expect_true(in_morton_order(t))
    r <- read_tile(tile_dirs(b$r$root)[1])
    # with one tile and one sequential scan, R's row order is deterministic
    expect_identical(t$tid, r$tid)
})

test_that("a nearly degenerate source (most x equal) matches R", {
    skip_if_no_reference()
    df <- data.frame(x = c(rep(5, 50), 7), y = c(seq(0, 49), 3), gene = "g", qv = 30,
        tid = 1:51)
    b <- both_stores(df, threshold = Inf)
    expect_same_store(b)
    t <- read_tile(tile_dirs(b$s$root)[1]); r <- read_tile(tile_dirs(b$r$root)[1])
    expect_identical(t$tid, r$tid)
})

test_that("degenerate tiles (zero-width extent, one row) are in R's Morton order", {
    skip_if_no_reference()
    # GiottoDisk's seed grid cannot be built for a zero-width source extent
    # (round(max / min) * 4 with min = 0), so give the tiles directly; the
    # Morton keys are then NA (0 / 0), which R's order() puts last.
    df <- data.frame(x = c(rep(5, 40), 20), y = c(seq(0, 39), 0), gene = "g", qv = 30,
        tid = 1:41)
    src <- gi_source(write_pq(df), "x", "y")
    plan <- list(bounds = rbind(c(0, 10, 0, 40), c(10, 30, 0, 40)), n = c(40, 1))
    s <- rs_point_store(src, plan, tempfile(), "gene")
    for (d in tile_dirs(s$root)) {
        expect_true(in_morton_order(read_tile(d)), label = basename(d))
    }
    expect_identical(s$counts, c(40, 1))
})

test_that("duplicate points and points on tile edges match R", {
    skip_if_no_reference()
    g <- expand.grid(x = seq(0, 100, 2.5), y = seq(0, 50, 2.5))
    df <- data.frame(x = rep(g$x, 2), y = rep(g$y, 2), gene = "g", qv = 30)
    df$tid <- seq_len(nrow(df))
    expect_same_store(both_stores(df, threshold = 600))
})

test_that("a small spill_rows (many spill parts) writes the same store", {
    skip_if_no_reference()
    df <- synth_points(6000, digits = 4)
    f <- write_pq(df, chunk_size = 100L)
    a <- both_stores(f)
    src <- gi_source(f, "x", "y")
    s2 <- rs_point_store(src, a$plan, tempfile("rs_spill_"), "gene", spill_rows = 37)
    expect_identical(compare_tile_roots(a$s$root, s2$root, "tid", skip = character()), character(0))
    # the spill directory is removed after a successful write
    expect_false(any(grepl("^\\.spill-", list.files(s2$root, all.files = TRUE))))
})

test_that("the written store is the same for 1 and 8 threads (deterministic)", {
    skip_if_no_reference()
    f <- write_pq(synth_points(6000, digits = 4), chunk_size = 100L)
    src <- gi_source(f, "x", "y", flip_y = TRUE)
    plan <- rs_plan(src, 400, 1L)
    expect_identical(rs_plan(src, 400, 8L)$bounds, plan$bounds)
    s1 <- rs_point_store(src, plan, tempfile(), "gene", threads = 1L)
    s8 <- rs_point_store(src, plan, tempfile(), "gene", threads = 8L)
    # identical including row_index and row order
    for (k in seq_along(tile_dirs(s1$root))) {
        expect_identical(read_tile(tile_dirs(s8$root)[k]), read_tile(tile_dirs(s1$root)[k]))
    }
})

test_that("geom = FALSE drops only the WKB column and the geo metadata", {
    skip_if_no_reference()
    b <- both_stores(synth_points(4000, digits = 4), geom = FALSE)
    expect_identical(compare_tile_roots(b$r$root, b$s$root, "tid", skip = c("row_index", "geom")),
        character(0))
    t <- read_tile(tile_dirs(b$s$root)[1])
    expect_false("geom" %in% names(t))
    expect_null(geo_meta(b$s$root))
})

test_that("`drop` leaves the columns out", {
    skip_if_no_reference()
    f <- write_pq(synth_points(1000))
    src <- gi_source(f, "x", "y")
    plan <- rs_plan(src, Inf)
    s <- rs_point_store(src, plan, tempfile(), "gene", drop = "qv")
    expect_identical(names(read_tile(tile_dirs(s$root)[1])),
        c("row_index", "x_index", "y_index", "geom", "feat_ID", "tid"))
})

test_that("rows with null coordinates are not written, as in R", {
    skip_if_no_reference()
    df <- synth_points(3000, digits = 4)
    df$x[c(3, 400, 401)] <- NA
    df$y[c(10, 700)] <- NA
    expect_same_store(both_stores(df))
})

# column types R writers meet in vendor files ------------------------------

test_that("int64 columns are written as R writes them", {
    skip_if_no_reference()
    # MERSCOPE detected_transcripts: barcode_id, fov, global_z are int64
    df <- synth_points(3000, digits = 4)
    tbl <- arrow::arrow_table(df, big = df$tid * 10)
    tbl$big <- tbl$big$cast(arrow::int64())
    expect_same_store(both_stores(write_pq(tbl)))
})

test_that("dictionary (factor) columns keep R's values and index type", {
    skip_if_no_reference()
    # R's tiles carry R's own "r" schema metadata, from which arrow in R
    # restores the factor's original levels() order; a parquet dictionary
    # (all a non-R writer can store) only holds levels in order of first
    # appearance. So compare the values as strings and the arrow types.
    df <- synth_points(3000, digits = 4)
    df$gene <- factor(df$gene)
    df$lab <- factor(sprintf("L%03d", df$tid %% 200))  # >= 127 levels: int16 indices
    b <- both_stores(df)
    expect_identical(b$plan$bounds, b$r$bounds)
    as_chr <- function(root) {
        lapply(tile_dirs(root), function(d) {
            t <- read_tile(d)
            t <- t[order(t$tid), c("tid", "feat_ID", "lab")]
            t$feat_ID <- as.character(t$feat_ID); t$lab <- as.character(t$lab)
            rownames(t) <- NULL
            t
        })
    }
    expect_identical(as_chr(b$s$root), as_chr(b$r$root))
    type_of <- function(root, col) {
        f <- list.files(tile_dirs(root)[1], pattern = "parquet$", full.names = TRUE)[1]
        arrow::ParquetFileReader$create(f)$GetSchema()[[col]]$type$ToString()
    }
    for (col in c("feat_ID", "lab")) {
        expect_identical(type_of(b$s$root, col), type_of(b$r$root, col))
    }
})

test_that("binary id columns cast to string match R", {
    skip_if_no_reference()
    # older Xenium transcripts.parquet store feature_name as binary
    df <- synth_points(3000, digits = 4)
    tbl <- arrow::arrow_table(df)
    tbl$gene <- arrow::Array$create(lapply(df$gene, charToRaw), type = arrow::binary())
    expect_same_store(both_stores(write_pq(tbl), to_string = "gene"))
})

test_that("large_utf8 columns are written as R writes them", {
    skip_if_no_reference()
    df <- synth_points(2000, digits = 4)
    tbl <- arrow::arrow_table(df)
    tbl$gene <- tbl$gene$cast(arrow::large_utf8())
    expect_same_store(both_stores(write_pq(tbl)))
})

# several files --------------------------------------------------------------

test_that("several files with the same schema match R", {
    skip_if_no_reference()
    df <- synth_points(6000, digits = 4)
    fs <- c(write_pq(df[1:2500, ], chunk_size = 300L), write_pq(df[2501:6000, ], chunk_size = 700L))
    expect_same_store(both_stores(fs, flip = TRUE))
})

test_that("a later file with float32 coordinates (first file float64) matches R", {
    skip_if_no_reference()
    df <- synth_points(4000, digits = 3)
    t2 <- arrow::arrow_table(df[2001:4000, ])
    t2$x <- t2$x$cast(arrow::float32()); t2$y <- t2$y$cast(arrow::float32())
    fs <- c(write_pq(df[1:2000, ]), write_pq(t2))
    expect_same_store(both_stores(fs))
})

test_that("a later file missing a metadata column matches R", {
    skip_if_no_reference()
    df <- synth_points(4000, digits = 3)
    df$cell_id <- sprintf("c%03d", df$tid %% 50)
    fs <- c(write_pq(df[1:2000, ]), write_pq(df[2001:4000, names(df) != "cell_id"]))
    expect_same_store(both_stores(fs))
})

test_that("geo metadata is in the file's key-value metadata, where DuckDB looks", {
    skip_if_not_installed("duckdb")
    skip_if_not_installed("DBI")
    f <- write_pq(synth_points(500))
    out <- tempfile()
    gi_write_points(gi_source(f, "x", "y"), list(bounds = rbind(c(0, 100, 0, 50)), n = 1),
        out, "gene")
    p <- list.files(out, pattern = "parquet$", recursive = TRUE, full.names = TRUE)[1]
    con <- DBI::dbConnect(duckdb::duckdb())
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
    kv <- DBI::dbGetQuery(con, sprintf(
        "SELECT decode(key) AS key, decode(value) AS value FROM parquet_kv_metadata(%s)",
        DBI::dbQuoteString(con, p)))
    expect_true("geo" %in% kv$key)
    expect_identical(kv$value[kv$key == "geo"],
        arrow::ParquetFileReader$create(p)$GetSchema()$metadata$geo)
})
