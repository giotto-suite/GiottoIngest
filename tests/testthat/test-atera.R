# The zarr input: Atera v1 transcripts.zarr.zip grid cells, against
# GiottoDisk's own R decoder (.atera_tx_cell_worker) on the same plan.
# Needs the 10x tiny Atera dataset under $GI_DATA (default ~/Documents/Datasets).

atera_tiny <- function() {
    d <- Sys.getenv("GI_DATA", file.path("~", "Documents", "Datasets"))
    z <- path.expand(file.path(d, "Atera", "tiny_atera_dataset_outs", "transcripts.zarr.zip"))
    if (!file.exists(z)) skip("tiny Atera dataset not found (set GI_DATA)")
    z
}

atera_parts <- function(root) {
    d <- list.dirs(root, recursive = TRUE)
    d <- sort(d[grepl("^tile_index=", basename(d))])
    names(d) <- basename(d)
    d
}

test_that("Atera zarr grid cells decode and write as GiottoDisk's R decoder", {
    skip_if_no_reference()
    z <- atera_tiny()
    ns <- gd()
    src <- ns$.zarr_open(z)
    n_total <- as.numeric(ns$.zarr_attrs(src, "")$number_rnas)
    plan <- ns$.atera_tx_tile_plan(src, ns$.auto_threshold(n_total, type = "points"))
    lk <- ns$.atera_tx_lookups(src)
    ns$.zarr_close(src)
    cells <- names(sort(attr(plan, "rows"), decreasing = TRUE))
    ords <- match(cells, unique(plan$tile)) - 1L
    s <- vapply(cells, function(t) plan$s[plan$tile == t][1L], numeric(1))
    ti <- lapply(cells, function(t) plan$tile_index[plan$tile == t])

    # whole buffers (one part per tile), then small ones: many parts per
    # tile, so row_index continues across parts
    for (rows in c(4194304, 97)) {
        r_root <- tempfile("r_"); rs_root <- tempfile("rs_")
        r_res <- lapply(cells, function(t) ns$.atera_tx_cell_worker(z, t, plan, r_root, "uid0",
            qv_threshold = 20, flip = TRUE, window_rows = rows, buffer_rows = rows,
            compression = ns$.parquet_compression()))
        rs_res <- gi_atera_tx_cells(z, cells, ords, s, ti, lk$gene_names, lk$cw_gene,
            root = rs_root, qv_threshold = 20, flip = TRUE, window_rows = rows,
            buffer_rows = rows, threads = 2L)
        for (k in seq_along(cells)) {
            expect_identical(rs_res[[k]]$tile_index, as.integer(r_res[[k]]$tile_index))
            expect_identical(rs_res[[k]]$counts, as.numeric(r_res[[k]]$counts))
            expect_identical(unname(rs_res[[k]]$bounds), unname(r_res[[k]]$bounds))
        }
        rp <- atera_parts(r_root); sp <- atera_parts(rs_root)
        expect_identical(names(sp), names(rp))
        for (d in names(rp)) {
            rf <- sort(list.files(rp[[d]], pattern = "parquet$"))
            expect_identical(sort(list.files(sp[[d]], pattern = "parquet$")), rf, label = d)
            for (f in rf) {
                a <- arrow::read_parquet(file.path(rp[[d]], f))
                b <- arrow::read_parquet(file.path(sp[[d]], f))
                # same rows in the same (Morton) order, every column
                expect_identical(as.data.frame(b), as.data.frame(a), label = paste(d, f))
                geo <- function(p) arrow::ParquetFileReader$create(p)$GetSchema()$metadata$geo
                expect_identical(geo(file.path(sp[[d]], f)), geo(file.path(rp[[d]], f)))
            }
        }
    }
})
