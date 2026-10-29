# Shared harness: synthetic parquet inputs, the GiottoDisk R writer as the
# reference, and the same plan / write done through giottoingest.

`%||%` <- function(a, b) if (is.null(a)) b else a

skip_if_no_reference <- function() {
    skip_if_not_installed("arrow")
    skip_if_not_installed("GiottoDisk")
    skip_if_not_installed("tilework")
    # arrow first: see README (reticulate pyarrow collision)
    suppressPackageStartupMessages({
        library(arrow)
        library(GiottoDisk)
    })
    options(tilework.warn_sequential = FALSE)
    future::plan(future::sequential)
}

gd <- function() asNamespace("GiottoDisk")
tw <- function() asNamespace("tilework")

# Write an arrow Table / data.frame to parquet; small chunk_size gives
# several row groups, so the parallel per-row-group passes are exercised.
write_pq <- function(x, path = tempfile(fileext = ".parquet"), chunk_size = 250L) {
    arrow::write_parquet(x, path, chunk_size = chunk_size)
    path
}

# The lazy query a GiottoDisk reader builds for the same source spec.
read_fun_for <- function(sx, sy, flip = FALSE, filters = list(), to_string = character()) {
    force(sx); force(sy); force(flip); force(filters); force(to_string)
    function(x, schema = NULL) {
        sc <- arrow::schema(arrow::open_dataset(x))
        if (length(to_string)) {
            schema <- do.call(arrow::schema, lapply(names(sc), function(nm) {
                if (nm %in% to_string) arrow::field(nm, arrow::string()) else sc$GetFieldByName(nm)
            }))
        }
        q <- arrow::open_dataset(x, schema = schema)
        for (f in filters) {
            q <- dplyr::filter(q, !!rlang::call2(f$op, as.name(f$col), f$value))
        }
        if (flip) q <- dplyr::mutate(q, !!rlang::sym(sy) := -!!rlang::sym(sy))
        q
    }
}

make_qs <- function(files, sx, sy, flip = FALSE, filters = list(), to_string = character()) {
    as(GiottoDisk::fileStore(path = files,
        read_fun = read_fun_for(sx, sy, flip, filters, to_string)), "queryableStore")
}

# R counts per tile for a set of bounds, via GiottoDisk's contiguous getTile
# over a freeTilePlan -- what each quadtree pass measures.
r_count_tiles <- function(qs, bounds, sx, sy) {
    fp <- tilework::freeTilePlan()
    fp$bounds <- bounds
    as.numeric(unlist(tilework::tileApply(qs, tiles = fp, FUN = nrow,
        sdimx = sx, sdimy = sy, contiguous = TRUE)))
}

# Seed grid, as storeWrite(parquetGeomTileStore, *) builds it.
seed_tiles <- function(e) {
    tiles <- tilework::tilePlan("spatial")
    terra::ext(tiles) <- terra::ext(e[1:2], e[3:4])
    erange <- range(terra::ext(tiles))
    length(tiles) <- round(max(erange) / min(erange)) * 4L
    tiles
}

# tilework::quadtreePlan's control flow with Rust counts (as the bench does).
rs_plan <- function(src, threshold, threads = 2L, max_depth = 10L) {
    ex <- gi_extent(src, threads)
    pending <- seed_tiles(ex$ext)[]
    leaves <- list(); leaf_n <- numeric(0); depth <- 0L
    measure <- function(ts) gi_count_tiles(src, do.call(rbind, lapply(ts, tw()$.ext_to_num_vec)), threads)
    while (length(pending) > 0L && depth < max_depth) {
        depth <- depth + 1L
        m <- measure(pending)
        nxt <- list()
        for (k in seq_along(pending)) {
            if (m[[k]] <= threshold) {
                leaves <- c(leaves, list(pending[[k]])); leaf_n <- c(leaf_n, m[[k]])
            } else {
                nxt <- c(nxt, tw()$.subdivide_extent(pending[[k]]))
            }
        }
        pending <- nxt
    }
    if (length(pending) > 0L) {
        leaves <- c(leaves, pending); leaf_n <- c(leaf_n, measure(pending))
    }
    merged <- tw()$.merge_rect_tiles(do.call(rbind, lapply(leaves, tw()$.ext_to_num_vec)),
        leaf_n, threshold)
    list(ext = ex, bounds = merged$bounds, n = merged$measures)
}

# GiottoDisk's own point tile store for the source.
r_point_store <- function(qs, out, threshold, id_col, sx, sy) {
    st <- GiottoDisk::storeCreate(path = out, type = "parquetGeomTile")
    st <- GiottoDisk::storeWrite(st, qs, threshold = threshold, id_col = id_col,
        sdimx = sx, sdimy = sy, type = "points")
    list(store = st, root = gd()$.idpath(st@path, st@uid),
        bounds = unname(st@tiles@bounds), n = as.numeric(st@tiles@metadata$n_records))
}

# giottoingest's write of the same store, on the given plan.
rs_point_store <- function(src, plan, out, id_col, to_string = character(), drop = character(),
                           geom = TRUE, threads = 2L, spill_rows = 4e6) {
    keep <- which(plan$n > 0)
    b <- plan$bounds
    geo <- gd()$.geoparquet_metadata(geom_col = "geom", geomtype = "points", crs = "",
        extent = c(min(b[, 1]), max(b[, 2]), min(b[, 3]), max(b[, 4])))
    dirs <- basename(gd()$.tilepath(out, seq_along(keep)))
    counts <- gi_write_point_tiles(src, b[keep, , drop = FALSE], dirs, id_col = id_col,
        to_string = to_string, drop = drop, root = out, geo_json = geo,
        spill_rows = spill_rows, geom = geom, threads = threads)
    list(root = out, counts = counts, dirs = dirs)
}

tile_dirs <- function(root) {
    d <- sort(list.dirs(root, recursive = FALSE))
    d[grepl("tile_index=", basename(d))]
}

read_tile <- function(d) as.data.frame(arrow::open_dataset(d))

# Compare two tile roots by key, all columns except `skip`. Returns a
# character vector of differences (empty when identical).
compare_tile_roots <- function(r_root, s_root, key, skip = "row_index") {
    rd <- tile_dirs(r_root); sd <- tile_dirs(s_root)
    if (!identical(basename(rd), basename(sd))) {
        return(sprintf("tile dirs differ: R [%s] vs Rust [%s]",
            paste(basename(rd), collapse = ","), paste(basename(sd), collapse = ",")))
    }
    out <- character()
    for (k in seq_along(rd)) {
        x <- read_tile(rd[k]); y <- read_tile(sd[k])
        x <- x[order(x[[key]]), !names(x) %in% skip, drop = FALSE]
        y <- y[order(y[[key]]), !names(y) %in% skip, drop = FALSE]
        rownames(x) <- rownames(y) <- NULL
        if (!identical(x, y)) {
            what <- if (!identical(names(x), names(y))) {
                sprintf("columns R [%s] vs Rust [%s]", paste(names(x), collapse = ","),
                    paste(names(y), collapse = ","))
            } else if (nrow(x) != nrow(y)) {
                sprintf("rows R %d vs Rust %d", nrow(x), nrow(y))
            } else {
                paste("values differ in", paste(names(x)[!mapply(identical, x, y)], collapse = ", "))
            }
            out <- c(out, paste0(basename(rd[k]), ": ", what))
        }
    }
    out
}

# Morton order check: rows sorted by (morton key over the tile's own
# extent, row_index), exactly what R's order() over .geom_parquet_frame gives.
in_morton_order <- function(tile) {
    gext <- c(min(tile$x_index), max(tile$x_index), min(tile$y_index), max(tile$y_index))
    key <- gd()$.morton_encode(tile$x_index, tile$y_index, gext)
    identical(order(key, tile$row_index), seq_len(nrow(tile)))
}

# Synthetic transcript-like points: float64 x/y, a gene id, a qv filter
# column, and an integer key.
synth_points <- function(n = 4000L, seed = 1L, w = 100, h = 50, digits = 1) {
    set.seed(seed)
    data.frame(
        x = round(stats::runif(n, 0, w), digits),
        y = round(stats::runif(n, 0, h), digits),
        gene = sample(c(letters, LETTERS), n, TRUE),
        qv = round(stats::runif(n, 0, 40), 2),
        tid = seq_len(n)
    )
}
