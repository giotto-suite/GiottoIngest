#!/usr/bin/env Rscript
# Run a polygon reader job sequentially and keep the polygon tile writer's
# inputs (the parquetStore of vertex rows and every writer argument), taken
# when it plans its tiles.
#   capture_polys.R JOB OUT_DIR
suppressPackageStartupMessages({ library(arrow); library(GiottoDisk) })
options(giotto.has_conda = FALSE, giotto.no_python_warn = TRUE,
        giotto.use_conda = FALSE, progressr.enable = FALSE, tilework.warn_sequential = FALSE)
a <- commandArgs(TRUE); job <- a[[1]]; out <- a[[2]]
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "jobs.R"))
unlink(out, recursive = TRUE); dir.create(out, recursive = TRUE)
keep <- file.path(out, "poly_args.rds")
trace(".tile_plan_quadtree", where = asNamespace("GiottoDisk"), print = FALSE,
      tracer = bquote({
        fr <- Filter(function(e) exists("group_col", envir = e, inherits = FALSE) &&
                       exists("split_geom", envir = e, inherits = FALSE), sys.frames())
        if (length(fr)) {
          pf <- fr[[length(fr)]]
          nm <- c("data", "store", "threshold", "id_col", "sdimx", "sdimy", "part_col",
                  "group_col", "contiguous", "split_geom", "split_geom_fmt",
                  "split_geom_sourcename", "flip_vertical", "write_param")
          saveRDS(mget(nm, envir = pf), .(keep))
        }
      }))
future::plan(future::sequential)
b <- GiottoDisk::sourceCreate(file.path(out, "project"), type = "gDirSource")
GiottoDisk::setArtifactDumpDir(b)
t0 <- Sys.time(); jobs[[job]](b)
cat(sprintf("%s: %.1f s, captured %s\n", job, as.numeric(Sys.time() - t0, units = "secs"),
            file.exists(keep)))
