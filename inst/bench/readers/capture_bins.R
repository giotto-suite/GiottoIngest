#!/usr/bin/env Rscript
# Run a bin-points reader job sequentially and keep a copy of the bins parquet
# that .binpoints_store() writes before tiling. The trace is set at the step
# after the bins file is written (`at = 11`); recheck it if GiottoDisk changes.
#   capture_bins.R JOB OUT_DIR
suppressPackageStartupMessages({ library(arrow); library(GiottoDisk) })
options(giotto.has_conda = FALSE, giotto.no_python_warn = TRUE,
        giotto.use_conda = FALSE, progressr.enable = FALSE, tilework.warn_sequential = FALSE)
a <- commandArgs(TRUE); job <- a[[1]]; out <- a[[2]]
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "jobs.R"))
unlink(out, recursive = TRUE); dir.create(out, recursive = TRUE)
keep <- file.path(out, "bins.parquet")
trace(".binpoints_store", where = asNamespace("GiottoDisk"), at = 11, print = FALSE,
      tracer = bquote(file.copy(f, .(keep), overwrite = TRUE)))
future::plan(future::sequential)
b <- GiottoDisk::sourceCreate(file.path(out, "project"), type = "gDirSource")
GiottoDisk::setArtifactDumpDir(b)
t0 <- Sys.time(); jobs[[job]](b)
cat(sprintf("%s: %.1f s, bins parquet %s (%s rows)\n", job,
    as.numeric(Sys.time() - t0, units = "secs"), file.exists(keep),
    format(nrow(open_dataset(keep)), big.mark = ",")))
