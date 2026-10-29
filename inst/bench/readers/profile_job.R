#!/usr/bin/env Rscript
# One profiling job in a fresh process.
#   profile_job.R JOB OUT_DIR MODE
# MODE "seq": sequential plan, Rprof (elapsed) -> stage breakdown.
# MODE "w<N>": N mirai workers, timing only; worker pids are written for the
# footprint sampler (fp_sampler.sh, started by run_jobs.sh).
# GD_USE_INGEST=0/1 picks GiottoDisk's R or Rust (GiottoIngest) write path.
suppressPackageStartupMessages({
  library(arrow)
  library(GiottoDisk)
  library(data.table)
})
options(giotto.has_conda = FALSE, giotto.no_python_warn = TRUE,
        giotto.use_conda = FALSE, progressr.enable = FALSE)
# GD_USE_INGEST=0 / 1: GiottoDisk's R path or the GiottoIngest (Rust) path,
# for readers built on GiottoDisk with the integration; unset: its default
if (nzchar(Sys.getenv("GD_USE_INGEST"))) {
  options(giottodisk.use_ingest = Sys.getenv("GD_USE_INGEST") == "1")
}
a <- commandArgs(TRUE)
job <- a[[1]]; out <- a[[2]]; mode <- a[[3]]
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "jobs.R"))
stopifnot(job %in% names(jobs))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
proj <- file.path(out, "project")
unlink(proj, recursive = TRUE)

if (mode == "seq") {
  future::plan(future::sequential)
} else {
  n <- as.integer(sub("^w", "", mode))
  future::plan(future.mirai::mirai_multisession, workers = n)
  wp <- unique(unlist(mirai::collect_mirai(mirai::mirai_map(1:64, function(i) {
    Sys.sleep(0.05); Sys.getpid()
  }))))
  writeLines(as.character(wp), file.path(out, "workers.txt"))
}

b <- GiottoDisk::sourceCreate(proj, type = "gDirSource")
GiottoDisk::setArtifactDumpDir(b)
prof <- file.path(out, "Rprof.out")
if (mode == "seq") Rprof(prof, interval = 0.01, event = "elapsed")
t0 <- Sys.time()
res <- tryCatch({ jobs[[job]](b); "ok" }, error = function(e) conditionMessage(e))
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (mode == "seq") Rprof(NULL)
sz <- sum(file.size(list.files(proj, recursive = TRUE, full.names = TRUE))) / 1e6
cat(sprintf("JOB %s MODE %s: %s, wall %.1f s, project %.0f MB\n", job, mode, res, secs, sz))
saveRDS(list(job = job, mode = mode, res = res, secs = secs, mb = sz),
        file.path(out, "summary.rds"))

if (mode == "seq" && file.exists(prof)) {
  s <- summaryRprof(prof)
  tot <- s$by.total
  st <- vapply(stages, function(f) {
    k <- paste0('"', f, '"'); if (k %in% rownames(tot)) tot[k, "total.time"] else 0
  }, numeric(1))
  tab <- data.frame(stage = names(stages), sec = round(st, 1),
                    share = round(st / s$sampling.time, 3))
  tab <- tab[tab$sec > 0, ]
  print(tab, row.names = FALSE)
  saveRDS(list(stages = tab, by_self = head(s$by.self, 30), total = s$sampling.time),
          file.path(out, "stages.rds"))
  cat("\ntop self:\n"); print(head(s$by.self, 12))
}
