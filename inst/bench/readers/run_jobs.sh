#!/bin/bash
# Serialized profiling runs with footprint sampling of R + workers.
#   run_jobs.sh OUT_ROOT MODE JOB [JOB ...]
# MODE: seq (profiled) or wN (N mirai workers). Uses the caller's R_LIBS.
# The footprint sampler (fp_sampler.sh) is macOS only.
P="$(cd "$(dirname "$0")" && pwd)"
FP="$P/fp_sampler.sh"
root=$1; mode=$2; shift 2
for job in "$@"; do
  out="$root/$job.$mode"
  mkdir -p "$out"; rm -f "$out/fp.tsv" "$out/fp.tsv.rss" "$out/workers.txt"
  LANG=en_US.UTF-8 Rscript "$P/profile_job.R" "$job" "$out" "$mode" > "$out/log.txt" 2>&1 &
  rpid=$!
  "$FP" "$rpid" "$out/fp.tsv" "$out/workers.txt" &
  wait "$rpid"; wait
  peak=$(awk 'BEGIN{m=0}{if($2>m)m=$2}END{printf "%.1f", m/1048576}' "$out/fp.tsv")
  echo "$(grep '^JOB' "$out/log.txt" | head -1) | peak footprint ${peak} GB"
done
