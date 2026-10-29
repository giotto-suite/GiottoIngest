#!/bin/sh
# Sum macOS physical footprint (KB) of a process tree once per second -- the
# memory macOS charges a process (Activity Monitor's "Memory"), which unlike
# RSS leaves out freed pages the system can take back.
#   fp_sampler.sh ROOT_PID OUT_TSV [EXTRA_PID_FILE]
# Writes "epoch footprint_kb n_processes" to OUT_TSV and the RSS sum to
# OUT_TSV.rss for comparison.
root=$1; out=$2; extra=$3
tree() {
  all=$1; frontier=$1
  while [ -n "$frontier" ]; do
    next=""
    for p in $frontier; do next="$next $(pgrep -P "$p" | tr '\n' ' ')"; done
    next=$(echo $next)
    [ -n "$next" ] && all="$all $next"
    frontier=$next
  done
  echo $all
}
while kill -0 "$root" 2>/dev/null; do
  roots=$root
  [ -n "$extra" ] && [ -f "$extra" ] && roots="$roots $(tr '\n' ' ' < "$extra")"
  pids=""
  for r in $roots; do pids="$pids $(tree $r)"; done
  pids=$(echo $pids | tr ' ' '\n' | sort -u)
  fp=0; n=0
  for p in $pids; do
    kb=$(/usr/bin/footprint -p "$p" 2>/dev/null | awk '/Footprint:/{
      for (i = 1; i <= NF; i++) if ($i == "Footprint:") { v = $(i+1); u = $(i+2) }
      m = (u == "GB") ? 1048576 : (u == "MB") ? 1024 : (u == "KB") ? 1 : 1/1024
      printf "%d", v * m; exit }')
    [ -n "$kb" ] && { fp=$((fp + kb)); n=$((n + 1)); }
  done
  rss=$(ps -o rss= -p "$(echo $pids | tr ' ' ',')" 2>/dev/null | awk '{s+=$1} END{print s+0}')
  t=$(date +%s)
  printf '%s\t%s\t%s\n' "$t" "$fp" "$n" >> "$out"
  printf '%s\t%s\n' "$t" "$rss" >> "$out.rss"
  sleep 1
done
