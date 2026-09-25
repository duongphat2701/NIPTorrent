#!/usr/bin/env bash
set -euo pipefail
INDIR=$(readlink -f "${1:?Usage: merge.sh <subset_dir>}")

shopt -s nullglob

for i in $(seq 1 22); do
  echo "Working on Chr: $i"

  files=( "$INDIR"/*/*."$i".start.fwd "$INDIR"/*/*."$i".start.rev )

  if (( ${#files[@]} == 0 )); then
    : > "$INDIR/merge.$i"
    continue
  fi

  # numeric merge-sort (input files already sorted individually by prepSamples)
  LC_ALL=C sort -n -m "${files[@]}" > "$INDIR/merge.$i"
done
