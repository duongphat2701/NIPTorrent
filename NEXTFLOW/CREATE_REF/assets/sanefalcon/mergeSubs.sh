#!/usr/bin/env bash
set -euo pipefail
INDIR=$(readlink -f "${1:?Usage: mergeSubs.sh <train_dir>}")

shopt -s nullglob

for i in $(seq 1 22); do
  echo "Working on Chr: $i"

  files=( "$INDIR"/*/merge."$i" )

  if (( ${#files[@]} == 0 )); then
    : > "$INDIR/merge.$i"
    continue
  fi

  LC_ALL=C sort -n -m "${files[@]}" > "$INDIR/merge.$i"
done
