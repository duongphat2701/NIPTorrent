#!/usr/bin/env bash
set -euo pipefail
INDIR=$(readlink -f "${1:?Usage: mergeAntiSub.sh <train_dir> <subset>}")
SUB="${2:?Usage: mergeAntiSub.sh <train_dir> <subset>}"

mkdir -p "$INDIR/$SUB"

shopt -s nullglob

for CHROM in $(seq 22 -1 1); do
  echo "$CHROM $SUB"

  files=( "$INDIR"/*/merge."$CHROM" )

  # remove the current subset’s own merge file from the list
  filtered=()
  for f in "${files[@]}"; do
    [[ "$f" == "$INDIR/$SUB/merge.$CHROM" ]] && continue
    filtered+=( "$f" )
  done

  if (( ${#filtered[@]} == 0 )); then
    : > "$INDIR/$SUB/anti.$CHROM"
    continue
  fi

  LC_ALL=C sort -n -m "${filtered[@]}" > "$INDIR/$SUB/anti.$CHROM"
done
