#!/usr/bin/env bash
set -euo pipefail

INDIR="${1:?Usage: prepSamples.sh <INDIR> <OUTDIR> [THREADS]}"
OUTDIR="${2:?Usage: prepSamples.sh <INDIR> <OUTDIR> [THREADS]}"
THREADS="${3:-8}"

mkdir -p "$OUTDIR"

# Use system tools from PATH (your conda env etc.)
SAMTOOLS="${SAMTOOLS:-samtools}"
PYTHON_BIN="${PYTHON_BIN:-python}"

# Use retro.py if it exists next to this script (recommended for SANEFALCON)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RETRO_PY="$SCRIPT_DIR/retro.py"
USE_RETRO=0
if [[ -f "$RETRO_PY" ]]; then
  USE_RETRO=1
fi

# Chromosomes to extract
CHROMS=( {1..22} X Y )

# Remove extensions to get clean sample ID:
#   24PC0738_010.bam        -> 24PC0738_010
#   24PC0738_010.sort.bam   -> 24PC0738_010
#   24PC0738_010.sorted.bam -> 24PC0738_010
sample_id_from_path() {
  local p="$1"
  local b
  b="$(basename "$p")"
  b="${b%.bam}"
  b="${b%.cram}"
  b="${b%.sort}"
  b="${b%.sorted}"
  b="${b%.dedup}"
  echo "$b"
}

# Detect if contigs use "chr" prefix (chr1 vs 1)
detect_chr_prefix() {
  local bam="$1"
  # idxstats is fast if BAM has an index; fallback will just assume "chr"
  if "$SAMTOOLS" idxstats "$bam" 2>/dev/null | head -n 5 | cut -f1 | grep -q '^chr1$'; then
    echo "chr"
  else
    echo ""
  fi
}

# Semaphore for background jobs
wait_for_slot() {
  local max_jobs="$1"
  while true; do
    local nj
    nj="$(jobs -rp | wc -l)"
    if (( nj < max_jobs )); then
      break
    fi
    sleep 0.1
  done
}

process_one_bam() {
  local bam="$1"
  local sid
  sid="$(sample_id_from_path "$bam")"

  local prefix
  prefix="$(detect_chr_prefix "$bam")"

  echo "IN: $bam"
  echo "OUT base: $OUTDIR/$sid"
  echo "CHR prefix: '${prefix}'"

  for chr in "${CHROMS[@]}"; do
    local region="${prefix}${chr}"

    # Forward (exclude reverse bit, unmapped, secondary, supplementary, dup, QCfail)
    # Reverse (require reverse bit, exclude unmapped, secondary, supplementary, dup, QCfail)
    # You can tighten mapping criteria by increasing -q (MAPQ).
    wait_for_slot "$THREADS"
    (
      if (( USE_RETRO )); then
        "$SAMTOOLS" view "$bam" "$region" -F 3852 -q 1 \
          | "$PYTHON_BIN" "$RETRO_PY" \
          | awk '{print $4}' > "$OUTDIR/$sid.$chr.start.fwd"
      else
        "$SAMTOOLS" view "$bam" "$region" -F 3852 -q 1 \
          | awk '{print $4}' > "$OUTDIR/$sid.$chr.start.fwd"
      fi
    ) &

    wait_for_slot "$THREADS"
    (
      if (( USE_RETRO )); then
        "$SAMTOOLS" view "$bam" "$region" -f 16 -F 3844 -q 1 \
          | "$PYTHON_BIN" "$RETRO_PY" \
          | awk '{print ($4 + length($10) - 1)}' > "$OUTDIR/$sid.$chr.start.rev"
      else
        "$SAMTOOLS" view "$bam" "$region" -f 16 -F 3844 -q 1 \
          | awk '{print ($4 + length($10) - 1)}' > "$OUTDIR/$sid.$chr.start.rev"
      fi
    ) &
  done

  wait
  echo "DONE sample: $sid"
}

export -f sample_id_from_path detect_chr_prefix wait_for_slot process_one_bam
export OUTDIR THREADS SAMTOOLS PYTHON_BIN RETRO_PY USE_RETRO
export -a CHROMS

# Find BAM/CRAM recursively
mapfile -t BAMS < <(find "$INDIR" -type f \( -name "*.bam" -o -name "*.cram" \) | sort)
if [[ "${#BAMS[@]}" -eq 0 ]]; then
  echo "No BAM/CRAM found under: $INDIR" >&2
  exit 1
fi

for bam in "${BAMS[@]}"; do
  process_one_bam "$bam"
done
