#!/usr/bin/env bash
set -euo pipefail

# Local version of the original qsub script:
#   python nuclDetector.py INDIR/anti.CHROM  INDIR/nucl_ex3.CHROM

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_NUCDEC="${SCRIPT_NUCDEC:-$SCRIPT_DIR/nuclDetector.py}"
SCRIPT_PYTHON="${SCRIPT_PYTHON:-python}"

INDIR="$(readlink -f "${1:?Usage: nuclDetectorAnti.sh <subset_dir>}")"

for CHROM in $(seq 22 -1 1); do
  IN="$INDIR/anti.$CHROM"
  OUT="$INDIR/nucl_ex3.$CHROM"

  echo "Running: $SCRIPT_PYTHON $SCRIPT_NUCDEC $IN $OUT"

  if [[ ! -s "$IN" ]]; then
    echo "WARNING: missing/empty $IN -> creating empty $OUT"
    : > "$OUT"
    continue
  fi

  "$SCRIPT_PYTHON" "$SCRIPT_NUCDEC" "$IN" "$OUT"
done
