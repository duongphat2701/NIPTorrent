#!/usr/bin/env bash
set -euo pipefail

########################################
# SANEFALCON training reference builder
#
# Usage (standalone):
#   ./create_ref.sh /path/to/BAM_DIR /path/to/WORK_DIR
#
# If arguments are omitted, defaults:
#   BAM_DIR = ./BAM
#   WORK_DIR = .
#
# LOG_LEVEL:
#   0 = silent (only fatal errors)
#   1 = normal (step messages)
########################################

BAM_DIR="${1:-$PWD/BAM}"
WORK_DIR="${2:-$PWD}"

SANEFALCON_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PYTHON="${PYTHON:-python}"
SAMTOOLS="${SAMTOOLS:-samtools}"

LOG_LEVEL="${LOG_LEVEL:-1}"

log() {
  if [[ "$LOG_LEVEL" -ge 1 ]]; then
    echo "$@" >&2
  fi
}

log "=== SANEFALCON training pipeline ==="
log "BAM dir : $BAM_DIR"
log "Work dir: $WORK_DIR"
log "Repo dir: $SANEFALCON_DIR"
log

mkdir -p "$WORK_DIR"

##############################
# 1) Extract read start positions
##############################
START_DIR="$WORK_DIR/start_positions"
mkdir -p "$START_DIR"

log "[1/6] Extracting read start positions from BAMs..."
shopt -s nullglob
first_bam_found=false

for bam in "$BAM_DIR"/*.bam; do
  if [[ -e "$bam" ]]; then
    first_bam_found=true
  fi

  base="$(basename "$bam")"
  sample="${base%.bam}"
  sample="${sample%.sort}"   # handle *.sort.bam if present
  out_prefix="$START_DIR/$sample"

  log "  -> Sample: $sample"

  for chr in $(seq 1 22); do
    region="chr${chr}"

    # Forward strand reads (not reverse)
    $SAMTOOLS view "$bam" "$region" -F 16 \
      | $PYTHON "$SANEFALCON_DIR/retro.py" \
      | awk '{print $4}' \
      > "${out_prefix}.${chr}.start.fwd"

    # Reverse strand reads (reverse)
    $SAMTOOLS view "$bam" "$region" -f 16 \
      | $PYTHON "$SANEFALCON_DIR/retro.py" \
      | awk '{print ($4 + length($10) - 1)}' \
      > "${out_prefix}.${chr}.start.rev"
  done
done

if [[ "$first_bam_found" = false ]]; then
  echo "ERROR: No .bam files found in $BAM_DIR" >&2
  exit 1
fi

log "[1/6] Done."
log

##############################
# 2) Build simple train structure (all samples in one 'run')
##############################
log "[2/6] Preparing train directory layout..."

TRAIN_RUN_DIR="$WORK_DIR/train/all/run1"
mkdir -p "$TRAIN_RUN_DIR"

find "$START_DIR" -maxdepth 1 -type f -name "*.start.*" -print0 \
  | while IFS= read -r -d '' f; do
      ln -sf "$f" "$TRAIN_RUN_DIR/"
    done

log "[2/6] Done."
log

##############################
# 3) Merge read starts per chromosome across all training samples
##############################
log "[3/6] Merging read starts per chromosome..."

NUCL_DIR="$WORK_DIR/nucl"
mkdir -p "$NUCL_DIR"

for chr in $(seq 1 22); do
  sort -n -m \
    "$TRAIN_RUN_DIR"/*.${chr}.start.fwd \
    "$TRAIN_RUN_DIR"/*.${chr}.start.rev \
    > "$NUCL_DIR/merge.${chr}"
done

log "[3/6] Done."
log

##############################
# 4) Call nucleosomes for merged training data
##############################
log "[4/6] Calling nucleosomes (nuclDetector.py)..."

for chr in $(seq 1 22); do
  in_file="$NUCL_DIR/merge.${chr}"
  out_file="$NUCL_DIR/nucl_all.${chr}"
  $PYTHON "$SANEFALCON_DIR/nuclDetector.py" "$in_file" "$out_file"
done

log "[4/6] Done."
log

##############################
# 5) Build per-sample nucleosome profiles (getProfile.py)
##############################
log "[5/6] Building nucleosome profiles per sample..."

PROFILE_DIR="$WORK_DIR/profiles"
mkdir -p "$PROFILE_DIR"

for fwd_chr1 in "$START_DIR"/*.1.start.fwd; do
  [[ -e "$fwd_chr1" ]] || continue

  base_nochr="${fwd_chr1%.1.start.fwd}"
  sample_name="$(basename "$base_nochr")"

  log "  -> Sample: $sample_name"

  for chr in $(seq 1 22); do
    nucl_file="$NUCL_DIR/nucl_all.${chr}"
    fwd_file="$START_DIR/${sample_name}.${chr}.start.fwd"
    rev_file="$START_DIR/${sample_name}.${chr}.start.rev"
    out_base="$PROFILE_DIR/${sample_name}.${chr}"

    $PYTHON "$SANEFALCON_DIR/getProfile.py" "$nucl_file" "$fwd_file" 0 "${out_base}.fwd"
    $PYTHON "$SANEFALCON_DIR/getProfile.py" "$nucl_file" "$fwd_file" 1 "${out_base}.ifwd"
    $PYTHON "$SANEFALCON_DIR/getProfile.py" "$nucl_file" "$rev_file" 1 "${out_base}.rev"
    $PYTHON "$SANEFALCON_DIR/getProfile.py" "$nucl_file" "$rev_file" 0 "${out_base}.irev"
  done
done

log "[5/6] Done."
log

##############################
# 6) Combine and fuse profiles → trainNucl.csv
##############################
log "[6/6] Combining upstream/downstream profiles and fusing..."

UPSTREAM_CSV="$WORK_DIR/upstreamProfs.csv"
DOWNSTREAM_CSV="$WORK_DIR/downstreamProfs.csv"
TRAIN_NUCL="$WORK_DIR/trainNucl.csv"

"$SANEFALCON_DIR/combProf.sh"  "$PROFILE_DIR" > "$DOWNSTREAM_CSV"
"$SANEFALCON_DIR/combProfI.sh" "$PROFILE_DIR" > "$UPSTREAM_CSV"

$PYTHON "$SANEFALCON_DIR/fuseProfiles.py" \
  -u "$UPSTREAM_CSV" \
  -d "$DOWNSTREAM_CSV" \
  > "$TRAIN_NUCL"

log
log "=== Done! ==="
log "Training nucleosome profiles written to:"
log "  $TRAIN_NUCL"
log
