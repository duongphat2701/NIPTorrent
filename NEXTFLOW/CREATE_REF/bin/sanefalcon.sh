#!/usr/bin/env bash
set -euo pipefail

log(){ echo "[$(date '+%F %T')] $*" >&2; }

usage() {
  cat >&2 <<'EOF'
Usage:
  sanefalcon.sh --sanefalcon-dir DIR --bamdir DIR [--outdir DIR] [--threads INT] [--ff_ref_tsv FILE]
                (--all | --starts | --split | --merge | --anti | --nucl | --profiles | --combine | --train)

Notes:
  - OUTDIR can be provided via --outdir or environment variable OUTDIR.
  - Output structure:
      $OUTDIR/ref_build/{logs,starts_raw,split,profiles,final}
EOF
}

need_cmd(){ command -v "$1" >/dev/null 2>&1 || { echo "Missing command: $1" >&2; exit 1; }; }

# ----------------------------
# Defaults / inputs
# ----------------------------
SANEFALCON_DIR=""
BAM_DIR=""
OUTDIR="${OUTDIR:-}"
THREADS="${THREADS:-16}"
FF_REF_TSV="${FF_REF_TSV:-}"

STEP=""

# split tuning (README used ~25 per batch; adjust if you want)
TARGET_PER_SUBSET="${TARGET_PER_SUBSET:-25}"

# tools
SAMTOOLS="${SAMTOOLS:-samtools}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

# ----------------------------
# Parse args
# ----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sanefalcon-dir) SANEFALCON_DIR="${2:?}"; shift 2;;
    --bamdir)         BAM_DIR="${2:?}"; shift 2;;
    --outdir)         OUTDIR="${2:?}"; shift 2;;
    --threads)        THREADS="${2:?}"; shift 2;;
    --ff_ref_tsv)     FF_REF_TSV="${2:?}"; shift 2;;

    --all|--starts|--split|--merge|--anti|--nucl|--profiles|--combine|--train)
      STEP="$1"; shift;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

[[ -n "$SANEFALCON_DIR" ]] || { echo "ERROR: --sanefalcon-dir is required" >&2; usage; exit 1; }
[[ -n "$BAM_DIR" ]]        || { echo "ERROR: --bamdir is required" >&2; usage; exit 1; }
[[ -n "$OUTDIR" ]]         || { echo "ERROR: --outdir or env OUTDIR is required" >&2; usage; exit 1; }
[[ -n "$STEP" ]]           || { echo "ERROR: one of --all/--starts/... is required" >&2; usage; exit 1; }

SANEFALCON_DIR="$(readlink -f "$SANEFALCON_DIR")"
BAM_DIR="$(readlink -f "$BAM_DIR")"
OUTDIR="$(readlink -f "$OUTDIR")"

need_cmd "$SAMTOOLS"
need_cmd "$PYTHON_BIN"
need_cmd find
need_cmd awk
need_cmd sort

WORK="$OUTDIR/ref_build"
LOGDIR="$WORK/logs"
STARTS_RAW="$WORK/starts_raw"
SPLIT_DIR="$WORK/split"
TRAIN_DIR="$SPLIT_DIR/train"
PROF_DIR="$WORK/profiles"
FINAL_DIR="$WORK/final"

mkdir -p "$LOGDIR" "$FINAL_DIR"

# Always run in SANEFALCON_DIR per your request
cd_sanefalcon() { cd "$SANEFALCON_DIR"; }

# Semaphore
wait_for_slot() {
  local max_jobs="$1"
  while true; do
    local nj
    nj="$(jobs -rp | wc -l | tr -d ' ')"
    if (( nj < max_jobs )); then break; fi
    sleep 0.1
  done
}

# Sample ID extraction must match prepSamples.sh logic
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

# -----------------------------------
# Step 1: Extract read starts
# -----------------------------------
step_starts() {
  log "STEP starts: prepSamples.sh"
  mkdir -p "$STARTS_RAW"
  cd_sanefalcon

  # Uses SANEFALCON_DIR/prepSamples.sh exactly
  # prepSamples.sh <INDIR> <OUTDIR> [THREADS]
  PYTHON_BIN="$PYTHON_BIN" SAMTOOLS="$SAMTOOLS" \
    bash "./prepSamples.sh" "$BAM_DIR" "$STARTS_RAW" "$THREADS" \
    >"$LOGDIR/01_prepSamples.log" 2>&1

  log "  Output: $STARTS_RAW"
}

# -----------------------------------
# Step 2: Split into subsets (symlinks only)
# -----------------------------------
step_split() {
  log "STEP split: create train/<subset>/<run>/ symlink layout"
  mkdir -p "$TRAIN_DIR"

  # 2.1 gather BAM/CRAM and map sample -> run (first path element under BAM_DIR)
  mapfile -t BAMLIST < <(find "$BAM_DIR" -type f \( -name "*.bam" -o -name "*.cram" \) | sort)
  if [[ "${#BAMLIST[@]}" -eq 0 ]]; then
    echo "ERROR: No BAM/CRAM under $BAM_DIR" >&2
    exit 1
  fi

  # Create a TSV mapping: sample_id <tab> run_id
  MAP_TSV="$WORK/sample_to_run.tsv"
  : > "$MAP_TSV"
  for bam in "${BAMLIST[@]}"; do
    sid="$(sample_id_from_path "$bam")"
    rel="${bam#"$BAM_DIR"/}"
    if [[ "$rel" == "$bam" ]]; then
      run="run1"
    else
      run="$(echo "$rel" | cut -d'/' -f1)"
      [[ -n "$run" ]] || run="run1"
    fi
    printf "%s\t%s\n" "$sid" "$run" >> "$MAP_TSV"
  done

  # 2.2 count samples per run
  RUN_COUNTS="$WORK/run_counts.tsv"
  awk -F'\t' '{c[$2]++} END{for(r in c) print r "\t" c[r]}' "$MAP_TSV" \
    | sort -k2,2nr -k1,1 > "$RUN_COUNTS"

  # 2.3 assign runs to subsets (greedy; keeps whole runs together)
  # subsets: a b c ... (as many as needed)
  ASSIGN_TSV="$WORK/run_to_subset.tsv"
  : > "$ASSIGN_TSV"

  # Create subset names list
  # (supports up to 26 subsets: a..z)
  alphabet=( {a..z} )

  subset_idx=0
  subset_samples=0
  while IFS=$'\t' read -r run cnt; do
    if (( subset_samples > 0 && subset_samples + cnt > TARGET_PER_SUBSET )); then
      subset_idx=$((subset_idx+1))
      subset_samples=0
    fi
    if (( subset_idx >= ${#alphabet[@]} )); then
      echo "ERROR: Need more than 26 subsets; increase TARGET_PER_SUBSET or customize." >&2
      exit 1
    fi
    subset="${alphabet[$subset_idx]}"
    subset_samples=$((subset_samples+cnt))
    printf "%s\t%s\n" "$run" "$subset" >> "$ASSIGN_TSV"
  done < "$RUN_COUNTS"

  # 2.4 make dirs and symlink start files into train/<subset>/<run>/
  # We symlink *all* chr start files for each sample.
  while IFS=$'\t' read -r sid run; do
    subset="$(awk -F'\t' -v r="$run" '$1==r{print $2; exit}' "$ASSIGN_TSV")"
    [[ -n "$subset" ]] || subset="a"
    out_run_dir="$TRAIN_DIR/$subset/$run"
    mkdir -p "$out_run_dir"

    # link starts (whatever chromosomes exist)
    for f in "$STARTS_RAW/${sid}."*.start.*; do
      [[ -e "$f" ]] || continue
      ln -sf "$f" "$out_run_dir/"
    done
  done < "$MAP_TSV"

  log "  Split dir: $TRAIN_DIR"
  log "  Run->subset: $ASSIGN_TSV"
}

# -----------------------------------
# Step 3: Merge per subset + all
# -----------------------------------
step_merge() {
  log "STEP merge: merge.sh per subset + mergeSubs.sh"

  cd_sanefalcon

  # per subset
  for subset_dir in "$TRAIN_DIR"/*; do
    [[ -d "$subset_dir" ]] || continue
    sub="$(basename "$subset_dir")"
    log "  merge subset: $sub"
    bash "./merge.sh" "$subset_dir" >"$LOGDIR/03_merge_${sub}.log" 2>&1
  done

  # all subsets
  bash "./mergeSubs.sh" "$TRAIN_DIR" >"$LOGDIR/03_mergeSubs.log" 2>&1
}

# -----------------------------------
# Step 4: Anti-merge per subset
# -----------------------------------
step_anti() {
  log "STEP anti: mergeAntiSub.sh for each subset"

  cd_sanefalcon

  for subset_dir in "$TRAIN_DIR"/*; do
    [[ -d "$subset_dir" ]] || continue
    sub="$(basename "$subset_dir")"
    log "  anti-merge subset: $sub"
    bash "./mergeAntiSub.sh" "$TRAIN_DIR" "$sub" >"$LOGDIR/04_mergeAnti_${sub}.log" 2>&1
  done
}

# -----------------------------------
# Step 5: Nucleosome tracks
# -----------------------------------
step_nucl() {
  log "STEP nucl: nuclDetectorAnti.sh per subset + nuclDetector.sh for all"

  cd_sanefalcon

  for subset_dir in "$TRAIN_DIR"/*; do
    [[ -d "$subset_dir" ]] || continue
    sub="$(basename "$subset_dir")"
    log "  nucl anti subset: $sub"
    bash "./nuclDetectorAnti.sh" "$subset_dir" >"$LOGDIR/05_nuclAnti_${sub}.log" 2>&1
  done

  log "  nucl for all training (merge.* -> nucl_ex4.*)"
  bash "./nuclDetector.sh" "$TRAIN_DIR" >"$LOGDIR/05_nuclAll.log" 2>&1
}

# -----------------------------------
# Step 6: Profiles (README 4× calls)
# -----------------------------------

run_profiles_one_sample() {
  local subset="$1"
  local sample_base="$2"   # full path without ".1.start.fwd"
  local sample_name
  sample_name="$(basename "$sample_base")"

  local subset_dir="$TRAIN_DIR/$subset"
  local out_subset_dir="$PROF_DIR/$subset"
  mkdir -p "$out_subset_dir"

  cd_sanefalcon

  for chr in $(seq 1 22); do
    # Prefer anti-track for this subset; fallback to global track (all training)
    local nucl_anti="$subset_dir/nucl_ex3.$chr"
    local nucl_all="$TRAIN_DIR/nucl_ex4.$chr"
    local nucl=""

    if [[ -s "$nucl_anti" ]]; then
      nucl="$nucl_anti"
    elif [[ -s "$nucl_all" ]]; then
      nucl="$nucl_all"
    else
      continue
    fi

    local fwd="${sample_base}.${chr}.start.fwd"
    local rev="${sample_base}.${chr}.start.rev"
    local out_base="$out_subset_dir/${sample_name}.${chr}"

    # Skip if missing
    [[ -s "$fwd" ]] || continue
    [[ -s "$rev" ]] || continue

    # Exactly as README (4 runs)
    "$PYTHON_BIN" "./getProfile.py" "$nucl" "$fwd" 0 "${out_base}.fwd"  >/dev/null
    "$PYTHON_BIN" "./getProfile.py" "$nucl" "$fwd" 1 "${out_base}.ifwd" >/dev/null
    "$PYTHON_BIN" "./getProfile.py" "$nucl" "$rev" 1 "${out_base}.rev"  >/dev/null
    "$PYTHON_BIN" "./getProfile.py" "$nucl" "$rev" 0 "${out_base}.irev" >/dev/null
  done
}

step_profiles() {
  log "STEP profiles: getProfile.py for each sample (anti-track per subset; fallback to global)"
  mkdir -p "$PROF_DIR"

  for subset_dir in "$TRAIN_DIR"/*; do
    [[ -d "$subset_dir" ]] || continue

    sub="$(basename "$subset_dir")"
    log "  subset: $sub"

    # Do NOT restrict to -type f; symlinks are common in split/
    mapfile -t sample_chr1 < <(find "$subset_dir" -name "*.1.start.fwd" | sort)

    log "  Found ${#sample_chr1[@]} samples in subset $sub"
    if [[ "${#sample_chr1[@]}" -eq 0 ]]; then
      echo "ERROR: No *.1.start.fwd found under $subset_dir" >&2
      exit 1
    fi

    for f in "${sample_chr1[@]}"; do
      base="${f%.1.start.fwd}"
      wait_for_slot "$THREADS"
      (
        run_profiles_one_sample "$sub" "$base" \
          >"$LOGDIR/06_profiles_${sub}_$(basename "$base").log" 2>&1
      ) &
    done
    wait
  done

  log "  Profiles output: $PROF_DIR"
}

# -----------------------------------
# Step 7: Combine profiles -> trainNucl.csv
# -----------------------------------
step_combine() {
  log "STEP combine: combProf.sh + combProfI.sh + fuseProfiles.py"

  cd_sanefalcon
  mkdir -p "$FINAL_DIR"

  UP_CSV="$FINAL_DIR/upstreamProfs.csv"
  DOWN_CSV="$FINAL_DIR/downstreamProfs.csv"
  TRAIN_NUCL="$FINAL_DIR/trainNucl.csv"

  # Follow README naming (combProf -> upstream, combProfI -> downstream)
  bash "./combProf.sh"  "$PROF_DIR" > "$UP_CSV"
  bash "./combProfI.sh" "$PROF_DIR" > "$DOWN_CSV"

  # fuseProfiles.py must be Py3-compatible if PYTHON_BIN=python3
  "$PYTHON_BIN" "./fuseProfiles.py" -u "$UP_CSV" -d "$DOWN_CSV" > "$TRAIN_NUCL"

  log "  Output: $TRAIN_NUCL"
}

# -----------------------------------
# Step 8: Train model (predictor.py)
# -----------------------------------
step_train() {
  log "STEP train: predictor.py (requires FF_REF_TSV)"

  if [[ -z "$FF_REF_TSV" ]]; then
    echo "ERROR: --ff_ref_tsv is required for --train/--all" >&2
    exit 1
  fi
  FF_REF_TSV="$(readlink -f "$FF_REF_TSV")"
  [[ -s "$FF_REF_TSV" ]] || { echo "ERROR: FF_REF_TSV not found or empty: $FF_REF_TSV" >&2; exit 1; }

  TRAIN_NUCL="$FINAL_DIR/trainNucl.csv"
  [[ -s "$TRAIN_NUCL" ]] || { echo "ERROR: Missing trainNucl.csv (run --combine first): $TRAIN_NUCL" >&2; exit 1; }

  mkdir -p "$FINAL_DIR" "$LOGDIR"

  # Build a predictor-compatible reference file (NO HEADER!)
  REF_FIXED="$FINAL_DIR/trainRef.sanefalcon.txt"
  awk '
    BEGIN{OFS=" "}
    {
      gsub(/\r$/,"");
      if ($0 ~ /^[[:space:]]*$/) next;
      if ($0 ~ /^[[:space:]]*#/) next;
      if ($0 ~ /(sample_id|ff_reference|sex)/) next;

      if (NF < 2) next;

      # path sample ff [sex]
      if ($1 ~ /\// && NF >= 3) {
        sex = (NF>=4 ? $4 : "Male");
        ff  = $3;
        if (ff ~ /^-?[0-9]+(\.[0-9]+)?$/) print $2, ff, sex;
        next;
      }

      # sample ff sex
      if (NF >= 3) {
        ff = $2;
        if (ff ~ /^-?[0-9]+(\.[0-9]+)?$/) print $1, ff, $3;
        next;
      }

      # sample ff
      ff = $2;
      if (ff ~ /^-?[0-9]+(\.[0-9]+)?$/) print $1, ff, "Male";
    }
  ' "$FF_REF_TSV" > "$REF_FIXED"

  # Hard guarantee: strip any lingering header/comment lines
  sed -i '/ff_reference/d; /sample_id/d; /^#/d; /^[[:space:]]*$/d' "$REF_FIXED"

  [[ -s "$REF_FIXED" ]] || { echo "ERROR: REF_FIXED empty after normalization: $REF_FIXED" >&2; exit 1; }

  # Show first lines for debugging
  log "  Using REF_FIXED=$REF_FIXED"
  head -n 3 "$REF_FIXED" >&2 || true

  cd_sanefalcon
  OUTBASE="$FINAL_DIR/sanefalcon_model"

  # Your predictor.py variant requires extra args (testNucl, testRef, startCol, endCol).
  # We use train as test to satisfy argv[4]/argv[5] and slice from col 1 onward.
  set +e
  "$PYTHON_BIN" "./predictor.py" \
    "$TRAIN_NUCL" \
    "$REF_FIXED" \
    "$OUTBASE" \
    "$TRAIN_NUCL" \
    "$REF_FIXED" \
    1 1000000 \
    >"$LOGDIR/08_train.log" 2>&1
  rc=$?
  set -e

  if (( rc != 0 )); then
    echo "ERROR: predictor.py failed (see $LOGDIR/08_train.log)." >&2
    tail -n 80 "$LOGDIR/08_train.log" >&2 || true
    exit $rc
  fi

  log "  Model outputs basename: $OUTBASE"
}

# ----------------------------
# Run steps
# ----------------------------
log "SANEFALCON_DIR=$SANEFALCON_DIR"
log "BAM_DIR=$BAM_DIR"
log "OUTDIR=$OUTDIR"
log "THREADS=$THREADS"
log "FF_REF_TSV=${FF_REF_TSV:-<none>}"
log "STEP=$STEP"
log "WORK=$WORK"

case "$STEP" in
  --starts)   step_starts ;;
  --split)    step_split ;;
  --merge)    step_merge ;;
  --anti)     step_anti ;;
  --nucl)     step_nucl ;;
  --profiles) step_profiles ;;
  --combine)  step_combine ;;
  --train)    step_train ;;
  --all)
    step_starts
    step_split
    step_merge
    step_anti
    step_nucl
    step_profiles
    step_combine
    step_train
    ;;
  *) echo "ERROR: unknown STEP: $STEP" >&2; exit 1;;
esac

log "DONE."