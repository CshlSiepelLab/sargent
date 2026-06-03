#!/usr/bin/env bash
# Block 1.1 v2 — SARGENT cross-ancestry "rescue" variants.
#
# For each existing v2 cell ($BASE/sp<split>_mix<m>_typed<t>_s<seed>/case.*),
# re-fit one or more variants (NO re-simulation) and re-run predict.py with
# --sargent-only (BayesR/GBLUP/ARG-uni already exist from the v2 run, so only
# the SARGENT arm K_sg·α is recomputed — ~half the Python kinship cost):
#
#   W  = `sargent --arg-only`  -> est_W  : per-tree ŵ², NO β fixed effect.
#        predict.py's "sargent" arm becomes K_sg·α only (β=0). The "SARGENT-W"
#        rescue: does per-tree weighting help once the β liability is removed?
#   EM = `sargent --em 5`      -> est_EM : joint coordinate descent on (β,ŵ²,σ²).
#        (expensive — ~5× a sequential fit; opt in via VARIANTS.)
#
# Which variants run is set by VARIANTS (default "W"; e.g. "W EM" for both).
# vf-cap stays active (VAR_FIXED_VFCAP, default 0.6). Resumable per variant:
# skips a variant whose est_<V>.predict.json already exists.
#
# Env: BASE, VARIANTS ("W" | "EM" | "W EM"), JOBS (xargs -P, default 12),
#      PYBIN (python w/ msprime), SARGENT_S2_VFCAP.
set -uo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
SARGENT=${SARGENT:-$ROOT/build/src/progs/sargent}
PRED=${PRED:-$ROOT/test/eval/predict.py}
PYBIN=${PYBIN:-python}
BASE=${BASE:-$HOME/test-sargent/block1_1_v2}
VARIANTS=${VARIANTS:-W}
JOBS=${JOBS:-12}

[ -x "$SARGENT" ] || { echo "ERROR: sargent not at $SARGENT (build it)"; exit 1; }
export SARGENT PRED PYBIN VARIANTS
# Per-fit thread cap. Default 1 (so JOBS-way doesn't oversubscribe with the old
# single-threaded binary). The new EMMA binary benefits from threading; set
# OMP_THREADS>1 when it's memory-bound (idle cores) to use them.
OMP_THREADS=${OMP_THREADS:-1}
export OMP_NUM_THREADS=$OMP_THREADS OPENBLAS_NUM_THREADS=$OMP_THREADS \
       MKL_NUM_THREADS=$OMP_THREADS NUMEXPR_NUM_THREADS=$OMP_THREADS \
       VECLIB_MAXIMUM_THREADS=$OMP_THREADS

# sargent flags per variant
sg_flags() { case "$1" in W) echo "--arg-only";; EM) echo "--em 5";; *) echo "";; esac; }

run_cell() {
  local dir="$1" c v t0; c=$(basename "$dir"); t0=$(date +%s)
  [ -s "$dir/case.trees" ] || { echo "SKIP  $c (no case.trees)"; return 0; }
  for v in $VARIANTS; do
    [ -s "$dir/est_${v}.predict.json" ] && continue
    if "$SARGENT" --silent --nprobe 10 $(sg_flags "$v") --outroot "$dir/est_${v}" \
         "$dir/case.phen" "$dir/case.trees" "$dir/case" > "$dir/sargent_${v}.log" 2>&1 \
       && ! grep -qi 'out of memory' "$dir/sargent_${v}.log"; then
      "$PYBIN" "$PRED" --sargent-only --case-prefix "$dir/case" \
         --sargent-prefix "$dir/est_${v}" > "$dir/predict_${v}.log" 2>&1 \
         || echo "WARN  $c predict_${v} rc=$?"
    else
      echo "FAIL  $c sargent $v (err/OOM)"
    fi
  done
  local st=""; for v in $VARIANTS; do st+=$([ -s "$dir/est_${v}.predict.json" ] && echo "$v " || echo "- "); done
  echo "DONE  $c  [$st]  $(( $(date +%s)-t0 ))s"
}
export -f run_cell sg_flags

mapfile -t CELLS < <(ls -d "$BASE"/sp*_mix*_typed*_s* 2>/dev/null)
echo "SARGENT variants v2: ${#CELLS[@]} cells × {$VARIANTS}, JOBS=$JOBS, base=$BASE"
echo "  vfcap=${SARGENT_S2_VFCAP:-0.6(default)}   start $(date '+%F %T')"
printf '%s\n' "${CELLS[@]}" | xargs -P "$JOBS" -I{} bash -c 'run_cell "$@"' _ {}
echo "ALL_VARIANTS_DONE $(date '+%F %T')"
