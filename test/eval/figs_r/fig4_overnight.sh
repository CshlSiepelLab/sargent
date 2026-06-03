#!/usr/bin/env bash
# Fig 4 overnight driver: re-run all cells used by Panels A, B, C with
# `--ws2-l2 0.5` and --arg-only.  This locks in the recommended
# localization estimator (TV smoothing over per-tree weights) for the
# main-text figure.
#
# Picked from smoke test on c1000_dpanmictic_s2 (2026-06-02):
#   default: AUC 0.837 (--arg-only)
#   l2=0.5:  AUC 0.858  ← +0.021, ratio 20× active/floor
#   l2=1.0:  AUC 0.857  (saturated)
#
# Cells:
#   - Block 1.2 v1 (30 cells: c{10,100,1000} × d{admixed,panmictic} × 5 seeds)
#   - Phase 10 admixed (10 cells: 5 sparse + 5 diffuse)
#   - Panel A (c1000_dpanmictic_s2 with l2_0.5) reuses the smoke-test
#     output at ~/test-sargent/fig4_smoke/, no new run.
#
# Output: ~/test-sargent/fig4_l2/<cell-name>/{est.w2.tsv, est.summary.json, sargent.log}
# Driver is resumable: cells whose est.w2.tsv exists are skipped.

set -euo pipefail

ROOT=${ROOT:-/Users/acs/git-src/sargent}
SARGENT=${SARGENT:-$ROOT/build/src/progs/sargent}
BASE=${BASE:-$HOME/test-sargent/fig4_l2}
N_PAR=${N_PAR:-4}
WS2_L2=${WS2_L2:-0.5}

mkdir -p "$BASE"

if [ ! -x "$SARGENT" ]; then
  echo "ERROR: sargent binary not found at $SARGENT" >&2
  exit 1
fi

JOBLIST="$BASE/.joblist.txt"
> "$JOBLIST"

# Phase 10 cells
for src in $HOME/test-sargent/phase10_admix/sparse_s? $HOME/test-sargent/phase10_admix/diffuse_s?; do
  [ -d "$src" ] || continue
  name=$(basename "$src")
  echo "phase10  $src  $name" >> "$JOBLIST"
done

# Block 1.2 v1 cells
for src in $HOME/test-sargent/block1_2_v1/c*_d*_s*; do
  [ -d "$src" ] || continue
  name=$(basename "$src")
  echo "block12  $src  $name" >> "$JOBLIST"
done

N_JOBS=$(wc -l < "$JOBLIST" | tr -d ' ')
echo "Fig 4 overnight L2=$WS2_L2 sweep"
echo "  base   : $BASE"
echo "  cells  : $N_JOBS"
echo "  parallel: $N_PAR cells at a time"
echo

run_cell() {
  local group="$1"; local src="$2"; local name="$3"
  local dir="$BASE/${group}_${name}"
  mkdir -p "$dir"
  if [ -s "$dir/est.w2.tsv" ] && [ -s "$dir/est.summary.json" ]; then
    echo "[skip] $name"
    return 0
  fi
  local t0=$(date +%s)
  $SARGENT --silent --nprobe 10 --arg-only --ws2-l2 "$WS2_L2" \
      --outroot "$dir/est" \
      "$src/case.phen" "$src/case.trees" "$src/case" \
      > "$dir/sargent.log" 2>&1 || { echo "  [FAIL sargent] $name"; tail -3 "$dir/sargent.log"; return 1; }
  local t1=$(date +%s)
  echo "[done] $name in $((t1-t0))s"
}
export -f run_cell
export BASE SARGENT WS2_L2

t_start=$(date +%s)
cat "$JOBLIST" | xargs -n 3 -P "$N_PAR" bash -c 'run_cell "$@"' _

t_end=$(date +%s)
echo
echo "Fig 4 overnight done in $(( (t_end - t_start) / 60 )) min."
echo "Output root: $BASE"
