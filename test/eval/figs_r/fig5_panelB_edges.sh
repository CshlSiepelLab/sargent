#!/usr/bin/env bash
# Fig 5 Panel B (per-edge version) — generate cells under per-edge age-stratified
# architectures (Block 2.7.1).  Three architectures × 3 seeds × small n.
#
# Per cell: simulate → sargent (--arg-only) → done.  edge_h2 runs in a
# separate post-pass via test/eval/figs_r/preprocess_fig5_edge.py.
#
# Small n + short L so edge_h2 V-build is fast (≤ 1 min/cell).

set -euo pipefail

ROOT=${ROOT:-/Users/acs/git-src/sargent}
SARGENT=${SARGENT:-$ROOT/build/src/progs/sargent}
SIM=${SIM:-$ROOT/scripts/simulate.py}
BASE=${BASE:-$HOME/test-sargent/fig5_panelB_edges}

N_TRAIN=${N_TRAIN:-300}
N_TEST=${N_TEST:-0}
LENGTH=${LENGTH:-200000}      # 200 kb → ~1500 trees → fast edge_h2
NE=${NE:-1e4}
H2_ARG=${H2_ARG:-0.85}
SEEDS=${SEEDS:-"1 2 3"}
N_PAR=${N_PAR:-3}

PY="mamba run -n sargent python"
mkdir -p "$BASE"

ARCHS=("age_recent_edges" "age_middle_edges" "age_ancient_edges")
JOBLIST="$BASE/.joblist.txt"
> "$JOBLIST"
for arch in "${ARCHS[@]}"; do
  for s in $SEEDS; do
    echo "$arch $s" >> "$JOBLIST"
  done
done

echo "Fig 5 Panel B (per-edge) cells: $(wc -l <"$JOBLIST" | tr -d ' ')"
echo "  base : $BASE"
echo "  n=$N_TRAIN  L=$LENGTH  h²_arg=$H2_ARG"
echo

run_cell() {
  local arch="$1"; local seed="$2"
  local dir="$BASE/${arch}_s${seed}"
  mkdir -p "$dir"
  if [ -s "$dir/est.w2.tsv" ]; then
    echo "[skip] $arch s$seed"
    return 0
  fi
  local t0=$(date +%s)
  $PY "$SIM" --demography panmictic \
      --n-indiv "$N_TRAIN" --n-test "$N_TEST" \
      --length "$LENGTH" --Ne "$NE" \
      --h2-snp 0.0 --h2-arg "$H2_ARG" \
      --ws2-architecture "$arch" \
      --seed "$seed" --out "$dir/case" \
      > "$dir/sim.log" 2>&1 || { echo "[FAIL sim] $arch s$seed"; tail -3 "$dir/sim.log"; return 1; }
  $SARGENT --silent --nprobe 10 --arg-only \
      --outroot "$dir/est" \
      "$dir/case.phen" "$dir/case.trees" "$dir/case" \
      > "$dir/sargent.log" 2>&1 || { echo "[FAIL sargent] $arch s$seed"; tail -3 "$dir/sargent.log"; return 1; }
  local t1=$(date +%s)
  echo "[done] $arch s$seed in $((t1-t0))s"
}
export -f run_cell
export BASE PY SIM SARGENT N_TRAIN N_TEST LENGTH NE H2_ARG

t_start=$(date +%s)
cat "$JOBLIST" | xargs -n 2 -P "$N_PAR" bash -c 'run_cell "$@"' _
t_end=$(date +%s)

echo
echo "Sweep done in $(( (t_end - t_start) / 60 )) min."
echo "Next: python test/eval/figs_r/preprocess_fig5_edge.py"
