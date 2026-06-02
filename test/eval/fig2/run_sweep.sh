#!/usr/bin/env bash
# Fig 2 Panel A — mixture sweep at full typed-SNP coverage.
#
# 5 mixtures of (h²_SNP, h²_ARG) at total h² = 0.6:
#   (0.6, 0.0) all-SNP
#   (0.45, 0.15) heavy-SNP
#   (0.30, 0.30) balanced
#   (0.15, 0.45) heavy-ARG
#   (0.0, 0.6) all-ARG
#
# 10 seeds × 5 mixtures = 50 cells.  Panmictic, n_train=300, L=300kb.
# Bed coverage = 1.0 (full typing).
#
# Per cell: sim → sargent (full pipeline: BayesR + ARG random).
# Sargent's est.summary.json records h2_snp / h2_arg / h2 — the values
# the plot will use.
#
# Parallelism: xargs -P $N_PAR.  Default 6 (leaves 2 cores for the
# other background sweep).  Resumable: cells with est.summary.json
# already present are skipped.

set -euo pipefail

ROOT=${ROOT:-/Users/acs/git-src/sargent}
SARGENT=${SARGENT:-$ROOT/build/src/progs/sargent}
SIM=${SIM:-$ROOT/scripts/simulate.py}

BASE=${BASE_FIG2:-$HOME/test-sargent/fig2}

N_TRAIN=${N_TRAIN:-300}
N_TEST=${N_TEST:-100}
LENGTH=${LENGTH:-300000}
NE=${NE:-1e4}
CAUSAL_FRAC=${CAUSAL_FRAC:-0.02}
BED_COV=${BED_COV:-1.0}
H2_TOTAL=${H2_TOTAL:-0.6}
N_SEEDS=${N_SEEDS:-10}
N_PAR=${N_PAR:-6}

if [ -n "${PYBIN:-}" ]; then
  PY="$PYBIN"
elif command -v mamba >/dev/null 2>&1 && mamba env list 2>/dev/null | grep -q '^sargent'; then
  PY="mamba run -n sargent python"
else
  PY="python"
fi

mkdir -p "$BASE"

# Mixtures encoded as <name> <h2_snp> <h2_arg>.  Stay synced with the
# legend labels used in plot.py.
MIXTURES=(
  "allSNP 0.60 0.00"
  "heavySNP 0.45 0.15"
  "balanced 0.30 0.30"
  "heavyARG 0.15 0.45"
  "allARG 0.00 0.60"
)

N_MIX=${#MIXTURES[@]}
N_CELLS=$(( N_MIX * N_SEEDS ))

echo "Fig 2 Panel A: mixture sweep"
echo "  base   : $BASE"
echo "  cells  : ${N_MIX} mixtures × ${N_SEEDS} seeds = ${N_CELLS}"
echo "  scale  : n_train=$N_TRAIN n_test=$N_TEST L=$LENGTH Ne=$NE causal_frac=$CAUSAL_FRAC"
echo "  cov    : $BED_COV   total h² = $H2_TOTAL   demog = panmictic"
echo "  parallel: $N_PAR cells at a time"
echo

if [ ! -x "$SARGENT" ]; then
  echo "ERROR: sargent binary not found at $SARGENT" >&2
  exit 1
fi

# Generate (mix_name h2_snp h2_arg seed) tuples for xargs.
JOBLIST="$BASE/.joblist.txt"
> "$JOBLIST"
for mix in "${MIXTURES[@]}"; do
  set -- $mix
  for seed in $(seq 1 $N_SEEDS); do
    echo "$1 $2 $3 $seed" >> "$JOBLIST"
  done
done

# Per-cell worker function (called by xargs).
run_cell() {
  local mix_name="$1"; local h2_snp="$2"; local h2_arg="$3"; local seed="$4"
  local dir="$BASE/${mix_name}_s${seed}"
  mkdir -p "$dir"
  if [ -s "$dir/est.summary.json" ]; then
    echo "[skip] $mix_name s$seed"
    return 0
  fi
  local t0=$(date +%s)
  $PY "$SIM" --demography panmictic \
      --n-indiv "$N_TRAIN" --n-test "$N_TEST" \
      --length "$LENGTH" --Ne "$NE" \
      --h2-snp "$h2_snp" --h2-arg "$h2_arg" \
      --causal-fraction "$CAUSAL_FRAC" \
      --bed-coverage "$BED_COV" \
      --seed "$seed" --out "$dir/case" \
      > "$dir/sim.log" 2>&1 || { echo "  [FAIL sim] $mix_name s$seed"; tail -3 "$dir/sim.log"; return 1; }
  $SARGENT --silent --nprobe 10 --outroot "$dir/est" \
      "$dir/case.phen" "$dir/case.trees" "$dir/case" \
      > "$dir/sargent.log" 2>&1 || { echo "  [FAIL sargent] $mix_name s$seed"; tail -3 "$dir/sargent.log"; return 1; }
  local t1=$(date +%s)
  echo "[done] $mix_name s$seed in $((t1-t0))s"
}
export -f run_cell
export BASE PY SARGENT SIM N_TRAIN N_TEST LENGTH NE CAUSAL_FRAC BED_COV

# Parallel execution.  ETA reported by xargs' total throughput.
t_start=$(date +%s)
cat "$JOBLIST" | xargs -n 4 -P "$N_PAR" bash -c 'run_cell "$@"' _

t_end=$(date +%s)
echo
echo "All cells processed in $(( (t_end - t_start) / 60 )) min."
echo "Run: python test/eval/fig2/plot.py"
