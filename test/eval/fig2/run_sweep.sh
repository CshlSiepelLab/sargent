#!/usr/bin/env bash
# Fig 2 sweep (revised 2026-06-02): rare-variant ablation at n=600 with
# extended coverage grid including bed_cov=0.05.  Sparse architecture
# only, panmictic, single population.  Drives both Fig 2 Panel A
# (h² recovery) and Panel B (r²_calib) — same 4 arms, same cells.
#
# Replaces the earlier mixture-sweep driver (mixture story dropped per
# Adam 2026-06-02).
#
# Per cell: sim → sargent → predict.py.  est.predict.json gives the
# four arms' r2_calib + h2 estimates for plotting.
#
# Tunable env (with defaults):
#   BASE        = output root           (~/test-sargent/fig2_n600)
#   N_TRAIN     = train indivs          (600)
#   N_TEST      = test  indivs          (200)
#   LENGTH      = sequence length       (300000)
#   NE          = effective pop size    (1e4)
#   H2_SNP      = target h²_snp         (0.3)
#   H2_ARG      = target h²_arg         (0.3)
#   CAUSAL_FRAC = fraction causal       (0.02)
#   COVERAGES   = bed-coverage values   ("1.0 0.7 0.3 0.1 0.05")
#   ARCH        = ws² architecture      (sparse)
#   BLOCKFRAC   = active fraction       (0.05)
#   CONTRAST    = w² contrast           (1000)
#   N_BLOCKS    = contiguous blocks     (3)
#   N_SEEDS     = number of seeds       (10)
#   N_PAR       = parallel workers      (4)
#
# Resumable: cells with est.predict.json present are skipped.

set -euo pipefail

ROOT=${ROOT:-/Users/acs/git-src/sargent}
SARGENT=${SARGENT:-$ROOT/build/src/progs/sargent}
SIM=${SIM:-$ROOT/scripts/simulate.py}
PRED=${PRED:-$ROOT/test/eval/predict.py}

BASE=${BASE:-$HOME/test-sargent/fig2_n600}

N_TRAIN=${N_TRAIN:-600}
N_TEST=${N_TEST:-200}
LENGTH=${LENGTH:-300000}
NE=${NE:-1e4}
H2_SNP=${H2_SNP:-0.3}
H2_ARG=${H2_ARG:-0.3}
CAUSAL_FRAC=${CAUSAL_FRAC:-0.02}
COVERAGES=${COVERAGES:-"1.0 0.7 0.3 0.1 0.05"}
ARCH=${ARCH:-sparse}
BLOCKFRAC=${BLOCKFRAC:-0.05}
CONTRAST=${CONTRAST:-1000}
N_BLOCKS=${N_BLOCKS:-3}
N_SEEDS=${N_SEEDS:-10}
N_PAR=${N_PAR:-4}

if [ -n "${PYBIN:-}" ]; then
  PY="$PYBIN"
elif command -v mamba >/dev/null 2>&1 && mamba env list 2>/dev/null | grep -q '^sargent'; then
  PY="mamba run -n sargent python"
else
  PY="python"
fi

mkdir -p "$BASE"

N_CV=$(echo $COVERAGES | wc -w | tr -d ' ')
N_CELLS=$(( N_CV * N_SEEDS ))

echo "Fig 2 sweep (n=$N_TRAIN, sparse, coverage axis)"
echo "  base   : $BASE"
echo "  cells  : ${N_CV} cov × ${N_SEEDS} seeds = ${N_CELLS}"
echo "  scale  : n_train=$N_TRAIN n_test=$N_TEST L=$LENGTH Ne=$NE"
echo "  trait  : h2_snp=$H2_SNP h2_arg=$H2_ARG causal_frac=$CAUSAL_FRAC"
echo "  arch   : $ARCH  blocks=$N_BLOCKS frac=$BLOCKFRAC contrast=$CONTRAST"
echo "  parallel: $N_PAR cells at a time"
echo

if [ ! -x "$SARGENT" ]; then
  echo "ERROR: sargent binary not found at $SARGENT" >&2
  exit 1
fi

JOBLIST="$BASE/.joblist.txt"
> "$JOBLIST"
for cov in $COVERAGES; do
  for seed in $(seq 1 $N_SEEDS); do
    echo "$cov $seed" >> "$JOBLIST"
  done
done

run_cell() {
  local cov="$1"; local seed="$2"
  local covtag=$(printf "%g" "$cov")
  local dir="$BASE/a${ARCH}_cv${covtag}_dpanmictic_s${seed}"
  mkdir -p "$dir"
  if [ -s "$dir/est.predict.json" ]; then
    echo "[skip] cov=$covtag s$seed"
    return 0
  fi
  local t0=$(date +%s)
  $PY "$SIM" --demography panmictic \
      --n-indiv "$N_TRAIN" --n-test "$N_TEST" \
      --length "$LENGTH" --Ne "$NE" \
      --h2-snp "$H2_SNP" --h2-arg "$H2_ARG" \
      --causal-fraction "$CAUSAL_FRAC" \
      --bed-coverage "$cov" \
      --ws2-architecture "$ARCH" \
      --ws2-n-blocks "$N_BLOCKS" --ws2-block-frac "$BLOCKFRAC" \
      --ws2-block-contrast "$CONTRAST" \
      --seed "$seed" --out "$dir/case" \
      > "$dir/sim.log" 2>&1 || { echo "  [FAIL sim] cov=$covtag s$seed"; tail -3 "$dir/sim.log"; return 1; }
  $SARGENT --silent --nprobe 10 --outroot "$dir/est" \
      "$dir/case.phen" "$dir/case.trees" "$dir/case" \
      > "$dir/sargent.log" 2>&1 || { echo "  [FAIL sargent] cov=$covtag s$seed"; tail -3 "$dir/sargent.log"; return 1; }
  $PY "$PRED" \
      --case-prefix "$dir/case" --sargent-prefix "$dir/est" \
      > "$dir/predict.log" 2>&1 || { echo "  [FAIL predict] cov=$covtag s$seed"; tail -3 "$dir/predict.log"; return 1; }
  local t1=$(date +%s)
  echo "[done] cov=$covtag s$seed in $((t1-t0))s"
}
export -f run_cell
export BASE PY SARGENT SIM PRED N_TRAIN N_TEST LENGTH NE H2_SNP H2_ARG \
       CAUSAL_FRAC ARCH BLOCKFRAC CONTRAST N_BLOCKS

t_start=$(date +%s)
cat "$JOBLIST" | xargs -n 2 -P "$N_PAR" bash -c 'run_cell "$@"' _

t_end=$(date +%s)
echo
echo "All cells processed in $(( (t_end - t_start) / 60 )) min."
echo "Run: python test/eval/fig2/plot.py --base-rva $BASE"
