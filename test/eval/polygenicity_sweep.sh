#!/usr/bin/env bash
# Block 2.5.2 — controlled-polygenicity sweep driver.
#
# Sweeps ws²-block-fraction (= effective polygenicity) at fixed contrast +
# block count + demography, to test whether SARGENT's M_e,ARG / variance-
# capture curve can separate the three regimes:
#
#   BLOCKFRAC = 0.005  -> oligogenic    (~ 5 active trees out of 1000)
#   BLOCKFRAC = 0.05   -> moderate      (~ 50 active trees; Phase 10 / Block 1.2 baseline)
#   BLOCKFRAC = 0.5    -> highly polygenic (~ 500 active trees)
#
# Truth M_e,ARG ≈ BLOCKFRAC × n_trees ≈ {5, 50, 500} at S~1000 (closed form).
#
# Decision rule:
#   If estimated M_e,ARG separates the three regimes monotonically (even
#     compressed): qualitative architecture story is the manuscript answer.
#   If estimated M_e,ARG is flat across regimes: estimator inadequate, the
#     fused-lasso (Phase 11) revival is justified.
#
# Resumable: skips any cell whose est.summary.json + est.w2.tsv exist.
#
# Tunable via env:
#   BASE        = output root           (default: ~/test-sargent/block2_5_2)
#   N_TRAIN     = train indivs          (300)
#   N_TEST      = test  indivs          (100)
#   LENGTH      = sequence length       (300000)
#   NE          = effective pop size    (1e4)
#   H2_ARG      = target h²_arg         (0.5)
#   BLOCKFRACS  = polygenicity sweep    ("0.005 0.05 0.5")
#   CONTRAST    = block:floor w² ratio  (1000; fixed)
#   N_BLOCKS    = contiguous blocks     (3; fixed)
#   DEMOGS      = demographies          ("admixed panmictic")
#   ARCHS       = ws² architectures     ("sparse scattered")
#                 sparse = n_blocks contiguous active runs (default fused-lasso-friendly)
#                 scattered = active trees individually distributed (fused-lasso
#                             robustness check; ~biology-realistic distribution)
#                 Per Adam (2026-06-01): scattered MUST be included in every
#                 polygenicity / estimator-amplitude validation going forward.
#   SPLIT, ADMIXT, ALPHA = admixed params (4000, 20, 0.5)
#   SEEDS       = seed list             ("1 2 3 4 5")

set -euo pipefail

ROOT=${ROOT:-/Users/acs/git-src/sargent}
SARGENT=${SARGENT:-$ROOT/build/src/progs/sargent}
SIM=${SIM:-$ROOT/scripts/simulate.py}

BASE=${BASE:-$HOME/test-sargent/block2_5_2}

N_TRAIN=${N_TRAIN:-300}
N_TEST=${N_TEST:-100}
LENGTH=${LENGTH:-300000}
NE=${NE:-1e4}
H2_ARG=${H2_ARG:-0.5}
BLOCKFRACS=${BLOCKFRACS:-"0.005 0.05 0.5"}
CONTRAST=${CONTRAST:-1000}
N_BLOCKS=${N_BLOCKS:-3}
DEMOGS=${DEMOGS:-"admixed panmictic"}
ARCHS=${ARCHS:-"sparse scattered"}
SPLIT=${SPLIT:-4000}
ADMIXT=${ADMIXT:-20}
ALPHA=${ALPHA:-0.5}
SEEDS=${SEEDS:-"1 2 3 4 5"}

if [ -n "${PYBIN:-}" ]; then
  PY="$PYBIN"
elif command -v mamba >/dev/null 2>&1 && mamba env list 2>/dev/null | grep -q '^sargent'; then
  PY="mamba run -n sargent python"
else
  PY="python"
fi

mkdir -p "$BASE"

N_BF=$(echo $BLOCKFRACS | wc -w | tr -d ' ')
N_DEM=$(echo $DEMOGS | wc -w | tr -d ' ')
N_AR=$(echo $ARCHS | wc -w | tr -d ' ')
N_SD=$(echo $SEEDS | wc -w | tr -d ' ')
N_CELLS=$(( N_BF * N_DEM * N_AR * N_SD ))

echo "Block 2.5.2: polygenicity sweep"
echo "  base   : $BASE"
echo "  cells  : ${N_BF} block_frac × ${N_DEM} demog × ${N_AR} arch × ${N_SD} seeds = ${N_CELLS}"
echo "  scale  : n_train=$N_TRAIN n_test=$N_TEST L=$LENGTH Ne=$NE h2_arg=$H2_ARG"
echo "  arch   : ARCHS={$ARCHS} contrast=$CONTRAST n_blocks=$N_BLOCKS"
echo "  python : $PY"
echo "  sargent: $SARGENT"
echo

if [ ! -x "$SARGENT" ]; then
  echo "ERROR: sargent binary not found at $SARGENT" >&2
  exit 1
fi

t_start=$(date +%s)
cell=0
skipped=0

for arch in $ARCHS; do
  for demog in $DEMOGS; do
    for bf in $BLOCKFRACS; do
      bftag=$(printf "%g" "$bf")
      for seed in $SEEDS; do
        cell=$((cell + 1))
        # Cell-name convention includes arch so 'sparse' and 'scattered'
        # cells live side-by-side in the same BASE without collision.
        dir="$BASE/a${arch}_bf${bftag}_d${demog}_s${seed}"
        mkdir -p "$dir"

        if [ -s "$dir/est.w2.tsv" ] && [ -s "$dir/est.summary.json" ]; then
          echo "[$cell/$N_CELLS] SKIP $dir (already done)"
          skipped=$((skipped + 1))
          continue
        fi

        echo "[$cell/$N_CELLS] arch=$arch bf=$bftag demog=$demog seed=$seed -> $dir"
        cstart=$(date +%s)

        DEMOG_ARGS="--demography $demog"
        if [ "$demog" = "admixed" ]; then
          DEMOG_ARGS="$DEMOG_ARGS --admix-split-time $SPLIT --admix-time $ADMIXT --admix-alpha $ALPHA"
        fi

        # 1. simulate (pure random; --ws2-architecture varied per ARCHS)
        $PY "$SIM" $DEMOG_ARGS \
            --n-indiv "$N_TRAIN" --n-test "$N_TEST" \
            --length "$LENGTH" --Ne "$NE" \
            --h2-snp 0 --h2-arg "$H2_ARG" \
            --ws2-architecture "$arch" \
            --ws2-n-blocks "$N_BLOCKS" --ws2-block-frac "$bf" \
            --ws2-block-contrast "$CONTRAST" \
            --seed "$seed" --out "$dir/case" \
            > "$dir/sim.log" 2>&1 || { echo "  sim FAILED"; tail -5 "$dir/sim.log"; continue; }

        # 2. sargent --arg-only (Stage 2)
        $SARGENT --silent --nprobe 10 --arg-only --outroot "$dir/est" \
            "$dir/case.phen" "$dir/case.trees" "$dir/case" \
            > "$dir/sargent.log" 2>&1 || { echo "  sargent FAILED"; tail -5 "$dir/sargent.log"; continue; }

        cend=$(date +%s)
        celapsed=$(( cend - cstart ))
        total_elapsed=$(( cend - t_start ))
        done_count=$(( cell - skipped ))
        avg=$(( done_count > 0 ? total_elapsed / done_count : 0 ))
        eta=$(( (N_CELLS - cell) * avg ))
        printf "    done in %ds  (avg %ds/cell, ETA %dm)\n" "$celapsed" "$avg" "$((eta / 60))"
      done
    done
  done
done

echo
echo "All cells processed ($cell total, $skipped skipped)."
echo
echo "Run arch_metrics + plot:"
echo "  cd $BASE && \\"
echo "    $PY $ROOT/test/eval/arch_metrics.py --batch . --no-branches --out arch_b252.json && \\"
echo "    $PY $ROOT/test/eval/plot_arch.py arch_b252.json --group-key ws2_block_frac --out-prefix arch_"
echo

t_end=$(date +%s)
echo "Total wall time: $(( (t_end - t_start) / 60 )) min"
