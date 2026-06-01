#!/usr/bin/env bash
# Rare-variant ablation — backs CLAIMS B6 (ARG kinship recovers prediction
# value from variants missing from the typed-SNP panel).
#
# Sweeps --bed-coverage at fixed h²; compares BayesR vs GBLUP vs
# ARG-BLUP-uniform vs SARGENT prediction R² and h² recovery. The
# rare-variant story: as f_typed drops, typed-SNP-only methods (BayesR,
# GBLUP) lose, but ARG-based methods (ARG-BLUP-uniform, SARGENT) keep
# picking up the missing-variant LD via the ARG.
#
# Per the project rule (Adam 2026-06-01), runs both --ws2-architecture
# sparse AND scattered.
#
# Resumable: skips cells with est.predict.json present.
#
# Tunable env:
#   BASE        = output root           (default: ~/test-sargent/rare_variant)
#   N_TRAIN     = train indivs          (300)
#   N_TEST      = test  indivs          (100)
#   LENGTH      = sequence length       (300000)
#   NE          = effective pop size    (1e4)
#   H2_SNP      = target h²_snp         (0.3)
#   H2_ARG      = target h²_arg         (0.3)
#   CAUSAL_FRAC = fraction causal       (0.02)
#   COVERAGES   = bed-coverage values   ("1.0 0.7 0.3 0.1")
#   ARCHS       = ws² architectures     ("sparse scattered")
#                 (REQUIRED both, per project rule.)
#   BLOCKFRAC   = active fraction       (0.05; fixed)
#   CONTRAST    = w² contrast           (1000; fixed)
#   N_BLOCKS    = contiguous blocks     (3; fixed; ignored for scattered)
#   DEMOGS      = demographies          ("panmictic")
#   SEEDS       = seed list             ("1 2 3 4 5")

set -euo pipefail

ROOT=${ROOT:-/Users/acs/git-src/sargent}
SARGENT=${SARGENT:-$ROOT/build/src/progs/sargent}
SIM=${SIM:-$ROOT/scripts/simulate.py}
PRED=${PRED:-$ROOT/test/eval/predict.py}

BASE=${BASE:-$HOME/test-sargent/rare_variant}

N_TRAIN=${N_TRAIN:-300}
N_TEST=${N_TEST:-100}
LENGTH=${LENGTH:-300000}
NE=${NE:-1e4}
H2_SNP=${H2_SNP:-0.3}
H2_ARG=${H2_ARG:-0.3}
CAUSAL_FRAC=${CAUSAL_FRAC:-0.02}
COVERAGES=${COVERAGES:-"1.0 0.7 0.3 0.1"}
ARCHS=${ARCHS:-"sparse scattered"}
BLOCKFRAC=${BLOCKFRAC:-0.05}
CONTRAST=${CONTRAST:-1000}
N_BLOCKS=${N_BLOCKS:-3}
DEMOGS=${DEMOGS:-"panmictic"}
SEEDS=${SEEDS:-"1 2 3 4 5"}

if [ -n "${PYBIN:-}" ]; then
  PY="$PYBIN"
elif command -v mamba >/dev/null 2>&1 && mamba env list 2>/dev/null | grep -q '^sargent'; then
  PY="mamba run -n sargent python"
else
  PY="python"
fi

mkdir -p "$BASE"

N_CV=$(echo $COVERAGES | wc -w | tr -d ' ')
N_AR=$(echo $ARCHS | wc -w | tr -d ' ')
N_DM=$(echo $DEMOGS | wc -w | tr -d ' ')
N_SD=$(echo $SEEDS | wc -w | tr -d ' ')
N_CELLS=$(( N_CV * N_AR * N_DM * N_SD ))

echo "Rare-variant ablation"
echo "  base   : $BASE"
echo "  cells  : ${N_CV} cov × ${N_AR} arch × ${N_DM} demog × ${N_SD} seeds = ${N_CELLS}"
echo "  scale  : n_train=$N_TRAIN n_test=$N_TEST L=$LENGTH"
echo "  trait  : h2_snp=$H2_SNP h2_arg=$H2_ARG causal_frac=$CAUSAL_FRAC"
echo "  arch   : ARCHS={$ARCHS} block_frac=$BLOCKFRAC contrast=$CONTRAST"
echo "  python : $PY"
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
    for cov in $COVERAGES; do
      covtag=$(printf "%g" "$cov")
      for seed in $SEEDS; do
        cell=$((cell + 1))
        dir="$BASE/a${arch}_cv${covtag}_d${demog}_s${seed}"
        mkdir -p "$dir"

        if [ -s "$dir/est.predict.json" ]; then
          echo "[$cell/$N_CELLS] SKIP $dir (already done)"
          skipped=$((skipped + 1))
          continue
        fi

        echo "[$cell/$N_CELLS] arch=$arch cov=$covtag demog=$demog seed=$seed -> $dir"
        cstart=$(date +%s)

        DEMOG_ARGS="--demography $demog"

        $PY "$SIM" $DEMOG_ARGS \
            --n-indiv "$N_TRAIN" --n-test "$N_TEST" \
            --length "$LENGTH" --Ne "$NE" \
            --h2-snp "$H2_SNP" --h2-arg "$H2_ARG" \
            --causal-fraction "$CAUSAL_FRAC" \
            --bed-coverage "$cov" \
            --ws2-architecture "$arch" \
            --ws2-n-blocks "$N_BLOCKS" --ws2-block-frac "$BLOCKFRAC" \
            --ws2-block-contrast "$CONTRAST" \
            --seed "$seed" --out "$dir/case" \
            > "$dir/sim.log" 2>&1 || { echo "  sim FAILED"; tail -5 "$dir/sim.log"; continue; }

        $SARGENT --silent --nprobe 10 --outroot "$dir/est" \
            "$dir/case.phen" "$dir/case.trees" "$dir/case" \
            > "$dir/sargent.log" 2>&1 || { echo "  sargent FAILED"; tail -5 "$dir/sargent.log"; continue; }

        $PY "$PRED" \
            --case-prefix "$dir/case" --sargent-prefix "$dir/est" \
            > "$dir/predict.log" 2>&1 || { echo "  predict FAILED"; tail -5 "$dir/predict.log"; continue; }

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
echo "Total wall: $(( ($(date +%s) - t_start) / 60 )) min"
