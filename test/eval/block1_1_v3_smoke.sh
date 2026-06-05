#!/usr/bin/env bash
# Block 1.1 v3 — smoke test for stdpopsim integration in simulate.py.
#
# Tests three models end-to-end (simulate → sargent → predict.py):
#   OOA:   OutOfAfrica_3G09   train=CEU  test=YRI   pure transfer
#   ADM:   AmericanAdmixture_4B18  train=EUR  test=ADMIX  admixed test
#   ARC:   OutOfAfricaArchaicAdmixture_5R19  train=CEU  test=YRI
#
# Runs at deliberately small n / L so the full pipeline finishes in a few
# minutes on a laptop. Not a stats validation — confirms the plumbing works.
#
# Usage:
#   bash test/eval/block1_1_v3_smoke.sh
#
# Optional env overrides:
#   BASE    output root (default: ~/test-sargent/block1_1_v3_smoke)
#   SEED    random seed (default: 1)

set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}
SARGENT=${SARGENT:-$ROOT/build/src/progs/sargent}
SIM=${SIM:-$ROOT/scripts/simulate.py}
PRED=${PRED:-$ROOT/test/eval/predict.py}

BASE=${BASE:-$HOME/test-sargent/block1_1_v3_smoke}
SEED=${SEED:-1}

# Small n/L for speed
N_TRAIN=100
N_TEST=50
LENGTH=300000
H2_SNP=0.3
H2_ARG=0.3
CAUSAL_FRAC=0.01
TYPED=0.7

if [ -n "${PYBIN:-}" ]; then
  PY="$PYBIN"
elif command -v mamba >/dev/null 2>&1 && mamba env list 2>/dev/null | grep -q '^sargent'; then
  PY="mamba run -n sargent python"
else
  PY="python"
fi

if [ ! -x "$SARGENT" ]; then
  echo "ERROR: sargent binary not found at $SARGENT" >&2
  echo "  Build with: cd $ROOT && cmake -B build && cmake --build build -j" >&2
  exit 1
fi

# Verify stdpopsim is importable
if ! $PY -c "import stdpopsim" 2>/dev/null; then
  echo "ERROR: stdpopsim not installed in the Python env. Run:" >&2
  echo "  pip install stdpopsim" >&2
  exit 1
fi

mkdir -p "$BASE"
echo "Block 1.1 v3 smoke test"
echo "  base   : $BASE"
echo "  scale  : n_train=$N_TRAIN n_test=$N_TEST L=$LENGTH"
echo "  sargent: $SARGENT"
echo "  python : $PY"
echo

PASS=0
FAIL=0
RESULTS=()

run_cell() {
  local tag="$1"; shift    # e.g. "OOA_CEU_YRI"
  local dir="$BASE/$tag"
  mkdir -p "$dir"

  echo "=== $tag ==="

  # 1. simulate
  $PY "$SIM" "$@" \
      --n-indiv "$N_TRAIN" --n-test "$N_TEST" \
      --length "$LENGTH" \
      --h2-snp "$H2_SNP" --h2-arg "$H2_ARG" \
      --causal-fraction "$CAUSAL_FRAC" \
      --bed-coverage "$TYPED" \
      --seed "$SEED" --out "$dir/case" \
      > "$dir/sim.log" 2>&1 \
    || { echo "  FAIL: simulate. See $dir/sim.log"; tail -5 "$dir/sim.log"; FAIL=$((FAIL+1)); return; }

  # 2. sargent
  $SARGENT --silent --nprobe 10 --outroot "$dir/est" \
      "$dir/case.phen" "$dir/case.trees" "$dir/case" \
      > "$dir/sargent.log" 2>&1 \
    || { echo "  FAIL: sargent. See $dir/sargent.log"; tail -5 "$dir/sargent.log"; FAIL=$((FAIL+1)); return; }

  # 3. predict
  $PY "$PRED" \
      --case-prefix "$dir/case" --sargent-prefix "$dir/est" \
      > "$dir/predict.log" 2>&1 \
    || { echo "  FAIL: predict. See $dir/predict.log"; tail -5 "$dir/predict.log"; FAIL=$((FAIL+1)); return; }

  # Report key fields from params.json and predict.json
  local fst
  fst=$($PY -c "import json; p=json.load(open('$dir/case.params.json')); print(p.get('fst_train_test','N/A'))")
  local r2c_arguni
  r2c_arguni=$($PY -c "
import json
d=json.load(open('$dir/est.predict.json'))
a=d.get('arg_blup_uniform',{})
v=a.get('r2_calib',a.get('r2','N/A'))
print(f'{v:.4f}' if isinstance(v,(int,float)) else v)
")
  echo "  PASS  F_ST=$fst  ARG-BLUP-uni r2_calib=$r2c_arguni"
  PASS=$((PASS+1))
  RESULTS+=("$tag: F_ST=$fst  r2_calib=$r2c_arguni")
}

# ---- Cell 1: OutOfAfrica_3G09, CEU train → YRI test (pure transfer, mix=0) ----
run_cell "OOA_3G09_CEU_YRI" \
  --demography stdpopsim \
  --stdpopsim-model OutOfAfrica_3G09 \
  --train-pop CEU --test-pop YRI

# ---- Cell 2: OutOfAfrica_3G09, CEU train → YRI test (5% YRI in training) ----
run_cell "OOA_3G09_CEU_YRI_mix05" \
  --demography stdpopsim \
  --stdpopsim-model OutOfAfrica_3G09 \
  --train-pop CEU --test-pop YRI \
  --mix-fraction-B 0.05

# ---- Cell 3: AmericanAdmixture_4B18, EUR train → ADMIX test ----
run_cell "ADM_4B18_EUR_ADMIX" \
  --demography stdpopsim \
  --stdpopsim-model AmericanAdmixture_4B18 \
  --train-pop EUR --test-pop ADMIX

# ---- Cell 4: OutOfAfricaArchaicAdmixture_5R19, CEU train → YRI test ----
run_cell "ARC_5R19_CEU_YRI" \
  --demography stdpopsim \
  --stdpopsim-model OutOfAfricaArchaicAdmixture_5R19 \
  --train-pop CEU --test-pop YRI

echo
echo "============================================================"
echo "Smoke test complete:  PASS=$PASS  FAIL=$FAIL"
for r in "${RESULTS[@]}"; do echo "  $r"; done
echo "============================================================"

if [ "$FAIL" -gt 0 ]; then exit 1; fi
