#!/usr/bin/env bash
# Phase 10 oracle diagnostic driver.
#
# For each ws² architecture {diffuse, sparse} × seed:
#   1. simulate a PURE-RANDOM dataset (--h2-snp 0) with a known per-tree ws²
#   2. run `sargent --arg-only` to estimate ŵ²
#   3. score oracle / uniform / estimated ARG-BLUP prediction + ws² recovery
# Then aggregate into a markdown table.
#
# Override scale via env, e.g.:
#   N_TRAIN=200 N_TEST=50 LENGTH=300000 SEEDS="42" ARCHS="diffuse sparse" \
#     bash test/eval/oracle_ws2.sh
set -euo pipefail

ROOT=/Users/acs/git-src/sargent
SARGENT=$ROOT/build/src/progs/sargent
SIM=$ROOT/scripts/simulate.py
ORACLE=$ROOT/test/eval/oracle_ws2.py
BASE=${BASE:-/Users/acs/test-sargent/phase10}

N_TRAIN=${N_TRAIN:-300}
N_TEST=${N_TEST:-100}
LENGTH=${LENGTH:-300000}
NE=${NE:-1e4}
H2_ARG=${H2_ARG:-0.5}
SEEDS=${SEEDS:-"42"}
ARCHS=${ARCHS:-"diffuse sparse"}
# sparse architecture parameters
NBLOCKS=${NBLOCKS:-3}
BLOCKFRAC=${BLOCKFRAC:-0.05}
CONTRAST=${CONTRAST:-1000.0}
# demography: panmictic | admixed (admixed = optimistic regime for ws²)
DEMOG=${DEMOG:-panmictic}
SPLIT=${SPLIT:-4000}
ADMIXT=${ADMIXT:-20}
ALPHA=${ALPHA:-0.5}

mkdir -p "$BASE"
echo "Phase 10: BASE=$BASE  n=$N_TRAIN+$N_TEST  L=$LENGTH  Ne=$NE  h2_arg=$H2_ARG"
echo "          archs=[$ARCHS]  seeds=[$SEEDS]"

for arch in $ARCHS; do
  for seed in $SEEDS; do
    dir="$BASE/${arch}_s${seed}"
    mkdir -p "$dir"
    echo "=== arch=$arch seed=$seed -> $dir ==="

    mamba run -n sargent python "$SIM" \
        --n-indiv "$N_TRAIN" --n-test "$N_TEST" \
        --length "$LENGTH" --Ne "$NE" \
        --demography "$DEMOG" --admix-split-time "$SPLIT" \
        --admix-time "$ADMIXT" --admix-alpha "$ALPHA" \
        --h2-snp 0 --h2-arg "$H2_ARG" \
        --ws2-architecture "$arch" \
        --ws2-n-blocks "$NBLOCKS" --ws2-block-frac "$BLOCKFRAC" \
        --ws2-block-contrast "$CONTRAST" \
        --seed "$seed" --out "$dir/case" \
        > "$dir/sim.log" 2>&1

    # Stage-2-only estimate (pure random; no SNP fixed effects)
    $SARGENT --silent --nprobe 10 --arg-only --outroot "$dir/est" \
        "$dir/case.phen" "$dir/case.trees" "$dir/case" \
        > "$dir/sargent.log" 2>&1

    mamba run -n sargent python "$ORACLE" \
        --case-prefix "$dir/case" --sargent-prefix "$dir/est" \
        > "$dir/oracle.log" 2>&1 || { echo "  oracle FAILED — see $dir/oracle.log"; tail -5 "$dir/oracle.log"; }
  done
done

# --- aggregate ---------------------------------------------------------------
mamba run -n sargent python - "$BASE" <<'PY' | tee "$BASE/phase10_table.md"
import json, os, sys, glob
base = sys.argv[1]
rows = []
for jf in sorted(glob.glob(os.path.join(base, "*", "est.oracle.json"))):
    d = json.load(open(jf))
    rows.append(d)

def g(d, k, sub="r2"):
    return d.get(k, {}).get(sub, float("nan")) if isinstance(d.get(k), dict) else float("nan")

print("# Phase 10 — oracle diagnostic\n")
print("| arch | seed-dir | R²_uni | R²_oracle | R²_est | **oracle−uni** | est−uni | loc-ceil AUC | est AUC |")
print("|---|---|---|---|---|---|---|---|---|")
for d in rows:
    case = d.get("case","")
    seeddir = os.path.basename(os.path.dirname(case)) if case else "?"
    rec = d.get("ws2_recovery", {})
    print(f"| {d.get('ws2_architecture')} | {seeddir} "
          f"| {g(d,'uniform'):+.3f} | {g(d,'oracle'):+.3f} | {g(d,'estimated'):+.3f} "
          f"| {d.get('oracle_minus_uniform_r2', float('nan')):+.3f} "
          f"| {d.get('estimated_minus_uniform_r2', float('nan')):+.3f} "
          f"| {d.get('loc_ceiling_auc', float('nan')):.3f} "
          f"| {rec.get('auc_active', float('nan')):.3f} |")
PY
echo "Wrote $BASE/phase10_table.md"
