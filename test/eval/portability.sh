#!/usr/bin/env bash
# Block 1.1 — cross-ancestry PRS portability driver.
#
# Loops over (mix_fraction_B × typed_fraction × seed); each cell runs
#   simulate.py (transfer demography) → sargent → predict.py
# and emits a `<cell>/est.predict.json` with all 4 arms' OOS R² on B-test.
# Resumable: if `est.predict.json` already exists, the cell is skipped.
#
# Aggregates results into `<BASE>/portability_table.md` at the end.
#
# Tunable via environment:
#   BASE             = output root directory (default: ~/test-sargent/block1_1_v1)
#   N_TRAIN, N_TEST  = sample sizes per cell                   (500, 200)
#   LENGTH           = sequence length in bp                   (1000000)
#   NE               = effective population size               (1e4)
#   SPLIT            = transfer split time, generations        (3500 → F_ST ≈ 0.15)
#   H2_SNP, H2_ARG   = target heritabilities                   (0.3, 0.3)
#   CAUSAL_FRAC      = fraction of variants that are causal    (0.01)
#   MIX_FRAC_B       = space-separated train-B fractions       ("0.00 0.05 0.20")
#   TYPED_FRACS      = space-separated typed-SNP fractions     ("0.7 0.1")
#   SEEDS            = space-separated seeds                   ("1 2 3 4 5")
#   SARGENT_EXTRA    = extra args to pass to sargent           ("")
#
# Usage examples:
#   # full v1 (120 fits):
#   BASE=~/test-sargent/block1_1_v1 bash test/eval/portability.sh
#
#   # smoke test (1 cell, fast):
#   BASE=/tmp/b11_test N_TRAIN=80 N_TEST=40 LENGTH=200000 \
#     MIX_FRAC_B=0.05 TYPED_FRACS=0.7 SEEDS=1 \
#     bash test/eval/portability.sh

set -euo pipefail

ROOT=${ROOT:-/Users/acs/git-src/sargent}
SARGENT=${SARGENT:-$ROOT/build/src/progs/sargent}
SIM=${SIM:-$ROOT/scripts/simulate.py}
PRED=${PRED:-$ROOT/test/eval/predict.py}

BASE=${BASE:-$HOME/test-sargent/block1_1_v1}

N_TRAIN=${N_TRAIN:-500}
N_TEST=${N_TEST:-200}
LENGTH=${LENGTH:-1000000}
NE=${NE:-1e4}
SPLIT=${SPLIT:-3500}
H2_SNP=${H2_SNP:-0.3}
H2_ARG=${H2_ARG:-0.3}
CAUSAL_FRAC=${CAUSAL_FRAC:-0.01}
MIX_FRAC_B=${MIX_FRAC_B:-"0.00 0.05 0.20"}
TYPED_FRACS=${TYPED_FRACS:-"0.7 0.1"}
SEEDS=${SEEDS:-"1 2 3 4 5"}
SARGENT_EXTRA=${SARGENT_EXTRA:-""}

# Resolve a `python` that has msprime/tskit. Order of preference:
#   1. $PYBIN if set
#   2. `mamba run -n sargent python`
#   3. plain `python` (assume already in the right env)
if [ -n "${PYBIN:-}" ]; then
  PY="$PYBIN"
elif command -v mamba >/dev/null 2>&1 && mamba env list 2>/dev/null | grep -q '^sargent'; then
  PY="mamba run -n sargent python"
else
  PY="python"
fi

mkdir -p "$BASE"

# Count total cells for progress reporting.
N_MIX=$(echo $MIX_FRAC_B | wc -w | tr -d ' ')
N_TYP=$(echo $TYPED_FRACS | wc -w | tr -d ' ')
N_SD=$(echo $SEEDS | wc -w | tr -d ' ')
N_CELLS=$(( N_MIX * N_TYP * N_SD ))

echo "Block 1.1 v1: portability driver"
echo "  base   : $BASE"
echo "  cells  : ${N_MIX} mix × ${N_TYP} typed × ${N_SD} seeds = ${N_CELLS}"
echo "  scale  : n_train=$N_TRAIN n_test=$N_TEST L=$LENGTH Ne=$NE"
echo "  trait  : h2_snp=$H2_SNP h2_arg=$H2_ARG causal_frac=$CAUSAL_FRAC"
echo "  demog  : transfer split=$SPLIT gen (F_ST ≈ 0.15 at Ne=1e4)"
echo "  python : $PY"
echo "  sargent: $SARGENT $SARGENT_EXTRA"
echo

# Sanity: sargent binary present?
if [ ! -x "$SARGENT" ]; then
  echo "ERROR: sargent binary not found at $SARGENT" >&2
  echo "  Build with: cd $ROOT && cmake -B build && cmake --build build -j" >&2
  exit 1
fi

t_start=$(date +%s)
cell=0
skipped=0

for mix in $MIX_FRAC_B; do
  mtag=$(printf "%.2f" "$mix")
  for typed in $TYPED_FRACS; do
    ttag=$(printf "%.2f" "$typed")
    for seed in $SEEDS; do
      cell=$((cell + 1))
      dir="$BASE/mix${mtag}_typed${ttag}_s${seed}"
      mkdir -p "$dir"

      # Resumability: skip if the final predict.json is already written.
      if [ -s "$dir/est.predict.json" ]; then
        echo "[$cell/$N_CELLS] SKIP $dir (already done)"
        skipped=$((skipped + 1))
        continue
      fi

      echo "[$cell/$N_CELLS] mix=$mtag typed=$ttag seed=$seed -> $dir"
      cstart=$(date +%s)

      # 1. simulate
      $PY "$SIM" \
          --demography transfer --transfer-split-time "$SPLIT" \
          --mix-fraction-B "$mix" \
          --n-indiv "$N_TRAIN" --n-test "$N_TEST" \
          --length "$LENGTH" --Ne "$NE" \
          --h2-snp "$H2_SNP" --h2-arg "$H2_ARG" \
          --causal-fraction "$CAUSAL_FRAC" \
          --bed-coverage "$typed" \
          --seed "$seed" --out "$dir/case" \
          > "$dir/sim.log" 2>&1 || { echo "  sim FAILED, see $dir/sim.log"; tail -10 "$dir/sim.log"; continue; }

      # 2. sargent (full pipeline: BayesR fixed + per-segment-weighted ARG random)
      $SARGENT --silent --nprobe 10 --outroot "$dir/est" $SARGENT_EXTRA \
          "$dir/case.phen" "$dir/case.trees" "$dir/case" \
          > "$dir/sargent.log" 2>&1 || { echo "  sargent FAILED, see $dir/sargent.log"; tail -10 "$dir/sargent.log"; continue; }

      # 3. predict (runs all 4 arms: BayesR / GBLUP / ARG-BLUP-uniform / SARGENT)
      $PY "$PRED" \
          --case-prefix "$dir/case" --sargent-prefix "$dir/est" \
          > "$dir/predict.log" 2>&1 || { echo "  predict FAILED, see $dir/predict.log"; tail -10 "$dir/predict.log"; continue; }

      cend=$(date +%s)
      celapsed=$(( cend - cstart ))
      total_elapsed=$(( cend - t_start ))
      remaining=$(( N_CELLS - cell ))
      avg=$(( cell > skipped ? total_elapsed / (cell - skipped) : 0 ))
      eta=$(( remaining * avg ))
      printf "    done in %ds  (avg %ds/cell, ETA %dm)\n" "$celapsed" "$avg" "$((eta / 60))"
    done
  done
done

echo
echo "All cells processed (${cell} total, ${skipped} skipped)."

# ----- aggregate ---------------------------------------------------------------
echo "Aggregating into $BASE/portability_table.md ..."
$PY - "$BASE" <<'PY' | tee "$BASE/portability_table.md"
import json, os, sys, glob, math, re
import numpy as np

base = sys.argv[1]
rows = []
for jf in sorted(glob.glob(os.path.join(base, "mix*_typed*_s*", "est.predict.json"))):
    d = json.load(open(jf))
    cell = os.path.basename(os.path.dirname(jf))
    # parse mix/typed/seed from cell name
    m = re.match(r"mix([\d.]+)_typed([\d.]+)_s(\d+)", cell)
    if not m:
        continue
    mix, typed, seed = float(m.group(1)), float(m.group(2)), int(m.group(3))
    # try to recover realised h2 / composition from params
    pjf = os.path.join(os.path.dirname(jf), "case.params.json")
    params = json.load(open(pjf)) if os.path.exists(pjf) else {}
    rows.append({
        "cell": cell, "mix": mix, "typed": typed, "seed": seed,
        "n_train_A": params.get("n_train_A"), "n_train_B": params.get("n_train_B"),
        "h2_snp_real": params.get("h2_snp_realised"),
        "h2_arg_real": params.get("h2_arg_realised"),
        "bayesR_r2":   d.get("bayesR", {}).get("r2"),
        "gblup_r2":    d.get("gblup", {}).get("r2"),
        "argunif_r2":  d.get("arg_blup_uniform", {}).get("r2"),
        "sargent_r2":  d.get("sargent", {}).get("r2"),
        "sargent_u_r": d.get("sargent_u_r"),
    })

print("# Block 1.1 — cross-ancestry PRS portability\n")
print(f"Aggregated {len(rows)} cells from `{base}`.\n")

# Per-cell table
print("## Per-seed table\n")
print("| mix_B | typed | seed | BayesR R² | GBLUP R² | ARG-BLUP-uni R² | SARGENT R² | SARGENT u-r |")
print("|---|---|---|---|---|---|---|---|")
def fnum(x, fmt="{:+.3f}"):
    return fmt.format(x) if isinstance(x, (int, float)) and not (isinstance(x, float) and math.isnan(x)) else "—"
for r in rows:
    print(f"| {r['mix']:.2f} | {r['typed']:.2f} | {r['seed']} "
          f"| {fnum(r['bayesR_r2'])} | {fnum(r['gblup_r2'])} "
          f"| {fnum(r['argunif_r2'])} | {fnum(r['sargent_r2'])} "
          f"| {fnum(r['sargent_u_r'])} |")

# Summary by (mix, typed) — mean ± sd over seeds, with paired SARGENT - ARG-BLUP-uni gap
print("\n## Summary (mean ± sd over seeds)\n")
print("| mix_B | typed | n_seeds | BayesR | GBLUP | ARG-BLUP-uni | SARGENT | SARGENT − ARG-BLUP-uni |")
print("|---|---|---|---|---|---|---|---|")
keys = sorted({(r["mix"], r["typed"]) for r in rows})
for (mix, typed) in keys:
    sub = [r for r in rows if r["mix"] == mix and r["typed"] == typed]
    def stat(field):
        vals = [r[field] for r in sub if isinstance(r[field], (int, float))]
        if not vals:
            return "—"
        return f"{np.mean(vals):+.3f} ± {np.std(vals, ddof=1):.3f}" if len(vals) > 1 else f"{vals[0]:+.3f}"
    # paired SARGENT − ARG-BLUP-uni gap
    pairs = [(r["sargent_r2"], r["argunif_r2"]) for r in sub
             if isinstance(r["sargent_r2"], (int, float)) and isinstance(r["argunif_r2"], (int, float))]
    if pairs:
        diffs = np.array([s - u for (s, u) in pairs])
        if len(diffs) > 1:
            gap_str = f"{diffs.mean():+.3f} ± {diffs.std(ddof=1):.3f}"
        else:
            gap_str = f"{diffs[0]:+.3f}"
    else:
        gap_str = "—"
    print(f"| {mix:.2f} | {typed:.2f} | {len(sub)} "
          f"| {stat('bayesR_r2')} | {stat('gblup_r2')} "
          f"| {stat('argunif_r2')} | {stat('sargent_r2')} | {gap_str} |")

# Headline portability curve: how does SARGENT−uniform gap scale with mix_B?
print("\n## Portability story\n")
print("SARGENT − ARG-BLUP-uni paired gap by mix fraction (averaged over typed × seeds):")
for mix in sorted({r["mix"] for r in rows}):
    pairs = [(r["sargent_r2"], r["argunif_r2"]) for r in rows
             if r["mix"] == mix and isinstance(r["sargent_r2"], (int, float))
             and isinstance(r["argunif_r2"], (int, float))]
    if not pairs:
        continue
    diffs = np.array([s - u for (s, u) in pairs])
    print(f"  mix_B = {mix:.2f}:  {diffs.mean():+.4f}  (n={len(diffs)}, sd={diffs.std(ddof=1) if len(diffs)>1 else 0:.4f})")
PY

t_end=$(date +%s)
echo
echo "Total wall time: $(( (t_end - t_start) / 60 )) min"
echo "Wrote $BASE/portability_table.md"
