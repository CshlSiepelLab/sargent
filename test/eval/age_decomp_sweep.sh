#!/usr/bin/env bash
# Block 2.1 v1 — controlled allele-age decomposition sweep.
#
# Tests: when traits are simulated with age-stratified causal effect sizes
# (via --age-tau-coef), does SARGENT's BayesR fit recover per-age-bin h²?
#
# v1 scope: FIXED-effect side only (h²_snp > 0, h²_arg = 0).  Age
# stratification on the random-effect side requires extending simulate.py
# (a v2 extension; the current age-tau only modifies fixed-effect β²).
#
# Sweep:
#   3 tau values × 5 seeds = 15 cells
#   n_train = 300, L = 300 kb, h²_snp = 0.5, h²_arg = 0
#   --bed-coverage = 1.0 (all causal SNPs typed; clean test of age-bin recovery
#                         without coverage confounding)
#
# Per cell: simulate → sargent --snp-only → age_partition.py.
# Aggregator: per-tau truth-vs-est, 3 equal-quantile bins.
#
# Resumable: skips cells with the partition JSON present.

set -euo pipefail

ROOT=${ROOT:-/Users/acs/git-src/sargent}
SARGENT=${SARGENT:-$ROOT/build/src/progs/sargent}
SIM=${SIM:-$ROOT/scripts/simulate.py}
PARTITION=${PARTITION:-$ROOT/test/eval/age_partition.py}

BASE=${BASE:-$HOME/test-sargent/block2_1_v1}

N_TRAIN=${N_TRAIN:-300}
N_TEST=${N_TEST:-100}
LENGTH=${LENGTH:-300000}
NE=${NE:-1e4}
H2_SNP=${H2_SNP:-0.5}
H2_ARG=${H2_ARG:-0.0}
CAUSAL_FRAC=${CAUSAL_FRAC:-0.02}
BED_COV=${BED_COV:-1.0}
AGE_REF=${AGE_REF:-10000}
TAUS=${TAUS:-"-1.0 0.0 1.0"}
SEEDS=${SEEDS:-"1 2 3 4 5"}

if [ -n "${PYBIN:-}" ]; then
  PY="$PYBIN"
elif command -v mamba >/dev/null 2>&1 && mamba env list 2>/dev/null | grep -q '^sargent'; then
  PY="mamba run -n sargent python"
else
  PY="python"
fi

mkdir -p "$BASE"

N_T=$(echo $TAUS | wc -w | tr -d ' ')
N_SD=$(echo $SEEDS | wc -w | tr -d ' ')
N_CELLS=$(( N_T * N_SD ))

echo "Block 2.1 v1: allele-age decomposition sweep"
echo "  base   : $BASE"
echo "  cells  : ${N_T} tau × ${N_SD} seeds = ${N_CELLS}"
echo "  scale  : n_train=$N_TRAIN n_test=$N_TEST L=$LENGTH Ne=$NE"
echo "  trait  : h2_snp=$H2_SNP h2_arg=$H2_ARG causal_frac=$CAUSAL_FRAC bed_cov=$BED_COV"
echo "  tau    : TAUS={$TAUS} age_ref=$AGE_REF gen"
echo

if [ ! -x "$SARGENT" ]; then
  echo "ERROR: sargent binary not found at $SARGENT" >&2
  exit 1
fi

t_start=$(date +%s)
cell=0
skipped=0

for tau in $TAUS; do
  ttag=$(printf "%g" "$tau")
  for seed in $SEEDS; do
    cell=$((cell + 1))
    dir="$BASE/tau${ttag}_s${seed}"
    mkdir -p "$dir"

    if [ -s "$dir/case.age_partition.json" ]; then
      echo "[$cell/$N_CELLS] SKIP $dir (already done)"
      skipped=$((skipped + 1))
      continue
    fi

    echo "[$cell/$N_CELLS] tau=$ttag seed=$seed -> $dir"
    cstart=$(date +%s)

    $PY "$SIM" \
        --n-indiv "$N_TRAIN" --n-test "$N_TEST" \
        --length "$LENGTH" --Ne "$NE" \
        --h2-snp "$H2_SNP" --h2-arg "$H2_ARG" \
        --causal-fraction "$CAUSAL_FRAC" \
        --bed-coverage "$BED_COV" \
        --age-tau-coef "$tau" --age-tau-ref "$AGE_REF" \
        --seed "$seed" --out "$dir/case" \
        > "$dir/sim.log" 2>&1 || { echo "  sim FAILED"; tail -5 "$dir/sim.log"; continue; }

    $SARGENT --silent --nprobe 10 --snp-only --outroot "$dir/est" \
        "$dir/case.phen" "$dir/case.trees" "$dir/case" \
        > "$dir/sargent.log" 2>&1 || { echo "  sargent FAILED"; tail -5 "$dir/sargent.log"; continue; }

    $PY "$PARTITION" --case-prefix "$dir/case" --sargent-prefix "$dir/est" \
        --bins 3 --out "$dir/case.age_partition.json" \
        > "$dir/partition.log" 2>&1 || { echo "  partition FAILED"; tail -5 "$dir/partition.log"; continue; }

    cend=$(date +%s)
    celapsed=$(( cend - cstart ))
    total_elapsed=$(( cend - t_start ))
    done_count=$(( cell - skipped ))
    avg=$(( done_count > 0 ? total_elapsed / done_count : 0 ))
    eta=$(( (N_CELLS - cell) * avg ))
    printf "    done in %ds  (avg %ds/cell, ETA %dm)\n" "$celapsed" "$avg" "$((eta / 60))"
  done
done

echo
echo "All cells processed ($cell total, $skipped skipped)."
echo "Aggregating ..."

$PY - "$BASE" <<'PY' | tee "$BASE/age_decomp_table.md"
import json, os, sys, glob, re
import numpy as np
from collections import defaultdict

base = sys.argv[1]
rows = []
for jf in sorted(glob.glob(os.path.join(base, "tau*_s*", "case.age_partition.json"))):
    d = json.load(open(jf))
    cell = os.path.basename(os.path.dirname(jf))
    m = re.match(r"tau([-\d.]+)_s(\d+)", cell)
    if not m: continue
    tau, seed = float(m.group(1)), int(m.group(2))
    rows.append({"tau": tau, "seed": seed, "data": d, "cell": cell})

# Stratify by tau, average per-bin h²-fraction across seeds.
print("# Block 2.1 v1 — allele-age decomposition\n")
print(f"Aggregated {len(rows)} cells from `{base}`.\n")
print("Trait: h²_snp = 0.5, h²_arg = 0, causal_fraction = 0.02, n_train = 300, L = 300 kb, bed_coverage = 1.0.\n")
print("Three equal-quantile age bins (per cell; edges differ across cells but the relative shift in variance is comparable).\n")
print("## Per-tau truth vs estimated h² fraction per bin (mean ± sd over seeds)\n")
print("| tau | bin | truth frac | est frac | bias (est − truth) | n_in_bin (mean) |")
print("|---|---|---|---|---|---|")
by_tau = defaultdict(lambda: [[], [], []])
by_tau_est = defaultdict(lambda: [[], [], []])
by_tau_nin = defaultdict(lambda: [[], [], []])
for r in rows:
    d = r["data"]
    if "truth" not in d: continue
    for c in range(3):
        by_tau[r["tau"]][c].append(d["truth"]["h2_frac"][c])
        by_tau_nin[r["tau"]][c].append(d["truth"]["n_in_bin"][c])
        if "est" in d:
            by_tau_est[r["tau"]][c].append(d["est"]["h2_frac"][c])

for tau in sorted(by_tau):
    for c, label in enumerate(["recent", "middle", "ancient"]):
        t = np.array(by_tau[tau][c])
        e = np.array(by_tau_est[tau][c]) if by_tau_est[tau][c] else np.array([])
        n = np.array(by_tau_nin[tau][c])
        tm, ts = (t.mean(), t.std(ddof=1) if len(t)>1 else 0)
        if len(e):
            em, es = (e.mean(), e.std(ddof=1) if len(e)>1 else 0)
            bias = em - tm
        else:
            em, es, bias = (float('nan'),)*3
        print(f"| {tau:+.1f} | {label} | {tm:.3f} ± {ts:.3f} | {em:.3f} ± {es:.3f} | {bias:+.3f} | {n.mean():.1f} |")

# Trend by tau across the recent bin (cleanest single-axis summary)
print("\n## Trend: recent-bin h² fraction as a function of tau\n")
print("(Lower tau → more variance on recent variants; should rise with negative tau.)\n")
print("| tau | truth recent | est recent | difference |")
print("|---|---|---|---|")
for tau in sorted(by_tau):
    tm = np.mean(by_tau[tau][0])
    em = np.mean(by_tau_est[tau][0]) if by_tau_est[tau][0] else float('nan')
    print(f"| {tau:+.1f} | {tm:.3f} | {em:.3f} | {em-tm:+.3f} |")
PY

t_end=$(date +%s)
echo
echo "Total wall: $(( (t_end - t_start) / 60 )) min"
