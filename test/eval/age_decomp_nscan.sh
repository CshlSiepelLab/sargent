#!/usr/bin/env bash
# Block 2.1 v2 — n-scaling test for ARG-based age decomposition.
#
# Tests whether SARGENT's per-tree ŵ²_s recovers the natural genealogy's
# per-tree-age-bin variance distribution better with larger n.
#
# Sim: h²_arg = 0.5, h²_snp = 0, default uniform ws² (no artificial age
# stratification). The TRUE per-tree-age decomposition is determined by
# the coalescent (branch lengths × subtree sizes binned by mean branch
# age); SARGENT's per-tree REML EM should approach this truth as n grows.
#
# Grid: n ∈ {300, 600, 1000}, 3 seeds each = 9 cells.
# L = 100 kb (smaller than Block 1.2's 300 kb to keep n=1000 feasible
# overnight). Panmictic (no demography complexity).
#
# Per cell: simulate → sargent --arg-only → arch_metrics post-hoc to get
# per-tree-age-bin h² for BOTH truth and est.

set -euo pipefail

ROOT=${ROOT:-/Users/acs/git-src/sargent}
SARGENT=${SARGENT:-$ROOT/build/src/progs/sargent}
SIM=${SIM:-$ROOT/scripts/simulate.py}

BASE=${BASE:-$HOME/test-sargent/block2_1_v2_nscan}

LENGTH=${LENGTH:-100000}
NE=${NE:-1e4}
H2_ARG=${H2_ARG:-0.5}
NS=${NS:-"300 600 1000"}
SEEDS=${SEEDS:-"1 2 3"}

if [ -n "${PYBIN:-}" ]; then
  PY="$PYBIN"
elif command -v mamba >/dev/null 2>&1 && mamba env list 2>/dev/null | grep -q '^sargent'; then
  PY="mamba run -n sargent python"
else
  PY="python"
fi

mkdir -p "$BASE"

N_N=$(echo $NS | wc -w | tr -d ' ')
N_SD=$(echo $SEEDS | wc -w | tr -d ' ')
N_CELLS=$(( N_N * N_SD ))

echo "Block 2.1 v2: n-scaling test for ARG-based age decomposition"
echo "  base   : $BASE"
echo "  cells  : ${N_N} n × ${N_SD} seeds = ${N_CELLS}"
echo "  scale  : L=$LENGTH Ne=$NE"
echo "  trait  : h2_arg=$H2_ARG  h2_snp=0  uniform ws² (no artificial age stratification)"
echo

if [ ! -x "$SARGENT" ]; then
  echo "ERROR: sargent binary not found at $SARGENT" >&2
  exit 1
fi

t_start=$(date +%s)
cell=0
skipped=0

for n in $NS; do
  for seed in $SEEDS; do
    cell=$((cell + 1))
    n_test=$(( n / 3 ))
    dir="$BASE/n${n}_s${seed}"
    mkdir -p "$dir"

    if [ -s "$dir/est.w2.tsv" ] && [ -s "$dir/case.true_ws2.tsv" ]; then
      echo "[$cell/$N_CELLS] SKIP $dir (already done)"
      skipped=$((skipped + 1))
      continue
    fi

    echo "[$cell/$N_CELLS] n=$n seed=$seed -> $dir"
    cstart=$(date +%s)

    $PY "$SIM" --demography panmictic \
        --n-indiv "$n" --n-test "$n_test" \
        --length "$LENGTH" --Ne "$NE" \
        --h2-snp 0 --h2-arg "$H2_ARG" \
        --seed "$seed" --out "$dir/case" \
        > "$dir/sim.log" 2>&1 || { echo "  sim FAILED"; tail -5 "$dir/sim.log"; continue; }

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

echo
echo "All cells processed ($cell total, $skipped skipped)."
echo "Aggregating per-tree-age-binned h² (truth vs est) ..."

$PY - "$BASE" "$ROOT" <<'PY' | tee "$BASE/age_decomp_nscan_table.md"
import json, os, sys, glob, re
import numpy as np
sys.path.insert(0, os.path.join(sys.argv[2], "test", "eval"))
import tskit
from arch_metrics import per_tree_age_binned_h2

base = sys.argv[1]
# Define common age edges across all cells so the bins are comparable.
# Use coalescent-typical generations at Ne=1e4: split [0, 5000, 20000, ∞]
COMMON_EDGES = [-1.0, 5000.0, 20000.0, 1e9]

rows = []
for d in sorted(glob.glob(os.path.join(base, "n*_s*"))):
    cell = os.path.basename(d)
    m = re.match(r"n(\d+)_s(\d+)", cell)
    if not m: continue
    n, seed = int(m.group(1)), int(m.group(2))
    case_truth = os.path.join(d, "case.full.true_ws2.tsv")
    if not os.path.exists(case_truth): continue
    w_truth = np.loadtxt(case_truth, skiprows=1)[:,-1]
    ts_full = tskit.load(os.path.join(d, "case.full.trees"))
    truth_r = per_tree_age_binned_h2(w_truth, ts_full, age_edges=COMMON_EDGES)

    est_path = os.path.join(d, "est.w2.tsv")
    if not os.path.exists(est_path): continue
    w_est = np.loadtxt(est_path, skiprows=1)[:,-1]
    ts_train = tskit.load(os.path.join(d, "case.trees"))
    est_r = per_tree_age_binned_h2(w_est, ts_train, age_edges=COMMON_EDGES)

    rows.append({"n": n, "seed": seed,
                 "truth": truth_r["fraction"], "est": est_r["fraction"],
                 "truth_ntrees": truth_r["n_trees_per_bin"],
                 "est_ntrees":   est_r["n_trees_per_bin"]})

print("# Block 2.1 v2 — n-scaling test for ARG-based age decomposition\n")
print(f"Aggregated {len(rows)} cells.\n")
print(f"Common age bins (gen): {COMMON_EDGES}\n")
print(f"Trait: h²_arg = 0.5, h²_snp = 0, default uniform ws².\n")
print("## Per-n truth vs estimated h² fraction per age bin (mean ± sd over seeds)\n")
print("| n | bin | truth frac | est frac | bias (est − truth) | mean trees |")
print("|---|---|---|---|---|---|")
from collections import defaultdict
by_n = defaultdict(lambda: {"truth": [[],[],[]], "est": [[],[],[]], "nt": [[],[],[]]})
for r in rows:
    for c in range(3):
        by_n[r["n"]]["truth"][c].append(r["truth"][c])
        by_n[r["n"]]["est"][c].append(r["est"][c])
        by_n[r["n"]]["nt"][c].append(r["truth_ntrees"][c])
labels = ["recent (≤5k gen)", "middle (5k–20k)", "ancient (>20k)"]
for n in sorted(by_n):
    for c, lbl in enumerate(labels):
        t = np.array(by_n[n]["truth"][c]); e = np.array(by_n[n]["est"][c]); nt = np.array(by_n[n]["nt"][c])
        tm, ts_ = (t.mean(), t.std(ddof=1) if len(t)>1 else 0)
        em, es_ = (e.mean(), e.std(ddof=1) if len(e)>1 else 0)
        print(f"| {n} | {lbl} | {tm:.3f} ± {ts_:.3f} | {em:.3f} ± {es_:.3f} | {em-tm:+.3f} | {nt.mean():.1f} |")

print("\n## Bias trend (est − truth, mean over seeds) by n\n")
print("| n | recent | middle | ancient | total |L1 bias| |")
print("|---|---|---|---|---|")
for n in sorted(by_n):
    biases = []
    for c in range(3):
        t = np.mean(by_n[n]["truth"][c]); e = np.mean(by_n[n]["est"][c])
        biases.append(e - t)
    l1 = sum(abs(b) for b in biases)
    print(f"| {n} | {biases[0]:+.3f} | {biases[1]:+.3f} | {biases[2]:+.3f} | {l1:.3f} |")
PY

t_end=$(date +%s)
echo
echo "Total wall: $(( (t_end - t_start) / 60 )) min"
