#!/usr/bin/env bash
# Block 1.2 — localization robustness sweep driver.
#
# For each (demography × block-contrast × seed):
#   1. simulate.py --ws2-architecture sparse with the cell's BLOCKCONTRAST
#   2. sargent --arg-only (estimate ŵ²)
#   3. oracle_ws2.py (compute estimator AUC, ceiling AUC, ws² recovery)
# Then aggregate.
#
# Resumable: skips any cell whose est.oracle.json is non-empty.
#
# Tunable via environment:
#   BASE              = output root (default: ~/test-sargent/block1_2_v1)
#   N_TRAIN, N_TEST   = sample sizes per cell        (300, 100)
#   LENGTH            = sequence length              (300000)
#   NE                = effective population size    (1e4)
#   H2_ARG            = target h²_arg                (0.5)
#   BLOCKFRAC         = active block fraction        (0.05; pinned in v1)
#   NBLOCKS           = number of active blocks      (3)
#   CONTRASTS         = block:floor ws² ratios       ("1000 100 10")
#   DEMOGS            = demographies                 ("admixed panmictic")
#   SPLIT, ADMIXT, ALPHA = admixed parameters        (4000, 20, 0.5)
#   SEEDS             = seed list                    ("1 2 3 4 5")
#
# Usage:
#   # Full v1 (30 cells, ~2h):
#   BASE=~/test-sargent/block1_2_v1 bash test/eval/localization_sweep.sh
#
#   # Smoke test (1 cell, fast):
#   BASE=/tmp/b12_test N_TRAIN=80 LENGTH=80000 \
#     DEMOGS=admixed CONTRASTS=1000 SEEDS=1 \
#     bash test/eval/localization_sweep.sh

set -euo pipefail

ROOT=${ROOT:-/Users/acs/git-src/sargent}
SARGENT=${SARGENT:-$ROOT/build/src/progs/sargent}
SIM=${SIM:-$ROOT/scripts/simulate.py}
ORACLE=${ORACLE:-$ROOT/test/eval/oracle_ws2.py}

BASE=${BASE:-$HOME/test-sargent/block1_2_v1}

N_TRAIN=${N_TRAIN:-300}
N_TEST=${N_TEST:-100}
LENGTH=${LENGTH:-300000}
NE=${NE:-1e4}
H2_ARG=${H2_ARG:-0.5}
BLOCKFRAC=${BLOCKFRAC:-0.05}
NBLOCKS=${NBLOCKS:-3}
CONTRASTS=${CONTRASTS:-"1000 100 10"}
DEMOGS=${DEMOGS:-"admixed panmictic"}
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

N_DEM=$(echo $DEMOGS | wc -w | tr -d ' ')
N_CON=$(echo $CONTRASTS | wc -w | tr -d ' ')
N_SD=$(echo $SEEDS | wc -w | tr -d ' ')
N_CELLS=$(( N_DEM * N_CON * N_SD ))

echo "Block 1.2 v1: localization robustness sweep"
echo "  base   : $BASE"
echo "  cells  : ${N_DEM} demog × ${N_CON} contrast × ${N_SD} seeds = ${N_CELLS}"
echo "  scale  : n_train=$N_TRAIN n_test=$N_TEST L=$LENGTH Ne=$NE h2_arg=$H2_ARG"
echo "  arch   : block_frac=$BLOCKFRAC n_blocks=$NBLOCKS (architecture pinned)"
echo "  python : $PY"
echo "  sargent: $SARGENT"
echo

if [ ! -x "$SARGENT" ]; then
  echo "ERROR: sargent binary not found at $SARGENT" >&2
  echo "  Build with: cd $ROOT && cmake -B build && cmake --build build -j" >&2
  exit 1
fi

t_start=$(date +%s)
cell=0
skipped=0

for demog in $DEMOGS; do
  for contrast in $CONTRASTS; do
    ctag=$(printf "%g" "$contrast")
    for seed in $SEEDS; do
      cell=$((cell + 1))
      dir="$BASE/c${ctag}_d${demog}_s${seed}"
      mkdir -p "$dir"

      if [ -s "$dir/est.oracle.json" ]; then
        echo "[$cell/$N_CELLS] SKIP $dir (already done)"
        skipped=$((skipped + 1))
        continue
      fi

      echo "[$cell/$N_CELLS] contrast=$ctag demog=$demog seed=$seed -> $dir"
      cstart=$(date +%s)

      # admixed needs split/admix args; panmictic ignores them
      DEMOG_ARGS="--demography $demog"
      if [ "$demog" = "admixed" ]; then
        DEMOG_ARGS="$DEMOG_ARGS --admix-split-time $SPLIT --admix-time $ADMIXT --admix-alpha $ALPHA"
      fi

      # 1. simulate
      $PY "$SIM" $DEMOG_ARGS \
          --n-indiv "$N_TRAIN" --n-test "$N_TEST" \
          --length "$LENGTH" --Ne "$NE" \
          --h2-snp 0 --h2-arg "$H2_ARG" \
          --ws2-architecture sparse \
          --ws2-n-blocks "$NBLOCKS" --ws2-block-frac "$BLOCKFRAC" \
          --ws2-block-contrast "$contrast" \
          --seed "$seed" --out "$dir/case" \
          > "$dir/sim.log" 2>&1 || { echo "  sim FAILED, see $dir/sim.log"; tail -10 "$dir/sim.log"; continue; }

      # 2. sargent --arg-only (Stage 2 only; pure random)
      $SARGENT --silent --nprobe 10 --arg-only --outroot "$dir/est" \
          "$dir/case.phen" "$dir/case.trees" "$dir/case" \
          > "$dir/sargent.log" 2>&1 || { echo "  sargent FAILED, see $dir/sargent.log"; tail -10 "$dir/sargent.log"; continue; }

      # 3. oracle_ws2.py (estimator AUC, ceiling, recovery)
      $PY "$ORACLE" --case-prefix "$dir/case" --sargent-prefix "$dir/est" \
          > "$dir/oracle.log" 2>&1 || { echo "  oracle FAILED, see $dir/oracle.log"; tail -10 "$dir/oracle.log"; continue; }

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

echo
echo "All cells processed ($cell total, $skipped skipped)."

# ----- aggregate ---------------------------------------------------------------
echo "Aggregating into $BASE/localization_table.md ..."
$PY - "$BASE" <<'PY' | tee "$BASE/localization_table.md"
import json, os, sys, glob, math, re
import numpy as np

base = sys.argv[1]
rows = []
for jf in sorted(glob.glob(os.path.join(base, "c*_d*_s*", "est.oracle.json"))):
    d = json.load(open(jf))
    cell = os.path.basename(os.path.dirname(jf))
    m = re.match(r"c([\d.]+)_d(\w+)_s(\d+)", cell)
    if not m:
        continue
    contrast, demog, seed = float(m.group(1)), m.group(2), int(m.group(3))
    pjf = os.path.join(os.path.dirname(jf), "case.params.json")
    params = json.load(open(pjf)) if os.path.exists(pjf) else {}
    rec = d.get("ws2_recovery", {}) or {}
    rows.append({
        "cell": cell, "contrast": contrast, "demog": demog, "seed": seed,
        "h2_arg_real": params.get("h2_arg_realised"),
        "est_auc":  rec.get("auc_active"),
        "pearson":  rec.get("pearson"),
        "spearman": rec.get("spearman"),
        "ceiling_auc": d.get("loc_ceiling_auc"),
        "r2_uni": (d.get("uniform") or {}).get("r2"),
        "r2_est": (d.get("estimated") or {}).get("r2"),
        "n_active": rec.get("n_active"),
    })

print("# Block 1.2 — localization robustness sweep\n")
print(f"Aggregated {len(rows)} cells from `{base}`.\n")

def fnum(x, fmt="{:+.3f}"):
    return fmt.format(x) if isinstance(x, (int, float)) and not (isinstance(x, float) and math.isnan(x)) else "—"

# Per-seed table
print("## Per-seed table\n")
print("| demog | contrast | seed | est AUC | ceiling AUC | gap | spearman | R²_uni | R²_est |")
print("|---|---|---|---|---|---|---|---|---|")
for r in rows:
    gap = r["ceiling_auc"] - r["est_auc"] if (isinstance(r["est_auc"], (int, float))
                                              and isinstance(r["ceiling_auc"], (int, float))) else None
    print(f"| {r['demog']} | {r['contrast']:.0f} | {r['seed']} "
          f"| {fnum(r['est_auc'], '{:.3f}')} | {fnum(r['ceiling_auc'], '{:.3f}')} "
          f"| {fnum(gap, '{:.3f}')} | {fnum(r['spearman'])} "
          f"| {fnum(r['r2_uni'])} | {fnum(r['r2_est'])} |")

# Summary by (demog × contrast)
print("\n## Summary (mean ± sd over seeds)\n")
print("| demog | contrast | n_seeds | est AUC | ceiling AUC | gap (ceiling−est) | spearman |")
print("|---|---|---|---|---|---|---|")
keys = sorted({(r["demog"], r["contrast"]) for r in rows}, key=lambda x: (x[0], -x[1]))
def stat(sub, field):
    vals = [r[field] for r in sub if isinstance(r[field], (int, float))]
    if not vals:
        return "—"
    return f"{np.mean(vals):.3f} ± {np.std(vals, ddof=1):.3f}" if len(vals) > 1 else f"{vals[0]:.3f}"
for demog, contrast in keys:
    sub = [r for r in rows if r["demog"] == demog and r["contrast"] == contrast]
    gaps = [(r["ceiling_auc"] - r["est_auc"]) for r in sub
            if isinstance(r["est_auc"], (int, float)) and isinstance(r["ceiling_auc"], (int, float))]
    if gaps and len(gaps) > 1:
        gap_str = f"{np.mean(gaps):.3f} ± {np.std(gaps, ddof=1):.3f}"
    elif gaps:
        gap_str = f"{gaps[0]:.3f}"
    else:
        gap_str = "—"
    print(f"| {demog} | {contrast:.0f} | {len(sub)} | {stat(sub, 'est_auc')} | "
          f"{stat(sub, 'ceiling_auc')} | {gap_str} | {stat(sub, 'spearman')} |")

# Trend per demog
print("\n## Robustness trend\n")
for demog in sorted({r["demog"] for r in rows}):
    print(f"\n**{demog}** — estimator AUC vs contrast (mean over seeds):")
    for contrast in sorted({r["contrast"] for r in rows}, reverse=True):
        vals = [r["est_auc"] for r in rows
                if r["demog"] == demog and r["contrast"] == contrast
                and isinstance(r["est_auc"], (int, float))]
        if vals:
            label = "above-chance" if np.mean(vals) > 0.55 else "near-chance"
            print(f"  contrast={contrast:.0f}: {np.mean(vals):.3f}  ({label}, n={len(vals)})")
PY

t_end=$(date +%s)
echo
echo "Total wall time: $(( (t_end - t_start) / 60 )) min"
echo "Wrote $BASE/localization_table.md"
