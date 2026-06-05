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
#   SPLIT            = transfer split time, generations (legacy single value;
#                       overridden by SPLITS if SPLITS is set)
#   SPLITS           = space-separated split times for F_ST ladder
#                       (Block 1.1 v2; default "3500", i.e. legacy single-F_ST)
#                       e.g. "1500 3500 7000" → F_ST ≈ {0.05, 0.15, 0.25}
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

# Auto-detect the repo root from this script's own location (portable across
# the laptop and the server); override with ROOT=... if needed.
ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}
SARGENT=${SARGENT:-$ROOT/build/src/progs/sargent}
SIM=${SIM:-$ROOT/scripts/simulate.py}
PRED=${PRED:-$ROOT/test/eval/predict.py}

BASE=${BASE:-$HOME/test-sargent/block1_1_v1}

N_TRAIN=${N_TRAIN:-500}
N_TEST=${N_TEST:-200}
LENGTH=${LENGTH:-1000000}
NE=${NE:-1e4}
SPLIT=${SPLIT:-3500}
# Block 1.1 v2 — F_ST ladder. If SPLITS is unset, use the legacy single SPLIT.
SPLITS=${SPLITS:-"$SPLIT"}
H2_SNP=${H2_SNP:-0.3}
H2_ARG=${H2_ARG:-0.3}
CAUSAL_FRAC=${CAUSAL_FRAC:-0.01}
MIX_FRAC_B=${MIX_FRAC_B:-"0.00 0.05 0.20"}
TYPED_FRACS=${TYPED_FRACS:-"0.7 0.1"}
SEEDS=${SEEDS:-"1 2 3 4 5"}
SARGENT_EXTRA=${SARGENT_EXTRA:-""}
# Extra flags passed to simulate.py (NOT the sargent binary).  Used e.g. for
# Block 1.1 v3a to set the per-tree ws² truth architecture
# (--ws2-architecture sparse --ws2-block-frac ... --ws2-block-contrast ...),
# which is a SIMULATION property, not an inference flag.
SIM_EXTRA=${SIM_EXTRA:-""}

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
N_SP=$(echo $SPLITS | wc -w | tr -d ' ')
N_MIX=$(echo $MIX_FRAC_B | wc -w | tr -d ' ')
N_TYP=$(echo $TYPED_FRACS | wc -w | tr -d ' ')
N_SD=$(echo $SEEDS | wc -w | tr -d ' ')
N_CELLS=$(( N_SP * N_MIX * N_TYP * N_SD ))

echo "Block 1.1 portability driver"
echo "  base   : $BASE"
echo "  cells  : ${N_SP} split × ${N_MIX} mix × ${N_TYP} typed × ${N_SD} seeds = ${N_CELLS}"
echo "  scale  : n_train=$N_TRAIN n_test=$N_TEST L=$LENGTH Ne=$NE"
echo "  trait  : h2_snp=$H2_SNP h2_arg=$H2_ARG causal_frac=$CAUSAL_FRAC"
echo "  demog  : transfer SPLITS={$SPLITS} gen"
echo "  python : $PY"
echo "  sargent: $SARGENT $SARGENT_EXTRA"
echo "  sim-extra: ${SIM_EXTRA:-(none)}"
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

for split in $SPLITS; do
 # Cell-name encoding: when only one SPLIT is in play, keep the legacy
 # name (no sp${split}_ prefix) so v1 outputs remain resumable. When
 # multiple SPLITS are swept, prefix with sp${split}_ to keep cells
 # separate per F_ST.
 sptag_prefix=""
 if [ "$N_SP" -gt 1 ]; then sptag_prefix="sp${split}_"; fi
 for mix in $MIX_FRAC_B; do
  mtag=$(printf "%.2f" "$mix")
  for typed in $TYPED_FRACS; do
    ttag=$(printf "%.2f" "$typed")
    for seed in $SEEDS; do
      cell=$((cell + 1))
      dir="$BASE/${sptag_prefix}mix${mtag}_typed${ttag}_s${seed}"
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
          --demography transfer --transfer-split-time "$split" \
          --mix-fraction-B "$mix" \
          --n-indiv "$N_TRAIN" --n-test "$N_TEST" \
          --length "$LENGTH" --Ne "$NE" \
          --h2-snp "$H2_SNP" --h2-arg "$H2_ARG" \
          --causal-fraction "$CAUSAL_FRAC" \
          --bed-coverage "$typed" $SIM_EXTRA \
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
# Accept both legacy (mix...) and SPLITS-prefixed (sp<split>_mix...) cell names.
for jf in sorted(glob.glob(os.path.join(base, "*mix*_typed*_s*", "est.predict.json"))):
    d = json.load(open(jf))
    cell = os.path.basename(os.path.dirname(jf))
    # Optional 'sp<split>_' prefix for SPLITS-sweep mode (Block 1.1 v2).
    m = re.match(r"(?:sp([\d.]+)_)?mix([\d.]+)_typed([\d.]+)_s(\d+)", cell)
    if not m:
        continue
    split = float(m.group(1)) if m.group(1) else None
    mix, typed, seed = float(m.group(2)), float(m.group(3)), int(m.group(4))
    # try to recover realised h2 / composition from params
    pjf = os.path.join(os.path.dirname(jf), "case.params.json")
    params = json.load(open(pjf)) if os.path.exists(pjf) else {}
    def r2c(arm):
        """Calibration-invariant R² = squared Pearson r. Robust to the
        cross-ancestry scale blowup that makes raw r2 misleading for the
        fixed-effect arm. Uses stored r2_calib (new predict.py) when present,
        else derives r² from the stored r (older JSONs)."""
        a = d.get(arm, {})
        if isinstance(a.get("r2_calib"), (int, float)):
            return a["r2_calib"]
        rr = a.get("r")
        return rr * rr if isinstance(rr, (int, float)) else None
    rows.append({
        "cell": cell, "split": split, "mix": mix, "typed": typed, "seed": seed,
        "n_train_A": params.get("n_train_A"), "n_train_B": params.get("n_train_B"),
        "fst":         params.get("fst_train_test"),   # present for stdpopsim cells
        "demography":  params.get("demography"),
        "train_pop":   params.get("train_pop"),
        "test_pop":    params.get("test_pop"),
        "h2_snp_real": params.get("h2_snp_realised"),
        "h2_arg_real": params.get("h2_arg_realised"),
        "bayesR_r2":   d.get("bayesR", {}).get("r2"),
        "gblup_r2":    d.get("gblup", {}).get("r2"),
        "argunif_r2":  d.get("arg_blup_uniform", {}).get("r2"),
        "sargent_r2":  d.get("sargent", {}).get("r2"),
        "sargent_u_r": d.get("sargent_u_r"),
        "bayesR_r2c":  r2c("bayesR"),
        "gblup_r2c":   r2c("gblup"),
        "argunif_r2c": r2c("arg_blup_uniform"),
        "sargent_r2c": r2c("sargent"),
    })

print("# Block 1.1 — cross-ancestry PRS portability\n")
print(f"Aggregated {len(rows)} cells from `{base}`.\n")

# Per-cell table
print("## Per-seed table\n")
has_fst = any(isinstance(r["fst"], (int, float)) for r in rows)
def fnum(x, fmt="{:+.3f}"):
    return fmt.format(x) if isinstance(x, (int, float)) and not (isinstance(x, float) and math.isnan(x)) else "—"
if has_fst:
    print("| mix_B | typed | seed | F_ST | BayesR R² | GBLUP R² | ARG-BLUP-uni R² | SARGENT R² | SARGENT u-r |")
    print("|---|---|---|---|---|---|---|---|---|")
else:
    print("| mix_B | typed | seed | BayesR R² | GBLUP R² | ARG-BLUP-uni R² | SARGENT R² | SARGENT u-r |")
    print("|---|---|---|---|---|---|---|---|")
for r in rows:
    fst_col = f" {fnum(r['fst'], '{:.3f}')} |" if has_fst else ""
    print(f"| {r['mix']:.2f} | {r['typed']:.2f} | {r['seed']} |{fst_col}"
          f" {fnum(r['bayesR_r2'])} | {fnum(r['gblup_r2'])} "
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
print("\n## Portability story (raw R²)\n")
print("SARGENT − ARG-BLUP-uni paired gap by mix fraction (averaged over typed × seeds):")
for mix in sorted({r["mix"] for r in rows}):
    pairs = [(r["sargent_r2"], r["argunif_r2"]) for r in rows
             if r["mix"] == mix and isinstance(r["sargent_r2"], (int, float))
             and isinstance(r["argunif_r2"], (int, float))]
    if not pairs:
        continue
    diffs = np.array([s - u for (s, u) in pairs])
    print(f"  mix_B = {mix:.2f}:  {diffs.mean():+.4f}  (n={len(diffs)}, sd={diffs.std(ddof=1) if len(diffs)>1 else 0:.4f})")

# Calibration-invariant headline: raw R² is destroyed by the cross-ancestry
# scale blowup of the fixed-effect arm (a single rare-in-train/common-in-test
# SNP can send BayesR/SARGENT R² to −3); the squared-correlation gap is the
# fair discrimination comparison. Prefer this for the verdict.
print("\n## Portability story — CALIBRATION-INVARIANT (squared correlation r²)\n")
print("SARGENT − ARG-BLUP-uni paired gap on r²_calib by mix fraction:")
for mix in sorted({r["mix"] for r in rows}):
    pairs = [(r["sargent_r2c"], r["argunif_r2c"]) for r in rows
             if r["mix"] == mix and isinstance(r["sargent_r2c"], (int, float))
             and isinstance(r["argunif_r2c"], (int, float))]
    if not pairs:
        continue
    diffs = np.array([s - u for (s, u) in pairs])
    print(f"  mix_B = {mix:.2f}:  {diffs.mean():+.4f}  (n={len(diffs)}, sd={diffs.std(ddof=1) if len(diffs)>1 else 0:.4f})")

# ---- F_ST-stratified portability surface (Block 1.1 v2) -----------------------
# Only meaningful when SPLITS were swept (cells carry an sp<split>_ prefix).
splits = sorted({r["split"] for r in rows if r["split"] is not None})
if splits:
    FST = {1500.0: 0.07, 2500.0: 0.12, 3500.0: 0.16, 7000.0: 0.30}
    def amean(rs, k):
        v = [r[k] for r in rs if isinstance(r[k], (int, float))]
        return float(np.mean(v)) if v else float("nan")
    def fst(sp):
        return FST.get(sp, float("nan"))

    print("\n## Portability surface by F_ST — CALIBRATION-INVARIANT r²\n")
    print("Mean r²_calib over (mix × typed × seed) per split:\n")
    print("| split | F_ST | BayesR | GBLUP | ARG-uni | SARGENT | ARG-uni−GBLUP | ARG-uni−BayesR |")
    print("|---|---|---|---|---|---|---|---|")
    for sp in splits:
        rs = [r for r in rows if r["split"] == sp]
        b, g, a, s = (amean(rs, "bayesR_r2c"), amean(rs, "gblup_r2c"),
                      amean(rs, "argunif_r2c"), amean(rs, "sargent_r2c"))
        print(f"| {int(sp)} | {fst(sp):.2f} | {b:+.3f} | {g:+.3f} | {a:+.3f} | {s:+.3f} "
              f"| {a-g:+.3f} | {a-b:+.3f} |")

    print("\nFramework−typed gap by F_ST "
          "(framework=max(ARG-uni,SARGENT), typed=max(BayesR,GBLUP); per-cell paired):")
    for sp in splits:
        gaps = []
        for r in (x for x in rows if x["split"] == sp):
            fw = [v for v in (r["argunif_r2c"], r["sargent_r2c"]) if isinstance(v, (int, float))]
            tp = [v for v in (r["bayesR_r2c"], r["gblup_r2c"]) if isinstance(v, (int, float))]
            if fw and tp:
                gaps.append(max(fw) - max(tp))
        if gaps:
            ga = np.array(gaps)
            sd = ga.std(ddof=1) if len(ga) > 1 else 0.0
            print(f"  F_ST={fst(sp):.2f}: {ga.mean():+.4f} ± {sd:.4f}  (n={len(ga)})")

    print("\nARG-uni vs GBLUP by F_ST × typed coverage "
          "(ARG kinship is coverage-invariant; does low coverage widen the gap?):")
    for sp in splits:
        for tc in sorted({r["typed"] for r in rows if r["split"] == sp}):
            rs = [r for r in rows if r["split"] == sp and r["typed"] == tc]
            gg, aa = amean(rs, "gblup_r2c"), amean(rs, "argunif_r2c")
            print(f"  F_ST={fst(sp):.2f} typed={tc:.2f}: GBLUP={gg:+.3f}  ARG-uni={aa:+.3f}  (gap {aa-gg:+.3f})")
PY

t_end=$(date +%s)
echo
echo "Total wall time: $(( (t_end - t_start) / 60 )) min"
echo "Wrote $BASE/portability_table.md"
