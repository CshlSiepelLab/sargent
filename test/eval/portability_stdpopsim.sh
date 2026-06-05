#!/usr/bin/env bash
# Block 1.1 v3 — cross-ancestry portability driver for stdpopsim demographies.
#
# Mirrors portability.sh but drives the stdpopsim simulation path:
#   simulate.py (--demography stdpopsim) → sargent → predict.py
# Resumable: skips any cell whose est.predict.json already exists.
#
# Supported sweep types (SWEEP env var):
#   OOA     OutOfAfrica_3G09    — 3 pop-pairs × 2 typed × 2 mix_B × 5 seeds (60)
#   ADM-tr  AmericanAdmixture_4B18 — EUR train → ADMIX test,
#                                    2 typed × 2 mix_B × 5 seeds (20)
#   ADM-te  AmericanAdmixture_4B18 — EUR train → ADMIX test mix=0,
#                                    2 typed × 5 seeds (10, admixed-test focus)
#   ARC     OutOfAfricaArchaicAdmixture_5R19 — CEU → Nean@t=1000,
#                                    2 typed × 5 seeds (10; Neanderthal probe)
#   ALL     run all four sweeps (100 cells total)
#
# Output: <BASE>/<model>_<trainpop>x<testpop>_mix<m>_typed<t>_s<seed>/
#
# Tunable via environment:
#   BASE         output root                      (~/test-sargent/block1_1_v3)
#   SWEEP        which sub-sweep(s)               (ALL)
#   N_TRAIN      training sample size             (500)
#   N_TEST       test sample size                 (200)
#   LENGTH       sequence length in bp            (1000000)
#   H2_SNP       target h²_snp                   (0.3)
#   H2_ARG       target h²_arg                   (0.3)
#   CAUSAL_FRAC  fraction causal SNPs             (0.01)
#   MIX_FRAC_B   space-sep train-B fractions      ("0.00 0.05")
#   TYPED_FRACS  space-sep typed-SNP fractions    ("0.7 0.1")
#   SEEDS        space-sep seeds                  ("1 2 3 4 5")
#   NEAN_TIME    Neanderthal sampling time (gen)  (1000)
#   SARGENT_EXTRA  extra args to sargent          ("")
#
# Usage:
#   # Full v3 (all sweeps):
#   BASE=~/test-sargent/block1_1_v3 bash test/eval/portability_stdpopsim.sh
#
#   # OOA only:
#   SWEEP=OOA BASE=~/test-sargent/block1_1_v3_ooa bash test/eval/portability_stdpopsim.sh
#
#   # Pilot (1 cell, fast):
#   SWEEP=OOA N_TRAIN=100 N_TEST=50 LENGTH=300000 MIX_FRAC_B=0.00 \
#     TYPED_FRACS=0.7 SEEDS=1 BASE=/tmp/v3_pilot \
#     bash test/eval/portability_stdpopsim.sh

set -euo pipefail

ROOT=${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}
SARGENT=${SARGENT:-$ROOT/build/src/progs/sargent}
SIM=${SIM:-$ROOT/scripts/simulate.py}
PRED=${PRED:-$ROOT/test/eval/predict.py}

BASE=${BASE:-$HOME/test-sargent/block1_1_v3}
SWEEP=${SWEEP:-ALL}

N_TRAIN=${N_TRAIN:-500}
N_TEST=${N_TEST:-200}
LENGTH=${LENGTH:-1000000}
H2_SNP=${H2_SNP:-0.3}
H2_ARG=${H2_ARG:-0.3}
CAUSAL_FRAC=${CAUSAL_FRAC:-0.01}
MIX_FRAC_B=${MIX_FRAC_B:-"0.00 0.05"}
TYPED_FRACS=${TYPED_FRACS:-"0.7 0.1"}
SEEDS=${SEEDS:-"1 2 3 4 5"}
NEAN_TIME=${NEAN_TIME:-1000}
SARGENT_EXTRA=${SARGENT_EXTRA:-""}

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

if ! $PY -c "import stdpopsim" 2>/dev/null; then
  echo "ERROR: stdpopsim not installed. Run: pip install stdpopsim" >&2
  exit 1
fi

mkdir -p "$BASE"

# ---- build cell list -----------------------------------------------------------
# Each entry: "model|train_pop|test_pop|test_pop_time|mix|typed|seed"
CELLS=()

add_cells() {
  local model="$1" train="$2" test="$3" test_time="${4:-0}"
  local use_mix="${5:-yes}"   # "no" → only mix=0.00
  for typed in $TYPED_FRACS; do
    for seed in $SEEDS; do
      if [ "$use_mix" = "yes" ]; then
        for mix in $MIX_FRAC_B; do
          CELLS+=("${model}|${train}|${test}|${test_time}|${mix}|${typed}|${seed}")
        done
      else
        CELLS+=("${model}|${train}|${test}|${test_time}|0.00|${typed}|${seed}")
      fi
    done
  done
}

if [[ "$SWEEP" == "OOA" || "$SWEEP" == "ALL" ]]; then
  # Three pop-pairs from OutOfAfrica_3G09
  add_cells OutOfAfrica_3G09      CEU YRI 0
  add_cells OutOfAfrica_3G09      CEU CHB 0
  add_cells OutOfAfrica_3G09      YRI CHB 0
fi
if [[ "$SWEEP" == "ADM-tr" || "$SWEEP" == "ALL" ]]; then
  # EUR → ADMIX with mix_B fraction of ADMIX in training
  add_cells AmericanAdmixture_4B18  EUR ADMIX 0
fi
if [[ "$SWEEP" == "ADM-te" || "$SWEEP" == "ALL" ]]; then
  # EUR → ADMIX, pure transfer (mix=0): admixed-test focus
  add_cells AmericanAdmixture_4B18  EUR ADMIX 0 no
fi
if [[ "$SWEEP" == "ARC" || "$SWEEP" == "ALL" ]]; then
  # CEU → Neanderthal@t=NEAN_TIME (no mix — train is always CEU)
  add_cells OutOfAfricaArchaicAdmixture_5R19  CEU Neanderthal "$NEAN_TIME" no
fi

N_CELLS=${#CELLS[@]}
if [ "$N_CELLS" -eq 0 ]; then
  echo "ERROR: SWEEP=$SWEEP produced no cells. Valid values: OOA ADM-tr ADM-te ARC ALL" >&2
  exit 1
fi

echo "Block 1.1 v3 stdpopsim portability driver"
echo "  base   : $BASE"
echo "  sweep  : $SWEEP  ($N_CELLS cells)"
echo "  scale  : n_train=$N_TRAIN n_test=$N_TEST L=$LENGTH"
echo "  trait  : h2_snp=$H2_SNP h2_arg=$H2_ARG causal_frac=$CAUSAL_FRAC"
echo "  python : $PY"
echo "  sargent: $SARGENT $SARGENT_EXTRA"
echo

t_start=$(date +%s)
cell_n=0
skipped=0

for entry in "${CELLS[@]}"; do
  IFS='|' read -r model train_pop test_pop test_time mix typed seed <<< "$entry"

  cell_n=$((cell_n + 1))

  # Short model tag for directory name
  case "$model" in
    OutOfAfrica_3G09)                      mtag="OOA3G09" ;;
    AmericanAdmixture_4B18)                mtag="ADM4B18" ;;
    OutOfAfricaArchaicAdmixture_5R19)      mtag="ARC5R19" ;;
    *)                                     mtag="$model"  ;;
  esac

  mtag_full="${mtag}_${train_pop}x${test_pop}"
  [ "$test_time" != "0" ] && [ "$test_time" != "0.0" ] \
    && mtag_full="${mtag_full}@t${test_time}"

  mix_fmt=$(printf "%.2f" "$mix")
  typed_fmt=$(printf "%.2f" "$typed")
  dir="$BASE/${mtag_full}_mix${mix_fmt}_typed${typed_fmt}_s${seed}"
  mkdir -p "$dir"

  if [ -s "$dir/est.predict.json" ]; then
    echo "[$cell_n/$N_CELLS] SKIP $dir"
    skipped=$((skipped + 1))
    continue
  fi

  echo "[$cell_n/$N_CELLS] $mtag_full mix=$mix_fmt typed=$typed_fmt seed=$seed"
  cstart=$(date +%s)

  # Build optional test-pop-time arg
  time_arg=""
  if [ "$test_time" != "0" ] && [ "$test_time" != "0.0" ]; then
    time_arg="--test-pop-time $test_time"
  fi

  # 1. simulate
  # shellcheck disable=SC2086
  $PY "$SIM" \
      --demography stdpopsim \
      --stdpopsim-model "$model" \
      --train-pop "$train_pop" \
      --test-pop  "$test_pop" \
      $time_arg \
      --mix-fraction-B "$mix" \
      --n-indiv "$N_TRAIN" --n-test "$N_TEST" \
      --length "$LENGTH" \
      --h2-snp "$H2_SNP" --h2-arg "$H2_ARG" \
      --causal-fraction "$CAUSAL_FRAC" \
      --bed-coverage "$typed" \
      --seed "$seed" --out "$dir/case" \
      > "$dir/sim.log" 2>&1 \
    || { echo "  sim FAILED  →  $dir/sim.log"; tail -6 "$dir/sim.log"; continue; }

  # 2. sargent
  # shellcheck disable=SC2086
  $SARGENT --silent --nprobe 10 --outroot "$dir/est" $SARGENT_EXTRA \
      "$dir/case.phen" "$dir/case.trees" "$dir/case" \
      > "$dir/sargent.log" 2>&1 \
    || { echo "  sargent FAILED  →  $dir/sargent.log"; tail -6 "$dir/sargent.log"; continue; }

  # 3. predict
  $PY "$PRED" \
      --case-prefix "$dir/case" --sargent-prefix "$dir/est" \
      > "$dir/predict.log" 2>&1 \
    || { echo "  predict FAILED  →  $dir/predict.log"; tail -6 "$dir/predict.log"; continue; }

  cend=$(date +%s)
  celapsed=$(( cend - cstart ))
  total_elapsed=$(( cend - t_start ))
  done_cells=$(( cell_n - skipped ))
  avg=$(( done_cells > 0 ? total_elapsed / done_cells : 0 ))
  eta=$(( (N_CELLS - cell_n) * avg ))
  printf "    done in %ds  (avg %ds/cell, ETA %dm)\n" \
    "$celapsed" "$avg" "$(( eta / 60 ))"
done

echo
echo "All cells processed (${cell_n} total, ${skipped} skipped)."

# ---- aggregate -----------------------------------------------------------------
echo "Aggregating into $BASE/portability_table_v3.md ..."

$PY - "$BASE" "$SWEEP" <<'PY' | tee "$BASE/portability_table_v3.md"
import json, os, sys, glob, math, re
import numpy as np

base  = sys.argv[1]
sweep = sys.argv[2] if len(sys.argv) > 2 else "ALL"

rows = []
for jf in sorted(glob.glob(os.path.join(base, "*_s*", "est.predict.json"))):
    d    = json.load(open(jf))
    cell = os.path.basename(os.path.dirname(jf))

    # Parse cell name: {model}_{trainpop}x{testpop}[@t{time}]_mix{m}_typed{t}_s{seed}
    pat = (r"([^_]+)_([^x]+)x([^_@]+)(?:@t(\d+))?_mix([\d.]+)_typed([\d.]+)_s(\d+)")
    m = re.match(pat, cell)
    if not m:
        continue
    mtag, train_pop, test_pop = m.group(1), m.group(2), m.group(3)
    test_time  = int(m.group(4)) if m.group(4) else 0
    mix, typed, seed = float(m.group(5)), float(m.group(6)), int(m.group(7))

    pjf    = os.path.join(os.path.dirname(jf), "case.params.json")
    params = json.load(open(pjf)) if os.path.exists(pjf) else {}

    def r2c(arm):
        a = d.get(arm, {})
        if isinstance(a.get("r2_calib"), (int, float)):
            return a["r2_calib"]
        rr = a.get("r")
        return rr * rr if isinstance(rr, (int, float)) else None

    rows.append({
        "cell": cell, "mtag": mtag,
        "train_pop": train_pop, "test_pop": test_pop,
        "test_time": test_time,
        "mix": mix, "typed": typed, "seed": seed,
        "fst":          params.get("fst_train_test"),
        "n_train_A":    params.get("n_train_A"),
        "n_train_B":    params.get("n_train_B"),
        "h2_snp_real":  params.get("h2_snp_realised"),
        "h2_arg_real":  params.get("h2_arg_realised"),
        "bayesR_r2c":   r2c("bayesR"),
        "gblup_r2c":    r2c("gblup"),
        "argunif_r2c":  r2c("arg_blup_uniform"),
        "sargent_r2c":  r2c("sargent"),
    })

def fnum(x, fmt="{:.4f}"):
    return fmt.format(x) if isinstance(x, (int, float)) and not math.isnan(float(x)) else "—"
def amean(rs, k):
    v = [r[k] for r in rs if isinstance(r[k], (int, float))]
    return float(np.mean(v)) if v else float("nan")

print(f"# Block 1.1 v3 — Realistic demographies ({sweep})\n")
print(f"Aggregated {len(rows)} cells from `{base}`.\n")

# Per-cell table
print("## Per-cell results (r²_calib)\n")
print("| model | train→test | F_ST | mix_B | typed | seed | BayesR | GBLUP | ARG-uni | SARGENT |")
print("|---|---|---|---|---|---|---|---|---|---|")
for r in rows:
    pop_pair = f"{r['train_pop']}→{r['test_pop']}"
    if r['test_time']:
        pop_pair += f"@t{r['test_time']}"
    print(f"| {r['mtag']} | {pop_pair} | {fnum(r['fst'])} "
          f"| {r['mix']:.2f} | {r['typed']:.2f} | {r['seed']} "
          f"| {fnum(r['bayesR_r2c'])} | {fnum(r['gblup_r2c'])} "
          f"| {fnum(r['argunif_r2c'])} | {fnum(r['sargent_r2c'])} |")

# Summary by (model, pop-pair, mix, typed) over seeds
print("\n## Summary by (model × pop-pair × mix × typed) — mean over seeds\n")
print("| model | train→test | mean F_ST | mix_B | typed | n | BayesR | GBLUP | ARG-uni | SARGENT | ARG-uni−GBLUP |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
keys = sorted({(r["mtag"], r["train_pop"], r["test_pop"], r["test_time"],
                r["mix"], r["typed"])
               for r in rows})
for (mt, tr, te, tt, mix, typed) in keys:
    sub = [r for r in rows
           if r["mtag"] == mt and r["train_pop"] == tr and r["test_pop"] == te
           and r["test_time"] == tt and r["mix"] == mix and r["typed"] == typed]
    pop_pair = f"{tr}→{te}" + (f"@t{tt}" if tt else "")

    fst_vals = [r["fst"] for r in sub if isinstance(r["fst"], (int, float))]
    fst_mean = fnum(float(np.mean(fst_vals))) if fst_vals else "—"

    b  = amean(sub, "bayesR_r2c");  g  = amean(sub, "gblup_r2c")
    a  = amean(sub, "argunif_r2c"); s_ = amean(sub, "sargent_r2c")
    print(f"| {mt} | {pop_pair} | {fst_mean} | {mix:.2f} | {typed:.2f} | {len(sub)} "
          f"| {fnum(b)} | {fnum(g)} | {fnum(a)} | {fnum(s_)} | {fnum(a-g)} |")

# Framework vs typed-SNP headline by F_ST
print("\n## Framework advantage over typed-SNP methods by F_ST (mean over mix × typed × seeds)\n")
fst_bins = {}
for r in rows:
    fst_val = r["fst"]
    if not isinstance(fst_val, (int, float)):
        continue
    # bin by 0.05 steps
    b_key = round(round(fst_val / 0.05) * 0.05, 2)
    fst_bins.setdefault(b_key, []).append(r)

if fst_bins:
    print("| F_ST (bin) | n | BayesR | GBLUP | ARG-uni | SARGENT | ARG-uni−GBLUP | ARG-uni−BayesR |")
    print("|---|---|---|---|---|---|---|---|")
    for bk in sorted(fst_bins):
        rs = fst_bins[bk]
        b  = amean(rs, "bayesR_r2c"); g  = amean(rs, "gblup_r2c")
        a  = amean(rs, "argunif_r2c"); s_ = amean(rs, "sargent_r2c")
        print(f"| ≈{bk:.2f} | {len(rs)} "
              f"| {fnum(b)} | {fnum(g)} | {fnum(a)} | {fnum(s_)} "
              f"| {fnum(a-g)} | {fnum(a-b)} |")
PY

t_end=$(date +%s)
echo
echo "Total wall time: $(( (t_end - t_start) / 60 )) min"
echo "Wrote $BASE/portability_table_v3.md"
