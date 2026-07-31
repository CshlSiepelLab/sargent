#!/usr/bin/env bash
# native_predict_test.sh — merge gate for the native --predict path.
#
# Generates a small case, runs `sargent --predict` (native JSON) and the
# reference analysis/predict.py, and asserts every metric agrees to < 1e-6.
# Also checks that omitting --predict leaves the fit outputs untouched (T5).
#
# Env:
#   BASE      work dir            (default: mktemp under TMPDIR)
#   PAPER     paper repo w/ predict.py  (default: ~/git-src/sargent-paper)
#   TOL       numeric tolerance   (default: 1e-6)
#   KEEP      1 = keep BASE       (default: 0)
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
SARGENT=${SARGENT:-$ROOT/build/src/progs/sargent}
SIM=${SIM:-$ROOT/scripts/simulate.py}
PAPER=${PAPER:-$HOME/git-src/sargent-paper}
PRED=${PRED:-$PAPER/analysis/predict.py}
TOL=${TOL:-1e-6}
BASE=${BASE:-$(mktemp -d "${TMPDIR:-/tmp}/native_predict.XXXXXX")}
KEEP=${KEEP:-0}

if [ -n "${PYBIN:-}" ]; then PY="$PYBIN"
elif command -v mamba >/dev/null 2>&1 && mamba env list 2>/dev/null | grep -q '^sargent'; then
  PY="mamba run -n sargent python"
else PY="python"; fi

[ -x "$SARGENT" ] || { echo "ERROR: sargent binary not found at $SARGENT (build first)"; exit 1; }
[ -f "$PRED" ]    || { echo "ERROR: predict.py not found at $PRED (set PAPER=)"; exit 1; }

echo "native_predict_test"
echo "  base    : $BASE"
echo "  sargent : $SARGENT"
echo "  predict : $PRED"
mkdir -p "$BASE"

# --- 1. Simulate a small case with a test split ----------------------------
echo "[1/4] simulate ..."
$PY "$SIM" --n-indiv 150 --n-test 60 --length 300000 \
  --h2-snp 0.3 --h2-arg 0.3 --n-causal 40 --bed-coverage 0.7 \
  --seed 7 --out "$BASE/case" >"$BASE/sim.log" 2>&1

# --- 2. sargent fit + native predict ---------------------------------------
echo "[2/4] sargent --predict ..."
$SARGENT --predict -o "$BASE/est" \
  "$BASE/case.phen" "$BASE/case.trees" "$BASE/case" >"$BASE/sargent.log" 2>&1
test -f "$BASE/est.predict.json" || { echo "FAIL: no native est.predict.json"; cat "$BASE/sargent.log"; exit 1; }
cp "$BASE/est.beta.tsv" "$BASE/est.beta.ref.tsv"
cp "$BASE/est.w2.tsv"   "$BASE/est.w2.ref.tsv"

# --- 3. reference predict.py -----------------------------------------------
echo "[3/4] predict.py (reference) ..."
$PY "$PRED" --case-prefix "$BASE/case" --sargent-prefix "$BASE/est" \
  --out "$BASE/ref.predict.json" >"$BASE/predpy.log" 2>&1

# --- 4. compare -------------------------------------------------------------
echo "[4/4] compare (tol=$TOL) ..."
$PY - "$BASE/est.predict.json" "$BASE/ref.predict.json" "$TOL" <<'PYEOF'
import json, sys, math
native, ref, tol = sys.argv[1], sys.argv[2], float(sys.argv[3])
A = json.load(open(native)); B = json.load(open(ref))

def flat(d, pre=""):
    out = {}
    for k, v in d.items():
        key = f"{pre}{k}"
        if isinstance(v, dict): out.update(flat(v, key + "."))
        else: out[key] = v
    return out

fa, fb = flat(A), flat(B)
keys = sorted(set(fa) & set(fb))
worst = 0.0; worst_k = None; fails = []
print(f"{'field':28s} {'native':>14s} {'ref':>14s} {'absdiff':>12s}")
for k in keys:
    a, b = fa[k], fb[k]
    if not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
        if a != b: fails.append((k, a, b, 'non-numeric mismatch'))
        continue
    if isinstance(a, float) and math.isnan(a) and isinstance(b, float) and math.isnan(b):
        d = 0.0
    else:
        d = abs(a - b)
    rel = d / max(1.0, abs(b))
    flag = "" if rel <= tol else "  <-- FAIL"
    if flag: fails.append((k, a, b, f"rel={rel:.2e}"))
    if rel > worst: worst, worst_k = rel, k
    print(f"{k:28s} {a:14.7g} {b:14.7g} {d:12.3e}{flag}")

only = sorted(set(fa) ^ set(fb))
if only: print("\nkeys only in one file:", only)
print(f"\nworst relative diff: {worst:.3e} at {worst_k}")
if fails:
    print(f"\nRESULT: FAIL ({len(fails)} field(s) exceed tol={tol})")
    sys.exit(1)
print(f"\nRESULT: PASS (all {len(keys)} shared metrics agree within {tol})")
PYEOF
CMP=$?

# --- T5: fit outputs unchanged without --predict ---------------------------
echo "[T5] no-op check (fit unchanged without --predict) ..."
$SARGENT -o "$BASE/est_nopred" \
  "$BASE/case.phen" "$BASE/case.trees" "$BASE/case" >"$BASE/sargent_nopred.log" 2>&1
if [ -f "$BASE/est_nopred.predict.json" ]; then
  echo "T5 FAIL: predict.json written without --predict"; exit 1; fi
if ! diff -q "$BASE/est.beta.ref.tsv" "$BASE/est_nopred.beta.tsv" >/dev/null; then
  echo "T5 FAIL: beta differs with/without --predict"; exit 1; fi
if ! diff -q "$BASE/est.w2.ref.tsv" "$BASE/est_nopred.w2.tsv" >/dev/null; then
  echo "T5 FAIL: w2 differs with/without --predict"; exit 1; fi
echo "T5 PASS: fit outputs identical, no predict.json without --predict"

[ "$KEEP" = "1" ] || rm -rf "$BASE"
exit $CMP
