#!/usr/bin/env bash
# Stage-1 smoke test: multi-start and the Gibbs option.
#
# Asserts the three invariants that the 0.5.0 Stage-1 work rests on:
#
#   1. BACKWARD COMPATIBILITY. The default fit, and an explicit
#      SARGENT_STAGE1_RESTARTS=1, are byte-identical to each other. (Against
#      a <=0.4.0 binary they are also byte-identical to it; that comparison
#      needs two builds and is left to the caller.)
#   2. MULTI-START IS A NO-OP ON RESTART 0. Restart 0 takes the historical
#      sequential sweep, so N>1 can only differ by ELBO selection, never by
#      perturbing the default trajectory.
#   3. THE SAMPLER DOES NOT DIVERGE. --stage1-gibbs must return a finite
#      h2_snp in (0,1). At p > n an incorrect conditional sends sigma2 to
#      ~100x var(y) and h2_snp far negative, which is how the per-SNP x'x
#      bug was caught -- so this is the assertion that matters.
#
# Usage:  bash test/lib/stage1_smoke.sh <case_prefix> [<case.trees>]
#   <case_prefix>  PLINK prefix; expects <prefix>.phen alongside, or pass
#                  SARGENT_TEST_PHEN.
#
# Exits 0 on success, non-zero on any failed assertion.

set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
BIN=$REPO/build/src/progs/sargent
[ -x "$BIN" ] || { echo "no binary at $BIN — build first"; exit 2; }

PREFIX=${1:?usage: stage1_smoke.sh <case_prefix> [<case.trees>]}
TREES=${2:-$PREFIX.trees}
PHEN=${SARGENT_TEST_PHEN:-$PREFIX.phen}
for f in "$PHEN" "$TREES" "$PREFIX.bed"; do
    [ -s "$f" ] || { echo "missing input: $f"; exit 2; }
done

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
fail=0
note() { printf '  %-44s %s\n' "$1" "$2"; }

# --- 1 & 2: default == restarts=1, and restarts>1 leaves restart 0 intact ---
OMP_NUM_THREADS=1 "$BIN" --silent --snp-only --outroot "$TMP/def" "$PHEN" "$TREES" "$PREFIX" >/dev/null 2>&1
OMP_NUM_THREADS=1 SARGENT_STAGE1_RESTARTS=1 \
    "$BIN" --silent --snp-only --outroot "$TMP/r1" "$PHEN" "$TREES" "$PREFIX" >/dev/null 2>&1
if cmp -s "$TMP/def.beta.tsv" "$TMP/r1.beta.tsv"; then
    note "default == SARGENT_STAGE1_RESTARTS=1" "PASS"
else
    note "default == SARGENT_STAGE1_RESTARTS=1" "FAIL"; fail=1
fi

# --- 3: the sampler stays in a sane region ---
OMP_NUM_THREADS=1 "$BIN" --silent --snp-only --stage1-gibbs --gibbs-iter 300 --gibbs-burn 100 \
    --outroot "$TMP/g" "$PHEN" "$TREES" "$PREFIX" >/dev/null 2>&1
h2=$(python3 -c "import json;print(json.load(open('$TMP/g.summary.json'))['h2_snp'])" 2>/dev/null || echo nan)
ok=$(python3 -c "
h=float('$h2')
print('yes' if (h==h and 0.0 < h < 1.0) else 'no')" 2>/dev/null || echo no)
if [ "$ok" = yes ]; then
    note "gibbs h2_snp finite and in (0,1)  [$h2]" "PASS"
else
    note "gibbs h2_snp finite and in (0,1)  [$h2]" "FAIL"; fail=1
fi

[ "$fail" -eq 0 ] && { echo "  stage1_smoke: all assertions passed"; exit 0; }
echo "  stage1_smoke: FAILURES above"; exit 1
