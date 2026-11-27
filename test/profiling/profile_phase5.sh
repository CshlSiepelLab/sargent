#!/usr/bin/env bash
# profile_phase5.sh — Phase 5 chunk 4 re-profile driver.
#
# Phase 4 profiled n=300, L=300kb, BEFORE chunks 2-4 wired incremental
# matvec into production.  This script re-profiles AFTER chunk 4 on a
# representative medium-scale case (n=1000, L=1Mb) so we can see where
# wall time actually goes now that tsk_tree_postorder/preorder are out
# of the hot path.
#
# Prerequisites (desktop):
#   - mamba env "sargent" activated (for simulate.py dependencies)
#   - samply installed (brew install samply)
#   - sargent built with -DCMAKE_BUILD_TYPE=RelWithDebInfo from
#     branch lap/phase5-chunk4 (PR #8) or main once #8 lands
#
# Output profile files (not committed):
#   ~/test-sargent/profiling/phase5/profile_joint_em.json.gz
#   ~/test-sargent/profiling/phase5/profile_sequential.json.gz
#
# Load in Firefox Profiler: https://profiler.firefox.com/

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SARGENT="$REPO_ROOT/build/src/progs/sargent"
SIMULATE="$REPO_ROOT/scripts/simulate.py"
DATA_DIR="$HOME/test-sargent/profiling/phase5"
CASE="$DATA_DIR/medium"

mkdir -p "$DATA_DIR"

if [[ ! -x "$SARGENT" ]]; then
    echo "ERROR: sargent not built.  Build with:" >&2
    echo "  mamba run -n sargent cmake -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo" >&2
    echo "  mamba run -n sargent cmake --build build" >&2
    exit 2
fi

# --- Step 1: simulate medium-scale case (n=1000, L=1Mb, balanced h², seed=42) ---
if [[ ! -f "$CASE.trees" ]]; then
    echo "==> Simulating medium case (n=1000, L=1Mb, h2_snp=h2_arg=0.4, seed=42)..."
    conda run -n sargent python3 "$SIMULATE" \
        --n-indiv 1000 --n-test 0 --length 1000000 \
        --h2-snp 0.4 --h2-arg 0.4 \
        --bed-coverage 1.0 --seed 42 \
        --out "$CASE"
else
    echo "==> Case already exists at $CASE.*, skipping simulation."
fi

# Sanity: print case stats
echo ""
echo "Case parameters:"
conda run -n sargent python3 -c "
import json
p = json.load(open('$CASE.params.json'))
for k in ['n_indiv','n_trees','n_variants','length','h2_snp_target','h2_arg_target','sigma2']:
    print(f'  {k:20s} = {p[k]}')
"

# --- Step 2: profile sequential mode (faster; ~30-40s on n=1000) ---
echo ""
echo "==> Profiling sequential (--sequential)..."
time samply record --save-only \
    -o "$DATA_DIR/profile_sequential.json.gz" \
    "$SARGENT" \
        --silent --sequential --nprobe 10 \
        --outroot "$DATA_DIR/run_sequential" \
        "$CASE.phen" "$CASE.trees" "$CASE"
echo "    Profile saved: $DATA_DIR/profile_sequential.json.gz"

# --- Step 3: profile joint EM 50 (slower; ~1.5-3 min on n=1000) ---
echo ""
echo "==> Profiling joint EM (--em 50)..."
time samply record --save-only \
    -o "$DATA_DIR/profile_joint_em.json.gz" \
    "$SARGENT" \
        --silent --em 50 --nprobe 10 \
        --outroot "$DATA_DIR/run_joint_em" \
        "$CASE.phen" "$CASE.trees" "$CASE"
echo "    Profile saved: $DATA_DIR/profile_joint_em.json.gz"

echo ""
echo "==> Done.  Load profiles in Firefox Profiler:"
echo "    https://profiler.firefox.com/"
echo "    File → Load from file → select the .json.gz"
echo ""
echo "Expected: tskit traversal cost dropped from ~70% to <5%; new"
echo "hotspots likely include tsk_tree_next (edge_diff bookkeeping),"
echo "the per-sample walk-up in our_parent[], geno_read_variant in"
echo "the BayesR sweep, or cg_solve."
