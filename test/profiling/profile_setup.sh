#!/usr/bin/env bash
# profile_setup.sh — Phase 4 profiling driver
#
# Simulates a fixed medium-scale case and profiles sargent under samply
# for joint-EM and sequential modes.
#
# Prerequisites:
#   - mamba env "sargent" activated (for simulate.py dependencies)
#   - samply installed (brew install samply)
#   - sargent built with -DCMAKE_BUILD_TYPE=RelWithDebInfo
#
# Output profile files (not committed):
#   ~/test-sargent/profiling/profile_joint_em.json.gz
#   ~/test-sargent/profiling/profile_sequential.json.gz

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SARGENT="$REPO_ROOT/build/src/progs/sargent"
SIMULATE="$REPO_ROOT/scripts/simulate.py"
DATA_DIR="$HOME/test-sargent/profiling"
CASE="$DATA_DIR/case"

mkdir -p "$DATA_DIR"

# --- Step 1: simulate medium-scale case (n=300, L=300kb, balanced h², seed=42) ---
if [[ ! -f "$CASE.trees" ]]; then
    echo "==> Simulating case..."
    conda run -n sargent python3 "$SIMULATE" \
        --n-indiv 300 --n-test 0 --length 300000 \
        --h2-snp 0.4 --h2-arg 0.4 \
        --bed-coverage 1.0 --seed 42 \
        --out "$CASE"
else
    echo "==> Case already exists, skipping simulation."
fi

# --- Step 2: profile joint EM (default, --em 50) ---
echo "==> Profiling joint EM (--em 50)..."
samply record --save-only \
    -o "$DATA_DIR/profile_joint_em.json.gz" \
    "$SARGENT" \
        "$CASE.phen" "$CASE.trees" "$CASE" \
        -o "$DATA_DIR/run_joint_em" \
        --em 50
echo "    Profile saved: $DATA_DIR/profile_joint_em.json.gz"

# --- Step 3: profile sequential mode ---
echo "==> Profiling sequential (--sequential)..."
samply record --save-only \
    -o "$DATA_DIR/profile_sequential.json.gz" \
    "$SARGENT" \
        "$CASE.phen" "$CASE.trees" "$CASE" \
        -o "$DATA_DIR/run_sequential" \
        --sequential
echo "    Profile saved: $DATA_DIR/profile_sequential.json.gz"

echo ""
echo "==> Done. Load profiles in Firefox Profiler:"
echo "    https://profiler.firefox.com/"
echo "    File → Load from file → select the .json.gz"
