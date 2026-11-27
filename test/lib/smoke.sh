#!/usr/bin/env bash
# Fast smoke test for the tree-matvec paranoia harness.
#
# Runs check_tree_matvec with light parameters (10 probes × 20 trees) on a
# committed debug-scale tree sequence.  Designed for routine iteration during
# Phase 5 development; finishes in <10 s on a laptop.
#
# Full validation (100 probes × all trees) takes ~4 min and is reserved for
# milestone runs:
#   build/test/lib/check_tree_matvec <case.trees>
#
# Usage:
#   bash test/lib/smoke.sh
#
# Exits 0 on success, non-zero on harness failure or missing binary.

set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
BIN=${BIN:-$REPO/build/test/lib/check_tree_matvec}

# Smoke case: a *.trees file passed as $1, or $SARGENT_SMOKE_CASE.
CASE=${1:-${SARGENT_SMOKE_CASE:-}}
if [[ -z $CASE ]]; then
    echo "usage: $0 case.trees   (or set SARGENT_SMOKE_CASE)"
    exit 2
fi

if [[ ! -x $BIN ]]; then
    echo "ERROR: $BIN not built.  Run:"
    echo "    mamba run -n sargent cmake --build build"
    exit 2
fi

if [[ ! -e $CASE ]]; then
    echo "ERROR: tree sequence not found: $CASE"
    exit 2
fi

# Light load: 10 probes × first 20 trees.  Tighter tol (1e-10) since fewer
# accumulations → less worst-case roundoff.
exec "$BIN" "$CASE" 10 20 1e-10 42
