#!/usr/bin/env python3
"""
preprocess_fig5.py — Build TSV caches for Fig 5 Panels B and C.

Panel B (WHEN axis): variance-capture curves stratified by allele-age bin
  (recent ≤5k, middle 5k–20k, ancient >20k generations) for Block 2.1 v2
  ARG-side cells at n=1000, 3 seeds.

Panel C (HOW axis): variance-capture curves by architecture (block-fraction
  bf=0.005/0.05/0.5) across Block 2.5.2 polygenicity sweep (30 cells).

Outputs → test/eval/figs_r/_cache/:
  fig5_panelB.tsv   columns: age_bin, kind, k_frac, captured_frac, seed, S_bin
  fig5_panelC.tsv   columns: bf, kind, k_frac, captured_frac, seed, demog, S_total
"""

import os, sys, csv
import numpy as np
import tskit

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
EVAL_DIR   = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, EVAL_DIR)

from arch_metrics import per_tree_mean_branch_age, variance_capture_curve
from predict import load_ws2

CACHE = os.path.join(SCRIPT_DIR, "_cache")
os.makedirs(CACHE, exist_ok=True)

# ---------------------------------------------------------------------------
# Age bin edges (generations)
# ---------------------------------------------------------------------------
AGE_EDGES  = [0.0, 5_000.0, 20_000.0, np.inf]
AGE_LABELS = ["recent", "middle", "ancient"]
N_CURVE_POINTS_B = 150   # linear k_frac grid for Panel B (linear x-axis)
N_CURVE_POINTS_C = 200   # log k_frac grid for Panel C (log x-axis)


def sample_curve_linear(k_arr, capt, S_bin, n_points):
    """Sample (k_frac, captured_frac) pairs on a linear k_frac grid."""
    rows = []
    if len(k_arr) == 0:
        return rows
    idx = np.unique(np.round(np.linspace(0, len(k_arr) - 1, n_points)).astype(int))
    for i in idx:
        rows.append((float(k_arr[i]) / S_bin, float(capt[i])))
    # guarantee endpoint
    if (float(k_arr[-1]) / S_bin, 1.0) not in rows:
        rows.append((1.0, 1.0))
    return rows


def sample_curve_log(k_arr, capt, S_total, n_points):
    """Sample (k_frac, captured_frac) pairs on a log k_frac grid."""
    rows = []
    if len(k_arr) == 0:
        return rows
    log_min = np.log10(1.0 / S_total)
    log_grid = np.logspace(log_min, 0, n_points)
    k_norm = np.array(k_arr, dtype=float) / S_total
    idx = np.unique(np.searchsorted(k_norm, log_grid).clip(0, len(k_arr) - 1))
    for i in idx:
        rows.append((float(k_arr[i]) / S_total, float(capt[i])))
    rows.append((1.0, 1.0))
    return rows


# ---------------------------------------------------------------------------
# Panel B — WHEN axis (age-stratified curves)
# ---------------------------------------------------------------------------

def process_panelB():
    # New data: age_* architecture sweep (fig5_panelB cells).
    # ws2 is assigned by the simulator to the target age tertile only, so we
    # compute the within-bin capture curves using the ACTUAL tree-sequence ages
    # (not the simulator's assignment tertiles) to reproduce what the estimator sees.
    base = os.path.expanduser("~/test-sargent/fig5_panelB")
    archs = ["age_recent", "age_middle", "age_ancient"]
    seeds = [1, 2, 3, 4, 5]
    out_path = os.path.join(CACHE, "fig5_panelB.tsv")

    writer_rows = []

    for arch in archs:
        for seed in seeds:
            cell = os.path.join(base, f"{arch}_s{seed}")
            ts_path   = os.path.join(cell, "case.trees")
            est_path  = os.path.join(cell, "est.w2.tsv")
            true_path = os.path.join(cell, "case.true_ws2.tsv")

            if not os.path.exists(ts_path):
                print(f"  [B] {arch} s{seed}: case.trees not found, skipping")
                continue

            print(f"  [B] {arch} s{seed}: loading ts (S=?)", end="", flush=True)
            ts = tskit.load(ts_path)
            S  = ts.num_trees
            print(f" S={S}", flush=True)

            print(f"  [B] {arch} s{seed}: computing ages …", flush=True)
            ages = per_tree_mean_branch_age(ts)

            sources = {}
            for kind, path in [("est", est_path), ("truth", true_path)]:
                if not os.path.exists(path):
                    continue
                w2 = load_ws2(path)
                if len(w2) == S:
                    sources[kind] = w2
                else:
                    print(f"  [B] {arch} s{seed}: {kind} length {len(w2)} != S={S}, skipping")

            for kind, w2 in sources.items():
                for bi, label in enumerate(AGE_LABELS):
                    lo = AGE_EDGES[bi]
                    hi = AGE_EDGES[bi + 1]
                    mask = (ages >= lo) & (ages < hi) if hi < np.inf else (ages >= lo)
                    valid_mask = mask & ~np.isnan(ages)
                    S_bin = int(valid_mask.sum())

                    if S_bin == 0:
                        print(f"  [B] {arch} s{seed} kind={kind} bin={label}: 0 trees")
                        continue

                    bin_w2 = w2[valid_mask]
                    k_arr, capt = variance_capture_curve(bin_w2)
                    pts = sample_curve_linear(k_arr, capt, S_bin, N_CURVE_POINTS_B)
                    for k_frac, captured_frac in pts:
                        writer_rows.append({
                            "arch": arch,
                            "age_bin": label,
                            "kind": kind,
                            "k_frac": f"{k_frac:.6f}",
                            "captured_frac": f"{captured_frac:.6f}",
                            "seed": seed,
                            "S_bin": S_bin,
                        })

            print(f"  [B] {arch} s{seed}: done", flush=True)

    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["arch", "age_bin", "kind", "k_frac",
                                          "captured_frac", "seed", "S_bin"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(writer_rows)
    print(f"Panel B → {out_path}  ({len(writer_rows)} rows)")


# ---------------------------------------------------------------------------
# Panel C — HOW axis (architecture sweep)
# ---------------------------------------------------------------------------

def process_panelC():
    base  = os.path.expanduser("~/test-sargent/block2_5_2")
    bfs   = ["0.005", "0.05", "0.5"]
    demogs = ["admixed", "panmictic"]
    seeds  = [1, 2, 3, 4, 5]
    out_path = os.path.join(CACHE, "fig5_panelC.tsv")

    writer_rows = []

    for bf in bfs:
        for demog in demogs:
            for seed in seeds:
                cell = os.path.join(base, f"bf{bf}_d{demog}_s{seed}")
                est_path    = os.path.join(cell, "est_fl05.w2.tsv")
                est_nfl_path = os.path.join(cell, "est.w2.tsv")
                true_path   = os.path.join(cell, "case.true_ws2.tsv")

                sources = {}
                for kind, path in [("est_fl", est_path), ("est", est_nfl_path),
                                   ("truth", true_path)]:
                    if not os.path.exists(path):
                        print(f"  [C] bf={bf} {demog} s={seed}: {kind} missing")
                        continue
                    sources[kind] = load_ws2(path)

                if not sources:
                    continue

                S_total = len(next(iter(sources.values())))

                for kind, w2 in sources.items():
                    k_arr, capt = variance_capture_curve(w2)
                    pts = sample_curve_log(k_arr, capt, S_total, N_CURVE_POINTS_C)
                    for k_frac, captured_frac in pts:
                        writer_rows.append({
                            "bf": bf,
                            "kind": kind,
                            "k_frac": f"{k_frac:.8f}",
                            "captured_frac": f"{captured_frac:.6f}",
                            "seed": seed,
                            "demog": demog,
                            "S_total": S_total,
                        })

    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["bf", "kind", "k_frac",
                                          "captured_frac", "seed", "demog", "S_total"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(writer_rows)
    print(f"Panel C → {out_path}  ({len(writer_rows)} rows)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", choices=["B", "C", "both"], default="both")
    args = ap.parse_args()

    if args.panel in ("B", "both"):
        print("=== Panel B (WHEN axis, age bins) ===")
        process_panelB()
    if args.panel in ("C", "both"):
        print("=== Panel C (HOW axis, architectures) ===")
        process_panelC()
