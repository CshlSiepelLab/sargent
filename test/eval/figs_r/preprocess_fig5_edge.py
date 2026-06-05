#!/usr/bin/env python3
"""
preprocess_fig5_edge.py — Fig 5 Panel B from per-edge architectures
(Block 2.7.1).

For each cell under $BASE (fig5_panelB_edges):
  - Run edge_h2.py if est.edge_h2.tsv not already present.
  - Read per-edge TSV; bin edges by midpoint age (recent / middle /
    ancient).
  - Sum Δh²_total per (architecture, observed bin, seed).

Output a single tidy TSV the R plotter can consume:
  _cache/fig5_panelB_edges.tsv
    columns: true_arch, observed_bin, seed, dh2_total, dh2_arg, dh2_snp,
             dh2_total_norm  (normalised so the 3 bins sum to 1.0 per seed)

Also writes truth-side per-bin sums for an idealised-target overlay:
  _cache/fig5_panelB_edges_truth.tsv

CLI:
  python preprocess_fig5_edge.py [--base BASE] [--rerun-edge-h2]
"""
import argparse
import csv
import glob
import os
import subprocess
import sys

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
EVAL_DIR   = os.path.dirname(SCRIPT_DIR)
CACHE      = os.path.join(SCRIPT_DIR, "_cache")
os.makedirs(CACHE, exist_ok=True)

AGE_EDGES  = [0.0, 5_000.0, 20_000.0, np.inf]
AGE_LABELS = ["recent", "middle", "ancient"]


def bin_age(mid):
    for i, label in enumerate(AGE_LABELS):
        lo, hi = AGE_EDGES[i], AGE_EDGES[i + 1]
        if lo <= mid < hi:
            return label
    return AGE_LABELS[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=os.path.expanduser(
        "~/test-sargent/fig5_panelB_edges"))
    ap.add_argument("--rerun-edge-h2", action="store_true",
                    help="re-run edge_h2.py even if est.edge_h2.tsv exists.")
    args = ap.parse_args()

    cells = sorted(glob.glob(os.path.join(args.base, "*_s*")))
    if not cells:
        print(f"No cells found under {args.base}", file=sys.stderr)
        sys.exit(1)

    edge_h2_py = os.path.join(EVAL_DIR, "edge_h2.py")
    if not os.path.exists(edge_h2_py):
        print(f"edge_h2.py not found at {edge_h2_py}", file=sys.stderr)
        sys.exit(1)

    out_rows       = []
    out_truth_rows = []

    for cell in cells:
        name = os.path.basename(cell)
        # name format: age_<tertile>_edges_s<seed>
        # tertile ∈ {recent, middle, ancient}
        try:
            arch = name.rsplit("_s", 1)[0]   # 'age_recent_edges'
            seed = int(name.rsplit("_s", 1)[1])
        except (ValueError, IndexError):
            continue

        edge_tsv = os.path.join(cell, "est.edge_h2.tsv")
        case     = os.path.join(cell, "case")
        est      = os.path.join(cell, "est")

        if args.rerun_edge_h2 or not os.path.exists(edge_tsv):
            print(f"[edge_h2] {name} …", flush=True)
            r = subprocess.run(
                ["mamba", "run", "-n", "sargent", "python", edge_h2_py,
                 "--case-prefix", case, "--sargent-prefix", est,
                 "--no-bayesR-fit",
                 "--out", edge_tsv],
                stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            if r.returncode != 0:
                print(f"  FAIL: {r.stderr.decode()[:200]}", file=sys.stderr)
                continue

        # Read per-edge data and aggregate
        per_bin = {b: {"dh2_total": 0.0, "dh2_arg": 0.0, "dh2_snp": 0.0}
                   for b in AGE_LABELS}
        with open(edge_tsv) as f:
            rd = csv.DictReader(f, delimiter="\t")
            for row in rd:
                age = 0.5 * (float(row["child_age"]) +
                             float(row["child_age"]) + 2 * float(row["length"]) / 2.0)
                # midpoint = child_age + length/2 (parent_age = child_age + length)
                age = float(row["child_age"]) + 0.5 * float(row["length"])
                b   = bin_age(age)
                per_bin[b]["dh2_total"] += float(row["dh2_total"])
                per_bin[b]["dh2_arg"]   += float(row["dh2_arg"])
                per_bin[b]["dh2_snp"]   += float(row["dh2_snp"])

        # Normalise across bins
        total = sum(per_bin[b]["dh2_total"] for b in AGE_LABELS)
        for b in AGE_LABELS:
            out_rows.append({
                "true_arch":         arch,
                "observed_bin":      b,
                "seed":              seed,
                "dh2_total":         per_bin[b]["dh2_total"],
                "dh2_arg":           per_bin[b]["dh2_arg"],
                "dh2_snp":           per_bin[b]["dh2_snp"],
                "dh2_total_norm":    (per_bin[b]["dh2_total"] / total) if total > 0 else 0.0,
            })

        # Truth-side per-bin sums (from .true_ws2_per_edge.tsv)
        truth_path = os.path.join(cell, "case.true_ws2_per_edge.tsv")
        if os.path.exists(truth_path):
            truth_data = np.loadtxt(truth_path, skiprows=1, ndmin=2)
            mids       = 0.5 * (truth_data[:, 2] + truth_data[:, 3])
            ws2_e      = truth_data[:, 5]
            ell        = truth_data[:, 4]
            truth_per_bin = {b: 0.0 for b in AGE_LABELS}
            for i, m in enumerate(mids):
                truth_per_bin[bin_age(float(m))] += ws2_e[i] * ell[i]
            t_total = sum(truth_per_bin.values())
            for b in AGE_LABELS:
                out_truth_rows.append({
                    "true_arch":      arch,
                    "observed_bin":   b,
                    "seed":           seed,
                    "truth_weighted": truth_per_bin[b],
                    "truth_norm":     (truth_per_bin[b] / t_total) if t_total > 0 else 0.0,
                })

    # Write per-cell aggregated TSVs
    out_path = os.path.join(CACHE, "fig5_panelB_edges.tsv")
    fields  = ["true_arch", "observed_bin", "seed", "dh2_total",
               "dh2_arg", "dh2_snp", "dh2_total_norm"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for row in out_rows:
            w.writerow(row)
    print(f"  wrote {out_path} ({len(out_rows)} rows)")

    truth_path = os.path.join(CACHE, "fig5_panelB_edges_truth.tsv")
    tfields = ["true_arch", "observed_bin", "seed",
               "truth_weighted", "truth_norm"]
    with open(truth_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=tfields, delimiter="\t")
        w.writeheader()
        for row in out_truth_rows:
            w.writerow(row)
    print(f"  wrote {truth_path} ({len(out_truth_rows)} rows)")

    # Quick summary print
    if out_rows:
        print()
        print("Estimator-side mean Δh²_total_norm by (true_arch, observed_bin):")
        archs = sorted({r["true_arch"] for r in out_rows})
        labels = AGE_LABELS
        print(f"  {'true':>22s}  " + "  ".join(f"{b:>8s}" for b in labels))
        for a in archs:
            row_means = []
            for b in labels:
                vals = [r["dh2_total_norm"] for r in out_rows
                        if r["true_arch"] == a and r["observed_bin"] == b]
                row_means.append(float(np.mean(vals)) if vals else float("nan"))
            print(f"  {a:>22s}  " + "  ".join(f"{v:>8.3f}" for v in row_means))


if __name__ == "__main__":
    main()
