#!/usr/bin/env python3
"""
plot.py — Fig 4 generator (manuscript localization headline).

Three panels:

  Panel A: Manhattan-style estimated per-tree ŵ²_s along the genome,
           with the true active blocks shaded.  Single representative
           cell from Block 1.2 (sparse + admixed + contrast 1000).
  Panel B: Variance-capture curves cumulative-vs-top-k for two
           architectures (sparse vs diffuse), truth + estimator.
  Panel C: Localization AUC across (contrast × demography), bar groups.

Reads: ~/test-sargent/block1_2_v1/ + ~/test-sargent/phase10_admix/.
"""

import argparse
import glob
import json
import os
import re
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import tskit


# ---------- Panel A: Manhattan per-tree ŵ² ----------

def load_cell_ws2(case_prefix, est_prefix):
    """Load per-tree truth + estimated ŵ² plus the train tree-sequence
    breakpoints for x-axis positioning."""
    truth_path = case_prefix + ".true_ws2.tsv"  # train-tree-indexed
    if not os.path.exists(truth_path):
        truth_path = case_prefix + ".full.true_ws2.tsv"
    w2_truth = np.loadtxt(truth_path, skiprows=1)[:, -1]
    w2_est   = np.loadtxt(est_prefix + ".w2.tsv", skiprows=1)[:, -1]
    ts_train = tskit.load(case_prefix + ".trees")
    # Tree midpoints in genomic coordinates for x-axis
    breaks = list(ts_train.breakpoints())
    mids = np.array([0.5 * (breaks[i] + breaks[i+1])
                     for i in range(ts_train.num_trees)])
    return w2_truth, w2_est, mids


def fig4_panel_a(ax, cell_dir):
    case = os.path.join(cell_dir, "case")
    est  = os.path.join(cell_dir, "est")
    if not os.path.exists(est + ".w2.tsv"):
        ax.text(0.5, 0.5, f"No data:\n{est}.w2.tsv missing",
                ha="center", va="center", transform=ax.transAxes)
        return
    w2_truth, w2_est, mids = load_cell_ws2(case, est)
    # Map both arrays to the same length (the train tree count) if needed
    n = min(len(w2_truth), len(w2_est), len(mids))
    w2_truth = w2_truth[:n]; w2_est = w2_est[:n]; mids = mids[:n]

    # Active blocks (truth has elevated ŵ² above the floor); shade their
    # extent on the plot.
    is_active = w2_truth > w2_truth.min() + 1e-12
    # Build contiguous blocks
    blocks = []
    in_block = False
    lo = 0
    for i, a in enumerate(is_active):
        if a and not in_block:
            in_block = True
            lo = i
        elif not a and in_block:
            in_block = False
            blocks.append((lo, i - 1))
    if in_block:
        blocks.append((lo, len(is_active) - 1))

    # Panel
    mb = mids / 1e3   # kb
    # Shade true active blocks first (background)
    ymax = max(np.max(w2_est) * 1.1, np.max(w2_truth) * 1.1, 1e-30)
    for (a, b) in blocks:
        ax.axvspan(mb[a], mb[b], color="#fee08b", alpha=0.6, lw=0, zorder=1)
    if blocks:
        ax.axvspan(np.nan, np.nan, color="#fee08b", alpha=0.6, lw=0,
                   label="true active blocks")

    # Estimated ŵ² as a stem-like plot
    ax.plot(mb, w2_est, color="#4393c3", lw=0.6, zorder=2,
            label=r"$\hat{w}^2_s$ (SARGENT estimate)")
    ax.scatter(mb, w2_est, s=4, color="#4393c3", alpha=0.6, zorder=2)

    ax.set_yscale("log")
    ax.set_xlabel("Genomic position (kb)")
    ax.set_ylabel(r"Per-tree variance weight $\hat{w}^2_s$  (log scale)")
    ax.set_title("A — SARGENT identifies trait-relevant genomic regions "
                 "(Block 1.2 cell: admixed, contrast 1000, sparse)",
                 loc="left", fontsize=11)
    ax.grid(True, alpha=0.25, which="both")
    ax.legend(loc="lower right", fontsize=9, frameon=True)


# ---------- Panel B: Variance-capture curves ----------

def load_capture(case_dir):
    case = os.path.join(case_dir, "case")
    est  = os.path.join(case_dir, "est")
    truth_path = case + ".full.true_ws2.tsv"
    if not os.path.exists(truth_path):
        truth_path = case + ".true_ws2.tsv"
    w2_truth = np.loadtxt(truth_path, skiprows=1)[:, -1]
    w2_est   = np.loadtxt(est + ".w2.tsv", skiprows=1)[:, -1]
    return w2_truth, w2_est


def capture_curve(w):
    w = np.sort(np.asarray(w, dtype=float))[::-1]
    cum = np.cumsum(w)
    total = cum[-1] if cum[-1] > 0 else 1.0
    k = np.arange(1, len(w) + 1) / len(w)   # normalised top-k fraction
    return k, cum / total


def fig4_panel_b(ax, base_phase10, sparse_cells=None, diffuse_cells=None):
    if sparse_cells is None:
        sparse_cells = sorted(glob.glob(os.path.join(base_phase10, "sparse_s*")))
    if diffuse_cells is None:
        diffuse_cells = sorted(glob.glob(os.path.join(base_phase10, "diffuse_s*")))

    cmap = plt.get_cmap("tab10")
    c_sparse = cmap(1)   # orange
    c_diffuse = cmap(0)  # blue

    # Plot truth + est curves for each cell in each arch (alpha-blended)
    def plot_arch(cells, color, label_truth, label_est):
        first = True
        for cd in cells:
            try:
                w_t, w_e = load_capture(cd)
            except Exception:
                continue
            for w, ls, alpha, label in [
                (w_t, "-",  0.8, label_truth if first else None),
                (w_e, "--", 0.6, label_est   if first else None),
            ]:
                k, f = capture_curve(w)
                ax.plot(k, f, color=color, lw=1.2 if ls == "-" else 1.0,
                        ls=ls, alpha=alpha, label=label)
            first = False

    plot_arch(sparse_cells, c_sparse,
              "sparse truth", "sparse estimated (SARGENT)")
    plot_arch(diffuse_cells, c_diffuse,
              "diffuse truth", "diffuse estimated (SARGENT)")

    ax.set_xscale("log")
    ax.set_xlim(1e-3, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Top-k trees, fraction of total (k/S)")
    ax.set_ylabel(r"Cumulative fraction of $\sum_s \hat{w}^2_s$")
    ax.set_title("B — Variance-capture curves: sparse vs diffuse "
                 "architectures (Phase 10 admixed)", loc="left", fontsize=11)
    ax.grid(True, alpha=0.25, which="both")
    ax.legend(loc="lower right", fontsize=8, frameon=True)


# ---------- Panel C: AUC across contrast × demography ----------

def load_auc_grid(base_b12):
    rows = []
    for jf in sorted(glob.glob(os.path.join(base_b12, "c*_d*_s*", "est.oracle.json"))):
        cell = os.path.basename(os.path.dirname(jf))
        m = re.match(r"c([\d.]+)_d(\w+)_s(\d+)", cell)
        if not m:
            continue
        contrast, demog, seed = float(m.group(1)), m.group(2), int(m.group(3))
        d = json.load(open(jf))
        auc = d.get("ws2_recovery", {}).get("auc_active")
        ceil = d.get("loc_ceiling_auc")
        if isinstance(auc, (int, float)):
            rows.append({"contrast": contrast, "demog": demog, "seed": seed,
                         "auc": auc, "ceiling": ceil})
    return rows


def fig4_panel_c(ax, base_b12):
    rows = load_auc_grid(base_b12)
    if not rows:
        ax.text(0.5, 0.5, "no Block 1.2 oracle JSONs found",
                ha="center", va="center", transform=ax.transAxes)
        return
    contrasts = sorted({r["contrast"] for r in rows}, reverse=True)   # 1000, 100, 10
    demogs    = sorted({r["demog"] for r in rows})                    # admixed, panmictic
    means = {}
    sems = {}
    means_ceil = {}
    for d in demogs:
        means[d]      = []
        sems[d]       = []
        means_ceil[d] = []
        for c in contrasts:
            vs = [r["auc"] for r in rows if r["demog"] == d and r["contrast"] == c]
            vs = np.array(vs)
            cs = np.array([r["ceiling"] for r in rows
                           if r["demog"] == d and r["contrast"] == c
                           and isinstance(r["ceiling"], (int, float))])
            if len(vs):
                means[d].append(vs.mean())
                sems[d].append(vs.std(ddof=1) / np.sqrt(len(vs)) if len(vs) > 1 else 0)
            else:
                means[d].append(float("nan"))
                sems[d].append(0)
            means_ceil[d].append(cs.mean() if len(cs) else float("nan"))

    x = np.arange(len(contrasts))
    width = 0.36
    color_admix    = "#762a83"
    color_panmict  = "#1a9850"
    # Bars: estimator AUC
    ax.bar(x - width/2, means["admixed"], width, yerr=sems["admixed"],
           color=color_admix, capsize=3,
           edgecolor="white", linewidth=0.5,
           label="admixed: SARGENT estimator AUC")
    ax.bar(x + width/2, means["panmictic"], width, yerr=sems["panmictic"],
           color=color_panmict, capsize=3,
           edgecolor="white", linewidth=0.5,
           label="panmictic: SARGENT estimator AUC")
    # Overlay ceiling as hollow markers
    ax.scatter(x - width/2, means_ceil["admixed"], s=60, marker="D",
               facecolor="none", edgecolor=color_admix, linewidth=1.4,
               label="admixed: idealized ceiling (V at uniform null)")
    ax.scatter(x + width/2, means_ceil["panmictic"], s=60, marker="D",
               facecolor="none", edgecolor=color_panmict, linewidth=1.4,
               label="panmictic: idealized ceiling")

    ax.axhline(0.5, color="0.6", lw=0.8, linestyle=":", label="chance (AUC = 0.5)")
    ax.set_xticks(x, [f"{int(c)}" for c in contrasts])
    ax.set_xlabel("Block:floor ws² contrast")
    ax.set_ylabel("Localization AUC")
    ax.set_ylim(0.40, 1.02)
    ax.set_title("C — Localization AUC across demography and contrast "
                 "(Block 1.2 v1; 5 seeds per cell)", loc="left", fontsize=11)
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="lower left", fontsize=8, frameon=True)


# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cell", default=os.path.expanduser(
        "~/test-sargent/block1_2_v1/c1000_dadmixed_s1"),
        help="Representative cell for Panel A.")
    ap.add_argument("--base-b12", default=os.path.expanduser(
        "~/test-sargent/block1_2_v1"))
    ap.add_argument("--base-phase10", default=os.path.expanduser(
        "~/test-sargent/phase10_admix"))
    ap.add_argument("--out", default=None,
                    help="Output prefix (default: ./fig4).")
    args = ap.parse_args()

    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.35, wspace=0.28)
    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1, 0])
    axC = fig.add_subplot(gs[1, 1])

    fig4_panel_a(axA, args.cell)
    fig4_panel_b(axB, args.base_phase10)
    fig4_panel_c(axC, args.base_b12)

    fig.suptitle(
        "Figure 4 — SARGENT localizes trait variance to genomic regions at "
        "recombination-aware resolution",
        fontsize=12, y=1.00)

    out_prefix = args.out or "fig4"
    for ext in ("pdf", "png"):
        path = f"{out_prefix}.{ext}"
        fig.savefig(path, dpi=160, bbox_inches="tight")
        print(f"Wrote {path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
