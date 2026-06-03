#!/usr/bin/env python3
"""
plot.py — Fig 2 generator (revised 2026-06-02, coverage-axis layout).

Two-panel within-population figure for the manuscript.  Same 50 cells,
same 4 arms in both panels:

  Panel A: ĥ² recovery vs typed-SNP coverage.  Truth dashed at 0.60.
           BayesR's h² (h2_snp from summary.json) shows what the typed-
           SNP fixed effect alone explains; GBLUP / SARGENT-U / SARGENT
           h² shown alongside.  Error bars are ±1 SE across seeds.

  Panel B: r²_calib vs typed-SNP coverage.  Same 4 arms, same cells.
           Error bars are ±1 SE across seeds.

Sparse architecture, panmictic, single population, h²_SNP = h²_ARG = 0.3,
n_train = 600, L = 300 kb, 10 seeds per coverage level.

Naming:
  - "ARG-BLUP-uniform" → **SARGENT-U** (a special case of SARGENT: no
    fixed effect, uniform per-tree weights).
  - SARGENT = full pipeline (BayesR fixed effect + per-tree-weighted
    ARG random effect).

Reads:
  $BASE_RVA/a<arch>_cv<cov>_d<demog>_s<seed>/est.predict.json   (R², h²)
  $BASE_RVA/a<arch>_cv<cov>_d<demog>_s<seed>/est.summary.json   (BayesR h²)
"""

import argparse
import glob
import json
import os
import re
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ARMS = ("sargent", "argunif", "gblup", "bayesR")
ARM_COL = {
    "sargent": "#762a83",  # purple
    "argunif": "#1a9850",  # green
    "gblup":   "#777777",  # grey
    "bayesR":  "#d6604d",  # warm
}
ARM_LBL = {
    "sargent": "SARGENT (full)",
    "argunif": "SARGENT-U",
    "gblup":   "GBLUP",
    "bayesR":  "BayesR",
}
ARM_MARK = {
    "sargent": "o",
    "argunif": "s",
    "gblup":   "^",
    "bayesR":  "D",
}
TRUE_H2 = 0.60


def _load(base, arch="sparse"):
    """coverage → list of {arm: r2_calib, arm_h2: float} dicts per seed."""
    rows = defaultdict(list)
    for jf in sorted(glob.glob(os.path.join(base, f"a{arch}_cv*_d*_s*",
                                            "est.predict.json"))):
        cell_dir = os.path.dirname(jf)
        cell = os.path.basename(cell_dir)
        m = re.match(rf"a{arch}_cv([\d.]+)_d\w+_s(\d+)", cell)
        if not m: continue
        cov = float(m.group(1))
        seed = int(m.group(2))
        d = json.load(open(jf))
        sj = os.path.join(cell_dir, "est.summary.json")
        s = json.load(open(sj)) if os.path.exists(sj) else {}
        rows[cov].append({
            "seed":       seed,
            "bayesR_r2":  d.get("bayesR",           {}).get("r2_calib"),
            "gblup_r2":   d.get("gblup",            {}).get("r2_calib"),
            "argunif_r2": d.get("arg_blup_uniform", {}).get("r2_calib"),
            "sargent_r2": d.get("sargent",          {}).get("r2_calib"),
            "bayesR_h2":  s.get("h2_snp"),
            "gblup_h2":   d.get("gblup_h2"),
            "argunif_h2": d.get("arg_blup_uniform_h2"),
            "sargent_h2": s.get("h2"),
        })
    return rows


def _mean_se(vs):
    vs = np.asarray([v for v in vs if isinstance(v, (int, float))],
                    dtype=float)
    if len(vs) == 0: return float("nan"), 0
    return float(vs.mean()), float(vs.std(ddof=1) / np.sqrt(len(vs)) if len(vs) > 1 else 0)


def _plot_panel(ax, rows, field_suffix, ylabel, title, ylim_max=None,
                truth_line=None):
    covs = sorted(rows.keys(), reverse=True)
    if not covs:
        ax.set_title(f"no data found for {title}", loc="left"); return

    # tiny x-offset per arm so error bars don't overlap
    offsets = np.linspace(-0.020, 0.020, len(ARMS))

    max_y_seen = 0.0
    for i, arm in enumerate(ARMS):
        means, ses = [], []
        for cv in covs:
            m, se = _mean_se([r[f"{arm}_{field_suffix}"] for r in rows[cv]])
            means.append(m); ses.append(se)
        means = np.array(means); ses = np.array(ses)
        if np.isfinite(means).any():
            max_y_seen = max(max_y_seen, float(np.nanmax(means + ses)))
        xs = np.array(covs) + offsets[i]
        ax.errorbar(xs, means, yerr=ses,
                    color=ARM_COL[arm], lw=1.6,
                    marker=ARM_MARK[arm], markersize=7,
                    markerfacecolor=ARM_COL[arm], markeredgecolor="white",
                    capsize=4, elinewidth=1.3,
                    linestyle="-" if arm in ("sargent", "argunif") else "--",
                    label=ARM_LBL[arm])

    if truth_line is not None:
        ax.axhline(truth_line, color="black", lw=1.0, linestyle="--",
                   label=f"True $h^2 = {truth_line}$")

    ax.set_xlim(1.07, 0.0)
    if ylim_max is None:
        ylim_max = max(0.65, max_y_seen + 0.05)
    ax.set_ylim(0, ylim_max)
    ax.set_xlabel(r"Typed-SNP coverage $f_\mathrm{typed}$")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=10.5)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower right", fontsize=9, frameon=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-rva",  default=os.path.expanduser("~/test-sargent/fig2_n600"))
    ap.add_argument("--arch",      default="sparse",
                    help="ws²-architecture to plot (default: sparse).")
    ap.add_argument("--out",       default=None, help="Output prefix (default: ./fig2).")
    args = ap.parse_args()

    rows = _load(args.base_rva, arch=args.arch)
    n_cells = sum(len(v) for v in rows.values())
    print(f"Loaded {n_cells} cells across {len(rows)} coverage levels "
          f"(arch={args.arch}) from {args.base_rva}.")

    fig, axes = plt.subplots(1, 2, figsize=(14, 4.8))
    _plot_panel(axes[0], rows, "h2",
                ylabel=r"Estimated $\hat{h}^2$",
                title="A — $\\hat{h}^2$ recovery (sparse, balanced mix)",
                ylim_max=1.0, truth_line=TRUE_H2)
    _plot_panel(axes[1], rows, "r2",
                ylabel=r"Held-out $r^2_\mathrm{calib}$",
                title="B — Prediction accuracy",
                ylim_max=None)
    fig.suptitle(
        "Figure 2 — Within-population: SARGENT and SARGENT-U recover $h^2$ "
        "and held-out prediction $r^2$ where BayesR / GBLUP fall off "
        "(sparse architecture, $n=600$, $L=300$ kb, $h^2_\\mathrm{total}=0.60$)",
        fontsize=11.0, y=1.02)
    fig.tight_layout()

    out_prefix = args.out or "fig2"
    for ext in ("pdf", "png"):
        path = f"{out_prefix}.{ext}"
        fig.savefig(path, dpi=160, bbox_inches="tight")
        print(f"Wrote {path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
