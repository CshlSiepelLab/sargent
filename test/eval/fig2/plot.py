#!/usr/bin/env python3
"""
plot.py — Fig 2 generator.

Two-panel figure for the manuscript:

  Panel A: mixture sweep at full coverage.  Stacked bars of estimated
           h²_SNP and h²_ARG vs the simulated mixture; reference lines
           for the simulated h²_ARG share and the total h².
  Panel B: coverage sweep at the balanced mixture.  Lines vs typed-SNP
           coverage for SARGENT (h²_SNP / h²_ARG / total), BayesR
           total, and ARG-BLUP-uniform total.  Shaded ribbons = ±1 SE
           across seeds.

All sims panmictic, single population, n_train=300, L=300 kb, total h²=0.6.

Reads:
  Panel A:  $BASE_FIG2/<mixname>_s<seed>/est.summary.json
  Panel B:  $BASE_RVA/a<arch>_cv<cov>_d<demog>_s<seed>/
             est.summary.json + est.predict.json
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


# ---------- Panel A: mixture sweep ----------

MIX_LABELS = ["allSNP", "heavySNP", "balanced", "heavyARG", "allARG"]
MIX_H2_ARG = {"allSNP": 0.0, "heavySNP": 0.15, "balanced": 0.30,
              "heavyARG": 0.45, "allARG": 0.60}
MIX_H2_SNP = {"allSNP": 0.60, "heavySNP": 0.45, "balanced": 0.30,
              "heavyARG": 0.15, "allARG": 0.0}
MIX_DISPLAY = {"allSNP": "all SNP\n(0.60/0)",
               "heavySNP": "heavy SNP\n(0.45/0.15)",
               "balanced": "balanced\n(0.30/0.30)",
               "heavyARG": "heavy ARG\n(0.15/0.45)",
               "allARG":  "all ARG\n(0/0.60)"}


def load_panel_a(base):
    rows = defaultdict(list)
    for jf in sorted(glob.glob(os.path.join(base, "*_s*", "est.summary.json"))):
        cell = os.path.basename(os.path.dirname(jf))
        m = re.match(r"(\w+)_s(\d+)", cell)
        if not m:
            continue
        mix, seed = m.group(1), int(m.group(2))
        if mix not in MIX_LABELS:
            continue
        d = json.load(open(jf))
        rows[mix].append({
            "seed": seed,
            "h2_snp": d.get("h2_snp"),
            "h2_arg": d.get("h2_arg"),
            "h2":     d.get("h2"),
        })
    return rows


def fig2_panel_a(ax, rows):
    x = np.arange(len(MIX_LABELS))
    width = 0.62

    means_snp, errs_snp, means_arg, errs_arg, means_total, errs_total = ([] for _ in range(6))
    for mix in MIX_LABELS:
        rs = rows.get(mix, [])
        snps = np.array([r["h2_snp"] for r in rs if isinstance(r["h2_snp"], (int, float))])
        args_ = np.array([r["h2_arg"] for r in rs if isinstance(r["h2_arg"], (int, float))])
        tots = np.array([r["h2"]      for r in rs if isinstance(r["h2"], (int, float))])
        def msm(v):
            return (v.mean(), v.std(ddof=1)/np.sqrt(len(v))) if len(v) > 1 else (v.mean() if len(v) else float("nan"), 0)
        m_s, e_s = msm(snps)
        m_a, e_a = msm(args_)
        m_t, e_t = msm(tots)
        means_snp.append(m_s); errs_snp.append(e_s)
        means_arg.append(m_a); errs_arg.append(e_a)
        means_total.append(m_t); errs_total.append(e_t)

    means_snp = np.array(means_snp); errs_snp = np.array(errs_snp)
    means_arg = np.array(means_arg); errs_arg = np.array(errs_arg)
    means_total = np.array(means_total); errs_total = np.array(errs_total)

    c_snp = "#d6604d"   # warm
    c_arg = "#4393c3"   # cool

    ax.bar(x, means_snp, width, label=r"$\hat{h}^2_\mathrm{SNP}$ (typed)",
           color=c_snp, edgecolor="white", linewidth=0.5)
    ax.bar(x, means_arg, width, bottom=means_snp,
           label=r"$\hat{h}^2_\mathrm{ARG}$ (excess via genealogy)",
           color=c_arg, edgecolor="white", linewidth=0.5)
    # Error bars on the bar tops (total h² SE)
    ax.errorbar(x, means_snp + means_arg, yerr=errs_total, fmt="none",
                ecolor="black", capsize=3, lw=1)

    # True ARG share marker on each bar (line at h²_SNP_true, above which
    # truth is ARG).  Annotates the "right answer" the bar should display.
    truth_x_left  = x - width/2
    truth_x_right = x + width/2
    truth_snp = [MIX_H2_SNP[m] for m in MIX_LABELS]
    truth_total = [MIX_H2_SNP[m] + MIX_H2_ARG[m] for m in MIX_LABELS]
    for i in range(len(x)):
        ax.plot([truth_x_left[i], truth_x_right[i]],
                [truth_snp[i], truth_snp[i]],
                color="black", lw=1.4, linestyle="--", zorder=5)
    ax.plot([], [], color="black", lw=1.4, linestyle="--",
            label=r"True $h^2_\mathrm{SNP}$ / $h^2_\mathrm{ARG}$ boundary")
    # Total truth (0.6) reference
    ax.axhline(0.60, color="0.4", lw=0.8, linestyle=":", label=r"True $h^2_\mathrm{total}=0.6$")

    ax.set_xticks(x, [MIX_DISPLAY[m] for m in MIX_LABELS])
    ax.set_xlabel(r"True mixture ($h^2_\mathrm{SNP}$ / $h^2_\mathrm{ARG}$)")
    ax.set_ylabel(r"Estimated $h^2$ components")
    ax.set_ylim(0, 1.0)
    ax.set_title("A — Variance decomposition under five true mixtures (full typed-SNP coverage)",
                 loc="left", fontsize=11)
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="upper right", fontsize=9, frameon=True)


# ---------- Panel B: coverage sweep ----------

def load_panel_b(base):
    rows = defaultdict(list)
    for jf in sorted(glob.glob(os.path.join(base, "a*_cv*_d*_s*", "est.summary.json"))):
        cell = os.path.basename(os.path.dirname(jf))
        m = re.match(r"a(\w+)_cv([\d.]+)_d(\w+)_s(\d+)", cell)
        if not m:
            continue
        arch, cov, demog, seed = m.group(1), float(m.group(2)), m.group(3), int(m.group(4))
        d = json.load(open(jf))
        # The predict.json has the per-arm fitted h² (gblup_h2, arg_blup_uniform_h2)
        # but they come from the eval-side fit_blup, not from sargent itself.
        pj = os.path.join(os.path.dirname(jf), "est.predict.json")
        p = json.load(open(pj)) if os.path.exists(pj) else {}
        rows[(arch, cov)].append({
            "seed": seed,
            "h2_snp":  d.get("h2_snp"),
            "h2_arg":  d.get("h2_arg"),
            "h2":      d.get("h2"),
            "gblup_h2":  p.get("gblup_h2"),
            "argunif_h2": p.get("arg_blup_uniform_h2"),
        })
    return rows


def fig2_panel_b(ax, rows, arch="scattered"):
    covs = sorted({k[1] for k in rows if k[0] == arch}, reverse=True)
    if not covs:
        ax.set_title(f"B — no Panel B data found for arch={arch}", loc="left")
        return

    def line(field, color, label, ls="-"):
        m_arr, se_arr = [], []
        for cv in covs:
            vs = [r[field] for r in rows[(arch, cv)] if isinstance(r[field], (int, float))]
            if vs:
                vs = np.array(vs)
                m_arr.append(vs.mean())
                se_arr.append(vs.std(ddof=1)/np.sqrt(len(vs)) if len(vs) > 1 else 0)
            else:
                m_arr.append(float("nan"))
                se_arr.append(0)
        m_arr = np.array(m_arr); se_arr = np.array(se_arr)
        ax.plot(covs, m_arr, color=color, lw=1.8, label=label, linestyle=ls,
                marker="o", markersize=5)
        ax.fill_between(covs, m_arr - se_arr, m_arr + se_arr, color=color, alpha=0.15)

    # SARGENT decomposition
    line("h2",       "#762a83", r"SARGENT $\hat{h}^2_\mathrm{total}$", ls="-")
    line("h2_snp",   "#d6604d", r"SARGENT $\hat{h}^2_\mathrm{SNP}$",   ls="-")
    line("h2_arg",   "#4393c3", r"SARGENT $\hat{h}^2_\mathrm{ARG}$",   ls="-")
    # External baselines
    line("gblup_h2", "#777777", r"GBLUP $\hat{h}^2$ (typed-SNP GRM)",  ls="--")
    line("argunif_h2", "#1a9850", r"ARG-BLUP-uniform $\hat{h}^2$ (ARG only)", ls="--")
    ax.axhline(0.60, color="0.4", lw=0.8, linestyle=":", label=r"True $h^2_\mathrm{total}=0.6$")

    ax.set_xlabel(r"Typed-SNP coverage $f_\mathrm{typed}$")
    ax.set_ylabel(r"Estimated $h^2$")
    ax.set_xlim(1.05, 0.05)   # reverse: 1.0 on left, 0.1 on right
    ax.set_ylim(0, 1.0)
    ax.set_title(f"B — Variance attribution shifts as typed coverage drops "
                 f"(balanced mixture; arch = {arch})", loc="left", fontsize=11)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right", fontsize=8, frameon=True)


# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-fig2", default=os.path.expanduser("~/test-sargent/fig2"),
                    help="Directory with Panel A mixture-sweep cells.")
    ap.add_argument("--base-rva",  default=os.path.expanduser("~/test-sargent/rare_variant"),
                    help="Directory with Panel B rare-variant ablation cells.")
    ap.add_argument("--arch-b", default="scattered",
                    help="Architecture to plot in Panel B (sparse|scattered).")
    ap.add_argument("--out", default=None,
                    help="Output prefix (default: ./fig2).")
    args = ap.parse_args()

    rows_a = load_panel_a(args.base_fig2)
    n_a = sum(len(v) for v in rows_a.values())
    rows_b = load_panel_b(args.base_rva)
    n_b = sum(len(v) for v in rows_b.values())
    print(f"Loaded Panel A: {n_a} cells across {len(rows_a)} mixtures.")
    print(f"Loaded Panel B: {n_b} cells across {len(rows_b)} (arch, cov) groups.")

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    fig2_panel_a(axes[0], rows_a)
    fig2_panel_b(axes[1], rows_b, arch=args.arch_b)
    fig.suptitle(
        "Figure 2 — SARGENT recovers total $h^2$; the SNP–ARG split shifts "
        "with true mixture and typed-SNP coverage",
        fontsize=11.5, y=1.01)
    fig.tight_layout()

    out_prefix = args.out or "fig2"
    for ext in ("pdf", "png"):
        path = f"{out_prefix}.{ext}"
        fig.savefig(path, dpi=160, bbox_inches="tight")
        print(f"Wrote {path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
