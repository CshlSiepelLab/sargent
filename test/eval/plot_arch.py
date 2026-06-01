#!/usr/bin/env python3
"""
plot_arch.py — Block 2.5.1 visualisation.

Reads one or more `arch_metrics_batch.json` (or single-case .arch_metrics.json)
files and emits:

  1. variance_capture.{pdf,png}   — overlay of variance-capture curves
                                    (truth + estimated) per case-group.
  2. me_arg_vs_arch.{pdf,png}     — M_e,ARG as a function of architecture
                                    (truth vs estimated, bar chart).
  3. branches_hist.{pdf,png}      — per-branch effect-size distribution
                                    histogram per architecture (if branches
                                    were captured).

Usage:
    python plot_arch.py <BATCH_JSON> [--out-prefix PREFIX] [--format pdf,png]
        [--group-key ws2_architecture]    # how to colour/group cases
"""

import argparse
import json
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def fig_variance_capture(rows, group_key, out_paths):
    """Overlay variance-capture curves, one colour per group."""
    by_group = defaultdict(list)
    for r in rows:
        g = r.get(group_key, "?")
        by_group[g].append(r)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    cmap = plt.get_cmap("tab10")
    colours = {g: cmap(i) for i, g in enumerate(sorted(by_group))}

    for which, title, ax in [("truth", "Truth", axes[0]),
                              ("est",   "Estimated (sargent)", axes[1])]:
        for g, rs in by_group.items():
            for r in rs:
                cc = r.get(which, {}).get("capture_curve")
                if not cc:
                    continue
                k = np.array(cc["k"], dtype=float)
                f = np.array(cc["f"], dtype=float)
                # k normalised by total trees for cross-arch comparability
                S = float(r.get(f"S_{'truth' if which=='truth' else 'est'}", k.max()))
                ax.plot(k / S, f, color=colours[g], alpha=0.6, lw=1.0)
        # legend handles per group (one per group)
        handles = [plt.Line2D([], [], color=c, label=str(g)) for g, c in colours.items()]
        ax.legend(handles=handles, title=group_key, loc="lower right",
                  fontsize=9, frameon=False)
        ax.set_xlabel("Top-k trees, fraction of total (k/S)")
        ax.set_ylabel("Captured fraction of ARG-variance")
        ax.set_xscale("log")
        ax.set_xlim(1e-3, 1.0)
        ax.set_ylim(0, 1.02)
        ax.grid(True, alpha=0.3)
        ax.set_title(title)

    fig.suptitle("Variance-capture curve per architecture (Block 2.5.1)", y=1.02)
    fig.tight_layout()
    for p in out_paths:
        fig.savefig(p, dpi=160, bbox_inches="tight")
    plt.close(fig)


def fig_me_arg(rows, group_key, out_paths):
    """Bar chart of M_e,ARG (truth vs estimated) per case-group."""
    by_group = defaultdict(lambda: {"truth": [], "est": []})
    for r in rows:
        g = r.get(group_key, "?")
        if isinstance(r.get("truth", {}).get("me_arg"), (int, float)):
            by_group[g]["truth"].append(r["truth"]["me_arg"])
        if isinstance(r.get("est", {}).get("me_arg"), (int, float)):
            by_group[g]["est"].append(r["est"]["me_arg"])

    groups = sorted(by_group)
    x = np.arange(len(groups))
    width = 0.36

    fig, ax = plt.subplots(figsize=(max(5, len(groups) * 0.8), 4.0))
    t_means = [np.mean(by_group[g]["truth"]) if by_group[g]["truth"] else float("nan")
               for g in groups]
    t_errs  = [np.std(by_group[g]["truth"], ddof=1) if len(by_group[g]["truth"]) > 1 else 0
               for g in groups]
    e_means = [np.mean(by_group[g]["est"])   if by_group[g]["est"]   else float("nan")
               for g in groups]
    e_errs  = [np.std(by_group[g]["est"], ddof=1) if len(by_group[g]["est"]) > 1 else 0
               for g in groups]
    ax.bar(x - width/2, t_means, width, yerr=t_errs, capsize=3, label="truth", color="C0")
    ax.bar(x + width/2, e_means, width, yerr=e_errs, capsize=3, label="estimated", color="C1")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=0)
    ax.set_xlabel(group_key)
    ax.set_ylabel("M_e,ARG  (effective # of contributing trees)")
    ax.set_title("Polygenicity index per architecture")
    ax.legend(frameon=False)
    ax.grid(True, axis="y", alpha=0.3, which="both")
    fig.tight_layout()
    for p in out_paths:
        fig.savefig(p, dpi=160, bbox_inches="tight")
    plt.close(fig)


def fig_age_bin_curves(rows, out_paths):
    """Variance-capture curves stratified by per-tree-age bin.

    Expects --age-bins arch_metrics output. Plots one panel per group
    (e.g., per architecture); within each panel, one curve per age bin.

    Returns True if anything was plotted, False otherwise.
    """
    # collect groups
    grouped = defaultdict(list)
    for r in rows:
        if r.get("truth", {}).get("age_bin_curves") or r.get("est", {}).get("age_bin_curves"):
            grouped[r.get("ws2_architecture", "?")].append(r)
    if not grouped:
        return False

    n_panels = len(grouped)
    fig, axes = plt.subplots(1, n_panels, figsize=(5.0 * n_panels, 4.0),
                             sharey=True)
    if n_panels == 1:
        axes = [axes]
    cmap = plt.get_cmap("plasma")

    for ax, (arch, rs) in zip(axes, grouped.items()):
        # Collect all unique age-bin labels (one per case, may differ)
        # Plot each case's curves with a colour per bin.
        for r in rs:
            ab = r.get("truth", {}).get("age_bin_curves") or r.get("est", {}).get("age_bin_curves")
            if not ab:
                continue
            for i, (lbl, c) in enumerate(ab.items()):
                k = np.array(c["k"], dtype=float)
                f = np.array(c["f"], dtype=float)
                S = c["S_bin"]
                col = cmap(i / max(1, len(ab) - 1))
                ax.plot(k / S, f, color=col, alpha=0.6, lw=1.0,
                        label=lbl if rs.index(r) == 0 else None)
        ax.set_xscale("log")
        ax.set_xlim(1e-3, 1.0)
        ax.set_ylim(0, 1.02)
        ax.grid(True, alpha=0.3)
        ax.set_xlabel("Top-k trees, fraction of bin total (k/S_bin)")
        if ax is axes[0]:
            ax.set_ylabel("Captured fraction of in-bin variance")
        ax.set_title(f"arch = {arch}")
        ax.legend(title="age bin (gen)", loc="lower right", fontsize=8, frameon=False)

    fig.suptitle("Variance-capture curve by tree-age tercile (Block 2.5.3)", y=1.02)
    fig.tight_layout()
    for p in out_paths:
        fig.savefig(p, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return True


def fig_branches_hist(rows, group_key, out_paths):
    """Histogram of per-branch variance contributions per group."""
    by_group = defaultdict(list)
    for r in rows:
        g = r.get(group_key, "?")
        # Prefer truth-branches when present, else estimated.
        for src in ("truth", "est"):
            if r.get(src, {}).get("branches", {}).get("sample"):
                by_group[g].extend(r[src]["branches"]["sample"])
                break

    if not any(by_group.values()):
        return False  # nothing to plot

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    cmap = plt.get_cmap("tab10")
    groups = sorted(by_group)
    colours = {g: cmap(i) for i, g in enumerate(groups)}
    for g in groups:
        vals = np.array(by_group[g], dtype=float)
        if len(vals) == 0:
            continue
        vals = vals[vals > 0]
        if len(vals) == 0:
            continue
        ax.hist(np.log10(vals), bins=60, density=True, histtype="step",
                color=colours[g], label=f"{g} (n={len(vals)})")
    ax.set_xlabel(r"log$_{10}$ (per-branch variance contribution = $\hat{w}^2_s \cdot bl(b)$)")
    ax.set_ylabel("density")
    ax.set_title("Per-branch effect-size distribution per architecture")
    ax.legend(frameon=False, title=group_key, fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    for p in out_paths:
        fig.savefig(p, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="arch_metrics_batch.json (or single-case JSON)")
    ap.add_argument("--out-prefix", default=None,
                    help="Output filename prefix (default: <input_dir>/arch_)")
    ap.add_argument("--format", default="pdf,png",
                    help="Comma-separated formats (default: pdf,png)")
    ap.add_argument("--group-key", default="ws2_architecture",
                    help="Field to group/colour cases by (default: ws2_architecture)")
    args = ap.parse_args()

    with open(args.input) as f:
        data = json.load(f)
    rows = data if isinstance(data, list) else [data]
    print(f"Loaded {len(rows)} cases from {args.input}")

    pref = args.out_prefix or os.path.join(os.path.dirname(args.input) or ".", "arch_")
    fmts = [f.strip() for f in args.format.split(",") if f.strip()]
    out_vcap = [f"{pref}variance_capture.{f}" for f in fmts]
    out_me   = [f"{pref}me_arg.{f}"            for f in fmts]
    out_br   = [f"{pref}branches_hist.{f}"     for f in fmts]

    fig_variance_capture(rows, args.group_key, out_vcap)
    fig_me_arg(rows, args.group_key, out_me)
    had_branches = fig_branches_hist(rows, args.group_key, out_br)
    out_age = [f"{pref}age_bins.{f}" for f in fmts]
    had_age = fig_age_bin_curves(rows, out_age)

    for p in out_vcap + out_me:
        print(f"Wrote {p}")
    if had_branches:
        for p in out_br:
            print(f"Wrote {p}")
    else:
        print("(no branches in input — skipped branches_hist; rerun arch_metrics without --no-branches)")
    if had_age:
        for p in out_age:
            print(f"Wrote {p}")
    else:
        print("(no age-bin curves in input — skipped; rerun arch_metrics with --age-bins)")


if __name__ == "__main__":
    main()
