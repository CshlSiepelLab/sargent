#!/usr/bin/env python3
"""
plot_portability.py — Block 1.1 visualisation.

Reads `<BASE>/mix*_typed*_s*/est.predict.json` produced by `portability.sh`
and emits two figures:

  1. <BASE>/portability_curve.{pdf,png}
     The headline: SARGENT − ARG-BLUP-uniform gap as a function of f_B
     (training mix fraction from population B), one line per typed fraction.
     Errorbars are the within-cell SEM across seeds.

  2. <BASE>/portability_arms.{pdf,png}
     Per-arm comparison bars: BayesR / GBLUP / ARG-BLUP-uni / SARGENT R²
     on the B-held-out test set, faceted by (mix_B × typed_fraction).
     Errorbars are within-cell SEM across seeds.

Usage:
    python plot_portability.py <BASE>
        [--out-prefix PREFIX]    # default: <BASE>/portability_
        [--format pdf,png]       # default: both
"""

import argparse
import glob
import json
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CELL_RE = re.compile(r"mix([\d.]+)_typed([\d.]+)_s(\d+)")


def load_cells(base):
    rows = []
    for jf in sorted(glob.glob(os.path.join(base, "mix*_typed*_s*", "est.predict.json"))):
        cell = os.path.basename(os.path.dirname(jf))
        m = CELL_RE.match(cell)
        if not m:
            continue
        mix, typed, seed = float(m.group(1)), float(m.group(2)), int(m.group(3))
        d = json.load(open(jf))
        rows.append({
            "mix": mix, "typed": typed, "seed": seed,
            "bayesR":  (d.get("bayesR")           or {}).get("r2"),
            "gblup":   (d.get("gblup")            or {}).get("r2"),
            "argunif": (d.get("arg_blup_uniform") or {}).get("r2"),
            "sargent": (d.get("sargent")          or {}).get("r2"),
        })
    return rows


def cell_stats(rows, mix, typed, field):
    vals = [r[field] for r in rows
            if r["mix"] == mix and r["typed"] == typed
            and isinstance(r[field], (int, float))]
    if not vals:
        return float("nan"), float("nan"), 0
    if len(vals) == 1:
        return vals[0], 0.0, 1
    return float(np.mean(vals)), float(np.std(vals, ddof=1) / np.sqrt(len(vals))), len(vals)


def paired_gap_stats(rows, mix, typed, a, b):
    """SEM of paired difference R²(a) − R²(b) within (mix, typed) cell."""
    pairs = [(r[a], r[b]) for r in rows
             if r["mix"] == mix and r["typed"] == typed
             and isinstance(r[a], (int, float)) and isinstance(r[b], (int, float))]
    if not pairs:
        return float("nan"), float("nan"), 0
    diffs = np.array([x - y for (x, y) in pairs])
    if len(diffs) == 1:
        return float(diffs[0]), 0.0, 1
    return float(diffs.mean()), float(diffs.std(ddof=1) / np.sqrt(len(diffs))), len(diffs)


def fig_portability_curve(rows, out_paths):
    mixes  = sorted({r["mix"]   for r in rows})
    typeds = sorted({r["typed"] for r in rows})
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    for typed in typeds:
        xs, ys, errs = [], [], []
        for mix in mixes:
            mean, sem, n = paired_gap_stats(rows, mix, typed, "sargent", "argunif")
            xs.append(mix); ys.append(mean); errs.append(sem)
        ax.errorbar(xs, ys, yerr=errs, marker="o", capsize=3,
                    label=f"f_typed = {typed:.2f}")
    ax.axhline(0, color="0.6", lw=0.8, ls=":")
    ax.set_xlabel("Training B-fraction (f_B)")
    ax.set_ylabel("R²  SARGENT − ARG-BLUP-uniform  (B-held-out)")
    ax.set_title("Block 1.1 — cross-ancestry portability gain over ARG-BLUP-uniform")
    ax.legend(loc="best", fontsize=9, frameon=False)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    for p in out_paths:
        fig.savefig(p, dpi=160)
    plt.close(fig)


def fig_arms(rows, out_paths):
    mixes  = sorted({r["mix"]   for r in rows})
    typeds = sorted({r["typed"] for r in rows})
    arms = [("bayesR",  "BayesR"),
            ("gblup",   "GBLUP"),
            ("argunif", "ARG-BLUP-uni"),
            ("sargent", "SARGENT")]
    arm_keys = [a[0] for a in arms]
    n_arms = len(arms)
    width = 0.18

    nrows = len(typeds)
    fig, axes = plt.subplots(nrows, 1, figsize=(7.0, 2.6 * nrows),
                             sharex=True, sharey=True)
    if nrows == 1:
        axes = [axes]
    for ax, typed in zip(axes, typeds):
        x = np.arange(len(mixes))
        for j, (k, label) in enumerate(arms):
            ys, errs = [], []
            for mix in mixes:
                mean, sem, _ = cell_stats(rows, mix, typed, k)
                ys.append(mean); errs.append(sem)
            ax.bar(x + (j - n_arms / 2 + 0.5) * width, ys, width,
                   yerr=errs, capsize=2, label=label)
        ax.axhline(0, color="0.6", lw=0.8, ls=":")
        ax.set_xticks(x, [f"{m:.2f}" for m in mixes])
        ax.set_ylabel("R²  (B-held-out)")
        ax.set_title(f"f_typed = {typed:.2f}", fontsize=10, loc="left")
        ax.grid(True, axis="y", alpha=0.25)
    axes[-1].set_xlabel("Training B-fraction (f_B)")
    axes[0].legend(loc="upper left", fontsize=8, frameon=False, ncol=4)
    fig.suptitle("Block 1.1 — per-arm prediction R² on B test", y=1.00)
    fig.tight_layout()
    for p in out_paths:
        fig.savefig(p, dpi=160)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("base", help="Directory containing mix*_typed*_s*/ cells")
    ap.add_argument("--out-prefix", default=None,
                    help="Output filename prefix (default: <base>/portability_)")
    ap.add_argument("--format", default="pdf,png",
                    help="Comma-separated formats (default: pdf,png)")
    args = ap.parse_args()

    rows = load_cells(args.base)
    if not rows:
        sys.exit(f"No cells found under {args.base}/mix*_typed*_s*/est.predict.json")
    print(f"Loaded {len(rows)} cells from {args.base}")

    pref = args.out_prefix or os.path.join(args.base, "portability_")
    fmts = [f.strip() for f in args.format.split(",") if f.strip()]

    p_curve = [f"{pref}curve.{f}" for f in fmts]
    p_arms  = [f"{pref}arms.{f}" for f in fmts]
    fig_portability_curve(rows, p_curve)
    fig_arms(rows, p_arms)
    for p in p_curve + p_arms:
        print(f"Wrote {p}")


if __name__ == "__main__":
    main()
