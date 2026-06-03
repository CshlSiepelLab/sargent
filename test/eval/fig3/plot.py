#!/usr/bin/env python3
"""
plot.py — Fig 3 generator (cross-ancestry portability, revised 2026-06-02).

Two-panel cross-ancestry figure for the manuscript:

  Panel A: Target-injection axis at F_ST=0.16, typed=0.70.
           Grouped bars over mix_B ∈ {0.0, 0.05, 0.20} for 4 arms.  Light
           error bars (±1 SE across seeds).  Shows that including target-
           population individuals in the training set is what unlocks
           cross-ancestry prediction: SARGENT-U and GBLUP lift dramatically;
           BayesR (and SARGENT, which inherits the BayesR fixed effect)
           lift much less.

  Panel B: F_ST ladder at typed=0.70, mix_B=0.20 (the predictive regime).
           Lines for 4 arms with light ±1 SE error bars.  All solid; arms
           are distinguished by marker shape, not line style.  Shows
           SARGENT-U dominates across the entire ladder; SARGENT collapses
           with BayesR at F_ST ≥ 0.16 due to the BayesR fixed-effect
           cross-ancestry liability.

Data: `test/eval/block1_1_v2_results.csv` (server-side Block 1.1 v2,
2026-06-02; 4 F_ST × 3 mix × 2 typed × 5 seeds = 120 cells).
"""

import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ARMS = ("sargent", "argunif", "gblup", "bayesR")
ARM_COL = {
    "sargent": "#762a83",  # purple
    "argunif": "#1a9850",  # green
    "gblup":   "#3690c0",  # blue
    "bayesR":  "#d6604d",  # warm red
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
ARM_CSV_COL = {
    "sargent": "sargent_r2c",
    "argunif": "argunif_r2c",
    "gblup":   "gblup_r2c",
    "bayesR":  "bayesR_r2c",
}


def _load_rows(csv_path):
    out = []
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            try:
                out.append({
                    "fst":   float(r["fst"]),
                    "mix":   float(r["mix"]),
                    "typed": float(r["typed"]),
                    "seed":  int(r["seed"]),
                    **{a: float(r[ARM_CSV_COL[a]]) for a in ARMS},
                })
            except (ValueError, KeyError):
                continue
    return out


def _filter(rows, **kw):
    out = []
    for r in rows:
        ok = True
        for k, v in kw.items():
            if isinstance(v, float):
                if abs(r[k] - v) > 1e-9: ok = False; break
            else:
                if r[k] != v: ok = False; break
        if ok: out.append(r)
    return out


def _msm(vs):
    vs = np.asarray(vs, dtype=float)
    if len(vs) == 0: return float("nan"), 0
    m = float(vs.mean())
    se = float(vs.std(ddof=1)/np.sqrt(len(vs)) if len(vs) > 1 else 0)
    return m, se


# ---------- Panel A: target-injection axis at F_ST=0.16, typed=0.70 ----------

def fig3_panel_a(ax, rows, fst=0.16, typed=0.70):
    mix_levels = [0.0, 0.05, 0.20]
    x = np.arange(len(mix_levels))
    width = 0.20
    arm_offsets = (np.arange(len(ARMS)) - (len(ARMS) - 1) / 2) * width

    max_y = 0.0
    for i, arm in enumerate(ARMS):
        means, ses = [], []
        for mx in mix_levels:
            subset = _filter(rows, fst=fst, mix=mx, typed=typed)
            m, se = _msm([r[arm] for r in subset])
            means.append(m); ses.append(se)
        means = np.array(means); ses = np.array(ses)
        if np.isfinite(means).any():
            max_y = max(max_y, float(np.nanmax(means + ses)))
        ax.bar(x + arm_offsets[i], means, width,
               yerr=ses,
               color=ARM_COL[arm], edgecolor="white", linewidth=0.6,
               capsize=2.5, ecolor="0.35", error_kw={"alpha": 0.7, "lw": 0.9},
               label=ARM_LBL[arm])

    ax.set_xticks(x, [r"$\mathrm{mix}_B=0$",
                      r"$\mathrm{mix}_B=0.05$",
                      r"$\mathrm{mix}_B=0.20$"])
    ax.set_ylabel(r"Held-out $r^2_\mathrm{calib}$ (cross-ancestry)")
    ax.set_ylim(0, max(0.32, max_y + 0.04))
    ax.set_title(f"A — Target injection unlocks cross-ancestry prediction "
                 f"($F_{{ST}}={fst}$, typed$=${typed:.2f})",
                 loc="left", fontsize=10.5)
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="upper left", fontsize=9, frameon=True, ncol=2)


# ---------- Panel B: F_ST ladder at typed=0.70, mix=0.20 ----------

def fig3_panel_b(ax, rows, typed=0.70, mix=0.20):
    fst_levels = sorted({r["fst"] for r in rows})
    max_y = 0.0
    for arm in ARMS:
        means, ses = [], []
        for f in fst_levels:
            subset = _filter(rows, fst=f, typed=typed, mix=mix)
            m, se = _msm([r[arm] for r in subset])
            means.append(m); ses.append(se)
        means = np.array(means); ses = np.array(ses)
        max_y = max(max_y, float(np.nanmax(means + ses)))
        ax.errorbar(fst_levels, means, yerr=ses,
                    color=ARM_COL[arm], lw=1.8,
                    marker=ARM_MARK[arm], markersize=8,
                    markerfacecolor=ARM_COL[arm], markeredgecolor="white",
                    capsize=3, elinewidth=0.9, ecolor="0.45",
                    linestyle="-", label=ARM_LBL[arm])

    ax.set_xlabel(r"$F_{ST}$ (training pop A → test pop B)")
    ax.set_ylabel(r"Held-out $r^2_\mathrm{calib}$ (cross-ancestry)")
    ax.set_xlim(min(fst_levels) - 0.02, max(fst_levels) + 0.02)
    ax.set_ylim(0, max(0.32, max_y + 0.03))
    ax.set_title(f"B — Across the $F_{{ST}}$ ladder with target injection "
                 f"(typed$=${typed:.2f}, $\\mathrm{{mix}}_B={mix}$): "
                 f"SARGENT-U holds; SARGENT collapses with BayesR",
                 loc="left", fontsize=10.5)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right", fontsize=9, frameon=True)


# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=os.path.join(
        os.path.dirname(__file__), "..", "block1_1_v2_results.csv"))
    ap.add_argument("--out", default=None, help="Output prefix (default: ./fig3).")
    ap.add_argument("--panel-a-fst",   type=float, default=0.16)
    ap.add_argument("--panel-a-typed", type=float, default=0.70)
    ap.add_argument("--panel-b-typed", type=float, default=0.70)
    ap.add_argument("--panel-b-mix",   type=float, default=0.20)
    args = ap.parse_args()

    rows = _load_rows(args.csv)
    print(f"Loaded {len(rows)} cells from {args.csv}")

    fig, axes = plt.subplots(1, 2, figsize=(14, 4.8))
    fig3_panel_a(axes[0], rows, fst=args.panel_a_fst, typed=args.panel_a_typed)
    fig3_panel_b(axes[1], rows, typed=args.panel_b_typed, mix=args.panel_b_mix)
    fig.suptitle(
        "Figure 3 — Cross-ancestry portability: target injection is the "
        "lever; SARGENT-U dominates across $F_{ST}$ when target individuals "
        "are present in training",
        fontsize=11.0, y=1.02)
    fig.tight_layout()

    out_prefix = args.out or "fig3"
    for ext in ("pdf", "png"):
        path = f"{out_prefix}.{ext}"
        fig.savefig(path, dpi=160, bbox_inches="tight")
        print(f"Wrote {path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
