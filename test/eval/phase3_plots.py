#!/usr/bin/env python3
"""Generate Phase 3 diagnostic plots from /Users/acs/test-sargent/overnight/."""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path("/Users/acs/test-sargent/overnight")
FIG_DIR = ROOT / "figs"
FIG_DIR.mkdir(exist_ok=True)


def load(path):
    with open(path) as f:
        return json.load(f)


def gather_multiseed():
    """Returns {(seed, f): {method: r2}}."""
    out = {}
    for d in sorted((ROOT / "multiseed").iterdir()):
        seed = int(d.name.split("_")[0][1:])
        f = float(d.name.split("_f")[1])
        seq = load(d / "seq.predict.json")
        em = load(d / "em.predict.json")
        out[(seed, f)] = {
            "BayesR":         seq["bayesR"]["r2"],
            "GBLUP":          seq["gblup"]["r2"],
            "ARG-BLUP-uni":   seq["arg_blup_uniform"]["r2"],
            "Seq SARGENT":    seq["sargent"]["r2"],
            "Joint EM":       em["sargent"]["r2"],
        }
    return out


def gather_h2sweep():
    """Returns {(seed, case): {method: {h2_snp,h2_arg}, true: {h2_snp,h2_arg}}}."""
    out = {}
    for d in sorted((ROOT / "h2sweep").iterdir()):
        seed = int(d.name.split("_")[0][1:])
        case = "_".join(d.name.split("_")[1:])
        p = load(d / "case.params.json")
        s = load(d / "seq.summary.json")
        e = load(d / "em.summary.json")
        out[(seed, case)] = {
            "true":  (p["h2_snp_realised"], p["h2_arg_realised"]),
            "seq":   (s["h2_snp"], s["h2_arg"]),
            "joint": (e["h2_snp"], e["h2_arg"]),
        }
    return out


def plot_bed_coverage(data):
    fs = sorted({k[1] for k in data}, reverse=True)
    methods = ["BayesR", "GBLUP", "ARG-BLUP-uni", "Seq SARGENT", "Joint EM"]
    colors  = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

    fig, ax = plt.subplots(figsize=(7, 5))
    for method, color in zip(methods, colors):
        means, sds = [], []
        for f in fs:
            vals = [data[(s, f)][method] for s in sorted({k[0] for k in data})]
            means.append(np.mean(vals))
            sds.append(np.std(vals, ddof=1))
        means, sds = np.array(means), np.array(sds)
        ax.plot(fs, means, "o-", label=method, color=color, linewidth=2)
        ax.fill_between(fs, means - sds, means + sds, alpha=0.15, color=color)

    ax.set_xscale("log")
    ax.set_xlabel("Typed-SNP coverage f")
    ax.set_ylabel("Test-set prediction R²")
    ax.set_title("Prediction R² vs BED coverage\n"
                 "(balanced h²_snp = h²_arg = 0.4, n=200, 5 seeds, ±SD)")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower left")
    ax.set_xticks(fs)
    ax.set_xticklabels([f"{f:g}" for f in fs])
    fig.tight_layout()
    out = FIG_DIR / "phase3_bed_coverage.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


def plot_h2_decomposition(data):
    cases = ["snp_only", "snp_heavy", "balanced", "arg_heavy", "arg_only"]
    seeds = sorted({k[0] for k in data})

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    width = 0.27
    x = np.arange(len(cases))

    for ax, method, label in [(axes[0], "seq", "Sequential SARGENT"),
                              (axes[1], "joint", "Joint EM SARGENT")]:
        h2_snp_mean = []; h2_snp_sd = []
        h2_arg_mean = []; h2_arg_sd = []
        h2_tot_mean = []; h2_tot_sd = []
        h2_snp_true = []; h2_arg_true = []; h2_tot_true = []
        for c in cases:
            snp_vals = [data[(s, c)][method][0] for s in seeds]
            arg_vals = [data[(s, c)][method][1] for s in seeds]
            tot_vals = [a + b for a, b in zip(snp_vals, arg_vals)]
            h2_snp_mean.append(np.mean(snp_vals)); h2_snp_sd.append(np.std(snp_vals, ddof=1))
            h2_arg_mean.append(np.mean(arg_vals)); h2_arg_sd.append(np.std(arg_vals, ddof=1))
            h2_tot_mean.append(np.mean(tot_vals)); h2_tot_sd.append(np.std(tot_vals, ddof=1))
            true_snp = [data[(s, c)]["true"][0] for s in seeds]
            true_arg = [data[(s, c)]["true"][1] for s in seeds]
            h2_snp_true.append(np.mean(true_snp))
            h2_arg_true.append(np.mean(true_arg))
            h2_tot_true.append(np.mean(true_snp) + np.mean(true_arg))

        ax.bar(x - width, h2_snp_mean, width, yerr=h2_snp_sd,
               label="h²_snp inferred", color="#1f77b4", capsize=3)
        ax.bar(x,         h2_arg_mean, width, yerr=h2_arg_sd,
               label="h²_arg inferred", color="#d62728", capsize=3)
        ax.bar(x + width, h2_tot_mean, width, yerr=h2_tot_sd,
               label="h²_total inferred", color="#7f7f7f", capsize=3)
        # Truth as black caps over each bar
        for i in range(len(cases)):
            ax.hlines(h2_snp_true[i], x[i] - 1.5*width, x[i] - 0.5*width, colors="black", linewidth=2)
            ax.hlines(h2_arg_true[i], x[i] - 0.5*width, x[i] + 0.5*width, colors="black", linewidth=2)
            ax.hlines(h2_tot_true[i], x[i] + 0.5*width, x[i] + 1.5*width, colors="black", linewidth=2)

        ax.set_xticks(x)
        ax.set_xticklabels(cases, rotation=20)
        ax.set_ylabel("h²")
        ax.set_title(label)
        ax.set_ylim(-0.05, 1.1)
        ax.grid(alpha=0.3, axis="y")
        ax.legend(loc="upper center", fontsize=9)

    fig.suptitle("h² decomposition at f=0.5 (total h²=0.8, 3 seeds, ±SD).  "
                 "Black bars = realised truth.")
    fig.tight_layout()
    out = FIG_DIR / "phase3_h2_decomposition.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    ms = gather_multiseed()
    print(f"multiseed cells: {len(ms)}")
    plot_bed_coverage(ms)

    hs = gather_h2sweep()
    print(f"h2sweep cells: {len(hs)}")
    plot_h2_decomposition(hs)
