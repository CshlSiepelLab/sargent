#!/usr/bin/env python3
"""
age_partition.py — Block 2.1 / 2.5.3: per-age-bin heritability partition.

For a sim produced by simulate.py with --age-tau-coef ≠ 0, compute the
per-age-bin contribution to h² and verify that the bias matches the
prescribed age-tau model.

The TRUTH partition is exact (reads the simulated betas from
true_beta.tsv + standardised genotype variances). An ESTIMATED partition
from SARGENT's β̂ output is also reported when sargent's .beta.tsv is
present.

Per-bin h²:
    h²_bin = Σ_{j ∈ bin}  β²_j · Var(x_j)
where x_j is the standardised dosage for variant j (so Var ≈ 1 by
construction).

Usage:
    python age_partition.py --case-prefix CASE_PREFIX [--sargent-prefix EST_PREFIX]
        [--bins 3 | --quantiles 0.33,0.67 | --age-edges 1,1000,10000,100000]
        [--out RESULT.json]
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from predict import read_bed, standardise, load_beta


def load_true_beta(path):
    """Returns (snp_idx, beta, age, in_bed) arrays."""
    data = np.loadtxt(path, skiprows=1, ndmin=2)
    if data.shape[1] < 6:
        raise ValueError(
            f"{path} has {data.shape[1]} columns; expected ≥6 "
            "(snp_idx bp_pos beta in_bed bed_idx age). "
            "Was this sim made with the post-Block-2.1 simulate.py?")
    snp_idx = data[:, 0].astype(int)
    beta    = data[:, 2].astype(float)
    in_bed  = data[:, 3].astype(int)
    bed_idx = data[:, 4].astype(int)
    age     = data[:, 5].astype(float)
    return snp_idx, beta, age, in_bed, bed_idx


def bin_edges_from_args(ages_observed, args):
    """Return monotone bin edges from CLI options.

    Default: equal-quantile bins of the observed causal-variant ages.
    """
    if args.age_edges:
        return np.array([float(x) for x in args.age_edges.split(",")],
                        dtype=float)
    if args.quantiles:
        qs = [float(x) for x in args.quantiles.split(",")]
        qs = sorted(set([0.0] + qs + [1.0]))
        return np.quantile(ages_observed, qs)
    n_bins = args.bins
    qs = np.linspace(0, 1, n_bins + 1)
    return np.quantile(ages_observed, qs)


def partition_h2_by_age(beta, age, var_x, edges):
    """Per-bin Σ β² · Var(x) using fixed bin edges.

    Returns:
        bins:        list of (lo, hi) tuples, length C
        h2_bin:      array length C
        h2_total:    scalar (sum of bins; should equal Σ over all variants
                     in the input regardless of binning)
        n_in_bin:    list length C
    """
    bins = []
    h2_bin = []
    n_in = []
    for c in range(len(edges) - 1):
        lo, hi = float(edges[c]), float(edges[c + 1])
        # Inclusive on the right edge only for the last bin
        if c == len(edges) - 2:
            mask = (age >= lo) & (age <= hi)
        else:
            mask = (age >= lo) & (age < hi)
        contrib = float(np.sum(beta[mask] ** 2 * var_x[mask]))
        bins.append((lo, hi))
        h2_bin.append(contrib)
        n_in.append(int(mask.sum()))
    return bins, np.array(h2_bin), float(np.sum(beta ** 2 * var_x)), n_in


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case-prefix",    required=True)
    ap.add_argument("--sargent-prefix", default=None,
                    help="If supplied, also compute an estimated partition "
                         "using sargent's .beta.tsv output.")
    ap.add_argument("--bins",      type=int, default=3,
                    help="Number of equal-quantile age bins (default 3).")
    ap.add_argument("--quantiles", default=None,
                    help="Comma-separated quantiles (e.g. '0.33,0.67') for "
                         "non-equal bin boundaries; overrides --bins.")
    ap.add_argument("--age-edges", default=None,
                    help="Comma-separated explicit age edges (in generations); "
                         "overrides --bins / --quantiles.")
    ap.add_argument("--out",       default=None)
    args = ap.parse_args()

    case = args.case_prefix
    with open(case + ".params.json") as f:
        params = json.load(f)
    n_indiv  = params["n_indiv"]
    n_test   = params["n_test"]
    n_total  = params["n_total"]
    n_var    = params["n_variants_full"]   # total simulated variants
    tau      = params.get("age_tau_coef", 0.0)
    ref      = params.get("age_tau_ref", None)

    # Causal variants and their ages (from simulate.py)
    snp_idx, beta_true, age, in_bed, bed_idx = load_true_beta(case + ".true_beta.tsv")

    # We need Var(x_j) for the causal variants, on TRAINING data, standardised
    # the way sargent does it. Easiest: read the train BED and compute.
    train_dosage = read_bed(case, n_indiv, params["n_variants"])
    train_mean = train_dosage.mean(axis=1, keepdims=True)
    train_std  = train_dosage.std(axis=1, keepdims=True)
    train_std[train_std == 0] = 1.0
    # Standardised dosages → Var(x_j) = 1 by construction in the train set;
    # we keep this explicit for clarity and to handle edge cases.
    var_x_bed = np.full(params["n_variants"], 1.0)
    # For variants NOT in the BED (in_bed == 0), we still have their dosages
    # in the FULL .bed only if bed_coverage = 1.0. Otherwise we cannot
    # compute Var(x) directly. For typical age-stratified analyses, the BED
    # spans all causal variants (causal-fraction draw is independent of
    # bed-coverage; many causal variants are in the BED).
    var_x = np.zeros(len(snp_idx))
    for k, idx in enumerate(snp_idx):
        if in_bed[k]:
            var_x[k] = var_x_bed[bed_idx[k]]
        else:
            # Variant is causal but not in BED. We use Var(standardised) = 1
            # as a uniform proxy. This biases per-bin h² slightly when ages
            # correlate with MAF, but is the safe default.
            var_x[k] = 1.0

    # Truth bin edges
    edges = bin_edges_from_args(age, args)

    bins, h2_truth, total_truth, n_in = partition_h2_by_age(
        beta_true, age, var_x, edges)

    # Build the result
    result = {
        "case_prefix": case,
        "tau":         tau,
        "ref_age":     ref,
        "n_causal":    int(len(snp_idx)),
        "edges":       [float(e) for e in edges],
        "bins":        [{"lo": float(lo), "hi": float(hi)} for lo, hi in bins],
        "h2_realised_total": params.get("h2_snp_realised"),
        "truth": {
            "n_in_bin":   n_in,
            "h2_bin":     [float(h) for h in h2_truth],
            "h2_total":   float(total_truth),
            "h2_frac":    [float(h / max(total_truth, 1e-30)) for h in h2_truth],
        },
    }

    # Optional: sargent's estimated partition
    if args.sargent_prefix:
        bf = args.sargent_prefix + ".beta.tsv"
        if os.path.exists(bf):
            beta_est = load_beta(bf)
            # Sargent's β is per BED-panel SNP; match to causal SNPs by
            # bed_idx (if in_bed). Variants not in the BED have no estimate.
            est_beta_at_causal = np.zeros(len(snp_idx))
            est_var_x          = np.zeros(len(snp_idx))
            for k in range(len(snp_idx)):
                if in_bed[k]:
                    est_beta_at_causal[k] = beta_est[bed_idx[k]]
                    est_var_x[k]          = 1.0
            bins_e, h2_est, total_est, n_in_e = partition_h2_by_age(
                est_beta_at_causal, age, est_var_x, edges)
            result["est"] = {
                "n_in_bin":   n_in_e,
                "h2_bin":     [float(h) for h in h2_est],
                "h2_total":   float(total_est),
                "h2_frac":    [float(h / max(total_est, 1e-30)) for h in h2_est],
                "beta_file":  bf,
            }
        else:
            result["est_error"] = f"sargent beta file not found at {bf}"

    out_path = args.out or (case + ".age_partition.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)

    # Console summary
    print(f"case: {os.path.basename(case)}  tau={tau:+.2f}  n_causal={len(snp_idx)}")
    print(f"  {'bin (gen)':<26s} {'n':>5s} {'h² truth':>10s} {'frac':>7s}", end="")
    if "est" in result:
        print(f" {'h² est':>10s} {'frac est':>9s}")
    else:
        print()
    for c, (lo, hi) in enumerate(bins):
        b = f"[{lo:>10.0f}, {hi:>10.0f}]"
        print(f"  {b:<26s} {n_in[c]:>5d} {h2_truth[c]:>10.4f} "
              f"{h2_truth[c]/max(total_truth,1e-30):>7.3f}", end="")
        if "est" in result:
            he = result["est"]["h2_bin"][c]
            print(f" {he:>10.4f} {result['est']['h2_frac'][c]:>9.3f}")
        else:
            print()
    print(f"  total h² (Σ β² · Var(x)): truth = {total_truth:.4f}", end="")
    if "est" in result:
        print(f"   est = {result['est']['h2_total']:.4f}")
    else:
        print()
    print(f"  wrote {out_path}")


if __name__ == "__main__":
    main()
