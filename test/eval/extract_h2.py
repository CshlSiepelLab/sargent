#!/usr/bin/env python3
"""
extract_h2.py — read sargent's β̂ and ŵ² outputs and recompute h²_fixed,
h²_random, σ²_inferred for each case in the overnight + medium sweeps.

For each case, h²_fixed = var(X_train · β̂) / var(y_train), and
h²_random ≈ tr(K_train) / n_train / var(y_train) where
K_train = Σ_s ŵ²[s] · G_s_train_train.
"""

import json
import os
import sys

import numpy as np
import tskit

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from predict import read_bed, streaming_combined_kinship, hap_to_indiv_2d


def compute_h2(case_prefix, sargent_prefix):
    """Return dict with h2_fixed_inf, h2_random_inf, sigma2_inf, plus targets."""
    with open(case_prefix + ".params.json") as f:
        params = json.load(f)
    n_train = params["n_indiv"]
    p_bed   = params["n_variants"]

    # Phenotype
    y_train = np.loadtxt(case_prefix + ".phen")
    var_y   = float(y_train.var())

    # Fixed-effect contribution: var(X β̂)
    h2_fixed_inf = float("nan")
    if os.path.exists(sargent_prefix + ".beta.tsv"):
        beta = np.loadtxt(sargent_prefix + ".beta.tsv", skiprows=1)
        if beta.ndim == 1:
            beta = beta[np.newaxis, :]
        beta_vals = beta[:, -1]
        # Read and standardise X_train (same as predict.py)
        dosage = read_bed(case_prefix, n_train, p_bed)
        # filter to polymorphic-in-train (matches sargent's filter)
        mn = dosage.min(axis=1)
        mx = dosage.max(axis=1)
        poly = mx > mn
        dosage = dosage[poly]
        beta_keep = beta_vals[:dosage.shape[0]]
        # standardise
        mu = dosage.mean(axis=1, keepdims=True)
        sd = dosage.std(axis=1, keepdims=True)
        sd[sd == 0] = 1.0
        X = (dosage - mu) / sd       # shape (p, n)
        Xb = X.T @ beta_keep
        h2_fixed_inf = float(Xb.var()) / max(var_y, 1e-15)

    # Random-effect contribution: tr(K_train) / n / var(y)
    h2_random_inf = float("nan")
    if os.path.exists(sargent_prefix + ".w2.tsv"):
        ws2 = np.loadtxt(sargent_prefix + ".w2.tsv", skiprows=1)
        if ws2.ndim == 1:
            ws2 = ws2[np.newaxis, :]
        ws2_vals = ws2[:, -1]
        ts_train = tskit.load(case_prefix + ".trees")
        # train sample nodes
        train_nodes = []
        for i in range(n_train):
            train_nodes.extend(ts_train.individual(i).nodes)
        K_hap = streaming_combined_kinship(ts_train, train_nodes, ws2_vals)
        K_ind = hap_to_indiv_2d(K_hap, n_train)
        tr_K = float(np.trace(K_ind))
        h2_random_inf = tr_K / n_train / max(var_y, 1e-15)

    # σ²_inferred ≈ var_y * (1 - h2_fixed_inf - h2_random_inf), but clamp
    s2_inf = var_y * max(0.0, 1.0 - h2_fixed_inf - h2_random_inf)

    return {
        "h2_fixed_target":  params["h2_fixed_target"],
        "h2_random_target": params["h2_random_target"],
        "h2_fixed_real":    params.get("h2_fixed_realised", float("nan")),
        "h2_random_real":   params.get("h2_random_realised", float("nan")),
        "h2_fixed_inf":     h2_fixed_inf,
        "h2_random_inf":    h2_random_inf,
        "sigma2_target":    params.get("sigma2", float("nan")),
        "sigma2_inf":       s2_inf,
        "var_y":            var_y,
    }


def main():
    base = "/Users/acs/test-sargent/overnight"

    # --- Debug multi-seed: average over seeds per f, per method ---
    print("## Debug — multi-seed bed-coverage sweep (balanced h², 5 seeds)")
    print("    target h²_fixed = 0.40, h²_random = 0.40 (sum 0.80; some realised drift)")
    print()
    print(f"{'f':>6} | {'method':>7} | {'h²_f_real':>10} | {'h²_r_real':>10} | "
          f"{'h²_f_inf':>10} | {'h²_r_inf':>10} | {'sum_inf':>8}")
    print("-" * 80)
    seeds = [42, 43, 44, 45, 46]
    for f in ["1.0", "0.6", "0.3", "0.1", "0.03"]:
        for method in ["seq", "em"]:
            r = []
            for seed in seeds:
                d = os.path.join(base, "multiseed", f"s{seed}_f{f}")
                case   = os.path.join(d, "case")
                sarg   = os.path.join(d, method)
                if not os.path.exists(sarg + ".beta.tsv"):
                    continue
                r.append(compute_h2(case, sarg))
            if not r:
                continue
            hfr = np.mean([x["h2_fixed_real"]  for x in r])
            hrr = np.mean([x["h2_random_real"] for x in r])
            hfi = np.mean([x["h2_fixed_inf"]   for x in r])
            hri = np.mean([x["h2_random_inf"]  for x in r])
            print(f"{f:>6} | {method:>7} | {hfr:>10.3f} | {hrr:>10.3f} | "
                  f"{hfi:>10.3f} | {hri:>10.3f} | {hfi+hri:>8.3f}")

    # --- Debug h² sweep: average over seeds per case, per method ---
    print()
    print("## Debug — h² sweep at f=0.5 (3 seeds each)")
    print()
    print(f"{'case':>14} | {'method':>7} | {'h²_f_real':>10} | {'h²_r_real':>10} | "
          f"{'h²_f_inf':>10} | {'h²_r_inf':>10} | {'sum_inf':>8}")
    print("-" * 90)
    for case_name in ["fixed_only", "fixed_heavy", "balanced", "rand_heavy", "rand_only"]:
        for method in ["seq", "em"]:
            r = []
            for seed in [42, 43, 44]:
                d = os.path.join(base, "h2sweep", f"s{seed}_{case_name}")
                case   = os.path.join(d, "case")
                sarg   = os.path.join(d, method)
                if not os.path.exists(sarg + ".beta.tsv"):
                    continue
                r.append(compute_h2(case, sarg))
            if not r:
                continue
            hfr = np.mean([x["h2_fixed_real"]  for x in r])
            hrr = np.mean([x["h2_random_real"] for x in r])
            hfi = np.mean([x["h2_fixed_inf"]   for x in r])
            hri = np.mean([x["h2_random_inf"]  for x in r])
            print(f"{case_name:>14} | {method:>7} | {hfr:>10.3f} | {hrr:>10.3f} | "
                  f"{hfi:>10.3f} | {hri:>10.3f} | {hfi+hri:>8.3f}")

    # --- Medium: single seed=42 ---
    print()
    print("## Medium — single seed=42 (balanced, n_train=300, L=300kb)")
    print()
    print(f"{'f':>6} | {'method':>7} | {'h²_f_real':>10} | {'h²_r_real':>10} | "
          f"{'h²_f_inf':>10} | {'h²_r_inf':>10} | {'sum_inf':>8}")
    print("-" * 80)
    for f in ["1.0", "0.3", "0.1"]:
        for method in ["seq", "em"]:
            d = os.path.join(base, "medium", f"f{f}")
            case   = os.path.join(d, "case")
            sarg   = os.path.join(d, method)
            if not os.path.exists(sarg + ".beta.tsv"):
                continue
            x = compute_h2(case, sarg)
            print(f"{f:>6} | {method:>7} | {x['h2_fixed_real']:>10.3f} | "
                  f"{x['h2_random_real']:>10.3f} | {x['h2_fixed_inf']:>10.3f} | "
                  f"{x['h2_random_inf']:>10.3f} | "
                  f"{x['h2_fixed_inf']+x['h2_random_inf']:>8.3f}")


if __name__ == "__main__":
    main()
