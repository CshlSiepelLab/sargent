#!/usr/bin/env python3
"""
overnight_aggregate.py — read all predict.json files from the overnight sweep
and produce a clean markdown summary with mean ± SD across seeds.
"""

import argparse
import json
import os
import sys
from collections import defaultdict

import numpy as np


def load_predict_json(path):
    if not os.path.exists(path):
        return None
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:
        return None


def fmt_mean_sd(values, prec=3):
    if not values:
        return "-"
    arr = np.asarray(values, float)
    arr = arr[~np.isnan(arr)]
    if len(arr) == 0:
        return "nan"
    if len(arr) == 1:
        return f"{arr[0]:+.{prec}f}"
    return f"{arr.mean():+.{prec}f} ± {arr.std(ddof=1):.{prec}f}"


def collect_r2(records, method):
    return [r[method]["r2"] for r in records if r is not None and method in r]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="/Users/acs/test-sargent/overnight")
    args = ap.parse_args()

    print("# Overnight sweep results")
    print()
    print(f"Generated from `{args.base}`.")
    print()

    # ===== Multi-seed bed-coverage sweep =====
    print("## Multi-seed bed-coverage sweep")
    print()
    print("Balanced h² case (h²_fixed = h²_random = 0.40, σ²_target = 0.20). "
          "n_train=200, n_test=100, length=100kb, regime=sparse.")
    print()

    f_values = [1.0, 0.6, 0.3, 0.1, 0.03]
    seeds_ms = [42, 43, 44, 45, 46]

    seq_by_f = defaultdict(list)
    em_by_f  = defaultdict(list)
    seq_var_by_f = defaultdict(list)  # (h2f_inf, h2r_inf, s2_inf) tuples
    em_var_by_f  = defaultdict(list)

    for seed in seeds_ms:
        for f in f_values:
            d = os.path.join(args.base, "multiseed", f"s{seed}_f{f}")
            seq = load_predict_json(os.path.join(d, "seq.predict.json"))
            em  = load_predict_json(os.path.join(d, "em.predict.json"))
            if seq: seq_by_f[f].append(seq)
            if em:  em_by_f[f].append(em)

    print("### Prediction R² (mean ± SD across 5 seeds)")
    print()
    print("| f | BayesR | GBLUP | ARG-BLUP-uni | Seq SARGENT | **Joint EM SARGENT** |")
    print("|---:|:---:|:---:|:---:|:---:|:---:|")
    for f in f_values:
        seq = seq_by_f.get(f, [])
        em  = em_by_f.get(f, [])
        b_seq  = collect_r2(seq, "bayesR")
        g_seq  = collect_r2(seq, "gblup")
        a_seq  = collect_r2(seq, "arg_blup_uniform")
        s_seq  = collect_r2(seq, "sargent")
        s_em   = collect_r2(em,  "sargent")
        print(f"| {f:.2f} | {fmt_mean_sd(b_seq)} | {fmt_mean_sd(g_seq)} | "
              f"{fmt_mean_sd(a_seq)} | {fmt_mean_sd(s_seq)} | "
              f"**{fmt_mean_sd(s_em)}** |")
    print()
    print("(Note: BayesR/GBLUP/ARG-BLUP-uniform are computed inside each "
          "predict.py run; the same values appear under seq-rows and em-rows "
          "because they don't depend on the SARGENT fit. Reported once.)")
    print()

    # Per-seed individual values for the most informative f
    print("### Per-seed sequential vs joint EM SARGENT R²")
    print()
    print("| seed | f=1.00 (seq / em) | f=0.60 (seq / em) | f=0.30 (seq / em) | f=0.10 (seq / em) | f=0.03 (seq / em) |")
    print("|---|---|---|---|---|---|")
    for seed in seeds_ms:
        row = [str(seed)]
        for f in f_values:
            d = os.path.join(args.base, "multiseed", f"s{seed}_f{f}")
            seq = load_predict_json(os.path.join(d, "seq.predict.json"))
            em  = load_predict_json(os.path.join(d, "em.predict.json"))
            s_r2 = seq.get("sargent", {}).get("r2") if seq else None
            e_r2 = em.get("sargent", {}).get("r2") if em else None
            s_str = f"{s_r2:+.3f}" if s_r2 is not None else "-"
            e_str = f"{e_r2:+.3f}" if e_r2 is not None else "-"
            row.append(f"{s_str} / {e_str}")
        print("| " + " | ".join(row) + " |")
    print()

    # ===== h² sweep at f=0.5 =====
    print("## h² sweep at f=0.5")
    print()
    print("All cases: n_train=200, n_test=100, length=100kb, bed_coverage=0.5, "
          "regime=sparse, total h²=0.8.")
    print()

    case_names = ["fixed_only", "fixed_heavy", "balanced", "rand_heavy", "rand_only"]
    seeds_h2 = [42, 43, 44]

    print("### Prediction R² (mean ± SD across 3 seeds)")
    print()
    print("| Case | h²_f / h²_r | BayesR | GBLUP | ARG-BLUP-uni | Seq SARGENT | **Joint EM SARGENT** |")
    print("|---|---|:---:|:---:|:---:|:---:|:---:|")
    h2_targets = {
        "fixed_only":  ("0.80", "0.00"),
        "fixed_heavy": ("0.60", "0.20"),
        "balanced":    ("0.40", "0.40"),
        "rand_heavy":  ("0.20", "0.60"),
        "rand_only":   ("0.00", "0.80"),
    }
    for case_name in case_names:
        seq, em = [], []
        for seed in seeds_h2:
            d = os.path.join(args.base, "h2sweep", f"s{seed}_{case_name}")
            s = load_predict_json(os.path.join(d, "seq.predict.json"))
            e = load_predict_json(os.path.join(d, "em.predict.json"))
            if s: seq.append(s)
            if e: em.append(e)
        h2f, h2r = h2_targets[case_name]
        print(f"| {case_name} | {h2f} / {h2r} | "
              f"{fmt_mean_sd(collect_r2(seq, 'bayesR'))} | "
              f"{fmt_mean_sd(collect_r2(seq, 'gblup'))} | "
              f"{fmt_mean_sd(collect_r2(seq, 'arg_blup_uniform'))} | "
              f"{fmt_mean_sd(collect_r2(seq, 'sargent'))} | "
              f"**{fmt_mean_sd(collect_r2(em, 'sargent'))}** |")
    print()

    # ===== Status =====
    print("## Status")
    print()
    n_ms_total = len(seeds_ms) * len(f_values)
    n_h2_total = len(seeds_h2) * len(case_names)
    n_ms_seq = sum(len(seq_by_f[f]) for f in f_values)
    n_ms_em  = sum(len(em_by_f[f])  for f in f_values)
    print(f"- Multi-seed bed-coverage: {n_ms_seq}/{n_ms_total} sequential, "
          f"{n_ms_em}/{n_ms_total} joint EM completed.")
    # Count h2 sweep
    n_h2_seq = 0; n_h2_em = 0
    for seed in seeds_h2:
        for case_name in case_names:
            d = os.path.join(args.base, "h2sweep", f"s{seed}_{case_name}")
            if os.path.exists(os.path.join(d, "seq.predict.json")): n_h2_seq += 1
            if os.path.exists(os.path.join(d, "em.predict.json")):  n_h2_em  += 1
    print(f"- h² sweep at f=0.5: {n_h2_seq}/{n_h2_total} sequential, "
          f"{n_h2_em}/{n_h2_total} joint EM completed.")


if __name__ == "__main__":
    main()
