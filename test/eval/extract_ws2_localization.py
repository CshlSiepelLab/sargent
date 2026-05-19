#!/usr/bin/env python3
"""
extract_ws2_localization.py — compute ws² localization metrics for both
sequential SARGENT and joint EM SARGENT outputs across all overnight cases.

For each case:
  - load true_ws2.tsv (true per-tree ws²)
  - load .w2.tsv (rho_active and E_ws2 columns)
  - compute AUC of rho_active for classifying truly-active trees
  - compute AUC of E_ws2 for the same classification
  - compute log-ws² correlation among truly-active trees
  - precision/recall at rho_active >= 0.5
"""

import json
import os
import sys

import numpy as np


def roc_auc(scores, labels):
    scores = np.asarray(scores, float)
    labels = np.asarray(labels, int)
    pos = labels == 1
    neg = labels == 0
    n_pos, n_neg = int(pos.sum()), int(neg.sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(scores)
    ranks = np.empty_like(order, dtype=float)
    sorted_scores = scores[order]
    i = 0
    while i < len(sorted_scores):
        j = i
        while j + 1 < len(sorted_scores) and sorted_scores[j + 1] == sorted_scores[i]:
            j += 1
        avg = (i + j) / 2.0 + 1
        ranks[order[i:j + 1]] = avg
        i = j + 1
    sum_ranks_pos = ranks[pos].sum()
    u = sum_ranks_pos - n_pos * (n_pos + 1) / 2.0
    return float(u / (n_pos * n_neg))


def pearson(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.std() == 0 or b.std() == 0: return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def load_true_ws2(prefix):
    rows = np.loadtxt(prefix + ".true_ws2.tsv", skiprows=1, usecols=(0, 1))
    if rows.ndim == 1: rows = rows[np.newaxis, :]
    idx  = rows[:, 0].astype(int)
    vals = rows[:, 1]
    out = np.zeros(int(idx.max()) + 1)
    out[idx] = vals
    return out


def load_w2(out_prefix):
    """Return (rho_active, E_ws2) arrays from a sargent output."""
    with open(out_prefix + ".w2.tsv") as f:
        header = f.readline().strip().split("\t")
    rho_cols = [i for i, h in enumerate(header) if h.startswith("rho_")]
    if not rho_cols:
        return None, None
    last_rho = rho_cols[-1]
    data = np.loadtxt(out_prefix + ".w2.tsv", skiprows=1)
    if data.ndim == 1: data = data[np.newaxis, :]
    return data[:, last_rho], data[:, -1]


def metrics_for_case(case_prefix, sargent_prefix):
    if not os.path.exists(sargent_prefix + ".w2.tsv"):
        return None
    if not os.path.exists(case_prefix + ".true_ws2.tsv"):
        return None
    true_ws2 = load_true_ws2(case_prefix)
    rho, e_ws2 = load_w2(sargent_prefix)
    if rho is None:
        return None
    n_true = len(true_ws2)
    n_inf  = len(rho)
    if n_true != n_inf:
        return None   # tree-count mismatch (shouldn't happen for train)

    labels = (true_ws2 > 0).astype(int)
    if labels.sum() == 0:    # no active trees (fixed_only case)
        return {"n_active_true": 0, "auc_rho": float("nan"),
                "auc_ws2": float("nan"), "precision_05": float("nan"),
                "recall_05": float("nan"), "logws2_corr_active": float("nan"),
                "n_called_active": int((rho >= 0.5).sum())}

    auc_rho = roc_auc(rho, labels)
    auc_ws2 = roc_auc(e_ws2, labels)

    # log(ws²) correlation among truly-active trees
    mask = labels == 1
    if mask.sum() >= 2:
        ltrue = np.log(true_ws2[mask])
        linf  = np.log(np.maximum(e_ws2[mask], 1e-300))
        corr_active = pearson(ltrue, linf)
    else:
        corr_active = float("nan")

    called = (rho >= 0.5).astype(int)
    tp = int(((called == 1) & (labels == 1)).sum())
    fp = int(((called == 1) & (labels == 0)).sum())
    fn = int(((called == 0) & (labels == 1)).sum())
    prec = tp / max(tp + fp, 1) if (tp + fp) > 0 else float("nan")
    rec  = tp / max(tp + fn, 1) if labels.sum() > 0 else float("nan")

    return {
        "n_active_true": int(labels.sum()),
        "n_called_active": int(called.sum()),
        "auc_rho":  auc_rho,
        "auc_ws2":  auc_ws2,
        "precision_05": prec,
        "recall_05":    rec,
        "logws2_corr_active": corr_active,
    }


def fmt(x, prec=3):
    if x is None: return "-"
    if isinstance(x, float):
        if x != x: return "nan"
        return f"{x:+.{prec}f}"
    return str(x)


def aggregate(values_list, prec=3):
    vs = [v for v in values_list if v is not None and not (isinstance(v, float) and v != v)]
    if not vs:
        return "-"
    a = np.asarray(vs, float)
    if len(a) == 1:
        return f"{a[0]:+.{prec}f}"
    return f"{a.mean():+.{prec}f} ± {a.std(ddof=1):.{prec}f}"


def main():
    base = "/Users/acs/test-sargent/overnight"

    # Multi-seed bed-coverage (balanced h2, 5 seeds, 5 f levels)
    print("## Debug — multi-seed ws² localization (5 seeds per f)")
    print()
    print("| f | method | AUC_rho | AUC_ws² | precision@0.5 | recall@0.5 | log(ws²)_corr_active | n_active_true | n_called |")
    print("|---:|:---:|:---:|:---:|:---:|:---:|:---:|---:|---:|")
    for f in ["1.0", "0.6", "0.3", "0.1", "0.03"]:
        for method in ["seq", "em"]:
            recs = []
            for seed in [42, 43, 44, 45, 46]:
                d = os.path.join(base, "multiseed", f"s{seed}_f{f}")
                m = metrics_for_case(os.path.join(d, "case"), os.path.join(d, method))
                if m: recs.append(m)
            if not recs:
                print(f"| {f} | {method} | - | - | - | - | - | - | - |")
                continue
            print(f"| {f} | {method} | "
                  f"{aggregate([r['auc_rho']  for r in recs])} | "
                  f"{aggregate([r['auc_ws2']  for r in recs])} | "
                  f"{aggregate([r['precision_05'] for r in recs])} | "
                  f"{aggregate([r['recall_05']    for r in recs])} | "
                  f"{aggregate([r['logws2_corr_active'] for r in recs])} | "
                  f"{np.mean([r['n_active_true']    for r in recs]):.1f} | "
                  f"{np.mean([r['n_called_active']  for r in recs]):.1f} |")
    print()

    # h2 sweep at f=0.5
    print("## Debug — h² sweep at f=0.5 ws² localization (3 seeds per case)")
    print()
    print("| case | method | AUC_rho | AUC_ws² | precision@0.5 | recall@0.5 | log(ws²)_corr_active | n_active_true | n_called |")
    print("|---|:---:|:---:|:---:|:---:|:---:|:---:|---:|---:|")
    for case_name in ["fixed_only", "fixed_heavy", "balanced", "rand_heavy", "rand_only"]:
        for method in ["seq", "em"]:
            recs = []
            for seed in [42, 43, 44]:
                d = os.path.join(base, "h2sweep", f"s{seed}_{case_name}")
                m = metrics_for_case(os.path.join(d, "case"), os.path.join(d, method))
                if m: recs.append(m)
            if not recs:
                print(f"| {case_name} | {method} | - | - | - | - | - | - | - |")
                continue
            print(f"| {case_name} | {method} | "
                  f"{aggregate([r['auc_rho']  for r in recs])} | "
                  f"{aggregate([r['auc_ws2']  for r in recs])} | "
                  f"{aggregate([r['precision_05'] for r in recs])} | "
                  f"{aggregate([r['recall_05']    for r in recs])} | "
                  f"{aggregate([r['logws2_corr_active'] for r in recs])} | "
                  f"{np.mean([r['n_active_true']    for r in recs]):.1f} | "
                  f"{np.mean([r['n_called_active']  for r in recs]):.1f} |")
    print()

    # Medium (single seed)
    print("## Medium — ws² localization (single seed=42)")
    print()
    print("| f | method | AUC_rho | AUC_ws² | precision@0.5 | recall@0.5 | log(ws²)_corr_active | n_active_true | n_called |")
    print("|---:|:---:|:---:|:---:|:---:|:---:|:---:|---:|---:|")
    for f in ["1.0", "0.3", "0.1"]:
        for method in ["seq", "em"]:
            d = os.path.join(base, "medium", f"f{f}")
            m = metrics_for_case(os.path.join(d, "case"), os.path.join(d, method))
            if not m:
                print(f"| {f} | {method} | - | - | - | - | - | - | - |")
                continue
            def s(v, prec=3):
                if isinstance(v, float):
                    return "nan" if v != v else f"{v:+.{prec}f}"
                return str(v)
            print(f"| {f} | {method} | {s(m['auc_rho'])} | {s(m['auc_ws2'])} | "
                  f"{s(m['precision_05'])} | {s(m['recall_05'])} | "
                  f"{s(m['logws2_corr_active'])} | {m['n_active_true']} | "
                  f"{m['n_called_active']} |")


if __name__ == "__main__":
    main()
