#!/usr/bin/env python3
"""
oracle_ws2.py — Phase 10 oracle diagnostic for per-segment ws².

Given a PURE-RANDOM dataset from simulate.py (--h2-snp 0, --n-test > 0, and a
--ws2-architecture), this compares out-of-sample ARG-BLUP prediction under three
choices of per-tree weights, all sharing the same kinship machinery so that only
the weight vector differs:

  - uniform   : wₛ = 1/S            (ARG-BLUP-uniform)
  - oracle    : wₛ = true ws²        (from <case>.full.true_ws2.tsv)
  - estimated : wₛ = sargent ŵ²      (from <sargent>.w2.tsv, --arg-only)

It also reports how well the estimated profile recovers the true one (local-h²
recovery): Pearson/Spearman correlation and the AUC for detecting active
(in-block) trees.

Decision rule (Phase 10):
  oracle ≈ uniform even on the sparse architecture  → fundamental ceiling.
  oracle ≫ uniform but estimated ≈ uniform          → estimation gap → Phase 11.

Usage:
    python oracle_ws2.py --case-prefix DIR/CASE --sargent-prefix DIR/OUT \
        [--out RESULT.json]
"""

import argparse
import json
import os
import sys

import numpy as np
import tskit

# Reuse the tested predict.py machinery (same directory).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from predict import (read_phen, load_ws2, hap_to_indiv_2d, fit_blup,
                     predict_blup, metrics, map_train_to_full,
                     center_kinship)


def streaming_multi_kinship(ts, sample_nodes_a, weight_vectors,
                            sample_nodes_b=None):
    """Σ_s w[s]·G_s for SEVERAL weight vectors in a single pass over trees.

    weight_vectors: dict {name: array length ts.num_trees}.
    Returns: dict {name: (n_a × n_b) array}.

    Same G_s convention as predict.streaming_combined_kinship:
      G_s[i,j] = 0.5 · Σ_b bl(b) · I_b(i) · I_b(j).
    """
    if sample_nodes_b is None:
        sample_nodes_b = sample_nodes_a
    n_a, n_b = len(sample_nodes_a), len(sample_nodes_b)
    rank_a = {nd: r for r, nd in enumerate(sample_nodes_a)}
    rank_b = {nd: r for r, nd in enumerate(sample_nodes_b)}
    set_a, set_b = set(sample_nodes_a), set(sample_nodes_b)

    Ks = {name: np.zeros((n_a, n_b), dtype=np.float64) for name in weight_vectors}

    for s, tree in enumerate(ts.trees()):
        ws = {name: float(w[s]) for name, w in weight_vectors.items()}
        if all(v == 0.0 for v in ws.values()):
            continue
        for nd in tree.nodes(order="preorder"):
            parent = tree.parent(nd)
            if parent == tskit.NULL:
                continue
            bl = tree.branch_length(nd)
            below = list(tree.samples(nd))
            ia = [rank_a[x] for x in below if x in set_a]
            ib = [rank_b[x] for x in below if x in set_b]
            if not ia or not ib:
                continue
            ia_arr = np.array(ia, dtype=int)
            ib_arr = np.array(ib, dtype=int)
            block = np.ix_(ia_arr, ib_arr)
            for name, w in ws.items():
                if w != 0.0:
                    Ks[name][block] += 0.5 * w * bl
    return Ks


def auc(scores, labels):
    """ROC-AUC via the Mann-Whitney U statistic (no sklearn dependency)."""
    scores = np.asarray(scores, float)
    labels = np.asarray(labels, int)
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores), float)
    ranks[order] = np.arange(1, len(scores) + 1)
    # average ranks for ties
    _, inv, counts = np.unique(scores, return_inverse=True, return_counts=True)
    tie_mean = np.zeros(len(counts))
    np.add.at(tie_mean, inv, ranks)
    tie_mean /= counts
    ranks = tie_mean[inv]
    r_pos = ranks[labels == 1].sum()
    u = r_pos - len(pos) * (len(pos) + 1) / 2.0
    return float(u / (len(pos) * len(neg)))


def spearman(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ra = np.argsort(np.argsort(a))
    rb = np.argsort(np.argsort(b))
    if ra.std() == 0 or rb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def idealized_loc_auc(ts, train_nodes, r_tilde, active_labels):
    """Localization CEILING: AUC of the idealized per-tree signal A_s/tr(G_s),
    computed at the true/oracle V via r_tilde = V^{-1}(y - ȳ), for detecting
    active trees.  This is the best a score-based localizer could do with no
    estimation error in V — it upper-bounds the real estimator's AUC.

    Per-branch decomposition (individual-level, matching predict.py's G_s):
      branch b contributes 0.25·bl·c_i·c_j to G_s[i,j], where c_i = #{i's
      haplotypes below b}.  Hence
        A_s   = r̃ᵀ G_s r̃ = Σ_b 0.25·bl·(Σ_i c_i r̃_i)²,
        tr Gs = Σ_b 0.25·bl·Σ_i c_i².
    Detector_s = A_s / tr(G_s) (scale-normalised so deep trees aren't favoured).
    """
    rank_of = {nd: r for r, nd in enumerate(train_nodes)}   # rank//2 = individual
    set_tr  = set(train_nodes)
    S = ts.num_trees
    detector = np.full(S, np.nan)
    for s, tree in enumerate(ts.trees()):
        A = 0.0
        trG = 0.0
        for nd in tree.nodes(order="preorder"):
            if tree.parent(nd) == tskit.NULL:
                continue
            bl = tree.branch_length(nd)
            counts = {}
            for h in tree.samples(nd):
                if h in set_tr:
                    i = rank_of[h] // 2
                    counts[i] = counts.get(i, 0) + 1
            if not counts:
                continue
            s1 = sum(c * r_tilde[i] for i, c in counts.items())
            s2 = sum(c * c for c in counts.values())
            A   += 0.25 * bl * s1 * s1
            trG += 0.25 * bl * s2
        if trG > 0:
            detector[s] = A / trG
    ok = ~np.isnan(detector)
    return auc(detector[ok], np.asarray(active_labels)[ok])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case-prefix",    required=True)
    ap.add_argument("--sargent-prefix", required=True)
    ap.add_argument("--out",            default=None)
    args = ap.parse_args()

    case, out = args.case_prefix, args.sargent_prefix

    with open(case + ".params.json") as f:
        params = json.load(f)
    n_train, n_test, n_total = params["n_indiv"], params["n_test"], params["n_total"]
    arch = params.get("ws2_architecture")
    print(f"case={case}  arch={arch}  n_train={n_train} n_test={n_test}", flush=True)

    y_train = read_phen(case + ".phen")
    y_test  = read_phen(case + ".test.phen")

    ts_full = tskit.load(case + ".full.trees")
    ts_train_only = tskit.load(case + ".trees")
    S_full = ts_full.num_trees

    train_nodes = []
    for i in range(n_train):
        train_nodes.extend(ts_full.individual(i).nodes)
    test_nodes = []
    for i in range(n_train, n_total):
        test_nodes.extend(ts_full.individual(i).nodes)

    # --- weight vectors over the FULL tree sequence ---------------------------
    # uniform
    w_uniform = np.ones(S_full) / S_full

    # oracle (true ws², full-ARG indexed)
    true_path = case + ".full.true_ws2.tsv"
    if not os.path.exists(true_path):
        true_path = case + ".true_ws2.tsv"   # no-test fallback
    true_full = np.loadtxt(true_path, skiprows=1)[:, -1]
    if len(true_full) != S_full:
        sys.exit(f"ERROR: true ws² length {len(true_full)} != {S_full} full trees")
    w_oracle = true_full.copy()

    # estimated (sargent ŵ², train-tree indexed → map to full)
    w_est_full = None
    if os.path.exists(out + ".w2.tsv"):
        ws2_est = load_ws2(out + ".w2.tsv")
        if len(ws2_est) == ts_train_only.num_trees:
            idx_map = map_train_to_full(ts_train_only, ts_full)
            w_est_full = ws2_est[idx_map]
        elif len(ws2_est) == S_full:
            w_est_full = ws2_est
        else:
            print(f"  WARNING: ŵ² length {len(ws2_est)} matches neither train "
                  f"({ts_train_only.num_trees}) nor full ({S_full}) trees; "
                  f"skipping estimated arm.", flush=True)

    weight_vectors = {"uniform": w_uniform, "oracle": w_oracle}
    if w_est_full is not None:
        weight_vectors["estimated"] = w_est_full

    # --- kinship (single pass for train, single pass for cross) ---------------
    print(f"Computing kinship over {S_full} trees for {len(weight_vectors)} "
          f"weightings ...", flush=True)
    K_train_hap = streaming_multi_kinship(ts_full, train_nodes, weight_vectors)
    K_cross_hap = streaming_multi_kinship(ts_full, test_nodes, weight_vectors,
                                          train_nodes)

    results = {"case": case, "ws2_architecture": arch,
               "n_train": n_train, "n_test": n_test, "n_trees_full": S_full,
               "h2_arg_realised": params.get("h2_arg_realised")}

    print("\n=== OOS prediction (pure ARG-BLUP; only the weights differ) ===")
    uniform_fit = None
    for name in ["uniform", "oracle", "estimated"]:
        if name not in weight_vectors:
            continue
        K_tr = hap_to_indiv_2d(K_train_hap[name], n_train)
        K_cr = hap_to_indiv_2d(K_cross_hap[name], n_test, n_train)
        # Double-center to absorb the deep-coalescent constant offset (Block 1.3
        # audit; previously this was just mean-diag normalisation, which fit_blup
        # then biased h² upward by ~+0.25). Same centering applied to K_cr so
        # the BLUP prediction stays calibrated.
        K_tr, K_cr = center_kinship(K_tr, K_cr)
        s2, su2, alpha, mean_y = fit_blup(y_train, K_tr)
        y_hat = predict_blup(K_cr, alpha, mean_y)
        m = metrics(y_hat, y_test)
        m["h2_fit"] = su2 / (su2 + s2)
        results[name] = m
        if name == "uniform":
            uniform_fit = {"K": K_tr, "su2": su2, "s2": s2}
        print(f"  {name:10s}  R²={m['r2']:+.4f}  r={m['r']:+.4f}  "
              f"h2_fit={m['h2_fit']:.3f}")

    # paired gaps
    if "oracle" in results and "uniform" in results:
        results["oracle_minus_uniform_r2"] = results["oracle"]["r2"] - results["uniform"]["r2"]
    if "estimated" in results and "uniform" in results:
        results["estimated_minus_uniform_r2"] = results["estimated"]["r2"] - results["uniform"]["r2"]

    # --- localization CEILING: idealized per-segment score at the homogeneous
    # (uniform) null V.  Active trees carry excess signal beyond the genome-wide
    # average; the score r̃ᵀG_sr̃ / tr(G_s) (r̃ = V_uniform^{-1}(y-ȳ)) detects them.
    # Evaluated with the exact V (no Hutchinson/estimation noise), so it
    # upper-bounds any real estimator's localization AUC.
    active = (true_full > true_full.min() + 1e-12).astype(int)
    if uniform_fit is not None and active.sum() > 0:
        n = n_train
        V = uniform_fit["su2"] * uniform_fit["K"] + uniform_fit["s2"] * np.eye(n)
        r_tilde = np.linalg.solve(V, y_train - y_train.mean())
        ceil = idealized_loc_auc(ts_full, train_nodes, r_tilde, active)
        results["loc_ceiling_auc"] = ceil
        print(f"\n=== localization ceiling (idealized score at uniform-null V) ===")
        print(f"  AUC(active) ceiling = {ceil:.3f}   (n_active={int(active.sum())})")

    # --- local-h² recovery (true vs estimated profile) ------------------------
    if w_est_full is not None:
        active = (true_full > true_full.min() + 1e-12).astype(int)   # in-block trees
        rec = {
            "pearson":  float(np.corrcoef(true_full, w_est_full)[0, 1])
                        if true_full.std() > 0 and w_est_full.std() > 0 else float("nan"),
            "spearman": spearman(true_full, w_est_full),
            "auc_active": auc(w_est_full, active),
            "n_active": int(active.sum()),
        }
        results["ws2_recovery"] = rec
        print(f"\n=== ws² recovery (true vs estimated) ===")
        print(f"  pearson={rec['pearson']:+.3f}  spearman={rec['spearman']:+.3f}  "
              f"AUC(active)={rec['auc_active']:.3f}  n_active={rec['n_active']}")

    out_path = args.out or (out + ".oracle.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
