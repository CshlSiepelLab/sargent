#!/usr/bin/env python3
"""
predict.py — out-of-sample prediction for SARGENT.

Given a dataset produced by simulate.py with --n-test > 0, this script:
  1. Loads sargent's trained outputs (β̂, ŵ²[s], σ̂² implicit)
  2. Computes per-tree kinship matrices G_s from the full ARG (.full.trees)
  3. Predicts ŷ_test under four methods:
        a. BayesR-only      — ŷ = X_test β̂
        b. GBLUP            — single genome-wide GRM, h² fit by REML on training
        c. ARG-BLUP-uniform — Σ_s G_s with uniform ws², h² fit by REML
        d. Full SARGENT     — X_test β̂ + Σ_s ŵ²[s] G_s_test_train V_train^{-1} (y - X β̂)
  4. Reports OOS R²(ŷ, y) and Pearson r for each method.

Usage:
    python predict.py --case-prefix DIR/CASE --sargent-prefix DIR/OUT [--out RESULT.json]
"""

import argparse
import json
import os
import sys

import numpy as np
import tskit


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def read_phen(path):
    return np.loadtxt(path).astype(float)


def read_bed(prefix, n_indiv, n_variants):
    """Standardise dosages (mean 0, var 1) as sargent does internally."""
    bytes_per_snp = (n_indiv + 3) // 4
    dosage = np.zeros((n_variants, n_indiv), dtype=float)
    with open(prefix + ".bed", "rb") as f:
        magic = f.read(3)
        if magic != bytes([0x6c, 0x1b, 0x01]):
            raise ValueError(f"bad BED magic in {prefix}.bed")
        for j in range(n_variants):
            buf = f.read(bytes_per_snp)
            if len(buf) < bytes_per_snp:
                raise ValueError(f"BED truncated at SNP {j}")
            for i in range(n_indiv):
                code = (buf[i // 4] >> ((i % 4) * 2)) & 0x3
                if code == 0:
                    v = 0.0
                elif code == 2:
                    v = 1.0
                elif code == 3:
                    v = 2.0
                else:
                    v = float("nan")
                dosage[j, i] = v
    return dosage


def standardise(dosage_train, dosage_test):
    """Standardise both train and test using TRAIN mean/std (sargent convention)."""
    mean = dosage_train.mean(axis=1, keepdims=True)
    std  = dosage_train.std(axis=1, keepdims=True)
    std[std == 0] = 1.0
    return (dosage_train - mean) / std, (dosage_test - mean) / std


def load_beta(beta_tsv):
    """Read E_beta column from sargent's .beta.tsv output."""
    data = np.loadtxt(beta_tsv, skiprows=1)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return data[:, -1]


def load_ws2(w2_tsv):
    data = np.loadtxt(w2_tsv, skiprows=1)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return data[:, -1]


# ---------------------------------------------------------------------------
# Per-tree kinship from a tskit tree sequence
# ---------------------------------------------------------------------------

def streaming_combined_kinship(ts, sample_nodes_a, weights_per_tree,
                               sample_nodes_b=None):
    """Compute Σ_s weights_per_tree[s] · G_s without storing per-tree kinship.

    Kinship between samples i, j on tree s is 0.5 · Σ_b bl(b) · I_b(i) · I_b(j)
    where I_b(k)=1 iff sample k is in the subtree below branch b.

    weights_per_tree: array of length ts.num_trees (or None for uniform sum).

    Returns: (n_a × n_b) numpy array containing the weighted sum.
    """
    if sample_nodes_b is None:
        sample_nodes_b = sample_nodes_a
    n_a, n_b = len(sample_nodes_a), len(sample_nodes_b)
    rank_a = {nd: r for r, nd in enumerate(sample_nodes_a)}
    rank_b = {nd: r for r, nd in enumerate(sample_nodes_b)}
    sample_set_a = set(sample_nodes_a)
    sample_set_b = set(sample_nodes_b)

    K = np.zeros((n_a, n_b), dtype=np.float64)
    have_weights = weights_per_tree is not None

    for s, tree in enumerate(ts.trees()):
        w = float(weights_per_tree[s]) if have_weights else 1.0
        if w == 0.0:
            continue
        scale = 0.5 * w
        for nd in tree.nodes(order="preorder"):
            parent = tree.parent(nd)
            if parent == tskit.NULL:
                continue
            bl = tree.branch_length(nd)
            below = list(tree.samples(nd))
            ia = [rank_a[x] for x in below if x in sample_set_a]
            ib = [rank_b[x] for x in below if x in sample_set_b]
            if not ia or not ib:
                continue
            ia_arr = np.array(ia, dtype=int)
            ib_arr = np.array(ib, dtype=int)
            K[np.ix_(ia_arr, ib_arr)] += scale * bl

    return K


def hap_to_indiv_2d(K_hap, n_a, n_b=None):
    """Collapse haplotype-level (2n_a x 2n_b) to individual-level (n_a x n_b).

    Each individual contributes two haplotypes; K_ind[i,j] = (1/2) Σ over the
    four (h_i, h_j) pairs.
    """
    if n_b is None:
        n_b = n_a
    return K_hap.reshape(n_a, 2, n_b, 2).sum(axis=(1, 3)) * 0.5


def center_kinship(K_train, K_cross=None):
    """Double-center K_train (H K H, H = I − 11ᵀ/n) and apply the consistent
    centering to K_cross if supplied. Then mean-diag = 1 normalisation.

    Why: the aggregate ARG kinship has a rank-1 constant offset
    (deep-coalescent shared-branch baseline; |cos(v₁,1)|=0.998 and λ₁/tr ≈ 0.65
    on test cases) that fit_blup otherwise absorbs into σ_u², inflating h̅²_fit
    by ≈+0.25. Removing the all-ones eigenmode brings the bias to ≈+0.05, the
    same floor as SARGENT's per-tree C-code estimator and as fit_blup on a
    SNP-GRM. See test/eval/fit_blup_audit.md.

    Args:
        K_train: (n_train × n_train) symmetric train block.
        K_cross: (n_cross × n_train) optional cross block (test → train).

    Returns:
        K_train_centered, K_cross_centered (None if K_cross is None).
        The same row/column statistics are used for both, so the cross
        block is on the same scale as the train block — required for
        BLUP prediction.
    """
    K_train = np.asarray(K_train, float)
    row_mean_train  = K_train.mean(axis=1, keepdims=True)   # (n_train, 1)
    col_mean_train  = K_train.mean(axis=0, keepdims=True)   # (1, n_train)
    grand_mean_train = float(K_train.mean())
    K_train_c = K_train - row_mean_train - col_mean_train + grand_mean_train
    d = float(np.mean(np.diag(K_train_c)))
    if d > 0:
        K_train_c = K_train_c / d
    if K_cross is None:
        return K_train_c, None
    K_cross = np.asarray(K_cross, float)
    row_mean_cross = K_cross.mean(axis=1, keepdims=True)    # (n_cross, 1)
    K_cross_c = K_cross - row_mean_cross - col_mean_train + grand_mean_train
    if d > 0:
        K_cross_c = K_cross_c / d
    return K_train_c, K_cross_c


# ---------------------------------------------------------------------------
# Predictors
# ---------------------------------------------------------------------------

def predict_bayesR(X_test, beta, mean_y):
    """Predict via fixed effects only, adding back the training mean."""
    return X_test.T @ beta + mean_y


def fit_blup(y_train, K_train_train, h2_init=0.5, max_iter=50, tol=1e-5):
    """Simple REML for a single random-effect model y = u + ε with K_train_train.

    Estimates h2 and σ² by iterative variance-component updates.
    Returns: (sigma2, sigma_u2, alpha)
      where alpha = (K + λI)^{-1} y_centered, λ = σ²/σ_u² (i.e., BLUP solution).
    """
    n = len(y_train)
    y_c = y_train - y_train.mean()
    var_y = y_c.var()

    sigma_u2 = h2_init * var_y
    sigma2   = (1 - h2_init) * var_y

    for _ in range(max_iter):
        V = sigma_u2 * K_train_train + sigma2 * np.eye(n)
        try:
            L = np.linalg.cholesky(V)
        except np.linalg.LinAlgError:
            V = V + 1e-6 * np.eye(n)
            L = np.linalg.cholesky(V)
        Vinv_y = np.linalg.solve(L.T, np.linalg.solve(L, y_c))
        # gradient-style updates (REML-ish)
        K_Vinv_y = K_train_train @ Vinv_y
        # tr(V^{-1} K)
        Vinv = np.linalg.solve(L.T, np.linalg.solve(L, np.eye(n)))
        tr_VinvK = np.einsum("ij,ji->", Vinv, K_train_train)
        tr_Vinv  = np.trace(Vinv)
        # Multiplicative update (REML EM style)
        a_u = float(K_Vinv_y @ Vinv_y)
        b_u = float(tr_VinvK)
        a_e = float(Vinv_y @ Vinv_y)
        b_e = float(tr_Vinv)
        new_su2 = sigma_u2 * a_u / max(b_u, 1e-30)
        new_s2  = sigma2   * a_e / max(b_e, 1e-30)
        # Clamp
        new_su2 = max(new_su2, 1e-8 * var_y)
        new_s2  = max(new_s2, 1e-8 * var_y)
        if abs(new_su2 - sigma_u2) / max(sigma_u2, 1e-30) < tol and \
           abs(new_s2 - sigma2) / max(sigma2, 1e-30) < tol:
            sigma_u2, sigma2 = new_su2, new_s2
            break
        sigma_u2, sigma2 = new_su2, new_s2

    # Solve for BLUP alpha: u = K alpha, where alpha = (K + (σ²/σ_u²) I)^{-1} y_c
    lam = sigma2 / max(sigma_u2, 1e-30)
    A = K_train_train + lam * np.eye(n)
    alpha = np.linalg.solve(A, y_c)
    return sigma2, sigma_u2, alpha, y_train.mean()


def predict_blup(K_test_train, alpha, mean_y):
    return K_test_train @ alpha + mean_y


def predict_sargent(X_test, beta, K_train, K_cross,
                    y_train, X_train):
    """Full SARGENT prediction with pre-combined kinship.

    Fits σ², σ_u² on the residual r = y - X β̂ by REML (same machinery as
    GBLUP / ARG-BLUP-uniform) so that prediction quality across methods is
    comparable on equal footing — only the kinship matrix differs.

    Centers K (and the cross block consistently) before fitting; see
    `center_kinship` for rationale.

    y_test_hat = X_test β̂ + mean_r + K_cross · α,  α = (K_train + λI)^{-1} (r − mean_r)
    """
    r = y_train - X_train.T @ beta
    K_train_c, K_cross_c = center_kinship(K_train, K_cross)
    sigma2, _sigma_u2, alpha, mean_r = fit_blup(r, K_train_c)
    u_test = K_cross_c @ alpha
    y_test_hat = X_test.T @ beta + mean_r + u_test
    return y_test_hat, sigma2


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def metrics(y_hat, y_true):
    """Out-of-sample prediction metrics.

    Reports two R² flavours because they answer different questions, and they
    diverge sharply in the cross-ancestry regime:

      r2       — raw coefficient of determination, 1 − SS_res/SS_tot. NOT
                 scale-invariant: it penalises any mean/variance miscalibration
                 of ŷ. The fixed-effect (BayesR) PRS has no built-in calibration,
                 so when a SNP is rare in the training population but common in
                 the test population, standardising test dosages by *train*
                 mean/std inflates ŷ's variance and r2 can go sharply negative
                 (observed down to −3 in Block 1.1) even when ŷ is positively
                 correlated with y. The BLUP arms self-calibrate via shrinkage,
                 so this asymmetry makes raw r2 an unfair cross-method metric.

      r2_calib — scale-invariant R² = squared Pearson correlation = the r2 you
                 would obtain after the optimal linear recalibration ŷ → a+bŷ.
                 This isolates *discrimination* from *calibration* and is the
                 field-standard cross-ancestry PRS metric. Use this for the
                 SARGENT-vs-baseline comparison; keep raw r2 for transparency.
    """
    y_hat, y_true = np.asarray(y_hat, float), np.asarray(y_true, float)
    ss_res = np.sum((y_true - y_hat) ** 2)
    ss_tot = np.sum((y_true - y_true.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    if y_hat.std() > 0 and y_true.std() > 0:
        r = float(np.corrcoef(y_hat, y_true)[0, 1])
    else:
        r = float("nan")
    r2_calib = r * r if r == r else float("nan")   # r==r guards against NaN
    return {"r": r, "r2": r2, "r2_calib": r2_calib, "mse": ss_res / len(y_true)}


# ---------------------------------------------------------------------------
# Tree-index mapping: train.trees → full.trees by position
# ---------------------------------------------------------------------------

def map_train_to_full(ts_train, ts_full):
    """For each tree in ts_full, return the index of the train tree containing its midpoint."""
    train_breaks = list(ts_train.breakpoints())
    out = np.zeros(ts_full.num_trees, dtype=int)
    for t, tree in enumerate(ts_full.trees()):
        mid = 0.5 * (tree.interval.left + tree.interval.right)
        idx = max(0, min(ts_train.num_trees - 1,
                         np.searchsorted(train_breaks, mid, side="right") - 1))
        out[t] = idx
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case-prefix",    required=True,
                    help="path prefix used by simulate.py (e.g. /path/case)")
    ap.add_argument("--sargent-prefix", required=True,
                    help="prefix of sargent's output files (.beta.tsv, .w2.tsv)")
    ap.add_argument("--out",            default=None,
                    help="output JSON path (default: <sargent-prefix>.predict.json)")
    args = ap.parse_args()

    case = args.case_prefix
    out  = args.sargent_prefix

    with open(case + ".params.json") as f:
        params = json.load(f)
    n_train = params["n_indiv"]
    n_test  = params["n_test"]
    n_total = params["n_total"]
    p       = params["n_variants"]

    print(f"Loading data: n_train={n_train}, n_test={n_test}, p={p}", flush=True)

    # Phenotypes
    y_train = read_phen(case + ".phen")
    y_test  = read_phen(case + ".test.phen")
    assert len(y_train) == n_train and len(y_test) == n_test

    # Genotypes
    dosage_train = read_bed(case,           n_train, p)
    dosage_test  = read_bed(case + ".test", n_test,  p)

    # Filter to SNPs that are polymorphic in the TRAINING set (matches sargent's filter)
    train_min = dosage_train.min(axis=1)
    train_max = dosage_train.max(axis=1)
    poly_mask = train_max > train_min
    n_poly    = int(poly_mask.sum())
    if n_poly != p:
        print(f"  Filtered {p - n_poly} monomorphic-in-train SNPs (kept {n_poly}/{p})",
              flush=True)
    dosage_train = dosage_train[poly_mask]
    dosage_test  = dosage_test[poly_mask]
    p = n_poly
    X_train, X_test = standardise(dosage_train, dosage_test)   # shape (p, n)

    # Sargent's β̂ — only meaningful if .beta.tsv exists (Stage 1 was run)
    has_beta = os.path.exists(out + ".beta.tsv")
    beta = load_beta(out + ".beta.tsv") if has_beta else np.zeros(p)
    if len(beta) != p:
        sys.exit(f"ERROR: β has length {len(beta)} but expected {p} polymorphic SNPs")

    # Sargent's ŵ² — only meaningful if .w2.tsv exists (Stage 2 was run)
    has_ws2 = os.path.exists(out + ".w2.tsv")
    ws2 = load_ws2(out + ".w2.tsv") if has_ws2 else None

    # ARG: load full and train, map trees
    print("Loading tree sequences ...", flush=True)
    ts_full  = tskit.load(case + ".full.trees")
    ts_train_only = tskit.load(case + ".trees")

    # Sample nodes for train and test in the full ARG
    train_nodes = []
    for i in range(n_train):
        train_nodes.extend(ts_full.individual(i).nodes)
    test_nodes = []
    for i in range(n_train, n_total):
        test_nodes.extend(ts_full.individual(i).nodes)

    # Map ws² (indexed against train.trees) to full-ARG tree indices, since
    # all our kinship computation is done in the full-ARG framework.
    if ws2 is not None and len(ws2) == ts_train_only.num_trees:
        idx_map = map_train_to_full(ts_train_only, ts_full)
        ws2_full = ws2[idx_map]
    elif ws2 is not None and len(ws2) == ts_full.num_trees:
        ws2_full = ws2
    else:
        ws2_full = None

    # Streaming kinship computation (avoids storing per-tree (n×n) matrices):
    #   K_uniform_train: Σ_s G_s / S_full       — for ARG-BLUP-uniform
    #   K_uniform_cross
    #   K_sargent_train: Σ_s ws²[s] · G_s        — for full SARGENT
    #   K_sargent_cross
    S_full = ts_full.num_trees
    ones = np.ones(S_full) / S_full

    print(f"Computing combined kinship over {S_full} full-ARG trees ...", flush=True)
    K_uni_train_hap  = streaming_combined_kinship(ts_full, train_nodes, ones)
    K_uni_cross_hap  = streaming_combined_kinship(ts_full, test_nodes, ones, train_nodes)
    K_uni_train = hap_to_indiv_2d(K_uni_train_hap, n_train)
    K_uni_cross = hap_to_indiv_2d(K_uni_cross_hap, n_test, n_train)

    if ws2_full is not None:
        K_sg_train_hap = streaming_combined_kinship(ts_full, train_nodes, ws2_full)
        K_sg_cross_hap = streaming_combined_kinship(ts_full, test_nodes, ws2_full, train_nodes)
        K_sg_train = hap_to_indiv_2d(K_sg_train_hap, n_train)
        K_sg_cross = hap_to_indiv_2d(K_sg_cross_hap, n_test, n_train)
    else:
        K_sg_train = K_sg_cross = None

    # -----------------------------------------------------------------------
    # Predictors
    # -----------------------------------------------------------------------
    results = {
        "case": case,
        "n_train": n_train,
        "n_test":  n_test,
        "n_trees_full":  ts_full.num_trees,
        "n_trees_train": ts_train_only.num_trees,
    }

    y_true_test = y_test
    y_true_components = None
    tcfile = case + ".test.true_components.tsv"
    if os.path.exists(tcfile):
        tc = np.loadtxt(tcfile, skiprows=1)
        if tc.ndim == 1:
            tc = tc[np.newaxis, :]
        y_true_components = {"fixed": tc[:, 1], "random": tc[:, 2], "total": tc[:, 3]}

    # (a) BayesR-only
    print("Predicting: BayesR-only ...", flush=True)
    y_hat_b = predict_bayesR(X_test, beta, y_train.mean())
    results["bayesR"] = metrics(y_hat_b, y_true_test)

    # (b) GBLUP — single genome-wide GRM, fitted h² on training
    print("Predicting: GBLUP ...", flush=True)
    K_g_train = (X_train.T @ X_train) / p          # (n_train, n_train)
    K_g_cross = (X_test.T  @ X_train) / p          # (n_test, n_train)
    K_g_train_c, K_g_cross_c = center_kinship(K_g_train, K_g_cross)
    s2, su2, alpha, mean_y = fit_blup(y_train, K_g_train_c)
    y_hat_g = predict_blup(K_g_cross_c, alpha, mean_y)
    results["gblup"] = metrics(y_hat_g, y_true_test)
    results["gblup_h2"] = su2 / (su2 + s2)

    # (c) ARG-BLUP uniform — average of per-tree kinship
    print("Predicting: ARG-BLUP-uniform ...", flush=True)
    K_uni_train_c, K_uni_cross_c = center_kinship(K_uni_train, K_uni_cross)
    s2, su2, alpha, mean_y = fit_blup(y_train, K_uni_train_c)
    y_hat_a = predict_blup(K_uni_cross_c, alpha, mean_y)
    results["arg_blup_uniform"] = metrics(y_hat_a, y_true_test)
    results["arg_blup_uniform_h2"] = su2 / (su2 + s2)

    # (d) Full SARGENT
    y_hat_s = None
    if ws2_full is not None and has_beta and K_sg_train is not None:
        print("Predicting: full SARGENT ...", flush=True)
        y_hat_s, sigma2_used = predict_sargent(
            X_test, beta, K_sg_train, K_sg_cross,
            y_train=y_train, X_train=X_train,
        )
        results["sargent"] = metrics(y_hat_s, y_true_test)
        results["sargent_sigma2"] = sigma2_used
    elif has_beta and ws2_full is None:
        # Stage 1 only — equivalent to BayesR
        results["sargent"] = results["bayesR"]

    # Decomposition metric: how well does SARGENT predict the random component alone?
    if y_true_components is not None and y_hat_s is not None:
        u_hat_s = y_hat_s - X_test.T @ beta
        results["sargent_u_r"] = float(
            np.corrcoef(u_hat_s, y_true_components["random"])[0, 1]
        )

    # -----------------------------------------------------------------------
    # Write
    # -----------------------------------------------------------------------
    out_path = args.out or (out + ".predict.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {out_path}", flush=True)

    print("\n=== Out-of-sample prediction ===")
    for k in ["bayesR", "gblup", "arg_blup_uniform", "sargent"]:
        if k in results:
            m = results[k]
            print(f"  {k:20s}  R²={m['r2']:+.3f}  r={m['r']:+.3f}  MSE={m['mse']:.3f}")


if __name__ == "__main__":
    main()
