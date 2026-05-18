"""
metrics.py — shared evaluation metrics for SARGENT test suite.
"""

import json
import numpy as np


def pearson_r(a, b):
    """Pearson correlation; returns NaN if either array is constant."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.std() == 0 or b.std() == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def load_params(prefix):
    with open(prefix + ".params.json") as f:
        return json.load(f)


def load_true_ws2(prefix):
    """Returns array of true ws² values indexed by tree_idx (0-based)."""
    rows = np.loadtxt(prefix + ".true_ws2.tsv", skiprows=1, usecols=(0, 1))
    if rows.ndim == 1:
        rows = rows[np.newaxis, :]
    idx = rows[:, 0].astype(int)
    ws2 = rows[:, 1]
    n_trees = int(idx.max()) + 1
    out = np.zeros(n_trees)
    out[idx] = ws2
    return out


def load_true_beta(prefix, n_variants):
    """Returns array of true β values (length n_variants, 0 for non-causal)."""
    beta = np.zeros(n_variants)
    rows = np.loadtxt(prefix + ".true_beta.tsv", skiprows=1, usecols=(0, 2))
    if rows.size == 0:
        return beta
    if rows.ndim == 1:
        rows = rows[np.newaxis, :]
    idx  = rows[:, 0].astype(int)
    vals = rows[:, 1]
    beta[idx] = vals
    return beta


def load_inferred_ws2(out_prefix):
    """Load E_ws2 column from sargent .w2.tsv output."""
    data = np.loadtxt(out_prefix + ".w2.tsv", skiprows=1)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return data[:, -1]      # last column is E_ws2


def load_inferred_rho(out_prefix):
    """Load rho_active column (last rho before E_ws2) from sargent .w2.tsv output.

    Returns posterior activation probability per tree (length n_trees).
    """
    with open(out_prefix + ".w2.tsv") as f:
        header = f.readline().strip().split("\t")
    # Find the last rho_* column index
    rho_cols = [i for i, h in enumerate(header) if h.startswith("rho_")]
    if not rho_cols:
        return None
    last_rho = rho_cols[-1]   # rho_M = posterior probability of the active slab
    data = np.loadtxt(out_prefix + ".w2.tsv", skiprows=1)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return data[:, last_rho]


def roc_auc(scores, labels):
    """Compute ROC AUC of scores against binary labels (0/1).
    NaN if labels are all-same."""
    scores = np.asarray(scores, float)
    labels = np.asarray(labels, int)
    pos = labels == 1
    neg = labels == 0
    n_pos, n_neg = pos.sum(), neg.sum()
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    # Mann-Whitney U formulation
    order = np.argsort(scores)
    ranks = np.empty_like(order, dtype=float)
    # Handle ties by averaging ranks
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


def load_inferred_beta(out_prefix):
    """Load E_beta column from sargent .beta.tsv output."""
    data = np.loadtxt(out_prefix + ".beta.tsv", skiprows=1)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return data[:, -1]      # last column is E_beta


def ws2_metrics(true_ws2, inferred_ws2, inferred_rho=None):
    """
    Returns dict with:
      ws2_corr_active   — Pearson r of log(true) vs log(inferred) among truly active trees
      ws2_corr_all      — Pearson r across all trees (treating null as min)
      auc_rho           — ROC AUC of rho_active for classifying truly-active trees
      auc_ws2           — ROC AUC of inferred ws² for classifying truly-active trees
      n_active_true     — # trees with true ws² > 0
      n_called_active   — # trees with rho_active >= 0.5 (if rho available)
    """
    n = len(true_ws2)
    assert len(inferred_ws2) == n, f"ws2 length mismatch: {len(inferred_ws2)} vs {n}"

    active_true = (true_ws2 > 0).astype(int)

    # log-correlation among truly-active trees
    if active_true.sum() >= 2:
        mask = active_true == 1
        log_true = np.log(true_ws2[mask])
        log_inf  = np.log(np.maximum(inferred_ws2[mask], 1e-300))
        corr_active = pearson_r(log_true, log_inf)
    else:
        corr_active = float("nan")

    floor_for_null = max(inferred_ws2[inferred_ws2 > 0].min() if (inferred_ws2 > 0).any() else 1e-30,
                         1e-30)
    log_inf_all  = np.log(np.where(inferred_ws2 > 0, inferred_ws2, floor_for_null))
    log_true_all = np.log(np.where(true_ws2 > 0, true_ws2,
                                   true_ws2[true_ws2 > 0].min() / 10
                                   if (true_ws2 > 0).any() else 1.0))
    corr_all = pearson_r(log_inf_all, log_true_all)

    out = {
        "ws2_corr_active":   corr_active,
        "ws2_corr_all":      corr_all,
        "n_trees":           n,
        "n_active_true":     int(active_true.sum()),
        "frac_active_true":  float(active_true.mean()),
        "auc_ws2":           roc_auc(inferred_ws2, active_true),
    }

    if inferred_rho is not None and len(inferred_rho) == n:
        out["auc_rho"] = roc_auc(inferred_rho, active_true)
        called = (inferred_rho >= 0.5).astype(int)
        tp = int(((called == 1) & (active_true == 1)).sum())
        fp = int(((called == 1) & (active_true == 0)).sum())
        fn = int(((called == 0) & (active_true == 1)).sum())
        out["n_called_active"] = int(called.sum())
        out["precision_rho"]   = tp / max(tp + fp, 1) if (tp + fp) > 0 else float("nan")
        out["recall_rho"]      = tp / max(tp + fn, 1) if active_true.sum() > 0 else float("nan")

    return out


def beta_metrics(true_beta, inferred_beta):
    """
    Returns dict with:
      corr_causal   — Pearson r among causal SNPs (true_beta != 0)
      corr_all      — Pearson r across all SNPs
      mse_causal    — MSE among causal SNPs
    """
    assert len(true_beta) == len(inferred_beta), "beta length mismatch"
    causal = true_beta != 0

    if causal.sum() >= 2:
        corr_causal = pearson_r(true_beta[causal], inferred_beta[causal])
        mse_causal  = float(np.mean((true_beta[causal] - inferred_beta[causal]) ** 2))
    else:
        corr_causal = float("nan")
        mse_causal  = float("nan")

    return {
        "beta_corr_causal": corr_causal,
        "beta_corr_all":    pearson_r(true_beta, inferred_beta),
        "beta_mse_causal":  mse_causal,
        "n_causal":         int(causal.sum()),
    }


def variance_metrics(sigma2_inferred, sigma2_true, h2_fixed_inferred, h2_fixed_true,
                     h2_rand_inferred, h2_rand_true):
    """Returns dict with relative errors for sigma2, h2_fixed, h2_random."""
    def rel_err(inferred, true):
        if true == 0:
            return float("nan") if inferred == 0 else float("inf")
        return (inferred - true) / true

    return {
        "sigma2_inferred":  sigma2_inferred,
        "sigma2_true":      sigma2_true,
        "sigma2_relerr":    rel_err(sigma2_inferred, sigma2_true),
        "h2_fixed_inferred":  h2_fixed_inferred,
        "h2_fixed_true":      h2_fixed_true,
        "h2_fixed_relerr":    rel_err(h2_fixed_inferred, h2_fixed_true),
        "h2_rand_inferred":   h2_rand_inferred,
        "h2_rand_true":       h2_rand_true,
        "h2_rand_relerr":     rel_err(h2_rand_inferred, h2_rand_true),
    }
