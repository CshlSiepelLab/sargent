#!/usr/bin/env python3
"""
arch_metrics.py — Block 2.5.1: genetic-architecture summary metrics from
per-tree ŵ²_s output.

Three metrics, all computable from existing SARGENT outputs (no refitting):

  (1) Tree-level polygenicity index `M_e,ARG`
      M_e,ARG = (Σ_s ŵ²_s)² / Σ_s (ŵ²_s)²
      Analog of Wray's effective number of contributing variants, at the
      tree (≈ recombination block) level.  Sparse architecture → small
      `M_e`; diffuse → large.

  (2) Variance-capture curve
      Fraction of total ARG-variance captured by the top-k trees, plotted
      against k.  Shape characterises the architecture at a glance.

  (3) Per-branch effect-size distribution (untyped piece)
      Per branch b in tree s: variance contribution ∝ ŵ²_s · bl(b).
      Aggregated across the genome → non-parametric distribution of
      branch-level effect sizes.

Usage:
    # Single case
    python arch_metrics.py --case-prefix DIR/case --sargent-prefix DIR/est \
        --out DIR/arch_metrics.json

    # Batch (one row per case)
    python arch_metrics.py --batch ROOT_DIR --pattern "*/case" \
        --out ROOT_DIR/arch_metrics_batch.json
"""

import argparse
import glob
import json
import os
import sys

import numpy as np
import tskit

# Reuse existing helpers (load_ws2 reads a .w2.tsv file).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from predict import load_ws2


# ---------------------------------------------------------------------------
# Core metrics
# ---------------------------------------------------------------------------

def me_arg(w2):
    """Effective number of contributing trees (analog of Wray's M_e).

    M_e,ARG = (Σ_s w²_s)² / Σ_s (w²_s)²
            = 1 / Σ_s p_s²,   p_s = w²_s / Σ_t w²_t

    Returns NaN if Σ w² = 0.
    """
    w2 = np.asarray(w2, dtype=float)
    s1 = float(w2.sum())
    if s1 == 0.0:
        return float("nan")
    s2 = float((w2 * w2).sum())
    return (s1 * s1) / s2 if s2 > 0 else float("nan")


def variance_capture_curve(w2, sort_descending=True):
    """Cumulative-fraction-of-variance vs top-k trees.

    Returns:
        k          : array of length S, [1, 2, ..., S]
        captured   : array of length S, cumulative variance fraction in [0, 1]

    The captured-fraction at k = S is 1 by construction.
    """
    w2 = np.asarray(w2, dtype=float)
    if sort_descending:
        w2 = np.sort(w2)[::-1]
    cum = np.cumsum(w2)
    total = float(cum[-1]) if cum[-1] > 0 else 1.0
    return np.arange(1, len(w2) + 1), cum / total


def variance_capture_summary(w2, fractions=(0.1, 0.25, 0.5, 0.75, 0.9, 0.95)):
    """For each target fraction f, return the *smallest* k such that the
    top-k trees capture at least fraction f of the total variance.

    Equivalent shape summary: "k_at_50pct" is the median-variance count.
    """
    k_arr, capt = variance_capture_curve(w2)
    out = {}
    for f in fractions:
        idx = int(np.searchsorted(capt, f, side="left"))
        # idx is 0-based position in capt; the k value is k_arr[idx].
        if idx >= len(k_arr):
            out[f"k_at_{int(f*100)}pct"] = int(k_arr[-1])
        else:
            out[f"k_at_{int(f*100)}pct"] = int(k_arr[idx])
    out["S_total"] = int(len(w2))
    return out


def gini(w2):
    """Gini coefficient of ŵ²_s.  0 = perfectly uniform; 1 = one tree
    carries everything.  Complement of M_e,ARG / S — both summarise
    concentration but in different units."""
    w2 = np.asarray(w2, dtype=float)
    if w2.size == 0:
        return float("nan")
    if (w2 < 0).any():
        w2 = np.clip(w2, 0, None)
    s = w2.sum()
    if s == 0:
        return float("nan")
    sorted_w = np.sort(w2)
    n = len(sorted_w)
    cum_w = np.cumsum(sorted_w)
    return float((n + 1 - 2 * (cum_w.sum() / cum_w[-1])) / n)


def per_branch_effect_distribution(w2, ts, max_n=None):
    """Per-branch variance contribution: w²_s · bl(b) for branch b in tree s.

    Returns a 1-D array of variance values, one entry per (tree, branch).
    Skips trees with w² = 0 for efficiency.

    Args:
        w2: array of length ts.num_trees.
        ts: tskit TreeSequence (typically the full ARG).
        max_n: if set, downsample the returned array to at most max_n entries
               (for memory in plotting).

    Returns:
        np.ndarray of variance contributions.
    """
    w2 = np.asarray(w2, dtype=float)
    if len(w2) != ts.num_trees:
        raise ValueError(f"w2 length {len(w2)} != ts.num_trees {ts.num_trees}")
    contribs = []
    for s, tree in enumerate(ts.trees()):
        w_s = float(w2[s])
        if w_s == 0.0:
            continue
        for nd in tree.nodes():
            if tree.parent(nd) == tskit.NULL:
                continue
            bl = tree.branch_length(nd)
            if bl > 0:
                contribs.append(w_s * bl)
    arr = np.array(contribs, dtype=float)
    if max_n is not None and len(arr) > max_n:
        rng = np.random.default_rng(0)
        arr = rng.choice(arr, size=max_n, replace=False)
    return arr


def per_branch_summary(branches):
    """Summary stats for the per-branch variance distribution."""
    if len(branches) == 0:
        return {"n_branches": 0}
    return {
        "n_branches": int(len(branches)),
        "mean":   float(np.mean(branches)),
        "median": float(np.median(branches)),
        "p90":    float(np.quantile(branches, 0.90)),
        "p99":    float(np.quantile(branches, 0.99)),
        "max":    float(np.max(branches)),
        "gini_branches": gini(branches),
    }


def per_tree_mean_branch_age(ts):
    """For each tree, return the mean of (parent_time of branch) over all
    internal branches, weighted by branch length.

    This is the "age" of the tree in the same units as variant ages
    (generations before present, assuming the underlying simulation set up
    times that way). Used by Block 2.5.3 joint distributions: variance-
    capture curves stratified by tree-age bin.

    Returns:
        ages: np.ndarray of length ts.num_trees, in same time units as ts.
    """
    ages = np.zeros(ts.num_trees)
    for s, tree in enumerate(ts.trees()):
        total_bl = 0.0
        total_weighted_time = 0.0
        for nd in tree.nodes():
            parent = tree.parent(nd)
            if parent == tskit.NULL:
                continue
            bl = tree.branch_length(nd)
            if bl <= 0:
                continue
            total_bl += bl
            total_weighted_time += bl * tree.time(parent)
        ages[s] = total_weighted_time / total_bl if total_bl > 0 else float("nan")
    return ages


def capture_curve_by_age_bin(w2, ts, age_quantiles=(0.33, 0.67)):
    """Variance-capture curves stratified by tree-age quantile bins.

    Returns:
        bins: list of (lo, hi) age ranges
        curves: dict mapping (lo, hi) → (k_arr_normalised, captured_fraction)
                Each curve uses log-spaced sampling of k like the main
                variance_capture_curve.
        per_bin_sum: dict (lo, hi) → total Σ w² in bin (for normalisation
                     across bins later)
    """
    if len(w2) != ts.num_trees:
        return None
    ages = per_tree_mean_branch_age(ts)
    qs = sorted(set([0.0] + list(age_quantiles) + [1.0]))
    edges = np.quantile(ages, qs)
    curves = {}
    sums = {}
    bins = []
    for c in range(len(edges) - 1):
        lo, hi = float(edges[c]), float(edges[c + 1])
        if c == len(edges) - 2:
            mask = (ages >= lo) & (ages <= hi)
        else:
            mask = (ages >= lo) & (ages < hi)
        if not mask.any():
            continue
        bin_w2 = w2[mask]
        k_arr, capt = variance_capture_curve(bin_w2)
        # downsample k for storage
        idx = _logspace_sample(len(k_arr))
        curves[(lo, hi)] = {
            "k":      [int(k_arr[i])  for i in idx],
            "f":      [float(capt[i]) for i in idx],
            "S_bin":  int(len(bin_w2)),
            "sum_w2": float(bin_w2.sum()),
            "me_arg": me_arg(bin_w2),
        }
        sums[(lo, hi)] = float(bin_w2.sum())
        bins.append((lo, hi))
    return {"bins": bins, "curves": curves, "sums": sums}


# ---------------------------------------------------------------------------
# I/O for a single case
# ---------------------------------------------------------------------------

def find_w2_file(sargent_prefix):
    """Locate the .w2.tsv file (sargent's per-tree weight estimates)."""
    for suffix in [".w2.tsv"]:
        p = sargent_prefix + suffix
        if os.path.exists(p):
            return p
    return None


def find_true_w2_file(case_prefix):
    """Locate the .full.true_ws2.tsv (preferred, full-ARG indexed) or fall
    back to .true_ws2.tsv (train-tree indexed when only one ARG was saved)."""
    for suffix in [".full.true_ws2.tsv", ".true_ws2.tsv"]:
        p = case_prefix + suffix
        if os.path.exists(p):
            return p
    return None


def load_ts_for_branches(case_prefix):
    """Load the tree sequence to use for per-branch metrics.

    Prefer the FULL ARG (case.full.trees) since the per-branch distribution
    should be computed on all trees, not just the simplified-train subset.
    Falls back to case.trees if the full ARG isn't available.
    """
    for suffix in [".full.trees", ".trees"]:
        p = case_prefix + suffix
        if os.path.exists(p):
            return tskit.load(p), suffix
    return None, None


def compute_case_metrics(case_prefix, sargent_prefix, want_branches=True,
                         want_age_bins=False):
    """Compute the three metrics for a single case.

    Returns a dict suitable for serialisation as JSON.
    """
    result = {
        "case_prefix":   case_prefix,
        "sargent_prefix": sargent_prefix,
    }

    # Parameters of the case (for downstream stratification by architecture)
    pjf = case_prefix + ".params.json"
    params = json.load(open(pjf)) if os.path.exists(pjf) else {}
    for k in ("ws2_architecture", "ws2_block_frac", "ws2_block_contrast",
              "ws2_n_blocks", "demography", "h2_arg_realised",
              "h2_snp_realised", "n_trees", "n_trees_full", "seed"):
        if k in params:
            result[k] = params[k]

    # ---- Estimated metrics from sargent ŵ² -----------------------------------
    w2_file = find_w2_file(sargent_prefix)
    if w2_file is None:
        result["error"] = f"sargent .w2.tsv not found at {sargent_prefix}.w2.tsv"
        return result

    w2_est = load_ws2(w2_file)
    result["w2_file"] = w2_file
    result["S_est"]   = int(len(w2_est))
    result["est"] = {
        "me_arg":           me_arg(w2_est),
        "gini":             gini(w2_est),
        "capture_quantiles": variance_capture_summary(w2_est),
        "sum_w2":           float(w2_est.sum()),
    }
    # Save the curve sparsely (log-spaced k) so JSON doesn't blow up.
    k_arr, capt = variance_capture_curve(w2_est)
    sample_k = _logspace_sample(len(k_arr))
    result["est"]["capture_curve"] = {
        "k":  [int(k_arr[i])  for i in sample_k],
        "f":  [float(capt[i]) for i in sample_k],
    }

    # ---- Truth metrics (where available) --------------------------------------
    true_file = find_true_w2_file(case_prefix)
    if true_file is not None:
        w2_true = load_ws2(true_file)
        result["truth_file"] = true_file
        result["S_truth"]    = int(len(w2_true))
        result["truth"] = {
            "me_arg":           me_arg(w2_true),
            "gini":             gini(w2_true),
            "capture_quantiles": variance_capture_summary(w2_true),
            "sum_w2":           float(w2_true.sum()),
        }
        k_arr_t, capt_t = variance_capture_curve(w2_true)
        sample_k_t = _logspace_sample(len(k_arr_t))
        result["truth"]["capture_curve"] = {
            "k": [int(k_arr_t[i])  for i in sample_k_t],
            "f": [float(capt_t[i]) for i in sample_k_t],
        }

    # ---- Per-branch effect-size distribution ---------------------------------
    if want_branches or want_age_bins:
        ts, which = load_ts_for_branches(case_prefix)
    else:
        ts, which = None, None

    if want_age_bins and ts is not None:
        # Use estimated w² if it matches; else truth.
        for src_name, src_w2, src_S in (("est",   w2_est, result["S_est"]),
                                        ("truth", (w2_true if 'w2_true' in locals() else None),
                                                  result.get("S_truth", 0))):
            if src_w2 is None or src_S != ts.num_trees:
                continue
            ab = capture_curve_by_age_bin(src_w2, ts)
            if ab is not None:
                result.setdefault(src_name, {})
                result[src_name]["age_bin_curves"] = {
                    f"[{lo:.0f},{hi:.0f}]": c for (lo, hi), c in ab["curves"].items()
                }
                break

    if want_branches and ts is not None:
            # Choose which w2 vector to use: prefer truth if available and S
            # matches; else use estimated (and tolerate length mismatch).
            for src_name, src_w2, src_S in (("truth", result.get("truth") and w2_true, result.get("S_truth", 0)),
                                            ("est",   w2_est,                              result["S_est"])):
                if src_w2 is None:
                    continue
                if src_S != ts.num_trees:
                    # Length mismatch (e.g., train-only w2 vs full ARG):
                    # only attempt branches when they match.
                    continue
                branches = per_branch_effect_distribution(src_w2, ts, max_n=200_000)
                result.setdefault(src_name, {})
                result[src_name]["branches"] = per_branch_summary(branches)
                # Also save a downsampled raw branch array for plotting
                if len(branches) > 5_000:
                    rng = np.random.default_rng(0)
                    sampled = rng.choice(branches, size=5_000, replace=False)
                else:
                    sampled = branches
                result[src_name]["branches"]["sample"] = sampled.tolist()
                break   # one source is enough

            result["ts_source_for_branches"] = which

    return result


def _logspace_sample(n, n_points=64):
    """Return roughly log-spaced indices into [0, n-1] (no duplicates)."""
    if n <= n_points:
        return list(range(n))
    raw = np.unique(np.logspace(0, np.log10(n), n_points, dtype=int)) - 1
    raw = np.clip(raw, 0, n - 1)
    if raw[-1] != n - 1:
        raw = np.append(raw, n - 1)
    return list(np.unique(raw))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case-prefix",    default=None)
    ap.add_argument("--sargent-prefix", default=None)
    ap.add_argument("--batch",          default=None,
                    help="Root directory; will iterate matching subdirs.")
    ap.add_argument("--pattern",        default="*/case.params.json",
                    help="Glob (relative to --batch) for case .params.json "
                         "files; the case prefix is the matched path minus "
                         "'.params.json', and the sargent prefix is the "
                         "same directory's 'est'. Default: '*/case.params.json'.")
    ap.add_argument("--out",            default=None,
                    help="Output JSON path. For batch mode, a single JSON "
                         "with one entry per case is written.")
    ap.add_argument("--no-branches",    action="store_true",
                    help="Skip the per-branch effect-size distribution "
                         "(faster; useful for batch runs).")
    ap.add_argument("--age-bins",       action="store_true",
                    help="Block 2.5.3: also compute variance-capture curves "
                         "stratified by per-tree mean-branch-age tercile "
                         "(slower; requires the tree sequence).")
    args = ap.parse_args()

    want_branches = not args.no_branches
    want_age_bins = args.age_bins

    if args.batch:
        rows = []
        matches = sorted(glob.glob(os.path.join(args.batch, args.pattern)))
        for match in matches:
            # Derive case prefix from the matched .params.json path.
            if match.endswith(".params.json"):
                case = match[:-len(".params.json")]
            else:
                case = match
            if not os.path.exists(case + ".params.json"):
                continue
            dir_ = os.path.dirname(case)
            est_prefix = os.path.join(dir_, "est")
            print(f"[arch_metrics] {os.path.relpath(case, args.batch)}",
                  flush=True)
            rows.append(compute_case_metrics(case, est_prefix,
                                             want_branches=want_branches,
                                             want_age_bins=want_age_bins))
        out_path = args.out or os.path.join(args.batch, "arch_metrics_batch.json")
        with open(out_path, "w") as f:
            json.dump(rows, f, indent=2)
        print(f"Wrote {out_path} ({len(rows)} cases)")
    else:
        if not (args.case_prefix and args.sargent_prefix):
            sys.exit("ERROR: --case-prefix and --sargent-prefix required "
                     "(or use --batch).")
        result = compute_case_metrics(args.case_prefix, args.sargent_prefix,
                                      want_branches=want_branches,
                                      want_age_bins=want_age_bins)
        out_path = args.out or (args.sargent_prefix + ".arch_metrics.json")
        with open(out_path, "w") as f:
            json.dump(result, f, indent=2)
        print(f"Wrote {out_path}")
        if "est" in result:
            print(f"  M_e,ARG (est)  = {result['est']['me_arg']:.2f}")
            print(f"  Gini    (est)  = {result['est']['gini']:.3f}")
            print(f"  k_at_50pct      = {result['est']['capture_quantiles']['k_at_50pct']}")
        if "truth" in result:
            print(f"  M_e,ARG (truth)= {result['truth']['me_arg']:.2f}")
            print(f"  Gini    (truth)= {result['truth']['gini']:.3f}")


if __name__ == "__main__":
    main()
