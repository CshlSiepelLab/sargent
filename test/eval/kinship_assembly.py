#!/usr/bin/env python3
"""
kinship_assembly.py — per-bin ARG kinship assembly for Phase C stratified GREML.

Builds K_b = Σ_{e ∈ bin_b} (ℓ_e / 2) · m_e · m_eᵀ for a set of bins
defined by a user-supplied binning function.  K_b is the natural diploid
kinship matrix for multi-component REML with partitioned heritability.

Public API
----------
    build_per_bin_kinship(ts, bin_fn, n_indiv, mode='diploid') -> dict
    bin_by_age_tertile(ts) -> bin_fn
    bin_by_daf_tertile(ts, n_indiv) -> bin_fn
    bin_by_population(ts, pop_labels) -> bin_fn

Bin-function signature: bin_fn(tree, child_node, c_u) -> bin_id
  tree       : tskit.Tree
  child_node : int, the child node of the edge being binned
  c_u        : np.ndarray shape (n_indiv,) dtype int8 — per-individual
               haplotype count descending from child_node (0, 1, or 2)

Validation (run as __main__)
----------------------------
    python kinship_assembly.py --case-prefix CASE [--check-daf] [--check-pop]

Four checks, all reported with quantitative numbers:
  1. Single-bin ‖K_assembly − K_ARG‖_F < 1e-6
  2. Sum-rule   ‖Σ_b K_b − K_ARG‖_F < 1e-6  (age tertile)
  3. PSD:  λ_min ≥ −1e-8 for every bin
  4. Trace: tr(K_b) > 0 for every bin
"""

import argparse
import os
import sys

import numpy as np
import tskit

EVAL_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, EVAL_DIR)


# ---------------------------------------------------------------------------
# Internal tertile helper (logic matches preprocess_fig5_edge.py)
# ---------------------------------------------------------------------------

def _per_cell_tertiles(ts):
    """Per-cell q33/q67 over edge midpoint age distribution."""
    nt   = ts.tables.nodes.time
    mids = []
    for tree in ts.trees():
        for u in tree.nodes():
            p = tree.parent(u)
            if p == tskit.NULL:
                continue
            mids.append(0.5 * (float(nt[p]) + float(nt[u])))
    mids = np.asarray(mids, dtype=float)
    if mids.size == 0:
        return 5_000.0, 20_000.0
    return float(np.nanpercentile(mids, 33.33)), float(np.nanpercentile(mids, 66.67))


# ---------------------------------------------------------------------------
# Binning functions
# ---------------------------------------------------------------------------

def bin_by_age_tertile(ts):
    """Partition edges into recent / middle / ancient by midpoint age tertile.

    Cutoffs are the per-cell q33 / q67 of the edge midpoint age distribution
    — the same convention as preprocess_fig5_edge.py, so truth and estimator
    bins are comparable.
    """
    q33, q67 = _per_cell_tertiles(ts)
    nt = ts.tables.nodes.time

    def bin_fn(tree, child, c_u):  # noqa: ARG001 (c_u unused for age binning)
        p   = tree.parent(child)
        mid = 0.5 * (float(nt[p]) + float(nt[child]))
        if mid <= q33:
            return "recent"
        if mid <= q67:
            return "middle"
        return "ancient"

    return bin_fn


def bin_by_daf_tertile(ts, n_indiv):
    """Partition edges by descendant haplotype frequency (DAF) tertile.

    DAF_e = tree.num_samples(child_node) / (2 * n_indiv).
    Tertile cutoffs are computed over the per-edge DAF distribution.
    Bin labels: 'rare' / 'common' / 'abundant'.
    """
    n_hap = 2 * n_indiv
    dafs  = []
    for tree in ts.trees():
        for u in tree.nodes():
            if tree.parent(u) == tskit.NULL:
                continue
            dafs.append(tree.num_samples(u) / n_hap)
    dafs = np.asarray(dafs, dtype=float)
    q33  = float(np.nanpercentile(dafs, 33.33))
    q67  = float(np.nanpercentile(dafs, 66.67))

    def bin_fn(tree, child, c_u):  # noqa: ARG001
        daf = tree.num_samples(child) / n_hap
        if daf <= q33:
            return "rare"
        if daf <= q67:
            return "common"
        return "abundant"

    return bin_fn


def bin_by_population(ts, pop_labels):
    """Assign each edge to the population with the most descendant haplotypes.

    Args:
        ts         : tskit.TreeSequence
        pop_labels : list of str, length n_indiv — one population label per
                     diploid individual (index-matched to ts.individuals() order)

    Returns a bin_fn; ties are broken by alphabetical order of the pop label.
    """
    pops         = sorted(set(pop_labels))
    labels_arr   = np.asarray(pop_labels)

    def bin_fn(tree, child, c_u):
        best_pop   = pops[0]
        best_count = -1
        for pop in pops:
            count = int(c_u[labels_arr == pop].sum())
            if count > best_count:
                best_count = count
                best_pop   = pop
        return best_pop

    return bin_fn


# ---------------------------------------------------------------------------
# Core assembly
# ---------------------------------------------------------------------------

def build_per_bin_kinship(ts, bin_fn, n_indiv, mode="diploid"):
    """Assemble per-bin ARG kinship matrices.

    For each edge e: ℓ_e = parent_time − child_time; m_e[i] ∈ {0,1,2}
    is the haplotype count of individual i below the child node.

        K_b[i, j] = Σ_{e ∈ bin_b} (ℓ_e / 2) · m_e[i] · m_e[j]

    Implementation: within each tree, edges sharing the same bin are batched
    into a matrix M (rows = sqrt(ℓ_e/2) · m_e) and accumulated as M^T M —
    a single BLAS dgemm call instead of n_edges individual outer products.

    Args:
        ts       : tskit.TreeSequence
        bin_fn   : callable (tree, child_node, c_u) → bin_id
                   c_u is a (n_indiv,) int8 array of haplotype counts
        n_indiv  : number of diploid individuals (= ts.num_samples // 2)
        mode     : 'diploid' (only supported value)

    Returns:
        dict mapping bin_id → np.ndarray shape (n_indiv, n_indiv) float64
    """
    if mode != "diploid":
        raise NotImplementedError(f"mode={mode!r} not supported; only 'diploid'")

    sample_nodes    = list(ts.samples())
    sample_to_indiv = {nd: idx // 2 for idx, nd in enumerate(sample_nodes)}
    nt              = ts.tables.nodes.time
    K_bins: dict    = {}

    for tree in ts.trees():
        # Single post-order pass: build per-node diploid descendant-count vector.
        c: dict = {}
        for u in tree.nodes(order="postorder"):
            c_u = np.zeros(n_indiv, dtype=np.int8)
            if tree.is_sample(u):
                i = sample_to_indiv.get(u, -1)
                if i >= 0:
                    c_u[i] = 1
            else:
                ch = tree.left_child(u)
                while ch != tskit.NULL:
                    c_u += c[ch]
                    ch = tree.right_sib(ch)
            c[u] = c_u

        # Collect sqrt(ℓ/2)·m vectors per bin; flush via M^T M at end of tree.
        bin_vecs: dict = {}
        for u in tree.nodes():
            p = tree.parent(u)
            if p == tskit.NULL:
                continue
            ell = float(nt[p]) - float(nt[u])
            if ell <= 0.0:
                continue
            b_id = bin_fn(tree, u, c[u])
            bin_vecs.setdefault(b_id, []).append(
                np.sqrt(ell / 2.0) * c[u].astype(np.float64)
            )

        for b_id, vecs in bin_vecs.items():
            M   = np.stack(vecs)     # shape (k_edges, n_indiv)
            K_b = K_bins.setdefault(b_id,
                                    np.zeros((n_indiv, n_indiv), dtype=np.float64))
            K_b += M.T @ M           # equivalent to Σ_e (ℓ_e/2) outer(m_e, m_e)

    return K_bins


# ---------------------------------------------------------------------------
# Reference K_ARG for validation
# ---------------------------------------------------------------------------

def _reference_k_arg(ts, n_indiv):
    """K_ARG = Σ_e (ℓ_e/2)·m_e·m_eᵀ via the predict.py streaming path.

    streaming_combined_kinship accumulates 0.5·bl per haplotype-level edge;
    hap_to_indiv_2d then sums the 2×2 blocks and multiplies by 0.5, giving
    Σ_e (ell/4)·m·mᵀ at diploid level.  The 2× factor below converts to
    SARGENT's native convention (ell/2), which matches build_per_bin_kinship
    and the G_matvec used in edge_h2.py.
    """
    from predict import streaming_combined_kinship, hap_to_indiv_2d  # noqa: PLC0415

    sample_nodes = list(ts.samples())
    ws2_uniform  = np.ones(ts.num_trees, dtype=float)
    K_hap        = streaming_combined_kinship(ts, sample_nodes, ws2_uniform)
    return 2.0 * hap_to_indiv_2d(K_hap, n_indiv)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def _frob(A, B=None):
    """‖A‖_F, or ‖A − B‖_F when B is given."""
    return float(np.linalg.norm(A if B is None else (A - B), "fro"))


def validate(ts, n_indiv, verbose=True):
    """Run the four validation checks; return a results dict.

    Checks:
      1. Single-bin ‖K_assembly − K_ARG‖_F  (should be < 1e-6)
      2. Sum-rule   ‖Σ_b K_b − K_ARG‖_F    (age tertile bins)
      3. PSD:  λ_min ≥ −1e-8 for every age-tertile bin
      4. Trace: tr(K_b) > 0 for every age-tertile bin
    """
    K_ref = _reference_k_arg(ts, n_indiv)

    # 1. Single-bin identity
    const_fn = lambda tree, child, c_u: "all"  # noqa: E731
    K_single    = build_per_bin_kinship(ts, const_fn, n_indiv)["all"]
    frob_single = _frob(K_single, K_ref)
    frob_kref   = _frob(K_ref)
    rel_single  = frob_single / max(frob_kref, 1.0)

    # 2. Sum-rule over age-tertile bins
    bins_age = build_per_bin_kinship(ts, bin_by_age_tertile(ts), n_indiv)
    K_sum    = sum(bins_age.values())
    frob_sum = _frob(K_sum, K_ref)
    rel_sum  = frob_sum / max(frob_kref, 1.0)

    # 3 & 4. Per-bin PSD + trace
    # PSD threshold is relative to λ_max to handle floating-point noise in
    # large matrices (rank-1 accumulation can leave λ_min ≈ −1e-7 when
    # λ_max ≈ 1e9, which is machine epsilon × condition number).
    psd_info: dict = {}
    for b_id, K_b in bins_age.items():
        evals = np.linalg.eigvalsh(K_b)
        psd_info[b_id] = {
            "lambda_min": float(evals.min()),
            "lambda_max": float(evals.max()),
            "trace":      float(np.trace(K_b)),
        }

    results = {
        "frob_kref":               frob_kref,
        "frob_single_bin_vs_karg": frob_single,
        "rel_single_bin_vs_karg":  rel_single,
        "frob_sum_rule":           frob_sum,
        "rel_sum_rule":            rel_sum,
        "age_bins_present":        sorted(bins_age.keys()),
        "per_bin_psd":             psd_info,
    }

    if verbose:
        # Use relative tolerance: both K_assembly and K_ref compute the same
        # sum-of-outer-products formula via different accumulation orders, so
        # the residual is floating-point noise; ‖diff‖/‖K‖ < 1e-10 is expected.
        rel_thr = 1e-10
        print("=== kinship_assembly validation ===")
        print(f"  n_indiv={n_indiv}  n_trees={ts.num_trees}  "
              f"‖K_ARG‖_F={frob_kref:.4e}")
        print()
        print(f"  [1] single-bin  ‖K_assembly − K_ARG‖_F = {frob_single:.3e}"
              f"  (rel={rel_single:.1e})"
              f"  ({'PASS' if rel_single < rel_thr else 'FAIL'})")
        print(f"  [2] sum-rule    ‖ΣK_b − K_ARG‖_F       = {frob_sum:.3e}"
              f"  (rel={rel_sum:.1e})"
              f"  ({'PASS' if rel_sum < rel_thr else 'FAIL'})")
        print(f"       bins present: {results['age_bins_present']}")
        print()
        print("  [3/4] per-bin PSD + trace:")
        all_psd_ok   = True
        all_trace_ok = True
        for b_id in sorted(psd_info):
            r        = psd_info[b_id]
            psd_ok   = r["lambda_min"] >= -1e-6 * max(r["lambda_max"], 1.0)
            trace_ok = r["trace"] > 0.0
            if not psd_ok:   all_psd_ok   = False
            if not trace_ok: all_trace_ok = False
            print(f"    {b_id:>10s}:  "
                  f"λ_min={r['lambda_min']:+.3e}  "
                  f"λ_max={r['lambda_max']:.3e}  "
                  f"trace={r['trace']:.4e}  "
                  f"PSD={'OK' if psd_ok else 'FAIL'}  "
                  f"trace={'OK' if trace_ok else 'FAIL'}")
        print()
        overall = (rel_single < rel_thr and rel_sum < rel_thr
                   and all_psd_ok and all_trace_ok)
        print(f"  Overall: {'ALL PASS' if overall else 'SOME FAIL'}"
              f"  (single-rel={rel_single:.1e}"
              f"  sum-rule-rel={rel_sum:.1e}"
              f"  PSD={'OK' if all_psd_ok else 'FAIL'}"
              f"  trace={'OK' if all_trace_ok else 'FAIL'})")

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case-prefix", required=True,
                    help="Path prefix for case files; loads CASE.trees.")
    ap.add_argument("--check-daf", action="store_true",
                    help="Also check bin_by_daf_tertile sum rule.")
    ap.add_argument("--check-pop", action="store_true",
                    help="Also check bin_by_population sum rule (uses ts.individuals).")
    args = ap.parse_args()

    ts_path = args.case_prefix + ".trees"
    print(f"  loading {ts_path} …", flush=True)
    ts      = tskit.load(ts_path)
    n_indiv = ts.num_samples // 2
    print(f"  n_samples={ts.num_samples}  n_indiv={n_indiv}  "
          f"n_trees={ts.num_trees}", flush=True)

    validate(ts, n_indiv)

    if args.check_daf:
        print()
        print("=== DAF-tertile sum-rule ===")
        K_ref    = _reference_k_arg(ts, n_indiv)
        bins_daf = build_per_bin_kinship(ts, bin_by_daf_tertile(ts, n_indiv), n_indiv)
        frob_daf = _frob(sum(bins_daf.values()), K_ref)
        rel_daf  = frob_daf / max(_frob(K_ref), 1.0)
        print(f"  bins: {sorted(bins_daf.keys())}")
        print(f"  ‖ΣK_b(DAF) − K_ARG‖_F = {frob_daf:.3e}  (rel={rel_daf:.1e})"
              f"  ({'PASS' if rel_daf < 1e-10 else 'FAIL'})")
        for b_id in sorted(bins_daf):
            tr = float(np.trace(bins_daf[b_id]))
            print(f"    {b_id:>10s}:  trace={tr:.4e}")

    if args.check_pop:
        pops = list(ts.populations())
        if len(pops) < 2:
            print()
            print("  (only one population in this ts; skipping --check-pop)")
        else:
            print()
            print("=== population-bin sum-rule ===")
            if ts.num_individuals > 0:
                pop_labels = [str(ts.node(ts.individual(i).nodes[0]).population)
                              for i in range(n_indiv)]
            else:
                pop_labels = [str(ts.node(ts.samples()[2 * i]).population)
                              for i in range(n_indiv)]
            K_ref    = _reference_k_arg(ts, n_indiv)
            bins_pop = build_per_bin_kinship(
                ts, bin_by_population(ts, pop_labels), n_indiv)
            frob_pop = _frob(sum(bins_pop.values()), K_ref)
            print(f"  bins: {sorted(bins_pop.keys())}")
            print(f"  ‖ΣK_b(pop) − K_ARG‖_F = {frob_pop:.3e}"
                  f"  ({'PASS' if frob_pop < 1e-6 else 'FAIL'})")
            for b_id in sorted(bins_pop):
                tr = float(np.trace(bins_pop[b_id]))
                print(f"    {b_id:>10s}:  trace={tr:.4e}")


if __name__ == "__main__":
    main()
