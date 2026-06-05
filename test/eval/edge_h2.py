#!/usr/bin/env python3
"""
edge_h2.py — per-edge posterior contribution to phenotypic variance.

Implements the calculation specified in test/eval/edge_h2_design.md.

Model (as implemented in src/lib/variational.c):
    V       =  σ² I + Σ_s  ŵ²_s · G_s
    u_s    ~~  N(0, ŵ²_s · G_s)                   (no extra σ² factor on per-tree var)
    per-edge prior  α_e ~ N(0, ŵ²_s · ℓ_e / 2)    (matvec decomposition, 1/2 from
                                                    diploid hap-splitting)

For each cell (one tree sequence + one sargent fit), compute:

    z   :=  V^{-1} (y - X β̂)               (BLUP residual)
    E[α_e | y]  =  ŵ²_s · ℓ_e · S_z[u(e)]
                                            (S_z is the SARGENT-matvec analog of S
                                             with z instead of `in`; the 1/2 prior
                                             factor is absorbed by the matvec's
                                             2·S_z descendant-sum identity)
    β_e         :=  Σ_{j: SNP j sits on e} β̂_j
                                            (cumulative typed-SNP effect on edge e)
    Δh²_e^ARG   =  (E[α_e | y])² · Σ_i (c_i^{u(e)})² / (n · Var(y))
    Δh²_e^SNP   =  β_e² · Σ_i (c_i^{u(e)})² / (n · Var(y))
    Δh²_e^total =  (E[α_e | y] + β_e)² · Σ_i (c_i^{u(e)})² / (n · Var(y))

Same descendant-count structure on both sides, so SNP and ARG are on the same basis.

This file lives next to predict.py and reuses its helpers (read_bed,
standardise, streaming_combined_kinship, load_ws2, load_beta).

CLI:

    python edge_h2.py --case-prefix CASE --sargent-prefix EST [--out PATH.tsv]
                      [--top-k-trees K]  [--top-k-edges K]
                      [--age-bin-edges 0,5000,20000,inf]
                      [--summary-json SUMMARY.json]

Outputs per-edge TSV plus optional aggregated summaries to stdout.
"""

import argparse
import csv
import json
import os
import sys

import numpy as np
import tskit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from predict import (
    read_bed,
    standardise,
    load_beta,
    load_ws2,
    streaming_combined_kinship,
    hap_to_indiv_2d,
)


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_phen(phen_path):
    """Read phenotype vector (last whitespace-separated column per row)."""
    vals = []
    with open(phen_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            vals.append(float(line.split()[-1]))
    return np.asarray(vals, dtype=float)


def load_summary(summary_path):
    return json.load(open(summary_path))


def compute_X_standardised(case_prefix, n_indiv, n_variants):
    """Load BED, standardise as sargent does (mean 0, var 1, train moments)."""
    dosage = read_bed(case_prefix, n_indiv, n_variants)        # (p, n)
    # use the cell's own moments (we have only one cohort here)
    dosage_std, _ = standardise(dosage, dosage)
    return dosage_std                                          # (p, n)


# ---------------------------------------------------------------------------
# V z = r solver  (direct n×n solve; appropriate for n ≲ 5000)
# ---------------------------------------------------------------------------

def compute_V_and_solve(ts, ws2, sigma2, residual, sample_nodes=None):
    """
    Build V = σ² I + Σ_s ŵ²_s · G_s   (diploid kinship)  and solve V z = residual.

    Uses predict.py::streaming_combined_kinship to accumulate the weighted ARG
    kinship at the haplotype level, then collapses to diploid level.

    Args:
        ts          : tskit tree sequence
        ws2         : per-tree variance weights (length S)
        sigma2      : sargent's σ²
        residual    : (n,) vector to solve against
        sample_nodes: list of haplotype sample node IDs (default: all samples)

    Returns:
        z           : V^{-1} residual           (length n)
        V           : the assembled n×n matrix  (returned for self-consistency checks)
    """
    if sample_nodes is None:
        sample_nodes = list(ts.samples())
    n_hap = len(sample_nodes)
    if n_hap % 2 != 0:
        raise ValueError(f"odd haplotype count {n_hap}: not diploid")
    n = n_hap // 2

    K_hap = streaming_combined_kinship(ts, sample_nodes, ws2)
    K     = hap_to_indiv_2d(K_hap, n)                          # (n, n)
    # streaming_combined_kinship + hap_to_indiv_2d returns half the matvec-scale
    # G_s (the python helper uses 0.5·bl at branches AND 0.5 at hap→indiv collapse
    # → factor 0.25 of shared_length; sargent's matvec gives 0.5).  Scale up to
    # match SARGENT's V = σ²I + Σ ŵ²_s · G_s_matvec parameterisation.
    G_matvec_combined = 2.0 * K
    V     = sigma2 * np.eye(n) + G_matvec_combined
    z     = np.linalg.solve(V, residual)
    return z, V, G_matvec_combined


# ---------------------------------------------------------------------------
# Per-tree single-pass descendant accumulators
# ---------------------------------------------------------------------------

def per_node_descendant_arrays(tree, sample_to_indiv, n_indiv):
    """
    For one tree, compute per-node c[u] (length n_indiv): the count of
    haplotypes of each individual descending from node u.  c[u][i] ∈ {0, 1, 2}.

    Single post-order pass.  O(n_nodes × n_indiv) memory; OK for n_indiv ≲ 5000.

    Returns:
        c           : dict {node_id: np.ndarray of shape (n_indiv,) int}
        S_one_arr   : dict {node_id: n_haps_in_subtree(u)} = sum_i c[u][i]
    """
    c       = {}
    S_one   = {}
    # Post-order: tree.nodes(order="postorder") iterates over nodes of this tree
    for u in tree.nodes(order="postorder"):
        c_u = np.zeros(n_indiv, dtype=np.int8)
        if tree.is_sample(u):
            # leaf == one haplotype
            i = sample_to_indiv.get(u, -1)
            if i >= 0:
                c_u[i] = 1
        else:
            # internal: sum from children
            child = tree.left_child(u)
            while child != tskit.NULL:
                c_u += c[child]
                child = tree.right_sib(child)
        c[u] = c_u
        S_one[u] = int(c_u.sum())
    return c, S_one


def per_node_S_z(c_per_node, z):
    """S_z[u] = 0.5 · sum_i c[u][i] · z[i] for each node u."""
    return {u: 0.5 * float(c.dot(z)) for u, c in c_per_node.items()}


def per_node_sum_csquared(c_per_node):
    """Σ_i (c[u][i])²  for each node u.  Used in ‖u_e^posterior‖²."""
    return {u: int((c.astype(int) ** 2).sum()) for u, c in c_per_node.items()}


# ---------------------------------------------------------------------------
# Per-edge SNP coefficient accumulator
# ---------------------------------------------------------------------------

def per_edge_beta_sum(ts, beta, snp_to_mutation=None):
    """
    For each branch (edge), accumulate β̂_j over the typed SNPs whose mutation
    sits on that branch.  Returns dict {(tree_idx, child_node): β_e}.

    By default we use a position-based lookup: SNP j has a coordinate from
    case.bim; the BED order is assumed to match the ts.sites() order so that
    site index j → mutation in tree containing position pos.  When the
    `snp_to_mutation` argument is provided (length p mapping into ts.mutations
    indices), we use it directly.

    For each mutation m: tree_idx = tree containing m.position; edge_child = m.node.
    """
    sums = {}

    if snp_to_mutation is not None:
        # explicit mapping
        for j, m_idx in enumerate(snp_to_mutation):
            if m_idx < 0:
                continue
            m = ts.mutation(m_idx)
            tree_idx = ts.at(m.position).index  # may be slow but correct
            key = (tree_idx, m.node)
            sums[key] = sums.get(key, 0.0) + float(beta[j])
        return sums

    # Default: assume SNPs in BED order match the sites in ts.sites() order
    # (which is what simulate.py does — it writes one BED record per site).
    if ts.num_sites != len(beta):
        # Not a strict requirement, but warn.
        print(f"  warning: ts.num_sites={ts.num_sites} != len(beta)={len(beta)}; "
              "SNP order assumption may be wrong", file=sys.stderr)

    for tree_idx, tree in enumerate(ts.trees()):
        for site in tree.sites():
            # site.id is the index into ts.sites()
            j = site.id
            if j >= len(beta):
                continue
            bj = float(beta[j])
            if bj == 0.0:
                continue
            for m in site.mutations:
                key = (tree_idx, m.node)
                sums[key] = sums.get(key, 0.0) + bj
    return sums


# ---------------------------------------------------------------------------
# Main per-edge computation
# ---------------------------------------------------------------------------

def compute_per_edge(ts, ws2, sigma2, beta, X_standardised, y, n_indiv,
                    top_k_trees=None, pop_assignment_fn=None):
    """
    Build the per-edge dataframe.

    Args:
        ts                 : tree sequence
        ws2                : per-tree variance weights (length S)
        sigma2             : σ²
        beta               : β̂ per SNP (length p), or None
        X_standardised     : (p, n) standardised genotype matrix, or None
        y                  : phenotype vector (length n)
        n_indiv            : number of diploid individuals
        top_k_trees        : restrict to trees with K largest ŵ²_s (None = all
                             trees with ŵ² > 0)
        pop_assignment_fn  : callable (tree, child_node) → pop_label, or None

    Returns:
        rows: list of dicts with per-edge fields.
    """
    n = n_indiv
    var_y = float(np.var(y))

    # residual = y - X β  (or just y if no β)
    if beta is not None and X_standardised is not None:
        Xbeta = X_standardised.T @ beta            # (n,)
        residual = y - Xbeta
    else:
        residual = y - y.mean()

    # Solve V z = residual
    print(f"  building V and solving V z = (y - Xβ)  (n={n})", file=sys.stderr)
    z, V, K = compute_V_and_solve(ts, ws2, sigma2, residual)

    # Per-edge β_e (SNP side)
    beta_e_sums = (
        per_edge_beta_sum(ts, beta) if beta is not None else {}
    )

    # Sample-to-individual map: tskit's default sample numbering pairs
    # haplotype 2i and 2i+1 as individual i.
    sample_nodes = list(ts.samples())
    sample_to_indiv = {nd: idx // 2 for idx, nd in enumerate(sample_nodes)}

    # Decide which trees to process
    S = ts.num_trees
    tree_order = np.argsort(ws2)[::-1]
    if top_k_trees is not None:
        tree_order = tree_order[:top_k_trees]
    keep = set(int(t) for t in tree_order if ws2[t] > 0)

    rows = []
    print(f"  processing {len(keep)}/{S} trees …", file=sys.stderr)

    for tree_idx, tree in enumerate(ts.trees()):
        if tree_idx not in keep:
            continue
        w = float(ws2[tree_idx])
        if w == 0.0:
            continue

        c_per_node       = {}
        S_z_per_node     = {}
        sum_csq_per_node = {}

        # post-order: build c arrays
        c_per_node, _ = per_node_descendant_arrays(tree, sample_to_indiv, n)
        S_z_per_node       = per_node_S_z(c_per_node, z)
        sum_csq_per_node   = per_node_sum_csquared(c_per_node)

        node_times = ts.tables.nodes.time
        pos_lo = tree.interval.left
        pos_hi = tree.interval.right

        # iterate edges = (child, parent) pairs in this tree
        for u in tree.nodes():
            p = tree.parent(u)
            if p == tskit.NULL:
                continue
            length_t   = float(node_times[p] - node_times[u])   # in time units
            child_age  = float(node_times[u])
            sum_csq    = sum_csq_per_node[u]
            S_z        = S_z_per_node[u]
            beta_e     = float(beta_e_sums.get((tree_idx, u), 0.0))
            # Per implementation (V = σ²I + Σ ŵ²·G_s): E[α_e|y] = ŵ²_s · ℓ_e · S_z[u].
            # No extra σ² factor — that lives only in the residual variance, not in
            # the per-tree variance scaling.
            E_alpha    = w * length_t * S_z

            dh2_arg   = (E_alpha ** 2) * sum_csq / (n * var_y)
            dh2_snp   = (beta_e ** 2) * sum_csq / (n * var_y)
            dh2_total = ((E_alpha + beta_e) ** 2) * sum_csq / (n * var_y)

            row = {
                "tree_idx":      tree_idx,
                "child_node":    int(u),
                "length":        length_t,
                "child_age":     child_age,
                "n_haps":        int(c_per_node[u].sum()),
                "sum_csq":       int(sum_csq),
                "beta_e":        beta_e,
                "E_alpha":       E_alpha,
                "dh2_snp":       dh2_snp,
                "dh2_arg":       dh2_arg,
                "dh2_total":     dh2_total,
                "position_lo":   float(pos_lo),
                "position_hi":   float(pos_hi),
            }
            if pop_assignment_fn is not None:
                row["pop_label"] = pop_assignment_fn(tree, u)
            rows.append(row)

    return rows, z, V


# ---------------------------------------------------------------------------
# Self-consistency checks (per design doc § Validation)
# ---------------------------------------------------------------------------

def self_consistency_arg(rows, ts, ws2, sigma2, z):  # noqa: ARG001 (sigma2 unused — kept for symmetry with caller)
    """
    Check that Σ_e u_e^ARG[i] across all edges in a tree matches the per-tree
    ARG prediction σ² · ŵ²_s · G_s · z (the per-tree posterior contribution).

    Strict per-tree identity that should hold to machine precision.
    """
    by_tree = {}
    for r in rows:
        by_tree.setdefault(r["tree_idx"], []).append(r)
    sample_nodes = list(ts.samples())
    sample_to_indiv = {nd: idx // 2 for idx, nd in enumerate(sample_nodes)}
    n_indiv = max(sample_to_indiv.values()) + 1

    max_residual_norm = 0.0
    for tree_idx, tree in enumerate(ts.trees()):
        if tree_idx not in by_tree:
            continue
        w = ws2[tree_idx]
        if w == 0.0:
            continue
        c_per_node, _ = per_node_descendant_arrays(tree, sample_to_indiv, n_indiv)
        # Build u_e^ARG[i] = c[u(e)][i] · E[α_e | y]; sum over edges
        u_arg = np.zeros(n_indiv, dtype=float)
        for r in by_tree[tree_idx]:
            u_arg += r["E_alpha"] * c_per_node[r["child_node"]].astype(float)
        # Compare to σ² · w · (G_s · z) using streaming_combined_kinship for this tree
        # G_s · z = streaming_combined_kinship(ts, samples, [1 if t==tree_idx else 0]) · z
        single = np.zeros(ts.num_trees)
        single[tree_idx] = 1.0
        Gs_python = hap_to_indiv_2d(
            streaming_combined_kinship(ts, sample_nodes, single), n_indiv
        )
        # streaming + hap_to_indiv_2d returns half the matvec G_s (see
        # compute_V_and_solve for the same fix).  Scale up to match.
        Gs_matvec = 2.0 * Gs_python
        expected = w * (Gs_matvec @ z)         # ws²_s · G_s_matvec · z
        residual_norm = float(np.linalg.norm(u_arg - expected))
        max_residual_norm = max(max_residual_norm, residual_norm)
    return max_residual_norm


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case-prefix",    required=True,
                    help="Path prefix to case.* simulation files.")
    ap.add_argument("--sargent-prefix", required=True,
                    help="Path prefix to sargent's est.* output files.")
    ap.add_argument("--out", default=None,
                    help="Output TSV for per-edge data (default: <sargent>.edge_h2.tsv).")
    ap.add_argument("--top-k-trees",  type=int, default=None,
                    help="Restrict to trees with the K largest ŵ²_s (default: all).")
    ap.add_argument("--top-k-edges",  type=int, default=None,
                    help="In the output TSV keep only the K largest Δh²_total edges "
                         "(default: keep all).")
    ap.add_argument("--age-bin-edges", default="0,5000,20000,inf",
                    help="Comma-separated edges in generations for per-age aggregation "
                         "(default: 0,5000,20000,inf).")
    ap.add_argument("--no-bayesR-fit", action="store_true",
                    help="Skip loading β̂ and X (SNP side); compute ARG side only.")
    ap.add_argument("--check-self-consistency", action="store_true",
                    help="Run the per-tree ARG self-consistency check and report.")
    args = ap.parse_args()

    case = args.case_prefix
    est  = args.sargent_prefix
    out  = args.out or (est + ".edge_h2.tsv")

    print(f"== edge_h2 ==")
    print(f"  case: {case}")
    print(f"  est:  {est}")

    # Tree sequence
    ts = tskit.load(case + ".trees")
    print(f"  ts:   n_trees={ts.num_trees}  n_sites={ts.num_sites}  "
          f"n_samples={ts.num_samples}")

    # Phenotype
    y = load_phen(case + ".phen")
    n_indiv = len(y)
    print(f"  y:    n={n_indiv}  var={float(np.var(y)):.4g}")

    # sargent outputs
    ws2 = load_ws2(est + ".w2.tsv")
    if len(ws2) != ts.num_trees:
        print(f"  WARNING: len(ws2)={len(ws2)} != n_trees={ts.num_trees}")
    summary = load_summary(est + ".summary.json")
    sigma2  = float(summary.get("sigma2", 1.0))
    print(f"  ws2:  S={len(ws2)}  σ²={sigma2:.4g}")

    # Optional BayesR side
    if args.no_bayesR_fit or not os.path.exists(est + ".beta.tsv"):
        beta = None
        X    = None
        print(f"  β:    skipped (no .beta.tsv or --no-bayesR-fit)")
    else:
        beta = load_beta(est + ".beta.tsv")
        # Need n_variants from the .bim file
        bim = case + ".bim"
        if not os.path.exists(bim):
            print(f"  WARNING: {bim} missing; cannot compute Xβ; SNP side disabled")
            beta = None
            X    = None
        else:
            with open(bim) as f:
                n_variants = sum(1 for _ in f)
            X = compute_X_standardised(case, n_indiv, n_variants)
            print(f"  β:    p={len(beta)}  ‖β‖={float(np.linalg.norm(beta)):.4g}")

    # Build per-edge data
    rows, z, _ = compute_per_edge(
        ts, ws2, sigma2, beta, X, y, n_indiv,
        top_k_trees=args.top_k_trees,
    )

    print(f"  rows: {len(rows)} edges")

    # Optional top-K filter
    if args.top_k_edges is not None and len(rows) > args.top_k_edges:
        rows.sort(key=lambda r: r["dh2_total"], reverse=True)
        rows = rows[: args.top_k_edges]

    # Write per-edge TSV
    fieldnames = ["tree_idx", "child_node", "length", "child_age",
                  "n_haps", "sum_csq", "beta_e", "E_alpha",
                  "dh2_snp", "dh2_arg", "dh2_total",
                  "position_lo", "position_hi"]
    if rows and "pop_label" in rows[0]:
        fieldnames.append("pop_label")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})
    print(f"  wrote {out}")

    # Per-age-bin aggregation
    if args.age_bin_edges:
        edges = [float("inf") if e.strip().lower() == "inf" else float(e)
                 for e in args.age_bin_edges.split(",")]
        if len(edges) >= 2:
            bins = list(zip(edges[:-1], edges[1:]))
            print()
            print(f"  Per-age-bin Δh² (SNP / ARG / total):")
            for lo, hi in bins:
                # Δh²_b = ‖Σ_e in b  E[α_e|y]·c[u(e)] + β_e·c[u(e)]‖² / (n · Var(y))
                # but the rows have per-edge u_e^SNP / u_e^ARG only as scalars;
                # to get the bin-level squared norm we need the per-individual sums.
                # For the printed quick summary, just sum the per-edge Δh² (a lower
                # bound on the true bin-level total — cross-terms between edges
                # within the bin are dropped).  Document this in the table.
                dh2_arg   = sum(r["dh2_arg"]   for r in rows
                                if lo <= r["child_age"] < hi)
                dh2_snp   = sum(r["dh2_snp"]   for r in rows
                                if lo <= r["child_age"] < hi)
                dh2_total = sum(r["dh2_total"] for r in rows
                                if lo <= r["child_age"] < hi)
                lo_s = f"{lo:.0f}" if lo != float("inf") else "inf"
                hi_s = f"{hi:.0f}" if hi != float("inf") else "inf"
                print(f"    [{lo_s:>5s}, {hi_s:>5s}):  "
                      f"snp={dh2_snp:.4f}  arg={dh2_arg:.4f}  total={dh2_total:.4f}")
            print(f"    (these are sums of per-edge Δh² within bin; cross-edge"
                  f" terms within a bin are not included — see design doc § Open extensions.)")

    # Optional self-consistency check
    if args.check_self_consistency:
        print()
        print(f"  Self-consistency: per-tree ARG sum vs σ²·ŵ²·G_s·z")
        max_res = self_consistency_arg(rows, ts, ws2, sigma2, z)
        print(f"    max ‖per-edge sum − σ²·ŵ²·G_s·z‖ = {max_res:.3e}")
        print(f"    (should be small — machine epsilon up to single-precision noise)")


if __name__ == "__main__":
    main()
