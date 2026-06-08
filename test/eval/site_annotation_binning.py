#!/usr/bin/env python
"""Block 2.2 Phase C — site-annotation edge binning for Acts I-II.

Provides `bin_by_site_annotation`, a bin_fn factory compatible with
kinship_assembly.build_per_bin_kinship, for SITE (genomic-position) annotations
like exon / promoter / intergenic.  Kept in a separate module (not the iMac's
canonical kinship_assembly.py) until the convention below is signed off, to
avoid divergence.

CONVENTIONS (first cut — flagged for Adam, see design tensions):

  * One bin per edge (library constraint: bin_fn returns a single bin_id, and
    build_per_bin_kinship weights every edge by ℓ_e/2).  An edge is assigned to
    the **majority annotation of the mutations it carries** in that tree.
  * Mutation-free edges → a dedicated `unannotated` bin.  KEPT (not dropped):
    these branches are genealogical background the ARG sees but typed SNPs do
    not — the very signal that should make ARG-REML beat SNP-REML in Regimes
    2/3.  A driver can include or exclude this bin as needed (exclude → reduce
    toward SNP-REML; include → full ARG).

DESIGN TENSIONS (need Adam's call before this goes into kinship_assembly.py):

  1. Weighting.  The library uses the eGRM (ℓ_e/2) per-edge weight, NOT a
     per-mutation count.  The design doc's "per-mutation rank-1 decomposition"
     and "reduces to SNP-REML when every mutation is typed" only coincide with
     (ℓ_e/2) in a limit.  For an exact Regime-1 calibration to SNP-REML we may
     need a per-mutation-count variant (weight ∝ #mutations on the edge in bin
     b).  This factory uses the eGRM weight (one bin/edge) as the default.
  2. Majority vs fractional.  The design rejected edge-fraction (double-counts
     shared edges) and edge-majority (ambiguous) but chose per-mutation, which
     the one-bin-per-edge interface cannot express.  Majority is the closest
     single-bin approximation; almost all edges carry 0-1 mutations so the
     approximation only bites on the rare multi-annotation edge.
"""
from collections import Counter, defaultdict

import numpy as np
import tskit


def read_annot_intervals(case_prefix):
    """Read a Phase-A <case>.annot_map.tsv into sorted (start, end, ann) arrays."""
    s, e, a = [], [], []
    with open(case_prefix + ".annot_map.tsv") as f:
        next(f)
        for line in f:
            x, y, z = line.rstrip("\n").split("\t")
            s.append(float(x)); e.append(float(y)); a.append(z)
    return np.array(s), np.array(e), np.array(a)


def make_position_annotation(starts, ends, anns, default="unannotated"):
    """Return pos -> annotation lookup over a sorted, contiguous interval map."""
    def pos_ann(pos):
        i = int(np.searchsorted(starts, pos, side="right") - 1)
        if 0 <= i < len(ends) and pos < ends[i]:
            return anns[i]
        return default
    return pos_ann


def bin_by_site_annotation(ts, starts, ends, anns, unannotated_label="unannotated"):
    """bin_fn factory: assign each edge to its mutations' majority annotation.

    Args:
        ts                : tskit.TreeSequence
        starts/ends/anns  : sorted contiguous annotation interval map (as from
                            read_annot_intervals on a <case>.annot_map.tsv)
        unannotated_label : bin for edges with no mutation (kept separate)

    Returns a bin_fn(tree, child, c_u) -> bin_id compatible with
    kinship_assembly.build_per_bin_kinship.
    """
    pos_ann = make_position_annotation(starts, ends, anns, default=unannotated_label)
    breaks = ts.breakpoints(as_array=True)
    site_pos = ts.tables.sites.position
    mut_site = ts.tables.mutations.site
    mut_node = ts.tables.mutations.node

    # (tree_index, node) -> Counter of mutation annotations
    edge_ann = defaultdict(Counter)
    for m in range(ts.num_mutations):
        pos = float(site_pos[mut_site[m]])
        ti = int(np.searchsorted(breaks, pos, side="right") - 1)
        edge_ann[(ti, int(mut_node[m]))][pos_ann(pos)] += 1

    def bin_fn(tree, child, c_u):  # noqa: ARG001 (c_u unused; binning is positional)
        c = edge_ann.get((tree.index, child))
        if not c:
            return unannotated_label
        return c.most_common(1)[0][0]

    return bin_fn


if __name__ == "__main__":
    # Smoke test on an existing Phase-A annotation sim (block2_2_v1/pilot).
    import os
    import kinship_assembly as ka
    cp = os.path.expanduser("~/test-sargent/block2_2_v1/pilot/case")
    if not os.path.exists(cp + ".annot_map.tsv"):
        raise SystemExit(f"no annot_map at {cp}; point at a --annotation-genes case")
    ts = tskit.load(cp + ".trees")
    n = ts.num_individuals
    s, e, a = read_annot_intervals(cp)
    bin_fn = bin_by_site_annotation(ts, s, e, a)
    Kb = ka.build_per_bin_kinship(ts, bin_fn, n)
    print(f"bins: {sorted(Kb)}")
    for b in sorted(Kb):
        K = Kb[b]
        print(f"  {b:>14}: meandiag={np.mean(np.diag(K)):.3g} "
              f"trace={np.trace(K):.3g} PSD_mineig={np.linalg.eigvalsh(K).min():.2e}")
    Ksum = sum(Kb.values())
    Kall = ka.build_per_bin_kinship(ts, lambda *_: "all", n)["all"]   # fast K_ARG
    rel = np.linalg.norm(Ksum - Kall) / np.linalg.norm(Kall)
    print(f"sum-rule (Σ_b K_b incl. unannotated == K_ARG) rel err = {rel:.2e}")
