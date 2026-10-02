#!/usr/bin/env python3
"""
simulate.py — generate synthetic ARG + genotype + phenotype data for SARGENT testing.

Model:  y = X β + u + ε
  X β : fixed effects (BayesR discrete mixture)
   u  : random effects drawn from tree-branch Brownian motion with per-tree ws²,
        weighted by each tree's genomic span (the span-weighted eGRM)
   ε  : i.i.d. N(0, σ²)
   V  = Σ_s ws²[s] · span_s · G_s + σ² I

Outputs (all with --out prefix):
  .trees          tskit tree sequence
  .phen           phenotype file (one value per line, individual order matches tree seq)
  .bed/.bim/.fam  PLINK genotype files for fixed-effect SNPs
  .true_beta.tsv  causal SNP indices and true betas (after scaling)
  .true_ws2.tsv   per-tree true ws² values (after scaling)
  .params.json    all simulation parameters and realised h²
"""

import argparse
import json
import math
import os
import sys

import msprime
import numpy as np
import tskit


# ---------------------------------------------------------------------------
# BED writer (PLINK SNP-major format, no external dependency)
# ---------------------------------------------------------------------------

def write_plink(prefix, dosage, positions, n_indiv, n_variants):
    """Write PLINK BED/BIM/FAM files.

    dosage: shape (n_variants, n_indiv), integer 0/1/2 or NaN for missing
    """
    bytes_per_snp = (n_indiv + 3) // 4

    with open(prefix + ".bed", "wb") as f:
        f.write(bytes([0x6c, 0x1b, 0x01]))   # magic + SNP-major mode
        for j in range(n_variants):
            row = dosage[j]
            buf = bytearray(bytes_per_snp)
            for i in range(n_indiv):
                v = row[i]
                if np.isnan(v):
                    code = 1   # missing: 01
                elif v <= 0:
                    code = 0   # hom ref: 00
                elif v == 1:
                    code = 2   # het: 10
                else:
                    code = 3   # hom alt: 11
                buf[i // 4] |= code << ((i % 4) * 2)
            f.write(buf)

    with open(prefix + ".bim", "w") as f:
        for j in range(n_variants):
            f.write(f"1\tsnp{j}\t0\t{positions[j]}\tA\tC\n")

    with open(prefix + ".fam", "w") as f:
        for i in range(n_indiv):
            f.write(f"ind{i} ind{i} 0 0 0 -9\n")


# ---------------------------------------------------------------------------
# Random-effect simulation via Brownian motion on tree branches
#
# Covariance: Cov(u_s[i], u_s[j]) = ws²[s] · G_s[i,j]
# where G_s[i,j] = 0.5 · Σ_b bl(b) · A_b(i) · A_b(j)
# and A_b(k) = number of haplotypes of k in the subtree below branch b.
#
# Draw z_b ~ N(0, ws² · bl(b)) independently for each branch b, then
#   node_val[leaf] = Σ_{b ancestor of leaf} z_b   (cumulative path sum)
#   u_s[i] = (1/√2) · (node_val[hap_a(i)] + node_val[hap_b(i)])
# ---------------------------------------------------------------------------

def sim_random_effect(ts, hap_ranks, n_indiv, ws2_input, rng):
    """Return u vector (length n_indiv) for all trees.

    ws2_input: either
      (a) array of ws² values, one per tree (length ts.num_trees), in which case
          every edge in tree s uses the same per-tree ws²_s; or
      (b) dict {(tree_idx, child_node) → ws²_e} for per-edge architectures
          (Block 2.7.1).  Edges absent from the dict get ws²_e = 0.

    hap_ranks: shape (n_indiv, 2) — sample ranks for each individual's two haplotypes.
    """
    u = np.zeros(n_indiv)
    node_val = np.zeros(ts.num_nodes)
    samples = ts.samples()
    per_edge = isinstance(ws2_input, dict)

    for tree_idx, tree in enumerate(ts.trees()):
        # Span-weight each tree's variance contribution so the generative kernel
        # is the span-weighted eGRM (V = Sum_s ws2_s * span_s * G_s), matching the
        # estimator (gall.c / variational.c weight G_s by span_s).  A branch's
        # increment variance scales by the tree's genomic span.  The caller
        # renormalises ws2 so realised var(u) hits the target h2, so only the
        # RELATIVE per-tree weighting matters here (raw span is fine).
        span = tree.interval.right - tree.interval.left
        if per_edge:
            # Skip trees with no active edges
            has_any = any(
                ws2_input.get((tree_idx, int(u_)), 0.0) > 0.0
                for u_ in tree.nodes() if tree.parent(u_) != tskit.NULL
            )
            if not has_any:
                continue
            ws2_tree_lookup = ws2_input  # dict lookup per edge
        else:
            ws2_tree = ws2_input[tree_idx]
            if ws2_tree == 0.0:
                continue

        node_val[:] = 0.0
        # Preorder: root first, then children — parent value is ready before children
        for nd in tree.nodes(order="preorder"):
            parent = tree.parent(nd)
            if parent == tskit.NULL:
                node_val[nd] = 0.0
            else:
                if per_edge:
                    ws2_e = ws2_tree_lookup.get((tree_idx, int(nd)), 0.0)
                else:
                    ws2_e = ws2_tree
                if ws2_e == 0.0:
                    node_val[nd] = node_val[parent]
                    continue
                bl = tree.branch_length(nd)
                node_val[nd] = node_val[parent] + rng.normal(0.0, math.sqrt(ws2_e * bl * span))

        # Accumulate into u: (1/√2) · (val[hap0] + val[hap1])
        h0 = hap_ranks[:, 0]
        h1 = hap_ranks[:, 1]
        u += (1.0 / math.sqrt(2.0)) * (node_val[samples[h0]] + node_val[samples[h1]])

    return u


def build_ws2_architecture(kind, n_trees, n_blocks, block_frac, contrast,
                           rng=None, ts=None):
    """Phase 10 / Block 2.6 — construct a per-tree ws² profile.

    'diffuse'    : constant ws² across all trees (flat landscape).
    'sparse'     : floor + `n_blocks` evenly-spaced, contiguous runs of trees set
                   high. `block_frac` of all trees are active; block:floor ratio
                   is `contrast`.
    'scattered'  : floor + `block_frac · n_trees` individually scattered active
                   trees (Bernoulli-style draw, no contiguity). `n_blocks` ignored.
    'age_recent' : binary spike on the youngest third of trees by mean branch age.
    'age_middle' : binary spike on the middle third of trees by mean branch age.
    'age_ancient': binary spike on the oldest third of trees by mean branch age.
                   Requires `ts` (tskit.TreeSequence).

    Per-edge variants (Block 2.7.1; 2026-06-05).  Return a *dict* keyed by
    (tree_idx, child_node) → ws² rather than a 1-D per-tree array.  Downstream
    sim_random_effect dispatches on the return type.  Intended for the WHEN
    axis where per-edge age is the natural binning, since trees span many
    ages and per-tree-mean-age is too coarse a summary.

    'age_recent_edges' : binary spike on edges whose midpoint age is in the
                         youngest tertile of edge ages (across the whole
                         tree sequence).
    'age_middle_edges' : same, for the middle tertile of edge ages.
    'age_ancient_edges': same, for the oldest tertile of edge ages.

    Absolute scale is arbitrary — caller renormalises ws2_vec so realised var(u)
    matches h2_arg; only the relative profile matters downstream.
    """
    if n_trees <= 0:
        return np.zeros(0, dtype=float)
    if kind == "diffuse":
        return np.ones(n_trees, dtype=float)

    # Per-edge age tertile architectures must be checked BEFORE the per-tree
    # age_* branch because both share the "age_" prefix.
    if kind.startswith("age_") and kind.endswith("_edges"):
        pass  # handled below, after the per-tree age block
    elif kind.startswith("age_"):
        if ts is None:
            raise ValueError("age_* architectures require ts=<TreeSequence>")
        # Import here to avoid adding tskit as a top-level dependency for callers
        # that don't use age modes.
        import os as _os, sys as _sys
        _eval = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                              "..", "test", "eval")
        if _eval not in _sys.path:
            _sys.path.insert(0, _eval)
        from arch_metrics import per_tree_mean_branch_age
        ages = per_tree_mean_branch_age(ts)
        # Use equal-count tertile edges (robust to any age distribution).
        q33, q67 = float(np.nanpercentile(ages, 33.33)), float(np.nanpercentile(ages, 66.67))
        if kind == "age_recent":
            active = ages <= q33
        elif kind == "age_middle":
            active = (ages > q33) & (ages <= q67)
        elif kind == "age_ancient":
            active = ages > q67
        else:
            raise ValueError(f"Unknown age architecture: {kind!r}")
        ws2 = np.where(active & ~np.isnan(ages), 1.0, 0.0)
        return ws2

    if kind.startswith("age_") and kind.endswith("_edges"):
        # Per-edge age tertile architectures (Block 2.7.1 design).  Return a
        # dict {(tree_idx, child_node) → ws²_e}.  sim_random_effect dispatches.
        if ts is None:
            raise ValueError(f"{kind!r} requires ts=<TreeSequence>")
        node_times = ts.tables.nodes.time
        # Collect (tree_idx, child_node, midpoint_age) for every edge.
        edge_entries = []
        for ti, tree in enumerate(ts.trees()):
            for u in tree.nodes():
                p = tree.parent(u)
                if p == tskit.NULL:
                    continue
                age_child  = float(node_times[u])
                age_parent = float(node_times[p])
                mid        = 0.5 * (age_child + age_parent)
                edge_entries.append((ti, int(u), mid))
        if not edge_entries:
            return {}
        mids = np.array([e[2] for e in edge_entries], dtype=float)
        q33  = float(np.nanpercentile(mids, 33.33))
        q67  = float(np.nanpercentile(mids, 66.67))
        if kind == "age_recent_edges":
            active_fn = lambda m: m <= q33
        elif kind == "age_middle_edges":
            active_fn = lambda m: (m > q33) & (m <= q67)
        elif kind == "age_ancient_edges":
            active_fn = lambda m: m > q67
        else:
            raise ValueError(f"Unknown per-edge age architecture: {kind!r}")
        ws2_per_edge = {}
        for ti, u, mid in edge_entries:
            ws2_per_edge[(ti, u)] = 1.0 if active_fn(mid) else 0.0
        return ws2_per_edge

    floor = 1.0 / max(contrast, 1.0)
    ws2 = np.full(n_trees, floor, dtype=float)
    if kind == "scattered":
        if rng is None:
            rng = np.random.default_rng(0)
        n_active = max(1, int(round(block_frac * n_trees)))
        n_active = min(n_active, n_trees)
        active_idx = rng.choice(n_trees, size=n_active, replace=False)
        ws2[active_idx] = 1.0
        return ws2
    # sparse (n_blocks contiguous, evenly spaced)
    n_active = max(n_blocks, int(round(block_frac * n_trees)))
    width = max(1, n_active // max(n_blocks, 1))
    centres = np.linspace(0, n_trees, n_blocks + 2)[1:-1]   # avoid genome edges
    for c in centres:
        lo = int(round(c - width / 2.0))
        lo = min(max(lo, 0), n_trees - width)
        ws2[lo:lo + width] = 1.0
    return ws2


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_gene_annotations(length, exon_frac, mean_exon_bp, promoter_bp,
                           tau_exon, tau_promoter_ratio, rng):
    """Generate a simple gene-structured annotation partition (Block 2.2).

    Single-exon "genes" of length ~lognormal(mean_exon_bp) are placed along
    [0, length) with exponential intergenic gaps sized so exonic bases occupy
    ~exon_frac of the genome.  Each exon has a `promoter_bp` promoter immediately
    upstream; everything not exon/promoter is intergenic.  Per-base causal
    enrichment: exon = tau_exon, promoter = tau_exon * tau_promoter_ratio,
    intergenic = 1.0 (baseline).

    Returns:
      ann_ids   : ['exon', 'promoter', 'intergenic']  (indices 0,1,2)
      ann_tau   : np.ndarray of per-annotation tau (same order)
      intervals : list of (start_bp, end_bp, ann_idx) for exon/promoter blocks,
                  sorted by start (intergenic is the implicit complement).
    """
    ann_ids = ["exon", "promoter", "intergenic"]
    ann_tau = np.array([tau_exon, tau_exon * tau_promoter_ratio, 1.0], dtype=float)
    EXON, PROMOTER = 0, 1
    # Mean bp consumed per gene so exonic bases average to exon_frac of the genome.
    period   = mean_exon_bp / max(exon_frac, 1e-9)
    mean_gap = max(period - promoter_bp - mean_exon_bp, 1.0)
    sigma    = 0.6                                   # lognormal spread of exon length
    mu       = math.log(max(mean_exon_bp, 1.0)) - 0.5 * sigma * sigma
    intervals = []
    x = 0.0
    while True:
        x += rng.exponential(mean_gap)               # intergenic gap before promoter
        p_start = int(round(x))
        e_start = p_start + int(round(promoter_bp))
        exon_len = int(round(min(max(rng.lognormal(mu, sigma), 30.0),
                                 10.0 * mean_exon_bp)))
        e_end = e_start + exon_len
        if e_end >= length:
            break
        intervals.append((p_start, e_start, PROMOTER))
        intervals.append((e_start, e_end,  EXON))
        x = float(e_end)
    return ann_ids, ann_tau, intervals


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-indiv",     type=int,   default=200,  help="diploid individuals (training set)")
    ap.add_argument("--n-test",      type=int,   default=0,    help="additional held-out test individuals (0 = no test split)")
    ap.add_argument("--length",      type=float, default=5e6,  help="sequence length (bp)")
    ap.add_argument("--recomb-rate", type=float, default=1e-8, help="recombination rate per bp per gen")
    ap.add_argument("--mut-rate",    type=float, default=1e-8, help="mutation rate per bp per gen")
    ap.add_argument("--Ne",          type=float, default=1e4,  help="effective population size")
    ap.add_argument("--demography",  type=str,   default="panmictic",
                    choices=["panmictic", "admixed", "transfer", "stdpopsim"],
                    help="population model. 'panmictic' = single constant-size pop; "
                         "'admixed' = two diverged sources + recent admixture pulse, "
                         "giving local-ancestry heterogeneity along the genome (the "
                         "Phase 10 optimistic regime, where different regions carry "
                         "genuinely different relatedness so per-segment ws² can help); "
                         "'transfer' = two diverged pops A and B (no admixture pulse) "
                         "with training drawn from A + optional B fraction (--mix-fraction-B) "
                         "and test drawn entirely from B (Block 1.1 cross-ancestry PRS); "
                         "'stdpopsim' = use a named stdpopsim HomSap demographic model "
                         "(requires --stdpopsim-model, --train-pop, --test-pop).")
    ap.add_argument("--stdpopsim-model", type=str, default=None,
                    help="stdpopsim HomSap demographic model name "
                         "(only with --demography stdpopsim). "
                         "E.g. 'OutOfAfrica_3G09', 'AmericanAdmixture_4B18', "
                         "'OutOfAfricaArchaicAdmixture_5R19'.")
    ap.add_argument("--train-pop", type=str, default=None,
                    help="population name for training samples "
                         "(only with --demography stdpopsim). "
                         "Must match a population ID in the stdpopsim model "
                         "(e.g. 'CEU', 'YRI', 'EUR'). "
                         "--mix-fraction-B controls target-pop injection into training.")
    ap.add_argument("--test-pop", type=str, default=None,
                    help="population name for test samples "
                         "(only with --demography stdpopsim). "
                         "For admixed-cohort test, use an admixed deme name "
                         "(e.g. 'ADMIX', 'ASW', 'MXL', 'PEL').")
    ap.add_argument("--train-pop-time", type=float, default=0.0,
                    help="sampling time in generations before present for train "
                         "individuals (--demography stdpopsim; default 0 = present). "
                         "Non-zero enables sampling ancient or archaic haplotypes as "
                         "training data (scenario 2: archaic-ancestry mapping).")
    ap.add_argument("--test-pop-time", type=float, default=0.0,
                    help="sampling time in generations before present for test "
                         "individuals (--demography stdpopsim; default 0 = present). "
                         "Set to ~1000 to sample Neanderthal haplotypes near "
                         "Vindija timing (~28 kya) from OutOfAfricaArchaicAdmixture_5R19; "
                         "gives F_ST ≈ 0.54 vs CEU, well above the modern-human range. "
                         "Must be < the population's split/merge time in the model.")
    ap.add_argument("--admix-split-time", type=float, default=4000.0,
                    help="generations since the two source pops split (--demography admixed; "
                         "larger ⇒ stronger structure, F_ST ≈ 1-exp(-T/2Ne))")
    ap.add_argument("--admix-time",  type=float, default=20.0,
                    help="generations since the admixture pulse (--demography admixed)")
    ap.add_argument("--admix-alpha", type=float, default=0.5,
                    help="ancestry fraction from source A (--demography admixed)")
    ap.add_argument("--transfer-split-time", type=float, default=3500.0,
                    help="generations since A and B split (--demography transfer; "
                         "3500 gen with Ne=1e4 gives F_ST ≈ 0.15, the EUR-vs-EAS regime).")
    ap.add_argument("--mix-fraction-B", type=float, default=0.0,
                    help="fraction of TRAIN sampled from population B in --demography transfer "
                         "(test is always 100%% B). 0.0 = pure A→B transfer; 0.05 = trace mix; "
                         "0.20 = moderate mix. Block 1.1 ladder = 0/0.05/0.20.")
    ap.add_argument("--h2-snp",    type=float, default=0.3,  help="target fixed-effect heritability")
    ap.add_argument("--h2-arg",   type=float, default=0.3,  help="target random-effect heritability")
    # --- Frequency architecture of the causal set (sargent-paper SIMULATION_SPECS §5) ---
    # Both default to the values that reproduce prior behaviour exactly.
    ap.add_argument("--causal-maf-min", type=float, default=0.0,
                    help="minimum pooled MAF for a variant to be ELIGIBLE as causal "
                         "(default 0.0 = no floor, prior behaviour).  Without a floor, "
                         "causal variants are drawn uniformly from the site-frequency "
                         "spectrum and are therefore mostly singletons/doubletons: their "
                         "realised contribution is dominated by which few individuals "
                         "happen to carry them, which destabilises the realised h2 and "
                         "makes train and test heritability diverge.  Published analyses "
                         "impose a floor even when rare variants are the subject "
                         "(Schoech 2019: MAF>0.07%%; Gazal 2018: 0.5%%; Zeng 2018: 1%%).  "
                         "NOTE this changes the RNG stream: the eligible set feeding "
                         "rng.choice differs, so floor-on and floor-off runs are "
                         "different draws, not the same draw perturbed.")
    ap.add_argument("--alpha", type=float, default=-1.0,
                    help="frequency dependence of PER-ALLELE effect variance: "
                         "Var(beta_j) proportional to [2p(1-p)]^alpha, so the variance a "
                         "variant EXPLAINS scales as [2p(1-p)]^(alpha+1).  "
                         "alpha=-1 (default) is the standardised-genotype model -- every "
                         "causal variant explains the same variance regardless of "
                         "frequency -- and is what you get implicitly from 'standardise "
                         "genotypes, draw effects i.i.d.'.  alpha=0 means no frequency "
                         "dependence of per-allele effects.  Empirical estimates cluster "
                         "near -0.4 (Schoech 2019: -0.38; Gazal 2018: -0.37; Zeng 2018 "
                         "corroborating), under which MAF<1%% variants explain <10%% of "
                         "SNP-heritability.  Applied as a deterministic reweighting of "
                         "the drawn effects, so it consumes NO random numbers and "
                         "alpha=-1 is bit-identical to prior behaviour.  Composes "
                         "multiplicatively with --age-tau-coef; both act before the "
                         "h2_snp normalisation, so order is irrelevant.")
    ap.add_argument("--n-causal",    type=int,   default=None,
                    help="number of causal SNPs (fixed effects).  If unset, defaults to "
                         "max(20, int(--causal-fraction * n_variants)).  An explicit value "
                         "takes priority over --causal-fraction.")
    ap.add_argument("--causal-fraction", type=float, default=0.01,
                    help="fraction of polymorphic variants that are causal SNPs "
                         "(default 0.01 = 1%%, matches SBayesR's canonical scenario for human "
                         "polygenic traits; typical GWAS-polygenic range is 0.01-0.10).  "
                         "Overridden by an explicit --n-causal.")
    ap.add_argument("--min-causal", type=int, default=20,
                    help="floor on the number of causal SNPs when deriving the count "
                         "from --causal-fraction (default 20, so very small simulations "
                         "don't get pathologically few causal variants).  When this floor "
                         "binds, the realised causal fraction exceeds --causal-fraction; a "
                         "warning is printed and the realised fraction is recorded in "
                         "params.json.  Set --min-causal 0 to disable the floor.")
    # Block 2.1 — age-stratified effect-size architecture.
    # Scale per-variant β² by (age_j / age_ref)^tau before the h²_snp
    # normalisation. tau=0 disables (default, current behaviour).
    # tau<0 puts more variance on RECENT variants (rare-variant burden /
    # purifying-selection-driven architecture; e.g., O'Connor / Zeng / Schoech).
    # tau>0 puts more variance on ANCIENT variants (drift-dominated background).
    # The total h²_snp is preserved by the downstream scale_f normalisation
    # (line 488); tau redistributes the variance among causal SNPs by age.
    ap.add_argument("--age-tau-coef", type=float, default=0.0,
                    help="age-dependent effect-size exponent. β² scales as "
                         "(age / age_ref)^tau.  Default 0 (no age dependence). "
                         "Negative: young variants larger (purifying-selection "
                         "model). Positive: ancient variants larger (drift "
                         "background).")
    ap.add_argument("--age-tau-ref", type=float, default=10000.0,
                    help="reference age (generations) for age-tau scaling. "
                         "Choose near the median variant age to keep magnitudes "
                         "comparable across runs. Default 10000 generations "
                         "(near 4·Ne for Ne=1e4).")
    ap.add_argument("--bed-coverage", type=float, default=1.0,
                    help="fraction of polymorphic SNPs to keep in the BED file (0 < f <= 1). "
                         "Hidden SNPs still contribute to y via their true β and remain in the "
                         "ARG/.trees output. Default: 1.0 (no hiding).")
    ap.add_argument("--bed-maf-min", type=float, default=0.0,
                    help="minimum MAF for a variant to enter the BED (typed) panel "
                         "(default 0 = no filter).  Set ~0.05 to make the panel "
                         "common-SNP / genotyping-array-like, so rare causal variants "
                         "are untyped and only the ARG can see them (Block 2.2: this is "
                         "what creates the flavour-A-vs-B gap).  Composes with "
                         "--bed-coverage as a plain AND, so f_typed retains its meaning "
                         "as the fraction of the ELIGIBLE pool that is typed.")
    ap.add_argument("--bed-freq-weight", type=float, default=0.0,
                    help="frequency-dependent ascertainment of the typed panel: variant "
                         "j is sampled with probability proportional to [2p(1-p)]^w "
                         "(default 0 = uniform, prior behaviour).  A uniform draw leaves "
                         "the panel's frequency spectrum INVARIANT to --bed-coverage, so "
                         "a simulated 5%% 'array' has the same MAF profile as full "
                         "sequencing -- the opposite of a real genotyping array, which is "
                         "designed from common variation and becomes progressively more "
                         "common-biased as it thins.  A single fixed w reproduces that "
                         "behaviour automatically, because weighted sampling from a "
                         "rare-dominated pool concentrates on the common tail as the "
                         "panel shrinks: at w=1 (weight proportional to genotypic "
                         "variance) the median typed MAF runs ~0.09 at f_typed=0.3 and "
                         "~0.29 at f_typed=0.05.  Unlike the uniform draw, which keeps a "
                         "Binomial number of variants, the weighted draw keeps an EXACT "
                         "count of round(f_typed * p).")
    ap.add_argument("--bed-snps-per-mb", type=float, default=0.0,
                    help="size the typed panel in SNPs per Mb instead of as a fraction "
                         "of the simulated variant count (default 0 = use --bed-coverage). "
                         "Overrides --bed-coverage when > 0.  --bed-coverage's meaning "
                         "depends on how many variants the coalescent happened to produce, "
                         "so it is comparable neither across region lengths nor to a real "
                         "platform; SNPs/Mb is both.  Reference points, genome-wide counts "
                         "over ~3 Gb: ~1400/Mb = Omni5 (4.3M), ~830/Mb = Omni2.5 (2.5M), "
                         "~250/Mb = the common arrays (GSA 650K, OmniExpress 700K, UK "
                         "Biobank Axiom 820K), ~80/Mb = 240K, ~30/Mb = 90K.  These "
                         "simulations produce ~3800 variants/Mb, so a common array is "
                         "f_typed ~ 0.066 -- far sparser than the axis Fig 2 used.")
    ap.add_argument("--causal-not-typed", action="store_true",
                    help="force causal variants OUT of the typed panel, so the coverage "
                         "axis measures LD tagging rather than direct genotyping.  Without "
                         "this, --causal-maf-min and --bed-freq-weight both select on "
                         "frequency and the panel preferentially contains exactly the "
                         "variants made causal: measured at f_typed=0.3, the panel directly "
                         "genotyped causal variants carrying 99%% of the genetic variance, "
                         "leaving nothing for tagging to explain.  Mirrors the treatment "
                         "the --h2-unobs class already gets.")
    ap.add_argument("--seed",        type=int,   default=42,   help="random seed")
    ap.add_argument("--out",         type=str,   default="sim", help="output prefix")
    # Phase 6b: SV-flavoured unobserved-variants class.  Variants are sampled
    # from the polymorphic-variant pool and EXCLUDED from the BED, then assigned
    # effect sizes from a heavier-tailed distribution.  Their contribution to y
    # is missing heritability that SARGENT's ARG random component should recover
    # but the SNP fit cannot.
    ap.add_argument("--h2-unobs",    type=float, default=0.0,
                    help="target heritability for unobserved (SV-style) variants "
                         "(default 0 = disabled)")
    ap.add_argument("--n-unobs-causal", type=int, default=10,
                    help="number of unobserved-class causal variants (only used "
                         "when --h2-unobs > 0)")
    ap.add_argument("--unobs-effect-dist", type=str, default="gaussian",
                    choices=["gaussian", "laplace"],
                    help="effect-size distribution for unobserved variants (default: gaussian)")
    ap.add_argument("--unobs-effect-scale", type=float, default=5.0,
                    help="effect-size SD multiplier (relative to BayesR slab) for "
                         "unobserved variants -- larger = more SV-like (default: 5.0)")
    ap.add_argument("--unobs-maf-min", type=float, default=0.0,
                    help="minimum MAF for SV-class carve-out (default 0.0 = no lower filter). "
                         "Restricting to rare variants weakens LD with common BED SNPs and "
                         "makes SV variance genuinely 'missing' from BayesR.")
    ap.add_argument("--unobs-maf-max", type=float, default=0.5,
                    help="maximum MAF for SV-class carve-out (default 0.5 = no upper filter). "
                         "Set e.g. 0.05 to restrict to rare-variant SV class.")
    # ---- structure of the unobserved class --------------------------------
    # WHY THESE EXIST.  Under infinite sites an untyped variant sits on exactly
    # ONE branch, and a COMMON variant's branch is old and widely carried, so
    # nearby typed SNPs tag it almost perfectly.  Measured on these cells: at 830
    # SNPs/Mb, of 1384 untyped variants exactly ONE has MAF > 5% -- and it is well
    # tagged.  A carve-out selected by FREQUENCY is therefore untyped but not
    # untagged: the SNP arms lose nothing by not typing it and the ARG has no
    # advantage to demonstrate.  That is what the first SV experiment measured
    # (+0.010) and what a 20-cell repeat at h2_arg = 0 reproduced
    # (+0.0095 +/- 0.0080).  Both knobs below break the one-variant-one-branch
    # correspondence that makes "untyped" imply "tagged".
    #
    # Defaults (1 / 0.0) are OFF and the legacy code path is kept literally, so
    # existing seeds reproduce bit for bit.
    # ---- the unobserved class as a SEPARATE MUTATIONAL PROCESS -------------
    # WHY.  The --h2-unobs machinery above draws the unobserved class from the
    # SNV pool and then hides those variants from the panel.  But hiding an SNV
    # does not make it untaggable: OTHER SNVs on the SAME BRANCH are still in the
    # panel and tag it exactly, which is why the carve-out was only ever worth
    # ~+0.01 (20 cells, +0.0095 +/- 0.0080).
    #
    # A structural variant is not a hidden SNV.  It arises by its own mutational
    # process, at its own (much lower) rate, and is absent from an SNV callset by
    # construction.  Whether it is tagged is then EMERGENT: it depends on whether
    # an SNV happened to land on the same branch -- a Poisson hit with mean
    # mu_snv * l_b * span_b.  Measured on these cells, ~80% of branch weight
    # carries no SNV, so most SVs land untagged without anything being imposed.
    #
    # --sv-rate turns the class on.  Human de novo SV estimates (~0.16 per genome
    # per generation over ~3.1 Gb) give ~5e-11 /bp/gen, about 1/250 of the SNV
    # rate; at 1.25e-8 that is --sv-rate 5e-11.  Default 0 = off, and the legacy
    # code path is kept literally so existing seeds reproduce bit for bit.
    # ---- RECURRENCE: one variant, k independent origins ---------------------
    # WHY.  Validation of --sv-rate showed a single-origin SV is usually WELL
    # TAGGED (median max r2 0.88-1.00): mutations land in proportion to
    # branch_length x span, so they concentrate on LONG branches -- exactly the
    # ones that also caught an SNV.  A carrier set that is a UNION of branches is
    # the only kind no single mutation can match, because no mutation sits on
    # more than one branch.
    #
    # Biologically this is ordinary: NAHR-recurrent CNVs, VNTR length classes and
    # STR expansions all arise repeatedly at one locus on different haplotype
    # backgrounds.  The composite variant is built from k branches drawn in
    # proportion to branch length at one locus -- i.e. exactly where independent
    # mutations would have fallen -- and its carriers are the union of their
    # descendants.  It is never in the BED panel: an SNV callset contains no
    # VNTR length alleles.
    #
    # This is a SIMULATION construct only.  The inference methods are blind to
    # it; they see a phenotype and a panel, which is what keeps the comparison
    # fair.  Default 1 = off, legacy path kept literally.
    ap.add_argument("--unobs-recurrence-zipf", type=float, default=0.0,
                    help="draw each unobserved variant's number of origins K from "
                         "a Zipf(a) distribution with MODE 1 instead of fixing it "
                         "(default 0 = off; try a=2.0). More defensible than a "
                         "fixed K: most structural variants are single-origin and "
                         "a minority (VNTRs, NAHR-hotspot CNVs) recur heavily, so "
                         "phi emerges from the mixture rather than being dialled, "
                         "and there is no arbitrary K to justify. Capped by "
                         "--unobs-recurrence-kmax.")
    ap.add_argument("--unobs-recurrence-kmax", type=int, default=20,
                    help="cap on K when --unobs-recurrence-zipf is used (default 20). "
                         "Zipf has an infinite mean for a<=2, so a cap is required.")
    ap.add_argument("--unobs-recurrence-k", type=int, default=1,
                    help="build each unobserved causal variant as a RECURRENT "
                         "variant with K independent origins at one locus, its "
                         "carriers being the union of K clades (default 1 = off). "
                         "Models NAHR-recurrent CNVs, VNTRs and STR expansions. "
                         "K > 1 is the only construct that makes a variant "
                         "untaggable by ANY single marker.")
    ap.add_argument("--sv-rate", type=float, default=0.0,
                    help="per-bp per-generation mutation rate for a SEPARATE "
                         "structural-variant class, overlaid on the same ARG and "
                         "never placed in the BED panel (default 0 = off). When "
                         "> 0 the --h2-unobs class is drawn from these variants "
                         "rather than from hidden SNVs, so 'untyped' is no longer "
                         "the same as 'tagged'. Suggested: 5e-11 (~1/250 of the "
                         "SNV rate, from human de novo SV estimates).")
    ap.add_argument("--unobs-series-k", type=int, default=1,
                    help="allelic series: place the unobserved class as groups of K "
                         "DISTINCT causal variants within one locus, with same-sign effects, "
                         "instead of independent singletons (default 1 = off).  The aggregate "
                         "carrier set is then a union of K branches, which no single marker "
                         "can tag -- the reason gene-based burden tests exist -- while a "
                         "branch-additive ARG random effect still spans it.")
    ap.add_argument("--unobs-series-kb", type=float, default=25.0,
                    help="locus half-width in kb for --unobs-series-k (default 25): the K "
                         "variants of a series are drawn within this distance of the anchor.")
    ap.add_argument("--unobs-mask-kb", type=float, default=0.0,
                    help="locally-sparse markers: exclude typed SNPs within this many kb of an "
                         "unobserved causal variant (default 0 = off).  Models SVs arising in "
                         "repetitive regions where arrays and short-read calling under-perform. "
                         "The panel SIZE is held fixed -- markers are redistributed elsewhere -- "
                         "so overall density is unchanged and only the LOCAL depletion differs.")
    # BayesR-style mixture for fixed effects (sigma² relative to var_pheno=1)
    ap.add_argument("--sigma2-mix",  type=float, nargs="+",
                    default=[0.0, 1e-4, 1e-3, 1e-2],
                    help="fixed-effect prior variances (first must be 0)")
    ap.add_argument("--pi-mix",      type=float, nargs="+", default=None,
                    help="fixed-effect mixture proportions (default: uniform)")
    # Segment-weight mixture (tau² relative to var_pheno/S)
    # Two presets: 'sparse' (default, spike-and-slab style, ~15% active)
    #              'diffuse' (uniform 4-component, ~75% active — legacy behaviour)
    # If --tau2-frac / --phi-mix are given explicitly, --regime is ignored.
    ap.add_argument("--regime",      type=str,   default="sparse",
                    choices=["sparse", "diffuse"],
                    help="ws² mixture regime (default: sparse)")
    ap.add_argument("--tau2-frac",   type=float, nargs="+", default=None,
                    help="segment-weight fractions relative to var_pheno/S (first must be 0; overrides --regime)")
    ap.add_argument("--phi-mix",     type=float, nargs="+", default=None,
                    help="segment-weight mixture proportions (overrides --regime)")
    # Phase 10 oracle diagnostic: directly construct the per-tree ws² profile,
    # overriding the --regime mixture draw.  Lets us contrast a flat landscape
    # (control: oracle ≈ uniform) with a sparse contiguous large-effect landscape
    # (the regime that should favour per-segment weights, if anything does).
    ap.add_argument("--ws2-architecture", type=str, default=None,
                    choices=["diffuse", "sparse", "scattered",
                             "age_recent", "age_middle", "age_ancient",
                             "age_recent_edges", "age_middle_edges",
                             "age_ancient_edges"],
                    help="directly build the ws² profile (Phase 10 / Block 2.6), "
                         "overriding the --regime mixture draw.  "
                         "'diffuse' = flat constant ws²; "
                         "'sparse' = floor + n_blocks evenly-spaced contiguous blocks; "
                         "'scattered' = floor + block_frac·n_trees individually scattered "
                         "active trees (n_blocks ignored). The 'scattered' option is the "
                         "fused-lasso robustness check: confirms that any per-segment "
                         "advantage isn't an artifact of the prior matching contiguous truth. "
                         "'age_{recent,middle,ancient}' = binary spike on the tertile of "
                         "TREES by mean branch age. "
                         "'age_{recent,middle,ancient}_edges' = binary spike on the tertile "
                         "of EDGES by midpoint age (Block 2.7.1; finer-grained per-edge "
                         "ground truth used by Fig 5 Panel B).")
    ap.add_argument("--ws2-n-blocks", type=int, default=3,
                    help="number of contiguous large-effect blocks (--ws2-architecture sparse)")
    ap.add_argument("--ws2-block-frac", type=float, default=0.05,
                    help="fraction of trees that are active/in-blocks (--ws2-architecture sparse)")
    ap.add_argument("--ws2-block-contrast", type=float, default=1000.0,
                    help="block:floor ws² ratio (--ws2-architecture sparse; large ⇒ "
                         "nearly all variance in the blocks)")
    # Block 2.2 — gene-structured functional annotations + per-category causal
    # enrichment.  Single-exon genes + upstream promoters + intergenic.  Causal
    # beta^2 is scaled by tau_{annotation} (mirrors --age-tau-coef), so the
    # realised per-annotation h^2 is proportional to base_frac * tau.
    ap.add_argument("--annotation-genes", action="store_true",
                    help="lay down a gene-structured annotation partition "
                         "(single-exon genes + upstream promoters + intergenic) "
                         "and enrich causal variants per category.  Writes "
                         "<out>.annotations.tsv (per-variant labels) and "
                         "<out>.annot_map.tsv (genomic partition).")
    ap.add_argument("--gene-exon-frac", type=float, default=0.015,
                    help="target fraction of bases that are exonic (default 0.015).")
    ap.add_argument("--gene-mean-exon-bp", type=float, default=150.0,
                    help="mean exon length in bp (lognormal spread; default 150).")
    ap.add_argument("--gene-promoter-bp", type=float, default=100.0,
                    help="promoter length in bp immediately upstream of each exon "
                         "(default 100; promoter base fraction follows as "
                         "~exon_frac * promoter_bp / mean_exon_bp).")
    ap.add_argument("--gene-tau-exon", type=float, default=10.0,
                    help="per-base causal enrichment of exons (default 10).")
    ap.add_argument("--gene-tau-promoter-ratio", type=float, default=1.0/3.0,
                    help="promoter enrichment as a fraction of the exon enrichment "
                         "(default 1/3); intergenic is baseline (tau=1).")
    args = ap.parse_args()

    # Apply --regime presets if mixture not given explicitly.
    # Sparse defaults are designed to be identifiable at n=100: only ~5% of trees
    # are active, so each active tree carries a large per-tree variance share,
    # large enough to rise above Hutchinson noise.
    if args.tau2_frac is None and args.phi_mix is None:
        if args.regime == "sparse":
            # 3-component spike-and-slab: 95% null, 4% medium slab, 1% large slab
            args.tau2_frac = [0.0, 0.1, 1.0]
            args.phi_mix   = [0.95, 0.04, 0.01]
        else:   # diffuse
            args.tau2_frac = [0.0, 0.01, 0.1, 1.0]
            args.phi_mix   = [0.25, 0.25, 0.25, 0.25]
    elif args.tau2_frac is None or args.phi_mix is None:
        sys.exit("ERROR: --tau2-frac and --phi-mix must be given together")

    if not (0.0 < args.causal_fraction <= 1.0):
        sys.exit("ERROR: --causal-fraction must satisfy 0 < f <= 1")
    if args.h2_snp + args.h2_arg + args.h2_unobs > 1.0:
        sys.exit("ERROR: h2_snp + h2_arg + h2_unobs must be ≤ 1")
    if args.h2_snp < 0 or args.h2_arg < 0 or args.h2_unobs < 0:
        sys.exit("ERROR: heritabilities must be non-negative")
    if not (0.0 < args.bed_coverage <= 1.0):
        sys.exit("ERROR: --bed-coverage must satisfy 0 < f <= 1")
    if args.h2_unobs > 0 and args.n_unobs_causal <= 0:
        sys.exit("ERROR: --h2-unobs > 0 requires --n-unobs-causal > 0")
    if not (0.0 <= args.unobs_maf_min <= args.unobs_maf_max <= 0.5):
        sys.exit("ERROR: require 0 <= --unobs-maf-min <= --unobs-maf-max <= 0.5")

    rng = np.random.default_rng(args.seed)

    # ---- 1. Simulate ARG and mutations -------------------------------------
    n_total = args.n_indiv + args.n_test
    # Transfer demography: split n_train between A and B; n_test is all B.
    # Layout in the sampled ts: [n_train_A A indivs, n_train_B B indivs, n_test B indivs],
    # which preserves the existing convention that train = ts[:n_train], test = ts[n_train:].
    n_train_A = args.n_indiv
    n_train_B = 0
    if args.demography in ("transfer", "stdpopsim"):
        n_train_B = int(round(args.mix_fraction_B * args.n_indiv))
        n_train_A = args.n_indiv - n_train_B

    if args.demography == "stdpopsim":
        if not args.stdpopsim_model:
            sys.exit("ERROR: --stdpopsim-model is required with --demography stdpopsim")
        if not args.train_pop:
            sys.exit("ERROR: --train-pop is required with --demography stdpopsim")
        if not args.test_pop and (args.n_test > 0 or n_train_B > 0):
            sys.exit("ERROR: --test-pop is required with --demography stdpopsim "
                     "when n_test > 0 or mix-fraction-B > 0")

    print(f"Simulating ARG: n_train={args.n_indiv} + n_test={args.n_test} "
          f"= {n_total} individuals, L={args.length:.0f} bp ...", flush=True)
    if args.demography == "admixed":
        dem = msprime.Demography()
        dem.add_population(name="ANC", initial_size=args.Ne)
        dem.add_population(name="A",   initial_size=args.Ne)
        dem.add_population(name="B",   initial_size=args.Ne)
        dem.add_population(name="ADM", initial_size=args.Ne)
        dem.add_admixture(time=args.admix_time, derived="ADM",
                          ancestral=["A", "B"],
                          proportions=[args.admix_alpha, 1.0 - args.admix_alpha])
        dem.add_population_split(time=args.admix_split_time,
                                 derived=["A", "B"], ancestral="ANC")
        ts = msprime.sim_ancestry(
            samples={"ADM": n_total},
            sequence_length=args.length,
            recombination_rate=args.recomb_rate,
            demography=dem,
            random_seed=args.seed,
        )
    elif args.demography == "transfer":
        dem = msprime.Demography()
        dem.add_population(name="ANC", initial_size=args.Ne)
        dem.add_population(name="A",   initial_size=args.Ne)
        dem.add_population(name="B",   initial_size=args.Ne)
        dem.add_population_split(time=args.transfer_split_time,
                                 derived=["A", "B"], ancestral="ANC")
        # msprime returns samples in population-dict order: A first, then B.
        # Sampling A=n_train_A, B=n_train_B+n_test gives the layout above.
        sample_spec = {"A": n_train_A, "B": n_train_B + args.n_test}
        # Drop pops with zero samples (msprime rejects 0-count entries on some versions)
        sample_spec = {p: n for p, n in sample_spec.items() if n > 0}
        ts = msprime.sim_ancestry(
            samples=sample_spec,
            sequence_length=args.length,
            recombination_rate=args.recomb_rate,
            demography=dem,
            random_seed=args.seed,
        )
        print(f"  transfer composition: train A={n_train_A} B={n_train_B}, "
              f"test B={args.n_test}", flush=True)
    elif args.demography == "stdpopsim":
        try:
            import stdpopsim as _stdpopsim
        except ImportError:
            sys.exit("ERROR: stdpopsim is not installed. "
                     "Install with: pip install stdpopsim")
        _species = _stdpopsim.get_species("HomSap")
        try:
            _sps_model = _species.get_demographic_model(args.stdpopsim_model)
        except Exception as e:
            avail = [m.id for m in _species.demographic_models]
            sys.exit(f"ERROR: unknown stdpopsim model {args.stdpopsim_model!r}: {e}\n"
                     f"Available HomSap models: {avail}")
        # Validate population names against the model's population list, and
        # warn when an archaic population (default_sampling_time=None) is
        # requested at time 0 — msprime will sample from a near-empty lineage.
        _pop_names = {p.name for p in _sps_model.populations}
        for _pn in filter(None, [args.train_pop, args.test_pop]):
            if _pn not in _pop_names:
                sys.exit(f"ERROR: population {_pn!r} not in {args.stdpopsim_model!r}. "
                         f"Available: {sorted(_pop_names)}")
        for _pn, _pt, _label in [
                (args.train_pop, args.train_pop_time, "train"),
                (args.test_pop,  args.test_pop_time,  "test")]:
            if _pn is None:
                continue
            _pop_obj = next((p for p in _sps_model.populations
                             if p.name == _pn), None)
            if (_pop_obj is not None
                    and getattr(_pop_obj, "default_sampling_time", 0) is None
                    and _pt == 0.0):
                print(
                    f"  WARNING: {_pn!r} has no defined present-day sampling "
                    f"time (default_sampling_time=None in stdpopsim). "
                    f"Sampling at t=0 may produce empty or unexpected results. "
                    f"Consider --{_label}-pop-time N "
                    f"(e.g. --{_label}-pop-time 1000 for Neanderthal ~28 kya).",
                    file=sys.stderr,
                )
        # Layout (same convention as transfer):
        #   [0 .. n_train_A-1]          from train_pop  (at train_pop_time)
        #   [n_train_A .. n_total-1]    from test_pop   (at test_pop_time)
        #
        # Use msprime.SampleSet unconditionally so non-zero sampling times
        # work for archaic populations (e.g. Neanderthal@t=1000 gen).
        _n_test_pop = n_train_B + args.n_test
        _sample_sets = [msprime.SampleSet(
            n_train_A, population=args.train_pop,
            time=args.train_pop_time)]
        if _n_test_pop > 0:
            _sample_sets.append(msprime.SampleSet(
                _n_test_pop, population=args.test_pop,
                time=args.test_pop_time))
        ts = msprime.sim_ancestry(
            samples=_sample_sets,
            demography=_sps_model.model,
            sequence_length=args.length,
            recombination_rate=args.recomb_rate,
            random_seed=args.seed,
        )
        _train_time_str = (f"@t={args.train_pop_time:.0f}gen"
                           if args.train_pop_time > 0 else "")
        _test_time_str  = (f"@t={args.test_pop_time:.0f}gen"
                           if args.test_pop_time  > 0 else "")
        print(f"  stdpopsim {args.stdpopsim_model}: "
              f"train_pop={args.train_pop}{_train_time_str} ({n_train_A}), "
              f"test_pop={args.test_pop}{_test_time_str} ({_n_test_pop})",
              flush=True)
        if n_train_B > 0:
            print(f"  mixed training: {n_train_A} from {args.train_pop} + "
                  f"{n_train_B} from {args.test_pop} "
                  f"({args.mix_fraction_B:.0%})", flush=True)
    else:
        ts = msprime.sim_ancestry(
            samples=n_total,
            sequence_length=args.length,
            recombination_rate=args.recomb_rate,
            population_size=args.Ne,
            random_seed=args.seed,
        )
    ts = msprime.sim_mutations(ts, rate=args.mut_rate, random_seed=args.seed + 1)

    # Second, INDEPENDENT mutational process for the structural-variant class.
    # Overlaid on the same genealogy, so an SV is a real mutation on a real
    # branch -- not an effect attached to a branch by fiat. Identified by the
    # site positions that did not exist before the overlay; a collision with an
    # existing SNV position makes the site multiallelic and it is dropped by the
    # biallelic filter below, which is the conservative outcome.
    sv_positions = set()
    if args.sv_rate > 0.0:
        _pre = set(ts.tables.sites.position.tolist())
        ts = msprime.sim_mutations(ts, rate=args.sv_rate,
                                   random_seed=args.seed + 7, keep=True)
        sv_positions = set(ts.tables.sites.position.tolist()) - _pre
        print(f"  SV class: {len(sv_positions)} sites at rate {args.sv_rate:g} "
              f"({args.sv_rate/args.mut_rate:.4f} x the SNV rate)", flush=True)

    # F_ST between train_pop and test_pop (stdpopsim only, computed on the
    # full tree sequence using tskit's Hudson estimator).
    fst_train_test = None
    if args.demography == "stdpopsim" and args.train_pop and args.test_pop:
        _pop_name_to_id = {}
        for _p in ts.populations():
            _pname = (_p.metadata.get("name", str(_p.id))
                      if _p.metadata else str(_p.id))
            _pop_name_to_id[_pname] = _p.id
        _train_pid = _pop_name_to_id.get(args.train_pop, -1)
        _test_pid  = _pop_name_to_id.get(args.test_pop,  -1)
        _train_nodes, _test_nodes = [], []
        for _ind in ts.individuals():
            _pid = ts.node(_ind.nodes[0]).population
            if _pid == _train_pid:
                _train_nodes.extend(_ind.nodes)
            elif _pid == _test_pid:
                _test_nodes.extend(_ind.nodes)
        if _train_nodes and _test_nodes:
            try:
                fst_train_test = round(float(ts.Fst(
                    [_train_nodes, _test_nodes])), 6)
                print(f"  F_ST ({args.train_pop} vs {args.test_pop}): "
                      f"{fst_train_test:.4f}", flush=True)
            except Exception as _e:
                print(f"  F_ST computation failed: {_e}", flush=True)
    print(f"  {ts.num_individuals} individuals, {ts.num_trees} trees, "
          f"{ts.num_sites} sites", flush=True)

    n_indiv  = ts.num_individuals
    n_trees  = ts.num_trees
    samples  = ts.samples()          # node IDs for all sample haplotypes

    # Per-individual haplotype ranks (sample indices, not node IDs)
    # ts.individual(i).nodes gives node IDs; map to sample ranks
    sample_rank = {node: rank for rank, node in enumerate(samples)}
    hap_ranks = np.array(
        [[sample_rank[ts.individual(i).nodes[0]],
          sample_rank[ts.individual(i).nodes[1]]]
         for i in range(n_indiv)],
        dtype=int,
    )  # shape (n_indiv, 2)

    # ---- 2. Collect polymorphic biallelic SNPs ----------------------------
    print("Collecting polymorphic variants ...", flush=True)
    dosage_list = []   # each entry: dosage vector (n_indiv,)
    pos_list    = []
    age_list    = []   # per-variant age (generations), from mutation time
    is_sv_list  = []   # True when the variant came from the SV overlay

    # Pre-build a tree iterator so we can resolve mutation.node → branch parent
    # time when the mutation's own `time` field is unknown (older tskit).
    for var in ts.variants():
        alleles = [a for a in var.alleles if a is not None]
        if len(alleles) != 2:
            continue
        g = var.genotypes           # length n_hap, sample-rank indexed
        d = g[hap_ranks[:, 0]] + g[hap_ranks[:, 1]]   # 0/1/2 dosage
        # require polymorphic
        if d.min() == d.max():
            continue
        dosage_list.append(d.astype(float))
        pos_list.append(int(var.site.position))
        # Age of the variant = mutation time (generations before present).
        # Use the first mutation at the site (biallelic ⇒ exactly one).
        mut = var.site.mutations[0]
        age = float(mut.time) if not math.isnan(mut.time) else float("nan")
        if math.isnan(age):
            # Fallback: time of the parent of the mutation node in the local tree.
            tree = ts.at(var.site.position)
            parent_node = tree.parent(mut.node)
            age = float(tree.time(parent_node)) if parent_node != tskit.NULL else 0.0
        age_list.append(age)
        is_sv_list.append(var.site.position in sv_positions)

    p = len(dosage_list)
    if p == 0:
        sys.exit("ERROR: no polymorphic biallelic variants found")
    print(f"  {p} polymorphic variants", flush=True)

    dosage = np.stack(dosage_list)   # (p, n_indiv)
    positions = pos_list
    ages = np.array(age_list, dtype=float)  # per-variant age in generations
    is_sv = np.array(is_sv_list, dtype=bool)
    if is_sv.any():
        print(f"  {int(is_sv.sum())} of {p} polymorphic variants are SV-class",
              flush=True)

    # ---- Block 2.2: gene-structured annotation partition -------------------
    # Build the exon/promoter/intergenic partition (dedicated rng so the
    # phenotype rng is unperturbed), label each variant by the interval it falls
    # in, and set var_tau[j] = the enrichment multiplier applied to beta_j^2.
    if args.annotation_genes:
        ann_ids, ann_tau, gene_intervals = build_gene_annotations(
            args.length, args.gene_exon_frac, args.gene_mean_exon_bp,
            args.gene_promoter_bp, args.gene_tau_exon,
            args.gene_tau_promoter_ratio,
            np.random.default_rng(args.seed + 2_000_003))
        n_ann = len(ann_ids)
        INTERGENIC = n_ann - 1
        starts = np.array([s for (s, e, a) in gene_intervals], dtype=float)
        ends   = np.array([e for (s, e, a) in gene_intervals], dtype=float)
        anns   = np.array([a for (s, e, a) in gene_intervals], dtype=int)
        posarr = np.array(positions, dtype=float)
        var_ann = np.full(p, INTERGENIC, dtype=int)
        if len(starts):
            k = np.searchsorted(starts, posarr, side="right") - 1
            inblk = (k >= 0) & (posarr < ends[np.clip(k, 0, len(ends) - 1)])
            var_ann[inblk] = anns[k[inblk]]
        var_tau = ann_tau[var_ann]
        ann_base = np.zeros(n_ann)
        for (s, e, a) in gene_intervals:
            ann_base[a] += (e - s)
        ann_base /= args.length
        ann_base[INTERGENIC] = max(0.0, 1.0 - ann_base.sum())
        print(f"  annotation: {len(gene_intervals)//2} genes; base fracs "
              f"exon={ann_base[0]:.4f} promoter={ann_base[1]:.4f} "
              f"intergenic={ann_base[2]:.4f}", flush=True)
    else:
        ann_ids = None
        var_ann = None
        var_tau = None

    # Decide which variants will be written to the BED file (the "observed panel").
    # Variants are still in the ARG and can carry true β; this controls what Stage 1
    # sees as X.
    #
    # Use a SEPARATE rng stream (derived from args.seed but a different generator)
    # so the bed-coverage draw does not perturb the phenotype-generation rng state.
    # This guarantees that fixing --seed produces identical β, ws², and y across
    # different --bed-coverage values — only the BED panel itself changes.
    # Phase 6b: when h2_unobs > 0, carve out unobserved-causal variants FIRST so
    # they are guaranteed disjoint from the SNP causal set and force-excluded
    # from the BED below.  Done here (before bed_mask) so the bed_mask draw
    # can exclude them.  Default (h2_unobs=0) keeps the legacy behaviour.
    #
    # MAF-conditional sampling: when --unobs-maf-max < 0.5, restrict the SV
    # pool to rare variants.  Rare SVs have weak LD with common BED SNPs,
    # so the carved-out variance is genuinely "missing" from BayesR rather
    # than leaked back through tag-SNP LD.  Default range [0, 0.5] keeps
    # legacy behaviour (no MAF filter).
    # ---- recurrent (union-of-branches) unobserved variants -----------------
    recurrence_built = False
    if args.h2_unobs > 0 and (args.unobs_recurrence_k > 1 or args.unobs_recurrence_zipf > 0):
        K = args.unobs_recurrence_k
        rec_rng = np.random.default_rng(args.seed + 99_991)
        node_to_indiv = np.full(ts.num_nodes, -1, dtype=np.int64)
        for _i in range(n_indiv):
            for _nd in ts.individual(_i).nodes:
                node_to_indiv[_nd] = _i
        new_dos, new_pos, new_age, k_used = [], [], [], []
        tries = 0
        while len(new_dos) < args.n_unobs_causal and tries < 200 * args.n_unobs_causal:
            tries += 1
            if args.unobs_recurrence_zipf > 0:
                K = int(min(rec_rng.zipf(args.unobs_recurrence_zipf),
                            args.unobs_recurrence_kmax))
            locus = float(rec_rng.uniform(0, args.length))
            tree = ts.at(locus)
            # Draw a target frequency for the composite, then take K origins of
            # COMPARABLE size (~target/K each).  Sampling origins by branch
            # length alone picks long branches, so one large clade dominates the
            # union and a marker tagging that clade tags the union -- measured
            # median max r2 0.68-0.90 even at K=8.  Real recurrent CNVs / VNTRs
            # arise repeatedly on modest, recent backgrounds, so the union is
            # several comparable rare clades; that is both the realistic model
            # and the one no single marker can match.
            f_target = float(rec_rng.uniform(max(args.unobs_maf_min, 0.02),
                                             min(args.unobs_maf_max, 0.35)))
            f_each = f_target / K
            n_hap = 2 * n_indiv
            cand, wts = [], []
            for u in tree.nodes():
                par = tree.parent(u)
                if par == tskit.NULL:
                    continue
                bl = tree.time(par) - tree.time(u)
                if bl <= 0:
                    continue
                fu = tree.num_samples(u) / n_hap
                if not (0.4 * f_each <= fu <= 1.8 * f_each):
                    continue
                cand.append(u); wts.append(bl)
            if len(cand) < K:
                continue
            wts = np.asarray(wts, dtype=float); wts /= wts.sum()
            pick = rec_rng.choice(len(cand), size=K, replace=False, p=wts)
            hap = set()
            for j in pick:
                hap.update(tree.samples(cand[j]))          # union of the K clades
            d = np.zeros(n_indiv)
            for h in hap:
                ii = node_to_indiv[h]
                if ii >= 0:
                    d[ii] += 1.0                            # per-haplotype count
            f = d.sum() / (2.0 * n_indiv)
            m = min(f, 1.0 - f)
            if m < max(args.unobs_maf_min, 1e-9) or m > args.unobs_maf_max:
                continue
            if d.min() == d.max():
                continue
            new_dos.append(d); new_pos.append(int(locus)); k_used.append(K)
            new_age.append(float(np.mean([tree.time(tree.parent(cand[j])) for j in pick])))
        if len(new_dos) < args.n_unobs_causal:
            sys.exit(f"ERROR: could only build {len(new_dos)} of "
                     f"{args.n_unobs_causal} recurrent variants at K={K} in MAF "
                     f"[{args.unobs_maf_min}, {args.unobs_maf_max}]; widen the MAF "
                     f"range, lower K, or lengthen --length.")
        dosage    = np.vstack([dosage, np.stack(new_dos)])
        positions = list(positions) + new_pos
        ages      = np.concatenate([ages, np.asarray(new_age, dtype=float)])
        is_sv     = np.concatenate([is_sv, np.ones(len(new_dos), dtype=bool)])
        unobs_idx = np.arange(p, p + len(new_dos), dtype=int)
        p         = dosage.shape[0]
        n_unobs_causal = len(unobs_idx)
        unobs_series   = np.arange(n_unobs_causal)
        unobs_set      = set(unobs_idx.tolist())
        recurrence_built = True
        # The typed-causal pool excludes the composites, exactly as the
        # non-recurrent path excludes its unobs_idx.
        remaining = np.array([j for j in range(p) if j not in unobs_set])
        _f = dosage[unobs_idx].mean(axis=1) / 2.0
        _kdesc = (f"K~Zipf({args.unobs_recurrence_zipf}) capped {args.unobs_recurrence_kmax}, "
                  f"mean {np.mean(k_used):.1f}, max {int(np.max(k_used))}, "
                  f"{100*np.mean(np.asarray(k_used)==1):.0f}% single-origin"
                  if args.unobs_recurrence_zipf > 0 else f"K={K}")
        print(f"  recurrence: {n_unobs_causal} variants, {_kdesc}, "
              f"MAF {np.min(np.minimum(_f,1-_f)):.3f}-{np.max(np.minimum(_f,1-_f)):.3f}",
              flush=True)

    if args.h2_unobs > 0 and not recurrence_built:
        n_unobs_causal = min(args.n_unobs_causal, p)
        freq = dosage.mean(axis=1) / 2.0
        maf  = np.minimum(freq, 1.0 - freq)
        eligible = np.where((maf >= args.unobs_maf_min) &
                            (maf <= args.unobs_maf_max))[0]
        # With a real SV class present, the unobserved variants ARE the SVs.
        # Their untaggability is then emergent (did an SNV land on the same
        # branch?) rather than imposed by hiding an SNV that has twins.
        if is_sv.any():
            eligible = np.array([j for j in eligible if is_sv[j]], dtype=int)
            print(f"  unobserved class drawn from the SV pool: "
                  f"{len(eligible)} eligible SV variants", flush=True)
        if len(eligible) < n_unobs_causal:
            sys.exit(f"ERROR: only {len(eligible)} variants in MAF "
                     f"[{args.unobs_maf_min}, {args.unobs_maf_max}], need "
                     f"{n_unobs_causal}.  Lengthen --length or widen MAF range.")
        if args.unobs_series_k > 1:
            # ---- allelic series -------------------------------------------
            # Draw n_loci anchors, then take K-1 further ELIGIBLE variants within
            # +/- unobs_series_kb of each.  The K variants are distinct sites (an
            # allelic series, not one recurrent mutation), so each still sits on
            # its own branch; what changes is that the class's variance is carried
            # by a UNION of K branches per locus.  No single marker correlates
            # with that union, while a branch-additive random effect spans it.
            # unobs_series[] records group membership so the effect draw can give
            # one locus a common sign (see the beta_unobs block).
            K = args.unobs_series_k
            n_loci = max(1, n_unobs_causal // K)
            win = args.unobs_series_kb * 1e3
            posv = np.asarray(positions, dtype=float)
            picked, groups, used = [], [], set()
            anchors = rng.permutation(eligible)
            for a in anchors:
                if len(picked) >= n_loci * K:
                    break
                if a in used:
                    continue
                near = eligible[(np.abs(posv[eligible] - posv[a]) <= win)]
                near = np.array([j for j in near if j not in used])
                if len(near) < K:
                    continue          # locus too sparse to host a series
                sel = rng.choice(near, size=K, replace=False)
                gid = len(groups)
                for j in sel:
                    picked.append(int(j)); groups.append(gid); used.add(int(j))
                groups_last = gid
            if len(picked) < K:
                sys.exit(f"ERROR: --unobs-series-k {K} could not find any locus with "
                         f"{K} eligible variants within +/-{args.unobs_series_kb} kb.  "
                         f"Widen --unobs-series-kb, lower --unobs-series-k, or widen "
                         f"the MAF range.")
            order = np.argsort(picked)
            unobs_idx    = np.array(picked, dtype=int)[order]
            unobs_series = np.array(groups, dtype=int)[order]
            n_unobs_causal = len(unobs_idx)
            print(f"  Allelic series: {len(set(groups))} loci x K={K} "
                  f"= {n_unobs_causal} unobserved causal variants", flush=True)
        else:
            # Legacy path, kept literally so existing seeds reproduce bit for bit.
            if var_tau is not None:
                ew = var_tau[eligible].astype(float)
                unobs_idx = rng.choice(eligible, size=n_unobs_causal, replace=False,
                                       p=ew / ew.sum())
            else:
                unobs_idx = rng.choice(eligible, size=n_unobs_causal, replace=False)
            unobs_idx.sort()
            unobs_series = np.arange(len(unobs_idx))   # each variant its own group
        unobs_set = set(unobs_idx.tolist())
        remaining = np.array([j for j in range(p) if j not in unobs_set])
    elif not recurrence_built:
        # Only clear when there is genuinely no unobserved class. Guarding on
        # `not recurrence_built` matters: the recurrence branch above has already
        # populated unobs_idx / n_unobs_causal, and falling through to here would
        # silently wipe them -- the class would be built and then never given
        # effects, leaving the phenotype with no unassayed component at all.
        n_unobs_causal = 0
        unobs_idx = np.zeros(0, dtype=int)
        unobs_series = np.zeros(0, dtype=int)
        remaining = np.arange(p)

    # --bed-snps-per-mb targets the panel SIZE directly.  --bed-coverage's meaning
    # depends on how many variants the coalescent happened to produce, so it is
    # comparable neither across region lengths nor to a real array; SNPs/Mb is
    # both.  For scale: a real array is ~333 SNPs/Mb (1M SNPs / 3 Gb), which is
    # f_typed ~ 0.087 at the ~3800 variants/Mb these simulations produce.
    _target_k = None
    if args.bed_snps_per_mb > 0.0:
        _target_k = int(round(args.bed_snps_per_mb * args.length / 1e6))

    if _target_k is not None or args.bed_coverage < 1.0:
        bed_rng = np.random.default_rng(args.seed + 1_000_003)
        if _target_k is None and args.bed_freq_weight == 0.0:
            # Uniform thinning (prior behaviour).  Kept as the literal original
            # call so the bed RNG stream is untouched when the new knobs are off.
            bed_mask = bed_rng.random(p) < args.bed_coverage
        else:
            # Fixed-size panel sampled without replacement, weighted by
            # [2p(1-p)]^w when --bed-freq-weight > 0 and uniform otherwise.
            _af  = dosage.mean(axis=1) / 2.0
            _maf = np.minimum(_af, 1.0 - _af)
            _het = 2.0 * _maf * (1.0 - _maf)
            _elig = np.flatnonzero(_het > 0.0)   # monomorphic sites are ineligible
            # Locally-sparse markers: drop candidates near an unobserved causal
            # variant BEFORE the draw, not after.  Doing it before keeps _k -- and
            # hence the panel size and overall density -- unchanged: the markers
            # are redistributed to the rest of the region rather than deleted.
            # That is what makes the contrast attributable to LOCAL depletion
            # rather than to a thinner panel, which would confound it with the
            # density axis this experiment is trying to escape.
            if args.unobs_mask_kb > 0.0 and len(unobs_idx):
                _posv = np.asarray(positions, dtype=float)
                _win  = args.unobs_mask_kb * 1e3
                _near = np.abs(_posv[_elig][:, None] - _posv[unobs_idx][None, :]).min(axis=1)
                _keep = _near > _win
                _dropped = int((~_keep).sum())
                # The panel size that WOULD have been achieved without masking:
                # the request, clamped by how many eligible sites exist at all.
                # Comparing against the unclamped request would reject regions
                # that were already variant-limited, which masking did not cause.
                _want = _target_k if _target_k is not None else int(round(args.bed_coverage * p))
                _want = max(1, min(_want, len(_elig)))
                if _keep.sum() < _want:
                    # The panel could only be filled by re-using masked sites, so
                    # the "same density, locally depleted" property this knob is
                    # for would silently become "thinner panel" -- which is the
                    # density axis, i.e. exactly the confound the experiment
                    # exists to avoid.  Refuse rather than quietly conflate them.
                    sys.exit(
                        f"ERROR: --unobs-mask-kb {args.unobs_mask_kb} leaves only "
                        f"{int(_keep.sum())} eligible sites but the panel needs {_want}. "
                        f"The masked windows ({len(unobs_idx)} variants x +/-"
                        f"{args.unobs_mask_kb} kb) cover too much of a {args.length/1e3:.0f} kb "
                        f"region. Lengthen --length, lower --unobs-mask-kb or "
                        f"--n-unobs-causal, or reduce --bed-snps-per-mb.")
                _elig = _elig[_keep]
                print(f"  Locally-sparse markers: {_dropped} candidate sites within "
                      f"{args.unobs_mask_kb} kb of an unobserved causal variant removed "
                      f"from the panel pool ({len(_elig)} remain, panel size {_want} "
                      f"unchanged)", flush=True)
            if args.bed_freq_weight == 0.0:
                _pw = None
            else:
                _w  = _het ** args.bed_freq_weight
                _pw = _w[_elig] / _w[_elig].sum()
            _k = _target_k if _target_k is not None else int(round(args.bed_coverage * p))
            _k = max(1, min(_k, len(_elig)))
            bed_mask = np.zeros(p, dtype=bool)
            bed_mask[bed_rng.choice(_elig, size=_k, replace=False, p=_pw)] = True
    else:
        bed_mask = np.ones(p, dtype=bool)

    # Phase 6b: unobserved-class variants are never in the BED, regardless of
    # coverage f.  Force-exclude them after the f-coverage draw.
    if args.h2_unobs > 0:
        bed_mask[unobs_idx] = False

    # An SNV callset contains no structural variants, causal or not. Excluding
    # only the causal ones would leak the class into the panel and understate
    # how much of it is unassayed.
    if is_sv.any():
        bed_mask &= ~is_sv

    # Block 2.2: common-SNP / array-like typed panel.  Drop rare variants from
    # the BED so they are untyped (only the ARG can see them).
    if args.bed_maf_min > 0.0:
        freq_all = dosage.mean(axis=1) / 2.0
        maf_all  = np.minimum(freq_all, 1.0 - freq_all)
        bed_mask &= (maf_all >= args.bed_maf_min)

    # ---- SV-class validation -------------------------------------------------
    # The mechanism claims an SV is untaggable when no SNV landed on its branch.
    # Two checks, both reported rather than assumed:
    #   exact  -- does any PANEL SNV have precisely this carrier set? (a mutation
    #             on the same branch produces exactly that vector)
    #   max r2 -- realised best tagging, so PARTIAL tagging is visible too
    # Independent branch-weight calculation predicts ~80% untagged; a large
    # departure means the overlay is not behaving as intended.
    sv_diag = None
    if is_sv.any():
        _panel = np.flatnonzero(bed_mask)
        _svv   = np.flatnonzero(is_sv)
        _keys  = {dosage[j].astype(np.int8).tobytes() for j in _panel}
        _exact = np.array([dosage[j].astype(np.int8).tobytes() in _keys for j in _svv])
        _P = dosage[_panel].T
        _Pc = _P - _P.mean(0); _Pn = np.sqrt((_Pc * _Pc).sum(0)); _Pn[_Pn < 1e-12] = 1.0
        _r2 = []
        for j in _svv:
            t = dosage[j] - dosage[j].mean()
            tn = math.sqrt(float(t @ t))
            if tn < 1e-12:
                _r2.append(0.0); continue
            r = (t @ _Pc) / (tn * _Pn)
            _r2.append(float(np.max(r * r)))
        _r2 = np.array(_r2)
        sv_diag = {"n_sv_variants": int(len(_svv)),
                   "frac_no_exact_snv_match": float(1.0 - _exact.mean()),
                   "mean_max_r2_to_panel": float(_r2.mean()),
                   "median_max_r2_to_panel": float(np.median(_r2)),
                   "frac_max_r2_below_0.2": float((_r2 < 0.2).mean())}
        print(f"  SV validation: {100*sv_diag['frac_no_exact_snv_match']:.0f}% have no exact "
              f"panel match (branch-weight prediction ~80%); "
              f"median max r2 to panel {sv_diag['median_max_r2_to_panel']:.3f}, "
              f"{100*sv_diag['frac_max_r2_below_0.2']:.0f}% below r2=0.2", flush=True)

    # ---- phi, MEASURED -------------------------------------------------------
    # phi = sum_j beta_j^2 (1 - rho_j^2) / sum_j beta_j^2, rho_j = max r^2 to the
    # panel.  It is not binary: a causal variant tagged at r2 = 0.88 still leaves
    # 12% of its variance unassayed.  Computed over ALL causal variants (assayed
    # and unassayed) after effects are drawn, so phi is an OUTPUT of a cell, not
    # a knob -- which is how it should be quoted.  Filled in below once betas
    # exist; the panel-side machinery is set up here.
    _phi_panel = np.flatnonzero(bed_mask)

    if not bed_mask.any():
        sys.exit("ERROR: bed mask produced an empty BED panel; raise --bed-coverage")
    n_bed = int(bed_mask.sum())
    if args.bed_coverage < 1.0 or args.h2_unobs > 0:
        print(f"  Keeping {n_bed}/{p} variants in BED ({100*n_bed/p:.1f}%)", flush=True)

    # Standardize dosage to match SARGENT's GENO_STANDARDIZE (zero-mean, unit-variance).
    # betas are therefore in units of per-s.d. of genotype, consistent with inference.
    dmean = dosage.mean(axis=1, keepdims=True)
    dstd  = dosage.std(axis=1, keepdims=True)
    dstd[dstd == 0] = 1.0   # avoid divide-by-zero for monomorphic sites
    dosage_std = (dosage - dmean) / dstd

    # ---- 3. Fixed effects --------------------------------------------------
    sigma2_mix = np.array(args.sigma2_mix, dtype=float)
    K = len(sigma2_mix)
    if sigma2_mix[0] != 0.0:
        sys.exit("ERROR: first sigma2_mix component must be 0 (null)")

    pi_mix = np.array(args.pi_mix if args.pi_mix else [1.0/K]*K, dtype=float)
    if len(pi_mix) != K:
        sys.exit(f"ERROR: pi_mix length {len(pi_mix)} != sigma2_mix length {K}")
    pi_mix /= pi_mix.sum()

    # Resolve n_causal: explicit --n-causal wins; otherwise derive from
    # --causal-fraction.  Floor at 20 SNPs so very small simulations
    # don't get pathologically few causal variants.  unobs_idx / remaining
    # are computed above (before bed_mask) so the SV-class carve-out is
    # disjoint from the SNP causal set.
    n_causal_floored = False
    if args.n_causal is not None:
        n_causal_target = min(args.n_causal, len(remaining))
    else:
        raw = int(round(args.causal_fraction * p))
        floor = max(0, args.min_causal)
        n_causal_target = min(max(floor, raw), len(remaining))
        n_causal_floored = raw < floor
        if n_causal_floored:
            sys.stderr.write(
                f"WARNING: --causal-fraction {args.causal_fraction} x {p} polymorphic "
                f"variants = {raw} < --min-causal {floor}; using {n_causal_target} causal "
                f"SNPs (realised causal fraction {n_causal_target / p:.4f}, "
                f"intended {args.causal_fraction}).\n")
    n_causal = n_causal_target
    # Block 2.2 — annotation enrichment via causal DENSITY: a variant's causal
    # probability is proportional to tau_{ann(j)}, so enriched annotations carry
    # tau-fold more causal variants per base (constant effect sizes).  This makes
    # per-annotation h^2 proportional to base_frac * tau while keeping the signal
    # polygenic (many modest-effect SNPs), which the downstream estimators need.
    # SIMULATION_SPECS §5: restrict the ELIGIBLE set before the draw, so both the
    # annotation-weighted and the plain branch below see the same pool.  This
    # necessarily changes what rng.choice returns, hence the RNG-stream caveat in
    # --causal-maf-min's help.
    if args.causal_maf_min > 0.0:
        _af  = dmean.ravel() / 2.0            # dosage is 0/1/2, so its mean is 2p
        _maf = np.minimum(_af, 1.0 - _af)
        # var_tau below is indexed by variant id, so filtering `remaining` keeps
        # the annotation weights aligned with no extra bookkeeping.
        remaining = remaining[_maf[remaining] >= args.causal_maf_min]
        if len(remaining) < n_causal:
            sys.exit(
                f"ERROR: only {len(remaining)} variants have MAF >= "
                f"{args.causal_maf_min}, but --n-causal {n_causal} are required.\n"
                f"  Lower --causal-maf-min, lower --n-causal, or simulate a longer "
                f"sequence (of {p} polymorphic variants, "
                f"{int((_maf >= args.causal_maf_min).sum())} clear the floor).")

    if var_tau is not None:
        w = var_tau[remaining].astype(float)
        causal_idx = rng.choice(remaining, size=n_causal, replace=False,
                                p=w / w.sum())
    else:
        causal_idx = rng.choice(remaining, size=n_causal, replace=False)
    causal_idx.sort()

    # --causal-not-typed: force causal variants out of the BED, mirroring what
    # Phase 6b already does for the unobserved class.  Without this, a coverage
    # axis does not measure LD TAGGING: --causal-maf-min and --bed-freq-weight
    # both select on frequency, so ascertainment preferentially types exactly the
    # variants made causal (measured: at f_typed=0.3 the panel directly genotyped
    # causal variants carrying 99% of the genetic variance).  Applied here, after
    # causal selection and before bed_mask is used, so in_bed and n_causal_in_bed
    # report the truth.
    if args.causal_not_typed:
        _n_before = n_bed
        bed_mask[causal_idx] = False
        if not bed_mask.any():
            sys.exit("ERROR: --causal-not-typed emptied the BED panel; raise "
                     "--bed-snps-per-mb / --bed-coverage")
        # n_bed was computed from bed_mask above, before causal_idx existed, and
        # write_plink iterates range(n_bed) over bed_dosage -- so it MUST be
        # refreshed here or the panel and its declared size disagree.
        # Dropping causal variants shrinks the panel below the size that was
        # requested, so without a top-up each density point would silently mean a
        # different density (at 1400 SNPs/Mb the panel came back 11% short).  Refill
        # from the eligible non-causal variants using the same ascertainment
        # weights, so --bed-snps-per-mb stays exact and comparable across the axis.
        _short = _n_before - int(bed_mask.sum())
        if _short > 0:
            _af2  = dosage.mean(axis=1) / 2.0
            _maf2 = np.minimum(_af2, 1.0 - _af2)
            _het2 = 2.0 * _maf2 * (1.0 - _maf2)
            _avail = np.flatnonzero((_het2 > 0.0) & (~bed_mask))
            _avail = np.setdiff1d(_avail, causal_idx)
            if args.h2_unobs > 0:
                _avail = np.setdiff1d(_avail, unobs_idx)
            if args.bed_maf_min > 0.0:
                _avail = _avail[_maf2[_avail] >= args.bed_maf_min]
            _take = min(_short, len(_avail))
            if _take > 0:
                _rng2 = np.random.default_rng(args.seed + 1_000_009)
                if args.bed_freq_weight == 0.0:
                    _pw2 = None
                else:
                    _w2  = _het2[_avail] ** args.bed_freq_weight
                    _pw2 = _w2 / _w2.sum()
                bed_mask[_rng2.choice(_avail, size=_take, replace=False, p=_pw2)] = True
        n_bed = int(bed_mask.sum())
        print(f"  --causal-not-typed: dropped {_n_before - int(bed_mask.sum()) + _short} "
              f"causal variant(s), refilled {_short}; {n_bed} typed", flush=True)

    beta = np.zeros(p)
    for j in causal_idx:
        k = rng.choice(K, p=pi_mix)
        if sigma2_mix[k] > 0:
            beta[j] = rng.normal(0.0, math.sqrt(sigma2_mix[k]))

    # SIMULATION_SPECS §5: impose the frequency dependence of per-allele effects.
    # beta above is on the STANDARDISED-genotype scale, i.e. alpha = -1 already.
    # Per-allele beta_j = beta_std_j / dstd_j, so Var(beta_per_allele) ∝ dstd^-2 *
    # Var(beta_std); to reach Var(beta_per_allele) ∝ (dstd^2)^alpha we need
    # Var(beta_std) ∝ (dstd^2)^(alpha+1), i.e. scale beta_std by dstd^(alpha+1).
    # dstd is used rather than a recomputed sqrt(2pq) because it is exactly what
    # the standardisation divided by, so no HWE assumption enters.  Deterministic:
    # consumes no RNG, and the exponent is 0 (a no-op) at the alpha=-1 default.
    if args.alpha != -1.0:
        beta[causal_idx] *= dstd.ravel()[causal_idx] ** (args.alpha + 1.0)

    # Block 2.1 — age-stratified effect-size scaling.
    # β_j² is scaled by (age_j / age_ref)^tau before h²_snp normalisation.
    # The downstream scale_f rescaling preserves total h²_snp; tau only
    # changes the SHAPE of the effect-size distribution across variant ages.
    if args.age_tau_coef != 0.0 and args.h2_snp > 0:
        ref = max(args.age_tau_ref, 1.0)
        # Clip ages to avoid 0^negative = inf for variants at the leaves.
        # Use 1 generation as the floor.
        eff_ages = np.clip(ages, 1.0, None) / ref
        # β scales as sqrt of variance scaling → exponent is tau/2.
        beta_scale = eff_ages ** (args.age_tau_coef / 2.0)
        beta *= beta_scale

    # Fixed-effect phenotype contribution (uses standardized genotypes)
    y_fixed = dosage_std.T @ beta           # (n_indiv,)

    # Trait-scale normalization reference. We pin realised heritability on the
    # TEST population (the last n_test individuals = pure target pop) so the
    # prediction target's h² is held at target across mixture levels and
    # cross-ancestry scenarios; otherwise pooled-sample normalization lets the
    # test-pop h² drift with sample composition, biasing r2_calib's ceiling
    # (see SIMULATION_SPECS.md §4.3 in the sargent-paper repo). With no test set
    # (n_test == 0) fall back to the full sample (legacy behaviour).
    norm_idx = slice(args.n_indiv, n_indiv) if args.n_test > 0 else slice(None)

    var_fixed = float(np.var(y_fixed[norm_idx]))
    if var_fixed > 0 and args.h2_snp > 0:
        scale_f = math.sqrt(args.h2_snp / var_fixed)
        beta   *= scale_f
        y_fixed = dosage_std.T @ beta
    elif args.h2_snp == 0:
        beta[:] = 0.0
        y_fixed[:] = 0.0

    # ---- 4. Random effects (tree-branch Brownian motion) ------------------
    tau2_frac = np.array(args.tau2_frac, dtype=float)
    M = len(tau2_frac)
    if tau2_frac[0] != 0.0:
        sys.exit("ERROR: first tau2_frac component must be 0 (null)")

    tau2_mix = tau2_frac / (n_trees if n_trees > 0 else 1)

    phi_mix = np.array(args.phi_mix if args.phi_mix else [1.0/M]*M, dtype=float)
    if len(phi_mix) != M:
        sys.exit(f"ERROR: phi_mix length {len(phi_mix)} != tau2_frac length {M}")
    phi_mix /= phi_mix.sum()

    # Per-tree ws²: either a directly-constructed architecture (Phase 10) or the
    # discrete-mixture draw (legacy / --regime).
    if args.ws2_architecture is not None:
        ws2_vec = build_ws2_architecture(
            args.ws2_architecture, n_trees, args.ws2_n_blocks,
            args.ws2_block_frac, args.ws2_block_contrast,
            rng=rng, ts=ts)
    else:
        mix_draw = rng.choice(M, size=n_trees, p=phi_mix)
        ws2_vec  = tau2_mix[mix_draw]       # shape (n_trees,)

    print("Simulating random effects ...", flush=True)
    y_random = sim_random_effect(ts, hap_ranks, n_indiv, ws2_vec, rng)

    var_random = float(np.var(y_random[norm_idx]))
    if var_random > 0 and args.h2_arg > 0:
        scale_r   = math.sqrt(args.h2_arg / var_random)
        y_random *= scale_r
        # scale the variance parameters accordingly — array or per-edge dict
        if isinstance(ws2_vec, dict):
            ws2_vec = {k: v * (scale_r ** 2) for k, v in ws2_vec.items()}
        else:
            ws2_vec  *= scale_r ** 2
    elif args.h2_arg == 0:
        y_random[:] = 0.0
        if isinstance(ws2_vec, dict):
            ws2_vec = {k: 0.0 for k in ws2_vec}
        else:
            ws2_vec[:]  = 0.0

    # ---- 4b. Unobserved (SV-style) class (Phase 6b; default off) ----------
    # Variants are carved out of the BED and contribute to y via heavier-tail
    # effect sizes (gaussian or laplace scaled by --unobs-effect-scale).  The
    # contribution looks like missing heritability to SARGENT's SNP fit but
    # should be at least partially recoverable by the ARG random component
    # via the linkage structure shared between unobs variants and BED variants.
    y_unobs   = np.zeros(n_indiv)
    beta_unobs = np.zeros(p)
    if args.h2_unobs > 0 and n_unobs_causal > 0:
        # Scale relative to the BayesR-slab SD (= sqrt of largest sigma2_mix).
        slab_sd = math.sqrt(max(sigma2_mix))
        scale   = args.unobs_effect_scale * slab_sd
        if args.unobs_effect_dist == "gaussian":
            eff = rng.normal(0.0, scale, size=n_unobs_causal)
        else:  # laplace
            # Laplace(0, b) has var = 2b² ⇒ b = scale/sqrt(2) to match Gaussian variance.
            eff = rng.laplace(0.0, scale / math.sqrt(2), size=n_unobs_causal)
        if args.unobs_series_k > 1:
            # SAME SIGN WITHIN A LOCUS -- this is the point of an allelic series.
            # Distinct alleles of one gene converge on the same functional
            # consequence, so their effects point the same way; that is why the
            # BURDEN carries signal a single marker does not.  With random signs
            # the K effects would partly cancel and the series would look like
            # noise to every method, ARG included.  Magnitudes keep their draw;
            # only the sign is imposed, one Rademacher draw per locus.
            _sgn = rng.choice([-1.0, 1.0], size=int(unobs_series.max()) + 1)
            eff = np.abs(eff) * _sgn[unobs_series]
        beta_unobs[unobs_idx] = eff
        y_unobs = dosage_std.T @ beta_unobs

        var_unobs = float(np.var(y_unobs[norm_idx]))
        if var_unobs > 0:
            scale_u = math.sqrt(args.h2_unobs / var_unobs)
            beta_unobs *= scale_u
            y_unobs    *= scale_u

    # ---- 4c. Joint rescale: make the TOTAL genetic variance hit its target --
    #
    # Each component above is scaled to its OWN target variance on norm_idx, but
    # the components are correlated -- typed and untyped causal variants sit in
    # the same window and are therefore in LD -- so their covariance is genetic
    # variance that nothing budgeted for.  var(y) then lands above 1 and realised
    # h2 above target.  Measured on the phi grid before this fix: var(y) ran
    # 1.031 to 1.154 and true h2 (= 1 - sigma2/var(y)) ran 0.806 to 0.827
    # against a 0.800 target, with the excess largest where the typed and
    # untyped components were closest in size -- the signature of a cross term.
    #
    # Rescaling the SUM by one common factor makes the total exact while leaving
    # every RATIO between components untouched, so phi_src (a ratio of component
    # variances) is unchanged and only the overall trait scale moves.  ws2_vec is
    # a variance, hence scale_g**2.
    h2_total_target = args.h2_snp + args.h2_arg + args.h2_unobs
    var_g_total = float(np.var((y_fixed + y_random + y_unobs)[norm_idx]))
    if var_g_total > 0 and h2_total_target > 0:
        scale_g = math.sqrt(h2_total_target / var_g_total)
        beta       *= scale_g
        beta_unobs *= scale_g
        y_fixed    *= scale_g
        y_random   *= scale_g
        y_unobs    *= scale_g
        ws2_vec    *= scale_g ** 2

    # ---- 5. Residual noise and final phenotype ----------------------------
    sigma2 = max(0.0, 1.0 - args.h2_snp - args.h2_arg - args.h2_unobs)
    y_noise = rng.normal(0.0, math.sqrt(sigma2) if sigma2 > 0 else 0.0, size=n_indiv)
    y = y_fixed + y_random + y_unobs + y_noise

    # Realised heritabilities. Pooled (all sampled individuals) kept for
    # backward-compat; per-population added because the scale is pinned on the
    # TEST population (§4.3), so train- and test-population h² generally differ.
    def _realised_h2(idx):
        vy = float(np.var(y[idx]))
        if vy <= 0:
            return 0.0, 0.0, 0.0
        return (float(np.var(y_fixed[idx]))  / vy,
                float(np.var(y_random[idx])) / vy,
                float(np.var(y_unobs[idx]))  / vy)

    h2_snp_realised, h2_arg_realised, h2_unobs_realised = _realised_h2(slice(None))
    if args.n_test > 0:
        (h2_snp_realised_test, h2_arg_realised_test,
         h2_unobs_realised_test)  = _realised_h2(slice(args.n_indiv, n_indiv))
        (h2_snp_realised_train, h2_arg_realised_train,
         h2_unobs_realised_train) = _realised_h2(slice(0, args.n_indiv))
    else:
        h2_snp_realised_test = h2_arg_realised_test = h2_unobs_realised_test = None
        h2_snp_realised_train = h2_arg_realised_train = h2_unobs_realised_train = None

    # ---- Block 2.2: realised per-annotation heritability (ground truth) -----
    # Per annotation, the genetic variance contributed by its causal variants,
    # split into typed (in BED) and untyped (not in BED, incl. the unobs class).
    # These are the truth columns the three estimators are scored against:
    #   flavour A / C-typed  vs  h2_typed
    #   flavour B            vs  h2_total
    #   C-untyped (= B - A)  vs  h2_untyped
    # Cross-annotation covariance is ~0 for scattered annotations, so per-anno
    # totals approximately sum to total genetic h^2 (residual reported in writeup).
    h2_by_annotation = None
    if var_ann is not None:
        beta_all = beta + beta_unobs            # disjoint supports (unobs ∉ causal)
        is_causal = beta_all != 0.0
        # LD-free additive partition: per causal variant, the variance it explains
        # alone is beta_j^2 (genotypes standardised to unit variance).  Summing
        # beta^2 per annotation is annotation-additive (sums to total h^2 with no
        # within-annotation LD cross-terms) and is the quantity S-LDSC targets.
        vy = float(np.var(y)) if float(np.var(y)) > 0 else 1.0
        b2 = beta_all ** 2
        h2_by_annotation = {}
        for c, name in enumerate(ann_ids):
            in_c      = (var_ann == c)
            typed_c   = in_c & bed_mask  & is_causal
            untyped_c = in_c & (~bed_mask) & is_causal
            h2_typed   = float(b2[typed_c].sum()   / vy)
            h2_untyped = float(b2[untyped_c].sum() / vy)
            h2_by_annotation[name] = {
                "base_frac":  float(ann_base[c]),
                "tau_c":      float(ann_tau[c]),
                "n_causal":   int((in_c & is_causal).sum()),
                "h2_typed":   h2_typed,
                "h2_untyped": h2_untyped,
                "h2_total":   h2_typed + h2_untyped,
            }

    # ---- 6. Write outputs -------------------------------------------------
    print(f"Writing outputs with prefix '{args.out}' ...", flush=True)

    n_train = args.n_indiv
    n_test  = args.n_test
    has_test = n_test > 0
    ts_train = None

    # Restrict the BED panel to variants chosen by --bed-coverage.
    bed_dosage   = dosage[bed_mask]
    bed_positions = [positions[j] for j in range(p) if bed_mask[j]]

    if has_test:
        # Save full ARG (needed for cross-kinship between train and test in prediction)
        ts.dump(args.out + ".full.trees")
        # Build the train tree sequence by simplifying to train sample nodes
        train_sample_nodes = []
        for i in range(n_train):
            train_sample_nodes.extend(ts.individual(i).nodes)
        ts_train = ts.simplify(samples=train_sample_nodes)
        ts_train.dump(args.out + ".trees")

        # Train: phenotypes + BED-panel genotypes for first n_train individuals
        np.savetxt(args.out + ".phen", y[:n_train], fmt="%.8g")
        write_plink(args.out, bed_dosage[:, :n_train], bed_positions, n_train, n_bed)

        # Test: phenotypes + BED-panel genotypes for last n_test individuals
        np.savetxt(args.out + ".test.phen", y[n_train:], fmt="%.8g")
        write_plink(args.out + ".test", bed_dosage[:, n_train:], bed_positions,
                    n_test, n_bed)

        # Also save full phen for completeness
        np.savetxt(args.out + ".full.phen", y, fmt="%.8g")
        # Per-individual true components (for evaluation of prediction quality)
        with open(args.out + ".test.true_components.tsv", "w") as f:
            f.write("idx\ty_fixed\ty_random\ty_total\n")
            for i in range(n_test):
                gi = n_train + i
                f.write(f"{i}\t{y_fixed[gi]:.8g}\t{y_random[gi]:.8g}\t{y[gi]:.8g}\n")

        # Per-individual ancestry labels (Block 1.1 transfer demography).
        # For 'transfer': first n_train_A train indivs are A, remaining n_train_B
        # train indivs are B, all n_test indivs are B. For other demographies,
        # all indivs share the same population label.
        if args.demography == "transfer":
            with open(args.out + ".train.anc.tsv", "w") as f:
                f.write("idx\tpop\n")
                for i in range(n_train):
                    pop = "A" if i < n_train_A else "B"
                    f.write(f"{i}\t{pop}\n")
            with open(args.out + ".test.anc.tsv", "w") as f:
                f.write("idx\tpop\n")
                for i in range(n_test):
                    f.write(f"{i}\tB\n")
        elif args.demography == "stdpopsim":
            with open(args.out + ".train.anc.tsv", "w") as f:
                f.write("idx\tpop\n")
                for i in range(n_train):
                    pop = args.train_pop if i < n_train_A else args.test_pop
                    f.write(f"{i}\t{pop}\n")
            with open(args.out + ".test.anc.tsv", "w") as f:
                f.write("idx\tpop\n")
                for i in range(n_test):
                    f.write(f"{i}\t{args.test_pop}\n")
    else:
        ts.dump(args.out + ".trees")
        np.savetxt(args.out + ".phen", y, fmt="%.8g")
        write_plink(args.out, bed_dosage, bed_positions, n_train, n_bed)

    # Map full-panel indices to BED-panel indices (for SNPs that are in the BED)
    bed_idx_of = np.full(p, -1, dtype=int)
    bed_idx_of[bed_mask] = np.arange(n_bed)

    # `maf` is APPENDED as the last column on purpose: metrics.py reads this file
    # by header name, but age_partition.py and e5_localized_age.py use positional
    # np.loadtxt, so inserting a column would break them while appending is safe.
    _af_out  = dmean.ravel() / 2.0
    _maf_out = np.minimum(_af_out, 1.0 - _af_out)
    with open(args.out + ".true_beta.tsv", "w") as f:
        f.write("snp_idx\tbp_pos\tbeta\tin_bed\tbed_idx\tage\tmaf\n")
        for j in causal_idx:
            in_bed = 1 if bed_mask[j] else 0
            f.write(f"{j}\t{positions[j]}\t{beta[j]:.8g}\t{in_bed}\t"
                    f"{bed_idx_of[j] if in_bed else -1}\t{ages[j]:.6g}\t"
                    f"{_maf_out[j]:.6g}\n")

    # UNOBSERVED causal variants (2026-09-09).  These are the --h2-unobs class:
    # real causal variants deliberately kept OUT of the typed panel, so a
    # tag-based predictor has to reach them through LD.  Until now they existed
    # only as aggregates in params.json (h2_unobs_realised, n_unobs_causal,
    # phi_measured) and their identities were never written, which made every
    # variant-level analysis at phi > 0 impossible -- the analysis could see only
    # the easy, fully typed half of the causal set.
    #
    # This writer is purely ADDITIVE: it emits a file and consumes no RNG, so the
    # same seed reproduces byte-identical phenotypes and existing fits stay valid.
    if len(unobs_idx):
        with open(args.out + ".true_unobs.tsv", "w") as f:
            f.write("snp_idx\tbp_pos\tbeta\tage\tmaf\n")
            for j in unobs_idx:
                # beta_unobs, NOT beta: the two arrays have disjoint supports
                # (see beta_all = beta + beta_unobs above), so beta[j] is 0 here.
                f.write(f"{j}\t{positions[j]}\t{beta_unobs[j]:.8g}\t"
                        f"{ages[j]:.6g}\t{_maf_out[j]:.6g}\n")

    # PER-INDIVIDUAL TRUE GENETIC VALUE, every individual (2026-10-02).
    #
    # WHY.  Variant-level analyses score a predictor against g, the true genetic
    # value, and rebuilt it from .true_beta.tsv + .true_unobs.tsv by looking each
    # causal variant up in the tree sequence BY POSITION.  That works only when
    # every causal variant is a site.  Under --unobs-recurrence-zipf the
    # --h2-unobs class is a set of COMPOSITES: K independent origins summed into
    # one pseudo-variant and APPENDED to the variant matrix after the real sites
    # (unobs_idx = arange(p, p + len(new_dos)) above).  A composite therefore has
    # no site of its own, and the bp_pos written for it maps to nothing --
    # measured on a Fig 3 cell, 1 of 40 untyped causal variants mapped, by
    # coincidence.  Reconstruction silently returned a g missing the entire
    # untyped class on exactly the cells where that class is the whole point, and
    # a guard keyed on "did any map" passed because of that one coincidence.
    #
    # .test.true_components.tsv did not cover this: its y_fixed is the TYPED class
    # alone, y_unobs was never broken out, and nothing at all was written for the
    # training individuals.
    #
    # Written in .full.phen order (train then test) and AFTER scale_g, so g_total
    # is on the phenotype's scale.  Purely ADDITIVE: it emits arrays that already
    # exist and consumes no RNG, so the same seed reproduces byte-identical
    # phenotypes, panels and trees, and existing fits stay valid.
    with open(args.out + ".true_g.tsv", "w") as f:
        f.write("idx\tg_typed\tg_arg\tg_unobs\tg_total\n")
        for i in range(n_indiv):
            _gt = float(y_fixed[i]) + float(y_random[i]) + float(y_unobs[i])
            f.write(f"{i}\t{float(y_fixed[i]):.8g}\t{float(y_random[i]):.8g}\t"
                    f"{float(y_unobs[i]):.8g}\t{_gt:.8g}\n")

    # ws2 ground truth.  For per-tree architectures (ws2_vec is a numpy array),
    # writes .true_ws2.tsv as before.  For per-edge architectures (ws2_vec is a
    # dict {(tree_idx, child_node) → ws²_e}, Block 2.7.1), writes
    # .true_ws2_per_edge.tsv (sparse per-edge ground truth) AND a derived
    # per-tree summary (sum of per-edge ws²_e · ℓ_e / sum ℓ_e — the
    # average ws² weight in tree s) as .true_ws2.tsv so downstream tools
    # that expect the per-tree format still get a sensible aggregate.
    if isinstance(ws2_vec, dict):
        node_times_full = ts.tables.nodes.time
        # Per-edge file: tree_idx, child_node, parent_age, child_age, length, ws2_e
        with open(args.out + ".true_ws2_per_edge.tsv", "w") as f:
            f.write("tree_idx\tchild_node\tparent_age\tchild_age\tlength\tws2_e\n")
            for ti, tree in enumerate(ts.trees()):
                for u_node in tree.nodes():
                    p = tree.parent(u_node)
                    if p == tskit.NULL:
                        continue
                    w_e = float(ws2_vec.get((ti, int(u_node)), 0.0))
                    if w_e == 0.0:
                        continue
                    f.write(f"{ti}\t{int(u_node)}\t"
                            f"{float(node_times_full[p]):.6g}\t"
                            f"{float(node_times_full[u_node]):.6g}\t"
                            f"{float(node_times_full[p]-node_times_full[u_node]):.6g}\t"
                            f"{w_e:.8g}\n")
        # Derived per-tree summary: ws²_tree = (sum ws²_e · ℓ_e) / (sum ℓ_e)
        # over edges in tree s; gives a single per-tree value matching the
        # per-tree-architecture output convention.
        def per_tree_summary(ts_in):
            out = np.zeros(ts_in.num_trees, dtype=float)
            full_breaks = list(ts.breakpoints())
            for ti, tree in enumerate(ts_in.trees()):
                if ts_in is ts:
                    src_ti = ti
                else:
                    mid = 0.5 * (tree.interval.left + tree.interval.right)
                    src_ti = max(0, min(ts.num_trees - 1,
                                np.searchsorted(full_breaks, mid, side="right") - 1))
                num = 0.0; denom = 0.0
                for u_node in tree.nodes():
                    p = tree.parent(u_node)
                    if p == tskit.NULL:
                        continue
                    ell = float(node_times_full[p] - node_times_full[u_node])
                    w_e = float(ws2_vec.get((src_ti, int(u_node)), 0.0))
                    num += w_e * ell
                    denom += ell
                out[ti] = (num / denom) if denom > 0 else 0.0
            return out
        train_summary = per_tree_summary(ts_train) if has_test else per_tree_summary(ts)
        with open(args.out + ".true_ws2.tsv", "w") as f:
            f.write("tree_idx\tws2\n")
            for s, w in enumerate(train_summary):
                f.write(f"{s}\t{w:.8g}\n")
        if has_test:
            full_summary = per_tree_summary(ts)
            with open(args.out + ".full.true_ws2.tsv", "w") as g:
                g.write("tree_idx\tws2\n")
                for s, w in enumerate(full_summary):
                    g.write(f"{s}\t{w:.8g}\n")
    else:
        # Per-tree architecture (legacy):
        # ws2_vec is defined on the FULL ARG's trees (or on the only ARG if no
        # test split).  If there is a test split, we also record train-tree ws²
        # indexed by the simplified train tree sequence; each train tree inherits
        # ws² from the full-ARG tree containing its midpoint.
        with open(args.out + ".true_ws2.tsv", "w") as f:
            f.write("tree_idx\tws2\n")
            if has_test:
                full_breaks = list(ts.breakpoints())
                for s, tr in enumerate(ts_train.trees()):
                    mid = 0.5 * (tr.interval.left + tr.interval.right)
                    full_s = max(0, min(len(ws2_vec) - 1,
                                        np.searchsorted(full_breaks, mid, side="right") - 1))
                    f.write(f"{s}\t{ws2_vec[full_s]:.8g}\n")
                with open(args.out + ".full.true_ws2.tsv", "w") as g:
                    g.write("tree_idx\tws2\n")
                    for s, w in enumerate(ws2_vec):
                        g.write(f"{s}\t{w:.8g}\n")
            else:
                for s, w in enumerate(ws2_vec):
                    f.write(f"{s}\t{w:.8g}\n")

    # ---- Block 2.2: annotation outputs -------------------------------------
    if var_ann is not None:
        with open(args.out + ".annotations.tsv", "w") as f:
            f.write("snp_idx\tbp_pos\tannotation\tin_bed\tis_causal\n")
            for j in range(p):
                f.write(f"{j}\t{positions[j]}\t{ann_ids[var_ann[j]]}\t"
                        f"{1 if bed_mask[j] else 0}\t"
                        f"{1 if beta_all[j] != 0.0 else 0}\n")
        # Full genomic partition (exon/promoter blocks + intergenic complement):
        # downstream flavour A (per-SNP annotation) and flavour B (per-tree map).
        with open(args.out + ".annot_map.tsv", "w") as f:
            f.write("start_bp\tend_bp\tannotation\n")
            cur = 0
            for (s, e, a) in gene_intervals:
                if s > cur:
                    f.write(f"{cur}\t{s}\tintergenic\n")
                f.write(f"{s}\t{e}\t{ann_ids[a]}\n")
                cur = e
            if cur < int(args.length):
                f.write(f"{cur}\t{int(args.length)}\tintergenic\n")

    # Realised MAF of the causal set, for the provenance block below.  Computed
    # here rather than reused from the .true_beta.tsv writer so this does not
    # depend on that block having run.
    _af_params  = dmean.ravel() / 2.0
    _maf_params = np.minimum(_af_params, 1.0 - _af_params)[causal_idx]
    if _maf_params.size == 0:                     # h2_snp == 0: no causal SNPs
        _maf_params = np.zeros(1)

    # ---- measured phi --------------------------------------------------------
    phi_measured = None
    # beta and beta_unobs are both in PER-SD units (y is built from dosage_std),
    # so the variance a causal variant explains is beta^2 -- not beta^2 * Var(x).
    # beta_all is defined only inside a conditional upstream, so combine the two
    # arrays directly rather than relying on it existing.
    _b = np.asarray(beta) + np.asarray(beta_unobs)
    _cidx = np.flatnonzero(_b != 0.0)
    if os.environ.get("SIM_PHI_DEBUG"):
        _u = np.flatnonzero(np.asarray(beta_unobs) != 0.0)
        print(f"  [phi-dbg] p={p} len(beta)={len(beta)} len(beta_unobs)={len(beta_unobs)} "
              f"n_causal_nonzero={_cidx.size} n_unobs_nonzero={_u.size}", flush=True)
        print(f"  [phi-dbg] unobs_idx[:5]={list(unobs_idx[:5]) if args.h2_unobs>0 else []} "
              f"beta_unobs nonzero idx[:5]={list(_u[:5])}", flush=True)
        print(f"  [phi-dbg] sum b^2 total={float((_b**2).sum()):.4f} "
              f"unobs={float((np.asarray(beta_unobs)**2).sum()):.4f} "
              f"typed={float((np.asarray(beta)**2).sum()):.4f}", flush=True)
    if _cidx.size and _phi_panel.size:
        _P = dosage[_phi_panel].T
        _Pc = _P - _P.mean(0); _Pn = np.sqrt((_Pc * _Pc).sum(0)); _Pn[_Pn < 1e-12] = 1.0
        _num = _den = 0.0
        _rows = []
        for j in _cidx:
            t = dosage[j] - dosage[j].mean()
            tn = math.sqrt(float(t @ t))
            rho2 = 0.0
            if tn >= 1e-12:
                r = (t @ _Pc) / (tn * _Pn)
                rho2 = float(np.max(r * r))
            v = float(_b[j]) ** 2
            _num += v * (1.0 - rho2); _den += v
            _rows.append((int(j), rho2, v))
        phi_measured = float(_num / _den) if _den > 0 else None
        _unob = set(unobs_idx.tolist()) if args.h2_unobs > 0 else set()
        _nu = sum(v * (1.0 - r2) for j, r2, v in _rows if j in _unob)
        _na = sum(v * (1.0 - r2) for j, r2, v in _rows if j not in _unob)
        print(f"  phi (measured) = {phi_measured:.4f}   "
              f"[unassayed-class share {_nu/max(_num,1e-12):.2f}, "
              f"typed-class share {_na/max(_num,1e-12):.2f}]", flush=True)

    params = {
        "n_indiv":         n_train,
        "n_test":          n_test,
        "n_total":         n_train + n_test,
        "n_trees":         (ts_train.num_trees if (has_test and ts_train is not None)
                            else n_trees),
        "n_trees_full":    n_trees,
        "n_variants":      n_bed,
        "n_variants_full": p,
        "bed_coverage":    args.bed_coverage,
        "bed_maf_min":     args.bed_maf_min,
        "bed_freq_weight": args.bed_freq_weight,
        "n_causal":        n_causal,
        "n_causal_in_bed": int(bed_mask[causal_idx].sum()),
        # Frequency architecture (SIMULATION_SPECS §5).  The realised quantiles
        # make it checkable after the fact that the floor did what was asked,
        # and record how rare the causal set actually was when it did not.
        "causal_maf_min":  args.causal_maf_min,
        "alpha":           args.alpha,
        "sv_rate":         args.sv_rate,
        "unobs_recurrence_k": args.unobs_recurrence_k,
        "unobs_recurrence_zipf": args.unobs_recurrence_zipf,
        "phi_measured":    phi_measured,
        "n_sv_variants":   int(is_sv.sum()) if is_sv.size else 0,
        "sv_validation":   sv_diag,
        "causal_maf_realised": {
            "min":    float(np.min(_maf_params)),
            "q25":    float(np.percentile(_maf_params, 25)),
            "median": float(np.median(_maf_params)),
            "q75":    float(np.percentile(_maf_params, 75)),
            "max":    float(np.max(_maf_params)),
            "frac_below_0.01": float((_maf_params < 0.01).mean()),
        },
        "n_polymorphic":            int(p),
        "causal_fraction_target":   float(args.causal_fraction),
        "causal_fraction_realised": (float(n_causal) / p) if p else 0.0,
        "n_causal_floored":         bool(n_causal_floored),
        "length":          args.length,
        "recomb_rate":     args.recomb_rate,
        "mut_rate":        args.mut_rate,
        "Ne":              args.Ne,
        "demography":      args.demography,
        "admix_split_time": args.admix_split_time if args.demography == "admixed" else None,
        "admix_time":       args.admix_time       if args.demography == "admixed" else None,
        "admix_alpha":      args.admix_alpha      if args.demography == "admixed" else None,
        "transfer_split_time": args.transfer_split_time if args.demography == "transfer" else None,
        "mix_fraction_B":     args.mix_fraction_B     if args.demography in ("transfer", "stdpopsim") else None,
        "n_train_A":          n_train_A               if args.demography in ("transfer", "stdpopsim") else None,
        "n_train_B":          n_train_B               if args.demography in ("transfer", "stdpopsim") else None,
        "stdpopsim_model":    args.stdpopsim_model    if args.demography == "stdpopsim" else None,
        "train_pop":          args.train_pop           if args.demography == "stdpopsim" else None,
        "test_pop":           args.test_pop            if args.demography == "stdpopsim" else None,
        "train_pop_time":     args.train_pop_time      if args.demography == "stdpopsim" else None,
        "test_pop_time":      args.test_pop_time       if args.demography == "stdpopsim" else None,
        "fst_train_test":     fst_train_test,
        "seed":            args.seed,
        "h2_snp_target":   args.h2_snp,
        "h2_arg_target":   args.h2_arg,
        "h2_unobs_target": args.h2_unobs,
        "age_tau_coef":    args.age_tau_coef,
        "age_tau_ref":     args.age_tau_ref if args.age_tau_coef != 0 else None,
        "h2_snp_realised":   h2_snp_realised,
        "h2_arg_realised":   h2_arg_realised,
        "h2_unobs_realised": h2_unobs_realised,
        "h2_norm_population":     "test" if args.n_test > 0 else "pooled",
        "h2_snp_realised_test":   h2_snp_realised_test,
        "h2_arg_realised_test":   h2_arg_realised_test,
        "h2_unobs_realised_test": h2_unobs_realised_test,
        "h2_snp_realised_train":  h2_snp_realised_train,
        "h2_arg_realised_train":  h2_arg_realised_train,
        "h2_unobs_realised_train": h2_unobs_realised_train,
        "n_unobs_causal":    n_unobs_causal,
        "unobs_effect_dist":  args.unobs_effect_dist if args.h2_unobs > 0 else None,
        "unobs_effect_scale": args.unobs_effect_scale if args.h2_unobs > 0 else None,
        "unobs_series_k":  args.unobs_series_k  if args.h2_unobs > 0 else None,
        "unobs_series_kb": args.unobs_series_kb if args.h2_unobs > 0
                                                  and args.unobs_series_k > 1 else None,
        "unobs_mask_kb":   args.unobs_mask_kb   if args.h2_unobs > 0 else None,
        "n_unobs_loci":    (int(len(set(unobs_series.tolist()))) if args.h2_unobs > 0
                            and len(unobs_series) else None),
        "unobs_maf_min":      args.unobs_maf_min      if args.h2_unobs > 0 else None,
        "unobs_maf_max":      args.unobs_maf_max      if args.h2_unobs > 0 else None,
        "sigma2":          sigma2,
        "regime":          args.regime,
        "ws2_architecture":   args.ws2_architecture,
        "ws2_n_blocks":       args.ws2_n_blocks       if args.ws2_architecture == "sparse" else None,
        "ws2_block_frac":     args.ws2_block_frac     if args.ws2_architecture == "sparse" else None,
        "ws2_block_contrast": args.ws2_block_contrast if args.ws2_architecture == "sparse" else None,
        "sigma2_mix":      list(sigma2_mix),
        "pi_mix":          list(pi_mix),
        "tau2_frac":       list(tau2_frac),
        "tau2_mix":        list(tau2_mix),
        "phi_mix":         list(phi_mix),
        "annotation_genes":        bool(args.annotation_genes),
        "gene_exon_frac":          args.gene_exon_frac          if args.annotation_genes else None,
        "gene_mean_exon_bp":       args.gene_mean_exon_bp       if args.annotation_genes else None,
        "gene_promoter_bp":        args.gene_promoter_bp        if args.annotation_genes else None,
        "gene_tau_exon":           args.gene_tau_exon           if args.annotation_genes else None,
        "gene_tau_promoter_ratio": args.gene_tau_promoter_ratio if args.annotation_genes else None,
        "annotation_ids":          (list(ann_ids) if var_ann is not None else None),
        "h2_by_annotation":        h2_by_annotation,
    }
    with open(args.out + ".params.json", "w") as f:
        json.dump(params, f, indent=2)

    print(f"\nDone.")
    if has_test:
        print(f"  n_train / n_test : {n_train} / {n_test}")
        print(f"  n_trees (train)  : {ts_train.num_trees}")
        print(f"  n_trees (full)   : {n_trees}")
    else:
        print(f"  n_indiv          : {n_indiv}")
        print(f"  n_trees          : {n_trees}")
    print(f"  n_variants       : {p}")
    print(f"  n_causal (fixed) : {n_causal}")
    print(f"  h2_snp  target / realised : {args.h2_snp:.3f} / {h2_snp_realised:.4f}")
    print(f"  h2_arg  target / realised : {args.h2_arg:.3f} / {h2_arg_realised:.4f}")
    if args.h2_unobs > 0:
        print(f"  h2_unobs target / realised : {args.h2_unobs:.3f} / {h2_unobs_realised:.4f}  "
              f"(n_unobs_causal={n_unobs_causal}, dist={args.unobs_effect_dist}, "
              f"scale={args.unobs_effect_scale}, "
              f"MAF=[{args.unobs_maf_min}, {args.unobs_maf_max}])")
    print(f"  sigma2           : {sigma2:.4f}")


if __name__ == "__main__":
    main()
